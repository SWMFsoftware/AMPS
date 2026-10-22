"""Activity field on the device (kernels/activity.cu; SCALING_PLAN 3.5, invariant I7, decision Q3).

* the device act_lin is within the summation bound of the host expression on every facet
  (``CudaBackend.activity_check``), checked at the first evaluation and every 100th;
* device-activity pressures and Jacobian columns equal the host-activity ones to rounding;
* the device activity is deterministic and partition-independent: with and without face
  compaction, resident or streamed, on two GPUs, after a basis rebind, the results are BITWISE
  those of one resident GPU with device activity (FP64 and FP32 records);
* a corrupted device basis is caught by the check (BackendError) instead of being used.
"""

import numpy as np
import pytest

from cops_cuda import backend as bk
from cops_cuda.errors import BackendError
import test_backend_contract as contract

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import evaluator, masks  # noqa: E402
from test_gpu_model import make_gpu_model  # noqa: E402

cpu_model = contract.cpu_model
NC31 = contract.CONFIGS[-1]


def same(u, v):
    return np.asarray(u).tobytes() == np.asarray(v).tobytes()


def problem(env):
    G = make_gpu_model(env)
    P = contract.make_problem(dict(env, model=G), NC31)
    request = bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True)
    states = [bk.ParameterState.from_problem(P, contract.perturbed(P, seed)) for seed in (1, 2, 3)]
    return G, P, request, states


def run_all(be, states, request, rebind=None):
    """Jacobians and pressures; then, with ``rebind=(basis, states, request)``, the same after
    a basis rebind."""
    out = [be.evaluate(s, request) for s in states]
    out += [be.evaluate(s) for s in states]                     # pressures only
    if rebind is not None:
        basis, states2, request2 = rebind
        be.rebind_basis(basis, 0.7)
        out += [be.evaluate(s, request2) for s in states2]
    return out


def check_same(ref, got):
    for x, y in zip(ref, got):
        assert same(x.pNG, y.pNG) and same(x.pRG, y.pRG)
        for key in x.columns:
            assert same(x.columns[key][0], y.columns[key][0]), key
            assert same(x.columns[key][1], y.columns[key][1]), key


def test_device_activity_is_within_the_bound(cpu_model, clean_pool):
    G, P, request, states = problem(cpu_model)
    run = cpu_model["synthetic"]["run"]
    out = {}
    for activity in ("host", "device"):
        with evaluator.CudaBackend(G, run, P.enhancement_basis, P.enhancement_scale,
                                   records_per_block=17, activity=activity) as be:
            out[activity] = [be.evaluate(s, request) for s in states]
            if activity == "device":
                checks = [be.activity_check(s) for s in states]
                ident = be.identity()["activity"]
                assert ident["mode"] == "device" and ident["checks"]["checks"] == 1 + len(states)
                assert ident["checks"]["evaluations"] == len(states)
    for c in checks:
        assert c["faces"] > 0 and c["max_error_eps"] <= evaluator.ACTIVITY_CHECK_EPS
        assert c["active_flips"] == 0
    for h, d in zip(out["host"], out["device"]):
        for u, v in [(h.pNG, d.pNG), (h.pRG, d.pRG)] + [(h.columns[k][i], d.columns[k][i])
                                                         for k in h.columns for i in (0, 1)]:
            u, v = np.asarray(u), np.asarray(v)
            assert np.all(np.abs(u - v) <= 1e-12 * np.maximum(np.abs(u), np.abs(v)) + 1e-300)


@pytest.mark.parametrize("precision", ["fp64", "fp32"])
def test_device_activity_is_layout_and_partition_invariant(cpu_model, clean_pool, monkeypatch, precision):
    G, P, request, states = problem(cpu_model)
    run = cpu_model["synthetic"]["run"]
    nf = len(cpu_model["elements"]["area"])
    basis2 = np.random.default_rng(9).normal(size=(nf, 3))
    request2 = bk.ColumnRequest(P.fit_j, range(3), temperature=True, kappa=True)
    states2 = [bk.ParameterState(s.a, s.T, s.kappa, np.array([0.3, -0.2, 0.5])) for s in states]
    kw = dict(records_per_block=17, activity="device", record_precision=precision)

    def evaluate(**extra):
        with evaluator.CudaBackend(G, run, P.enhancement_basis, P.enhancement_scale,
                                   **dict(kw, **extra)) as be:
            res = run_all(be, states, request, (basis2, states2, request2))
            return res, masks.seen_positive(be.dataset), be.dataset.streamed

    ref, seen, streamed = evaluate()
    assert not streamed
    got, seen2, _ = evaluate(compact_faces=False)
    check_same(ref, got)
    assert same(seen, seen2)
    if precision == "fp32":
        monkeypatch.setenv(evaluator.STREAM_CHUNK_ENVIRONMENT, "40")
        monkeypatch.setenv(evaluator.RESIDENT_ENVIRONMENT, "120")
        got, seen2, streamed = evaluate()
        assert streamed
        check_same(ref, got)
        monkeypatch.delenv(evaluator.RESIDENT_ENVIRONMENT)
    if cp.cuda.runtime.getDeviceCount() >= 2:
        from cops_cuda.multi_gpu import MultiGpuBackend
        monkeypatch.setenv("COPS_CUDA_RECORD_PRECISION", precision)
        with MultiGpuBackend(G._spec, G, P.enhancement_basis, P.enhancement_scale, devices=(0, 1),
                             records_per_block=17, activity="device") as two:
            assert two.shm is None and two.act is None
            got = run_all(two, states, request, (basis2, states2, request2))
            assert [w["activity"] for w in two.ready] == ["device", "device"]
        check_same(ref, got)


def test_corrupted_device_basis_is_caught(cpu_model, clean_pool):
    G, P, request, states = problem(cpu_model)
    run = cpu_model["synthetic"]["run"]
    with evaluator.CudaBackend(G, run, P.enhancement_basis, P.enhancement_scale,
                               records_per_block=17, activity="device") as be:
        be.arrays["actbasis"][5] += 1.0                          # one facet's Y_0 changed
        with pytest.raises(BackendError, match="device activity differs"):
            be.evaluate(states[0], request)


def test_registry_device_activity(cpu_model, clean_pool):
    from cops_cuda import gpu_model
    reg = gpu_model.REGISTRY
    reg.close()
    reg.__init__()
    try:
        reg.configure(0, "tree", 17, activity="device")
        G = make_gpu_model(cpu_model)
        pNG, pRG = G.evaluate(cpu_model["coef"][0], cpu_model["coef"][1][:3], 1.0)
        want = cpu_model["model"].evaluate(cpu_model["coef"][0], cpu_model["coef"][1][:3], 1.0)
        assert np.allclose(pNG, want[0], rtol=1e-12, atol=0) and np.allclose(pRG, want[1], rtol=1e-12, atol=0)
        (backend, _), = reg.backends.values()
        assert backend.activity_mode == "device" and reg.stats()["activity"] == "device"
    finally:
        reg.close()
        reg.__init__()

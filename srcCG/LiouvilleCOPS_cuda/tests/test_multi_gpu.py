"""Multi-GPU backend (cops_cuda/multi_gpu.py): one vs two devices (needs 2 GPUs).

* partition: contiguous, complete, balanced by records;
* two devices give BITWISE the one-device results (pressures and every column),
  in both modes, through the contract and through gpu_model and an unmodified
  SciPy TRF trajectory;
* seen_faces is the union over workers and equals the reference;
* failure: a dead worker becomes a BackendError, and close() is clean.
"""

import numpy as np
import pytest

from cops_cuda import backend as bk
from cops_cuda.errors import BackendError
from cops_cuda.multi_gpu import partition
import test_backend_contract as contract
from test_gpu_model import make_gpu_model, same

cpu_model = contract.cpu_model


def test_partition_is_contiguous_complete_and_balanced():
    rng = np.random.default_rng(2)
    for n, ndev in ((10, 2), (1150, 2), (7, 3), (2, 2)):
        counts = rng.integers(0, 400000, n)
        blocks = partition(counts, ndev)
        assert blocks[0][0] == 0 and blocks[-1][1] == n
        assert all(a < b for a, b in blocks)
        assert all(blocks[k][1] == blocks[k + 1][0] for k in range(ndev - 1))
        loads = [counts[a:b].sum() for a, b in blocks]
        assert max(loads) - min(loads) <= 2 * counts.max()
    with pytest.raises(ValueError):
        partition([1], 2)


cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() < 2:
    pytest.skip("needs two CUDA devices", allow_module_level=True)
from cops_cuda import evaluator, gpu_model  # noqa: E402
from cops_cuda.multi_gpu import MultiGpuBackend  # noqa: E402


@pytest.fixture
def registry(clean_pool):
    reg = gpu_model.REGISTRY
    reg.close()
    reg.__init__()
    yield reg
    reg.close()
    reg.__init__()


@pytest.mark.parametrize("mode", ["sequential", "tree"])
def test_two_devices_are_bitwise_one_device(cpu_model, registry, mode):
    G = make_gpu_model(cpu_model)
    run = cpu_model["synthetic"]["run"]
    for config in (contract.CONFIGS[0], contract.CONFIGS[-1]):
        P = contract.make_problem(dict(cpu_model, model=G), config)
        request = bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True)
        basis = P.enhancement_basis if P.nb else None
        with evaluator.CudaBackend(G, run, basis, P.enhancement_scale, mode=mode,
                                   records_per_block=17) as one, \
                MultiGpuBackend(G._spec, G, basis, P.enhancement_scale, devices=(0, 1),
                                mode=mode, records_per_block=17) as two:
            assert len(two.blocks) == 2 and sum(two.block_records) == int(G.indptr[-1])
            for seed in (1, 2, 3):
                state = bk.ParameterState.from_problem(P, contract.perturbed(P, seed))
                a, b = one.evaluate(state, request), two.evaluate(state, request)
                assert same(a.pNG, b.pNG) and same(a.pRG, b.pRG)
                for key in request.keys():
                    assert same(a.columns[key][0], b.columns[key][0]), key
                    assert same(a.columns[key][1], b.columns[key][1]), key
            ident = two.identity()
            assert [w["device_id"] for w in ident["workers_ready"]] == [0, 1]
        assert two.closed and two.shm is None


def test_trf_trajectory_on_two_devices_is_bitwise(cpu_model, registry):
    fit_cops = cpu_model["fit_cops"]
    results = []
    for devices in ([0], [0, 1]):
        registry.close()
        registry.__init__()
        registry.configure(devices=devices, mode="tree", records_per_block=17)
        G = make_gpu_model(cpu_model)
        P = contract.make_problem(dict(cpu_model, model=G), dict(contract.CONFIGS[0]))
        P.args.max_nfev = 6
        saved = getattr(fit_cops, "P", None)
        fit_cops.P = P
        try:
            results.append(fit_cops.fit_one(0.0))
        finally:
            fit_cops.P = saved
        if devices == [0, 1]:
            assert same(G.seen_faces(), cpu_model["model"].seen_faces())
    one, two = results
    assert same(one["x"], two["x"]) and one["nfev"] == two["nfev"] > 2
    assert same(one["cost"], two["cost"])


def test_dead_worker_is_a_clean_backend_error(cpu_model, registry):
    G = make_gpu_model(cpu_model)
    P = contract.make_problem(dict(cpu_model, model=G), contract.CONFIGS[0])
    two = MultiGpuBackend(G._spec, G, None, 1.0, devices=(0, 1), timeout_s=60)
    try:
        state = bk.ParameterState.from_problem(P, P.x0)
        two.evaluate(state)
        two.workers[1]["proc"].terminate()
        two.workers[1]["proc"].join(10)
        with pytest.raises(BackendError, match="device 1"):
            two.evaluate(state)
        with pytest.raises(BackendError, match="closed or broken"):
            two.evaluate(state)
    finally:
        two.close()
    assert all(not w["proc"].is_alive() for w in two.workers)
    two.close()                                       # idempotent

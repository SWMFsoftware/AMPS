"""GPU-backed fit_cops.Model (cops_cuda/gpu_model.py) against the reference Model (GPU).

* sequential mode: ``evaluate`` (pressures and every Jacobian entry the
  reference returns) and ``seen_faces`` are bitwise the reference's, and a
  whole SciPy TRF trajectory of ``fit_cops.fit_one`` is therefore bitwise
  identical;
* tree mode: within the approved budgets;
* backends are shared between models with equal rows and rebound only for
  a new enhancement basis, keeping the device records (bitwise a new backend).
"""

import numpy as np
import pytest

from cops_cuda import backend as bk
import test_backend_contract as contract

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import gpu_model  # noqa: E402

cpu_model = contract.cpu_model
EPS = np.finfo(float).eps


@pytest.fixture
def registry(clean_pool):
    reg = gpu_model.REGISTRY
    reg.close()
    reg.__init__()
    yield reg
    reg.close()
    reg.__init__()


def make_gpu_model(env):
    fit_cops = env["fit_cops"]
    run = str(env["synthetic"]["run"])
    header, names, rows = fit_cops.read_model_file(run + "/COPS_model.dat")
    cols = {name: index for index, name in enumerate(names)}
    M = gpu_model.gpu_model_class()(run, header, cols, rows, env["elements"], env["coef"])
    M.load_kernels(M.sel)
    return M


def same(u, v):
    return np.asarray(u).tobytes() == np.asarray(v).tobytes()


def test_sequential_evaluate_and_seen_faces_are_bitwise(cpu_model, registry):
    registry.configure(0, "sequential")
    ref, G = cpu_model["model"], make_gpu_model(cpu_model)
    assert same(G.seen_faces(), ref.seen_faces())
    for config in contract.CONFIGS:
        P = contract.make_problem(cpu_model, config)
        rng = np.random.default_rng(7)
        for _ in range(2):
            x = contract.perturbed(P, int(rng.integers(1000)))
            a, b, T, kappa = P.unpack(x)
            kw = dict(b=b if P.nb else None, enhancement_basis=P.enhancement_basis,
                      enhancement_scale=P.enhancement_scale)
            for got, want in zip(G.evaluate(a, T, kappa, **kw), ref.evaluate(a, T, kappa, **kw)):
                assert same(got, want)
            kw.update(jac=True, basis_index=P.fit_j, enhancement_index=list(range(P.nb)))
            gN, gR, gJ = G.evaluate(a, T, kappa, **kw)
            rN, rR, rJ = ref.evaluate(a, T, kappa, **kw)
            assert same(gN, rN) and same(gR, rR)
            assert set(gJ) == set(rJ)
            for key in rJ:
                assert same(gJ[key][0], rJ[key][0]) and same(gJ[key][1], rJ[key][1]), key


def test_sequential_trf_trajectory_is_bitwise(cpu_model, registry):
    registry.configure(0, "sequential")
    fit_cops = cpu_model["fit_cops"]
    G = make_gpu_model(cpu_model)
    results = []
    for model in (cpu_model["model"], G):
        P = contract.make_problem(dict(cpu_model, model=model), dict(contract.CONFIGS[0]))
        P.args.max_nfev = 6
        P.lam = 0.0
        saved = getattr(fit_cops, "P", None)
        fit_cops.P = P                 # fit_one reads the module-level problem, as main() sets it
        try:
            results.append(fit_cops.fit_one(0.0))
        finally:
            fit_cops.P = saved
    cpu, gpu = results
    assert same(cpu["x"], gpu["x"])
    assert cpu["nfev"] == gpu["nfev"] and cpu["status"] == gpu["status"]
    assert same(cpu["cost"], gpu["cost"])
    assert cpu["nfev"] > 2                       # a real trajectory, not a single step
    assert registry.calls["jacobian"] > 0 and registry.calls["pressures"] > 0


def test_tree_mode_within_budget(cpu_model, registry):
    registry.configure(0, "tree", records_per_block=17)
    ref, G = cpu_model["model"], make_gpu_model(cpu_model)
    P = contract.make_problem(cpu_model, contract.CONFIGS[0])
    x = contract.perturbed(P, 11)
    a, b, T, kappa = P.unpack(x)
    import test_evaluator
    mNG, mRG = test_evaluator.magnitudes(P, x)
    gN, gR = G.evaluate(a, T, kappa)
    rN, rR = ref.evaluate(a, T, kappa)
    assert np.all(np.abs(gN - rN) <= 512 * EPS * mNG)
    assert np.all(np.abs(gR - rR) <= 512 * EPS * mRG)


def test_backends_are_shared_and_rebound_for_a_new_basis(cpu_model, registry):
    registry.configure(0, "sequential")
    G1, G2 = make_gpu_model(cpu_model), make_gpu_model(cpu_model)
    a, T = cpu_model["coef"][0], G1.T0
    G1.evaluate(a, T, 1.0)
    G2.evaluate(a, T, 1.0)
    assert [e["event"] for e in registry.events] == ["bind"]          # shared dataset
    nf = len(cpu_model["elements"]["area"])
    basis = np.ones((nf, 1))
    G2.evaluate(a, T, 1.0, b=np.array([0.1]), enhancement_basis=basis, enhancement_scale=2.0)
    G1.evaluate(a, T, 1.0, b=np.array([0.1]), enhancement_basis=basis.copy(), enhancement_scale=2.0)
    assert [e["event"] for e in registry.events] == ["bind", "rebind_basis_only"]   # records kept
    with pytest.raises(ValueError):
        G1.evaluate(a, T, 1.0, b=np.array([0.1]))                     # reference refusal
    stats = registry.stats()
    assert stats["calls"]["pressures"] == 4 and stats["backends"]


@pytest.mark.parametrize("precision", ["fp64", "fp32"])
def test_basis_rebind_is_bitwise_a_new_backend(cpu_model, clean_pool, precision):
    """SCALING_PLAN P1: a template basis change keeps the records; results equal a fresh bind."""
    from cops_cuda import evaluator
    P = contract.make_problem(cpu_model, next(c for c in contract.CONFIGS if c.get("templates")))
    run = cpu_model["synthetic"]["run"]
    x = contract.perturbed(P, 5)
    state = bk.ParameterState.from_problem(P, x)
    request = bk.ColumnRequest.from_problem(P)
    with evaluator.from_problem(P, run, mode="tree", record_precision=precision) as fresh:
        want, want_j = fresh.evaluate(state), fresh.evaluate(state, request)
    other = np.ones_like(P.enhancement_basis[:, :1])
    with evaluator.CudaBackend(P.M, run, other, 3.0, mode="tree", record_precision=precision) as be:
        dataset = be.dataset
        be.rebind_basis(None)
        assert be.nbasis_enhancement == 0
        be.rebind_basis(P.enhancement_basis, P.enhancement_scale)
        assert be.dataset is dataset and be.timing["rebind_basis_s"] >= 0
        got, got_j = be.evaluate(state), be.evaluate(state, request)
    assert same(got.pNG, want.pNG) and same(got.pRG, want.pRG)
    for key in request.keys():
        assert same(got_j.columns[key][0], want_j.columns[key][0])
        assert same(got_j.columns[key][1], want_j.columns[key][1])

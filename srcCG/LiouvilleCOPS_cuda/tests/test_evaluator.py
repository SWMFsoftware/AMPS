"""C06 CUDA evaluator against the CPU reference through the C02 contract (GPU).

* "sequential" mode must reproduce ``fit_cops.Model.evaluate`` BITWISE: pressures
  and every requested Jacobian column (activity, template, temperature, κ).
  That proves each per-contribution expression and its order are the
  reference's.
* "tree" mode (production) must satisfy the approved budget: pressure
  |Δp| <= 512 eps M (M = sum of contribution magnitudes) and, through the
  unmodified ``fit_cops.Problem``, the derived residual bound.
"""

import numpy as np
import pytest

from cops_cuda import backend as bk
import test_backend_contract as contract

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import evaluator  # noqa: E402
from c06_jacobian_budget import jacobian_magnitudes  # noqa: E402

EPS = np.finfo(float).eps
cpu_model = contract.cpu_model
CONFIGS = contract.CONFIGS


def fresh_pool():
    pool = cp.get_default_memory_pool()
    pool.free_all_blocks()
    pool.set_limit(size=0)


def magnitudes(P, x):
    """Per-point sums of |contribution| for NG and RG at x (the M of the budget)."""
    M = P.M
    a, b, T, kappa = P.unpack(x)
    act = M.Y @ a
    if P.nb:
        act = act + P.enhancement_scale * (P.enhancement_basis @ b)
    af = np.maximum(act, 0.0)[M.face]
    Tmax, T75, Tmin = T
    Tl = Tmax + M.__class__.__init__.__globals__["KT"] * (Tmax - T75) * M.xlit
    u = 1.0 / np.sqrt(np.where(M.lit & (Tl > Tmin), Tl, Tmin))
    bc = lambda w: np.bincount(M.pt, weights=np.abs(w), minlength=M.npt)
    c, r, bm = M.cNG, M.cRG, M.bm
    mNG = abs(c["H2O"] * bm["H2O"]) * bc(M.wn["H2O"] * af * u) + \
        abs(kappa * c["CO2"] * bm["CO2"]) * bc(M.wn["CO2"] * u)
    mRG = np.abs(r["H2O"]) * (bc(M.wp["H2O"] * af) + bm["H2O"] * bc(M.wv["H2O"] * af * u)) + \
        np.abs(kappa * r["CO2"]) * (bc(M.Q["CO2"] * M.ACO2[M.face] * M.KPhi) + bm["CO2"] * bc(M.wv["CO2"] * u))
    return mNG, mRG


def test_host_constants_equal_the_model(cpu_model):
    M = cpu_model["model"]
    c = evaluator.HostConstants(M)
    for s in ("H2O", "CO2"):
        qmin, qmax = c.qmin[s][M.pt], c.qmax[s][M.pt]
        expected = np.where(M.night, qmin, qmin + (qmax - qmin) * np.minimum(M.cos, 1.0))
        assert expected.tobytes() == M.Q[s].tobytes()


@pytest.mark.parametrize("config", CONFIGS, ids=[c["name"] for c in CONFIGS])
def test_sequential_mode_is_bitwise_the_cpu_reference(cpu_model, config):
    fresh_pool()
    P = contract.make_problem(cpu_model, config)
    ref = bk.ReferenceBackend.from_problem(P)
    request = bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True)
    with evaluator.from_problem(P, cpu_model["synthetic"]["run"], mode="sequential") as gpu:
        assert gpu.PhiCO2.tobytes() == gpu.constants.PhiCO2_cpu.tobytes()
        for seed in (1, 2, 3):
            state = bk.ParameterState.from_problem(P, contract.perturbed(P, seed))
            want = ref.evaluate(state, request)
            got = gpu.evaluate(state, request)
            assert got.pNG.tobytes() == want.pNG.tobytes()
            assert got.pRG.tobytes() == want.pRG.tobytes()
            for key in request.keys():
                for g, w in zip(got.columns[key], want.columns[key]):
                    assert g.tobytes() == w.tobytes(), key


@pytest.mark.parametrize("config", [CONFIGS[0], CONFIGS[-1]], ids=["B29", "NC31"])
def test_tree_mode_within_budget_and_drop_in_for_problem(cpu_model, config):
    fresh_pool()
    P = contract.make_problem(cpu_model, config)
    x = contract.perturbed(P, 4)
    r_cpu, J_cpu = P.data_residuals(x), P.data_jacobian(x)
    mNG, mRG = magnitudes(P, x)
    with evaluator.from_problem(P, cpu_model["synthetic"]["run"], mode="tree",
                                records_per_block=17) as gpu:
        state = bk.ParameterState.from_problem(P, x)
        ref = bk.ReferenceBackend.from_problem(P).evaluate(state)
        got = gpu.evaluate(state)
        assert np.all(np.abs(got.pNG - ref.pNG) <= 512 * EPS * mNG)
        assert np.all(np.abs(got.pRG - ref.pRG) <= 512 * EPS * mRG)
        repeat = gpu.evaluate(state)
        assert repeat.pNG.tobytes() == got.pNG.tobytes()      # deterministic
        # Approved Jacobian budget (VALIDATION.md): |dJ| <= 1024 eps M_J per entry.
        request = bk.ColumnRequest.from_problem(P)
        ref_j = bk.ReferenceBackend.from_problem(P).evaluate(state, request)
        got_j = gpu.evaluate(state, request)
        fit_cops = cpu_model["fit_cops"]
        for key, mags in jacobian_magnitudes(P, P.M, x, request, fit_cops).items():
            for g, w, m in zip(got_j.columns[key], ref_j.columns[key], mags):
                assert np.all(np.abs(g - w) <= 1024 * EPS * m), key
        bk.install(P, gpu)
        r_gpu, J_gpu = P.data_residuals(x), P.data_jacobian(x)
    # Derived residual bound (VALIDATION.md): |dr_i| <= w * 512 eps M_i / p_i.
    # A point that sees no facet has p = M = 0 (the pressure floor applies on
    # both sides identically), so its bound is 0, not 0/0.
    floor = 1e-300
    pNG, pRG = np.maximum(ref.pNG[P.useNG], floor), np.maximum(ref.pRG[P.useRG], floor)
    delta = np.concatenate([P.wNG * 512 * EPS * mNG[P.useNG] / pNG,
                            P.wRG * 512 * EPS * mRG[P.useRG] / pRG])
    assert np.all(np.abs(r_gpu - r_cpu) <= delta)
    assert J_gpu.shape == J_cpu.shape
    scale = np.max(np.abs(J_cpu), axis=0)
    assert np.all(np.abs(J_gpu - J_cpu) <= 1e-11 * scale)    # diagnostic (former provisional gate)


def test_unsupported_sza_mode_is_refused(cpu_model):
    M = cpu_model["model"]
    old = M.SZAMode
    try:
        M.SZAMode = "srccg_table"
        with pytest.raises(bk.UnsupportedConfiguration):
            evaluator.HostConstants(M)
    finally:
        M.SZAMode = old


@pytest.mark.parametrize("records_per_block", [17, 4096])
def test_fused_kernels_are_bitwise_the_generic_tree(cpu_model, records_per_block):
    """C08: the compile-time-specialized kernels change no bit of any sum."""
    fresh_pool()
    from cops_cuda.memory import DeviceBudget
    run = cpu_model["synthetic"]["run"]
    budget = DeviceBudget(0)
    try:
        for config in CONFIGS:
            P = contract.make_problem(cpu_model, config)
            requests = [bk.PRESSURES_ONLY,
                        bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True),
                        bk.ColumnRequest(P.fit_j, range(P.nb), temperature=False, kappa=True),
                        bk.ColumnRequest((), (), temperature=True, kappa=False)]
            with evaluator.from_problem(P, run, mode="tree", fused=False, budget=budget,
                                        records_per_block=records_per_block) as generic, \
                    evaluator.from_problem(P, run, mode="tree", fused=True, budget=budget,
                                           records_per_block=records_per_block) as fused:
                assert fused.identity()["fused_kernels"] and not generic.identity()["fused_kernels"]
                assert fused.PhiCO2.tobytes() == generic.PhiCO2.tobytes()
                for seed in (1, 2):
                    state = bk.ParameterState.from_problem(P, contract.perturbed(P, seed))
                    for request in requests:
                        a, b = generic.evaluate(state, request), fused.evaluate(state, request)
                        assert a.pNG.tobytes() == b.pNG.tobytes(), config["name"]
                        assert a.pRG.tobytes() == b.pRG.tobytes(), config["name"]
                        for key in request.keys():
                            for u, v in zip(a.columns[key], b.columns[key]):
                                assert u.tobytes() == v.tobytes(), (config["name"], key)
    finally:
        budget.close()


def test_fused_plan_covers_the_standard_requests_only():
    plan = evaluator.CudaBackend._fused_plan
    fake = type("F", (), {})()
    codes = list(range(5)) + list(range(5, 14)) + [1000 + c for c in range(25)] + [2000 + c for c in range(25)]
    steps = plan(fake, codes)
    assert [s[0] for s in steps] == ["eval_fused_p14", "eval_fused_c16", "eval_fused_c9"]
    assert sorted(c for s in steps for c in s[2]) == sorted(codes)
    assert plan(fake, [14]) is None                       # PhiCO2 alone -> generic
    assert plan(fake, [0, 1, 2, 3, 4, 1000]) is None      # NG column without its RG pair
    assert [s[0] for s in plan(fake, list(range(5)))] == ["eval_fused_p5"]

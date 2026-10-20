"""Header-only model (cops_cuda/host_model.py) against the reference Model.

* CPU: every per-point attribute ``fit_cops.Problem`` and ``CudaBackend`` read
  is bitwise the reference's; the Problem built on it is the same problem.
* CPU: record-dependent methods refuse loudly.
* GPU: a backend bound from the header model evaluates bitwise like one bound
  from the full model, and the unmodified Problem works through the facade.
"""

import numpy as np
import pytest

from cops_cuda import backend as bk
from cops_cuda import host_model
from cops_cuda.errors import UnsupportedConfiguration
import test_backend_contract as contract

cpu_model = contract.cpu_model
PER_POINT = ("ipoint", "useNG", "useRG", "refNG", "refRG", "TRG", "pNG_run", "pRG_run",
             "T0", "S", "sel", "rows", "indptr", "VdotB", "rh")


@pytest.fixture
def light(cpu_model):
    fit_cops = cpu_model["fit_cops"]
    run = str(cpu_model["synthetic"]["run"])
    header, names, rows = fit_cops.read_model_file(run + "/COPS_model.dat")
    cols = {name: index for index, name in enumerate(names)}
    M = host_model.header_model(run, header, cols, rows, cpu_model["elements"], cpu_model["coef"])
    M.load_kernels(M.sel)
    return M


def test_per_point_data_is_bitwise_the_reference(cpu_model, light):
    ref = cpu_model["model"]
    for name in PER_POINT:
        assert np.asarray(getattr(light, name)).tobytes() == \
            np.asarray(getattr(ref, name)).tobytes(), name
    for name in ("SZAMode", "helio", "TNG", "RGModel", "nfaces", "npt"):
        assert getattr(light, name) == getattr(ref, name), name
    assert light.Y is cpu_model["elements"]["Y"]
    for s in ("H2O", "CO2"):
        assert light.bm[s] == ref.bm[s] and light.cNG[s] == ref.cNG[s]
        assert light.cRG[s].tobytes() == ref.cRG[s].tobytes()
    assert bk.ModelSettings.from_model(light) == bk.ModelSettings.from_model(ref)
    assert host_model.is_header_model(light) and not host_model.is_header_model(ref)
    assert light.PhiCO2 is None and light.face.size == 0


def test_record_methods_refuse(light):
    with pytest.raises(UnsupportedConfiguration):
        light.evaluate(np.zeros(25), (200.0, 150.0, 100.0), 1.0)
    with pytest.raises(UnsupportedConfiguration):
        light.seen_faces()
    with pytest.raises(UnsupportedConfiguration):
        light.csr(np.zeros(0))


@pytest.mark.parametrize("config", contract.CONFIGS, ids=[c["name"] for c in contract.CONFIGS])
def test_problem_on_header_model_is_the_same_problem(cpu_model, light, config):
    P_ref = contract.make_problem(cpu_model, config)
    env = dict(cpu_model, model=light)
    P = contract.make_problem(env, config)
    assert P.names == P_ref.names
    for name in ("x0", "lo", "hi", "Lt", "wNG", "wRG", "useNG", "useRG", "enhancement_basis"):
        if hasattr(P_ref, name):
            assert np.asarray(getattr(P, name)).tobytes() == \
                np.asarray(getattr(P_ref, name)).tobytes(), name
    assert bk.problem_contract(P) == bk.problem_contract(P_ref)


def test_backend_from_header_model_is_bitwise(cpu_model, light, clean_pool):
    from cops_cuda import evaluator
    from cops_cuda.memory import DeviceBudget
    run = cpu_model["synthetic"]["run"]
    budget = DeviceBudget(0)            # one pool owner for both backends
    for config in contract.CONFIGS:
        P_ref = contract.make_problem(cpu_model, config)
        P = contract.make_problem(dict(cpu_model, model=light), config)
        request = bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True)
        ref = bk.ReferenceBackend.from_problem(P_ref)
        for mode in ("sequential", "tree"):
            with evaluator.from_problem(P_ref, run, mode=mode, budget=budget) as full, \
                    evaluator.from_problem(P, run, mode=mode, budget=budget) as lean:
                assert lean.constants.PhiCO2_cpu is None
                assert lean.PhiCO2.tobytes() == full.PhiCO2.tobytes()
                for seed in (1, 2):
                    state = bk.ParameterState.from_problem(P, contract.perturbed(P, seed))
                    a, b = full.evaluate(state, request), lean.evaluate(state, request)
                    assert a.pNG.tobytes() == b.pNG.tobytes()
                    assert a.pRG.tobytes() == b.pRG.tobytes()
                    for key in request.keys():
                        for u, v in zip(a.columns[key], b.columns[key]):
                            assert u.tobytes() == v.tobytes(), key
                    if mode == "sequential":
                        want = ref.evaluate(state, request)
                        assert b.pNG.tobytes() == want.pNG.tobytes()
                if mode == "tree":
                    # The unmodified Problem through the facade, on the header model.
                    bk.install(P, lean)
                    bk.install(P_ref, full)
                    x = contract.perturbed(P, 3)
                    P.lam = P_ref.lam = 0.0
                    assert P.data_residuals(x).tobytes() == P_ref.data_residuals(x).tobytes()
                    assert P.data_jacobian(x).tobytes() == P_ref.data_jacobian(x).tobytes()
                    P.M, P_ref.M = light, cpu_model["model"]
    budget.close()

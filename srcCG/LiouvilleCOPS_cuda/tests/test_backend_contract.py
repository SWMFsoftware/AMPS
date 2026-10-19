"""C02 backend contract: exact column order and drop-in use by fit_cops.Problem.

Runs on the CPU (both interpreters): the contract is proven with the
reference backend so that C06 only has to show CUDA == reference through the
same interface.
"""

import argparse
import json

import numpy as np
import pytest

from cops_cuda import backend as bk
from cops_cuda import reference
from cops_cuda.errors import IdentityError, UnsupportedConfiguration

CONFIGS = (
    dict(name="legacy_d4_free_T_kappa", degree=4, enhancement_degree=4,
         fix_legacy_activity=False, fix_temperature=False, fix_kappa=False,
         ordered_temperature=False, templates=False),
    dict(name="legacy_d3_ordered_T", degree=3, enhancement_degree=4,
         fix_legacy_activity=False, fix_temperature=False, fix_kappa=False,
         ordered_temperature=True, templates=False),
    dict(name="legacy_d2_fixed_T", degree=2, enhancement_degree=4,
         fix_legacy_activity=False, fix_temperature=True, fix_kappa=False,
         ordered_temperature=False, templates=False),
    dict(name="enhancement_l4_control", degree=4, enhancement_degree=4,
         fix_legacy_activity=True, fix_temperature=True, fix_kappa=True,
         ordered_temperature=False, templates=False),
    dict(name="enhancement_l5_templates", degree=4, enhancement_degree=5,
         fix_legacy_activity=True, fix_temperature=True, fix_kappa=True,
         ordered_temperature=False, templates=True),
    # AMMENDED_ROADMAP primary family: B29 is legacy_d4_free_T_kappa above;
    # NC31 = B29 + two nonnegative template amplitudes, nothing frozen, one clip.
    dict(name="joint_nc31_free_baseline_templates", degree=4, enhancement_degree=4,
         fix_legacy_activity=False, fix_temperature=False, fix_kappa=False,
         ordered_temperature=False, templates=True),
)


@pytest.fixture
def cpu_model(synthetic):
    fit_cops = reference.load("fit_cops")
    run = str(synthetic["run"])
    header, names, rows = fit_cops.read_model_file(run + "/COPS_model.dat")
    cols = {name: index for index, name in enumerate(names)}
    elements = fit_cops.read_elements(run + "/elements.bin")
    rng = np.random.default_rng(3)
    coef = np.zeros((3, 25))
    coef[0] = rng.normal(0.0, 0.4, 25)
    coef[0][0] = 0.8                  # mostly positive activity, some clipping
    coef[1:] = 0.5
    model = fit_cops.Model(run, header, cols, rows, elements, coef)
    model.load_kernels(model.sel)
    return {"fit_cops": fit_cops, "model": model, "elements": elements, "coef": coef,
            "header": header, "synthetic": synthetic}


def make_problem(env, config):
    fit_cops = env["fit_cops"]
    args = argparse.Namespace(
        degree=config["degree"], enhancement_degree=config["enhancement_degree"],
        fix_legacy_activity=config["fix_legacy_activity"],
        fix_temperature=config["fix_temperature"], fix_kappa=config["fix_kappa"],
        ordered_temperature=config["ordered_temperature"], tbounds=(30.0, 400.0),
        enhancement_scale=None, max_nfev=50, ftol=1e-10, xtol=1e-10, gtol=1e-10, jobs=1)
    if config["templates"]:
        centroid = env["elements"]["centroid"]
        values = np.column_stack(((centroid[:, 2] > 0).astype(float),
                                  (centroid[:, 0] > 1000.0).astype(float)))
        args.activity_template_names = ["syn_north", "syn_east"]
        args.activity_template_sources = [{"name": n, "synthetic": True}
                                          for n in args.activity_template_names]
        args.activity_template_basis = values
    return fit_cops.Problem(env["model"], env["elements"], env["coef"], args)


def perturbed(problem, seed):
    rng = np.random.default_rng(seed)
    x = problem.x0.copy()
    if x.size:
        x = x + 0.05 * np.maximum(np.abs(x), 1.0) * rng.normal(size=x.size)
        x = np.clip(x, problem.lo, problem.hi)
    return x


def test_settings_parse_like_model(cpu_model):
    assert bk.ModelSettings.from_header(cpu_model["header"]) == \
        bk.ModelSettings.from_model(cpu_model["model"])


def test_settings_refuse_unsupported_modes(cpu_model):
    header = [h.replace("TemperatureModel=analytical", "TemperatureModel=table")
              for h in cpu_model["header"]]
    with pytest.raises(UnsupportedConfiguration):
        bk.ModelSettings.from_header(header)


def test_column_order_equals_reference_insertion_order(cpu_model):
    M = cpu_model["model"]
    request = bk.ColumnRequest([3, 0, 7], (), temperature=True, kappa=True)
    _, _, J = M.evaluate(cpu_model["coef"][0], M.T0, 1.0, jac=True, basis_index=[3, 0, 7])
    assert list(J.keys()) == request.keys()


@pytest.mark.parametrize("config", CONFIGS, ids=[c["name"] for c in CONFIGS])
def test_reference_backend_is_model_evaluate(cpu_model, config):
    P = make_problem(cpu_model, config)
    backend = bk.ReferenceBackend(P.M, P.enhancement_basis if P.nb else None,
                                  P.enhancement_scale)
    x = perturbed(P, 11)
    state = bk.ParameterState.from_problem(P, x)
    request = bk.ColumnRequest.from_problem(P)
    result = backend.evaluate(state, request)
    a, b, T, kappa = P.unpack(x)
    kwargs = ({"b": b, "enhancement_basis": P.enhancement_basis,
               "enhancement_scale": P.enhancement_scale} if P.nb else {})
    pNG, pRG, J = P.M.evaluate(a, T, kappa, jac=True, basis_index=P.fit_j,
                               enhancement_index=list(range(P.nb)), **kwargs)
    assert result.pNG.tobytes() == pNG.tobytes() and result.pRG.tobytes() == pRG.tobytes()
    assert list(result.columns) == request.keys()
    for key in request.keys():
        assert result.columns[key][0].tobytes() == J[key][0].tobytes()
        assert result.columns[key][1].tobytes() == J[key][1].tobytes()
    assert result.matrix("NG").shape == (len(P.M.rows), len(request.keys()))


@pytest.mark.parametrize("config", CONFIGS, ids=[c["name"] for c in CONFIGS])
def test_problem_through_facade_is_bitwise_identical(cpu_model, config):
    P = make_problem(cpu_model, config)
    x = perturbed(P, 5)
    r0, J0 = P.data_residuals(x), P.data_jacobian(x)
    lam = 1e-3
    P.lam = lam
    full_r0, full_J0 = P.residuals(x), P.jacobian(x)
    backend = bk.ReferenceBackend(P.M, P.enhancement_basis if P.nb else None,
                                  P.enhancement_scale)
    facade = bk.install(P, backend)
    assert P.M is facade
    r1, J1 = P.data_residuals(x), P.data_jacobian(x)
    assert r1.tobytes() == r0.tobytes()
    assert J1.shape == J0.shape and J1.tobytes() == J0.tobytes()
    assert P.residuals(x).tobytes() == full_r0.tobytes()
    assert P.jacobian(x).tobytes() == full_J0.tobytes()
    assert facade.calls["pressures"] == 2 and facade.calls["jacobian"] == 2


def test_facade_rejects_foreign_basis_and_scale(cpu_model):
    P = make_problem(cpu_model, CONFIGS[4])
    backend = bk.ReferenceBackend(P.M, P.enhancement_basis, P.enhancement_scale)
    facade = bk.install(P, backend)
    a, b, T, kappa = P.unpack(perturbed(P, 2))
    other = P.enhancement_basis.copy()
    facade.evaluate(a, T, kappa, b=b, enhancement_basis=other,
                    enhancement_scale=P.enhancement_scale)       # identical bytes: accepted
    other[0, 0] += 1.0
    with pytest.raises(IdentityError, match="not the one bound"):
        facade.evaluate(a, T, kappa, b=b, enhancement_basis=other,
                        enhancement_scale=P.enhancement_scale)
    with pytest.raises(IdentityError, match="scale"):
        facade.evaluate(a, T, kappa, b=b, enhancement_basis=P.enhancement_basis,
                        enhancement_scale=P.enhancement_scale * 2)


def test_facade_refuses_reordered_rows(cpu_model):
    P = make_problem(cpu_model, CONFIGS[0])
    backend = bk.ReferenceBackend(P.M)
    backend.rows = tuple(reversed(backend.rows))
    with pytest.raises(IdentityError, match="row order"):
        bk.install(P, backend)


class LimitedBackend(bk.ReferenceBackend):
    name = "limited"

    def capabilities(self):
        caps = super().capabilities()
        caps["column_kinds"] = {"a", "kappa"}
        caps["rg_models"] = {"calibration_srccg"}
        return caps


def test_unsupported_configuration_is_refused_before_evaluation(cpu_model):
    P = make_problem(cpu_model, CONFIGS[0])
    backend = LimitedBackend(P.M)
    state = bk.ParameterState.from_problem(P, P.x0)
    with pytest.raises(UnsupportedConfiguration) as excinfo:
        backend.evaluate(state, bk.ColumnRequest.from_problem(P))
    message = str(excinfo.value)
    assert "rg_model=equilibrium_sphere" in message and "['T']" in message


def test_parameter_state_validation(cpu_model):
    with pytest.raises(ValueError):
        bk.ParameterState(np.zeros(24), [1, 2, 3], 1.0)
    with pytest.raises(ValueError):
        bk.ParameterState(np.zeros(25), [1, np.nan, 3], 1.0)
    state = bk.ParameterState(np.zeros(25), [150, 140, 133], 1.0, b=np.zeros(4))
    assert not state.a.flags.writeable
    backend = bk.ReferenceBackend(cpu_model["model"])
    with pytest.raises(ValueError, match="enhancement coefficients"):
        backend.evaluate(state)
    assert state.digest() != bk.ParameterState(np.zeros(25), [150, 140, 133], 1.0).digest()


# ------------------------------------------------- C02 joint-mode extensions
NC31 = CONFIGS[-1]


def test_parameter_state_and_request_round_trip_exactly(cpu_model):
    P = make_problem(cpu_model, NC31)
    state = bk.ParameterState.from_problem(P, perturbed(P, 9))
    again = bk.ParameterState.from_dict(json.loads(json.dumps(state.to_dict())))
    assert again.digest() == state.digest()
    assert again.a.tobytes() == state.a.tobytes() and again.b.tobytes() == state.b.tobytes()
    request = bk.ColumnRequest.from_problem(P)
    assert bk.ColumnRequest.from_dict(json.loads(json.dumps(request.as_dict()))) == request
    assert len(request.keys()) == 31


def test_problem_contract_identifies_templates_and_prior(cpu_model):
    P = make_problem(cpu_model, NC31)
    contract = bk.problem_contract(P)
    assert contract["names"][25:27] == ["c_syn_north", "c_syn_east"]
    assert contract["prior"]["reg_slice"] == [0, 27]
    assert contract == bk.problem_contract(make_problem(cpu_model, NC31))
    # Reordering the templates changes the contract (column order is identity).
    swapped = dict(NC31)
    P2 = make_problem(cpu_model, swapped)
    P2_basis = P2.enhancement_basis[:, ::-1].copy()
    assert bk.array_digest(P2_basis) != contract["enhancement_basis_sha256"]
    # The prior is part of the identity: a different run prior is detected.
    B29 = make_problem(cpu_model, CONFIGS[0])
    assert bk.problem_contract(B29)["prior"]["coefficients_sha256"] == contract["prior"]["coefficients_sha256"]
    assert bk.problem_contract(B29)["sha256"] != contract["sha256"]


def test_install_refuses_stale_or_reordered_templates(cpu_model):
    P = make_problem(cpu_model, NC31)
    reordered = bk.ReferenceBackend(P.M, P.enhancement_basis[:, ::-1].copy(), P.enhancement_scale)
    with pytest.raises(IdentityError, match="stale or reordered"):
        bk.install(P, reordered)
    rescaled = bk.ReferenceBackend(P.M, P.enhancement_basis, P.enhancement_scale * 1.5)
    with pytest.raises(IdentityError, match="scale"):
        bk.install(P, rescaled)
    one_template = bk.ReferenceBackend(P.M, P.enhancement_basis[:, :1].copy(), P.enhancement_scale)
    with pytest.raises(IdentityError, match="enhancement columns"):
        bk.install(P, one_template)
    facade = bk.install(P, bk.ReferenceBackend.from_problem(P))
    assert facade._request.keys()[-1] == ("kappa", 0)

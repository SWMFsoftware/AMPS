"""C12 fit restart files (cops_cuda/restart.py) on a synthetic run.

* write/read round trip, with corruption detection;
* the optimizer start rebuilt from a restart equals the fitted x bitwise, for the same
  problem; for a problem in the ordered-temperature parameterization the same physical
  temperatures are reproduced; names missing from the restart keep the problem's start;
* GPU: the restart state's pressures (StateEvaluator) equal the CPU reference within the
  approved pressure budget.
"""

import json
import shutil

import numpy as np
import pytest

from cops_cuda import restart as rs
import test_backend_contract as contract

cpu_model = contract.cpu_model
EPS = np.finfo(float).eps


def template_config():
    return next(c for c in contract.CONFIGS if c["templates"])


def make_restart(env, tmp_path, config, seed=5, at_prior=False):
    P = contract.make_problem(env, config)
    x = P.x0.copy() if at_prior else contract.perturbed(P, seed)
    run = env["synthetic"]["run"]
    fit_dir = tmp_path / "fit"
    fit_dir.mkdir()
    shutil.copy(run / "COPS_model.dat", fit_dir / "COPS_model.dat")     # same table format
    (fit_dir / "fit_manifest.json").write_text(json.dumps({"inputs": {}, "selection": {"sha256": "x"}}))
    result = {"x": x, "lam": 0.0, "message": "test", "status": 0, "nfev": 3, "optimality": 1.0,
              "cost": 1.0, "success": False}
    r = rs.from_fit(P, result, fit_dir, {"inputs": {}}, env["elements"])
    return P, x, r


def test_round_trip_and_corruption(cpu_model, tmp_path):
    P, x, r = make_restart(cpu_model, tmp_path, template_config())
    path = tmp_path / "fit_restart.npz"
    rs.write(path, r)
    back = rs.read(path)
    for k in rs.ARRAY_KEYS:
        assert np.array_equal(back.arrays[k], r.arrays[k]), k
    assert back.template_names == P.template_names
    assert back.mesh_identity() == r.mesh_identity()
    with pytest.raises(FileExistsError):
        rs.write(path, r)
    # Corrupt one array: the stored digest catches it.
    bad = r.copy()
    bad.arrays["a"] = bad.arrays["a"].copy()
    path2 = tmp_path / "bad.npz"
    rs.write(path2, bad)
    with np.load(str(path2)) as z:
        data = {k: np.array(z[k]) for k in z.files}
    data["a"][0] += 1.0
    path3 = tmp_path / "bad3.npz"
    with open(path3, "xb") as f:
        np.savez_compressed(f, **data)
    with pytest.raises(ValueError, match="corrupt"):
        rs.read(path3)


def test_start_vector_reproduces_the_fit_state(cpu_model, tmp_path):
    for config in contract.CONFIGS:
        sub = tmp_path / config["name"]
        sub.mkdir()
        P, x, r = make_restart(cpu_model, sub, config)
        xs, report = rs.start_vector(P, r)
        assert xs.tobytes() == x.tobytes(), config["name"]
        assert set(report["from_restart"]) == set(P.names)
        assert not report["kept_problem_start"]


def test_start_vector_across_parameterizations(cpu_model, tmp_path):
    base = dict(contract.CONFIGS[0])            # free T; the run's T0 is ordered (Tmax >= T75 >= Tmin)
    P, x, r = make_restart(cpu_model, tmp_path, base, at_prior=True)
    ordered = dict(base, ordered_temperature=True)
    P2 = contract.make_problem(cpu_model, ordered)
    xs, report = rs.start_vector(P2, r)
    _, _, T2, kappa2 = P2.unpack(xs)
    _, _, T1, kappa1 = P.unpack(x)
    np.testing.assert_allclose(T2, T1, rtol=1e-14)
    assert kappa2 == kappa1
    # A problem with templates the restart lacks keeps its own start for them.
    P3 = contract.make_problem(cpu_model, template_config())
    xs3, rep3 = rs.start_vector(P3, r)
    assert all(n.startswith("c_") or n.startswith("b_") for n in rep3["kept_problem_start"])


cp = pytest.importorskip("cupy")


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() == 0, reason="no CUDA device")
def test_state_pressures_on_the_gpu(cpu_model, tmp_path, clean_pool):
    from cops_cuda import backend as bk
    config = template_config()
    P, x, r = make_restart(cpu_model, tmp_path, config)
    run = cpu_model["synthetic"]["run"]
    ev = rs.StateEvaluator(run, 0)
    try:
        png, prg = ev.pressures(r)
    finally:
        ev.release()
    # CPU reference: the full-row model, the same state.
    fit_cops = cpu_model["fit_cops"]
    header, names, rows = fit_cops.read_model_file(str(run / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    M = fit_cops.Model(str(run), header, cols, rows, cpu_model["elements"], cpu_model["coef"])
    basis = rs.enhancement_basis(r)
    b = rs.b_vector(r)
    rN, rR = M.evaluate(r.a, r.T, float(r.kappa[0]), b=b if b.size else None,
                        enhancement_basis=basis, enhancement_scale=r.enhancement_scale)
    scale = np.maximum(np.abs(rN), 1e-300)
    assert np.all(np.abs(png - rN) <= 1e-12 * scale + 1e-300)
    assert np.all(np.abs(prg - rR) <= 1e-12 * np.maximum(np.abs(rR), 1e-300) + 1e-300)
    del bk


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() == 0, reason="no CUDA device")
def test_point_changes_are_the_exact_pressure_difference(cpu_model, tmp_path, clean_pool):
    """The pressures are linear in each facet's clipped activity: the summed per-record changes
    (kernels/restart.cu) equal the difference of two full evaluations."""
    config = template_config()
    P, x, r = make_restart(cpu_model, tmp_path, config)
    run = cpu_model["synthetic"]["run"]
    act = np.asarray(r.activity_linear)
    region = np.zeros(len(act), bool)
    region[::7] = True
    column = region.astype(float)
    delta_af = np.maximum(act + 0.3 * r.enhancement_scale * column, 0.0) - np.maximum(act, 0.0)
    ev = rs.StateEvaluator(run, 0)
    try:
        share_ng, share_rg, p0n, p0r = ev.facet_point_shares(r, delta_af, ev.point_ids)
        dng, drg = ev.point_changes(r, delta_af)
        basis = np.column_stack((rs.enhancement_basis(r), column))
        b1 = np.concatenate((rs.b_vector(r), [0.3]))
        p1n, p1r = ev.pressures(r, basis=basis, b=b1)
    finally:
        ev.release()
    scale = np.maximum(np.abs(p0n), 1e-300)
    assert np.all(np.abs((p1n - p0n) - dng) <= 1e-12 * scale)
    assert np.all(np.abs((p1r - p0r) - drg) <= 1e-12 * np.maximum(np.abs(p0r), 1e-300))
    assert np.all(share_ng[~region] == 0.0) and np.any(share_ng[region] > 0.0)

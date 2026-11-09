#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_positive_fit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_positive_fit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 7809a9cfef80dafd05dc4355c2d52716cf1f26eab7b8a776ebfced71cb9ee59c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. CPU tests of positive_problem.py (projection, Tikhonov incl. the held term, W10x5, grid rows, Jacobian,
#   start restart, resume refusal with another ratio or held term). user decisions 2026-11-07/08.
# USED BY: jobs/positive_unit_12m.sbatch.
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
"""CPU tests of model_ext/positive_problem.py (positive_fit.py's fit-level part; 2026-11-07).

    .venv/bin/python test_positive_fit.py [--scratch DIR]       light tests (login node: < 1 min, < 1 GB)
    .venv/bin/python test_positive_fit.py --12m [--out JSON]    start projection + Tikhonov on the 12M base (batch job)

Light tests:
 1. W10x5 operator == late_response_eval.responses() on synthetic series (random grid offset, random stops), and on the
    122289 forward run (runs/shap8_012m_flat_earlier/late_response/122289/eval/series_dz1.5.npz: W @ grid series ==
    stored W10x5 of every model, NG and RG);
 2. start projection (emission-weighted since 2026-11-07: exact recovery of a synthetic exponent field, the
    negative-cliff rule) and the log Tikhonov
    operator (Lt^T Lt == G, D(x0) = 0, D == area-weighted rms of s - s0);
 3. argument rewriting; the selection wrapper with a real cops_selection profile on synthetic rows;
 4. end to end on a synthetic base with the UNCHANGED fit_cops.main(): all layers of positive_fit.py (selection,
    exponent Problem, grain + crossing, W10x5 on M and MF, grain start, writers, Jacobian check) with a CPU stand-in
    for agent B's GPU patch (ExpModel: af = exp(act_lin), reference night ratio, no floor), lambda 0 and 0.3, each in
    a subprocess; checks: Jacobian check, D, COPS_model.dat late rows == W10x5(surface) + grain, the patched activity
    products, the manifest record.
 5. (2026-11-08, CO2 flux floor) the variable handling: set_positive_environment sets COPS_POS_FLOOR_CO2=1 for the
    ratio 1e-5 only and refuses it (any value) for 0.02 / 0.002; positive_activity._read_environment: the canonical
    configuration of every earlier setting equals the 2026-11-07 module's (frozen c1 snapshot), refusals; the restart
    marker's co2_floor and the resume refusal; the host assembly of positive_activity's CO2-floor evaluation
    (_evaluate_co2_floor) with an emulation of the documented kernel sums against an independent NumPy model of
    F_CO2 = max(kappa Q ACO2, r qmax ACO2) (pressures; kappa, T and activity columns vs the analytic kappa derivative
    and central differences; kappa ~ 0: kappa column exactly 0; kappa = 2: bitwise the original _evaluate); e2e:
    --ratio 1e-5 from the lambda-0.3 restart (start check reports the CO2-floor change; marker, manifest, text, npz)
    and its resume (reproduces).
 6. (2026-11-08, --start-set) rewrite_arguments: absent option -> the configuration of before (no key), refusals
    (without --start-restart, bad syntax, duplicate, grain names with --grain-start-uPa); e2e: --start-restart +
    --start-set kappa=1 (start x kappa = 1, the rest the source's; start check "NOT the source's state"; manifest,
    marker, fit_chain, ACTIVITY_PARAMETERIZATION.txt), its resume (chain keeps start_set, reproduces), an unknown name
    and a value outside the bounds refused; the fits without the option record no start_set anywhere.
"""

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
REF = APP.parent / "LiouvilleCOPS"
for p in (str(HERE), str(APP / "tools"), str(APP)):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("COPS_OUTPUT_ROOT", "/data001/highresolu_vt/cops_output")
SERIES = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_flat_earlier/late_response/122289/eval/"
              "series_dz1.5.npz")
RESULTS = []


def report(name, ok, detail):
    RESULTS.append((name, ok, detail))
    print("%-58s %s  %s" % (name, "PASS" if ok else "FAIL", detail), flush=True)


# ------------------------------------------------------------------------------------------------ 1 window
def test_window():
    import positive_problem as pp
    sys.path.insert(0, str(HERE))
    import late_response_eval as lre
    rng = np.random.default_rng(1)
    worst = 0.0
    for trial in range(5):
        t0 = 528503600.0 + rng.uniform(0, 1)
        tg = t0 + np.arange(200, dtype=float) + (rng.uniform(-0.05, 0.05, 200) if trial == 4 else 0.0)
        stops = np.sort(rng.uniform(tg[0] + 18.5, tg[-1] - 0.01, 40))
        stops[-1] = tg[-1]                                   # the last observation ends on the last grid point
        W = pp.w10x5_operator(tg, stops)
        for _ in range(3):
            y = np.exp(rng.normal(-14, 1.0, len(tg)))
            ref = np.array([lre.responses(tg, y, T)["W10x5"] for T in stops])
            worst = max(worst, float(np.max(np.abs(W @ y - ref) / np.abs(ref))))
    report("W10x5 operator == late_response_eval (synthetic)", worst < 1e-13, "max rel diff %.2e" % worst)
    # the copies in positive_problem equal the originals
    tg = 100.0 + np.arange(50.0)
    y = rng.normal(size=50)
    d = max(abs(pp.responses(tg, y, 140.3)[k] - lre.responses(tg, y, 140.3)[k]) for k in ("W4", "W10", "W10x5"))
    report("positive_problem.responses == late_response_eval.responses", d == 0.0, "diff %.1e" % d)
    if not SERIES.is_file():
        report("W10x5 on the 122289 series", False, "missing %s" % SERIES)
        return
    z = np.load(SERIES)
    grid = np.nonzero(z["kind"] == "grid")[0]
    tg, stops = z["et"][grid], z["stop_et"]
    W = pp.w10x5_operator(tg, stops)
    worst, n = 0.0, 0
    for key in z.files:
        if key.endswith("|W10x5"):
            series = z[key.rsplit("|", 1)[0]]
            worst = max(worst, float(np.max(np.abs(W @ series[grid] - z[key]) / np.abs(z[key]))))
            n += 1
    report("W10x5 operator on 122289 series_dz1.5 (%d series)" % n, worst < 1e-12 and n > 0,
           "max rel diff %.2e; W %s, nnz %d" % (worst, W.shape, W.nnz))
    # the layout class on the 122289 point table (late-only table, base rows = table rows)
    table = pp.read_point_table(SERIES.parent.parent / "late_grid_points.txt")
    points = {i: (str(u), float(e), 0.0, 0.0) for i, (u, e) in enumerate(zip(z["utc"], z["et"]))}
    lay = pp.LateLayout(table, points, np.arange(len(points)))
    d = float(np.max(np.abs(lay.stop_et - stops)))
    report("LateLayout on the 122289 table (stops, grid)", d < 1e-6 and np.array_equal(lay.grid_rows, grid),
           "max |stop - stored| %.1e s; %d late, %d grid" % (d, len(lay.late_rows), len(lay.grid_rows)))


# ------------------------------------------------------------------------------------------------ 2 start, Tikhonov
def test_projection():
    import positive_problem as pp
    rng = np.random.default_rng(2)
    n = 20000
    Y = rng.normal(size=(n, 25))
    Y[:, 0] = 1.0
    area = rng.uniform(0.5, 2.0, n)
    T = (rng.random((n, 1)) < 0.15).astype(float)
    a_true = rng.normal(0, 0.3, 25)
    a_true[0] = 40.0
    s_true = Y @ a_true + 0.7 * T[:, 0]
    A = np.exp(s_true)
    A[rng.random(n) < 0.01] = 0.0                                  # zero facets are excluded
    pr = pp.start_projection(Y, area, A, T)
    err = max(float(np.max(np.abs(pr["a"] - a_true))), abs(float(pr["c"][0]) - 0.7))
    st = pr["stats"]
    report("projection recovers an exact exponent field", err < 1e-9 and st["rms_ln_residual_over_A_orig_gt_0"] < 1e-9,
           "max coef err %.1e, rms %.1e, zero area %.4f" % (err, st["rms_ln_residual_over_A_orig_gt_0"],
                                                           st["area_fraction_A_orig_zero"]))
    A2 = np.exp(Y @ a_true - 0.5 * T[:, 0] + rng.normal(0, 0.2, n))
    pr2 = pp.start_projection(Y, area, A2, T)
    w = area * A2                                                  # emission weights (user 2026-11-07)
    a_y = np.linalg.solve(Y.T @ (w[:, None] * Y), Y.T @ (w * np.log(A2)))
    ok = (pr2["stats"]["template_coefficient_clamped_to_0_and_Y_reprojected"] and pr2["c"][0] == 0.0 and
          np.allclose(pr2["a"], a_y, rtol=0, atol=1e-10) and pr2["stats"]["projected_template_coefficients"][0] < 0)
    report("negative cliff coefficient -> 0, Y re-projected", ok,
           "projected c %.4f; rms %.4f" % (pr2["stats"]["projected_template_coefficients"][0],
                                          pr2["stats"]["rms_ln_residual_over_A_orig_gt_0"]))
    ascale = pr["ascale"]
    Lt, G, info = pp.log_tikhonov(pr["G_all"], float(area.sum()), ascale, 25)
    x0 = np.concatenate([pr["a"] / ascale, pr["c"]])
    x = x0 + rng.normal(0, 1e-2, 26)
    s0 = Y @ (x0[:25] * ascale) + T[:, 0] * x0[25]
    s = Y @ (x[:25] * ascale) + T[:, 0] * x[25]
    D_direct = math.sqrt(float(np.sum(area * (s - s0) ** 2) / area.sum()))
    D_op = float(np.linalg.norm(Lt @ (x - x0)))
    ok = info["LtT_Lt_minus_G_rel"] < 1e-12 and np.linalg.norm(Lt @ (x0 - x0)) == 0 and abs(D_op / D_direct - 1) < 1e-10
    report("log Tikhonov: Lt^T Lt == G, D(x0) = 0, D == rms(s - s0)", ok,
           "Lt^T Lt rel %.1e; D op/direct - 1 = %.1e" % (info["LtT_Lt_minus_G_rel"], D_op / D_direct - 1))


# ------------------------------------------------------------------------------------------------ 3 arguments, selection
def test_arguments(scratch):
    import positive_problem as pp
    table = Path(scratch) / "pts.txt"
    table.write_text("# index UTC kind original_index L2_start L2_stop\n0 2016-09-30T10:33:10 obs 1 a b\n")
    cfg, argv, notes = pp.rewrite_arguments(["--point-table", str(table), "--check-jacobian", "--devices", "0,1",
                                             "--", "/r/base", "--profile", "all_available", "--lambda", "0",
                                             "--check-tolerance", "1e-6", "--activity-template", "B", "cliff_ge_60deg"])
    ok = (cfg["check_jacobian"] and cfg["grain_start_uPa"] is None and cfg["grain_start_B"] is None and
          cfg["start_restart"] is None and cfg["ratio"] is None and
          argv[:4] == ["--devices", "0,1", "--template-start", "0"] and "--fix-temperature" in argv and
          argv[argv.index("--check-tolerance") + 1] == "inf" and "--check-jacobian" not in argv and cfg["run"] == "/r/base")
    report("rewrite_arguments: inserts/replaces", ok, " | ".join(notes))
    tmpl = ["--activity-template", "B", "cliff_ge_60deg"]
    cases = ((["--template-start", "1"], tmpl), ([], ["--fix-kappa"] + tmpl), ([], ["--degree", "3"] + tmpl),
             ([], ["--ordered-temperature"] + tmpl), ([], ["--check-jacobian"] + tmpl), ([], []),
             ([], ["--activity-template", "B", "cliff_ge_50deg"]), ([], tmpl + ["--activity-template", "B", "ft_c1"]))
    refused = 0
    for head, tail in cases:
        try:
            pp.rewrite_arguments(["--point-table", str(table)] + head + ["--", "x"] + tail)
        except SystemExit:
            refused += 1
    try:
        pp.rewrite_arguments(["--point-table", str(table), "--", "x"] + tmpl)
        accepted = True
    except SystemExit:
        accepted = False
    report("rewrite_arguments: refusals (incl. template != one cliff_ge_60deg)", refused == len(cases) and accepted,
           "%d of %d refused; the valid command accepted: %s" % (refused, len(cases), accepted))


def test_selection():
    import positive_problem as pp
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    cols = {k: i for i, k in enumerate(("iPoint", "UTC", "inside", "refNG", "refRG", "useNG", "useRG"))}
    rows, gkeys = [], set()
    for k in range(10):                                      # 6 observations, 4 grid rows
        grid = k in (2, 5, 6, 9)
        u = "2016-09-30T10:%02d:%02d%s" % (30 + k // 6, (k * 7) % 60, ".282" if grid else "")
        rows.append([str(k), u, "0", "0" if grid else "1e-9", "0" if grid else "2e-9", "0", "0"])
        if grid:
            gkeys.add(pp.time_key(u))
    profile = {"schema": "liouville_cops.selection_profile", "schema_version": 1, "name": "t",
               "candidate_policy": "positive_outside_shape", "apply_fit_interval": False,
               "apply_run_exclusions": False, "exclude_intervals": {"NG": [], "RG": []}, "saturated_ng": "include",
               "expected_counts": {"rows": 6, "used_ng": 6, "used_rg": 6}}
    rec = {}
    original = fit_cops.apply_selection_profile
    try:
        pp.install_selection(fit_cops, gkeys, rec)
        res, report_ = fit_cops.apply_selection_profile(rows, cols, profile, {"FitInterval": []}, [False] * 10)
    finally:
        fit_cops.apply_selection_profile = original
    use = [(r[cols["useNG"]], r[cols["useRG"]]) for r in res]
    ok = (all(use[k] == ("0", "0") for k in (2, 5, 6, 9)) and all(use[k] == ("1", "1") for k in (0, 1, 3, 4, 7, 8))
          and report_["reasons"]["NG"][5] == pp.GRID_REASON and report_["reasons"]["RG"][0] == "selected"
          and report_["reason_counts"]["NG"] == {"late_grid_point": 4, "selected": 6}
          and report_["counts"]["rows"] == 6 and report_["counts"]["late_grid_rows"] == 4)
    report("selection wrapper: audited profile on observation rows only", ok, str(report_["reason_counts"]["NG"]))


# ------------------------------------------------------------------------------------------------ 4 end to end
HEADER_LINES = (
    "LiouvilleCOPS model of the COPS gauge pressures; synthetic test run (positive_fit tests)",
    "ProductKind=forward TableSchema=liouville_cops_forward_v1",
    "ShapeModel=synthetic FlightTimeRotation=none FitInterval=all",
    "AnalyticalTmax=182.09922 AnalyticalT75=136.11149 AnalyticalTmin=133 HeliocentricModel=shap5 ActivityFile=none",
    "SurfaceBC=flat TemperatureModel=analytical ActivityModel=shap5 SZAMode=continuous CastShadows=on "
    "SpacecraftVelocity=on",
    "NGAcceptance=srccg_shield NGHeadOnWeight=on RGAcceptance=hemisphere RGPressureModel=equilibrium_sphere "
    "TemperatureNG=293 TemperatureRGMode=hk_table SpeciesSensitivity=on (H2O 0.893, CO2 0.704) ReferenceDataUnit=mbar",
)
COLUMNS = ("iPoint", "UTC", "et", "x", "y", "z", "altitude", "inside", "rh_AU", "VdotB_RG", "nVisible", "Q_H2O",
           "Q_CO2", "n_H2O", "n_CO2", "PhiRG_H2O", "PhiRG_CO2", "pNG_H2O", "pNG_CO2", "pNG", "pRG_H2O", "pRG_CO2",
           "pRG", "refNG", "refRG", "NGBelowDeckFraction", "T_RG", "useNG", "useRG")
RECORD = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"), ("KV", "<f8"), ("cos", "<f8"), ("shadow", "<i4")])
ET0 = 528503598.182348            # et of 2016-09-30T10:32:10 (et - UTC offset as in the 12M base: t_end 10:33:10)


def exp_model_class(fit_cops):
    """CPU stand-in for agent B's GPU patch: the reference Model.evaluate with af = exp(act_lin) (dA/da = A Y);
    reference night ratio, no flux floor."""
    np_ = np

    class ExpModel(fit_cops.Model):
        def evaluate(self, a, T, kappa, b=None, enhancement_basis=None, enhancement_scale=1.0, jac=False,
                     basis_index=None, enhancement_index=None):
            import positive_problem as pp_
            Tmax, T75, Tmin = T
            act_lin = self.Y @ a
            if b is not None:
                act_lin = act_lin + enhancement_scale * (enhancement_basis @ b)
            if pp_.HELD["h"] is not None:              # as positive_problem.install_held_activity on the GPU host path
                pp_.HELD["calls"] += 1
                act_lin = act_lin + pp_.HELD["h"]
            with np_.errstate(over="ignore", invalid="ignore"):
                A = np_.exp(act_lin)
            af = A[self.face]
            Tl = Tmax + fit_cops.KT * (Tmax - T75) * self.xlit
            use_lit = self.lit & (Tl > Tmin)
            Te = np_.where(use_lit, Tl, Tmin)
            u = 1.0 / np_.sqrt(Te)
            bc = lambda w: np_.bincount(self.pt, weights=w, minlength=self.npt)     # noqa: E731
            cNG, cRG, bm = self.cNG, self.cRG, self.bm
            with np_.errstate(over="ignore", invalid="ignore"):
                pNG_H = cNG["H2O"] * bm["H2O"] * bc(self.wn["H2O"] * af * u)
                pRG_H = cRG["H2O"] * (bc(self.wp["H2O"] * af) + bm["H2O"] * bc(self.wv["H2O"] * af * u))
            pNG_C = cNG["CO2"] * bm["CO2"] * bc(self.wn["CO2"] * u)
            pRG_C = cRG["CO2"] * (self.PhiCO2 + bm["CO2"] * bc(self.wv["CO2"] * u))
            pNG, pRG = pNG_H + kappa * pNG_C, pRG_H + kappa * pRG_C
            if not jac:
                return pNG, pRG
            J = {}

            def columns(Bm):
                m = A[:, None] * Bm
                return (cNG["H2O"] * bm["H2O"] * (self.csr(self.wn["H2O"] * u) @ m),
                        cRG["H2O"][:, None] * (self.csr(self.wp["H2O"] + bm["H2O"] * self.wv["H2O"] * u) @ m))
            if basis_index:
                JN, JR = columns(self.Y[:, basis_index])
                for k, j in enumerate(basis_index):
                    J[("a", j)] = (JN[:, k], JR[:, k])
            if enhancement_index:
                JN, JR = columns(enhancement_scale * enhancement_basis[:, enhancement_index])
                for k, j in enumerate(enhancement_index):
                    J[("b", j)] = (JN[:, k], JR[:, k])
            du = -0.5 * u / Te
            dl, dn_ = np_.where(use_lit, du, 0.0), np_.where(use_lit, 0.0, du)
            dlx = dl * self.xlit

            def sums(w):
                s_l, s_lx, s_n = bc(w * dl), bc(w * dlx), bc(w * dn_)
                return fit_cops.KT * s_lx + s_l, -fit_cops.KT * s_lx, s_n
            gN = cNG["H2O"] * bm["H2O"] * self.wn["H2O"] * af + kappa * cNG["CO2"] * bm["CO2"] * self.wn["CO2"]
            tN, tRH, tRC = sums(gN), sums(bm["H2O"] * self.wv["H2O"] * af), sums(bm["CO2"] * self.wv["CO2"])
            for k in range(3):
                J[("T", k)] = (tN[k], cRG["H2O"] * tRH[k] + kappa * cRG["CO2"] * tRC[k])
            J[("kappa", 0)] = (pNG_C, pRG_C)
            return pNG, pRG, J
    return ExpModel


def utc_of(sec):
    """'2016-09-30Th:mm:ss[.fff]' for seconds after 2016-09-30T00:00."""
    h, r = divmod(sec, 3600.0)
    m, s = divmod(r, 60.0)
    frac = round(s - int(s), 3)
    return "2016-09-30T%d:%02d:%02d%s" % (int(h), int(m), int(s), (".%03d" % round(frac * 1000)) if frac else "")


def make_base(root, seed=5):
    """A synthetic base with observations around an NG gap and a 1-s late grid; refs from the exp model at the start."""
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    at = reference.load("activity_templates")
    import positive_problem as pp
    rng = np.random.default_rng(seed)
    root = Path(root)
    run = root / "base"
    (run / "kernels").mkdir(parents=True)
    nf = 400
    k = np.arange(nf) + 0.5
    pol, azi = np.arccos(1.0 - 2.0 * k / nf), np.pi * (1.0 + 5.0 ** 0.5) * k
    unit = np.column_stack((np.sin(pol) * np.cos(azi), np.sin(pol) * np.sin(azi), np.cos(pol)))
    coef = fit_cops.read_activity_coefficients(str(REF), list(HEADER_LINES))
    el = np.empty((nf, 34))
    el[:, 0:3], el[:, 3:6] = 2000.0 * unit, unit
    el[:, 6] = rng.uniform(1e3, 3e3, nf)
    el[:, 7:32] = rng.normal(0.0, 0.3, (nf, 25))
    el[:, 7] = 1.0
    el[:, 32] = np.maximum(el[:, 7:32] @ coef[0], 0.0)
    el[:, 33] = np.maximum(el[:, 7:32] @ coef[1], 0.0)
    with (run / "elements.bin").open("wb") as f:
        f.write(np.array([1, nf, 25], dtype="<i4").tobytes())
        f.write(np.ascontiguousarray(el, dtype="<f8").tobytes())
    tcol = (rng.random(nf) < 0.2).astype(float)
    gc, gt = reference.load("gravity_cache"), reference.load("gravity_topology")
    at.write_bundle(str(root / "templates.npz"), ["cliff_ge_60deg"], tcol[:, None], el[:, 0:3], el[:, 6],
                    {"templates": {"cliff_ge_60deg": {"definition": "synthetic test"}},
                     "gravity": {"cache": {"definition": {"algorithm": gc.ALGORITHM,
                                                          "orientation_policy": gt.ORIENTATION_POLICY}},
                                 "integrity_check": "DISABLE after topology audit"}})
    # points: observations 10:20:00..10:24:30 /30 s (NG+RG), 10:25..10:33 /60 s (RG only), 10:33:10..10:34:26 /4 s
    # (NG+RG, the late window); grid 10:32:46.282 + k s up to 10:34:30.282
    obs = [36000 + 1200 + 30 * i for i in range(10)] + [36000 + 1500 + 60 * i for i in range(9)] + \
          [36000 + 1990 + 4 * i for i in range(20)]
    grid = [36000 + 1966.282 + i for i in range(105)]
    pts = sorted([(t, "obs") for t in obs] + [(t, "grid") for t in grid])
    stops = {}
    rows, geo = [], []
    for j, (t, kind) in enumerate(pts):
        et = ET0 + (t - (36000 + 1930))
        sc = np.array([2600.0 + 3.0 * j, 400.0 - 2.0 * j, 1800.0 - 1.5 * j])
        rh = 3.831 * fit_cops.AU
        alt = float(np.linalg.norm(sc) - 2000.0)
        vb = 0.3 + 0.001 * j
        g = np.array([et, vb, sc[0], sc[1], sc[2], 0.6, 0.0, 0.8, rh, alt])
        cnt = int(rng.integers(150, 300))
        rec = np.zeros(cnt, dtype=RECORD)
        rec["face"] = rng.permutation(nf)[:cnt]
        rec["Kn"] = rng.uniform(0.0, 1e-6, cnt) * (1 + 0.002 * j)
        rec["KPhi"] = rng.uniform(0.0, 1e-6, cnt)
        rec["KV"] = rng.uniform(0.0, 1e-7, cnt)
        rec["cos"] = rng.uniform(-1.0, 1.0, cnt)
        rec["shadow"] = (rng.random(cnt) < 0.1).astype("<i4")
        hdr = np.zeros(1, dtype=[("version", "<i4"), ("point", "<i4"), ("geometry", "<f8", (10,)), ("count", "<i4")])
        hdr["version"], hdr["point"], hdr["geometry"], hdr["count"] = 1, j, g, cnt
        with (run / "kernels" / ("point=%04d.bin" % j)).open("wb") as f:
            f.write(hdr.tobytes())
            f.write(rec.tobytes())
        rows.append([str(j), utc_of(t), "%.6f" % et, "%.6e" % sc[0], "%.6e" % sc[1], "%.6e" % sc[2], "%.6e" % alt, "0",
                     "%.8f" % (rh / fit_cops.AU), "%.6e" % vb, str(cnt)] + ["nan"] * 8 + ["1e-7"] + ["nan"] * 2 +
                    ["1e-7", "0", "0", "nan", "250.082", "0", "0"])
        if kind == "obs" and t >= 36000 + 1990:
            stops[j] = utc_of(t + 0.224)
    write_table(run, rows)
    (run / "input.snapshot.in").write_text("ExcludeNG = 2016-09-30T10:24:40 2016-09-30T10:40:00\n")
    with (root / "late_grid_points.txt").open("w") as f:
        f.write("# index UTC kind original_index L2_start L2_stop\n")
        for j, (t, kind) in enumerate(pts):
            f.write("%d %s %s %d %s %s\n" % (j, utc_of(t), kind, j if kind == "obs" else -1,
                                             utc_of(t - 3.776) if j in stops else "-", stops.get(j, "-")))
    # reference data: the exp model at the projected start (+ a grain after the gap), x lognormal noise
    cols = {c: i for i, c in enumerate(COLUMNS)}
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    M = exp_model_class(fit_cops)(str(run), list(HEADER_LINES), cols, rows, elements, coef)
    proj = pp.start_projection(elements["Y"], elements["area"], elements["AH2O"], tcol[:, None])
    pNG, pRG = M.evaluate(proj["a"], M.T0, 1.0, b=proj["c"], enhancement_basis=tcol[:, None], enhancement_scale=1.0)
    for j, (t, kind) in enumerate(pts):
        if kind == "grid":
            continue
        noise = np.exp(rng.normal(0.0, 0.15, 2))
        ng = pNG[j] * noise[0] * 0.8 + (3e-1 * pNG[j] if t >= 36000 + 1990 else 0.0)
        rows[j][cols["refRG"]] = "%.6e" % (pRG[j] * noise[1])
        rows[j][cols["useRG"]] = "1"
        if not 36000 + 1480 <= t < 36000 + 1990:
            rows[j][cols["refNG"]] = "%.6e" % ng
            rows[j][cols["useNG"]] = "1"
    write_table(run, rows)
    return run, root / "templates.npz", root / "late_grid_points.txt"


def write_table(run, rows):
    with (Path(run) / "COPS_model.dat").open("w") as f:
        for line in HEADER_LINES:
            f.write("# %s\n" % line)
        f.write("#" + " ".join(COLUMNS) + "\n")
        for r in rows:
            f.write(" ".join(r) + "\n")


def e2e_child(run, bundle, points, out, lam, extra=(), max_nfev=4):
    """One fit_cops.main() with every layer of positive_fit.py (child process). ``extra``: more driver options (before
    '--'), e.g. --start-restart R, --grain-start-uPa A, --ratio R; --restart-from R is emulated with gpu_fit's
    install_template_start (the GPU driver's resume path)."""
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    import grain_crossing
    import positive_problem as pp
    extra = list(extra)
    ratio = extra[extra.index("--ratio") + 1] if "--ratio" in extra else None
    pp.set_positive_environment(ratio=ratio)
    fit_cops.Model = exp_model_class(fit_cops)
    cfg, argv, notes = pp.rewrite_arguments(
        ["--point-table", str(points), "--check-jacobian", "--check-parameters",
         "a00,a02,a07,c_cliff_ge_60deg,kappa,grain_A_uPa,grain_B_uPa_per_100s"] + extra + ["--",
         str(run), "--allow-legacy-base", "--degree", "4", "--enhancement-degree", "4", "--lambda", str(lam),
         "--jobs", "1", "--tbounds", "30", "400", "--max-nfev", str(max_nfev), "--no-save-restart", "--out", str(out),
         "--activity-template", str(bundle), "cliff_ge_60deg"])
    if cfg.get("held_exponent"):                      # as positive_fit.main (the GPU patch is emulated by ExpModel)
        pp.load_held(cfg["held_exponent"], pp.elements_faces(run))
    if cfg.get("restart_from"):
        pp.check_restart_marker(cfg["restart_from"], cfg)
    if cfg.get("start_restart"):
        pp.check_restart_marker(cfg["start_restart"], cfg)
    record = {"command": ["test_positive_fit.py e2e"]}
    sub = {"config": cfg, "notes": notes}
    record["positive"] = sub
    pp.install_core(fit_cops, cfg, sub)
    grain_crossing.install(fit_cops, 10.0, record)
    pp.install_outer(fit_cops, cfg, sub)
    original_save = fit_cops.save_active_manifest

    def save():
        if fit_cops.ACTIVE_MANIFEST is not None:
            fit_cops.ACTIVE_MANIFEST["model_extension"] = record
            fit_cops.ACTIVE_MANIFEST["compute_backend"] = {"kind": "CPU test stand-in (ExpModel)"}
        original_save()
    fit_cops.save_active_manifest = save
    captured = []
    original_fit_one = fit_cops.fit_one

    def fit_one(lam_, x_start=None):
        r = original_fit_one(lam_, x_start)
        captured.append(r)
        return r
    fit_cops.fit_one = fit_one
    if cfg.get("restart_from"):                     # gpu_fit's resume path (drive -> install_template_start)
        import gpu_fit
        from cops_cuda import restart as rs_
        gpu_fit.install_template_start(fit_cops, 0.0, {}, rs_.read(cfg["restart_from"]), [])
    sys.argv = [str(REF / "tools" / "fit_cops.py")] + argv[argv.index("--") + 1:]
    fit_cops.main()
    # gpu_fit's restart writer path with the marker
    from cops_cuda import restart as rs
    elements = fit_cops.read_elements(str(Path(run) / "elements.bin"))
    state = rs.from_fit(fit_cops.P, captured[-1], out, fit_cops.ACTIVE_MANIFEST, elements)
    rs.write(Path(out) / "fit_restart.npz", state)
    P = fit_cops.P
    from cops_cuda import figures                     # gpu_fit.write_figures' reader, with the grid filter
    nfig = [len(figures.read_table(Path(out))["utc"]), len(figures.read_table(Path(run))["utc"])]
    np.savez(Path(out) / "test_state.npz", x=captured[-1]["x"], x0=P.x0, D0=P.deviation(P.x0),
             Lt=P.Lt, ascale=P.ascale, et_all=P.et_all, et_end=P.et_end, nfig=np.array(nfig))
    absorb = os.environ.get("POSITIVE_TEST_ABSORB")
    if absorb and P.held is not None:
        # h = B z (in the span): residuals and Jacobian with h at x == without h at x + h_B (exact absorption)
        x = np.array(captured[-1]["x"], float)
        rH, JH, DH = P.residuals(x), P.jacobian(x), P.deviation(x)
        xs = x.copy()
        xs[P.reg_slice] += P.held_hB
        h_saved, held_saved = pp.HELD["h"], P.held
        pp.HELD["h"], P.held = None, None
        try:
            rN, JN, DN = P.residuals(xs), P.jacobian(xs), P.deviation(xs)
        finally:
            pp.HELD["h"], P.held = h_saved, held_saved
        n = min(len(rH), len(rN))
        out_ = {"residual_max_abs_diff": float(np.max(np.abs(rH[:n] - rN[:n]))), "rows": [len(rH), len(rN)],
                "residual_scale": float(np.max(np.abs(rN))), "D_held": DH, "D_shifted": DN,
                "jacobian_max_rel_diff": float(np.max(np.abs(JH[:JN.shape[0]] - JN)) / np.max(np.abs(JN))),
                "extra_rows_held": [float(v) for v in rH[n:]], "hB": [float(v) for v in P.held_hB],
                "perp2": float(P.held_perp2)}
        Path(absorb).write_text(json.dumps(out_, indent=1))


def test_e2e(scratch):
    import positive_problem as pp
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    root = Path(scratch) / "e2e"
    run, bundle, points = make_base(root)
    for lam in (0.0, 0.3):
        out = root / ("fit_lam%g" % lam)
        log = root / ("fit_lam%g.log" % lam)
        cmd = [sys.executable, "-W", "ignore", __file__, "--e2e-child", str(run), str(bundle), str(points), str(out),
               str(lam)]
        rc = subprocess.call(cmd, stdout=log.open("w"), stderr=subprocess.STDOUT)
        text = log.read_text()
        if rc != 0:
            report("e2e lambda %g: fit_cops.main with all layers" % lam, False, "exit %d; log %s\n%s" % (rc, log, text[-3000:]))
            continue
        man = json.loads((out / "fit_manifest.json").read_text())
        pos = man["model_extension"]["positive"]
        jc = pos["jacobian_check"]
        cols = {k: v for k, v in jc["columns"].items() if isinstance(v, dict)}
        worst = max(v["data"] for v in cols.values())
        worst_late = max(v.get("data_late_window", 0.0) for v in cols.values())
        reg = max(v.get("regularization", 0.0) for v in cols.values())
        report("e2e lambda %g: fit ran, Jacobian check (data/late/reg)" % lam,
               worst < 1e-6 and worst_late < 1e-6 and reg < 1e-6,
               "data %.1e late %.1e reg %.1e; crossing %s; rows %s" % (
                   worst, worst_late, reg, ["%.1e" % v.get("crossing", float("nan")) for v in cols.values()], jc["rows"]))
        st = np.load(out / "test_state.npz")
        z = np.load(out / "activity_exp_facets.npz")
        el = fit_cops.read_elements(str(run / "elements.bin"))
        D_direct = math.sqrt(float(np.sum(el["area"] * (z["s"] - z["s0"]) ** 2) / np.sum(el["area"])))
        D_fit = man["optimizer_results"][0]["deviation_from_prior"]
        report("e2e lambda %g: D(x0) = 0, D == rms(s - s0)" % lam,
               float(st["D0"]) == 0.0 and abs(D_fit - D_direct) <= 1e-9 * max(D_direct, 1e-300) + 1e-14,
               "D fitter %.6e direct %.6e" % (D_fit, D_direct))
        # COPS_model.dat late rows == W10x5 of the instantaneous surface model + grain at the time stamp
        names, rows, cols_ = None, [], None
        for line in open(out / "COPS_model.dat"):
            if line.startswith("#iPoint"):
                names = line[1:].split()
            elif not line.startswith("#"):
                rows.append(line.split())
        c = {n: i for i, n in enumerate(names)}
        pNG_tab = np.array([float(r[c["pNG"]]) for r in rows])
        pRG_tab = np.array([float(r[c["pRG"]]) for r in rows])
        header, cn, rws = fit_cops.read_model_file(str(run / "COPS_model.dat"))
        cmap = {n: i for i, n in enumerate(cn)}
        Mraw = exp_model_class(fit_cops)(str(run), header, cmap, rws, el, fit_cops.read_activity_coefficients(
            str(REF), header))
        P_names = man["parameters"]["names"]
        x = st["x"]
        a = x[:25] * float(st["ascale"])
        cc = np.array([x[P_names.index("c_cliff_ge_60deg")]])
        kappa = x[P_names.index("kappa")]
        tcol = np.load(bundle)["values"]
        pNG_i, pRG_i = Mraw.evaluate(a, Mraw.T0, kappa, b=cc, enhancement_basis=tcol, enhancement_scale=1.0)
        gj = json.loads((out / "grain_term.json").read_text())
        A_, B_ = gj["result"]["A_Pa"], gj["result"]["B_Pa_per_s"]
        et = st["et_all"]
        grain = np.where(et >= float(st["et_end"]), np.maximum(A_ - B_ * (et - float(st["et_end"])), 0.0), 0.0)
        table = pp.read_point_table(points)
        kind = np.array([e["kind"] for e in table])
        gi = np.nonzero(kind == "grid")[0]
        tg = et[gi]
        worst_ng = worst_rg = worst_inst = 0.0
        late = [i for i, e in enumerate(table) if e["kind"] == "obs" and e["l2_stop"] != "-"]
        for i in late:
            T = et[i] + (pp.parse_utc_seconds(table[i]["l2_stop"]) - pp.parse_utc_seconds(table[i]["utc"]))
            import late_response_eval as lre
            eng = lre.responses(tg, pNG_i[gi], T)["W10x5"] + grain[i]
            erg = lre.responses(tg, pRG_i[gi], T)["W10x5"]
            worst_ng = max(worst_ng, abs(pNG_tab[i] / eng - 1))
            worst_rg = max(worst_rg, abs(pRG_tab[i] / erg - 1))
        early = [i for i, e in enumerate(table) if e["kind"] == "obs" and e["l2_stop"] == "-"] + list(gi)
        for i in early:
            worst_inst = max(worst_inst, abs(pNG_tab[i] / (pNG_i[i] + grain[i]) - 1), abs(pRG_tab[i] / pRG_i[i] - 1))
        report("e2e lambda %g: COPS_model.dat late = W10x5 + grain, else inst." % lam,
               max(worst_ng, worst_rg, worst_inst) < 2e-6,
               "late NG %.1e RG %.1e; other rows %.1e (%d late, %d other; table %%.6e)" % (
                   worst_ng, worst_rg, worst_inst, len(late), len(early)))
        # patched activity products
        fac = np.loadtxt(out / "activity_facets.dat", comments="#", usecols=(6, 7, 8, 9))
        ok = (np.allclose(fac[:, 2], z["A"], rtol=1e-12, atol=0) and np.allclose(fac[:, 3], z["A"], rtol=1e-12, atol=0)
              and np.allclose(fac[:, 0], el["Y"] @ fit_cops.read_activity_coefficients(str(REF), header)[0], rtol=1e-12))
        res_txt = (out / "fit_result.txt").read_text()
        ratio = float([l.split()[1] for l in res_txt.splitlines() if l.startswith("activity_ratio ")][0])
        ok_ratio = abs(ratio - np.sum(el["area"] * z["A"]) / np.sum(el["area"] * el["AH2O"])) <= 6e-7   # %.6f
        rmeta = json.loads(str(np.load(out / "fit_restart.npz")["meta_json"]))
        ok_marker = rmeta["model_extension"]["parameterization"] == pp.PARAMETERIZATION_ID
        try:
            pp.check_restart_marker(out / "fit_restart.npz")
            ok_marker &= True
        except SystemExit:
            ok_marker = False
        report("e2e lambda %g: activity_facets.dat, fit_result ratio, restart marker" % lam,
               ok and ok_ratio and ok_marker and (out / "ACTIVITY_PARAMETERIZATION.txt").is_file(),
               "facets %s ratio %s marker %s" % (ok, ok_ratio, ok_marker))
        # grid rows: flagged in COPS_model.dat, listed in late_window_model.dat, absent from plots/residual tables
        n_obs = len([e for e in table if e["kind"] == "obs"])
        gutc = {e["utc"] for e in table if e["kind"] == "grid"}
        note = [l for l in open(out / "COPS_model.dat") if l.startswith("# PositiveFitGridRows")]
        lwm = [l.split() for l in open(out / "late_window_model.dat") if not l.startswith("#")]
        kinds = [v[3] for v in lwm]
        res_rows = [l.split() for l in open(out / "pressure_residuals.dat") if not l.startswith("#")]
        surf_ok = all(abs(float(v[8]) / (float(v[5]) + float(v[7])) - 1) < 1e-6 and
                      abs(float(v[9]) / float(v[6]) - 1) < 1e-6 for v in lwm)
        ok = (len(note) == 1 and kinds.count("grid") == 105 and kinds.count("late_obs") == 20 and surf_ok and
              len(res_rows) == 2 * n_obs and not any(v[1] in gutc for v in res_rows) and
              list(st["nfig"]) == [n_obs, n_obs] and
              pos["late_window_products"]["table_vs_surface_plus_grain_max_rel"] < 1e-6)
        report("e2e lambda %g: grid rows flagged/listed; plots+residuals obs only" % lam, ok,
               "note %d; late_window_model grid %d late_obs %d (surface+grain ok %s); residual rows %d (2 x %d obs); "
               "figure tables %s; table vs surface+grain %.1e" % (
                   len(note), kinds.count("grid"), kinds.count("late_obs"), surf_ok, len(res_rows), n_obs,
                   list(st["nfig"]), pos["late_window_products"]["table_vs_surface_plus_grain_max_rel"]))
        sel = pos.get("selection")
        lw = pos["late_window"]
        sp_ = pos["start_projection"]
        report("e2e lambda %g: record (late window, start, grain start)" % lam,
               lw["late_observations"] == 20 and lw["grid_rows"] == 105 and lw["fit_rows_equal_all_rows"] and
               lw["grain_t_end_equals_late_start"] and pos["grain_start"]["grain_A_uPa"] == pp.DEFAULT_GRAIN_START_UPA
               and sel is None,
               "late %d grid %d shared %s; rms ln %.3f; cliff c %s; t_end ok %s" % (
                   lw["late_observations"], lw["grid_rows"], lw["fit_rows_equal_all_rows"],
                   sp_["rms_ln_residual_over_A_orig_gt_0"], sp_["c"], lw["grain_t_end_equals_late_start"]))
        print("   log: %s" % log)
        for line in text.splitlines():
            if line.startswith("positive_fit:") or line.startswith("consistency") or "cost" in line[:40]:
                print("   | " + line[:220])


def run_child(root, name, run, bundle, points, lam, extra=(), env=None, max_nfev=4):
    out = Path(root) / name
    log = Path(root) / (name + ".log")
    cmd = [sys.executable, "-W", "ignore", __file__, "--e2e-child", str(run), str(bundle), str(points), str(out),
           str(lam), str(max_nfev)] + list(extra)
    environ = dict(os.environ, **(env or {}))
    rc = subprocess.call(cmd, stdout=log.open("w"), stderr=subprocess.STDOUT, env=environ)
    return rc, out, log.read_text()


def start_lines(text):
    return [l for l in text.splitlines() if l.startswith("positive_fit: start check") or
            l.startswith("positive_fit: grain start") or "WARNING" in l]


def test_e2e_start(scratch):
    """--start-restart (cycle 1) and the gpu_fit resume path on the synthetic base (CPU stand-in): reproduction of the
    source's final cost (same problem), a ratio change reported (not a failure), an explicit grain start (wins over the
    restart's grain; reproduction not expected), a resume (--restart-from) that must reproduce."""
    root = Path(scratch) / "e2e"
    run, bundle, points = root / "base", root / "templates.npz", root / "late_grid_points.txt"
    src = root / "fit_lam0.3"
    restart = src / "fit_restart.npz"
    if not restart.is_file():
        report("e2e start restart: source restart", False, "missing %s (test_e2e first)" % restart)
        return
    sman = json.loads((src / "fit_manifest.json").read_text())
    src_cost = sman["optimizer_results"][0]["cost"]
    src_x = np.load(restart)["x"]
    # 1 same problem: the start must reproduce the source's final cost
    rc, out, text = run_child(root, "start_same", run, bundle, points, 0.3, ["--start-restart", str(restart)])
    ok = rc == 0
    if ok:
        man = json.loads((out / "fit_manifest.json").read_text())
        pos = man["model_extension"]["positive"]
        sc = pos["start_check"]
        g = pos["grain_start"]
        ok = (sc["reproduction_expected"] and "failed" not in sc and abs(sc["cost_rel_diff"]) <= 1e-12 and
              sc["x_start_equals_source_x"] and abs(sc["start"]["cost"] / src_cost - 1) <= 1e-12 and
              g["grain_A_uPa"] == float(src_x[-2]) and g["grain_B_uPa_per_100s"] == float(src_x[-1]) and
              pos["start_restart"]["x_equals_source_x"] and
              man["model_extension"]["positive"]["fit_chain"]["start_restart"]["path"] == str(restart.resolve()))
        rmeta = json.loads(str(np.load(out / "fit_restart.npz")["meta_json"]))["model_extension"]
        ok &= rmeta["fit_chain"]["start_restart"]["sha256"] == sc["source"]["sha256"] and rmeta["night_ratio"] == 0.02
        ok &= "Optimizer start: the restart" in (out / "ACTIVITY_PARAMETERIZATION.txt").read_text()
        detail = "rel %.1e; start cost %.12e source %.12e; grain %s" % (sc["cost_rel_diff"], sc["start"]["cost"],
                                                                         src_cost, [g["grain_A_uPa"], g["grain_B_uPa_per_100s"]])
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e start restart: same problem -> reproduces the source cost", ok, detail)
    for line in start_lines(text):
        print("   | " + line[:240])
    # 2 ratio 0.002: reported change (the CPU stand-in ignores the ratio, so the cost is unchanged here)
    rc, out, text = run_child(root, "start_r0002", run, bundle, points, 0.3, ["--start-restart", str(restart),
                                                                             "--ratio", "0.002"])
    ok = rc == 0
    if ok:
        sc = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]["start_check"]
        rmeta = json.loads(str(np.load(out / "fit_restart.npz")["meta_json"]))["model_extension"]
        ok = (not sc["reproduction_expected"] and "failed" not in sc and sc["changes"][0].startswith("ratio 0.02/0.02")
              and rmeta["night_ratio"] == 0.002 and rmeta["floor_ratio"] == 0.002 and
              "0.002 -- user 2026-11-07 (second" in (out / "ACTIVITY_PARAMETERIZATION.txt").read_text())
        detail = "changes %s; marker ratios %s/%s" % (sc["changes"], rmeta["night_ratio"], rmeta["floor_ratio"])
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e start restart: ratio 0.002 reported, recorded in marker/text", ok, detail)
    for line in start_lines(text):
        print("   | " + line[:240])
    # 3 explicit grain start (fine tuning): wins over the restart's grain; reproduction not expected
    rc, out, text = run_child(root, "start_grain", run, bundle, points, 0.3, ["--start-restart", str(restart),
                                                                             "--grain-start-uPa", "5.5",
                                                                             "--grain-start-B", "0.25"])
    ok = rc == 0
    if ok:
        pos = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]
        sc, g = pos["start_check"], pos["grain_start"]
        ok = (g["grain_A_uPa"] == 5.5 and g["grain_B_uPa_per_100s"] == 0.25 and not sc["reproduction_expected"] and
              any("grain_A_uPa" in c for c in sc["changes"]))
        detail = "grain %s; changes %s" % ([g["grain_A_uPa"], g["grain_B_uPa_per_100s"]], sc["changes"])
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e start restart: explicit grain start wins, reported", ok, detail)
    # 4 resume (gpu_fit --restart-from path) of the same problem: must reproduce; another ratio is refused
    rs1 = root / "start_same" / "fit_restart.npz"
    rc, out, text = run_child(root, "resume_same", run, bundle, points, 0.3, ["--restart-from", str(rs1)])
    ok = rc == 0
    if ok:
        pos = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]
        sc = pos["start_check"]
        rmeta = json.loads(str(np.load(out / "fit_restart.npz")["meta_json"]))["model_extension"]
        ok = (sc["reproduction_expected"] and "failed" not in sc and abs(sc["cost_rel_diff"]) <= 1e-12 and
              sc["source"]["role"] == "restart_from" and
              rmeta["fit_chain"]["start_restart"]["path"] == str(restart.resolve()) and
              rmeta["fit_chain"]["resumed_from"][-1]["path"] == str(rs1.resolve()))
        detail = "rel %.1e; chain start %s, resumed from %s" % (sc["cost_rel_diff"], rmeta["fit_chain"]["start_restart"]["path"][-30:],
                                                                rmeta["fit_chain"]["resumed_from"][-1]["path"][-35:])
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e resume (--restart-from): reproduces, chain recorded", ok, detail)
    rc, out, text = run_child(root, "resume_other_ratio", run, bundle, points, 0.3,
                              ["--restart-from", str(rs1), "--ratio", "0.002"])
    report("e2e resume with another ratio is refused", rc != 0 and "a resume must continue the same problem" in text,
           "exit %d: %s" % (rc, [l for l in text.splitlines() if "resume must" in l][:1]))


# ------------------------------------------------------------------------------------------------ held exponent
def test_held_units():
    """held_projection (D decomposition), install_held_activity (act_lin + h on the host path, device refused)."""
    import positive_problem as pp
    rng = np.random.default_rng(7)
    n = 30000
    Y = rng.normal(size=(n, 25))
    Y[:, 0] = 1.0
    T = (rng.random((n, 1)) < 0.2).astype(float)
    area = rng.uniform(0.5, 2.0, n)
    A = np.exp(Y @ rng.normal(0, 0.3, 25) + 0.5 * T[:, 0])
    pr = pp.start_projection(Y, area, A, T)
    ascale = pr["ascale"]
    Lt, G, _ = pp.log_tikhonov(pr["G_all"], float(area.sum()), ascale, 25)
    x0 = np.concatenate([pr["a"] / ascale, pr["c"]])
    S = np.concatenate([np.full(25, ascale), [1.0]])
    worst = 0.0
    for trial in range(4):
        h = rng.normal(0, 0.4, n) if trial < 3 else np.zeros(n)
        hp = pp.held_projection(Y, T, area, h, pr["G_all"], ascale, chunk=7000)
        x = x0 + rng.normal(0, 0.05, 26)
        s = np.column_stack((Y, T)) @ (x * S) + h
        s0 = np.column_stack((Y, T)) @ (x0 * S)
        D_direct = math.sqrt(float(np.sum(area * (s - s0) ** 2) / area.sum()))
        v = Lt @ (x - x0 + hp["hB"])
        D_fit = math.sqrt(float(v @ v) + hp["perp2"])
        worst = max(worst, abs(D_fit / D_direct - 1))
        if trial == 3:
            zero_ok = np.all(hp["hB"] == 0.0) and hp["perp2"] == 0.0
    report("held: D direct == sqrt(|Lt(x-x0+h_B)|^2 + |h_perp|^2) (random h)", worst < 1e-12 and zero_ok,
           "max rel %.1e; h = 0 -> h_B = 0, perp 0: %s; orthogonality %.1e" % (worst, zero_ok, hp["orthogonality_rel"]))
    z = rng.normal(0, 0.2, 26)
    h = np.column_stack((Y, T)) @ z
    hp = pp.held_projection(Y, T, area, h, pr["G_all"], ascale)
    err = float(np.max(np.abs(hp["hB"] * S - z)))
    report("held: h in span(B) -> h_B = z, |h_perp| ~ 0", err < 1e-12 and hp["perp2"] < 1e-24 * hp["h2"],
           "max |h_B S - z| %.1e; perp2/h2 %.1e" % (err, hp["perp2"] / hp["h2"]))
    # the GPU host-path patch (on the classes of this process; restored afterwards)
    from types import SimpleNamespace
    from cops_cuda import evaluator, multi_gpu
    saved = {(c, k): getattr(c, k) for c in (evaluator.CudaBackend, multi_gpu.MultiGpuBackend)
             for k in ("activity_linear", "_upload_activity", "_evaluate") if k in c.__dict__}
    saved_held = dict(pp.HELD)
    try:
        hh = rng.normal(size=50)
        pp.HELD.update(h=hh, info={"path": "-"}, calls=0)
        rec = {}
        pp.install_held_activity(rec)
        fake = SimpleNamespace(Y=rng.normal(size=(50, 25)), enhancement_scale=1.0,
                               enhancement_basis=rng.normal(size=(50, 1)), activity_mode="device", activity="device")
        st = SimpleNamespace(a=rng.normal(size=25), b=np.array([0.3]))
        ref = fake.Y @ st.a + fake.enhancement_basis @ st.b + hh
        got_c = evaluator.CudaBackend.activity_linear(fake, st)
        got_m = multi_gpu.MultiGpuBackend.activity_linear(fake, st)
        ok = np.array_equal(got_c, ref) and np.array_equal(got_m, ref) and pp.HELD["calls"] == 2
        refused = 0
        for f in (lambda: evaluator.CudaBackend._upload_activity(fake, st),
                  lambda: multi_gpu.MultiGpuBackend._evaluate(fake, st, None)):
            try:
                f()
            except Exception as exc:
                refused += "host activity" in str(exc)
        report("held: GPU host-path patch act_lin + h; device activity refused", ok and refused == 2,
               "act_lin == Y a + B b + h: %s; calls %d; device refused %d of 2" % (ok, pp.HELD["calls"], refused))
    finally:
        for (c, k), f in saved.items():
            setattr(c, k, f)
        pp.HELD.clear()
        pp.HELD.update(saved_held)


def test_e2e_held(scratch):
    """The held term end to end on the synthetic base (CPU stand-in adds h to act_lin as the GPU patch does)."""
    import positive_problem as pp
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    root = Path(scratch) / "e2e"
    run, bundle, points = root / "base", root / "templates.npz", root / "late_grid_points.txt"
    el = fit_cops.read_elements(str(run / "elements.bin"))
    nf = len(el["area"])
    tcol = np.load(bundle)["values"][:, 0]
    rng = np.random.default_rng(11)
    files = {"zeros": root / "held_zeros.npz", "random": root / "held_random.npz", "span": root / "held_span.npz",
             "ft": root / "finetune_state.npz", "random2": root / "held_random2.npz"}
    np.savez(files["zeros"], held=np.zeros(nf))
    h_rand = rng.normal(0.0, 0.3, nf)
    np.savez(files["random"], held=h_rand)
    np.savez(files["random2"], held=h_rand * 0.5)
    z = np.concatenate([rng.normal(0, 0.05, 25), [0.2]])
    np.savez(files["span"], held=el["Y"] @ z[:25] + tcol * z[25])
    np.savez(files["ft"], held_next=h_rand, grain_A_uPa=np.float64(4.0), grain_B_uPa_per_100s=np.float64(0.5))
    base = root / "fit_lam0.3"
    # 1 held = 0: bitwise the run without a held term (costs, x, D)
    rc, out, text = run_child(root, "held_zeros", run, bundle, points, 0.3, ["--held-exponent", str(files["zeros"])])
    if rc == 0:
        m0 = json.loads((base / "fit_manifest.json").read_text())["optimizer_results"][0]
        m1 = json.loads((out / "fit_manifest.json").read_text())["optimizer_results"][0]
        x0_, x1_ = np.load(base / "test_state.npz")["x"], np.load(out / "test_state.npz")["x"]
        keys = ("cost", "data_cost", "regularization_cost", "deviation_from_prior")
        ok = all(m0[k] == m1[k] for k in keys) and np.array_equal(x0_, x1_)
        detail = "; ".join("%s %r vs %r" % (k, m0[k], m1[k]) for k in keys) + "; x equal %s" % np.array_equal(x0_, x1_)
    else:
        ok, detail = False, "exit %d\n%s" % (rc, text[-2500:])
    report("e2e held = 0: costs, x, D bitwise as without held", ok, detail)
    # 2 random h: D direct == fitter D (1e-9), Jacobian check, products, marker
    rc, out, text = run_child(root, "held_random", run, bundle, points, 0.3, ["--held-exponent", str(files["random"])])
    if rc == 0:
        man = json.loads((out / "fit_manifest.json").read_text())
        pos = man["model_extension"]["positive"]
        zf = np.load(out / "activity_exp_facets.npz")
        meta = json.loads(str(zf["meta_json"]))
        D_dir = math.sqrt(float(np.sum(el["area"] * (zf["s"] - zf["s0"]) ** 2) / np.sum(el["area"])))
        D_fit = man["optimizer_results"][0]["deviation_from_prior"]
        st = np.load(out / "test_state.npz")
        names = man["parameters"]["names"]
        x = st["x"]
        s_fit = el["Y"] @ (x[:25] * float(st["ascale"])) + tcol * x[names.index("c_cliff_ge_60deg")]
        jc = pos["jacobian_check"]
        cols = {k: v for k, v in jc["columns"].items() if isinstance(v, dict)}
        jac_ok = max(max(v["data"], v.get("data_late_window", 0.0), v.get("regularization", 0.0))
                     for v in cols.values()) < 1e-6
        rz = np.load(out / "fit_restart.npz")
        rmeta = json.loads(str(rz["meta_json"]))["model_extension"]
        txt = (out / "ACTIVITY_PARAMETERIZATION.txt").read_text()
        cost_parts = man["optimizer_results"][0]
        ok = (abs(D_dir / D_fit - 1) <= 1e-9 and abs(meta["D_direct_vs_fitter_rel"]) <= 1e-9 and jac_ok and
              np.array_equal(zf["held"], h_rand) and np.allclose(zf["s"] - zf["held"], s_fit, rtol=0, atol=1e-12) and
              rmeta["held_exponent"]["sha256"] == pp.sha256(files["random"]) and rmeta["activity_linear_includes_held"]
              and np.allclose(rz["activity_linear"], zf["s"], rtol=0, atol=1e-12) and
              pos["held_exponent"]["D_perp_nepers"] > 0 and "Held exponent term: h from" in txt and
              pp.sha256(files["random"]) in txt and pos["held_activity_calls"] > 0 and
              jc["rows"]["regularization"] == 27 and
              abs(cost_parts["cost"] - (cost_parts["data_cost"] + cost_parts["regularization_cost"])) <=
              1e-12 * cost_parts["cost"])
        detail = ("D direct %.12f fitter %.12f (rel %.1e; meta %.1e); D_perp %.4f; Jacobian ok %s (reg rows %d); "
                  "cost = data + lam D^2/2: %.1e; held calls %d" % (
                      D_dir, D_fit, D_dir / D_fit - 1, meta["D_direct_vs_fitter_rel"],
                      pos["held_exponent"]["D_perp_nepers"], jac_ok, jc["rows"]["regularization"],
                      cost_parts["cost"] - cost_parts["data_cost"] - cost_parts["regularization_cost"],
                      pos["held_activity_calls"]))
    else:
        ok, detail = False, "exit %d\n%s" % (rc, text[-2500:])
    report("e2e held random: D direct == fitter, Jacobian, products, marker", ok, detail)
    for line in text.splitlines():
        if line.startswith("positive_fit: held") or line.startswith("positive_fit: wrote") and "D direct" in line:
            print("   | " + line[:240])
    # 3 h in span(B): absorbed exactly by the shift x + h_B (residuals, Jacobian, D)
    absorb = root / "held_span_absorb.json"
    rc, out, text = run_child(root, "held_span", run, bundle, points, 0.3, ["--held-exponent", str(files["span"]),
                                                                           "--check-held-span"],
                              env={"POSITIVE_TEST_ABSORB": str(absorb)}, max_nfev=2)
    if rc == 0 and absorb.is_file():
        ab = json.loads(absorb.read_text())
        hs = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]["held_span_check"]
        print("   | held_span_check (start x): residuals %.1e, Jacobian %.1e, D %.9f vs %.9f, calls with/without %d/%d"
              % (hs["residual_max_abs_diff"], hs["jacobian_max_rel_diff"], hs["D_held"], hs["D_shifted"],
                 hs["held_calls_with"], hs["held_calls_without"]))
        span_ok = (hs["residual_max_abs_diff"] <= 1e-12 * max(1.0, hs["residual_scale"]) and
                   hs["jacobian_max_rel_diff"] <= 1e-12 and hs["held_calls_without"] == 0 and
                   abs(hs["D_held"] / hs["D_shifted"] - 1) <= 1e-12)
        S = np.concatenate([np.full(25, float(np.load(out / "test_state.npz")["ascale"])), [1.0]])
        ok = (ab["residual_max_abs_diff"] <= 1e-12 * max(1.0, ab["residual_scale"]) and
              ab["jacobian_max_rel_diff"] <= 1e-12 and abs(ab["D_held"] / ab["D_shifted"] - 1) <= 1e-12 and
              np.allclose(np.array(ab["hB"]) * S, z, rtol=0, atol=1e-10) and ab["perp2"] < 1e-20 and span_ok)
        detail = ("residuals %.1e (scale %.2e), Jacobian rel %.1e, D %.9f vs %.9f, extra rows %s, h_B S - z %.1e"
                  % (ab["residual_max_abs_diff"], ab["residual_scale"], ab["jacobian_max_rel_diff"], ab["D_held"],
                     ab["D_shifted"], ab["extra_rows_held"], float(np.max(np.abs(np.array(ab["hB"]) * S - z)))))
    else:
        ok, detail = False, "exit %d\n%s" % (rc, text[-2500:])
    report("e2e held in span(B): absorbed exactly by x + h_B", ok, detail)
    # 4 cycle-2 style: start restart + finetune_state.npz (held_next, grain) ; resume with the same held / another
    src = root / "fit_lam0.3" / "fit_restart.npz"
    rc, out, text = run_child(root, "held_cycle2", run, bundle, points, 0.3,
                              ["--start-restart", str(src), "--held-exponent", str(files["ft"]),
                               "--grain-start-from", str(files["ft"])])
    if rc == 0:
        pos = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]
        g, sc = pos["grain_start"], pos["start_check"]
        ok = (g["grain_A_uPa"] == 4.0 and g["grain_B_uPa_per_100s"] == 0.5 and str(files["ft"]) in g["source_A"] and
              pos["held_exponent"]["key"] == "held_next" and not sc["reproduction_expected"] and
              any(c.startswith("held term") for c in sc["changes"]))
        detail = "grain %s from %s; changes %s" % ([g["grain_A_uPa"], g["grain_B_uPa_per_100s"]],
                                                  g["source_A"][-40:], sc["changes"])
    else:
        ok, detail = False, "exit %d\n%s" % (rc, text[-2500:])
    report("e2e cycle-2 start: start restart + held_next + grain from file", ok, detail)
    rs2 = out / "fit_restart.npz"
    rc, out, text = run_child(root, "held_resume", run, bundle, points, 0.3,
                              ["--restart-from", str(rs2), "--held-exponent", str(files["ft"])])
    if rc == 0:
        sc = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]["start_check"]
        ok = sc["reproduction_expected"] and "failed" not in sc and abs(sc["cost_rel_diff"]) <= 1e-12
        detail = "rel %.1e" % sc["cost_rel_diff"]
    else:
        ok, detail = False, "exit %d\n%s" % (rc, text[-2500:])
    report("e2e resume with the same held term reproduces", ok, detail)
    rc, out, text = run_child(root, "held_resume_other", run, bundle, points, 0.3,
                              ["--restart-from", str(rs2), "--held-exponent", str(files["random2"])])
    report("e2e resume with another held term is refused", rc != 0 and "a resume must continue the same problem" in text,
           "exit %d" % rc)
    rc, out, text = run_child(root, "held_resume_none", run, bundle, points, 0.3, ["--restart-from", str(rs2)])
    report("e2e resume without the held term is refused", rc != 0 and "a resume must continue the same problem" in text,
           "exit %d" % rc)


# ------------------------------------------------------------------------------------------------ 12M
# ------------------------------------------------------------------------------------------------ CO2 flux floor
OLD_C1 = Path("/data001/highresolu_vt/cops_output/model_ext_frozen/c1_20261108T000327Z")


def _load_module(path, name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_co2_environment(scratch):
    """COPS_POS_FLOOR_CO2 handling (2026-11-08): positive_problem.set_positive_environment and the restart marker,
    positive_activity._read_environment and its canonical form."""
    import positive_problem as pp
    import positive_activity as pa
    res = {}
    for r in ("0.02", "0.002", "1e-5"):
        env = {}
        out = pp.set_positive_environment(env, ratio=r)
        res[r] = (env.get("COPS_POS_FLOOR_CO2"), out.get("COPS_POS_FLOOR_CO2"), pp.current_co2_floor(env))
    ok = (res["0.02"] == (None, None, False) and res["0.002"] == (None, None, False) and res["1e-5"] == ("1", "1", True))
    env = {"COPS_POS_NIGHT_RATIO": "1e-5", "COPS_POS_FLOOR_RATIO": "1e-5"}       # preset ratio, no --ratio (workers)
    pp.set_positive_environment(env)
    ok &= env.get("COPS_POS_FLOOR_CO2") == "1"
    refused = []
    for preset, r in ((("COPS_POS_FLOOR_CO2", "1"), "0.002"), (("COPS_POS_FLOOR_CO2", "0"), "0.002"),
                      (("COPS_POS_FLOOR_CO2", "1"), "0.02"), (("COPS_POS_FLOOR_CO2", "0"), "1e-5"),
                      (("COPS_POS_FLOOR_CO2", "yes"), "1e-5")):
        try:
            pp.set_positive_environment(dict([preset]), ratio=r)
            refused.append(False)
        except SystemExit:
            refused.append(True)
    ok &= all(refused)
    report("CO2 floor: set_positive_environment (1e-5 only; refusals)", ok,
           "%s; preset ratio 1e-5 -> %s; refused %s" % (res, env.get("COPS_POS_FLOOR_CO2"), refused))
    # positive_activity: canonical form of the earlier settings == the 2026-11-07 module's
    old = _load_module(OLD_C1 / "positive_activity.py", "positive_activity_c1")
    base = {"COPS_POS_ENABLE": "1", "COPS_POS_ACTIVITY": "exp"}
    same = []
    for e in ({"COPS_POS_ENABLE": "1"}, dict(base, COPS_POS_NIGHT_RATIO="0.02", COPS_POS_FLOOR_RATIO="0.02"),
              dict(base, COPS_POS_NIGHT_RATIO="0.002", COPS_POS_FLOOR_RATIO="0.002"),
              dict(base, COPS_POS_NIGHT_RATIO="1e-5", COPS_POS_FLOOR_RATIO="1e-5"),
              {"COPS_POS_ENABLE": "1", "COPS_POS_NIGHT_RATIO": "0.02", "COPS_POS_FLOOR_RATIO": "0.02"}):
        c_new, c_old = pa._read_environment(e), old._read_environment(e)
        same.append(pa._canonical(c_new) == old._canonical(c_old) and c_new["co2_floor"] is False)
    c = pa._read_environment(dict(base, COPS_POS_NIGHT_RATIO="1e-5", COPS_POS_FLOOR_RATIO="1e-5", COPS_POS_FLOOR_CO2="1"))
    can = json.loads(pa._canonical(c))
    ok2 = all(same) and can.get("co2_floor") is True and c["modified"] and c["co2_floor_ratio"] == 1e-5
    bad = []
    for e in (dict(base, COPS_POS_FLOOR_CO2="1"), dict(base, COPS_POS_FLOOR_RATIO="0", COPS_POS_FLOOR_CO2="1"),
              dict(base, COPS_POS_FLOOR_RATIO="1e-5", COPS_POS_FLOOR_CO2="2"), {"COPS_POS_FLOOR_CO2": "1"}):
        try:
            pa._read_environment(e)
            bad.append(False)
        except ValueError:
            bad.append(True)
    ok2 &= all(bad)
    report("CO2 floor: positive_activity config (canonical of earlier settings unchanged, refusals)", ok2,
           "earlier settings equal the c1 module's canonical form: %s; CO2 canonical %s; refused %s" % (same, can, bad))
    # restart marker: co2_floor field, resume refusal
    root = Path(scratch) / "co2_marker"
    root.mkdir(parents=True, exist_ok=True)
    marker = {"model_extension": {"parameterization": pp.PARAMETERIZATION_ID, "night_ratio": 1e-5, "floor_ratio": 1e-5}}
    paths = {}
    for tag, co2 in (("old", None), ("off", False), ("on", True)):
        m = json.loads(json.dumps(marker))
        if co2 is not None:
            m["model_extension"]["co2_floor"] = co2
        paths[tag] = root / ("marker_%s.npz" % tag)
        np.savez(paths[tag], meta_json=np.array(json.dumps(m)))
    saved = {k: os.environ.get(k) for k in ("COPS_POS_NIGHT_RATIO", "COPS_POS_FLOOR_RATIO", "COPS_POS_FLOOR_CO2")}
    try:
        os.environ.update(COPS_POS_NIGHT_RATIO="1e-5", COPS_POS_FLOOR_RATIO="1e-5", COPS_POS_FLOOR_CO2="1")
        outc = {}
        for tag, path in paths.items():
            try:
                info = pp.check_restart_marker(path, {"restart_from": str(path.resolve())})
                outc[tag] = "accepted co2_floor=%s" % info["co2_floor"]
            except SystemExit as exc:
                outc[tag] = "refused" if "CO2 floor" in str(exc) else "refused (other: %s)" % exc
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    ok3 = outc == {"old": "refused", "off": "refused", "on": "accepted co2_floor=True"}
    report("CO2 floor: restart marker co2_floor, resume refused without it", ok3, str(outc))


class _FakeDeviceArray:
    """A host array with the slice .set() of a CuPy array (the kappa slot of mesh.ACO2)."""

    def __init__(self, values):
        self.a = np.array(values, dtype=np.float64)
        self.size = self.a.size

    def __getitem__(self, sl):
        parent = self

        class View:
            def set(self, v):
                parent.a[sl] = np.asarray(v, dtype=np.float64)
        return View()


def _co2_case(rng, npt=7, nf=60, nrec=900, ratio=0.02):
    """Synthetic records: per point qmin/qmax (H2O, CO2), VdotB; per record face, Kn, KPhi, KV, cos (<= 0: night)."""
    import types
    case = types.SimpleNamespace(npt=npt, nf=nf, ratio=ratio)
    case.pt = rng.integers(0, npt, nrec)
    case.pt.sort()
    case.face = rng.integers(0, nf, nrec)
    case.kn, case.kphi, case.kv = (rng.uniform(0.1, 1.0, nrec) for _ in range(3))
    case.cos = rng.uniform(-0.5, 1.0, nrec)
    case.qmaxH, case.qmaxC = rng.uniform(1.0, 2.0, npt), rng.uniform(5.0, 9.0, npt)
    case.qminH, case.qminC = ratio * case.qmaxH, ratio * case.qmaxC
    case.vb = rng.uniform(0.5, 1.5, npt)
    case.aco2 = np.where(rng.uniform(size=nf) < 0.1, 0.0, rng.uniform(0.2, 1.0, nf))     # some ACO2 = 0
    case.ah2o = rng.uniform(0.2, 1.0, nf)
    case.Y = 0.3 * rng.normal(size=(nf, 25))                    # 25 exponent columns (ParameterState needs 25)
    case.a = 0.3 * rng.normal(size=25)
    case.KT, case.COS75 = 0.3, float(np.cos(np.radians(75.0)))
    case.cNG = {"H2O": 1.3, "CO2": 0.7}
    case.bm = {"H2O": 0.11, "CO2": 0.21}
    case.cRG = {"H2O": rng.uniform(0.5, 1.0, npt), "CO2": rng.uniform(0.5, 1.0, npt)}
    return case


def _co2_records(case, a, T):
    """Per-record factors of the model (fit_cops Model.evaluate + the H2O floor): q's, u, dl, dlx, dn, af."""
    night = case.cos <= 0.0
    cz = np.minimum(case.cos, 1.0)
    p = case.pt
    qh = np.where(night, case.qminH[p], case.qminH[p] + (case.qmaxH[p] - case.qminH[p]) * cz)
    qc = np.where(night, case.qminC[p], case.qminC[p] + (case.qmaxC[p] - case.qminC[p]) * cz) * case.aco2[case.face]
    c = np.where(night, -1.0, np.clip(case.cos, -1.0, 1.0))
    lit = c >= case.COS75
    x = np.where(lit, 1.0 / np.where(lit, c, 1.0) - 1.0, 0.0)
    Tmax, T75, Tmin = T
    tl = Tmax + case.KT * (Tmax - T75) * x
    use = lit & (tl > Tmin)
    te = np.where(use, tl, Tmin)
    u = 1.0 / np.sqrt(te)
    du = -0.5 * u / te
    dl, dn = np.where(use, du, 0.0), np.where(use, 0.0, du)
    af = np.exp(case.Y @ a)[case.face]
    fminH = case.qmaxH[p] * (case.ratio * case.ah2o[case.face])
    return {"qh": qh, "qc": qc, "u": u, "dl": dl, "dlx": dl * x, "dn": dn, "af": af, "fminH": fminH,
            "fminC": case.qmaxC[p] * (case.ratio * case.aco2[case.face]), "VB": case.vb[p]}


def _co2_reference(case, a, T, kappa):
    """Independent model: F_H = max(q_H A, r qmax_H A_orig), F_C = max(kappa qc, r qmax_C ACO2)."""
    r = _co2_records(case, a, T)
    FH = np.maximum(r["qh"] * r["af"], r["fminH"])
    FC = np.maximum(kappa * r["qc"], r["fminC"])
    bc = lambda w: np.bincount(case.pt, weights=w, minlength=case.npt)     # noqa: E731
    cN, cR, bm = case.cNG, case.cRG, case.bm
    pNG = cN["H2O"] * bm["H2O"] * bc(FH * case.kn * r["u"]) + cN["CO2"] * bm["CO2"] * bc(FC * case.kn * r["u"])
    pRG = (cR["H2O"] * (bc(FH * case.kphi) + bm["H2O"] * bc(FH * r["VB"] * case.kv * r["u"]))
           + cR["CO2"] * (bc(FC * case.kphi) + bm["CO2"] * bc(FC * r["VB"] * case.kv * r["u"])))
    act = ~(kappa * r["qc"] < r["fminC"])
    dk = np.where(act, r["qc"], 0.0)
    dkNG = cN["CO2"] * bm["CO2"] * bc(dk * case.kn * r["u"])
    dkRG = cR["CO2"] * (bc(dk * case.kphi) + bm["CO2"] * bc(dk * r["VB"] * case.kv * r["u"]))
    return pNG, pRG, dkNG, dkRG, int((~act).sum())


def _co2_emulated_backend(case, cf):
    """A stand-in CudaBackend whose _run_fused emulates the documented sums of pressure.cu (cf: the CO2-floor kernels
    eval_fused_m_p5f / p14f, header item 5; otherwise eval_fused_m_p5 / p14 with the bind-time PhiCO2)."""
    import contextlib
    import types
    from cops_cuda import evaluator
    be = types.SimpleNamespace()
    be.mixed, be.device_id, be.timing = True, 0, {}
    be.constants = types.SimpleNamespace(cNG=case.cNG, cRG=case.cRG, bm=case.bm, KT=case.KT, COS75=case.COS75)
    tail = [np.nan, case.ratio] if cf else []
    be.dataset = types.SimpleNamespace(nfaces_device=case.nf, npoints=case.npt, streamed=False, arrays={
        "mesh.ACO2": _FakeDeviceArray(np.concatenate([case.aco2, case.ratio * case.ah2o, tail]))})
    be.cupy = types.SimpleNamespace(cuda=types.SimpleNamespace(Device=lambda d: contextlib.nullcontext()))
    be.arrays = {}
    be._positive_activity = {}
    be.state = {}
    be._upload_activity = lambda state: be.state.update(a=np.asarray(state.a, float))
    be._upload_columns = lambda request: be.state.update(cols=list(request.basis_index))
    be._fused_plan = lambda codes: evaluator.CudaBackend._fused_plan(be, codes)
    be._params = lambda state: evaluator.CudaBackend._params(be, state)

    def run_fused(plan, codes, params):
        T = (params["Tmax"], params["T75"], params["Tmin"])
        r = _co2_records(case, be.state["a"], T)
        d = be.dataset.arrays["mesh.ACO2"].a
        kappa = d[2 * case.nf] if cf else np.nan
        u, af = r["u"], r["af"]
        atH = r["qh"] * af < r["fminH"]
        afr = np.where(atH, 1.0, af)
        qH = np.where(atH, r["fminH"], r["qh"])
        wn_h, wp_h, wv_h = qH * case.kn, qH * case.kphi, qH * r["VB"] * case.kv
        wn_c, wp_c, wv_c = r["qc"] * case.kn, r["qc"] * case.kphi, r["qc"] * r["VB"] * case.kv
        atC = (kappa * r["qc"] < r["fminC"]) if cf else np.zeros(len(u), bool)
        on = ~atC
        fn, fp, fv = r["fminC"] * case.kn, r["fminC"] * case.kphi, r["fminC"] * r["VB"] * case.kv
        w0 = params["cNbmH"] * wn_h * afr + np.where(on, params["kcNbmC"] * wn_c, 0.0)
        w1 = params["bmH"] * wv_h * afr
        w2 = params["bmC"] * wv_c
        terms = {0: wn_h * afr * u, 1: wp_h * afr, 2: wv_h * afr * u, 3: np.where(on, wn_c * u, 0.0),
                 4: np.where(on, wv_c * u, 0.0)}
        for k, w in enumerate((w0, w1, np.where(on, w2, 0.0))):
            terms[5 + 3 * k], terms[6 + 3 * k], terms[7 + 3 * k] = w * r["dl"], w * r["dlx"], w * r["dn"]
        terms.update({15: np.where(on, wp_c, 0.0), 16: np.where(atC, fn * u, 0.0), 17: np.where(atC, fv * u, 0.0),
                      18: np.where(atC, fp, 0.0)})
        w3 = params["bmC"] * fv
        for k, w in enumerate((fn, w3)):
            terms[19 + 3 * k], terms[20 + 3 * k], terms[21 + 3 * k] = (np.where(atC, w * r[q], 0.0)
                                                                       for q in ("dl", "dlx", "dn"))
        mb_on = ~atH
        for c, j in enumerate(be.state.get("cols", [])):
            mb = np.where(mb_on, af * case.Y[case.face, j], 0.0)
            terms[1000 + c] = wn_h * u * mb
            terms[2000 + c] = (wp_h + params["bmH"] * wv_h * u) * mb
        names = [e[0] for e in plan]
        if cf != any(n.endswith("f") for n in names):
            raise AssertionError("plan %s for cf=%s" % (names, cf))
        return np.array([np.bincount(case.pt, weights=terms[k], minlength=case.npt) for k in codes])
    be._run_fused = run_fused
    be._run_sums = lambda codes, params: run_fused(evaluator.CudaBackend._fused_plan(be, codes), codes, params)
    be.PhiCO2 = None
    return be


def test_co2_assembly():
    """positive_activity._evaluate_co2_floor (the host part of the CO2 floor) with emulated kernel sums against the
    independent NumPy model; kappa ~ 0 (all at the floor) and kappa = 2 (none at the floor: bitwise the original)."""
    import positive_activity as pa
    from cops_cuda import backend as bk, evaluator
    from cops_cuda.errors import BackendError
    case = _co2_case(np.random.default_rng(11))
    T = np.array([180.0, 140.0, 120.0])
    req = bk.ColumnRequest(basis_index=(0, 1, 2), temperature=True, kappa=True)
    worst = {"p": 0.0, "kappa_analytic": 0.0, "fd": 0.0}
    detail = []
    for kappa in (0.3, 0.05, 1.7):
        be = _co2_emulated_backend(case, cf=True)
        st = bk.ParameterState(case.a, T, kappa, None)
        res = pa._evaluate_co2_floor(be, st, req, evaluator, bk, BackendError)
        pNG, pRG, dkNG, dkRG, nfl = _co2_reference(case, case.a, T, kappa)
        rel = lambda x, y: float(np.max(np.abs(x - y)) / np.max(np.abs(y)))    # noqa: E731
        worst["p"] = max(worst["p"], rel(res.pNG, pNG), rel(res.pRG, pRG))
        kNG, kRG = res.columns[("kappa", 0)]
        worst["kappa_analytic"] = max(worst["kappa_analytic"], rel(kNG, dkNG), rel(kRG, dkRG))
        # central differences of the independent model (a, T, kappa); kappa only where no record is within 1e-4 of
        # its CO2-floor kink (the step is 1e-6 relative)
        r0 = _co2_records(case, case.a, T)
        pos = r0["fminC"] > 0
        kink_free = bool(np.min(np.abs(kappa * r0["qc"][pos] / r0["fminC"][pos] - 1.0)) > 1e-4)
        for key, (ng, rg) in res.columns.items():
            kind, k = key
            if kind == "kappa" and not kink_free:
                continue
            h = 1e-6 * (abs(kappa) if kind == "kappa" else (100.0 if kind == "T" else 1.0))
            args = []
            for sgn in (1.0, -1.0):
                a_, T_, k_ = case.a.copy(), T.copy(), kappa
                if kind == "a":
                    a_[k] += sgn * h
                elif kind == "T":
                    T_[k] += sgn * h
                else:
                    k_ += sgn * h
                args.append(_co2_reference(case, a_, T_, k_)[:2])
            fdNG, fdRG = ((args[0][i] - args[1][i]) / (2 * h) for i in (0, 1))
            worst["fd"] = max(worst["fd"], rel(ng, fdNG), rel(rg, fdRG))
        detail.append("kappa %g: %d of %d records at the CO2 floor%s" % (kappa, nfl, len(case.pt),
                                                                         "" if kink_free else " (kappa FD skipped)"))
    ok = worst["p"] < 1e-13 and worst["kappa_analytic"] < 1e-13 and worst["fd"] < 1e-6
    report("CO2 floor: host assembly vs independent model (emulated sums)", ok,
           "max rel: pressures %.1e, kappa column vs analytic %.1e, columns vs central FD %.1e; %s"
           % (worst["p"], worst["kappa_analytic"], worst["fd"], "; ".join(detail)))
    # kappa ~ 0: every record with ACO2 > 0 at the CO2 floor: kappa column exactly 0
    be = _co2_emulated_backend(case, cf=True)
    res = pa._evaluate_co2_floor(be, bk.ParameterState(case.a, T, 4.04e-19, None), req, evaluator, bk, BackendError)
    pNG, pRG, _, _, nfl = _co2_reference(case, case.a, T, 4.04e-19)
    zero = not np.any(res.columns[("kappa", 0)][0]) and not np.any(res.columns[("kappa", 0)][1])
    kslot = be.dataset.arrays["mesh.ACO2"].a[2 * case.nf]
    ok0 = zero and kslot == 4.04e-19 and float(np.max(np.abs(res.pNG / pNG - 1))) < 1e-13
    report("CO2 floor: kappa ~ 0 -> kappa column exactly 0, kappa slot written", ok0,
           "%d of %d records at the floor (%d with ACO2 = 0); kappa slot %r" % (
               nfl, len(case.pt), int(np.sum(case.aco2[case.face] == 0)), kslot))
    # kappa = 2 > 1: no record at the CO2 floor -> bitwise the original CudaBackend._evaluate (no-floor sums)
    be_cf = _co2_emulated_backend(case, cf=True)
    st = bk.ParameterState(case.a, T, 2.0, None)
    r_cf = pa._evaluate_co2_floor(be_cf, st, req, evaluator, bk, BackendError)
    be_or = _co2_emulated_backend(case, cf=False)
    be_or.PhiCO2 = np.bincount(case.pt, weights=_co2_records(case, case.a, T)["qc"] * case.kphi, minlength=case.npt)
    r_or = evaluator.CudaBackend._evaluate(be_or, st, req)
    same = (r_cf.pNG.tobytes() == r_or.pNG.tobytes() and r_cf.pRG.tobytes() == r_or.pRG.tobytes() and
            list(r_cf.columns) == list(r_or.columns) and
            all(np.asarray(r_cf.columns[k][i]).tobytes() == np.asarray(r_or.columns[k][i]).tobytes()
                for k in r_cf.columns for i in (0, 1)))
    report("CO2 floor: kappa 2 (none at the floor) == original _evaluate, bitwise", same,
           "columns %s" % list(r_cf.columns))


def test_e2e_co2(scratch):
    """--ratio 1e-5 (the CO2 floor decision) from the lambda-0.3 restart and its resume, CPU stand-in (it ignores the
    GPU environment: here only the recording and the start-check logic are tested)."""
    root = Path(scratch) / "e2e"
    run, bundle, points = root / "base", root / "templates.npz", root / "late_grid_points.txt"
    restart = root / "fit_lam0.3" / "fit_restart.npz"
    if not restart.is_file():
        report("e2e CO2 floor: source restart", False, "missing %s (test_e2e first)" % restart)
        return
    rc, out, text = run_child(root, "start_r1e-5", run, bundle, points, 0.3, ["--start-restart", str(restart),
                                                                             "--ratio", "1e-5"])
    ok = rc == 0
    if ok:
        pos = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]
        sc = pos["start_check"]
        rmeta = json.loads(str(np.load(out / "fit_restart.npz")["meta_json"]))["model_extension"]
        txt = (out / "ACTIVITY_PARAMETERIZATION.txt").read_text()
        z = np.load(out / "activity_exp_facets.npz")
        zmeta = json.loads(str(z["meta_json"]))
        ok = (not sc["reproduction_expected"] and "failed" not in sc and "CO2 floor off -> on" in sc["changes"][0]
              and sc["start"]["co2_floor"] is True and rmeta["co2_floor"] is True and rmeta["night_ratio"] == 1e-5
              and pos["ratio"]["co2_floor"] is True and "use the same floor for CO2" in pos["ratio"]["decision"]
              and "CO2 flux floor (COPS_POS_FLOOR_CO2=1)" in txt and "H2O floor = CO2 floor" in txt
              and "F_min_CO2_over_qmax_CO2" in z.files and zmeta.get("co2_floor") is True
              and "use the same floor for CO2" in rmeta["policy"])
        detail = "changes %s; marker co2_floor %s; npz keys %s" % (sc["changes"], rmeta["co2_floor"], z.files)
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e CO2 floor: ratio 1e-5 start (change reported, recorded everywhere)", ok, detail)
    rs1 = root / "start_r1e-5" / "fit_restart.npz"
    rc, out, text = run_child(root, "resume_r1e-5", run, bundle, points, 0.3, ["--restart-from", str(rs1),
                                                                              "--ratio", "1e-5"])
    ok = rc == 0
    if ok:
        sc = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]["start_check"]
        ok = sc["reproduction_expected"] and "failed" not in sc and abs(sc["cost_rel_diff"]) <= 1e-12
        detail = "rel %.1e" % sc["cost_rel_diff"]
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e CO2 floor: resume of the 1e-5 fit reproduces", ok, detail)
    rc, out, text = run_child(root, "resume_r1e-5_as_0002", run, bundle, points, 0.3, ["--restart-from", str(rs1),
                                                                                      "--ratio", "0.002"])
    report("e2e CO2 floor: resume of the 1e-5 fit with ratio 0.002 is refused",
           rc != 0 and "a resume must continue the same problem" in text, "exit %d" % rc)



def test_start_set_arguments(scratch):
    import positive_problem as pp
    table = Path(scratch) / "pts_ss.txt"
    table.write_text("# index UTC kind original_index L2_start L2_stop\n0 2016-09-30T10:33:10 obs 1 a b\n")
    rs = Path(scratch) / "fake_restart.npz"
    np.savez(rs, x=np.zeros(1))
    tail = ["--", "/r/base", "--activity-template", "B", "cliff_ge_60deg"]
    cfg, argv, notes = pp.rewrite_arguments(["--point-table", str(table), "--start-restart", str(rs)] + tail)
    ok = "start_set" not in cfg and not any("start-set" in n for n in notes)
    cfg2, argv2, notes2 = pp.rewrite_arguments(["--point-table", str(table), "--start-restart", str(rs),
                                                "--start-set", "kappa=1"] + tail)
    ok &= cfg2["start_set"] == [["kappa", 1.0]] and "--start-set" not in argv2 and argv2 == argv
    bad = []
    for extra in (["--start-set", "kappa=1"], ["--start-restart", str(rs), "--start-set", "kappa"],
                  ["--start-restart", str(rs), "--start-set", "kappa=x"], ["--start-restart", str(rs), "--start-set",
                                                                           "kappa=1,kappa=2"],
                  ["--start-restart", str(rs), "--start-set", "kappa=inf"],
                  ["--start-restart", str(rs), "--grain-start-uPa", "5", "--start-set", "grain_A_uPa=4"]):
        try:
            pp.rewrite_arguments(["--point-table", str(table)] + extra + tail)
            bad.append(False)
        except SystemExit:
            bad.append(True)
    ok &= all(bad)
    report("start-set: arguments (absent = before; refusals)", ok, "parsed %s; refused %s" % (cfg2["start_set"], bad))


def test_e2e_start_set(scratch):
    root = Path(scratch) / "e2e"
    run, bundle, points = root / "base", root / "templates.npz", root / "late_grid_points.txt"
    restart = root / "fit_lam0.3" / "fit_restart.npz"
    if not restart.is_file():
        report("e2e start-set: source restart", False, "missing %s (test_e2e first)" % restart)
        return
    src = np.load(restart)
    src_x = np.asarray(src["x"], float)
    src_names = json.loads(str(src["meta_json"]))["names"]
    kk = src_names.index("kappa")
    rc, out, text = run_child(root, "start_set_k1", run, bundle, points, 0.3, ["--start-restart", str(restart),
                                                                              "--start-set", "kappa=1"])
    ok = rc == 0
    if ok:
        man = json.loads((out / "fit_manifest.json").read_text())
        pos = man["model_extension"]["positive"]
        sc = pos["start_check"]
        rmeta = json.loads(str(np.load(out / "fit_restart.npz")["meta_json"]))["model_extension"]
        txt = (out / "ACTIVITY_PARAMETERIZATION.txt").read_text()
        a = pos["start_set"]["assignments"][0]
        ok = (not sc["reproduction_expected"] and "failed" not in sc and sc["same"].get("start_set") is False
              and any(c.startswith("--start-set kappa=1.0") for c in sc["changes"])
              and sc["start_set"][0]["start_value"] == 1.0 and sc["start_set"][0]["source_value"] == src_x[kk]
              and a["name"] == "kappa" and a["value"] == 1.0 and a["start_restart_value"] == src_x[kk]
              and "kappa = 1 where stuck" in a["decision"]
              and rmeta["start_set"]["assignments"][0]["value"] == 1.0
              and rmeta["fit_chain"]["start_set"][0]["value"] == 1.0
              and "then --start-set kappa=1.0" in txt
              and "the start is NOT the source's state: --start-set kappa = 1.0" in text
              and pos["config"]["start_set"] == [["kappa", 1.0]])
        detail = "changes %s; start kappa %r (source %r)" % (sc["changes"], sc["start_set"][0]["start_value"], src_x[kk])
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e start-set: kappa=1 applied after --start-restart, recorded everywhere", ok, detail)
    for line in start_lines(text):
        print("   | " + line[:240])
    rs1 = root / "start_set_k1" / "fit_restart.npz"
    rc, out, text = run_child(root, "start_set_resume", run, bundle, points, 0.3, ["--restart-from", str(rs1)])
    ok = rc == 0
    if ok:
        sc = json.loads((out / "fit_manifest.json").read_text())["model_extension"]["positive"]["start_check"]
        rmeta = json.loads(str(np.load(out / "fit_restart.npz")["meta_json"]))["model_extension"]
        txt = (out / "ACTIVITY_PARAMETERIZATION.txt").read_text()
        ok = (sc["reproduction_expected"] and "failed" not in sc and abs(sc["cost_rel_diff"]) <= 1e-12 and
              rmeta["fit_chain"]["start_set"][0]["value"] == 1.0 and "then --start-set kappa=1.0" in txt and
              "start_set" not in sc["same"])
        detail = "rel %.1e; chain start_set %s" % (sc["cost_rel_diff"], rmeta["fit_chain"].get("start_set"))
    else:
        detail = "exit %d\n%s" % (rc, text[-2500:])
    report("e2e start-set: resume keeps the chain's start_set, reproduces", ok, detail)
    res = []
    for name, extra in (("start_set_unknown", ["--start-set", "kapa=1"]), ("start_set_bounds", ["--start-set",
                                                                                              "kappa=-1"])):
        rc, out, text = run_child(root, name, run, bundle, points, 0.3, ["--start-restart", str(restart)] + extra)
        res.append(rc != 0 and ("not a parameter of this fit" in text or "outside its bounds" in text))
    report("e2e start-set: unknown name and out-of-bounds value refused", all(res), str(res))
    # the fits without the option: no start_set recorded anywhere
    clean = []
    for name in ("start_same", "fit_lam0.3", "resume_same"):
        d = root / name
        if not (d / "fit_manifest.json").is_file():
            clean.append("%s missing" % name)
            continue
        pos = json.loads((d / "fit_manifest.json").read_text())["model_extension"]["positive"]
        rmeta = json.loads(str(np.load(d / "fit_restart.npz")["meta_json"]))["model_extension"]
        clean.append("start_set" not in pos and "start_set" not in (pos.get("config") or {}) and
                     "start_set" not in rmeta and "start_set" not in (rmeta.get("fit_chain") or {}) and
                     "start_set" not in ((pos.get("start_check") or {}).get("same") or {}) and
                     "--start-set" not in (d / "ACTIVITY_PARAMETERIZATION.txt").read_text())
    report("e2e start-set: fits without the option record nothing new", all(c is True for c in clean), str(clean))


def test_12m(out_json):
    """Start projection (emission-weighted) and log Tikhonov on the 12M base (elements.bin of shap8_012m_flat_earlier; the lategrid base
    has the same mesh) with cliff_ge_60deg of restart_templates.npz."""
    import time
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    at = reference.load("activity_templates")
    import positive_problem as pp
    R = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_flat_earlier")
    t0 = time.time()
    el = fit_cops.read_elements(str(R / "base" / "elements.bin"))
    bundle = at.load_bundle(str(R / "bundle" / "restart_templates.npz"), el["centroid"], el["area"])
    T = np.ascontiguousarray(at.select_templates(bundle, ["cliff_ge_60deg"]))
    print("loaded %d facets in %.0f s" % (len(el["area"]), time.time() - t0), flush=True)
    header, _, _ = fit_cops.read_model_file(str(R / "base" / "COPS_model.dat"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    e = float(np.max(np.abs(np.maximum(el["Y"] @ coef[0], 0.0) - el["AH2O"])) / np.max(el["AH2O"]))
    grams = pp.gram_matrices(el["Y"], el["area"], el["AH2O"], T)
    proj = pp.start_projection(el["Y"], el["area"], el["AH2O"], T, grams)
    Lt, G, info = pp.log_tikhonov(grams[0], float(np.sum(el["area"])), proj["ascale"], 25)
    x0 = np.concatenate([proj["a"] / proj["ascale"], proj["c"]])
    rng = np.random.default_rng(3)
    dx = rng.normal(0, 1e-3, 26)
    ds = pp.exponent_field(el["Y"], dx[:25] * proj["ascale"], T, dx[25:])
    D_direct = math.sqrt(float(np.sum(el["area"] * ds ** 2) / np.sum(el["area"])))
    D_op = float(np.linalg.norm(Lt @ dx))
    res = {"facets": int(len(el["area"])), "A_orig_check_rel": e, "seconds": time.time() - t0,
           "projection": proj["stats"], "ascale": proj["ascale"], "a": [float(v) for v in proj["a"]],
           "c_cliff": [float(v) for v in proj["c"]], "x0_scaled": [float(v) for v in x0],
           "tikhonov": info, "D_operator_over_direct_minus_1": D_op / D_direct - 1.0,
           "cliff_area_fraction": float(np.sum(el["area"] * T[:, 0]) / np.sum(el["area"]))}
    res["alternatives"] = projection_alternatives(el["Y"], el["area"], el["AH2O"], T, pp)
    ind = np.array(res["alternatives"]["activity_weighted_log_ls"]["theta"])
    mod = np.concatenate([proj["a"], proj["c"]])
    res["module_vs_independent_emission_projection_max_rel"] = float(np.max(np.abs(mod - ind)) / np.max(np.abs(ind)))
    print(json.dumps(res, indent=1), flush=True)
    if out_json:
        Path(out_json).write_text(json.dumps(res, indent=1) + "\n")
    ok = (info["LtT_Lt_minus_G_rel"] < 1e-12 and abs(D_op / D_direct - 1.0) < 1e-8 and e < 1e-12 and
          res["module_vs_independent_emission_projection_max_rel"] < 1e-8 and
          proj["stats"]["weighting"] == "emission_weighted_log_ls")
    report("12M: projection, Lt^T Lt == G, D == rms(s - s0)", ok, "see %s" % out_json)


def projection_alternatives(Y, area, A, T, pp, chunk=1 << 20):
    """Diagnostics: independent computations of other maps of A_orig onto exp([Y, T] theta).
    activity_weighted_log_ls: weights area x A_orig over A_orig > 0 (unnormalized; = the module's emission-weighted
    projection, the user's rule since 2026-11-07: cross-check); area_weighted_log_ls: weights area (the first rule);
    linear_space_nls: min sum area (exp(s) - A)^2 over all facets (Gauss-Newton from the activity-weighted solution)."""
    n, m = len(area), Y.shape[1] + T.shape[1]

    def blocks():
        for i in range(0, n, chunk):
            sl = slice(i, min(i + chunk, n))
            yield sl, np.column_stack((Y[sl], T[sl]))

    def field(theta):
        s = np.empty(n)
        for sl, B in blocks():
            s[sl] = B @ theta
        return s
    out = {}
    G, r = np.zeros((m, m)), np.zeros(m)
    for sl, B in blocks():
        w = np.where(A[sl] > 0, area[sl] * A[sl], 0.0)
        G += B.T @ (w[:, None] * B)
        r += B.T @ (w * np.log(np.where(A[sl] > 0, A[sl], 1.0)))
    th = np.linalg.solve(G, r)
    s = field(th)
    out["activity_weighted_log_ls"] = dict(pp.projection_stats(s, area, A), c_cliff=float(th[-1]),
                                           theta=[float(v) for v in th])
    Ga, ra = np.zeros((m, m)), np.zeros(m)
    for sl, B in blocks():
        w = np.where(A[sl] > 0, area[sl], 0.0)
        Ga += B.T @ (w[:, None] * B)
        ra += B.T @ (w * np.log(np.where(A[sl] > 0, A[sl], 1.0)))
    tha = np.linalg.solve(Ga, ra)
    out["area_weighted_log_ls"] = dict(pp.projection_stats(field(tha), area, A), c_cliff=float(tha[-1]))
    cost = float(np.sum(area * (np.exp(s) - A) ** 2))
    it = 0
    for it in range(60):
        H, g = np.zeros((m, m)), np.zeros(m)
        for sl, B in blocks():
            e = np.exp(s[sl])
            H += B.T @ ((area[sl] * e * e)[:, None] * B)
            g += B.T @ (area[sl] * e * (A[sl] - e))
        d = np.linalg.solve(H, g)
        step = 1.0
        while step > 1e-6:
            s_try = field(th + step * d)
            c_try = float(np.sum(area * (np.exp(s_try) - A) ** 2))
            if c_try < cost:
                break
            step *= 0.5
        if step <= 1e-6 or (cost - c_try) < 1e-12 * cost:
            break
        th, s, cost = th + step * d, s_try, c_try
    out["linear_space_nls"] = dict(pp.projection_stats(s, area, A), c_cliff=float(th[-1]), iterations=it + 1,
                                   relative_rms_activity_residual=float(np.sqrt(cost / np.sum(area * A ** 2))))
    return out


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--e2e-child":
        # internal: --e2e-child run bundle points out lambda [max_nfev [driver options ...]]
        child = sys.argv[2:]
        run, bundle, points, out, lam = child[:5]
        max_nfev = int(child[5]) if len(child) > 5 else 4
        e2e_child(run, bundle, points, out, float(lam), child[6:], max_nfev)
        return 0
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=os.environ.get("POSITIVE_TEST_SCRATCH"))
    ap.add_argument("--12m", dest="m12", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.m12:
        test_12m(a.out)
    else:
        if not a.scratch:
            sys.exit("give --scratch DIR (a new directory for the synthetic base)")
        Path(a.scratch).mkdir(parents=True, exist_ok=False)
        test_window()
        test_projection()
        test_arguments(a.scratch)
        test_selection()
        test_held_units()
        test_co2_environment(a.scratch)
        test_co2_assembly()
        test_e2e(a.scratch)
        test_e2e_start(a.scratch)
        test_e2e_co2(a.scratch)
        test_start_set_arguments(a.scratch)
        test_e2e_start_set(a.scratch)
        test_e2e_held(a.scratch)
    bad = [r for r in RESULTS if not r[1]]
    print("\n%d tests, %d failed" % (len(RESULTS), len(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

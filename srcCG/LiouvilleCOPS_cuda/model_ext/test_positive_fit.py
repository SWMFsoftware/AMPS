#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_positive_fit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_positive_fit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: eb52c3990e70e3dcc2a305cec4866e4bbb34f8d9c56c091b29bba39c94f6d623
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. Unit tests of positive_problem.py (projection, Tikhonov, W10x5 operator, grid rows, Jacobian). user decisions 2026-11-07.
# USED BY: jobs/positive_unit_12m.sbatch.
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
    ok = (cfg["check_jacobian"] and cfg["grain_start_uPa"] == pp.DEFAULT_GRAIN_START_UPA and
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
            Tmax, T75, Tmin = T
            act_lin = self.Y @ a
            if b is not None:
                act_lin = act_lin + enhancement_scale * (enhancement_basis @ b)
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


def e2e_child(run, bundle, points, out, lam):
    """One fit_cops.main() with every layer of positive_fit.py (child process)."""
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    import grain_crossing
    import positive_problem as pp
    pp.set_positive_environment()
    fit_cops.Model = exp_model_class(fit_cops)
    cfg, argv, notes = pp.rewrite_arguments(
        ["--point-table", str(points), "--check-jacobian", "--check-parameters",
         "a00,a02,a07,c_cliff_ge_60deg,kappa,grain_A_uPa,grain_B_uPa_per_100s", "--",
         str(run), "--allow-legacy-base", "--degree", "4", "--enhancement-degree", "4", "--lambda", str(lam),
         "--jobs", "1", "--tbounds", "30", "400", "--max-nfev", "4", "--no-save-restart", "--out", str(out),
         "--activity-template", str(bundle), "cliff_ge_60deg"])
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


# ------------------------------------------------------------------------------------------------ 12M
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=os.environ.get("POSITIVE_TEST_SCRATCH"))
    ap.add_argument("--12m", dest="m12", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--e2e-child", nargs=5, default=None)
    a = ap.parse_args()
    if a.e2e_child:
        run, bundle, points, out, lam = a.e2e_child
        e2e_child(run, bundle, points, out, float(lam))
        return 0
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
        test_e2e(a.scratch)
    bad = [r for r in RESULTS if not r[1]]
    print("\n%d tests, %d failed" % (len(RESULTS), len(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

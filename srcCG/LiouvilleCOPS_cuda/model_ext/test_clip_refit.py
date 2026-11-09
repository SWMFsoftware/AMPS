#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_clip_refit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_clip_refit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1034089955e6d68ca1007db15a62ca424c85ede19f46ebf73621f6c1519c662c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. CPU tests of clip_refit.py / clip_problem.py (14 tests; option off == grain_fit_fixed path bitwise). 2026-11-08.
# ==== end of header ====
"""CPU tests of model_ext/clip_refit.py + clip_problem.py (refit of the clipped c6 model on the lategrid base, D1/D2;
2026-11-08).

    .venv/bin/python -W ignore test_clip_refit.py --scratch DIR      (login node: ~1-2 min, < 2 GB)

 1. arguments: rewrite_arguments (variants, inserted --template-start 0, refusals);
 2. real tables (COPS_model.dat only, no kernels): the all_available profile on the NEW base (1551 rows) through
    positive_problem's selection wrapper == the unwrapped profile on the OLD base (1150 rows, c6's) for every
    observation row (use flags, reasons; matched by UTC); grid rows use 0 / late_grid_point; 934 NG / 1150 RG;
    references and inside flags of the observation rows equal; the new base's LateLayout (grid rows, late observations,
    STOP offsets);
 3. end to end on synthetic bases (test_positive_fit.make_base: 400 facets, observations around the NG gap, a 1-s late
    grid) with the REAL driver chain clip_refit.main -> grain_fit_fixed.main -> grain_fit.main -> gpu_fit.main; only
    gpu_fit.drive is replaced by a CPU stand-in (the reference fit_cops.Model; gpu_fit's restart, template-start,
    segment and restart-writer steps); each fit in a subprocess:
    a. option off == the grain_fit_fixed path, bitwise: on a base WITHOUT grid rows, grain_fit_fixed.py started with
       gpu_fit --restart-from (c6's own mechanism) vs clip_refit.py --variant D1 --start-restart (same restart; + start
       check + full Jacobian check): x, cost, COPS_model.dat, grain_term.txt, ActivityFit.txt identical;
    b. D1 on the base WITH grid rows (same observation rows and kernels + the grid): x and cost bitwise equal to (a)
       (grid rows enter neither the cost nor the crossing pairs nor the grain term); the start check reproduces the
       source's final cost; grid rows flagged in COPS_model.dat and listed in late_window_model.dat, absent from
       pressure_residuals.dat and the figure tables; positive_table.py output rows == (a)'s table rows;
    c. D2 on the grid base: full Jacobian check (data, late window, crossing); COPS_model.dat late rows ==
       late_response_eval W10x5(instantaneous surface) + grain(t_obs), other rows instantaneous; the start model: rows
       before 10:33:10 bitwise equal to D1's, late rows == W @ D2's own grid rows;
    d. resume: D2 segment 2 from segment 1's restart (gpu_fit --restart-from) reproduces its cost (start check
       CONFIRMED); a D1 restart is refused for D2; held values different from the source's are refused.
"""

import argparse
import json
import os
import shutil
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
os.environ.setdefault("COPS_DATA_ROOT", "/data001/highresolu_vt/cops_cluster_transfer_20261026/"
                                        "rosetta_data_end_of_mission_cops")
OUT_ROOT = Path("/data001/highresolu_vt/cops_output")
NEW_RUN = OUT_ROOT / "runs/shap8_012m_lategrid_dz1.5"
OLD_BASE = OUT_ROOT / "runs/shap8_012m_flat_earlier/base"
RESULTS = []
FIX = ["--fix", "Tmax=182.09922", "--fix", "T75=136.11149", "--fix", "Tmin=133.0", "--fix", "c_ft_c1=1"]


def report(name, ok, detail):
    RESULTS.append((name, ok, detail))
    print("%-62s %s  %s" % (name, "PASS" if ok else "FAIL", detail), flush=True)


# ------------------------------------------------------------------------------------------------ 1 arguments
def test_arguments(scratch):
    import clip_problem as cp
    pt = Path(scratch) / "pt.txt"
    pt.write_text("# index UTC kind original_index L2_start L2_stop\n")
    rs = Path(scratch) / "r.npz"
    rs.write_bytes(b"x")
    tail = ["--", "/run", "--lambda", "0", "--activity-template", "b.npz", "cliff_ge_60deg"]
    cfg, argv, notes = cp.rewrite_arguments(["--variant", "D2", "--point-table", str(pt), "--start-restart", str(rs),
                                             "--check-full-jacobian"] + FIX + ["--devices", "0"] + tail)
    ok = (cfg["variant"] == "D2" and cfg["check_full_jacobian"] and cfg["start_restart"] == str(rs.resolve()) and
          argv[:argv.index("--")] == FIX + ["--devices", "0", "--template-start", "0"] and
          cfg["fixed"] == {"Tmax": 182.09922, "T75": 136.11149, "Tmin": 133.0, "c_ft_c1": 1.0} and
          cfg["templates"] == ["cliff_ge_60deg"] and cfg["run"] == "/run" and len(notes) == 1)
    report("arguments: D2 config, --template-start 0 inserted", ok, "argv %s" % argv[:12])
    bad = [(["--point-table", str(pt)] + tail, "variant"),
           (["--variant", "D2"] + tail, "point-table"),
           (["--variant", "D1", "--start-restart", str(rs), "--restart-from", str(rs)] + tail, "exclude"),
           (["--variant", "D1", "--template-start", "1"] + tail, "template-start"),
           (["--variant", "D1"] + tail + ["--fix-temperature"], "fix-temperature"),
           (["--variant", "D1"] + tail + ["--restart-from", "x"], "restart-from"),
           (["--variant", "D1"] + tail + ["--check-jacobian"], "check-jacobian")]
    fails = []
    for argv_, why in bad:
        try:
            cp.rewrite_arguments(argv_)
            fails.append(why)
        except SystemExit as exc:
            if why not in str(exc):
                fails.append("%s: %s" % (why, exc))
    report("arguments: refusals", not fails, "not refused / wrong message: %s" % fails if fails else "%d refused" % len(bad))


# ------------------------------------------------------------------------------------------------ 2 real tables
def select(fit_cops, run, gkeys=None):
    import positive_problem as pp
    header, names, rows = fit_cops.read_model_file(str(Path(run) / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    srcdir = str(REF)
    _, iv = fit_cops.run_intervals(header, srcdir, str(run))
    sat = fit_cops.ng_saturated(fit_cops.product_input("data/COPS_HK_20160930.txt"), [r[cols["UTC"]] for r in rows])
    profile, _ = fit_cops.load_profile("all_available", fit_cops.profile_directory(srcdir))
    original = fit_cops.apply_selection_profile
    rec = {}
    try:
        if gkeys is not None:
            pp.install_selection(fit_cops, gkeys, rec)
        res, rep = fit_cops.apply_selection_profile(rows, cols, profile, iv, sat)
    finally:
        fit_cops.apply_selection_profile = original
    return rows, res, rep, cols, rec


def test_real_tables():
    import positive_problem as pp
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    table = pp.read_point_table(NEW_RUN / "late_grid_points.txt")
    gkeys = pp.grid_keys(table)
    rows_o, res_o, rep_o, co, _ = select(fit_cops, OLD_BASE)
    rows_n, res_n, rep_n, cn, rec = select(fit_cops, NEW_RUN / "base", gkeys)
    old = {pp.time_key(r[co["UTC"]]): (j, r) for j, r in enumerate(res_o)}
    bad, ngrid, nobs = [], 0, 0
    for j, r in enumerate(res_n):
        key = pp.time_key(r[cn["UTC"]])
        if key in gkeys:
            ngrid += 1
            if (r[cn["useNG"]], r[cn["useRG"]]) != ("0", "0") or rep_n["reasons"]["NG"][j] != pp.GRID_REASON \
                    or rep_n["reasons"]["RG"][j] != pp.GRID_REASON or float(r[cn["refNG"]]) != 0 \
                    or float(r[cn["refRG"]]) != 0:
                bad.append("grid %s" % r[cn["UTC"]])
            continue
        nobs += 1
        if key not in old:
            bad.append("obs %s not in the old base" % r[cn["UTC"]])
            continue
        jo, ro = old[key]
        same = all(r[cn[c]] == ro[co[c]] for c in ("useNG", "useRG", "refNG", "refRG", "inside"))
        same &= rep_n["reasons"]["NG"][j] == rep_o["reasons"]["NG"][jo] and rep_n["reasons"]["RG"][j] == rep_o["reasons"]["RG"][jo]
        if not same:
            bad.append("obs %s differs" % r[cn["UTC"]])
    used_ng = sum(int(r[cn["useNG"]]) for r in res_n)
    used_rg = sum(int(r[cn["useRG"]]) for r in res_n)
    ok = (not bad and ngrid == 401 and nobs == 1150 == len(res_o) and used_ng == 934 and used_rg == 1150 and
          rep_n["reason_sha256_observation_rows"] == rep_o["reason_sha256"])
    report("real tables: profile on new-base obs rows == old base (c6)", ok,
           "%d obs rows equal (use, ref, inside, reasons), %d grid rows use 0 / %s; used NG %d RG %d; reason sha "
           "obs-rows %s; differences %s" % (nobs, ngrid, pp.GRID_REASON, used_ng, used_rg,
                                           rep_n["reason_sha256_observation_rows"] == rep_o["reason_sha256"], bad[:3]))
    points = pp.read_run_points(NEW_RUN / "base")
    ipoint = np.array([int(r[cn["iPoint"]]) for r in rows_n])
    lay = pp.LateLayout(table, points, ipoint)
    s = lay.summary
    ok = (s["grid_rows"] == 401 and s["first_late_utc"] == "2016-09-30T10:33:10" and s["grid_step_s"] == [1.0, 1.0]
          and s["table_index_equals_row"] and abs(s["operator"]["row_sum_range"][0] - 1) < 1e-12
          and abs(s["operator"]["row_sum_range"][1] - 1) < 1e-12)
    report("real tables: LateLayout of the new base", ok,
           "late obs %d (%s .. %s), grid %d (%s .. %s), STOP - time %s s, W %s nnz %d" % (
               s["late_observations"], s["first_late_utc"], s["last_late_utc"], s["grid_rows"], s["grid_first_utc"],
               s["grid_last_utc"], ["%.3f" % v for v in s["stop_minus_time_s"]], s["operator"]["shape"],
               s["operator"]["nnz"]))


# ------------------------------------------------------------------------------------------------ 3 end to end
def make_bases(root):
    import test_positive_fit as tpf
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    at = reference.load("activity_templates")
    gc, gt = reference.load("gravity_cache"), reference.load("gravity_topology")
    root = Path(root)
    run_g, bundle, points = tpf.make_base(root)
    el = fit_cops.read_elements(str(run_g / "elements.bin"))
    rng = np.random.default_rng(11)
    n = len(el["area"])
    col = np.where(rng.random(n) < 0.1, rng.uniform(0.05, 0.3, n), 0.0)
    ft = root / "ft_c1.templates.npz"
    at.write_bundle(str(ft), ["ft_c1"], col[:, None], el["centroid"], el["area"],
                    {"templates": {"ft_c1": {"definition": "synthetic fine-tuning template (test)"}},
                     "gravity": {"cache": {"definition": {"algorithm": gc.ALGORITHM,
                                                          "orientation_policy": gt.ORIENTATION_POLICY}},
                                 "integrity_check": "DISABLE after topology audit"}})
    # the base without grid rows: the observation rows of run_g (same iPoints, kernels, elements)
    import positive_problem as pp
    table = pp.read_point_table(points)
    grid_utc = {e["utc"] for e in table if e["kind"] == "grid"}
    run_n = root / "base_nogrid"
    (run_n / "kernels").mkdir(parents=True)
    shutil.copy2(run_g / "elements.bin", run_n / "elements.bin")
    shutil.copy2(run_g / "input.snapshot.in", run_n / "input.snapshot.in")
    keep = []
    for line in open(run_g / "COPS_model.dat"):
        if line.startswith("#"):
            keep.append(line)
            continue
        v = line.split()
        if v[1] in grid_utc:
            continue
        keep.append(line)
        shutil.copy2(run_g / "kernels" / ("point=%04d.bin" % int(v[0])), run_n / "kernels" / ("point=%04d.bin" % int(v[0])))
    (run_n / "COPS_model.dat").write_text("".join(keep))
    return run_g, run_n, bundle, ft, points


def fit_args(run, out, bundle, ft):
    return [str(run), "--allow-legacy-base", "--degree", "4", "--enhancement-degree", "4", "--lambda", "0", "--jobs",
            "1", "--tbounds", "30", "400", "--ftol", "1e-10", "--xtol", "1e-10", "--gtol", "1e-10",
            "--check-tolerance", "inf", "--out", str(out), "--activity-template", str(bundle), "cliff_ge_60deg",
            "--activity-template", str(ft), "ft_c1"]


def cpu_drive(args, argv, fit_args_, devices, max_nfev, template_start, comm, deadline=None):
    """gpu_fit.drive without a GPU: the reference CPU fit_cops.Model; gpu_fit's restart (--restart-from: templates,
    mesh check, install_template_start), segment and restart-writer steps; the figure tables (cops_cuda.figures)."""
    import gpu_fit
    from cops_cuda import provenance, reference
    from cops_cuda import restart as rs
    reference.bind()
    fit_cops = reference.load("fit_cops")
    restart_state, restart_info = None, {"save": True}
    if args.restart_from:
        restart_state = rs.read(args.restart_from)
        restart_info["loaded"] = provenance.file_identity(args.restart_from, "fit_restart")
        fit_args_ = gpu_fit.restart_template_args(restart_state, fit_args_, restart_info)
        gpu_fit.check_restart_mesh(restart_state, fit_args_, fit_cops)
    fit_start = {"template_start": template_start, "applied": False, "max_nfev": max_nfev}
    original_save = fit_cops.save_active_manifest

    def save():
        if fit_cops.ACTIVE_MANIFEST is not None:
            fit_cops.ACTIVE_MANIFEST["compute_backend"] = {"kind": "CPU test stand-in for gpu_fit.drive",
                                                           "fit_start": fit_start, "restart": restart_info,
                                                           "driver_command": list(argv)}
        original_save()
    fit_cops.save_active_manifest = save
    captured = []
    gpu_fit.install_template_start(fit_cops, template_start, fit_start, restart_state, captured)
    previous = list(restart_state.meta.get("segments", [])) if restart_state is not None else []
    segment = args.segment if args.segment is not None else len(previous) + 1
    sys.argv = [str(Path(fit_cops.__file__))] + fit_args_
    fit_cops.main()
    out = Path(fit_cops.ACTIVE_MANIFEST_PATH).parent
    elements = fit_cops.read_elements(fit_cops.ACTIVE_MANIFEST["inputs"]["base_elements"]["path"])
    state = rs.from_fit(fit_cops.P, captured[-1], out, fit_cops.ACTIVE_MANIFEST, elements)
    if restart_state is not None:
        state.meta["history"] = list(restart_state.meta.get("history", [])) + state.meta["history"]
        state.meta["started_from"] = restart_info["loaded"]
    r = captured[-1]
    state.meta["segments"] = previous + [{"index": int(segment), "status": "complete", "nfev": int(r["nfev"]),
                                          "max_nfev": int(max_nfev), "cost": float(r["cost"])}]
    state.meta["status"] = "complete"
    rs.write(out / "fit_restart.npz", state)
    from cops_cuda import figures
    base_dir = Path(fit_cops.ACTIVE_MANIFEST["inputs"]["base_elements"]["path"]).parent
    (out / "test_figure_rows.json").write_text(json.dumps([len(figures.read_table(out)["utc"]),
                                                           len(figures.read_table(base_dir)["utc"])]))
    save()
    return 0


def child(spec):
    import gpu_fit
    gpu_fit.drive = cpu_drive
    if spec["mode"] == "driver":
        import clip_refit
        sys.argv = [str(HERE / "clip_refit.py")] + spec["argv"]
        return clip_refit.main()
    import grain_fit_fixed
    sys.argv = [str(HERE / "grain_fit_fixed.py")] + spec["argv"]
    return grain_fit_fixed.main()


def run_child(root, name, mode, argv):
    log = Path(root) / (name + ".log")
    spec = json.dumps({"mode": mode, "argv": [str(a) for a in argv]})
    rc = subprocess.call([sys.executable, "-W", "ignore", __file__, "--child", spec], stdout=log.open("w"),
                         stderr=subprocess.STDOUT)
    return rc, log.read_text()


def restart_arrays(path):
    with np.load(str(path), allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        return {k: np.array(z[k]) for k in z.files if k != "meta_json"}, meta


def table_rows(path):
    return [l for l in open(path) if l.strip() and not l.startswith("#")]


def lines_of(text, prefix="clip_refit:"):
    return [l for l in text.splitlines() if l.startswith(prefix) or "Error" in l or "Traceback" in l]


def test_e2e(scratch):
    import positive_problem as pp
    root = Path(scratch) / "e2e"
    run_g, run_n, bundle, ft, points = make_bases(root)
    gff_head = FIX + ["--crossing-weight", "10", "--devices", "0", "--template-start", "0"]
    # source ("c6" analog): grain_fit_fixed on the no-grid base from the default start
    S = root / "S0"
    rc, text = run_child(root, "S0", "gff", gff_head + ["--max-nfev", "8", "--"] + fit_args(run_n, S, bundle, ft))
    if rc != 0:
        report("e2e source fit (grain_fit_fixed, no-grid base)", False, "exit %d\n%s" % (rc, text[-3000:]))
        return
    s_arr, s_meta = restart_arrays(S / "fit_restart.npz")
    report("e2e source fit (grain_fit_fixed, no-grid base)", True, "cost %.9e nfev %s" % (
        s_meta["optimizer"]["cost"], s_meta["optimizer"]["nfev"]))
    # a. option off == grain_fit_fixed path (c6's restart mechanism), no-grid base
    R1, C1 = root / "a_gff_restart_from", root / "a_clip_D1"
    rc1, t1 = run_child(root, R1.name, "gff", gff_head + ["--max-nfev", "4", "--restart-from", S / "fit_restart.npz",
                                                          "--"] + fit_args(run_n, R1, bundle, ft))
    rc2, t2 = run_child(root, C1.name, "driver", ["--variant", "D1", "--start-restart", S / "fit_restart.npz",
                                                   "--check-full-jacobian"] + gff_head +
                        ["--max-nfev", "4", "--"] + fit_args(run_n, C1, bundle, ft))
    if rc1 or rc2:
        report("e2e a: option off == grain_fit_fixed path (bitwise)", False,
               "exit %d / %d\n%s\n%s" % (rc1, rc2, t1[-2000:], t2[-2500:]))
        return
    a1, m1 = restart_arrays(R1 / "fit_restart.npz")
    a2, m2 = restart_arrays(C1 / "fit_restart.npz")
    same_files = {f: (R1 / f).read_bytes() == (C1 / f).read_bytes()
                  for f in ("COPS_model.dat", "grain_term.txt", "ActivityFit.txt", "pressure_residuals.dat")}
    ok = (np.array_equal(a1["x"], a2["x"]) and m1["optimizer"]["cost"] == m2["optimizer"]["cost"] and
          np.array_equal(a1["pred_ng"], a2["pred_ng"]) and all(same_files.values()) and
          m1["optimizer"]["nfev"] == m2["optimizer"]["nfev"])
    man2 = json.loads((C1 / "fit_manifest.json").read_text())
    me = man2["model_extension"]
    cr = me["clip_refit"]
    ok_rec = (me["command"][0].endswith("clip_refit.py") and me["grain_fit_fixed_command"][0].endswith("grain_fit_fixed.py")
              and cr["start_restart"]["x_equals_source_x"] and cr["start_check"]["cost_rel_diff_vs_source_final"] == 0.0
              and cr["jacobian_check"]["worst_data"] < 1e-5 and m2["model_extension"]["parameterization"] ==
              "clip_refit_c6_D1_instantaneous_v1" and m2["segments"][-1]["index"] == 1 and len(m2["segments"]) == 1)
    report("e2e a: option off == grain_fit_fixed path (bitwise)", ok and ok_rec,
           "x equal %s, cost %.12e / %.12e, files %s; start cost rel %s, Jacobian worst data %.1e, segments %d (gff "
           "restart-from path carries %d)" % (np.array_equal(a1["x"], a2["x"]), m1["optimizer"]["cost"],
                                               m2["optimizer"]["cost"], same_files,
                                               cr["start_check"]["cost_rel_diff_vs_source_final"],
                                               cr["jacobian_check"]["worst_data"], len(m2["segments"]),
                                               len(m1["segments"])))
    for line in lines_of(t2)[:8]:
        print("   | " + line[:230])
    # b. D1 on the grid base
    D1 = root / "b_clip_D1_grid"
    rc, t3 = run_child(root, D1.name, "driver", ["--variant", "D1", "--point-table", points, "--start-restart",
                                                  S / "fit_restart.npz", "--check-full-jacobian"] + gff_head +
                       ["--max-nfev", "4", "--"] + fit_args(run_g, D1, bundle, ft))
    if rc:
        report("e2e b: D1 with grid rows == without (bitwise)", False, "exit %d\n%s" % (rc, t3[-3000:]))
        return
    a3, m3 = restart_arrays(D1 / "fit_restart.npz")
    man3 = json.loads((D1 / "fit_manifest.json").read_text())
    cr3 = man3["model_extension"]["clip_refit"]
    gw1, gw3 = cr["grain_window"], cr3["grain_window"]
    ok = (np.array_equal(a3["x"], a2["x"]) and m3["optimizer"]["cost"] == m2["optimizer"]["cost"] and
          cr3["start_check"]["cost_rel_diff_vs_source_final"] == 0.0 and
          all(gw1[k] == gw3[k] for k in ("t_end_et", "gap_seconds", "pairs_before_gap", "pairs_after_gap", "NNG", "NRG"))
          and cr3["late_window"]["grid_rows_in_fit_model"] == 0 and cr3["late_window"]["grid_data_entries_used"] == 0)
    report("e2e b: D1 with grid rows == without (bitwise)", ok,
           "x equal %s, cost %.12e / %.12e; start cost rel %s; grain window %s; late_window %s" % (
               np.array_equal(a3["x"], a2["x"]), m3["optimizer"]["cost"], m2["optimizer"]["cost"],
               cr3["start_check"]["cost_rel_diff_vs_source_final"],
               {k: gw3[k] for k in ("gap_seconds", "pairs_before_gap", "pairs_after_gap", "NNG", "NRG")},
               {k: cr3["late_window"].get(k) for k in ("fit_model_rows", "all_rows_model_rows", "grid_rows_in_fit_model")}))
    table = pp.read_point_table(points)
    n_obs = sum(1 for e in table if e["kind"] == "obs")
    n_grid = sum(1 for e in table if e["kind"] == "grid")
    gutc = {e["utc"] for e in table if e["kind"] == "grid"}

    def products(out):
        note = [l for l in open(out / "COPS_model.dat") if l.startswith("# ClipRefitGridRows")]
        lwm = [l.split() for l in open(out / "late_window_model.dat") if not l.startswith("#")]
        kinds = [v[3] for v in lwm]
        res_rows = [l.split() for l in open(out / "pressure_residuals.dat") if not l.startswith("#")]
        nfig = json.loads((out / "test_figure_rows.json").read_text())
        pt_out = out / "positive_table_out.dat"
        rc_ = subprocess.call([sys.executable, str(OUT_ROOT / "figures_cluster/scripts/positive_table.py"), str(out),
                               str(pt_out)], stdout=subprocess.DEVNULL)
        return note, kinds, res_rows, nfig, rc_, pt_out
    note, kinds, res_rows, nfig, rc_pt, pt_out = products(D1)
    ok = (len(note) == 1 and kinds.count("grid") == n_grid and kinds.count("late_obs") > 0 and
          len(res_rows) == 2 * n_obs and not any(v[1] in gutc for v in res_rows) and nfig == [n_obs, n_obs] and
          rc_pt == 0 and table_rows(pt_out) == table_rows(C1 / "COPS_model.dat") and
          len(table_rows(D1 / "COPS_model.dat")) == n_obs + n_grid and
          cr3["late_window_products"]["table_vs_surface_plus_grain_max_rel"] < 1e-6 and (D1 / "CLIP_REFIT.txt").is_file())
    report("e2e b: D1 grid rows flagged/listed; residuals+figures obs only", ok,
           "note %d; late_window_model grid %d late_obs %d; residual rows %d (2 x %d obs); figure tables %s; "
           "positive_table.py exit %d, rows == no-grid table %s" % (
               len(note), kinds.count("grid"), kinds.count("late_obs"), len(res_rows), n_obs, nfig, rc_pt,
               table_rows(pt_out) == table_rows(C1 / "COPS_model.dat") if pt_out.is_file() else None))
    # c. D2 on the grid base
    D2 = root / "c_clip_D2_grid"
    rc, t4 = run_child(root, D2.name, "driver", ["--variant", "D2", "--point-table", points, "--start-restart",
                                                  S / "fit_restart.npz", "--check-full-jacobian",
                                                  "--check-parameters",
                                                  "a00,a02,a07,c_cliff_ge_60deg,kappa,grain_A_uPa,grain_B_uPa_per_100s"]
                       + gff_head + ["--max-nfev", "4", "--"] + fit_args(run_g, D2, bundle, ft))
    if rc:
        report("e2e c: D2 W10x5", False, "exit %d\n%s" % (rc, t4[-3000:]))
        return
    man4 = json.loads((D2 / "fit_manifest.json").read_text())
    cr4 = man4["model_extension"]["clip_refit"]
    jc = cr4["jacobian_check"]["columns"]
    worst = max(v["data"] for v in jc.values() if isinstance(v, dict))
    worst_late = max(v.get("data_late_window", 0.0) for v in jc.values() if isinstance(v, dict))
    report("e2e c: D2 full Jacobian check (data / late window / crossing)", worst < 1e-5 and worst_late < 1e-5,
           "data %.1e late %.1e; crossing %s; rows %s" % (worst, worst_late, ["%.1e" % v.get("crossing", float("nan"))
                                                                              for v in jc.values()],
                                                           cr4["jacobian_check"]["rows"]))
    # COPS_model.dat late rows vs the independent late_response_eval definition
    import late_response_eval as lre
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    header, cn, rws = fit_cops.read_model_file(str(run_g / "COPS_model.dat"))
    el = fit_cops.read_elements(str(run_g / "elements.bin"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    Mraw = fit_cops.Model(str(run_g), header, {n: i for i, n in enumerate(cn)}, rws, el, coef)
    a4, m4 = restart_arrays(D2 / "fit_restart.npz")
    tcols = a4["template_columns"]
    b = np.asarray(a4["template_values"], float)
    pNG_i, pRG_i = Mraw.evaluate(a4["a"], a4["T"], float(a4["kappa"][0]), b=b, enhancement_basis=tcols,
                                 enhancement_scale=float(m4["enhancement_scale"]))
    gj = json.loads((D2 / "grain_term.json").read_text())
    A_, B_, t_end = gj["result"]["A_Pa"], gj["result"]["B_Pa_per_s"], gj["t_end_et"]
    et = np.array([float(r[cn.index("et")]) for r in rws])
    grain = np.where(et >= t_end, np.maximum(A_ - B_ * (et - t_end), 0.0), 0.0)
    rows4 = table_rows(D2 / "COPS_model.dat")
    ci = {n: i for i, n in enumerate(cn)}
    pNG_tab = np.array([float(l.split()[ci["pNG"]]) for l in rows4])
    pRG_tab = np.array([float(l.split()[ci["pRG"]]) for l in rows4])
    kind = np.array([e["kind"] for e in table])
    gi = np.nonzero(kind == "grid")[0]
    late = [i for i, e in enumerate(table) if e["kind"] == "obs" and e["l2_stop"] != "-"]
    wng = wrg = winst = 0.0
    for i in late:
        T = et[i] + (pp.parse_utc_seconds(table[i]["l2_stop"]) - pp.parse_utc_seconds(table[i]["utc"]))
        eng = lre.responses(et[gi], pNG_i[gi], T)["W10x5"] + grain[i]
        erg = lre.responses(et[gi], pRG_i[gi], T)["W10x5"]
        wng, wrg = max(wng, abs(pNG_tab[i] / eng - 1)), max(wrg, abs(pRG_tab[i] / erg - 1))
    other = [i for i in range(len(table)) if i not in set(late)]
    for i in other:
        winst = max(winst, abs(pNG_tab[i] / (pNG_i[i] + grain[i]) - 1), abs(pRG_tab[i] / pRG_i[i] - 1))
    report("e2e c: D2 COPS_model.dat late = W10x5 + grain, else inst.", max(wng, wrg, winst) < 2e-6,
           "late NG %.1e RG %.1e (%d rows); other rows %.1e (%d); table %%.6e" % (wng, wrg, len(late), winst, len(other)))
    # start models: D1 vs D2 before the late window bitwise; D2 late rows = W @ D2's grid rows
    z1, z2 = np.load(D1 / "start_check_model.npz"), np.load(D2 / "start_check_model.npz")
    late_key = pp.time_key(pp.LATE_START_UTC)
    k1 = {int(k): j for j, k in enumerate(z1["time_key_ms"])}
    early = [j for j, k in enumerate(z2["time_key_ms"]) if int(k) < late_key and int(k) in k1]
    e1 = [k1[int(z2["time_key_ms"][j])] for j in early]
    eq_early = (np.array_equal(z2["surface_pNG"][early], z1["surface_pNG"][e1]) and
                np.array_equal(z2["surface_pRG"][early], z1["surface_pRG"][e1]))
    lay = pp.LateLayout(table, pp.read_run_points(run_g), np.array([int(l.split()[0]) for l in table_rows(run_g / "COPS_model.dat")]))
    pos2 = {int(r): j for j, r in enumerate(z2["rows"])}
    lp = [pos2[int(r)] for r in lay.late_rows]
    gp = [pos2[int(r)] for r in lay.grid_rows]
    dW = float(np.max(np.abs(z2["surface_pNG"][lp] - lay.W @ z2["surface_pNG"][gp]) / np.abs(z2["surface_pNG"][lp])))
    late_differs = not np.array_equal(z2["surface_pNG"][lp], z1["surface_pNG"][[k1[int(z2["time_key_ms"][j])] for j in lp]])
    report("e2e c: start model D2 vs D1 (early bitwise, late = W @ grid)", eq_early and dW < 1e-14 and late_differs,
           "%d early rows bitwise equal %s; late |D2 - W grid|/D2 %.1e (%d late); late differs from D1 %s; start "
           "costs D1 %.9e D2 %.9e" % (len(early), eq_early, dW, len(lp), late_differs,
                                      cr3["start_check"]["start"]["cost"], cr4["start_check"]["start"]["cost"]))
    note, kinds, res_rows, nfig, rc_pt, pt_out = products(D2)
    ok = (len(note) == 1 and kinds.count("grid") == n_grid and len(res_rows) == 2 * n_obs and nfig == [n_obs, n_obs]
          and rc_pt == 0 and len(table_rows(pt_out)) == n_obs and cr4["late_window"]["fit_rows_equal_all_rows"] and
          cr4["late_window"]["grain_t_end_equals_late_start"] and m4["model_extension"]["variant"] == "D2")
    report("e2e c: D2 products (grid rows flagged, obs-only residuals/figures)", ok,
           "note %d; grid %d late_obs %d; residual rows %d; figures %s; positive_table.py exit %d rows %d" % (
               len(note), kinds.count("grid"), kinds.count("late_obs"), len(res_rows), nfig, rc_pt,
               len(table_rows(pt_out)) if pt_out.is_file() else -1))
    for line in lines_of(t4)[:10]:
        print("   | " + line[:230])
    # d. resume of D2 (segment 2), refusals
    R2 = root / "d_clip_D2_seg02"
    rc, t5 = run_child(root, R2.name, "driver", ["--variant", "D2", "--point-table", points, "--restart-from",
                                                  D2 / "fit_restart.npz", "--segment", "2"] + gff_head +
                       ["--max-nfev", "3", "--"] + fit_args(run_g, R2, bundle, ft))
    ok, detail = rc == 0, "exit %d" % rc
    if ok:
        cr5 = json.loads((R2 / "fit_manifest.json").read_text())["model_extension"]["clip_refit"]
        a5, m5 = restart_arrays(R2 / "fit_restart.npz")
        sc = cr5["start_check"]
        ok = (sc["resume_reproduction"]["ok"] and len(m5["segments"]) == 2 and m5["segments"][-1]["index"] == 2 and
              cr5["fit_chain"]["resumed_from"][-1]["path"] == str((D2 / "fit_restart.npz").resolve()) and
              m5["model_extension"]["fit_chain"]["start_restart"]["path"] == str((S / "fit_restart.npz").resolve()))
        detail = "start cost rel %.1e; segments %s; chain start %s" % (
            sc["cost_rel_diff_vs_source_final"], [s_["index"] for s_ in m5["segments"]],
            m5["model_extension"]["fit_chain"].get("start_restart", {}).get("path"))
    else:
        detail += "\n" + t5[-2500:]
    report("e2e d: D2 resume reproduces segment 1 (start check)", ok, detail)
    rc6, t6 = run_child(root, "d_refuse_D1_restart_for_D2", "driver",
                        ["--variant", "D2", "--point-table", points, "--restart-from", D1 / "fit_restart.npz"] + gff_head +
                        ["--max-nfev", "3", "--"] + fit_args(run_g, root / "d_refuse1", bundle, ft))
    fix_bad = ["--fix", "Tmax=180", "--fix", "T75=136.11149", "--fix", "Tmin=133.0", "--fix", "c_ft_c1=1"]
    rc7, t7 = run_child(root, "d_refuse_held", "driver",
                        ["--variant", "D1", "--point-table", points, "--start-restart", S / "fit_restart.npz"] + fix_bad +
                        ["--crossing-weight", "10", "--devices", "0", "--template-start", "0", "--max-nfev", "3", "--"]
                        + fit_args(run_g, root / "d_refuse2", bundle, ft))
    ok = (rc6 != 0 and "must continue the same problem" in t6 and rc7 != 0 and "held values" in t7)
    report("e2e d: refusals (D1 restart for D2; other held values)", ok,
           "exit %d / %d; %s | %s" % (rc6, rc7, [l for l in t6.splitlines() if "clip_refit" in l][-1:],
                                      [l for l in t7.splitlines() if "clip_refit" in l][-1:]))


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "--child":
        return child(json.loads(sys.argv[2]))
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", required=True, help="a new directory for the synthetic bases and fits")
    a = ap.parse_args()
    Path(a.scratch).mkdir(parents=True, exist_ok=False)
    test_arguments(a.scratch)
    test_real_tables()
    test_e2e(a.scratch)
    bad = [r for r in RESULTS if not r[1]]
    print("\n%d tests, %d failed" % (len(RESULTS), len(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

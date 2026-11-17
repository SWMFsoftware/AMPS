#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_tools.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_tools.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 669e27c5473308632289efc9df2a24d6ca5abb8338931492306967a9260497dc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: adapter (option C): export folder in the layout of the nonneg chains (marked: no fine tuning) and the frozen
#   nonneg_maps.py / nonneg_contrib.py run unchanged with state_activity replaced by the localnn one (A re-derived bitwise).
# RUN: jobs/localnn_tools.sbatch (frozen model_ext_frozen/localnn_post_20261114T085037Z).
# 2026-11-15: stage 2: export --ft-dir (the chain's final state = A_final of the fine tuning, checked against A_fit
#   exp(delta) < 1e-12; links to the fine tuning's tables) with a fit_state/ subfolder = the fit's own files and
#   tables, used as --fit-dir by maps / contrib (2026-11-15 fix: with one folder nonneg_contrib compared the fit
#   state with the fine tuning's late_window_model.dat, max rel 2.6e-2, and failed, post job 132045); state_activity
#   with the held term.
# RUN 2026-11-15: post 132045 (failed, before the fix); 132812 (f10 / 0.01, fixed, frozen
#   localnn_s2post_20261115T195616Z); posts 132814-132816.
# ==== end of header ====
"""Option C (local non-negative basis): adapter so that the maps and contribution tools of the softplus / squared
chains (the frozen nonneg_maps.py and nonneg_contrib.py, model_ext_frozen/nonneg_contrib_20261112T210012Z, UNCHANGED)
run on a localnn fit's final state.

Those tools read a fit directory and a fine-tuning directory of a nonneg chain and check the softplus / squared
activity's reproduction (nonneg_maps.state_activity). A localnn fit has no fine tuning and another activity model, so:

  export   --seg SEG --base BASE --d2 D2 --point-table PT --design DESIGN --out DIR
           DIR = one folder used as BOTH --fit-dir and --ft-dir, with the files the tools read:
             fit_result.txt, nonneg_facets.npz, late_window_model.dat, COPS_model.dat, fit_manifest.json  (symlinks to SEG)
             fit_restart.npz     NOT a gpu_fit restart: T, point_ids, use_ng, use_rg, template_values = the 7 multipliers,
                                 kappa (meta kind localnn_export)
             finetune_state.npz  NO fine tuning: A_final = A, held_prev = held_next = delta = 0, seen = the fit's seen,
                                 late = 0, grain A, B of the fit (meta kind localnn_export_no_fine_tuning)
             finetune_result.txt ("complete yes", "parameters 0"), fine_tuning.json (fit directory, point table and the
                                 positive_problem.py of the localnn snapshot with sha256), grain_term.json (t_end_et)
             pressure_comparison.txt  every base row, fit = tuned = the localnn model (10 digits, the format of
                                 nonneg_finetune), kinds obs | late_obs | grid
             README.txt
           with --ft-dir FT (a localnn fine tuning of SEG; stage 2, 2026-11-15): DIR is the fine-tuning folder (A_final;
             late_window_model.dat, pressure_comparison.txt, finetune_state.npz, finetune_result.txt of FT) and
             DIR/fit_state the fit folder (fit_result.txt, nonneg_facets.npz, late_window_model.dat, COPS_model.dat,
             fit_manifest.json, localnn_state.npz of SEG; fit_restart.npz, grain_term.json): the tools verify the fit
             state against the FIT's tables (2026-11-15 fix: with one folder the fit state was compared with the fine
             tuning's late_window_model.dat, max rel 2.6e-2, and nonneg_contrib's verification failed, job 132045)
  maps     --state DIR [nonneg_maps 'maps' arguments: --base --geometry --out --procs]
  contrib  --state DIR [nonneg_contrib arguments: --base --geometry --out --procs]
           run the frozen tool with nonneg_maps.state_activity replaced by localnn_state_activity: A = the fit's A,
           re-derived from localnn_state.npz (c, m), the grid's hat functions and TS, bitwise; kappa, seen, bins from the
           fit; model name "localnn"; the checks the tools require (passed) are the localnn ones.
"""

import argparse
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402

CONTRIB_DEFAULT = Path("/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z")
LABELS = ("[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]")


def cmd_export(args):
    seg = Path(args.seg).resolve()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    stage = out.parent / ("%s.incomplete_%s" % (out.name, lp.utc_now().replace(":", "").replace("-", "")[:15]))
    design = lp.load_design(args.design)
    D = lp.Data(args.base, args.d2, args.point_table, design)
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        m = np.asarray(z["m"], float)
        kappa, gA, gB = float(z["kappa"]), float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        A = np.asarray(z["A"], np.float64)
        seen = np.asarray(z["seen"], bool)
    ftd = Path(args.ft_dir).resolve() if args.ft_dir else None
    if ftd is not None:
        import localnn_gpu_eval as lge
        A, gA, gB = lge.ft_state(ftd, seg)
    stage.mkdir(parents=True)
    names = ["fit_result.txt", "nonneg_facets.npz", "COPS_model.dat", "fit_manifest.json", "localnn_state.npz"]
    if ftd is None:
        names.append("late_window_model.dat")
    for name in names:
        os.symlink(seg / name, stage / name)
    if ftd is not None:                      # the fine tuning's own tables and state (fit = tuned at delta = 0)
        for name in ("late_window_model.dat", "pressure_comparison.txt", "finetune_state.npz", "finetune_result.txt"):
            os.symlink(ftd / name, stage / name)
    tag = {"kind": "localnn_export", "segment": str(seg), "note": "written by localnn_tools.py export for the maps / "
           "contribution tools; not a gpu_fit restart and not a fine tuning"}
    np.savez(stage / "fit_restart.npz", T=np.array(lp.USER_T0), point_ids=D.ipoint, use_ng=D.use_ng, use_rg=D.use_rg,
             template_values=m, kappa=np.array([kappa]), meta_json=json.dumps(tag))
    if ftd is None:
        z0 = np.zeros(len(A))
        np.savez(stage / "finetune_state.npz", A_final=A, held_prev=z0, held_next=z0, delta=z0, seen=seen,
                 late=np.zeros(len(A), bool), grain_A_uPa=np.array(gA), grain_B_uPa_per_100s=np.array(gB),
                 meta_json=json.dumps(dict(tag, kind="localnn_export_no_fine_tuning")))
        (stage / "finetune_result.txt").write_text("# localnn export: NO fine tuning (localnn_tools.py)\ncomplete yes\n"
                                                   "parameters 0\n")
    pp_path = HERE / "positive_problem.py"
    ftj = {"kind": "localnn_export", "fit": {"directory": str(out.resolve()),
                                             "point_table": {"path": D.point_table, "sha256": lp.sha256_file(D.point_table)}},
           "modules": {"positive_problem.py": {"path": str(pp_path), "sha256": lp.sha256_file(pp_path)}},
           "segment": str(seg), "ft_dir": str(ftd) if ftd else None}
    if ftd is not None:
        ftj["fit"]["directory"] = str((out / "fit_state").resolve())
    (stage / "fine_tuning.json").write_text(json.dumps(ftj, indent=1))
    gt = json.dumps({"t_end_et": D.et_end, "last_ng_before_gap_et": D.et_before,
                     "source": "localnn_problem.Data (grain_crossing definition)"})
    (stage / "grain_term.json").write_text(gt)
    if ftd is not None:
        # the fit folder of a fine-tuned state: the fit's own files and tables (2026-11-15)
        fs = stage / "fit_state"
        fs.mkdir()
        for name in ("fit_result.txt", "nonneg_facets.npz", "late_window_model.dat", "COPS_model.dat", "fit_manifest.json",
                     "localnn_state.npz"):
            os.symlink(seg / name, fs / name)
        shutil.copyfile(stage / "fit_restart.npz", fs / "fit_restart.npz")
        (fs / "grain_term.json").write_text(gt)
    c = meta["costs"]
    sNG = None
    lw = lp.read_table_cols(seg / "late_window_model.dat", ("kind",))
    cm = lp.read_table_cols(seg / "COPS_model.dat", ("pNG", "pRG"))
    # pNG / pRG of every row from the state (design route, 10 digits)
    st = np.load(seg / "localnn_state.npz", allow_pickle=False)
    import localnn_gpu_eval as lge
    sN, sR, _ = lge.state_rows_surface(design, seg, ftd)
    rows = D.model_rows(sN, sR, gA, gB)
    late = set(int(v) for v in D.late_rows)
    grid = set(int(v) for v in D.grid_rows)
    L = ["# per row of the base: localnn export (localnn_tools.py) -- NO fine tuning: fit = tuned = the localnn model of "
         "%s; kind obs | late_obs (W10x5 + grain) | grid (model only); ln_* = ln(model / reference) where the reference "
         "is > 0 (nan otherwise); pNG includes the grain" % seg,
         "#iPoint UTC et kind useNG useRG refNG_Pa refRG_Pa pNG_fit_Pa pNG_tuned_Pa pRG_fit_Pa pRG_tuned_Pa ln_NG_fit "
         "ln_NG_tuned ln_RG_fit ln_RG_tuned"]
    for i in range(D.nrow):
        kind = "late_obs" if i in late else ("grid" if i in grid else "obs")
        pn, pr = rows["pNG"][i], rows["pRG"][i]
        ln_n = "%.6f" % np.log(pn / D.ref_ng[i]) if D.ref_ng[i] > 0 else "nan"
        ln_r = "%.6f" % np.log(pr / D.ref_rg[i]) if D.ref_rg[i] > 0 else "nan"
        L.append("%d %s %.6f %s %d %d %.6e %.6e %.9e %.9e %.9e %.9e %s %s %s %s" % (
            D.ipoint[i], D.utc[i], D.et[i], kind, D.use_ng[i], D.use_rg[i], D.ref_ng[i], D.ref_rg[i], pn, pn, pr, pr,
            ln_n, ln_n, ln_r, ln_r))
    if ftd is None:
        (stage / "pressure_comparison.txt").write_text("\n".join(L) + "\n")
    if ftd is not None:
        cm = lp.read_table_cols(ftd / "COPS_model.dat", ("pNG", "pRG"))
    chk = lp.compare(rows["pNG"], {k: v["pNG"] for k, v in cm.items()}, D.ipoint, "%.6e")
    (stage / "README.txt").write_text(
        "localnn export for the maps / contribution tools (localnn_tools.py export, %s)\n"
        "Segment: %s (data cost %.12e). localnn_tools.py maps | contrib --state <this folder> (--ft-dir = this folder,\n"
        "--fit-dir = this folder, or its fit_state/ for a fine-tuned state)\n"
        "(the frozen nonneg_maps.py / nonneg_contrib.py with state_activity replaced). fit_restart.npz, finetune_state.npz,\n"
        "finetune_result.txt, fine_tuning.json are NOT a gpu_fit restart and NOT a fine tuning: they carry the localnn\n"
        "state (multipliers as template_values; A_final = A; held = delta = 0) in the layout the tools read.\n"
        "pressure_comparison.txt vs COPS_model.dat (printed digits): %s\n" % (lp.utc_now(), seg, c["data_cost"],
                                                                             json.dumps(chk)))
    stage.rename(out)
    lp.say("export -> %s (%s)" % (out, json.dumps(chk)))
    return 0


def load_frozen(cdir):
    cdir = Path(cdir)
    if str(cdir) not in sys.path:
        sys.path.insert(0, str(cdir))
    spec = importlib.util.spec_from_file_location("nonneg_maps", str(cdir / "nonneg_maps.py"))
    NM = importlib.util.module_from_spec(spec)
    sys.modules["nonneg_maps"] = NM          # nonneg_contrib's "import nonneg_maps as NM" gets this module
    spec.loader.exec_module(NM)
    spec = importlib.util.spec_from_file_location("nonneg_contrib", str(cdir / "nonneg_contrib.py"))
    nc = importlib.util.module_from_spec(spec)
    sys.modules["nonneg_contrib"] = nc       # its pool workers are pickled by module name
    spec.loader.exec_module(nc)
    if nc.NM is not NM:
        raise SystemExit("nonneg_contrib did not import the patched nonneg_maps")
    return NM, nc


def make_state_activity(NM):
    def localnn_state_activity(ft_dir, fit_dir, el, pool_threads):
        ft_dir, fit_dir = Path(ft_dir), Path(fit_dir)
        ftj = json.load(open(ft_dir / "fine_tuning.json"))
        if ftj.get("kind") != "localnn_export":
            raise SystemExit("%s is not a localnn export" % ft_dir)
        want_fit = (ft_dir / "fit_state") if ftj.get("ft_dir") and (ft_dir / "fit_state").is_dir() else ft_dir
        if fit_dir.resolve() != want_fit.resolve():
            raise SystemExit("localnn: --fit-dir must be %s for the export %s" % (want_fit, ft_dir))
        seg = Path(ftj["segment"])
        with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
            meta = json.loads(str(z["meta_json"]))
            c, m, kappa = np.asarray(z["c"]), np.asarray(z["m"]), float(z["kappa"])
        with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
            A, TS, bins, seen = (np.asarray(z[k]) for k in ("A", "TS", "bin", "seen"))
            held = np.asarray(z["held"], float) if "held" in z.files else np.zeros(len(A))
        grid = lp.load_grid(meta["grid"]["path"])
        if lp.sha256_file(Path(meta["grid"]["path"]) / "grid.npz") != meta["grid"]["sha256"]:
            raise SystemExit("localnn: the grid changed since the fit")
        u = lp.facet_field(grid, c)
        A_re = m[bins.astype(np.intp)] * np.maximum(u + TS, 0.0)
        hz = held != 0.0
        A_re[hz] = A_re[hz] * np.exp(held[hz])          # the fit's operations (localnn_fit.write_products)
        n = el["nfaces"]
        chk = {"model": "localnn", "kappa": kappa, "multipliers": [float(v) for v in m], "bin_names": list(lp.BIN_NAMES),
               "bin_labels_deg": list(LABELS), "segment": str(seg), "A_reproduced_bitwise": bool(np.array_equal(A_re, A)),
               "held_facets": int(hz.sum()), "facets": int(len(A)), "seen_fit_facets": int(seen.sum()),
               "kappa_matches_fit_result": NM.parse_result_kappa(seg / "fit_result.txt") is not None}
        late = np.zeros(n, bool)
        A_out = A
        if ftj.get("ft_dir"):
            with np.load(Path(ftj["ft_dir"]) / "finetune_state.npz", allow_pickle=False) as z:
                A_fin, delta, late, params = (np.asarray(z[k]) for k in ("A_final", "delta", "late", "seen"))
            tuned = delta != 0.0
            with np.errstate(divide="ignore", invalid="ignore"):
                relt = np.abs(A_fin[tuned] / (A[tuned] * np.exp(delta[tuned])) - 1.0)
            chk.update({"note": "localnn chain state: A_final of the fine tuning %s of this fit" % ftj["ft_dir"],
                        "A_fit_reproduced_bitwise": chk["A_reproduced_bitwise"], "tuned_facets": int(tuned.sum()),
                        "A_final_equals_A_fit_at_untuned": bool(np.array_equal(A_fin[~tuned], A[~tuned])),
                        "A_final_equals_recomputed_at_untuned": bool(np.array_equal(A_fin[~tuned], A_re[~tuned])),
                        "A_final_vs_recomputed_max_rel_at_tuned": float(relt.max()) if relt.size else 0.0,
                        "ft_parameter_facets": int(params.sum()),
                        "ft_parameters_subset_of_fit_seen": bool(np.all(seen[params.astype(bool)]))})
            A_out = A_fin
            ok_ft = (chk["A_final_equals_A_fit_at_untuned"] and chk["A_final_vs_recomputed_max_rel_at_tuned"] < 1e-12)
        else:
            # the keys the tools' README writers read, with their values for a state WITHOUT fine tuning
            chk.update({"note": "localnn state, no fine tuning: A_final = the fit's A (re-derived from c, m, the grid, "
                                "TS and the held term)",
                        "A_fit_reproduced_bitwise": chk["A_reproduced_bitwise"], "tuned_facets": 0,
                        "A_final_equals_A_fit_at_untuned": True,
                        "A_final_equals_recomputed_at_untuned": chk["A_reproduced_bitwise"],
                        "A_final_vs_recomputed_max_rel_at_tuned": 0.0, "ft_parameter_facets": 0,
                        "ft_parameters_subset_of_fit_seen": True})
            ok_ft = True
        chk["A_min"] = float(A_out.min())
        chk["passed"] = bool(chk["A_reproduced_bitwise"] and len(A) == n and chk["A_min"] >= 0.0 and ok_ft)
        man = json.load(open(seg / "fit_manifest.json"))
        return {"A": A_out, "kappa": kappa, "seen": seen.astype(bool), "bins": bins.astype(np.uint8),
                "late": late.astype(bool), "model": "localnn", "checks": chk, "manifest": man, "ft_json": ftj}
    return localnn_state_activity


def cmd_tool(args, rest):
    NM, nc = load_frozen(args.contrib_code)
    NM.state_activity = make_state_activity(NM)
    fit_dir = Path(args.state) / "fit_state" if (Path(args.state) / "fit_state").is_dir() else Path(args.state)
    argv = list(rest) + ["--ft-dir", args.state, "--fit-dir", str(fit_dir)]
    if args.cmd == "maps":
        return NM.main(["maps"] + argv)
    return nc.main(argv)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("export")
    for k in ("seg", "base", "d2", "point-table", "design", "out"):
        e.add_argument("--" + k, required=True)
    e.add_argument("--ft-dir", help="a localnn fine tuning of SEG: export its final state (A_final, its tables)")
    for name in ("maps", "contrib"):
        p = sub.add_parser(name)
        p.add_argument("--state", required=True)
        p.add_argument("--contrib-code", default=str(CONTRIB_DEFAULT))
    args, rest = ap.parse_known_args(argv)
    if args.cmd == "export":
        if rest:
            raise SystemExit("unknown arguments %s" % rest)
        return cmd_export(args)
    return cmd_tool(args, rest)


if __name__ == "__main__":
    sys.exit(main())

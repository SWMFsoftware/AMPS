#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/d2_tools.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/d2_tools.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 056b4b4ba46b3db5ae5e229f29c1e41c3682926457870a2e085dc1827489080b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the clipped refit D2: export adapter (nonneg layout, A = max(0, activity_linear) re-derived bitwise), maps /
#   contrib with state_activity replaced, verify (forward model vs D2's model and data cost) (2026-11-17)
# RUN: jobs/d2_tools.sbatch (frozen model_ext_frozen/d2_tools_20261117T225752Z)
# ==== end of header ====
"""The clipped refit D2 (runs/shap8_012m_lategrid_dz1.5/clip_refit/D2/final): adapter so that the maps and contribution
tools of the nonneg chains (the frozen nonneg_maps.py / nonneg_contrib.py, UNCHANGED) and the figure scripts run on D2's
final state (2026-11-17, lead / user: full figure sets of the six cases of 12M_stage2_reg_best, incl. D2; pattern of
localnn_tools.py of option C).

D2 is a fit_cops fit with the clipped linear activity ("clipped linear: max(0, activity_linear) (the fitter's semantics)",
its fit_manifest model_extension) and no fine tuning. Its activity, per facet e:
    activity_linear_e = Y_e . a + A_s sum_k T_ek c_k      (a: 25 harmonic coefficients, T: cliff_ge_60deg + ft_c1..ft_c6,
                                                           c: their values, A_s = the enhancement scale; fit_restart.npz)
    A_e = max(0, activity_linear_e)
(re-derived bitwise from the restart: tests/fullset_stage2_feasibility_20261117/d2_probe.log).

  export   --d2 D2 --base BASE --bins-from <a nonneg fit's nonneg_facets.npz on the same base> --out DIR
           DIR = one folder used as BOTH --fit-dir and --ft-dir, with the files the tools read:
             fit_result.txt, COPS_model.dat, fit_manifest.json, late_window_model.dat, grain_term.json (symlinks to D2)
             nonneg_facets.npz  A, P = A, u = activity_linear, TS = 0, bin (slope bins of the nonneg fits on this base,
                                from --bins-from), seen (D2's own activity_facets.dat 'seen'), meta (model "clipped_D2",
                                bin_names, multipliers 1)
             fit_restart.npz    NOT a gpu_fit restart: T, point_ids, use_ng, use_rg (D2's), template_values = 7 ones (no
                                slope multipliers), kappa (D2's) (meta kind d2_export)
             finetune_state.npz NO fine tuning: A_final = A, held = delta = 0, seen, late = 0, grain A, B of D2
             finetune_result.txt ("complete yes", "parameters 0"), fine_tuning.json (kind d2_export; D2's point table and
                                its positive_problem.py with sha256)
             pressure_comparison.txt  every base row: fit = tuned = D2's model per row, 10 significant digits
                                (observation rows: pressure_residuals.dat model_Pa, full precision; COPS_model.dat and the
                                restart's pred_ng / pred_rg carry only its 7 printed digits; kinds from late_window_model.dat)
             README.txt, export_checks.json (A vs D2's activity_facets.dat, seen, bins, pred vs COPS_model.dat)
  maps     --state DIR [nonneg_maps 'maps' arguments: --base --geometry --out --procs]
  contrib  --state DIR [nonneg_contrib arguments: --base --geometry --out --procs]
           the frozen tool with nonneg_maps.state_activity replaced by d2_state_activity (A re-derived bitwise from
           D2's restart and compared with the export's nonneg_facets.npz A and finetune_state.npz A_final).
  verify   --state DIR --contrib DIR2 --out FILE: the forward model of the exported state (the contribution tool's
           per-row sums, slope_bins.npz, + D2's grain term) against D2's COPS_model.dat and data cost.
Also importable: make_state_activity(NM) (contrib_phase_data.py --state-adapter).
"""

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CONTRIB_DEFAULT = Path("/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z")
BIN_NAMES = ["slope_00_40", "slope_40_50", "slope_50_60", "slope_60_70", "slope_70_80", "slope_80_90", "slope_90_180"]
LABELS = ("[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]")
T0 = time.time()


def say(t):
    print("d2_tools [%5.0f s]: %s" % (time.time() - T0, t), flush=True)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def read_Y(base):
    with open(Path(base) / "elements.bin", "rb") as f:
        f.read(12)
        d = np.fromfile(f, dtype="<f8").reshape(-1, 34)
    Y = np.ascontiguousarray(d[:, 7:32])
    del d
    return Y


def d2_activity(d2, base):
    """(A, activity_linear stored, re-derived bitwise?, restart dict)."""
    with np.load(Path(d2) / "fit_restart.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        lin, a, tc, tv = (np.asarray(z[k], np.float64) for k in ("activity_linear", "a", "template_columns", "template_values"))
        r = {k: np.asarray(z[k]) for k in ("T", "point_ids", "use_ng", "use_rg", "kappa", "pred_ng", "pred_rg", "x")}
    r["names"] = [str(n) for n in meta["names"]]
    es = float(meta["enhancement_scale"])
    Y = read_Y(base)
    re_lin = Y @ a + es * (tc @ tv)
    del Y, tc
    return np.maximum(lin, 0.0), lin, bool(np.array_equal(re_lin, lin)), r


def read_activity_facets(d2):
    with open(Path(d2) / "activity_facets.dat") as f:
        for _ in range(3):
            f.readline()
        v = np.fromfile(f, sep=" ")
    return v.reshape(-1, 10)


def residual_models(d2):
    """D2's model per observation row and gauge in full precision: pressure_residuals.dat model_Pa (COPS_model.dat prints 7
    significant digits; this table 16). Returns {"NG": {iPoint: model}, "RG": {...}}."""
    out = {"NG": {}, "RG": {}}
    for line in (Path(d2) / "pressure_residuals.dat").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        out[t[3]][int(t[0])] = float(t[7])
    return out


def cmd_export(args):
    d2, base, out = Path(args.d2).resolve(), Path(args.base).resolve(), Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    stage = out.parent / ("%s.incomplete_%s" % (out.name, utc_now().replace(":", "").replace("-", "")))
    A, lin, rederived, r = d2_activity(d2, base)
    say("A = max(0, activity_linear): %d facets > 0; activity_linear re-derived bitwise: %s" % (int((A > 0).sum()), rederived))
    v = read_activity_facets(d2)
    seen = v[:, 5].astype(bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        relA = np.abs(v[:, 9] - A) / np.abs(A)
    chk = {"activity_linear_rederived_bitwise": rederived,
           "activity_facets_rows": int(v.shape[0]),
           "activity_facets_zero_pattern_equal": bool(np.array_equal(v[:, 9] == 0.0, A == 0.0)),
           "activity_facets_clipped_activity_max_rel_at_A_gt_0": float(np.nanmax(relA[A > 0])),
           "activity_facets_total_linear_max_abs_over_max_abs_linear": float(np.max(np.abs(v[:, 8] - lin)) / np.max(np.abs(lin)))}
    del v
    with np.load(args.bins_from, allow_pickle=False) as z:
        bins = np.asarray(z["bin"]).astype(np.int8)
        seen_b = np.asarray(z["seen"]).astype(bool)
    chk["seen_equal_bins_from_fit_seen"] = bool(np.array_equal(seen, seen_b))
    chk["bins_from"] = {"path": str(Path(args.bins_from).resolve()), "sha256": sha256_file(args.bins_from)}
    names = r["names"]
    ik = names.index("kappa")
    kappa = float(np.asarray(r["kappa"]).reshape(-1)[0])
    chk["kappa_restart_equals_x"] = bool(kappa == float(r["x"][ik]))
    gA = float(r["x"][names.index("grain_A_uPa")])
    gB = float(r["x"][names.index("grain_B_uPa_per_100s")])
    # per-row model of D2 (pred_*): vs COPS_model.dat (7 significant digits)
    cols, rows = None, []
    for line in (d2 / "COPS_model.dat").read_text().splitlines():
        if line.startswith("#iPoint"):
            cols = line[1:].split()
        elif not line.startswith("#") and line.strip():
            rows.append(line.split())
    ci = {c: k for k, c in enumerate(cols)}
    ipt = np.array([int(t[ci["iPoint"]]) for t in rows])
    if not np.array_equal(ipt, np.asarray(r["point_ids"])):
        raise SystemExit("COPS_model.dat rows differ from the restart's point ids")
    pNG = np.array([float(t[ci["pNG"]]) for t in rows])
    pRG = np.array([float(t[ci["pRG"]]) for t in rows])
    rm = residual_models(d2)
    full = {"NG": np.array([rm["NG"].get(int(i), np.nan) for i in ipt]), "RG": np.array([rm["RG"].get(int(i), np.nan) for i in ipt])}
    with np.errstate(divide="ignore", invalid="ignore"):
        chk["pred_vs_COPS_model_max_rel"] = {"NG": float(np.nanmax(np.abs(r["pred_ng"] / pNG - 1.0))),
                                             "RG": float(np.nanmax(np.abs(r["pred_rg"] / pRG - 1.0)))}
        chk["pressure_residuals_model_vs_COPS_model_max_rel"] = {
            "NG": float(np.nanmax(np.abs(full["NG"] / pNG - 1.0))), "RG": float(np.nanmax(np.abs(full["RG"] / pRG - 1.0)))}
    chk["pressure_residuals_rows"] = {"NG": int(np.isfinite(full["NG"]).sum()), "RG": int(np.isfinite(full["RG"]).sum())}
    kinds = {}
    for line in (d2 / "late_window_model.dat").read_text().splitlines():
        if not line.startswith("#") and line.strip():
            t = line.split()
            kinds[int(t[0])] = t[3]
    stage.mkdir(parents=True)
    for name in ("fit_result.txt", "COPS_model.dat", "fit_manifest.json", "late_window_model.dat", "grain_term.json"):
        os.symlink(d2 / name, stage / name)
    meta = {"created_utc": utc_now(), "kind": "d2_export", "model": "clipped_D2", "bin_names": BIN_NAMES,
            "multipliers": [1.0] * 7, "kappa": kappa, "d2": str(d2),
            "definitions": {"A": "max(0, activity_linear) (D2, the fitter's semantics)", "u": "activity_linear of D2's restart",
                            "TS": "0 (D2's templates are inside activity_linear)", "P": "A",
                            "bin": "the slope bins of the nonneg fits on this base (%s)" % args.bins_from,
                            "seen": "D2's activity_facets.dat seen (a facet contributes to a loaded selected-data kernel)"}}
    np.savez_compressed(stage / "nonneg_facets.npz", u=lin, TS=np.zeros(len(A)), A=A, P=A, bin=bins, seen=seen,
                        meta_json=np.array(json.dumps(meta)))
    tag = {"kind": "d2_export", "d2": str(d2), "note": "written by d2_tools.py export for the maps / contribution tools; "
           "not a gpu_fit restart and not a fine tuning"}
    np.savez(stage / "fit_restart.npz", T=np.asarray(r["T"]), point_ids=np.asarray(r["point_ids"]),
             use_ng=np.asarray(r["use_ng"]), use_rg=np.asarray(r["use_rg"]), template_values=np.ones(7),
             kappa=np.array([kappa]), meta_json=json.dumps(tag))
    z0 = np.zeros(len(A))
    np.savez(stage / "finetune_state.npz", A_final=A, held_prev=z0, held_next=z0, delta=z0, seen=seen,
             late=np.zeros(len(A), bool), grain_A_uPa=np.array(gA), grain_B_uPa_per_100s=np.array(gB),
             meta_json=json.dumps(dict(tag, kind="d2_export_no_fine_tuning")))
    (stage / "finetune_result.txt").write_text("# D2 export: NO fine tuning (d2_tools.py)\ncomplete yes\nparameters 0\n")
    man = json.load(open(d2 / "fit_manifest.json"))
    cr = man["model_extension"]["clip_refit"]
    pp = cr["uses"]["positive_problem"]
    ftj = {"kind": "d2_export", "d2": str(d2),
           "fit": {"directory": str(out.resolve()), "point_table": {"path": cr["point_table"]["path"],
                                                                     "sha256": cr["point_table"]["sha256"]}},
           "modules": {"positive_problem.py": {"path": pp["module"], "sha256": pp["sha256"]}}}
    if sha256_file(pp["module"]) != pp["sha256"] or sha256_file(cr["point_table"]["path"]) != cr["point_table"]["sha256"]:
        raise SystemExit("D2's positive_problem.py or point table changed since the fit")
    (stage / "fine_tuning.json").write_text(json.dumps(ftj, indent=1))
    L = ["# per row of the base: D2 export (d2_tools.py) -- NO fine tuning: fit = tuned = D2's model per row "
         "(pressure_residuals.dat model_Pa, 16 digits, at the observation rows; fit_restart.npz pred_ng / pred_rg = the "
         "7-digit COPS_model.dat values at the grid rows); kind obs | late_obs | grid; ln_* = ln(model / reference) where the "
         "reference is > 0 (nan otherwise); pNG includes the grain",
         "#iPoint UTC et kind useNG useRG refNG_Pa refRG_Pa pNG_fit_Pa pNG_tuned_Pa pRG_fit_Pa pRG_tuned_Pa ln_NG_fit "
         "ln_NG_tuned ln_RG_fit ln_RG_tuned"]
    for k, t in enumerate(rows):
        ip = int(ipt[k])
        rn, rr = float(t[ci["refNG"]]), float(t[ci["refRG"]])
        # full precision where pressure_residuals.dat has the row (every observation row), else pred (grid rows: the
        # tools read late_window_model.dat there)
        pn = float(full["NG"][k]) if np.isfinite(full["NG"][k]) else float(r["pred_ng"][k])
        pr = float(full["RG"][k]) if np.isfinite(full["RG"][k]) else float(r["pred_rg"][k])
        ln_n = "%.6f" % np.log(pn / rn) if rn > 0 else "nan"
        ln_r = "%.6f" % np.log(pr / rr) if rr > 0 else "nan"
        L.append("%d %s %s %s %d %d %.6e %.6e %.9e %.9e %.9e %.9e %s %s %s %s" % (
            ip, t[ci["UTC"]], t[ci["et"]], kinds.get(ip, "obs"), int(r["use_ng"][k]), int(r["use_rg"][k]), rn, rr,
            pn, pn, pr, pr, ln_n, ln_n, ln_r, ln_r))
    (stage / "pressure_comparison.txt").write_text("\n".join(L) + "\n")
    json.dump(chk, open(stage / "export_checks.json", "x"), indent=1)
    (stage / "README.txt").write_text(
        "D2 export for the maps / contribution tools (d2_tools.py export, %s)\nD2: %s\n"
        "d2_tools.py maps | contrib --state <this folder> (--ft-dir = --fit-dir = this folder): the frozen nonneg_maps.py /\n"
        "nonneg_contrib.py with state_activity replaced by d2_state_activity. fit_restart.npz, finetune_state.npz,\n"
        "finetune_result.txt, fine_tuning.json are NOT a gpu_fit restart and NOT a fine tuning: they carry D2's state\n"
        "(A = max(0, activity_linear); multipliers 1; held = delta = 0) in the layout the tools read.\nChecks: %s\n" % (
            utc_now(), d2, json.dumps(chk)))
    stage.rename(out)
    say("export -> %s %s" % (out, json.dumps(chk)))
    return 0


def make_state_activity(NM):
    def d2_state_activity(ft_dir, fit_dir, el, pool_threads):
        ft_dir, fit_dir = Path(ft_dir), Path(fit_dir)
        ftj = json.load(open(ft_dir / "fine_tuning.json"))
        if ftj.get("kind") != "d2_export":
            raise SystemExit("%s is not a D2 export" % ft_dir)
        if fit_dir.resolve() != ft_dir.resolve():
            raise SystemExit("d2: --fit-dir must be the export %s" % ft_dir)
        d2 = Path(ftj["d2"])
        mrun = NM.re.search(r"^# run: (\S+)$", (d2 / "fit_result.txt").read_text(), NM.re.M)
        A_re, lin, rederived, r = d2_activity(d2, Path(mrun.group(1)))
        with np.load(ft_dir / "nonneg_facets.npz", allow_pickle=False) as z:
            A, bins, seen = (np.asarray(z[k]) for k in ("A", "bin", "seen"))
        with np.load(ft_dir / "finetune_state.npz", allow_pickle=False) as z:
            A_fin = np.asarray(z["A_final"])
        kappa = float(np.asarray(r["kappa"]).reshape(-1)[0])
        n = el["nfaces"]
        chk = {"model": "clipped_D2", "kappa": kappa, "multipliers": [1.0] * 7, "bin_names": BIN_NAMES,
               "bin_labels_deg": list(LABELS), "d2": str(d2), "activity_linear_rederived_bitwise": rederived,
               "A_reproduced_bitwise": bool(np.array_equal(A_re, A)), "A_final_equals_A": bool(np.array_equal(A_fin, A)),
               "note": "D2 state, no fine tuning: A_final = A = max(0, activity_linear), re-derived from D2's restart",
               "A_fit_reproduced_bitwise": bool(rederived and np.array_equal(A_re, A)), "tuned_facets": 0,
               "held_fit_equals_held_prev": True, "held_next_equals_held_prev_plus_delta": True,
               "A_final_equals_A_fit_at_untuned": bool(np.array_equal(A_fin, A)),
               "A_final_equals_recomputed_at_untuned": bool(np.array_equal(A_fin, A_re)),
               "A_final_vs_recomputed_max_rel_at_tuned": 0.0, "ft_parameter_facets": 0,
               "ft_parameter_facets_result_txt": 0, "ft_parameters_subset_of_fit_seen": True, "delta_only_at_parameters": True,
               "seen_fit_facets": int(seen.sum()), "seen_ok": True,
               "kappa_fit_result_txt": NM.parse_result_kappa(d2 / "fit_result.txt"),
               "kappa_matches_fit_result": NM.parse_result_kappa(d2 / "fit_result.txt") is not None,
               "A_final_min": float(A_fin.min()), "A_final_zero_facets": int(np.sum(A_fin == 0.0))}
        chk["passed"] = bool(chk["A_fit_reproduced_bitwise"] and chk["A_final_equals_A"] and len(A) == n and
                             chk["A_final_min"] >= 0.0)
        man = json.load(open(d2 / "fit_manifest.json"))
        return {"A": A_fin, "kappa": kappa, "seen": seen.astype(bool), "bins": bins.astype(np.uint8),
                "late": np.zeros(n, bool), "model": "clipped_D2", "checks": chk, "manifest": man, "ft_json": ftj}
    return d2_state_activity


def load_frozen(cdir):
    cdir = Path(cdir)
    if str(cdir) not in sys.path:
        sys.path.insert(0, str(cdir))
    spec = importlib.util.spec_from_file_location("nonneg_maps", str(cdir / "nonneg_maps.py"))
    NM = importlib.util.module_from_spec(spec)
    sys.modules["nonneg_maps"] = NM
    spec.loader.exec_module(NM)
    spec = importlib.util.spec_from_file_location("nonneg_contrib", str(cdir / "nonneg_contrib.py"))
    nc = importlib.util.module_from_spec(spec)
    sys.modules["nonneg_contrib"] = nc
    spec.loader.exec_module(nc)
    if nc.NM is not NM:
        raise SystemExit("nonneg_contrib did not import the patched nonneg_maps")
    return NM, nc


def cmd_tool(args, rest):
    NM, nc = load_frozen(args.contrib_code)
    NM.state_activity = make_state_activity(NM)
    argv = list(rest) + ["--ft-dir", args.state, "--fit-dir", args.state]
    if args.cmd == "maps":
        return NM.main(["maps"] + argv)
    return nc.main(argv)


def cmd_verify(args):
    """Forward model of the exported state = the contribution tool's per-row model (slope_bins.npz model_surface_NG_Pa /
    model_surface_RG_Pa: kernel sums with the exported A and kappa; W10x5 at the late observations) + D2's grain term at
    the late NG observations (late_window_model.dat grain_NG_Pa), against D2's per-row model (fit_restart.npz pred_ng /
    pred_rg, full precision; COPS_model.dat prints 7 digits) and D2's data cost (0.5 sum over the used rows of
    (ln(model/obs))^2 / N_gauge, the normalized log residuals of pressure_residuals.dat; crossing cost 0 for D2)."""
    import re
    st, cd = Path(args.state), Path(args.contrib)
    d2 = Path(json.load(open(st / "fine_tuning.json"))["d2"])
    with np.load(cd / "slope_bins.npz", allow_pickle=False) as z:
        ip_c, kind = np.asarray(z["iPoint"]), [str(k) for k in z["kind"]]
        surf = {"NG": np.asarray(z["model_surface_NG_Pa"], np.float64), "RG": np.asarray(z["model_surface_RG_Pa"], np.float64)}
        tab = {"NG": np.asarray(z["table_surface_NG_Pa"], np.float64), "RG": np.asarray(z["table_surface_RG_Pa"], np.float64)}
    rm = residual_models(d2)
    with np.load(d2 / "fit_restart.npz", allow_pickle=False) as z:
        pred = {g: np.array([rm[g].get(int(i), np.nan) for i in np.asarray(z["point_ids"])]) for g in ("NG", "RG")}
        ipt = np.asarray(z["point_ids"])
        use = {"NG": np.asarray(z["use_ng"]) > 0, "RG": np.asarray(z["use_rg"]) > 0}
        ref = {"NG": np.asarray(z["ref_ng"], np.float64), "RG": np.asarray(z["ref_rg"], np.float64)}
    if not np.array_equal(ip_c, ipt):
        raise SystemExit("the contribution rows differ from D2's rows")
    grain = {}
    for line in (d2 / "late_window_model.dat").read_text().splitlines():
        if not line.startswith("#") and line.strip():
            t = line.split()
            grain[int(t[0])] = float(t[7])
    late = np.array([k == "late_obs" for k in kind])
    obs = np.array([k != "grid" for k in kind])
    g_ng = np.array([grain.get(int(i), 0.0) if late[r] else 0.0 for r, i in enumerate(ipt)])
    model = {"NG": surf["NG"] + g_ng, "RG": surf["RG"]}
    out = {"rows": int(len(ipt)), "observation_rows": int(obs.sum()), "late_obs_rows": int(late.sum())}
    for g in ("NG", "RG"):
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.abs(model[g] / pred[g] - 1.0)
            dt = np.abs(surf[g] / tab[g] - 1.0)
        out["%s_model_vs_D2_pred_max_rel_observation_rows" % g] = float(np.nanmax(d[obs]))
        out["%s_model_vs_D2_pred_bitwise_rows" % g] = int(np.sum(model[g][obs] == pred[g][obs]))
        out["%s_surface_vs_contrib_table_max_rel" % g] = float(np.nanmax(dt))
    cost, cost_pred = 0.0, 0.0
    for g in ("NG", "RG"):
        u = use[g] & (ref[g] > 0)
        lr = np.log(model[g][u] / ref[g][u])
        lp = np.log(pred[g][u] / ref[g][u])
        cost += 0.5 * float(np.sum(lr ** 2)) / int(u.sum())
        cost_pred += 0.5 * float(np.sum(lp ** 2)) / int(u.sum())
        out["rms_ln_%s" % g] = float(np.sqrt(np.mean(lr ** 2)))
        out["used_rows_%s" % g] = int(u.sum())
    t = (d2 / "fit_result.txt").read_text()
    out.update({"data_cost_forward_model": cost, "data_cost_from_D2_pred": cost_pred,
                "data_cost_fit_result": float(re.search(r"^final_data_cost (\S+)", t, re.M).group(1))})
    out["data_cost_forward_vs_fit_result_rel"] = abs(cost / out["data_cost_fit_result"] - 1.0)
    out["data_cost_pred_vs_fit_result_rel"] = abs(cost_pred / out["data_cost_fit_result"] - 1.0)
    json.dump(out, open(args.out, "x"), indent=1)
    say("verify: %s" % json.dumps(out))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("export")
    for k in ("d2", "base", "bins-from", "out"):
        e.add_argument("--" + k, required=True)
    for name in ("maps", "contrib"):
        p = sub.add_parser(name)
        p.add_argument("--state", required=True)
        p.add_argument("--contrib-code", default=str(CONTRIB_DEFAULT))
    v = sub.add_parser("verify")
    for k in ("state", "contrib", "out"):
        v.add_argument("--" + k, required=True)
    args, rest = ap.parse_known_args(argv)
    if args.cmd == "export":
        if rest:
            raise SystemExit("unknown arguments %s" % rest)
        return cmd_export(args)
    if args.cmd == "verify":
        return cmd_verify(args)
    return cmd_tool(args, rest)


if __name__ == "__main__":
    sys.exit(main())

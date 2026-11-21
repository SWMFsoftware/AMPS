#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_fit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_fit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e553bcd74384cf83fcbcb372a08edfb65de65945058a9ea2e470825098531092
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver (option C): verify (D2 and c9 reproduced), start (NNLS projection of D2 or --start state: another
#   fit's end point; lambda rule), fit (scipy TRF, segments with wall-time stop / exit 75, products fit_result.txt,
#   localnn_state.npz, nonneg_facets.npz, COPS_model.dat, late_window_model.dat, fit_manifest.json, README.txt,
#   HOW-TO-RERUN.md); --smooth-rule-file (the rule that gave --smooth-lambda); --mult-penalty-f (P_m, user 2026-11-14
#   ~17:30 CST; recorded separately); jacobian_check's node step scaled by the state's median x_c.
# RUN: jobs/localnn_build.sbatch, jobs/localnn_fit.sbatch (anchored stage 1: localnn_anch_20261114T090437Z; P_m check: model_ext_frozen/localnn_mp_20261115T020233Z).
# 2026-11-18: stage 2: start --start state --grain-from --held --late-kernels (--state-drop-templates); fit --held
#   --late-kernels --checkpoint-every and --device gpu / --gpu-svd-mode; jacobian_check (2026-11-18): several FD
#   steps per column, the smallest discrepancy counts (1e-4 max(|x|, s), s = 1 and the state scale, each also noise-
#   floored), the one-sided second-order difference whenever a bound clips; regression fixed: the state-scaled steps
#   of 2026-11-15..17 failed test_localnn check 5 at the D2-projection start (6.0e-3 / 9.4e-3), the analytic
#   Jacobian is unchanged and right.
# RUN 2026-11-18: localnn_s2jc3_20261118T081235Z, localnn_s2law_20261118T081952Z; tests 136504, 136507.
# ==== end of header ====
"""Option C (local non-negative basis, user 2026-11-14): driver. CPU only; the model and the cost are in
localnn_problem.py, the design matrices come from localnn_design.py, the grid from localnn_mesh.py.

Subcommands
  verify  D2's final activity and the bounded softplus c9 final state through the evaluator (direct facet sums of the
          design pass + W10x5 + grain + the fit's cost): data costs vs the reported ones, pNG / pRG of every row vs
          COPS_model.dat / late_window_model.dat / pressure_comparison.txt (printed digits), and the design route vs
          the direct facet sums of the two test activities -> <out>/verification.json (exit 1 if a check fails)
  start   the start (--start d2: NNLS projection of max(0, D2) [- TS with --templates on], m = 1, kappa and grain
          from D2), its cost terms, x_c0^T Lg x_c0 and lambda_C of the proposed rule for f = 0.1, 1, 10
          -> <out>/start.npz, start_report.json
  fit     one fit: --start-file <start dir>, --smooth-lambda X [--smooth-rule-file R.json: the rule that gave X,
          checked and recorded; 2026-11-14 the solution-anchored rule] or --smooth-factor F (lambda = F * data cost at the
          start / x_c0^T Lg x_c0), --templates on|off (must equal the start's), --slope-mult-min 1, --max-nfev N,
          --wall-time S -> <out>/segNN (segNN.incomplete_<UTC> renamed when written), <out>/final -> the segment
          that ended by an optimizer condition; a wall-time stop writes the best state so far and exits 75, the
          same command continues from it (new TRF run from that x; not bitwise equal to an uninterrupted run)

Products of a fit segment: fit_result.txt (costs, termination, nfev, parameters, active bounds), localnn_state.npz
(x, names, bounds, c, m, kappa, grain, start), nonneg_facets.npz (u = sum c phi, TS, A, P = A / m_bin, bin, seen,
held = 0, area; the layout of nonneg_problem.py plus area), COPS_model.dat and late_window_model.dat (formats of D2;
pNG / pRG of every row from the design route), fit_manifest.json, README.txt, HOW-TO-RERUN.md.
"""

import argparse
import datetime as dt
import json
import os
import re
import shlex
import socket
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp       # noqa: E402

EXIT_WALLTIME = 75


def say(text):
    lp.say(text)


def stage_dir(out):
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return out.parent / ("%s.incomplete_%s" % (out.name, stamp))


def code_record():
    return {p.name: lp.sha256_file(p) for p in sorted(HERE.glob("*.py"))
            if p.name.startswith(("localnn_", "positive_problem", "test_localnn"))}


def common(args):
    design = lp.load_design(args.design)
    data = lp.Data(args.base, args.d2, args.point_table, design)
    return design, data


def d2_value(data, name):
    return float(data.d2_x[data.d2_names.index(name)])


# ------------------------------------------------------------------------------------------------ verify
def cmd_verify(args):
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    res = verification(args)
    stage = stage_dir(out)
    stage.mkdir(parents=True)
    with open(stage / "verification.json", "x") as f:
        json.dump(res, f, indent=1, sort_keys=True, default=float)
    stage.rename(out)
    say("verification %s -> %s" % ("PASSED" if res["passed"] else "FAILED", out))
    return 0 if res["passed"] else 1


def verification(args):
    t0 = time.time()
    design, D = common(args)
    res = {"created_utc": lp.utc_now(), "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "command": sys.argv, "data": D.summary, "design_meta": {k: design["meta"][k] for k in
                                                                   ("paths", "inputs_sha256", "code", "checks")}}
    ok = True
    d2 = Path(args.d2)
    with np.load(d2 / "fit_restart.npz", allow_pickle=False) as z:
        m2 = json.loads(str(z["meta_json"]))
        kappa2 = float(np.asarray(z["kappa"]).reshape(-1)[0])
    rep2 = float(m2["optimizer"]["cost"])
    gt = json.load(open(d2 / "grain_term.json"))
    states = {"D2": {"sum": design["sum_D2"], "kappa": kappa2, "A": lp_d2_grain(D, "grain_A_uPa"),
                     "B": lp_d2_grain(D, "grain_B_uPa_per_100s"), "reported_data_cost": rep2,
                     "reported_crossing_cost": float(gt["final"]["penalty_cost"]),
                     "t_end_reported": float(gt["t_end_et"])}}
    c9ft, c9fit = Path(args.c9_ft), Path(args.c9_fit)
    with np.load(c9ft / "finetune_state.npz", allow_pickle=False) as z:
        mft = json.loads(str(z["meta_json"]))
        A9, B9 = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    with np.load(c9fit / "fit_restart.npz", allow_pickle=False) as z:
        kappa9 = float(np.asarray(z["kappa"]).reshape(-1)[0])
    states["c9"] = {"sum": design["sum_c9"], "kappa": kappa9, "A": A9, "B": B9,
                    "reported_data_cost": float(mft["costs_final"]["cost_data"]),
                    "reported_crossing_cost": float(mft["costs_final"]["cost_crossing"])}
    for name, s in states.items():
        sNG = s["sum"][:, 0] + s["kappa"] * design["C_NG"]
        sRG = s["sum"][:, 1] + s["kappa"] * design["C_RG"]
        rows = D.model_rows(sNG, sRG, s["A"], s["B"])
        cst = D.costs(rows["pNG"], rows["pRG"])
        v = {"data_cost": cst["data_cost"], "reported_data_cost": s["reported_data_cost"],
             "data_cost_rel_diff": abs(cst["data_cost"] / s["reported_data_cost"] - 1.0),
             "crossing_cost": cst["crossing_cost"], "reported_crossing_cost": s["reported_crossing_cost"],
             "rms_ln_NG": cst["rms_ln_NG"], "rms_ln_RG": cst["rms_ln_RG"], "kappa": s["kappa"],
             "grain_A_uPa": s["A"], "grain_B_uPa_per_100s": s["B"]}
        if "t_end_reported" in s:
            v["t_end_equal"] = bool(D.et_end == s["t_end_reported"])
        if name == "D2":
            cm = lp.read_table_cols(d2 / "COPS_model.dat", ("pNG", "pRG"))
            lw = lp.read_table_cols(d2 / "late_window_model.dat", ("surface_pNG_Pa", "surface_pRG_Pa", "grain_NG_Pa"))
            v["COPS_model_pNG"] = lp.compare(rows["pNG"], {k: x["pNG"] for k, x in cm.items()}, D.ipoint, "%.6e")
            v["COPS_model_pRG"] = lp.compare(rows["pRG"], {k: x["pRG"] for k, x in cm.items()}, D.ipoint, "%.6e")
        else:
            pc = lp.read_table_cols(c9ft / "pressure_comparison.txt", ("pNG_tuned_Pa", "pRG_tuned_Pa"))
            lw = lp.read_table_cols(c9ft / "late_window_model.dat", ("surface_pNG_Pa", "surface_pRG_Pa", "grain_NG_Pa"))
            cm = lp.read_table_cols(c9ft / "COPS_model.dat", ("pNG", "pRG"))
            v["pressure_comparison_pNG"] = lp.compare(rows["pNG"], {k: x["pNG_tuned_Pa"] for k, x in pc.items()},
                                                      D.ipoint, "%.9e")
            v["pressure_comparison_pRG"] = lp.compare(rows["pRG"], {k: x["pRG_tuned_Pa"] for k, x in pc.items()},
                                                      D.ipoint, "%.9e")
            v["COPS_model_pNG"] = lp.compare(rows["pNG"], {k: x["pNG"] for k, x in cm.items()}, D.ipoint, "%.6e")
            v["COPS_model_pRG"] = lp.compare(rows["pRG"], {k: x["pRG"] for k, x in cm.items()}, D.ipoint, "%.6e")
        v["late_window_surface_pNG"] = lp.compare(rows["surface_NG"], {k: x["surface_pNG_Pa"] for k, x in lw.items()},
                                                  D.ipoint, "%.9e")
        v["late_window_surface_pRG"] = lp.compare(rows["surface_RG"], {k: x["surface_pRG_Pa"] for k, x in lw.items()},
                                                  D.ipoint, "%.9e")
        v["late_window_grain"] = lp.compare(rows["grain"], {k: x["grain_NG_Pa"] for k, x in lw.items()}, D.ipoint, "%.9e")
        passed = (v["data_cost_rel_diff"] < 1e-9 and abs(v["crossing_cost"] - v["reported_crossing_cost"]) < 1e-15
                  and all(v[k]["max_rel"] <= 5.0e-7 * (1 + 1e-6) for k in v if k.startswith("COPS_model"))
                  and all(v[k]["max_rel"] <= 5.0e-10 * (1 + 1e-6) + 1e-14 for k in v
                          if k.startswith(("late_window", "pressure_comparison"))))
        v["passed"] = bool(passed)
        ok &= passed
        res[name] = v
        say("%s: data cost %.13e (reported %.13e, rel %.2e), crossing %.3e; tables: %s" % (
            name, v["data_cost"], v["reported_data_cost"], v["data_cost_rel_diff"], v["crossing_cost"],
            json.dumps({k: (v[k]["max_rel"], v[k]["printed_differences"]) for k in v if isinstance(v[k], dict)})))
    tc = design["meta"]["design_vs_direct_test"]
    res["design_vs_direct_test_activities"] = tc
    res["design_vs_direct_passed"] = bool(all(val < 1e-12 for val in tc.values()))
    ok &= res["design_vs_direct_passed"]
    res["passed"] = bool(ok)
    res["code_sha256"] = code_record()
    res["wall_seconds"] = time.time() - t0
    return res


def lp_d2_grain(D, name):
    return d2_value(D, name)


# ------------------------------------------------------------------------------------------------ start
def held_setup(args, design):
    """(held dict for localnn_problem.Problem, record) from --held FINETUNE_STATE (held_next) and --late-kernels, or
    (None, None)."""
    if not getattr(args, "held", None):
        return None, None
    if not args.late_kernels:
        raise SystemExit("--held needs --late-kernels")
    import localnn_late as llate
    with np.load(args.held, allow_pickle=False) as z:
        h = np.asarray(z["held_next"], float)
    with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
        TS, bins = np.asarray(z["TS"], np.float64), np.asarray(z["bin"])
    LK = llate.load(args.late_kernels, len(h))
    import hashlib
    rec = {"path": str(Path(args.held).resolve()), "key": "held_next", "file_sha256": lp.sha256_file(args.held),
           "sha256": hashlib.sha256(np.ascontiguousarray(h).tobytes()).hexdigest(),
           "nonzero_facets": int(np.sum(h != 0.0)), "min": float(h.min()), "max": float(h.max()),
           "late_kernels": {"path": str(Path(args.late_kernels).resolve()), "sha256": lp.sha256_file(args.late_kernels)},
           "definition": "A_e = m_b (sum_j c_j phi_j + TS_e) exp(h_e); h = the previous fine tuning's held_next"}
    return {"h": h, "LK": LK, "bins": bins, "TS": TS}, rec


def make_problem(args, design, D, templates, lam=0.0):
    grid = lp.load_grid(args.grid)
    fm = float(getattr(args, "mult_penalty_f", 0.0) or 0.0)
    held, hrec = held_setup(args, design)
    kw = {}
    if fm > 0:
        kw["mult_penalty_f"] = fm
    if held is not None:
        kw["held"] = held
    P = lp.Problem(D, design, grid, templates, args.slope_mult_min, design["meta"]["A_scale"], lam, **kw)
    P.held_record = hrec
    return P, grid


def cmd_start(args):
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    templates = args.templates == "on"
    design, D = common(args)
    P, grid = make_problem(args, design, D, templates)
    import localnn_mesh as lm
    el = lm.read_elements(Path(args.base) / "elements.bin")
    area = el["area"]
    with np.load(Path(args.d2) / "fit_restart.npz", allow_pickle=False) as z:
        A_D2 = np.maximum(np.asarray(z["activity_linear"], np.float64), 0.0)
        kappa2 = float(np.asarray(z["kappa"]).reshape(-1)[0])
    with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
        TS = np.asarray(z["TS"], np.float64)
    if args.start == "state":
        # another localnn fit's end point (user 2026-11-14 ~17:30 CST: continuation of the f = 10 fit with P_m on)
        sf = Path(args.state_file).resolve()
        with np.load(sf, allow_pickle=False) as z:
            x0 = np.asarray(z["x"], float)
            smeta = json.loads(str(z["meta_json"]))
        dropped = bool(getattr(args, "state_drop_templates", False))
        if dropped and not (smeta["templates"] == "on" and args.templates == "off"):
            raise SystemExit("--state-drop-templates needs a templates-on state and --templates off")
        if (smeta["templates"] != args.templates and not dropped) or list(np.load(sf)["names"]) != list(P.names):
            raise SystemExit("--state-file %s: templates or parameter names differ" % sf)
        if smeta["grid"]["sha256"] != lp.sha256_file(Path(args.grid) / "grid.npz"):
            raise SystemExit("--state-file %s was fitted on another grid" % sf)
        if args.grain_from:
            with np.load(args.grain_from, allow_pickle=False) as z:
                x0[P.iA], x0[P.iB] = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        c0 = x0[P.ic] * P.A_s
        prep = {"start": "state", "state_file": str(sf), "state_sha256": lp.sha256_file(sf),
                "grain_from": str(Path(args.grain_from).resolve()) if args.grain_from else None,
                "held": P.held_record,
                "state_segment": smeta.get("segment"), "state_costs": smeta.get("costs"),
                "state_lambda_C": smeta.get("lambda_C"), "state_run": str(sf.parent),
                "state_templates_dropped": dropped,
                "note": ("the state's x (c, m, kappa, grain) unchanged, D2's templates removed from the model "
                         "(user / lead 2026-11-16: templates OFF test)") if dropped else None}
        say("start from the state %s: %s" % (sf, json.dumps(prep["state_costs"])))
    else:
        target = A_D2 - TS if templates else A_D2
        c0, prep = lp.nnls_projection(grid, area, target)
        say("NNLS projection: " + json.dumps(prep))
        x0 = np.concatenate([c0 / P.A_s, np.ones(lp.NBINS - 1), [kappa2, d2_value(D, "grain_A_uPa"),
                                                                 d2_value(D, "grain_B_uPa_per_100s")]])
    cst, rows = P.evaluate(x0)
    xc0 = x0[P.ic]
    e0 = float(xc0 @ (P.Lg @ xc0))
    lam = {k: lp.lambda_from_factor(f, cst["data_cost"], xc0, P.Lg)[0] for k, f in lp.LAMBDA_FACTORS.items()}
    # activity statistics of the start
    u = lp.facet_field(grid, c0)
    A0 = u + TS if templates else u
    with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
        seen = np.asarray(z["seen"], bool)
    sarea = float(np.sum(area[seen]))
    stats = {"seen_area_fraction_with_A_gt_0": float(np.sum(area[seen & (A0 > 0)]) / sarea),
             "D2_seen_area_fraction_with_A_gt_0": float(np.sum(area[seen & (A_D2 > 0)]) / sarea),
             "total_area_activity_start": float(np.sum(area * A0)), "total_area_activity_D2": float(np.sum(area * A_D2))}
    rep = {"created_utc": lp.utc_now(), "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "command": sys.argv, "templates": args.templates, "slope_mult_min": args.slope_mult_min,
           "projection": prep, "activity": stats,
           "costs_at_start": {k: cst[k] for k in ("data_cost", "crossing_cost", "rms_ln_NG", "rms_ln_RG", "wrong_side",
                                                   "late_rms_ln_NG", "late_rms_ln_RG")},
           "xc0_Lg_xc0": e0, "lambda_rule": "lambda_C = f * data_cost(start) / (x_c0^T Lg x_c0), x_c = c / A_s",
           "lambda_C": lam, "penalty_at_start": {k: lam[k] * e0 for k in lam}, "factors": lp.LAMBDA_FACTORS,
           "smoothness_factor": P.smooth, "data": D.summary, "A_scale": P.A_s, "code_sha256": code_record(),
           "wall_seconds": time.time() - t0}
    say("start: " + json.dumps({k: rep[k] for k in ("costs_at_start", "xc0_Lg_xc0", "lambda_C", "activity")}))
    stage = stage_dir(out)
    stage.mkdir(parents=True)
    np.savez(stage / "start.npz", x0=x0, c0=c0, names=np.array(P.names), templates=np.array(args.templates),
             meta_json=json.dumps(rep, sort_keys=True, default=float))
    with open(stage / "start_report.json", "x") as f:
        json.dump(rep, f, indent=1, sort_keys=True, default=float)
    stage.rename(out)
    return 0


# ------------------------------------------------------------------------------------------------ fit
def cmd_fit(args):
    t0 = time.time()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "final").exists():
        say("%s/final exists: the fit is complete" % out)
        return 0
    templates = args.templates == "on"
    with np.load(Path(args.start_file) / "start.npz", allow_pickle=False) as z:
        x_start = np.asarray(z["x0"], float)
        st_meta = json.loads(str(z["meta_json"]))
        st_templates = str(z["templates"])
    if st_templates != args.templates:
        raise SystemExit("the start was made with --templates %s" % st_templates)
    if (args.smooth_lambda is None) == (args.smooth_factor is None):
        raise SystemExit("give exactly one of --smooth-lambda and --smooth-factor")
    if args.smooth_factor is not None:
        if args.smooth_rule_file:
            raise SystemExit("--smooth-rule-file goes with --smooth-lambda")
        lam = args.smooth_factor * st_meta["costs_at_start"]["data_cost"] / st_meta["xc0_Lg_xc0"]
        lam_src = {"factor": args.smooth_factor, "rule": st_meta["lambda_rule"]}
    else:
        lam, lam_src = float(args.smooth_lambda), {"explicit": float(args.smooth_lambda)}
        if args.smooth_rule_file:
            # a rule decided elsewhere (e.g. the solution-anchored rule, user 2026-11-14): lambda must be one of its values
            rule = json.loads(Path(args.smooth_rule_file).read_text())
            hit = [k for k, v in rule["lambda_C"].items() if float(v) == lam]
            if not hit:
                raise SystemExit("--smooth-lambda %r is not one of the rule's lambda_C %s" % (lam, rule["lambda_C"]))
            lam_src.update({"factor": float(hit[0]), "rule": rule["rule"], "rule_file": str(Path(args.smooth_rule_file).resolve()),
                            "rule_sha256": lp.sha256_file(args.smooth_rule_file), "rule_record": rule})
    design, D = common(args)
    P, grid = make_problem(args, design, D, templates, lam)
    if list(st_meta.get("data", {}).get("late_layout", {}).items()) and st_meta["data"]["t_end_et"] != D.et_end:
        raise SystemExit("the start's data layout differs")
    segs = sorted(int(m.group(1)) for m in (re.match(r"^seg(\d\d)$", p.name) for p in out.iterdir()) if m)
    prev = None
    nfev_used = 0
    x0 = x_start.copy()
    if segs:
        prev = out / ("seg%02d" % segs[-1])
        with np.load(prev / "localnn_state.npz", allow_pickle=False) as z:
            pm = json.loads(str(z["meta_json"]))
            x0 = np.asarray(z["x"], float)
        if pm["stop"] != "wall_time":
            raise SystemExit("%s ended with %s but final is missing" % (prev, pm["stop"]))
        for key in ("lambda_C", "templates", "slope_mult_min", "start_file"):
            if pm[key] != {"lambda_C": lam, "templates": args.templates, "slope_mult_min": args.slope_mult_min,
                           "start_file": str(Path(args.start_file).resolve())}[key]:
                raise SystemExit("resume: %s differs from %s" % (key, prev))
        if float((pm.get("mult_penalty") or {}).get("f_m", 0.0)) != float(args.mult_penalty_f):
            raise SystemExit("resume: --mult-penalty-f differs from %s" % prev)
        if (pm.get("held") or {}).get("sha256") != (P.held_record or {}).get("sha256"):
            raise SystemExit("resume: the held term differs from %s" % prev)
        nfev_used = int(pm["nfev_total"])
    seg = (segs[-1] + 1) if segs else 1
    settings_key = {"lambda_C": lam, "templates": args.templates, "mult_penalty_f": float(args.mult_penalty_f),
                    "start_file": str(Path(args.start_file).resolve()),
                    "held_sha256": (P.held_record or {}).get("sha256"), "max_nfev": args.max_nfev}
    ckpath = out / "checkpoint.npz"
    resumed_ck = None
    if ckpath.exists():
        with np.load(ckpath, allow_pickle=False) as z:
            ck = {"x": np.asarray(z["x"], float), "nfev_total": int(z["nfev_total"]),
                  "settings": json.loads(str(z["settings_json"])), "objective": float(z["objective"]),
                  "segment": int(z["segment"])}
        if ck["settings"] != settings_key:
            raise SystemExit("checkpoint %s was written with other settings: %s" % (ckpath, ck["settings"]))
        if ck["nfev_total"] > nfev_used:
            x0, nfev_used = ck["x"], ck["nfev_total"]
            resumed_ck = {"path": str(ckpath), "sha256": lp.sha256_file(ckpath), "nfev_total": ck["nfev_total"],
                          "objective": ck["objective"], "segment_of_checkpoint": ck["segment"]}
            say("resume from the checkpoint %s (%d evaluations, objective %.10e)" % (ckpath, ck["nfev_total"],
                                                                                    ck["objective"]))
            moved = out / ("checkpoint_resumed_seg%02d_%s.npz" % (seg, dt.datetime.now(dt.timezone.utc).strftime(
                "%Y%m%dT%H%M%SZ")))
            ckpath.rename(moved)
            resumed_ck["kept_as"] = str(moved)
    budget = args.max_nfev - nfev_used
    if budget <= 0:
        raise SystemExit("no evaluations left")
    x0 = np.minimum(np.maximum(x0, P.lo), P.hi)
    c_start, _ = P.evaluate(x_start)
    c_seg0, _ = P.evaluate(x0)
    if P.fm > 0:
        say("multiplier penalty: f_m %r, C_ref %r: P_m at the fit start %.6e, at the segment start %.6e" % (
            P.fm, P.cref_m, c_start["mult_penalty"], c_seg0["mult_penalty"]))
    say("segment %d: lambda_C %.6e (%s), %d parameters, %d residual rows, budget %d evaluations; cost at the segment "
        "start: data %.10e crossing %.3e penalty %.6e objective %.10e" % (
            seg, lam, json.dumps(lam_src), len(x0), D.NNG + D.NRG + len(D.cross["sign"]) + P.R.shape[0], budget,
            c_seg0["data_cost"], c_seg0["crossing_cost"], c_seg0["penalty"], c_seg0["objective"]))
    if args.check_jacobian:
        chk = jacobian_check(P, x0)
        say("Jacobian check: " + json.dumps(chk))
        if not chk["passed"]:
            raise SystemExit("Jacobian check failed")
    else:
        chk = None
    P.deadline = (t0 + args.wall_time) if args.wall_time else None
    dev_rec = None
    if getattr(args, "device", "cpu") == "gpu":
        import localnn_gpu_solver as lgs
        dev_rec = lgs.install(P, args.gpu_svd_mode)
        say("device gpu: " + json.dumps(dev_rec))
    if args.checkpoint_every and args.checkpoint_every > 0:
        def write_ck(xb, ob, nfev_seg):
            tmp = out / "checkpoint.npz.tmp.npz"
            np.savez(tmp, x=xb, objective=np.array(ob), nfev_total=np.array(nfev_used + nfev_seg),
                     segment=np.array(seg), utc=np.array(lp.utc_now()),
                     settings_json=np.array(json.dumps(settings_key, sort_keys=True, default=float)))
            os.replace(tmp, ckpath)
        P.checkpoint = {"every": int(args.checkpoint_every), "write": write_ck}
    stop, res = None, None
    t1 = time.time()
    try:
        res = least_squares(P.fun, x0, jac=P.jac, bounds=(P.lo, P.hi), x_scale="jac", method="trf",
                            tr_solver="exact", ftol=lp.TOL, xtol=lp.TOL, gtol=lp.TOL, max_nfev=budget, verbose=0)
        x = res.x
        stop = "optimizer"
    except lp.WallTime:
        x = P.best[1]
        stop = "wall_time"
    finally:
        if dev_rec is not None:
            dev_rec = lgs.stats()
            lgs.uninstall()
    tfit = time.time() - t1
    cst, rows = P.evaluate(x)
    nfev_total = nfev_used + P.nfev
    say("segment %d ended (%s) after %d evaluations (%.0f s, %.2f s per evaluation): data %.10e crossing %.3e "
        "penalty %.6e objective %.10e" % (seg, stop if stop == "wall_time" else res.message, P.nfev, tfit,
                                          tfit / max(P.nfev, 1), cst["data_cost"], cst["crossing_cost"],
                                          cst["penalty"], cst["objective"]))
    if P.fm > 0:
        say("multiplier penalty at the end: P_m %.6e" % cst["mult_penalty"])
    segdir = out / ("seg%02d" % seg)
    stage = stage_dir(segdir)
    stage.mkdir(parents=True)
    info = {"stop": stop, "segment": seg, "nfev": P.nfev, "nfev_total": nfev_total, "fit_seconds": tfit,
            "seconds_per_evaluation": tfit / max(P.nfev, 1), "lambda_C": lam, "lambda_source": lam_src,
            "templates": args.templates, "slope_mult_min": args.slope_mult_min,
            "start_file": str(Path(args.start_file).resolve()), "max_nfev": args.max_nfev,
            "previous_segment": str(prev) if prev else None, "jacobian_check": chk,
            "resumed_from_checkpoint": resumed_ck, "checkpoint_every": int(args.checkpoint_every or 0),
            "cost_at_fit_start": {k: c_start[k] for k in ("data_cost", "crossing_cost", "penalty", "objective",
                                                           "mult_penalty") if k in c_start},
            "cost_at_segment_start": {k: c_seg0[k] for k in ("data_cost", "crossing_cost", "penalty", "objective",
                                                              "mult_penalty") if k in c_seg0}}
    if P.held_record is not None:
        info["held"] = P.held_record
    if dev_rec is not None:
        info["device"] = dev_rec
        say("device gpu: %d SVDs in %.1f s" % (dev_rec["svd_calls"], dev_rec["svd_seconds"]))
    if P.fm > 0:
        info["mult_penalty"] = {"f_m": P.fm, "C_ref": P.cref_m, "definition": lp.MULT_PENALTY_DEFINITION,
                                "value_final": cst["mult_penalty"]}
    if res is not None:
        info.update(message=res.message, status=int(res.status), njev=int(res.njev) if res.njev is not None else None,
                    optimality=float(res.optimality), lsq_cost=float(res.cost))
    else:
        info.update(message="wall-time stop at the best evaluated state", status=None, njev=None, optimality=None)
    write_products(stage, args, P, D, design, grid, x, cst, rows, info, templates)
    stage.rename(segdir)
    if ckpath.exists():
        ckpath.rename(segdir / "checkpoint_last.npz")
    if stop == "optimizer":
        os.symlink(segdir.name, out / "final")
        say("final -> %s" % segdir.name)
        return 0
    say("wall-time stop: resubmit the same command to continue (exit %d)" % EXIT_WALLTIME)
    return EXIT_WALLTIME


def node_fd_scale(xc):
    """Typical node coefficient for the FD steps of the Jacobian checks: min(1, 75th percentile of the positive x_c)."""
    xc = np.asarray(xc, float)
    return min(1.0, float(np.percentile(xc[xc > 0], 75))) if np.any(xc > 0) else 1.0


def jacobian_check(P, x, ncols=24, seed=1):
    """Analytic vs central-difference columns of the full residual vector (data + crossing + smoothness)."""
    rng = np.random.default_rng(seed)
    J = P.jac(x)
    seen_cols = np.flatnonzero(np.any(J[:, P.ic] != 0, axis=0))
    cols = sorted(set(rng.choice(seen_cols, size=min(ncols, len(seen_cols)), replace=False).tolist()) |
                  set(range(P.im.start, len(x))))
    worst, rep = 0.0, {}
    nf0 = P.nfev
    typ = np.ones(len(x))
    typ[P.ik:] = 1e-3                     # kappa, grain A, B (start values ~0 for kappa and B)
    # node columns: the state's typical x_c (2026-11-14: at a converged state x_c ~ 1e-3, a step of 1e-4 x 1 is ~5 % of
    # the field and the FD truncation error of ln p reaches 5e-4; at the D2-projection starts x_c is O(1 .. 100)).
    # 2026-11-15: the 75th percentile of the positive x_c, not the median: with P_m half of the nodes sit at the bound
    # 0 as 1e-24 .. 1e-19 (TRF interior), the median 8e-20 gave steps far below the rounding of the residuals
    # (FD column 0, rel. error 1.0); the 75th percentile is 1.1e-3 .. 1.5e-3 at the three stage-1 end points
    xc = x[P.ic]
    typ[P.ic] = node_fd_scale(xc)
    # 2026-11-16: a noise floor of the step. At the f = 1, f_m = 1e-4 c4 start (chain stage2_f1_pm1e-4) the 75th
    # percentile of the positive x_c is 3.8e-5, so h = 3.8e-9 at the bound nodes and the FD rounding noise
    # eps max|f| / h reached 2.9e-5 of the column (it fell as 1/h: 3.9e-6, 9.5e-7, 2.9e-7, 3.0e-8 at 10, 30, 100,
    # 1000 h; the analytic columns are right). The floor keeps eps max|f| / h <= 1e-8 of the column's max |J|
    # (capped at 1000 x the base step; 2026-11-17: 100 x left 6.7e-6 at the stage2r f1 / 1e-4 c4 start; at 1000 h the
    # truncation error of those columns was 3e-8).
    # 2026-11-18: SEVERAL steps per column, the smallest discrepancy counts (regression found by the night-fix tests:
    # at the D2-projection start of test_localnn check 5 the state-scaled steps (median / 75th percentile, from
    # localnn_mp_20261115T020233Z on) and their noise floor gave 1.08 / 5.7e-2 / 6.2e-3 where the 2026-11-14 step
    # 1e-4 max(|x|, 1) gives 6.8e-6; at the converged stage-2 states the 1e-4 x 1 step is too large and the scaled /
    # noise-floored steps are needed). Candidates: 1e-4 max(|x|, s) for s = 1 and s = the scale typ (75th percentile of
    # the positive x_c for the nodes; 1 / 1e-3 for the other parameters), each also with the noise floor
    # eps max|f| / (1e-8 max|J_col|) (at most 1000 x). A wrong analytic column disagrees at every step.
    fmax = float(np.max(np.abs(P.fun(x))))
    eps = np.finfo(float).eps
    steps_used = {}
    for k in cols:
        colmax = float(np.max(np.abs(J[:, k])))
        den = max(colmax, 1e-300)
        cand = set()
        for s_ in (1.0, typ[k]):
            h0 = 1e-4 * max(abs(x[k]), s_)
            cand.add(h0)
            if colmax > 0:
                cand.add(max(h0, min(eps * fmax / (1e-8 * colmax), 1000.0 * h0)))
        best_e, best_h = np.inf, None
        for h in sorted(cand):
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] = max(x[k] - h, P.lo[k]) if np.isfinite(P.lo[k]) else x[k] - h
            if xm[k] > x[k] - h:
                # the lower point clipped by the bound: second-order one-sided difference from x (2026-11-15 / 17: the
                # asymmetric two-point difference is first-order; at m = 1 + 5e-5 its error reached 2.2e-5 / 2.5e-5)
                xp2 = x.copy()
                xp2[k] += 2.0 * h
                num = (-3.0 * P.fun(x) + 4.0 * P.fun(xp) - P.fun(xp2)) / (2.0 * h)
            else:
                num = (P.fun(xp) - P.fun(xm)) / (xp[k] - xm[k])
            e = float(np.max(np.abs(num - J[:, k])) / den)
            if e < best_e:
                best_e, best_h = e, h
        rep[P.names[k]] = best_e
        steps_used[P.names[k]] = best_h
        worst = max(worst, best_e)
    P.nfev = nf0
    P.best = None
    return {"columns": len(cols), "max_rel_error": worst, "per_column": rep, "passed": bool(worst < 1e-5),
            "steps": steps_used}


# ------------------------------------------------------------------------------------------------ products
def write_products(stage, args, P, D, design, grid, x, cst, rows, info, templates):
    import localnn_mesh as lm
    c, m, kappa, A, B = P.unpack(x)
    el = lm.read_elements(Path(args.base) / "elements.bin")
    area = el["area"]
    with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
        TS = np.asarray(z["TS"], np.float64) if templates else np.zeros(len(area))
        seen = np.asarray(z["seen"], bool)
        bins = np.asarray(z["bin"]).astype(np.int8)
    u = lp.facet_field(grid, c)
    Pf = np.maximum(u + TS, 0.0)
    held_arr = np.zeros(len(area))
    if P.held is not None:
        held_arr[P.held["idx"]] = np.asarray(np.load(args.held)["held_next"], float)[P.held["idx"]]
    clipped = int(np.sum(u + TS < 0))
    Af = m[bins.astype(np.intp)] * Pf
    if P.held is not None:                   # the same operations as localnn_ft.LocalnnFitContext's reproduction
        hz = held_arr != 0.0
        Af[hz] = Af[hz] * np.exp(held_arr[hz])
        Pf = Pf.copy()
        Pf[hz] = Pf[hz] * np.exp(held_arr[hz])
    sarea = float(np.sum(area[seen]))
    act = {"seen_area_fraction_A_gt_0": float(np.sum(area[seen & (Af > 0)]) / sarea),
           "area_fraction_A_gt_0": float(np.sum(area[Af > 0]) / np.sum(area)),
           "nodes_at_zero": int(np.sum(x[P.ic] <= 1e-8)), "nodes": P.n, "clipped_facets": clipped,
           "total_area_activity": float(np.sum(area * Af)), "seen_area_activity": float(np.sum(area[seen] * Af[seen]))}
    tol = 1e-8
    act_bounds = {"c_at_lower_bound": int(np.sum(x[P.ic] - P.lo[P.ic] <= tol)),
                  "multipliers_at_lower_bound": [lp.BIN_NAMES[1 + k] for k in range(lp.NBINS - 1)
                                                 if x[P.im][k] - P.lo[P.im][k] <= tol * max(1.0, P.lo[P.im][k])],
                  "kappa_at_0": bool(x[P.ik] <= tol), "grain_A_at_0": bool(x[P.iA] <= tol),
                  "grain_B_at_0": bool(x[P.iB] <= tol)}
    if np.any(x < P.lo) or np.any(x > P.hi):
        raise SystemExit("a parameter is outside its bounds")
    meta = {"schema": lp.SCHEMA, "parameterization": lp.PARAMETERIZATION, "created_utc": lp.utc_now(),
            "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
            "env": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                                                    "SLURM_CPUS_PER_TASK", "SLURM_JOB_NODELIST")},
            "grid": {"path": str(Path(args.grid).resolve()), "sha256": lp.sha256_file(Path(args.grid) / "grid.npz"),
                     "nodes": P.n, "meta": grid["meta"]},
            "design": {"path": str(Path(args.design).resolve()), "sha256": lp.sha256_file(Path(args.design) / "design.npz")},
            "d2": str(Path(args.d2).resolve()), "base": str(Path(args.base).resolve()), "data": D.summary,
            "A_scale": P.A_s, "T": list(lp.USER_T0), "smoothness": P.smooth, "costs": {k: cst[k] for k in (
                "data_cost", "crossing_cost", "penalty", "objective", "rms_ln_NG", "rms_ln_RG", "late_rms_ln_NG",
                "late_rms_ln_RG", "wrong_side", "mult_penalty") if k in cst}, "multipliers": list(m), "kappa": kappa,
            "grain_A_uPa": A,
            "grain_B_uPa_per_100s": B, "activity": act, "active_bounds": act_bounds, "code_sha256": code_record()}
    meta.update(info)
    np.savez(stage / "localnn_state.npz", x=x, names=np.array(P.names), lo=P.lo, hi=P.hi, c=c, m=m, kappa=kappa,
             grain_A_uPa=A, grain_B_uPa_per_100s=B, x_start=np.asarray(np.load(Path(args.start_file) / "start.npz")["x0"]),
             meta_json=json.dumps(meta, sort_keys=True, default=float))
    fmeta = {"created_utc": meta["created_utc"], "parameterization": lp.PARAMETERIZATION, "model": "localnn",
             "templates": args.templates, "multipliers": list(m), "kappa": kappa, "A_scale": P.A_s, "T": list(lp.USER_T0),
             "bin_names": list(lp.BIN_NAMES),
             "definitions": {"u": "sum_j c_j phi_j(e) (coarse hat functions, localnn grid)",
                             "TS": "A_s sum_k ft_ck (D2's templates held at 1; 0 with --templates off)",
                             "A": "m[bin] max(0, u + TS)", "P": "A / m[bin] = max(0, u + TS)", "bin": "slope bin 0..6",
                             "seen": "facet contributes to a fit point (the c9 fit's nonneg_facets.npz seen)",
                             "held": "the held fine-tuning term h (A = m max(0, u + TS) exp(h)); 0 without --held",
                             "area": "facet area [m^2]"}, "held": P.held_record}
    with open(stage / "nonneg_facets.npz", "xb") as f:
        np.savez_compressed(f, u=u, TS=TS, A=Af, P=Pf, bin=bins, seen=seen, held=held_arr, area=area,
                            meta_json=json.dumps(fmeta, sort_keys=True, default=float))
    write_tables(stage, args, D, rows, x, meta)
    write_result(stage, args, P, D, x, cst, meta)
    with open(stage / "fit_manifest.json", "x") as f:
        json.dump(meta, f, indent=1, sort_keys=True, default=float)
    write_readme(stage, meta)
    write_howto(stage, args, meta)


def write_tables(stage, args, D, rows, x, meta):
    d2 = Path(args.d2)
    lines = (d2 / "COPS_model.dat").read_text().splitlines()
    k_hdr = next(i for i, line in enumerate(lines) if line.startswith("#iPoint"))
    names = lines[k_hdr][1:].split()
    c = {n: i for i, n in enumerate(names)}
    hdr = ["# LocalnnFit (model_ext/localnn_fit.py, %s, %s): pNG/pRG of EVERY row recomputed on CPU from the design "
           "matrices of the local non-negative basis (localnn_design.py; FP32 records, FP64 arithmetic) for the state of "
           "localnn_state.npz; every other column is D2's (same base, same selection). Not the fitter's GPU output."
           % (lp.PARAMETERIZATION, meta["created_utc"])]
    body = []
    pos = {int(ip): k for k, ip in enumerate(D.ipoint)}
    for line in lines[:k_hdr]:
        if line.startswith("# LiouvilleCOPS FIT PREDICTION"):
            hdr.append(line.replace("(tools/fit_cops.py", "(D2's header, tools/fit_cops.py"))
        else:
            hdr.append(line)
    hdr.append(lines[k_hdr])
    for line in lines[k_hdr + 1:]:
        if line.startswith("#") or not line.strip():
            body.append(line)
            continue
        v = line.split()
        k = pos[int(v[c["iPoint"]])]
        v[c["pNG"]] = "%.6e" % rows["pNG"][k]
        v[c["pRG"]] = "%.6e" % rows["pRG"][k]
        body.append(" ".join(v))
    with open(stage / "COPS_model.dat", "x") as f:
        f.write("\n".join(hdr + body) + "\n")
    lay = D.layout
    stop = {int(r): float(t) for r, t in zip(lay.late_rows, lay.stop_et)}
    rws = sorted([int(r) for r in lay.grid_rows] + [int(r) for r in lay.late_rows], key=lambda i: lay.et[i])
    with open(stage / "late_window_model.dat", "x") as f:
        f.write("# late window of this fit (model_ext/localnn_fit.py, %s): every 1-s grid row (model only) and every "
                "observation from %s UTC on\n" % (meta["created_utc"], lp.pp.LATE_START_UTC))
        f.write("# kind grid: surface_* = instantaneous surface model at et; kind late_obs: surface_* = W10x5 of the "
                "surface model ending at the L2 STOP time (stop_et); grain_NG_Pa = grain(et) (NG only, unaveraged);\n"
                "# table_pNG/pRG = COPS_model.dat (grid: instantaneous surface + grain; late_obs: W10x5 + grain); "
                "re-average only surface_* of grid rows, never table_pNG\n")
        f.write("#iPoint UTC et kind stop_et surface_pNG_Pa surface_pRG_Pa grain_NG_Pa table_pNG_Pa table_pRG_Pa "
                "refNG_Pa refRG_Pa fitNG fitRG\n")
        for i in rws:
            kind = "late_obs" if i in stop else "grid"
            f.write("%d %s %.6f %s %s %.9e %.9e %.9e %.6e %.6e %.6e %.6e %d %d\n" % (
                int(D.ipoint[i]), D.utc[i], D.et[i], kind, "%.6f" % stop[i] if i in stop else "-",
                rows["surface_NG"][i], rows["surface_RG"][i], rows["grain"][i], rows["pNG"][i], rows["pRG"][i],
                D.ref_ng[i], D.ref_rg[i], int(D.useNG[i]), int(D.useRG[i])))


def write_result(stage, args, P, D, x, cst, meta):
    c, m, kappa, A, B = P.unpack(x)
    L = ["# LiouvilleCOPS local non-negative basis fit (option C; model_ext/localnn_fit.py, %s)" % lp.PARAMETERIZATION,
         "# run: %s" % Path(args.base).resolve(),
         "# A_e = m_b(e) (sum_j c_j phi_j(e) + TS_e), c_j >= 0; %d hat functions (grid %s); templates %s; multipliers: "
         "slope_00_40 held at 1, six free >= %g; T held %s; kappa >= 0; grain A, B >= 0; crossing weight %g; W10x5 from %s" % (
             P.n, meta["grid"]["path"], args.templates, args.slope_mult_min, list(lp.USER_T0), lp.CROSSING_WEIGHT,
             lp.pp.LATE_START_UTC),
         "# data: %d NG, %d RG points; misfit: log residuals, equal weight per gauge (fit_cops); costs = 0.5 sum r^2" % (
             D.NNG, D.NRG),
         "# smoothness: lambda_C x_c^T Lg x_c, x_c = c / A_s (A_s %.17g), Lg = P1 stiffness of the coarse mesh" % P.A_s,
         "lambda_C %.9e" % meta["lambda_C"],
         "lambda_source %s" % json.dumps(meta["lambda_source"]),
         "initial_data_cost %.12e" % meta["cost_at_fit_start"]["data_cost"],
         "initial_penalty %.12e" % meta["cost_at_fit_start"]["penalty"],
         "initial_objective %.12e" % meta["cost_at_fit_start"]["objective"],
         "final_cost %.12e" % cst["data_cost"],
         "final_data_cost %.12e" % cst["data_cost"],
         "final_crossing_cost %.12e" % cst["crossing_cost"],
         "final_penalty %.12e" % cst["penalty"],
         "final_objective %.12e" % cst["objective"]] + ([
         "# multiplier penalty: %s; f_m %r" % (lp.MULT_PENALTY_DEFINITION, P.fm),
         "initial_mult_penalty %.12e" % meta["cost_at_fit_start"]["mult_penalty"],
         "final_mult_penalty %.12e" % cst["mult_penalty"]] if P.fm > 0 else []) + [
         "rms_ln_residual_NG %.4f" % cst["rms_ln_NG"], "rms_ln_residual_RG %.4f" % cst["rms_ln_RG"],
         "late_window_rms_ln_NG %.4f" % cst["late_rms_ln_NG"], "late_window_rms_ln_RG %.4f" % cst["late_rms_ln_RG"],
         "crossing_wrong_side_pairs %d" % cst["wrong_side"],
         "optimizer: %s (status %s, stop %s, nfev %d, nfev_total %d, njev %s, first-order optimality %s)" % (
             meta["message"], meta["status"], meta["stop"], meta["nfev"], meta["nfev_total"], meta["njev"],
             "%.3e" % meta["optimality"] if meta["optimality"] is not None else "-"),
         "segment %d  fit_seconds %.1f  seconds_per_evaluation %.3f" % (meta["segment"], meta["fit_seconds"],
                                                                         meta["seconds_per_evaluation"]),
         "active_bounds %s" % json.dumps(meta["active_bounds"]),
         "# activity: %s" % json.dumps(meta["activity"]),
         "# parameter  start  fitted  lower_bound"]
    xs = np.asarray(np.load(Path(args.start_file) / "start.npz")["x0"])
    for k in list(range(P.im.start, len(x))):
        L.append("%-22s %+.9e %+.9e %g" % (P.names[k], xs[k], x[k], P.lo[k]))
    L.append("# c_j (activity units, = A_s x_c): sum %.6e, max %.6e, nodes > 0: %d of %d" % (
        float(np.sum(c)), float(np.max(c)), int(np.sum(x[P.ic] > 1e-8)), P.n))
    L += ["Tmax %.6f" % lp.USER_T0[0], "T75 %.6f" % lp.USER_T0[1], "Tmin %.6f" % lp.USER_T0[2], "kappa %.6e" % kappa]
    with open(stage / "fit_result.txt", "x") as f:
        f.write("\n".join(L) + "\n")


def write_readme(stage, meta):
    L = ["Local non-negative basis fit (option C, user 2026-11-14) -- %s" % lp.PARAMETERIZATION,
         "Written %s by localnn_fit.py on %s (Slurm job %s)." % (meta["created_utc"], meta["host"], meta["slurm_job_id"]),
         "", "MODEL: see localnn_problem.py (docstring). Grid: %s (%d nodes)." % (meta["grid"]["path"], meta["grid"]["nodes"]),
         "Templates %s, slope-multiplier lower bound %g, lambda_C %.6e (%s)." % (
             meta["templates"], meta["slope_mult_min"], meta["lambda_C"], json.dumps(meta["lambda_source"])),
         "Stop: %s (%s), nfev %d (total %d)." % (meta["stop"], meta["message"], meta["nfev"], meta["nfev_total"]),
         "Costs: data %.10e, crossing %.3e, penalty %.6e, objective %.10e." % (
             meta["costs"]["data_cost"], meta["costs"]["crossing_cost"], meta["costs"]["penalty"], meta["costs"]["objective"]),
         "", "FILES",
         "  fit_result.txt     costs, termination, multipliers, kappa, grain, active bounds",
         "  localnn_state.npz  x (optimizer coordinates: c/A_s, m_1..m_6, kappa, grain A [uPa], B [uPa/100 s]), names,",
         "                     bounds, c [activity units], m (7, m_0 = 1), the start x",
         "  nonneg_facets.npz  per facet u, TS, A, P, bin, seen, held (0), area (the nonneg_problem layout + area)",
         "  COPS_model.dat     D2's table with pNG / pRG of every row from the design route (grid rows model-only)",
         "  late_window_model.dat  the late window row by row (D2's format)",
         "  fit_manifest.json  everything above + inputs, code sha256",
         "  HOW-TO-RERUN.md    the exact commands",
         "NOTE: the maps / contribution tools (nonneg_maps.py, nonneg_contrib.py) check the softplus/squared model's",
         "  reproduction and need an adapter for this model; COPS_model.dat here is the CPU design route (verified",
         "  against the GPU tables of D2 and c9 by localnn_fit.py verify), not a GPU forward run."]
    with open(stage / "README.txt", "x") as f:
        f.write("\n".join(L) + "\n")


def write_howto(stage, args, meta):
    env = ["export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda",
           "export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output",
           "export OMP_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s MKL_NUM_THREADS=%s PYTHONUNBUFFERED=1" % (
               os.environ.get("OMP_NUM_THREADS", "1"), os.environ.get("OPENBLAS_NUM_THREADS", "1"),
               os.environ.get("MKL_NUM_THREADS", "1")), "cd $APP"]
    L = ["# HOW TO RERUN this fit segment (%s)" % meta["created_utc"], "",
         "Code: %s (sha256 in fit_manifest.json code_sha256); interpreter $APP/.venv/bin/python." % HERE, "",
         "Environment:", "```bash"] + env + ["```", "", "This segment (exact command):", "```bash",
                                             " ".join(shlex.quote(a) for a in [sys.executable] + sys.argv), "```", "",
         "Inputs: grid %s, design %s, start %s, base %s, D2 %s." % (meta["grid"]["path"], meta["design"]["path"],
                                                                 meta["start_file"], meta["base"], meta["d2"]),
         "From scratch: localnn_mesh.py (grid) -> localnn_design.py (design) -> localnn_fit.py start -> localnn_fit.py fit;",
         "the job scripts jobs/localnn_build.sbatch and jobs/localnn_fit.sbatch run them (see their headers and",
         "srcCG/LiouvilleCOPS_cuda/model_ext/README.md, section 'Option C: local non-negative basis')."]
    with open(stage / "HOW-TO-RERUN.md", "x") as f:
        f.write("\n".join(L) + "\n")


# ------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def base_args(p):
        p.add_argument("--design", required=True)
        p.add_argument("--base", required=True)
        p.add_argument("--d2", required=True, help="D2 fit directory (selection, references, kappa, grain)")
        p.add_argument("--point-table", required=True)
        p.add_argument("--out", required=True)

    pv = sub.add_parser("verify")
    base_args(pv)
    pv.add_argument("--c9-ft", required=True)
    pv.add_argument("--c9-fit", required=True)
    for name in ("start", "fit"):
        p = sub.add_parser(name)
        base_args(p)
        p.add_argument("--grid", required=True)
        p.add_argument("--templates", choices=("on", "off"), required=True)
        p.add_argument("--slope-mult-min", type=float, default=1.0)
        p.add_argument("--held", help="finetune_state.npz whose held_next is held (stage 2; needs --late-kernels)")
        p.add_argument("--late-kernels", help="localnn_late.py output (G = dp/dA of the late facets)")
    ps = sub.choices["start"]
    ps.add_argument("--start", choices=("d2", "state"), default="d2")
    ps.add_argument("--state-file", help="--start state: a localnn_state.npz (its x is the start)")
    ps.add_argument("--grain-from", help="--start state: grain A, B from this finetune_state.npz")
    ps.add_argument("--state-drop-templates", action="store_true",
                    help="--start state with --templates off from a templates-on state: the same x, the template part "
                         "TS removed (2026-11-16, the templates OFF test)")
    pf = sub.choices["fit"]
    pf.add_argument("--start-file", required=True, help="a 'start' output folder")
    pf.add_argument("--smooth-lambda", type=float)
    pf.add_argument("--smooth-factor", type=float)
    pf.add_argument("--mult-penalty-f", type=float, default=0.0,
                    help="f_m of P_m = f_m C_ref sum_b (ln m_b)^2 (user 2026-11-14 ~17:30 CST; 0 = off, the default)")
    pf.add_argument("--smooth-rule-file", help="JSON of the rule that gave --smooth-lambda (keys rule, lambda_C {f: value}, "
                    "...): checked (lambda must be one of its values) and recorded in the products")
    pf.add_argument("--max-nfev", type=int, required=True)
    pf.add_argument("--wall-time", type=float, default=0.0, help="seconds after which the segment stops (0: none)")
    pf.add_argument("--check-jacobian", action="store_true")
    pf.add_argument("--checkpoint-every", type=int, default=0,
                    help="write <out>/checkpoint.npz (the best x so far) every N evaluations; a later run of the same "
                         "command resumes from it (0: off)")
    pf.add_argument("--device", choices=("cpu", "gpu"), default="cpu",
                    help="gpu (lead / user 2026-11-15 'option 3'): the design products and TRF's SVD on one GPU "
                         "(localnn_gpu_solver.py, CuPy FP64); cpu (default) = the previous path, bitwise")
    pf.add_argument("--gpu-svd-mode", choices=("qr", "gesvd"), default="qr",
                    help="with --device gpu: SVD of J_augmented as QR + SVD(R) (qr) or cuSOLVER gesvd directly")
    args = ap.parse_args(argv)
    if getattr(args, "slope_mult_min", 1.0) != 1.0:
        raise SystemExit("--slope-mult-min: the decision is 1 (user 2026-11-12 / the lead 2026-11-14)")
    return {"verify": cmd_verify, "start": cmd_start, "fit": cmd_fit}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())

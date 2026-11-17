#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_localnn_stage2.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_localnn_stage2.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 7885ab3ff3fb03242c80fec49b1edeb9a55d9c1120cbd42302e1dd515a90b10f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tests of option C stage 2 (CPU): 1 held h = 0 bitwise; 2 held term vs direct facet sums; 3 Jacobian vs FD
#   with held + P_m; 4 checkpoint resume after SIGKILL; 5 localnn_ft run --validate-only (validation checks 1-8).
# RUN 2026-11-15: 131980 (test 3 failed: FD step from the median x_c 8e-20, fixed), 132009 (frozen
#   localnn_s2_20261115T185842Z: 5 of 5 passed).
# ==== end of header ====
"""Tests of option C stage 2 (CPU): the held fine-tuning term, the checkpoints and the CPU fine tuning of a localnn fit.

    test_localnn_stage2.py --grid GRID --base BASE --d2 D2 --point-table PT --late LK --contrib-code FROZEN_CONTRIB
        --state SEG (a localnn fit segment) --rule RULE --work DIR [--out report.json]

 1 held h = 0: the Problem with held (all zeros) = the Problem without held, bitwise (fun, jac)
 2 held term vs direct facet sums: h = 0.3 / -0.2 on 40 late facets; the all-rows design route + held correction vs
   the frozen record arithmetic of A = m (u + TS) exp(h) at 4 rows (2 late grid rows, 1 late observation's grid row
   inside the window, 1 early row) (<= 1e-12)
 3 Jacobian vs finite differences with the held term and P_m (localnn_fit.jacobian_check + every held facet's nodes)
 4 checkpoint resume: localnn_fit.py fit --max-nfev 14 --checkpoint-every 4 killed (SIGKILL) after the checkpoint at 8
   evaluations, then the same command: resumes from the checkpoint (nfev_total continues at 8, the start x = the
   checkpoint's x), the segment ends with nfev_total 14; an uninterrupted run with the same budget for reference
 5 localnn_ft.py run --validate-only on SEG: validation checks 1-8 of the frozen nonneg_finetune pass with the CPU engine
"""

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_ft as lft          # noqa: E402  (imports the frozen nonneg modules first)
import localnn_problem as lp      # noqa: E402
import localnn_fit as lf          # noqa: E402
import localnn_design as ld       # noqa: E402
import localnn_late as llate       # noqa: E402

RES = []
PY = sys.executable


def check(name, ok, **info):
    RES.append({"test": name, "passed": bool(ok), **info})
    print("%s %s %s" % ("PASS" if ok else "FAIL", name, json.dumps(info, default=float)), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    for k in ("grid", "base", "d2", "point-table", "late", "contrib-code", "state", "rule", "work"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=False)
    grid = lp.load_grid(args.grid)
    design = lp.load_design(Path(args.grid) / "design")
    D = lp.Data(args.base, args.d2, args.point_table, design)
    seg = Path(args.state).resolve()
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        x = np.asarray(z["x"], float)
        meta = json.loads(str(z["meta_json"]))
    nf = len(grid["facet_node"])
    LK = llate.load(args.late, nf)
    with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
        TS, bins = np.asarray(z["TS"], np.float64), np.asarray(z["bin"])
    lam, A_s = float(meta["lambda_C"]), design["meta"]["A_scale"]
    # 1
    P0 = lp.Problem(D, design, grid, True, 1.0, A_s, lam, mult_penalty_f=0.01)
    Pz = lp.Problem(D, design, grid, True, 1.0, A_s, lam, mult_penalty_f=0.01,
                    held={"h": np.zeros(nf), "LK": LK, "bins": bins, "TS": TS})
    check("1 held h = 0 is bitwise the Problem without held", np.array_equal(P0.fun(x), Pz.fun(x)) and
          np.array_equal(P0.jac(x), Pz.jac(x)), held_nonzero=Pz.held["nonzero"])
    # 2
    rng = np.random.default_rng(11)
    li = LK["late_idx"]
    sel = rng.choice(li, size=min(40, len(li)), replace=False)
    h = np.zeros(nf)
    h[sel] = np.where(rng.uniform(size=len(sel)) > 0.5, 0.3, -0.2)
    Ph = lp.Problem(D, design, grid, True, 1.0, A_s, lam, mult_penalty_f=0.01, held={"h": h, "LK": LK, "bins": bins,
                                                                                  "TS": TS})
    c, m, kappa, gA, gB = Ph.unpack(x)
    cst, rows = Ph.evaluate(x)
    A = m[bins.astype(np.intp)] * np.maximum(lp.facet_field(grid, c) + TS, 0.0)
    hz = h != 0.0
    A[hz] = A[hz] * np.exp(h[hz])
    S = ld.contrib_setup(args.base, args.contrib_code, lp.USER_T0)
    test_rows = [int(D.grid_rows[-1]), int(D.grid_rows[-60]), int(D.grid_rows[150]), 0]
    worst, det = 0.0, {}
    for r in test_rows:
        dn, dr = ld.direct_row(S, args.base, r, int(D.ipoint[r]), A, kappa)
        en, er = abs(rows["inst_NG"][r] / dn - 1), abs(rows["inst_RG"][r] / dr - 1)
        det["row%d" % r] = [en, er]
        worst = max(worst, en, er)
    # the held facets do change those late rows (the test is not empty)
    st0, rows0 = P0.evaluate(x)
    ch = float(np.max(np.abs(rows["inst_NG"][test_rows[:3]] / rows0["inst_NG"][test_rows[:3]] - 1)))
    check("2 held term vs direct facet sums (4 rows)", worst <= 1e-12 and ch > 1e-12, max_rel=worst, rows=det,
          held_change_of_late_rows=ch)
    # 3
    chk = lf.jacobian_check(Ph, x)
    # plus the nodes of the held facets
    nodes = sorted(set(grid["facet_node"][sel].ravel().tolist()))[:20]
    Jh = Ph.jac(x)
    w3 = 0.0
    xc_typ = lf.node_fd_scale(x[Ph.ic])
    for j in nodes:
        hh = 1e-4 * max(abs(x[j]), xc_typ)
        xp, xm = x.copy(), x.copy()
        xp[j] += hh
        xm[j] = max(x[j] - hh, 0.0)
        num = (Ph.fun(xp) - Ph.fun(xm)) / (xp[j] - xm[j])
        w3 = max(w3, float(np.max(np.abs(num - Jh[:, j])) / max(np.max(np.abs(Jh[:, j])), 1e-300)))
    check("3 Jacobian vs finite differences with held + P_m", chk["passed"] and w3 < 1e-5, max_rel=chk["max_rel_error"],
          held_node_columns=len(nodes), held_node_max_rel=w3)
    # 4 checkpoint resume (subprocesses: the real CLI)
    held_file = work / "held_test.npz"
    np.savez(held_file, held_next=h)
    com = ["--design", str(Path(args.grid) / "design"), "--base", args.base, "--d2", args.d2, "--point-table",
           args.point_table, "--grid", args.grid, "--templates", "on"]
    hel = ["--held", str(held_file), "--late-kernels", args.late]
    env = dict(os.environ)
    r0 = subprocess.run([PY, "-W", "ignore", str(HERE / "localnn_fit.py"), "start"] + com + hel +
                        ["--start", "state", "--state-file", str(seg / "localnn_state.npz"), "--out",
                         str(work / "start")], capture_output=True, text=True, env=env)
    fitcmd = [PY, "-W", "ignore", str(HERE / "localnn_fit.py"), "fit"] + com + hel + [
        "--slope-mult-min", "1", "--start-file", str(work / "start"), "--smooth-lambda", repr(lam), "--smooth-rule-file",
        args.rule, "--mult-penalty-f", "0.01", "--max-nfev", "14", "--checkpoint-every", "4"]
    ref = subprocess.run(fitcmd + ["--out", str(work / "fit_ref")], capture_output=True, text=True, env=env)
    out_k = work / "fit_kill"
    pr = subprocess.Popen(fitcmd + ["--out", str(out_k)], stdout=open(work / "kill_run.log", "w"),
                          stderr=subprocess.STDOUT, env=env)
    killed_at = None
    t0 = time.time()
    while pr.poll() is None and time.time() - t0 < 1500:
        ck = out_k / "checkpoint.npz"
        if ck.exists():
            try:
                with np.load(ck) as z:
                    if int(z["nfev_total"]) >= 8:
                        pr.send_signal(signal.SIGKILL)
                        killed_at = int(z["nfev_total"])
                        x_ck = np.asarray(z["x"], float)
                        break
            except Exception:
                pass
        time.sleep(1.0)
    pr.wait()
    res = subprocess.run(fitcmd + ["--out", str(out_k)], capture_output=True, text=True, env=env)
    ok4, info4 = False, {"killed_at_checkpoint": killed_at, "start_rc": r0.returncode, "ref_rc": ref.returncode,
                         "resume_rc": res.returncode}
    try:
        segs = sorted(p.name for p in out_k.iterdir() if p.name.startswith("seg"))
        with np.load(out_k / "final" / "localnn_state.npz", allow_pickle=False) as z:
            mk = json.loads(str(z["meta_json"]))
        with np.load(work / "fit_ref" / "final" / "localnn_state.npz", allow_pickle=False) as z:
            mr = json.loads(str(z["meta_json"]))
        rc_ = mk.get("resumed_from_checkpoint") or {}
        info4.update(segments=segs, nfev_total=mk["nfev_total"], resumed=rc_, ref_nfev_total=mr["nfev_total"],
                     resume_log="resume from the checkpoint" in res.stdout,
                     segment_start_objective=mk["cost_at_segment_start"]["objective"],
                     ref_final_objective=mr["costs"]["objective"], resumed_final_objective=mk["costs"]["objective"],
                     held_nonzero=(mk.get("held") or {}).get("nonzero_facets"))
        ok4 = (killed_at is not None and killed_at >= 8 and res.returncode == 0 and mk["nfev_total"] == 14 and
               rc_.get("nfev_total") == killed_at and info4["resume_log"] and
               abs(mk["cost_at_segment_start"]["objective"] / rc_["objective"] - 1) < 1e-12 and segs == ["seg01"])
    except Exception as exc:
        info4["error"] = "%s: %s" % (exc.__class__.__name__, exc)
        info4["resume_stdout_tail"] = res.stdout[-2000:] + res.stderr[-2000:]
    check("4 checkpoint resume after SIGKILL", ok4, **info4)
    # 5 the CPU fine tuning of a localnn fit: validation checks 1-8
    v = subprocess.run([PY, "-W", "ignore", str(HERE / "localnn_ft.py"), "run", "--fit", str(seg), "--late", args.late,
                        "--design", str(Path(args.grid) / "design"), "--base", args.base, "--d2", args.d2, "--out",
                        str(work / "ft_validate"), "--facets", "late_only", "--window-min", "5", "--validate-only"],
                       capture_output=True, text=True, env=env)
    vj = {}
    try:
        vj = json.load(open(work / "ft_validate" / "validation.json"))
    except Exception as exc:
        vj = {"error": str(exc)}
    check("5 localnn fine tuning: validation checks 1-8 (frozen nonneg_finetune, CPU engine)",
          v.returncode == 0 and vj.get("passed") is True, rc=v.returncode,
          checks={k: (vv.get("passed") if isinstance(vv, dict) else vv) for k, vv in vj.items()},
          start=vj.get("1_start_cost"), tail=(v.stdout[-600:] if v.returncode else ""))
    ok = all(r["passed"] for r in RES)
    print("%d tests, %d failed" % (len(RES), sum(not r["passed"] for r in RES)))
    if args.out:
        with open(args.out, "x") as f:
            json.dump({"results": RES, "passed": ok, "argv": sys.argv}, f, indent=1, default=float)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

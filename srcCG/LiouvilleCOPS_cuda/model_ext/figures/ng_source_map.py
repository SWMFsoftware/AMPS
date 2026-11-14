#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/ng_source_map.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/ng_source_map.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 96223ea54fdded990129f4922f7ffda6d3d59b7bd3398365eb3a7b6edf46acec
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Per-facet sums of the modelled NG surface signal (frozen nonneg_contrib record_contributions): 95 late
#   observations (W10x5) and all 934 NG observations, Pa and per-observation shares.
# RUN: jobs/ng_source_map.sbatch (normal128cores, CPU)
# ==== end of header ====
"""Where does the MODELLED NG surface signal come from? (lead / user 2026-11-13 ~20:05 CST, for the camera path of the
detach-300 m panorama movie: "the other comet's lobe that contributes most to the NG data").

Per facet e, from the fitted forward model term by term (the frozen nonneg_contrib.py record_contributions, as
fov_share.py and the contribution data: the fit's FP32 records, F = q A with the night / shadow rule, the fit's T and kappa,
H2O + CO2; the NG grain term is not a surface signal and is excluded):
  late_Pa[e]   = sum over the 95 late observations k of c_k,NG,e (W10x5 of the 1-s grid rows, as the fit)
  late_sh[e]   = sum_k c_k,NG,e / C_k (C_k = the observation's total surface contribution): sum of per-observation shares
  all_Pa[e], all_sh[e] = the same over ALL NG observations of the descent (the fit's NG data: observation rows with
                 refNG > 0 used by the fit, before the late window, plus the 95 late observations)
Writes <out>.npz (the four arrays, float64, per facet) and <out>_rows.json (the rows used, totals, checks).
Run (CPU job): ng_source_map.py --contrib-code <frozen nonneg_contrib dir> --ft-dir <c9/ft/final> --fit-dir <c9/fit/final>
    --base <base> --late-npz <contributions_c9/source_areas_NG.npz> --out <prefix> [--procs N]
"""

import argparse
import json
import multiprocessing as mp
import re
import sys
import time
from pathlib import Path

import numpy as np

_G = {}


def chunk_task(task):
    rows_w = task                     # list of (row, ip, w_Pa, w_share or None (= 1 / own total))
    NC = _G["NC"]
    n = _G["nfaces"]
    acc_Pa = np.zeros(n)
    acc_sh = np.zeros(n)
    tot = {}
    for row, ip, wPa, wsh in rows_w:
        kp, h, rec = NC.NM.read_kernel_full(Path(_G["base"]) / "kernels" / ("point=%04d.bin" % ip))
        if kp != ip:
            raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
        face, ((ng, rg),), _ = NC.record_contributions(rec, h, (_G["A"],), "fp32", row)
        c = np.bincount(face, weights=ng, minlength=n)
        t = float(c.sum())
        tot[row] = t
        acc_Pa += wPa * c
        acc_sh += (wsh if wsh is not None else (1.0 / t if t > 0 else 0.0)) * c
    return acc_Pa, acc_sh, tot


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("--contrib-code", "--ft-dir", "--fit-dir", "--base", "--late-npz", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--procs", type=int, default=32)
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out)
    for suf in (".npz", "_rows.json"):
        if Path(str(out) + suf).exists():
            raise SystemExit("%s%s exists" % (out, suf))
    sys.path.insert(0, a.contrib_code)
    import nonneg_contrib as NC
    NM = NC.NM
    base, ft_dir, fit_dir = Path(a.base), Path(a.ft_dir), Path(a.fit_dir)
    fc, header, cols, rows = NM.base_table(base)                 # the same setup as fov_share.py / nonneg_contrib
    helio = fc.header_value(header, "HeliocentricModel")
    mrun = re.search(r"^# run: (\S+)$", (fit_dir / "fit_result.txt").read_text(), re.M)
    if mrun is None or Path(mrun.group(1)).resolve() != base.resolve():
        raise SystemExit("the fit did not run on %s" % base)
    el = NM.read_elements(base / "elements.bin")
    S = NM.state_activity(ft_dir, fit_dir, el, a.procs)
    if not S["checks"]["passed"]:
        raise SystemExit("state checks failed")
    with np.load(fit_dir / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
        point_ids = np.asarray(z["point_ids"], np.int64)
        use_ng = np.asarray(z["use_ng"])
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("restart rows differ from the base rows")
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    utc = [r[cols["UTC"]] for r in rows]
    et = np.array([float(r[cols["et"]]) for r in rows])
    refng = np.array([float(r[cols["refNG"]]) for r in rows])
    bm, cNG, cRG, gconst = NC.gauge_constants(fc, header, TRG)
    ftj = S["ft_json"]
    pprec = ftj["modules"]["positive_problem.py"]
    pp = NM.load_module(pprec["path"], "positive_problem_of_state", pprec["sha256"])
    ptab = ftj["fit"]["point_table"]
    layout = pp.LateLayout(pp.read_point_table(ptab["path"]), pp.read_run_points(base), ipoint)
    Wm = layout.W.tocsr()
    from cops_cuda import device_store
    NC._G.update({"base": str(base), "A": S["A"], "A_fit": S["A"], "ACO2": el["ACO2"], "bins": S["bins"].astype(np.intp),
                  "T": T, "kappa": S["kappa"], "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio, "AU": float(fc.AU),
                  "COS75": float(fc.COS75), "KT": float(fc.KT), "fp32_record_fields": device_store.fp32_record_fields,
                  "nfaces": el["nfaces"], "area": el["area"], "seen": S["seen"]})
    late, gr = list(layout.late_rows), list(layout.grid_rows)
    # late observation totals C_k: from the contribution data (source_areas_NG.npz), aligned by iPoint
    L = np.load(a.late_npz, allow_pickle=False)
    ck = {int(ip): float(t) for ip, t in zip(L["iPoint"], L["total_facets_Pa"])}
    C_late = np.array([ck[int(ipoint[r])] for r in late])
    Wd = Wm.toarray()                                   # (95, n_grid)
    w_Pa = Wd.sum(axis=0)
    w_sh = (Wd / C_late[:, None]).sum(axis=0)
    grid = np.zeros(len(rows), bool)
    grid[gr] = True
    is_late = np.zeros(len(rows), bool)
    is_late[late] = True
    pre = [r for r in range(len(rows)) if refng[r] > 0 and use_ng[r] and not grid[r] and not is_late[r]]
    _G.update({"NC": NC, "nfaces": int(el["nfaces"]), "base": str(base), "A": S["A"]})
    tasks_late = [(gr[j], int(ipoint[gr[j]]), float(w_Pa[j]), float(w_sh[j])) for j in range(len(gr)) if w_Pa[j] > 0]
    tasks_pre = [(r, int(ipoint[r]), 1.0, None) for r in pre]
    print("ng_source_map: %d late grid rows (W10x5 -> 95 observations), %d NG observation rows before the late window" % (
        len(tasks_late), len(tasks_pre)), flush=True)
    res = {}
    for name, tasks in (("late", tasks_late), ("pre", tasks_pre)):
        chunks = [tasks[i::a.procs] for i in range(a.procs) if tasks[i::a.procs]]
        accP, accS, tot = np.zeros(_G["nfaces"]), np.zeros(_G["nfaces"]), {}
        with mp.get_context("fork").Pool(min(a.procs, len(chunks))) as pool:
            for P_, S_, t_ in pool.imap_unordered(chunk_task, chunks):
                accP += P_
                accS += S_
                tot.update(t_)
        res[name] = (accP, accS, tot)
        print("  %s done (%.0f s)" % (name, time.time() - t0), flush=True)
    late_Pa, late_sh, tot_late = res["late"]
    pre_Pa, pre_sh, tot_pre = res["pre"]
    # checks: W10x5 of the grid-row totals vs the contribution data's late totals; sum of shares = 95
    tg = np.array([tot_late.get(r, 0.0) for r in gr])
    C_chk = Wd @ tg
    chk = {"late_totals_vs_contribution_data_max_rel": float(np.max(np.abs(C_chk / C_late - 1))),
           "late_share_sum": float(late_sh.sum()), "pre_share_sum": float(pre_sh.sum()),
           "late_Pa_sum_vs_sum_C_k_rel": float(late_Pa.sum() / C_late.sum() - 1)}
    print("checks:", json.dumps(chk), flush=True)
    np.savez(str(out) + ".npz", late_Pa=late_Pa, late_share=late_sh, all_Pa=late_Pa + pre_Pa, all_share=late_sh + pre_sh,
             pre_Pa=pre_Pa, pre_share=pre_sh)
    json.dump({"definition": __doc__, "inputs": vars(a), "late_observations": len(late), "late_grid_rows": len(tasks_late),
               "pre_rows": len(pre), "pre_first_last_utc": [utc[pre[0]], utc[pre[-1]]] if pre else None,
               "late_totals_Pa": C_late.tolist(), "checks": chk, "seconds": time.time() - t0},
              open(str(out) + "_rows.json", "x"), indent=1)
    print("ng_source_map done in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

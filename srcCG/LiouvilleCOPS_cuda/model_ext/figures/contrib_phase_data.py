#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/contrib_phase_data.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/contrib_phase_data.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 70d91a5a797f6ad061fbfc9a04aba1ba21f7ecd25a43c2a20e3ec5188d695e7b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Per-facet contribution per unit area (NG density, RG flux; the fit's kernels and arithmetic) at snapshot base
#   rows, phase time means (before the NG gap, gap, late W10x5) and per-row region sums (lobe, distance shells), with
#   checks.
# RUN: jobs/contrib_phase_data.sbatch (normal128cores, CPU)
# ==== end of header ====
"""Data of the contribution-density maps and the region-share figures of one fitted state (2026-11-14, user: "make all
plots for the current output"; lead: contribution per unit surface area exactly as the fit computes it; shares by lobe and
by distance shells from the touchdown point). CPU.

Per base row (the fitting base's kernels, exact; frozen nonneg_contrib arithmetic as frame_contrib.py: the fit's FP32
record fields, FP64 math, the state's A, kappa, T; night / shadow rule; H2O + CO2) and record (facet e):
  NG density n_e = sum_s F_s,e (sqrt(pi) beta_s,e / 2) K_n,e [m^-3]  (the fit's NG weight incl. the shield side in K_n)
  RG flux  Phi_e = sum_s F_s,e [K_Phi,e + (V.b_RG)(sqrt(pi) beta_s,e / 2) K_V,e] [m^-2 s^-1]
  per unit area: n_e / A_e, Phi_e / A_e; pressures c_NG,e, c_RG,e [Pa] (record_contributions) for the shares.
Outputs (<out>/):
  snap_<HHMMSS>.npz      per snapshot row: face, ng, rg (per-area), ng_Pa, rg_Pa (per facet), totals
  phase_means.npz        dense per-facet time means of the per-area values: NG_pre, RG_pre (observation rows before the
                         NG gap), RG_gap (RG rows in the gap), NG_late, RG_late (the 95 late observations, W10x5 of the
                         grid rows as the fit); and the same for the pressure contributions (*_Pa) for the source areas
  regions.npz / .json    per base row: NG / RG pressure sums by lobe (small / big, neck plane of the A10 analysis) and by
                         distance shell from the touchdown point (0-30, 30-100, 100-300, 300-1000, 1000-3000, >3000 m);
                         late observations = W10x5 of the grid rows; checks against the contribution data's per-row model
  manifest.json          inputs, row classes, checks
"""

import argparse
import json
import multiprocessing as mp
import re
import sys
import time
from pathlib import Path

import numpy as np

SHELLS = (0.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, np.inf)
_G = {}


def chunk_task(rows):
    NC = _G["NC"]
    G = NC._G
    FCM = _G["FCM"]
    n = _G["nf"]
    acc = {k: np.zeros(n) for k in _G["acc_keys"]}
    res = {}
    for row, ip, wts in rows:
        kp, h, rec = NC.NM.read_kernel_full(Path(_G["base"]) / "kernels" / ("point=%04d.bin" % ip))
        if kp != ip:
            raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
        face, nH, nC, pH, pC = FCM.decompose(rec, h, G["A"], G["ACO2"], G["kappa"], G["T"], "fp32", NC)
        _, ((ng, rg),), _ = NC.record_contributions(rec, h, (G["A"],), "fp32", row)
        area = _G["area"][face]
        vng, vrg = (nH + nC) / area, (pH + pC) / area
        for key, w in wts.items():
            if key.endswith("_Pa"):
                acc[key][face] += w * (ng if key.startswith("NG") else rg)
            else:
                acc[key][face] += w * (vng if key.startswith("NG") else vrg)
        lb = _G["lobe"][face]
        sh = _G["shell"][face]
        r = {"ng": float(ng.sum()), "rg": float(rg.sum()), "n": float((nH + nC).sum()), "phi": float((pH + pC).sum()),
             "lobe_ng": np.bincount(lb, weights=ng, minlength=2), "lobe_rg": np.bincount(lb, weights=rg, minlength=2),
             "shell_ng": np.bincount(sh, weights=ng, minlength=len(SHELLS) - 1),
             "shell_rg": np.bincount(sh, weights=rg, minlength=len(SHELLS) - 1)}
        res[row] = r
        if row in _G["snap"]:
            np.savez(Path(_G["out"]) / ("snap_%s.npz" % _G["snap"][row]), face=face.astype(np.int32), ng=vng, rg=vrg,
                     ng_Pa=ng, rg_Pa=rg, total_n=r["n"], total_phi=r["phi"], total_ng_Pa=r["ng"], total_rg_Pa=r["rg"],
                     row=row, ipoint=ip)
    return acc, res


def main():
    ap = argparse.ArgumentParser()
    for k in ("--contrib-code", "--ft-dir", "--fit-dir", "--base", "--lobe", "--frame", "--slope-bins", "--out",
              "--frame-contrib-py"):
        ap.add_argument(k, required=True)
    ap.add_argument("--snap-rows", default="1079,1406,1528,1549")
    ap.add_argument("--procs", type=int, default=32)
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, a.contrib_code)
    import nonneg_contrib as NC
    import importlib.util
    spec = importlib.util.spec_from_file_location("frame_contrib_mod", a.frame_contrib_py)
    FCM = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(FCM)
    NM = NC.NM
    base, ft_dir, fit_dir = Path(a.base), Path(a.ft_dir), Path(a.fit_dir)
    fc, header, cols, rows = NM.base_table(base)
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
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("restart rows differ from the base rows")
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    utc = [r[cols["UTC"]] for r in rows]
    et = np.array([float(r[cols["et"]]) for r in rows])
    alt = np.array([float(r[cols["altitude"]]) for r in rows])
    refng = np.array([float(r[cols["refNG"]]) for r in rows])
    refrg = np.array([float(r[cols["refRG"]]) for r in rows])
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
    # row classes (as fov_share.py)
    nrow = len(rows)
    late, gr = list(layout.late_rows), list(layout.grid_rows)
    grid = np.zeros(nrow, bool)
    grid[gr] = True
    is_late = np.zeros(nrow, bool)
    is_late[late] = True
    obs = ((refng > 0) | (refrg > 0)) & ~grid
    ng_obs = [r for r in range(nrow) if obs[r] and refng[r] > 0 and not is_late[r]]
    r_gap = max(ng_obs, key=lambda r: et[r])
    pre = obs & ~is_late & (et <= et[r_gap])
    gap = obs & ~is_late & (et > et[r_gap])
    pre_ng = np.flatnonzero(pre & (refng > 0))
    pre_rg = np.flatnonzero(pre & (refrg > 0))
    gap_rg = np.flatnonzero(gap & (refrg > 0))
    wl = np.asarray(Wm.sum(axis=0)).ravel() / len(late)          # weight of grid row j in the late mean
    weights = {}
    for r in pre_ng:
        weights.setdefault(int(r), {}).update({"NG_pre": 1.0 / len(pre_ng), "NG_pre_Pa": 1.0 / len(pre_ng)})
    for r in pre_rg:
        weights.setdefault(int(r), {}).update({"RG_pre": 1.0 / len(pre_rg), "RG_pre_Pa": 1.0 / len(pre_rg)})
    for r in gap_rg:
        weights.setdefault(int(r), {}).update({"RG_gap": 1.0 / len(gap_rg), "RG_gap_Pa": 1.0 / len(gap_rg)})
    for j, r in enumerate(gr):
        if wl[j] > 0:
            weights.setdefault(int(r), {}).update({"NG_late": wl[j], "RG_late": wl[j], "NG_late_Pa": wl[j], "RG_late_Pa": wl[j]})
    acc_keys = ["NG_pre", "RG_pre", "RG_gap", "NG_late", "RG_late", "NG_pre_Pa", "RG_pre_Pa", "RG_gap_Pa", "NG_late_Pa",
                "RG_late_Pa"]
    # lobes and shells
    lobe = np.load(a.lobe).astype(np.intp)
    ptd = np.load(a.frame)["ptd"]
    d = np.linalg.norm(el["centroid"] - ptd, axis=1)
    shell = (np.searchsorted(np.asarray(SHELLS), d, side="right") - 1).astype(np.intp)
    snap_rows = [int(v) for v in a.snap_rows.split(",")]
    snap = {r: utc[r].split("T")[1].replace(":", "").replace(".", "_") for r in snap_rows}
    _G.update({"NC": NC, "FCM": FCM, "nf": int(el["nfaces"]), "acc_keys": acc_keys, "base": str(base), "area": el["area"],
               "lobe": lobe, "shell": shell, "snap": snap, "out": str(out)})
    tasks = [(r, int(ipoint[r]), weights.get(r, {})) for r in range(nrow)]
    tasks.sort(key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    chunks = [tasks[i::a.procs] for i in range(a.procs)]
    print("contrib_phase_data: %d rows (%d NG before the gap, %d RG before, %d RG in the gap, %d late from %d grid rows); "
          "snapshots %s (%.0f s)" % (nrow, len(pre_ng), len(pre_rg), len(gap_rg), len(late), int((wl > 0).sum()),
                                     snap, time.time() - t0), flush=True)
    acc = {k: np.zeros(int(el["nfaces"])) for k in acc_keys}
    R = {}
    with mp.get_context("fork").Pool(a.procs) as pool:
        for n, (acc_c, res) in enumerate(pool.imap_unordered(chunk_task, chunks), start=1):
            for k in acc_keys:
                acc[k] += acc_c[k]
            R.update(res)
            print("  chunk %d / %d (%.0f s)" % (n, len(chunks), time.time() - t0), flush=True)
    np.savez(out / "phase_means.npz", **acc)
    # regions per row; late observations = W10x5 of the grid rows
    tot = {k: np.array([R[r][k] for r in range(nrow)]) for k in ("ng", "rg", "n", "phi")}
    lob = {g: np.array([R[r]["lobe_" + g] for r in range(nrow)]) for g in ("ng", "rg")}
    shl = {g: np.array([R[r]["shell_" + g] for r in range(nrow)]) for g in ("ng", "rg")}
    mod = {k: v.copy() for k, v in tot.items()}
    for k in mod:
        mod[k][late] = Wm @ tot[k][gr]
    for g in ("ng", "rg"):
        lob[g][late] = Wm @ lob[g][gr]
        shl[g][late] = Wm @ shl[g][gr]
    # checks against the contribution data (per-row model surface pressure)
    SB = np.load(a.slope_bins, allow_pickle=False)
    pos = {int(i): j for j, i in enumerate(SB["iPoint"])}
    chk = {}
    for g, key in (("ng", "model_surface_NG_Pa"), ("rg", "model_surface_RG_Pa")):
        rr = [r for r in range(nrow) if int(ipoint[r]) in pos and (obs[r] or is_late[r])]
        ref = np.array([SB[key][pos[int(ipoint[r])]] for r in rr])
        mine = mod[g][rr]
        ok = ref != 0
        chk[g] = {"rows": len(rr), "max_rel_vs_contribution_data": float(np.max(np.abs(mine[ok] / ref[ok] - 1)))}
        nz = mod[g] != 0
        chk[g]["lobe_sum_vs_total_max_rel"] = float(np.max(np.abs(lob[g].sum(axis=1)[nz] / mod[g][nz] - 1)))
        chk[g]["shell_sum_vs_total_max_rel"] = float(np.max(np.abs(shl[g].sum(axis=1)[nz] / mod[g][nz] - 1)))
    for k in acc_keys:
        chk["phase_" + k] = {"facets_nonzero": int((acc[k] > 0).sum())}
    np.savez(out / "regions.npz", ipoint=ipoint, utc=np.array(utc), et=et, alt=alt, obs=obs, late=is_late, grid=grid,
             pre=pre, gap=gap, refng=refng, refrg=refrg, total_ng_Pa=mod["ng"], total_rg_Pa=mod["rg"],
             lobe_ng_Pa=lob["ng"], lobe_rg_Pa=lob["rg"], shell_ng_Pa=shl["ng"], shell_rg_Pa=shl["rg"],
             shells_m=np.array(SHELLS), lobe_names=np.array(["small lobe (head)", "big lobe (body)"]))
    json.dump({"definition": __doc__, "inputs": vars(a), "T": T, "kappa": float(S["kappa"]), "gauge_constants": gconst,
               "classes": {"NG_pre_rows": len(pre_ng), "RG_pre_rows": len(pre_rg), "RG_gap_rows": len(gap_rg),
                           "late_observations": len(late), "late_grid_rows_used": int((wl > 0).sum()),
                           "last_NG_row_before_gap": utc[r_gap]},
               "snapshots": {str(r): {"utc": utc[r], "altitude_m": float(alt[r]), "file": "snap_%s.npz" % snap[r],
                                      "NG_Pa": R[r]["ng"], "RG_Pa": R[r]["rg"], "n_m3": R[r]["n"], "phi_m2s": R[r]["phi"]}
                             for r in snap_rows},
               "checks": chk, "seconds": time.time() - t0}, open(out / "manifest.json", "x"), indent=1, default=float)
    print(json.dumps(chk, indent=1))
    print("contrib_phase_data done in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

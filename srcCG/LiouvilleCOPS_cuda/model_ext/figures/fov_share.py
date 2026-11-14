#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/fov_share.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/fov_share.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e07d3707c5d0be4aeb2f57735a65a1bce170026e6323f744631a3de259815e49
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: share of the modelled NG / RG surface signal from directions outside each gauge's IK field of view (frozen
#   nonneg_contrib record contributions; late rows W10x5); 2026-11-13 --cone-rg / --cone-ng: circular cones instead of
#   the IK rectangles.
# RUN: jobs/fov_share.sbatch (CONE_RG, CONE_NG); jobs 124731 (c5), 124789 (c9), 124885 (c9, circles 45 / 77.5 deg)
# ==== end of header ====
"""Share of the MODELLED NG and RG surface signal that comes from directions outside each gauge's IK field of view
(lead / user 2026-11-13, with the movies of the observable part: "for each of the 95 post-gap observations and for the
pre-gap rows, the share of the modelled RG signal (K_Phi-weighted, as the fit computes it) that comes from directions
outside the RG IK FOV (45 x 45 deg rectangle), and the share of the modelled NG signal outside the NG IK FOV. The fit
uses RGAcceptance = hemisphere, so this tells how different an IK-FOV movie is from what the model's RG integrates").

Per kernel record (point p, facet e) the NG and RG contributions are computed EXACTLY as the fitted model and the
contribution data (frozen model_ext nonneg_contrib.py record_contributions: the fit's FP32 records, F = q A with the
night/shadow rule, the fit's T and kappa, the fit_cops gauge constants; H2O + CO2); the record counts as inside a FOV
when the direction from the gauge (kernel header position) to the facet centroid passes LiouvilleCOPS
src/Geometry.cpp InsideFOV with that time's body-frame FOV bounds (cops_movies.py fov: fov.npz rows_*). Facets near the
gauge are subdivided by the model; this test uses the facet centroid for the whole record (stated in the output).
Late observations (95, >= 10:33:10 UTC): W10x5 of the 1-s grid rows, as the fit (out_k = sum_j W_kj out_j, the same
for the totals). The NG grain term is not a surface signal and is excluded.

Run (CPU job, jobs/fov_share.sbatch):
  fov_share.py --contrib-code <frozen nonneg_contrib dir> --ft-dir <c5/ft/final> --fit-dir <c5/fit/final> --base <base>
      --geometry <maps geometry of the chain> --fov <movies_work/fov.npz> --out <prefix> [--procs N]
Writes <prefix>.txt (per row) and <prefix>_summary.json / .txt (medians, ranges: before the NG gap, in the gap, the 95
late observations).
Added 2026-11-13 (detach-200 m movies; user: circular fields of view): --cone-rg DEG / --cone-ng DEG replace that gauge's
IK rectangle by a circular cone of DEG half-angle around its boresight (unit direction d inside when d.b >= cos DEG; the
boresight of the row time from fov.npz); default: the IK rectangles (unchanged).
"""

import argparse
import json
import multiprocessing as mp
import re
import sys
import time
from pathlib import Path

import numpy as np

_FOV = {}


def inside(d, bs, bounds):
    """Geometry.cpp InsideFOV for unit directions d (n, 3)."""
    m = (d @ bs) > 0.0
    b = np.asarray(bounds, np.float64)
    for i in range(len(b)):
        mi = np.cross(b[i], b[(i + 1) % len(b)])
        m &= (d @ mi) * float(np.dot(mi, bs)) >= 0.0
    return m


def row_task(task):
    row, ip = task
    NC = _FOV["NC"]
    G = NC._G
    kp, h, rec = NC.NM.read_kernel_full(Path(G["base"]) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
    face, ((ng, rg),), _ = NC.record_contributions(rec, h, (G["A"],), "fp32", row)
    x = np.asarray(h["sc"], np.float64)
    d = _FOV["cent"][face] - x[None, :]
    d /= np.linalg.norm(d, axis=1)[:, None]
    out = {"row": row, "ip": ip, "ng": float(ng.sum()), "rg": float(rg.sum())}
    for g, c in (("NG", ng), ("RG", rg)):
        if _FOV.get("cone_" + g) is not None:          # circular cone (2026-11-13 option)
            m = (d @ np.asarray(_FOV["bs_" + g][row], np.float64)) >= np.cos(np.radians(_FOV["cone_" + g]))
        else:
            m = inside(d, _FOV["bs_" + g][row], _FOV["bounds_" + g][row])
        out["in_" + g] = float(c[m].sum())
        out["recs_in_" + g] = int(m.sum())
    out["records"] = int(len(rec))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("--contrib-code", "--ft-dir", "--fit-dir", "--base", "--geometry", "--fov", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--cone-rg", type=float, default=None, help="RG: circular cone of this half-angle [deg] instead of the IK rectangle")
    ap.add_argument("--cone-ng", type=float, default=None, help="NG: circular cone of this half-angle [deg] instead of the IK rectangle")
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out)
    for suf in (".txt", "_summary.json", "_summary.txt"):
        if Path(str(out) + suf).exists():
            raise SystemExit("%s%s exists" % (out, suf))
    sys.path.insert(0, a.contrib_code)
    import nonneg_contrib as NC            # frozen copy (its nonneg_maps from the same folder)
    NM = NC.NM
    base, ft_dir, fit_dir = Path(a.base), Path(a.ft_dir), Path(a.fit_dir)
    # ---- the same setup as nonneg_contrib.main (constants, activity, late layout)
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
        use_ng, use_rg = np.asarray(z["use_ng"]), np.asarray(z["use_rg"])
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("restart rows differ from the base rows")
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    utc = [r[cols["UTC"]] for r in rows]
    et = np.array([float(r[cols["et"]]) for r in rows])            # time order (the UTC strings are not zero-padded)
    alt = np.array([float(r[cols["altitude"]]) for r in rows])
    refng = np.array([float(r[cols["refNG"]]) for r in rows])
    refrg = np.array([float(r[cols["refRG"]]) for r in rows])
    bm, cNG, cRG, gconst = NC.gauge_constants(fc, header, TRG)
    ftj = S["ft_json"]
    pprec = ftj["modules"]["positive_problem.py"]
    pp = NM.load_module(pprec["path"], "positive_problem_of_state", pprec["sha256"])
    ptab = ftj["fit"]["point_table"]
    layout = pp.LateLayout(pp.read_point_table(ptab["path"]), pp.read_run_points(base), ipoint)
    W = layout.W.tocsr()
    W.sort_indices()
    from cops_cuda import device_store
    NC._G.update({"base": str(base), "A": S["A"], "A_fit": S["A"], "ACO2": el["ACO2"], "bins": S["bins"].astype(np.intp),
                  "T": T, "kappa": S["kappa"], "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio, "AU": float(fc.AU),
                  "COS75": float(fc.COS75), "KT": float(fc.KT), "fp32_record_fields": device_store.fp32_record_fields,
                  "nfaces": el["nfaces"], "area": el["area"], "seen": S["seen"]})
    F = np.load(a.fov, allow_pickle=False)
    fov_utc = None
    _FOV.update({"NC": NC, "cent": el["centroid"],
                 "bs_NG": F["rows_bs_ng"], "bounds_NG": F["rows_bounds_ng"],
                 "bs_RG": F["rows_bs_rg"], "bounds_RG": F["rows_bounds_rg"],
                 "cone_RG": a.cone_rg, "cone_NG": a.cone_ng})
    if len(F["rows_bs_ng"]) != len(ipoint):
        raise SystemExit("fov.npz rows (%d) differ from the base rows (%d)" % (len(F["rows_bs_ng"]), len(ipoint)))
    tasks = sorted(((r, int(ipoint[r])) for r in range(len(ipoint))),
                   key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    print("fov_share: %d rows, %d processes" % (len(tasks), a.procs), flush=True)
    R = {}
    with mp.get_context("fork").Pool(a.procs) as pool:
        for n, res in enumerate(pool.imap_unordered(row_task, tasks, chunksize=1), start=1):
            R[res["row"]] = res
            if n % 300 == 0:
                print("  %d / %d rows (%.0f s)" % (n, len(tasks), time.time() - t0), flush=True)
    nrow = len(ipoint)
    tot = {k: np.array([R[r][k] for r in range(nrow)]) for k in ("ng", "rg", "in_NG", "in_RG")}
    late, gr = layout.late_rows, layout.grid_rows
    mod = {k: v.copy() for k, v in tot.items()}
    for k in mod:
        mod[k][late] = W @ tot[k][gr]
    out_ng = 1.0 - mod["in_NG"] / mod["ng"]
    out_rg = 1.0 - mod["in_RG"] / mod["rg"]
    # row classes: observation rows (refNG or refRG > 0, not grid), late = the 95 W10x5 rows; the NG gap
    grid = np.zeros(nrow, bool)
    grid[gr] = True
    is_late = np.zeros(nrow, bool)
    is_late[late] = True
    obs = ((refng > 0) | (refrg > 0)) & ~grid
    r_late0 = min(late, key=lambda r: et[r])
    ng_obs = [r for r in range(nrow) if obs[r] and refng[r] > 0 and not is_late[r]]
    r_gap = max(ng_obs, key=lambda r: et[r])           # last NG observation before the late window
    classes = {"before the NG gap (observation rows up to %s)" % utc[r_gap]:
               obs & ~is_late & (et <= et[r_gap]),
               "in the NG gap (RG observation rows after %s, before %s)" % (utc[r_gap], utc[r_late0]):
               obs & ~is_late & (et > et[r_gap]),
               "late window (95 observations, W10x5)": is_late}
    summ = {"definition": {
        "share_outside": "1 - (modelled surface signal from records whose facet-centroid direction is inside the gauge's "
                         "IK FOV) / (total modelled surface signal); NG grain excluded; late = W10x5 of the grid rows",
        "contributions": "frozen nonneg_contrib.py record_contributions (the fit's FP32 records, H2O + CO2)",
        "fov": ("LiouvilleCOPS Geometry.cpp InsideFOV on the getfov bounds rotated to 67P/C-G_CK at the row time"
                if a.cone_rg is None and a.cone_ng is None else
                "RG: %s; NG: %s (circular cone: d.b >= cos(half-angle), the getfov boresight at the row time)" % (
                    "IK rectangle (InsideFOV)" if a.cone_rg is None else "circular cone of %g deg half-angle" % a.cone_rg,
                    "IK rectangle (InsideFOV)" if a.cone_ng is None else "circular cone of %g deg half-angle" % a.cone_ng))},
        "inputs": vars(a), "classes": {}}
    lines = ["# fov_share.py: modelled surface-signal share from OUTSIDE the IK FOV, per row (iPoint UTC altitude_m "
             "class share_out_NG share_out_RG model_NG_Pa model_RG_Pa); late rows = W10x5", ""]
    for name, m in classes.items():
        ent = {"rows": int(m.sum())}
        for g, v in (("NG", out_ng), ("RG", out_rg)):
            sel = m & np.isfinite(v) & ((refng > 0) if g == "NG" else (refrg > 0))
            if sel.any():
                ent[g] = {"rows": int(sel.sum()), "median": float(np.median(v[sel])), "min": float(v[sel].min()),
                          "max": float(v[sel].max()), "mean": float(v[sel].mean())}
        summ["classes"][name] = ent
    for r in range(nrow):
        if not obs[r] and not is_late[r]:
            continue
        cl = [n for n, m in classes.items() if m[r]]
        lines.append("%d %s %.3f %s %.6f %.6f %.6e %.6e" % (ipoint[r], utc[r], alt[r], (cl[0].split(" (")[0] if cl else
                                                                                           "-").replace(" ", "_"),
                                                             out_ng[r], out_rg[r], mod["ng"][r], mod["rg"][r]))
    Path(str(out) + ".txt").write_text("\n".join(lines) + "\n")
    json.dump(summ, open(str(out) + "_summary.json", "x"), indent=1, default=str)
    T_ = ["# share of the modelled surface signal from directions OUTSIDE the IK FOV (RG: 45 x 45 deg rectangle; NG: "
          "77.5 x 77.5 deg); the fit integrates RG over the hemisphere and NG with the srccg_shield acceptance", ""]
    if a.cone_rg is not None or a.cone_ng is not None:
        T_[0] = ("# share of the modelled surface signal from directions OUTSIDE the field of view (RG: %s; NG: %s); the fit "
                 "integrates RG over the hemisphere and NG with the srccg_shield acceptance" % (
                     "45 x 45 deg IK rectangle" if a.cone_rg is None else "circular cone, %g deg half-angle" % a.cone_rg,
                     "77.5 x 77.5 deg IK rectangle" if a.cone_ng is None else "circular cone, %g deg half-angle" % a.cone_ng))
    for name, ent in summ["classes"].items():
        T_.append("%s: %d rows" % (name, ent["rows"]))
        for g in ("RG", "NG"):
            if g in ent:
                e = ent[g]
                T_.append("  %s: median %.4f, range %.4f-%.4f, mean %.4f (%d rows)" % (g, e["median"], e["min"], e["max"],
                                                                                    e["mean"], e["rows"]))
    Path(str(out) + "_summary.txt").write_text("\n".join(T_) + "\n")
    print("\n".join(T_))
    print("fov_share done in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/cliff_effect_data.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/cliff_effect_data.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b12ae130a18fc0ad22c37b01ab780aedecb71241c2c24c737233a0cda1190ae1
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: data of the slope-at-equal-illumination analysis of one state: per-record F, illumination class, slope bin,
#   observation-weight split and contributions with the fit's arithmetic, snapshot rows and phase means (2026-11-17)
# RUN: CPU Slurm job jobs/cliff_effect_data.sbatch, one per case (HOW-TO-RERUN.md of
#   figures_cluster/12M/cliff_effect_stage2/)
# ==== end of header ====
"""Does the source rate depend on the slope angle at the same illumination? Data of ONE state (2026-11-17, user: "test
several definitions: >40, >60, >80, as well as 0-40,40-50,50-60,60-70,70-80, 80-90; for the six cases individually. the
task is to find whether at the same illuminations condition the source rate depends on the slope angle"; lead's
specification and defaults).

    cliff_effect_data.py --spec <full_sets json> --key <set key> --out <dir> [--procs N]      (CPU job)

Model (LiouvilleCOPS src/Source.cpp, Liouville.h; base: SurfaceBC flat, SZAMode continuous, CastShadows on): per facet i
    F_i = A_i q(mu_i),  q = Qmin + (Qmax - Qmin) min(cos SZA, 1) on the lit side, Qmin at night (cos SZA <= 0) or in a
    ray-traced cast shadow;  NG contribution c_i = cNG bm F_i u_i Kn_i (u = 1/sqrt(T(mu)), T the fit's law), RG c_i =
    cRG F_i (KPhi_i + bm (V.b) u_i KV_i); Kn = sum over the visible sub-elements of dS/r^2 cos(theta)/pi W_NG(l).
Per kernel record the arithmetic is the fit's (nonneg_contrib.record_contributions of the frozen contrib code, FP32
records; checked equal); the H2O part is decomposed exactly: c/dS = const x A x I x G with
    NG: I = q u, G = Kn/dS;   RG: I = q, G = cRG (KPhi + bm (V.b) u KV)/dS / cRG-constant (the temperature enters G).
Illumination class of a record (FP64 file fields): lit SZA in 10 deg bins (0-10 ... 80-90) when cos SZA > 0 and not
shadowed, "shadowed" (cos SZA > 0, cast shadow), "night" (cos SZA <= 0). Slope: the fit's slope bins (nonneg_facets.npz
bin: 0-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-180); categories (unions) are formed by the figure script.
Surfaces: the gauge-seen facets of a base row (NG: Kn > 0; RG: KV > 0), at the snapshot rows 1079 (10:33:10), 1406
(10:37:31.282), 1528 (10:39:09.282), 1549 (10:39:26) and pooled over the phases (before the NG gap, the gap (RG), the late
window = W10x5 of the 1-s grid rows, weights = the mean W row of the late observations); the whole surface and the fit-seen
surface at the map times 10:33:10 / 10:39:26 from the maps' flux .vtu (F_H2O_per_m2_s, cos_SZA, illumination: the code
path of the production-rate maps).
Geometric split (P2) at the facet centroid: r = |x_gauge - centroid| (x_gauge = the kernel file's gauge position), l =
(x_gauge - centroid)/r, cos theta = l.n, W_NG = srccg_shield acceptance with the head-on weight (1 - l.b_NG) (base
settings; shield table of the movie preparation), W_RG = 1 for l.b_RG < 0 (hemisphere); the boresights b of the row from
the movie preparation's prep.npz (SPICE); visible fraction f = K/dS / (cos theta W / (pi r^2)) (K = Kn, KV): the remainder,
including partial visibility and the sub-element sampling (can exceed 1); records with cos theta <= 0 or W = 0 at the
centroid are counted and excluded from the split.
Writes <out>/cliff_data.npz (+ .json: inputs, checks). Refuses an existing --out.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
RUN = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5")
CX = Path("/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z")
MOVIE_WORK = RUN / "nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/movies_work_c9_detach300_pan"
SNAP_ROWS = (1079, 1406, 1528, 1549)
NJ, NB = 11, 7                     # illumination classes, slope bins
ND = 130                           # distance bins: log(r / 1 m) / log(1.1), r in 1 m ... 2e5 m
NT = 10                            # cos theta bins of 0.1
NCELL = NJ * ND * NT
GAUGES = ("NG", "RG")
T0 = time.time()
_G = {}


def say(t):
    print("cliff_effect_data [%6.0f s]: %s" % (time.time() - T0, t), flush=True)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def illum_class(cos64, shadow):
    lit = (cos64 > 0.0) & (shadow == 0)
    sza = np.degrees(np.arccos(np.clip(cos64, -1.0, 1.0)))
    j = np.where(lit, np.minimum((sza // 10.0).astype(np.int64), 8), np.where(cos64 > 0.0, 9, 10))
    return j


def new_acc():
    a = {}
    for g in GAUGES:
        a[g] = {"p1": np.zeros((4, NJ, NB)),        # area, A area, F area, active area
                "p2": np.zeros((6, NJ, NB)),        # area (split defined), inv_r2, cos, W, fvis, g   (area weighted)
                "p3": np.zeros((2, NB, NCELL)),     # area, c (H2O contribution) per cell
                "p4": np.zeros((6, NB)),            # c, area, A area, I area, G area, F area
                "n": np.zeros(4)}                    # records seen, split excluded (cos<=0), (W=0), records total
    a["co2"] = np.zeros(2)                          # CO2 contribution NG, RG (weighted)
    a["tot"] = np.zeros(2)                          # H2O contribution NG, RG (weighted)
    return a


def add(acc, b, w):
    for g in GAUGES:
        for k in acc[g]:
            acc[g][k] += w * b[g][k]
    acc["co2"] += w * b["co2"]
    acc["tot"] += w * b["tot"]


def row_stats(row):
    """Accumulators of one base row (weight 1) and the check data."""
    G = _G
    NC, NM = G["NC"], G["NM"]
    ip = int(G["ipoint"][row])
    kp, h, rec = NM.read_kernel_full(Path(G["base"]) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
    face, ((ng_ref, rg_ref),), _ = NC.record_contributions(rec, h, (G["A"],), "fp32", row)
    # the same arithmetic, H2O part separated (copy of record_contributions' lines; checked below)
    conv = {}
    f32 = NC._G["fp32_record_fields"](rec, NC._G["COS75"], conv)
    kn, kphi, kv, c0 = (f32[k].astype(np.float64) for k in ("kn", "kphi", "kv", "cos"))
    night = c0 <= 0.0
    qf = NM.q_factors(h["rh"], NC._G["helio"], NC._G["AU"])
    (qmin_h, qmax_h), (qmin_c, qmax_c) = qf["H2O"], qf["CO2"]
    vb = h["VdotB"]
    cz = np.where(c0 < 1.0, c0, 1.0)
    qh = np.where(night, qmin_h, qmin_h + (qmax_h - qmin_h) * cz)
    qc = np.where(night, qmin_c, qmin_c + (qmax_c - qmin_c) * cz) * NC._G["ACO2"][face]
    clipped = np.clip(c0, -1.0, 1.0)
    cc = np.where(night, -1.0, clipped)
    lit = cc >= NC._G["COS75"]
    x = np.where(lit, 1.0 / np.where(lit, cc, 1.0) - 1.0, 0.0)
    Tmax, T75, Tmin = NC._G["T"]
    tl = Tmax + NC._G["KT"] * (Tmax - T75) * x
    te = np.where(lit & (tl > Tmin), tl, Tmin)
    u = 1.0 / np.sqrt(te)
    kappa, bm, cNG = NC._G["kappa"], NC._G["bm"], NC._G["cNG"]
    cRGH, cRGC = NC._G["cRG"]["H2O"][row], NC._G["cRG"]["CO2"][row]
    A = G["A"][face]
    wnH, wvH, wpH = qh * kn, qh * vb * kv, qh * kphi
    wnC, wvC, phiC = qc * kn, qc * vb * kv, qc * kphi
    ng_h = cNG["H2O"] * bm["H2O"] * (wnH * A * u)
    ng_c = (kappa * cNG["CO2"] * bm["CO2"]) * (wnC * u)
    rg_h = cRGH * (wpH * A + bm["H2O"] * (wvH * A * u))
    rg_c = kappa * cRGC * (phiC + bm["CO2"] * (wvC * u))
    chk = {"row": row, "ip": ip, "records": int(len(rec)),
           "gauge_position_vs_prep_max_abs_m": float(np.max(np.abs(np.asarray(h["sc"]) - G["x_prep"][row]))),
           "ng_equal_record_contributions": bool(np.array_equal(ng_h + ng_c, ng_ref)),
           "rg_equal_record_contributions": bool(np.array_equal(rg_h + rg_c, rg_ref)),
           "ng_total": float(np.sum(ng_ref)), "rg_total": float(np.sum(rg_ref))}
    dS = G["area"][face]
    F = A * qh
    jcl = illum_class(np.asarray(rec["cos"], np.float64), np.asarray(rec["shadow"]))
    b = G["bins"][face]
    # geometry at the centroid
    xg = np.asarray(h["sc"], np.float64)
    d = xg[None, :] - G["centroid"][face]
    r = np.sqrt(np.einsum("ij,ij->i", d, d))
    l = d / r[:, None]
    cth = np.einsum("ij,ij->i", l, G["normal"][face])
    bng, brg = G["b_ng"][row], G["b_rg"][row]
    mu = l @ bng
    lrg = l @ brg
    sf = G["shield"]
    xs = np.clip(mu, 0.0, 1.0) * (len(sf) - 1)
    i0 = np.minimum(xs.astype(np.int64), len(sf) - 2)
    wsf = xs - i0
    shield = (1.0 - wsf) * sf[i0] + wsf * sf[i0 + 1]
    W_ng = np.where(mu < 0.0, 1.0 - mu, np.where(lrg < 0.0, shield, 0.0))
    W_rg = np.where(lrg < 0.0, 1.0, 0.0)
    out = {"co2": np.array([np.sum(ng_c), np.sum(rg_c)]), "tot": np.array([np.sum(ng_h), np.sum(rg_h)])}
    jb = jcl * NB + b
    ldist = np.clip(np.floor(np.log(np.maximum(r, 1.0)) / np.log(1.1)), 0, ND - 1).astype(np.int64)
    for g, K, W, c, I, Gm in (("NG", kn, W_ng, ng_h, qh * u, kn / dS),
                              ("RG", kv, W_rg, rg_h, qh, cRGH * (kphi + bm["H2O"] * vb * kv * u) / dS)):
        s = K > 0.0
        a = {"n": np.zeros(4)}
        js, bs, dss = jb[s], b[s], dS[s]
        p1 = np.zeros((4, NJ * NB))
        p1[0] = np.bincount(js, dss, NJ * NB)
        p1[1] = np.bincount(js, dss * A[s], NJ * NB)
        p1[2] = np.bincount(js, dss * F[s], NJ * NB)
        p1[3] = np.bincount(js, dss * (A[s] > 0.0), NJ * NB)
        a["p1"] = p1.reshape(4, NJ, NB)
        ideal = cth[s] * W[s] / (np.pi * r[s] ** 2)
        ok = (cth[s] > 0.0) & (W[s] > 0.0)
        a["n"][:] = [s.sum(), np.sum(cth[s] <= 0.0), np.sum((cth[s] > 0.0) & (W[s] <= 0.0)), len(rec)]
        gg = K[s] / dss
        fv = np.where(ok, gg / np.where(ok, ideal, 1.0), 0.0)
        jo, do = js[ok], dss[ok]
        p2 = np.zeros((6, NJ * NB))
        p2[0] = np.bincount(jo, do, NJ * NB)
        p2[1] = np.bincount(jo, do / r[s][ok] ** 2, NJ * NB)
        p2[2] = np.bincount(jo, do * cth[s][ok], NJ * NB)
        p2[3] = np.bincount(jo, do * W[s][ok], NJ * NB)
        p2[4] = np.bincount(jo, do * fv[ok], NJ * NB)
        p2[5] = np.bincount(jo, do * gg[ok], NJ * NB)
        a["p2"] = p2.reshape(6, NJ, NB)
        tb = np.clip(np.floor(cth[s] * NT), 0, NT - 1).astype(np.int64)
        cell = (jcl[s] * ND + ldist[s]) * NT + tb
        idx = bs * NCELL + cell
        p3 = np.zeros((2, NB * NCELL))
        p3[0] = np.bincount(idx, dss, NB * NCELL)
        p3[1] = np.bincount(idx, c[s], NB * NCELL)
        a["p3"] = p3.reshape(2, NB, NCELL)
        p4 = np.zeros((6, NB))
        p4[0] = np.bincount(bs, c[s], NB)
        p4[1] = np.bincount(bs, dss, NB)
        p4[2] = np.bincount(bs, dss * A[s], NB)
        p4[3] = np.bincount(bs, dss * I[s], NB)
        p4[4] = np.bincount(bs, dss * Gm[s], NB)
        p4[5] = np.bincount(bs, dss * F[s], NB)
        a["p4"] = p4
        out[g] = a
    # F of the records vs the maps' F (snapshot rows at the map times)
    fm = G["maps_F"].get(row)
    if fm is not None:
        Fm, Jm = fm
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.abs(F - Fm[face]) / np.abs(Fm[face])
        both = (Fm[face] > 0) | (F > 0)
        chk["F_records_vs_maps_max_rel"] = float(np.nanmax(rel[both])) if both.any() else 0.0
        chk["F_records_vs_maps_zero_pattern_equal"] = bool(np.array_equal(F == 0, Fm[face] == 0))
        jm = Jm[face]
        jr = np.where(jcl <= 8, 0, np.where(jcl == 9, 2, 1))       # day 0, night 1, cast shadow 2 (maps codes)
        chk["illumination_class_vs_maps_mismatch_records"] = int(np.sum(jr != jm))
    return row, out, chk


def chunk_work(task):
    rows, weights = task
    acc = new_acc()
    chks = []
    for r, w in zip(rows, weights):
        _, b, chk = row_stats(r)
        add(acc, b, w)
        chks.append({k: chk[k] for k in ("row", "ng_equal_record_contributions", "rg_equal_record_contributions")})
    return acc, chks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    ap.add_argument("--phases", default="all", choices=("all", "none"), help="none: snapshot rows only (test)")
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    sp = json.load(open(a.spec))["sets"][a.key]
    ft, fit, maps, phase = Path(sp["ft"]), Path(sp["fit"]), Path(sp["maps"]), Path(sp["phase"])
    sys.path.insert(0, str(CX))
    sys.path.insert(0, str(APP))
    import nonneg_contrib as NC
    NM = NC.NM
    base = RUN / "base"
    fc, header, cols, rows = NM.base_table(base)
    helio = fc.header_value(header, "HeliocentricModel")
    el = NM.read_elements(base / "elements.bin")
    with np.load(ft / "finetune_state.npz", allow_pickle=False) as z:
        A = np.asarray(z["A_final"], np.float64)
    with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
        seen_fit = np.asarray(z["seen"]).astype(bool)
        bins = np.asarray(z["bin"]).astype(np.int64)
    with np.load(fit / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
        kappa = float(np.asarray(z["kappa"]).reshape(-1)[0])
        point_ids = np.asarray(z["point_ids"], np.int64)
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("restart rows differ from the base rows")
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    bm, cNG, cRG, gconst = NC.gauge_constants(fc, header, TRG)
    from cops_cuda import device_store
    NC._G.update({"base": str(base), "A": A, "A_fit": A, "ACO2": el["ACO2"], "bins": bins.astype(np.intp), "T": T,
                  "kappa": kappa, "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio, "AU": float(fc.AU),
                  "COS75": float(fc.COS75), "KT": float(fc.KT), "fp32_record_fields": device_store.fp32_record_fields,
                  "nfaces": el["nfaces"], "area": el["area"], "seen": seen_fit})
    prep = np.load(MOVIE_WORK / "prep.npz", allow_pickle=False)
    if not np.array_equal(np.asarray(prep["rows_ipoint"]).astype(np.int64), ipoint):
        raise SystemExit("prep.npz rows differ from the base rows")
    shield = np.loadtxt(MOVIE_WORK / "shield_table.txt")
    # whole surface (maps' .vtu at the map times)
    man = json.load(open(maps / "maps_manifest.json"))
    sys.path.insert(0, man["code_dir"])
    import importlib.util
    spec_ = importlib.util.spec_from_file_location("nm_of_maps", str(Path(man["code_dir"]) / "nonneg_maps.py"))
    nmm = importlib.util.module_from_spec(spec_)
    spec_.loader.exec_module(nmm)
    area = el["area"]
    whole, maps_F = {}, {}
    for t in man["times"]:
        lab = t["utc"][-8:]
        if lab not in ("10:33:10", "10:39:26"):
            continue
        rd, _, _ = nmm.read_vtu(maps / t["vtu"])
        F = np.asarray(rd["F_H2O_per_m2_s"], np.float64)
        cs = np.asarray(rd["cos_SZA"], np.float64)
        il = np.asarray(rd["illumination"]).astype(np.int64)     # 0 day, 1 night, 2 cast shadow
        sza = np.degrees(np.arccos(np.clip(cs, -1.0, 1.0)))
        j = np.where(il == 0, np.minimum((sza // 10.0).astype(np.int64), 8), np.where(il == 2, 9, 10))
        jb = j * NB + bins
        tabs = {}
        for name, m in (("all", np.ones(len(A), bool)), ("fit_seen", seen_fit)):
            p = np.zeros((4, NJ * NB))
            p[0] = np.bincount(jb[m], area[m], NJ * NB)
            p[1] = np.bincount(jb[m], area[m] * A[m], NJ * NB)
            p[2] = np.bincount(jb[m], area[m] * F[m], NJ * NB)
            p[3] = np.bincount(jb[m], area[m] * (A[m] > 0), NJ * NB)
            tabs[name] = p.reshape(4, NJ, NB)
        with np.errstate(divide="ignore", invalid="ignore"):
            whole[lab] = tabs
        rowi = {"10:33:10": 1079, "10:39:26": 1549}[lab]
        if not rows[rowi][cols["UTC"]].endswith(lab):
            raise SystemExit("base row %d is %s, not the map time %s" % (rowi, rows[rowi][cols["UTC"]], lab))
        maps_F[rowi] = (F, il)
        del rd
    say("whole surface tables: %s" % list(whole))
    _G.update({"NC": NC, "NM": NM, "base": str(base), "ipoint": ipoint, "A": A, "bins": bins, "area": area,
               "centroid": el["centroid"], "normal": el["normal"], "b_ng": np.asarray(prep["rows_b_ng"]),
               "b_rg": np.asarray(prep["rows_b_rg"]), "x_prep": np.asarray(prep["rows_x"]), "shield": shield,
               "maps_F": maps_F})
    # gauge position check: kernel header vs prep rows_x
    R = np.load(phase / "regions.npz", allow_pickle=False)
    if not np.array_equal(np.asarray(R["ipoint"]).astype(np.int64), ipoint):
        raise SystemExit("regions.npz rows differ from the base rows")
    ftj = json.load(open(ft / "fine_tuning.json"))
    pprec = ftj["modules"]["positive_problem.py"]
    pp = NM.load_module(pprec["path"], "positive_problem_of_state", pprec["sha256"])
    layout = pp.LateLayout(pp.read_point_table(ftj["fit"]["point_table"]["path"]), pp.read_run_points(base), ipoint)
    Wm = layout.W.tocsr()
    late_rows, grid_rows = np.asarray(layout.late_rows), np.asarray(layout.grid_rows)
    wgrid = np.asarray(Wm.sum(axis=0)).ravel() / len(late_rows)
    if len(wgrid) != len(grid_rows):
        raise SystemExit("W has %d columns, %d grid rows" % (len(wgrid), len(grid_rows)))
    phases = {"pre_gap": (np.nonzero(R["pre"] & R["obs"])[0], None),
              "gap": (np.nonzero(R["gap"] & R["obs"])[0], None),
              "late_W10x5": (grid_rows, wgrid)}
    if a.phases == "none":
        phases = {}
    say("phases: %s; late W weights sum %.6f" % ({k: len(v[0]) for k, v in phases.items()}, float(np.sum(wgrid))))
    pool = mp.get_context("fork").Pool(a.procs)
    res = {"snap": {}, "phase": {}}
    checks = {"rows": {}}
    try:
        for r, b, chk in pool.imap_unordered(row_stats, SNAP_ROWS):
            res["snap"][r] = b
            checks["rows"][r] = chk
            say("snapshot row %d: %s" % (r, json.dumps({k: v for k, v in chk.items() if k != "row"})))
        for name, (rws, wts) in phases.items():
            if wts is None:
                wts = np.full(len(rws), 1.0 / len(rws))
            n = max(1, min(len(rws), a.procs * 2))
            tasks = [(list(rws[i::n]), list(wts[i::n])) for i in range(n)]
            acc = new_acc()
            eq = []
            for ac, ch in pool.imap_unordered(chunk_work, tasks):
                add(acc, ac, 1.0)
                eq += ch
            res["phase"][name] = acc
            checks["phase_%s_records_equal_record_contributions" % name] = bool(all(
                c["ng_equal_record_contributions"] and c["rg_equal_record_contributions"] for c in eq))
            say("phase %s: %d rows, contributions equal: %s" % (name, len(eq), checks["phase_%s_records_equal_record_contributions" % name]))
    finally:
        pool.close()
        pool.join()
    # per-row totals vs the contribution-density snapshots of the full set (contrib_phase_data.py, instantaneous rows)
    for r in SNAP_ROWS:
        tag = rows[r][cols["UTC"]].split("T")[1].replace(":", "").replace(".", "_")
        with np.load(phase / ("snap_%s.npz" % tag), allow_pickle=False) as z:
            tn, tr_ = float(z["total_ng_Pa"]), float(z["total_rg_Pa"])
            if int(z["row"]) != r:
                raise SystemExit("snap_%s.npz is row %d, not %d" % (tag, int(z["row"]), r))
        tot = checks["rows"][r]
        checks["rows"][r]["ng_vs_contrib_phase_snapshot_rel"] = abs(tot["ng_total"] / tn - 1.0)
        checks["rows"][r]["rg_vs_contrib_phase_snapshot_rel"] = abs(tot["rg_total"] / tr_ - 1.0)
    out.mkdir(parents=True)
    arrs = {}
    for lab, tabs in whole.items():
        for name, p in tabs.items():
            arrs["whole_%s_%s_p1" % (lab.replace(":", ""), name)] = p
    for r, b in res["snap"].items():
        for g in GAUGES:
            for k, v in b[g].items():
                arrs["snap%d_%s_%s" % (r, g, k)] = v
        arrs["snap%d_tot" % r], arrs["snap%d_co2" % r] = b["tot"], b["co2"]
    for name, b in res["phase"].items():
        for g in GAUGES:
            for k, v in b[g].items():
                arrs["phase_%s_%s_%s" % (name, g, k)] = v
        arrs["phase_%s_tot" % name], arrs["phase_%s_co2" % name] = b["tot"], b["co2"]
    np.savez_compressed(out / "cliff_data.npz", **arrs)
    meta = {"key": a.key, "spec": sp, "inputs": {"ft": str(ft), "fit": str(fit), "maps": str(maps), "phase": str(phase),
            "contrib_code": str(CX), "prep": str(MOVIE_WORK / "prep.npz"), "shield_table": str(MOVIE_WORK / "shield_table.txt"),
            "shield_table_bins": int(len(shield))},
            "T": T, "kappa": kappa, "snap_rows": list(SNAP_ROWS),
            "snap_utc": {int(r): rows[r][cols["UTC"]] for r in SNAP_ROWS},
            "phases": {k: int(len(v[0])) for k, v in phases.items()}, "NJ": NJ, "NB": NB, "ND": ND, "NT": NT,
            "illumination_classes": ["SZA %d-%d" % (10 * k, 10 * k + 10) for k in range(9)] + ["shadowed", "night"],
            "bins": ["0-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-180"],
            "distance_bin": "floor(ln(max(r,1 m))/ln 1.1), 0..%d" % (ND - 1), "checks": checks,
            "command": sys.argv, "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "wall_seconds": time.time() - T0}
    json.dump(meta, open(out / "cliff_data.json", "x"), indent=1, default=str)
    say("done -> %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

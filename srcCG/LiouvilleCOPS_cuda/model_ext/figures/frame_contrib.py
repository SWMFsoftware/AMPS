#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/frame_contrib.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/frame_contrib.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5fa84399495c2145e9b55af839c7a4d5eca812b535745eb8433ede92f6ef0776
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Per-facet contribution per unit area (NG density, RG flux) of the fitted model at every movie frame from the
#   frames base kernels (frozen nonneg_contrib arithmetic, FP32 records); verification against record_contributions,
#   the frames base COPS_model.dat (base source) and the base rows; colour range.
# RUN: jobs/frame_contrib.sbatch (normal128cores, CPU)
# ==== end of header ====
"""Per-facet contribution per unit surface area at every movie frame, exactly as the fitted model computes it
(2026-11-14, user: "colour = contribution per unit surface area, exactly as the fitted model computes it").

Kernels: the fitting base built at the frame times (frames_base.sbatch: the lategrid base's boundary.in and executable,
ObservationFile = the frame times; point k = frame k); --frames-base-tail: a second build for the frames after the end
of the RG temperature table (TemperatureRGMode constant; frames >= its first point are read from it; the overlap frame
is compared record by record). Per record (frame k, facet e) with the frozen nonneg_contrib
arithmetic (the fit's FP32 record fields, FP64 math, the c9 state: activity A, CO2 scale kappa, temperature law T; the
night / shadow rule; H2O + CO2):
  NG  n_e   = sum_s F_s,e (sqrt(pi) beta_s,e / 2) K_n,e                          [m^-3]   (density at NG; the fit's
             NG weight incl. 1 - l.b and the shield fraction SF are inside K_n)
  RG  Phi_e = sum_s F_s,e [K_Phi,e + (V.b_RG)(sqrt(pi) beta_s,e / 2) K_V,e]      [m^-2 s^-1] (flux into RG)
  per unit area: n_e / A_e [m^-5], Phi_e / A_e [m^-4 s^-1] (A_e = facet area of elements.bin).
Facets without a record (not visible from the gauge, facing away) contribute 0, as in the kernels. The pressures of the
fit follow as p_NG = sum_s cNG_s n_s, p_RG = sum_s cRG_s Phi_s (fit_cops gauge constants, T_RG of the frame point).
Writes <out>/fNNNNN.npz (face int32, ng, rg float32 per-area values; totals float64), <out>/frame_totals.tsv,
<out>/verification.json (checks below), <out>/color_range.json (the fixed log range per gauge: top = 99.99th percentile
of the per-area values of the facets with a contribution in the last --top-frames frames; bottom = the per-area level
above which the facets carry --bottom-share of the gauge's surface signal at frame 0).
Checks: (1) for every frame the sums of the decomposition times the gauge constants equal the frozen
record_contributions Pa sums; (2) with the base's own source (elements.bin AH2O / ACO2, kappa 1, the analytical
temperature law of the boundary file) and the FP64 records, the sums equal the frames base's COPS_model.dat columns
n_H2O, n_CO2, PhiRG_H2O, PhiRG_CO2 (and pNG, pRG); (3) frames at the times of base rows (frame 0 = row 0 within 1 us,
frame 1799 = row 1549): the kernel records equal the base's, and the c9 sums equal record_contributions on the base row
and the contribution data's model_surface (row 0).
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


def decompose(rec, h, A, ACO2, kappa, T, prec, NC):
    """the record_contributions arithmetic split into species densities and fluxes (no gauge constants)."""
    G = NC._G
    NM = NC.NM
    face = rec["face"].astype(np.int64)
    if prec == "fp32":
        conv = {}
        f32 = G["fp32_record_fields"](rec, G["COS75"], conv)
        kn, kphi, kv, c0 = (f32[k].astype(np.float64) for k in ("kn", "kphi", "kv", "cos"))
        night = c0 <= 0.0
    else:
        kn, kphi, kv, c0 = (np.asarray(rec[k], np.float64) for k in ("Kn", "KPhi", "KV", "cos"))
        night = (c0 <= 0.0) | (rec["shadow"] != 0)
    qf = NM.q_factors(h["rh"], G["helio"], G["AU"])
    (qmin_h, qmax_h), (qmin_c, qmax_c) = qf["H2O"], qf["CO2"]
    vb = h["VdotB"]
    cz = np.where(c0 < 1.0, c0, 1.0)
    qh = np.where(night, qmin_h, qmin_h + (qmax_h - qmin_h) * cz)
    qc = np.where(night, qmin_c, qmin_c + (qmax_c - qmin_c) * cz) * ACO2[face]
    wnH, wvH, wpH = qh * kn, qh * vb * kv, qh * kphi
    wnC, wvC, phiC = qc * kn, qc * vb * kv, qc * kphi
    clipped = np.clip(c0, -1.0, 1.0)
    c = np.where(night, -1.0, clipped)
    lit = c >= G["COS75"]
    x = np.where(lit, 1.0 / np.where(lit, c, 1.0) - 1.0, 0.0)
    Tmax, T75, Tmin = T
    tl = Tmax + G["KT"] * (Tmax - T75) * x
    te = np.where(lit & (tl > Tmin), tl, Tmin)
    u = 1.0 / np.sqrt(te)
    bm = G["bm"]
    af = A[face]
    nH = bm["H2O"] * (wnH * af * u)
    nC = kappa * bm["CO2"] * (wnC * u)
    pH = wpH * af + bm["H2O"] * (wvH * af * u)
    pC = kappa * (phiC + bm["CO2"] * (wvC * u))
    return face, nH, nC, pH, pC


def fb_of(k):
    return _G["fb_tail"] if (_G.get("fb_tail") and k >= _G["tail_first"]) else _G["fb"]


def frame_task(k):
    NC = _G["NC"]
    G = NC._G
    kp, h, rec = NC.NM.read_kernel_full(Path(fb_of(k)) / "kernels" / ("point=%04d.bin" % k))
    if kp != k:
        raise RuntimeError("kernel file of point %d stores %d" % (k, kp))
    face, nH, nC, pH, pC = decompose(rec, h, G["A"], G["ACO2"], G["kappa"], G["T"], "fp32", NC)
    _, ((ng, rg),), _ = NC.record_contributions(rec, h, (G["A"],), "fp32", k)
    cNG, cRG = G["cNG"], G["cRG"]
    ng_dec = cNG["H2O"] * nH.sum() + cNG["CO2"] * nC.sum()
    rg_dec = cRG["H2O"][k] * pH.sum() + cRG["CO2"][k] * pC.sum()
    area = _G["area"][face]
    vng = (nH + nC) / area
    vrg = (pH + pC) / area
    out = {"k": k, "records": int(len(rec)), "n": float((nH + nC).sum()), "phi": float((pH + pC).sum()),
           "nH": float(nH.sum()), "nC": float(nC.sum()), "pH": float(pH.sum()), "pC": float(pC.sum()),
           "ng_Pa": float(ng.sum()), "rg_Pa": float(rg.sum()), "ng_Pa_dec": float(ng_dec), "rg_Pa_dec": float(rg_dec)}
    # check (2): the base's own source, FP64 records
    _, bH, bC, bpH, bpC = decompose(rec, h, _G["AH2O"], _G["ACO2"], 1.0, _G["T_base"], "fp64", NC)
    out.update({"base_nH": float(bH.sum()), "base_nC": float(bC.sum()), "base_pH": float(bpH.sum()),
                "base_pC": float(bpC.sum())})
    np.savez(Path(_G["out"]) / ("f%05d.npz" % k), face=face.astype(np.int32), ng=vng.astype(np.float32),
             rg=vrg.astype(np.float32), total_ng=out["n"], total_rg=out["phi"], total_ng_Pa=out["ng_Pa"],
             total_rg_Pa=out["rg_Pa"])
    if k in _G["keep"]:
        _G_keep = {"face": face, "vng": vng, "vrg": vrg, "area": area}
        np.savez(Path(_G["out"]) / ("full_f%05d.npz" % k), **_G_keep)
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("--contrib-code", "--ft-dir", "--fit-dir", "--base", "--frames-base", "--out", "--slope-bins"):
        ap.add_argument(k, required=True)
    ap.add_argument("--procs", type=int, default=48)
    ap.add_argument("--frames-base-tail", default="")
    ap.add_argument("--top-frames", type=int, default=10)
    ap.add_argument("--top-percentile", type=float, default=99.99)
    ap.add_argument("--bottom-share", type=float, default=0.90)
    ap.add_argument("--row-pairs", default="0:0,1799:1549", help="frame:base row pairs at the same time")
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, a.contrib_code)
    import nonneg_contrib as NC
    NM = NC.NM
    base, fb, ft_dir, fit_dir = Path(a.base), Path(a.frames_base), Path(a.ft_dir), Path(a.fit_dir)
    fc, header, cols, rows = NM.base_table(base)
    helio = fc.header_value(header, "HeliocentricModel")
    mrun = re.search(r"^# run: (\S+)$", (fit_dir / "fit_result.txt").read_text(), re.M)
    if mrun is None or Path(mrun.group(1)).resolve() != base.resolve():
        raise SystemExit("the fit did not run on %s" % base)
    el = NM.read_elements(base / "elements.bin")
    import hashlib
    sh0 = hashlib.sha256((base / "elements.bin").read_bytes()).hexdigest()
    for b_ in [fb] + ([Path(a.frames_base_tail)] if a.frames_base_tail else []):
        if hashlib.sha256((b_ / "elements.bin").read_bytes()).hexdigest() != sh0:
            raise SystemExit("elements.bin of %s differs from the base's" % b_)
    S = NM.state_activity(ft_dir, fit_dir, el, a.procs)
    if not S["checks"]["passed"]:
        raise SystemExit("state checks failed")
    with np.load(fit_dir / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
    fcb, header_f, cols_f, rows_f = NM.base_table(fb)
    tail_first = None
    if a.frames_base_tail:                     # rows of the tail build replace the frames from its first point on
        _, header_t, cols_t, rows_t = NM.base_table(Path(a.frames_base_tail))
        if cols_t != cols_f:
            raise SystemExit("the tail build's table columns differ")
        ip_t = [int(float(r[cols_t["iPoint"]])) for r in rows_t]
        tail_first = ip_t[0] + 1                # the first point is the overlap frame (compared, taken from the first build)
        ip_f = [int(float(r[cols_f["iPoint"]])) for r in rows_f]
        if ip_f[-1] < ip_t[0]:
            raise SystemExit("no overlap frame between the builds")
        overlap_row_t = rows_t[0]
        rows_f = [r for r, i in zip(rows_f, ip_f) if i < tail_first] + [r for r, i in zip(rows_t, ip_t) if i >= tail_first]
    for key in ("HeliocentricModel", "TemperatureNG", "RGPressureModel"):
        if fc.header_value(header, key) != fc.header_value(header_f, key):
            raise SystemExit("header %s differs between the bases" % key)
    TRG = np.array([float(r[cols_f["T_RG"]]) for r in rows_f])
    bm, cNG, cRG, gconst = NC.gauge_constants(fc, header_f, TRG)
    T_base = tuple(float(fc.header_value(header_f, k)) for k in ("AnalyticalTmax", "AnalyticalT75", "AnalyticalTmin"))
    from cops_cuda import device_store
    NC._G.update({"base": str(fb), "A": S["A"], "A_fit": S["A"], "ACO2": el["ACO2"], "bins": S["bins"].astype(np.intp),
                  "T": T, "kappa": S["kappa"], "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio, "AU": float(fc.AU),
                  "COS75": float(fc.COS75), "KT": float(fc.KT), "fp32_record_fields": device_store.fp32_record_fields,
                  "nfaces": el["nfaces"], "area": el["area"], "seen": S["seen"]})
    nfr = len(rows_f)
    keep = {0} | set(range(nfr - a.top_frames, nfr))
    _G.update({"NC": NC, "fb": str(fb), "out": str(out), "area": el["area"], "AH2O": el["AH2O"], "ACO2": el["ACO2"],
               "T_base": T_base, "keep": keep, "fb_tail": a.frames_base_tail, "tail_first": tail_first})
    if [int(float(r[cols_f["iPoint"]])) for r in rows_f] != list(range(len(rows_f))):
        raise SystemExit("frame rows are not 0..N-1")
    print("frame_contrib: %d frames, c9 T %s kappa %.6g; base T %s (%.0f s)" % (nfr, T, S["kappa"], T_base, time.time() - t0),
          flush=True)
    R = {}
    order = sorted(range(nfr), key=lambda k: -int(float(rows_f[k][cols_f["nVisible"]])))
    with mp.get_context("fork").Pool(a.procs) as pool:
        for n, res in enumerate(pool.imap_unordered(frame_task, order, chunksize=1), start=1):
            R[res["k"]] = res
            if n % 300 == 0:
                print("  %d / %d frames (%.0f s)" % (n, nfr, time.time() - t0), flush=True)
    # verification
    col = lambda name: np.array([float(r[cols_f[name]]) for r in rows_f])   # noqa: E731
    ver = {}
    rel = lambda a_, b_: float(np.max(np.abs(a_ / b_ - 1.0)))   # noqa: E731
    ng_Pa = np.array([R[k]["ng_Pa"] for k in range(nfr)])
    rg_Pa = np.array([R[k]["rg_Pa"] for k in range(nfr)])
    ver["1_decomposition_vs_record_contributions_max_rel"] = {
        "NG": rel(np.array([R[k]["ng_Pa_dec"] for k in range(nfr)]), ng_Pa),
        "RG": rel(np.array([R[k]["rg_Pa_dec"] for k in range(nfr)]), rg_Pa)}
    v2 = {}
    for key, cname in (("base_nH", "n_H2O"), ("base_nC", "n_CO2"), ("base_pH", "PhiRG_H2O"), ("base_pC", "PhiRG_CO2")):
        mine = np.array([R[k][key] for k in range(nfr)])
        tab = col(cname)
        ok = tab != 0
        v2[cname] = {"max_rel": rel(mine[ok], tab[ok]), "median_rel": float(np.median(np.abs(mine[ok] / tab[ok] - 1))),
                     "rows": int(ok.sum())}
    pNG_b = np.array([cNG["H2O"] * R[k]["base_nH"] + cNG["CO2"] * R[k]["base_nC"] for k in range(nfr)])
    pRG_b = np.array([cRG["H2O"][k] * R[k]["base_pH"] + cRG["CO2"][k] * R[k]["base_pC"] for k in range(nfr)])
    v2["pNG"] = {"max_rel": rel(pNG_b, col("pNG"))}
    v2["pRG"] = {"max_rel": rel(pRG_b, col("pRG"))}
    v2["note"] = ("COPS_model.dat prints 7 significant digits (format .6e); max_rel at that level means agreement to "
                  "print precision; base source: elements.bin AH2O, ACO2, kappa 1, temperature %s, FP64 records" % (T_base,))
    ver["2_base_source_vs_frames_base_COPS_model"] = v2
    # (3) frames at base-row times
    SB = np.load(a.slope_bins, allow_pickle=False)
    sb_ip = {int(i): j for j, i in enumerate(SB["iPoint"])}
    v3 = []
    for pair in a.row_pairs.split(","):
        kf, rb = (int(x) for x in pair.split(":"))
        kp, hb, recb = NM.read_kernel_full(base / "kernels" / ("point=%04d.bin" % rb))
        _, hf, recf = NM.read_kernel_full(fb / "kernels" / ("point=%04d.bin" % kf))
        same = (len(recb) == len(recf)) and all(np.array_equal(recb[n], recf[n]) for n in recb.dtype.names)
        Gs = NC._G
        cRG_save = Gs["cRG"]
        TRG_b = np.array([float(r[cols["T_RG"]]) for r in rows])
        _, _, cRG_b, _ = NC.gauge_constants(fc, header, TRG_b)
        Gs["cRG"] = cRG_b
        _, ((ngb, rgb),), _ = NC.record_contributions(recb, hb, (Gs["A"],), "fp32", rb)
        Gs["cRG"] = cRG_save
        ent = {"frame": kf, "base_row": rb, "frame_et": hf["et"], "row_et": hb["et"], "dt_s": hf["et"] - hb["et"],
               "gauge_position_diff_m": float(np.linalg.norm(np.asarray(hf["sc"]) - np.asarray(hb["sc"]))),
               "records_identical": bool(same), "records": [len(recf), len(recb)],
               "c9_NG_Pa_frame": R[kf]["ng_Pa"], "c9_NG_Pa_base_row": float(ngb.sum()),
               "c9_RG_Pa_frame": R[kf]["rg_Pa"], "c9_RG_Pa_base_row": float(rgb.sum()),
               "density_n_frame_m3": R[kf]["n"], "flux_phi_frame_m2s": R[kf]["phi"]}
        ent["rel_NG"] = ent["c9_NG_Pa_frame"] / ent["c9_NG_Pa_base_row"] - 1
        ent["rel_RG"] = ent["c9_RG_Pa_frame"] / ent["c9_RG_Pa_base_row"] - 1
        ip = int(float(rows[rb][cols["iPoint"]]))
        if ip in sb_ip:
            j = sb_ip[ip]
            ent["contribution_data_model_surface_NG_Pa"] = float(SB["model_surface_NG_Pa"][j])
            ent["contribution_data_model_surface_RG_Pa"] = float(SB["model_surface_RG_Pa"][j])
            ent["contribution_data_kind"] = str(SB["kind"][j])
        v3.append(ent)
    ver["3_frames_at_base_row_times"] = v3
    if a.frames_base_tail:                     # the overlap frame in both builds
        ko = tail_first - 1
        _, h1, r1 = NM.read_kernel_full(fb / "kernels" / ("point=%04d.bin" % ko))
        _, h2, r2 = NM.read_kernel_full(Path(a.frames_base_tail) / "kernels" / ("point=%04d.bin" % ko))
        same = (len(r1) == len(r2)) and all(np.array_equal(r1[n], r2[n]) for n in r1.dtype.names)
        row1 = [r for r in NM.base_table(fb)[3] if int(float(r[cols_f["iPoint"]])) == ko][0]
        ver["4_tail_overlap_frame"] = {"frame": ko, "records_identical": bool(same), "headers_identical": h1 == h2,
                                       "T_RG_first_build": row1[cols_f["T_RG"]], "T_RG_tail_build": overlap_row_t[cols_f["T_RG"]],
                                       "n_H2O": [row1[cols_f["n_H2O"]], overlap_row_t[cols_f["n_H2O"]]],
                                       "PhiRG_H2O": [row1[cols_f["PhiRG_H2O"]], overlap_row_t[cols_f["PhiRG_H2O"]]],
                                       "pRG": [row1[cols_f["pRG"]], overlap_row_t[cols_f["pRG"]]]}
    # colour range
    rng = {}
    for g, key in (("NG", "vng"), ("RG", "vrg")):
        z0 = np.load(out / "full_f00000.npz")
        v, aa = z0[key], z0["area"]
        o = np.argsort(-v)
        cs = np.cumsum((v * aa)[o])
        j = int(np.searchsorted(cs, a.bottom_share * cs[-1]))
        bottom = float(v[o][min(j, len(o) - 1)])
        pool = np.concatenate([np.load(out / ("full_f%05d.npz" % k))[key] for k in range(nfr - a.top_frames, nfr)])
        top = float(np.percentile(pool[pool > 0], a.top_percentile))
        rng[g] = {"bottom": bottom, "top": top, "decades": float(np.log10(top / bottom)),
                  "units": "m^-3 per m^2 (m^-5)" if g == "NG" else "m^-2 s^-1 per m^2 (m^-4 s^-1)",
                  "label": "NG: density per m² of surface [m⁻³ m⁻²]" if g == "NG" else "RG: flux into RG per m² [m⁻² s⁻¹ m⁻²]",
                  "bottom_rule": "per-area level above which the facets carry %g %% of the gauge's surface signal at frame 0"
                                 % (100 * a.bottom_share),
                  "top_rule": "%g th percentile of the per-area values of the facets with a contribution in frames %d-%d"
                              % (a.top_percentile, nfr - a.top_frames, nfr - 1)}
    json.dump(rng, open(out / "color_range.json", "x"), indent=1)
    with open(out / "frame_totals.tsv", "x") as f:
        f.write("frame\tutc\taltitude_m\trecords\tn_m3\tn_H2O\tn_CO2\tphi_m2s\tphi_H2O\tphi_CO2\tNG_surface_Pa\tRG_surface_Pa\n")
        for k in range(nfr):
            r = R[k]
            f.write("%d\t%s\t%s\t%d\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\n" % (
                k, rows_f[k][cols_f["UTC"]], rows_f[k][cols_f["altitude"]], r["records"], r["n"], r["nH"], r["nC"],
                r["phi"], r["pH"], r["pC"], r["ng_Pa"], r["rg_Pa"]))
    json.dump({"definition": __doc__, "inputs": vars(a), "c9_T": T, "c9_kappa": float(S["kappa"]), "base_T": T_base,
               "gauge_constants": gconst, "verification": ver, "color_range": rng, "seconds": time.time() - t0},
              open(out / "verification.json", "x"), indent=1, default=float)
    print(json.dumps({"verification": ver, "color_range": rng}, indent=1, default=float))
    print("frame_contrib done in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

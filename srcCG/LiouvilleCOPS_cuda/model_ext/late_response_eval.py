#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/late_response_eval.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/late_response_eval.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 58158de1b4c1f0785aa9e9d8b6fb5f01b96ac490a9ed6568b44102b1ac0279f5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Forward evaluation (no fit) of finished grain fits and their fine tunings on late-window bases built
#   on a fine time grid (late_time_grid.py + native/liouville_cops_cuda, optional GaugeOffsetSC), with the COPS time
#   response: inst (model at the time stamp), W4 (L2 interval), W10 (10-s running average), W10x5 (10 s + DPU
#   5-readout average); ROSINA COPS instrument modes RO-ROS-MAN-1019 sect. 4.2. Validates the zero-offset base against
#   each fit's own values. User 2026-11-07.
# USED BY: jobs/late_response.sbatch (step 3); figures/late_response_figures.py reads its output.
# RUN: .venv/bin/python late_response_eval.py --points late_grid_points.txt --base dz0=<base> [--base dz1.5=<base> ...]
#   --model NAME=<fit dir>,<fine-tuning dir> [...] --out <new dir>
# ==== end of header ====
"""Forward evaluation of finished 12M fits on late-window bases built on a fine time grid, with the COPS time response
(user 2026-11-07: "1: both" = 10-s running average and 10 s + DPU 5-readout running average; "3: make the sensitivity
check" = gauge offsets along +Z of the spacecraft; "for the forward run use two 12M best fits with and without extra
weight of the last point observations"). No fitting: every model keeps its fitted per-facet activity, temperature law,
CO2 scale kappa and dust-grain term; only the geometry (time grid, gauge offset) and the measurement model change.

    .venv/bin/python late_response_eval.py --points <late_grid_points.txt> \
        --base dz0=<base dir> [--base dz1.5=<dir> ...] \
        --model NAME=<fit dir>,<fine-tuning dir> [...] --out <dir>

Pressures (the formulas of tools/fit_cops.py Model.evaluate, SZAMode continuous, flat surface BC; checked on the fitting
base to 1e-6 relative): per kernel record of a point, night = (cos SZA <= 0 or cast shadow);
  q_H2O = 600 r_h^-4.2143229 (0.2 + 0.8 cos SZA | 0.2 night),  q_CO2 = 600 r_h^-2 (0.1 + 0.9 cos SZA | 0.1 night);
  T: analytical law of the fit (Tmax, T75, Tmin), u = T^-1/2;
  pNG = kT_NG [ (sqrt(m/2k) sqrt(pi)/2) sum q A Kn u ]_H2O / 0.893 + kappa [ ... A_CO2 ... ]_CO2 / 0.704  + grain(t)
  pRG = sqrt(2 pi m k T_RG) [ sum q A (KPhi + VdotB sqrt(m/2k) sqrt(pi)/2 u KV) ] / S, the same species sum;
  A = the fit's clipped H2O activity per facet (fine_tuning_state.npz activity_frozen; equal to the fit's
  activity_facets.dat) or the fine-tuned one (activity_final); A_CO2 = the base's elements.bin (run activity);
  grain(t) = max(0, A_g - B_g (t - t_end)) for t >= t_end (fit: grain_term.json result; fine-tuned: grain_finetune.json).
Measurement models for an observation whose L2 interval ends at T (the HK table's STOP_TIME; the observation table time
is that end, rounded; ROSINA COPS instrument modes RO-ROS-MAN-1019 sect. 4.2: readout every 2 s, electrometer time
constant ~10 s = "a running average over the last 10 s"; in monitoring modes the DPU averages the last 5 readouts):
  inst   : the model at the observation time (what every fit so far used);
  W4     : mean over the L2 interval [T-4, T] (reference only);
  W10    : mean over [T-10, T];
  W10x5  : mean of W10 ending at T, T-2, T-4, T-6, T-8 (a trapezoid over [T-18, T]).
Window means use the grid samples (linear in time, 0.05-s sub-samples, trapezoidal rule).
Outputs in --out: validation.txt (the dz0 base at the original observation times against each fit's own model values and
each fine tuning's prediction), late_response_<base>.dat (every observation point, every model and window),
series_<base>.npz (the 1-s model series), summary.txt (rms of ln(model/obs) over the late window and the last points).
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS/tools")
import fit_cops as fc                        # noqa: E402  (constants and the kernel reader of the fitter)
from fit_cops import read_kernel             # noqa: E402

KT = 1.0 / (1.0 - 1.0 / np.cos(np.radians(75.0)))
COS75 = np.cos(np.radians(75.0))
M_H, M_C = fc.MASS["H2O"], fc.MASS["CO2"]
BM_H, BM_C = (np.sqrt(m / (2.0 * fc.KB)) * fc.SQRTPI / 2.0 for m in (M_H, M_C))
S_H, S_C, T_NG = 0.893, 0.704, 293.0
WINDOWS = ("inst", "W4", "W10", "W10x5")


def fit_parameters(fit):
    vals = {}
    for line in open(fit / "fit_result.txt"):
        f = line.split()
        if len(f) == 2 and f[0] in ("Tmax", "T75", "Tmin", "kappa"):
            vals[f[0]] = float(f[1])
    return vals


def load_model(spec):
    name, rest = spec.split("=", 1)
    fit, ft = (Path(p) for p in rest.split(","))
    z = np.load(ft / "fine_tuning_state.npz", allow_pickle=True)
    meta = json.loads(str(z["meta_json"]))
    if Path(meta["source_restart"]["path"]).resolve() != (fit / "fit_restart.npz").resolve():
        raise SystemExit("%s: fine tuning %s was not made from %s" % (name, ft, fit))
    par = fit_parameters(fit)
    g_fit = json.load(open(fit / "grain_term.json"))
    g_ft = json.load(open(ft / "grain_finetune.json"))
    pred = {int(p): z["prediction"][i] for i, p in enumerate(z["point_ids"])}
    common = dict(Tmax=par["Tmax"], T75=par["T75"], Tmin=par["Tmin"], kappa=par["kappa"], t_end=float(g_fit["t_end_et"]))
    return [dict(common, name=name + ":fit", act=np.asarray(z["activity_frozen"]), A=float(g_fit["result"]["A_Pa"]),
                 B=float(g_fit["result"]["B_Pa_per_s"]), ref=fit / "pressure_residuals.dat", pred=None),
            dict(common, name=name + ":finetuned", act=np.asarray(z["activity_final"]), A=float(g_ft["A_Pa"]),
                 B=float(g_ft["B_Pa_per_s"]), ref=None, pred=pred)]


def read_model_table(path):
    cols, rows = None, []
    for line in open(path):
        if line.startswith("#iPoint"):
            cols = line[1:].split()
        elif not line.startswith("#"):
            rows.append(line.split())
    c = {k: i for i, k in enumerate(cols)}
    get = lambda k, f=float: np.array([f(r[c[k]]) for r in rows])
    return dict(ip=get("iPoint", int), utc=[r[c["UTC"]] for r in rows], et=get("et"), alt=get("altitude"),
                TRG=get("T_RG"), refNG=get("refNG"), refRG=get("refRG"), inside=get("inside", int))


def evaluate_base(base, models):
    """pNG, pRG of every model at every point of the base."""
    tab = read_model_table(base / "COPS_model.dat")
    el = np.memmap(base / "elements.bin", dtype="<f8", mode="r", offset=12).reshape(-1, 34)
    n = len(tab["ip"])
    out = {m["name"]: (np.zeros(n), np.zeros(n)) for m in models}
    for j, ip in enumerate(tab["ip"]):
        _, h, rec = read_kernel(str(base / "kernels" / ("point=%04d.bin" % ip)))
        rhA = h["rh"] / fc.AU
        c = np.minimum(rec["cos"], 1.0)
        night = (rec["cos"] <= 0.0) | (rec["shadow"] == 1)
        qH = np.where(night, 0.2, 0.2 + 0.8 * c) * 600.0 / rhA ** 4.2143229
        qC = np.where(night, 0.1, 0.1 + 0.9 * c) * 600.0 / rhA ** 2
        cc = np.where(night, -1.0, np.clip(rec["cos"], -1.0, 1.0))
        lit = cc >= COS75
        xlit = np.where(lit, 1.0 / np.where(lit, cc, 1.0) - 1.0, 0.0)
        aC = np.asarray(el[rec["face"], 33])
        cRG_H = np.sqrt(2 * np.pi * M_H * fc.KB * tab["TRG"][j]) / S_H
        cRG_C = np.sqrt(2 * np.pi * M_C * fc.KB * tab["TRG"][j]) / S_C
        VB = h["VdotB"]
        for m in models:
            Tl = m["Tmax"] + KT * (m["Tmax"] - m["T75"]) * xlit
            u = 1.0 / np.sqrt(np.where(lit & (Tl > m["Tmin"]), Tl, m["Tmin"]))
            a = m["act"][rec["face"]]
            nH, nC = (qH * a * BM_H * u * rec["Kn"]).sum(), (qC * aC * BM_C * u * rec["Kn"]).sum()
            fH = (qH * a * (rec["KPhi"] + VB * BM_H * u * rec["KV"])).sum()
            fC = (qC * aC * (rec["KPhi"] + VB * BM_C * u * rec["KV"])).sum()
            et = tab["et"][j]
            grain = max(0.0, m["A"] - m["B"] * (et - m["t_end"])) if et >= m["t_end"] else 0.0
            out[m["name"]][0][j] = fc.KB * T_NG * (nH / S_H + m["kappa"] * nC / S_C) + grain
            out[m["name"]][1][j] = cRG_H * fH + m["kappa"] * cRG_C * fC
        if j % 50 == 0:
            print("  %s: point %d/%d (%s, %d records)" % (base.name, j, n, tab["utc"][j], len(rec)), flush=True)
    return tab, out


def window_mean(tg, yg, t_end, width, dt=0.05):
    s = np.arange(t_end - width, t_end + 0.5 * dt, dt)
    s[-1] = t_end
    if s[0] < tg[0] - 1e-6 or s[-1] > tg[-1] + 1e-6:
        return np.nan
    return np.trapz(np.interp(s, tg, yg), s) / (s[-1] - s[0])


def responses(tg, yg, T):
    w10 = lambda e: window_mean(tg, yg, e, 10.0)
    return {"W4": window_mean(tg, yg, T, 4.0), "W10": w10(T), "W10x5": np.mean([w10(T - 2.0 * k) for k in range(5)])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", required=True, type=Path)
    ap.add_argument("--base", action="append", required=True)
    ap.add_argument("--model", action="append", required=True)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    pts = [l.split() for l in open(a.points) if not l.startswith("#")]     # index UTC kind orig L2_start L2_stop
    kind = np.array([p[2] for p in pts])
    orig = np.array([int(p[3]) for p in pts])
    models = [v for spec in a.model for v in load_model(spec)]
    names = [m["name"] for m in models]
    summary, validation = [], []
    for spec in a.base:
        label, path = spec.split("=", 1)
        base = Path(path)
        tab, val = evaluate_base(base, models)
        if len(tab["ip"]) != len(pts) or any(tab["utc"][i] != pts[i][1] for i in range(len(pts))):
            raise SystemExit("%s: points differ from %s" % (base, a.points))
        grid = np.nonzero(kind == "grid")[0]
        obs = np.nonzero(kind == "obs")[0]
        tg = tab["et"][grid]
        # L2 interval end of each observation as ephemeris time: its own et plus (STOP - observation time) in UTC seconds
        stop = np.array([tab["et"][i] + (utc_seconds(pts[i][5]) - utc_seconds(pts[i][1])) for i in obs])
        res = {}
        for nm in names:
            pNG, pRG = val[nm]
            r = {"inst": (pNG[obs], pRG[obs])}
            for w in ("W4", "W10", "W10x5"):
                r[w] = tuple(np.array([responses(tg, y[grid], T)[w] for T in stop]) for y in (pNG, pRG))
            res[nm] = r
        refNG, refRG = tab["refNG"][obs], tab["refRG"][obs]
        # validation: the base with the fitting geometry (dz0) at the original observation times
        if label == "dz0":
            for m in models:
                pNG, pRG = res[m["name"]]["inst"]
                if m["ref"] is not None:
                    ref = {}
                    for line in open(m["ref"]):
                        if not line.startswith("#"):
                            f = line.split()
                            ref[(int(f[0]), f[3])] = float(f[7])
                    eNG = max(abs(pNG[k] / ref[(orig[i], "NG")] - 1) for k, i in enumerate(obs))
                    eRG = max(abs(pRG[k] / ref[(orig[i], "RG")] - 1) for k, i in enumerate(obs))
                    src = "the fit's pressure_residuals.dat"
                else:
                    eNG = max(abs(pNG[k] / m["pred"][orig[i]][0] - 1) for k, i in enumerate(obs))
                    eRG = max(abs(pRG[k] / m["pred"][orig[i]][1] - 1) for k, i in enumerate(obs))
                    src = "the fine tuning's prediction"
                validation.append("%-28s max |rel diff| vs %s: NG %.2e, RG %.2e" % (m["name"], src, eNG, eRG))
                if max(eNG, eRG) > 1e-4:
                    raise SystemExit("validation failed: " + validation[-1])
        with open(a.out / ("late_response_%s.dat" % label), "w") as f:
            f.write("# %s: %s; gauge altitude = altitude column of the base; obs = observation (Pa)\n" % (label, base))
            f.write("# UTC L2_stop altitude_m obsNG obsRG " + " ".join("%s:%s:%s" % (nm, g, w) for nm in names for g in ("NG", "RG") for w in WINDOWS) + "\n")
            for k, i in enumerate(obs):
                f.write("%s %s %.3f %.4e %.4e " % (pts[i][1], pts[i][5], tab["alt"][i], refNG[k], refRG[k]))
                f.write(" ".join("%.5e" % res[nm][w][gi][k] for nm in names for gi in (0, 1) for w in WINDOWS) + "\n")
        np.savez(a.out / ("series_%s.npz" % label), et=tab["et"], alt=tab["alt"], kind=kind, utc=np.array(tab["utc"]),
                 obs_index=obs, stop_et=stop, refNG=refNG, refRG=refRG,
                 **{"%s|%s" % (nm, g): val[nm][gi] for nm in names for gi, g in enumerate(("NG", "RG"))},
                 **{"%s|%s|%s" % (nm, g, w): res[nm][w][gi] for nm in names for gi, g in enumerate(("NG", "RG")) for w in WINDOWS})
        last5 = np.arange(len(obs))[-5:]
        for nm in names:
            for w in WINDOWS:
                lnNG = np.log(res[nm][w][0] / refNG)
                lnRG = np.log(res[nm][w][1] / refRG)
                summary.append("%-6s %-30s %-6s %7.4f %7.4f | %+7.3f %+7.3f | %s | %s" % (
                    label, nm, w, np.sqrt(np.nanmean(lnNG ** 2)), np.sqrt(np.nanmean(lnRG ** 2)), lnNG[-1], lnRG[-1],
                    " ".join("%+6.3f" % v for v in lnNG[last5]), " ".join("%+6.3f" % v for v in lnRG[last5])))
        print("%s done" % label, flush=True)
    with open(a.out / "validation.txt", "w") as f:
        f.write("\n".join(validation) + "\n")
    with open(a.out / "summary.txt", "w") as f:
        f.write("# late window %s .. %s (%d observations); ln(model/obs)\n" % (pts[obs[0]][1], pts[obs[-1]][1], len(obs)))
        f.write("# base   model                          window rms_NG  rms_RG |  last NG  last RG | NG last 5 | RG last 5\n")
        f.write("\n".join(summary) + "\n")
    print(open(a.out / "validation.txt").read())
    print(open(a.out / "summary.txt").read())


def utc_seconds(s):
    d, t = s.split("T")
    h, m, sec = t.split(":")
    return int(h) * 3600 + int(m) * 60 + float(sec)

if __name__ == "__main__":
    main()

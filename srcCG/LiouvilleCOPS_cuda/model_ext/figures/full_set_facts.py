#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/full_set_facts.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/full_set_facts.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 38453c8a2e17f927ce05fb667d2a0c56e8e78a5526390f548e30000f28ea77d2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Facts of one full figure set computed from its data (2026-11-17)
# RUN: the full sets of the stage-2 states (jobs/full_set_submit.sh, jobs/finish_full_set.sh)
# ==== end of header ====
"""Facts of ONE full figure set, computed from its data (2026-11-17, full sets of the stage-2 states; the "Facts" section
of the c9 set's README.txt was computed by hand from the same products). CPU (reads the final A of the state: srun).

    full_set_facts.py --set <full_set folder> --ft <FT dir> --fit <FIT dir> --maps <maps> --contrib <contributions>
                      --phase <contrib_phase data> --base <fit base> [--chain <chain root> --cycles N] [--caveat TEXT ...]
Writes <set>/checks/facts.json and prints the Facts text (README.txt is written by the caller).
Numbers: data cost after each cycle (fit | fine tuning) and the last terminations (chains); multipliers and their bounds;
kappa (value, at its lower bound 0); A = 0 on facets / area / seen area; Q_H2O, Q_CO2 at 10:33:10 and 10:39:26 UTC (maps
manifest); production share by slope bin (contributions slope_bins.npz); the modelled surface-pressure share by slope
bin before the NG gap, in the late window, below 20 m and at the last observation (slope_bins.npz fractions); P50 / P90
source-area unions (checks/source_area_tables.txt); lobe and distance-shell shares (contrib_phase regions.npz).
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stage2_figures as S2  # noqa: E402  (kv, fit_result)
import pub_maps_figures as M  # noqa: E402  (multiplier_bounds, bin_edges_from_meta)


def med(x):
    return float(np.median(x)) if len(x) else float("nan")


def main():
    ap = argparse.ArgumentParser()
    for k in ("--set", "--ft", "--fit", "--maps", "--contrib", "--phase", "--base"):
        ap.add_argument(k, required=True)
    ap.add_argument("--chain", default=None)
    ap.add_argument("--cycles", type=int, default=5)
    ap.add_argument("--caveat", action="append", default=[])
    a = ap.parse_args()
    O, ft, fit = Path(a.set), Path(a.ft), Path(a.fit)
    F = {}
    L = []
    # chain history
    if a.chain:
        ch = Path(a.chain)
        hist = []
        for j in range(1, a.cycles + 1):
            fr = ch / ("c%d" % j) / "fit" / "final" / "fit_result.txt"
            fn = ch / ("c%d" % j) / "ft" / "final" / "finetune_result.txt"
            f1 = S2.fit_result(fr) if fr.exists() else {}
            k1 = S2.kv(fn) if fn.exists() else {}
            hist.append({"cycle": j, "fit_data": f1.get("final_data_cost"), "fit_stop": f1.get("termination"),
                         "fit_nfev": f1.get("nfev"), "ft_data": float(k1["cost_data"]) if "cost_data" in k1 else None,
                         "ft_stop": k1.get("termination"), "ft_nfev": k1.get("nfev"),
                         "ft_delta_penalty": float(k1["cost_delta_penalty"]) if "cost_delta_penalty" in k1 else None})
        F["history"] = hist
        L.append("* Cycles: data cost after each fit | fine tuning: " + "; ".join(
            "c%d %.4e | %.4e" % (h["cycle"], h["fit_data"], h["ft_data"]) for h in hist) +
            ". Last fit %s after %s evaluations, last fine tuning %s after %s (P_delta %s). Reaching max_nfev is not "
            "convergence." % (S2.short(hist[-1]["fit_stop"]), hist[-1]["fit_nfev"], hist[-1]["ft_stop"], hist[-1]["ft_nfev"],
                              ("%.2e" % hist[-1]["ft_delta_penalty"]) if hist[-1]["ft_delta_penalty"] is not None else "-"))
    else:
        fr = S2.fit_result(Path(json.load(open(ft / "fine_tuning.json")).get("d2", str(fit))) / "fit_result.txt")
        F["fit"] = fr
        L.append("* Final state = the fit (no fine tuning): data cost %.4e, %s after %s evaluations (not converged)." % (
            fr["final_data_cost"], S2.short(fr["termination"]), fr["nfev"]))
    # multipliers, kappa
    mb = M.multiplier_bounds(fit)
    man = json.load(open(Path(a.maps) / "maps_manifest.json"))
    sc = man["state_checks"]
    kappa = float(sc["kappa"])
    F["multipliers"] = [(n, float(x), float(lo)) for n, x, lo in mb]
    F["kappa"] = kappa
    if mb:
        L.append("* Multipliers (0-40 deg held at 1): " + ", ".join(
            "%s %.4g%s" % (n.replace("c_slope_", "").replace("_", "-"), x, " (at its bound %g)" % lo
                           if x - lo <= 1e-9 * max(1.0, abs(lo)) else "") for n, x, lo in mb) + ".")
    else:
        L.append("* No slope multipliers (all bins 1).")
    L.append("* kappa (CO2 scale on the srcCG CO2 map; lower bound 0) = %.3e: at its lower bound (the fits prefer the H2O "
             "law for the whole signal; total pressure only, so this does not measure CO2 as absent)." % kappa)
    # A = 0
    with np.load(ft / "finetune_state.npz", allow_pickle=False) as z:
        A = np.asarray(z["A_final"], np.float64)
    with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
        seen = np.asarray(z["seen"]).astype(bool)
    with open(Path(a.base) / "elements.bin", "rb") as f:
        f.read(12)
        d = np.fromfile(f, dtype="<f8").reshape(-1, 34)
    area = np.ascontiguousarray(d[:, 6])
    del d
    z0 = A <= 0.0
    F["A_zero"] = {"facets": float(z0.mean()), "area": float(area[z0].sum() / area.sum()),
                   "seen_area": float(area[z0 & seen].sum() / area[seen].sum()), "nonzero_seen_facets": int((~z0 & seen).sum())}
    L.append("* A = 0 on %.1f %% of the facets, %.1f %% of the area, %.1f %% of the seen area; %d seen facets have A > 0." % (
        100 * F["A_zero"]["facets"], 100 * F["A_zero"]["area"], 100 * F["A_zero"]["seen_area"], F["A_zero"]["nonzero_seen_facets"]))
    # Q
    Q = {t["utc"][-8:]: (t["Q_H2O_s-1"], t["Q_CO2_s-1"]) for t in man["times"]}
    F["Q"] = Q
    # seen / unseen share of the H2O production (the maps' flux files, as the model computes it)
    nmc = M.load_module(Path(man["code_dir"]) / "nonneg_maps.py", "nonneg_maps_of_maps_facts")
    with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
        fbins = np.asarray(z["bin"]).astype(np.intp)
    share = {}
    for t in man["times"]:
        lab = t["utc"][-8:]
        if lab not in ("10:33:10", "10:39:26"):
            continue
        rd, _, _ = nmc.read_vtu(Path(a.maps) / t["vtu"])
        P = np.asarray(rd["F_H2O_per_m2_s"], np.float64) * area
        tot = P.sum()
        share[lab] = {"Q_H2O_from_vtu": float(tot), "Q_H2O_manifest": float(t["Q_H2O_s-1"]),
                      "seen_share": float(P[seen].sum() / tot), "unseen_share": float(P[~seen].sum() / tot),
                      "unseen_ge60deg_share": float(P[(~seen) & (fbins >= 3)].sum() / tot),
                      "seen_ge60deg_share": float(P[seen & (fbins >= 3)].sum() / tot)}
        del rd, P
    F["production_seen_unseen"] = share
    L.append("* Q_H2O = %.3g s^-1 at 10:33:10 UTC, %.3g at 10:39:26 (Q_CO2 %.2g / %.2g); at 10:33:10 %.1f %% of Q_H2O comes "
             "from facets no fitted observation sees (%.1f %% from unseen facets >= 60 deg), %.1f %% from the seen ones." % (
                 Q["10:33:10"][0], Q["10:39:26"][0], Q["10:33:10"][1], Q["10:39:26"][1],
                 100 * share["10:33:10"]["unseen_share"], 100 * share["10:33:10"]["unseen_ge60deg_share"],
                 100 * share["10:33:10"]["seen_share"]))
    # slope bins
    with np.load(Path(a.contrib) / "slope_bins.npz", allow_pickle=False) as z:
        sb = {k: np.asarray(z[k]) for k in z.files}
    labs = [str(x) for x in sb["bin_labels_deg"]]
    steep = np.array([k >= 3 for k in range(len(labs))])
    F["production_share_seen_103310"] = dict(zip(labs, [float(v) for v in sb["production_share_seen_103310"]]))
    F["seen_area_share"] = dict(zip(labs, [float(v) for v in sb["bin_share_of_seen_area"]]))
    L.append("* Slope bins >= 60 deg: %.1f %% of the seen area, %.1f %% of the seen H2O production at 10:33:10 (%s)." % (
        100 * sb["bin_share_of_seen_area"][steep].sum(), 100 * sb["production_share_seen_103310"][steep].sum(),
        ", ".join("%s %.1f %%" % (l, 100 * v) for l, v in zip(labs, sb["production_share_seen_103310"]))))
    kind = np.array([str(k) for k in sb["kind"]])
    utc = np.array([str(u) for u in sb["utc"]])
    with np.load(Path(a.phase) / "regions.npz", allow_pickle=False) as z:
        R = {k: np.asarray(z[k]) for k in z.files}
    alt = R["alt"]
    if not np.array_equal(R["ipoint"], sb["iPoint"]):
        raise SystemExit("regions and slope-bin rows differ")
    pre = (kind == "obs") & np.array([u.split("T")[1] < "10:33:10" for u in utc])
    late = kind == "late_obs"
    low = late & (alt < 20.0)
    last = np.nonzero(late)[0][np.argmin(alt[late])]
    fr = {g: sb["fraction_%s" % g] for g in ("NG", "RG")}
    use = {"NG": sb["use_ng"] > 0, "RG": sb["use_rg"] > 0}
    sh = {}
    for g in ("NG", "RG"):
        sh[g] = {"pre_mean": [float(v) for v in fr[g][pre & use[g]].mean(axis=0)] if (pre & use[g]).any() else None,
                 "late_mean": [float(v) for v in fr[g][late].mean(axis=0)],
                 "below20_mean": [float(v) for v in fr[g][low].mean(axis=0)] if low.any() else None,
                 "last": [float(v) for v in fr[g][last]], "last_alt_m": float(alt[last])}
    F["pressure_share_by_bin"] = sh
    L.append("* Share of the modelled surface pressure by slope bin (contributions, means over rows): before the NG gap NG "
             ">= 60 deg %.0f %%, 0-40 %.0f %%; late window >= 60 deg NG %.0f %%, RG %.0f %% (90-180 deg %.0f %% / %.0f %%); "
             "below 20 m (%d rows) 90-180 deg NG %.0f %%, RG %.0f %%; at the last observation (%.1f m) 90-180 deg NG %.0f %%, "
             "RG %.0f %%." % (
                 100 * sum(sh["NG"]["pre_mean"][3:]), 100 * sh["NG"]["pre_mean"][0], 100 * sum(sh["NG"]["late_mean"][3:]),
                 100 * sum(sh["RG"]["late_mean"][3:]), 100 * sh["NG"]["late_mean"][6], 100 * sh["RG"]["late_mean"][6],
                 int(low.sum()), 100 * sh["NG"]["below20_mean"][6], 100 * sh["RG"]["below20_mean"][6], alt[last],
                 100 * sh["NG"]["last"][6], 100 * sh["RG"]["last"][6]))
    # source areas
    t = (O / "checks" / "source_area_tables.txt").read_text()
    sa = {}
    for g in ("NG", "RG"):
        for p in (50, 90):
            m = re.search(r"\[%s, P = %d %%\] facets with N > 0: (\d+) .*?area (\S+) m\^2" % (g, p), t)
            sa["%s_P%d" % (g, p)] = {"facets": int(m.group(1)), "area_m2": float(m.group(2))}
    F["source_area_unions"] = sa
    L.append("* P50 source-area unions over the 95 post-gap observations: NG %d facets, %.3g m^2; RG %d facets, %.3g m^2 "
             "(P90: NG %.3g m^2, RG %.3g m^2)." % (sa["NG_P50"]["facets"], sa["NG_P50"]["area_m2"], sa["RG_P50"]["facets"],
                                                 sa["RG_P50"]["area_m2"], sa["NG_P90"]["area_m2"], sa["RG_P90"]["area_m2"]))
    # regions
    lob = [str(x) for x in R["lobe_names"]]
    ibig = [k for k, n in enumerate(lob) if "big" in n.lower() or "body" in n.lower()]
    ismall = [k for k, n in enumerate(lob) if k not in ibig]
    with np.errstate(divide="ignore", invalid="ignore"):
        big_ng = R["lobe_ng_Pa"][:, ibig[0]] / R["total_ng_Pa"] if ibig else None
        sm_ng = R["lobe_ng_Pa"][:, ismall[0]] / R["total_ng_Pa"] if ismall else None
        sm_rg = R["lobe_rg_Pa"][:, ismall[0]] / R["total_rg_Pa"] if ismall else None
        s30 = (R["shell_ng_Pa"][:, 1] / R["total_ng_Pa"], R["shell_rg_Pa"][:, 1] / R["total_rg_Pa"])
    pre_r = R["pre"] & R["obs"]
    F["regions"] = {"lobe_names": lob, "shells_m": [float(v) for v in R["shells_m"]],
                    "big_lobe_NG_median_pre": med(big_ng[pre_r]) if big_ng is not None else None,
                    "small_lobe_median_late": [med(sm_ng[R["late"]]), med(sm_rg[R["late"]])] if sm_ng is not None else None,
                    "shell_30_100m_median_late": [med(s30[0][R["late"]]), med(s30[1][R["late"]])]}
    L.append("* Regions (contrib_phase): NG before the late window from the %s median %.0f %%; late window from the %s "
             "median NG %.0f %%, RG %.0f %%; from 30-100 m around the touchdown point median NG %.0f %%, RG %.0f %%." % (
                 lob[ibig[0]] if ibig else "-", 100 * F["regions"]["big_lobe_NG_median_pre"],
                 lob[ismall[0]] if ismall else "-", 100 * F["regions"]["small_lobe_median_late"][0],
                 100 * F["regions"]["small_lobe_median_late"][1], 100 * F["regions"]["shell_30_100m_median_late"][0],
                 100 * F["regions"]["shell_30_100m_median_late"][1]))
    for c in a.caveat:
        L.append("* Caveat: " + c)
    F["text"] = L
    json.dump(F, open(O / "checks" / "facts.json", "x"), indent=1, default=str)
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())

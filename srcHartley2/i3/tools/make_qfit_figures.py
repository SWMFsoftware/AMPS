#!/usr/bin/env python3
"""Comparison figures of the OH-fitted SA1982_VH campaign (runs NOT converged): OH (0-0),
H Lyman-alpha and O I 130.4 nm of all six cases against MAVEN/IUVS, and case 5 O I 130.4 nm
against the digitized MAVEN/IUVS profile.

Copy of validation_runs/i3_qfit_campaign_20261123/report_20261124/scripts/make_qfit_figures.py
(2026-11-25, user request), made to run from git: the helpers come from i3_figure_helpers.py and the
small input files from inputs/ (both copied here, see inputs/README.md); repository-relative defaults and
the options --campaign, --report-dir and --predicted. Paths below are relative to the AMPS root.

Runs
- <campaign>/sa1982_vh/caseN_np128/production, default campaign
  validation_runs/i3_qfit_campaign_20261123: inner-boundary fix + SA1982_VH OH loss, Q(H2O) refitted
  to MAVEN/IUVS OH from output 4 of the cancelled runs (README.md of the campaign); the newest output
  written by all six runs unless --output is given. Shown as run, NOT scaled.
- "predicted": the values that the fit expected, read from inputs/ohfit_numbers_output4.json (copy of
  validation_runs/i3_fixed_campaign_20261122/report_interim_20261123/data/ohfit_numbers.json:
  cancelled SA1982_VH runs at output 4 x fit factor; "fill-corr." = their fill-corrected estimate).
  Numbers only. The script also checks that each run's input Q(H2O) is the fitted Q of that file.
The default-chemistry runs are not shown (user request, 2026-11-24).

Brightness is edge-corrected (srcHartley2/visualization/comet_edge_correction.py, through
make_interim_figures.load) and shown to 2/3 of the domain radius. O I 130.4 nm = direct (KVB2022
solar line) + Ly-beta Bowen; the SUMER2001 direct model is not shown (user request). Model/obs for
OH and H uses the MAVEN/IUVS points inside the tabulated range; case 5 O uses edge-corrected means
over each IUVS bin, as in the earlier reports. No optical-depth band is drawn (user request); the
range where the O line-centre tau exceeds 0.1 is given in the panel text.

Inputs in git: i3_figure_helpers.py, inputs/ (case-5 IUVS O profile, the slide's model curve, the fit
numbers), srcHartley2/visualization (edge correction; MAVEN/IUVS OH and H points in
references/i3_brightness_MAVEN). Outside git: only the AMPS run output (--campaign).

Figures (300 dpi, <report-dir>/plots/): qfit_OH_H_O_vs_maven_out<K>.png,
qfit_o_case5_vs_iuvs_out<K>.png. Data: <report-dir>/data/qfit_numbers_out<K>.json,
<report-dir>/data/qfit_o_case5_vs_iuvs_out<K>.csv. Default report directory:
<campaign>/report_20261124. Run from any directory:
    python3 srcHartley2/i3/tools/make_qfit_figures.py [--output K] [--campaign DIR] [--report-dir DIR]
Reads run output only; writes only under <report-dir>/plots and <report-dir>/data.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent        # srcHartley2/i3/tools
ROOT = HERE.parents[2]                        # AMPS repository root
VR = ROOT / "validation_runs"
CAMPAIGN = VR / "i3_qfit_campaign_20261123"   # default; --campaign
REPORT = CAMPAIGN / "report_20261124"         # default; --report-dir
sys.path.insert(0, str(HERE))
import i3_figure_helpers as H  # noqa: E402  (loaders, edge correction, palette, case-5 observation)

np, plt, Line2D = H.np, H.plt, H.Line2D
P = H.P
S1, INK, MUTED = H.S1, H.INK, H.MUTED
KVB, BOWEN = H.KVB, H.BOWEN
OKEY = f"O total {KVB}"                       # direct (KVB2022) + Ly-beta Bowen, R
CASES = H.CASES
QFIT = {c: CAMPAIGN / f"sa1982_vh/case{c}_np128/production" for c in CASES}
PRED = HERE / "inputs/ohfit_numbers_output4.json"   # default; --predicted
X_LIM = (0.2, 25.0)                           # 10^5 km, all rows (as model_vs_obs_Q.png)


def rel(path):
    """Path relative to the AMPS root when inside it, else absolute (for the JSON)."""
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def corrected(d, key):
    return H.bright(d, key) * H.edge_of(d, key)


def geo(v):
    return float(np.exp(np.mean(np.log(v))))


def maven_ratio(c, sp, d):
    """Model/obs at the MAVEN/IUVS points inside the tabulated range (edge-corrected model)."""
    _, obs = P.read_reference(P.DEFAULT_REFERENCE_DIR, c, sp)
    y = corrected(d, sp)
    ok = np.isfinite(y) & (y > 0)
    inside = obs["x"] * 1e8 <= d["b"][ok].max()
    r = H.loglog_interp(obs["x"][inside] * 1e8, d["b"], y) / obs["y"][inside]
    return {"geomean": geo(r), "median": float(np.median(r)), "min": float(r.min()), "max": float(r.max()),
            "n_points": int(len(r)), "n_points_total": int(len(obs["x"]))}


def at_km(d, key, km):
    return float(H.loglog_interp(np.array([km * 1e3]), d["b"], corrected(d, key))[0])


def title(c, what, d):
    return (f"Case {c}: {what}\nr$_h$ = {d['rh']:.2f} AU, $\\dot r$ = {d['rdot']:+.1f} km/s, "
            f"{d['iterations'] // 1000}k it{'' if d['converged'] else ' (running)'}")


def main():
    global CAMPAIGN, REPORT, QFIT, PRED
    parser = argparse.ArgumentParser(description=" ".join(__doc__.split("\n\n")[0].split()))
    parser.add_argument("--output", type=int, default=None,
                        help="output index of the OH-fitted runs (default: newest written by all six)")
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN,
                        help="campaign directory with sa1982_vh/caseN_np128/production (default: %(default)s)")
    parser.add_argument("--report-dir", type=Path, default=None,
                        help="where plots/ and data/ are written (default: <campaign>/report_20261124)")
    parser.add_argument("--predicted", type=Path, default=PRED,
                        help="fit numbers of the cancelled runs (default: %(default)s)")
    args = parser.parse_args()
    PRED = args.predicted.resolve()
    CAMPAIGN = args.campaign.resolve()
    REPORT = (args.report_dir or CAMPAIGN / "report_20261124").resolve()
    QFIT = {c: CAMPAIGN / f"sa1982_vh/case{c}_np128/production" for c in CASES}
    out = args.output if args.output is not None else min(max(H.outputs_of(QFIT[c])) for c in CASES)
    (REPORT / "plots").mkdir(parents=True, exist_ok=True)
    (REPORT / "data").mkdir(exist_ok=True)
    pred = json.loads(PRED.read_text())

    runs = {}
    for c in CASES:
        if out not in H.outputs_of(QFIT[c]):
            raise SystemExit(f"case {c}: output {out} of the OH-fitted run is incomplete")
        n = H.load(QFIT[c], c, out)
        if n["variants"].get("OH_LOSS") != "SA1982_VH":
            raise SystemExit(f"case {c}: the run is not SA1982_VH ({n['variants']})")
        q_fit = pred["cases"][str(c)]["now"]["q_h2o_fitted_s-1"]
        if abs(n["q_h2o"] / q_fit - 1) > 1e-3:
            raise SystemExit(f"case {c}: input Q(H2O) {n['q_h2o']:.4e} is not the fitted {q_fit:.4e}")
        runs[c] = n
        print(f"case {c}: output {out} ({H.status(n)}), Q(H2O) {n['q_h2o']:.4g}")

    numbers = {"output": out, "window_iterations": [out * H.SAMPLE_LENGTH + 1, (out + 1) * H.SAMPLE_LENGTH],
               "o_model": "direct OI_1304_TRIPLET_DIRECT_KVB2022 + OI_1304_TRIPLET_LYB_BOWEN_SUMER2001, edge-corrected",
               "predicted_from": rel(PRED), "cases": {}}
    for c in CASES:
        n = runs[c]
        p = pred["cases"][str(c)]
        tau = n["cols"]["O"][f"Line Center Tau {KVB}"]
        spans = H.flagged_spans(n["b"], n["cols"]["O"][f"Tau Threshold Exceeded {KVB}"])
        e = numbers["cases"][str(c)] = {
            "run": rel(QFIT[c]), "output": out, "iterations": n["iterations"],
            "converged": n["converged"], "q_h2o_s-1": n["q_h2o"], "variants": n["variants"], "sha256": n["sha256"],
            "rh_au": n["rh"], "rdot_km_s": n["rdot"],
            "predicted": {"h_model_over_obs": p["now"]["h_model_over_obs_after_fit_geomean"],
                          "h_model_over_obs_fill_corrected": p["now_est"]["h_model_over_obs_after_fit_geomean"]},
            "O_KVB2022_bowen_edge_corrected_R": {f"{km:.0e} km": at_km(n, OKEY, km) for km in (1e3, 1e4, 1e5, 1e6)},
            "O_line_centre_tau_max": float(np.max(tau)),
            "O_tau_above_0.1_up_to_km": max((h for _, h in spans), default=0.0) / 1e3}
        for sp in ("OH", "H"):
            e[f"{sp}_model_over_obs"] = maven_ratio(c, sp, n)

    # ---------------------------------------------------------- case 5 bin means
    c5 = H.OBS_CASE
    n5 = runs[c5]
    b5 = n5["b"]
    obs = H.read_commented_csv(H.OBS_CSV)
    slide = H.read_commented_csv(H.SLIDE_CSV)
    xo = np.array([float(r["distance_1e5_km"]) for r in obs])
    lo = np.array([float(r["distance_lower_1e5_km"]) for r in obs])
    hi = np.array([float(r["distance_upper_1e5_km"]) for r in obs])
    yo = np.array([float(r["radiance_R"]) for r in obs])
    if hi.max() * 1e8 > 2.0 / 3.0 * n5["r_edge"]:
        raise SystemExit("an observed bin lies beyond 2/3 of the domain radius")
    mean5 = np.array([H.bin_average(b5, corrected(n5, OKEY), l * 1e8, h * 1e8) for l, h in zip(lo, hi)])
    r5 = yo / mean5
    pc5 = pred["cases"][str(c5)]["o_case5"]
    numbers["cases"][str(c5)]["o_case5_obs_over_model_binmean"] = {
        "per_bin": [float(v) for v in r5], "geomean": geo(r5),
        "predicted_geomean": pc5["now_KVB2022"]["geomean"],
        "predicted_geomean_fill_corrected": pc5["now_est_KVB2022"]["geomean"]}
    with open(REPORT / f"data/qfit_o_case5_vs_iuvs_out{out}.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["bin", "distance_1e5_km", "observed_R",
                                                    "model_KVB2022_bowen_R_binmean_edge", "obs_over_model"])
        writer.writeheader()
        for k in range(len(xo)):
            writer.writerow({"bin": obs[k]["bin"], "distance_1e5_km": f"{xo[k]:.4f}", "observed_R": f"{yo[k]:.3f}",
                             "model_KVB2022_bowen_R_binmean_edge": f"{mean5[k]:.4f}", "obs_over_model": f"{r5[k]:.4f}"})

    # ------------------------------------- figure 1: OH, H and O of all six cases
    fig, axes = plt.subplots(3, 6, figsize=(19, 11.4), dpi=300)
    for j, c in enumerate(CASES):
        n = runs[c]
        b = n["b"]
        e = numbers["cases"][str(c)]
        r_plot = 2.0 / 3.0 * n["r_edge"]
        for i, sp in enumerate(("OH", "H")):
            ax = axes[i, j]
            _, ob = P.read_reference(P.DEFAULT_REFERENCE_DIR, c, sp)
            y = corrected(n, sp)
            ok = np.isfinite(y) & (y > 0) & (b <= r_plot)
            ax.plot(b[ok] / 1e8, y[ok], color=S1, lw=1.6, zorder=3)
            win = ok & (b > 2e7)
            ylo, yhi = min(ob["y"].min(), y[win & (b < 2e9)].min()), max(ob["y"].max(), y[win].max())
            yerr = np.vstack([ob["y"] - np.where(np.isnan(ob["y_lower"]), ob["y"], ob["y_lower"]),
                              np.where(np.isnan(ob["y_upper"]), ob["y"], ob["y_upper"]) - ob["y"]])
            ax.errorbar(ob["x"], ob["y"], yerr=yerr, xerr=[ob["x"] - ob["x_lower"], ob["x_upper"] - ob["x"]],
                        fmt="o", ms=3.5, color=INK, ecolor=MUTED, elinewidth=0.6, capsize=0, zorder=5)
            r = e[f"{sp}_model_over_obs"]
            text = f"model/obs: {r['geomean']:.2f} (geo. mean),\n{r['min']:.2f}-{r['max']:.2f} over the points"
            if sp == "OH":
                text = f"Q(H$_2$O) = {H.sci(n['q_h2o'])} s$^{{-1}}$ (fitted to OH)\n" + text
            else:
                text += (f"\npredicted by the fit: {e['predicted']['h_model_over_obs']:.2f} "
                         f"(fill-corr. {e['predicted']['h_model_over_obs_fill_corrected']:.2f})")
            ax.text(0.97, 0.96, text, transform=ax.transAxes, ha="right", va="top", color=INK, fontsize=7.0,
                    linespacing=1.45)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(*X_LIM)
            ax.set_ylim(max(1e-3, 0.4 * ylo), 8.0 * yhi)
            H.style(ax)
            ax.set_title(title(c, "OH (0-0)" if sp == "OH" else "H Lyman-alpha", n), color=INK, loc="left",
                         fontsize=8.4, linespacing=1.35)
            if j == 0:
                ax.set_ylabel(f"{sp} brightness (kR)")
        # ---- O row: direct (KVB2022) + Ly-beta Bowen
        ax = axes[2, j]
        y = corrected(n, OKEY)
        ok = np.isfinite(y) & (y > 0) & (b <= r_plot)
        ax.plot(b[ok] / 1e8, y[ok], color=S1, lw=1.6, zorder=3)
        shown = y[ok & (b / 1e8 >= X_LIM[0])]
        ylo, yhi = shown.min(), shown.max()
        oc = e["O_KVB2022_bowen_edge_corrected_R"]
        tau_km = e["O_tau_above_0.1_up_to_km"]
        tau_text = (f"line-centre $\\tau_0$ > 0.1 for b $\\leq$ {tau_km / 1e5:.2g}$\\times$10$^5$ km"
                    if tau_km > 0 else "line-centre $\\tau_0$ $\\leq$ 0.1")
        text = f"{oc['1e+04 km']:.3g} R at 10$^4$ km, {oc['1e+05 km']:.3g} R at 10$^5$ km\n{tau_text}"
        if c == c5:
            ylo, yhi = min(ylo, yo.min()), max(yhi, yo.max())
            # the model bin means are not drawn here (user request, 2026-11-24); they are in the text and
            # in qfit_o_case5_vs_iuvs_out<K>.png
            ax.errorbar(xo, yo, xerr=[xo - np.maximum(lo, X_LIM[0]), hi - xo], fmt="o", ms=3.8, color=INK,
                        elinewidth=0.9, capsize=0, zorder=5)
            g = numbers["cases"][str(c5)]["o_case5_obs_over_model_binmean"]
            text = (f"IUVS / model (bin means): {g['geomean']:.2f} (geo. mean),\n"
                    f"{min(r5):.2f}-{max(r5):.2f} over the bins; predicted by the fit {g['predicted_geomean']:.2f}\n"
                    + text)
        ax.text(0.97, 0.96, text, transform=ax.transAxes, ha="right", va="top", color=INK, fontsize=7.0,
                linespacing=1.45)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(*X_LIM)
        ax.set_ylim(0.4 * ylo, 12.0 * yhi)
        H.style(ax)
        what = "O I 130.4 nm vs MAVEN/IUVS" if c == c5 else "O I 130.4 nm"
        ax.set_title(title(c, what, n), color=INK, loc="left", fontsize=8.4, linespacing=1.35)
        ax.set_xlabel("Distance from comet (10$^5$ km)")
        if j == 0:
            ax.set_ylabel("O I 130.4 nm brightness (R)")
    handles = [Line2D([], [], color=INK, ls="none", marker="o", ms=3.5, label="MAVEN/IUVS"),
               Line2D([], [], color=S1, lw=1.6, label=f"AMPS inner-boundary fix + SA1982_VH OH loss, Q(H$_2$O) fitted to OH "
                                                      f"(output {out}, NOT converged); O: direct (KVB2022) + Ly$\\beta$ Bowen")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.03), ncol=2, frameon=False, fontsize=8.0)
    fig.text(0.5, 0.004,
             f"Output {out} = iterations {out * H.SAMPLE_LENGTH + 1}-{(out + 1) * H.SAMPLE_LENGTH} of "
             f"{CAMPAIGN.name} (NOT converged), shown as run (not scaled); Q(H$_2$O) was fitted to OH from the "
             "cancelled runs at output 4 ('predicted by the fit' = those runs x fit factor).\nAll curves edge-corrected, "
             "shown to 2/3 of the domain radius. Model O comes from H$_2$O only (no CO/CO$_2$ photolysis or electron "
             "impact in I3); optically thin brightness.",
             ha="center", va="bottom", color=MUTED, fontsize=7.8, linespacing=1.4)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(REPORT / f"plots/qfit_OH_H_O_vs_maven_out{out}.png")
    plt.close(fig)

    # ------------------------------------------------ figure 2: case 5 O vs IUVS
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(14.2, 6.0), dpi=300, gridspec_kw={"width_ratios": [1.35, 1]})
    xs = b5 / 1e8
    y = corrected(n5, OKEY)
    ok = (xs >= H.OBS_X[0]) & (xs <= H.OBS_X[1]) & (b5 <= 2.0 / 3.0 * n5["r_edge"]) & np.isfinite(y) & (y > 0)
    ax.plot(xs[ok], y[ok], color=S1, lw=1.6,
            label=f"AMPS fix + SA1982_VH, Q(H$_2$O) = {H.sci(n5['q_h2o'])} s$^{{-1}}$ (fitted to OH): "
                  "direct (KVB2022) + Ly$\\beta$ Bowen")
    ax.plot(xo, mean5, ls="none", marker="s", ms=5.2, mfc="white", mec=S1, mew=1.3, label="   same, mean over each bin (1-D)")
    ax.plot([float(r["distance_1e5_km"]) for r in slide], [float(r["radiance_R"]) for r in slide],
            color=H.SLOTS[7], lw=1.3, label=f"IUVS team model on the slide (Q(H$_2$O) = {H.sci(H.SLIDE_Q_H2O)} s$^{{-1}}$)")
    ax.errorbar(xo, yo, xerr=[xo - np.maximum(lo, H.OBS_X[0]), hi - xo], fmt="o", color=INK, ms=4.6, lw=1.1,
                capsize=0, zorder=5, label=f"MAVEN/IUVS {H.OBS_TIME} (digitized; bars: bins)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*H.OBS_X); ax.set_ylim(*H.OBS_Y)
    ax.set_xlabel("Distance from comet (10$^5$ km)")
    ax.set_ylabel("O I 130.4 nm brightness (R)")
    ax.set_title(f"(a) Case {c5} O I 130.4 nm vs MAVEN/IUVS: r$_h$ = {n5['rh']:.2f} AU, "
                 f"$\\dot r$ = {n5['rdot']:+.1f} km/s; output {out}, NOT converged",
                 loc="left", color=INK, fontsize=8.8)
    ax.text(0.035, 0.035,
            "AMPS edge-corrected, shown as run (not scaled), optically thin. Model O comes from H$_2$O only\n"
            "(no CO/CO$_2$ photolysis or electron impact in I3), which would make the model brighter.",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.4, color=INK, linespacing=1.45)
    H.style(ax)
    ax.legend(loc="upper right", frameon=False, fontsize=6.8, handlelength=2.6)
    axr.plot(xo, r5, ls="-", lw=0.9, marker="s", ms=4.8, mfc="white", mec=S1, color=S1,
             label=f"direct (KVB2022) + Bowen (geo. mean {geo(r5):.2f})")
    axr.axhline(pc5["now_KVB2022"]["geomean"], color=S1, lw=0.8, ls=(0, (1, 1.5)),
                label=f"predicted by the fit: {pc5['now_KVB2022']['geomean']:.2f}")
    axr.axhline(1.0, color=MUTED, lw=0.8)
    axr.set_xscale("log"); axr.set_yscale("log")
    axr.set_xlim(0.3, 15)
    H.ratio_axis(axr, min(0.4, 0.8 * r5.min()), max(2.0, 1.25 * r5.max()))
    axr.set_xlabel("Distance from comet (10$^5$ km)")
    axr.set_ylabel("observed / AMPS")
    axr.set_title("(b) IUVS / AMPS per bin, edge-corrected bin means", loc="left", color=INK, fontsize=8.8)
    H.style(axr)
    axr.legend(loc="lower left", frameon=False, fontsize=6.8)
    fig.text(0.5, 0.005,
             f"Output {out} = iterations {out * H.SAMPLE_LENGTH + 1}-{(out + 1) * H.SAMPLE_LENGTH} of "
             f"{CAMPAIGN.name}, NOT converged. 'Predicted by the fit' = cancelled runs at output 4 x fit factor.\n"
             "Observation and emission model as in o_case5_vs_iuvs.png of the O report (digitized slide; resonant scattering "
             "+ Ly$\\beta$ Bowen only; bin means along the slit).",
             ha="center", va="bottom", color=MUTED, fontsize=7.3, linespacing=1.4)
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    fig.savefig(REPORT / f"plots/qfit_o_case5_vs_iuvs_out{out}.png")
    plt.close(fig)

    (REPORT / f"data/qfit_numbers_out{out}.json").write_text(json.dumps(numbers, indent=1))
    for c in CASES:
        e = numbers["cases"][str(c)]
        print(f"case {c}: model/obs OH {e['OH_model_over_obs']['geomean']:.3f}, H {e['H_model_over_obs']['geomean']:.3f} "
              f"(predicted {e['predicted']['h_model_over_obs']:.3f}); O(KVB2022+Bowen) at 1e5 km "
              f"{e['O_KVB2022_bowen_edge_corrected_R']['1e+05 km']:.3g} R")
    print(f"case 5 O obs/model (bin means): {', '.join(f'{v:.2f}' for v in r5)}; geo. mean {geo(r5):.3f}")
    print(f"wrote plots/qfit_OH_H_O_vs_maven_out{out}.png, plots/qfit_o_case5_vs_iuvs_out{out}.png, "
          f"data/qfit_numbers_out{out}.json, data/qfit_o_case5_vs_iuvs_out{out}.csv")


if __name__ == "__main__":
    main()

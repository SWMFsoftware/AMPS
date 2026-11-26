#!/usr/bin/env python3
"""O I 130.4-nm g-factor of the I3 O brightness comparison: direct (KVB2022) + Ly-beta Bowen.

2026-11-25, user request: "prepare a plot of the g-factor used for O brightness comparison". Copied
2026-11-26 from validation_runs/i3_qfit_campaign_20261123/report_20261124/scripts/ to run from git (options
--campaign, --report-dir, --output K|final, --np; the computation is unchanged).

(a) g at 1 AU against heliocentric radial velocity, from the tables compiled into AMPS
    (src/models/resonance_fluorescence/generated/reference_tables.inc), with the case-5 solar
    fluxes; the six cases are marked at their r-dot with their own fluxes.
(b) g of each case at the comet (r-dot, r_h, case solar fluxes): the value AMPS evaluates for the
    comet's bulk motion (CometBrightness.cpp, nominal g).
(c) g the run applied / g of (b), against distance from the comet. AMPS evaluates g for every
    sampled O particle at its own heliocentric radial velocity (comet r-dot + particle velocity);
    the brightness table records the column-weighted mean, "Effective Photon g Factor".

The g evaluation repeats resonance_fluorescence.cpp: linear interpolation of the uniform tables,
direct = table * F_triplet / F_reference, Bowen = table * F_Lyb, both * (1 AU / r_h)^2.
Runs: <campaign>/sa1982_vh/caseN_np<NP>/production, default campaign validation_runs/i3_qfit_campaign_20261123
(SA1982_VH, Q(H2O) fitted to OH); --output K (default 4) or final (last output of each run).
Writes <report-dir>/plots/o1304_gfactor_<tag>.png and <report-dir>/data/o1304_gfactor_{cases_<tag>.csv,<tag>.json,
1au_case5_fluxes.csv}, tag = out<K> or final; default report directory <campaign>/report_20261124.
    python3 srcHartley2/i3/tools/make_o_gfactor_figure.py [--output K|final] [--campaign DIR] [--report-dir DIR]
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent        # srcHartley2/i3/tools
ROOT = HERE.parents[2]                        # AMPS repository root
CAMPAIGN = ROOT / "validation_runs/i3_qfit_campaign_20261123"   # default; --campaign
REPORT = CAMPAIGN / "report_20261124"                           # default; --report-dir
sys.path.insert(0, str(HERE))
import i3_figure_helpers as H  # noqa: E402

np, plt = H.np, H.plt
INC = ROOT / "src/models/resonance_fluorescence/generated/reference_tables.inc"
KEY_K = "OI_1304_TRIPLET_DIRECT_KVB2022"
KEY_B = "OI_1304_TRIPLET_LYB_BOWEN_SUMER2001"


def inc_array(text, name):
    m = re.search(rf"constexpr double {name}\[\] = \{{(.*?)\}};", text, re.S)
    return np.array([float(x) for x in m.group(1).split(",") if x.strip()])


def inc_const(text, name):
    return float(re.search(rf"constexpr \S+ {name} = ([^;]+);", text).group(1))


def tables():
    text = INC.read_text()
    t = {"kvb": inc_array(text, "OKVB2022DirectPhotonRate1AUS1"),
         "bowen": inc_array(text, "OSUMER2001BowenRatePerLymanBetaFlux1AU"),
         "f_ref": inc_const(text, "OKVB2022ReferenceTripletFlux1AUPhotonM2S")}
    for key, prefix in (("kvb", "OKVB2022"), ("bowen", "OSUMER2001")):
        v0 = inc_const(text, f"{prefix}VelocityMinimumMPS") / 1e3
        dv = inc_const(text, f"{prefix}VelocityStepMPS") / 1e3
        t[f"v_{key}"] = v0 + dv * np.arange(len(t[key]))
    return t


def g_direct(t, v, f_oi_m2, rh=1.0):
    v = np.asarray(v, float)
    out = np.interp(v, t["v_kvb"], t["kvb"]) * f_oi_m2 / t["f_ref"] / rh ** 2
    return np.where((v >= t["v_kvb"][0]) & (v <= t["v_kvb"][-1]), out, np.nan)


def g_bowen(t, v, f_lyb_m2, rh=1.0):
    return np.interp(v, t["v_bowen"], t["bowen"]) * f_lyb_m2 / rh ** 2


def rel(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def main():
    global CAMPAIGN, REPORT
    ap = argparse.ArgumentParser(description=" ".join(__doc__.split("\n\n")[0].split()))
    ap.add_argument("--output", default="4", help="output index K of the runs, or 'final' (default 4)")
    ap.add_argument("--campaign", type=Path, default=CAMPAIGN, help="default: %(default)s")
    ap.add_argument("--report-dir", type=Path, default=None, help="default: <campaign>/report_20261124")
    ap.add_argument("--np", type=int, default=128, help="ranks in the run-directory names (default 128)")
    args = ap.parse_args()
    CAMPAIGN = args.campaign.resolve()
    REPORT = (args.report_dir or CAMPAIGN / "report_20261124").resolve()
    final = args.output == "final"
    out = "final" if final else int(args.output)
    tag = "final" if final else f"out{out}"
    t = tables()
    rows, effective = [], {}
    for case in H.CASES:
        run = CAMPAIGN / f"sa1982_vh/case{case}_np{args.np}/production"
        k = max(H.outputs_of(run)) if final else out
        d = H.load(run, case, k)
        inp = H.edge.parse_runtime_input(run / "i3.in", case)
        f_oi = float(inp["brightness.solar_oi1304_triplet_flux_1au_photon_cm-2_s-1"]) * 1e4
        f_lyb = float(inp["brightness.solar_lyb_flux_1au_photon_cm-2_s-1"]) * 1e4
        gk = float(g_direct(t, d["rdot"], f_oi, d["rh"]))
        gb = float(g_bowen(t, d["rdot"], f_lyb, d["rh"]))
        c = d["cols"]["O"]
        if not np.array_equal(c[f"Column Density {KEY_K} [m^-2]"], c[f"Column Density {KEY_B} [m^-2]"]):
            raise SystemExit(f"case {case}: direct and Bowen O columns differ")
        if max(c[f"Excluded Column Fraction {k}"].max() for k in (KEY_K, KEY_B)) > 0:
            raise SystemExit(f"case {case}: particles outside the g-factor velocity domain")
        ek, eb = c[f"Effective Photon g Factor {KEY_K} [s^-1]"], c[f"Effective Photon g Factor {KEY_B} [s^-1]"]
        keep = (d["b"] <= 2.0 / 3.0 * d["r_edge"]) & (c[f"Column Density {KEY_K} [m^-2]"] > 0)
        effective[case] = (d["b"][keep], (ek + eb)[keep] / (gk + gb), ek[keep] / gk, eb[keep] / gb)
        at = {f"{x:g}": float(np.interp(np.log(x * 1e3), np.log(d["b"][keep]), ((ek + eb) / (gk + gb))[keep]))
              for x in (1e3, 1e4, 1e5, 1e6)}
        rows.append({"case": case, "r_h_au": d["rh"], "rdot_km_s": d["rdot"],
                     "f_oi1304_triplet_1au_ph_cm-2_s-1": f_oi / 1e4, "f_lyb_1au_ph_cm-2_s-1": f_lyb / 1e4,
                     "g_direct_kvb2022_s-1": gk, "g_bowen_s-1": gb, "g_total_s-1": gk + gb,
                     "bowen_share": gb / (gk + gb),
                     "g_direct_1au_s-1": gk * d["rh"] ** 2, "g_bowen_1au_s-1": gb * d["rh"] ** 2,
                     "effective_over_bulk_total_at_km": at,
                     "effective_over_bulk_direct_range": [float(effective[case][2].min()), float(effective[case][2].max())],
                     "effective_over_bulk_bowen_range": [float(effective[case][3].min()), float(effective[case][3].max())],
                     "run": rel(run), "output": k, "iterations": d["iterations"],
                     "converged": d["converged"], "brightness_O_sha256": d["sha256"]["O"]})
        print(f"case {case}: r_h={d['rh']:.2f} au rdot={d['rdot']:+.1f} km/s g_direct={gk:.3e} g_bowen={gb:.3e} "
              f"total={gk + gb:.3e} s^-1 (Bowen {100 * gb / (gk + gb):.0f}%); effective/bulk at 1e3..1e6 km "
              + " ".join(f"{v:.3f}" for v in at.values()))

    fig, (ax, axb, axc) = plt.subplots(1, 3, figsize=(16.5, 5.4), dpi=300,
                                       gridspec_kw={"width_ratios": [1.25, 0.95, 1.1]})
    r5 = rows[H.OBS_CASE - 1]
    f5 = (r5["f_oi1304_triplet_1au_ph_cm-2_s-1"] * 1e4, r5["f_lyb_1au_ph_cm-2_s-1"] * 1e4)
    v = np.round(np.arange(-55.0, 55.0001, 0.1), 1)
    gd, gbw = g_direct(t, v, f5[0]), g_bowen(t, v, f5[1])
    ax.plot(v, gd + gbw, color=H.INK, lw=2.0, label="total (used)", zorder=3)
    ax.plot(v, gd, color=H.MUTED, lw=1.3, ls="--", label="direct, KVB2022 (130.217 nm)")
    ax.plot(v, gbw, color=H.MUTED, lw=1.3, ls=":", label="Ly$\\beta$ Bowen (102.576 nm)")
    for edge_v in t["v_kvb"][[0, -1]]:
        ax.axvline(edge_v, color=H.GRID, lw=1.0, zorder=0)
    ax.text(t["v_kvb"][-1] - 0.8, 2e-5, "end of\nKVB2022\ntable", ha="right", va="top", fontsize=7, color=H.MUTED)
    for r in rows:
        k = r["case"]
        g1 = r["g_direct_1au_s-1"] + r["g_bowen_1au_s-1"]
        ax.plot(r["rdot_km_s"], g1, "o", ms=7, color=H.SLOTS[k - 1], mec="white", mew=1.0, zorder=4)
        dy = {4: 0.62, 5: 1.55, 6: 0.62}.get(k, 1.6)
        dx = {4: -2.2, 5: 0.0, 6: 2.2}.get(k, 0.0)
        ax.text(r["rdot_km_s"] + dx, g1 * dy, str(k), ha="center", va="center", fontsize=8, color=H.INK)
    ax.set_yscale("log")
    ax.set_xlim(-55, 55); ax.set_ylim(1e-9, 3e-5)
    ax.set_xlabel("Heliocentric radial velocity (km s$^{-1}$)")
    ax.set_ylabel("g at 1 AU (130.4-nm triplet photons s$^{-1}$ per O atom)")
    ax.set_title("(a) g at 1 AU, case-5 solar fluxes; points: cases 1-6 at their $\\dot r$", loc="left")
    ax.legend(frameon=False, loc="lower center", fontsize=7.5)
    H.style(ax)

    x = np.arange(1, 7)
    for r in rows:
        k, col = r["case"], H.SLOTS[r["case"] - 1]
        # direct left, total centre, Bowen right: in cases 4-6 total and direct differ by 3-4%
        axb.plot(k - 0.2, r["g_direct_kvb2022_s-1"], "s", ms=6.0, mfc="white", mec=col, mew=1.4, zorder=3)
        axb.plot(k, r["g_total_s-1"], "o", ms=8, color=col, mec="white", mew=1.0, zorder=3)
        axb.plot(k + 0.2, r["g_bowen_s-1"], "^", ms=6.5, mfc="white", mec=col, mew=1.4, zorder=3)
        axb.text(k, r["g_total_s-1"] * 1.32, f"{r['g_total_s-1']:.1e}".replace("e-0", "e-"),
                 ha="center", va="bottom", fontsize=7, color=H.INK)
    axb.set_yscale("log")
    axb.set_xlim(0.5, 6.5); axb.set_ylim(5e-9, 3e-5)
    axb.set_xticks(x)
    axb.set_xticklabels([f"{r['case']}\n{r['r_h_au']:.2f} au\n{r['rdot_km_s']:+.1f}" for r in rows])
    axb.set_xlabel("case, r$_h$, $\\dot r$ (km s$^{-1}$)")
    axb.set_ylabel("g at the comet (photons s$^{-1}$ per O atom)")
    axb.set_title("(b) g of each case at $\\dot r$ and r$_h$", loc="left")
    shapes = [H.Line2D([], [], ls="", marker="o", ms=7, color=H.MUTED, mec="white", label="total (used)"),
              H.Line2D([], [], ls="", marker="s", ms=5.5, mfc="white", mec=H.MUTED, mew=1.3, label="direct, KVB2022"),
              H.Line2D([], [], ls="", marker="^", ms=6, mfc="white", mec=H.MUTED, mew=1.3, label="Ly$\\beta$ Bowen")]
    axb.legend(handles=shapes, frameon=False, loc="upper left", fontsize=7.5)
    H.style(axb)

    for case, (b, ratio, _, _) in effective.items():
        axc.plot(b / 1e8, ratio, color=H.SLOTS[case - 1], lw=1.6, label=f"case {case}")
    axc.axhline(1.0, color=H.MUTED, lw=0.8, zorder=0)
    axc.set_xscale("log")
    axc.set_xlim(0.01, 25); axc.set_ylim(0.85, 1.6)
    axc.set_xlabel("Distance from comet (10$^5$ km)")
    axc.set_ylabel("g applied in the run / g of (b)")
    axc.set_title("(c) column-weighted g of the sampled O atoms, " + ("last output of each run" if final else f"output {out}"),
                  loc="left")
    axc.legend(frameon=False, loc="upper left", fontsize=7.5, ncol=2)
    H.style(axc)

    fig.text(0.006, 0.012,
             "O I 130.4-nm triplet photons per ground-state (2p$^4$ $^3$P$_2$) O atom. Direct: Killen, Vervack & Burger "
             "(2022) Fig. 14 g of 130.217 nm (Gladstone 1992 line shape, SORCE SOLSTICE 2009-03-01), x 0.352$^2$ to 1 AU, "
             "x case FISM2 triplet flux / 8.41e9 cm$^{-2}$ s$^{-1}$. Bowen: SOHO/SUMER quiet-Sun Ly$\\beta$ profile "
             "(Curdt et al. 2001)\nnormalized to the case FISM2 102.4-102.8 nm flux, NIST ASD f and branching to 1128.7 nm "
             "(0.288). Both x (1 AU / r$_h$)$^2$. Runs: "
             f"{CAMPAIGN.name} (SA1982_VH, Q(H$_2$O) fitted to OH), "
             + ("last output of each run" if final else f"output {out}") + "; (c) to 2/3 of the domain radius.",
             fontsize=6.8, color=H.MUTED, va="bottom")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    (REPORT / "plots").mkdir(parents=True, exist_ok=True)
    (REPORT / "data").mkdir(exist_ok=True)
    png = REPORT / f"plots/o1304_gfactor_{tag}.png"
    fig.savefig(png)
    with open(REPORT / f"data/o1304_gfactor_cases_{tag}.csv", "w", newline="") as stream:
        w = csv.writer(stream)
        w.writerow(["case", "r_h_au", "rdot_km_s", "f_oi1304_triplet_1au_ph_cm-2_s-1", "f_lyb_1au_ph_cm-2_s-1",
                    "g_direct_kvb2022_s-1", "g_bowen_s-1", "g_total_s-1", "bowen_share",
                    "applied_over_bulk_at_1e3km", "applied_over_bulk_at_1e4km", "applied_over_bulk_at_1e5km",
                    "applied_over_bulk_at_1e6km"])
        for r in rows:
            w.writerow([r["case"], r["r_h_au"], r["rdot_km_s"], f"{r['f_oi1304_triplet_1au_ph_cm-2_s-1']:.10e}",
                        f"{r['f_lyb_1au_ph_cm-2_s-1']:.10e}", f"{r['g_direct_kvb2022_s-1']:.6e}",
                        f"{r['g_bowen_s-1']:.6e}", f"{r['g_total_s-1']:.6e}", f"{r['bowen_share']:.4f}"]
                       + [f"{x:.4f}" for x in r["effective_over_bulk_total_at_km"].values()])
    with open(REPORT / "data/o1304_gfactor_1au_case5_fluxes.csv", "w", newline="") as stream:
        w = csv.writer(stream)
        w.writerow(["# g at 1 AU with the case-5 solar fluxes; photons s^-1 per 2p4 3P2 O atom; "
                    "direct undefined beyond the KVB2022 table (+-48 km/s)"])
        w.writerow(["v_km_s", "g_direct_kvb2022_1au_s-1", "g_bowen_1au_s-1", "g_total_1au_s-1"])
        for vi, a, bw in zip(v, gd, gbw):
            w.writerow([f"{vi:.1f}", "" if np.isnan(a) else f"{a:.6e}", f"{bw:.6e}",
                        "" if np.isnan(a) else f"{a + bw:.6e}"])
    meta = {"tables": rel(INC), "tables_sha256": H.sha256(INC),
            "f_reference_triplet_1au_ph_cm-2_s-1": t["f_ref"] / 1e4, "cases": rows,
            "figure": rel(png)}
    (REPORT / f"data/o1304_gfactor_{tag}.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(png)


if __name__ == "__main__":
    main()

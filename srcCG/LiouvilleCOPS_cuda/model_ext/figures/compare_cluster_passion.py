#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/compare_cluster_passion.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/compare_cluster_passion.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c1e4968f073ab26cd8e8c8e0839a45d5e5ffd184f95c7f7860a3d2b97ee2532b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: early validation figures: the cluster (Voyager, 8 x A100) 12M fits against the passion (2 x A6000) ones.
# RUN (from $APP): .venv/bin/python <this> --run <cluster run dir> --ref-root <passion products> --out <dir> [--set ...]
# ==== end of header ====
"""Cluster (Voyager, 8 x A100) against passion (2 x A6000): the 12M earlier-set or C16 fits.

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python (uses cops_cuda.figures):

    .venv/bin/python <this> --run <restart dir> --ref-root <T>/reference_fits --out <dir> [--set earlier|c16] [--overwrite]

Writes, in the publication format of cops_cuda.figures (.eps and .png, 300 dpi):
* <out>/model_pressure_difference_{time,altitude}: |p_cluster / p_passion - 1| of the fitted model
  pressures (pressure_residuals.dat, 16 digits), NG and RG, one line per fit; late window shaded;
* <out>/parameter_shift_sigma: |fitted_cluster - fitted_passion| / sigma_passion per parameter;
* <out>/cost_table.{txt,json}: costs, optimizer, rms, active bounds, wall and GPU times.
Exact zeros are drawn at the floor FLOOR and listed in the table.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

APP = Path.cwd()
sys.path.insert(0, str(APP))
from cops_cuda import figures  # noqa: E402

SETS = {   # fits, cluster dir under --run, passion dir under --ref-root, table title
    "earlier": (("B29", "N30", "C30", "NC31"), "fits/%s_n2000/seg01", "fit_012m_gpu_%s_n2000_mixed", "earlier set"),
    "c16": (("B29", "N30", "S30", "NS31"), "c16/fits/%s_n2000_mono10_ordered",
            "fit_012m_gpu_%s_n2000_mixed_mono10_ordered", "C16 set (--monotonic-weight 10 --ordered-temperature)"),
}
FLOOR = 1e-17
TBOUNDS = (30.0, 400.0)          # the --tbounds of every fit compared here


def residual_rows(fit_dir):
    """{(iPoint, gauge): (utc, altitude, model pressure)} from pressure_residuals.dat."""
    out = {}
    for line in open(Path(fit_dir) / "pressure_residuals.dat"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        out[(int(f[0]), f[3])] = (f[1], float(f[2]), float(f[7]))
    return out


def fit_result(fit_dir):
    """Scalars and the parameter table of fit_result.txt."""
    text = (Path(fit_dir) / "fit_result.txt").read_text()
    res = {"params": {}}
    for key in ("initial_cost", "final_data_cost", "rms_ln_residual_NG", "rms_ln_residual_RG"):
        m = re.search(r"^%s\s+(\S+)" % key, text, re.M)
        res[key] = float(m.group(1)) if m else None
    m = re.search(r"^optimizer: (.*)\(status (\-?\d+), nfev (\d+), first-order optimality (\S+)\)", text, re.M)
    res["optimizer"] = {"message": m.group(1).strip(), "status": int(m.group(2)), "nfev": int(m.group(3)),
                        "optimality": float(m.group(4))}
    m = re.search(r"^data_jacobian_rank (\d+) of (\d+)", text, re.M)
    res["rank"] = "%s/%s" % (m.group(1), m.group(2))
    block = text.split("# parameter", 1)[1].splitlines()[1:]
    for line in block:
        f = line.split()
        if len(f) != 4:
            break
        res["params"][f[0]] = {"fitted": float(f[2]), "sigma": float(f[3])}
    # the temperatures as printed after the table (with --ordered-temperature the parameters are Tmin, dT75, dTmax)
    t = {k: float(re.search(r"^%s\s+(\S+)\s*$" % k, text, re.M).group(1)) for k in ("Tmax", "T75", "Tmin")}
    res["active_bounds"] = [k for k, v in t.items() if abs(v - TBOUNDS[0]) < 1e-6 or abs(v - TBOUNDS[1]) < 1e-6]
    res["active_bounds"] += [k for k, v in res["params"].items() if k.startswith("c_") and v["fitted"] == 0.0]
    return res


def backend_times(fit_dir):
    try:
        cb = json.load(open(Path(fit_dir) / "fit_manifest.json"))["compute_backend"]
    except (OSError, KeyError, ValueError):
        return {}
    keep = {}
    for key in ("wall_seconds", "call_seconds", "call_counts", "devices"):
        if key in cb:
            keep[key] = cb[key]
    return keep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="runner restart dir (base/, bundle/, fits/<FIT>_n2000/seg01)")
    ap.add_argument("--ref-root", required=True, help="$T/reference_fits")
    ap.add_argument("--out", required=True)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--set", choices=sorted(SETS), default="earlier",
                    help="earlier: fits/<FIT>_n2000/seg01; c16: c16/fits/<FIT>_n2000_mono10_ordered")
    args = ap.parse_args()
    fits, cluster_dir, passion_dir, title = SETS[args.set]
    run, ref, out = Path(args.run), Path(args.ref_root), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cluster = {f: run / (cluster_dir % f) for f in fits}
    passion = {f: ref / (passion_dir % f) for f in fits}

    plt = figures.style()
    base = figures.read_table(run / "base")
    window = figures.late_window(base, bundle=run / "bundle" / "restart_templates.npz")
    colors = {f: figures.PALETTE[i] for i, f in enumerate(fits)}

    # ---- model pressure differences
    diff, zeros = {}, {}
    for f in fits:
        a, b = residual_rows(cluster[f]), residual_rows(passion[f])
        if set(a) != set(b):
            raise SystemExit("%s: point sets differ" % f)
        for gauge in ("NG", "RG"):
            keys = sorted(k for k in a if k[1] == gauge and np.isfinite(a[k][2]) and np.isfinite(b[k][2])
                          and b[k][2] != 0)
            t = np.array([figures.parse_utc(a[k][0]) for k in keys])
            alt = np.array([a[k][1] for k in keys])
            d = np.array([abs(a[k][2] / b[k][2] - 1.0) for k in keys])
            zeros[(f, gauge)] = int((d == 0).sum())
            diff[(f, gauge)] = (t, alt, np.maximum(d, FLOOR), d)
    for xname in ("time", "altitude"):
        fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.6), sharex=True)
        for ax, gauge in zip(axes, ("NG", "RG")):
            for f in fits:
                t, alt, dplot, _ = diff[(f, gauge)]
                x = t if xname == "time" else alt
                ax.plot(x, dplot, lw=1.0, color=colors[f], label=f)
            ax.set_yscale("log")
            nonzero = [d[d > 0].min() for (f2, g2), (_, _, _, d) in diff.items() if g2 == gauge and (d > 0).any()]
            exact = any(zeros[(f2, gauge)] for f2 in fits)
            ax.set_ylim(FLOOR / 2 if exact else min(nonzero) / 3, 1e-1)
            ax.set_ylabel(r"|$p_\mathrm{cluster}/p_\mathrm{passion}-1$|, %s" % gauge)
            if xname == "time":
                ax.axvspan(window["t0"], window["t1"], facecolor=figures.LATE_FILL, edgecolor=figures.LATE_EDGE,
                           lw=0.5, zorder=0)
            else:
                ax.axvspan(window["altitude_min_m"], window["altitude_max_m"], facecolor=figures.LATE_FILL,
                           edgecolor=figures.LATE_EDGE, lw=0.5, zorder=0)
                ax.set_xscale("log")
        axes[0].text(0.01, 0.97, "shaded: %s" % figures.window_label(window).replace("\n", ", "),
                     transform=axes[0].transAxes, fontsize=7, va="top", color="0.15")
        if xname == "time":
            axes[1].set_xlabel("UTC, 2016-09-30")
            import matplotlib.dates as mdates
            axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        else:
            axes[1].set_xlabel("altitude [m]")
            axes[1].invert_xaxis()
        axes[0].legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False)
        fig.subplots_adjust(left=0.13, right=0.97, top=0.90, bottom=0.13)
        fig.text(0.5, 0.01, "fitted model pressures (pressure_residuals.dat, 16 digits), 12M SHAP8, mixed precision, "
                 "2000 evaluations" + ("; exact equality drawn at %.0e" % FLOOR if any(zeros.values()) else ""), ha="center", fontsize=6, color="0.35")
        figures.save(fig, str(out / ("model_pressure_difference_%s" % xname)), overwrite=args.overwrite)
        plt.close(fig)

    # ---- parameter shifts in sigma units
    res = {f: (fit_result(cluster[f]), fit_result(passion[f])) for f in fits}
    names = []
    for f in fits:
        for n in res[f][1]["params"]:
            if n not in names:
                names.append(n)
    equal_level = 1e-7                   # fit_result.txt prints 7 significant digits
    no_sigma = {}
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    for k, f in enumerate(fits):
        c, p = res[f]
        xs, ys = [], []
        for i, n in enumerate(names):
            if n in p["params"] and n in c["params"]:
                s = p["params"][n]["sigma"]
                if s == 0.0:             # no gradient (e.g. an inert template amplitude): no sigma scale
                    no_sigma.setdefault(f, []).append(n)
                    continue
                xs.append(i + (k - 1.5) * 0.18)
                ys.append(max(abs(c["params"][n]["fitted"] - p["params"][n]["fitted"]) / s, equal_level))
        ax.plot(xs, ys, "o", ms=3.2, mew=0, color=colors[f], label=f)
    ax.set_yscale("log")
    ax.set_ylim(equal_level / 3, 1.0)
    ax.axhline(equal_level * 1.8, color="0.6", lw=0.5, ls=(0, (2, 2)))
    ax.text(len(names) - 0.3, equal_level * 2.2, "equal at the 7 printed digits (fit_result.txt)",
            ha="right", va="bottom", fontsize=6.5, color="0.35")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylabel(r"|$\Delta$ fitted| / $\sigma_\mathrm{passion}$")
    ax.set_xlim(-0.8, len(names) - 0.2)
    ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.tight_layout()
    figures.save(fig, str(out / "parameter_shift_sigma"), overwrite=args.overwrite)
    plt.close(fig)

    # ---- cost table
    rows = []
    for f in fits:
        c, p = res[f]
        rows.append({"fit": f,
                     "initial_cost": {"cluster": c["initial_cost"], "passion": p["initial_cost"]},
                     "final_data_cost": {"cluster": c["final_data_cost"], "passion": p["final_data_cost"],
                                         "relative_difference": c["final_data_cost"] / p["final_data_cost"] - 1},
                     "optimizer": {"cluster": c["optimizer"], "passion": p["optimizer"]},
                     "rms_ln_NG": {"cluster": c["rms_ln_residual_NG"], "passion": p["rms_ln_residual_NG"]},
                     "rms_ln_RG": {"cluster": c["rms_ln_residual_RG"], "passion": p["rms_ln_residual_RG"]},
                     "rank": {"cluster": c["rank"], "passion": p["rank"]},
                     "active_bounds": {"cluster": c["active_bounds"], "passion": p["active_bounds"]},
                     "max_model_pressure_rel_diff": {g: float(diff[(f, g)][3].max()) for g in ("NG", "RG")},
                     "exactly_equal_model_pressures": {g: zeros[(f, g)] for g in ("NG", "RG")},
                     "not_in_sigma_plot_sigma_passion_zero": {n: {"cluster": c["params"][n]["fitted"], "passion": p["params"][n]["fitted"]}
                                                              for n in no_sigma.get(f, [])},
                     "compute_backend": {"cluster": backend_times(cluster[f]), "passion": backend_times(passion[f])}})
    (out / "cost_table.json").write_text(json.dumps(rows, indent=1) + "\n")
    lines = ["# 12M SHAP8 %s, mixed precision, --tbounds 30 400, max_nfev 2000: cluster" % title + " (Voyager 8 x A100 80GB)"
             " vs passion (2 x RTX A6000)",
             "# fit  initial_cost(cl/pa)      final_data_cost cluster  passion        rel.diff    nfev/status(cl,pa)"
             "  rmsNG(cl/pa)   rmsRG(cl/pa)   active bounds (cl | pa)      max|dp/p| NG, RG"]
    for r in rows:
        lines.append("%-5s %.6e/%.6e  %.6e  %.6e  %+.2e  %d/%d, %d/%d   %.4f/%.4f  %.4f/%.4f  %-12s | %-12s  %.1e, %.1e" % (
            r["fit"], r["initial_cost"]["cluster"], r["initial_cost"]["passion"],
            r["final_data_cost"]["cluster"], r["final_data_cost"]["passion"], r["final_data_cost"]["relative_difference"],
            r["optimizer"]["cluster"]["nfev"], r["optimizer"]["cluster"]["status"],
            r["optimizer"]["passion"]["nfev"], r["optimizer"]["passion"]["status"],
            r["rms_ln_NG"]["cluster"], r["rms_ln_NG"]["passion"], r["rms_ln_RG"]["cluster"], r["rms_ln_RG"]["passion"],
            ",".join(r["active_bounds"]["cluster"]) or "-", ",".join(r["active_bounds"]["passion"]) or "-",
            r["max_model_pressure_rel_diff"]["NG"], r["max_model_pressure_rel_diff"]["RG"]))
    if all(r["optimizer"][h]["nfev"] >= 2000 and r["optimizer"][h]["status"] == 0 for r in rows for h in ("cluster", "passion")):
        lines.append("# none converged: every fit stopped at max_nfev 2000 (status 0)")
    for r in rows:
        for n, v in r["not_in_sigma_plot_sigma_passion_zero"].items():
            lines.append("# %s %s: sigma_passion = 0 (no gradient), not in parameter_shift_sigma; fitted cluster %.7e, passion %.7e"
                         % (r["fit"], n, v["cluster"], v["passion"]))
    (out / "cost_table.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/analysis/positive_lambda_scan_12m/scan_summary.py (on Voyager: /data001/highresolu_vt/cops_output/analysis/positive_lambda_scan_12m/scan_summary.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b428b99b6b93235b980d571cea4a05245575d750514b1bbeaab68b431098dd3e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: analysis/figure script. Table (costs, rms, late-window rms, D, kappa, cliff factor, grain, floor and
#   production fractions) and L-curve of the 12M positive-activity lambda scan. 2026-11-07.
# RUN (from <AMPS>/srcCG/LiouvilleCOPS_cuda): .venv/bin/python scan_summary.py --runs <positive dir> --area area.txt
#   --out <dir> [--overwrite]   (area.txt: awk '!/^#/{print $5}' <seg>/activity_facets.dat)
# ==== end of header ====
"""Summary and L-curve of the 12M positive-activity lambda scan (user decisions 2026-11-07: A = exp(Y a + c T_cliff),
c >= 0, night/shadow ratio 0.02, H2O floor 0.02 q_max A_orig, prior = emission-weighted projection of ln A_orig,
area-weighted log Tikhonov D, lambda 0, 1e-5 ... 1e-1 "the user picks from the L-curve"; original temperature law held,
kappa free, grain, crossing weight 10, W10x5 from 10:33:10, gauge offset +1.5 m, max_nfev 4000).

Reads only the fit products (fit_result.txt, grain_term.txt, activity_exp_facets.npz, late_window_model.dat) of
    runs/shap8_012m_lategrid_dz1.5/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam<L>/seg<NN>/  (last segment)
and a facet-area column (awk '!/^#/{print $5}' <any seg>/activity_facets.dat > area.txt; the mesh is the same in all).

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda (cops_cuda.figures):
    .venv/bin/python -W ignore scan_summary.py --runs <positive dir> --area area.txt --out <dir> [--overwrite]
Writes <out>/scan_table.txt and <out>/lcurve.{eps,png} (300 dpi).
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd()))
from cops_cuda import figures  # noqa: E402

LATE = "2016-09-30T10:33:10"


def val(text, key):
    m = re.search(r"^%s\s+(\S+)" % re.escape(key), text, re.M)
    return float(m.group(1)) if m else np.nan


def param(text, name):
    """fitted value and sigma of a parameter row 'name initial fitted sigma'"""
    m = re.search(r"^%s\s+\S+\s+(\S+)\s+(\S+)" % re.escape(name), text, re.M)
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, type=Path)
    ap.add_argument("--area", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    area = np.loadtxt(a.area)
    rows = []
    for d in sorted(a.runs.glob("C30_n4000_exp_cliff_w10x5_grain_cross10_lam*")):
        lam = d.name.rsplit("_lam", 1)[1]
        seg = sorted(d.glob("seg[0-9][0-9]"))[-1]
        t = (seg / "fit_result.txt").read_text()
        g = (seg / "grain_term.txt").read_text()
        opt = re.search(r"^optimizer: (.*?) \(status (\S+), nfev (\d+), first-order optimality (\S+)\)", t, re.M)
        z = np.load(seg / "activity_exp_facets.npz")
        A, seen, A0 = z["A"], z["seen"], z["F_min_over_qmax"] / 0.02
        w_s = area * seen
        # late window: observation rows (kind late_obs) used in the fit, model = table column (W10x5 surface + grain)
        L = [l.split() for l in (seg / "late_window_model.dat").read_text().splitlines() if not l.startswith("#")]
        L = [r for r in L if r[3] == "late_obs"]
        late = {}
        et0 = min(float(x[2]) for x in L)              # first late observation (10:33:10, the grain onset)
        for gauge, im, ir, iu in (("NG", 8, 10, 12), ("RG", 9, 11, 13)):
            r = np.array([np.log(float(x[im]) / float(x[ir])) for x in L if x[iu] == "1" and float(x[ir]) > 0])
            late[gauge] = float(np.sqrt(np.mean(r ** 2)))
            # the set of the late-response forward check (summary_after_grain_onset.txt): 18-s window after the onset
            r = np.array([np.log(float(x[im]) / float(x[ir])) for x in L
                          if x[iu] == "1" and float(x[ir]) > 0 and float(x[4]) - 18.0 >= et0])
            late[gauge + "30"] = float(np.sqrt(np.mean(r ** 2)))
            if gauge == "NG":
                r = np.array([np.log(float(x[im]) / float(x[ir])) for x in L
                              if x[iu] == "1" and float(x[ir]) >= 8e-6 and float(x[4]) - 18.0 >= et0])
                late["NG30hi"] = float(np.sqrt(np.mean(r ** 2)))
                late["n30"] = int(sum(1 for x in L if float(x[4]) - 18.0 >= et0))
        c, cs = param(t, "c_cliff_ge_60deg")
        k, ks = param(t, "kappa")
        rows.append(dict(
            lam=float(lam), name=lam, seg=seg.name, cost=val(t, "final_cost"), data=val(t, "final_data_cost"),
            reg=val(t, "final_regularization_cost"), rmsNG=val(t, "rms_ln_residual_NG"), rmsRG=val(t, "rms_ln_residual_RG"),
            lateNG=late["NG"], lateRG=late["RG"], nlate=len(L), late30NG=late["NG30"], late30NGhi=late["NG30hi"],
            late30RG=late["RG30"], n30=late["n30"], D=val(t, "deviation_from_prior_D"),
            peff=val(t, "effective_number_of_parameters"), term=opt.group(1), status=int(opt.group(2)),
            nfev=int(opt.group(3)), opt=float(opt.group(4)), c=c, c_sig=cs, kappa=k, kappa_sig=ks,
            grainA=val(g, "A_Pa") * 1e6, grainB=val(g, "B_Pa_per_s") * 1e8,
            floor_seen=float(np.sum(w_s[A <= 0.02 * A0]) / np.sum(w_s)),
            below_seen=float(np.sum(w_s[A < A0]) / np.sum(w_s)),
            prod_seen=float(np.sum((area * A)[seen]) / np.sum((area * A0)[seen])),
            prod_unseen=float(np.sum((area * A)[~seen]) / np.sum((area * A0)[~seen])),
            smax=float(np.max(z["s"])), smin=float(np.min(z["s"]))))
    rows.sort(key=lambda r: r["lam"])
    a.out.mkdir(parents=True, exist_ok=True)
    tab = a.out / "scan_table.txt"
    if tab.exists() and not a.overwrite:
        sys.exit("exists: %s" % tab)
    with open(tab, "w") as f:
        f.write("# 12M positive-activity lambda scan, %s (scan_summary.py)\n" % a.runs)
        f.write("# cost = data + regularization (fitter units); rms = rms ln(model/obs) over all selected points; late = over the "
                "late observations >= %s UTC (W10x5 surface + grain);\n# late30 = over the late observations whose 18-s window starts after the first late observation (from 10:33:30, the set of figures_cluster/12M/model_vs_observations_late_response/summary_after_grain_onset.txt), NGhi without the NG readings < 8e-6 Pa\n" % LATE)
        f.write("# D = rms ln(A / A_prior) area weighted (nepers); c = ln(cliff factor); grain A [uPa], B [uPa/100 s];\n"
                "# floor_seen = seen-area fraction with A <= 0.02 A_orig (the H2O floor binds for every record); below_seen = "
                "seen area with A < A_orig;\n# prod_seen/unseen = sum area A / sum area A_orig on the seen / unseen facets "
                "(unseen facets do not reach the gauges)\n")
        keys = ["name", "seg", "cost", "data", "reg", "rmsNG", "rmsRG", "lateNG", "lateRG", "n30", "late30NG", "late30NGhi", "late30RG", "D", "peff", "nfev", "status",
                "opt", "c", "c_sig", "kappa", "kappa_sig", "grainA", "grainB", "floor_seen", "below_seen", "prod_seen",
                "prod_unseen", "smin", "smax", "term"]
        f.write("#" + " ".join(keys) + "\n")
        for r in rows:
            f.write(" ".join(r[k] if isinstance(r[k], str) and k != "term" else ("'%s'" % r[k] if k == "term" else "%.6g" % r[k])
                             for k in keys) + "\n")
    print(tab.read_text())
    plt = figures.style()
    fig, axs = plt.subplots(1, 2, figsize=(figures.DOUBLE_COLUMN_IN, 3.0))
    pos = [r for r in rows if r["lam"] > 0]
    zero = [r for r in rows if r["lam"] == 0]
    ax = axs[0]
    ax.plot([r["D"] for r in pos], [r["data"] for r in pos], "-o", color=figures.PALETTE[0], ms=4, lw=1.2, label="lambda > 0")
    for r in pos:
        ax.annotate("%g" % r["lam"], (r["D"], r["data"]), textcoords="offset points", xytext=(5, 3), fontsize=6.5)
    for r in zero:
        ax.plot([r["D"]], [r["data"]], "s", color=figures.PALETTE[1], ms=4,
                label="lambda = 0 (max_nfev, not converged)" if r["status"] == 0 else "lambda = 0")
        ax.annotate("0", (r["D"], r["data"]), textcoords="offset points", xytext=(5, 3), fontsize=6.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.7, 1.2e4)
    ax.set_xlabel("D = rms ln(A / A_prior), area weighted [nepers]")
    ax.set_ylabel("data cost")
    ax.set_title("L-curve", fontsize=8)
    ax.legend(fontsize=6, loc="upper right")
    ax = axs[1]
    lx = [r["lam"] for r in pos]
    for key, lab, c, mk in (("rmsNG", "NG, all points", figures.PALETTE[0], "o"), ("rmsRG", "RG, all points", figures.PALETTE[1], "s"),
                            ("lateNG", "NG, late window", figures.PALETTE[2], "^"), ("lateRG", "RG, late window", figures.PALETTE[3], "v")):
        ax.plot(lx, [r[key] for r in pos], "-", marker=mk, color=c, ms=3.5, lw=1.0, label=lab)
        for r in zero:
            ax.plot([min(lx) / 10], [r[key]], marker=mk, color=c, ms=3.5, ls="none")
    ax.set_xscale("log")
    ticks = [min(lx) / 10] + lx
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0"] + ["%g" % x for x in lx], fontsize=6.5)
    ax.set_xlabel("lambda (leftmost point: lambda = 0)")
    ax.set_ylabel("rms ln(model / obs)")
    ax.set_ylim(0.06, 0.34)
    ax.set_title("residuals against the COPS observations", fontsize=8)
    ax.legend(fontsize=6, loc="upper right", ncol=2, title="late window: >= 10:33:10 UTC", title_fontsize=6)
    fig.suptitle("SHAP8 12M, positive activity exp(Y a + c T_cliff), floor 0.02, W10x5 late window, gauge +1.5 m Z", fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    print("wrote", figures.save(fig, str(a.out / "lcurve"), overwrite=a.overwrite))


if __name__ == "__main__":
    main()

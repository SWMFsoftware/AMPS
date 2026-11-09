#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-08, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/positive_cycles_table.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/positive_cycles_table.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: bb5b79867cebe23356e0c6ae939f09f50c7f0a64cb8ef11fb95d2b61aff10211
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure helper. Table of every stage of the positive-activity cycles (first-scan fit, cN fit, cN fine tuning):
#   costs, rms, late-window rms, D, kappa, grain, nfev, stop. 2026-11-08.
# USED BY: figures/figures_12m_positive_cycles.sh.
# ==== end of header ====
"""Table of the 12M positive-activity fit / fine-tuning cycles (user decisions 2026-11-07/08: ratio 0.002, three cycles
of fit + per-element fine tuning, fine tuning only on the facets seen exclusively in the last 5 min).

    positive_cycles_table.py <positive_r0.002 dir> <first-scan positive dir> <output file>

One row per stage and lambda: first-scan fit (ratio 0.02), c1 fit, c1 ft, c2 fit, c2 ft, c3 fit, c3 ft (whatever
exists and is complete). Fits: fit_result.txt of <fit>/final; fine tunings: finetune_result.txt of <ft> (top level,
written on completion). Read only.
"""
import re
import sys
from pathlib import Path

run, scan, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
LAMS = ["0", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1"]


def fit_row(d):
    t = (d / "fit_result.txt").read_text()
    v = lambda k: float(re.search(r"^%s\s+(\S+)" % k, t, re.M).group(1))  # noqa: E731
    p = lambda k: float(re.search(r"^%s\s+\S+\s+(\S+)" % re.escape(k), t, re.M).group(1))  # noqa: E731
    opt = re.search(r"^optimizer: (.*?) \(status \S+, nfev (\d+)", t, re.M)
    lw = [l.split() for l in (d / "late_window_model.dat").read_text().splitlines() if not l.startswith("#")]
    import math
    late = {}
    for g, im, ir, iu in (("NG", 8, 10, 12), ("RG", 9, 11, 13)):
        r = [math.log(float(x[im]) / float(x[ir])) for x in lw if x[3] == "late_obs" and x[iu] == "1" and float(x[ir]) > 0]
        late[g] = (sum(a * a for a in r) / len(r)) ** 0.5
    return dict(total=v("final_cost"), data=v("final_data_cost"), reg=v("final_regularization_cost"),
                ng=v("rms_ln_residual_NG"), rg=v("rms_ln_residual_RG"), lng=late["NG"], lrg=late["RG"],
                D=v("deviation_from_prior_D"), kappa=p("kappa"), grain=p("grain_A_uPa"), nfev=int(opt.group(2)),
                stop="ftol" if "ftol" in opt.group(1) else ("max_nfev" if "maximum" in opt.group(1) else opt.group(1)[:20]))


def ft_row(d):
    t = (d / "finetune_result.txt").read_text()
    f = lambda k: float(re.search(r"^%s\s+initial\s+\S+\s+final\s+(\S+)" % k, t, re.M).group(1))  # noqa: E731
    stop = re.search(r"^stop (\S+)", t, re.M).group(1)
    nf = re.search(r"^nfev\S*\s+(\d+)", t, re.M) or re.search(r"(\d+) evaluations", t)
    D = re.search(r"^D\S*\s+initial\s+\S+\s+final\s+(\S+)", t, re.M)
    return dict(total=f("cost_total"), data=f("cost_data"), reg=f("cost_regularization"), ng=f("rms_ln_NG"),
                rg=f("rms_ln_RG"), lng=f("late_window_rms_ln_NG"), lrg=f("late_window_rms_ln_RG"),
                D=float(D.group(1)) if D else float("nan"), kappa=float("nan"), grain=f("grain_A_uPa"),
                nfev=int(nf.group(1)) if nf else -1, stop=stop)


rows = []
for L in LAMS:
    s = scan / ("C30_n4000_exp_cliff_w10x5_grain_cross10_lam%s/seg01" % L)
    if (s / "fit_result.txt").exists():
        rows.append((L, "scan fit (ratio 0.02)", fit_row(s)))
    for c in (1, 2, 3):
        f = run / ("lam%s/c%d/fit/final" % (L, c))
        if (f / "fit_result.txt").exists():
            rows.append((L, "c%d fit" % c, fit_row(f)))
        ft = run / ("lam%s/c%d/ft" % (L, c))
        if (ft / "finetune_result.txt").exists():
            rows.append((L, "c%d fine tuning" % c, ft_row(ft)))
with open(out, "w") as fh:
    fh.write("# 12M positive-activity cycles: %s (first scan %s)\n" % (run, scan))
    fh.write("# costs: fitter units (total = data + regularization; crossing 0 in all); rms = rms ln(model/obs) over the "
             "selected points; late = over the late observations >= 10:33:10 UTC; D nepers; kappa of fits (fine tunings "
             "hold it); grain A uPa; nfev; stop\n")
    fh.write("#%-5s %-22s %11s %11s %11s %7s %7s %7s %7s %8s %8s %7s %5s %s\n" % (
        "lam", "stage", "total", "data", "reg", "rmsNG", "rmsRG", "lateNG", "lateRG", "D", "kappa", "grainA", "nfev", "stop"))
    for L, st, r in rows:
        fh.write("%-6s %-22s %11.5e %11.5e %11.5e %7.4f %7.4f %7.4f %7.4f %8.3f %8.3g %7.3f %5d %s\n" % (
            L, st, r["total"], r["data"], r["reg"], r["ng"], r["rg"], r["lng"], r["lrg"], r["D"], r["kappa"], r["grain"],
            r["nfev"], r["stop"]))
print(out.read_text())

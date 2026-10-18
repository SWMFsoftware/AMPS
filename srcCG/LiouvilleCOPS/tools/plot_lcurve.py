#!/usr/bin/env python3
"""Plot the L-curve of tools/fit_cops.py (<out>/lcurve.txt): misfit vs. deviation from the prior, and the
fitted temperature, CO2 scale and activity diagnostics as functions of lambda.

DATA SOURCES
  <out>/lcurve.txt written by tools/fit_cops.py with several --lambda values: one row per lambda with
  rmsNG, rmsRG, rho = sqrt(msNG + msRG) (log misfit of the COPS NG and RG data selected in the base run), D (relative
  deviation of the H2O activity map from the srcCG prior), p_eff, Tmax, T75, Tmin, kappa, activity diagnostics,
  nfev, status (0: maximum number of evaluations reached), optimality; the header line gives lambda at the
  maximum curvature of (ln rho, ln D).

usage: plot_lcurve.py <L-curve directory or lcurve.txt> [--save file.png]
"""
import argparse
import os
import re

from cops_paths import prepare_output, product_input

import numpy as np
import matplotlib

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("path")
ap.add_argument("--save", default=None, help="save under cops_output (or absolute path outside LiouvilleCOPS) instead of opening a window")
args = ap.parse_args()
if args.save:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

source = product_input(args.path)
fname = source if source.endswith(".txt") else os.path.join(source, "lcurve.txt")
names, lam_c = None, None
for line in open(fname):
    if line.startswith("#lambda"):
        names = line[1:].split()
    m = re.search(r"maximum curvature .* lambda = (\S+)", line)
    if m:
        lam_c = float(m.group(1))
d = np.loadtxt(fname)
c = {n: d[:, i] for i, n in enumerate(names)}
lam = c["lambda"]
conv = c["status"] > 0

fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))

a = ax[0]
a.loglog(c["rho"], c["D"], "-", color="0.6", lw=1)
a.loglog(c["rho"][conv], c["D"][conv], "o", color="C0", label="converged")
a.loglog(c["rho"][~conv], c["D"][~conv], "o", mfc="none", color="C0", label="max. evaluations reached")
for x, y, l in zip(c["rho"], c["D"], lam):
    a.annotate("%.0e" % l, (x, y), textcoords="offset points", xytext=(5, 3), fontsize=8)
if lam_c is not None:
    k = np.argmin(np.abs(np.log(lam / lam_c)))
    a.loglog(c["rho"][k], c["D"][k], "*", ms=16, color="C3", label="max. curvature, lambda=%.1e" % lam_c)
a.set_xlabel("data misfit rho = sqrt(msNG + msRG)  [ln units]")
a.set_ylabel("D = rms(A - A_srcCG) / rms(A_srcCG)")
a.set_title("L-curve")
a.legend(fontsize=8)
a.grid(True, which="both", alpha=0.3)

a = ax[1]
for n, col in (("Tmax", "C3"), ("T75", "C1"), ("Tmin", "C0")):
    a.semilogx(lam, c[n], "o-", color=col, label=n)
a.set_xlabel("lambda")
a.set_ylabel("temperature parameter [K]")
a.set_title("analytical temperature model")
a.legend(fontsize=8)
a.grid(True, which="both", alpha=0.3)

a = ax[2]
for n, lab, col in (("kappa", "kappa (CO2 / srcCG)", "C2"), ("activity_ratio", "H2O activity / srcCG (all)", "C0"),
                    ("activity_ratio_seen", "H2O activity / srcCG (seen elements)", "C4"),
                    ("unseen_activity_fraction", "fraction of H2O activity on unseen elements", "C7"),
                    ("rmsNG", "rms ln residual NG", "C3"), ("rmsRG", "rms ln residual RG", "C1")):
    a.semilogx(lam, c[n], "o-", color=col, label=lab, ms=4)
a.set_xlabel("lambda")
a.set_title("CO2 scale, activity, misfit")
a.legend(fontsize=8)
a.grid(True, which="both", alpha=0.3)

for a in ax[1:]:
    if lam_c is not None:
        a.axvline(lam_c, color="C3", ls=":", lw=1)
fig.suptitle(os.path.abspath(fname), fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.95))
if args.save:
    saved = prepare_output(args.save)
    fig.savefig(saved, dpi=110)
    print("Figure written: %s" % saved)
else:
    plt.show()

#!/usr/bin/env python3
"""Regenerate the manual's analytical teaching figure, not simulation output.

The formulas are independent elementary limits used to explain observable and
boundary conventions. This optional editorial script needs NumPy/matplotlib;
the PDF figure is supplied locally and normal TeX builds do not run it.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
# The older constrained_layout spelling also supports the distribution's
# matplotlib 3.1; do not require the newer layout= keyword just for figures.
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)

# Constant first-order parent loss: normalized time x=beta*t.
x = np.linspace(0, 4, 201)
axes[0].plot(x, np.exp(-x), color="#1565a8", label="survival")
axes[0].plot(x, -np.expm1(-x), "--", color="#bd4b27", label="loss")
axes[0].set(xlabel=r"$\beta t$", ylabel="probability", ylim=(0, 1.02))
axes[0].legend(frameon=False, fontsize=7)
axes[0].set_title("Constant parent loss", fontsize=9)

# Outgoing cosine/Knudsen law versus uniform solid angle. These are CDFs
# in mu=cos(theta), not angular PDFs per dtheta or per steradian.
mu = np.linspace(0, 1, 201)
axes[1].plot(mu, mu**2, color="#1565a8", label="cosine")
axes[1].plot(mu, mu, "--", color="#bd4b27", label="uniform solid angle")
axes[1].set(xlabel=r"$\mu=\cos\theta$", ylabel=r"$P(M\leq\mu)$",
            xlim=(0, 1), ylim=(0, 1))
axes[1].legend(frameon=False, fontsize=6.5)
axes[1].set_title("Emission CDFs", fontsize=9)

# Steady spherically symmetric constant-speed outward flow, normalized at R.
# Excludes gravity, returning molecules, collisions and distributed sources.
r = np.geomspace(1, 10, 201)
axes[2].loglog(r, r**-2, color="#1565a8", label="no loss")
axes[2].loglog(r, np.exp(-0.3*(r-1))/r**2, "--", color="#bd4b27",
               label=r"$\beta R/u=0.3$")
axes[2].set(xlabel=r"$r/R$", ylabel=r"$n(r)/n(R)$")
axes[2].legend(frameon=False, fontsize=7)
axes[2].set_title("Radial reference density", fontsize=9)
for ax in axes:
    ax.grid(alpha=0.2)
    ax.tick_params(labelsize=7)
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)
fig.savefig(OUT / "kinetic_reference.pdf", metadata={
    "Title": "Analytical kinetic reference limits",
    "Subject": "Teaching reference curves; not AMPS simulation output",
    "CreationDate": None,
})
fig.savefig(OUT / "kinetic_reference.png", dpi=180)
plt.close(fig)
(OUT / "kinetic_reference.json").write_text(json.dumps({
    "kind": "analytical_reference_not_simulation",
    "parent_loss": {"normalized_time": x.tolist(), "survival": np.exp(-x).tolist()},
    "emission": {"mu": mu.tolist(), "knudsen_cdf": (mu**2).tolist()},
    "radial": {"r_over_R": r.tolist(), "no_loss": (r**-2).tolist(),
               "beta_R_over_u": 0.3,
               "with_loss": (np.exp(-0.3*(r-1))/r**2).tolist()},
}, indent=2) + "\n")
print("Wrote analytical reference PDF, PNG and numeric provenance")

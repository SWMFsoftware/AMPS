#!/usr/bin/env python3
"""Draw AMPS ownership/interfaces, not a compulsory execution sequence.

Applications, shared physical libraries, core algorithms and supporting
services retain different responsibilities. The local PDF is a normal manual
input; this optional NumPy/matplotlib editorial tool is not run by LaTeX.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)
fig, ax = plt.subplots(figsize=(10.2, 7.8))
fig.subplots_adjust(left=0.015, right=0.985, bottom=0.015, top=0.985)
ax.set(xlim=(0, 14), ylim=(0, 10))
ax.axis("off")
navy, teal = "#07254b", "#008e99"


def box(x, y, width, height, title, detail, fill="#edf5f8", edge=teal):
    """Use separate heading/body sizes for legibility at page width."""
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.035,rounding_size=0.08",
        linewidth=1.1, facecolor=fill, edgecolor=edge))
    ax.text(x+width/2, y+height*0.72, title, ha="center", va="center",
            fontsize=10.8, weight="bold", color=navy)
    ax.text(x+width/2, y+height*0.32, detail, ha="center", va="center",
            fontsize=8.7, linespacing=1.35, color="#172d42")


def arrow(start, end, both=False):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="<->" if both else "-|>",
        mutation_scale=12, linewidth=1.1, color=navy,
        shrinkA=3, shrinkB=3))


ax.text(7, 9.73, "AMPS architecture: responsibilities and interfaces",
        ha="center", va="center", fontsize=15, weight="bold", color=navy)
box(0.3, 8.6, 13.4, 0.85, "Configuration and external data",
    "Compile-time deck  |  runtime case  |  source / field tables  |  SPICE geometry",
    fill="#fff6e4", edge="#c68d20")
box(0.3, 7.15, 8.0, 1.05, "Applications and thin adapters",
    "Moon / Enceladus / comets / geospace / SEP\nStandalone driver or SWMF PT/PC wrapper")
box(8.65, 7.15, 5.05, 1.05, "Shared physical models",
    "Rates, source laws, diffusion, backgrounds\nExplicit units, validity and typed outcomes")
arrow((4.3, 8.6), (4.3, 8.2))
arrow((11.1, 8.6), (11.1, 8.2))
arrow((8.3, 7.67), (8.65, 7.67), both=True)

# The containing panel is the core boundary. Arrows communicate dependency,
# not the numerical splitting order chosen by an application's driver.
ax.add_patch(FancyBboxPatch(
    (0.3, 2.5), 13.4, 4.05, boxstyle="round,pad=0.04",
    facecolor="#f7f9fc", edgecolor=navy, linewidth=1.35))
ax.text(7, 6.25, "AMPS core", ha="center", va="center",
        fontsize=13, weight="bold", color=navy)
arrow((4.3, 7.15), (4.3, 6.55))
arrow((11.1, 7.15), (11.1, 6.55))
ax.text(7.1, 6.85, "registered hooks and accepted physical state",
        ha="center", va="center", fontsize=8.3, color=navy)
box(0.6, 4.85, 3.95, 1.05, "AMR mesh and geometry",
    "Tree / cells / surfaces / cut volumes\nOwner and interpolation identities")
box(5.0, 4.85, 3.8, 1.05, "Particle evolution",
    "Records / weights / time steps\nMovers and boundary events")
box(9.25, 4.85, 4.1, 1.05, "Processes and fields",
    "Collisions / chemistry / internal energy\nPrescribed or supported PIC fields")
box(0.6, 3.05, 3.95, 1.05, "Sampling and observers",
    "Collect / complete / normalize\nTrajectories, rays and diagnostics")
box(5.0, 3.05, 3.8, 1.05, "Parallel ownership",
    "Halo exchange / migration\nCommunicators and coherent epochs")
box(9.25, 3.05, 4.1, 1.05, "Persistence and output",
    "Restart / layout / checksums\nPhysical and numerical provenance")
arrow((4.55, 5.37), (5.0, 5.37), both=True)
arrow((8.8, 5.37), (9.25, 5.37), both=True)
for centre in (2.57, 6.9, 11.3):
    arrow((centre, 4.85), (centre, 4.1))
box(0.3, 1.35, 13.4, 0.72, "Common support",
    "Units  |  random streams  |  containers  |  numerical utilities  |  logging and timing")
ax.text(7, 2.27, "support services used throughout the core and its adapters",
        ha="center", va="center", fontsize=8.2, color=navy)
box(0.3, 0.18, 13.4, 0.72, "Products and independent evidence",
    "Mesh / particle observables  |  synthetic measurements  |  restart records  |  test receipts",
    fill="#fff6e4", edge="#c68d20")
arrow((0.18, 3.55), (0.18, 0.55))
ax.text(0.05, 2.8, "output", fontsize=7, rotation=90, va="center")
fig.savefig(OUT / "architecture.pdf", metadata={
    "Title": "AMPS architecture", "CreationDate": None})
fig.savefig(OUT / "architecture.png", dpi=180)
plt.close(fig)
print("Wrote local architecture vector PDF and PNG")

"""Helpers of make_qfit_figures.py: reading AMPS brightness tables, edge correction, the case-5
MAVEN/IUVS O I 130.4-nm profile, palette and axis style.

Copied 2026-11-25 (user request: make the figure script run from git) from two scripts outside
git, without changing what they compute:
- validation_runs/i3_master_campaign3_oxygen_20261119/report_oxygen/scripts/make_o_brightness_figures.py:
  constants, palette and rcParams, sci, style, sha256, loglog_interp, flagged_spans, run_state,
  read_commented_csv, bin_average;
- validation_runs/i3_fixed_campaign_20261122/report_interim_20261123/scripts/make_interim_figures.py:
  SPEC, LINE, bfile, read_tecplot, outputs_of, load, bright, edge_of, status, ratio_axis.
Changed: the case-5 observation files are read from ./inputs (copies), and the modules of
srcHartley2/visualization are found from this file's location.

Units: impact parameter b in m; H and OH brightness in kR, O I 130.4 nm in R; distances on the
plots in 10^5 km.
"""
import csv
import hashlib
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent        # srcHartley2/i3/tools
ROOT = HERE.parents[2]                        # AMPS repository root
DATA = HERE / "inputs"
sys.path.insert(0, str(ROOT / "srcHartley2/visualization"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402,F401  (re-exported for the figure script)

import comet_edge_correction as edge  # noqa: E402
import plot_i3_brightness_comparison as P  # noqa: E402,F401  (MAVEN/IUVS OH and H reference points)

CASES = range(1, 7)
SAMPLE_LENGTH = 10000                     # iterations per output window (i3.in)
MODELS = (
    ("OI_1304_TRIPLET_DIRECT_KVB2022", "O I 130.4 nm, direct (KVB2022)"),
    ("OI_1304_TRIPLET_DIRECT_SUMER2001", "O I 130.4 nm, direct (SUMER2001)"),
    ("OI_1304_TRIPLET_LYB_BOWEN_SUMER2001", "O I 130.4 nm, Ly$\\beta$ Bowen"),
)
KVB, SUMER, BOWEN = (m for m, _ in MODELS)
# Case-5 comparison with the digitized MAVEN/IUVS profile (inputs/README.md).
OBS_CASE = 5
OBS_CSV = DATA / "iuvs_o1304_case5_digitized.csv"
SLIDE_CSV = DATA / "iuvs_o1304_case5_slide_model.csv"
SLIDE_Q_H2O = 1.89e29                     # label of the slide's model curve, s^-1
OBS_TIME = "2025-10-29 17:32"             # time label of the slide panel
OBS_X, OBS_Y = (0.03, 20.0), (1.0, 3000.0)  # axes of the slide panel

# Reference categorical palette, fixed order (as the earlier I3 report figures).
SLOTS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
         "#6250d6", "#e34948")
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
S1, S2 = SLOTS[0], SLOTS[1]
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 9,
                     "font.family": "DejaVu Sans"})

SPEC = {"H": 0, "O": 2, "OH": 3}          # species index of pic.<sp>.brightness.s=<i>
LINE = {"H": "H_LYMAN_ALPHA", "OH": "OH_A2SIGMA_X2PI_00"}
RATIO_TICKS = (0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0)


def sci(x, digits=2):
    """3.16e+28 -> '3.16×10$^{28}$' for matplotlib text."""
    mantissa, exponent = f"{x:.{digits}e}".split("e")
    return f"{mantissa}$\\times$10$^{{{int(exponent)}}}$"


def style(ax):
    ax.grid(True, which="major", color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def loglog_interp(x, xp, fp):
    ok = np.isfinite(fp) & (fp > 0)
    return np.exp(np.interp(np.log(x), np.log(xp[ok]), np.log(fp[ok])))


def flagged_spans(b, flag):
    """Contiguous [b_lo, b_hi] ranges (m) where the tau flag is set."""
    spans, start = [], None
    for i, f in enumerate(flag > 0):
        if f and start is None:
            start = i
        if start is not None and (not f or i == len(flag) - 1):
            stop = i if f else i - 1
            spans.append((b[start], b[stop]))
            start = None
    return spans


def run_state(run, case, out):
    """(last iteration of output window `out`, True if the run converged exactly there)."""
    log = (run / f"case{case}.run.log").read_text(errors="replace")
    m = re.search(r"converged after (\d+) iterations", log)
    window_end = (out + 1) * SAMPLE_LENGTH
    if m and int(m.group(1)) == window_end:
        return int(m.group(1)), True
    return window_end, False


def read_commented_csv(path):
    with open(path) as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith("#")))


def bin_average(b, B, lo_m, hi_m):
    """Mean of B(b) over lo <= b <= hi along the slit (1-D), log-log interpolation
    of the table; below the smallest tabulated b the first value is used."""
    b0 = b[np.isfinite(B) & (B > 0)][0]
    lo = max(lo_m, b0)
    grid = np.geomspace(lo, hi_m, 4001)
    val = loglog_interp(grid, b, B)
    integral = float(np.sum(0.5 * (val[1:] + val[:-1]) * np.diff(grid)))
    if lo_m < b0:
        integral += float(loglog_interp(np.array([b0]), b, B)[0]) * (b0 - lo_m)
    return integral / (hi_m - lo_m)


def bfile(run, sp, out):
    return run / f"PT/plots/pic.{sp}.brightness.s={SPEC[sp]}.out={out}.dat"


def read_tecplot(path):
    """(column dict, aux dict) of an AMPS brightness table, rows sorted by b. The aux keys may
    contain dots (e.g. Chemistry.Variant.OH_LOSS in the SA1982_VH outputs)."""
    names, rows, aux = None, [], {}
    for line in path.read_text().splitlines():
        if line.startswith("VARIABLES"):
            names = re.findall(r'"([^"]+)"', line)
        elif line.startswith("DATASETAUXDATA"):
            m = re.match(r'DATASETAUXDATA ([\w.]+)="([^"]*)"', line)
            if m is None:
                raise SystemExit(f"{path}: unreadable line {line[:80]!r}")
            aux[m.group(1)] = m.group(2)
        elif line.strip() and line.lstrip()[0] in "0123456789-+.":
            rows.append(line.split())
    data = np.array(rows, dtype=float)
    if names is None or data.shape[1] != len(names):
        raise SystemExit(f"{path}: VARIABLES and data columns disagree")
    data = data[np.argsort(data[:, 0])]
    return {n: data[:, i] for i, n in enumerate(names)}, aux


def outputs_of(run):
    """Output indices for which all three brightness files exist."""
    sets = [{int(p.name.rsplit("=", 1)[1][:-4]) for p in (run / "PT/plots").glob(f"pic.{sp}.brightness.s={i}.out=*.dat")}
            for sp, i in SPEC.items()]
    return set.intersection(*sets)


def load(run, case, out, with_edge=True):
    """Brightness tables of output `out` of a run, its runtime inputs and (with_edge) the
    edge-correction factors of srcHartley2/visualization/comet_edge_correction.py."""
    cols = {sp: read_tecplot(bfile(run, sp, out))[0] for sp in SPEC}
    b = cols["O"]["Impact Parameter [m]"]
    for sp in SPEC:
        if not np.array_equal(cols[sp]["Impact Parameter [m]"], b):
            raise SystemExit(f"{run} output {out}: impact parameters differ between species")
    inp = edge.parse_runtime_input(run / "i3.in", case)
    iterations, converged = run_state(run, case, out)
    d = {"run": run, "out": out, "cols": cols, "b": b, "iterations": iterations, "converged": converged,
         "q_h2o": float(inp["source.H2O.rate_s-1"]), "rh": float(inp["heliocentric_distance_au"]),
         "rdot": float(inp["brightness.heliocentric_radial_velocity_km_s"]),
         "variants": edge.variants_from_output(bfile(run, "O", out)),
         "sha256": {sp: sha256(bfile(run, sp, out)) for sp in SPEC}}
    if with_edge:
        factors, inputs = edge.correction_for_run(run, case, out, b, ("H", "OH", "O"))
        d.update(edge=factors, r_edge=inputs.outer_radius_m, tau_oh_s=1.0 / inputs.loss_rate("OH"))
    return d


def bright(d, key):
    """Brightness array: 'H' or 'OH' in kR, an O model name in R, 'O total <direct>' = direct + Bowen in R."""
    if key in LINE:
        return d["cols"][key][f"Brightness {LINE[key]} [kR]"]
    if key.startswith("O total "):
        return d["cols"]["O"][f"Brightness {key[8:]} [R]"] + d["cols"]["O"][f"Brightness {BOWEN} [R]"]
    return d["cols"]["O"][f"Brightness {key} [R]"]


def edge_of(d, key):
    return d["edge"]["O" if key not in LINE else key]


def status(d):
    return f"{d['iterations'] // 1000}k it" + ("" if d["converged"] else ", running")


def ratio_axis(ax, lo, hi):
    """Log ratio axis with plain labels at fixed values (no crowded minor labels)."""
    ax.set_ylim(lo, hi)
    ticks = [t for t in RATIO_TICKS if lo <= t <= hi]
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter([f"{t:g}" for t in ticks]))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())

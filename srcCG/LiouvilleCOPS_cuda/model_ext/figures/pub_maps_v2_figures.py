#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/pub_maps_v2_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/pub_maps_v2_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1125ad0628aa5e19b69dfe54fc009591ad578a251321a3684aaf674252d2e524
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script (version 2): activity or H2O / total production-rate template and ground-track figures coloured
#   only on the seen facets, discrete log classes per panel (--low-cut-area, --template-relief); 2026-11-13: --model-
#   name, multiplier-bound note from fit_restart.npz.
# RUN: jobs/pub_maps_v2.sbatch (GPU; QUANTITY, TIME, LOWCUT, TEMPLATE_RELIEF, MODEL_NAME)
# ==== end of header ====
"""Version 2 of the activity and production-rate template / ground-track figures of ONE fitted state (softplus model
after 5 cycles), with the surface coloured only where the gauges saw it during the fit.

User request 2026-11-12 (second round): "redo the ground track source rate and activity maps: 1) plot final activity;
2) scan the range of the final activity within the shown window and plot the entire range (use log scale, 15-20
intervals in the log scale -- make sure that it looks nice and the plots are correct representation of the model
results); put the colour only on that part of the nucleus that has been seen by RG or NG during the fitting intervals
and keep grey for the rest of the surface; 3) do the same for the ground tracks with the source rate; as CO2 is very
minor, do not specify H2O or CO2 but plot the total; ... add captions to all plots in a separate file; in the template
activity plot do not add the text on the top and also plots in each of the panel in a separate figure file".
Coordinator's specification: seen = the fit's nonneg_facets.npz seen (facet contributes to at least one fit point);
on the seen facets a DISCRETE log10 scale with 15-20 intervals covering the entire range of the non-zero values on
the seen facets DISPLAYED in each panel (scanned per panel), round boundaries; seen facets with value 0 their own
labelled colour; unseen facets relief-shaded grey; class counts checked; total production rate only.

WHAT IS SHOWN
  quantity  activity: A = finetune_state.npz A_final; production: F_total = F_H2O + F_CO2 [m^-2 s^-1] at 10:33:10 /
            10:39:26 UTC (pub_maps_figures.load_fields: the maps law, verified bitwise against the maps' .vtu files).
  surface   not seen during the fit: grey, relief-shaded (cops_cuda.figures.SURFACE_COLOR, the background of the
            late-visibility maps); seen, value exactly 0: ZERO_SEEN (light green), relief-shaded; seen, value > 0:
            the class colour of its log10 value, unshaded (the colour equals the colour bar's).
  classes   per panel: the facets drawn in the panel (face index per pixel of the rasterizer), their seen non-zero
            values -> [min, max] of log10 (with --low-cut-area f > 0: [area-weighted f quantile, max]; user decision
            2026-11-12 for the activity: f = 0.001, the range holding 99.9 % of the seen non-zero area; values below
            the lowest boundary get the lowest class colour) -> boundaries at multiples of a round step (STEPS; integer
            decades when that gives 15-20 classes) -> 15-20 classes; colours: matplotlib 'plasma' at the class centres.
            Class tables
            (facets, area, value range per class, colour) in checks/class_tables.txt; every facet value is checked to
            lie inside its class.
  track     ground track, late window and impact point exactly as in the late-visibility set: groundtrack_data of
            --track-base (the observation-only base), the fit's fields from --base; elements.bin fields checked equal.

CODE REUSE: cops_cuda.figures (AMPS bad75ec, figures.py last changed 47fccd5) unchanged: rasterize, render_ortho
(return_faces), _shade, lonlat, OrthoView, _colored_track, _minute_marks, _scale_bar, gravity_slope, bundle_metadata,
groundtrack_data, save, style. Copied with changes (marked "copy of"): the projection part of render_lonlat (returns
the face image), the panel layouts of templates_figure and groundtrack_figures (per-panel discrete colour bars, no
header text, single-panel figures). Fields, checks and caption helpers: figures_cluster/scripts/pub_maps_figures.py.

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python on a GPU node (jobs/pub_maps_v2.sbatch):
    pub_maps_v2_figures.py --base <fit base> --track-base <observation-only base> --bundle <bundle .npz>
        --ft-dir <c5/ft/final> --fit-dir <c5/fit/final> --maps-dir <chain/maps> --geometry <maps geometry>
        --maps-code <frozen nonneg_maps.py> --out-dir <v2 folder> --quantity activity | --quantity total --time HH:MM:SS
        [--classes 15,20]
Writes into <out-dir> (new files only; staged in <out-dir>/v2_<quantity>.incomplete_<UTC>, moved into place):
  activity/ or production_rate/:  12M_templates_<q>, panels/12M_templates_<q>_{a..f}_*, 12M_groundtrack_<q>_{3d,latlon}
  (.eps, .png, _caption.txt each); checks/<q>_class_tables.txt, checks/<q>_fields_verification.json,
  checks/<q>_manifest.json, checks/slope_bins_<q>.json (per-bin area, seen area, multiplier, seen production share)
"""

import argparse
import datetime as dt
import json
import math
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import pub_maps_figures as M  # noqa: E402

FIG = M.FIG
SCHEMA = "pub_maps_v2_v1"
CMAP_NAME = "plasma"
ZERO_SEEN = np.array([0xA1, 0xD7, 0x6A]) / 255.0           # seen during the fit, value exactly 0
NOT_SEEN = FIG.SURFACE_COLOR                                # grey: not seen during the fit
NMIN, NMAX = 15, 20
STEPS = (0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 100)
ROUND_STEPS = (0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20, 25, 50, 100)   # preferred; the others only when none fits
GLOBAL_W, ZOOM_W = 1800, 1600                               # lon/lat image widths (pixels, before 2x supersampling)
T0 = time.time()


def say(text):
    print("pub_maps_v2 [%6.0f s]: %s" % (time.time() - T0, text), flush=True)


def hexcol(c):
    return "#%02X%02X%02X" % tuple(int(round(255 * x)) for x in c[:3])


# ============================================================================================ classes
def class_edges(lgmin, lgmax):
    """Round class boundaries covering [lgmin, lgmax] with NMIN..NMAX classes (integer-decade steps preferred)."""
    cands = []
    for st in STEPS:
        lo = math.floor(lgmin / st + 1e-9) * st
        hi = math.ceil(lgmax / st - 1e-9) * st
        if hi <= lo:
            hi = lo + st
        n = int(round((hi - lo) / st))
        cands.append((0 if NMIN <= n <= NMAX else 1, 0 if st >= 1 else 1, 0 if st in ROUND_STEPS else 1,
                      abs(n - 18), st, lo, n))
    _, _, _, _, st, lo, n = min(cands)
    return np.array([round(lo + st * k, 6) for k in range(n + 1)]), st


class Scale:
    """Discrete log10 classes of the seen non-zero values displayed in one panel."""

    def __init__(self, name, values, areas=None, low_cut=0.0):
        """low_cut > 0 (user decision 2026-11-12, activity): the class range starts at the area-weighted low_cut
        quantile of the values (the smallest value at which the cumulative seen non-zero area, values ascending,
        reaches low_cut of the total) and ends at the maximum (one-sided cut, the top is kept); values below the lowest
        class boundary get the lowest class colour."""
        import matplotlib
        self.name = name
        self.vmin, self.vmax = float(values.min()), float(values.max())
        self.low_cut = float(low_cut)
        if self.low_cut > 0.0:
            o = np.argsort(values, kind="stable")
            cw = np.cumsum(areas[o]) / areas.sum()
            self.cut_value = float(values[o][min(int(np.searchsorted(cw, self.low_cut)), len(o) - 1)])
        else:
            self.cut_value = self.vmin
        lg = np.log10(values)
        self.edges, self.step = class_edges(float(np.log10(self.cut_value)), float(lg.max()))
        below = lg < self.edges[0] - 1e-9
        self.n_below = int(below.sum())
        self.area_below = float(areas[below].sum()) if areas is not None else float("nan")
        self.area_total = float(areas.sum()) if areas is not None else float("nan")
        self.n = len(self.edges) - 1
        self.colors = matplotlib.colormaps[CMAP_NAME]((np.arange(self.n) + 0.5) / self.n)[:, :3]
        self.denormal = self.vmin < np.finfo(np.float64).tiny
        self.count = int(values.size)

    def index(self, values):
        k = np.floor((np.log10(values) - self.edges[0]) / self.step + 1e-12).astype(np.int64)
        return np.clip(k, 0, self.n - 1)

    def text(self, unit=""):
        return ("%d classes of %g decade%s from $10^{%g}$ to $10^{%g}$%s" % (
            self.n, self.step, "" if self.step == 1 else "s", self.edges[0], self.edges[-1], unit))

    def below_text(self):
        """Caption sentence on the values below the class range (only with a low cut)."""
        if self.low_cut <= 0.0:
            return ""
        return (" Class range: from the area-weighted %g %% quantile of these values (%.3e: the seen non-zero facets "
                "below it hold %g %% of their area) to their maximum, rounded outward to the class step (one-sided cut "
                "at the bottom, the top kept); the %s facets below the lowest class boundary 10^%g (%.3g %% of the seen "
                "non-zero facets drawn here, %.3g %% of their area; down to %.3e) have the colour of the lowest class "
                "(triangle at the low end of the colour bar)." % (
                    100 * self.low_cut, self.cut_value, 100 * self.low_cut, "{:,}".format(self.n_below), self.edges[0],
                    100.0 * self.n_below / self.count, 100.0 * self.area_below / self.area_total, self.vmin))


# relief modulation of the class colours (user request 2026-11-12, production-rate template maps: "plot the colour map
# such that the topology of the nucleus is still seen"): k = strength * (I - I_FLAT) / 0.7 with I the relief intensity of
# the lon/lat renderer (cops_cuda.figures.render_lonlat at shade weight 1: I = 1 - 0.7 (1 - max(0, n.L)), light from the
# north-west 40 deg above the local horizon) and I_FLAT its value on level ground (n.L = sin 40 deg): k < 0 (slopes facing
# away from the light) darkens the class colour (c (1 + k)), k > 0 lightens it toward white (c + k (1 - c)); on level
# ground the colour equals the colour bar's. Hue is kept.
I_FLAT = 1.0 - 0.70 + 0.70 * math.sin(math.radians(40.0))


def relief_modulate(rgb, intensity, strength):
    k = strength * (np.asarray(intensity) - I_FLAT) / 0.70
    k = k[:, None]
    return np.clip(np.where(k < 0, rgb * (1.0 + k), rgb + k * (1.0 - rgb)), 0.0, 1.0)


def relief_check(scale, strength):
    """CIELAB colour differences between adjacent classes at the same illumination (darkest slope, level ground,
    brightest slope) and the lightness range of each class (the shaded class must stay distinguishable)."""
    def lab(c):
        c = np.asarray(c, float)
        lin = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
        M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
        xyz = lin @ M.T / np.array([0.95047, 1.0, 1.08883])
        f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16.0 / 116)
        return np.stack((116 * f[:, 1] - 16, 500 * (f[:, 0] - f[:, 1]), 200 * (f[:, 1] - f[:, 2])), axis=1)
    out = {}
    for name, I in (("darkest_slope", 0.30), ("level_ground", I_FLAT), ("brightest_slope", 1.0)):
        L = lab(relief_modulate(scale.colors, np.full(scale.n, I), strength))
        d = np.linalg.norm(np.diff(L, axis=0), axis=1)
        out[name] = {"min_deltaE_adjacent": float(d.min()), "median_deltaE_adjacent": float(np.median(d))}
    out["strength"] = strength
    return out


def facet_rgb(values, seen, scale, intensity, relief=0.0):
    """Displayed colour per facet (shading included): grey unseen, ZERO_SEEN for seen zeros, class colours."""
    rgb = NOT_SEEN[None, :] * intensity[:, None]
    z = seen & (values == 0)
    rgb[z] = ZERO_SEEN[None, :] * intensity[z, None]
    c = seen & (values > 0)
    rgb[c] = scale.colors[scale.index(values[c])]
    if relief > 0.0:
        rgb[c] = relief_modulate(rgb[c], intensity[c], relief)
    return np.clip(rgb, 0.0, 1.0)


def displayed(face):
    f = np.unique(face)
    return f[f >= 0]


def class_table(label, values, seen, area, disp, scale):
    """Text table + dict: facets displayed in a panel per class; checks every value lies inside its class."""
    s = disp[seen[disp] & (values[disp] > 0)]
    z = disp[seen[disp] & (values[disp] == 0)]
    u = disp[~seen[disp]]
    k = scale.index(values[s])
    lg = np.log10(values[s])
    below = lg < scale.edges[0] - 1e-9
    inside = bool(np.all(((lg >= scale.edges[k] - 1e-9) & (lg < scale.edges[k + 1] + 1e-9)) | (below & (k == 0))))
    rows = []
    for j in range(scale.n):
        m = k == j
        rows.append({"class": j, "log10_from": float(scale.edges[j]), "log10_to": float(scale.edges[j + 1]),
                     "colour": hexcol(scale.colors[j]), "facets": int(m.sum()),
                     "area_m2": float(area[s][m].sum()),
                     "value_min": float(values[s][m].min()) if m.any() else None,
                     "value_max": float(values[s][m].max()) if m.any() else None})
    out = {"panel": label, "displayed_facets": int(disp.size), "seen_nonzero": int(s.size), "seen_zero": int(z.size),
           "unseen": int(u.size), "values_inside_their_class": inside,
           "sum_check": int(sum(r["facets"] for r in rows)) == int(s.size) and
           int(s.size + z.size + u.size) == int(disp.size),
           "scale": {"edges_log10": [float(e) for e in scale.edges], "step": scale.step, "classes": scale.n,
                     "min": scale.vmin, "max": scale.vmax, "denormal_min": bool(scale.denormal),
                     "low_cut_area_fraction": scale.low_cut, "cut_value": scale.cut_value},
           "below_range": {"facets": int(below.sum()), "area_m2": float(area[s][below].sum()),
                           "facet_fraction": float(below.mean()) if below.size else 0.0,
                           "area_fraction": float(area[s][below].sum() / area[s].sum()) if s.size else 0.0,
                           "coloured_as_class": 0},
           "classes": rows}
    L = ["[%s]  displayed facets %d: seen non-zero %d, seen zero %d (colour %s), not seen %d (grey); values inside their "
         "class: %s; counts add up: %s" % (label, disp.size, s.size, z.size, hexcol(ZERO_SEEN), u.size, inside,
                                          out["sum_check"]),
         "  range of the seen non-zero values in the panel: %.6e .. %.6e%s; %d classes of %g decade(s)" % (
             scale.vmin, scale.vmax, " (minimum is a denormal float64 value: below 2.2e-308)" if scale.denormal else "",
             scale.n, scale.step),
         "  class range rule: %s; below the lowest boundary: %d facets (%.4g %% of the seen non-zero facets, %.4g %% of "
         "their area), coloured as class 0 (counted in class 0 below)" % (
             ("area-weighted %g %% quantile %.6e .. maximum (one-sided)" % (100 * scale.low_cut, scale.cut_value))
             if scale.low_cut > 0 else "the entire range (minimum .. maximum)", out["below_range"]["facets"],
             100 * out["below_range"]["facet_fraction"], 100 * out["below_range"]["area_fraction"]),
         "  %-4s %-22s %-8s %10s %14s %14s %14s" % ("cls", "log10 interval", "colour", "facets", "area_m2", "min", "max")]
    for r in rows:
        L.append("  %-4d [%8g, %8g)%s %-8s %10d %14.6e %14s %14s" % (
            r["class"], r["log10_from"], r["log10_to"], "]" if r["class"] == scale.n - 1 else " ", r["colour"],
            r["facets"], r["area_m2"], "-" if r["value_min"] is None else "%.4e" % r["value_min"],
            "-" if r["value_max"] is None else "%.4e" % r["value_max"]))
    if not (inside and out["sum_check"]):
        raise SystemExit("class check failed for %s" % label)
    return out, "\n".join(L)


# ============================================================================================ rendering
def lonlat_faces(V, F, lon_range, lat_range, width, device_id=0, ss=2):
    """Copy of the projection part of cops_cuda.figures.render_lonlat (47fccd5), unchanged: the facet drawn at every
    supersampled pixel (-1: background) of the equirectangular image; the colours are applied by FIG._shade."""
    (lo0, lo1), (la0, la1) = lon_range, lat_range
    W = int(width) * ss
    H = int(round(width * (la1 - la0) / (lo1 - lo0))) * ss
    tri = V[F]
    lon, lat = FIG.lonlat(tri)
    r = np.linalg.norm(tri, axis=-1)
    span = lon.max(axis=1) - lon.min(axis=1)
    seam = span > 180.0
    lon_w = np.where(seam[:, None] & (lon < 0), lon + 360.0, lon)
    lon_all = np.concatenate((lon_w, lon_w[seam] - 360.0))
    lat_all = np.concatenate((lat, lat[seam]))
    r_all = np.concatenate((r, r[seam]))
    index = np.concatenate((np.arange(len(F)), np.nonzero(seam)[0]))
    u = (lon_all - lo0) / (lo1 - lo0) * W
    v = (la1 - lat_all) / (la1 - la0) * H
    w = r_all.max() - r_all
    corner = np.arange(3 * len(index), dtype=np.int64).reshape(-1, 3)
    face_k, _ = FIG.rasterize(u.ravel(), v.ravel(), w.ravel(), corner, W, H, device_id=device_id,
                              max_box=max(W, H) // 2)
    return np.where(face_k < 0, -1, index[np.maximum(face_k, 0)]), ss


def shade(face, rgb, ss):
    return FIG._shade(face, rgb, np.ones(len(rgb)), ss)


# ============================================================================================ drawing helpers
def discrete_colorbar(fig, cax, scale, label, orientation="horizontal", max_labels=7):
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cm = ListedColormap(scale.colors)
    norm = BoundaryNorm(scale.edges, cm.N)
    ext = "neither"
    if getattr(scale, "n_below", 0) > 0:
        cm.set_under(scale.colors[0])             # values below the range: the lowest class colour (user 2026-11-12)
        ext = "min"
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation=orientation, ticks=scale.edges,
                      spacing="uniform", extend=ext)
    n = len(scale.edges) - 1
    every = next((e for e in range(1, n + 1) if n % e == 0 and n // e + 1 <= max_labels), None)
    if every is None:
        every = int(math.ceil(n / (max_labels - 1)))
    lab_idx = set(range(0, n + 1, every))
    if len(lab_idx) < 3:          # e.g. a prime class count: evenly spaced labels and the last boundary
        every = next(e for e in range(1, n + 1) if len(range(0, n + 1, e)) + 1 <= max_labels)
        lab_idx = set(range(0, n + 1, every))
        if n - max(lab_idx) >= 0.6 * every:
            lab_idx.add(n)
        else:
            lab_idx.discard(max(lab_idx))
            lab_idx.add(n)
    labels = ["$10^{%g}$" % e if i in lab_idx else "" for i, e in enumerate(scale.edges)]
    (cb.ax.set_xticklabels if orientation == "horizontal" else cb.ax.set_yticklabels)(labels)
    cb.ax.tick_params(labelsize=6.3, length=2.5)
    cb.set_label(label, fontsize=6.8)
    return cb


def under_axes(fig, ax, gap_in, h_in, width_frac=1.0):
    """Axes just below an axes' drawn box (after its aspect is applied): gap and height in inches."""
    ax.apply_aspect()
    pos = ax.get_position()
    W, H = fig.get_size_inches()
    w = pos.width * width_frac
    return fig.add_axes([pos.x0 + 0.5 * (pos.width - w), pos.y0 - (gap_in + h_in) / H, w, h_in / H])


def latlon_track_templates(ax, gt, late, glon, glat, ok_seg, impact, c, zl, za, annotate):
    """The track / impact / start / zoom-box drawing of cops_cuda.figures.templates_figure (47fccd5), copied."""
    import matplotlib.patheffects as pe
    from matplotlib.collections import LineCollection
    halo = [pe.Stroke(linewidth=1.5, foreground="white"), pe.Normal()]
    halo_thick = [pe.Stroke(linewidth=2.9, foreground="white"), pe.Normal()]
    pts = np.column_stack((glon, glat))
    seg = np.stack((pts[:-1], pts[1:]), axis=1)
    keep = gt["outside"][:-1] & gt["outside"][1:] & ok_seg
    ax.add_collection(LineCollection(seg[keep], colors="k", linewidths=0.8, zorder=5, path_effects=halo))
    lk = keep & late[:-1] & late[1:]
    ax.add_collection(LineCollection(seg[lk], colors="k", linewidths=2.0, zorder=6, path_effects=halo_thick))
    ax.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=8 if c == 0 else 11, mfc="#F0E442", mec="k",
            mew=0.6, zorder=9)
    if c == 0:
        ax.plot([glon[0]], [glat[0]], marker="s", ms=3.6, mfc="white", mec="k", mew=0.6, zorder=9)
        ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k", lw=0.6, zorder=7,
                path_effects=halo)
    else:
        FIG._minute_marks(ax, gt, late, lambda p: FIG.lonlat(p), label=False)
        if annotate:
            ax.annotate("impact %s UTC" % impact["utc"].split("T")[1], xy=(impact["lon_deg"], impact["lat_deg"]),
                        xytext=(18, -48), textcoords="offset points", fontsize=6.0, ha="center",
                        arrowprops={"arrowstyle": "-|>", "lw": 0.6, "color": "k", "shrinkB": 5},
                        bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.4, "boxstyle": "round,pad=0.2"},
                        zorder=10)


def map_axes_setup(ax, c, lon_r, lat_r, xlabel=True, ylabel=True):
    ax.set_xlim(*lon_r)
    ax.set_ylim(*lat_r)
    if c == 0:
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-90, 91, 45))
    if ylabel:
        ax.set_ylabel("latitude [°]")
    if xlabel:
        ax.set_xlabel("longitude [°]")


def surface_handles(qshort):
    from matplotlib.patches import Patch
    return [Patch(facecolor=ZERO_SEEN * 0.92, edgecolor="none", label="seen during the fit, %s = 0" % qshort),
            Patch(facecolor=NOT_SEEN * 0.85, edgecolor="none", label="not seen during the fit")]


def track_handles_templates(gt):
    from matplotlib.lines import Line2D
    table, window, impact = gt["table"], gt["window"], gt["impact"]
    return [Line2D([], [], color="k", lw=0.8, label="ground track (sub-spacecraft intercept)"),
            Line2D([], [], color="k", lw=2.0, label="ground track in the late-visibility window"),
            Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", mew=0.5,
                   label="1-min marks from %s UTC" % window["start"].split("T")[1]),
            Line2D([], [], marker="s", ls="none", ms=3.6, mfc="white", mec="k", mew=0.6,
                   label="track start %s UTC, %.1f km" % (table["time"][0].strftime("%H:%M:%S"),
                                                         table["altitude"][0] / 1000)),
            Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                   label="Rosetta impact, lon %.2f°, lat %.2f°" % (impact["lon_deg"], impact["lat_deg"]))]


def cliff_handles(meta):
    from matplotlib.patches import Patch
    h = []
    for name, col, thr in (("cliff_ge_45deg", FIG.CLIFF45_COLOR, "45–60°"), ("cliff_ge_60deg", FIG.CLIFF60_COLOR, "≥ 60°")):
        frac = meta.get("templates", {}).get(name, {}).get("area_fraction", np.nan)
        h.append(Patch(facecolor=col, label="slope %s (%s: %.1f%% of area)" % (thr, name, 100 * frac)))
    h[0].set_label(h[0].get_label().replace("of area)", "of area, incl. ≥ 60°)"))
    return h


# ============================================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for a in ("--base", "--track-base", "--bundle", "--ft-dir", "--fit-dir", "--maps-dir", "--geometry", "--maps-code",
              "--out-dir"):
        ap.add_argument(a, required=True)
    ap.add_argument("--quantity", required=True, choices=("activity", "total", "h2o"),
                    help="activity: A_final; total: F_H2O + F_CO2; h2o: F_H2O = q_H2O A (user decision 2026-11-12: "
                         "'plot H2O only for the rate')")
    ap.add_argument("--time", default=None, help="HH:MM:SS (UTC, 2016-09-30) of the total production rate: one of the "
                                                 "map times of the maps job's geometry (times.txt)")
    ap.add_argument("--classes", default="15,20", help="MIN,MAX number of colour classes per panel")
    ap.add_argument("--template-relief", type=float, default=0.0,
                    help="strength of the relief modulation of the class colours in the TEMPLATE figure and its value "
                         "panels (0: flat colours; user request 2026-11-12 for the production rate); ground tracks keep "
                         "flat colours")
    ap.add_argument("--low-cut-area", type=float, default=0.0,
                    help="area fraction cut at the bottom of each panel's class range (0: entire range; 0.001: the "
                         "range holding 99.9 %% of the seen non-zero area, user decision 2026-11-12 for the activity)")
    ap.add_argument("--gravity-cache", default=None)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--pixels3d", type=int, default=800)
    ap.add_argument("--model-name", default=M.MODEL_NAME,
                    help="name of the fitted state in the captions (default: %(default)s; 2026-11-13: set per fit)")
    args = ap.parse_args()
    M.MODEL_NAME = args.model_name
    global NMIN, NMAX
    NMIN, NMAX = (int(v) for v in args.classes.split(","))
    out = Path(args.out_dir)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    q = args.quantity
    if q in ("total", "h2o"):
        if not args.time:
            raise SystemExit("--quantity %s needs --time HH:MM:SS" % q)
        utc = "2016-09-30T%s" % args.time
        labs = [ln.split()[0] for ln in open(Path(args.geometry) / "times.txt")
                if not ln.startswith("#") and ln.split() and ln.split()[0].startswith("t") and ln.split()[1] == utc]
        if len(labs) != 1:
            raise SystemExit("%s is not a map time of %s/times.txt" % (utc, args.geometry))
        M.TIMES = ((labs[0], utc),)              # load_fields computes the flux at this time only
        q = "%s_%s" % (q, args.time.replace(":", ""))
    elif args.time:
        raise SystemExit("--time applies to --quantity total / h2o")
    sub = "activity" if q == "activity" else "production_rate"
    stem = ("activity" if q == "activity" else
            "production_rate_%sUTC" % q.split("_")[1] if q.startswith("total") else
            "h2o_production_rate_%sUTC" % q.split("_")[1])
    stage = out / ("v2_%s.incomplete_%s" % (q, stamp))
    finals = [out / sub / ("12M_templates_%s.png" % stem), out / "checks" / ("%s_class_tables.txt" % q)]
    for f in finals:
        if f.exists():
            raise SystemExit("%s exists (never overwritten; move it aside first)" % f)
    stage.mkdir(parents=True)
    say("stage %s, quantity %s" % (stage, q))

    # ---- track (observation-only base) and fields (fit base), as pub_maps_figures.main
    base_tab = FIG.read_table(args.track_base)
    window = FIG.late_window(base_tab, args.bundle, None, None)
    gt = FIG.groundtrack_data(args.track_base, None, args.bundle, window, args.device)
    say("track: %d rows, impact %s UTC; window %d rows" % (len(gt["table"]["utc"]), gt["impact"]["utc"],
                                                           len(window["points"])))
    fields = M.load_fields(args)
    if not fields["checks"]["passed"]:
        raise SystemExit("field checks failed")
    # bounds of the slope multipliers from the fit (data-driven caption note, 2026-11-13: refit with multipliers >= 1)
    with np.load(Path(args.fit_dir) / "fit_restart.npz", allow_pickle=False) as z:
        rnames = [str(n) for n in json.loads(str(z["meta_json"]))["names"]]
        rx, rlo = np.asarray(z["x"], np.float64), np.asarray(z["lo"], np.float64)
    mb = [(n[2:], rx[i], rlo[i]) for i, n in enumerate(rnames) if n.startswith("c_slope_")]
    lows = sorted({lo for _, _, lo in mb})
    at = [n for n, x, lo in mb if x - lo <= 1e-9 * max(1.0, abs(lo))]
    if mb and lows != [0.0]:
        M.MULT_NOTE = ("; the %d free multipliers bounded below by %s, %s at that bound" % (
            len(mb), " / ".join("%g" % v for v in lows),
            ("none" if not at else ", ".join(n.replace("slope_", "").replace("_", "-") + " deg" for n in at))))
    fields["checks"]["multiplier_bounds"] = {"lower": {n: lo for n, _, lo in mb}, "at_lower_bound": at,
                                            "values": {n: x for n, x, _ in mb}}
    say("multipliers: %s; lower bounds %s; at the bound: %s" % ({n: "%.6g" % x for n, x, _ in mb}, lows, at))
    if not np.array_equal(fields["normal"], gt["normals"]):
        raise SystemExit("normals differ between the track mesh and the fields")
    nm = M.load_module(args.maps_code, "nonneg_maps_frozen_v2")
    e1 = nm.read_elements(Path(args.base) / "elements.bin")
    e2 = nm.read_elements(Path(args.track_base) / "elements.bin")
    same = {k: bool(np.array_equal(e1[k], e2[k])) for k in ("centroid", "normal", "area", "ACO2")}
    del e1, e2
    if not all(same.values()):
        raise SystemExit("track base mesh differs: %s" % same)
    fields["checks"]["track_base_elements_equal_fit_base"] = same
    say("field checks passed; track base mesh equal %s" % same)
    seen, area = fields["seen"], fields["area"]
    if q == "activity":
        values = fields["A"]
        qshort, qlabel, unit = "A", "H₂O activity factor A", ""
        cb_label = "A (seen facets; log classes)"
    else:
        lab = M.TIMES[0][0]
        tinfo = fields["times"][lab]
        if q.startswith("total"):
            values = fields["times"][lab]["total"]
            qshort, qlabel, unit = "F", "total production rate F (H₂O + CO₂)", " m⁻² s⁻¹"
            cb_label = "F [m⁻² s⁻¹] (seen facets; log classes)"
        else:
            values = fields["times"][lab]["F_H2O"]
            qshort, qlabel, unit = "F_H2O", "H₂O production rate", " m⁻² s⁻¹"
            cb_label = "H₂O production rate [m⁻² s⁻¹] (seen facets; log classes)"
    meta = FIG.bundle_metadata(args.bundle)
    slope, cache_path = FIG.gravity_slope(args.bundle, args.gravity_cache, len(gt["faces"]))
    ctx = {"gt": gt, "fields": fields, "meta": meta, "cache_path": cache_path, "args": args,
           "stats": {"A": M.value_stats(fields["A"], seen, area)},
           "track_base": str(Path(args.track_base).resolve()), "fit_base": str(Path(args.base).resolve()),
           "track_offset": M.gauge_offset(args.track_base), "fit_offset": M.gauge_offset(args.base)}
    plt = FIG.style()
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    I_ll = M.lonlat_intensity(V, F, N)
    table, impact = gt["table"], gt["impact"]
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = FIG.lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    zl, za = M.zoom_box(gt, late)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    wins = {"global": ((-180, 180), (-90, 90), GLOBAL_W), "zoom": (zl, za, ZOOM_W)}

    # ---- lon/lat windows: face images, displayed facets, scales, images
    tables, ttext, scales, imgs = {}, [], {}, {}
    from matplotlib.colors import Normalize
    snorm = Normalize(0.0, 90.0)
    magma = plt.get_cmap("magma")
    slope_rgb = magma(snorm(np.clip(slope, 0, 90)))[:, :3]
    cliff_rgb = FIG.face_colors(len(F), None)
    masks = gt["masks"] or {}
    cliff_rgb[masks.get("cliff_ge_45deg", slope >= 45.0)] = FIG.CLIFF45_COLOR
    cliff_rgb[masks.get("cliff_ge_60deg", slope >= 60.0)] = FIG.CLIFF60_COLOR
    cliff_rgb = cliff_rgb * I_ll[:, None]
    for wname, (lon_r, lat_r, width) in wins.items():
        face, ss = lonlat_faces(V, F, lon_r, lat_r, width, args.device)
        disp = displayed(face)
        sel = disp[seen[disp] & (values[disp] > 0)]
        sc = Scale("%s/%s" % (q, wname), values[sel], area[sel], args.low_cut_area)
        scales[wname] = sc
        tab, txt = class_table("%s, lon/lat %s window" % (q, wname), values, seen, area, disp, sc)
        if args.template_relief > 0:
            rc = relief_check(sc, args.template_relief)
            tab["template_relief_check"] = rc
            txt += ("\n  template relief modulation (strength %g): CIELAB deltaE between adjacent classes at the same "
                    "illumination, min / median: darkest slope %.1f / %.1f, level ground %.1f / %.1f, brightest slope "
                    "%.1f / %.1f" % (rc["strength"], rc["darkest_slope"]["min_deltaE_adjacent"],
                                     rc["darkest_slope"]["median_deltaE_adjacent"],
                                     rc["level_ground"]["min_deltaE_adjacent"], rc["level_ground"]["median_deltaE_adjacent"],
                                     rc["brightest_slope"]["min_deltaE_adjacent"],
                                     rc["brightest_slope"]["median_deltaE_adjacent"]))
        tables["latlon_" + wname] = tab
        ttext.append(txt)
        imgs[("value", wname)] = shade(face, facet_rgb(values, seen, sc, I_ll), ss)
        imgs[("value_tpl", wname)] = (shade(face, facet_rgb(values, seen, sc, I_ll, args.template_relief), ss)
                                      if args.template_relief > 0 else imgs[("value", wname)])
        imgs[("slope", wname)] = shade(face, slope_rgb, ss)
        imgs[("cliff", wname)] = shade(face, cliff_rgb, ss)
        del face
        say("window %s: %d displayed facets, %d seen non-zero -> %s" % (wname, disp.size, sel.size, sc.text()))

    def scale_caption(sc, where):
        return ("%s: the seen non-zero values drawn in this panel range from %.3e to %.3e%s; %s." % (
            where, sc.vmin, sc.vmax, unit.replace("⁻²", "^-2").replace("⁻¹", "^-1") if unit else "",
            sc.text().replace("$10^{", "10^").replace("}$", "")) + sc.below_text() +
            (" The smallest value is a denormal (subnormal) float64 number, i.e. an underflow of the model's "
             "softplus / multiplier product below 2.2e-308, not a physical value." if sc.denormal else ""))

    surface_def = (
        "Colour only on the facets seen during the fit (the fit's nonneg_facets.npz 'seen': the facet contributes to at "
        "least one fitted NG or RG observation; %s of the %s facets, %s of the area); facets not seen during the fit "
        "are grey with relief shading. Seen facets with %s = 0 exactly are light green (relief-shaded). Seen facets with "
        "%s > 0: discrete log10 classes (colours 'plasma', drawn without relief shading so that the colour equals the "
        "colour bar's); the classes cover %s with 15-20 intervals at round boundaries (checks/%s_class_tables.txt "
        "lists every class with its facet count, area and value range)." % (
            M.thousands(seen.sum()), M.thousands(len(seen)), M.pct(np.sum(area[seen]) / np.sum(area)), qshort, qshort,
            ("the range of the seen non-zero values drawn in each panel (scanned per panel) that holds %g %% of their "
             "area (from their area-weighted %g %% quantile to their maximum; smaller values have the colour of the "
             "lowest class; user decision 2026-11-12)" % (100 - 100 * args.low_cut_area, 100 * args.low_cut_area))
            if args.low_cut_area > 0 else
            "the entire range of the seen non-zero values drawn in each panel (scanned per panel)", q))
    surface_def_tpl = surface_def
    if args.template_relief > 0:
        surface_def_tpl = surface_def.replace(
            "drawn without relief shading so that the colour equals the colour bar's",
            "colours modulated by the relief shading (strength %g: darkened on slopes facing away from the light, "
            "lightened toward white on slopes facing it, unchanged on level ground lit from 40 deg above the horizon; "
            "hue kept; adjacent classes stay distinguishable under the same illumination, checks/%s_class_tables.txt); "
            "the colour bar shows the unshaded class colours" % (args.template_relief, q))
        assert surface_def_tpl != surface_def
    if q == "activity":
        what = ("H2O activity factor A of the surface in %s. %s" % (args.model_name, M.model_text(ctx)))
    elif q.startswith("total"):
        what = ("Total production rate F = F_H2O + F_CO2 (instantaneous emitted molecular flux per facet, m^-2 s^-1) of "
                "%s, at %s UTC on 2016-09-30. "
                "%s %s" % (args.model_name, tinfo["utc"].split("T")[1], M.time_text(ctx, lab), M.model_text(ctx)))
    else:
        n_ = tinfo["facets"]
        what = ("H2O production rate F_H2O (instantaneous emitted H2O molecular flux per facet, m^-2 s^-1) of "
                "%s, at %s UTC on 2016-09-30. "
                "Time %s UTC (r_h = %.6f AU). F_H2O = q_H2O A as the model computes it: q_H2O = q_min on the night side "
                "(cos SZA <= 0) and in cast shadow (Sun ray from the facet blocked by the nucleus), else q_min + (q_max - "
                "q_min) cos SZA, with q_max = 600 / r_h^%.7g = %.4f and q_min = 0.2 q_max (r_h in AU); facets at this "
                "time: %s day, %s night, %s in cast shadow; total H2O production Q_H2O = sum F_H2O area = %.4e s^-1. "
                "F_H2O = 0 exactly where A = 0. %s" % (
                    args.model_name, tinfo["utc"].split("T")[1], tinfo["utc"].split("T")[1], tinfo["rh_AU"], tinfo["p"],
                    tinfo["q"]["H2O"][1], M.thousands(n_["day"]), M.thousands(n_["night"]),
                    M.thousands(n_["cast_shadow"]), tinfo["Q_H2O"], M.model_text(ctx)))
    trk = M.track_text(ctx)
    slp = M.slope_text(ctx)
    surf0 = what.split(". ")[0].replace("Total production rate", "total production rate", 1)
    mesh = ("on the %s shape model (%s facets), longitude/latitude, outermost surface where overhangs put several "
            "surface layers at one longitude/latitude" % (Path(gt["shape"]).name, M.thousands(len(F))))
    written, captions = [], []

    def save(fig, prefix, caption):
        written.extend(FIG.save(fig, prefix))
        plt.close(fig)
        written.append(FIG.write_caption(prefix, caption))
        captions.append((str(Path(prefix).relative_to(stage)), caption))

    title_row = {"value": "%s, %s" % (qlabel, "final state" if q == "activity" else tinfo["utc"].split("T")[1] + " UTC"),
                 "slope": "effective gravitational slope", "cliff": "cliff templates"}

    # ---- templates figure (copy of cops_cuda.figures.templates_figure layout; changed: value row with per-panel
    #      discrete colour bars under the panels, no header text)
    def draw_map(ax, row, wname, c, ylabel):
        lon_r, lat_r, _ = wins[wname]
        img = imgs[("value_tpl", wname)] if row == "value" else imgs[(row, wname)]   # templates: relief option
        ax.imshow(img, extent=(lon_r[0], lon_r[1], lat_r[0], lat_r[1]), interpolation="lanczos",
                  aspect="equal" if c == 0 else 1.0 / coslat)
        latlon_track_templates(ax, gt, late, glon, glat, ok_seg, impact, c, zl, za, annotate=(row == "value"))
        map_axes_setup(ax, c, lon_r, lat_r, xlabel=False, ylabel=ylabel)

    row_h = 2.05
    hsp = 0.30
    rows_in = 3 * row_h + 0.05
    height = (rows_in + 3 * hsp * rows_in / 4 + 1.35) / 0.995      # rows + gaps + legend; panels keep row_h
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, height))
    gs = fig.add_gridspec(4, 3, width_ratios=[2.0, 1.0, 0.05], height_ratios=[row_h, 0.05, row_h, row_h],
                          wspace=0.10, hspace=hsp, left=0.07, right=0.93, top=0.995, bottom=1.35 / height)
    letters = iter("abcdef")
    axes = {}
    for r, row in ((0, "value"), (2, "slope"), (3, "cliff")):
        for c, wname in enumerate(("global", "zoom")):
            ax = fig.add_subplot(gs[r, c])
            draw_map(ax, row, wname, c, ylabel=(c == 0))
            if c == 0:
                ax.set_title(title_row[row], fontsize=8, loc="left", pad=3)
            if r == 3:
                ax.set_xlabel("longitude [°]")
            ax.text(0.015, 0.96, "(%s)" % next(letters), transform=ax.transAxes, va="top", fontsize=7.5,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8})
            axes[(row, wname)] = ax
        cax = fig.add_subplot(gs[r, 2])
        if row == "slope":
            from matplotlib.cm import ScalarMappable
            cb = fig.colorbar(ScalarMappable(norm=snorm, cmap=magma), cax=cax, extend="max")
            cb.set_label("slope [°] (> 90: overhang)", fontsize=6.5)
            cb.set_ticks([0, 15, 30, 45, 60, 75, 90])
        else:
            cax.axis("off")
    for wname, ml in (("global", 7), ("zoom", 5)):
        discrete_colorbar(fig, under_axes(fig, axes[("value", wname)], 0.33, 0.11), scales[wname],
                          cb_label if wname == "global" else cb_label.split(" (")[0] + " (impact area)", max_labels=ml)
    fig.legend(handles=surface_handles(qshort) + cliff_handles(meta) + track_handles_templates(gt), loc="lower left",
               bbox_to_anchor=(0.07, 0.0), ncol=2, fontsize=6.3, frameon=False, handlelength=2.0, columnspacing=1.2)
    cap_t = ("%s, with the effective gravitational slope and the cliff templates, %s. Left: global; right: the impact "
             "area (box in the left panels). (a, b) %s %s %s %s (c, d) Effective gravitational slope and (e, f) cliff "
             "templates (as in the late-visibility figure). %s %s" % (
                 what.split(". ")[0], mesh, surface_def_tpl, scale_caption(scales["global"], "(a)"),
                 scale_caption(scales["zoom"], "(b)"), what[len(what.split(". ")[0]) + 2:], slp, trk))
    save(fig, stage / sub / ("12M_templates_%s" % stem), cap_t)
    say("wrote templates figure")

    # ---- each panel of the templates figure as its own figure
    vname = "activity" if q == "activity" else "production_rate" if q.startswith("total") else "h2o_production_rate"
    pan_names = {("value", "global"): "a_%s_global" % vname, ("value", "zoom"): "b_%s_impact_area" % vname,
                 ("slope", "global"): "c_slope_global", ("slope", "zoom"): "d_slope_impact_area",
                 ("cliff", "global"): "e_cliff_global", ("cliff", "zoom"): "f_cliff_impact_area"}
    for (row, wname), pname in pan_names.items():
        c = 0 if wname == "global" else 1
        figsize = (FIG.DOUBLE_COLUMN_IN, 4.6) if c == 0 else (4.4, 5.2)
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.10, 0.36, 0.86, 0.60] if c == 0 else [0.15, 0.36, 0.80, 0.60])
        draw_map(ax, row, wname, c, ylabel=True)
        ax.set_xlabel("longitude [°]")
        ax.set_title(title_row[row] + ("" if c == 0 else ", impact area"), fontsize=8, loc="left", pad=3)
        if row == "value":
            discrete_colorbar(fig, under_axes(fig, ax, 0.42, 0.12, 0.9), scales[wname], cb_label,
                              max_labels=7 if c == 0 else 5)
            hs = surface_handles(qshort) + track_handles_templates(gt)
        elif row == "slope":
            from matplotlib.cm import ScalarMappable
            cb = fig.colorbar(ScalarMappable(norm=snorm, cmap=magma), cax=under_axes(fig, ax, 0.42, 0.12, 0.9),
                              orientation="horizontal", extend="max")
            cb.set_label("slope [°] (> 90: overhang)", fontsize=6.8)
            cb.set_ticks([0, 15, 30, 45, 60, 75, 90])
            hs = track_handles_templates(gt)
        else:
            hs = cliff_handles(meta) + track_handles_templates(gt)
        fig.legend(handles=hs, loc="lower left", bbox_to_anchor=(0.02, 0.0), ncol=2, fontsize=6.0, frameon=False,
                   handlelength=1.8, columnspacing=1.0)
        if row == "value":
            pc = "%s %s, %s. %s %s %s %s" % (
                what.split(". ")[0], "global map" if c == 0 else "the impact area", mesh, surface_def_tpl,
                scale_caption(scales[wname], "Colour classes"), what[len(what.split(". ")[0]) + 2:], trk)
        elif row == "slope":
            pc = "Effective gravitational slope, %s, %s. %s %s" % ("global map" if c == 0 else "the impact area", mesh,
                                                                   slp.split(" Cliff templates")[0], trk)
        else:
            pc = "Cliff templates, %s, %s. %s %s" % ("global map" if c == 0 else "the impact area", mesh, slp, trk)
        save(fig, stage / sub / "panels" / ("12M_templates_%s_%s" % (stem, pname)), pc)
    say("wrote 6 panel figures")

    # ---- ground track, lon/lat (copy of cops_cuda.figures.groundtrack_figures layout; changed: surface colours,
    #      per-panel discrete colour bars under the maps, legend surface entries)
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    alt = np.maximum(table["altitude"], 0.1)
    anorm = LogNorm(vmin=max(alt.min(), 0.1), vmax=alt.max())
    vir = plt.get_cmap("viridis")
    no_jump = step <= max(150.0, 5.0 * np.median(step))
    no_seam = (np.abs(np.diff(glon)) <= 180.0) & no_jump
    legend_extra = [Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                           label="Rosetta impact (lon %.2f°, lat %.2f°)" % (impact["lon_deg"], impact["lat_deg"])),
                    Line2D([], [], marker="s", ls="none", ms=4, mfc="white", mec="k", mew=0.6,
                           label="track start, %s UTC, %.1f km" % (table["time"][0].strftime("%H:%M:%S"),
                                                                    table["altitude"][0] / 1000.0))]

    def gt_legend(fig, anchor, hidden):
        hs = surface_handles(qshort)
        hs.append(Line2D([], [], color="#21918C", lw=2.6, label="ground track, late-visibility window"))
        hs.append(Line2D([], [], color="#21918C", lw=1.1, label="ground track (sub-spacecraft intercept)"))
        if hidden:
            hs.append(Line2D([], [], color="0.45", lw=0.6, ls=(0, (1, 1.5)), label="ground track, hidden"))
        hs.append(Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", mew=0.5,
                         label="1-min marks from window start"))
        hs.extend(legend_extra)
        fig.legend(handles=hs, loc="lower left", bbox_to_anchor=anchor, fontsize=6.3, ncol=2, frameon=False,
                   handlelength=2.2, columnspacing=1.0)

    Hll = 9.6
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, Hll))
    ax = fig.add_axes([0.08, 0.665, 0.90, 0.315])
    ax.imshow(imgs[("value", "global")], extent=(-180, 180, -90, 90), interpolation="lanczos", aspect="equal")
    FIG._colored_track(ax, glon, glat, alt, anorm, vir, 2.6, mask=late & gt["outside"], seg_mask=no_seam, zorder=5)
    lc = FIG._colored_track(ax, glon, glat, alt, anorm, vir, 1.0, mask=gt["outside"], seg_mask=no_seam, zorder=6)
    ax.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=9, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
    ax.plot([glon[0]], [glat[0]], marker="s", ms=4, mfc="white", mec="k", mew=0.6, zorder=9)
    ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k", lw=0.6, zorder=7)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-90, 91, 30))
    ax.set_xlabel("longitude [°]")
    ax.set_ylabel("latitude [°]")
    ax.text(0.01, 0.97, "(a)", transform=ax.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    discrete_colorbar(fig, under_axes(fig, ax, 0.48, 0.11, 0.8), scales["global"], cb_label, max_labels=7)
    axz = fig.add_axes([0.08, 0.17, 0.90, 0.31])
    axz.imshow(imgs[("value", "zoom")], extent=(zl[0], zl[1], za[0], za[1]), interpolation="lanczos",
               aspect=1.0 / coslat)
    FIG._colored_track(axz, glon, glat, alt, anorm, vir, 2.6, mask=late, seg_mask=no_seam, zorder=5)
    FIG._colored_track(axz, glon, glat, alt, anorm, vir, 1.1, mask=gt["outside"], seg_mask=no_seam, zorder=6)
    FIG._minute_marks(axz, gt, late, lambda p: FIG.lonlat(p), label=False)
    axz.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=11, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
    axz.annotate("Rosetta impact\n%s UTC\nlon %.2f°, lat %.2f°" % (
        impact["utc"].split("T")[1], impact["lon_deg"], impact["lat_deg"]),
        xy=(impact["lon_deg"], impact["lat_deg"]), xytext=(-75, 45), textcoords="offset points",
        fontsize=6.3, ha="center", arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": "k", "shrinkB": 6},
        bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.5, "boxstyle": "round,pad=0.25"}, zorder=10)
    axz.set_xlim(*zl)
    axz.set_ylim(*za)
    axz.set_xlabel("longitude [°]")
    axz.set_ylabel("latitude [°]")
    axz.text(0.01, 0.97, "(b)", transform=axz.transAxes, va="top", fontsize=8,
             bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    discrete_colorbar(fig, under_axes(fig, axz, 0.48, 0.11, 1.0), scales["zoom"],
                      cb_label.split(" (")[0] + " (impact area)", max_labels=5)
    gt_legend(fig, (0.08, 0.0), hidden=False)
    cb = fig.colorbar(lc, cax=fig.add_axes([0.70, 0.03, 0.28, 0.011]), orientation="horizontal")
    cb.set_label("spacecraft altitude [m]", fontsize=7)
    cap_ll = ("End-of-mission ground track of Rosetta and its impact point on the %s shape model (%s facets), "
              "longitude/latitude (equirectangular), outermost surface where overhangs put several surface layers at "
              "one longitude/latitude. (a) Global map; (b) the impact area (box in a). Surface: %s. %s %s %s %s Track: "
              "coloured by altitude (log colour bar 'spacecraft altitude'), thick within the late-visibility window. "
              "%s" % (Path(gt["shape"]).name, M.thousands(len(F)), surf0, surface_def,
                      scale_caption(scales["global"], "(a)"), scale_caption(scales["zoom"], "(b)"),
                      what[len(what.split(". ")[0]) + 2:], trk))
    save(fig, stage / sub / ("12M_groundtrack_%s_latlon" % stem), cap_ll)
    say("wrote ground-track lon/lat figure")

    # ---- ground track, 3-D views (copy of cops_cuda.figures.groundtrack_figures; changed: surface colours from the
    #      face images, per-panel discrete colour bars under the panels, legend surface entries)
    xyz = np.column_stack((table["x"], table["y"], table["z"]))
    ip = np.asarray(impact["point_m"])
    ground = gt["ground"]
    view = FIG.OrthoView(xyz[late][0] - ip)
    I3, front = M.ortho_intensity(N, view)
    body_half = 1.04 * np.max(np.abs(np.column_stack(view.project(V - V.mean(axis=0))[:2])))
    zoom_half = 1.15 * max(np.max(np.linalg.norm(ground[late] - ip, axis=1)), 1.0)
    for name in ("late_ng_not_rg", "late_rg_not_ng"):
        if name in masks and masks[name].any():
            cen = V[F[masks[name]]].mean(axis=1)
            zoom_half = max(zoom_half, 1.08 * np.max(np.abs(np.column_stack(view.project(cen - ip)[:2]))))
    close_half = 1.3 * max(np.max(np.abs(np.column_stack(view.project(ground[late] - ip)[:2]))), 20.0)
    panels = ((V.mean(axis=0), body_half), (ip, zoom_half), (ip, close_half))
    pixels = args.pixels3d
    H3 = 4.05
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, H3))
    width = 0.315
    hpan = width * FIG.DOUBLE_COLUMN_IN / H3
    axs = [fig.add_axes([0.005 + k * (width + 0.0175), 1.0 - hpan - 0.005, width, hpan]) for k in range(3)]
    lc = None
    dummy = np.tile(NOT_SEEN, (len(F), 1))
    sc3 = []
    for k, (ax, (center, half)) in enumerate(zip(axs, panels)):
        _, to_image, visible, face, ss = FIG.render_ortho(V, F, N, view, center, half, pixels, dummy, args.device,
                                                          return_faces=True)
        disp = displayed(face)
        sel = disp[seen[disp] & (values[disp] > 0)]
        sc = Scale("%s/3d_%s" % (q, "abc"[k]), values[sel], area[sel], args.low_cut_area)
        sc3.append(sc)
        tab, txt = class_table("%s, 3-D panel (%s)" % (q, "abc"[k]), values, seen, area, disp, sc)
        tables["3d_" + "abc"[k]] = tab
        ttext.append(txt)
        image = shade(face, facet_rgb(values, seen, sc, I3), ss)
        del face
        ax.imshow(image, interpolation="lanczos", origin="upper", extent=(0, pixels, pixels, 0))
        gx, gy, _ = to_image(ground)
        tol = max(3.0, 4.0 * 2.0 * half / pixels)
        vis = visible(ground, tol) & gt["outside"]
        hidden = ~vis & gt["outside"]
        FIG._colored_track(ax, gx, gy, alt, anorm, vir, 0.6, mask=hidden, seg_mask=no_jump, zorder=4,
                           solid_color="0.45", dotted=True)
        FIG._colored_track(ax, gx, gy, alt, anorm, vir, 2.6, mask=vis & late, seg_mask=no_jump, zorder=5)
        lc = FIG._colored_track(ax, gx, gy, alt, anorm, vir, 1.1, mask=vis, seg_mask=no_jump, zorder=6)
        if vis[0]:
            ax.plot([gx[0]], [gy[0]], marker="s", ms=4, mfc="white", mec="k", mew=0.6, zorder=9)
        ix, iy, _ = to_image(ip[None, :])
        ax.plot([ix[0]], [iy[0]], marker="*", ms=9 if k < 2 else 11, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
        if k < 2:
            s = panels[k + 1][1] / half * pixels / 2
            ax.plot([ix[0] - s, ix[0] + s, ix[0] + s, ix[0] - s, ix[0] - s],
                    [iy[0] - s, iy[0] - s, iy[0] + s, iy[0] + s, iy[0] - s], color="k", lw=0.6, zorder=7)
        if k == 2:
            FIG._minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=True)
            ax.annotate("Rosetta impact\n%s UTC" % impact["utc"].split("T")[1], xy=(ix[0], iy[0]),
                        xytext=(-45, -55), textcoords="offset points", fontsize=6.3, ha="center",
                        arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": "k", "shrinkB": 6},
                        bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.5, "boxstyle": "round,pad=0.25"},
                        zorder=10)
        else:
            FIG._minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=False)
        ax.set_xlim(0, pixels)
        ax.set_ylim(pixels, 0)
        ax.set_xticks([])
        ax.set_yticks([])
        FIG._scale_bar(ax, half, pixels)
        ax.text(0.03, 0.97, "(%s)" % "abc"[k], transform=ax.transAxes, va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
        discrete_colorbar(fig, under_axes(fig, ax, 0.08, 0.10, 0.78), sc,
                          "%s, panel (%s)" % (cb_label.split(" (")[0], "abc"[k]), max_labels=4)
        say("3-D panel %s: %d displayed facets -> %s" % ("abc"[k], disp.size, sc.text()))
    gt_legend(fig, (0.005, 0.0), hidden=True)
    cb = fig.colorbar(lc, cax=fig.add_axes([0.70, 0.07, 0.29, 0.03]), orientation="horizontal")
    cb.set_label("spacecraft altitude [m]", fontsize=7)
    cap_3d = ("End-of-mission ground track of Rosetta and its impact point on the %s shape model (%s facets), "
              "orthographic views rendered on the GPU (z-buffer). (a) The whole nucleus seen down the final approach "
              "(view direction: from the start of the late-visibility window toward the impact point); (b) the impact "
              "region (box in a); (c) close-up of the late descent (box in b), with labelled 1-minute marks (UTC, "
              "altitude). Scale bars in each panel. Surface: %s. %s %s %s %s %s Track: coloured by altitude (log colour "
              "bar 'spacecraft altitude'), thick within the late-visibility window, dotted grey where hidden behind the "
              "nucleus. %s" % (Path(gt["shape"]).name, M.thousands(len(F)), surf0, surface_def,
                               scale_caption(sc3[0], "(a)"), scale_caption(sc3[1], "(b)"), scale_caption(sc3[2], "(c)"),
                               what[len(what.split(". ")[0]) + 2:], trk))
    save(fig, stage / sub / ("12M_groundtrack_%s_3d" % stem), cap_3d)
    say("wrote ground-track 3-D figure")

    # ---- per slope bin: area, seen area, multiplier, share of the seen production (input of the cliff figure)
    bins = fields["bins"]
    edges = fields["edges"]
    binrows = []
    for b in range(len(edges) - 1):
        m = bins == b
        row = {"bin": b, "interval_deg": [edges[b], edges[b + 1]], "multiplier": float(fields["mult"][b]),
               "facets": int(m.sum()), "area_fraction": float(area[m].sum() / area.sum()),
               "seen_area_fraction_of_bin": float(area[m & seen].sum() / area[m].sum()),
               "bin_share_of_seen_area": float(area[m & seen].sum() / area[seen].sum()),
               "activity_area_share_seen": float(fields["A"][m & seen].dot(area[m & seen]) /
                                                 fields["A"][seen].dot(area[seen]))}
        if q != "activity":
            Ft = values
            row["production_share_seen"] = float(Ft[m & seen].dot(area[m & seen]) / Ft[seen].dot(area[seen]))
            row["production_share_all"] = float(Ft[m].dot(area[m]) / Ft.dot(area))
        binrows.append(row)
    (stage / "checks").mkdir(exist_ok=True)
    json.dump({"quantity": q, "bins": binrows, "seen": "nonneg_facets.npz seen", "time_utc": None if q == "activity"
               else tinfo["utc"]}, open(stage / "checks" / ("slope_bins_%s.json" % q), "w"), indent=1)
    with open(stage / "checks" / ("%s_class_tables.txt" % q), "w") as f:
        f.write("# pub_maps_v2_figures.py (%s): colour classes of every panel showing %s; facets = the facets drawn in "
                "the panel (rasterizer face index per supersampled pixel); seen = the fit's nonneg_facets.npz seen\n\n"
                % (SCHEMA, q))
        f.write("\n\n".join(ttext) + "\n")
    json.dump(fields["checks"], open(stage / "checks" / ("%s_fields_verification.json" % q), "w"), indent=1)
    man = {"schema": SCHEMA, "created_utc": stamp, "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "command": sys.argv, "script_sha256": M.sha256_file(__file__), "helper_sha256": M.sha256_file(M.__file__),
           "figures_py_sha256": M.sha256_file(FIG.__file__), "quantity": q, "impact": impact,
           "window": {k: window[k] for k in ("start", "stop", "source", "altitude_min_m", "altitude_max_m")},
           "track_base": ctx["track_base"], "fit_base": ctx["fit_base"], "class_tables": tables,
           "colormap": CMAP_NAME, "zero_seen_colour": hexcol(ZERO_SEEN), "not_seen_colour": hexcol(NOT_SEEN),
           "template_relief": args.template_relief, "low_cut_area": args.low_cut_area,
           "files": [str(Path(p).relative_to(stage)) for p in written], "wall_seconds": time.time() - T0}
    json.dump(man, open(stage / "checks" / ("%s_manifest.json" % q), "w"), indent=1, default=float)
    with open(stage / "checks" / ("%s_captions.json" % q), "w") as f:
        json.dump([{"figure": n, "caption": c} for n, c in captions], f, indent=1, ensure_ascii=False)
    # ---- move into place (file by file; the folders are shared by the quantity jobs)
    for p in sorted(stage.rglob("*")):
        if p.is_file():
            dest = out / p.relative_to(stage)
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                raise SystemExit("%s exists" % dest)
            p.rename(dest)
    for d in sorted((p for p in stage.rglob("*") if p.is_dir()), key=lambda p: -len(p.parts)):
        d.rmdir()
    stage.rmdir()
    say("done: %d files" % (len(written) + 5))
    return 0


if __name__ == "__main__":
    sys.exit(main())

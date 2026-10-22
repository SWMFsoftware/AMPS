#!/usr/bin/env python3
"""3D views of Rosetta's final approach and descent where the NG pressure rose (publication format, .eps + .png).

    .venv/bin/python tools/descent3d_figures.py --base <base run> --bundle <restart template bundle> \\
        --visibility-bundle <visibility bundle> --out <prefix> [--impact <groundtrack impact.json>]

The whole nucleus (nothing cropped) seen broadside to the trajectory: the view looks along the normal of
the plane that best fits the spacecraft positions from ``--start`` to the last row, tilted ``--elevation``
upward about the local vertical of the impact site, from the side on which the spacecraft descends from
left to right. The trajectory passes across the view; parts hidden behind the nucleus are dotted. Every
figure pairs the whole nucleus with a zoom on the impact area in the same direction (the last
``--zoom-altitude`` metres of the descent), with altitude drop lines to the ground track (nearest surface
point, the run's altitude) in the zoom. The "NG jump" is the first NG measurement after the NG data
gap with the largest rise of the measured NG pressure (10:24:38 -> 10:33:10 UTC, x4.1, in the end-of-mission data;
``--jump`` sets it explicitly). Figures (``<prefix>_<name>``):

* ``terrain``: relief shading; trajectory coloured by the measured NG pressure (grey: no NG data);
* ``cos_sza``: cosine of the solar zenith angle at the NG jump, 0 in night or cast shadow (shadow map from
  the Sun direction of that observation, checked against the forward model's own shadow flags of the
  facets in that point's kernel file);
* ``slope``: effective gravitational slope (self-gravity + centrifugal; the bundle's gravity cache);
* ``visibility``: facets seen by NG and/or RG at the NG jump (non-zero transport kernels of that point),
  and the late-window NG-only and RG-only templates;
* ``overview``: terrain, cos(SZA), slope and visibility at the NG jump, whole nucleus;
* ``classes``: facets seen by NG only, RG only, or both, (1) at every observation point from the NG jump to
  the last point and (2) at one or more of them (zoom); ``classes_global``: the same classes in
  longitude/latitude over the whole nucleus.

Every descent figure also marks the surface seen by only one gauge at one or more observation points
since the NG jump (NG only: orange, RG only: blue): filled on the plain terrain, outlined on the colour maps.

``<prefix>_data.json`` records the view, the NG gap, the shadow-map check and the class statistics.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))

from cops_cuda import figures as fg  # noqa: E402
from cops_cuda import formats  # noqa: E402

COS75 = np.cos(np.radians(75.0))
BOTH_COLOR = np.array([0x9A, 0xD0, 0x8E]) / 255.0
BOTH_ALL_COLOR = np.array([0x1B, 0x7A, 0x3E]) / 255.0
NIGHT_COLOR = np.array([0x45, 0x55, 0x75]) / 255.0
CLASS_COLORS = {"ng_only": fg.NG_ONLY_COLOR, "rg_only": fg.RG_ONLY_COLOR, "both": BOTH_COLOR}
CLASS_NAMES = {"ng_only": "NG only", "rg_only": "RG only", "both": "NG and RG"}


def read_point(base, point):
    """(header geometry dict, records) of one kernel file."""
    path = Path(base) / "kernels" / ("point=%04d.bin" % point)
    with path.open("rb") as stream:
        header = np.frombuffer(stream.read(formats.KERNEL_HEADER_DTYPE.itemsize), dtype=formats.KERNEL_HEADER_DTYPE)[0]
        records = np.frombuffer(stream.read(), dtype=formats.KERNEL_RECORD_DTYPE)
    if int(header["point"]) != point:
        raise ValueError("kernel point mismatch in %s" % path)
    return dict(zip(formats.KERNEL_GEOMETRY_NAMES, (float(x) for x in header["geometry"]))), records


def class_counts(base, points, nfaces):
    """Per facet, the number of points at which it is seen by NG only, RG only, both, or either."""
    out = {k: np.zeros(nfaces, np.uint16) for k in ("ng_only", "rg_only", "both", "seen")}
    for point in points:
        _, rec = read_point(base, point)
        ng = rec["Kn"] > 0
        rg = (rec["KPhi"] != 0) | (rec["KV"] != 0)
        face = rec["face"]
        out["ng_only"][face[ng & ~rg]] += 1
        out["rg_only"][face[rg & ~ng]] += 1
        out["both"][face[ng & rg]] += 1
        out["seen"][face[ng | rg]] += 1
    return out


def ng_jump(table, rows, at=None, min_gap=60.0):
    """Row indices (last NG before, first NG after) of the NG jump among ``rows``, and the gap length [s].

    The jump is the NG data gap (longer than ``min_gap`` seconds) across which the measured NG pressure rises
    most, or, with ``at`` (UTC), the first NG measurement at or after that time.
    """
    ng = [i for i in rows if table["refNG"][i] > 0]
    if at is not None:
        t = fg.parse_utc(at)
        after = next(i for i in ng if table["time"][i] >= t)
        before = ng[ng.index(after) - 1]
    else:
        gaps = [(table["refNG"][b] / table["refNG"][a], a, b) for a, b in zip(ng[:-1], ng[1:])
                if (table["time"][b] - table["time"][a]).total_seconds() > min_gap]
        _, before, after = max(gaps)
    return before, after, (table["time"][after] - table["time"][before]).total_seconds()


def shadow_map(V, F, N, sun, max_pixels=8192, device_id=0):
    """Lit flag of every facet: facing the Sun and not hidden from it (orthographic depth map along the Sun)."""
    view = fg.OrthoView(sun)
    sx, sy, sz = view.project(V)
    pad = 5.0
    x0, x1, y0, y1 = sx.min() - pad, sx.max() + pad, sy.min() - pad, sy.max() + pad
    pixel = max((x1 - x0) / max_pixels, (y1 - y0) / max_pixels)
    W, H = int(np.ceil((x1 - x0) / pixel)), int(np.ceil((y1 - y0) / pixel))
    face_img, depth = fg.rasterize((sx - x0) / pixel, (y1 - sy) / pixel, sz - sz.min(), F, W, H, device_id=device_id)
    c = V[F].mean(axis=1)
    cx, cy, cz = view.project(c)
    col = np.clip(((cx - x0) / pixel).astype(int), 0, W - 1)
    row = np.clip(((y1 - cy) / pixel).astype(int), 0, H - 1)
    mu = N @ sun
    tol = 0.3 + 2.0 * pixel / np.maximum(mu, 0.1)
    own = face_img[row, col] == np.arange(len(F))
    lit = (mu > 0) & (own | (cz - sz.min() <= depth[row, col] + tol))
    return lit, mu, {"pixel_m": pixel, "width": W, "height": H}


def frame(view, points, margin=0.53):
    """(center, half width) of a square orthographic window around ``points``."""
    px, py, _ = view.project(points)
    center = 0.5 * (px.min() + px.max()) * view.right + 0.5 * (py.min() + py.max()) * view.up
    return center, margin * max(np.ptp(px), np.ptp(py))


def build(args):
    base = Path(args.base)
    table = fg.read_table(base)
    shape = args.shape or fg.shape_of_run(base)
    print("reading %s" % shape, flush=True)
    V, F = fg.load_mesh(shape)
    el = formats.map_elements(base / "elements.bin")
    N = np.asarray(el["normal"], float)
    area = np.asarray(el["area"], float)
    if len(F) != len(N):
        raise ValueError("mesh and elements.bin differ in facet count")
    impact = json.load(open(args.impact)) if args.impact else fg.impact_point(V, F, table)
    ip = np.asarray(impact["point_m"], float)
    dist = np.linalg.norm(np.asarray(el["centroid"], float) - ip, axis=1)
    near = dist <= 100.0
    up = (N[near] * area[near][:, None]).sum(axis=0)
    up /= np.linalg.norm(up)

    t0 = fg.parse_utc(args.start)
    rows = [i for i in range(len(table["time"])) if table["time"][i] >= t0 and table["inside"][i] == 0]
    P = np.column_stack((table["x"], table["y"], table["z"]))[rows]
    G, gd, _ = fg.ground_track(V, F, P)                          # nearest surface point: the altitude
    print("nearest-point distance vs the table's altitude: max |difference| %.3f m" %
          np.max(np.abs(gd - table["altitude"][rows])), flush=True)
    before, after, gap = ng_jump(table, rows, args.jump)
    jump_point = int(table["iPoint"][after])
    print("NG gap %.0f s: %s -> %s (point %d)" % (gap, table["utc"][before], table["utc"][after], jump_point), flush=True)

    # broadside view: the normal of the trajectory's best-fit plane, tilted up about the impact's vertical
    _, sv, vt = np.linalg.svd(P - P.mean(axis=0), full_matrices=False)
    normal = vt[2] - (vt[2] @ up) * up
    normal /= np.linalg.norm(normal)
    el_ = np.radians(args.elevation)
    choice = None
    for sign in (1.0, -1.0):
        view = fg.OrthoView(np.cos(el_) * sign * normal + np.sin(el_) * up, up_hint=up)
        sx = view.project(P)[0]
        if choice is None or sx[0] < sx[-1]:                    # descending from left to right
            choice = view
        if sx[0] < sx[-1]:
            break
    view = choice
    print("trajectory plane from %s: singular values %s; view %s" % (args.start, np.round(sv / sv[0], 3),
                                                                     np.round(view.d, 3)), flush=True)
    whole = frame(view, np.vstack((V, P)))
    zsel = table["altitude"][rows] <= args.zoom_altitude
    geom, rec = read_point(base, jump_point)
    masks = fg.bundle_masks(args.visibility_bundle) if args.visibility_bundle else {}
    late = masks.get("late_ng_not_rg", np.zeros(len(F), bool)) | masks.get("late_rg_not_ng", np.zeros(len(F), bool))
    zpts = [P[zsel], ip[None, :]]
    if late.any():
        zpts.append(V[F[late].ravel()])
    zoom = frame(view, np.vstack(zpts), margin=0.58)

    sun = np.array([geom["sun_x"], geom["sun_y"], geom["sun_z"]])
    sun /= np.linalg.norm(sun)
    lit, mu, smap = shadow_map(V, F, N, sun, device_id=args.device)
    model_lit = (rec["cos"] > 0) & (rec["shadow"] == 0)
    agree = float(np.mean(model_lit == lit[rec["face"]]))
    print("shadow map vs forward model (%d facets in the kernel): %.4f agree" % (len(rec), agree), flush=True)
    slope, cache = fg.gravity_slope(args.bundle, args.gravity_cache, nfaces=len(F))
    ng_seen = np.zeros(len(F), bool)
    rg_seen = np.zeros(len(F), bool)
    ng_seen[rec["face"][rec["Kn"] > 0]] = True
    rg_seen[rec["face"][(rec["KPhi"] != 0) | (rec["KV"] != 0)]] = True
    since = [int(table["iPoint"][i]) for i in rows if i >= after]
    counts = class_counts(base, since, len(F))
    print("visibility classes over %d points since the NG jump" % len(since), flush=True)
    excl_ng, excl_rg = counts["ng_only"] >= 1, counts["rg_only"] >= 1
    n_since = len(since)
    class_faces = excl_ng | excl_rg | (counts["both"] == n_since)
    jsel = np.array([i >= after for i in rows])
    cframe = frame(view, np.vstack((P[jsel], ip[None, :], V[F[class_faces].ravel()])), margin=0.56)
    return {"table": table, "rows": rows, "P": P, "G": G, "V": V, "F": F, "N": N, "area": area, "dist": dist,
            "view": view, "frames": {"whole": whole, "zoom": zoom, "classes": cframe}, "zoom_rows": zsel, "plane_singular": sv / sv[0],
            "impact": impact, "ip": ip, "before": before, "after": after, "gap": gap, "jump_point": jump_point,
            "sun": sun, "cos_sza": np.where(lit, mu, 0.0), "lit_check": {"facets": int(len(rec)), "agreement": agree},
            "shadow_map": smap, "slope": slope, "cache": cache, "ng_seen": ng_seen, "rg_seen": rg_seen,
            "late_ng": masks.get("late_ng_not_rg", np.zeros(len(F), bool)),
            "late_rg": masks.get("late_rg_not_ng", np.zeros(len(F), bool)), "counts": counts, "class_points": since,
            "excl_ng": excl_ng, "excl_rg": excl_rg}


def outline(ax, face, ss, mask, color):
    """Boundary of the facets in ``mask`` as seen in a rendered face image (white halo)."""
    m = np.zeros(face.shape)
    valid = face >= 0
    m[valid] = mask[face[valid]]
    if not m.any():
        return
    x = (np.arange(face.shape[1]) + 0.5) / ss
    y = (np.arange(face.shape[0]) + 0.5) / ss
    ax.contour(x, y, m, levels=[0.5], colors=["w"], linewidths=1.7, zorder=4)
    ax.contour(x, y, m, levels=[0.5], colors=[color], linewidths=0.85, zorder=4)


def render(d, colors, ax, pixels, device_id, shade_weight, where="whole", trajectory="ng", labels=True,
           cmap_ng=None, norm_ng=None, exclusive=None):
    """One panel of the whole nucleus (``where`` = whole or zoom window): trajectory, ground track, marks.

    ``exclusive``: "fill" colours, "outline" outlines the surface seen by NG only / RG only since the NG jump.
    """
    from matplotlib.collections import LineCollection
    center, half = d["frames"][where]
    if exclusive == "fill":
        colors = np.array(colors, copy=True)
        colors[d["excl_rg"]] = fg.RG_ONLY_COLOR
        colors[d["excl_ng"]] = fg.NG_ONLY_COLOR
    image, to_image, visible, face, ss = fg.render_ortho(d["V"], d["F"], d["N"], d["view"], center, half, pixels,
                                                         colors, device_id, shade_weight=shade_weight, return_faces=True)
    ax.imshow(image, interpolation="lanczos", origin="upper", extent=(0, pixels, pixels, 0))
    if exclusive == "outline":
        outline(ax, face, ss, d["excl_rg"], fg.RG_ONLY_COLOR)
        outline(ax, face, ss, d["excl_ng"], fg.NG_ONLY_COLOR)
    table, rows, P, G = d["table"], d["rows"], d["P"], d["G"]
    tol = max(3.0, 8.0 * half / pixels)
    px, py, _ = to_image(P)
    gx, gy, _ = to_image(G)
    pvis = visible(P, tol)
    gvis = visible(G, tol)
    seg_g = np.stack((np.column_stack((gx, gy))[:-1], np.column_stack((gx, gy))[1:]), axis=1)
    keep = gvis[:-1] & gvis[1:]
    ax.add_collection(LineCollection(seg_g[keep], colors="k", linewidths=0.9, zorder=5))
    ax.add_collection(LineCollection(seg_g[keep], colors="w", linewidths=0.35, zorder=5))
    times = table["time"][rows]
    j, b = rows.index(d["after"]), rows.index(d["before"])
    step = 600 if where == "whole" else (120 if where == "zoom" else 60)
    marks = [0]
    for i in range(1, len(rows)):
        if (times[i] - times[0]).total_seconds() >= step * len(marks):
            marks.append(i)
    if where != "whole":
        zr = np.nonzero(d["zoom_rows"])[0]
        marks = [m for m in marks if m >= zr[0]]
        for i in sorted(set(marks[::2] + [b, j])):
            ax.plot([px[i], gx[i]], [py[i], gy[i]], color="0.25", lw=0.5, ls=(0, (2, 1.5)), zorder=6)
    seg = np.stack((np.column_stack((px, py))[:-1], np.column_stack((px, py))[1:]), axis=1)
    shown = pvis[:-1] & pvis[1:]
    ax.add_collection(LineCollection(seg[~shown], colors="0.45", linewidths=0.8, linestyles=(0, (1, 1.5)), zorder=6))
    lc = None
    if trajectory == "ng":
        ng = table["refNG"][rows]
        has = (ng[:-1] > 0) & (ng[1:] > 0)
        ax.add_collection(LineCollection(seg[shown & ~has], colors="0.55", linewidths=2.0, zorder=7))
        lc = LineCollection(seg[shown & has], cmap=cmap_ng, norm=norm_ng, linewidths=2.4, zorder=8)
        lc.set_array(0.5 * (ng[:-1] + ng[1:])[shown & has])
        ax.add_collection(lc)
    else:
        ax.add_collection(LineCollection(seg[shown], colors="w", linewidths=2.6, zorder=7))
        ax.add_collection(LineCollection(seg[shown], colors="k", linewidths=1.3, zorder=8))
    ax.plot(px[marks], py[marks], "o", ms=2.6, mfc="white", mec="k", mew=0.5, zorder=9)
    ix, iy, _ = to_image(d["ip"][None, :])
    ax.plot([ix[0]], [iy[0]], marker="*", ms=10, mfc="#F0E442", mec="k", mew=0.6, zorder=10)
    ax.plot([px[j]], [py[j]], marker="D", ms=4.5, mfc="#D55E00", mec="k", mew=0.6, zorder=10)
    ax.plot([px[b]], [py[b]], marker="D", ms=3.5, mfc="white", mec="k", mew=0.6, zorder=10)
    if labels:
        box = {"facecolor": "white", "edgecolor": "0.3", "lw": 0.4, "boxstyle": "round,pad=0.2"}
        arrow = {"arrowstyle": "-", "lw": 0.5, "color": "0.2"}
        alt = table["altitude"]
        if where == "whole":
            ax.annotate("%s UTC, %.1f km" % (times[0].strftime("%H:%M"), alt[rows[0]] / 1000.0), xy=(px[0], py[0]),
                        xytext=(-10, -14), textcoords="offset points", fontsize=5.6, ha="right", va="top", bbox=box,
                        arrowprops=arrow, zorder=11)
            for m in marks[1:]:
                ax.annotate(times[m].strftime("%H:%M"), xy=(px[m], py[m]), xytext=(-6, -2), textcoords="offset points",
                            fontsize=5.0, ha="right", va="center", zorder=11)
        else:
            ax.annotate("NG jump %s UTC\n%.0f m, NG %.1e Pa" % (times[j].strftime("%H:%M:%S"), alt[d["after"]],
                                                               table["refNG"][d["after"]]),
                        xy=(px[j], py[j]), xytext=(38, 4), textcoords="offset points", fontsize=5.6, va="center",
                        bbox=box, arrowprops=arrow, zorder=11)
            ax.annotate("last NG before the gap\n%s UTC, %.0f m, NG %.1e Pa" % (
                times[b].strftime("%H:%M:%S"), alt[d["before"]], table["refNG"][d["before"]]),
                xy=(px[b], py[b]), xytext=(38, 6), textcoords="offset points", fontsize=5.6, va="center",
                bbox=box, arrowprops=arrow, zorder=11)
            ax.annotate("impact %s UTC" % d["impact"]["utc"].split("T")[1], xy=(ix[0], iy[0]), xytext=(-72, 22),
                        textcoords="offset points", fontsize=5.6, ha="center", bbox=box, arrowprops=arrow, zorder=11)
    ax.set_xlim(0, pixels)
    ax.set_ylim(pixels, 0)
    ax.set_xticks([])
    ax.set_yticks([])
    fg._scale_bar(ax, half, pixels)
    return lc, to_image


def zoom_box(ax, d, to_image, pixels):
    """Outline of the zoom window in a whole-nucleus panel."""
    zc, zh = d["frames"]["zoom"]
    view = d["view"]
    corners = [zc + sx * zh * view.right + sy * zh * view.up for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))]
    x, y, _ = to_image(np.array(corners))
    ax.plot(x, y, color="k", lw=0.6, zorder=9)


def colormap(name):
    import matplotlib
    return matplotlib.colormaps[name] if hasattr(matplotlib, "colormaps") else matplotlib.cm.get_cmap(name)


def colour_fields(d):
    n = len(d["F"])
    out = {"terrain": (np.tile(fg.SURFACE_COLOR, (n, 1)), 1.0)}
    c = colormap("inferno")(0.08 + 0.92 * d["cos_sza"])[:, :3]
    c[d["cos_sza"] <= 0] = NIGHT_COLOR
    out["cos_sza"] = (c, 0.35)
    if d["slope"] is not None:
        out["slope"] = (colormap("viridis")(np.clip(d["slope"], 0, 90) / 90.0)[:, :3], 0.35)
    v = np.tile(fg.SURFACE_COLOR, (n, 1))
    v[d["ng_seen"] & d["rg_seen"]] = BOTH_COLOR
    v[d["ng_seen"] & ~d["rg_seen"]] = fg.NG_ONLY_COLOR
    v[~d["ng_seen"] & d["rg_seen"]] = fg.RG_ONLY_COLOR
    out["visibility_jump"] = (v, 0.6)
    u = np.tile(fg.SURFACE_COLOR, (n, 1))
    u[d["late_rg"]] = fg.RG_ONLY_COLOR
    u[d["late_ng"]] = fg.NG_ONLY_COLOR
    out["visibility_late"] = (u, 0.6)
    return out


def legend_patches(kind, d):
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    t = d["table"]["time"]
    j, b = t[d["after"]].strftime("%H:%M:%S"), t[d["before"]].strftime("%H:%M:%S")
    common = [Line2D([], [], color="k", lw=1.3, label="spacecraft trajectory from %s UTC" % t[d["rows"][0]].strftime("%H:%M")
                     + {"terrain": " (colour: measured NG pressure)",
                        "overview": " (in (a) coloured by measured NG pressure)"}.get(kind, "")),
              Line2D([], [], color="0.45", lw=0.8, ls=(0, (1, 1.5)), label="trajectory hidden behind the nucleus"),
              Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", label="every 10 min (whole), 2 min (zoom)"),
              Line2D([], [], color="k", lw=0.8, label="ground track (nearest surface point)"),
              Line2D([], [], color="0.25", lw=0.6, ls=(0, (2, 1.5)), label="altitude drop lines (zoom)"),
              Line2D([], [], marker="D", ls="none", ms=3.5, mfc="white", mec="k", label="last NG before the gap (%s UTC)" % b),
              Line2D([], [], marker="D", ls="none", ms=4.5, mfc="#D55E00", mec="k", label="NG jump (%s UTC)" % j),
              Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", label="Rosetta impact")]
    excl = []
    if kind in ("terrain", "cos_sza", "slope", "overview", "visibility"):
        style = "fill" if kind == "terrain" else "outline"
        for name, mask, colr in (("NG", d["excl_ng"], fg.NG_ONLY_COLOR), ("RG", d["excl_rg"], fg.RG_ONLY_COLOR)):
            label = "seen only by %s at one or more points since the NG jump (%d facets)" % (name, mask.sum())
            excl.append(Patch(color=colr, label=label) if style == "fill" else
                        Line2D([], [], color=colr, lw=1.2, label=label + (" (filled in (a))" if kind == "overview" else "")))
    common = excl + common
    if kind == "visibility":
        ng, rg = d["ng_seen"], d["rg_seen"]
        return [Patch(color=BOTH_COLOR, label="seen by NG and RG at %s UTC (%d facets)" % (j, (ng & rg).sum())),
                Patch(color=fg.NG_ONLY_COLOR, label="NG only: at %s UTC %d facets; late window %d" % (
                    j, (ng & ~rg).sum(), d["late_ng"].sum())),
                Patch(color=fg.RG_ONLY_COLOR, label="RG only: at %s UTC %d facets; late window %d" % (
                    j, (~ng & rg).sum(), d["late_rg"].sum()))] + common
    if kind == "cos_sza":
        return [Patch(color=NIGHT_COLOR, label="night side or cast shadow at %s UTC" % j)] + common
    return common


def ticks(cb):
    """Readable ticks on the logarithmic NG-pressure bar."""
    lo, hi = cb.mappable.norm.vmin, cb.mappable.norm.vmax
    values = [v for v in (1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5) if lo <= v <= hi]
    cb.set_ticks(values)
    cb.set_ticklabels(["%g" % v for v in values])
    cb.ax.minorticks_off()


def figures(d, prefix, pixels, device_id, overwrite):
    plt = fg.style()
    from matplotlib.colors import LogNorm, Normalize
    import matplotlib
    fields = colour_fields(d)
    table, rows = d["table"], d["rows"]
    ng = table["refNG"][rows]
    norm_ng = LogNorm(vmin=ng[ng > 0].min(), vmax=ng.max())
    cmap_ng = plt.get_cmap("viridis")
    j = table["time"][d["after"]].strftime("%H:%M:%S")
    bars = {"terrain": (cmap_ng, norm_ng, "measured NG pressure along the trajectory [Pa] (grey: no data)"),
            "cos_sza": (matplotlib.colors.ListedColormap(colormap("inferno")(np.linspace(0.08, 1, 256))),
                        Normalize(0, 1), "cos(solar zenith angle) at %s UTC (white line: SZA 75°, the T75 node)" % j),
            "slope": (colormap("viridis"), Normalize(0, 90), "effective gravitational slope [°] (white line: 60°)")}
    W = fg.DOUBLE_COLUMN_IN
    written = []

    def colorbar(fig, kind, cax, lc, short=False):
        cmap, norm, label = bars[kind]
        if short:
            label = {"terrain": "measured NG [Pa]", "cos_sza": "cos(SZA) at %s UTC" % j,
                     "slope": "gravitational slope [°]"}[kind]
        sm = lc if kind == "terrain" else matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(label, fontsize=6.3, labelpad=2)
        cb.ax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=5.8)
        if kind == "terrain":
            ticks(cb)
        if kind == "cos_sza":
            cb.ax.axvline(COS75, color="w", lw=1.0)
        if kind == "slope":
            cb.ax.axvline(60, color="w", lw=1.0)

    def pair(kind, title, field=None, legend_kind=None):
        """Whole nucleus (a) and the zoom (b) side by side."""
        H = 5.7
        fig = plt.figure(figsize=(W, H))
        pw = 0.49
        ph = pw * W / H
        colors, sw = fields[field or kind]
        lc = None
        for k, where in enumerate(("whole", "zoom")):
            ax = fig.add_axes([0.005 + k * 0.50, 0.30, pw, ph])
            out, to_image = render(d, colors, ax, pixels, device_id, sw, where=where,
                                   trajectory="ng" if kind == "terrain" else "plain", cmap_ng=cmap_ng, norm_ng=norm_ng,
                                   exclusive="fill" if kind == "terrain" else "outline")
            lc = lc or out
            if where == "whole":
                zoom_box(ax, d, to_image, pixels)
            ax.text(0.03, 0.97, "(%s)" % "ab"[k], transform=ax.transAxes, va="top", fontsize=8,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
        if kind in bars:
            colorbar(fig, kind, fig.add_axes([0.53, 0.21, 0.44, 0.02]), lc)
        fig.legend(handles=legend_patches(legend_kind or kind, d), loc="lower left", bbox_to_anchor=(0.005, 0.0),
                   fontsize=5.4, frameon=False, ncol=2 if kind not in bars else 1)
        fig.suptitle(title, fontsize=7.5, y=0.995)
        out = fg.save(fig, "%s_%s" % (prefix, kind), overwrite)
        plt.close(fig)
        return out

    alt = table["altitude"]
    span = "%s–%s UTC, %.1f km–%.0f m" % (table["time"][rows[0]].strftime("%H:%M"),
                                         table["time"][rows[-1]].strftime("%H:%M:%S"), alt[rows[0]] / 1000.0, alt[rows[-1]])
    written += pair("terrain", "Final approach and descent onto the impact site (%s): (a) whole nucleus, (b) impact area" % span)
    written += pair("cos_sza", "Illumination at the NG jump (%s UTC): (a) whole nucleus, (b) impact area" % j)
    if "slope" in fields:
        written += pair("slope", "Effective gravitational slope: (a) whole nucleus, (b) impact area")

    # visibility: rows at the NG jump / late-window templates; columns whole / zoom
    H = 9.0
    fig = plt.figure(figsize=(W, H))
    pw = 0.49
    ph = pw * W / H
    for r, field in enumerate(("visibility_jump", "visibility_late")):
        colors, sw = fields[field]
        for k, where in enumerate(("whole", "zoom")):
            ax = fig.add_axes([0.005 + k * 0.50, 0.565 - r * (ph + 0.035), pw, ph])
            _, to_image = render(d, colors, ax, pixels, device_id, sw, where=where, trajectory="plain",
                                 labels=(r == 0), exclusive="outline" if r == 0 else None)
            if where == "whole":
                zoom_box(ax, d, to_image, pixels)
            ax.text(0.03, 0.97, "(%s)" % "abcd"[2 * r + k], transform=ax.transAxes, va="top", fontsize=8,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    fig.legend(handles=legend_patches("visibility", d), loc="lower left", bbox_to_anchor=(0.005, 0.0), fontsize=5.6,
               frameon=False, ncol=2)
    fig.suptitle("Gauge visibility: (a, b) at the NG jump (%s UTC); (c, d) late-window templates %s–%s UTC" % (
        j, d["window"]["start"].split("T")[1], d["window"]["stop"].split("T")[1]), fontsize=7.5, y=0.997)
    written += fg.save(fig, "%s_visibility" % prefix, overwrite)
    plt.close(fig)

    # overview: whole nucleus, four colourings
    kinds = [k for k in ("terrain", "cos_sza", "slope", "visibility_jump") if k in fields]
    H = W * 1.13
    fig = plt.figure(figsize=(W, H))
    pw = 0.49
    ph = pw * W / H
    for k, kind in enumerate(kinds):
        r, c = divmod(k, 2)
        ax = fig.add_axes([0.005 + c * 0.50, 0.53 - r * (ph + 0.015), pw, ph])
        colors, sw = fields[kind]
        lc, _ = render(d, colors, ax, pixels, device_id, sw, where="whole",
                       trajectory="ng" if kind == "terrain" else "plain", labels=(kind == "terrain"),
                       cmap_ng=cmap_ng, norm_ng=norm_ng, exclusive="fill" if kind == "terrain" else "outline")
        ax.text(0.03, 0.97, "(%s)" % "abcd"[k], transform=ax.transAxes, va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
        if kind in bars:
            colorbar(fig, kind, ax.inset_axes([0.53, 0.05, 0.44, 0.025]), lc, short=True)
        else:
            ng_, rg_ = d["ng_seen"], d["rg_seen"]
            ax.text(0.97, 0.04, "green: seen by NG and RG at %s UTC\n(NG only %d, RG only %d facets)" % (
                j, (ng_ & ~rg_).sum(), (~ng_ & rg_).sum()), transform=ax.transAxes, ha="right", fontsize=5.4,
                bbox={"facecolor": "white", "edgecolor": "0.4", "lw": 0.4, "pad": 1.5})
    fig.legend(handles=legend_patches("overview", d), loc="lower left", bbox_to_anchor=(0.005, 0.0), fontsize=5.6,
               frameon=False, ncol=2)
    fig.suptitle("Final approach and descent (%s): terrain, illumination, gravitational slope, gauge visibility" % span,
                 fontsize=7.5, y=0.998)
    written += fg.save(fig, "%s_overview" % prefix, overwrite)
    plt.close(fig)
    return written


def class_masks(d):
    n = len(d["class_points"])
    c = d["counts"]
    return n, {k: c[k] == n for k in CLASS_COLORS}, {k: c[k] >= 1 for k in CLASS_COLORS}


def class_figures(d, prefix, pixels, device_id, overwrite):
    """``<prefix>_classes`` (zoom, 2 x 3, whole nucleus rendered) and ``<prefix>_classes_global`` (lon/lat)."""
    plt = fg.style()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    n, every, once = class_masks(d)
    area = d["area"]
    table = d["table"]
    t_jump = table["time"][d["after"]].strftime("%H:%M:%S")
    t_last = table["time"][d["rows"][-1]].strftime("%H:%M:%S")
    seen_all = d["counts"]["seen"] == n
    d["class_stats"] = {
        "points": n, "interval_utc": [table["utc"][d["after"]], table["utc"][d["rows"][-1]]],
        "every_point": {k: {"facets": int(m.sum()), "area_m2": float(area[m].sum())} for k, m in every.items()},
        "at_least_one_point": {k: {"facets": int(m.sum()), "area_m2": float(area[m].sum())} for k, m in once.items()},
        "seen_every_point_any_class": int(seen_all.sum()),
        "seen_every_point_switching_class": int((seen_all & ~every["both"]).sum())}
    written = []
    W, H = fg.DOUBLE_COLUMN_IN, 6.9
    fig = plt.figure(figsize=(W, H))
    pw = 0.325
    ph = pw * W / H
    rows_def = (("every", every, "(1) seen at every one of the %d points, %s–%s UTC" % (n, t_jump, t_last)),
                ("once", once, "(2) seen at one or more of these points"))
    k = 0
    for r, (key, masks, label) in enumerate(rows_def):
        bottom = 0.57 - r * (ph + 0.08)
        fig.text(0.01, bottom + ph + 0.036, label, fontsize=7.2, fontweight="bold")
        for c, name in enumerate(("ng_only", "rg_only", "both")):
            ax = fig.add_axes([0.01 + c * (pw + 0.005), bottom, pw, ph])
            colors = np.tile(fg.SURFACE_COLOR, (len(d["F"]), 1))
            m = masks[name]
            colors[m] = BOTH_ALL_COLOR if (key == "every" and name == "both") else CLASS_COLORS[name]
            render(d, colors, ax, pixels, device_id, 0.6, where="classes", trajectory="plain", labels=False)
            ax.set_title("%s: %d facets, %.3g km²" % (CLASS_NAMES[name], m.sum(), area[m].sum() / 1e6), fontsize=6.6, pad=2)
            ax.text(0.03, 0.97, "(%s)" % "abcdef"[k], transform=ax.transAxes, va="top", fontsize=7.5,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
            if not m.any():
                ax.text(0.76, 0.80, "no facet", transform=ax.transAxes, ha="center", fontsize=8,
                        bbox={"facecolor": "white", "edgecolor": "0.4", "lw": 0.5, "pad": 2})
            k += 1
    handles = [Patch(color=fg.NG_ONLY_COLOR, label="NG sees the facet, RG does not (Kn > 0, KΦ = KV = 0)"),
               Patch(color=fg.RG_ONLY_COLOR, label="RG sees the facet, NG does not"),
               Patch(color=BOTH_ALL_COLOR, label="both gauges see it at every point"),
               Patch(color=BOTH_COLOR, label="both gauges see it at one or more points"),
               Line2D([], [], color="k", lw=1.3, label="spacecraft trajectory"),
               Line2D([], [], marker="D", ls="none", ms=4, mfc="#D55E00", mec="k", label="NG jump (%s UTC)" % t_jump),
               Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", label="Rosetta impact")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.01, 0.0), fontsize=5.8, frameon=False, ncol=2)
    fig.text(0.01, 0.092, "Impact area seen broadside to the trajectory, with marks every minute (the whole nucleus is "
             "rendered; %s_classes_global maps all of it). %d facets are seen at every point but switch class "
             "(%d NG only and %d RG only at some points)." % (
                 Path(prefix).name, int((seen_all & ~every["both"]).sum()), int((seen_all & once["ng_only"]).sum()),
                 int((seen_all & once["rg_only"]).sum())), ha="left", va="bottom", fontsize=5.6, wrap=True)
    fig.suptitle("Which facets each gauge sees from the NG jump to the end (%s–%s UTC, %.0f–%.0f m)" % (
        t_jump, t_last, table["altitude"][d["after"]], table["altitude"][d["rows"][-1]]), fontsize=7.5, y=0.999)
    written += fg.save(fig, "%s_classes" % prefix, overwrite)
    plt.close(fig)

    # whole nucleus, longitude/latitude: the classes layered so the small ones stay visible
    V, F, N = d["V"], d["F"], d["N"]
    colors = np.tile(fg.SURFACE_COLOR, (len(F), 1))
    colors[once["both"]] = BOTH_COLOR
    colors[every["both"]] = BOTH_ALL_COLOR
    colors[once["rg_only"]] = fg.RG_ONLY_COLOR
    colors[once["ng_only"]] = fg.NG_ONLY_COLOR
    local = once["ng_only"] | once["rg_only"] | every["both"] | (once["both"] & (d["dist"] <= 1000.0))
    lon, lat = fg.lonlat(V[F[local]].mean(axis=1))
    zl, za = (lon.min() - 4.0, lon.max() + 4.0), (lat.min() - 4.0, lat.max() + 4.0)
    ilon, ilat = d["impact"]["lon_deg"], d["impact"]["lat_deg"]
    gl, ga = fg.lonlat(d["G"])
    since = np.array([i >= d["after"] for i in d["rows"]])
    fig = plt.figure(figsize=(W, 6.6))
    ax = fig.add_axes([0.07, 0.60, 0.90, 0.37])
    image, _ = fg.render_lonlat(V, F, N, (-180, 180), (-90, 90), 2 * pixels, colors, device_id, shade_weight=0.6)
    ax.imshow(image, extent=(-180, 180, -90, 90), interpolation="lanczos", aspect="equal")
    ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k", lw=0.6)
    ax.plot([ilon], [ilat], marker="*", ms=8, mfc="#F0E442", mec="k", mew=0.6)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-90, 91, 30))
    ax.set_xlabel("longitude [°]")
    ax.set_ylabel("latitude [°]")
    ax.text(0.01, 0.97, "(a)", transform=ax.transAxes, va="top", fontsize=8, bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    axz = fig.add_axes([0.07, 0.19, 0.90, 0.36])
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    image, _ = fg.render_lonlat(V, F, N, zl, za, 2 * pixels, colors, device_id, shade_weight=0.6)
    axz.imshow(image, extent=(zl[0], zl[1], za[0], za[1]), interpolation="lanczos", aspect=1.0 / coslat)
    axz.plot(gl[since], ga[since], color="k", lw=1.0)
    axz.plot([ilon], [ilat], marker="*", ms=10, mfc="#F0E442", mec="k", mew=0.6)
    axz.set_xlim(*zl)
    axz.set_ylim(*za)
    axz.set_xlabel("longitude [°]")
    axz.set_ylabel("latitude [°]")
    axz.text(0.01, 0.97, "(b)", transform=axz.transAxes, va="top", fontsize=8, bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})

    def lab(name, m):
        return "%s (%d facets, %.3g km²)" % (name, m.sum(), area[m].sum() / 1e6)
    handles = [Patch(color=fg.NG_ONLY_COLOR, label=lab("NG only at one or more points", once["ng_only"])),
               Patch(color=fg.RG_ONLY_COLOR, label=lab("RG only at one or more points", once["rg_only"])),
               Patch(color=BOTH_ALL_COLOR, label=lab("NG and RG at every point", every["both"])),
               Patch(color=BOTH_COLOR, label=lab("NG and RG at one or more points", once["both"])),
               Line2D([], [], color="k", lw=1.0, label="ground track since the NG jump (nearest surface point)"),
               Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", label="Rosetta impact")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.07, 0.0), fontsize=6.0, frameon=False, ncol=2)
    fig.text(0.07, 0.092, "Drawn in this order, each over the previous: NG and RG at one or more points, at every point, "
             "RG only, NG only. No facet is NG only or RG only at every point.\n(b): the area around the impact that holds "
             "every NG-only, RG-only and every-point facet and the NG-and-RG facets within 1000 m.", fontsize=5.8)
    fig.suptitle("Gauge visibility over the whole nucleus, %d points %s–%s UTC" % (n, t_jump, t_last), fontsize=7.5, y=0.997)
    written += fg.save(fig, "%s_classes_global" % prefix, overwrite)
    plt.close(fig)
    return written


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--base", required=True)
    p.add_argument("--bundle", required=True, help="restart template bundle (gravity cache path, cliff templates)")
    p.add_argument("--visibility-bundle", help="bundle with late_ng_not_rg / late_rg_not_ng (default: --bundle)")
    p.add_argument("--gravity-cache")
    p.add_argument("--impact", help="ground-track impact.json (default: recomputed)")
    p.add_argument("--shape")
    p.add_argument("--start", default="2016-09-30T09:00:00", help="trajectory start (UTC); also fixes the view plane")
    p.add_argument("--jump", help="UTC of the NG jump (default: the NG data gap with the largest pressure rise)")
    p.add_argument("--elevation", type=float, default=20.0,
                   help="tilt of the broadside view above the trajectory plane, about the impact's vertical [deg]")
    p.add_argument("--zoom-altitude", type=float, default=1100.0, help="the zoom shows the descent below this altitude [m]")
    p.add_argument("--pixels", type=int, default=900)
    p.add_argument("--device", type=int, default=0)
    p.add_argument("--out", required=True)
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    args.visibility_bundle = args.visibility_bundle or args.bundle
    d = build(args)
    d["window"] = fg.late_window(d["table"], args.visibility_bundle)
    written = figures(d, args.out, args.pixels, args.device, args.overwrite)
    written += class_figures(d, args.out, args.pixels, args.device, args.overwrite)
    info = {"base": str(args.base), "bundle": str(args.bundle), "visibility_bundle": str(args.visibility_bundle),
            "gravity_cache": d["cache"], "start": args.start,
            "view": {"elevation_deg": args.elevation, "towards_viewer": d["view"].d.tolist(),
                     "trajectory_plane_singular_values": d["plane_singular"].tolist(),
                     "zoom_altitude_m": args.zoom_altitude},
            "ng_gap": {"last_before": d["table"]["utc"][d["before"]], "first_after": d["table"]["utc"][d["after"]],
                       "seconds": d["gap"], "jump_point": d["jump_point"],
                       "ng_before_pa": float(d["table"]["refNG"][d["before"]]),
                       "ng_after_pa": float(d["table"]["refNG"][d["after"]])},
            "sun_direction_body": d["sun"].tolist(), "shadow_map": d["shadow_map"],
            "shadow_check_vs_forward_model": d["lit_check"], "visibility_classes": d.get("class_stats"),
            "figures": written}
    Path(args.out + "_data.json").write_text(json.dumps(info, indent=2) + "\n")
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    sys.exit(main())

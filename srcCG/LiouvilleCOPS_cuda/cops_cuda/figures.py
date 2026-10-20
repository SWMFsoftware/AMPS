"""Publication figures: pressures with the late-visibility window, and ground-track maps.

This is the figure format of LiouvilleCOPS_cuda (user request 2026-10-20). Every figure is
written twice, as ``<prefix>.eps`` (vector; images embedded as rasters) and ``<prefix>.png``
(300 dpi). ``tools/pub_figures.py`` is the command line; ``gpu_fit.py``,
``c05_write_bundle.py`` and ``restart_tool.py`` call the same functions, so future runs
produce this format automatically.

Pressure figures (``pressure_figures``):

* rows: NG pressure, RG pressure, NG/RG;
* the time figure has a full-day column and a zoom on the end of the descent;
* the altitude figure has a logarithmic altitude axis, with the descent running left to right;
* the **late-visibility window** is shaded in every panel, and annotated with its UTC and
  altitude range. This is the interval whose observation points define the
  ``late_ng_not_rg`` and ``late_rg_not_ng`` templates (``c05_write_bundle.py``,
  ``--late-start/--late-stop``).

Ground-track figures (``groundtrack_figures``) show the end-of-mission ground track of
Rosetta and its impact point (user decision 2026-10-20: the Rosetta impact only):

* ``<prefix>_3d``: shaded orthographic views of the shape model, rendered on the GPU
  (``kernels/render.cu``, z-buffer, correct for the concave nucleus). The global view looks
  down the final approach, and a zoom shows the late descent;
* ``<prefix>_latlon``: an equirectangular longitude/latitude map, where the outermost surface
  wins at overhangs, with a zoom on the impact area.

Definitions:

* **Ground track** (plotted): the sub-spacecraft intercept point, i.e. the first surface hit on
  the line from the spacecraft toward the body-frame origin (SPICE ``subpnt`` INTERCEPT on the
  shape model). It is continuous along the descent.
* **Nearest surface point** (in the track table): the definition behind the run's ``altitude``
  column. It is not plotted. Far from the nucleus it jumps between the bumps of the
  nucleus's outer envelope, and near the surface it differs from the intercept by up to about
  altitude × tan(slope).
* **Impact point**: the body-fixed trajectory of the last two observation rows, extrapolated
  in a straight line to its first intersection with the shape model. The last row is 0.38 m
  above the surface (SHAP5 1M).

Positions are the ``x y z`` columns of the base run's ``COPS_model.dat``, which are body-fixed
metres from SPICE. Longitude and latitude are planetocentric in that frame.

The figures contain no transparency, because PostScript has none: shading uses solid light
colours.
"""

import datetime as dt
import json
from pathlib import Path

import numpy as np

# c05_write_bundle.py defaults (--late-start/--late-stop); a bundle's metadata takes precedence.
DEFAULT_LATE_WINDOW = ("2016-09-30T10:33:10", "2016-09-30T10:39:26")

# Okabe–Ito colour-blind safe palette.
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
LATE_FILL = "#FBE3C8"          # late-visibility shading (opaque; EPS has no alpha)
LATE_EDGE = "#E69F00"
NG_ONLY_COLOR = np.array([0xD5, 0x5E, 0x00]) / 255.0
RG_ONLY_COLOR = np.array([0x00, 0x72, 0xB2]) / 255.0
SURFACE_COLOR = np.array([0.80, 0.80, 0.80])
BACKGROUND = np.array([1.0, 1.0, 1.0])

SINGLE_COLUMN_IN = 3.5
DOUBLE_COLUMN_IN = 7.2


def style():
    """Apply the publication rcParams and return matplotlib.pyplot (Agg backend)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.0, "lines.markersize": 2.5,
        "legend.frameon": False, "legend.framealpha": 1.0, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02, "ps.fonttype": 42, "pdf.fonttype": 42, "ps.usedistiller": None,
        "axes.grid": False, "figure.dpi": 100})
    return plt


def save(fig, prefix, overwrite=False):
    """Write ``<prefix>.eps`` and ``<prefix>.png`` (300 dpi); return both paths."""
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = [prefix.with_name(prefix.name + ".eps"), prefix.with_name(prefix.name + ".png")]
    if not overwrite:
        existing = [str(p) for p in paths if p.exists()]
        if existing:
            raise FileExistsError("figure files exist (use overwrite): %s" % ", ".join(existing))
    fig.savefig(str(paths[0]), format="eps", dpi=300)
    fig.savefig(str(paths[1]), format="png", dpi=300)
    return [str(p) for p in paths]


# ----------------------------------------------------------------------------- tables
def parse_utc(text):
    """``YYYY-MM-DDTh:mm:ss[.sss]`` (hours may lack the leading zero) -> datetime."""
    day, clock = text.strip().rstrip("Z").split("T")
    h, m, s = clock.split(":")
    return dt.datetime.fromisoformat(day) + dt.timedelta(hours=int(h), minutes=int(m), seconds=float(s))


def read_table(path):
    """Columns of a LiouvilleCOPS or fit ``COPS_model.dat`` (a directory names its table)."""
    path = Path(path)
    if path.is_dir():
        path = path / "COPS_model.dat"
    names, rows = None, []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#iPoint"):
                names = line[1:].split()
                continue
            if line.startswith("#") or not line.strip() or names is None:
                continue
            rows.append(line.split())
    if names is None:
        raise ValueError("%s has no '#iPoint' header" % path)
    col = {n: i for i, n in enumerate(names)}
    out = {"path": str(path), "utc": [r[col["UTC"]] for r in rows]}
    out["time"] = np.array([parse_utc(u) for u in out["utc"]])
    for name in ("iPoint", "et", "x", "y", "z", "altitude", "pNG", "pRG", "refNG", "refRG",
                 "useNG", "useRG", "inside"):
        if name in col:
            out[name] = np.array([float(r[col[name]]) for r in rows])
    return out


def late_window(table, bundle=None, start=None, stop=None):
    """The late-visibility interval and the altitude range of its observation rows.

    The interval comes from ``start/stop`` if given, else from the bundle's
    ``late_ng_not_rg`` metadata, else ``DEFAULT_LATE_WINDOW``. Rows are selected as in
    ``masks.late_points``: start <= UTC <= stop and not inside the nucleus.
    """
    source = "explicit"
    if start is None or stop is None:
        start, stop, source = DEFAULT_LATE_WINDOW[0], DEFAULT_LATE_WINDOW[1], "default (c05_write_bundle)"
        if bundle is not None:
            meta = bundle_metadata(bundle)
            interval = meta.get("templates", {}).get("late_ng_not_rg", {}).get("interval_utc")
            if interval:
                start, stop, source = interval[0], interval[1], "bundle %s" % Path(bundle).name
    t0, t1 = parse_utc(start), parse_utc(stop)
    sel = (table["time"] >= t0) & (table["time"] <= t1)
    if "inside" in table:
        sel &= table["inside"] == 0
    if not sel.any():
        raise ValueError("late window %s .. %s selects no rows" % (start, stop))
    alt = table["altitude"][sel]
    return {"start": start, "stop": stop, "t0": t0, "t1": t1, "source": source,
            "points": [int(p) for p in table["iPoint"][sel]],
            "altitude_min_m": float(alt.min()), "altitude_max_m": float(alt.max())}


def bundle_metadata(path):
    """The metadata of a template bundle (.npz with ``metadata_json``, or its .json)."""
    path = Path(path)
    if path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    with np.load(str(path), allow_pickle=False) as z:
        return json.loads(str(z["metadata_json"])) if "metadata_json" in z.files else {}


def bundle_masks(path):
    """{name: bool facet mask} of a template bundle."""
    with np.load(str(path), allow_pickle=False) as z:
        names = [str(n) for n in z["names"]]
        values = np.asarray(z["values"])
    return {n: values[:, k] > 0 for k, n in enumerate(names)}


def _hm(seconds):
    return "%.0f" % seconds if seconds >= 10 else "%.1f" % seconds


def window_label(window):
    return "late visibility\n%s–%s UTC\n%s–%s m" % (
        window["start"].split("T")[1], window["stop"].split("T")[1],
        _hm(window["altitude_max_m"]), _hm(window["altitude_min_m"]))


# -------------------------------------------------------------------------- pressures
def pressure_figures(base, series, ref, window, prefix, title=None, overwrite=False,
                     zoom_minutes=25.0):
    """Time and altitude figures of measured and modelled NG/RG pressures.

    * ``base``: table of the unfitted base run, or None. Drawn dashed grey.
    * ``series``: list of ``(label, table)`` or ``(label, table, linestyle)``; tables give
      pNG and pRG. Dashed and dotted series are drawn above solid ones, and within a style
      earlier series above later ones, so nearly coincident curves all stay visible.
    * ``ref``: the table whose ``refNG/refRG`` are the data and whose ``useNG/useRG`` mark
      the selection. Hollow markers are points the fit did not use.
    * ``window``: from ``late_window``.

    Writes ``<prefix>_time`` and ``<prefix>_altitude`` (.eps and .png); returns the paths.
    """
    plt = style()
    import matplotlib.dates as mdates
    import matplotlib.transforms as mtransforms
    lines = []
    if base is not None:
        lines.append(("unfitted base", base, (0, (4, 2)), "0.45", 0.9))
    for k, item in enumerate(series):
        label, table = item[0], item[1]
        ls = item[2] if len(item) > 2 else "-"
        lines.append((label, table, ls, PALETTE[k % len(PALETTE)], 1.0))
    # dashed/dotted lines above solid ones; within a style, earlier series on top
    order = {id(item): 4 + len(lines) - k + (0 if item[2] == "-" else 20) for k, item in enumerate(lines)}
    t, alt = ref["time"], ref["altitude"]
    data = {"NG": (ref["refNG"], ref["useNG"]), "RG": (ref["refRG"], ref["useRG"])}
    written = []

    def panels(ax, x):
        for row, gauge in enumerate(("NG", "RG")):
            obs, use = data[gauge]
            ok = obs > 0
            ax[row].plot(x[ok & (use > 0)], obs[ok & (use > 0)], "o", color="k", ms=1.6, mew=0,
                         label="COPS data (fitted)", zorder=3)
            if (ok & (use == 0)).any():
                ax[row].plot(x[ok & (use == 0)], obs[ok & (use == 0)], "o", mfc="none", mec="k", ms=2.2,
                             mew=0.4, label="COPS data (not fitted)", zorder=3)
            for item in lines:
                label, table, ls, color, lw = item
                ax[row].plot(x, table["p" + gauge], ls=ls, color=color, lw=lw, label=label, zorder=order[id(item)])
            ax[row].set_yscale("log")
        ok = (ref["refNG"] > 0) & (ref["refRG"] > 0)
        ax[2].plot(x[ok], ref["refNG"][ok] / ref["refRG"][ok], "o", color="k", ms=1.6, mew=0, zorder=3)
        for item in lines:
            label, table, ls, color, lw = item
            ax[2].plot(x, table["pNG"] / table["pRG"], ls=ls, color=color, lw=lw, zorder=order[id(item)])
        ax[2].set_yscale("log")

    def shade(a, lo, hi):
        a.axvspan(lo, hi, facecolor=LATE_FILL, edgecolor=LATE_EDGE, lw=0.5, zorder=0)

    def window_title(a, xmid, text):
        tr = mtransforms.blended_transform_factory(a.transData, a.transAxes)
        a.text(xmid, 1.02, text, transform=tr, ha="center", va="bottom", fontsize=6.5, color="0.15",
               bbox={"facecolor": LATE_FILL, "edgecolor": LATE_EDGE, "lw": 0.5, "pad": 1.5})

    ylabels = ("NG pressure [Pa]", "RG pressure [Pa]", "NG / RG")
    span_text = "late visibility: %s–%s UTC, %s–%s m" % (
        window["start"].split("T")[1], window["stop"].split("T")[1],
        _hm(window["altitude_max_m"]), _hm(window["altitude_min_m"]))

    # --- time: full descent (left) and the last minutes (right)
    fig, ax = plt.subplots(3, 2, figsize=(DOUBLE_COLUMN_IN, 6.9), sharex="col", sharey="row",
                           gridspec_kw={"width_ratios": [2.2, 1.0], "wspace": 0.06, "hspace": 0.08})
    fig.subplots_adjust(top=0.88)
    tz0 = window["t1"] - dt.timedelta(minutes=zoom_minutes)
    tz1 = window["t1"] + dt.timedelta(minutes=0.6)
    for col in (0, 1):
        a = ax[:, col]
        panels(a, t)
        for r in range(3):
            shade(a[r], window["t0"], window["t1"])
    for r in range(3):
        ax[r, 0].set_ylabel(ylabels[r])
        ax[r, 1].set_xlim(tz0, tz1)
    ax[0, 0].set_xlim(t.min() - dt.timedelta(minutes=10), tz1)
    ax[2, 0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax[2, 1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax[2, 1].xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, 5)))
    ax[2, 0].set_xlabel("UTC, %s" % t[0].strftime("%Y-%m-%d"))
    ax[2, 1].set_xlabel("UTC (last %.0f min)" % zoom_minutes)
    ax[1, 1].text(0.04, 0.05, "shaded: " + window_label(window), transform=ax[1, 1].transAxes, va="bottom",
                  ha="left", fontsize=6.3, color="0.15",
                  bbox={"facecolor": LATE_FILL, "edgecolor": LATE_EDGE, "lw": 0.5, "pad": 1.5})
    for k, a in enumerate(ax[:, 0]):
        a.text(0.01, 0.95, "(%s)" % "ace"[k], transform=a.transAxes, va="top", fontsize=8)
    for k, a in enumerate(ax[:, 1]):
        a.text(0.04, 0.95, "(%s)" % "bdf"[k], transform=a.transAxes, va="top", fontsize=8)
    _legend(fig, ax[0, 0], title, y=0.965)
    written += save(fig, str(prefix) + "_time", overwrite)
    plt.close(fig)

    # --- altitude (log, descent left to right)
    fig, ax = plt.subplots(3, 1, figsize=(DOUBLE_COLUMN_IN * 0.72, 6.9), sharex=True,
                           gridspec_kw={"hspace": 0.08})
    fig.subplots_adjust(top=0.88)
    panels(ax, alt)
    for a in ax:
        shade(a, window["altitude_min_m"], window["altitude_max_m"])
        a.set_xscale("log")
    ax[0].set_xlim(alt.max() * 1.3, max(alt.min(), 1e-3) * 0.7)
    for r in range(3):
        ax[r].set_ylabel(ylabels[r])
        ax[r].text(0.01, 0.95, "(%s)" % "abc"[r], transform=ax[r].transAxes, va="top", fontsize=8)
    ax[2].set_xlabel("altitude above the shape model [m] (descent →)")
    window_title(ax[0], np.sqrt(window["altitude_min_m"] * window["altitude_max_m"]), span_text)
    _legend(fig, ax[0], title, y=0.965)
    written += save(fig, str(prefix) + "_altitude", overwrite)
    plt.close(fig)
    return written


def _legend(fig, ax, title, y=0.96):
    handles, labels = ax.get_legend_handles_labels()
    seen, h2, l2 = set(), [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l)
            h2.append(h)
            l2.append(l)
    fig.legend(h2, l2, loc="upper center", bbox_to_anchor=(0.5, y), ncol=min(3, len(l2)),
               fontsize=7, handlelength=2.6, columnspacing=1.4)
    if title:
        fig.suptitle(title, y=y + 0.025, fontsize=8)


# ------------------------------------------------------------------------ geometry
def load_mesh(bdf):
    """(vertices [m], faces) of a NASTRAN surface, in the face order of ``elements.bin``.

    The reader is the reference builder's (vertices in sorted GRID-id order, CTRIA3 order).
    """
    from . import reference
    builder = reference.load("build_physical_activity_templates")
    vertices, faces = builder.read_nastran_surface(str(bdf))
    return np.ascontiguousarray(vertices, float), np.ascontiguousarray(faces, np.int64)


def shape_of_run(run):
    """The ShapeModel path of a base run (from its ``input.snapshot.in``)."""
    for line in (Path(run) / "input.snapshot.in").read_text(encoding="utf-8").splitlines():
        key, _, value = line.partition("=")
        if key.strip() == "ShapeModel":
            return value.split("#")[0].strip()
    raise ValueError("no ShapeModel in %s/input.snapshot.in" % run)


def _closest_on_triangles(p, a, b, c):
    """Closest points to p (n,3) on triangles (n,k,3) a, b, c; returns points and distances."""
    p = p[:, None, :]
    ab, ac = b - a, c - a
    n = np.cross(ab, ac)
    nn = np.einsum("...i,...i", n, n)
    t = np.einsum("...i,...i", p - a, n) / np.where(nn > 0, nn, 1.0)
    q = p - t[..., None] * n
    # barycentric coordinates of the plane projection
    v0, v1, v2 = ab, ac, q - a
    d00, d01, d11 = (np.einsum("...i,...i", v0, v0), np.einsum("...i,...i", v0, v1),
                     np.einsum("...i,...i", v1, v1))
    d20, d21 = np.einsum("...i,...i", v2, v0), np.einsum("...i,...i", v2, v1)
    den = d00 * d11 - d01 * d01
    den = np.where(den != 0, den, 1.0)
    bv = (d11 * d20 - d01 * d21) / den
    bw = (d00 * d21 - d01 * d20) / den
    inside = (bv >= 0) & (bw >= 0) & (bv + bw <= 1)
    best = np.where(inside[..., None], q, np.nan)
    dist = np.where(inside, np.linalg.norm(p - q, axis=-1), np.inf)
    for s0, s1 in ((a, b), (b, c), (c, a)):
        e = s1 - s0
        ee = np.einsum("...i,...i", e, e)
        u = np.clip(np.einsum("...i,...i", p - s0, e) / np.where(ee > 0, ee, 1.0), 0.0, 1.0)
        r = s0 + u[..., None] * e
        d = np.linalg.norm(p - r, axis=-1)
        take = d < dist
        best = np.where(take[..., None], r, best)
        dist = np.where(take, d, dist)
    return best, dist


def ground_track(vertices, faces, xyz, k=48):
    """Nearest surface points of positions xyz (n,3): (points, distances, face indices)."""
    from scipy.spatial import cKDTree
    centroid = vertices[faces].mean(axis=1)
    _, idx = cKDTree(centroid).query(xyz, k=min(k, len(faces)))
    idx = idx.reshape(len(xyz), -1)
    tri = vertices[faces[idx]]                                          # (n, k, 3, 3)
    pts, dist = _closest_on_triangles(xyz, tri[:, :, 0], tri[:, :, 1], tri[:, :, 2])
    j = np.argmin(dist, axis=1)
    rows = np.arange(len(xyz))
    return pts[rows, j], dist[rows, j], idx[rows, j]


def intercept_track(vertices, faces, xyz, device_id=0):
    """Sub-spacecraft intercept points: the first surface hit on the line from each position
    toward the body-frame origin (SPICE ``subpnt`` INTERCEPT on the shape model).

    Returns (points, face indices); rows without a hit (inside the body) are NaN / -1.
    The nearest hit is found on the GPU (``kernels/render.cu`` ``ray_first_hit``) and its
    intersection recomputed here in FP64.
    """
    import cupy as cp
    from . import cuda_sources
    xyz = np.asarray(xyz, float)
    dirs = -xyz / np.linalg.norm(xyz, axis=1)[:, None]
    nf = len(faces)
    key = np.empty(len(xyz), np.uint64)
    with cp.cuda.Device(int(device_id)):
        kernel = cuda_sources.function("render", "ray_first_hit", device_id)
        dv = cp.asarray(np.ascontiguousarray(vertices, float))
        dfc = cp.asarray(np.ascontiguousarray(faces, np.int64))
        for r0 in range(0, len(xyz), 512):
            r1 = min(r0 + 512, len(xyz))
            best = cp.full(r1 - r0, np.iinfo(np.uint64).max, dtype=cp.uint64)
            kernel((256, r1 - r0), (256,), (dv, dfc, np.int64(nf), cp.asarray(xyz[r0:r1]),
                                            cp.asarray(dirs[r0:r1]), best))
            key[r0:r1] = cp.asnumpy(best)
        del dv, dfc
        cp.get_default_memory_pool().free_all_blocks()
    hit = key != np.iinfo(np.uint64).max
    face = np.where(hit, (key & np.uint64(0xffffffff)).astype(np.int64), -1)
    points = np.full(xyz.shape, np.nan)
    for i in np.nonzero(hit)[0]:
        a, b, c = vertices[faces[face[i]]]
        e1, e2 = b - a, c - a
        h = np.cross(dirs[i], e2)
        s = xyz[i] - a
        q = np.cross(s, e1)
        t = np.dot(e2, q) / np.dot(e1, h)
        points[i] = xyz[i] + t * dirs[i]
    return points, face


def impact_point(vertices, faces, table, radius=60.0):
    """Straight-line extrapolation of the last two rows to the surface (Möller–Trumbore)."""
    from scipy.spatial import cKDTree
    xyz = np.column_stack((table["x"], table["y"], table["z"]))
    p0, p1 = xyz[-2], xyz[-1]
    dt_s = float(table["et"][-1] - table["et"][-2]) if "et" in table else \
        (table["time"][-1] - table["time"][-2]).total_seconds()
    vel = (p1 - p0) / dt_s
    centroid = vertices[faces].mean(axis=1)
    cand = np.array(cKDTree(centroid).query_ball_point(p1, radius), dtype=np.int64)
    a, b, c = (vertices[faces[cand, i]] for i in range(3))
    e1, e2 = b - a, c - a
    h = np.cross(vel, e2)
    det = np.einsum("ij,ij->i", e1, h)
    ok = np.abs(det) > 1e-30
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    s = p1 - a
    u = np.einsum("ij,ij->i", s, h) * inv
    q = np.cross(s, e1)
    v = (q @ vel) * inv
    tt = np.einsum("ij,ij->i", e2, q) * inv
    hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (tt > 0)
    gp, gd, gf = ground_track(vertices, faces, p1[None, :])
    result = {"last_row_utc": table["utc"][-1], "last_row_altitude_m": float(table["altitude"][-1]),
              "speed_m_s": float(np.linalg.norm(vel)), "nearest_point_m": gp[0].tolist(),
              "nearest_distance_m": float(gd[0]), "nearest_face": int(gf[0])}
    if hit.any():
        k = int(np.argmin(np.where(hit, tt, np.inf)))
        point = p1 + tt[k] * vel
        when = table["time"][-1] + dt.timedelta(seconds=float(tt[k]))
        result.update({"method": "straight-line extrapolation of the last two rows",
                       "point_m": point.tolist(), "face": int(cand[k]),
                       "seconds_after_last_row": float(tt[k]),
                       "utc": when.strftime("%Y-%m-%dT%H:%M:%S.") + "%01d" % (when.microsecond // 100000)})
    else:
        when = table["time"][-1]
        result.update({"method": "nearest surface point of the last row (no ray hit)",
                       "point_m": gp[0].tolist(), "face": int(gf[0]), "seconds_after_last_row": 0.0,
                       "utc": table["utc"][-1]})
    r = np.asarray(result["point_m"])
    result["lon_deg"] = float(np.degrees(np.arctan2(r[1], r[0])))
    result["lat_deg"] = float(np.degrees(np.arcsin(r[2] / np.linalg.norm(r))))
    return result


def lonlat(points):
    r = np.linalg.norm(points, axis=-1)
    return (np.degrees(np.arctan2(points[..., 1], points[..., 0])),
            np.degrees(np.arcsin(np.clip(points[..., 2] / r, -1, 1))))


# --------------------------------------------------------------------------- rendering
def rasterize(u, v, w, faces, width, height, draw=None, device_id=0, max_box=None):
    """Face index per pixel (-1 = background) and its depth, via kernels/render.cu."""
    import cupy as cp
    from . import cuda_sources
    nf = len(faces)
    draw = np.ones(nf, np.uint8) if draw is None else np.asarray(draw, np.uint8)
    max_box = int(max_box if max_box is not None else max(width, height))
    with cp.cuda.Device(int(device_id)):
        zbuf = cp.full(width * height, np.iinfo(np.uint64).max, dtype=cp.uint64)
        args = [cp.asarray(np.ascontiguousarray(u, float)), cp.asarray(np.ascontiguousarray(v, float)),
                cp.asarray(np.ascontiguousarray(w, float)), cp.asarray(np.ascontiguousarray(faces, np.int64)),
                cp.asarray(draw), np.int64(nf), np.int32(width), np.int32(height), np.int32(max_box), zbuf]
        kernel = cuda_sources.function("render", "raster_faces", device_id)
        threads = 256
        blocks = int(min((nf + threads - 1) // threads, 65535 * 4))
        kernel((max(blocks, 1),), (threads,), tuple(args))
        key = cp.asnumpy(zbuf).reshape(height, width)
        del args, zbuf
        cp.get_default_memory_pool().free_all_blocks()
    empty = key == np.iinfo(np.uint64).max
    face = np.where(empty, -1, (key & np.uint64(0xffffffff)).astype(np.int64))
    depth = np.where(empty, np.inf, (key >> np.uint64(32)).astype(np.uint32).view(np.float32).astype(float))
    return face, depth


def _shade(face_img, colors, intensity, ss):
    """RGB image from per-face colours and shading, averaged over ss x ss supersamples."""
    rgb = np.empty(face_img.shape + (3,))
    bg = face_img < 0
    f = np.where(bg, 0, face_img)
    rgb[...] = colors[f] * intensity[f][..., None]
    rgb[bg] = BACKGROUND
    if ss > 1:
        h, w = face_img.shape
        rgb = rgb.reshape(h // ss, ss, w // ss, ss, 3).mean(axis=(1, 3))
    return np.clip(rgb, 0, 1)


def face_colors(nfaces, masks):
    colors = np.tile(SURFACE_COLOR, (nfaces, 1))
    if masks is not None:
        if "late_rg_not_ng" in masks:
            colors[masks["late_rg_not_ng"]] = RG_ONLY_COLOR
        if "late_ng_not_rg" in masks:
            colors[masks["late_ng_not_rg"]] = NG_ONLY_COLOR
    return colors


class OrthoView:
    """Orthographic camera: looks along -``towards`` (from the viewer to the scene)."""

    def __init__(self, towards_viewer, up_hint=(0.0, 0.0, 1.0)):
        d = np.asarray(towards_viewer, float)
        self.d = d / np.linalg.norm(d)
        up = np.asarray(up_hint, float) - np.dot(up_hint, self.d) * self.d
        if np.linalg.norm(up) < 1e-6:
            up = np.array([1.0, 0, 0]) - self.d[0] * self.d
        self.up = up / np.linalg.norm(up)
        self.right = np.cross(self.up, self.d)

    def project(self, p):
        return p @ self.right, p @ self.up, -(p @ self.d)       # screen x, screen y, depth


def render_ortho(vertices, faces, normals, view, center, half_width, pixels, colors, device_id=0, ss=2,
                 light=None):
    """Shaded orthographic image (H x W x 3) of the mesh around ``center``; returns image + mapper."""
    W = H = int(pixels) * ss
    sx, sy, sz = view.project(vertices - center)
    scale = W / (2.0 * half_width)
    u = (sx + half_width) * scale
    v = (half_width - sy) * scale                                  # image rows go down
    w = sz - sz.min()
    face, depth = rasterize(u, v, w, faces, W, H, device_id=device_id)
    L = light if light is not None else (0.75 * view.d - 0.35 * view.right + 0.45 * view.up)
    L = L / np.linalg.norm(L)
    lam = np.clip(normals @ L, 0.0, 1.0)
    facing = normals @ view.d
    intensity = 0.30 + 0.70 * lam
    intensity = np.where(facing < 0, 0.25, intensity)              # back faces seen through gaps
    image = _shade(face, colors, intensity, ss)

    def to_image(points):
        px, py, pz = view.project(np.atleast_2d(points) - center)
        return (px + half_width) * scale / ss, (half_width - py) * scale / ss, pz - sz.min()

    def visible(points, tol_m):
        ix, iy, pz = to_image(points)
        col = np.clip((ix * ss).astype(int), 0, W - 1)
        row = np.clip((iy * ss).astype(int), 0, H - 1)
        inside = (ix >= 0) & (ix < W / ss) & (iy >= 0) & (iy < H / ss)
        return inside & (pz <= depth[row, col] + tol_m)

    return image, to_image, visible


def render_lonlat(vertices, faces, normals, lon_range, lat_range, width, colors, device_id=0, ss=2,
                  shade_weight=1.0):
    """Equirectangular image of the outermost surface over a longitude/latitude box.

    ``shade_weight`` blends the hill shading into the facet colours: 1 = full relief shading,
    0 = the colours exactly (for quantitative colour maps).
    """
    (lo0, lo1), (la0, la1) = lon_range, lat_range
    W = int(width) * ss
    H = int(round(width * (la1 - la0) / (lo1 - lo0))) * ss
    tri = vertices[faces]                                            # (nf, 3, 3)
    lon, lat = lonlat(tri)
    r = np.linalg.norm(tri, axis=-1)
    span = lon.max(axis=1) - lon.min(axis=1)
    seam = span > 180.0
    # seam faces: unwrap their negative longitudes by +360 and draw them twice (at +0 and -360)
    lon_w = np.where(seam[:, None] & (lon < 0), lon + 360.0, lon)
    lon_all = np.concatenate((lon_w, lon_w[seam] - 360.0))
    lat_all = np.concatenate((lat, lat[seam]))
    r_all = np.concatenate((r, r[seam]))
    index = np.concatenate((np.arange(len(faces)), np.nonzero(seam)[0]))
    u = (lon_all - lo0) / (lo1 - lo0) * W
    v = (la1 - lat_all) / (la1 - la0) * H
    w = r_all.max() - r_all
    corner = np.arange(3 * len(index), dtype=np.int64).reshape(-1, 3)
    face_k, _ = rasterize(u.ravel(), v.ravel(), w.ravel(), corner, W, H, device_id=device_id,
                          max_box=max(W, H) // 2)
    face = np.where(face_k < 0, -1, index[np.maximum(face_k, 0)])
    # local light from the north-west, 40 degrees above the horizon
    c = tri.mean(axis=1)
    up = c / np.linalg.norm(c, axis=1)[:, None]
    east = np.cross([0.0, 0.0, 1.0], up)
    east /= np.maximum(np.linalg.norm(east, axis=1), 1e-12)[:, None]
    north = np.cross(up, east)
    L = np.cos(np.radians(40)) * (north - east) / np.sqrt(2) + np.sin(np.radians(40)) * up
    lam = np.clip(np.einsum("ij,ij->i", normals, L), 0.0, 1.0)
    intensity = 1.0 - shade_weight * (0.70 - 0.70 * lam)
    image = _shade(face, colors, intensity, ss)

    def to_image(lo, la):
        return (np.asarray(lo) - lo0) / (lo1 - lo0) * (W / ss), (la1 - np.asarray(la)) / (la1 - la0) * (H / ss)

    return image, to_image


# ------------------------------------------------------------------ ground-track figures
def groundtrack_data(base_run, shape=None, bundle=None, window=None, device_id=0):
    """Mesh, tracks (sub-spacecraft intercept and nearest point), impact point and template
    masks for the ground-track figures."""
    from . import formats
    base_run = Path(base_run)
    table = read_table(base_run)
    shape = shape or shape_of_run(base_run)
    vertices, faces = load_mesh(shape)
    elements, _ = formats.load_elements(base_run / "elements.bin")
    centroid = vertices[faces].mean(axis=1)
    if len(centroid) != len(elements["centroid"]) or np.abs(centroid - elements["centroid"]).max() > 1e-6:
        raise ValueError("shape %s does not match %s/elements.bin (face count or order)" % (shape, base_run))
    xyz = np.column_stack((table["x"], table["y"], table["z"]))
    outside = table["inside"] == 0 if "inside" in table else np.ones(len(xyz), bool)
    near, near_d, near_f = ground_track(vertices, faces, xyz)
    sub, sub_f = intercept_track(vertices, faces, xyz, device_id)
    outside &= sub_f >= 0
    window = window or late_window(table, bundle)
    return {"table": table, "shape": str(shape), "vertices": vertices, "faces": faces,
            "normals": elements["normal"], "ground": np.where(np.isnan(sub), near, sub), "ground_face": sub_f,
            "nearest": near, "nearest_distance": near_d, "nearest_face": near_f,
            "outside": outside, "impact": impact_point(vertices, faces, table),
            "masks": bundle_masks(bundle) if bundle else None, "bundle": str(bundle) if bundle else None,
            "window": window}


def _colored_track(ax, x, y, values, norm, cmap, lw, mask=None, seg_mask=None, zorder=5,
                   solid_color=None, dotted=False):
    """Track coloured by ``values``; a segment is drawn when both ends pass ``mask`` and it
    passes ``seg_mask`` (e.g. no longitude seam crossing)."""
    from matplotlib.collections import LineCollection
    pts = np.column_stack((x, y))
    seg = np.stack((pts[:-1], pts[1:]), axis=1)
    keep = np.ones(len(seg), bool) if mask is None else (mask[:-1] & mask[1:])
    if seg_mask is not None:
        keep &= seg_mask
    if solid_color is not None:
        lc = LineCollection(seg[keep], colors=solid_color, linewidths=lw, zorder=zorder,
                            linestyles=(0, (1, 1.5)) if dotted else "solid")
    else:
        lc = LineCollection(seg[keep], cmap=cmap, norm=norm, linewidths=lw, zorder=zorder, capstyle="round")
        lc.set_array(0.5 * (values[:-1] + values[1:])[keep])
    ax.add_collection(lc)
    return lc


def groundtrack_figures(gt, prefix, device_id=0, overwrite=False, title=None, pixels=800):
    """Write ``<prefix>_3d`` and ``<prefix>_latlon`` (.eps, .png) and ``<prefix>_track.dat``.

    ``<prefix>_3d``: (a) the whole nucleus seen down the final approach, (b) the template
    regions around the impact, (c) a close-up of the late-visibility descent with minute marks.
    ``<prefix>_latlon``: (a) the global longitude/latitude map, (b) the impact area.
    The track is the nearest-surface-point track coloured by altitude. It is thick inside the
    late-visibility window, and dotted grey where hidden behind the nucleus (3D views). Where the
    nearest point jumps between distant surface parts (lobe switches at high altitude), the line
    is interrupted instead of drawn across.
    """
    plt = style()
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    table, window, impact = gt["table"], gt["window"], gt["impact"]
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    alt = np.maximum(table["altitude"], 0.1)
    norm = LogNorm(vmin=max(alt.min(), 0.1), vmax=alt.max())
    cmap = plt.get_cmap("viridis")
    late = np.zeros(len(alt), bool)
    late[[int(np.nonzero(table["iPoint"] == p)[0][0]) for p in window["points"]]] = True
    colors = face_colors(len(F), gt["masks"])
    ip = np.asarray(impact["point_m"])
    ground = gt["ground"]
    step = np.linalg.norm(np.diff(ground, axis=0), axis=1)
    no_jump = step <= max(150.0, 5.0 * np.median(step))
    written = []
    legend_extra = [Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                           label="Rosetta impact (lon %.2f°, lat %.2f°)" % (impact["lon_deg"], impact["lat_deg"])),
                    Line2D([], [], marker="s", ls="none", ms=4, mfc="white", mec="k", mew=0.6,
                           label="track start, %s UTC, %.1f km" % (table["time"][0].strftime("%H:%M:%S"),
                                                                    table["altitude"][0] / 1000.0))]

    # ---------------- 3D orthographic views
    approach = np.column_stack((table["x"], table["y"], table["z"]))[late][0] - ip
    view = OrthoView(approach)
    body_half = 1.04 * np.max(np.abs(np.column_stack(view.project(V - V.mean(axis=0))[:2])))
    zoom_half = 1.15 * max(np.max(np.linalg.norm(ground[late] - ip, axis=1)), 1.0)
    if gt["masks"] is not None:
        for name in ("late_ng_not_rg", "late_rg_not_ng"):
            if name in gt["masks"] and gt["masks"][name].any():
                cen = V[F[gt["masks"][name]]].mean(axis=1)
                zoom_half = max(zoom_half, 1.08 * np.max(np.abs(np.column_stack(view.project(cen - ip)[:2]))))
    close_half = 1.3 * max(np.max(np.abs(np.column_stack(view.project(ground[late] - ip)[:2]))), 20.0)
    panels = ((V.mean(axis=0), body_half), (ip, zoom_half), (ip, close_half))
    fig = plt.figure(figsize=(DOUBLE_COLUMN_IN, 3.35))
    width = 0.315
    axs = [fig.add_axes([0.005 + k * (width + 0.0175), 0.27, width, width * DOUBLE_COLUMN_IN / 3.35])
           for k in range(3)]
    lc = None
    for k, (ax, (center, half)) in enumerate(zip(axs, panels)):
        image, to_image, visible = render_ortho(V, F, N, view, center, half, pixels, colors, device_id)
        ax.imshow(image, interpolation="lanczos", origin="upper", extent=(0, pixels, pixels, 0))
        gx, gy, _ = to_image(ground)
        tol = max(3.0, 4.0 * 2.0 * half / pixels)
        vis = visible(ground, tol) & gt["outside"]
        hidden = ~vis & gt["outside"]
        _colored_track(ax, gx, gy, alt, norm, cmap, 0.6, mask=hidden, seg_mask=no_jump, zorder=4,
                       solid_color="0.45", dotted=True)
        _colored_track(ax, gx, gy, alt, norm, cmap, 2.6, mask=vis & late, seg_mask=no_jump, zorder=5)
        lc = _colored_track(ax, gx, gy, alt, norm, cmap, 1.1, mask=vis, seg_mask=no_jump, zorder=6)
        if vis[0]:
            ax.plot([gx[0]], [gy[0]], marker="s", ms=4, mfc="white", mec="k", mew=0.6, zorder=9)
        ix, iy, _ = to_image(ip[None, :])
        ax.plot([ix[0]], [iy[0]], marker="*", ms=9 if k < 2 else 11, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
        if k < 2:                                   # outline of the next panel's area
            s = panels[k + 1][1] / half * pixels / 2
            ax.plot([ix[0] - s, ix[0] + s, ix[0] + s, ix[0] - s, ix[0] - s],
                    [iy[0] - s, iy[0] - s, iy[0] + s, iy[0] + s, iy[0] - s], color="k", lw=0.6, zorder=7)
        if k == 2:
            _minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=True)
            ax.annotate("Rosetta impact\n%s UTC" % impact["utc"].split("T")[1], xy=(ix[0], iy[0]),
                        xytext=(-45, -55), textcoords="offset points", fontsize=6.3, ha="center",
                        arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": "k", "shrinkB": 6},
                        bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.5, "boxstyle": "round,pad=0.25"},
                        zorder=10)
        else:
            _minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=False)
        ax.set_xlim(0, pixels)
        ax.set_ylim(pixels, 0)
        ax.set_xticks([])
        ax.set_yticks([])
        _scale_bar(ax, half, pixels)
        ax.text(0.03, 0.97, "(%s)" % "abc"[k], transform=ax.transAxes, va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    _region_legend(fig, gt["masks"], legend_extra, anchor=(0.005, 0.0), loc="lower left", ncol=2)
    cax = fig.add_axes([0.70, 0.10, 0.29, 0.035])
    cb = fig.colorbar(lc, cax=cax, orientation="horizontal")
    cb.set_label("spacecraft altitude [m]", fontsize=7)
    if title:
        fig.suptitle(title, fontsize=8, y=1.0)
    written += save(fig, str(prefix) + "_3d", overwrite)
    plt.close(fig)

    # ---------------- longitude/latitude map
    glon, glat = lonlat(ground)
    no_seam = (np.abs(np.diff(glon)) <= 180.0) & no_jump
    lons, lats = [glon[late]], [glat[late]]
    if gt["masks"] is not None:
        for name in ("late_ng_not_rg", "late_rg_not_ng"):
            if name in gt["masks"] and gt["masks"][name].any():
                lo, la = lonlat(V[F[gt["masks"][name]]].mean(axis=1))
                lons.append(lo)
                lats.append(la)
    lons, lats = np.concatenate(lons), np.concatenate(lats)
    pad = 0.10 * max(np.ptp(lons), np.ptp(lats), 1.0)
    zl = (lons.min() - pad, lons.max() + pad)
    za = (lats.min() - pad, lats.max() + pad)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    fig = plt.figure(figsize=(DOUBLE_COLUMN_IN, 7.0))
    ax = fig.add_axes([0.08, 0.60, 0.90, 0.38])
    image, _ = render_lonlat(V, F, N, (-180, 180), (-90, 90), 2 * pixels, colors, device_id)
    ax.imshow(image, extent=(-180, 180, -90, 90), interpolation="lanczos", aspect="equal")
    _colored_track(ax, glon, glat, alt, norm, cmap, 2.6, mask=late & gt["outside"], seg_mask=no_seam, zorder=5)
    lc = _colored_track(ax, glon, glat, alt, norm, cmap, 1.0, mask=gt["outside"], seg_mask=no_seam, zorder=6)
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
    axz = fig.add_axes([0.08, 0.13, 0.90, 0.40])
    image, _ = render_lonlat(V, F, N, zl, za, 2 * pixels, colors, device_id)
    axz.imshow(image, extent=(zl[0], zl[1], za[0], za[1]), interpolation="lanczos", aspect=1.0 / coslat)
    _colored_track(axz, glon, glat, alt, norm, cmap, 2.6, mask=late, seg_mask=no_seam, zorder=5)
    _colored_track(axz, glon, glat, alt, norm, cmap, 1.1, mask=gt["outside"], seg_mask=no_seam, zorder=6)
    _minute_marks(axz, gt, late, lambda p: lonlat(p), label=False)
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
    _region_legend(fig, gt["masks"], legend_extra, anchor=(0.08, 0.0), loc="lower left", ncol=2, hidden=False)
    cax = fig.add_axes([0.70, 0.035, 0.28, 0.014])
    cb = fig.colorbar(lc, cax=cax, orientation="horizontal")
    cb.set_label("spacecraft altitude [m]", fontsize=7)
    if title:
        fig.suptitle(title, fontsize=8, y=1.0)
    written += save(fig, str(prefix) + "_latlon", overwrite)
    plt.close(fig)
    written.append(write_track_table(gt, str(prefix) + "_track.dat", overwrite))
    return written


REGION_COLOR = np.array([0xCC, 0x79, 0xA7]) / 255.0


def region_map(gt, mask, prefix, label, device_id=0, overwrite=False, title=None, pixels=800):
    """Longitude/latitude map of a facet region (e.g. a restart_tool enhancement region).

    (a) global map, (b) zoom on the region; ground track and impact as in
    ``groundtrack_figures``. Writes ``<prefix>`` (.eps, .png); returns the paths.
    """
    plt = style()
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    table, impact = gt["table"], gt["impact"]
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    mask = np.asarray(mask, bool)
    colors = face_colors(len(F), None)
    colors[mask] = REGION_COLOR
    alt = np.maximum(table["altitude"], 0.1)
    norm = LogNorm(vmin=max(alt.min(), 0.1), vmax=alt.max())
    cmap = plt.get_cmap("viridis")
    glon, glat = lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    fig = plt.figure(figsize=(DOUBLE_COLUMN_IN, 7.0))
    ax = fig.add_axes([0.08, 0.60, 0.90, 0.38])
    image, _ = render_lonlat(V, F, N, (-180, 180), (-90, 90), 2 * pixels, colors, device_id)
    ax.imshow(image, extent=(-180, 180, -90, 90), interpolation="lanczos", aspect="equal")
    lc = _colored_track(ax, glon, glat, alt, norm, cmap, 1.0, mask=gt["outside"], seg_mask=ok_seg, zorder=6)
    ax.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=9, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-90, 91, 30))
    ax.set_xlabel("longitude [°]")
    ax.set_ylabel("latitude [°]")
    ax.text(0.01, 0.97, "(a)", transform=ax.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    if mask.any():
        lo, la = lonlat(V[F[mask]].mean(axis=1))
        if np.ptp(lo) > 180:                       # region across the seam: zoom on the whole map
            zl, za = (-180.0, 180.0), (-90.0, 90.0)
        else:
            pad = 0.10 * max(np.ptp(lo), np.ptp(la), 2.0)
            zl, za = (lo.min() - pad, lo.max() + pad), (max(la.min() - pad, -90), min(la.max() + pad, 90))
        ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k", lw=0.6, zorder=7)
        axz = fig.add_axes([0.08, 0.13, 0.90, 0.40])
        image, _ = render_lonlat(V, F, N, zl, za, 2 * pixels, colors, device_id)
        coslat = max(np.cos(np.radians(0.5 * (za[0] + za[1]))), 0.2)
        axz.imshow(image, extent=(zl[0], zl[1], za[0], za[1]), interpolation="lanczos", aspect=1.0 / coslat)
        _colored_track(axz, glon, glat, alt, norm, cmap, 1.1, mask=gt["outside"], seg_mask=ok_seg, zorder=6)
        axz.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=11, mfc="#F0E442", mec="k", mew=0.6,
                 zorder=9)
        axz.set_xlim(*zl)
        axz.set_ylim(*za)
        axz.set_xlabel("longitude [°]")
        axz.set_ylabel("latitude [°]")
        axz.text(0.01, 0.97, "(b)", transform=axz.transAxes, va="top", fontsize=8,
                 bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    handles = [Patch(facecolor=REGION_COLOR, edgecolor="none", label="%s (%d facets)" % (label, int(mask.sum()))),
               Line2D([], [], color="#21918C", lw=1.1, label="ground track (sub-spacecraft intercept)"),
               Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                      label="Rosetta impact (lon %.2f°, lat %.2f°)" % (impact["lon_deg"], impact["lat_deg"]))]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.08, 0.0), fontsize=6.3, ncol=1, frameon=False)
    cax = fig.add_axes([0.70, 0.035, 0.28, 0.014])
    cb = fig.colorbar(lc, cax=cax, orientation="horizontal")
    cb.set_label("spacecraft altitude [m]", fontsize=7)
    if title:
        fig.suptitle(title, fontsize=8, y=1.0)
    paths = save(fig, prefix, overwrite)
    plt.close(fig)
    return paths


def _region_legend(fig, masks, extra=(), anchor=(0.0, 0.0), loc="lower left", ncol=2, hidden=True):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = []
    if masks is not None and "late_ng_not_rg" in masks:
        handles.append(Patch(facecolor=NG_ONLY_COLOR, edgecolor="none",
                             label="NG-only facets (%d)" % int(masks["late_ng_not_rg"].sum())))
    if masks is not None and "late_rg_not_ng" in masks:
        handles.append(Patch(facecolor=RG_ONLY_COLOR, edgecolor="none",
                             label="RG-only facets (%d)" % int(masks["late_rg_not_ng"].sum())))
    handles.append(Line2D([], [], color="#21918C", lw=2.6, label="ground track, late-visibility window"))
    handles.append(Line2D([], [], color="#21918C", lw=1.1, label="ground track (sub-spacecraft intercept)"))
    if hidden:
        handles.append(Line2D([], [], color="0.45", lw=0.6, ls=(0, (1, 1.5)), label="ground track, hidden"))
    handles.append(Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", mew=0.5,
                          label="1-min marks from %s" % "window start"))
    handles.extend(extra)
    fig.legend(handles=handles, loc=loc, bbox_to_anchor=anchor, fontsize=6.3, ncol=ncol, frameon=False,
               handlelength=2.2, columnspacing=1.0)


def _scale_bar(ax, half_width, pixels):
    """A scale bar of a round length about one fifth of the panel width."""
    target = 0.2 * 2 * half_width
    length = 10 ** np.floor(np.log10(target))
    for m in (5, 2, 1):
        if m * length <= target:
            length *= m
            break
    px = length / (2 * half_width) * pixels
    x0, y0 = 0.06 * pixels, 0.93 * pixels
    ax.plot([x0, x0 + px], [y0, y0], color="k", lw=1.6, solid_capstyle="butt", zorder=8)
    ax.text(x0 + px / 2, y0 - 0.015 * pixels, "%g %s" % ((length / 1000, "km") if length >= 1000 else (length, "m")),
            ha="center", va="bottom", fontsize=6.5, zorder=8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6})


def _minute_marks(ax, gt, late, project, label):
    """Marks along the late track at the window start and every full minute after it.

    With ``label``, a mark is labelled (UTC, altitude) when it lies at least 5.5% of the panel
    width from the previously labelled one; labels sit perpendicular to the track, alternating
    sides, with a leader line.
    """
    table = gt["table"]
    idx = np.nonzero(late)[0]
    t0 = table["time"][idx[0]]
    marks = [idx[0]]
    for i in idx[1:]:
        if (table["time"][i] - t0).total_seconds() >= 60 * len(marks):
            marks.append(i)
    x, y = project(gt["ground"][marks])
    x, y = np.atleast_1d(x).astype(float), np.atleast_1d(y).astype(float)
    ax.plot(x, y, "o", ms=2.8, mfc="white", mec="k", mew=0.5, ls="none", zorder=8)
    if not label:
        return
    xl = ax.get_xlim()
    min_gap = 0.055 * abs(xl[1] - xl[0])
    last, side = None, 1
    for k, (xi, yi, i) in enumerate(zip(x, y, marks)):
        if last is not None and np.hypot(xi - last[0], yi - last[1]) < min_gap:
            continue
        j0, j1 = max(k - 1, 0), min(k + 1, len(x) - 1)
        tx, ty = x[j1] - x[j0], y[j1] - y[j0]
        n = np.hypot(tx, ty)
        px, py = (-ty / n, tx / n) if n > 0 else (0.0, 1.0)
        ax.annotate("%s UTC\n%.0f m" % (table["time"][i].strftime("%H:%M:%S"), table["altitude"][i]),
                    xy=(xi, yi), xytext=(26 * side * px, -26 * side * py), textcoords="offset points",
                    fontsize=5.4, ha="center", va="center", zorder=8,
                    arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "0.2"},
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4})
        last, side = (xi, yi), -side


# ------------------------------------------------------------------- template maps
CLIFF45_COLOR = np.array([0x56, 0xB4, 0xE9]) / 255.0
CLIFF60_COLOR = np.array([0x00, 0x5A, 0x8C]) / 255.0


def gravity_slope(bundle, cache=None, nfaces=None):
    """Per-facet effective gravitational slope [deg] from the bundle's admitted gravity cache."""
    meta = bundle_metadata(bundle)
    path = cache or (((meta.get("gravity") or {}).get("cache") or {}).get("file") or {}).get("path")
    if not path or not Path(path).exists():
        return None, None
    with np.load(str(path), allow_pickle=False) as z:
        slope = np.asarray(z["reference_slope_deg"], float)
    if nfaces is not None and len(slope) != nfaces:
        raise ValueError("gravity cache %s has %d facets, the mesh %d" % (path, len(slope), nfaces))
    return slope, str(path)


def late_visibility_definition(window, meta=None):
    """Plain-language definition of the late-visibility templates (figure text and captions)."""
    n = len(window["points"])
    if meta:
        n = meta.get("templates", {}).get("late_ng_not_rg", {}).get("selected_point_count", n)
    return ("Late visibility: the %d COPS observation points of %s–%s UTC (spacecraft altitude %s–%s m). "
            "NG-only facets are seen by the nude gauge (NG) and not by the ram gauge (RG) from at least one "
            "of these points: Kn > 0 and KΦ = KV = 0 at the same point. RG-only facets are the reverse: "
            "KΦ ≠ 0 or KV ≠ 0, and Kn = 0, at the same point. Kn = Σ (dS/r²) A(θ) W_NG and "
            "KΦ, KV = Σ (dS/r²) A(θ) Acc_RG (−l·b_RG, resp. 1) are the collisionless (Liouville) transport "
            "kernels of the facet at the gauge, summed over its parts visible from the spacecraft "
            "(LiouvilleCOPS src/Liouville.h), so 'seen' means a nonzero instrument response, not only a line "
            "of sight. The union is taken over the points." % (
                n, window["start"].split("T")[1], window["stop"].split("T")[1],
                _hm(window["altitude_max_m"]), _hm(window["altitude_min_m"])))


def write_caption(prefix, text, overwrite=False):
    path = Path(str(prefix) + "_caption.txt")
    if path.exists() and not overwrite:
        raise FileExistsError(str(path))
    path.write_text(text.strip() + "\n", encoding="utf-8")
    return str(path)


def templates_figure(gt, prefix, gravity_cache=None, device_id=0, overwrite=False, title=None, pixels=900):
    """Template maps with the ground track: the publication form of the builder's bundle map.

    Rows (left: global longitude/latitude map; right: zoom on the impact area):
    (a, b) late-visibility templates (NG-only, RG-only);
    (c, d) effective gravitational slope from the admitted gravity cache;
    (e, f) cliff templates (slope >= 45 deg and >= 60 deg).
    Rows (c-f) are omitted when the bundle has no gravity cache (visibility-only bundles).
    Writes ``<prefix>`` (.eps, .png) and ``<prefix>_caption.txt``; returns the paths.
    """
    plt = style()
    import matplotlib.patheffects as pe
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    masks, window, impact, table = gt["masks"] or {}, gt["window"], gt["impact"], gt["table"]
    meta = bundle_metadata(gt["bundle"]) if gt["bundle"] else {}
    slope, cache_path = gravity_slope(gt["bundle"], gravity_cache, len(F)) if gt["bundle"] else (None, None)
    rows = ["late"] + (["slope", "cliff"] if slope is not None else [])
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    # zoom box: the late track and the late-visibility regions
    lons, lats = [glon[late]], [glat[late]]
    for name in ("late_ng_not_rg", "late_rg_not_ng"):
        if name in masks and masks[name].any():
            lo, la = lonlat(V[F[masks[name]]].mean(axis=1))
            lons.append(lo)
            lats.append(la)
    lons, lats = np.concatenate(lons), np.concatenate(lats)
    pad = 0.10 * max(np.ptp(lons), np.ptp(lats), 1.0)
    zl, za = (lons.min() - pad, lons.max() + pad), (lats.min() - pad, lats.max() + pad)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))

    cmap = plt.get_cmap("magma")
    snorm = Normalize(0.0, 90.0)
    layers = {}
    base = face_colors(len(F), None)
    lc = base.copy()
    if "late_rg_not_ng" in masks:
        lc[masks["late_rg_not_ng"]] = RG_ONLY_COLOR
    if "late_ng_not_rg" in masks:
        lc[masks["late_ng_not_rg"]] = NG_ONLY_COLOR
    layers["late"] = (lc, 1.0)
    if slope is not None:
        layers["slope"] = (cmap(snorm(np.clip(slope, 0, 90)))[:, :3], 0.0)
        cc = base.copy()
        c45 = masks.get("cliff_ge_45deg", slope >= 45.0)
        c60 = masks.get("cliff_ge_60deg", slope >= 60.0)
        cc[c45] = CLIFF45_COLOR
        cc[c60] = CLIFF60_COLOR
        layers["cliff"] = (cc, 1.0)

    nrow = len(rows)
    row_h = 2.05
    height = row_h * nrow + 1.55
    fig = plt.figure(figsize=(DOUBLE_COLUMN_IN, height))
    gs = fig.add_gridspec(nrow, 3, width_ratios=[2.0, 1.0, 0.05], wspace=0.10, hspace=0.30,
                          left=0.07, right=0.93, top=1.0 - 0.10 / nrow, bottom=1.35 / height)
    halo = [pe.Stroke(linewidth=1.5, foreground="white"), pe.Normal()]
    halo_thick = [pe.Stroke(linewidth=2.9, foreground="white"), pe.Normal()]
    letters = iter("abcdef")
    row_titles = {"late": "late-visibility templates",
                  "slope": "effective gravitational slope",
                  "cliff": "cliff templates"}
    for r, row in enumerate(rows):
        colors, shade_w = layers[row]
        for c, (lon_r, lat_r) in enumerate((((-180, 180), (-90, 90)), (zl, za))):
            ax = fig.add_subplot(gs[r, c])
            width = 2 * pixels if c == 0 else pixels
            image, _ = render_lonlat(V, F, N, lon_r, lat_r, width, colors, device_id, shade_weight=shade_w)
            ax.imshow(image, extent=(lon_r[0], lon_r[1], lat_r[0], lat_r[1]), interpolation="lanczos",
                      aspect="equal" if c == 0 else 1.0 / coslat)
            pts = np.column_stack((glon, glat))
            seg = np.stack((pts[:-1], pts[1:]), axis=1)
            keep = gt["outside"][:-1] & gt["outside"][1:] & ok_seg
            from matplotlib.collections import LineCollection
            ax.add_collection(LineCollection(seg[keep], colors="k", linewidths=0.8, zorder=5, path_effects=halo))
            lk = keep & late[:-1] & late[1:]
            ax.add_collection(LineCollection(seg[lk], colors="k", linewidths=2.0, zorder=6, path_effects=halo_thick))
            ax.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=8 if c == 0 else 11,
                    mfc="#F0E442", mec="k", mew=0.6, zorder=9)
            if c == 0:
                ax.plot([glon[0]], [glat[0]], marker="s", ms=3.6, mfc="white", mec="k", mew=0.6, zorder=9)
                ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k",
                        lw=0.6, zorder=7, path_effects=halo)
                ax.set_xticks(range(-180, 181, 60))
                ax.set_yticks(range(-90, 91, 45))
                ax.set_ylabel("latitude [°]")
                ax.set_title(row_titles[row], fontsize=8, loc="left", pad=3)
            else:
                _minute_marks(ax, gt, late, lambda p: lonlat(p), label=False)
                if r == 0:
                    ax.annotate("impact %s UTC" % impact["utc"].split("T")[1],
                                xy=(impact["lon_deg"], impact["lat_deg"]), xytext=(18, -48),
                                textcoords="offset points", fontsize=6.0, ha="center",
                                arrowprops={"arrowstyle": "-|>", "lw": 0.6, "color": "k", "shrinkB": 5},
                                bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.4,
                                      "boxstyle": "round,pad=0.2"}, zorder=10)
            ax.set_xlim(*lon_r)
            ax.set_ylim(*lat_r)
            if r == nrow - 1:
                ax.set_xlabel("longitude [°]")
            ax.text(0.015, 0.96, "(%s)" % next(letters), transform=ax.transAxes, va="top", fontsize=7.5,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8})
        cax = fig.add_subplot(gs[r, 2])
        if row == "slope":
            from matplotlib.cm import ScalarMappable
            cb = fig.colorbar(ScalarMappable(norm=snorm, cmap=cmap), cax=cax, extend="max")
            cb.set_label("slope [°] (> 90: overhang)", fontsize=6.5)
            cb.set_ticks([0, 15, 30, 45, 60, 75, 90])
        else:
            cax.axis("off")
    handles = []
    if "late_ng_not_rg" in masks:
        handles.append(Patch(facecolor=NG_ONLY_COLOR, label="NG-only (%d facets, %.3f%% of area)" % (
            int(masks["late_ng_not_rg"].sum()),
            100 * meta.get("templates", {}).get("late_ng_not_rg", {}).get("area_fraction", np.nan))))
    if "late_rg_not_ng" in masks:
        handles.append(Patch(facecolor=RG_ONLY_COLOR, label="RG-only (%d facets, %.3f%% of area)" % (
            int(masks["late_rg_not_ng"].sum()),
            100 * meta.get("templates", {}).get("late_rg_not_ng", {}).get("area_fraction", np.nan))))
    if slope is not None:
        for name, col, thr in (("cliff_ge_45deg", CLIFF45_COLOR, "45–60°"), ("cliff_ge_60deg", CLIFF60_COLOR, "≥ 60°")):
            frac = meta.get("templates", {}).get(name, {}).get("area_fraction", np.nan)
            handles.append(Patch(facecolor=col, label="slope %s (%s: %.1f%% of area)" % (thr, name, 100 * frac)))
        handles[-2].set_label(handles[-2].get_label().replace("of area)", "of area, incl. ≥ 60°)"))
    handles += [Line2D([], [], color="k", lw=0.8, label="ground track (sub-spacecraft intercept)"),
                Line2D([], [], color="k", lw=2.0, label="ground track in the late-visibility window"),
                Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", mew=0.5,
                       label="1-min marks from %s UTC" % window["start"].split("T")[1]),
                Line2D([], [], marker="s", ls="none", ms=3.6, mfc="white", mec="k", mew=0.6,
                       label="track start %s UTC, %.1f km" % (table["time"][0].strftime("%H:%M:%S"),
                                                             table["altitude"][0] / 1000)),
                Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                       label="Rosetta impact, lon %.2f°, lat %.2f°" % (impact["lon_deg"], impact["lat_deg"]))]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.07, 0.0), ncol=2, fontsize=6.3,
               frameon=False, handlelength=2.0, columnspacing=1.2)
    defn = late_visibility_definition(window, meta)
    fig.text(0.07, 1.0 + 0.02 / nrow,
             "Late visibility: the %d observation points of %s–%s UTC (altitude %s–%s m).\n"
             "NG-only = facets seen by the nude gauge (Kn > 0) and not by the ram gauge (KΦ = KV = 0) from at least "
             "one of these points;\nRG-only = the reverse. 'Seen' = nonzero collisionless transport kernel "
             "(visible and within the gauge acceptance)." % (
                 len(window["points"]), window["start"].split("T")[1], window["stop"].split("T")[1],
                 _hm(window["altitude_max_m"]), _hm(window["altitude_min_m"])),
             fontsize=6.8, va="bottom", ha="left", linespacing=1.35)
    if title:
        fig.suptitle(title, fontsize=8, y=1.0 + 0.06 / nrow)
    paths = save(fig, prefix, overwrite)
    plt.close(fig)
    g = meta.get("gravity", {})
    caption = (
        "Activity-template maps of the %s shape model (%d facets), longitude/latitude, outermost surface "
        "where overhangs put several surface layers at one longitude/latitude. Left: global; right: the "
        "impact area (box in the left panels). (a, b) Late-visibility templates. %s "
        % (Path(gt["shape"]).name, len(F), defn))
    if slope is not None:
        caption += (
            "(c, d) Effective gravitational slope: the angle between the inward facet normal and the effective "
            "gravity (self-gravity of a homogeneous polyhedron, density %s kg/m³, evaluated %s m above each "
            "facet centroid, minus the centrifugal term ω×(ω×r) of the SPICE rotation at %s). Values above "
            "90° are overhangs. (e, f) Cliff templates: slope ≥ 45° (light) and ≥ 60° (dark). Gravity cache: %s. "
            % (g.get("density_kg_m3", "?"), g.get("evaluation_offset_m", "?"),
               ", ".join(g.get("omega_sample_epochs_utc", [])) or "the late interval", cache_path))
    caption += (
        "Line: ground track of Rosetta on 2016-09-30, the sub-spacecraft intercept (first surface hit on the "
        "line from the spacecraft toward the body-frame origin), from %s UTC at %.1f km altitude (square) to "
        "impact; thick within the late-visibility window, with 1-minute marks from its start. Star: Rosetta "
        "impact point, the straight-line extrapolation of the last two observation rows (last row %s UTC, "
        "%.2f m above the surface, %.2f m/s) to the shape model: %s UTC, longitude %.3f°, latitude %.3f° "
        "(planetocentric, body-fixed frame of the SPICE kernels)." % (
            table["time"][0].strftime("%H:%M:%S"), table["altitude"][0] / 1000, table["utc"][-1].split("T")[1],
            table["altitude"][-1], impact["speed_m_s"], impact["utc"].split("T")[1], impact["lon_deg"],
            impact["lat_deg"]))
    paths.append(write_caption(prefix, caption, overwrite))
    return paths


def write_track_table(gt, path, overwrite=False):
    """Ground track table: per observation row, UTC, altitude, ground point, lon/lat, late flag."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(str(path))
    table, impact = gt["table"], gt["impact"]
    lon, lat = lonlat(gt["ground"])
    late = set(gt["window"]["points"])
    with path.open("w", encoding="utf-8") as f:
        f.write("# End-of-mission ground track of Rosetta on %s\n" % gt["shape"])
        f.write("# sub_*: sub-spacecraft intercept (first surface hit toward the body-frame origin; plotted track)\n")
        f.write("# near_*: nearest surface point (the run's altitude definition)\n")
        f.write("# late-visibility window %s .. %s (%s)\n" % (gt["window"]["start"], gt["window"]["stop"],
                                                             gt["window"]["source"]))
        f.write("# impact: %s\n" % json.dumps(impact))
        f.write("#iPoint UTC altitude_m sc_x_m sc_y_m sc_z_m sub_x_m sub_y_m sub_z_m sub_lon_deg sub_lat_deg "
                "near_x_m near_y_m near_z_m near_lon_deg near_lat_deg late\n")
        nlon, nlat = lonlat(gt["nearest"])
        for i in range(len(lon)):
            p = int(table["iPoint"][i])
            f.write("%d %s %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.5f %.5f %.6e %.6e %.6e %.5f %.5f %d\n" % (
                p, table["utc"][i], table["altitude"][i], table["x"][i], table["y"][i], table["z"][i],
                gt["ground"][i, 0], gt["ground"][i, 1], gt["ground"][i, 2], lon[i], lat[i],
                gt["nearest"][i, 0], gt["nearest"][i, 1], gt["nearest"][i, 2], nlon[i], nlat[i], int(p in late)))
    return str(path)

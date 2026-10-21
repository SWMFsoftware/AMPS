"""CPU publication pressure figures, adapted from the CUDA branch figure layout.

No CUDA imports and no fitting or interpolation. Callers supply evaluated pressure
arrays; the standalone reader retains the CPU manifest/provenance checks. This
module intentionally covers fit-result pressures, not the CUDA mesh renderer.
"""

import datetime as dt
import functools
import json
from pathlib import Path

import numpy as np

# c05_write_bundle.py defaults (--late-start/--late-stop); a bundle's metadata takes precedence.
DEFAULT_LATE_WINDOW = ("2016-09-30T10:33:10", "2016-09-30T10:39:26")

# Okabe–Ito colour-blind safe palette.
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
LATE_FILL = "#FBE3C8"          # late-visibility shading (opaque; EPS has no alpha)
LATE_EDGE = "#E69F00"
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
    # Older deployed Matplotlib versions need a draw before tight bounding
    # boxes are measured, especially for a figure-level multi-column legend.
    fig.canvas.draw()
    fig.savefig(str(paths[1]), format="png", dpi=300)
    fig.savefig(str(paths[0]), format="eps", dpi=300)
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
    # Preserve CPU schema/manifest/report hash checks for existing results.
    from plot_cops_model import Run
    verified = Run(str(path))
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
    out = {"path": str(path), "utc": [r[col["UTC"]] for r in rows],
           "provenance": verified.display_name()}
    out["time"] = np.array([parse_utc(u) for u in out["utc"]])
    for name in ("iPoint", "et", "x", "y", "z", "altitude", "pNG", "pRG", "refNG", "refRG",
                 "useNG", "useRG", "inside"):
        if name in col:
            out[name] = np.array([float(r[col[name]]) for r in rows])
    validate_table(out)
    return out


def late_window(table, bundle=None, start=None, stop=None):
    """The late-visibility interval and the altitude range of its observation rows.

    The interval comes from ``start/stop`` if given, else from the bundle's
    ``late_ng_not_rg`` metadata, else ``DEFAULT_LATE_WINDOW``. Rows are selected as in
    ``masks.late_points``: start <= UTC <= stop and not inside the nucleus.
    """
    source = "explicit"
    if (start is None) != (stop is None):
        raise ValueError("late-start and late-stop must be supplied together")
    if start is None:
        start, stop, source = DEFAULT_LATE_WINDOW[0], DEFAULT_LATE_WINDOW[1], "default (c05_write_bundle)"
        if bundle is not None:
            meta = bundle_metadata(bundle)
            interval = meta.get("templates", {}).get("late_ng_not_rg", {}).get("interval_utc")
            if interval:
                start, stop, source = interval[0], interval[1], "bundle %s" % Path(bundle).name
    t0, t1 = parse_utc(start), parse_utc(stop)
    if t0 >= t1:
        raise ValueError("late-start must precede late-stop")
    sel = (table["time"] >= t0) & (table["time"] <= t1)
    if "inside" in table:
        sel &= table["inside"] == 0
    sel &= np.isfinite(table["altitude"]) & (table["altitude"] > 0)
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


def validate_table(table):
    """Fail on malformed rows rather than pairing unrelated observations."""
    n = len(table["time"])
    for key in ("iPoint", "time", "altitude", "pNG", "pRG", "refNG", "refRG", "useNG", "useRG"):
        if np.shape(table[key]) != (n,) or not n:
            raise ValueError("empty or inconsistent column: " + key)
    if len(np.unique(table["iPoint"])) != n:
        raise ValueError("duplicate iPoint")
    if any(b <= a for a, b in zip(table["time"][:-1], table["time"][1:])):
        raise ValueError("UTC rows must be strictly increasing")
    for key in ("useNG", "useRG"):
        if not np.isin(table[key], [0, 1]).all():
            raise ValueError("selection flags must be zero or one")


def aligned(reference, table):
    validate_table(table)
    if not (np.array_equal(reference["iPoint"], table["iPoint"]) and
            np.array_equal(reference["time"], table["time"])):
        raise ValueError("comparison requires matching iPoint/UTC rows; no interpolation is performed")


def positive(values, table):
    """NaN breaks lines at missing/nonphysical points, without filling gaps."""
    values = np.asarray(values, float)
    valid = np.isfinite(values) & (values > 0)
    valid &= np.asarray(table.get("inside", np.zeros(len(values)))) == 0
    return np.where(valid, values, np.nan)


def ratio(ng, rg, table):
    ng, rg = positive(ng, table), positive(rg, table)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return positive(ng / rg, table)


def isolated_style(function):
    """Publication styling must not leak into existing A5 diagnostics."""
    @functools.wraps(function)
    def wrapped(*args, **kwargs):
        import matplotlib
        import matplotlib.pyplot as plt
        existing = set(plt.get_fignums())
        with matplotlib.rc_context():
            try:
                return function(*args, **kwargs)
            finally:
                for number in set(plt.get_fignums()) - existing:
                    plt.close(number)
    return wrapped


def _hm(seconds):
    return "%.0f" % seconds if seconds >= 10 else "%.1f" % seconds


def window_label(window):
    """Shaded-window label; ``window["name"]`` (default "late visibility") names the interval."""
    return "%s\n%s–%s UTC\n%s–%s m" % (
        window.get("name", "late visibility"), window["start"].split("T")[1], window["stop"].split("T")[1],
        _hm(window["altitude_max_m"]), _hm(window["altitude_min_m"]))


# -------------------------------------------------------------------------- pressures
@isolated_style
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
    from cops_paths import output_path
    prefix = output_path(prefix)
    validate_table(ref)
    if not np.isfinite(zoom_minutes) or zoom_minutes <= 0:
        raise ValueError("zoom-minutes must be finite and positive")
    targets = [Path(str(prefix) + suffix + ext) for suffix in ("_time", "_altitude")
               for ext in (".eps", ".png")]
    if not overwrite and any(p.exists() for p in targets):
        raise FileExistsError("figure files exist (use --overwrite)")
    for table in ([base] if base is not None else []) + [s[1] for s in series]:
        aligned(ref, table)
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
    t, alt = ref["time"], positive(ref["altitude"], ref)
    if not np.isfinite(alt).any():
        raise ValueError("no positive outside-surface altitude")
    data = {"NG": (ref["refNG"], ref["useNG"]), "RG": (ref["refRG"], ref["useRG"])}
    written = []

    def panels(ax, x):
        for row, gauge in enumerate(("NG", "RG")):
            obs, use = data[gauge]
            ok = np.isfinite(positive(obs, ref))
            ax[row].plot(x[ok & (use > 0)], obs[ok & (use > 0)], "o", color="k", ms=1.6, mew=0,
                         label="COPS data (fitted)", zorder=3)
            if (ok & (use == 0)).any():
                ax[row].plot(x[ok & (use == 0)], obs[ok & (use == 0)], "o", mfc="none", mec="k", ms=2.2,
                             mew=0.4, label="COPS data (not fitted)", zorder=3)
            for item in lines:
                label, table, ls, color, lw = item
                ax[row].plot(x, positive(table["p" + gauge], table), ls=ls, color=color, lw=lw, label=label, zorder=order[id(item)])
            ax[row].set_yscale("log")
        observed_ratio = ratio(ref["refNG"], ref["refRG"], ref)
        ok = np.isfinite(observed_ratio)
        used = (ref["useNG"] > 0) & (ref["useRG"] > 0)
        ax[2].plot(x[ok & used], observed_ratio[ok & used], "o", color="k", ms=1.6, mew=0, zorder=3)
        ax[2].plot(x[ok & ~used], observed_ratio[ok & ~used], "o", mfc="none", mec="k", ms=2.2, mew=0.4, zorder=3)
        for item in lines:
            label, table, ls, color, lw = item
            ax[2].plot(x, ratio(table["pNG"], table["pRG"], table), ls=ls, color=color, lw=lw, zorder=order[id(item)])
        ax[2].set_yscale("log")

    def shade(a, lo, hi):
        a.axvspan(lo, hi, facecolor=LATE_FILL, edgecolor=LATE_EDGE, lw=0.5, zorder=0)

    def window_title(a, xmid, text):
        tr = mtransforms.blended_transform_factory(a.transData, a.transAxes)
        a.text(xmid, 1.02, text, transform=tr, ha="center", va="bottom", fontsize=6.5, color="0.15",
               bbox={"facecolor": LATE_FILL, "edgecolor": LATE_EDGE, "lw": 0.5, "pad": 1.5})

    ylabels = ("NG pressure [Pa]", "RG pressure [Pa]", "NG / RG")
    span_text = "%s: %s–%s UTC, %s–%s m" % (
        window.get("name", "late visibility"), window["start"].split("T")[1], window["stop"].split("T")[1],
        _hm(window["altitude_max_m"]), _hm(window["altitude_min_m"]))

    # --- time: full descent (left) and the last minutes (right)
    fig, ax = plt.subplots(3, 2, figsize=(DOUBLE_COLUMN_IN, 6.9), sharex="col", sharey="row",
                           gridspec_kw={"width_ratios": [2.2, 1.0], "wspace": 0.06, "hspace": 0.08})
    fig.subplots_adjust(top=0.88)
    # Zoom is anchored to the final observation, not a possibly earlier
    # visibility window. The full panel must never truncate later data.
    tz0 = t[-1] - dt.timedelta(minutes=zoom_minutes)
    tz1 = t[-1] + dt.timedelta(minutes=0.6)
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
    ax[0].set_xlim(np.nanmax(alt) * 1.3, np.nanmin(alt) * 0.7)
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
    # Anchor an axes legend in figure coordinates. This also avoids the
    # figure-legend tight-bbox bug in the host's Matplotlib 3.1 series.
    ax.legend(h2, l2, loc="upper center", bbox_to_anchor=(0.5, y),
              bbox_transform=fig.transFigure, ncol=min(3, len(l2)),
              fontsize=7, handlelength=2.6, columnspacing=1.4)
    if title:
        fig.suptitle(title, y=y + 0.025, fontsize=8)


def evaluated_table(base, utc, point_ids, predictions, references, used):
    """Attach evaluated pressures to audited base geometry, never refit data."""
    table = dict(base, utc=list(utc), time=np.array([parse_utc(str(u)) for u in utc]),
                 iPoint=np.asarray(point_ids))
    for g, name in enumerate(("NG", "RG")):
        table["p" + name] = np.asarray(predictions)[:, g]
        table["ref" + name] = np.asarray(references)[:, g]
        table["use" + name] = np.asarray(used)[:, g]
    aligned(base, table)
    return table


def table_from_rows(rows, columns):
    """In-memory adapter for rows already audited by the running fitter.

    The fit's output manifest is not finalized yet. Reopening its new table
    through the standalone provenance reader would incorrectly require the
    final output hash before the report has even been written.
    """
    utc = [row[columns["UTC"]] for row in rows]
    table = dict(utc=utc, time=np.array([parse_utc(u) for u in utc]))
    for name in ("iPoint", "altitude", "inside", "pNG", "pRG", "refNG", "refRG", "useNG", "useRG"):
        table[name] = np.array([float(row[columns[name]]) for row in rows])
    return table


def read_fine_tuning(directory, base):
    """Read existing fine-tuning pressure products with their manifest hashes.

    These files deliberately are not COPS_model.dat: geometry comes from the
    original base. The base-table hash and exact row identity are both checked.
    No multi-GB kernels or facet fields need to be loaded to draw a figure.
    """
    from cops_io import sha256_file
    directory = Path(directory)
    meta = json.loads((directory / "fine_tuning.json").read_text())
    if meta.get("schema") != "liouville_cops.late_interval_fine_tuning.v1" or meta.get("state") != "complete":
        raise ValueError("fine-tuning result is not a completed supported product")
    path = directory / "pressure_comparison.txt"
    if meta.get("outputs", {}).get(path.name, {}).get("sha256") != sha256_file(path):
        raise ValueError("fine-tuning pressure table does not match manifest")
    if meta["definition"]["table_sha256"] != sha256_file(base["path"]):
        raise ValueError("fine-tuning result requires its original base geometry")
    rows = [line.split() for line in path.read_text().splitlines() if line and not line.startswith("#")]
    if not rows or any(len(r) != 12 for r in rows):
        raise ValueError("malformed fine-tuning pressure table")
    ids = np.array([int(r[0]) for r in rows])
    utc = [r[1] for r in rows]
    refs = np.array([[float(r[3]), float(r[4])] for r in rows])
    used = np.array([[int(r[9]), int(r[10])] for r in rows])
    before, after = [evaluated_table(base, utc, ids,
        np.array([[float(r[k]), float(r[k+1])] for r in rows]), refs, used) for k in (5, 7)]
    for table in (before, after):
        table["saturatedNG"] = np.array([int(r[11]) for r in rows], bool)
        table["provenance"] = "fine tuning; optimizer success=%s; scientific review=%s" % (
            meta["optimizer"]["success"], meta.get("scientific_review", "unknown"))
    return before, after, meta


def fine_interval(end_utc, minutes):
    end = parse_utc(end_utc)
    return (end - dt.timedelta(minutes=minutes)).isoformat(), end.isoformat()


def write_result_figures(directory, base, fit, series=None, interval=None,
                         window_name="late visibility", title=None):
    """Nonfatal visualization hook for a newly created scientific result.

    Failures are explicit in the small status JSON, never reclassified as
    optimizer failure. On success the caller includes all returned relative
    filenames in its authoritative output inventory. A short/test trajectory
    outside the mission window therefore still produces a usable fit.
    """
    from cops_io import file_identity, write_json_atomic
    directory = Path(directory)
    status = {"schema": "liouville_cops.publication_figures.v1", "status": "pending"}
    products = []
    try:
        window = late_window(fit, start=interval[0] if interval else None,
                             stop=interval[1] if interval else None)
        window["name"] = window_name
        paths = pressure_figures(base, series or [("fitted model", fit)], fit, window,
                                 directory / "figures" / "pressures", title=title)
        products = [str(Path(p).relative_to(directory)) for p in paths]
        status.update(status="complete", window={k: v for k, v in window.items() if k not in ("t0", "t1")},
                      outputs={p: file_identity(directory / p) for p in products})
    except Exception as exc:
        status.update(status="failed", error="%s: %s" % (type(exc).__name__, exc))
        print("Publication figures not written: " + status["error"], flush=True)
    status_path = directory / "figures" / "pressures.json"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(status_path, status)
    return products + ["figures/pressures.json"]

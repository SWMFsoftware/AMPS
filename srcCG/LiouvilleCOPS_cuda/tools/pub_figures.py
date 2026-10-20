#!/usr/bin/env python3
"""Publication figures (.eps and .png, 300 dpi) of LiouvilleCOPS runs and fits (CUDA env).

This is the figure format of LiouvilleCOPS_cuda (user request 2026-10-20). The functions are
in ``cops_cuda/figures.py``. ``gpu_fit.py``, ``c05_write_bundle.py`` and ``restart_tool.py``
call them automatically; this tool regenerates or combines figures.

Pressures of one or more fits against the COPS data, with the late-visibility window shaded::

    .venv/bin/python tools/pub_figures.py pressures --base $BASE \
        --fit "N30, 788k=$O/fit_788k_gpu_N30_n2000_r2:dashed" --fit "N30, 1M=$O/fit_001m_gpu_N30_n2000" \
        [--ref <fit dir>] [--bundle <template bundle .npz>] --out <dir>/N30

writes ``<out>_time`` and ``<out>_altitude``.
* ``--fit LABEL=DIR[:dashed|dotted|solid]``: drawn in the given order; earlier fits are drawn
  on top.
* ``--ref``: the fit whose data selection (useNG/useRG) the markers show. Defaults to the
  last ``--fit``.
* ``--bundle``: the late-visibility window is read from its ``late_ng_not_rg`` metadata.
  Without it, ``--late-start/--late-stop`` apply, defaulting to c05_write_bundle's.
* ``--no-base``: omit the unfitted base model.

End-of-mission ground track and the Rosetta impact point on the shape model::

    .venv/bin/python tools/pub_figures.py groundtrack --base $BASE [--shape <BDF>] \
        [--bundle <template bundle .npz>] --out <dir>/groundtrack [--device 0]

writes ``<out>_3d``, ``<out>_latlon`` and ``<out>_track.dat``, plus ``<out>_impact.json``.
* The shape defaults to the base run's ShapeModel.
* With ``--bundle``, the NG-only and RG-only template facets are coloured.

``--overwrite`` replaces existing figure files; by default existing files are refused.
"""

import argparse
import json
from pathlib import Path
import sys

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))

from cops_cuda import figures  # noqa: E402

STYLES = {"solid": "-", "dashed": (0, (3, 1.2)), "dotted": (0, (1, 1.2)), "dashdot": (0, (4, 1.2, 1, 1.2))}


def parse_fit(text):
    label, _, rest = text.partition("=")
    if not rest:
        raise SystemExit("--fit expects LABEL=DIR[:style], got %r" % text)
    path, _, style = rest.rpartition(":") if rest.rsplit(":", 1)[-1] in STYLES else (rest, "", "solid")
    return label, path, STYLES[style or "solid"]


def cmd_pressures(args):
    fits = [parse_fit(f) for f in args.fit]
    base = figures.read_table(args.base)
    series = [(label, figures.read_table(path), style) for label, path, style in fits]
    ref = figures.read_table(args.ref) if args.ref else series[-1][1]
    window = figures.late_window(base, args.bundle, args.late_start, args.late_stop)
    paths = figures.pressure_figures(None if args.no_base else base, series, ref, window, args.out,
                                     title=args.title, overwrite=args.overwrite)
    print("\n".join(paths))
    return 0


def cmd_groundtrack(args):
    base = figures.read_table(args.base)
    window = figures.late_window(base, args.bundle, args.late_start, args.late_stop)
    gt = figures.groundtrack_data(args.base, args.shape, args.bundle, window, args.device)
    paths = figures.groundtrack_figures(gt, args.out, device_id=args.device, overwrite=args.overwrite,
                                        title=args.title)
    impact_path = Path(str(args.out) + "_impact.json")
    if impact_path.exists() and not args.overwrite:
        raise FileExistsError(str(impact_path))
    impact_path.write_text(json.dumps(dict(gt["impact"], shape=gt["shape"], base_run=str(Path(args.base).resolve()),
                                           late_window={k: window[k] for k in ("start", "stop", "source",
                                                                              "altitude_min_m", "altitude_max_m")}),
                                      indent=1), encoding="utf-8")
    paths.append(str(impact_path))
    print("\n".join(paths))
    print("impact: %s UTC, lon %.3f°, lat %.3f° (%s)" % (gt["impact"]["utc"], gt["impact"]["lon_deg"],
                                                      gt["impact"]["lat_deg"], gt["impact"]["method"]))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--base", required=True, help="base run directory (COPS_model.dat, elements.bin)")
    common.add_argument("--bundle", help="template bundle .npz (late window, NG-only/RG-only facets)")
    common.add_argument("--late-start", default=None)
    common.add_argument("--late-stop", default=None)
    common.add_argument("--out", required=True, help="output prefix")
    common.add_argument("--title", default=None)
    common.add_argument("--overwrite", action="store_true")
    p = sub.add_parser("pressures", parents=[common])
    p.add_argument("--fit", action="append", required=True, help="LABEL=DIR[:solid|dashed|dotted|dashdot]")
    p.add_argument("--ref", help="fit directory whose selection the data markers show")
    p.add_argument("--no-base", action="store_true")
    p.set_defaults(func=cmd_pressures)
    g = sub.add_parser("groundtrack", parents=[common])
    g.add_argument("--shape", help="BDF of the base's mesh (default: the base's ShapeModel)")
    g.add_argument("--device", type=int, default=0)
    g.set_defaults(func=cmd_groundtrack)
    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

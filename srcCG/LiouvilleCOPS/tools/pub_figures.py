#!/usr/bin/env python3
"""GPU-style pressure figures from existing CPU fits, without CUDA or refitting.

Run ``pub_figures.py pressures --help``. Result tables and fine-tuning tables
retain their native provenance checks. All output stays outside this application.
The CUDA branch's groundtrack/templates mesh-renderer commands are not provided.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import cops_publication as figures
from cops_io import file_identity, write_json_atomic
from cops_paths import output_path, product_input

STYLES = {"solid": "-", "dashed": (0, (3, 1.2)),
          "dotted": (0, (1, 1.2)), "dashdot": (0, (4, 1.2, 1, 1.2))}


def parse_fit(text):
    label, sep, value = text.partition("=")
    if not sep or not label or not value:
        raise argparse.ArgumentTypeError("expected LABEL=DIR[:solid|dashed|dotted|dashdot]")
    path, colon, suffix = value.rpartition(":")
    return label, product_input(path if colon and suffix in STYLES else value), STYLES.get(suffix, "-")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("pressures", help="NG, RG and NG/RG vs UTC and altitude; PNG+EPS")
    p.add_argument("--base", required=True, help="original base run directory or COPS_model.dat")
    p.add_argument("--fit", action="append", required=True, type=parse_fit,
                   help="LABEL=DIR[:style]; repeat for fits or fine-tuning directories")
    p.add_argument("--ref", help="reference observations/selection; default last --fit")
    p.add_argument("--no-base", action="store_true", help="omit unfitted-base curve")
    p.add_argument("--bundle", help="template bundle .npz/.json defining shaded late window")
    p.add_argument("--late-start", help="explicit window UTC start; requires --late-stop")
    p.add_argument("--late-stop", help="explicit window UTC stop; requires --late-start")
    p.add_argument("--zoom-minutes", type=float, default=25.)
    p.add_argument("--out", required=True, help="output prefix under COPS_OUTPUT_ROOT, or absolute")
    p.add_argument("--title")
    p.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        base = figures.read_table(product_input(args.base))
        inputs = [base["path"]]
        def read(path):
            path = Path(path)
            if (path / "fine_tuning.json").is_file():
                _, table, fine_meta = figures.read_fine_tuning(path, base)
                table["fine_metadata"] = fine_meta
                inputs.extend([str(path / "pressure_comparison.txt"), str(path / "fine_tuning.json")])
                return table
            table = figures.read_table(path)
            inputs.append(table["path"])
            return table

        series = [(label, read(path), style) for label, path, style in args.fit]
        ref = read(product_input(args.ref)) if args.ref else series[-1][1]
        fine_meta = ref.get("fine_metadata")
        start, stop = args.late_start, args.late_stop
        window_name = "late visibility"
        if start is None and stop is None and not args.bundle and fine_meta is not None:
            start, stop = figures.fine_interval(fine_meta["end_utc"], fine_meta["definition"]["minutes"])
            window_name = "fine-tuning window"
        bundle = product_input(args.bundle) if args.bundle else None
        if bundle:
            inputs.append(bundle)
        window = figures.late_window(ref, bundle, start, stop)
        window["name"] = window_name
        prefix = output_path(args.out)
        metadata = Path(str(prefix) + "_figures.json")
        if metadata.exists() and not args.overwrite:
            raise FileExistsError("figure metadata exists (use --overwrite): " + str(metadata))
        paths = figures.pressure_figures(None if args.no_base else base, series, ref, window,
                                        prefix, title=args.title, overwrite=args.overwrite,
                                        zoom_minutes=args.zoom_minutes)
        write_json_atomic(metadata, {
            "schema": "liouville_cops.publication_figures.v1",
            "inputs": [file_identity(path) for path in dict.fromkeys(inputs)],
            "series": [{"label": label, "provenance": table.get("provenance")} for label, table, _ in series],
            "reference_provenance": ref.get("provenance"),
            "altitude_coordinate": "reference table altitude; common iPoint/UTC required",
            "window": {k: v for k, v in window.items() if k not in ("t0", "t1")},
            "zoom_minutes": args.zoom_minutes,
            "outputs": [file_identity(path) for path in paths]})
        print("\n".join(paths + [str(metadata)]))
    except (ValueError, OSError, KeyError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

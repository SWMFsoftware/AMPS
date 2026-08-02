#!/usr/bin/env python3
"""Build the real-measurement C10 POES/MetOp cutoff-boundary reference.

This command consumes a local mirror of the NOAA/NCEI historical SEM-2 Level-2
16-second files.  It never substitutes an empirical or synthetic reference.
The output is split into three products:

``C10_poes_boundary_crossings.csv``
    The primary observational product.  Each row is one inbound or outbound
    satellite-pass crossing of 50 percent of the polar-cap channel flux.

``reference_C10_poes_meped_boundary.csv.gz``
    A two-hour-window, one-hour-step, 3-hour-MLT-bin table used by ``run_C10.py``.
    Empty cells are retained and explicitly marked ``missing=TRUE``.

``C10_reference_manifest.json``
    Source-file SHA-256 values, extraction settings, channel-to-rigidity mapping,
    and row counts.  This file is required for reproducibility and should be
    archived together with the compressed reference CSV.

The December-2006 literature provides the method, fitted model, and figures, but
not the complete pass-level numerical crossing list.  Consequently, a real C10
reference must be regenerated from the NCEI archive; digitizing a paper figure
would not be an equivalent observational reference.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence

from poes_sem2 import (
    aggregate_crossings,
    extract_boundary_crossings,
    load_observations,
    sha256_file,
    write_crossings_csv,
    write_manifest,
    write_reference_csv,
)

SUPPORTED_SUFFIXES = (".txt", ".txt.gz", ".asc", ".asc.gz", ".dat", ".dat.gz", ".cdf")


def parse_utc(text: str) -> datetime:
    """Parse a command-line UTC timestamp and normalize it to UTC."""

    value = text.strip().replace("Z", "+00:00")
    result = datetime.fromisoformat(value)
    if result.tzinfo is None:
        result = result.replace(tzinfo=timezone.utc)
    return result.astimezone(timezone.utc)


def discover_source_files(input_dir: Path) -> List[Path]:
    """Recursively find supported Level-2 files while excluding manifests."""

    files = [
        path for path in input_dir.rglob("*")
        if path.is_file()
        and path.name.lower() != "download_manifest.json"
        and any(path.name.lower().endswith(suffix) for suffix in SUPPORTED_SUFFIXES)
    ]
    return sorted(files)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", type=Path, default=Path("data/reference_source"),
                        help="Directory containing downloaded NCEI Level-2 files")
    parser.add_argument("--event-start", type=parse_utc, default=parse_utc("2006-12-14T00:00:00Z"))
    parser.add_argument("--event-end", type=parse_utc, default=parse_utc("2006-12-16T12:00:00Z"))
    parser.add_argument("--default-altitude-km", type=float, default=850.0,
                        help="Used only when a Level-2 product lacks an altitude column")
    parser.add_argument("--minimum-pass-abs-lat-deg", type=float, default=45.0)
    parser.add_argument("--polar-plateau-abs-lat-deg", type=float, default=75.0)
    parser.add_argument("--minimum-polar-samples", type=int, default=4)
    parser.add_argument("--minimum-plateau-to-low-ratio", type=float, default=2.0)
    parser.add_argument("--minimum-leg-samples", type=int, default=4)
    parser.add_argument("--window-hours", type=float, default=2.0,
                        help="Published Dmitriev-style fit/aggregation window")
    parser.add_argument("--step-hours", type=float, default=1.0,
                        help="Step between adjacent reference-window centers")
    parser.add_argument("--mlt-bin-hours", type=float, default=3.0)
    parser.add_argument("--minimum-crossings-per-cell", type=int, default=1)
    parser.add_argument("--crossings-output", type=Path, default=Path("C10_poes_boundary_crossings.csv"))
    parser.add_argument(
        "--reference-output", type=Path,
        default=Path("reference_C10_poes_meped_boundary.csv.gz"),
        help="gzip-compressed CSV consumed directly by run_C10.py; must end in .csv.gz",
    )
    parser.add_argument("--manifest-output", type=Path, default=Path("C10_reference_manifest.json"))
    parser.add_argument("--summary-output", type=Path, default=Path("C10_reference_summary.json"))
    parser.add_argument("--validate-only", action="store_true",
                        help="Read and extract data but do not overwrite output files")
    args = parser.parse_args(argv)
    if args.event_end <= args.event_start:
        parser.error("--event-end must be later than --event-start")
    if args.mlt_bin_hours <= 0.0 or 24.0 / args.mlt_bin_hours < 1.0:
        parser.error("--mlt-bin-hours must define at least one bin")
    if not args.reference_output.name.lower().endswith(".csv.gz"):
        parser.error(
            "--reference-output must end in .csv.gz; "
            "the production C10 reference is kept gzip-compressed"
        )
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    source_files = discover_source_files(args.input_dir)
    if not source_files:
        print(
            f"No supported Level-2 files were found below {args.input_dir}.\n"
            "Run download_poes_sem2.py first or copy the NCEI daily files into that directory.",
            file=sys.stderr,
        )
        return 2

    print(f"Reading {len(source_files)} NCEI Level-2 files...")
    observations = load_observations(source_files, args.default_altitude_km)
    observations = [row for row in observations if args.event_start <= row.time_utc <= args.event_end]
    if not observations:
        print("No valid observations fall inside the requested event interval.", file=sys.stderr)
        return 2

    crossings = extract_boundary_crossings(
        observations,
        minimum_abs_lat_deg=args.minimum_pass_abs_lat_deg,
        polar_plateau_abs_lat_deg=args.polar_plateau_abs_lat_deg,
        minimum_polar_samples=args.minimum_polar_samples,
        minimum_plateau_to_low_ratio=args.minimum_plateau_to_low_ratio,
        minimum_leg_samples=args.minimum_leg_samples,
    )
    if not crossings:
        print(
            "No boundary crossings were extracted. Inspect the source columns and relax quality "
            "settings only after confirming that the selected days contain an SEP enhancement.",
            file=sys.stderr,
        )
        return 3

    n_bins = max(1, round(24.0 / args.mlt_bin_hours))
    mlt_bins = tuple(round(index * 24.0 / n_bins, 6) for index in range(n_bins))
    cells = aggregate_crossings(
        crossings,
        event_start=args.event_start,
        event_end=args.event_end,
        window_hours=args.window_hours,
        step_hours=args.step_hours,
        mlt_bin_centers=mlt_bins,
        minimum_crossings_per_cell=args.minimum_crossings_per_cell,
    )

    configuration = {
        "event_start_utc": args.event_start.isoformat(),
        "event_end_utc": args.event_end.isoformat(),
        "default_altitude_km": args.default_altitude_km,
        "minimum_pass_abs_lat_deg": args.minimum_pass_abs_lat_deg,
        "polar_plateau_abs_lat_deg": args.polar_plateau_abs_lat_deg,
        "minimum_polar_samples": args.minimum_polar_samples,
        "minimum_plateau_to_low_ratio": args.minimum_plateau_to_low_ratio,
        "minimum_leg_samples": args.minimum_leg_samples,
        "window_hours": args.window_hours,
        "step_hours": args.step_hours,
        "mlt_bins": mlt_bins,
        "minimum_crossings_per_cell": args.minimum_crossings_per_cell,
        "boundary_definition": "50_percent_of_median_flux_at_abs_AACGM_latitude_ge_75_deg",
        "papers_are_not_numerical_reference": True,
    }

    summary = {
        "n_source_files": len(source_files),
        "n_observations": len(observations),
        "n_crossings": len(crossings),
        "n_reference_cells": len(cells),
        "n_nonmissing_reference_cells": sum(not cell.missing for cell in cells),
        "reference_compression": "gzip",
        "satellites": sorted({row.satellite for row in observations}),
        "channels": sorted({row.channel for row in crossings}),
        "first_observation_utc": min(row.time_utc for row in observations).isoformat(),
        "last_observation_utc": max(row.time_utc for row in observations).isoformat(),
    }
    print(json.dumps(summary, indent=2))

    if args.validate_only:
        return 0

    write_crossings_csv(crossings, args.crossings_output)
    manifest_sha = write_manifest(source_files, crossings, cells, configuration, args.manifest_output)
    write_reference_csv(cells, args.reference_output, manifest_sha)
    summary["crossings_output"] = str(args.crossings_output)
    summary["reference_output"] = str(args.reference_output)
    summary["manifest_output"] = str(args.manifest_output)
    summary["manifest_sha256"] = manifest_sha
    summary["reference_sha256"] = sha256_file(args.reference_output)
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"crossings: {args.crossings_output}")
    print(f"reference: {args.reference_output}")
    print(f"manifest:  {args.manifest_output}")
    print(f"summary:   {args.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

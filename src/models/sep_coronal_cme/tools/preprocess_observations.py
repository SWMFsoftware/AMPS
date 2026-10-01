#!/usr/bin/env python3
"""Stage-12 offline preprocessing and immutable campaign assets.

Examples (run from the model directory):
  python3 tools/preprocess_observations.py preprocess --job examples/stage12/job.json --output build/observations
  python3 tools/preprocess_observations.py verify --bundle build/observations
  python3 tools/preprocess_observations.py preregister --bundle build/observations --input preregistration.json --output frozen-preregistration.json
  python3 tools/preprocess_observations.py select --bundle build/observations --preregistration frozen-preregistration.json --input realized-candidates.json --output selection.json

All output files are immutable. No network access, per-cell/particle Python,
implicit coordinate conversion, or SEP-dependent background selection occurs.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
from preprocessing.core import (PreprocessingError, canonical, load_json, read_bundle, verify_frozen)
from preprocessing.inference import preprocess_job
from preprocessing.campaign import preregister_campaign, select_candidates
from preprocessing.protocols import (field_line_requests, freeze_transfer_protocol,
                                     bind_transfer_run, compare_backgrounds)


def write_new(path, record):
    # Exclusive create protects a frozen output against a silent overwrite.
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    data = canonical(record)
    with destination.open("xb") as handle:
        handle.write(data)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("preprocess", help="infer and publish a checksummed source-owned observation bundle")
    p.add_argument("--job", required=True, type=Path); p.add_argument("--output", required=True, type=Path)
    p = commands.add_parser("verify", help="verify complete bundle bytes and metadata/data-use roles")
    p.add_argument("--bundle", required=True, type=Path)
    for name in ("preregister", "select", "field-line-requests"):
        p = commands.add_parser(name)
        p.add_argument("--bundle", required=True, type=Path)
        p.add_argument("--input", required=True, type=Path); p.add_argument("--output", required=True, type=Path)
        if name == "select": p.add_argument("--preregistration", required=True, type=Path)
        if name == "field-line-requests":
            p.add_argument("--ephemeris-asset", required=True); p.add_argument("--reviewed-by", required=True)
    p = commands.add_parser("freeze-transfer")
    p.add_argument("--input", required=True, type=Path); p.add_argument("--output", required=True, type=Path)
    p = commands.add_parser("bind-transfer")
    p.add_argument("--protocol", required=True, type=Path); p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p = commands.add_parser("compare-mhd")
    p.add_argument("--analytic", required=True, type=Path); p.add_argument("--mhd", required=True, type=Path)
    p.add_argument("--allow-combined-discrepancy", action="store_true"); p.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.command == "preprocess":
        record = preprocess_job(load_json(args.job), args.job.resolve().parent, args.output)
    elif args.command == "verify":
        record, _ = read_bundle(args.bundle)
    elif args.command in {"preregister", "select", "field-line-requests"}:
        _, assets = read_bundle(args.bundle)
        if args.command == "preregister": record = preregister_campaign(load_json(args.input), assets)
        elif args.command == "select": record = select_candidates(load_json(args.preregistration), load_json(args.input), assets)
        else:
            matches = [a for a in assets if a["metadata"]["asset_id"] == args.ephemeris_asset]
            if len(matches) != 1: raise PreprocessingError("ephemeris asset missing/ambiguous")
            record = field_line_requests(matches[0], load_json(args.input), args.reviewed_by)
        write_new(args.output, record)
    elif args.command == "freeze-transfer":
        record = freeze_transfer_protocol(load_json(args.input)); write_new(args.output, record)
    elif args.command == "bind-transfer":
        inputs = load_json(args.input)
        record = bind_transfer_run(load_json(args.protocol), inputs["event_inputs"], inputs["release_authorities"])
        write_new(args.output, record)
    else:
        record = compare_backgrounds(load_json(args.analytic), load_json(args.mhd), args.allow_combined_discrepancy)
        write_new(args.output, record)
    verify_frozen(record)
    print("stage12 " + args.command + " PASS identity=" + record["identity"])
    return 0 if record.get("status") != "no-qualified-candidate" else 1

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PreprocessingError, OSError, KeyError, TypeError, ValueError, IndexError) as error:
        print("stage12 ERROR: " + str(error), file=sys.stderr)
        raise SystemExit(2)

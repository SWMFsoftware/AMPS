#!/usr/bin/env python3
"""Audit immutable observations; optionally reproduce bounded archive subsets.

This is reference provenance/integrity software, not a transport simulation or
an observational comparison. It never changes case.json, acceptance criteria,
model forcing or a checked-in reference. Source updates require explicit review.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent

def digest(data):
    return hashlib.sha256(data).hexdigest()

def extract(raw, rule):
    # Keep original header/channel/uncertainty columns and exact numeric bytes.
    # Missing intervals stay absent; no interpolation, averaging or fill repair.
    lines=raw.splitlines(keepends=True)
    if rule["kind"]=="voyager-ascii-event-window":
        count=int(lines[3]);header=lines[:4+count]
        rows=[line for line in lines[4+count:] if rule["date_start_inclusive"]<=line.decode("ascii")[:10]<rule["date_stop_exclusive"]]
        return b"".join(header+rows)
    if rule["kind"]=="ulysses-pha-year-doy":
        header=[];rows=[]
        for line in lines:
            parts=line.split()
            try:
                numeric=len(parts)==11 and 90<=float(parts[0])<110 and 1<=float(parts[1])<367
            except ValueError:
                numeric=False
            if numeric:
                if rule["fractional_year_after1900"]<=float(parts[0])<rule["fractional_year_after1900"]+1 and rule["doy_min_inclusive"]<=float(parts[1])<rule["doy_max_exclusive"]:
                    rows.append(line)
            elif not rows:
                header.append(line)
        return b"".join(header+rows)
    raise ValueError("unknown reviewed extraction rule")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch",action="store_true",help="fetch/check original archive bytes and reproduce into a separate destination")
    parser.add_argument("--source-cache",type=Path,help="offline original source-N.bin files; original SHA256 checks remain mandatory")
    parser.add_argument("--output",type=Path,help="fresh output directory required with --fetch or --source-cache")
    args=parser.parse_args()
    case=json.loads((ROOT/"case.json").read_text())
    for artifact in case["artifacts"]:
        data=(ROOT/artifact["path"]).read_bytes()
        if len(data)!=artifact["bytes"] or digest(data)!=artifact["sha256"]:
            raise ValueError("reference artifact integrity mismatch: "+artifact["path"])
    if case["comparison"]["status"]!="NOT_RUN" or case["comparison"]["validation_claim"]:
        raise ValueError("an integrity audit cannot promote observational qualification")
    if args.fetch and args.source_cache:
        raise ValueError("choose live fetch or hash-enforced source cache")
    if args.fetch or args.source_cache:
        if args.output is None or args.output.exists():
            raise ValueError("reproduction requires a fresh separate --output directory")
        output=args.output.resolve();output.mkdir(parents=True)
        retrieval=json.loads((ROOT/"retrieval.json").read_text())
        for index, record in enumerate(retrieval["downloads"]):
            raw_path=(args.source_cache if args.source_cache else output)/("source-"+str(index)+".bin")
            # curl is invoked as argv with a strict2MB bound. This archive
            # audit never downloads bulk ERNE exports or entire mission data.
            if not args.source_cache:
                subprocess.run(["curl","-fsSL","--max-time","45","--max-filesize","2000000",record["url"],"-o",str(raw_path)],check=True)
            raw=raw_path.read_bytes()
            expected=record.get("source_sha256",record.get("sha256"))
            if digest(raw)!=expected:
                raise ValueError("upstream source changed; retain new bytes for review, do not update reference automatically")
            if "extract" in record:
                generated=extract(raw,record["extract"])
                if digest(generated)!=record["output_sha256"]:
                    raise ValueError("reproduced subset differs from reviewed reference")
                (output/record["extract"]["path"]).write_bytes(generated)
    print(case["case_id"]+": reference integrity verified; observational comparison NOT_RUN; model readiness "+case["readiness"]["model_inputs"])
    return 0

if __name__=="__main__":
    try:
        raise SystemExit(main())
    except (ValueError,OSError,subprocess.CalledProcessError) as error:
        print("ERROR: "+str(error),file=sys.stderr)
        raise SystemExit(2)

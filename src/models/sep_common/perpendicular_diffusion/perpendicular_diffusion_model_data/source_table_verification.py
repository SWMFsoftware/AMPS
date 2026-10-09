#!/usr/bin/env python3
"""Compare all 576 CSV coefficients with primary Qin-Zhang Tables 3/4.

Retrieve the source PDF listed in NLGCE_source_audit.json, and extract it with
`pdftotext -layout source.pdf source.txt`. Pass that text as --primary-text.
The script reads only local files and never regenerates missing coefficients.
"""
import argparse
import csv
import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def verify(primary_text,primary_pdf=None):
    text = primary_text.read_text()
    records = []
    for kind,number,next_number in [("parallel",3,4),("perpendicular",4,None)]:
        section = text.split(f"Table {number}:")[1]
        if next_number is not None:
            section = section.split(f"Table {next_number}:")[0]
        source = {}
        for line in section.splitlines():
            match = re.match(r"^\s*\((\d)\s+(\d)\s+(\d)\)\s+(.+)$",line)
            if match:
                j,k,l,values = match.groups()
                key = (int(j),int(k),int(l))
                if key in source:
                    raise ValueError(f"Duplicate source row {key}")
                tokens = values.split()
                if len(tokens)!=6:
                    raise ValueError(f"Incorrect source row length {key}")
                source[key] = [Decimal(value) for value in tokens]
        expected = {(j,k,l) for j in range(4) for k in range(4) for l in range(3)}
        if set(source)!=expected:
            raise AssertionError("Missing/out-of-range primary coefficient rows")
        path = ROOT/f"NLGCE_F_2014_{kind}.csv"
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows)!=48:
            raise AssertionError("Incorrect CSV row count")
        seen = set()
        for row in rows:
            key = tuple(int(row[n]) for n in ("j","k","l"))
            if key in seen:
                raise AssertionError("Duplicate CSV tuple")
            seen.add(key)
            actual = [Decimal(row[f"d_i{i}"]) for i in range(6)]
            if actual!=source[key]:
                raise AssertionError(f"Coefficient mismatch {kind}, {key}")
        if seen!=expected:
            raise AssertionError("CSV index set mismatch")
        records.append({"asset":path.name,"source_table":number,"rows":48,"coefficients":288,
                        "sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
    if primary_pdf is not None:
        expected_hash = json.loads((ROOT/"NLGCE_source_audit.json").read_text())["source_pdf_sha256"]
        if hashlib.sha256(primary_pdf.read_bytes()).hexdigest()!=expected_hash:
            raise AssertionError("Source PDF bytes differ from the audited version")
    return {"state":"passed","all_coefficients_compared_by_explicit_index":576,
            "source_pdf_hash_compared":primary_pdf is not None,"tables":records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-text",required=True,type=Path)
    parser.add_argument("--primary-pdf",type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.primary_text,args.primary_pdf),indent=2))


if __name__ == "__main__":
    main()

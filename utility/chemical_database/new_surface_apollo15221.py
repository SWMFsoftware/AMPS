#!/usr/bin/env python3
"""Lossless CSV views of the four Zenodo5044539 Apollo15221 source tables.

The source README declares temperature in kelvin and desorption in monolayers
per second. Native energy headers specify J/mol. ``Probability`` is retained
as the producer's label: this reader does not decide whether it is a mass,
density, or normalized distribution. TPD files mix measured and model columns;
their columns are never pooled or converted into a universal surface law.

All numerical fields remain strings, including signed values and repeated
coverage coordinates. A successful projection means lexical access only. It
does not qualify desorption order, surface site distributions, accommodation,
sticking, an Arrhenius prefactor, or use on any other material/specimen.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


SCHEMAS = {
    "hydrogen-tpd": ["Temperature (K)", "Void TPD", "Vacuum TPD", "Sum", "Experimental TPD"],
    "water-tpd": ["Temperature (K)", "Void TPD", "Vacuum TPD", "Sum", "Experimental TPD"],
    "coverage-energy": ["Fractional Coverage", "Energy (J/mol)"],
    "site-distribution": ["Energy (J/mol)", "Probability"],
}


def sha256(path):
    """Hash the exact source/recipe bytes; never normalize their text first."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def decode(raw, kind):
    """Read only the explicitly selected native TSV schema, transactionally.

    Tab boundaries matter: whitespace splitting would erase empty fields and
    move values between measured/model columns. Empty represented fields are
    preserved. A changed header, ragged row or wholly missing data is an ERROR,
    rather than an opportunity to infer a column or silently skip a source row.
    """
    if kind not in SCHEMAS:
        raise ValueError("unknown explicit Apollo15221 table kind")
    if not isinstance(raw, bytes):
        raise TypeError("native table must be bytes")
    text = raw.decode("utf-8")
    rows = list(csv.reader(io.StringIO(text, newline=""), delimiter="\t", strict=True))
    if not rows or rows[0] != SCHEMAS[kind]:
        raise ValueError("native header changed; no column/unit inference permitted")
    if len(rows) < 2:
        raise ValueError("native table has no data rows")
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError("native table has a ragged or blank source row")
    return rows


def project(source, output, *, kind, expected_sha256):
    """Authenticate then export all tokens; refuse replacement of old evidence.

    The caller supplies a reviewed source SHA and a NEW destination. All input
    validation completes before output publication. Failed IO may leave partial
    evidence for diagnosis; it is never reused automatically as success.
    """
    source, output = Path(source), Path(output)
    if output.exists():
        raise FileExistsError("CSV output already exists")
    raw = source.read_bytes()
    actual_sha = hashlib.sha256(raw).hexdigest()
    if actual_sha != expected_sha256:
        raise ValueError("source SHA256 mismatch")
    rows = decode(raw, kind)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(rows)
    with output.open(encoding="utf-8", newline="") as stream:
        if list(csv.reader(stream, strict=True)) != rows:
            raise ValueError("CSV lexical round-trip changed native fields")
    return dict(schema="amps-apollo15221-lexical-projection-v1", status="PASS",
                scope="source lexical preservation only", physical_activation=False,
                source_sha256=actual_sha, recipe_sha256=sha256(__file__),
                kind=kind, columns=rows[0], row_count=len(rows)-1,
                output_sha256=sha256(output), original_source_modified=False)


def main(argv=None):
    """Return 0/PASS, 77/missing prerequisite, or 2/parser/IO ERROR honestly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-file", type=Path, required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--kind", choices=sorted(SCHEMAS), required=True)
    parser.add_argument("--output-file", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = project(args.source_file, args.output_file, kind=args.kind,
                         expected_sha256=args.source_sha256)
    except FileNotFoundError as error:
        print(json.dumps(dict(status="SKIPPED", reason=str(error)))); return 77
    except (OSError, ValueError, TypeError, csv.Error) as error:
        print(json.dumps(dict(status="ERROR", reason=str(error)))); return 2
    print(json.dumps(result, indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Validate frozen data, input decks, and vendored observation operators."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

from study_common import load_config, read_driver, sha256


def count_csv(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        content = (line for line in stream if not line.startswith("#"))
        return sum(1 for _ in csv.DictReader(content))


def required_columns(path: Path, names: set[str]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(line for line in stream if not line.startswith("#"))
        missing = names - set(reader.fieldnames or ())
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")


def validate_input(path: Path) -> None:
    """Check parser-critical directives and reject previously problematic keys."""
    text = path.read_text(encoding="utf-8")
    required = {
        "CALC_TARGET": "CUTOFF_RIGIDITY",
        "CUTOFF_SAMPLING": "VERTICAL",
        "CUTOFF_SEARCH_ALGORITHM": "RIGIDITY_LIST",
        "FIELD_MODEL": "T05",
        "OUTPUT_MODE": "SHELLS",
        "SHELL_GEOMETRY": "GEODETIC",
    }
    active: Dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(("!", "#")):
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            active[parts[0]] = parts[1].strip()
    for key, value in required.items():
        if active.get(key) != value:
            raise ValueError(f"{path}: expected {key} {value}, found {active.get(key)!r}")
    forbidden = {
        "CUTOFF_UNRESOLVED_EXTENSION_PASSES",
        "CUTOFF_UNRESOLVED_EXTENSION_FACTOR",
    }
    present = forbidden.intersection(active)
    if present:
        raise ValueError(f"{path}: unsupported parser keyword(s): {sorted(present)}")


def run_check(command: List[str], cwd: Path) -> None:
    completed = subprocess.run(command, cwd=cwd, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if completed.returncode != 0:
        raise ValueError(f"command failed: {' '.join(command)}\n{completed.stdout}")
    print(completed.stdout.strip())


def main() -> int:
    root, config = load_config()
    provenance = json.loads((root / "data" / "provenance.json").read_text())
    paths = {
        "pamela_table_s1.csv": root / config["data"]["pamela_reference"],
        "poes_metop_meped_boundaries.csv.gz": root / config["data"]["poes_reference"],
        "ts05_dec2006_5min.txt": root / config["data"]["driver"],
    }
    problems: List[str] = []
    for name, path in paths.items():
        expected = provenance[name]
        try:
            actual = sha256(path)
            if actual != expected["sha256"]:
                raise ValueError(f"digest {actual} != {expected['sha256']}")
            if name.endswith((".csv", ".csv.gz")):
                rows = count_csv(path)
            else:
                rows = len(read_driver(path))
            if rows != int(expected["rows"]):
                raise ValueError(f"row count {rows} != {expected['rows']}")
            print(f"PASS {name}: {rows} rows, sha256={actual}")
        except Exception as exc:
            problems.append(f"{name}: {exc}")

    try:
        required_columns(paths["pamela_table_s1.csv"], {
            "interval_midpoint_utc", "rigidity_geometric_center_gv",
            "pamela_cutoff_aacgm_deg", "sigma_plus_deg", "sigma_minus_deg",
        })
        required_columns(paths["poes_metop_meped_boundaries.csv.gz"], {
            "interval_midpoint_utc", "rigidity_gv", "channel", "hemisphere",
            "mlt_hour", "boundary_aacgm_lat_deg", "validation_role",
            "acceptance_eligible",
        })
    except Exception as exc:
        problems.append(str(exc))

    for altitude in (475, 850):
        try:
            validate_input(root / "inputs" / f"AMPS_PARAM_DEC2006_{altitude}km.in")
            print(f"PASS AMPS_PARAM_DEC2006_{altitude}km.in")
        except Exception as exc:
            problems.append(str(exc))

    # Exercise the original validators as a defense against a study-level check
    # accidentally becoming less strict than C9 or C10.
    try:
        run_check([sys.executable, "run_C9.py", "--validate-references"], root / "vendor" / "C9")
        run_check([sys.executable, "run_C9.py", "--validate-driver"], root / "vendor" / "C9")
        run_check([sys.executable, "run_C10.py", "--validate-references"], root / "vendor" / "C10")
        run_check([sys.executable, "run_C10.py", "--validate-driver"], root / "vendor" / "C10")
        run_check([sys.executable, "run_C10.py", "--self-test"], root / "vendor" / "C10")
    except Exception as exc:
        problems.append(str(exc))

    if problems:
        print("PACKAGE VALIDATION FAILED", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        return 1
    print("PACKAGE VALIDATION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

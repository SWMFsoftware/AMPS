#!/usr/bin/env python3
"""Compare a live SWMF cutoff/access artifact with replay of the same snapshot.

This is the Step-10 cross-path comparator, not a replacement for C/F physics gates.
Both inputs must come from the Mode3D production writer and must declare the same
snapshot ID, mesh revision, content fingerprint, variable schema, and row topology.
Numeric comparison defaults to exact equality because live and replay use the same
solver and identical frozen compact field.  A validation campaign may supply an
explicit, preregistered tolerance; the runner never silently widens it.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


AUX_RE = re.compile(r'^AUXDATA\s+([A-Za-z0-9_]+)\s*=\s*"(.*)"\s*$')
QUOTED_RE = re.compile(r'"([^"]+)"')
REQUIRED_PROVENANCE = (
    "SNAPSHOT_ID",
    "SNAPSHOT_EPOCH_UTC",
    "SNAPSHOT_MESH_REVISION",
    "SNAPSHOT_CONTENT_FINGERPRINT",
    "OUTER_BOUNDARY_POLICY",
)
SHUE_PROVENANCE = (
    "SHUE_R0_RE",
    "SHUE_ALPHA",
    "SHUE_TAIL_CAP_X_M",
)


class ComparisonFailure(RuntimeError):
    """Raised for a scientific/schema mismatch, not a runner implementation error."""


@dataclass
class Product:
    path: Path
    variables: List[str]
    aux: Dict[str, str]
    rows: List[List[str]]
    zones: List[str]


def parse_product(path: Path) -> Product:
    """Read the append-only Tecplot contract without assuming a fixed product width."""
    if not path.is_file():
        raise ComparisonFailure(f"artifact does not exist: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    variables: List[str] = []
    aux: Dict[str, str] = {}
    rows: List[List[str]] = []
    zones: List[str] = []
    collecting_variables = False

    for line_number, raw in enumerate(lines, start=1):
        text = raw.strip()
        if not text or text.startswith("#"):
            continue
        upper = text.upper()
        if upper.startswith("TITLE="):
            collecting_variables = False
            continue
        if upper.startswith("VARIABLES="):
            collecting_variables = True
            variables.extend(QUOTED_RE.findall(text.split("=", 1)[1]))
            continue
        if upper.startswith("ZONE"):
            collecting_variables = False
            zones.append(text)
            continue
        match = AUX_RE.match(text)
        if match:
            collecting_variables = False
            key, value = match.groups()
            if key in aux and aux[key] != value:
                raise ComparisonFailure(
                    f"{path}:{line_number}: conflicting AUXDATA {key}"
                )
            aux[key] = value
            continue
        if collecting_variables:
            names = QUOTED_RE.findall(text)
            if names:
                variables.extend(names)
                continue
            collecting_variables = False

        try:
            tokens = shlex.split(text, comments=False, posix=True)
        except ValueError as error:
            raise ComparisonFailure(f"{path}:{line_number}: {error}") from error
        if tokens:
            rows.append(tokens)

    if not variables:
        raise ComparisonFailure(f"{path}: missing VARIABLES declaration")
    if len(set(variables)) != len(variables):
        raise ComparisonFailure(f"{path}: VARIABLES contains duplicate names")
    if not zones:
        raise ComparisonFailure(f"{path}: missing ZONE declaration")
    if not rows:
        raise ComparisonFailure(f"{path}: contains no product rows")
    for row_number, row in enumerate(rows, start=1):
        if len(row) != len(variables):
            raise ComparisonFailure(
                f"{path}: row {row_number} has {len(row)} columns; "
                f"VARIABLES declares {len(variables)}"
            )
    return Product(path=path, variables=variables, aux=aux, rows=rows, zones=zones)


def as_float(token: str) -> Tuple[bool, float]:
    try:
        return True, float(token)
    except ValueError:
        return False, 0.0


def compare_products(live: Product, replay: Product, *, rtol: float, atol: float) -> dict:
    if not (math.isfinite(rtol) and rtol >= 0.0 and math.isfinite(atol) and atol >= 0.0):
        raise ComparisonFailure("rtol and atol must be finite and non-negative")
    if live.variables != replay.variables:
        raise ComparisonFailure("VARIABLES schema differs between live and replay")
    if live.zones != replay.zones:
        raise ComparisonFailure("ZONE declarations/topology differ between live and replay")
    if len(live.rows) != len(replay.rows):
        raise ComparisonFailure("row count differs between live and replay")

    for key in REQUIRED_PROVENANCE:
        if not live.aux.get(key):
            raise ComparisonFailure(f"live artifact is missing AUXDATA {key}")
        if not replay.aux.get(key):
            raise ComparisonFailure(f"replay artifact is missing AUXDATA {key}")
        if live.aux[key] != replay.aux[key]:
            raise ComparisonFailure(f"provenance mismatch for {key}")

    # A common snapshot ID is not enough when the physical escape surface is supplied
    # by input rather than by the field CSV. For SHUE products require the fully
    # resolved surface on both paths, so AUTO resolution or a tail-cap mismatch fails
    # before numeric rows are considered.
    if live.aux["OUTER_BOUNDARY_POLICY"].upper() == "SHUE":
        for key in SHUE_PROVENANCE:
            if not live.aux.get(key):
                raise ComparisonFailure(f"live SHUE artifact is missing AUXDATA {key}")
            if not replay.aux.get(key):
                raise ComparisonFailure(f"replay SHUE artifact is missing AUXDATA {key}")
            if live.aux[key] != replay.aux[key]:
                raise ComparisonFailure(f"provenance mismatch for {key}")

    max_abs = 0.0
    max_rel = 0.0
    numeric_values = 0
    for row_index, (left, right) in enumerate(zip(live.rows, replay.rows), start=1):
        for column_index, (a_token, b_token) in enumerate(zip(left, right), start=1):
            a_is_number, a = as_float(a_token)
            b_is_number, b = as_float(b_token)
            if a_is_number != b_is_number:
                raise ComparisonFailure(
                    f"row {row_index} column {column_index}: token type differs"
                )
            if not a_is_number:
                if a_token != b_token:
                    raise ComparisonFailure(
                        f"row {row_index} column {column_index}: text differs"
                    )
                continue
            numeric_values += 1
            if math.isnan(a) or math.isnan(b):
                if not (math.isnan(a) and math.isnan(b)):
                    raise ComparisonFailure(
                        f"row {row_index} column {column_index}: NaN mismatch"
                    )
                continue
            if math.isinf(a) or math.isinf(b):
                if a != b:
                    raise ComparisonFailure(
                        f"row {row_index} column {column_index}: infinity mismatch"
                    )
                continue
            difference = abs(a - b)
            scale = max(abs(a), abs(b))
            relative = 0.0 if scale == 0.0 else difference / scale
            max_abs = max(max_abs, difference)
            max_rel = max(max_rel, relative)
            if difference > atol + rtol * scale:
                variable = live.variables[column_index - 1]
                raise ComparisonFailure(
                    f"row {row_index} variable {variable!r}: live={a:.17g}, "
                    f"replay={b:.17g}, abs={difference:.3g} exceeds "
                    f"atol+rtol*scale={atol + rtol * scale:.3g}"
                )

    return {
        "schema": "sep-in-geospace/swmf-coupled-access-comparison/v1",
        "RESULT": "PASS",
        "snapshot_id": live.aux["SNAPSHOT_ID"],
        "content_fingerprint": live.aux["SNAPSHOT_CONTENT_FINGERPRINT"],
        "mesh_revision": live.aux["SNAPSHOT_MESH_REVISION"],
        "epoch_utc": live.aux["SNAPSHOT_EPOCH_UTC"],
        "outer_boundary_policy": live.aux["OUTER_BOUNDARY_POLICY"],
        "rows": len(live.rows),
        "columns": len(live.variables),
        "numeric_values": numeric_values,
        "rtol": rtol,
        "atol": atol,
        "max_abs_difference": max_abs,
        "max_relative_difference": max_rel,
        "live": str(live.path),
        "replay": str(replay.path),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare coupled cutoff/access output with identical-snapshot replay."
    )
    parser.add_argument("--live", required=True, type=Path,
                        help="Tecplot artifact written during live SWMF coupling")
    parser.add_argument("--replay", required=True, type=Path,
                        help="corresponding standalone SWMF_SNAPSHOT replay artifact")
    parser.add_argument("--rtol", type=float, default=0.0,
                        help="explicit numeric relative tolerance (default: exact)")
    parser.add_argument("--atol", type=float, default=0.0,
                        help="explicit numeric absolute tolerance (default: exact)")
    parser.add_argument("--json", type=Path,
                        help="optional machine-readable comparison summary")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = compare_products(
            parse_product(args.live), parse_product(args.replay),
            rtol=args.rtol, atol=args.atol,
        )
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
        print(
            f"PASS: identical snapshot {summary['snapshot_id']}; "
            f"{summary['rows']} rows, max_abs={summary['max_abs_difference']:.3g}"
        )
        print("RESULT: PASS")
        return 0
    except (OSError, ComparisonFailure) as error:
        failure = {"RESULT": "FAIL", "message": str(error),
                   "live": str(args.live), "replay": str(args.replay)}
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
        print(f"FAIL: {error}", file=sys.stderr)
        print("RESULT: FAIL", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

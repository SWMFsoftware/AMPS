#!/usr/bin/env python3
"""Compare Step-12 standalone/mesh and live/replay products row by row.

The input is the ``cross_path_campaigns`` section of a Step-12 release manifest.
Every CSV byte stream is SHA-256 pinned; every column is declared as a key, an exact
identity/outcome, or a finite numeric observable.  This avoids permissive comparisons
that accidentally ignore a new product column or treat NaN as agreement.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from release_validation.contract import (  # type: ignore
        COMMON_PHYSICS_TAG,
        PHASE_1_SCOPE,
        REQUIRED_CAMPAIGN_KINDS,
        REQUIRED_PARITY_ROLES,
        TOLERANCE_CEILINGS,
        ReleaseValidationError,
        ResourceResolver,
        atomic_write_json,
        manifest_sha256,
        read_json,
        require,
        require_mapping,
        require_nonempty_string,
        require_sequence,
        require_sha256,
        validate_release_manifest_shape,
    )
else:
    from .contract import (
        COMMON_PHYSICS_TAG,
        PHASE_1_SCOPE,
        REQUIRED_CAMPAIGN_KINDS,
        REQUIRED_PARITY_ROLES,
        TOLERANCE_CEILINGS,
        ReleaseValidationError,
        ResourceResolver,
        atomic_write_json,
        manifest_sha256,
        read_json,
        require,
        require_mapping,
        require_nonempty_string,
        require_sequence,
        require_sha256,
        validate_release_manifest_shape,
    )


COMPARISON_SCHEMA = "sep-in-geospace/cross-path-comparison/v1"


def _read_csv(path: Path, label: str) -> Tuple[List[str], List[Dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None, "%s has no CSV header" % label)
            header = [str(name).strip() for name in reader.fieldnames]
            require(all(header), "%s contains an empty CSV column name" % label)
            require(len(header) == len(set(header)), "%s contains duplicate CSV columns" % label)
            rows = []
            for row_number, raw in enumerate(reader, 2):
                require(None not in raw, "%s row %d has too many fields" % (label, row_number))
                row = {name: str(raw.get(name, "")) for name in header}
                rows.append(row)
    except (OSError, csv.Error) as exc:
        raise ReleaseValidationError("cannot parse %s %s: %s" % (label, path, exc)) from exc
    require(bool(rows), "%s contains no comparison rows" % label)
    return header, rows


def _column_list(value: Any, label: str, allow_empty: bool = False) -> List[str]:
    raw = require_sequence(value, label)
    columns = [require_nonempty_string(item, "%s[%d]" % (label, index)) for index, item in enumerate(raw)]
    require(allow_empty or bool(columns), "%s must not be empty" % label)
    require(len(columns) == len(set(columns)), "%s contains duplicate columns" % label)
    return columns


def _indexed_rows(
    rows: Sequence[Mapping[str, str]], key_columns: Sequence[str], label: str
) -> Dict[Tuple[str, ...], Mapping[str, str]]:
    output: Dict[Tuple[str, ...], Mapping[str, str]] = {}
    for number, row in enumerate(rows, 2):
        key = tuple(row[column] for column in key_columns)
        require(all(part != "" for part in key), "%s row %d has an empty key" % (label, number))
        require(key not in output, "%s contains duplicate key %r" % (label, key))
        output[key] = row
    return output


def _finite_number(token: str, label: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise ReleaseValidationError("%s is not numeric: %r" % (label, token)) from exc
    require(math.isfinite(value), "%s must be finite; NaN/Inf cannot pass parity" % label)
    return value


def _validate_tolerance(pair: Mapping[str, Any], campaign_kind: str, role: str) -> Tuple[str, float, float]:
    tolerance_class = require_nonempty_string(pair.get("tolerance_class"), "pair.tolerance_class")
    require(tolerance_class in TOLERANCE_CEILINGS, "unknown tolerance class %s" % tolerance_class)
    rtol = pair.get("rtol")
    atol = pair.get("atol")
    require(isinstance(rtol, (int, float)) and math.isfinite(float(rtol)) and float(rtol) >= 0.0,
            "%s rtol is invalid" % role)
    require(isinstance(atol, (int, float)) and math.isfinite(float(atol)) and float(atol) >= 0.0,
            "%s atol is invalid" % role)
    rtol = float(rtol)
    atol = float(atol)
    ceiling_rtol, ceiling_atol = TOLERANCE_CEILINGS[tolerance_class]
    require(rtol <= ceiling_rtol and atol <= ceiling_atol,
            "%s relaxes the %s ceiling (rtol %.17g, atol %.17g)" % (role, tolerance_class, ceiling_rtol, ceiling_atol))

    # A frozen SWMF snapshot replay has the same cell values and must reproduce the
    # live calculation to kernel precision.  The sampled analytic-field campaign is
    # permitted the roadmap's documented five-percent mesh representation envelope.
    if campaign_kind == "SWMF_OFFLINE_REPLAY":
        require(rtol <= TOLERANCE_CEILINGS["KERNEL"][0] and atol <= TOLERANCE_CEILINGS["KERNEL"][1],
                "%s SWMF replay comparison exceeds kernel precision" % role)
    else:
        maximum = TOLERANCE_CEILINGS["DETECTOR" if role == "DETECTOR" else "MESH"]
        require(rtol <= maximum[0] and atol <= maximum[1],
                "%s sampled-mesh comparison exceeds the roadmap ceiling" % role)
    return tolerance_class, rtol, atol


def compare_pair(
    pair: Mapping[str, Any], campaign_kind: str, resolver: ResourceResolver
) -> Dict[str, Any]:
    role = require_nonempty_string(pair.get("role"), "pair.role")
    require(role in REQUIRED_PARITY_ROLES, "unknown cross-path role %s" % role)
    reference = resolver.resolve(pair.get("reference"), "%s reference" % role)
    candidate = resolver.resolve(pair.get("candidate"), "%s candidate" % role)
    require(reference != candidate, "%s reference and candidate are the same file" % role)
    tolerance_class, rtol, atol = _validate_tolerance(pair, campaign_kind, role)

    reference_header, reference_rows = _read_csv(reference, "%s reference" % role)
    candidate_header, candidate_rows = _read_csv(candidate, "%s candidate" % role)
    require(reference_header == candidate_header, "%s CSV schemas differ" % role)
    key_columns = _column_list(pair.get("key_columns"), "%s.key_columns" % role)
    exact_columns = _column_list(pair.get("exact_columns"), "%s.exact_columns" % role)
    numeric_columns = _column_list(pair.get("numeric_columns"), "%s.numeric_columns" % role)
    declared = key_columns + exact_columns + numeric_columns
    require(len(declared) == len(set(declared)), "%s column groups overlap" % role)
    require(set(declared) == set(reference_header),
            "%s must classify every CSV column exactly once" % role)

    expected = _indexed_rows(reference_rows, key_columns, "%s reference" % role)
    actual = _indexed_rows(candidate_rows, key_columns, "%s candidate" % role)
    missing = sorted(set(expected) - set(actual))
    extra = sorted(set(actual) - set(expected))
    failures: List[Dict[str, Any]] = []
    if missing:
        failures.append({"kind": "missing_keys", "count": len(missing), "sample": missing[:5]})
    if extra:
        failures.append({"kind": "extra_keys", "count": len(extra), "sample": extra[:5]})

    max_abs = 0.0
    max_rel = 0.0
    compared_values = 0
    for key in sorted(set(expected) & set(actual)):
        left = expected[key]
        right = actual[key]
        for column in exact_columns:
            if left[column] != right[column] and len(failures) < 25:
                failures.append({"kind": "exact", "key": key, "column": column,
                                 "reference": left[column], "candidate": right[column]})
        for column in numeric_columns:
            a = _finite_number(left[column], "%s reference %r %s" % (role, key, column))
            b = _finite_number(right[column], "%s candidate %r %s" % (role, key, column))
            difference = abs(a - b)
            scale = max(abs(a), abs(b))
            relative = difference / scale if scale > 0.0 else 0.0
            max_abs = max(max_abs, difference)
            max_rel = max(max_rel, relative)
            compared_values += 1
            if difference > atol + rtol * scale and len(failures) < 25:
                failures.append({"kind": "numeric", "key": key, "column": column,
                                 "reference": a, "candidate": b, "absolute_error": difference,
                                 "relative_error": relative,
                                 "allowed_error": atol + rtol * scale})

    return {
        "role": role,
        "status": "PASS" if not failures else "FAIL",
        "reference": str(reference),
        "candidate": str(candidate),
        "row_count": len(expected),
        "numeric_value_count": compared_values,
        "tolerance_class": tolerance_class,
        "rtol": rtol,
        "atol": atol,
        "max_absolute_error": max_abs,
        "max_relative_error": max_rel,
        "failures": failures,
    }


def compare_campaigns(
    manifest: Mapping[str, Any], resolver: ResourceResolver
) -> Dict[str, Any]:
    campaigns = require_sequence(manifest.get("cross_path_campaigns"), "cross_path_campaigns")
    by_kind: Dict[str, Mapping[str, Any]] = {}
    results: List[Dict[str, Any]] = []
    for index, raw in enumerate(campaigns):
        campaign = require_mapping(raw, "cross_path_campaigns[%d]" % index)
        kind = require_nonempty_string(campaign.get("kind"), "campaign.kind")
        require(kind in REQUIRED_CAMPAIGN_KINDS, "unknown cross-path campaign kind %s" % kind)
        require(kind not in by_kind, "duplicate cross-path campaign kind %s" % kind)
        by_kind[kind] = campaign
        campaign_id = require_nonempty_string(campaign.get("campaign_id"), "%s.campaign_id" % kind)
        require(campaign.get("common_physics_tag") == COMMON_PHYSICS_TAG,
                "%s has a different common physics tag" % kind)
        require(campaign.get("phase_1_scope") == PHASE_1_SCOPE, "%s has the wrong scope" % kind)
        invariants = require_mapping(campaign.get("invariants"), "%s.invariants" % kind)
        for field in ("field_state_id", "species", "epoch_utc", "frame"):
            require_nonempty_string(invariants.get(field), "%s.invariants.%s" % (kind, field))
        for field in ("boundary_sha256", "integration_sha256", "output_schema_sha256"):
            require_sha256(invariants.get(field), "%s.invariants.%s" % (kind, field))
        expected_sources = (
            ("ANALYTIC_OR_PHENOMENOLOGICAL", "SAMPLED_MESH")
            if kind == "SAMPLED_ANALYTIC_MESH"
            else ("LIVE_SWMF", "OFFLINE_SNAPSHOT_REPLAY")
        )
        require(campaign.get("reference_field_source") == expected_sources[0],
                "%s has wrong reference field source" % kind)
        require(campaign.get("candidate_field_source") == expected_sources[1],
                "%s has wrong candidate field source" % kind)
        pairs = require_sequence(campaign.get("pairs"), "%s.pairs" % kind)
        roles = [require_mapping(pair, "%s pair" % kind).get("role") for pair in pairs]
        require(len(roles) == len(set(roles)), "%s contains duplicate product roles" % kind)
        require(set(roles) == set(REQUIRED_PARITY_ROLES),
                "%s must compare all required roles; missing=%s extra=%s" %
                (kind, sorted(set(REQUIRED_PARITY_ROLES) - set(roles)),
                 sorted(set(roles) - set(REQUIRED_PARITY_ROLES))))
        pair_results = [compare_pair(require_mapping(pair, "%s pair" % kind), kind, resolver) for pair in pairs]
        results.append({"campaign_id": campaign_id, "kind": kind,
                        "status": "PASS" if all(item["status"] == "PASS" for item in pair_results) else "FAIL",
                        "pairs": pair_results})

    missing = set(REQUIRED_CAMPAIGN_KINDS) - set(by_kind)
    require(not missing, "missing cross-path campaigns: %s" % ", ".join(sorted(missing)))
    return {"schema": COMPARISON_SCHEMA,
            "status": "PASS" if all(item["status"] == "PASS" for item in results) else "FAIL",
            "campaigns": results}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="Step-12 release manifest")
    parser.add_argument("--output", required=True, type=Path, help="atomic JSON comparison report")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    manifest_path = args.manifest.expanduser().resolve()
    try:
        manifest = read_json(manifest_path, "release manifest")
        validate_release_manifest_shape(manifest)
        resolver = ResourceResolver(manifest_path)
        report = compare_campaigns(manifest, resolver)
        report["manifest_sha256"] = manifest_sha256(manifest_path)
    except ReleaseValidationError as exc:
        report = {"schema": COMPARISON_SCHEMA, "status": "FAIL", "error": str(exc)}
    atomic_write_json(args.output.expanduser().resolve(), report)
    print("RESULT: %s" % report["status"], flush=True)
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())

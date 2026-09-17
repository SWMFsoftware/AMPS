#!/usr/bin/env python3
"""Run Phase-V linked, cross-model, and observational validation cases.

The runner intentionally owns evidence mechanics, not transport physics.  It
never searches for or invokes the independent ``srcSEP`` application.  A
cross-model producer exports a checksum-identified CSV bundle, and srcSEP3D
consumes that immutable bundle exactly like an observational dataset.  This
keeps the two applications independently buildable while making provenance and
the comparison boundary explicit.

Missing external prerequisites are ``SKIP``.  Malformed evidence is ``ERROR``;
valid evidence outside its declared tolerance is ``FAIL``.  Those categories
must not be collapsed because a release campaign may accept a diagnostic miss
but can never treat absent observations as scientific success.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import statistics
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape as xml_escape


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = Path(__file__).with_name("case_registry.json")


class EvidenceError(RuntimeError):
    """Evidence is present but unusable; this maps to ERROR, not FAIL."""


def _read_json(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise EvidenceError(f"cannot read JSON object {path}: {error}") from error
    if not isinstance(value, dict):
        raise EvidenceError(f"JSON root must be an object: {path}")
    return value


def _write_json_atomic(path: Path, value: Dict[str, Any]) -> None:
    """Publish reports with rename so an interrupted write is never evidence."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise EvidenceError(f"cannot hash {path}: {error}") from error
    return digest.hexdigest()


def _registry() -> List[Dict[str, Any]]:
    payload = _read_json(REGISTRY_PATH)
    if payload.get("schema") != "srcsep3d-validation-case-registry-v1":
        raise EvidenceError(f"unsupported registry schema: {REGISTRY_PATH}")
    cases = payload.get("cases")
    if not isinstance(cases, list) or not cases:
        raise EvidenceError("Phase-V registry contains no cases")
    required = ("id", "name", "evidence_class", "role", "description")
    seen = set()
    result: List[Dict[str, Any]] = []
    for original in cases:
        if not isinstance(original, dict) or any(not original.get(k) for k in required):
            raise EvidenceError("every Phase-V descriptor requires complete metadata")
        case = dict(original)
        case_id = str(case["id"]).upper()
        if case_id in seen:
            raise EvidenceError(f"duplicate Phase-V case ID: {case_id}")
        if case["evidence_class"] not in ("linked", "series", "convergence",
                                          "reserved"):
            raise EvidenceError(f"unsupported evidence class for {case_id}")
        seen.add(case_id)
        case["id"] = case_id
        result.append(case)
    return sorted(result, key=lambda item: str(item["id"]))


def _contained_path(directory: Path, relative: Any, label: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise EvidenceError(f"{label} path is absent")
    path = (directory / relative).resolve()
    root = directory.resolve()
    if path != root and root not in path.parents:
        raise EvidenceError(f"{label} path escapes its evidence directory: {path}")
    if not path.is_file():
        raise EvidenceError(f"{label} file is absent: {path}")
    return path


def _verify_declared_hash(path: Path, declared: Any, label: str) -> str:
    if not isinstance(declared, str) or len(declared) != 64:
        raise EvidenceError(f"{label} requires a complete SHA-256 declaration")
    observed = _sha256(path)
    if observed.lower() != declared.lower():
        raise EvidenceError(
            f"{label} checksum mismatch: declared {declared}, observed {observed}")
    return observed


def _load_series(path: Path) -> List[Tuple[float, float]]:
    """Read the intentionally narrow interchange grammar.

    Unit meanings live in the manifest; the CSV uses canonical column names so
    a time profile cannot accidentally be parsed as a spectrum with guessed
    columns.  Positive values are required because Phase-V profile errors are
    logarithmic and no hidden numerical floor is allowed.
    """
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != ["coordinate_si", "value_si"]:
                raise EvidenceError(
                    f"{path} header must be exactly coordinate_si,value_si")
            result = []
            for line_number, row in enumerate(reader, start=2):
                try:
                    coordinate = float(row["coordinate_si"])
                    value = float(row["value_si"])
                except (TypeError, ValueError) as error:
                    raise EvidenceError(
                        f"non-numeric scientific row {line_number} in {path}") from error
                if not math.isfinite(coordinate) or not math.isfinite(value) or value <= 0.0:
                    raise EvidenceError(
                        f"invalid scientific row {line_number} in {path}")
                if result and coordinate <= result[-1][0]:
                    raise EvidenceError(f"coordinates are not strictly increasing in {path}")
                result.append((coordinate, value))
    except OSError as error:
        raise EvidenceError(f"cannot read scientific CSV {path}: {error}") from error
    if len(result) < 2:
        raise EvidenceError(f"scientific CSV needs at least two rows: {path}")
    return result


def _interpolate_log(series: Sequence[Tuple[float, float]], x: float) -> Optional[float]:
    if x < series[0][0] or x > series[-1][0]:
        return None
    for index, point in enumerate(series):
        if point[0] == x:
            return point[1]
        if point[0] > x:
            left, right = series[index - 1], point
            fraction = (x - left[0]) / (right[0] - left[0])
            return math.exp(math.log(left[1]) +
                            fraction * (math.log(right[1]) - math.log(left[1])))
    return series[-1][1]


def _correlation(left: Sequence[float], right: Sequence[float]) -> float:
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    dl = [item - left_mean for item in left]
    dr = [item - right_mean for item in right]
    left_var = sum(item * item for item in dl)
    right_var = sum(item * item for item in dr)
    if left_var == 0.0 or right_var == 0.0:
        return 1.0 if list(left) == list(right) else 0.0
    return sum(a * b for a, b in zip(dl, dr)) / math.sqrt(left_var * right_var)


def _onset(coordinates: Sequence[float], values: Sequence[float], fraction: float) -> float:
    threshold = max(values) * fraction
    return next(x for x, value in zip(coordinates, values) if value >= threshold)


def _trapezoid(coordinates: Sequence[float], values: Sequence[float]) -> float:
    return sum(0.5 * (values[i - 1] + values[i]) *
               (coordinates[i] - coordinates[i - 1])
               for i in range(1, len(values)))


def _compare_series(model: Sequence[Tuple[float, float]],
                    reference: Sequence[Tuple[float, float]],
                    acceptance: Dict[str, Any]) -> Tuple[Dict[str, float], List[str]]:
    coordinates: List[float] = []
    model_values: List[float] = []
    reference_values: List[float] = []
    for coordinate, value in reference:
        interpolated = _interpolate_log(model, coordinate)
        if interpolated is None:
            continue
        coordinates.append(coordinate)
        model_values.append(interpolated)
        reference_values.append(value)
    if len(model_values) < 2:
        raise EvidenceError("model/reference overlap contains fewer than two points")

    normalization = acceptance.get("normalization", "absolute")
    if normalization == "absolute":
        scale = 1.0
    elif normalization == "one-global-amplitude":
        scale = math.exp(statistics.fmean(
            math.log(ref / value) for value, ref in zip(model_values, reference_values)))
    elif normalization == "unit-peak":
        scale = max(reference_values) / max(model_values)
    else:
        raise EvidenceError(f"unknown comparison normalization: {normalization}")
    model_values = [value * scale for value in model_values]

    model_logs = [math.log10(value) for value in model_values]
    reference_logs = [math.log10(value) for value in reference_values]
    errors = [model_value - reference_value
              for model_value, reference_value in zip(model_logs, reference_logs)]
    model_peak = max(range(len(model_values)), key=model_values.__getitem__)
    reference_peak = max(range(len(reference_values)), key=reference_values.__getitem__)
    model_fluence = _trapezoid(coordinates, model_values)
    reference_fluence = _trapezoid(coordinates, reference_values)
    metrics = {
        "coverage": len(model_values) / len(reference),
        "amplitude_scale": scale,
        "log10_rmse": math.sqrt(statistics.fmean(error * error for error in errors)),
        "median_absolute_log10_error": statistics.median(abs(error) for error in errors),
        "log10_correlation": _correlation(model_logs, reference_logs),
        "onset_coordinate_error": abs(
            _onset(coordinates, model_values, 0.01) -
            _onset(coordinates, reference_values, 0.01)),
        "peak_coordinate_error": abs(
            coordinates[model_peak] - coordinates[reference_peak]),
        "log10_peak_ratio": math.log10(
            model_values[model_peak] / reference_values[reference_peak]),
        "log10_fluence_ratio": math.log10(model_fluence / reference_fluence),
    }
    violations: List[str] = []
    minimums = {"coverage": "minimum_coverage",
                "log10_correlation": "minimum_log10_correlation"}
    maximums = {
        "log10_rmse": "maximum_log10_rmse",
        "median_absolute_log10_error": "maximum_median_absolute_log10_error",
        "onset_coordinate_error": "maximum_onset_coordinate_error",
        "peak_coordinate_error": "maximum_peak_coordinate_error",
    }
    for metric, threshold_name in minimums.items():
        threshold = float(acceptance[threshold_name])
        if metrics[metric] < threshold:
            violations.append(f"{metric}={metrics[metric]:.17g} below {threshold:.17g}")
    for metric, threshold_name in maximums.items():
        threshold = float(acceptance[threshold_name])
        if metrics[metric] > threshold:
            violations.append(f"{metric}={metrics[metric]:.17g} exceeds {threshold:.17g}")
    for metric, threshold_name in (
            ("log10_peak_ratio", "maximum_absolute_log10_peak_ratio"),
            ("log10_fluence_ratio", "maximum_absolute_log10_fluence_ratio")):
        threshold = float(acceptance[threshold_name])
        if abs(metrics[metric]) > threshold:
            violations.append(f"abs({metric})={abs(metrics[metric]):.17g} "
                              f"exceeds {threshold:.17g}")
    return metrics, violations


def _verify_provenance(manifest: Dict[str, Any]) -> None:
    provenance = manifest.get("provenance")
    if not isinstance(provenance, dict):
        raise EvidenceError("scientific evidence has no provenance object")
    for key in ("model", "reference"):
        if not isinstance(provenance.get(key), str) or not provenance[key].strip():
            raise EvidenceError(f"scientific evidence provenance.{key} is empty")


def _run_series(case: Dict[str, Any], evidence_root: Optional[Path]) -> Dict[str, Any]:
    case_id = str(case["id"])
    if evidence_root is None:
        return _result(case, "SKIP", "provide --evidence-root with reviewed evidence")
    case_dir = evidence_root.resolve() / case_id
    manifest_path = case_dir / "manifest.json"
    if not manifest_path.is_file():
        return _result(case, "SKIP", f"evidence manifest is absent: {manifest_path}")
    manifest = _read_json(manifest_path)
    if manifest.get("schema") != "srcsep3d-scientific-evidence-v1" or \
            str(manifest.get("case_id", "")).upper() != case_id:
        raise EvidenceError(f"scientific manifest schema/case mismatch: {manifest_path}")
    for key in ("coordinate_name", "coordinate_units", "value_name", "value_units"):
        if not isinstance(manifest.get(key), str) or not manifest[key].strip():
            raise EvidenceError(f"scientific manifest field {key} is empty")
    _verify_provenance(manifest)
    model_path = _contained_path(case_dir, manifest.get("model_csv"), "model CSV")
    reference_path = _contained_path(
        case_dir, manifest.get("reference_csv"), "reference CSV")
    _verify_declared_hash(model_path, manifest.get("model_sha256"), "model CSV")
    _verify_declared_hash(reference_path, manifest.get("reference_sha256"),
                          "reference CSV")

    exact_pairs = manifest.get("exact_pairs", [])
    if not isinstance(exact_pairs, list):
        raise EvidenceError("exact_pairs must be an array")
    if case.get("requires_exact_pairs") and not exact_pairs:
        raise EvidenceError(f"{case_id} requires at least one exact source pair")
    for index, pair in enumerate(exact_pairs):
        if not isinstance(pair, dict):
            raise EvidenceError(f"exact_pairs[{index}] is not an object")
        left = _contained_path(case_dir, pair.get("model"), "exact model artifact")
        right = _contained_path(case_dir, pair.get("reference"),
                                "exact reference artifact")
        left_hash = _verify_declared_hash(
            left, pair.get("model_sha256"), "exact model artifact")
        right_hash = _verify_declared_hash(
            right, pair.get("reference_sha256"), "exact reference artifact")
        if left_hash != right_hash:
            raise EvidenceError(
                f"exact source pair {index} is not byte-identical for {case_id}")

    metrics, violations = _compare_series(
        _load_series(model_path), _load_series(reference_path),
        dict(case.get("acceptance", {})))
    status = "PASS" if not violations else "FAIL"
    message = ("all declared scientific thresholds satisfied" if not violations
               else "; ".join(violations))
    result = _result(case, status, message)
    result["metrics"] = metrics
    result["artifacts"] = [str(manifest_path), str(model_path), str(reference_path)]
    result["manifest_sha256"] = _sha256(manifest_path)
    return result


def _run_convergence(case: Dict[str, Any], evidence_root: Optional[Path]) -> Dict[str, Any]:
    case_id = str(case["id"])
    if evidence_root is None:
        return _result(case, "SKIP", "provide --evidence-root with convergence evidence")
    manifest_path = evidence_root.resolve() / case_id / "manifest.json"
    if not manifest_path.is_file():
        return _result(case, "SKIP", f"evidence manifest is absent: {manifest_path}")
    manifest = _read_json(manifest_path)
    if manifest.get("schema") != "srcsep3d-convergence-evidence-v1" or \
            str(manifest.get("case_id", "")).upper() != case_id:
        raise EvidenceError(f"convergence manifest schema/case mismatch: {manifest_path}")
    _verify_provenance(manifest)
    levels = manifest.get("levels")
    if not isinstance(levels, list) or len(levels) < 3:
        raise EvidenceError("convergence evidence requires at least three levels")
    parsed: List[Tuple[float, float]] = []
    for item in levels:
        if not isinstance(item, dict):
            raise EvidenceError("convergence level is not an object")
        resolution = float(item.get("resolution", math.nan))
        error = float(item.get("error", math.nan))
        if not math.isfinite(resolution) or resolution <= 0.0 or \
                not math.isfinite(error) or error <= 0.0:
            raise EvidenceError("convergence resolution/error must be positive and finite")
        parsed.append((resolution, error))
    parsed.sort(reverse=True)  # coarse spacing to fine spacing
    orders = [math.log(parsed[i][1] / parsed[i + 1][1]) /
              math.log(parsed[i][0] / parsed[i + 1][0])
              for i in range(len(parsed) - 1)]
    observed_order = min(orders)
    default_difference = float(manifest.get("default_relative_difference", math.nan))
    if not math.isfinite(default_difference) or default_difference < 0.0:
        raise EvidenceError("default_relative_difference is invalid")
    acceptance = dict(case.get("acceptance", {}))
    violations = []
    if observed_order < float(acceptance["minimum_observed_order"]):
        violations.append(f"observed_order={observed_order:.17g} below "
                          f"{float(acceptance['minimum_observed_order']):.17g}")
    if default_difference > float(acceptance["maximum_default_relative_difference"]):
        violations.append(f"default_relative_difference={default_difference:.17g} exceeds "
                          f"{float(acceptance['maximum_default_relative_difference']):.17g}")
    result = _result(case, "PASS" if not violations else "FAIL",
                     "convergence thresholds satisfied" if not violations
                     else "; ".join(violations))
    result["metrics"] = {
        "minimum_observed_order": observed_order,
        "default_relative_difference": default_difference,
    }
    result["artifacts"] = [str(manifest_path)]
    result["manifest_sha256"] = _sha256(manifest_path)
    return result


def _run_linked(case: Dict[str, Any], executable: Optional[Path],
                launch_prefix: str, output_dir: Path,
                timeout: float) -> Dict[str, Any]:
    case_id = str(case["id"])
    if executable is None:
        return _result(case, "SKIP", "provide --amps with a configured linked executable")
    executable = executable.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise EvidenceError(f"linked executable is missing or not executable: {executable}")
    prefix = shlex.split(launch_prefix)

    def execute(arguments: Sequence[str]) -> subprocess.CompletedProcess[str]:
        try:
            return subprocess.run([*prefix, str(executable), *arguments],
                                  cwd=str(ROOT), text=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, check=False,
                                  timeout=timeout)
        except (OSError, subprocess.TimeoutExpired) as error:
            raise EvidenceError(f"linked command failed to execute: {error}") from error

    listed = execute(["--list-tests"])
    if listed.returncode != 0:
        raise EvidenceError(
            f"linked executable rejected --list-tests: {listed.stdout[-2000:]}")
    advertised = {line.split("|", 1)[0].strip().upper()
                  for line in listed.stdout.splitlines() if "|" in line}
    if case_id not in advertised:
        raise EvidenceError(
            f"linked executable does not advertise {case_id}; rebuild it from this Phase-V tree")
    native_report = output_dir / case_id / "native.json"
    native_report.parent.mkdir(parents=True, exist_ok=True)
    completed = execute(["--test", case_id, "--test-json", str(native_report),
                         "--artifact-directory", str(native_report.parent)])
    if not native_report.is_file():
        raise EvidenceError(
            f"linked case wrote no JSON report (exit {completed.returncode}): "
            f"{completed.stdout[-2000:]}")
    payload = _read_json(native_report)
    records = payload.get("results")
    if not isinstance(records, list):
        raise EvidenceError("linked JSON report contains no results array")
    record = next((item for item in records if isinstance(item, dict) and
                   str(item.get("id", "")).upper() == case_id), None)
    if record is None:
        raise EvidenceError(f"linked JSON report does not contain {case_id}")
    status = str(record.get("status", "ERROR")).upper()
    if status not in ("PASS", "FAIL", "SKIP", "ERROR"):
        status = "ERROR"
    result = _result(case, status, str(record.get("message", "")))
    result["metrics"] = record.get("metrics", [])
    result["artifacts"] = [str(native_report), *record.get("artifacts", [])]
    result["command"] = [*prefix, str(executable), "--test", case_id]
    return result


def _result(case: Dict[str, Any], status: str, message: str) -> Dict[str, Any]:
    return {
        "id": str(case["id"]),
        "name": str(case["name"]),
        "evidence_class": str(case["evidence_class"]),
        "role": str(case["role"]),
        "status": status,
        "message": message,
        "metrics": {},
        "artifacts": [],
        "elapsed_seconds": 0.0,
    }


def _totals(results: Iterable[Dict[str, Any]]) -> Dict[str, int]:
    totals = {"passed": 0, "failed": 0, "skipped": 0, "errors": 0}
    mapping = {"PASS": "passed", "FAIL": "failed", "SKIP": "skipped",
               "ERROR": "errors"}
    for result in results:
        totals[mapping.get(str(result.get("status", "ERROR")), "errors")] += 1
    return totals


def _write_junit(path: Path, results: Sequence[Dict[str, Any]],
                 totals: Dict[str, int]) -> None:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (f'<testsuite name="srcSEP3D Phase V" tests="{len(results)}" '
         f'failures="{totals["failed"]}" errors="{totals["errors"]}" '
         f'skipped="{totals["skipped"]}">'),
    ]
    for result in results:
        case_id = xml_escape(str(result["id"]))
        message = xml_escape(str(result["message"]))
        lines.append(f'  <testcase classname="srcSEP3D.validation" name="{case_id}" '
                     f'time="{float(result["elapsed_seconds"]):.9f}">')
        if result["status"] == "FAIL":
            lines.append(f'    <failure message="{message}"/>')
        elif result["status"] == "ERROR":
            lines.append(f'    <error message="{message}"/>')
        elif result["status"] == "SKIP":
            lines.append(f'    <skipped message="{message}"/>')
        lines.append("  </testcase>")
    lines.extend(("</testsuite>", ""))
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text("\n".join(lines), encoding="utf-8")
    os.replace(temporary, path)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run srcSEP3D Phase-V external validation evidence")
    selectors = parser.add_mutually_exclusive_group(required=True)
    selectors.add_argument("--list", action="store_true")
    selectors.add_argument("--case", action="append", default=[])
    selectors.add_argument("--all", action="store_true")
    parser.add_argument("--amps", type=Path,
                        help="configured linked srcSEP3D/AMPS executable")
    parser.add_argument("--launch-prefix", default="",
                        help="optional direct argv prefix, e.g. 'mpiexec -n 8'")
    parser.add_argument("--evidence-root", type=Path,
                        help="directory containing CASE_ID/manifest.json bundles")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "test_output" / "phase-v")
    parser.add_argument("--timeout", type=float, default=3600.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.timeout <= 0.0:
        raise EvidenceError("--timeout must be positive")
    cases = _registry()
    if args.list:
        print(f"{'ID':<10} {'CLASS':<12} {'ROLE':<13} NAME")
        for case in cases:
            print(f"{case['id']:<10} {case['evidence_class']:<12} "
                  f"{case['role']:<13} {case['name']}")
        return 0
    by_id = {str(case["id"]): case for case in cases}
    requested = sorted(by_id) if args.all else sorted({item.upper() for item in args.case})
    unknown = sorted(set(requested) - set(by_id))
    if unknown:
        raise EvidenceError("unknown Phase-V case(s): " + ", ".join(unknown))
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_root = (args.evidence_root.expanduser().resolve()
                     if args.evidence_root is not None else None)
    results: List[Dict[str, Any]] = []
    for case_id in requested:
        case = by_id[case_id]
        started = time.monotonic()
        try:
            kind = case["evidence_class"]
            if kind == "linked":
                result = _run_linked(case, args.amps, args.launch_prefix,
                                     output_dir, args.timeout)
            elif kind == "series":
                result = _run_series(case, evidence_root)
            elif kind == "convergence":
                result = _run_convergence(case, evidence_root)
            else:
                result = _result(case, "SKIP", str(case["skip_reason"]))
        except EvidenceError as error:
            result = _result(case, "ERROR", str(error))
        result["elapsed_seconds"] = time.monotonic() - started
        _write_json_atomic(output_dir / case_id / "result.json", result)
        results.append(result)
        print(f"[{case_id}] {result['status']} {result['message']}", flush=True)

    totals = _totals(results)
    summary = {
        "schema": "srcsep3d-validation-summary-v1",
        "totals": totals,
        "results": results,
    }
    _write_json_atomic(output_dir / "validation-summary.json", summary)
    _write_junit(output_dir / "validation-summary.xml", results, totals)
    return 2 if totals["errors"] else (1 if totals["failed"] else 0)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except EvidenceError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)

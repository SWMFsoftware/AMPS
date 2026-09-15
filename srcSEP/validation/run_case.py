#!/usr/bin/env python3
"""Execute registered end-to-end cases through linked srcSEP/AMPS.

This runner owns only common campaign mechanics: deterministic case selection,
configuration-path resolution, linked-executable verification, isolated output
directories, report aggregation, and exit-status propagation. Each registered
case owns its native invocation, independent reference, metrics, and scientific
acceptance rule. That boundary lets CV02 and later cases use a different
reference method without creating a second command-line convention or weakening
the shared evidence contract.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Dict, Iterable, List, Optional, Sequence
from xml.sax.saxutils import escape as xml_escape


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "validation" / "case_registry.json"

# Cross-model publication cases are immutable scientific definitions.  The
# public test/run_tests.py front end enforces the same rule, while duplicating
# the guard here protects operators who invoke validation/run_case.py directly.
FIXED_PUBLICATION_INPUT_CASES = {
    "XM02", "XM03", "OV01", "OV02", "OV03", "OV04", "OV05", "EV01", "EV02"
}


class CaseRunnerError(RuntimeError):
    """Describe a configuration or execution failure visible to operators."""


def _load_json(path: Path) -> Dict[str, Any]:
    """Read one JSON object and retain the filename in any diagnostic."""
    try:
        with path.open("r", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise CaseRunnerError(f"cannot read JSON object {path}: {error}") from error
    if not isinstance(value, dict):
        raise CaseRunnerError(f"JSON root must be an object: {path}")
    return value


def _write_json(path: Path, value: Dict[str, Any]) -> None:
    """Commit JSON through a same-directory temporary file.

    Validation reports are release evidence.  A partially written report must
    never be mistaken for a completed case after a process or filesystem error.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_junit(path: Path, results: Sequence[Dict[str, Any]],
                 totals: Dict[str, int]) -> None:
    """Write a compact JUnit mirror of the authoritative JSON result.

    XML consumers receive the same categorical status and diagnostic, while
    the full configuration/metrics/artifact record remains in JSON.  Building
    the document in memory permits the same atomic commit used for JSON.
    """
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (f'<testsuite name="srcSEP validation cases" tests="{len(results)}" '
         f'failures="{totals["failed"]}" errors="{totals["errors"]}" '
         f'skipped="{totals["skipped"]}">'),
    ]
    for result in results:
        identifier = xml_escape(str(result.get("id", "unknown")))
        elapsed = float(result.get("elapsed_seconds", 0.0))
        message = xml_escape(str(result.get("message", "")))
        status = str(result.get("status", "ERROR")).upper()
        lines.append(f'  <testcase classname="srcSEP.validation" name="{identifier}" '
                     f'time="{elapsed:.17g}">')
        if status == "FAIL":
            lines.append(f'    <failure message="{message}"/>')
        elif status == "ERROR":
            lines.append(f'    <error message="{message}"/>')
        elif status == "SKIP":
            lines.append(f'    <skipped message="{message}"/>')
        lines.append(f'    <system-out>{message}</system-out>')
        lines.append("  </testcase>")
    lines.append("</testsuite>")
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        stream.write("\n".join(lines) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _registry() -> List[Dict[str, Any]]:
    payload = _load_json(REGISTRY_PATH)
    if payload.get("schema") != "srcsep-validation-case-registry-v1":
        raise CaseRunnerError(f"unsupported case registry schema in {REGISTRY_PATH}")
    cases = payload.get("cases")
    if not isinstance(cases, list) or not cases:
        raise CaseRunnerError("validation case registry contains no cases")

    seen = set()
    validated: List[Dict[str, Any]] = []
    required = ("id", "name", "group", "runtime_class", "entrypoint",
                "default_input", "description")
    for case in cases:
        if not isinstance(case, dict) or any(not case.get(key) for key in required):
            raise CaseRunnerError("each validation descriptor requires complete metadata")
        canonical = str(case["id"]).upper()
        if canonical in seen:
            raise CaseRunnerError(f"duplicate validation case ID: {canonical}")
        seen.add(canonical)
        descriptor = dict(case)
        descriptor["id"] = canonical
        for key in ("entrypoint", "default_input"):
            candidate = (ROOT / str(descriptor[key])).resolve()
            if ROOT not in candidate.parents or not candidate.is_file():
                raise CaseRunnerError(
                    f"{canonical} {key} is missing or outside the source tree: {candidate}")
        validated.append(descriptor)
    return sorted(validated, key=lambda item: str(item["id"]))


def _load_entrypoint(descriptor: Dict[str, Any]):
    """Load a case module by registered path, not by user-controlled import name."""
    path = (ROOT / str(descriptor["entrypoint"])).resolve()
    module_name = f"srcsep_validation_{descriptor['id'].casefold()}"
    specification = importlib.util.spec_from_file_location(module_name, path)
    if specification is None or specification.loader is None:
        raise CaseRunnerError(f"cannot load validation entrypoint {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    if not callable(getattr(module, "run_case", None)):
        raise CaseRunnerError(f"validation entrypoint has no run_case(): {path}")
    return module


def _totals(results: Iterable[Dict[str, Any]]) -> Dict[str, int]:
    totals = {"passed": 0, "failed": 0, "skipped": 0, "errors": 0}
    keys = {"PASS": "passed", "FAIL": "failed", "SKIP": "skipped",
            "ERROR": "errors"}
    for result in results:
        key = keys.get(str(result.get("status", "ERROR")).upper(), "errors")
        totals[key] += 1
    return totals


def _verify_linked_registry(executable: Path,
                            selected_ids: Sequence[str],
                            timeout: Optional[float]) -> None:
    """Fail before case setup unless the requested linked binary advertises it.

    Executing ``--list-tests`` is deliberately part of the gate: an arbitrary
    executable at the requested path, or a stale srcSEP build without CV01,
    must not be accepted merely because it is executable.
    """
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise CaseRunnerError(
            f"linked srcSEP/AMPS executable is missing or not executable: {executable}")
    try:
        completed = subprocess.run(
            [str(executable), "--list-tests"], cwd=str(ROOT), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            timeout=timeout)
    except subprocess.TimeoutExpired as error:
        raise CaseRunnerError(
            "linked srcSEP/AMPS --list-tests exceeded the configured timeout") from error
    except OSError as error:
        raise CaseRunnerError(
            f"cannot execute linked srcSEP/AMPS application {executable}: {error}") from error
    if completed.returncode != 0:
        raise CaseRunnerError(
            "linked srcSEP/AMPS executable rejected --list-tests "
            f"with status {completed.returncode}: {completed.stdout.strip()}")
    advertised = {
        line.split("|", 1)[0].strip().upper()
        for line in completed.stdout.splitlines() if "|" in line
    }
    missing = sorted(set(selected_ids) - advertised)
    if missing:
        raise CaseRunnerError(
            "linked srcSEP/AMPS executable does not advertise validation "
            "case(s): " + ", ".join(missing) +
            "; rebuild the application from this source tree")


def run(selected_ids: Sequence[str], output_dir: Path,
        input_override: Optional[Path], executable: Path,
        timeout: Optional[float]) -> int:
    """Run selected descriptors in registry order and write one common report."""
    descriptors = _registry()
    by_id = {str(item["id"]): item for item in descriptors}
    requested = [item.upper() for item in selected_ids]
    unknown = sorted(set(requested) - set(by_id))
    if unknown:
        raise CaseRunnerError("unknown validation case ID(s): " + ", ".join(unknown))
    if input_override is not None and len(requested) != 1:
        raise CaseRunnerError("--input may be used only when exactly one case is selected")
    if input_override is not None and requested[0] in FIXED_PUBLICATION_INPUT_CASES:
        raise CaseRunnerError(
            f"{requested[0]} uses its single registered publication-derived input; "
            "remove --input")

    # The validation portfolio is explicitly an application-level gate. Every
    # case must be present in the linked registry before Python generates any
    # evidence or calls a reference solution.
    _verify_linked_registry(executable, requested, timeout)

    output_dir.mkdir(parents=True, exist_ok=True)
    results: List[Dict[str, Any]] = []
    case_records: List[Dict[str, Any]] = []
    for case_id in sorted(set(requested)):
        descriptor = by_id[case_id]
        case_output = output_dir / case_id
        case_output.mkdir(parents=True, exist_ok=True)
        input_path = (input_override if input_override is not None else
                      ROOT / str(descriptor["default_input"])).resolve()
        if not input_path.is_file():
            raise CaseRunnerError(f"input for {case_id} does not exist: {input_path}")

        module = _load_entrypoint(descriptor)
        try:
            result = module.run_case(
                source_root=ROOT,
                input_path=input_path,
                output_dir=case_output,
                executable=executable,
                timeout=timeout,
            )
        except Exception as error:  # Preserve later cases after one local error.
            result = {
                "id": case_id,
                "status": "ERROR",
                "message": f"case entrypoint raised {type(error).__name__}: {error}",
                "elapsed_seconds": 0.0,
                "seed": None,
                "configuration": [f"input={input_path}"],
                "metrics": [],
                "artifacts": [],
            }
        if not isinstance(result, dict):
            raise CaseRunnerError(f"{case_id} run_case() did not return a result object")
        result["id"] = case_id
        results.append(result)
        case_records.append({
            "id": case_id,
            "entrypoint": str(descriptor["entrypoint"]),
            "input": str(input_path),
            "input_sha256": _sha256(input_path),
            "output_directory": str(case_output),
        })
        print(f"{str(result.get('status', 'ERROR')).upper()} {case_id}: "
              f"{result.get('message', 'no diagnostic')}")

    totals = _totals(results)
    exit_code = 2 if totals["errors"] else (1 if totals["failed"] else 0)
    report = {
        "schema": "srcsep-component-tests-v1",
        "exit_code": exit_code,
        "totals": totals,
        "results": results,
    }
    _write_json(output_dir / "srcsep-tests.json", report)
    _write_junit(output_dir / "srcsep-tests.xml", results, totals)
    _write_json(output_dir / "validation-run-manifest.json", {
        "schema": "srcsep-validation-case-run-v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "registry": str(REGISTRY_PATH),
        "registry_sha256": _sha256(REGISTRY_PATH),
        "source_root": str(ROOT),
        "python": sys.version,
        "platform": platform.platform(),
        "linked_executable": str(executable),
        "linked_executable_sha256": _sha256(executable),
        "cases": case_records,
        "report": "srcsep-tests.json",
        "report_exit_code": exit_code,
    })
    print(f"Validation summary: PASS={totals['passed']} FAIL={totals['failed']} "
          f"SKIP={totals['skipped']} ERROR={totals['errors']}")
    return exit_code


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run registered end-to-end srcSEP validation cases.")
    parser.add_argument("--list", action="store_true",
                        help="list case IDs without compiling or running a model")
    parser.add_argument("--case", action="append", default=[], metavar="ID",
                        help="case ID to execute; repeatable")
    parser.add_argument("--all", action="store_true",
                        help="execute every case in deterministic registry order")
    parser.add_argument("--input", type=Path,
                        help=("override one selected CV/IV input; XM02/XM03 and "
                              "OV01-OV05 always use their registered publication-derived input"))
    parser.add_argument("--output-dir", type=Path, required=False,
                        default=Path("test_output") / "validation-cases",
                        help="directory for model, reference, metrics, and manifests")
    parser.add_argument(
        "--amps", default=os.environ.get("SEP_EXECUTABLE", str(ROOT.parent / "amps")),
        help="linked srcSEP/AMPS executable that owns the numerical model run")
    parser.add_argument("--timeout", type=float,
                        help="timeout in seconds for each linked/reference command")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    descriptors = _registry()
    if args.list:
        for item in descriptors:
            print(f"{item['id']} | {item['group']} | {item['runtime_class']} | "
                  f"{item['name']}")
        return 0
    selected = [str(item["id"]) for item in descriptors] if args.all else args.case
    if not selected:
        raise CaseRunnerError("select --case ID, --all, or --list")
    if args.timeout is not None and args.timeout <= 0.0:
        raise CaseRunnerError("--timeout must be positive")
    return run(selected, args.output_dir.expanduser().resolve(),
               args.input.expanduser().resolve() if args.input else None,
               Path(args.amps).expanduser().resolve(), args.timeout)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CaseRunnerError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)

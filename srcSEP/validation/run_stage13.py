#!/usr/bin/env python3
"""Run cases from the spherical Stage 13 AMPS .in deck and save evidence."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
from run_case import CaseRunnerError, _load_json, _sha256, _totals, _write_json, _write_junit, run


ROOT = Path(__file__).resolve().parents[1]
DECK = ROOT.parent / "input/sep_stage13.in"
REFERENCE_SCHEMA = "sep-stage13-observational-case-v1"


def native(executable, deck, arguments, work, log, timeout):
    command = [str(executable), "-i", str(deck), *arguments]
    process = subprocess.run(command, cwd=work, text=True, errors="replace",
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             timeout=timeout, check=False)
    log.write_text(f"COMMAND: {command!r}\n{process.stdout}\nEXIT: {process.returncode}\n",
                   encoding="utf-8")
    return process


def selection(lines):
    """Read the resolved values printed by the generic C++ case selector."""
    fields = {}
    for line in lines.splitlines():
        if line.startswith("stage13."):
            key, value = line[len("stage13."):].split("=", 1)
            if key in fields:
                raise CaseRunnerError(f"duplicate native result field {key}")
            fields[key] = value
    return fields


def verify_reference(path, name, continuity):
    """A present artifact with wrong bytes is an error, never a missing-data skip."""
    package = _load_json(path)
    if package.get("schema") != REFERENCE_SCHEMA or package.get("case_id") != name:
        raise CaseRunnerError(f"observation manifest schema or case mismatch in {path}")
    if package.get("existing_case") != continuity:
        raise CaseRunnerError(f"continuity reference disagrees with the selected deck in {path}")
    comparison = package.get("comparison", {})
    if (not isinstance(comparison, dict) or comparison.get("status") != "NOT_RUN" or
            comparison.get("validation_claim") is not False):
        raise CaseRunnerError(f"an observation package cannot supply its own comparison verdict: {path}")
    missing = package.get("missing_required_inputs")
    if not isinstance(missing, list) or any(not isinstance(item, str) or not item for item in missing):
        raise CaseRunnerError(f"invalid missing_required_inputs in {path}")
    artifacts = package.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise CaseRunnerError(f"invalid artifacts list in {path}")
    verified = []
    seen = set()
    for artifact in artifacts:
        recorded_path = Path(artifact["path"])
        file = (path.parent / recorded_path).resolve()
        if recorded_path.is_absolute() or not file.is_file() or file in seen:
            raise CaseRunnerError(f"unresolved or duplicate observation artifact {file}")
        seen.add(file)
        if file.stat().st_size != artifact["bytes"] or _sha256(file) != artifact["sha256"]:
            raise CaseRunnerError(f"observation artifact differs from its recorded bytes: {file}")
        verified.append({**artifact, "resolved_path": str(file)})
    return package, verified


def run_selected(executable, deck, name, work, timeout, input_only):
    start = time.monotonic()
    result = {"id": name, "status": "ERROR", "metrics": [], "artifacts": []}
    try:
        arguments = ["-case", name]
        if input_only:
            arguments.append("--validate-runtime-input")
        log = work / "native.log"
        result["artifacts"].append(str(log))
        process = native(executable, deck, arguments, work, log, timeout)
        fields = selection(process.stdout)
        result["native_exit_code"] = process.returncode
        result["configuration"] = [f"{key}={value}" for key, value in fields.items()]
        expected = 0 if input_only else 77
        expected_terminal = ("runtime_input", "PASS") if input_only else ("comparison", "SKIP")
        if (fields.get("case") != name or process.returncode != expected or
                fields.get(expected_terminal[0]) != expected_terminal[1]):
            raise CaseRunnerError(f"native validation returned unexpected status; see {log}")
        manifest = Path(fields["manifest"])
        continuity = fields.get("continuity_reference")
        package, verified = verify_reference(manifest, name, continuity)
        result["observation_package"] = {
            "path": str(manifest), "sha256": _sha256(manifest),
            "readiness": package["readiness"], "verified_artifacts": verified,
            "missing_required_inputs": package["missing_required_inputs"],
        }
        # Keep the existing OV01 run and its scoring contract authoritative.
        # Its linked wind/source/boundaries are different, so its report cannot
        # become a spherical-model validation PASS even if OV01 itself passes.
        if continuity and not input_only:
            baseline_dir = work / "continuity"
            status = run([continuity], baseline_dir, None, executable, timeout)
            report = baseline_dir / "srcsep-tests.json"
            baseline = _load_json(report)["results"]
            if len(baseline) != 1 or baseline[0]["id"] != continuity:
                raise CaseRunnerError(f"unexpected linked continuity result in {report}")
            result["continuity"] = {
                "id": continuity, "status": baseline[0]["status"],
                "exit_code": status, "report": str(report),
            }
            result["artifacts"].append(str(report))
            if status != 0 or baseline[0]["status"] != "PASS":
                result["status"] = baseline[0]["status"]
                result["message"] = "The maintained OV01 continuity gate did not pass."
                return result
        result["status"] = "PASS" if input_only else "SKIP"
        result["message"] = ("Resolved runtime selection and verified observation package integrity."
                             if input_only else fields["reason"])
    except (CaseRunnerError, OSError, subprocess.SubprocessError, KeyError, TypeError, ValueError) as error:
        result["message"] = str(error)
    finally:
        result["spherical_observation_validation"] = "NOT_RUN"
        result["elapsed_seconds"] = time.monotonic() - start
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, default=ROOT.parent / "amps")
    parser.add_argument("-i", "--input", type=Path, default=DECK)
    selected = parser.add_mutually_exclusive_group(required=True)
    selected.add_argument("-case", "--case", action="append", metavar="NAME")
    selected.add_argument("--all", action="store_true")
    selected.add_argument("--list", action="store_true")
    parser.add_argument("--validate-runtime-input", action="store_true",
                        help="verify the selection and reference hashes without scoring")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be finite and positive")
    executable, deck = args.amps.resolve(), args.input.resolve()
    if args.list:
        return subprocess.run([str(executable), "-i", str(deck), "--list-validation-cases"],
                              timeout=args.timeout, check=False).returncode
    if args.output_dir is None:
        parser.error("--output-dir is required for a case or campaign run")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    discovery = native(executable, deck, ["--list-validation-cases"], output,
                       output / "discovery.log", args.timeout)
    available = discovery.stdout.splitlines()
    if discovery.returncode != 0 or not available or len(set(available)) != len(available):
        raise CaseRunnerError(f"cannot enumerate linked runtime cases; see {output / 'discovery.log'}")
    requested = available if args.all else args.case
    unknown = set(requested) - set(available)
    if unknown:
        raise CaseRunnerError(f"unknown cases {', '.join(sorted(unknown))}; available: {', '.join(available)}")
    results = []
    for name in available:
        if name not in requested:
            continue
        work = output / name
        work.mkdir()
        result = run_selected(executable, deck, name, work, args.timeout, args.validate_runtime_input)
        results.append(result)
        print(f"{result['status']} {name} {result['message']}", flush=True)
    totals = _totals(results)
    code = 2 if totals["errors"] else 1 if totals["failed"] else 77 if totals["skipped"] else 0
    _write_json(output / "stage13-validation.json", {
        "schema": "srcsep-stage13-validation-run-v1", "exit_code": code,
        "mode": "input-integrity" if args.validate_runtime_input else "observation",
        "executable": str(executable), "executable_sha256": _sha256(executable),
        "deck": str(deck), "deck_sha256": _sha256(deck), "totals": totals, "results": results,
    })
    _write_junit(output / "stage13-validation.xml", results, totals)
    return code


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CaseRunnerError, OSError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)

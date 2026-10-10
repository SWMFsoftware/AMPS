#!/usr/bin/env python3
"""Run every registered srcMoon U-series test and summarize the campaign.

Each test is launched through its own ``UXX_*/test.py`` command-line runner.
The aggregate runner preserves that process's combined stdout/stderr and reads
the structured ``UXX/result.json`` artifact to determine scientific status.
It never infers PASS from a zero process exit alone.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any


TEST_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TEST_ROOT.parents[1]
DEFAULT_BUILD = REPO_ROOT / "build"
DEFAULT_OUTPUT_PARENT = REPO_ROOT / "test_output" / "srcMoon"

# This is the campaign registry, not a directory-discovery convenience.  A new
# test must be added here deliberately so code review can see that the aggregate
# campaign changed.  validate_registry() also rejects an unregistered UXX
# directory, preventing a newly added test from being silently omitted.
TESTS: tuple[tuple[str, str], ...] = (
    ("U01", "U01_los_geometry/test.py"),
    ("U02", "U02_target_isolation/test.py"),
    ("U03", "U03_lunar_gravity/test.py"),
    ("U04", "U04_rotating_frame/test.py"),
    ("U05", "U05_na_radiation_shadow/test.py"),
    ("U06", "U06_lorentz/test.py"),
    ("U07", "U07_surface_refinement/test.py"),
    ("U08", "U08_lola_geometry/test.py"),
    ("U09", "U09_terrain_shadow/test.py"),
    ("U10", "U10_diviner_temperature/test.py"),
    ("U11", "U11_thermal_inertia/test.py"),
    ("U12", "U12_surface_interaction/test.py"),
    ("U13", "U13_cold_trapping/test.py"),
    ("U14", "U14_na_photoionization/test.py"),
    ("U15", "U15_multispecies_photoionization/test.py"),
    ("U16", "U16_electron_impact/test.py"),
    ("U17", "U17_charge_exchange/test.py"),
    ("U18", "U18_plasma_driver/test.py"),
    ("U19", "U19_helium_source/test.py"),
    ("U20", "U20_neon_source/test.py"),
    ("U21", "U21_argon_source/test.py"),
    ("U22", "U22_na_sources/test.py"),
    ("U23", "U23_meteoroid_driver/test.py"),
    ("U24", "U24_water_chemistry/test.py"),
    ("U25", "U25_data_provenance/test.py"),
    ("U26", "U26_observation_geometry/test.py"),
)

VALID_STATUSES = ("PASS", "FAIL", "SKIPPED", "ERROR")
EXPECTED_EXIT = {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}


def validate_registry() -> list[str]:
    """Return registry errors without running any test or modifying outputs."""
    errors: list[str] = []
    ids = [test_id for test_id, _ in TESTS]
    if len(ids) != len(set(ids)):
        errors.append("the run_all_test.py registry contains duplicate test IDs")

    registered_directories: set[str] = set()
    for test_id, runner_name in TESTS:
        runner = TEST_ROOT / runner_name
        registered_directories.add(runner.parts[-2])
        if not runner.is_file():
            errors.append(f"{test_id}: registered runner does not exist: {runner}")
        if runner.parts[-2][:3] != test_id:
            errors.append(
                f"{test_id}: runner directory does not match its ID: {runner_name}"
            )

    discovered_directories = {
        path.name
        for path in TEST_ROOT.glob("U[0-9][0-9]_*")
        if path.is_dir()
    }
    for directory in sorted(discovered_directories - registered_directories):
        errors.append(
            f"unregistered test directory {directory}; add its runner to TESTS"
        )
    for directory in sorted(registered_directories - discovered_directories):
        errors.append(f"registered test directory is missing: {directory}")
    return errors


def git_revision() -> str:
    """Return the source revision recorded in the aggregate result."""
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def default_output_directory() -> Path:
    """Create a run-specific default name so prior campaign evidence survives."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return DEFAULT_OUTPUT_PARENT / f"run-all-{timestamp}"


def run_logged(command: list[str], log_path: Path, stream: bool) -> int:
    """Run one command with a durable combined stdout/stderr transcript."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("command: " + " ".join(command) + "\n")
        log.flush()
        try:
            if not stream:
                completed = subprocess.run(
                    command,
                    cwd=REPO_ROOT,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
                return completed.returncode

            # Streaming is optional because most tests are short.  When
            # enabled, each line is copied to the terminal and durable log.
            process = subprocess.Popen(
                command,
                cwd=REPO_ROOT,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)
            return process.wait()
        except OSError as error:
            # Treat an unlaunchable compiler/interpreter/runner as
            # infrastructure ERROR and retain the operating-system diagnostic.
            message = f"ERROR: could not execute command: {error}\n"
            log.write(message)
            log.flush()
            if stream:
                print(message, end="", file=sys.stderr, flush=True)
            return 127


def clean_build(output_root: Path, stream: bool) -> dict[str, Any]:
    """Perform the repository-required clean ``make -j`` build."""
    if DEFAULT_BUILD.exists():
        # AMPS generates configured source beneath build/.  A clean tree is
        # required so the test campaign cannot link against stale copies.
        shutil.rmtree(DEFAULT_BUILD)
    log_path = output_root / "build.log"
    started = time.monotonic()
    returncode = run_logged(["make", "-j"], log_path, stream)
    return {
        "status": "PASS" if returncode == 0 else "ERROR",
        "command": ["make", "-j"],
        "returncode": returncode,
        "log": str(log_path),
        "duration_seconds": time.monotonic() - started,
    }


def load_runner_result(
    test_id: str, runner: Path, result_path: Path, log_path: Path, returncode: int
) -> dict[str, Any]:
    """Validate one runner's JSON/exit-code protocol and normalize its record."""
    base = {
        "id": test_id,
        "runner": str(runner),
        "log": str(log_path),
        "result_json": str(result_path),
        "returncode": returncode,
    }
    try:
        with result_path.open(encoding="utf-8") as stream:
            payload = json.load(stream)
    except FileNotFoundError:
        return {
            **base,
            "status": "ERROR",
            "reason": "individual runner did not create its required result.json",
        }
    except (OSError, json.JSONDecodeError) as error:
        return {
            **base,
            "status": "ERROR",
            "reason": f"individual runner produced unreadable result.json: {error}",
        }

    status = payload.get("status")
    if payload.get("id") != test_id or status not in VALID_STATUSES:
        return {
            **base,
            "status": "ERROR",
            "reason": (
                "result.json must contain the registered id and one of "
                + ", ".join(VALID_STATUSES)
            ),
            "reported_result": payload,
        }

    expected_exit = EXPECTED_EXIT[status]
    if returncode != expected_exit:
        return {
            **base,
            "status": "ERROR",
            "reason": (
                f"runner exit {returncode} disagrees with {status} "
                f"result (expected {expected_exit})"
            ),
            "reported_status": status,
            "name": payload.get("name", test_id),
        }

    return {
        **base,
        "status": status,
        "name": payload.get("name", test_id),
        "reason": payload.get("reason"),
    }


def run_test(
    test_id: str,
    runner_name: str,
    build_dir: Path,
    output_root: Path,
    stream: bool,
) -> dict[str, Any]:
    """Execute one registered CLI runner and consume its structured result."""
    runner = (TEST_ROOT / runner_name).resolve()
    result_path = output_root / test_id / "result.json"
    if result_path.exists():
        # A stale artifact must never let a crashed runner appear successful.
        result_path.unlink()
    log_path = output_root / "runner_logs" / f"{test_id}.log"
    command = [
        sys.executable,
        str(runner),
        "--build-dir",
        str(build_dir),
        "--output-dir",
        str(output_root),
    ]
    started = time.monotonic()
    returncode = run_logged(command, log_path, stream)
    result = load_runner_result(test_id, runner, result_path, log_path, returncode)
    result["duration_seconds"] = time.monotonic() - started
    return result


def print_summary(counts: dict[str, int], results: list[dict[str, Any]]) -> None:
    """Print totals and mandatory diagnostics for FAIL/ERROR outcomes."""
    print("\nsrcMoon test summary")
    print(f"  executed: {counts['executed']}")
    print(f"  passed:   {counts['PASS']}")
    print(f"  failed:   {counts['FAIL']}")
    print(f"  skipped:  {counts['SKIPPED']}")
    print(f"  errors:   {counts['ERROR']}")

    problems = [result for result in results if result["status"] in {"FAIL", "ERROR"}]
    if problems:
        print("\nFailed/error tests:")
        for result in problems:
            print(f"  {result['status']} {result['id']} {result.get('name', '')}".rstrip())
            print(f"    runner: {result['runner']}")
            print(f"    log:    {result['log']}")
            if result.get("reason"):
                print(f"    reason: {result['reason']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list the explicit registry")
    parser.add_argument(
        "--only", action="append", metavar="UXX", help="run only a registered ID"
    )
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="new/empty artifact directory (default: UTC-stamped under test_output)",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="reuse an existing generated build instead of clean make -j",
    )
    parser.add_argument("--stream", action="store_true", help="also stream runner logs")
    parser.add_argument(
        "--strict-skips",
        action="store_true",
        help="return failure when any test is SKIPPED",
    )
    args = parser.parse_args()

    registry_errors = validate_registry()
    if registry_errors:
        for error in registry_errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2

    registry = dict(TESTS)
    if args.list:
        for test_id, runner in TESTS:
            print(f"{test_id}  {runner}")
        return 0

    selected = list(registry) if args.only is None else args.only
    unknown = [test_id for test_id in selected if test_id not in registry]
    if unknown:
        parser.error("unknown test ID(s): " + ", ".join(unknown))
    if len(selected) != len(set(selected)):
        parser.error("a test ID may be selected only once")

    output_root = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else default_output_directory().resolve()
    )
    build_dir = args.build_dir.resolve()
    if output_root.exists() and any(output_root.iterdir()):
        parser.error(
            f"output directory is not empty: {output_root}; choose a new directory "
            "so prior test evidence is not overwritten"
        )
    output_root.mkdir(parents=True, exist_ok=True)
    campaign_started = time.monotonic()

    build_result: dict[str, Any] | None = None
    if not args.skip_build:
        if build_dir != DEFAULT_BUILD.resolve():
            parser.error("a clean build can target only the repository build directory")
        build_result = clean_build(output_root, args.stream)
        if build_result["status"] == "ERROR":
            summary = {
                "source_revision": git_revision(),
                "build": build_result,
                "counts": {
                    "executed": 0,
                    "PASS": 0,
                    "FAIL": 0,
                    "SKIPPED": 0,
                    "ERROR": 0,
                },
                "campaign_status": "ERROR",
                "campaign_error": "clean build failed; no test runner was executed",
                "results": [],
            }
            with (output_root / "run_all_summary.json").open(
                "w", encoding="utf-8"
            ) as stream:
                json.dump(summary, stream, indent=2, sort_keys=True)
                stream.write("\n")
            print(f"ERROR: clean build failed; see {build_result['log']}", file=sys.stderr)
            return 2

    results: list[dict[str, Any]] = []
    for test_id in selected:
        result = run_test(
            test_id, registry[test_id], build_dir, output_root, args.stream
        )
        results.append(result)
        print(f"{result['status']:7s} {test_id}  log={result['log']}")

    counts = {status: 0 for status in VALID_STATUSES}
    for result in results:
        counts[result["status"]] += 1
    counts["executed"] = len(results)
    print_summary(counts, results)

    campaign_status = "PASS"
    if counts["ERROR"]:
        campaign_status = "ERROR"
    elif counts["FAIL"] or (args.strict_skips and counts["SKIPPED"]):
        campaign_status = "FAIL"

    summary = {
        "source_revision": git_revision(),
        "build": build_result,
        "build_reused": args.skip_build,
        "build_directory": str(build_dir),
        "output_directory": str(output_root),
        "duration_seconds": time.monotonic() - campaign_started,
        "counts": counts,
        "campaign_status": campaign_status,
        "strict_skips": args.strict_skips,
        "results": results,
    }
    summary_path = output_root / "run_all_summary.json"
    with summary_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(f"\nAggregate result: {summary_path}")

    if counts["ERROR"]:
        return 2
    if counts["FAIL"] or (args.strict_skips and counts["SKIPPED"]):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Global runner for AMPS core unit and module tests.

The runner is intentionally standard-library only and Python-3.8 compatible.
It validates the explicit registry and component inventory, resolves portable
machine configuration, gives each test an isolated work/artifact directory,
and writes JSON plus JUnit evidence.  Individual runners communicate through
a small versioned context/result protocol; they never need host-specific paths.
"""

from __future__ import annotations

import argparse
import datetime as dt
import fnmatch
import hashlib
import json
import os
import queue
import random
import signal
import shutil
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


TEST_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TEST_ROOT.parents[1]
sys.path.insert(0, str(TEST_ROOT))
# The source checkout is read-mostly and may be mounted read-only on cluster
# workers. Prevent the control plane itself from leaving __pycache__ state;
# child adapters receive the matching environment setting in run_one().
sys.dont_write_bytecode = True

from lib.configuration import (  # noqa: E402
    ConfigurationError, command_available, load_configuration)
from lib.manifest import (  # noqa: E402
    CAMPAIGNS, ManifestError, audit_inventory, incomplete_components,
    load_manifest)


STATUSES = ("PASS", "FAIL", "SKIP", "ERROR")
# Severity is intentionally not the process-exit-code order: SKIP is a valid
# outcome and exits zero, but it still dominates PASS when a suite summarizes
# individual cases.  FAIL means the tested behavior was observed and was
# wrong; ERROR means the harness could not make that observation reliably.
SEVERITY = {"PASS": 0, "SKIP": 1, "FAIL": 2, "ERROR": 3}
EXPECTED_EXIT = {"PASS": 0, "SKIP": 0, "FAIL": 1, "ERROR": 2}


def stamp() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def one_line(value: Any, limit: int = 500) -> str:
    text = " ".join(str(value or "").split())
    return text if len(text) <= limit else text[:limit - 3] + "..."


def write_json(path: Path, value: Any) -> None:
    """Atomically publish JSON so a killed nightly never leaves valid-looking partial evidence."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(str(temporary), str(path))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run registered AMPS core unit/module tests with portable settings")
    selection = parser.add_argument_group("selection")
    selection.add_argument("--campaign", choices=CAMPAIGNS,
                           help="test campaign (default: smoke)")
    selection.add_argument("--tier", choices=CAMPAIGNS,
                           help="backward-compatible alias for --campaign")
    selection.add_argument("--all", action="store_true", help="select every registered test")
    selection.add_argument("--test", action="append", default=[], metavar="ID",
                           help="select an exact test ID; may be repeated")
    selection.add_argument("--match", action="append", default=[], metavar="GLOB",
                           help="select test IDs matching a shell glob")
    selection.add_argument("--component", action="append", default=[], metavar="PATH",
                           help="select a component or component prefix")
    selection.add_argument("--case", action="append", default=[], metavar="ID",
                           help="select an exact individual case ID; may be repeated")
    selection.add_argument("--order", choices=("normal", "reverse", "shuffle"),
                           default="normal", help="suite and case execution order")
    selection.add_argument("--shuffle-seed", type=int, default=0, metavar="N",
                           help="deterministic seed for --order shuffle (default: 0)")
    selection.add_argument("--list", action="store_true", help="list selected tests and exit")
    selection.add_argument("--coverage", action="store_true",
                           help="print the component coverage inventory and exit")

    machine = parser.add_argument_group("machine configuration")
    machine.add_argument("--config", action="append", type=Path, default=[], metavar="JSON",
                         help="overlay machine configuration; may be repeated")
    machine.add_argument("--output-dir", metavar="PATH",
                         help="exact output directory for this run")
    machine.add_argument("--build-jobs", type=int, metavar="N")
    machine.add_argument("--mpi-ranks", type=int, metavar="N")

    execution = parser.add_argument_group("execution and evidence")
    execution.add_argument(
        "--jobs", type=int, default=1, metavar="N",
        help=("maximum concurrently executing nonexclusive suites "
              "(default: 1; exclusive native suites always run alone)"))
    execution.add_argument("--fail-fast", action="store_true")
    execution.add_argument("--dry-run", action="store_true")
    execution.add_argument("--require-complete-coverage", action="store_true",
                           help="fail while any component is partial or planned")
    execution.add_argument("--json", type=Path, metavar="FILE",
                           help="additional JSON report path")
    execution.add_argument("--junit", type=Path, metavar="FILE",
                           help="additional JUnit report path")
    return parser.parse_args()


def select_tests(tests: Sequence[Mapping[str, Any]], args: argparse.Namespace) -> List[Dict[str, Any]]:
    """Select in manifest order; simultaneous explicit selectors are intersections.

    Manifest order is stable evidence order.  For example, ``--component
    src/pic --match '*BUFFER*'`` narrows to PIC buffer tests; it does not form
    a union that unexpectedly runs every test matching either selector.
    """

    requested_tests = list(getattr(args, "test", []))
    requested_matches = list(getattr(args, "match", []))
    requested_components = list(getattr(args, "component", []))
    requested_cases = list(getattr(args, "case", []))
    campaign = getattr(args, "campaign", None) or getattr(args, "tier", None) or "smoke"
    tier = getattr(args, "tier", None)
    if getattr(args, "campaign", None) and tier and campaign != tier:
        raise ManifestError("--campaign and --tier name different campaigns")
    if len(set(requested_cases)) != len(requested_cases):
        raise ManifestError("duplicate requested case IDs are not allowed")
    explicit_suite = bool(requested_tests or requested_matches or requested_components)
    selected: List[Dict[str, Any]] = []
    known_tests = set()
    known_cases = set()
    for item in tests:
        known_tests.add(str(item["id"]))
        known_tests.update(str(value) for value in item.get("aliases", []))
        known_cases.update(str(value) for value in item["expected_cases"])
    unknown = set(requested_tests) - known_tests
    if unknown:
        raise ManifestError("unknown requested test IDs: {}".format(", ".join(sorted(unknown))))
    unknown_cases = set(requested_cases) - known_cases
    if unknown_cases:
        raise ManifestError("unknown requested case IDs: {}".format(
            ", ".join(sorted(unknown_cases))))
    for item in tests:
        identifier = str(item["id"])
        suite_names = {identifier}
        suite_names.update(str(value) for value in item.get("aliases", []))
        alias_cases = item.get("alias_cases", {})
        component = str(item["component"])
        owned_cases = set(str(value) for value in item["expected_cases"])
        alias_subset: Optional[set] = None
        if requested_tests and identifier not in requested_tests:
            matched_aliases = suite_names.intersection(requested_tests) - {identifier}
            if matched_aliases and all(alias in alias_cases for alias in matched_aliases):
                alias_subset = set(
                    value for alias in matched_aliases for value in alias_cases[alias])
        elif requested_matches and not any(
                fnmatch.fnmatchcase(identifier, pattern)
                for pattern in requested_matches):
            matched_aliases = {
                alias for alias in item.get("aliases", [])
                if any(fnmatch.fnmatchcase(alias, pattern)
                       for pattern in requested_matches)}
            if matched_aliases and all(alias in alias_cases for alias in matched_aliases):
                alias_subset = set(
                    value for alias in matched_aliases for value in alias_cases[alias])
        if getattr(args, "all", False):
            include = not requested_cases or bool(owned_cases & set(requested_cases))
        elif explicit_suite:
            include = True
            if requested_tests:
                include = include and bool(suite_names & set(requested_tests))
            if requested_matches:
                include = include and any(
                    fnmatch.fnmatchcase(name, pattern)
                    for name in suite_names for pattern in requested_matches)
            if requested_components:
                include = include and any(
                    component == prefix.rstrip("/") or
                    component.startswith(prefix.rstrip("/") + "/")
                    for prefix in requested_components)
            if requested_cases:
                include = include and bool(owned_cases & set(requested_cases))
                if alias_subset is not None:
                    include = include and bool(alias_subset & set(requested_cases))
        elif requested_cases:
            include = bool(owned_cases & set(requested_cases))
        else:
            include = campaign in item.get("campaigns", item["tiers"])
        if include:
            selected_item = dict(item)
            if requested_cases:
                selected_item["selected_cases"] = [
                    value for value in requested_cases
                    if value in owned_cases and
                    (alias_subset is None or value in alias_subset)]
            elif alias_subset is not None:
                selected_item["selected_cases"] = [
                    value for value in item["expected_cases"]
                    if value in alias_subset]
            else:
                selected_item["selected_cases"] = list(item["expected_cases"])
            selected.append(selected_item)
    return selected


def order_tests(tests: Sequence[Mapping[str, Any]], order: str,
                shuffle_seed: int = 0) -> List[Dict[str, Any]]:
    """Return a non-mutating deterministic suite/case execution order."""

    ordered = [dict(item) for item in tests]
    for item in ordered:
        item["selected_cases"] = list(item.get("selected_cases", item["expected_cases"]))
        item["case_order"] = order
        item["shuffle_seed"] = int(shuffle_seed)
    if order == "normal":
        return ordered
    if order == "reverse":
        for item in ordered:
            item["selected_cases"].reverse()
        ordered.reverse()
        return ordered
    if order != "shuffle":
        raise ManifestError("unsupported execution order {!r}".format(order))
    generator = random.Random(int(shuffle_seed))
    for item in ordered:
        generator.shuffle(item["selected_cases"])
    generator.shuffle(ordered)
    return ordered


def print_tests(tests: Sequence[Mapping[str, Any]]) -> None:
    for item in tests:
        print("{:<30} {:<12} {:<30} cases={:<3} {}".format(
            item["id"], item["level"], item["component"],
            len(item.get("selected_cases", item["expected_cases"])), item["description"]))


def print_coverage(entries: Sequence[Mapping[str, Any]]) -> None:
    for item in entries:
        print("{:<38} {:<8} {}".format(
            item.get("component", "?"), item.get("status", "?"),
            ",".join(item.get("test_ids", [])) or "-"))
    incomplete = incomplete_components(entries)
    print("coverage: components={} complete={} incomplete={}".format(
        len(entries), len(entries) - len(incomplete), len(incomplete)))


def make_output_directory(config: Mapping[str, Any], requested: Optional[str]) -> Path:
    """Resolve the evidence root without requiring it to live in the repository."""

    if requested:
        candidate = Path(os.path.expandvars(os.path.expanduser(requested)))
        output = candidate if candidate.is_absolute() else REPO_ROOT / candidate
    else:
        output = Path(config["paths"]["output_root"]) / ("run-{}-{}".format(stamp(), os.getpid()))
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    return output


def missing_requirements(test: Mapping[str, Any], config: Mapping[str, Any]) -> List[str]:
    """Probe declared prerequisites before creating a misleading failed test.

    An unavailable compiler/reference fixture is SKIP, not a scientific FAIL.
    Paths are constrained to the repository because machine-local datasets are
    passed through the configured data root rather than manifest path entries.
    """

    missing: List[str] = []
    for name in test["requires"]["tools"]:
        configured = config["tools"].get(name)
        if configured is None:
            missing.append("tool '{}' is not configured".format(name))
        elif not command_available(configured, REPO_ROOT):
            missing.append("tool '{}' is unavailable: {}".format(name, configured[0]))
    for relative in test["requires"]["paths"]:
        candidate = (REPO_ROOT / relative).resolve()
        try:
            candidate.relative_to(REPO_ROOT)
        except ValueError:
            missing.append("required path escapes repository: {}".format(relative))
            continue
        if not candidate.exists():
            missing.append("required path is absent: {}".format(relative))
    capabilities = config.get("capabilities", {})
    for name in test["requires"].get("capabilities", []):
        declaration = capabilities.get(name, {})
        if not bool(declaration.get("available", False)):
            policy = declaration.get("policy", "optional")
            missing.append("{} capability '{}' is unavailable".format(policy, name))
    return missing


def required_capability_missing(reasons: Sequence[str]) -> bool:
    return any(reason.startswith("required capability ") for reason in reasons)


def skipped_result(test: Mapping[str, Any], reason: str, log: Path) -> Dict[str, Any]:
    return {
        "id": test["id"], "component": test["component"], "level": test["level"],
        "status": "SKIP", "message": one_line(reason), "seconds": 0.0,
        "runner": test["runner"], "log": str(log), "work_dir": "",
        "artifact_dir": "", "cases": [], "evidence": [], "metrics": {},
    }


def prerequisite_error_result(test: Mapping[str, Any], reason: str,
                              log: Path) -> Dict[str, Any]:
    result = skipped_result(test, reason, log)
    result["status"] = "ERROR"
    return result


def validate_result(test: Mapping[str, Any], payload: Any, exit_code: int) -> Tuple[str, str, List[Dict[str, Any]]]:
    """Enforce the versioned runner protocol and prevent false-green suites.

    The suite status must be at least as severe as every reported case, IDs
    must be unique, and status must agree with the child process exit code.
    Thus neither a buggy runner nor a wrapper script can hide a failed case.
    """

    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        return "ERROR", "runner result is not a schema_version 1 object", []
    if payload.get("id") != test["id"]:
        return "ERROR", "runner result ID {!r} does not match {}".format(
            payload.get("id"), test["id"]), []
    status = str(payload.get("status", "")).upper()
    if status not in STATUSES:
        return "ERROR", "runner returned invalid status {!r}".format(status), []
    raw_cases = payload.get("cases", [])
    if not isinstance(raw_cases, list):
        return "ERROR", "runner cases must be a list", []
    cases: List[Dict[str, Any]] = []
    case_ids = set()
    for index, raw in enumerate(raw_cases):
        if not isinstance(raw, dict) or not isinstance(raw.get("id"), str):
            return "ERROR", "runner case {} is malformed".format(index), []
        case_status = str(raw.get("status", "")).upper()
        if case_status not in STATUSES:
            return "ERROR", "runner case {} has invalid status".format(raw.get("id")), []
        if raw["id"] in case_ids:
            return "ERROR", "runner returned duplicate case {}".format(raw["id"]), []
        case_ids.add(raw["id"])
        cases.append({"id": raw["id"], "status": case_status,
                      "message": one_line(raw.get("message", ""))})
    if status == "PASS" and not cases:
        return "ERROR", "runner returned PASS without any executed cases", []
    expected_cases = list(test.get("selected_cases", test.get("expected_cases", [])))
    actual_cases = [item["id"] for item in cases]
    if expected_cases and actual_cases != expected_cases:
        return "ERROR", "runner case registry mismatch: expected {}, observed {}".format(
            expected_cases, actual_cases), cases
    if cases:
        worst = max(cases, key=lambda item: SEVERITY[item["status"]])["status"]
        if SEVERITY[worst] > SEVERITY[status]:
            return "ERROR", "suite status {} hides a {} case".format(status, worst), cases
    expected = EXPECTED_EXIT[status]
    if exit_code != expected:
        return "ERROR", "runner status {} requires exit {}, got {}".format(
            status, expected, exit_code), cases
    return status, one_line(payload.get("message", "")), cases


def run_one(test: Mapping[str, Any], config: Mapping[str, Any], output: Path,
            dry_run: bool, case_prefix: str = "    ") -> Dict[str, Any]:
    """Execute one runner in isolated state and normalize its retained evidence.

    Runners get paths and commands only through the context receipt.  They may
    create data below ``work_dir`` and ``artifact_dir`` but cannot choose the
    global report location.  This separation lets work use node-local scratch
    while artifacts remain on durable CI storage.
    """

    identifier = str(test["id"])
    log = output / "logs" / (identifier + ".log")
    log.parent.mkdir(parents=True, exist_ok=True)
    missing = missing_requirements(test, config)
    if missing:
        message = "; ".join(missing)
        if required_capability_missing(missing):
            log.write_text("ERROR: {}\n".format(message), encoding="utf-8")
            return prerequisite_error_result(test, message, log)
        log.write_text("SKIP: {}\n".format(message), encoding="utf-8")
        return skipped_result(test, message, log)

    temporary_root = str(config["paths"].get("temporary_root", ""))
    if temporary_root:
        # output.name alone is not unique: two simultaneous CI jobs can use
        # different output parents with the same basename. Bind scratch state
        # to the complete resolved evidence path while keeping names readable.
        run_key = hashlib.sha256(str(output).encode("utf-8")).hexdigest()[:12]
        work = (Path(temporary_root) / "amps-core-tests" /
                (output.name + "-" + run_key) / identifier)
    else:
        work = output / "work" / identifier
    artifacts = output / "artifacts" / identifier
    work.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)
    result_file = work / "result.json"
    context_file = work / "context.json"
    # An explicit --output-dir may intentionally be reused. Invalidate the
    # runner-owned receipt before launch so a crashed/buggy adapter cannot have
    # an old PASS result accepted as evidence for this invocation.
    try:
        result_file.unlink()
    except FileNotFoundError:
        pass
    # The context is the only supported host/application interface for an
    # individual runner.  Keeping it versioned makes old runners fail clearly
    # after an incompatible schema change rather than guessing new defaults.
    context = {
        "schema_version": 1,
        "test": {key: value for key, value in test.items() if key != "runner_path"},
        "repo_root": str(REPO_ROOT),
        "test_root": str(TEST_ROOT),
        "work_dir": str(work),
        "artifact_dir": str(artifacts),
        "result_file": str(result_file),
        "paths": config["paths"],
        "tools": config["tools"],
        "execution": config["execution"],
        "capabilities": config.get("capabilities", {}),
    }
    write_json(context_file, context)
    runner_command = list(config["tools"]["python"]) + [
        str(test["runner_path"]), "--context", str(context_file)]
    if dry_run:
        log.write_text("DRY-RUN: {}\n".format(" ".join(runner_command)), encoding="utf-8")
        result = skipped_result(test, "dry run", log)
        result["work_dir"] = str(work)
        result["artifact_dir"] = str(artifacts)
        return result

    timeout = int(test.get("timeout_seconds") or config["execution"]["default_timeout_seconds"])
    started = time.monotonic()
    environment = dict(os.environ)
    environment["AMPS_TEST_CONTEXT"] = str(context_file)
    # Bytecode is disabled so executing a read-mostly source checkout does not
    # add __pycache__ state.  The explicit context variable also supports
    # diagnostics and future non-Python adapters; --context remains canonical.
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    # Stream the adapter so case-start records reach the user before their
    # assertions execute.  Other adapter/build output remains in the complete
    # retained log: native rebuilds can emit thousands of compiler lines that
    # would otherwise obscure the case progression requested by the user.
    process = subprocess.Popen(
        runner_command, cwd=str(REPO_ROOT), env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, bufsize=1,
        start_new_session=(os.name == "posix"))
    output_lines: List[str] = []
    events: "queue.Queue[object]" = queue.Queue()
    completed_marker = object()

    def read_runner_output() -> None:
        try:
            assert process.stdout is not None
            for line in process.stdout:
                events.put(line)
        finally:
            events.put(completed_marker)

    reader = threading.Thread(target=read_runner_output,
                              name="amps-global-test-output-reader")
    reader.daemon = True
    reader.start()
    deadline = started + timeout
    timed_out = False
    try:
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                timed_out = True
                break
            try:
                item = events.get(timeout=min(0.1, remaining))
            except queue.Empty:
                # An adapter may exit while an incorrectly detached descendant
                # still holds its stdout descriptor. Kill the adapter's
                # dedicated group so the reader reaches EOF instead of waiting
                # for timeout.
                if os.name == "posix" and process.poll() is not None:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                continue
            if item is completed_marker:
                break
            line = str(item)
            output_lines.append(line)
            if line.startswith("[CASE] RUN"):
                print("{}{}".format(case_prefix, line.rstrip()), flush=True)
    except BaseException:
        # KeyboardInterrupt and scheduler cancellation are not normal result
        # states, but they still must not orphan make/mpiexec descendants or
        # leave the native repository lock held. The adapter owns this process
        # group by construction; terminate it before re-raising to the caller.
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
        reader.join(timeout=1.0)
        raise

    if timed_out:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
        reader.join(timeout=1.0)
        while True:
            try:
                item = events.get_nowait()
            except queue.Empty:
                break
            if item is not completed_marker:
                output_lines.append(str(item))
        output_lines.append("\nERROR: runner timed out after {} seconds\n".format(timeout))
        exit_code = 124
    else:
        try:
            # EOF is not proof that a process exited: a faulty adapter can
            # close stdout and continue running. Preserve the same absolute
            # deadline while waiting for its status.
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                raise subprocess.TimeoutExpired(runner_command, timeout)
            exit_code = process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            timed_out = True
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            process.wait()
            output_lines.append(
                "\nERROR: runner timed out after {} seconds\n".format(timeout))
            exit_code = 124
        except BaseException:
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            process.wait()
            reader.join(timeout=1.0)
            raise
        if os.name == "posix":
            # Normally this group is already empty. This defensive cleanup
            # prevents a faulty adapter from leaving compiler/MPI descendants.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        reader.join(timeout=1.0)
    output_text = "".join(output_lines)
    seconds = time.monotonic() - started
    with log.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("command: {}\n".format(" ".join(runner_command)))
        stream.write("timeout_seconds: {}\n\n".format(timeout))
        stream.write(output_text)
        if output_text and not output_text.endswith("\n"):
            stream.write("\n")

    if timed_out:
        status, message, cases, payload = "ERROR", "runner timed out", [], {}
    elif not result_file.is_file():
        status, message, cases, payload = (
            "ERROR", "runner did not write {}".format(result_file), [], {})
    else:
        try:
            with result_file.open(encoding="utf-8") as stream:
                payload = json.load(stream)
            status, message, cases = validate_result(test, payload, exit_code)
        except (OSError, ValueError) as error:
            status, message, cases, payload = "ERROR", "cannot read runner result: {}".format(error), [], {}

    result = {
        "id": identifier, "component": test["component"], "level": test["level"],
        "status": status, "message": message, "seconds": seconds,
        "runner": test["runner"], "log": str(log), "work_dir": str(work),
        "artifact_dir": str(artifacts), "cases": cases,
        "evidence": payload.get("evidence", []) if isinstance(payload, dict) else [],
        "metrics": payload.get("metrics", {}) if isinstance(payload, dict) else {},
    }
    keep = config["execution"]["keep_work"]
    if keep == "never" or (keep == "failures" and status in ("PASS", "SKIP")):
        # The target was created by this invocation below the resolved run
        # directory; never accept a runner-supplied deletion path.
        shutil.rmtree(str(work))
        result["work_dir"] = ""
    return result


def counts(results: Iterable[Mapping[str, Any]]) -> Dict[str, int]:
    totals = {status: 0 for status in STATUSES}
    for item in results:
        totals[str(item["status"])] += 1
    totals["TOTAL"] = sum(totals.values())
    return totals


def individual_results(results: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """Flatten suite receipts into the final human/JSON result inventory.

    Successful adapters normally publish explicit cases. If a suite cannot
    start or cannot produce trustworthy case evidence, retain one suite-level
    row instead of silently omitting that selected test from the final report.
    This makes the terminal inventory reconcile exactly with its summary even
    for prerequisite skips, compile errors, and fail-fast runs.
    """

    rows: List[Dict[str, Any]] = []
    for result in results:
        cases = result.get("cases", [])
        if cases:
            for child in cases:
                rows.append({
                    "id": child["id"],
                    "suite": result["id"],
                    "component": result["component"],
                    "status": child["status"],
                    "message": one_line(child.get("message", "")),
                    "kind": "case",
                })
        else:
            rows.append({
                "id": result["id"],
                "suite": result["id"],
                "component": result["component"],
                "status": result["status"],
                "message": one_line(result.get("message", "")),
                "kind": "suite",
            })
    return rows


def print_final_results(rows: Sequence[Mapping[str, Any]]) -> None:
    """Print every final outcome after all interleaved execution has ended."""

    print("FINAL INDIVIDUAL TEST RESULTS")
    for row in rows:
        suffix = ""
        if row.get("kind") == "suite":
            suffix = " [suite-level: no individual case evidence]"
        if row["status"] != "PASS" and row.get("message"):
            suffix += " {}".format(row["message"])
        print("{:<5} {:<42} suite={}{}".format(
            row["status"], row["id"], row["suite"], suffix))


def execution_batches(
    tests: Sequence[Mapping[str, Any]],
) -> List[List[Tuple[int, Mapping[str, Any]]]]:
    """Place exclusive suites behind barriers while retaining manifest order.

    A native AMPS adapter changes the checkout-wide generated configuration,
    so it must neither overlap another native adapter nor a dependency-light
    suite that could inspect the tree. Consecutive nonexclusive suites form a
    batch that the caller may execute concurrently. Each exclusive suite is a
    one-element batch and therefore has both a leading and trailing barrier.
    """

    batches: List[List[Tuple[int, Mapping[str, Any]]]] = []
    parallel: List[Tuple[int, Mapping[str, Any]]] = []
    for index, test in enumerate(tests):
        if bool(test.get("exclusive", False)):
            if parallel:
                batches.append(parallel)
                parallel = []
            batches.append([(index, test)])
        else:
            parallel.append((index, test))
    if parallel:
        batches.append(parallel)
    return batches


def run_batch(
    batch: Sequence[Tuple[int, Mapping[str, Any]]],
    selected_count: int,
    config: Mapping[str, Any],
    output: Path,
    dry_run: bool,
    jobs: int,
    fail_fast: bool,
) -> List[Tuple[int, Dict[str, Any]]]:
    """Run one barrier-delimited batch and return results in manifest order."""

    workers = min(jobs, len(batch))
    if workers <= 1:
        completed = []
        for index, test in batch:
            number = index + 1
            print("[{}/{}] RUN  {} ({})".format(
                number, selected_count, test["id"], test["component"]), flush=True)
            result = run_one(test, config, output, dry_run)
            print("[{}/{}] {:5} {} {:.3f}s {}".format(
                number, selected_count, result["status"], result["id"],
                result["seconds"], result["message"]), flush=True)
            completed.append((index, result))
            if fail_fast and result["status"] in ("FAIL", "ERROR"):
                break
        return completed

    futures: Dict[Future, Tuple[int, Mapping[str, Any]]] = {}
    with ThreadPoolExecutor(max_workers=workers,
                            thread_name_prefix="amps-core-suite") as executor:
        pending = iter(batch)

        def submit(index: int, test: Mapping[str, Any]) -> None:
            number = index + 1
            print("[{}/{}] RUN  {} ({})".format(
                number, selected_count, test["id"], test["component"]), flush=True)
            # Prefix streamed case starts with the owning suite because output
            # from independent adapters may legitimately interleave.
            case_prefix = "    [{}] ".format(test["id"])
            future = executor.submit(
                run_one, test, config, output, dry_run, case_prefix)
            futures[future] = (index, test)

        for _ in range(workers):
            try:
                index, test = next(pending)
            except StopIteration:
                break
            submit(index, test)

        completed = []
        stop_submitting = False
        while futures:
            done, _ = wait(list(futures), return_when=FIRST_COMPLETED)
            for future in done:
                index, test = futures.pop(future)
                number = index + 1
                try:
                    result = future.result()
                except Exception as error:
                    # run_one normally translates child failures into evidence.
                    # This boundary catches control-plane defects without losing
                    # the other concurrently completed receipts.
                    result = {
                        "id": test["id"], "component": test["component"],
                        "level": test["level"], "status": "ERROR",
                        "message": one_line("scheduler worker failed: {}".format(error)),
                        "seconds": 0.0, "runner": test["runner"], "log": "",
                        "work_dir": "", "artifact_dir": "", "cases": [],
                        "evidence": [], "metrics": {},
                    }
                print("[{}/{}] {:5} {} {:.3f}s {}".format(
                    number, selected_count, result["status"], result["id"],
                    result["seconds"], result["message"]), flush=True)
                completed.append((index, result))
                if fail_fast and result["status"] in ("FAIL", "ERROR"):
                    stop_submitting = True
            if not stop_submitting:
                for _ in range(len(done)):
                    try:
                        next_index, next_test = next(pending)
                    except StopIteration:
                        break
                    submit(next_index, next_test)
    return sorted(completed, key=lambda item: item[0])


def write_junit(path: Path, results: Sequence[Mapping[str, Any]], elapsed: float) -> None:
    """Flatten suites/cases into portable JUnit while retaining the full log path."""

    flat_count = sum(max(1, len(item.get("cases", []))) for item in results)
    failures = sum(1 for item in results for case in (item.get("cases") or [item])
                   if case["status"] == "FAIL")
    errors = sum(1 for item in results for case in (item.get("cases") or [item])
                 if case["status"] == "ERROR")
    skipped = sum(1 for item in results for case in (item.get("cases") or [item])
                  if case["status"] == "SKIP")
    suite = ET.Element("testsuite", {
        "name": "amps-core", "tests": str(flat_count), "failures": str(failures),
        "errors": str(errors), "skipped": str(skipped), "time": "{:.6f}".format(elapsed),
    })
    for result in results:
        children = result.get("cases") or [{
            "id": result["id"], "status": result["status"], "message": result["message"]}]
        per_case_time = float(result["seconds"]) / max(1, len(children))
        for child in children:
            case_node = ET.SubElement(suite, "testcase", {
                "classname": str(result["component"]),
                "name": "{}/{}".format(result["id"], child["id"]),
                "time": "{:.6f}".format(per_case_time),
            })
            status = child["status"]
            message = one_line(child.get("message", ""))
            if status == "FAIL":
                ET.SubElement(case_node, "failure", {"message": message}).text = message
            elif status == "ERROR":
                ET.SubElement(case_node, "error", {"message": message}).text = message
            elif status == "SKIP":
                ET.SubElement(case_node, "skipped", {"message": message})
            ET.SubElement(case_node, "system-out").text = "log: {}".format(result["log"])
    path.parent.mkdir(parents=True, exist_ok=True)
    tree = ET.ElementTree(suite)
    tree.write(str(path), encoding="utf-8", xml_declaration=True)


def copy_report(source: Path, requested: Optional[Path]) -> None:
    if requested is None:
        return
    target = requested if requested.is_absolute() else REPO_ROOT / requested
    target = target.resolve()
    if target == source.resolve():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(str(source), str(target))


def main() -> int:
    args = parse_arguments()
    if args.jobs <= 0:
        print("configuration error: --jobs must be a positive integer", file=sys.stderr)
        return 2
    try:
        config = load_configuration(
            REPO_ROOT, TEST_ROOT / "config" / "default.json", args.config,
            output_root=args.output_dir, build_jobs=args.build_jobs,
            mpi_ranks=args.mpi_ranks)
        if args.campaign is None and args.tier is None:
            args.campaign = config["execution"]["default_campaign"]
        tests = load_manifest(TEST_ROOT / "manifest.json", REPO_ROOT, TEST_ROOT)
        inventory, inventory_errors = audit_inventory(
            TEST_ROOT / "component_inventory.json", REPO_ROOT, tests)
        if inventory_errors:
            raise ManifestError("; ".join(inventory_errors))
        selected = order_tests(
            select_tests(tests, args), args.order, args.shuffle_seed)
    except (ConfigurationError, ManifestError) as error:
        print("configuration error: {}".format(error), file=sys.stderr)
        return 2

    if args.coverage:
        print_coverage(inventory)
        if args.require_complete_coverage and incomplete_components(inventory):
            return 1
        return 0
    if args.list:
        print_tests(selected)
        return 0
    if not selected:
        print("no tests selected", file=sys.stderr)
        return 2

    output = make_output_directory(config, args.output_dir)
    write_json(output / "effective-config.json", config)
    campaign = args.campaign or args.tier or "smoke"
    selected_case_ids = [
        case_id for item in selected for case_id in item["selected_cases"]]
    write_json(output / "selection.json", {
        "schema_version": 1,
        "campaign": campaign,
        "order": args.order,
        "shuffle_seed": args.shuffle_seed,
        "case_ids": selected_case_ids,
        "tests": selected,
    })
    print("AMPS core tests: campaign={} selected_suites={} selected_cases={} order={} output={}".format(
        campaign, len(selected), len(selected_case_ids), args.order, output))
    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    started = time.monotonic()
    indexed_results: List[Tuple[int, Dict[str, Any]]] = []
    for batch in execution_batches(selected):
        batch_results = run_batch(
            batch, len(selected), config, output, args.dry_run, args.jobs,
            args.fail_fast)
        indexed_results.extend(batch_results)
        # Concurrent work already submitted in this barrier batch is allowed
        # to finish. Fail-fast prevents scheduling subsequent batches; it
        # cannot safely revoke already-running adapters and their evidence.
        if (args.fail_fast and any(
                result["status"] in ("FAIL", "ERROR")
                for _, result in batch_results)):
            break
    results = [result for _, result in sorted(indexed_results, key=lambda item: item[0])]

    elapsed = time.monotonic() - started
    totals = counts(results)
    final_individual_results = individual_results(results)
    individual_totals = counts(final_individual_results)
    incomplete = incomplete_components(inventory)
    report = {
        "schema_version": 1,
        "started_utc": started_utc,
        "elapsed_seconds": elapsed,
        "output_directory": str(output),
        "campaign": campaign,
        "order": args.order,
        "shuffle_seed": args.shuffle_seed,
        "selection": [item["id"] for item in selected],
        "case_selection": selected_case_ids,
        "not_executed_suites": [
            item["id"] for item in selected
            if item["id"] not in {result["id"] for result in results}],
        "counts": totals,
        "individual_counts": individual_totals,
        "individual_results": final_individual_results,
        "coverage": {
            "component_count": len(inventory),
            "complete": not incomplete,
            "incomplete_components": incomplete,
        },
        "results": results,
    }
    json_report = output / "amps-core-tests.json"
    junit_report = output / "amps-core-tests.xml"
    write_json(json_report, report)
    write_junit(junit_report, results, elapsed)
    copy_report(json_report, args.json)
    copy_report(junit_report, args.junit)

    print_final_results(final_individual_results)
    print("INDIVIDUAL SUMMARY total={TOTAL} pass={PASS} fail={FAIL} skip={SKIP} error={ERROR}".format(
        **individual_totals))
    print("SUITE SUMMARY total={TOTAL} pass={PASS} fail={FAIL} skip={SKIP} error={ERROR} seconds={:.3f}".format(
        elapsed, **totals))
    print("JSON:  {}".format(json_report))
    print("JUnit: {}".format(junit_report))
    # SKIP is visible in both reports but does not make a heterogeneous-machine
    # run fail.  A nightly policy that requires a capability can detect the
    # nonzero skip count or provision that prerequisite explicitly.
    failed = totals["FAIL"] or totals["ERROR"]
    if args.require_complete_coverage and incomplete:
        print("COVERAGE FAIL: incomplete components: {}".format(", ".join(incomplete)))
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

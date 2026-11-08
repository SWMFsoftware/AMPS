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
import socket
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
# Stream workers can finish on different threads.  Serialize complete terminal
# lines while leaving suite logs independent so a user never sees two stream
# status records spliced together.
OUTPUT_LOCK = threading.Lock()


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
        description="Run registered AMPS core and application tests with portable settings")
    selection = parser.add_argument_group("selection")
    selection.add_argument("--campaign", choices=CAMPAIGNS,
                           help="test campaign (default: smoke)")
    selection.add_argument("--tier", choices=CAMPAIGNS,
                           help="backward-compatible alias for --campaign")
    selection.add_argument("--all", action="store_true", help="select every registered test")
    selection.add_argument(
        "--full-validation", action="store_true",
        help=("select the complete registry and request expensive, portability, "
              "distributed, and application qualification execution"))
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
    selection.add_argument(
        "--list", action="store_true",
        help="list the complete registry, or an explicitly narrowed selection, by stream")
    selection.add_argument("--coverage", action="store_true",
                           help="print the component coverage inventory and exit")

    machine = parser.add_argument_group("machine configuration")
    machine.add_argument(
        "--config", action="append", type=Path, default=[], metavar="JSON",
        help=("required machine configuration overlay; may be repeated "
              "(profiles: test/qualification-machine/)"))
    machine.add_argument("--output-dir", metavar="PATH",
                         help="exact output directory for this run")
    machine.add_argument("--build-jobs", type=int, metavar="N")
    machine.add_argument("--mpi-ranks", type=int, metavar="N")

    execution = parser.add_argument_group("execution and evidence")
    execution.add_argument(
        "--jobs", type=int, default=1, metavar="N",
        help=("maximum concurrently executing nonexclusive suites "
              "inside each stream (default: 1)"))
    execution.add_argument(
        "--streams", type=int, default=1, metavar="N",
        help=("maximum independently executing core/application test streams "
              "(default: 1; checkout-mutating suites remain globally serialized)"))
    execution.add_argument("--fail-fast", action="store_true")
    execution.add_argument("--dry-run", action="store_true")
    execution.add_argument("--require-complete-coverage", action="store_true",
                           help="fail while any component is partial or planned")
    execution.add_argument("--json", type=Path, metavar="FILE",
                           help="additional JSON report path")
    execution.add_argument("--junit", type=Path, metavar="FILE",
                           help="additional JUnit report path")
    arguments = parser.parse_args()
    if not arguments.config:
        profiles = REPO_ROOT / "test" / "qualification-machine"
        # Match the make wrapper's ``hostname -s`` convention.  The expected
        # filename is diagnostic guidance rather than an implicit fallback:
        # an operator must still name the reviewed profile explicitly, which
        # prevents a costly campaign from silently using another host's tools,
        # MPI launcher, scratch filesystem, or qualification capabilities.
        short_hostname = socket.gethostname().split(".", 1)[0]
        parser.error(
            "--config is required; select a machine profile from {} "
            "(for this host use {}, or copy and edit {})".format(
                profiles, profiles / (short_hostname + ".json"),
                profiles / "template.json"))
    return arguments


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
    full_validation = bool(getattr(args, "full_validation", False))
    campaign = getattr(args, "campaign", None) or getattr(args, "tier", None) or "smoke"
    tier = getattr(args, "tier", None)
    if getattr(args, "campaign", None) and tier and campaign != tier:
        raise ManifestError("--campaign and --tier name different campaigns")
    if full_validation and (
            getattr(args, "all", False) or requested_tests or
            requested_matches or requested_components or requested_cases or
            getattr(args, "campaign", None) or tier):
        raise ManifestError(
            "--full-validation cannot be combined with campaign or narrowing selectors")
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
        if full_validation:
            # Full validation is deliberately not a synonym for a particular
            # campaign.  It is the union of the complete reviewed registry and
            # also propagates qualification intent through the execution
            # context below.  Keeping selection here in manifest order gives
            # the final combined report the same stable evidence ordering as
            # every other global run.
            include = True
        elif getattr(args, "all", False):
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


def stream_id_for_test(test: Mapping[str, Any]) -> str:
    """Return the stable execution-stream identity for a registered suite.

    Core ownership paths use the ``src/...`` form validated by the manifest.
    Explicitly allow-listed applications use top-level names such as
    ``srcEnceladus`` and ``srcMoon``.  Deriving streams from that already
    validated distinction prevents a second application registry from drifting
    away from the manifest while keeping each application's logs independent.
    """

    component = str(test["component"])
    return "core" if component.startswith("src/") else "application-" + component


def group_test_streams(
    tests: Sequence[Mapping[str, Any]],
) -> List[Dict[str, Any]]:
    """Group selected suites into core and per-application streams.

    Stream order is deliberately stable: core first, followed by application
    components in lexical order.  Suite order inside each stream remains the
    exact selected order (normal, reverse, or seeded shuffle).  Original
    selection indices are retained so the combined report can be restored to
    the requested logical order regardless of concurrent completion order.
    """

    grouped: Dict[str, Dict[str, Any]] = {}
    for index, test in enumerate(tests):
        stream_id = stream_id_for_test(test)
        if stream_id not in grouped:
            component = "" if stream_id == "core" else str(test["component"])
            grouped[stream_id] = {
                "id": stream_id,
                "kind": "core" if stream_id == "core" else "application",
                "component": component,
                "tests": [],
            }
        grouped[stream_id]["tests"].append((index, test))
    order = (["core"] if "core" in grouped else []) + sorted(
        identifier for identifier in grouped if identifier != "core")
    return [grouped[identifier] for identifier in order]


def print_tests(
    tests: Sequence[Mapping[str, Any]],
    registered_tests: Optional[Sequence[Mapping[str, Any]]] = None,
    stream_jobs: int = 1,
    suite_jobs: int = 1,
) -> None:
    """Print registered applications and the exact selected stream plan."""

    registry = list(registered_tests if registered_tests is not None else tests)
    applications: Dict[str, List[Mapping[str, Any]]] = {}
    selected_ids = {str(item["id"]) for item in tests}
    for item in registry:
        component = str(item["component"])
        if not component.startswith("src/"):
            applications.setdefault(component, []).append(item)
    print("REGISTERED APPLICATIONS count={}".format(len(applications)))
    for component in sorted(applications):
        suites = applications[component]
        selected_count = sum(
            1 for item in suites if str(item["id"]) in selected_ids)
        print("  {:<20} registered_suites={:<3} registered_cases={:<4} selected_suites={}".format(
            component, len(suites),
            sum(len(item["expected_cases"]) for item in suites), selected_count))

    streams = group_test_streams(tests)
    print("EXECUTION STREAMS selected={} concurrent_limit={} suites_per_stream={}".format(
        len(streams), min(stream_jobs, len(streams)) if streams else 0,
        suite_jobs))
    for number, stream in enumerate(streams, 1):
        indexed = stream["tests"]
        case_count = sum(len(item.get("selected_cases", item["expected_cases"]))
                         for _, item in indexed)
        application = (" component={}".format(stream["component"])
                       if stream["kind"] == "application" else "")
        print("\nSTREAM {}/{} id={} kind={}{} suites={} cases={}".format(
            number, len(streams), stream["id"], stream["kind"], application,
            len(indexed), case_count))
        for _, item in indexed:
            print("  {:<30} {:<12} {:<30} cases={:<3} {}".format(
                item["id"], item["level"], item["component"],
                len(item.get("selected_cases", item["expected_cases"])),
                item["description"]))


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
    """Place one stream's exclusive suites behind order-preserving barriers."""

    return indexed_execution_batches(list(enumerate(tests)))


def indexed_execution_batches(
    tests: Sequence[Tuple[int, Mapping[str, Any]]],
) -> List[List[Tuple[int, Mapping[str, Any]]]]:
    """Batch indexed suites without losing their combined-selection indices.

    A native AMPS adapter changes the checkout-wide generated configuration,
    so it remains a one-element barrier inside its stream.  ``CheckoutGate``
    below extends that barrier across all concurrently executing streams.
    """

    batches: List[List[Tuple[int, Mapping[str, Any]]]] = []
    parallel: List[Tuple[int, Mapping[str, Any]]] = []
    for index, test in tests:
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


class CheckoutGate:
    """Writer-priority gate protecting checkout-mutating native suites.

    Nonexclusive suites enter as readers and may overlap across streams.
    An exclusive suite enters as a writer only after every reader finishes;
    once a writer is waiting, new readers pause so a busy core stream cannot
    starve an application build.  This is stricter than merely placing a mutex
    around builds: no dependency-light reader may inspect generated headers
    while another stream changes the repository-wide AMPS configuration.
    """

    def __init__(self) -> None:
        self._condition = threading.Condition()
        self._readers = 0
        self._writer = False
        self._waiting_writers = 0

    def acquire(self, exclusive: bool) -> None:
        with self._condition:
            if exclusive:
                self._waiting_writers += 1
                try:
                    while self._writer or self._readers:
                        self._condition.wait()
                    self._writer = True
                finally:
                    self._waiting_writers -= 1
            else:
                while self._writer or self._waiting_writers:
                    self._condition.wait()
                self._readers += 1

    def release(self, exclusive: bool) -> None:
        with self._condition:
            if exclusive:
                if not self._writer:
                    raise RuntimeError("exclusive checkout gate released without ownership")
                self._writer = False
            else:
                if self._readers <= 0:
                    raise RuntimeError("shared checkout gate released without ownership")
                self._readers -= 1
            self._condition.notify_all()


class StreamRecorder:
    """Publish one stream's concise terminal events and retained text log."""

    def __init__(self, identifier: str, path: Path, reset: bool = True) -> None:
        self.identifier = identifier
        self.path = path
        self._lock = threading.Lock()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if reset or not self.path.exists():
            self.path.write_text("", encoding="utf-8")

    def emit(self, message: str) -> None:
        rendered = "[stream {}] {}".format(self.identifier, message)
        with OUTPUT_LOCK:
            print(rendered, flush=True)
        with self._lock:
            with self.path.open("a", encoding="utf-8", newline="\n") as stream:
                stream.write(rendered + "\n")


def run_batch(
    batch: Sequence[Tuple[int, Mapping[str, Any]]],
    selected_count: int,
    config: Mapping[str, Any],
    output: Path,
    dry_run: bool,
    jobs: int,
    fail_fast: bool,
    gate: Optional[CheckoutGate] = None,
    recorder: Optional[StreamRecorder] = None,
) -> List[Tuple[int, Dict[str, Any]]]:
    """Run one barrier-delimited batch and return results in manifest order."""

    def emit(message: str) -> None:
        if recorder is None:
            with OUTPUT_LOCK:
                print(message, flush=True)
        else:
            recorder.emit(message)

    def execute(test: Mapping[str, Any], case_prefix: str = "    ",
                start_message: str = "") -> Dict[str, Any]:
        exclusive = bool(test.get("exclusive", False))
        if gate is not None:
            gate.acquire(exclusive)
        try:
            # Announce RUN only after the cross-stream checkout gate is held.
            # Otherwise several queued native builds would look simultaneous
            # even though the safety gate was correctly serializing them.
            if start_message:
                emit(start_message)
            return run_one(test, config, output, dry_run, case_prefix)
        finally:
            if gate is not None:
                gate.release(exclusive)

    workers = min(jobs, len(batch))
    if workers <= 1:
        completed = []
        for index, test in batch:
            number = index + 1
            start_message = "[{}/{}] RUN  {} ({})".format(
                number, selected_count, test["id"], test["component"])
            case_prefix = ("    [stream {}] ".format(recorder.identifier)
                           if recorder is not None else "    ")
            result = execute(test, case_prefix, start_message)
            emit("[{}/{}] {:5} {} {:.3f}s {}".format(
                number, selected_count, result["status"], result["id"],
                result["seconds"], result["message"]))
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
            start_message = "[{}/{}] RUN  {} ({})".format(
                number, selected_count, test["id"], test["component"])
            # Prefix streamed case starts with the owning suite because output
            # from independent adapters may legitimately interleave.
            case_prefix = "    [{}{}] ".format(
                ("stream {} / ".format(recorder.identifier)
                 if recorder is not None else ""), test["id"])
            future = executor.submit(
                execute, test, case_prefix, start_message)
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
                emit("[{}/{}] {:5} {} {:.3f}s {}".format(
                    number, selected_count, result["status"], result["id"],
                    result["seconds"], result["message"]))
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


def run_stream(
    stream: Mapping[str, Any],
    selected_count: int,
    config: Mapping[str, Any],
    output: Path,
    dry_run: bool,
    jobs: int,
    fail_fast: bool,
    gate: CheckoutGate,
    stop_event: threading.Event,
) -> Tuple[str, List[Tuple[int, Dict[str, Any]]]]:
    """Execute one stream sequentially across its internal barrier batches."""

    identifier = str(stream["id"])
    recorder = StreamRecorder(identifier, output / "streams" / (identifier + ".log"))
    indexed_tests = list(stream["tests"])
    case_count = sum(len(item.get("selected_cases", item["expected_cases"]))
                     for _, item in indexed_tests)
    recorder.emit("START kind={} suites={} cases={} suite_jobs={}".format(
        stream["kind"], len(indexed_tests), case_count, jobs))
    completed: List[Tuple[int, Dict[str, Any]]] = []
    for batch in indexed_execution_batches(indexed_tests):
        if fail_fast and stop_event.is_set():
            break
        batch_results = run_batch(
            batch, selected_count, config, output, dry_run, jobs, fail_fast,
            gate=gate, recorder=recorder)
        completed.extend(batch_results)
        if fail_fast and any(
                result["status"] in ("FAIL", "ERROR")
                for _, result in batch_results):
            stop_event.set()
            break
    stream_counts = counts(result for _, result in completed)
    recorder.emit(
        "END executed_suites={} pass={} fail={} skip={} error={}".format(
            len(completed), stream_counts["PASS"], stream_counts["FAIL"],
            stream_counts["SKIP"], stream_counts["ERROR"]))
    return identifier, completed


def run_test_streams(
    streams: Sequence[Mapping[str, Any]],
    selected_count: int,
    config: Mapping[str, Any],
    output: Path,
    dry_run: bool,
    jobs: int,
    stream_jobs: int,
    fail_fast: bool,
) -> Tuple[List[Tuple[int, Dict[str, Any]]], Dict[str, List[Tuple[int, Dict[str, Any]]]]]:
    """Run independent streams concurrently under one global checkout gate."""

    if not streams:
        return [], {}
    gate = CheckoutGate()
    stop_event = threading.Event()
    workers = min(stream_jobs, len(streams))
    by_stream: Dict[str, List[Tuple[int, Dict[str, Any]]]] = {}
    streams_by_id = {str(stream["id"]): stream for stream in streams}

    def record_stream_failure(identifier: str, error: BaseException) -> None:
        """Convert a stream-control defect into visible combined evidence."""

        recorder = StreamRecorder(
            identifier, output / "streams" / (identifier + ".log"),
            reset=False)
        recorder.emit("ERROR stream scheduler failed: {}".format(error))
        first_index, first_test = streams_by_id[identifier]["tests"][0]
        by_stream[identifier] = [(first_index, {
            "id": first_test["id"],
            "component": first_test["component"],
            "level": first_test["level"],
            "status": "ERROR",
            "message": one_line("stream scheduler failed: {}".format(error)),
            "seconds": 0.0,
            "runner": first_test["runner"],
            "log": str(recorder.path),
            "work_dir": "", "artifact_dir": "", "cases": [],
            "evidence": [], "metrics": {},
        })]
        stop_event.set()

    if workers <= 1:
        for stream in streams:
            if fail_fast and stop_event.is_set():
                break
            identifier = str(stream["id"])
            try:
                returned_identifier, completed = run_stream(
                    stream, selected_count, config, output, dry_run, jobs,
                    fail_fast, gate, stop_event)
                if returned_identifier != identifier:
                    raise RuntimeError(
                        "stream identity mismatch: expected {}, got {}".format(
                            identifier, returned_identifier))
                by_stream[identifier] = completed
            except BaseException as error:
                record_stream_failure(identifier, error)
                if isinstance(error, (KeyboardInterrupt, SystemExit)):
                    raise
    else:
        with ThreadPoolExecutor(max_workers=workers,
                                thread_name_prefix="amps-test-stream") as executor:
            futures = {
                executor.submit(
                    run_stream, stream, selected_count, config, output,
                    dry_run, jobs, fail_fast, gate, stop_event): str(stream["id"])
                for stream in streams
            }
            for future in futures:
                identifier = futures[future]
                try:
                    returned_identifier, completed = future.result()
                    if returned_identifier != identifier:
                        raise RuntimeError(
                            "stream identity mismatch: expected {}, got {}".format(
                                identifier, returned_identifier))
                    by_stream[identifier] = completed
                except Exception as error:
                    # A stream-level scheduler defect has no trustworthy suite
                    # identity to attach. Preserve it in the stream log and
                    # let the missing selected suites remain visible in the
                    # combined report's not_executed_suites inventory.
                    record_stream_failure(identifier, error)
    combined = [item for stream in streams
                for item in by_stream.get(str(stream["id"]), [])]
    return sorted(combined, key=lambda item: item[0]), by_stream


def write_stream_reports(
    output: Path,
    streams: Sequence[Mapping[str, Any]],
    by_stream: Mapping[str, Sequence[Tuple[int, Mapping[str, Any]]]],
) -> List[Dict[str, Any]]:
    """Write one JSON receipt per stream and return combined-report metadata."""

    metadata: List[Dict[str, Any]] = []
    for stream in streams:
        identifier = str(stream["id"])
        indexed_tests = list(stream["tests"])
        completed = list(by_stream.get(identifier, []))
        results = [result for _, result in sorted(completed, key=lambda item: item[0])]
        selected_ids = [str(test["id"]) for _, test in indexed_tests]
        executed_ids = {str(result["id"]) for result in results}
        report_path = output / "streams" / (identifier + ".json")
        log_path = output / "streams" / (identifier + ".log")
        payload = {
            "schema_version": 1,
            "id": identifier,
            "kind": stream["kind"],
            "component": stream["component"],
            "selection": selected_ids,
            "selected_case_count": sum(
                len(test.get("selected_cases", test["expected_cases"]))
                for _, test in indexed_tests),
            "not_executed_suites": [
                identifier for identifier in selected_ids
                if identifier not in executed_ids],
            "counts": counts(results),
            "results": results,
            "log": str(log_path),
        }
        write_json(report_path, payload)
        metadata.append({
            "id": identifier,
            "kind": stream["kind"],
            "component": stream["component"],
            "selected_suite_count": len(selected_ids),
            "selected_case_count": payload["selected_case_count"],
            "counts": payload["counts"],
            "report": str(report_path),
            "log": str(log_path),
        })
    return metadata


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
    if args.streams <= 0:
        print("configuration error: --streams must be a positive integer", file=sys.stderr)
        return 2
    try:
        config = load_configuration(
            REPO_ROOT, TEST_ROOT / "config" / "default.json", args.config,
            output_root=args.output_dir, build_jobs=args.build_jobs,
            mpi_ranks=args.mpi_ranks)
        # A bare --list is a registry inventory, so it must expose every
        # registered application as well as core. Explicit selectors still
        # narrow listing exactly as they narrow execution.
        unfiltered_list = bool(
            args.list and not args.full_validation and not args.all and
            args.campaign is None and args.tier is None and not args.test and
            not args.match and not args.component and not args.case)
        if unfiltered_list:
            args.all = True
        elif (args.campaign is None and args.tier is None and
              not args.full_validation):
            args.campaign = config["execution"]["default_campaign"]
        # This flag is an execution request, not a machine-property setting.
        # Adapters use it to select expensive/qualification modes while still
        # validating real paths, tools, data packages, and physical gates.
        config["execution"]["full_validation"] = bool(args.full_validation)
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
        print_tests(selected, registered_tests=tests,
                    stream_jobs=args.streams, suite_jobs=args.jobs)
        return 0
    if not selected:
        print("no tests selected", file=sys.stderr)
        return 2

    output = make_output_directory(config, args.output_dir)
    write_json(output / "effective-config.json", config)
    campaign = ("full-validation" if args.full_validation else
                (args.campaign or args.tier or "smoke"))
    selected_case_ids = [
        case_id for item in selected for case_id in item["selected_cases"]]
    streams = group_test_streams(selected)
    stream_selection = [{
        "id": stream["id"],
        "kind": stream["kind"],
        "component": stream["component"],
        "suites": [test["id"] for _, test in stream["tests"]],
        "case_count": sum(
            len(test.get("selected_cases", test["expected_cases"]))
            for _, test in stream["tests"]),
    } for stream in streams]
    write_json(output / "selection.json", {
        "schema_version": 1,
        "campaign": campaign,
        "full_validation": bool(args.full_validation),
        "included_campaigns": list(CAMPAIGNS) if args.full_validation else [campaign],
        "order": args.order,
        "shuffle_seed": args.shuffle_seed,
        "stream_jobs": args.streams,
        "suite_jobs_per_stream": args.jobs,
        "streams": stream_selection,
        "case_ids": selected_case_ids,
        "tests": selected,
    })
    print("AMPS tests: campaign={} selected_streams={} stream_jobs={} suite_jobs_per_stream={} selected_suites={} selected_cases={} order={} output={}".format(
        campaign, len(streams), min(args.streams, len(streams)), args.jobs,
        len(selected), len(selected_case_ids), args.order, output))
    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    started = time.monotonic()
    indexed_results, results_by_stream = run_test_streams(
        streams, len(selected), config, output, args.dry_run, args.jobs,
        args.streams, args.fail_fast)
    results = [result for _, result in sorted(indexed_results, key=lambda item: item[0])]
    stream_reports = write_stream_reports(output, streams, results_by_stream)

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
        "full_validation": bool(args.full_validation),
        "included_campaigns": list(CAMPAIGNS) if args.full_validation else [campaign],
        "order": args.order,
        "shuffle_seed": args.shuffle_seed,
        "stream_jobs": args.streams,
        "suite_jobs_per_stream": args.jobs,
        "streams": stream_reports,
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
    print("STREAM REPORTS")
    for item in stream_reports:
        print("  {} log={} json={}".format(
            item["id"], item["log"], item["report"]))
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

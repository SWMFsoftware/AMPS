#!/usr/bin/env python3
"""Self-test the portable AMPS core-test control plane."""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import os
import signal
import sys
import threading
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, List


HERE = Path(__file__).resolve()
TEST_ROOT = HERE.parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.configuration import ConfigurationError, load_configuration  # noqa: E402
from lib.manifest import ManifestError, audit_inventory, load_manifest  # noqa: E402
from lib.runner_support import case, load_context, main_guard, write_result  # noqa: E402


def _global_module():
    """Import the executable runner as a module without starting its CLI."""

    spec = importlib.util.spec_from_file_location("amps_core_global_runner", TEST_ROOT / "run_tests.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import global runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_layered_configuration(context: Dict) -> None:
    """Check file < environment < explicit-CLI precedence and path anchoring."""

    repo = Path(context["repo_root"])
    with TemporaryDirectory(prefix="amps-test-config-") as directory:
        overlay = Path(directory) / "machine.json"
        overlay.write_text(json.dumps({
            "schema_version": 1,
            "paths": {"output_root": "relative-results"},
            "tools": {"cxx": ["machine-cxx", "--driver"]},
            "execution": {"build_jobs": 3},
        }), encoding="utf-8")
        resolved = load_configuration(
            repo, TEST_ROOT / "config" / "default.json", [overlay],
            environ={"AMPS_TEST_BUILD_JOBS": "7", "AMPS_TEST_MPI_RANKS": "2"},
            build_jobs=11)
        assert resolved["tools"]["cxx"] == ["machine-cxx", "--driver"]
        assert resolved["execution"]["build_jobs"] == 11
        assert resolved["execution"]["mpi_ranks"] == 2
        assert resolved["paths"]["output_root"] == str((repo / "relative-results").resolve())


def test_schema_extensions(context: Dict) -> None:
    """Exercise backward-compatible capability and schema-version handling."""

    repo = Path(context["repo_root"])
    with TemporaryDirectory(prefix="amps-test-schema-") as directory:
        root = Path(directory)
        capability_file = root / "capabilities.json"
        capability_file.write_text(json.dumps({
            "schema_version": 1,
            "capabilities": {
                "cuda": {"available": False, "policy": "required"},
                "spice": True,
            },
        }), encoding="utf-8")
        resolved = load_configuration(
            repo, TEST_ROOT / "config" / "default.json", [capability_file],
            environ={})
        assert resolved["capabilities"]["cuda"] == {
            "available": False, "policy": "required"}
        assert resolved["capabilities"]["spice"] == {
            "available": True, "policy": "optional"}
        assert resolved["execution"]["default_campaign"] == "smoke"
        profiled = load_configuration(
            repo, TEST_ROOT / "config" / "default.json", [],
            environ={"AMPS_TEST_CAMPAIGN": "profile"})
        assert profiled["execution"]["default_campaign"] == "profile"
        module = _global_module()
        capability_test = {
            "id": "CAPABILITY-FIXTURE", "component": "src/test",
            "level": "unit", "runner": "never-started.py",
            "requires": {"tools": [], "paths": [], "capabilities": ["cuda"]}}
        reasons = module.missing_requirements(capability_test, resolved)
        assert reasons == ["required capability 'cuda' is unavailable"]
        assert module.required_capability_missing(reasons)
        policy_result = module.run_one(
            capability_test, resolved, root / "capability-output", False)
        assert policy_result["status"] == "ERROR"

        bad_version = root / "bad-version.json"
        bad_version.write_text('{"schema_version": 2}', encoding="utf-8")
        try:
            load_configuration(repo, TEST_ROOT / "config" / "default.json",
                               [bad_version], environ={})
        except ConfigurationError:
            pass
        else:
            raise AssertionError("configuration schema_version 2 was accepted")

        manifest_payload = json.loads((TEST_ROOT / "manifest.json").read_text(
            encoding="utf-8"))
        manifest_payload["tests"][1]["expected_cases"][0] = \
            manifest_payload["tests"][0]["expected_cases"][0]
        bad_manifest = root / "duplicate-case.json"
        bad_manifest.write_text(json.dumps(manifest_payload), encoding="utf-8")
        try:
            load_manifest(bad_manifest, repo, TEST_ROOT)
        except ManifestError:
            pass
        else:
            raise AssertionError("globally duplicated case ID was accepted")

        manifest_payload = json.loads((TEST_ROOT / "manifest.json").read_text(
            encoding="utf-8"))
        manifest_payload["schema_version"] = 2
        bad_manifest_version = root / "bad-manifest-version.json"
        bad_manifest_version.write_text(
            json.dumps(manifest_payload), encoding="utf-8")
        try:
            load_manifest(bad_manifest_version, repo, TEST_ROOT)
        except ManifestError:
            pass
        else:
            raise AssertionError("manifest schema_version 2 was accepted")


def test_registry_and_inventory(context: Dict) -> None:
    """Ensure the committed registry and coverage ledger agree with this tree."""

    repo = Path(context["repo_root"])
    tests = load_manifest(TEST_ROOT / "manifest.json", repo, TEST_ROOT)
    entries, errors = audit_inventory(TEST_ROOT / "component_inventory.json", repo, tests)
    assert not errors, "; ".join(errors)
    identifiers = {item["id"] for item in tests}
    assert "CORE-HARNESS-001" in identifiers
    assert len(entries) >= 20


def test_selection_contract(context: Dict) -> None:
    """Pin exact-ID selection so future scheduler changes cannot broaden runs."""

    module = _global_module()
    repo = Path(context["repo_root"])
    tests = load_manifest(TEST_ROOT / "manifest.json", repo, TEST_ROOT)
    args = argparse.Namespace(
        test=["CORE-GENERAL-INPUT-001"], match=[], component=[], all=False,
        case=[], campaign="smoke", tier=None)
    selected = module.select_tests(tests, args)
    assert [item["id"] for item in selected] == ["CORE-GENERAL-INPUT-001"]
    case_args = argparse.Namespace(
        test=[], match=[], component=[], all=False, case=["INP-03"],
        campaign="smoke", tier=None)
    selected = module.select_tests(tests, case_args)
    assert [item["id"] for item in selected] == ["CORE-GENERAL-INPUT-001"]
    assert selected[0]["selected_cases"] == ["INP-03"]
    alias_args = argparse.Namespace(
        test=["CORE-PIC-PHYSICS-001"], match=[], component=[], all=False,
        case=[], campaign="smoke", tier=None)
    selected = module.select_tests(tests, alias_args)
    assert [item["id"] for item in selected] == ["CORE-PIC-NATIVE-001"]
    assert len(selected[0]["selected_cases"]) == 19
    assert all(value.startswith("PHYS-") for value in selected[0]["selected_cases"])
    legacy_tier_args = argparse.Namespace(
        test=[], match=[], component=[], all=False, case=[],
        campaign=None, tier="smoke")
    assert module.select_tests(tests, legacy_tier_args)

    # Campaign names are operational contracts, not aliases for identical
    # selections. Smoke is dependency-light, nightly is the complete current
    # deterministic registry, distributed selects the two-rank channel suite,
    # extended adds the optional-provider/interface module to the configured
    # native applications, and profile selects only the clean native profile
    # matrix.  T3's mesh owner performs three serial DIM builds internally.
    # Pinning these complete sets makes a campaign-membership edit deliberate:
    # a new suite cannot silently introduce an MPI launch or AMPS rebuild into
    # machine bring-up, nor disappear from the routine nightly campaign.
    def campaign_ids(name: str) -> set:
        campaign_args = argparse.Namespace(
            test=[], match=[], component=[], all=False, case=[],
            campaign=name, tier=None)
        return {item["id"] for item in module.select_tests(tests, campaign_args)}

    assert campaign_ids("smoke") == {
        "CORE-HARNESS-001", "CORE-GENERAL-INPUT-001",
        "CORE-GENERAL-CONTAINERS-001", "CORE-GENERAL-NUMERICS-001",
        "CORE-GENERAL-RNG-001", "CORE-GENERAL-LOGGER-001",
        "CORE-UNITS-001", "CORE-SPECIES-RADIATION-001"}
    assert campaign_ids("nightly") == {
        "CORE-HARNESS-001", "CORE-GENERAL-INPUT-001",
        "CORE-GENERAL-CONTAINERS-001", "CORE-GENERAL-NUMERICS-001",
        "CORE-GENERAL-RNG-001", "CORE-GENERAL-MPI-001",
        "CORE-GENERAL-LOGGER-001", "CORE-UNITS-001",
        "CORE-SPECIES-RADIATION-001", "CORE-INTERFACE-ABI-001",
        "CORE-PIC-NATIVE-001", "CORE-PIC-WEIGHT-TIME-LOCAL-001",
        "CORE-PIC-WEIGHT-TIME-SINGLE-GLOBAL-001",
        "CORE-PIC-SPLIT-MERGE-001",
        "CORE-PIC-PBUFFER-LAYOUT-ALIGNED-FIXED-001",
        "CORE-PIC-PBUFFER-LAYOUT-ALIGNED-INDIVIDUAL-001",
        "CORE-MESH-NATIVE-001"}
    assert campaign_ids("distributed") == {"CORE-GENERAL-MPI-001"}
    assert campaign_ids("extended") == {
        "CORE-INTERFACE-ABI-001", "CORE-PIC-NATIVE-001",
        "CORE-PIC-SPLIT-MERGE-001", "CORE-MESH-NATIVE-001"}
    assert campaign_ids("profile") == {
        "CORE-PIC-NATIVE-001", "CORE-PIC-WEIGHT-TIME-LOCAL-001",
        "CORE-PIC-WEIGHT-TIME-SINGLE-GLOBAL-001",
        "CORE-PIC-SPLIT-MERGE-001",
        "CORE-PIC-PBUFFER-LAYOUT-ALIGNED-FIXED-001",
        "CORE-PIC-PBUFFER-LAYOUT-ALIGNED-INDIVIDUAL-001",
        "CORE-MESH-NATIVE-001"}

    unknown_case_args = argparse.Namespace(
        test=[], match=[], component=[], all=False, case=["UNKNOWN-CASE-001"],
        campaign="smoke", tier=None)
    try:
        module.select_tests(tests, unknown_case_args)
    except ManifestError:
        pass
    else:
        raise AssertionError("unknown exact case ID was accepted")

    ordered = module.order_tests(tests[:3], "shuffle", 1729)
    repeated = module.order_tests(tests[:3], "shuffle", 1729)
    assert [(item["id"], item["selected_cases"]) for item in ordered] == [
        (item["id"], item["selected_cases"]) for item in repeated]
    reversed_tests = module.order_tests(tests[:2], "reverse", 0)
    assert [item["id"] for item in reversed_tests] == [
        tests[1]["id"], tests[0]["id"]]
    assert reversed_tests[-1]["selected_cases"] == list(
        reversed(tests[0]["expected_cases"]))


def test_result_protocol(context: Dict) -> None:
    """Prove that a PASS wrapper cannot conceal a failed individual case."""

    module = _global_module()
    test = {
        "id": "PROTOCOL-TEST", "component": "src/test", "level": "unit",
        "expected_cases": ["CASE-A"],
    }
    payload = {
        "schema_version": 1, "id": "PROTOCOL-TEST", "status": "PASS",
        "message": "ok", "cases": [{"id": "CASE-A", "status": "PASS", "message": "ok"}],
    }
    status, message, cases = module.validate_result(test, payload, 0)
    assert status == "PASS" and len(cases) == 1
    hidden = dict(payload)
    hidden["cases"] = [{"id": "CASE-A", "status": "FAIL", "message": "bad"}]
    status, message, cases = module.validate_result(test, hidden, 0)
    assert status == "ERROR" and "hides" in message
    empty = dict(payload)
    empty["cases"] = []
    status, message, cases = module.validate_result(test, empty, 0)
    assert status == "ERROR" and "without any executed cases" in message
    wrong = dict(payload)
    wrong["cases"] = [{"id": "CASE-B", "status": "PASS", "message": "ok"}]
    status, message, cases = module.validate_result(test, wrong, 0)
    assert status == "ERROR" and "registry mismatch" in message


def test_scheduler_barriers(context: Dict) -> None:
    """Prove exclusive suites are isolated without reordering the registry."""

    module = _global_module()
    tests = [
        {"id": "A", "exclusive": False},
        {"id": "B", "exclusive": False},
        {"id": "C", "exclusive": True},
        {"id": "D", "exclusive": False},
        {"id": "E", "exclusive": True},
    ]
    batches = module.execution_batches(tests)
    observed = [[item[1]["id"] for item in batch] for batch in batches]
    assert observed == [["A", "B"], ["C"], ["D"], ["E"]]


def _fixture_config(context: Dict, directory: Path) -> Dict:
    """Return a real normalized profile redirected into a disposable root."""

    config = load_configuration(
        Path(context["repo_root"]), TEST_ROOT / "config" / "default.json",
        environ={})
    config["paths"]["output_root"] = str(directory)
    config["paths"]["temporary_root"] = ""
    config["tools"]["python"] = [sys.executable]
    config["execution"]["keep_work"] = "always"
    return config


def _fixture_test(identifier: str, runner: Path, mode: str,
                  timeout: int = 5) -> Dict:
    return {
        "id": identifier,
        "component": "src/test",
        "level": "unit",
        "tiers": ["smoke"],
        "campaigns": ["smoke"],
        "runner": "fixture.py",
        "runner_path": str(runner),
        "description": "process-control fixture",
        "requires": {"tools": ["python"], "paths": [], "capabilities": []},
        "expected_cases": ["FIXTURE-CASE-001"],
        "selected_cases": ["FIXTURE-CASE-001"],
        "fixture_mode": mode,
        "timeout_seconds": timeout,
        "exclusive": False,
    }


def _write_fixture_runner(path: Path) -> None:
    """Create a hostile child adapter used to exercise the real subprocess path."""

    path.write_text(r'''#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--context", required=True)
args = parser.parse_args()
context = json.loads(Path(args.context).read_text(encoding="utf-8"))
mode = context["test"].get("fixture_mode")
result = Path(context["result_file"])
if mode == "malformed":
    result.write_text("{not-json", encoding="utf-8")
elif mode == "partial":
    Path(str(result) + ".tmp").write_text('{"schema_version": 1', encoding="utf-8")
elif mode == "missing" or mode == "crash":
    if mode == "crash":
        os._exit(7)
elif mode == "timeout":
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    artifact = Path(context["artifact_dir"])
    artifact.mkdir(parents=True, exist_ok=True)
    (artifact / "child.pid").write_text(str(child.pid), encoding="utf-8")
    time.sleep(30)
else:
    payload = {
        "schema_version": 1,
        "id": context["test"]["id"],
        "status": "PASS",
        "message": "fixture pass",
        "cases": [{
            "id": "FIXTURE-CASE-001",
            "status": "FAIL" if mode == "hidden-failure" else "PASS",
            "message": "injected hidden failure" if mode == "hidden-failure" else "ok",
        }],
    }
    result.write_text(json.dumps(payload), encoding="utf-8")
    if mode == "bad-exit":
        sys.exit(1)
''', encoding="utf-8")


def test_subprocess_receipts(context: Dict) -> None:
    """Reject stale, missing, partial, malformed, crashed, and exit-mismatched receipts."""

    module = _global_module()
    with TemporaryDirectory(prefix="amps-test-subprocess-") as directory:
        root = Path(directory)
        runner = root / "fixture_runner.py"
        _write_fixture_runner(runner)
        config = _fixture_config(context, root)
        for mode in ("malformed", "partial", "missing", "crash", "bad-exit",
                     "hidden-failure"):
            output = root / ("output-" + mode)
            test = _fixture_test("FIXTURE-" + mode.upper(), runner, mode)
            if mode == "missing":
                stale = output / "work" / test["id"] / "result.json"
                stale.parent.mkdir(parents=True, exist_ok=True)
                stale.write_text(json.dumps({
                    "schema_version": 1, "id": test["id"], "status": "PASS",
                    "message": "stale", "cases": [{
                        "id": "FIXTURE-CASE-001", "status": "PASS", "message": "stale"}],
                }), encoding="utf-8")
            result = module.run_one(test, config, output, False)
            assert result["status"] == "ERROR", (mode, result)


def test_timeout_cleanup(context: Dict) -> None:
    """Timeout and an actual SIGINT must both reap the complete child tree."""

    module = _global_module()
    with TemporaryDirectory(prefix="amps-test-timeout-") as directory:
        root = Path(directory)
        runner = root / "fixture_runner.py"
        _write_fixture_runner(runner)
        config = _fixture_config(context, root)
        output = root / "output"
        test = _fixture_test("FIXTURE-TIMEOUT", runner, "timeout", timeout=1)
        result = module.run_one(test, config, output, False)
        assert result["status"] == "ERROR" and "timed out" in result["message"]
        pid_file = output / "artifacts" / test["id"] / "child.pid"
        child_pid = int(pid_file.read_text(encoding="utf-8"))
        def descendant_is_gone(child_pid: int) -> bool:
            """Treat a reaped process or a transient zombie as no longer active."""

            proc_stat = Path("/proc") / str(child_pid) / "stat"
            if not proc_stat.exists():
                return True
            fields = proc_stat.read_text(encoding="utf-8", errors="replace").split()
            if len(fields) > 2 and fields[2] == "Z":
                return True
            return False

        deadline = time.monotonic() + 3.0
        while time.monotonic() < deadline and not descendant_is_gone(child_pid):
            time.sleep(0.05)
        if not descendant_is_gone(child_pid):
            raise AssertionError("timed-out descendant {} remains alive".format(child_pid))

        # Exercise run_one's BaseException boundary with a real interrupt, not
        # a monkey-patched exception.  The forked controller is deliberately
        # outside the adapter's new session.  Sending SIGINT only to that
        # controller therefore proves that run_one catches KeyboardInterrupt,
        # kills the adapter-owned process group, waits for it, and then
        # re-raises.  This is the same control path used by Ctrl-C in a nightly
        # driver and is distinct from the timeout path above.
        if os.name != "posix" or not hasattr(os, "fork"):
            raise AssertionError("interrupt cleanup test requires POSIX process groups")
        interrupt_output = root / "interrupt-output"
        interrupt_test = _fixture_test(
            "FIXTURE-INTERRUPT", runner, "timeout", timeout=20)
        controller_pid = os.fork()
        if controller_pid == 0:  # pragma: no cover - assertions run in parent
            signal.signal(signal.SIGINT, signal.default_int_handler)
            try:
                module.run_one(
                    interrupt_test, config, interrupt_output, False)
            except KeyboardInterrupt:
                os._exit(130)
            except BaseException:
                os._exit(2)
            os._exit(3)

        interrupt_child_file = (
            interrupt_output / "artifacts" / interrupt_test["id"] / "child.pid")
        waited = False
        try:
            deadline = time.monotonic() + 5.0
            while time.monotonic() < deadline and not interrupt_child_file.is_file():
                time.sleep(0.05)
            assert interrupt_child_file.is_file(), "interrupt fixture did not start descendant"
            interrupt_child_pid = int(interrupt_child_file.read_text(encoding="utf-8"))
            os.kill(controller_pid, signal.SIGINT)
            _, controller_status = os.waitpid(controller_pid, 0)
            waited = True
            assert os.WIFEXITED(controller_status), controller_status
            assert os.WEXITSTATUS(controller_status) == 130, controller_status
            deadline = time.monotonic() + 3.0
            while (time.monotonic() < deadline and
                   not descendant_is_gone(interrupt_child_pid)):
                time.sleep(0.05)
            if not descendant_is_gone(interrupt_child_pid):
                raise AssertionError(
                    "interrupted descendant {} remains alive".format(
                        interrupt_child_pid))
        finally:
            if not waited:
                try:
                    os.kill(controller_pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                os.waitpid(controller_pid, 0)


def test_scheduler_concurrency(context: Dict) -> None:
    """Bound concurrency and stop scheduling after fail-fast observes failure."""

    module = _global_module()
    lock = threading.Lock()
    active = [0]
    maximum = [0]
    started: List[str] = []

    def fake_run(test, config, output, dry_run, case_prefix="    "):
        with lock:
            active[0] += 1
            maximum[0] = max(maximum[0], active[0])
            started.append(test["id"])
        time.sleep(0.01 if test.get("fail") else 0.04)
        with lock:
            active[0] -= 1
        status = "FAIL" if test.get("fail") else "PASS"
        return {
            "id": test["id"], "component": "src/test", "level": "unit",
            "status": status, "message": "fixture", "seconds": 0.04,
            "runner": "fixture", "log": "", "work_dir": "",
            "artifact_dir": "", "cases": [], "evidence": [], "metrics": {},
        }

    original = module.run_one
    module.run_one = fake_run
    try:
        tests = [(index, {"id": name, "component": "src/test", "level": "unit",
                          "runner": "fixture", "fail": False})
                 for index, name in enumerate(("A", "B", "C", "D"))]
        completed = module.run_batch(
            tests, 4, {}, Path("."), False, jobs=2, fail_fast=False)
        assert [index for index, _ in completed] == [0, 1, 2, 3]
        assert maximum[0] == 2

        started[:] = []
        active[0] = 0
        maximum[0] = 0
        fail_tests = [
            (0, {"id": "FAIL-FIRST", "component": "src/test", "level": "unit",
                 "runner": "fixture", "fail": True}),
            (1, {"id": "IN-FLIGHT", "component": "src/test", "level": "unit",
                 "runner": "fixture", "fail": False}),
            (2, {"id": "MUST-NOT-START", "component": "src/test", "level": "unit",
                 "runner": "fixture", "fail": False}),
        ]
        module.run_batch(
            fail_tests, 3, {}, Path("."), False, jobs=2, fail_fast=True)
        assert set(started) == {"FAIL-FIRST", "IN-FLIGHT"}
    finally:
        module.run_one = original


def test_final_result_inventory(context: Dict) -> None:
    """Require terminal, JSON, and JUnit views to reconcile exactly."""

    module = _global_module()
    results = [
        {
            "id": "SUITE-A", "component": "src/test", "status": "PASS",
            "message": "ok", "seconds": 0.1, "log": "suite-a.log", "cases": [
                {"id": "CASE-A1", "status": "PASS", "message": "ok"},
                {"id": "CASE-A2", "status": "FAIL", "message": "bad"},
            ],
        },
        {
            "id": "SUITE-B", "component": "src/test", "status": "SKIP",
            "message": "missing tool", "seconds": 0.025,
            "log": "suite-b.log", "cases": [],
        },
    ]
    rows = module.individual_results(results)
    assert [row["id"] for row in rows] == ["CASE-A1", "CASE-A2", "SUITE-B"]
    assert [row["kind"] for row in rows] == ["case", "case", "suite"]
    assert module.counts(rows) == {
        "PASS": 1, "FAIL": 1, "SKIP": 1, "ERROR": 0, "TOTAL": 3}

    # The user's contract is the terminal report, not merely an internal list.
    # Capture the formatter directly so removal of the heading, a case row, or
    # the typed suite-level fallback becomes a deterministic harness failure.
    terminal = io.StringIO()
    with contextlib.redirect_stdout(terminal):
        module.print_final_results(rows)
    report = terminal.getvalue()
    assert report.startswith("FINAL INDIVIDUAL TEST RESULTS\n")
    assert "PASS  CASE-A1" in report
    assert "FAIL  CASE-A2" in report
    assert "SKIP  SUITE-B" in report
    assert "suite-level: no individual case evidence" in report

    # Exercise the actual report writers, not a second implementation of their
    # counting rules.  The JSON normalized rows and JUnit attributes/case names
    # must describe the same three outcomes printed above.
    with TemporaryDirectory(prefix="amps-test-report-") as directory:
        report_root = Path(directory)
        json_path = report_root / "result.json"
        junit_path = report_root / "result.xml"
        module.write_json(json_path, {
            "individual_results": rows,
            "individual_counts": module.counts(rows),
        })
        module.write_junit(junit_path, results, 0.125)
        json_receipt = json.loads(json_path.read_text(encoding="utf-8"))
        assert json_receipt["individual_results"] == rows
        assert json_receipt["individual_counts"] == {
            "PASS": 1, "FAIL": 1, "SKIP": 1, "ERROR": 0, "TOTAL": 3}
        junit = ET.parse(str(junit_path)).getroot()
        assert junit.attrib["tests"] == "3"
        assert junit.attrib["failures"] == "1"
        assert junit.attrib["errors"] == "0"
        assert junit.attrib["skipped"] == "1"
        assert [node.attrib["name"] for node in junit.findall("testcase")] == [
            "SUITE-A/CASE-A1", "SUITE-A/CASE-A2", "SUITE-B/SUITE-B"]


def test_traceability_contract(context: Dict) -> None:
    """Require exhaustive source/mode classification and known case owners."""

    repo = Path(context["repo_root"])
    tests = load_manifest(TEST_ROOT / "manifest.json", repo, TEST_ROOT)
    case_ids = [case_id for test in tests for case_id in test["expected_cases"]]
    assert case_ids and len(case_ids) == len(set(case_ids))
    for test in tests:
        assert test["expected_cases"], test["id"]

    spec = importlib.util.spec_from_file_location(
        "amps_traceability_validator", TEST_ROOT / "validate_traceability.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import traceability validator")
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    validation = validator.Validation()
    trace = validator.load_json(TEST_ROOT / "source_traceability.json", validation)
    contracts = validator.load_json(
        TEST_ROOT / "configuration_contracts.json", validation)
    manifest = validator.load_json(TEST_ROOT / "manifest.json", validation)
    inventory = validator.load_json(
        TEST_ROOT / "component_inventory.json", validation)
    expanded = validator.validate_traceability(
        trace, manifest, inventory, validation)
    contract_counts = validator.validate_contracts(
        contracts, manifest, trace, validation)
    assert not validation.errors, "; ".join(validation.errors)
    assert expanded["translation_unit_count"] > 0
    assert expanded["manifest_case_count"] == len(case_ids)
    assert sum(contract_counts.values()) > 0


def main() -> int:
    context = load_context("Self-test the AMPS core-test runner")
    tests: List[tuple] = [
        ("HARNESS-CONFIG-001", "configuration precedence and path resolution",
         test_layered_configuration),
        ("HARNESS-MANIFEST-001", "manifest and component-inventory consistency",
         test_registry_and_inventory),
        ("HARNESS-SELECTION-001", "exact test-selection contract",
         test_selection_contract),
        ("HARNESS-PROTOCOL-001", "result protocol rejects hidden failures",
         test_result_protocol),
        ("HARNESS-SCHEDULER-001", "concurrent scheduler preserves exclusive barriers",
         test_scheduler_barriers),
        ("HARNESS-REPORT-001", "final individual-result inventory and counts reconcile",
         test_final_result_inventory),
        ("HARNESS-SCHEMA-001", "campaign/capability schema and version rejection",
         test_schema_extensions),
        ("HARNESS-SUBPROCESS-001", "hostile subprocess receipt rejection",
         test_subprocess_receipts),
        ("HARNESS-TIMEOUT-001", "timeout/interrupt kill complete process groups",
         test_timeout_cleanup),
        ("HARNESS-CONCURRENCY-001", "bounded concurrency and fail-fast scheduling",
         test_scheduler_concurrency),
        ("HARNESS-TRACEABILITY-001", "global case ownership registry",
         test_traceability_contract),
    ]
    by_id = {item[0]: item for item in tests}
    selected = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", by_id)))
    results = []
    for identifier in selected:
        if identifier not in by_id:
            results.append(case(identifier, "ERROR", "unknown harness case"))
            continue
        identifier, description, function = by_id[identifier]
        print("[CASE] RUN  {} {}".format(identifier, description), flush=True)
        try:
            function(context)
            results.append(case(identifier, "PASS", "ok"))
        except Exception as error:  # every self-test failure belongs in structured evidence
            results.append(case(identifier, "FAIL", "{}: {}".format(type(error).__name__, error)))
    status = "PASS" if all(item["status"] == "PASS" for item in results) else "FAIL"
    return write_result(context, status,
                        "{} of {} harness checks passed".format(
                            sum(item["status"] == "PASS" for item in results), len(results)),
                        cases=results)


if __name__ == "__main__":
    main_guard(main)

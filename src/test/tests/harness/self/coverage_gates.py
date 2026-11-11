"""Adversarial coverage declarations and the real standalone/strict CLIs.

These are governance tests. Their synthetic complete tree cannot establish
coverage of any AMPS component, and a truthful incomplete real inventory is
expected to remain nonpassing under the strict CLI.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

from lib.manifest import audit_inventory, incomplete_components


def test_coverage_gates(context: dict) -> None:
    """Require real CLI execution and reject contradictory/fake coverage."""

    repo = Path(context["repo_root"])
    test_root = repo / "src/test"
    artifacts = Path(context["artifact_dir"])
    work = Path(context["work_dir"])
    python = list(context["tools"]["python"])
    expanded_path = artifacts / "standalone-traceability.json"
    completed = subprocess.run(
        python + ["-B", str(test_root / "validate_traceability.py"),
                  "--output", str(expanded_path)],
        cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=120, check=False)
    (artifacts / "standalone-traceability.txt").write_text(completed.stdout)
    assert completed.returncode == 0, completed.stdout
    expanded = json.loads(expanded_path.read_text())
    assert expanded["manifest_case_count"] > 0
    assert expanded["translation_unit_count"] > 0

    # The same validator is also an importable read-only authority in the
    # global strict CLI. Importing it here exercises the actual owner/gap
    # validation; no test-local substitute decides whether a declaration is
    # internally consistent.
    spec = importlib.util.spec_from_file_location(
        "amps_coverage_gate_validator", test_root / "validate_traceability.py")
    assert spec is not None and spec.loader is not None
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)

    scenarios = []
    with TemporaryDirectory(prefix="amps-coverage-gates-", dir=work) as directory:
        fixture = Path(directory)
        components = ["src/general", "src/interface", "src/meshAMR",
                      "src/pic", "src/species", "src/test"]
        for component in components:
            (fixture / component).mkdir(parents=True)
        entries = [{"component": component, "status": "covered",
                    "test_ids": ["FIXTURE-OWNER-001"]}
                   for component in components]
        manifest = [{"id": "FIXTURE-OWNER-001"}]
        inventory_file = fixture / "inventory.json"

        def audit(rows: list) -> tuple:
            inventory_file.write_text(json.dumps({"schema_version": 1,
                                                  "components": rows}))
            return audit_inventory(inventory_file, fixture, manifest)

        observed, errors = audit(entries)
        assert not errors and not incomplete_components(observed)
        scenarios.append("synthetic complete declarations accepted only as fixture")
        partial = copy.deepcopy(entries)
        partial[0]["status"] = "partial"
        observed, errors = audit(partial)
        assert not errors and incomplete_components(observed) == ["src/general"]
        scenarios.append("partial component cannot satisfy strict completeness")
        missing_owner = copy.deepcopy(entries)
        missing_owner[0]["test_ids"] = []
        assert any("covered but has no" in message
                   for message in audit(missing_owner)[1])
        scenarios.append("covered component requires registered owners")
        unknown_owner = copy.deepcopy(entries)
        unknown_owner[0]["test_ids"] = ["UNKNOWN-OWNER-001"]
        assert any("unknown tests" in message
                   for message in audit(unknown_owner)[1])
        scenarios.append("unknown owner is rejected")
        (fixture / "src/models/unregistered").mkdir(parents=True)
        assert any("missing from inventory" in message
                   for message in audit(entries)[1])
        scenarios.append("new discovered production component is rejected")

    owner = {"status": "tested", "reason": "independently evidenced fixture",
             "suite_ids": ["FIXTURE-OWNER-001"],
             "case_ids": ["FIXTURE-CASE-001"], "gaps": []}

    def validate_owner(value: dict) -> list:
        validation = validator.Validation()
        validator.validate_owner(
            value, "fixture", {"FIXTURE-OWNER-001"},
            {"FIXTURE-CASE-001"}, {"FIXTURE-CASE-001"}, validation)
        return validation.errors

    assert not validate_owner(owner)
    with_gap = copy.deepcopy(owner)
    with_gap["gaps"] = ["required physics has no independent reference"]
    assert any("cannot retain gaps" in message
               for message in validate_owner(with_gap))
    scenarios.append("tested declaration retaining a required gap is rejected")
    without_cases = copy.deepcopy(owner)
    without_cases["case_ids"] = []
    assert any("requires direct cases" in message
               for message in validate_owner(without_cases))
    scenarios.append("tested declaration without direct assertions is rejected")
    stale_case = copy.deepcopy(owner)
    stale_case["case_ids"] = ["STALE-CASE-001"]
    assert any("absent from" in message
               for message in validate_owner(stale_case))
    scenarios.append("stale case identity is rejected")

    # Use an explicit portable profile even for a read-only nested CLI. Only
    # Python is used here; its argv comes from the resolved parent profile.
    # No hostname path, MPI installation, or global output default is assumed.
    profile = work / "coverage-gates-machine.json"
    profile.write_text(json.dumps({
        "schema_version": 1,
        "tools": {"python": python},
        "paths": {"output_root": str(work / "unused-results"),
                  "temporary_root": str(work / "unused-work")},
    }, indent=2) + "\n")
    strict = subprocess.run(
        python + ["-B", str(test_root / "run_tests.py"), "--config", str(profile),
                  "--coverage", "--require-complete-coverage"],
        cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=120, check=False)
    (artifacts / "strict-coverage-cli.txt").write_text(strict.stdout)
    real_inventory = json.loads((test_root / "component_inventory.json").read_text())
    incomplete = incomplete_components(real_inventory["components"])
    assert strict.returncode == (1 if incomplete else 0), strict.stdout
    assert "coverage: components=" in strict.stdout
    scenarios.append("real strict CLI preserves incomplete-inventory exit status")

    # Exercise the strict CLI's FAILURE control flow too. Merely checking an
    # incomplete inventory's exit 1 cannot detect removal of the cross-ledger
    # audit: the inventory would still produce exit 1. In this private child
    # process only, inject one unmistakable governance defect through the real
    # validator interface, then execute the unmodified CLI. It must surface
    # that defect as configuration ERROR=2, never ordinary incomplete FAIL=1.
    # The successful standalone invocation above separately exercises the real
    # validator itself; this fixture does not replace that positive evidence
    # or manufacture complete coverage for any production component.
    fault = "INJECTED-COVERAGE-LEDGER-CONTRADICTION-001"
    bootstrap = work / "strict-coverage-fault-control.py"
    argv = [str(test_root / "run_tests.py"), "--config", str(profile),
            "--coverage", "--require-complete-coverage"]
    bootstrap.write_text(
        "import runpy, sys\n"
        f"sys.path.insert(0, {str(test_root)!r})\n"
        "import validate_traceability as validator\n"
        "def injected_audit():\n"
        "    validation = validator.Validation()\n"
        f"    validation.errors.append({fault!r})\n"
        "    return validation, {}\n"
        "validator.audit_repository = injected_audit\n"
        f"sys.argv = {argv!r}\n"
        "runpy.run_path(sys.argv[0], run_name='__main__')\n")
    rejected = subprocess.run(
        python + ["-B", str(bootstrap)], cwd=repo, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120, check=False)
    (artifacts / "strict-coverage-injected-defect.txt").write_text(rejected.stdout)
    assert rejected.returncode == 2 and fault in rejected.stdout, rejected.stdout
    scenarios.append("real strict CLI surfaces injected ledger contradiction as ERROR=2")
    (artifacts / "coverage-gate-controls.json").write_text(json.dumps({
        "scenarios": scenarios, "scenario_count": len(scenarios),
        "standalone_cli_exit": completed.returncode,
        "strict_cli_exit": strict.returncode,
        "injected_defect_cli_exit": rejected.returncode,
        "real_incomplete_components": incomplete,
        "qualification": "governance controls only; no component promoted",
    }, indent=2) + "\n")

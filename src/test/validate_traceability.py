#!/usr/bin/env python3
"""Validate and expand the AMPS T0 source/configuration traceability ledgers.

This program deliberately uses only the Python standard library.  It is a
structural governance check, not a scientific test: discovering a source file,
compiling it, or naming a case never upgrades behavioral coverage.  Status is
reviewed explicitly in the two ledgers and the validator only proves that the
declared ownership is complete, internally consistent, and anchored to files
and stable registered case IDs that exist in the current checkout.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Set, Tuple


TEST_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TEST_ROOT.parents[1]
sys.path.insert(0, str(TEST_ROOT))

from lib.manifest import ManifestError, load_manifest  # noqa: E402
from generate_test_catalog import render_catalog  # noqa: E402
ALLOWED_STATUSES = {"tested", "partial", "planned", "exempt"}
INVENTORY_STATUS = {
    "tested": "covered", "partial": "partial", "planned": "planned",
    # src/test is intentionally exempt from production-source coverage, while
    # its own remaining harness work is represented as partial in inventory.
    "exempt": "partial",
}
STABLE_ID = re.compile(r"(?<![A-Z0-9])([A-Z][A-Z0-9]*(?:-[A-Z0-9]+)+)(?![A-Z0-9])")


class Validation:
    """Accumulate all structural faults so one run gives actionable output."""

    def __init__(self) -> None:
        self.errors: List[str] = []

    def require(self, condition: bool, message: str) -> None:
        if not condition:
            self.errors.append(message)


def load_json(path: Path, validation: Validation) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        validation.errors.append("{}: {}".format(path.relative_to(REPO_ROOT), error))
        return {}
    validation.require(isinstance(value, dict), "{}: root must be an object".format(
        path.relative_to(REPO_ROOT)))
    return value if isinstance(value, dict) else {}


def repo_path(raw: str) -> Path:
    """Resolve a repository-relative authority without accepting escapes."""

    candidate = (REPO_ROOT / raw).resolve()
    try:
        candidate.relative_to(REPO_ROOT)
    except ValueError as error:
        raise ValueError("path escapes repository: {}".format(raw)) from error
    return candidate


def excluded(relative: str, patterns: Sequence[str]) -> bool:
    return any(fnmatch.fnmatchcase(relative, pattern) for pattern in patterns)


def discover_files(discovery: Mapping[str, Any]) -> Tuple[List[str], List[str]]:
    """Return production translation units and compile-time public contracts."""

    translation_extensions = set(discovery.get("translation_unit_extensions", []))
    contract_extensions = set(discovery.get("public_contract_extensions", []))
    excludes = list(discovery.get("exclude_globs", []))
    translation_units: List[str] = []
    public_contracts: List[str] = []
    for path in (REPO_ROOT / "src").rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(REPO_ROOT).as_posix()
        if excluded(relative, excludes):
            continue
        suffix = path.suffix.lower()
        if suffix in translation_extensions:
            translation_units.append(relative)
        if suffix in contract_extensions:
            public_contracts.append(relative)
    return sorted(translation_units), sorted(public_contracts)


def owning_area(relative: str, areas: Sequence[Mapping[str, Any]]) -> List[str]:
    candidates = []
    for area in areas:
        prefix = str(area.get("path", "")).rstrip("/")
        if relative == prefix or relative.startswith(prefix + "/"):
            candidates.append(prefix)
    if not candidates:
        return []
    longest = max(len(prefix) for prefix in candidates)
    return sorted(prefix for prefix in candidates if len(prefix) == longest)


def case_ids_in_sources(sources: Iterable[str], validation: Validation) -> Set[str]:
    identifiers: Set[str] = set()
    for raw in sources:
        try:
            path = repo_path(raw)
        except ValueError as error:
            validation.errors.append(str(error))
            continue
        validation.require(path.is_file(), "case registry source is missing: {}".format(raw))
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        identifiers.update(match.group(1) for match in STABLE_ID.finditer(text))
    return identifiers


def manifest_catalog(manifest: Mapping[str, Any], validation: Validation) -> Tuple[Set[str], Set[str]]:
    suites: Set[str] = set()
    cases: Set[str] = set()
    for index, suite in enumerate(manifest.get("tests", [])):
        suite_id = suite.get("id")
        validation.require(isinstance(suite_id, str) and bool(suite_id),
                           "manifest tests[{}] has no stable id".format(index))
        if not isinstance(suite_id, str) or not suite_id:
            continue
        validation.require(suite_id not in suites, "duplicate manifest suite: {}".format(suite_id))
        suites.add(suite_id)
        for case_id in suite.get("expected_cases", []):
            validation.require(case_id not in cases, "duplicate manifest case: {}".format(case_id))
            cases.add(case_id)
    return suites, cases


def validate_owner(owner: Mapping[str, Any], label: str, suites: Set[str],
                   cases: Set[str], source_cases: Set[str], validation: Validation) -> None:
    status = owner.get("status")
    validation.require(status in ALLOWED_STATUSES,
                       "{} has invalid status {!r}".format(label, status))
    validation.require(isinstance(owner.get("reason"), str) and bool(owner.get("reason", "").strip()),
                       "{} needs a nonempty reason".format(label))
    gaps = owner.get("gaps")
    validation.require(isinstance(gaps, list), "{} gaps must be an array".format(label))
    suite_ids = owner.get("suite_ids", [])
    case_ids = owner.get("case_ids", [])
    validation.require(isinstance(suite_ids, list), "{} suite_ids must be an array".format(label))
    validation.require(isinstance(case_ids, list), "{} case_ids must be an array".format(label))
    for suite_id in suite_ids if isinstance(suite_ids, list) else []:
        validation.require(suite_id in suites,
                           "{} names unknown manifest suite {}".format(label, suite_id))
    for case_id in case_ids if isinstance(case_ids, list) else []:
        validation.require(case_id in cases,
                           "{} names case absent from manifest {}".format(label, case_id))
        validation.require(case_id in source_cases,
                           "{} names case absent from declared registry sources {}".format(
                               label, case_id))
    if status == "tested":
        validation.require(bool(case_ids), "{} tested status requires direct cases".format(label))
        validation.require(not gaps, "{} tested status cannot retain gaps".format(label))
    elif status in {"partial", "planned"}:
        validation.require(bool(gaps), "{} {} status requires explicit gaps".format(label, status))


def validate_traceability(trace: Mapping[str, Any], manifest: Mapping[str, Any],
                          inventory: Mapping[str, Any], validation: Validation) -> Dict[str, Any]:
    validation.require(trace.get("schema_version") == 1,
                       "source_traceability.json schema_version must be 1")
    validation.require(set(trace.get("statuses", [])) == ALLOWED_STATUSES,
                       "source_traceability.json statuses differ from the validator contract")
    suites, cases = manifest_catalog(manifest, validation)
    source_cases = case_ids_in_sources(trace.get("case_registry_sources", []), validation)
    for case_id in sorted(cases):
        validation.require(case_id in source_cases,
                           "manifest case is absent from declared registry sources: {}".format(case_id))

    areas = trace.get("areas", [])
    validation.require(isinstance(areas, list) and bool(areas),
                       "source traceability areas must be a nonempty array")
    area_ids: Set[str] = set()
    area_paths: Set[str] = set()
    components: Set[str] = set()
    area_by_path: Dict[str, Mapping[str, Any]] = {}
    for index, area in enumerate(areas if isinstance(areas, list) else []):
        label = "trace area[{}]".format(index)
        area_id = area.get("area_id")
        path = area.get("path")
        component = area.get("component")
        validation.require(isinstance(area_id, str) and bool(area_id),
                           "{} needs area_id".format(label))
        validation.require(isinstance(path, str) and bool(path), "{} needs path".format(label))
        validation.require(isinstance(component, str) and bool(component),
                           "{} needs component".format(label))
        if isinstance(area_id, str):
            validation.require(area_id not in area_ids, "duplicate trace area_id: {}".format(area_id))
            area_ids.add(area_id)
        if isinstance(path, str):
            validation.require(path not in area_paths, "duplicate trace path: {}".format(path))
            area_paths.add(path)
            area_by_path[path] = area
            try:
                validation.require(repo_path(path).is_dir(), "trace area path is missing: {}".format(path))
            except ValueError as error:
                validation.errors.append(str(error))
        if isinstance(component, str):
            validation.require(component not in components,
                               "duplicate trace component: {}".format(component))
            components.add(component)
        validate_owner(area, label, suites, cases, source_cases, validation)
        for registry in area.get("local_registries", []):
            try:
                validation.require(repo_path(registry).is_file(),
                                   "{} local registry is missing: {}".format(label, registry))
            except ValueError as error:
                validation.errors.append(str(error))

    inventory_entries = inventory.get("components", [])
    inventory_by_component = {entry.get("component"): entry for entry in inventory_entries}
    validation.require(len(inventory_by_component) == len(inventory_entries),
                       "component_inventory.json has duplicate component entries")
    validation.require(components == set(inventory_by_component),
                       "trace components and component inventory differ: trace-only={} inventory-only={}".format(
                           sorted(components - set(inventory_by_component)),
                           sorted(set(inventory_by_component) - components)))
    for area in areas if isinstance(areas, list) else []:
        component = area.get("component")
        entry = inventory_by_component.get(component, {})
        expected_status = INVENTORY_STATUS.get(area.get("status"))
        validation.require(entry.get("status") == expected_status,
                           "{} status {} disagrees with inventory {}".format(
                               component, area.get("status"), entry.get("status")))
        validation.require(set(entry.get("test_ids", [])) == set(area.get("suite_ids", [])),
                           "{} suite ownership disagrees with component inventory".format(component))

    translation_units, public_contracts = discover_files(trace.get("discovery", {}))
    expanded = []
    counts: Counter = Counter()
    for kind, paths in (("translation_unit", translation_units),
                        ("public_contract", public_contracts)):
        for relative in paths:
            owners = owning_area(relative, areas)
            validation.require(len(owners) == 1,
                               "{} has {} owners: {}".format(relative, len(owners), owners))
            if len(owners) != 1:
                continue
            area = area_by_path[owners[0]]
            status = area["status"]
            counts[(kind, status)] += 1
            expanded.append({
                "path": relative,
                "kind": kind,
                "area_id": area["area_id"],
                "status": status,
            })

    return {
        "schema_version": 1,
        "translation_unit_count": len(translation_units),
        "public_contract_count": len(public_contracts),
        "area_count": len(areas) if isinstance(areas, list) else 0,
        "manifest_suite_count": len(suites),
        "manifest_case_count": len(cases),
        "counts": {
            kind: {status: counts[(kind, status)] for status in sorted(ALLOWED_STATUSES)}
            for kind in ("translation_unit", "public_contract")
        },
        "files": expanded,
    }


def validate_contracts(contracts: Mapping[str, Any], manifest: Mapping[str, Any],
                       trace: Mapping[str, Any], validation: Validation) -> Dict[str, int]:
    validation.require(contracts.get("schema_version") == 1,
                       "configuration_contracts.json schema_version must be 1")
    validation.require(set(contracts.get("statuses", [])) == ALLOWED_STATUSES,
                       "configuration_contracts.json statuses differ from validator contract")
    suites, cases = manifest_catalog(manifest, Validation())
    source_cases = case_ids_in_sources(trace.get("case_registry_sources", []), validation)
    known_components = {area.get("component") for area in trace.get("areas", [])}
    identifiers: Set[str] = set()
    status_counts: Counter = Counter()
    for index, contract in enumerate(contracts.get("contracts", [])):
        label = "configuration contract[{}]".format(index)
        mode_id = contract.get("mode_id")
        validation.require(isinstance(mode_id, str) and bool(mode_id),
                           "{} needs mode_id".format(label))
        if isinstance(mode_id, str):
            validation.require(mode_id not in identifiers,
                               "duplicate configuration mode_id: {}".format(mode_id))
            identifiers.add(mode_id)
        validate_owner(contract, label, suites, cases, source_cases, validation)
        status_counts[contract.get("status")] += 1
        components = contract.get("components", [])
        validation.require(isinstance(components, list) and bool(components),
                           "{} needs at least one owning component".format(label))
        for component in components if isinstance(components, list) else []:
            validation.require(component in known_components,
                               "{} names unknown component {}".format(label, component))
        authorities = contract.get("authority", [])
        validation.require(isinstance(authorities, list) and bool(authorities),
                           "{} needs at least one authority path".format(label))
        for authority in authorities if isinstance(authorities, list) else []:
            try:
                validation.require(repo_path(authority).exists(),
                                   "{} authority is missing: {}".format(label, authority))
            except ValueError as error:
                validation.errors.append(str(error))
        validation.require(isinstance(contract.get("profile"), str) and
                           bool(contract.get("profile", "").strip()),
                           "{} needs a profile".format(label))
        validation.require(isinstance(contract.get("selectors"), list),
                           "{} selectors must be an array".format(label))
        validation.require(isinstance(contract.get("generated_outputs"), list),
                           "{} generated_outputs must be an array".format(label))
    return {status: status_counts[status] for status in sorted(ALLOWED_STATUSES)}


def validate_test_catalog(
    tests: Sequence[Mapping[str, Any]],
    inventory: Mapping[str, Any],
    validation: Validation,
) -> None:
    """Require generated human documentation to match the live registry.

    The manifest remains authoritative. This comparison makes the contributor
    catalog an enforced projection of that authority, so an agent cannot add a
    case while leaving the documented what/how/why/where inventory stale.
    """

    catalog = TEST_ROOT / "TEST_CATALOG.md"
    try:
        expected = render_catalog(REPO_ROOT, TEST_ROOT, tests, inventory)
        actual = catalog.read_text(encoding="utf-8")
    except (OSError, ValueError) as error:
        validation.errors.append("test catalog: {}".format(error))
        return
    validation.require(
        actual == expected,
        "TEST_CATALOG.md is stale; run python3 src/test/generate_test_catalog.py")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        help="write expanded machine-readable ownership to this file")
    return parser.parse_args(argv)


def audit_repository() -> Tuple[Validation, Dict[str, Any]]:
    """Read and reconcile every governance ledger without executing tests.

    Both the standalone validator and the strict-coverage CLI use this exact
    path. An inventory row alone cannot establish coverage: the owning source
    area must agree, its named cases must exist, and a tested declaration may
    not retain unresolved gaps. This is declaration consistency, not a fresh
    physics qualification or a replacement for retained execution receipts.
    The function only reads repository files; the caller owns output and exit
    classification so malformed governance remains distinct from valid but
    incomplete coverage.
    """

    validation = Validation()
    trace = load_json(TEST_ROOT / "source_traceability.json", validation)
    contracts = load_json(TEST_ROOT / "configuration_contracts.json", validation)
    # Use the same manifest loader as the global runner.  T11 introduced
    # reviewable out-of-line case registries for the several-hundred-row model
    # imports; reading raw manifest JSON here would falsely treat those suites
    # as case-less and let the two governance paths disagree.
    try:
        loaded_tests = load_manifest(
            TEST_ROOT / "manifest.json", REPO_ROOT, TEST_ROOT)
        manifest = {"tests": loaded_tests}
    except ManifestError as error:
        validation.errors.append("manifest: {}".format(error))
        manifest = {"tests": []}
    inventory = load_json(TEST_ROOT / "component_inventory.json", validation)
    validate_test_catalog(manifest["tests"], inventory, validation)
    expanded = validate_traceability(trace, manifest, inventory, validation)
    contract_counts = validate_contracts(contracts, manifest, trace, validation)
    expanded["configuration_contract_count"] = sum(contract_counts.values())
    expanded["configuration_contract_status_counts"] = contract_counts
    return validation, expanded


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    validation, expanded = audit_repository()

    if validation.errors:
        for error in validation.errors:
            print("ERROR: {}".format(error), file=sys.stderr)
        print("traceability validation failed with {} error(s)".format(
            len(validation.errors)), file=sys.stderr)
        return 1

    if args.output:
        output = args.output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_name(output.name + ".tmp")
        temporary.write_text(json.dumps(expanded, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
        temporary.replace(output)
    print("PASS traceability: {} areas, {} production translation units, {} public contracts, "
          "{} configuration contracts, {} suites, {} cases".format(
              expanded["area_count"], expanded["translation_unit_count"],
              expanded["public_contract_count"], expanded["configuration_contract_count"],
              expanded["manifest_suite_count"], expanded["manifest_case_count"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

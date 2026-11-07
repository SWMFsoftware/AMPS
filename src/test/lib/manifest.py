"""Manifest and component-inventory validation for AMPS core tests."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Set, Tuple


class ManifestError(ValueError):
    """The test registry or coverage inventory is inconsistent."""


_ID = re.compile(r"^[A-Z][A-Z0-9-]*$")
_LEVELS = {"unit", "module", "integration", "system"}
CAMPAIGNS = (
    "smoke", "nightly", "profile", "extended", "distributed",
    "application-smoke",
)
_CAMPAIGNS = set(CAMPAIGNS)
_COVERAGE = {"covered", "partial", "planned", "exempt"}


def _load(path: Path) -> Dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as stream:
            value = json.load(stream)
    except OSError as error:
        raise ManifestError("cannot read {}: {}".format(path, error))
    except ValueError as error:
        raise ManifestError("invalid JSON in {}: {}".format(path, error))
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise ManifestError("{} must be a schema_version 1 object".format(path))
    return value


def load_manifest(path: Path, repo_root: Path, test_root: Path) -> List[Dict[str, Any]]:
    """Load the execution registry and resolve runners within ``test_root``.

    The containment checks are both correctness and safety boundaries: a
    committed manifest cannot execute an arbitrary file outside the test tree
    or claim coverage for a path outside production ``src``.
    """

    payload = _load(path)
    unknown_top = set(payload) - {"schema_version", "campaigns", "tests"}
    if unknown_top:
        raise ManifestError("{} has unknown top-level keys: {}".format(
            path, ", ".join(sorted(unknown_top))))
    raw_campaigns = payload.get("campaigns", list(CAMPAIGNS))
    if (not isinstance(raw_campaigns, list) or not raw_campaigns or
            any(value not in _CAMPAIGNS for value in raw_campaigns) or
            len(set(raw_campaigns)) != len(raw_campaigns)):
        raise ManifestError("{}: campaigns must be a unique nonempty subset of {}".format(
            path, sorted(_CAMPAIGNS)))
    enabled_campaigns = set(raw_campaigns)
    raw_tests = payload.get("tests")
    if not isinstance(raw_tests, list):
        raise ManifestError("{}: 'tests' must be a list".format(path))
    tests: List[Dict[str, Any]] = []
    identifiers: Set[str] = set()
    all_case_ids: Set[str] = set()
    runner_root = test_root.resolve()
    for index, raw in enumerate(raw_tests):
        label = "{}: tests[{}]".format(path, index)
        if not isinstance(raw, dict):
            raise ManifestError("{} must be an object".format(label))
        required = {
            "id", "component", "level", "tiers", "runner", "description",
            "expected_cases",
        }
        missing = required - set(raw)
        if missing:
            raise ManifestError("{} missing {}".format(label, ", ".join(sorted(missing))))
        allowed = required | {
            "aliases", "alias_cases", "requires", "timeout_seconds",
            "mpi_ranks", "exclusive",
        }
        unknown_keys = set(raw) - allowed
        if unknown_keys:
            raise ManifestError("{} has unknown keys: {}".format(
                label, ", ".join(sorted(unknown_keys))))
        identifier = raw["id"]
        if not isinstance(identifier, str) or not _ID.match(identifier):
            raise ManifestError("{} has invalid id {!r}".format(label, identifier))
        if identifier in identifiers:
            raise ManifestError("duplicate test id '{}'".format(identifier))
        identifiers.add(identifier)
        aliases = raw.get("aliases", [])
        if (not isinstance(aliases, list) or
                any(not isinstance(value, str) or not _ID.match(value)
                    for value in aliases) or
                len(set(aliases)) != len(aliases)):
            raise ManifestError("{} aliases must be unique test IDs".format(label))
        collisions = set(aliases) & identifiers
        if collisions:
            raise ManifestError("duplicate test ID or alias '{}'".format(
                sorted(collisions)[0]))
        identifiers.update(aliases)
        alias_cases = raw.get("alias_cases", {})
        if not isinstance(alias_cases, dict):
            raise ManifestError("{} alias_cases must be an object".format(label))
        unknown_alias_cases = set(alias_cases) - set(aliases)
        if unknown_alias_cases:
            raise ManifestError("{} alias_cases names undeclared aliases: {}".format(
                label, ", ".join(sorted(unknown_alias_cases))))
        component = raw["component"]
        if not isinstance(component, str) or not component.startswith("src/"):
            raise ManifestError("{} component must be below src/".format(label))
        component_path = (repo_root / component).resolve()
        try:
            component_path.relative_to((repo_root / "src").resolve())
        except ValueError:
            raise ManifestError("{} component escapes src/".format(label))
        if not component_path.exists():
            raise ManifestError("{} component does not exist: {}".format(label, component_path))
        if raw["level"] not in _LEVELS:
            raise ManifestError("{} has invalid level {!r}".format(label, raw["level"]))
        tiers = raw["tiers"]
        if (not isinstance(tiers, list) or not tiers or
                any(tier not in enabled_campaigns for tier in tiers) or
                len(set(tiers)) != len(tiers)):
            raise ManifestError("{} tiers must be a unique nonempty subset of {}".format(
                label, sorted(enabled_campaigns)))
        expected_cases = raw["expected_cases"]
        if (not isinstance(expected_cases, list) or not expected_cases or
                any(not isinstance(value, str) or not _ID.match(value)
                    for value in expected_cases) or
                len(set(expected_cases)) != len(expected_cases)):
            raise ManifestError("{} expected_cases must be unique nonempty case IDs".format(
                label))
        duplicate_cases = set(expected_cases) & all_case_ids
        if duplicate_cases:
            raise ManifestError("case IDs must be globally unique; duplicate {}".format(
                ", ".join(sorted(duplicate_cases))))
        all_case_ids.update(expected_cases)
        for alias, alias_expected in alias_cases.items():
            if (not isinstance(alias_expected, list) or not alias_expected or
                    len(set(alias_expected)) != len(alias_expected) or
                    any(value not in expected_cases for value in alias_expected)):
                raise ManifestError(
                    "{} alias_cases.{} must be a unique nonempty subset of expected_cases".format(
                        label, alias))
        runner = (test_root / str(raw["runner"])).resolve()
        try:
            runner.relative_to(runner_root)
        except ValueError:
            raise ManifestError("{} runner escapes {}".format(label, runner_root))
        if not runner.is_file():
            raise ManifestError("{} runner not found: {}".format(label, runner))
        requires = raw.get("requires", {})
        if not isinstance(requires, dict):
            raise ManifestError("{} requires must be an object".format(label))
        tool_names = requires.get("tools", [])
        path_names = requires.get("paths", [])
        capability_names = requires.get("capabilities", [])
        unknown_requirements = set(requires) - {"tools", "paths", "capabilities"}
        if unknown_requirements:
            raise ManifestError("{} requires has unknown keys: {}".format(
                label, ", ".join(sorted(unknown_requirements))))
        if (not isinstance(tool_names, list) or
                not all(isinstance(v, str) and v for v in tool_names) or
                len(set(tool_names)) != len(tool_names)):
            raise ManifestError("{} requires.tools must be a unique list of names".format(label))
        if (not isinstance(path_names, list) or
                not all(isinstance(v, str) and v for v in path_names) or
                len(set(path_names)) != len(path_names)):
            raise ManifestError(
                "{} requires.paths must be a unique list of repository paths".format(label))
        if (not isinstance(capability_names, list) or
                not all(isinstance(v, str) and v for v in capability_names) or
                len(set(capability_names)) != len(capability_names)):
            raise ManifestError(
                "{} requires.capabilities must be a unique list of names".format(label))
        timeout = raw.get("timeout_seconds")
        if timeout is not None and (not isinstance(timeout, int) or timeout <= 0):
            raise ManifestError("{} timeout_seconds must be positive".format(label))
        mpi_ranks = raw.get("mpi_ranks")
        if mpi_ranks is not None and (not isinstance(mpi_ranks, int) or mpi_ranks <= 0):
            raise ManifestError("{} mpi_ranks must be positive".format(label))
        if "exclusive" in raw and not isinstance(raw["exclusive"], bool):
            raise ManifestError("{} exclusive must be boolean".format(label))
        item = dict(raw)
        item["runner_path"] = str(runner)
        item["requires"] = {
            "tools": tool_names, "paths": path_names,
            "capabilities": capability_names,
        }
        item["aliases"] = list(aliases)
        item["alias_cases"] = {
            str(alias): list(values) for alias, values in alias_cases.items()}
        item["campaigns"] = list(tiers)
        item["expected_cases"] = list(expected_cases)
        item["exclusive"] = bool(raw.get("exclusive", False))
        tests.append(item)
    return tests


def discover_components(repo_root: Path) -> Set[str]:
    """Return production component directories governed by the inventory.

    Discovery is deliberately one level deep for the major ownership roots.
    That granularity matches independently testable AMPS libraries without
    turning every implementation subdirectory into a separate coverage gate.
    """

    src = repo_root / "src"
    components = {
        "src/general", "src/interface", "src/meshAMR", "src/pic", "src/species",
        "src/test",
    }
    for parent in (src / "general", src / "pic", src / "models"):
        if parent.is_dir():
            for child in parent.iterdir():
                if child.is_dir() and child.name not in {"build", "test", "__pycache__"}:
                    components.add(child.relative_to(repo_root).as_posix())
    return components


def audit_inventory(
    path: Path,
    repo_root: Path,
    manifest_tests: Sequence[Mapping[str, Any]],
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Validate inventory structure and return ``(entries, errors)``.

    Coverage may legitimately remain partial during migration.  Structural
    errors (missing/new component, unknown test ID, inconsistent covered
    state) are always errors; incomplete states are reported separately by
    :func:`incomplete_components`.
    """

    payload = _load(path)
    raw_entries = payload.get("components")
    if not isinstance(raw_entries, list):
        raise ManifestError("{}: 'components' must be a list".format(path))
    tests = {str(item["id"]) for item in manifest_tests}
    for item in manifest_tests:
        tests.update(str(value) for value in item.get("aliases", []))
    entries: List[Dict[str, Any]] = []
    errors: List[str] = []
    seen: Set[str] = set()
    for index, raw in enumerate(raw_entries):
        if not isinstance(raw, dict):
            errors.append("components[{}] is not an object".format(index))
            continue
        component = raw.get("component")
        status = raw.get("status")
        ids = raw.get("test_ids", [])
        if not isinstance(component, str):
            errors.append("components[{}] has no string component".format(index))
            continue
        if component in seen:
            errors.append("duplicate inventory component {}".format(component))
        seen.add(component)
        if status not in _COVERAGE:
            errors.append("{} has invalid coverage state {!r}".format(component, status))
        if not isinstance(ids, list) or not all(isinstance(value, str) for value in ids):
            errors.append("{} test_ids must be a list of IDs".format(component))
            ids = []
        unknown = set(ids) - tests
        if unknown:
            errors.append("{} names unknown tests: {}".format(component, ", ".join(sorted(unknown))))
        if status == "covered" and not ids:
            errors.append("{} is covered but has no registered tests".format(component))
        entries.append(dict(raw))

    discovered = discover_components(repo_root)
    missing = discovered - seen
    stale = seen - discovered
    if missing:
        errors.append("components missing from inventory: {}".format(", ".join(sorted(missing))))
    if stale:
        errors.append("inventory paths not discovered in src: {}".format(", ".join(sorted(stale))))
    return entries, errors


def incomplete_components(entries: Sequence[Mapping[str, Any]]) -> List[str]:
    """Return every component that cannot yet satisfy strict coverage."""

    return [str(entry.get("component")) for entry in entries
            if entry.get("status") not in ("covered", "exempt")]

#!/usr/bin/env python3
"""Fail-closed data contract for the Roadmap Step-12 release decision.

This module intentionally contains validation and provenance rules only.  It never
runs a solver and never changes an acceptance threshold.  Production calculations can
therefore be expensive and site-specific while the final audit remains deterministic,
portable, and suitable for an independent reviewer.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


SCHEMA = "sep-in-geospace/step12-release/v1"
SUMMARY_SCHEMA = "sep-in-geospace/step12-release-summary/v1"
EVIDENCE_SCHEMA = "sep-in-geospace/gate-evidence/v1"
BUILD_SCHEMA = "sep-in-geospace/build-provenance/v1"
CAPABILITY_SCHEMA = "sep-in-geospace/capabilities/v1"
RESOURCE_SCHEMA = "sep-in-geospace/resource-estimates/v1"
COMMON_PHYSICS_TAG = "sep-in-geospace-phase1-static-characteristics-v1"
PHASE_1_SCOPE = "INSTANTANEOUS_QUASI_STATIC_MAGNETIC"
NORMALIZATION_POLICY = "SHARED_EVENT_BOUNDARY_NO_PLATFORM_SCALE"

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
UTC_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}(:\d{2}(?:\.\d+)?)?Z?$")
PLACEHOLDER_RE = re.compile(
    r"(?:<[^>]+>|\b(?:REPLACE|TODO|TBD|INSERT|CHANGEME)(?:_|\b))", re.IGNORECASE
)

# Step 12 closes only the Phase-1 static/quasi-static product.  Dynamic electric-field
# characteristics (U-F11/I-F09) and forward/trapping/loss development are explicitly
# outside that scope.  The registry is fixed in code so deleting a difficult gate from
# a site manifest cannot make an incomplete campaign pass.
REQUIRED_GATES: Mapping[str, Tuple[str, ...]] = {
    "U": tuple("U-F%02d" % n for n in list(range(1, 11)) + [12, 13]),
    "I": tuple("I-F%02d" % n for n in list(range(1, 8)) + [10]),
    "C": tuple("C%d" % n for n in range(1, 20)),
    "F": tuple("F%d" % n for n in range(1, 18)),
    "O": ("O1", "O2", "O3", "O4"),
}
REQUIRED_GATE_IDS = tuple(
    gate for family in ("U", "I", "C", "F", "O") for gate in REQUIRED_GATES[family]
)

REQUIRED_PARITY_ROLES = (
    "FIELD_SAMPLES",
    "TRAJECTORIES",
    "DIRECTIONAL_ACCESS",
    "CUTOFF",
    "SPECTRUM",
    "DENSITY",
    "DETECTOR",
)
REQUIRED_CAMPAIGN_KINDS = ("SAMPLED_ANALYTIC_MESH", "SWMF_OFFLINE_REPLAY")

# These are ceilings, not defaults.  A manifest may be stricter, but cannot raise a
# threshold above the validation plan.  Absolute tolerances remain at the numerical
# floor and cannot be used to hide a large relative discrepancy near a physical scale.
TOLERANCE_CEILINGS: Mapping[str, Tuple[float, float]] = {
    "EXACT": (0.0, 0.0),
    "KERNEL": (1.0e-10, 1.0e-12),
    "INTEGRATED": (0.02, 1.0e-12),
    "MESH": (0.05, 1.0e-12),
    "DETECTOR": (0.05, 1.0e-12),
}

REQUIRED_SUPPORTED_CAPABILITIES = {
    "STANDALONE_CUTOFF",
    "STANDALONE_FLUX_SPECTRUM",
    "SWMF_CUTOFF",
    "SWMF_FLUX_SPECTRUM",
    "OFFLINE_SWMF_REPLAY",
}
REQUIRED_UNSUPPORTED_CAPABILITIES = {
    "DYNAMIC_EB_CHARACTERISTICS",
    "LONG_DURATION_TRAPPING",
    "LOCAL_ACCELERATION",
    "LOSS_PHYSICS",
    "FORWARD_TRAPPED_FLUX",
}
REQUIRED_RESOURCE_CASES = {
    "STANDALONE_CUTOFF",
    "STANDALONE_COMBINED",
    "SWMF_CUTOFF",
    "SWMF_COMBINED",
    "CROSS_PATH_PARITY",
}
REQUIRED_COMMANDS = {
    "standalone_cutoff_only",
    "standalone_combined",
    "swmf_cutoff_only",
    "swmf_combined",
    "offline_replay",
}
REQUIRED_HOOKS = {
    "ccmc_archive",
    "swmf_snapshot_export",
    "swmf_offline_replay",
}


class ReleaseValidationError(RuntimeError):
    """A missing, malformed, mutable, or scientifically invalid release input."""


def require(condition: bool, message: str) -> None:
    """Raise one uniform exception so every CLI can produce a structured failure."""

    if not condition:
        raise ReleaseValidationError(message)


def require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    require(isinstance(value, dict), "%s must be a JSON object" % label)
    return value


def require_sequence(value: Any, label: str) -> Sequence[Any]:
    require(isinstance(value, list), "%s must be a JSON array" % label)
    return value


def require_nonempty_string(value: Any, label: str, allow_placeholder: bool = False) -> str:
    require(isinstance(value, str) and bool(value.strip()), "%s must be a non-empty string" % label)
    text = value.strip()
    if not allow_placeholder:
        require(not PLACEHOLDER_RE.search(text), "%s contains an unresolved placeholder" % label)
    return text


def require_sha256(value: Any, label: str) -> str:
    text = require_nonempty_string(value, label)
    require(bool(SHA256_RE.fullmatch(text)), "%s must be a lowercase SHA-256" % label)
    return text


def read_json(path: Path, label: str) -> Mapping[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise ReleaseValidationError("cannot read %s %s: %s" % (label, path, exc)) from exc
    return require_mapping(value, label)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise ReleaseValidationError("cannot hash %s: %s" % (path, exc)) from exc
    return digest.hexdigest()


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    """Publish a complete report or leave the previous report untouched."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


class ResourceResolver:
    """Resolve only local, hash-pinned inputs relative to the release manifest.

    Network locations and unpinned paths are rejected.  Every verified digest is kept
    so the caller can build a restart fingerprint covering all consumed bytes.
    """

    def __init__(self, manifest_path: Path):
        self.manifest_path = manifest_path.expanduser().resolve()
        self.base = self.manifest_path.parent
        self.verified: MutableMapping[str, str] = {}

    def resolve(self, reference: Any, label: str) -> Path:
        item = require_mapping(reference, label)
        raw_path = require_nonempty_string(item.get("path"), label + ".path")
        require("://" not in raw_path, "%s.path must be a local file" % label)
        expected = require_sha256(item.get("sha256"), label + ".sha256")
        candidate = Path(raw_path).expanduser()
        if not candidate.is_absolute():
            candidate = self.base / candidate
        candidate = candidate.resolve()
        require(candidate.is_file(), "%s does not exist: %s" % (label, candidate))
        # A single release manifest often references the same immutable artifact from
        # several gates.  Reuse a digest already verified in this process, while still
        # requiring every declaration to carry the identical expected hash.
        cached = self.verified.get(str(candidate))
        if cached is not None:
            require(cached == expected, "%s declares a conflicting SHA-256 for %s" % (label, candidate))
            return candidate
        observed = sha256_file(candidate)
        require(observed == expected, "%s SHA-256 mismatch for %s" % (label, candidate))
        self.verified[str(candidate)] = observed
        return candidate

    def fingerprint(self, manifest_sha256: str) -> str:
        material = {
            "manifest_sha256": manifest_sha256,
            "resources": sorted(self.verified.items()),
        }
        return hashlib.sha256(canonical_json_bytes(material)).hexdigest()


def validate_build_manifest(
    value: Mapping[str, Any], expected_kind: str, tag: str, source_sha256: str
) -> Dict[str, Any]:
    label = expected_kind + " build manifest"
    require(value.get("schema") == BUILD_SCHEMA, "%s has wrong schema" % label)
    require(value.get("build_kind") == expected_kind, "%s has wrong build_kind" % label)
    require(value.get("common_physics_tag") == tag, "%s uses a different common physics tag" % label)
    require(
        value.get("common_physics_source_sha256") == source_sha256,
        "%s uses a different common physics source digest" % label,
    )
    revision = require_nonempty_string(value.get("source_revision"), label + ".source_revision")
    require(value.get("dirty") is False, "%s must be built from a clean recorded revision" % label)
    compiler = require_nonempty_string(value.get("compiler"), label + ".compiler")
    dependencies = require_mapping(value.get("dependencies"), label + ".dependencies")
    require(bool(dependencies), "%s.dependencies must not be empty" % label)
    for name, version in dependencies.items():
        require_nonempty_string(name, label + ".dependency name")
        require_nonempty_string(version, "%s.dependencies.%s" % (label, name))
    executable = require_sha256(value.get("executable_sha256"), label + ".executable_sha256")
    return {"kind": expected_kind, "revision": revision, "compiler": compiler, "executable_sha256": executable}


def _validate_artifact_references(
    artifacts: Any, resolver: ResourceResolver, label: str
) -> List[Dict[str, str]]:
    values = require_sequence(artifacts, label)
    require(bool(values), "%s must contain at least one closed artifact" % label)
    seen = set()
    output: List[Dict[str, str]] = []
    for index, reference in enumerate(values):
        path = resolver.resolve(reference, "%s[%d]" % (label, index))
        require(str(path) not in seen, "%s contains a duplicate artifact" % label)
        seen.add(str(path))
        output.append({"path": str(path), "sha256": resolver.verified[str(path)]})
    return output


def validate_gate_evidence(
    value: Mapping[str, Any], expected_gate: str, resolver: ResourceResolver, tag: str
) -> Dict[str, Any]:
    """Validate a completed runner record without interpreting a dry run as PASS."""

    label = "gate evidence %s" % expected_gate
    require(value.get("schema") == EVIDENCE_SCHEMA, "%s has wrong schema" % label)
    require(value.get("gate_id") == expected_gate, "%s has the wrong gate_id" % label)
    require(value.get("status") == "PASS", "%s is not PASS" % label)
    require(value.get("result_line") == "RESULT: PASS", "%s lacks the exact RESULT: PASS line" % label)
    require(value.get("exit_code") == 0, "%s did not exit zero" % label)
    require(value.get("dry_run") is False, "%s is a dry run and cannot count as physics evidence" % label)
    require(value.get("common_physics_tag") == tag, "%s uses a different common physics tag" % label)
    require_nonempty_string(value.get("command"), label + ".command")
    require_nonempty_string(value.get("source_revision"), label + ".source_revision")
    require(value.get("dirty") is False, "%s must record dirty=false" % label)
    require_nonempty_string(value.get("compiler"), label + ".compiler")
    dependencies = require_mapping(value.get("dependencies"), label + ".dependencies")
    require(bool(dependencies), "%s.dependencies must not be empty" % label)
    for key in ("executable_sha256", "input_sha256", "boundary_sha256", "response_sha256"):
        require_sha256(value.get(key), label + "." + key)
    for key in ("field_id", "snapshot_id", "species", "frame", "mover", "scheduler"):
        require_nonempty_string(value.get(key), label + "." + key)
    epoch = require_nonempty_string(value.get("epoch_utc"), label + ".epoch_utc")
    require(bool(UTC_RE.fullmatch(epoch)), "%s.epoch_utc is not an absolute ISO-like UTC" % label)
    require_mapping(value.get("integrator_tolerances"), label + ".integrator_tolerances")
    grid_history = require_sequence(value.get("grid_history"), label + ".grid_history")
    require(bool(grid_history), "%s.grid_history must retain every refinement level" % label)
    require(isinstance(value.get("mpi_ranks"), int) and value["mpi_ranks"] > 0, "%s.mpi_ranks must be positive" % label)
    require(isinstance(value.get("threads"), int) and value["threads"] > 0, "%s.threads must be positive" % label)
    seeds = require_sequence(value.get("random_seeds"), label + ".random_seeds")
    require(all(isinstance(seed, int) for seed in seeds), "%s.random_seeds must contain integers" % label)
    wall = value.get("wall_time_s")
    require(isinstance(wall, (int, float)) and math.isfinite(float(wall)) and float(wall) >= 0.0, "%s.wall_time_s is invalid" % label)
    counts = require_mapping(value.get("termination_counts"), label + ".termination_counts")
    require(bool(counts), "%s.termination_counts must not be empty" % label)
    require(all(isinstance(count, int) and count >= 0 for count in counts.values()), "%s termination counts must be non-negative integers" % label)
    artifacts = _validate_artifact_references(value.get("artifacts"), resolver, label + ".artifacts")

    # Observation-facing gates must prove that a comparison was produced.  This rule
    # directly prevents the historic C19/C9/C10 failure mode where the solver exited
    # but no comparison table existed.
    if expected_gate.startswith("O") or expected_gate in {"F8", "F9", "F10", "F17"}:
        require(isinstance(value.get("observation_comparison_count"), int) and value["observation_comparison_count"] > 0,
                "%s produced no observation comparison" % label)
        require(value.get("normalization_policy") == NORMALIZATION_POLICY,
                "%s uses an unapproved platform normalization" % label)
    if expected_gate == "O4":
        require(value.get("frozen") is True, "O4 evidence must be frozen")
        require(value.get("retuned") is False, "O4 evidence shows post-freeze retuning")
        require_nonempty_string(value.get("reported_outcome"), "O4.reported_outcome")

    return {"gate_id": expected_gate, "status": "PASS", "artifact_count": len(artifacts)}


def validate_capabilities(value: Mapping[str, Any]) -> Dict[str, Any]:
    require(value.get("schema") == CAPABILITY_SCHEMA, "capability table has wrong schema")
    rows = require_sequence(value.get("capabilities"), "capabilities")
    by_id: Dict[str, Mapping[str, Any]] = {}
    for index, raw in enumerate(rows):
        row = require_mapping(raw, "capabilities[%d]" % index)
        capability_id = require_nonempty_string(row.get("id"), "capabilities[%d].id" % index)
        require(capability_id not in by_id, "duplicate capability %s" % capability_id)
        status = row.get("status")
        require(status in ("SUPPORTED", "UNSUPPORTED"), "capability %s has invalid status" % capability_id)
        require_nonempty_string(row.get("basis"), "capability %s basis" % capability_id)
        by_id[capability_id] = row
    for capability_id in REQUIRED_SUPPORTED_CAPABILITIES:
        require(capability_id in by_id and by_id[capability_id]["status"] == "SUPPORTED",
                "required supported capability is missing: %s" % capability_id)
    for capability_id in REQUIRED_UNSUPPORTED_CAPABILITIES:
        require(capability_id in by_id and by_id[capability_id]["status"] == "UNSUPPORTED",
                "Phase-1 limitation is not explicit: %s" % capability_id)
    return {"supported": sorted(REQUIRED_SUPPORTED_CAPABILITIES), "unsupported": sorted(REQUIRED_UNSUPPORTED_CAPABILITIES)}


def validate_resource_estimates(value: Mapping[str, Any]) -> Dict[str, Any]:
    require(value.get("schema") == RESOURCE_SCHEMA, "resource estimate table has wrong schema")
    rows = require_sequence(value.get("estimates"), "resource estimates")
    found = set()
    for index, raw in enumerate(rows):
        row = require_mapping(raw, "estimates[%d]" % index)
        case = require_nonempty_string(row.get("case"), "estimates[%d].case" % index)
        require(case not in found, "duplicate resource estimate %s" % case)
        found.add(case)
        for key in ("mpi_ranks", "threads", "memory_gib", "wall_minutes", "disk_gib"):
            number = row.get(key)
            require(isinstance(number, (int, float)) and math.isfinite(float(number)) and float(number) > 0.0,
                    "resource estimate %s.%s must be positive" % (case, key))
        require_nonempty_string(row.get("basis"), "resource estimate %s.basis" % case)
    missing = REQUIRED_RESOURCE_CASES - found
    require(not missing, "resource estimates missing cases: %s" % ", ".join(sorted(missing)))
    return {"cases": sorted(found)}


def validate_named_commands(value: Any, required: Iterable[str], label: str) -> Dict[str, str]:
    commands = require_mapping(value, label)
    output: Dict[str, str] = {}
    for name in sorted(required):
        require(name in commands, "%s is missing %s" % (label, name))
        output[name] = require_nonempty_string(commands[name], "%s.%s" % (label, name))
    return output


def validate_release_manifest_shape(manifest: Mapping[str, Any]) -> None:
    """Check immutable top-level semantics before opening referenced files."""

    require(manifest.get("schema") == SCHEMA, "release manifest has wrong schema")
    require_nonempty_string(manifest.get("release_id"), "release_id")
    require(manifest.get("phase_1_scope") == PHASE_1_SCOPE, "release manifest has wrong Phase-1 scope")
    common = require_mapping(manifest.get("common_physics"), "common_physics")
    require(common.get("tag") == COMMON_PHYSICS_TAG, "common_physics.tag is not the released tag")
    require_sha256(common.get("source_sha256"), "common_physics.source_sha256")
    require_mapping(common.get("standalone_build"), "common_physics.standalone_build")
    require_mapping(common.get("coupled_build"), "common_physics.coupled_build")
    require_sequence(manifest.get("gate_evidence"), "gate_evidence")
    require_sequence(manifest.get("cross_path_campaigns"), "cross_path_campaigns")
    require_mapping(manifest.get("holdout"), "holdout")
    require_mapping(manifest.get("capabilities"), "capabilities reference")
    require_mapping(manifest.get("resource_estimates"), "resource_estimates reference")
    validate_named_commands(manifest.get("hooks"), REQUIRED_HOOKS, "hooks")
    validate_named_commands(manifest.get("commands"), REQUIRED_COMMANDS, "commands")


def manifest_sha256(path: Path) -> str:
    return sha256_file(path.expanduser().resolve())

#!/usr/bin/env python3
"""Shared policy, hashing, and archive helpers for the Stage-3 science gate."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List


POLICY_SCHEMA = "amps.comet.stage3-science-policy.v1"
MANIFEST_SCHEMA = "amps.comet.stage3-science-manifest.v1"
IMPORTANT_MODES = ("BERNOULLI", "SYSTEMATIC")


class GateError(RuntimeError):
    """A policy, archive, or evidence contract was violated."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_json(value: object) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def atomic_write_json(path: Path, value: object) -> None:
    """Write JSON with rename semantics so a reboot cannot leave a half file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _finite_number(policy: Dict[str, Any], name: str, *, positive: bool = True) -> float:
    value = policy.get(name)
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise GateError(f"policy {name} must be numeric")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0.0):
        raise GateError(f"policy {name} must be finite and positive")
    return result


def validate_policy(policy: Dict[str, Any]) -> Dict[str, Any]:
    """Validate the approved policy without silently supplying defaults."""
    if policy.get("schema") != POLICY_SCHEMA:
        raise GateError("unknown Stage-3 science policy schema")
    if policy.get("modes") != list(IMPORTANT_MODES):
        raise GateError("policy modes must be BERNOULLI then SYSTEMATIC")
    seeds = policy.get("seeds")
    if (
        not isinstance(seeds, list)
        or len(seeds) != 4
        or any(not isinstance(seed, int) or isinstance(seed, bool) for seed in seeds)
        or any(seed < 1 or seed > 2147483000 for seed in seeds)
        or len(set(seeds)) != 4
    ):
        raise GateError("policy requires four distinct valid integer seeds")
    if policy.get("mpi_ranks") != 12 or policy.get("omp_threads") != 1:
        raise GateError("approved science gate requires 12 MPI ranks and one OpenMP thread")
    if policy.get("sampling_length") != 10000:
        raise GateError("approved sampling length is 10000")
    if policy.get("burn_in_iterations") != 10000:
        raise GateError("approved burn-in is 10000 iterations")
    if policy.get("maximum_iterations") != 40000:
        raise GateError("approved hard limit is 40000 iterations")
    if policy.get("required_converged_windows") != 3:
        raise GateError("approved convergence history requires three comparisons")
    if not math.isclose(_finite_number(policy, "convergence_tolerance"), 0.05):
        raise GateError("approved convergence tolerance is 0.05")
    altitude = policy.get("comparison_altitude_km")
    if altitude != [10.0, 100000.0]:
        raise GateError("approved altitude interval is 10--100000 km")
    if policy.get("require_fully_contained_cells") is not True:
        raise GateError("science gate must retain only fully contained exact cells")
    if policy.get("common_support_rule") != "all-runs-minimum":
        raise GateError("science gate requires all-runs-minimum exposure support")
    species = policy.get("important_species")
    if species != ["C2", "C2H", "CN", "C2H5", "HCO"]:
        raise GateError("important-species list differs from the approved gate")
    for name in (
        "nucleus_radius_m",
        "minimum_particle_iteration_exposure",
        "systematic_cv_target",
        "systematic_cv_required_bin_fraction",
        "systematic_cv_hard_maximum",
        "systematic_raw_zero_cell_fraction_maximum",
        "systematic_population_required_block_fraction",
        "mode_mean_relative_tolerance",
        "mode_mean_combined_standard_error_multiplier",
        "direct_source_relative_tolerance",
        "rubin_maximum_regression_fraction",
    ):
        _finite_number(policy, name)
    for name in (
        "systematic_cv_target",
        "systematic_cv_required_bin_fraction",
        "systematic_cv_hard_maximum",
        "systematic_raw_zero_cell_fraction_maximum",
        "systematic_population_required_block_fraction",
        "mode_mean_relative_tolerance",
        "direct_source_relative_tolerance",
        "rubin_maximum_regression_fraction",
    ):
        if float(policy[name]) > 1.0:
            raise GateError(f"policy {name} must not exceed one")

    # These numbers are the user-approved scientific decision boundaries, not
    # plotting knobs.  The manifest hash proves internal consistency, but by
    # itself would also bless a newly hashed, relaxed policy.  Pin the values
    # here so changing a gate requires an explicit source/test review.
    approved_numeric = {
        "nucleus_radius_m": 6000.0,
        "minimum_particle_iteration_exposure": 100000.0,
        "systematic_cv_target": 0.25,
        "systematic_cv_required_bin_fraction": 0.90,
        "systematic_cv_hard_maximum": 0.50,
        "systematic_raw_zero_cell_fraction_maximum": 0.05,
        "systematic_population_required_block_fraction": 0.90,
        "mode_mean_relative_tolerance": 0.05,
        "mode_mean_combined_standard_error_multiplier": 2.0,
        "direct_source_relative_tolerance": 0.01,
        "rubin_maximum_regression_fraction": 0.10,
    }
    for name, approved in approved_numeric.items():
        if not math.isclose(float(policy[name]), approved, rel_tol=0.0, abs_tol=1e-15):
            raise GateError(f"policy {name} differs from the approved value {approved:g}")
    if policy.get("chemical_only_population_target") != [50, 200]:
        raise GateError("approved chemical-only population target is 50--200")
    if policy.get("population_metric") != "completed-window-sampled-mean-occupancy":
        raise GateError("population metric must use completed-window sampled occupancy")
    if policy.get("require_zero_emergency_events_in_accepted_windows") is not True:
        raise GateError("approved gate requires zero emergency executions")
    return policy


def load_policy(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GateError(f"cannot read science policy {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise GateError("science policy root must be an object")
    return validate_policy(value)


def expected_run_order(policy: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Alternate mode order between seeds without changing the sample set."""
    result: List[Dict[str, Any]] = []
    for index, seed in enumerate(policy["seeds"]):
        modes = IMPORTANT_MODES if index % 2 == 0 else tuple(reversed(IMPORTANT_MODES))
        for mode in modes:
            result.append({"seed": seed, "mode": mode})
    return result


def parse_scalar_input_records(path: Path) -> Dict[str, str]:
    """Read simple NAME=VALUE records while ignoring AMPS `!` comments."""
    records: Dict[str, str] = {}
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*(?:\([^)]*\))?)\s*=\s*(.*?)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        body = line.split("!", 1)[0]
        match = pattern.match(body)
        if match:
            name, value = match.groups()
            # Species-qualified source/reaction records legitimately repeat.
            # Preserve a duplicate sentinel so a scalar control that should be
            # unique still fails its exact-value check without rejecting the
            # canonical multi-species input grammar.
            if name in records:
                records[name] = "<DUPLICATE>"
            else:
                records[name] = value
    return records


def validate_canonical_input(path: Path, policy: Dict[str, Any]) -> None:
    records = parse_scalar_input_records(path)
    expected = {
        "RadialSamplingLength": str(policy["sampling_length"]),
        "RadialMaxIterations": str(policy["maximum_iterations"]),
        "RadialConvergenceWindows": str(policy["required_converged_windows"]),
        "CometDiagnostics": "on",
        "ChemistryMode": "on",
    }
    for name, value in expected.items():
        if records.get(name, "").lower() != value.lower():
            raise GateError(f"canonical input {name} does not match policy")
    numeric = {
        "NucleusRadius": policy["nucleus_radius_m"],
        "RadialConvergenceTolerance": policy["convergence_tolerance"],
    }
    for name, expected_value in numeric.items():
        try:
            actual = float(records[name])
        except (KeyError, ValueError) as exc:
            raise GateError(f"canonical input has invalid {name}") from exc
        if not math.isclose(actual, float(expected_value), rel_tol=0.0, abs_tol=1e-12):
            raise GateError(f"canonical input {name} does not match policy")
    compact = records.get("RadialParticleNumberPerCellRange(CHEMICAL_ONLY)", "")
    values = [part.strip() for part in compact.split("|")]
    if values != ["50", "200"]:
        raise GateError("canonical chemical-only population range is not 50--200")


def ensure_exact_keys(records: Iterable[Dict[str, Any]], policy: Dict[str, Any]) -> None:
    records = list(records)
    expected = {(item["mode"], item["seed"]) for item in expected_run_order(policy)}
    actual = {(item.get("mode"), item.get("seed")) for item in records}
    if actual != expected or len(records) != 8:
        raise GateError("manifest does not contain exactly the approved eight runs")

#!/usr/bin/env python3
"""Strict independent-seed science gate for comet chemistry Stage 9.

This module intentionally consumes compact, checksum-protected *run summaries*
instead of scraping human-readable logs.  A summary must be generated from a
completed transported run and must carry every quantity needed by the Stage-9
decision.  Missing data are invalid evidence: the gate never substitutes a
default, infers convergence from the process return code, or treats a smooth
profile as a statistically resolved result.

The schema is deliberately independent of a particular comet species list.
The manifest policy names primary and daughter species and supplies every
scientific threshold.  This keeps thresholds out of the parser and prevents a
future run from silently inheriting values chosen for the synthetic tests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple


MANIFEST_SCHEMA = "amps.comet.stage9-science-manifest.v1"
SUMMARY_SCHEMA = "amps.comet.stage9-run-summary.v1"
METHODS = ("ANALOG", "FORCED")


class EvidenceError(ValueError):
    """Raised when evidence is malformed, incomplete, or unauthenticated."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_object(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EvidenceError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise EvidenceError(f"JSON root is not an object: {path}")
    return value


def number(value: object, label: str, *, minimum: float | None = None) -> float:
    # ``bool`` is an ``int`` in Python.  Rejecting it explicitly avoids
    # accepting JSON ``true`` as a scientific value of one.
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EvidenceError(f"{label} is not numeric")
    result = float(value)
    if not math.isfinite(result):
        raise EvidenceError(f"{label} is not finite")
    if minimum is not None and result < minimum:
        raise EvidenceError(f"{label} is below {minimum}")
    return result


def integer(value: object, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise EvidenceError(f"{label} is not an integer >= {minimum}")
    return value


def string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise EvidenceError(f"{label} is not a nonempty string")
    return value


def boolean(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise EvidenceError(f"{label} is not boolean")
    return value


def exact_keys(value: Mapping[str, Any], expected: Iterable[str], label: str) -> None:
    expected_set = set(expected)
    actual = set(value)
    if actual != expected_set:
        missing = sorted(expected_set - actual)
        extra = sorted(actual - expected_set)
        raise EvidenceError(f"{label} keys differ: missing={missing} extra={extra}")


def relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0)


def close(left: float, right: float, tolerance: float = 5.0e-13) -> bool:
    return relative_error(left, right) <= tolerance


def sample_statistics(values: Sequence[float], multiplier: float) -> Dict[str, float]:
    if len(values) < 2:
        raise EvidenceError("an independent-seed statistic needs at least two values")
    mean = statistics.fmean(values)
    variance = statistics.variance(values)
    standard_error = math.sqrt(variance / len(values))
    return {
        "mean": mean,
        "variance": variance,
        "standard_error": standard_error,
        "confidence_low": mean - multiplier * standard_error,
        "confidence_high": mean + multiplier * standard_error,
    }


def validate_policy(raw: object) -> Dict[str, Any]:
    if not isinstance(raw, dict):
        raise EvidenceError("manifest policy is not an object")
    required = {
        "minimum_seeds",
        "mpi_ranks",
        "omp_threads",
        "burn_in_iterations",
        "minimum_steady_windows",
        "maximum_buffer_utilization",
        "confidence_multiplier",
        "minimum_forced_daughter_collisions",
        "flux_closure",
        "primary_species",
        "daughter_species",
    }
    exact_keys(raw, required, "policy")
    policy = dict(raw)
    policy["minimum_seeds"] = integer(raw["minimum_seeds"], "minimum_seeds", minimum=4)
    policy["mpi_ranks"] = integer(raw["mpi_ranks"], "mpi_ranks", minimum=1)
    policy["omp_threads"] = integer(raw["omp_threads"], "omp_threads", minimum=1)
    policy["burn_in_iterations"] = integer(
        raw["burn_in_iterations"], "burn_in_iterations", minimum=0
    )
    policy["minimum_steady_windows"] = integer(
        raw["minimum_steady_windows"], "minimum_steady_windows", minimum=2
    )
    maximum_buffer = number(
        raw["maximum_buffer_utilization"], "maximum_buffer_utilization", minimum=0.0
    )
    if not 0.0 < maximum_buffer < 1.0:
        raise EvidenceError("maximum_buffer_utilization must be in (0,1)")
    policy["maximum_buffer_utilization"] = maximum_buffer
    policy["confidence_multiplier"] = number(
        raw["confidence_multiplier"], "confidence_multiplier", minimum=0.0
    )
    if policy["confidence_multiplier"] == 0.0:
        raise EvidenceError("confidence_multiplier must be positive")
    policy["minimum_forced_daughter_collisions"] = integer(
        raw["minimum_forced_daughter_collisions"],
        "minimum_forced_daughter_collisions",
        minimum=1,
    )

    # Flux closure is a scientific decision, not an evaluator default.  The
    # independent-seed method treats the signed residual as an estimator and
    # asks whether its confidence interval contains zero.  The alternative is
    # deliberately stricter operationally: every accepted run must satisfy an
    # explicitly supplied absolute relative-residual bound.  Keeping the
    # method and its optional threshold in the authenticated manifest makes it
    # impossible to reinterpret an archive after the run.
    flux = raw["flux_closure"]
    if not isinstance(flux, dict):
        raise EvidenceError("flux_closure policy is not an object")
    exact_keys(flux, {"method", "relative_tolerance"}, "flux_closure policy")
    method = string(flux["method"], "flux_closure method")
    if method == "independent_seed_ci_contains_zero":
        if flux["relative_tolerance"] is not None:
            raise EvidenceError(
                "flux_closure relative_tolerance must be null for independent_seed_ci_contains_zero"
            )
        policy["flux_closure"] = {
            "method": method,
            "relative_tolerance": None,
        }
    elif method == "relative_tolerance":
        tolerance = number(
            flux["relative_tolerance"],
            "flux_closure relative_tolerance",
            minimum=0.0,
        )
        if tolerance <= 0.0:
            raise EvidenceError("flux_closure relative_tolerance must be positive")
        policy["flux_closure"] = {
            "method": method,
            "relative_tolerance": tolerance,
        }
    else:
        raise EvidenceError(f"unknown flux_closure method: {method}")

    primary = raw["primary_species"]
    daughters = raw["daughter_species"]
    if not isinstance(primary, dict) or not primary:
        raise EvidenceError("primary_species must be a nonempty object")
    if not isinstance(daughters, dict) or not daughters:
        raise EvidenceError("daughter_species must be a nonempty object")
    overlap = set(primary).intersection(daughters)
    if overlap:
        raise EvidenceError(f"primary/daughter species overlap: {sorted(overlap)}")

    primary_result: Dict[str, Dict[str, float]] = {}
    for species, settings in primary.items():
        string(species, "primary species")
        if not isinstance(settings, dict):
            raise EvidenceError(f"primary policy for {species} is not an object")
        exact_keys(settings, {"relative_tolerance", "combined_se_multiplier"}, f"primary {species}")
        primary_result[species] = {
            "relative_tolerance": number(
                settings["relative_tolerance"], f"{species} relative_tolerance", minimum=0.0
            ),
            "combined_se_multiplier": number(
                settings["combined_se_multiplier"],
                f"{species} combined_se_multiplier",
                minimum=0.0,
            ),
        }

    daughter_result: Dict[str, Dict[str, float]] = {}
    daughter_keys = {
        "relative_tolerance",
        "combined_se_multiplier",
        "minimum_ess_fraction",
        "maximum_zero_cell_fraction",
        "minimum_variance_cost_gain",
        "minimum_gain_bin_fraction",
    }
    for species, settings in daughters.items():
        string(species, "daughter species")
        if not isinstance(settings, dict):
            raise EvidenceError(f"daughter policy for {species} is not an object")
        exact_keys(settings, daughter_keys, f"daughter {species}")
        converted = {
            key: number(settings[key], f"{species} {key}", minimum=0.0)
            for key in daughter_keys
        }
        for key in ("minimum_ess_fraction", "maximum_zero_cell_fraction", "minimum_gain_bin_fraction"):
            if converted[key] > 1.0:
                raise EvidenceError(f"{species} {key} exceeds one")
        if converted["minimum_ess_fraction"] <= 0.0:
            raise EvidenceError(f"{species} minimum_ess_fraction must be positive")
        if converted["minimum_gain_bin_fraction"] <= 0.0:
            raise EvidenceError(f"{species} minimum_gain_bin_fraction must be positive")
        if converted["minimum_variance_cost_gain"] <= 1.0:
            raise EvidenceError(
                f"{species} minimum_variance_cost_gain must demonstrate improvement (>1)"
            )
        daughter_result[species] = converted
    policy["primary_species"] = primary_result
    policy["daughter_species"] = daughter_result
    return policy


def _summary_path(manifest_path: Path, relative: object) -> Path:
    relative_text = string(relative, "summary path")
    root = manifest_path.parent.resolve()
    candidate = (root / relative_text).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise EvidenceError(f"summary path escapes manifest directory: {relative_text}") from exc
    if not candidate.is_file():
        raise EvidenceError(f"summary is missing: {candidate}")
    return candidate


def _profile(raw: object, label: str) -> List[Tuple[float, float]]:
    if not isinstance(raw, list) or not raw:
        raise EvidenceError(f"{label} profile is empty or not a list")
    result: List[Tuple[float, float]] = []
    last_radius = -math.inf
    for index, row in enumerate(raw):
        if not isinstance(row, dict):
            raise EvidenceError(f"{label} profile row {index} is not an object")
        exact_keys(row, {"radius_m", "value"}, f"{label} profile row {index}")
        radius = number(row["radius_m"], f"{label} radius[{index}]", minimum=0.0)
        value = number(row["value"], f"{label} value[{index}]", minimum=0.0)
        if radius <= last_radius:
            raise EvidenceError(f"{label} radii are not strictly increasing")
        last_radius = radius
        result.append((radius, value))
    return result


def _validate_daughter_record(raw: object, label: str) -> Dict[str, float]:
    keys = {
        "zero_cell_fraction",
        "model_carriers",
        "sum_weights",
        "sum_squared_weights",
        "surface_source_s-1",
        "chemical_production_s-1",
        "inner_reemitted_s-1",
        "chemical_destruction_s-1",
        "inner_incident_s-1",
        "outer_escape_s-1",
        "net_flux_s-1",
    }
    if not isinstance(raw, dict):
        raise EvidenceError(f"{label} daughter record is not an object")
    exact_keys(raw, keys, f"{label} daughter record")
    # All physical rates are non-negative except the signed net residual.
    # Keeping that exception explicit is important: a statistically ordinary
    # small negative imbalance must reach the closure test, not be rejected as
    # malformed evidence.
    result = {
        key: number(
            raw[key],
            f"{label} {key}",
            minimum=None if key == "net_flux_s-1" else 0.0,
        )
        for key in keys
    }
    if result["zero_cell_fraction"] > 1.0:
        raise EvidenceError(f"{label} zero_cell_fraction exceeds one")
    carriers = result["model_carriers"]
    if carriers <= 0.0 or result["sum_weights"] <= 0.0 or result["sum_squared_weights"] <= 0.0:
        raise EvidenceError(f"{label} has no positive carrier/weight evidence")
    ess = result["sum_weights"] ** 2 / result["sum_squared_weights"]
    if ess > carriers * (1.0 + 5.0e-13):
        raise EvidenceError(f"{label} ESS exceeds model carrier count")
    result["effective_sample_size"] = ess
    result["ess_fraction"] = ess / carriers

    # Recompute both the accounting identity and the steady-state residual.
    # The former detects a malformed diagnostic.  The latter is evaluated
    # against the species-specific science tolerance below.
    recomputed_net = (
        result["surface_source_s-1"]
        + result["chemical_production_s-1"]
        + result["inner_reemitted_s-1"]
        - result["chemical_destruction_s-1"]
        - result["inner_incident_s-1"]
        - result["outer_escape_s-1"]
    )
    if not close(recomputed_net, result["net_flux_s-1"]):
        raise EvidenceError(f"{label} flux identity does not close")
    gross = sum(
        abs(result[key])
        for key in (
            "surface_source_s-1",
            "chemical_production_s-1",
            "inner_reemitted_s-1",
            "chemical_destruction_s-1",
            "inner_incident_s-1",
            "outer_escape_s-1",
        )
    )
    result["flux_closure_fraction"] = abs(recomputed_net) / max(gross, 1.0)
    result["signed_flux_fraction"] = recomputed_net / max(gross, 1.0)
    return result


def validate_summary(raw: Dict[str, Any], policy: Dict[str, Any], label: str) -> Dict[str, Any]:
    required = {
        "schema",
        "method",
        "seed",
        "mpi_ranks",
        "omp_threads",
        "physical_input_sha256",
        "mesh_sha256",
        "reaction_network_sha256",
        "started_from_iteration_zero",
        "restart_used",
        "convergence",
        "operations",
        "collision_audit",
        "profiles",
        "daughters",
    }
    exact_keys(raw, required, label)
    if raw["schema"] != SUMMARY_SCHEMA:
        raise EvidenceError(f"{label} has an unknown summary schema")
    method = string(raw["method"], f"{label} method")
    if method not in METHODS:
        raise EvidenceError(f"{label} has unknown method {method}")
    result: Dict[str, Any] = {
        "method": method,
        "seed": integer(raw["seed"], f"{label} seed", minimum=1),
        "mpi_ranks": integer(raw["mpi_ranks"], f"{label} mpi_ranks", minimum=1),
        "omp_threads": integer(raw["omp_threads"], f"{label} omp_threads", minimum=1),
    }
    for key in ("physical_input_sha256", "mesh_sha256", "reaction_network_sha256"):
        value = string(raw[key], f"{label} {key}")
        if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
            raise EvidenceError(f"{label} {key} is not a lowercase SHA-256")
        result[key] = value
    if result["mpi_ranks"] != policy["mpi_ranks"] or result["omp_threads"] != policy["omp_threads"]:
        raise EvidenceError(f"{label} MPI/OpenMP decomposition differs from policy")
    if not boolean(raw["started_from_iteration_zero"], f"{label} started_from_iteration_zero"):
        raise EvidenceError(f"{label} did not start from iteration zero")
    if boolean(raw["restart_used"], f"{label} restart_used"):
        raise EvidenceError(f"{label} used a restart without a Stage-9 restart contract")

    convergence = raw["convergence"]
    if not isinstance(convergence, dict):
        raise EvidenceError(f"{label} convergence is not an object")
    exact_keys(
        convergence,
        {"status", "hard_limit_reached", "complete_steady_windows", "discarded_burn_in_iterations"},
        f"{label} convergence",
    )
    result["convergence"] = {
        "status": string(convergence["status"], f"{label} convergence status"),
        "hard_limit_reached": boolean(
            convergence["hard_limit_reached"], f"{label} hard_limit_reached"
        ),
        "complete_steady_windows": integer(
            convergence["complete_steady_windows"], f"{label} complete_steady_windows"
        ),
        "discarded_burn_in_iterations": integer(
            convergence["discarded_burn_in_iterations"],
            f"{label} discarded_burn_in_iterations",
        ),
    }

    operations = raw["operations"]
    if not isinstance(operations, dict):
        raise EvidenceError(f"{label} operations is not an object")
    operation_keys = {
        "wall_seconds",
        "peak_buffer_utilization",
        "buffer_exhaustions",
        "emergency_events",
        "unfinished_trajectories",
    }
    exact_keys(operations, operation_keys, f"{label} operations")
    result["operations"] = {
        "wall_seconds": number(operations["wall_seconds"], f"{label} wall_seconds", minimum=0.0),
        "peak_buffer_utilization": number(
            operations["peak_buffer_utilization"], f"{label} peak_buffer_utilization", minimum=0.0
        ),
        "buffer_exhaustions": integer(operations["buffer_exhaustions"], f"{label} buffer_exhaustions"),
        "emergency_events": integer(operations["emergency_events"], f"{label} emergency_events"),
        "unfinished_trajectories": integer(
            operations["unfinished_trajectories"], f"{label} unfinished_trajectories"
        ),
    }
    if result["operations"]["wall_seconds"] <= 0.0:
        raise EvidenceError(f"{label} wall_seconds must be positive")
    if result["operations"]["peak_buffer_utilization"] > 1.0:
        raise EvidenceError(f"{label} peak_buffer_utilization exceeds one")
    result["operations"]["mpi_rank_seconds"] = (
        result["operations"]["wall_seconds"] * result["mpi_ranks"]
    )

    collision = raw["collision_audit"]
    if not isinstance(collision, dict):
        raise EvidenceError(f"{label} collision_audit is not an object")
    exact_keys(
        collision,
        {"transported_forced_daughter_events", "correction_mismatch_events"},
        f"{label} collision_audit",
    )
    result["collision_audit"] = {
        "transported_forced_daughter_events": integer(
            collision["transported_forced_daughter_events"],
            f"{label} transported_forced_daughter_events",
        ),
        "correction_mismatch_events": integer(
            collision["correction_mismatch_events"], f"{label} correction_mismatch_events"
        ),
    }

    profiles = raw["profiles"]
    expected_species = set(policy["primary_species"]) | set(policy["daughter_species"])
    if not isinstance(profiles, dict):
        raise EvidenceError(f"{label} profiles is not an object")
    exact_keys(profiles, expected_species, f"{label} profiles")
    result["profiles"] = {
        species: _profile(profiles[species], f"{label} {species}")
        for species in sorted(expected_species)
    }

    daughters = raw["daughters"]
    if not isinstance(daughters, dict):
        raise EvidenceError(f"{label} daughters is not an object")
    exact_keys(daughters, policy["daughter_species"], f"{label} daughters")
    result["daughters"] = {
        species: _validate_daughter_record(daughters[species], f"{label} {species}")
        for species in sorted(daughters)
    }
    return result


def _common_profile_grid(runs: Sequence[Dict[str, Any]], species: str) -> List[float]:
    reference = [radius for radius, _value in runs[0]["profiles"][species]]
    for run in runs[1:]:
        radii = [radius for radius, _value in run["profiles"][species]]
        if radii != reference:
            raise EvidenceError(f"{species} profiles do not share an exact radial grid")
    return reference


def _profile_by_method(
    runs: Sequence[Dict[str, Any]], species: str, multiplier: float
) -> Dict[str, List[Dict[str, float]]]:
    radii = _common_profile_grid(runs, species)
    result: Dict[str, List[Dict[str, float]]] = {}
    for method in METHODS:
        selected = [run for run in runs if run["method"] == method]
        rows: List[Dict[str, float]] = []
        for index, radius in enumerate(radii):
            values = [run["profiles"][species][index][1] for run in selected]
            row = {"radius_m": radius, **sample_statistics(values, multiplier)}
            rows.append(row)
        result[method] = rows
    return result


def evaluate(manifest_path: Path) -> Dict[str, Any]:
    manifest_path = manifest_path.resolve()
    manifest = load_object(manifest_path)
    exact_keys(manifest, {"schema", "policy", "runs"}, "manifest")
    if manifest["schema"] != MANIFEST_SCHEMA:
        raise EvidenceError("unknown Stage-9 manifest schema")
    policy = validate_policy(manifest["policy"])
    records = manifest["runs"]
    if not isinstance(records, list):
        raise EvidenceError("manifest runs is not a list")

    minimum_count = 2 * policy["minimum_seeds"]
    if len(records) < minimum_count:
        raise EvidenceError(
            f"manifest needs at least {minimum_count} balanced ANALOG/FORCED records"
        )
    runs: List[Dict[str, Any]] = []
    seen: set[Tuple[str, int]] = set()
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            raise EvidenceError(f"run record {index} is not an object")
        exact_keys(record, {"method", "seed", "summary", "sha256"}, f"run record {index}")
        path = _summary_path(manifest_path, record["summary"])
        expected_hash = string(record["sha256"], f"run record {index} sha256")
        if sha256_file(path) != expected_hash:
            raise EvidenceError(f"summary hash mismatch: {path}")
        summary = validate_summary(load_object(path), policy, str(path))
        if summary["method"] != record["method"] or summary["seed"] != record["seed"]:
            raise EvidenceError(f"run record metadata differs from summary: {path}")
        key = (summary["method"], summary["seed"])
        if key in seen:
            raise EvidenceError(f"duplicate run record {key}")
        seen.add(key)
        summary["path"] = str(path)
        runs.append(summary)

    seeds_by_method = {
        method: sorted(run["seed"] for run in runs if run["method"] == method)
        for method in METHODS
    }
    if seeds_by_method["ANALOG"] != seeds_by_method["FORCED"]:
        raise EvidenceError("ANALOG and FORCED seed sets differ")
    if len(seeds_by_method["ANALOG"]) < policy["minimum_seeds"]:
        raise EvidenceError("independent-seed count is below policy")

    # These hashes describe the physical problem rather than estimator
    # controls.  Requiring exact equality prevents a mesh, chemistry network,
    # or physical source change from masquerading as a variance improvement.
    for key in ("physical_input_sha256", "mesh_sha256", "reaction_network_sha256"):
        values = {run[key] for run in runs}
        if len(values) != 1:
            raise EvidenceError(f"runs do not share one {key}")

    failures: List[str] = []
    for run in runs:
        label = f"{run['method']}/seed={run['seed']}"
        convergence = run["convergence"]
        if convergence["status"] != "converged" or convergence["hard_limit_reached"]:
            failures.append(f"{label}: run is not converged before the hard limit")
        if convergence["complete_steady_windows"] < policy["minimum_steady_windows"]:
            failures.append(f"{label}: too few complete steady windows")
        if convergence["discarded_burn_in_iterations"] < policy["burn_in_iterations"]:
            failures.append(f"{label}: insufficient discarded burn-in")
        operations = run["operations"]
        if operations["peak_buffer_utilization"] > policy["maximum_buffer_utilization"]:
            failures.append(f"{label}: peak particle-buffer utilization exceeds policy")
        if operations["buffer_exhaustions"]:
            failures.append(f"{label}: particle-buffer exhaustion occurred")
        if operations["emergency_events"]:
            failures.append(f"{label}: unplanned emergency adaptation occurred")
        if operations["unfinished_trajectories"]:
            failures.append(f"{label}: unfinished local-timestep trajectory remains")
        collision = run["collision_audit"]
        if collision["correction_mismatch_events"]:
            failures.append(f"{label}: a collision changed a carrier correction")
        if (
            run["method"] == "FORCED"
            and collision["transported_forced_daughter_events"]
            < policy["minimum_forced_daughter_collisions"]
        ):
            failures.append(f"{label}: no qualifying transported forced-daughter collision")
        for species, daughter in run["daughters"].items():
            settings = policy["daughter_species"][species]
            if daughter["ess_fraction"] < settings["minimum_ess_fraction"]:
                failures.append(f"{label}/{species}: ESS fraction is below policy")
            if daughter["zero_cell_fraction"] > settings["maximum_zero_cell_fraction"]:
                failures.append(f"{label}/{species}: zero-cell fraction exceeds policy")

    # Apply exactly the closure rule frozen in the authenticated policy.  The
    # evaluator has no implicit preferred method: CI and per-run relative
    # tolerance are both supported, but the archive must choose one.
    flux_closure: Dict[str, Any] = {}
    flux_policy = policy["flux_closure"]
    for species in policy["daughter_species"]:
        flux_closure[species] = {}
        for method in METHODS:
            values = [
                run["daughters"][species]["signed_flux_fraction"]
                for run in runs
                if run["method"] == method
            ]
            if flux_policy["method"] == "independent_seed_ci_contains_zero":
                row = sample_statistics(values, policy["confidence_multiplier"])
                accepted = row["confidence_low"] <= 0.0 <= row["confidence_high"]
                row.update(
                    policy_method=flux_policy["method"],
                    accepted=accepted,
                )
                failure = (
                    f"{method}/{species}: independent-seed flux-residual CI excludes zero"
                )
            else:
                tolerance = flux_policy["relative_tolerance"]
                absolute_values = [abs(value) for value in values]
                accepted = all(value <= tolerance for value in absolute_values)
                row = {
                    "policy_method": flux_policy["method"],
                    "relative_tolerance": tolerance,
                    "maximum_absolute_residual": max(absolute_values),
                    "accepted_run_count": sum(
                        value <= tolerance for value in absolute_values
                    ),
                    "total_run_count": len(absolute_values),
                    "accepted": accepted,
                }
                failure = (
                    f"{method}/{species}: a run exceeds the flux relative-tolerance policy"
                )
            flux_closure[species][method] = row
            if not accepted:
                failures.append(failure)

    multiplier = policy["confidence_multiplier"]
    profiles: Dict[str, Any] = {}
    all_species = list(policy["primary_species"]) + list(policy["daughter_species"])
    for species in all_species:
        method_rows = _profile_by_method(runs, species, multiplier)
        settings = (
            policy["primary_species"].get(species)
            or policy["daughter_species"][species]
        )
        comparisons: List[Dict[str, float | bool]] = []
        for analog, forced in zip(method_rows["ANALOG"], method_rows["FORCED"]):
            combined_se = math.hypot(
                analog["standard_error"], forced["standard_error"]
            )
            scale = max(abs(analog["mean"]), abs(forced["mean"]), 1.0)
            allowed = max(
                settings["relative_tolerance"] * scale,
                settings["combined_se_multiplier"] * combined_se,
            )
            difference = abs(analog["mean"] - forced["mean"])
            accepted = difference <= allowed
            comparisons.append(
                {
                    "radius_m": analog["radius_m"],
                    "difference": difference,
                    "allowed_difference": allowed,
                    "accepted": accepted,
                }
            )
            if not accepted:
                failures.append(
                    f"{species} profile differs between ANALOG and FORCED at "
                    f"r={analog['radius_m']:.17g} m"
                )
        profiles[species] = {
            "ANALOG": method_rows["ANALOG"],
            "FORCED": method_rows["FORCED"],
            "agreement": comparisons,
        }

    efficiency: Dict[str, Any] = {}
    mean_cost = {
        method: statistics.fmean(
            run["operations"]["mpi_rank_seconds"]
            for run in runs
            if run["method"] == method
        )
        for method in METHODS
    }
    for species, settings in policy["daughter_species"].items():
        rows: List[Dict[str, float | bool]] = []
        accepted_bins = 0
        for analog, forced in zip(
            profiles[species]["ANALOG"], profiles[species]["FORCED"]
        ):
            analog_variance = analog["variance"]
            forced_variance = forced["variance"]
            unbounded_gain = forced_variance == 0.0 and analog_variance > 0.0
            if analog_variance == 0.0:
                gain = 1.0 if forced_variance == 0.0 else 0.0
            elif unbounded_gain:
                # JSON has no representation for infinity.  Preserve the
                # mathematical distinction explicitly instead of emitting a
                # non-standard Infinity token or an arbitrary large number.
                gain = None
            else:
                gain = (
                    analog_variance * mean_cost["ANALOG"]
                    / (forced_variance * mean_cost["FORCED"])
                )
            accepted = unbounded_gain or (
                gain is not None and gain >= settings["minimum_variance_cost_gain"]
            )
            accepted_bins += int(accepted)
            rows.append(
                {
                    "radius_m": analog["radius_m"],
                    "analog_variance": analog_variance,
                    "forced_variance": forced_variance,
                    "variance_cost_gain": gain,
                    "unbounded_gain": unbounded_gain,
                    "accepted": accepted,
                }
            )
        fraction = accepted_bins / len(rows)
        if fraction < settings["minimum_gain_bin_fraction"]:
            failures.append(
                f"{species}: variance-per-cost gain passes only {fraction:.6g} of bins"
            )
        efficiency[species] = {
            "bins": rows,
            "accepted_bin_fraction": fraction,
            "required_bin_fraction": settings["minimum_gain_bin_fraction"],
        }

    return {
        "schema": "amps.comet.stage9-science-result.v1",
        "outcome": "PASS" if not failures else "FAIL_SCIENCE",
        "accepted": not failures,
        "seeds": seeds_by_method["ANALOG"],
        "run_count": len(runs),
        "failures": failures,
        "cost": {
            "metric": "wall_seconds_times_mpi_ranks",
            "mean_mpi_rank_seconds": mean_cost,
        },
        "profiles": profiles,
        "efficiency": efficiency,
        "flux_closure": flux_closure,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate a checksum-protected Stage-9 ANALOG/FORCED seed matrix."
    )
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, help="optional JSON result path")
    arguments = parser.parse_args(argv)
    try:
        result = evaluate(arguments.manifest)
    except EvidenceError as exc:
        result = {
            "schema": "amps.comet.stage9-science-result.v1",
            "outcome": "INVALID_EVIDENCE",
            "accepted": False,
            "error": str(exc),
        }
        exit_code = 12
    else:
        exit_code = 0 if result["accepted"] else 10
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    sys.stdout.write(rendered)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

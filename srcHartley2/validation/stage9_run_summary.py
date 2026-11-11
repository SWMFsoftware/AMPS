#!/usr/bin/env python3
"""Build one strict ``amps.comet.stage9-run-summary.v1`` record.

The builder is intentionally a verifier, not a data-repair utility.  It reads
the completed Stage-0 chemistry reports, Stage-9 checkpoint/carrier tables,
and AMPS plot profiles from one fresh campaign attempt.  Missing fields,
inconsistent provenance, ambiguous profile files, non-closing carrier sums,
or absent collision evidence terminate the build.  No science value is
interpolated or inferred from a successful process return code.

``run_stage9_science_campaign.py`` imports :func:`build_summary` after AMPS
returns zero.  The standalone CLI is useful for auditing a preserved attempt;
it never modifies the run's raw diagnostic files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple


SUMMARY_SCHEMA = "amps.comet.stage9-run-summary.v1"
CHECKPOINT_SCHEMA = "amps.comet.stage9-checkpoint.v1"


class SummaryError(RuntimeError):
    """A run lacks complete, internally consistent Stage-9 evidence."""


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
        raise SummaryError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SummaryError(f"JSON root is not an object: {path}")
    return value


def atomic_write_json(path: Path, value: object) -> None:
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


def number(value: object, label: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SummaryError(f"{label} is not numeric")
    result = float(value)
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        raise SummaryError(f"{label} is not a valid finite value")
    return result


def integer(value: object, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise SummaryError(f"{label} is not an integer >= {minimum}")
    return value


def exact_keys(value: Mapping[str, Any], expected: Iterable[str], label: str) -> None:
    actual = set(value)
    wanted = set(expected)
    if actual != wanted:
        raise SummaryError(
            f"{label} keys differ: missing={sorted(wanted-actual)} extra={sorted(actual-wanted)}"
        )


def close(left: float, right: float, tolerance: float = 5.0e-12) -> bool:
    return abs(left - right) <= tolerance * max(abs(left), abs(right), 1.0)


def normalize_mode(value: object, label: str) -> str:
    if not isinstance(value, str):
        raise SummaryError(f"{label} estimator mode is not text")
    mode = value.strip().upper()
    aliases = {"ANALOG": "ANALOG", "ENHANCED": "FORCED", "FORCED": "FORCED"}
    try:
        return aliases[mode]
    except KeyError as exc:
        raise SummaryError(f"{label} has unknown estimator mode {value!r}") from exc


def _indexed_files(directory: Path, pattern: str) -> Dict[int, Path]:
    regex = re.compile(pattern)
    result: Dict[int, Path] = {}
    for path in directory.iterdir() if directory.is_dir() else ():
        match = regex.fullmatch(path.name)
        if match is None:
            continue
        index = int(match.group(1))
        if index in result:
            raise SummaryError(f"duplicate output index {index} in {directory}")
        result[index] = path
    return result


def _read_checkpoint(
    path: Path, output_index: int, method: str, seed: int, ranks: int
) -> Dict[str, Any]:
    value = load_object(path)
    required = {
        "schema",
        "iteration",
        "output_index",
        "sampling_length",
        "estimator_mode",
        "random_seed",
        "mpi_ranks",
        "window_wall_seconds",
        "peak_rank_live_particles",
        "total_live_particles",
        "rank_capacity_min",
        "rank_capacity_max",
        "peak_buffer_utilization",
        "safe_boundary_audits_min",
        "safe_boundary_audits_max",
        "pending_remainders",
        "pending_reaction_events",
        "species",
    }
    exact_keys(value, required, str(path))
    if value["schema"] != CHECKPOINT_SCHEMA:
        raise SummaryError(f"unknown checkpoint schema: {path}")
    if integer(value["output_index"], f"{path} output_index") != output_index:
        raise SummaryError(f"checkpoint index disagrees with filename: {path}")
    if normalize_mode(value["estimator_mode"], str(path)) != method:
        raise SummaryError(f"checkpoint estimator mode mismatch: {path}")
    if integer(value["random_seed"], f"{path} random_seed", minimum=1) != seed:
        raise SummaryError(f"checkpoint seed mismatch: {path}")
    if integer(value["mpi_ranks"], f"{path} mpi_ranks", minimum=1) != ranks:
        raise SummaryError(f"checkpoint MPI size mismatch: {path}")
    sampling = integer(value["sampling_length"], f"{path} sampling_length", minimum=1)
    iteration = integer(value["iteration"], f"{path} iteration", minimum=1)
    if iteration != (output_index + 1) * sampling:
        raise SummaryError(f"checkpoint iteration/index/sampling identity fails: {path}")
    for key in (
        "window_wall_seconds",
        "peak_rank_live_particles",
        "total_live_particles",
        "rank_capacity_min",
        "rank_capacity_max",
        "peak_buffer_utilization",
        "safe_boundary_audits_min",
        "safe_boundary_audits_max",
        "pending_remainders",
        "pending_reaction_events",
    ):
        number(value[key], f"{path} {key}", minimum=0.0)
    if value["peak_buffer_utilization"] > 1.0:
        raise SummaryError(f"checkpoint buffer utilization exceeds one: {path}")
    if value["pending_remainders"] != 0 or value["pending_reaction_events"] != 0:
        raise SummaryError(f"checkpoint contains unfinished trajectories: {path}")
    if value["safe_boundary_audits_min"] <= 0:
        raise SummaryError(f"checkpoint has no safe-boundary audit: {path}")
    species = value["species"]
    if not isinstance(species, list) or not species:
        raise SummaryError(f"checkpoint has no species records: {path}")
    by_symbol: Dict[str, Dict[str, Any]] = {}
    for row in species:
        if not isinstance(row, dict) or not isinstance(row.get("symbol"), str):
            raise SummaryError(f"checkpoint species record is malformed: {path}")
        symbol = row["symbol"]
        if symbol in by_symbol:
            raise SummaryError(f"duplicate checkpoint species {symbol}: {path}")
        by_symbol[symbol] = row
    value["species_by_symbol"] = by_symbol
    return value


def _read_carrier_cells(
    path: Path, output: int, method: str, seed: int, ranks: int
) -> Dict[str, float]:
    required = {
        "cell_index",
        "model_carriers",
        "sum_weights",
        "sum_squared_weights",
        "carrier_ess",
        "correction_min",
        "correction_max",
        "output_index",
        "sampling_length",
        "estimator_mode",
        "random_seed",
        "mpi_ranks",
    }
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or set(reader.fieldnames) != required:
                raise SummaryError(f"carrier cell schema differs: {path}")
            rows = list(reader)
    except OSError as exc:
        raise SummaryError(f"cannot read carrier table {path}: {exc}") from exc
    if not rows:
        raise SummaryError(f"carrier table is empty: {path}")
    count = 0
    sum_weights = 0.0
    sum_squares = 0.0
    occupied = 0
    for index, row in enumerate(rows):
        try:
            if int(row["cell_index"]) != index:
                raise SummaryError(f"noncontiguous carrier cell index: {path}")
            if int(row["output_index"]) != output or int(row["random_seed"]) != seed:
                raise SummaryError(f"carrier provenance mismatch: {path}")
            if int(row["mpi_ranks"]) != ranks:
                raise SummaryError(f"carrier MPI provenance mismatch: {path}")
            if normalize_mode(row["estimator_mode"], str(path)) != method:
                raise SummaryError(f"carrier estimator mode mismatch: {path}")
            carriers = int(row["model_carriers"])
            weight = float(row["sum_weights"])
            weight2 = float(row["sum_squared_weights"])
            ess = float(row["carrier_ess"])
        except (KeyError, ValueError) as exc:
            raise SummaryError(f"malformed carrier value: {path}") from exc
        if carriers < 0 or min(weight, weight2, ess) < 0.0 or not all(
            math.isfinite(item) for item in (weight, weight2, ess)
        ):
            raise SummaryError(f"invalid carrier value: {path}")
        expected_ess = weight * weight / weight2 if weight2 > 0.0 else 0.0
        if not close(ess, expected_ess):
            raise SummaryError(f"carrier ESS identity fails: {path}")
        if ess > carriers * (1.0 + 5.0e-12):
            raise SummaryError(f"carrier ESS exceeds count: {path}")
        if carriers == 0:
            if weight != 0.0 or weight2 != 0.0 or row["correction_min"] or row["correction_max"]:
                raise SummaryError(f"empty carrier cell has nonempty moments: {path}")
        else:
            occupied += 1
            if not row["correction_min"] or not row["correction_max"]:
                raise SummaryError(f"occupied cell lacks correction bounds: {path}")
        count += carriers
        sum_weights += weight
        sum_squares += weight2
    return {
        "cell_count": float(len(rows)),
        "occupied_cells": float(occupied),
        "model_carriers": float(count),
        "sum_weights": sum_weights,
        "sum_squared_weights": sum_squares,
    }


def _parse_variables(line: str) -> List[str]:
    return re.findall(r'"([^"]+)"', line)


def _read_plot_profile(path: Path) -> List[Tuple[float, float]]:
    try:
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    except OSError as exc:
        raise SummaryError(f"cannot read plot profile {path}: {exc}") from exc
    variable_line = next((line for line in lines if line.lstrip().upper().startswith("VARIABLES")), None)
    if variable_line is None:
        raise SummaryError(f"plot profile lacks VARIABLES: {path}")
    variables = _parse_variables(variable_line)
    try:
        radius_column = variables.index("X")
        density_column = variables.index("Number Density[1/m^3]")
    except ValueError as exc:
        raise SummaryError(f"plot profile lacks radius/density columns: {path}") from exc
    rows: List[Tuple[float, float]] = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped[0].isalpha() or stripped[0] in {'#', '"'}:
            continue
        try:
            values = [float(field) for field in stripped.replace(",", " ").split()]
        except ValueError:
            continue
        if len(values) != len(variables):
            raise SummaryError(f"plot data width differs from VARIABLES: {path}")
        radius, density = values[radius_column], values[density_column]
        if not math.isfinite(radius) or radius < 0.0 or not math.isfinite(density) or density < 0.0:
            raise SummaryError(f"plot profile has invalid radius/density: {path}")
        rows.append((radius, density))
    if not rows:
        raise SummaryError(f"plot profile has no numeric data: {path}")

    # AMPS can repeat a block-face node in adjacent Tecplot zones.  Such rows
    # describe the same physical point and may be collapsed only if their
    # sampled density agrees to roundoff; averaging disagreeing duplicates
    # would hide a decomposition/output defect.
    result: List[Tuple[float, float]] = []
    for radius, density in rows:
        if result and radius < result[-1][0]:
            raise SummaryError(f"plot profile is not radially ordered: {path}")
        if result and radius == result[-1][0]:
            if not close(density, result[-1][1]):
                raise SummaryError(f"duplicate plot radius has conflicting density: {path}")
            continue
        result.append((radius, density))
    return result


def _find_plot(plots: Path, species: str, output: int) -> Path:
    matches = sorted(plots.glob(f"pic.{species}.s=*.out={output}.dat"))
    if len(matches) != 1:
        raise SummaryError(
            f"expected one plot profile for {species} out={output}, found {len(matches)}"
        )
    return matches[0]


def _mean_profiles(profiles: Sequence[List[Tuple[float, float]]], label: str) -> List[Dict[str, float]]:
    reference = [radius for radius, _value in profiles[0]]
    for profile in profiles[1:]:
        if [radius for radius, _value in profile] != reference:
            raise SummaryError(f"{label} profile grid changes between accepted windows")
    return [
        {
            "radius_m": radius,
            "value": sum(profile[index][1] for profile in profiles) / len(profiles),
        }
        for index, radius in enumerate(reference)
    ]


def _mean(values: Sequence[float], label: str) -> float:
    if not values:
        raise SummaryError(f"cannot average empty {label}")
    return sum(values) / len(values)


def _flux_entry(report: Dict[str, Any], species: str, path: Path) -> Dict[str, Any]:
    rows = report.get("species_flux")
    if not isinstance(rows, list):
        raise SummaryError(f"Stage-0 report has no species_flux: {path}")
    matches = [row for row in rows if isinstance(row, dict) and row.get("species") == species]
    if len(matches) != 1:
        raise SummaryError(f"Stage-0 report lacks unique {species} flux: {path}")
    return matches[0]


def _collision_audit(
    attempt: Path, method: str, seed: int, ranks: int
) -> Tuple[Dict[str, int], List[Path]]:
    """Authenticate the focused production-NTC collision prerequisite.

    Adding a per-collision observer to the general NTC hot path would change
    the cost and risk profile of every long science run.  The campaign runner
    instead executes the frozen executable's focused transported-daughter
    fixture immediately before each ENHANCED case.  This parser accepts only
    its complete PASS record and independently checks the conservation and
    ownership fields needed by Stage 9.  ANALOG has no forced daughters by
    definition and therefore reports an exact zero without a fixture.
    """
    if method == "ANALOG":
        return {
            "transported_forced_daughter_events": 0,
            "correction_mismatch_events": 0,
        }, []
    root = attempt / "collision-preflight"
    run_log = root / "run.log"
    error_log = root / "error.log"
    if not run_log.is_file() or not error_log.is_file():
        raise SummaryError("ENHANCED attempt lacks focused collision-preflight logs")
    text = run_log.read_text(encoding="utf-8", errors="replace")
    if f"comet random seed={seed}" not in text:
        raise SummaryError("focused collision preflight does not prove the requested seed")
    if f"The total number of the MPI processes={ranks}" not in text:
        raise SummaryError("focused collision preflight does not prove its MPI size")
    records = [
        line.strip()
        for line in text.splitlines()
        if " transported-daughter-collision " in line
    ]
    if len(records) != 1:
        raise SummaryError("focused collision preflight lacks one transported-daughter record")
    fields = dict(re.findall(r"([a-z][a-z0-9-]*)=([^ ]+)", records[0]))
    required = {
        "ranks",
        "reaction",
        "daughters",
        "partners",
        "daughter-velocity-changes",
        "partner-velocity-changes",
        "unfinished-trajectories",
        "wrong-owner",
        "cycles",
        "trials",
        "accepted",
        "updates",
        "invalid-probabilities",
        "roundoff-clamps",
        "correction-error",
        "represented-rate-error",
        "momentum-error",
    }
    if not records[0].startswith("PASS ") or not required.issubset(fields):
        raise SummaryError("focused collision preflight did not emit its complete PASS schema")
    try:
        integers = {
            key: int(fields[key])
            for key in required
            if key
            not in {
                "reaction",
                "correction-error",
                "represented-rate-error",
                "momentum-error",
            }
        }
        correction_error = float(fields["correction-error"])
        rate_error = float(fields["represented-rate-error"])
        momentum_error = float(fields["momentum-error"])
    except ValueError as exc:
        raise SummaryError("focused collision preflight contains malformed values") from exc
    if (
        integers["ranks"] != ranks
        or fields["reaction"] != "C2_2C"
        or integers["daughters"] != 16
        or integers["partners"] != 64
        or integers["daughter-velocity-changes"] < 8
        or integers["partner-velocity-changes"] < 8
        or integers["unfinished-trajectories"] != 0
        or integers["wrong-owner"] != 0
        or integers["cycles"] <= 0
        or integers["trials"] <= 0
        or integers["accepted"] <= 0
        or integers["updates"] < 64
        or integers["invalid-probabilities"] != 0
        or integers["roundoff-clamps"] != 0
        or correction_error != 0.0
        or rate_error != 0.0
        or not math.isfinite(momentum_error)
        or momentum_error > 2.0e-12
    ):
        raise SummaryError("focused collision preflight violates Stage-9 invariants")
    return {
        "transported_forced_daughter_events": integers["accepted"],
        "correction_mismatch_events": 0,
    }, [run_log, error_log]


def build_summary(archive: Path, attempt: Path) -> Tuple[Dict[str, Any], List[Path]]:
    """Validate an attempt and return its compact summary plus raw artifacts."""
    archive = archive.resolve()
    attempt = attempt.resolve()
    try:
        attempt.relative_to(archive)
    except ValueError as exc:
        raise SummaryError("attempt is outside the Stage-9 archive") from exc
    campaign = load_object(archive / "CAMPAIGN_STATE.json")
    policy = load_object(archive / "frozen" / "stage9_policy.json")
    status = load_object(attempt / "status.json")
    if campaign.get("schema") != "amps.comet.stage9-campaign-state.v1":
        raise SummaryError("unknown Stage-9 campaign-state schema")
    if policy.get("schema") != "amps.comet.stage9-campaign-policy.v1":
        raise SummaryError("unknown Stage-9 campaign policy schema")
    if status.get("state") not in {"summarizing", "completed"} or status.get("return_code") != 0:
        raise SummaryError("attempt is not a successful completed process")
    arm = status.get("estimator_mode")
    method = normalize_mode(arm, str(attempt / "status.json"))
    seed = integer(status.get("seed"), "attempt seed", minimum=1)
    ranks = integer(policy.get("mpi_ranks"), "policy mpi_ranks", minimum=1)
    if ranks != 12 or policy.get("omp_threads") != 1:
        raise SummaryError("Stage-9 campaign must use 12 MPI ranks and one thread")

    diagnostics = attempt / "PT" / "diagnostics"
    plots = attempt / "PT" / "plots"
    checkpoints = _indexed_files(
        diagnostics, r"comet_stage9_checkpoint\.out=(\d+)\.json"
    )
    reports = _indexed_files(diagnostics, r"comet_stage0\.out=(\d+)\.json")
    if not checkpoints or sorted(checkpoints) != sorted(reports):
        raise SummaryError("Stage-0 and Stage-9 checkpoint index sets differ or are empty")
    checkpoint_values = {
        index: _read_checkpoint(path, index, method, seed, ranks)
        for index, path in checkpoints.items()
    }
    sampling_lengths = {value["sampling_length"] for value in checkpoint_values.values()}
    if len(sampling_lengths) != 1:
        raise SummaryError("sampling length changes between Stage-9 checkpoints")
    sampling_length = int(next(iter(sampling_lengths)))
    if sampling_length != policy.get("sampling_length"):
        raise SummaryError("checkpoint sampling length differs from policy")
    burn_in = integer(policy.get("burn_in_iterations"), "policy burn_in_iterations")
    accepted = [
        index
        for index in sorted(checkpoints)
        if checkpoint_values[index]["iteration"] > burn_in
    ]
    minimum_windows = integer(
        policy.get("minimum_steady_windows"), "policy minimum_steady_windows", minimum=2
    )
    if len(accepted) < minimum_windows:
        raise SummaryError("too few complete post-burn-in Stage-9 windows")

    report_values: Dict[int, Dict[str, Any]] = {}
    artifacts: List[Path] = [attempt / "run.log", attempt / "error.log"]
    for index in sorted(reports):
        report = load_object(reports[index])
        if report.get("schema") != "amps.comet.stage0-diagnostics.v1":
            raise SummaryError(f"unknown Stage-0 diagnostic schema: {reports[index]}")
        if (
            report.get("output_index") != index
            or report.get("iteration") != checkpoint_values[index]["iteration"]
            or report.get("sampling_length") != sampling_length
            or report.get("random_seed") != seed
            or report.get("random_seed_configured") is not True
            or report.get("mpi_ranks") != ranks
        ):
            raise SummaryError(f"Stage-0 provenance mismatch: {reports[index]}")
        report_values[index] = report
        artifacts.extend(
            [
                reports[index],
                diagnostics / f"comet_stage0.out={index}.csv",
                checkpoints[index],
                diagnostics / f"comet_stage9_checkpoint.out={index}.csv",
            ]
        )
    if any(not path.is_file() for path in artifacts):
        raise SummaryError("a required paired CSV/log artifact is missing")

    gate_policy = policy.get("gate_policy")
    if not isinstance(gate_policy, dict):
        raise SummaryError("campaign policy has no gate_policy")
    primary = gate_policy.get("primary_species")
    daughters = gate_policy.get("daughter_species")
    if not isinstance(primary, dict) or not isinstance(daughters, dict) or not primary or not daughters:
        raise SummaryError("campaign gate_policy has no primary/daughter species")
    requested_species = list(primary) + list(daughters)
    profiles: Dict[str, List[Dict[str, float]]] = {}
    daughter_summary: Dict[str, Dict[str, float]] = {}

    for species in requested_species:
        plot_profiles: List[List[Tuple[float, float]]] = []
        aggregate_count = 0.0
        aggregate_weight = 0.0
        aggregate_weight2 = 0.0
        zero_cells = 0
        total_cells = 0
        for index in accepted:
            plot_path = _find_plot(plots, species, index)
            profile = _read_plot_profile(plot_path)
            plot_profiles.append(profile)
            artifacts.append(plot_path)
            zero_cells += sum(value == 0.0 for _radius, value in profile)
            total_cells += len(profile)
            carrier_path = diagnostics / f"stage9_carriers.{species}.out={index}.csv"
            cells = _read_carrier_cells(carrier_path, index, method, seed, ranks)
            artifacts.append(carrier_path)
            checkpoint_species = checkpoint_values[index]["species_by_symbol"].get(species)
            if checkpoint_species is None:
                raise SummaryError(f"checkpoint lacks requested species {species}")
            for key, checkpoint_key in (
                ("model_carriers", "model_carriers"),
                ("sum_weights", "physical_inventory"),
                ("sum_squared_weights", "sum_squared_weights"),
            ):
                expected = number(
                    checkpoint_species.get(checkpoint_key),
                    f"checkpoint {species} {checkpoint_key}",
                    minimum=0.0,
                )
                if not close(cells[key], expected):
                    raise SummaryError(
                        f"carrier-cell and checkpoint {key} differ for {species} out={index}"
                    )
            aggregate_count += cells["model_carriers"]
            aggregate_weight += cells["sum_weights"]
            aggregate_weight2 += cells["sum_squared_weights"]
        profiles[species] = _mean_profiles(plot_profiles, species)
        if species not in daughters:
            continue
        flux_rows = [
            _flux_entry(report_values[index], species, reports[index])
            for index in accepted
        ]
        flux_mapping = {
            "surface_source_s-1": "sampled_surface_source_s-1",
            "chemical_production_s-1": "chemical_represented_production_s-1",
            "inner_reemitted_s-1": "inner_reemitted_s-1",
            "chemical_destruction_s-1": "chemical_destruction_s-1",
            "inner_incident_s-1": "inner_incident_s-1",
            "outer_escape_s-1": "outer_escape_s-1",
            "net_flux_s-1": "net_represented_flux_s-1",
        }
        daughter_summary[species] = {
            "zero_cell_fraction": zero_cells / total_cells,
            "model_carriers": aggregate_count,
            "sum_weights": aggregate_weight,
            "sum_squared_weights": aggregate_weight2,
            **{
                destination: _mean(
                    [
                        number(
                            row.get(source),
                            f"Stage-0 {species} {source}",
                            minimum=None if destination == "net_flux_s-1" else 0.0,
                        )
                        for row in flux_rows
                    ],
                    f"{species} {source}",
                )
                for destination, source in flux_mapping.items()
            },
        }

    collision, collision_artifacts = _collision_audit(attempt, method, seed, ranks)
    artifacts.extend(collision_artifacts)

    run_log = (attempt / "run.log").read_text(encoding="utf-8", errors="replace")
    converged = re.search(
        r"radial simulation converged after (\d+) iterations and (\d+) complete sampling windows",
        run_log,
    )
    hard_limit = "radial simulation reached the hard limit" in run_log
    convergence_status = "converged" if converged is not None and not hard_limit else "nonconverged"
    if converged is not None and int(converged.group(1)) != checkpoint_values[max(checkpoints)]["iteration"]:
        raise SummaryError("final convergence iteration differs from final checkpoint")
    emergency_total = 0
    for index in sorted(report_values):
        events = report_values[index].get("emergency_events")
        if not isinstance(events, dict) or not isinstance(events.get("total"), int):
            raise SummaryError(f"Stage-0 report lacks emergency counters: {reports[index]}")
        emergency_total += events["total"]

    provenance = campaign.get("frozen")
    if not isinstance(provenance, dict):
        raise SummaryError("campaign state lacks frozen provenance")
    input_hash = provenance.get("input", {}).get("sha256")
    chemistry_hash = provenance.get("headers", {}).get(
        "CometChemistry.generated.h", {}
    ).get("sha256")
    if not isinstance(input_hash, str) or not isinstance(chemistry_hash, str):
        raise SummaryError("campaign state lacks input/reaction-network hashes")
    mesh_path = attempt / "comet.radial.mesh.msh"
    if not mesh_path.is_file():
        raise SummaryError("run has no radial mesh artifact")
    artifacts.append(mesh_path)
    pending = sum(
        int(checkpoint_values[index]["pending_remainders"])
        + int(checkpoint_values[index]["pending_reaction_events"])
        for index in checkpoint_values
    )
    peak_buffer = max(
        float(checkpoint_values[index]["peak_buffer_utilization"])
        for index in checkpoint_values
    )
    wall_seconds = number(status.get("wall_seconds"), "attempt wall_seconds", minimum=0.0)
    if wall_seconds <= 0.0:
        raise SummaryError("attempt wall time is not positive")

    summary = {
        "schema": SUMMARY_SCHEMA,
        "method": method,
        "seed": seed,
        "mpi_ranks": ranks,
        "omp_threads": 1,
        "physical_input_sha256": input_hash,
        "mesh_sha256": sha256_file(mesh_path),
        "reaction_network_sha256": chemistry_hash,
        "started_from_iteration_zero": status.get("fresh_start_verified") is True,
        "restart_used": status.get("restart_used") is True,
        "convergence": {
            "status": convergence_status,
            "hard_limit_reached": hard_limit,
            "complete_steady_windows": len(accepted),
            "discarded_burn_in_iterations": burn_in,
        },
        "operations": {
            "wall_seconds": wall_seconds,
            "peak_buffer_utilization": peak_buffer,
            # Buffer exhaustion terminates AMPS.  This zero is admitted only
            # after a zero process status and a continuously sampled peak
            # strictly below one; otherwise summary construction fails.
            "buffer_exhaustions": 0,
            "emergency_events": emergency_total,
            "unfinished_trajectories": pending,
        },
        "collision_audit": collision,
        "profiles": profiles,
        "daughters": daughter_summary,
    }
    if peak_buffer >= 1.0:
        raise SummaryError("particle buffer reached full utilization")
    return summary, sorted(set(artifacts))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--attempt", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args(argv)
    try:
        summary, _artifacts = build_summary(arguments.archive, arguments.attempt)
    except SummaryError as exc:
        print(f"stage9_run_summary.py: {exc}", file=os.sys.stderr)
        return 12
    output = arguments.output or arguments.attempt / "stage9_summary.json"
    atomic_write_json(output, summary)
    print(f"PASS Stage-9 run summary output={output} sha256={sha256_file(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Evaluate the permanent eight-run Stage-3 comet science gate.

The evaluator is deliberately read-only with respect to AMPS run directories.
It writes compact JSON/CSV/Markdown summaries under ``evaluation/`` and never
repairs, interpolates, or silently drops invalid evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
sys.path.insert(0, str(THIS_DIR))
sys.path.insert(0, str(THIS_DIR.parent / "visualization"))

from adaptive_radial_bins import Cell, load_cells_with_identity  # noqa: E402
from plot_1d import PlotError  # noqa: E402
from stage3_gate_common import (  # noqa: E402
    GateError,
    MANIFEST_SCHEMA,
    atomic_write_json,
    load_policy,
    sha256_file,
    sha256_json,
)


OUTCOME_EXIT = {
    "PASS_SKIP_STAGE4": 0,
    "FAIL_CONSIDER_STAGE4": 10,
    "INCONCLUSIVE_RERUN": 11,
    "INVALID_EVIDENCE": 12,
}


def _load_json(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GateError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise GateError(f"JSON root is not an object: {path}")
    return value


def _archive_file(archive: Path, relative: object, label: str) -> Path:
    """Resolve one manifest path without permitting archive traversal."""
    if not isinstance(relative, str) or not relative:
        raise GateError(f"manifest is missing {label}")
    candidate = (archive / relative).resolve()
    try:
        candidate.relative_to(archive.resolve())
    except ValueError as exc:
        raise GateError(f"manifest {label} escapes the archive") from exc
    if not candidate.is_file():
        raise GateError(f"manifest {label} is missing: {candidate}")
    return candidate


def _verify_artifact_manifest(
    run_root: Path, run: Dict[str, Any], policy: Dict[str, Any]
) -> None:
    """Re-hash the production evidence instead of trusting runner state.

    The operational runner performs this check before it declares a case
    complete, but a table can be damaged or edited between that operation and
    scientific evaluation.  The evaluator is independently responsible for
    rejecting such evidence; embedded seed/mode metadata is not a checksum.
    """
    checksum_path = run_root / "artifacts.sha256"
    expected_manifest_hash = run.get("artifact_manifest_sha256")
    if (
        not checksum_path.is_file()
        or not isinstance(expected_manifest_hash, str)
        or sha256_file(checksum_path) != expected_manifest_hash
    ):
        raise GateError(f"{run_root}: artifact manifest is missing or changed")

    records: Dict[str, str] = {}
    try:
        lines = checksum_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise GateError(f"cannot read artifact manifest {checksum_path}: {exc}") from exc
    for line in lines:
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            raise GateError(f"{checksum_path}: malformed checksum record")
        expected_hash, relative = match.groups()
        if relative in records:
            raise GateError(f"{checksum_path}: duplicate artifact {relative}")
        path = (run_root / relative).resolve()
        try:
            path.relative_to(run_root.resolve())
        except ValueError as exc:
            raise GateError(f"{checksum_path}: artifact escapes attempt: {relative}") from exc
        if not path.is_file() or sha256_file(path) != expected_hash:
            raise GateError(f"{run_root}: artifact changed: {relative}")
        records[relative] = expected_hash

    diagnostics = Path("PT") / "diagnostics"
    required = {
        "run.log",
        "error.log",
        str(diagnostics / "cell_measures.csv"),
    }
    output_count = policy["required_converged_windows"] + 1
    for output in range(output_count):
        required.add(str(diagnostics / f"comet_stage0.out={output}.json"))
        required.add(str(diagnostics / f"comet_stage0.out={output}.csv"))
        for species in policy["important_species"]:
            required.add(str(diagnostics / f"stage3_cells.{species}.out={output}.csv"))
    missing = sorted(required.difference(records))
    if missing:
        raise GateError(
            f"{checksum_path}: required artifact is not protected: {missing[0]}"
        )


def _manifest_runs(archive: Path, policy: Dict[str, Any]) -> List[Dict[str, Any]]:
    manifest = _load_json(archive / "RUN_MANIFEST.json")
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise GateError("unknown science-run manifest schema")
    if manifest.get("status") != "completed":
        raise GateError("science-run manifest is not operationally complete")
    if manifest.get("policy_sha256") != sha256_json(policy):
        raise GateError("manifest policy hash does not match the selected policy")
    frozen_policy = _archive_file(
        archive, manifest.get("policy_path"), "frozen policy"
    )
    if sha256_json(load_policy(frozen_policy)) != sha256_json(policy):
        raise GateError("frozen archive policy differs from the selected policy")
    for label, path_key, hash_key in (
        ("executable", "frozen_executable", "frozen_executable_sha256"),
        ("input", "frozen_input", "frozen_input_sha256"),
        ("picGlobal.dfn", "frozen_pic_global", "frozen_pic_global_sha256"),
    ):
        path = _archive_file(archive, manifest.get(path_key), f"frozen {label}")
        expected_hash = manifest.get(hash_key)
        if not isinstance(expected_hash, str) or sha256_file(path) != expected_hash:
            raise GateError(f"frozen {label} differs from the manifest")
    runs = manifest.get("runs")
    if not isinstance(runs, list) or len(runs) != 8:
        raise GateError("manifest must contain exactly eight run records")
    expected = {(mode, seed) for mode in policy["modes"] for seed in policy["seeds"]}
    actual = {(run.get("mode"), run.get("seed")) for run in runs}
    if actual != expected or len(actual) != 8:
        raise GateError("manifest mode/seed matrix differs from the policy")
    binary_hash = manifest.get("frozen_executable_sha256")
    input_hash = manifest.get("frozen_input_sha256")
    pic_global_hash = manifest.get("frozen_pic_global_sha256")
    for run in runs:
        if run.get("status") != "completed" or run.get("return_code") != 0:
            raise GateError(
                f"run {run.get('mode')}/{run.get('seed')} is not operationally complete"
            )
        if run.get("mpi_ranks") != policy["mpi_ranks"]:
            raise GateError("run manifest contains the wrong MPI rank count")
        if run.get("executable_sha256") != binary_hash:
            raise GateError("run manifest mixes executable hashes")
        if run.get("input_sha256") != input_hash:
            raise GateError("run manifest mixes input hashes")
        if run.get("pic_global_sha256") != pic_global_hash:
            raise GateError("run manifest mixes picGlobal.dfn hashes")
        relative = run.get("attempt_path")
        if not isinstance(relative, str):
            raise GateError("run manifest is missing attempt_path")
        root = (archive / relative).resolve()
        try:
            root.relative_to(archive.resolve())
        except ValueError as exc:
            raise GateError("run attempt escapes the archive") from exc
        _verify_artifact_manifest(root, run, policy)
        run["_root"] = root
    return runs


def _load_complete_windows(
    run: Dict[str, Any], policy: Dict[str, Any]
) -> Tuple[Dict[int, Dict[str, Any]], List[str]]:
    diagnostics = run["_root"] / "PT" / "diagnostics"
    paths = sorted(diagnostics.glob("comet_stage0.out=*.json"))
    if not paths:
        raise GateError(f"{run['_root']}: no complete diagnostic windows")
    windows: Dict[int, Dict[str, Any]] = {}
    for path in paths:
        report = _load_json(path)
        index = report.get("output_index")
        if not isinstance(index, int) or index in windows:
            raise GateError(f"{path}: invalid or duplicate output index")
        windows[index] = report
    expected_indices = list(range(policy["required_converged_windows"] + 1))
    if sorted(windows) != expected_indices:
        raise GateError(
            f"{run['_root']}: complete outputs must be contiguous {expected_indices}"
        )

    inconclusive: List[str] = []
    for index in expected_indices:
        report = windows[index]
        expected_iteration = (index + 1) * policy["sampling_length"]
        if report.get("iteration") != expected_iteration:
            raise GateError(f"{run['_root']}: output {index} iteration mismatch")
        if report.get("sampling_length") != policy["sampling_length"]:
            raise GateError(f"{run['_root']}: partial or wrong-length output {index}")
        if report.get("random_seed") != run["seed"] or report.get(
            "random_seed_configured"
        ) is not True:
            raise GateError(f"{run['_root']}: seed provenance mismatch")
        if report.get("mpi_ranks") != policy["mpi_ranks"]:
            raise GateError(f"{run['_root']}: MPI provenance mismatch")
        if report.get("residual_sampling") != run["mode"]:
            raise GateError(f"{run['_root']}: residual-mode provenance mismatch")

    accepted = expected_indices[1:]
    for index in accepted:
        report = windows[index]
        events = report.get("emergency_events")
        if not isinstance(events, dict) or set(events) != {
            "total",
            "buffer_utilization",
            "sampling_population",
            "allocation_burst",
        }:
            raise GateError(f"{run['_root']}: output {index} lacks emergency evidence")
        if any(not isinstance(value, int) or value < 0 for value in events.values()):
            raise GateError(f"{run['_root']}: invalid emergency counters")
        if policy["require_zero_emergency_events_in_accepted_windows"] and any(
            events.values()
        ):
            inconclusive.append(f"output {index} executed emergency adaptation")
        requests = report.get("emergency_requests")
        if not isinstance(requests, int) or requests < 0:
            raise GateError(f"{run['_root']}: invalid emergency-request counter")
        if policy["require_zero_emergency_events_in_accepted_windows"] and requests:
            inconclusive.append(f"output {index} requested emergency adaptation")

        flux = report.get("species_flux")
        reactions = report.get("reactions")
        if not isinstance(flux, list) or not isinstance(reactions, list):
            raise GateError(f"{run['_root']}: output {index} has incomplete chemistry")
        changes: List[float] = []
        species_seen = set()
        for entry in flux:
            species = entry.get("species")
            species_seen.add(species)
            change = entry.get("relative_l1_change")
            if not isinstance(change, (int, float)) or not math.isfinite(change):
                raise GateError(f"{run['_root']}: output {index} lacks convergence data")
            changes.append(float(change))
            if entry.get("production_agreement") is not True:
                inconclusive.append(
                    f"output {index} species {species} production closure failed"
                )
            configured = float(entry.get("configured_surface_source_s-1", 0.0))
            sampled = float(entry.get("sampled_surface_source_s-1", 0.0))
            if configured > 0.0:
                relative = abs(sampled - configured) / configured
                if relative > policy["direct_source_relative_tolerance"]:
                    inconclusive.append(
                        f"output {index} species {species} direct source differs by {relative:.6g}"
                    )
        if not set(policy["important_species"]).issubset(species_seen):
            raise GateError(f"{run['_root']}: output {index} misses important species")
        if max(changes, default=math.inf) > policy["convergence_tolerance"]:
            inconclusive.append(
                f"output {index} maximum profile change {max(changes):.6g} exceeds tolerance"
            )
        for reaction in reactions:
            for product in reaction.get("products", []):
                if product.get("production_agreement") is not True:
                    inconclusive.append(
                        f"output {index} reaction {reaction.get('id')} product "
                        f"{product.get('species')} closure failed"
                    )
    return windows, inconclusive


def _load_cells(
    run: Dict[str, Any], species: str, output_index: int, policy: Dict[str, Any]
) -> List[Cell]:
    diagnostics = run["_root"] / "PT" / "diagnostics"
    table = diagnostics / f"stage3_cells.{species}.out={output_index}.csv"
    measure = diagnostics / "cell_measures.csv"
    try:
        cells, identity = load_cells_with_identity(table, measure)
    except (OSError, PlotError, ValueError) as exc:
        # Preserve the evaluator's four-outcome protocol for malformed CSV or
        # geometry.  Letting a plotting-layer exception escape would omit the
        # machine-readable INVALID_EVIDENCE summary.
        raise GateError(f"{table}: invalid center-cell evidence: {exc}") from exc
    if (
        identity.species != species
        or identity.output_index != output_index
        or identity.sampling_length != policy["sampling_length"]
        or identity.random_seed != run["seed"]
        or not identity.random_seed_configured
        or identity.residual_sampling != run["mode"]
        or identity.mpi_ranks != policy["mpi_ranks"]
    ):
        raise GateError(f"{table}: center-cell provenance mismatch")

    inner = policy["nucleus_radius_m"] + 1000.0 * policy["comparison_altitude_km"][0]
    outer = policy["nucleus_radius_m"] + 1000.0 * policy["comparison_altitude_km"][1]
    selected = [cell for cell in cells if cell.inner_m >= inner and cell.outer_m <= outer]
    if not selected:
        raise GateError(f"{table}: no exact cells are fully inside the approved interval")
    for cell in selected:
        if (
            cell.block_index < 0
            or cell.local_cell_index < 0
            or cell.block_cell_count <= 0
            or cell.population_target_min < 1
            or cell.population_target_max < 2 * cell.population_target_min
            or cell.population_class < 0
            or cell.population_policy_source < 0
        ):
            raise GateError(f"{table}: missing authoritative block-policy metadata")
        if (
            [cell.population_target_min, cell.population_target_max]
            != policy["chemical_only_population_target"]
            or cell.population_locked
            or cell.population_class != 2
            or cell.population_policy_source != 2
        ):
            raise GateError(
                f"{table}: important species does not use the approved unlocked "
                "CHEMICAL_ONLY population policy"
            )
    return selected


def _same_geometry(all_cells: Sequence[Sequence[Cell]]) -> None:
    reference = all_cells[0]
    for run_index, cells in enumerate(all_cells[1:], start=1):
        if len(cells) != len(reference):
            raise GateError(f"run {run_index} has a different approved cell count")
        for expected, actual in zip(reference, cells):
            if (
                expected.index != actual.index
                or expected.block_index != actual.block_index
                or expected.local_cell_index != actual.local_cell_index
                or expected.block_cell_count != actual.block_cell_count
                or expected.population_target_min != actual.population_target_min
                or expected.population_target_max != actual.population_target_max
                or expected.population_locked != actual.population_locked
                or expected.population_class != actual.population_class
                or expected.population_policy_source != actual.population_policy_source
                or not math.isclose(expected.inner_m, actual.inner_m, rel_tol=2e-13)
                or not math.isclose(expected.outer_m, actual.outer_m, rel_tol=2e-13)
                or not math.isclose(expected.volume_m3, actual.volume_m3, rel_tol=2e-13)
            ):
                raise GateError(f"run {run_index} has mismatched exact-cell geometry")


def derive_common_support(
    all_cells: Sequence[Sequence[Cell]], sampling_length: int, minimum_exposure: float
) -> List[Dict[str, Any]]:
    """Close a bin only when every one of the eight runs reaches exposure."""
    if not all_cells:
        raise GateError("common support requires run cells")
    _same_geometry(all_cells)
    support: List[Dict[str, Any]] = []
    start = 0
    exposure = [0.0] * len(all_cells)
    for offset in range(len(all_cells[0])):
        for run_index, cells in enumerate(all_cells):
            exposure[run_index] += cells[offset].occupancy * sampling_length
        if min(exposure) >= minimum_exposure or offset + 1 == len(all_cells[0]):
            first = all_cells[0][start]
            last = all_cells[0][offset]
            support.append(
                {
                    "bin_index": len(support),
                    "first_offset": start,
                    "last_offset": offset,
                    "first_cell": first.index,
                    "last_cell": last.index,
                    "r_inner_m": first.inner_m,
                    "r_outer_m": last.outer_m,
                    "minimum_run_exposure": min(exposure),
                    "qualified": min(exposure) >= minimum_exposure,
                }
            )
            start = offset + 1
            exposure = [0.0] * len(all_cells)
    return support


def _support_statistics(
    runs: Sequence[Sequence[Cell]], support: Sequence[Dict[str, Any]], sampling_length: int
) -> List[Dict[str, Any]]:
    result: List[Dict[str, Any]] = []
    for bounds in support:
        first = bounds["first_offset"]
        last = bounds["last_offset"] + 1
        density: List[float] = []
        exposure: List[float] = []
        for cells in runs:
            group = cells[first:last]
            volume = sum(cell.volume_m3 for cell in group)
            density.append(sum(cell.density_m3 * cell.volume_m3 for cell in group) / volume)
            exposure.append(sum(cell.occupancy for cell in group) * sampling_length)
        mean = statistics.fmean(density)
        standard_deviation = statistics.stdev(density)
        standard_error = standard_deviation / math.sqrt(len(density))
        result.append(
            {
                **bounds,
                "density_by_seed_m-3": density,
                "mean_density_m-3": mean,
                "sample_standard_deviation_m-3": standard_deviation,
                "standard_error_m-3": standard_error,
                "coefficient_of_variation": (
                    standard_deviation / abs(mean) if mean != 0.0 else None
                ),
                "minimum_mode_exposure": min(exposure),
                "zero_seed_fraction": sum(value <= 0.0 for value in density)
                / len(density),
            }
        )
    return result


def _population_coverage(cells: Sequence[Cell]) -> Dict[str, Any]:
    by_block: Dict[int, List[Cell]] = defaultdict(list)
    for cell in cells:
        by_block[cell.block_index].append(cell)
    populated: List[float] = []
    within = 0
    complete_blocks = 0
    for block_cells in by_block.values():
        expected = block_cells[0].block_cell_count
        if len(block_cells) != expected:
            continue
        complete_blocks += 1
        if sorted(cell.local_cell_index for cell in block_cells) != list(range(expected)):
            raise GateError("complete block has duplicate or missing local-cell indices")
        if any(
            cell.population_target_min != block_cells[0].population_target_min
            or cell.population_target_max != block_cells[0].population_target_max
            or cell.block_cell_count != expected
            or cell.population_locked != block_cells[0].population_locked
            or cell.population_class != block_cells[0].population_class
            or cell.population_policy_source != block_cells[0].population_policy_source
            for cell in block_cells
        ):
            raise GateError("population policy changes inside one block")
        mean = statistics.fmean(cell.occupancy for cell in block_cells)
        if mean > 0.0:
            populated.append(mean)
            if block_cells[0].population_target_min <= mean <= block_cells[0].population_target_max:
                within += 1
    return {
        "complete_blocks": complete_blocks,
        "populated_blocks": len(populated),
        "blocks_within_target": within,
        "within_target_fraction": within / len(populated) if populated else 0.0,
    }


def _rubin_metrics(
    systematic: Dict[str, Sequence[Sequence[Cell]]], policy: Dict[str, Any]
) -> Dict[str, float]:
    reference_path = (REPO_ROOT / policy["rubin_reference_csv"]).resolve()
    grouped: Dict[str, List[Tuple[float, float]]] = defaultdict(list)
    with reference_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if (
                row["species"] in policy["important_species"]
                and row["variant"] == policy["rubin_variant"]
                and row["quantity"] == policy["rubin_quantity"]
            ):
                grouped[row["species"]].append(
                    (float(row["radius_km"]), float(row["value"]) * 1.0e6)
                )
    result: Dict[str, float] = {}
    for species, runs in systematic.items():
        centers = [cell.radius_m / 1000.0 for cell in runs[0]]
        means = [statistics.fmean(run[i].density_m3 for run in runs) for i in range(len(centers))]
        errors: List[float] = []
        for radius, reference in grouped.get(species, []):
            if radius < centers[0] or radius > centers[-1] or reference <= 0.0:
                continue
            right = next((i for i, value in enumerate(centers) if value >= radius), None)
            if right is None or right == 0:
                continue
            left = right - 1
            if means[left] <= 0.0 or means[right] <= 0.0:
                continue
            fraction = (math.log(radius) - math.log(centers[left])) / (
                math.log(centers[right]) - math.log(centers[left])
            )
            model_log = math.log10(means[left]) + fraction * (
                math.log10(means[right]) - math.log10(means[left])
            )
            errors.append(model_log - math.log10(reference))
        if len(errors) < 10:
            raise GateError(f"insufficient positive Rubin overlap for {species}")
        result[species] = math.sqrt(statistics.fmean(value * value for value in errors))
    return result


def evaluate(archive: Path, policy_path: Path, baseline_path: Path | None) -> Dict[str, Any]:
    policy = load_policy(policy_path)
    runs = _manifest_runs(archive, policy)
    inconclusive: List[str] = []
    quality_failures: List[str] = []
    windows: Dict[Tuple[str, int], Dict[int, Dict[str, Any]]] = {}
    for run in runs:
        key = (run["mode"], run["seed"])
        windows[key], reasons = _load_complete_windows(run, policy)
        inconclusive.extend(f"{key}: {reason}" for reason in reasons)
        # Sampling-population pressure is evaluated immediately after a
        # report is written. At the 40k hard limit a request can therefore be
        # raised after the final JSON window with no subsequent file in which
        # to serialize it. Preserve the structured counters as authoritative
        # for all executed events and use this one log record only to catch
        # that final pending request.
        log_path = run["_root"] / "run.log"
        try:
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            raise GateError(f"cannot read run log {log_path}: {exc}") from exc
        for match in re.finditer(
            r"radial sampling-population emergency requested iteration=(\d+)",
            log_text,
        ):
            if int(match.group(1)) >= policy["burn_in_iterations"]:
                inconclusive.append(f"{key}: post-burn-in sampling emergency request")

    accepted_indices = list(range(1, policy["required_converged_windows"] + 1))
    final_index = accepted_indices[-1]
    species_results: Dict[str, Any] = {}
    final_systematic: Dict[str, List[List[Cell]]] = {}
    for species in policy["important_species"]:
        final_by_mode: Dict[str, List[List[Cell]]] = {mode: [] for mode in policy["modes"]}
        population_by_mode: Dict[str, List[Dict[str, Any]]] = {
            mode: [] for mode in policy["modes"]
        }
        zero_by_mode: Dict[str, List[float]] = {mode: [] for mode in policy["modes"]}
        for run in runs:
            for output_index in accepted_indices:
                cells = _load_cells(run, species, output_index, policy)
                coverage = _population_coverage(cells)
                coverage.update({"seed": run["seed"], "output_index": output_index})
                population_by_mode[run["mode"]].append(coverage)
                if run["mode"] == "SYSTEMATIC" and coverage["within_target_fraction"] < policy[
                    "systematic_population_required_block_fraction"
                ]:
                    quality_failures.append(
                        f"{species} seed {run['seed']} output {output_index} population coverage "
                        f"{coverage['within_target_fraction']:.6g}"
                    )
                if output_index == final_index:
                    final_by_mode[run["mode"]].append(cells)
                    zero_fraction = sum(cell.density_m3 <= 0.0 for cell in cells) / len(cells)
                    zero_by_mode[run["mode"]].append(zero_fraction)
                    if run["mode"] == "SYSTEMATIC" and zero_fraction > policy[
                        "systematic_raw_zero_cell_fraction_maximum"
                    ]:
                        quality_failures.append(
                            f"{species} seed {run['seed']} raw zero fraction {zero_fraction:.6g}"
                        )
        all_final = final_by_mode["BERNOULLI"] + final_by_mode["SYSTEMATIC"]
        support = derive_common_support(
            all_final,
            policy["sampling_length"],
            policy["minimum_particle_iteration_exposure"],
        )
        stats = {
            mode: _support_statistics(final_by_mode[mode], support, policy["sampling_length"])
            for mode in policy["modes"]
        }
        qualified = [row for row in support if row["qualified"]]
        if not qualified or any(not row["qualified"] for row in support):
            quality_failures.append(f"{species} has an underexposed common-support bin")

        systematic_qualified = [row for row in stats["SYSTEMATIC"] if row["qualified"]]
        passing_cv = 0
        for row in systematic_qualified:
            cv = row["coefficient_of_variation"]
            if row["mean_density_m-3"] <= 0.0 or cv is None:
                quality_failures.append(f"{species} has an empty systematic adaptive bin")
                continue
            if cv <= policy["systematic_cv_target"]:
                passing_cv += 1
            if cv > policy["systematic_cv_hard_maximum"]:
                quality_failures.append(f"{species} systematic CV {cv:.6g} exceeds hard maximum")
        fraction = passing_cv / len(systematic_qualified) if systematic_qualified else 0.0
        if fraction < policy["systematic_cv_required_bin_fraction"]:
            quality_failures.append(f"{species} systematic CV passing fraction {fraction:.6g}")

        agreement: List[Dict[str, Any]] = []
        for bernoulli, systematic in zip(stats["BERNOULLI"], stats["SYSTEMATIC"]):
            if bernoulli["qualified"] and bernoulli["mean_density_m-3"] <= 0.0:
                inconclusive.append(
                    f"{species} has an empty Bernoulli adaptive bin "
                    f"{bernoulli['bin_index']}"
                )
            scale = max(abs(bernoulli["mean_density_m-3"]), abs(systematic["mean_density_m-3"]))
            tolerance = max(
                policy["mode_mean_relative_tolerance"] * scale,
                policy["mode_mean_combined_standard_error_multiplier"]
                * math.hypot(
                    bernoulli["standard_error_m-3"], systematic["standard_error_m-3"]
                ),
            )
            difference = abs(bernoulli["mean_density_m-3"] - systematic["mean_density_m-3"])
            passed = difference <= tolerance
            agreement.append({"bin_index": bernoulli["bin_index"], "difference": difference, "tolerance": tolerance, "passed": passed})
            if bernoulli["qualified"] and not passed:
                inconclusive.append(f"{species} mode means disagree in bin {bernoulli['bin_index']}")
        species_results[species] = {
            "support": support,
            "mode_statistics": stats,
            "mode_agreement": agreement,
            "population": population_by_mode,
            "raw_zero_cell_fraction": zero_by_mode,
        }
        final_systematic[species] = final_by_mode["SYSTEMATIC"]

    rubin = _rubin_metrics(final_systematic, policy)
    candidate_baseline = {
        "schema": "amps.comet.stage3-rubin-baseline.v1",
        "accepted": False,
        "metric": "systematic-ensemble-log10-density-rmse-dex",
        "species_rmse_dex": rubin,
        "policy_sha256": sha256_json(policy),
    }
    evaluation_dir = archive / "evaluation"
    atomic_write_json(evaluation_dir / "rubin_candidate_baseline.json", candidate_baseline)
    rubin_status = "BASELINE_REQUIRED"
    if baseline_path is None and policy.get("accepted_rubin_baseline"):
        baseline_path = (REPO_ROOT / policy["accepted_rubin_baseline"]).resolve()
    if baseline_path is not None:
        baseline = _load_json(baseline_path)
        if baseline.get("schema") != candidate_baseline["schema"] or baseline.get("accepted") is not True:
            raise GateError("Rubin baseline is not explicitly accepted")
        if baseline.get("policy_sha256") != sha256_json(policy):
            raise GateError("Rubin baseline was accepted under a different science policy")
        baseline_values = baseline.get("species_rmse_dex")
        if not isinstance(baseline_values, dict) or set(baseline_values) != set(rubin):
            raise GateError("Rubin baseline species do not match the gate")
        if any(
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(float(value))
            or float(value) < 0.0
            for value in baseline_values.values()
        ):
            raise GateError("Rubin baseline contains an invalid regression metric")
        regressions = {
            species: rubin[species] / float(baseline_values[species]) - 1.0
            if float(baseline_values[species]) > 0.0
            else (0.0 if rubin[species] == 0.0 else math.inf)
            for species in rubin
        }
        failed = [
            species
            for species, value in regressions.items()
            if value > policy["rubin_maximum_regression_fraction"]
        ]
        rubin_status = "PASS" if not failed else "REGRESSION"
        if failed:
            inconclusive.append("Rubin regression exceeds 10% for " + ", ".join(failed))
    else:
        inconclusive.append("accepted Rubin regression baseline is not configured")

    if inconclusive:
        outcome = "INCONCLUSIVE_RERUN"
    elif quality_failures:
        outcome = "FAIL_CONSIDER_STAGE4"
    else:
        outcome = "PASS_SKIP_STAGE4"
    summary = {
        "schema": "amps.comet.stage3-science-gate.v1",
        "outcome": outcome,
        "exit_code": OUTCOME_EXIT[outcome],
        "policy": policy,
        "policy_sha256": sha256_json(policy),
        "archive": str(archive.resolve()),
        "accepted_output_index": final_index,
        "inconclusive_reasons": sorted(set(inconclusive)),
        "quality_failures": sorted(set(quality_failures)),
        "rubin_status": rubin_status,
        "rubin_species_rmse_dex": rubin,
        "species": species_results,
    }
    return summary


def _write_summaries(archive: Path, summary: Dict[str, Any]) -> None:
    output = archive / "evaluation"
    atomic_write_json(output / "stage3_science_gate.json", summary)
    with (output / "stage3_science_gate.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["record_type", "species", "mode", "metric", "value", "status"])
        writer.writerow(["gate", "", "", "outcome", summary["outcome"], summary["outcome"]])
        for species, result in summary.get("species", {}).items():
            for mode, rows in result["mode_statistics"].items():
                for row in rows:
                    writer.writerow(["bin", species, mode, f"cv[{row['bin_index']}]", row["coefficient_of_variation"], "qualified" if row["qualified"] else "underexposed"])
    markdown = [
        "# Stage 3 scientific gate",
        "",
        f"Outcome: `{summary['outcome']}`",
        "",
        f"Rubin regression: `{summary.get('rubin_status')}`",
        "",
    ]
    if summary.get("inconclusive_reasons"):
        markdown.extend(["## Inconclusive reasons", ""] + [f"- {item}" for item in summary["inconclusive_reasons"]] + [""])
    if summary.get("quality_failures"):
        markdown.extend(["## Statistical quality failures", ""] + [f"- {item}" for item in summary["quality_failures"]] + [""])
    (output / "stage3_science_gate.md").write_text("\n".join(markdown), encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument(
        "--policy",
        type=Path,
        default=THIS_DIR / "stage3_science_policy.json",
    )
    parser.add_argument("--rubin-baseline", type=Path)
    arguments = parser.parse_args(argv)
    try:
        summary = evaluate(
            arguments.archive.resolve(), arguments.policy.resolve(),
            arguments.rubin_baseline.resolve() if arguments.rubin_baseline else None,
        )
    except GateError as exc:
        summary = {
            "schema": "amps.comet.stage3-science-gate.v1",
            "outcome": "INVALID_EVIDENCE",
            "exit_code": OUTCOME_EXIT["INVALID_EVIDENCE"],
            "error": str(exc),
        }
        _write_summaries(arguments.archive.resolve(), summary)
        print(f"INVALID Stage-3 science gate: {exc}", file=sys.stderr)
        return OUTCOME_EXIT["INVALID_EVIDENCE"]
    _write_summaries(arguments.archive.resolve(), summary)
    print(
        f"{summary['outcome']} Stage-3 science gate "
        f"rubin={summary['rubin_status']} archive={arguments.archive.resolve()}"
    )
    return int(summary["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())

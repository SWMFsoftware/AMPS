#!/usr/bin/env python3
"""Build exposure-controlled radial diagnostics without changing AMPS output.

The original cell tables remain authoritative and are never edited.  This
tool joins the Stage-3 center-cell diagnostic to ``cell_measures.csv``, then merges
contiguous inner-to-outer cells for each species until the requested raw
model-particle-iteration exposure is reached.  The final bin is retained and
flagged when it remains below target; zero and non-positive estimates are
reported explicitly rather than removed for plotting convenience.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from plot_1d import PlotError


FILE_PATTERN = re.compile(
    r"^stage3_cells\.(?P<species>.+)\.out=(?P<output>\d+)\.csv$"
)


@dataclass(frozen=True)
class Cell:
    index: int
    radius_m: float
    inner_m: float
    outer_m: float
    volume_m3: float
    density_m3: float
    occupancy: float
    # Stage-3 science-gate metadata was added after the original adaptive-bin
    # format.  Defaults keep historical focused archives readable; the science
    # gate itself requires nonnegative values and never infers a block size.
    block_index: int = -1
    local_cell_index: int = -1
    block_cell_count: int = -1
    population_target_min: int = -1
    population_target_max: int = -1
    population_locked: bool = False
    population_class: int = -1
    population_policy_source: int = -1


@dataclass(frozen=True)
class TableIdentity:
    """Provenance that must agree across statistically combined runs."""

    species: str
    output_index: int
    sampling_length: int
    random_seed: int
    random_seed_configured: bool
    residual_sampling: str
    mpi_ranks: int


def read_measure_rows(path: Path) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"cell_index", "radius_m", "r_inner_m", "r_outer_m", "measure_m3"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise PlotError(
                f"{path}: Stage-3 bins require {', '.join(sorted(required))}"
            )
        for row in reader:
            rows.append({name: float(row[name]) for name in required})
    if not rows:
        raise PlotError(f"{path}: empty cell-measure table")
    for expected, row in enumerate(rows):
        if int(row["cell_index"]) != expected:
            raise PlotError(f"{path}: cell indices are not contiguous")
        if not (
            math.isfinite(row["radius_m"])
            and math.isfinite(row["r_inner_m"])
            and math.isfinite(row["r_outer_m"])
            and math.isfinite(row["measure_m3"])
            and row["r_inner_m"] < row["radius_m"] < row["r_outer_m"]
            and row["measure_m3"] > 0.0
        ):
            raise PlotError(f"{path}: invalid radial support at cell {expected}")
        if expected and not math.isclose(
            rows[expected - 1]["r_outer_m"], row["r_inner_m"], rel_tol=2e-13
        ):
            raise PlotError(f"{path}: radial supports have a gap or overlap")
    return rows


def load_cells_with_identity(
    table_path: Path, measure_path: Path
) -> Tuple[List[Cell], TableIdentity]:
    measures = read_measure_rows(measure_path)
    samples: List[Dict[str, float]] = []
    identity: TableIdentity | None = None
    name_match = FILE_PATTERN.match(table_path.name)
    if name_match is None:
        raise PlotError(
            f"{table_path}: expected stage3_cells.SPECIES.out=N.csv naming"
        )
    try:
        with table_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            required = {
                "cell_index",
                "radius_m",
                "r_inner_m",
                "r_outer_m",
                "measure_m3",
                "number_density_m-3",
                "mean_model_particles",
                "sampling_length",
                "output_index",
                "random_seed",
                "random_seed_configured",
                "residual_sampling",
                "mpi_ranks",
            }
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise PlotError(
                    f"{table_path}: expected a Stage-3 center-cell CSV; "
                    "Tecplot pic.*.dat coordinates are node-interpolated and "
                    "cannot be joined exactly to shell exposure"
                )
            policy_fields = {
                "block_index",
                "local_cell_index",
                "block_cell_count",
                "population_target_min_per_cell",
                "population_target_max_per_cell",
                "population_locked",
                "population_class",
                "population_policy_source",
            }
            has_policy_fields = policy_fields.issubset(reader.fieldnames)
            if policy_fields.intersection(reader.fieldnames) and not has_policy_fields:
                raise PlotError(
                    f"{table_path}: incomplete Stage-3 block-policy metadata"
                )
            for row_number, row in enumerate(reader, start=2):
                try:
                    row_identity = TableIdentity(
                        species=name_match.group("species"),
                        output_index=int(row["output_index"]),
                        sampling_length=int(row["sampling_length"]),
                        random_seed=int(row["random_seed"]),
                        random_seed_configured=(
                            int(row["random_seed_configured"]) == 1
                        ),
                        residual_sampling=row["residual_sampling"].upper(),
                        mpi_ranks=int(row["mpi_ranks"]),
                    )
                    sample = {
                        "cell_index": float(row["cell_index"]),
                        "radius_m": float(row["radius_m"]),
                        "r_inner_m": float(row["r_inner_m"]),
                        "r_outer_m": float(row["r_outer_m"]),
                        "measure_m3": float(row["measure_m3"]),
                        "number_density_m-3": float(row["number_density_m-3"]),
                        "mean_model_particles": float(row["mean_model_particles"]),
                    }
                    if has_policy_fields:
                        sample.update(
                            {
                                "block_index": int(row["block_index"]),
                                "local_cell_index": int(row["local_cell_index"]),
                                "block_cell_count": int(row["block_cell_count"]),
                                "population_target_min": int(
                                    row["population_target_min_per_cell"]
                                ),
                                "population_target_max": int(
                                    row["population_target_max_per_cell"]
                                ),
                                "population_locked": int(row["population_locked"]),
                                "population_class": int(row["population_class"]),
                                "population_policy_source": int(
                                    row["population_policy_source"]
                                ),
                            }
                        )
                except (TypeError, ValueError) as exc:
                    raise PlotError(
                        f"{table_path}:{row_number}: invalid Stage-3 sample metadata"
                    ) from exc
                if row_identity.output_index != int(name_match.group("output")):
                    raise PlotError(
                        f"{table_path}:{row_number}: filename/output_index mismatch"
                    )
                if (
                    row_identity.sampling_length <= 0
                    or row_identity.random_seed < 0
                    or int(row["random_seed_configured"]) not in (0, 1)
                    or (
                        row_identity.random_seed_configured
                        and row_identity.random_seed <= 0
                    )
                    or (
                        not row_identity.random_seed_configured
                        and row_identity.random_seed != 0
                    )
                    or row_identity.mpi_ranks <= 0
                    or row_identity.residual_sampling not in {"BERNOULLI", "SYSTEMATIC"}
                ):
                    raise PlotError(
                        f"{table_path}:{row_number}: invalid Stage-3 provenance"
                    )
                if identity is None:
                    identity = row_identity
                elif identity != row_identity:
                    raise PlotError(
                        f"{table_path}:{row_number}: metadata changes within one table"
                    )
                samples.append(sample)
    except OSError as exc:
        raise PlotError(f"cannot read Stage-3 cell table {table_path}: {exc}") from exc
    if len(samples) != len(measures):
        raise PlotError(
            f"{table_path}: {len(samples)} samples do not match "
            f"{len(measures)} exact cell measures"
        )
    for index, sample in enumerate(samples):
        if int(sample["cell_index"]) != index:
            raise PlotError(f"{table_path}: cell indices are not contiguous")
        if not math.isfinite(sample["number_density_m-3"]):
            raise PlotError(f"{table_path}: density contains a non-finite value")
        embedded_geometry = (
            sample["radius_m"],
            sample["r_inner_m"],
            sample["r_outer_m"],
            sample["measure_m3"],
        )
        reference_geometry = (
            measures[index]["radius_m"],
            measures[index]["r_inner_m"],
            measures[index]["r_outer_m"],
            measures[index]["measure_m3"],
        )
        if any(
            not math.isclose(actual, expected, rel_tol=2.0e-13, abs_tol=1.0e-10)
            for actual, expected in zip(embedded_geometry, reference_geometry)
        ):
            raise PlotError(
                f"{table_path}: embedded geometry disagrees with cell measures "
                f"at cell {index}"
            )
        if not (
            math.isfinite(sample["mean_model_particles"])
            and sample["mean_model_particles"] >= 0.0
        ):
            raise PlotError(
                f"{table_path}: mean model-particle occupancy must be nonnegative"
            )
        if "block_index" in sample:
            if not (
                sample["block_index"] >= 0
                and sample["local_cell_index"] >= 0
                and sample["block_cell_count"] > 0
                and sample["local_cell_index"] < sample["block_cell_count"]
                and sample["population_target_min"] > 0
                and sample["population_target_max"]
                >= 2 * sample["population_target_min"]
                and sample["population_locked"] in (0, 1)
                and sample["population_class"] >= 0
                and sample["population_policy_source"] >= 0
            ):
                raise PlotError(
                    f"{table_path}: invalid block-policy metadata at cell {index}"
                )
            # DIM1 AMPS leaves all contain the same compile-time cell count and
            # the C++ writer numbers blocks in radial traversal order.  Verify
            # that contract explicitly so a corrupted policy table cannot be
            # mistaken for authoritative block-mean population evidence.
            block_cell_count = int(sample["block_cell_count"])
            if (
                int(sample["block_index"]) != index // block_cell_count
                or int(sample["local_cell_index"]) != index % block_cell_count
            ):
                raise PlotError(
                    f"{table_path}: block/local-cell indices disagree at cell {index}"
                )
            if index:
                first = samples[0]
                invariant_policy = (
                    "block_cell_count",
                    "population_target_min",
                    "population_target_max",
                    "population_locked",
                    "population_class",
                    "population_policy_source",
                )
                if any(sample[name] != first[name] for name in invariant_policy):
                    raise PlotError(
                        f"{table_path}: species policy changes at cell {index}"
                    )
    if identity is None:
        raise PlotError(f"{table_path}: empty Stage-3 center-cell table")

    cells = [
        Cell(
            index=i,
            radius_m=samples[i]["radius_m"],
            inner_m=samples[i]["r_inner_m"],
            outer_m=samples[i]["r_outer_m"],
            volume_m3=samples[i]["measure_m3"],
            density_m3=samples[i]["number_density_m-3"],
            occupancy=samples[i]["mean_model_particles"],
            block_index=int(samples[i].get("block_index", -1)),
            local_cell_index=int(samples[i].get("local_cell_index", -1)),
            block_cell_count=int(samples[i].get("block_cell_count", -1)),
            population_target_min=int(
                samples[i].get("population_target_min", -1)
            ),
            population_target_max=int(
                samples[i].get("population_target_max", -1)
            ),
            population_locked=bool(samples[i].get("population_locked", 0)),
            population_class=int(samples[i].get("population_class", -1)),
            population_policy_source=int(
                samples[i].get("population_policy_source", -1)
            ),
        )
        for i, row in enumerate(measures)
    ]
    return cells, identity


def load_cells(table_path: Path, measure_path: Path) -> List[Cell]:
    """Compatibility wrapper used by unit tests and single-table clients."""
    return load_cells_with_identity(table_path, measure_path)[0]


def merge_cells(
    cells: Sequence[Cell], sampling_length: int, minimum_exposure: float
) -> List[Dict[str, object]]:
    if sampling_length <= 0 or not math.isfinite(minimum_exposure) or minimum_exposure <= 0:
        raise ValueError("sampling length and minimum exposure must be positive")
    if not cells:
        raise ValueError("at least one cell is required")

    result: List[Dict[str, object]] = []
    start = 0
    exposure = 0.0
    for index, cell in enumerate(cells):
        if cell.occupancy < 0.0 or not math.isfinite(cell.occupancy):
            raise ValueError("cell occupancy must be finite and nonnegative")
        exposure += cell.occupancy * sampling_length
        if exposure >= minimum_exposure or index + 1 == len(cells):
            group = cells[start : index + 1]
            volume = sum(entry.volume_m3 for entry in group)
            integrated = sum(entry.density_m3 * entry.volume_m3 for entry in group)
            nonpositive = sum(entry.density_m3 <= 0.0 for entry in group)
            status = "complete"
            if exposure == 0.0:
                status = "zero_exposure"
            elif exposure < minimum_exposure:
                status = "final_underexposed"
            result.append(
                {
                    "bin_index": len(result),
                    "first_cell": group[0].index,
                    "last_cell": group[-1].index,
                    "cell_count": len(group),
                    "r_inner_m": group[0].inner_m,
                    "r_outer_m": group[-1].outer_m,
                    "representative_radius_m": math.sqrt(
                        group[0].inner_m * group[-1].outer_m
                    ),
                    "volume_m3": volume,
                    "particle_time_exposure": exposure,
                    "mean_model_particles_per_cell": (
                        exposure / sampling_length / len(group)
                    ),
                    "minimum_requested_exposure": minimum_exposure,
                    "number_density_m-3": integrated / volume,
                    "integrated_number": integrated,
                    "nonpositive_cell_count": nonpositive,
                    "nonpositive_cell_fraction": nonpositive / len(group),
                    "status": status,
                    "uncertainty_method": "independent-seed-or-window-required",
                    "uncertainty_sample_count": 0,
                }
            )
            start = index + 1
            exposure = 0.0
    return result


def merge_replicate_cells(
    replicates: Sequence[Sequence[Cell]],
    sampling_length: int,
    minimum_exposure: float,
) -> List[Dict[str, object]]:
    """Merge independent runs on one common, exposure-controlled support.

    Bin boundaries are selected from the arithmetic mean model-particle
    occupancy across replicas.  This prevents a noisy minor-species run from
    choosing a different support for each seed and thereby making the
    between-seed variance meaningless.  Density remains volume weighted
    inside each run; uncertainty is then calculated across those independent
    run estimates (never across the already-correlated cells of one run).
    """
    if not replicates:
        raise ValueError("at least one replicate is required")
    reference = replicates[0]
    if not reference:
        raise ValueError("replicates must contain at least one cell")
    for run_index, cells in enumerate(replicates):
        if len(cells) != len(reference):
            raise ValueError(f"replicate {run_index} has a different cell count")
        for expected, actual in zip(reference, cells):
            geometry = (
                (expected.index, actual.index),
                (expected.radius_m, actual.radius_m),
                (expected.inner_m, actual.inner_m),
                (expected.outer_m, actual.outer_m),
                (expected.volume_m3, actual.volume_m3),
            )
            if geometry[0][0] != geometry[0][1] or any(
                not math.isclose(a, b, rel_tol=2.0e-13, abs_tol=1.0e-10)
                for a, b in geometry[1:]
            ):
                raise ValueError(
                    f"replicate {run_index} does not share the reference mesh"
                )
            policy = (
                "block_index",
                "local_cell_index",
                "block_cell_count",
                "population_target_min",
                "population_target_max",
                "population_locked",
                "population_class",
                "population_policy_source",
            )
            if any(getattr(expected, name) != getattr(actual, name) for name in policy):
                raise ValueError(
                    f"replicate {run_index} does not share block-policy metadata"
                )

    average_cells = [
        Cell(
            index=cell.index,
            radius_m=cell.radius_m,
            inner_m=cell.inner_m,
            outer_m=cell.outer_m,
            volume_m3=cell.volume_m3,
            density_m3=sum(run[i].density_m3 for run in replicates) / len(replicates),
            occupancy=sum(run[i].occupancy for run in replicates) / len(replicates),
            block_index=cell.block_index,
            local_cell_index=cell.local_cell_index,
            block_cell_count=cell.block_cell_count,
            population_target_min=cell.population_target_min,
            population_target_max=cell.population_target_max,
            population_locked=cell.population_locked,
            population_class=cell.population_class,
            population_policy_source=cell.population_policy_source,
        )
        for i, cell in enumerate(reference)
    ]
    rows = merge_cells(average_cells, sampling_length, minimum_exposure)
    for row in rows:
        first = int(row["first_cell"])
        last = int(row["last_cell"]) + 1
        run_density: List[float] = []
        run_exposure: List[float] = []
        for cells in replicates:
            group = cells[first:last]
            volume = sum(cell.volume_m3 for cell in group)
            run_density.append(
                sum(cell.density_m3 * cell.volume_m3 for cell in group) / volume
            )
            run_exposure.append(
                sum(cell.occupancy for cell in group) * sampling_length
            )

        count = len(run_density)
        mean = statistics.fmean(run_density)
        sample_std = statistics.stdev(run_density) if count >= 2 else None
        standard_error = (
            sample_std / math.sqrt(count) if sample_std is not None else None
        )
        cv = (
            sample_std / abs(mean)
            if sample_std is not None and mean != 0.0
            else None
        )
        zeros = sum(value <= 0.0 for value in run_density)
        row.update(
            {
                "number_density_m-3": mean,
                "integrated_number": mean * float(row["volume_m3"]),
                "density_sample_standard_deviation_m-3": sample_std,
                "density_standard_error_m-3": standard_error,
                "density_coefficient_of_variation": cv,
                "zero_or_negative_replicate_count": zeros,
                "zero_or_negative_replicate_fraction": zeros / count,
                "particle_time_exposure": statistics.fmean(run_exposure),
                "particle_time_exposure_min": min(run_exposure),
                "particle_time_exposure_max": max(run_exposure),
                "mean_model_particles_per_cell": (
                    statistics.fmean(run_exposure)
                    / sampling_length
                    / int(row["cell_count"])
                ),
                "uncertainty_method": (
                    "independent-seed-sample-standard-deviation"
                    if count >= 2
                    else "independent-seed-or-window-required"
                ),
                "uncertainty_sample_count": count,
                "decision_sample_status": (
                    "ready" if count >= 4 else "requires_at_least_four_independent_seeds"
                ),
            }
        )
    return rows


def write_outputs(prefix: Path, metadata: Dict[str, object], bins: Iterable[Dict[str, object]]) -> None:
    rows = list(bins)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")
    json_path.write_text(
        json.dumps({**metadata, "bins": rows}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "table", type=Path,
        help="one PT/diagnostics/stage3_cells.SPECIES.out=N.csv table",
    )
    parser.add_argument(
        "--replicate",
        action="append",
        default=[],
        type=Path,
        metavar="TABLE",
        help=(
            "an independent-seed table for the same species and mesh; repeat "
            "the option, with at least four total tables for decision-ready statistics"
        ),
    )
    parser.add_argument("--cell-measures", type=Path, required=True)
    parser.add_argument("--sampling-length", type=int, required=True)
    parser.add_argument("--minimum-exposure", type=float, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    arguments = parser.parse_args(argv)

    table_paths = [arguments.table, *arguments.replicate]
    resolved_paths = [path.resolve() for path in table_paths]
    if len(set(resolved_paths)) != len(resolved_paths):
        raise PlotError("replicate tables must be distinct files")
    loaded = [
        load_cells_with_identity(path, arguments.cell_measures)
        for path in table_paths
    ]
    replicas = [entry[0] for entry in loaded]
    identities = [entry[1] for entry in loaded]
    reference_identity = identities[0]
    for identity in identities[1:]:
        if (
            identity.species != reference_identity.species
            or identity.output_index != reference_identity.output_index
            or identity.sampling_length != reference_identity.sampling_length
            or identity.residual_sampling != reference_identity.residual_sampling
            or identity.mpi_ranks != reference_identity.mpi_ranks
            or identity.random_seed_configured
            != reference_identity.random_seed_configured
        ):
            raise PlotError(
                "replicates must share species, output index, sampling length, "
                "residual mode, and MPI rank count"
            )
    seeds = [identity.random_seed for identity in identities]
    if len(identities) > 1 and not reference_identity.random_seed_configured:
        raise PlotError(
            "independent-replicate statistics require an explicit random seed"
        )
    if len(set(seeds)) != len(seeds):
        raise PlotError("replicates must declare distinct random seeds")
    if arguments.sampling_length != reference_identity.sampling_length:
        raise PlotError(
            "--sampling-length does not match the value recorded in the table"
        )
    species = reference_identity.species
    bins = merge_replicate_cells(
        replicas, arguments.sampling_length, arguments.minimum_exposure
    )
    write_outputs(
        arguments.output,
        {
            "schema": "amps.comet.stage3-adaptive-radial.v1",
            "species": species,
            "source_table": str(arguments.table.resolve()),
            "source_tables": [str(path.resolve()) for path in table_paths],
            "source_cell_measures": str(arguments.cell_measures.resolve()),
            "sampling_length": arguments.sampling_length,
            "output_index": reference_identity.output_index,
            "random_seeds": seeds,
            "random_seed_configured": reference_identity.random_seed_configured,
            "residual_sampling": reference_identity.residual_sampling,
            "mpi_ranks": reference_identity.mpi_ranks,
            "replicate_provenance": "validated-distinct-seed-common-window",
            "exposure_definition": "mean model particle occupancy times completed sampling iterations",
            "merge_direction": "inner-to-outer contiguous cells",
            "bin_support_exposure_statistic": "arithmetic mean across independent seeds",
            "independent_sample_count": len(replicas),
            "decision_sample_status": (
                "ready"
                if len(replicas) >= 4
                else "requires_at_least_four_independent_seeds"
            ),
        },
        bins,
    )
    print(
        f"PASS Stage-3 adaptive bins species={species} cells={len(replicas[0])} "
        f"bins={len(bins)} replicas={len(replicas)} "
        f"minimum-exposure={arguments.minimum_exposure:g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

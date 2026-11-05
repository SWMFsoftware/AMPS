#!/usr/bin/env python3
"""Create machine-readable statistical diagnostics for DIM1 minor species.

The analyzer is read-only.  It never treats a missing legacy field as if it
were known: in particular, the convergence-weighted L1 change is reported as
unavailable unless exact AMPS cell measures are present in the table or are
provided explicitly with ``--cell-measures``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from plot_1d import PlotError, read_amps_table


FILE_PATTERN = re.compile(
    r"^pic\.(?P<species>.+)\.s=(?P<species_index>\d+)\.out=(?P<output_index>\d+)\.dat$"
)
DENSITY_COLUMN = "Number Density[1/m^3]"
OCCUPANCY_COLUMN = "Particle Number"
TIMESTEP_COLUMN = "Local Time Step"
WEIGHT_COLUMN = "log10(Local Particle Weight)"
MEASURE_COLUMN = "Cell Measure [m^3]"
DEFAULT_RADII_KM = (10.0, 100.0, 1_000.0, 10_000.0, 100_000.0)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_range(values: np.ndarray, positive: bool = False) -> Optional[List[float]]:
    valid = np.isfinite(values)
    if positive:
        valid &= values > 0.0
    if not np.any(valid):
        return None
    return [float(np.min(values[valid])), float(np.max(values[valid]))]


def column(variables: Sequence[str], table: np.ndarray, name: str, path: Path) -> np.ndarray:
    try:
        return table[:, variables.index(name)]
    except ValueError as exc:
        raise PlotError(f"{path}: required column {name!r} is absent") from exc


def discover(plots_dir: Path) -> Dict[int, Dict[str, Path]]:
    if not plots_dir.is_dir():
        raise PlotError(f"plot directory does not exist: {plots_dir}")
    outputs: Dict[int, Dict[str, Path]] = {}
    for path in sorted(plots_dir.iterdir()):
        match = FILE_PATTERN.match(path.name)
        if match is None:
            continue
        output = int(match.group("output_index"))
        species = match.group("species")
        if species in outputs.setdefault(output, {}):
            raise PlotError(f"duplicate {species} output at index {output}")
        outputs[output][species] = path
    if not outputs:
        raise PlotError(f"no AMPS species tables were found in {plots_dir}")
    return outputs


def select_output(
    outputs: Dict[int, Dict[str, Path]], species: Sequence[str], requested: Optional[int]
) -> int:
    if requested is not None:
        missing = [name for name in species if name not in outputs.get(requested, {})]
        if missing:
            raise PlotError(
                f"output {requested} is missing requested species: {', '.join(missing)}"
            )
        return requested
    complete = [
        index for index, files in outputs.items() if all(name in files for name in species)
    ]
    if not complete:
        raise PlotError("no output window contains every requested species")
    return max(complete)


def positive_log_interpolate(radius_m: np.ndarray, values: np.ndarray, target_m: float) -> Optional[float]:
    valid = np.isfinite(radius_m) & np.isfinite(values) & (radius_m > 0.0) & (values > 0.0)
    x = radius_m[valid]
    y = values[valid]
    if len(x) < 2 or target_m < x[0] or target_m > x[-1]:
        return None
    return float(np.exp(np.interp(np.log(target_m), np.log(x), np.log(y))))


def load_measure_table(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    radii: List[float] = []
    measures: List[float] = []
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not {"radius_m", "measure_m3"}.issubset(reader.fieldnames):
                raise PlotError(
                    f"{path}: measure CSV requires radius_m and measure_m3 columns"
                )
            for row in reader:
                radii.append(float(row["radius_m"]))
                measures.append(float(row["measure_m3"]))
    except OSError as exc:
        raise PlotError(f"cannot read cell measures {path}: {exc}") from exc
    radius = np.asarray(radii, dtype=float)
    measure = np.asarray(measures, dtype=float)
    if len(radius) == 0 or np.any(~np.isfinite(radius)) or np.any(~np.isfinite(measure)):
        raise PlotError(f"{path}: cell measures must be finite and nonempty")
    if np.any(measure <= 0.0) or np.any(np.diff(radius) <= 0.0):
        raise PlotError(f"{path}: radii must increase and cell measures must be positive")
    return radius, measure


def exact_measures(
    variables: Sequence[str], table: np.ndarray, radius: np.ndarray, path: Path,
    external: Optional[Tuple[np.ndarray, np.ndarray]],
) -> Tuple[Optional[np.ndarray], str]:
    if MEASURE_COLUMN in variables:
        measure = table[:, variables.index(MEASURE_COLUMN)]
        if np.any(~np.isfinite(measure)) or np.any(measure <= 0.0):
            raise PlotError(f"{path}: {MEASURE_COLUMN} must be finite and positive")
        return measure, "amps_output_column"
    if external is None:
        return None, "unavailable: exact AMPS cell measures were not supplied"
    measure_radius, measure = external
    if len(radius) != len(measure_radius) or not np.allclose(
        radius, measure_radius, rtol=0.0, atol=1.0e-8
    ):
        raise PlotError(f"{path}: external cell-measure radii do not match output radii")
    return measure, "external_exact_table"


def relative_l1(
    current: np.ndarray, previous: np.ndarray, measure: np.ndarray
) -> float:
    # This is algebraically identical to RadialSimulation.cpp:
    # sum |n_k V - n_(k-1) V| / max(sum |n_k V|, sum |n_(k-1) V|).
    current_integrated = current * measure
    previous_integrated = previous * measure
    difference = float(np.sum(np.abs(current_integrated - previous_integrated)))
    scale = max(
        float(np.sum(np.abs(current_integrated))),
        float(np.sum(np.abs(previous_integrated))),
    )
    return difference / scale if scale > 0.0 else 0.0


def convergence_metadata(run_log: Optional[Path], iteration: Optional[int]) -> Dict[str, object]:
    result: Dict[str, object] = {"record": None, "status": "not_provided"}
    if run_log is None:
        return result
    try:
        lines = run_log.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise PlotError(f"cannot read run log {run_log}: {exc}") from exc
    candidates: List[Tuple[int, str]] = []
    for line in lines:
        if "radial convergence" not in line:
            continue
        match = re.search(r"\biteration=(\d+)\b", line)
        if match is None:
            continue
        record_iteration = int(match.group(1))
        if iteration is None or record_iteration <= iteration:
            candidates.append((record_iteration, line.strip()))
    if not candidates:
        return result
    record_iteration, record = max(candidates)
    result.update({"record": record, "iteration": record_iteration})
    consecutive = re.search(r"\bconsecutive=(\d+)/(\d+)\b", record)
    if consecutive is None:
        result["status"] = "reported_but_unparsed"
    else:
        completed, required = (int(value) for value in consecutive.groups())
        result.update({"consecutive_completed": completed, "consecutive_required": required})
        result["status"] = "converged" if required > 0 and completed >= required else "nonconverged"
    return result


def load_comparison(path: Optional[Path]) -> Dict[float, Dict[str, str]]:
    if path is None:
        return {}
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            return {float(row["radius_km"]): row for row in csv.DictReader(stream)}
    except (OSError, KeyError, ValueError) as exc:
        raise PlotError(f"cannot read comparison table {path}: {exc}") from exc


def optional_float(value: Optional[str]) -> Optional[float]:
    if value is None or value.strip() == "":
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def summarize_species(
    species: str, path: Path, previous_path: Optional[Path], radii_km: Sequence[float],
    external_measures: Optional[Tuple[np.ndarray, np.ndarray]],
    comparison: Dict[float, Dict[str, str]], output_index: int,
) -> Dict[str, object]:
    variables, table = read_amps_table(path)
    radius = column(variables, table, "X", path)
    density = column(variables, table, DENSITY_COLUMN, path)
    occupancy = column(variables, table, OCCUPANCY_COLUMN, path)
    timestep = column(variables, table, TIMESTEP_COLUMN, path)
    weight = column(variables, table, WEIGHT_COLUMN, path)
    # AMPS repeats a block-face plotting node in both adjacent Tecplot zones.
    # Preserve those rows for occupancy statistics and accept nondecreasing
    # radius, while still rejecting a genuinely unsorted table.
    if np.any(~np.isfinite(radius)) or np.any(np.diff(radius) < 0.0):
        raise PlotError(f"{path}: radial coordinates must be finite and nondecreasing")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise PlotError(f"{path}: density must be finite and nonnegative")
    if np.any(~np.isfinite(occupancy)) or np.any(occupancy < 0.0):
        raise PlotError(f"{path}: particle occupancy must be finite and nonnegative")

    positive = occupancy > 0.0
    quantiles = np.quantile(occupancy, [0.25, 0.5, 0.75, 0.9, 0.99])
    positive_quantiles = (
        np.quantile(occupancy[positive], [0.25, 0.5, 0.75, 0.9, 0.99])
        if np.any(positive)
        else np.full(5, np.nan)
    )
    carrier = np.divide(weight, timestep, out=np.full_like(weight, np.nan), where=timestep > 0.0)
    result: Dict[str, object] = {
        "species": species,
        "output_index": output_index,
        "source_file": str(path.resolve()),
        "source_sha256": file_sha256(path),
        "row_count": int(len(table)),
        "occupancy": {
            "mean_all_cells": float(np.mean(occupancy)),
            "maximum": float(np.max(occupancy)),
            "zero_cell_count": int(np.count_nonzero(~positive)),
            "zero_cell_fraction": float(np.mean(~positive)),
            "quantiles_all_cells": dict(zip(("p25", "p50", "p75", "p90", "p99"), map(float, quantiles))),
            "quantiles_positive_cells": {
                name: None if math.isnan(value) else float(value)
                for name, value in zip(("p25", "p50", "p75", "p90", "p99"), positive_quantiles)
            },
            "first_nonzero_radius_m": float(radius[positive][0]) if np.any(positive) else None,
            "last_nonzero_radius_m": float(radius[positive][-1]) if np.any(positive) else None,
        },
        "carrier": {
            "base_weight_column": WEIGHT_COLUMN,
            "base_weight_range": finite_range(weight, positive=True),
            "local_timestep_s_range": finite_range(timestep, positive=True),
            "weight_over_timestep_s-1_range": finite_range(carrier, positive=True),
        },
        "requested_radii": [],
    }

    previous_density: Optional[np.ndarray] = None
    previous_radius: Optional[np.ndarray] = None
    previous_table: Optional[np.ndarray] = None
    previous_variables: Optional[Sequence[str]] = None
    if previous_path is not None:
        previous_variables, previous_table = read_amps_table(previous_path)
        previous_radius = column(previous_variables, previous_table, "X", previous_path)
        previous_density = column(previous_variables, previous_table, DENSITY_COLUMN, previous_path)

    for requested_radius_km in radii_km:
        current_value = positive_log_interpolate(radius, density, requested_radius_km * 1.0e3)
        previous_value = (
            positive_log_interpolate(previous_radius, previous_density, requested_radius_km * 1.0e3)
            if previous_radius is not None and previous_density is not None
            else None
        )
        comparison_row = comparison.get(float(requested_radius_km), {})
        amps_paper = optional_float(comparison_row.get(f"AMPS_{species}_density_cm-3"))
        reference = optional_float(comparison_row.get(f"Rubin2011_{species}_density_cm-3"))
        result["requested_radii"].append(
            {
                "radius_km": float(requested_radius_km),
                "density_m-3": current_value,
                "previous_density_m-3": previous_value,
                "relative_window_change": (
                    abs(current_value - previous_value) / max(abs(current_value), abs(previous_value))
                    if current_value is not None and previous_value is not None and max(abs(current_value), abs(previous_value)) > 0.0
                    else None
                ),
                "comparison_amps_density_cm-3": amps_paper,
                "comparison_reference_density_cm-3": reference,
                "amps_to_reference_ratio": (
                    amps_paper / reference
                    if amps_paper is not None and reference is not None and amps_paper > 0.0 and reference > 0.0
                    else None
                ),
            }
        )

    integrated: Dict[str, object] = {"value": None, "status": "previous window unavailable"}
    if previous_density is not None and previous_radius is not None and previous_table is not None and previous_variables is not None:
        if len(radius) != len(previous_radius) or not np.allclose(radius, previous_radius, rtol=0.0, atol=1.0e-8):
            integrated["status"] = "unavailable: window meshes differ"
        else:
            measures, measure_status = exact_measures(variables, table, radius, path, external_measures)
            if measures is None:
                integrated["status"] = measure_status
            else:
                integrated = {
                    "value": relative_l1(density, previous_density, measures),
                    "status": "available",
                    "measure_source": measure_status,
                    "definition": "sum(abs((n_current-n_previous)*V))/max(sum(abs(n_current*V)),sum(abs(n_previous*V)))",
                }
    result["integrated_relative_l1_change"] = integrated
    return result


def write_csv(path: Path, report: Dict[str, object]) -> None:
    fields = [
        "record_type", "species", "output_index", "radius_km", "metric", "value",
        "previous_value", "reference_value", "ratio", "status",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for species_report in report["species"]:
            occupancy = species_report["occupancy"]
            for metric in ("mean_all_cells", "maximum", "zero_cell_fraction"):
                writer.writerow({
                    "record_type": "summary", "species": species_report["species"],
                    "output_index": species_report["output_index"], "metric": metric,
                    "value": occupancy[metric],
                })
            integrated = species_report["integrated_relative_l1_change"]
            writer.writerow({
                "record_type": "integrated_change", "species": species_report["species"],
                "output_index": species_report["output_index"], "metric": "relative_l1",
                "value": integrated["value"], "status": integrated["status"],
            })
            for radius in species_report["requested_radii"]:
                writer.writerow({
                    "record_type": "radius", "species": species_report["species"],
                    "output_index": species_report["output_index"],
                    "radius_km": radius["radius_km"], "metric": "density_m-3",
                    "value": radius["density_m-3"],
                    "previous_value": radius["previous_density_m-3"],
                    "reference_value": radius["comparison_reference_density_cm-3"],
                    "ratio": radius["amps_to_reference_ratio"],
                })


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize DIM1 minor-species occupancy, carriers, convergence, and paper comparison.",
        epilog=(
            "example:\n  python3 srcHartley2/visualization/analyze_minor_species.py "
            "PT/plots --species C2 C2H --run-log run.log --comparison-csv "
            "PT/plots/rubin2011_comparison/comparison_profiles.csv"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("plots_dir", type=Path)
    parser.add_argument("--species", nargs="+", default=["C2"])
    parser.add_argument("--output-index", type=int)
    parser.add_argument("--sampling-length", type=int)
    parser.add_argument("--radii-km", type=float, nargs="+", default=DEFAULT_RADII_KM)
    parser.add_argument("--run-log", type=Path)
    parser.add_argument("--comparison-csv", type=Path)
    parser.add_argument("--cell-measures", type=Path, help="exact CSV with radius_m,measure_m3")
    parser.add_argument("--output-prefix", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.sampling_length is not None and args.sampling_length <= 0:
            raise PlotError("--sampling-length must be positive")
        if any(not math.isfinite(value) or value <= 0.0 for value in args.radii_km):
            raise PlotError("--radii-km values must be finite and positive")
        outputs = discover(args.plots_dir)
        output_index = select_output(outputs, args.species, args.output_index)
        previous_index = output_index - 1
        comparison = load_comparison(args.comparison_csv)
        measures = load_measure_table(args.cell_measures) if args.cell_measures else None
        iteration = (
            (output_index + 1) * args.sampling_length
            if args.sampling_length is not None else None
        )
        report: Dict[str, object] = {
            "schema": "amps.comet.minor-species.stage0.v1",
            "plots_directory": str(args.plots_dir.resolve()),
            "output_index": output_index,
            "sampling_length": args.sampling_length,
            "nominal_iteration": iteration,
            "convergence": convergence_metadata(args.run_log, iteration),
            "species": [],
        }
        for species in args.species:
            previous = outputs.get(previous_index, {}).get(species)
            report["species"].append(
                summarize_species(
                    species, outputs[output_index][species], previous,
                    args.radii_km, measures, comparison, output_index,
                )
            )
        prefix = args.output_prefix or args.plots_dir / "minor_species_stage0"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        json_path = prefix.with_suffix(".json")
        csv_path = prefix.with_suffix(".csv")
        json_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        write_csv(csv_path, report)
        print(f"wrote {json_path}")
        print(f"wrote {csv_path}")
        return 0
    except (PlotError, OSError, ValueError) as exc:
        print(f"analyze_minor_species.py: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

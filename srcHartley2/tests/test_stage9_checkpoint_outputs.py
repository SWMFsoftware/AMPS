#!/usr/bin/env python3
"""Validate paired one-window Stage-9 checkpoint smoke artifacts.

The production writer deliberately emits the same information at three
levels: a JSON checkpoint, a per-species row in the checkpoint CSV, and a
per-cell carrier CSV.  This checker does not trust one level to validate
another.  It recomputes all counts, sums, Kish effective sample sizes, occupied
cell counts, and correction extrema from the cell rows, then requires both
summaries to agree.  Missing, duplicate, non-finite, or ambiguous evidence is
an error.

This is a focused output-contract test.  It proves that one analog and one
enhanced smoke run used the same seed, decomposition, sampling window, and
species/cell layout.  It does not claim convergence or scientific acceptance.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


SCHEMA = "amps.comet.stage9-checkpoint.v1"
MODES = ("analog", "enhanced")
CHECKPOINT_KEYS = {
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
SPECIES_KEYS = {
    "symbol",
    "model_carriers",
    "physical_inventory",
    "sum_squared_weights",
    "carrier_ess",
    "occupied_cells",
    "minimum_positive_cell_ess",
    "correction_min",
    "correction_max",
}
SUMMARY_FIELDS = (
    "iteration",
    "output_index",
    "sampling_length",
    "estimator_mode",
    "random_seed",
    "mpi_ranks",
    "species",
    "model_carriers",
    "physical_inventory",
    "sum_squared_weights",
    "carrier_ess",
    "occupied_cells",
    "minimum_positive_cell_ess",
    "correction_min",
    "correction_max",
    "window_wall_seconds",
    "peak_rank_live_particles",
    "total_live_particles",
    "rank_capacity_min",
    "rank_capacity_max",
    "peak_buffer_utilization",
    "safe_boundary_audits",
    "pending_remainders",
    "pending_reaction_events",
)
CELL_FIELDS = (
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
)


class EvidenceError(ValueError):
    """Malformed, incomplete, or internally inconsistent checkpoint evidence."""


def exact_keys(value: Mapping[str, object], expected: Iterable[str], label: str) -> None:
    actual, wanted = set(value), set(expected)
    if actual != wanted:
        raise EvidenceError(
            f"{label} keys differ: missing={sorted(wanted-actual)} "
            f"extra={sorted(actual-wanted)}"
        )


def integer(value: object, label: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise EvidenceError(f"{label} is not an integer >= {minimum}")
    return value


def csv_integer(value: str, label: str, minimum: int = 0) -> int:
    try:
        result = int(value, 10)
    except (TypeError, ValueError) as exc:
        raise EvidenceError(f"{label} is not an integer") from exc
    if str(result) != value.strip() or result < minimum:
        raise EvidenceError(f"{label} is not a canonical integer >= {minimum}")
    return result


def number(value: object, label: str, minimum: float = 0.0) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EvidenceError(f"{label} is not numeric")
    result = float(value)
    if not math.isfinite(result) or result < minimum:
        raise EvidenceError(f"{label} is not finite and >= {minimum}")
    return result


def csv_number(value: str, label: str, minimum: float = 0.0) -> float:
    if value is None or value.strip() == "":
        raise EvidenceError(f"{label} is blank")
    try:
        result = float(value)
    except ValueError as exc:
        raise EvidenceError(f"{label} is not numeric") from exc
    if not math.isfinite(result) or result < minimum:
        raise EvidenceError(f"{label} is not finite and >= {minimum}")
    return result


def close(left: float, right: float, label: str) -> None:
    scale = max(abs(left), abs(right), 1.0)
    if abs(left - right) > 5.0e-13 * scale:
        raise EvidenceError(f"{label} differs: {left:.17e} != {right:.17e}")


def read_json(path: Path) -> Dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EvidenceError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise EvidenceError(f"JSON root is not an object: {path}")
    return value


def read_csv(path: Path, fields: Sequence[str]) -> List[Dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != list(fields):
                raise EvidenceError(
                    f"CSV header differs in {path}: {reader.fieldnames}"
                )
            rows = list(reader)
    except OSError as exc:
        raise EvidenceError(f"cannot read CSV {path}: {exc}") from exc
    if not rows:
        raise EvidenceError(f"CSV contains no data rows: {path}")
    return rows


def diagnostic_directory(case_root: Path) -> Path:
    nested = case_root / "PT" / "diagnostics"
    if nested.is_dir():
        return nested
    if case_root.name == "diagnostics" and case_root.is_dir():
        return case_root
    raise EvidenceError(f"case has no PT/diagnostics directory: {case_root}")


@dataclass(frozen=True)
class CaseEvidence:
    mode: str
    seed: int
    ranks: int
    sampling_length: int
    iteration: int
    output_index: int
    species: Tuple[str, ...]
    cell_count: int


def validate_case(mode: str, case_root: Path) -> CaseEvidence:
    diagnostics = diagnostic_directory(case_root)
    json_paths = sorted(diagnostics.glob("comet_stage9_checkpoint.out=*.json"))
    csv_paths = sorted(diagnostics.glob("comet_stage9_checkpoint.out=*.csv"))
    if len(json_paths) != 1 or len(csv_paths) != 1:
        raise EvidenceError(
            f"{mode}: expected exactly one checkpoint JSON and CSV, found "
            f"{len(json_paths)} and {len(csv_paths)}"
        )
    json_path, summary_path = json_paths[0], csv_paths[0]
    if json_path.stem != summary_path.stem:
        raise EvidenceError(f"{mode}: checkpoint JSON/CSV output indices differ")

    checkpoint = read_json(json_path)
    exact_keys(checkpoint, CHECKPOINT_KEYS, f"{mode} checkpoint")
    if checkpoint["schema"] != SCHEMA:
        raise EvidenceError(f"{mode}: unknown checkpoint schema")
    iteration = integer(checkpoint["iteration"], f"{mode} iteration", 1)
    output_index = integer(checkpoint["output_index"], f"{mode} output_index")
    sampling = integer(checkpoint["sampling_length"], f"{mode} sampling_length", 1)
    if checkpoint["estimator_mode"] != mode:
        raise EvidenceError(f"{mode}: estimator_mode provenance differs")
    seed = integer(checkpoint["random_seed"], f"{mode} random_seed")
    ranks = integer(checkpoint["mpi_ranks"], f"{mode} mpi_ranks", 1)
    wall = number(checkpoint["window_wall_seconds"], f"{mode} wall seconds")
    peak = integer(checkpoint["peak_rank_live_particles"], f"{mode} peak live")
    live = integer(checkpoint["total_live_particles"], f"{mode} total live")
    cap_min = integer(checkpoint["rank_capacity_min"], f"{mode} capacity min", 1)
    cap_max = integer(checkpoint["rank_capacity_max"], f"{mode} capacity max", 1)
    utilization = number(
        checkpoint["peak_buffer_utilization"], f"{mode} utilization"
    )
    audits_min = integer(
        checkpoint["safe_boundary_audits_min"], f"{mode} audits min", 1
    )
    audits_max = integer(
        checkpoint["safe_boundary_audits_max"], f"{mode} audits max", 1
    )
    if cap_min > cap_max or peak > cap_max or live > ranks * cap_max:
        raise EvidenceError(f"{mode}: impossible particle-buffer totals")
    if utilization > 1.0 or audits_min > audits_max:
        raise EvidenceError(f"{mode}: invalid utilization or audit range")
    # The rank producing the maximum occupancy can have any capacity between
    # the reported extrema, so only these enclosing bounds are justified.
    if utilization + 5.0e-13 < peak / cap_max or (
        peak > 0 and utilization - 5.0e-13 > peak / cap_min
    ):
        raise EvidenceError(f"{mode}: peak utilization is inconsistent with capacity")
    if integer(checkpoint["pending_remainders"], f"{mode} pending remainders") != 0:
        raise EvidenceError(f"{mode}: checkpoint has pending trajectory remainders")
    if integer(checkpoint["pending_reaction_events"], f"{mode} pending events") != 0:
        raise EvidenceError(f"{mode}: checkpoint has pending reaction events")

    raw_species = checkpoint["species"]
    if not isinstance(raw_species, list) or not raw_species:
        raise EvidenceError(f"{mode}: species list is empty or malformed")
    json_species: Dict[str, Dict[str, object]] = {}
    species_order: List[str] = []
    for index, raw in enumerate(raw_species):
        if not isinstance(raw, dict):
            raise EvidenceError(f"{mode}: species row {index} is not an object")
        exact_keys(raw, SPECIES_KEYS, f"{mode} species row {index}")
        symbol = raw["symbol"]
        if not isinstance(symbol, str) or not symbol or symbol in json_species:
            raise EvidenceError(f"{mode}: duplicate or invalid species symbol")
        json_species[symbol] = raw
        species_order.append(symbol)

    summary_rows = read_csv(summary_path, SUMMARY_FIELDS)
    if len(summary_rows) != len(species_order):
        raise EvidenceError(f"{mode}: checkpoint CSV species count differs from JSON")
    summary_by_species: Dict[str, Dict[str, str]] = {}
    for row in summary_rows:
        symbol = row["species"]
        if symbol not in json_species or symbol in summary_by_species:
            raise EvidenceError(f"{mode}: duplicate or unknown summary species {symbol}")
        summary_by_species[symbol] = row
        integer_fields = {
            "iteration": iteration,
            "output_index": output_index,
            "sampling_length": sampling,
            "random_seed": seed,
            "mpi_ranks": ranks,
            "peak_rank_live_particles": peak,
            "total_live_particles": live,
            "rank_capacity_min": cap_min,
            "rank_capacity_max": cap_max,
            "safe_boundary_audits": audits_min,
            "pending_remainders": 0,
            "pending_reaction_events": 0,
        }
        for field, expected in integer_fields.items():
            if csv_integer(row[field], f"{mode}/{symbol} {field}") != expected:
                raise EvidenceError(f"{mode}/{symbol}: {field} provenance differs")
        if row["estimator_mode"] != mode:
            raise EvidenceError(f"{mode}/{symbol}: summary mode differs")
        close(csv_number(row["window_wall_seconds"], f"{mode}/{symbol} wall"), wall,
              f"{mode}/{symbol} wall")
        close(csv_number(row["peak_buffer_utilization"], f"{mode}/{symbol} util"),
              utilization, f"{mode}/{symbol} utilization")

    carrier_paths = sorted(diagnostics.glob("stage9_carriers.*.out=*.csv"))
    if len(carrier_paths) != len(species_order):
        raise EvidenceError(f"{mode}: carrier-file count differs from species count")
    carrier_by_species: Dict[str, Path] = {}
    suffix = f".out={output_index}.csv"
    for path in carrier_paths:
        prefix = "stage9_carriers."
        if not path.name.startswith(prefix) or not path.name.endswith(suffix):
            raise EvidenceError(f"{mode}: malformed carrier filename {path.name}")
        symbol = path.name[len(prefix) : -len(suffix)]
        if symbol not in json_species or symbol in carrier_by_species:
            raise EvidenceError(f"{mode}: duplicate or unknown carrier species {symbol}")
        carrier_by_species[symbol] = path

    common_cell_count: int | None = None
    case_carriers = 0
    for symbol in species_order:
        rows = read_csv(carrier_by_species[symbol], CELL_FIELDS)
        if common_cell_count is None:
            common_cell_count = len(rows)
        elif len(rows) != common_cell_count:
            raise EvidenceError(f"{mode}: per-species cell counts differ")

        count = 0
        occupied = 0
        weight = 0.0
        weight2 = 0.0
        min_cell_ess = math.inf
        correction_min = math.inf
        correction_max = 0.0
        for expected_cell, row in enumerate(rows):
            if csv_integer(row["cell_index"], f"{mode}/{symbol} cell index") != expected_cell:
                raise EvidenceError(f"{mode}/{symbol}: cell indices are not contiguous")
            metadata = {
                "output_index": output_index,
                "sampling_length": sampling,
                "random_seed": seed,
                "mpi_ranks": ranks,
            }
            for field, expected in metadata.items():
                if csv_integer(row[field], f"{mode}/{symbol}/{expected_cell} {field}") != expected:
                    raise EvidenceError(f"{mode}/{symbol}: cell provenance differs")
            if row["estimator_mode"] != mode:
                raise EvidenceError(f"{mode}/{symbol}: cell mode differs")

            cell_count = csv_integer(
                row["model_carriers"], f"{mode}/{symbol}/{expected_cell} carriers"
            )
            cell_weight = csv_number(
                row["sum_weights"], f"{mode}/{symbol}/{expected_cell} sum_weights"
            )
            cell_weight2 = csv_number(
                row["sum_squared_weights"],
                f"{mode}/{symbol}/{expected_cell} sum_squared_weights",
            )
            cell_ess = csv_number(
                row["carrier_ess"], f"{mode}/{symbol}/{expected_cell} ESS"
            )
            if cell_count == 0:
                if cell_weight != 0.0 or cell_weight2 != 0.0 or cell_ess != 0.0:
                    raise EvidenceError(f"{mode}/{symbol}: empty cell has carrier mass")
                if row["correction_min"] != "" or row["correction_max"] != "":
                    raise EvidenceError(f"{mode}/{symbol}: empty cell has correction extrema")
            else:
                if not (cell_weight > 0.0 and cell_weight2 > 0.0):
                    raise EvidenceError(f"{mode}/{symbol}: occupied cell has nonpositive mass")
                expected_ess = cell_weight * cell_weight / cell_weight2
                close(cell_ess, expected_ess, f"{mode}/{symbol} cell ESS")
                if cell_ess > cell_count * (1.0 + 5.0e-13):
                    raise EvidenceError(f"{mode}/{symbol}: cell ESS exceeds carrier count")
                cmin = csv_number(
                    row["correction_min"], f"{mode}/{symbol}/{expected_cell} correction min",
                    minimum=0.0,
                )
                cmax = csv_number(
                    row["correction_max"], f"{mode}/{symbol}/{expected_cell} correction max",
                    minimum=0.0,
                )
                if not (cmin > 0.0 and cmax > 0.0) or cmin > cmax:
                    raise EvidenceError(f"{mode}/{symbol}: reversed correction extrema")
                occupied += 1
                min_cell_ess = min(min_cell_ess, cell_ess)
                correction_min = min(correction_min, cmin)
                correction_max = max(correction_max, cmax)
            count += cell_count
            weight += cell_weight
            weight2 += cell_weight2

        case_carriers += count
        summary = summary_by_species[symbol]
        js = json_species[symbol]
        json_count = integer(js["model_carriers"], f"{mode}/{symbol} JSON carriers")
        json_occupied = integer(js["occupied_cells"], f"{mode}/{symbol} JSON occupied")
        summary_count = csv_integer(summary["model_carriers"], f"{mode}/{symbol} summary carriers")
        summary_occupied = csv_integer(summary["occupied_cells"], f"{mode}/{symbol} summary occupied")
        if count != json_count or count != summary_count:
            raise EvidenceError(f"{mode}/{symbol}: aggregate carrier count differs")
        if occupied != json_occupied or occupied != summary_occupied:
            raise EvidenceError(f"{mode}/{symbol}: occupied-cell count differs")
        for value, json_key, csv_key in (
            (weight, "physical_inventory", "physical_inventory"),
            (weight2, "sum_squared_weights", "sum_squared_weights"),
        ):
            close(number(js[json_key], f"{mode}/{symbol} JSON {json_key}"), value,
                  f"{mode}/{symbol} JSON {json_key}")
            close(csv_number(summary[csv_key], f"{mode}/{symbol} CSV {csv_key}"), value,
                  f"{mode}/{symbol} CSV {csv_key}")
        aggregate_ess = weight * weight / weight2 if weight2 > 0.0 else 0.0
        close(number(js["carrier_ess"], f"{mode}/{symbol} JSON ESS"), aggregate_ess,
              f"{mode}/{symbol} JSON ESS")
        close(csv_number(summary["carrier_ess"], f"{mode}/{symbol} CSV ESS"), aggregate_ess,
              f"{mode}/{symbol} CSV ESS")

        if occupied == 0:
            if any(js[key] is not None for key in (
                "minimum_positive_cell_ess", "correction_min", "correction_max"
            )):
                raise EvidenceError(f"{mode}/{symbol}: empty JSON species has extrema")
            if any(summary[key] != "" for key in (
                "minimum_positive_cell_ess", "correction_min", "correction_max"
            )):
                raise EvidenceError(f"{mode}/{symbol}: empty CSV species has extrema")
        else:
            for key, expected in (
                ("minimum_positive_cell_ess", min_cell_ess),
                ("correction_min", correction_min),
                ("correction_max", correction_max),
            ):
                close(number(js[key], f"{mode}/{symbol} JSON {key}",
                             minimum=0.0), expected,
                      f"{mode}/{symbol} JSON {key}")
                close(csv_number(summary[key], f"{mode}/{symbol} CSV {key}",
                                 minimum=0.0), expected,
                      f"{mode}/{symbol} CSV {key}")
                if expected <= 0.0:
                    raise EvidenceError(f"{mode}/{symbol}: nonpositive {key}")

    if case_carriers <= 0:
        raise EvidenceError(f"{mode}: smoke checkpoint contains no carriers")
    assert common_cell_count is not None
    return CaseEvidence(mode, seed, ranks, sampling, iteration, output_index,
                        tuple(species_order), common_cell_count)


def validate_cases(cases: Mapping[str, Path]) -> Dict[str, CaseEvidence]:
    if set(cases) != set(MODES):
        raise EvidenceError(f"cases must be exactly {list(MODES)}")
    evidence = {mode: validate_case(mode, cases[mode]) for mode in MODES}
    analog, enhanced = evidence["analog"], evidence["enhanced"]
    for field in ("seed", "ranks", "sampling_length", "iteration", "output_index",
                  "species", "cell_count"):
        if getattr(analog, field) != getattr(enhanced, field):
            raise EvidenceError(f"analog/enhanced {field} provenance differs")
    return evidence


def parse_case(text: str) -> Tuple[str, Path]:
    mode, separator, raw_path = text.partition("=")
    if not separator or mode not in MODES or not raw_path:
        raise argparse.ArgumentTypeError(
            "case must be analog=/path or enhanced=/path"
        )
    return mode, Path(raw_path).resolve()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", action="append", type=parse_case, required=True,
        metavar="MODE=PATH",
        help="one-window case root; supply analog=PATH and enhanced=PATH",
    )
    arguments = parser.parse_args(argv)
    cases: Dict[str, Path] = {}
    for mode, path in arguments.case:
        if mode in cases:
            parser.error(f"duplicate --case mode: {mode}")
        cases[mode] = path
    try:
        evidence = validate_cases(cases)
    except EvidenceError as exc:
        print(f"FAIL Stage-9 checkpoint smoke evidence: {exc}", file=sys.stderr)
        return 1
    reference = evidence["analog"]
    print(
        "PASS Stage-9 checkpoint smoke "
        f"modes=analog,enhanced seed={reference.seed} ranks={reference.ranks} "
        f"sampling-length={reference.sampling_length} "
        f"species={len(reference.species)} cells={reference.cell_count}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

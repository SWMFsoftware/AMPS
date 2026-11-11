#!/usr/bin/env python3
"""Synthetic mutation tests for the strict Stage-9 checkpoint validator."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List


MODULE_PATH = Path(__file__).with_name("test_stage9_checkpoint_outputs.py")
SPEC = importlib.util.spec_from_file_location("stage9_checkpoint", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
CHECKER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CHECKER
SPEC.loader.exec_module(CHECKER)


def write_csv(path: Path, fields: tuple[str, ...], rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_case(root: Path, mode: str) -> Path:
    diagnostics = root / mode / "PT" / "diagnostics"
    diagnostics.mkdir(parents=True)
    output = 0
    h2o_cells: List[Dict[str, object]] = [
        {
            "cell_index": 0, "model_carriers": 2, "sum_weights": "5.00000000000000000e+00",
            "sum_squared_weights": "1.30000000000000000e+01",
            "carrier_ess": f"{25.0/13.0:.17e}", "correction_min": "1.00000000000000000e+00",
            "correction_max": "1.50000000000000000e+00", "output_index": output,
            "sampling_length": 1, "estimator_mode": mode, "random_seed": 104729,
            "mpi_ranks": 1,
        },
        {
            "cell_index": 1, "model_carriers": 1, "sum_weights": "4.00000000000000000e+00",
            "sum_squared_weights": "1.60000000000000000e+01",
            "carrier_ess": "1.00000000000000000e+00", "correction_min": "2.00000000000000000e+00",
            "correction_max": "2.00000000000000000e+00", "output_index": output,
            "sampling_length": 1, "estimator_mode": mode, "random_seed": 104729,
            "mpi_ranks": 1,
        },
    ]
    empty_cells: List[Dict[str, object]] = []
    for index in range(2):
        empty_cells.append({
            "cell_index": index, "model_carriers": 0,
            "sum_weights": "0.00000000000000000e+00",
            "sum_squared_weights": "0.00000000000000000e+00",
            "carrier_ess": "0.00000000000000000e+00", "correction_min": "",
            "correction_max": "", "output_index": output, "sampling_length": 1,
            "estimator_mode": mode, "random_seed": 104729, "mpi_ranks": 1,
        })
    write_csv(diagnostics / "stage9_carriers.H2O.out=0.csv", CHECKER.CELL_FIELDS, h2o_cells)
    write_csv(diagnostics / "stage9_carriers.C.out=0.csv", CHECKER.CELL_FIELDS, empty_cells)

    ess = 81.0 / 29.0
    minimum_cell_ess = min(25.0 / 13.0, 1.0)
    common = {
        "iteration": 1, "output_index": 0, "sampling_length": 1,
        "estimator_mode": mode, "random_seed": 104729, "mpi_ranks": 1,
        "window_wall_seconds": "5.00000000000000000e-01",
        "peak_rank_live_particles": 3, "total_live_particles": 3,
        "rank_capacity_min": 100, "rank_capacity_max": 100,
        "peak_buffer_utilization": "3.00000000000000000e-02",
        "safe_boundary_audits": 2, "pending_remainders": 0,
        "pending_reaction_events": 0,
    }
    summary_rows = [
        {
            **common, "species": "H2O", "model_carriers": 3,
            "physical_inventory": "9.00000000000000000e+00",
            "sum_squared_weights": "2.90000000000000000e+01",
            "carrier_ess": f"{ess:.17e}", "occupied_cells": 2,
            "minimum_positive_cell_ess": f"{minimum_cell_ess:.17e}",
            "correction_min": "1.00000000000000000e+00",
            "correction_max": "2.00000000000000000e+00",
        },
        {
            **common, "species": "C", "model_carriers": 0,
            "physical_inventory": "0.00000000000000000e+00",
            "sum_squared_weights": "0.00000000000000000e+00",
            "carrier_ess": "0.00000000000000000e+00", "occupied_cells": 0,
            "minimum_positive_cell_ess": "", "correction_min": "",
            "correction_max": "",
        },
    ]
    write_csv(diagnostics / "comet_stage9_checkpoint.out=0.csv", CHECKER.SUMMARY_FIELDS, summary_rows)
    checkpoint = {
        "schema": CHECKER.SCHEMA, "iteration": 1, "output_index": 0,
        "sampling_length": 1, "estimator_mode": mode, "random_seed": 104729,
        "mpi_ranks": 1, "window_wall_seconds": 0.5,
        "peak_rank_live_particles": 3, "total_live_particles": 3,
        "rank_capacity_min": 100, "rank_capacity_max": 100,
        "peak_buffer_utilization": 0.03, "safe_boundary_audits_min": 2,
        "safe_boundary_audits_max": 2, "pending_remainders": 0,
        "pending_reaction_events": 0,
        "species": [
            {
                "symbol": "H2O", "model_carriers": 3, "physical_inventory": 9.0,
                "sum_squared_weights": 29.0, "carrier_ess": ess,
                "occupied_cells": 2, "minimum_positive_cell_ess": minimum_cell_ess,
                "correction_min": 1.0, "correction_max": 2.0,
            },
            {
                "symbol": "C", "model_carriers": 0, "physical_inventory": 0.0,
                "sum_squared_weights": 0.0, "carrier_ess": 0.0,
                "occupied_cells": 0, "minimum_positive_cell_ess": None,
                "correction_min": None, "correction_max": None,
            },
        ],
    }
    (diagnostics / "comet_stage9_checkpoint.out=0.json").write_text(
        json.dumps(checkpoint), encoding="utf-8"
    )
    return root / mode


def mutate_json(case: Path, callback) -> None:
    path = case / "PT" / "diagnostics" / "comet_stage9_checkpoint.out=0.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    callback(value)
    path.write_text(json.dumps(value), encoding="utf-8")


def mutate_csv(path: Path, callback) -> None:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fields = tuple(reader.fieldnames or ())
        rows = list(reader)
    callback(rows)
    write_csv(path, fields, rows)


class CheckpointValidatorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.analog = make_case(root, "analog")
        self.enhanced = make_case(root, "enhanced")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def validate(self) -> None:
        CHECKER.validate_cases({"analog": self.analog, "enhanced": self.enhanced})

    def assert_rejected(self) -> None:
        with self.assertRaises(CHECKER.EvidenceError):
            self.validate()

    def test_valid_pair(self) -> None:
        self.validate()

    def test_missing_or_duplicate_checkpoint_is_rejected(self) -> None:
        path = self.analog / "PT" / "diagnostics" / "comet_stage9_checkpoint.out=0.json"
        path.unlink()
        self.assert_rejected()

    def test_mode_and_cross_arm_provenance_are_strict(self) -> None:
        mutate_json(self.enhanced, lambda value: value.__setitem__("random_seed", 7))
        self.assert_rejected()

    def test_pending_state_and_safe_boundary_are_strict(self) -> None:
        mutate_json(self.analog, lambda value: value.__setitem__("pending_remainders", 1))
        self.assert_rejected()

    def test_buffer_bounds_are_strict(self) -> None:
        mutate_json(self.analog, lambda value: value.__setitem__("peak_buffer_utilization", 1.1))
        self.assert_rejected()

    def test_summary_aggregate_corruption_is_rejected(self) -> None:
        path = self.enhanced / "PT" / "diagnostics" / "comet_stage9_checkpoint.out=0.csv"
        mutate_csv(path, lambda rows: rows[0].__setitem__("physical_inventory", "8.0"))
        self.assert_rejected()

    def test_cell_ess_and_correction_corruption_is_rejected(self) -> None:
        path = self.analog / "PT" / "diagnostics" / "stage9_carriers.H2O.out=0.csv"
        mutate_csv(path, lambda rows: rows[0].__setitem__("carrier_ess", "2.0"))
        self.assert_rejected()

    def test_empty_cell_semantics_are_strict(self) -> None:
        path = self.analog / "PT" / "diagnostics" / "stage9_carriers.C.out=0.csv"
        mutate_csv(path, lambda rows: rows[0].__setitem__("correction_min", "1.0"))
        self.assert_rejected()

    def test_species_and_cell_layout_must_match(self) -> None:
        path = self.enhanced / "PT" / "diagnostics" / "stage9_carriers.C.out=0.csv"
        mutate_csv(path, lambda rows: rows.pop())
        self.assert_rejected()

    def test_nonfinite_and_noncontiguous_cells_are_rejected(self) -> None:
        path = self.enhanced / "PT" / "diagnostics" / "stage9_carriers.H2O.out=0.csv"
        mutate_csv(path, lambda rows: rows[1].__setitem__("sum_weights", "nan"))
        self.assert_rejected()


if __name__ == "__main__":
    unittest.main(verbosity=2)

#!/usr/bin/env python3
"""Mutation tests for the Stage-9 transported-collision log parser.

These records are synthetic parser fixtures, not AMPS science evidence.  Each
mutation proves that a producer PASS token cannot conceal one broken physical
or ownership invariant.
"""

from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path

import test_stage9_transported_collision_logs as stage9


VALID = {
    "status": "PASS",
    "ranks": 1,
    "reaction": "C2_2C",
    "product": 2,
    "partner": 13,
    "owner": 0,
    "daughters": 16,
    "partners": 64,
    "daughter_changes": 16,
    "partner_changes": 48,
    "unchanged_records": 80,
    "unfinished": 0,
    "wrong_owner": 0,
    "cycles": 1,
    "trials": 128,
    "accepted": 64,
    "updates": 128,
    "probability_checks": 128,
    "invalid_probabilities": 0,
    "roundoff_clamps": 0,
    "probability_minimum": 1.0,
    "probability_maximum": 1.0,
    "correction_error": 0.0,
    "represented_rate_error": 0.0,
    "momentum_error": 1.0e-15,
}


def record(row: dict[str, object]) -> str:
    return (
        "{status} transported-daughter-collision ranks={ranks} "
        "reaction={reaction} product={product} partner={partner} owner={owner} "
        "daughters={daughters} partners={partners} "
        "daughter-velocity-changes={daughter_changes} "
        "partner-velocity-changes={partner_changes} "
        "unchanged-records={unchanged_records} "
        "unfinished-trajectories={unfinished} wrong-owner={wrong_owner} "
        "cycles={cycles} trials={trials} accepted={accepted} updates={updates} "
        "probability-checks={probability_checks} "
        "invalid-probabilities={invalid_probabilities} "
        "roundoff-clamps={roundoff_clamps} "
        "probability-min={probability_minimum:.17e} "
        "probability-max={probability_maximum:.17e} "
        "correction-error={correction_error:.17e} "
        "represented-rate-error={represented_rate_error:.17e} "
        "momentum-error={momentum_error:.17e}\n"
    ).format(**row)


class ParserTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.path = Path(self.temporary.name) / "run.log"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def parse(self, changes: dict[str, object] | None = None, ranks: int = 1):
        row = copy.deepcopy(VALID)
        row["ranks"] = ranks
        if ranks > 1:
            row["owner"] = ranks - 1
        if changes:
            row.update(changes)
        self.path.write_text(record(row), encoding="utf-8")
        return stage9.parse_log(self.path, ranks)

    def assert_rejected(self, changes: dict[str, object], text: str) -> None:
        with self.assertRaisesRegex(ValueError, text):
            self.parse(changes)

    def test_valid_record(self) -> None:
        self.assertEqual(self.parse()["daughters"], 16)

    def test_status_and_unique_record(self) -> None:
        self.assert_rejected({"status": "FAIL"}, "not PASS")
        self.path.write_text("no record\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            stage9.parse_log(self.path, 1)
        self.path.write_text(record(VALID) + record(VALID), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            stage9.parse_log(self.path, 1)

    def test_identity_and_owner(self) -> None:
        self.assert_rejected({"reaction": "H2O_OH_H"}, "unexpected reaction")
        self.assert_rejected({"partner": 2}, "product/partner")
        self.assert_rejected({"owner": 1}, "outside the MPI")
        row = copy.deepcopy(VALID)
        row["ranks"] = 4
        row["owner"] = 3
        self.path.write_text(record(row), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "does not match"):
            stage9.parse_log(self.path, 1)

    def test_particle_counts_and_feedback(self) -> None:
        self.assert_rejected({"daughters": 15}, "expected 16")
        self.assert_rejected({"partners": 63}, "expected 16")
        self.assert_rejected({"unchanged_records": 79}, "invariants")
        self.assert_rejected({"daughter_changes": 7}, "daughter velocity")
        self.assert_rejected({"partner_changes": 65}, "partner velocity")

    def test_trajectory_and_owner_failures(self) -> None:
        self.assert_rejected({"unfinished": 1}, "unfinished")
        self.assert_rejected({"wrong_owner": 1}, "wrong-owner")

    def test_collision_work(self) -> None:
        self.assert_rejected({"cycles": 0}, "cycle count")
        self.assert_rejected({"cycles": 17}, "cycle count")
        self.assert_rejected({"trials": 0}, "no collision work")
        self.assert_rejected({"accepted": 129}, "exceed NTC trials")
        self.assert_rejected({"updates": 63}, "too few")
        self.assert_rejected({"probability_checks": 127}, "check two")
        self.assert_rejected({"updates": 127}, "did not update both")

    def test_probability_guards(self) -> None:
        self.assert_rejected({"invalid_probabilities": 1}, "invalid or clamped")
        self.assert_rejected({"roundoff_clamps": 1}, "invalid or clamped")
        self.assert_rejected({"probability_minimum": -1.0e-6}, "outside")
        self.assert_rejected(
            {"probability_minimum": 0.9, "probability_maximum": 0.9},
            "not one",
        )

    def test_carrier_and_momentum_closure(self) -> None:
        self.assert_rejected({"correction_error": 1.0e-16}, "correction")
        self.assert_rejected({"represented_rate_error": 1.0e-16}, "rate changed")
        self.assert_rejected({"momentum_error": 3.0e-12}, "momentum")

    def test_cross_rank_invariants(self) -> None:
        baseline = self.parse()
        current = copy.deepcopy(baseline)
        current.update(ranks=4, owner=3, trials=256, accepted=100,
                       updates=200, probability_checks=200,
                       daughter_changes=15, partner_changes=50,
                       momentum_error=2.0e-15)
        stage9.compare_invariants(baseline, current, "valid")
        current["partner"] = 12
        with self.assertRaisesRegex(ValueError, "partner"):
            stage9.compare_invariants(baseline, current, "invalid")


if __name__ == "__main__":
    unittest.main(verbosity=2)

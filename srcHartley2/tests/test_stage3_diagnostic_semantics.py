#!/usr/bin/env python3
"""Pure tests for Bernoulli versus systematic uncertainty report semantics."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "srcHartley2" / "tests"))

from test_stage0_production_report import expected_production_tolerance  # noqa: E402


class DiagnosticSemanticsTest(unittest.TestCase):
    def test_bernoulli_uses_five_sigma(self) -> None:
        record = {
            "chemical_production_sigma_s-1": 3.0,
            "systematic_error_bound_s-1": None,
            "uncertainty_method": "independent-bernoulli-analytic",
        }
        self.assertEqual(
            expected_production_tolerance(
                record, 100.0, "BERNOULLI", "chemical_production_sigma_s-1"
            ),
            15.0,
        )

    def test_systematic_uses_fixed_epoch_bound(self) -> None:
        record = {
            "chemical_production_sigma_s-1": None,
            "systematic_error_bound_s-1": 7.0,
            "uncertainty_method": "fixed-carrier-epoch-bound",
        }
        self.assertEqual(
            expected_production_tolerance(
                record, 100.0, "SYSTEMATIC", "chemical_production_sigma_s-1"
            ),
            7.0,
        )

    def test_legacy_bernoulli_archive_remains_readable(self) -> None:
        record = {
            "chemical_production_sigma_s-1": 3.0,
        }
        self.assertEqual(
            expected_production_tolerance(
                record,
                100.0,
                "LEGACY_BERNOULLI",
                "chemical_production_sigma_s-1",
            ),
            15.0,
        )

    def test_systematic_rejects_bernoulli_sigma(self) -> None:
        record = {
            "represented_number_sigma_s-1": 1.0,
            "systematic_error_bound_s-1": 7.0,
            "uncertainty_method": "fixed-carrier-epoch-bound",
        }
        with self.assertRaises(AssertionError):
            expected_production_tolerance(
                record, 100.0, "SYSTEMATIC", "represented_number_sigma_s-1"
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)

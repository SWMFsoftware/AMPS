#!/usr/bin/env python3
"""Analytical statistics and corrupted-ledger fixtures; no observations."""
import csv
from pathlib import Path
import tempfile
import unittest

from analyze_diagnostics import statistics, read_blocks
from run_budget_check import verify_ledger, GAS_SIGNS


class IndependentDiagnostics(unittest.TestCase):
    def test_constant_positive_blocks_have_exact_zero_spread(self):
        result = statistics([2.0]*9, 1)
        self.assertEqual(result["lag1"], 0.0)
        self.assertEqual(result["equilibrium_failed_checks"], 0)
        self.assertEqual(result["relative_ci95_half_width"], 0.0)
        self.assertFalse(result["zero_science_support"])

    def test_linear_transient_has_analytical_lag_two_thirds(self):
        result = statistics(list(range(1, 10)), 1)
        self.assertAlmostEqual(result["lag1"], 2/3)
        self.assertFalse(result["independence_pass"])
        self.assertEqual(result["equilibrium_failed_checks"], 3)

    def test_ten_base_samples_cannot_supply_nine_doubled_batches(self):
        result = statistics([2.0]*10, 2)
        self.assertEqual(result["batch_count"], 5)
        self.assertFalse(result["eligible"])

    def test_incomplete_block_is_parser_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"blocks.csv"
            path.write_text("base_block,metric,value_si\n1,a,1\n1,b,2\n2,a,1\n")
            with self.assertRaisesRegex(ValueError, "incomplete"):
                read_blocks(path)

    def test_unexplained_inventory_is_fail_and_false_receipt_is_error(self):
        # Structural fixture values deliberately violate conservation; these
        # are parser/identity tests, not model or observational evidence.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"ledger.csv"
            base = {"checkpoint": "initial", "gas_model_particles": 0,
                    "gas_atoms": "0", "model_count_residual": 0,
                    "weighted_residual_atoms": "0", "weighted_roundoff_bound_atoms": "0"}
            for name in GAS_SIGNS:
                base[name+"_count"] = 0
                base[name+"_atoms"] = "0"
            def write(reported):
                with path.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=base)
                    writer.writeheader();writer.writerow(base)
                    writer.writerow(dict(base, checkpoint="final", gas_model_particles=1,
                                         model_count_residual=reported))
            write(1)
            self.assertEqual(verify_ledger(path)["status"], "FAIL")
            write(0)
            with self.assertRaisesRegex(ValueError, "disagrees"):
                verify_ledger(path)


if __name__ == "__main__":
    unittest.main()

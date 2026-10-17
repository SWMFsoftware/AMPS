import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from compare_extra_jet_studies import compare_family, summarize  # noqa: E402


class ExtraJetComparisonTests(unittest.TestCase):
    def test_comparison_sign_rate_and_training_summary(self):
        labels = ["jet:porco_001", "background:b"]
        fit = {"total_source_rate_s": 10.0,
               "source_shape_fractions": {"jet:porco_001": 1.0},
               "chi2": 8.0, "samples": 4}
        base = {"component_labels": labels,
                "source_rate_strategies": {
                    "shared_shape_e3_e5": {"target_fits": {"E3": fit, "E14": fit}}}}
        new_fit = {"total_source_rate_s": 12.0,
                   "source_shape_fractions": {
                       "jet:porco_001": .75,
                       "jet:e17_candidate_baghdad": .25},
                   "chi2": 6.0, "samples": 4}
        augmented = {"component_labels": [
            "jet:porco_001", "jet:e17_candidate_baghdad", "background:b"],
            "source_rate_strategies": {
                "shared_shape_e3_e5": {
                    "target_fits": {"E3": new_fit, "E14": new_fit}}}}
        rows = compare_family(base, augmented, "CO2")
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["delta_chi2_base_minus_augmented"], 2.0)
        self.assertEqual(rows[0]["candidate_source_rate_s-1"], 3.0)
        self.assertTrue(rows[0]["used_to_train_shape"])
        self.assertFalse(rows[1]["used_to_train_shape"])
        summary = summarize(rows)[0]
        self.assertEqual(summary["training_delta_objective"], .5)
        self.assertEqual(summary["all_targets_delta_sum_chi2_per_sample"], 1.0)


if __name__ == "__main__":
    unittest.main()

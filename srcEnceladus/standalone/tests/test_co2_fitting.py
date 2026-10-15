import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "co2_fitting_study", ROOT / "co2_fitting_study.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TestCO2Study(unittest.TestCase):
    def test_co2_mass_and_thermal_speed_relation(self):
        self.assertAlmostEqual(MODULE.CO2_MASS_KG / 1.66053906660e-27,
                               44.0095, places=8)
        h2o_mass = 18.01528 * 1.66053906660e-27
        ratio = np.sqrt(h2o_mass / MODULE.CO2_MASS_KG)
        self.assertLess(ratio, 1.0)
        self.assertAlmostEqual(ratio, np.sqrt(18.01528 / 44.0095), places=12)

    def test_two_epoch_ratio_is_proportional_and_bracketed(self):
        h2o = {
            "source_rate_shape_strategies": {
                name: {"target_fits": {
                    "E3": {"total_source_rate_s": 20.0},
                    "E5": {"total_source_rate_s": 60.0},
                }} for name in ("independent_each_observation",
                                "shared_shape_inms_only")
            }
        }
        co2_targets = {
            "E3": {"total_source_rate_s": 2.0},
            "E5": {"total_source_rate_s": 6.0},
        }
        co2 = {
            name: {"target_fits": co2_targets}
            for name in ("independent_each_observation", "shared_shape_e3_e5",
                         "shared_shape_all_co2")
        }
        result = MODULE.pair_h2o_co2(h2o, co2)
        for item in result.values():
            self.assertEqual(item["proportional_through_origin_apparent_ratio"], 10.0)
            self.assertEqual(item["actual_h2o_to_co2_ratio_fragmentation_bracket"],
                             [20.0, 50.0])
            self.assertIn("n=2", item["correlation_status"])


if __name__ == "__main__":
    unittest.main()

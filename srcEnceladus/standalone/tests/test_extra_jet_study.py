import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "extra_jet_e17_study", ROOT / "extra_jet_e17_study.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TestExtraJetStudy(unittest.TestCase):
    def test_small_nonnegative_fit_recovers_active_set(self):
        matrix = np.array([[1.0, 0.0, 1.0],
                           [0.0, 1.0, 1.0],
                           [1.0, 1.0, 0.0],
                           [2.0, 1.0, 1.0]])
        expected = np.array([2.0, 0.0, 0.5])
        values = matrix @ expected
        result = MODULE.small_nonnegative_fit(
            matrix, values, np.ones(len(values)))
        np.testing.assert_allclose(result["coefficients"], expected,
                                   rtol=0, atol=1.0e-12)
        self.assertLess(result["chi2"], 1.0e-24)

    def test_negative_unconstrained_amplitude_is_clamped(self):
        matrix = np.eye(2)
        result = MODULE.small_nonnegative_fit(
            matrix, np.array([3.0, -2.0]), np.ones(2))
        np.testing.assert_allclose(result["coefficients"], [3.0, 0.0])
        self.assertAlmostEqual(result["chi2"], 4.0)

    def test_angular_distance_wraps_west_longitude(self):
        distance = MODULE.angular_distance_deg(-80.0, 359.0, -80.0, 1.0)
        self.assertLess(distance, 1.0)
        self.assertGreater(distance, 0.0)

    def test_stripe_sampling_preserves_named_polylines(self):
        samples = MODULE.stripe_candidates(
            ROOT / "data" / "published_stripes_2014.json", 1.0)
        names = {item["stripe"] for item in samples}
        self.assertIn("baghdad", names)
        self.assertIn("damascus", names)
        self.assertTrue(all(-90 <= item["latitude_deg"] <= 90
                            and 0 <= item["west_longitude_deg"] < 360
                            for item in samples))


if __name__ == "__main__":
    unittest.main()

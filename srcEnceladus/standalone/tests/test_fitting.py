from pathlib import Path
import csv
import json
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fitting import solve_amplitudes  # noqa: E402
from plume_fit import load_density  # noqa: E402


class FittingTests(unittest.TestCase):
    def test_nonnegative_linear_recovery(self):
        matrix = np.array([[1., 0., 1.], [0., 1., 1.],
                           [2., 1., 1.], [1., 3., 1.], [4., 2., 1.]])
        truth = np.array([2.5, 4.0, 0.75])
        sigma = np.full(5, 0.2)
        result = solve_amplitudes(matrix, matrix@truth, sigma)
        np.testing.assert_allclose(result["coefficients"], truth, rtol=1e-10)
        self.assertFalse(result["rank_deficient"])
        self.assertEqual(result["effective_rank_by_relative_singular_cutoff"]["1e-4"], 3)
        self.assertEqual(result["fit_warnings"], [])

    def test_unobserved_column_is_flagged(self):
        matrix = np.array([[1., 0.], [2., 0.], [3., 0.]])
        result = solve_amplitudes(matrix, np.array([1., 2., 3.]), np.ones(3))
        self.assertEqual(result["unobserved_columns"].tolist(), [1])
        self.assertTrue(result["rank_deficient"])
        self.assertTrue(result["fit_warnings"])

    def test_raw_counts_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            observations = directory/"data.csv"
            with observations.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["x_m", "y_m", "z_m", "density_m3", "sigma_m3"])
                writer.writerow([1, 2, 3, 4, 1])
            metadata = directory/"metadata.json"
            metadata.write_text(json.dumps({
                "observable": "counts", "frame": "IAU_ENCELADUS",
                "position_units": "m", "density_units": "m^-3",
                "uncertainty_units": "m^-3",
                "uncertainty_model": "independent_gaussian",
                "calibration_reference": "none",
                "sticking_treatment": "none", "selection_reference": "test",
                "data_kind": "calibrated_observation",
                "time_sampling": "instantaneous"}))
            with self.assertRaises(ValueError):
                load_density(observations, metadata)

    def test_digitized_profiles_have_published_peaks(self):
        location = os.environ.get("ENCELADUS_PLUME_FIT_DATA")
        if not location and os.environ.get("ENCELADUS_VALIDATION_DATA_ROOT"):
            location = str(Path(os.environ["ENCELADUS_VALIDATION_DATA_ROOT"])
                           / "data" / "processed" / "plume_fit_v1")
        if not location:
            self.skipTest(
                "set ENCELADUS_VALIDATION_DATA_ROOT or ENCELADUS_PLUME_FIT_DATA")
        data_root = Path(location)
        for flyby, peak in (("e3", 7.6e13), ("e5", 3.4e14)):
            observations = data_root/(flyby+"_h2o_sticking_corrected.csv")
            metadata = data_root/(flyby+"_h2o_sticking_corrected.metadata.json")
            with observations.open() as stream:
                rows = list(csv.DictReader(stream))
            values = np.array([float(row["density_m3"]) for row in rows])
            self.assertEqual(len(values), 79)
            self.assertAlmostEqual(values.max()/peak, 1.0, places=9)
            self.assertTrue(np.all(values > 0))
            x, y, sigma, qualified = load_density(observations, metadata)
            self.assertEqual(x.shape, (79, 3))
            np.testing.assert_array_equal(y, values)
            self.assertTrue(np.all(sigma > 0))
            self.assertEqual(qualified["fit_role"], "absolute_density_calibration")


if __name__ == "__main__":
    unittest.main()

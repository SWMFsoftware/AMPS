import csv
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from plot_reconstructed_figures import polar_xy  # noqa: E402


class FigureReconstructionTests(unittest.TestCase):
    def test_polar_orientation_matches_paper_layout(self):
        x, y = polar_xy([-80, -80, -80, -80], [0, 90, 180, 270], shift=False)
        np.testing.assert_allclose(x, [0, 10, 0, -10], atol=1e-12)
        np.testing.assert_allclose(y, [-10, 0, 10, 0], atol=1e-12)

    def test_current_iau_overlay_shift_is_exact(self):
        x_shifted, y_shifted = polar_xy(-80, 0, shift=True)
        x_expected, y_expected = polar_xy(-80, 3.5, shift=False)
        self.assertAlmostEqual(float(x_shifted), float(x_expected), places=13)
        self.assertAlmostEqual(float(y_shifted), float(y_expected), places=13)

    def test_published_figure3_archive_has_all_six_curves(self):
        path = ROOT / "figure_data/figure3_published_curves.csv"
        with open(path, newline="") as stream:
            rows = list(csv.DictReader(stream))
        groups = {(row["flyby"], row["series"]) for row in rows}
        expected = {(flyby, series) for flyby in ("E3", "E5") for series in
                    ("INMS", "with_tiger_stripes", "without_tiger_stripes")}
        self.assertEqual(groups, expected)
        self.assertGreater(len(rows), 1000)
        for row in rows:
            self.assertGreater(float(row["density_m3"]), 0.0)
            self.assertGreaterEqual(float(row["time_ca_s"]), -401.0)
            self.assertLessEqual(float(row["time_ca_s"]), 401.0)


if __name__ == "__main__":
    unittest.main()

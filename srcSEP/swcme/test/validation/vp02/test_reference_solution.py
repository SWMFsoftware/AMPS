#!/usr/bin/env python3
"""Focused deterministic tests for the VP02 independent equations."""

import unittest

from reference_solution import fold_observed_angle_deg, parker_angle_deg


class ParkerReferenceTests(unittest.TestCase):
    def test_angle_grows_with_radius_and_falls_with_speed(self) -> None:
        self.assertGreater(parker_angle_deg(1.0, 0.0, 400.0), parker_angle_deg(0.3, 0.0, 400.0))
        self.assertLess(parker_angle_deg(1.0, 0.0, 800.0), parker_angle_deg(1.0, 0.0, 400.0))

    def test_sector_folding_is_invariant_under_field_reversal(self) -> None:
        self.assertAlmostEqual(fold_observed_angle_deg(4.0, -3.0), fold_observed_angle_deg(-4.0, 3.0))
        self.assertAlmostEqual(fold_observed_angle_deg(1.0, -1.0), 45.0)

    def test_polar_limit_is_radial(self) -> None:
        self.assertAlmostEqual(parker_angle_deg(1.0, 90.0, 400.0), 0.0, places=12)


if __name__ == "__main__":
    unittest.main(verbosity=2)

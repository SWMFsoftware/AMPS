#!/usr/bin/env python3
"""Focused tests for VP03's independent thermodynamic oracle."""

import math
import unittest

from reference_solution import background_state, leblanc_density_cm3, observed_state


class ReferenceTests(unittest.TestCase):
    def test_leblanc_is_normalized_at_one_au(self) -> None:
        self.assertAlmostEqual(leblanc_density_cm3(1.0), 5.0, places=13)
        self.assertGreater(leblanc_density_cm3(0.3), leblanc_density_cm3(0.8))

    def test_default_sound_speed_is_radius_independent(self) -> None:
        inner = background_state(0.3, 0.0)
        outer = background_state(1.0, 0.0)
        self.assertAlmostEqual(inner["sound_speed_km_s"], outer["sound_speed_km_s"], places=13)
        self.assertAlmostEqual(outer["speed_km_s"], 400.0)

    def test_observed_closure_is_positive_and_ordered(self) -> None:
        state = observed_state(5.0, 1.0e5, 5.0)
        self.assertTrue(all(math.isfinite(value) and value > 0.0 for value in state.values()))
        self.assertGreater(state["fast_speed_km_s"], state["sound_speed_km_s"])
        self.assertGreater(state["fast_speed_km_s"], state["alfven_speed_km_s"])


if __name__ == "__main__":
    unittest.main(verbosity=2)

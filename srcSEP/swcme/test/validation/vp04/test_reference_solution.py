#!/usr/bin/env python3
"""Focused mathematical tests for VP04's independent PCHIP oracle."""

import unittest

from reference_solution import evaluate


class PchipTests(unittest.TestCase):
    def test_linear_data_are_exact(self) -> None:
        times, radii = [0.0, 2.0, 5.0], [1.0, 5.0, 11.0]
        radius, derivative = evaluate(times, radii, 3.0)
        self.assertAlmostEqual(radius, 7.0, places=13)
        self.assertAlmostEqual(derivative, 2.0, places=13)

    def test_knots_are_interpolated_and_curve_is_monotone(self) -> None:
        times, radii = [0.0, 1.0, 3.0, 4.0], [2.0, 2.5, 6.0, 6.2]
        for time, expected in zip(times, radii):
            self.assertAlmostEqual(evaluate(times, radii, time)[0], expected, places=13)
        samples = [evaluate(times, radii, index / 20.0)[0] for index in range(81)]
        self.assertTrue(all(samples[index + 1] >= samples[index] for index in range(len(samples) - 1)))

    def test_outside_queries_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            evaluate([0.0, 1.0], [2.0, 3.0], -0.1)


if __name__ == "__main__":
    unittest.main(verbosity=2)

"""Stage-A2 mathematical and Python/C++ real-harmonic agreement tests."""

import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

from activity_basis import (BASIS_ID, BasisError, coefficient_count,
                            coefficient_index, evaluate_angles,
                            evaluate_cartesian, modes)


class ActivityBasisMathTests(unittest.TestCase):
    def test_identifier_count_order_and_index(self):
        self.assertEqual(BASIS_ID, "liouville_cops.real_orthonormal_cs_v1")
        self.assertEqual([coefficient_count(5, degree) for degree in range(5, 9)],
                         [11, 24, 39, 56])
        ordered = modes(5, 8)
        self.assertEqual(ordered[0], (5, -5))
        self.assertEqual(ordered[10], (5, 5))
        self.assertEqual(ordered[11], (6, -6))
        self.assertEqual(ordered[-1], (8, 8))
        for index, (degree, order) in enumerate(ordered):
            self.assertEqual(coefficient_index(degree, order, 5, 8), index)

    def test_known_normalization_phase_and_poles(self):
        root = math.sqrt
        y00 = evaluate_cartesian([1.0, 0.0, 0.0], 0, 1)
        self.assertAlmostEqual(y00[coefficient_index(0, 0, 0, 1)],
                               1.0 / root(4.0 * math.pi), places=15)
        north = evaluate_cartesian([0.0, 0.0, 2.0], 0, 1)
        south = evaluate_cartesian([0.0, 0.0, -3.0], 0, 1)
        self.assertAlmostEqual(north[coefficient_index(1, 0, 0, 1)],
                               root(3.0 / (4.0 * math.pi)), places=15)
        self.assertAlmostEqual(south[coefficient_index(1, 0, 0, 1)],
                               -root(3.0 / (4.0 * math.pi)), places=15)
        for order in (-1, 1):
            self.assertEqual(north[coefficient_index(1, order, 0, 1)], 0.0)
            self.assertEqual(south[coefficient_index(1, order, 0, 1)], 0.0)
        plus_x = evaluate_cartesian([1.0, 0.0, 0.0], 1, 1)
        plus_y = evaluate_cartesian([0.0, 1.0, 0.0], 1, 1)
        expected = -root(3.0 / (4.0 * math.pi))
        self.assertAlmostEqual(plus_x[coefficient_index(1, 1, 1, 1)],
                               expected, places=15)
        self.assertAlmostEqual(plus_y[coefficient_index(1, -1, 1, 1)],
                               expected, places=15)

    def test_scale_invariance_and_scalar_table_shapes(self):
        point = np.array([0.31, -0.47, 0.82])
        one = evaluate_cartesian(point, 5, 8)
        table = evaluate_cartesian([point, 1e6 * point, 1e-6 * point], 5, 8)
        self.assertEqual(one.shape, (56,))
        self.assertEqual(table.shape, (3, 56))
        np.testing.assert_allclose(table, np.tile(one, (3, 1)), rtol=2e-15,
                                   atol=2e-15)

    def test_orthonormality_through_degree_eight(self):
        # Gauss-Legendre integrates mu polynomials through degree 2n-1; 20
        # nodes exceed the degree-16 products here.  A 36-point periodic
        # longitude grid exactly resolves all Fourier products through |m|=8.
        mu, weight_mu = np.polynomial.legendre.leggauss(20)
        phi = 2.0 * math.pi * np.arange(36) / 36.0
        theta, longitude = np.meshgrid(np.arccos(mu), phi, indexing="ij")
        values = evaluate_angles(theta, longitude, 0, 8).reshape(-1, 81)
        weights = np.repeat(weight_mu, len(phi)) * (2.0 * math.pi / len(phi))
        gram = values.T @ (weights[:, None] * values)
        error = np.max(np.abs(gram - np.eye(81)))
        self.assertLess(error, 2e-13)

    def test_invalid_inputs(self):
        for bounds in ((-1, 8), (5, 4), (0, 9)):
            with self.assertRaises(BasisError):
                coefficient_count(*bounds)
        for point in ([0, 0, 0], [math.nan, 0, 1], [math.inf, 0, 1], [1, 2]):
            with self.assertRaises(BasisError):
                evaluate_cartesian(point, 5, 8)
        with self.assertRaises(BasisError):
            coefficient_index(4, 0, 5, 8)
        with self.assertRaises(BasisError):
            evaluate_angles(-0.1, 0.0, 5, 8)


class ActivityBasisCrossLanguageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("g++ is unavailable")
        cls.temp = tempfile.TemporaryDirectory(prefix="cops_a2_basis_")
        cls.probe = Path(cls.temp.name) / "activity_basis_probe"
        command = [compiler, "-std=c++17", "-O2", "-Wall", "-Wextra",
                   "-pedantic", "-I", str(APP / "src"),
                   str(APP / "tests/activity_basis_probe.cpp"),
                   str(APP / "src/ActivityEnhancementBasis.cpp"),
                   "-o", str(cls.probe)]
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError("C++ basis probe failed to compile:\n" + completed.stdout)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "temp"):
            cls.temp.cleanup()

    def cpp_values(self, points, min_degree, max_degree):
        text = "".join("%.17g %.17g %.17g\n" % tuple(point) for point in points)
        completed = subprocess.run(
            [str(self.probe), str(min_degree), str(max_degree)], input=text,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return np.array([[float(value) for value in line.split()]
                         for line in completed.stdout.splitlines()])

    def test_python_cpp_all_modes_reference_coordinates(self):
        points = np.array([
            [0.0, 0.0, 1.0], [0.0, 0.0, -2.0],
            [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
            [0.31, -0.47, 0.82], [-2.3, 0.17, -0.91],
            [1e-12, -2e-12, 1.0], [1e-12, 2e-12, -1.0],
        ])
        for minimum, maximum in ((0, 8), (5, 5), (5, 6), (5, 7), (5, 8)):
            python = evaluate_cartesian(points, minimum, maximum)
            cpp = self.cpp_values(points, minimum, maximum)
            difference = np.max(np.abs(python - cpp))
            self.assertLessEqual(difference, 2e-14,
                                 "range %d..%d differs by %.3e" %
                                 (minimum, maximum, difference))

    def test_cpp_rejects_zero_vector_and_bad_degrees(self):
        for arguments, text in ((["5", "8"], "0 0 0\n"),
                                (["5", "4"], "1 0 0\n"),
                                (["0", "9"], "1 0 0\n")):
            completed = subprocess.run([str(self.probe)] + arguments, input=text,
                                       text=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Error:", completed.stderr)


if __name__ == "__main__":
    unittest.main()


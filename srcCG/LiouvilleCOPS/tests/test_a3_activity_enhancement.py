"""Stage-A3 strict-file, residual, clipping, and compatibility tests."""

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

from activity_basis import coefficient_index, evaluate_cartesian


def enhancement_text(max_degree=5, coefficients=None, **replace):
    """Create one complete canonical file, allowing one-field mutations."""
    if coefficients is None:
        coefficients = {(degree, order): 0.0
                        for degree in range(5, max_degree + 1)
                        for order in range(-degree, degree + 1)}
    fields = {
        "Format": "liouville_cops.activity_enhancement",
        "FormatVersion": "1",
        "BasisID": "liouville_cops.real_orthonormal_cs_v1",
        "Species": "H2O",
        "BodyFrame": "67P/C-G_CK",
        "AngleConvention": "theta=acos(z/r);phi=atan2(y,x)",
        "Ordering": "increasing_l_then_m_minus_l_to_plus_l",
        "Clipping": "max(0,legacy_linear+enhancement)",
        "MinimumDegree": "5",
        "MaximumDegree": str(max_degree),
        "Scale": "2.0",
        "ScaleInterpretation": "dimensionless_multiplier_of_normalized_basis",
        "LegacyActivityIdentity": "builtin:shap5",
        "Regularization": "not_applicable",
        "Lambda": "0.0",
    }
    fields.update({key: str(value) for key, value in replace.items()})
    lines = ["%s = %s" % item for item in fields.items()]
    lines.extend("Coefficient = %d %d %.17g" % (degree, order, value)
                 for (degree, order), value in coefficients.items())
    return "\n".join(lines) + "\n"


class ActivityEnhancementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("g++ is unavailable")
        cls.temp = tempfile.TemporaryDirectory(prefix="cops_a3_enhancement_")
        cls.directory = Path(cls.temp.name)
        cls.probe = cls.directory / "activity_enhancement_probe"
        command = [compiler, "-std=c++17", "-O2", "-Wall", "-Wextra",
                   "-pedantic", "-I", str(APP / "src"),
                   str(APP / "tests/activity_enhancement_probe.cpp"),
                   str(APP / "src/ActivityEnhancement.cpp"),
                   str(APP / "src/ActivityEnhancementBasis.cpp"),
                   "-o", str(cls.probe)]
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError("C++ enhancement probe failed to compile:\n" + completed.stdout)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_file(self, text, point=(0.31, -0.47, 0.82), legacy=0.7,
                 frame="67P/C-G_CK", identity="builtin:shap5"):
        path = self.directory / ("case_%d.txt" % abs(hash(text)))
        path.write_text(text, encoding="utf-8")
        return subprocess.run(
            [str(self.probe), str(path), frame, identity,
             *("%.17g" % value for value in point), "%.17g" % legacy],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def values(self, completed):
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return [float(value) for value in completed.stdout.split()]

    def test_zero_coefficients_and_degree_four_control_reproduce_baseline(self):
        for maximum in (4, 5, 8):
            result = self.values(self.run_file(enhancement_text(maximum), legacy=0.7))
            self.assertEqual(int(result[0]), maximum)
            self.assertEqual(int(result[1]), 0 if maximum == 4 else (maximum + 1) ** 2 - 25)
            self.assertEqual(result[2], 0.0)
            self.assertEqual(result[3], 0.7)
            self.assertEqual(result[4], 0.7)
        clipped = self.values(self.run_file(enhancement_text(5), legacy=-0.25))
        self.assertEqual(clipped[3], -0.25)
        self.assertEqual(clipped[4], 0.0)

    def test_single_mode_matches_normative_python_basis_and_clips_once(self):
        point = np.array([0.31, -0.47, 0.82])
        coefficient = 1.25
        coefficients = {(5, order): 0.0 for order in range(-5, 6)}
        coefficients[(5, -2)] = coefficient
        completed = self.run_file(enhancement_text(5, coefficients), point=point,
                                  legacy=-0.1)
        maximum, count, residual, total, applied = self.values(completed)
        expected = 2.0 * coefficient * evaluate_cartesian(point, 5, 5)[
            coefficient_index(5, -2, 5, 5)]
        self.assertEqual((int(maximum), int(count)), (5, 11))
        self.assertAlmostEqual(residual, expected, places=14)
        self.assertAlmostEqual(total, -0.1 + expected, places=14)
        self.assertEqual(applied, max(0.0, total))

    def test_rejects_wrong_frame_and_legacy_identity(self):
        text = enhancement_text(5)
        for frame, identity, fragment in (
                ("WRONG_FRAME", "builtin:shap5", "BodyFrame"),
                ("67P/C-G_CK", "builtin:shap5_1", "LegacyActivityIdentity")):
            completed = self.run_file(text, frame=frame, identity=identity)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn(fragment, completed.stderr)

    def test_rejects_incompatible_metadata_and_nonfinite_values(self):
        mutations = (
            ({"Species": "CO2"}, "Species"),
            ({"BasisID": "wrong"}, "BasisID"),
            ({"Ordering": "wrong"}, "Ordering"),
            ({"Clipping": "max(0,legacy)+enhancement"}, "Clipping"),
            ({"Scale": "nan"}, "Scale"),
            ({"Scale": "0"}, "Scale"),
            ({"Lambda": "-1"}, "Lambda"),
            ({"MaximumDegree": "9"}, "MaximumDegree"),
        )
        for mutation, fragment in mutations:
            with self.subTest(mutation=mutation):
                completed = self.run_file(enhancement_text(5, **mutation))
                self.assertNotEqual(completed.returncode, 0)
                self.assertIn(fragment, completed.stderr)

    def test_rejects_missing_duplicate_unexpected_and_unknown_content(self):
        base = enhancement_text(5)
        cases = (
            (base.replace("Coefficient = 5 -5 0\n", ""), "missing coefficient"),
            (base + "Coefficient = 5 -5 0\n", "duplicate coefficient"),
            (base + "Coefficient = 6 0 0\n", "unexpected coefficient"),
            (base + "Mystery = value\n", "unknown field"),
            (base + "Scale = 2\n", "duplicate field"),
        )
        for text, fragment in cases:
            with self.subTest(fragment=fragment):
                completed = self.run_file(text)
                self.assertNotEqual(completed.returncode, 0)
                self.assertIn(fragment, completed.stderr)


if __name__ == "__main__":
    unittest.main()

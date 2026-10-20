"""Mesh identity, template basis, bounds, and exact C++ replay tests."""

from argparse import Namespace
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

import fit_cops
from activity_templates import load_bundle, select_templates, write_bundle
from gravity_cache import (derive_effective_fields, load_cache, make_definition,
                           write_cache)


class ActivityTemplateTests(unittest.TestCase):
    def test_bundle_is_mesh_bound_and_pickle_free(self):
        centroid = np.array([[1.0, 2.0, 3.0], [-2.0, 1.0, 0.5]])
        area = np.array([4.0, 5.0])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "templates.npz"
            write_bundle(path, ["ng_only", "cliff"],
                         np.array([[1.0, 0.0], [0.0, 1.0]]),
                         centroid, area, {"purpose": "unit_test"})
            bundle = load_bundle(path, centroid, area)
            np.testing.assert_array_equal(
                select_templates(bundle, ["cliff", "ng_only"]),
                np.array([[0.0, 1.0], [1.0, 0.0]]))
            changed = centroid.copy()
            changed[0, 0] = np.nextafter(changed[0, 0], np.inf)
            with self.assertRaisesRegex(ValueError, "mesh identity"):
                load_bundle(path, changed, area)

    def test_template_parameter_is_nonnegative_and_in_combined_jacobian(self):
        # Minimal positive synthetic model using the same sparse evaluator as
        # the production fit.  Degree 4 leaves only the localized parameter.
        nface, npoint = 6, 3
        model = fit_cops.Model.__new__(fit_cops.Model)
        model.nfaces, model.npt = nface, npoint
        model.Y = np.zeros((nface, 25)); model.Y[:, 0] = 1.0
        model.face = np.arange(nface); model.pt = np.arange(nface) % npoint
        order = np.argsort(model.pt, kind="stable")
        model.face, model.pt = model.face[order], model.pt[order]
        model.indptr = np.searchsorted(model.pt, np.arange(npoint + 1))
        model.lit = np.ones(nface, dtype=bool); model.xlit = np.zeros(nface)
        model.bm = {species: 1.0 for species in fit_cops.SPECIES}
        model.cNG = {species: 1.0 for species in fit_cops.SPECIES}
        model.cRG = {species: np.ones(npoint) for species in fit_cops.SPECIES}
        model.wn = {species: np.ones(nface) for species in fit_cops.SPECIES}
        model.wv = {species: np.ones(nface) for species in fit_cops.SPECIES}
        model.wp = {"H2O": np.ones(nface)}
        model.PhiCO2 = np.ones(npoint)
        model.T0 = np.array([180.0, 150.0, 120.0])
        model.rows = np.arange(npoint)
        model.useNG = model.useRG = np.ones(npoint, dtype=int)
        coef = np.zeros((3, 25)); coef[0, 0] = 2.0; coef[1, 0] = 0.1
        base_ng, base_rg = model.evaluate(coef[0], model.T0, 1.0)
        model.refNG, model.refRG = base_ng, base_rg
        elements = {"centroid": np.column_stack((np.ones(nface),
                                                  np.arange(nface),
                                                  np.ones(nface))),
                    "area": np.ones(nface), "Y": model.Y}
        template = np.array([[1.0], [1.0], [0.0], [0.0], [0.0], [0.0]])
        args = Namespace(enhancement_degree=4, fix_legacy_activity=True,
                         enhancement_scale=None, degree=4,
                         ordered_temperature=False, fix_temperature=True,
                         fix_kappa=True, tbounds=(30.0, 400.0),
                         body_frame="67P/C-G_CK",
                         activity_template_names=["ng_only"],
                         activity_template_sources=[{"name": "ng_only"}],
                         activity_template_basis=template)
        problem = fit_cops.Problem(model, elements, coef, args)
        self.assertEqual((problem.nh, problem.nt, problem.nb), (0, 1, 1))
        self.assertEqual(problem.names, ["c_ng_only"])
        self.assertEqual(problem.lo[0], 0.0)
        self.assertTrue(np.isposinf(problem.hi[0]))
        self.assertEqual(problem.data_jacobian(problem.x0).shape,
                         (problem.ndata, 1))


class GravityCacheTests(unittest.TestCase):
    """The expensive self-gravity field is reusable only on the exact mesh."""

    def setUp(self):
        self.centroid = np.array([[2.0, 0.0, 0.0], [0.0, 3.0, 0.5]])
        self.normal = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        self.area = np.array([4.0, 5.0])
        self.gravity = np.array([[-1.0e-4, 0.0, 0.0],
                                 [0.0, -2.0e-4, -1.0e-5]])
        self.omega = np.array([0.0, 0.0, 1.0e-4])
        self.definition = make_definition(
            self.centroid, self.normal, self.area,
            {"size_bytes": 123, "sha256": "a" * 64}, 8, 12, 533.0, 0.1)

    def test_round_trip_and_new_omega_reuse_self_gravity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gravity.npz"
            written = write_cache(
                path, self.centroid, self.normal, self.area, self.definition,
                self.gravity, self.omega,
                {"offset_halving_slope_difference_deg_max": 0.01})
            loaded = load_cache(
                path, self.centroid, self.normal, self.area, self.definition,
                self.omega)
            self.assertTrue(loaded["reference_omega_reused"])
            np.testing.assert_array_equal(loaded["self_gravity"], self.gravity)
            np.testing.assert_array_equal(loaded["slope"], written["slope"])

            new_omega = self.omega * 1.1
            changed = load_cache(
                path, self.centroid, self.normal, self.area, self.definition,
                new_omega)
            self.assertFalse(changed["reference_omega_reused"])
            expected = derive_effective_fields(
                self.centroid, self.normal, self.gravity, new_omega)
            np.testing.assert_allclose(changed["centrifugal"], expected[0])
            np.testing.assert_allclose(changed["effective"], expected[1])
            np.testing.assert_allclose(changed["slope"], expected[2])

    def test_rejects_changed_mesh_or_physical_input(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gravity.npz"
            write_cache(path, self.centroid, self.normal, self.area,
                        self.definition, self.gravity, self.omega, {})
            changed_normal = self.normal.copy()
            changed_normal[0, 0] = np.nextafter(1.0, 0.0)
            with self.assertRaisesRegex(ValueError, "geometry identity"):
                load_cache(path, self.centroid, changed_normal, self.area,
                           self.definition, self.omega)
            changed_definition = dict(self.definition, density_kg_m3=500.0)
            with self.assertRaisesRegex(ValueError, "physical inputs"):
                load_cache(path, self.centroid, self.normal, self.area,
                           changed_definition, self.omega)


class FacetReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("g++ is unavailable")
        cls.temp = tempfile.TemporaryDirectory(prefix="cops_facet_replay_")
        cls.directory = Path(cls.temp.name)
        cls.probe = cls.directory / "probe"
        command = [compiler, "-std=c++17", "-O2", "-I", str(APP / "src"),
                   str(APP / "tests/activity_facet_enhancement_probe.cpp"),
                   str(APP / "src/ActivityEnhancement.cpp"),
                   str(APP / "src/ActivityEnhancementBasis.cpp"),
                   "-o", str(cls.probe)]
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError("facet replay probe failed to compile:\n" + completed.stdout)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def text(self, x="1", residual="0.35"):
        return "\n".join((
            "Format = liouville_cops.activity_enhancement_facets",
            "FormatVersion = 1",
            "BasisID = real_spherical_harmonics_plus_declared_facet_templates_v1",
            "Species = H2O", "BodyFrame = 67P/C-G_CK",
            "Ordering = increasing_face_index",
            "Clipping = max(0,legacy_linear+enhancement)",
            "HarmonicMinimumDegree = 5", "HarmonicMaximumDegree = 8",
            "Scale = 1",
            "ScaleInterpretation = explicit_dimensionless_residual_per_facet",
            "LegacyActivityIdentity = builtin:shap5",
            "Regularization = unit_test", "Lambda = 0", "FacetCount = 1",
            "TemplateCoefficient = cliff 0.5",
            "Facet = 0 %s 2 3 4 %s" % (x, residual), ""))

    def run_text(self, text):
        path = self.directory / (hashlib.sha256(text.encode()).hexdigest() + ".txt")
        path.write_text(text)
        return subprocess.run([str(self.probe), str(path), "67P/C-G_CK",
                               "builtin:shap5", "1", "2", "3", "4", "-0.1"],
                              text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

    def test_exact_field_replay_and_single_clip(self):
        completed = self.run_text(self.text())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        count, residual, total, applied = map(float, completed.stdout.split())
        self.assertEqual(count, 1.0)
        self.assertAlmostEqual(residual, 0.35)
        self.assertAlmostEqual(total, 0.25)
        self.assertAlmostEqual(applied, 0.25)

    def test_rejects_different_mesh_and_negative_template_amplitude(self):
        for text, fragment in ((self.text(x="1.0000000000000002"), "active mesh"),
                               (self.text() + "TemplateCoefficient = bad -1\n",
                                "nonnegative")):
            completed = self.run_text(text)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn(fragment, completed.stderr)


if __name__ == "__main__":
    unittest.main()

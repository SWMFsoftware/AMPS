"""Stage-A4 fitter tests that do not require mission kernels or SPICE."""

from argparse import Namespace
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

import fit_cops
from activity_basis import coefficient_count, evaluate_cartesian


def synthetic_problem(max_degree):
    """Build a deterministic positive-field pressure problem in memory."""
    rng = np.random.default_rng(92341)
    nface, npoint = 80, 8
    centroid = rng.normal(size=(nface, 3))
    centroid /= np.linalg.norm(centroid, axis=1)[:, None]
    legacy_basis = rng.normal(scale=0.05, size=(nface, 25))
    legacy_basis[:, 0] = 1.0
    area = rng.uniform(0.3, 2.0, nface)
    coef = np.zeros((3, 25))
    coef[0, 0] = 1.7
    coef[1, 0] = 0.2

    # One contribution from every face.  The numerical values need only be
    # positive and non-degenerate; Model.evaluate supplies the real production
    # algebra and sparse accumulation used by the mission fit.
    model = fit_cops.Model.__new__(fit_cops.Model)
    model.nfaces, model.npt = nface, npoint
    model.Y = legacy_basis
    model.face = np.arange(nface)
    model.pt = np.arange(nface) % npoint
    order = np.argsort(model.pt, kind="stable")
    model.face, model.pt = model.face[order], model.pt[order]
    model.indptr = np.searchsorted(model.pt, np.arange(npoint + 1))
    ncontrib = len(model.face)
    model.lit = np.ones(ncontrib, dtype=bool)
    model.xlit = rng.uniform(0.0, 0.4, ncontrib)
    model.bm = {species: 0.7 + 0.2 * i
                for i, species in enumerate(fit_cops.SPECIES)}
    model.cNG = {species: 0.8 + 0.2 * i
                 for i, species in enumerate(fit_cops.SPECIES)}
    model.cRG = {species: rng.uniform(0.7, 1.3, npoint)
                 for species in fit_cops.SPECIES}
    model.wn = {species: rng.uniform(0.2, 1.0, ncontrib)
                for species in fit_cops.SPECIES}
    model.wv = {species: rng.uniform(0.2, 1.0, ncontrib)
                for species in fit_cops.SPECIES}
    model.wp = {"H2O": rng.uniform(0.2, 1.0, ncontrib)}
    model.PhiCO2 = rng.uniform(0.2, 1.0, npoint)
    model.T0 = np.array([205.0, 175.0, 95.0])
    model.rows = np.arange(npoint)
    model.useNG = np.ones(npoint, dtype=int)
    model.useRG = np.ones(npoint, dtype=int)

    enh = (evaluate_cartesian(centroid, 5, max_degree)
           if max_degree >= 5 else np.empty((nface, 0)))
    scale = np.sqrt(np.sum(area * (legacy_basis @ coef[0]) ** 2) /
                    np.sum(area))
    base_ng, base_rg = model.evaluate(coef[0], model.T0, 1.0)
    model.refNG = base_ng * np.linspace(0.93, 1.07, npoint)
    model.refRG = base_rg * np.linspace(1.06, 0.94, npoint)

    elements = {"centroid": centroid, "area": area, "Y": legacy_basis}
    args = Namespace(enhancement_degree=max_degree,
                     fix_legacy_activity=True,
                     enhancement_scale=scale,
                     degree=4, ordered_temperature=False,
                     fix_temperature=True, fix_kappa=True,
                     tbounds=(30.0, 400.0), body_frame="67P/C-G_CK")
    return fit_cops.Problem(model, elements, coef, args), elements, coef


class EnhancementFitTests(unittest.TestCase):
    def test_parameter_counts_scale_and_degree_four_control(self):
        for degree, expected in ((4, 0), (5, 11), (6, 24), (7, 39), (8, 56)):
            problem, elements, coef = synthetic_problem(degree)
            self.assertEqual(problem.nb, expected)
            self.assertEqual(problem.nb,
                             0 if degree == 4 else coefficient_count(5, degree))
            expected_scale = np.sqrt(
                np.sum(elements["area"] * (elements["Y"] @ coef[0]) ** 2) /
                np.sum(elements["area"]))
            self.assertEqual(problem.enhancement_scale, expected_scale)
            if degree == 4:
                self.assertEqual(len(problem.x0), 0)
                np.testing.assert_array_equal(problem.data_jacobian(problem.x0),
                                              np.empty((problem.ndata, 0)))
                self.assertEqual(problem.deviation(problem.x0), 0.0)

    def test_analytic_data_and_regularization_jacobian(self):
        problem, _, _ = synthetic_problem(6)
        problem.lam = 3.2e-4
        x = np.linspace(-0.03, 0.04, len(problem.x0))
        analytic = problem.jacobian(x)
        for column in range(len(x)):
            step = 1.0e-6
            plus, minus = x.copy(), x.copy()
            plus[column] += step
            minus[column] -= step
            finite_difference = ((problem.residuals(plus) -
                                  problem.residuals(minus)) / (2.0 * step))
            np.testing.assert_allclose(analytic[:, column], finite_difference,
                                       rtol=3e-6, atol=2e-9)

    def test_regularization_is_unclipped_area_weighted_residual(self):
        problem, elements, _ = synthetic_problem(5)
        problem.lam = 1.0
        b = np.linspace(-0.2, 0.15, problem.nb)
        delta = problem.enhancement_scale * (problem.enhancement_basis @ b)
        expected = np.sqrt(np.sum(elements["area"] * delta ** 2) / problem.N0)
        self.assertAlmostEqual(np.linalg.norm(problem.regularization(b)), expected,
                               places=13)

    def test_one_sided_derivative_uses_inactive_branch_near_clip(self):
        """At a clip kink, approach from the declared inactive side."""
        problem, _, _ = synthetic_problem(5)
        column = int(np.argmax(np.abs(problem.enhancement_basis[0])))
        derivative = (problem.enhancement_scale *
                      problem.enhancement_basis[0, column])
        x = problem.x0.copy()
        # Place facet zero just below the clip, then perturb in the direction
        # that makes it more negative.  The analytic convention is dA/db=0
        # for total_linear <= 0, so this is its meaningful one-sided check.
        target = -1.0e-10
        x[column] = (target - problem.legacy_linear[0]) / derivative
        total, field_derivative = problem.activity_linear_and_derivative(x, column)
        self.assertLessEqual(total[0], 0.0)
        direction = -1.0 if field_derivative[0] > 0 else 1.0
        step = 1.0e-8
        active = total > 0.0
        while np.any((total + direction * step * field_derivative > 0.0) != active):
            step *= 0.25
            self.assertGreater(step, 1.0e-14)
        trial = x.copy()
        trial[column] += direction * step
        finite_difference = direction * (
            problem.residuals(trial) - problem.residuals(x)) / step
        np.testing.assert_allclose(problem.jacobian(x)[:, column],
                                   finite_difference, rtol=2e-6, atol=2e-8)

    def test_written_exchange_file_binds_legacy_hash_and_all_modes(self):
        problem, elements, _ = synthetic_problem(5)
        b = np.linspace(-0.1, 0.1, problem.nb)
        with tempfile.TemporaryDirectory() as directory:
            legacy = Path(directory) / "ActivityFit.txt"
            legacy.write_text("legacy snapshot\n")
            fit_cops.write_enhancement_products(
                directory, problem, b, elements,
                np.arange(len(elements["area"])) % 2 == 0,
                1e-4, legacy)
            text = (Path(directory) / "ActivityEnhancementFit.txt").read_text()
            digest = hashlib.sha256(legacy.read_bytes()).hexdigest()
            self.assertIn("LegacyActivityIdentity = sha256:" + digest, text)
            self.assertEqual(text.count("Coefficient ="), 11)
            facets = (Path(directory) / "ActivityEnhancementFacets.dat").read_text()
            self.assertIn("seen legacy_linear enhancement_linear", facets)


if __name__ == "__main__":
    unittest.main()

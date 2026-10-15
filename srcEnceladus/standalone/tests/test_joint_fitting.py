import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "joint_fitting_study", ROOT/"joint_fitting_study.py")
JOINT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(JOINT)


class JointFittingTests(unittest.TestCase):
    def test_signed_nuisance_with_more_sources_than_samples(self):
        coordinate = np.linspace(-1.0, 1.0, 7)
        sources = np.column_stack([
            np.exp(-((coordinate-center)/0.35)**2)
            for center in np.linspace(-1.2, 1.2, 10)])
        matrix = np.column_stack((sources, np.ones(len(coordinate))))
        truth = np.zeros(matrix.shape[1])
        truth[[1, 7, -1]] = [2.0, 3.0, -0.4]
        values = matrix@truth
        sigma = np.linspace(0.5, 1.1, len(values))

        fit = JOINT.bounded_weighted_fit(
            matrix, values, sigma, signed_last=True)

        self.assertEqual(
            fit["solver"],
            "projected_nnls_with_analytic_signed_nuisance")
        self.assertTrue(np.all(fit["coefficient"][:-1] >= 0))
        np.testing.assert_allclose(fit["prediction"], values,
                                   rtol=0, atol=2.0e-14)

    def test_signed_nuisance_must_be_observable(self):
        with self.assertRaisesRegex(ValueError, "finite and nonzero"):
            JOINT.bounded_weighted_fit(
                np.zeros((4, 3)), np.ones(4), np.ones(4), signed_last=True)

    def test_zero_residual_dof_has_no_reduced_chi_square(self):
        matrix = np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        fit = JOINT.fit_uvis(
            matrix, np.asarray([0.8, 1.8, 2.8]), np.ones(3),
            1.0e20, 1.0e19, ["source:1", "source:2"])
        self.assertEqual(fit["shape_dof"], 0)
        self.assertIsNone(fit["shape_reduced_chi2"])

    def test_combined_fit_recovers_one_vector_and_balances_datasets(self):
        truth = np.asarray([2.0, 3.0])
        first = np.asarray([[1.0, 0.0], [0.0, 1.0]])
        second = np.asarray([
            [1.0, 1.0], [2.0, 1.0], [1.0, 2.0], [3.0, 1.0]])
        fit = JOINT.solve_combined_blocks([
            ("short", first, first@truth, np.ones(2)),
            ("long", second, second@truth, np.ones(4)),
        ])

        np.testing.assert_allclose(fit["coefficients"], truth,
                                   rtol=0, atol=2.0e-14)
        self.assertEqual(fit["datasets"], ["short", "long"])
        self.assertAlmostEqual(
            fit["dataset_weighting"]["short"]["sigma_multiplier"],
            np.sqrt(2.0))
        self.assertAlmostEqual(
            fit["dataset_weighting"]["long"]["sigma_multiplier"], 2.0)
        self.assertTrue(fit["balance_by_dataset"])
        self.assertAlmostEqual(
            fit["per_dataset"]["short"]["chi2_per_sample"], 0.0,
            places=25)
        self.assertAlmostEqual(
            fit["per_dataset"]["long"]["chi2_per_sample"], 0.0,
            places=25)
        self.assertAlmostEqual(
            fit["chi2"], sum(item["balanced_objective_contribution"]
                             for item in fit["per_dataset"].values()),
            places=25)

    def test_fixed_shape_fit_recovers_total_rate_and_background(self):
        labels = ["jet:a", "stripe:b", "background:c0_m-3"]
        shape_vector = np.asarray([2.0, 6.0, 999.0])
        matrix = np.asarray([
            [1.0, 0.0, 1.0],
            [0.0, 2.0, 1.0],
            [1.0, 1.0, 1.0],
            [2.0, 0.5, 1.0],
        ])
        normalized = np.asarray([0.25, 0.75, 0.0])
        truth_total = 4.5
        truth_background = 0.7
        values = matrix@(truth_total*normalized+np.asarray([0.0, 0.0, truth_background]))

        fit = JOINT.fit_rate_for_fixed_source_shape(
            matrix, values, np.ones(len(values)), labels, shape_vector)

        self.assertAlmostEqual(fit["total_source_rate_s"]/1.0e28, 4.5)
        self.assertAlmostEqual(fit["background_coefficients"]["background:c0_m-3"], 0.7)
        self.assertAlmostEqual(sum(fit["source_shape_fractions"].values()), 1.0)
        np.testing.assert_allclose(
            fit["full_coefficients"], [1.125, 3.375, 0.7], atol=1.0e-12)

    def test_fixed_shape_rejects_zero_source(self):
        with self.assertRaisesRegex(ValueError, "zero total source"):
            JOINT.fit_rate_for_fixed_source_shape(
                np.eye(2), np.ones(2), np.ones(2),
                ["jet:a", "background:c0_m-3"], np.zeros(2))


if __name__ == "__main__":
    unittest.main()

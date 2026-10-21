"""Joint baseline/template tests using the actual sparse pressure evaluator.

The synthetic transport geometry avoids mission/SPICE dependencies but retains
temperature dependence, both gauges, CO2, clipping, and analytic derivatives.
"""
from copy import deepcopy
import unittest
from unittest.mock import patch
from pathlib import Path
import tempfile
import shutil
import subprocess
import numpy as np
from test_a4_enhancement_fit import synthetic_problem
import fit_cops


def joint_problem():
    base, elements, coef = synthetic_problem(4)
    args = deepcopy(base.args)
    args.fix_legacy_activity = args.fix_temperature = args.fix_kappa = False
    args.activity_template_names = ["ng_only", "cliff"]
    xyz = elements["centroid"]
    args.activity_template_basis = np.column_stack((xyz[:, 0] > 0.4,
                                                   xyz[:, 1] < -0.2)).astype(float)
    # Include both illuminated and night facets to exercise Tmin derivatives.
    base.M.lit[::3] = False
    return fit_cops.Problem(base.M, elements, coef, args), elements, coef


class JointTemplateTests(unittest.TestCase):
    def test_complete_joint_output_uses_fitted_baseline(self):
        p, elements, coef = joint_problem()
        p.args.profile = "all_available"
        p.args.data = None
        p.saturated = np.zeros(p.M.npt, dtype=bool)
        p.selection_reasons = {g: ["selected"]*p.M.npt for g in ("NG", "RG")}
        p.transition_intervals = []
        p.M.helio = "analytical"
        x = p.x0.copy(); x[0] += 0.2; x[p.bslice] = [0.1, 0.4]
        a, b, _, _ = p.unpack(x)
        result = dict(x=x, lam=0., msNG=.1, msRG=.1, cost=.1,
                      D=p.deviation(x), peff=3., message="synthetic", status=1,
                      nfev=1, optimality=0., sv=np.ones(31), sv_data=np.ones(16),
                      rank_data=16, rank_augmented=16, sig=np.zeros(31), success=True)
        names = ("iPoint UTC altitude inside pNG pRG refNG refRG useNG useRG "
                 "Q_H2O Q_CO2 n_H2O n_CO2 PhiRG_H2O PhiRG_CO2 pNG_H2O "
                 "pNG_CO2 pRG_H2O pRG_CO2 NGBelowDeckFraction").split()
        cols = {n:i for i,n in enumerate(names)}
        rows = []
        for i in range(p.M.npt):
            row = ["1"]*len(names)
            row[0] = str(i); row[1] = "2016-09-30T10:35:%02d" % i
            row[cols["inside"]] = "0"; rows.append(row)
        manifest = dict(inputs={"housekeeping_table":{"sha256":"a"*64}},
                        selection={"sha256":"b"*64},
                        prior={"activity_coefficients_sha256":"c"*64},
                        objective={"selection_policy":{},
                                   "regularization":{"identifier":fit_cops.JOINT_REGULARIZATION}})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/"input.snapshot.in").write_text("ActivityModel = shap5\n")
            with patch.multiple(fit_cops, P=p, MF=p.M, ACTIVE_MANIFEST=manifest):
                fit_cops.write_fit(directory, result, directory, ["synthetic"],
                                   names, rows, cols, directory, elements,
                                   np.ones(len(elements["area"]), dtype=bool), .2)
            np.testing.assert_allclose(np.loadtxt(root/"ActivityFit.txt")[0], a,
                                       atol=0, rtol=0)
            fields = np.loadtxt(root/"ActivityEnhancementFacets.dat")
            expected = elements["Y"]@a + p.enhancement_scale*p.enhancement_basis@b
            np.testing.assert_allclose(fields[:,8], expected, rtol=1e-14)
            report_fields = np.loadtxt(root/"activity_facets.dat")
            np.testing.assert_allclose(report_fields[:,8], expected, rtol=1e-14)
            np.testing.assert_allclose(report_fields[:,7],
                                       expected-elements["Y"]@coef[0], rtol=1e-14)
            self.assertTrue((root/"pressure_comparison_time.png").is_file())
            # The automatic publication hook uses audited in-memory rows:
            # no completed manifest exists yet at this point in the writer.
            self.assertTrue((root/"figures/pressures_time.eps").is_file())
            self.assertTrue((root/"figures/pressures_altitude.png").is_file())
            export = (root/"ActivityEnhancementFit.txt").read_text()
            self.assertIn("sha256:" + fit_cops.file_identity(root/"ActivityFit.txt")["sha256"], export)
            self.assertIn(fit_cops.JOINT_REGULARIZATION, export)
            # The independent C++ probe accepts one facet. Export a slice of
            # this very same fitted state, preserving the fitted baseline hash,
            # and verify C++'s single-clip result numerically.
            compiler = shutil.which("g++")
            if compiler:
                app = Path(fit_cops.__file__).resolve().parents[1]
                probe = root/"probe"
                subprocess.run([compiler, "-std=c++17", "-O2", "-I", str(app/"src"),
                                str(app/"tests/activity_facet_enhancement_probe.cpp"),
                                str(app/"src/ActivityEnhancement.cpp"),
                                str(app/"src/ActivityEnhancementBasis.cpp"),
                                "-o", str(probe)], check=True, capture_output=True)
                one = root/"one"; one.mkdir()
                p.enhancement_basis = p.enhancement_basis[:1]
                sliced = {k:v[:1] for k,v in elements.items()}
                path = fit_cops.write_enhancement_products(
                    str(one), p, b, sliced, np.ones(1, dtype=bool), 0.,
                    root/"ActivityFit.txt", fitted_activity=a)
                identity = "sha256:"+fit_cops.file_identity(root/"ActivityFit.txt")["sha256"]
                command = [str(probe), path, p.args.body_frame, identity]
                command += [repr(float(v)) for v in sliced["centroid"][0]]
                command += [repr(float(sliced["area"][0])), repr(float(sliced["Y"][0]@a))]
                replay = subprocess.run(command, check=True, capture_output=True, text=True)
                values = np.fromstring(replay.stdout, sep=" ")
                self.assertAlmostEqual(values[2], expected[0], places=13)
                self.assertAlmostEqual(values[3], max(0., expected[0]), places=13)

    def test_parameter_bounds_and_all_analytic_columns(self):
        p, _, _ = joint_problem()
        self.assertEqual(len(p.x0), 31)
        self.assertTrue(p.enhancement_mode)
        np.testing.assert_array_equal(p.lo[p.bslice], [0., 0.])
        x = p.x0.copy(); x[p.bslice] = [0.2, 0.3]
        p.lam = 0.7
        analytic = p.jacobian(x)
        for j in range(len(x)):
            h = 1e-5 * max(1., abs(x[j]))
            xp, xm = x.copy(), x.copy()
            xp[j] += h; xm[j] -= h
            numerical = (p.residuals(xp) - p.residuals(xm)) / (2*h)
            np.testing.assert_allclose(analytic[:, j], numerical,
                                       rtol=2e-6, atol=2e-9, err_msg=p.names[j])

    def test_zero_template_matches_free_legacy_problem(self):
        p, elements, coef = joint_problem()
        args = deepcopy(p.args)
        args.activity_template_names = []
        args.activity_template_basis = np.empty((len(elements["area"]), 0))
        legacy = fit_cops.Problem(p.M, elements, coef, args)
        np.testing.assert_allclose(p.data_residuals(p.x0),
                                   legacy.data_residuals(legacy.x0), atol=0, rtol=0)
        indices = list(range(p.na)) + list(range(p.na+p.nb, len(p.x0)))
        np.testing.assert_allclose(p.data_jacobian(p.x0)[:, indices],
                                   legacy.data_jacobian(legacy.x0), atol=0, rtol=0)

    def test_regularizer_is_total_map_with_cross_terms(self):
        p, elements, coef = joint_problem()
        x = p.x0.copy(); x[0] += 0.25; x[p.bslice] = [0.4, 0.3]
        a, b, _, _ = p.unpack(x)
        delta = elements["Y"] @ (a-coef[0]) + p.enhancement_scale*p.enhancement_basis@b
        expected = np.sum(elements["area"]*delta**2)/p.N0
        self.assertAlmostEqual(p.deviation(x)**2, expected, places=12)
        # A constant template duplicates the constant baseline direction;
        # this must remain a valid semidefinite map norm, without a ridge.
        args = deepcopy(p.args)
        args.activity_template_basis[:, 0] = elements["Y"][:, 0]
        dependent = fit_cops.Problem(p.M, elements, coef, args)
        y = dependent.x0.copy()
        y[0] -= dependent.enhancement_scale/dependent.ascale
        y[dependent.na] = 1.
        self.assertLess(dependent.deviation(y), 1e-7)


if __name__ == "__main__":
    unittest.main()

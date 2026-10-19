"""Stage-A6 synthetic recovery and cross-language validation tests.

These tests deliberately use the production ``Model.evaluate`` and ``Problem``
classes, but replace mission kernels with small deterministic transport
operators.  That separation has two purposes:

* recovery tests can inject an exactly known high-degree activity field; and
* visibility can be controlled explicitly, including a genuinely unsupported
  coefficient direction that must *not* be reported as recovered.

The tests are not a physical model of the Rosetta descent.  They validate the
algebra, parameter ordering, derivatives, regularization, and observability
logic used by the fitter before those components are applied to flight data.
"""

from argparse import Namespace
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.optimize import least_squares


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

import fit_cops
from activity_basis import coefficient_index, evaluate_cartesian, modes


def _unit_vectors(rng, count):
    """Return deterministic, nonzero directions distributed over the sphere."""
    vectors = rng.normal(size=(count, 3))
    return vectors / np.linalg.norm(vectors, axis=1)[:, None]


def synthetic_recovery_problem(max_degree=5, injected=None, noise_sigma=0.0,
                               geometry="full"):
    """Construct an in-memory enhancement-only fit with known truth.

    ``geometry='full'`` gives every observation an independent positive
    weighting of every facet.  It is intentionally overdetermined and is used
    to test recoverable modes.  ``geometry='north_cap'`` exposes only repeated
    north-pole facet directions.  At the pole every non-zonal (m != 0) real
    harmonic is exactly zero, creating a documented null space without using
    a numerical near-zero as a surrogate for invisibility.

    Temperature, legacy degree-0--4 activity, and CO2 kappa are fixed.  Thus
    the only free parameters are the normalized degree-5--``max_degree`` H2O
    residual coefficients introduced in Stage A4.
    """
    rng = np.random.default_rng(6102026 + max_degree)
    nface = 180
    npoint = 72 if max_degree == 5 else 120
    centroid = _unit_vectors(rng, nface)

    if geometry == "north_cap":
        # Duplicate directions are intentional: all contributing facets see
        # precisely the same harmonic row.  Remaining facets retain general
        # directions so the whole-surface regularization Gram matrix stays
        # positive definite.
        nseen = 12
        centroid[:nseen] = np.array([0.0, 0.0, 1.0])
        contributing_faces = np.arange(nseen)
    elif geometry == "full":
        contributing_faces = np.arange(nface)
    else:
        raise ValueError("unknown synthetic geometry: " + geometry)

    # The legacy field is a strictly positive constant.  This keeps every
    # recovery experiment away from the clipping kink; clipping itself and its
    # one-sided derivative convention are tested separately in A3/A4.
    legacy_basis = np.zeros((nface, 25), dtype=float)
    legacy_basis[:, 0] = 1.0
    area = rng.uniform(0.4, 1.8, nface)
    coefficients = np.zeros((3, 25), dtype=float)
    coefficients[0, 0] = 2.0
    coefficients[1, 0] = 0.1

    model = fit_cops.Model.__new__(fit_cops.Model)
    model.nfaces, model.npt = nface, npoint
    model.Y = legacy_basis
    model.face = np.tile(contributing_faces, npoint)
    model.pt = np.repeat(np.arange(npoint), len(contributing_faces))
    model.indptr = np.arange(npoint + 1, dtype=np.int64) * len(contributing_faces)
    ncontribution = len(model.face)
    model.lit = np.ones(ncontribution, dtype=bool)
    model.xlit = np.zeros(ncontribution, dtype=float)
    model.bm = {"H2O": 0.9, "CO2": 1.1}
    model.cNG = {"H2O": 0.8, "CO2": 1.0}
    model.cRG = {species: np.ones(npoint) * value
                 for species, value in (("H2O", 0.75), ("CO2", 0.95))}

    # Independent positive transport weights make the full geometry resolve
    # all requested modes.  CO2 is exactly zero because kappa is fixed and the
    # synthetic experiment is specifically about the H2O activity residual.
    model.wn = {
        "H2O": rng.lognormal(mean=0.0, sigma=0.7, size=ncontribution),
        "CO2": np.zeros(ncontribution),
    }
    model.wv = {
        "H2O": rng.lognormal(mean=-0.1, sigma=0.6, size=ncontribution),
        "CO2": np.zeros(ncontribution),
    }
    model.wp = {"H2O": rng.lognormal(mean=-0.2, sigma=0.5,
                                      size=ncontribution)}
    model.PhiCO2 = np.zeros(npoint)
    model.T0 = np.array([190.0, 145.0, 110.0])
    model.rows = np.arange(npoint)
    model.useNG = np.ones(npoint, dtype=int)
    model.useRG = np.ones(npoint, dtype=int)

    elements = {"centroid": centroid, "area": area, "Y": legacy_basis}
    basis = evaluate_cartesian(centroid, 5, max_degree)
    scale = 2.0                         # RMS of the constant legacy field
    truth = np.zeros(basis.shape[1], dtype=float)
    if injected:
        for mode, value in injected.items():
            truth[coefficient_index(mode[0], mode[1], 5, max_degree)] = value

    p_ng, p_rg = model.evaluate(
        coefficients[0], model.T0, 1.0, b=truth,
        enhancement_basis=basis, enhancement_scale=scale)
    if noise_sigma:
        # Pressure errors are multiplicative because the objective is in log
        # pressure.  Separate deterministic sequences prevent the two gauges
        # from receiving artificially identical noise.
        p_ng = p_ng * np.exp(rng.normal(scale=noise_sigma, size=npoint))
        p_rg = p_rg * np.exp(rng.normal(scale=noise_sigma, size=npoint))
    model.refNG, model.refRG = p_ng, p_rg

    args = Namespace(
        enhancement_degree=max_degree, fix_legacy_activity=True,
        enhancement_scale=scale, degree=4, ordered_temperature=False,
        fix_temperature=True, fix_kappa=True, tbounds=(30.0, 400.0),
        body_frame="67P/C-G_CK", max_nfev=1000, ftol=1e-12, xtol=1e-12,
        gtol=1e-12, jobs=1)
    problem = fit_cops.Problem(model, elements, coefficients, args)
    return problem, elements, coefficients, truth, contributing_faces


def solve(problem, lam=0.0):
    """Run the same bounded TRF least-squares method used by ``fit_cops``."""
    problem.lam = float(lam)
    return least_squares(
        problem.residuals, problem.x0, jac=problem.jacobian,
        bounds=(problem.lo, problem.hi), x_scale="jac", method="trf",
        max_nfev=problem.args.max_nfev, ftol=problem.args.ftol,
        xtol=problem.args.xtol, gtol=problem.args.gtol)


class SyntheticRecoveryTests(unittest.TestCase):
    """A6 requirements 5--9: derivatives, recovery, visibility, and lambda."""

    def test_finite_difference_jacobian_at_nonzero_multimode_state(self):
        injected = {(5, -3): 0.06, (5, 2): -0.04, (6, 0): 0.03,
                    (6, 5): -0.02}
        problem, _, _, truth, _ = synthetic_recovery_problem(6, injected)
        analytic = problem.jacobian(truth)
        worst_relative = 0.0
        for column in range(problem.nb):
            step = 2.0e-7
            plus, minus = truth.copy(), truth.copy()
            plus[column] += step
            minus[column] -= step
            finite = ((problem.residuals(plus) - problem.residuals(minus)) /
                      (2.0 * step))
            scale = max(np.max(np.abs(finite)), 1.0e-300)
            worst_relative = max(
                worst_relative,
                float(np.max(np.abs(analytic[:, column] - finite)) / scale))
            np.testing.assert_allclose(analytic[:, column], finite,
                                       rtol=4e-6, atol=3e-9)
        print("A6 synthetic Jacobian: %d columns, worst relative difference %.3e" %
              (problem.nb, worst_relative))

    def test_exact_single_degree_five_mode_is_recovered(self):
        problem, _, _, truth, _ = synthetic_recovery_problem(
            5, {(5, -2): 0.12})
        result = solve(problem, lam=0.0)
        rank = np.linalg.matrix_rank(problem.data_jacobian(result.x))
        self.assertTrue(result.success, result.message)
        self.assertEqual(rank, problem.nb)
        np.testing.assert_allclose(result.x, truth, rtol=2e-8, atol=2e-10)
        residual_norm = np.linalg.norm(problem.data_residuals(result.x))
        coefficient_error = np.max(np.abs(result.x - truth))
        self.assertLess(residual_norm, 1e-10)
        print("A6 exact degree-5 recovery: rank %d/%d, max |db| %.3e, residual norm %.3e" %
              (rank, problem.nb, coefficient_error, residual_norm))

    def test_several_modes_with_controlled_noise_recover_supported_field(self):
        injected = {(5, -4): 0.08, (5, 1): -0.055, (6, -2): 0.045,
                    (6, 0): -0.035, (6, 6): 0.025}
        problem, _, _, truth, _ = synthetic_recovery_problem(
            6, injected, noise_sigma=2.0e-4)
        result = solve(problem, lam=0.0)
        self.assertTrue(result.success, result.message)
        jacobian = problem.data_jacobian(result.x)
        self.assertEqual(np.linalg.matrix_rank(jacobian), problem.nb)

        true_field = problem.enhancement_basis @ truth
        recovered_field = problem.enhancement_basis @ result.x
        relative_field_error = (np.linalg.norm(recovered_field - true_field) /
                                np.linalg.norm(true_field))
        correlation = np.corrcoef(true_field, recovered_field)[0, 1]
        self.assertLess(relative_field_error, 0.08)
        self.assertGreater(correlation, 0.997)
        print("A6 noisy multimode recovery: field relative error %.3e, correlation %.9f" %
              (relative_field_error, correlation))

    def test_unseen_non_zonal_signal_is_explicitly_unsupported(self):
        # At the north pole Y_5,+1 is exactly zero.  The injected field is
        # therefore nonzero over the unseen surface but produces no pressure
        # change in this controlled observation geometry.
        problem, elements, _, truth, seen_faces = synthetic_recovery_problem(
            5, {(5, 1): 0.2}, geometry="north_cap")
        seen = np.zeros(len(elements["area"]), dtype=bool)
        seen[seen_faces] = True
        field = problem.enhancement_basis @ truth
        self.assertLess(np.max(np.abs(field[seen])), 1e-14)
        self.assertGreater(np.sqrt(np.mean(field[~seen] ** 2)), 1e-2)

        jacobian = problem.data_jacobian(problem.x0)
        rank = np.linalg.matrix_rank(jacobian)
        self.assertLess(rank, problem.nb)
        result = solve(problem, lam=0.0)
        self.assertTrue(result.success, result.message)
        self.assertLess(np.linalg.norm(result.x), 1e-12)
        self.assertGreater(np.linalg.norm(result.x - truth), 0.19)
        print("A6 unsupported geometry: rank %d/%d, seen max %.3e, unseen RMS %.3e, recovered norm %.3e" %
              (rank, problem.nb, np.max(np.abs(field[seen])),
               np.sqrt(np.mean(field[~seen] ** 2)), np.linalg.norm(result.x)))

    def test_regularization_limits_and_zero_lambda_rank_diagnostic(self):
        problem, _, _, truth, _ = synthetic_recovery_problem(
            5, {(5, -2): 0.12, (5, 3): -0.07})
        unregularized = solve(problem, lam=0.0)
        strong = solve(problem, lam=1.0e10)
        self.assertTrue(unregularized.success, unregularized.message)
        self.assertTrue(strong.success, strong.message)
        self.assertEqual(np.linalg.matrix_rank(
            problem.data_jacobian(unregularized.x)), problem.nb)
        np.testing.assert_allclose(unregularized.x, truth,
                                   rtol=2e-8, atol=2e-10)
        self.assertLess(np.linalg.norm(strong.x),
                        1e-4 * np.linalg.norm(unregularized.x))
        self.assertGreater(np.linalg.norm(problem.data_residuals(strong.x)),
                           np.linalg.norm(problem.data_residuals(unregularized.x)))
        print("A6 lambda limits: |b_strong|/|b_lambda0| %.3e" %
              (np.linalg.norm(strong.x) / np.linalg.norm(unregularized.x)))


class CoefficientRoundTripTests(unittest.TestCase):
    """A6 requirement 3: Python writer -> strict C++ reader/evaluator."""

    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("g++ is unavailable")
        cls.temporary = tempfile.TemporaryDirectory(prefix="cops_a6_roundtrip_")
        cls.root = Path(cls.temporary.name)
        cls.probe = cls.root / "activity_enhancement_probe"
        command = [
            compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
            "-I", str(APP / "src"),
            str(APP / "tests/activity_enhancement_probe.cpp"),
            str(APP / "src/ActivityEnhancement.cpp"),
            str(APP / "src/ActivityEnhancementBasis.cpp"),
            "-o", str(cls.probe),
        ]
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError("C++ A6 round-trip probe failed to compile:\n" +
                               completed.stdout)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_writer_order_and_binary64_values_round_trip_through_cpp(self):
        problem, elements, _, _, _ = synthetic_recovery_problem(8)
        values = np.linspace(-0.125, 0.175, problem.nb)
        legacy = self.root / "ActivityFit.txt"
        legacy.write_text("deterministic legacy identity\n", encoding="utf-8")
        fit_cops.write_enhancement_products(
            self.root, problem, values, elements,
            np.ones(len(elements["area"]), dtype=bool), 3.0e-4, legacy)
        exchange = self.root / "ActivityEnhancementFit.txt"

        records = [line.split() for line in exchange.read_text().splitlines()
                   if line.startswith("Coefficient =")]
        written_modes = [(int(record[2]), int(record[3])) for record in records]
        written_values = np.array([float(record[4]) for record in records])
        self.assertEqual(written_modes, modes(5, 8))
        np.testing.assert_array_equal(written_values, values)

        identity = "sha256:" + hashlib.sha256(legacy.read_bytes()).hexdigest()
        point = np.array([0.31, -0.47, 0.82])
        completed = subprocess.run(
            [str(self.probe), str(exchange), "67P/C-G_CK", identity,
             *("%.17g" % value for value in point), "2.0"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        maximum, count, residual, total, applied = (
            float(value) for value in completed.stdout.split())
        expected = problem.enhancement_scale * (
            evaluate_cartesian(point, 5, 8) @ values)
        self.assertEqual((int(maximum), int(count)), (8, 56))
        self.assertAlmostEqual(residual, expected, places=13)
        self.assertAlmostEqual(total, 2.0 + expected, places=13)
        self.assertAlmostEqual(applied, max(0.0, 2.0 + expected), places=13)
        print("A6 coefficient round trip: 56 modes, residual absolute difference %.3e" %
              abs(residual - expected))


if __name__ == "__main__":
    unittest.main()

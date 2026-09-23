#!/usr/bin/env python3
"""Deterministic unit tests for the F4 output validator.

No AMPS executable is required.  The fixture writes the same three Tecplot-style
files as the solver, reconstructs every nominal/lower/upper product independently,
and then injects targeted corruptions.  The important regression case contains a
node with N_resolved=0: its nominal values are NaN while both uncertainty bounds
remain finite and fully gated.
"""

from __future__ import print_function

import importlib.util
import io
import json
import math
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path


F4_DIR = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("run_F4", str(F4_DIR / "run_F4.py"))
F4 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F4)


def _number(value):
    if math.isnan(value):
        return "nan"
    return "%.17g" % value


def _write_zone(path, title, variables, rows):
    with path.open("w") as stream:
        stream.write('TITLE="%s"\n' % title)
        stream.write("VARIABLES=" + " ".join('"%s"' % name for name in variables) + "\n")
        stream.write('ZONE T="fixture" I=%d F=POINT\n' % len(rows))
        for row in rows:
            stream.write(" ".join(_number(value) for value in row) + "\n")


def _make_fixture(workdir, all_resolved=False, corruption=None):
    """Write one-point F4 output with exact independently computed products."""
    E = [1.0, 10.0, 100.0, 1000.0]
    sampled = [16.0] * len(E)
    if all_resolved:
        resolved = [16.0] * len(E)
        nominal = [0.25, 0.50, 0.75, 1.00]
    else:
        resolved = [16.0, 0.0, 8.0, 16.0]
        nominal = [0.25, float("nan"), 0.75, 1.00]

    lower = []
    upper = []
    unresolved = []
    for t, n_sampled, n_resolved in zip(nominal, sampled, resolved):
        fraction = (n_sampled - n_resolved) / n_sampled
        unresolved.append(fraction)
        t_lower = t * n_resolved / n_sampled if n_resolved > 0.0 else 0.0
        lower.append(t_lower)
        upper.append(t_lower + fraction)

    Jb = [F4.j_boundary(energy) for energy in E]
    Jl = [t * jb for t, jb in zip(nominal, Jb)]
    Jl_lower = [t * jb for t, jb in zip(lower, Jb)]
    Jl_upper = [t * jb for t, jb in zip(upper, Jb)]

    # Compute all file-level products before injecting a spectrum corruption.
    # That guarantees each negative test exercises a real closure gate rather
    # than regenerating the expected output from the same bad value.
    density = F4.density_from_local_spectrum(E, Jl, F4.EMIN, F4.EMAX)
    density_lower = F4.density_from_local_spectrum(E, Jl_lower, F4.EMIN, F4.EMAX)
    density_upper = F4.density_from_local_spectrum(E, Jl_upper, F4.EMIN, F4.EMAX)
    flux_total = F4.flux_from_local_spectrum(E, Jl, F4.EMIN, F4.EMAX)
    flux_total_lower = F4.flux_from_local_spectrum(E, Jl_lower, F4.EMIN, F4.EMAX)
    flux_total_upper = F4.flux_from_local_spectrum(E, Jl_upper, F4.EMIN, F4.EMAX)

    channel_values = []
    for _name, e1, e2 in F4.ENERGY_BINS:
        channel_values.extend([
            F4.flux_channel_as_amps(E, nominal, e1, e2),
            F4.flux_channel_as_amps(E, lower, e1, e2),
            F4.flux_channel_as_amps(E, upper, e1, e2),
        ])

    if corruption == "finite_empty_nominal":
        nominal[1] = 0.0
        Jl[1] = 0.0
    elif corruption == "lower_spectrum":
        Jl_lower[0] *= 1.01
    elif corruption == "reversed_interval":
        lower[2] = upper[2] + 0.1

    spectrum_variables = [
        "E_MeV", "T", "T_lower", "T_upper", "unresolved_fraction",
        "N_sampled", "N_resolved", "J_boundary_perMeV",
        "J_local_perMeV", "J_local_lower_perMeV", "J_local_upper_perMeV",
    ]
    spectrum_rows = list(zip(
        E, nominal, lower, upper, unresolved, sampled, resolved,
        Jb, Jl, Jl_lower, Jl_upper,
    ))
    _write_zone(workdir / "gridless_points_spectrum.dat", "F4 fixture spectrum",
                spectrum_variables, spectrum_rows)

    density_variables = [
        "X_km", "Y_km", "Z_km", "N_m^-3", "N_lower_m^-3",
        "N_upper_m^-3", "N_cm^-3", "N_lower_cm^-3", "N_upper_cm^-3",
        "Rc_lower_GV", "Rc_effective_GV", "Rc_upper_GV", "PenumbraWidth_GV",
        "T_high",
    ]
    density_row = [
        F4.RE_KM + F4.ALT_KM_DEFAULT, 0.0, 0.0,
        density, density_lower, density_upper,
        density * 1.0e-6, density_lower * 1.0e-6, density_upper * 1.0e-6,
        0.0, 0.0, 0.0, 0.0, 0.0,
    ]
    _write_zone(workdir / "gridless_points_density.dat", "F4 fixture density",
                density_variables, [density_row])

    flux_variables = [
        "X_km", "Y_km", "Z_km", "F_tot_m2s1",
        "F_tot_lower_m2s1", "F_tot_upper_m2s1",
    ]
    for name, _e1, _e2 in F4.ENERGY_BINS:
        flux_variables.extend([
            "F_%s_m2s1" % name,
            "F_%s_lower_m2s1" % name,
            "F_%s_upper_m2s1" % name,
        ])
    flux_row = [
        F4.RE_KM + F4.ALT_KM_DEFAULT, 0.0, 0.0,
        flux_total, flux_total_lower, flux_total_upper,
    ] + channel_values
    _write_zone(workdir / "gridless_points_flux.dat", "F4 fixture flux",
                flux_variables, [flux_row])


class F4ValidatorTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="f4_validator_"))
        self.points = F4.build_points(F4.ALT_KM_DEFAULT, [0.0], [0.0])
        self.args = F4.parse_args([
            "--skip-run", "--workdir", str(self.root),
            "--scan-n", "4", "--max-particles", "64",
            "--lons", "0", "--lats", "0",
        ])
        # main() normally materializes these parsed list values before analyze().
        self.args.lons_values = [0.0]
        self.args.lats_values = [0.0]

    def tearDown(self):
        shutil.rmtree(str(self.root))

    def _failed_check_names(self):
        with (self.root / "F4_result.json").open() as stream:
            result = json.load(stream)
        return set(row["check"] for row in result["failed_checks"])

    def _analyze(self):
        # Negative fixtures intentionally print a detailed failure report in
        # production.  Suppress it here so the unit-test output stays focused on
        # the assertion names while still exercising the identical code path.
        with redirect_stdout(io.StringIO()):
            return F4.analyze(self.root, self.points, self.args)

    def test_numeric_gate_defaults_are_unchanged(self):
        self.assertEqual(self.args.closure_tol, 2.0e-5)
        self.assertEqual(self.args.integral_tol, 2.0e-5)
        self.assertEqual(self.args.differential_tol, 2.0e-5)
        self.assertEqual(self.args.t_bounds_tol, 1.0e-12)

    def test_all_resolved_nominal_products_pass(self):
        _make_fixture(self.root, all_resolved=True)
        self.assertTrue(self._analyze())

    def test_empty_resolved_node_uses_finite_gated_bounds(self):
        _make_fixture(self.root, all_resolved=False)
        self.assertTrue(self._analyze())

    def test_skip_run_entrypoint_accepts_the_same_bounded_fixture(self):
        _make_fixture(self.root, all_resolved=False)
        argv = [
            "--skip-run", "--workdir", str(self.root),
            "--scan-n", "4", "--max-particles", "64",
            "--lons", "0", "--lats", "0",
        ]
        with redirect_stdout(io.StringIO()):
            self.assertEqual(F4.main(argv), 0)

    def test_empty_resolved_node_cannot_be_silently_mapped_to_zero(self):
        _make_fixture(self.root, corruption="finite_empty_nominal")
        self.assertFalse(self._analyze())
        self.assertIn("nominal_definedness_contract", self._failed_check_names())

    def test_corrupted_lower_spectrum_fails_existing_closure_tolerance(self):
        _make_fixture(self.root, corruption="lower_spectrum")
        self.assertFalse(self._analyze())
        self.assertIn("local_spectrum_lower_closure", self._failed_check_names())

    def test_reversed_uncertainty_interval_fails(self):
        _make_fixture(self.root, corruption="reversed_interval")
        self.assertFalse(self._analyze())
        self.assertIn("transmission_bounds", self._failed_check_names())


if __name__ == "__main__":
    unittest.main(verbosity=2)

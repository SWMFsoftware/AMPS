#!/usr/bin/env python3
"""Independent mathematical/physical limits for the validation-only warm solver."""
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest import mock

import mpmath as mp

MODULE = Path(__file__).resolve().parents[1] / "validation/physics/warm_parallel_maxwellian.py"
SPEC = importlib.util.spec_from_file_location("srcsep_warm_parallel", str(MODULE))
WARM = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = WARM
SPEC.loader.exec_module(WARM)


def controlled_case():
    # These are explicitly declared verification inputs, not inferred production
    # composition/thermodynamics or a new production initialization default.
    return {"provenance": "controlled uniform ep verification fixture; no observational or production qualification",
            "distribution": "isotropic-nonrelativistic-Maxwellian", "propagation": "parallel",
            "proton_density_m3": "5000000", "electron_density_m3": "5000000",
            "magnetic_field_t": "5e-9", "proton_temperature_k": "120000",
            "electron_temperature_k": "120000", "proton_mass_kg": "1.67262192369e-27",
            "electron_mass_kg": "9.1093837015e-31", "charge_magnitude_c": "1.602176634e-19",
            "light_speed_mps": "299792458", "vacuum_permeability_hpm": "1.256637062e-6",
            "boltzmann_jpk": "1.380649e-23", "proton_drift_mps": "0", "electron_drift_mps": "0"}


class WarmParallelTests(unittest.TestCase):
    def test_landau_contour_independent_of_complex_erfc(self):
        with mp.workdps(50):
            for point in (mp.mpc(0), mp.mpc(1, ".5"), mp.mpc(2, "-.3"), mp.mpc(0, "-.5"),
                          mp.mpc(100, ".1"), mp.mpc(100000, "-.1"), mp.mpc(1, "-3")):
                offset = max(mp.mpf(0), -mp.im(point))+mp.mpf(".5")
                actual = WARM.plasma_z(point)
                for distance in (offset, offset+mp.mpf(".25")):
                    independent = WARM.z_contour(point, distance)
                    self.assertLess(abs(actual-independent), mp.mpf("1e-40")*max(1, abs(actual)))
            self.assertLess(abs(WARM.plasma_z(0)-mp.j*mp.sqrt(mp.pi)), mp.mpf("1e-45"))

    def test_cold_plasma_limit_and_signed_resonance(self):
        with mp.workdps(65):
            case = controlled_case()
            case["proton_temperature_k"] = "0.000001"
            case["electron_temperature_k"] = "0.000001"
            plasma = WARM.Plasma(case)
            scales = plasma.scales()
            mass, a2 = scales["mass_ratio"], scales["alfven_ion_over_c"]**2
            z, k = mp.mpc(".3", ".1"), mp.mpf(".4")
            operator = WARM.Dispersion(plasma)
            for h in (-1, 1):
                cold = z*z*a2-k*k-z/(z+h)-mass*z/(z-h*mass)
                actual = operator.evaluate(z, k, h)[0]
                self.assertLess(abs(actual-cold), mp.mpf("1e-11"))
            # Isotropic no-drift counterpropagating symmetry, including the
            # signed ion/electron gyrofrequencies and lower-half-plane residue.
            operator = WARM.Dispersion(WARM.Plasma(controlled_case()))
            z = mp.mpc(".3", "-.1")
            left = operator.evaluate(-mp.conj(z), k, 1)[0]
            right = mp.conj(operator.evaluate(z, k, -1)[0])
            self.assertLess(abs(left-right), mp.mpf("1e-55"))

    def test_first_warm_circular_mhd_correction(self):
        with mp.workdps(60):
            plasma = WARM.Plasma(controlled_case())
            s = plasma.scales()
            mass, beta, beta_e = s["mass_ratio"], s["beta_proton"], s["beta_electron"]
            aa = 1+1/mass+s["alfven_ion_over_c"]**2
            leading = 1/mp.sqrt(aa)
            coefficient = ((1-1/mass**2)/aa+(beta-beta_e/mass)/2)/(2*aa)
            # Independent expansion of the full kinetic susceptibility:
            # z=K/sqrt(A)+h*K^2*((1-M^-2)/A+(beta_p-beta_e/M)/2)/(2*A)+O(K^3).
            for h in (-1, 1):
                k = mp.mpf(".00001")
                root, _ = WARM.newton(WARM.Dispersion(plasma), k, h, leading*k, 40)
                measured = (mp.re(root)/k-leading)/k
                self.assertLess(abs(measured-h*coefficient), mp.mpf("0.00003"))
                energy = WARM.Dispersion(plasma).weak_damping_energy_ratio(mp.re(root), k, h)
                self.assertLess(abs(energy-2), 10*k)

    def test_precisions_and_numerical_root_isolation(self):
        plasma = WARM.Plasma(controlled_case())
        for h in (-1, 1):
            result = WARM.solve_point(plasma, ".1", h, 30, 50)
            self.assertEqual(result["status"], "RESOLVED")
            self.assertEqual(result["numerical_isolation"]["winding_count"], 1)
            self.assertLess(float(result["precision_refinement_relative"]), 1e-15)
            self.assertLess(float(result["continuation_refinement_relative"]), 1e-15)
            self.assertFalse(result["damping_resolved"])
            self.assertTrue(float(result["phase_speed_over_model_va_ion"]) > 1 if h == 1 else
                            float(result["phase_speed_over_model_va_ion"]) < 1)

    def test_small_k_neutral_cancellation_increases_precision(self):
        result = WARM.solve_point(WARM.Plasma(controlled_case()), "1e-30", 1, 25, 45)
        self.assertEqual(result["status"], "RESOLVED")
        self.assertGreater(result["additional_cancellation_guard_digits"], 50)
        self.assertEqual(result["numerical_isolation"]["winding_count"], 1)
        self.assertLess(float(result["continuation_refinement_relative"]), 1e-10)

    def test_default_precision_serialized_frequency_roundtrip(self):
        plasma = WARM.Plasma(controlled_case())
        point = ".1000000000000000000000000000123456789012345"
        result = WARM.solve_point(plasma, point, 1, 40, 70)
        self.assertEqual(result["requested_k_di"], point)
        with mp.workdps(120):
            omega = mp.mpc(result["omega_real_radps"], result["omega_imaginary_raw_radps"])
            gyro = plasma.scales()["ion_gyrofrequency_radps"]
            value, derivative, _ = WARM.Dispersion(plasma).evaluate(omega/gyro, mp.mpf(result["k_di"]), 1)
            correction = abs(value/derivative)*gyro
            self.assertLessEqual(correction, mp.mpf(result["frequency_numerical_error_radps"]))
            self.assertLessEqual(mp.mpf(result["frequency_serialization_roundoff_radps"]),
                                 mp.mpf(result["frequency_numerical_error_radps"]))

    def test_winding_alias_negative_control(self):
        class ManyRoots:
            def evaluate(self, value, _k, _h):
                w = value-1
                a = mp.mpf(1)/16
                return w*(w**64-a**64), 65*w**64-a**64, mp.mpf(0)
        with mp.workdps(70):
            # This polynomial has65 zeros inside r=1/8. A64-point phase
            # unwrap aliases them to1; the independent F'/F integral must not.
            operator = ManyRoots()
            measured = WARM.winding(operator, mp.mpf(1), mp.mpf(1), 1, mp.mpf(1)/8, 64)
            self.assertEqual(measured[0], 1)
            self.assertLess(abs(measured[3]-65), mp.mpf("1e-15"))
            isolated = WARM.isolate(operator, mp.mpf(1), mp.mpf(1), 1)
            self.assertLess(mp.mpf(isolated["contour_radius"]), mp.mpf(1)/16)

    def test_singular_predictor_is_explicit_unresolved(self):
        with mp.workdps(50), \
                mock.patch.object(WARM, "newton", return_value=(mp.mpc(".0001"), 1)), \
                mock.patch.object(WARM.Dispersion, "evaluate", return_value=(mp.mpc(1), mp.mpc(0), mp.mpc(1))):
            with self.assertRaises(WARM.Unresolved):
                WARM.trace_branch(WARM.Plasma(controlled_case()), mp.mpf(".1"), 1, 30, mp.mpf(".1"))

    def test_no_implicit_species_temperature_or_distribution(self):
        for field in ("electron_temperature_k", "electron_density_m3", "electron_mass_kg", "proton_drift_mps"):
            case = controlled_case()
            del case[field]
            with self.assertRaises((ValueError, KeyError)):
                WARM.Plasma(case).numbers()
        case = controlled_case()
        case["electron_density_m3"] = "4000000"
        with self.assertRaises(ValueError):
            WARM.Plasma(case).numbers()
        case = controlled_case()
        case["distribution"] = "unspecified"
        with self.assertRaises(ValueError):
            WARM.Plasma(case).numbers()


if __name__ == "__main__":
    unittest.main(verbosity=2)

"""Independent source-export contracts, using synthetic rates and geometry.

These fixtures are mathematical test inputs, not replacements for published
catalog coordinates. Real campaigns resolve the external catalog and record
its hash. The deliberately unequal rates expose accidental equal partition,
normalization of absolute rates, dropped zeros or nuisance-background export.
"""
import copy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fit_source_export import export_vectors, export_single
from fitting import load_config, solve_amplitudes
from source_contract import CASES, HYBRID, SOLAR, source_ids
import numpy as np
from joint_fitting_study import solve_combined_blocks, strict_diagnostic_json


def fixture_config(case):
    if case == HYBRID:
        return load_config(ROOT / "config/tenishev2014_stage0.json",
            stripe_surface_mapping="local_vertical_to_reference_sphere")
    return {"schema_version": 1, "source_case": case, "frame": "IAU_ENCELADUS",
        "reference_radius_m": 252100.0, "source_law": "tenishev2010_point",
        "discrete_jets": [{"id": label[4:], "latitude_deg": -80.0,
            "west_longitude_deg": float(i), "tilt_deg": 0.0 if case == SOLAR else 30.0,
            "azimuth_deg": 0.0, "temperature_K": 180.0, "bulk_speed_m_s": 1500.0}
            for i, label in enumerate(source_ids(case))], "tiger_stripes": []}


def fixture_rates(case):
    return {name: (0.0 if i == 0 else (i+1)*1.0e22)
            for i, name in enumerate(source_ids(case))}


class AutomaticExportTests(unittest.TestCase):
    def test_stagnated_trf_recovers_same_nnls_problem_with_kkt_gate(self):
        from scipy.optimize import lsq_linear as original
        def stagnated(*args,**kwargs):
            result=original(*args,**kwargs)
            result.success=False; result.status=-1; result.message="forced stagnation"
            return result
        matrix=np.column_stack(([1.,2.,4.],np.ones(3),[1.,2.,4.]))
        truth=matrix@np.array([2.,3.,0.])
        with patch("fitting.lsq_linear",side_effect=stagnated):
            fitted=solve_amplitudes(matrix,truth,np.ones(3))
        np.testing.assert_allclose(fitted["predicted"],truth,atol=1e-12)
        diagnostic=fitted["solver_diagnostics"]
        self.assertEqual(diagnostic["method"],"scipy_nnls_recovery_after_trf_failure")
        self.assertLessEqual(diagnostic["kkt_negative_gradient_over_target_norm"],1e-12)

    def test_full_turn_catalog_azimuth_is_same_axis_in_half_open_receipt(self):
        config=fixture_config(CASES[1])
        config["discrete_jets"][84]["azimuth_deg"]=360.0
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"export"
            export_vectors(config,[("test","E3","H2O",fixture_rates(CASES[1]),"H2O")],out,{})
            snapshot=json.loads((out/"resolved_geometry.json").read_text())
            self.assertEqual(snapshot["discrete_jets"][84]["azimuth_deg"],0.0)
            self.assertEqual(snapshot["_angle_endpoint_canonicalization"][0]["fit_value"],360.0)
            self.assertEqual(config["discrete_jets"][84]["azimuth_deg"],360.0)

    def test_singular_diagnostic_is_valid_json_but_invalid_physics_is_rejected(self):
        encoded = json.dumps(strict_diagnostic_json({
            "condition_number_scaled": float("inf"), "rate": 1e28}), allow_nan=False)
        self.assertEqual(json.loads(encoded)["condition_number_scaled"], "Infinity")
        for key in ("rate", "chi2", "condition_number_scaled"):
            with self.assertRaises(ValueError):
                strict_diagnostic_json({key: float("nan")})
        with self.assertRaises(ValueError):
            strict_diagnostic_json({"rate": float("inf")})

    def test_separate_observation_backgrounds_recover_independent_truth(self):
        # One common physical coefficient and two genuinely different offsets
        # give an overdetermined, identifiable block fixture. The shared-offset
        # alternative cannot reproduce both profiles, so the control detects
        # accidentally reconnecting the two nuisance columns.
        matrix = np.column_stack(([1.0,2.0,4.0], np.ones(3)))
        blocks = [("E3", matrix, matrix@np.array([2.0,3.0]), np.ones(3)),
                  ("E5", matrix, matrix@np.array([2.0,11.0]), np.ones(3))]
        fit = solve_combined_blocks(blocks, labels=["jet:a", "background:c0"],
                                    separate_backgrounds=True)
        np.testing.assert_allclose(fit["dataset_coefficients"]["E3"], [2.0,3.0], atol=1e-12)
        np.testing.assert_allclose(fit["dataset_coefficients"]["E5"], [2.0,11.0], atol=1e-12)
        np.testing.assert_allclose(fit["coefficients"], [2.0,0.0], atol=1e-12)
        self.assertLess(fit["chi2"], 1e-20)
        self.assertGreater(solve_combined_blocks(blocks)["chi2"], 1.0)

    def test_all_four_cases_preserve_absolute_and_shape_rates(self):
        for case in CASES:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as tmp:
                destination = Path(tmp) / "export"
                rates = fixture_rates(case)
                index = export_vectors(fixture_config(case), [(
                    "independent_each_observation", "E3", "H2O", rates, "H2O")],
                    destination, {"fixture": "synthetic independent nonuniform rates"})
                self.assertEqual(index["source_case"], case)
                for mode in ("absolute", "shape"):
                    files = index["artifacts"][0]["files"][mode]
                    m = json.loads((destination / files["manifest"]).read_text())
                    expected = {key: value/math.fsum(rates.values()) if mode == "shape" else value
                                for key, value in rates.items()}
                    self.assertEqual(m["species"]["H2O"], expected)
                    text = (destination / files["source_file"]).read_text()
                    self.assertEqual(text.count("#begin source: "), len(rates))
                    self.assertNotIn("background:", text)
                    self.assertIn("source_case="+case, text)

    def test_fitted_and_exported_stripe_vertices_are_identical(self):
        config = fixture_config(HYBRID)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "export"
            export_vectors(config, [("test", "E3", "H2O", fixture_rates(HYBRID), "H2O")], out, {})
            snapshot = json.loads((out / "resolved_geometry.json").read_text())
            self.assertEqual([s["points_m"] for s in snapshot["tiger_stripes"]],
                             [s["points_m"] for s in config["tiger_stripes"]])
            for stripe in snapshot["tiger_stripes"]:
                for point in stripe["points_m"]:
                    self.assertAlmostEqual(math.sqrt(math.fsum(x*x for x in point)),252100.0,places=7)

    def test_apparent_co2_absolute_export_is_explicitly_deferred(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "export"
            index = export_vectors(fixture_config(HYBRID), [(
                "test", "E5", "CO2", fixture_rates(HYBRID), "apparent_CO2")], out, {})
            item = index["artifacts"][0]
            self.assertEqual(set(item["files"]), {"shape"})
            self.assertEqual(item["absolute_export_status"], "DEFERRED_USER_FRAGMENTATION_CHOICE")
            self.assertFalse(list(out.glob("*absolute*")))

    def test_single_e3_co2_fit_has_same_absolute_deferral(self):
        config=fixture_config(HYBRID);config["molecular_species"]="CO2"
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"fit.json";path.write_text("{}")
            index=export_single({"config":config,"source_rate_s":fixture_rates(HYBRID)},path,"E3")
            self.assertEqual(set(index["artifacts"][0]["files"]),{"shape"})

    def test_raw_legacy_fit_is_not_silently_projected_after_fitting(self):
        config = load_config(ROOT / "config/tenishev2014_stage0.json")
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "fit/export stripe geometry differs"):
                export_vectors(config, [("test", "E3", "H2O", fixture_rates(HYBRID), "H2O")],
                               Path(tmp) / "export", {})

    def test_single_fit_excludes_background_and_keeps_fitted_kinematics(self):
        config = fixture_config(HYBRID)
        config["discrete_jets"][0]["bulk_speed_m_s"] = 812.5
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fit.json"
            result = {"config": config, "source_rate_s": fixture_rates(HYBRID),
                      "component_amplitudes": {"background:c0_m-3": 8e13}}
            path.write_text(json.dumps(result))
            export_single(result, path, "E3")
            geometry = json.loads((Path(tmp)/"fit_amps_sources/resolved_geometry.json").read_text())
            self.assertEqual(geometry["discrete_jets"][0]["bulk_speed_m_s"], 812.5)

    def test_zero_missing_and_overwrite_inputs_fail_closed(self):
        config = fixture_config(HYBRID)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/"export"
            zero = {key: 0.0 for key in source_ids(HYBRID)}
            index = export_vectors(config, [("test", "E3", "H2O", zero, "H2O")], out, {})
            self.assertEqual(index["artifacts"][0]["status"], "DEFERRED")
            self.assertFalse(list(out.glob("*.in")))
            with self.assertRaises(FileExistsError):
                export_vectors(config, [], out, {})
            bad = fixture_rates(HYBRID);bad.pop(next(iter(bad)))
            with self.assertRaisesRegex(ValueError, "incomplete"):
                export_vectors(config, [("test", "E3", "H2O", bad, "H2O")],Path(tmp)/"bad",{})


if __name__ == "__main__":
    unittest.main()

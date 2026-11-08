#!/usr/bin/env python3
"""Focused deterministic tests for the combined 3I/ATLAS runtime input."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[3]
MODULE_PATH = REPOSITORY / "srcHartley2" / "i3" / "prepare_i3_cases.py"
SPEC = importlib.util.spec_from_file_location("prepare_i3_cases", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class I3CaseGeneratorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.cases = MODULE._load_csv(
            REPOSITORY / "srcHartley2" / "i3" / "reference" / "i3_cases.csv"
        )
        self.reactions = MODULE._load_reactions(
            REPOSITORY
            / "srcHartley2"
            / "i3"
            / "reference"
            / "photochemistry_channels.csv"
        )

    def test_catalog_and_photochemistry_are_complete(self) -> None:
        self.assertEqual([row["case_id"] for row in self.cases], list("123456"))
        self.assertEqual(len(self.reactions), 17)
        self.assertEqual(
            {reaction.reactants for reaction in self.reactions},
            {"H2O", "OH", "H2", "H", "O"},
        )
        self.assertEqual(
            [float(row["exposure_mean_f107p_sfu"]) for row in self.cases],
            [
                136.194566908,
                142.874791643,
                139.711261811,
                143.568741662,
                150.444624050,
                138.071807780,
            ],
        )
        self.assertEqual(
            [float(row["oh00_photon_g_at_rh_s-1"]) for row in self.cases],
            [
                1.801508954861e-4,
                1.879160202915e-4,
                2.003402199802e-4,
                1.438101113966e-4,
                1.261056268403e-4,
                2.236355943965e-4,
            ],
        )

    def test_surface_resolution_uses_minimum_parent_mixture_mfp(self) -> None:
        expected = [
            1.7526150930392894,
            1.2978528764154647,
            0.5313123870760047,
            0.13144233686591544,
            0.12273476681866985,
            0.13144233686591544,
        ]
        for case, expected_value in zip(self.cases, expected):
            paths, limiting, selected = MODULE.derive_surface_resolution(case)
            self.assertEqual(limiting, "CO2")
            self.assertEqual(selected, min(paths.values()))
            self.assertAlmostEqual(selected, expected_value, places=12)

    def test_injection_and_recoil_controls(self) -> None:
        self.assertEqual(MODULE.derive_reference_injection_count(), 114)
        self.assertEqual(
            MODULE.derive_characteristic_speeds(self.reactions),
            {"H": 46945, "H2": 18783, "O": 4660, "OH": 1997, "H2O": 400, "CO": 400, "CO2": 400},
        )

    def test_observation_geometry_is_explicit_and_consistent(self) -> None:
        for case in self.cases:
            MODULE._validate_observation_geometry(case)

        # A copied LOS is nonzero and unit length, so this specifically proves
        # that the generator rejects a physically inconsistent sky-plane axis
        # rather than checking only that nine numerical fields were supplied.
        invalid = dict(self.cases[0])
        invalid["impact_x"] = invalid["los_x"]
        invalid["impact_y"] = invalid["los_y"]
        invalid["impact_z"] = invalid["los_z"]
        with self.assertRaisesRegex(ValueError, "not orthogonal"):
            MODULE._validate_observation_geometry(invalid)

    def test_generated_runtime_input_matches_tracked_file(self) -> None:
        with tempfile.TemporaryDirectory(prefix="amps-i3-input-test-") as directory:
            outputs = MODULE.generate(REPOSITORY, Path(directory))
            self.assertEqual([output.name for output in outputs], ["i3.in"])
            output = outputs[0]
            tracked = REPOSITORY / "input" / "i3.in"
            text = output.read_text(encoding="utf-8")
            self.assertEqual(text, tracked.read_text(encoding="utf-8"))
            self.assertNotIn("@", text)
            self.assertEqual(text.count("#begin case:"), 6)
            self.assertEqual(text.count("#end case"), 6)
            self.assertEqual(text.count("chemistry.rate."), 6 * 17)
            self.assertEqual(text.count("source.H2O.rate_s-1"), 6)
            self.assertEqual(text.count("mesh.surface_resolution_m"), 6)
            self.assertIn("brightness.enabled = true", text)
            self.assertIn(
                "brightness.impact_parameter_grid = LOG | 2500.0 | 1.0e9 | 200",
                text,
            )
            self.assertEqual(
                text.count("brightness.heliocentric_radial_velocity_km_s"), 6
            )
            self.assertEqual(text.count("brightness.observation_geometry"), 6)
            self.assertEqual(
                text.count("brightness.solar_lya_flux_1au_photon_cm-2_s-1"), 6
            )
            self.assertIn("H_LYA_IKL2020", text)
            self.assertIn("OH_AX00_SA1988_LOWELL_UNQUENCHED", text)
            self.assertIn("mesh.outer_radius_m = 1000000000", text)
            self.assertIn("mesh.radial_cell_count = 2000", text)
            self.assertIn("simulation.sampling_length = 10000", text)
            self.assertIn("simulation.max_iterations = 30000", text)

    def test_compile_configuration_is_layout_only(self) -> None:
        text = (REPOSITORY / "input" / "i3.conf").read_text(encoding="utf-8")
        self.assertIn("RuntimeInputMode=on", text)
        self.assertIn("MaxMeshRefinmentLevel=30", text)
        self.assertEqual(text.count("ChemistryReaction="), 17)
        self.assertEqual(text.count(" | HELIOCENTRIC | 0 |"), 17)
        self.assertIn("BrightnessMode=off", text)
        self.assertNotIn("2025-", text)
        self.assertNotIn("8.0e27", text)


if __name__ == "__main__":
    unittest.main()

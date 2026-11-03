import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "export_forward_sources", ROOT / "export_forward_sources.py")
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)


class ForwardSourceExportTests(unittest.TestCase):
    def manifest(self, mode="shape"):
        value = 1.0 / len(EXPORT.SOURCE_IDS) if mode == "shape" else 2.0e25
        return {
            "schema_version": 1,
            "source_case": "tenishev2014_hybrid_8",
            "distribution_name": "unit_fixture",
            "mode": mode,
            "species": {
                "H2O": {name: value for name in EXPORT.SOURCE_IDS},
                "CO2": {name: value for name in EXPORT.SOURCE_IDS},
            },
            "fit_provenance": {"strategy": "unit-test"},
        }

    def test_shape_export_has_complete_hybrid_geometry(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "manifest.json"
            manifest.write_text(json.dumps(self.manifest()), encoding="utf-8")
            text = EXPORT.render(manifest, EXPORT.DEFAULT_CONFIG)
        self.assertEqual(text.count("#begin source:"), 15)
        self.assertIn("#begin source: I_baghdad", text)
        self.assertIn("#begin source: damascus_extra", text)
        self.assertIn("mode=shape", text)
        self.assertIn("schema_version=2", text)
        self.assertIn(
            "stripe_surface_mapping=local_vertical_to_reference_sphere", text)
        self.assertIn("fit_provenance={\"strategy\":\"unit-test\"}", text)
        coordinates = []
        for line in text.splitlines():
            if line.startswith(("x_m=", "y_m=", "z_m=")):
                coordinates.append(float(line.split("=", 1)[1]))
        self.assertEqual(len(coordinates) % 3, 0)
        for offset in range(0, len(coordinates), 3):
            radius = math.sqrt(math.fsum(
                value * value for value in coordinates[offset:offset + 3]))
            self.assertAlmostEqual(radius, 252100.0, places=7)

        config = EXPORT.load_json(EXPORT.DEFAULT_CONFIG)
        stripe_path = (EXPORT.DEFAULT_CONFIG.parent /
                       config["stripe_geometry_file"]).resolve()
        raw = EXPORT.load_json(stripe_path)["stripes"]
        raw_coordinates = [
            point
            for item in config["tiger_stripes"]
            for point in raw[item["geometry_id"]]
        ]
        projected_coordinates = [coordinates[offset:offset + 3]
                                 for offset in range(0, len(coordinates), 3)]
        self.assertEqual(len(raw_coordinates), len(projected_coordinates))
        for source, projected in zip(raw_coordinates, projected_coordinates):
            source_radius = math.sqrt(math.fsum(value * value for value in source))
            projected_radius = math.sqrt(math.fsum(
                value * value for value in projected))
            for source_value, projected_value in zip(source, projected):
                self.assertAlmostEqual(source_value / source_radius,
                                       projected_value / projected_radius,
                                       places=15)

    def test_incomplete_strength_map_is_rejected(self):
        manifest = self.manifest()
        del manifest["species"]["H2O"][EXPORT.SOURCE_IDS[-1]]
        with self.assertRaisesRegex(ValueError, "canonical 15"):
            EXPORT.validate_manifest(manifest)

    def test_absolute_rates_are_not_normalized(self):
        name, mode, values = EXPORT.validate_manifest(self.manifest("absolute"))
        self.assertEqual(name, "unit_fixture")
        self.assertEqual(mode, "absolute")
        self.assertEqual(sum(values["H2O"].values()), 3.0e26)

    def test_e3_e5_equal_startup_manifests_preserve_reported_totals(self):
        expected = {
            "e3_equal_h2o_co2_startup_manifest.json": (
                "e3_h2o_co2_equal_startup",
                "enceladus_e3_equal_h2o_co2_sources.in",
                2.060909950673409e28,
                5.356e25,
                [1.0712e25, 2.678e25],
            ),
            "e5_equal_h2o_co2_startup_manifest.json": (
                "e5_h2o_co2_equal_startup",
                "enceladus_e5_equal_h2o_co2_sources.in",
                4.520015885781454e28,
                1.013e26,
                [2.026e25, 5.065e25],
            ),
        }
        for filename, (expected_name, receipt_name, h2o_total, co2_total,
                       physical_bracket) in expected.items():
            with self.subTest(filename=filename):
                manifest_path = ROOT / "results" / filename
                manifest = EXPORT.load_json(manifest_path)
                name, mode, values = EXPORT.validate_manifest(manifest)
                self.assertEqual(name, expected_name)
                self.assertEqual(mode, "absolute")
                self.assertEqual(set(values), {"H2O", "CO2"})
                self.assertEqual(len(set(values["H2O"].values())), 1)
                self.assertEqual(len(set(values["CO2"].values())), 1)
                self.assertEqual(math.fsum(values["H2O"].values()), h2o_total)
                self.assertEqual(math.fsum(values["CO2"].values()), co2_total)
                co2_provenance = manifest["fit_provenance"]["species"]["CO2"]
                self.assertEqual(
                    co2_provenance["impact_fragmentation_enhancement_factor_range"],
                    [2.0, 5.0])
                self.assertEqual(
                    co2_provenance["physical_co2_total_bracket_s-1"],
                    physical_bracket)
                receipt = ROOT.parents[1] / "input" / receipt_name
                self.assertEqual(
                    receipt.read_text(encoding="utf-8"),
                    EXPORT.render(manifest_path, EXPORT.DEFAULT_CONFIG))

    def test_unsafe_distribution_and_species_names_are_rejected(self):
        manifest = self.manifest()
        manifest["distribution_name"] = "bad: parser block"
        with self.assertRaisesRegex(ValueError, "distribution_name"):
            EXPORT.validate_manifest(manifest)
        manifest = self.manifest()
        manifest["species"]["H2O;bad"] = manifest["species"].pop("H2O")
        with self.assertRaisesRegex(ValueError, "species names"):
            EXPORT.validate_manifest(manifest)

    def test_duplicate_canonical_geometry_id_is_rejected(self):
        config = EXPORT.load_json(EXPORT.DEFAULT_CONFIG)
        stripe_path = (EXPORT.DEFAULT_CONFIG.parent /
                       config["stripe_geometry_file"]).resolve()
        config["stripe_geometry_file"] = str(stripe_path)
        config["discrete_jets"][1]["id"] = config["discrete_jets"][0]["id"]
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = Path(directory) / "manifest.json"
            config_path = Path(directory) / "config.json"
            manifest_path.write_text(json.dumps(self.manifest()), encoding="utf-8")
            config_path.write_text(json.dumps(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical jet"):
                EXPORT.render(manifest_path, config_path)

    def test_nonfinite_kinematics_are_rejected_before_rendering(self):
        config = EXPORT.load_json(EXPORT.DEFAULT_CONFIG)
        stripe_path = (EXPORT.DEFAULT_CONFIG.parent /
                       config["stripe_geometry_file"]).resolve()
        config["stripe_geometry_file"] = str(stripe_path)
        config["discrete_jets"][0]["bulk_speed_m_s"] = "nan"
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = Path(directory) / "manifest.json"
            config_path = Path(directory) / "config.json"
            manifest_path.write_text(json.dumps(self.manifest()), encoding="utf-8")
            config_path.write_text(json.dumps(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                EXPORT.render(manifest_path, config_path)


if __name__ == "__main__":
    unittest.main()

import json
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fitting import (SOURCE_CASES, append_extra_jet, load_config,  # noqa: E402
                     sources_from_config)


EXPECTED_41 = [
    1, 4, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 21, 24, 25,
    27, 28, 32, 35, 36, 37, 40, 54, 60, 61, 62, 63, 66, 67, 68,
    73, 84, 85, 86, 89, 90, 91, 92, 93, 94,
]


class SourceCaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        name = os.environ.get("ENCELADUS_JET_CATALOG")
        if not name and os.environ.get("ENCELADUS_VALIDATION_DATA_ROOT"):
            name = str(Path(os.environ["ENCELADUS_VALIDATION_DATA_ROOT"])
                       / "data" / "derived" / "porco2014_jet_catalog.json")
        if not name:
            raise unittest.SkipTest(
                "set ENCELADUS_VALIDATION_DATA_ROOT or ENCELADUS_JET_CATALOG")
        cls.catalog_path = Path(name)
        cls.catalog = json.loads(cls.catalog_path.read_text())
        cls.config_path = ROOT/"config"/"tenishev2014_stage0.json"

    def test_all_three_source_cases_are_declared(self):
        self.assertEqual(SOURCE_CASES, (
            "tenishev2014_hybrid_8",
            "porco2014_98",
            "portyankina2022_solar_41",
        ))

    def test_legacy_case_retains_eight_jets_and_seven_stripes(self):
        config = load_config(self.config_path)
        sources, labels = sources_from_config(config)
        self.assertEqual(config["source_case"], "tenishev2014_hybrid_8")
        self.assertEqual(len(config["discrete_jets"]), 8)
        self.assertEqual(len(config["tiger_stripes"]), 7)
        self.assertEqual(len(sources), 15)
        self.assertEqual(len(labels), 17)

    def test_porco_case_has_all_98_and_no_stripes(self):
        config = load_config(
            self.config_path, "porco2014_98", self.catalog_path)
        sources, labels = sources_from_config(config)
        self.assertEqual([item["porco_id"] for item in config["discrete_jets"]],
                         list(range(1, 99)))
        self.assertEqual(config["tiger_stripes"], [])
        self.assertEqual(len(sources), 98)
        self.assertEqual(len(labels), 100)  # jets plus two empirical backgrounds
        first = sources[0]
        self.assertAlmostEqual(
            float(np.dot(first.normal, first.axis)), math.cos(math.radians(28)),
            places=12)

    def test_provenance_tracked_extra_jet_makes_99_source_variant(self):
        config = load_config(
            self.config_path, "porco2014_98", self.catalog_path)
        document = {
            "schema_version": 1,
            "frame": "IAU_ENCELADUS",
            "base_source_case": "porco2014_98",
            "variant_name": "porco2014_98_plus_test_candidate",
            "source": {
                "id": "test_candidate", "latitude_deg": -88.0,
                "west_longitude_deg": 296.0, "tilt_deg": 0.0,
                "azimuth_deg": 0.0, "temperature_K": 180.0,
                "bulk_speed_m_s": 1500.0,
                "source_law": "tenishev2010_point",
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"candidate.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            append_extra_jet(config, path)
        sources, labels = sources_from_config(config)
        self.assertEqual(config["source_variant"], document["variant_name"])
        self.assertEqual(config["resolved_source_case"]["jet_count"], 99)
        self.assertEqual(len(sources), 99)
        self.assertIn("jet:test_candidate", labels)

    def test_packaged_e17_candidate_is_opt_in_and_makes_99_sources(self):
        config = load_config(
            self.config_path, "porco2014_98", self.catalog_path)
        self.assertEqual(len(config["discrete_jets"]), 98)
        append_extra_jet(config, ROOT/"e17_extra_jet_candidate.json")
        sources, labels = sources_from_config(config)
        self.assertEqual(len(sources), 99)
        self.assertIn("jet:e17_candidate_baghdad", labels)
        self.assertEqual(config["source_variant"],
                         "porco2014_98_plus_e17_candidate")

    def test_portyankina_case_uses_exact_2010_selection_and_vertical_axes(self):
        config = load_config(
            self.config_path, "portyankina2022_solar_41", self.catalog_path)
        sources, labels = sources_from_config(config)
        self.assertEqual([item["porco_id"] for item in config["discrete_jets"]],
                         EXPECTED_41)
        self.assertEqual(config["tiger_stripes"], [])
        self.assertEqual(len(sources), 41)
        self.assertEqual(len(labels), 43)
        for source in sources:
            np.testing.assert_allclose(source.axis, source.normal, atol=2e-15)
            self.assertEqual(source.temperature_K, 180.0)
            self.assertEqual(source.bulk_speed_m_s, 1500.0)

    def test_external_catalog_has_expected_provenance_and_distributions(self):
        self.assertEqual(self.catalog["jet_count"], 98)
        self.assertEqual(self.catalog["selections"]["portyankina2022_solar_41"],
                         EXPECTED_41)
        all_counts = {name: 0 for name in
                      ("cairo", "baghdad", "damascus", "alexandria")}
        selected_counts = dict(all_counts)
        selected = set(EXPECTED_41)
        for item in self.catalog["jets"]:
            all_counts[item["sulcus"]] += 1
            if item["porco_id"] in selected:
                selected_counts[item["sulcus"]] += 1
        self.assertEqual(all_counts, {
            "cairo": 26, "baghdad": 37, "damascus": 29, "alexandria": 6})
        self.assertEqual(selected_counts, {
            "cairo": 16, "baghdad": 9, "damascus": 11, "alexandria": 5})


if __name__ == "__main__":
    unittest.main()

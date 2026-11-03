"""Validate the checked-in E7/E21 observation release without raw archives."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import unittest


OBSERVATIONS = Path(__file__).resolve().parents[1] / "observations"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ForwardObservationReleaseTest(unittest.TestCase):
    def test_manifest_products_are_portable_and_self_consistent(self) -> None:
        manifest = json.loads(
            (OBSERVATIONS / "release_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema_version"], 1)
        self.assertIn("must not be scored", manifest["contract"])
        self.assertEqual(len(manifest["products"]), 3)

        expected_rows = {
            "e7_co2_density.csv": 43,
            "e21_h2o_osnb_counts.csv": 704,
            "e21_co2_osnb_counts.csv": 481,
        }
        for product in manifest["products"]:
            relative = Path(product["path"])
            metadata_relative = Path(product["metadata"])
            self.assertFalse(relative.is_absolute())
            self.assertFalse(metadata_relative.is_absolute())
            self.assertEqual(relative.parent, Path("."))
            self.assertEqual(metadata_relative.parent, Path("."))
            table = OBSERVATIONS / relative
            metadata = OBSERVATIONS / metadata_relative
            self.assertEqual(product["rows"], expected_rows[relative.name])
            self.assertEqual(product["sha256"], sha256(table))
            self.assertTrue(metadata.is_file())

    def test_tables_have_lf_records_and_physical_geometry(self) -> None:
        for path in sorted(OBSERVATIONS.glob("*.csv")):
            raw = path.read_bytes()
            self.assertNotIn(b"\r", raw)
            self.assertTrue(raw.endswith(b"\n"))
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertTrue(rows)
            times = [float(row["time_ca_s"]) for row in rows]
            self.assertTrue(all(b > a for a, b in zip(times, times[1:])))
            for row in rows:
                radius = math.sqrt(sum(float(row[name]) ** 2 for name in (
                    "spacecraft_x_m", "spacecraft_y_m", "spacecraft_z_m")))
                self.assertGreater(radius, 252100.0)

            if path.name == "e7_co2_density.csv":
                self.assertTrue(all(float(row["density_m3"]) >= 0.0 for row in rows))
                self.assertTrue(all(float(row["sigma_counting_m3"]) > 0.0
                                    for row in rows))
            else:
                self.assertTrue(all(float(row["counts_per_integration"]) >= 0.0
                                    for row in rows))
                self.assertTrue(all(float(row["poisson_sigma_counts"]) > 0.0
                                    for row in rows))


if __name__ == "__main__":
    unittest.main()

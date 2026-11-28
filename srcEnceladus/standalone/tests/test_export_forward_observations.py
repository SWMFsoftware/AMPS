"""Validate the checked-in E7/E21 observation release without raw archives."""

from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
OBSERVATIONS = ROOT / "observations"
sys.path.insert(0, str(ROOT))

from export_forward_observations import parse_arguments  # noqa: E402
from validation_paths import VALIDATION_ROOT_ENV  # noqa: E402


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


class ExporterValidationRootTest(unittest.TestCase):
    """The data root is an explicit input; no machine path is assumed."""

    def parse_without_env(self, argv: list[str]):
        environment = {key: value for key, value in os.environ.items()
                       if key != VALIDATION_ROOT_ENV}
        with patch.dict(os.environ, environment, clear=True):
            return parse_arguments(argv)

    def assert_usage_error(self, argv: list[str]) -> str:
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), \
                self.assertRaises(SystemExit) as raised:
            self.parse_without_env(argv)
        self.assertEqual(raised.exception.code, 2)
        return stderr.getvalue()

    def test_missing_root_is_a_usage_error(self) -> None:
        message = self.assert_usage_error([])
        self.assertIn("--validation-data-root", message)
        self.assertIn(VALIDATION_ROOT_ENV, message)

    def test_nonexistent_root_is_a_usage_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            missing = Path(temporary) / "absent"
            message = self.assert_usage_error(
                ["--validation-data-root", str(missing)])
        self.assertIn("not a directory", message)

    def test_environment_root_and_explicit_precedence(self) -> None:
        with tempfile.TemporaryDirectory() as first, \
                tempfile.TemporaryDirectory() as second:
            with patch.dict(os.environ, {VALIDATION_ROOT_ENV: first}):
                self.assertEqual(parse_arguments([]).validation_data_root,
                                 Path(first).resolve())
                explicit = parse_arguments(["--validation-data-root", second])
                self.assertEqual(explicit.validation_data_root,
                                 Path(second).resolve())

    def test_original_option_spelling_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = self.parse_without_env(["--validation-root", temporary])
        self.assertEqual(args.validation_data_root, Path(temporary).resolve())


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression tests for the extracted MAVEN/IUVS brightness references."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
VISUALIZATION = ROOT / "srcHartley2/visualization"
REFERENCE_DIR = VISUALIZATION / "references/i3_brightness_MAVEN"
REFERENCE_TABLE = REFERENCE_DIR / "maven_iuvs_brightness.csv"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


PLOT = load_module("i3_brightness_plot", VISUALIZATION / "plot_i3_brightness_comparison.py")
EXTRACT = load_module("i3_brightness_extract", VISUALIZATION / "extract_i3_brightness_maven.py")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class I3BrightnessReferenceTest(unittest.TestCase):
    def test_all_case_species_subsets_are_complete_and_valid(self) -> None:
        expected_counts = {
            "H": (11, 14, 14, 7, 8, 9),
            "OH": (4, 8, 5, 3, 2, 2),
        }
        total = 0
        for species, counts in expected_counts.items():
            for case_id, expected in enumerate(counts, start=1):
                table, values = PLOT.read_reference(REFERENCE_DIR, case_id, species)
                self.assertEqual(table, REFERENCE_TABLE)
                self.assertEqual(values["x"].size, expected)
                self.assertTrue((values["x_lower"] < values["x"]).all())
                self.assertTrue((values["x"] < values["x_upper"]).all())
                self.assertTrue((values["y"] > 0.0).all())
                total += expected
        self.assertEqual(total, 87)

    def test_metadata_authenticates_committed_output(self) -> None:
        metadata = json.loads((REFERENCE_DIR / "metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["schema_version"], 1)
        self.assertEqual(metadata["row_count"], 87)
        self.assertEqual(metadata["output_table"], REFERENCE_TABLE.name)
        self.assertEqual(metadata["output_table_sha256"], sha256(REFERENCE_TABLE))
        self.assertEqual(
            metadata["anchor_table_sha256"], sha256(REFERENCE_DIR / "digitization_pixels.csv")
        )
        self.assertEqual(len(metadata["panels"]), 12)

    def test_alternate_reference_directory(self) -> None:
        with tempfile.TemporaryDirectory(prefix="i3-brightness-reference-") as temporary:
            alternate = Path(temporary)
            shutil.copy2(REFERENCE_TABLE, alternate / REFERENCE_TABLE.name)
            _, original = PLOT.read_reference(REFERENCE_DIR, 1, "H")
            _, copied = PLOT.read_reference(alternate, 1, "H")
            for name in original:
                self.assertTrue(np.array_equal(original[name], copied[name], equal_nan=True))

    def test_archive_regeneration_is_byte_exact(self) -> None:
        # Keep the reviewed raster source beside the test that verifies its
        # deterministic extraction. Do not fall back to a repository-root
        # developer copy: omitting the tracked artifact must fail on every
        # machine instead of becoming a machine-dependent pass or skip.
        archive = Path(__file__).resolve().parent / "data" / "i3" / "Archive.zip"
        self.assertTrue(archive.is_file(), f"missing bundled I3 archive: {archive}")
        with tempfile.TemporaryDirectory(prefix="i3-brightness-extract-") as temporary:
            destination = Path(temporary)
            status = EXTRACT.main(["--archive", str(archive), "--output-dir", str(destination)])
            self.assertEqual(status, 0)
            self.assertEqual(
                (destination / REFERENCE_TABLE.name).read_bytes(), REFERENCE_TABLE.read_bytes()
            )
            self.assertEqual(
                (destination / "metadata.json").read_bytes(),
                (REFERENCE_DIR / "metadata.json").read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()

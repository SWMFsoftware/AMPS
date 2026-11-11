#!/usr/bin/env python3
"""Synthetic fail-closed tests for Stage-9 ensemble comparison publication."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


THIS_DIR = Path(__file__).resolve().parent
TESTS_DIR = THIS_DIR.parent / "tests"
sys.path.insert(0, str(THIS_DIR))
sys.path.insert(0, str(TESTS_DIR))

import generate_stage9_rubin_comparison as comparison  # noqa: E402
from test_stage9_science_gate_selftest import Archive, SEEDS  # noqa: E402


def write_reference(root: Path) -> tuple[Path, Path]:
    data = root / "paper.csv"
    metadata = root / "paper.metadata.json"
    fields = (
        "figure",
        "pdf_page",
        "species",
        "variant",
        "quantity",
        "radius_km",
        "value",
        "unit",
    )
    rows = []
    curves = []
    for figure, species, values in (
        (8, "H2O", ((7.0, 1.0e8), (17.0, 5.0e7), (107.0, 1.0e7))),
        (9, "C2", ((7.0, 1.0e2), (17.0, 5.0e1), (107.0, 1.0e1))),
    ):
        for radius, value in values:
            rows.append(
                {
                    "figure": figure,
                    "pdf_page": figure + 1,
                    "species": species,
                    "variant": "nominal",
                    "quantity": "number_density",
                    "radius_km": radius,
                    "value": value,
                    "unit": "cm^-3",
                }
            )
        curves.append(
            {
                "figure": figure,
                "species": species,
                "variant": "nominal",
                "quantity": "number_density",
                "samples": len(values),
            }
        )
    with data.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    metadata.write_text(
        json.dumps(
            {
                "paper": "synthetic Rubin fixture; not science evidence",
                "source_pdf_sha256": "a" * 64,
                "method": "synthetic",
                "curves": curves,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return data, metadata


class Stage9ComparisonTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.archive = Archive(self.root / "archive")
        self.paper, self.paper_metadata = write_reference(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_publishes_all_formats_with_ensemble_metadata(self) -> None:
        output = self.root / "comparison"
        product = comparison.publish(
            self.archive.manifest_path, output, self.paper, self.paper_metadata, dpi=72
        )
        base = output / "stage9_rubin2011_ensemble_comparison"
        for suffix in (".png", ".eps", ".pdf", ".csv", ".json"):
            path = base.with_suffix(suffix)
            self.assertTrue(path.is_file(), path)
            self.assertGreater(path.stat().st_size, 10)
        self.assertEqual(product["schema"], comparison.SCHEMA)
        self.assertEqual(product["seeds"], list(SEEDS))
        self.assertEqual(product["run_count"], 8)
        self.assertTrue(product["convergence"]["all_runs_converged_before_hard_limit"])
        self.assertEqual(len(product["effective_exposure"]["runs"]), 8)
        self.assertEqual(len(product["profiles"]["C2"]["RUBIN2011"]), 1)
        with base.with_suffix(".csv").open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual({row["source"] for row in rows}, {"ANALOG", "FORCED", "RUBIN2011"})

    def test_altered_summary_is_rejected_before_output(self) -> None:
        record = self.archive.manifest["runs"][0]
        path = self.archive.root / record["summary"]
        path.write_text(path.read_text(encoding="utf-8") + " ", encoding="utf-8")
        output = self.root / "comparison"
        with self.assertRaisesRegex(comparison.ComparisonError, "hash mismatch"):
            comparison.publish(
                self.archive.manifest_path, output, self.paper, self.paper_metadata
            )
        self.assertFalse((output / "stage9_rubin2011_ensemble_comparison.json").exists())

    def test_failed_science_is_not_published(self) -> None:
        self.archive.mutate_summary(
            "ANALOG",
            SEEDS[0],
            lambda value: value["operations"].__setitem__("emergency_events", 1),
        )
        with self.assertRaisesRegex(comparison.ComparisonError, "did not pass"):
            comparison.publish(
                self.archive.manifest_path,
                self.root / "comparison",
                self.paper,
                self.paper_metadata,
            )

    def test_reference_inventory_mismatch_is_rejected(self) -> None:
        metadata = json.loads(self.paper_metadata.read_text(encoding="utf-8"))
        metadata["curves"].pop()
        self.paper_metadata.write_text(json.dumps(metadata), encoding="utf-8")
        with self.assertRaisesRegex(comparison.ComparisonError, "inventories differ"):
            comparison.publish(
                self.archive.manifest_path,
                self.root / "comparison",
                self.paper,
                self.paper_metadata,
            )

    def test_existing_evidence_is_never_overwritten(self) -> None:
        output = self.root / "comparison"
        output.mkdir()
        (output / "stage9_rubin2011_ensemble_comparison.json").write_text("keep\n")
        with self.assertRaisesRegex(comparison.ComparisonError, "overwrite"):
            comparison.publish(
                self.archive.manifest_path, output, self.paper, self.paper_metadata
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)

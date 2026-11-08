#!/usr/bin/env python3
"""Verify provenance, normalized grids, and generated core reference data."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import unittest


MODULE = Path(__file__).resolve().parents[1]
REFERENCE = MODULE / "reference" / "g_factors"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class ReferenceDataTest(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads((REFERENCE / "manifest.json").read_text())

    def test_manifest_hashes_match_every_authoritative_file(self) -> None:
        models = self.manifest["models"]
        h = models["H_LYA_IKL2020"]
        oh = models["OH_AX00_SA1988_LOWELL_UNQUENCHED"]
        self.assertEqual(sha256(REFERENCE / h["source_file"]), h["source_sha256"])
        self.assertEqual(sha256(REFERENCE / oh["raw_file"]), oh["raw_sha256"])
        self.assertEqual(
            sha256(REFERENCE / oh["normalized_file"]), oh["normalized_sha256"]
        )
        generated = self.manifest["generated"]
        self.assertEqual(
            sha256((REFERENCE / generated["file"]).resolve()), generated["sha256"]
        )

    def test_lowell_table_is_complete_and_not_extrapolated(self) -> None:
        oh = self.manifest["models"]["OH_AX00_SA1988_LOWELL_UNQUENCHED"]
        with (REFERENCE / oh["normalized_file"]).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 1201)
        velocity = [float(row["heliocentric_radial_velocity_km_s"]) for row in rows]
        self.assertEqual((velocity[0], velocity[-1]), (-60.0, 60.0))
        self.assertTrue(
            all(abs(right - left - 0.1) < 1.0e-10 for left, right in zip(velocity, velocity[1:]))
        )
        self.assertTrue(all(row["api_extrapolated"] == "false" for row in rows))
        self.assertTrue(
            all(float(row["oh00_energy_g_1au_erg_s-1_molecule-1"]) > 0 for row in rows)
        )
        with (REFERENCE / oh["raw_file"]).open() as stream:
            raw = [json.loads(line) for line in stream]
        self.assertEqual(len(raw), 1201)
        self.assertTrue(all(row["response"]["extrapolated"] is False for row in raw))

    def test_generated_include_is_current(self) -> None:
        subprocess.run(
            [sys.executable, str(MODULE / "tools/generate_reference_tables.py"), "--check"],
            check=True,
        )


if __name__ == "__main__":
    unittest.main()

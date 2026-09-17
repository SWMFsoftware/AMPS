#!/usr/bin/env python3
"""Self-tests for the Phase-V external-evidence runner.

Temporary bundles exercise the parser and checksum boundary only.  They prove
that the campaign machinery classifies valid, missing, and corrupt evidence
correctly; they are never registered as scientific event evidence themselves.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "validation" / "run_validation.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ValidationRunnerTests(unittest.TestCase):
    def run_runner(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run([sys.executable, str(RUNNER), *arguments],
                              cwd=str(ROOT), text=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, check=False)

    def test_registry_lists_all_evidence_classes(self) -> None:
        completed = self.run_runner("--list")
        self.assertEqual(completed.returncode, 0, completed.stdout)
        for case_id in ("NAT3D01", "MPI3D01", "XM3D01", "XM3D05",
                        "OV3D01", "OV3D04"):
            self.assertIn(case_id, completed.stdout)

    def test_missing_external_evidence_is_skip(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            completed = self.run_runner("--case", "XM3D01", "--output-dir",
                                        str(output))
            self.assertEqual(completed.returncode, 0, completed.stdout)
            result = json.loads((output / "XM3D01" / "result.json").read_text())
            self.assertEqual(result["status"], "SKIP")

    def make_series_bundle(self, root: Path) -> Path:
        case_dir = root / "XM3D01"
        case_dir.mkdir(parents=True)
        model = case_dir / "model.csv"
        reference = case_dir / "reference.csv"
        content = ("coordinate_si,value_si\n"
                   "0,1\n"
                   "1,10\n"
                   "2,100\n"
                   "3,10\n"
                   "4,1\n")
        model.write_text(content, encoding="utf-8")
        reference.write_text(content, encoding="utf-8")
        manifest = {
            "schema": "srcsep3d-scientific-evidence-v1",
            "case_id": "XM3D01",
            "coordinate_name": "time",
            "coordinate_units": "s",
            "value_name": "intensity",
            "value_units": "m^-2 s^-1 sr^-1 J^-1",
            "model_csv": model.name,
            "model_sha256": sha256(model),
            "reference_csv": reference.name,
            "reference_sha256": sha256(reference),
            "provenance": {
                "model": "unit-test fixture for parser mechanics",
                "reference": "independent unit-test fixture for parser mechanics",
            },
            "exact_pairs": [],
        }
        (case_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return model

    def test_valid_bundle_passes_and_hash_change_errors(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            evidence = root / "evidence"
            model = self.make_series_bundle(evidence)
            output = root / "pass"
            completed = self.run_runner(
                "--case", "XM3D01", "--evidence-root", str(evidence),
                "--output-dir", str(output))
            self.assertEqual(completed.returncode, 0, completed.stdout)
            result = json.loads((output / "XM3D01" / "result.json").read_text())
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["metrics"]["coverage"], 1.0)

            # Change one byte after the manifest was signed.  The second run
            # must be ERROR before any scientific metric is calculated.
            model.write_text(model.read_text(encoding="utf-8") + "5,1\n",
                             encoding="utf-8")
            corrupt_output = root / "corrupt"
            corrupt = self.run_runner(
                "--case", "XM3D01", "--evidence-root", str(evidence),
                "--output-dir", str(corrupt_output))
            self.assertEqual(corrupt.returncode, 2, corrupt.stdout)
            result = json.loads(
                (corrupt_output / "XM3D01" / "result.json").read_text())
            self.assertEqual(result["status"], "ERROR")
            self.assertIn("checksum mismatch", result["message"])

    def test_convergence_bundle_reports_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "evidence" / "XM3D05"
            case_dir.mkdir(parents=True)
            manifest = {
                "schema": "srcsep3d-convergence-evidence-v1",
                "case_id": "XM3D05",
                "quantity": "L1 error",
                "levels": [
                    {"resolution": 4.0, "error": 0.16},
                    {"resolution": 2.0, "error": 0.04},
                    {"resolution": 1.0, "error": 0.01},
                ],
                "default_relative_difference": 0.04,
                "provenance": {
                    "model": "unit-test convergence fixture",
                    "reference": "unit-test finest-grid fixture",
                },
            }
            (case_dir / "manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            output = root / "output"
            completed = self.run_runner(
                "--case", "XM3D05", "--evidence-root",
                str(root / "evidence"), "--output-dir", str(output))
            self.assertEqual(completed.returncode, 0, completed.stdout)
            result = json.loads((output / "XM3D05" / "result.json").read_text())
            self.assertEqual(result["status"], "PASS")
            self.assertAlmostEqual(result["metrics"]["minimum_observed_order"], 2.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)

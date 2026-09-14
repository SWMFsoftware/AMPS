#!/usr/bin/env python3
"""Focused tests for the dependency-light Python test orchestrator.

The fixture deliberately contains both supported analytical evidence forms:
an ordinary registry metric and a CSV solution series.  The analytical values
are literal independent fixture data, not values computed by a srcSEP kernel,
so this test can detect plotting/report regressions without validating a
production numerical routine against itself.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "test" / "run_tests.py"


def _load_runner_module():
    """Load the runner by path without requiring test/ to be a package."""
    specification = importlib.util.spec_from_file_location(
        "srcsep_python_test_runner", RUNNER)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot import runner from {RUNNER}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class PythonTestRunnerTests(unittest.TestCase):
    def test_help_contains_annotated_selection_and_plot_examples(self):
        """Keep the command-line quick-start available to archive users."""
        completed = subprocess.run(
            [sys.executable, str(RUNNER), "--help"], cwd=str(ROOT), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        self.assertEqual(completed.returncode, 0, completed.stdout)
        # These anchors cover the distinct native, source-only, report-only,
        # MPI, plotting, and forwarded-model-argument workflows.  Exact spacing
        # is intentionally not asserted because argparse may wrap descriptions
        # differently when a Python distribution changes terminal width.
        for expected in (
                "Discover the tests registered",
                "--test PARK07",
                "--group parker --group fte-dmumu",
                "--routine",
                "--all",
                "--suite controlled-analytical",
                "--mpi-np 4",
                "--from-json results.json",
                "--plot none",
                "--mover parker --turbulence-source prescribed",
                "Outputs:"):
            self.assertIn(expected, completed.stdout)

    def test_registry_list_parser_ignores_headers_and_deduplicates(self):
        runner = _load_runner_module()
        listing = """Available tests:\nPARK01 | parker | deterministic\nFTED08 | fte-dmumu | extended\nPARK01 | parker | duplicate\n"""
        self.assertEqual(runner._parse_list_output(listing),
                         ["FTED08", "PARK01"])

    def test_existing_report_produces_png_eps_and_manifests(self):
        with tempfile.TemporaryDirectory(prefix="srcsep-python-runner-") as tmp:
            temporary = Path(tmp)
            solution_csv = temporary / "FTED08_solution.csv"
            # Conventional column names exercise the pointwise overlay path.
            # Values are fixed fixture evidence for y=x^2 at three coordinates.
            solution_csv.write_text(
                "x,numerical,analytical\n"
                "0,0,0\n"
                "1,1.02,1\n"
                "2,3.96,4\n",
                encoding="utf-8")
            report = {
                "schema": "srcsep-component-tests-v1",
                "exit_code": 0,
                "totals": {"passed": 2, "failed": 0,
                           "skipped": 0, "errors": 0},
                "results": [
                    {
                        "id": "PARK01",
                        "status": "PASS",
                        "message": "controlled analytical moment",
                        "elapsed_seconds": 0.01,
                        "seed": 1234,
                        "configuration": ["fixture=true"],
                        "metrics": [{
                            "name": "relative_error",
                            "value": 0.01,
                            "tolerance": 0.05,
                            "comparison": "relative_error <= tolerance",
                            "units": "1",
                        }],
                        "artifacts": [],
                    },
                    {
                        "id": "FTED08",
                        "status": "PASS",
                        "message": "controlled solution series",
                        "elapsed_seconds": 0.01,
                        "seed": 5678,
                        "configuration": ["fixture=true"],
                        "metrics": [],
                        "artifacts": [solution_csv.name],
                    },
                ],
            }
            report_path = temporary / "fixture-results.json"
            report_path.write_text(json.dumps(report), encoding="utf-8")
            output = temporary / "output"

            completed = subprocess.run(
                [sys.executable, str(RUNNER), "--from-json", str(report_path),
                 "--output-dir", str(output), "--formats", "png,eps"],
                cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=False)
            self.assertEqual(completed.returncode, 0, completed.stdout)

            manifest = json.loads(
                (output / "analytical_plot_manifest.json").read_text(
                    encoding="utf-8"))
            by_id = {entry["test_id"]: entry for entry in manifest["plots"]}
            self.assertEqual(by_id["PARK01"]["kind"], "metric-acceptance")
            self.assertEqual(by_id["FTED08"]["kind"], "solution-series")
            for test_id in ("PARK01", "FTED08"):
                for suffix in ("png", "eps"):
                    figure = output / "plots" / f"{test_id}_comparison.{suffix}"
                    self.assertTrue(figure.is_file(), figure)
                    self.assertGreater(figure.stat().st_size, 0)

            run_manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(run_manifest["exit_code"], 0)
            self.assertTrue(run_manifest["report_sha256"])


if __name__ == "__main__":
    unittest.main(verbosity=2)

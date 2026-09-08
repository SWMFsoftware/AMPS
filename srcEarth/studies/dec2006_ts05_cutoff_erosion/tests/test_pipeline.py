"""Small integration tests for generated inputs and comparison normalization."""

from __future__ import annotations

import csv
import contextlib
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from compare_observations import normalize_pamela, normalize_poes
from run_study import execute, independent_validation_remains
import run_study
from study_common import read_driver


class PipelineTests(unittest.TestCase):
    def test_independent_observation_validation_is_not_suppressed(self):
        """A PAMELA failure must not prevent the independent POES check."""

        self.assertTrue(independent_validation_remains(
            "pamela", ["poes", "morphology", "compare"]
        ))
        self.assertFalse(independent_validation_remains(
            "poes", ["morphology", "compare"]
        ))
        self.assertFalse(independent_validation_remains(
            "validate", ["pamela", "poes"]
        ))

    def test_pamela_failure_allows_poes_then_blocks_morphology(self):
        """Default policy collects both observation checks before stopping."""

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "study"

            def fake_execute(command, cwd, log_path, dry_run, **keywords):
                return 1 if keywords["stage"] == "pamela" else 0

            arguments = [
                "run_study.py", "--prepare-only", "--output-root", str(output),
                "--stage", "validate", "--stage", "pamela",
                "--stage", "poes", "--stage", "morphology",
            ]
            with mock.patch.object(sys, "argv", arguments), \
                    mock.patch.object(run_study, "execute", side_effect=fake_execute), \
                    contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                return_code = run_study.main()

            self.assertEqual(return_code, 1)
            manifest = json.loads((output / "study_run_manifest.json").read_text())
            self.assertEqual(manifest["return_codes"]["validate"], 0)
            self.assertEqual(manifest["return_codes"]["pamela"], 1)
            self.assertEqual(manifest["return_codes"]["poes"], 0)
            self.assertNotIn("morphology", manifest["return_codes"])

    def test_stage_execution_tees_live_output_and_records_status(self):
        """The orchestrator must show child output without sacrificing logs."""

        with tempfile.TemporaryDirectory() as temporary:
            log_path = Path(temporary) / "logs" / "diagnostic.log"
            terminal = io.StringIO()
            command = [
                sys.executable, "-c",
                "print('synthetic AMPS progress', flush=True)",
            ]
            with contextlib.redirect_stdout(terminal):
                return_code = execute(
                    command, ROOT, log_path, False,
                    stage="pamela", stage_index=2, stage_count=7,
                )

            self.assertEqual(return_code, 0)
            screen = terminal.getvalue()
            self.assertIn("[2/7] START PAMELA", screen)
            self.assertIn("[PAMELA] synthetic AMPS progress", screen)
            self.assertIn("[2/7] PASS PAMELA", screen)
            saved = log_path.read_text()
            self.assertIn("Stage: pamela", saved)
            self.assertIn("synthetic AMPS progress", saved)

    def test_normalized_comparison_schema(self):
        pamela = normalize_pamela([{
            "interval_midpoint_utc": "2006-12-14T00:00:00Z",
            "rigidity_center_gv": "0.5", "pamela_cutoff_aacgm_deg": "60",
            "amps_cutoff_aacgm_deg": "59", "pamela_sigma_plus_deg": "0.5",
            "pamela_sigma_minus_deg": "0.6",
        }])[0]
        self.assertEqual(pamela["model_minus_observation_deg"], -1.0)
        self.assertTrue(pamela["used_for_primary_metrics"])
        poes = normalize_poes([{
            "interval_midpoint_utc": "2006-12-14T00:00:00Z",
            "rigidity_gv": "0.174013525", "channel": "P6", "hemisphere": "N",
            "mlt_hour": "3", "observed_boundary_aacgm_deg": "62",
            "modeled_boundary_aacgm_deg": "63", "used_for_acceptance": "True",
            "validation_role": "PRIMARY", "sigma_deg": "1",
        }])[0]
        self.assertEqual(poes["model_minus_observation_deg"], 1.0)
        self.assertTrue(poes["used_for_primary_metrics"])

    def test_smoke_prepare_renders_eight_parser_safe_cases(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "morphology"
            completed = subprocess.run([
                sys.executable, str(ROOT / "scripts" / "run_morphology.py"),
                "--profile", "SMOKE", "--prepare-only", "--output-root", str(output),
            ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            result = json.loads((output / "morphology_result.json").read_text())
            self.assertEqual(result["n_cases"], 8)
            inputs = list(output.rglob("AMPS_PARAM_C10.in"))
            self.assertEqual(len(inputs), 8)
            for path in inputs:
                text = path.read_text()
                self.assertIn("CUTOFF_SAMPLING         VERTICAL", text)
                self.assertIn("CUTOFF_SEARCH_ALGORITHM RIGIDITY_LIST", text)
                self.assertNotIn("CUTOFF_UNRESOLVED_EXTENSION_PASSES", text)

    def test_sensitivity_driver_generation_preserves_cadence(self):
        with tempfile.TemporaryDirectory() as temporary:
            completed = subprocess.run([
                sys.executable,
                str(ROOT / "scripts" / "make_ts05_sensitivity_drivers.py"),
                "--output-root", temporary,
            ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            for name in ("history_frozen", "instantaneous_frozen"):
                rows = read_driver(Path(temporary) / f"ts05_dec2006_{name}.txt")
                self.assertEqual(len(rows), 865)
                self.assertEqual((rows[1].epoch - rows[0].epoch).total_seconds(), 300.0)


if __name__ == "__main__":
    unittest.main()

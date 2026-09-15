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
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "test" / "run_tests.py"
CASE_RUNNER = ROOT / "validation" / "run_case.py"


def _load_runner_module():
    """Load the runner by path without requiring test/ to be a package."""
    specification = importlib.util.spec_from_file_location(
        "srcsep_python_test_runner", RUNNER)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot import runner from {RUNNER}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _load_cross_model_runner():
    """Load the case scorer while preserving its sibling-module imports."""
    case_directory = ROOT / "validation" / "cases"
    sys.path.insert(0, str(case_directory))
    try:
        specification = importlib.util.spec_from_file_location(
            "srcsep_cross_model_case_runner",
            case_directory / "cross_model_case_runner.py")
        if specification is None or specification.loader is None:
            raise RuntimeError("cannot import cross-model case runner")
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


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
                "--amps ../amps --validation-case CV01",
                "--validation-case XM01 --validation-case XM02",
                "validation/cases/XM02/publication_input.json",
                "validation/cases/XM03/publication_input.json",
                "validation/cases/XM03/input/earth_shock_thermal_source.csv",
                "validation/cases/XM03/reference/liu_figure12_earth_observations.csv",
                "Do not add --case-input",
                "No external production CSV is used by XM02",
                "XM03 needs no external model CSV",
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
        # This reproduces the real C++ registry preamble and column heading.
        # The historical parser accepted ``ID`` and later asked AMPS to run an
        # unknown test called "id", terminating a monolithic --all command.
        listing = """Registered srcSEP standalone component tests
ID | group | class | initialization | build modes | description
PARK01 | parker | deterministic
FTED08 | fte-dmumu | extended
PARK01 | parker | duplicate
"""
        self.assertEqual(runner._parse_list_output(listing),
                         ["FTED08", "PARK01"])

    def test_isolated_all_continues_and_lists_failed_and_error_tests(self):
        """Retain GOOD01 after BAD02 errors and FAIL03 reports failure.

        A fake process layer keeps this unit test independent of AMPS while
        exercising the real per-ID paths, synthetic ERROR record, aggregate
        JSON/JUnit generation, printed results, and actionable final lists.
        """
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-isolated-all-") as tmp:
            temporary = Path(tmp)
            executable = temporary / "amps"
            executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            executable.chmod(0o755)
            output = temporary / "output"
            output.mkdir()
            arguments = SimpleNamespace(
                amps=str(executable), model_args=[], mpi_np=None,
                mpiexec="mpiexec", timeout=None)
            called: list[list[str]] = []

            def fake_run(command, cwd, log_path, environment=None,
                         timeout_s=None):
                del cwd, log_path, environment, timeout_s
                command = list(command)
                called.append(command)
                identifier = command[command.index("--test") + 1]
                if identifier == "BAD02":
                    # Model a fatal signal/abort: no JSON exists for this ID.
                    return 134
                result_status = "FAIL" if identifier == "FAIL03" else "PASS"
                result_message = ("analytical tolerance exceeded"
                                  if result_status == "FAIL" else
                                  "fixture passed")
                report_path = Path(command[command.index("--test-json") + 1])
                report_path.write_text(json.dumps({
                    "schema": "srcsep-component-tests-v1",
                    "exit_code": 1 if result_status == "FAIL" else 0,
                    "totals": {"passed": int(result_status == "PASS"),
                               "failed": int(result_status == "FAIL"),
                               "skipped": 0, "errors": 0},
                    "results": [{
                        "id": identifier, "status": result_status,
                        "message": result_message, "elapsed_seconds": 0.01,
                        "seed": 7, "configuration": [], "metrics": [],
                        "artifacts": [],
                    }],
                }), encoding="utf-8")
                return 1 if result_status == "FAIL" else 0

            stdout = io.StringIO()
            with (mock.patch.object(runner, "_discover_all_ids",
                                    return_value=["BAD02", "FAIL03", "GOOD01"]),
                  mock.patch.object(runner, "_validation_case_ids",
                                    return_value=set()),
                  mock.patch.object(runner, "_run_streaming",
                                    side_effect=fake_run),
                  mock.patch("sys.stdout", stdout)):
                exit_code, report_path, commands = runner._run_all_tests(
                    arguments, output, output / "test-run.log")
                # main() calls the final summary only after plot/manifest output;
                # invoke the same helper here because this focused unit calls
                # the lower-level isolated campaign function directly.
                runner._print_summary(json.loads(
                    report_path.read_text(encoding="utf-8")))

            self.assertEqual(exit_code, 2)
            self.assertEqual(len(called), 3)
            self.assertEqual(commands, called)
            self.assertEqual([command[command.index("--test") + 1]
                              for command in called],
                             ["BAD02", "FAIL03", "GOOD01"])
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["totals"], {
                "passed": 1, "failed": 1, "skipped": 0, "errors": 1})
            self.assertTrue((output / "srcsep-tests.xml").is_file())
            rendered = stdout.getvalue()
            self.assertIn("RESULT BAD02: ERROR", rendered)
            self.assertIn("RESULT FAIL03: FAIL", rendered)
            self.assertIn("RESULT GOOD01: PASS", rendered)
            self.assertIn(
                "Overall test summary: TOTAL=3 PASS=1 FAIL=1 SKIP=0 ERROR=1",
                rendered)
            self.assertIn("Failed tests (1):\n  - FAIL03: analytical tolerance exceeded",
                          rendered)
            self.assertIn("Error tests (1):\n  - BAD02: process exited 134",
                          rendered)

    def test_publication_plot_labels_name_source_and_exact_figures(self):
        """Ensure detached XM02/XM03 figures remain self-attributing."""
        runner = _load_cross_model_runner()
        for case_id, expected in (
                ("XM02", "Figure 7"),
                ("XM03", "Figure 12(a), Figure 12(b), Figure 12(c)")):
            case = json.loads((ROOT / "validation" / "cases" / case_id /
                               "input.json").read_text(encoding="utf-8"))
            label = runner._publication_plot_label(case)
            self.assertIn("Reference:", label)
            self.assertIn(case["reference"]["plot_citation"], label)
            self.assertIn(expected, label)

    def test_xm03_uses_one_global_amplitude_for_all_observations(self):
        """Prevent accidental per-time or per-instrument normalization.

        The synthetic linked spectrum is exactly one fifth of every reference
        value. A correct global log-space nuisance fit must recover five while
        preserving zero residual and full coverage across all three times.
        """
        runner = _load_cross_model_runner()
        case = json.loads((ROOT / "validation" / "cases" / "XM03" /
                           "input.json").read_text(encoding="utf-8"))
        model = []
        reference = []
        for elapsed, time_factor in ((4.0, 1.0), (12.0, 10.0), (36.0, 100.0)):
            for energy, spectral_factor in ((1.0, 2.0), (10.0, 0.2)):
                observed = time_factor * spectral_factor
                model.append({
                    "elapsed_hours": str(elapsed),
                    "energy_mev": str(energy),
                    "relative_differential_intensity": str(observed / 5.0),
                    "effective_sample_count": "1000",
                })
                reference.append({
                    "elapsed_hours": str(elapsed),
                    "instrument": "fixture",
                    "energy_low_mev": str(0.9 * energy),
                    "energy_high_mev": str(1.1 * energy),
                    "effective_energy_mev": str(energy),
                    "differential_intensity_pfu_per_mev": str(observed),
                })
        metrics, scale, comparison = runner._xm03_score(case, model, reference)
        self.assertAlmostEqual(scale, 5.0, places=12)
        by_name = {item["name"]: item["value"] for item in metrics}
        self.assertEqual(by_name["observation_point_coverage"], 1.0)
        self.assertAlmostEqual(by_name["global_log10_intensity_rmse"], 0.0)
        self.assertAlmostEqual(by_name["log10_intensity_correlation"], 1.0)
        for row in comparison:
            self.assertAlmostEqual(
                float(row["model_scaled_pfu_per_mev"]),
                float(row["differential_intensity_pfu_per_mev"]), places=12)

    def test_xm02_and_xm03_duration_keys_match_their_registered_inputs(self):
        """Exercise argument construction before a costly linked execution.

        XM02 retains the historical ``duration_hours`` field. XM03 names its
        06:00-UTC origin explicitly because its spectral clock begins at the
        later CME launch. Constructing both native vectors here detects a key
        rename applied to the wrong case before AMPS is started.
        """
        runner = _load_cross_model_runner()
        xm02_path = ROOT / "validation" / "cases" / "XM02" / "input.json"
        xm03_path = ROOT / "validation" / "cases" / "XM03" / "input.json"
        xm02 = json.loads(xm02_path.read_text(encoding="utf-8"))
        xm03 = json.loads(xm03_path.read_text(encoding="utf-8"))
        xm02_arguments = runner._xm02_arguments(xm02)
        xm03_arguments = runner._xm03_arguments(xm03, xm03_path)
        self.assertEqual(
            xm02_arguments[xm02_arguments.index("--duration-s") + 1],
            "158400.0")
        self.assertEqual(
            xm03_arguments[xm03_arguments.index("--duration-s") + 1],
            "158400.0")

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
            # The actionable lists intentionally close the transcript, after
            # plot generation and the output-directory announcement. A clean
            # report still names both categories and marks them empty.
            self.assertGreater(completed.stdout.rfind("Failed tests (0):"),
                               completed.stdout.rfind("Results:"))
            self.assertGreater(completed.stdout.rfind("Error tests (0):"),
                               completed.stdout.rfind("Failed tests (0):"))
            self.assertTrue(completed.stdout.rstrip().endswith(
                "Error tests (0):\n  none"), completed.stdout)

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

    def test_validation_mode_forwards_the_selected_linked_executable(self):
        """Prevent CV01 from regressing to a Python-compiled model driver."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-validation-command-") as tmp:
            output = Path(tmp)
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps",
                validation_all=False,
                validation_cases=["CV01"],
                case_input=None,
                timeout=42.0,
            )
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[:4], [
                sys.executable,
                str(ROOT / "validation" / "run_case.py"),
                "--amps",
                "/opt/amps/bin/srcsep-amps",
            ])
            self.assertIn("CV01", command)
            self.assertIn("--timeout", command)
            run.assert_called_once()

    def test_cv02_cv05_share_the_linked_application_selection_path(self):
        """Require every new controlled case to use the requested AMPS binary."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-validation-command-") as tmp:
            output = Path(tmp)
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps",
                validation_all=False,
                validation_cases=["CV02", "CV03", "CV04", "CV05"],
                case_input=None,
                timeout=120.0,
            )
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in ("CV02", "CV03", "CV04", "CV05"):
                self.assertIn(case_id, command)
            self.assertNotIn("controlled_transport_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_cv06_cv12_share_the_linked_application_selection_path(self):
        """Keep advanced physics cases attached to the selected application.

        The case-local source harness is useful for compile testing, but it is
        not scientific evidence. This assertion prevents future convenience
        changes from bypassing validation/run_case.py or replacing ``--amps``
        with advanced_validation_models.cpp.
        """
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-validation-command-") as tmp:
            output = Path(tmp)
            selected = [f"CV{index:02d}" for index in range(6, 13)]
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps",
                validation_all=False,
                validation_cases=selected,
                case_input=None,
                timeout=600.0,
            )
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in selected:
                self.assertIn(case_id, command)
            self.assertNotIn("advanced_validation_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_iv01_iv06_use_the_selected_linked_application(self):
        """Integrated cases must never be relabelled standalone evidence."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-iv-command-") as tmp:
            output = Path(tmp)
            selected = [f"IV{index:02d}" for index in range(1, 7)]
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps", validation_all=False,
                validation_cases=selected, case_input=None, timeout=600.0)
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in selected:
                self.assertIn(case_id, command)
            self.assertNotIn("integrated_validation_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_xm01_xm03_use_the_selected_linked_application(self):
        """Cross-model orchestration must retain the production binary gate."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-xm-command-") as tmp:
            output = Path(tmp)
            selected = ["XM01", "XM02", "XM03"]
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps", validation_all=False,
                validation_cases=selected, case_input=None, timeout=900.0)
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in selected:
                self.assertIn(case_id, command)
            # XM02/XM03 resolve their only input from the registry.  The
            # top-level runner must not manufacture an --input override while
            # assembling a multi-case command.
            self.assertNotIn("--input", command)
            self.assertNotIn("cross_model_validation_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_xm02_xm03_reject_command_line_input_overrides(self):
        """Keep each publication comparison tied to its registered input.

        This check occurs before executable validation, so a deliberately
        missing binary can be used to prove that input-policy diagnostics are
        deterministic and do not depend on the local AMPS installation.
        """
        with tempfile.TemporaryDirectory(prefix="srcsep-xm-fixed-input-") as tmp:
            temporary = Path(tmp)
            alternate = temporary / "alternate.json"
            alternate.write_text("{}\n", encoding="utf-8")
            for case_id in ("XM02", "XM03"):
                completed = subprocess.run(
                    [sys.executable, str(RUNNER), "--amps",
                     str(temporary / "missing-amps"), "--validation-case", case_id,
                     "--case-input", str(alternate), "--output-dir",
                     str(temporary / case_id)],
                    cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, check=False)
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn(
                    f"{case_id} uses its single registered publication-derived input",
                    completed.stdout)

    def test_direct_case_runner_enforces_fixed_xm_input(self):
        """Apply the same no-override rule below the convenience front end."""
        with tempfile.TemporaryDirectory(prefix="srcsep-xm-direct-input-") as tmp:
            temporary = Path(tmp)
            alternate = temporary / "alternate.json"
            alternate.write_text("{}\n", encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(CASE_RUNNER), "--amps",
                 str(temporary / "missing-amps"), "--case", "XM02", "--input",
                 str(alternate), "--output-dir", str(temporary / "output")],
                cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=False)
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn(
                "XM02 uses its single registered publication-derived input",
                completed.stdout)

    def test_validation_case_rejects_a_missing_linked_executable(self):
        """Never replace unavailable application evidence with a local build."""
        with tempfile.TemporaryDirectory(prefix="srcsep-missing-amps-") as tmp:
            missing = Path(tmp) / "amps-does-not-exist"
            completed = subprocess.run(
                [sys.executable, str(CASE_RUNNER), "--amps", str(missing),
                 "--case", "CV01", "--output-dir", str(Path(tmp) / "out")],
                cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=False)
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn("linked srcSEP/AMPS executable is missing", completed.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)

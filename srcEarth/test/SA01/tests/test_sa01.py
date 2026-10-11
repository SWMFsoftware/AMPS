#!/usr/bin/env python3
"""Numerical, schema, gate, preparation, and safety tests for SA01."""

from __future__ import annotations

import json
import math
import shlex
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

TEST_DIR = Path(__file__).resolve().parent
SA01_DIR = TEST_DIR.parent
AMPS_ROOT = SA01_DIR.parents[2]
sys.path.insert(0, str(SA01_DIR))

from prepare_SA01 import prepare_trajectory, render_deck  # noqa: E402
from run_SA01 import (  # noqa: E402
    parser as run_parser,
    record_workflow_completion,
    reserve_campaign_output,
    safe_rebuild,
)
from sa01_common import (  # noqa: E402
    CHANNELS, LIGHT_SPEED_CM_S, PROTON_REST_MEV, bin_average,
    fractional_change, isotropic_density_center_cm3, ratio_metrics,
    read_point_table, read_spectrum_zones, spectrum_integral,
)


class SA01Tests(unittest.TestCase):
    def test_unit_conversion_and_relativistic_density(self) -> None:
        values = [2.0] * 7
        expected = 0.0
        for value, (_, low, high) in zip(values, CHANNELS):
            center = 0.5 * (low + high)
            gamma = 1.0 + center / PROTON_REST_MEV
            speed = LIGHT_SPEED_CM_S * math.sqrt(1.0 - gamma ** -2)
            expected += 4.0 * math.pi * value * (high - low) / speed
        self.assertAlmostEqual(isotropic_density_center_cm3(values), expected, places=22)
        self.assertEqual(1.0e4 * 1.0e-4, 1.0)

    def test_power_law_channel_integral(self) -> None:
        # J=10000/E^2 integrates to 300 over [20,50], average 10.
        self.assertAlmostEqual(bin_average([10, 100], [100, 1], 20, 50), 10.0, places=11)

    def test_physical_zero_is_not_deleted_or_bridged(self) -> None:
        # Piecewise linear transitions around a zero give two triangles, not a
        # positive log interpolation spanning from E=1 directly to E=3.
        self.assertAlmostEqual(spectrum_integral([1, 2, 3], [1, 0, 1], 1, 3), 1.0)
        metrics = ratio_metrics([{"model": 0.0, "reference": 1.0},
                                 {"model": 1.0, "reference": 1.0}])
        self.assertEqual(metrics["count"], 2)
        self.assertEqual(metrics["zero_model_count"], 1)
        self.assertEqual(metrics["factor_two_fraction"], 0.5)

    def test_nonfinite_and_support_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            spectrum_integral([1, 2], [1, np.nan], 1, 2)
        with self.assertRaises(ValueError):
            spectrum_integral([1, 2], [1, 1], 0.5, 2)

    def test_duplicate_tecplot_times_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "spectrum.dat"
            path.write_text(
                'VARIABLES="E_MeV" "J_local_perMeV" "unresolved_fraction"\n'
                'ZONE T="T0 2014-01-06T09:00:00Z" I=2 F=POINT\n18.5 1 0\n116 1 0\n'
                'ZONE T="T1 2014-01-06T09:00:00Z" I=2 F=POINT\n18.5 1 0\n116 1 0\n',
                encoding="ascii",
            )
            with self.assertRaises(ValueError):
                read_spectrum_zones(path)

    def test_malformed_point_table_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "points.dat"
            path.write_text('VARIABLES="TimeUTC" "N_cm^-3"\nZONE I=1 F=POINT\n"2014-01-01T00:00:00Z"\n')
            with self.assertRaises(ValueError):
                read_point_table(path)

    def test_boundary_availability_selection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trajectory, boundary, output = root / "trajectory", root / "boundary", root / "out"
            trajectory.write_text(
                "2014-01-06T09:00:00Z 0 0 1000\n2014-01-06T09:15:00Z 0 0 1000\n"
                + "\n".join(f"2014-01-06T{9 + (30 + i * 15) // 60:02d}:{(30 + i * 15) % 60:02d}:00Z 0 0 1000" for i in range(12)) + "\n"
            )
            epochs = ["2014-01-06T09:00:00Z"] + [f"2014-01-06T{9 + (30 + i * 15) // 60:02d}:{(30 + i * 15) % 60:02d}:00Z" for i in range(12)]
            boundary.write_text("\n".join(epoch + " " + " ".join(["1"] * 10) for epoch in epochs) + "\n")
            result = prepare_trajectory(trajectory, boundary, output, 15)
            self.assertEqual(result["kept"], 13)
            self.assertEqual(result["excluded_by_reason"]["missing_exact_boundary_spectrum"], 1)

    def test_template_uses_absolute_paths_and_strict_gates(self) -> None:
        replacements = {
            "RUN_ID": "test", "FIRST_EPOCH": "2014-01-06T09:00:00",
            "DRIVER_FILE": "/data/driver", "SPECTRUM_FILE": "/data/spectrum",
            "TRAJECTORY_FILE": "/data/trajectory", "CUTOFF_NODES": 24,
            "DS_INTERVALS": 24, "SPECTRUM_NODES": 25,
            "TOTAL_PARTICLES": 2400, "THREADS": 2,
            "FIELD_EVAL_METHOD": "GRIDLESS",
            "MODE3D_MESH_BLOCK": "! unused",
            "TEMPORAL_BLOCK": "! unused",
            "PARALLEL_BLOCK": "  GRIDLESS_PARALLEL THREADS\n  GRIDLESS_THREADS 2",
        }
        deck = render_deck(replacements)
        self.assertIn("DRIVER_FILE /data/driver", deck)
        self.assertIn("DS_FAIL_ON_UNRESOLVED T", deck)
        self.assertIn("DS_TRANSMISSION_SCAN_N 25", deck)
        self.assertIn("DS_MAX_PARTICLES 2400", deck)
        self.assertIn("DS_UNRESOLVED_TOL 0.01", deck)
        self.assertIn("DS_UNRESOLVED_EXTENSION_PASSES 3", deck)
        self.assertIn("DS_UNRESOLVED_EXTENSION_FACTOR 2.0", deck)
        self.assertIn("TRAP_DRIFT_DETECTION T", deck)
        self.assertIn("TRAP_MIN_DRIFT_REVOLUTIONS 3", deck)
        self.assertIn("SPEC_TIME_MAX_GAP_S 900", deck)
        self.assertNotIn("{{", deck)

    def test_template_renders_explicit_gridded_snapshot_contract(self) -> None:
        replacements = {
            "RUN_ID": "test_gridded", "FIRST_EPOCH": "2014-01-06T09:00:00",
            "DRIVER_FILE": "/data/driver", "SPECTRUM_FILE": "/data/spectrum",
            "TRAJECTORY_FILE": "/data/trajectory", "CUTOFF_NODES": 24,
            "DS_INTERVALS": 24, "SPECTRUM_NODES": 25,
            "TOTAL_PARTICLES": 2400, "THREADS": 8,
            "FIELD_EVAL_METHOD": "GRID_3D",
            "MODE3D_MESH_BLOCK": (
                "#MODE3D_MESH\n  MODE3D_MESH_RES_EARTH_RE 0.05\n"
                "  MODE3D_MESH_RES_BOUNDARY_RE 2.0\n"
                "  MODE3D_MESH_COARSENING LINEAR\n  MODE3D_MESH_EXPONENT 1.0"
            ),
            "TEMPORAL_BLOCK": (
                "#TEMPORAL\n  TEMPORAL_MODE SNAPSHOT_LIST\n"
                "  SNAPSHOT_LIST_FILE /data/snapshot_epochs.txt"
            ),
            "PARALLEL_BLOCK": (
                "  DENSITY_PARALLEL OPENMP\n  DENSITY_THREADS 8\n"
                "  MODE3D_TRAJECTORY_SNAPSHOT_SCHEDULER EPOCH_PARALLEL\n"
                "  MODE3D_EPOCH_PARALLEL_MEMORY_FRACTION 0.80\n"
                "  MODE3D_MPI_SCHEDULER DYNAMIC\n  MODE3D_MPI_DYNAMIC_CHUNK 1"
            ),
        }
        deck = render_deck(replacements)
        self.assertIn("FIELD_EVAL_METHOD GRID_3D", deck)
        self.assertIn("TEMPORAL_MODE SNAPSHOT_LIST", deck)
        self.assertIn("MODE3D_MESH_RES_EARTH_RE 0.05", deck)
        self.assertIn("DENSITY_PARALLEL OPENMP", deck)
        self.assertIn(
            "MODE3D_TRAJECTORY_SNAPSHOT_SCHEDULER EPOCH_PARALLEL", deck
        )
        self.assertIn("DS_FAIL_ON_UNRESOLVED T", deck)
        self.assertIn("DS_UNRESOLVED_EXTENSION_PASSES 3", deck)
        self.assertIn("TRAP_DRIFT_DETECTION T", deck)
        self.assertNotIn("{{", deck)

    def test_gate_denominator_and_near_zero_convergence(self) -> None:
        rows = [{"model": value, "reference": 1.0} for value in (1.0, 1.5, 0.0, 3.0)]
        metric = ratio_metrics(rows)
        self.assertEqual(metric["count"], 4)
        self.assertEqual(metric["factor_two_fraction"], 0.5)
        self.assertEqual(fractional_change(0.0, 0.0), 0.0)
        self.assertGreater(fractional_change(0.0, 1e-20), 0.02)

    def test_rebuild_refuses_unexpected_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                safe_rebuild(Path(temporary))

    def test_list_driven_output_reservation_is_unique_and_non_destructive(self) -> None:
        """Repeated list invocations must never reuse prior physical evidence."""

        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary) / "campaigns"
            first = reserve_campaign_output(parent, "QUICK")
            second = reserve_campaign_output(parent, "QUICK")
            self.assertNotEqual(first, second)
            self.assertEqual(list(first.iterdir()), [])
            self.assertEqual(list(second.iterdir()), [])
            self.assertTrue(first.name.startswith("sa01_quick_"))

    def test_workflow_exit_contract_preserves_scientific_failure(self) -> None:
        """Code-validation PASS must not rewrite a failed observation result."""

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            physical = output / "physical_execution_summary.json"
            comparison_path = output / "comparison_summary.json"
            plot_manifest = output / "plot_manifest.json"
            qualification_path = output / "qualification_summary.json"
            report = output / "SA01_validation_report.md"
            figure = output / "figures" / "comparison.png"
            figure.parent.mkdir()
            physical.write_text('{"status":"PASS"}\n', encoding="utf-8")
            comparison_path.write_text('{"status":"FAIL"}\n', encoding="utf-8")
            plot_manifest.write_text('{"figure_sha256":{}}\n', encoding="utf-8")
            qualification_path.write_text('{"status":"FAIL"}\n', encoding="utf-8")
            report.write_text("# failed scientific qualification\n", encoding="utf-8")
            figure.write_bytes(b"not-an-image-but-a-stable-test-artifact")

            result = record_workflow_completion(
                output, {"status": "FAIL"}, [str(figure)], report,
                qualification="FAIL", profile="QUICK",
            )
            saved = json.loads((output / "workflow_completion.json").read_text())
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(saved["scientific_qualification"], "FAIL")
            self.assertEqual(saved["observational_status"], "FAIL")
            self.assertTrue(saved["qualification_is_unchanged"])

    def test_quick_and_standard_are_independent_repository_commands(self) -> None:
        """Freeze the requested profile ownership and complete-workflow options."""

        test_list = (AMPS_ROOT / "srcEarth/test/list").read_text(encoding="utf-8")
        studies_list = (AMPS_ROOT / "srcEarth/studies/list").read_text(encoding="utf-8")
        quick = [line for line in test_list.splitlines()
                 if line.startswith("P ") and "test/SA01/run_SA01.py" in line]
        standard = [line for line in studies_list.splitlines()
                    if "test/SA01/run_SA01.py" in line and not line.startswith("!")]
        self.assertEqual(len(quick), 1)
        self.assertEqual(len(standard), 1)
        for command, profile, parent in (
                (quick[0], "QUICK", "quick_test_runs"),
                (standard[0], "STANDARD", "standard_study_runs")):
            self.assertIn(f"--profile {profile}", command)
            self.assertIn(f"--output-parent test_output/SA01/{parent}", command)
            self.assertIn("--formats png eps", command)
            self.assertIn("--exit-code-policy WORKFLOW", command)
            self.assertIn("--solver GRIDDED", command)
            self.assertIn("--mover HC4", command)
        # Exclusivity is a test-runner scheduling property, not a private CLI
        # prerequisite. Keeping the exact QUICK command standalone-executable lets
        # operators reproduce the registered workflow without fabricating runner
        # environment state.
        self.assertNotIn("--require-runner-exclusive", quick[0])
        self.assertNotIn("--require-runner-exclusive", standard[0])
        quick_index = test_list.splitlines().index(quick[0])
        self.assertEqual(test_list.splitlines()[quick_index - 1], "! runner: exclusive")

        # The parser itself must accept a reusable output parent without weakening
        # the default qualification-based exit behavior for manual invocations.
        parsed = run_parser().parse_args([
            "--data-root", "/data/package", "--output-parent", "/tmp/runs",
            "--profile", "QUICK",
        ])
        self.assertEqual(parsed.profile, "QUICK")
        self.assertEqual(parsed.exit_code_policy, "QUALIFICATION")

    def test_registered_quick_command_is_directly_parseable(self) -> None:
        """The exact master-list payload must also be a valid standalone CLI."""

        test_list = (AMPS_ROOT / "srcEarth/test/list").read_text(encoding="utf-8")
        command = next(line for line in test_list.splitlines()
                       if line.startswith("P ") and "test/SA01/run_SA01.py" in line)
        tokens = shlex.split(command)
        # Remove only the repository expectation marker, interpreter, and script.
        # Parsing the untouched remainder prevents a future list-only prerequisite
        # from recreating the direct-execution failure reported by the operator.
        self.assertEqual(tokens[:3], ["P", "python3", "srcEarth/test/SA01/run_SA01.py"])
        parsed = run_parser().parse_args(tokens[3:])
        self.assertEqual(parsed.profile, "QUICK")
        self.assertEqual(parsed.exit_code_policy, "WORKFLOW")
        self.assertEqual(parsed.output_parent, Path("test_output/SA01/quick_test_runs"))

    def test_mode3d_scheduler_preserves_epoch_ownership_boundary(self) -> None:
        """Guard the communicator split that makes concurrent T05 epochs safe.

        This is deliberately a source-contract test in addition to linked runtime
        tests: replacing either singleton physics collectives or the world epoch queue
        with one communicator can compile cleanly but mix fields or deadlock only in a
        multi-epoch production run.
        """
        mode3d = (AMPS_ROOT / "srcEarth/3d/Mode3D.cpp").read_text(encoding="utf-8")
        self.assertIn('snapshotScheduler=="EPOCH_PARALLEL"', mode3d)
        self.assertIn(
            'outputMode=="TRAJECTORY" && Mode3DSnapshotListRequested(prm)',
            mode3d,
        )
        self.assertIn("PIC::InitMPI(/*independentDomainMode=*/trajectoryEpochParallel)", mode3d)
        self.assertIn("MPI_COMM_WORLD,static_cast<long long>(snapshotEpochs.size()),1", mode3d)
        self.assertIn("MPI_Barrier(MPI_COMM_WORLD);", mode3d)
        self.assertIn("epochScheduler.reset();", mode3d)

    def test_mode3d_density_uses_fine_tasks_only_for_shared_epoch_modes(self) -> None:
        """Ensure POINTS/SHELLS and TRAJECTORY cannot silently swap work units."""
        density = (AMPS_ROOT / "srcEarth/3d/DensityMode3D.cpp").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            '(outputMode == "POINTS" || outputMode == "SHELLS") && mpiSize > 1',
            density,
        )
        self.assertIn(
            "const long long nEnergyTasks=static_cast<long long>(nLoc)*",
            density,
        )
        self.assertIn("computeOneSharedEpochEnergy", density)
        self.assertIn("if (sharedEpochFineTasks)", density)
        self.assertIn("EvaluateProductSet_(prm,E_MeV,T,TLower,TUpper)", density)


if __name__ == "__main__":
    unittest.main(verbosity=2)

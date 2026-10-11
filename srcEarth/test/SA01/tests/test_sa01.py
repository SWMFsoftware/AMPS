#!/usr/bin/env python3
"""Numerical, schema, gate, preparation, and safety tests for SA01."""

from __future__ import annotations

import math
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
from run_SA01 import safe_rebuild  # noqa: E402
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

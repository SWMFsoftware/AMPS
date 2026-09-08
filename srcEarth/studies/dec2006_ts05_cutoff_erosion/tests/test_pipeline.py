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
from make_figures import cutoff_degradation_figures
from run_morphology import load_c10_module, snapshot_suffix, split_multishell_access
from run_study import execute, independent_validation_remains
import run_study
from study_common import read_driver


class PipelineTests(unittest.TestCase):
    def test_snapshot_suffix_matches_native_mode3d_filename_contract(self):
        """The postprocessor must address one exact file per batch epoch."""

        from datetime import datetime, timezone
        epoch = datetime(2006, 12, 15, 0, 50, tzinfo=timezone.utc)
        self.assertEqual(
            snapshot_suffix(2, epoch),
            "_snapshot_000002_2006_12_15T00_50_00",
        )

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

    def test_shared_morphology_then_pamela_failure_still_allows_poes(self):
        """Shared production runs first; both independent comparisons are collected."""

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
            self.assertEqual(manifest["return_codes"]["morphology"], 0)
            self.assertEqual(manifest["return_codes"]["pamela"], 1)
            self.assertEqual(manifest["return_codes"]["poes"], 0)

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

    def test_smoke_prepare_renders_one_native_multi_epoch_batch(self):
        """BATCHED must mean one process/mesh, not directory-only grouping."""

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "morphology"
            completed = subprocess.run([
                sys.executable, str(ROOT / "scripts" / "run_morphology.py"),
                "--profile", "SMOKE", "--prepare-only", "--output-root", str(output),
            ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            result = json.loads((output / "morphology_result.json").read_text())
            self.assertEqual(result["n_cases"], 8)
            self.assertEqual(result["n_amps_launches"], 1)
            self.assertEqual(result["mesh_layout"], "BATCHED")
            self.assertTrue(result["mesh_reused_across_shells"])
            self.assertTrue(result["mesh_reused_across_epochs"])
            inputs = list(output.rglob("AMPS_PARAM_C10.in"))
            self.assertEqual(len(inputs), 1)
            for path in inputs:
                text = path.read_text()
                self.assertIn("CUTOFF_SAMPLING         VERTICAL", text)
                self.assertIn("CUTOFF_SEARCH_ALGORITHM RIGIDITY_LIST", text)
                self.assertIn("TEMPORAL_MODE          SNAPSHOT_LIST", text)
                self.assertIn("SNAPSHOT_LIST_FILE     snapshot_epochs.txt", text)
                self.assertIn("SHELL_COUNT 2", text)
                self.assertIn("SHELL_ALTS_KM 475 850", text)
                self.assertNotIn("CUTOFF_UNRESOLVED_EXTENSION_PASSES", text)
                epochs = (path.parent / "snapshot_epochs.txt").read_text()
                self.assertEqual(len([
                    line for line in epochs.splitlines()
                    if line and not line.startswith("#")
                ]), 4)
            inventory = json.loads((output / "command_inventory.json").read_text())
            self.assertEqual(len(inventory), 1)
            for item in inventory:
                self.assertNotIn("--epoch", item["command"])
                self.assertTrue(item["reuses_mesh_across_epochs"])
                self.assertEqual(item["command"][item["command"].index("-mover") + 1],
                                 "RK4")

    def test_native_mode3d_supports_snapshot_list_shells_with_one_mesh(self):
        """Guard new SHELLS support and the pre-existing dispatch boundaries."""

        source_path = ROOT.parents[1] / "3d" / "Mode3D.cpp"
        source = source_path.read_text(encoding="utf-8")
        self.assertGreaterEqual(source.count('if (outputMode=="SHELLS") return snap;'), 1)
        self.assertGreaterEqual(source.count('if (outputMode=="SHELLS") return;'), 1)

        # The earliest return is the compatibility firewall: every ordinary
        # SNAPSHOT/TIME_SERIES calculation exits before the new SHELLS branch.
        # Within SNAPSHOT_LIST, TRAJECTORY must retain its historical timestamp
        # filtering, flattened-point rebuild, and aperture-index remapping.
        build_start = source.index("EarthUtil::AmpsParam Mode3DBuildSnapshotWorkParam(")
        build_end = source.index("void Mode3DValidateSnapshotListCoverage(", build_start)
        build = source[build_start:build_end]
        self.assertLess(
            build.index("if (!Mode3DSnapshotListRequested(snap)) return snap;"),
            build.index('if (outputMode=="SHELLS") return snap;'),
        )
        for legacy_marker in (
            "globalToLocal", "RebuildFlattenedPointsFromTrajectories",
            "remapped.locationIndex=found->second",
        ):
            self.assertIn(legacy_marker, build)

        validate_start = build_end
        validate_end = source.index("EarthUtil::AmpsParam Mode3DBuildSnapshotParam(",
                                    validate_start)
        validation = source[validate_start:validate_end]
        self.assertLess(
            validation.index("if (!Mode3DSnapshotListRequested(prm)) return;"),
            validation.index('if (outputMode=="SHELLS") return;'),
        )
        self.assertIn("every location must match exactly one epoch", validation)

        mesh_position = source.index("  amps_init_mesh();   // build")
        snapshot_loop_position = source.index(
            "for (std::size_t iSnapshot=0; iSnapshot<snapshotEpochs.size();"
        )
        self.assertLess(mesh_position, snapshot_loop_position)

    def test_multishell_split_requires_altitude_labeled_zones(self):
        """The optimizer must never infer shell identity from row order."""

        c10 = load_c10_module(ROOT)
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = work / "combined.dat"
            source.write_text(
                'TITLE="fixture"\n'
                'VARIABLES="lon_deg" "lat_deg" "rigidity_gv" '
                '"access_state" "allowed" "unresolved"\n'
                'ZONE T="Alt_km=475"\n0 40 0.5 1 1 0\n'
                'ZONE T="Alt_km=850"\n0 50 0.5 0 0 0\n'
            )
            destinations = {475.0: work / "475.dat", 850.0: work / "850.dat"}
            counts = split_multishell_access(
                source, (475.0, 850.0), destinations, c10
            )
            self.assertEqual(counts, {475.0: 1, 850.0: 1})
            self.assertEqual(
                c10.parse_tecplot_shell_access(destinations[475.0])[0].access_state,
                1,
            )

    def test_multishell_split_accepts_explicit_zero_based_shell_indices(self):
        """AMPS Shell_0/Shell_1 zone labels follow configured altitude order."""

        c10 = load_c10_module(ROOT)
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = work / "combined.dat"
            source.write_text(
                'TITLE="fixture"\n'
                'VARIABLES="lon_deg" "lat_deg" "rigidity_gv" '
                '"access_state" "allowed" "unresolved"\n'
                'ZONE T="Shell_0"\n0 40 0.5 1 1 0\n'
                'ZONE T="Shell_1"\n0 50 0.5 0 0 0\n'
            )
            destinations = {475.0: work / "475.dat", 850.0: work / "850.dat"}
            counts = split_multishell_access(
                source, (475.0, 850.0), destinations, c10
            )
            self.assertEqual(counts, {475.0: 1, 850.0: 1})

    def test_publication_degradation_figures_include_png_and_eps(self):
        """The required science visualization must include raster and vector files."""

        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = work / "cutoff_dynamics_timeseries.csv"
            with source.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=(
                    "epoch_utc", "altitude_km", "rigidity_gv", "hemisphere",
                    "cutoff_erosion_deg",
                ))
                writer.writeheader()
                for altitude in (475.0, 850.0):
                    for rigidity in (0.4, 0.7):
                        for hemisphere in ("N", "S"):
                            for epoch, erosion in (
                                ("2006-12-14T00:00:00Z", 0.0),
                                ("2006-12-15T00:00:00Z", -2.0),
                            ):
                                writer.writerow({
                                    "epoch_utc": epoch, "altitude_km": altitude,
                                    "rigidity_gv": rigidity,
                                    "hemisphere": hemisphere,
                                    "cutoff_erosion_deg": erosion,
                                })
            products = cutoff_degradation_figures(source, work / "figures")
            self.assertEqual(len(products), 6)
            for stem in ("figure_cutoff_degradation",
                         "figure_peak_cutoff_degradation"):
                self.assertTrue((work / "figures" / f"{stem}.png").is_file())
                self.assertTrue((work / "figures" / f"{stem}.eps").is_file())

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

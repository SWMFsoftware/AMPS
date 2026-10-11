"""Regression checks for saved-product coverage, refinement and plot compatibility."""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_SA01 import compare_branch, compare_convergence, validate_branch_evidence
from plot_SA01 import _representative_spectra, _time_energy_maps, _time_series_channels, plot_all
from prepare_SA01 import prepare_branch
from run_SA01 import write_report, qualification_status, diagnostic_evidence, wait_for_evidence
from sa01_common import CHANNELS, atomic_json, sha256_file, read_point_table, read_spectrum_zones, write_csv


def plot_fixture(root: Path) -> Path:
    """Manufactured plotting inputs: never represent these as physical evidence."""
    output = root / "synthetic_output"
    processed = root / "validation_data/processed/SA01"
    epochs = [f"2014-01-06T{hour:02d}:00:00Z" for hour in (9, 10, 14)]
    rows = [{"spacecraft": "RBSPA", "utc": utc, "channel": name,
             "energy_min_mev": low, "energy_max_mev": high, "model": 0 if i == 0 else 0.5,
             "reference": 1.0, "l_shell": 4.0, "mlt_hour": 6.0,
             "max_unresolved_fraction_at_epoch": 0.0}
            for utc in epochs for i, (name, low, high) in enumerate(CHANNELS)]
    for product in ("spectrum", "channel_flux", "total_flux", "density"):
        selected = rows if product in ("spectrum", "channel_flux") else rows[::7]
        write_csv(output / "RBSPA/comparison" / (product + "_comparison.csv"), selected)
    write_csv(output / "RBSPA/comparison/native_reconstruction.csv", [
        {"utc": utc, "total_flux_fractional_residual": 0.01, "density_fractional_residual": 0.005}
        for utc in epochs])
    write_csv(processed / "soho_erne_5min.csv", [
        {"utc": utc, "energy_center_mev": 57, "flux_cm2_s_sr_mev": 1} for utc in epochs])
    write_csv(processed / "qin_denton_5min.csv", [{"utc": utc, "dst_nt": -20} for utc in epochs])
    write_csv(output / "convergence_summary.csv", [
        {"spacecraft": "RBSPA", "utc": utc, "product": product, "fractional_change": 0.03}
        for utc in epochs for product in ("channel_flux", "total_flux", "density")])
    atomic_json(output / "comparison_summary.json", {"status": "SYNTHETIC FIXTURE", "pooled": {}})
    return output


def fixture(root: Path, nodes: int, last_channel: float = 1.0) -> Path:
    """Create a hash-bound twelve-epoch fixture with independent channel values."""
    branch = root / "RBSPA"
    branch.mkdir(parents=True)
    epochs = [f"2014-01-06T{9 + i:02d}:00:00Z" for i in range(12)]
    trajectory = branch / "trajectory.txt"
    trajectory.write_text("\n".join(t + " 0 0 1000" for t in epochs) + "\n")
    prepared = {"solver": "GRIDLESS", "trajectory": {
        "output": str(trajectory), "kept": 12, "sha256": sha256_file(trajectory)},
        "profile_parameters": {"spectrum_nodes": nodes, "requested_directions_per_energy": 96}}
    for key in ("input", "driver", "spectrum"):
        path = branch / (key + ".txt")
        path.write_text("immutable fixture\n")
        prepared[key + "_file"] = str(path)
        prepared[key + "_sha256"] = sha256_file(path)
    atomic_json(branch / "prepared_run.json", prepared)
    flux = branch / "gridless_points_flux.dat"
    flux.write_text('VARIABLES="TimeUTC" "F_tot_m2s1" ' +
                    " ".join(f'"F_{name}_m2s1"' for name, _, _ in CHANNELS) +
                    '\nZONE I=12 F=POINT\n' +
                    "\n".join(f'"{t}" 7 1 1 1 1 1 1 {last_channel}' for t in epochs) + "\n")
    density = branch / "gridless_points_density.dat"
    density.write_text('VARIABLES="TimeUTC" "N_cm^-3"\nZONE I=12 F=POINT\n' +
                       "\n".join(f'"{t}" 1' for t in epochs) + "\n")
    spectrum = branch / "gridless_points_spectrum.dat"
    spectrum.write_text('VARIABLES="E_MeV" "J_local_perMeV" "unresolved_fraction"\n' +
                        "".join(f'ZONE T="point {t}" I=2 F=POINT\n18.5 1 0\n116 1 0\n' for t in epochs))
    termination = branch / "gridless_termination_summary.dat"
    termination.write_text("fixture termination\n")
    manifest = branch / "standalone_run_manifest.json"
    manifest.write_text("{}\n")
    atomic_json(branch / "run_status.json", {**prepared, "status": "PASS", "return_code": 0,
        "output_sha256": {p.name: sha256_file(p) for p in (flux, density, spectrum, termination, manifest)}})
    return branch


class CompletionTests(unittest.TestCase):
    def test_wait_requires_terminal_evidence_not_process_disappearance(self):
        """An exited/reused runner cannot manufacture successful campaign output."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atomic_json(root / "active_run.json", {"runner_pid": 123})
            with patch("run_SA01.process_identity", side_effect=["start-one", None]):
                with self.assertRaisesRegex(RuntimeError, "without terminal evidence"):
                    wait_for_evidence(root)
            # An actual failure is a completed prerequisite, but stays a failure;
            # independent next runs may proceed without crediting it as a PASS.
            atomic_json(root / "physical_execution_summary.json", {"status": "ERROR"})
            self.assertEqual(wait_for_evidence(root)["status"], "ERROR")
            atomic_json(root / "diagnostic_status.json", {"status": "COMPLETE"})
            self.assertEqual(wait_for_evidence(root, diagnostic=True)["status"], "COMPLETE")

    def test_wait_observes_new_atomic_completion(self):
        """The wait loop notices completion without killing or restarting a job."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atomic_json(root / "active_run.json", {"runner_pid": 123})
            def finish_on_sleep(seconds):
                """Replace real waiting with a deterministic fake runner finish."""
                self.assertEqual(seconds, 5.0)
                atomic_json(root / "physical_execution_summary.json", {"status": "PASS"})
            with patch("run_SA01.process_identity", return_value="same-start"), \
                    patch("run_SA01.time.sleep", side_effect=finish_on_sleep):
                self.assertEqual(wait_for_evidence(root)["status"], "PASS")

    def test_qualification_never_promotes_missing_or_partial_evidence(self):
        """Successful execution/observation scores alone must not qualify SA01."""
        observation = {"status": "PASS", "spacecraft": {"RBSPA": {}, "RBSPB": {}}}
        self.assertEqual(qualification_status(observation, None), "INCOMPLETE")
        convergence = {"status": "PASS", "evidence": {"RBSPA": {}}}
        self.assertEqual(qualification_status(observation, convergence), "INCOMPLETE")
        convergence["evidence"]["RBSPB"] = {}
        self.assertEqual(qualification_status(observation, convergence), "PASS")
        convergence["status"] = "ERROR"
        self.assertEqual(qualification_status(observation, convergence), "INCOMPLETE")
        convergence["status"] = "FAIL"
        self.assertEqual(qualification_status(observation, convergence), "FAIL")

    def test_bounded_failure_cannot_be_hidden_by_other_refinement(self):
        """A sampling PASS cannot override a registered timestep failure or gap."""
        observation = {"status": "PASS", "spacecraft": {"RBSPA": {}, "RBSPB": {}}}
        convergence = {"status": "PASS", "evidence": {"RBSPA": {}, "RBSPB": {}}}
        self.assertEqual(qualification_status(observation, convergence, [{"status": "FAIL"}]), "FAIL")
        self.assertEqual(qualification_status(observation, convergence, [{"status": "INCOMPLETE"}]), "INCOMPLETE")
        self.assertEqual(qualification_status(observation, None, [{"status": "PASS"}]), "INCOMPLETE")

    def test_registered_diagnostics_persist_and_reject_changed_evidence(self):
        """Registration survives omitted CLI flags and detects result/native drift.

        Mock only native-run validation here (covered separately by diagnostics
        tests); exercise real on-disk comparison and campaign fingerprint checks.
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atomic_json(root / "RBSPA/prepared_run.json", {"fixture": True})
            comparison = root / "bounded.json"
            payload = {"check": "bounded_numerical_difference", "status": "FAIL",
                       "coarse": str(root / "a"), "fine": str(root / "b"),
                       "numerical_controls": {"coarse": {"control": 1}, "fine": {"control": 1}},
                       "artifact_sha256": {}, "maximum_fractional_change": .1, "threshold": .02}
            atomic_json(comparison, payload)
            with patch("analyze_diagnostics.completed", return_value={"control": 1}):
                registered = diagnostic_evidence(root, [comparison])
                self.assertEqual(diagnostic_evidence(root), registered)
                payload["status"] = "PASS"
                atomic_json(comparison, payload)
                with self.assertRaisesRegex(ValueError, "evidence changed"):
                    diagnostic_evidence(root)
                payload["status"] = "FAIL"
                atomic_json(comparison, payload)
            with patch("analyze_diagnostics.completed", return_value={"control": 2}):
                with self.assertRaisesRegex(ValueError, "stale coarse"):
                    diagnostic_evidence(root)
            atomic_json(root / "RBSPA/prepared_run.json", {"fixture": "changed"})
            with self.assertRaisesRegex(ValueError, "different prepared campaign"):
                diagnostic_evidence(root)

    def test_report_rejects_stale_convergence_evidence(self):
        """A saved PASS cannot be reused after either underlying run record changes."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture(root / "coarse", 2)
            fine = fixture(root / "fine", 3)
            compare_convergence(root / "coarse", root / "fine", ("RBSPA",),
                                root / "coarse/convergence_summary")
            comparison = {"status": "PASS", "spacecraft": {}, "pooled": {}, "gates": []}
            write_report(root / "coarse", comparison, [])
            status_path = fine / "run_status.json"
            status = json.loads(status_path.read_text())
            status["changed_after_comparison"] = True
            atomic_json(status_path, status)
            with self.assertRaisesRegex(ValueError, "stale fine convergence"):
                write_report(root / "coarse", comparison, [])

    def test_complete_comparison_uses_native_flux_and_checks_units(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            branch = fixture(root / "output", 2)
            path = branch / "gridless_points_spectrum.dat"
            text = path.read_text().replace('"unresolved_fraction"', '"unresolved_fraction" "N_sampled" "N_resolved"')
            text = text.replace("18.5 1 0", "18.5 1 0 96 96").replace("116 1 0", "116 1 0 96 96")
            path.write_text('AUXDATA SPECTRUM_INTENSITY_UNIT="m^-2 s^-1 sr^-1 MeV^-1"\n'
                            'AUXDATA SPECTRUM_ENERGY_BASIS="PER_PARTICLE"\n'
                            'AUXDATA STEP6_CHARACTERISTIC_MAPPING="STATIC_MAGNETIC"\n' + text)
            status_path = branch / "run_status.json"
            status = json.loads(status_path.read_text())
            status["output_sha256"][path.name] = sha256_file(path)
            atomic_json(status_path, status)
            epochs = [f"2014-01-06T{9 + i:02d}:00:00Z" for i in range(12)]
            processed = root / "data/processed/SA01"
            write_csv(processed / "rbspa_rept_spectrum_5min.csv", [
                dict(utc=utc, energy_min_mev=low, energy_max_mev=high,
                     differential_flux_cm2_s_sr_mev=1e-4, l_shell=4, mlt_hour=6)
                for utc in epochs for _, low, high in CHANNELS])
            write_csv(processed / "rbspa_rept_density_5min.csv", [
                dict(utc=utc, isotropic_omnidirectional_flux_cm2_s=4*math.pi*97.5e-4,
                     isotropic_equivalent_density_cm3=1) for utc in epochs])
            metrics, rows = compare_branch(root / "data", root / "output", "RBSPA")
            self.assertEqual(metrics["products"]["spectrum"]["factor_two_fraction"], 1)
            self.assertEqual(metrics["products"]["channel_flux"]["factor_two_fraction"], 0)
            self.assertEqual(rows["channel_flux"][0]["model"], 1e-4)
            path.write_text(path.read_text().replace('"m^-2 s^-1 sr^-1 MeV^-1"', '"cm^-2 s^-1 sr^-1 MeV^-1"'))
            status["output_sha256"][path.name] = sha256_file(path)
            atomic_json(status_path, status)
            with self.assertRaisesRegex(ValueError, "declared SPECTRUM_INTENSITY_UNIT"):
                compare_branch(root / "data", root / "output", "RBSPA")

    def test_all_plot_families_including_convergence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = plot_fixture(root)
            paths = plot_all(root, output, ("RBSPA",), ("png",))
            self.assertEqual(len(paths), 11)
            self.assertTrue(all(Path(path).stat().st_size > 1000 for path in paths))
            # Figures remain traceable when generated independently of execution.
            manifest = json.loads((output / "plot_manifest.json").read_text())
            self.assertEqual(len(manifest["figure_sha256"]), 11)
            self.assertEqual(set(manifest["source_sha256"]), {"plot_SA01.py", "sa01_common.py"})

    def test_mode3d_epoch_metadata_and_ambiguous_locations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = 'AUXDATA SNAPSHOT_EPOCH_UTC="2014-01-06T10:00:00Z"\n'
            point = root / "flux.dat"
            point.write_text(metadata + 'VARIABLES="F_tot_m2s1"\nZONE I=1 F=POINT\n1\n')
            self.assertEqual(read_point_table(point)[0]["TimeUTC"], "2014-01-06T10:00:00Z")
            point.write_text(point.read_text().replace("I=1", "I=2") + "2\n")
            with self.assertRaisesRegex(ValueError, "missing or duplicate"):
                read_point_table(point)
            spectrum = root / "spectrum.dat"
            spectrum.write_text(metadata + 'VARIABLES="E_MeV" "J_local_perMeV"\n'
                'ZONE T="loc_000000" I=2 F=POINT\n18.5 1\n116 1\n')
            self.assertEqual(read_spectrum_zones(spectrum)[0].utc, "2014-01-06T10:00:00Z")

    def test_publication_directions_and_explicit_mover(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            processed = root / "data/processed/SA01"
            processed.mkdir(parents=True)
            epochs = [f"2014-01-06T{9 + i:02d}:00:00Z" for i in range(12)]
            (processed / "rbspa_trajectory_geo_5min.txt").write_text(
                "\n".join(t + " 0 0 1000" for t in epochs) + "\n")
            (processed / "boundary_soho_erne_5min.dat").write_text(
                "\n".join(t + " 1 1" for t in epochs) + "\n")
            (processed / "ts05_driver_5min.txt").write_text("fixture\n")
            prepared = prepare_branch(root / "data", root / "output", "RBSPA", "PUBLICATION", 8,
                "GRIDDED", numerical_controls={"mover": "BORIS", "trace_dt": 0.001,
                                               "cadence_minutes": 60})
            deck = Path(prepared["input_file"]).read_text()
            self.assertIn("DS_DIRECTION_NMU 32", deck)
            self.assertIn("DS_DIRECTION_NPHI 48", deck)
            self.assertEqual(prepared["profile_parameters"]["requested_directions_per_energy"], 1536)
            self.assertEqual(prepared["numerical_controls"]["mover"], "BORIS")
            self.assertEqual(prepared["trajectory"]["kept"], 12)

    def test_refinement_catches_failed_seventh_channel(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture(root / "coarse", 25)
            fixture(root / "fine", 49, last_channel=1.1)
            result = compare_convergence(root / "coarse", root / "fine", ("RBSPA",), root / "convergence")
            self.assertEqual(result["status"], "FAIL")
            self.assertEqual(result["common_population_rows"], 12 * 9)

    def test_cadence_only_is_not_numerical_refinement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture(root / "coarse", 25)
            fixture(root / "fine", 25)
            with self.assertRaisesRegex(ValueError, "cadence alone"):
                compare_convergence(root / "coarse", root / "fine", ("RBSPA",), root / "convergence")

    def test_changed_product_cannot_reuse_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            branch = fixture(Path(directory), 25)
            self.assertEqual(len(validate_branch_evidence(branch)), 12)
            (branch / "gridless_points_density.dat").write_text("truncated\n")
            with self.assertRaisesRegex(ValueError, "changed output"):
                validate_branch_evidence(branch)

    def test_common_missing_epoch_does_not_hide_coverage_loss(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            branch = fixture(root, 25)
            # Every physical product loses the same epoch; update output hashes as
            # if the producer accepted the incomplete set. Ephemeris binding must
            # still reject it before observation lookup.
            for name in ("flux", "density", "spectrum"):
                path = branch / f"gridless_points_{name}.dat"
                text = path.read_text()
                if name == "spectrum":
                    text = text.replace('ZONE T="point 2014-01-06T09:00:00Z" I=2 F=POINT\n18.5 1 0\n116 1 0\n', "")
                else:
                    text = "\n".join(line for line in text.splitlines() if "2014-01-06T09:00:00Z" not in line) + "\n"
                    text = text.replace("I=12", "I=11")
                path.write_text(text)
            status_path = branch / "run_status.json"
            status = json.loads(status_path.read_text())
            status["output_sha256"] = {name: sha256_file(branch / name) for name in status["output_sha256"]}
            atomic_json(status_path, status)
            with self.assertRaisesRegex(ValueError, "prepared epochs"):
                compare_branch(root / "unused_reference", root, "RBSPA")

    def test_single_spacecraft_plot_families_on_supported_matplotlib(self):
        rows = [{"spacecraft": "RBSPA", "utc": f"2014-01-06T{hour:02d}:00:00Z",
                 "channel": name, "energy_min_mev": low, "energy_max_mev": high,
                 "reference": 1.0, "model": 0.0 if index == 0 else 0.5}
                for hour in (9, 10, 14) for index, (name, low, high) in enumerate(CHANNELS)]
        frame = pd.DataFrame(rows)
        frame["time"] = pd.to_datetime(frame.utc, utc=True)
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for render in (_representative_spectra, _time_energy_maps, _time_series_channels):
                paths += render(frame, Path(directory), ("png",))
            self.assertEqual(len(paths), 3)
            self.assertTrue(all(Path(path).stat().st_size > 1000 for path in paths))


if __name__ == "__main__":
    unittest.main()

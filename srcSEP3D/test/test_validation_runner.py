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
NATIVE_MATRIX = ROOT / "validation" / "run_native_matrix.py"
OV3D01_PARAMETERS = (ROOT / "validation" / "cases" / "2013-04-11" /
                     "parameters.json")
OV3D01_PREPARER = (ROOT / "validation" / "cases" / "2013-04-11" /
                   "prepare_case.py")
OV3D01_FIT_TEMPLATE = (ROOT / "validation" / "cases" / "2013-04-11" /
                       "event_fit.template.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ValidationRunnerTests(unittest.TestCase):
    def run_runner(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run([sys.executable, str(RUNNER), *arguments],
                              cwd=str(ROOT), text=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, check=False)

    def run_preparer(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run([sys.executable, str(OV3D01_PREPARER),
                               *arguments], cwd=str(ROOT), text=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, check=False)

    def run_native_matrix(
            self, *arguments: str) -> subprocess.CompletedProcess[str]:
        """Invoke the matrix runner exactly as an operator would.

        The test below supplies a process-level launcher and executable rather
        than importing runner internals. This exercises argv construction,
        discovery, environment propagation, report parsing, and hashing across
        the same subprocess boundary used by an actual MPI campaign.
        """
        return subprocess.run([sys.executable, str(NATIVE_MATRIX), *arguments],
                              cwd=str(ROOT), text=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, check=False)

    def complete_ov3d01_fit(self) -> dict:
        """Return a structurally complete synthetic fit for mechanics tests.

        Values intentionally reproduce the generic example where possible;
        this object tests derivation/refusal/checksum mechanics and is never
        registered or presented as observational event evidence.
        """
        au = 149597870700.0
        payload = json.loads(OV3D01_FIT_TEMPLATE.read_text(encoding="utf-8"))
        payload["coordinate_transform"] = {
            "name": "event-frozen-carrington-fixture",
            "model_from_event_frozen_carrington": [
                [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "provenance": "unit-test identity rotation",
        }
        payload["ambient"] = {
            "number_density_at_one_au_m3": 5.0e6,
            "magnetic_field_magnitude_at_one_au_t": 5.0e-9,
            "proton_temperature_k": 1.0e5,
            "magnetic_polarity": 1,
            "averaging_start_utc": "2013-04-10T18:55:00Z",
            "averaging_end_utc": "2013-04-11T06:55:00Z",
            "source_file": "unit-test-omni.csv",
            "source_sha256": "0" * 64,
        }
        payload["swcme"] = {
            "solar_rotation_rate_rad_per_s": 2.865e-6,
            "adiabatic_index": 5.0 / 3.0,
            "thermodynamic_closure": "proton-only",
            "alpha_to_proton_ratio": 0.0,
            "electron_temperature_k": 1.0e5,
            "alpha_temperature_k": 1.0e5,
            "drag_coefficient_per_km": 8.0e-8,
            "sheath_thickness_at_one_au_m": 0.12 * au,
            "ejecta_thickness_at_one_au_m": 0.22 * au,
            "shock_smoothing_width_at_one_au_m": 0.01 * au,
            "leading_edge_smoothing_width_at_one_au_m": 0.02 * au,
            "trailing_edge_smoothing_width_at_one_au_m": 0.03 * au,
            "sheath_ramp_power": 2.0,
            "leading_edge_speed_factor": 1.12,
            "ejecta_density_factor": 0.5,
            "ejecta_speed_factor": 0.8,
            "valid_from_s": 3600.0,
            "valid_until_s": 86400.0,
            "fit_provenance": "unit-test fixture",
        }
        payload["source"] = {
            "physical_particle_rate_per_s": 1.0e30,
            "injection_efficiency": 1.0e-4,
            "maximum_energy_j": 1.602176634e-11,
            "reference_energy_mev": 10.0,
            "samples_per_step_per_species": 1000,
            "relative_source_weight_per_area": 1.0,
            "normalization_provenance": "unit-test fixture",
        }
        payload["mesh"] = {
            "global_cell_size_m": 3.7399467675e10,
            "minimum_cell_size_m": 1.495978707e9,
            "maximum_level": 7,
            "memory_budget_bytes": 68719476736,
            "solar_surface_cell_size_m": 1.495978707e9,
            "solar_transition_outer_radius_m": 3.7399467675e10,
            "refinement_tube_radius_at_reference_m": 4.487936121e9,
            "refinement_tube_center_cell_size_m": 1.495978707e9,
            "active_tube_radius_at_reference_m": 7.479893535e9,
            "active_tube_buffer_blocks": 1,
            "observer_collection_radius_m": 1.495978707e9,
            "convergence_provenance": "unit-test fixture",
        }
        payload["population_control"] = {
            "mode": "split-merge",
            "minimum_particles_per_cell_per_species": 16,
            "target_particles_per_cell_per_species": 24,
            "maximum_particles_per_cell_per_species": 32,
            "cadence_steps": 1,
            "convergence_provenance": "unit-test fixture",
        }
        payload["turbulence"] = {
            "delta_b_over_b": 0.3,
            "correlation_length_m": 4.487936121e9,
            "provenance": "unit-test fixture",
        }
        payload["run"] = {
            "time_step_s": 1.0,
            "duration_s": 259200.0,
            "campaign_seed": 4701,
            "observer_cadence_s": 60.0,
            "observer_energy_bins": 48,
            "provenance": "unit-test fixture",
        }
        return payload

    def test_registry_lists_all_evidence_classes(self) -> None:
        completed = self.run_runner("--list")
        self.assertEqual(completed.returncode, 0, completed.stdout)
        for case_id in ("NAT3D01", "MPI3D01", "XM3D01", "XM3D05",
                        "OV3D01", "OV3D04"):
            self.assertIn(case_id, completed.stdout)

    def test_ov3d01_blueprint_separates_known_and_unresolved_physics(self) -> None:
        """Guard the event law and the explicit no-guesses campaign boundary."""
        payload = json.loads(OV3D01_PARAMETERS.read_text(encoding="utf-8"))
        self.assertEqual(payload["case_id"], "OV3D01")
        self.assertEqual(payload["execution_status"],
                         "blueprint-needs-event-fit")
        values = {item["name"]: item["value"]
                  for item in payload["published_parameters"]}
        self.assertEqual(values["mean_free_path_reference_m"], 44879361210.0)
        self.assertEqual(values["mean_free_path_reference_rigidity_v"], 1.0e9)
        self.assertEqual(values["mean_free_path_radial_exponent"], 1.0)
        self.assertAlmostEqual(values["mean_free_path_rigidity_exponent"],
                               1.0 / 3.0)
        unresolved = {item["name"] for item in payload["unresolved_inputs"]}
        self.assertIn("dbm_drag_and_shock_arrival", unresolved)
        self.assertIn("source_rate_and_macroparticle_weight", unresolved)
        self.assertIn("seed_boundary_rate_mapping", unresolved)
        mappings = {item["name"]: item["srcsep3d_mapping"]
                    for item in payload["published_parameters"]}
        self.assertIn("fixed-phase-space-power-law",
                      mappings["seed_momentum_power_law_index"])
        self.assertIn("phase_space_power_index=5",
                      mappings["seed_momentum_power_law_index"])

    def test_ov3d01_preparer_refuses_missing_fit_and_builds_hashed_artifacts(self) -> None:
        incomplete = self.run_preparer(
            "status", "--fit", str(OV3D01_FIT_TEMPLATE))
        self.assertEqual(incomplete.returncode, 1, incomplete.stdout)
        self.assertIn("ambient.number_density_at_one_au_m3", incomplete.stdout)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fit = root / "fit.json"
            fit.write_text(json.dumps(
                self.complete_ov3d01_fit(), indent=2) + "\n",
                encoding="utf-8")
            rendered = root / "ov3d01.in"
            completed = self.run_preparer(
                "render-input", "--fit", str(fit), "--output", str(rendered))
            self.assertEqual(completed.returncode, 0, completed.stdout)
            deck = rendered.read_text(encoding="utf-8")
            self.assertIn("start_mode = explicit", deck)
            self.assertIn("mode = parker-tube", deck)
            self.assertIn("spectrum_model = fixed-phase-space-power-law", deck)
            self.assertIn("phase_space_power_index = 5", deck)
            self.assertIn("mean_free_path_model = radial-rigidity-power-law", deck)
            self.assertNotIn("[observer.inner]", deck)

            model = root / "prepared-model.csv"
            reference = root / "prepared-reference.csv"
            series = ("coordinate_si,value_si\n"
                      "0,1\n60,2\n120,4\n180,2\n240,1\n")
            model.write_text(series, encoding="utf-8")
            reference.write_text(series, encoding="utf-8")
            evidence_root = root / "evidence"
            case_dir = evidence_root / "OV3D01"
            bundled = self.run_preparer(
                "build-evidence", "--model-csv", str(model),
                "--reference-csv", str(reference), "--output-dir",
                str(case_dir), "--coordinate-name", "time",
                "--coordinate-units", "s", "--value-name",
                "differential_intensity", "--value-units",
                "m^-2 s^-1 sr^-1 J^-1", "--model-provenance",
                "unit-test model fixture", "--reference-provenance",
                "unit-test observation fixture")
            self.assertEqual(bundled.returncode, 0, bundled.stdout)
            manifest = json.loads(
                (case_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["model_sha256"],
                             sha256(case_dir / "model.csv"))
            self.assertEqual(manifest["reference_sha256"],
                             sha256(case_dir / "reference.csv"))

            output = root / "validation-output"
            validated = self.run_runner(
                "--case", "OV3D01", "--evidence-root", str(evidence_root),
                "--output-dir", str(output))
            self.assertEqual(validated.returncode, 0, validated.stdout)
            result = json.loads(
                (output / "OV3D01" / "result.json").read_text())
            self.assertEqual(result["status"], "PASS")

    def test_native_matrix_is_capability_gated_and_hashes_its_input(self) -> None:
        """Prove that native evidence cannot be fabricated by the runner.

        These tiny scripts only emulate the *test-host protocol*. They do not
        emulate AMPS or count as native physics evidence. Their purpose is to
        make the runner's fail-closed contract testable on a machine without
        MPI: the launcher explicitly forwards the requested rank count and the
        host writes the same JSON schema as the shared C++ registry.
        """
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            launcher = root / "launcher.py"
            launcher.write_text(
                "#!/usr/bin/env python3\n"
                "import os\n"
                "import subprocess\n"
                "import sys\n"
                "if len(sys.argv) < 4 or sys.argv[1] != '--ranks':\n"
                "    raise SystemExit(91)\n"
                "environment = dict(os.environ)\n"
                "environment['FAKE_MPI_RANKS'] = sys.argv[2]\n"
                "raise SystemExit(subprocess.call(sys.argv[3:], "
                "env=environment))\n",
                encoding="utf-8")
            launcher.chmod(0o755)

            host = root / "native-host.py"

            def write_host(advertised: tuple[str, ...]) -> None:
                listing = "".join(
                    f"print('{case_id}|linked|fixture')\n"
                    for case_id in advertised)
                host.write_text(
                    "#!/usr/bin/env python3\n"
                    "import json\n"
                    "import os\n"
                    "from pathlib import Path\n"
                    "import sys\n"
                    "if '--list-tests' in sys.argv:\n" +
                    "    " + listing.replace("\n", "\n    ").rstrip() + "\n"
                    "    raise SystemExit(0)\n"
                    "def value(option):\n"
                    "    return sys.argv[sys.argv.index(option) + 1]\n"
                    "case_id = value('--test')\n"
                    "input_path = Path(value('--test-input'))\n"
                    "report = Path(value('--test-json'))\n"
                    "artifact_dir = Path(value('--artifact-directory'))\n"
                    "if not input_path.is_file() or report.parent != "
                    "artifact_dir:\n"
                    "    raise SystemExit(92)\n"
                    "message = ('fixture ranks=' + "
                    "os.environ.get('FAKE_MPI_RANKS', '') + ' threads=' + "
                    "os.environ.get('OMP_NUM_THREADS', ''))\n"
                    "payload = {\n"
                    "    'schema': 'srcsep-component-tests-v1',\n"
                    "    'exit_code': 0,\n"
                    "    'totals': {'passed': 1, 'failed': 0, "
                    "'skipped': 0, 'errors': 0},\n"
                    "    'results': [{'id': case_id, 'status': 'PASS', "
                    "'message': message, 'metrics': [], 'artifacts': []}],\n"
                    "}\n"
                    "report.write_text(json.dumps(payload) + '\\n', "
                    "encoding='utf-8')\n"
                    "raise SystemExit(0)\n",
                    encoding="utf-8")
                host.chmod(0o755)

            test_input = root / "qualification.in"
            test_input.write_text("schema_version = 4\n", encoding="utf-8")
            launcher_template = (
                f"{sys.executable} {launcher} --ranks {{ranks}}")

            # A host missing one required callback must be rejected before any
            # matrix cell is accepted or a summary is published.
            write_host(("MPI3D01",))
            rejected_output = root / "rejected"
            rejected = self.run_native_matrix(
                "--amps", str(host), "--test-input", str(test_input),
                "--profile", "small", "--launcher", launcher_template,
                "--output-dir", str(rejected_output), "--timeout", "10")
            self.assertEqual(rejected.returncode, 2, rejected.stdout)
            self.assertIn("does not advertise required native case(s): MPI3D02",
                          rejected.stdout)
            self.assertFalse((rejected_output / "native-matrix.json").exists())

            # A protocol-complete host exercises all 2-rank x 2-thread x
            # 2-case cells in the small profile. The resulting hashes prove
            # the summary owns the exact input and executable it launched.
            write_host(("MPI3D01", "MPI3D02"))
            accepted_output = root / "accepted"
            accepted = self.run_native_matrix(
                "--amps", str(host), "--test-input", str(test_input),
                "--profile", "small", "--launcher", launcher_template,
                "--output-dir", str(accepted_output), "--timeout", "10")
            self.assertEqual(accepted.returncode, 0, accepted.stdout)
            summary = json.loads(
                (accepted_output / "native-matrix.json").read_text(
                    encoding="utf-8"))
            self.assertEqual(summary["schema"],
                             "srcsep3d-native-matrix-v2")
            self.assertEqual(summary["test_input_sha256"], sha256(test_input))
            self.assertEqual(summary["executable_sha256"], sha256(host))
            self.assertEqual(summary["totals"]["PASS"], 8)
            self.assertEqual(len(summary["records"]), 8)
            self.assertEqual(
                {row["ranks"] for row in summary["records"]}, {1, 2})
            self.assertEqual(
                {row["threads"] for row in summary["records"]}, {1, 2})
            for row in summary["records"]:
                self.assertIn("--test-input", row["command"])
                self.assertTrue(Path(row["native_report"]).is_file())

    def test_missing_external_evidence_is_skip(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            completed = self.run_runner("--case", "XM3D01", "--output-dir",
                                        str(output))
            self.assertEqual(completed.returncode, 0, completed.stdout)
            result = json.loads((output / "XM3D01" / "result.json").read_text())
            self.assertEqual(result["status"], "SKIP")

    def test_linked_case_requires_an_explicit_input_deck(self) -> None:
        """A linked executable alone is insufficient native provenance."""
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            completed = self.run_runner(
                "--case", "NAT3D01", "--amps", sys.executable,
                "--output-dir", str(output))
            self.assertEqual(completed.returncode, 2, completed.stdout)
            result = json.loads(
                (output / "NAT3D01" / "result.json").read_text())
            self.assertEqual(result["status"], "ERROR")
            self.assertIn("require --test-input", result["message"])

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

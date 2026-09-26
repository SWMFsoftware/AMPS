#!/usr/bin/env python3
"""Reference and end-to-end tests for the Roadmap Step-8 campaign layer."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import stat
import sys
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
SRC_EARTH = HERE.parents[1]
sys.path.insert(0, str(SRC_EARTH))

from standalone_campaign.adapters import (  # noqa: E402
    adapt_observations,
    cadence_average,
    compare_with_observations,
    validate_instrument_response,
)
from standalone_campaign.campaign import (  # noqa: E402
    CampaignError,
    cache_resources,
    load_and_validate_manifest,
    read_prediction_csv,
    sha256_file,
)
from standalone_campaign.run_campaign import (  # noqa: E402
    _write_extracted_predictions,
    main as run_campaign,
)
from standalone_campaign.release_gate import evaluate_release_gate  # noqa: E402


UTC = timezone.utc


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def resource(identifier: str, kind: str, path: Path) -> dict:
    return {
        "id": identifier,
        "kind": kind,
        "path": path.name,
        "sha256": sha256_file(path),
        # The unit fixture is generated locally from the source text below.
        # Production manifests should name the archive, DOI, or processing
        # record from which each immutable file was obtained.
        "provenance": "UStandaloneCampaign deterministic local fixture",
    }


def make_fixture(
    root: Path,
    unresolved: float = 0.005,
    prediction: bool = True,
    event_id: str = "O1",
    role: str = "VALIDATION",
):
    """Create a three-epoch campaign with a deterministic fake AMPS executable."""

    template = write(root / "input.in", "MODEL {{FIELD_MODEL}}\nEPOCH {{EPOCH_UTC}}\nDRIVER {{DRIVER_PATH}}\n")
    driver = write(root / "driver.dat", "# deterministic driver descriptor\n")
    boundary = write(root / "boundary.csv", "energy_mev,intensity\n10,1\n100,0.1\n")
    ephemeris = write(root / "ephemeris.csv", "utc,x_km,y_km,z_km\n2020-01-01T00:00:00Z,42164,0,0\n")
    attitude = write(root / "attitude.csv", "utc,q0,q1,q2,q3\n2020-01-01T00:00:00Z,1,0,0,0\n")
    response = write(root / "response.csv", "channel,energy_min_mev,energy_max_mev,response\nP5,38,82,1\n")
    observation = write(
        root / "goes.csv",
        "utc,spacecraft,channel,east_west_ratio,east_quality_flag,west_quality_flag,quality_status\n"
        "2020-01-01T00:05:00Z,GOES13,P5,2.0,0,0,VALID\n",
    )

    resources = [
        resource("template", "input_template", template),
        resource("driver", "driver", driver),
        resource("boundary", "boundary_spectrum", boundary),
        resource("ephemeris", "ephemeris", ephemeris),
        resource("attitude", "attitude", attitude),
        resource("response", "instrument_response", response),
        resource("observation", "observation", observation),
    ]
    manifest = {
        "schema_version": "earth-standalone-campaign/v1",
        "campaign_id": "STEP8_TEST",
        "event_id": event_id,
        "role": role,
        "frozen": True,
        "normalization_policy": "SHARED_EVENT_BOUNDARY_NO_PLATFORM_SCALE",
        "event": {
            "start_utc": "2020-01-01T00:00:00Z",
            "end_utc": "2020-01-01T00:10:00Z",
            "field_cadence_seconds": 300,
        },
        "models": ["DIPOLE"],
        "products": ["cutoff", "directional_access", "spectrum", "density", "flux", "detector_rate"],
        "exclusions": [],
        "resources": resources,
        "instruments": [{
            "id": "GOES_EPEAD",
            "adapter": "GOES_EPEAD_DIRECTIONAL",
            "observation_resource": "observation",
            "response_resource": "response",
            "cadence_seconds": 300,
            "quantity": "directional_ratio",
            "units": "1",
            "comparison_required": True,
        }],
        "execution": {
            "input_template_resource": "template",
            "command": ["{amps}", "{input}"],
            "mpi_ranks": 1,
            "threads": 1,
            "expected_artifacts": [
                "standalone_run_manifest.json",
                "gridless_termination_summary.dat",
                "convergence_metrics.json",
                "campaign_prediction.csv",
            ],
            "prediction_artifact": "campaign_prediction.csv",
        },
        "validation": {
            "numerical_gates": [
                {
                    "id": "unresolved_support",
                    "source": "termination_summary",
                    "metric": "max_unresolved_fraction",
                    "operator": "<=",
                    "threshold": 0.01,
                    "required": True,
                },
                {
                    "id": "energy_convergence",
                    "source": "json_artifact",
                    "path": "convergence_metrics.json",
                    "metric": "energy_relative_change",
                    "operator": "<=",
                    "threshold": 0.02,
                    "required": True,
                },
                {
                    "id": "angular_convergence",
                    "source": "json_artifact",
                    "path": "convergence_metrics.json",
                    "metric": "angular_relative_change",
                    "operator": "<=",
                    "threshold": 0.02,
                    "required": True,
                },
            ]
        },
    }
    manifest_path = root / "campaign.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    fake = root / "fake_amps.py"
    prediction_code = """
with open('campaign_prediction.csv', 'w', encoding='utf-8') as out:
    out.write('utc,instrument_id,platform,channel,direction,quantity,units,value,lower,upper\\n')
    out.write('%s,GOES_EPEAD,GOES13,P5,EAST_OVER_WEST,directional_ratio,1,2.0,,\\n' % epoch)
""" if prediction else ""
    fake.write_text(
        "#!/usr/bin/env python3\n"
        "import json, pathlib, sys\n"
        "text = pathlib.Path(sys.argv[1]).read_text(encoding='utf-8')\n"
        "epoch = [line.split()[1] for line in text.splitlines() if line.startswith('EPOCH ')][0]\n"
        "pathlib.Path('standalone_run_manifest.json').write_text(json.dumps({'epoch': epoch})+'\\n')\n"
        "pathlib.Path('convergence_metrics.json').write_text(json.dumps({"
        "'energy_relative_change': 0.01, 'angular_relative_change': 0.01})+'\\n')\n"
        "pathlib.Path('gridless_termination_summary.dat').write_text(\n"
        "    'VARIABLES=\\\"N_sampled\\\" \\\"N_resolved\\\" \\\"unresolved_fraction\\\" \\\"N_OUTER_BOUNDARY_ALLOWED\\\"\\n'\n"
        "    'ZONE T=\\\"termination\\\" I=1 F=POINT\\n'\n"
        "    '100 99 %.17g 50\\n' % " + repr(unresolved) + ")\n"
        + prediction_code,
        encoding="utf-8",
    )
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    return manifest_path, fake


class ManifestTests(unittest.TestCase):
    def test_manifest_hash_and_cache_are_strict(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, _ = make_fixture(root)
            manifest = load_and_validate_manifest(manifest_path)
            cached = cache_resources(manifest["_verified_resources"], root / "cache")
            self.assertEqual(set(cached), {item["id"] for item in manifest["resources"]})
            self.assertTrue(all(sha256_file(cached[item["id"]]) == item["sha256"] for item in manifest["resources"]))

            # A changed reference must block the campaign before execution.
            write(root / "driver.dat", "changed\n")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

    def test_network_resource_and_platform_scale_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["resources"][1]["path"] = "https://example.invalid/driver.dat"
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["execution"]["prediction_artifact"] = "../../escaped.csv"
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["normalization_policy"] = "FIT_EACH_PLATFORM"
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

    def test_frozen_numerical_contract_cannot_be_weakened_or_omitted(self):
        """Step 8 may orchestrate Step 7, but may not loosen its gates."""

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            gates = {item["id"]: item for item in data["validation"]["numerical_gates"]}
            gates["unresolved_support"]["threshold"] = 0.0100001
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["validation"]["numerical_gates"] = [
                item for item in data["validation"]["numerical_gates"]
                if item["id"] != "angular_convergence"
            ]
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

    def test_provenance_and_preregistered_exclusions_are_mandatory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            del data["resources"][0]["provenance"]
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["instruments"][0]["quantity"] = "integral_flux"
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)

            manifest_path, _ = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["exclusions"] = [{
                "id": "late_drop", "reason": "added after inspection",
                "preregistered": False,
            }]
            manifest_path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(CampaignError):
                load_and_validate_manifest(manifest_path)


class AdapterAndCadenceTests(unittest.TestCase):
    def test_released_o1_o2_references_are_accepted(self):
        """Exercise the committed observation assets, not only small fixtures."""

        pamela = adapt_observations(
            SRC_EARTH / "test" / "C9" / "reference_C9_pamela_table_s1.csv",
            {"id": "PAMELA", "adapter": "PAMELA_CUTOFF", "units": "deg_AACGM"},
        )
        poes = adapt_observations(
            SRC_EARTH / "test" / "C10" / "reference_C10_poes_meped_boundary.csv.gz",
            {"id": "POES", "adapter": "POES_METOP_MEPED_CUTOFF", "units": "deg_AACGM"},
        )
        goes = adapt_observations(
            SRC_EARTH / "test" / "C19" / "data" / "reference_C19_goes_epead_ew.csv.gz",
            {
                "id": "GOES", "adapter": "GOES_EPEAD_DIRECTIONAL",
                "cadence_seconds": 300, "units": "1",
            },
        )
        self.assertGreater(len(pamela), 100)
        self.assertGreater(len(poes), 100)
        self.assertGreater(len(goes), 100)
        self.assertAlmostEqual(pamela[0]["value"], 60.72)
        self.assertEqual(poes[0]["quality"], "VALID")
        self.assertEqual(goes[0]["quantity"], "directional_ratio")

    def test_rept_adapter_reference_row(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = write(
                Path(temporary) / "rept.csv",
                "time_utc,spacecraft,energy_mev,proton_flux,lower,upper,quality_status,units\n"
                "2017-09-08T00:00:00Z,RBSP-A,20,4.0,3.0,5.0,VALID,cm^-2_s^-1_sr^-1_MeV^-1\n",
            )
            rows = adapt_observations(
                path,
                {
                    "id": "REPT", "adapter": "REPT_PROTON_SPECTRUM",
                    "cadence_seconds": 300,
                    "units": "cm^-2_s^-1_sr^-1_MeV^-1",
                },
            )
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["platform"], "RBSP-A")
            self.assertEqual(rows[0]["quantity"], "differential_flux")
            self.assertAlmostEqual(rows[0]["value"], 4.0)

    def test_detector_response_schema_and_bounds(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid = write(
                root / "response.csv",
                "channel,energy_min_mev,energy_max_mev,relative_response\n"
                "P5,38,82,1.0\n",
            )
            config = {
                "id": "GOES", "adapter": "GOES_EPEAD_DIRECTIONAL", "units": "1",
            }
            summary = validate_instrument_response(valid, config)
            self.assertEqual(summary["status"], "PASS")
            self.assertEqual(summary["row_count"], 1)

            invalid = write(
                root / "invalid.csv",
                "channel,energy_min_mev,energy_max_mev,relative_response\n"
                "P5,82,38,1.0\n",
            )
            with self.assertRaises(CampaignError):
                validate_instrument_response(invalid, config)

    def test_goes_adapter_and_linear_exposure_reference(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, _ = make_fixture(root)
            manifest = load_and_validate_manifest(manifest_path)
            instrument = manifest["instruments"][0]
            rows = adapt_observations(root / "goes.csv", instrument)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["quantity"], "directional_ratio")
            self.assertAlmostEqual(rows[0]["value"], 2.0)

            # y(t)=1+t/300 over [0,600] has exact mean 2.0.  This is an
            # independent analytic reference for the trapezoidal exposure fold.
            samples = [
                (datetime(2020, 1, 1, 0, 0, tzinfo=UTC), 1.0),
                (datetime(2020, 1, 1, 0, 5, tzinfo=UTC), 2.0),
                (datetime(2020, 1, 1, 0, 10, tzinfo=UTC), 3.0),
            ]
            average = cadence_average(samples, samples[0][0], samples[-1][0])
            self.assertAlmostEqual(average, 2.0, places=14)

            # An observation exposure outside the simulated cadence must not
            # be filled by nearest-neighbor or linear extrapolation.
            with self.assertRaises(CampaignError):
                cadence_average(
                    samples,
                    samples[0][0] - (samples[1][0] - samples[0][0]),
                    samples[-1][0],
                )

    def test_comparison_reference_metrics(self):
        observations = [{
            "instrument_id": "I", "platform": "P", "channel": "C",
            "direction": "EAST_OVER_WEST", "quantity": "directional_ratio",
            "units": "1",
            "start_utc": "2020-01-01T00:00:00Z", "end_utc": "2020-01-01T00:10:00Z",
            "value": 2.0,
        }]
        predictions = [
            {"instrument_id": "I", "platform": "P", "channel": "C",
             "direction": "EAST_OVER_WEST", "quantity": "directional_ratio",
             "units": "1",
             "utc": "2020-01-01T00:00:00Z", "value": 1.0, "lower": None, "upper": None},
            {"instrument_id": "I", "platform": "P", "channel": "C",
             "direction": "EAST_OVER_WEST", "quantity": "directional_ratio",
             "units": "1",
             "utc": "2020-01-01T00:10:00Z", "value": 3.0, "lower": None, "upper": None},
        ]
        rows, metrics = compare_with_observations(predictions, observations)
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0]["modeled"], 2.0, places=14)
        self.assertAlmostEqual(metrics["directional_log10_rmse"], 0.0, places=14)
        self.assertEqual(metrics["directional_sign_fraction"], 1.0)
        self.assertEqual(metrics["comparison_fraction"], 1.0)

        # Units are part of the identity key.  A numerically equal value with
        # incompatible units is missing evidence, never an implicit conversion.
        incompatible = [dict(predictions[0], units="percent")]
        rows, metrics = compare_with_observations(incompatible, observations)
        self.assertEqual(rows, [])
        self.assertEqual(metrics["missing_model_coverage_count"], 1)


class RunnerTests(unittest.TestCase):
    def test_tecplot_ratio_extractor_reference(self):
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            write(
                run_dir / "gridless_points_flux.dat",
                'TITLE="rates"\n'
                'VARIABLES="X_km" "R_EAST_s1" "R_WEST_s1"\n'
                'ZONE T="flux" I=1 F=POINT\n'
                '42164 6 3\n',
            )
            execution = {
                "prediction_artifact": "campaign_prediction.csv",
                "prediction_extractors": [{
                    "instrument_id": "GOES", "artifact": "gridless_points_flux.dat",
                    "row_index": 0, "operation": "RATIO",
                    "numerator_column": "R_EAST_s1", "denominator_column": "R_WEST_s1",
                    "platform": "GOES13", "channel": "P5",
                    "direction": "EAST_OVER_WEST", "quantity": "directional_ratio",
                    "units": "1",
                    "expected_variables": ["X_km", "R_EAST_s1", "R_WEST_s1"],
                    "expected_row_count": 1,
                }],
            }
            epoch = datetime(2020, 1, 1, 0, 0, tzinfo=UTC)
            _write_extracted_predictions(run_dir, epoch, execution)
            rows = read_prediction_csv(run_dir / "campaign_prediction.csv")
            self.assertEqual(len(rows), 1)
            self.assertAlmostEqual(rows[0]["value"], 2.0)
            self.assertEqual(rows[0]["utc"], "2020-01-01T00:00:00Z")

            # A changed producer schema or point count is a scientific-input
            # change.  The extractor must not continue by column-name guessing.
            execution["prediction_extractors"][0]["expected_variables"] = [
                "X_km", "R_WEST_s1", "R_EAST_s1",
            ]
            with self.assertRaises(CampaignError):
                _write_extracted_predictions(run_dir, epoch, execution)
            execution["prediction_extractors"][0]["expected_variables"] = [
                "X_km", "R_EAST_s1", "R_WEST_s1",
            ]
            execution["prediction_extractors"][0]["expected_row_count"] = 2
            with self.assertRaises(CampaignError):
                _write_extracted_predictions(run_dir, epoch, execution)

    def test_end_to_end_pass_and_verified_restart(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, fake = make_fixture(root)
            output = root / "output"
            arguments = [
                "--manifest", str(manifest), "--amps", str(fake),
                "--output-dir", str(output),
            ]
            self.assertEqual(run_campaign(arguments), 0)
            summary = json.loads((output / "campaign_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "PASS")
            self.assertEqual(summary["run_pass_count"], 3)
            self.assertEqual(summary["observation_comparison_fail_count"], 0)
            comparison = output / "comparisons" / "DIPOLE" / "GOES_EPEAD" / "comparison.csv"
            self.assertGreater(len(comparison.read_text(encoding="utf-8").splitlines()), 1)

            self.assertEqual(run_campaign(arguments + ["--restart"]), 0)
            summary = json.loads((output / "campaign_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["restart_skip_count"], 3)

            # Artifact hashes, not only a previous PASS label, control restart.
            # Corrupt one run and prove that only that run is recomputed.
            first_manifest = sorted(output.glob("runs/*/*/standalone_run_manifest.json"))[0]
            first_manifest.write_text("tampered\n", encoding="utf-8")
            self.assertEqual(run_campaign(arguments + ["--restart"]), 0)
            summary = json.loads((output / "campaign_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["restart_skip_count"], 2)
            repaired_status = json.loads(
                (first_manifest.parent / "run_status.json").read_text(encoding="utf-8")
            )
            self.assertEqual(repaired_status["restart_action"], "EXECUTED")

    def test_dry_run_is_never_reported_as_scientific_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, fake = make_fixture(root)
            output = root / "dry"
            self.assertEqual(run_campaign([
                "--manifest", str(manifest), "--amps", str(fake),
                "--output-dir", str(output), "--dry-run",
            ]), 0)
            summary = json.loads((output / "campaign_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "NOT_RUN")
            self.assertEqual(summary["run_pass_count"], 0)
            self.assertEqual(summary["run_not_executed_count"], 3)
            self.assertFalse(summary["all_required_numerical_gates_pass"])

    def test_process_launch_failure_writes_fail_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, fake = make_fixture(root)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["execution"]["command"] = [
                "step8-command-that-does-not-exist", "{amps}", "{input}",
            ]
            manifest_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
            output = root / "launch_failure"
            self.assertEqual(run_campaign([
                "--manifest", str(manifest_path), "--amps", str(fake),
                "--output-dir", str(output),
            ]), 2)
            summary = json.loads((output / "campaign_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "FAIL")
            status_files = sorted(output.glob("runs/*/*/run_status.json"))
            self.assertEqual(len(status_files), 3)
            for status_path in status_files:
                status = json.loads(status_path.read_text(encoding="utf-8"))
                self.assertEqual(status["status"], "FAIL")
                self.assertEqual(status["return_code"], 127)
                self.assertIn("FileNotFoundError", status["launch_error"])

    def test_unresolved_gate_is_not_relaxed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, fake = make_fixture(root, unresolved=0.02)
            output = root / "output"
            code = run_campaign([
                "--manifest", str(manifest), "--amps", str(fake),
                "--output-dir", str(output),
            ])
            self.assertEqual(code, 2)
            status_files = sorted(output.glob("runs/*/*/run_status.json"))
            self.assertEqual(len(status_files), 3)
            for path in status_files:
                status = json.loads(path.read_text(encoding="utf-8"))
                gate = {item["id"]: item for item in status["numerical_gates"]}
                self.assertEqual(gate["unresolved_support"]["status"], "FAIL")
                self.assertEqual(gate["unresolved_support"]["threshold"], 0.01)

    def test_missing_comparison_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, fake = make_fixture(root, prediction=False)
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            data["execution"]["expected_artifacts"].remove("campaign_prediction.csv")
            manifest_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
            output = root / "output"
            code = run_campaign([
                "--manifest", str(manifest_path), "--amps", str(fake),
                "--output-dir", str(output),
            ])
            self.assertEqual(code, 2)
            summary = json.loads((output / "campaign_summary.json").read_text(encoding="utf-8"))
            self.assertGreater(summary["observation_comparison_fail_count"], 0)

    def test_o1_o2_o4_release_gate_requires_complete_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            holdout_manifest, _ = make_fixture(
                root, event_id="O4", role="HOLDOUT"
            )

            def completed(event_id: str) -> dict:
                return {
                    "campaign_id": "STEP8_" + event_id,
                    "event_id": event_id,
                    "status": "PASS",
                    "role": "VALIDATION",
                    "frozen": True,
                    "normalization_policy": "SHARED_EVENT_BOUNDARY_NO_PLATFORM_SCALE",
                    "all_required_numerical_gates_pass": True,
                    "observation_comparison_count": 1,
                    "observation_comparison_fail_count": 0,
                    "observation_comparison_row_count": 1,
                    "required_observation_comparison_count": 1,
                    "required_observation_comparison_fail_count": 0,
                    "run_count": 3,
                    "run_pass_count": 3,
                    "run_fail_count": 0,
                    "run_not_executed_count": 0,
                    "manifest_sha256": "a" * 64,
                    "executable_sha256": "b" * 64,
                    "resource_count": 7,
                }

            o1 = write(root / "o1.json", json.dumps(completed("O1")) + "\n")
            o2 = write(root / "o2.json", json.dumps(completed("O2")) + "\n")
            decision = evaluate_release_gate(o1, o2, holdout_manifest)
            self.assertEqual(decision["status"], "PASS")
            self.assertEqual(
                [item["status"] for item in decision["decisions"]],
                ["PASS", "PASS", "READY"],
            )

            missing = completed("O2")
            missing["observation_comparison_count"] = 0
            o2.write_text(json.dumps(missing) + "\n", encoding="utf-8")
            decision = evaluate_release_gate(o1, o2, holdout_manifest)
            self.assertEqual(decision["status"], "FAIL")
            self.assertFalse(decision["decisions"][1]["checks"]["observation_comparisons"])


if __name__ == "__main__":
    unittest.main(verbosity=2)

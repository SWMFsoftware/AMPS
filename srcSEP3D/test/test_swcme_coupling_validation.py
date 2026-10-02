#!/usr/bin/env python3
"""Independent exact-crossing and rejection checks for the evidence consumer.

Fixtures exercise report mechanics only. No test here claims to run SWCME,
MPI or an observed CME. Actual native verification remains a campaign gate.
"""
import csv
import io
from types import SimpleNamespace
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cme_validation", ROOT / "validation/run_swcme_coupling_validation.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
prepare_spec=importlib.util.spec_from_file_location("cme_reference_preparer",ROOT/"validation/cases/2012-07-12/prepare_reference_data.py")
preparer=importlib.util.module_from_spec(prepare_spec)
prepare_spec.loader.exec_module(preparer)


class CouplingConsumerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle, self.output = self.root / "bundle", self.root / "output"
        runner.write_demo(self.bundle)

    def manifest(self, **changes):
        path = self.bundle / "manifest.json"
        data = runner.read_json(path)
        data.update(changes)
        runner.atomic_json(path, data)

    def rewrite(self, file, mutator, descriptor):
        path = self.bundle / file
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            columns, rows = reader.fieldnames, list(reader)
        mutator(rows)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader(); writer.writerows(rows)
        data = runner.read_json(self.bundle / "manifest.json")
        data[descriptor]["sha256"] = runner.sha256(path)
        self.manifest(**data)

    def evaluate(self):
        # Tests below check evidence logic independently of renderer speed.
        # One separate test renders both publication formats through real Agg.
        with patch.object(runner, "plot_comparison", return_value=[]):
            return runner.evaluate(self.bundle, self.output)

    def test_exact_ballistic_crossing_and_no_synthetic_pass(self):
        result = self.evaluate()
        exact = (149597870700.0 - 20 * 695700000.0) / 700000.0
        self.assertAlmostEqual(result["metrics"]["modeled_one_au_time_s"], exact, places=8)
        self.assertLess(abs(result["metrics"]["arrival_error_s"]), 1e-5)
        self.assertEqual(result["status"], "SKIP")
        self.assertFalse(result["production_release_qualified"])
        self.assertFalse(result["full_plasma_coupling_validated"])

    def test_particle_and_rank_failures_are_explicit(self):
        def change(rows):
            rows[10]["particle_count"] = "1"
            rows[100]["mpi_radius_spread_m"] = "2"
        self.rewrite("history.csv", change, "model")
        result = self.evaluate()
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("particle population/injection was nonzero", result["violations"])
        self.assertIn("MPI ranks disagree", result["violations"])

    def test_early_termination_does_not_extrapolate_arrival(self):
        self.rewrite("history.csv", lambda rows: rows.__delitem__(slice(1800, None)), "model")
        result = self.evaluate()
        self.assertEqual(result["status"], "FAIL")
        self.assertIsNone(result["metrics"]["modeled_one_au_time_s"])
        self.assertIsNone(result["metrics"]["arrival_error_s"])
        self.assertLess(result["metrics"]["holdout_coverage"], 0.9)

    def test_bad_track_misses_threshold_without_posthoc_time_shift(self):
        def change(rows):
            for row in rows[2:]:
                row["radius_m"] = str(float(row["radius_m"]) + 10 * 695700000.0)
        self.rewrite("observations.csv", change, "observations")
        result = self.evaluate()
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("shock-track residual threshold missed", result["violations"])

    def test_no_holdout_and_front_substitution_rejected(self):
        self.manifest(fit_used_holdout=True)
        with self.assertRaisesRegex(runner.EvidenceError, "used to fit"):
            self.evaluate()
        self.manifest(fit_used_holdout=False, feature="ejecta-front")
        with self.assertRaisesRegex(runner.EvidenceError, "ejecta tracks"):
            self.evaluate()

    def test_overlapping_fit_and_observation_intervals_rejected(self):
        self.rewrite("observations.csv", lambda rows: rows[1].update(role="validation"), "observations")
        with self.assertRaisesRegex(runner.EvidenceError, "overlap"):
            self.evaluate()

    def test_hash_and_clock_corruption_rejected(self):
        with (self.bundle / "history.csv").open("a") as stream:
            stream.write("corrupt\n")
        with self.assertRaisesRegex(runner.EvidenceError, "SHA256 mismatch"):
            self.evaluate()
        runner.write_demo(self.bundle)
        self.rewrite("history.csv", lambda rows: rows[10].update(tick="999"), "model")
        with self.assertRaisesRegex(runner.EvidenceError, "integer clock"):
            self.evaluate()

    def test_missing_evidence_skip(self):
        self.assertEqual(runner.evaluate(self.root / "missing", self.output)["status"], "SKIP")

    def test_cli_requires_real_evidence_and_does_not_reuse_outputs(self):
        command = [sys.executable, str(ROOT / "validation/run_swcme_coupling_validation.py"),
                   "--bundle", str(self.root / "missing"), "--output-dir", str(self.output),
                   "--require-evidence"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        record = json.loads((self.output / "result.json").read_text())
        self.assertEqual(record["status"], "SKIP")
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("must be empty", result.stderr)

    def test_real_renderer_produces_png_and_vector_eps(self):
        result = runner.evaluate(self.bundle, self.output)
        png = self.output / "swcme-coupling-comparison.png"
        eps = self.output / "swcme-coupling-comparison.eps"
        self.assertEqual(png.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
        self.assertTrue(eps.read_bytes().startswith(b"%!PS-Adobe-3.0 EPSF-3.0"))
        self.assertIn("SYNTHETIC", (self.output / "figure-caption.txt").read_text())
        self.assertEqual(result["status"], "SKIP")

    def make_reference(self):
        # Tiny parser fixtures, not an observed event or native validation run.
        columns={
            "wind_protons.csv": "time_utc,proton_bulk_speed_m_s,proton_density_m3,proton_temperature_k,quality\n2000-01-01T00:00:00Z,400000,5000000,100000,good\n",
            "wind_magnetic_1min.csv": "time_utc,magnitude_t,bx_gse_t,by_gse_t,bz_gse_t,quality\n2000-01-01T00:00:00Z,5e-9,3e-9,4e-9,0,good\n",
            "wind_magnetic_3sec_shock.csv": "time_utc,magnitude_t,bx_gse_t,by_gse_t,bz_gse_t,quality\n2000-01-01T00:00:00Z,5e-9,3e-9,4e-9,0,good\n",
            "wind_orbit.csv": "time_utc,x_hec_m,y_hec_m,z_hec_m,heliocentric_radius_m\n2000-01-01T00:00:00Z,149597870700,0,0,149597870700\n",
            "helcats_time_elongation.csv": "time_utc,event_id,trace_id,elongation_deg,position_angle_deg,feature,role\n2000-01-01T00:00:00Z,fixture,0,10,75,unclassified CME brightness,diagnostic\n",
        }
        for name,text in columns.items():(self.bundle/name).write_text(text)
        original=runner.read_json(self.bundle/"manifest.json")
        arrival=dict(original["arrival"],spacecraft="Wind")
        runner.atomic_json(self.bundle/"arrival.json",arrival)
        reference=dict(schema=runner.REFERENCE_SCHEMA,case_id="CME3D02",reference_kind="observations",
                       event_name="Parser fixture",comparison_scope="arrival-only",arrival=arrival,
                       native_manifest="native-manifest.json",
                       files=[dict(file=name,sha256=runner.sha256(self.bundle/name)) for name in [*columns,"arrival.json"]])
        runner.atomic_json(self.bundle/"manifest.json",reference)
        return original

    def test_reference_only_data_are_found_without_inventing_native_evidence(self):
        self.make_reference()
        result=self.evaluate()
        self.assertEqual(result["status"],"SKIP")
        self.assertTrue(result["reference_ready"])
        self.assertIn("native-manifest.json",result["message"])
        self.assertEqual(result["metrics"]["reference_counts"]["wind_protons.csv"],1)
        self.assertFalse(result["shock_radius_track_validated"])

    def test_reference_corruption_errors_before_native_skip(self):
        self.make_reference()
        with (self.bundle/"wind_protons.csv").open("a") as stream:stream.write("corruption\n")
        with self.assertRaisesRegex(runner.EvidenceError,"SHA256 mismatch"):self.evaluate()

    def test_attached_native_history_uses_arrival_diagnostic_not_unclassified_tracks(self):
        original=self.make_reference()
        # Retain the explicitly synthetic producer in the native fixture, and
        # prove it cannot pass as an installed native AMPS history.
        runner.atomic_json(self.bundle/"native-manifest.json",original)
        with self.assertRaisesRegex(runner.EvidenceError,"native producer"):self.evaluate()

    def test_hapi_fill_masking_and_si_conversion(self):
        metadata={"density":dict(fill="99999.9"),"field":dict(fill="-1e31")}
        self.assertEqual(preparer.value(dict(density=5),metadata,"density",1e6),5e6)
        self.assertIsNone(preparer.value(dict(density=99999.9),metadata,"density",1e6))
        converted=preparer.value(dict(field=[3,-1e31,4]),metadata,"field",1e-9)
        self.assertAlmostEqual(converted[0],3e-9,places=16)
        self.assertIsNone(converted[1])
        self.assertAlmostEqual(converted[2],4e-9,places=16)

    def test_arrival_only_renderer_keeps_synthetic_fixture_unqualified(self):
        self.manifest(comparison_scope="arrival-only")
        result=runner.evaluate(self.bundle,self.output)
        self.assertEqual(result["status"],"SKIP")
        self.assertFalse(result["shock_radius_track_validated"])
        self.assertEqual(result["metrics"]["holdout_points"],0)
        self.assertIn("SYNTHETIC",(self.output/"figure-caption.txt").read_text())
        self.assertTrue((self.output/"swcme-coupling-comparison.eps").read_bytes().endswith(b"showpage\n"))

    def test_empty_and_nonmonotonic_hapi_data_are_errors(self):
        path=self.root/"hapi.json"
        path.write_text(json.dumps(dict(status=dict(code=1200),parameters=[dict(name="Time")],data=[])))
        with self.assertRaisesRegex(ValueError,"successful data response"):preparer.hapi_records(path)
        payload=dict(status=dict(code=1200),parameters=[dict(name="Time")],data=[["2000-01-02T00:00:00Z"],["2000-01-01T00:00:00Z"]])
        path.write_text(json.dumps(payload))
        with self.assertRaisesRegex(ValueError,"nonmonotonic"):preparer.hapi_records(path)


class NativeOrchestrationTests(unittest.TestCase):
    """Software protocol fixtures; these do not simulate MPI or qualify AMPS."""
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name); self.bundle=self.root/"fixture"
        runner.write_demo(self.bundle)
        self.model=runner.csv_rows(self.bundle/"history.csv",runner.MODEL_COLUMNS)
        self.runtime=dict(schema="srcsep3d-native-shock-runtime-v1",producer="srcSEP3D-native",
            run_intent="shock-propagation",source_enabled=False,time_step_s=60,
            maximum_time_steps=len(self.model)-1,completed_steps=len(self.model)-1,final_time_s=float(self.model[-1]["time_s"]),
            final_radius_m=float(self.model[-1]["shock_radius_m"]),stop_shock_radius_m=0,
            stop_reason="step-budget",history_file="shock-history.csv",mpi_ranks=4,
            provider_identity="canonical-swcme3d-standalone",provider_configuration_fingerprint="fixture-only",
            application_configuration_fingerprint="fixture-fingerprint")

    def test_native_completion_and_missing_ticks_are_checked(self):
        rows=runner.read_native_history(self.bundle/"history.csv",self.runtime)
        self.assertEqual(len(rows),4001)
        self.runtime["completed_steps"]-=1
        with self.assertRaisesRegex(runner.EvidenceError,"completion count"):
            runner.read_native_history(self.bundle/"history.csv",self.runtime)
        self.runtime["completed_steps"]+=1
        self.model[20]["tick"]=21
        with (self.bundle/"history.csv").open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=runner.MODEL_COLUMNS);writer.writeheader();writer.writerows(self.model)
        with self.assertRaisesRegex(runner.EvidenceError,"missing ticks"):
            runner.read_native_history(self.bundle/"history.csv",self.runtime)

    def test_radius_stop_must_have_a_real_crossing_bracket(self):
        self.runtime.update(stop_reason="shock-radius",stop_shock_radius_m=runner.AU_M)
        with self.assertRaisesRegex(runner.EvidenceError,"bracket"):
            runner.read_native_history(self.bundle/"history.csv",self.runtime)

    def launch_fixture(self, **changes):
        # Mock only the OS launch. The actual orchestration freezes the shipped
        # deck, preserves non-UTF8 raw diagnostics, checks reports and renders.
        executable=self.root/"fixture-executable";executable.write_text("TEST FIXTURE ONLY\n");executable.chmod(0o755)
        output=self.root/"output";output.mkdir()
        args=SimpleNamespace(amps=executable,input=ROOT/"examples/sep3d_swcme_20rs_1au.in",
             ranks=4,launcher="fixture-launch -n {ranks}",event_fit=None,bundle=self.root/"missing-reference")
        for key,value in changes.items():setattr(args,key,value)
        def launch(command,**kwargs):
            native=Path(command[command.index("--output-dir")+1])
            (native/"shock-history.csv").write_bytes((self.bundle/"history.csv").read_bytes())
            runner.atomic_json(native/"native-runtime.json",self.runtime)
            return SimpleNamespace(stdout=io.BytesIO(b"fixture log \xc1 preserved\n"),wait=lambda:0)
        completed=SimpleNamespace(returncode=0,stdout=b"run_intent=shock-propagation\nsource_enabled=false\nphysics_fingerprint=fixture-fingerprint\n")
        with patch.object(runner.subprocess,"run",return_value=completed),patch.object(runner.subprocess,"Popen",side_effect=launch):
            result=runner.native_campaign(args,output)
        return result,output

    def test_native_control_stays_unqualified_and_keeps_provenance(self):
        result,output=self.launch_fixture()
        self.assertEqual(result["status"],"SKIP")
        self.assertTrue(result["native_mechanics_pass"])
        self.assertFalse(result["production_release_qualified"])
        launch=runner.read_json(output/"launch.json")
        self.assertEqual(launch["mpi_ranks"],4)
        self.assertEqual(launch["input_sha256"],runner.sha256(output/"native/input.in"))
        self.assertIn(b"\xc1",(output/"native.log").read_bytes())
        self.assertTrue((output/"native-propagation-control.png").is_file())
        self.assertIn(str(output/"native"),(output/"native/input.in").read_text())

    def test_native_source_violation_fails_before_observations(self):
        self.model[0]["particle_count"]=1
        with (self.bundle/"history.csv").open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=runner.MODEL_COLUMNS);writer.writeheader();writer.writerows(self.model)
        result,_=self.launch_fixture()
        self.assertEqual(result["status"],"FAIL")
        self.assertIn("nonzero",result["message"])

    def test_event_fit_must_match_original_input(self):
        fit=self.root/"event-fit.json"
        runner.atomic_json(fit,dict(schema="srcsep3d-swcme-event-fit-v1",fit_used_holdout=False,input_sha256="0"*64))
        with self.assertRaisesRegex(runner.EvidenceError,"original input"):
            self.launch_fixture(event_fit=fit)


if __name__ == "__main__":
    unittest.main(verbosity=2)

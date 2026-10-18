#!/usr/bin/env python3
"""Deterministic Stage-9 UI/compiler boundary and workflow tests.

These tests exercise the interface contract without spending Geant4 histories.
When a shieldSim executable is supplied, every generated case is additionally
passed through the production C++ parser in ``--validate-case`` mode.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time
import unittest


if len(sys.argv) < 3:
    raise SystemExit("usage: test_stage9_interface.py SOURCE_ROOT OUTPUT_DIR [SHIELDSIM]")
ROOT = pathlib.Path(sys.argv.pop(1)).resolve()
OUTPUT = pathlib.Path(sys.argv.pop(1)).resolve()
SHIELDSIM = pathlib.Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else None
# Test execution must not create source-tree artifacts.
sys.dont_write_bytecode = True
SPEC = importlib.util.spec_from_file_location("stage9_interface",
                                              ROOT / "scripts/stage9_interface.py")
stage9 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(stage9)


def load_example() -> dict:
    return json.loads((ROOT / "examples/stage9_ui_request.json").read_text())


FAKE_RUNNER = """#!/usr/bin/env python3
import json,pathlib,re,sys
case=pathlib.Path(sys.argv[1].split('=',1)[1])
definition=json.loads(case.read_text())
case_hash='synthetic-test-case-hash'
for absorber in definition['geometry']['absorbers']:
    safe=re.sub(r'[^A-Za-z0-9_-]','_',absorber['absorber_id'])
    prefix=definition['output']['prefix']+'__'+safe
    pathlib.Path(prefix+'_resolved_case.json').write_text(json.dumps({
        'case_id':definition['case_id'],'active_absorber_id':absorber['absorber_id'],
        'case_hash':case_hash}))
    pathlib.Path(prefix+'_result.json').write_text(json.dumps({
        'absorber_id':absorber['absorber_id'],'case_hash':case_hash}))
"""


class Stage9InterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        OUTPUT.mkdir(parents=True, exist_ok=True)

    def test_matrix_traceability_round_trip_and_reproducibility(self) -> None:
        request = load_example()
        first = stage9.compile_request(request)
        second = stage9.compile_request(json.loads(json.dumps(request)))
        manifest = first["manifest"]
        self.assertEqual(manifest["case_matrix"]["shield_count"], 2)
        self.assertEqual(manifest["case_matrix"]["absorber_count"], 2)
        self.assertEqual(manifest["case_matrix"]["transport_case_count"], 4)
        self.assertEqual(manifest["cost_preflight"]["requested_quantity_series_count"], 4)
        self.assertEqual(len(manifest["case_matrix"]["mappings"]), 4)
        self.assertEqual(manifest["ui_request"], request)
        self.assertEqual(manifest["manifest_content_hash"],
                         second["manifest"]["manifest_content_hash"])
        self.assertEqual([x["case_hash"] for x in manifest["generated_cases"]],
                         [x["case_hash"] for x in second["manifest"]["generated_cases"]])
        for mapping in manifest["case_matrix"]["mappings"]:
            self.assertIn(mapping["shield_id"], {"al_hdpe", "al_only"})
            self.assertIn(mapping["absorber_id"], {"silicon_proxy", "water_proxy"})
            self.assertEqual(mapping["case_json_paths"]["quantities"], "$.quantities")

    def test_order_density_and_thickness_conversion(self) -> None:
        manifest = stage9.compile_request(load_example())["manifest"]
        layers = manifest["ordered_shield_layers"][0]["layers"]
        self.assertEqual([layer["source_to_absorber_order"] for layer in layers], [1, 2])
        self.assertAlmostEqual(layers[0]["density_g_cm3"], 2.699)
        self.assertAlmostEqual(layers[0]["mass_thickness_g_cm2"], 0.2699)
        self.assertAlmostEqual(layers[1]["geometric_thickness_mm"], 2.0)
        self.assertIn("density_uncertainty", layers[0])

        # Case-local material properties remain changeable through the input.
        request = load_example()
        request["material_overrides"] = [{
            "material_id": "Al_project_density", "revision": "ui-test-1",
            "extends": "Al_natural_G4_11_4@1", "status": "project_specific",
            "density": {"value": 2.5, "unit": "g/cm3",
                        "uncertainty": {"type": "standard", "value": 0.01,
                                        "unit": "g/cm3"},
                        "reference_condition": "synthetic UI contract test"},
            "provenance": [{"properties": ["density"], "source": "test fixture",
                            "version": "1"}]
        }]
        request["shield_configurations"][0]["layers"][0]["material"] = "Al_project_density"
        resolved = stage9.compile_request(request)["manifest"]["ordered_shield_layers"][0]["layers"][0]
        self.assertEqual(resolved["density_g_cm3"], 2.5)
        self.assertAlmostEqual(resolved["mass_thickness_g_cm2"], 0.25)

    def assertRejected(self, mutate, phrase: str) -> None:  # noqa: N802
        request = load_example()
        mutate(request)
        with self.assertRaisesRegex(stage9.Stage9Error, phrase):
            stage9.compile_request(request)

    def test_scientific_and_geometry_guards(self) -> None:
        def make_six_layers(request: dict) -> None:
            prototype = request["shield_configurations"][0]["layers"][0]
            request["shield_configurations"][0]["layers"] = []
            for index in range(6):
                layer = copy.deepcopy(prototype)
                layer["layer_id"] = f"layer_{index}"
                request["shield_configurations"][0]["layers"].append(layer)

        self.assertRejected(lambda r: r["quantities"].append("DDD"), "cannot combine")
        self.assertRejected(
            lambda r: r.update(purpose="spacecraft_anomaly_environment_effect_assessment"),
            "disables purpose")
        self.assertRejected(make_six_layers, "one through five")
        self.assertRejected(
            lambda r: r["shield_configurations"][0]["layers"][0].update(
                thickness={"type": "geometric", "value": 101, "unit": "mm"}),
            "exceeds the approved 100 mm")
        self.assertRejected(
            lambda r: r["absorbers"][0]["geometry"].pop("size_x"), "size_x")
        self.assertRejected(
            lambda r: r["absorbers"][0]["geometry"].update(type="mesh"),
            "arbitrary shapes remain Stage 10")
        self.assertRejected(
            lambda r: r["shield_configurations"][0]["layers"][0].update(
                material="LunarRegolith"), "disabled template")

    def test_cost_and_uncertainty_preflight(self) -> None:
        request = load_example()
        manifest = stage9.compile_request(request)["manifest"]
        self.assertEqual(manifest["cost_preflight"]["wall_time_status"],
                         "not_assessable_no_measured_throughput")
        self.assertEqual(manifest["uncertainty_preflight"][0]["status"],
                         "assessable_at_runtime_with_remaining_statistical_budget")
        request["runtime_calibration"] = {
            "status": "measured", "seconds_per_history": 0.25,
            "provenance": "test host timing fixture"
        }
        cost = stage9.compile_request(request)["manifest"]["cost_preflight"]
        self.assertEqual(cost["minimum_wall_seconds"], 200.0)
        self.assertEqual(cost["maximum_wall_seconds"], 2000.0)

    def test_every_exposed_profile_compiles(self) -> None:
        selections = {
            "tid_slab_verification_v1": ["TID"],
            "silicon_damage_slab_verification_v1": ["DDD", "n_eq"],
            "let_hardness_slab_verification_v1": ["LET", "H100/10"],
        }
        for profile, quantities in selections.items():
            with self.subTest(profile=profile):
                request = load_example()
                request["request_id"] = "profile_" + profile
                request["profile_id"] = profile
                request["quantities"] = quantities
                # Silicon is inside all three bounded response domains.  One
                # absorber keeps this a profile-contract check rather than a
                # repetition of the matrix-cardinality test.
                request["absorbers"] = request["absorbers"][:1]
                compiled = stage9.compile_request(request)
                self.assertEqual(compiled["manifest"]["quantities"], quantities)
                self.assertEqual({badge["quantity"] for badge in
                                  compiled["manifest"]["validation_badges"]},
                                 set(quantities))
                if SHIELDSIM is not None:
                    bundle = OUTPUT / ("profile_" + profile)
                    if bundle.exists():
                        shutil.rmtree(bundle)
                    stage9.write_bundle(compiled, bundle)
                    for case in compiled["manifest"]["generated_cases"]:
                        validation = bundle / "validate" / case["case_id"]
                        validation.mkdir(parents=True, exist_ok=True)
                        result = subprocess.run(
                            [str(SHIELDSIM), "--validate-case=" +
                             str(bundle / case["relative_path"])],
                            cwd=validation, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, check=False)
                        self.assertEqual(result.returncode, 0, result.stdout)

    def test_candidate_regression_snapshot(self) -> None:
        baseline = json.loads((ROOT / "tests/reference/stage9_candidate_baseline.v1.json").read_text())
        manifest = stage9.compile_request(load_example())["manifest"]
        self.assertEqual(baseline["status"], "candidate_unreviewed_not_acceptance")
        self.assertEqual(baseline["portable_regression_contract_sha256"],
                         manifest["reproducibility"]["portable_regression_contract_sha256"])
        self.assertEqual(baseline["expected_case_matrix"], {
            "shield_count": 2, "absorber_count": 2, "transport_case_count": 4})

    def test_bundle_cpp_validation_and_cancel_resume(self) -> None:
        work = OUTPUT / "bundle_workflow"
        if work.exists():
            shutil.rmtree(work)
        compiled = stage9.compile_request(load_example())
        stage9.write_bundle(compiled, work)

        if SHIELDSIM is not None:
            for case in compiled["manifest"]["generated_cases"]:
                validation = OUTPUT / (case["case_id"] + "_validation")
                validation.mkdir(parents=True, exist_ok=True)
                result = subprocess.run(
                    [str(SHIELDSIM), "--validate-case=" + str(work / case["relative_path"])],
                    cwd=validation, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, check=False)
                self.assertEqual(result.returncode, 0, result.stdout)

        fake = OUTPUT / "fake_shieldsim.py"
        fake.write_text(FAKE_RUNNER, encoding="utf-8")
        fake.chmod(0o755)
        canceled = stage9.run_bundle(work, fake, max_new_cases=1)
        self.assertEqual(canceled["status"], "canceled")
        self.assertEqual(len(canceled["completed_case_ids"]), 1)
        resumed = stage9.run_bundle(work, fake, resume=True)
        self.assertEqual(resumed["status"], "complete")
        self.assertEqual(len(resumed["completed_case_ids"]), 2)
        for case in compiled["manifest"]["generated_cases"]:
            provenance = json.loads((work / "runs" / case["case_id"] /
                                     "stage9_provenance.json").read_text())
            self.assertEqual(provenance["manifest_content_hash"],
                             compiled["manifest"]["manifest_content_hash"])
            self.assertTrue(provenance["validation_badges"])

    def test_bundle_tampering_is_rejected_before_transport(self) -> None:
        work = OUTPUT / "tamper_guard"
        if work.exists():
            shutil.rmtree(work)
        compiled = stage9.compile_request(load_example())
        stage9.write_bundle(compiled, work)
        case_path = work / compiled["manifest"]["generated_cases"][0]["relative_path"]
        case = json.loads(case_path.read_text())
        case["accuracy"]["target"] = 0.99
        case_path.write_text(json.dumps(case), encoding="utf-8")
        fake = OUTPUT / "must_not_run.py"
        fake.write_text("#!/usr/bin/env python3\nraise RuntimeError('must not run')\n",
                        encoding="utf-8")
        fake.chmod(0o755)
        with self.assertRaisesRegex(stage9.Stage9Error, "changed"):
            stage9.run_bundle(work, fake)
        self.assertFalse((work / "job_state.json").exists())

    def test_active_case_cancellation_restarts_partial_case(self) -> None:
        work = OUTPUT / "active_cancel_workflow"
        if work.exists():
            shutil.rmtree(work)
        stage9.write_bundle(stage9.compile_request(load_example()), work)
        fake = OUTPUT / "slow_fake_shieldsim.py"
        fake.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib,time\n"
            "pathlib.Path('started').write_text('active')\n"
            "time.sleep(30)\n",
            encoding="utf-8")
        fake.chmod(0o755)
        runner = subprocess.Popen(
            [sys.executable, str(ROOT / "scripts/stage9_interface.py"), "run",
             "--bundle", str(work), "--executable", str(fake)],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            state_path = work / "job_state.json"
            if state_path.exists() and json.loads(state_path.read_text()).get("current_case_id"):
                break
            time.sleep(0.05)
        else:
            runner.kill()
            self.fail("runner did not enter an active case")
        stage9.request_cancel(work)
        output, _ = runner.communicate(timeout=10)
        self.assertEqual(runner.returncode, 0, output)
        canceled = json.loads((work / "job_state.json").read_text())
        self.assertEqual(canceled["status"], "canceled")
        self.assertEqual(canceled["completed_case_ids"], [])
        current_attempts = next(iter(canceled["case_attempts"].values()))
        self.assertEqual(current_attempts[-1]["status"], "interrupted_restart_required")

        # Replace only the test double, then resume.  The interrupted case must
        # execute again because within-case Geant4 histories are never reused.
        fake.write_text(FAKE_RUNNER, encoding="utf-8")
        fake.chmod(0o755)
        resumed = stage9.run_bundle(work, fake, resume=True)
        self.assertEqual(resumed["status"], "complete")
        first_case_attempts = next(iter(resumed["case_attempts"].values()))
        self.assertEqual(len(first_case_attempts), 2)
        self.assertEqual(first_case_attempts[-1]["status"], "complete")

    def test_browser_assets_expose_required_controls(self) -> None:
        html = (ROOT / "ui/stage9/index.html").read_text()
        javascript = (ROOT / "ui/stage9/app.js").read_text()
        for control in ("add-shield", "add-absorber", "preview", "compile",
                        "run", "cancel", "resume", "status"):
            self.assertIn(f'id="{control}"', html)
        self.assertIn("source_to_absorber", (ROOT / "scripts/stage9_interface.py").read_text())
        self.assertIn("/api/preview", javascript)
        self.assertIn("transport_case_count", javascript)
        node = shutil.which("node")
        if node:
            result = subprocess.run([node, "--check", str(ROOT / "ui/stage9/app.js")],
                                    capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            # The CI host may provide Chromium rather than Node. Loading the
            # page from file still parses and starts app.js; API fetches are
            # expected to fail because no HTTP server is opened in this test.
            chrome = shutil.which("google-chrome") or shutil.which("chromium")
            if chrome:
                chrome_profile = OUTPUT / "chrome_profile"
                if chrome_profile.exists():
                    shutil.rmtree(chrome_profile)
                result = subprocess.run(
                    [chrome, "--headless=new", "--no-sandbox", "--disable-gpu",
                     "--user-data-dir=" + str(chrome_profile),
                     "--dump-dom", (ROOT / "ui/stage9/index.html").as_uri()],
                    capture_output=True, text=True, check=False)
                if result.returncode != 0 and ("Operation not permitted" in result.stderr or
                                                "unique user data directory" in result.stderr):
                    return  # sandbox forbids Chrome's helper sockets/directories
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertNotIn("SyntaxError", result.stdout + result.stderr)
                self.assertIn("Initialization failed: Failed to fetch", result.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)

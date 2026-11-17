#!/usr/bin/env python3
"""I17 release/scorer protocol tests; temporary values are NOT observations.

All CSV/JSON inputs below live in TemporaryDirectory and carry protocol-only
identity. They test refusal, hashing, selection and arithmetic. No production
simulation or canonical observational product is opened, scored or qualified.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

TEST_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TEST_ROOT))
import i17_holdout as gate


class HoldoutContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="i17-protocol-only-")
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)
        self.native = {"protocol-only.cpp": "0"*64}
        self.patch = patch("build_m2_qualification.production_source_identity",
                           return_value=self.native)
        self.patch.start()
        self.addCleanup(self.patch.stop)
        self.git_patch = patch("i17_holdout.reviewed_git_identity", return_value={
            "source_revision": "protocol-fixture-only", "blob_oid": "a"*40})
        self.git_patch.start()
        self.addCleanup(self.git_patch.stop)
        self.freeze = self.root/"reviewed.json"
        self.campaign = self.root/"approved-protocol-campaign"
        self.exe = self.write("executable", "PROTOCOL ONLY, NEVER EXECUTED\n")
        self.history = self.write("history.csv", "selection\ntraining\nholdout\n")
        self.contract = self.write("contract.json", '{"protocol_only":true}\n')
        self.processor = self.write("processor.py", "# protocol-only bytes\n")
        common = ("#begin application: moon\n"
                  "surface_geometry = lola\n"
                  "spice_path = /protocol-only/spice\n"
                  "lola_image_file = /protocol-only/lola.img\n"
                  "lola_label_file = /protocol-only/lola.lbl\n"
                  "surface_cea_file = /protocol-only/surface.cea\n"
                  "surface_tecplot_file = /protocol-only/surface.dat\n"
                  "surface_inventory_checkpoint_file = none\n"
                  "diviner_eclipse_history_file = none\n"
                  "diviner_eclipse_selection_contract_file = none\n"
                  "thermal_infrared_emissivity = 0.95\n"
                  "thermal_time_step_s = 300\n")
        self.training_input = self.write("training.in", common+
            "simulation_start_utc = 2011-06-15T16:50:00Z\n"
            "thermal_spinup_cycle_start_utc = 2011-05-17T04:05:57.024Z\n"
            "thermal_spinup_cycle_end_utc = 2011-06-15T16:50:00Z\n#end application\n")
        self.template = self.write("holdout.in", common+
            "simulation_start_utc = 2011-12-10T11:00:00Z\n"
            "thermal_spinup_cycle_start_utc = 2011-11-10T00:00:00Z\n"
            "thermal_spinup_cycle_end_utc = 2011-12-10T11:00:00Z\n"
            "diviner_eclipse_selection = holdout\n"
            f"diviner_eclipse_freeze_file = {self.freeze}\n#end application\n")
        candidate = self.write_json("candidate.json", {
            "status": "REVIEW_REQUIRED_NOT_FROZEN", "selection": "training", "comparison_count": 1})
        result = self.write_json("training_result.json", {"status": "SKIPPED", "metrics": {
            "thermal_convergence": {"failed_bins": 0, "compared_bins": 1,
                                    "maximum_update_changes_K": [1, .5, .25]}}})
        manifest = self.write_json("training_manifest.json", {"executable": {
            "sha256": gate.sha256(self.exe)}})
        model = self.write("training_model.csv", "selection,site_time_bin_id\ntraining,protocol-only\n")
        convergence = self.write("convergence.csv", "protocol_fixture_only\n")
        scope = self.write("scope.txt", "PROTOCOL ONLY: NOT MODEL QUALIFICATION\n")
        self.doc = {"schema": "srcMoon-I17-reviewed-holdout-v1",
            "status": "APPROVED_FROZEN", "approved_by": "protocol fixture reviewer",
            "approved_utc": "2026-11-01T00:00:00Z", "frozen_before_holdout_residuals": True,
            "holdout_campaign_id": self.campaign.name,
            "initialization": "independent-observation-free-periodic-spinup",
            "processed_history_sha256": gate.sha256(self.history),
            "selection_contract_sha256": gate.sha256(self.contract),
            "processor_sha256": gate.sha256(self.processor),
            "executable_sha256": gate.sha256(self.exe),
            "production_source_identity_sha256": hashlib.sha256(json.dumps(
                self.native, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "reviewed_holdout_input": self.record(self.template),
            "training_evidence": {key: self.record(path) for key, path in (
                ("input", self.training_input), ("candidate", candidate),
                ("result", result), ("manifest", manifest),
                ("model_csv", model), ("convergence_csv", convergence))},
            "acceptance": {"weighted_rmse_k_max": 0,
                "absolute_weighted_bias_k_max": 0, "normalized_residual_rms_max": 0,
                "coverage_sigma": 1, "minimum_coverage_fraction": 1,
                "phase_amplitude_operator": "sampled-extrema-on-frozen-bin-groups",
                "sampled_phase_error_s_max": 0, "sampled_amplitude_error_k_max": 0},
            "phase_amplitude_groups": [{"id": "protocol-only-site",
                "site_time_bin_ids": ["a0", "a1", "a2"]}],
            "model_scope_qualification": {"status": "READY", "evidence": [self.record(scope)]}}
        self.save()

    def write(self, name, text):
        path = self.root/name
        path.write_text(text, encoding="utf-8")
        return path

    def write_json(self, name, document):
        return self.write(name, json.dumps(document, indent=2)+"\n")

    def record(self, path):
        return {"path": path.name, "sha256": gate.sha256(path)}

    def save(self):
        self.freeze.write_text(json.dumps(self.doc, indent=2)+"\n")

    def release(self):
        return gate.validate_release(self.template, self.exe, history=self.history,
            contract=self.contract, processor=self.processor, campaign=self.campaign)

    def model_rows(self):
        # Exact arithmetic fixtures test the gate; they are not synthetic
        # observations, are never plotted and cannot enter a formal campaign.
        return [{"selection": "holdout", "site_time_bin_id": f"a{i}",
            "observed_bolometric_brightness_temperature_K": str(value),
            "model_bolometric_brightness_temperature_K": str(value),
            "observed_uncertainty_K": "2", "latitude_deg": "80",
            "longitude_east_deg": "0", "ephemeris_time_tdb_s": str(60*i)}
            for i, value in enumerate((100, 80, 95))]

    def test_valid_review_is_read_only(self):
        before = {p: gate.sha256(p) for p in self.root.iterdir() if p.is_file()}
        release = self.release()
        self.assertEqual(release["sha256"], gate.sha256(self.freeze))
        self.assertEqual(before, {p: gate.sha256(p) for p in before})
        self.assertFalse(self.campaign.exists())

    def test_legacy_default_and_unknown_selection(self):
        self.assertEqual(gate.selection(self.training_input), "training")
        self.template.write_text(self.template.read_text().replace("selection = holdout", "selection = other"))
        with self.assertRaises(ValueError): gate.selection(self.template)

    def test_candidate_never_authorizes_holdout(self):
        self.doc["status"] = "REVIEW_REQUIRED_NOT_FROZEN"
        self.save()
        with self.assertRaises(gate.HoldoutUnavailable): self.release()

    def test_missing_review_is_unavailable(self):
        self.freeze.unlink()
        with self.assertRaises(gate.HoldoutUnavailable): self.release()

    def test_named_utc_preholdout_review_required(self):
        for key, value in (("approved_by", ""), ("approved_utc", "2026-11-01"),
                           ("frozen_before_holdout_residuals", False)):
            old = self.doc[key]
            self.doc[key] = value; self.save()
            with self.subTest(key=key), self.assertRaises(ValueError): self.release()
            self.doc[key] = old

    def test_changed_science_or_binary_bytes_are_errors(self):
        for key in ("processed_history_sha256", "selection_contract_sha256",
                    "processor_sha256", "executable_sha256", "production_source_identity_sha256"):
            old = self.doc[key]; self.doc[key] = "1"*64; self.save()
            with self.subTest(key=key), self.assertRaises(ValueError): self.release()
            self.doc[key] = old

    def test_changed_training_record_is_error(self):
        self.doc["training_evidence"]["model_csv"]["sha256"] = "1"*64; self.save()
        with self.assertRaises(ValueError): self.release()

    def test_missing_training_record_is_unavailable(self):
        (self.root/"training_model.csv").unlink()
        with self.assertRaises(gate.HoldoutUnavailable): self.release()

    def test_material_change_cannot_be_hidden_by_new_holdout_input_hash(self):
        self.template.write_text(self.template.read_text().replace("emissivity = 0.95", "emissivity = 0.90"))
        self.doc["reviewed_holdout_input"] = self.record(self.template); self.save()
        with self.assertRaises(ValueError): self.release()

    def test_june_spinup_cannot_be_reused_for_december(self):
        self.template.write_text(self.template.read_text().replace("2011-12-10T11:00:00Z", "2011-06-15T16:50:00Z"))
        self.doc["reviewed_holdout_input"] = self.record(self.template); self.save()
        with self.assertRaises(ValueError): self.release()

    def test_incomplete_numerical_gate_cannot_be_reviewed_as_pass(self):
        path = self.root/"training_result.json"
        path.write_text('{"metrics":{"thermal_convergence":{"failed_bins":1}}}')
        self.doc["training_evidence"]["result"] = self.record(path); self.save()
        with self.assertRaises(ValueError): self.release()

    def test_phase_or_scope_gap_remains_unavailable(self):
        for key, value in (("phase_amplitude_groups", []),
                           ("model_scope_qualification", {"status": "PARTIAL"})):
            old = self.doc[key]; self.doc[key] = value; self.save()
            with self.subTest(key=key), self.assertRaises(gate.HoldoutUnavailable): self.release()
            self.doc[key] = old

    def test_invalid_limits_are_errors(self):
        for value in (None, True, -1, float("nan"), float("inf")):
            self.doc["acceptance"]["weighted_rmse_k_max"] = value; self.save()
            with self.subTest(value=value), self.assertRaises(ValueError): self.release()

    def test_approval_is_bound_to_one_campaign(self):
        self.doc["holdout_campaign_id"] = "another-campaign"; self.save()
        with self.assertRaises(ValueError): self.release()

    def test_repeated_claim_requires_review(self):
        self.campaign.mkdir()
        release = self.release()
        gate.claim_release(self.campaign, release)
        with self.assertRaises(gate.HoldoutUnavailable): gate.claim_release(self.campaign, release)

    def test_exact_arithmetic_contract(self):
        status, metrics, reason = gate.score_holdout(self.model_rows(), self.release())
        self.assertEqual(status, "PASS")
        self.assertTrue(all(metrics["acceptance_checks"].values()))

    def test_resolved_amplitude_and_phase_misses_fail(self):
        for predictions, failed in (((100, 80, 105), "sampled_amplitude"),
                                    ((80, 100, 95), "sampled_phase")):
            rows = self.model_rows()
            for row, value in zip(rows, predictions):
                row["model_bolometric_brightness_temperature_K"] = str(value)
            status, metrics, _ = gate.score_holdout(rows, self.release())
            self.assertEqual(status, "FAIL")
            self.assertFalse(metrics["acceptance_checks"][failed])

    def test_ambiguous_observed_phase_is_not_fabricated(self):
        rows = self.model_rows()
        rows[2]["observed_bolometric_brightness_temperature_K"] = "80"
        self.assertEqual(gate.score_holdout(rows, self.release())[0], "SKIPPED")

    def test_missing_or_changed_group_support_is_error(self):
        rows = self.model_rows()
        with self.assertRaises(ValueError): gate.score_holdout(rows[:-1], self.release())
        rows[1]["latitude_deg"] = "81"
        with self.assertRaises(ValueError): gate.score_holdout(rows, self.release())

    def test_unselected_rows_cannot_be_scored(self):
        rows = self.model_rows(); rows[0]["selection"] = "training"
        with self.assertRaises(ValueError): gate.score_holdout(rows, self.release())

    def canonical_protocol_package(self):
        """Create metadata/filter fixtures, with NO temperature observations."""
        import integration_common as integration
        package = self.root/"canonical-protocol-only/diviner/processed/i17"
        package.mkdir(parents=True)
        history = package/"eclipse_thermal_history.csv"
        history.write_text("selection,protocol_fixture\ntraining,june-filter-token\nholdout,december-filter-token\n")
        contract = package/"selection_contract.json"
        contract.write_bytes((TEST_ROOT/"integration-tests/I17/reference/selection_contract.json").read_bytes())
        processor = TEST_ROOT/"data_processing/prepare_d03_rdr_eclipse.py"
        (package/"d03_rdr_qa.json").write_text(json.dumps({"status": "READY",
            "processed": {"sha256": gate.sha256(history)},
            "contract": {"sha256": gate.sha256(contract)},
            "processor": {"sha256": gate.sha256(processor)}}))
        self.doc.update(processed_history_sha256=gate.sha256(history),
            selection_contract_sha256=gate.sha256(contract), processor_sha256=gate.sha256(processor))
        self.save()
        return integration, package.parents[2]

    def test_missing_review_stops_before_csv_parser(self):
        integration, canonical = self.canonical_protocol_package()
        self.freeze.unlink()
        self.campaign.mkdir()
        destination = self.campaign/"I17/application.in"
        destination.parent.mkdir()
        with patch.dict(os.environ, {"MOON_VALIDATION_DATA_ROOT": str(canonical)}), \
             patch.object(integration.csv, "DictReader", side_effect=AssertionError("science CSV parser was reached")):
            with self.assertRaises(integration.MissingPrerequisite):
                integration.prepare_i17_input(self.template, destination,
                    self.campaign, executable=self.exe)
        self.assertFalse(destination.exists())
        self.assertFalse((self.campaign/"_m2_validation_data").exists())

    def test_approved_release_stages_selected_tokens_and_detects_tamper(self):
        integration, canonical = self.canonical_protocol_package()
        self.campaign.mkdir()
        destination = self.campaign/"I17/application.in"
        destination.parent.mkdir()
        with patch.dict(os.environ, {"MOON_VALIDATION_DATA_ROOT": str(canonical)}):
            integration.prepare_i17_input(self.template, destination,
                self.campaign, executable=self.exe)
            release = gate.staged_release(destination)
            selected = Path(integration.moon_setting(destination,
                "diviner_eclipse_history_file", preserve_case=True))
            self.assertIn("december-filter-token", selected.read_text())
            self.assertNotIn("june-filter-token", selected.read_text())
            self.assertEqual(release["document"], self.doc)
            reader = csv.DictReader
            def refuse_second_holdout_read(stream, *args, **kwargs):
                # Revalidating the already visible June evidence is allowed;
                # the claim must stop a second December history inspection.
                if Path(stream.name).name in ("eclipse_thermal_history.csv",
                                               "holdout_eclipse_thermal_history.csv"):
                    raise AssertionError("second December science read")
                return reader(stream, *args, **kwargs)
            with patch.object(integration.csv, "DictReader", side_effect=refuse_second_holdout_read):
                with self.assertRaises(integration.MissingPrerequisite):
                    integration.prepare_i17_input(self.template, destination,
                        self.campaign, executable=self.exe)
        frozen = gate.freeze_path(destination)
        frozen.write_text(frozen.read_text()+" ")
        with self.assertRaises(ValueError): gate.staged_release(destination)

    def test_holdout_numerical_gate_preserves_selection(self):
        from i17_convergence import compare_time_steps
        directories = [self.root/str(index) for index in range(3)]
        for directory, model in zip(directories, (100, 101, 101.05)):
            directory.mkdir()
            row = self.model_rows()[0]
            row["utc"] = "protocol-only"
            row["model_bolometric_brightness_temperature_K"] = str(model)
            with (directory/"receipt.txt.I17_holdout_model.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, list(row)); writer.writeheader(); writer.writerow(row)
        passed, metrics = compare_time_steps(directories, self.root/"dt.csv", "holdout")
        self.assertTrue(passed)
        self.assertEqual(metrics["selection"], "holdout")
        path = directories[-1]/"receipt.txt.I17_holdout_model.csv"
        path.write_text(path.read_text().replace("holdout", "training"))
        with self.assertRaises(ValueError): compare_time_steps(directories, self.root/"dt.csv", "holdout")

    def test_direct_scorer_cannot_bypass_release_before_csv_or_plot(self):
        import integration_common as integration
        work = self.root/"unreleased-score"; work.mkdir()
        (work/"application.in").write_bytes(self.template.read_bytes())
        with patch.object(integration.csv, "DictReader", side_effect=AssertionError("unreleased science read")), \
             patch.object(integration, "render_comparison") as render:
            with self.assertRaises(integration.MissingPrerequisite):
                integration.check("I17", {"comparison_selection": "holdout"}, work)
            render.assert_not_called()

    def test_unidentifiable_phase_cannot_publish_a_comparison_plot(self):
        import integration_common as integration
        work = self.root/"ambiguous-protocol-score"; work.mkdir()
        rows = self.model_rows()
        rows[2]["observed_bolometric_brightness_temperature_K"] = "80"
        with (work/"receipt.txt.I17_holdout_model.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, list(rows[0])); writer.writeheader(); writer.writerows(rows)
        receipt = {"comparison_selection": "holdout", "holdout_comparison_row_count": 3,
            "holdout_rows_accessed": 3, "initialization_mode": "independent-observation-free-periodic-spinup",
            "capability_status": "EXPERIMENTAL / NOT VALIDATED"}
        with patch.object(gate, "staged_release", return_value=self.release()), \
             patch.object(integration, "render_comparison") as render:
            self.assertEqual(integration.check("I17", receipt, work)[0], "SKIPPED")
            render.assert_not_called()
        self.assertFalse((work/"comparison.csv").exists())

    def test_review_outside_committed_tree_cannot_be_promoted(self):
        self.git_patch.stop()
        with self.assertRaises(gate.HoldoutUnavailable):
            gate.reviewed_git_identity(self.freeze)


class NativeGuardContract(unittest.TestCase):
    """Compile a thin adapter around the actual production native guard."""
    @classmethod
    def setUpClass(cls):
        compiler = shlex.split(os.environ.get("CXX", "c++"))
        if not compiler or shutil.which(compiler[0]) is None:
            raise unittest.SkipTest("C++ compiler unavailable for native authorization protocol")
        cls.temporary = tempfile.TemporaryDirectory(prefix="i17-native-guard-protocol-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.build = Path(cls.temporary.name)
        source = cls.build/"native_guard.cpp"
        source.write_text('#include "MoonI17Authorization.h"\n#include <iostream>\n'
            'int main(int argc,char** argv){ if(argc!=4) return 2; std::string error;'
            'bool ok=Moon::I17Authorization::Validate(argv[1],argv[2],argv[3],&error);'
            'if(!ok) std::cerr<<error<<"\\n"; return ok?0:77;}\n')
        cls.executable = cls.build/"native_guard"
        header = TEST_ROOT.parent/"MoonI17Authorization.h"
        command = compiler+["-std=c++17", "-O2", "-I", str(TEST_ROOT.parent),
                            str(source), "-o", str(cls.executable)]
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=False)
        if result.returncode:
            raise RuntimeError("production native authorization guard did not compile:\n"+result.stdout)
        print("I17 native authorization protocol: production header SHA-256="+
              gate.sha256(header)+" compiler="+compiler[0], flush=True)

    def setUp(self):
        self.fixture = tempfile.TemporaryDirectory(prefix="i17-native-guard-data-")
        self.addCleanup(self.fixture.cleanup)
        self.root = Path(self.fixture.name)
        self.history = self.root/"holdout_history.csv"
        self.freeze = self.root/"reviewed.json"
        self.validation = self.root/"holdout_release_validation.json"
        for path in (self.history, self.freeze, self.validation):
            path.write_text("PROTOCOL ONLY: NO SCIENCE VALUES\n")
        self.auth = self.root/"authorization.txt"
        gate.write_native_authorization(self.auth, history=self.history,
            freeze=self.freeze, validation=self.validation,
            release={"sha256": gate.sha256(self.freeze), "source_identity_sha256": "0"*64,
                     "physics_identity_sha256": "1"*64,
                     "reviewed_git_identity": {"blob_oid": "a"*40}})

    def run_guard(self):
        return subprocess.run([str(self.executable), str(self.auth), str(self.history),
            str(self.freeze)], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)

    def test_actual_native_guard_accepts_only_matching_staged_context(self):
        self.assertEqual(self.run_guard().returncode, 0)

    def test_actual_native_guard_rejects_missing_release(self):
        self.auth.unlink()
        result = self.run_guard()
        self.assertEqual(result.returncode, 77)
        self.assertIn("sidecar unavailable", result.stderr)

    def test_actual_native_guard_rejects_unknown_or_duplicate_keys(self):
        original = self.auth.read_text()
        for suffix in ("unknown=value\n", "selection=holdout\n"):
            self.auth.write_text(original+suffix)
            with self.subTest(suffix=suffix): self.assertEqual(self.run_guard().returncode, 77)

    def test_actual_native_guard_rejects_wrong_selection_or_review(self):
        original = self.auth.read_text()
        for old, new in (("selection=holdout", "selection=training"),
                         ("review_status=APPROVED_FROZEN", "review_status=REVIEW_REQUIRED")):
            self.auth.write_text(original.replace(old, new))
            with self.subTest(new=new): self.assertEqual(self.run_guard().returncode, 77)

    def test_actual_native_guard_rejects_other_history_or_validation(self):
        original = self.auth.read_text()
        for old, new in ((str(self.history), str(self.root/"another.csv")),
                         (str(self.validation), str(self.root/"another.json"))):
            self.auth.write_text(original.replace(old, new))
            with self.subTest(new=new): self.assertEqual(self.run_guard().returncode, 77)

    def test_actual_native_guard_rejects_bad_hash_or_review_token(self):
        original = self.auth.read_text()
        for old, new in (("history_sha256="+gate.sha256(self.history), "history_sha256=invalid"),
                         ("review_token="+gate.sha256(self.freeze), "review_token="+"f"*64)):
            self.auth.write_text(original.replace(old, new))
            with self.subTest(new=new): self.assertEqual(self.run_guard().returncode, 77)

    def test_actual_native_guard_requires_real_staged_files(self):
        self.history.unlink()
        self.assertEqual(self.run_guard().returncode, 77)


if __name__ == "__main__":
    unittest.main()

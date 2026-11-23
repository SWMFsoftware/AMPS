#!/usr/bin/env python3
"""Software-only negative controls for routine membership/report/protocol routing.

No fake result in this file represents a particle/model calculation. The fake
execution boundary is explicit and labelled, letting malformed child reports
and crashes be tested without a missing-physics fallback or a native gate edit.
"""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("prototype", HERE / "routine_fixtures.py")
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)


class ProtocolTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="sep-routing-software-fixture-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.executable = self.root / "labelled-software-executable"
        self.executable.write_text("SOFTWARE PROTOCOL FIXTURE ONLY; NOT AN AMPS MODEL\n")
        self.manifest = json.loads((HERE / "routine_fixture_manifest.json").read_text())
        self.rows = [{"id": i, "group": "software-fixture", "class": "routine",
                      "initialization": "none", "build_modes": "software-only",
                      "description": "labelled software protocol, no model output"}
                     for i in self.manifest["expected_ids"]]
        # Use actual membership, but self-contained artificial descriptor and
        # route data. The complete catalog is asserted independently below.
        self.manifest["fixtures"] = {row["id"]: {"kind": "native"} for row in self.rows}
        self.refresh()
        self.framework = P.framework_module()
        self.framework._validation_case_ids = lambda: []
        self.calls, self.behavior = [], {}
        self.catalog = self.catalog_text()
        self.framework._run_streaming = self.execute

    def refresh(self):
        self.manifest["routine_inventory_sha256"] = P.inventory_digest(self.rows)
        self.manifest["expected_ids"] = [row["id"] for row in self.rows]

    def catalog_text(self):
        return "ID | group | class | initialization | build_modes | description\n" + "\n".join(
            " | ".join(row[k] for k in ("id", "group", "class", "initialization", "build_modes", "description"))
            for row in self.rows)

    def execute(self, command, cwd, log, environment, timeout):
        self.calls.append(command)
        if "--list-tests" in command:
            log.write_text(self.catalog)
            return 0
        identity = command[command.index("--test") + 1]
        path = Path(command[command.index("--test-json") + 1])
        action = self.behavior.get(identity, "PASS")
        if action == "missing":
            return 2
        row = {"id": identity, "status": "PASS" if action in ("duplicate", "foreign", "crash") else action,
               "message": "SOFTWARE PROTOCOL FIXTURE ONLY", "seed": 917,
               "configuration": ["software_only=true"], "metrics": [{"name": "protocol_marker", "value": 7}],
               "artifacts": [], "duration_seconds": 0.0}
        rows = [row, dict(row)] if action == "duplicate" else [row]
        if action == "foreign":
            rows[0]["id"] = "FOREIGN999"
        path.write_text(json.dumps({"results": rows}))
        return 9 if action == "crash" else 1 if action in ("FAIL", "ERROR") else 0

    def run_driver(self, forwarded=()):
        self.manifest_path = self.root / "fixture-manifest.json"
        self.manifest_path.write_text(json.dumps(self.manifest))
        self.output = self.root / ("run-" + str(len(list(self.root.glob("run-*")))))
        args = argparse.Namespace(amps=self.executable, manifest=self.manifest_path,
                                  output=self.output, timeout=9, native_args=list(forwarded))
        code = P.run(args, self.framework)
        return code, json.loads((self.output / "srcsep-tests.json").read_text())

    def test_exact_223_membership_and_no_extended_leak(self):
        self.assertEqual(len(self.rows), 223)
        self.assertIn("S13AD01", self.manifest["expected_ids"])
        extra = dict(self.rows[0], id="EXTENDED999", **{"class": "extended"})
        self.catalog += "\n" + " | ".join(extra[k] for k in ("id", "group", "class", "initialization", "build_modes", "description"))
        code, report = self.run_driver()
        self.assertEqual(code, 0)
        self.assertEqual(report["totals"]["passed"], 223)
        self.assertEqual([r["id"] for r in report["results"]], self.manifest["expected_ids"])
        self.assertTrue(all("EXTENDED999" not in c for c in self.calls))

    def test_missing_manifest_id_retains_all_required_errors(self):
        self.manifest["expected_ids"].pop()
        code, report = self.run_driver()
        self.assertEqual(code, 2)
        self.assertEqual(len(report["results"]), len(self.rows))
        self.assertTrue(all(r["status"] == "ERROR" for r in report["results"]))
        self.assertEqual(len(self.calls), 1)

    def test_metadata_mutation_does_not_run_any_callback(self):
        self.catalog = self.catalog.replace("software-only", "changed-mode", 1)
        code, report = self.run_driver()
        self.assertEqual((code, len(report["results"]), len(self.calls)), (2, len(self.rows), 1))

    def test_unavailable_is_required_error_with_continuation(self):
        identity = self.rows[20]["id"]
        self.manifest["fixtures"][identity] = {"kind": "unavailable", "reason": "declared fixture missing"}
        code, report = self.run_driver()
        self.assertEqual((code, report["totals"]["passed"], report["totals"]["errors"]), (2, len(self.rows) - 1, 1))
        self.assertEqual(next(r for r in report["results"] if r["id"] == identity)["status"], "ERROR")

    def test_child_fail_skip_crash_missing_duplicate_foreign(self):
        actions = ("FAIL", "SKIP", "crash", "missing", "duplicate", "foreign")
        for row, action in zip(self.rows, actions):
            self.behavior[row["id"]] = action
        code, report = self.run_driver()
        self.assertEqual(code, 2)
        self.assertEqual([r["status"] for r in report["results"][:6]], ["FAIL", "SKIP", "ERROR", "ERROR", "ERROR", "ERROR"])
        self.assertEqual(len(report["results"]), len(self.rows))
        raw = self.output / self.rows[2]["id"] / "native.json"
        self.assertEqual(json.loads(raw.read_text())["results"][0]["status"], "PASS")
        self.assertEqual(json.loads((raw.parent / "process.json").read_text())["actual_exit_code"], 9)

    def test_seed_repeat_forwarded_unchanged_and_metrics_preserved(self):
        forwarded = ["--test-seed", "1234567", "--test-repeat", "2"]
        code, report = self.run_driver(forwarded)
        self.assertEqual(code, 0)
        self.assertTrue(all(c[-4:] == forwarded for c in self.calls[1:]))
        self.assertTrue(all(r["seed"] == 917 and r["metrics"][0]["value"] == 7 for r in report["results"]))

    def test_reserved_override_rejected_and_existing_output_preserved(self):
        with self.assertRaisesRegex(ValueError, "reserved"):
            self.run_driver(["--test-output=foreign-directory"])
        marker = self.output / "untouched"
        marker.write_text("preserved")
        args = argparse.Namespace(amps=self.executable, manifest=self.manifest_path,
                                  output=self.output, timeout=9, native_args=[])
        with self.assertRaisesRegex(ValueError, "fresh"):
            P.run(args, self.framework)
        self.assertEqual(marker.read_text(), "preserved")

    def test_checkpoint_protocol_scope_hashes_and_quote(self):
        route = {"kind": "material-checkpoint-v2"}
        for key in ("startup_input", "event"):
            path = self.root / ('quoted "' + key + '.in')
            path.write_text("IMMUTABLE SOFTWARE INPUT\n")
            route[key] = {"path": str(path), "sha256": P.digest(path)}
        request, inputs = P.checkpoint_request(route, self.root, self.root)
        self.assertIn('SEP_CHECKPOINT_FRESH_PROCESS 2\nrole "routine"', request.read_text())
        self.assertIn('\\"startup_input.in', request.read_text())
        inputs["event"].write_text("CHANGED")
        with self.assertRaisesRegex(ValueError, "changed/missing"):
            P.checkpoint_request(route, self.root, self.root)
        foreign = copy.deepcopy(self.manifest)
        foreign["fixtures"][self.rows[0]["id"]] = route
        with self.assertRaisesRegex(ValueError, "another test ID"):
            P.validate_manifest(foreign, self.rows)

    def test_storage_metadata_and_explicit_startup_cannot_be_guessed(self):
        row = next(r for r in self.rows if r["id"] == "FTE01")
        row["initialization"] = "native-storage"
        self.refresh()
        self.manifest["fixtures"]["FTE01"] = {"kind": "native-storage", "startup_args": ["--shock-model", "analytical"]}
        P.validate_manifest(self.manifest, self.rows)
        self.manifest["fixtures"]["FTE01"]["startup_args"] = []
        with self.assertRaisesRegex(ValueError, "startup differs"):
            P.validate_manifest(self.manifest, self.rows)
        self.manifest["fixtures"]["FTE01"] = {"kind": "native"}
        with self.assertRaisesRegex(ValueError, "cannot be assumed"):
            P.validate_manifest(self.manifest, self.rows)

    def test_manifest_argument_override_and_midrun_mutation_fail_closed(self):
        for kind in ("native", "material-checkpoint-v2", "validation", "birth-operator", "unavailable"):
            altered=copy.deepcopy(self.manifest)
            altered["fixtures"][self.rows[0]["id"]] = {"kind":kind, "startup_args":["--test-input", "foreign"]}
            with self.assertRaisesRegex(ValueError, "startup_args are forbidden"):
                P.validate_manifest(altered,self.rows)
        execution=self.framework._run_streaming
        def mutate(command,cwd,log,environment,timeout):
            code=execution(command,cwd,log,environment,timeout)
            if "--test" in command:
                self.manifest_path.write_text(self.manifest_path.read_text()+" ")
            return code
        self.framework._run_streaming=mutate
        code,report=self.run_driver()
        self.assertEqual((code,len(report["results"]),report["totals"]["errors"]),(2,len(self.rows),len(self.rows)))
        self.assertEqual(len(self.calls),2) # discovery + firstnative; remainingblocked
        provenance=json.loads((self.output/"run_manifest.json").read_text())
        self.assertNotEqual(provenance["fixture_manifest_sha256"],provenance["fixture_manifest_final_sha256"])
        self.assertEqual(P.digest(self.output/"parsed-fixture-manifest.json"),provenance["fixture_manifest_sha256"])

    def test_catalog_retains_seed_repeat_isolation_and_tags(self):
        text=self.catalog+"\n  seed: explicit seed policy\n  state/isolation: private software fixture\n  seed-capability: context; repeat-supported: yes\n  tags: repeat:process-isolated mover:parker\n"
        parsed=P.parse_catalog(text)[-1]
        self.assertEqual(parsed["seed_policy"], "explicit seed policy")
        self.assertTrue(parsed["repeat_supported"])
        self.assertEqual(parsed["tags"], ["repeat:process-isolated", "mover:parker"])
        self.assertNotEqual(P.inventory_digest(P.parse_catalog(text)), P.inventory_digest(self.rows))

    def test_birth_seed_repeat_groups_preserve_existing_controls(self):
        self.assertEqual(P.birth_campaign_arguments(["--test-seed", "42", "--test-repeat=2"]),
                         ["--seed", "42", "--repeat-count", "2"])
        self.assertEqual(P.birth_campaign_arguments(["--test-seed=41", "--test-seed", "42"]), ["--seed", "42"])
        with self.assertRaisesRegex(ValueError, "does not support"):
            P.birth_campaign_arguments(["--particle-mover", "parker"])
        with self.assertRaisesRegex(ValueError, "missing"):
            P.birth_campaign_arguments(["--test-repeat"])

    def test_duplicate_and_unknown_catalog_metadata_reject(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            P.parse_catalog(self.catalog + "\n" + self.catalog.splitlines()[1])
        with self.assertRaisesRegex(ValueError, "initialization"):
            P.parse_catalog(self.catalog.replace("| none |", "| guessed-owner |", 1))


if __name__ == "__main__":
    unittest.main(verbosity=2)

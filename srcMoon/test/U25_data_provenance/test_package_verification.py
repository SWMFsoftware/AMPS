#!/usr/bin/env python3
"""Planted-defect controls for the U25 package verifier.

Every case builds a small fixture package in a temporary directory, plants one
defect, and checks that U25 assigns the status required by its contract
(acceptance.json status_mapping).  The fixture contains an arbitrary list of
integers and their sum: it is test infrastructure for U25 itself and is never
validation data, an observation, or a model reference.  Nothing here reads or
writes the canonical validation-data root.

Run from the AMPS root:  python3 -B srcMoon/test/U25_data_provenance/test_package_verification.py
"""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import package_verification as pv  # noqa: E402

# The fixture's processing script.  It writes provenance/QA files with a fresh
# generation time on every run, exactly like the real packages.  @VARIANT@
# selects one planted defect; defects that must appear only during the
# reproduction run check U25_FIXTURE_PHASE, which the tests set around it.
FIXTURE_SCRIPT = r'''
import argparse, hashlib, json, os, sys, time
from datetime import datetime, timezone
from pathlib import Path
VARIANT = "@VARIANT@"
REPRODUCING = os.environ.get("U25_FIXTURE_PHASE") == "reproduce"
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def now(): return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
parser = argparse.ArgumentParser(); parser.add_argument("--root", type=Path, required=True)
root = parser.parse_args().root.resolve()
if VARIANT == "exit_failure" and REPRODUCING: sys.exit(3)
if VARIANT == "escape_write" and REPRODUCING: Path(os.environ["U25_FIXTURE_ESCAPE"]).write_text("escaped\n")
values = [int(item) for item in (root / "raw" / "input.txt").read_text().split()]
(root / "processed").mkdir(exist_ok=True)
text = "sum=%d\n" % sum(values)
if VARIANT == "nondeterministic_product": text += "wall=%r\n" % time.time()
(root / "processed" / "sum.txt").write_text(text)
raw = root / "raw" / "input.txt"
inventory = {"files": [{"path": "raw/input.txt", "sha256": sha(raw), "byte_count": raw.stat().st_size, "immutable": True}]}
(root / "raw_inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
product = root / "processed" / "sum.txt"
provenance = {"generated_utc": now(),
              "outputs": [{"path": "processed/sum.txt", "sha256": sha(product), "byte_count": product.stat().st_size}],
              "processing_scripts": [{"path": "scripts/process.py", "sha256": sha(root / "scripts" / "process.py")}],
              "raw_inventory": {"path": "raw_inventory.json", "sha256": sha(root / "raw_inventory.json")},
              "prohibited_methods": {"synthetic_replacement_data": False}}
if VARIANT == "unregistered_time_field": provenance["run_utc"] = now()
if VARIANT == "manual_edit_declared": provenance["manual_edits"] = True
(root / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
qa = {"status": "PASS", "generated_utc": now()}
(root / "qa_report.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
(root / "qa_report.txt").write_text("Fixture QA\nGenerated UTC: %s\nStatus: PASS\n" % qa["generated_utc"])
lines = []
for path in sorted(root.rglob("*")):
    relative = path.relative_to(root)
    if not path.is_file() or path.name == "SHA256SUMS.txt" or "__pycache__" in relative.parts: continue
    lines.append("%s  %s" % (sha(path), relative.as_posix()))
(root / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
'''

TIMESTAMPS = {"provenance.json": {"json_keys": [["generated_utc"]]},
              "qa_report.json": {"json_keys": [["generated_utc"]]},
              "qa_report.txt": {"line_prefixes": ["Generated UTC: "]}}


def fixture_spec() -> dict:
    return {
        "required_artifacts": ["README.md", "raw_inventory.json", "provenance.json", "qa_report.json", "SHA256SUMS.txt"],
        "ignored_names": ["__pycache__"],
        "hash_manifests": [{"file": "SHA256SUMS.txt", "base": "."}],
        "complete_manifest_directories": ["raw", "processed"],
        "raw_directory": "raw",
        "raw_inventory": {"file": "raw_inventory.json", "list_key": ["files"], "path_key": "path",
                          "sha256_key": "sha256", "bytes_key": "byte_count", "immutable_key": "immutable"},
        "provenance": {
            "file": "provenance.json",
            "outputs": {"list_key": ["outputs"], "path_key": "path", "sha256_key": "sha256",
                        "bytes_key": "byte_count", "covers_directory": "processed"},
            "scripts": {"list_key": ["processing_scripts"], "path_key": "path", "sha256_key": "sha256",
                        "covers_directory": "scripts", "covers_suffix": ".py"},
            "self_hashes": [{"path_key": ["raw_inventory", "path"], "sha256_key": ["raw_inventory", "sha256"]}],
            "must_be_false": [["prohibited_methods", "synthetic_replacement_data"]]},
        "reproduction": {"interpreter": sys.executable,
                         "command": ["{python}", "-B", "{copy}/scripts/process.py", "--root", "{copy}"],
                         "network": False, "timeout_seconds": 60, "copy_exclude": ["__pycache__"],
                         "timestamp_fields": TIMESTAMPS}}


def build_fixture(root: Path, variant: str = "clean", extra_files: dict | None = None) -> Path:
    package = root / "package"
    (package / "raw").mkdir(parents=True)
    (package / "scripts").mkdir()
    (package / "raw" / "input.txt").write_text("1 2 3\n")
    (package / "README.md").write_text("U25 control fixture: test infrastructure, not validation data.\n")
    for name, text in (extra_files or {}).items():
        (package / name).write_text(text)
    (package / "scripts" / "process.py").write_text(FIXTURE_SCRIPT.replace("@VARIANT@", variant))
    environment = {key: value for key, value in os.environ.items() if key != "U25_FIXTURE_PHASE"}
    subprocess.run([sys.executable, "-B", str(package / "scripts" / "process.py"), "--root", str(package)],
                   check=True, env=environment)
    return package


def rewrite_manifest(package: Path) -> None:
    """Re-hash consistently, as a hand edit that also updated the records would."""
    provenance = json.loads((package / "provenance.json").read_text())
    for item in provenance["outputs"]:
        path = package / item["path"]
        item["sha256"] = pv.sha256_file(path)
        item["byte_count"] = path.stat().st_size
    (package / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    lines = []
    for path in sorted(package.rglob("*")):
        relative = path.relative_to(package)
        if path.is_file() and path.name != "SHA256SUMS.txt" and "__pycache__" not in relative.parts:
            lines.append(f"{pv.sha256_file(path)}  {relative.as_posix()}")
    (package / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")


class PackageControls(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory(prefix="u25-controls-")
        self.root = Path(self.directory.name)
        os.environ.pop("U25_FIXTURE_PHASE", None)

    def tearDown(self) -> None:
        os.environ.pop("U25_FIXTURE_PHASE", None)
        os.environ.pop("U25_FIXTURE_ESCAPE", None)
        self.directory.cleanup()

    def evaluate(self, package: Path, spec: dict | None = None) -> dict:
        os.environ["U25_FIXTURE_PHASE"] = "reproduce"
        try:
            return pv.evaluate_package(spec or fixture_spec(), package, self.root / "work")
        finally:
            os.environ.pop("U25_FIXTURE_PHASE", None)

    def test_clean_fixture_passes_with_only_timestamps_normalized(self) -> None:
        package = build_fixture(self.root)
        before = {key: pv.sha256_file(path) for key, path in pv.package_files(package, []).items()}
        result = self.evaluate(package)
        self.assertEqual(result["status"], "PASS", result.get("reason"))
        self.assertEqual(result["reproduction"]["normalized_files"],
                         ["provenance.json", "qa_report.json", "qa_report.txt"])
        after = {key: pv.sha256_file(path) for key, path in pv.package_files(package, []).items()}
        self.assertEqual(before, after, "the canonical fixture must not change")
        self.assertIn("manual_edits", result["integrity"]["metadata_completeness_c6"]["absent"])

    def test_changed_raw_byte_is_error(self) -> None:
        package = build_fixture(self.root)
        (package / "raw" / "input.txt").write_text("1 2 4\n")
        self.assertEqual(self.evaluate(package)["status"], "ERROR")

    def test_unlisted_raw_file_is_error(self) -> None:
        package = build_fixture(self.root)
        (package / "raw" / "unlisted.txt").write_text("x\n")
        result = self.evaluate(package)
        self.assertEqual(result["status"], "ERROR")
        self.assertIn("not in any hash manifest", result["reason"])

    def test_consistent_hand_edit_fails_reproduction(self) -> None:
        package = build_fixture(self.root)
        (package / "processed" / "sum.txt").write_text("sum=7\n")
        rewrite_manifest(package)
        result = self.evaluate(package)
        self.assertEqual(result["status"], "FAIL", result.get("reason"))
        self.assertIn("processed/sum.txt", [item["path"] for item in result["reproduction"]["differences"]])

    def test_nondeterministic_product_fails(self) -> None:
        self.assertEqual(self.evaluate(build_fixture(self.root, "nondeterministic_product"))["status"], "FAIL")

    def test_unregistered_time_field_fails(self) -> None:
        result = self.evaluate(build_fixture(self.root, "unregistered_time_field"))
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("provenance.json", [item["path"] for item in result["reproduction"]["differences"]])

    def test_declared_manual_edit_fails(self) -> None:
        result = self.evaluate(build_fixture(self.root, "manual_edit_declared"))
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("manual_edits", result["reason"])

    def test_reproduction_exit_is_error(self) -> None:
        self.assertEqual(self.evaluate(build_fixture(self.root, "exit_failure"))["status"], "ERROR")

    def test_write_into_canonical_package_is_error(self) -> None:
        package = build_fixture(self.root, "escape_write")
        os.environ["U25_FIXTURE_ESCAPE"] = str(package / "processed" / "escaped.txt")
        result = self.evaluate(package)
        self.assertEqual(result["status"], "ERROR")
        self.assertIn("canonical package changed", result["reason"])

    def test_missing_required_artifact_is_skipped(self) -> None:
        package = build_fixture(self.root)
        (package / "README.md").unlink()
        self.assertEqual(self.evaluate(package)["status"], "SKIPPED")

    def test_symlink_is_error(self) -> None:
        package = build_fixture(self.root)
        (package / "processed" / "link.txt").symlink_to(package / "processed" / "sum.txt")
        self.assertEqual(self.evaluate(package)["status"], "ERROR")

    def test_unregistered_reproduction_is_skipped(self) -> None:
        spec = fixture_spec()
        spec["reproduction"] = None
        result = self.evaluate(build_fixture(self.root), spec)
        self.assertEqual(result["status"], "SKIPPED")
        self.assertTrue(result["integrity_evaluated"])


class DatasetGates(unittest.TestCase):
    """Readiness and observation-selection gates through evaluate_dataset()."""

    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory(prefix="u25-gates-")
        self.root = Path(self.directory.name)
        self.contract = {"accepted_readiness_values": ["READY", "PASS"]}

    def tearDown(self) -> None:
        self.directory.cleanup()

    def registry(self, **entry: object) -> dict:
        base = {"package_directory": "package", "readiness": {"file": "qa_report.json", "key": ["status"]},
                "observation_dataset": False, "verification": "fixture"}
        base.update(entry)
        return {"packages": {"fixture": fixture_spec()}, "datasets": {"D90": base}}

    def run_gate(self, registry: dict, dataset: str = "D90") -> dict:
        return pv.evaluate_dataset(dataset, registry, self.contract, self.root, self.root / "work", {})

    def test_ready_dataset_passes(self) -> None:
        build_fixture(self.root)
        self.assertEqual(self.run_gate(self.registry())["status"], "PASS")

    def test_partial_qa_is_skipped_without_integrity(self) -> None:
        package = build_fixture(self.root)
        (package / "qa_report.json").write_text(json.dumps({"status": "PARTIAL"}))
        result = self.run_gate(self.registry())
        self.assertEqual(result["status"], "SKIPPED")
        self.assertIn("status=PARTIAL", result["reason"])
        self.assertFalse(result["integrity_evaluated"])

    def test_corrupt_qa_is_error(self) -> None:
        package = build_fixture(self.root)
        (package / "qa_report.json").write_text("{not json")
        self.assertEqual(self.run_gate(self.registry())["status"], "ERROR")

    def test_absent_registered_key_is_error(self) -> None:
        build_fixture(self.root)
        registry = self.registry(readiness={"file": "qa_report.json", "key": ["overall_status"]})
        self.assertEqual(self.run_gate(registry)["status"], "ERROR")

    def test_unfrozen_selection_is_skipped(self) -> None:
        build_fixture(self.root, extra_files={"selection.json": json.dumps({"status": "PROVISIONAL_NOT_FROZEN"})})
        registry = self.registry(observation_dataset=True, selection_files=[
            {"file": "selection.json", "status_key": ["status"], "frozen_status_values": []}])
        result = self.run_gate(registry)
        self.assertEqual(result["status"], "SKIPPED")
        self.assertIn("not frozen", result["reason"])

    def test_observation_dataset_without_selection_is_skipped(self) -> None:
        build_fixture(self.root)
        self.assertEqual(self.run_gate(self.registry(observation_dataset=True, selection_files=[]))["status"], "SKIPPED")

    def test_registered_frozen_selection_passes(self) -> None:
        build_fixture(self.root, extra_files={"selection.json": json.dumps({"status": "FIXTURE_FROZEN"})})
        registry = self.registry(observation_dataset=True, selection_files=[
            {"file": "selection.json", "status_key": ["status"], "frozen_status_values": ["FIXTURE_FROZEN"]}])
        self.assertEqual(self.run_gate(registry)["status"], "PASS")

    def test_unregistered_verification_and_dataset_are_skipped(self) -> None:
        build_fixture(self.root)
        self.assertEqual(self.run_gate(self.registry(verification=None))["status"], "SKIPPED")
        self.assertEqual(self.run_gate(self.registry(), dataset="D91")["status"], "SKIPPED")


class Rules(unittest.TestCase):
    def test_timestamp_line_rule(self) -> None:
        a = '{\n  "generated_utc": "2026-10-07T06:14:47.967418Z",\n  "x": 1\n}\n'
        b = a.replace("2026-10-07T06:14:47.967418Z", "2026-11-23T23:53:20Z")
        self.assertIsNone(pv.timestamp_only_lines(a, a, ["generated_utc"], []))
        self.assertIsNone(pv.timestamp_only_lines(a, b, ["generated_utc"], []))
        self.assertIsNotNone(pv.timestamp_only_lines(a, b, ["processing_utc"], []))
        self.assertIsNotNone(pv.timestamp_only_lines(a, a.replace("2026-10-07T06:14:47.967418Z", "later"), ["generated_utc"], []))
        self.assertIsNotNone(pv.timestamp_only_lines(a, b.replace('  "generated', '   "generated'), ["generated_utc"], []))
        self.assertIsNotNone(pv.timestamp_only_lines(a, b + "\n", ["generated_utc"], []))
        self.assertIsNotNone(pv.timestamp_only_lines(a, b.replace('"x": 1', '"x": 2'), ["generated_utc"], []))
        t = "Generated UTC: 2026-10-07T06:14:48.165570Z\n"
        self.assertIsNone(pv.timestamp_only_lines(t, "Generated UTC: 2026-11-23T23:53:20.967251Z\n", [], ["Generated UTC: "]))
        self.assertIsNotNone(pv.timestamp_only_lines(t, "Generated UTC: yesterday\n", [], ["Generated UTC: "]))

    def test_same_name_outside_registered_path_is_not_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            left, right = Path(name) / "a.json", Path(name) / "b.json"
            document = {"generated_utc": "2026-10-07T06:14:47Z", "inner": {"generated_utc": "2026-10-07T06:14:47Z"}}
            left.write_text(json.dumps(document, indent=2) + "\n")
            changed = copy.deepcopy(document)
            changed["inner"]["generated_utc"] = "2026-11-23T00:00:00Z"
            right.write_text(json.dumps(changed, indent=2) + "\n")
            self.assertIsNotNone(pv.compare_normalized(left, right, {"json_keys": [["generated_utc"]]}))
            changed = copy.deepcopy(document)
            changed["generated_utc"] = "2026-11-23T00:00:00Z"
            right.write_text(json.dumps(changed, indent=2) + "\n")
            self.assertIsNone(pv.compare_normalized(left, right, {"json_keys": [["generated_utc"]]}))

    def test_dataset_token_parser(self) -> None:
        self.assertEqual(pv.expand_dataset_tokens(["instrument-specific D07-D09/D12-D15 metadata as applicable"]),
                         ["D07", "D08", "D09", "D12", "D13", "D14", "D15"])
        self.assertEqual(pv.expand_dataset_tokens(["D10 global QA-clean potassium geography", "Killen (2002)"]), ["D10"])
        self.assertEqual(pv.expand_dataset_tokens(["target-specific Dxx package"]), [])

    def test_committed_contracts_and_registry_agree(self) -> None:
        """Every dataset declared by a requesting test has a registry entry,
        and the registry uses only the contract's frozen timestamp names."""
        contract, registry = pv.load_contract(), pv.load_registry()
        pv.check_registry(contract, registry)
        declared = pv.declared_datasets()
        self.assertIn("U15", declared)
        for test_id, datasets in declared.items():
            for dataset in datasets:
                self.assertIn(dataset, registry["datasets"], f"{test_id} declares unregistered {dataset}")


if __name__ == "__main__":
    unittest.main(verbosity=2)

#!/usr/bin/env python3
"""Exercise I04 reference identity/error handling, not new physical evidence."""
from pathlib import Path
import importlib.util
import json
import tempfile
import sys
import unittest
from unittest.mock import patch
from subprocess import CompletedProcess

path = Path(__file__).with_name("check_reference.py")
spec = importlib.util.spec_from_file_location("reference_checker", path)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class ReferenceIdentity(unittest.TestCase):
    def test_canonical_manifest_does_not_claim_unused_spice_dependencies(self):
        """Check metadata with a protocol file, not a model/physics fixture."""
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        import integration_common as integration
        with tempfile.TemporaryDirectory(prefix="amps-i04-manifest-protocol-") as directory:
            executable = Path(directory) / "not-a-model-executable"
            executable.write_text("manifest protocol fixture; never executed\n")
            record = integration.manifest("I04", executable,
                ["canonical-test_Moon"], None, mpi_ranks=4)
            self.assertEqual(record["configuration"], "canonical-test_Moon")
            self.assertFalse(record["uses_spice_kernels"])
            self.assertEqual(record["external_data"], {})
            self.assertEqual(record["mpi_ranks"], 4)

    def prepare_layout(self, root):
        """Create only rank/file-layout fixtures, never validation values."""
        (root / "PT/plots").mkdir(parents=True)
        for rank in range(4):
            (root / f"PT/thread={rank}.log").touch()

    def test_bundled_hash_is_the_existing_reviewed_contract(self):
        self.assertEqual(checker.hashlib.sha256(checker.REFERENCE.read_bytes()).hexdigest(),
                         checker.expected_hash())

    def test_corrupt_reference_is_error_before_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            corrupt = root / "corrupt.ref"
            corrupt.write_bytes(b"corrupt storage fixture; not science values")
            diff = root / "check.diff"
            with patch.object(checker, "REFERENCE", corrupt):
                self.assertEqual(checker.check(root, diff, root / "absent-diffnum"), 2)
            evidence = json.loads(diff.with_suffix(".reference.json").read_text())
            self.assertIn("checksum", evidence["reason"])
            self.assertGreater(diff.stat().st_size, 0)

    def test_missing_rank_layout_is_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            diff = root / "check.diff"
            self.assertEqual(checker.check(root, diff, root / "absent-diffnum"), 2)
            evidence = json.loads(diff.with_suffix(".reference.json").read_text())
            self.assertIn("four MPI ranks", evidence["reason"])

    def test_corrupt_provenance_cannot_bless_matching_reference_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            corrupt = root / "provenance.json"
            corrupt.write_text('{"status":"not reviewed"}\n')
            diff = root / "check.diff"
            with patch.object(checker, "PROVENANCE", corrupt):
                self.assertEqual(checker.check(root, diff, root / "absent-diffnum"), 2)
            evidence = json.loads(diff.with_suffix(".reference.json").read_text())
            self.assertIn("provenance checksum", evidence["reason"])

    def test_missing_production_output_is_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.prepare_layout(root)
            diff = root / "check.diff"
            self.assertEqual(checker.check(root, diff, root / "absent-diffnum"), 2)
            evidence = json.loads(diff.with_suffix(".reference.json").read_text())
            self.assertIn("production regression output is missing", evidence["reason"])

    def test_comparator_failure_and_evaluated_difference_remain_distinct(self):
        # These are comparator protocol messages, not physical science data.
        # The same Perl exit code represents an I/O error and a real mismatch.
        # Check that the runner preserves the distinction required by AGENTS.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.prepare_layout(root)
            (root / "PT/plots/test_Moon.dat").touch()
            diff = root / "check.diff"
            for diagnostic, expected in [
                    (b"ERROR in DiffNum.pl: cannot open result\n", 2),
                    (b"ERROR in DiffNum.pl: number of significant numerical differences is 4\n", 1)]:
                with self.subTest(diagnostic=diagnostic), patch.object(
                        checker.subprocess, "run", return_value=CompletedProcess(
                            args=[], returncode=255, stdout=diagnostic, stderr=b"")):
                    self.assertEqual(checker.check(root, diff, root / "diffnum"), expected)
                    self.assertEqual(diff.read_bytes(), diagnostic)

    def test_m0_driver_recognizes_only_this_runs_evaluated_mismatch(self):
        # Protocol-only fixtures verify the aggregate's FAIL/ERROR boundary.
        # No fixture data are compared or treated as model/observational truth.
        test_root=Path(__file__).resolve().parents[2]
        driver_spec=importlib.util.spec_from_file_location(
            "i04_m0_driver",test_root/"run_integration_tests.py")
        driver=importlib.util.module_from_spec(driver_spec)
        with patch.object(sys,"path",[str(test_root)]+sys.path):
            driver_spec.loader.exec_module(driver)
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);nightly=root/"nightly"
            self.prepare_layout(nightly)
            (nightly/"amps").touch()
            result=nightly/"PT/plots/test_Moon.dat";result.touch()
            receipt=root/"test_Moon.reference.json"
            def write(status,compared):
                receipt.write_text(json.dumps({"status":status,
                    "command":["DiffNum.pl",str(compared),"reviewed.ref"]}))
            write("FAIL",result)
            self.assertTrue(driver.evaluated_reference_failure(nightly,receipt))
            write("ERROR",result)
            self.assertFalse(driver.evaluated_reference_failure(nightly,receipt))
            write("FAIL",root/"previous-run.dat")
            self.assertFalse(driver.evaluated_reference_failure(nightly,receipt))
            write("FAIL",result);result.unlink()
            self.assertFalse(driver.evaluated_reference_failure(nightly,receipt))


if __name__ == "__main__":
    # These infrastructure tests deliberately exercise corrupt references,
    # absent files and evaluated differences. The checker correctly prints
    # ERROR/FAIL for those fixtures, but successful self-tests must not make
    # the surrounding production regression look broken. unittest buffers
    # each test's stdout/stderr and replays it if that test actually fails;
    # its assertion traceback and nonzero exit status remain visible too.
    # The separate production checker command in MakefileTest/Table is not
    # buffered here, so its real PASS/FAIL/ERROR is always printed normally.
    unittest.main(buffer=True)

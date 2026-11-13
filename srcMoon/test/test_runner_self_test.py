#!/usr/bin/env python3
"""Infrastructure self-checks for the global srcMoon campaign runner.

This file deliberately has no U/I identifier: it verifies orchestration and
contract discovery, not lunar physics.  In particular, it cannot turn a
missing production capability or dataset into a PASS for a roadmap test.
"""

from __future__ import annotations

import importlib.util
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


TEST_ROOT = Path(__file__).resolve().parent


def load_runner():
    """Load the adjacent executable module without changing Python search paths."""
    path = TEST_ROOT / "test_runner.py"
    specification = importlib.util.spec_from_file_location(
        "srcmoon_global_test_runner", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load global runner from {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()


class GlobalRunnerContractTest(unittest.TestCase):
    """Check the frozen registry and every implemented per-ID contract."""

    def test_implemented_ids_are_exact(self):
        expected = ({f"U{number:02d}" for number in range(1, 27)} |
                    {f"I{number:02d}" for number in range(1, 21)})
        self.assertEqual(set(RUNNER.implemented_layout()), expected)

    def test_planned_ids_are_exact_and_not_implemented(self):
        expected = {f"I{number:02d}" for number in range(21, 47)}
        self.assertEqual(set(RUNNER.PLANNED_I_TESTS), expected)
        self.assertTrue(expected.isdisjoint(RUNNER.implemented_layout()))

    def test_architecture_contract_is_complete(self):
        self.assertEqual(RUNNER.audit_architecture(), [])

    def test_missing_capability_record_is_skipped(self):
        record = RUNNER.skipped_record("I21", "planned", "not implemented")
        self.assertEqual(record["status"], "SKIPPED")
        self.assertEqual(record["execution"], "not-run")
        self.assertIsNone(record["runner"])

    def test_supporting_results_stay_separate_and_control_exit(self):
        """Exercise orchestration only; fixture statuses are not science evidence.

        A fake child writes its protocol receipt instead of compiling/running
        AMPS. Verify that a supporting FAIL/ERROR cannot disappear behind a
        formal predecessor SKIPPED, and supporting PASS never becomes I19 PASS.
        """
        for status, expected in (("PASS", 0), ("FAIL", 1), ("ERROR", 2)):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as name:
                output=Path(name)/"campaign"
                def fake_child(command, log, stream):
                    child=Path(command[command.index("--output-dir")+1])
                    child.mkdir(parents=True)
                    (child/"development_summary.json").write_text(json.dumps({
                        "results":[{"id":"I19","status":"SKIPPED"}],
                        "supporting_checks":[{"check":"I19-production-support",
                            "status":status,"runner":"fixture-runner","log":"fixture.log"}]}))
                    return expected
                with patch.object(sys,"argv",["test_runner.py","--phase","m2-development",
                        "--output-dir",str(output)]), patch.object(RUNNER,"run_logged",fake_child), \
                        contextlib.redirect_stdout(io.StringIO()) as stream:
                    self.assertEqual(RUNNER.main(),expected)
                report=json.loads((output/"test_summary.json").read_text())
                self.assertEqual(report["counts"]["SKIPPED"],1)
                self.assertEqual(report["counts"]["PASS"],0)
                self.assertEqual(report["development_counts"][status],1)
                if expected:
                    self.assertIn("fixture-runner",stream.getvalue())
                    self.assertIn("fixture.log",stream.getvalue())


if __name__ == "__main__":
    unittest.main()

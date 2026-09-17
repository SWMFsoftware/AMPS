#!/usr/bin/env python3
"""Unit tests for the relocated SWCME R1 runner's public selection contract."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


RUNNER_PATH = Path(__file__).resolve().parents[1] / "run_tests.py"
SPEC = importlib.util.spec_from_file_location("swcme_r1_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
RUNNER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


def arguments(**overrides):
    values = dict(list=False, tests=[], groups=[], routine=False, all=False,
                  output_dir=Path("unused"), jobs=1, rebuild=False,
                  no_build=True, verbose=False)
    values.update(overrides)
    return argparse.Namespace(**values)


class SelectionTests(unittest.TestCase):
    def test_all_is_stable_and_complete(self):
        selected = RUNNER.select(arguments(all=True))
        self.assertEqual([item.test_id for item in selected],
                         [item.test_id for item in RUNNER.TESTS])

    def test_ids_are_case_insensitive_and_deduplicated(self):
        selected = RUNNER.select(
            arguments(tests=["r1cfg01", "R1CFG01", "r3d01"]))
        self.assertEqual([item.test_id for item in selected],
                         ["R1CFG01", "R3D01"])

    def test_groups_expand_in_registry_order(self):
        selected = RUNNER.select(arguments(groups=["common", "3d"]))
        self.assertEqual([item.test_id for item in selected],
                         ["R1CFG01", "R3D01"])

    def test_missing_or_conflicting_mode_is_usage_error(self):
        with self.assertRaises(RUNNER.UsageError):
            RUNNER.select(arguments())
        with self.assertRaises(RUNNER.UsageError):
            RUNNER.select(arguments(routine=True, all=True))
        with self.assertRaises(RUNNER.UsageError):
            RUNNER.select(arguments(tests=["missing"]))


class ReportTests(unittest.TestCase):
    def test_json_and_junit_are_written_from_same_results(self):
        records = [
            RUNNER.Result("R1CFG01", "COMMON", "PASS", "ok", 0.1, ["x"]),
            RUNNER.Result("R3D01", "3D", "FAIL", "bad", 0.2, ["y"]),
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            RUNNER.write_reports(root, records)
            json_text = (root / "swcme-r1-tests.json").read_text(encoding="utf-8")
            xml_text = (root / "swcme-r1-tests.xml").read_text(encoding="utf-8")
            self.assertIn('"PASS": 1', json_text)
            self.assertIn('"FAIL": 1', json_text)
            self.assertIn('tests="2"', xml_text)
            self.assertIn('failures="1"', xml_text)


if __name__ == "__main__":
    unittest.main(verbosity=2)

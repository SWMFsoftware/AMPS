#!/usr/bin/env python3
"""Regression tests for shared C9/C10 executable and status provenance guards."""

from __future__ import annotations

import json
import os
import shlex
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, Tuple


TEST_ROOT = Path(__file__).resolve().parents[1]
if str(TEST_ROOT) not in sys.path:
    sys.path.insert(0, str(TEST_ROOT))

from earth_test_runtime import (  # noqa: E402
    WORKDIR_KEY_IGNORED_OPTIONS,
    EarthExecutableError,
    archive_previous_result,
    atomic_write_json,
    configuration_workdir,
    finish_run_status,
    initial_run_status,
    probe_earth_executable,
    record_workdir_configuration,
)


class EarthTestRuntimeTests(unittest.TestCase):
    def _make_driver(self, root: Path, body: str) -> Path:
        root.mkdir(parents=True, exist_ok=True)
        path = root / "amps"
        path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
        path.chmod(0o755)
        return path

    def test_earth_help_identity_and_hash_are_required(self) -> None:
        with tempfile.TemporaryDirectory(prefix="earth-runtime-") as directory:
            root = Path(directory)
            earth = self._make_driver(
                root,
                "printf '%s\\n' 'AMPS Earth: Gridless Energetic Particle Solver' "
                "'-mode <3d|gridless>' 'No solver is run.'\n",
            )
            identity = probe_earth_executable(earth, cwd=root)
            self.assertEqual(identity["path"], str(earth.resolve()))
            self.assertRegex(identity["sha256"], r"^[0-9a-f]{64}$")

            sep = self._make_driver(root / "sep", "printf '%s\\n' 'SEP standalone driver'\n")
            with self.assertRaisesRegex(EarthExecutableError, "not the required Earth driver"):
                probe_earth_executable(sep, cwd=root)

    def test_status_write_is_atomic_and_result_archive_is_recoverable(self) -> None:
        with tempfile.TemporaryDirectory(prefix="earth-status-") as directory:
            root = Path(directory)
            status = root / "run_status.json"
            payload = initial_run_status(
                schema="test-status-v1",
                test_id="C9",
                output_root=root,
                argv=["--dry-run"],
                execution_mode="DRY_RUN",
            )
            payload["status"] = "RUNNING"
            atomic_write_json(status, payload)
            self.assertIn('"status": "RUNNING"', status.read_text(encoding="utf-8"))
            self.assertIn('"result_current": false', status.read_text(encoding="utf-8"))
            self.assertEqual(list(root.glob(".*.tmp.*")), [])

            finished = finish_run_status(
                payload,
                state="DRY_RUN_COMPLETE",
                return_code=0,
                message="commands generated",
                result_current=False,
            )
            atomic_write_json(status, finished)
            self.assertIn(
                '"status": "DRY_RUN_COMPLETE"', status.read_text(encoding="utf-8")
            )
            self.assertIsNotNone(finished["finished_utc"])

            result = root / "C9_result.json"
            result.write_text('{"passed": true}\n', encoding="utf-8")
            archived = archive_previous_result(
                result, run_started_utc="2026-10-08T22:00:00Z"
            )
            self.assertIsNotNone(archived)
            assert archived is not None
            self.assertFalse(result.exists())
            self.assertEqual(archived.read_text(encoding="utf-8"), '{"passed": true}\n')

    def test_c9_and_c10_wire_the_shared_guard_before_physical_results(self) -> None:
        # Keep runner integration in the permanent gate: testing the helper in
        # isolation would not catch a later edit that bypassed it in one of the
        # two expensive observational validations.
        for test_id in ("C9", "C10"):
            source = (TEST_ROOT / test_id / ("run_%s.py" % test_id)).read_text(
                encoding="utf-8"
            )
            self.assertIn("probe_earth_executable(", source, test_id)
            self.assertIn("archive_previous_result(", source, test_id)
            self.assertIn("%s_run_status.json" % test_id, source, test_id)
            self.assertIn('result_current=True', source, test_id)

    def test_configuration_workdir_separates_option_sets(self) -> None:
        base = "test_output/C4_exit"
        # A bare invocation keeps the documented default directory.
        self.assertEqual(configuration_workdir(base, []), base)
        hc4 = configuration_workdir(base, ["--adaptive-dt", "T", "--movers", "HC4"])
        boris = configuration_workdir(base, ["--adaptive-dt", "T", "--movers", "BORIS"])
        self.assertTrue(hc4.startswith(base + "_"), hc4)
        self.assertNotEqual(hc4, boris)
        # The same options give the same directory, and options that only say how
        # to treat an existing directory do not change it, so --skip-run finds it.
        self.assertEqual(
            configuration_workdir(base, ["--adaptive-dt", "T", "--movers", "HC4"]), hc4
        )
        self.assertEqual(
            configuration_workdir(
                base, ["--adaptive-dt", "T", "--skip-run", "--movers", "HC4", "--keep"]
            ),
            hc4,
        )
        for option in WORKDIR_KEY_IGNORED_OPTIONS:
            self.assertEqual(configuration_workdir(base, [option]), base)

    def test_configuration_workdir_is_recorded(self) -> None:
        with tempfile.TemporaryDirectory(prefix="earth-runtime-") as directory:
            workdir = Path(directory) / "C4_exit_x"
            record_workdir_configuration(workdir, ["--movers", "HC4", "--keep"])
            payload = json.loads((workdir / "workdir_configuration.json").read_text())
            self.assertEqual(payload["argv"], ["--movers", "HC4", "--keep"])
            self.assertEqual(payload["key_options"], ["--movers", "HC4"])

    def test_list_entries_of_one_wrapper_get_separate_workdirs(self) -> None:
        # The wrappers delete their default work directory when they start, and
        # test_runner.py -j N runs list entries in parallel.  Every wrapper that
        # the list runs more than once must take its default directory from
        # configuration_workdir(), and its entries must differ in the options
        # that form the key; otherwise two entries share one directory and
        # delete each other's files (the C4, C5 and C14 failures at -j 6).
        wrappers = ("C3", "C4", "C5", "C14")
        for test_id in wrappers:
            source = (TEST_ROOT / test_id / ("run_%s.py" % test_id)).read_text(
                encoding="utf-8"
            )
            self.assertIn("configuration_workdir(", source, test_id)

        keys: Dict[str, Dict[Tuple[str, ...], int]] = {w: {} for w in wrappers}
        list_path = TEST_ROOT / "list"
        for number, line in enumerate(list_path.read_text(encoding="utf-8").splitlines(), 1):
            if line[:2] not in ("P ", "F "):
                continue
            tokens = shlex.split(line[2:])
            if tokens and tokens[0] == "python3":
                tokens = tokens[1:]
            for test_id in wrappers:
                if tokens and tokens[0].endswith("/%s/run_%s.py" % (test_id, test_id)):
                    if "--workdir" in tokens or any(t.startswith("--workdir=") for t in tokens):
                        continue
                    key = tuple(
                        t for t in tokens[1:] if t not in WORKDIR_KEY_IGNORED_OPTIONS
                    )
                    self.assertNotIn(
                        key,
                        keys[test_id],
                        "%s list lines %s and %d share a default work directory"
                        % (test_id, keys[test_id].get(key), number),
                    )
                    keys[test_id][key] = number


if __name__ == "__main__":
    unittest.main(verbosity=2)

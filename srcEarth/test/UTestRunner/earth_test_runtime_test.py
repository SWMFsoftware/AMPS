#!/usr/bin/env python3
"""Regression tests for shared C9/C10 executable and status provenance guards."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[1]
if str(TEST_ROOT) not in sys.path:
    sys.path.insert(0, str(TEST_ROOT))

from earth_test_runtime import (  # noqa: E402
    EarthExecutableError,
    archive_previous_result,
    atomic_write_json,
    finish_run_status,
    initial_run_status,
    probe_earth_executable,
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


if __name__ == "__main__":
    unittest.main(verbosity=2)

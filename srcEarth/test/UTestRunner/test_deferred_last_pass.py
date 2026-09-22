#!/usr/bin/env python3
"""End-to-end regression tests for test_runner.py last-pass persistence."""

from __future__ import annotations

import json
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TEST_DIRECTORY = Path(__file__).resolve().parents[1]
RUNNER = TEST_DIRECTORY / "test_runner.py"


class DeferredLastPassTests(unittest.TestCase):
    """Exercise the public CLI so these tests also guard against reruns."""

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.runner = self.root / "test_runner.py"
        self.test_list = self.root / "list"
        self.counter = self.root / "executions.txt"
        self.worker = self.root / "worker.py"
        self.pending = self.root / ".list.last-pass-results.json"
        shutil.copy2(RUNNER, self.runner)

        # Each launched command appends one label.  The exact number of lines
        # lets the deferred-commit test prove that --commit-last-pass did not
        # silently schedule either command again.
        self.worker.write_text(
            "from pathlib import Path\n"
            "import sys\n"
            "counter = Path(sys.argv[1])\n"
            "old = counter.read_text(encoding='utf-8') if counter.exists() else ''\n"
            "counter.write_text(old + sys.argv[2] + '\\n', encoding='utf-8')\n"
            "raise SystemExit(0 if sys.argv[3] == 'pass' else 1)\n",
            encoding="utf-8",
        )

        python = shlex.quote(sys.executable)
        worker = shlex.quote(str(self.worker))
        counter = shlex.quote(str(self.counter))
        self.original_list = (
            f"P {python} {worker} {counter} scalar pass\n"
            "last pass: old-scalar\n"
            "for $outcome={pass,fail}\n"
            f"{{P,F}} {python} {worker} {counter} $outcome $outcome\n"
            "last pass: {old-loop-pass,old-loop-fail}\n"
        )
        self.test_list.write_text(self.original_list, encoding="utf-8")

    def run_runner(self, *arguments: object) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(self.runner),
                *(str(argument) for argument in arguments),
            ],
            cwd=self.root,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    def execute_test_list(self, *extra_arguments: object) -> subprocess.CompletedProcess[str]:
        return self.run_runner(
            self.test_list,
            "--no-memory-gate",
            "--no-start-progress",
            "--commit-id",
            "tested-commit",
            "--log-dir",
            self.root / "logs",
            "--report-prefix",
            self.root / "report",
            *extra_arguments,
        )

    def assert_run_succeeded(self, completed: subprocess.CompletedProcess[str]) -> None:
        self.assertEqual(completed.returncode, 0, completed.stdout)

    def test_deferred_commit_updates_passes_without_rerunning(self) -> None:
        run = self.execute_test_list()
        self.assert_run_succeeded(run)
        self.assertEqual(self.test_list.read_text(encoding="utf-8"), self.original_list)
        self.assertEqual(
            self.counter.read_text(encoding="utf-8").splitlines(),
            ["scalar", "pass", "fail"],
        )
        self.assertTrue(self.pending.is_file())

        payload = json.loads(self.pending.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(payload["commit_id"], "tested-commit")
        self.assertEqual(payload["test_count"], 3)
        self.assertEqual(payload["passing_count"], 2)

        # The runner copy and list share a directory, exactly reproducing the
        # requested no-positional `test_runner.py --commit-last-pass` form.
        commit = self.run_runner("--commit-last-pass")
        self.assert_run_succeeded(commit)
        self.assertIn("Updated 2 last-pass entries", commit.stdout)
        self.assertFalse(self.pending.exists())
        self.assertEqual(
            self.counter.read_text(encoding="utf-8").splitlines(),
            ["scalar", "pass", "fail"],
            "--commit-last-pass executed test commands again",
        )
        self.assertEqual(
            self.test_list.read_text(encoding="utf-8"),
            self.original_list
            .replace("last pass: old-scalar", "last pass: tested-commit")
            .replace(
                "last pass: {old-loop-pass,old-loop-fail}",
                "last pass: {tested-commit,old-loop-fail}",
            ),
        )

    def test_commit_rejects_results_if_list_changed(self) -> None:
        run = self.execute_test_list()
        self.assert_run_succeeded(run)
        self.test_list.write_text(
            self.original_list + "! edited after the run\n",
            encoding="utf-8",
        )

        commit = self.run_runner("--commit-last-pass")
        self.assertEqual(commit.returncode, 2, commit.stdout)
        self.assertIn("changed after the saved run", commit.stdout)
        self.assertTrue(self.pending.is_file(), "rejected cache should be retained")
        self.assertIn("last pass: old-scalar", self.test_list.read_text(encoding="utf-8"))
        self.assertEqual(len(self.counter.read_text(encoding="utf-8").splitlines()), 3)

    def test_update_last_pass_remains_immediate_and_leaves_no_cache(self) -> None:
        run = self.execute_test_list("--update-last-pass")
        self.assert_run_succeeded(run)
        updated = self.test_list.read_text(encoding="utf-8")
        self.assertIn("last pass: tested-commit", updated)
        self.assertIn("last pass: {tested-commit,old-loop-fail}", updated)
        self.assertFalse(self.pending.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)

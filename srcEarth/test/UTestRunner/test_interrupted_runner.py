#!/usr/bin/env python3
"""UTR-F08: an interrupted runner still writes its reports.

The AMPS nightly stops a runner at its time limit with SIGINT (SIGKILL 60 s
later).  On 2026-11-18 the Earth runner was stopped this way while its last
entry ran, printed "Interrupted by user" and wrote no report, so none of the
63 finished entries reached the nightly results.  These tests run the
production runner as a subprocess on a three-entry list with ``-j 1``, signal
it while the second entry runs, and check the report: the first entry
finished, the second was stopped (``interrupted``) and its command is gone,
and the third never started (``not_started``).  Only the standard library and
``sh``/``sleep``/``true`` are needed.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


RUNNER_PATH = Path(__file__).resolve().parents[1] / "test_runner.py"

# Test adapter for the UTR-F08 signal race (nightly 2026-11-22): load the
# production runner unchanged and wrap its wait_for_process_completion().  Once the
# second entry's command has written its pid, the wrapper sends the signal to the
# runner's own process from inside that entry's coroutine.  This is the situation of
# the 2026-11-22 failure, made deterministic: a runner whose signal handler raises an
# exception in the executing Python frame (before 2026-11-23) loses the entry's
# "# Interrupted" record and leaves its command running; the event-loop handler stops
# the command and writes the record.
SIGNAL_INSIDE_ENTRY_LAUNCHER = """
import asyncio, importlib.util, os, signal, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("earth_test_runner", sys.argv[1])
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
signum, pid_file = int(sys.argv[2]), Path(sys.argv[3])
original_wait = module.wait_for_process_completion
fired = []

async def wait_and_signal_inside_entry(proc, **kwargs):
    if not fired:
        while proc.returncode is None and not (pid_file.exists() and pid_file.read_text().strip()):
            await asyncio.sleep(0.01)
        if proc.returncode is None:
            fired.append(True)
            os.kill(os.getpid(), signum)
    return await original_wait(proc, **kwargs)

module.wait_for_process_completion = wait_and_signal_inside_entry
sys.exit(module.main(sys.argv[4:]))
"""

def process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    try:   # a zombie waiting for its parent counts as gone
        with open(f"/proc/{pid}/stat", encoding="utf-8") as stream:
            return stream.read().split(") ", 1)[1].split()[0] != "Z"
    except (OSError, IndexError):
        return False


class InterruptedRunTests(unittest.TestCase):
    """UTR-F08: SIGINT and SIGTERM stop the run, and every entry is reported."""

    def interrupt(self, sig: int, signal_inside_entry: bool = False):
        with tempfile.TemporaryDirectory(prefix="earth-runner-interrupt-") as tmp:
            root = Path(tmp)
            pid_file = root / "sleeper.pid"
            (root / "list").write_text(
                "P true\n"
                f"P sh -c 'echo $$ > {pid_file}; exec sleep 60'\n"
                "P true\n",
                encoding="utf-8",
            )
            runner_args = [str(root / "list"), "-j", "1", "--no-memory-gate",
                           "--workdir", str(root), "--log-dir", str(root / "logs"),
                           "--report-prefix", str(root / "report")]
            if signal_inside_entry:
                command = [sys.executable, "-c", SIGNAL_INSIDE_ENTRY_LAUNCHER, str(RUNNER_PATH),
                           str(int(sig)), str(pid_file)] + runner_args
            else:
                command = [sys.executable, str(RUNNER_PATH)] + runner_args
            proc = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            )
            try:
                if signal_inside_entry:
                    # The launcher signals the runner itself; only wait for the end.
                    output, _ = proc.communicate(timeout=60)
                    sleeper = int(pid_file.read_text())
                else:
                    deadline = time.monotonic() + 30
                    while not (pid_file.exists() and pid_file.read_text().strip()):
                        self.assertLess(time.monotonic(), deadline, "the second entry did not start")
                        self.assertIsNone(proc.poll(), "the runner ended before the second entry started")
                        time.sleep(0.05)
                    sleeper = int(pid_file.read_text())
                    time.sleep(0.2)
                    proc.send_signal(sig)
                    output, _ = proc.communicate(timeout=60)
            finally:
                if proc.poll() is None:
                    proc.kill()
                    proc.communicate()
            deadline = time.monotonic() + 10
            while process_exists(sleeper) and time.monotonic() < deadline:
                time.sleep(0.05)
            report = json.loads((root / "report.json").read_text(encoding="utf-8"))
            stopped_log = (Path(report[1]["log_file"]).read_text(encoding="utf-8")
                           if report[1]["log_file"] else "")
            csv_header = (root / "report.csv").read_text(encoding="utf-8").splitlines()[0]
            return proc.returncode, output, report, stopped_log, csv_header, process_exists(sleeper)

    def check(self, sig: int, signal_inside_entry: bool = False):
        code, output, report, stopped_log, csv_header, alive = self.interrupt(sig, signal_inside_entry)
        self.assertEqual(code, 130, output)
        self.assertEqual([e["line_no"] for e in report], [1, 2, 3])
        first, second, third = report
        self.assertEqual((first["actual"], first["interrupted"], first["not_started"]), ("P", False, False))
        self.assertEqual((second["actual"], second["interrupted"], second["not_started"]), ("F", True, False))
        self.assertGreater(second["elapsed_s"], 0.0)
        self.assertIn("# Interrupted: the runner was stopped", stopped_log)
        self.assertEqual((third["actual"], third["interrupted"], third["not_started"]), ("F", False, True))
        self.assertEqual(third["log_file"], "")
        self.assertFalse(alive, "the stopped command outlived the runner")
        self.assertIn("[INTERRUPTED] #002 line 2", output)
        self.assertIn("1 stopped while running, 1 not started", output)
        # UTR-F09 appends the decoded status and failure cause after these columns.
        self.assertTrue(csv_header.endswith(
            "exclusive,interrupted,not_started,exit_status,failure_cause"))
        self.assertIn("the runner was interrupted", second["failure_cause"])
        self.assertIn("not started", third["failure_cause"])
        self.assertIn("stopped while running (the runner was interrupted; no result)", output)
        self.assertIn("not started (the runner was interrupted; no result)", output)

    def test_sigint_reports_every_entry(self):
        self.check(signal.SIGINT)

    def test_sigterm_reports_every_entry(self):
        self.check(signal.SIGTERM)

    def test_sigint_inside_running_entry_stops_and_records_it(self):
        self.check(signal.SIGINT, signal_inside_entry=True)

    def test_sigterm_inside_running_entry_stops_and_records_it(self):
        self.check(signal.SIGTERM, signal_inside_entry=True)

if __name__ == "__main__":
    unittest.main(verbosity=2)

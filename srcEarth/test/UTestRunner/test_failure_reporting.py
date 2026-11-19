#!/usr/bin/env python3
"""UTR-F09: a failed entry reports why it failed, not only its exit status.

A test that dies with ``exit=139`` must say that this is 128 + 11, a SIGSEGV
in a process the wrapper ran, and quote the crash message from the log (for
Open MPI: "mpirun noticed that process rank N ... exited on signal 11").  These
tests exercise the production runner (``test_runner.py``), the shared wrapper
helper (``earth_test_runtime.report_amps_failure``), and the C1-C4 wrappers
with a stand-in ``mpirun``.  None of them needs AMPS, MPI or a compiler.

Reporting only: actual P/F still comes from the exit status and the timeout.
Wrappers still return the unchanged nonzero AMPS status.
"""

from __future__ import annotations

import asyncio
import contextlib
import importlib.util
import io
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = TEST_ROOT / "test_runner.py"
SPEC = importlib.util.spec_from_file_location("earth_test_runner_f09", RUNNER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load production test_runner.py")
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)

if str(TEST_ROOT) not in sys.path:
    sys.path.insert(0, str(TEST_ROOT))
import earth_test_runtime as runtime  # noqa: E402


SEP = "=" * 80
# What Open MPI 4.1.7 with UCX wrote when rank 1 of a 4-rank job dereferenced
# NULL (reproduced on the validation host, hostname and PIDs replaced).
OPEN_MPI_SEGV = [
    "rank 0: starting field initialization",
    "[node:101:0:101] Caught signal 11 (Segmentation fault: address not mapped to object at address (nil))",
    "==== backtrace (tid: 101) ====",
    " 3 0x00000000004011e0 main()  ???:0",
    "=================================",
    "[node:101] *** Process received signal ***",
    "[node:101] Signal: Segmentation fault (11)",
    "[node:101] *** End of error message ***",
    "--------------------------------------------------------------------------",
    "Primary job  terminated normally, but 1 process returned",
    "a non-zero exit code. Per user-direction, the job has been aborted.",
    "--------------------------------------------------------------------------",
    "mpirun noticed that process rank 1 with PID 0 on node node exited on signal 11 (Segmentation fault).",
    "--------------------------------------------------------------------------",
]
MPIRUN_LINE = OPEN_MPI_SEGV[-2]


def runner_log(path: Path, body, *, command: str = "srcEarth/test/C1/run_C1.py --mode gridless") -> Path:
    """A log laid out exactly as run_one_test writes it."""

    header = [
        "# Test 14 from line 43",
        "# Expected: P",
        f"# Command: {command}",
        "# Workdir: /tmp",
        "# np=4, nt=16",
        "# Started: 2026-11-19 10:00:00",
        SEP,
    ]
    footer = ["", SEP, "# Finished: 2026-11-19 10:00:03", "# Elapsed seconds: 3.000",
              "# Timed out: False", "# Exit code: 139"]
    path.write_text("\n".join(header + list(body) + footer) + "\n", encoding="utf-8")
    return path


class ExitStatusTests(unittest.TestCase):
    def test_known_statuses_are_decoded_exactly(self):
        cases = {
            None: "no exit status",
            0: "exit 0",
            1: "exit 1",
            139: "exit 139 = 128 + 11 (SIGSEGV: Segmentation fault)",
            134: "exit 134 = 128 + 6 (SIGABRT: Aborted)",
            137: "exit 137 = 128 + 9 (SIGKILL: Killed)",
            -11: "killed by signal 11 (SIGSEGV: Segmentation fault)",
            -15: "killed by signal 15 (SIGTERM: Terminated)",
            255: "exit 255",
        }
        for code, text in cases.items():
            with self.subTest(code=code):
                self.assertEqual(runner.describe_exit_code(code, compact=True), text)
        self.assertIn("by the shell/MPI-launcher convention", runner.describe_exit_code(139))

    def test_runner_and_wrapper_helper_decode_identically(self):
        codes = [None] + list(range(-64, 0)) + list(range(0, 260))
        for compact in (False, True):
            for code in codes:
                with self.subTest(code=code, compact=compact):
                    self.assertEqual(runner.describe_exit_code(code, compact=compact),
                                     runtime.describe_exit_status(code, compact=compact))
        self.assertEqual(runner.FAILURE_CAUSE_PREFIX, runtime.FAILURE_CAUSE_PREFIX)
        self.assertEqual([(p, r.pattern) for p, r in runner.CRASH_MESSAGE_PATTERNS],
                         [(p, r.pattern) for p, r in runtime.CRASH_MESSAGE_PATTERNS])

    def test_signal_statuses_are_recognised(self):
        for code in (-11, -9, 129, 139, 192):
            self.assertTrue(runner.exit_code_is_signal(code), code)
        for code in (None, 0, 1, 2, 128, 193, 255):
            self.assertFalse(runner.exit_code_is_signal(code), code)


class ExtractFailureCauseTests(unittest.TestCase):
    def cause(self, body, **kwargs) -> str:
        with tempfile.TemporaryDirectory(prefix="earth-f09-") as tmp:
            return runner.extract_failure_cause(str(runner_log(Path(tmp) / "t.log", body, **kwargs)))

    def test_wrapper_cause_line_wins_over_crash_messages(self):
        report = [SEP, "AMPS FAILURE REPORT: C1", SEP,
                  "FAILURE CAUSE: C1: AMPS exit 139 = 128 + 11 (SIGSEGV: Segmentation fault); " + MPIRUN_LINE]
        cause = self.cause(OPEN_MPI_SEGV + report)
        self.assertTrue(cause.startswith("C1: AMPS exit 139 = 128 + 11 (SIGSEGV"), cause)
        self.assertIn("exited on signal 11", cause)

    def test_launcher_message_naming_the_rank_is_preferred(self):
        self.assertEqual(self.cause(OPEN_MPI_SEGV + ["AMPS command failed with exit code 139"]), MPIRUN_LINE)

    def test_process_message_is_used_without_a_launcher_message(self):
        # The first message of the most specific kind: the failing process's own.
        self.assertEqual(self.cause(OPEN_MPI_SEGV[:8]), OPEN_MPI_SEGV[1])
        self.assertEqual(self.cause(OPEN_MPI_SEGV[:8] + ["rank 3: Killed"]), OPEN_MPI_SEGV[1])

    def test_python_traceback_gives_the_exception_line(self):
        body = ["test_x ... FAIL", "Traceback (most recent call last):",
                '  File "t.py", line 3, in test_x', "    self.assertEqual(a, b)",
                "AssertionError: '! runner: timeout 9000' != '! runner: exclusive'",
                "FAILED (failures=1)", "RESULT: FAIL"]
        self.assertEqual(self.cause(body),
                         "Python exception: AssertionError: '! runner: timeout 9000' != '! runner: exclusive'")

    def test_fail_word_then_last_line_fallbacks(self):
        self.assertEqual(self.cause(["step 1", "C2 FAIL: max |dR|/R = 0.3 > 0.05", "wrote plot"]),
                         "C2 FAIL: max |dR|/R = 0.3 > 0.05")
        # The header is not command output: "--fail-on-paxis" there must not count.
        self.assertEqual(self.cause(["computing"], command="run_C4.py --fail-on-paxis"),
                         "last output: computing")

    def test_launch_error_empty_and_missing_logs(self):
        self.assertEqual(self.cause(["", "ERROR launching command: [Errno 2] No such file"]),
                         "ERROR launching command: [Errno 2] No such file")
        self.assertEqual(self.cause([]), "the command wrote no output")
        self.assertEqual(runner.extract_failure_cause(""), "")
        self.assertEqual(runner.extract_failure_cause("/nonexistent/earth-f09.log"), "")


class RunnerEndToEndTests(unittest.TestCase):
    @staticmethod
    def run_entry(script: str, expected: str = "P"):
        command = "%s -c %s" % (shlex.quote(sys.executable), shlex.quote(script))
        test = runner.TestCase(1, 1, expected, command)
        with tempfile.TemporaryDirectory(prefix="earth-f09-run-") as tmp:
            root = Path(tmp)
            (root / "logs").mkdir()
            plan = runner.build_test_plan(test, default_np=1, default_nt=1)
            return asyncio.run(runner.run_one_test(
                plan, workdir=root, log_dir=root / "logs", timeout_s=30.0, use_shell=False,
                base_env=os.environ.copy(), set_thread_env_vars=False, preserve_thread_env=True,
                print_start=False, mem_total_bytes=None, mem_available_bytes=None,
            ))

    def test_mpirun_style_139_reports_signal_and_rank(self):
        lines = OPEN_MPI_SEGV + ["AMPS command failed with exit code 139"]
        result = self.run_entry("import sys; print(%r); sys.exit(139)" % "\n".join(lines))
        self.assertEqual((result.actual, result.exit_code), ("F", 139))
        self.assertEqual(result.exit_status, "exit 139 = 128 + 11 (SIGSEGV: Segmentation fault)")
        self.assertEqual(result.failure_cause, MPIRUN_LINE)

    def test_command_killed_by_a_signal_is_named(self):
        result = self.run_entry("import os, signal; print('working', flush=True); "
                                "os.kill(os.getpid(), signal.SIGTERM)")
        self.assertEqual((result.actual, result.exit_code), ("F", -15))
        self.assertEqual(result.exit_status, "killed by signal 15 (SIGTERM: Terminated)")
        self.assertEqual(result.failure_cause, "last output: working")

    def test_a_pass_has_no_cause_and_p_f_is_unchanged(self):
        result = self.run_entry("print('ok')")
        self.assertEqual((result.actual, result.matched_reference), ("P", True))
        self.assertEqual((result.exit_status, result.failure_cause), ("exit 0", ""))
        result = self.run_entry("import sys; sys.exit(3)", expected="F")
        self.assertEqual((result.actual, result.matched_reference), ("F", True))
        self.assertEqual(result.failure_cause, "the command wrote no output")


class ReportLabelTests(unittest.TestCase):
    @staticmethod
    def result(index, expected, actual, **kwargs):
        defaults = dict(index=index, line_no=index, expected=expected, actual=actual,
                        matched_reference=(expected == actual), exit_code=None, timed_out=False,
                        elapsed_s=1.0, command="cmd%d" % index, log_file="")
        defaults.update(kwargs)
        return runner.TestResult(**defaults)

    def test_reasons_are_never_unexpected_pass_for_a_failure(self):
        cases = [
            (self.result(1, "F", "F", matched_reference=False, interrupted=True),
             "stopped while running (the runner was interrupted; no result)"),
            (self.result(2, "F", "F", matched_reference=False, not_started=True),
             "not started (the runner was interrupted; no result)"),
            (self.result(3, "P", "F", exit_code=139), "unexpected failure"),
            (self.result(4, "F", "P", exit_code=0), "unexpected pass"),
        ]
        for result, reason in cases:
            with self.subTest(index=result.index):
                self.assertEqual(runner.mismatch_reason(result), reason)

    def test_reports_carry_status_and_cause_and_flag_crashed_expected_f(self):
        results = [
            self.result(1, "P", "F", exit_code=139, exit_status="exit 139 = 128 + 11 (SIGSEGV: Segmentation fault)",
                        failure_cause=MPIRUN_LINE, last_pass="a" * 40),
            self.result(2, "F", "F", exit_code=139, exit_status="exit 139 = 128 + 11 (SIGSEGV: Segmentation fault)",
                        failure_cause="C4 case X: AMPS exit 139"),
            self.result(3, "F", "F", exit_code=1, exit_status="exit 1", failure_cause="C4 FAIL: gate"),
        ]
        with tempfile.TemporaryDirectory(prefix="earth-f09-report-") as tmp:
            paths = runner.write_reports(results, Path(tmp) / "report")
            text = Path(paths[2]).read_text(encoding="utf-8")
            self.assertIn("  status:   exit 139 = 128 + 11 (SIGSEGV: Segmentation fault)\n", text)
            self.assertIn("  cause:    " + MPIRUN_LINE + "\n", text)
            header = Path(paths[1]).read_text(encoding="utf-8").splitlines()[0]
            self.assertTrue(header.endswith(",exit_status,failure_cause"), header)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                runner.print_final_summary(results, json_path=paths[0], csv_path=paths[1],
                                           to_address_txt_path=paths[2], to_address_csv_path=paths[3],
                                           actual_failed_txt_path=paths[4], actual_failed_csv_path=paths[5],
                                           last_pass_failed_txt_path=paths[6], last_pass_failed_csv_path=paths[7],
                                           new_or_unexpected_passed_txt_path=paths[8],
                                           new_or_unexpected_passed_csv_path=paths[9])
        summary = out.getvalue()
        self.assertIn("#001 line 1: unexpected failure", summary)
        self.assertIn("status: exit 139 = 128 + 11 (SIGSEGV: Segmentation fault); cause: " + MPIRUN_LINE, summary)
        crashed = summary.split("Expected-F entries that ended by a signal", 1)
        self.assertEqual(len(crashed), 2, summary)
        self.assertIn("#002 line 2", crashed[1])
        self.assertNotIn("#003 line 3", crashed[1])


class WrapperReportTests(unittest.TestCase):
    def test_report_ends_with_the_cause_and_is_appended_to_the_log(self):
        with tempfile.TemporaryDirectory(prefix="earth-f09-wrap-") as tmp:
            log = Path(tmp) / "C1_amps.log"
            log.write_text("\n".join(["Command:", "  mpirun -np 4 ./amps", ""] + OPEN_MPI_SEGV) + "\n",
                           encoding="utf-8")
            stream = io.StringIO()
            cause = runtime.report_amps_failure(139, ["mpirun", "-np", "4", "./amps"], Path(tmp), log,
                                                label="C1", stream=stream)
            self.assertEqual(stream.getvalue().splitlines()[-1], cause)
            self.assertEqual(log.read_text(encoding="utf-8").splitlines()[-1], cause)
        self.assertEqual(cause, "FAILURE CAUSE: C1: AMPS exit 139 = 128 + 11 (SIGSEGV: Segmentation fault); "
                         + MPIRUN_LINE)
        report = stream.getvalue()
        for token in ("status:   exit status 139 = 128 + 11", "stack size limit (ulimit -s):",
                      "OMP_NUM_THREADS=", "L10: [node:101] Signal: Segmentation fault (11)",
                      "last 17 lines of the AMPS log:"):
            self.assertIn(token, report)

    def test_missing_log_and_log_without_crash_message(self):
        with tempfile.TemporaryDirectory(prefix="earth-f09-wrap-") as tmp:
            missing = Path(tmp) / "absent" / "amps.log"
            cause = runtime.report_amps_failure(1, ["amps"], Path(tmp), missing, label="C2",
                                                stream=io.StringIO())
            self.assertEqual(cause, "FAILURE CAUSE: C2: AMPS exit 1; the AMPS log %s does not exist" % missing)
            plain = Path(tmp) / "amps.log"
            plain.write_text("reading input\n", encoding="utf-8")
            cause = runtime.report_amps_failure(2, ["amps"], Path(tmp), plain, label="C3/UPPER_SCAN",
                                                stream=io.StringIO())
            self.assertEqual(cause, "FAILURE CAUSE: C3/UPPER_SCAN: AMPS exit 2; no crash message "
                             "(signal, rank, backtrace, abort) in the AMPS log")


class CWrapperWiringTests(unittest.TestCase):
    WRAPPERS = ("C1", "C2", "C3", "C4")

    def test_every_amps_failure_branch_reports(self):
        for case in self.WRAPPERS:
            source = (TEST_ROOT / case / ("run_%s.py" % case)).read_text(encoding="utf-8")
            with self.subTest(case=case):
                self.assertGreater(source.count("rc = run_command("), 0)
                self.assertEqual(source.count("rc = run_command("), source.count("report_amps_failure(rc"))
                self.assertNotIn("AMPS command failed with exit code", source)

    def test_c1_with_a_crashing_launcher_returns_139_and_names_the_cause(self):
        with tempfile.TemporaryDirectory(prefix="earth-f09-c1-") as tmp:
            root = Path(tmp)
            mpirun = root / "fake_mpirun"
            mpirun.write_text("#!/bin/sh\ncat <<'EOF'\n%s\nEOF\nexit 139\n" % "\n".join(OPEN_MPI_SEGV),
                              encoding="utf-8")
            mpirun.chmod(0o755)
            completed = subprocess.run(
                [sys.executable, str(TEST_ROOT / "C1" / "run_C1.py"), "--mode", "gridless",
                 "-np", "4", "-nt", "1", "--mpirun", str(mpirun), "--amps", str(mpirun),
                 "--workdir", str(root / "work")],
                cwd=str(root), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120,
            )
        self.assertEqual(completed.returncode, 139, completed.stdout[-2000:])
        last = completed.stdout.splitlines()[-1]
        self.assertEqual(last, "FAILURE CAUSE: C1: AMPS exit 139 = 128 + 11 (SIGSEGV: Segmentation fault); "
                         + MPIRUN_LINE)
        # AMPS output (stdout) precedes the report (stderr) in the merged stream.
        self.assertLess(completed.stdout.index("Caught signal 11"),
                        completed.stdout.index("AMPS FAILURE REPORT: C1"))


if __name__ == "__main__":
    unittest.main(verbosity=2)

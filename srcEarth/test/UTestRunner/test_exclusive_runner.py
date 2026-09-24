#!/usr/bin/env python3
"""Dependency-free regressions for test_runner exclusive scheduling.

The production runner is imported directly from the adjacent test directory.
No AMPS executable, MPI runtime, or third-party Python module is needed.  The
scheduler test replaces only the subprocess coroutine with a timed in-memory
stand-in; it therefore exercises the real dispatcher state machine without
starting artificial shell processes or depending on host load.
"""

from __future__ import annotations

import asyncio
import contextlib
import io
import importlib.util
import os
import shlex
import sys
import tempfile
import time
import unittest
from pathlib import Path


RUNNER_PATH = Path(__file__).resolve().parents[1] / "test_runner.py"
TEST_LIST_PATH = RUNNER_PATH.with_name("list")
SPEC = importlib.util.spec_from_file_location("earth_test_runner", RUNNER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load production test_runner.py")
runner = importlib.util.module_from_spec(SPEC)
# dataclasses resolves the defining module through sys.modules while classes are
# created, so register the temporary import name before executing the module.
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def parse_text(text: str):
    """Parse one temporary list using the production parser."""

    with tempfile.TemporaryDirectory(prefix="earth-runner-parse-") as tmp:
        path = Path(tmp) / "list"
        path.write_text(text, encoding="utf-8")
        return runner.parse_test_file(path)


class ExclusiveParserTests(unittest.TestCase):
    """UTR-F01: directive parsing is strict and command-preserving."""

    def test_scalar_directive_preserves_command_and_provenance(self):
        command = "python3 -c \"print('unchanged')\""
        tests = parse_text(
            "! runner: exclusive\n"
            f"P {command}\n"
            "last pass: old-commit\n"
        )
        self.assertEqual(len(tests), 1)
        self.assertTrue(tests[0].exclusive)
        self.assertEqual(tests[0].command, command)
        self.assertEqual(tests[0].expected, "P")
        self.assertEqual(tests[0].last_pass, "old-commit")

    def test_directive_before_loop_applies_to_every_expansion(self):
        tests = parse_text(
            "! runner: exclusive\n"
            "for $m={RK4,HC4}\n"
            "! comments remain legal between a loop and its command\n"
            "{P,F} command --mover $m\n"
            "last pass: {rk4-pass,hc4-pass}\n"
        )
        self.assertEqual([test.command for test in tests], [
            "command --mover RK4", "command --mover HC4"])
        self.assertEqual([test.expected for test in tests], ["P", "F"])
        self.assertTrue(all(test.exclusive for test in tests))
        self.assertEqual(
            [test.last_pass for test in tests], ["rk4-pass", "hc4-pass"])

    def test_directive_may_immediately_precede_loop_command(self):
        tests = parse_text(
            "for $m={RK4,HC4}\n"
            "! runner: exclusive\n"
            "P command --mover $m\n"
        )
        self.assertEqual(len(tests), 2)
        self.assertTrue(all(test.exclusive for test in tests))

    def test_ordinary_entry_remains_nonexclusive(self):
        tests = parse_text("! ordinary comment\nP command --value 1\n")
        self.assertEqual(len(tests), 1)
        self.assertFalse(tests[0].exclusive)

    def test_invalid_directives_fail_loudly(self):
        with self.assertRaisesRegex(ValueError, "unsupported runner directive"):
            parse_text("! runner: exclusve\nP command\n")
        with self.assertRaisesRegex(ValueError, "no following test command"):
            parse_text("! runner: exclusive\n")
        with self.assertRaisesRegex(ValueError, "comment separates"):
            parse_text(
                "! runner: exclusive\n"
                "! this must not silently retarget the directive\n"
                "P command\n"
            )

    def test_last_pass_update_retains_directive_and_command(self):
        source = (
            "! runner: exclusive\n"
            "P command --scientific-gate strict\n"
            "last pass: old-commit\n"
        )
        with tempfile.TemporaryDirectory(prefix="earth-runner-update-") as tmp:
            path = Path(tmp) / "list"
            path.write_text(source, encoding="utf-8")
            tests = runner.parse_test_file(path)
            result = runner.TestResult(
                index=tests[0].index,
                line_no=tests[0].line_no,
                expected="P",
                actual="P",
                matched_reference=True,
                exit_code=0,
                timed_out=False,
                elapsed_s=0.0,
                command=tests[0].command,
                log_file="unused.log",
                last_pass=tests[0].last_pass,
                last_pass_line_no=tests[0].last_pass_line_no,
                exclusive=True,
            )
            changed = runner.update_last_pass_entries(
                path, tests, [result], "new-commit")
            self.assertEqual(changed, 1)
            self.assertEqual(
                path.read_text(encoding="utf-8"),
                source.replace("old-commit", "new-commit"),
            )


class ExclusiveSchedulerTests(unittest.TestCase):
    """UTR-F02: exclusive work drains and owns the runner pool."""

    def test_exclusive_plan_never_overlaps_neighboring_plans(self):
        tests = [
            runner.TestCase(1, 1, "P", "ordinary-a"),
            runner.TestCase(2, 2, "P", "ordinary-b"),
            runner.TestCase(3, 3, "P", "exclusive-x", exclusive=True),
            runner.TestCase(4, 4, "P", "ordinary-c"),
        ]
        delays = {
            "ordinary-a": 0.030,
            "ordinary-b": 0.050,
            "exclusive-x": 0.020,
            "ordinary-c": 0.010,
        }
        starts = {}
        finishes = {}

        async def fake_run_one_test(plan, **_kwargs):
            command = plan.test.command
            starts[command] = time.monotonic()
            await asyncio.sleep(delays[command])
            finishes[command] = time.monotonic()
            return runner.TestResult(
                index=plan.test.index,
                line_no=plan.test.line_no,
                expected=plan.test.expected,
                actual="P",
                matched_reference=True,
                exit_code=0,
                timed_out=False,
                elapsed_s=delays[command],
                command=command,
                log_file="unused.log",
                exclusive=plan.test.exclusive,
            )

        original = runner.run_one_test
        runner.run_one_test = fake_run_one_test
        try:
            with tempfile.TemporaryDirectory(prefix="earth-runner-schedule-") as tmp:
                # Completion lines are useful in the production runner but only
                # add noise to this unit-test transcript.
                with contextlib.redirect_stdout(io.StringIO()):
                    results = asyncio.run(runner.run_all_tests(
                        tests,
                        jobs=3,
                        workdir=Path(tmp),
                        log_dir=Path(tmp) / "logs",
                        timeout_s=None,
                        use_shell=False,
                        default_np=1,
                        default_nt=1,
                        set_thread_env_vars=False,
                        preserve_thread_env=True,
                        print_start=False,
                        memory_gate_enabled=False,
                        min_free_memory_fraction=0.5,
                        memory_check_interval=0.0,
                    ))
        finally:
            runner.run_one_test = original

        self.assertEqual([result.index for result in results], [1, 2, 3, 4])
        # The first two ordinary tests should really overlap; otherwise this
        # test could pass accidentally under a fully serial dispatcher.
        self.assertLess(starts["ordinary-b"], finishes["ordinary-a"])
        # Exclusive X cannot start until both earlier jobs have drained.
        self.assertGreaterEqual(
            starts["exclusive-x"],
            max(finishes["ordinary-a"], finishes["ordinary-b"]),
        )
        # No later work may start while exclusive X owns the pool.
        self.assertGreaterEqual(starts["ordinary-c"], finishes["exclusive-x"])


class ExclusiveEnvironmentTests(unittest.TestCase):
    """UTR-F03: exclusive provenance reaches only the protected child."""

    def test_run_one_test_propagates_and_scrubs_marker(self):
        marker = runner.RUNNER_EXCLUSIVE_ENV
        expected = runner.RUNNER_EXCLUSIVE_ENV_VALUE

        def marker_command(wanted):
            check = (
                "import os,sys; "
                "sys.exit(0 if os.environ.get(%r) == %r else 19)" %
                (marker, wanted))
            return "%s -c %s" % (
                shlex.quote(sys.executable), shlex.quote(check))

        async def execute(plan, root, base_env):
            return await runner.run_one_test(
                plan,
                workdir=root,
                log_dir=root / "logs",
                timeout_s=10.0,
                use_shell=False,
                base_env=base_env,
                set_thread_env_vars=False,
                preserve_thread_env=True,
                print_start=False,
                mem_total_bytes=None,
                mem_available_bytes=None,
            )

        with tempfile.TemporaryDirectory(prefix="earth-runner-env-") as tmp:
            root = Path(tmp)
            (root / "logs").mkdir()
            # Seed a bogus parent value.  The runner must overwrite it for the
            # exclusive child and remove it completely for the ordinary child.
            base_env = os.environ.copy()
            base_env[marker] = "stale-parent-value"

            exclusive_test = runner.TestCase(
                1, 1, "P", marker_command(expected), exclusive=True)
            exclusive_plan = runner.build_test_plan(
                exclusive_test, default_np=1, default_nt=1)
            exclusive_result = asyncio.run(
                execute(exclusive_plan, root, base_env))
            self.assertEqual(exclusive_result.exit_code, 0)
            exclusive_log = Path(exclusive_result.log_file).read_text(
                encoding="utf-8")
            self.assertIn("# Exclusive scheduling: yes", exclusive_log)
            self.assertIn("%s=%s" % (marker, expected), exclusive_log)

            ordinary_test = runner.TestCase(
                2, 2, "P", marker_command(None), exclusive=False)
            ordinary_plan = runner.build_test_plan(
                ordinary_test, default_np=1, default_nt=1)
            ordinary_result = asyncio.run(execute(ordinary_plan, root, base_env))
            self.assertEqual(ordinary_result.exit_code, 0)
            ordinary_log = Path(ordinary_result.log_file).read_text(
                encoding="utf-8")
            self.assertIn("# Exclusive scheduling: no", ordinary_log)
            self.assertIn("%s=<unset>" % marker, ordinary_log)

        # The shared base environment is never mutated by either launch.
        self.assertEqual(base_env[marker], "stale-parent-value")


class ActiveC19EntryTests(unittest.TestCase):
    """UTR-F04: C19 has isolation provenance and bounded CPU fan-out."""

    def test_active_c19_entry_requires_exclusive_bounded_parallelism(self):
        tests = runner.parse_test_file(TEST_LIST_PATH)
        matches = [
            test for test in tests
            if "srcEarth/test/C19/run_C19.py" in test.command
            and "--require-runner-exclusive" in test.command
        ]
        self.assertEqual(len(matches), 1)
        test = matches[0]
        self.assertTrue(test.exclusive)

        argv = shlex.split(test.command)
        self.assertEqual(argv.count("--mode3d-parallel-field-init"), 1)
        self.assertEqual(argv.count("--require-runner-exclusive"), 1)
        self.assertEqual(int(argv[argv.index("-np") + 1]), 4)
        self.assertEqual(int(argv[argv.index("-nt") + 1]), 8)

        # Mode3D field initialization uses nt temporary workers plus the
        # original caller on every MPI rank: 4*(8+1)=36 participants, matching
        # the 36-CPU allocation in the failing production log instead of the
        # previous oversubscribed 4*(33+1)=136 participants.
        np_value = int(argv[argv.index("-np") + 1])
        nt_value = int(argv[argv.index("-nt") + 1])
        self.assertEqual(np_value * (nt_value + 1), 36)


if __name__ == "__main__":
    unittest.main(verbosity=2)

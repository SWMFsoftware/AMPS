#!/usr/bin/env python3
"""Self-test for tools/sep_test_orchestrator.py (run_all_test.sh backend).

Each result reader is exercised on a small synthetic fixture written in the
exact format its runner produces (JSON report fields, shared-log markers,
unittest verbose text, ``[ID] STATUS`` lines).  Expected counts, statuses,
messages and extracted log contents are written into the fixtures and
asserted literally; nothing is computed by the code under test.  The
step-level guards (missing report, nonzero exit without a failing record),
the status vocabulary, de-duplication into "unique tests", and the summary's
exit status and report files are covered as well.

Run:  python3 tools/test_sep_test_orchestrator.py      (verbosity 2)
Needs no AMPS build; it is registered as the ``orchestrator-selftest`` step.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_SPEC = importlib.util.spec_from_file_location("sep_test_orchestrator",
                                               _HERE / "sep_test_orchestrator.py")
orch = importlib.util.module_from_spec(_SPEC)
# dataclasses resolves string annotations through sys.modules, so the module
# must be registered before it executes.
sys.modules[_SPEC.name] = orch
_SPEC.loader.exec_module(orch)


class Fixture(unittest.TestCase):
    """Temporary output directory plus a context factory for one step."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def ctx(self, name="step", log_text="", exit_code=0, kind="test",
            runner="runner.py", namespace=""):
        step = orch.Step(name, "fixture step", [], kind=kind, runner=runner,
                         namespace=namespace)
        log = self.out / f"01-{name}.log"
        log.write_text(log_text, encoding="utf-8")
        return orch.StepContext(step, self.out, log, exit_code, time.time() - 5,
                                self.out / "failed-tests" / name)

    def write_json(self, relative, payload):
        path = self.out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path


class StatusVocabulary(unittest.TestCase):
    def test_runner_spellings_map_to_four_statuses(self):
        cases = {"PASS": "PASS", "passed": "PASS", "ok": "PASS", "FAIL": "FAIL",
                 "failed": "FAIL", "SKIP": "SKIP", "skipped": "SKIP",
                 "ERROR": "ERROR", "errors": "ERROR", "MISSING": "ERROR",
                 "weird": "ERROR"}
        for raw, expected in cases.items():
            self.assertEqual(orch.normalize_status(raw), expected, raw)


class SrcSepRunTests(Fixture):
    def test_failed_test_section_is_copied_from_shared_log(self):
        report = self.out / "run-tests"
        self.write_json("run-tests/srcsep-tests.json", {"results": [
            {"id": "PARK01", "status": "PASS", "message": "ok"},
            {"id": "FTE01", "status": "ERROR", "message": "process exited 0"},
            {"id": "OV02", "status": "ERROR", "message": "TypeError gating"},
            {"id": "SCAT01", "status": "SKIP", "message": "no rule"}]})
        self.write_json("run-tests/individual/OV02/OV02.json", {})
        (report / "test-run.log").write_text(
            "header\n"
            "[2026-10-10T06:17:00+00:00] RUN /amps --test PARK01 --test-json x\n"
            "park output\n"
            "[2026-10-10T06:17:04+00:00] RUN /amps --test FTE01 --test-json y\n"
            "fte line 1\nAMPS: exit: focused_transport_dmumu.cpp:13\n"
            "[2026-10-10T06:18:00+00:00] RUN python3 run_case.py --case OV02 --output-dir z\n"
            "Traceback: TypeError gating\n", encoding="utf-8")
        records = orch.reader_srcsep_run_tests(report)(self.ctx("run-tests"))
        by_id = {r.test_id: r for r in records}
        self.assertEqual(orch.count(records),
                         {"PASS": 1, "FAIL": 0, "SKIP": 1, "ERROR": 2, "TOTAL": 4})
        fte = Path(by_id["FTE01"].log).read_text()
        self.assertIn("# source: ", fte)
        self.assertIn("focused_transport_dmumu.cpp:13", fte)
        self.assertNotIn("park output", fte)
        self.assertNotIn("TypeError", fte)
        self.assertIn("--case OV02", Path(by_id["OV02"].log).read_text())
        self.assertTrue(by_id["OV02"].evidence[0].endswith("individual/OV02/OV02.json"))
        self.assertEqual(by_id["PARK01"].log, str(report / "test-run.log"))

    def test_missing_report_is_unavailable(self):
        with self.assertRaises(orch.ResultsUnavailable):
            orch.reader_srcsep_run_tests(self.out / "absent")(self.ctx())


class SrcSep3dRunTests(Fixture):
    def test_runner_log_is_used_and_fallback_locates_line(self):
        report = self.out / "rt"
        self.write_json("rt/srcsep3d-tests.json", {"results": [
            {"test_id": "R3D01", "status": "FAIL", "message": "deck", "log": "/logs/R3D01.log"},
            {"test_id": "LAY01", "status": "ERROR", "message": "crash"}]})
        ctx = self.ctx("run-tests-all", "[R3D01] FAIL\nnoise\n[LAY01] ERROR (0.1s) crash\n")
        by_id = {r.test_id: r for r in orch.reader_srcsep3d_run_tests(report)(ctx)}
        self.assertEqual(by_id["R3D01"].log, "/logs/R3D01.log")
        self.assertEqual((by_id["LAY01"].log, by_id["LAY01"].line), (str(ctx.step_log), 3))


class PdMover(Fixture):
    def test_missing_maps_to_error_and_suite_namespaces(self):
        report = self.write_json("pd.json", {"suites": [
            {"name": "library", "results": [
                {"test_id": "PD13-HOST-GATE", "status": "PASS", "suite": "library",
                 "detail": "ok", "log": "/l/library.log", "line": 9}]},
            {"name": "stage1", "results": [
                {"test_id": "CFG3D17", "status": "MISSING", "suite": "stage1",
                 "detail": "expected test did not report", "log": "/l/stage1.log", "line": 0}]}]})
        records = orch.reader_pd_mover(report, "pd.py", {"library": "lib"})(self.ctx("pd-mover"))
        self.assertEqual([r.status for r in records], ["PASS", "ERROR"])
        self.assertEqual(records[0].namespace, "lib")
        self.assertEqual(records[1].namespace, "pd-mover/stage1")
        self.assertEqual(records[0].runner, "pd.py [library]")
        self.assertEqual((records[0].log, records[0].line), ("/l/library.log", 9))


class NativeJson(Fixture):
    def test_results_point_to_step_log_line_and_artifacts(self):
        report = self.write_json("n/native.json", {"results": [
            {"id": "SCCM3D04", "status": "FAIL", "message": "epoch", "artifacts": ["/a/x"]}]})
        ctx = self.ctx("sccm3d", "init\n[SCCM3D04] FAIL epoch\n")
        record = orch.reader_native_json(report)(ctx)[0]
        self.assertEqual((record.status, record.line), ("FAIL", 2))
        self.assertEqual(record.evidence, [str(report), "/a/x"])


class Coupled(Fixture):
    def test_existing_diagnostic_used_else_output_copied(self):
        run_dir = self.out / "cr" / "runs" / "r1"
        diag = run_dir / "failures" / "native-amps" / "N01.log"
        diag.parent.mkdir(parents=True)
        diag.write_text("diagnostic")
        self.write_json("cr/summary.json", {"run_directory": str(run_dir), "results": [
            {"id": "N01", "scope": "native-amps", "status": "FAIL", "message": "m1"},
            {"id": "S01", "scope": "shared-model", "status": "ERROR", "message": "m2",
             "output": "captured traceback"},
            {"id": "S02", "scope": "shared-model", "status": "SKIP", "message": ""}]})
        records = {r.test_id: r for r in orch.reader_coupled(self.out / "cr")(self.ctx("coupled-runner"))}
        self.assertEqual(records["native-amps/N01"].log, str(diag))
        self.assertIn("captured traceback", Path(records["shared-model/S01"].log).read_text())
        self.assertEqual(records["shared-model/S02"].status, "SKIP")


class ReducedAndPositive(Fixture):
    def test_reduced_front_and_positive_shock_logs(self):
        self.write_json("rf/summary.json", {"results": [
            {"scope": "native", "test_id": "RSF07", "status": "FAIL", "message": "t",
             "log": "/rf/native.log", "report": "/rf/r.json"}]})
        self.write_json("ps/summary.json", {"checks": [
            {"run": "rank-4", "name": "arrival", "status": "FAIL", "detail": "late",
             "log": "/ps/execution.log"}]})
        rf = orch.reader_reduced_front(self.out / "rf")(self.ctx("reduced-front"))[0]
        ps = orch.reader_positive_shock(self.out / "ps")(self.ctx("positive-shock"))[0]
        self.assertEqual((rf.test_id, rf.log, rf.evidence), ("native/RSF07", "/rf/native.log", ["/rf/r.json"]))
        self.assertEqual((ps.test_id, ps.message, ps.log), ("rank-4/arrival", "late", "/ps/execution.log"))


class CoronalCme(Fixture):
    def test_fresh_report_copies_output_and_stale_report_is_rejected(self):
        report = self.write_json("cc/results.json", {"tests": [
            {"id": "SCCM01", "status": "FAIL", "message": "m", "output": "assert x"},
            {"id": "SCCM02", "status": "SKIP", "message": ""}]})
        records = orch.reader_coronal_cme(report)(self.ctx("coronal-cme"))
        self.assertIn("assert x", Path(records[0].log).read_text())
        old = time.time() - 3600
        os.utime(report, (old, old))
        with self.assertRaises(orch.ResultsUnavailable):
            orch.reader_coronal_cme(report)(self.ctx("coronal-cme"))


class UnittestText(Fixture):
    def test_verbose_lines_docstrings_and_failure_blocks(self):
        text = ("test_a (test_mod.Viewer) ... ok\n"
                "test_b (test_mod.Viewer)\nDocstring first line. ... FAIL\n"
                "test_c (test_mod.Viewer) ... skipped 'no display'\n"
                "\n" + "=" * 70 + "\n"
                "FAIL: test_b (test_mod.Viewer)\n" + "-" * 70 + "\n"
                "Traceback (most recent call last):\nAssertionError: 3 != 4\n\n"
                + "-" * 70 + "\nRan 3 tests in 0.1s\n\nFAILED (failures=1, skipped=1)\n")
        records = {r.test_id: r for r in orch.read_unittest(self.ctx("viewer-tests", text))}
        self.assertEqual({k: r.status for k, r in records.items()},
                         {"Viewer.test_a": "PASS", "Viewer.test_b": "FAIL",
                          "Viewer.test_c": "SKIP"})
        block = Path(records["Viewer.test_b"].log).read_text()
        self.assertIn("FAIL: test_b", block)
        self.assertIn("# source:", block)

    def test_no_result_lines_is_unavailable(self):
        with self.assertRaises(orch.ResultsUnavailable):
            orch.read_unittest(self.ctx("viewer-tests", "ImportError: x\n"))


class BracketText(Fixture):
    def test_bracket_and_named_verdicts_worst_status_wins(self):
        text = ("[CMBGU01] PASS config\n[CMBGU05-CONTACT] flux_bound=1e-17\n"
                "BG3D-4 FAIL: leakage exceeds bound\n[ARCH01] PASS ok\n"
                "[CMBGU01] FAIL later\n")
        records = {r.test_id: r for r in orch.read_bracket_text(self.ctx("corona-swcme", text))}
        self.assertEqual({k: r.status for k, r in records.items()},
                         {"CMBGU01": "FAIL", "BG3D-4": "FAIL", "ARCH01": "PASS"})
        self.assertEqual(records["BG3D-4"].line, 3)
        self.assertEqual(records["CMBGU01"].line, 5)


class MakeRecipe(Fixture):
    def test_every_binary_counted_and_hidden_failures_recovered(self):
        text = ("[ARCH01] PASS headers\n"
                "build/test_a\n[A1] PASS one\nBG3D-4 FAIL: leakage\n"
                "make: [makefile:122: test] Error 1 (ignored)\n"
                "build/test_b\nmetric=1\n"
                "build/test_c\nsegfault\n"
                "make: [makefile:122: test] Error 139 (ignored)\n"
                "build/test_d\n[D1] PASS ok\n"
                "make: [makefile:122: test] Error 2 (ignored)\n")
        reader = orch.reader_make_recipe(r"^build/(test_\w+)$")
        records = {r.test_id: r for r in reader(self.ctx("corona-swcme", text))}
        self.assertEqual({k: r.status for k, r in records.items()},
                         {"ARCH01": "PASS", "A1": "PASS", "BG3D-4": "FAIL",
                          "test_b": "PASS", "test_c": "ERROR", "D1": "PASS",
                          "test_d": "ERROR"})
        self.assertEqual(records["test_c"].message, "exited with status 139")
        self.assertEqual(records["test_c"].line, 8)


class StepGuards(Fixture):
    def run_step(self, command, reader):
        step = orch.Step("guard", "fixture", [command], reader=reader, runner="r")
        args = argparse.Namespace(stream=False)
        with contextlib.redirect_stdout(io.StringIO()):
            return orch.execute(step, 1, self.out, orch.Console(None), args)

    def test_nonzero_exit_without_failing_record_is_an_error(self):
        run = self.run_step(["sh", "-c", "exit 3"],
                            lambda ctx: [ctx.record("T1", "PASS")])
        self.assertEqual([r.test_id for r in run.records], ["T1", "guard:exit-status"])
        self.assertEqual(run.status, "FAIL")

    def test_missing_report_becomes_results_unavailable(self):
        def reader(ctx):
            raise orch.ResultsUnavailable("report not found: x")
        run = self.run_step(["sh", "-c", "echo hi"], reader)
        self.assertEqual([(r.test_id, r.status) for r in run.records],
                         [("guard:results-unavailable", "ERROR")])
        self.assertIn("hi", Path(run.log).read_text())


class Summary(Fixture):
    def test_counts_unique_dedupe_files_and_exit_status(self):
        def record(test_id, status, step, namespace, kind="test"):
            return orch.TestRecord(test_id, status, "msg", step, "r", log="/log",
                                   namespace=namespace, kind=kind)
        runs = [
            orch.StepRun(orch.Step("build", "", [], kind="check"), "PASS", 0, 1.0, "",
                         [record("build", "PASS", "build", "build", "check")]),
            orch.StepRun(orch.Step("a", "", []), "FAIL", 1, 1.0, "",
                         [record("X1", "PASS", "a", "shared"), record("X2", "FAIL", "a", "shared")]),
            orch.StepRun(orch.Step("b", "", []), "PASS", 0, 1.0, "",
                         [record("X1", "PASS", "b", "shared"), record("Y1", "SKIP", "b", "b")]),
        ]
        args = argparse.Namespace(dry_run=False, junit=self.out / "j.xml")
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = orch.summarize("App", runs, orch.Console(None), self.out, args, "t0")
        text = buffer.getvalue()
        self.assertEqual(code, 1)
        self.assertIn("test executions: total=4 pass=2 fail=1 skip=1 error=0", text)
        self.assertIn("unique tests:    total=3 pass=1 fail=1 skip=1 error=0", text)
        self.assertIn("checks:          total=1 pass=1 fail=0", text)
        self.assertIn("FAIL   X2", text)
        report = json.loads((self.out / "run_all_test.json").read_text())
        self.assertEqual(report["counts"]["unique_tests"]["TOTAL"], 3)
        self.assertIn("X2", (self.out / "failed_tests.txt").read_text())
        self.assertIn('failures="1"', (self.out / "j.xml").read_text())

    def test_all_passing_returns_zero(self):
        runs = [orch.StepRun(orch.Step("a", "", []), "PASS", 0, 1.0, "",
                             [orch.TestRecord("X1", "PASS", "", "a", "r", namespace="a")])]
        args = argparse.Namespace(dry_run=False, junit=None)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(orch.summarize("App", runs, orch.Console(None), self.out, args, "t0"), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)

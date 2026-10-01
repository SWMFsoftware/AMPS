#!/usr/bin/env python3
"""Process-boundary regressions for the aggregate runner, never MPI evidence.

The tiny executable/launcher/shared fixtures below exercise discovery, argv,
report completeness and failure classification. They deliberately do not
allocate AMPS cells or simulate physical/MPI qualification. Published test
fixtures have synthetic provider names and cannot establish coronal runtime
qualification; the real shared-model suite is exercised separately.
"""
from __future__ import annotations
import json
import contextlib
import importlib.util
import io
from pathlib import Path
import queue
import subprocess
import sys
import tempfile
import threading
import time
import unittest
import xml.etree.ElementTree as ET

RUNNER = Path(__file__).with_name("run_coupled_sep_corona.py")


class CoupledRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="runner-mechanics-only-")
        self.root = Path(self.temporary.name)
        self.model = self.root/"model"
        (self.model/"test").mkdir(parents=True)
        (self.model/"makefile").write_text("build/sep_coronal_cme_tests check-adapters:\n\t@true\n")
        self.fixture = {"shared": [{"id": "CFGDEMO01", "stage": 0, "passed": True, "seconds": 0.01, "output": "synthetic fixture"},
                                   {"id": "PROVDEMO01", "stage": 12, "passed": True, "seconds": 0.01, "output": "synthetic fixture"}],
                        "native": [{"id": "HOSTDEMO01", "status": "PASS", "metrics": [], "artifacts": [], "message": "synthetic fixture"},
                                   {"id": "FUTURECHECK9", "status": "PASS", "metrics": [], "artifacts": [], "message": "synthetic fixture"}],
                        "mode": "complete"}
        self.config = self.root/"fixture.json"
        self.save()
        self.deck = self.root/"synthetic.in"; self.deck.write_text("runner-mechanics-only\n")
        # Each fixture returns the current registry separately from its report;
        # mutations below can invalidate either without changing the driver.
        (self.model/"test/run_tests.py").write_text('''from pathlib import Path
import json,sys
config=json.loads(Path(__file__).resolve().parents[2].joinpath("fixture.json").read_text())
rows=config["shared"]
if "--list" in sys.argv:
    for row in rows: print(row["id"]+"\\tstage="+str(row["stage"])+"\\tpreprocessing")
    raise SystemExit(0)
output=Path(sys.argv[sys.argv.index("--output-dir")+1]); output.mkdir(parents=True)
reported=rows[:-1] if config["mode"] == "shared-missing" else rows
for index,row in enumerate(reported):
    print("["+row["id"]+"] "+("PASS" if row["passed"] else "FAIL"))
    if config["mode"]=="slow-progress" and index==0:
        # No explicit flush: the driver must make this visible before exit.
        import time
        release=Path(__file__).resolve().parents[2]/"release-progress"
        deadline=time.monotonic()+8
        while not release.exists() and time.monotonic()<deadline: time.sleep(0.02)
passed=sum(r["passed"] for r in reported)
(output/"results.json").write_text(json.dumps(dict(suite="sep_coronal_cme",total=len(reported),passed=passed,failed=len(reported)-passed,tests=reported)))
raise SystemExit(0 if passed==len(reported) else 1)
''')
        self.amps = self.root/"fixture-amps"
        self.amps.write_text("#!"+sys.executable+'''\nfrom pathlib import Path
import json,sys
config=json.loads(Path(__file__).with_name("fixture.json").read_text())
rows=config["native"]
if "--list-tests" in sys.argv:
    for row in rows: print(row["id"]+" | synthetic | mechanics only | suite=sep-corona")
    print("OTHER01 | unrelated | mechanics only | suite=amps")
    raise SystemExit(0)
assert sys.argv[sys.argv.index("--test-suite")+1]=="sep-corona" and "--test" not in sys.argv
Path(__file__).with_name("native-argv.json").write_text(json.dumps(sys.argv))
if config["mode"]=="missing-report": raise SystemExit(0)
if config["mode"]=="native-missing": rows=rows[:-1]
if config["mode"]=="native-duplicate": rows=rows+[rows[0]]
if config["mode"]=="native-extra": rows=rows+[dict(rows[0],id="EXTRA01")]
code=2 if any(r["status"]=="ERROR" for r in rows) else (1 if any(r["status"]=="FAIL" for r in rows) else 0)
ranks=int(sys.argv[sys.argv.index("--expect-mpi-ranks")+1])
data=dict(schema="srcsep-component-tests-v1",configuration_fingerprint="synthetic-mechanics-only",mpi_ranks=ranks,
          completed_steps=0,providers=dict(input_schema_version=4,background="synthetic-fixture",shock="synthetic-fixture"),results=rows)
if config["mode"]=="rank-mismatch": data["mpi_ranks"]=1
if config["mode"]=="obsolete-provider-report": del data["providers"]
report=Path(sys.argv[sys.argv.index("--test-json")+1]);report.parent.mkdir(parents=True)
report.write_text(json.dumps(data))
for row in rows: print("["+row["id"]+"] "+row["status"]+" - synthetic fixture")
raise SystemExit(0 if config["mode"]=="exit-mismatch" else code)
''')
        self.amps.chmod(0o755)
        self.launcher = self.root/"fixture-launcher.py"
        self.launcher.write_text("import subprocess,sys\nassert int(sys.argv[1])==3\nraise SystemExit(subprocess.run(sys.argv[2:]).returncode)\n")
        self.output = self.root/"output"

    def tearDown(self):
        self.temporary.cleanup()

    def save(self):
        self.config.write_text(json.dumps(self.fixture))

    def run_suite(self, *arguments):
        process = subprocess.run([sys.executable, str(RUNNER), "--model-dir", str(self.model),
            "--amps", str(self.amps), "--test-input", str(self.deck), "--ranks", "3",
            "--launcher", sys.executable+" "+str(self.launcher)+" {ranks}", "--output-dir", str(self.output),
            *arguments], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        report = json.loads((self.output/"summary.json").read_text())
        return process, report

    def test_complete_two_scopes_and_dynamic_future_registry(self):
        process, report = self.run_suite()
        self.assertEqual(process.returncode, 0, process.stdout)
        self.assertEqual(report["counts"], {"total": 4, "pass": 4, "fail": 0, "skip": 0, "error": 0})
        self.assertEqual({r["scope"] for r in report["results"]}, {"shared-model", "native-amps"})
        self.assertIn("not-established", report["coronal_runtime_provider_qualification"])
        self.assertIn("failed_test_summary: fail=0 error=0 runner_errors=0", process.stdout)
        self.assertIn("test_logs_directory="+report["run_directory"], process.stdout)
        argv = json.loads((self.root/"native-argv.json").read_text())
        self.assertIn("--test-suite", argv); self.assertNotIn("--test", argv)
        self.fixture["native"].append(dict(self.fixture["native"][0], id="UNNAMEDNEXT88"))
        self.save(); process, report = self.run_suite()
        self.assertEqual(process.returncode, 0, process.stdout)
        self.assertEqual(report["counts"]["total"], 5)
        self.assertIn("UNNAMEDNEXT88", [r["id"] for r in report["results"]])

    def test_mixed_native_pass_fail_reconciles_suite_exit(self):
        self.fixture["native"][1]["status"] = "FAIL"; self.save()
        process, report = self.run_suite()
        self.assertEqual(process.returncode, 1, process.stdout)
        self.assertEqual(report["counts"]["pass"], 3)
        self.assertEqual(report["counts"]["fail"], 1)
        self.assertEqual(report["counts"]["error"], 0)

    def test_failed_case_summary_retains_details_and_fresh_log_paths(self):
        # A traceback may be too long for a console summary. Its complete body
        # belongs to the per-case diagnostic, with the final line as the reason.
        captured = "debug context\n"+"detailed output "*150+"\nAssertionError: synthetic shared mismatch\n"
        self.fixture["shared"][0].update(passed=False, output=captured)
        self.fixture["native"][1].update(status="FAIL", message="synthetic native mismatch",
            metrics=[{"name": "synthetic-value", "value": 3}], artifacts=["synthetic-artifact.json"])
        self.save()
        process, report = self.run_suite()
        self.assertEqual(process.returncode, 1, process.stdout)
        self.assertIn("failed_test_summary: fail=2 error=0 runner_errors=0", process.stdout)
        self.assertIn("FAIL shared-model/CFGDEMO01: AssertionError: synthetic shared mismatch", process.stdout)
        self.assertIn("FAIL native-amps/FUTURECHECK9: synthetic native mismatch", process.stdout)
        run = Path(report["run_directory"])
        summary = (run/"failures.txt").read_text()
        self.assertEqual(summary, (self.output/"failures.txt").read_text())
        self.assertEqual(report["failure_summary"], str(run/"failures.txt"))
        self.assertIn(summary, process.stdout)
        for filename in ("shared.log", "native.log", "shared-build.log", "shared-discovery.log", "native-discovery.log"):
            self.assertIn(str(run/filename), summary)
            self.assertTrue((run/filename).is_file())
        for row in report["results"]:
            if row["status"] == "PASS":
                self.assertNotIn("diagnostic_log", row)
                continue
            diagnostic = Path(row["diagnostic_log"])
            self.assertTrue(diagnostic.is_file())
            self.assertIn(str(diagnostic), summary)
            if row["scope"] == "shared-model":
                self.assertIn(captured, diagnostic.read_text())
                self.assertIn(str(run/"shared.log"), diagnostic.read_text())
            else:
                self.assertIn("synthetic-value", diagnostic.read_text())
                self.assertIn("synthetic-artifact.json", diagnostic.read_text())
                self.assertIn(str(run/"native.log"), diagnostic.read_text())
        # Publishing another run must replace the convenient latest summary,
        # without changing the previous run's diagnostics or retaining its FAIL.
        self.fixture["shared"][0]["passed"] = True
        self.fixture["native"][1]["status"] = "PASS"; self.save()
        process, next_report = self.run_suite()
        self.assertEqual(process.returncode, 0, process.stdout)
        self.assertNotEqual(report["run_directory"], next_report["run_directory"])
        self.assertEqual((run/"failures.txt").read_text(), summary)
        self.assertNotIn("FAIL shared-model/CFGDEMO01", (self.output/"failures.txt").read_text())

    def test_missing_duplicate_extra_rank_and_provider_evidence(self):
        for mode in ("native-missing", "native-duplicate", "native-extra", "rank-mismatch", "obsolete-provider-report"):
            with self.subTest(mode=mode):
                self.fixture["mode"] = mode; self.save()
                process, report = self.run_suite()
                self.assertEqual(process.returncode, 2, process.stdout)
                self.assertEqual(report["counts"]["error"], 2)
                self.assertEqual(report["counts"]["pass"], 2)

    def test_stale_output_never_satisfies_a_later_crash(self):
        process, first = self.run_suite(); self.assertEqual(process.returncode, 0, process.stdout)
        self.fixture["mode"] = "missing-report"; self.save()
        process, second = self.run_suite()
        self.assertEqual(process.returncode, 2, process.stdout)
        self.assertNotEqual(first["run_directory"], second["run_directory"])
        self.assertEqual(second["counts"]["error"], 2)

    def test_process_exit_mismatch_and_shared_missing_rows(self):
        self.fixture["mode"] = "exit-mismatch"; self.fixture["native"][0]["status"] = "FAIL"; self.save()
        process, report = self.run_suite(); self.assertEqual(process.returncode, 2, process.stdout)
        self.assertEqual(report["counts"]["error"], 2)
        self.fixture["native"][0]["status"] = "PASS"; self.fixture["mode"] = "shared-missing"; self.save()
        process, report = self.run_suite(); self.assertEqual(process.returncode, 2, process.stdout)
        self.assertEqual(report["counts"]["error"], 2)
        self.assertEqual(report["counts"]["pass"], 2)

    def test_skips_remain_explicit_and_strict_policy_is_optional(self):
        self.fixture["native"][0]["status"] = "SKIP"; self.save()
        process, report = self.run_suite(); self.assertEqual(process.returncode, 0, process.stdout)
        self.assertEqual(report["counts"]["skip"], 1)
        process, report = self.run_suite("--require-no-skips")
        self.assertEqual(process.returncode, 1, process.stdout)
        self.assertEqual(report["counts"]["fail"], 0)
        self.assertEqual(ET.parse(self.output/"junit.xml").getroot().get("skipped"), "1")
        self.assertIn("strict_skip_summary: skip=1", process.stdout)
        self.assertIn("SKIP native-amps/HOSTDEMO01", process.stdout)

    def test_explicit_model_only_and_missing_executable_scope(self):
        self.amps.unlink()
        process, report = self.run_suite("--model-only")
        self.assertEqual(process.returncode, 0, process.stdout)
        self.assertEqual(report["scope"], "shared-model-only")
        self.assertEqual(report["counts"]["total"], 2)
        process, report = self.run_suite()
        self.assertEqual(process.returncode, 2, process.stdout)
        self.assertEqual(len(report["runner_errors"]), 1)
        self.assertEqual(ET.parse(self.output/"junit.xml").getroot().get("errors"), "1")
        self.assertIn("failed_test_summary: fail=0 error=0 runner_errors=1", process.stdout)
        self.assertIn("runner_error native-amps: linked --amps executable is missing/nonexecutable", process.stdout)
        self.assertNotIn("phase_log native-amps-tests=", process.stdout)

    def test_live_progress_precedes_child_exit_and_heartbeat(self):
        self.fixture["mode"] = "slow-progress"; self.save()
        command = [sys.executable, str(RUNNER), "--model-dir", str(self.model), "--model-only",
                   "--output-dir", str(self.output), "--progress-interval", "0.05", "--timeout", "6"]
        process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        lines, received = [], queue.Queue()
        def collect():
            for line in process.stdout:
                lines.append(line); received.put(line)
        reader = threading.Thread(target=collect, daemon=True); reader.start()
        try:
            found_result, found_heartbeat = False, False
            deadline = time.monotonic()+5
            while time.monotonic() < deadline and not (found_result and found_heartbeat):
                try:
                    line = received.get(timeout=0.1)
                except queue.Empty:
                    continue
                if "[progress] shared-model-tests completed=1/2" in line:
                    found_result |= "result=CFGDEMO01 PASS" in line
                    found_heartbeat |= "RUNNING" in line
            self.assertTrue(found_result, "result was buffered until child exit: " + "".join(lines))
            self.assertTrue(found_heartbeat, "quiet child produced no heartbeat: " + "".join(lines))
            self.assertIsNone(process.poll(), "progress arrived only after execution finished")
            self.assertFalse((self.output/"summary.json").exists(), "final publication occurred before release")
        finally:
            (self.root/"release-progress").write_text("release\n")
            process.wait(timeout=10); reader.join(timeout=2); process.stdout.close()
        self.assertEqual(process.returncode, 0, "".join(lines))
        self.assertIn("completed=2/2 (100.0%)", "".join(lines))
        report = json.loads((self.output/"summary.json").read_text())
        self.assertEqual(report["counts"]["pass"], 2)
        child_log = Path(report["run_directory"])/"shared.log"
        self.assertNotIn("[progress]", child_log.read_text(), "progress modified the original child log")
        self.assertEqual(child_log.read_text().count("[CFGDEMO01] PASS"), 1)

    def test_quiet_process_timeout_stays_bounded_and_preserves_output(self):
        spec = importlib.util.spec_from_file_location("aggregate_timeout_test", RUNNER)
        runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
        operations, console = [], io.StringIO()
        with contextlib.redirect_stdout(console):
            operation, text = runner.execute([sys.executable, "-c",
                "import time; print('before timeout',flush=True); time.sleep(4)"],
                self.root, self.root/"timeout.log", 0.5, operations, "quiet-test", 0.05)
        self.assertEqual(operation["returncode"], 2)
        self.assertIsNotNone(operation["error"])
        self.assertLess(operation["seconds"], 3, "timeout waited for the quiet child to finish")
        self.assertIn("before timeout", text)
        self.assertIn("RUNNING output_idle=", console.getvalue())
        self.assertIn("END exit=2 ERROR=", console.getvalue())


if __name__ == "__main__":
    unittest.main(verbosity=2)

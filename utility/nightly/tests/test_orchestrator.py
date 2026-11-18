import json
import os
import datetime as dt
import io
import subprocess
import sys
import tarfile
import tempfile
import time
import unittest
from pathlib import Path

from amps_nightly import procs, readers, transfer
from amps_nightly.machines import parse_machine
from amps_nightly.orchestrator import Night, trim_log
from amps_nightly.records import read_run_dir, write_run_dir
from amps_nightly.registry import Registry
from amps_nightly.runners import Routine
from amps_nightly.store import Store, ingest

GIT_ENV = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")

GOOD = """import json, os, sys
out = sys.argv[1]
open(os.path.join(out, "t2.log"), "w").write("details of t2\\n")
open(os.path.join(out, "extra.txt"), "w").write("x\\n")
with open(os.path.join(out, "r.jsonl"), "w") as stream:
    stream.write(json.dumps({"test": "t1", "status": "PASS"}) + "\\n")
    stream.write(json.dumps({"test": "t2", "status": "FAILED", "reason": "bad value",
                             "log": os.path.join(out, "t2.log")}) + "\\n")
"""

MAKEFILE_TEST = """test_manifest:
\t@echo "Alpha|pic-field-solver|||amps"
\t@echo "Beta|pic-field-solver|||amps"
\t@echo "Gamma|pic-field-solver||SKIP|amps"
\t@echo "Delta|pic-field-solver|||amps"
\t@echo "Epsilon|pic-field-solver|||amps"
test_compile:
\t@true
test_Alpha_compile:
\t@echo "test_Alpha_compile..." > test_Alpha.diff
\t@echo "test_Alpha_compile... done" >> test_Alpha.diff
test_Alpha_rundir:
\t@echo "test_Alpha_rundir..." >> test_Alpha.diff
\t@echo "test_Alpha_rundir... done" >> test_Alpha.diff
test_Alpha_run:
\t@echo "running Alpha with $(MPIRUN)"
test_Alpha_check:
\t@: > test_Alpha.diff
test_Epsilon_compile:
\t@echo "test_Epsilon_compile..." > test_Epsilon.diff
\t@echo "test_Epsilon_compile... done" >> test_Epsilon.diff
test_Epsilon_rundir:
\t@echo "test_Epsilon_rundir... done" >> test_Epsilon.diff
test_Epsilon_run:
\t@echo "AMPS: exit: line=8328, file=pic.h, message=Error: found a negative interpolation weight (error code=0.0)"
\t@echo "test_Epsilon_run... done" >> test_Epsilon.diff
test_Epsilon_check:
\t@echo "WARNING in DiffNum.pl: there are extra numbers in ref" >> test_Epsilon.diff
test_Delta_compile:
\t@echo "test_Delta_compile..." > test_Delta.diff
\t@echo "main.cpp:12:3: error: 'Coefficient' was not declared"
\t@echo "test_Delta_compile... done" >> test_Delta.diff
test_Delta_rundir:
\t@echo "rundir must not run after a failed build"; exit 7
test_Beta_compile:
\t@echo "test_Beta_compile..." > test_Beta.diff
\t@echo "compiler error in Beta"; exit 1
"""


def git(*args, cwd):
    subprocess.run(["git"] + list(args), cwd=str(cwd), check=True, env=GIT_ENV, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def repository(path: Path, files: dict) -> str:
    path.mkdir(parents=True)
    git("init", "-q", cwd=path)
    git("symbolic-ref", "HEAD", "refs/heads/master", cwd=path)
    for name, text in files.items():
        (path / name).write_text(text)
        if name.endswith((".pl", ".py")):
            (path / name).chmod(0o755)
    git("add", "-A", cwd=path)
    git("commit", "-q", "-m", "init", cwd=path)
    return str(path)


class Processes(unittest.TestCase):
    def test_time_limit_kills_children_in_other_sessions(self):
        old = procs.GRACE_S
        procs.GRACE_S = 1
        try:
            with tempfile.TemporaryDirectory() as tmp:
                tag = "test-{}".format(os.getpid())
                started = time.time()
                outcome = procs.run(["/bin/bash", "-c", "setsid sleep 300 & sleep 300"], Path(tmp),
                                    dict(os.environ), Path(tmp, "log"), 1.5, tag)
                log = Path(tmp, "log").read_text()
        finally:
            procs.GRACE_S = old
        self.assertTrue(outcome.timed_out)
        self.assertLess(time.time() - started, 30)
        self.assertEqual(procs._tagged(tag), [])            # the setsid grandchild is gone too
        self.assertIn("time limit of 2 s reached", log)


class Readers(unittest.TestCase):
    def read(self, name, data, filename="r.json"):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, filename)
            path.write_text(data if isinstance(data, str) else json.dumps(data))
            return readers.READERS[name](readers.ReadContext("x", str(path), Path(tmp), 1))

    def test_src_test_json(self):
        rows = self.read("src-test-json", {"output_directory": "/o", "results": [
            {"id": "CORE-A", "status": "PASS", "seconds": 3, "cases": [{"status": "PASS"}], "log": "logs/CORE-A.log"},
            {"id": "CORE-B", "status": "ERROR", "message": "runner timed out", "cases": []},
            {"id": "CORE-C", "status": "FAIL", "message": "native core-test application did not build",
             "cases": [{"status": "FAIL"}, {"status": "PASS"}]}]})
        a, b, c = rows
        self.assertEqual((a["status"], a["cases_total"], a["log_src"]), ("PASS", 1, str(Path("/o/logs/CORE-A.log"))))
        self.assertEqual((b["status"], b["stage"], b["reason"]), ("FAILED", "run", "harness: runner timed out"))
        self.assertEqual((c["stage"], c["cases_total"], c["cases_failed"]), ("build", 2, 1))

    def test_earth_report_marks_expected_failures(self):
        rows = self.read("earth-report", [
            {"command": "srcEarth/test/C5/run_C5.py --x", "line_no": 80, "expected": "F", "actual": "F",
             "exit_code": 1, "timed_out": False, "elapsed_s": 3},
            {"command": "srcEarth/test/C9/run_C9.py --a", "line_no": 90, "expected": "P", "actual": "F",
             "exit_code": None, "timed_out": True, "elapsed_s": 60},
            {"command": "srcEarth/test/C9/run_C9.py --b", "line_no": 92, "expected": "P", "actual": "P"}])
        self.assertEqual([r["test"] for r in rows], ["C5", "C9 (line 90)", "C9 (line 92)"])
        self.assertTrue(rows[0]["reason"].startswith("expected to fail"))
        self.assertEqual((rows[1]["stage"], rows[1]["reason"]), ("run", "timed out after 60 s"))

    def test_earth_report_names_the_script(self):
        commands = ["python3 srcEarth/test/SA01/run_SA01.py --data-root /d", "./srcEarth/test/C10/run_C10.py --x",
                    "env OMP_NUM_THREADS=2 python3.8 -B srcEarth/test/C3/run_C3.py", "srcEarth/test/UFluxNumerics/run_test.sh",
                    "python3 tools/check.py"]
        rows = self.read("earth-report", [{"command": c, "line_no": n, "expected": "P", "actual": "P"}
                                          for n, c in enumerate(commands)])
        self.assertEqual([r["test"] for r in rows], ["SA01", "C10", "C3", "UFluxNumerics", "check"])

    def test_earth_report_of_an_interrupted_runner(self):
        """test_runner.py stopped at the runner's time limit: every entry still gets a row."""
        rows = self.read("earth-report", [
            {"command": "srcEarth/test/C8/run_C8.py", "line_no": 223, "expected": "P", "actual": "P",
             "elapsed_s": 110},
            {"command": "python3 srcEarth/test/SA03/run_SA03.py", "line_no": 251, "expected": "P", "actual": "F",
             "exit_code": None, "timed_out": False, "elapsed_s": 2950.4, "interrupted": True, "not_started": False,
             "log_file": "/w/logs/test_064.log"},
            {"command": "python3 srcEarth/test/SA04/run_SA04.py", "line_no": 260, "expected": "P", "actual": "F",
             "elapsed_s": 0.0, "interrupted": False, "not_started": True, "log_file": ""}])
        self.assertEqual([(r["test"], r["status"]) for r in rows],
                         [("C8", "PASS"), ("SA03", "FAILED"), ("SA04", "SKIPPED")])
        self.assertEqual((rows[1]["stage"], rows[1]["reason"]), ("run", "stopped after 2950 s: the runner was stopped"))
        self.assertEqual(rows[2]["reason"], "not started: the runner was stopped before this entry")

    def test_sep_moon_junit_and_missing_report(self):
        rows = self.read("sep-orchestrator-json", {"records": [
            {"test_id": "build", "status": "PASS", "step": "build", "kind": "check"},
            {"test_id": "T1", "status": "ERROR", "step": "a", "message": "crash"},
            {"test_id": "T1", "status": "PASS", "step": "b"}]})
        self.assertEqual([r["test"] for r in rows], ["build", "a/T1", "b/T1"])
        self.assertEqual(rows[1]["reason"], "error: crash")
        rows = self.read("moon-summary", {"results": [{"id": "U01", "status": "SKIPPED", "reason": "no build"},
                                                      {"id": "U02", "status": "ERROR", "reason": "moon.input missing"}]},
                         "summary.json")
        self.assertEqual([(r["status"], r["reason"]) for r in rows],
                         [("SKIPPED", "no build"), ("FAILED", "error: moon.input missing")])
        rows = self.read("junit", '<testsuite><testcase classname="k" name="a" time="2"/>'
                                  '<testcase name="b"><failure message="boom"/></testcase>'
                                  '<testcase name="c"><skipped/></testcase></testsuite>', "r.xml")
        self.assertEqual([(r["test"], r["status"]) for r in rows], [("k.a", "PASS"), ("b", "FAILED"), ("c", "SKIPPED")])
        with self.assertRaises(readers.ResultsUnavailable):
            readers.READERS["src-test-json"](readers.ReadContext("x", "/nonexistent.json", Path("/"), 2))


@unittest.skipUnless(subprocess.run(["which", "make"], stdout=subprocess.DEVNULL).returncode == 0, "needs make")
class Night_(unittest.TestCase):
    """A whole night with a local git repository, three command runners and the Makefile.test runner."""

    def test_night(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            amps = repository(tmp / "src" / "AMPS", {
                "good.py": GOOD, "Makefile": "include Makefile.test\n", "Makefile.test": MAKEFILE_TEST,
                "Config.pl": "#!/bin/sh\necho \"configured $AMPS_TABLE_HOST $*\" > .configured\n"})
            share = repository(tmp / "src" / "share", {"README": "share\n"})
            routine = Routine({"schema_version": 1,
                               "checkout": {"repository": amps, "branch": "master",
                                            "companions": [{"path": "share", "repository": share}]},
                               "runners": [
                                   {"id": "good", "command": ["{python}", "good.py", "{out}"], "families": "core-pic",
                                    "results": {"reader": "amps-results", "file": "{out}/r.jsonl"},
                                    "logs": ["{out}/*.txt"]},
                                   {"id": "crash", "workspace": "plain", "families": "core-mesh",
                                    "command": "echo 'ERROR: boom in {table_host}'; exit 3",
                                    "results": {"reader": "amps-results", "file": "{out}/none.jsonl"}},
                                   {"id": "gpu", "requires": ["gpu"], "command": ["true"], "families": "core-pic",
                                    "results": {"reader": "exit-code"}},
                                   {"id": "legacy", "kind": "makefile"}],
                               "disabled": [{"test": "t1", "source": "good", "reason": "flaky"}]})
            machine = parse_machine({
                "schema_version": 1, "machine": "box", "workdir": str(tmp / "work"), "upload": str(tmp / "incoming"),
                "python": sys.executable, "lanes": 2, "cores": 2, "make_vars": {"MPIRUN": "mpirun -np 2"},
                "configs": [{"config": "gnu", "config_pl": ["-install -compiler=gcc"]}]}, Path("box.json"))
            (tmp / "incoming").mkdir()
            lines = []
            rundir, = Night(machine, routine, Registry.load(), machine.configs, echo=lines.append).run()

            run, rows = read_run_dir(rundir)
            got = {r["test"]: r for r in rows}
            self.assertEqual(len(run["commit"]), 40)
            self.assertEqual(got["t1"]["status"], "SKIPPED")
            self.assertEqual(got["t1"]["reason"], "disabled: flaky (its result was ignored)")
            self.assertEqual((got["t2"]["status"], got["t2"]["log"]), ("FAILED", "logs/good/tests/t2.log"))
            self.assertEqual((rundir / got["t2"]["log"]).read_text(), "details of t2\n")
            self.assertTrue((rundir / "logs/good/files/extra.txt").is_file())
            self.assertEqual(got["crash"]["status"], "FAILED")
            self.assertIn("exited with 3: ERROR: boom in box", got["crash"]["reason"])
            self.assertEqual((got["gpu"]["status"], got["gpu"]["reason"]), ("SKIPPED", "box does not declare capability gpu"))
            self.assertEqual(got["Alpha"]["status"], "PASS")
            self.assertEqual((got["Beta"]["status"], got["Beta"]["stage"]), ("FAILED", "compile"))
            self.assertEqual(got["Gamma"]["status"], "SKIPPED")
            self.assertEqual((got["Delta"]["status"], got["Delta"]["stage"]), ("FAILED", "compile"))
            self.assertEqual((got["Epsilon"]["status"], got["Epsilon"]["stage"]), ("FAILED", "run"))
            self.assertIn("the run aborted: AMPS: exit: line=8328", got["Epsilon"]["reason"])
            self.assertIn("error: 'Coefficient' was not declared", got["Delta"]["reason"])
            self.assertNotIn("rundir must not run", (rundir / "logs/legacy/tests/Delta.log").read_text())
            self.assertIn("running Alpha with mpirun -np 2", (rundir / "logs/legacy/tests/Alpha.log").read_text())
            runners = {r["runner"]: r for r in json.loads((rundir / "runners.json").read_text())}
            self.assertEqual(runners["crash"]["exit_code"], 3)
            self.assertTrue(runners["good"]["workspace_kept"])            # it has a failure
            ws = Path(runners["good"]["workspace"])
            self.assertEqual(ws, tmp / "work" / "GNU" / "good" / "AMPS")        # <workdir>/<COMPILER>/<runner>/AMPS
            self.assertTrue((ws.parent / ".amps-nightly").is_file())
            self.assertEqual((ws / ".configured").read_text().split(), ["configured", "box", "-install", "-compiler=gcc"])
            self.assertTrue((ws / "share" / "README").is_file())
            self.assertFalse((Path(runners["crash"]["workspace"]) / ".configured").exists())   # plain workspace
            self.assertFalse((rundir / transfer.PARTIAL_RUN_FILE).exists())

            finished, abandoned = transfer.runs_to_pack(tmp / "work")
            self.assertEqual((finished, abandoned), ([rundir], []))
            tarball = transfer.pack(rundir, tmp / "work" / "outbox")
            self.assertEqual(transfer.runs_to_pack(tmp / "work"), ([], []))            # PACKED
            with tarfile.open(str(tarball)) as tar:
                names = tar.getnames()
            self.assertEqual(names[-1], run["run_id"] + "/DONE")
            self.assertIn(run["run_id"] + "/logs/good/tests/t2.log", names)
            self.assertNotIn(run["run_id"] + "/PACKED", names)
            self.assertEqual(transfer.unsent(tmp / "work" / "outbox"), [tarball])
            transfer.send(tarball, str(tmp / "incoming"))
            self.assertEqual(transfer.unsent(tmp / "work" / "outbox"), [])
            self.assertTrue((tmp / "incoming" / tarball.name).is_file())

            with Store(tmp / "results.db") as st:
                accepted, rejected = ingest(tmp, tmp / "www", st)
            self.assertEqual((accepted, rejected), ([run["run_id"]], []))
            self.assertEqual(list((tmp / "incoming").iterdir()), [])

            # the next night installs into the same folders again, and a full disk stops runners from starting
            later = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=1)
            rundir2, = Night(machine, routine, Registry.load(), machine.configs, echo=lines.append, now=later).run()
            _, rows2 = read_run_dir(rundir2)
            self.assertEqual({r["test"]: r["status"] for r in rows2}["Alpha"], "PASS")
            machine.min_free_gb = 1e12
            rundir3, = Night(machine, routine, Registry.load(), machine.configs, echo=lines.append,
                             now=later + dt.timedelta(hours=1)).run()
            _, rows3 = read_run_dir(rundir3)
            # (a test that needs a capability box does not declare keeps that reason: the gpu runner, and the Table's
            # <test>--valgrind variants)
            self.assertTrue(rows3 and all(r["status"] == "SKIPPED" and "GB free" in r["reason"] for r in rows3
                                          if "does not declare capability" not in r["reason"]))
            night = list((tmp / "www" / "archive").glob("*/*/*/box-gnu"))[0]
            index = (night / "index.html").read_text()
            self.assertIn("Runners", index)
            self.assertIn('href="logs/good/tests/t2.log"', index)


class Transfer(unittest.TestCase):
    RUN = {"schema_version": 1, "run_id": "20261114T0100Z-box-gnu", "machine": "box", "config": "gnu",
           "started": "2026-11-14T01:00:00Z", "finished": "2026-11-14T02:00:00Z", "complete": True}

    def test_trim_log_keeps_each_step_head_and_tail(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp, "t.log")
            body = "".join("line {}\n".format(i) for i in range(5000))
            log.write_text("### compile\nok\n### exit 0\n### run\n" + body + "### exit 1\nbye\n")
            removed = trim_log(log, head=100, tail=200)
            text = log.read_text()
        self.assertGreater(removed, 40000)
        self.assertIn("### compile\nok\n### exit 0\n### run\nline 0\n", text)
        self.assertIn("line 4999\n### exit 1\nbye\n", text)
        self.assertIn("amps-nightly left out", text)
        self.assertLess(len(text), 600)

    def bad_tarball(self, path: Path, member: tarfile.TarInfo, data: bytes = b"x") -> Path:
        with tarfile.open(str(path), "w:gz") as tar:
            member.size = len(data) if member.isfile() else 0
            tar.addfile(member, io.BytesIO(data) if member.isfile() else None)
        return path

    def test_unpack_refuses_unsafe_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            link = tarfile.TarInfo("run/evil")
            link.type, link.linkname = tarfile.SYMTYPE, "/etc/passwd"
            for name, member in (("up", tarfile.TarInfo("run/../../x")), ("abs", tarfile.TarInfo("/tmp/x")),
                                 ("link", link), ("top", tarfile.TarInfo("x"))):
                tarball = self.bad_tarball(tmp / (name + ".tar.gz"), member)
                with self.assertRaises(transfer.TransferError, msg=name):
                    transfer.unpack(tarball, tmp / "staging")
            self.assertFalse((tmp / "x").exists())

    def test_bad_tarball_is_rejected_by_ingest(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "incoming").mkdir()
            (tmp / "incoming" / "broken.tar.gz").write_bytes(b"not a tarball")
            (tmp / "incoming" / ".arriving.tar.gz").write_bytes(b"partial")          # rsync temporary: ignored
            with Store(tmp / "results.db") as st:
                accepted, rejected = ingest(tmp, tmp / "www", st)
            self.assertEqual(accepted, [])
            self.assertEqual(len(rejected), 1)
            self.assertTrue(list((tmp / "rejected").glob("broken.tar.gz-*/REJECTED.txt")))
            self.assertTrue((tmp / "incoming" / ".arriving.tar.gz").exists())

    def test_finalize_completes_a_run_whose_orchestrator_died(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            rundir = work / "nights" / "20261113-2000" / "runs" / "box-gnu"
            (rundir / "logs" / "legacy").mkdir(parents=True)
            (rundir / "logs" / "legacy" / "runner.log").write_text("compiling ...\n")
            (rundir / transfer.PARTIAL_RUN_FILE).write_text(json.dumps(dict(self.RUN, complete=False)))
            (rundir / transfer.PROGRESS_FILE).write_text(json.dumps(
                {"family": "core-pic", "test": "t1", "status": "PASS", "source": "core", "stage": "", "reason": ""}) + "\n")
            (rundir / "runners.json").write_text(json.dumps([{"runner": "core", "state": "ran"}]))
            (rundir.parent.parent / "orchestrator.log").write_text(
                "20:00:01 --- box-gnu (x)\n20:00:02 core: start (limit 120 min, workspace w)\n"
                "20:10:00 core: 1 PASS, 0 FAILED, 0 SKIPPED in 10 min\n20:10:01 legacy: start (limit 600 min, workspace w)\n")
            finished, abandoned = transfer.runs_to_pack(work)          # no run in progress
            self.assertEqual((finished, abandoned), ([], [rundir]))
            transfer.finalize(rundir, lambda runner: "legacy-fam" if runner == "legacy" else "")
            run, rows = read_run_dir(rundir)
        self.assertFalse(run["complete"])
        self.assertIn("orchestrator stopped", run["notes"])
        got = {r["test"]: r for r in rows}
        self.assertEqual(got["t1"]["status"], "PASS")
        self.assertEqual((got["legacy"]["status"], got["legacy"]["family"], got["legacy"]["log"]),
                         ("FAILED", "legacy-fam", "logs/legacy/runner.log"))


if __name__ == "__main__":
    unittest.main()

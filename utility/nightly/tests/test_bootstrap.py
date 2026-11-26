"""bin/amps_nightly_test.sh (the ~/bin script) and `amps-nightly install`, with local git repositories."""

import contextlib
import fcntl
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

from amps_nightly import cli, install
from amps_nightly.runners import Checkout
from amps_nightly.workspace import Mirrors, WorkspaceError

SCRIPT = install.SCRIPT
GIT_ENV = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
STEP = '#!/bin/sh\necho "{} $*" >> "$(dirname "$0")/../../../../calls.txt"\n'


def git(*args, cwd):
    return subprocess.run(["git"] + list(args), cwd=str(cwd), check=True, env=GIT_ENV, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip()


def commit(repo: Path, files: dict) -> str:
    for name, text in files.items():
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        if name.endswith(".sh") or name.endswith("amps-nightly"):
            path.chmod(0o755)
    git("add", "-A", cwd=repo)
    git("commit", "-q", "-m", "c", cwd=repo)
    return git("rev-parse", "HEAD", cwd=repo)


def new_repo(path: Path, files: dict) -> str:
    path.mkdir(parents=True)
    git("init", "-q", cwd=path)
    git("symbolic-ref", "HEAD", "refs/heads/master", cwd=path)
    return commit(path, files)


@unittest.skipUnless(shutil.which("git") and shutil.which("sh"), "needs git and sh")
class Bootstrap(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        t = self.t = Path(self.tmp.name)
        self.root, self.data = t / "ampstest", t / "ampstest" / "SWMF_data"
        stub = "#!/usr/bin/env python3\nimport sys\nif sys.argv[1:3] == ['machines', '--var']:\n    print({!r})\n".format(
            str(self.data))
        self.files = {"utility/nightly/amps-nightly": stub,
                      "utility/nightly/bin/amps_nightly_test.sh": SCRIPT.read_text(),
                      "utility/nightly/nightly_run.sh": STEP.format("run"),
                      "utility/nightly/nightly_pack.sh": STEP.format("pack"),
                      "utility/nightly/nightly_update.sh": STEP.format("update")}
        self.amps = t / "src" / "AMPS"
        self.sha1 = new_repo(self.amps, self.files)
        self.swmf = t / "src" / "SWMF_data"
        new_repo(self.swmf, {"PT/ref.dat": "1\n"})
        self.bin = t / "bin"
        self.bin.mkdir()
        shutil.copy(str(SCRIPT), str(self.bin / install.SCRIPT_NAME))
        self.conf(REPOSITORY="file://" + str(self.amps))

    def tearDown(self):
        self.tmp.cleanup()

    def conf(self, **extra):
        settings = dict(AMPSTEST=str(self.root), TESTDIR=str(self.root / "Tmp_AMPS_test"), BRANCH="master", UPDATE_SITE="yes", PYTHON=sys.executable,
                        SWMF_DATA_REPOSITORY="file://" + str(self.swmf))
        settings.update(extra)
        (self.bin / install.CONF_NAME).write_text(install.conf_text(settings))

    def night(self, *args, pattern="nightly-*.log"):
        done = subprocess.run(["sh", str(self.bin / install.SCRIPT_NAME)] + list(args), stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, universal_newlines=True, timeout=120,
                              env={"HOME": str(self.t), "PATH": "/usr/bin:/bin"})
        self.output = done.stdout
        logs = sorted((self.root / "logs").glob(pattern))
        text = logs[-1].read_text() if logs else ""
        for log in logs:
            log.unlink()
        calls = [c.strip() for c in (self.root / "calls.txt").read_text().splitlines()] \
            if (self.root / "calls.txt").exists() else []
        if (self.root / "calls.txt").exists():
            (self.root / "calls.txt").unlink()
        return done.returncode, text, calls

    def test_nights(self):
        code, log, calls = self.night("--only", "moon")
        self.assertEqual(code, 0, log)
        self.assertEqual(calls, ["run --commit {} --only moon".format(self.sha1), "pack", "update"])
        self.assertEqual((self.data / "PT/ref.dat").read_text(), "1\n")
        self.assertNotIn("WARNING", log)

        # the next night: new AMPS commit with a different script, new SWMF_data commit, a test wrote into the data
        sha2 = commit(self.amps, {"utility/nightly/bin/amps_nightly_test.sh": SCRIPT.read_text() + "# changed\n"})
        commit(self.swmf, {"PT/ref.dat": "2\n"})
        (self.data / "PT/ref.dat").write_text("edited by a test\n")
        (self.data / "PT/ref.dat.tmp").write_text("left behind\n")
        code, log, calls = self.night()
        self.assertEqual(code, 0, log)
        self.assertEqual(calls[0].split()[:3], ["run", "--commit", sha2])
        self.assertIn("--note WARNING: {} differs from".format(self.bin / install.SCRIPT_NAME), calls[0])
        self.assertEqual((self.data / "PT/ref.dat").read_text(), "2\n")
        self.assertFalse((self.data / "PT/ref.dat.tmp").exists())
        self.assertTrue((self.root / "Tmp_AMPS_test/AMPS.prev").is_dir())         # TESTDIR/AMPS, as before

        # the download fails: the previous copy is used
        self.conf(REPOSITORY="file://" + str(self.t / "nowhere"))
        code, log, calls = self.night()
        self.assertEqual(calls[0].split()[:3], ["run", "--commit", sha2])
        self.assertIn("downloading AMPS failed", log)

    def test_lock_and_server_setting(self):
        self.conf(REPOSITORY="file://" + str(self.amps), UPDATE_SITE="no")
        self.root.mkdir(parents=True)
        with open(str(self.root / "nightly.lock"), "a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            code, log, calls = self.night()
        self.assertEqual((code, calls), (0, []))
        self.assertIn("another night is still running", log)
        code, log, calls = self.night()
        self.assertEqual([c.split()[0] for c in calls], ["run", "pack"])
        self.assertIn("UPDATE_SITE=no", log)

    def test_long_routine(self):
        """--routine long goes on to nightly_run.sh; the night's log is named after the routine."""
        code, log, calls = self.night("--routine", "long", pattern="long-*.log")
        self.assertEqual(code, 0, log)
        self.assertEqual(calls, ["run --commit {} --routine long".format(self.sha1), "pack", "update"])
        self.assertIn("AMPS nightly test, routine long, on", log)
        self.assertEqual(list((self.root / "logs").glob("nightly-*.log")), [])
        code, log, calls = self.night("--routine", "Bad/name", pattern="*.log")
        self.assertEqual((code, calls, log), (2, [], ""))
        self.assertIn("must be a routine name", self.output)

    def test_long_routine_waits_for_the_lock(self):
        """Two routines never run at once: the long routine waits for the lock (LOCK_WAIT_H), then gives up."""
        self.conf(REPOSITORY="file://" + str(self.amps), LOCK_WAIT_H="0.0005")       # 1.8 s
        self.root.mkdir(parents=True)
        with open(str(self.root / "nightly.lock"), "a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            code, log, calls = self.night("--routine", "long", pattern="long-*.log")
        self.assertEqual((code, calls), (0, []))
        self.assertIn("is taken (another routine is running); waiting for it up to 0.0005 h", log)
        self.assertIn("was still taken after 0.0005 h (another routine is running); routine long not run", log)
        # the nightly holding the lock finishes while the long routine waits
        self.conf(REPOSITORY="file://" + str(self.amps), LOCK_WAIT_H="0.05")         # 3 min
        lock = open(str(self.root / "nightly.lock"), "a")
        fcntl.flock(lock, fcntl.LOCK_EX)
        release = threading.Timer(1.5, lock.close)                                   # closing releases the flock
        release.start()
        try:
            code, log, calls = self.night("--routine", "long", pattern="long-*.log")
        finally:
            release.cancel()
            lock.close()
        self.assertEqual(code, 0, log)
        self.assertIn("waiting for it up to 0.05 h", log)
        self.assertEqual(calls[0], "run --commit {} --routine long".format(self.sha1))


@unittest.skipUnless(shutil.which("git"), "needs git")
class PinnedCommit(unittest.TestCase):
    def test_commit_must_be_on_the_branch(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            first = new_repo(tmp / "AMPS", {"a": "1"})
            second = commit(tmp / "AMPS", {"a": "2"})
            mirrors = Mirrors(tmp / "mirrors", Checkout(str(tmp / "AMPS"), "master", []))
            mirrors.update(lambda _: None, pin=first[:10])
            self.assertEqual(mirrors.commit, first)
            mirrors.update(lambda _: None)
            self.assertEqual(mirrors.commit, second)
            with self.assertRaises(WorkspaceError):
                mirrors.update(lambda _: None, pin="0" * 40)


class Install(unittest.TestCase):
    CRONTAB = "# old routine\n# 30 22 * * * /home/x/bin/run_test_amps.sh\nMAILTO=x\n"

    def test_install_reinstall_uninstall(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            settings = {"AMPSTEST": str(tmp / "ampstest"), "BRANCH": "master", "UPDATE_SITE": "yes",
                        "PYTHON": "/usr/bin/python3"}
            plan = install.plan_install(tmp / "bin", settings, "20:00", self.CRONTAB)
            written = []
            install.apply(plan, self.CRONTAB, writer=written.append)
            script, conf = tmp / "bin" / install.SCRIPT_NAME, tmp / "bin" / install.CONF_NAME
            self.assertEqual(script.read_bytes(), SCRIPT.read_bytes())
            self.assertTrue(os.access(str(script), os.X_OK))
            self.assertEqual(install.read_conf(conf), settings)
            line = "00 20 * * *  {} >> {} 2>&1".format(script, tmp / "ampstest" / "logs" / "cron.log")
            self.assertEqual(written[-1], self.CRONTAB + install.MARKER + "\n" + line + "\n")
            self.assertEqual(plan.backup.read_text(), self.CRONTAB)

            again = install.plan_install(tmp / "bin", settings, "20:00", written[-1])
            self.assertFalse(again.changes)
            later = install.plan_install(tmp / "bin", settings, "21:30", written[-1])
            self.assertEqual(later.crontab.count(install.SCRIPT_NAME), 1)
            self.assertIn("30 21 * * *", later.crontab)

            gone = install.plan_uninstall(tmp / "bin", written[-1], tmp / "backup")
            install.apply(gone, written[-1], writer=written.append)
            self.assertEqual(written[-1], self.CRONTAB)
            self.assertFalse(script.exists() or conf.exists())

    def test_routine_lines(self):
        """One crontab line per routine below the marker: the nightly's as before, the long routine's weekly."""
        self.assertEqual(install.cron_line("Sun 06:30", Path("/b/s.sh"), Path("/l/cron.log"), "long"),
                         "30 06 * * 0  /b/s.sh --routine long >> /l/cron.log 2>&1")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            settings = {"AMPSTEST": str(tmp / "ampstest"), "BRANCH": "master", "UPDATE_SITE": "yes",
                        "PYTHON": "/usr/bin/python3"}
            plan = install.plan_install(tmp / "bin", settings, "20:00", self.CRONTAB, routines=[("long", "Sat 10:00")])
            written = []
            install.apply(plan, self.CRONTAB, writer=written.append)
            script, log = tmp / "bin" / install.SCRIPT_NAME, tmp / "ampstest" / "logs" / "cron.log"
            nightly = "00 20 * * *  {} >> {} 2>&1".format(script, log)
            weekly = "00 10 * * 6  {} --routine long >> {} 2>&1".format(script, log)
            self.assertEqual(written[-1], self.CRONTAB + install.MARKER + "\n" + nightly + "\n" + weekly + "\n")
            self.assertFalse(install.plan_install(tmp / "bin", settings, "20:00", written[-1],
                                                  routines=[("long", "Sat 10:00")]).changes)
            only_nightly = install.plan_install(tmp / "bin", settings, "20:00", written[-1])
            self.assertEqual(only_nightly.crontab, self.CRONTAB + install.MARKER + "\n" + nightly + "\n")
            gone = install.plan_uninstall(tmp / "bin", written[-1], tmp / "backup")
            install.apply(gone, written[-1], writer=written.append)
            self.assertEqual(written[-1], self.CRONTAB)

    def test_install_command_routines(self):
        """amps-nightly install --dry-run: the lines of the machine file's routines; --routine picks what --time sets.
        The crontab is only read (mocked) and never written."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "machines").mkdir()
            (tmp / "machines" / "box.json").write_text(json.dumps({
                "schema_version": 1, "machine": "box", "upload": str(tmp / "incoming"), "start": "21:00",
                "workdir": str(tmp / "work"), "configs": [{"config": "gnu", "config_pl": ["-install"]}],
                "routines": {"long": {"start": "Sat 10:00"}, "weekly": {"start": "Sun 01:00", "enabled": False}}}))

            def install_(*extra):
                out = io.StringIO()
                with mock.patch.object(install, "read_crontab", return_value=""), \
                        mock.patch.object(install, "apply", side_effect=AssertionError("must not change anything")), \
                        contextlib.redirect_stdout(out):
                    code = cli.main(["install", "--machines", str(tmp / "machines"), "--machine", "box", "--bin",
                                     str(tmp / "bin"), "--root", str(tmp / "ampstest"), "--dry-run"] + list(extra))
                return code, [l.strip().split("  ")[0] for l in out.getvalue().splitlines() if "* *" in l]
            self.assertEqual(install_(), (0, ["00 21 * * *", "00 10 * * 6"]))
            self.assertEqual(install_("--time", "22:15"), (0, ["15 22 * * *", "00 10 * * 6"]))
            self.assertEqual(install_("--routine", "long", "--time", "Fri 23:00"), (0, ["00 21 * * *", "00 23 * * 5"]))
            self.assertEqual(install_("--time", "Sat 10:00")[0], 1)            # the nightly runs every day
            self.assertEqual(install_("--routine", "weekly", "--time", "Sun 02:00")[0], 1)    # not enabled

    def test_conf_values_with_spaces(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "c.conf")
            path.write_text(install.conf_text({"AMPSTEST": "/a dir/x", "PYTHON": "python3"}))
            self.assertEqual(install.read_conf(path), {"AMPSTEST": "/a dir/x", "PYTHON": "python3"})


if __name__ == "__main__":
    unittest.main()

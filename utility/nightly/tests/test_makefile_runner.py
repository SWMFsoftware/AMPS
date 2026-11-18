import os
import tempfile
import time
import unittest
from pathlib import Path

from amps_nightly import makefile_runner


class ConfigReset(unittest.TestCase):
    """The Config.pl files are put back to their post-install state between the tests of a night."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ws = Path(self.tmp.name)
        (self.ws / "Makefile").write_text("EXTRALINKEROPTIONS_F=-lmpi_cxx\n")
        (self.ws / "Makefile.local").write_text("SPICE=/opt/spice\n")
        (self.ws / ".amps.conf").write_text("SPICE=/opt/spice\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_options_set_by_a_test_are_undone(self):
        saved = makefile_runner.config_snapshot(self.ws)
        # what Config.pl -f-link-option=..., -compiler-option=... and -model-data-path=... of one test leave behind
        (self.ws / "Makefile").write_text("EXTRALINKEROPTIONS_F=-L/gtest/lib -lgtest -pthread\n")
        with open(str(self.ws / "Makefile.local"), "a") as stream:
            stream.write("EXTRACOMPILEROPTIONS+=-I/gtest/include\n")
        (self.ws / ".amps.conf").write_text("SPICE=/opt/spice\nMODELINPUTDATA=/data/input/const-field-bc\n")
        (self.ws / ".general.conf").write_text("created by the test\n")
        changed = makefile_runner.config_restore(self.ws, saved)
        self.assertEqual(sorted(changed), [".amps.conf", ".general.conf", "Makefile", "Makefile.local"])
        self.assertEqual((self.ws / "Makefile").read_text(), "EXTRALINKEROPTIONS_F=-lmpi_cxx\n")
        self.assertEqual((self.ws / "Makefile.local").read_text(), "SPICE=/opt/spice\n")
        self.assertEqual((self.ws / ".amps.conf").read_text(), "SPICE=/opt/spice\n")
        self.assertFalse((self.ws / ".general.conf").exists())     # did not exist after the install

    def test_unchanged_files_are_not_rewritten(self):
        saved = makefile_runner.config_snapshot(self.ws)
        os.utime(str(self.ws / "Makefile"), (1000000000, 1000000000))
        self.assertEqual(makefile_runner.config_restore(self.ws, saved), [])
        self.assertEqual((self.ws / "Makefile").stat().st_mtime, 1000000000)

    def test_a_removed_file_comes_back(self):
        saved = makefile_runner.config_snapshot(self.ws)
        (self.ws / "Makefile.local").unlink()
        self.assertEqual(makefile_runner.config_restore(self.ws, saved), ["Makefile.local"])
        self.assertEqual((self.ws / "Makefile.local").read_text(), "SPICE=/opt/spice\n")


class BuildError(unittest.TestCase):
    """Only real compiler and linker diagnostics make a test count as not built."""

    def error(self, *lines):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "X.log"
            log.write_text("\n### 2026-11-15 10:00:00 test_X_compile\n### cwd /ws\n### $ make test_X_compile\n" +
                           "\n".join(lines) + "\n### exit 0 after 1 s\n")
            return makefile_runner.build_error(log, "test_X_compile")

    def test_quoted_source_and_messages_are_not_errors(self):
        self.assertEqual(self.error(
            "Comet.cpp:130:12: warning: '%s' directive writing up to 999 bytes [-Wformat-overflow=]",
            '  130 | sprintf(msg,"Error: The neutral binary file (%s) does not exist.",fname);',
            "      |            ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "OH.cpp:1128:11: warning: passing argument 1 to 'restrict'-qualified parameter aliases with argument 3: "
            '"%s AMPS:Error: the file is not found"'), "")

    def test_compiler_errors(self):
        self.assertIn("geopackContext_", self.error(
            "make[1]: Entering directory '/ws'",
            "gridless/CutoffRigidityGridless.cpp:929:12: error: 'geopackContext_' was not declared in this scope"))
        self.assertIn("SpiceUsr.h", self.error(
            "build/models/exosphere/Exosphere.h:64:10: fatal error: SpiceUsr.h: No such file or directory"))
        self.assertIn("IMPLICIT", self.error("interface.f90:10:5:", "", "Error: Symbol 'x' at (1) has no IMPLICIT type"))

    def test_linker_errors(self):
        self.assertIn("undefined reference", self.error("main.o: in function `main':",
                                                        "main.cpp:(.text+0x1d): undefined reference to `foo()'",
                                                        "collect2: error: ld returned 1 exit status"))


if __name__ == "__main__":
    unittest.main()


class WorkspaceTable(unittest.TestCase):
    """Results are attributed with the tested workspace's Table, not this tool's."""

    def test_test_only_in_workspace_table_gets_its_result(self):
        from amps_nightly import legacy
        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp)
            (ws / "MakefileTest").mkdir()
            (ws / "MakefileTest" / "Table").write_text(
                "#>\nName=unit_test\nFamily=unit-tests\n<#\n#>\nName=unit_test_new_only_here\nFamily=unit-tests\n<#\n")
            (ws / "test_unit_test.diff").write_text("")
            (ws / "test_unit_test_new_only_here.diff").write_text("")
            results = {r["test"]: r for r in legacy.collect_results(legacy.read_diff_files(ws),
                                                                    makefile_runner.workspace_table(ws))}
            self.assertEqual(results["unit_test_new_only_here"]["status"], "PASS")
            self.assertEqual(results["unit_test"]["status"], "PASS")

    def test_falls_back_to_the_tool_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertTrue(len(makefile_runner.workspace_table(Path(tmp))) > 0)


DEPENDENT_TABLE = """#>
Name=Moon
Family=planets
<#
#>
Name=Moon-restart
Family=planets
Run=
>>>cp ${TESTMoonDIR}/ParticleData.restart* ${TESTMoon-restartDIR}
<#
"""

# Moon's run writes its restart file after a while; Moon-restart's run copies it, and its check reports whether it got it
DEPENDENT_MAKEFILE = """\
test_%_compile:
\t@echo "test_$*_compile..." > test_$*.diff; touch amps
test_%_rundir:
\t@mkdir -p "run_test_$*"; mv amps "run_test_$*"; echo "test_$*_rundir... done" >> test_$*.diff
test_Moon_run:
\t@sleep 1.5; touch run_test_Moon/ParticleData.restart
test_Moon-restart_run:
\tcp run_test_Moon/ParticleData.restart* run_test_Moon-restart
test_Moon_check:
\t@: > test_Moon.diff
test_Moon-restart_check:
\t@if [ -f run_test_Moon-restart/ParticleData.restart ]; then : > test_Moon-restart.diff; \\
\telse echo "no restart file" > test_Moon-restart.diff; fi
"""


class RunDependencies(unittest.TestCase):
    """A test that uses another test's run directory (${TEST<x>DIR}) runs after that test's run (Moon-restart)."""

    def test_table_records_the_run_directories_a_block_uses(self):
        from amps_nightly.table import read_table
        uses = {t.name: t.uses for t in read_table()}
        self.assertEqual(uses["Moon-restart"], {"Moon"})
        self.assertEqual(uses["Moon"], set())
        self.assertEqual(uses["Periodic-BC-NonUniform"], set())     # its check counts its own ranks

    def test_order_puts_the_tests_waited_for_first(self):
        from amps_nightly.table import TableTest

        def block(name, *uses):
            return TableTest(name=name, family="", selector="", line=0, uses=set(uses))

        table = [block("B", "A"), block("A"), block("C", "D"), block("D", "C"), block("E", "missing")]
        order, waits = makefile_runner.run_order(["B", "A", "C", "D", "E"], table)
        self.assertEqual(order, ["A", "B", "D", "C", "E"])       # A moves before B; otherwise the Table order
        self.assertEqual(waits["B"], ["A"])
        self.assertEqual(waits["A"], [])
        self.assertEqual(waits["E"], [])                         # a test that is not run is not waited for
        self.assertEqual(sorted(len(waits[n]) for n in ("C", "D")), [0, 1])   # a cycle: one of them does not wait

    def run_night(self):
        from amps_nightly import plan as plan_mod
        from amps_nightly.machines import parse_machine
        with tempfile.TemporaryDirectory() as tmp:
            ws, logs = Path(tmp, "ws"), Path(tmp, "logs")
            (ws / "MakefileTest").mkdir(parents=True)
            (ws / "MakefileTest" / "Table").write_text(DEPENDENT_TABLE)
            (ws / "Makefile").write_text(DEPENDENT_MAKEFILE)
            m = parse_machine({"schema_version": 1, "machine": "box", "upload": "/srv/incoming", "lanes": 2,
                               "configs": [{"config": "gnu", "config_pl": ["-install"]}]}, Path("box.json"))
            items = [plan_mod.makefile_item(n, "planets", "") for n in ("Moon", "Moon-restart")]
            env = {"PATH": "/usr/local/bin:/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp")}
            messages = []
            rows = makefile_runner.run(ws, m, items, env, logs, "t", time.time() + 300, messages.append)
            got = {r["test"]: r for r in rows}
            return got, Path(got["Moon-restart"]["log_src"]).read_text(), messages

    def test_dependent_run_waits(self):
        got, log, messages = self.run_night()
        self.assertEqual(got["Moon"]["status"], "PASS")
        self.assertEqual(got["Moon-restart"]["status"], "PASS")
        self.assertIn("waiting for test_Moon_run: this test uses run_test_Moon", log)
        self.assertIn("Makefile.test: Moon-restart runs after Moon (it uses their run directories)", messages)

    def test_without_the_wait_the_copy_finds_nothing(self):
        """The control: with two lanes and no waiting, Moon-restart starts beside Moon and fails (2026-11-18)."""
        from unittest import mock
        with mock.patch.object(makefile_runner, "run_order", lambda ready, table: (list(ready), {n: [] for n in ready})):
            got, log, _ = self.run_night()
        self.assertEqual(got["Moon-restart"]["status"], "FAILED")
        self.assertIn("cannot stat", log)

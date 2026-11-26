import os
import shutil
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


@unittest.skipUnless(shutil.which("g++") and shutil.which("make"), "needs g++ and make")
class LongManifest(unittest.TestCase):
    """A tested tree without a manifest target gets its list from this tool's table_parser on the tree's own Table:
    test_manifest the nightly tests, test_long_manifest the Tier=long ones."""

    def test_fallback_per_target(self):
        from amps_nightly.machines import parse_machine
        m = parse_machine({"schema_version": 1, "machine": "box", "upload": "/srv/incoming",
                           "configs": [{"config": "gnu", "config_pl": ["-install"]}]}, Path("box.json"))
        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp)
            (ws / "MakefileTest").mkdir()
            (ws / "MakefileTest" / "Table").write_text("#>\nName=Alpha\nFamily=pic-field-solver\n<#\n"
                                                       "#>\nName=Zeta-long\nFamily=sep-long\nTier=long\nRanks=4\n<#\n")
            (ws / "Makefile.conf").write_text("COMPILE.c=gcc\n")
            messages = []
            long_ = makefile_runner.manifest(ws, m, dict(os.environ), messages.append, "test_long_manifest")
            nightly = makefile_runner.manifest(ws, m, dict(os.environ), messages.append)
        self.assertEqual((long_, nightly), ([("Zeta-long", "sep-long", "")], [("Alpha", "pic-field-solver", "")]))
        self.assertIn("make test_long_manifest unavailable in the tested commit", messages[0])
        self.assertIn("make test_manifest unavailable in the tested commit", messages[1])


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


LANES_TABLE = """#>
Name=Coupled
Family=planets
Type=swmf
<#
#>
Name=A
Family=planets
<#
#>
Name=B
Family=planets
<#
#>
Name=C
Family=planets
<#
"""

# each compile records the tree it ran in (PTDIR of Makefile.def, as $(PTDIR)/data/input/... in the Table) and its -j
LANES_MAKEFILE = """\
include Makefile.def
test_%_compile:
\t@echo "test_$*_compile..." > test_$*.diff; sleep 1; echo "PTDIR=$(PTDIR) cwd=$$(pwd) jobs=$(TESTCOMPILEJOBS)" > amps
test_%_rundir:
\t@mkdir -p run_test_$*; mv amps run_test_$*; echo "test_$*_rundir... done" >> test_$*.diff
test_%_run:
\t@cp run_test_$*/amps run_test_$*/out
test_%_check:
\t@if [ -f run_test_$*/out ]; then : > test_$*.diff; else echo "no output" > test_$*.diff; fi
"""


class CompileLanes(unittest.TestCase):
    """compile_lanes > 1: tests compile at once in copies of the workspace; runs and checks stay in the workspace."""

    def night(self, lanes, cores=6, table=LANES_TABLE, makefile=LANES_MAKEFILE, names=("Coupled", "A", "B", "C")):
        from amps_nightly import plan as plan_mod
        from amps_nightly.machines import parse_machine
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        ws, logs = Path(tmp.name, "AMPS"), Path(tmp.name, "logs")
        (ws / "MakefileTest").mkdir(parents=True)
        (ws / "MakefileTest" / "Table").write_text(table)
        (ws / "Makefile").write_text(makefile)
        (ws / "Makefile.def").write_text("PTDIR = {0}\nMYDIR = {0}\n".format(ws))
        m = parse_machine({"schema_version": 1, "machine": "box", "upload": "/srv/incoming", "lanes": 2,
                           "cores": cores, "compile_lanes": lanes, "min_free_gb": 0,
                           "configs": [{"config": "gnu", "config_pl": ["-install"]}]}, Path("box.json"))
        items = [plan_mod.makefile_item(n, "planets", "") for n in names]
        env = {"PATH": "/usr/local/bin:/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp")}
        messages = []
        rows = makefile_runner.run(ws, m, items, env, logs, "t", time.time() + 300, messages.append, m.configs[0])
        built = {name: (ws / "run_test_{}".format(name) / "amps").read_text().split()
                 for name in names if (ws / "run_test_{}".format(name) / "amps").exists()}
        return ws, {r["test"]: r for r in rows}, built, messages

    def test_tests_compile_in_copies_and_run_in_the_workspace(self):
        ws, got, built, messages = self.night(3)
        self.assertEqual({n: r["status"] for n, r in got.items()}, dict.fromkeys(("Coupled", "A", "B", "C"), "PASS"))
        lanes = {str(ws.parent / "lane2" / "AMPS"), str(ws.parent / "lane3" / "AMPS")}
        # the SWMF test compiles in the workspace; A and B, the first tests the copies take, in the two copies
        self.assertEqual(built["Coupled"], ["PTDIR={}".format(ws), "cwd={}".format(ws), "jobs=2"])
        for name in ("A", "B"):
            ptdir, cwd, jobs = built[name]
            self.assertIn(cwd[len("cwd="):], lanes)
            self.assertEqual(ptdir, "PTDIR=" + cwd[len("cwd="):])      # Makefile.def of the copy names the copy
            self.assertEqual(jobs, "jobs=2")                          # 6 cores / 3 compile lanes
        self.assertNotEqual(built["A"][1], built["B"][1])
        for lane in lanes:                                            # what the copies built was handed over
            self.assertEqual(sorted(p.name for p in Path(lane).glob("run_test_*")), [])
            self.assertEqual(sorted(p.name for p in Path(lane).glob("*.diff")), [])
        self.assertIn("Makefile.test: 4 tests to run; compile: 3 lanes, make -j2 each; run: 2 lanes", messages)
        self.assertTrue(any(m.startswith("Makefile.test: compile lane 2: ") and "Makefile.def" in m for m in messages))

    def test_one_lane_compiles_in_the_workspace(self):
        ws, got, built, messages = self.night(1)
        self.assertEqual({r["status"] for r in got.values()}, {"PASS"})
        self.assertEqual({b[1] for b in built.values()}, {"cwd={}".format(ws)})
        self.assertEqual({b[2] for b in built.values()}, {"jobs=6"})
        self.assertFalse((ws.parent / "lane2").exists())

    def test_a_copy_that_cannot_be_made_leaves_fewer_lanes(self):
        from unittest import mock

        def fail(ws, lane):
            raise OSError("No space left on device")
        with mock.patch.object(makefile_runner, "compile_lane", fail):
            ws, got, built, messages = self.night(3)
        self.assertEqual({r["status"] for r in got.values()}, {"PASS"})
        self.assertEqual({b[1] for b in built.values()}, {"cwd={}".format(ws)})
        self.assertIn("Makefile.test: compile lane 2 not made (No space left on device): the tests compile in 1 lane",
                      messages)

    def test_copy(self):
        """The copy keeps links, leaves out what the tests make, and names its own path where the workspace's was."""
        with tempfile.TemporaryDirectory() as tmp:
            ws, lane = Path(tmp, "AMPS"), Path(tmp, "lane2", "AMPS")
            for d in ("data/output", "share/lib", "run_test_X", "SWMF/PT", "srcTemp", ".git", "src"):
                (ws / d).mkdir(parents=True)
            (ws / "Makefile.def").write_text("PTDIR = {0}\nMYDIR = {0}\n".format(ws))
            (ws / "Makefile.local").write_text("SPICE=/opt/spice\n")
            (ws / "share" / "lib" / "libSHARE.a").write_bytes(b"\0" + str(ws).encode())     # binary: left alone
            (ws / ".git" / "config").write_text(str(ws))
            (ws / "src" / "pic.cpp").write_text("int main(){}\n")
            (ws / "output").symlink_to("data/output")
            (ws / "SWMF_data").symlink_to(tmp)
            for name in ("amps", "test_X.diff", "run_test_X/amps", "SWMF/PT/x"):
                (ws / name).write_text("made by a test")
            rewritten = makefile_runner.compile_lane(ws, lane)
            self.assertEqual(rewritten, ["Makefile.def"])
            self.assertEqual((lane / "Makefile.def").read_text(), "PTDIR = {0}\nMYDIR = {0}\n".format(lane))
            self.assertEqual((lane / "Makefile.local").read_text(), "SPICE=/opt/spice\n")
            self.assertEqual((lane / "share" / "lib" / "libSHARE.a").read_bytes(), b"\0" + str(ws).encode())
            self.assertEqual((lane / ".git" / "config").read_text(), str(ws))
            self.assertEqual(os.readlink(str(lane / "output")), "data/output")
            self.assertEqual(os.readlink(str(lane / "SWMF_data")), tmp)
            self.assertEqual(sorted(p.name for p in lane.iterdir()),
                             [".git", "Makefile.def", "Makefile.local", "SWMF_data", "data", "output", "share", "src"])

    def test_a_test_using_another_compiles_after_it_in_its_tree(self):
        """Moon-restart builds only where the compile of Moon left moon.input.spicekernels (2026-11-18)."""
        table = LANES_TABLE + ("#>\nName=D\nFamily=planets\nRun=\n>>>cp ${TESTCDIR}/amps ${TESTDDIR}/from-C\n<#\n"
                               "#>\nName=E\nFamily=planets\n<#\n")
        makefile = LANES_MAKEFILE + (
            "test_C_compile:\n"
            "\t@echo \"test_C_compile...\" > test_C.diff; sleep 1; touch c.staged; echo \"cwd=$$(pwd)\" > amps\n"
            "test_D_compile:\n"
            "\t@echo \"test_D_compile...\" > test_D.diff; test -f c.staged; echo \"cwd=$$(pwd)\" > amps\n")
        names = ("Coupled", "A", "B", "C", "D", "E")
        ws, got, built, messages = self.night(4, table=table, makefile=makefile, names=names)
        self.assertEqual({n: r["status"] for n, r in got.items()}, dict.fromkeys(names, "PASS"))
        self.assertEqual(built["D"], built["C"])                     # the same tree, after C


class CompileGroups(unittest.TestCase):
    def test_units(self):
        from amps_nightly.table import TableTest

        def block(name, *uses):
            return TableTest(name=name, family="", selector="", line=0, uses=set(uses))
        table = [block("Moon"), block("Moon-restart", "Moon"), block("A"), block("B", "Z"), block("Y", "X"), block("X")]
        self.assertEqual(makefile_runner.compile_groups(["Moon-restart", "A", "Moon", "B", "Y", "X"], table),
                         [["Moon", "Moon-restart"], ["A"], ["B"], ["X", "Y"]])

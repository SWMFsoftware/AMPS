import os
import tempfile
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

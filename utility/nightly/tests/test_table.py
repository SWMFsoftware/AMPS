import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from amps_nightly.registry import Registry
from amps_nightly.table import DEFAULT_TABLE, check_table, read_table

REPO = Path(__file__).resolve().parents[3]
SNIPPET = """header line Name=<APP> is documentation, not a block
#>
Name=test/Alpha@@(0+passion)
Family=fam-a
Time=1
Check=
>>>@echo "check" > test_Alpha$(TESTAlpha-REF).diff
>>>-@DiffNum.pl a b > test_Alpha-extra.diff
>>>Family=not-a-keyword-inside-a-section
<#
#>
Name=Beta
Time=1
<#
"""


JOBS_SNIPPET = """#>
Name=Plain
Family=fam-a
Compile=
>>>@($(MAKE) -j8 amps)
>>>make -j$(TESTCOMPILEJOBS)
>>>-@($(MAKE) -j1 amps)
Rundir=
>>>cp -r SWMF_data/SP/x $(TESTPlainDIR)/; cp SWMF.exe $(TESTPlainDIR)
<#
#>
Name=Coupled
Family=fam-a
Type=swmf
<#
#>
Name=UsesTree
Family=fam-a
Compile=
>>>cd SWMF-IH/SWMF; ${MAKE} -j$(TESTCOMPILEJOBS) LIB
<#
"""


class TableReader(unittest.TestCase):
    def test_snippet(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "Table")
            path.write_text(SNIPPET)
            alpha, beta = read_table(path)
        self.assertEqual((alpha.name, alpha.family, alpha.selector, alpha.family_lines), ("Alpha", "fam-a", "(0+passion)", 1))
        self.assertEqual(alpha.result_files, {"Alpha-extra"})
        self.assertEqual((beta.name, beta.family, beta.family_lines), ("Beta", "", 0))

    def test_type_swmf_trees_and_fixed_jobs(self):
        """Tests that build in the workspace's SWMF trees cannot compile in a copy of it; make -j<N> is reported."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "Table")
            path.write_text(JOBS_SNIPPET)
            plain, coupled, tree = read_table(path)
        self.assertEqual([(t.type, t.swmf_tree, t.workspace_only) for t in (plain, coupled, tree)],
                         [("amps", False, False), ("swmf", False, True), ("amps", True, True)])
        self.assertEqual(plain.fixed_jobs, [(5, "@($(MAKE) -j8 amps)")])     # -j$(TESTCOMPILEJOBS) and -j1 are fine
        self.assertEqual(tree.fixed_jobs, [])
        jobs = [p for p in check_table([plain, coupled, tree], Registry.load()) if "make jobs" in p]
        self.assertEqual(jobs, ["Table line 5 (Plain): '@($(MAKE) -j8 amps)' fixes the number of make jobs: use "
                                "-j$(TESTCOMPILEJOBS) (the nightly sets it from compile_jobs of the machine file; "
                                "8 otherwise)"])

    def test_repository_table_has_registered_families(self):
        problems = check_table(read_table(DEFAULT_TABLE), Registry.load())
        self.assertEqual(problems, [])


@unittest.skipUnless(shutil.which("g++") and shutil.which("make"), "needs g++ and make")
class TableParserManifest(unittest.TestCase):
    """table_parser.cpp emits TEST<name>FAMILY and a test_manifest target listing every selected test."""

    def test_manifest(self):
        registry = Registry.load()
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run(["g++", "-std=c++14", str(REPO / "utility/TestScripts/table_parser.cpp"), "-o", "table_parser"],
                           cwd=tmp, check=True)
            shutil.copy(str(DEFAULT_TABLE), str(Path(tmp, "Table")))
            subprocess.run(["./table_parser", "4"], cwd=tmp, check=True, stdout=subprocess.DEVNULL)
            out = subprocess.run(["make", "-s", "-f", "Makefile.test", "test_manifest"], cwd=tmp, check=True,
                                 stdout=subprocess.PIPE, universal_newlines=True).stdout
        lines = [l.split("|") for l in out.splitlines() if l.strip()]
        self.assertGreater(len(lines), 50)
        for fields in lines:
            self.assertEqual(len(fields), 5, fields)
            name, family, _ref, _code, kind = fields
            self.assertIsNotNone(registry.family(family), "{} has unregistered family {!r}".format(name, family))
            self.assertIn(kind, ("amps", "swmf"))

    def test_compile_jobs_variable(self):
        """The compile rules run make -j$(TESTCOMPILEJOBS): 8 unless the command line or the environment sets it."""
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run(["g++", "-std=c++14", str(REPO / "utility/TestScripts/table_parser.cpp"), "-o", "table_parser"],
                           cwd=tmp, check=True)
            Path(tmp, "Table").write_text(SNIPPET)
            subprocess.run(["./table_parser", "2"], cwd=tmp, check=True, stdout=subprocess.DEVNULL,
                           env=dict(os.environ, AMPS_TABLE_HOST="passion"))

            def compile_lines(*args, **env):
                # make -n prints the recipe; it still runs the recursive $(MAKE) amps, which fails here (no amps rule)
                out = subprocess.run(["make", "-n", "-f", "Makefile.test", "test_Beta_compile"] + list(args), cwd=tmp,
                                     stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True,
                                     env=dict(os.environ, **env)).stdout
                return [line.strip() for line in out.splitlines() if "amps)" in line]
            self.assertEqual(compile_lines(), ["(make -j8 amps)"])
            self.assertEqual(compile_lines("TESTCOMPILEJOBS=9"), ["(make -j9 amps)"])
            self.assertEqual(compile_lines(TESTCOMPILEJOBS="12"), ["(make -j12 amps)"])

    def test_host_override(self):
        """AMPS_TABLE_HOST (the machine file's table_host) replaces the host name in the @@(...) selection."""
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run(["g++", "-std=c++14", str(REPO / "utility/TestScripts/table_parser.cpp"), "-o", "table_parser"],
                           cwd=tmp, check=True)
            Path(tmp, "Table").write_text(SNIPPET)
            names = {}
            for host in ("passion", "newbox"):
                subprocess.run(["./table_parser", "2"], cwd=tmp, check=True, stdout=subprocess.DEVNULL,
                               env=dict(os.environ, AMPS_TABLE_HOST=host))
                text = Path(tmp, "Makefile.test").read_text()
                self.assertIn("HOSTNAME={}\n".format(host), text)
                out = subprocess.run(["make", "-s", "-f", "Makefile.test", "test_manifest"], cwd=tmp, check=True,
                                     stdout=subprocess.PIPE, universal_newlines=True).stdout
                names[host] = [line.split("|")[0] for line in out.splitlines()]
        self.assertEqual(names, {"passion": ["Alpha", "Beta"], "newbox": ["Beta"]})


if __name__ == "__main__":
    unittest.main()

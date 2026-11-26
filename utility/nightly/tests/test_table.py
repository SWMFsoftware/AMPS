import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from amps_nightly.registry import Registry
from amps_nightly.table import DEFAULT_TABLE, check_table, read_table, tier_tests

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


TIER_SNIPPET = """#>
Name=Event
Family=sep
<#
#>
Name=Long-event--Parker
Family=sep-long
Tier=long
Ranks=32
Compile=
>>>Tier=nightly is not a keyword inside a section
<#
#>
Name=Weekly
Family=sep-long
Tier=weekly
<#
#>
Name=NoRanks
Family=sep-long
Tier=long
<#
#>
Name=ZeroRanks
Family=sep-long
Tier=long
Ranks=0
<#
#>
Name=NightlyRanks
Family=sep
Ranks=8
<#
#>
Name=LongInNightlyFamily
Family=sep
Tier=long
Ranks=4
<#
#>
Name=NightlyInLongFamily
Family=sep-long
<#
#>
Name=Event--long
Family=sep-long
Tier=long
Ranks=4
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

    def test_tier_and_ranks(self):
        """Tier= and Ranks= are keywords of the block's keyword part; a block without Tier= is nightly."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "Table")
            path.write_text(TIER_SNIPPET)
            tests = {t.name: t for t in read_table(path)}
        event, long_ = tests["Event"], tests["Long-event--Parker"]
        self.assertEqual((event.tier, event.ranks, event.long), ("nightly", None, False))
        self.assertEqual((long_.tier, long_.ranks, long_.long, long_.ranks_text), ("long", 32, True, "32"))
        self.assertEqual((tests["Weekly"].tier, tests["ZeroRanks"].ranks, tests["NoRanks"].ranks), ("weekly", None, None))
        self.assertEqual([t.name for t in tier_tests(list(tests.values()), "nightly")],
                         ["Event", "NightlyRanks", "NightlyInLongFamily"])

    def test_tier_rules(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "Table")
            path.write_text(TIER_SNIPPET)
            tests = read_table(path)
            problems = check_table(tests, Registry.load())
        line = {t.name: t.line for t in tests}
        self.assertEqual(problems, [p.format(**{k.replace("-", "_"): v for k, v in line.items()}) for p in [
            "Table line {Weekly} (Weekly): Tier=weekly is not nightly or long",
            "Table line {NoRanks} (NoRanks): Tier=long needs Ranks=<positive integer> (the MPI ranks of its run)",
            "Table line {ZeroRanks} (ZeroRanks): Tier=long needs Ranks=<positive integer> (the MPI ranks of its run; "
            "got Ranks=0)",
            "Table line {NightlyRanks} (NightlyRanks): Ranks= is only for Tier=long blocks (a nightly test's ranks are "
            "set by its run rule)",
            "Table line {LongInNightlyFamily} (LongInNightlyFamily): Family=sep belongs to group legacy, not long (the "
            "group of Tier=long tests)",
            "Table line {NightlyInLongFamily} (NightlyInLongFamily): Family=sep-long belongs to group long, not legacy "
            "(the group of nightly tests)",
            "Table line {Event__long} (Event--long): the name of a Tier=long test must not start with the name of the "
            "nightly test Event (its result files test_Event--long*.diff would be taken as Event's); rename it"]])

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

    def test_long_tier_goes_to_its_own_manifest(self):
        """The contract with table_parser: a Tier=long block is in Makefile.test.long and its test_long_manifest
        (same five fields), with TEST<name>RANKS; test_manifest never lists it."""
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run(["g++", "-std=c++14", str(REPO / "utility/TestScripts/table_parser.cpp"), "-o", "table_parser"],
                           cwd=tmp, check=True)
            Path(tmp, "Table").write_text(SNIPPET + "#>\nName=Zeta-long\nFamily=sep-long\nTier=long\nRanks=4\n<#\n")
            subprocess.run(["./table_parser", "2"], cwd=tmp, check=True, stdout=subprocess.DEVNULL,
                           env=dict(os.environ, AMPS_TABLE_HOST="passion"))

            def manifest(*files):
                args = [a for f in files for a in ("-f", f)]
                out = subprocess.run(["make", "-s"] + args + [files[-1].replace("Makefile.test", "test").replace(
                    ".long", "_long") + "_manifest"], cwd=tmp, check=True, stdout=subprocess.PIPE,
                    universal_newlines=True).stdout
                return [line.split("|") for line in out.splitlines()]
            self.assertEqual([f[0] for f in manifest("Makefile.test")], ["Alpha", "Beta"])
            self.assertEqual(manifest("Makefile.test", "Makefile.test.long"), [["Zeta-long", "sep-long", "", "", "amps"]])
            self.assertIn("TESTZeta-longRANKS=4", Path(tmp, "Makefile.test.long").read_text())

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

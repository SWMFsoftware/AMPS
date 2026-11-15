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


class TableReader(unittest.TestCase):
    def test_snippet(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "Table")
            path.write_text(SNIPPET)
            alpha, beta = read_table(path)
        self.assertEqual((alpha.name, alpha.family, alpha.selector, alpha.family_lines), ("Alpha", "fam-a", "(0+passion)", 1))
        self.assertEqual(alpha.result_files, {"Alpha-extra"})
        self.assertEqual((beta.name, beta.family, beta.family_lines), ("Beta", "", 0))

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

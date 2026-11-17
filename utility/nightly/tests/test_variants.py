"""Makefile.test variants: a Table test <test>--valgrind is shown as <test>(valgrind) and needs the capability valgrind."""

import datetime as dt
import json
import os
import tempfile
import time
import unittest
from pathlib import Path

from amps_nightly import legacy, makefile_runner
from amps_nightly import plan as plan_mod
from amps_nightly.machines import parse_machine
from amps_nightly.registry import Registry
from amps_nightly.runners import Routine
from amps_nightly.table import display_name, read_table, required_capabilities, table_name

TODAY = dt.date(2026, 11, 17)
TABLE = """#>
Name=test/Mars-ions@@(0+box)
Family=planets
<#
#>
Name=test/Mars-ions--valgrind@@(0+box)
Family=planets
<#
#>
Name=unit_test_field-lines--valgrind
Family=unit-tests
<#
"""


def machine(**extra):
    raw = {"schema_version": 1, "machine": "box", "upload": "/srv/incoming",
           "configs": [{"config": "gnu", "config_pl": ["-install -compiler=gfortran,gcc_mpicc"]}]}
    raw.update(extra)
    return parse_machine(raw, Path("box.json"))


class Names(unittest.TestCase):
    def test_display_and_table_names(self):
        self.assertEqual(display_name("Mars-ions--valgrind"), "Mars-ions(valgrind)")
        self.assertEqual(display_name("unit_test_field-lines--valgrind"), "unit_test_field-lines(valgrind)")
        for name in ("Mars-ions", "fast-wave--legacy-ecsim", "SEP--Parker_spiral--FTE", "--valgrind", "x--valgrind--y"):
            self.assertEqual(display_name(name), name)            # only a known variant suffix is mapped
        self.assertEqual(table_name("Mars-ions(valgrind)"), "Mars-ions--valgrind")
        self.assertEqual(table_name("Mars-ions"), "Mars-ions")

    def test_required_capabilities(self):
        self.assertEqual(required_capabilities("Mars-ions--valgrind"), ["valgrind"])
        self.assertEqual(required_capabilities("Mars-ions(valgrind)"), ["valgrind"])
        self.assertEqual(required_capabilities("Mars-ions"), [])

    def test_log_file_is_named_as_shown(self):
        self.assertEqual(makefile_runner.test_log(Path("logs/tests"), "Mars-ions--valgrind"),
                         Path("logs/tests/Mars-ions(valgrind).log"))
        self.assertEqual(makefile_runner.test_log(Path("logs/tests"), "Mars-ions"), Path("logs/tests/Mars-ions.log"))


class Gating(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.table = Path(self.tmp.name, "Table")
        self.table.write_text(TABLE)
        self.manifest = Path(self.tmp.name, "manifest.json")
        self.manifest.write_text(json.dumps({"campaigns": ["nightly"], "tests": []}))
        self.routine = Routine({"schema_version": 1, "runners": [{"id": "legacy", "kind": "makefile"}]})

    def tearDown(self):
        self.tmp.cleanup()

    def legacy_items(self, m):
        plans = plan_mod.plan(m, m.configs[0], self.routine, Registry.load(), TODAY, table=self.table,
                              manifest=self.manifest)
        return {i.test: i for sp in plans if sp.source.id == "legacy" for i in sp.items}

    def test_skipped_without_the_capability(self):
        items = self.legacy_items(machine())
        self.assertEqual(sorted(items), ["Mars-ions", "Mars-ions(valgrind)", "unit_test_field-lines(valgrind)"])
        self.assertEqual((items["Mars-ions"].action, items["Mars-ions"].make_name), ("run", "Mars-ions"))
        variant = items["Mars-ions(valgrind)"]
        self.assertEqual((variant.action, variant.make_name, variant.family), ("skip", "Mars-ions--valgrind", "planets"))
        self.assertEqual(variant.reason, "box does not declare capability valgrind")
        self.assertEqual(items["unit_test_field-lines(valgrind)"].family, "unit-tests")

    def test_run_with_the_capability(self):
        items = self.legacy_items(machine(capabilities=["valgrind"]))
        self.assertEqual({name: i.action for name, i in items.items()},
                         {"Mars-ions": "run", "Mars-ions(valgrind)": "run", "unit_test_field-lines(valgrind)": "run"})

    def test_repository_files(self):
        """The Table has the two variants in the families of their tests, and passion declares valgrind."""
        tests = {t.name: t.family for t in read_table()}
        for name in ("Mars-ions", "unit_test_field-lines"):
            self.assertEqual(tests.get(name + "--valgrind"), tests.get(name))
        passion = json.loads((Path(__file__).resolve().parents[1] / "machines" / "passion.json").read_text())
        self.assertIn("valgrind", passion["capabilities"])


class Results(unittest.TestCase):
    def test_variant_result_file_is_its_own_row(self):
        table = read_table_text(TABLE)
        files = {"test_Mars-ions.diff": (0, ""),
                 "test_Mars-ions--valgrind.diff": (60, "memcheck: 1 distinct error(s); amps exit status 0\n")}
        rows = {r["test"]: r for r in legacy.collect_results(files, table)}
        self.assertEqual(sorted(rows), ["Mars-ions", "Mars-ions(valgrind)"])     # not folded into Mars-ions
        self.assertEqual(rows["Mars-ions"]["status"], "PASS")
        self.assertEqual(rows["Mars-ions(valgrind)"]["status"], "FAILED")
        self.assertIn("memcheck: 1 distinct error(s)", rows["Mars-ions(valgrind)"]["reason"])
        self.assertEqual(rows["Mars-ions(valgrind)"]["family"], "planets")


def read_table_text(text):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp, "Table")
        path.write_text(text)
        return read_table(path)


MAKEFILE = """\
D=test_$(T).diff
test_%_compile:
\t@echo "test_$*_compile..." > test_$*.diff; touch amps; echo "test_$*_compile... done" >> test_$*.diff
test_%_rundir:
\t@echo "test_$*_rundir..." >> test_$*.diff; mkdir -p "run_test_$*"; mv amps "run_test_$*"; echo "test_$*_rundir... done" >> test_$*.diff
test_%_run:
\t@echo "test_$*_run..." >> test_$*.diff; echo "running $* in run_test_$*"; echo "test_$*_run... done" >> test_$*.diff
test_Mars-ions_check:
\t@: > test_Mars-ions.diff
test_Mars-ions--valgrind_check:
\t@echo "memcheck: 1 distinct error(s); amps exit status 0" > test_Mars-ions--valgrind.diff
"""


class Runner(unittest.TestCase):
    """makefile_runner.run on a workspace with a fake Makefile: make names, shown names, log files."""

    def test_rows_and_logs_of_a_variant(self):
        with tempfile.TemporaryDirectory() as tmp:
            ws, logs = Path(tmp, "ws"), Path(tmp, "logs")
            (ws / "MakefileTest").mkdir(parents=True)
            (ws / "MakefileTest" / "Table").write_text(TABLE)
            (ws / "Makefile").write_text(MAKEFILE)
            items = [plan_mod.makefile_item("Mars-ions", "planets", ""),
                     plan_mod.makefile_item("Mars-ions--valgrind", "planets", "")]
            env = {"PATH": "/usr/local/bin:/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp")}
            rows = makefile_runner.run(ws, machine(), items, env, logs, "t", time.time() + 300, lambda message: None)
            got = {r["test"]: r for r in rows}
            self.assertEqual(sorted(got), ["Mars-ions", "Mars-ions(valgrind)"])
            self.assertEqual(got["Mars-ions"]["status"], "PASS")
            self.assertEqual(got["Mars-ions(valgrind)"]["status"], "FAILED")
            self.assertIn("memcheck: 1 distinct error(s)", got["Mars-ions(valgrind)"]["reason"])
            log = Path(got["Mars-ions(valgrind)"]["log_src"])
            self.assertEqual(log.name, "Mars-ions(valgrind).log")
            text = log.read_text()
            self.assertIn("test_Mars-ions--valgrind_compile", text)              # make ran the Table name
            self.assertIn("running Mars-ions--valgrind in run_test_Mars-ions--valgrind", text)
            self.assertTrue((ws / "run_test_Mars-ions--valgrind" / "amps").exists())
            self.assertEqual(Path(got["Mars-ions"]["log_src"]).name, "Mars-ions.log")


if __name__ == "__main__":
    unittest.main()

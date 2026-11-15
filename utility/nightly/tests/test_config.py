import datetime as dt
import json
import tempfile
import unittest
from pathlib import Path

from amps_nightly import plan as plan_mod
from amps_nightly.checks import check_all
from amps_nightly.machines import MachineError, declared_columns, find_machine, load_machines, parse_machine
from amps_nightly.registry import Registry
from amps_nightly.runners import Routine, RoutineError

TODAY = dt.date(2026, 11, 14)
TABLE = """#>
Name=Alpha
Family=pic-field-solver
<#
#>
Name=Beta@@(0+passion)
Family=pic-field-solver
<#
#>
Name=Gamma@@(-passion)
Family=earth
<#
#>
Name=Retired@@(0+)
Family=earth
<#
"""


def machine(name="box", **extra):
    raw = {"schema_version": 1, "machine": name, "upload": "/srv/incoming",
           "configs": [{"config": "gnu", "config_pl": ["-install -compiler=gfortran,gcc_mpicc"]}]}
    raw.update(extra)
    return raw


def core(rid="core", match=(), **extra):
    """A runner entry for src/test/run_tests.py, as in runners.json."""
    raw = {"id": rid, "command": ["{python}", "src/test/run_tests.py", "--config", "{core_profile}"],
           "results": {"reader": "src-test-json", "file": "{out}/report.json"}, "families": "by-prefix",
           "expect": {"src_test_campaign": "nightly", "match": list(match)}, "requires_vars": ["core_profile"]}
    raw.update(extra)
    return raw


def routine(**extra):
    raw = {"schema_version": 1, "runners": [{"id": "legacy", "kind": "makefile"}, core()], "disabled": []}
    raw.update(extra)
    return Routine(raw)


class Machines(unittest.TestCase):
    def test_defaults_and_comments(self):
        m = parse_machine(machine(_comment="ignored"), Path("box.json"))
        self.assertEqual((m.table_host, m.hostnames, m.start, m.lanes), ("box", ["box"], "22:30", 1))
        self.assertEqual(m.configs[0].column, "box-gnu")
        self.assertEqual(m.configs[0].compiler, "GNU")
        self.assertIsNone(m.configs[0].sources)

    def test_errors(self):
        for raw, path, text in ((machine(), "other.json", "file is named"),
                                (machine(enable=True), "box.json", "unknown key enable"),
                                (machine(start="10pm"), "box.json", "HH:MM"),
                                (machine(timeouts_s={"link": 5}), "box.json", "timeouts_s"),
                                (machine(configs=[]), "box.json", "non-empty"),
                                (machine(name="Box"), "Box.json", "lower case")):
            with self.assertRaises(MachineError) as caught:
                parse_machine(raw, Path(path))
            self.assertIn(text, str(caught.exception))

    def test_switching_off_hides_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "a.json").write_text(json.dumps(machine("a", hostnames=["a.example.org"])))
            two = machine("b")
            two["configs"].append({"config": "intel", "config_pl": ["-install"], "enabled": False})
            Path(tmp, "b.json").write_text(json.dumps(two))
            Path(tmp, "c.json").write_text(json.dumps(machine("c", enabled=False)))
            Path(tmp, "bad.json").write_text("{")
            machines, problems = load_machines(Path(tmp))
        self.assertEqual([m.machine for m in machines], ["a", "b", "c"])
        self.assertEqual(len(problems), 1)
        self.assertEqual([c["column"] for c in declared_columns(machines)], ["a-gnu", "b-gnu"])
        self.assertEqual(find_machine(machines, hostname="a.example.org").machine, "a")
        self.assertEqual(find_machine(machines, hostname="b.example.org").machine, "b")   # short name


class RoutineFile(unittest.TestCase):
    def test_repository_files_are_consistent(self):
        errors, _warnings, summary = check_all(today=TODAY)
        self.assertEqual(errors, [])
        self.assertTrue(any("passion-gnu" in line for line in summary))

    def test_disabled_entries(self):
        r = routine(disabled=[{"test": "SEP-*", "reason": "new reference pending", "until": "2026-11-30",
                               "machines": ["passion"]},
                              {"test": "Alpha", "source": "legacy", "reason": "old", "until": "2026-11-01"}])
        hit = r.disabled_entry("legacy", "SEP-CME", "passion", "gnu", TODAY)
        self.assertEqual(hit.describe(), "disabled: new reference pending (until 2026-11-30)")
        self.assertIsNone(r.disabled_entry("legacy", "SEP-CME", "pleiades", "gnu", TODAY))
        self.assertIsNone(r.disabled_entry("legacy", "SEP-CME", "passion", "gnu", dt.date(2026, 12, 1)))
        self.assertIsNone(r.disabled_entry("legacy", "Alpha", "passion", "gnu", TODAY))   # expired

    def test_errors(self):
        bad_runners = ([{"id": "x", "kind": "shell"}],
                       [{"id": "x", "command": ["true"], "results": {"format": "exit-code"}, "families": "f"}],
                       [{"id": "x", "command": ["true"], "results": {"reader": "exit-code"}}],      # no families
                       [{"id": "x", "kind": "makefile"}, {"id": "x", "kind": "makefile"}],
                       [core(workspace="shared")],
                       [core(campaign="nightly")])                                                 # unknown key
        for runners in bad_runners:
            with self.assertRaises(RoutineError):
                Routine({"schema_version": 1, "runners": runners})
        with self.assertRaises(RoutineError):
            routine(disabled=[{"test": "A"}])                       # a reason is required
        with self.assertRaises(RoutineError):
            routine(disabled=[{"test": "A", "reason": "r", "source": "nope"}])


class Plan(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.table = Path(self.tmp.name, "Table")
        self.table.write_text(TABLE)
        self.manifest = Path(self.tmp.name, "manifest.json")
        self.manifest.write_text(json.dumps({"campaigns": ["smoke", "nightly"], "tests": [
            {"id": "CORE-PIC-001", "tiers": ["nightly"]}, {"id": "CORE-MESH-001", "tiers": ["smoke"]},
            {"id": "APPLICATION-MOON-001", "tiers": ["nightly"]}]}))
        self.registry = Registry.load()

    def tearDown(self):
        self.tmp.cleanup()

    def run_plan(self, m, r):
        plans = plan_mod.plan(m, m.configs[0], r, self.registry, TODAY, table=self.table, manifest=self.manifest)
        return {sp.source.id: sp for sp in plans}

    def test_host_selection_follows_table_parser(self):
        self.assertTrue(plan_mod.host_selected("", "box"))
        self.assertTrue(plan_mod.host_selected("(0+passion+box)", "box"))
        self.assertFalse(plan_mod.host_selected("(0+)", "box"))
        self.assertFalse(plan_mod.host_selected("(-box)", "box"))
        self.assertTrue(plan_mod.host_selected("(-pleiades)", "box"))

    def test_plan(self):
        m = parse_machine(machine("passion", core_profile="test/qualification-machine/passion.json"), Path("passion.json"))
        r = routine(disabled=[{"test": "Beta", "reason": "flaky"}])
        got = self.run_plan(m, r)
        legacy = {i.test: (i.action, i.reason) for i in got["legacy"].items}
        self.assertEqual(legacy, {"Alpha": ("run", ""), "Beta": ("skip", "disabled: flaky")})
        core = {i.test: i.family for i in got["core"].items}
        self.assertEqual(core, {"CORE-PIC-001": "core-pic", "APPLICATION-MOON-001": "app-moon"})

    def test_new_machine_without_profile_or_capability(self):
        m = parse_machine(machine("box", table_host="box"), Path("box.json"))
        r = Routine({"schema_version": 1, "runners": [
            {"id": "legacy", "kind": "makefile"},
            core(match=["CORE-*"]),
            core("moon", match=["APPLICATION-MOON-*"], requires=["spice"], requires_vars=[]),
            {"id": "sep", "families": "app-sep", "command": ["true"], "results": {"reader": "exit-code"},
             "enabled": False}]})
        got = self.run_plan(m, r)
        self.assertEqual(sorted(i.test for i in got["legacy"].items), ["Alpha", "Gamma"])
        self.assertEqual([(i.action, i.reason) for i in got["core"].items],
                         [("skip", "box does not define core_profile in its machine file")])
        self.assertIn("capability spice", got["moon"].items[0].reason)
        self.assertFalse(got["sep"].active)

    def test_config_sources_and_duplicates(self):
        m = parse_machine(machine("box", core_profile="p.json"), Path("box.json"))
        m.configs[0].sources = ["core"]
        r = Routine({"schema_version": 1, "runners": [
            {"id": "legacy", "kind": "makefile"}, core(), core("pic", match=["CORE-PIC-*"])]})
        m.configs[0].sources = ["core", "pic"]
        got = self.run_plan(m, r)
        self.assertFalse(got["legacy"].active)
        self.assertEqual(plan_mod.duplicates(list(got.values())), ["core-pic / CORE-PIC-001 comes from both core and pic"])


if __name__ == "__main__":
    unittest.main()

import contextlib
import datetime as dt
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from amps_nightly import cli
from amps_nightly import plan as plan_mod
from amps_nightly.checks import check_all
from amps_nightly.machines import (MachineError, declared_columns, declared_routines, find_machine, load_machines,
                                   parse_machine, parse_start)
from amps_nightly.registry import Registry
from amps_nightly.runners import Routine, RoutineError, routine_names, routine_path

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
#>
Name=Event-long--Parker
Family=sep-long
Tier=long
Ranks=32
<#
#>
Name=Wide-long@@(0+pleiades)
Family=sep-long
Tier=long
Ranks=64
<#
"""
LONG = {"enabled": True, "start": "Sat 10:00", "time_limit_h": 32, "lanes": 1, "compile_lanes": 1,
        "timeouts_s": {"compile": 3600, "rundir": 600, "run": 86400, "check": 1800},
        "make_vars": {"TESTMPILAUNCH": "mpirun -np"}, "max_ranks": 16, "stale_after_h": 192}


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
                                (machine(name="Box"), "Box.json", "lower case"),
                                (machine(compile_lanes=0), "box.json", "compile_lanes must be a positive integer"),
                                (machine(compile_jobs="all"), "box.json", "compile_jobs must be a positive integer"),
                                (machine(configs=[{"config": "gnu", "config_pl": ["-install"], "compile_jobs": 0}]),
                                 "box.json", "configs[1]: compile_jobs must be a positive integer"),
                                (machine(make_vars={"TESTCOMPILEJOBS": "8"}), "box.json", "compile_jobs")):
            with self.assertRaises(MachineError) as caught:
                parse_machine(raw, Path(path))
            self.assertIn(text, str(caught.exception))

    def test_compile_jobs(self):
        """make -j of a Makefile.test compile (TESTCOMPILEJOBS): the cores shared by the compile lanes and by the
        runners that may run at once, unless the machine or the configuration sets compile_jobs."""
        def jobs(config=None, **extra):
            m = parse_machine(machine(**extra), Path("box.json"))
            return m.compile_jobs_for(m.configs[0] if config is None else config)
        self.assertEqual(jobs(cores=36, compile_lanes=4), 9)
        self.assertEqual(jobs(cores=36), 36)
        self.assertEqual(jobs(cores=36, compile_lanes=4, max_parallel=2), 4)
        self.assertEqual(jobs(cores=2, compile_lanes=4), 1)
        self.assertEqual(jobs(cores=36, compile_lanes=4, compile_jobs=12), 12)
        m = parse_machine(machine(cores=36, compile_lanes=4, compile_jobs=12, configs=[
            {"config": "gnu", "config_pl": ["-install"]},
            {"config": "intel", "config_pl": ["-install"], "compile_jobs": 6}]), Path("box.json"))
        self.assertEqual([m.compile_jobs_for(c) for c in m.configs] + [m.compile_jobs_for()], [12, 6, 12])
        self.assertEqual(parse_machine(machine(), Path("box.json")).compile_lanes, 1)

    def test_routines(self):
        """A routine of "routines" replaces the machine's start, limit and lanes and merges timeouts and make_vars."""
        m = parse_machine(machine(cores=36, lanes=4, compile_lanes=4, start="20:00", night_limit_h=14,
                                  timeouts_s={"run": 7200}, make_vars={"MPIRUN": "mpirun -np 8"},
                                  routines={"long": LONG, "_note": "ignored"}), Path("box.json"))
        self.assertEqual((m.runs_routine("nightly"), m.runs_routine("long"), m.runs_routine("weekly")), (True, True, False))
        self.assertIs(m.for_routine("nightly"), m)
        long_ = m.for_routine("long")
        self.assertEqual((long_.routine, long_.start, long_.night_limit_h, long_.lanes, long_.compile_lanes),
                         ("long", "Sat 10:00", 32, 1, 1))
        self.assertEqual(long_.timeouts_s, {"compile": 3600, "rundir": 600, "run": 86400, "check": 1800})
        self.assertEqual(long_.make_vars, {"MPIRUN": "mpirun -np 8", "TESTMPILAUNCH": "mpirun -np"})
        self.assertEqual(long_.compile_jobs_for(), 36)                    # one compile lane gets the cores
        self.assertEqual((m.lanes, m.start, m.timeouts_s["run"], m.compile_jobs_for()), (4, "20:00", 7200, 9))
        self.assertEqual((m.max_ranks("long"), m.stale_after_h("long")), (16, 192))
        # defaults: the machine's values, max_ranks = cores, stale after 192 h
        d = parse_machine(machine(cores=36, routines={"long": {"start": "12:00"}}), Path("box.json")).for_routine("long")
        self.assertEqual((d.lanes, d.night_limit_h, d.max_ranks("long"), d.stale_after_h("long")), (1, 14, 36, 192))
        off = parse_machine(machine(routines={"long": dict(LONG, enabled=False)}), Path("box.json"))
        self.assertFalse(off.runs_routine("long"))
        self.assertEqual(off.enabled_routines(), [])
        with self.assertRaises(MachineError):
            off.for_routine("long")
        self.assertEqual(parse_machine(machine(), Path("box.json")).enabled_routines(), [])

    def test_routine_errors(self):
        for routines, text in (({"long": dict(LONG, lane=1)}, "routines.long: unknown key lane"),
                               ({"nightly": LONG}, "the nightly routine is set by the machine's own keys"),
                               ({"Long": LONG}, "a routine name must be lower case"),
                               ({"long": dict(LONG, start="Saturday 10:00")}, "must be HH:MM or <day> HH:MM"),
                               ({"long": {k: v for k, v in LONG.items() if k != "start"}}, "start must be a non-empty"),
                               ({"long": dict(LONG, time_limit_h=200)}, "time_limit_h must be a number of hours"),
                               ({"long": dict(LONG, lanes=0)}, "lanes must be a positive integer"),
                               ({"long": dict(LONG, max_ranks="32")}, "max_ranks must be a positive integer"),
                               ({"long": dict(LONG, stale_after_h=0)}, "stale_after_h must be a number of hours"),
                               ({"long": dict(LONG, timeouts_s={"link": 5})}, "timeouts_s may only set"),
                               ({"long": dict(LONG, make_vars={"TESTCOMPILEJOBS": "4"})}, "make_vars may not set"),
                               ({"long": dict(LONG, enabled="yes")}, "enabled must be true or false"),
                               ([LONG], "routines must map routine names to objects")):
            with self.assertRaises(MachineError) as caught:
                parse_machine(machine(routines=routines), Path("box.json"))
            self.assertIn(text, str(caught.exception))

    def test_start_times(self):
        self.assertEqual([parse_start(t) for t in ("Sat 10:00", "Sun 06:30", "20:00", "9:05")],
                         [(6, 10, 0), (0, 6, 30), (None, 20, 0), (None, 9, 5)])
        for bad in ("Saturday 10:00", "sat 10:00", "24:00", "10"):
            with self.assertRaises(ValueError):
                parse_start(bad)

    def test_repository_passion_routines(self):
        """passion runs the long routine on Saturday at 10:00, when the Friday nightly (20:00 + 14 h) has reached its
        limit; its 32 h end by Sunday 18:00, 2 h before the Sunday nightly (20:00)."""
        raw = json.loads((Path(__file__).resolve().parents[1] / "machines" / "passion.json").read_text())
        self.assertEqual(raw["routines"], {"long": LONG})
        self.assertEqual((raw["start"], raw["night_limit_h"], raw["cores"]), ("20:00", 14, 36))
        friday_end = dt.datetime(2026, 11, 27, 20, 0) + dt.timedelta(hours=raw["night_limit_h"])
        start = dt.datetime(2026, 11, 28, 10, 0)
        self.assertLessEqual(friday_end, start)
        self.assertEqual(start + dt.timedelta(hours=LONG["time_limit_h"]), dt.datetime(2026, 11, 29, 18, 0))
        machines, _ = load_machines()
        self.assertEqual([c["column"] for c in declared_routines(machines)["long"]], ["passion-gnu"])

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

    def test_routine_files(self):
        """runners.json is the nightly; runners.<routine>.json next to it is another routine, with the same schema."""
        with tempfile.TemporaryDirectory() as tmp:
            nightly = Path(tmp, "runners.json")
            nightly.write_text(json.dumps({"schema_version": 1, "runners": [{"id": "legacy", "kind": "makefile"}]}))
            Path(tmp, "runners.long.json").write_text(json.dumps({"schema_version": 1, "runners": [
                {"id": "legacy-long", "kind": "makefile", "manifest_target": "test_long_manifest"}]}))
            Path(tmp, "runners.Bad.json").write_text("{}")                 # not a routine name: ignored
            self.assertEqual(routine_names(nightly), ["nightly", "long"])
            self.assertEqual(routine_path("long", nightly), Path(tmp, "runners.long.json"))
            self.assertEqual(routine_path("nightly", nightly), nightly)
            long_ = Routine.load_routine("long", nightly)
            self.assertEqual((long_.name, long_.label), ("long", "runners.long.json"))
            self.assertEqual(long_.sources[0].manifest_target, "test_long_manifest")
            self.assertEqual(Routine.load_routine("nightly", nightly).sources[0].manifest_target, "test_manifest")
            with self.assertRaises(RoutineError) as caught:
                Routine.load_routine("weekly", nightly)
            self.assertIn("no routine 'weekly': {} does not exist (routines: nightly, long)".format(
                Path(tmp, "runners.weekly.json")), str(caught.exception))
            Path(tmp, "runners.long.json").write_text(json.dumps({"schema_version": 1, "runners": [
                {"id": "x", "kind": "makefile", "manifest_target": "test_all"}]}))
            with self.assertRaises(RoutineError) as caught:
                Routine.load_routine("long", nightly)
            self.assertIn("runners.long.json: runner x: manifest_target must be one of test_long_manifest, "
                          "test_manifest", str(caught.exception))
        with self.assertRaises(RoutineError) as caught:
            Routine({"schema_version": 1, "runners": [core(manifest_target="test_manifest")]})
        self.assertIn("manifest_target is only for kind makefile", str(caught.exception))

    def test_repository_long_routine(self):
        long_ = Routine.load_routine("long")
        self.assertEqual([(s.id, s.kind, s.enabled) for s in long_.sources],
                         [("sep-event-inputs-long", "command", True), ("legacy-long", "makefile", True)])
        self.assertEqual(long_.source("legacy-long").manifest_target, "test_long_manifest")
        self.assertEqual(long_.checkout.repository, Routine.load().checkout.repository)
        self.assertEqual(long_.source("legacy-long").time_limit_min, 32 * 60)   # passion's time_limit_h

    def test_repository_sep_event_input_check(self):
        """Both routines regenerate the 2013-04-11 SEP event inputs from the raw archives (prepare_event.py check),
        with the paths of passion's machine file; a dry run describes the file-less exit-code reader."""
        machines, _ = load_machines()
        passion = find_machine(machines, "passion")
        for name, rid, family in (("nightly", "sep-event-inputs", "app-sep"),
                                  ("long", "sep-event-inputs-long", "sep-long")):
            source = Routine.load_routine(name).source(rid)
            self.assertIsNotNone(source, rid)
            self.assertEqual((source.kind, source.workspace, source.family, source.results),
                             ("command", "plain", family, {"reader": "exit-code"}))
            self.assertEqual(source.command[2:5], ["srcSEP/test/event_2013_04_11/prepare_event.py", "check", "--raw"])
            self.assertEqual(source.env, {"PYTHONPATH": "{sep_python_path}"})
            for var in source.requires_vars:
                self.assertTrue(passion.vars.get(var), var)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = cli.main(["run", "--machine", "passion", "--routine", name, "--only", rid, "--dry-run"])
            self.assertEqual(code, 0)
            self.assertIn("prepare_event.py check --raw /data/vtenishe/sep_validation_data/event_2013_04_11",
                          out.getvalue())
            self.assertIn("results: exit-code\n", out.getvalue() + "\n")

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


class LongChecks(unittest.TestCase):
    """check: the routine files, the Ranks= of Tier=long tests against the machines that run them, and more."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        t = self.dir = Path(self.tmp.name)
        (t / "Table").write_text(TABLE)
        (t / "machines").mkdir()
        (t / "manifest.json").write_text(json.dumps({"campaigns": ["nightly"], "tests": []}))
        self.write_runners([{"id": "legacy", "kind": "makefile"}],
                           [{"id": "legacy-long", "kind": "makefile", "manifest_target": "test_long_manifest"}])

    def tearDown(self):
        self.tmp.cleanup()

    def write_runners(self, nightly, long_=None):
        (self.dir / "runners.json").write_text(json.dumps({"schema_version": 1, "runners": nightly}))
        if long_ is not None:
            (self.dir / "runners.long.json").write_text(json.dumps({"schema_version": 1, "runners": long_}))
        elif (self.dir / "runners.long.json").exists():
            (self.dir / "runners.long.json").unlink()

    def write_machine(self, name, **extra):
        raw = machine(name, table_host=name, **extra)
        (self.dir / "machines" / (name + ".json")).write_text(json.dumps(raw))

    def check(self):
        return check_all(self.dir / "Table", None, self.dir / "runners.json", self.dir / "machines",
                         self.dir / "manifest.json", today=TODAY)

    def test_ranks_against_max_ranks(self):
        self.write_machine("passion", cores=36, routines={"long": dict(LONG, max_ranks=16)})
        self.write_machine("pleiades", cores=8, routines={"long": {"start": "Sun 01:00"}})     # max_ranks = cores
        self.write_machine("small", cores=4)                                                 # no long routine
        errors, _, summary = self.check()
        self.assertEqual(sorted(e for e in errors if "max_ranks" in e), [
            "MakefileTest/Table line {} (Event-long--Parker): Ranks=32 exceeds max_ranks 16 of routine long on passion "
            "(machines/passion.json)".format(TABLE.splitlines().index("Name=Event-long--Parker")),
            "MakefileTest/Table line {} (Event-long--Parker): Ranks=32 exceeds max_ranks 8 of routine long on pleiades "
            "(machines/pleiades.json)".format(TABLE.splitlines().index("Name=Event-long--Parker")),
            "MakefileTest/Table line {} (Wide-long): Ranks=64 exceeds max_ranks 8 of routine long on pleiades "
            "(machines/pleiades.json)".format(TABLE.splitlines().index("Name=Wide-long@@(0+pleiades)"))])
        self.assertIn("runners.long.json: 1 runner, enabled: legacy-long", summary)
        self.assertIn("routine long (runners.long.json): runs on passion (Sat 10:00, 32 h), pleiades (Sun 01:00, 14 h)",
                      summary)
        self.assertIn("MakefileTest/Table: 6 blocks, 6 tests, 2 of them Tier=long", summary)

    def test_routine_files_and_machines(self):
        self.write_machine("passion", routines={"weekly": LONG})
        self.write_runners([{"id": "legacy", "kind": "makefile"}],
                           [{"id": "legacy", "kind": "makefile", "manifest_target": "test_long_manifest"}])
        errors, _, _ = self.check()
        self.assertIn("machines/passion.json: routines.weekly: there is no runners.weekly.json next to runners.json",
                      errors)
        self.assertIn("runners.long.json runner legacy: the id is also a runner of runners.json; runner ids must "
                      "differ between routines (a runner's workspace is <workdir>/<COMPILER>/<runner>)", errors)
        (self.dir / "runners.long.json").write_text('{"schema_version": 1, "runners": [{"id": "x", "kind": "y"}]}')
        errors, _, _ = self.check()
        self.assertIn("runners.long.json: runner x: kind must be one of makefile, command", errors)

    def test_long_tests_without_a_long_routine(self):
        self.write_runners([{"id": "legacy", "kind": "makefile"}])
        errors, warnings, _ = self.check()
        self.assertEqual(errors, [])
        self.assertIn("MakefileTest/Table has Tier=long tests, but no routine file has a makefile runner with "
                      "\"manifest_target\": \"test_long_manifest\": they never run", warnings)


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

    def configured_tree(self, host="passion"):
        """A real tiny manifest isolates precedence from the developer checkout.

        This is administrative selection data, not model output or observation
        data. make reads the declared manifest normally, so both manifest use
        and Table use are tested at the public planning boundary.
        """
        tree = Path(self.tmp.name, "configured")
        tree.mkdir()
        (tree / "Makefile").write_text("include Makefile.test\n")
        (tree / "Makefile.test").write_text(
            "HOSTNAME={}\ntest_manifest:\n\t@echo 'Configured|pic-field-solver|||amps'\n".format(host))
        return tree

    def test_explicit_table_wins_over_implicitly_configured_tree(self):
        machine_value = parse_machine(machine("passion"), Path("passion.json"))
        tree = self.configured_tree()
        self.assertEqual(plan_mod.configured_host(tree), "passion")
        with mock.patch.object(plan_mod, "REPO_ROOT", tree):
            items, note = plan_mod.legacy_items(machine_value, None, self.table)
        self.assertEqual([item.test for item in items], ["Alpha", "Beta"])
        self.assertIn("from MakefileTest/Table", note)

    def test_explicit_configured_tree_keeps_priority_over_table(self):
        machine_value = parse_machine(machine("passion"), Path("passion.json"))
        tree = self.configured_tree()
        items, note = plan_mod.legacy_items(machine_value, tree, self.table)
        self.assertEqual([item.test for item in items], ["Configured"])
        self.assertIn("exact list from 'make test_manifest'", note)

    def test_explicit_configured_tree_wrong_host_does_not_fall_back_to_table(self):
        machine_value = parse_machine(machine("passion"), Path("passion.json"))
        tree = self.configured_tree(host="another-host")
        with self.assertRaisesRegex(RuntimeError, "configured for host 'another-host', not 'passion'"):
            plan_mod.legacy_items(machine_value, tree, self.table)

    def test_implicit_configured_tree_keeps_exact_manifest_without_table(self):
        machine_value = parse_machine(machine("passion"), Path("passion.json"))
        tree = self.configured_tree()
        with mock.patch.object(plan_mod, "REPO_ROOT", tree):
            items, note = plan_mod.legacy_items(machine_value, None, None)
        self.assertEqual([item.test for item in items], ["Configured"])
        self.assertIn("exact list from 'make test_manifest'", note)

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

    def test_plan_of_each_tier(self):
        """The nightly makefile runner lists the nightly tier, the long routine's only the Tier=long blocks."""
        m = parse_machine(machine("passion", routines={"long": LONG}), Path("passion.json"))
        got = self.run_plan(m, routine())
        self.assertEqual(sorted(i.test for i in got["legacy"].items), ["Alpha", "Beta"])
        long_ = Routine({"schema_version": 1, "runners": [
            {"id": "legacy-long", "kind": "makefile", "manifest_target": "test_long_manifest", "enabled": True}]}, "long")
        got = self.run_plan(m.for_routine("long"), long_)
        self.assertEqual([(i.test, i.family) for i in got["legacy-long"].items], [("Event-long--Parker", "sep-long")])
        self.assertIn("Tier=long blocks", got["legacy-long"].note)
        off = Routine({"schema_version": 1, "runners": [{"id": "legacy-long", "kind": "makefile", "enabled": False}]},
                      "long")
        self.assertEqual(self.run_plan(m, off)["legacy-long"].note, "switched off in runners.long.json")

    def test_configured_tree_long_manifest(self):
        """With a configured tree, the long routine's list is 'make -s test_long_manifest' of Makefile.test.long."""
        m = parse_machine(machine("passion"), Path("passion.json"))
        tree = self.configured_tree()
        self.assertEqual(plan_mod.configured_host(tree, "test_long_manifest"), "")      # no Makefile.test.long yet
        (tree / "Makefile").write_text("include Makefile.test\ninclude Makefile.test.long\n")
        (tree / "Makefile.test.long").write_text(
            "TESTMPILAUNCH?=mpirun -np\n\ntest_long_manifest:\n\t@echo 'Configured-long|sep-long|||amps'\n")
        self.assertEqual(plan_mod.configured_host(tree, "test_long_manifest"), "passion")
        items, note = plan_mod.legacy_items(m, tree, self.table, "test_long_manifest")
        self.assertEqual([(i.test, i.family) for i in items], [("Configured-long", "sep-long")])
        self.assertIn("exact list from 'make test_long_manifest'", note)
        items, _ = plan_mod.legacy_items(m, tree, self.table)
        self.assertEqual([i.test for i in items], ["Configured"])

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

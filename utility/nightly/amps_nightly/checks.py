"""`amps-nightly check`: validate the Table, families.json, the routine files (runners.json, runners.<routine>.json)
and machines/*.json together."""

from __future__ import annotations

import datetime as dt
import fnmatch
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .machines import REPO_ROOT, load_machines
from .orchestrator import BUILTIN_PLACEHOLDERS, placeholders
from .plan import CORE_MANIFEST, host_selected, machine_vars
from .readers import READERS
from .registry import Registry, RegistryError
from .runners import NIGHTLY, Routine, RoutineError, file_name, routine_names
from .table import MANIFESTS, check_table, read_table


def check_all(table: Optional[Path] = None, families: Optional[Path] = None, routine_path: Optional[Path] = None,
              machines_dir: Optional[Path] = None, manifest: Optional[Path] = None,
              today: Optional[dt.date] = None) -> Tuple[List[str], List[str], List[str]]:
    """(errors, warnings, summary lines). routine_path is runners.json; every runners.<routine>.json next to it is
    checked as well."""
    today = today or dt.date.today()
    errors: List[str] = []
    warnings: List[str] = []
    summary: List[str] = []
    try:
        registry = Registry.load(families)
    except RegistryError as error:
        return [str(error)], [], []
    try:
        routines = [Routine.load(routine_path)]
    except RoutineError as error:
        return [str(error)], [], []
    for name in routine_names(routine_path)[1:]:
        try:
            routines.append(Routine.load_routine(name, routine_path))
        except RoutineError as error:
            errors.append(str(error))

    tests = read_table(table)
    errors += ["MakefileTest/Table: " + p for p in check_table(tests, registry)]
    legacy_names = {t.name for t in tests}
    long_tests = [t for t in tests if t.long]
    summary.append("MakefileTest/Table: {} blocks, {} tests{}".format(
        len(tests), len(legacy_names), ", {} of them Tier=long".format(len({t.name for t in long_tests}))
        if long_tests else ""))

    try:
        core = json.loads(Path(manifest or CORE_MANIFEST).read_text())
        campaigns: List[str] = core.get("campaigns", [])
        suites: Dict[str, List[str]] = {t["id"]: t.get("campaigns", t.get("tiers", [])) for t in core["tests"]}
    except (OSError, ValueError, KeyError) as error:
        warnings.append("cannot read the src/test manifest ({}); core sources not checked".format(error))
        campaigns, suites = [], {}

    owner: Dict[str, str] = {}               # runner id -> the routine file that has it
    for routine in routines:
        for source in routine.sources:
            where = "{} runner {}".format(routine.label, source.id)
            if source.id in owner:
                errors.append("{}: the id is also a runner of {}; runner ids must differ between routines (a "
                              "runner's workspace is <workdir>/<COMPILER>/<runner>)".format(where, owner[source.id]))
            owner.setdefault(source.id, routine.label)
            if source.family and registry.family(source.family) is None:
                errors.append("{}: family {!r} is not defined in families.json".format(where, source.family))
            if source.kind == "command" and source.results.get("reader") not in READERS:
                errors.append("{}: unknown reader {!r} (known: {})".format(where, source.results.get("reader"),
                                                                         ", ".join(sorted(READERS))))
            for name in sorted(placeholders(source) - BUILTIN_PLACEHOLDERS - set(source.requires_vars)):
                errors.append("{}: uses {{{}}}, which is not built in; list it in requires_vars so machines "
                              "without it skip the runner".format(where, name))
            if source.campaign and suites:
                if source.campaign not in campaigns:
                    errors.append("{}: src/test has no campaign {!r}".format(where, source.campaign))
                for glob in source.match:
                    if not any(fnmatch.fnmatchcase(s, glob) and source.campaign in c for s, c in suites.items()):
                        warnings.append("{}: match {!r} selects no suite of campaign {}".format(
                            where, glob, source.campaign))
        enabled_sources = [s.id for s in routine.sources if s.enabled]
        summary.append("{}: {} runner{}, enabled: {}".format(routine.label, len(routine.sources),
                                                              "" if len(routine.sources) == 1 else "s",
                                                              ", ".join(enabled_sources) or "none"))
    # the routines whose makefile runners list a tier of the Table (test_manifest: nightly, test_long_manifest: long)
    tier_routines = {(MANIFESTS[s.manifest_target][0], r.name) for r in routines for s in r.sources
                     if s.kind == "makefile"}
    if long_tests and not any(tier == "long" for tier, _ in tier_routines):
        warnings.append("MakefileTest/Table has Tier=long tests, but no routine file has a makefile runner with "
                        "\"manifest_target\": \"test_long_manifest\": they never run")

    machines, problems = load_machines(machines_dir)
    errors += ["machines/" + p for p in problems]
    known_hosts: Dict[str, str] = {}
    names = {r.name for r in routines}
    for machine in machines:
        where = "machines/{}".format(machine.path.name)
        for host in machine.hostnames:
            if host in known_hosts:
                errors.append("{}: hostname {} is also listed by {}".format(where, host, known_hosts[host]))
            known_hosts.setdefault(host, machine.path.name)
        for name in sorted(machine.routines or {}):
            if name not in names:
                errors.append("{}: routines.{}: there is no {} next to runners.json".format(where, name, file_name(name)))
        for routine in routines:
            if not machine.runs_routine(routine.name):
                continue
            for source in routine.sources:
                missing = [v for v in source.requires_vars if v not in machine_vars(machine)]
                if source.enabled and missing and machine.enabled:
                    warnings.append("{}: does not define {} (runner {}{} will be reported as skipped)".format(
                        where, ", ".join(missing), source.id,
                        "" if routine.name == NIGHTLY else " of " + routine.label))
        if machine.core_profile and not (REPO_ROOT / machine.core_profile).is_file():
            errors.append("{}: core_profile {} does not exist".format(where, machine.core_profile))
        for config in machine.configs:
            for sid in config.sources or []:
                if not any(r.source(sid) for r in routines):
                    errors.append("{} config {}: unknown source {!r}".format(where, config.config, sid))
        if machine.enabled and not any(machine.table_host in t.selector for t in tests):
            warnings.append("{}: no @@(...) selector in MakefileTest/Table names {!r}; tests limited to other "
                            "machines will not run there and machine-specific Ref=/Keys= lines will not apply"
                            .format(where, machine.table_host))
        # a test of a fixed rank count (Ranks= of a Tier=long block) must fit the machines whose routine runs it
        for tier, name in sorted(tier_routines):
            if not machine.enabled or not machine.runs_routine(name):
                continue
            limit = machine.max_ranks(name)
            for test in tests:
                if test.tier == tier and test.ranks and test.ranks > limit and \
                        host_selected(test.selector, machine.table_host):
                    errors.append("MakefileTest/Table line {} ({}): Ranks={} exceeds max_ranks {} of routine {} on {} "
                                  "({})".format(test.line, test.name, test.ranks, limit, name, machine.machine, where))
    columns = [c.column for m in machines for c in m.active_configs()]
    summary.append("machines: {} file{}, {} enabled, report columns: {}".format(
        len(machines), "" if len(machines) == 1 else "s", sum(1 for m in machines if m.enabled),
        ", ".join(columns) or "none"))
    for routine in routines[1:]:
        hosts = ["{} ({}, {:g} h)".format(m.machine, m.routines[routine.name].start,
                                          m.for_routine(routine.name).night_limit_h)
                 for m in machines if m.enabled and m.runs_routine(routine.name)]
        summary.append("routine {} ({}): runs on {}".format(routine.name, routine.label, ", ".join(hosts) or
                                                          "no machine (no machine file enables it in \"routines\")"))

    machine_names = {m.machine: m for m in machines}
    active = 0
    for routine in routines:
        for entry in routine.disabled:
            where = "{} disabled[{}] ({})".format(routine.label, entry.number, entry.test)
            if entry.expired(today):
                warnings.append("{}: expired on {}; the test runs again, delete the entry".format(where, entry.until))
            for name in entry.machines:
                if name not in machine_names:
                    warnings.append("{}: unknown machine {!r}".format(where, name))
            for config in entry.configs:
                if not any(m.config(config) for m in machines):
                    warnings.append("{}: no machine has configuration {!r}".format(where, config))
            source = routine.source(entry.source) if entry.source else None
            if source is None or source.kind == "makefile" or source.campaign:
                candidates = (legacy_names if source is None or source.kind == "makefile" else set()) | \
                             (set(suites) if source is None or source.campaign else set())
                if not any(fnmatch.fnmatchcase(name, entry.test) for name in candidates):
                    warnings.append("{}: matches no MakefileTest/Table test or src/test suite".format(where))
        active += sum(1 for d in routine.disabled if not d.expired(today))
    summary.append("disabled tests: {} entr{}".format(active, "y" if active == 1 else "ies"))
    return errors, warnings, summary

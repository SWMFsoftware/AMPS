"""`amps-nightly check`: validate the Table, families.json, runners.json and machines/*.json together."""

from __future__ import annotations

import datetime as dt
import fnmatch
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .machines import REPO_ROOT, load_machines
from .orchestrator import BUILTIN_PLACEHOLDERS, placeholders
from .plan import CORE_MANIFEST, machine_vars
from .readers import READERS
from .registry import Registry, RegistryError
from .runners import Routine, RoutineError
from .table import check_table, read_table


def check_all(table: Optional[Path] = None, families: Optional[Path] = None, routine_path: Optional[Path] = None,
              machines_dir: Optional[Path] = None, manifest: Optional[Path] = None,
              today: Optional[dt.date] = None) -> Tuple[List[str], List[str], List[str]]:
    """(errors, warnings, summary lines)."""
    today = today or dt.date.today()
    errors: List[str] = []
    warnings: List[str] = []
    summary: List[str] = []
    try:
        registry = Registry.load(families)
    except RegistryError as error:
        return [str(error)], [], []
    try:
        routine = Routine.load(routine_path)
    except RoutineError as error:
        return [str(error)], [], []

    tests = read_table(table)
    errors += ["MakefileTest/Table: " + p for p in check_table(tests, registry)]
    legacy_names = {t.name for t in tests}
    summary.append("MakefileTest/Table: {} blocks, {} tests".format(len(tests), len(legacy_names)))

    try:
        core = json.loads(Path(manifest or CORE_MANIFEST).read_text())
        campaigns: List[str] = core.get("campaigns", [])
        suites: Dict[str, List[str]] = {t["id"]: t.get("campaigns", t.get("tiers", [])) for t in core["tests"]}
    except (OSError, ValueError, KeyError) as error:
        warnings.append("cannot read the src/test manifest ({}); core sources not checked".format(error))
        campaigns, suites = [], {}

    for source in routine.sources:
        where = "runners.json runner {}".format(source.id)
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
                    warnings.append("{}: match {!r} selects no suite of campaign {}".format(where, glob, source.campaign))
    enabled_sources = [s.id for s in routine.sources if s.enabled]
    summary.append("runners.json: {} runners, enabled: {}".format(len(routine.sources), ", ".join(enabled_sources) or "none"))

    machines, problems = load_machines(machines_dir)
    errors += ["machines/" + p for p in problems]
    known_hosts: Dict[str, str] = {}
    for machine in machines:
        where = "machines/{}".format(machine.path.name)
        for host in machine.hostnames:
            if host in known_hosts:
                errors.append("{}: hostname {} is also listed by {}".format(where, host, known_hosts[host]))
            known_hosts.setdefault(host, machine.path.name)
        for source in routine.sources:
            missing = [v for v in source.requires_vars if v not in machine_vars(machine)]
            if source.enabled and missing and machine.enabled:
                warnings.append("{}: does not define {} (runner {} will be reported as skipped)".format(
                    where, ", ".join(missing), source.id))
        if machine.core_profile and not (REPO_ROOT / machine.core_profile).is_file():
            errors.append("{}: core_profile {} does not exist".format(where, machine.core_profile))
        for config in machine.configs:
            for sid in config.sources or []:
                if routine.source(sid) is None:
                    errors.append("{} config {}: unknown source {!r}".format(where, config.config, sid))
        if machine.enabled and not any(machine.table_host in t.selector for t in tests):
            warnings.append("{}: no @@(...) selector in MakefileTest/Table names {!r}; tests limited to other "
                            "machines will not run there and machine-specific Ref=/Keys= lines will not apply"
                            .format(where, machine.table_host))
    columns = [c.column for m in machines for c in m.active_configs()]
    summary.append("machines: {} file{}, {} enabled, report columns: {}".format(
        len(machines), "" if len(machines) == 1 else "s", sum(1 for m in machines if m.enabled),
        ", ".join(columns) or "none"))

    names = {m.machine: m for m in machines}
    for entry in routine.disabled:
        where = "runners.json disabled[{}] ({})".format(entry.number, entry.test)
        if entry.expired(today):
            warnings.append("{}: expired on {}; the test runs again, delete the entry".format(where, entry.until))
        for name in entry.machines:
            if name not in names:
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
    active = [d for d in routine.disabled if not d.expired(today)]
    summary.append("disabled tests: {} entr{}".format(len(active), "y" if len(active) == 1 else "ies"))
    return errors, warnings, summary

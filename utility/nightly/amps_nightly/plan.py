"""What the nightly routine runs on one machine configuration (`amps-nightly plan`).

The runner uses the same resolution, so the plan is what the report will show for that column:
tests that run, and tests reported as SKIPPED with a reason (switched off in runners.json, marked
SKIP in MakefileTest/Table, or needing a capability the machine does not declare).
"""

from __future__ import annotations

import datetime as dt
import fnmatch
import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .machines import REPO_ROOT, Config, Machine
from .registry import Registry, UNASSIGNED
from .runners import Routine, Source
from .table import display_name, read_table, required_capabilities

RUN, SKIP = "run", "skip"
CORE_MANIFEST = REPO_ROOT / "src" / "test" / "manifest.json"


@dataclass
class Item:
    family: str
    test: str
    action: str           # RUN or SKIP
    reason: str = ""      # why it is skipped
    listed: bool = True   # False: a placeholder for a runner whose tests are known only after it ran
    target: str = ""      # Makefile.test name when it differs from the shown name (test_<target>_compile, ...)

    @property
    def make_name(self) -> str:
        return self.target or self.test


def makefile_item(name: str, family: str, code: str) -> Item:
    """The plan item of a Makefile.test test; a variant <test>--valgrind is shown as <test>(valgrind)."""
    shown = display_name(name)
    return Item(family or UNASSIGNED, shown, SKIP if code == "SKIP" else RUN,
                "ExeptionCode=SKIP in MakefileTest/Table for this compiler" if code == "SKIP" else "",
                target=name if shown != name else "")


@dataclass
class SourcePlan:
    source: Source
    active: bool          # False: the source does not run in this configuration at all
    note: str             # how the list was obtained, or why the source is inactive
    items: List[Item] = field(default_factory=list)

    def count(self, action: str) -> int:
        return sum(1 for i in self.items if i.action == action)


# ----------------------------------------------------------------------------- Table selection
def host_selected(selector: str, host: str) -> bool:
    """The @@(...) rule of utility/TestScripts/table_parser.cpp (IsThisHost), including its substring match."""
    if not selector:
        return True
    text = selector.replace("(", " ").replace(")", " ").replace(",", " ")
    return (host in text) if "0+" in text else (host not in text)


def makefile_manifest(amps_dir: Path, timeout: int = 300) -> List[Tuple[str, str, str]]:
    """(name, family, exception code) from `make -s test_manifest` of a configured tree."""
    done = subprocess.run(["make", "-s", "test_manifest"], cwd=str(amps_dir), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True, timeout=timeout)
    if done.returncode != 0:
        raise RuntimeError("make test_manifest failed in {}: {}".format(
            amps_dir, (done.stderr.strip().splitlines() or ["exit code {}".format(done.returncode)])[-1]))
    tests = []
    for line in done.stdout.splitlines():
        parts = line.split("|")
        if len(parts) == 5 and parts[0]:
            tests.append((parts[0], parts[1], parts[3]))
    return tests


def configured_host(amps_dir: Path) -> str:
    """HOSTNAME= of Makefile.test in a configured tree that has the test_manifest target, else ""."""
    makefile = amps_dir / "Makefile.test"
    try:
        text = makefile.read_text(errors="replace")
    except OSError:
        return ""
    if "\ntest_manifest:" not in text:
        return ""
    for line in text.splitlines():
        if line.startswith("HOSTNAME="):
            return line[len("HOSTNAME="):].strip()
    return ""


def legacy_items(machine: Machine, amps_dir: Optional[Path], table: Optional[Path]) -> Tuple[List[Item], str]:
    tree = amps_dir or REPO_ROOT
    # An explicit configured tree requests its exact compiler-selected manifest
    # and retains the host check below. Otherwise an explicit Table is the
    # caller's source of selection, even when the default checkout happens to
    # contain a Makefile.test from an earlier configuration. Only the implicit
    # default may use that cached manifest. This keeps --table and private
    # planning fixtures independent of unrelated local build state.
    if amps_dir or (table is None and configured_host(tree) == machine.table_host):
        if configured_host(tree) != machine.table_host:
            raise RuntimeError("{} was configured for host {!r}, not {!r}".format(
                tree, configured_host(tree) or "?", machine.table_host))
        items = [makefile_item(name, family, code) for name, family, code in makefile_manifest(tree)]
        return items, "exact list from 'make test_manifest' in {}".format(tree)
    seen: Dict[str, Item] = {}
    for test in read_table(table):
        if host_selected(test.selector, machine.table_host) and test.name not in seen:
            seen[test.name] = makefile_item(test.name, test.family, "")
    return (list(seen.values()),
            "from MakefileTest/Table for host {!r}; compiler-specific SKIP rules are applied at run time".format(
                machine.table_host))


# ----------------------------------------------------------------------------- src/test selection
def machine_vars(machine: Machine) -> Dict[str, str]:
    """Machine values for {placeholders}; core_profile counts as a variable."""
    values = dict(machine.vars or {})
    if machine.core_profile:
        values["core_profile"] = machine.core_profile
    return values


def core_suites(source: Source, manifest: Optional[Path] = None) -> List[str]:
    data = json.loads(Path(manifest or CORE_MANIFEST).read_text())
    if source.campaign not in data.get("campaigns", []):
        raise RuntimeError("src/test has no campaign {!r} (it has: {})".format(
            source.campaign, ", ".join(data.get("campaigns", []))))
    suites = [t["id"] for t in data["tests"] if source.campaign in t.get("campaigns", t.get("tiers", []))]
    if source.match:
        suites = [s for s in suites if any(fnmatch.fnmatchcase(s, g) for g in source.match)]
    return suites


# ----------------------------------------------------------------------------- plan
def apply_rules(items: List[Item], source: Source, machine: Machine, config: Config, routine: Routine,
                today: dt.date) -> None:
    """Mark as SKIP, with the reason, the tests this machine cannot run or runners.json switches off.

    A capability is required by the runner (runners.json "requires") or by the test itself: a Makefile.test
    variant such as <test>--valgrind needs the capability of its variant (table.VARIANTS)."""
    values = machine_vars(machine)
    missing_vars = [v for v in source.requires_vars if not values.get(v)]
    for item in items:
        if item.action != RUN:
            continue
        missing = [c for c in list(source.requires) + required_capabilities(item.make_name) if c not in machine.capabilities]
        entry = routine.disabled_entry(source.id, item.test, machine.machine, config.config, today)
        if missing:
            item.action, item.reason = SKIP, "{} does not declare capability {}".format(machine.machine, ", ".join(missing))
        elif missing_vars:
            item.action, item.reason = SKIP, "{} does not define {} in its machine file".format(
                machine.machine, ", ".join(missing_vars))
        elif entry:
            item.action, item.reason = SKIP, entry.describe()


def plan(machine: Machine, config: Config, routine: Routine, registry: Registry, today: Optional[dt.date] = None,
         amps_dir: Optional[Path] = None, table: Optional[Path] = None,
         manifest: Optional[Path] = None) -> List[SourcePlan]:
    today = today or dt.date.today()
    plans: List[SourcePlan] = []
    for source in routine.sources:
        if not source.enabled:
            plans.append(SourcePlan(source, False, "switched off in runners.json" +
                                    (" ({})".format(source.note) if source.note else "")))
            continue
        if config.sources is not None and source.id not in config.sources:
            plans.append(SourcePlan(source, False, "not in the 'sources' of {} {}".format(machine.machine, config.config)))
            continue
        try:
            if source.kind == "makefile":
                items, note = legacy_items(machine, amps_dir, table)
            elif source.campaign:
                items = [Item(source.family or registry.family_for_suite(s), s, RUN)
                         for s in core_suites(source, manifest)]
                note = "src/test campaign {!r}{}".format(
                    source.campaign, ", suites matching " + " ".join(source.match) if source.match else "")
            else:
                items = [Item(source.family or UNASSIGNED, source.id, RUN, listed=False)]
                note = "tests are known once the runner has run; reader {}".format(source.results.get("reader", ""))
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
            plans.append(SourcePlan(source, False, "cannot list tests: {}".format(error)))
            continue

        apply_rules(items, source, machine, config, routine, today)
        plans.append(SourcePlan(source, True, note, items))
    return plans


def duplicates(plans: List[SourcePlan]) -> List[str]:
    """(family, test) pairs produced by more than one source; a run record cannot hold both."""
    owner: Dict[Tuple[str, str], str] = {}
    found = []
    for sp in plans:
        for item in sp.items:
            key = (item.family, item.test)
            if key in owner and owner[key] != sp.source.id:
                found.append("{} / {} comes from both {} and {}".format(item.family, item.test, owner[key], sp.source.id))
            owner.setdefault(key, sp.source.id)
    return found

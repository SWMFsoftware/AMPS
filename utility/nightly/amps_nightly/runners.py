"""The test procedure: utility/nightly/runners.json, the list of test-runner commands the orchestrator executes.

"checkout"  the repository and branch tested each night, and companion repositories placed inside it
"runners"   one entry per test runner (an existing test program), run in this order:
              kind "command"   prepare commands, then the runner command, in its own workspace; its report
                               is read by the named reader (amps_nightly/readers.py)
              kind "makefile"  the Makefile.test tests, run by the orchestrator itself (makefile_runner.py);
                               "manifest_target" names the make target that lists them (test_manifest, or
                               test_long_manifest for the Tier=long tests of Makefile.test.long)
"disabled"  tests switched off without being deleted: reported as SKIPPED with the reason until "until".
Commands are argv lists, or strings run by /bin/bash -c. {placeholders} are filled from the machine file
and the run (see ORCHESTRATOR.md). Keys starting with "_" are comments and are ignored.

Routines: runners.json is the routine "nightly", run every night. Another routine <name> (e.g. "long", run
periodically: `amps-nightly run --routine long`) is the file runners.<name>.json next to it, with the same
schema; a machine runs it only when its machine file enables it ("routines": {"<name>": {...}}).
"""

from __future__ import annotations

import datetime as dt
import fnmatch
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .table import MANIFESTS

DEFAULT_PATH = Path(__file__).resolve().parent.parent / "runners.json"
NIGHTLY = "nightly"                 # the routine of runners.json
ROUTINE_NAME = re.compile(r"^[a-z][a-z0-9-]{0,30}$")
KINDS = ("makefile", "command")
WORKSPACES = ("configured", "plain")
DEFAULT_MANIFEST_TARGET = "test_manifest"
_ID = re.compile(r"^[a-z0-9][a-z0-9-]{0,62}$")
_RUNNER_KEYS = {"id", "kind", "title", "enabled", "note", "requires", "requires_vars", "workspace", "prepare",
                "command", "env", "time_limit_min", "jobs", "results", "families", "logs", "expect",
                "resources", "locks", "after", "manifest_target"}
_DISABLED_KEYS = {"test", "reason", "source", "machines", "configs", "until"}
Command = Union[str, List[str]]


class RoutineError(ValueError):
    """runners.json is malformed."""


@dataclass
class Companion:
    path: str                  # directory inside the checkout
    repository: str = ""       # cloned from a local mirror of this repository ...
    branch: str = "master"
    link_var: str = ""         # ... or a symbolic link to the machine variable of this name (shared data)


@dataclass
class Checkout:
    repository: str
    branch: str
    companions: List[Companion]


@dataclass
class Source:
    """One test runner."""
    id: str
    kind: str
    title: str
    enabled: bool
    note: str
    requires: List[str]                       # machine capabilities
    requires_vars: List[str] = field(default_factory=list)   # machine vars the commands use
    workspace: str = "configured"             # configured: Config.pl calls of the configuration; plain: none
    prepare: List[Command] = field(default_factory=list)
    command: Command = ""
    env: Dict[str, str] = field(default_factory=dict)
    time_limit_min: int = 120
    jobs: int = 0                             # upper limit for {jobs}; 0 = all cores of the machine
    results: Dict[str, str] = field(default_factory=dict)    # reader, file
    families: str = ""                        # "by-prefix" or a family id
    logs: List[str] = field(default_factory=list)            # globs of log files to archive
    campaign: str = ""                        # expect: src/test campaign whose suites this runner reports
    match: List[str] = field(default_factory=list)           # expect: suite-ID globs within the campaign
    manifest_target: str = DEFAULT_MANIFEST_TARGET           # kind makefile: the make target listing its tests

    @property
    def family(self) -> str:
        return "" if self.families == "by-prefix" else self.families


@dataclass
class Disabled:
    test: str
    reason: str
    source: str = ""
    machines: List[str] = field(default_factory=list)
    configs: List[str] = field(default_factory=list)
    until: Optional[dt.date] = None
    number: int = 0

    def expired(self, today: dt.date) -> bool:
        return self.until is not None and today > self.until

    def applies(self, source: str, test: str, machine: str, config: str, today: dt.date) -> bool:
        return (not self.expired(today) and fnmatch.fnmatchcase(test, self.test)
                and (not self.source or self.source == source)
                and (not self.machines or machine in self.machines)
                and (not self.configs or config in self.configs))

    def describe(self) -> str:
        return "disabled: {}{}".format(self.reason, " (until {})".format(self.until) if self.until else "")


# ----------------------------------------------------------------------------- field readers
def _strings(raw: Dict[str, Any], key: str, where: str) -> List[str]:
    value = raw.get(key, [])
    if not isinstance(value, list) or not all(isinstance(v, str) and v.strip() for v in value):
        raise RoutineError("{}: {} must be a list of non-empty strings".format(where, key))
    return [v.strip() for v in value]


def _text(raw: Dict[str, Any], key: str, where: str, required: bool = False, default: str = "") -> str:
    value = raw.get(key, default)
    if not isinstance(value, str) or (required and not value.strip()):
        raise RoutineError("{}: {} must be a {}string".format(where, key, "non-empty " if required else ""))
    return value.strip()


def _unknown(raw: Dict[str, Any], allowed: set, where: str) -> None:
    unknown = sorted(k for k in raw if k not in allowed and not k.startswith("_"))
    if unknown:
        raise RoutineError("{}: unknown key{} {}".format(where, "" if len(unknown) == 1 else "s", ", ".join(unknown)))


def _command(value: Any, where: str) -> Command:
    if isinstance(value, str) and value.strip():
        return value.strip()
    if isinstance(value, list) and value and all(isinstance(v, str) and v for v in value):
        return list(value)
    raise RoutineError(where + " must be a non-empty string or list of strings")


def _positive(raw: Dict[str, Any], key: str, where: str, default: int, minimum: int = 1) -> int:
    value = raw.get(key, default)
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise RoutineError("{}: {} must be an integer of at least {}".format(where, key, minimum))
    return value


def _source(raw: Any, number: int) -> Source:
    where = "runners[{}]".format(number)
    if not isinstance(raw, dict):
        raise RoutineError(where + ": must be an object")
    sid = raw.get("id", "")
    if not isinstance(sid, str) or not _ID.match(sid):
        raise RoutineError("{}: id must be lower case letters, digits and -".format(where))
    where = "runner {}".format(sid)
    _unknown(raw, _RUNNER_KEYS, where)
    kind = raw.get("kind", "command")
    if kind not in KINDS:
        raise RoutineError("{}: kind must be one of {}".format(where, ", ".join(KINDS)))
    enabled = raw.get("enabled", True)
    if not isinstance(enabled, bool):
        raise RoutineError(where + ": enabled must be true or false")
    workspace = raw.get("workspace", "configured")
    if workspace not in WORKSPACES:
        raise RoutineError("{}: workspace must be one of {}".format(where, ", ".join(WORKSPACES)))
    env = raw.get("env", {})
    if not isinstance(env, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in env.items()):
        raise RoutineError(where + ": env must map names to strings")
    prepare = raw.get("prepare", [])
    if not isinstance(prepare, list):
        raise RoutineError(where + ": prepare must be a list of commands")
    expect = raw.get("expect", {})
    if not isinstance(expect, dict) or set(expect) - {"src_test_campaign", "match"}:
        raise RoutineError(where + ": expect may only set src_test_campaign and match")
    source = Source(
        id=sid, kind=kind, title=_text(raw, "title", where) or sid, enabled=enabled, note=_text(raw, "note", where),
        requires=_strings(raw, "requires", where), requires_vars=_strings(raw, "requires_vars", where),
        workspace=workspace, prepare=[_command(p, "{} prepare[{}]".format(where, i + 1)) for i, p in enumerate(prepare)],
        env=dict(env), time_limit_min=_positive(raw, "time_limit_min", where, 120),
        jobs=_positive(raw, "jobs", where, 0, minimum=0), families=_text(raw, "families", where),
        logs=_strings(raw, "logs", where), campaign=_text(expect, "src_test_campaign", where),
        match=_strings(expect, "match", where))
    if "manifest_target" in raw:
        target = raw["manifest_target"]
        if kind != "makefile":
            raise RoutineError(where + ": manifest_target is only for kind makefile")
        if target not in MANIFESTS:
            raise RoutineError("{}: manifest_target must be one of {} (the targets table_parser writes)".format(
                where, ", ".join(sorted(MANIFESTS))))
        source.manifest_target = target
    if kind == "command":
        source.command = _command(raw.get("command"), where + ": command")
        results = raw.get("results")
        if not isinstance(results, dict) or not isinstance(results.get("reader"), str):
            raise RoutineError(where + ": results must name a reader")
        source.results = {k: str(v) for k, v in results.items()}
        if not source.families:
            raise RoutineError(where + ": families must be 'by-prefix' or a family id")
    return source


def _disabled(raw: Any, number: int) -> Disabled:
    where = "disabled[{}]".format(number)
    if not isinstance(raw, dict):
        raise RoutineError(where + ": must be an object")
    _unknown(raw, _DISABLED_KEYS, where)
    until = None
    if raw.get("until") is not None:
        try:
            until = dt.datetime.strptime(str(raw["until"]), "%Y-%m-%d").date()
        except ValueError:
            raise RoutineError(where + ": until must be a date YYYY-MM-DD")
    return Disabled(test=_text(raw, "test", where, required=True), reason=_text(raw, "reason", where, required=True),
                    source=_text(raw, "source", where), machines=_strings(raw, "machines", where),
                    configs=_strings(raw, "configs", where), until=until, number=number)


def _checkout(raw: Any) -> Checkout:
    if raw is None:
        return Checkout("", "master", [])
    if not isinstance(raw, dict):
        raise RoutineError("checkout must be an object")
    _unknown(raw, {"repository", "branch", "companions"}, "checkout")
    companions = []
    for number, item in enumerate(raw.get("companions", []), 1):
        where = "checkout.companions[{}]".format(number)
        if not isinstance(item, dict):
            raise RoutineError(where + ": must be an object")
        _unknown(item, {"repository", "path", "branch", "link_var"}, where)
        path = _text(item, "path", where, required=True)
        if Path(path).is_absolute() or ".." in Path(path).parts:
            raise RoutineError(where + ": path must be inside the checkout")
        companion = Companion(path=path, repository=_text(item, "repository", where),
                              branch=_text(item, "branch", where, default="master") or "master",
                              link_var=_text(item, "link_var", where))
        if bool(companion.repository) == bool(companion.link_var):
            raise RoutineError(where + ": give exactly one of repository and link_var")
        companions.append(companion)
    return Checkout(_text(raw, "repository", "checkout"), _text(raw, "branch", "checkout", default="master"),
                    companions)


def file_name(routine: str) -> str:
    """runners.json for the nightly routine, runners.<routine>.json for another."""
    return "runners.json" if routine == NIGHTLY else "runners.{}.json".format(routine)


def routine_path(routine: str, nightly: Optional[Path] = None) -> Path:
    """The file of a routine: the nightly's runners.json (default: utility/nightly/runners.json), or
    runners.<routine>.json next to it."""
    nightly = Path(nightly) if nightly else DEFAULT_PATH
    return nightly if routine == NIGHTLY else nightly.with_name(file_name(routine))


def routine_names(nightly: Optional[Path] = None) -> List[str]:
    """nightly, then every routine that has a runners.<routine>.json next to runners.json, by name."""
    nightly = Path(nightly) if nightly else DEFAULT_PATH
    found = []
    for path in sorted(nightly.parent.glob("runners.*.json")) if nightly.parent.is_dir() else []:
        name = path.name[len("runners."):-len(".json")]
        if ROUTINE_NAME.match(name) and name != NIGHTLY:
            found.append(name)
    return [NIGHTLY] + found


class Routine:
    def __init__(self, data: Any, name: str = NIGHTLY):
        if not ROUTINE_NAME.match(name):
            raise RoutineError("routine name {!r} must be lower case letters, digits and -".format(name))
        self.name = name                    # the routine: nightly (runners.json), or <name> (runners.<name>.json)
        self.label = label = file_name(name)
        if not isinstance(data, dict) or data.get("schema_version") != 1:
            raise RoutineError("{} must be an object with schema_version 1".format(label))
        _unknown(data, {"schema_version", "checkout", "runners", "disabled"}, label)
        if not isinstance(data.get("runners"), list) or not isinstance(data.get("disabled", []), list):
            raise RoutineError("{}: runners and disabled must be lists".format(label))
        try:
            self.checkout = _checkout(data.get("checkout"))
            self.sources = [_source(raw, n) for n, raw in enumerate(data["runners"], 1)]
        except RoutineError as error:
            raise RoutineError(str(error) if name == NIGHTLY else "{}: {}".format(label, error))
        ids = [s.id for s in self.sources]
        duplicates = sorted({i for i in ids if ids.count(i) > 1})
        if duplicates:
            raise RoutineError("{}: runner id used twice: {}".format(label, ", ".join(duplicates)))
        try:
            self.disabled = [_disabled(raw, n) for n, raw in enumerate(data.get("disabled", []), 1)]
            for entry in self.disabled:
                if entry.source and entry.source not in ids:
                    raise RoutineError("disabled[{}]: unknown runner {!r}".format(entry.number, entry.source))
        except RoutineError as error:
            raise RoutineError(str(error) if name == NIGHTLY else "{}: {}".format(label, error))

    @classmethod
    def load(cls, path: Optional[Path] = None, name: str = NIGHTLY) -> "Routine":
        path = Path(path) if path else DEFAULT_PATH
        try:
            return cls(json.loads(path.read_text()), name)
        except RoutineError:
            raise
        except (OSError, ValueError) as error:
            raise RoutineError("cannot read {}: {}".format(path, error))

    @classmethod
    def load_routine(cls, name: str, nightly: Optional[Path] = None) -> "Routine":
        """The routine `name`, from runners.json (nightly) or the runners.<name>.json next to it."""
        path = routine_path(name, nightly)
        if name != NIGHTLY and not path.is_file():
            raise RoutineError("no routine {!r}: {} does not exist (routines: {})".format(
                name, path, ", ".join(routine_names(nightly))))
        return cls.load(path, name)

    def source(self, sid: str) -> Optional[Source]:
        return next((s for s in self.sources if s.id == sid), None)

    def disabled_entry(self, source: str, test: str, machine: str, config: str,
                       today: dt.date) -> Optional[Disabled]:
        return next((d for d in self.disabled if d.applies(source, test, machine, config, today)), None)

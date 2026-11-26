"""Machine files: utility/nightly/machines/<machine>.json, one file per test machine.

Each file names the machine, how to build AMPS there (one entry per compiler configuration; each
configuration is one column of the report), where to work and where to upload the results.
Adding a machine is adding a file; "enabled": false pauses it and hides its columns on the report
(its history stays in the results store); deleting the file removes it.
Keys starting with "_" are comments and are ignored.

"routines" enables routines other than the nightly on the machine (each one is runners.<routine>.json):
    "routines": {"long": {"enabled": true, "start": "Sat 10:00", "time_limit_h": 34, "lanes": 1, ...}}
While routine <routine> runs, its values replace the machine's (start, night_limit_h by time_limit_h, lanes,
compile_lanes) or are merged into them (timeouts_s, make_vars). A routine the file does not enable does not run
on the machine.
"""

from __future__ import annotations

import dataclasses
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .records import column_id

DEFAULT_DIR = Path(__file__).resolve().parent.parent / "machines"
REPO_ROOT = Path(__file__).resolve().parents[3]
STAGES = ("compile", "rundir", "run", "check")
DEFAULT_TIMEOUTS = {"compile": 3600, "rundir": 600, "run": 7200, "check": 900}
_NAME = re.compile(r"^[a-z0-9][a-z0-9_.+-]{0,63}$")
_TIME = re.compile(r"^([01][0-9]|2[0-3]):[0-5][0-9]$")
_MACHINE_KEYS = {"schema_version", "machine", "enabled", "description", "hostnames", "table_host", "workdir",
                 "upload", "start", "core_profile", "capabilities", "lanes", "timeouts_s", "configs",
                 "python", "search_path", "mpiexec", "vars", "make_vars", "cores", "memory_gb", "max_parallel",
                 "night_limit_h", "keep_workspaces", "server_root", "min_free_gb", "compile_lanes", "compile_jobs",
                 "routines"}
KEEP_WORKSPACES = ("failures", "always", "never")
_CONFIG_KEYS = {"config", "enabled", "compiler", "toolchain", "setup", "config_pl", "after_install", "sources",
                "compile_jobs"}
# The make variable of MakefileTest/Table that sets the parallel jobs (make -j) of a test's compile stage. The nightly
# sets it from compile_jobs (see Machine.compile_jobs_for); make_vars may not set it.
COMPILE_JOBS_VAR = "TESTCOMPILEJOBS"
NIGHTLY = "nightly"
_ROUTINE_NAME = re.compile(r"^[a-z][a-z0-9-]{0,30}$")
_ROUTINE_KEYS = {"enabled", "start", "time_limit_h", "lanes", "compile_lanes", "timeouts_s", "make_vars", "max_ranks",
                 "stale_after_h"}
DEFAULT_STALE_AFTER_H = 192          # the report flags a routine's results older than this (nightly columns: 36 h)
DAYS = ("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")       # cron's day-of-week numbers 0-6
_START = re.compile(r"^(?:(Sun|Mon|Tue|Wed|Thu|Fri|Sat) )?([01]?[0-9]|2[0-3]):([0-5][0-9])$")


def parse_start(start: str) -> Tuple[Optional[int], int, int]:
    """(cron day of week or None for every day, hour, minute) of "HH:MM" or "<Day> HH:MM" ("Sat 10:00")."""
    match = _START.match(start.strip())
    if not match:
        raise ValueError("start {!r} must be HH:MM or <day> HH:MM, e.g. 'Sat 10:00' (days: {})".format(
            start, ", ".join(DAYS)))
    day = DAYS.index(match.group(1)) if match.group(1) else None
    return day, int(match.group(2)), int(match.group(3))


class MachineError(ValueError):
    """A machine file is malformed."""


@dataclass
class Config:
    """One compiler configuration of a machine = one column of the report."""
    machine: str
    config: str
    compiler: str
    toolchain: str
    enabled: bool
    setup: List[str]            # shell lines run before building (module loads, environment)
    config_pl: List[str]        # Config.pl invocations, in order; the first one installs
    after_install: List[str]    # shell lines run after the Config.pl invocations
    sources: Optional[List[str]]  # runners.json source ids to run; None = every enabled source
    compile_jobs: int = 0       # make -j of a Makefile.test compile with this compiler; 0 = the machine's

    @property
    def column(self) -> str:
        return column_id(self.machine, self.config)

    @property
    def folder(self) -> str:
        """Directory of this configuration's workspaces in the workdir: the compiler name (GNU, Intel, ...),
        as the per-compiler folders of the old routine's ~/Tmp_AMPS_test."""
        return re.sub(r"[^A-Za-z0-9._+-]", "_", self.compiler or self.config)


@dataclass
class RoutineSettings:
    """A routine other than the nightly on one machine (machine file "routines")."""
    name: str
    start: str                                  # "<Day> HH:MM" (weekly) or "HH:MM" (daily), local time
    enabled: bool = True
    time_limit_h: float = 0                     # replaces night_limit_h; 0 = the machine's
    lanes: int = 0                              # replaces lanes; 0 = the machine's
    compile_lanes: int = 0                      # replaces compile_lanes; 0 = the machine's
    timeouts_s: Dict[str, int] = field(default_factory=dict)    # merged into the machine's timeouts_s
    make_vars: Dict[str, str] = field(default_factory=dict)     # merged into the machine's make_vars
    max_ranks: int = 0                          # most MPI ranks of a test (Ranks= of Tier=long); 0 = cores
    stale_after_h: float = DEFAULT_STALE_AFTER_H   # the report flags results older than this


@dataclass
class Machine:
    machine: str
    path: Path
    enabled: bool
    description: str
    hostnames: List[str]
    table_host: str             # host name used by the @@(...) selectors of MakefileTest/Table
    workdir: str
    upload: str                 # local directory, or [user@]host:directory (rsync over ssh)
    start: str                  # local start time HH:MM of the nightly run
    core_profile: str           # src/test machine profile (path relative to the repository), "" = none
    capabilities: List[str]
    lanes: int
    timeouts_s: Dict[str, int]
    configs: List[Config]
    python: str = "/usr/bin/python3"
    search_path: List[str] = None     # PATH of the runners (default /usr/local/bin, /usr/bin, /bin)
    mpiexec: str = "mpiexec"
    vars: Dict[str, str] = None       # machine-specific values for {placeholders} in runners.json
    make_vars: Dict[str, str] = None  # variables passed to the Makefile.test targets (MPIRUN=...)
    cores: int = 1
    memory_gb: int = 0
    max_parallel: int = 1             # runners at once; 1 = sequential (phase 1)
    night_limit_h: float = 14         # runners not started within this many hours are reported as skipped
    keep_workspaces: str = "failures"
    server_root: str = ""             # results server only: results.db, incoming/, www/ (amps-nightly update)
    min_free_gb: float = 5            # no runner or Makefile.test test starts with less free space in workdir
    compile_lanes: int = 1            # Makefile.test tests compiled at once, each in its own copy of the workspace
    compile_jobs: int = 0             # make -j of each Makefile.test compile; 0 = cores / (compile_lanes x max_parallel)
    routines: Dict[str, RoutineSettings] = None   # routines other than the nightly ("routines")
    routine: str = NIGHTLY            # the routine this Machine's values are for (for_routine)

    def runs_routine(self, name: str) -> bool:
        """The machine file enables this routine (the nightly always)."""
        settings = (self.routines or {}).get(name)
        return name == NIGHTLY or bool(settings and settings.enabled)

    def enabled_routines(self) -> List[str]:
        """The routines other than the nightly that run here, by name."""
        return sorted(name for name, s in (self.routines or {}).items() if s.enabled)

    def max_ranks(self, name: str) -> int:
        """The most MPI ranks a test of routine `name` may use here: its max_ranks, else the machine's cores."""
        settings = (self.routines or {}).get(name)
        return settings.max_ranks if settings and settings.max_ranks else self.cores

    def stale_after_h(self, name: str) -> float:
        settings = (self.routines or {}).get(name)
        return settings.stale_after_h if settings else DEFAULT_STALE_AFTER_H

    def for_routine(self, name: str) -> "Machine":
        """This machine as routine `name` sees it: the routine's values in place of the machine's."""
        if name == NIGHTLY:
            return self
        settings = (self.routines or {}).get(name)
        if settings is None or not settings.enabled:
            raise MachineError("{}: routine {} is not enabled (\"routines\" of the machine file)".format(
                self.path.name, name))
        timeouts = dict(self.timeouts_s)
        timeouts.update(settings.timeouts_s)
        make_vars = dict(self.make_vars or {})
        make_vars.update(settings.make_vars)
        return dataclasses.replace(
            self, routine=name, start=settings.start, night_limit_h=settings.time_limit_h or self.night_limit_h,
            lanes=settings.lanes or self.lanes, compile_lanes=settings.compile_lanes or self.compile_lanes,
            timeouts_s=timeouts, make_vars=make_vars)

    def active_configs(self) -> List[Config]:
        return [c for c in self.configs if c.enabled] if self.enabled else []

    def config(self, name: str) -> Optional[Config]:
        return next((c for c in self.configs if c.config == name), None)

    def compile_jobs_for(self, config: Optional[Config] = None) -> int:
        """make -j of one Makefile.test compile (TESTCOMPILEJOBS): the configuration's compile_jobs, else the
        machine's, else the cores shared by the compile lanes (and by the runners that may run at once), so the
        lanes together use the whole machine."""
        if config is not None and config.compile_jobs:
            return config.compile_jobs
        if self.compile_jobs:
            return self.compile_jobs
        return max(1, self.cores // (self.compile_lanes * self.max_parallel))


# ----------------------------------------------------------------------------- field readers
def _keys(raw: Dict[str, Any], allowed: set, where: str) -> None:
    unknown = sorted(k for k in raw if k not in allowed and not k.startswith("_"))
    if unknown:
        raise MachineError("{}: unknown key{} {}".format(where, "" if len(unknown) == 1 else "s", ", ".join(unknown)))


def _str(raw: Dict[str, Any], key: str, where: str, default: Optional[str] = None) -> str:
    value = raw.get(key, default)
    if not isinstance(value, str) or (default is None and not value.strip()):
        raise MachineError("{}: {} must be {}a string".format(where, key, "" if default is not None else "a non-empty "))
    return value.strip()


def _strings(raw: Dict[str, Any], key: str, where: str, default: Optional[List[str]] = None) -> List[str]:
    value = raw.get(key, default)
    if not isinstance(value, list) or not all(isinstance(v, str) and v.strip() for v in value):
        raise MachineError("{}: {} must be a list of non-empty strings".format(where, key))
    return [v.strip() for v in value]


def _bool(raw: Dict[str, Any], key: str, where: str, default: bool) -> bool:
    value = raw.get(key, default)
    if not isinstance(value, bool):
        raise MachineError("{}: {} must be true or false".format(where, key))
    return value


def _hours(value: Any, where: str, key: str = "night_limit_h", limit: float = 48) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 < value <= limit:
        raise MachineError("{}: {} must be a number of hours between 0 and {:g}".format(where, key, limit))
    return float(value)


def _choice(raw: Dict[str, Any], key: str, choices: tuple, where: str) -> str:
    value = raw.get(key, choices[0])
    if value not in choices:
        raise MachineError("{}: {} must be one of {}".format(where, key, ", ".join(choices)))
    return value


def _positive(value: Any, label: str, where: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise MachineError("{}: {} must be a positive integer".format(where, label))
    return value


# ----------------------------------------------------------------------------- loading
def parse_machine(raw: Any, path: Path) -> Machine:
    where = path.name
    if not isinstance(raw, dict) or raw.get("schema_version") != 1:
        raise MachineError("{}: must be a JSON object with schema_version 1".format(where))
    _keys(raw, _MACHINE_KEYS, where)
    name = _str(raw, "machine", where)
    if not _NAME.match(name):
        raise MachineError("{}: machine {!r} must be lower case letters, digits and . _ + -".format(where, name))
    if path.suffix == ".json" and name != path.stem:
        raise MachineError("{}: machine is {!r} but the file is named {}".format(where, name, path.name))
    start = _str(raw, "start", where, "22:30")
    if not _TIME.match(start):
        raise MachineError("{}: start must be HH:MM (24-hour local time)".format(where))
    core = raw.get("core_profile", "")
    if core is None:
        core = ""
    if not isinstance(core, str):
        raise MachineError("{}: core_profile must be a path or null".format(where))
    timeouts = dict(DEFAULT_TIMEOUTS)
    given = raw.get("timeouts_s", {})
    if not isinstance(given, dict) or set(given) - set(STAGES):
        raise MachineError("{}: timeouts_s may only set {}".format(where, ", ".join(STAGES)))
    for stage, seconds in given.items():
        timeouts[stage] = _positive(seconds, "timeouts_s." + stage, where)

    configs: List[Config] = []
    raw_configs = raw.get("configs")
    if not isinstance(raw_configs, list) or not raw_configs:
        raise MachineError("{}: configs must be a non-empty list".format(where))
    for number, item in enumerate(raw_configs, 1):
        at = "{} configs[{}]".format(where, number)
        if not isinstance(item, dict):
            raise MachineError(at + ": must be an object")
        _keys(item, _CONFIG_KEYS, at)
        config = _str(item, "config", at)
        if not _NAME.match(config):
            raise MachineError("{}: config {!r} must be lower case letters, digits and . _ + -".format(at, config))
        if any(c.config == config for c in configs):
            raise MachineError("{}: config {!r} appears twice".format(at, config))
        sources = item.get("sources")
        configs.append(Config(
            machine=name, config=config, compiler=_str(item, "compiler", at, config.upper()) or config.upper(),
            toolchain=_str(item, "toolchain", at, ""), enabled=_bool(item, "enabled", at, True),
            setup=_strings(item, "setup", at, []), config_pl=_strings(item, "config_pl", at),
            after_install=_strings(item, "after_install", at, []),
            sources=None if sources is None else _strings(item, "sources", at),
            compile_jobs=_positive(item["compile_jobs"], "compile_jobs", at) if "compile_jobs" in item else 0))

    folders = [c.folder for c in configs]
    for folder in set(folders):
        if folders.count(folder) > 1:
            raise MachineError("{}: configs {} have the same compiler {!r}, so their workspaces would share "
                               "{}/{}; give them different \"compiler\" names".format(
                                   where, ", ".join(c.config for c in configs if c.folder == folder), folder,
                                   "<workdir>", folder))
    make_vars = _mapping(raw, "make_vars", where)
    if COMPILE_JOBS_VAR in make_vars:
        raise MachineError("{}: make_vars may not set {}: compile_jobs (of the machine or of a configuration) "
                           "sets it".format(where, COMPILE_JOBS_VAR))
    routines = _routines(raw, where)
    return Machine(
        machine=name, path=path, enabled=_bool(raw, "enabled", where, True),
        description=_str(raw, "description", where, ""),
        hostnames=_strings(raw, "hostnames", where, [name]),
        table_host=_str(raw, "table_host", where, name) or name,
        workdir=_str(raw, "workdir", where, "~/AMPSTEST/work") or "~/AMPSTEST/work",
        upload=_str(raw, "upload", where), start=start, core_profile=core.strip(),
        capabilities=_strings(raw, "capabilities", where, []),
        lanes=_positive(raw.get("lanes", 1), "lanes", where), timeouts_s=timeouts, configs=configs,
        python=_str(raw, "python", where, "/usr/bin/python3") or "/usr/bin/python3",
        search_path=_strings(raw, "search_path", where, ["/usr/local/bin", "/usr/bin", "/bin"]),
        mpiexec=_str(raw, "mpiexec", where, "mpiexec") or "mpiexec",
        vars=_mapping(raw, "vars", where), make_vars=make_vars,
        cores=_positive(raw.get("cores", 1), "cores", where),
        memory_gb=0 if raw.get("memory_gb", 0) == 0 else _positive(raw.get("memory_gb"), "memory_gb", where),
        max_parallel=_positive(raw.get("max_parallel", 1), "max_parallel", where),
        night_limit_h=_hours(raw.get("night_limit_h", 14), where),
        keep_workspaces=_choice(raw, "keep_workspaces", KEEP_WORKSPACES, where),
        server_root=_str(raw, "server_root", where, ""), min_free_gb=_free(raw, where),
        compile_lanes=_positive(raw.get("compile_lanes", 1), "compile_lanes", where),
        compile_jobs=_positive(raw["compile_jobs"], "compile_jobs", where) if "compile_jobs" in raw else 0,
        routines=routines)


def _routines(raw: Dict[str, Any], where: str) -> Dict[str, RoutineSettings]:
    """The "routines" object: routine name -> its settings on this machine."""
    given = raw.get("routines", {})
    if not isinstance(given, dict):
        raise MachineError("{}: routines must map routine names to objects".format(where))
    routines: Dict[str, RoutineSettings] = {}
    for name, item in given.items():
        if name.startswith("_"):
            continue
        at = "{} routines.{}".format(where, name)
        if name == NIGHTLY:
            raise MachineError("{}: the nightly routine is set by the machine's own keys (start, night_limit_h, "
                               "lanes, timeouts_s, ...), not in routines".format(at))
        if not _ROUTINE_NAME.match(name):
            raise MachineError("{}: a routine name must be lower case letters, digits and - (its file is "
                               "runners.<name>.json)".format(at))
        if not isinstance(item, dict):
            raise MachineError(at + ": must be an object")
        _keys(item, _ROUTINE_KEYS, at)
        start = _str(item, "start", at)
        try:
            parse_start(start)
        except ValueError as error:
            raise MachineError("{}: {}".format(at, error))
        timeouts = item.get("timeouts_s", {})
        if not isinstance(timeouts, dict) or set(timeouts) - set(STAGES):
            raise MachineError("{}: timeouts_s may only set {}".format(at, ", ".join(STAGES)))
        make_vars = _mapping(item, "make_vars", at)
        if COMPILE_JOBS_VAR in make_vars:
            raise MachineError("{}: make_vars may not set {}: compile_jobs of the machine sets it".format(
                at, COMPILE_JOBS_VAR))
        routines[name] = RoutineSettings(
            name=name, start=start, enabled=_bool(item, "enabled", at, True),
            time_limit_h=_hours(item["time_limit_h"], at, "time_limit_h", 168) if "time_limit_h" in item else 0,
            lanes=_positive(item["lanes"], "lanes", at) if "lanes" in item else 0,
            compile_lanes=_positive(item["compile_lanes"], "compile_lanes", at) if "compile_lanes" in item else 0,
            timeouts_s={stage: _positive(seconds, "timeouts_s." + stage, at) for stage, seconds in timeouts.items()},
            make_vars=make_vars,
            max_ranks=_positive(item["max_ranks"], "max_ranks", at) if "max_ranks" in item else 0,
            stale_after_h=_hours(item["stale_after_h"], at, "stale_after_h", 24 * 366)
            if "stale_after_h" in item else DEFAULT_STALE_AFTER_H)
    return routines


def _free(raw: Dict[str, Any], where: str) -> float:
    value = raw.get("min_free_gb", 5)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
        raise MachineError("{}: min_free_gb must be a number of GB, 0 or more".format(where))
    return float(value)


def _mapping(raw: Dict[str, Any], key: str, where: str) -> Dict[str, str]:
    value = raw.get(key, {})
    if not isinstance(value, dict) or not all(isinstance(k, str) and k and isinstance(v, str) for k, v in value.items()):
        raise MachineError("{}: {} must map names to strings".format(where, key))
    return dict(value)


def load_machine(path: Path) -> Machine:
    path = Path(path)
    try:
        raw = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise MachineError("{}: {}".format(path.name, error))
    return parse_machine(raw, path)


def load_machines(directory: Optional[Path] = None) -> Tuple[List[Machine], List[str]]:
    """Every machines/*.json, sorted by machine name, and the problems of files that could not be read."""
    directory = Path(directory) if directory else DEFAULT_DIR
    machines, problems = [], []
    for path in sorted(directory.glob("*.json")) if directory.is_dir() else []:
        try:
            machines.append(load_machine(path))
        except MachineError as error:
            problems.append(str(error))
    return machines, problems


def find_machine(machines: List[Machine], name: str = "", hostname: str = "") -> Optional[Machine]:
    """The machine called name, or else the one listing hostname (full or short) in its hostnames."""
    if name:
        return next((m for m in machines if m.machine == name), None)
    short = hostname.split(".", 1)[0]
    return next((m for m in machines if hostname in m.hostnames or short in m.hostnames), None)


def declared_columns(machines: List[Machine]) -> List[Dict[str, str]]:
    """Columns the report shows: every enabled configuration of every enabled machine, in file order."""
    return [{"column": c.column, "machine": m.machine, "config": c.config, "compiler": c.compiler,
             "toolchain": c.toolchain} for m in machines for c in m.active_configs()]


def declared_routines(machines: List[Machine]) -> Dict[str, List[Dict[str, Any]]]:
    """The report's columns of each routine other than the nightly: the declared columns of the machines that enable
    it, each with the routine's stale_after_h there."""
    names = sorted({name for m in machines if m.enabled for name in m.enabled_routines()})
    return {name: [{"column": c.column, "machine": m.machine, "config": c.config, "compiler": c.compiler,
                    "toolchain": c.toolchain, "stale_after_h": m.stale_after_h(name)}
                   for m in machines if m.runs_routine(name) for c in m.active_configs()]
            for name in names}

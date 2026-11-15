"""Machine files: utility/nightly/machines/<machine>.json, one file per test machine.

Each file names the machine, how to build AMPS there (one entry per compiler configuration; each
configuration is one column of the report), where to work and where to upload the results.
Adding a machine is adding a file; "enabled": false pauses it and hides its columns on the report
(its history stays in the results store); deleting the file removes it.
Keys starting with "_" are comments and are ignored.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
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
                 "night_limit_h", "keep_workspaces"}
KEEP_WORKSPACES = ("failures", "always", "never")
_CONFIG_KEYS = {"config", "enabled", "compiler", "toolchain", "setup", "config_pl", "after_install", "sources"}


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

    @property
    def column(self) -> str:
        return column_id(self.machine, self.config)


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

    def active_configs(self) -> List[Config]:
        return [c for c in self.configs if c.enabled] if self.enabled else []

    def config(self, name: str) -> Optional[Config]:
        return next((c for c in self.configs if c.config == name), None)


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


def _hours(value: Any, where: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 < value <= 48:
        raise MachineError("{}: night_limit_h must be a number of hours between 0 and 48".format(where))
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
            sources=None if sources is None else _strings(item, "sources", at)))

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
        vars=_mapping(raw, "vars", where), make_vars=_mapping(raw, "make_vars", where),
        cores=_positive(raw.get("cores", 1), "cores", where),
        memory_gb=0 if raw.get("memory_gb", 0) == 0 else _positive(raw.get("memory_gb"), "memory_gb", where),
        max_parallel=_positive(raw.get("max_parallel", 1), "max_parallel", where),
        night_limit_h=_hours(raw.get("night_limit_h", 14), where),
        keep_workspaces=_choice(raw, "keep_workspaces", KEEP_WORKSPACES, where))


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

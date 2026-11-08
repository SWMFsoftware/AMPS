"""Portable, layered configuration for the AMPS core-test framework.

The framework deliberately keeps machine paths out of test implementations.
Committed defaults are overlaid by zero or more JSON files, environment
variables, and finally command-line overrides.  Relative paths always refer
to the AMPS repository root; environment variables and ``~`` are expanded.
"""

from __future__ import annotations

import copy
import json
import os
import shlex
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional


class ConfigurationError(ValueError):
    """A configuration file or override violates the public schema."""


_ENVIRONMENT_OVERRIDES = {
    "AMPS_TEST_OUTPUT_ROOT": ("paths", "output_root", "path"),
    "AMPS_TEST_TEMP_ROOT": ("paths", "temporary_root", "path"),
    "AMPS_TEST_DATA_ROOT": ("paths", "data_root", "path"),
    # SPICE installations are site resources and must never be embedded in an
    # application adapter.  These optional overrides are resolved with the
    # same rules as the ordinary test paths and are consumed only by adapters
    # that explicitly check their named machine capability before execution.
    "AMPS_TEST_SPICE_PATH": ("paths", "spice_path", "path"),
    "AMPS_TEST_SPICE_KERNELS": ("paths", "spice_kernels", "path"),
    "AMPS_TEST_PYTHON": ("tools", "python", "command"),
    "AMPS_TEST_CXX": ("tools", "cxx", "command"),
    "AMPS_TEST_MPI_CXX": ("tools", "mpi_cxx", "command"),
    "AMPS_TEST_MAKE": ("tools", "make", "command"),
    "AMPS_TEST_MPIEXEC": ("tools", "mpiexec", "command"),
    "AMPS_TEST_BUILD_JOBS": ("execution", "build_jobs", "integer"),
    "AMPS_TEST_MPI_RANKS": ("execution", "mpi_ranks", "integer"),
    "AMPS_TEST_MPI_RANK_FLAG": ("execution", "mpi_rank_flag", "text"),
    "AMPS_TEST_TIMEOUT": ("execution", "default_timeout_seconds", "integer"),
    "AMPS_TEST_CAMPAIGN": ("execution", "default_campaign", "text"),
}

_CAMPAIGNS = {
    "smoke", "nightly", "profile", "extended", "distributed",
    "application-smoke",
}


def _deep_merge(base: MutableMapping[str, Any], overlay: Mapping[str, Any]) -> None:
    """Recursively merge objects while replacing scalar/list leaves.

    Commands are lists and must be replaced as one argv vector.  Element-wise
    list merging would combine two compiler or launcher commands and produce a
    configuration that was never requested by either profile.
    """

    for key, value in overlay.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), MutableMapping):
            _deep_merge(base[key], value)
        else:
            base[key] = copy.deepcopy(value)


def _read_json(path: Path) -> Dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as stream:
            value = json.load(stream)
    except OSError as error:
        raise ConfigurationError("cannot read configuration {}: {}".format(path, error))
    except ValueError as error:
        raise ConfigurationError("invalid JSON in configuration {}: {}".format(path, error))
    if not isinstance(value, dict):
        raise ConfigurationError("configuration {} must contain a JSON object".format(path))
    if value.get("schema_version", 1) != 1:
        raise ConfigurationError("unsupported schema_version in {}".format(path))
    return value


def _command(value: Any, label: str) -> List[str]:
    """Normalize a tool into argv form; no runner later invokes a shell."""

    if isinstance(value, str):
        result = shlex.split(os.path.expandvars(os.path.expanduser(value)))
    elif isinstance(value, list) and all(isinstance(item, str) for item in value):
        result = [os.path.expandvars(os.path.expanduser(item)) for item in value]
    else:
        raise ConfigurationError("{} must be a command string or list of strings".format(label))
    if not result:
        raise ConfigurationError("{} must not be empty".format(label))
    return result


def _path(value: Any, repo_root: Path, label: str, allow_empty: bool = False) -> Optional[Path]:
    """Expand a portable path and anchor relative values at the repository root."""

    if value in (None, ""):
        if allow_empty:
            return None
        raise ConfigurationError("{} must not be empty".format(label))
    if not isinstance(value, str):
        raise ConfigurationError("{} must be a path string".format(label))
    expanded = Path(os.path.expandvars(os.path.expanduser(value)))
    return (expanded if expanded.is_absolute() else repo_root / expanded).resolve()


def _positive_integer(value: Any, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise ConfigurationError("{} must be an integer".format(label))
    if parsed <= 0:
        raise ConfigurationError("{} must be positive".format(label))
    return parsed


def _apply_environment(config: MutableMapping[str, Any], environ: Mapping[str, str]) -> None:
    """Apply the documented scalar environment layer with early type checks."""

    for variable, (section, key, kind) in _ENVIRONMENT_OVERRIDES.items():
        if variable not in environ:
            continue
        raw = environ[variable]
        if kind == "integer":
            value: Any = _positive_integer(raw, variable)
        elif kind == "command":
            value = _command(raw, variable)
        else:
            value = raw
        config.setdefault(section, {})[key] = value

    # Capability names are data rather than environment-variable suffixes:
    # sites may use names such as ``external-field`` that are not portable
    # shell identifiers. The available list describes the machine; the
    # required list upgrades an unavailable SKIP to a policy failure.
    available = environ.get("AMPS_TEST_CAPABILITIES")
    required = environ.get("AMPS_TEST_REQUIRED_CAPABILITIES")
    if available is not None or required is not None:
        capabilities = config.setdefault("capabilities", {})
        if not isinstance(capabilities, MutableMapping):
            raise ConfigurationError("capabilities must be an object")
        for name in [value.strip() for value in (available or "").split(",")
                     if value.strip()]:
            prior = capabilities.get(name, {})
            policy = prior.get("policy", "optional") if isinstance(prior, Mapping) else "optional"
            capabilities[name] = {"available": True, "policy": policy}
        for name in [value.strip() for value in (required or "").split(",")
                     if value.strip()]:
            prior = capabilities.get(name, {})
            if isinstance(prior, bool):
                is_available = prior
            elif isinstance(prior, Mapping):
                is_available = prior.get("available", False)
            else:
                is_available = False
            capabilities[name] = {"available": bool(is_available), "policy": "required"}


def _environment_config_paths(environ: Mapping[str, str]) -> List[Path]:
    raw = environ.get("AMPS_TEST_CONFIG", "").strip()
    if not raw:
        return []
    return [Path(item).expanduser().resolve() for item in raw.split(os.pathsep) if item]


def load_configuration(
    repo_root: Path,
    default_file: Path,
    config_files: Iterable[Path] = (),
    environ: Optional[Mapping[str, str]] = None,
    output_root: Optional[str] = None,
    build_jobs: Optional[int] = None,
    mpi_ranks: Optional[int] = None,
) -> Dict[str, Any]:
    """Load and validate the effective configuration.

    Precedence, low to high: committed defaults, files named by
    ``AMPS_TEST_CONFIG``, repeated ``--config`` files, environment scalar
    overrides, and explicit CLI overrides.
    """

    repo_root = repo_root.resolve()
    environment = dict(os.environ if environ is None else environ)
    config: Dict[str, Any] = _read_json(default_file)
    loaded_files = [str(default_file.resolve())]

    for path in _environment_config_paths(environment) + [Path(p).resolve() for p in config_files]:
        _deep_merge(config, _read_json(path))
        loaded_files.append(str(path))

    _apply_environment(config, environment)
    if output_root is not None:
        config.setdefault("paths", {})["output_root"] = output_root
    if build_jobs is not None:
        config.setdefault("execution", {})["build_jobs"] = build_jobs
    if mpi_ranks is not None:
        config.setdefault("execution", {})["mpi_ranks"] = mpi_ranks

    # Reject misspellings instead of silently ignoring them.  This matters on
    # clusters where a typo in a scratch or launcher key could redirect a long
    # test to the login-node filesystem or invoke the wrong MPI implementation.
    expected_sections = {"schema_version", "paths", "tools", "execution", "capabilities"}
    unknown = set(config) - expected_sections
    if unknown:
        raise ConfigurationError("unknown top-level configuration keys: {}".format(", ".join(sorted(unknown))))
    for section in ("paths", "tools", "execution"):
        if not isinstance(config.get(section), dict):
            raise ConfigurationError("configuration section '{}' must be an object".format(section))
    if not isinstance(config.get("capabilities", {}), dict):
        raise ConfigurationError("configuration section 'capabilities' must be an object")

    paths = config["paths"]
    tools = config["tools"]
    execution = config["execution"]
    raw_capabilities = config.get("capabilities", {})
    output = _path(paths.get("output_root"), repo_root, "paths.output_root")
    temporary = _path(paths.get("temporary_root"), repo_root,
                      "paths.temporary_root", allow_empty=True)
    data_root = _path(paths.get("data_root"), repo_root,
                      "paths.data_root", allow_empty=True)
    spice_path = _path(paths.get("spice_path"), repo_root,
                       "paths.spice_path", allow_empty=True)
    spice_kernels = _path(paths.get("spice_kernels"), repo_root,
                          "paths.spice_kernels", allow_empty=True)
    resolved_tools = {
        name: _command(value, "tools.{}".format(name)) for name, value in tools.items()
    }
    required_tools = {"python", "cxx", "mpi_cxx", "make", "mpiexec"}
    missing = required_tools - set(resolved_tools)
    if missing:
        raise ConfigurationError("missing tool definitions: {}".format(", ".join(sorted(missing))))

    keep_work = str(execution.get("keep_work", "failures")).lower()
    if keep_work not in ("always", "failures", "never"):
        raise ConfigurationError("execution.keep_work must be always, failures, or never")
    mpi_rank_flag = execution.get("mpi_rank_flag", "-n")
    if not isinstance(mpi_rank_flag, str) or not mpi_rank_flag.strip():
        raise ConfigurationError("execution.mpi_rank_flag must be a nonempty string")
    default_campaign = execution.get("default_campaign", "smoke")
    if default_campaign not in _CAMPAIGNS:
        raise ConfigurationError(
            "execution.default_campaign must be one of {}".format(
                ", ".join(sorted(_CAMPAIGNS))))

    capabilities: Dict[str, Dict[str, Any]] = {}
    for name, raw in raw_capabilities.items():
        if not isinstance(name, str) or not name:
            raise ConfigurationError("capability names must be nonempty strings")
        if isinstance(raw, bool):
            available = raw
            policy = "optional"
        elif isinstance(raw, Mapping):
            unknown_capability_keys = set(raw) - {"available", "policy"}
            if unknown_capability_keys:
                raise ConfigurationError(
                    "capabilities.{} has unknown keys: {}".format(
                        name, ", ".join(sorted(unknown_capability_keys))))
            available = raw.get("available", False)
            policy = raw.get("policy", "optional")
            if not isinstance(available, bool):
                raise ConfigurationError(
                    "capabilities.{}.available must be boolean".format(name))
        else:
            raise ConfigurationError(
                "capabilities.{} must be boolean or an object".format(name))
        if policy not in ("optional", "required"):
            raise ConfigurationError(
                "capabilities.{}.policy must be optional or required".format(name))
        capabilities[name] = {"available": bool(available), "policy": str(policy)}

    return {
        "schema_version": 1,
        "repo_root": str(repo_root),
        "loaded_files": loaded_files,
        "paths": {
            "output_root": str(output),
            "temporary_root": "" if temporary is None else str(temporary),
            "data_root": "" if data_root is None else str(data_root),
            "spice_path": "" if spice_path is None else str(spice_path),
            "spice_kernels": "" if spice_kernels is None else str(spice_kernels),
        },
        "tools": resolved_tools,
        "capabilities": capabilities,
        "execution": {
            "build_jobs": _positive_integer(execution.get("build_jobs", 16), "execution.build_jobs"),
            "mpi_ranks": _positive_integer(execution.get("mpi_ranks", 1), "execution.mpi_ranks"),
            "mpi_rank_flag": mpi_rank_flag,
            "default_timeout_seconds": _positive_integer(
                execution.get("default_timeout_seconds", 300),
                "execution.default_timeout_seconds"),
            "keep_work": keep_work,
            "default_campaign": str(default_campaign),
        },
    }


def command_available(command: List[str], repo_root: Path) -> bool:
    """Return whether the executable at the front of *command* is resolvable.

    Only argv[0] is probed: flags and wrapper arguments are not files.  A path
    containing a separator may be repository-relative; a bare executable is
    resolved through the current machine's PATH.
    """

    import shutil

    executable = os.path.expanduser(os.path.expandvars(command[0]))
    if os.sep in executable:
        path = Path(executable)
        if not path.is_absolute():
            path = repo_root / path
        return path.is_file() and os.access(str(path), os.X_OK)
    return shutil.which(executable) is not None

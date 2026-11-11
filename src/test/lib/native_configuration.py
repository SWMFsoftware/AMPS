"""Pin and verify a native fixture's global no-SPICE build selection.

An exosphere input-block ``SPICE=off`` controls component behavior; it does
not remove the global installation saved by Config.pl/Makefile.local. Native
fixtures with no target/frame/kernel contract must therefore select their
global mode explicitly. This helper runs only inside an adapter's existing
exclusive checkout lock, after its mandatory root/idle/build-removal preflight.
It never edits a site file or selects a mode for another application owner.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import shlex
from typing import Any, Dict, List, Mapping
import uuid

from .runner_support import RunnerError, command, run_command


class NativeConfigurationError(RunnerError):
    """A configuration failure with durable evidence before native execution."""

    def __init__(self, message: str, evidence: List[str]):
        super().__init__(message)
        self.evidence = list(evidence)


def _selection(path: Path) -> Dict[str, str]:
    # Config.pl writes simple key=value records. Preserve last-assignment
    # semantics rather than accepting an earlier, superseded SPICE=nospice.
    result: Dict[str, str] = {}
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines():
            match = re.fullmatch(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*", line)
            if match:
                result[match.group(1)] = match.group(2)
    return result


def _no_spice_define(options: str) -> bool:
    """Interpret -D/-U in compiler order, including split/combined spellings.

    These headers use #ifdef, so -D_NO_SPICE_CALLS_=0 still defines the macro;
    a later -U_NO_SPICE_CALLS_ removes it. Checking substring presence would
    accept contradictory compiler options and conceal the original defect.
    """
    tokens = iter(shlex.split(options))
    defined = False
    for token in tokens:
        if token in ("-D", "-U"):
            operation, value = token, next(tokens, "")
        elif token.startswith(("-D", "-U")):
            operation, value = token[:2], token[2:]
        else:
            continue
        if value.split("=", 1)[0] == "_NO_SPICE_CALLS_":
            defined = operation == "-D"
    return defined


@dataclass
class NativeConfiguration:
    command: List[str]
    stdout: str
    receipt_path: Path
    receipt: Dict[str, Any]
    evidence: List[str]

    def save(self) -> None:
        temporary = self.receipt_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(self.receipt, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
        temporary.replace(self.receipt_path)

    def reject(self, message: str) -> None:
        self.receipt.update(status="ERROR", message=message)
        self.save()
        raise NativeConfigurationError(
            message + "; configuration receipt=" + str(self.receipt_path), self.evidence)


def configure_no_spice(context: Mapping[str, Any], application: str,
                       input_deck: str, label: str) -> NativeConfiguration:
    """Use Config.pl, then inspect the actual make expansion without building.

    A private, dependency-free target is selected explicitly. GNU make -n amps
    is unsuitable here because recursive make recipes may still execute. The
    probe reads the same Makefile and tool/environment as the subsequent native
    build, while leaving the production Makefile and all site settings intact.
    """
    if not re.fullmatch(r"[A-Za-z0-9_-]+", label):
        raise RunnerError("invalid native configuration evidence label")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    artifacts.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    timeout = int(context["test"].get("timeout_seconds", 1800))
    argv = [str(repo / "Config.pl"), "-application=" + application,
            "-spice-path=nospice"]
    log = artifacts / (label + "-configuration.txt")
    receipt_path = artifacts / (label + "-configuration.json")
    result = NativeConfiguration(argv, "", receipt_path, {
        "schema_version": 1, "status": "STARTED", "configure_command": argv,
        "expected_application": application.rsplit("/", 1)[-1],
        "expected_input": Path(input_deck).name, "expected_spice": "nospice",
        "before_selection": _selection(repo / ".amps.conf"),
    }, [str(receipt_path), str(log)])
    result.save()
    log.write_text("", encoding="utf-8")
    try:
        configured = run_command(argv, repo, timeout=timeout)
        result.stdout = configured.stdout or ""
        log.write_text(result.stdout, encoding="utf-8")
        result.receipt["configure_returncode"] = configured.returncode
        if configured.returncode != 0:
            result.reject("Config.pl failed to select the native no-SPICE fixture")

        observed = _selection(repo / ".amps.conf")
        result.receipt["observed_selection"] = observed
        selected_input = repo / Path(input_deck).name
        checks = {
            "application": observed.get("APPLICATION") == result.receipt["expected_application"],
            "input": observed.get("InputFileAMPS") == Path(input_deck).name,
            "saved_spice": observed.get("SPICE") == "nospice",
            "staged_input": selected_input.is_file() and
                selected_input.read_bytes() == (repo / input_deck).read_bytes(),
        }
        result.receipt["configuration_sha256"] = {
            name: hashlib.sha256((repo / name).read_bytes()).hexdigest()
            for name in (".amps.conf", "Makefile.local", "Makefile.def", "Makefile.conf")
            if (repo / name).is_file()}

        target = "amps_test_native_configuration_" + uuid.uuid4().hex
        probe = work / (label + "-configuration.mk")
        probe.write_text(
            "$(info AMPS_NATIVE_CONFIG_SPICE=$(SPICE))\n"
            "$(info AMPS_NATIVE_CONFIG_SEARCH_C=$(SEARCH_C))\n"
            "$(info AMPS_NATIVE_CONFIG_FLAGCC=$(FLAGCC))\n"
            ".PHONY: " + target + "\n" + target + ":\n\t@:\n",
            encoding="utf-8")
        probe_argv = command(context, "make") + [
            "--no-print-directory", "-f", str(repo / "Makefile"),
            "-f", str(probe), target]
        result.receipt["probe_command"] = probe_argv
        inspected = run_command(probe_argv, repo, timeout=timeout,
                                env={"MAKEFLAGS": "-j16"})
        result.receipt["probe_returncode"] = inspected.returncode
        with log.open("a", encoding="utf-8") as stream:
            stream.write(inspected.stdout or "")
        effective: Dict[str, str] = {}
        for key in ("SPICE", "SEARCH_C", "FLAGCC"):
            values = re.findall(r"^AMPS_NATIVE_CONFIG_" + key + r"=(.*)$",
                                inspected.stdout or "", re.MULTILINE)
            if len(values) == 1:
                effective[key] = values[0].strip()
        result.receipt["effective_make"] = effective
        checks.update({
            "make_probe": inspected.returncode == 0,
            "effective_spice": effective.get("SPICE") == "nospice",
            "search_define": _no_spice_define(effective.get("SEARCH_C", "")),
            "compiler_define": _no_spice_define(effective.get("FLAGCC", "")),
        })
        result.receipt["checks"] = checks
        if not all(checks.values()):
            result.reject("native no-SPICE configuration mismatch: " +
                          ", ".join(key for key, passed in checks.items() if not passed))
        result.receipt["status"] = "PASS"
        result.save()
        print("[preflight] native_configuration=no-SPICE verified receipt=" +
              str(receipt_path), flush=True)
        return result
    except NativeConfigurationError:
        raise
    except Exception as error:
        result.reject("native configuration inspection failed: " + str(error))
        raise  # unreachable after reject; preserves static return analysis


def verify_no_spice_build(configuration: NativeConfiguration, output: str) -> None:
    """Check the real clean-build transcript before any MPI execution.

    The probe establishes the intended mode. Checking compiler commands as
    well catches a recursive-build override that disagrees with that probe.
    A clean native fixture always compiles C++ sources, so an empty collection
    is missing evidence rather than an implicit pass.
    """
    compiled = []
    libraries = []
    for line in output.splitlines():
        try:
            tokens = shlex.split(line)
        except ValueError:
            continue  # non-command diagnostics may contain unmatched quotes
        if "-c" in tokens and any(token.endswith((".cpp", ".cc", ".cxx"))
                                   for token in tokens):
            compiled.append(_no_spice_define(line))
        if any(token == "-lcspice" or token.endswith("cspice.a") for token in tokens):
            libraries.append(line)
    checks = {"cpp_commands_present": bool(compiled),
              "all_cpp_commands_no_spice": bool(compiled) and all(compiled),
              "no_cspice_archive_or_link": not libraries}
    configuration.receipt["build_checks"] = checks
    configuration.receipt["cpp_compile_command_count"] = len(compiled)
    configuration.receipt["unexpected_cspice_commands"] = libraries
    if not all(checks.values()):
        configuration.reject("native build contradicts the no-SPICE configuration: " +
                             ", ".join(key for key, passed in checks.items() if not passed))
    configuration.save()

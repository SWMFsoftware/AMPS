#!/usr/bin/env python3
"""Build the configured core-test application and run native PIC cases.

This is an extended, exclusive test because AMPS' established native build
workflow selects one repository-wide application.  It obeys the mandatory
preflight: confirm the repository, reject detected concurrent builds/runs,
remove only ``build``, select the application through Config.pl, and compile
with ``MAKEFLAGS=-j16``.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


PBUFFER_CASES = [
    "PIC-PBUFFER-INIT-001", "PIC-PBUFFER-STATE-001",
    "PIC-PBUFFER-LIST-001",
]

PHYSICS_CASES = [
    "PHYS-PHOTOCHEM-RATES-001",
    "PHYS-PHOTOCHEM-DECAY-001",
    "PHYS-PHOTOCHEM-KINEMATICS-001",
    "PHYS-PHOTOCHEM-BRANCHING-001",
    "PHYS-COLLISION-CONSERVATION-001",
    "PHYS-COLLISION-ISOTROPY-001",
    "PHYS-DISTRIBUTION-MAXWELLIAN-001",
    "PHYS-IDF-STORAGE-001",
    "PHYS-IDF-ROTATION-001",
    "PHYS-IDF-VIBRATION-001",
    "PHYS-IDF-PARTITION-001",
    "PHYS-IDF-RT-001",
    "PHYS-IDF-RT-RATE-001",
    "PHYS-MOVER-SIMPLE-001",
    "PHYS-MOVER-BORIS-001",
    "PHYS-MOVER-SURFACE-001",
    "PHYS-3D-COLLISION-RELAXATION-001",
    "PHYS-3D-COLLISION-RATE-001",
    "PHYS-3D-PHOTODECAY-001",
]

SPECIES_CASES = [
    "PIC-SPECIES-CONSTANTS-001",
    "PIC-SPECIES-GENERATED-MAPPING-001",
    "PIC-SPECIES-GENERATED-TABLES-001",
    "PIC-SPECIES-LOOKUP-DIAGNOSTICS-001",
]

UNIT_CASES = [
    "PIC-UNITS-GLOBAL-INIT-001",
    "PIC-UNITS-CONVERSION-TABLE-001",
]

WEIGHT_TIME_CASES = [
    "PIC-WEIGHT-MODE-ROUNDTRIP-001",
    "PIC-TIMESTEP-MODE-ROUNDTRIP-001",
    "PIC-WEIGHT-TIMESTEP-SPECIES-COPY-001",
    "PIC-WEIGHT-TIMESTEP-INVALID-VALUE-001",
    "PIC-WEIGHT-TIMESTEP-INVALID-SPECIES-001",
]

ADAPTER_CASES = [
    "CONFIG-GENERATOR-001", "CONFIG-INVALID-DECK-001",
    "CONFIG-DUPLICATE-SPECIES-001",
]
NATIVE_CASES = (PBUFFER_CASES + PHYSICS_CASES + SPECIES_CASES + UNIT_CASES +
                WEIGHT_TIME_CASES)

# The manifest aliases the historical family IDs to the combined suite, but
# accepting direct legacy contexts keeps this adapter backward-compatible.
SUITES: Dict[str, Dict[str, Any]] = {
    "CORE-PIC-NATIVE-001": {
        "artifact": "native-core-output.txt",
        "cases": ADAPTER_CASES + NATIVE_CASES,
    },
    "CORE-PIC-PARTICLE-BUFFER-001": {
        "artifact": "native-particle-buffer-output.txt",
        "cases": PBUFFER_CASES,
    },
    "CORE-PIC-PHYSICS-001": {
        "artifact": "native-physics-output.txt",
        "cases": PHYSICS_CASES,
    },
    "CORE-SPECIES-001": {
        "artifact": "native-species-output.txt",
        "cases": ["CONFIG-DUPLICATE-SPECIES-001"] + SPECIES_CASES,
    },
    "CORE-PIC-WEIGHT-TIME-001": {
        "artifact": "native-weight-time-output.txt",
        "cases": WEIGHT_TIME_CASES,
    },
}


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 identity for retained evidence."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(str(temporary), str(path))


def retain_build_identity(repo: Path, artifacts: Path) -> Tuple[Path, Path]:
    """Copy the executable/generated configuration and bind each to SHA-256."""

    retained = artifacts / "native-build"
    retained.mkdir(parents=True, exist_ok=True)
    candidates = [
        repo / "amps",
        repo / ".amps.conf",
        repo / "build" / "Makefile.conf",
        repo / "build" / "pic" / "picGlobal.dfn",
        repo / "build" / "pic" / "picSpeciesMacro.dfn",
        repo / "build" / "pic" / "pic.h",
        repo / "build" / "meshAMR" / "meshAMRdef.h",
        repo / "input" / "test" / "core-unit-tests.input",
    ]
    records: List[Dict[str, Any]] = []
    retained_executable = retained / "amps-core-unit-tests"
    for source in candidates:
        if not source.is_file():
            continue
        target = retained_executable if source == repo / "amps" else (
            retained / source.relative_to(repo).as_posix().replace("/", "__"))
        shutil.copy2(str(source), str(target))
        records.append({
            "repository_path": source.relative_to(repo).as_posix(),
            "retained_path": str(target),
            "size_bytes": target.stat().st_size,
            "sha256": sha256(target),
        })
    identity = artifacts / "native-build-identity.json"
    write_json(identity, {
        "schema_version": 1,
        "application": "test/core-unit-tests",
        "spice_path": "nospice",
        "configure_commands": [
            ["./Config.pl", "-application=test/core-unit-tests"],
            ["./Config.pl", "-spice-path=nospice"],
        ],
        "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
        "files": records,
    })
    return retained_executable, identity


def validate_generated_configuration(repo: Path) -> Tuple[bool, str]:
    """Check application identity and effective mutually exclusive modes.

    ``picGlobal.dfn`` intentionally contains defaults followed by generated
    ``#undef``/``#define`` overrides.  Searching for an expected token can
    therefore produce a false PASS when a later override selects a different
    mode.  Replaying those object-macro directives gives the effective value
    seen by the compiler and proves that each single-choice mode ends with the
    value selected by the input deck.
    """

    required = [repo / ".amps.conf", repo / "amps",
                repo / "build" / "pic" / "picGlobal.dfn",
                repo / "build" / "pic" / "pic.h",
                repo / "build" / "meshAMR" / "meshAMRdef.h"]
    missing = [path.relative_to(repo).as_posix() for path in required
               if not path.is_file()]
    if missing:
        return False, "missing generated artifacts: {}".format(", ".join(missing))
    selection = (repo / ".amps.conf").read_text(encoding="utf-8", errors="replace")
    definitions = (repo / "build" / "pic" / "picGlobal.dfn").read_text(
        encoding="utf-8", errors="replace")
    header = (repo / "build" / "pic" / "pic.h").read_text(
        encoding="utf-8", errors="replace")
    effective: Dict[str, str] = {}
    for raw in definitions.splitlines():
        undef = re.match(r"\s*#\s*undef\s+([A-Za-z_][A-Za-z0-9_]*)", raw)
        if undef:
            effective.pop(undef.group(1), None)
            continue
        define = re.match(
            r"\s*#\s*define\s+([A-Za-z_][A-Za-z0-9_]*)\s+([^/\s]+)", raw)
        if define:
            effective[define.group(1)] = define.group(2)

    expected_modes = {
        "_SIMULATION_TIME_STEP_MODE_":
            "_SPECIES_DEPENDENT_GLOBAL_TIME_STEP_",
        "_SIMULATION_PARTICLE_WEIGHT_MODE_":
            "_SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_",
        "_INDIVIDUAL_PARTICLE_WEIGHT_MODE_":
            "_INDIVIDUAL_PARTICLE_WEIGHT_OFF_",
        "_PIC_PHOTOLYTIC_REACTIONS_MODE_":
            "_PIC_PHOTOLYTIC_REACTIONS_MODE_ON_",
        "_PIC__PARTICLE_COLLISION_MODEL__MODE_": "_PIC_MODE_ON_",
        "_PIC__PARTICLE_COLLISION_MODEL_":
            "_PIC__PARTICLE_COLLISION_MODEL__NTC_",
        "_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_":
            "_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_",
    }
    checks = {
        "APPLICATION=core-unit-tests": "APPLICATION=core-unit-tests" in selection,
        "InputFileAMPS=core-unit-tests.input":
            "InputFileAMPS=core-unit-tests.input" in selection,
        "SPICE=nospice": "SPICE=nospice" in selection,
        "H2O species generated": bool(re.search(
            r"#define\s+_H2O_SPEC_\s+0(?:\s|$)", header)),
    }
    checks.update({
        "effective {}={}".format(name, value): effective.get(name) == value
        for name, value in expected_modes.items()
    })
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        return False, "generated configuration checks failed: {}".format(
            ", ".join(failed))
    return True, "selected application and generated physics macros verified"


def active_repo_processes(repo: Path) -> List[str]:
    """Best-effort Linux preflight, supplemented by the portable file lock.

    The /proc scan protects unrelated, non-framework AMPS runs that cannot
    know about our lock.  On systems without /proc, the advisory lock still
    serializes instances of this runner; site schedulers must serialize legacy
    builds that do not participate in that lock.
    """

    proc = Path("/proc")
    if not proc.is_dir():
        return []
    active = []
    own_pid = os.getpid()
    names = {"make", "gmake", "amps", "Config.pl", "ampsConfig.pl"}
    for entry in proc.iterdir():
        if not entry.name.isdigit() or int(entry.name) == own_pid:
            continue
        try:
            cwd = (entry / "cwd").resolve()
            cwd.relative_to(repo)
            raw = (entry / "cmdline").read_bytes().replace(b"\x00", b" ").decode(
                "utf-8", errors="replace").strip()
        except (OSError, ValueError):
            continue
        executable = Path(raw.split()[0]).name if raw else ""
        if (executable in names or executable.startswith("amps") or
                any(token in raw for token in ("Config.pl", "ampsConfig.pl"))):
            active.append("pid={} cwd={} command={}".format(entry.name, cwd, raw))
    return active


def main() -> int:
    context = load_context("Run configured native AMPS core tests")
    suite_id = str(context["test"]["id"])
    if suite_id not in SUITES:
        raise RunnerError("native core runner does not define suite {}".format(suite_id))
    suite = SUITES[suite_id]
    requested = list(context["test"].get("selected_cases", suite["cases"]))
    if (not requested or any(not isinstance(value, str) for value in requested) or
            len(set(requested)) != len(requested)):
        raise RunnerError("selected_cases must be a unique nonempty case list")
    unknown = [value for value in requested if value not in suite["cases"]]
    if unknown:
        raise RunnerError("suite {} does not own selected cases {}".format(
            suite_id, unknown))
    repo = Path(context["repo_root"]).resolve()
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 1800))
    if not (repo / "Config.pl").is_file() or not (repo / "src").is_dir():
        return write_result(context, "ERROR", "repository preflight failed: {}".format(repo))
    actual_cwd = Path.cwd().resolve()
    if actual_cwd != repo:
        return write_result(
            context, "ERROR",
            "native adapter cwd {} is not repository {}".format(actual_cwd, repo))
    print("[preflight] repository={} cwd={} status=confirmed".format(
        repo, actual_cwd), flush=True)

    # Lock a stable inode inside this checkout rather than a TMPDIR file.
    # Cluster nodes and site profiles can have different local temporary
    # directories while sharing the source tree; locking Config.pl therefore
    # coordinates every instance that can mutate this checkout and naturally
    # leaves independent checkouts unconstrained. The file is opened read-only
    # and its contents/metadata are never changed.
    lock_path = repo / "Config.pl"
    lock_stream = lock_path.open("r")
    try:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return write_result(context, "ERROR", "another native AMPS test holds {}".format(lock_path))
        active = active_repo_processes(repo)
        if active:
            return write_result(context, "ERROR", "native build/run already active: {}".format(
                "; ".join(active)))
        print("[preflight] active_native_build_or_test_processes=0", flush=True)

        # Exact, validated generated target required by the repository build
        # policy.  No source, configuration, library, or parent path is removed.
        build = (repo / "build").resolve()
        if build != repo / "build":
            raise RunnerError("resolved build target is not repository/build")
        print("[preflight] generated_removal_target={} status=validated".format(
            build), flush=True)
        removed = run_command(["rm", "-rf", "--", "build"], repo, timeout=60)
        if removed.returncode != 0:
            return write_result(context, "ERROR", "could not remove generated build directory")

        # Config.pl is the sole application-selection authority.  In
        # particular, this runner never edits .amps.conf or Makefile.local.
        configured = run_command(
            [str(repo / "Config.pl"), "-application=test/core-unit-tests"],
            repo, timeout=timeout)
        if configured.returncode != 0:
            return write_result(context, "FAIL", "Config.pl could not select core-unit-tests",
                                cases=[case(identifier, "FAIL", "configuration failed")
                                       for identifier in requested])
        spice = run_command([str(repo / "Config.pl"), "-spice-path=nospice"],
                            repo, timeout=timeout)
        if spice.returncode != 0:
            return write_result(context, "FAIL", "Config.pl could not disable SPICE",
                                cases=[case(identifier, "FAIL", "SPICE configuration failed")
                                       for identifier in requested])

        # Repository policy fixes native rebuild concurrency at -j16.  The
        # profile's build_jobs remains available to dependency-light adapters
        # but intentionally cannot alter this native build gate.
        make_command = command(context, "make") + ["amps"]
        print("[preflight] build_environment=MAKEFLAGS=-j16 command={}".format(
            " ".join(make_command)), flush=True)
        built = run_command(make_command, repo, timeout=timeout,
                            env={"MAKEFLAGS": "-j16"})
        if built.returncode != 0 or not (repo / "amps").is_file():
            return write_result(context, "FAIL", "native core-test application did not build",
                                cases=[case(identifier, "FAIL",
                                            "make exit {}".format(built.returncode))
                                       for identifier in requested])

        artifacts.mkdir(parents=True, exist_ok=True)
        executable, identity = retain_build_identity(repo, artifacts)
        evidence = [str(identity), str(executable)]
        observed: Dict[str, Dict[str, Any]] = {}

        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        # These are process-local numerical/module tests.  The mover fixture
        # constructs a single owner block and is intentionally a one-rank
        # analytic test; multi-rank ownership/ghost tests are separate future
        # suites.  The manifest records this override explicitly so a machine
        # profile requesting four ranks cannot silently change its semantics.
        ranks = str(context["test"].get(
            "mpi_ranks", context["execution"].get("mpi_ranks", 1)))
        native_requested = [identifier for identifier in requested
                            if identifier in NATIVE_CASES]
        transcript_chunks: List[str] = []
        # Execute the already-final context sequence exactly, including the
        # adapter-level configuration cases.  Each native case gets a fresh
        # process while sharing one retained build. This is hard isolation for
        # PIC's persistent mesh, spherical-boundary registration, random
        # generator, and other global singletons. Applying --order here would
        # double-transform reverse/shuffle selected_cases from the root runner.
        for identifier in requested:
            if identifier == "CONFIG-GENERATOR-001":
                print("[CASE] RUN  CONFIG-GENERATOR-001 generated application identity and physics macros",
                      flush=True)
                # Verify both the human-readable selection and the generated
                # headers used to build the retained executable. Hash evidence
                # already exists before this case publishes its result.
                passed, message = validate_generated_configuration(repo)
                observed[identifier] = case(
                    identifier, "PASS" if passed else "FAIL", message)
                continue

            if identifier == "CONFIG-INVALID-DECK-001":
                print("[CASE] RUN  CONFIG-INVALID-DECK-001 absent input deck rejection",
                      flush=True)
                suffix = 0
                while True:
                    missing = work / (
                        "guaranteed-absent-core-deck-{}-{}.input".format(
                            os.getpid(), suffix))
                    if not missing.exists():
                        break
                    suffix += 1
                invalid = run_command(
                    [str(repo / "ampsConfig.pl"), "-input", str(missing),
                     "-no-compile"], repo, timeout=min(timeout, 120))
                invalid_transcript = artifacts / "invalid-deck-output.txt"
                invalid_transcript.write_text(
                    invalid.stdout or "", encoding="utf-8")
                evidence.append(str(invalid_transcript))
                diagnostic = "Input file must be specified with -input and must exist"
                passed = (invalid.returncode != 0 and
                          diagnostic in (invalid.stdout or ""))
                message = (
                    "missing input rejected before generated-tree mutation"
                    if passed else
                    "expected nonzero exit and missing-input diagnostic; exit={}".format(
                        invalid.returncode))
                observed[identifier] = case(
                    identifier, "PASS" if passed else "FAIL", message)
                continue

            if identifier == "CONFIG-DUPLICATE-SPECIES-001":
                print("[CASE] RUN  CONFIG-DUPLICATE-SPECIES-001 unique compiled names and repeated DUST groups",
                      flush=True)
                source_deck = repo / "input" / "test" / "core-unit-tests.input"
                invalid_deck = work / "duplicate-species.input"
                dust_deck = work / "repeated-dust-groups.input"
                original = source_deck.read_text(encoding="utf-8")
                duplicate = original.replace(
                    "SpeciesList=H2O,OH,H,O,H2,O2,CO2",
                    "SpeciesList=H2O,OH,H2O,O,H2,O2,CO2", 1)
                repeated_dust = original.replace(
                    "SpeciesList=H2O,OH,H,O,H2,O2,CO2",
                    "SpeciesList=H2O,DUST,DUST,DUST", 1)
                if duplicate == original or repeated_dust == original:
                    observed[identifier] = case(
                        identifier, "FAIL", "could not construct SpeciesList fixtures")
                    continue
                invalid_deck.write_text(duplicate, encoding="utf-8")
                dust_deck.write_text(repeated_dust, encoding="utf-8")
                before = sha256(executable)
                invalid = run_command(
                    [str(repo / "ampsConfig.pl"), "-input", str(invalid_deck),
                     "-validate-species-list-only"], repo,
                    timeout=min(timeout, 120))
                dust = run_command(
                    [str(repo / "ampsConfig.pl"), "-input", str(dust_deck),
                     "-validate-species-list-only"], repo,
                    timeout=min(timeout, 120))
                duplicate_transcript = artifacts / "duplicate-species-output.txt"
                duplicate_transcript.write_text(
                    "===== duplicate ordinary species =====\n{}\n"
                    "===== repeated DUST groups =====\n{}".format(
                        invalid.stdout or "", dust.stdout or ""),
                    encoding="utf-8")
                evidence.append(str(duplicate_transcript))
                diagnostic = "compiled species names must be unique"
                passed = (invalid.returncode != 0 and
                          diagnostic in (invalid.stdout or "") and
                          dust.returncode == 0 and
                          "SpeciesList validation passed" in (dust.stdout or "") and
                          executable.is_file() and sha256(executable) == before)
                message = (
                    "ordinary duplicate rejected, repeated DUST groups accepted, "
                    "and configured executable unchanged"
                    if passed else
                    "expected duplicate rejection, repeated-DUST acceptance, and "
                    "unchanged executable; duplicate_exit={} dust_exit={}".format(
                        invalid.returncode, dust.returncode))
                observed[identifier] = case(
                    identifier, "PASS" if passed else "FAIL", message)
                continue

            run = run_command(command(context, "mpiexec") +
                              [rank_flag, ranks, str(executable),
                               "--case", identifier],
                              work, timeout=timeout)
            transcript_chunks.append(
                "===== {} =====\n{}".format(identifier, run.stdout or ""))
            found = re.findall(
                r"^\[([A-Z][A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                run.stdout or "", flags=re.MULTILINE)
            exact = [entry for entry in found if entry[0] == identifier]
            foreign = [entry[0] for entry in found if entry[0] != identifier]
            if len(exact) == 1 and not foreign:
                _, status, message = exact[0]
                if run.returncode != 0 and status == "PASS":
                    status = "FAIL"
                    message = "PASS record with mpiexec exit {}".format(run.returncode)
                observed[identifier] = case(identifier, status, message or "ok")
            else:
                observed[identifier] = case(
                    identifier, "FAIL",
                    "expected one exact record; observed exact={} foreign={} exit={}".format(
                        len(exact), foreign, run.returncode))

        transcript = artifacts / str(suite["artifact"])
        transcript.write_text("\n".join(transcript_chunks), encoding="utf-8")
        evidence.append(str(transcript))
        missing_results = [identifier for identifier in requested
                           if identifier not in observed]
        for identifier in missing_results:
            observed[identifier] = case(identifier, "FAIL", "case was not executed")
        cases = [observed[identifier] for identifier in requested]
        status = "PASS" if all(item["status"] == "PASS" for item in cases) else "FAIL"
        return write_result(
            context, status,
            "{} selected cases executed from one isolated native build".format(
                len(cases)),
            cases=cases, evidence=evidence,
            metrics={
                "mpi_ranks": int(ranks),
                "native_processes": len(native_requested),
                "case_order": str(context["test"].get("case_order", "normal")),
                "shuffle_seed": int(context["test"].get("shuffle_seed", 0)),
                "executable_sha256": sha256(executable),
            })
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)

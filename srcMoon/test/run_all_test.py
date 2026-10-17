#!/usr/bin/env python3
"""Run every registered srcMoon local and M1/M2 linked test and summarize it.

Each test is launched through its own command-line runner. U01--U26 consume the
generated build tree.  The default REGRESSION campaign uses a deliberately
coarse production receipt for one bounded I12 linked wiring run; I11 and I13
report SKIPPED because an ordinary regression must not silently launch the
costly M1 convergence/freeze campaign. QUALIFICATION explicitly runs the
reviewed, 5000-step-bounded I11 sequence before its dependent I12--I13 gates.
The aggregate runner preserves its outer
runner log and reads the individual structured result.json to determine
scientific status. It never infers PASS from a zero process exit alone.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any


TEST_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TEST_ROOT.parents[1]
DEFAULT_BUILD = REPO_ROOT / "build"
DEFAULT_OUTPUT_PARENT = REPO_ROOT / "test_output" / "srcMoon"
DEFAULT_SPICE_PATH = Path("/home/vtenishe/SPICE/cspice")
DEFAULT_SPICE_KERNELS = Path("/home/vtenishe/SPICE/Kernels")
DEFAULT_AMPS = REPO_ROOT / "amps"
DEFAULT_M2_BUILD_DIR_NAME = "_m2_na_ar_build"
DEFAULT_M1_INPUT = REPO_ROOT / "srcMoon" / "examples" / "m1_memory200.in"
DEFAULT_REGRESSION_INPUT = (
    REPO_ROOT / "srcMoon" / "examples" / "m1_regression.in"
)
DEFAULT_M2_INPUT = REPO_ROOT / "srcMoon" / "examples" / "lola_surface.in"

# This is the campaign registry, not a directory-discovery convenience.  A new
# test must be added here deliberately so code review can see that the aggregate
# campaign changed. validate_registry() also rejects an unregistered UXX
# directory and enforces the exact requested I11--I13 linked range, preventing
# a newly added test from being silently omitted.
U_TESTS: tuple[tuple[str, str], ...] = (
    ("U01", "U01_los_geometry/test.py"),
    ("U02", "U02_target_isolation/test.py"),
    ("U03", "U03_lunar_gravity/test.py"),
    ("U04", "U04_rotating_frame/test.py"),
    ("U05", "U05_na_radiation_shadow/test.py"),
    ("U06", "U06_lorentz/test.py"),
    ("U07", "U07_surface_refinement/test.py"),
    ("U08", "U08_lola_geometry/test.py"),
    ("U09", "U09_terrain_shadow/test.py"),
    ("U10", "U10_diviner_temperature/test.py"),
    ("U11", "U11_thermal_inertia/test.py"),
    ("U12", "U12_surface_interaction/test.py"),
    ("U13", "U13_cold_trapping/test.py"),
    ("U14", "U14_na_photoionization/test.py"),
    ("U15", "U15_multispecies_photoionization/test.py"),
    ("U16", "U16_electron_impact/test.py"),
    ("U17", "U17_charge_exchange/test.py"),
    ("U18", "U18_plasma_driver/test.py"),
    ("U19", "U19_helium_source/test.py"),
    ("U20", "U20_neon_source/test.py"),
    ("U21", "U21_argon_source/test.py"),
    ("U22", "U22_na_sources/test.py"),
    ("U23", "U23_meteoroid_driver/test.py"),
    ("U24", "U24_water_chemistry/test.py"),
    ("U25", "U25_data_provenance/test.py"),
    ("U26", "U26_observation_geometry/test.py"),
)

# M1 linked tests use the same executable and reviewed application receipt but
# retain separate runners and artifacts. They are registered explicitly here
# so adding a directory cannot silently omit it from the aggregate campaign.
# I11 has no statistical time/block ceiling in QUALIFICATION mode: unresolved
# confidence intervals keep its production campaign running. REGRESSION mode
# still invokes all three runners, but I11/I13 return explicit SKIPPED records
# without starting their expensive linked models. I12 performs the one bounded
# coarse production initialization/wiring run. The aggregate preserves every
# runner status and never converts successful process launch into PASS.
LINKED_M1_TESTS: tuple[tuple[str, str], ...] = (
    ("I11", "integration-tests/I11/test.py"),
    ("I12", "integration-tests/I12/test.py"),
    ("I13", "integration-tests/I13/test.py"),
)

# M2 IDs are always invoked after I13.  Their individual runners enforce the
# predecessor result and retain a structured SKIPPED artifact when the FULL M1
# baseline has not passed; local U09--U13 success never substitutes for that
# linked gate.
LINKED_M2_TESTS: tuple[tuple[str, str], ...] = tuple(
    (f"I{number:02d}", f"integration-tests/I{number:02d}/test.py")
    for number in range(14, 21)
)

TESTS = U_TESTS + LINKED_M1_TESTS + LINKED_M2_TESTS

VALID_STATUSES = ("PASS", "FAIL", "SKIPPED", "ERROR")
EXPECTED_EXIT = {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}
REQUIRED_LINKED_ARTIFACTS = {
    "stdout",
    "stderr",
    "receipt",
    "result",
    "run_manifest",
    "input_snapshot",
}


def validate_registry() -> list[str]:
    """Return registry errors without running any test or modifying outputs."""
    errors: list[str] = []
    ids = [test_id for test_id, _ in TESTS]
    if len(ids) != len(set(ids)):
        errors.append("the run_all_test.py registry contains duplicate test IDs")

    registered_u_directories: set[str] = set()
    for test_id, runner_name in TESTS:
        runner = TEST_ROOT / runner_name
        if not runner.is_file():
            errors.append(f"{test_id}: registered runner does not exist: {runner}")
        if runner.parts[-2] != test_id and not runner.parts[-2].startswith(
            test_id + "_"
        ):
            errors.append(
                f"{test_id}: runner directory does not match its ID: {runner_name}"
            )
        if test_id.startswith("U"):
            registered_u_directories.add(runner.parts[-2])

    discovered_directories = {
        path.name
        for path in TEST_ROOT.glob("U[0-9][0-9]_*")
        if path.is_dir()
    }
    for directory in sorted(discovered_directories - registered_u_directories):
        errors.append(
            f"unregistered test directory {directory}; add its runner to TESTS"
        )
    for directory in sorted(registered_u_directories - discovered_directories):
        errors.append(f"registered test directory is missing: {directory}")

    # I01--I10 belong to the separately configured M0 baseline campaign. Only
    # the requested M1/M2 range is checked here, avoiding any implication that a
    # single SPICE build can replace the dedicated M0 multi-build runner.
    expected_m1_directories = {f"I{number:02d}" for number in range(11, 14)}
    registered_m1_directories = {
        Path(runner).parts[-2] for _, runner in LINKED_M1_TESTS
    }
    if registered_m1_directories != expected_m1_directories:
        errors.append(
            "linked M1 registry must contain exactly I11, I12, and I13"
        )
    expected_m2_directories = {f"I{number:02d}" for number in range(14, 21)}
    registered_m2_directories = {
        Path(runner).parts[-2] for _, runner in LINKED_M2_TESTS
    }
    if registered_m2_directories != expected_m2_directories:
        errors.append("linked M2 registry must contain exactly I14 through I20")
    return errors


def git_revision() -> str:
    """Return the source revision recorded in the aggregate result."""
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def default_output_directory() -> Path:
    """Create a run-specific default name so prior campaign evidence survives."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return DEFAULT_OUTPUT_PARENT / f"run-all-{timestamp}"


def run_logged(
    command: list[str], log_path: Path, stream: bool, append: bool = False
) -> int:
    """Run one command with a durable combined stdout/stderr transcript.

    ``append`` is used only for a multi-step operation such as configuration
    followed by compilation.  Individual test-runner logs retain the safer
    default of replacing an old transcript rather than mixing two executions.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a" if append else "w", encoding="utf-8") as log:
        if append:
            log.write("\n")
        log.write("command: " + " ".join(command) + "\n")
        log.flush()
        try:
            if not stream:
                completed = subprocess.run(
                    command,
                    cwd=REPO_ROOT,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
                return completed.returncode

            # Streaming is optional because most tests are short.  When
            # enabled, each line is copied to the terminal and durable log.
            process = subprocess.Popen(
                command,
                cwd=REPO_ROOT,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)
            return process.wait()
        except OSError as error:
            # Treat an unlaunchable compiler/interpreter/runner as
            # infrastructure ERROR and retain the operating-system diagnostic.
            message = f"ERROR: could not execute command: {error}\n"
            log.write(message)
            log.flush()
            if stream:
                print(message, end="", file=sys.stderr, flush=True)
            return 127


def clean_build(
    output_root: Path,
    stream: bool,
    spice_path: Path,
    spice_kernels: Path,
) -> dict[str, Any]:
    """Clean, configure the Moon application with SPICE, and run ``make -j``.

    Config.pl changes the generated AMPS build contract, so configuration is a
    required part of this campaign rather than an undocumented prerequisite.
    The build tree is removed before Config.pl runs, exactly matching the
    repository build discipline.  Both commands share one durable log and are
    also recorded separately in the machine-readable campaign summary.
    """
    if DEFAULT_BUILD.exists():
        # AMPS generates configured source beneath build/.  A clean tree is
        # required so the test campaign cannot link against stale copies.
        shutil.rmtree(DEFAULT_BUILD)
    log_path = output_root / "build.log"
    started = time.monotonic()

    configure_command = [
        str(REPO_ROOT / "Config.pl"),
        "-application=moon",
        f"-spice-path={spice_path}",
        f"-spice-kernels={spice_kernels}",
    ]
    configure_started = time.monotonic()
    configure_returncode = run_logged(configure_command, log_path, stream)
    configure_result = {
        "name": "configure",
        "command": configure_command,
        "returncode": configure_returncode,
        "duration_seconds": time.monotonic() - configure_started,
    }
    if configure_returncode != 0:
        return {
            "status": "ERROR",
            "failed_step": "configure",
            "returncode": configure_returncode,
            "log": str(log_path),
            "steps": [configure_result],
            "duration_seconds": time.monotonic() - started,
        }

    # Debug symbols and compiler optimization are independent.  Passing OPT3
    # on make's command line overrides the generated Makefile.conf value and
    # propagates through recursive makes.  This keeps GDB-capable symbols while
    # avoiding the effectively -O0 production executable produced when OPT3 is
    # replaced by bare ``-g``.  Makefile.conf and Makefile remain generated,
    # local files and are never edited or proposed for commit by this runner.
    build_command = ["make", "-j", "OPT3=-O3 -g"]
    build_started = time.monotonic()
    build_returncode = run_logged(build_command, log_path, stream, append=True)
    build_result = {
        "name": "build",
        "command": build_command,
        "returncode": build_returncode,
        "duration_seconds": time.monotonic() - build_started,
    }
    return {
        "status": "PASS" if build_returncode == 0 else "ERROR",
        "failed_step": None if build_returncode == 0 else "build",
        "returncode": build_returncode,
        "log": str(log_path),
        "steps": [configure_result, build_result],
        "duration_seconds": time.monotonic() - started,
    }


def build_m2_qualification(
    output_root: Path,
    stream: bool,
    spice_path: Path,
    spice_kernels: Path,
) -> dict[str, Any]:
    """Create the separate Na+Ar executable before restoring canonical build.

    The helper performs its own clean Config.pl/build cycle and copies the
    executable into the immutable campaign tree.  ``clean_build`` is called
    immediately afterward, so U-series and M0/M1 tests still consume the
    canonical Na-only generated tree and ``amps`` executable.
    """
    artifact_dir = output_root / DEFAULT_M2_BUILD_DIR_NAME
    command = [
        sys.executable,
        str(TEST_ROOT / "build_m2_qualification.py"),
        "--output-dir", str(artifact_dir),
        "--spice-path", str(spice_path),
        "--spice-kernels", str(spice_kernels),
    ]
    outer_log = output_root / "build-m2-na-ar-runner.log"
    started = time.monotonic()
    returncode = run_logged(command, outer_log, stream)
    receipt = artifact_dir / "build-m2-na-ar.json"
    result: dict[str, Any] = {
        "status": "ERROR", "returncode": returncode,
        "command": command, "runner_log": str(outer_log),
        "receipt": str(receipt),
        "duration_seconds": time.monotonic() - started,
    }
    if returncode != 0 or not receipt.is_file():
        result["reason"] = (
            "dedicated Na+Ar Config.pl/build failed or omitted its receipt")
        return result
    try:
        payload = json.loads(receipt.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        result["reason"] = f"unreadable Na+Ar build receipt: {error}"
        return result
    executable = artifact_dir / "amps-moon-m2-na-ar"
    if payload.get("status") != "PASS" or not executable.is_file():
        result["reason"] = str(payload.get("error", "Na+Ar build did not pass"))
        result["reported_receipt"] = payload
        return result
    result.update(status="PASS", executable=str(executable),
                  executable_sha256=payload["executable"]["sha256"],
                  reported_receipt=payload)
    return result


def linked_artifact_index(result_path: Path) -> dict[str, str]:
    """Index the mandatory linked-run evidence that already exists on disk.

    The individual I-series runner owns these files. The aggregate runner only
    records paths; it neither copies nor rewrites scientific evidence. Checking
    existence here makes partial ERROR runs auditable without claiming that a
    missing mandatory artifact was produced successfully.
    """
    artifact_directory = result_path.parent
    expected_files = {
        "stdout": artifact_directory / "stdout.log",
        "stderr": artifact_directory / "stderr.log",
        "receipt": artifact_directory / "receipt.txt",
        "result": result_path,
        "run_manifest": artifact_directory / "run_manifest.json",
        "input_snapshot": artifact_directory / "application.in",
    }
    artifacts = {
        name: str(path) for name, path in expected_files.items() if path.is_file()
    }
    # Future linked scorers may add convergence tables or science diagnostics.
    # Index those files by their relative names so the aggregate result remains
    # complete without teaching this orchestrator their scientific meaning.
    if artifact_directory.is_dir():
        for path in sorted(artifact_directory.rglob("*")):
            if path.is_file() and str(path) not in artifacts.values():
                artifacts[str(path.relative_to(artifact_directory))] = str(path)
    return artifacts


def load_runner_result(
    test_id: str, runner: Path, result_path: Path, log_path: Path, returncode: int
) -> dict[str, Any]:
    """Validate one runner's JSON/exit-code protocol and normalize its record."""
    base = {
        "id": test_id,
        "runner": str(runner),
        "log": str(log_path),
        "result_json": str(result_path),
        "returncode": returncode,
    }
    indexed_artifacts: dict[str, str] = {}
    if test_id.startswith("I"):
        # The per-ID directory is itself evidence even when result.json is
        # absent or corrupt: integration_common.py may already have retained
        # stdout/stderr, the input snapshot, and a partial receipt there.
        # Record it before parsing so every ERROR path remains auditable.
        base["artifact_directory"] = str(result_path.parent)
        indexed_artifacts = linked_artifact_index(result_path)
        base["artifacts"] = indexed_artifacts
        if "run_manifest" in indexed_artifacts:
            base["run_manifest"] = indexed_artifacts["run_manifest"]
    try:
        with result_path.open(encoding="utf-8") as stream:
            payload = json.load(stream)
    except FileNotFoundError:
        return {
            **base,
            "status": "ERROR",
            "reason": "individual runner did not create its required result.json",
        }
    except (OSError, json.JSONDecodeError) as error:
        return {
            **base,
            "status": "ERROR",
            "reason": f"individual runner produced unreadable result.json: {error}",
        }

    status = payload.get("status")
    if payload.get("id") != test_id or status not in VALID_STATUSES:
        return {
            **base,
            "status": "ERROR",
            "reason": (
                "result.json must contain the registered id and one of "
                + ", ".join(VALID_STATUSES)
            ),
            "reported_result": payload,
            "artifacts": indexed_artifacts,
        }

    expected_exit = EXPECTED_EXIT[status]
    if returncode != expected_exit:
        return {
            **base,
            "status": "ERROR",
            "reason": (
                f"runner exit {returncode} disagrees with {status} "
                f"result (expected {expected_exit})"
            ),
            "reported_status": status,
            "name": payload.get("name", test_id),
            "artifacts": indexed_artifacts,
        }

    if test_id.startswith("I") and status != "ERROR":
        # REGRESSION deliberately does not launch the costly bounded I11 or the
        # prerequisite-blocked I13 model. Their manifests say execution=not-run
        # and no receipt can honestly exist. All actually linked runs retain
        # the complete evidence set, including the production receipt.
        not_run = False
        manifest_path = indexed_artifacts.get("run_manifest")
        if manifest_path is not None:
            try:
                not_run = json.loads(
                    Path(manifest_path).read_text(encoding="utf-8")
                ).get("execution") == "not-run"
            except (OSError, json.JSONDecodeError):
                not_run = False
        required = (
            REQUIRED_LINKED_ARTIFACTS - {"receipt"}
            if status == "SKIPPED" and not_run
            else REQUIRED_LINKED_ARTIFACTS
        )
        missing_artifacts = sorted(required - indexed_artifacts.keys())
        if missing_artifacts:
            return {
                **base,
                "status": "ERROR",
                "reason": (
                    f"{status} linked result is missing mandatory retained "
                    "artifact(s): " + ", ".join(missing_artifacts)
                ),
                "reported_status": status,
                "name": payload.get("name", test_id),
                "artifacts": indexed_artifacts,
            }

    normalized = {
        **base,
        "status": status,
        "name": payload.get("name", test_id),
        "reason": payload.get("reason"),
    }
    # Preserve discoverability of every linked artifact in the aggregate JSON
    # without copying or flattening the scientifically meaningful files. The
    # complete individual result remains authoritative at result_json.
    if test_id.startswith("I"):
        # Only files observed after runner completion are indexed. The older
        # individual schema can declare an intended receipt path before that
        # receipt exists, which must not be represented as retained evidence.
        normalized["artifacts"] = indexed_artifacts
    return normalized


def run_test(
    test_id: str,
    runner_name: str,
    build_dir: Path,
    output_root: Path,
    stream: bool,
    amps: Path,
    m1_input: Path,
    m2_input: Path,
    campaign_mode: str,
    mpi_ranks: int,
    m2_amps: Path | None,
) -> dict[str, Any]:
    """Execute one registered CLI runner and consume its structured result."""
    runner = (TEST_ROOT / runner_name).resolve()
    result_path = output_root / test_id / "result.json"
    if result_path.exists():
        # A stale artifact must never let a crashed runner appear successful.
        result_path.unlink()
    log_path = output_root / "runner_logs" / f"{test_id}.log"
    if test_id.startswith("I"):
        # integration_common.py owns the inner linked stdout/stderr, receipt,
        # manifest, and exact application-input snapshot. This outer log keeps
        # only the runner-level transcript and status JSON, so no artifact is
        # overwritten or mistaken for model output.
        linked_input = m2_input if int(test_id[1:]) >= 14 else m1_input
        linked_executable = (
            m2_amps if test_id in {"I18", "I19", "I20"} and
            m2_amps is not None else amps)
        command = [
            sys.executable,
            str(runner),
            "--amps",
            str(linked_executable),
            "--input",
            str(linked_input),
            "--output-dir",
            str(output_root),
            "--campaign-mode",
            campaign_mode,
            "--mpi-ranks",
            str(mpi_ranks),
        ]
    else:
        command = [
            sys.executable,
            str(runner),
            "--build-dir",
            str(build_dir),
            "--output-dir",
            str(output_root),
        ]
    started = time.monotonic()
    returncode = run_logged(command, log_path, stream)
    result = load_runner_result(test_id, runner, result_path, log_path, returncode)
    result["duration_seconds"] = time.monotonic() - started
    return result


def print_summary(counts: dict[str, int], results: list[dict[str, Any]]) -> None:
    """Print totals and mandatory diagnostics for FAIL/ERROR outcomes."""
    print("\nsrcMoon test summary")
    print(f"  executed: {counts['executed']}")
    print(f"  passed:   {counts['PASS']}")
    print(f"  failed:   {counts['FAIL']}")
    print(f"  skipped:  {counts['SKIPPED']}")
    print(f"  errors:   {counts['ERROR']}")

    problems = [result for result in results if result["status"] in {"FAIL", "ERROR"}]
    if problems:
        print("\nFailed/error tests:")
        for result in problems:
            print(f"  {result['status']} {result['id']} {result.get('name', '')}".rstrip())
            print(f"    runner: {result['runner']}")
            print(f"    log:    {result['log']}")
            if result.get("artifact_directory"):
                print(f"    linked artifacts: {result['artifact_directory']}")
            if result.get("reason"):
                print(f"    reason: {result['reason']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list the explicit registry")
    parser.add_argument(
        "--only", action="append", metavar="UXX|IXX",
        help="run only a registered ID",
    )
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument(
        "--amps",
        type=Path,
        default=DEFAULT_AMPS,
        help="SPICE Moon executable used by linked I11-I13 (default: %(default)s)",
    )
    parser.add_argument(
        "--m2-amps",
        type=Path,
        help=("prebuilt dedicated Na+Ar executable used by I18-I20 when "
              "--skip-build is selected; normal builds create it automatically"),
    )
    parser.add_argument(
        "--m1-input",
        type=Path,
        help=("Moon .in receipt used by linked I11-I13; the default is "
              "m1_regression.in for REGRESSION and m1_memory200.in for "
              "QUALIFICATION"),
    )
    parser.add_argument(
        "--m2-input", type=Path, default=DEFAULT_M2_INPUT,
        help="Moon LOLA/M2 .in receipt used by linked I14-I20",
    )
    parser.add_argument(
        "--campaign-mode",
        choices=("regression", "qualification"),
        default="regression",
        help=("bounded linked regression (default) or the 5000-step-bounded "
              "M1 qualification campaign"),
    )
    parser.add_argument(
        "--mpi-ranks",
        type=int,
        default=8,
        help=("positive MPI process count passed to every linked M1 runner "
              "and recorded in each manifest (default: 8)"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="new/empty artifact directory (default: UTC-stamped under test_output)",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help=(
            "reuse an existing generated build instead of clean Moon/SPICE "
            "configuration and make -j"
        ),
    )
    parser.add_argument(
        "--spice-path",
        type=Path,
        default=DEFAULT_SPICE_PATH,
        help="CSPICE installation used by Config.pl (default: %(default)s)",
    )
    parser.add_argument(
        "--spice-kernels",
        type=Path,
        default=DEFAULT_SPICE_KERNELS,
        help="SPICE kernel directory used by Config.pl (default: %(default)s)",
    )
    parser.add_argument("--stream", action="store_true", help="also stream runner logs")
    parser.add_argument(
        "--strict-skips",
        action="store_true",
        help="return failure when any test is SKIPPED",
    )
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("--mpi-ranks must be a positive integer")

    registry_errors = validate_registry()
    if registry_errors:
        for error in registry_errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2

    registry = dict(TESTS)
    if args.list:
        for test_id, runner in TESTS:
            print(f"{test_id}  {runner}")
        return 0

    selected = list(registry) if args.only is None else args.only
    unknown = [test_id for test_id in selected if test_id not in registry]
    if unknown:
        parser.error("unknown test ID(s): " + ", ".join(unknown))
    if len(selected) != len(set(selected)):
        parser.error("a test ID may be selected only once")

    output_root = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else default_output_directory().resolve()
    )
    build_dir = args.build_dir.resolve()
    amps = args.amps.resolve()
    selected_m1_input = args.m1_input
    if selected_m1_input is None:
        selected_m1_input = (
            DEFAULT_REGRESSION_INPUT
            if args.campaign_mode == "regression"
            else DEFAULT_M1_INPUT
        )
    m1_input = selected_m1_input.resolve()
    m2_input = args.m2_input.resolve()
    m2_amps = args.m2_amps.resolve() if args.m2_amps is not None else None
    if output_root.exists() and any(output_root.iterdir()):
        parser.error(
            f"output directory is not empty: {output_root}; choose a new directory "
            "so prior test evidence is not overwritten"
        )
    output_root.mkdir(parents=True, exist_ok=True)
    campaign_started = time.monotonic()

    # I18-I20 cannot launch unless this fresh campaign also executes the FULL
    # I13 predecessor.  In the default regression campaign, or when I13 is not
    # selected, their runners deliberately write not-run SKIPPED evidence.
    # Avoid spending two clean builds on an executable that the layered gate
    # makes unreachable; a qualification selection containing I13 still
    # builds Na+Ar before restoring the canonical Na-only tree.
    m2_executable_reachable = (
        args.campaign_mode == "qualification" and "I13" in selected and
        any(test_id in {"I18", "I19", "I20"} for test_id in selected)
    )

    build_result: dict[str, Any] | None = None
    m2_build_result: dict[str, Any] | None = None
    if not args.skip_build:
        if build_dir != DEFAULT_BUILD.resolve():
            parser.error("a clean build can target only the repository build directory")
        if m2_executable_reachable:
            m2_build_result = build_m2_qualification(
                output_root, args.stream, args.spice_path.resolve(),
                args.spice_kernels.resolve())
            if m2_build_result["status"] == "ERROR":
                summary = {
                    "source_revision": git_revision(),
                    "m2_build": m2_build_result,
                    "counts": {"executed": 0, "PASS": 0, "FAIL": 0,
                               "SKIPPED": 0, "ERROR": 0},
                    "campaign_status": "ERROR",
                    "campaign_error": "dedicated Na+Ar M2 build failed",
                    "results": [],
                }
                with (output_root / "run_all_summary.json").open(
                    "w", encoding="utf-8") as stream:
                    json.dump(summary, stream, indent=2, sort_keys=True)
                    stream.write("\n")
                print("ERROR: dedicated Na+Ar M2 build failed; see " +
                      str(m2_build_result["runner_log"]), file=sys.stderr)
                return 2
            m2_amps = Path(str(m2_build_result["executable"]))
        build_result = clean_build(
            output_root,
            args.stream,
            args.spice_path.resolve(),
            args.spice_kernels.resolve(),
        )
        if build_result["status"] == "ERROR":
            summary = {
                "source_revision": git_revision(),
                "build": build_result,
                "counts": {
                    "executed": 0,
                    "PASS": 0,
                    "FAIL": 0,
                    "SKIPPED": 0,
                    "ERROR": 0,
                },
                "campaign_status": "ERROR",
                "campaign_error": (
                    "Moon/SPICE configuration or clean build failed; "
                    "no test runner was executed"
                ),
                "results": [],
            }
            with (output_root / "run_all_summary.json").open(
                "w", encoding="utf-8"
            ) as stream:
                json.dump(summary, stream, indent=2, sort_keys=True)
                stream.write("\n")
            print(
                "ERROR: Moon/SPICE configuration or clean build failed; "
                f"see {build_result['log']}",
                file=sys.stderr,
            )
            return 2
    elif (m2_executable_reachable and
          (m2_amps is None or not m2_amps.is_file())):
        parser.error(
            "--skip-build with I18-I20 requires --m2-amps pointing to the "
            "dedicated Na+Ar executable")

    results: list[dict[str, Any]] = []
    for test_id in selected:
        result = run_test(
            test_id, registry[test_id], build_dir, output_root, args.stream,
            amps, m1_input, m2_input, args.campaign_mode, args.mpi_ranks,
            m2_amps,
        )
        results.append(result)
        print(f"{result['status']:7s} {test_id}  log={result['log']}")

    counts = {status: 0 for status in VALID_STATUSES}
    for result in results:
        counts[result["status"]] += 1
    counts["executed"] = len(results)
    print_summary(counts, results)

    # The aggregate scientific status is independent of the convenience exit
    # policy. A campaign with unresolved prerequisites remains SKIPPED even
    # when non-strict mode returns zero for developer workflows; it must never
    # be serialized as PASS merely because no evaluated criterion failed.
    campaign_status = "PASS"
    if counts["ERROR"]:
        campaign_status = "ERROR"
    elif counts["FAIL"] or (args.strict_skips and counts["SKIPPED"]):
        campaign_status = "FAIL"
    elif counts["SKIPPED"]:
        campaign_status = "SKIPPED"

    summary = {
        "source_revision": git_revision(),
        "build": build_result,
        "m2_build": m2_build_result,
        "build_reused": args.skip_build,
        "build_directory": str(build_dir),
        "linked_executable": str(amps),
        "linked_m2_na_ar_executable": str(m2_amps) if m2_amps else None,
        "linked_m1_input": str(m1_input),
        "linked_m2_input": str(m2_input),
        "campaign_mode": args.campaign_mode,
        "mpi_ranks": args.mpi_ranks,
        "output_directory": str(output_root),
        "duration_seconds": time.monotonic() - campaign_started,
        "counts": counts,
        "campaign_status": campaign_status,
        "strict_skips": args.strict_skips,
        "results": results,
    }
    summary_path = output_root / "run_all_summary.json"
    with summary_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(f"\nAggregate result: {summary_path}")

    if counts["ERROR"]:
        return 2
    if counts["FAIL"] or (args.strict_skips and counts["SKIPPED"]):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

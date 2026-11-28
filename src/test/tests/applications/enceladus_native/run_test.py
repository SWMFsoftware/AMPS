#!/usr/bin/env python3
"""Build Enceladus natively and import startup, case, MPI and mesh gates.

The adapter owns the checkout-wide build transition and follows the mandatory
preflight literally.  Scientific/configuration assertions remain in
``srcEnceladus/test/run_native_startup.py``; this file validates its exact
registry and translates the local receipt without upgrading any outcome.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import sys
import uuid
from pathlib import Path
from typing import List, Mapping, Tuple


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


LOCAL_CASES = [
    "N00_help", "N01_valid_default", "N02_stripe_surface_gate",
    "N03_plasma_cache_directory_gate", "N04_pointer_mover_generated",
    "N05_chemistry_processing_generated",
    "N06_case_local", "N07_case_mpi", "N08_case_transport",
    "N09_case_uvis", "N10_runtime_mesh", "N11_case_errors",
]
EXPECTED_CASES = [
    "ENCELADUS-N00-HELP", "ENCELADUS-N01-VALID-DEFAULT",
    "ENCELADUS-N02-STRIPE-SURFACE-GATE",
    "ENCELADUS-N03-PLASMA-CACHE-DIRECTORY-GATE",
    "ENCELADUS-N04-POINTER-MOVER-GENERATED",
    "ENCELADUS-N05-CHEMISTRY-PROCESSING-GENERATED",
    "ENCELADUS-N06-CASE-LOCAL", "ENCELADUS-N07-CASE-MPI",
    "ENCELADUS-N08-CASE-TRANSPORT", "ENCELADUS-N09-CASE-UVIS",
    "ENCELADUS-N10-RUNTIME-MESH", "ENCELADUS-N11-CASE-ERRORS",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def active_repo_processes(repo: Path) -> List[str]:
    """Return other native build/run processes rooted in this checkout."""

    proc = Path("/proc")
    if not proc.is_dir():
        return []
    active: List[str] = []
    native_names = {"make", "gmake", "amps", "Config.pl", "ampsConfig.pl",
                    "g++", "cc1plus", "mpicxx", "mpiexec", "mpirun"}
    for entry in proc.iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            cwd = (entry / "cwd").resolve()
            cwd.relative_to(repo)
            raw = (entry / "cmdline").read_bytes().replace(
                b"\x00", b" ").decode("utf-8", errors="replace").strip()
        except (OSError, ValueError):
            continue
        executable = Path(raw.split()[0]).name if raw else ""
        if (executable in native_names or executable.startswith("amps") or
                any(token in raw for token in ("Config.pl", "ampsConfig.pl"))):
            active.append("pid={} cwd={} command={}".format(entry.name, cwd, raw))
    return active


def build_native(context: Mapping, repo: Path, artifacts: Path,
                 timeout: int, model_data_root: Path = None) -> Tuple[Path, List[str]]:
    """Perform the exact Enceladus native preflight and retain its identity."""

    if Path.cwd().resolve() != repo or not (repo / "Config.pl").is_file():
        raise RunnerError("Enceladus native preflight is not at repository root")
    active = active_repo_processes(repo)
    if active:
        raise RunnerError("native build/run already active: {}".format("; ".join(active)))
    print("[preflight] repository={} cwd={} status=confirmed".format(
        repo, Path.cwd().resolve()), flush=True)
    print("[preflight] active_native_build_or_test_processes=0", flush=True)
    build = (repo / "build").resolve()
    if build != repo / "build":
        raise RunnerError("resolved build target is not repository/build")
    print("[preflight] generated_removal_target={} status=validated".format(build), flush=True)
    removed = run_command(["rm", "-rf", "--", "build"], repo, timeout=60)
    configure = [str(repo / "Config.pl"), "-application=enceladus"]
    if model_data_root is not None:
        configure.append("-model-data-path=" + str(model_data_root.resolve()))
    # Config copies the application tree into generated build storage. Old
    # application output can be many gigabytes and is never build input.
    # Hold it outside that tree on the SAME filesystem, and restore it even
    # when configuration fails. Nothing in the user's output is deleted.
    source_output = repo / "srcEnceladus" / "test_output"
    hold = repo / "test_output" / ("enceladus-build-hold-" + uuid.uuid4().hex)
    held = False
    try:
        if source_output.exists():
            hold.parent.mkdir(parents=True, exist_ok=True)
            source_output.rename(hold)
            held = True
        configured = run_command(configure, repo, timeout=timeout)
        if configured.returncode:
            raise RunnerError("Enceladus Config.pl failed")
        hooks = run_command([str(repo / "ampsConfig.pl"), "-input",
                             "enceladus.input", "-no-compile"], repo, timeout=timeout)
    finally:
        if held:
            hold.rename(source_output)
    make_argv = command(context, "make") + ["amps"]
    print("[preflight] build_environment=MAKEFLAGS=-j16 command={}".format(
        " ".join(make_argv)), flush=True)
    built = run_command(make_argv, repo, timeout=timeout,
                        env={"MAKEFLAGS": "-j16"})
    transcript = artifacts / "enceladus-native-build.txt"
    transcript.write_text(
        (removed.stdout or "") + (configured.stdout or "") +
        (hooks.stdout or "") + (built.stdout or ""),
        encoding="utf-8")
    if (removed.returncode != 0 or configured.returncode != 0 or
            hooks.returncode != 0 or built.returncode != 0 or not (repo / "amps").is_file()):
        raise RunnerError("native Enceladus configuration/build failed")
    retained = artifacts / "amps-enceladus"
    shutil.copy2(str(repo / "amps"), str(retained))
    identity = artifacts / "enceladus-build-identity.json"
    identity.write_text(json.dumps({
        "schema_version": 1,
        "application": "enceladus",
        "configure_command": configure,
        "production_hooks_command": ["./ampsConfig.pl", "-input", "enceladus.input", "-no-compile"],
        "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
        "executable_sha256": sha256(retained),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return retained, [str(transcript), str(identity)]


def load_receipt(path: Path) -> Mapping[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RunnerError("cannot read native startup receipt: {}".format(error))
    if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
        raise RunnerError("native startup receipt has no results array")
    return payload


def main() -> int:
    context = load_context("Run native Enceladus startup gates")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    if (not requested or len(set(requested)) != len(requested) or
            any(identifier not in EXPECTED_CASES for identifier in requested)):
        raise RunnerError("invalid Enceladus native case selection")
    repo = Path(str(context["repo_root"])).resolve()
    artifacts = Path(str(context["artifact_dir"]))
    artifacts.mkdir(parents=True, exist_ok=True)
    timeout = int(context["test"].get("timeout_seconds", 3600))
    lock_stream = (repo / "Config.pl").open("r")
    try:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            raise RunnerError("another native AMPS test holds Config.pl")
        # These fixtures contain declared synthetic rates/geometry, not fitted
        # Cassini physics. Keep the immutable publication beside the receipt
        # so a retained binary's compiled data-root remains inspectable.
        sys.path.insert(0, str(repo / "srcEnceladus" / "test"))
        from case_fixtures import make_fixture
        fixture_deck, data_root, _, _ = make_fixture(artifacts / "case-fixtures")
        retained, evidence = build_native(context, repo, artifacts, timeout, data_root)
        local_output = artifacts / "native-startup"
        for identifier in requested:
            print("[CASE] RUN  {} native startup gate".format(identifier), flush=True)
        completed = run_command(command(context, "python") + [
            str(repo / "srcEnceladus" / "test" / "run_native_startup.py"),
            "--executable", str(retained), "--repo-root", str(repo),
            "--output", str(local_output),
        ], repo, timeout=timeout, env={"PYTHONDONTWRITEBYTECODE": "1"})
        transcript = artifacts / "enceladus-native-startup.txt"
        transcript.write_text(completed.stdout or "", encoding="utf-8")
        receipt_path = local_output / "test_receipt.json"
        payload = load_receipt(receipt_path)
        rows = payload["results"]
        if completed.returncode != (0 if all(r.get("status") == "PASS" for r in rows) else 1):
            raise RunnerError("native startup exit/receipt disagreement")
        launcher = context["tools"].get("mpiexec", [])
        case_output = artifacts / "native-cases"
        case_completed = run_command(command(context, "python") + [
            str(repo / "srcEnceladus" / "test" / "run_case_native_tests.py"),
            "--executable", str(retained), "--fixture-deck", str(fixture_deck),
            "--output", str(case_output), "--launcher-json", json.dumps(launcher),
            "--rank-flag=" + str(context["execution"].get("mpi_rank_flag", "-n")),
        ], repo, timeout=timeout, env={"PYTHONDONTWRITEBYTECODE": "1"})
        case_transcript = artifacts / "enceladus-native-cases.txt"
        case_transcript.write_text(case_completed.stdout or "", encoding="utf-8")
        case_receipt = case_output / "test_receipt.json"
        case_rows = load_receipt(case_receipt)["results"]
        if case_completed.returncode != (1 if any(r.get("status") == "FAIL" for r in case_rows) else 0):
            raise RunnerError("native case exit/receipt disagreement")
        rows = rows + case_rows
        observed = [str(row.get("test")) for row in rows if isinstance(row, dict)]
        if observed != LOCAL_CASES or len(observed) != len(set(observed)):
            raise RunnerError("native startup registry mismatch: {}".format(observed))
        translated = {}
        for local_id, global_id, row in zip(LOCAL_CASES, EXPECTED_CASES, rows):
            status = str(row.get("status"))
            if status not in {"PASS", "FAIL", "SKIP"}:
                raise RunnerError("unsupported native startup status {!r}".format(status))
            translated[global_id] = case(
                global_id, status,
                " ".join(str(row.get("detail", "")).split()) or local_id)
        selected_rows = [translated[identifier] for identifier in requested]
        status = "FAIL" if any(row["status"] == "FAIL" for row in selected_rows) else "PASS"
        if all(row["status"] == "SKIP" for row in selected_rows):
            status = "SKIP"
        evidence.extend([str(transcript), str(receipt_path), str(case_transcript), str(case_receipt)])
        return write_result(
            context, status,
            "imported {} native Enceladus startup gates".format(len(selected_rows)),
            cases=selected_rows, evidence=evidence)
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)

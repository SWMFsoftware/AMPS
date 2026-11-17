#!/usr/bin/env python3
"""Run paired native checkpoint roles in genuinely separate MPI processes.

This is a software/API continuation gate on the documented fixed analytic
diagnostic scope. It does not add a production resume parser or qualify the
material-grid/four-polarization/relativistic-source physics. Each process owns
unique scratch, while only explicitly named checkpoint artifacts cross roles.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE / "validation"))
try:
    from sep_native_execution import create_native_work_directory
finally:
    sys.path.pop(0)


def quoted(value: str) -> str:
    """Match C++ std::quoted exactly, retaining UTF-8 and literal whitespace."""
    return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'


def binary_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_role(command: list[str], work: Path, log: Path, timeout: float) -> int:
    """Retain complete output and terminate the whole MPI job on timeout."""
    with log.open("w", encoding="utf-8") as stream:
        stream.write("COMMAND: " + json.dumps(command) + "\nWORK: " + str(work) + "\n")
        stream.flush()
        process = subprocess.Popen(command, cwd=work, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            return process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            stream.write("ERROR: complete MPI process group terminated after timeout\n")
            return 124


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpi-exec", default="mpiexec")
    parser.add_argument("--timeout", type=float, default=1200)
    parser.add_argument("--cases", nargs="+", choices=("CP05F01", "CP05F02", "CP05F03"),
                        default=("CP05F01", "CP05F02", "CP05F03"))
    args = parser.parse_args()
    amps, output = args.amps.resolve(), args.output.resolve()
    if not amps.is_file() or args.mpi_ranks < 1 or args.timeout <= 0:
        parser.error("a linked executable, positive rank count and timeout are required")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("choose a new or empty output directory; previous evidence is preserved")
    output.mkdir(parents=True, exist_ok=True)
    executable_hash = binary_hash(amps)
    results = []
    for case in dict.fromkeys(args.cases):
        case_dir = output / case
        case_dir.mkdir()
        prefix = case_dir / "physics"
        roles, failure = [], ""
        for role in ("write", "resume"):
            role_dir = case_dir / role
            role_dir.mkdir()
            work = create_native_work_directory(role_dir, role)
            request = role_dir / "fixture.args"
            # All nested paths are explicit absolute output paths. There is no
            # input-directory search, writable link to source, or env parser.
            request.write_text("SEP_CHECKPOINT_FRESH_PROCESS 1\nrole " + quoted(role) +
                               "\nprefix " + quoted(str(prefix)) + "\n", encoding="utf-8")
            report = role_dir / "srcsep-tests.json"
            command = [args.mpi_exec, "-n", str(args.mpi_ranks), str(amps),
                       "--test", case, "--test-input", str(request),
                       "--test-output-dir", str(role_dir / "artifacts"),
                       "--test-json", str(report), "--test-junit", str(role_dir / "srcsep-tests.xml")]
            start = time.monotonic()
            try:
                code = run_role(command, work, role_dir / "native.log", args.timeout)
                document = json.loads(report.read_text(encoding="utf-8")) if report.is_file() else {}
                matches = [row for row in document.get("results", []) if row.get("id") == case]
                passed = code == 0 and len(matches) == 1 and matches[0].get("status") == "PASS" and binary_hash(amps) == executable_hash
                if not passed:
                    failure = f"{role}: exit={code}, matching native PASS results={len(matches) if matches and matches[0].get('status') == 'PASS' else 0}"
            except (OSError, ValueError) as error:
                code, passed, failure = 2, False, f"{role}: {error}"
            roles.append({"role": role, "exit_code": code, "status": "PASS" if passed else "FAIL",
                          "duration_s": time.monotonic() - start, "report": str(report),
                          "log": str(role_dir / "native.log"), "working_directory": str(work)})
            if not passed:
                break
        expected, resumed = Path(str(prefix) + ".expected"), Path(str(prefix) + ".resumed")
        equal = not failure and expected.is_file() and resumed.is_file() and expected.read_bytes() == resumed.read_bytes()
        if not equal and not failure:
            failure = "final uninterrupted/resumed artifacts missing or differ byte-for-byte"
        results.append({"id": case, "status": "FAIL" if failure else "PASS", "message": failure,
                        "roles": roles, "final_artifact_equality": bool(equal)})
        print(case, results[-1]["status"], failure, flush=True)
    summary = {"schema": "srcsep-fresh-process-checkpoint-v1", "mpi_ranks": args.mpi_ranks,
               "executable": str(amps), "executable_sha256": executable_hash,
               "scope": "fixed-analytic-native-diagnostic-API;production-resume-CLI=false;output-continuation=false",
               "results": results}
    temporary = output / "summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, output / "summary.json")
    return 1 if any(row["status"] != "PASS" for row in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())

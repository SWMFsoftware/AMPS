#!/usr/bin/env python3
"""Run fresh native checkpoint roles with explicit material startup inputs.

The maintained CP05F IDs now use complete host steps on the approved moving
reduced background. --operator-only selects the separately linked frozen wave/
codec fixture, which never traces particles. Both paths preserve actual child
return codes, full reports and byte-exact final checkpoint comparison.
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
    parser.add_argument("--operator-only", action="store_true")
    parser.add_argument("--startup-input", type=Path,default=SOURCE/"examples/reduced-shock/field_line.in")
    parser.add_argument("--event", type=Path,default=SOURCE/"examples/reduced-shock/positive_1au.event")
    parser.add_argument("--cases", nargs="+",choices=("CP05F01","CP05F02","CP05F03","CP05OF01","CP05OF02","CP05OF03"))
    args = parser.parse_args()
    amps, output = args.amps.resolve(), args.output.resolve()
    if not amps.is_file() or args.mpi_ranks < 1 or args.timeout <= 0:
        parser.error("a linked executable, positive rank count and timeout are required")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("choose a new or empty output directory; previous evidence is preserved")
    output.mkdir(parents=True, exist_ok=True)
    executable_hash = binary_hash(amps)
    cases=args.cases or (["CP05OF01","CP05OF02","CP05OF03"] if args.operator_only else ["CP05F01","CP05F02","CP05F03"])
    if any(identifier.startswith("CP05OF")!=args.operator_only for identifier in cases):
        parser.error("case IDs and --operator-only scope disagree")
    results = []
    for case in dict.fromkeys(cases):
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
            version=1 if args.operator_only else 2
            text="SEP_CHECKPOINT_FRESH_PROCESS "+str(version)+"\nrole "+quoted(role)+"\nprefix "+quoted(str(prefix))+"\n"
            if not args.operator_only:
                text+="startup_input "+quoted(str(args.startup_input.resolve(strict=True)))+"\nevent "+quoted(str(args.event.resolve(strict=True)))+"\n"
            request.write_text(text,encoding="utf-8")
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
               "scope": "isolated-frozen-wave-codec;particle_motion=false" if args.operator_only else
                        "prescribed-reduced-material-native-host;physical_dt_unchanged=true",
               "results": results}
    temporary = output / "summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, output / "summary.json")
    return 1 if any(row["status"] != "PASS" for row in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())

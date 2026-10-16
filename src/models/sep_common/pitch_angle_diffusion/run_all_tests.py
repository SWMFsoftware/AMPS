#!/usr/bin/env python3
"""Run every test of the pitch-angle diffusion library and summarize them.

Stages (in order; each prints '[ID] STATUS message' lines):

  build        make: library, unit tests, reference driver, acceptance binary
  PAD01-PAD12  deterministic unit tests (test_pitch_angle_diffusion.cpp)
  PAR01-PAR11  independent reference checks (reference_checks.py: mpmath,
               companion bundles, supplement check_physics.py, registry
               digest, interface contract)
  AP-*         statistical acceptance campaign (test_pitch_angle_montecarlo
               .cpp; roadmap Section 26.5), 52 cases, ~2.3 h on 16 idle
               cores (3.9 h measured on a shared node); each case's
               result is printed as soon as it finishes

Acceptance evidence reuse.  The campaign is deterministic: every batch has
its own seeded random stream, so a run of the same binary reproduces the same
numbers bit for bit, independently of the thread count.  The runner
therefore reuses a recorded run (build/acceptance/results.txt) only when the
SHA-256 of the current acceptance binary equals the digest recorded with it
(build/acceptance/binary.sha256); any change to the kernels, coefficients,
quadrature or the campaign itself changes the binary and forces a rerun.
The certification matrix (pitch_angle_certification.cpp) is not linked into
that binary, so recording evidence there does not invalidate it.

Options:
  --rerun-acceptance   ignore a recorded run and execute the campaign
  --skip-acceptance    do not run or reuse it (reported as SKIP)
  --quick-acceptance   1/50 particle counts; reported as SKIP (not evidence)
  --threads N          campaign threads (default 16)
  --json PATH          machine-readable summary

Exit status: 0 when nothing failed, 1 otherwise.  Usage (any directory):
  python3 src/models/sep_common/pitch_angle_diffusion/run_all_tests.py
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
ACCEPTANCE_DIR = BUILD / "acceptance"
LINE = re.compile(r"^\[([^\]]+)\]\s+(\S+)\s*(.*)$")

records: list[dict] = []


def emit(test_id: str, status: str, message: str, stage: str) -> None:
    records.append({"id": test_id, "status": status, "message": message, "stage": stage})
    print(f"[{test_id}] {status} {message}", flush=True)


def parse_lines(text: str, stage: str, quick: bool = False) -> None:
    for line in text.splitlines():
        match = LINE.match(line.strip())
        if not match:
            continue
        test_id, status, message = match.groups()
        if status.startswith("QUICK-"):
            emit(test_id, "SKIP", f"quick mode (not acceptance evidence): {status[6:]} "
                                  f"{message}", stage)
        elif status in ("PASS", "FAIL", "SKIP"):
            emit(test_id, status, message, stage)
        else:
            emit(test_id, "FAIL", f"unrecognized status {status}: {message}", stage)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=HERE, capture_output=True, text=True, **kwargs)


def stage_build() -> bool:
    started = time.time()
    result = run(["make", "all", "tests"])
    ok = result.returncode == 0
    emit("build", "PASS" if ok else "FAIL",
         f"make all tests ({time.time() - started:.0f} s)" if ok
         else (result.stdout + result.stderr)[-1500:], "build")
    return ok


def stage_unit() -> None:
    result = run([str(BUILD / "test_pitch_angle_diffusion")])
    parse_lines(result.stdout, "unit")
    if "pitch-angle unit summary" not in result.stdout:
        emit("PAD", "FAIL", f"unit test binary exited {result.returncode} without a "
                            f"summary: {result.stderr[-500:]}", "unit")


def stage_reference() -> None:
    result = run([sys.executable, str(HERE / "reference_checks.py"),
                  "--driver", str(BUILD / "pitch_angle_driver")])
    parse_lines(result.stdout, "reference")
    if "pitch-angle reference summary" not in result.stdout:
        emit("PAR", "FAIL", f"reference checks exited {result.returncode}: "
                            f"{(result.stdout + result.stderr)[-500:]}", "reference")


def stage_acceptance(args) -> None:
    binary = BUILD / "test_pitch_angle_montecarlo"
    digest = sha256(binary)
    results = ACCEPTANCE_DIR / "results.txt"
    recorded = ACCEPTANCE_DIR / "binary.sha256"
    if args.skip_acceptance:
        emit("AP", "SKIP", "acceptance campaign skipped (--skip-acceptance)", "acceptance")
        return
    if args.quick_acceptance:
        result = run([str(binary), "--quick", "--threads", str(args.threads)])
        parse_lines(result.stdout, "acceptance")
        return
    reusable = (not args.rerun_acceptance and results.exists() and recorded.exists() and
                recorded.read_text().split()[0] == digest and
                "pitch-angle acceptance summary" in results.read_text())
    if reusable:
        stamp = datetime.datetime.fromtimestamp(results.stat().st_mtime).isoformat(
            timespec="seconds")
        print(f"acceptance: reusing the run recorded {stamp} for binary sha256 {digest[:16]}... "
              f"(identical binary; deterministic seeds)", flush=True)
        parse_lines(results.read_text(), "acceptance")
        return
    ACCEPTANCE_DIR.mkdir(parents=True, exist_ok=True)
    recorded.write_text(f"{digest}  build/test_pitch_angle_montecarlo\n")
    print(f"acceptance: running the campaign (binary sha256 {digest[:16]}...; ~2.3 h on "
          f"{args.threads} idle cores, longer on a shared node; results.txt is written "
          f"case by case)", flush=True)
    started = time.time()
    # Stream the campaign: every case line is written to results.txt and
    # judged as soon as the binary prints it, so progress is visible while
    # the campaign runs.  A file without the final summary line (interrupted
    # run) is never reused.
    output: list[str] = []
    with results.open("w") as record, subprocess.Popen(
            [str(binary), "--threads", str(args.threads)], cwd=HERE,
            stdout=subprocess.PIPE, text=True, bufsize=1) as process:
        for line in process.stdout:
            output.append(line)
            record.write(line)
            record.flush()
            parse_lines(line, "acceptance")
        code = process.wait()
        record.write(f"real {time.time() - started:.0f} s\nEXIT={code}\n")
    if "pitch-angle acceptance summary" not in "".join(output):
        emit("AP", "FAIL", f"acceptance binary exited {code} without a summary",
             "acceptance")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rerun-acceptance", action="store_true")
    parser.add_argument("--skip-acceptance", action="store_true")
    parser.add_argument("--quick-acceptance", action="store_true")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if stage_build():
        stage_unit()
        stage_reference()
        stage_acceptance(args)
    counts = {status: sum(r["status"] == status for r in records)
              for status in ("PASS", "FAIL", "SKIP")}
    print(f"pitch-angle library summary: {counts['PASS']} passed, {counts['FAIL']} failed, "
          f"{counts['SKIP']} skipped")
    failed = [r for r in records if r["status"] == "FAIL"]
    for r in failed:
        print(f"  FAILED [{r['id']}] ({r['stage']}): {r['message'][:200]}")
    if args.json:
        args.json.write_text(json.dumps({"counts": counts, "records": records}, indent=1))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Top-level orchestration for the December 2006 cutoff study.

Stages are intentionally explicit.  A failed validation command returns a
nonzero status and normally stops later scientific interpretation.  Use
``--continue-on-validation-failure`` only for debugging or to inspect partial
products; the final archive records that choice.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

from study_common import default_output_root, load_config, resolve_output_path


STAGE_ORDER = ("validate", "pamela", "poes", "morphology", "compare", "dynamics", "figures")


def execute(
    command: Sequence[str],
    cwd: Path,
    log_path: Path,
    dry_run: bool,
    *,
    stage: str = "command",
    stage_index: int = 1,
    stage_count: int = 1,
) -> int:
    """Execute one study stage while teeing its complete output to the terminal.

    The original orchestrator redirected the child process exclusively to its
    log file.  That was safe for archiving, but a long AMPS calculation then
    appeared completely idle.  The child now writes to a pipe which is copied
    immediately to both destinations.  ``PYTHONUNBUFFERED`` is set for Python
    stage runners so their own progress and AMPS pass-through output arrive as
    soon as they are produced rather than at process exit.

    The function intentionally returns the child's unmodified exit code.  This
    preserves the existing stop/continue policy and makes shell, MPI, and AMPS
    failures visible in both the manifest and the terminal summary.
    """

    command_text = " ".join(shlex.quote(token) for token in command)
    label = stage.upper()
    started = time.monotonic()
    print("=" * 78, flush=True)
    print(f"[{stage_index}/{stage_count}] START {label}", flush=True)
    print(f"[{label}] working directory: {cwd}", flush=True)
    print(f"[{label}] log: {log_path}", flush=True)
    print(f"[{label}] command: {command_text}", flush=True)
    if dry_run:
        print(f"[{stage_index}/{stage_count}] PREPARED {label} (not executed)", flush=True)
        return 0

    log_path.parent.mkdir(parents=True, exist_ok=True)
    child_environment = os.environ.copy()
    child_environment["PYTHONUNBUFFERED"] = "1"
    with log_path.open("w", encoding="utf-8", buffering=1) as log:
        log.write(f"Stage: {stage}\nWorking directory: {cwd}\nCommand: {command_text}\n\n")
        process = subprocess.Popen(
            list(command), cwd=cwd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1,
            env=child_environment,
        )
        assert process.stdout is not None
        with process.stdout:
            for line in process.stdout:
                # Prefixing every line identifies the responsible stage when
                # AMPS, MPI, and postprocessing messages are interleaved in
                # batch output.
                sys.stdout.write(f"[{label}] {line}")
                sys.stdout.flush()
                log.write(line)
        return_code = process.wait()

    elapsed = time.monotonic() - started
    outcome = "PASS" if return_code == 0 else "FAIL"
    stream = sys.stdout if return_code == 0 else sys.stderr
    print(
        f"[{stage_index}/{stage_count}] {outcome} {label}: "
        f"exit={return_code}, elapsed={elapsed:.1f} s, log={log_path}",
        file=stream, flush=True,
    )
    return return_code


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", action="append", choices=STAGE_ORDER + ("all",),
                        help="May be repeated; default is all stages")
    parser.add_argument("--profile", choices=("SMOKE", "ROUTINE", "FULL"), default="SMOKE")
    parser.add_argument("--amps", type=Path, default=Path("./amps"))
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("-np", type=int)
    parser.add_argument("-nt", type=int)
    parser.add_argument("--prepare-only", action="store_true",
                        help="Validate data and generate commands/inputs without running AMPS")
    parser.add_argument("--continue-on-validation-failure", action="store_true")
    parser.add_argument(
        "--output-root", type=Path, default=default_output_root(),
        help=("Study output root (default: repository-level "
              "test_output/dec2006_ts05_cutoff_erosion)"),
    )
    return parser.parse_args()


def main() -> int:
    run_started = time.monotonic()
    args = parse_args()
    root, config = load_config()
    execution = config["execution"]
    np_value = args.np if args.np is not None else int(execution["mpi_ranks"])
    nt_value = args.nt if args.nt is not None else int(execution["threads_per_rank"])
    output = resolve_output_path(args.output_root)
    output.mkdir(parents=True, exist_ok=True)
    amps_path = args.amps.expanduser()
    if not amps_path.is_absolute():
        amps_path = (Path.cwd() / amps_path).resolve()
    stages = list(args.stage or ["all"])
    if "all" in stages:
        stages = list(STAGE_ORDER)
    stages = [stage for stage in STAGE_ORDER if stage in stages]

    # Fail before starting the validation pipeline when a requested modeling
    # stage cannot possibly launch.  This catches common relative-path errors
    # (for example ``--amps ../amps`` from the repository root) without first
    # creating dozens of case directories or relying on an opaque MPI status.
    amps_stages = {"pamela", "poes", "morphology"}
    needs_amps = bool(amps_stages.intersection(stages)) and not args.prepare_only
    if needs_amps:
        if not amps_path.is_file():
            print(f"ERROR: AMPS executable does not exist: {amps_path}", file=sys.stderr)
            print("Pass --amps ./amps when launching from the AMPS repository root.",
                  file=sys.stderr)
            return 2
        if not os.access(amps_path, os.X_OK):
            print(f"ERROR: AMPS file is not executable: {amps_path}", file=sys.stderr)
            return 2
        if shutil.which(args.mpirun) is None:
            print(f"ERROR: MPI launcher is not available: {args.mpirun}", file=sys.stderr)
            return 2

    python = sys.executable
    commands: Dict[str, List[str]] = {
        "validate": [python, str(root / "scripts" / "validate_package.py")],
        "pamela": [
            python, "run_C9.py", "--profile", args.profile,
            "--solver", "GRIDDED", "--cutoff-evaluation", "DIRECT_ACCESS",
            "--comparison-observable", "PAMELA_T50", "--max-trace-time", "300",
            "--output-root", str(output / "C9"), "--amps", str(amps_path),
            "--mpirun", args.mpirun, "-np", str(np_value), "-nt", str(nt_value),
            "--mode3d-parallel-field-init",
        ],
        "poes": [
            python, "run_C10.py", "--profile", args.profile,
            "--solver", "GRIDDED", "--cutoff-evaluation", "DIRECT_ACCESS",
            "--comparison-observable", "ACCESS_T50", "--max-trace-time", "300",
            "--output-root", str(output / "C10"), "--amps", str(amps_path),
            "--mpirun", args.mpirun, "-np", str(np_value), "-nt", str(nt_value),
            "--mode3d-parallel-field-init",
        ],
        "morphology": [
            python, str(root / "scripts" / "run_morphology.py"),
            "--profile", args.profile, "--amps", str(amps_path),
            "--mpirun", args.mpirun, "-np", str(np_value), "-nt", str(nt_value),
            "--output-root", str(output / "morphology"),
        ],
        "compare": [
            python, str(root / "scripts" / "compare_observations.py"),
            "--c9-root", str(output / "C9" / "gridded"),
            "--c10-root", str(output / "C10" / "gridded"),
            "--output-root", str(output / "comparison"),
        ],
        "dynamics": [
            python, str(root / "scripts" / "analyze_dynamics.py"),
            "--morphology-root", str(output / "morphology"),
            "--output-root", str(output / "dynamics"),
        ],
        "figures": [
            python, str(root / "scripts" / "make_figures.py"),
            "--comparison-root", str(output / "comparison"),
            "--dynamics-root", str(output / "dynamics"),
            "--output-root", str(output / "figures"),
        ],
    }
    if args.prepare_only:
        commands["pamela"].append("--dry-run")
        commands["poes"].append("--dry-run")
        commands["morphology"].append("--prepare-only")

    cwd = {
        "validate": root,
        "pamela": root / "vendor" / "C9",
        "poes": root / "vendor" / "C10",
        "morphology": root,
        "compare": root,
        "dynamics": root,
        "figures": root,
    }
    record = {
        "study_id": config["study_id"], "profile": args.profile,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "prepare_only": args.prepare_only, "stages": stages,
        "mpi_ranks": np_value, "threads_per_rank": nt_value,
        "commands": {stage: commands[stage] for stage in stages},
        "return_codes": {},
        "stage_elapsed_seconds": {},
    }
    print("=" * 78, flush=True)
    print("December 2006 TS05 cutoff-erosion study", flush=True)
    print(f"profile: {args.profile}", flush=True)
    print(f"stages: {', '.join(stages)}", flush=True)
    print(f"MPI ranks / threads per rank: {np_value} / {nt_value}", flush=True)
    print(f"AMPS executable: {amps_path}", flush=True)
    print(f"output root: {output}", flush=True)
    print(f"prepare only: {args.prepare_only}", flush=True)

    final_rc = 0
    for stage_index, stage in enumerate(stages, start=1):
        # Comparison and inference cannot exist in preparation-only mode because
        # no model output was generated.  Their commands remain in the record if
        # explicitly requested, but are skipped rather than failing on absence.
        if args.prepare_only and stage in ("compare", "dynamics", "figures"):
            print(f"[{stage_index}/{len(stages)}] SKIP {stage.upper()}: --prepare-only",
                  flush=True)
            record["return_codes"][stage] = None
            record["stage_elapsed_seconds"][stage] = 0.0
            continue
        stage_started = time.monotonic()
        rc = execute(
            commands[stage], cwd[stage], output / "logs" / f"{stage}.log", False,
            stage=stage, stage_index=stage_index, stage_count=len(stages),
        )
        record["return_codes"][stage] = rc
        record["stage_elapsed_seconds"][stage] = round(
            time.monotonic() - stage_started, 3
        )
        final_rc = final_rc or rc
        if rc and not args.continue_on_validation_failure:
            print(
                f"Stopping after {stage.upper()} failure. "
                "Use --continue-on-validation-failure only for diagnostic runs.",
                file=sys.stderr, flush=True,
            )
            break
    record["passed"] = final_rc == 0
    record["elapsed_seconds"] = round(time.monotonic() - run_started, 3)
    (output / "study_run_manifest.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print("=" * 78, flush=True)
    print(
        f"STUDY {'PASS' if final_rc == 0 else 'FAIL'}: exit={final_rc}, "
        f"elapsed={record['elapsed_seconds']:.1f} s",
        flush=True,
    )
    for stage in stages:
        if stage not in record["return_codes"]:
            status = "NOT RUN"
        elif record["return_codes"][stage] is None:
            status = "SKIPPED"
        elif record["return_codes"][stage] == 0:
            status = "PASS"
        else:
            status = f"FAIL ({record['return_codes'][stage]})"
        elapsed = record["stage_elapsed_seconds"].get(stage)
        elapsed_text = "" if elapsed is None else f", {elapsed:.1f} s"
        print(f"  {stage:<10} {status}{elapsed_text}", flush=True)
    print(f"Run manifest: {output / 'study_run_manifest.json'}", flush=True)
    return final_rc


if __name__ == "__main__":
    raise SystemExit(main())

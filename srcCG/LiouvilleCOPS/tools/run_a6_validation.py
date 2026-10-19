#!/usr/bin/env python3
"""Run the Stage-A6 validation gate and optionally write a JSON record.

The default gate runs the complete fast Python suite (integration tests remain
explicitly skipped) and the C++ analytic/transport self-test.  ``--integration``
also enables the temporary-directory SHAP5 tests, including the zero-
enhancement regression and fitted Python/C++ pressure replay.

This driver does not fit the full flight dataset and does not turn same-model
synthetic recovery into an independent physical validation.  Its purpose is
to execute the implementation-consistency gate defined in
``amendement_development_plan.md`` in one reproducible command.
"""

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from cops_paths import prepare_output


APP = Path(__file__).resolve().parents[1]


def sha256(path):
    """Stream a file identity without loading a large executable into RAM."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_step(name, command, environment):
    """Run one gate step, echoing output while retaining it for the report."""
    print("\n=== %s ===" % name, flush=True)
    print("command: " + " ".join(command), flush=True)
    started = time.monotonic()
    process = subprocess.Popen(
        command, cwd=str(APP), env=environment, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        bufsize=1)
    lines = []
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="", flush=True)
        lines.append(line)
    returncode = process.wait()
    output = "".join(lines)
    elapsed = time.monotonic() - started
    match = re.search(r"Ran (\d+) tests? in", output)
    skipped = re.search(r"OK \(skipped=(\d+)\)", output)
    return {
        "name": name,
        "command": command,
        "returncode": returncode,
        "elapsed_seconds": elapsed,
        "test_count": int(match.group(1)) if match else None,
        "skipped_count": int(skipped.group(1)) if skipped else 0,
        "output_sha256": hashlib.sha256(output.encode("utf-8")).hexdigest(),
        # Preserve only the tail: it contains unittest/selftest summaries but
        # avoids turning the JSON record into a duplicate multi-megabyte log.
        "output_tail": "\n".join(output.splitlines()[-80:]),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--integration", action="store_true",
        help="enable local solver/mission-data integration tests")
    parser.add_argument(
        "--threads", type=int, default=2,
        help="OpenMP threads for C++ integration and self-tests (default: 2)")
    parser.add_argument(
        "--skip-selftest", action="store_true",
        help="run only the Python unittest gate")
    parser.add_argument(
        "--report", default=None,
        help="optional new JSON report path under cops_output, or an absolute path outside the application")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")

    # Resolve and validate the report name before starting the potentially
    # expensive integration suite.  This deliberately checks only existence:
    # prepare_output() creates the parent directory but never the report file,
    # so a failed validation can still leave a complete diagnostic record.
    report_target = None
    if args.report:
        report_target = Path(prepare_output(args.report))
        if report_target.exists():
            parser.error("report already exists: %s" % report_target)

    executable = APP / "liouville_cops"
    if not executable.is_file():
        parser.error("missing executable %s; compile before validation" % executable)

    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = str(args.threads)
    # The test suite controls process-level concurrency.  One BLAS thread per
    # test prevents accidental multiplication of OpenMP and BLAS workers.
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    if args.integration:
        environment["COPS_RUN_INTEGRATION"] = "1"
    else:
        environment.pop("COPS_RUN_INTEGRATION", None)

    steps = [run_step(
        "Python unit/synthetic%s suite" %
        ("/integration" if args.integration else ""),
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
        environment)]
    if not args.skip_selftest:
        steps.append(run_step(
            "C++ analytic and transport self-test",
            [str(executable), "selftest.in"], environment))

    status = "passed" if all(step["returncode"] == 0 for step in steps) else "failed"
    report = {
        "schema": "liouville_cops.a6_validation",
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": status,
        "integration_enabled": bool(args.integration),
        "threads": args.threads,
        "environment": {name: environment[name] for name in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                         "MKL_NUM_THREADS")},
        "software": {
            "python": sys.version.split()[0],
            "executable": {
                "path": str(executable.resolve()),
                "size_bytes": executable.stat().st_size,
                "sha256": sha256(executable),
            },
            "a6_test": {
                "path": str((APP / "tests/test_a6_validation.py").resolve()),
                "sha256": sha256(APP / "tests/test_a6_validation.py"),
            },
        },
        "scope": {
            "claim": "implementation consistency before flight-data experiments",
            "not_claimed": [
                "independent validation of coma physics",
                "instrument calibration validation",
                "flight-data source localization",
                "optimizer convergence for the full descent",
            ],
        },
        "steps": steps,
    }
    if report_target is not None:
        report_target.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        print("\nA6 report: %s" % report_target)
    print("\nStage A6: %s" % status.upper())
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Run the layered srcMoon stand-alone verification suite.

This is the local U-series campaign driver.  It does not run or summarize the
I-series linked/convergence/validation campaign, and a later test never masks
an earlier FAIL or ERROR.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from common.moon_testlib import (
    DEFAULT_BUILD,
    DEFAULT_OUTPUT,
    REPO_ROOT,
    load_acceptance,
    run_one,
    test_directories,
)


def rebuild(output_root: Path) -> int:
    """Use the repository's required clean AMPS build sequence.

    Removing build/ is intentional because Config.pl copies and rewrites a
    generated source tree.  Reusing that tree after source/configuration
    changes can otherwise test stale code.  The complete build stream is kept
    as an artifact rather than mixed into the JSON science result.
    """
    build_dir = REPO_ROOT / "build"
    if build_dir.exists():
        shutil.rmtree(build_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    log_path = output_root / "build.log"
    with log_path.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(
            ["make", "-j"],
            cwd=REPO_ROOT,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return completed.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--list", action="store_true")
    selection.add_argument("--all", action="store_true")
    selection.add_argument("--test", action="append", metavar="UXX")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--strict-skips", action="store_true")
    args = parser.parse_args()

    directories = test_directories()
    if args.list:
        for test_id in directories:
            acceptance = load_acceptance(test_id)
            print(f"{test_id}  {acceptance['name']}")
        return 0

    output_root = args.output_dir.resolve()
    if args.rebuild and rebuild(output_root) != 0:
        print(f"ERROR clean build failed; see {output_root / 'build.log'}", file=sys.stderr)
        return 2

    selected = list(directories) if args.all else args.test
    unknown = [test_id for test_id in selected if test_id not in directories]
    if unknown:
        parser.error(f"unknown test ID(s): {', '.join(unknown)}")

    # Execute in frozen ID order.  These tests are intentionally serial here:
    # several modes share a cached probe executable and deterministic output
    # paths, while campaign-level parallelism offers little benefit.
    results = [
        run_one(test_id, args.build_dir.resolve(), output_root)
        for test_id in selected
    ]
    counts = {status: 0 for status in ("PASS", "FAIL", "ERROR", "SKIPPED")}
    for result in results:
        counts[result["status"]] += 1
        print(f"{result['status']:7s} {result['id']} {result['name']}")

    summary = {"counts": counts, "results": results}
    output_root.mkdir(parents=True, exist_ok=True)
    with (output_root / "summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")

    # Infrastructure errors outrank evaluated failures.  SKIPPED is normally a
    # recorded incomplete capability; CI/release gates can opt into rejecting
    # it with --strict-skips.
    if counts["ERROR"]:
        return 2
    if counts["FAIL"]:
        return 1
    if args.strict_skips and counts["SKIPPED"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

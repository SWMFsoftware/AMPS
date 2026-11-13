#!/usr/bin/env python3
"""Prepare a review-only I04 candidate from two reproducible linked runs.

This script never installs/blesses a reference or changes an acceptance hash.
A candidate is a regression proposal, not an independent physical oracle.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def values(path):
    """Read the named nightly microscopic state, preserving its original file."""
    for line in path.read_text().splitlines():
        if line.startswith("AMPS: 0\t"):
            row = line.split()
            return dict(zip(("mean_speed_m_s", "mean_vx_m_s", "mean_vy_m_s", "mean_vz_m_s"),
                            map(float, row[2:6])))
    raise ValueError(f"missing neutral-Na microscopic output: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-run", type=Path, required=True)
    parser.add_argument("--second-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--budget-result", type=Path, required=True)
    parser.add_argument("--previous-reference", type=Path,
                        default=HERE/"reference/test_Moon.ref_np4",
                        help="use the archived prior reference to reproduce a historical review after approval")
    # Preserve reproducibility of the original chemistry review while letting
    # a later, independently justified core integration name its actual cause.
    # This metadata never changes output bytes or grants reference approval.
    parser.add_argument("--reason",default=(
        "corrected the configured neutral-Na photochemistry dispatch; physical lifetime unchanged"),
        help="documented mechanism requiring this particular reference review")
    args = parser.parse_args()
    if not args.reason.strip():
        parser.error("--reason must state the documented change")
    first, second = args.first_run.resolve(), args.second_run.resolve()
    for run in (first, second):
        if len(list((run/"PT").glob("thread=*.log"))) != 4:
            raise ValueError(f"reviewed nightly layout requires four ranks: {run}")
    files = [run/"PT/plots/test_Moon.dat" for run in (first, second)]
    if files[0].read_bytes() != files[1].read_bytes() or sha(first/"amps") != sha(second/"amps"):
        raise ValueError("candidate requires identical output from the same executable in two runs")
    budget = json.loads(args.budget_result.read_text())
    if budget["status"] != "PASS" or not budget["conditions"]["chemistry_survival"]:
        raise ValueError("candidate requires independent linked chemistry/ledger fixture PASS")
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite reference-review evidence: {output}")
    output.mkdir(parents=True, exist_ok=True)
    candidate = output/"test_Moon.ref_np4"
    shutil.copy2(files[0], candidate)
    frozen = args.previous_reference.resolve()
    report = {"status": "AWAITING REVIEW", "installed": False,
        "reason": args.reason.strip(),
        "generator": {"path": str(Path(__file__).resolve()),
                      "sha256": sha(Path(__file__))},
        "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "source_is_dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=REPO, text=True).strip()),
        "frozen_reference_sha256": sha(frozen),
        "candidate_sha256": sha(candidate), "nightly_executable_sha256": sha(first/"amps"),
        "first_run": str(first), "second_run": str(second),
        "budget_result": {"path": str(args.budget_result.resolve()), "sha256": sha(args.budget_result)},
        "previous_values": values(frozen), "candidate_values": values(candidate),
        "input_sha256": sha(REPO/"input/moon.input"),
        "note": "No tolerance, archived reference, shared SWMF_data file or I13 baseline was changed."}
    (output/"reference_review.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

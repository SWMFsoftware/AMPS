#!/usr/bin/env python3
"""Compare the nightly Moon output with I04's already-reviewed reference.

The checker pins the user-reviewed reference and its provenance. The external
SWMF_data tree is shared between AMPS checkouts; its HEAD cannot supersede
that contract. Historical import recovers only the superseded Git reference.
This checker never derives or approves a reference from current model output,
and the shared tree is read-only throughout.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "reference/test_Moon.ref_np4"
PROVENANCE = HERE / "reference/provenance.json"
# This publication of a previously reviewed regression reference is itself
# traceable to the source-data repository object. Scientific tolerances remain
# the normal DiffNum defaults and the exact checksum in acceptance.json.
REVISION = "7a1c437369a600bd248bcf8846b66b286e7109c1"
SOURCE_PATH = "PT/AMPS/data/output/test_Moon/test_Moon.ref_np4"
GIT_REFERENCE_HASH = "abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4"


def expected_hash() -> str:
    """Use the existing frozen I04 acceptance contract as the authority."""
    return json.loads((HERE / "reference/acceptance.json").read_text())["data_hashes"][
        "output/test_Moon/test_Moon.ref_np4"]


def import_reviewed(repository: Path) -> None:
    """Reproduce the exact reference bytes without touching external files."""
    raw = subprocess.run(["git", "-C", str(repository), "show",
                          f"{REVISION}:{SOURCE_PATH}"], check=True,
                         stdout=subprocess.PIPE).stdout
    if hashlib.sha256(raw).hexdigest() != GIT_REFERENCE_HASH:
        raise ValueError("authoritative historical Git reference checksum mismatch")
    # The corrected Na reference was reviewed in AMPS, not published at this
    # historical SWMF_data object. Reproducing old evidence must preserve the
    # new reference and put the old bytes in their explicitly superseded path.
    destination = REFERENCE if expected_hash() == GIT_REFERENCE_HASH else (
        HERE / "reference/superseded/test_Moon.ref_np4.pre-na-photochemistry")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.read_bytes() != raw:
        raise ValueError("refusing to overwrite a different existing reviewed reference")
    destination.write_bytes(raw)


def check(run_directory: Path, diff: Path, diffnum: Path) -> int:
    """Fail closed on missing/corrupt reference evidence or wrong MPI layout."""
    evidence = {"id": "I04", "reference_path": str(REFERENCE),
                "expected_sha256": expected_hash(), "status": "ERROR"}
    try:
        actual_hash = hashlib.sha256(REFERENCE.read_bytes()).hexdigest()
        evidence["actual_sha256"] = actual_hash
        if actual_hash != expected_hash():
            raise ValueError("reviewed reference checksum mismatch")
        hashes = json.loads((HERE / "reference/acceptance.json").read_text())["data_hashes"]
        provenance_hash = hashlib.sha256(PROVENANCE.read_bytes()).hexdigest()
        if provenance_hash != hashes["reference/provenance.json"]:
            raise ValueError("reviewed reference provenance checksum mismatch")
        provenance = json.loads(PROVENANCE.read_text())
        if provenance.get("status") != "REVIEWED" or provenance.get("reference_sha256") != actual_hash:
            raise ValueError("reference provenance is unreviewed or identifies different bytes")
        evidence["provenance_sha256"] = provenance_hash
        ranks = len(list((run_directory / "PT").glob("thread=*.log")))
        if ranks != 4:
            raise ValueError(f"the reviewed I04 reference requires four MPI ranks, got {ranks}")
        result = run_directory / "PT/plots/test_Moon.dat"
        if not result.is_file():
            raise ValueError(f"production regression output is missing: {result}")
        command = [str(diffnum.resolve()),
                   str(result.resolve()),
                   str(REFERENCE)]
        compared = subprocess.run(command, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, check=False)
        # DiffNum uses Perl die() for both evaluated disagreement and I/O
        # errors, so its nonzero exit status alone cannot distinguish FAIL
        # from ERROR. Only its completed difference-count summaries establish
        # an evaluated comparison. Preserve every diagnostic in either case.
        diagnostics = compared.stdout + compared.stderr
        diff.write_bytes(diagnostics)
        completed_difference = re.search(
            rb"ERROR in DiffNum\.pl: number of (?:significant numerical differences|differences in the text) is [1-9][0-9]*",
            diagnostics) is not None
        if compared.returncode < 0 or (compared.returncode != 0 and not completed_difference):
            status = "ERROR"
        else:
            status = "PASS" if compared.returncode == 0 and not diagnostics else "FAIL"
        evidence.update(command=command,
                        status=status, comparator_exit_code=compared.returncode,
                        diff_bytes=diff.stat().st_size)
    except (OSError, ValueError) as error:
        diff.write_text(f"ERROR I04: {error}\n")
        evidence["reason"] = str(error)
    diff.with_suffix(".reference.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(f"{evidence['status']} I04: reviewed reference {expected_hash()}")
    return {"PASS": 0, "FAIL": 1, "ERROR": 2}[evidence["status"]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--import-from-repository", type=Path)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--diff-file", type=Path, default=Path("test_Moon.diff"))
    parser.add_argument("--diffnum", type=Path, default=Path("share/Scripts/DiffNum.pl"))
    args = parser.parse_args()
    if args.import_from_repository:
        import_reviewed(args.import_from_repository)
    elif args.run_dir:
        raise SystemExit(check(args.run_dir, args.diff_file, args.diffnum))
    else:
        parser.error("supply --run-dir or --import-from-repository")

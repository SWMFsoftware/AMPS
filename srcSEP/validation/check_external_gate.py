#!/usr/bin/env python3
"""Validate one external evidence class without merging validation claims.

Controlled analytical tests, native AMPS execution, real SWMF coupling, and
held-out spacecraft comparisons answer different scientific questions.  This
small command keeps the two external-data classes independently executable:
each supplied manifest and each referenced byte stream is authenticated by the
Step 15 validator, while observational coverage is checked across all supplied
events.  It intentionally does not synthesize missing evidence or convert an
INCOMPLETE class into PASS.
"""

import argparse
import pathlib
import sys
from typing import Optional, Sequence

import run_campaign


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--evidence-class",
        required=True,
        choices=("COUPLED_INTEGRATION", "OBSERVATIONAL_VALIDATION"),
        help="validate real SWMF coupling or held-out spacecraft evidence")
    parser.add_argument(
        "manifest", nargs="+",
        help="one or more srcsep-external-evidence-v1 manifest paths")
    return parser.parse_args(argv)


def validate_gate(args: argparse.Namespace) -> None:
    observed_families = set()
    for manifest_name in args.manifest:
        manifest, _ = run_campaign.validate_external_manifest(
            pathlib.Path(manifest_name).resolve(), args.evidence_class)
        run_campaign.require(
            manifest["status"] == "PASS",
            "external gate requires PASS manifest: {}".format(manifest_name))
        observed_families.update(
            metric["metric_family"]
            for metric in manifest["acceptance_metrics"])

    if args.evidence_class == "OBSERVATIONAL_VALIDATION":
        # Multiple held-out events may jointly close the required metric
        # families.  Validate the union here while retaining every event's own
        # independently checked status, provenance, and input checksums.
        missing = run_campaign.OBSERVATIONAL_METRIC_FAMILIES - observed_families
        run_campaign.require(
            not missing,
            "observational gate is missing metric families: {}".format(
                ", ".join(sorted(missing))))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_arguments(argv)
    try:
        validate_gate(args)
    except run_campaign.EvidenceError as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        return 2
    label = "SWMF" if args.evidence_class == "COUPLED_INTEGRATION" \
        else "OBSERVATIONAL"
    print("PASS {}-GATE: {} checksum-verified manifest(s)".format(
        label, len(args.manifest)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

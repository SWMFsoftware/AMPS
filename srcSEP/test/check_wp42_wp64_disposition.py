#!/usr/bin/env python3
"""Fail closed when WP42--WP64 status drifts from code, build, or runner.

The advanced overlay is deliberately split between dependency-light
experimental APIs and external qualification gates.  This checker prevents a
third, ambiguous state: a document cannot call a package production-active
while its implementation object is absent from the production archive or its
test is unavailable from the unified runner.
"""

from __future__ import print_function

import argparse
import json
import pathlib
import re
from typing import List, Sequence


ALLOWED_STATUS = frozenset(("integrated", "experimental", "superseded",
                            "remove", "external-gate"))
FORBIDDEN_DOCUMENT_CLAIMS = (
    "Production seam active",
    "All returned exit status zero",
    "one PASS for every work package",
)


def _documentation_claim_errors(text: str, name: str) -> List[str]:
    return ["{} retains unqualified historical claim: {}".format(name, claim)
            for claim in FORBIDDEN_DOCUMENT_CLAIMS if claim in text]


def audit(root: pathlib.Path) -> List[str]:
    root = root.resolve()
    errors: List[str] = []
    disposition_path = root / "WP42_WP64_DISPOSITION.json"
    try:
        disposition = json.loads(disposition_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        return ["cannot read disposition manifest: {}".format(error)]
    if disposition.get("schema") != "srcsep-wp42-wp64-disposition-v1":
        errors.append("disposition manifest schema is missing or unsupported")

    packages = disposition.get("work_packages")
    if not isinstance(packages, list):
        return errors + ["disposition manifest has no work_packages array"]
    expected = ["WP{:02d}".format(number) for number in range(42, 65)]
    observed = [str(item.get("id", "")) for item in packages
                if isinstance(item, dict)]
    if observed != expected:
        errors.append("work-package IDs are not the exact ordered WP42--WP64 set")
    for item in packages:
        if not isinstance(item, dict):
            errors.append("work-package disposition is not an object")
            continue
        identifier = str(item.get("id", "unknown"))
        if item.get("status") not in ALLOWED_STATUS:
            errors.append(identifier + " has an unsupported status")
        for required in ("symbols", "reason", "promotion_gate"):
            if not item.get(required):
                errors.append(identifier + " omits " + required)

    makefile = (root / "makefile").read_text(encoding="utf-8")
    experimental = disposition.get("experimental_sources", [])
    if not isinstance(experimental, list) or not experimental:
        errors.append("experimental_sources is empty or invalid")
    else:
        for relative in experimental:
            source = root / str(relative)
            if not source.is_file():
                errors.append("experimental source is absent: " + str(relative))
                continue
            object_name = str(pathlib.PurePosixPath(str(relative)).with_suffix(".o"))
            if object_name in makefile:
                errors.append(
                    "experimental source entered the production object list: " +
                    object_name)

    common_header = root.parent / "src" / "models" / "sep_common" / \
        "sep_transport_common.h"
    try:
        common_text = common_header.read_text(encoding="utf-8")
    except OSError as error:
        errors.append("cannot read canonical transport header: {}".format(error))
        common_text = ""
    if "enum class ParkerMeasure" not in common_text:
        errors.append("accepted ParkerMeasure supporting contract is absent")

    wave_header = (root / "util" / "sep_focused_transport_core.h").read_text(
        encoding="utf-8")
    mover_source = (root / "focused_transport_mfp.cpp").read_text(
        encoding="utf-8")
    core_source = (root / "util" / "sep_focused_transport_mfp_core.cpp").read_text(
        encoding="utf-8")
    for token in ("resonantBranch", "preWaveMomentumKgMPerS",
                  "postWaveMomentumKgMPerS"):
        if token not in wave_header or ("emitted." + token) not in mover_source:
            errors.append("accepted event metadata is not wired: " + token)
    if core_source.find("waveAccumulator->Deposit(contribution)") < \
            core_source.find("contribution.postWaveMomentumKgMPerS"):
        errors.append("event contribution is deposited before event metadata commits")

    runner = (root / "test" / "run_tests.py").read_text(encoding="utf-8")
    for token in ("wp42-wp64-experimental", "wp59-wp64-native"):
        if token not in runner:
            errors.append("unified runner does not expose suite " + token)
    for token in ("test-wp42-wp64-experimental:", "test-wp59-wp64-native:"):
        if token not in makefile:
            errors.append("makefile does not expose target " + token[:-1])

    test_source = (root / "test" / "test_wp42_wp64.cpp").read_text(
        encoding="utf-8")
    for identifier in expected:
        if not re.search(r"\b{}\b".format(identifier), test_source):
            errors.append(identifier + " has no dependency-light contract test")

    for name in ("WP42_WP64_IMPLEMENTATION.md",
                 "WP42_WP64_VALIDATION_REPORT.md"):
        path = root / name
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as error:
            errors.append("cannot read {}: {}".format(name, error))
            continue
        errors.extend(_documentation_claim_errors(text, name))
        if "experimental" not in text.lower() or "external" not in text.lower():
            errors.append(name + " does not state experimental/external boundaries")

    return errors


def main(argv: Sequence[str] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path,
                        default=pathlib.Path(__file__).resolve().parents[1])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    errors = audit(args.root)
    if errors:
        for error in errors:
            print("FAIL B05-DISPOSITION: " + error)
        return 1
    print("PASS B05-DISPOSITION: WP42--WP64 have explicit build/test/evidence status")
    if args.self_test:
        negative = _documentation_claim_errors(
            "Production seam active", "negative-control.md")
        if not negative:
            print("FAIL B05-DISPOSITION: false production claim was accepted")
            return 1
        print("PASS B05-DISPOSITION-NEGATIVE: false production claim rejected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

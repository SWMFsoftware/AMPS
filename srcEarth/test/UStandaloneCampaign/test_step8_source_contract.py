#!/usr/bin/env python3
"""Audit Step-8 production wiring without requiring a linked AMPS binary.

This is a source-contract gate, not the numerical reference suite. It protects
orchestration features that can otherwise disappear while the small Python
fixtures continue to pass: immutable inputs, exact extraction schemas,
observation adapters, restart evidence, release aggregation, documentation, and
registration in the main Earth test list.
"""

import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    runner = (ROOT / "standalone_campaign" / "run_campaign.py").read_text(encoding="utf-8")
    campaign = (ROOT / "standalone_campaign" / "campaign.py").read_text(encoding="utf-8")
    adapters = (ROOT / "standalone_campaign" / "adapters.py").read_text(encoding="utf-8")
    release = (ROOT / "standalone_campaign" / "release_gate.py").read_text(encoding="utf-8")
    tests = (ROOT / "test" / "list").read_text(encoding="utf-8")
    root_readme = (ROOT / "README.md").read_text(encoding="utf-8")
    module_readme = (ROOT / "standalone_campaign" / "README.md").read_text(encoding="utf-8")
    test_readme = (HERE / "README.md").read_text(encoding="utf-8")
    example = json.loads(
        (ROOT / "examples" / "standalone_step8_campaign" / "campaign.json")
        .read_text(encoding="utf-8")
    )

    require("network/mutable resource is forbidden" in campaign,
            "mutable-network input rejection is missing")
    require("provenance statement" in campaign and "sha256" in campaign,
            "resource provenance/hash validation is missing")
    require("unresolved_support" in campaign and "angular_convergence" in campaign and
            "energy_convergence" in campaign,
            "frozen numerical-gate contract is missing")
    require("verified_restart" in runner and "SKIPPED_VERIFIED_PASS" in runner,
            "hash-verified restart wiring is missing")
    require("no model/observation comparison rows were produced" in runner and
            "comparison_fraction" in runner,
            "fail-closed comparison coverage gate is missing")
    require("expected_variables" in runner and "expected_row_count" in runner and
            '"units": extractor["units"]' in runner,
            "exact Tecplot extraction schema/unit gate is missing")
    require("GOES_EPEAD_DIRECTIONAL" in adapters and "PAMELA_CUTOFF" in adapters and
            "POES_METOP_MEPED_CUTOFF" in adapters and "REPT_PROTON_SPECTRUM" in adapters,
            "one or more required observation adapters are missing")
    require("cadence_average" in adapters and "model series does not cover" in adapters,
            "exposure averaging or its no-extrapolation gate is missing")
    require("validate_instrument_response" in adapters and "PRIMARY_ACCEPTANCE" in adapters,
            "response validation or the released C10 acceptance mapping is missing")
    require("O1" in release and "O2" in release and "O4" in release and
            "all_required_numerical_gates_pass" in release,
            "O1/O2/O4 release evidence aggregation is missing")
    require("P srcEarth/test/UStandaloneCampaign/run_test.sh" in tests,
            "Step-8 suite is not registered independently")
    require("Roadmap Step 8" in root_readme and
            "Existing gates remain unchanged" in module_readme and
            "Reference solutions" in test_readme,
            "Step-8 implementation/test documentation is incomplete")

    require(example.get("event_id") == "STEP8_SMOKE",
            "example manifest has no explicit event identity")
    require(all(item.get("provenance") for item in example["resources"]),
            "example manifest resource provenance is incomplete")
    extractor = example["execution"]["prediction_extractors"][0]
    require(extractor.get("units") == "1" and
            extractor.get("expected_row_count") == 1 and
            extractor.get("expected_variables"),
            "example extractor does not freeze units, variables, and row count")

    print("PASS: Step-8 source, evidence, adapters, release gate, and registration")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as exc:
        print("FAIL: %s" % exc, file=sys.stderr)
        raise SystemExit(1)

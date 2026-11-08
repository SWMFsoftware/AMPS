#!/usr/bin/env python3
"""Register the Moon I01--I10 M0 campaign without overstating execution.

The local M0 orchestrator performs several mutually incompatible native AMPS
rebuilds (SPICE neutral, no-SPICE neutral, ion, and canonical regression).
Its historical rebuild helper has not yet been converted to the mandatory
root/process preflight and fixed ``MAKEFLAGS=-j16`` contract.  Running it from
the global routine would therefore violate the repository build policy.

This adapter still audits the complete frozen I01--I10 runner/README inventory
and publishes one typed SKIP per scientific case.  A future change may replace
those SKIPs with imported outcomes only after the local multi-build runner is
made preflight-compliant and independently exercised.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, load_context, main_guard, write_result)


LOCAL_IDS = ["I{:02d}".format(index) for index in range(1, 11)]
EXPECTED_CASES = ["MOON-{}".format(identifier) for identifier in LOCAL_IDS]


def main() -> int:
    context = load_context("Audit the Moon I01-I10 M0 integration boundary")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    if (not requested or len(set(requested)) != len(requested) or
            any(identifier not in EXPECTED_CASES for identifier in requested)):
        raise RunnerError("invalid Moon M0 case selection")

    repo = Path(str(context["repo_root"]))
    local_root = repo / "srcMoon" / "test"
    missing = []
    for local_id in LOCAL_IDS:
        directory = local_root / "integration-tests" / local_id
        for name in ("test.py", "README.md"):
            path = directory / name
            if not path.is_file():
                missing.append(str(path))
    orchestrator = local_root / "run_integration_tests.py"
    if not orchestrator.is_file():
        missing.append(str(orchestrator))
    if missing:
        raise RunnerError("Moon M0 registry files are missing: {}".format(
            ", ".join(missing)))

    reason = (
        "not executed: local multi-build M0 orchestrator has not yet been "
        "converted to the mandatory native-build preflight")
    rows = [case(identifier, "SKIP", reason) for identifier in requested]
    for row in rows:
        print("[CASE] RUN  {} {}".format(row["id"], reason), flush=True)

    artifacts = Path(str(context["artifact_dir"]))
    artifacts.mkdir(parents=True, exist_ok=True)
    boundary = artifacts / "moon-m0-boundary.json"
    boundary.write_text(json.dumps({
        "schema_version": 1,
        "local_runner": str(orchestrator),
        "registered_local_ids": LOCAL_IDS,
        "execution": "not-run",
        "reason": reason,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return write_result(
        context, "SKIP", reason, cases=rows, evidence=[str(boundary)],
        metrics={"registered_case_count": len(LOCAL_IDS),
                 "selected_case_count": len(rows)})


if __name__ == "__main__":
    main_guard(main)

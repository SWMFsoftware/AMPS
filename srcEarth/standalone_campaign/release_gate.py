#!/usr/bin/env python3
"""Evaluate the Roadmap Step-8 standalone release gate.

Step 8 is complete only when the O1 and O2 campaign summaries are real PASS
records, every required numerical/observation gate passed without a
platform-specific scale, and the frozen O4 holdout manifest can be preflighted
without accessing the network.  This command deliberately does not rerun or
reinterpret either campaign; it consumes their machine-readable evidence and
fails closed when any required field is absent.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from standalone_campaign.campaign import (  # type: ignore
        CampaignError,
        atomic_write_json,
        load_and_validate_manifest,
    )
else:
    from .campaign import CampaignError, atomic_write_json, load_and_validate_manifest


NORMALIZATION_POLICY = "SHARED_EVENT_BOUNDARY_NO_PLATFORM_SCALE"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _read_json(path: Path, label: str) -> Mapping[str, object]:
    try:
        with path.expanduser().resolve().open(encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise CampaignError("cannot read %s %s: %s" % (label, path, exc)) from exc
    if not isinstance(value, dict):
        raise CampaignError("%s must be a JSON object: %s" % (label, path))
    return value


def _campaign_decision(
    path: Path, expected_event_id: str
) -> Dict[str, object]:
    """Validate one completed campaign summary without guessing defaults."""

    summary = _read_json(path, expected_event_id + " summary")
    checks = {
        "event_id": summary.get("event_id") == expected_event_id,
        "status": summary.get("status") == "PASS",
        "role": summary.get("role") == "VALIDATION",
        "frozen": summary.get("frozen") is True,
        "normalization_policy": (
            summary.get("normalization_policy") == NORMALIZATION_POLICY
        ),
        "numerical_gates": summary.get("all_required_numerical_gates_pass") is True,
        "immutable_identity": (
            bool(SHA256_RE.fullmatch(str(summary.get("manifest_sha256", ""))))
            and bool(SHA256_RE.fullmatch(str(summary.get("executable_sha256", ""))))
            and isinstance(summary.get("resource_count"), int)
            and int(summary["resource_count"]) >= 7
        ),
        "observation_comparisons": (
            isinstance(summary.get("observation_comparison_count"), int)
            and int(summary["observation_comparison_count"]) > 0
            and isinstance(summary.get("observation_comparison_row_count"), int)
            and int(summary["observation_comparison_row_count"]) > 0
            and isinstance(summary.get("required_observation_comparison_count"), int)
            and int(summary["required_observation_comparison_count"]) > 0
            and summary.get("observation_comparison_fail_count") == 0
            and summary.get("required_observation_comparison_fail_count") == 0
        ),
        "executed_runs": (
            isinstance(summary.get("run_count"), int)
            and int(summary["run_count"]) > 0
            and summary.get("run_pass_count") == summary.get("run_count")
            and summary.get("run_fail_count") == 0
            and summary.get("run_not_executed_count") == 0
        ),
    }
    return {
        "event_id": expected_event_id,
        "summary_path": str(path.expanduser().resolve()),
        "campaign_id": summary.get("campaign_id"),
        "manifest_sha256": summary.get("manifest_sha256"),
        "executable_sha256": summary.get("executable_sha256"),
        "checks": checks,
        "status": "PASS" if all(checks.values()) else "FAIL",
    }


def _holdout_decision(path: Path) -> Dict[str, object]:
    """Preflight the untouched O4 manifest and all of its hashed resources."""

    manifest = load_and_validate_manifest(path)
    checks = {
        "event_id": manifest.get("event_id") == "O4",
        "role": manifest.get("role") == "HOLDOUT",
        "frozen": manifest.get("frozen") is True,
        "normalization_policy": (
            manifest.get("normalization_policy") == NORMALIZATION_POLICY
        ),
        "resources_verified": bool(manifest.get("_verified_resources")),
    }
    return {
        "event_id": "O4",
        "manifest_path": str(path.expanduser().resolve()),
        "campaign_id": manifest.get("campaign_id"),
        "manifest_sha256": manifest.get("_manifest_sha256"),
        "checks": checks,
        "status": "READY" if all(checks.values()) else "FAIL",
    }


def evaluate_release_gate(
    o1_summary: Path, o2_summary: Path, holdout_manifest: Path
) -> Dict[str, object]:
    """Return complete evidence even when one input is malformed or missing."""

    decisions = []
    for label, callback in (
        ("O1", lambda: _campaign_decision(o1_summary, "O1")),
        ("O2", lambda: _campaign_decision(o2_summary, "O2")),
        ("O4", lambda: _holdout_decision(holdout_manifest)),
    ):
        try:
            decisions.append(callback())
        except CampaignError as exc:
            decisions.append({"event_id": label, "status": "FAIL", "error": str(exc)})
    passed = (
        decisions[0].get("status") == "PASS"
        and decisions[1].get("status") == "PASS"
        and decisions[2].get("status") == "READY"
    )
    return {
        "schema_version": "earth-standalone-release/v1",
        "status": "PASS" if passed else "FAIL",
        "normalization_policy": NORMALIZATION_POLICY,
        "decisions": decisions,
        "non_relaxation_rule": (
            "Failed numerical or observational gates remain failures; missing "
            "comparisons and unfrozen holdout evidence cannot be accepted."
        ),
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--o1-summary", required=True, type=Path)
    parser.add_argument("--o2-summary", required=True, type=Path)
    parser.add_argument("--holdout-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    result = evaluate_release_gate(
        args.o1_summary, args.o2_summary, args.holdout_manifest
    )
    atomic_write_json(args.output.expanduser().resolve(), result)
    print("STEP-8 STANDALONE RELEASE %s" % result["status"], flush=True)
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())

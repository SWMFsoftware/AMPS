#!/usr/bin/env python3
"""Create a complete, deliberately non-releasable Step-12 manifest skeleton.

The generator includes every fixed gate and both seven-product parity campaigns so a
site cannot accidentally start from an incomplete hand-written list.  Placeholder
paths and all-zero digests are intentional: ``run_release.py`` will fail until they
are replaced with real, hash-pinned artifacts and measured commands/resources.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from release_validation.contract import (  # type: ignore
        COMMON_PHYSICS_TAG,
        NORMALIZATION_POLICY,
        PHASE_1_SCOPE,
        REQUIRED_GATE_IDS,
        REQUIRED_PARITY_ROLES,
        SCHEMA,
        atomic_write_json,
    )
else:
    from .contract import (
        COMMON_PHYSICS_TAG,
        NORMALIZATION_POLICY,
        PHASE_1_SCOPE,
        REQUIRED_GATE_IDS,
        REQUIRED_PARITY_ROLES,
        SCHEMA,
        atomic_write_json,
    )


ZERO_SHA256 = "0" * 64


def _reference(path: str) -> Dict[str, str]:
    return {"path": path, "sha256": ZERO_SHA256}


def _campaign(kind: str) -> Dict[str, Any]:
    sampled = kind == "SAMPLED_ANALYTIC_MESH"
    pairs = []
    for role in REQUIRED_PARITY_ROLES:
        stem = role.lower()
        pairs.append({
            "role": role,
            "reference": _reference("evidence/%s/%s_reference.csv" % (kind.lower(), stem)),
            "candidate": _reference("evidence/%s/%s_candidate.csv" % (kind.lower(), stem)),
            "key_columns": ["REPLACE_KEY_COLUMN"],
            "exact_columns": ["REPLACE_OUTCOME_OR_IDENTITY_COLUMN"],
            "numeric_columns": ["REPLACE_NUMERIC_COLUMN"],
            "tolerance_class": "MESH" if sampled else "KERNEL",
            "rtol": 0.05 if sampled else 1.0e-10,
            "atol": 1.0e-12,
        })
    return {
        "campaign_id": "REPLACE_CAMPAIGN_ID",
        "kind": kind,
        "common_physics_tag": COMMON_PHYSICS_TAG,
        "phase_1_scope": PHASE_1_SCOPE,
        "reference_field_source": "ANALYTIC_OR_PHENOMENOLOGICAL" if sampled else "LIVE_SWMF",
        "candidate_field_source": "SAMPLED_MESH" if sampled else "OFFLINE_SNAPSHOT_REPLAY",
        "invariants": {
            "field_state_id": "REPLACE_FIELD_STATE_ID",
            "species": "PROTON",
            "epoch_utc": "REPLACE_ABSOLUTE_UTC",
            "frame": "GSM",
            "boundary_sha256": ZERO_SHA256,
            "integration_sha256": ZERO_SHA256,
            "output_schema_sha256": ZERO_SHA256,
        },
        "pairs": pairs,
    }


def build_skeleton() -> Dict[str, Any]:
    return {
        "schema": SCHEMA,
        "release_id": "REPLACE_RELEASE_ID",
        "phase_1_scope": PHASE_1_SCOPE,
        "common_physics": {
            "tag": COMMON_PHYSICS_TAG,
            "source_sha256": ZERO_SHA256,
            "standalone_build": _reference("evidence/standalone_build.json"),
            "coupled_build": _reference("evidence/swmf_coupled_build.json"),
        },
        "gate_evidence": [
            {"gate_id": gate_id,
             "evidence": _reference("evidence/gates/%s.json" % gate_id.replace("-", "_"))}
            for gate_id in REQUIRED_GATE_IDS
        ],
        "cross_path_campaigns": [
            _campaign("SAMPLED_ANALYTIC_MESH"),
            _campaign("SWMF_OFFLINE_REPLAY"),
        ],
        "holdout": {
            "gate_id": "O4", "frozen": True, "retuned": False,
            "selected_before_release": True, "outcome_reported": True,
            "normalization_policy": NORMALIZATION_POLICY,
            "preregistration_sha256": ZERO_SHA256,
        },
        "capabilities": _reference("capabilities.json"),
        "resource_estimates": _reference("resource_estimates.json"),
        "hooks": {
            "ccmc_archive": "REPLACE_CCMC_ARCHIVE_COMMAND",
            "swmf_snapshot_export": "REPLACE_SWMF_EXPORT_COMMAND",
            "swmf_offline_replay": "REPLACE_OFFLINE_REPLAY_COMMAND",
        },
        "commands": {
            "standalone_cutoff_only": "REPLACE_STANDALONE_CUTOFF_COMMAND",
            "standalone_combined": "REPLACE_STANDALONE_COMBINED_COMMAND",
            "swmf_cutoff_only": "REPLACE_SWMF_CUTOFF_COMMAND",
            "swmf_combined": "REPLACE_SWMF_COMBINED_COMMAND",
            "offline_replay": "REPLACE_OFFLINE_REPLAY_COMMAND",
        },
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path,
                        help="destination JSON skeleton; existing files are atomically replaced")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    atomic_write_json(args.output.expanduser().resolve(), build_skeleton())
    print("RESULT: SKELETON_WRITTEN", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

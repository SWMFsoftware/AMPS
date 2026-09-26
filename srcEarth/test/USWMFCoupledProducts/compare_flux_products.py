#!/usr/bin/env python3
"""Compare a complete live Step-11 product set with offline Mode3D replay.

The comparator consumes the live ``swmf_flux_spectrum_manifest*.json`` and derives
the corresponding replay names by replacing the live content/time suffix with an
optional standalone suffix.  Every Tecplot artifact is parsed with the strict Step-10
reader.  Snapshot identity, boundary-spectrum time, product controls, observation
state, units, response/channel schemas, VARIABLES, ZONE topology, and numeric rows
must agree before the product set passes.

Exact numeric comparison is the default because live and replay invoke the same
backward-characteristic solver on the same exported field.  Non-zero tolerances are
accepted only when explicitly provided by a preregistered validation campaign; this
script never widens them automatically.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Sequence


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "USWMFCoupledAccess"))
from compare_cutoff_access import (  # noqa: E402
    ComparisonFailure,
    compare_products,
    parse_product,
)


SCHEMA = "sep-in-geospace/swmf-coupled-products/v1"
REQUIRED_MANIFEST_KEYS = (
    "snapshot_id",
    "content_fingerprint",
    "mesh_revision",
    "epoch_utc",
    "boundary_spectrum_evaluation_epoch_utc",
    "output_suffix",
    "product_control_fingerprint",
    "boundary_spectrum_fingerprint",
    "channel_schema_fingerprint",
    "detector_response_fingerprint",
    "observation_state_fingerprint",
)
REQUIRED_PRODUCT_AUX = (
    "PHASE_1_INTERPRETATION",
    "STEP6_CHARACTERISTIC_MAPPING",
    "SNAPSHOT_ID",
    "SNAPSHOT_EPOCH_UTC",
    "SNAPSHOT_MESH_REVISION",
    "SNAPSHOT_CONTENT_FINGERPRINT",
    "OUTER_BOUNDARY_POLICY",
    "BOUNDARY_SPECTRUM_EVALUATION_EPOCH_UTC",
    "ACTIVE_SPECTRUM_TABLE_EPOCH_UTC",
    "PRODUCT_CONTROL_FINGERPRINT",
    "BOUNDARY_SPECTRUM_FINGERPRINT",
    "CHANNEL_SCHEMA_FINGERPRINT",
    "DETECTOR_RESPONSE_FINGERPRINT",
    "OBSERVATION_STATE_FINGERPRINT",
    "SPECTRUM_ENERGY_BASIS",
    "SPECTRUM_MASS_NUMBER",
    "SPECTRUM_INTENSITY_UNIT",
    "SPECTRUM_RELATIVE_UNCERTAINTY",
    "SPECTRUM_TEMPORAL_STATUS",
    "SPECTRUM_TEMPORAL_GAP",
    "SPECTRUM_TEMPORAL_FRACTION",
    "PLANAR_FLUX_CONVENTION",
    "MAXIMUM_UNRESOLVED_FRACTION",
    "RESPONSE_WEIGHTED_UNRESOLVED_UPPER_BOUND",
    "UNRESOLVED_TOLERANCE",
)
MANIFEST_TO_AUX = {
    "snapshot_id": "SNAPSHOT_ID",
    "epoch_utc": "SNAPSHOT_EPOCH_UTC",
    "mesh_revision": "SNAPSHOT_MESH_REVISION",
    "content_fingerprint": "SNAPSHOT_CONTENT_FINGERPRINT",
    "boundary_spectrum_evaluation_epoch_utc":
        "BOUNDARY_SPECTRUM_EVALUATION_EPOCH_UTC",
    "active_spectrum_table_epoch_utc": "ACTIVE_SPECTRUM_TABLE_EPOCH_UTC",
    "product_control_fingerprint": "PRODUCT_CONTROL_FINGERPRINT",
    "boundary_spectrum_fingerprint": "BOUNDARY_SPECTRUM_FINGERPRINT",
    "channel_schema_fingerprint": "CHANNEL_SCHEMA_FINGERPRINT",
    "detector_response_fingerprint": "DETECTOR_RESPONSE_FINGERPRINT",
    "observation_state_fingerprint": "OBSERVATION_STATE_FINGERPRINT",
    "outer_boundary_policy": "OUTER_BOUNDARY_POLICY",
    "spectrum_energy_basis": "SPECTRUM_ENERGY_BASIS",
    "spectrum_intensity_unit": "SPECTRUM_INTENSITY_UNIT",
    "spectrum_temporal_status": "SPECTRUM_TEMPORAL_STATUS",
}
ROLE_VARIABLES = {
    "DENSITY": {"N_m^-3", "N_lower_m^-3", "N_upper_m^-3"},
    "SPECTRUM": {
        "E_MeV",
        "T",
        "T_lower",
        "T_upper",
        "unresolved_fraction",
        "J_boundary_perMeV",
        "J_local_perMeV",
        "J_local_lower_perMeV",
        "J_local_upper_perMeV",
        "J_omni_perMeV",
        "J_omni_lower_perMeV",
        "J_omni_upper_perMeV",
        "J_planar_perMeV",
        "J_planar_lower_perMeV",
        "J_planar_upper_perMeV",
    },
    "FLUX": {
        "F_tot_m2s1",
        "F_tot_lower_m2s1",
        "F_tot_upper_m2s1",
        "F_planar_m2s1",
        "F_planar_lower_m2s1",
        "F_planar_upper_m2s1",
    },
    "DENSITY_FLUX": {
        "N_m^-3",
        "N_lower_m^-3",
        "N_upper_m^-3",
        "F_tot_m2s1",
        "F_tot_lower_m2s1",
        "F_tot_upper_m2s1",
        "F_planar_m2s1",
        "F_planar_lower_m2s1",
        "F_planar_upper_m2s1",
    },
    "TERMINATION": {
        "location_index",
        "E_MeV",
        "N_sampled",
        "N_retried",
        "N_resolved",
        "N_allowed",
        "T",
        "T_lower",
        "T_upper",
        "unresolved_fraction",
    },
}


def _require_manifest(manifest_path: Path) -> dict:
    if not manifest_path.is_file():
        raise ComparisonFailure(f"manifest does not exist: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ComparisonFailure(f"cannot parse {manifest_path}: {error}") from error
    if manifest.get("schema") != SCHEMA:
        raise ComparisonFailure("unsupported or missing Step-11 manifest schema")
    if manifest.get("RESULT") != "PASS":
        raise ComparisonFailure("live Step-11 manifest is not PASS")
    if manifest.get("phase_1_interpretation") != "INSTANTANEOUS_QUASI_STATIC":
        raise ComparisonFailure("manifest lacks the Phase-1 quasi-static limitation")
    if manifest.get("characteristic_mapping") != "STATIC_MAGNETIC":
        raise ComparisonFailure("manifest does not use the released static mapping")
    for key in REQUIRED_MANIFEST_KEYS:
        if not isinstance(manifest.get(key), str) or not manifest[key]:
            raise ComparisonFailure(f"manifest is missing {key}")
    if manifest["epoch_utc"] != manifest["boundary_spectrum_evaluation_epoch_utc"]:
        raise ComparisonFailure("field and boundary-spectrum epochs differ in manifest")

    termination = manifest.get("termination")
    if not isinstance(termination, dict):
        raise ComparisonFailure("manifest is missing termination accounting")
    required_counts = ("sampled", "retried", "resolved", "allowed")
    if any(not isinstance(termination.get(key), int) for key in required_counts):
        raise ComparisonFailure("manifest termination counters are not integers")
    sampled = termination["sampled"]
    resolved = termination["resolved"]
    allowed = termination["allowed"]
    counts = termination.get("counts")
    if (
        sampled <= 0
        or not 0 <= allowed <= resolved <= sampled
        or not isinstance(counts, list)
        or not counts
        or any(not isinstance(value, int) or value < 0 for value in counts)
        or sum(counts) != sampled
    ):
        raise ComparisonFailure("manifest termination counters do not close")
    try:
        maximum = float(termination["maximum_unresolved_fraction"])
        response_upper = float(termination["response_weighted_unresolved_upper_bound"])
        tolerance = float(termination["unresolved_tolerance"])
    except (KeyError, TypeError, ValueError) as error:
        raise ComparisonFailure("manifest unresolved gate is malformed") from error
    if not (
        math.isfinite(maximum)
        and math.isfinite(response_upper)
        and math.isfinite(tolerance)
        and 0.0 <= maximum <= tolerance <= 1.0
        and response_upper == maximum
    ):
        raise ComparisonFailure("manifest unresolved fraction exceeds its retained gate")

    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ComparisonFailure("manifest has no product artifacts")
    roles: List[str] = []
    paths: List[str] = []
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ComparisonFailure("manifest artifact is not an object")
        role, path = artifact.get("role"), artifact.get("path")
        if role not in {"DENSITY", "FLUX", "DENSITY_FLUX", "SPECTRUM", "TERMINATION"}:
            raise ComparisonFailure(f"unrecognized artifact role: {role!r}")
        if not isinstance(path, str) or not path or path in paths:
            raise ComparisonFailure("artifact path is empty or duplicated")
        roles.append(role)
        paths.append(path)
    if "SPECTRUM" not in roles or "TERMINATION" not in roles:
        raise ComparisonFailure("product set lacks spectrum or termination evidence")
    if "DENSITY_FLUX" not in roles and not ({"DENSITY", "FLUX"} <= set(roles)):
        raise ComparisonFailure("product set lacks density and integral flux")
    return manifest


def _replay_name(live_name: str, live_suffix: str, replay_suffix: str) -> str:
    ending = f"{live_suffix}.dat"
    if not live_name.endswith(ending):
        raise ComparisonFailure(
            f"live artifact {live_name!r} is not bound to manifest output_suffix"
        )
    return live_name[: -len(ending)] + replay_suffix + ".dat"


def compare_product_set(
    manifest_path: Path,
    replay_dir: Path,
    *,
    replay_suffix: str = "",
    rtol: float = 0.0,
    atol: float = 0.0,
) -> dict:
    """Validate provenance and compare every manifest-enrolled numeric artifact."""
    manifest = _require_manifest(manifest_path)
    if not replay_dir.is_dir():
        raise ComparisonFailure(f"replay directory does not exist: {replay_dir}")

    product_summaries: List[dict] = []
    manifest_dir = manifest_path.resolve().parent
    live_suffix = manifest["output_suffix"]
    channels = manifest.get("channels", [])
    responses = manifest.get("detector_responses", [])
    if not isinstance(channels, list) or not all(isinstance(x, str) for x in channels):
        raise ComparisonFailure("manifest channel schema is malformed")
    if not isinstance(responses, list) or not all(isinstance(x, str) for x in responses):
        raise ComparisonFailure("manifest detector schema is malformed")

    for artifact in manifest["artifacts"]:
        live_path = Path(artifact["path"])
        if not live_path.is_absolute():
            live_path = manifest_dir / live_path
        replay_path = replay_dir / _replay_name(
            Path(artifact["path"]).name, live_suffix, replay_suffix
        )
        live = parse_product(live_path)
        replay = parse_product(replay_path)

        for key in REQUIRED_PRODUCT_AUX:
            if not live.aux.get(key) or not replay.aux.get(key):
                raise ComparisonFailure(
                    f"{artifact['role']} product is missing required AUXDATA {key}"
                )
            if live.aux[key] != replay.aux[key]:
                raise ComparisonFailure(
                    f"{artifact['role']} live/replay metadata differs for {key}"
                )
        for manifest_key, aux_key in MANIFEST_TO_AUX.items():
            value = manifest.get(manifest_key)
            if value is not None and str(value) != live.aux.get(aux_key):
                raise ComparisonFailure(
                    f"manifest/artifact mismatch for {manifest_key} ({aux_key})"
                )
        if live.aux["PHASE_1_INTERPRETATION"] != "INSTANTANEOUS_QUASI_STATIC":
            raise ComparisonFailure("artifact lacks the Phase-1 interpretation")
        if live.aux["STEP6_CHARACTERISTIC_MAPPING"] != "STATIC_MAGNETIC":
            raise ComparisonFailure("artifact uses an unreleased characteristic mapping")
        missing_variables = ROLE_VARIABLES[artifact["role"]] - set(live.variables)
        if missing_variables:
            raise ComparisonFailure(
                f"{artifact['role']} artifact lacks required value/bound columns: "
                + ", ".join(sorted(missing_variables))
            )

        summary = compare_products(live, replay, rtol=rtol, atol=atol)
        summary["role"] = artifact["role"]
        product_summaries.append(summary)

        if artifact["role"] in {"FLUX", "DENSITY_FLUX"}:
            for name in channels:
                required = {
                    f"F_{name}_m2s1",
                    f"F_{name}_lower_m2s1",
                    f"F_{name}_upper_m2s1",
                }
                if not required <= set(live.variables):
                    raise ComparisonFailure(
                        f"flux artifact lacks nominal/lower/upper channel {name}"
                    )
            for name in responses:
                required = {
                    f"R_{name}_s1",
                    f"R_{name}_lower_s1",
                    f"R_{name}_upper_s1",
                }
                if not required <= set(live.variables):
                    raise ComparisonFailure(
                        f"flux artifact lacks nominal/lower/upper detector response {name}"
                    )

    return {
        "schema": "sep-in-geospace/swmf-coupled-products-comparison/v1",
        "RESULT": "PASS",
        "snapshot_id": manifest["snapshot_id"],
        "product_control_fingerprint": manifest["product_control_fingerprint"],
        "artifacts_compared": len(product_summaries),
        "rtol": rtol,
        "atol": atol,
        "maximum_absolute_difference": max(
            summary["max_abs_difference"] for summary in product_summaries
        ),
        "maximum_relative_difference": max(
            summary["max_relative_difference"] for summary in product_summaries
        ),
        "products": product_summaries,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live-manifest", required=True, type=Path)
    parser.add_argument("--replay-dir", required=True, type=Path)
    parser.add_argument("--replay-suffix", default="")
    parser.add_argument("--rtol", type=float, default=0.0)
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--json", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        summary = compare_product_set(
            args.live_manifest,
            args.replay_dir,
            replay_suffix=args.replay_suffix,
            rtol=args.rtol,
            atol=args.atol,
        )
        if args.json:
            args.json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(summary, sort_keys=True))
        print("RESULT: PASS")
        return 0
    except (ComparisonFailure, OSError) as error:
        failure: Dict[str, str] = {"RESULT": "FAIL", "message": str(error)}
        if args.json:
            args.json.write_text(json.dumps(failure, indent=2) + "\n", encoding="utf-8")
        print(f"RESULT: FAIL: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

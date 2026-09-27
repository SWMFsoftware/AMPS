#!/usr/bin/env python3
"""Evaluate the complete Roadmap Step-12 release contract.

The command consumes hash-pinned build, numerical, functional, observational, and
cross-path evidence.  It does not run missing tests, manufacture defaults, or turn a
dry run into a PASS.  A failure always produces a JSON summary, an explicit
``RESULT: FAIL`` line, and a non-zero exit status.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from release_validation.compare_cross_path import compare_campaigns  # type: ignore
    from release_validation.contract import (  # type: ignore
        COMMON_PHYSICS_TAG,
        NORMALIZATION_POLICY,
        REQUIRED_GATE_IDS,
        REQUIRED_HOOKS,
        REQUIRED_COMMANDS,
        SUMMARY_SCHEMA,
        ReleaseValidationError,
        ResourceResolver,
        atomic_write_json,
        manifest_sha256,
        read_json,
        require,
        require_mapping,
        require_nonempty_string,
        require_sequence,
        require_sha256,
        validate_build_manifest,
        validate_capabilities,
        validate_gate_evidence,
        validate_named_commands,
        validate_release_manifest_shape,
        validate_resource_estimates,
    )
else:
    from .compare_cross_path import compare_campaigns
    from .contract import (
        COMMON_PHYSICS_TAG,
        NORMALIZATION_POLICY,
        REQUIRED_GATE_IDS,
        REQUIRED_HOOKS,
        REQUIRED_COMMANDS,
        SUMMARY_SCHEMA,
        ReleaseValidationError,
        ResourceResolver,
        atomic_write_json,
        manifest_sha256,
        read_json,
        require,
        require_mapping,
        require_nonempty_string,
        require_sequence,
        require_sha256,
        validate_build_manifest,
        validate_capabilities,
        validate_gate_evidence,
        validate_named_commands,
        validate_release_manifest_shape,
        validate_resource_estimates,
    )


def _verify_pair_resources(manifest: Mapping[str, Any], resolver: ResourceResolver) -> None:
    """Hash parity inputs before accepting a restart shortcut.

    Restart may skip row-by-row comparisons only after every referenced byte stream is
    known to be unchanged.  Merely matching the manifest hash would miss an artifact
    modified in place after the original PASS.
    """

    for campaign_index, raw_campaign in enumerate(require_sequence(manifest.get("cross_path_campaigns"), "cross_path_campaigns")):
        campaign = require_mapping(raw_campaign, "cross_path_campaigns[%d]" % campaign_index)
        for pair_index, raw_pair in enumerate(require_sequence(campaign.get("pairs"), "campaign.pairs")):
            pair = require_mapping(raw_pair, "campaign pair %d" % pair_index)
            role = require_nonempty_string(pair.get("role"), "campaign pair role")
            resolver.resolve(pair.get("reference"), "%s reference" % role)
            resolver.resolve(pair.get("candidate"), "%s candidate" % role)


def _validate_holdout(value: Any) -> Dict[str, Any]:
    holdout = require_mapping(value, "holdout")
    require(holdout.get("gate_id") == "O4", "holdout.gate_id must be O4")
    require(holdout.get("frozen") is True, "O4 holdout must be frozen")
    require(holdout.get("retuned") is False, "O4 holdout must not be retuned")
    require(holdout.get("selected_before_release") is True,
            "O4 must have been selected before the release decision")
    require(holdout.get("outcome_reported") is True,
            "O4 outcome must be reported even when it prevents release")
    require(holdout.get("normalization_policy") == NORMALIZATION_POLICY,
            "O4 uses an unapproved normalization policy")
    freeze_digest = require_sha256(holdout.get("preregistration_sha256"),
                                   "holdout.preregistration_sha256")
    return {"gate_id": "O4", "frozen": True, "retuned": False,
            "preregistration_sha256": freeze_digest}


def _validate_gate_registry(
    manifest: Mapping[str, Any], resolver: ResourceResolver, revision: str
) -> List[Dict[str, Any]]:
    entries = require_sequence(manifest.get("gate_evidence"), "gate_evidence")
    by_gate: Dict[str, Mapping[str, Any]] = {}
    for index, raw in enumerate(entries):
        entry = require_mapping(raw, "gate_evidence[%d]" % index)
        gate_id = require_nonempty_string(entry.get("gate_id"), "gate_evidence[%d].gate_id" % index)
        require(gate_id not in by_gate, "duplicate gate evidence %s" % gate_id)
        by_gate[gate_id] = entry
    require(set(by_gate) == set(REQUIRED_GATE_IDS),
            "gate registry differs from the fixed Phase-1 matrix; missing=%s extra=%s" %
            (sorted(set(REQUIRED_GATE_IDS) - set(by_gate)),
             sorted(set(by_gate) - set(REQUIRED_GATE_IDS))))

    results: List[Dict[str, Any]] = []
    for gate_id in REQUIRED_GATE_IDS:
        path = resolver.resolve(by_gate[gate_id].get("evidence"), "%s evidence file" % gate_id)
        evidence = read_json(path, "%s evidence" % gate_id)
        require(evidence.get("source_revision") == revision,
                "%s evidence was produced from a different source revision" % gate_id)
        results.append(validate_gate_evidence(evidence, gate_id, resolver, COMMON_PHYSICS_TAG))
    return results


def _prepare_inputs(
    manifest_path: Path, manifest: Mapping[str, Any], resolver: ResourceResolver
) -> Dict[str, Any]:
    """Validate and hash every release input except numeric parity rows."""

    validate_release_manifest_shape(manifest)
    common = require_mapping(manifest.get("common_physics"), "common_physics")
    tag = require_nonempty_string(common.get("tag"), "common_physics.tag")
    source_digest = require_sha256(common.get("source_sha256"), "common_physics.source_sha256")
    standalone_path = resolver.resolve(common.get("standalone_build"), "standalone build manifest")
    coupled_path = resolver.resolve(common.get("coupled_build"), "coupled build manifest")
    standalone = validate_build_manifest(read_json(standalone_path, "standalone build manifest"),
                                         "STANDALONE", tag, source_digest)
    coupled = validate_build_manifest(read_json(coupled_path, "coupled build manifest"),
                                      "SWMF_COUPLED", tag, source_digest)
    require(standalone["revision"] == coupled["revision"],
            "standalone and coupled builds use different source revisions")
    gates = _validate_gate_registry(manifest, resolver, str(standalone["revision"]))

    capability_path = resolver.resolve(manifest.get("capabilities"), "capability table")
    resource_path = resolver.resolve(manifest.get("resource_estimates"), "resource estimate table")
    capabilities = validate_capabilities(read_json(capability_path, "capability table"))
    resources = validate_resource_estimates(read_json(resource_path, "resource estimate table"))
    holdout = _validate_holdout(manifest.get("holdout"))
    commands = validate_named_commands(manifest.get("commands"), REQUIRED_COMMANDS, "commands")
    hooks = validate_named_commands(manifest.get("hooks"), REQUIRED_HOOKS, "hooks")
    _verify_pair_resources(manifest, resolver)
    return {"builds": [standalone, coupled], "gates": gates, "capabilities": capabilities,
            "resource_estimates": resources, "holdout": holdout,
            "commands": commands, "hooks": hooks}


def _restart_matches(
    output_path: Path, release_id: str, manifest_digest: str, input_fingerprint: str
) -> bool:
    if not output_path.is_file():
        return False
    try:
        previous = read_json(output_path, "previous Step-12 summary")
    except ReleaseValidationError:
        return False
    return (
        previous.get("schema") == SUMMARY_SCHEMA
        and previous.get("status") == "PASS"
        and previous.get("release_id") == release_id
        and previous.get("manifest_sha256") == manifest_digest
        and previous.get("input_fingerprint") == input_fingerprint
    )


def evaluate_release(
    manifest_path: Path, output_path: Path, restart: bool = False,
    validate_only: bool = False
) -> Tuple[Dict[str, Any], bool]:
    """Return the report and whether a verified previous PASS was reused."""

    manifest_path = manifest_path.expanduser().resolve()
    output_path = output_path.expanduser().resolve()
    manifest = read_json(manifest_path, "release manifest")
    resolver = ResourceResolver(manifest_path)
    prepared = _prepare_inputs(manifest_path, manifest, resolver)
    manifest_digest = manifest_sha256(manifest_path)
    input_fingerprint = resolver.fingerprint(manifest_digest)
    release_id = require_nonempty_string(manifest.get("release_id"), "release_id")

    if restart and not validate_only and _restart_matches(
        output_path, release_id, manifest_digest, input_fingerprint
    ):
        # Do not rewrite a verified report: preserving its bytes makes the restart
        # decision itself independently auditable.
        return read_json(output_path, "previous Step-12 summary"), True

    if validate_only:
        report = {
            "schema": SUMMARY_SCHEMA,
            "status": "VALID_NOT_RELEASE",
            "release_id": release_id,
            "manifest_sha256": manifest_digest,
            "input_fingerprint": input_fingerprint,
            "message": "All configuration and hashes are valid; numeric parity was not evaluated.",
            "prepared": prepared,
        }
        atomic_write_json(output_path, report)
        return report, False

    parity = compare_campaigns(manifest, resolver)
    require(parity.get("status") == "PASS", "one or more cross-path comparisons failed")
    # compare_campaigns visits the same hash-pinned resources; recompute after it in
    # case a future schema adds an input not known to the preflight collector.
    input_fingerprint = resolver.fingerprint(manifest_digest)
    report = {
        "schema": SUMMARY_SCHEMA,
        "status": "PASS",
        "release_id": release_id,
        "manifest_sha256": manifest_digest,
        "input_fingerprint": input_fingerprint,
        "common_physics_tag": COMMON_PHYSICS_TAG,
        "prepared": prepared,
        "cross_path": parity,
        "non_relaxation_rule": (
            "Missing, failed, dry-run, retuned, unpinned, or threshold-relaxed evidence "
            "is a release failure."
        ),
    }
    atomic_write_json(output_path, report)
    return report, False


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path,
                        help="hash-pinned Step-12 release manifest")
    parser.add_argument("--output", type=Path, default=Path("step12_release_summary.json"),
                        help="atomic JSON summary (default: step12_release_summary.json)")
    parser.add_argument("--restart", action="store_true",
                        help="reuse a prior PASS only when manifest and every resource hash still match")
    parser.add_argument("--validate-only", action="store_true",
                        help="validate configuration/hashes without claiming a release PASS")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        report, reused = evaluate_release(
            args.manifest, args.output, restart=args.restart,
            validate_only=args.validate_only
        )
    except ReleaseValidationError as exc:
        report = {"schema": SUMMARY_SCHEMA, "status": "FAIL", "error": str(exc)}
        atomic_write_json(args.output.expanduser().resolve(), report)
        print("RESULT: FAIL", flush=True)
        print("Step-12 release failure: %s" % exc, file=sys.stderr, flush=True)
        return 2
    if report.get("status") == "VALID_NOT_RELEASE":
        print("RESULT: VALID_NOT_RELEASE", flush=True)
        return 0
    print("RESULT: PASS", flush=True)
    if reused:
        print("RESTART: SKIPPED_VERIFIED_PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Check Step 15 report separation and observational-label safeguards."""

import hashlib
import json
import pathlib
import sys
import tempfile
from types import SimpleNamespace


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: test_campaign_governance.py REPORT SOURCE_ROOT")
    report_path = pathlib.Path(sys.argv[1]).resolve()
    source_root = pathlib.Path(sys.argv[2]).resolve()
    sys.path.insert(0, str(source_root / "validation"))
    import run_campaign  # pylint: disable=import-error,import-outside-toplevel
    import check_external_gate  # pylint: disable=import-error,import-outside-toplevel

    with report_path.open("r", encoding="utf-8") as stream:
        report = json.load(stream)
    require(report["schema"] == "srcsep-validation-campaign-v1",
            "campaign schema identity changed")
    require(report["release_status"] == "INCOMPLETE",
            "a source-only campaign must not claim release readiness")
    require(report["sections"]["numerical_verification"]["status"] == "PASS",
            "VAL01 did not remain separate and passing")
    require(report["sections"]["cross_mover_verification"]["status"] == "PASS",
            "VAL02 did not remain separate and passing")
    require(report["sections"]["cross_model_verification"]["status"] == "PASS",
            "VAL03 did not remain separate and passing")
    require(report["sections"]["coupled_integration"]["status"] == "INCOMPLETE",
            "SWCME success must not conceal missing real SWMF evidence")
    require(report["sections"]["observational_validation"]["status"] == "INCOMPLETE",
            "missing spacecraft evidence must remain explicit")

    # Exercise the semantic barrier with a structurally complete record.  The
    # input bytes and checksum are real, but this record deliberately carries a
    # prohibited label and must therefore never acquire observational status.
    with tempfile.TemporaryDirectory(prefix="srcsep-step15-governance.") as temp_name:
        temp_dir = pathlib.Path(temp_name)
        input_path = temp_dir / "unit-input.txt"
        input_path.write_text("governance boundary input\n", encoding="utf-8")
        digest = hashlib.sha256(input_path.read_bytes()).hexdigest()
        manifest = {
            "schema": "srcsep-external-evidence-v1",
            "campaign_id": "governance-unit",
            "evidence_class": "OBSERVATIONAL_VALIDATION",
            "data_class": "SPACECRAFT_OBSERVATION",
            "status": "PASS",
            "configuration": {"transport_mover": "parker"},
            "inputs": [{
                "path": input_path.name,
                "sha256": digest,
                "role": "governance-input",
                "source_url_or_pid": "local-unit-record",
                "access_utc": "2000-01-01T00:00:00Z",
                "license_or_acknowledgment": "unit-only",
            }],
            "random_seeds": [],
            "compiler": {"command": "none", "version": "unit", "flags": []},
            "output_schema": {
                "name": "unit", "version": "1",
                "variables": [{"name": "intensity", "units": "arbitrary"}],
            },
            "event": {
                "id": "unit", "start_utc": "2000-01-01T00:00:00Z",
                "end_utc": "2000-01-01T01:00:00Z",
            },
            "observers": [{
                "mission": "unit", "instrument": "unit", "product_level": "unit",
                "product_version": "1", "variables": ["intensity"],
                "cadence_s": 60, "coordinate_system": "unit",
                "quality_flags": "unit",
            }],
            "uncertainty": {"method": "unit", "confidence_level": 0.95},
            "held_out": True,
            "acceptance_metrics": [{
                "name": "unit onset", "metric_family": "onset", "value": 0,
                "tolerance": 1, "comparison": "<=", "units": "s",
            }],
            "notes": "synthetic record used only to test rejection",
        }
        manifest_path = temp_dir / "observation.json"
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        rejected = False
        try:
            run_campaign.validate_external_manifest(
                manifest_path, "OBSERVATIONAL_VALIDATION")
        except run_campaign.EvidenceError:
            rejected = True
        require(rejected, "prohibited observational label was accepted")

        # Exercise the independently invocable external gates with genuine
        # bytes and complete schemas.  The observation fixture covers every
        # required family across one held-out record; a second invocation then
        # removes five families and must fail even though that single manifest
        # is structurally valid.  This proves class completeness is evaluated
        # by the observational gate, not inferred from another test class.
        manifest.pop("notes")
        # WP39 upgrades a free-form label to a complete instrument-space
        # transformation contract.  The governance fixture uses dimensionless
        # unit responses; it validates schema separation only and is not
        # promoted to a scientific observation in the campaign report.
        manifest["configuration"]["forward_operator"] = {
            "version": "srcsep-observation-forward-v1",
            "energy_response": "unit-response",
            "angular_response": "unit-aperture",
            "cadence_s": 60,
            "species": "proton",
            "dead_time_s": 0,
            "saturation_policy": "reject-saturated",
            "background_subtraction": "unit-background",
            "uncertainty_propagation": "poisson-plus-background",
        }
        manifest["observers"].append(dict(manifest["observers"][0]))
        manifest["observers"][1]["mission"] = "unit-second-observer"
        manifest["acceptance_metrics"] = [
            {
                "name": "unit {}".format(family),
                "metric_family": family,
                "value": 0,
                "tolerance": 1,
                "comparison": "<=",
                "units": "dimensionless",
            }
            for family in sorted(run_campaign.OBSERVATIONAL_METRIC_FAMILIES)
        ]
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        check_external_gate.validate_gate(SimpleNamespace(
            evidence_class="OBSERVATIONAL_VALIDATION",
            manifest=[str(manifest_path)]))

        partial_manifest = dict(manifest)
        partial_manifest["acceptance_metrics"] = [
            manifest["acceptance_metrics"][0]
        ]
        partial_path = temp_dir / "partial-observation.json"
        partial_path.write_text(json.dumps(partial_manifest), encoding="utf-8")
        rejected = False
        try:
            check_external_gate.validate_gate(SimpleNamespace(
                evidence_class="OBSERVATIONAL_VALIDATION",
                manifest=[str(partial_path)]))
        except run_campaign.EvidenceError:
            rejected = True
        require(rejected,
                "observational gate accepted incomplete metric-family coverage")

        # A real SWMF gate has a different data class and requires both replay
        # fidelity and transport-response metrics.  Keeping this fixture
        # separate prevents a spacecraft record from satisfying SWMF coupling.
        swmf_manifest = {
            "schema": "srcsep-external-evidence-v1",
            "campaign_id": "unit-swmf-gate",
            "evidence_class": "COUPLED_INTEGRATION",
            "data_class": "SWMF_OUTPUT",
            "status": "PASS",
            "configuration": {
                "srcsep_parameter_file": "unit.input",
                "swmf_run_identifier": "unit-run",
                "field_line_selection": "unit-line",
                "coupling_cadence_s": 60,
            },
            "inputs": manifest["inputs"],
            "random_seeds": [],
            "compiler": manifest["compiler"],
            "output_schema": manifest["output_schema"],
            "acceptance_metrics": [
                {
                    "name": "unit background replay",
                    "metric_family": "background_replay",
                    "value": 0,
                    "tolerance": 1,
                    "comparison": "<=",
                    "units": "dimensionless",
                },
                {
                    "name": "unit transport response",
                    "metric_family": "transport_response",
                    "value": 0,
                    "tolerance": 1,
                    "comparison": "<=",
                    "units": "dimensionless",
                },
            ],
        }
        swmf_path = temp_dir / "swmf.json"
        swmf_path.write_text(json.dumps(swmf_manifest), encoding="utf-8")
        check_external_gate.validate_gate(SimpleNamespace(
            evidence_class="COUPLED_INTEGRATION", manifest=[str(swmf_path)]))

        # A pre-existing destination is replaced only after a complete JSON
        # payload is ready; this checks the public helper's commit semantics.
        destination = temp_dir / "transactional.txt"
        destination.write_text("old", encoding="utf-8")
        run_campaign.write_transactional(destination, "new\n")
        require(destination.read_text(encoding="utf-8") == "new\n",
                "transactional writer did not publish the complete payload")

    print("PASS VAL-GOVERNANCE: evidence classes, external gates, and label barriers enforced")
    return 0


if __name__ == "__main__":
    sys.exit(main())

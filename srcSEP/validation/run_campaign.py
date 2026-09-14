#!/usr/bin/env python3
"""Assemble auditable Step 15 evidence without inflating scientific claims.

The source-only C++ runners provide numerical and actual SWCME integration
records.  Real SWMF outputs and spacecraft products necessarily live outside a
source distribution, so this program accepts explicit manifests for them,
checks every local input checksum, and keeps missing evidence as INCOMPLETE.
All output is written transactionally: a complete temporary file is flushed
and fsynced before os.replace publishes it at the requested destination.
"""

import argparse
import datetime
import hashlib
import json
import os
import pathlib
import platform
import re
import subprocess
import sys
import tempfile
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


SCHEMA = "srcsep-validation-campaign-v1"
EXTERNAL_SCHEMA = "srcsep-external-evidence-v1"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
FORBIDDEN_OBSERVATIONAL_LABELS = ("synthetic", "fixture", "mock")
OBSERVATIONAL_METRIC_FAMILIES = {
    "onset",
    "anisotropy",
    "spectra",
    "fluence",
    "decay",
    "multi_spacecraft_longitude",
}
COUPLED_METRIC_FAMILIES = {"background_replay", "transport_response"}


class EvidenceError(Exception):
    """Raised when supplied evidence is malformed or its bytes do not match."""


def sha256_file(path: pathlib.Path) -> str:
    """Hash a file in bounded chunks so large model products remain practical."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_tree_sha256(root: pathlib.Path) -> str:
    """Return a deterministic digest of source-controlled-style input bytes.

    Git metadata is unavailable in some delivered archives.  This digest uses
    relative path, file length, and content for every relevant source file,
    while excluding generated reports, caches, object files, and archives.
    Thus an archive-only campaign still has a stable, independently repeatable
    identity without accidentally hashing its own output.
    """
    ignored_parts = {".git", "__pycache__", "validation-output", "test_output"}
    ignored_suffixes = {
        ".a", ".o", ".pyc", ".tar", ".gz", ".zip", ".json.tmp", ".xml.tmp"
    }
    digest = hashlib.sha256()
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        relative = path.relative_to(root)
        if any(part in ignored_parts for part in relative.parts):
            continue
        if any(path.name.endswith(suffix) for suffix in ignored_suffixes):
            continue
        encoded_name = relative.as_posix().encode("utf-8")
        digest.update(len(encoded_name).to_bytes(8, "big"))
        digest.update(encoded_name)
        size = path.stat().st_size
        digest.update(size.to_bytes(8, "big"))
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def load_json(path: pathlib.Path) -> Dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise EvidenceError("cannot read JSON {}: {}".format(path, error))
    if not isinstance(value, dict):
        raise EvidenceError("{} must contain one JSON object".format(path))
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def validate_internal_report(path: pathlib.Path,
                             expected_ids: Sequence[str]) -> Dict[str, Any]:
    """Validate the complete registry record and return cases by stable ID."""
    report = load_json(path)
    require(report.get("schema") == "srcsep-component-tests-v1",
            "{} has an unsupported registry schema".format(path))
    results = report.get("results")
    require(isinstance(results, list), "{} omits results[]".format(path))
    by_id = {item.get("id"): item for item in results if isinstance(item, dict)}
    require(set(by_id) == set(expected_ids),
            "{} must contain exactly {}".format(path, ", ".join(expected_ids)))
    for case_id in expected_ids:
        case = by_id[case_id]
        require(case.get("status") in ("PASS", "FAIL", "SKIP", "ERROR"),
                "{} has invalid status".format(case_id))
        require(isinstance(case.get("configuration"), list) and
                bool(case["configuration"]),
                "{} must record effective configuration".format(case_id))
        require(isinstance(case.get("metrics"), list) and bool(case["metrics"]),
                "{} must record acceptance metrics".format(case_id))
        require("seed" in case, "{} must record its deterministic seed".format(case_id))
        for metric in case["metrics"]:
            require(all(key in metric for key in
                        ("name", "value", "tolerance", "comparison", "units")),
                    "{} has an incomplete metric".format(case_id))
    return by_id


def require_nonempty_mapping(value: Any, label: str) -> Dict[str, Any]:
    require(isinstance(value, dict) and bool(value), "{} must be a nonempty object".format(label))
    return value


def require_nonempty_string(value: Any, label: str) -> str:
    require(isinstance(value, str) and bool(value.strip()), "{} must be a nonempty string".format(label))
    return value


def metric_satisfies_acceptance(metric: Dict[str, Any], label: str) -> bool:
    """Evaluate the declared scalar acceptance rule without using truthiness."""
    value = metric.get("value")
    tolerance = metric.get("tolerance")
    comparison = metric.get("comparison")
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            "{}.value must be numeric".format(label))
    require(isinstance(tolerance, (int, float)) and not isinstance(tolerance, bool),
            "{}.tolerance must be numeric".format(label))
    operations = {
        "<=": lambda: value <= tolerance,
        "<": lambda: value < tolerance,
        ">=": lambda: value >= tolerance,
        ">": lambda: value > tolerance,
        "==": lambda: value == tolerance,
    }
    require(comparison in operations,
            "{}.comparison must be one of <=, <, >=, >, ==".format(label))
    return operations[comparison]()


def validate_external_manifest(path: pathlib.Path,
                               evidence_class: str) -> Tuple[Dict[str, Any], List[Dict[str, str]]]:
    """Authenticate an external record and every file on which it depends."""
    manifest = load_json(path)
    require(manifest.get("schema") == EXTERNAL_SCHEMA,
            "{} has an unsupported external-evidence schema".format(path))
    require(manifest.get("evidence_class") == evidence_class,
            "{} has evidence_class {}, expected {}".format(
                path, manifest.get("evidence_class"), evidence_class))
    expected_data_class = (
        "SWMF_OUTPUT" if evidence_class == "COUPLED_INTEGRATION"
        else "SPACECRAFT_OBSERVATION"
    )
    require(manifest.get("data_class") == expected_data_class,
            "{} must declare data_class={}".format(path, expected_data_class))
    require(manifest.get("status") in ("PASS", "FAIL"),
            "{} is not completed evidence (status must be PASS or FAIL)".format(path))
    require_nonempty_string(manifest.get("campaign_id"), "campaign_id")
    require_nonempty_mapping(manifest.get("configuration"), "configuration")
    require(isinstance(manifest.get("random_seeds"), list),
            "random_seeds must be an array (empty is valid for deterministic runs)")

    compiler = require_nonempty_mapping(manifest.get("compiler"), "compiler")
    for field in ("command", "version"):
        require_nonempty_string(compiler.get(field), "compiler.{}".format(field))
    require(isinstance(compiler.get("flags"), list), "compiler.flags must be an array")
    output_schema = require_nonempty_mapping(manifest.get("output_schema"), "output_schema")
    for field in ("name", "version"):
        require_nonempty_string(output_schema.get(field), "output_schema.{}".format(field))
    require(isinstance(output_schema.get("variables"), list) and
            bool(output_schema["variables"]),
            "output_schema.variables must be a nonempty array")
    for index, variable in enumerate(output_schema["variables"]):
        require(isinstance(variable, dict),
                "output variable {} must be an object".format(index))
        require_nonempty_string(variable.get("name"),
                                "output_schema.variables[{}].name".format(index))
        require_nonempty_string(variable.get("units"),
                                "output_schema.variables[{}].units".format(index))

    metrics = manifest.get("acceptance_metrics")
    require(isinstance(metrics, list) and bool(metrics),
            "acceptance_metrics must be a nonempty array")
    for index, metric in enumerate(metrics):
        require(isinstance(metric, dict), "acceptance metric {} must be an object".format(index))
        for field in ("name", "metric_family", "value", "tolerance", "comparison", "units"):
            require(field in metric, "acceptance metric {} omits {}".format(index, field))
        require_nonempty_string(metric.get("name"),
                                "acceptance_metrics[{}].name".format(index))
        require_nonempty_string(metric.get("metric_family"),
                                "acceptance_metrics[{}].metric_family".format(index))
        require_nonempty_string(metric.get("units"),
                                "acceptance_metrics[{}].units".format(index))
    metric_results = [
        metric_satisfies_acceptance(metric, "acceptance_metrics[{}]".format(index))
        for index, metric in enumerate(metrics)
    ]
    require((manifest["status"] == "PASS") == all(metric_results),
            "manifest status must equal the combined metric acceptance result")

    inputs = manifest.get("inputs")
    require(isinstance(inputs, list) and bool(inputs), "inputs must be a nonempty array")
    checksums: List[Dict[str, str]] = []
    for index, item in enumerate(inputs):
        require(isinstance(item, dict), "input {} must be an object".format(index))
        for field in ("path", "sha256", "role", "source_url_or_pid",
                      "access_utc", "license_or_acknowledgment"):
            require_nonempty_string(item.get(field), "inputs[{}].{}".format(index, field))
        expected_sha = item["sha256"].lower()
        require(bool(SHA256_RE.match(expected_sha)),
                "inputs[{}].sha256 is not lowercase SHA-256".format(index))
        input_path = (path.parent / item["path"]).resolve()
        require(input_path.is_file(), "external input does not exist: {}".format(input_path))
        actual_sha = sha256_file(input_path)
        require(actual_sha == expected_sha,
                "checksum mismatch for {}: expected {}, got {}".format(
                    input_path, expected_sha, actual_sha))
        checksums.append({"path": str(input_path), "sha256": actual_sha})

    if evidence_class == "COUPLED_INTEGRATION":
        configuration = manifest["configuration"]
        for field in ("srcsep_parameter_file", "swmf_run_identifier",
                      "field_line_selection"):
            require_nonempty_string(configuration.get(field),
                                    "configuration.{}".format(field))
        cadence = configuration.get("coupling_cadence_s")
        require(isinstance(cadence, (int, float)) and cadence > 0,
                "configuration.coupling_cadence_s must be positive")
        coupled_families = {metric["metric_family"] for metric in metrics}
        require(COUPLED_METRIC_FAMILIES.issubset(coupled_families),
                "SWMF evidence must include background_replay and transport_response metrics")

    if evidence_class == "OBSERVATIONAL_VALIDATION":
        # A second, deliberately simple label barrier prevents a human-edited
        # record from assigning observational authority to manufactured data.
        serialized = json.dumps(manifest, sort_keys=True).lower()
        for forbidden in FORBIDDEN_OBSERVATIONAL_LABELS:
            require(forbidden not in serialized,
                    "observational evidence contains forbidden label '{}'".format(forbidden))
        configuration = require_nonempty_mapping(
            manifest.get("configuration"), "configuration")
        # WP39 requires the evidence manifest to describe the same operations
        # performed by sep_observation_forward_model: response folding,
        # exposure/cadence, species, detector nonlinearity, background, and
        # uncertainty propagation.  A free-form label cannot establish that a
        # simulation spectrum was converted into instrument-space counts.
        forward = require_nonempty_mapping(
            configuration.get("forward_operator"),
            "configuration.forward_operator")
        for field in ("version", "energy_response", "angular_response",
                      "cadence_s", "species", "dead_time_s",
                      "saturation_policy", "background_subtraction",
                      "uncertainty_propagation"):
            require(field in forward,
                    "configuration.forward_operator omits {}".format(field))
        require_nonempty_string(forward.get("version"),
                                "configuration.forward_operator.version")
        require_nonempty_string(forward.get("energy_response"),
                                "configuration.forward_operator.energy_response")
        require_nonempty_string(forward.get("angular_response"),
                                "configuration.forward_operator.angular_response")
        require_nonempty_string(forward.get("species"),
                                "configuration.forward_operator.species")
        require(isinstance(forward.get("cadence_s"), (int, float)) and
                forward["cadence_s"] > 0,
                "configuration.forward_operator.cadence_s must be positive")
        require(isinstance(forward.get("dead_time_s"), (int, float)) and
                forward["dead_time_s"] >= 0,
                "configuration.forward_operator.dead_time_s must be non-negative")
        for field in ("saturation_policy", "background_subtraction",
                      "uncertainty_propagation"):
            require_nonempty_string(
                forward.get(field),
                "configuration.forward_operator.{}".format(field))
        event = require_nonempty_mapping(manifest.get("event"), "event")
        for field in ("id", "start_utc", "end_utc"):
            require_nonempty_string(event.get(field), "event.{}".format(field))
        observers = manifest.get("observers")
        require(isinstance(observers, list) and bool(observers),
                "observers must be a nonempty array")
        for index, observer in enumerate(observers):
            require(isinstance(observer, dict), "observer {} must be an object".format(index))
            for field in ("mission", "instrument", "product_level", "product_version",
                          "variables", "cadence_s", "coordinate_system", "quality_flags"):
                require(field in observer, "observer {} omits {}".format(index, field))
        require(manifest.get("held_out") is True,
                "observational event must be held_out=true")
        uncertainty = require_nonempty_mapping(manifest.get("uncertainty"), "uncertainty")
        require_nonempty_string(uncertainty.get("method"), "uncertainty.method")
        confidence = uncertainty.get("confidence_level")
        require(isinstance(confidence, (int, float)) and 0.0 < confidence < 1.0,
                "uncertainty.confidence_level must lie strictly between 0 and 1")
        for index, metric in enumerate(metrics):
            require(metric["metric_family"] in OBSERVATIONAL_METRIC_FAMILIES,
                    "observational metric {} has unknown metric_family {}".format(
                        index, metric["metric_family"]))
        if any(metric["metric_family"] == "multi_spacecraft_longitude"
               for metric in metrics):
            require(len(observers) >= 2,
                    "multi_spacecraft_longitude requires at least two observers")

    checksums.append({"path": str(path.resolve()), "sha256": sha256_file(path)})
    return manifest, checksums


def section(status: str, required: Sequence[str], cases: Sequence[Dict[str, Any]],
            limitations: Sequence[str]) -> Dict[str, Any]:
    return {
        "status": status,
        "required_cases": list(required),
        "cases": list(cases),
        "limitations": list(limitations),
    }


def case_status(case: Dict[str, Any]) -> str:
    return "PASS" if case.get("status") == "PASS" else "FAIL"


def compiler_provenance() -> Dict[str, Any]:
    command = os.environ.get("CXX", "c++")
    try:
        completed = subprocess.run(
            [command, "--version"], check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        version = completed.stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        version = "unavailable: {}".format(error)
    return {
        "command": command,
        "version": version,
        "flags": [
            "-std=c++11 (VAL01-VAL03)",
            "-std=c++17 (VAL04)",
            "-O1", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            "-fsanitize=address,undefined", "-fno-omit-frame-pointer",
        ],
    }


def write_transactional(path: pathlib.Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def markdown_report(report: Dict[str, Any]) -> str:
    lines = [
        "# srcSEP Step 15 validation campaign",
        "",
        "Release status: **{}**".format(report["release_status"]),
        "",
        "| Evidence class | Status | Cases |",
        "|---|---:|---|",
    ]
    labels = {
        "numerical_verification": "Numerical verification",
        "cross_mover_verification": "Cross-mover verification",
        "cross_model_verification": "Cross-model verification",
        "coupled_integration": "Coupled integration",
        "observational_validation": "Observational validation",
    }
    for key, value in report["sections"].items():
        case_ids = ", ".join(str(case.get("id", case.get("campaign_id", "external")))
                             for case in value["cases"]) or "none"
        lines.append("| {} | {} | {} |".format(labels[key], value["status"], case_ids))
    lines.extend(["", "## Limitations", ""])
    limitations = []
    for value in report["sections"].values():
        limitations.extend(value["limitations"])
    if limitations:
        lines.extend("- " + item for item in limitations)
    else:
        lines.append("- None recorded.")
    lines.extend([
        "",
        "The JSON companion is the authoritative machine-readable record and "
        "contains complete configurations, seeds, metrics, tolerances, compiler "
        "identity, output schema, and input checksums.",
        "",
    ])
    return "\n".join(lines)


def assemble(args: argparse.Namespace) -> Dict[str, Any]:
    numerical_path = pathlib.Path(args.numerical_json).resolve()
    swcme_path = pathlib.Path(args.swcme_json).resolve()
    source_root = pathlib.Path(args.source_root).resolve()
    numerical = validate_internal_report(numerical_path, ("VAL01", "VAL02", "VAL03"))
    swcme = validate_internal_report(swcme_path, ("VAL04",))

    checksums = [
        {"path": str(numerical_path), "sha256": sha256_file(numerical_path)},
        {"path": str(swcme_path), "sha256": sha256_file(swcme_path)},
    ]
    sections: Dict[str, Dict[str, Any]] = {}
    sections["numerical_verification"] = section(
        case_status(numerical["VAL01"]), ("VAL01",), (numerical["VAL01"],), ())
    sections["cross_mover_verification"] = section(
        case_status(numerical["VAL02"]), ("VAL02",), (numerical["VAL02"],), ())
    sections["cross_model_verification"] = section(
        case_status(numerical["VAL03"]), ("VAL03",), (numerical["VAL03"],),
        ("VAL03 uses an independent finite-volume solver under an identical "
         "one-zone coefficient history; an additional external solver may be "
         "recorded in a future campaign.",))

    coupled_cases: List[Dict[str, Any]] = [swcme["VAL04"]]
    coupled_limitations: List[str] = []
    swmf_status = "INCOMPLETE"
    if args.swmf_manifest:
        swmf, external_checksums = validate_external_manifest(
            pathlib.Path(args.swmf_manifest).resolve(), "COUPLED_INTEGRATION")
        checksums.extend(external_checksums)
        swmf = dict(swmf)
        swmf["id"] = "VAL04-SWMF"
        coupled_cases.append(swmf)
        swmf_status = swmf["status"]
    else:
        coupled_limitations.append(
            "No checksum-verified real SWMF-to-srcSEP replay manifest was supplied.")
    swcme_status = case_status(swcme["VAL04"])
    coupled_status = (
        "FAIL" if "FAIL" in (swcme_status, swmf_status)
        else "PASS" if (swcme_status, swmf_status) == ("PASS", "PASS")
        else "INCOMPLETE"
    )
    sections["coupled_integration"] = section(
        coupled_status, ("VAL04-SWCME", "VAL04-SWMF"),
        coupled_cases, coupled_limitations)

    observation_cases: List[Dict[str, Any]] = []
    observation_limitations: List[str] = []
    observed_families = set()
    observation_failed = False
    for manifest_name in args.observational_manifest:
        observation, external_checksums = validate_external_manifest(
            pathlib.Path(manifest_name).resolve(), "OBSERVATIONAL_VALIDATION")
        checksums.extend(external_checksums)
        observation_cases.append(observation)
        observation_failed = observation_failed or observation["status"] == "FAIL"
        observed_families.update(metric["metric_family"]
                                 for metric in observation["acceptance_metrics"])
    missing_families = sorted(OBSERVATIONAL_METRIC_FAMILIES - observed_families)
    if not observation_cases:
        observation_limitations.append(
            "No held-out, checksum-verified spacecraft event manifest was supplied.")
    if missing_families:
        observation_limitations.append(
            "Observational metrics not closed: {}.".format(", ".join(missing_families)))
    observation_status = (
        "FAIL" if observation_failed
        else "PASS" if observation_cases and not missing_families
        else "INCOMPLETE"
    )
    sections["observational_validation"] = section(
        observation_status,
        tuple(sorted(OBSERVATIONAL_METRIC_FAMILIES)),
        observation_cases, observation_limitations)

    statuses = [value["status"] for value in sections.values()]
    release_status = (
        "FAILED" if "FAIL" in statuses
        else "RELEASE_READY" if all(value == "PASS" for value in statuses)
        else "INCOMPLETE"
    )
    now = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)
    return {
        "schema": SCHEMA,
        "campaign_id": args.campaign_id,
        "generated_utc": now.isoformat().replace("+00:00", "Z"),
        "release_status": release_status,
        "provenance": {
            "source_tree_sha256": source_tree_sha256(source_root),
            "input_checksums": checksums,
            "compiler": compiler_provenance(),
            "platform": platform.platform(),
            "output_schema": SCHEMA,
        },
        "sections": sections,
    }


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numerical-json", required=True)
    parser.add_argument("--swcme-json", required=True)
    parser.add_argument("--source-root", default=str(pathlib.Path(__file__).resolve().parents[1]))
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--campaign-id", default="srcsep-step15-local")
    parser.add_argument("--swmf-manifest")
    parser.add_argument("--observational-manifest", action="append", default=[])
    parser.add_argument(
        "--release", action="store_true",
        help="return nonzero unless every evidence class is PASS")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_arguments(argv)
    try:
        report = assemble(args)
        output_dir = pathlib.Path(args.output_dir).resolve()
        write_transactional(
            output_dir / "step15-validation-campaign.json",
            json.dumps(report, indent=2, sort_keys=True) + "\n")
        write_transactional(
            output_dir / "step15-validation-campaign.md", markdown_report(report))
    except EvidenceError as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        return 2
    print("Step 15 validation campaign: {}".format(report["release_status"]))
    if args.release and report["release_status"] != "RELEASE_READY":
        print("Release gate refused incomplete or failed evidence.", file=sys.stderr)
        return 1
    return 0 if report["release_status"] != "FAILED" else 1


if __name__ == "__main__":
    sys.exit(main())

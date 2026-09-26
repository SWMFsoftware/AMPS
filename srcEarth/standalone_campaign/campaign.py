#!/usr/bin/env python3
"""Strict manifest, immutable-cache, rendering, and evidence helpers.

The functions in this module are intentionally dependency-free so that a
campaign can be preflighted on a login node before an AMPS allocation starts.
Scientific files are addressed by both a local path and SHA-256 digest.  A URL
is rejected: acquisition is a separate, provenance-recorded step and mutable
network content is never fetched while a validation score is being produced.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple


UTC = timezone.utc
SCHEMA_VERSION = "earth-standalone-campaign/v1"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
PLACEHOLDER_RE = re.compile(r"\{\{([A-Z][A-Z0-9_]*)\}\}")

# Every Step-8 campaign must state where each class of scientific input came
# from.  An analytic smoke case may use a small descriptor saying that a class
# is not physically used, but it may not silently omit the class.
REQUIRED_RESOURCE_KINDS = {
    "input_template",
    "driver",
    "boundary_spectrum",
    "ephemeris",
    "attitude",
    "instrument_response",
    "observation",
}

# Reject misspelled or ambiguous resource classes.  Extra campaign material is
# allowed only when it is deliberately classified as auxiliary evidence; this
# keeps a typo such as ``ephemeries`` from satisfying a human review while the
# runner silently uses a different file.
ALLOWED_RESOURCE_KINDS = REQUIRED_RESOURCE_KINDS | {
    "auxiliary",
    "calibration",
    "quality_mask",
}

ALLOWED_PRODUCTS = {
    "cutoff",
    "directional_access",
    "spectrum",
    "density",
    "flux",
    "detector_rate",
}

ALLOWED_ADAPTERS = {
    "GOES_EPEAD_DIRECTIONAL",
    "PAMELA_CUTOFF",
    "POES_METOP_MEPED_CUTOFF",
    "REPT_PROTON_SPECTRUM",
}

ALLOWED_OPERATORS = {"<=", "<", ">=", ">", "=="}


class CampaignError(RuntimeError):
    """Raised when campaign provenance or evidence is incomplete or invalid."""


def parse_utc(value: str) -> datetime:
    """Parse an ISO-8601 timestamp and require an explicit UTC offset.

    Naive timestamps are rejected because treating local time as UTC can shift a
    field snapshot relative to the observation by hours while leaving otherwise
    plausible output.
    """

    text = str(value).strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        result = datetime.fromisoformat(text)
    except ValueError as exc:
        raise CampaignError("invalid ISO-8601 timestamp: %s" % value) from exc
    if result.tzinfo is None:
        raise CampaignError("timestamp lacks an explicit UTC offset: %s" % value)
    result = result.astimezone(UTC)
    return result


def format_utc(value: datetime) -> str:
    """Return the one canonical UTC representation used in paths and records."""

    return value.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def epoch_token(value: datetime) -> str:
    """Return a filesystem-safe UTC token with second resolution."""

    return value.astimezone(UTC).strftime("%Y%m%dT%H%M%SZ")


def sha256_file(path: Path) -> str:
    """Hash a file as a stream so large observation tables need little memory."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    """Return a stable digest for rendered inputs and run fingerprints."""

    return hashlib.sha256(payload).hexdigest()


def atomic_write_json(path: Path, value: object) -> None:
    """Replace a JSON status file only after the complete document is durable.

    A killed batch job therefore leaves either the previous valid record or no
    record, never a truncated file that a restart could mistake for PASS.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, str(path))
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass


def _require_mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise CampaignError("%s must be a JSON object" % label)
    return value


def _require_sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, list):
        raise CampaignError("%s must be a JSON array" % label)
    return value


def _require_id(value: object, label: str) -> str:
    text = str(value)
    if not ID_RE.fullmatch(text):
        raise CampaignError("%s is not a portable identifier: %s" % (label, text))
    return text


def _resource_path(manifest_root: Path, text: str) -> Path:
    """Resolve one local resource and reject every network-style locator."""

    if "://" in text or text.startswith("data:"):
        raise CampaignError(
            "network/mutable resource is forbidden during scoring: %s" % text
        )
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = manifest_root / path
    return path.resolve()


def _require_safe_relative_path(text: str, label: str) -> str:
    """Reject absolute and parent-traversing run-directory paths/patterns."""

    path = Path(text)
    if not text or path.is_absolute() or ".." in path.parts:
        raise CampaignError("%s must stay inside the per-run directory" % label)
    return text


def _validate_resource_records(
    manifest: Mapping[str, object], manifest_root: Path
) -> List[Dict[str, object]]:
    resources_value = _require_sequence(manifest.get("resources"), "resources")
    if not resources_value:
        raise CampaignError("resources must not be empty")

    resources: List[Dict[str, object]] = []
    identifiers = set()
    kinds = set()
    for index, raw in enumerate(resources_value):
        item = _require_mapping(raw, "resources[%d]" % index)
        identifier = _require_id(item.get("id"), "resources[%d].id" % index)
        if identifier in identifiers:
            raise CampaignError("duplicate resource id: %s" % identifier)
        identifiers.add(identifier)
        kind = str(item.get("kind", "")).strip().lower()
        if kind not in ALLOWED_RESOURCE_KINDS:
            raise CampaignError(
                "resource %s has unsupported kind: %s" % (identifier, kind)
            )
        kinds.add(kind)
        provenance = str(item.get("provenance", "")).strip()
        if not provenance:
            raise CampaignError(
                "resource %s has no provenance statement" % identifier
            )
        digest = str(item.get("sha256", "")).strip().lower()
        if not SHA256_RE.fullmatch(digest):
            raise CampaignError("resource %s has invalid SHA-256" % identifier)
        source = _resource_path(manifest_root, str(item.get("path", "")))
        if not source.is_file():
            raise CampaignError("resource %s is missing: %s" % (identifier, source))
        actual = sha256_file(source)
        if actual != digest:
            raise CampaignError(
                "resource %s hash mismatch: expected %s, got %s"
                % (identifier, digest, actual)
            )
        resources.append(
            {
                "id": identifier,
                "kind": kind,
                "source_path": str(source),
                "sha256": actual,
                "size_bytes": source.stat().st_size,
                "provenance": provenance,
            }
        )

    missing = sorted(REQUIRED_RESOURCE_KINDS.difference(kinds))
    if missing:
        raise CampaignError("campaign omits required resource kinds: %s" % ", ".join(missing))
    return resources


def campaign_epochs(event: Mapping[str, object]) -> List[datetime]:
    """Construct an inclusive and exactly aligned campaign cadence.

    A partial final interval is rejected rather than silently shortening one
    exposure.  Observation adapters retain their native windows separately.
    """

    start = parse_utc(str(event.get("start_utc", "")))
    end = parse_utc(str(event.get("end_utc", "")))
    if end < start:
        raise CampaignError("event end_utc precedes start_utc")
    try:
        cadence = int(event.get("field_cadence_seconds", 0))
    except (TypeError, ValueError) as exc:
        raise CampaignError("field_cadence_seconds must be an integer") from exc
    if cadence <= 0:
        raise CampaignError("field_cadence_seconds must be positive")
    duration = int((end - start).total_seconds())
    if duration % cadence:
        raise CampaignError(
            "event interval is not an integer number of field cadences"
        )
    return [start + timedelta(seconds=offset) for offset in range(0, duration + 1, cadence)]


def _validate_execution(
    manifest: Mapping[str, object], resource_ids: Sequence[str]
) -> None:
    execution = _require_mapping(manifest.get("execution"), "execution")
    template_id = str(execution.get("input_template_resource", ""))
    if template_id not in resource_ids:
        raise CampaignError("execution input_template_resource is not a resource id")
    command = _require_sequence(execution.get("command"), "execution.command")
    if not command or not all(isinstance(token, str) and token for token in command):
        raise CampaignError("execution.command must contain nonempty string tokens")
    if not any("{input}" in token for token in command):
        raise CampaignError("execution.command must contain {input}")
    if not any("{amps}" in token for token in command):
        raise CampaignError("execution.command must contain {amps}")
    expected = _require_sequence(
        execution.get("expected_artifacts"), "execution.expected_artifacts"
    )
    if not expected or not all(isinstance(pattern, str) and pattern for pattern in expected):
        raise CampaignError("execution.expected_artifacts must not be empty")
    for pattern in expected:
        _require_safe_relative_path(str(pattern), "expected artifact pattern")
    prediction = execution.get("prediction_artifact")
    if prediction is not None and (not isinstance(prediction, str) or not prediction):
        raise CampaignError("execution.prediction_artifact must be a nonempty path")
    if prediction is not None:
        _require_safe_relative_path(str(prediction), "prediction_artifact")
    extractors = execution.get("prediction_extractors", [])
    if not isinstance(extractors, list):
        raise CampaignError("execution.prediction_extractors must be a JSON array")
    for index, raw in enumerate(extractors):
        extractor = _require_mapping(raw, "execution.prediction_extractors[%d]" % index)
        operation = str(extractor.get("operation", "COLUMN")).upper()
        if operation not in ("COLUMN", "RATIO"):
            raise CampaignError("prediction extractor has unsupported operation: %s" % operation)
        for key in (
            "instrument_id", "artifact", "platform", "channel", "direction",
            "quantity", "units",
        ):
            if not str(extractor.get(key, "")):
                raise CampaignError("prediction extractor %d lacks %s" % (index, key))
        _require_safe_relative_path(str(extractor["artifact"]), "prediction extractor artifact")
        expected_variables = _require_sequence(
            extractor.get("expected_variables"),
            "execution.prediction_extractors[%d].expected_variables" % index,
        )
        if not expected_variables or not all(
            isinstance(name, str) and name for name in expected_variables
        ):
            raise CampaignError(
                "prediction extractor %d requires the exact expected_variables list"
                % index
            )
        try:
            expected_rows = int(extractor.get("expected_row_count", 0))
            row_index = int(extractor.get("row_index", 0))
        except (TypeError, ValueError) as exc:
            raise CampaignError(
                "prediction extractor %d row counts must be integers" % index
            ) from exc
        if expected_rows <= 0 or row_index < 0 or row_index >= expected_rows:
            raise CampaignError(
                "prediction extractor %d has an invalid expected row count/index"
                % index
            )
        if operation == "COLUMN" and not str(extractor.get("column", "")):
            raise CampaignError("COLUMN prediction extractor lacks column")
        if operation == "RATIO" and (
            not str(extractor.get("numerator_column", ""))
            or not str(extractor.get("denominator_column", ""))
        ):
            raise CampaignError("RATIO prediction extractor requires numerator/denominator columns")


def _validate_instruments(
    manifest: Mapping[str, object], resource_kinds: Mapping[str, str], role: str
) -> bool:
    """Validate observation contracts and report whether scoring is required.

    Adapter quantity/unit pairs are part of the scientific interface.  Checking
    them here prevents a misspelled quantity from surviving preflight and later
    appearing merely as a missing comparison.  REPT unit strings vary with the
    released export, so that adapter requires an explicit nonempty unit and then
    checks every row against it in :mod:`adapters`.
    """

    instruments = _require_sequence(manifest.get("instruments"), "instruments")
    if not instruments:
        raise CampaignError("instruments must not be empty")
    seen = set()
    required_comparisons = 0
    adapter_contracts = {
        "PAMELA_CUTOFF": ("cutoff_latitude", "deg_AACGM"),
        "POES_METOP_MEPED_CUTOFF": ("cutoff_latitude", "deg_AACGM"),
        "GOES_EPEAD_DIRECTIONAL": ("directional_ratio", "1"),
        "REPT_PROTON_SPECTRUM": ("differential_flux", None),
    }
    for index, raw in enumerate(instruments):
        item = _require_mapping(raw, "instruments[%d]" % index)
        identifier = _require_id(item.get("id"), "instruments[%d].id" % index)
        if identifier in seen:
            raise CampaignError("duplicate instrument id: %s" % identifier)
        seen.add(identifier)
        adapter = str(item.get("adapter", "")).upper()
        if adapter not in ALLOWED_ADAPTERS:
            raise CampaignError("unsupported observation adapter: %s" % adapter)
        for key in ("observation_resource", "response_resource"):
            if str(item.get(key, "")) not in resource_kinds:
                raise CampaignError("instrument %s has invalid %s" % (identifier, key))
        if resource_kinds[str(item["observation_resource"])] != "observation":
            raise CampaignError(
                "instrument %s observation_resource is not an observation" % identifier
            )
        if resource_kinds[str(item["response_resource"])] not in (
            "instrument_response", "calibration"
        ):
            raise CampaignError(
                "instrument %s response_resource is not a response/calibration"
                % identifier
            )
        if not isinstance(item.get("comparison_required", True), bool):
            raise CampaignError(
                "instrument %s comparison_required must be boolean" % identifier
            )
        comparison_required = bool(item.get("comparison_required", True))
        required_comparisons += int(comparison_required)
        quantity = str(item.get("quantity", "")).strip()
        units = str(item.get("units", "")).strip()
        if comparison_required and (not quantity or not units):
            raise CampaignError(
                "instrument %s requires comparison quantity and units" % identifier
            )
        expected_quantity, expected_units = adapter_contracts[adapter]
        if quantity != expected_quantity:
            raise CampaignError(
                "instrument %s quantity %s does not match adapter quantity %s"
                % (identifier, quantity or "<empty>", expected_quantity)
            )
        if expected_units is not None and units != expected_units:
            raise CampaignError(
                "instrument %s units %s do not match adapter units %s"
                % (identifier, units or "<empty>", expected_units)
            )
        try:
            cadence = int(item.get("cadence_seconds", 0))
        except (TypeError, ValueError) as exc:
            raise CampaignError(
                "instrument %s cadence_seconds must be an integer" % identifier
            ) from exc
        if cadence <= 0:
            raise CampaignError(
                "instrument %s cadence_seconds must be positive" % identifier
            )
    if role in ("VALIDATION", "HOLDOUT") and required_comparisons == 0:
        raise CampaignError(
            "%s campaign must require at least one observation comparison" % role
        )
    return required_comparisons > 0


def _validate_gates(manifest: Mapping[str, object], role: str) -> None:
    validation = _require_mapping(manifest.get("validation"), "validation")
    gates = _require_sequence(validation.get("numerical_gates"), "validation.numerical_gates")
    if not gates:
        raise CampaignError("at least one numerical gate is required")
    identifiers = set()
    for index, raw in enumerate(gates):
        gate = _require_mapping(raw, "validation.numerical_gates[%d]" % index)
        identifier = _require_id(gate.get("id"), "numerical gate id")
        if identifier in identifiers:
            raise CampaignError("duplicate numerical gate id: %s" % identifier)
        identifiers.add(identifier)
        if str(gate.get("operator", "")) not in ALLOWED_OPERATORS:
            raise CampaignError("gate %s has unsupported operator" % identifier)
        try:
            threshold = float(gate.get("threshold"))
        except (TypeError, ValueError) as exc:
            raise CampaignError("gate %s threshold is not numeric" % identifier) from exc
        if not math.isfinite(threshold):
            raise CampaignError("gate %s threshold is not finite" % identifier)
        source = str(gate.get("source", ""))
        if source not in ("termination_summary", "json_artifact"):
            raise CampaignError("gate %s has unsupported source" % identifier)
        if source == "json_artifact":
            if not str(gate.get("path", "")) or not str(gate.get("metric", "")):
                raise CampaignError("JSON gate %s requires path and metric" % identifier)
            _require_safe_relative_path(str(gate["path"]), "JSON gate path")

    by_id = {str(item["id"]): item for item in gates}

    def require_not_weaker(
        identifier: str, source: str, metric: str, maximum: float
    ) -> None:
        gate = by_id.get(identifier)
        if gate is None:
            raise CampaignError("required numerical gate is missing: %s" % identifier)
        if (
            gate.get("source") != source
            or str(gate.get("metric", "")) != metric
            or str(gate.get("operator", "")) != "<="
            or not bool(gate.get("required", True))
            or float(gate["threshold"]) > maximum
        ):
            raise CampaignError(
                "gate %s weakens or changes the frozen %s <= %.12g contract"
                % (identifier, metric, maximum)
            )

    # Step 8 is orchestration, not a license to raise the Step-5 unresolved
    # limit.  Every campaign, including a smoke campaign, retains the 0.01
    # response-support gate.  Validation/holdout campaigns must additionally
    # retain the roadmap's two-percent energy and angular convergence gates.
    require_not_weaker(
        "unresolved_support", "termination_summary",
        "max_unresolved_fraction", 0.01,
    )
    if role in ("VALIDATION", "HOLDOUT"):
        require_not_weaker(
            "energy_convergence", "json_artifact",
            "energy_relative_change", 0.02,
        )
        require_not_weaker(
            "angular_convergence", "json_artifact",
            "angular_relative_change", 0.02,
        )


def _validate_exclusions(manifest: Mapping[str, object], role: str) -> None:
    """Require every exclusion to be explicit and preregistered for scoring."""

    exclusions = _require_sequence(manifest.get("exclusions"), "exclusions")
    seen = set()
    for index, raw in enumerate(exclusions):
        item = _require_mapping(raw, "exclusions[%d]" % index)
        identifier = _require_id(item.get("id"), "exclusions[%d].id" % index)
        if identifier in seen:
            raise CampaignError("duplicate exclusion id: %s" % identifier)
        seen.add(identifier)
        if not str(item.get("reason", "")).strip():
            raise CampaignError("exclusion %s has no reason" % identifier)
        if not isinstance(item.get("preregistered"), bool):
            raise CampaignError(
                "exclusion %s preregistered must be boolean" % identifier
            )
        if role in ("VALIDATION", "HOLDOUT") and not item["preregistered"]:
            raise CampaignError(
                "validation exclusion %s was not preregistered" % identifier
            )


def load_and_validate_manifest(path: Path) -> Dict[str, object]:
    """Load, structurally validate, and hash-check a Step-8 campaign manifest."""

    manifest_path = path.expanduser().resolve()
    try:
        with manifest_path.open(encoding="utf-8") as stream:
            raw = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise CampaignError("cannot read campaign manifest %s: %s" % (manifest_path, exc)) from exc
    manifest = dict(_require_mapping(raw, "manifest"))
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise CampaignError("unsupported schema_version: %s" % manifest.get("schema_version"))
    _require_id(manifest.get("campaign_id"), "campaign_id")
    _require_id(manifest.get("event_id"), "event_id")
    role = str(manifest.get("role", "")).upper()
    if role not in ("SMOKE", "DEVELOPMENT", "VALIDATION", "HOLDOUT"):
        raise CampaignError("unsupported campaign role: %s" % role)
    if not isinstance(manifest.get("frozen"), bool):
        raise CampaignError("frozen must be an explicit boolean")
    if role in ("VALIDATION", "HOLDOUT") and not manifest.get("frozen"):
        raise CampaignError("a VALIDATION or HOLDOUT campaign must be frozen")
    if manifest.get("normalization_policy") != "SHARED_EVENT_BOUNDARY_NO_PLATFORM_SCALE":
        raise CampaignError(
            "platform-specific normalization is forbidden; use the shared event boundary policy"
        )

    event = _require_mapping(manifest.get("event"), "event")
    campaign_epochs(event)
    products = _require_sequence(manifest.get("products"), "products")
    if not products:
        raise CampaignError("products must not be empty")
    unknown_products = sorted(set(str(item) for item in products).difference(ALLOWED_PRODUCTS))
    if unknown_products:
        raise CampaignError("unsupported products: %s" % ", ".join(unknown_products))
    _validate_exclusions(manifest, role)

    resources = _validate_resource_records(manifest, manifest_path.parent)
    resource_ids = [str(item["id"]) for item in resources]
    resource_kinds = {
        str(item["id"]): str(item["kind"]) for item in resources
    }
    _validate_execution(manifest, resource_ids)
    comparisons_required = _validate_instruments(manifest, resource_kinds, role)
    if comparisons_required and not manifest["execution"].get("prediction_artifact"):
        raise CampaignError(
            "required observation comparisons need execution.prediction_artifact"
        )
    _validate_gates(manifest, role)

    models = _require_sequence(manifest.get("models"), "models")
    if not models:
        raise CampaignError("models must not be empty")
    for index, model in enumerate(models):
        _require_id(model, "models[%d]" % index)

    # Store resolved, verified records under a private key used by the runner.
    # The original JSON structure remains intact for the frozen-manifest hash.
    manifest["_manifest_path"] = str(manifest_path)
    manifest["_manifest_sha256"] = sha256_file(manifest_path)
    manifest["_verified_resources"] = resources
    # Store the canonical spelling consumed by summaries/release checks while
    # retaining the exact file digest above for provenance.
    manifest["role"] = role
    return manifest


def cache_resources(
    verified_resources: Sequence[Mapping[str, object]], cache_root: Path
) -> Dict[str, Path]:
    """Copy references into a content-addressed cache and verify both copies."""

    cache_root.mkdir(parents=True, exist_ok=True)
    result: Dict[str, Path] = {}
    for item in verified_resources:
        source = Path(str(item["source_path"]))
        digest = str(item["sha256"])
        target_dir = cache_root / digest
        target = target_dir / source.name
        target_dir.mkdir(parents=True, exist_ok=True)
        if target.exists() and sha256_file(target) != digest:
            raise CampaignError("cached resource has wrong digest: %s" % target)
        if not target.exists():
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=target.name + ".", suffix=".copying", dir=str(target_dir)
            )
            os.close(descriptor)
            temporary = Path(temporary_name)
            try:
                shutil.copyfile(str(source), str(temporary))
                if sha256_file(temporary) != digest:
                    raise CampaignError(
                        "resource changed while being cached: %s" % source
                    )
                os.replace(str(temporary), str(target))
            finally:
                try:
                    temporary.unlink()
                except FileNotFoundError:
                    pass
        result[str(item["id"])] = target.resolve()
    return result


def render_template(template: str, values: Mapping[str, object]) -> str:
    """Render only explicit ``{{UPPER_CASE}}`` placeholders.

    Unknown placeholders fail.  This is deliberately smaller than a general
    template language: campaign inputs must be reviewable data, not executable
    code whose behavior depends on the rendering host.
    """

    normalized = {str(key): str(value) for key, value in values.items()}

    def replace(match: re.Match) -> str:
        key = match.group(1)
        if key not in normalized:
            raise CampaignError("input template contains unknown placeholder: %s" % key)
        return normalized[key]

    rendered = PLACEHOLDER_RE.sub(replace, template)
    remaining = PLACEHOLDER_RE.findall(rendered)
    if remaining:
        raise CampaignError("unresolved input placeholders: %s" % ", ".join(remaining))
    return rendered


def template_values(
    manifest: Mapping[str, object], cached: Mapping[str, Path], model: str, epoch: datetime,
    output_dir: Path
) -> Dict[str, object]:
    """Build deterministic renderer values, including every resource path."""

    values: Dict[str, object] = {
        "CAMPAIGN_ID": manifest["campaign_id"],
        "EVENT_ID": manifest["event_id"],
        "FIELD_MODEL": model,
        "EPOCH_UTC": format_utc(epoch),
        "EPOCH_TOKEN": epoch_token(epoch),
        "OUTPUT_DIR": str(output_dir.resolve()),
    }
    for identifier, path in cached.items():
        key = re.sub(r"[^A-Za-z0-9]", "_", identifier).upper() + "_PATH"
        values[key] = str(path)
    constants = _require_mapping(manifest.get("template_values", {}), "template_values")
    for key, value in constants.items():
        name = str(key)
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*", name):
            raise CampaignError("invalid template_values key: %s" % name)
        if name in values:
            raise CampaignError("template_values attempts to replace reserved key: %s" % name)
        values[name] = value
    return values


def parse_tecplot_table(path: Path) -> Tuple[List[str], List[List[float]]]:
    """Read numeric rows from the simple Tecplot POINT tables written by AMPS."""

    variables: List[str] = []
    rows: List[List[float]] = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            text = raw.strip()
            if not text:
                continue
            if text.upper().startswith("VARIABLES"):
                variables = re.findall(r'"([^"]+)"', text)
                continue
            if text.upper().startswith(("TITLE", "ZONE", "AUXDATA", "#")):
                continue
            if not variables:
                continue
            tokens = text.replace(",", " ").split()
            if len(tokens) != len(variables):
                raise CampaignError(
                    "%s: data row has %d columns; VARIABLES defines %d"
                    % (path, len(tokens), len(variables))
                )
            try:
                values = [float(token) for token in tokens]
            except ValueError as exc:
                raise CampaignError("%s contains a nonnumeric data row" % path) from exc
            if not all(math.isfinite(value) for value in values):
                raise CampaignError("%s contains a nonfinite value" % path)
            rows.append(values)
    if not variables or not rows:
        raise CampaignError("%s is not a nonempty Tecplot POINT table" % path)
    return variables, rows


def termination_evidence(run_dir: Path) -> Dict[str, object]:
    """Aggregate explicit trajectory accounting from every termination table."""

    paths = sorted(run_dir.glob("*termination_summary*.dat"))
    maximum_unresolved = None  # type: Optional[float]
    total_sampled = 0
    total_resolved = 0
    termination_counts: Dict[str, int] = {}
    for path in paths:
        variables, rows = parse_tecplot_table(path)
        index = {name: position for position, name in enumerate(variables)}
        if "unresolved_fraction" not in index or "N_sampled" not in index:
            raise CampaignError("termination table lacks required columns: %s" % path)
        for row in rows:
            value = row[index["unresolved_fraction"]]
            maximum_unresolved = value if maximum_unresolved is None else max(maximum_unresolved, value)
            total_sampled += int(round(row[index["N_sampled"]]))
            if "N_resolved" in index:
                total_resolved += int(round(row[index["N_resolved"]]))
            for name, position in index.items():
                if name.startswith("N_") and name not in ("N_sampled", "N_resolved", "N_retried"):
                    termination_counts[name[2:]] = termination_counts.get(name[2:], 0) + int(
                        round(row[position])
                    )
    return {
        "files": [str(path.relative_to(run_dir)) for path in paths],
        "table_count": len(paths),
        "max_unresolved_fraction": maximum_unresolved,
        "total_sampled": total_sampled,
        "total_resolved": total_resolved,
        "termination_counts": termination_counts,
    }


def _nested_metric(document: object, dotted_name: str) -> float:
    current = document
    for part in dotted_name.split("."):
        if not isinstance(current, Mapping) or part not in current:
            raise CampaignError("JSON metric is missing: %s" % dotted_name)
        current = current[part]
    try:
        value = float(current)
    except (TypeError, ValueError) as exc:
        raise CampaignError("JSON metric is not numeric: %s" % dotted_name) from exc
    if not math.isfinite(value):
        raise CampaignError("JSON metric is not finite: %s" % dotted_name)
    return value


def _compare(value: float, operator: str, threshold: float) -> bool:
    if operator == "<=":
        return value <= threshold
    if operator == "<":
        return value < threshold
    if operator == ">=":
        return value >= threshold
    if operator == ">":
        return value > threshold
    if operator == "==":
        return value == threshold
    raise CampaignError("unsupported gate operator: %s" % operator)


def evaluate_numerical_gates(
    run_dir: Path,
    termination: Mapping[str, object],
    gates: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    """Evaluate declared gates; missing required evidence is always a failure."""

    results: List[Dict[str, object]] = []
    for gate in gates:
        identifier = str(gate["id"])
        operator = str(gate["operator"])
        threshold = float(gate["threshold"])
        required = bool(gate.get("required", True))
        try:
            if gate["source"] == "termination_summary":
                metric = str(gate.get("metric", "max_unresolved_fraction"))
                value_object = termination.get(metric)
                if value_object is None:
                    raise CampaignError("termination metric is missing: %s" % metric)
                value = float(value_object)
                source_path = ",".join(str(item) for item in termination.get("files", []))
            else:
                relative = Path(str(gate["path"]))
                source_file = run_dir / relative
                if not source_file.is_file():
                    raise CampaignError("gate evidence file is missing: %s" % relative)
                with source_file.open(encoding="utf-8") as stream:
                    document = json.load(stream)
                value = _nested_metric(document, str(gate["metric"]))
                source_path = str(relative)
            passed = _compare(value, operator, threshold)
            results.append(
                {
                    "id": identifier,
                    "status": "PASS" if passed else "FAIL",
                    "value": value,
                    "operator": operator,
                    "threshold": threshold,
                    "required": required,
                    "source": source_path,
                }
            )
        except (CampaignError, OSError, json.JSONDecodeError, ValueError) as exc:
            results.append(
                {
                    "id": identifier,
                    "status": "FAIL" if required else "NOT_EVALUATED",
                    "value": None,
                    "operator": operator,
                    "threshold": threshold,
                    "required": required,
                    "error": str(exc),
                }
            )
    return results


def artifact_inventory(run_dir: Path, patterns: Sequence[str]) -> Tuple[List[Dict[str, object]], List[str]]:
    """Hash every expected artifact and report unmatched patterns separately."""

    found: Dict[str, Dict[str, object]] = {}
    missing: List[str] = []
    for pattern in patterns:
        matches = sorted(path for path in run_dir.glob(pattern) if path.is_file())
        if not matches:
            missing.append(pattern)
        for path in matches:
            relative = str(path.relative_to(run_dir))
            found[relative] = {
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    return [found[key] for key in sorted(found)], missing


def verified_restart(status_path: Path, fingerprint: str) -> bool:
    """Return true only for a PASS record whose artifacts still match hashes."""

    if not status_path.is_file():
        return False
    try:
        with status_path.open(encoding="utf-8") as stream:
            status = json.load(stream)
    except (OSError, json.JSONDecodeError):
        return False
    if status.get("status") != "PASS" or status.get("fingerprint") != fingerprint:
        return False
    run_dir = status_path.parent
    artifacts = status.get("artifacts", [])
    if not isinstance(artifacts, list) or not artifacts:
        return False
    for item in artifacts:
        if not isinstance(item, Mapping):
            return False
        path = run_dir / str(item.get("path", ""))
        if not path.is_file() or sha256_file(path) != item.get("sha256"):
            return False
    return True


def read_prediction_csv(path: Path) -> List[Dict[str, object]]:
    """Read the canonical per-snapshot prediction table used by Step 8.

    Required columns intentionally separate instrument identity from the model
    value.  This prevents a comparison script from guessing that, for example,
    two identically named P5 channels on different spacecraft are interchangeable.
    """

    required = {
        "utc", "instrument_id", "platform", "channel", "direction",
        "quantity", "units", "value",
    }
    rows: List[Dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise CampaignError("%s lacks prediction columns: %s" % (path, ", ".join(sorted(missing))))
        for line_number, raw in enumerate(reader, start=2):
            try:
                utc = format_utc(parse_utc(raw["utc"]))
                value = float(raw["value"])
                lower = float(raw["lower"]) if raw.get("lower", "").strip() else None
                upper = float(raw["upper"]) if raw.get("upper", "").strip() else None
            except (CampaignError, TypeError, ValueError) as exc:
                raise CampaignError("%s:%d invalid prediction row" % (path, line_number)) from exc
            if not math.isfinite(value) or (lower is not None and not math.isfinite(lower)) or (
                upper is not None and not math.isfinite(upper)
            ):
                raise CampaignError("%s:%d nonfinite prediction" % (path, line_number))
            if lower is not None and upper is not None and not lower <= value <= upper:
                raise CampaignError(
                    "%s:%d prediction is outside its declared bounds"
                    % (path, line_number)
                )
            rows.append(
                {
                    "utc": utc,
                    "instrument_id": raw["instrument_id"].strip(),
                    "platform": raw["platform"].strip(),
                    "channel": raw["channel"].strip(),
                    "direction": raw["direction"].strip().upper(),
                    "quantity": raw["quantity"].strip(),
                    "units": raw["units"].strip(),
                    "value": value,
                    "lower": lower,
                    "upper": upper,
                }
            )
            if not all(
                str(rows[-1][name]).strip()
                for name in (
                    "instrument_id", "platform", "channel", "direction",
                    "quantity", "units",
                )
            ):
                raise CampaignError(
                    "%s:%d prediction identity/unit field is empty"
                    % (path, line_number)
                )
    if not rows:
        raise CampaignError("prediction artifact is empty: %s" % path)
    return rows

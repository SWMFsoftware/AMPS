#!/usr/bin/env python3
"""Strict, versioned contract for SWCME observational campaigns.

This module deliberately validates evidence *before* any analysis command is
started.  Its job is governance rather than physics: a syntactically successful
empty campaign must never be reported as scientific validation.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List

SCHEMA_VERSION = 2
LAYER_IDS = ("V1", "V2", "V3", "V4", "V5")
OUTCOMES = {"PASS", "FAIL", "INCOMPLETE", "NOT_APPLICABLE", "ERROR"}
PURPOSES = {"RELEASE_VALIDATION", "SCHEMA_SELF_TEST"}
SOURCE_KINDS = {"OBSERVATION", "SYNTHETIC_REGRESSION"}
TOP_LEVEL_KEYS = {
    "schema_version", "event_id", "campaign_purpose", "description",
    "provenance", "events", "observers", "validation_layers", "sweeps", "convergence",
    "comparisons", "plots",
}


def _nonempty_list(value: Any) -> bool:
    return isinstance(value, list) and len(value) > 0


def validate_campaign_config(config: Any) -> Dict[str, Any]:
    """Return a deterministic schema/scientific-completeness classification.

    ERROR means malformed or unsupported input. INCOMPLETE means the structure
    is understood but required scientific evidence is absent. FAIL is reserved
    for a complete comparison that violates an acceptance threshold.
    """
    errors: List[str] = []
    missing: List[str] = []
    if not isinstance(config, dict):
        return {"status": "ERROR", "errors": ["campaign root must be an object"],
                "missing_requirements": []}

    unknown = sorted(set(config) - TOP_LEVEL_KEYS)
    if unknown:
        errors.append("unknown top-level keys: " + ", ".join(unknown))
    if config.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"schema_version must equal {SCHEMA_VERSION}")
    for field in ("event_id", "description"):
        if not isinstance(config.get(field), str) or not config[field].strip():
            errors.append(f"{field} must be a nonempty string")
    purpose = config.get("campaign_purpose")
    if purpose not in PURPOSES:
        errors.append("campaign_purpose must be RELEASE_VALIDATION or SCHEMA_SELF_TEST")

    # Event and observer inventories are first-class schema objects. Requiring
    # them prevents a metrics table from being detached from its time window or
    # spacecraft/ephemeris identity.
    for field,required_keys in (("events", ("id", "start_utc", "end_utc")),
                                ("observers", ("id", "mission", "position_frame"))):
        records=config.get(field)
        if not _nonempty_list(records):
            missing.append("at least one " + field[:-1] + " record")
            continue
        for index,record in enumerate(records):
            where=f"{field}[{index}]"
            if not isinstance(record,dict):
                errors.append(where + " must be an object")
                continue
            allowed=set(required_keys) | {"description", "provenance_ids"}
            extra=sorted(set(record)-allowed)
            if extra:
                errors.append(where + " has unknown keys: " + ", ".join(extra))
            for key in required_keys:
                if not isinstance(record.get(key),str) or not record[key].strip():
                    missing.append(where + "." + key)

    provenance = config.get("provenance")
    provenance_ids = set()
    observation_count = 0
    if not _nonempty_list(provenance):
        missing.append("at least one provenance record")
        provenance = []
    for index, record in enumerate(provenance):
        where = f"provenance[{index}]"
        if not isinstance(record, dict):
            errors.append(where + " must be an object")
            continue
        allowed = {"id", "source_kind", "source", "citation", "sha256",
                   "retrieved_utc", "processing"}
        extra = sorted(set(record) - allowed)
        if extra:
            errors.append(where + " has unknown keys: " + ", ".join(extra))
        rid = record.get("id")
        if not isinstance(rid, str) or not rid:
            errors.append(where + ".id must be nonempty")
        elif rid in provenance_ids:
            errors.append(where + ".id is duplicated")
        else:
            provenance_ids.add(rid)
        if record.get("source_kind") not in SOURCE_KINDS:
            errors.append(where + ".source_kind is invalid")
        elif record["source_kind"] == "OBSERVATION":
            observation_count += 1
        for field in ("source", "citation", "retrieved_utc", "processing"):
            if not isinstance(record.get(field), str) or not record[field].strip():
                missing.append(where + "." + field)
        digest = record.get("sha256")
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            errors.append(where + ".sha256 must be 64 lowercase hexadecimal digits")

    layers = config.get("validation_layers")
    if not isinstance(layers, list):
        errors.append("validation_layers must be an array")
        layers = []
    seen = set()
    layer_statuses: List[str] = []
    for index, layer in enumerate(layers):
        where = f"validation_layers[{index}]"
        if not isinstance(layer, dict):
            errors.append(where + " must be an object")
            continue
        allowed = {"id", "status", "description", "evidence", "metrics",
                   "uncertainties", "convergence", "provenance_ids", "rationale"}
        extra = sorted(set(layer) - allowed)
        if extra:
            errors.append(where + " has unknown keys: " + ", ".join(extra))
        layer_id = layer.get("id")
        if layer_id not in LAYER_IDS:
            errors.append(where + ".id must be one of " + ", ".join(LAYER_IDS))
            continue
        if layer_id in seen:
            errors.append(where + ".id is duplicated")
        seen.add(layer_id)
        status = layer.get("status")
        if status not in OUTCOMES:
            errors.append(where + ".status is invalid")
            continue
        layer_statuses.append(status)
        if not isinstance(layer.get("description"), str) or not layer["description"].strip():
            missing.append(where + ".description")
        if status in {"PASS", "FAIL"}:
            # A conclusion is auditable only when raw evidence, metrics,
            # uncertainties, convergence information, and data provenance are
            # all present. Arrays may contain one or many records, but cannot
            # be placeholders.
            for field in ("evidence", "metrics", "uncertainties", "convergence",
                          "provenance_ids"):
                if not _nonempty_list(layer.get(field)):
                    missing.append(where + "." + field)
            for rid in layer.get("provenance_ids", []):
                if rid not in provenance_ids:
                    errors.append(where + f" references unknown provenance id {rid!r}")
        elif status == "NOT_APPLICABLE":
            if not isinstance(layer.get("rationale"), str) or not layer["rationale"].strip():
                missing.append(where + ".rationale")
        elif status == "INCOMPLETE":
            if not isinstance(layer.get("rationale"), str) or not layer["rationale"].strip():
                missing.append(where + ".rationale")

    for layer_id in LAYER_IDS:
        if layer_id not in seen:
            missing.append("validation layer " + layer_id)
    for field in ("sweeps", "convergence", "comparisons", "plots"):
        if field in config and not isinstance(config[field], list):
            errors.append(field + " must be an array")

    if purpose == "RELEASE_VALIDATION":
        if observation_count == 0:
            missing.append("observational provenance for release validation")
        if "NOT_APPLICABLE" in layer_statuses:
            missing.append("all V1-V5 layers are mandatory for release validation")
        if not _nonempty_list(config.get("comparisons")):
            missing.append("at least one model/observation comparison")

    if errors:
        status = "ERROR"
    elif missing or "INCOMPLETE" in layer_statuses:
        status = "INCOMPLETE"
    elif "ERROR" in layer_statuses:
        status = "ERROR"
    elif "FAIL" in layer_statuses:
        status = "FAIL"
    elif len(seen) == len(LAYER_IDS) and all(s in {"PASS", "NOT_APPLICABLE"}
                                             for s in layer_statuses):
        status = "PASS"
    else:
        status = "INCOMPLETE"
    return {"status": status, "errors": sorted(errors),
            "missing_requirements": sorted(set(missing)),
            "observation_provenance_count": observation_count,
            "layers_present": sorted(seen)}

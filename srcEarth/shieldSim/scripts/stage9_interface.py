#!/usr/bin/env python3
"""Stage-9 request compiler, preflight reporter, and resumable case runner.

This module is intentionally independent of Geant4.  It translates a narrow,
versioned UI request into the already validated ``shieldSim.case.v2`` contract;
the C++ parser remains the final scientific authority.  The split gives the UI
fast, deterministic feedback without creating a second transport model.

Stage 9 is currently *verification only*.  The capability registry records why
mission/anomaly runs and individual quantities are not production-qualified.
Changing a label in the browser therefore cannot bypass a scientific gate.
"""

from __future__ import annotations

import argparse
import copy
import datetime as _datetime
import hashlib
import json
import math
import os
import pathlib
import re
import signal
import subprocess
import sys
import time
from statistics import NormalDist
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Tuple


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = ROOT / "data" / "stage9_capabilities.v1.json"
DEFAULT_MATERIALS = ROOT / "data" / "material_library.v1.json"

REQUEST_SCHEMA = "shieldSim.ui_request.v1"
MANIFEST_SCHEMA = "shieldSim.ui_run_manifest.v1"
JOB_SCHEMA = "shieldSim.ui_job_state.v1"
_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,95}$")


class Stage9Error(ValueError):
    """A user-facing request or bundle error with no hidden fallback."""


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def content_hash(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _read_json(path: pathlib.Path) -> Any:
    try:
        with path.open("r", encoding="utf-8") as stream:
            return json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise Stage9Error(f"cannot read JSON '{path}': {exc}") from exc


def _write_json(path: pathlib.Path, value: Any) -> None:
    """Write canonical, human-readable JSON using an atomic same-dir rename."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".new")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, ensure_ascii=False)
        stream.write("\n")
    temporary.replace(path)


def _object(value: Any, context: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise Stage9Error(f"{context} must be an object")
    return value


def _array(value: Any, context: str) -> List[Any]:
    if not isinstance(value, list):
        raise Stage9Error(f"{context} must be an array")
    return value


def _string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise Stage9Error(f"{context} must be a non-empty string")
    return value


def _number(value: Any, context: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise Stage9Error(f"{context} must be a number")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "finite and > 0" if positive else "finite"
        raise Stage9Error(f"{context} must be {qualifier}")
    return result


def _only(obj: Mapping[str, Any], allowed: Iterable[str], context: str) -> None:
    extra = sorted(set(obj) - set(allowed))
    if extra:
        raise Stage9Error(f"{context} has unknown field(s): {', '.join(extra)}")


def _identifier(value: Any, context: str) -> str:
    result = _string(value, context)
    if not _ID.fullmatch(result):
        raise Stage9Error(
            f"{context} must match {_ID.pattern!r}; spaces and path separators are forbidden")
    return result


def _safe_output_id(value: str) -> str:
    """Match shieldSim.cc SafeId for absorber-specific output prefixes."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", value)


def _measure_mm(value: Any, context: str, *, thickness: bool = False) -> float:
    obj = _object(value, context)
    _only(obj, {"value", "unit", "type"} if thickness else {"value", "unit"}, context)
    magnitude = _number(obj.get("value"), context + ".value", positive=True)
    unit = _string(obj.get("unit"), context + ".unit")
    factor = {"mm": 1.0, "cm": 10.0}.get(unit)
    if factor is None:
        raise Stage9Error(f"{context}.unit must be 'mm' or 'cm'")
    return magnitude * factor


def _measure_areal(value: Any, context: str) -> float:
    obj = _object(value, context)
    _only(obj, {"value", "unit", "type"}, context)
    magnitude = _number(obj.get("value"), context + ".value", positive=True)
    unit = _string(obj.get("unit"), context + ".unit")
    if unit != "g/cm2":
        raise Stage9Error(f"{context}.unit must be 'g/cm2'")
    return magnitude


def _path_from(base: pathlib.Path, text: str) -> pathlib.Path:
    path = pathlib.Path(text).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _registry(path: pathlib.Path) -> Dict[str, Any]:
    result = _object(_read_json(path), "capability registry")
    if result.get("schema_version") != "shieldSim.stage9_capabilities.v1":
        raise Stage9Error("unsupported Stage-9 capability registry schema")
    return result


def _profiles(registry: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    profiles = {}
    for index, raw in enumerate(_array(registry.get("profiles"), "registry.profiles")):
        profile = _object(raw, f"registry.profiles[{index}]")
        profile_id = _identifier(profile.get("profile_id"),
                                 f"registry.profiles[{index}].profile_id")
        if profile_id in profiles:
            raise Stage9Error(f"duplicate profile_id '{profile_id}'")
        profiles[profile_id] = profile
    return profiles


def _material_records(library_path: pathlib.Path,
                      overrides: List[Any]) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, str]]:
    """Build enough of the material index for UI geometry calculations.

    The C++ MaterialLibrary still performs complete composition, provenance,
    role, and uncertainty validation.  Here we deliberately resolve only IDs,
    aliases, roles, status, and density needed by the layer editor.
    """
    library = _object(_read_json(library_path), "material library")
    records: Dict[str, Dict[str, Any]] = {}
    aliases: Dict[str, str] = {}
    for index, raw in enumerate(_array(library.get("materials"),
                                       "material library.materials")):
        record = copy.deepcopy(_object(raw, f"materials[{index}]"))
        material_id = _identifier(record.get("material_id"),
                                  f"materials[{index}].material_id")
        records[material_id] = record
        for alias in record.get("aliases", []):
            aliases[_string(alias, f"materials[{index}].aliases")] = material_id

    # Mirror the important MergeOverride rule: an extended record gets a new
    # identity and never inherits aliases.  Full validity is checked by C++.
    for index, raw in enumerate(overrides):
        override = copy.deepcopy(_object(raw, f"material_overrides[{index}]"))
        material_id = _identifier(override.get("material_id"),
                                  f"material_overrides[{index}].material_id")
        if material_id in records or material_id in aliases:
            raise Stage9Error(f"material override ID '{material_id}' collides with the library")
        if "extends" in override:
            parent_ref = _string(override["extends"], f"material_overrides[{index}].extends")
            parent_id = parent_ref.rsplit("@", 1)[0]
            if parent_id not in records:
                raise Stage9Error(f"material override extends unknown '{parent_id}'")
            merged = copy.deepcopy(records[parent_id])
            merged.pop("aliases", None)
            merged.update(override)
            record = merged
        else:
            record = override
        records[material_id] = record
        for alias in record.get("aliases", []):
            alias = _string(alias, f"material_overrides[{index}].aliases")
            if alias in aliases or alias in records:
                raise Stage9Error(f"material override alias '{alias}' collides")
            aliases[alias] = material_id
    return records, aliases


def _resolve_material(token: Any, role: str, records: Mapping[str, Dict[str, Any]],
                      aliases: Mapping[str, str], context: str) -> Dict[str, Any]:
    selected = _string(token, context)
    material_id = aliases.get(selected, selected)
    if material_id not in records:
        raise Stage9Error(f"{context} selects unknown material '{selected}'")
    record = records[material_id]
    status = record.get("status")
    if status == "template":
        reason = record.get("applicability_and_limitations", "properties unresolved")
        raise Stage9Error(f"{context}: material '{selected}' is a disabled template: {reason}")
    if role not in record.get("roles", []):
        raise Stage9Error(f"{context}: material '{material_id}' is not approved for role '{role}'")
    density = _object(record.get("density"), context + ".resolved_density")
    if density.get("unit") != "g/cm3":
        raise Stage9Error(f"{context}: resolved density unit must be g/cm3")
    _number(density.get("value"), context + ".resolved_density.value", positive=True)
    return record


def _validate_transform(value: Any, context: str) -> None:
    transform = _object(value, context)
    _only(transform, {"translation", "rotation_quaternion_local_to_world"}, context)
    translation = _object(transform.get("translation"), context + ".translation")
    _only(translation, {"x", "y", "z"}, context + ".translation")
    for axis in ("x", "y", "z"):
        item = _object(translation.get(axis), f"{context}.translation.{axis}")
        _only(item, {"value", "unit"}, f"{context}.translation.{axis}")
        _number(item.get("value"), f"{context}.translation.{axis}.value")
        if item.get("unit") not in ("mm", "cm"):
            raise Stage9Error(f"{context}.translation.{axis}.unit must be mm or cm")
    quaternion = _object(transform.get("rotation_quaternion_local_to_world"),
                         context + ".rotation_quaternion_local_to_world")
    _only(quaternion, {"w", "x", "y", "z"},
          context + ".rotation_quaternion_local_to_world")
    components = [_number(quaternion.get(axis), f"{context}.quaternion.{axis}")
                  for axis in ("w", "x", "y", "z")]
    if abs(sum(x * x for x in components) - 1.0) > 1.0e-8:
        raise Stage9Error(f"{context} rotation quaternion must have unit norm within 1e-8")


def _absolute_dependency_paths(value: Any, base: pathlib.Path) -> None:
    """Relocate template file references so generated cases are self-consistent."""
    if isinstance(value, dict):
        for key, child in value.items():
            if key in {"path", "metadata_path"} and isinstance(child, str):
                path = pathlib.Path(child).expanduser()
                if not path.is_absolute():
                    value[key] = str((base / path).resolve())
            else:
                _absolute_dependency_paths(child, base)
    elif isinstance(value, list):
        for child in value:
            _absolute_dependency_paths(child, base)


def _accuracy_preflight(case: Mapping[str, Any], quantities: List[str]) -> List[Dict[str, Any]]:
    accuracy = _object(case.get("accuracy"), "generated_case.accuracy")
    overrides = _object(accuracy.get("quantity_overrides", {}),
                        "generated_case.accuracy.quantity_overrides")
    reports = []
    for quantity in quantities:
        effective = copy.deepcopy(accuracy)
        effective.update(_object(overrides.get(quantity, {}),
                                 f"accuracy.quantity_overrides.{quantity}"))
        target = _number(effective.get("target"), f"accuracy target for {quantity}", positive=True)
        scope = effective.get("scope")
        budget = _object(effective.get("systematic_budget", {}),
                         f"accuracy systematic budget for {quantity}")
        report: Dict[str, Any] = {"quantity": quantity, "target_relative": target,
                                  "scope": scope}
        known = budget.get("known_relative_floor")
        if known is not None and _number(known, "known_relative_floor") > target:
            report.update(status="unachievable_known_systematic_floor",
                          known_relative_floor=float(known))
        elif scope == "combined" and budget.get("status") != "complete":
            report.update(status="not_assessable_incomplete_systematic_budget")
        elif scope == "statistical":
            report.update(status="statistically_assessable_at_runtime",
                          statement="target applies to Monte Carlo statistical uncertainty only")
        else:
            # Structured budgets carry a parser-computed covariance result only
            # after C++ validation.  The UI reports their components and a
            # conservative normal-coverage estimate without claiming a pass.
            if "relative_standard_uncertainty" in budget:
                u = _number(budget["relative_standard_uncertainty"],
                            "relative_standard_uncertainty")
            else:
                components = budget.get("components", [])
                correlations = budget.get("correlations", [])
                by_id = {c["id"]: _number(c["relative_standard_uncertainty"],
                                           "systematic component uncertainty")
                         for c in components}
                variance = sum(u * u for u in by_id.values())
                for correlation in correlations:
                    i, j = correlation["component_i"], correlation["component_j"]
                    variance += 2.0 * by_id[i] * by_id[j] * _number(
                        correlation["correlation_coefficient"], "correlation coefficient")
                u = math.sqrt(max(0.0, variance))
            probability = effective.get("coverage", {}).get("interval_probability", 0.95)
            probability = _number(probability, "interval_probability")
            z = NormalDist().inv_cdf(0.5 + probability / 2.0)
            floor = z * u
            report.update(
                status=("unachievable_declared_systematic_floor" if floor >= target
                        else "assessable_at_runtime_with_remaining_statistical_budget"),
                declared_systematic_standard_relative=u,
                normal_coverage_systematic_floor_relative=floor,
                coverage_probability=probability)
        reports.append(report)
    return reports


def compile_request(request: Mapping[str, Any], registry_path: pathlib.Path = DEFAULT_REGISTRY,
                    material_path: pathlib.Path = DEFAULT_MATERIALS) -> Dict[str, Any]:
    """Validate a UI request and return cases plus a fully traceable manifest."""
    request = copy.deepcopy(_object(request, "request"))
    _only(request, {"schema_version", "request_id", "purpose", "profile_id",
                    "shield_configurations", "absorbers", "quantities", "accuracy_target",
                    "material_overrides", "runtime_calibration"}, "request")
    if request.get("schema_version") != REQUEST_SCHEMA:
        raise Stage9Error(f"request.schema_version must be '{REQUEST_SCHEMA}'")
    request_id = _identifier(request.get("request_id"), "request.request_id")
    purpose = _string(request.get("purpose"), "request.purpose")
    profile_id = _identifier(request.get("profile_id"), "request.profile_id")

    registry_path = registry_path.resolve()
    registry = _registry(registry_path)
    profiles = _profiles(registry)
    if profile_id not in profiles:
        raise Stage9Error(f"profile '{profile_id}' is unknown or disabled")
    profile = profiles[profile_id]
    enabled_purposes = profile.get("enabled_purposes", [])
    if purpose not in enabled_purposes:
        reason = profile.get("disabled_purpose_reasons", {}).get(
            purpose, "purpose has not passed the profile's scientific entry gate")
        raise Stage9Error(f"profile '{profile_id}' disables purpose '{purpose}': {reason}")

    quantities = [_string(item, f"request.quantities[{index}]") for index, item in
                  enumerate(_array(request.get("quantities"), "request.quantities"))]
    if not quantities or len(set(quantities)) != len(quantities):
        raise Stage9Error("request.quantities must contain unique values and must not be empty")
    allowed = set(profile.get("allowed_quantities", []))
    unsupported = [quantity for quantity in quantities if quantity not in allowed]
    if unsupported:
        raise Stage9Error(
            f"profile '{profile_id}' cannot combine quantity/quantities {unsupported}; "
            f"allowed values are {sorted(allowed)}")

    target = _number(request.get("accuracy_target", 0.10),
                     "request.accuracy_target", positive=True)
    if target >= 1.0:
        raise Stage9Error("request.accuracy_target must be a relative fraction in (0,1)")
    overrides = _array(request.get("material_overrides", []), "request.material_overrides")
    records, aliases = _material_records(material_path.resolve(), overrides)

    shields = _array(request.get("shield_configurations"),
                     "request.shield_configurations")
    absorbers = _array(request.get("absorbers"), "request.absorbers")
    if not shields:
        raise Stage9Error("request.shield_configurations must not be empty")
    if not absorbers:
        raise Stage9Error("request.absorbers must not be empty")

    resolved_shields: List[Dict[str, Any]] = []
    shield_ids = set()
    for shield_index, raw in enumerate(shields):
        context = f"request.shield_configurations[{shield_index}]"
        shield = _object(raw, context)
        _only(shield, {"shield_id", "label", "size_x", "size_y", "transform", "layers"},
              context)
        shield_id = _identifier(shield.get("shield_id"), context + ".shield_id")
        if shield_id in shield_ids:
            raise Stage9Error(f"duplicate shield_id '{shield_id}'")
        shield_ids.add(shield_id)
        _string(shield.get("label"), context + ".label")
        _measure_mm(shield.get("size_x"), context + ".size_x")
        _measure_mm(shield.get("size_y"), context + ".size_y")
        _validate_transform(shield.get("transform"), context + ".transform")
        layers = _array(shield.get("layers"), context + ".layers")
        if not 1 <= len(layers) <= 5:
            raise Stage9Error(f"{context}.layers must contain one through five layers")
        layer_ids = set()
        resolved_layers = []
        for layer_index, raw_layer in enumerate(layers):
            layer_context = f"{context}.layers[{layer_index}]"
            layer = _object(raw_layer, layer_context)
            _only(layer, {"layer_id", "material", "thickness"}, layer_context)
            layer_id = _identifier(layer.get("layer_id"), layer_context + ".layer_id")
            if layer_id in layer_ids:
                raise Stage9Error(f"{layer_context} duplicates layer_id '{layer_id}'")
            layer_ids.add(layer_id)
            material = _resolve_material(layer.get("material"), "shield", records, aliases,
                                         layer_context + ".material")
            density = float(material["density"]["value"])
            thickness = _object(layer.get("thickness"), layer_context + ".thickness")
            kind = thickness.get("type")
            if kind == "geometric":
                geometric_mm = _measure_mm(thickness, layer_context + ".thickness",
                                           thickness=True)
                areal = density * geometric_mm / 10.0
            elif kind == "areal_density":
                areal = _measure_areal(thickness, layer_context + ".thickness")
                geometric_mm = 10.0 * areal / density
            else:
                raise Stage9Error(
                    f"{layer_context}.thickness.type must be geometric or areal_density")
            if geometric_mm > 100.0 + 1.0e-12:
                raise Stage9Error(
                    f"{layer_context} derived geometric thickness {geometric_mm:g} mm "
                    "exceeds the approved 100 mm per-layer limit")
            resolved_layers.append({
                "source_to_absorber_order": layer_index + 1,
                "layer_id": layer_id,
                "material_id": material["material_id"],
                "density_g_cm3": density,
                "density_uncertainty": material["density"].get(
                    "uncertainty", {"status": "not_provided"}),
                "mass_thickness_g_cm2": areal,
                "geometric_thickness_mm": geometric_mm,
                "input_thickness": copy.deepcopy(thickness),
            })
        resolved_shields.append({"shield_id": shield_id, "label": shield["label"],
                                 "layers": resolved_layers, "case_geometry": copy.deepcopy(shield)})

    resolved_absorbers: List[Dict[str, Any]] = []
    absorber_ids = set()
    for absorber_index, raw in enumerate(absorbers):
        context = f"request.absorbers[{absorber_index}]"
        absorber = _object(raw, context)
        _only(absorber, {"absorber_id", "label", "material", "geometry"}, context)
        absorber_id = _identifier(absorber.get("absorber_id"), context + ".absorber_id")
        if absorber_id in absorber_ids:
            raise Stage9Error(f"duplicate absorber_id '{absorber_id}'")
        absorber_ids.add(absorber_id)
        _string(absorber.get("label"), context + ".label")
        material = _resolve_material(absorber.get("material"), "absorber", records, aliases,
                                     context + ".material")
        geometry = _object(absorber.get("geometry"), context + ".geometry")
        _only(geometry, {"type", "size_x", "size_y", "thickness", "transform"},
              context + ".geometry")
        if geometry.get("type") != "slab":
            raise Stage9Error(
                f"{context}.geometry.type must be 'slab'; arbitrary shapes remain Stage 10")
        # Every dimension is mandatory.  No material-name or UI-display default
        # is allowed to become geometry silently.
        for dimension in ("size_x", "size_y", "thickness"):
            _measure_mm(geometry.get(dimension), f"{context}.geometry.{dimension}")
        _validate_transform(geometry.get("transform"), context + ".geometry.transform")
        resolved_absorbers.append({
            "absorber_id": absorber_id, "label": absorber["label"],
            "provider": "slab", "material_id": material["material_id"],
            "density_g_cm3": float(material["density"]["value"]),
            "density_uncertainty": material["density"].get(
                "uncertainty", {"status": "not_provided"}),
            "case_geometry": copy.deepcopy(absorber),
        })

    template_path = _path_from(registry_path.parent,
                               _string(profile.get("case_template"),
                                       f"profile {profile_id}.case_template"))
    template = _object(_read_json(template_path), f"profile template {template_path}")
    if template.get("schema_version") != "shieldSim.case.v2":
        raise Stage9Error("Stage-9 profiles require a shieldSim.case.v2 template")
    template_file_sha256 = hashlib.sha256(template_path.read_bytes()).hexdigest()
    _absolute_dependency_paths(template, template_path.parent)
    # A caller-selected material library must be the one used by both the UI
    # preview and the C++ case, never merely the one recorded in the manifest.
    template["material_library"]["path"] = str(material_path.resolve())

    generated: List[Dict[str, Any]] = []
    mappings: List[Dict[str, Any]] = []
    for shield_index, shield in enumerate(resolved_shields):
        case = copy.deepcopy(template)
        case_id = f"{request_id}__{shield['shield_id']}"
        if len(case_id) > 160:
            raise Stage9Error("generated case_id exceeds 160 characters")
        case["case_id"] = case_id
        case["purpose"] = purpose
        case["material_overrides"] = copy.deepcopy(overrides)
        case["geometry"]["shield_stack"] = copy.deepcopy(shield["case_geometry"])
        case["geometry"]["shield_stack"].pop("shield_id", None)
        case["geometry"]["shield_stack"].pop("label", None)
        case_absorbers = []
        for absorber in resolved_absorbers:
            item = copy.deepcopy(absorber["case_geometry"])
            item.pop("label", None)
            case_absorbers.append(item)
            mappings.append({
                "matrix_index": len(mappings), "shield_id": shield["shield_id"],
                "absorber_id": absorber["absorber_id"], "case_id": case_id,
                "case_file": f"cases/{shield_index:03d}_{shield['shield_id']}.json",
                "case_json_paths": {
                    "shield": "$.geometry.shield_stack",
                    "absorber": ("$.geometry.absorbers[?(@.absorber_id=='" +
                                 absorber["absorber_id"] + "')]"),
                    "quantities": "$.quantities",
                    "accuracy": "$.accuracy",
                },
                "expected_output_prefix": case_id + "__" +
                                          _safe_output_id(absorber["absorber_id"]),
            })
        case["geometry"]["absorbers"] = case_absorbers
        case["quantities"] = list(quantities)
        case["accuracy"]["target"] = target
        for override in case["accuracy"].get("quantity_overrides", {}).values():
            override["target"] = target
        base_seed = int(case["transport"]["random_seed"])
        case["transport"]["random_seed"] = base_seed + shield_index
        case["output"]["prefix"] = case_id
        generated.append({"case_id": case_id,
                          "relative_path": f"cases/{shield_index:03d}_{shield['shield_id']}.json",
                          "case": case, "case_hash": content_hash(case)})

    sampling = template["accuracy"]["sampling"]
    case_count = len(resolved_shields) * len(resolved_absorbers)
    min_histories = case_count * int(sampling["min_histories"])
    max_histories = case_count * int(sampling["max_histories"])
    calibration = _object(request.get("runtime_calibration", {"status": "not_provided"}),
                          "request.runtime_calibration")
    _only(calibration, {"status", "seconds_per_history", "provenance"},
          "request.runtime_calibration")
    cost: Dict[str, Any] = {
        "transport_case_count": case_count,
        "generated_case_file_count": len(resolved_shields),
        "requested_quantity_series_count": case_count * len(quantities),
        "minimum_primary_histories": min_histories,
        "maximum_primary_histories": max_histories,
        "history_definition": "sum over independent shield x absorber transports",
    }
    if calibration.get("status") == "measured":
        seconds = _number(calibration.get("seconds_per_history"),
                          "runtime_calibration.seconds_per_history", positive=True)
        provenance = _string(calibration.get("provenance"),
                             "runtime_calibration.provenance")
        cost.update(wall_time_status="calibrated_estimate",
                    seconds_per_history=seconds, calibration_provenance=provenance,
                    minimum_wall_seconds=min_histories * seconds,
                    maximum_wall_seconds=max_histories * seconds)
    elif calibration.get("status") == "not_provided":
        cost.update(wall_time_status="not_assessable_no_measured_throughput")
    else:
        raise Stage9Error("runtime_calibration.status must be measured or not_provided")

    badges = [copy.deepcopy(badge) for badge in profile.get("validation_badges", [])
              if badge.get("quantity") in quantities]
    regression_contract = {
        "request_hash": content_hash(request), "profile_id": profile_id,
        "template_file_sha256": template_file_sha256,
        "case_matrix": {"shield_count": len(resolved_shields),
                        "absorber_count": len(resolved_absorbers),
                        "transport_case_count": case_count},
        "ordered_shield_layers": [{"shield_id": item["shield_id"],
                                   "layers": item["layers"]} for item in resolved_shields],
        "absorbers": [{key: value for key, value in item.items() if key != "case_geometry"}
                      for item in resolved_absorbers],
        "quantities": quantities, "accuracy_target_relative": target,
    }
    manifest: Dict[str, Any] = {
        "schema_version": MANIFEST_SCHEMA,
        "stage9_qualification_status": registry.get("stage9_qualification_status"),
        "request_id": request_id, "purpose": purpose, "profile_id": profile_id,
        "request_hash": content_hash(request),
        # Embedding the canonical request makes recovery possible from the
        # manifest alone; request.canonical.json is a convenience duplicate.
        "ui_request": copy.deepcopy(request),
        "capability_registry": {"path": str(registry_path),
                                "sha256": content_hash(registry)},
        "material_library": {"path": str(material_path.resolve()),
                             "sha256": hashlib.sha256(material_path.read_bytes()).hexdigest()},
        "case_matrix": {"shield_count": len(resolved_shields),
                        "absorber_count": len(resolved_absorbers),
                        "transport_case_count": case_count,
                        "formula": "shield_count * absorber_count",
                        "mappings": mappings},
        "ordered_shield_layers": [{"shield_id": item["shield_id"],
                                   "layers": item["layers"]} for item in resolved_shields],
        "resolved_absorbers": [{key: value for key, value in item.items()
                                if key != "case_geometry"} for item in resolved_absorbers],
        "quantities": quantities, "accuracy_target_relative": target,
        "uncertainty_preflight": _accuracy_preflight(generated[0]["case"], quantities),
        "cost_preflight": cost,
        "validation_badges": badges,
        "disabled_scientific_states": copy.deepcopy(registry.get("disabled_states", [])),
        "generated_cases": [{key: value for key, value in item.items() if key != "case"}
                            for item in generated],
        "reproducibility": {
            "canonical_json": "RFC-8259 JSON, UTF-8, keys sorted for SHA-256",
            "seed_policy": "profile template seed plus zero-based shield index",
            "request_round_trip_required": True,
            "portable_regression_contract_sha256": content_hash(regression_contract),
            "portable_regression_contract": regression_contract,
        },
        "runner_boundary": {
            "checkpoint_granularity": "generated shield case file",
            "resume_behavior": "completed case files skipped; interrupted case restarts",
            "within_case_geant4_checkpoint": False,
        },
    }
    manifest["manifest_content_hash"] = content_hash(manifest)
    return {"request": request, "manifest": manifest, "generated_cases": generated}


def write_bundle(compiled: Mapping[str, Any], output: pathlib.Path) -> pathlib.Path:
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise Stage9Error(f"bundle output '{output}' exists and is not empty")
    output.mkdir(parents=True, exist_ok=True)
    _write_json(output / "request.canonical.json", compiled["request"])
    for item in compiled["generated_cases"]:
        _write_json(output / item["relative_path"], item["case"])
    manifest = copy.deepcopy(compiled["manifest"])
    manifest["bundle_root"] = str(output)
    # bundle_root is operational metadata and deliberately not folded into the
    # location-independent manifest_content_hash.
    _write_json(output / "manifest.json", manifest)
    return output


def _utc_now() -> str:
    return _datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z")


def request_cancel(bundle: pathlib.Path) -> pathlib.Path:
    marker = bundle.resolve() / "CANCEL_REQUESTED"
    marker.write_text(_utc_now() + "\n", encoding="utf-8")
    return marker


def _verify_bundle(bundle: pathlib.Path, manifest: Mapping[str, Any]) -> None:
    """Authenticate every interface-controlled input before launching Geant4."""
    hash_input = copy.deepcopy(dict(manifest))
    expected_manifest_hash = hash_input.pop("manifest_content_hash", None)
    hash_input.pop("bundle_root", None)  # location metadata is intentionally unhashed
    if expected_manifest_hash != content_hash(hash_input):
        raise Stage9Error("bundle manifest content hash mismatch")
    request = _object(manifest.get("ui_request"), "manifest.ui_request")
    if content_hash(request) != manifest.get("request_hash"):
        raise Stage9Error("embedded UI request hash mismatch")
    request_file = bundle / "request.canonical.json"
    if not request_file.is_file() or content_hash(_read_json(request_file)) != manifest.get("request_hash"):
        raise Stage9Error("request.canonical.json is missing or does not match the manifest")
    material = _object(manifest.get("material_library"), "manifest.material_library")
    material_path = pathlib.Path(_string(material.get("path"), "material_library.path"))
    if (not material_path.is_file() or
            hashlib.sha256(material_path.read_bytes()).hexdigest() != material.get("sha256")):
        raise Stage9Error("material library is missing or has changed since bundle compilation")
    for entry in _array(manifest.get("generated_cases"), "manifest.generated_cases"):
        case_entry = _object(entry, "manifest.generated_cases[]")
        relative = pathlib.Path(_string(case_entry.get("relative_path"), "generated case path"))
        if relative.is_absolute() or ".." in relative.parts:
            raise Stage9Error("generated case path escapes the bundle")
        case_path = bundle / relative
        if (not case_path.is_file() or
                content_hash(_read_json(case_path)) != case_entry.get("case_hash")):
            raise Stage9Error(f"generated case is missing or changed: {relative}")


def run_bundle(bundle: pathlib.Path, executable: pathlib.Path, *, resume: bool = False,
               poll_seconds: float = 0.1, max_new_cases: Optional[int] = None) -> Dict[str, Any]:
    """Run generated files with durable state and case-boundary resume.

    Cancellation is cooperative through ``CANCEL_REQUESTED``.  If a case is
    active, SIGTERM is sent and escalated to SIGKILL after five seconds.  A
    later resume never pretends that partial Geant4 histories are reusable.
    """
    bundle = bundle.resolve()
    manifest = _object(_read_json(bundle / "manifest.json"), "bundle manifest")
    if manifest.get("schema_version") != MANIFEST_SCHEMA:
        raise Stage9Error("bundle has unsupported manifest schema")
    _verify_bundle(bundle, manifest)
    executable = executable.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise Stage9Error(f"runner executable is not executable: {executable}")
    state_path = bundle / "job_state.json"
    old = _read_json(state_path) if state_path.exists() else None
    if old and not resume:
        raise Stage9Error("job_state.json exists; use --resume to continue without overwriting it")
    completed = set(old.get("completed_case_ids", [])) if old else set()
    attempt = int(old.get("attempt", 0)) + 1 if old else 1
    cancel_marker = bundle / "CANCEL_REQUESTED"
    if resume and cancel_marker.exists():
        cancel_marker.unlink()
    state: Dict[str, Any] = {
        "schema_version": JOB_SCHEMA, "request_id": manifest["request_id"],
        "manifest_content_hash": manifest["manifest_content_hash"],
        "status": "running", "attempt": attempt, "started_utc": _utc_now(),
        "updated_utc": _utc_now(), "executable": str(executable),
        "completed_case_ids": sorted(completed), "current_case_id": None,
        "case_attempts": copy.deepcopy(old.get("case_attempts", {})) if old else {},
        "checkpoint_granularity": "generated shield case file",
    }
    _write_json(state_path, state)
    newly_completed = 0
    active: Optional[subprocess.Popen] = None

    def stop_active() -> None:
        nonlocal active
        if active is None or active.poll() is not None:
            return
        active.terminate()
        deadline = time.monotonic() + 5.0
        while active.poll() is None and time.monotonic() < deadline:
            time.sleep(0.05)
        if active.poll() is None:
            active.kill()
        active.wait()

    previous_handlers: Dict[int, Any] = {}

    def signal_handler(_signum: int, _frame: Any) -> None:
        request_cancel(bundle)

    for signum in (signal.SIGINT, signal.SIGTERM):
        previous_handlers[signum] = signal.signal(signum, signal_handler)
    try:
        for case_entry in manifest["generated_cases"]:
            case_id = case_entry["case_id"]
            if case_id in completed:
                continue
            if cancel_marker.exists() or (max_new_cases is not None and
                                          newly_completed >= max_new_cases):
                state["status"] = "canceled"
                state["cancellation_boundary"] = "between_generated_case_files"
                break
            case_file = bundle / case_entry["relative_path"]
            run_dir = bundle / "runs" / case_id
            run_dir.mkdir(parents=True, exist_ok=True)
            case_mappings = [item for item in manifest["case_matrix"]["mappings"]
                             if item["case_id"] == case_id]
            _write_json(run_dir / "stage9_provenance.json", {
                "schema_version": "shieldSim.ui_output_provenance.v1",
                "request_id": manifest["request_id"], "case_id": case_id,
                "manifest_content_hash": manifest["manifest_content_hash"],
                "case_hash": case_entry["case_hash"],
                "validation_badges": manifest["validation_badges"],
                "stage9_qualification_status": manifest["stage9_qualification_status"],
                "matrix_mappings": case_mappings,
            })
            log_path = run_dir / "run.log"
            state["current_case_id"] = case_id
            attempts = state["case_attempts"].setdefault(case_id, [])
            attempt_record = {"attempt": len(attempts) + 1, "started_utc": _utc_now(),
                              "case_sha256": case_entry["case_hash"],
                              "log": str(log_path.relative_to(bundle))}
            attempts.append(attempt_record)
            state["updated_utc"] = _utc_now()
            _write_json(state_path, state)
            with log_path.open("w", encoding="utf-8") as log:
                active = subprocess.Popen([str(executable), "--case=" + str(case_file)],
                                          cwd=run_dir, stdout=log,
                                          stderr=subprocess.STDOUT, text=True)
                while active.poll() is None:
                    if cancel_marker.exists():
                        stop_active()
                        break
                    time.sleep(poll_seconds)
                return_code = active.returncode
            attempt_record["finished_utc"] = _utc_now()
            attempt_record["return_code"] = return_code
            if cancel_marker.exists():
                attempt_record["status"] = "interrupted_restart_required"
                state["status"] = "canceled"
                state["cancellation_boundary"] = "active_case_terminated_partial_not_reused"
                break
            if return_code != 0:
                attempt_record["status"] = "failed"
                state["status"] = "failed"
                state["failure"] = {"case_id": case_id, "return_code": return_code,
                                    "log": str(log_path.relative_to(bundle))}
                break
            output_error = None
            for mapping in case_mappings:
                prefix = mapping["expected_output_prefix"]
                resolved_path = run_dir / (prefix + "_resolved_case.json")
                result_path = run_dir / (prefix + "_result.json")
                try:
                    resolved = _object(_read_json(resolved_path), str(resolved_path))
                    result = _object(_read_json(result_path), str(result_path))
                    if resolved.get("case_id") != case_id:
                        raise Stage9Error("resolved-case case_id does not match the manifest")
                    if resolved.get("active_absorber_id") != mapping["absorber_id"]:
                        raise Stage9Error("resolved-case absorber does not match the matrix cell")
                    if result.get("absorber_id") != mapping["absorber_id"]:
                        raise Stage9Error("result absorber does not match the matrix cell")
                    if result.get("case_hash") != resolved.get("case_hash"):
                        raise Stage9Error("result/resolved-case scientific hashes differ")
                except Stage9Error as exc:
                    output_error = str(exc)
                    break
            if output_error is not None:
                attempt_record["status"] = "failed_output_contract"
                state["status"] = "failed"
                state["failure"] = {"case_id": case_id,
                                    "reason": "missing_or_invalid_expected_output",
                                    "detail": output_error,
                                    "log": str(log_path.relative_to(bundle))}
                break
            attempt_record["status"] = "complete"
            completed.add(case_id)
            newly_completed += 1
            state["completed_case_ids"] = sorted(completed)
            state["current_case_id"] = None
            state["updated_utc"] = _utc_now()
            _write_json(state_path, state)
        else:
            state["status"] = "complete"
    finally:
        stop_active()
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
        state["current_case_id"] = None
        state["updated_utc"] = _utc_now()
        state["finished_utc"] = _utc_now()
        _write_json(state_path, state)
    return state


def _print(value: Any) -> None:
    json.dump(value, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=pathlib.Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--materials", type=pathlib.Path, default=DEFAULT_MATERIALS)
    commands = parser.add_subparsers(dest="command", required=True)
    preview = commands.add_parser("preview", help="validate and print preflight manifest")
    preview.add_argument("--request", required=True, type=pathlib.Path)
    compile_cmd = commands.add_parser("compile", help="write reproducible run bundle")
    compile_cmd.add_argument("--request", required=True, type=pathlib.Path)
    compile_cmd.add_argument("--output", required=True, type=pathlib.Path)
    run = commands.add_parser("run", help="execute or resume a compiled bundle")
    run.add_argument("--bundle", required=True, type=pathlib.Path)
    run.add_argument("--executable", required=True, type=pathlib.Path)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--max-new-cases", type=int, help=argparse.SUPPRESS)
    cancel = commands.add_parser("cancel", help="request cooperative cancellation")
    cancel.add_argument("--bundle", required=True, type=pathlib.Path)
    args = parser.parse_args(argv)
    try:
        if args.command in ("preview", "compile"):
            request = _read_json(args.request.resolve())
            compiled = compile_request(request, args.registry, args.materials)
            if args.command == "preview":
                _print(compiled["manifest"])
            else:
                path = write_bundle(compiled, args.output)
                _print({"status": "compiled", "bundle": str(path),
                        "manifest_content_hash": compiled["manifest"]["manifest_content_hash"]})
        elif args.command == "run":
            state = run_bundle(args.bundle, args.executable, resume=args.resume,
                               max_new_cases=args.max_new_cases)
            _print(state)
            return 0 if state["status"] in ("complete", "canceled") else 1
        else:
            _print({"status": "cancellation_requested",
                    "marker": str(request_cancel(args.bundle))})
    except Stage9Error as exc:
        print(f"Stage-9 error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Nonnegative variable-projection fitting for qualified density profiles."""
import copy
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, lsq_linear

from plume_model import (POINT_SOURCE_LAW, Plume, design_matrix, line_source,
                         spherical_source)


SOURCE_CASES = (
    "tenishev2014_hybrid_8",
    "porco2014_98",
    "portyankina2022_solar_41",
)


def _external_path(config_path, configured, override, environment,
                   validation_data_root=None):
    """Resolve an external data product without embedding it in the code tree."""
    name = override or configured or os.environ.get(environment)
    if not name:
        root = (validation_data_root
                if validation_data_root is not None
                else os.environ.get("ENCELADUS_VALIDATION_DATA_ROOT"))
        if root:
            name = (Path(root).expanduser()/"data"/"derived"
                    /"porco2014_jet_catalog.json")
    if not name:
        raise ValueError(
            "Source case requires an external jet catalog; set jet_catalog_file, "
            "pass --jet-catalog, set "+environment+", or supply the common "
            "validation root with --validation-data-root or "
            "ENCELADUS_VALIDATION_DATA_ROOT")
    path = Path(name).expanduser()
    if not path.is_absolute():
        path = config_path.parent/path
    return path.resolve()


def _materialize_catalog_case(config, catalog, source_case):
    """Replace the legacy source geometry with a checked catalog selection."""
    if (catalog.get("schema_version") != 1
            or catalog.get("frame") != "IAU_ENCELADUS"
            or catalog.get("jet_count") != 98):
        raise ValueError("Jet catalog must be schema 1 with 98 IAU_ENCELADUS jets")
    jets = catalog.get("jets")
    if not isinstance(jets, list) or len(jets) != 98:
        raise ValueError("Jet catalog must contain exactly 98 jet records")
    by_id = {item.get("porco_id"): item for item in jets}
    if set(by_id) != set(range(1, 99)):
        raise ValueError("Jet catalog porco_id values must be exactly 1..98")

    options = config.get("source_case_options", {}).get(source_case)
    if not isinstance(options, dict):
        raise ValueError("Missing source_case_options for "+source_case)
    expected_selection = {
        "porco2014_98": ("porco2014_98", 98),
        "portyankina2022_solar_41": ("portyankina2022_solar_41", 41),
    }[source_case]
    selection_name = options.get("catalog_selection")
    if selection_name != expected_selection[0]:
        raise ValueError(
            source_case+" must use catalog selection "+expected_selection[0])
    selected = catalog.get("selections", {}).get(selection_name)
    if (not isinstance(selected, list) or len(selected) != expected_selection[1]
            or len(set(selected)) != len(selected)
            or any(item not in by_id for item in selected)):
        raise ValueError("Invalid catalog selection "+str(selection_name))

    temperature = options.get("temperature_K")
    speed = options.get("bulk_speed_m_s")
    if (not isinstance(temperature, (int, float))
            or not isinstance(speed, (int, float))
            or not np.isfinite([temperature, speed]).all()
            or temperature <= 0 or speed < 0):
        raise ValueError("Catalog source case requires explicit finite T>0 and u>=0")
    direction_mode = options.get("direction_mode")
    if direction_mode not in ("published", "local_normal"):
        raise ValueError("direction_mode must be published or local_normal")
    if options.get("include_tiger_stripes", False):
        raise ValueError("The 98- and 41-jet cases must not also include stripes")

    materialized = []
    for porco_id in selected:
        item = by_id[porco_id]
        required = ("id", "latitude_deg", "west_longitude_deg",
                    "tilt_deg", "azimuth_deg", "sulcus")
        if any(item.get(key) is None for key in required):
            raise ValueError("Incomplete geometry for Porco jet "+str(porco_id))
        source = {
            "id": item["id"],
            "latitude_deg": item["latitude_deg"],
            "west_longitude_deg": item["west_longitude_deg"],
            "tilt_deg": (item["tilt_deg"] if direction_mode == "published" else 0.0),
            "azimuth_deg": (item["azimuth_deg"] if direction_mode == "published" else 0.0),
            "bulk_speed_m_s": speed,
            "temperature_K": temperature,
            "porco_id": porco_id,
            "sulcus": item["sulcus"],
            "geometry_quality": item.get("geometry_quality"),
            "published_tilt_deg": item["tilt_deg"],
            "published_azimuth_deg": item["azimuth_deg"],
        }
        if "portyankina_2010_solar_coefficient" in item:
            source["portyankina_2010_solar_coefficient"] = item[
                "portyankina_2010_solar_coefficient"]
        materialized.append(source)
    config["discrete_jets"] = materialized
    config["tiger_stripes"] = []
    config["resolved_source_case"] = {
        "name": source_case,
        "catalog_selection": selection_name,
        "jet_count": len(materialized),
        "direction_mode": direction_mode,
        "temperature_K": temperature,
        "bulk_speed_m_s": speed,
        "kinematics_basis": options.get("kinematics_basis"),
    }


def load_config(path, source_case=None, jet_catalog_path=None,
                validation_data_root=None):
    """Load JSON and resolve one of the three published source configurations."""
    config_path = Path(path).resolve()
    config = json.loads(config_path.read_text())
    selected_case = source_case or config.get(
        "source_case", "tenishev2014_hybrid_8")
    if selected_case not in SOURCE_CASES:
        raise ValueError("Unknown source_case: "+str(selected_case))
    config["source_case"] = selected_case
    if selected_case != "tenishev2014_hybrid_8":
        catalog_path = _external_path(
            config_path, config.get("jet_catalog_file"), jet_catalog_path,
            "ENCELADUS_JET_CATALOG", validation_data_root)
        catalog = json.loads(catalog_path.read_text())
        _materialize_catalog_case(config, catalog, selected_case)
        config["_jet_catalog_path"] = str(catalog_path)
    else:
        config["resolved_source_case"] = {
            "name": selected_case,
            "jet_count": len(config.get("discrete_jets", [])),
            "stripe_count": len(config.get("tiger_stripes", [])),
            "direction_mode": "published_legacy",
        }
    geometry_name = config.get("stripe_geometry_file")
    if geometry_name and config.get("tiger_stripes"):
        geometry_path = (config_path.parent/geometry_name).resolve()
        geometry = json.loads(geometry_path.read_text())
        if geometry.get("frame") != "IAU_ENCELADUS" or geometry.get("units") != "m":
            raise ValueError("Stripe geometry must be IAU_ENCELADUS in metres")
        table = geometry.get("stripes", {})
        for stripe in config.get("tiger_stripes", []):
            key = stripe.get("geometry_id")
            if key not in table:
                raise ValueError("Missing stripe geometry: "+str(key))
            stripe["points_m"] = table[key]
        config["_stripe_geometry_path"] = str(geometry_path)
    return config


def append_extra_jet(config, definition_path):
    """Append one provenance-tracked development jet to a resolved source case."""
    path = Path(definition_path).expanduser().resolve()
    document = json.loads(path.read_text(encoding="utf-8"))
    if (document.get("schema_version") != 1
            or document.get("frame") != "IAU_ENCELADUS"
            or document.get("base_source_case") != config.get("source_case")):
        raise ValueError(
            "Extra-jet definition must be schema 1 in IAU_ENCELADUS and match "
            "the resolved base source case")
    variant = document.get("variant_name")
    source = document.get("source")
    if not isinstance(variant, str) or not variant or not isinstance(source, dict):
        raise ValueError("Extra-jet definition requires variant_name and source")
    required = ("id", "latitude_deg", "west_longitude_deg", "tilt_deg",
                "azimuth_deg", "temperature_K", "bulk_speed_m_s")
    if any(source.get(key) is None for key in required):
        raise ValueError("Extra-jet source geometry or kinematics is incomplete")
    if any(item.get("id") == source["id"]
           for item in config.get("discrete_jets", [])):
        raise ValueError("Duplicate extra-jet source ID: "+str(source["id"]))
    numeric = [source[key] for key in required[1:]]
    if (not np.isfinite(numeric).all()
            or not -90 <= source["latitude_deg"] <= 90
            or not 0 <= source["west_longitude_deg"] < 360
            or not 0 <= source["tilt_deg"] <= 90
            or source["temperature_K"] <= 0
            or source["bulk_speed_m_s"] < 0):
        raise ValueError("Invalid extra-jet geometry or kinematics")
    candidate = copy.deepcopy(source)
    config.setdefault("discrete_jets", []).append(candidate)
    config["_extra_jet_path"] = str(path)
    config["source_variant"] = variant
    resolved = config.setdefault("resolved_source_case", {})
    resolved["base_name"] = config["source_case"]
    resolved["variant_name"] = variant
    resolved["extra_jet_id"] = candidate["id"]
    resolved["jet_count"] = len(config["discrete_jets"])
    return config


def sources_from_config(config):
    if config.get("physics") != "gravity_free_collisionless":
        raise ValueError("physics must explicitly be gravity_free_collisionless")
    if config.get("frame") != "IAU_ENCELADUS":
        raise ValueError("Model and observations must use IAU_ENCELADUS")
    radius = config.get("reference_radius_m")
    mass = config.get("molecular_mass_kg")
    if not np.isfinite([radius, mass]).all() or radius <= 0 or mass <= 0:
        raise ValueError("Positive reference radius and molecular mass required")
    source_law = config.get("source_law", POINT_SOURCE_LAW)

    result, labels = [], []
    for source in config.get("discrete_jets", []):
        label = "jet:"+source["id"]
        if label in labels:
            raise ValueError("Duplicate source ID: "+source["id"])
        position, normal, axis = spherical_source(
            radius, source["latitude_deg"], source["west_longitude_deg"],
            source["tilt_deg"], source["azimuth_deg"])
        result.append(Plume(tuple(position), tuple(normal), tuple(axis),
                            source["temperature_K"], source["bulk_speed_m_s"],
                            mass, source.get("source_law", source_law)))
        labels.append(label)

    for stripe in config.get("tiger_stripes", []):
        label = "stripe:"+stripe["id"]
        if label in labels:
            raise ValueError("Duplicate source ID: "+stripe["id"])
        if "points_m" not in stripe:
            raise ValueError("Unresolved stripe geometry for "+stripe["id"])
        result.append(line_source(
            stripe["points_m"], stripe["temperature_K"],
            stripe["bulk_speed_m_s"], mass,
            stripe.get("source_law", source_law)))
        labels.append(label)

    background = config.get("background", {})
    selection = background.get("model", "both")
    if selection in ("constant", "both"):
        labels.append("background:c0_m-3")
    if selection in ("power_law", "both"):
        labels.append("background:density_at_reference_radius_m-3")
    if not result and not labels:
        raise ValueError("No sources or background configured")
    return result, labels


def matrix_from_config(config, points):
    sources, labels = sources_from_config(config)
    background = config.get("background", {})
    matrix = design_matrix(
        sources, points, config["reference_radius_m"],
        background.get("model", "both"),
        background.get("power_law_exponent", -2.0))
    return matrix, labels


def solve_amplitudes(matrix, density, sigma):
    """Weighted nonnegative least squares with scaling and rank diagnostics."""
    a = np.asarray(matrix, dtype=float)
    y = np.asarray(density, dtype=float)
    sd = np.asarray(sigma, dtype=float)
    if a.ndim != 2 or y.shape != (len(a),) or sd.shape != y.shape or len(y) == 0:
        raise ValueError("Invalid design/data dimensions")
    if (not all(np.isfinite(v).all() for v in (a, y, sd))
            or np.any(y < 0) or np.any(sd <= 0)):
        raise ValueError("Finite nonnegative density and positive sigma required")
    weighted = a/sd[:, None]
    scale = np.linalg.norm(weighted, axis=0)
    observable = scale > 0
    safe_scale = np.where(observable, scale, 1.0)
    scaled = weighted/safe_scale
    fit = lsq_linear(scaled, y/sd, bounds=(0, np.inf),
                     tol=1e-12, max_iter=5000)
    if not fit.success:
        raise RuntimeError("Amplitude solve failed: "+fit.message)
    coefficients = fit.x/safe_scale
    coefficients[~observable] = 0.0
    predicted = a @ coefficients
    singular = np.linalg.svd(scaled, compute_uv=False)
    rank = int(np.linalg.matrix_rank(scaled))
    relative_singular = (singular/singular[0] if len(singular) and singular[0] > 0
                         else np.zeros_like(singular))
    # A formal floating-point rank alone is misleading for this inverse
    # problem: nearby jets and stripe nodes can produce nearly parallel
    # flyby signatures.  Report several explicit resolution thresholds rather
    # than silently selecting one as a physical prior.
    effective_rank = {
        "1e-2": int(np.count_nonzero(relative_singular >= 1.0e-2)),
        "1e-3": int(np.count_nonzero(relative_singular >= 1.0e-3)),
        "1e-4": int(np.count_nonzero(relative_singular >= 1.0e-4)),
    }
    condition = (float(singular[0]/singular[-1])
                 if len(singular) and singular[-1] > 0 else float("inf"))
    warnings = []
    if condition > 1.0e6:
        warnings.append(
            "The scaled design matrix is severely ill-conditioned; individual "
            "jet/stripe amplitudes are not uniquely identifiable from this "
            "profile. Total production and its partition are model-dependent.")
    if rank < a.shape[1]:
        warnings.append(
            "The design matrix is formally rank deficient; at least one linear "
            "combination of amplitudes is unobservable.")
    return {
        "coefficients": coefficients,
        "predicted": predicted,
        "residual": (predicted-y)/sd,
        "chi2": float(np.sum(((predicted-y)/sd)**2)),
        "rank": rank,
        "singular_values": singular,
        "relative_singular_values": relative_singular,
        "effective_rank_by_relative_singular_cutoff": effective_rank,
        "condition_number_scaled": condition,
        "unobserved_columns": np.where(~observable)[0],
        "rank_deficient": rank < a.shape[1],
        "fit_warnings": warnings,
    }


def fit_density(config, points, density, sigma):
    """Fit source/background amplitudes and optional bounded jet shapes."""
    initial = copy.deepcopy(config)
    sources_from_config(initial)
    definitions = initial.get("shape_parameters", [])
    allowed = {"temperature_K", "bulk_speed_m_s", "tilt_deg", "azimuth_deg"}
    jets = {source["id"]: index
            for index, source in enumerate(initial.get("discrete_jets", []))}
    keys = [(item["source_id"], item["field"]) for item in definitions]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate shape parameter")

    x0, lower, upper = [], [], []
    for item in definitions:
        if item["field"] not in allowed or item["source_id"] not in jets:
            raise ValueError("Unknown jet shape parameter")
        value = initial["discrete_jets"][jets[item["source_id"]]][item["field"]]
        lo, hi = item["lower"], item["upper"]
        if (not np.isfinite([lo, hi, value]).all() or not lo < hi
                or not lo <= value <= hi):
            raise ValueError("Finite explicit bounds must contain initial value")
        if item["field"] == "temperature_K" and lo <= 0:
            raise ValueError("Temperature lower bound must be positive")
        if item["field"] == "bulk_speed_m_s" and lo < 0:
            raise ValueError("Speed lower bound must be nonnegative")
        if item["field"] == "tilt_deg" and not 0 <= lo < hi <= 90:
            raise ValueError("Tilt bounds must remain in [0,90]")
        x0.append(value)
        lower.append(lo)
        upper.append(hi)

    def evaluate(values):
        candidate = copy.deepcopy(initial)
        for item, value in zip(definitions, values):
            candidate["discrete_jets"][jets[item["source_id"]]][item["field"]] = float(value)
        matrix, labels = matrix_from_config(candidate, points)
        return candidate, labels, solve_amplitudes(matrix, density, sigma)

    if definitions:
        optimum = least_squares(
            lambda values: evaluate(values)[2]["residual"], x0,
            bounds=(lower, upper), x_scale="jac", ftol=1e-10, xtol=1e-10,
            gtol=1e-10, max_nfev=1000)
        if not optimum.success:
            raise RuntimeError("Shape optimization failed: "+optimum.message)
        best, labels, result = evaluate(optimum.x)
        result["shape_jacobian_rank"] = int(np.linalg.matrix_rank(optimum.jac))
        result["shape_at_bound"] = optimum.active_mask
    else:
        best, labels, result = evaluate([])

    result["config"] = best
    result["component_labels"] = labels
    result["component_amplitudes"] = dict(zip(labels, result["coefficients"]))
    source_count = sum(label.startswith(("jet:", "stripe:")) for label in labels)
    result["source_rate_s"] = {
        label: float(result["coefficients"][index]*1e28)
        for index, label in enumerate(labels[:source_count])}
    background = best.get("background", {})
    if "background:density_at_reference_radius_m-3" in labels:
        index = labels.index("background:density_at_reference_radius_m-3")
        exponent = background.get("power_law_exponent", -2.0)
        result["background_power_law_coefficient"] = float(
            result["coefficients"][index]
            / best["reference_radius_m"]**exponent)
    result["nominal_dof"] = len(density)-result["rank"]-len(definitions)
    result["nominal_reduced_chi2"] = (
        result["chi2"]/result["nominal_dof"]
        if result["nominal_dof"] > 0 else None)
    result["uncertainty_status"] = (
        "input sigmas used only as weights; source covariance is not reported "
        "because nonnegative bounds and ill-conditioning invalidate a simple "
        "linear covariance interpretation")
    return result

"""Automatic fitted-source handoff, with no implicit scientific selection.

Each observation/strategy exports its own complete source vector. Background
coefficients remain in the fit report and never become surface production.
Absolute rates use molecules/s; shapes use independent unit-sum fractions per
species. E3/E5 apparent CO2 cannot become an absolute physical CO2 boundary
without an explicitly selected fragmentation correction, so only their shapes
are emitted and the deferred absolute artifact is recorded in the index.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path

from export_forward_sources import render, sha256, validate_manifest
from source_contract import CASES, HYBRID, source_ids


def add_export_arguments(parser):
    """New CLI fits use the user-selected surface-consistent geometry.

    The library loader still preserves old JSON behavior for historical
    reproduction. Reproducing raw buried-node fits is an explicit CLI choice
    and requires disabling AMPS export because those positions violate its
    surface contract. No reference/source input file is rewritten.
    """
    parser.add_argument("--stripe-surface-mapping", choices=(
        "local_vertical_to_reference_sphere", "legacy_coordinates"),
        default="local_vertical_to_reference_sphere")
    parser.add_argument("--no-amps-export", action="store_true",
                        help="explicit historical/unsupported-law diagnostic; suppress AMPS source receipts")


def json_default(value):
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "tolist"):
        return value.tolist()
    raise TypeError("not JSON serializable: " + type(value).__name__)


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, default=json_default, allow_nan=False)
        stream.write("\n")


def resolved_geometry(config, destination):
    """Snapshot the actual fit geometry, never reload reference kinematics.

    Points must already be on the declared AMPS sphere within 1 micrometre.
    They are serialized unchanged; AMPS forms the same trapezoidal nodal
    weights from those points. The unprojected authority is retained solely
    as hashed provenance, not used to shift the fit after optimization.
    """
    snapshot = copy.deepcopy(config)
    case = snapshot.get("source_variant", snapshot["source_case"])
    if case not in CASES:
        raise ValueError("AMPS does not support source case " + str(case))
    snapshot["source_case"] = case
    snapshot["_resolved_export_geometry"] = True
    # The published catalog uses 360 degrees for porco_085's azimuth.
    # 360 and 0 are the same local tangent direction, whereas the native
    # receipt deliberately uses the half-open [0,360) representation. Record
    # this coordinate-notation conversion; never change a physical angle or
    # silently wrap an out-of-range catalog value.
    endpoint_conversions = []
    for jet in snapshot["discrete_jets"]:
        for key in ("azimuth_deg", "west_longitude_deg"):
            if jet[key] == 360.0:
                jet[key] = 0.0
                endpoint_conversions.append({"source":jet["id"], "key":key,
                                             "fit_value":360.0,"export_value":0.0})
    snapshot["_angle_endpoint_canonicalization"] = endpoint_conversions
    radius = float(snapshot["reference_radius_m"])
    if snapshot.get("source_law") != "tenishev2010_point":
        raise ValueError("AMPS export currently requires the Tenishev point law")
    for source in snapshot["discrete_jets"] + snapshot["tiger_stripes"]:
        if source.get("source_law", snapshot["source_law"]) != "tenishev2010_point":
            raise ValueError("per-source law differs from the supported AMPS law")
    for stripe in snapshot["tiger_stripes"]:
        for vertex in stripe["points_m"]:
            norm = math.sqrt(math.fsum(float(x)**2 for x in vertex))
            if not math.isfinite(norm) or abs(norm-radius) > 1e-6:
                raise ValueError("fit/export stripe geometry differs: fit projected sources or use --no-amps-export")
    if snapshot["tiger_stripes"]:
        original = Path(snapshot["_stripe_geometry_path"])
        retained = destination / "stripe_geometry_authority.json"
        # Exclusive creation protects the exact raw authority used by this
        # output bundle. It is not a new source fit or an altered observation.
        with retained.open("xb") as stream:
            stream.write(original.read_bytes())
        snapshot["_stripe_geometry_path"] = str(retained.resolve())
        snapshot["stripe_geometry_file"] = str(retained.resolve())
    path = destination / "resolved_geometry.json"
    write_json(path, snapshot)
    return case, path


def export_vectors(config, records, destination, provenance):
    """Write all declared vectors plus an auditable per-artifact index.

    records contain (strategy,target,species,rates,interpretation). All canonical
    labels are required, including exact zeros. No preferred fit, total-rate
    averaging, regularization, or CO2 fragmentation value is chosen here.
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    case, geometry = resolved_geometry(config, destination)
    index = {"schema_version": 1, "source_case": case,
             "geometry_file": geometry.name, "geometry_sha256": sha256(geometry),
             "qualification": "conditional fitted source; not a unique measured vent map",
             "provenance": provenance, "artifacts": []}
    for strategy, target, species, rates, interpretation in records:
        if set(rates) != set(source_ids(case)):
            raise ValueError("fit has incomplete or noncanonical physical source labels")
        if not all(math.isfinite(float(v)) and v >= 0 for v in rates.values()):
            raise ValueError("fit source rates must be finite and nonnegative")
        total = math.fsum(float(v) for v in rates.values())
        item = {"strategy": strategy, "target": target, "species": species,
                "interpretation": interpretation, "reported_source_rate_s-1": total,
                "files": {}}
        if total <= 0:
            item.update(status="DEFERRED", reason="no positive fitted physical-source coefficient")
            index["artifacts"].append(item)
            continue
        modes = ("shape",) if interpretation == "apparent_CO2" else ("absolute", "shape")
        if interpretation == "apparent_CO2":
            item["absolute_export_status"] = "DEFERRED_USER_FRAGMENTATION_CHOICE"
            item["absolute_export_reason"] = "E3/E5 apparent mass-44 production is not physical CO2; no factor in [2,5] selected"
        for mode in modes:
            stem = "__".join((strategy, target, species, mode))
            strengths = {label: float(value)/total if mode == "shape" else float(value)
                         for label, value in rates.items()}
            manifest = {"schema_version": 1, "source_case": case,
                "distribution_name": case+"__"+stem, "mode": mode,
                "species": {species: strengths}, "fit_provenance": {
                    **provenance, "strategy": strategy, "target": target,
                    "interpretation": interpretation,
                    "source_total_includes_background": False,
                    "geometry_matches_fit": True,
                    "conditional_shape_only": interpretation == "apparent_CO2"}}
            validate_manifest(manifest)
            manifest_path = destination / (stem+".manifest.json")
            write_json(manifest_path, manifest)
            output = destination / (stem+".in")
            text = render(manifest_path, geometry)
            with output.open("x", encoding="utf-8", newline="\n") as stream:
                stream.write(text)
            item["files"][mode] = {"source_file": output.name,
                "source_sha256": sha256(output), "manifest": manifest_path.name,
                "manifest_sha256": sha256(manifest_path)}
        item["status"] = "EXPORTED"
        index["artifacts"].append(item)
    write_json(destination / "index.json", index)
    return index


def export_study(summary, config, output_directory, species):
    """Expand each fitted shape and target total into its actual component rates."""
    strategy_key = "source_rate_shape_strategies" if species == "H2O" else "source_rate_strategies"
    records = []
    for strategy, result in summary[strategy_key].items():
        for target, fit in result["target_fits"].items():
            total = float(fit["total_source_rate_s"])
            rates = {key: total*float(fraction)
                     for key, fraction in fit["source_shape_fractions"].items()}
            interpretation = "apparent_CO2" if species == "CO2" and target in ("E3", "E5") else species
            records.append((strategy, target, species, rates, interpretation))
    result_file = Path(output_directory) / "study_results.json"
    return export_vectors(config, records, Path(output_directory) / "amps_sources", {
        "study_result": str(result_file.resolve()), "study_sha256": sha256(result_file),
        "backgrounds_excluded": True,
        "limitations": summary.get("interpretation_limits", [])})


def export_single(result, result_path, target):
    """The single-profile CLI retains exact per-component rates from its fit."""
    path = Path(result_path)
    symbol = result["config"]["molecular_species"]
    interpretation = ("apparent_CO2" if symbol == "CO2" and str(target).upper() in ("E3","E5")
                      else symbol)
    return export_vectors(result["config"], [(
        "independent_each_observation", target, symbol, result["source_rate_s"], interpretation)],
        path.parent / (path.stem+"_amps_sources"),
        {"fit_result": str(path.resolve()), "fit_sha256": sha256(path),
         "qualification": result.get("qualification", "conditional fit")})

"""Reviewed observer requests, frozen second-event protocols and MHD comparisons.

These are offline comparison/provenance contracts. They do not create a runtime
MHD provider or label a complete comparison as truth validation. Event inputs
and spacecraft identities come exclusively from the supplied records.
"""
from __future__ import annotations
import math
import re
from datetime import datetime
from .core import require, finite, freeze_record, verify_frozen, digest, gaussian, validate_asset_roles, covariance


def field_line_requests(ephemeris_asset, requests, reviewer):
    m = ephemeris_asset["metadata"]
    require(m["data_use_role"] == "construction", "line seeds must be construction ephemeris")
    require(ephemeris_asset["product"]["schema"] == "sep-observer-ephemeris-v1", "wrong ephemeris schema")
    require(isinstance(reviewer, str) and bool(reviewer), "observer requests need an explicit review authority")
    times = ephemeris_asset["product"]["time_s"]
    records = []
    identities = set()
    for request in requests:
        require(request["stable_line_id"] not in identities, "duplicate stable field-line request")
        identities.add(request["stable_line_id"])
        require(request["observer_id"] and request["time_s"] in times, "observer epoch must be explicitly covered (no extrapolation)")
        index = times.index(request["time_s"])
        position = ephemeris_asset["product"]["position_m"][index]
        require(request["solar_radius_m"] > 0 and request["outer_radius_m"] > math.sqrt(sum(x*x for x in position)) > request["solar_radius_m"], "observer seed outside declared trace domain")
        require(request["nominal_step_m"] > 0 and request["maximum_steps_per_branch"] > 0 and request["unsigned_magnetic_flux_wb"] >= 0, "invalid trace controls/measure")
        records.append(dict(request, seed_m=position,
                            position_velocity_covariance=ephemeris_asset["product"]["covariance"][index]))
    require(records, "empty observer request list")
    return freeze_record({"schema": "sep-reviewed-field-line-requests-v1",
                          "frame": ephemeris_asset["product"]["transform"]["target_frame"],
                          "epoch_utc": m["epoch_utc"], "ephemeris_identity": digest(ephemeris_asset),
                          "reviewed_by": reviewer, "requests": sorted(records, key=lambda r: r["stable_line_id"])})

TRANSFER_AUTHORITIES = {"equations", "source_efficiency_form", "spectrum_form", "mean_free_path_form",
                       "inference_procedures", "reference_surface_policy", "return_policy", "background_methodology"}

def freeze_transfer_protocol(protocol):
    require(protocol["schema"] == "sep-event-transfer-protocol-v1", "unsupported transfer protocol")
    require(set(protocol["frozen_authorities"]) == TRANSFER_AUTHORITIES and
            all(protocol["frozen_authorities"].values()), "incomplete transferable physics/inference freeze")
    require(protocol["calibration_event"] != protocol["transfer_event"], "second-event validation reuses calibration event")
    require(protocol["held_out_group"] and protocol["event_specific_inputs"], "missing event/group holdout and permitted event inputs")
    require(len(protocol["event_specific_inputs"]) == len(set(protocol["event_specific_inputs"])), "duplicate event-specific inference authority")
    require(set(protocol["transferable_constants"]).isdisjoint(protocol["event_specific_inputs"]), "same field declared transferable and event specific")
    require(protocol["freeze_before_withheld_sep"] is True, "withheld SEP was loaded before protocol freeze")
    require(protocol["front_reconstruction_independent"] is True, "second event needs independent front reconstruction")
    require(protocol["observer_paths"], "missing transfer observer paths")
    for path in protocol["observer_paths"]:
        require(set(path) >= {"observer_id", "icme_screening_asset_sha256", "background_coverage_complete", "sep_coverage_complete", "prior_magnetic_cloud", "cloud_background_validated"}, "incomplete transfer path screening")
        require(isinstance(path["icme_screening_asset_sha256"], str) and
                re.fullmatch(r"[a-f0-9]{64}", path["icme_screening_asset_sha256"]) is not None and
                path["background_coverage_complete"] is True and path["sep_coverage_complete"] is True,
                "transfer background/SEP coverage or checksummed ICME screening missing")
        require(isinstance(path["prior_magnetic_cloud"], bool) and isinstance(path["cloud_background_validated"], bool),
                "cloud screening/qualification must be explicit boolean decisions")
    prior_cloud = any(path["prior_magnetic_cloud"] and not path["cloud_background_validated"] for path in protocol["observer_paths"])
    stress = prior_cloud or (protocol["transfer_event"] == "2012-05-17" and
                            not all(path["cloud_background_validated"] for path in protocol["observer_paths"]))
    require(not stress or protocol["classification"] == "stress-test", "unrepresented prior magnetic cloud requires declared stress test")
    require(protocol["classification"] in {"primary-transfer", "stress-test", "fixed-parameter-stress-test"}, "unknown transfer classification")
    return freeze_record(protocol)


def bind_transfer_run(frozen_protocol, event_inputs, release_authorities):
    verify_frozen(frozen_protocol)
    require(set(event_inputs) == set(frozen_protocol["event_specific_inputs"]), "unregistered or missing event-specific transfer inputs")
    for name, asset in event_inputs.items():
        require(isinstance(asset, dict) and set(asset) == {"asset_id", "content_sha256", "inference_procedure"}, "event input must bind an immutable asset and frozen inference")
        require(asset["asset_id"] and re.fullmatch(r"[a-f0-9]{64}", asset["content_sha256"]) is not None, "missing event-input content checksum")
        require(asset["inference_procedure"] == frozen_protocol["frozen_authorities"]["inference_procedures"], "event input uses an unregistered inference procedure")
    require(release_authorities["transferable_constants"] == frozen_protocol["transferable_constants"], "post hoc change to transferable calibration constants")
    require(release_authorities["frozen_authorities"] == frozen_protocol["frozen_authorities"], "transfer changes frozen equations/forms/inference")
    # This is per-run identity. Reusing campaign-D's fingerprint after changing
    # event-bound inputs is specifically forbidden by the specification.
    record = {"schema": "sep-transfer-run-v1", "protocol_identity": frozen_protocol["identity"],
              "event_inputs": event_inputs, "release_authorities": release_authorities}
    record["release_calibration_fingerprint"] = digest(record)
    return freeze_record(record)


MATCHED_FIELDS = {"epoch_utc", "frame", "units", "variable_definitions", "cadence_s", "spatial_support",
                  "interpolation", "masks", "front_history_sha256", "composition", "equation_of_state"}
BOUNDARY_FIELDS = {"magnetogram_sha256", "preprocessing_sha256", "magnetic_normalization", "open_flux_normalization"}
VARIABLES = {"number_density_m3", "mass_density_kg_m3", "temperature_K", "pressure_Pa", "magnetic_field_T",
             "velocity_m_per_s", "alfven_speed_m_per_s", "fast_speed_m_per_s", "signed_open_flux_wb",
             "unsigned_open_flux_wb", "D2_first_fast_height_m", "D2_first_supercritical_height_m"}


def compare_backgrounds(analytic, mhd, allow_combined_discrepancy=False):
    for sample in (analytic, mhd):
        verify_frozen(sample)
        require(sample["schema"] == "sep-offline-background-samples-v1", "unsupported matched sample schema")
        require(MATCHED_FIELDS|BOUNDARY_FIELDS <= set(sample["metadata"]), "missing matched-boundary/sample metadata")
        m = sample["metadata"]
        # A nonempty label is not checksum ownership. These bindings must be
        # actual SHA-256 digests before matching or permitting a combined
        # discrepancy; that flag cannot rescue missing/ambiguous provenance.
        for key in ("magnetogram_sha256", "preprocessing_sha256", "front_history_sha256"):
            require(isinstance(m[key], str) and re.fullmatch(r"[a-f0-9]{64}", m[key]) is not None,
                    "invalid comparison authority checksum: " + key)
        for key, unit in (("magnetic_normalization", "T"), ("open_flux_normalization", "Wb")):
            normalization = m[key]
            require(isinstance(normalization, dict) and set(normalization) == {"value", "units", "definition"} and
                    normalization["units"] == unit and normalization["definition"] and
                    isinstance(normalization["value"], (int, float)) and normalization["value"] > 0 and
                    math.isfinite(normalization["value"]), "missing absolute SI comparison normalization: " + key)
        try:
            datetime.strptime(m["epoch_utc"], "%Y-%m-%dT%H:%M:%SZ")
        except (ValueError, TypeError) as error:
            raise ValueError("comparison epoch must be an explicit valid UTC instant") from error
        require(m["frame"] and isinstance(m["cadence_s"], (int, float)) and m["cadence_s"] > 0, "comparison frame/cadence missing")
        require(m["front_history_sha256"] and m["interpolation"] in {"none-analytic-at-nodes", "bounded-linear", "original-node-values"}, "unsupported/missing comparison front or interpolation")
        require(isinstance(m["masks"], list) and len(m["masks"]) == len(sample["sample_coordinates"]) and all(isinstance(v, bool) for v in m["masks"]) and any(m["masks"]), "comparison mask/support is empty or incomplete")
        support = m["spatial_support"]
        require(set(support) == {"minimum_m", "maximum_m", "time_s"} and len(support["minimum_m"]) == len(support["maximum_m"]) == 3 and len(support["time_s"]) == 2, "incomplete comparison spatial/time support")
        finite(support["minimum_m"]+support["maximum_m"]+support["time_s"])
        require(all(a <= b for a, b in zip(support["minimum_m"], support["maximum_m"])) and support["time_s"][0] <= support["time_s"][1], "reversed sample support")
        for point in sample["sample_coordinates"]:
            require(len(point) == 4 and all(a <= x <= b for a, x, b in zip(support["minimum_m"], point[:3], support["maximum_m"])) and support["time_s"][0] <= point[3] <= support["time_s"][1], "comparison interpolated/extrapolated outside declared support")
        require(sample["residual_tuning_used"] is False, "comparison cannot tune a background from residuals")
        require(sample["source_uri"] and isinstance(sample["source_sha256"], str) and
                re.fullmatch(r"[a-f0-9]{64}", sample["source_sha256"]) is not None and
                sample["model_version"] and sample["run_id"], "missing checksummed offline model/run/source provenance")
        require(set(sample["variables"]) == VARIABLES, "incomplete offline n/B/u/vA/cf/open-flux/D2 comparison")
    differences = sorted(key for key in MATCHED_FIELDS if analytic["metadata"][key] != mhd["metadata"][key])
    require(not differences, "frame/epoch/units/support/definition/mask mismatch: " + ", ".join(differences))
    boundary_mismatch = sorted(key for key in BOUNDARY_FIELDS if analytic["metadata"][key] != mhd["metadata"][key])
    require(not boundary_mismatch or allow_combined_discrepancy, "unreported boundary/normalization discrepancy")
    definitions = analytic["metadata"]["variable_definitions"]
    units = analytic["metadata"]["units"]
    require(isinstance(definitions, dict) and set(definitions) == VARIABLES and all(definitions.values()), "missing variable definitions")
    expected_units = {"number_density_m3": "m^-3", "mass_density_kg_m3": "kg/m^3", "temperature_K": "K", "pressure_Pa": "Pa",
        "magnetic_field_T": "T", "velocity_m_per_s": "m/s", "alfven_speed_m_per_s": "m/s", "fast_speed_m_per_s": "m/s",
        "signed_open_flux_wb": "Wb", "unsigned_open_flux_wb": "Wb", "D2_first_fast_height_m": "m", "D2_first_supercritical_height_m": "m"}
    require(units == expected_units, "comparison variable units differ from the declared SI definitions")
    require(analytic["sample_coordinates"] == mhd["sample_coordinates"] and analytic["sample_coordinates"], "unmatched coordinate/time sampling")
    count = len(analytic["sample_coordinates"])
    finite(v for coordinate in analytic["sample_coordinates"] for v in coordinate)
    metrics = {}
    for variable in sorted(VARIABLES):
        a, b = analytic["variables"][variable], mhd["variables"][variable]
        require(len(a["values"]) == len(b["values"]), "sample variable dimensions differ")
        expected = count*3 if variable in {"magnetic_field_T", "velocity_m_per_s"} else count
        require(len(a["values"]) == expected, "variable does not cover declared spatial support")
        finite(a["values"]+b["values"])
        indices = [i for i in range(expected) if analytic["metadata"]["masks"][i//3 if expected == count*3 else i]]
        residual = [a["values"][i]-b["values"][i] for i in indices]
        covariance(a["covariance"], expected); covariance(b["covariance"], expected)
        c = [[a["covariance"][i][j]+b["covariance"][i][j] for j in indices] for i in indices]
        score = gaussian(residual, c)
        metrics[variable] = dict(score, mean_difference=math.fsum(residual)/len(residual),
                                 rms_difference=math.sqrt(math.fsum(x*x for x in residual)/len(residual)),
                                 compared_values=len(indices), masked_values=expected-len(indices))
    topology = {}
    for key in ("topology", "connectivity"):
        require(key in analytic and key in mhd, "missing topology/connectivity comparability declaration")
        a, b = analytic[key], mhd[key]
        if a["definition"] != b["definition"]:
            topology[key] = {"applicable": False, "reason": "different categorical definitions"}
        else:
            require(len(a["values"]) == len(b["values"]) == count, "categorical coverage mismatch")
            topology[key] = {"applicable": True, "agreement": sum(x == y for x, y, valid in zip(a["values"], b["values"], analytic["metadata"]["masks"]) if valid)/sum(analytic["metadata"]["masks"])}
    return freeze_record({"schema": "sep-offline-mhd-comparison-v1",
        "classification": "combined-boundary-plus-model-discrepancy" if boundary_mismatch else "matched-structural-cross-model-comparison",
        "boundary_mismatches": boundary_mismatch, "analytic_identity": digest(analytic), "mhd_identity": digest(mhd),
        "metrics": metrics, "categorical_metrics": topology,
        "truth_validation": False, "runtime_imported_provider_qualified": False, "residual_tuning_used": False})

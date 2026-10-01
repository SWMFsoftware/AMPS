"""Preregistered Cartesian formation-height comparisons, with SEP held out.

Candidate scoring uses independent qualification observations. The exact
candidate product, density-conditioned radio likelihood, weights, gates and
source identities are frozen together. No candidate may disappear after
inspection, and no universal formation-height window is used as a shortcut.
"""
from __future__ import annotations
import itertools
import math
import re
from .core import (require, finite, digest, freeze_record, verify_frozen,
                   validate_asset_roles, gaussian, covariance)
from .inference import evaluate_power_law

AXES = ("magnetogram", "field_scale", "coupling_radii", "wind_density", "front")

def preregister_campaign(record, assets):
    require(set(record) == {"schema", "axes", "hard_gates", "constraint", "weighting_rule", "density_construction_procedure", "density_construction_equations_sha256"}, "unknown/incomplete preregistration field or selector")
    require(record["density_construction_procedure"] and isinstance(record["density_construction_equations_sha256"], str) and
            re.fullmatch(r"[a-f0-9]{64}", record["density_construction_equations_sha256"]) is not None,
            "candidate density construction equations/procedure must be frozen and checksummed")
    require(record["schema"] == "sep-formation-preregistration-v1", "unsupported campaign schema")
    validate_asset_roles(assets)
    require(set(record["axes"]) == set(AXES), "candidate axes must be the complete five-factor product")
    known = {a["metadata"]["asset_id"]: a for a in assets}
    for axis in AXES:
        members = record["axes"][axis]
        require(bool(members) and len({m["id"] for m in members}) == len(members), "empty/duplicate Cartesian axis")
        for member in members:
            require(set(member) == {"id", "asset_id", "prior_weight"}, "unknown candidate member field or SEP-derived selector")
            require(member["asset_id"] in known, "candidate construction asset is absent")
            m = known[member["asset_id"]]["metadata"]
            require(m["data_use_role"] == "construction" and "sep" not in m["kind"].lower(), "candidate selector must use independent non-SEP construction inputs")
            require(member["prior_weight"] > 0 and math.isfinite(member["prior_weight"]), "invalid preregistered prior")
    require(set(record["hard_gates"]) == {"D6_maximum_relative_error", "topology_minimum_agreement", "D1_maximum_relative_error"}, "unknown or SEP/universal-height candidate gate")
    finite(record["hard_gates"].values())
    require(0 <= record["hard_gates"]["topology_minimum_agreement"] <= 1 and
            record["hard_gates"]["D6_maximum_relative_error"] >= 0 and record["hard_gates"]["D1_maximum_relative_error"] >= 0,
            "invalid preregistered diagnostic bounds")
    require(set(record["constraint"]) == {"kind", "asset_id"} and
            record["constraint"]["kind"] in {"radio-frequency-time", "joint-preinferred-height"}, "SEP output is forbidden as a candidate selector")
    constraint = known.get(record["constraint"]["asset_id"])
    require(constraint is not None and constraint["metadata"]["data_use_role"] == "qualification" and
            "sep" not in constraint["metadata"]["kind"].lower(), "constraint must be independent non-SEP qualification data")
    require(record["weighting_rule"] == "preregistered-prior-times-joint-likelihood", "unsupported/post-hoc weighting rule")
    # Bind every observation and inference product, not only its ID/path. The
    # preregistration does not contain any realized SEP values or intensities.
    return freeze_record(dict(record, asset_identities={key: digest(value) for key, value in sorted(known.items())}))

def _logsumexp(values):
    maximum = max(values)
    return maximum+math.log(math.fsum(math.exp(x-maximum) for x in values))

def radio_likelihood(observation, profile, front_radii):
    frequency = observation["frequency_hz"]
    finite(frequency)
    require(all(x > 0 for x in frequency), "radio frequency must be positive")
    require(len(frequency) == len(front_radii), "radio/front cadence mismatch")
    density, density_cov = evaluate_power_law(profile, front_radii)
    # f_pe=(1/2pi)sqrt(n e^2/(epsilon0 m_e)); SI density is m^-3.
    epsilon0 = 1/(1.25663706212e-6*299792458.0**2)
    factor = math.sqrt(1.602176634e-19**2/(epsilon0*9.1093837015e-31))/(2*math.pi)
    base = [factor*math.sqrt(n) for n in density]
    covariance(observation["frequency_covariance_hz2"], len(frequency))
    branches = observation["harmonic_hypotheses"]
    require(set(h["harmonic"] for h in branches) == {1, 2} and len(branches) == 2, "fundamental/harmonic hypotheses must both be explicit")
    require(all(h["prior_probability"] > 0 for h in branches) and
            abs(math.fsum(h["prior_probability"] for h in branches)-1) < 1e-12, "radio harmonic priors must normalize")
    likelihoods = []
    cross = observation.get("density_frequency_cross_covariance", [[0.0]*len(frequency) for _ in frequency])
    require(len(cross) == len(frequency) and all(len(row) == len(frequency) for row in cross), "radio/density cross covariance dimensions")
    finite(v for row in cross for v in row)
    for hypothesis in branches:
        h = hypothesis["harmonic"]
        predicted = [h*f for f in base]
        jacobian = [h*f/(2*n) for f, n in zip(base, density)]
        # Residual covariance includes density uncertainty and any declared
        # covariance with radio inference; it never treats a reused density
        # inversion as an independent observed height.
        c = [[observation["frequency_covariance_hz2"][i][j]+jacobian[i]*density_cov[i][j]*jacobian[j]
              -jacobian[i]*cross[i][j]-cross[j][i]*jacobian[j]
              for j in range(len(frequency))] for i in range(len(frequency))]
        score = gaussian([a-b for a, b in zip(frequency, predicted)], c)
        likelihoods.append(dict(score, harmonic=h, prior_probability=hypothesis["prior_probability"],
                                predicted_frequency_hz=predicted, density_at_front_m3=density))
    return {"log_likelihood": _logsumexp([item["log_likelihood"]+math.log(item["prior_probability"]) for item in likelihoods]),
            "branches": likelihoods, "density_conditioned": True}

def height_likelihood(observation, profile, front_radii, density_checksum):
    require(observation["density_model_sha256"] == density_checksum, "height inversion density-model checksum differs from candidate authority")
    require(observation["independence"] == "joint-density-inference", "pre-inferred height cannot be independent across shared density members")
    observed = observation["height_m"]
    require(len(observed) == len(front_radii), "height/front cadence mismatch")
    # The density normalization is an additional observed inference coordinate.
    # Its covariance with every inferred height is retained in the full joint
    # likelihood. Candidates share this one constraint, not independent copies.
    residual = [a-b for a, b in zip(observed, front_radii)]+[
        observation["inferred_log_reference_density"]-profile["log_reference_value"]]
    score = gaussian(residual, observation["joint_height_log_density_covariance"])
    return dict(score, density_conditioned=True, independence="joint-density-inference")

def select_candidates(preregistration, realized, assets):
    verify_frozen(preregistration)
    validate_asset_roles(assets)
    known = {a["metadata"]["asset_id"]: a for a in assets}
    require(preregistration["asset_identities"] == {key: digest(value) for key, value in sorted(known.items())}, "assets changed after candidate preregistration")
    product = list(itertools.product(*[[m["id"] for m in preregistration["axes"][axis]] for axis in AXES]))
    tuples = [tuple(row["tuple"]) for row in realized]
    require(len(tuples) == len(set(tuples)) and set(tuples) == set(product), "missing/duplicated/unregistered Cartesian candidate tuple")
    rows_by_tuple = {tuple(row["tuple"]): row for row in realized}
    indices = {axis: {m["id"]: m for m in preregistration["axes"][axis]} for axis in AXES}
    constraint = known[preregistration["constraint"]["asset_id"]]
    observation = constraint["product"]
    rows = []
    for key in product:
        original = rows_by_tuple[key]
        require(set(original) == {"tuple", "D6", "topology", "D1", "D2", "qualification_assets", "front_radii_m", "density_construction"}, "unknown/post-hoc or SEP-dependent realized candidate field")
        for diagnostic in ("D6", "topology", "D1", "D2"):
            require(diagnostic in original["qualification_assets"], "missing independent diagnostic evidence")
            asset = known.get(original["qualification_assets"][diagnostic])
            require(asset is not None and asset["metadata"]["data_use_role"] == "qualification" and "sep" not in asset["metadata"]["kind"].lower(), "construction/withheld/SEP asset reused for candidate qualification")
        finite([original["D6"], original["topology"], original["D1"]]+original["front_radii_m"])
        require(original["D6"] >= 0 and original["D1"] >= 0 and 0 <= original["topology"] <= 1, "invalid independent diagnostic value")
        require(set(original["D2"]) == {"first_fast_height_m", "first_fast_time_s", "first_supercritical_height_m", "first_supercritical_time_s"}, "D2 must distinguish fast/supercritical heights and times")
        finite(original["D2"].values())
        require(len(original["front_radii_m"]) == len(observation["time_s"]), "event constraint time support differs from front")
        wind = known[indices["wind_density"][key[3]]["asset_id"]]
        require(wind["product"]["schema"] == "sep-positive-power-law-profile-v1" and wind["metadata"]["units"] == "m^-3", "radio/height likelihood needs an explicit density profile")
        construction = original["density_construction"]
        require(set(construction) == {"tuple", "procedure", "equations_sha256", "source_checksums", "profile"}, "incomplete per-tuple density construction")
        require(tuple(construction["tuple"]) == key and construction["procedure"] == preregistration["density_construction_procedure"] and
                construction["equations_sha256"] == preregistration["density_construction_equations_sha256"], "candidate density belongs to another tuple/procedure")
        expected_sources = {axis: known[indices[axis][member]["asset_id"]]["metadata"]["source_sha256"] for axis, member in zip(AXES, key)}
        require(construction["source_checksums"] == expected_sources, "candidate density does not bind every construction authority")
        profile = construction["profile"]
        require(profile["schema"] in {"sep-positive-power-law-profile-v1", "sep-front-sampled-density-v1"}, "unsupported realized density reconstruction")
        if profile["schema"] == "sep-positive-power-law-profile-v1":
            finite([profile["reference_radius_m"], profile["log_reference_value"], profile["exponent"]]+profile["support_m"])
            require(profile["reference_radius_m"] > 0 and len(profile["support_m"]) == 2 and 0 < profile["support_m"][0] < profile["support_m"][1], "invalid realized density support")
            covariance(profile["log_parameter_covariance"], 2)
        else:
            require(preregistration["constraint"]["kind"] == "radio-frequency-time", "joint height inference requires an explicit density normalization coordinate")
            require(profile["geometry_uncertainty_model"] in {"included-in-covariance", "preregistered-front-ensemble"}, "sampled density lacks front/projection uncertainty authority")
        # Each tuple supplies its fully rebuilt density product. In particular,
        # changing field scale/radii is not assumed to leave density unchanged
        # or to scale Alfvén speed linearly. The external construction equations
        # and all five member sources are bound before the likelihood is used.
        if preregistration["constraint"]["kind"] == "radio-frequency-time":
            score = radio_likelihood(observation, profile, original["front_radii_m"])
        else:
            score = height_likelihood(observation, profile, original["front_radii_m"], wind["metadata"]["source_sha256"])
        gates = preregistration["hard_gates"]
        accepted = (original["D6"] <= gates["D6_maximum_relative_error"] and
                    original["topology"] >= gates["topology_minimum_agreement"] and
                    original["D1"] <= gates["D1_maximum_relative_error"])
        log_prior = math.fsum(math.log(indices[axis][member]["prior_weight"]) for axis, member in zip(AXES, key))
        rows.append(dict(original, accepted=accepted, likelihood=score, log_prior=log_prior))
    eligible = [row["log_prior"]+row["likelihood"]["log_likelihood"] for row in rows if row["accepted"]]
    normalizer = _logsumexp(eligible) if eligible else None
    for row in rows:
        row["weight"] = math.exp(row["log_prior"]+row["likelihood"]["log_likelihood"]-normalizer) if row["accepted"] else 0.0
    return freeze_record({"schema": "sep-formation-selection-v1", "preregistration_identity": preregistration["identity"],
                          "status": "qualified-selection" if eligible else "no-qualified-candidate",
                          "realized_source_identity": digest(realized), "rows": rows,
                          "withheld_sep_used": False})

def load_withheld_after_freeze(selection, asset):
    verify_frozen(selection)
    require(selection["schema"] == "sep-formation-selection-v1" and selection["status"] == "qualified-selection", "selection must be frozen and qualified before withheld validation")
    require(asset["metadata"]["data_use_role"] == "withheld-validation", "requested product is not withheld validation")
    # A caller may now compute all registered validation metrics. The frozen
    # selector is immutable; these values are never returned to its weights.
    return asset["product"]

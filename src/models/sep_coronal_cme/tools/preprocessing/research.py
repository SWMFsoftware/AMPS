"""Version-6 research assets and explicitly bounded offline wave/campaign solvers.

These products are separate from schema-5 event inputs. Manufactured numerical
verification is not an observed-event or production-host qualification. The
wave producer is a named 1-D, two-moment, small-amplitude Alfvén approximation;
it carries its approximation and every coupling/energy-sink authority.
"""
from __future__ import annotations
import copy
import math
import re
from .core import require, finite, digest, freeze_record, verify_frozen, covariance

CAPABILITIES = {"integrated-peclet-family", "return-renewal", "foreshock-distance-proxy",
                "self-generated-waves", "smooth-field-drift", "dynamic-attitude",
                "impulsive-source", "wind-envelope", "transfer-campaign", "mhd-comparison",
                "nonradial-winding", "streaming-limit"}


def content_hash(value):
    require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value), "invalid content SHA-256")
    return value


def research_configuration(capabilities, authorities):
    require(set(capabilities) == CAPABILITIES and all(type(v) is bool for v in capabilities.values()),
            "unknown/missing/nonboolean independent research capability")
    for name, enabled in capabilities.items():
        if enabled:
            require(name in authorities and set(authorities[name]) >= {"algorithm_version", "verification_gate", "domain", "validation_claim"},
                    "enabled research capability lacks version/gate/domain/claim")
            require(all(authorities[name].values()), "empty active research authority")
    return freeze_record(dict(schema="sccm-research-configuration-v6", capability_flags=capabilities,
                              authorities=authorities, schema5_selector_available=False,
                              production_adapter_qualified=False))


WAVE_COMMON_BINDINGS = {"producer_configuration", "producer_qualification", "background_generation",
    "shock_history_generation", "shock_frame", "wave_frame_signed_direction", "wave_number_sign_convention",
    "spectral_grid", "boundary_normalization", "interpolation", "coverage_mask", "resonance_mapping",
    "tolerance_profile", "residual_definitions", "residual_fields", "uncertainty", "table_content"}
WAVE_COUPLED_BINDINGS = {"source_geometry_timing", "source_spectrum_normalization", "particle_distribution",
    "species_weights", "transport_return_policy", "coupling_iteration_cadence"}


def validate_wave_asset(asset, active):
    verify_frozen(asset)
    require(asset["schema"] == "sccm-foreshock-wave-asset-v6", "wrong wave asset schema")
    bindings = asset["bindings"]
    require(WAVE_COMMON_BINDINGS|WAVE_COUPLED_BINDINGS|{"external_wave_field"} == set(bindings), "missing/extra wave replay authority")
    for key in WAVE_COMMON_BINDINGS:
        require(bindings[key] and bindings[key] != "not-applicable" and bindings[key] == active[key], "wave replay mismatch: "+key)
    for key in WAVE_COMMON_BINDINGS-{"background_generation","shock_history_generation","shock_frame","wave_frame_signed_direction","wave_number_sign_convention"}:
        content_hash(bindings[key])
    require(type(bindings["background_generation"]) is int and bindings["background_generation"]>0 and
            type(bindings["shock_history_generation"]) is int and bindings["shock_history_generation"]>0,
            "wave generation identity is missing/stale")
    if asset["coupling_kind"] == "self-consistent-particle-wave-iteration":
        require(all(bindings[k] and bindings[k] != "not-applicable" and bindings[k] == active[k] for k in WAVE_COUPLED_BINDINGS),
                "missing coupled particle/source/transport/cadence authority")
        for key in WAVE_COUPLED_BINDINGS: content_hash(bindings[key])
        require(bindings["external_wave_field"] == "not-applicable", "coupled solver cannot claim prescribed wave authority")
    else:
        require(asset["coupling_kind"] == "prescribed-external-wave-field" and bindings["external_wave_field"] == active["external_wave_field"], "unknown/missing external wave authority")
        content_hash(bindings["external_wave_field"])
        require(all(bindings[k] == "not-applicable" for k in WAVE_COUPLED_BINDINGS), "external-field replay forged coupled authority")
    require(asset["coverage_mask"] and all(type(x) is bool for x in asset["coverage_mask"]), "wave coverage mask missing")
    require(all(asset["coverage_mask"]), "requested resonant cell is unsupported; no extrapolation/fallback")
    require(len(asset["coverage_mask"])==len(asset["spectrum"]),"wave mask/table shape mismatch")
    finite(asset["spectrum"]); require(all(v >= 0 for v in asset["spectrum"]), "negative wave spectrum")
    require(digest(asset["spectrum"]) == bindings["table_content"] and digest(asset["coverage_mask"]) == bindings["coverage_mask"], "wave table/mask content identity mismatch")
    return asset


def advance_resonant_waves(state, controls, dt):
    """Positive explicit FV evolution in s and log|k| with coupled reservoirs.

    W is total Alfvén energy per volume per log|k|. The two directions carry
    pseudomomentum +/-W/v_A in the declared stationary wave frame. Streaming
    growth uses gamma=(pi/4) Omega n_cr/n_i max(+/-v_stream/v_A-1,0).
    Particle energy and parallel pseudomomentum lose exactly the wave gain.
    This two-moment/quasilinear closure does not claim nonlinear kinetic
    saturation. Refraction, focusing, damping and spectral cascade are separate
    explicit authorities; background work and every sink are ledgered.
    """
    finite([dt, controls["ds_m"], controls["dlogk"], controls["alfven_speed_m_per_s"]])
    require(dt > 0 and controls["ds_m"] > 0 and controls["dlogk"] > 0 and controls["alfven_speed_m_per_s"] > 0, "invalid wave grid/cadence")
    out = copy.deepcopy(state); w = state["wave_energy"]
    ns, nk = len(w), len(w[0])
    require(ns >= 2 and nk >= 2 and all(len(row) == nk and all(len(cell) == 2 for cell in row) for row in w), "wave grid shape")
    require(len(state["particle_energy"]) == len(state["particle_momentum"]) == ns, "particle grid/cadence mismatch")
    finite(x for row in w for cell in row for x in cell)
    finite(state["particle_energy"]+state["particle_momentum"])
    require(all(x >= 0 for row in w for cell in row for x in cell), "negative initial wave energy")
    require(all(x >= 0 for x in state["particle_energy"]), "negative initial particle energy")
    va, ds, dk = controls["alfven_speed_m_per_s"], controls["ds_m"], controls["dlogk"]
    area=controls["tube_area_m2"]; finite([area]); require(area>0,"invalid wave physical volume")
    volume=ds*area
    require(controls["front_coordinate_support_m"][0] >= controls["reference_distance_m"] > 0 and
            controls["shock_propagation_approximation"] == "frozen-shock-over-coupling-step",
            "wave particle support overlaps unresolved release layer/unqualified shock-frame propagation")
    u, refract, cascade = controls["flow_speed_m_per_s"], controls["logk_advection_per_s"], controls["cascade_per_s"]
    damping, focus = controls["damping_per_s"], controls["focusing_per_s"]
    finite([u, refract, cascade, damping, focus, controls["gyrofrequency_per_s"], controls["cr_to_ion_density"], controls["particle_mass_density"]])
    require(cascade >= 0 and damping >= 0 and controls["particle_mass_density"] > 0 and controls["cr_to_ion_density"] >= 0, "invalid wave closure coefficients")
    require(dt*(abs(u)+va)/ds+dt*abs(refract)/dk+2*dt*cascade/(dk*dk)+dt*damping < 1,
            "wave positivity CFL/cadence violated")
    ledger = {"particle_wave_energy": 0.0, "particle_wave_momentum": 0.0, "damping_heat": 0.0,
              "background_work": 0.0, "spectral_boundary_sink": 0.0,
              "damping_momentum":0.0,"background_momentum_work":0.0,"spectral_boundary_momentum":0.0}
    # Spatial boundaries are explicitly periodic manufactured/tube boundaries.
    # This producer must not be reused for an open shock boundary without a new
    # boundary authority and its flux ledger. Spectral boundaries are outflow.
    require(controls["spatial_boundary"] == "periodic" and controls["spectral_boundary"] == "outflow", "unqualified wave boundary closure")
    for i in range(ns):
        streaming = state["particle_momentum"][i]/controls["particle_mass_density"]
        for j in range(nk):
            for direction, sign in enumerate((1, -1)):
                value = w[i][j][direction]; speed = u+sign*va
                spatial = -(abs(speed)*dt/ds)*(value-w[(i-1 if speed >= 0 else i+1)%ns][j][direction])
                neighbor = j-1 if refract >= 0 else j+1
                spectral = -(abs(refract)*dt/dk)*(value-(w[i][neighbor][direction] if 0 <= neighbor < nk else 0))
                if (refract > 0 and j == nk-1) or (refract < 0 and j == 0):
                    escaped=abs(refract)*dt/dk*value*dk*volume
                    ledger["spectral_boundary_sink"] += escaped
                    ledger["spectral_boundary_momentum"] += sign*escaped/va
                left = w[i][max(0, j-1)][direction]; right = w[i][min(nk-1, j+1)][direction]
                diffuse = dt*cascade/(dk*dk)*(left-2*value+right)
                gamma = 0.0 if not controls["growth_enabled"] else math.pi/4*controls["gyrofrequency_per_s"]*controls["cr_to_ion_density"]*max(sign*streaming/va-1, 0)
                gain = dt*2*gamma*value
                loss = dt*damping*value; work = dt*focus*value
                out["wave_energy"][i][j][direction] = value+spatial+spectral+diffuse+gain-loss+work
                out["particle_energy"][i] -= gain*dk
                out["particle_momentum"][i] -= sign*gain*dk/va
                ledger["particle_wave_energy"] += gain*dk*volume; ledger["particle_wave_momentum"] += sign*gain*dk*volume/va
                ledger["damping_heat"] += loss*dk*volume; ledger["background_work"] += work*dk*volume
                ledger["damping_momentum"] += sign*loss*dk*volume/va
                ledger["background_momentum_work"] += sign*work*dk*volume/va
    require(all(x >= 0 and math.isfinite(x) for row in out["wave_energy"] for cell in row for x in cell)
            and all(e >= 0 for e in out["particle_energy"]), "coupled wave/particle positivity failed; reduce cadence")
    def total_energy(data):
        return volume*(math.fsum(data["particle_energy"])+dk*math.fsum(v for row in data["wave_energy"] for cell in row for v in cell))
    def total_momentum(data):
        return volume*(math.fsum(data["particle_momentum"])+dk/va*math.fsum(cell[0]-cell[1] for row in data["wave_energy"] for cell in row))
    ledger["energy_residual_j"]=total_energy(out)-total_energy(state)-ledger["background_work"]+ledger["damping_heat"]+ledger["spectral_boundary_sink"]
    ledger["momentum_residual_kg_m_per_s"]=total_momentum(out)-total_momentum(state)-ledger["background_momentum_work"]+ledger["damping_momentum"]+ledger["spectral_boundary_momentum"]
    tolerance=controls["conservation_relative_tolerance"]; finite([tolerance]); require(tolerance>0,"missing conservation tolerance")
    # Counterpropagating waves may have zero *net* momentum. Normalize by the
    # absolute inventory, never by that cancellation, so symmetry does not
    # manufacture an impossible zero-roundoff acceptance tolerance.
    momentum_inventory=volume*(math.fsum(abs(x) for x in state["particle_momentum"])+
        dk/va*math.fsum(v for row in w for cell in row for v in cell))
    require(abs(ledger["energy_residual_j"])<=tolerance*max(abs(total_energy(state)),1e-300) and
            abs(ledger["momentum_residual_kg_m_per_s"])<=tolerance*max(momentum_inventory,1e-300),
            "particle-wave/background/sink energy-momentum accounting failed")
    out["ledger"] = ledger
    return out


def resonant_diffusion(mu, speed, gyrofrequency, field_t, logk, directional_power, broadening):
    """Declared broadened gyroresonance, bounded linear interpolation, no gap fill."""
    finite([mu, speed, gyrofrequency, field_t, broadening]); finite(logk+directional_power)
    require(-1 < mu < 1 and speed > 0 and gyrofrequency > 0 and field_t > 0 and broadening > 0,
            "invalid resonant particle/wave domain")
    require(len(logk) == len(directional_power) >= 2 and all(a < b for a,b in zip(logk,logk[1:])), "invalid spectral coordinate")
    coordinate = math.log(gyrofrequency/(speed*math.sqrt(mu*mu+broadening*broadening)))
    require(logk[0] <= coordinate <= logk[-1], "resonance outside asset support")
    for j in range(len(logk)-1):
        if logk[j] <= coordinate <= logk[j+1]:
            fraction = (coordinate-logk[j])/(logk[j+1]-logk[j])
            power = (1-fraction)*directional_power[j]+fraction*directional_power[j+1]
            require(power > 0, "zero/unsupported resonant scattering power")
            nu = math.pi/2*gyrofrequency*(4*math.pi*1e-7)*power/(field_t*field_t)
            return 0.5*nu*(1-mu*mu)
    raise ValueError("uncovered resonance")


def streaming_limit(rows, registration):
    verify_frozen(registration)
    require(registration["schema"] == "sccm-streaming-limit-registration-v6" and registration["frozen_before_sep"] is True,
            "post-hoc/missing streaming-limit protocol")
    for key in ("source_asset_sha256", "calibration_asset_sha256", "response_sha256"): content_hash(registration[key])
    require(registration["source_asset_sha256"] != registration["calibration_asset_sha256"], "reused limit calibration provenance")
    require(set(row["stratum"] for row in rows) == set(registration["required_strata"]), "incomplete streaming-limit strata")
    result = []
    radial=registration["radial_transform"]
    require(radial["law"] == "registered-power-law-sensitivity" and radial["target_radius_m"]>0 and radial["limit_reference_radius_m"]>0,
            "missing preregistered radial transform")
    content_hash(radial["authority_sha256"])
    finite([radial["exponent"],radial["target_radius_m"],radial["limit_reference_radius_m"]])
    response = registration["response_weights"]
    finite(response); require(all(v >= 0 for v in response) and sum(response) > 0, "invalid channel response")
    for row in rows:
        require(row["units"] == registration["units"] and row["window"] == registration["window"], "response/unit/window mismatch")
        require(registration["radius_support_m"][0] <= row["radius_m"] <= registration["radius_support_m"][1], "silent radial extrapolation")
        require(len(row["model"]) == len(row["limit"]) == len(response), "response channel mismatch")
        finite(row["model"]+row["limit"]+[row["variance_log_model"],row["variance_log_limit"],row["log_covariance"]])
        jm = sum(x*w for x,w in zip(row["model"],response))*(radial["target_radius_m"]/row["radius_m"])**radial["exponent"]
        jl = sum(x*w for x,w in zip(row["limit"],response))*(radial["target_radius_m"]/radial["limit_reference_radius_m"])**radial["exponent"]
        state = row["state"]
        require(state in {"valid","zero","censored","unsupported","unmapped"}, "unknown streaming-limit state")
        require(jm >= 0 and jl >= 0, "negative response-folded intensity")
        if state == "valid" and (jm == 0 or jl == 0): state = "zero"
        applicable = state == "valid" and jm > 0 and jl > 0
        variance = row["variance_log_model"]+row["variance_log_limit"]-2*row["log_covariance"]
        require(variance >= 0, "invalid propagated intensity covariance")
        q = math.log10(jm/jl) if applicable else None
        result.append(dict(row, transformed_model=jm, transformed_limit=jl, applicable=applicable,
            q_SL=q, q_interval=None if not applicable else [q-1.96*math.sqrt(variance)/math.log(10),q+1.96*math.sqrt(variance)/math.log(10)],
            state=state if not applicable else "valid", feedback_qualified=False))
    return freeze_record(dict(schema="sccm-streaming-limit-product-v6", registration_identity=registration["identity"],
                              rows=result, observational_agreement_gate=False,
                              applicable_count=sum(r["applicable"] for r in result),
                              typed_inapplicable_count=sum(not r["applicable"] for r in result)))


def family_bundle(members, calibration, target_peclet):
    """New major: a momentum-dependent joint source cannot be a schema-3 column."""
    verify_frozen(calibration)
    require(calibration["schema"] == "sccm-reference-calibration-v6" and target_peclet > 0,
            "missing family calibration/common target")
    required = {"coefficient_model", "mean_free_path_authority", "mover", "return_policy", "calibration_horizon", "geometry"}
    require(required <= set(calibration["authorities"]) and all(calibration["authorities"][key] for key in required), "incomplete reference calibration")
    require(members and all(member["source_is_separable"] is False for member in members), "family cannot use separable source measure")
    identities = set()
    for member in members:
        key = (member["species_id"], member["patch_id"], member["time_lower_s"], member["time_upper_s"], member["momentum_lower_si"])
        require(key not in identities, "duplicate family stratum"); identities.add(key)
        require(member["geometry_generation"] > 0 and member["joint_measure"] > 0 and
                abs(member["integrated_peclet"]-target_peclet) <= 1e-8*target_peclet,
                "family generation/joint measure/common target mismatch")
    from .research_products import validate_family
    return validate_family(freeze_record(dict(schema="sep-field-line-family-bundle-v4", calibration=calibration,
                              target_peclet=target_peclet,members=members,legacy_reader_must_reject=True)))


def transfer_campaign(protocol, run, observations, predictions, metric_registration):
    """Frozen transfer evaluation retains failed metrics and compound-event grouping."""
    verify_frozen(protocol); verify_frozen(run); verify_frozen(observations); verify_frozen(metric_registration)
    require(run["protocol_identity"] == protocol["identity"] and
            observations["data_use_role"] == "withheld-validation" and observations["held_out_group"] == protocol["held_out_group"],
            "transfer observation reuse/wrong compound-event holdout")
    require(metric_registration["frozen_before_observations"] is True and metric_registration["protocol_identity"] == protocol["identity"], "post-hoc transfer metrics")
    names = {"onset_time", "peak_time", "log_intensity", "fluence_ratio", "spectral_index", "anisotropy", "uncertainty_coverage"}
    require(set(metric_registration["tolerances"]) == names, "incomplete transfer metrics/tolerances")
    require(run["release_calibration_fingerprint"] != metric_registration["calibration_event_fingerprint"], "reused campaign-D run fingerprint")
    if protocol["transfer_event"] == "2020-05-29":
        require(protocol["held_out_group"] == "2020-05-27/2020-06-02" and observations["pure_radial_same_flux_tube_claim"] is False,
                "compound interval/pure-radial transfer interpretation")
    expected = {row["observer_id"] for row in observations["rows"]}
    require(expected == {row["observer_id"] for row in predictions} and len(expected) == len(predictions), "missing/duplicate transfer observer")
    rows = []
    for observation in observations["rows"]:
        prediction = next(row for row in predictions if row["observer_id"] == observation["observer_id"])
        require(set(observation["metrics"]) == set(prediction["metrics"]) == names, "missing onset/peak/fluence/spectrum/anisotropy/uncertainty metric")
        metrics = {}
        for name in sorted(names):
            a,b=prediction["metrics"][name],observation["metrics"][name]; finite([a,b])
            if name == "log_intensity":
                value = math.log10(a/b) if a > 0 and b > 0 else None
            elif name == "fluence_ratio":
                value = a/b if b > 0 else None
            else: value=a-b
            tolerance=metric_registration["tolerances"][name]
            residual=None if value is None else (value-1 if name=="fluence_ratio" else value)
            metrics[name]=dict(value=value,applicable=value is not None,passed=residual is not None and abs(residual)<=tolerance)
        rows.append(dict(observer_id=observation["observer_id"],metrics=metrics,passed=all(row["passed"] for row in metrics.values())))
    return freeze_record(dict(schema="sccm-transfer-campaign-result-v6",protocol_identity=protocol["identity"],
        run_identity=run["identity"],observations_identity=observations["identity"],rows=rows,
        independent_event_count=1,passed=all(row["passed"] for row in rows),classification=protocol["classification"]))

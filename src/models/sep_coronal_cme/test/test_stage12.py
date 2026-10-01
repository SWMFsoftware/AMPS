#!/usr/bin/env python3
"""Canonical stage-12 gates with independent synthetic references and mutations.

The global release runner selects one of these four classes by its canonical
ID. Standard unittest discovery runs the same checks. No downloaded event data
or simulated SEP output is used as a substitute for independent observations.
"""
from __future__ import annotations
import argparse
import copy
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"tools"))
from preprocessing.core import (PreprocessingError, canonical, digest, covariance,
    publish_bundle, read_bundle, gaussian, validate_asset_roles, verify_frozen, freeze_record)
from preprocessing.inference import (fit_magnetogram, fit_image_ellipse,
    transform_ephemeris, fit_power_law, fit_ellipsoid, fit_plasma_sheet,
    shock_compression, fold_response, fit_wsa_speeds, critical_mach_table, preprocess_job, fit_kinematics, validate_table)
from preprocessing.campaign import (AXES, preregister_campaign, select_candidates,
    radio_likelihood, height_likelihood, load_withheld_after_freeze)
from preprocessing.protocols import (field_line_requests, freeze_transfer_protocol,
    bind_transfer_run, compare_backgrounds, VARIABLES, MATCHED_FIELDS, BOUNDARY_FIELDS,
    TRANSFER_AUTHORITIES)


def diagonal(size, variance=1.0):
    return [[variance if i == j else 0.0 for j in range(size)] for i in range(size)]

def metadata(identifier, role="construction", kind="synthetic", unit="m^-3", source=None):
    return {"asset_id": identifier, "kind": kind, "units": unit, "frame": "synthetic-HCI",
        "epoch_utc": "2000-01-01T00:00:00Z", "processing_version": "synthetic-stage12-v1",
        "data_use_role": role, "source_uri": "synthetic://"+identifier,
        "source_sha256": digest(source if source is not None else {"fixture": identifier}),
        "coordinate_definition": "declared Cartesian/SI synthetic reference", "support": {"time_s": [0, 10]},
        "independence_id": identifier}

def map_source():
    theta = [math.pi*(i+0.5)/5 for i in range(5) for _ in range(8)]
    phi = [2*math.pi*j/8 for _ in range(5) for j in range(8)]
    # Independent low-order analytic spherical harmonics, NOT the fit basis.
    field = [2e-5*math.sqrt(3/(4*math.pi))*math.cos(t)
             -3e-6*math.sqrt(3/(4*math.pi))*math.sin(t)*math.cos(p)
             +5e-6*math.sqrt(5/(16*math.pi))*(3*math.cos(t)**2-1)
             +1e-8/math.sqrt(4*math.pi) for t, p in zip(theta, phi)]
    return {"theta_rad": theta, "phi_rad": phi, "radial_field_T": field, "covariance_T2": diagonal(len(field), 1e-14)}

def density_source(n0=1e6):
    radius = [1.0, 2.0, 3.0, 4.0]
    values = [n0/r**2 for r in radius]
    return {"radius_m": radius, "values": values, "covariance": [[(1e-3*values[i])**2 if i == j else 0.0 for j in range(4)] for i in range(4)]}

def ephemeris_source():
    return {"time_s": [0.0, 1.0], "position_m": [[10, 0, 0], [10, 1, 0]],
        "velocity_m_per_s": [[0, 1, 0], [0, 1, 0]], "covariance": [diagonal(6, 0.01), diagonal(6, 0.01)]}

def ephemeris_parameters():
    return {"rotation": [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            "source_frame": "synthetic-HCI", "target_frame": "synthetic-rotated", "transform_epoch_utc": "2000-01-01T00:00:00Z"}

def pipeline_fixture(base):
    source = density_source()
    (base/"density.json").write_bytes(canonical(source))
    return {"schema": "sep-observation-preprocess-job-v1", "epoch_utc": "2000-01-01T00:00:00Z",
        "output_frame": "synthetic-HCI", "processing_version": "synthetic-stage12-v1",
        "assets": [{"metadata": metadata("density", source=source), "source_file": "density.json",
                    "operator": "power-law-profile", "parameters": {"reference_radius_m": 1.0}}]}


class PROV3D01(unittest.TestCase):
    def test_magnetogram_parameters_and_covariance(self):
        fit = fit_magnetogram(map_source(), {"maximum_degree": 2, "maximum_removed_monopole_T": 1e-7})
        values = {(r["degree"], r["order"], r["sine"]): r["coefficient_T"] for r in fit["modes"]}
        for mode, expected in { (1, 0, False): 2e-5, (1, 1, False): 3e-6, (2, 0, False): 5e-6 }.items():
            self.assertAlmostEqual(values[mode], expected, delta=1e-15)
        self.assertAlmostEqual(fit["removed_monopole_T"], 1e-8, delta=1e-15)
        covariance(fit["covariance_T2"], 8)

    def test_image_ellipse_and_reconstructed_ellipsoid(self):
        angle = 0.3
        points = [[3+4*math.cos(t)*math.cos(angle)-2*math.sin(t)*math.sin(angle),
                   2+4*math.cos(t)*math.sin(angle)+2*math.sin(t)*math.cos(angle)]
                  for t in [2*math.pi*i/32 for i in range(32)]]
        fit = fit_image_ellipse({"contour_xy": points, "implicit_equation_covariance": diagonal(32, 1e-8)},
            {"coordinate_scale_m": 1.0, "projection": "synthetic-known-plane"})
        for actual, expected in zip(fit["parameters"], [3, 2, 4, 2, angle]):
            self.assertAlmostEqual(actual, expected, delta=1e-9)
        covariance(fit["covariance"], 5)
        surface = [[4, 0, 0], [-4, 0, 0], [0, 3, 0], [0, -3, 0], [0, 0, 2], [0, 0, -2]]
        ellipsoid = fit_ellipsoid({"surface_points_m": surface, "implicit_equation_covariance": diagonal(6, 1e-6)},
            {"body_to_frame_rotation": diagonal(3), "center_m": [0, 0, 0], "normalization_m": 1.0,
             "reconstruction_asset_sha256": digest(surface)})
        for actual, expected in zip(ellipsoid["axes_m"], [4, 3, 2]):
            self.assertAlmostEqual(actual, expected, delta=1e-12)
        covariance(ellipsoid["axis_covariance_m2"], 3)

    def test_kinematic_center_axes_joint_inference(self):
        times = [0, 1, 2, 3]
        source = {"time_s": times, "center_m": [[t, 2*t, -t] for t in times],
            "axes_m": [[2+t, 3+2*t, 4+0.5*t] for t in times], "joint_center_axes_covariance_m2": diagonal(24, 0.01)}
        fit = fit_kinematics(source, {"degree": 1, "reference_time_s": 0, "time_scale_s": 1,
            "attitude_asset_sha256": "a"*64})
        for row, expected in zip(fit["center_axes_coefficients_m"], [[0, 1], [0, 2], [0, -1], [2, 1], [3, 2], [4, 0.5]]):
            for value, reference in zip(row, expected): self.assertAlmostEqual(value, reference, delta=1e-12)
        covariance(fit["joint_coefficient_covariance"], 12)

    def test_ephemeris_frame_covariance_and_reviewed_requests(self):
        fit = transform_ephemeris(ephemeris_source(), ephemeris_parameters())
        self.assertEqual(fit["position_m"][0], [0, 10, 0])
        self.assertEqual(fit["velocity_m_per_s"][0], [-1, 0, 0])
        covariance(fit["covariance"][0], 6)
        requests = field_line_requests({"metadata": metadata("ephemeris"), "product": fit},
            [{"stable_line_id": "observer-line", "observer_id": "synthetic-observer", "time_s": 0.0,
              "solar_radius_m": 1, "outer_radius_m": 20, "nominal_step_m": 0.1,
              "maximum_steps_per_branch": 10000, "unsigned_magnetic_flux_wb": 0.1}], "synthetic-reference-review")
        verify_frozen(requests)
        self.assertEqual(requests["requests"][0]["seed_m"], [0, 10, 0])

    def test_plasma_wind_compression_and_response_uncertainty(self):
        profile = fit_power_law(density_source(), {"reference_radius_m": 1})
        self.assertAlmostEqual(profile["exponent"], -2, delta=1e-12)
        self.assertAlmostEqual(math.exp(profile["log_reference_value"]), 1e6, delta=1e-6)
        distance = [-2, -1, 0, 1, 2]
        enhancement = [5*math.exp(-d*d/2) for d in distance]
        sheet = fit_plasma_sheet({"distance_m": distance, "density_m3": [2+x for x in enhancement],
            "covariance_m6": diagonal(5, 1e-4)}, {"background_density_m3": 2, "distance_scale_m": 1})
        self.assertAlmostEqual(sheet["width_m"], 1, delta=1e-12)
        compression = shock_compression({"upstream_density_m3": 2, "downstream_density_m3": 6,
            "joint_covariance_m6": [[0.04, 0.02], [0.02, 0.09]]}, {"definition": "number-density-ratio"})
        self.assertEqual(compression["compression"], 3)
        # Analytic J C J^T, including the covariance cross term.
        self.assertAlmostEqual(compression["variance"], 0.0825, delta=1e-12)
        response = fold_response({"values": [2, 4], "covariance": [[1, 0.2], [0.2, 2]]},
            {"response_matrix": [[0.25, 0.75]], "response_sha256": digest([[0.25, 0.75]]),
             "input_definition": "differential intensity", "output_definition": "channel intensity"})
        self.assertEqual(response["values"], [3.5]); self.assertAlmostEqual(response["covariance"][0][0], 1.2625)
        parameters = {"alpha": 0.5, "beta": 0.5, "width_rad": 1, "distance_power": 2, "outer_exponent": 1}
        f, d = [1, 2, 3, 4], [0.0, 0.5, 1, 2]
        v = [300000+400000*(1+x)**(-0.5)*(1-0.5*math.exp(-y*y)) for x, y in zip(f, d)]
        speeds = fit_wsa_speeds({"expansion_factor": f, "footpoint_distance_rad": d, "speed_m_per_s": v,
            "speed_covariance": diagonal(4, 100)}, parameters)
        self.assertAlmostEqual(speeds["slow_speed_m_per_s"], 300000, delta=1e-6)
        self.assertAlmostEqual(speeds["fast_speed_m_per_s"], 700000, delta=1e-6)
        table = critical_mach_table({"beta": [0, 1], "obliquity_rad": [0, math.pi/2],
            "values": [[1.2, 2.5], [1.1, 2.0]], "value_covariance": diagonal(4, 0.01)},
            {"mach_convention": "fast", "gamma_ad": 5/3, "interpolation": "bounded-bilinear"})
        self.assertEqual(table["schema"], "sep-critical-mach-table-v1")


class PROV3D02(unittest.TestCase):
    def test_byte_determinism_source_and_processing_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary); job = pipeline_fixture(base)
            first = preprocess_job(job, base, base/"a")
            second = preprocess_job(copy.deepcopy(job), base, base/"b")
            self.assertEqual(first, second)
            for member in first["members"]:
                self.assertEqual((base/"a"/member["path"]).read_bytes(), (base/"b"/member["path"]).read_bytes())
            changed = copy.deepcopy(job); changed["processing_version"] = "synthetic-stage12-v2"
            changed["assets"][0]["metadata"]["processing_version"] = changed["processing_version"]
            newer = preprocess_job(changed, base, base/"v2")
            self.assertNotEqual(first["identity"], newer["identity"])
            source = density_source(1.1e6); (base/"density.json").write_bytes(canonical(source))
            job["assets"][0]["metadata"]["source_sha256"] = digest(source)
            altered = preprocess_job(job, base, base/"changed-source")
            self.assertNotEqual(first["identity"], altered["identity"])
            read_bundle(base/"a")
            with self.assertRaises(PreprocessingError): preprocess_job(job, base, base/"a")
            (base/"a"/first["members"][0]["path"]).write_text("{}\n")
            with self.assertRaises(PreprocessingError): read_bundle(base/"a")

    def test_actual_cli_round_trip(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary); job = pipeline_fixture(base)
            (base/"job.json").write_bytes(canonical(job))
            for args in (["preprocess", "--job", str(base/"job.json"), "--output", str(base/"bundle")],
                         ["verify", "--bundle", str(base/"bundle")]):
                result = subprocess.run([sys.executable, str(ROOT/"tools/preprocess_observations.py")]+args,
                    text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                self.assertEqual(result.returncode, 0, result.stdout)
            result = subprocess.run([sys.executable, str(ROOT/"tools/preprocess_observations.py"), "preprocess",
                "--job", str(base/"job.json"), "--output", str(base/"bundle")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 2)


class PROV3D03(unittest.TestCase):
    def test_metadata_epoch_unit_frame_and_covariance_rejections(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary); original = pipeline_fixture(base)
            for key in ("epoch_utc", "frame", "source_sha256", "units", "independence_id", "support"):
                job = copy.deepcopy(original); del job["assets"][0]["metadata"][key]
                with self.assertRaises(PreprocessingError): preprocess_job(job, base, base/key)
            for key, value in [("epoch_utc", "2001-01-01T00:00:00Z"), ("frame", "wrong-frame"),
                               ("units", "unknown"), ("source_sha256", "0"*64)]:
                job = copy.deepcopy(original); job["assets"][0]["metadata"][key] = value
                with self.assertRaises(PreprocessingError): preprocess_job(job, base, base/(key+"-bad"))
        with self.assertRaises(PreprocessingError): covariance([[1, 2], [2, 1]], 2)
        with self.assertRaises(PreprocessingError): covariance([[1, 0], [0, -1]], 2)
        with self.assertRaises(PreprocessingError): transform_ephemeris(ephemeris_source(), dict(ephemeris_parameters(), rotation=[[-1, 0, 0], [0, 1, 0], [0, 0, 1]]))

    def test_role_response_and_validation_reuse(self):
        a = {"metadata": metadata("construction"), "product": {}}
        b = {"metadata": metadata("qualification", "qualification"), "product": {}}
        for field in ("source_sha256", "independence_id", "response_id"):
            reused = copy.deepcopy([a, b]); reused[0]["metadata"][field] = reused[1]["metadata"][field] = "shared" if field != "source_sha256" else "f"*64
            with self.assertRaises(PreprocessingError): validate_asset_roles(reused)
        with self.assertRaises(PreprocessingError): validate_asset_roles([a, copy.deepcopy(a)])

    def test_transfer_and_matched_mhd_manifest_boundaries(self):
        protocol = {"schema": "sep-event-transfer-protocol-v1", "calibration_event": "event-D", "transfer_event": "event-E",
            "held_out_group": "compound-episode-E", "frozen_authorities": {key: "frozen-v1" for key in TRANSFER_AUTHORITIES},
            "transferable_constants": {"source_efficiency": 1e-4}, "event_specific_inputs": ["magnetogram", "front"],
            "freeze_before_withheld_sep": True, "front_reconstruction_independent": True, "classification": "primary-transfer",
            "observer_paths": [{"observer_id": "A", "icme_screening_asset_sha256": "a"*64,
                "background_coverage_complete": True, "sep_coverage_complete": True, "prior_magnetic_cloud": False, "cloud_background_validated": False}]}
        frozen = freeze_transfer_protocol(protocol)
        release = {"transferable_constants": protocol["transferable_constants"], "frozen_authorities": protocol["frozen_authorities"]}
        event_inputs = {key: {"asset_id": key+"-E", "content_sha256": digest(key), "inference_procedure": "frozen-v1"} for key in ("magnetogram", "front")}
        first = bind_transfer_run(frozen, event_inputs, release)
        changed_inputs = copy.deepcopy(event_inputs); changed_inputs["magnetogram"]["content_sha256"] = digest("changed")
        second = bind_transfer_run(frozen, changed_inputs, release)
        self.assertNotEqual(first["release_calibration_fingerprint"], second["release_calibration_fingerprint"])
        with self.assertRaises(PreprocessingError): bind_transfer_run(frozen, {"magnetogram": "a", "front": "b", "mfp": "retuned"}, release)
        stress = copy.deepcopy(protocol); stress["transfer_event"] = "2012-05-17"
        with self.assertRaises(PreprocessingError): freeze_transfer_protocol(stress)
        stress["classification"] = "stress-test"; freeze_transfer_protocol(stress)
        malformed = copy.deepcopy(protocol); malformed["observer_paths"][0]["icme_screening_asset_sha256"] = "named-but-unhashed"
        with self.assertRaises(PreprocessingError): freeze_transfer_protocol(malformed)
        metadata_ = {key: "same" for key in MATCHED_FIELDS|BOUNDARY_FIELDS}
        metadata_.update(epoch_utc="2000-01-01T00:00:00Z", frame="synthetic-HCI", cadence_s=1,
            spatial_support={"minimum_m": [0, -1, -1], "maximum_m": [2, 1, 1], "time_s": [0, 1]},
            masks=[True], interpolation="none-analytic-at-nodes")
        metadata_.update(magnetogram_sha256=digest("magnetogram"), preprocessing_sha256=digest("preprocessing"),
            front_history_sha256=digest("front-history"),
            magnetic_normalization={"value": 1.0, "units": "T", "definition": "synthetic-inner-field-amplitude"},
            open_flux_normalization={"value": 1.0, "units": "Wb", "definition": "synthetic-unsigned-open-flux"})
        metadata_["variable_definitions"] = {key: "synthetic-"+key for key in VARIABLES}
        metadata_["units"] = {"number_density_m3": "m^-3", "mass_density_kg_m3": "kg/m^3", "temperature_K": "K", "pressure_Pa": "Pa",
            "magnetic_field_T": "T", "velocity_m_per_s": "m/s", "alfven_speed_m_per_s": "m/s", "fast_speed_m_per_s": "m/s",
            "signed_open_flux_wb": "Wb", "unsigned_open_flux_wb": "Wb", "D2_first_fast_height_m": "m", "D2_first_supercritical_height_m": "m"}
        variables = {key: {"values": [1, 2, 3] if key in {"magnetic_field_T", "velocity_m_per_s"} else [1],
                          "covariance": diagonal(3 if key in {"magnetic_field_T", "velocity_m_per_s"} else 1, 0.01)} for key in VARIABLES}
        sample = {"schema": "sep-offline-background-samples-v1", "metadata": metadata_, "variables": variables,
            "sample_coordinates": [[1, 0, 0, 0]], "source_uri": "synthetic://model", "source_sha256": "b"*64,
            "model_version": "synthetic-model-v1", "run_id": "synthetic-run", "residual_tuning_used": False,
            "topology": {"definition": "open-closed", "values": ["open"]},
            "connectivity": {"definition": "stable-footpoint", "values": ["A"]}}
        sample = freeze_record(sample)
        comparison = compare_backgrounds(sample, copy.deepcopy(sample))
        self.assertFalse(comparison["truth_validation"])
        self.assertFalse(comparison["runtime_imported_provider_qualified"])
        self.assertTrue(all(value["rms_difference"] == 0 for value in comparison["metrics"].values()))
        changed = copy.deepcopy(sample); del changed["identity"]
        changed["metadata"]["magnetogram_sha256"] = digest("other")
        changed = freeze_record(changed)
        with self.assertRaises(PreprocessingError): compare_backgrounds(sample, changed)
        self.assertEqual(compare_backgrounds(sample, changed, True)["classification"], "combined-boundary-plus-model-discrepancy")
        del changed["identity"]; changed["metadata"]["frame"] = "other"
        changed = freeze_record(changed)
        with self.assertRaises(PreprocessingError): compare_backgrounds(sample, changed, True)
        for key in ("magnetogram_sha256", "preprocessing_sha256", "front_history_sha256", "magnetic_normalization", "open_flux_normalization"):
            malformed = copy.deepcopy(sample); del malformed["identity"]
            malformed["metadata"][key] = "named-but-unqualified"
            with self.assertRaises(PreprocessingError): compare_backgrounds(sample, freeze_record(malformed), True)


def candidate_fixture():
    assets = []
    def add(identifier, role="construction", product=None, kind="background", units="m^-3"):
        asset = {"metadata": metadata(identifier, role, kind, units), "product": product or {"value": identifier}}
        assets.append(asset); return asset
    axes = {}
    for axis in AXES:
        count = 2 if axis in {"magnetogram", "field_scale", "wind_density", "front"} else 1
        members = []
        for i in range(count):
            identifier = axis+str(i)
            product = fit_power_law(density_source(1e6*(1+0.5*i)), {"reference_radius_m": 1.0}) if axis == "wind_density" else None
            add(identifier, product=product)
            members.append({"id": identifier, "asset_id": identifier, "prior_weight": 1.0/count})
        axes[axis] = members
    for name in ("D6", "topology", "D1", "D2"): add("qual"+name, "qualification", kind="diagnostic-reference")
    # Manufactured radio frequency from a known electron density, independently
    # evaluated here using the SI constants, with correlated observing errors.
    epsilon0 = 1/(1.25663706212e-6*299792458.0**2)
    factor = math.sqrt(1.602176634e-19**2/(epsilon0*9.1093837015e-31))/(2*math.pi)
    observation = {"time_s": [1, 2], "frequency_hz": [factor*math.sqrt(1e6/r**2) for r in [2, 3]],
        "frequency_covariance_hz2": [[100, 20], [20, 100]],
        "harmonic_hypotheses": [{"harmonic": 1, "prior_probability": 0.8}, {"harmonic": 2, "prior_probability": 0.2}]}
    add("radio", "qualification", observation, "type-II-radio", "Hz")
    add("withheld", "withheld-validation", {"intensity": [100]}, "SEP-product")
    preregistration = {"schema": "sep-formation-preregistration-v1", "axes": axes,
        "hard_gates": {"D6_maximum_relative_error": 0.1, "topology_minimum_agreement": 0.8, "D1_maximum_relative_error": 0.1},
        "constraint": {"kind": "radio-frequency-time", "asset_id": "radio"}, "weighting_rule": "preregistered-prior-times-joint-likelihood",
        "density_construction_procedure": "synthetic-rebuild-every-tuple-v1", "density_construction_equations_sha256": digest("synthetic-n0-r^-2") }
    realized = [{"tuple": list(key), "D6": 0.01, "topology": 0.95, "D1": 0.02,
        "D2": {"first_fast_height_m": 1.5, "first_fast_time_s": 0.5, "first_supercritical_height_m": 1.8, "first_supercritical_time_s": 0.8},
        "qualification_assets": {name: "qual"+name for name in ("D6", "topology", "D1", "D2")},
        "front_radii_m": [2, 3]} for key in itertools.product(*[[m["id"] for m in axes[axis]] for axis in AXES])]
    attach_density_construction(preregistration, realized, assets)
    return preregistration, realized, assets


def attach_density_construction(preregistration, realized, assets):
    known = {a["metadata"]["asset_id"]: a for a in assets}
    indices = {axis: {m["id"]: m for m in preregistration["axes"][axis]} for axis in AXES}
    for row in realized:
        key = row["tuple"]
        wind = known[indices["wind_density"][key[3]]["asset_id"]]
        row["density_construction"] = {"tuple": list(key), "procedure": preregistration["density_construction_procedure"],
            "equations_sha256": preregistration["density_construction_equations_sha256"],
            "source_checksums": {axis: known[indices[axis][member]["asset_id"]]["metadata"]["source_sha256"] for axis, member in zip(AXES, key)},
            "profile": copy.deepcopy(wind["product"])}


class CAL3D01(unittest.TestCase):
    def test_complete_density_conditioned_radio_product_and_freeze(self):
        preregistration, realized, assets = candidate_fixture()
        frozen = preregister_campaign(preregistration, assets)
        selection = select_candidates(frozen, realized, assets)
        self.assertEqual(len(selection["rows"]), 16)
        self.assertAlmostEqual(sum(r["weight"] for r in selection["rows"]), 1, delta=1e-12)
        branch = selection["rows"][0]["likelihood"]["branches"][0]
        self.assertAlmostEqual(branch["chi2"], 0, delta=1e-16)
        first_density = branch["density_at_front_m3"]
        higher = next(r for r in selection["rows"] if r["tuple"][3] == "wind_density1")
        self.assertAlmostEqual(higher["likelihood"]["branches"][0]["density_at_front_m3"][0]/first_density[0], 1.5, delta=1e-12)
        self.assertGreater(selection["rows"][0]["weight"], higher["weight"])
        self.assertEqual(load_withheld_after_freeze(selection, assets[-1]), {"intensity": [100]})
        changed = copy.deepcopy(selection); changed["rows"][0]["weight"] += 0.01
        with self.assertRaises(PreprocessingError): verify_frozen(changed)
        sampled = copy.deepcopy(realized)
        for row in sampled:
            old = row["density_construction"]["profile"]
            from preprocessing.inference import evaluate_power_law
            n, c = evaluate_power_law(old, row["front_radii_m"])
            row["density_construction"]["profile"] = {"schema": "sep-front-sampled-density-v1", "front_radii_m": row["front_radii_m"],
                "density_m3": n, "covariance_m6": c, "geometry_uncertainty_model": "included-in-covariance"}
        sampled_selection = select_candidates(frozen, sampled, assets)
        for exact, tabulated in zip(selection["rows"], sampled_selection["rows"]):
            self.assertAlmostEqual(exact["weight"], tabulated["weight"], delta=1e-12)
        # Rejecting one tuple retains the full product and exactly zero weight.
        rejected = copy.deepcopy(realized); rejected[0]["D6"] = 0.5
        result = select_candidates(frozen, rejected, assets)
        self.assertEqual(len(result["rows"]), 16); self.assertEqual(result["rows"][0]["weight"], 0)

    def test_cartesian_leakage_and_posthoc_mutations(self):
        preregistration, realized, assets = candidate_fixture(); frozen = preregister_campaign(preregistration, assets)
        for bad in (realized[:-1], realized+[realized[0]]):
            with self.assertRaises(PreprocessingError): select_candidates(frozen, bad, assets)
        changed = copy.deepcopy(frozen); changed["axes"]["field_scale"][0]["prior_weight"] = 0.8
        with self.assertRaises(PreprocessingError): select_candidates(changed, realized, assets)
        for replacement in ("magnetogram0", "withheld"):
            bad = copy.deepcopy(realized); bad[0]["qualification_assets"]["D6"] = replacement
            with self.assertRaises(PreprocessingError): select_candidates(frozen, bad, assets)
        bad = copy.deepcopy(realized); bad[0]["density_construction"]["tuple"][0] = "another-magnetogram"
        with self.assertRaises(PreprocessingError): select_candidates(frozen, bad, assets)
        bad = copy.deepcopy(realized); bad[0]["density_construction"]["equations_sha256"] = "b"*64
        with self.assertRaises(PreprocessingError): select_candidates(frozen, bad, assets)
        bad = copy.deepcopy(realized); bad[0]["SEP_peak"] = 123
        with self.assertRaises(PreprocessingError): select_candidates(frozen, bad, assets)
        for kind in ("SEP-onset", "universal-height-window"):
            bad = copy.deepcopy(preregistration); bad["constraint"]["kind"] = kind
            with self.assertRaises(PreprocessingError): preregister_campaign(bad, assets)
        bad = copy.deepcopy(preregistration); bad["constraint"]["asset_id"] = "withheld"
        with self.assertRaises(PreprocessingError): preregister_campaign(bad, assets)
        bad_assets = copy.deepcopy(assets); bad_assets[-1]["metadata"]["source_sha256"] = bad_assets[0]["metadata"]["source_sha256"]
        with self.assertRaises(PreprocessingError): preregister_campaign(preregistration, bad_assets)

    def test_joint_height_covariance_and_shared_density_identity(self):
        profile = fit_power_law(density_source(), {"reference_radius_m": 1})
        observation = {"height_m": [3], "density_model_sha256": "a"*64,
            "inferred_log_reference_density": profile["log_reference_value"]+1,
            "joint_height_log_density_covariance": [[4, 1], [1, 1]], "independence": "joint-density-inference"}
        score = height_likelihood(observation, profile, [2], "a"*64)
        # Inverse [[4,1],[1,1]] = [[1,-1],[-1,4]]/3; [1,1] has chi2=1.
        self.assertAlmostEqual(score["chi2"], 1, delta=1e-12)
        independent_wrong = gaussian([1, 1], [[4, 0], [0, 1]])["chi2"]
        self.assertNotEqual(score["chi2"], independent_wrong)
        for field, value in (("independence", "independent"), ("density_model_sha256", "b"*64)):
            bad = dict(observation); bad[field] = value
            with self.assertRaises(PreprocessingError): height_likelihood(bad, profile, [2], "a"*64)
        preregistration, realized, assets = candidate_fixture()
        # All related wind members must declare the same inversion source and
        # use its FULL joint covariance; no duplicate independent likelihoods.
        wind = next(a for a in assets if a["metadata"]["asset_id"] == "wind_density0")
        for a in assets:
            if a["metadata"]["asset_id"].startswith("wind_density"):
                a["metadata"]["source_sha256"] = wind["metadata"]["source_sha256"]
            if a["metadata"]["asset_id"] == "radio":
                a["product"] = dict(observation, time_s=[1], density_model_sha256=wind["metadata"]["source_sha256"])
        preregistration["constraint"]["kind"] = "joint-preinferred-height"
        for row in realized: row["front_radii_m"] = [2]
        attach_density_construction(preregistration, realized, assets)
        selection = select_candidates(preregister_campaign(preregistration, assets), realized, assets)
        self.assertEqual(len(selection["rows"]), 16)
        self.assertTrue(all(row["likelihood"]["independence"] == "joint-density-inference" for row in selection["rows"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", choices=("PROV3D01", "PROV3D02", "PROV3D03", "CAL3D01"))
    args = parser.parse_args()
    classes = [globals()[args.test]] if args.test else [PROV3D01, PROV3D02, PROV3D03, CAL3D01]
    suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(case) for case in classes)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if args.test: print("["+args.test+"] "+("PASS" if result.wasSuccessful() else "FAIL"))
    return 0 if result.wasSuccessful() else 1

if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate a deterministic, explicitly synthetic stage-12 workflow example.

This is a maintenance generator for the checked-in examples, not an event data
fetcher. The manufactured references do not qualify any observational profile.
"""
from pathlib import Path
import math
from preprocessing.core import canonical, digest
from preprocessing.inference import fit_power_law
from preprocessing.campaign import AXES

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/"examples/stage12"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT/"sources").mkdir(exist_ok=True)
    assets = []
    sources_by_id = {}
    def add(identifier, source, operator, parameters, units, kind, role="construction"):
        path = "sources/"+identifier+".json"
        (OUTPUT/path).write_bytes(canonical(source))
        sources_by_id[identifier] = source
        assets.append({"source_file": path, "operator": operator, "parameters": parameters,
            "metadata": {"asset_id": identifier, "kind": kind, "units": units, "frame": "synthetic-HCI",
                "epoch_utc": "2000-01-01T00:00:00Z", "processing_version": "synthetic-stage12-v1",
                "data_use_role": role, "source_uri": "synthetic://stage12/"+identifier,
                "source_sha256": digest(source), "coordinate_definition": "declared SI synthetic reference",
                "support": {"time_s": [0, 2]}, "independence_id": "independent-"+identifier}})
    theta = [math.pi*(i+0.5)/5 for i in range(5) for _ in range(8)]
    phi = [2*math.pi*j/8 for _ in range(5) for j in range(8)]
    field = [2e-5*math.sqrt(3/(4*math.pi))*math.cos(t) for t in theta]
    add("magnetogram", {"theta_rad": theta, "phi_rad": phi, "radial_field_T": field,
        "covariance_T2": [[1e-14 if i == j else 0 for j in range(40)] for i in range(40)]},
        "magnetogram-harmonics", {"maximum_degree": 1, "maximum_removed_monopole_T": 1e-7}, "T", "magnetogram")
    table_parameters = {"column_units": {"time_s": "s", "value": "1"}, "positive_columns": ["value"],
        "variable_definitions": {"time_s": "elapsed since synthetic epoch", "value": "manufactured normalized reference"},
        "interpolation": "none", "output_schema": "sep-synthetic-reference-v1", "uncertainty_model": "covariance", "uncertainty_columns": ["value"]}
    for identifier, role in ([("scale", "construction"), ("radii", "construction"), ("front", "construction")]+
                             [("qual"+name, "qualification") for name in ("D6", "topology", "D1", "D2")]+[("withheld", "withheld-validation")]):
        source = {"source_id": identifier, "columns": {"time_s": [0, 1], "value": [1, 1.1]}, "covariance": [[0.01, 0], [0, 0.01]]}
        add(identifier, source, "declared-table", table_parameters, "explicit-column-SI",
            "SEP-product" if identifier == "withheld" else "independent-reference", role)
    for i, normalization in enumerate((1e6, 1.5e6)):
        radius = [1, 2, 3, 4]; values = [normalization/r**2 for r in radius]
        source = {"radius_m": radius, "values": values,
            "covariance": [[(1e-3*values[j])**2 if j == k else 0 for k in range(4)] for j in range(4)]}
        add("density"+str(i), source, "power-law-profile", {"reference_radius_m": 1}, "m^-3", "density")
    epsilon0 = 1/(1.25663706212e-6*299792458.0**2)
    plasma = math.sqrt(1.602176634e-19**2/(epsilon0*9.1093837015e-31))/(2*math.pi)
    add("radio", {"time_s": [1, 2], "frequency_hz": [plasma*math.sqrt(1e6/r**2) for r in (2, 3)],
        "frequency_covariance_hz2": [[100, 20], [20, 100]],
        "harmonic_hypotheses": [{"harmonic": 1, "prior_probability": 0.8}, {"harmonic": 2, "prior_probability": 0.2}]},
        "formation-constraint", {"kind": "radio-frequency-time"}, "explicit-radio-or-height-SI", "type-II-radio", "qualification")
    add("ephemeris", {"time_s": [0, 1], "position_m": [[10, 0, 0], [10, 1, 0]], "velocity_m_per_s": [[0, 1, 0], [0, 1, 0]],
        "covariance": [[[0.01 if i == j else 0 for j in range(6)] for i in range(6)] for _ in range(2)]},
        "ephemeris-transform", {"rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "source_frame": "synthetic-HCI",
            "target_frame": "synthetic-HCI", "transform_epoch_utc": "2000-01-01T00:00:00Z"}, "SI-position-velocity", "observer-ephemeris")
    job = {"schema": "sep-observation-preprocess-job-v1", "epoch_utc": "2000-01-01T00:00:00Z", "output_frame": "synthetic-HCI",
        "processing_version": "synthetic-stage12-v1", "assets": assets}
    (OUTPUT/"job.json").write_bytes(canonical(job))
    ids = {"magnetogram": ["magnetogram"], "field_scale": ["scale"], "coupling_radii": ["radii"],
           "wind_density": ["density0", "density1"], "front": ["front"]}
    axes = {axis: [{"id": name, "asset_id": name, "prior_weight": 1/len(names)} for name in names] for axis, names in ids.items()}
    preregistration = {"schema": "sep-formation-preregistration-v1", "axes": axes,
        "hard_gates": {"D6_maximum_relative_error": 0.1, "topology_minimum_agreement": 0.8, "D1_maximum_relative_error": 0.1},
        "constraint": {"kind": "radio-frequency-time", "asset_id": "radio"}, "weighting_rule": "preregistered-prior-times-joint-likelihood",
        "density_construction_procedure": "synthetic-rebuild-every-tuple-v1", "density_construction_equations_sha256": digest("synthetic-n0-r^-2")}
    (OUTPUT/"preregistration.json").write_bytes(canonical(preregistration))
    realized = []
    for wind in ids["wind_density"]:
        key = ["magnetogram", "scale", "radii", wind, "front"]
        realized.append({"tuple": key, "D6": 0.01, "topology": 0.95, "D1": 0.02,
            "D2": {"first_fast_height_m": 1.5, "first_fast_time_s": 0.5, "first_supercritical_height_m": 1.8, "first_supercritical_time_s": 0.8},
            "qualification_assets": {name: "qual"+name for name in ("D6", "topology", "D1", "D2")}, "front_radii_m": [2, 3],
            "density_construction": {"tuple": key, "procedure": preregistration["density_construction_procedure"],
                "equations_sha256": preregistration["density_construction_equations_sha256"],
                "source_checksums": {axis: digest(sources_by_id[member]) for axis, member in zip(AXES, key)},
                "profile": fit_power_law(sources_by_id[wind], {"reference_radius_m": 1})}})
    (OUTPUT/"realized-candidates.json").write_bytes(canonical(realized))
    (OUTPUT/"field-line-requests.json").write_bytes(canonical([{"stable_line_id": "synthetic-observer-line", "observer_id": "synthetic-observer",
        "time_s": 0, "solar_radius_m": 1, "outer_radius_m": 20, "nominal_step_m": 0.1,
        "maximum_steps_per_branch": 10000, "unsigned_magnetic_flux_wb": 0.1}]))
    print("generated synthetic stage12 inputs (not an observational qualification)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

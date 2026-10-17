#!/usr/bin/env python3
"""Stage-5 absorber-entry DDD/n_eq transport and strict-domain acceptance.

The exact-node cases traverse the normal case parser, primary generator,
Geant4 geometry, upstream-face scorer, NIEL lookup, per-event statistics, and
result writer. A negligible but nonzero aluminum layer is retained so these are
full shield transport cases while their independently calculated reference is
still the incident exact table node to the registered 0.5% tolerance.
"""

from __future__ import annotations

import copy
import json
import math
import pathlib
import subprocess
import sys


if len(sys.argv) != 5:
    raise SystemExit(
        "usage: test_stage5_transport.py SHIELDSIM SOURCE_ROOT EVIDENCE_DIR RUN_DIR"
    )

EXECUTABLE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
PLAN = json.loads(
    (ROOT / "tests/reference_ddd/stage5_validation_plan.v1.json").read_text(encoding="utf-8")
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def set_measure(container: dict, name: str, value: float, unit: str = "mm") -> None:
    container[name] = {"value": value, "unit": unit}


def configure_sampling(case: dict, histories: int, seed: int) -> None:
    case["accuracy"]["sampling"] = {
        "min_histories": histories,
        "max_histories": histories,
        "batch_size": histories,
    }
    case["transport"]["random_seed"] = seed


def run_case(name: str, case: dict) -> tuple[dict, dict]:
    directory = RUN_DIR / name
    directory.mkdir(parents=True, exist_ok=True)
    document = copy.deepcopy(case)
    document["case_id"] = f"stage5_{name}"
    document["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    document["damage_response"]["path"] = str(
        ROOT / "data/niel/damage_response_library.v1.json"
    )
    document["output"]["prefix"] = "stage5"
    case_path = directory / "case.json"
    case_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(EXECUTABLE), f"--case={case_path}"],
        cwd=directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (directory / "run.log").write_text(completed.stdout, encoding="utf-8")
    require(completed.returncode == 0, f"transport case {name} failed; see {directory/'run.log'}")
    absorber = document["geometry"]["absorbers"][0]["absorber_id"]
    result_path = directory / f"stage5__{absorber}_result.json"
    response_path = directory / f"stage5__{absorber}_response_cell.json"
    require(result_path.is_file() and response_path.is_file(),
            f"transport case {name} omitted result/response output")
    return (json.loads(result_path.read_text(encoding="utf-8")),
            json.loads(response_path.read_text(encoding="utf-8")))


def output(result: dict, name: str) -> dict:
    matches = [item for item in result["outputs"] if item["name"] == name]
    require(len(matches) == 1, f"expected exactly one {name} output")
    return matches[0]


def scores(payload: dict) -> dict[str, dict]:
    key = "raw_transport_scores" if "raw_transport_scores" in payload else "scores"
    return {item["score_id"]: item for item in payload[key]}


def main() -> int:
    require(EXECUTABLE.is_file(), f"missing shieldSim executable: {EXECUTABLE}")
    require(PLAN["registration_status"].startswith("fixed_before_automated"),
            "Stage-5 acceptance plan registration is missing")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    base = json.loads((ROOT / "examples/stage5_ddd_neq_case.json").read_text(encoding="utf-8"))
    area_cm2 = (
        base["geometry"]["absorbers"][0]["geometry"]["size_x"]["value"] / 10.0
    ) * (
        base["geometry"]["absorbers"][0]["geometry"]["size_y"]["value"] / 10.0
    )
    reference_niel = 95.0e-27 * 6.02214076e23 / 28.0855
    comparisons: list[dict] = []
    matrix: list[dict] = []

    for index, registered in enumerate(PLAN["transport_cases"]):
        case = copy.deepcopy(base)
        case["source"]["spectrum"]["incident_id"] = registered["id"]
        case["source"]["spectrum"]["particle"] = registered["particle"]
        case["source"]["spectrum"]["energy"] = {
            "coordinate": "total_kinetic_energy",
            "value": registered["energy_MeV"],
            "unit": "MeV",
        }
        case["geometry"]["shield_stack"]["layers"][0]["thickness"]["value"] = \
            PLAN["negligible_shield_thickness_mm"]
        case["quantities"] = ["DDD", "n_eq"] if registered["particle"]["name"] == "neutron" else ["DDD"]
        case["transport"]["physics_list"] = (
            "FTFP_BERT_HP" if registered["particle"]["name"] == "neutron" else "FTFP_BERT"
        )
        configure_sampling(case, PLAN["transport_histories"], 520000 + index)
        result, response = run_case(registered["id"], case)
        ddd = output(result, "DDD")
        expected_ddd = registered["niel_MeV_cm2_g"] / area_cm2
        require(ddd["value"] is not None and ddd["table_domain_status"] == "complete",
                f"{registered['id']} DDD unexpectedly unavailable")
        relative = abs(float(ddd["value"]) - expected_ddd) / expected_ddd
        require(relative <= PLAN["transport_reference_relative_tolerance"],
                f"{registered['id']} entry-fold DDD differs from the exact source node by {relative:.4g}")
        raw = scores(result)
        response_raw = scores(response)
        particle = registered["particle"]["name"]
        prefix = f"absorber/{result['absorber_id']}/entry/by_particle/{particle}"
        count = raw[f"{prefix}/crossing_count"]["mean_per_primary"]
        require(count >= PLAN["entry_crossing_fraction_min"],
                f"{registered['id']} has too few upstream-face crossings: {count}")
        raw_ddd = raw[f"{prefix}/DDD_MeV_g"]["mean_per_primary"]
        require(math.isclose(raw_ddd, float(ddd["value"]), rel_tol=1.0e-12),
                f"{registered['id']} result does not close to event-level entry score")
        require(math.isclose(response_raw[f"{prefix}/DDD_MeV_g"]["mean_per_primary"],
                             raw_ddd, rel_tol=1.0e-12),
                f"{registered['id']} response-cell primitive differs from result primitive")
        damage_key = response["response_key_document"].get("damage_response", {})
        require(damage_key.get("manifest_hash") == result["damage_response"]["manifest_hash"] and
                response["response_key_document"]["scoring_contract"] ==
                "shieldSim.transport_raw_plus_silicon_damage.v1",
                f"{registered['id']} response cache key omits Stage-5 damage provenance")
        comparisons.append({
            "id": f"{registered['id']}_absorber_entry_DDD",
            "unit": "MeV/g/primary",
            "expected": expected_ddd,
            "actual": float(ddd["value"]),
            "relative_difference": relative,
            "acceptance_relative": PLAN["transport_reference_relative_tolerance"],
        })
        if particle == "neutron":
            neq = output(result, "n_eq")
            expected_neq = expected_ddd / reference_niel
            require(neq["value"] is not None and
                    math.isclose(float(neq["value"]), expected_neq,
                                 rel_tol=PLAN["transport_reference_relative_tolerance"]),
                    "neutron n_eq does not use the prescribed E722 normalization")
            require("proton and He-4" in neq["excluded_from_definition"],
                    "n_eq output does not state its neutron-only boundary")
            comparisons.append({
                "id": "neutron_1p05MeV_ASTM_n_eq",
                "unit": "cm-2/primary",
                "expected": expected_neq,
                "actual": float(neq["value"]),
                "relative_difference": abs(float(neq["value"]) - expected_neq) / expected_neq,
                "acceptance_relative": PLAN["transport_reference_relative_tolerance"],
            })
        diagnostic = float(ddd["diagnostic_Geant4_nonionizing_deposition"]["value"])
        require(ddd["authoritative_method"].startswith("absorber-entry scalar fluence") and
                "never added" in ddd["diagnostic_Geant4_nonionizing_deposition"]["definition"],
                "table fold and Geant4 nonionizing-deposition diagnostic are not separated")
        matrix.append({
            "id": registered["id"], "DDD_MeV_g_per_primary": ddd["value"],
            "diagnostic_nonionizing_MeV_g_per_primary": diagnostic,
            "entry_crossings_per_primary": count,
        })
        require("pending_accelerator_device_validation_and_independent_review" in
                ddd["validation_status"], "Stage-5 completion boundary is missing")

    # Fluence is defined per area normal to the particle direction. A crossing
    # of a fixed slab face at cosine mu therefore represents 1/(A*mu), not 1/A.
    # Keep the source close to the shield and compact enough that every oblique
    # primary intersects both the finite shield and absorber; this isolates the
    # projected-area identity from lateral acceptance losses.
    angular = PLAN["oblique_fluence_case"]
    oblique = copy.deepcopy(base)
    oblique["source"]["spectrum"]["incident_id"] = "proton_100MeV_mu0p5"
    oblique["source"]["spectrum"]["particle"] = angular["particle"]
    oblique["source"]["spectrum"]["energy"]["value"] = angular["energy_MeV"]
    oblique["source"]["spectrum"]["direction"] = {
        "mu": angular["mu"], "phi": {"value": angular["phi_deg"], "unit": "deg"}
    }
    oblique["quantities"] = ["DDD"]
    oblique["transport"]["physics_list"] = "FTFP_BERT"
    oblique["geometry"]["source_plane"]["size_x"]["value"] = 10.0
    oblique["geometry"]["source_plane"]["size_y"]["value"] = 10.0
    oblique["geometry"]["source_plane"]["transform"]["translation"]["z"] = {
        "value": -22.0, "unit": "mm"
    }
    oblique["geometry"]["absorbers"][0]["geometry"]["transform"]["translation"]["z"] = {
        "value": -19.9, "unit": "mm"
    }
    oblique["geometry"]["shield_stack"]["layers"][0]["thickness"]["value"] = \
        PLAN["negligible_shield_thickness_mm"]
    configure_sampling(oblique, angular["histories"], 520050)
    oblique_result, _ = run_case("proton_100MeV_mu0p5", oblique)
    oblique_ddd = output(oblique_result, "DDD")
    expected_oblique = angular["niel_MeV_cm2_g"] / (area_cm2 * angular["mu"])
    oblique_relative = abs(float(oblique_ddd["value"]) - expected_oblique) / expected_oblique
    require(oblique_relative <= angular["acceptance_relative"],
            f"oblique scalar-fluence DDD omitted the 1/mu weight: {oblique_relative:.4g}")
    oblique_raw = scores(oblique_result)
    oblique_prefix = f"absorber/{oblique_result['absorber_id']}/entry/by_particle/proton"
    expected_fluence = 1.0 / (area_cm2 * angular["mu"])
    require(math.isclose(oblique_raw[f"{oblique_prefix}/scalar_fluence_cm2"]["mean_per_primary"],
                         expected_fluence, rel_tol=angular["acceptance_relative"]),
            "oblique raw scalar fluence does not equal crossings/(A*mu)")
    comparisons.append({
        "id": "proton_100MeV_oblique_projected_fluence_DDD",
        "unit": "MeV/g/primary",
        "expected": expected_oblique,
        "actual": float(oblique_ddd["value"]),
        "relative_difference": oblique_relative,
        "acceptance_relative": angular["acceptance_relative"],
    })
    matrix.append({
        "id": "proton_100MeV_mu0p5",
        "DDD_MeV_g_per_primary": oblique_ddd["value"],
        "scalar_fluence_cm2_per_primary":
            oblique_raw[f"{oblique_prefix}/scalar_fluence_cm2"]["mean_per_primary"],
    })

    # A transported primary above 20 MeV must reach the scoring surface and make
    # both quantities explicitly not assessable. It may never be omitted,
    # evaluated at 20 MeV, or extrapolated.
    registered = PLAN["out_of_domain_case"]
    outside = copy.deepcopy(base)
    outside["source"]["spectrum"]["incident_id"] = "neutron_25MeV_out_of_domain"
    outside["source"]["spectrum"]["particle"] = registered["particle"]
    outside["source"]["spectrum"]["energy"]["value"] = registered["energy_MeV"]
    outside["transport"]["physics_list"] = "FTFP_BERT"
    configure_sampling(outside, registered["histories"], 520100)
    outside_result, _ = run_case("neutron_25MeV_out_of_domain", outside)
    for quantity in ("DDD", "n_eq"):
        unavailable = output(outside_result, quantity)
        require(unavailable["value"] is None and
                unavailable["table_domain_status"] == "not_assessable" and
                unavailable["accuracy"]["status"] == "not_assessable",
                f"out-of-domain {quantity} was not explicitly unavailable")
    require(output(outside_result, "DDD")["out_of_domain_crossings"] > 0 and
            output(outside_result, "n_eq")["out_of_domain_neutron_crossings"] > 0,
            "out-of-domain transported neutron was not counted")

    # Secondary-neutron sensitivity: a thick 10-cm Al layer under a 1-GeV
    # proton beam must exercise secondary neutron scoring. Neutrons above the
    # 20-MeV E722 boundary are expected to make n_eq not assessable; lower-energy
    # contributions remain visible as explicitly labeled partial raw scores.
    sensitivity = PLAN["secondary_sensitivity"]
    secondary = copy.deepcopy(base)
    secondary["source"]["spectrum"]["incident_id"] = "proton_1GeV_secondary_neutrons"
    secondary["source"]["spectrum"]["particle"] = sensitivity["particle"]
    secondary["source"]["spectrum"]["energy"]["value"] = sensitivity["energy_MeV"]
    secondary["quantities"] = ["n_eq"]
    secondary["transport"]["physics_list"] = "FTFP_BERT"
    secondary["geometry"]["shield_stack"]["layers"][0]["thickness"]["value"] = \
        sensitivity["aluminum_thickness_mm"]
    secondary["geometry"]["absorbers"][0]["geometry"]["transform"]["translation"]["z"] = {
        "value": 85.0, "unit": "mm"
    }
    configure_sampling(secondary, sensitivity["histories"], 520200)
    secondary_result, _ = run_case("secondary_neutron_sensitivity", secondary)
    neutron_scores = [item for item in secondary_result["raw_transport_scores"]
                      if "/entry/by_particle/neutron/crossing_count" in item["score_id"]]
    require(neutron_scores and neutron_scores[0]["sum"] > 0,
            "thick-shield sensitivity did not produce absorber-entry secondary neutrons")
    secondary_neq = output(secondary_result, "n_eq")
    matrix.append({
        "id": "secondary_neutron_sensitivity",
        "entry_neutron_crossings": neutron_scores[0]["sum"],
        "n_eq_domain_status": secondary_neq["table_domain_status"],
        "out_of_domain_neutron_crossings": secondary_neq["out_of_domain_neutron_crossings"],
    })

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_case": "stage5_transport",
        "status": "passed",
        "scope": "full Geant4 absorber-entry DDD/n_eq path; formal external device validation remains pending",
        "comparisons": comparisons,
        "transport_matrix": matrix,
        "out_of_domain_policy_verified": True,
        "secondary_neutron_sensitivity": matrix[-1],
    }
    (EVIDENCE / "stage5_transport.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Stage-5 transport checks passed: {len(comparisons)} numeric comparisons")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

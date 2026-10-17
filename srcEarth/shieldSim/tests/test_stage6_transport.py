#!/usr/bin/env python3
"""Stage-6 LET-step and weighted-H100/10 end-to-end acceptance tests."""

from __future__ import annotations

import copy
import json
import math
import pathlib
import subprocess
import sys

if len(sys.argv) != 6:
    raise SystemExit(
        "usage: test_stage6_transport.py SHIELDSIM REFERENCE_TOOL SOURCE_ROOT EVIDENCE_DIR RUN_DIR"
    )

SHIELDSIM = pathlib.Path(sys.argv[1]).resolve()
REFERENCE = pathlib.Path(sys.argv[2]).resolve()
ROOT = pathlib.Path(sys.argv[3]).resolve()
EVIDENCE = pathlib.Path(sys.argv[4]).resolve()
RUN_DIR = pathlib.Path(sys.argv[5]).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def stage6_base() -> dict:
    case = json.loads((ROOT / "examples/stage4_tid_case.json").read_text(encoding="utf-8"))
    case["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    case["geometry"]["shield_stack"]["layers"] = [
        {"layer_id": "negligible_al", "material": "Al",
         "thickness": {"type": "geometric", "value": 1.0e-6, "unit": "mm"}}
    ]
    case["geometry"]["absorbers"][0]["geometry"]["thickness"] = {
        "value": 0.01, "unit": "mm"
    }
    case["quantities"] = ["LET"]
    case["let_response"] = {
        "definition_id": "shieldSim.electronic_let.v1",
        "linear_let_bins": {
            "scale": "logarithmic",
            "minimum": {"value": 0.001, "unit": "keV/um"},
            "maximum": {"value": 10000.0, "unit": "keV/um"},
            "bins_per_decade": 20,
        },
        "out_of_range_policy": "not_assessable_no_clamping",
    }
    case["accuracy"] = {
        "target": 0.10, "scope": "statistical",
        "coverage": {"standard_k": 1.0, "standard_probability": 0.68,
                     "interval_probability": 0.95,
                     "interval_method": "fieller_event_covariance"},
        "systematic_budget": {"status": "incomplete"},
        "quantity_overrides": {},
        "sampling": {"min_histories": 100, "max_histories": 100,
                     "batch_size": 100},
    }
    case["transport"]["production_cut"] = {"value": 0.01, "unit": "mm"}
    case["transport"]["max_step"] = {"value": 0.001, "unit": "mm"}
    return case


def fixed_histories(case: dict, histories: int, seed: int) -> None:
    case["accuracy"]["sampling"] = {
        "min_histories": histories, "max_histories": histories, "batch_size": histories
    }
    case["transport"]["random_seed"] = seed


def run_case(name: str, case: dict) -> dict:
    directory = RUN_DIR / name
    directory.mkdir(parents=True, exist_ok=True)
    document = copy.deepcopy(case)
    document["case_id"] = f"stage6_{name}"
    document["output"]["prefix"] = "stage6"
    path = directory / "case.json"
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    completed = subprocess.run([str(SHIELDSIM), f"--case={path}"], cwd=directory,
                               text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               check=False)
    (directory / "run.log").write_text(completed.stdout, encoding="utf-8")
    require(completed.returncode == 0, f"{name} failed; see {directory/'run.log'}")
    absorber = document["geometry"]["absorbers"][0]["absorber_id"]
    result = directory / f"stage6__{absorber}_result.json"
    require(result.is_file(), f"{name} omitted result JSON")
    return json.loads(result.read_text(encoding="utf-8"))


def output(result: dict, name: str) -> dict:
    matches = [item for item in result["outputs"] if item["name"] == name]
    require(len(matches) == 1, f"expected one {name} output")
    return matches[0]


def reference_values(cases: list[tuple[str, str, int, int, float]]) -> dict[str, float]:
    grouped: dict[str, list[dict]] = {
        "reference_list_default": [], "G4EmStandardPhysics_option4": []
    }
    for identifier, particle, z, a, energy in cases:
        query = {"id": identifier, "particle": particle, "material": "G4_Si",
                 "energy_MeV": energy, "quantities": ["electronic_stopping_power"]}
        if particle == "ion":
            query.update({"Z": z, "A": a})
        grouped["G4EmStandardPhysics_option4" if particle == "ion"
                else "reference_list_default"].append(query)
    values: dict[str, float] = {}
    for index, (em_physics, queries) in enumerate(grouped.items()):
        request = {"schema_version": "shieldSim.stage6_reference_request.v1",
                   "physics_list": "FTFP_BERT", "em_physics": em_physics,
                   "ion_physics": "reference_list_default", "queries": queries}
        path = RUN_DIR / f"transport_reference_{index}_request.json"
        result = RUN_DIR / f"transport_reference_{index}_result.json"
        path.write_text(json.dumps(request, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        completed = subprocess.run([str(REFERENCE), str(path), str(result)], cwd=RUN_DIR,
                                   text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   check=False)
        (RUN_DIR / f"transport_reference_{index}.log").write_text(completed.stdout,
                                                                   encoding="utf-8")
        require(completed.returncode == 0, "Stage-6 transport reference probe failed")
        payload = json.loads(result.read_text(encoding="utf-8"))
        values.update({row["id"]: float(row["values"]["electronic_stopping_power"]["value"])
                       for row in payload["results"]})
    return values


def main() -> int:
    require(SHIELDSIM.is_file() and REFERENCE.is_file(), "missing Stage-6 executables")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    comparisons: list[dict] = []
    references = reference_values([
        ("proton", "proton", 1, 1, 100.0),
        ("alpha", "alpha", 2, 4, 100.0),
        ("c12", "ion", 6, 12, 1200.0),
    ])

    species = [
        ("proton", {"name": "proton", "Z": 1, "A": 1},
         {"coordinate": "total_kinetic_energy", "value": 100.0, "unit": "MeV"}),
        ("alpha", {"name": "alpha", "Z": 2, "A": 4},
         {"coordinate": "total_kinetic_energy", "value": 100.0, "unit": "MeV"}),
        ("c12", {"name": "ion", "Z": 6, "A": 12},
         {"coordinate": "kinetic_energy_per_nucleon", "value": 100.0, "unit": "MeV/n"}),
    ]
    baseline_results: dict[str, float] = {}
    aggregate_baseline: dict[str, float] = {}
    for index, (identifier, particle, energy) in enumerate(species):
        case = stage6_base()
        case["source"]["spectrum"]["incident_id"] = f"{identifier}_thin_si"
        case["source"]["spectrum"]["particle"] = particle
        case["source"]["spectrum"]["energy"] = energy
        if identifier == "c12":
            case["transport"]["em_physics"] = "G4EmStandardPhysics_option4"
        fixed_histories(case, 80, 610000 + index)
        result = run_case(f"thin_{identifier}", case)
        let = output(result, "LET")
        track = let["track_mean"]
        primary_components = [component for component in let["components"]
                              if component["primary"] and component["generation"] == 0 and
                              component["Z"] == particle["Z"] and component["A"] == particle["A"]]
        require(len(primary_components) == 1,
                f"{identifier}: unique primary LET component missing")
        primary_track = primary_components[0]["track_mean"]
        require(track["status"] == "available" and primary_track["status"] == "available" and
                let["histogram_domain_status"] == "complete",
                f"{identifier}: LET result unavailable")
        # The reference describes the incident primary stopping. The aggregate
        # spectrum deliberately also includes transported delta electrons and
        # therefore is not substituted for the primary component comparison.
        actual = float(primary_track["value_MeV_cm2_per_g"])
        expected = references[identifier]
        relative = abs(actual - expected) / expected
        require(relative <= 0.03,
                f"{identifier}: thin-slab track mean differs from midpoint reference by {relative:.4g}")
        baseline_results[identifier] = actual
        aggregate_baseline[identifier] = float(track["value_MeV_cm2_per_g"])
        bins = sum(float(row["fluence_cm2_per_primary"]) for row in let["bins"])
        denominator = float(track["event_moments"]["denominator_mean"])
        require(math.isclose(bins + float(let["underflow_fluence_cm2_per_primary"]) +
                             float(let["overflow_fluence_cm2_per_primary"]),
                             denominator, rel_tol=2.0e-12, abs_tol=1.0e-15),
                f"{identifier}: LET histogram does not conserve track fluence")
        require("pending_independent_review" in let["validation_status"],
                "Stage-6 completion boundary missing")
        comparisons.append({
            "id": f"thin_{identifier}_track_mean_vs_midpoint_table",
            "unit": "MeV cm2/g", "expected": expected, "actual": actual,
            "relative_difference": relative, "acceptance_relative": 0.03,
        })

    # Max-step, production-cut, and finite-geometry perturbations are separate
    # convergence axes. They are compared with the registered proton baseline,
    # not pooled into a single undocumented numerical allowance.
    variants = {
        "max_step_half": ("max_step", 0.0005),
        "production_cut_half": ("production_cut", 0.005),
        "absorber_thickness_half": ("absorber_thickness", 0.005),
    }
    for index, (name, (kind, value)) in enumerate(variants.items()):
        case = stage6_base()
        if kind in ("max_step", "production_cut"):
            case["transport"][kind] = {"value": value, "unit": "mm"}
        else:
            case["geometry"]["absorbers"][0]["geometry"]["thickness"] = {
                "value": value, "unit": "mm"
            }
        fixed_histories(case, 80, 610100 + index)
        actual = float(output(run_case(name, case), "LET")["track_mean"]
                       ["value_MeV_cm2_per_g"])
        expected = aggregate_baseline["proton"]
        relative = abs(actual - expected) / expected
        require(relative <= 0.10, f"{name}: LET convergence residual {relative:.4g}")
        comparisons.append({"id": f"LET_convergence_{name}", "unit": "MeV cm2/g",
                            "expected": expected, "actual": actual,
                            "relative_difference": relative, "acceptance_relative": 0.10})

    # The exact two-node mixture has P(E>100 | E>10)=2/3. Surface weights are
    # equal at mu=1, so the transport ratio must recover that analytic result.
    hardness_case = stage6_base()
    hardness_case["quantities"] = ["H100/10"]
    hardness_case.pop("let_response")
    hardness_case["source"]["spectrum"] = {
        "format": "discrete_response_nodes", "nodes": [
            {"incident_id": "p50", "particle": {"name": "proton", "Z": 1, "A": 1},
             "energy": {"coordinate": "total_kinetic_energy", "value": 50.0, "unit": "MeV"},
             "direction": {"mu": 1.0, "phi": {"value": 0.0, "unit": "rad"}},
             "weight": 1.0 / 3.0},
            {"incident_id": "p150", "particle": {"name": "proton", "Z": 1, "A": 1},
             "energy": {"coordinate": "total_kinetic_energy", "value": 150.0, "unit": "MeV"},
             "direction": {"mu": 1.0, "phi": {"value": 0.0, "unit": "rad"}},
             "weight": 2.0 / 3.0},
        ]}
    fixed_histories(hardness_case, 3000, 610200)
    hardness = output(run_case("hardness_two_node", hardness_case), "H100/10")
    require(hardness["ratio_status"] == "available" and hardness["value"] is not None,
            "two-node hardness ratio unavailable")
    actual_h = float(hardness["value"])
    require(abs(actual_h - 2.0 / 3.0) <= 0.04,
            f"weighted H100/10 differs from analytic 2/3: {actual_h:g}")
    require(hardness["interval_95"]["low"] <= 2.0 / 3.0 <= hardness["interval_95"]["high"],
            "Fieller interval does not cover analytic two-node hardness")
    comparisons.append({"id": "weighted_H100_10_two_node_mixture", "unit": "1",
                        "expected": 2.0 / 3.0, "actual": actual_h,
                        "absolute_acceptance": 0.04})

    zero_case = stage6_base()
    zero_case["quantities"] = ["H100/10"]
    zero_case.pop("let_response")
    zero_case["source"]["spectrum"]["particle"] = {"name": "alpha", "Z": 2, "A": 4}
    fixed_histories(zero_case, 20, 610201)
    zero = output(run_case("hardness_zero_denominator", zero_case), "H100/10")
    require(zero["value"] is None and zero["ratio_status"] == "zero_denominator" and
            zero["accuracy"]["status"] == "not_assessable",
            "zero-denominator hardness was not explicitly undefined")

    evidence = {"schema_version": "shieldSim.validation_evidence.v1",
                "test_case": "stage6_transport", "comparisons": comparisons}
    (EVIDENCE / "stage6_transport.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Stage-6 transport checks passed: {len(comparisons)} comparisons")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

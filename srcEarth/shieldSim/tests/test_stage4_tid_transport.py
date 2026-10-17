#!/usr/bin/env python3
"""Stage-4 finite-slab TID transport, closure, and convergence gate.

All generated cases and Geant4 output stay in the caller-provided build-tree
directory. The material/species matrix exercises full transport/scoring across
the validated slabs. Quantity-specific cut, step, finite-size, history, and EM-
model sensitivity are tested without assigning an unreferenced absolute dose.
PSTAR/ASTAR are compared with the deterministic unrestricted-stopping probe;
they are not misused as finite-volume deposited-energy expected values.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import pathlib
import statistics
import subprocess
import sys


if len(sys.argv) != 5:
    raise SystemExit(
        "usage: test_stage4_tid_transport.py SHIELDSIM SOURCE_ROOT EVIDENCE_DIR RUN_DIR"
    )

EXECUTABLE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
PLAN_PATH = ROOT / "tests" / "reference_tid" / "stage4_transport_plan.v1.json"
NIST = ROOT / "tests" / "reference_tid" / "nist"
MEV_J = 1.602176634e-13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def exact_nist_row(table: str, material: str, energy: float) -> dict[str, str]:
    with (NIST / f"{table}_{material}.csv").open(newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream)
                if math.isclose(float(row["energy_MeV"]), energy,
                                rel_tol=0.0, abs_tol=1.0e-12)]
    require(len(rows) == 1, f"ambiguous NIST row: {table}/{material}/{energy}")
    return rows[0]


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
    case = copy.deepcopy(case)
    case["case_id"] = f"stage4_{name}"
    case["material_library"]["path"] = str(ROOT / "data" / "material_library.v1.json")
    case["output"]["prefix"] = "stage4"
    case_path = directory / "case.json"
    case_path.write_text(json.dumps(case, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(EXECUTABLE), f"--case={case_path}"],
        cwd=directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    (directory / "run.log").write_text(completed.stdout, encoding="utf-8")
    require(completed.returncode == 0, f"transport case {name} failed; see {directory/'run.log'}")
    absorber_id = case["geometry"]["absorbers"][0]["absorber_id"]
    result_path = directory / f"stage4__{absorber_id}_result.json"
    resolved_path = directory / f"stage4__{absorber_id}_resolved_case.json"
    require(result_path.is_file() and resolved_path.is_file(), f"case {name} omitted canonical outputs")
    return (json.loads(result_path.read_text(encoding="utf-8")),
            json.loads(resolved_path.read_text(encoding="utf-8")))


def tid_output(result: dict) -> dict:
    matches = [entry for entry in result["outputs"] if entry["name"] == "TID"]
    require(len(matches) == 1, "result must have exactly one TID output")
    return matches[0]


def scores(result: dict) -> dict[str, dict]:
    return {entry["score_id"]: entry for entry in result["raw_transport_scores"]}


def assert_exact_dose_and_mass(name: str, result: dict, resolved: dict,
                               comparisons: list[dict]) -> None:
    output = tid_output(result)
    absorber_id = result["absorber_id"]
    raw = scores(result)
    ion = raw[f"absorber/{absorber_id}/ionizing_deposited_energy_MeV"]
    mean_ion_mev = ion["sum"] / result["events"]
    mass_kg = output["absorber_mass"]["value"]
    expected_gy = mean_ion_mev * MEV_J / mass_kg
    require(math.isclose(output["value"], expected_gy, rel_tol=2.0e-12),
            f"{name}: TID is not ionizing event-ledger energy divided by mass")
    require(math.isclose(output["rad_per_primary"], 100.0 * output["value"], rel_tol=2.0e-12),
            f"{name}: rad/Gy conversion is not exact")
    geometry = resolved["geometry"]["absorbers"][0]
    volume_cm3 = (geometry["size_x"]["value"] / 10.0) * (geometry["size_y"]["value"] / 10.0) * (geometry["thickness"]["value"] / 10.0)
    expected_mass = volume_cm3 * geometry["density_g_cm3"] / 1000.0
    require(math.isclose(mass_kg, expected_mass, rel_tol=2.0e-12),
            f"{name}: Geant4 absorber mass differs from analytic slab mass")
    comparisons.extend([
        {"id": f"{name}_dose_ledger_conversion", "unit": "Gy/primary",
         "expected": expected_gy, "actual": output["value"], "passed": True},
        {"id": f"{name}_analytic_slab_mass", "unit": "kg",
         "expected": expected_mass, "actual": mass_kg, "passed": True},
    ])


def assert_breakdown_closure(result: dict, comparisons: list[dict]) -> None:
    raw = scores(result)
    absorber = result["absorber_id"]
    for ledger in ("total", "ionizing", "nonionizing"):
        suffix = f"/{ledger}_deposited_energy_MeV"
        aggregate = raw[f"absorber/{absorber}/{ledger}_deposited_energy_MeV"]["sum"]
        groups = {
            "particle": f"absorber/{absorber}/deposition_authoritative/by_particle/",
            "ancestry": f"absorber/{absorber}/deposition_authoritative/by_ancestry/",
            "creator_diagnostic": f"absorber/{absorber}/deposition_diagnostic_by_creator_process/",
            "step_process_diagnostic": f"absorber/{absorber}/deposition_diagnostic_by_step_defining_process/",
        }
        for group, prefix in groups.items():
            subtotal = sum(entry["sum"] for score_id, entry in raw.items()
                           if score_id.startswith(prefix) and score_id.endswith(suffix))
            require(math.isclose(subtotal, aggregate, rel_tol=2.0e-12, abs_tol=1.0e-14),
                    f"{group} {ledger} breakdown does not close to aggregate")
            comparisons.append({
                "id": f"breakdown_{group}_{ledger}_closure",
                "unit": "MeV",
                "expected": aggregate,
                "actual": subtotal,
                "passed": True,
            })
    scope = result["energy_ledger_scope"]
    require("authoritative_particle_breakdown" in scope and
            "authoritative" in scope["authoritative_particle_breakdown"],
            "particle/ancestry breakdown is not labeled authoritative")
    require("diagnostics only" in scope["diagnostic_process_breakdown"],
            "process breakdown is not explicitly labeled diagnostic")


def relative_change(a: float, b: float) -> float:
    return abs(a - b) / abs(b)


def main() -> int:
    require(EXECUTABLE.is_file(), f"missing shieldSim executable: {EXECUTABLE}")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    require(plan["registration_status"].startswith("registered_before_first"),
            "transport acceptance plan is not preregistered")
    base = json.loads((ROOT / "examples" / "stage4_tid_case.json").read_text(encoding="utf-8"))
    comparisons: list[dict] = []
    material_species_results: list[dict] = []
    all_outputs: list[tuple[str, dict]] = []
    material_aliases = {"aluminum": "Al", "silicon": "Si", "water": "Water"}
    species = {"proton": {"name": "proton", "Z": 1, "A": 1},
               "alpha": {"name": "alpha", "Z": 2, "A": 4}}

    # Full Geant4 primary generation, transport, event scoring, mass conversion,
    # and result serialization are exercised for every validated particle and
    # material. No PSTAR/ASTAR number is used as finite-volume absorbed energy:
    # escaping delta electrons and nuclear interactions make it a different
    # observable from unrestricted stopping power.
    matrix_plan=plan["material_species_matrix"]
    for case_index, registered in enumerate(matrix_plan["cases"]):
        for material_index, (material, alias) in enumerate(material_aliases.items()):
            case = copy.deepcopy(base)
            case["geometry"]["shield_stack"]["layers"] = [{
                "layer_id": "negligible_upstream_al",
                "material": "Al",
                "thickness": {"type": "geometric", "value": 0.001, "unit": "mm"},
            }]
            stack = case["geometry"]["shield_stack"]
            set_measure(stack, "size_x", matrix_plan["transverse_size_mm"])
            set_measure(stack, "size_y", matrix_plan["transverse_size_mm"])
            source_plane = case["geometry"]["source_plane"]
            set_measure(source_plane, "size_x", matrix_plan["source_size_mm"])
            set_measure(source_plane, "size_y", matrix_plan["source_size_mm"])
            absorber = case["geometry"]["absorbers"][0]
            absorber["absorber_id"] = f"{material}_slab_proxy"
            absorber["material"] = alias
            set_measure(absorber["geometry"], "size_x", matrix_plan["transverse_size_mm"])
            set_measure(absorber["geometry"], "size_y", matrix_plan["transverse_size_mm"])
            set_measure(absorber["geometry"], "thickness", matrix_plan["thickness_mm"])
            absorber["geometry"]["transform"]["translation"]["z"] = {"value": -10.0, "unit": "mm"}
            spectrum = case["source"]["spectrum"]
            spectrum["particle"] = species[registered["particle"]]
            spectrum["incident_id"] = f"{registered['particle']}_{registered['energy_MeV']:g}MeV"
            spectrum["energy"] = {"coordinate": "total_kinetic_energy",
                                  "value": registered["energy_MeV"], "unit": "MeV"}
            case["transport"]["production_cut"]["value"] = matrix_plan["production_cut_mm"]
            case["transport"]["max_step"]["value"] = matrix_plan["max_step_mm"]
            seed = 411000 + 100 * case_index + material_index
            configure_sampling(case, matrix_plan["histories"], seed)
            name = f"material_{registered['particle']}_{material}"
            result, resolved = run_case(name, case)
            output = tid_output(result)
            all_outputs.append((name, output))
            assert_exact_dose_and_mass(name, result, resolved, comparisons)
            material_species_results.append({
                "id": name,
                "particle": registered["particle"],
                "material": material,
                "energy_MeV": registered["energy_MeV"],
                "TID_Gy_per_primary": output["value"],
                "statistical_relative_standard_uncertainty": output["accuracy"]["statistical_relative_standard_uncertainty"],
                "passed": output["accuracy"]["statistical_relative_standard_uncertainty"] <= plan["monte_carlo_relative_standard_uncertainty_max"],
            })

    # Quantity-specific convergence uses the realistic example shield and a
    # one-millimetre silicon proxy. Only the final two refinements determine
    # acceptance; the coarsest setting documents the approach to the plateau.
    convergence: dict[str, list[dict]] = {}
    convergence_runs: dict[str, tuple[dict, dict]] = {}
    histories = plan["convergence"]["histories_per_setting"]
    for family, values in (("production_cut_mm", plan["convergence"]["production_cut_mm"]),
                           ("max_step_mm", plan["convergence"]["max_step_mm"]),
                           ("transverse_size_mm", plan["convergence"]["transverse_size_mm"])):
        points: list[dict] = []
        for index, value in enumerate(values):
            case = copy.deepcopy(base)
            configure_sampling(case, histories, 412000)
            # Keep the source compact and centered while expanding the material
            # boundaries; this isolates lateral-escape convergence from source
            # sampling changes.
            set_measure(case["geometry"]["source_plane"], "size_x", 10.0)
            set_measure(case["geometry"]["source_plane"], "size_y", 10.0)
            if family == "production_cut_mm":
                case["transport"]["production_cut"]["value"] = value
            elif family == "max_step_mm":
                case["transport"]["production_cut"]["value"] = 0.03
                case["transport"]["max_step"]["value"] = value
            else:
                case["transport"]["production_cut"]["value"] = 0.03
                case["transport"]["max_step"]["value"] = 0.05
                for key in ("size_x", "size_y"):
                    set_measure(case["geometry"]["shield_stack"], key, value)
                    set_measure(case["geometry"]["absorbers"][0]["geometry"], key, value)
            name = f"convergence_{family}_{value:g}"
            result, resolved = run_case(name, case)
            convergence_runs[name] = (result, resolved)
            output = tid_output(result)
            all_outputs.append((name, output))
            assert_exact_dose_and_mass(name, result, resolved, comparisons)
            aggregate=scores(result)[f"absorber/{result['absorber_id']}/ionizing_deposited_energy_MeV"]
            mean_ionizing_mev=aggregate["sum"]/result["events"]
            metric_value=mean_ionizing_mev if family=="transverse_size_mm" else output["value"]
            metric_name=(plan["convergence"]["transverse_size_metric"]
                         if family=="transverse_size_mm" else "TID_Gy_per_primary")
            points.append({"setting": value, "TID_Gy_per_primary": output["value"],
                           "ionizing_deposited_energy_MeV_per_primary":mean_ionizing_mev,
                           "convergence_metric":metric_name,
                           "convergence_value":metric_value,
                           "statistical_relative_standard_uncertainty": output["accuracy"]["statistical_relative_standard_uncertainty"]})
        final_change = relative_change(points[-2]["convergence_value"], points[-1]["convergence_value"])
        points[-1]["final_refinement_relative_change"] = final_change
        points[-1]["passed"] = final_change <= plan["final_refinement_relative_change_max"]
        require(points[-1]["passed"], f"{family}: final TID refinement changed by {final_change:.4g}")
        convergence[family] = points

    final_name = f"convergence_transverse_size_mm_{plan['convergence']['transverse_size_mm'][-1]:g}"
    final_result, final_resolved = convergence_runs[final_name]
    assert_breakdown_closure(final_result, comparisons)
    require(tid_output(final_result)["validation_status"].endswith("pending_independent_review"),
            "Stage-4 result status does not preserve the independent-review boundary")

    # Compare the approved production baseline against option 4 with identical
    # seed, histories, geometry, cuts, and step limit.
    option_case = json.loads((RUN_DIR / final_name / "case.json").read_text(encoding="utf-8"))
    option_case["transport"]["em_physics"] = "G4EmStandardPhysics_option4"
    option_result, option_resolved = run_case("em_option4_sensitivity", option_case)
    option_output = tid_output(option_result)
    all_outputs.append(("em_option4_sensitivity", option_output))
    assert_exact_dose_and_mass("em_option4_sensitivity", option_result, option_resolved, comparisons)
    baseline_value = tid_output(final_result)["value"]
    em_sensitivity = relative_change(option_output["value"], baseline_value)
    require(em_sensitivity <= plan["em_option4_transport_sensitivity_relative_max"],
            f"EM option-4 TID sensitivity {em_sensitivity:.4g} exceeds registered allocation")

    # Independent seeds verify the expected 1/sqrt(N) reduction of the reported
    # event-level standard error. The ratio is evaluated from medians so a
    # single unusually quiet realization cannot decide the gate.
    scaling = plan["convergence"]["history_scaling"]
    relative_se: dict[str, list[float]] = {"low": [], "high": []}
    scaling_points: list[dict] = []
    for label, count in (("low", scaling["low_histories"]),
                         ("high", scaling["high_histories"])):
        for seed in scaling["seeds"]:
            case = copy.deepcopy(base)
            configure_sampling(case, count, seed)
            case["transport"]["production_cut"]["value"] = 0.03
            case["transport"]["max_step"]["value"] = 0.05
            set_measure(case["geometry"]["source_plane"], "size_x", 10.0)
            set_measure(case["geometry"]["source_plane"], "size_y", 10.0)
            for key in ("size_x", "size_y"):
                set_measure(case["geometry"]["shield_stack"], key, 80.0)
                set_measure(case["geometry"]["absorbers"][0]["geometry"], key, 80.0)
            name = f"history_{label}_{count}_seed_{seed}"
            result, _ = run_case(name, case)
            output = tid_output(result)
            value = output["accuracy"]["statistical_relative_standard_uncertainty"]
            relative_se[label].append(value)
            scaling_points.append({"group": label, "histories": count, "seed": seed,
                                   "statistical_relative_standard_uncertainty": value})
    se_ratio = statistics.median(relative_se["high"]) / statistics.median(relative_se["low"])
    ratio_interval = scaling["accepted_ratio_interval"]
    require(ratio_interval[0] <= se_ratio <= ratio_interval[1],
            f"history scaling ratio {se_ratio:.4g} is outside {ratio_interval}")

    # The Stage-4 acceptance population must meet the approved <=2% MC standard
    # uncertainty independently of the broader default 10% combined target.
    mc_limit = plan["monte_carlo_relative_standard_uncertainty_max"]
    for name, output in all_outputs:
        statistical = output["accuracy"]["statistical_relative_standard_uncertainty"]
        require(statistical <= mc_limit,
                f"{name}: MC relative standard uncertainty {statistical:.4g} exceeds {mc_limit}")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_id": "stage4_tid_transport",
        "status": "passed",
        "scope": "proton/helium-4 finite material slabs; not device, organ, or crew dose",
        "comparisons": comparisons,
        "material_species_matrix": material_species_results,
        "convergence": convergence,
        "history_scaling": {"points": scaling_points, "median_high_to_low_se_ratio": se_ratio,
                            "accepted_ratio_interval": ratio_interval},
        "em_option4_transport_sensitivity_relative": em_sensitivity,
        "em_option4_transport_sensitivity_relative_max": plan["em_option4_transport_sensitivity_relative_max"],
        "monte_carlo_relative_standard_uncertainty_max": mc_limit,
    }
    (EVIDENCE / "stage4_tid_transport.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Stage-4 TID transport passed: {len(all_outputs)} runs, {len(comparisons)} exact/reference comparisons")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Strict Stage-6 capability and LET-response schema boundaries."""

from __future__ import annotations

import copy
import json
import pathlib
import subprocess
import sys

if len(sys.argv) != 4:
    raise SystemExit("usage: test_stage6_case_boundaries.py SHIELDSIM SOURCE_ROOT RUN_DIR")

EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
RUN_DIR = pathlib.Path(sys.argv[3]).resolve()


def base() -> dict:
    case = json.loads((ROOT / "examples/stage4_tid_case.json").read_text(encoding="utf-8"))
    case["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    case["quantities"] = ["LET"]
    case["let_response"] = {
        "definition_id": "shieldSim.electronic_let.v1",
        "linear_let_bins": {"scale": "logarithmic",
                            "minimum": {"value": 0.001, "unit": "keV/um"},
                            "maximum": {"value": 10000.0, "unit": "keV/um"},
                            "bins_per_decade": 20},
        "out_of_range_policy": "not_assessable_no_clamping",
    }
    case["accuracy"]["scope"] = "statistical"
    case["accuracy"]["coverage"]["interval_method"] = "fieller_event_covariance"
    case["accuracy"]["systematic_budget"] = {"status": "incomplete"}
    return case


def validate(name: str, document: dict, accepted: bool) -> None:
    path = RUN_DIR / f"{name}.json"
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    completed = subprocess.run([str(EXE), f"--validate-case={path}"], cwd=RUN_DIR,
                               text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               check=False)
    (RUN_DIR / f"{name}.log").write_text(completed.stdout, encoding="utf-8")
    if (completed.returncode == 0) != accepted:
        raise AssertionError(f"{name}: expected accepted={accepted}; see {name}.log")


def main() -> int:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    valid = base()
    validate("valid_let", valid, True)
    missing = copy.deepcopy(valid)
    missing.pop("let_response")
    validate("missing_let_response", missing, False)
    bad_unit = copy.deepcopy(valid)
    bad_unit["let_response"]["linear_let_bins"]["minimum"]["unit"] = "MeV cm2/g"
    validate("wrong_primary_let_unit", bad_unit, False)
    noninteger = copy.deepcopy(valid)
    noninteger["let_response"]["linear_let_bins"]["maximum"]["value"] = 9000.0
    validate("noninteger_log_bin_span", noninteger, False)
    unused = copy.deepcopy(valid)
    unused["quantities"] = ["H100/10"]
    validate("unused_let_response", unused, False)
    hardness = copy.deepcopy(unused)
    hardness.pop("let_response")
    validate("valid_hardness_without_let_response", hardness, True)
    gamma = copy.deepcopy(valid)
    gamma["source"]["spectrum"]["particle"] = {"name": "gamma", "Z": 0, "A": 0}
    validate("gamma_outside_let_incident_domain", gamma, False)
    carbon_al = copy.deepcopy(valid)
    carbon_al["source"]["spectrum"]["particle"] = {"name": "ion", "Z": 6, "A": 12}
    carbon_al["source"]["spectrum"]["energy"] = {
        "coordinate": "kinetic_energy_per_nucleon", "value": 100.0, "unit": "MeV/n"
    }
    carbon_al["geometry"]["absorbers"][0]["material"] = "Al"
    validate("c12_non_silicon_reference_domain", carbon_al, False)
    carbon_default = copy.deepcopy(valid)
    carbon_default["source"]["spectrum"]["particle"] = {"name": "ion", "Z": 6, "A": 12}
    carbon_default["source"]["spectrum"]["energy"] = {
        "coordinate": "kinetic_energy_per_nucleon", "value": 100.0, "unit": "MeV/n"
    }
    validate("c12_default_em_model_outside_reference_screen", carbon_default, False)
    carbon_option4 = copy.deepcopy(carbon_default)
    carbon_option4["transport"]["em_physics"] = "G4EmStandardPhysics_option4"
    validate("c12_option4_validated_model", carbon_option4, True)
    print("All Stage-6 case boundaries passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

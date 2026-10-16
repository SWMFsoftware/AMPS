#!/usr/bin/env python3
"""Schema and boundary tests for the Stage-1 case-file executable path.

The test derives every case from the documented example, changes one contract
at a time, and calls --validate-case.  No particle transport is needed, so this
suite is quick enough for every build while still exercising the exact CLI and
production parser used by users.
"""

from __future__ import annotations

import copy
import json
import pathlib
import subprocess
import sys
import tempfile


if len(sys.argv) not in (3, 4):
    raise SystemExit("usage: test_stage1_case_boundaries.py SHIELDSIM SOURCE_ROOT [EVIDENCE_DIR]")
EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE_DIR = pathlib.Path(sys.argv[3]).resolve() if len(sys.argv) == 4 else None
BASE = json.loads((ROOT / "examples/stage1_case.json").read_text())
BASE["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
BASE["source"]["spectrum"]["path"] = str(ROOT / "tests/data/mono_100MeV_proton.dat")
BASE["source"]["spectrum"]["metadata_path"] = str(
    ROOT / "tests/data/mono_spectrum_metadata.v1.json"
)


def run_case(work: pathlib.Path, name: str, case: dict, should_pass: bool) -> dict[str, dict]:
    path = work / f"{name}.json"
    case = copy.deepcopy(case)
    case["case_id"] = name
    case["output"]["prefix"] = str(work / name)
    path.write_text(json.dumps(case, indent=2) + "\n")
    completed = subprocess.run(
        [str(EXE), f"--validate-case={path}"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if (completed.returncode == 0) != should_pass:
        raise AssertionError(
            f"{name}: expected {'success' if should_pass else 'failure'}, "
            f"got rc={completed.returncode}\n{completed.stdout}"
        )
    if not should_pass:
        return {}
    result: dict[str, dict] = {}
    for absorber in case["geometry"]["absorbers"]:
        manifest = work / f"{name}__{absorber['absorber_id']}_resolved_case.json"
        result[absorber["absorber_id"]] = json.loads(manifest.read_text())
    return result


def main() -> int:
    comparisons: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="shieldsim-stage1-") as directory:
        work = pathlib.Path(directory)

        # Inclusive layer-count boundaries: one and five pass; six fails.
        one = copy.deepcopy(BASE)
        one["geometry"]["shield_stack"]["layers"] = one["geometry"]["shield_stack"]["layers"][:1]
        run_case(work, "one_layer", one, True)

        five = copy.deepcopy(BASE)
        prototype = five["geometry"]["shield_stack"]["layers"][0]
        five["geometry"]["shield_stack"]["layers"] = []
        for index in range(5):
            layer = copy.deepcopy(prototype)
            layer["layer_id"] = f"layer_{index}"
            layer["thickness"] = {"type": "geometric", "value": 1.0, "unit": "mm"}
            five["geometry"]["shield_stack"]["layers"].append(layer)
        for absorber in five["geometry"]["absorbers"]:
            absorber["geometry"]["front_z"]["value"] = -14.0
        run_case(work, "five_layers", five, True)

        six = copy.deepcopy(five)
        extra = copy.deepcopy(prototype)
        extra["layer_id"] = "layer_5"
        extra["thickness"] = {"type": "geometric", "value": 1.0, "unit": "mm"}
        six["geometry"]["shield_stack"]["layers"].append(extra)
        run_case(work, "six_layers_rejected", six, False)

        # A homogeneous layer and the same material split into adjacent layers
        # have exactly equal analytic thickness and areal-density ledgers.
        single_al = copy.deepcopy(BASE)
        single_al["geometry"]["shield_stack"]["layers"] = [{
            "layer_id": "al_all", "material": "Al",
            "thickness": {"type": "geometric", "value": 2.0, "unit": "mm"},
        }]
        for absorber in single_al["geometry"]["absorbers"]:
            absorber["geometry"]["front_z"]["value"] = -17.0
        split_al = copy.deepcopy(single_al)
        split_al["geometry"]["shield_stack"]["layers"] = [
            {"layer_id": "al_a", "material": "Al",
             "thickness": {"type": "geometric", "value": 1.0, "unit": "mm"}},
            {"layer_id": "al_b", "material": "Al",
             "thickness": {"type": "geometric", "value": 1.0, "unit": "mm"}},
        ]
        single_result = run_case(work, "single_al", single_al, True)["water_proxy"]
        split_result = run_case(work, "split_al", split_al, True)["water_proxy"]
        single_layers = single_result["geometry"]["shield_layers"]
        split_layers = split_result["geometry"]["shield_layers"]
        assert sum(x["thickness"]["value"] for x in single_layers) == \
               sum(x["thickness"]["value"] for x in split_layers)
        assert sum(x["areal_density"]["value"] for x in single_layers) == \
               sum(x["areal_density"]["value"] for x in split_layers)
        comparisons.extend([
            {
                "id": "single_vs_split_total_thickness",
                "expected": sum(x["thickness"]["value"] for x in single_layers),
                "actual": sum(x["thickness"]["value"] for x in split_layers),
                "unit": "mm",
            },
            {
                "id": "single_vs_split_areal_density",
                "expected": sum(x["areal_density"]["value"] for x in single_layers),
                "actual": sum(x["areal_density"]["value"] for x in split_layers),
                "unit": "g/cm2",
            },
        ])

        # Layer order is first-class geometry, not a material-set operation.
        order_a = run_case(work, "layer_order_a", BASE, True)["water_proxy"]
        order_case = copy.deepcopy(BASE)
        order_case["geometry"]["shield_stack"]["layers"].reverse()
        order_b = run_case(work, "layer_order_b", order_case, True)["water_proxy"]
        assert order_a["geometry_hash"] != order_b["geometry_hash"]

        # Inclusive 10-cm per-layer boundary. Enlarge the explicitly declared
        # world and move the absorber so geometry bounds do not mask this test.
        at_limit = copy.deepcopy(BASE)
        at_limit["geometry"]["world"]["size_z"]["value"] = 300.0
        at_limit["geometry"]["source_plane"]["z"]["value"] = -140.0
        at_limit["geometry"]["shield_stack"]["front_z"]["value"] = -100.0
        at_limit["geometry"]["shield_stack"]["layers"] = [copy.deepcopy(prototype)]
        at_limit["geometry"]["shield_stack"]["layers"][0]["thickness"] = {
            "type": "geometric", "value": 100.0, "unit": "mm"
        }
        for absorber in at_limit["geometry"]["absorbers"]:
            absorber["geometry"]["front_z"]["value"] = 1.0
        limit_result = run_case(work, "thickness_100mm", at_limit, True)["water_proxy"]
        comparisons.append({
            "id": "inclusive_layer_thickness_boundary",
            "expected": 100.0,
            "actual": limit_result["geometry"]["shield_layers"][0]["thickness"]["value"],
            "unit": "mm",
        })
        over_limit = copy.deepcopy(at_limit)
        over_limit["geometry"]["shield_stack"]["layers"][0]["thickness"]["value"] = 100.001
        run_case(work, "thickness_over_100mm", over_limit, False)

        missing_dimension = copy.deepcopy(BASE)
        del missing_dimension["geometry"]["world"]["size_x"]
        run_case(work, "missing_world_size", missing_dimension, False)
        missing_placement = copy.deepcopy(BASE)
        del missing_placement["geometry"]["shield_stack"]["center_x"]
        run_case(work, "missing_shield_placement", missing_placement, False)
        missing_absorber_placement = copy.deepcopy(BASE)
        del missing_absorber_placement["geometry"]["absorbers"][0]["geometry"]["center_y"]
        run_case(work, "missing_absorber_placement", missing_absorber_placement, False)
        bad_unit = copy.deepcopy(BASE)
        bad_unit["geometry"]["source_plane"]["size_x"]["unit"] = "inch"
        run_case(work, "unknown_length_unit", bad_unit, False)
        unknown_field = copy.deepcopy(BASE)
        unknown_field["geometry"]["shield_stack"]["typo_thickness"] = 2.0
        run_case(work, "unknown_field", unknown_field, False)

        reversed_range = copy.deepcopy(BASE)
        reversed_range["source"]["energy_limits"]["proton"]["min"]["value"] = 101.0
        run_case(work, "reversed_energy_range", reversed_range, False)
        unsupported = copy.deepcopy(BASE)
        unsupported["quantities"] = ["n_eq"]
        run_case(work, "stage1_unsupported_quantity", unsupported, False)

        # Default/override precedence is visible in serialized provenance.
        default_result = run_case(work, "accuracy_default", BASE, True)["water_proxy"]
        assert default_result["accuracy"]["target"] == 0.10
        assert default_result["accuracy"]["target_source"] == "default_0.10"
        overridden = copy.deepcopy(BASE)
        overridden["accuracy"]["target"] = 0.20
        overridden["accuracy"]["quantity_overrides"] = {
            "TID": {
                "target": 0.05,
                "scope": "statistical",
                "coverage": copy.deepcopy(BASE["accuracy"]["coverage"]),
            }
        }
        override_result = run_case(work, "accuracy_override", overridden, True)["water_proxy"]
        assert override_result["accuracy"]["target"] == 0.20
        assert override_result["accuracy"]["quantity_overrides"]["TID"]["target"] == 0.05

        invalid_accuracy = copy.deepcopy(BASE)
        invalid_accuracy["accuracy"]["target"] = 0.0
        run_case(work, "zero_accuracy", invalid_accuracy, False)

        # Absorber permutation changes the enclosing case serialization but
        # must not change either absorber's independent geometry identity.
        normal = run_case(work, "absorber_order_a", BASE, True)
        permuted = copy.deepcopy(BASE)
        permuted["geometry"]["absorbers"].reverse()
        reversed_result = run_case(work, "absorber_order_b", permuted, True)
        for absorber_id in normal:
            assert normal[absorber_id]["active_absorber_geometry_hash"] == \
                   reversed_result[absorber_id]["active_absorber_geometry_hash"]
            assert normal[absorber_id]["geometry_hash"] == \
                   reversed_result[absorber_id]["geometry_hash"]

        # A density override changes the mass-thickness conversion and both the
        # material and geometry/case hashes; the parent library remains intact.
        density = copy.deepcopy(BASE)
        density["material_overrides"] = [{
            "material_id": "Al_density_test", "revision": "case-1",
            "extends": "Al_natural_G4_11_4@1", "status": "project_specific",
            "density": {"value": 2.70, "unit": "g/cm3",
                        "uncertainty": {"type": "standard", "value": 0.01, "unit": "g/cm3"},
                        "reference_condition": "test measurement"},
            "provenance": [{"properties": ["density"], "source": "test record", "version": "1"}]
        }]
        density["geometry"]["shield_stack"]["layers"][0]["material"] = "Al_density_test"
        density_result = run_case(work, "density_override", density, True)["water_proxy"]
        default_layer = default_result["geometry"]["shield_layers"][0]
        density_layer = density_result["geometry"]["shield_layers"][0]
        assert density_layer["thickness"]["value"] != default_layer["thickness"]["value"]
        assert density_layer["material_hash"] != default_layer["material_hash"]
        assert density_result["geometry_hash"] != default_result["geometry_hash"]

        # The input areal density is unchanged, so the expected geometric
        # thickness follows t[mm] = 10*x[g/cm2]/rho[g/cm3].
        comparisons.extend([
            {
                "id": "default_density_areal_to_geometric",
                "expected": 10.0 * default_layer["areal_density"]["value"] / 2.699,
                "actual": default_layer["thickness"]["value"],
                "unit": "mm",
            },
            {
                "id": "override_density_areal_to_geometric",
                "expected": 10.0 * density_layer["areal_density"]["value"] / 2.70,
                "actual": density_layer["thickness"]["value"],
                "unit": "mm",
            },
        ])

    if EVIDENCE_DIR is not None:
        EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
        (EVIDENCE_DIR / "stage1_case_boundaries.json").write_text(
            json.dumps(
                {
                    "schema_version": "shieldSim.validation_evidence.v1",
                    "test_case": "stage1_case_boundaries",
                    "comparisons": comparisons,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )

    print("All Stage-1 case boundary tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Execute every user-facing material/species/shield/shape example.

This is a configuration-plumbing gate, not a precision or physical-validation
comparison.  Each committed JSON is run through the normal parser, Geant4
transport, and quantity writers.  It also checks the important semantic split:
proton all-quantity cases report the defined zero neutron-only n_eq, while the
dedicated neutron cases record a positive n_eq entry response.
"""

from __future__ import annotations

import hashlib
import json
import pathlib
import subprocess
import sys


if len(sys.argv) != 5:
    raise SystemExit(
        "usage: test_configuration_examples.py SHIELDSIM SOURCE_ROOT "
        "EVIDENCE_DIR RUN_DIR"
    )

EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
EXAMPLES = ROOT / "examples" / "configuration_cases"


EXPECTED = {
    "01_si_box_proton_all_quantities": {
        "species": {"proton"},
        "material": "Si",
        "provider": "analytic_csg",
        "layers": 3,
        "quantities": {"TID", "DDD", "n_eq", "LET", "H100/10"},
        "n_eq": "zero",
    },
    "02_si_cylinder_neutron_damage": {
        "species": {"neutron"},
        "material": "Si",
        "provider": "analytic_csg",
        "layers": 1,
        "quantities": {"DDD", "n_eq"},
        "n_eq": "positive",
    },
    "03_si_sphere_alpha_tid_ddd_let": {
        "species": {"alpha"},
        "material": "Si",
        "provider": "analytic_csg",
        "layers": 5,
        "quantities": {"TID", "DDD", "LET"},
    },
    "04_water_cylinder_proton_tid_let_h100": {
        "species": {"proton"},
        "material": "Water",
        "provider": "analytic_csg",
        "layers": 2,
        "quantities": {"TID", "LET", "H100/10"},
    },
    "05_al_mesh_box_proton_tid_let_h100": {
        "species": {"proton"},
        "material": "Al",
        "provider": "tessellated_mesh",
        "layers": 4,
        "quantities": {"TID", "LET", "H100/10"},
    },
    "06_si_mesh_tetrahedron_neutron_damage": {
        "species": {"neutron"},
        "material": "Si",
        "provider": "tessellated_mesh",
        "layers": 2,
        "quantities": {"DDD", "n_eq"},
        "n_eq": "positive",
    },
    "07_si_mesh_octahedron_case_material_override": {
        "species": {"proton"},
        "material": "Si",
        "provider": "tessellated_mesh",
        "layers": 2,
        "quantities": {"TID", "DDD", "n_eq", "LET", "H100/10"},
        "n_eq": "zero",
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sensitive_region(document: dict) -> dict:
    regions = document["geometry"]["absorbers"][0]["geometry"]["regions"]
    selected = [region for region in regions if region["role"] == "sensitive"]
    require(len(selected) == 1, "example must define exactly one sensitive region")
    return selected[0]


def quantity(payload: dict, name: str) -> dict:
    selected = [item for item in payload["outputs"] if item["name"] == name]
    require(len(selected) == 1, f"expected exactly one {name} output")
    return selected[0]


def main() -> int:
    require(EXE.is_file(), f"missing shieldSim executable {EXE}")
    paths = sorted(EXAMPLES.glob("*.json"))
    require({path.stem for path in paths} == set(EXPECTED),
            "configuration gallery is incomplete or contains an untested JSON file")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    evidence_cases = []

    for path in paths:
        name = path.stem
        expected = EXPECTED[name]
        document = json.loads(path.read_text(encoding="utf-8"))
        absorber = document["geometry"]["absorbers"][0]
        region = sensitive_region(document)
        stack = document["geometry"]["shield_stack"]
        nodes = document["source"]["spectrum"]["nodes"]

        # Check the committed user-facing selections before executing the case.
        require({node["particle"]["name"] for node in nodes} == expected["species"],
                f"{name}: incident-species selection changed")
        require(region["material"] == expected["material"],
                f"{name}: absorber material selection changed")
        require(absorber["geometry"]["type"] == expected["provider"],
                f"{name}: shape-provider selection changed")
        require(len(stack["layers"]) == expected["layers"],
                f"{name}: ordered layer count changed")
        require(1 <= len(stack["layers"]) <= 5,
                f"{name}: layer count is outside the model boundary")
        require(stack["size_x"]["value"] > 0 and stack["size_y"]["value"] > 0,
                f"{name}: shield footprint must be explicit and positive")
        require({layer["thickness"]["type"] for layer in stack["layers"]}
                <= {"geometric", "areal_density"},
                f"{name}: unknown layer thickness representation")
        require(set(document["quantities"]) == expected["quantities"],
                f"{name}: calculation selection changed")
        require(abs(sum(node["weight"] for node in nodes) - 1.0) <= 1.0e-12,
                f"{name}: source weights are not normalized")

        directory = RUN_DIR / name
        directory.mkdir(parents=True, exist_ok=True)
        completed = subprocess.run(
            [str(EXE), f"--case={path}"], cwd=directory, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
        )
        (directory / "run.log").write_text(completed.stdout, encoding="utf-8")
        require(completed.returncode == 0,
                f"{name}: transport failed; see {directory / 'run.log'}")

        prefix = document["output"]["prefix"]
        absorber_id = absorber["absorber_id"]
        result_path = directory / f"{prefix}__{absorber_id}_result.json"
        resolved_path = directory / f"{prefix}__{absorber_id}_resolved_case.json"
        require(result_path.is_file(), f"{name}: result JSON is missing")
        require(resolved_path.is_file(), f"{name}: resolved case is missing")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        resolved = json.loads(resolved_path.read_text(encoding="utf-8"))
        require({item["name"] for item in result["outputs"]} == expected["quantities"],
                f"{name}: runtime output set differs from input selection")
        require(set(resolved["quantities"]) == expected["quantities"],
                f"{name}: resolved case lost the quantity selection")

        if expected.get("n_eq") == "zero":
            require(quantity(result, "n_eq")["value"] == 0.0,
                    f"{name}: proton-only neutron-equivalent fluence must be zero")
        elif expected.get("n_eq") == "positive":
            require(quantity(result, "n_eq")["value"] > 0.0,
                    f"{name}: neutron entry case did not produce positive n_eq")

        if name == "07_si_mesh_octahedron_case_material_override":
            layers = resolved["geometry"]["shield_layers"]
            require(layers[0]["material_id"] == "Al_project_density_example",
                    "case-local material ID was not used by the shield")
            require(abs(layers[0]["thickness"]["value"] - 1.0) <= 1.0e-12,
                    "0.27 g/cm2 at 2.70 g/cm3 did not resolve to 1.0 mm")

        evidence_cases.append({
            "case": name,
            "incident_species": sorted(expected["species"]),
            "absorber_material": expected["material"],
            "geometry_provider": expected["provider"],
            "shield_size_mm": [stack["size_x"]["value"], stack["size_y"]["value"]],
            "ordered_layer_count": len(stack["layers"]),
            "quantities": sorted(expected["quantities"]),
            "result_sha256": "sha256:" + hashlib.sha256(result_path.read_bytes()).hexdigest(),
            "passed": True,
        })

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "scope": "working_input_configuration_examples",
        "physical_validation_claim": False,
        "note": (
            "This gate proves input-to-transport/output plumbing. Physical and "
            "analytic comparison claims remain in the stage-specific gates."
        ),
        "cases": evidence_cases,
    }
    evidence_path = EVIDENCE / "configuration_examples.json"
    evidence_path.write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"Configuration gallery passed {len(evidence_cases)} complete transport runs; "
        f"evidence={evidence_path}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL configuration examples: {error}", file=sys.stderr)
        raise SystemExit(1)

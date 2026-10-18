#!/usr/bin/env python3
"""Stage-10 analytic-CSG parser, geometry, transport, and cache-key gate.

The first arbitrary-shape provider is intentionally verification-only.  This
test compares Geant4 primitive masses with exact analytic values, proves that a
CSG box is transport-equivalent to the legacy slab, exercises a nontrivial
Boolean enclosure, checks rigid-rotation invariance for a square box, and
enforces fatal schema/overlap boundaries.  It does not supply the independent
nontrivial external reference required to complete Stage 10.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import pathlib
import subprocess
import sys


if len(sys.argv) != 6:
    raise SystemExit(
        "usage: test_stage10_analytic_csg.py SHIELDSIM SOURCE_ROOT EVIDENCE_DIR "
        "RUN_DIR FIGURE_DIR"
    )

EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
FIGURES = pathlib.Path(sys.argv[5]).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def transform(x: float, y: float, z: float,
              quaternion: tuple[float, float, float, float] = (1.0, 0.0, 0.0, 0.0)) -> dict:
    return {
        "translation": {
            "x": {"value": x, "unit": "mm"},
            "y": {"value": y, "unit": "mm"},
            "z": {"value": z, "unit": "mm"},
        },
        "rotation_quaternion_local_to_world": dict(
            zip(("w", "x", "y", "z"), quaternion)
        ),
    }


def mass_properties(volume_cm3: float, provenance: str,
                    tolerance: float = 1.0e-9,
                    centroid_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)) -> dict:
    return {
        "volume": {"value": volume_cm3, "unit": "cm3"},
        "centroid": {
            "x": {"value": centroid_mm[0], "unit": "mm"},
            "y": {"value": centroid_mm[1], "unit": "mm"},
            "z": {"value": centroid_mm[2], "unit": "mm"},
        },
        "relative_tolerance": tolerance,
        "provenance": provenance,
    }


def csg_absorber(solid: dict, volume_cm3: float, z_mm: float = -15.5,
                 quaternion: tuple[float, float, float, float] = (1.0, 0.0, 0.0, 0.0),
                 absorber_id: str = "silicon_proxy",
                 entry_surface: str = "local_negative_z_face") -> dict:
    return {
        "absorber_id": absorber_id,
        "geometry": {
            "type": "analytic_csg",
            "provider_version": "1.0.0",
            "coordinate_frame": {
                "name": "absorber_local_cartesian",
                "handedness": "right",
                "origin": "provider_origin",
                "length_unit": "mm",
            },
            "transform": transform(0.0, 0.0, z_mm, quaternion),
            "source": {
                "kind": "inline",
                "source_id": "stage10_synthetic_analytic_fixture",
                "source_version": "1.0.0",
                "provenance": "Synthetic verification primitive; not an as-built spacecraft model",
            },
            "regions": [{
                "region_id": "detector",
                "role": "sensitive",
                "material": "Si",
                "transform": transform(0.0, 0.0, 0.0),
                "solid": solid,
                "entry_surface": entry_surface,
                "mass_properties": mass_properties(
                    volume_cm3, "Exact closed-form analytic primitive volume and centroid"
                ),
            }],
        },
    }


def base_case() -> dict:
    case = json.loads((ROOT / "examples/stage4_tid_case.json").read_text(encoding="utf-8"))
    case["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    case["accuracy"]["sampling"] = {
        "min_histories": 40, "max_histories": 40, "batch_size": 40
    }
    case["transport"]["random_seed"] = 1010010
    return case


def write_case(directory: pathlib.Path, case: dict) -> pathlib.Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "case.json"
    path.write_text(json.dumps(case, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def invoke(name: str, case: dict, validate_only: bool = False,
           expect_success: bool = True) -> tuple[subprocess.CompletedProcess[str], pathlib.Path]:
    directory = RUN_DIR / name
    document = copy.deepcopy(case)
    document["case_id"] = f"stage10_{name}"
    document["output"]["prefix"] = name
    path = write_case(directory, document)
    option = f"--validate-case={path}" if validate_only else f"--case={path}"
    completed = subprocess.run(
        [str(EXE), option], cwd=directory, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    (directory / "run.log").write_text(completed.stdout, encoding="utf-8")
    require((completed.returncode == 0) == expect_success,
            f"{name}: expected success={expect_success}; see {directory/'run.log'}")
    return completed, directory


def result(directory: pathlib.Path, prefix: str, absorber_id: str = "silicon_proxy") -> dict:
    path = directory / f"{prefix}__{absorber_id}_result.json"
    require(path.is_file(), f"missing result: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def resolved(directory: pathlib.Path, prefix: str, absorber_id: str = "silicon_proxy") -> dict:
    path = directory / f"{prefix}__{absorber_id}_resolved_case.json"
    require(path.is_file(), f"missing resolved case: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def tid(payload: dict) -> dict:
    values = [entry for entry in payload["outputs"] if entry["name"] == "TID"]
    require(len(values) == 1, "expected exactly one TID output")
    return values[0]


def quantity(payload: dict, name: str) -> dict:
    values = [entry for entry in payload["outputs"] if entry["name"] == name]
    require(len(values) == 1, f"expected exactly one {name} output")
    return values[0]


def add_comparison(comparisons: list[dict], name: str, expected: float,
                   actual: float, unit: str, tolerance: float) -> None:
    relative = abs(actual - expected) / max(abs(expected), 1.0e-300)
    require(relative <= tolerance,
            f"{name}: relative difference {relative} exceeds {tolerance}")
    comparisons.append({
        "id": name, "expected": expected, "actual": actual, "unit": unit,
        "relative_tolerance": tolerance, "passed": True,
    })


def render_figure(comparisons: list[dict]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise AssertionError(f"Stage-10 comparison figure requires matplotlib: {error}")
    FIGURES.mkdir(parents=True, exist_ok=True)
    selected = [item for item in comparisons if "mass" in item["id"] or
                "equivalence" in item["id"] or "invariance" in item["id"]]
    labels = [item["id"].replace("_", "\n") for item in selected]
    normalized_expected = []
    normalized_actual = []
    for item in selected:
        scale = max(abs(item["expected"]), abs(item["actual"]), 1.0e-300)
        normalized_expected.append(item["expected"] / scale)
        normalized_actual.append(item["actual"] / scale)
    positions = range(len(selected))
    figure, axis = plt.subplots(figsize=(12, 6.5))
    axis.bar([value - 0.2 for value in positions], normalized_expected,
             width=0.4, label="analytic / legacy reference")
    axis.bar([value + 0.2 for value in positions], normalized_actual,
             width=0.4, label="analytic_csg.v1 / Geant4")
    axis.set_xticks(list(positions))
    axis.set_xticklabels(labels, fontsize=8)
    axis.set_ylabel("value normalized within each comparison")
    axis.set_title("Stage 10 analytic-CSG verification comparisons")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.text(
        0.01, 0.01,
        "Verification-only: primitive mass, slab/box transport equivalence, and rigid-rotation "
        "invariance. No independent nontrivial spacecraft-shape reference is claimed.",
        fontsize=8,
    )
    figure.tight_layout(rect=(0.0, 0.06, 1.0, 1.0))
    for extension in ("png", "svg"):
        path = FIGURES / f"stage10_analytic_csg_comparisons.{extension}"
        figure.savefig(path, dpi=170 if extension == "png" else None,
                       bbox_inches="tight")
        require(path.stat().st_size > 1000, f"unexpectedly small figure: {path}")
    plt.close(figure)


def main() -> int:
    require(EXE.is_file(), f"missing shieldSim executable: {EXE}")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    comparisons: list[dict] = []
    base = base_case()

    # Legacy slab and an exactly coextensive analytic box must generate the
    # same transport sequence with the same seed. This isolates provider
    # representation from source and physics-model differences.
    _, slab_dir = invoke("slab_reference", base)
    slab_result = result(slab_dir, "slab_reference")
    box = {"node_id": "detector_box", "type": "box",
           "size_x": {"value": 50.0, "unit": "mm"},
           "size_y": {"value": 50.0, "unit": "mm"},
           "size_z": {"value": 1.0, "unit": "mm"}}
    box_case = copy.deepcopy(base)
    box_case["geometry"]["absorbers"] = [csg_absorber(box, 2.5)]
    _, box_dir = invoke("csg_box", box_case)
    box_result = result(box_dir, "csg_box")
    add_comparison(comparisons, "box_analytic_mass_kg", 2.5 * 2.33 / 1000.0,
                   tid(box_result)["absorber_mass"]["value"], "kg", 2.0e-12)
    add_comparison(comparisons, "slab_csg_box_mass_equivalence_kg",
                   tid(slab_result)["absorber_mass"]["value"],
                   tid(box_result)["absorber_mass"]["value"], "kg", 2.0e-12)
    add_comparison(comparisons, "slab_csg_box_tid_equivalence_Gy_per_primary",
                   tid(slab_result)["value"], tid(box_result)["value"],
                   "Gy/primary", 2.0e-12)
    require("stage10_geometry_providers_verification_only" in
            tid(box_result)["validation_status"],
            "analytic CSG output overclaims a completed validation stage")

    # DDD/n_eq use the provider-defined entry surface, not deposited-energy
    # scoring. Exact slab/box agreement therefore independently exercises the
    # entire-exterior normal, area, and projected-cosine implementation.
    damage = json.loads((ROOT / "examples/stage5_ddd_neq_case.json").read_text(
        encoding="utf-8"
    ))
    damage["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    damage["damage_response"]["path"] = str(
        ROOT / "data/niel/damage_response_library.v1.json"
    )
    damage["accuracy"]["sampling"] = {
        "min_histories": 40, "max_histories": 40, "batch_size": 40
    }
    damage["transport"]["random_seed"] = 1010011
    _, damage_slab_dir = invoke("damage_slab_reference", damage)
    damage_slab = result(damage_slab_dir, "damage_slab_reference")
    damage_csg_case = copy.deepcopy(damage)
    damage_csg_case["geometry"]["absorbers"] = [csg_absorber(box, 2.5)]
    _, damage_csg_dir = invoke("damage_csg_box", damage_csg_case)
    damage_csg = result(damage_csg_dir, "damage_csg_box")
    add_comparison(comparisons, "slab_csg_entry_surface_DDD_equivalence_MeV_g_primary",
                   quantity(damage_slab, "DDD")["value"],
                   quantity(damage_csg, "DDD")["value"],
                   "MeV/g/primary", 2.0e-12)
    add_comparison(comparisons, "slab_csg_entry_surface_neq_equivalence_cm2_primary",
                   quantity(damage_slab, "n_eq")["value"],
                   quantity(damage_csg, "n_eq")["value"],
                   "cm^-2/primary", 2.0e-12)
    require("stage10_geometry_providers_verification_only" in
            quantity(damage_csg, "DDD")["validation_status"],
            "analytic CSG DDD output overclaims Stage-5 validation")

    sphere_radius = 5.0
    sphere_volume = 4.0 * math.pi * sphere_radius ** 3 / 3.0 / 1000.0
    sphere = {"node_id": "detector_sphere", "type": "sphere",
              "inner_radius": {"value": 0.0, "unit": "mm"},
              "outer_radius": {"value": sphere_radius, "unit": "mm"}}
    sphere_case = copy.deepcopy(base)
    sphere_case["geometry"]["absorbers"] = [csg_absorber(
        sphere, sphere_volume, z_mm=-10.0, entry_surface="entire_exterior"
    )]
    _, sphere_dir = invoke("csg_sphere", sphere_case)
    sphere_result = result(sphere_dir, "csg_sphere")
    add_comparison(comparisons, "sphere_analytic_mass_kg",
                   sphere_volume * 2.33 / 1000.0,
                   tid(sphere_result)["absorber_mass"]["value"], "kg", 2.0e-12)

    cylinder_volume = math.pi * 5.0 ** 2 * 2.0 / 1000.0
    cylinder = {"node_id": "detector_cylinder", "type": "cylinder",
                "inner_radius": {"value": 0.0, "unit": "mm"},
                "outer_radius": {"value": 5.0, "unit": "mm"},
                "length": {"value": 2.0, "unit": "mm"}}
    cylinder_case = copy.deepcopy(base)
    cylinder_case["geometry"]["absorbers"] = [csg_absorber(
        cylinder, cylinder_volume, z_mm=-14.0
    )]
    _, cylinder_dir = invoke("csg_cylinder", cylinder_case)
    cylinder_result = result(cylinder_dir, "csg_cylinder")
    add_comparison(comparisons, "cylinder_analytic_mass_kg",
                   cylinder_volume * 2.33 / 1000.0,
                   tid(cylinder_result)["absorber_mass"]["value"], "kg", 2.0e-12)

    # A z-axis rotation preserves this square target and directional source.
    # The geometry hash must change because placement changed, while the exact
    # physical score sequence remains invariant.
    half = math.sqrt(0.5)
    rotated_case = copy.deepcopy(base)
    rotated_case["geometry"]["absorbers"] = [csg_absorber(
        box, 2.5, quaternion=(half, 0.0, 0.0, half)
    )]
    _, rotated_dir = invoke("csg_box_rotated", rotated_case)
    rotated_result = result(rotated_dir, "csg_box_rotated")
    add_comparison(comparisons, "rigid_rotation_tid_invariance_Gy_per_primary",
                   tid(box_result)["value"], tid(rotated_result)["value"],
                   "Gy/primary", 2.0e-12)
    box_hash = resolved(box_dir, "csg_box")["active_absorber_geometry_hash"]
    rotated_hash = resolved(rotated_dir, "csg_box_rotated")["active_absorber_geometry_hash"]
    require(box_hash != rotated_hash, "placement rotation failed to invalidate geometry hash")

    # Exercise a passive aluminum Boolean shell plus a separate silicon sensor.
    shell_case = copy.deepcopy(base)
    shell = csg_absorber(box, 2.5, z_mm=-5.0, absorber_id="enclosed_sensor")
    outer = {"node_id": "outer_box", "type": "box",
             "size_x": {"value": 30.0, "unit": "mm"},
             "size_y": {"value": 30.0, "unit": "mm"},
             "size_z": {"value": 20.0, "unit": "mm"}}
    inner = {"node_id": "inner_void", "type": "box",
             "size_x": {"value": 26.0, "unit": "mm"},
             "size_y": {"value": 26.0, "unit": "mm"},
             "size_z": {"value": 16.0, "unit": "mm"}}
    shell_region = {
        "region_id": "aluminum_enclosure", "role": "passive", "material": "Al",
        "transform": transform(0.0, 0.0, 0.0),
        "solid": {"node_id": "enclosure_shell", "type": "subtraction",
                  "relationship": "right_fully_contained", "left": outer,
                  "right": inner, "right_transform": transform(0.0, 0.0, 0.0)},
        "mass_properties": mass_properties(
            (30.0 * 30.0 * 20.0 - 26.0 * 26.0 * 16.0) / 1000.0,
            "Exact difference of concentric rectangular-prism volumes", 0.05
        ),
    }
    sensor = shell["geometry"]["regions"][0]
    sensor["solid"] = {"node_id": "sensor_box", "type": "box",
                       "size_x": {"value": 5.0, "unit": "mm"},
                       "size_y": {"value": 5.0, "unit": "mm"},
                       "size_z": {"value": 1.0, "unit": "mm"}}
    sensor["mass_properties"] = mass_properties(
        0.025, "Exact rectangular-prism sensor volume"
    )
    shell["geometry"]["regions"] = [shell_region, sensor]
    shell_case["geometry"]["absorbers"] = [shell]
    _, shell_dir = invoke("csg_boolean_enclosure", shell_case)
    shell_result = result(shell_dir, "csg_boolean_enclosure", "enclosed_sensor")
    add_comparison(comparisons, "boolean_enclosure_sensor_mass_kg",
                   0.025 * 2.33 / 1000.0,
                   tid(shell_result)["absorber_mass"]["value"], "kg", 2.0e-12)

    # The separately committed Boolean example is also a user-facing input,
    # rather than only a geometry fixture constructed inside this test. Run it
    # through the normal case path and require that its JSON-selected shield,
    # absorber, and complete calculation list survive into the result. Make
    # reference-library paths absolute because invoke() intentionally copies
    # cases into the isolated test_output run tree.
    documented = json.loads(
        (ROOT / "examples/stage10_analytic_csg_case.json").read_text(encoding="utf-8")
    )
    documented["material_library"]["path"] = str(
        ROOT / "data/material_library.v1.json"
    )
    documented["damage_response"]["path"] = str(
        ROOT / "data/niel/damage_response_library.v1.json"
    )
    require(
        [layer["material"] for layer in documented["geometry"]["shield_stack"]["layers"]]
        == ["Al", "HDPE"],
        "documented Boolean example lost its input-defined shield selection",
    )
    documented_regions = documented["geometry"]["absorbers"][0]["geometry"]["regions"]
    documented_sensitive = [
        region for region in documented_regions if region["role"] == "sensitive"
    ]
    require(
        len(documented_sensitive) == 1 and documented_sensitive[0]["material"] == "Si",
        "documented Boolean example lost its input-defined absorber selection",
    )
    required_functions = ["TID", "DDD", "n_eq", "LET", "H100/10"]
    require(documented["quantities"] == required_functions,
            "documented Boolean example lost its calculation selection")
    _, documented_dir = invoke("documented_boolean_example", documented)
    documented_result = result(
        documented_dir, "documented_boolean_example", "enclosed_silicon_proxy"
    )
    require(
        [entry["name"] for entry in documented_result["outputs"]] == required_functions,
        "documented Boolean result does not match its input calculation selection",
    )

    # Exercise the remaining Boolean constructors with spatially separated
    # passive regions. The intersection reference is an exact overlap of two
    # translated boxes (2 x 4 x 4 mm); its imported centroid is explicit even
    # though general-intersection centroid qualification remains open.
    all_boolean_case = copy.deepcopy(base)
    tiny_sensor = {"node_id": "boolean_case_sensor", "type": "box",
                   "size_x": {"value": 2.0, "unit": "mm"},
                   "size_y": {"value": 2.0, "unit": "mm"},
                   "size_z": {"value": 1.0, "unit": "mm"}}
    all_boolean = csg_absorber(tiny_sensor, 0.004, z_mm=-10.0,
                               quaternion=(half, 0.0, 0.0, half),
                               absorber_id="all_boolean_proxy")
    union_region = {
        "region_id": "disjoint_union", "role": "passive", "material": "Al",
        "transform": transform(-15.0, 0.0, 0.0),
        "solid": {
            "node_id": "union_root", "type": "union", "relationship": "disjoint",
            "left": {"node_id": "union_left", "type": "box",
                     "size_x": {"value": 2.0, "unit": "mm"},
                     "size_y": {"value": 2.0, "unit": "mm"},
                     "size_z": {"value": 2.0, "unit": "mm"}},
            "right": {"node_id": "union_right", "type": "box",
                      "size_x": {"value": 2.0, "unit": "mm"},
                      "size_y": {"value": 2.0, "unit": "mm"},
                      "size_z": {"value": 2.0, "unit": "mm"}},
            "right_transform": transform(6.0, 0.0, 0.0),
        },
        "mass_properties": mass_properties(
            0.016, "Exact sum of two disjoint 2 mm cube volumes",
            centroid_mm=(3.0, 0.0, 0.0)
        ),
    }
    intersection_region = {
        "region_id": "box_intersection", "role": "passive", "material": "Al",
        "transform": transform(15.0, 0.0, 0.0),
        "solid": {
            "node_id": "intersection_root", "type": "intersection",
            "relationship": "general",
            "left": {"node_id": "intersection_left", "type": "box",
                     "size_x": {"value": 4.0, "unit": "mm"},
                     "size_y": {"value": 4.0, "unit": "mm"},
                     "size_z": {"value": 4.0, "unit": "mm"}},
            "right": {"node_id": "intersection_right", "type": "box",
                      "size_x": {"value": 4.0, "unit": "mm"},
                      "size_y": {"value": 4.0, "unit": "mm"},
                      "size_z": {"value": 4.0, "unit": "mm"}},
            "right_transform": transform(2.0, 0.0, 0.0),
        },
        "mass_properties": mass_properties(
            0.032, "Exact overlap of two 4 mm cubes shifted 2 mm along x",
            0.05, centroid_mm=(1.0, 0.0, 0.0)
        ),
    }
    all_boolean["geometry"]["regions"] = [
        union_region, intersection_region, all_boolean["geometry"]["regions"][0]
    ]
    all_boolean_case["geometry"]["absorbers"] = [all_boolean]
    invoke("csg_union_intersection", all_boolean_case)

    # Parser and construction failures exercise ambiguity and silent-repair
    # boundaries. Validation-only cases must fail before Geant4 initialization.
    invalid = copy.deepcopy(box_case)
    invalid["geometry"]["absorbers"][0]["geometry"]["provider_version"] = "2.0.0"
    invoke("reject_provider_version", invalid, validate_only=True, expect_success=False)
    invalid = copy.deepcopy(box_case)
    invalid["geometry"]["absorbers"][0]["geometry"]["coordinate_frame"]["length_unit"] = "unspecified"
    invoke("reject_ambiguous_units", invalid, validate_only=True, expect_success=False)
    invalid = copy.deepcopy(box_case)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["role"] = "passive"
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0].pop("entry_surface")
    invoke("reject_no_sensitive_region", invalid, validate_only=True, expect_success=False)
    invalid = copy.deepcopy(box_case)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["mass_properties"]["volume"]["value"] = 2.0
    invoke("reject_wrong_analytic_volume", invalid, validate_only=True, expect_success=False)
    invalid = copy.deepcopy(box_case)
    invalid["purpose"] = "preliminary_mission_environment_effect_assessment"
    invoke("reject_science_facing_use", invalid, validate_only=True, expect_success=False)
    invalid = copy.deepcopy(all_boolean_case)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["solid"][
        "right_transform"
    ] = transform(1.0, 0.0, 0.0)
    invoke("reject_false_disjoint_union", invalid, validate_only=True,
           expect_success=False)
    overlap = copy.deepcopy(shell_case)
    overlap["geometry"]["absorbers"][0]["geometry"]["regions"][1]["transform"] = \
        transform(0.0, 0.0, -8.0)
    invoke("reject_overlapping_regions", overlap, expect_success=False)

    # Response identity must distinguish providers even for coextensive solids.
    slab_cell = slab_dir / "slab_reference__silicon_proxy_response_cell.json"
    box_cell = box_dir / "csg_box__silicon_proxy_response_cell.json"
    require(slab_cell.is_file() and box_cell.is_file(), "response-cell outputs are missing")
    slab_key = json.loads(slab_cell.read_text(encoding="utf-8"))["response_key"]
    box_key = json.loads(box_cell.read_text(encoding="utf-8"))["response_key"]
    require(slab_key != box_key, "analytic provider failed to invalidate the response key")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "stage": 10,
        "scope": "analytic_csg_v1_verification_only",
        "physical_validation_claim": False,
        "completion_status": "incomplete_pending_independent_nontrivial_shape_reference",
        "comparisons": comparisons,
    }
    evidence_path = EVIDENCE / "stage10_analytic_csg.json"
    evidence_path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    render_figure(comparisons)
    digest = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
    print(f"Stage-10 analytic CSG verification passed; evidence sha256={digest}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL Stage-10 analytic CSG verification: {error}", file=sys.stderr)
        raise SystemExit(1)

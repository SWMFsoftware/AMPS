#!/usr/bin/env python3
"""Generate the directly runnable Stage-10 shape/configuration galleries.

The committed JSON files are the user-facing examples.  This generator is the
single auditable source for their repeated case boilerplate, exact mass
properties, and the deliberately varied source/shield/material selections;
rerunning it must reproduce those files byte for byte.

All mesh vertices are explicit millimetres and all triangular facets use
outward winding.  Provider v1 deliberately accepts convex, closed meshes only;
it is not an STL repair/import utility.
"""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "examples" / "stage10_shapes"
CONFIGURATION_OUTPUT = ROOT / "examples" / "configuration_cases"


def measure(value: float, unit: str = "mm") -> dict:
    return {"value": value, "unit": unit}


def transform(x: float = 0.0, y: float = 0.0, z: float = 0.0) -> dict:
    return {
        "translation": {axis: measure(value) for axis, value in zip("xyz", (x, y, z))},
        "rotation_quaternion_local_to_world": {
            "w": 1.0, "x": 0.0, "y": 0.0, "z": 0.0
        },
    }


def mass_properties(volume_cm3: float, centroid_mm=(0.0, 0.0, 0.0),
                    provenance: str = "Exact closed-form geometry") -> dict:
    return {
        "volume": measure(volume_cm3, "cm3"),
        "centroid": {
            axis: measure(value) for axis, value in zip("xyz", centroid_mm)
        },
        "relative_tolerance": 1.0e-9,
        "provenance": provenance,
    }


def common_geometry(provider_type: str, provider_z_mm: float, source_id: str,
                    regions: list[dict]) -> dict:
    return {
        "type": provider_type,
        "provider_version": "1.0.0",
        "coordinate_frame": {
            "name": "absorber_local_cartesian",
            "handedness": "right",
            "origin": "provider_origin",
            "length_unit": "mm",
        },
        "transform": transform(z=provider_z_mm),
        "source": {
            "kind": "inline",
            "source_id": source_id,
            "source_version": "1.0.0",
            "provenance": (
                "Synthetic Stage-10 shape example with exact dimensions; "
                "not an as-built spacecraft model"
            ),
        },
        "regions": regions,
    }


def analytic_region(solid: dict, volume_cm3: float, entry_surface: str) -> dict:
    return {
        "region_id": "detector",
        "role": "sensitive",
        "material": "Si",
        "transform": transform(),
        "solid": solid,
        "entry_surface": entry_surface,
        "mass_properties": mass_properties(volume_cm3),
    }


def mesh_region(vertices: list[list[float]], facets: list[tuple[str, int, int, int]],
                volume_cm3: float, centroid_mm: tuple[float, float, float],
                entry_surface: dict, provenance: str) -> dict:
    return {
        "region_id": "detector",
        "role": "sensitive",
        "material": "Si",
        "transform": transform(),
        "mesh": {
            "length_unit": "mm",
            "validation_tolerance": measure(1.0e-9),
            "vertices": vertices,
            "facets": [
                {"facet_id": facet_id, "vertex_indices": [a, b, c]}
                for facet_id, a, b, c in facets
            ],
        },
        "entry_surface": entry_surface,
        "mass_properties": mass_properties(volume_cm3, centroid_mm, provenance),
    }


def case_document(name: str, absorber: dict) -> dict:
    # Start from the Stage-6 two-node proton case because its 50/150 MeV mix
    # gives a nonzero, finite H100/10 example and already contains the explicit
    # LET binning contract. Add the registered silicon damage-response library
    # below so all five screenshot quantities can be selected in this same
    # ordinary case input. The shape, shield, absorber material, and quantities
    # are not command-line defaults: each remains visible and editable in the
    # generated JSON file.
    case = json.loads((ROOT / "examples" / "stage6_let_hardness_case.json").read_text())
    case["case_id"] = f"stage10_shape_example_{name}"
    # These cases live one directory below the other examples.
    case["material_library"]["path"] = "../../data/material_library.v1.json"
    case["geometry"]["absorbers"] = [absorber]
    case["quantities"] = ["TID", "DDD", "n_eq", "LET", "H100/10"]
    case["damage_response"] = {
        "path": "../../data/niel/damage_response_library.v1.json",
        "library_id": "shieldSim.silicon_displacement_damage",
        "library_version": "1.0.0",
    }
    case["accuracy"]["sampling"] = {
        "min_histories": 20, "max_histories": 20, "batch_size": 20
    }
    # The gallery is a short verification demonstration. A complete
    # application-specific systematic budget is intentionally not invented.
    case["accuracy"]["systematic_budget"] = {"status": "incomplete"}
    case["transport"]["random_seed"] = 1010100
    case["output"]["prefix"] = f"stage10_{name}"
    return case


def source_node(incident_id: str, particle_name: str, atomic_number: int,
                mass_number: int, energy_value: float, energy_unit: str = "MeV",
                energy_coordinate: str = "total_kinetic_energy",
                weight: float = 1.0) -> dict:
    """Return one exact verification source node with no hidden defaults."""
    return {
        "incident_id": incident_id,
        "particle": {
            "name": particle_name, "Z": atomic_number, "A": mass_number,
        },
        "energy": {
            "coordinate": energy_coordinate,
            "value": energy_value,
            "unit": energy_unit,
        },
        "direction": {"mu": 1.0, "phi": measure(0.0, "rad")},
        "weight": weight,
    }


def set_source(case: dict, nodes: list[dict], size_x_mm: float,
               size_y_mm: float, center_x_mm: float = 0.0,
               center_y_mm: float = 0.0) -> None:
    """Install a normalized discrete source and its explicit plane footprint."""
    if abs(sum(node["weight"] for node in nodes) - 1.0) > 1.0e-12:
        raise ValueError("configuration-example source weights must sum to one")
    case["source"]["spectrum"] = {
        "format": "discrete_response_nodes", "nodes": nodes,
    }
    plane = case["geometry"]["source_plane"]
    plane["size_x"] = measure(size_x_mm)
    plane["size_y"] = measure(size_y_mm)
    plane["transform"] = transform(center_x_mm, center_y_mm, -80.0)


def set_shield(case: dict, size_x_mm: float, size_y_mm: float,
               layers: list[dict]) -> None:
    """Install the ordered upstream-to-downstream shield definition."""
    if not 1 <= len(layers) <= 5:
        raise ValueError("configuration examples require one through five layers")
    stack = case["geometry"]["shield_stack"]
    stack["size_x"] = measure(size_x_mm)
    stack["size_y"] = measure(size_y_mm)
    stack["layers"] = layers


def geometric_layer(layer_id: str, material: str, value: float,
                    unit: str = "mm") -> dict:
    return {
        "layer_id": layer_id,
        "material": material,
        "thickness": {"type": "geometric", "value": value, "unit": unit},
    }


def areal_layer(layer_id: str, material: str, value: float) -> dict:
    return {
        "layer_id": layer_id,
        "material": material,
        "thickness": {
            "type": "areal_density", "value": value, "unit": "g/cm2",
        },
    }


def select_quantities(case: dict, quantities: list[str]) -> None:
    """Keep mandatory response blocks synchronized with quantity selection."""
    case["quantities"] = quantities
    if "DDD" in quantities or "n_eq" in quantities:
        case["damage_response"] = {
            "path": "../../data/niel/damage_response_library.v1.json",
            "library_id": "shieldSim.silicon_displacement_damage",
            "library_version": "1.0.0",
        }
    else:
        case.pop("damage_response", None)
    if "LET" not in quantities:
        case.pop("let_response", None)


def set_sensitive_material(case: dict, material: str) -> None:
    """Change the one scored region without altering passive regions."""
    regions = case["geometry"]["absorbers"][0]["geometry"]["regions"]
    sensitive = [region for region in regions if region["role"] == "sensitive"]
    if len(sensitive) != 1:
        raise ValueError("configuration examples require one sensitive region")
    sensitive[0]["material"] = material


def configuration_case(name: str, absorber: dict) -> dict:
    """Create short-run boilerplate shared by the selection examples."""
    case = case_document(name, copy.deepcopy(absorber))
    case["case_id"] = f"configuration_example_{name}"
    case["output"]["prefix"] = f"configuration_{name}"
    # These are executable configuration demonstrations, not converged science
    # calculations.  The 10% request remains visible while the fixed 20-history
    # bounds keep regression execution short and intentionally unconverged.
    case["accuracy"]["sampling"] = {
        "min_histories": 20, "max_histories": 20, "batch_size": 20,
    }
    case["accuracy"]["systematic_budget"] = {"status": "incomplete"}
    return case


def analytic_examples() -> dict[str, dict]:
    box = {
        "node_id": "detector_box", "type": "box",
        "size_x": measure(20.0), "size_y": measure(12.0), "size_z": measure(2.0),
    }
    sphere_radius = 5.0
    sphere = {
        "node_id": "detector_sphere", "type": "sphere",
        "inner_radius": measure(0.0), "outer_radius": measure(sphere_radius),
    }
    cylinder_radius, cylinder_length = 5.0, 4.0
    cylinder = {
        "node_id": "detector_cylinder", "type": "cylinder",
        "inner_radius": measure(0.0), "outer_radius": measure(cylinder_radius),
        "length": measure(cylinder_length),
    }
    return {
        "analytic_box": {
            "absorber_id": "analytic_box_proxy",
            "geometry": common_geometry(
                "analytic_csg", -14.0, "analytic_box_example",
                [analytic_region(box, 20.0 * 12.0 * 2.0 / 1000.0,
                                 "local_negative_z_face")],
            ),
        },
        "analytic_sphere": {
            "absorber_id": "analytic_sphere_proxy",
            "geometry": common_geometry(
                "analytic_csg", -10.0, "analytic_sphere_example",
                [analytic_region(sphere, 4.0 * math.pi * sphere_radius**3 / 3.0 / 1000.0,
                                 "entire_exterior")],
            ),
        },
        "analytic_cylinder": {
            "absorber_id": "analytic_cylinder_proxy",
            "geometry": common_geometry(
                "analytic_csg", -12.0, "analytic_cylinder_example",
                [analytic_region(cylinder,
                                 math.pi * cylinder_radius**2 * cylinder_length / 1000.0,
                                 "local_negative_z_face")],
            ),
        },
    }


def mesh_examples() -> dict[str, dict]:
    # A triangulated rectangular solid. The two z-negative facets together are
    # the same 20 x 12 mm planar entry patch as the analytic-box example.
    box_vertices = [
        [-10, -6, -1], [10, -6, -1], [10, 6, -1], [-10, 6, -1],
        [-10, -6, 1], [10, -6, 1], [10, 6, 1], [-10, 6, 1],
    ]
    box_facets = [
        ("zneg_a", 0, 2, 1), ("zneg_b", 0, 3, 2),
        ("zpos_a", 4, 5, 6), ("zpos_b", 4, 6, 7),
        ("yneg_a", 0, 1, 5), ("yneg_b", 0, 5, 4),
        ("xpos_a", 1, 2, 6), ("xpos_b", 1, 6, 5),
        ("ypos_a", 2, 3, 7), ("ypos_b", 2, 7, 6),
        ("xneg_a", 3, 0, 4), ("xneg_b", 3, 4, 7),
    ]

    # Right tetrahedron: x,y,z >= 0 and x+y+z <= 10 mm. The order below gives
    # outward normals; volume is abc/6 and centroid is (a,b,c)/4.
    tetra_vertices = [[0, 0, 0], [10, 0, 0], [0, 10, 0], [0, 0, 10]]
    tetra_facets = [
        ("z0_entry", 0, 2, 1), ("y0", 0, 1, 3),
        ("x0", 0, 3, 2), ("slanted", 1, 2, 3),
    ]

    # Regular axis-aligned octahedron |x|+|y|+|z| <= r. It is a useful curved-
    # body mesh surrogate while retaining an exact 4*r^3/3 volume.
    radius = 5.0
    octa_vertices = [
        [radius, 0, 0], [-radius, 0, 0], [0, radius, 0],
        [0, -radius, 0], [0, 0, radius], [0, 0, -radius],
    ]
    octa_facets = [
        ("top_pp", 0, 2, 4), ("top_np", 2, 1, 4),
        ("top_nn", 1, 3, 4), ("top_pn", 3, 0, 4),
        ("bottom_pp", 2, 0, 5), ("bottom_np", 1, 2, 5),
        ("bottom_nn", 3, 1, 5), ("bottom_pn", 0, 3, 5),
    ]

    return {
        "mesh_box": {
            "absorber_id": "mesh_box_proxy",
            "geometry": common_geometry(
                "tessellated_mesh", -14.0, "triangulated_box_example",
                [mesh_region(
                    box_vertices, box_facets, 0.48, (0.0, 0.0, 0.0),
                    {"type": "tagged_facets", "facet_ids": ["zneg_a", "zneg_b"]},
                    "Exact rectangular-prism volume and symmetry centroid",
                )],
            ),
        },
        "mesh_tetrahedron": {
            "absorber_id": "mesh_tetrahedron_proxy",
            "geometry": common_geometry(
                "tessellated_mesh", -12.0, "right_tetrahedron_example",
                [mesh_region(
                    tetra_vertices, tetra_facets, 1.0 / 6.0, (2.5, 2.5, 2.5),
                    {"type": "tagged_facets", "facet_ids": ["z0_entry"]},
                    "Exact right-tetrahedron abc/6 volume and vertex-average centroid",
                )],
            ),
        },
        "mesh_octahedron": {
            "absorber_id": "mesh_octahedron_proxy",
            "geometry": common_geometry(
                "tessellated_mesh", -10.0, "regular_octahedron_example",
                [mesh_region(
                    octa_vertices, octa_facets, 4.0 * radius**3 / 3.0 / 1000.0,
                    (0.0, 0.0, 0.0), {"type": "entire_exterior"},
                    "Exact |x|+|y|+|z|<=r octahedron volume and symmetry centroid",
                )],
            ),
        },
    }


def configuration_examples(absorbers: dict[str, dict]) -> dict[str, dict]:
    """Build runnable examples spanning the supported selection boundaries.

    No single source can make all five quantities nonzero: Stage-4 TID accepts
    proton/He-4 incidents, while Stage-5 n_eq is explicitly neutron-only.  The
    first case therefore shows the complete five-name input selection and an
    intentional zero n_eq; the neutron cases provide the complementary nonzero
    n_eq calculation.  Keeping those semantics in separate files prevents a
    synthetic mixed field from being mistaken for an environmental spectrum.
    """
    examples: dict[str, dict] = {}

    # Natural-silicon analytic box: one case selects every implemented output.
    # Proton nodes exercise TID, proton DDD, LET, and the two-threshold H100/10
    # ratio.  n_eq is present but exactly zero by its neutron-only definition.
    name = "01_si_box_proton_all_quantities"
    case = configuration_case(name, absorbers["analytic_box"])
    set_source(case, [
        source_node("proton_50MeV", "proton", 1, 1, 50.0, weight=0.5),
        source_node("proton_150MeV", "proton", 1, 1, 150.0, weight=0.5),
    ], 18.0, 10.0)
    set_shield(case, 60.0, 40.0, [
        geometric_layer("outer_al", "Al", 1.5),
        areal_layer("middle_hdpe", "HDPE", 0.19),
        geometric_layer("inner_ta", "Ta", 0.1),
    ])
    select_quantities(case, ["TID", "DDD", "n_eq", "LET", "H100/10"])
    case["transport"]["random_seed"] = 1020101
    examples[name] = case

    # A monoenergetic neutron in the ASTM E722 1.0--1.1 MeV group makes both
    # silicon DDD and n_eq directly inspectable.  The small source plane lies
    # inside the cylinder projection so this 20-history example records entry
    # crossings instead of relying on a rare geometric hit.
    name = "02_si_cylinder_neutron_damage"
    case = configuration_case(name, absorbers["analytic_cylinder"])
    set_source(case, [
        source_node("neutron_1p05MeV", "neutron", 0, 1, 1.05),
    ], 6.0, 6.0)
    set_shield(case, 45.0, 45.0, [
        geometric_layer("outer_al", "Al", 0.05),
    ])
    select_quantities(case, ["DDD", "n_eq"])
    case["transport"].update({
        "physics_list": "FTFP_BERT_HP",
        "production_cut": measure(0.1),
        "max_step": measure(0.1),
        "random_seed": 1020102,
    })
    examples[name] = case

    # A five-layer maximum-count example with He-4 demonstrates that source
    # species, layer count/order, and absorber shape remain independent input
    # selections.  Thin layers preserve a transmitted alpha component in this
    # quick verification fixture; they are not a recommended flight stack.
    name = "03_si_sphere_alpha_tid_ddd_let"
    case = configuration_case(name, absorbers["analytic_sphere"])
    set_source(case, [
        source_node("alpha_100MeV_total", "alpha", 2, 4, 100.0),
    ], 6.0, 6.0)
    set_shield(case, 50.0, 50.0, [
        geometric_layer("outer_kapton", "Kapton", 0.05),
        geometric_layer("al_structure", "Al", 0.10),
        geometric_layer("hdpe_liner", "HDPE", 0.10),
        geometric_layer("ta_spot_shield", "Ta", 0.01),
        geometric_layer("water_layer", "Water", 0.05),
    ])
    select_quantities(case, ["TID", "DDD", "LET"])
    case["transport"]["random_seed"] = 1020103
    examples[name] = case

    # Water is inside the qualified TID/LET absorber domain.  The cylinder's
    # exact dimensions remain in the analytic provider while H100/10 is scored
    # independently at the rear face of the explicitly sized Cu/HDPE shield.
    name = "04_water_cylinder_proton_tid_let_h100"
    case = configuration_case(name, absorbers["analytic_cylinder"])
    set_sensitive_material(case, "Water")
    set_source(case, [
        source_node("proton_50MeV", "proton", 1, 1, 50.0, weight=0.5),
        source_node("proton_150MeV", "proton", 1, 1, 150.0, weight=0.5),
    ], 6.0, 6.0)
    set_shield(case, 40.0, 30.0, [
        geometric_layer("outer_cu", "Cu", 0.5),
        geometric_layer("inner_hdpe", "HDPE", 4.0),
    ])
    select_quantities(case, ["TID", "LET", "H100/10"])
    case["transport"]["random_seed"] = 1020104
    examples[name] = case

    # The same proton calculation is legal for an aluminum absorber and a
    # triangulated box.  Four shielding materials and a non-square footprint
    # make the editable transverse dimensions and physical layer order clear.
    name = "05_al_mesh_box_proton_tid_let_h100"
    case = configuration_case(name, absorbers["mesh_box"])
    set_sensitive_material(case, "Al")
    set_source(case, [
        source_node("proton_50MeV", "proton", 1, 1, 50.0, weight=0.5),
        source_node("proton_150MeV", "proton", 1, 1, 150.0, weight=0.5),
    ], 18.0, 10.0)
    set_shield(case, 70.0, 45.0, [
        geometric_layer("outer_kapton", "Kapton", 0.10),
        geometric_layer("al_structure", "Al", 0.75),
        geometric_layer("kevlar_layer", "Kevlar", 0.50),
        areal_layer("inner_hdpe", "HDPE", 0.095),
    ])
    select_quantities(case, ["TID", "LET", "H100/10"])
    case["transport"]["random_seed"] = 1020105
    examples[name] = case

    # A tagged triangular entry patch demonstrates neutron damage scoring on a
    # non-box mesh.  Centering a 4x4 mm source at (2,2) mm keeps all unscattered
    # rays inside x>=0, y>=0, x+y<=10 on the tetrahedron's tagged z=0 face.
    name = "06_si_mesh_tetrahedron_neutron_damage"
    case = configuration_case(name, absorbers["mesh_tetrahedron"])
    set_source(case, [
        source_node("neutron_1p05MeV", "neutron", 0, 1, 1.05),
    ], 4.0, 4.0, 2.0, 2.0)
    set_shield(case, 55.0, 35.0, [
        geometric_layer("outer_al", "Al", 0.05),
        geometric_layer("inner_hdpe", "HDPE", 0.05),
    ])
    select_quantities(case, ["DDD", "n_eq"])
    case["transport"].update({
        "physics_list": "FTFP_BERT_HP",
        "production_cut": measure(0.1),
        "max_step": measure(0.1),
        "random_seed": 1020106,
    })
    examples[name] = case

    # Case-local properties can be selected without editing the shared
    # material library.  The deliberately synthetic density value is labeled
    # as such and affects only this verification input; a science case must
    # replace it and its provenance with project measurements.
    name = "07_si_mesh_octahedron_case_material_override"
    case = configuration_case(name, absorbers["mesh_octahedron"])
    case["material_overrides"] = [{
        "material_id": "Al_project_density_example",
        "revision": "verification-1",
        "extends": "Al_natural_G4_11_4@1",
        "status": "project_specific",
        "density": {
            "value": 2.70,
            "unit": "g/cm3",
            "uncertainty": {
                "type": "standard", "value": 0.01, "unit": "g/cm3",
            },
            "reference_condition": "synthetic configuration example",
        },
        "provenance": [{
            "properties": ["density"],
            "source": "synthetic verification value; replace for science use",
            "version": "configuration-example-1",
        }],
        "applicability_and_limitations": (
            "Input-schema demonstration only; not a flight alloy or lot."
        ),
    }]
    set_source(case, [
        source_node("proton_50MeV", "proton", 1, 1, 50.0, weight=0.5),
        source_node("proton_150MeV", "proton", 1, 1, 150.0, weight=0.5),
    ], 6.0, 6.0)
    set_shield(case, 65.0, 65.0, [
        areal_layer("project_al", "Al_project_density_example", 0.27),
        geometric_layer("inner_hdpe", "HDPE", 2.0),
    ])
    select_quantities(case, ["TID", "DDD", "n_eq", "LET", "H100/10"])
    case["transport"]["random_seed"] = 1020107
    examples[name] = case

    return examples


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    absorbers = {**analytic_examples(), **mesh_examples()}
    for name, absorber in sorted(absorbers.items()):
        document = case_document(name, copy.deepcopy(absorber))
        (OUTPUT / f"{name}.json").write_text(
            json.dumps(document, indent=2, sort_keys=False) + "\n", encoding="utf-8"
        )
    configurations = configuration_examples(absorbers)
    CONFIGURATION_OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, document in configurations.items():
        (CONFIGURATION_OUTPUT / f"{name}.json").write_text(
            json.dumps(document, indent=2, sort_keys=False) + "\n", encoding="utf-8"
        )
    print(f"wrote {len(absorbers)} Stage-10 shape examples to {OUTPUT}")
    print(
        f"wrote {len(configurations)} selection examples to "
        f"{CONFIGURATION_OUTPUT}"
    )


if __name__ == "__main__":
    main()

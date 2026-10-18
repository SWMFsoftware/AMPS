#!/usr/bin/env python3
"""Stage-10 tessellated-mesh geometry, transport, and rejection gate.

Runs every committed analytic/mesh gallery case, compares exact polyhedral
mass properties with parser and Geant4 results, proves a triangulated box is
transport-equivalent to its analytic counterpart, and tests the strict mesh
boundary.  The gate is deliberately verification-only: it does not qualify
general concave CAD, STL import, or mission-facing arbitrary geometry.
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
        "usage: test_stage10_tessellated_mesh.py SHIELDSIM SOURCE_ROOT "
        "EVIDENCE_DIR RUN_DIR FIGURE_DIR"
    )

EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
FIGURES = pathlib.Path(sys.argv[5]).resolve()
EXAMPLES = ROOT / "examples" / "stage10_shapes"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def invoke(name: str, document: dict, validate_only: bool = False,
           expect_success: bool = True) -> tuple[subprocess.CompletedProcess[str], pathlib.Path]:
    directory = RUN_DIR / name
    directory.mkdir(parents=True, exist_ok=True)
    case = copy.deepcopy(document)
    case["case_id"] = f"stage10_mesh_test_{name}"
    case["material_library"]["path"] = str(ROOT / "data" / "material_library.v1.json")
    case["output"]["prefix"] = name
    path = directory / "case.json"
    path.write_text(json.dumps(case, indent=2) + "\n", encoding="utf-8")
    option = f"--validate-case={path}" if validate_only else f"--case={path}"
    completed = subprocess.run(
        [str(EXE), option], cwd=directory, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    (directory / "run.log").write_text(completed.stdout, encoding="utf-8")
    require((completed.returncode == 0) == expect_success,
            f"{name}: expected success={expect_success}; see {directory / 'run.log'}")
    return completed, directory


def load_example(name: str) -> dict:
    return json.loads((EXAMPLES / f"{name}.json").read_text(encoding="utf-8"))


def result(directory: pathlib.Path, prefix: str, absorber_id: str) -> dict:
    path = directory / f"{prefix}__{absorber_id}_result.json"
    require(path.is_file(), f"missing result {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def resolved(directory: pathlib.Path, prefix: str, absorber_id: str) -> dict:
    path = directory / f"{prefix}__{absorber_id}_resolved_case.json"
    require(path.is_file(), f"missing resolved case {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def tid(payload: dict) -> dict:
    selected = [item for item in payload["outputs"] if item["name"] == "TID"]
    require(len(selected) == 1, "expected exactly one TID quantity")
    return selected[0]


def quantity(payload: dict, name: str) -> dict:
    selected = [item for item in payload["outputs"] if item["name"] == name]
    require(len(selected) == 1, f"expected exactly one {name} quantity")
    return selected[0]


def add_comparison(comparisons: list[dict], comparison_id: str, expected: float,
                   actual: float, unit: str, tolerance: float) -> None:
    relative = abs(actual - expected) / max(abs(expected), 1.0e-300)
    require(relative <= tolerance,
            f"{comparison_id}: relative difference {relative} exceeds {tolerance}")
    comparisons.append({
        "id": comparison_id,
        "expected": expected,
        "actual": actual,
        "unit": unit,
        "relative_tolerance": tolerance,
        "passed": True,
    })


def canonical_mesh_properties(manifest: dict) -> dict:
    absorber = manifest["geometry"]["absorbers"][0]
    require(absorber["geometry_provider"] == "tessellated_mesh.v1",
            "resolved manifest lost mesh provider identity")
    regions = absorber["canonical_provider_geometry"]["regions"]
    sensitive = [region for region in regions if region["role"] == "sensitive"]
    require(len(sensitive) == 1, "expected one canonical sensitive mesh region")
    return sensitive[0]["mass_properties"]


def render_figure(comparisons: list[dict]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise AssertionError(f"Stage-10 mesh comparison figure requires matplotlib: {error}")

    mass = [item for item in comparisons
            if item["id"].startswith("mesh_") and item["id"].endswith("_exact_mass_kg")]
    transport = [item for item in comparisons if "equivalence" in item["id"]]
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.8))

    labels = [item["id"].replace("_exact_mass_kg", "") for item in mass]
    positions = list(range(len(mass)))
    axes[0].bar([x - 0.18 for x in positions], [x["expected"] for x in mass],
                width=0.36, label="analytic polyhedron")
    axes[0].bar([x + 0.18 for x in positions], [x["actual"] for x in mass],
                width=0.36, label="Geant4 tessellation")
    axes[0].set_xticks(positions)
    axes[0].set_xticklabels(labels, rotation=18, ha="right")
    axes[0].set_ylabel("sensitive mass (kg)")
    axes[0].set_title("Exact mass-property checks")
    axes[0].legend(fontsize=8)
    axes[0].grid(axis="y", alpha=0.25)

    tlabels = [item["id"].replace("analytic_mesh_box_", "").replace("_equivalence", "")
               for item in transport]
    normalized_expected, normalized_actual = [], []
    for item in transport:
        scale = max(abs(item["expected"]), abs(item["actual"]), 1.0e-300)
        normalized_expected.append(item["expected"] / scale)
        normalized_actual.append(item["actual"] / scale)
    tpositions = list(range(len(transport)))
    axes[1].bar([x - 0.18 for x in tpositions], normalized_expected,
                width=0.36, label="analytic box")
    axes[1].bar([x + 0.18 for x in tpositions], normalized_actual,
                width=0.36, label="triangulated box")
    axes[1].set_xticks(tpositions)
    axes[1].set_xticklabels(tlabels, rotation=18, ha="right")
    axes[1].set_ylabel("value normalized per comparison")
    axes[1].set_title("Representation-equivalent transport")
    axes[1].legend(fontsize=8)
    axes[1].grid(axis="y", alpha=0.25)

    figure.suptitle("Stage 10 tessellated_mesh.v1 verification")
    figure.text(
        0.01, 0.01,
        "Closed outward convex meshes only; verification examples do not qualify arbitrary "
        "concave CAD/STL or mission-facing geometry.",
        fontsize=8,
    )
    figure.tight_layout(rect=(0.0, 0.05, 1.0, 0.95))
    for extension in ("png", "svg"):
        path = FIGURES / f"stage10_tessellated_mesh_comparisons.{extension}"
        figure.savefig(path, dpi=170 if extension == "png" else None,
                       bbox_inches="tight")
        require(path.stat().st_size > 1000, f"unexpectedly small figure {path}")
    plt.close(figure)


def main() -> int:
    require(EXE.is_file(), f"missing shieldSim executable {EXE}")
    expected_examples = {
        "analytic_box", "analytic_cylinder", "analytic_sphere",
        "mesh_box", "mesh_tetrahedron", "mesh_octahedron",
    }
    require({path.stem for path in EXAMPLES.glob("*.json")} == expected_examples,
            "the committed Stage-10 shape gallery is incomplete or has an untested case")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    comparisons: list[dict] = []

    runs: dict[str, tuple[dict, pathlib.Path, dict]] = {}
    for name in sorted(expected_examples):
        document = load_example(name)
        absorber_id = document["geometry"]["absorbers"][0]["absorber_id"]
        _, directory = invoke(name, document)
        payload = result(directory, name, absorber_id)
        manifest = resolved(directory, name, absorber_id)
        runs[name] = (payload, directory, manifest)
        require({item["name"] for item in payload["outputs"]} ==
                {"TID", "DDD", "n_eq", "LET", "H100/10"},
                f"{name} did not execute every input-selected calculation")
        require("stage10_geometry_providers_verification_only" in
                tid(payload)["validation_status"],
                f"{name} overclaims a completed Stage-10 validation")

    # Every shape case exposes the three screenshot selections through the
    # strict case schema. Validate a second configuration for every shape:
    # copper replaces the first shield layer, water replaces silicon, and the
    # calculation list drops the silicon-only DDD/n_eq functions. This proves
    # these choices are read from each JSON file rather than fixed in C++ or in
    # the example generator. Build and transport one analytic and one mesh box
    # variant as the end-to-end representative pair.
    selected_runs: dict[str, dict] = {}
    for name in sorted(expected_examples):
        selected = load_example(name)
        layers = selected["geometry"]["shield_stack"]["layers"]
        require(1 <= len(layers) <= 5, f"{name} has no input-defined shield selection")
        layers[0]["material"] = "Cu"
        geometry = selected["geometry"]["absorbers"][0]["geometry"]
        sensitive = [region for region in geometry["regions"]
                     if region["role"] == "sensitive"]
        require(len(sensitive) == 1, f"{name} has no unique input-defined absorber selection")
        sensitive[0]["material"] = "Water"
        selected["quantities"] = ["TID", "LET", "H100/10"]
        selected.pop("damage_response")
        run_transport = name in {"analytic_box", "mesh_box"}
        _, selected_directory = invoke(
            f"selected_{name}", selected, validate_only=not run_transport
        )
        if run_transport:
            absorber_id = selected["geometry"]["absorbers"][0]["absorber_id"]
            selected_payload = result(selected_directory, f"selected_{name}", absorber_id)
            require({item["name"] for item in selected_payload["outputs"]} ==
                    {"TID", "LET", "H100/10"},
                    f"{name} ignored its modified calculation selection")
            selected_runs[name] = selected_payload

    selected_analytic_tid = tid(selected_runs["analytic_box"])
    selected_mesh_tid = tid(selected_runs["mesh_box"])
    add_comparison(comparisons, "selected_water_analytic_mass_kg",
                   0.48 / 1000.0,
                   selected_analytic_tid["absorber_mass"]["value"], "kg", 2.0e-9)
    add_comparison(comparisons, "selected_water_mesh_mass_kg",
                   0.48 / 1000.0,
                   selected_mesh_tid["absorber_mass"]["value"], "kg", 2.0e-9)
    add_comparison(comparisons, "selected_water_analytic_mesh_TID_equivalence",
                   selected_analytic_tid["value"], selected_mesh_tid["value"],
                   "Gy/primary", 1.0e-8)

    # The parser independently derives oriented-polyhedron mass properties;
    # Geant4 then independently reports the constructed sensitive mass.
    volumes = {
        "mesh_box": 20.0 * 12.0 * 2.0 / 1000.0,
        "mesh_tetrahedron": 10.0**3 / 6.0 / 1000.0,
        "mesh_octahedron": 4.0 * 5.0**3 / 3.0 / 1000.0,
    }
    areas = {
        "mesh_box": 2.0 * (20.0 * 12.0 + 20.0 * 2.0 + 12.0 * 2.0) / 100.0,
        "mesh_tetrahedron": (150.0 + 50.0 * math.sqrt(3.0)) / 100.0,
        "mesh_octahedron": 4.0 * math.sqrt(3.0) * 5.0**2 / 100.0,
    }
    centroids = {
        "mesh_box": (0.0, 0.0, 0.0),
        "mesh_tetrahedron": (2.5, 2.5, 2.5),
        "mesh_octahedron": (0.0, 0.0, 0.0),
    }
    density = 2.33
    for name in volumes:
        payload, _, manifest = runs[name]
        props = canonical_mesh_properties(manifest)
        add_comparison(comparisons, f"{name}_exact_volume_cm3", volumes[name],
                       props["exact_volume_cm3"], "cm3", 2.0e-12)
        add_comparison(comparisons, f"{name}_exact_surface_area_cm2", areas[name],
                       props["exact_surface_area_cm2"], "cm2", 2.0e-12)
        for axis, expected, actual in zip("xyz", centroids[name], props["exact_centroid_mm"]):
            absolute = abs(actual - expected)
            require(absolute <= 1.0e-10,
                    f"{name} centroid {axis} error {absolute} mm exceeds tolerance")
            comparisons.append({
                "id": f"{name}_exact_centroid_{axis}_mm",
                "expected": expected, "actual": actual, "unit": "mm",
                "absolute_tolerance": 1.0e-10, "passed": True,
            })
        add_comparison(comparisons, f"{name}_exact_mass_kg",
                       volumes[name] * density / 1000.0,
                       tid(payload)["absorber_mass"]["value"], "kg", 2.0e-9)

    analytic_tid = tid(runs["analytic_box"][0])
    mesh_tid = tid(runs["mesh_box"][0])
    add_comparison(comparisons, "analytic_mesh_box_mass_kg_equivalence",
                   analytic_tid["absorber_mass"]["value"],
                   mesh_tid["absorber_mass"]["value"], "kg", 2.0e-9)
    add_comparison(comparisons, "analytic_mesh_box_TID_Gy_primary_equivalence",
                   analytic_tid["value"], mesh_tid["value"], "Gy/primary", 1.0e-8)

    # TID exercises sensitive mass and energy deposition but does not consume
    # the provider entry-surface estimator. The same two gallery inputs also
    # request DDD/n_eq, so these comparisons execute both tagged mesh triangles,
    # their summed area, outward normal, point-in-triangle test, and cosine
    # weight without changing the source, shield, or absorber selection.
    add_comparison(comparisons, "analytic_mesh_box_DDD_equivalence",
                   quantity(runs["analytic_box"][0], "DDD")["value"],
                   quantity(runs["mesh_box"][0], "DDD")["value"],
                   "MeV/g/primary", 1.0e-8)
    add_comparison(comparisons, "analytic_mesh_box_n_eq_equivalence",
                   quantity(runs["analytic_box"][0], "n_eq")["value"],
                   quantity(runs["mesh_box"][0], "n_eq")["value"],
                   "cm^-2/primary", 1.0e-8)

    analytic_cell = runs["analytic_box"][1] / \
        "analytic_box__analytic_box_proxy_response_cell.json"
    mesh_cell = runs["mesh_box"][1] / "mesh_box__mesh_box_proxy_response_cell.json"
    require(analytic_cell.is_file() and mesh_cell.is_file(), "response cells are missing")
    require(json.loads(analytic_cell.read_text())["response_key"] !=
            json.loads(mesh_cell.read_text())["response_key"],
            "analytic and mesh providers collided in the response cache key")

    # Fatal input boundaries: provider v1 never guesses units, closes holes,
    # flips faces, repairs concavity, or accepts undocumented entry patches.
    valid = load_example("mesh_box")
    mesh = valid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["mesh"]
    invalid = copy.deepcopy(valid)
    invalid["geometry"]["absorbers"][0]["geometry"]["provider_version"] = "2.0.0"
    invoke("reject_mesh_provider_version", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["mesh"]["length_unit"] = "m"
    invoke("reject_mesh_ambiguous_units", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid_mesh = invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["mesh"]
    invalid_mesh["facets"].pop()
    invoke("reject_open_mesh", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid_facet = invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["mesh"]["facets"][0]
    invalid_facet["vertex_indices"][1], invalid_facet["vertex_indices"][2] = \
        invalid_facet["vertex_indices"][2], invalid_facet["vertex_indices"][1]
    invoke("reject_reversed_facet", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["mesh"]["vertices"][6] = [0, 0, 0]
    invoke("reject_nonconvex_mesh", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["entry_surface"]["facet_ids"] = ["missing"]
    invoke("reject_unknown_entry_facet", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["entry_surface"]["facet_ids"] = ["zneg_a", "zpos_a"]
    invoke("reject_nonplanar_entry_patch", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid["geometry"]["absorbers"][0]["geometry"]["regions"][0]["mass_properties"]["volume"]["value"] = 0.40
    invoke("reject_wrong_mesh_volume", invalid, True, False)
    invalid = copy.deepcopy(valid)
    invalid["purpose"] = "preliminary_mission_environment_effect_assessment"
    invoke("reject_mesh_science_facing_use", invalid, True, False)
    # Keep this local variable used: it documents that the original fixture
    # remained unchanged across deep-copy rejection mutations.
    require(len(mesh["facets"]) == 12, "invalid-case mutation leaked into the valid fixture")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "stage": 10,
        "scope": "tessellated_mesh_v1_convex_verification_only",
        "physical_validation_claim": False,
        "completion_status": "incomplete_pending_independent_nontrivial_shape_reference",
        "tested_examples": sorted(expected_examples),
        "provider_boundary": (
            "one connected closed consistently outward-oriented convex triangular shell per "
            "region; inline millimetre coordinates; no STL/CAD repair or concave meshes"
        ),
        "comparisons": comparisons,
    }
    evidence_path = EVIDENCE / "stage10_tessellated_mesh.json"
    evidence_path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    render_figure(comparisons)
    digest = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
    print(f"Stage-10 tessellated-mesh verification passed; evidence sha256={digest}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL Stage-10 tessellated-mesh verification: {error}", file=sys.stderr)
        raise SystemExit(1)

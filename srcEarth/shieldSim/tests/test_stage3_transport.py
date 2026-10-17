#!/usr/bin/env python3
"""Stage-3 Geant4 transport, isolation, response, and convergence acceptance.

This test deliberately writes every generated case, log, response cell, matrix,
fold, and evidence record below the caller-supplied build directory.  Nothing is
written to ``/tmp``.  The numerical scope is the Stage-3 raw transport engine;
the test does not enable or physically validate TID, DDD, n_eq, LET, H100/10,
instrument response, or a mission environment.

Five independent seeds are pooled for comparisons that make a statistical
claim.  The pooling uses event-level sums and sums of squares written by the
production scorer, so histories with zero score remain part of the estimator.
The preregistered acceptance limits are:

* relative Monte Carlo standard uncertainty <= 2% for each compared nonzero
  primitive response;
* the full two-sided 95% confidence interval for a difference lies within
  +/-5% of the folded/reference response;
* the final cut, step, transverse-size, and response-grid refinements change
  the designated raw transmission primitive by <= 2%; and
* quantities whose expected value is exactly zero are compared absolutely,
  never by an undefined relative error.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import pathlib
import shutil
import subprocess
import sys
from typing import Callable


if len(sys.argv) != 8:
    raise SystemExit(
        "usage: test_stage3_transport.py SHIELDSIM RESPONSE_TOOL WEIGHT_TOOL "
        "SOURCE_ROOT EVIDENCE_DIR RUN_DIR MODE"
    )

EXE = pathlib.Path(sys.argv[1]).resolve()
RESPONSE_TOOL = pathlib.Path(sys.argv[2]).resolve()
WEIGHT_TOOL = pathlib.Path(sys.argv[3]).resolve()
ROOT = pathlib.Path(sys.argv[4]).resolve()
EVIDENCE_DIR = pathlib.Path(sys.argv[5]).resolve()
RUN_DIR = pathlib.Path(sys.argv[6]).resolve()
MODE = sys.argv[7]
if MODE not in {"full", "smoke"}:
    raise SystemExit("MODE must be full or smoke")

BASE = json.loads((ROOT / "examples/stage3_response_cell_case.json").read_text(encoding="utf-8"))
BASE["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
SEEDS = [310001, 310019, 310043, 310081, 310111]
ABSORBER_SCORE = "absorber/silicon_proxy/total_deposited_energy_MeV"
EXIT_SCORE = "shield_exit/primary/proton/crossing_count"
EVENT_CONVERGENCE_SCORE = "shield_layer/outer_al/step_count"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def canonical_write(path: pathlib.Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def command(argv: list[str], cwd: pathlib.Path, log_name: str) -> None:
    """Run a production command and retain complete output beside its inputs."""
    completed = subprocess.run(
        argv,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (cwd / log_name).write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise AssertionError(
            f"command failed with rc={completed.returncode}: {' '.join(argv)}\n"
            f"complete output: {cwd / log_name}\n{completed.stdout[-6000:]}"
        )


def fixed_histories(case: dict, histories: int) -> None:
    # Equal min/max values disable adaptive early stopping while still driving
    # the normal production batch loop and accuracy reporting path.
    case["accuracy"]["target"] = 0.02
    case["accuracy"]["scope"] = "statistical"
    case["accuracy"]["sampling"] = {
        "min_histories": histories,
        "max_histories": histories,
        "batch_size": histories,
    }


def set_incident(case: dict, incident_id: str, energy_mev: float, particle: str = "proton",
                 z: int = 1, a: int = 1, per_nucleon: bool = False) -> None:
    spectrum = case["source"]["spectrum"]
    spectrum["format"] = "monoenergetic_response_cell"
    spectrum["incident_id"] = incident_id
    spectrum["particle"] = {"name": particle, "Z": z, "A": a}
    spectrum["energy"] = {
        "coordinate": "kinetic_energy_per_nucleon" if per_nucleon else "total_kinetic_energy",
        "value": energy_mev,
        "unit": "MeV/n" if per_nucleon else "MeV",
    }
    spectrum["direction"] = {"mu": 1.0, "phi": {"value": 0.0, "unit": "rad"}}


def result_prefix(directory: pathlib.Path, absorber_id: str = "silicon_proxy") -> pathlib.Path:
    return directory / f"run__{absorber_id}"


def run_case(case: dict, name: str, histories: int, *, validate_only: bool = False) -> tuple[dict, pathlib.Path]:
    """Write and run one isolated case, returning its result and output prefix."""
    directory = RUN_DIR / name
    if directory.exists():
        # This directory is wholly owned by this CTest case.  Recreating it
        # prevents stale response cells from being accepted after a failed run.
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    document = copy.deepcopy(case)
    document["case_id"] = name
    document["output"]["prefix"] = str(directory / "run")
    fixed_histories(document, histories)
    case_path = directory / "case.json"
    canonical_write(case_path, document)
    option = "--validate-case" if validate_only else "--case"
    command([str(EXE), f"{option}={case_path}"], directory, "shieldSim.log")
    absorber_id = document["geometry"]["absorbers"][0]["absorber_id"]
    prefix = result_prefix(directory, absorber_id)
    result_path = prefix.with_name(prefix.name + ("_resolved_case.json" if validate_only else "_result.json"))
    require(result_path.is_file(), f"missing production output {result_path}")
    return json.loads(result_path.read_text(encoding="utf-8")), prefix


def scores(document: dict, field: str = "raw_transport_scores") -> dict[str, dict]:
    return {entry["score_id"]: entry for entry in document[field]}


def combine_moments(documents: list[dict], score_id: str,
                    field: str = "raw_transport_scores") -> dict[str, float]:
    """Pool independent event samples, with absent dynamic scores counted as zero."""
    n = sum(float(document.get("events", document.get("histories", 0.0))) for document in documents)
    total = 0.0
    squares = 0.0
    unit = None
    for document in documents:
        match = scores(document, field).get(score_id)
        if match is not None:
            total += float(match["sum"])
            squares += float(match["sum_squares"])
            unit = match["unit"] if unit is None else unit
            require(unit == match["unit"], f"unit changed for {score_id}")
    require(n > 1.0, f"insufficient pooled histories for {score_id}")
    mean = total / n
    variance = max(0.0, (squares - total * total / n) / (n - 1.0))
    return {
        "histories": n,
        "sum": total,
        "sum_squares": squares,
        "mean": mean,
        "sample_variance": variance,
        "standard_error": math.sqrt(variance / n),
        "unit": unit or "structural_zero_per_primary",
    }


def combine_response_cells(paths: list[pathlib.Path], output: pathlib.Path) -> dict:
    """Pool seed replicas of one response node into one matrix-ready cell."""
    cells = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    require(cells, "response-cell replica list is empty")
    key = cells[0]["response_key"]
    key_document = cells[0]["response_key_document"]
    incident = cells[0]["incident"]
    require(all(cell["response_key"] == key for cell in cells), "replica response keys differ")
    require(all(cell["response_key_document"] == key_document for cell in cells),
            "replica response-key documents differ")
    require(all(cell["incident"] == incident for cell in cells), "replica incident nodes differ")
    score_ids = sorted({entry["score_id"] for cell in cells for entry in cell["scores"]})
    pooled_scores = []
    pooled_moments: dict[str, dict[str, float]] = {}
    for score_id in score_ids:
        moment = combine_moments(cells, score_id, "scores")
        pooled_moments[score_id] = moment
        pooled_scores.append({
            "histories": moment["histories"],
            "mean_per_primary": moment["mean"],
            "sample_variance": moment["sample_variance"],
            "score_id": score_id,
            "standard_error": moment["standard_error"],
            "sum": moment["sum"],
            "sum_squares": moment["sum_squares"],
            "unit": moment["unit"],
        })
    cell_covariance = [
        {(entry["score_id_i"], entry["score_id_j"]): entry
         for entry in cell["score_covariance"]}
        for cell in cells
    ]
    pooled_covariance = []
    total_histories = sum(float(cell["histories"]) for cell in cells)
    for first_index, first in enumerate(score_ids):
        for second in score_ids[first_index:]:
            sum_products = sum(
                float(covariance.get((first, second), {"sum_products": 0.0})["sum_products"])
                for covariance in cell_covariance
            )
            first_moment, second_moment = pooled_moments[first], pooled_moments[second]
            sample_covariance = (
                sum_products - first_moment["sum"] * second_moment["sum"] / total_histories
            ) / (total_histories - 1.0)
            if first == second:
                sample_covariance = max(0.0, sample_covariance)
            pooled_covariance.append({
                "covariance_of_mean": sample_covariance / total_histories,
                "histories": total_histories,
                "sample_covariance": sample_covariance,
                "score_id_i": first,
                "score_id_j": second,
                "sum_products": sum_products,
                "unit": first_moment["unit"] + " * " + second_moment["unit"],
            })
    pooled = {
        "absorber_id": cells[0]["absorber_id"],
        "case_hash": cells[0]["case_hash"],
        "histories": total_histories,
        "incident": incident,
        "normalization": cells[0]["normalization"],
        # Zero denotes an explicitly pooled realization. Individual seed IDs
        # remain archived in the input cells beside this derived cell.
        "random_seed": 0,
        "response_key": key,
        "response_key_document": key_document,
        "schema_version": "shieldSim.response_cell.v1",
        "score_covariance": pooled_covariance,
        "scores": pooled_scores,
        # Per-seed secondary ledgers are retained in the source cells.  Their
        # dynamic categories are already represented by the pooled score rows.
        "secondary_production_by_layer": [],
        "validation_status": "stage3_five_seed_pooled_transport_raw_pending_independent_review",
    }
    canonical_write(output, pooled)
    return pooled


def assert_precision(moment: dict[str, float], label: str) -> None:
    require(moment["mean"] != 0.0, f"{label}: nonzero comparison unexpectedly has zero mean")
    relative = moment["standard_error"] / abs(moment["mean"])
    require(relative <= 0.02, f"{label}: relative MC standard uncertainty {relative:.4%} > 2%")


def equivalence(reference: dict[str, float], candidate: dict[str, float], label: str,
                margin: float = 0.05) -> dict[str, float]:
    """Require an independent-sample 95% difference CI inside +/-margin."""
    assert_precision(reference, label + " reference")
    assert_precision(candidate, label + " candidate")
    delta = candidate["mean"] - reference["mean"]
    standard_error = math.hypot(reference["standard_error"], candidate["standard_error"])
    low, high = delta - 1.96 * standard_error, delta + 1.96 * standard_error
    bound = margin * abs(reference["mean"])
    require(low >= -bound and high <= bound,
            f"{label}: 95% difference CI [{low:.6g},{high:.6g}] is outside +/-{bound:.6g}")
    return {"delta": delta, "ci95_low": low, "ci95_high": high, "bound": bound}


def exact_zero(value: float, label: str) -> None:
    require(value == 0.0, f"{label}: expected exact zero, got {value!r}")


def run_ensemble(factory: Callable[[int], dict], stem: str, histories: int,
                 response_cells: bool = False) -> tuple[list[dict], list[pathlib.Path]]:
    documents: list[dict] = []
    cell_paths: list[pathlib.Path] = []
    for index, seed in enumerate(SEEDS):
        case = factory(seed)
        result, prefix = run_case(case, f"{stem}_seed{index + 1}", histories)
        documents.append(result)
        if response_cells:
            path = prefix.with_name(prefix.name + "_response_cell.json")
            require(path.is_file(), f"missing response cell {path}")
            cell_paths.append(path)
    return documents, cell_paths


def transform_point(q: tuple[float, float, float, float], p: list[float], t: list[float]) -> list[float]:
    """Independent quaternion transform used to construct a rotated fixture."""
    w, x, y, z = q
    px, py, pz = p
    return [
        (1 - 2 * (y * y + z * z)) * px + 2 * (x * y - w * z) * py
        + 2 * (x * z + w * y) * pz + t[0],
        2 * (x * y + w * z) * px + (1 - 2 * (x * x + z * z)) * py
        + 2 * (y * z - w * x) * pz + t[1],
        2 * (x * z - w * y) * px + 2 * (y * z + w * x) * py
        + (1 - 2 * (x * x + y * y)) * pz + t[2],
    ]


def set_transform(item: dict, position: list[float], q: tuple[float, float, float, float]) -> None:
    item["transform"] = {
        "translation": {
            "x": {"value": position[0], "unit": "mm"},
            "y": {"value": position[1], "unit": "mm"},
            "z": {"value": position[2], "unit": "mm"},
        },
        "rotation_quaternion_local_to_world": dict(zip(("w", "x", "y", "z"), q)),
    }


def main() -> int:
    if RUN_DIR.exists():
        # RUN_DIR is a dedicated CTest artifact directory supplied by CMake.
        shutil.rmtree(RUN_DIR)
    RUN_DIR.mkdir(parents=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    comparisons: list[dict] = []
    acceptance: list[dict] = []

    # --- Geometry isolation and local-frame boundary scoring -----------------
    angle = math.radians(20.0)
    q = (math.cos(angle / 2.0), 0.0, math.sin(angle / 2.0), 0.0)
    translation = [7.0, -5.0, 3.0]
    rotated = copy.deepcopy(BASE)
    for object_, position in (
        (rotated["geometry"]["source_plane"], [0.0, 0.0, -80.0]),
        (rotated["geometry"]["shield_stack"], [0.0, 0.0, -20.0]),
        (rotated["geometry"]["absorbers"][0]["geometry"], [0.0, 0.0, -16.0]),
    ):
        set_transform(object_, transform_point(q, position, translation), q)
    rotated_result, _ = run_case(rotated, "rotated_translated", 512)
    rotated_scores = scores(rotated_result)
    require(rotated_result["schema_version"] == "shieldSim.result.v2", "wrong Stage-3 result schema")
    require(rotated_scores[EXIT_SCORE]["mean_per_primary"] > 0.95,
            "rotated local downstream-face scorer lost normally incident primaries")
    for layer_id in ("outer_al", "inner_hdpe"):
        require(any(key.startswith(f"shield_layer/{layer_id}/") for key in rotated_scores),
                f"rotated case lacks layer-resolved scoring for {layer_id}")

    # Two absorbers in one input are alternative independent cases. Reversing
    # their selection order must preserve both geometry identities and common-
    # seed tallies.  Geant4 navigation can change accumulated floating path
    # lengths at roundoff level after a geometry rebuild, so continuous sums
    # and their squared moments use a 1e-4 numerical comparison while
    # integer/zero tallies remain exact.
    alternatives = copy.deepcopy(BASE)
    water = copy.deepcopy(alternatives["geometry"]["absorbers"][0])
    water["absorber_id"] = "water_proxy"
    water["material"] = "Water"
    water["geometry"]["thickness"] = {"value": 2.0, "unit": "mm"}
    alternatives["geometry"]["absorbers"].append(water)
    first_result, first_prefix = run_case(alternatives, "absorber_order_a", 256)
    reversed_case = copy.deepcopy(alternatives)
    reversed_case["geometry"]["absorbers"].reverse()
    second_result, second_prefix = run_case(reversed_case, "absorber_order_b", 256)
    # run_case returns the first absorber output; load every ID explicitly.
    order_results: list[dict[str, dict]] = []
    for directory, case in ((first_prefix.parent, alternatives), (second_prefix.parent, reversed_case)):
        by_id = {}
        for absorber in case["geometry"]["absorbers"]:
            prefix = result_prefix(directory, absorber["absorber_id"])
            by_id[absorber["absorber_id"]] = json.loads(
                prefix.with_name(prefix.name + "_result.json").read_text(encoding="utf-8")
            )
        order_results.append(by_id)
    for absorber_id in ("silicon_proxy", "water_proxy"):
        left, right = order_results[0][absorber_id], order_results[1][absorber_id]
        require(left["absorber_geometry_hash"] == right["absorber_geometry_hash"],
                f"{absorber_id}: geometry hash depends on absorber selection order")
        left_scores, right_scores = scores(left), scores(right)
        for score_id in sorted(set(left_scores) | set(right_scores)):
            l = left_scores.get(score_id, {"sum": 0.0, "sum_squares": 0.0})
            r = right_scores.get(score_id, {"sum": 0.0, "sum_squares": 0.0})
            for field in ("sum", "sum_squares"):
                scale = max(abs(float(l[field])), abs(float(r[field])))
                difference = abs(float(l[field]) - float(r[field]))
                if scale == 0.0:
                    require(difference == 0.0,
                            f"{absorber_id}: zero score {score_id}/{field} changed with order")
                else:
                    require(difference / scale <= 1.0e-4,
                            f"{absorber_id}: score {score_id}/{field} changed with order")

    # --- Five-seed same-material split equivalence ---------------------------
    # A thousand histories per replica gives the relatively noisy finite-slab
    # deposition score enough precision when five replicas are pooled.  Fast
    # smoke mode keeps the same preregistered statistical criterion; it reduces
    # the response grid, not the evidence quality.
    histories = 1000

    def homogeneous(seed: int, split: bool) -> dict:
        case = copy.deepcopy(BASE)
        case["transport"]["random_seed"] = seed
        case["geometry"]["shield_stack"]["layers"] = (
            [
                {"layer_id": "al_a", "material": "Al",
                 "thickness": {"type": "geometric", "value": 1.0, "unit": "mm"}},
                {"layer_id": "al_b", "material": "Al",
                 "thickness": {"type": "geometric", "value": 1.0, "unit": "mm"}},
            ] if split else [
                {"layer_id": "al_all", "material": "Al",
                 "thickness": {"type": "geometric", "value": 2.0, "unit": "mm"}}
            ]
        )
        case["geometry"]["absorbers"][0]["geometry"]["transform"]["translation"]["z"]["value"] = -17.0
        return case

    single_docs, _ = run_ensemble(lambda seed: homogeneous(seed, False), "single_al", histories)
    split_docs, _ = run_ensemble(lambda seed: homogeneous(seed, True), "split_al", histories)
    single = combine_moments(single_docs, ABSORBER_SCORE)
    split = combine_moments(split_docs, ABSORBER_SCORE)
    split_acceptance = equivalence(single, split, "single-vs-split Al absorber deposition")
    comparisons.append({"id": "single_vs_split_al_transport", "expected": single["mean"],
                        "actual": split["mean"], "unit": single["unit"]})
    acceptance.append({"id": "single_vs_split_al_transport", **split_acceptance})

    # --- Response cells, five-seed direct mixtures, and both material orders --
    pooled_by_order: dict[str, dict[float, tuple[dict, pathlib.Path]]] = {}
    energy_nodes = [80.0, 100.0, 110.0, 120.0, 150.0] if MODE == "full" else [80.0, 150.0]
    for order_name, materials in (("w_hdpe", ("W", "HDPE")), ("hdpe_w", ("HDPE", "W"))):
        pooled_by_order[order_name] = {}

        def response_case(seed: int, energy: float) -> dict:
            case = copy.deepcopy(BASE)
            case["transport"]["random_seed"] = seed
            case["geometry"]["shield_stack"]["layers"] = [
                {"layer_id": "layer_0", "material": materials[0],
                 "thickness": {"type": "geometric", "value": 0.5, "unit": "mm"}},
                {"layer_id": "layer_1", "material": materials[1],
                 "thickness": {"type": "geometric", "value": 2.0, "unit": "mm"}},
            ]
            case["geometry"]["absorbers"][0]["geometry"]["transform"]["translation"]["z"]["value"] = -16.0
            set_incident(case, f"proton_{energy:g}MeV_mu1_phi0", energy)
            return case

        for energy in energy_nodes:
            documents, paths = run_ensemble(
                lambda seed, e=energy: response_case(seed, e),
                f"{order_name}_{energy:g}MeV", histories, response_cells=True,
            )
            pooled_path = RUN_DIR / order_name / f"pooled_{energy:g}MeV.json"
            pooled_path.parent.mkdir(parents=True, exist_ok=True)
            pooled = combine_response_cells(paths, pooled_path)
            moment = combine_moments(documents, ABSORBER_SCORE)
            assert_precision(moment, f"{order_name} {energy:g} MeV response cell")
            pooled_by_order[order_name][energy] = (pooled, pooled_path)

        selected = [80.0, 150.0]
        matrix_path = RUN_DIR / order_name / "matrix.json"
        command([str(RESPONSE_TOOL), "assemble", str(matrix_path)] +
                [str(pooled_by_order[order_name][energy][1]) for energy in selected],
                matrix_path.parent, "assemble.log")
        matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
        weights_path = matrix_path.parent / "weights.json"
        canonical_write(weights_path, {
            "schema_version": "shieldSim.response_weights.v1",
            "matrix_hash": matrix["matrix_hash"],
            "normalization": "unit-probability validation mixture",
            "quadrature_policy": "exact two-node discrete distribution",
            "interpolation_policy": "exact_nodes",
            "weights": [
                {"incident_id": "proton_80MeV_mu1_phi0", "weight": 0.4},
                {"incident_id": "proton_150MeV_mu1_phi0", "weight": 0.6},
            ],
        })
        fold_path = matrix_path.parent / "fold.json"
        command([str(RESPONSE_TOOL), "fold", str(matrix_path), str(weights_path), str(fold_path)],
                matrix_path.parent, "fold.log")
        fold = json.loads(fold_path.read_text(encoding="utf-8"))
        folded = {item["score_id"]: item for item in fold["results"]}[ABSORBER_SCORE]

        def direct_case(seed: int) -> dict:
            case = response_case(seed, 80.0)
            case["source"]["spectrum"] = {
                "format": "discrete_response_nodes",
                "nodes": [
                    {"incident_id": "proton_80MeV_mu1_phi0",
                     "particle": {"name": "proton", "Z": 1, "A": 1},
                     "energy": {"coordinate": "total_kinetic_energy", "value": 80.0, "unit": "MeV"},
                     "direction": {"mu": 1.0, "phi": {"value": 0.0, "unit": "rad"}},
                     "weight": 0.4},
                    {"incident_id": "proton_150MeV_mu1_phi0",
                     "particle": {"name": "proton", "Z": 1, "A": 1},
                     "energy": {"coordinate": "total_kinetic_energy", "value": 150.0, "unit": "MeV"},
                     "direction": {"mu": 1.0, "phi": {"value": 0.0, "unit": "rad"}},
                     "weight": 0.6},
                ],
            }
            return case

        direct_docs, _ = run_ensemble(direct_case, f"{order_name}_direct", histories)
        direct = combine_moments(direct_docs, ABSORBER_SCORE)
        folded_moment = {
            "mean": float(folded["value"]),
            "standard_error": float(folded["matrix_statistical_standard_uncertainty"]),
            "unit": folded["unit"],
        }
        direct_acceptance = equivalence(folded_moment, direct,
                                        f"{order_name} direct-vs-response fold")
        comparisons.append({"id": f"{order_name}_direct_vs_fold", "expected": folded_moment["mean"],
                            "actual": direct["mean"], "unit": direct["unit"]})
        acceptance.append({"id": f"{order_name}_direct_vs_fold", **direct_acceptance})

        # Every physical layer must retain its own production diagnostics.  We
        # assert ledger structure, not an unreferenced expected ordering effect.
        source_cell = pooled_by_order[order_name][80.0][0]
        score_ids = {item["score_id"] for item in source_cell["scores"]}
        for layer_id in ("layer_0", "layer_1"):
            require(f"shield_layer/{layer_id}/step_count" in score_ids,
                    f"{order_name}: missing step ledger for {layer_id}")
            require(any(item.startswith(f"shield_layer/{layer_id}/secondary_birth/")
                        for item in score_ids),
                    f"{order_name}: no layer-resolved secondary-production score for {layer_id}")

    # --- Actual response-grid refinement and weight-builder contract ---------
    if MODE == "full":
        order_dir = RUN_DIR / "w_hdpe"

        def fold_grid(label: str, energies: list[float], query: float) -> float:
            matrix_path = order_dir / f"grid_{label}.json"
            command([str(RESPONSE_TOOL), "assemble", str(matrix_path)] +
                    [str(pooled_by_order["w_hdpe"][energy][1]) for energy in energies],
                    order_dir, f"grid_{label}_assemble.log")
            source_path = order_dir / f"grid_{label}_source.json"
            canonical_write(source_path, {
                "schema_version": "shieldSim.response_source_nodes.v1",
                "normalization": "one validation primary",
                "quadrature_policy": "single exact query node",
                "nodes": [{"particle": "proton", "Z": 1, "A": 1,
                           "total_kinetic_energy_MeV": query, "mu": 1.0,
                           "phi_rad": 0.0, "probability": 1.0}],
            })
            weights = order_dir / f"grid_{label}_weights.json"
            command([sys.executable, str(WEIGHT_TOOL), str(matrix_path), str(source_path), str(weights)],
                    order_dir, f"grid_{label}_weights.log")
            fold_path = order_dir / f"grid_{label}_fold.json"
            command([str(RESPONSE_TOOL), "fold", str(matrix_path), str(weights), str(fold_path)],
                    order_dir, f"grid_{label}_fold.log")
            fold = json.loads(fold_path.read_text(encoding="utf-8"))
            return float({item["score_id"]: item for item in fold["results"]}[ABSORBER_SCORE]["value"])

        # The fine grid contains the query exactly.  The last refinement is the
        # medium [100,120] interpolation to that exact 110-MeV response.
        coarse_value = fold_grid("coarse", [80.0, 150.0], 110.0)
        medium_value = fold_grid("medium", [100.0, 120.0], 110.0)
        fine_value = fold_grid("fine", [100.0, 110.0, 120.0], 110.0)
        grid_change = abs(medium_value - fine_value) / abs(fine_value)
        require(grid_change <= 0.02, f"last response-grid refinement changed score by {grid_change:.4%}")
        comparisons.append({"id": "response_grid_medium_vs_fine", "expected": fine_value,
                            "actual": medium_value, "unit": "MeV/primary"})
        acceptance.append({"id": "response_grid_medium_vs_fine",
                           "coarse_value": coarse_value, "relative_change": grid_change,
                           "limit": 0.02})

    # --- Cut, step, and finite-transverse-size plateau drivers ---------------
    def convergence_case(name: str, value: float) -> dict:
        case = copy.deepcopy(BASE)
        case["transport"]["random_seed"] = 911337
        if name == "cut":
            case["transport"]["production_cut"] = {"value": value, "unit": "mm"}
        elif name == "step":
            case["transport"]["max_step"] = {"value": value, "unit": "mm"}
        elif name == "size":
            for item in (case["geometry"]["source_plane"], case["geometry"]["shield_stack"]):
                item["size_x"] = {"value": value, "unit": "mm"}
                item["size_y"] = {"value": value, "unit": "mm"}
            absorber = case["geometry"]["absorbers"][0]["geometry"]
            absorber["size_x"] = {"value": value, "unit": "mm"}
            absorber["size_y"] = {"value": value, "unit": "mm"}
        elif name == "source_plane":
            # Hold the shield and absorber at 100 mm while expanding only the
            # sampled plane.  This separates source-plane convergence from the
            # common transverse-size expansion above.
            case["geometry"]["source_plane"]["size_x"] = {"value": value, "unit": "mm"}
            case["geometry"]["source_plane"]["size_y"] = {"value": value, "unit": "mm"}
            for item in (case["geometry"]["shield_stack"],
                         case["geometry"]["absorbers"][0]["geometry"]):
                item["size_x"] = {"value": 100.0, "unit": "mm"}
                item["size_y"] = {"value": 100.0, "unit": "mm"}
        else:
            raise AssertionError(name)
        return case

    refinement_values = {
        "cut": [1.0, 0.1, 0.03],
        "step": [1.0, 0.2, 0.05],
        "size": [30.0, 50.0, 80.0],
        "source_plane": [30.0, 50.0, 80.0],
    }
    if MODE == "smoke":
        refinement_values = {key: values[-2:] for key, values in refinement_values.items()}
    for parameter, values in refinement_values.items():
        moments = []
        for value in values:
            result, _ = run_case(convergence_case(parameter, value),
                                 f"convergence_{parameter}_{value:g}", 1000)
            moments.append(combine_moments([result], EXIT_SCORE))
        final_change = abs(moments[-2]["mean"] - moments[-1]["mean"]) / abs(moments[-1]["mean"])
        require(final_change <= 0.02,
                f"last {parameter} refinement changed raw primary transmission by {final_change:.4%}")
        comparisons.append({"id": f"{parameter}_last_refinement", "expected": moments[-1]["mean"],
                            "actual": moments[-2]["mean"], "unit": moments[-1]["unit"]})
        acceptance.append({"id": f"{parameter}_last_refinement",
                           "tested_values": values, "relative_change": final_change, "limit": 0.02})

    # Event-count convergence is evaluated with the well-populated outer-layer
    # step-count primitive and independent five-seed ensembles. Quadrupling N
    # should reduce its standard error by approximately two. Deposited energy
    # contains rare nuclear tails and has its own 5,000-history precision gate
    # above; convergence criteria are intentionally quantity-specific.
    def event_case(seed: int) -> dict:
        case = copy.deepcopy(BASE)
        case["transport"]["random_seed"] = seed
        return case

    event_small_docs, _ = run_ensemble(event_case, "event_convergence_250", 250)
    event_large_docs, _ = run_ensemble(event_case, "event_convergence_1000", 1000)
    event_small = combine_moments(event_small_docs, EVENT_CONVERGENCE_SCORE)
    event_large = combine_moments(event_large_docs, EVENT_CONVERGENCE_SCORE)
    assert_precision(event_large, "event-convergence accepted ensemble")
    se_ratio = event_small["standard_error"] / event_large["standard_error"]
    require(1.5 <= se_ratio <= 2.5,
            f"quadrupling event count gave SE ratio {se_ratio:.3f}, outside [1.5,2.5]")
    comparisons.append({"id": "event_standard_error_sqrtN", "expected": 2.0,
                        "actual": se_ratio, "unit": "1"})
    acceptance.append({"id": "event_standard_error_sqrtN", "expected_ratio": 2.0,
                       "observed_ratio": se_ratio, "acceptance_interval": [1.5, 2.5]})

    # --- Species/material/physics configuration execution -------------------
    smoke_cases: list[tuple[str, dict]] = []
    alpha = copy.deepcopy(BASE)
    set_incident(alpha, "alpha_100MeVn", 100.0, "alpha", 2, 4, per_nucleon=True)
    smoke_cases.append(("species_alpha4", alpha))
    iron = copy.deepcopy(BASE)
    set_incident(iron, "fe56_100MeVn", 100.0, "ion", 26, 56, per_nucleon=True)
    smoke_cases.append(("species_fe56_default", iron))
    iron_qmd = copy.deepcopy(iron)
    iron_qmd["transport"]["ion_physics"] = "G4IonQMDPhysics"
    smoke_cases.append(("species_fe56_qmd", iron_qmd))
    hp = copy.deepcopy(BASE)
    hp["transport"]["physics_list"] = "FTFP_BERT_HP"
    smoke_cases.append(("physics_ftfp_bert_hp", hp))
    tungsten_water = copy.deepcopy(BASE)
    tungsten_water["geometry"]["shield_stack"]["layers"] = [{
        "layer_id": "tungsten", "material": "W",
        "thickness": {"type": "geometric", "value": 0.5, "unit": "mm"},
    }]
    tungsten_water["geometry"]["absorbers"][0]["material"] = "Water"
    tungsten_water["geometry"]["absorbers"][0]["absorber_id"] = "water_proxy"
    tungsten_water["geometry"]["absorbers"][0]["geometry"]["transform"]["translation"]["z"]["value"] = -18.0
    smoke_cases.append(("materials_w_water", tungsten_water))
    for index, (name, case) in enumerate(smoke_cases):
        case["transport"]["random_seed"] = 700001 + index * 20
        result, _ = run_case(case, name, 32)
        require(result["raw_transport_scores"], f"{name}: transport produced no raw scores")
        require(result["physics_list"] == case["transport"]["physics_list"],
                f"{name}: physics-list provenance mismatch")
        require(result["ion_physics"] == case["transport"]["ion_physics"],
                f"{name}: ion-physics provenance mismatch")
        require(result["energy_ledger_scope"]["nuclear_rest_mass_and_q_values"].startswith("not included"),
                f"{name}: raw kinetic/deposition ledger scope is not explicit")

    # A stopped-primary extreme supplies an exact zero-transmission case.  The
    # absence of a dynamic score row is defined as zero over all histories by
    # the event-level accumulator; it is not assigned a relative uncertainty.
    stopped = copy.deepcopy(BASE)
    stopped["geometry"]["shield_stack"]["layers"] = [{
        "layer_id": "thick_w", "material": "W",
        "thickness": {"type": "geometric", "value": 10.0, "unit": "mm"},
    }]
    stopped["geometry"]["absorbers"][0]["geometry"]["transform"]["translation"]["z"]["value"] = -9.0
    set_incident(stopped, "stopped_proton_2MeV", 2.0)
    stopped_result, _ = run_case(stopped, "extreme_stopped_primary", 256)
    stopped_exit = scores(stopped_result).get(EXIT_SCORE, {"sum": 0.0})["sum"]
    exact_zero(float(stopped_exit), "thick-W primary transmission")
    comparisons.append({"id": "stopped_primary_zero_transmission", "expected": 0.0,
                        "actual": float(stopped_exit), "unit": "crossings"})

    # A mandatory score with no nonionizing deposit in the controlled proton
    # fixture exercises the exact-zero path; the response fold must preserve it
    # exactly rather than applying a relative tolerance to 0/0.
    zero_score = "absorber/silicon_proxy/nonionizing_deposited_energy_MeV"
    for order_name in pooled_by_order:
        for energy, (cell, _) in pooled_by_order[order_name].items():
            match = scores(cell, "scores").get(zero_score)
            if match is not None and match["sum"] == 0.0:
                exact_zero(float(match["standard_error"]), f"{order_name} {energy:g} zero-score uncertainty")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_case": "stage3_transport",
        "scope": "raw Stage-3 transport/geometry/response verification only",
        "stage_status": "implemented_and_acceptance_tested_pending_independent_review",
        "independent_seed_count": len(SEEDS),
        "acceptance_policy": {
            "mc_relative_standard_uncertainty_max": 0.02,
            "equivalence_95_percent_relative_margin": 0.05,
            "last_refinement_relative_change_max": 0.02,
            "zero_expected_policy": "exact_absolute_comparison",
        },
        "comparisons": comparisons,
        "acceptance_details": acceptance,
        "run_directory": str(RUN_DIR),
    }
    canonical_write(EVIDENCE_DIR / "stage3_transport.json", evidence)
    canonical_write(RUN_DIR / "acceptance_summary.json", evidence)
    print(f"Stage-3 transport acceptance passed; artifacts: {RUN_DIR}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, KeyError, ValueError, OSError, json.JSONDecodeError) as error:
        print(f"FAIL Stage-3 transport acceptance: {error}", file=sys.stderr)
        raise SystemExit(1)

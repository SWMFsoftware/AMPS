#!/usr/bin/env python3
"""Exact tests for Stage-3 energy/mu/periodic-phi response weights.

Synthetic response values are affine in the interpolation coordinates, so the
multilinear answer is analytic.  The production matrix assembler, Python
weight builder, and production folder are all exercised; no transport result
is used as its own oracle.
"""

from __future__ import annotations

import hashlib
import json
import math
import pathlib
import shutil
import subprocess
import sys


if len(sys.argv) != 6:
    raise SystemExit(
        "usage: test_stage3_response_weights.py RESPONSE_TOOL WEIGHT_TOOL WORK_DIR EVIDENCE_DIR PYTHON"
    )
RESPONSE_TOOL = pathlib.Path(sys.argv[1]).resolve()
WEIGHT_TOOL = pathlib.Path(sys.argv[2]).resolve()
WORK = pathlib.Path(sys.argv[3]).resolve()
EVIDENCE = pathlib.Path(sys.argv[4]).resolve()
PYTHON = sys.argv[5]


def canonical(document: dict) -> str:
    return json.dumps(document, sort_keys=True, separators=(",", ":"))


def write(path: pathlib.Path, document: dict) -> None:
    path.write_text(canonical(document) + "\n", encoding="utf-8")


def run(argv: list[str], *, succeeds: bool, log: str) -> None:
    completed = subprocess.run(
        argv, cwd=WORK, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False
    )
    (WORK / log).write_text(completed.stdout, encoding="utf-8")
    if (completed.returncode == 0) != succeeds:
        raise AssertionError(
            f"unexpected rc={completed.returncode}: {' '.join(argv)}\n{completed.stdout}"
        )


def main() -> int:
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    key_document = {"synthetic_contract": "stage3_energy_mu_phi_tensor_v1"}
    response_key = "sha256:" + hashlib.sha256(canonical(key_document).encode()).hexdigest()
    cells: list[pathlib.Path] = []
    expected_corner_values: dict[str, float] = {}
    for energy in (100.0, 200.0):
        for mu in (0.5, 1.0):
            for phi in (0.0, math.pi):
                incident_id = f"e{energy:g}_m{mu:g}_p{phi:.12g}"
                # The declared policy interpolates response linearly in log(E),
                # mu, and phi between the selected periodic nodes.
                value = 3.0 * math.log(energy) + 2.0 * mu + 0.5 * phi
                expected_corner_values[incident_id] = value
                cell = {
                    "absorber_id": "synthetic",
                    "case_hash": "synthetic",
                    "histories": 1,
                    "incident": {
                        "A": 1, "Z": 1,
                        "energy_coordinate": "total_kinetic_energy",
                        "energy_coordinate_unit": "MeV",
                        "energy_coordinate_value": energy,
                        "incident_id": incident_id,
                        "mu": mu, "particle": "proton", "phi_rad": phi,
                        "total_kinetic_energy_MeV": energy,
                    },
                    "normalization": "synthetic per primary",
                    "random_seed": 1,
                    "response_key": response_key,
                    "response_key_document": key_document,
                    "schema_version": "shieldSim.response_cell.v1",
                    "score_covariance": [{
                        "covariance_of_mean": 0.0,
                        "histories": 1,
                        "sample_covariance": 0.0,
                        "score_id_i": "analytic",
                        "score_id_j": "analytic",
                        "sum_products": value * value,
                        "unit": "1 * 1",
                    }],
                    "scores": [{
                        "histories": 1,
                        "mean_per_primary": value,
                        "sample_variance": 0.0,
                        "score_id": "analytic",
                        "standard_error": 0.0,
                        "sum": value,
                        "sum_squares": value * value,
                        "unit": "1",
                    }],
                    "secondary_production_by_layer": [],
                    "validation_status": "synthetic_analytic_test",
                }
                path = WORK / f"{incident_id}.json"
                write(path, cell)
                cells.append(path)

    matrix = WORK / "matrix.json"
    run([str(RESPONSE_TOOL), "assemble", str(matrix)] + [str(path) for path in cells],
        succeeds=True, log="assemble.log")
    query_energy = math.sqrt(100.0 * 200.0)
    source = {
        "schema_version": "shieldSim.response_source_nodes.v1",
        "normalization": "one analytic query",
        "quadrature_policy": "single exact query",
        "nodes": [{"particle": "proton", "Z": 1, "A": 1,
                   "total_kinetic_energy_MeV": query_energy,
                   "mu": 0.75, "phi_rad": math.pi / 2.0, "probability": 1.0}],
    }
    source_path, weights_path = WORK / "source.json", WORK / "weights.json"
    write(source_path, source)
    run([PYTHON, str(WEIGHT_TOOL), str(matrix), str(source_path), str(weights_path)],
        succeeds=True, log="weights.log")
    weights = json.loads(weights_path.read_text(encoding="utf-8"))
    assert len(weights["weights"]) == 8
    assert math.isclose(sum(item["weight"] for item in weights["weights"]), 1.0,
                        rel_tol=0.0, abs_tol=1.0e-12)
    assert all(math.isclose(item["weight"], 0.125, rel_tol=0.0, abs_tol=1.0e-12)
               for item in weights["weights"])
    fold_path = WORK / "fold.json"
    run([str(RESPONSE_TOOL), "fold", str(matrix), str(weights_path), str(fold_path)],
        succeeds=True, log="fold.log")
    folded = json.loads(fold_path.read_text(encoding="utf-8"))["results"][0]["value"]
    expected = 3.0 * math.log(query_energy) + 2.0 * 0.75 + 0.5 * math.pi / 2.0
    assert math.isclose(folded, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)

    # Exact-node mode must reject this midpoint, and interpolation mode must
    # reject energy extrapolation and an incomplete surrounding tensor grid.
    exact_path = WORK / "exact.json"
    run([PYTHON, str(WEIGHT_TOOL), str(matrix), str(source_path), str(exact_path),
         "--policy", "exact_nodes"], succeeds=False, log="reject_exact.log")
    outside = json.loads(canonical(source))
    outside["nodes"][0]["total_kinetic_energy_MeV"] = 400.0
    outside_path = WORK / "outside.json"
    write(outside_path, outside)
    run([PYTHON, str(WEIGHT_TOOL), str(matrix), str(outside_path), str(WORK / "outside_weights.json")],
        succeeds=False, log="reject_extrapolation.log")
    incomplete = WORK / "incomplete_matrix.json"
    run([str(RESPONSE_TOOL), "assemble", str(incomplete)] + [str(path) for path in cells[:-1]],
        succeeds=True, log="assemble_incomplete.log")
    run([PYTHON, str(WEIGHT_TOOL), str(incomplete), str(source_path),
         str(WORK / "incomplete_weights.json")], succeeds=False, log="reject_incomplete.log")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_case": "stage3_response_weights",
        "comparisons": [{"id": "logE_mu_periodic_phi_multilinear_fold",
                         "expected": expected, "actual": folded, "unit": "1"}],
    }
    write(EVIDENCE / "stage3_response_weights.json", evidence)
    print("Stage-3 response-weight tests passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL Stage-3 response-weight test: {error}", file=sys.stderr)
        raise SystemExit(1)

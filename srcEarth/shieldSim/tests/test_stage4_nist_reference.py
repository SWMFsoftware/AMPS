#!/usr/bin/env python3
"""Stage-4 deterministic Geant4-to-NIST reference acceptance test.

The acceptance plan is committed independently of executable output.  This
driver first authenticates every pinned fixture, then runs the reference-table
probe with the production EM configuration and the option-4 sensitivity.  All
comparisons occur at exact NIST table nodes, so interpolation error is exactly
zero for this gate.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import pathlib
import subprocess
import sys


if len(sys.argv) != 5:
    raise SystemExit(
        "usage: test_stage4_nist_reference.py REFERENCE_TOOL SOURCE_ROOT EVIDENCE_DIR RUN_DIR"
    )

TOOL = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
FIXTURE = ROOT / "tests" / "reference_tid" / "nist"
PLAN_PATH = ROOT / "tests" / "reference_tid" / "stage4_reference_plan.v1.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_row(path: pathlib.Path, energy: float) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        matches = [row for row in csv.DictReader(stream)
                   if math.isclose(float(row["energy_MeV"]), energy,
                                   rel_tol=0.0, abs_tol=1.0e-12)]
    require(len(matches) == 1, f"{path.name} does not have one unambiguous row at {energy} MeV")
    return matches[0]


def run_probe(request: dict, stem: str) -> dict:
    request_path = RUN_DIR / f"{stem}_request.json"
    result_path = RUN_DIR / f"{stem}_result.json"
    request_path.write_text(json.dumps(request, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(TOOL), str(request_path), str(result_path)],
        cwd=RUN_DIR,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (RUN_DIR / f"{stem}.log").write_text(completed.stdout, encoding="utf-8")
    require(completed.returncode == 0, f"{stem} reference probe failed; see {stem}.log")
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    require(payload["schema_version"] == "shieldSim.stage4_reference_result.v1",
            "reference tool returned an unexpected schema")
    require(payload["csda_range_method"].startswith("composite_Simpson_quadrature"),
            "reference tool did not archive its CSDA construction method")
    return payload


def main() -> int:
    require(TOOL.is_file(), f"missing reference tool: {TOOL}")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    require(plan["registration_status"] == "registered_before_first_model_reference_comparison",
            "acceptance plan is not preregistered")

    manifest = json.loads((FIXTURE / "manifest.sha256.json").read_text(encoding="utf-8"))
    for relative, record in manifest["files"].items():
        path = FIXTURE / relative
        require(path.is_file(), f"pinned NIST file is missing: {relative}")
        require(path.stat().st_size == record["size_bytes"], f"size mismatch: {relative}")
        require(sha256(path) == record["sha256"], f"SHA-256 mismatch: {relative}")

    material_names = {"aluminum": "G4_Al", "silicon": "G4_Si", "water": "G4_WATER"}
    queries: list[dict] = []
    metadata: dict[str, dict] = {}
    for case in plan["charged_particle_cases"]:
        for material, geant_name in material_names.items():
            identifier = f"{case['particle']}_{material}_{case['energy_MeV']:g}MeV"
            queries.append({
                "id": identifier,
                "particle": case["particle"],
                "material": geant_name,
                "energy_MeV": case["energy_MeV"],
                "quantities": ["electronic_stopping_power", "total_stopping_power", "csda_range"],
            })
            metadata[identifier] = {"kind": "charged", "material": material, "case": case}
    for case in plan["photon_cases"]:
        for material, geant_name in material_names.items():
            identifier = f"gamma_{material}_{case['energy_MeV']:g}MeV"
            queries.append({
                "id": identifier,
                "particle": "gamma",
                "material": geant_name,
                "energy_MeV": case["energy_MeV"],
                "quantities": ["mass_attenuation_coefficient"],
            })
            metadata[identifier] = {"kind": "photon", "material": material, "case": case}

    common = {
        "schema_version": "shieldSim.stage4_reference_request.v1",
        "physics_list": "FTFP_BERT",
        "ion_physics": "reference_list_default",
        "queries": queries,
    }
    baseline = run_probe({**common, "em_physics": "reference_list_default"}, "baseline")
    option4 = run_probe({**common, "em_physics": "G4EmStandardPhysics_option4"}, "option4")
    baseline_by_id = {row["id"]: row for row in baseline["results"]}
    option4_by_id = {row["id"]: row for row in option4["results"]}
    require(set(baseline_by_id) == set(metadata) == set(option4_by_id),
            "reference probe omitted or added queries")

    comparisons: list[dict] = []
    sensitivity: list[dict] = []
    for identifier, info in metadata.items():
        base = baseline_by_id[identifier]
        alt = option4_by_id[identifier]
        material = info["material"]
        case = info["case"]
        if info["kind"] == "charged":
            row = exact_row(FIXTURE / f"{case['table']}_{material}.csv", case["energy_MeV"])
            fields = [
                ("electronic_stopping_power", "electronic_stopping_MeV_cm2_g",
                 case["electronic_acceptance_relative"][material]),
                # Total stopping is governed by the same predeclared envelope;
                # the nuclear term is retained in the reference CSV and is
                # small at these Stage-4 TID validation energies.
                ("total_stopping_power", "total_stopping_MeV_cm2_g",
                 case["electronic_acceptance_relative"][material]),
                ("csda_range", "csda_range_g_cm2",
                 case["csda_acceptance_relative"][material]),
            ]
        else:
            row = exact_row(FIXTURE / f"xcom_{material}.csv", case["energy_MeV"])
            fields = [("mass_attenuation_coefficient", "mass_attenuation_with_coherent_cm2_g",
                       case["mass_attenuation_acceptance_relative"][material])]
        for quantity, reference_field, tolerance in fields:
            expected = float(row[reference_field])
            actual = float(base["values"][quantity]["value"])
            relative = abs(actual - expected) / expected
            comparisons.append({
                "id": f"{identifier}_{quantity}",
                "unit": base["values"][quantity]["unit"],
                "expected": expected,
                "actual": actual,
                "relative_difference": relative,
                "acceptance_relative": tolerance,
                "interpolation_relative_uncertainty": 0.0,
                "passed": relative <= tolerance,
                "reference": f"NIST {case.get('table', 'XCOM')} exact table node",
            })
            alt_value = float(alt["values"][quantity]["value"])
            sensitivity.append({
                "id": f"{identifier}_{quantity}_option4_sensitivity",
                "unit": base["values"][quantity]["unit"],
                "baseline": actual,
                "option4": alt_value,
                "relative_difference": abs(alt_value - actual) / abs(actual) if actual else 0.0,
            })

    failed = [item for item in comparisons if not item["passed"]]
    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_id": "stage4_nist_reference",
        "status": "passed" if not failed else "failed",
        "scope": "deterministic Geant4 EM table comparison; not a finite-slab dose comparison",
        "validated_domain": plan["validated_domain"],
        "acceptance_basis": plan["acceptance_basis"],
        "fixture_manifest_sha256": sha256(FIXTURE / "manifest.sha256.json"),
        "comparisons": comparisons,
        "em_option4_sensitivity": sensitivity,
    }
    (EVIDENCE / "stage4_nist_reference.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if failed:
        details = ", ".join(
            f"{item['id']}={item['relative_difference']:.3g}>{item['acceptance_relative']:.3g}"
            for item in failed
        )
        raise AssertionError(f"NIST reference acceptance failed: {details}")
    print(f"Stage-4 NIST reference checks passed: {len(comparisons)} comparisons")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

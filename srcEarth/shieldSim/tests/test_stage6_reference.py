#!/usr/bin/env python3
"""Stage-6 electronic-LET reference comparisons.

Proton/He-4 comparisons reuse exact nodes from the authenticated NIST
PSTAR/ASTAR package. The selected heavy-ion implementation screen compares
fully stripped C-12 in natural silicon with checksum-pinned GSI ATIMA v1.41
calculator responses. ATIMA's 15% screen is not interpreted as an experimental
uncertainty or as validation outside the three registered energies.
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
        "usage: test_stage6_reference.py REFERENCE_TOOL SOURCE_ROOT EVIDENCE_DIR RUN_DIR"
    )

TOOL = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
NIST = ROOT / "tests/reference_tid/nist"
ATIMA = ROOT / "tests/reference_let/atima"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def rows(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def exact_row(path: pathlib.Path, energy: float) -> dict[str, str]:
    selected = [row for row in rows(path)
                if math.isclose(float(row["energy_MeV"]), energy,
                                rel_tol=0.0, abs_tol=1.0e-12)]
    require(len(selected) == 1, f"{path.name}: missing exact {energy:g} MeV node")
    return selected[0]


def main() -> int:
    require(TOOL.is_file(), f"missing reference tool: {TOOL}")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((NIST / "manifest.sha256.json").read_text(encoding="utf-8"))
    for relative, record in manifest["files"].items():
        path = NIST / relative
        require(path.stat().st_size == record["size_bytes"] and
                sha256(path) == record["sha256"], f"NIST fixture changed: {relative}")
    for line in (ATIMA / "manifest.sha256").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split(maxsplit=1)
        require(sha256(ATIMA / relative) == expected,
                f"ATIMA fixture changed: {relative}")

    plan = json.loads(
        (ROOT / "tests/reference_tid/stage4_reference_plan.v1.json").read_text(encoding="utf-8")
    )
    material_names = {"aluminum": "G4_Al", "silicon": "G4_Si", "water": "G4_WATER"}
    queries: list[dict] = []
    expected: dict[str, tuple[float, float, str]] = {}
    for case in plan["charged_particle_cases"]:
        for material, geant_material in material_names.items():
            identifier = f"{case['particle']}_{material}_{case['energy_MeV']:g}MeV"
            queries.append({"id": identifier, "particle": case["particle"],
                            "material": geant_material, "energy_MeV": case["energy_MeV"],
                            "quantities": ["electronic_stopping_power"]})
            reference = exact_row(NIST / f"{case['table']}_{material}.csv",
                                  float(case["energy_MeV"]))
            expected[identifier] = (
                float(reference["electronic_stopping_MeV_cm2_g"]),
                float(case["electronic_acceptance_relative"][material]),
                f"NIST {case['table'].upper()} exact node",
            )
    for row in rows(ATIMA / "c12_silicon_electronic_stopping.csv"):
        identifier = f"c12_silicon_{float(row['energy_MeV_per_u']):g}MeV_per_u"
        queries.append({"id": identifier, "particle": "ion", "Z": 6, "A": 12,
                        "material": "G4_Si", "energy_MeV": float(row["total_energy_MeV"]),
                        "quantities": ["electronic_stopping_power"]})
        expected[identifier] = (
            float(row["electronic_stopping_MeV_cm2_per_g"]), 0.15,
            "GSI ATIMA v1.41 web response; implementation/reference screen",
        )

    def probe(stem: str, selected: list[dict], em_physics: str) -> dict:
        request = {
            "schema_version": "shieldSim.stage6_reference_request.v1",
            "physics_list": "FTFP_BERT", "em_physics": em_physics,
            "ion_physics": "reference_list_default", "queries": selected,
        }
        request_path = RUN_DIR / f"{stem}_request.json"
        result_path = RUN_DIR / f"{stem}_result.json"
        request_path.write_text(json.dumps(request, indent=2, sort_keys=True) + "\n",
                                encoding="utf-8")
        completed = subprocess.run([str(TOOL), str(request_path), str(result_path)], cwd=RUN_DIR,
                                   text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   check=False)
        (RUN_DIR / f"{stem}.log").write_text(completed.stdout, encoding="utf-8")
        require(completed.returncode == 0,
                f"Stage-6 reference tool failed; see {RUN_DIR/f'{stem}.log'}")
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        require(payload["schema_version"] == "shieldSim.stage6_reference_result.v1",
                "reference result has the wrong schema")
        return payload

    # The fixed screen selects the reference-list default for p/He-4 and
    # option 4 for C-12. The latter is an explicit domain boundary: the default
    # model misses the 10 MeV/u ATIMA row by more than the predeclared 15%.
    baseline = probe("stage6_p_he_baseline",
                     [query for query in queries if query["particle"] != "ion"],
                     "reference_list_default")
    carbon = probe("stage6_c12_option4",
                   [query for query in queries if query["particle"] == "ion"],
                   "G4EmStandardPhysics_option4")
    result_rows = baseline["results"] + carbon["results"]

    comparisons: list[dict] = []
    for result in result_rows:
        reference, tolerance, provenance = expected.pop(result["id"])
        actual = float(result["values"]["electronic_stopping_power"]["value"])
        relative = abs(actual - reference) / reference
        comparisons.append({
            "id": result["id"], "unit": "MeV cm2/g", "expected": reference,
            "actual": actual, "relative_difference": relative,
            "acceptance_relative": tolerance, "passed": relative <= tolerance,
            "reference": provenance,
        })
    require(not expected, "reference tool omitted registered queries")
    failed = [item for item in comparisons if not item["passed"]]
    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_case": "stage6_electronic_let_reference",
        "scope": "deterministic unrestricted electronic stopping comparison",
        "heavy_ion_interpretation": "15% implementation/reference screen; not uncertainty or experimental validation",
        "comparisons": comparisons,
    }
    (EVIDENCE / "stage6_reference.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if failed:
        raise AssertionError("Stage-6 reference failures: " + ", ".join(
            f"{item['id']}={item['relative_difference']:.4g}>{item['acceptance_relative']:.4g}"
            for item in failed
        ))
    print(f"Stage-6 electronic stopping checks passed: {len(comparisons)} comparisons")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

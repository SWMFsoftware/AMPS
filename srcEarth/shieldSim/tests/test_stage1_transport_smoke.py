#!/usr/bin/env python3
"""Short Geant4 integration test for Stage-1 geometry and score ledgers."""

from __future__ import annotations

import copy
import json
import math
import pathlib
import subprocess
import sys
import tempfile


if len(sys.argv) not in (3, 4):
    raise SystemExit("usage: test_stage1_transport_smoke.py SHIELDSIM SOURCE_ROOT [EVIDENCE_DIR]")
EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE_DIR = pathlib.Path(sys.argv[3]).resolve() if len(sys.argv) == 4 else None


def parse_summary(path: pathlib.Path) -> dict:
    result: dict[str, object] = {"layers": [], "targets": [], "ledgers": []}
    for raw in path.read_text().splitlines():
        fields = raw.split()
        if not fields or fields[0] == "#":
            continue
        if fields[0] == "scalar":
            result[fields[1]] = float(fields[2])
        elif fields[0] == "count":
            result[fields[1]] = float(fields[2])
        elif fields[0] == "shield_layer_material":
            result["layers"].append((int(fields[1]), fields[2], float(fields[3])))
        elif fields[0] == "target_material":
            result["targets"].append(
                (int(fields[1]), fields[2], float(fields[3]), float(fields[4]))
            )
        elif fields[0] == "deposit_ledger":
            result["ledgers"].append(tuple(float(x) for x in fields[1:]))
    return result


def make_case(prefix: pathlib.Path, reverse: bool) -> dict:
    case = json.loads((ROOT / "examples/stage1_case.json").read_text())
    case["case_id"] = "stage1_transport_smoke_reversed" if reverse else "stage1_transport_smoke"
    case["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    case["source"]["spectrum"]["path"] = str(ROOT / "tests/data/mono_100MeV_proton.dat")
    case["source"]["spectrum"]["metadata_path"] = str(
        ROOT / "tests/data/mono_spectrum_metadata.v1.json"
    )
    case["material_overrides"] = [{
        "material_id": "B_enriched_smoke", "revision": "test-1",
        "extends": "B_natural_G4_11_4@1", "status": "project_specific",
        "composition": {
            "basis": "elemental_mass_fraction", "isotopes": "explicit",
            "residual_policy": "exact",
            "components": [
                {"isotope": "B-10", "value": 0.95},
                {"isotope": "B-11", "value": 0.05},
            ],
        },
        "provenance": [{
            "properties": ["composition", "isotopes"],
            "source": "synthetic Stage-1 Geant4 isotope-construction test",
            "version": "1",
        }],
    }]
    case["geometry"]["shield_stack"]["layers"][1]["material"] = "B_enriched_smoke"
    case["quantities"] = ["TID", "H100/10"]
    case["accuracy"]["scope"] = "combined"
    case["accuracy"]["target"] = 0.10
    case["accuracy"]["systematic_budget"] = {
        "status": "incomplete", "known_relative_floor": 0.20
    }
    case["accuracy"]["sampling"] = {
        "min_histories": 20,
        "max_histories": 60,
        "batch_size": 20,
    }
    case["output"]["prefix"] = str(prefix)
    if reverse:
        case["geometry"]["absorbers"].reverse()
    return case


def run(work: pathlib.Path, label: str, reverse: bool) -> dict[str, dict]:
    prefix = work / label
    case_path = work / f"{label}.json"
    case_path.write_text(json.dumps(make_case(prefix, reverse), indent=2) + "\n")
    completed = subprocess.run(
        [str(EXE), f"--case={case_path}"], cwd=work, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(f"transport smoke failed ({completed.returncode})\n{completed.stdout}")
    return {
        absorber: {
            "summary": parse_summary(work / f"{label}__{absorber}_run_summary.dat"),
            "result": json.loads((work / f"{label}__{absorber}_result.json").read_text()),
            "manifest": json.loads((work / f"{label}__{absorber}_resolved_case.json").read_text()),
        }
        for absorber in ("water_proxy", "silicon_proxy")
    }


def main() -> int:
    comparisons: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="shieldsim-transport-") as directory:
        work = pathlib.Path(directory)
        normal = run(work, "normal", False)
        reversed_result = run(work, "reversed", True)

        expected_density = {"water_proxy": 1.0, "silicon_proxy": 2.33}
        expected_mass = {"water_proxy": 0.0025, "silicon_proxy": 0.005825}
        for absorber, bundle in normal.items():
            result = bundle["summary"]
            assert result["events"] == 20
            assert [round(row[2], 6) for row in result["layers"]] == [2.699, 2.34]
            target = result["targets"][0]
            assert math.isclose(target[2], expected_density[absorber], rel_tol=1e-12)
            assert math.isclose(target[3], expected_mass[absorber], rel_tol=1e-12)
            ledger = result["ledgers"][0]
            _, total, ionizing, nonionizing = ledger
            assert math.isclose(total, ionizing + nonionizing, rel_tol=2e-10, abs_tol=1e-30)
            comparisons.extend([
                {
                    "id": f"{absorber}_density",
                    "expected": expected_density[absorber],
                    "actual": target[2],
                    "unit": "g/cm3",
                },
                {
                    "id": f"{absorber}_analytic_slab_mass",
                    "expected": expected_mass[absorber],
                    "actual": target[3],
                    "unit": "kg",
                },
                {
                    "id": f"{absorber}_deposition_ledger_closure",
                    "expected": ionizing + nonionizing,
                    "actual": total,
                    "unit": "MeV",
                },
            ])

            # Resetting the seed and rebuilding an independent world makes an
            # absorber's result exactly invariant to absorber list order.
            other = reversed_result[absorber]["summary"]
            assert result["ledgers"] == other["ledgers"]
            assert result["targets"] == other["targets"]

            result_json = bundle["result"]
            assert result_json["schema_version"] == "shieldSim.result.v1"
            assert result_json["outputs"][0]["validation_status"] == \
                   "implemented_not_physically_validated_stage1"
            assert result_json["outputs"][0]["accuracy"]["status"] == "not_met"
            assert "systematic floor" in result_json["outputs"][0]["accuracy"]["reason"]
            hardness = result_json["outputs"][1]
            assert hardness["name"] == "H100/10"
            assert 0 <= hardness["numerator_count"] <= hardness["denominator_count"] <= 20
            assert "indicator_sample_covariance" in hardness
            if hardness["value"] is not None:
                assert 0.0 <= hardness["interval_95"]["low"] <= hardness["value"] <= 1.0
                assert hardness["value"] <= hardness["interval_95"]["high"] <= 1.0
            manifest = bundle["manifest"]
            for field in ("code_revision", "compiler", "platform", "geant4_version",
                          "geant4_data_libraries", "material_library"):
                assert field in manifest

    if EVIDENCE_DIR is not None:
        EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
        (EVIDENCE_DIR / "stage1_transport_smoke.json").write_text(
            json.dumps(
                {
                    "schema_version": "shieldSim.validation_evidence.v1",
                    "test_case": "stage1_transport_smoke",
                    "comparisons": comparisons,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )

    print("Stage-1 transport smoke and ledger checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

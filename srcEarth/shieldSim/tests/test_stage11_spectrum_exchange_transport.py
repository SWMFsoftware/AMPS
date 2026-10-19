#!/usr/bin/env python3
"""Drive the strict Stage-11 exchange file through Geant4 and proxy output.

The expected source normalization is integrated independently from the JSON
nodes.  Positive-positive segments use the declared power law and the segment
touching zero uses the declared linear interpolation.  This test therefore
detects center/bin reinterpretation, omitted SI-area conversion, endpoint
extrapolation, and a lost isotropic angular factor.
"""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys


if len(sys.argv) != 5:
    raise SystemExit(
        "usage: test_stage11_spectrum_exchange_transport.py "
        "SHIELDSIM SOURCE_ROOT EVIDENCE_DIR RUN_DIR"
    )

EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE = pathlib.Path(sys.argv[3]).resolve()
RUN_DIR = pathlib.Path(sys.argv[4]).resolve()
FIXTURE = ROOT / "tests/data/stage11_spectrum_exchange.v1.json"


def parse_summary(path: pathlib.Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) == 3 and fields[0] in {"scalar", "count"}:
            values[fields[1]] = float(fields[2])
    return values


def invoke(case_path: pathlib.Path, validate_only: bool = False) -> subprocess.CompletedProcess[str]:
    option = "--validate-case" if validate_only else "--case"
    return subprocess.run(
        [str(EXE), f"{option}={case_path}"], cwd=RUN_DIR, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )


def make_case(prefix: pathlib.Path) -> dict:
    case = json.loads(
        (ROOT / "examples/stage11_spectrum_exchange_case.json").read_text(encoding="utf-8")
    )
    case["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    case["source"]["spectrum"]["path"] = str(FIXTURE)
    case["accuracy"]["sampling"] = {
        "min_histories": 256, "max_histories": 256, "batch_size": 256,
    }
    case["output"]["prefix"] = str(prefix)
    return case


def require_rejected(case: dict, label: str, expected: str) -> None:
    path = RUN_DIR / f"rejected_{label}.json"
    path.write_text(json.dumps(case, indent=2) + "\n", encoding="utf-8")
    completed = invoke(path, validate_only=True)
    assert completed.returncode != 0, f"invalid {label} case was accepted"
    assert expected in completed.stdout, completed.stdout


def main() -> int:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    prefix = RUN_DIR / "stage11_exchange"
    case = make_case(prefix)
    case_path = RUN_DIR / "stage11_exchange_case.json"
    case_path.write_text(json.dumps(case, indent=2) + "\n", encoding="utf-8")

    completed = invoke(case_path)
    if completed.returncode != 0:
        raise AssertionError(f"Stage-11 transport failed\n{completed.stdout}")

    stem = pathlib.Path(str(prefix) + "__silicon_proxy")
    summary = parse_summary(pathlib.Path(str(stem) + "_run_summary.dat"))
    result = json.loads(
        pathlib.Path(str(stem) + "_result.json").read_text(encoding="utf-8")
    )
    resolved = json.loads(
        pathlib.Path(str(stem) + "_resolved_case.json").read_text(encoding="utf-8")
    )

    # Segment 10--100 MeV: J=100000/E^2 -> integral 9000 m^-2 s^-1 sr^-1.
    # Segment 100--1000 MeV: linear 10 -> 0 -> integral 4500 in the same unit.
    # shieldSim's source-plane normalization is cm2, hence the exact 1e-4 factor.
    expected_integral_cm2 = (9000.0 + 4500.0) * 1.0e-4
    expected_rate = expected_integral_cm2 * math.pi * 25.0
    assert summary["events"] == 256
    assert summary["input_proton"] == 256 and summary["input_alpha"] == 0
    assert math.isclose(summary["source_norm_no_angular"], expected_integral_cm2,
                        rel_tol=2e-10)
    assert math.isclose(summary["source_angular_factor"], math.pi, rel_tol=2e-6)
    assert math.isclose(summary["incident_particle_rate_s"], expected_rate, rel_tol=2e-6)
    assert result["source_format"] == "spectrum_exchange_v1"
    assert result["source_component_id"] == "local_proton_points"
    assert result["source_science_status"] == "verification_fixture"
    assert result["source_hash"] == resolved["source"]["spectrum_hash"]
    assert result["source_normalization"]["integral_no_angular_unit"] == "1/(cm2 sr s)"
    assert resolved["source"]["metadata"]["context"]["location"]["frame"] == "GSM"
    assert resolved["source"]["metadata"]["quality"]["status"] == "accepted"
    assert resolved["run_status"] == \
        "stage11_exchange_consumer_implemented_and_tested_AMPS_producer_integration_pending"

    beam = json.loads(json.dumps(case))
    beam["source"]["mode"] = "beam"
    require_rejected(beam, "beam", "requires case.source.mode='isotropic'")
    missing = json.loads(json.dumps(case))
    missing["source"]["spectrum"]["component_id"] = "not_present"
    require_rejected(missing, "component", "component_id is absent")
    mission = json.loads(json.dumps(case))
    mission["purpose"] = "preliminary_mission_environment_effect_assessment"
    require_rejected(mission, "fixture_mission", "science-facing spectrum_exchange_v1")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "status": "PASS",
        "boundary": "exchange JSON -> strict case -> exact spectrum sampler -> Geant4 -> proxies",
        "fixture": str(FIXTURE),
        "component_id": "local_proton_points",
        "comparisons": [
            {
                "id": "declared_curve_integral",
                "expected": expected_integral_cm2,
                "actual": summary["source_norm_no_angular"],
                "unit": "1/(cm2 sr s)",
                "method": "independent analytic power-law plus linear-zero integration",
            },
            {
                "id": "isotropic_plane_rate",
                "expected": expected_rate,
                "actual": summary["incident_particle_rate_s"],
                "unit": "1/s",
                "method": "integral times pi times 25 cm2 source plane",
            },
        ],
    }
    (EVIDENCE / "stage11_spectrum_exchange_transport.json").write_text(
        json.dumps(evidence, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    print("PASS Stage-11 spectrum exchange drove the full Geant4/proxy path")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

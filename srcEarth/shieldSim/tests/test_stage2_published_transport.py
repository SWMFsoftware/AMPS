#!/usr/bin/env python3
"""End-to-end Geant4 test driven directly by a downloaded AMS-02 spectrum.

This test closes a boundary that adapter tests cannot: the checksum-pinned
published proton bins must pass through CaseDefinition, primary sampling,
Geant4 transport, and the production TID/H100/10 writers.  The expected source
integral is calculated independently from the pinned rigidity edges and CSV
values.  The test does not claim that one near-Earth Bartels rotation is a
complete mission GCR environment or that its transported dose matches an
observation.
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import subprocess
import sys
import tempfile


if len(sys.argv) not in (3, 4):
    raise SystemExit(
        "usage: test_stage2_published_transport.py SHIELDSIM SOURCE_ROOT [EVIDENCE_DIR]"
    )
EXE = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
EVIDENCE_DIR = pathlib.Path(sys.argv[3]).resolve() if len(sys.argv) == 4 else None
FIXTURE = ROOT / "tests/reference_spectra/gcr"


def parse_summary(path: pathlib.Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        fields = raw.split()
        if len(fields) == 3 and fields[0] in {"scalar", "count"}:
            values[fields[1]] = float(fields[2])
    return values


def expected_proton_integral() -> float:
    """Reproduce sum(J_i * Delta E_i) without using production adapter output."""
    interpretation = json.loads(
        (FIXTURE / "ams02_20131118_interpretation.v1.json").read_text(encoding="utf-8")
    )
    rigidity = interpretation["proton_rigidity_edges_GV"]
    mass = interpretation["nucleon_rest_mass_GeV"]

    def kinetic_mev(rigidity_gv: float) -> float:
        return 1000.0 * (math.sqrt(rigidity_gv**2 + mass**2) - mass)

    edges = [kinetic_mev(value) for value in rigidity]
    with (FIXTURE / "ams02_20131118_hapi_data.csv").open(newline="", encoding="utf-8") as stream:
        row = next(csv.reader(line for line in stream if not line.startswith("#")))
    # CSV fields 2..46 are the 45 proton flux bins. Convert
    # m^-2 GeV^-1 to cm^-2 MeV^-1 with the documented factor 1e-7.
    flux = [float(value) * 1.0e-7 for value in row[2:47]]
    assert len(flux) == len(edges) - 1 == 45
    return sum(value * (upper - lower) for value, lower, upper in zip(flux, edges, edges[1:]))


def make_case(prefix: pathlib.Path) -> dict:
    case = json.loads((ROOT / "examples/stage1_case.json").read_text(encoding="utf-8"))
    case["case_id"] = "stage2_downloaded_ams02_proton_transport"
    case["purpose"] = "verification_test"
    case["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    case["geometry"]["absorbers"] = [case["geometry"]["absorbers"][1]]
    case["source"] = {
        "mode": "isotropic",
        "spectrum": {
            "format": "ams02_hapi_component",
            "info_path": str(FIXTURE / "ams02_20131118_hapi_info.json"),
            "data_path": str(FIXTURE / "ams02_20131118_hapi_data.csv"),
            "interpretation_path": str(FIXTURE / "ams02_20131118_interpretation.v1.json"),
            "component_id": "proton",
        },
        "energy_limits": {
            "proton": {
                "min": {"value": 400.0, "unit": "MeV"},
                "max": {"value": 100.0, "unit": "GeV"},
            },
            # Required by the v1 case schema but unused because the admitted
            # published component contains no alpha population.
            "alpha": {
                "min": {"value": 1.0, "unit": "MeV"},
                "max": {"value": 100.0, "unit": "GeV"},
            },
        },
    }
    case["accuracy"]["scope"] = "combined"
    case["accuracy"]["target"] = 0.10
    case["accuracy"]["systematic_budget"] = {
        "status": "incomplete",
        "known_relative_floor": 0.20,
    }
    # A fixed history count makes this an exact deterministic regression while
    # the deliberately unachievable combined target prevents early stopping.
    case["accuracy"]["sampling"] = {
        "min_histories": 256,
        "max_histories": 256,
        "batch_size": 256,
    }
    case["output"]["prefix"] = str(prefix)
    return case


def invoke(case_path: pathlib.Path, validate_only: bool = False) -> subprocess.CompletedProcess[str]:
    option = "--validate-case" if validate_only else "--case"
    return subprocess.run(
        [str(EXE), f"{option}={case_path}"],
        cwd=case_path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def require_rejected(work: pathlib.Path, base: dict, label: str, expected: str) -> None:
    path = work / f"rejected_{label}.json"
    path.write_text(json.dumps(base, indent=2) + "\n", encoding="utf-8")
    completed = invoke(path, validate_only=True)
    assert completed.returncode != 0, f"invalid {label} case was accepted"
    assert expected in completed.stdout, completed.stdout


def main() -> int:
    expected_norm = expected_proton_integral()
    observed_norm = math.nan
    with tempfile.TemporaryDirectory(prefix="shieldsim-stage2-transport-") as directory:
        work = pathlib.Path(directory)
        prefix = work / "published"
        case = make_case(prefix)
        case_path = work / "published_case.json"
        case_path.write_text(json.dumps(case, indent=2) + "\n", encoding="utf-8")
        completed = invoke(case_path)
        if completed.returncode != 0:
            raise AssertionError(
                f"published-spectrum transport failed ({completed.returncode})\n{completed.stdout}"
            )

        stem = work / "published__silicon_proxy"
        summary = parse_summary(pathlib.Path(str(stem) + "_run_summary.dat"))
        result = json.loads(pathlib.Path(str(stem) + "_result.json").read_text(encoding="utf-8"))
        resolved = json.loads(
            pathlib.Path(str(stem) + "_resolved_case.json").read_text(encoding="utf-8")
        )

        assert summary["events"] == 256
        assert summary["input_proton"] == 256
        assert summary["input_alpha"] == 0
        observed_norm = summary["source_norm_no_angular"]
        assert math.isclose(summary["source_norm_no_angular"], expected_norm, rel_tol=2e-6)
        assert math.isclose(summary["source_angular_factor"], math.pi, rel_tol=2e-6)
        assert math.isclose(summary["source_plane_area_cm2"], 25.0, rel_tol=1e-12)
        expected_rate = expected_norm * math.pi * 25.0
        assert math.isclose(summary["incident_particle_rate_s"], expected_rate, rel_tol=2e-6)

        assert result["events"] == 256
        assert result["source_format"] == "ams02_hapi_component"
        assert result["source_component_id"] == "proton"
        assert result["source_id"] == "AMS02_BR2460_2013-11-18"
        assert result["source_hash"] == resolved["source"]["spectrum_hash"]
        normalization = result["source_normalization"]
        assert normalization["integral_no_angular_unit"] == "1/(cm2 sr s)"
        assert normalization["incident_particle_rate_unit"] == "1/s"
        assert math.isclose(normalization["integral_no_angular"], expected_norm, rel_tol=1e-12)
        assert math.isclose(normalization["angular_factor"], math.pi, rel_tol=1e-12)
        assert math.isclose(normalization["source_plane_area_cm2"], 25.0, rel_tol=1e-12)
        assert math.isclose(normalization["incident_particle_rate"], expected_rate, rel_tol=1e-12)
        assert resolved["source"]["transport_component_id"] == "proton"
        assert resolved["run_status"] == \
            "stage2_published_spectrum_transport_test_not_physical_validation"

        outputs = {item["name"]: item for item in result["outputs"]}
        tid = outputs["TID"]
        assert tid["validation_status"] == \
            "stage2_published_spectrum_transport_test_not_physical_validation"
        assert tid["value"] > 0.0
        assert tid["source_normalized_value"] > 0.0
        assert tid["source_normalized_unit"] == "Gy/s"
        assert math.isclose(
            tid["source_normalized_value"], tid["value"] * expected_rate, rel_tol=3e-6
        )
        hardness = outputs["H100/10"]
        assert hardness["validation_status"] == tid["validation_status"]
        assert 0 <= hardness["numerator_count"] <= hardness["denominator_count"] <= 256

        helium = json.loads(json.dumps(case))
        helium["source"]["spectrum"]["component_id"] = "helium"
        require_rejected(work, helium, "inclusive_helium", "admits only the explicit A=1 proton")
        mission = json.loads(json.dumps(case))
        mission["purpose"] = "preliminary_mission_environment_effect_assessment"
        require_rejected(work, mission, "mission_use", "admitted only for verification_test")
        beam = json.loads(json.dumps(case))
        beam["source"]["mode"] = "beam"
        require_rejected(work, beam, "beam_mapping", "requires case.source.mode='isotropic'")

    if EVIDENCE_DIR is not None:
        EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
        evidence = {
            "schema_version": "shieldSim.validation_evidence.v1",
            "test_case": "stage2_downloaded_ams02_proton_transport",
            "boundary": "published bins -> case -> source -> Geant4 -> TID/H100/10",
            "comparisons": [
                {
                    "id": "published_source_integral",
                    "expected": expected_norm,
                    "actual": observed_norm,
                    "unit": "1/(cm2 sr s)",
                    "method": "independent rigidity transform and rectangular bin integration",
                }
            ],
        }
        (EVIDENCE_DIR / "stage2_published_transport.json").write_text(
            json.dumps(evidence, sort_keys=True, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )

    print("PASS downloaded AMS-02 proton spectrum drove full Geant4/proxy path")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

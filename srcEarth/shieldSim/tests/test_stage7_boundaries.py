#!/usr/bin/env python3
"""Strict Stage-7 CLI, time-series, quadrature, and allocation boundaries.

All generated documents live below the caller-owned build directory.  The
test uses only synthetic values: published AMS/GOES qualification boundaries
are exercised separately and no synthetic value is presented as physics.
"""

from __future__ import annotations

import copy
import json
import math
import pathlib
import shutil
import subprocess
import sys


if len(sys.argv) != 5:
    raise SystemExit(
        "usage: test_stage7_boundaries.py EXPOSURE_TOOL SOURCE_ROOT RUN_DIR EVIDENCE_DIR"
    )
TOOL = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(sys.argv[2]).resolve()
RUN_DIR = pathlib.Path(sys.argv[3]).resolve()
EVIDENCE_DIR = pathlib.Path(sys.argv[4]).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def write(path: pathlib.Path, document: dict) -> None:
    path.write_text(json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n",
                    encoding="utf-8")


def run(command: str, document: dict, name: str, success: bool = True) -> tuple[subprocess.CompletedProcess[str], dict | None]:
    input_path = RUN_DIR / f"{name}.input.json"
    output_path = RUN_DIR / f"{name}.output.json"
    write(input_path, document)
    completed = subprocess.run(
        [str(TOOL), command, str(input_path), str(output_path)],
        text=True, capture_output=True, check=False,
    )
    (RUN_DIR / f"{name}.log").write_text(completed.stdout + completed.stderr, encoding="utf-8")
    if success:
        require(completed.returncode == 0, f"{name} failed: {completed.stderr}")
        require(output_path.is_file(), f"{name} did not write output")
        return completed, json.loads(output_path.read_text(encoding="utf-8"))
    require(completed.returncode != 0, f"{name} unexpectedly succeeded")
    return completed, None


def main() -> int:
    if RUN_DIR.exists():
        shutil.rmtree(RUN_DIR)
    RUN_DIR.mkdir(parents=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    example = json.loads(
        (ROOT / "examples/stage7_time_resolved_exposure.json").read_text(encoding="utf-8")
    )
    _, first = run("integrate", example, "valid_a")
    _, second = run("integrate", example, "valid_b")
    require(first == second, "Stage-7 integration is not deterministic")
    provenance = first.get("build_provenance", {})
    require(all(provenance.get(field) for field in
                ("code_revision", "source_tree_state", "compiler", "platform")),
            "Stage-7 CLI result lacks build/source provenance")
    exposure = first["integrated_exposures"][0]
    # 2 units/s for 300 s plus 1 unit/s for 600 s; the explicit 300 s gap
    # contributes exactly zero and is never interpolated.
    require(exposure["value"] == 1200.0, "example explicit-gap integral")
    require(exposure["standard_uncertainty"] is None and
            exposure["accuracy_status"] == "not_assessable",
            "incomplete uncertainty acquired a false accuracy claim")
    require(len(first["configuration_segments"]) == 2, "configuration split was averaged away")

    unknown = copy.deepcopy(example)
    unknown["typo"] = True
    failed, _ = run("integrate", unknown, "unknown_field", False)
    require("unknown field" in failed.stderr, "unknown Stage-7 field diagnosis")

    local_time = copy.deepcopy(example)
    local_time["intervals"][0]["start_utc"] = "2017-09-10T11:25:00-05:00"
    failed, _ = run("integrate", local_time, "local_time", False)
    require("YYYY-MM-DD" in failed.stderr, "local timestamp was not rejected")

    hidden_gap = copy.deepcopy(example)
    hidden_gap["intervals"].pop(1)
    failed, _ = run("integrate", hidden_gap, "hidden_gap", False)
    require("explicit gap" in failed.stderr, "hidden gap was not rejected")

    false_gcr = copy.deepcopy(example)
    false_gcr["source_qualification"]["environment"] = "ams02_a1_proton_component"
    false_gcr["source_qualification"]["environment_completeness"] = "complete_gcr"
    failed, _ = run("integrate", false_gcr, "false_full_gcr", False)
    require("incomplete_not_full_gcr" in failed.stderr, "false full-GCR label was not rejected")

    false_goes = copy.deepcopy(example)
    false_goes["source_qualification"].update({
        "environment": "goes16_sgps_directional_ingestion",
        "environment_completeness": "directional_sep",
        "transport_admission": "verification_test_only",
    })
    failed, _ = run("integrate", false_goes, "false_goes_transport", False)
    require("GOES transport" in failed.stderr, "GOES without boundary mapping was admitted")

    quadrature = {
        "schema_version": "shieldSim.response_quadrature.v1",
        "source_bin_edges_MeV": [1.0, 4.0],
        "differential_source": [2.0],
        "response_energy_nodes_MeV": [1.0, 2.0, 4.0],
        "response_values": [3.0, 3.0 + 0.5 * math.log(2.0), 3.0 + 0.5 * math.log(4.0)],
        "source_density_unit": "cm-2 MeV-1",
        "response_unit": "proxy cm2",
        "output_unit": "proxy",
        "relative_tolerance": 1.0e-6,
        "maximum_refinements": 16,
        "extrapolation_policy": "refuse",
    }
    _, folded = run("fold-spectrum", quadrature, "quadrature")
    expected_fold = 2.0 * (3.0 * 3.0 + 0.5 * ((4.0 * math.log(4.0) - 4.0) + 1.0))
    require(folded["converged"] and folded["extrapolation_used"] is False,
            "adaptive quadrature convergence/provenance")
    require(abs(folded["value"] - expected_fold) / expected_fold < 3.0e-6,
            "adaptive quadrature analytic result")
    extrapolated = copy.deepcopy(quadrature)
    extrapolated["source_bin_edges_MeV"][0] = 0.5
    failed, _ = run("fold-spectrum", extrapolated, "quadrature_extrapolation", False)
    require("extrapolation" in failed.stderr, "quadrature extrapolation was not refused")

    allocation_input = {
        "schema_version": "shieldSim.neyman_allocation.v1",
        "policy": "pilot_neyman_no_transport_biasing",
        "total_histories": 100,
        "minimum_per_stratum": 10,
        "strata": [
            {"stratum_id": "common", "source_weight": 0.5, "pilot_standard_deviation": 1.0},
            {"stratum_id": "rare", "source_weight": 0.3, "pilot_standard_deviation": 2.0},
            {"stratum_id": "zero_pilot", "source_weight": 0.2, "pilot_standard_deviation": 0.0},
        ],
    }
    _, allocation = run("allocate", allocation_input, "allocation")
    histories = [int(item["histories"]) for item in allocation["allocation"]]
    require(histories == [42, 48, 10] and sum(histories) == 100,
            "Neyman exact-total/minimum allocation")
    require(allocation["physics_biasing"] is False, "allocation silently enabled physics biasing")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_case": "stage7_boundaries",
        "validation_status": "stage7_software_framework_incomplete_physical_gcr_goes_qualification",
        "comparisons": [
            {"id": "explicit_gap_time_integral", "expected": 1200.0,
             "actual": exposure["value"], "unit": "verification_unit"},
            {"id": "adaptive_response_quadrature", "expected": expected_fold,
             "actual": folded["value"], "unit": "proxy"},
            {"id": "neyman_exact_history_total", "expected": 100.0,
             "actual": float(sum(histories)), "unit": "histories"},
        ],
        "rejections": ["unknown_field", "local_time", "hidden_gap", "false_full_gcr",
                       "false_goes_transport", "quadrature_extrapolation"],
    }
    write(EVIDENCE_DIR / "stage7_boundaries.json", evidence)
    print("PASS Stage-7 strict time, gap, qualification, quadrature, and allocation boundaries")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)

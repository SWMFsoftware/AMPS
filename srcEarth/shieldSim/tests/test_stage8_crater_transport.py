#!/usr/bin/env python3
"""End-to-end smoke test for the foreground-only axial CRaTER benchmark."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys

if len(sys.argv) != 7:
    raise SystemExit(
        "usage: test_stage8_crater_transport.py BENCHMARK PARALLEL_RUNNER SOURCE_ROOT "
        "EVIDENCE_DIR RUN_DIR PYTHON"
    )

BENCHMARK = pathlib.Path(sys.argv[1]).resolve()
RUNNER = pathlib.Path(sys.argv[2]).resolve()
ROOT = pathlib.Path(sys.argv[3]).resolve()
EVIDENCE = pathlib.Path(sys.argv[4]).resolve()
RUN_DIR = pathlib.Path(sys.argv[5]).resolve()
PYTHON = pathlib.Path(sys.argv[6]).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    case = json.loads(
        (ROOT / "examples/stage8_crater_foreground_axial.json").read_text(encoding="utf-8")
    )
    case["case_id"] = "stage8_crater_transport_smoke"
    # Twenty histories per elemental charge exercises the complete random
    # source path in expectation while keeping the routine CTest gate bounded.  This run is
    # intentionally not expected to meet the 10% per-bin target; production
    # comparison runs use the larger example history count or more.
    case["histories"] = 560
    case["random_seed"] = 82560
    case["source"]["path"] = str(
        ROOT / "tests/reference_crater/glace/solmin_gcr.json"
    )
    case["observation"]["path"] = str(
        ROOT / "tests/reference_spectra/observations/crater_gcr_let_spectra_v1_3.csv"
    )
    case["output"]["directory"] = str(RUN_DIR)
    case_path = RUN_DIR / "case.json"
    case_path.write_text(json.dumps(case, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    completed = subprocess.run(
        [str(PYTHON), str(RUNNER), "--benchmark", str(BENCHMARK),
         "--case", str(case_path), "--workers", "2", "--python", str(PYTHON)],
        cwd=RUN_DIR, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    (RUN_DIR / "run.log").write_text(completed.stdout, encoding="utf-8")
    require(completed.returncode == 0, f"CRaTER transport failed; see {RUN_DIR/'run.log'}")
    result_path = RUN_DIR / "stage8_crater_transport_smoke__foreground_only_axial_result.json"
    require(result_path.is_file(), "CRaTER result JSON is missing")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    require(result["benchmark_scope"] == "foreground_only_axial", "scope label missing")
    require("not_stage8" in result["validation_status"], "non-validation status missing")
    require(result["result_type"] == "complete_parallel_merge", "shards were not merged")
    require(result["parallel"]["merge_status"] == "complete_all_shards_present" and
            result["parallel"]["shard_count"] == 2,
            "parallel completeness provenance is missing")
    require(len(result["geometry"]["axial_elements"]) == 19, "axial stack changed")
    require(len(result["primary_ions"]) == 28, "H-Ni mapping changed")
    require(len(result["bins"]) == 291, "observation-aligned histogram changed")
    require(result["event_counts"]["coincidence"] > 0, "no telescope coincidences")
    diagnostic = result["diagnostic_comparison"]
    ideal = float(diagnostic["ideal_integral_per_cm2_sr_day"])
    cut = float(diagnostic["after_cuts_integral_per_cm2_sr_day"])
    observed = float(diagnostic["observed_integral_per_cm2_sr_day"])
    require(math.isfinite(ideal) and ideal > 0.0, "ideal integral is unavailable")
    require(math.isfinite(cut) and 0.0 <= cut <= ideal, "cut spectrum is not a subset")
    require(math.isfinite(observed) and observed > 0.0, "observation was not loaded")
    for row in result["bins"]:
        require(row["ideal_foreground_per_cm2_sr_day_keV_per_um"] >= 0.0,
                "negative ideal bin")
        require(row["foreground_after_cuts_per_cm2_sr_day_keV_per_um"] >= 0.0,
                "negative cut bin")

    figure_prefix = RUN_DIR / "stage8_crater_transport_smoke__foreground_only_axial_comparison"
    for extension in ("png", "svg"):
        figure = figure_prefix.with_suffix(f".{extension}")
        require(figure.is_file() and figure.stat().st_size > 1000,
                f"missing/empty comparison figure {figure}")

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_case": "stage8_crater_foreground_axial_transport_smoke",
        "benchmark_scope": "foreground_only_axial",
        "scientific_status": "transport_and_plot_smoke_not_observational_acceptance",
        "comparisons": [
            {
                "id": "published_D1_D6_geometry_factor_provenance",
                "unit": "cm2 sr",
                "expected": 0.569,
                "actual": result["geometry"]["D1_D6_geometry_factor_cm2_sr"],
            },
            {
                "id": "transport_history_accounting",
                "unit": "histories",
                "expected": float(case["histories"]),
                "actual": float(result["transport"]["histories"]),
            },
        ],
        "diagnostic_only": {
            "ideal_integral_per_cm2_sr_day": ideal,
            "after_cuts_integral_per_cm2_sr_day": cut,
            "observed_integral_per_cm2_sr_day": observed,
            "acceptance_decision": "not_defined",
        },
        "figures": [
            figure_prefix.with_suffix(".png").name,
            figure_prefix.with_suffix(".svg").name,
        ],
    }
    (EVIDENCE / "stage8_crater_transport.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("Foreground-only axial CRaTER transport smoke test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

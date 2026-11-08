#!/usr/bin/env python3
"""Import the complete dependency-light Enceladus forward-model campaign.

The application-local runner remains the scientific authority: it compiles
``ForwardModel.cpp``, checks the source/runtime grammars and chemistry/source
laws, and compares 31 calculated values with the independently maintained
reference JSON.  This adapter only supplies the globally configured compiler,
validates the exact local registry, and translates each outcome into the
versioned core-harness protocol.  It deliberately does not claim that this
dependency-light campaign is a native AMPS trajectory qualification.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Dict, Mapping


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, load_context, main_guard,
    run_command, write_result)


LOCAL_CASES = [
    "U00_exporter", "U01_source_export", "U02_core_compile",
    "U03b_canonical_source_geometry", "U03_source_parser",
    "U04_runtime_parser", "U04b_transport_only_parser", "U05_chemistry",
    "U06_source_law", "U07_numerics", "U08_validation",
]

# These names are copied from the reviewed, hand-maintained reference file,
# not discovered after execution.  Consequently, adding or removing a local
# scalar fails the global registry check until both test owners are reviewed.
REFERENCE_KEYS = [
    "minimum_cell_size_m", "global_time_step_s", "primary_weight_H2O",
    "primary_weight_CO2", "absolute_rate_H2O_s-1",
    "absolute_rate_CO2_s-1", "absolute_primary_weight_H2O",
    "planned_daughter_species", "photo_h2o_oh_h_quiet_s-1",
    "photo_h2o_oh_h_active_s-1", "cx_hplus_h2o_sigma_20kms_m2",
    "cx_hplus_h2o_rate_20kms_s-1", "cx_zero_reactant_rate_s-1",
    "tenishev_k3_c-0.35_s4.8", "tenishev_q73_c-0.35_s4.8",
    "stripe_weight_sum", "source_count", "point_normal_norm",
    "point_axis_norm", "stripe_minimum_radius_m", "stripe_maximum_radius_m",
    "stripe_reference_radius_m", "parsed_sample_length_steps",
    "parsed_total_iterations", "parsed_flyby_refinement_target_cell_size_m",
    "parsed_flyby_refinement_radius_m", "parsed_observation_requests",
    "parsed_split_position_cell_diagonal_fraction",
    "parsed_split_velocity_thermal_speed_fraction",
    "parsed_nonuniform_split_position_shift", "negative_validation_checks",
]


def normalize(value: str) -> str:
    """Convert a local mixed-case/punctuation identifier to manifest grammar."""

    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", value)).strip("-").upper()


def global_local_id(local_id: str) -> str:
    return "ENCELADUS-{}".format(normalize(local_id))


def global_reference_id(key: str) -> str:
    return "ENCELADUS-REF-{}".format(normalize(key))


EXPECTED_CASES = (
    [global_local_id(value) for value in LOCAL_CASES]
    + [global_reference_id(value) for value in REFERENCE_KEYS]
)


def load_receipt(path: Path) -> Mapping[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RunnerError("cannot read Enceladus receipt {}: {}".format(path, error))
    if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
        raise RunnerError("Enceladus receipt is not an object with a results array")
    return payload


def main() -> int:
    context = load_context("Run the Enceladus forward-model verification campaign")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    if (not requested or len(set(requested)) != len(requested) or
            any(identifier not in EXPECTED_CASES for identifier in requested)):
        raise RunnerError("invalid Enceladus case selection")

    repo = Path(str(context["repo_root"]))
    artifacts = Path(str(context["artifact_dir"]))
    artifacts.mkdir(parents=True, exist_ok=True)
    local_output = artifacts / "local-campaign"
    transcript = artifacts / "enceladus-forward-model.txt"
    timeout = int(context["test"].get("timeout_seconds", 600))

    argv = command(context, "python") + [
        str(repo / "srcEnceladus" / "test" / "run_tests.py"),
        "--output", str(local_output), "--cxx",
    ] + command(context, "cxx")
    # The local campaign is intentionally one compiled fixture, so it cannot
    # be split safely into 42 subprocesses.  Announce every selected identity
    # before entering that fixture; the exact post-run receipt below remains
    # the authority for which assertions actually completed.
    for identifier in requested:
        print("[CASE] RUN  {} application-local forward assertion".format(
            identifier), flush=True)
    completed = run_command(
        argv, repo, timeout=timeout, env={"PYTHONDONTWRITEBYTECODE": "1"})
    transcript.write_text(completed.stdout or "", encoding="utf-8")

    receipt_path = local_output / "test_receipt.json"
    payload = load_receipt(receipt_path)
    rows = payload["results"]
    observed_names = [str(row.get("test")) for row in rows if isinstance(row, dict)]
    expected_names = LOCAL_CASES + ["R_{}".format(key) for key in REFERENCE_KEYS]
    if observed_names != expected_names or len(observed_names) != len(set(observed_names)):
        raise RunnerError(
            "Enceladus local registry mismatch: expected {} observed {}".format(
                expected_names, observed_names))

    translated: Dict[str, Dict[str, str]] = {}
    for row in rows:
        local_name = str(row["test"])
        identifier = (global_reference_id(local_name[2:])
                      if local_name.startswith("R_")
                      else global_local_id(local_name))
        local_status = str(row.get("status"))
        if local_status not in {"PASS", "FAIL"}:
            raise RunnerError("unsupported Enceladus status {!r}".format(local_status))
        detail = " ".join(str(row.get("detail", "")).split())
        translated[identifier] = case(
            identifier, local_status,
            detail or "application-local {} result".format(local_name))

    selected_rows = [translated[identifier] for identifier in requested]
    suite_status = "FAIL" if any(row["status"] == "FAIL" for row in selected_rows) else "PASS"
    # A nonzero local exit is required when any imported assertion failed and
    # forbidden when all assertions passed.  Checking both data and process
    # status prevents a stale/partial receipt from being accepted.
    expected_exit = 0 if all(row["status"] == "PASS" for row in translated.values()) else 1
    if completed.returncode != expected_exit:
        raise RunnerError(
            "Enceladus exit {} disagrees with imported status (expected {})".format(
                completed.returncode, expected_exit))

    evidence = [str(receipt_path), str(transcript)]
    comparison = local_output / "reference_comparison.svg"
    if comparison.is_file():
        evidence.append(str(comparison))
    return write_result(
        context, suite_status,
        "imported {} exact Enceladus outcomes".format(len(selected_rows)),
        cases=selected_rows, evidence=evidence,
        metrics={"local_case_count": len(rows), "selected_case_count": len(selected_rows)})


if __name__ == "__main__":
    main_guard(main)

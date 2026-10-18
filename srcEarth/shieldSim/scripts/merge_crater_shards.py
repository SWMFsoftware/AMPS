#!/usr/bin/env python3
"""Merge complete CRaTER shard raw moments into one statistically valid result.

This tool does not average already-normalized spectra.  It first verifies that
every expected shard is present and scientifically identical, then adds the raw
importance-weight sums and sums of squares.  Rates and standard errors are
recomputed once from the global history count, which is the same estimator used
by a serial run.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import pathlib

SECONDS_PER_DAY = 86400.0
CM2_PER_M2 = 10000.0


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def standard_error(total: float, squares: float, histories: int) -> float:
    if histories < 2:
        return 0.0
    centered = max(0.0, squares - total * total / histories)
    return math.sqrt(centered / (histories * (histories - 1)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("shards", nargs="+", type=pathlib.Path)
    args = parser.parse_args()
    require(len(args.shards) > 1, "parallel merge requires at least two shard files")
    documents = [json.loads(path.read_text(encoding="utf-8")) for path in args.shards]
    first = documents[0]
    require(first.get("schema_version") == "shieldSim.crater_foreground_axial_result.v1",
            "unexpected shard schema")
    require(first.get("benchmark_scope") == "foreground_only_axial", "scope label missing")

    expected_count = int(first["parallel"]["shard_count"])
    global_requested = int(first["parallel"]["global_requested_histories"])
    base_seed = int(first["parallel"]["base_random_seed"])
    require(expected_count == len(documents),
            f"expected {expected_count} shards but received {len(documents)}")
    stable_fields = (
        "schema_version", "case_id", "benchmark_scope", "observable",
        "comparison_units", "source", "observation", "geometry", "primary_ions",
        "event_selection", "limitations",
    )
    transport_fields = (
        "geant4_version", "physics_list", "em_physics", "production_cut_um",
        "code_revision", "source_tree_state", "global_requested_histories",
    )
    indexes: set[int] = set()
    seeds: list[int] = []
    histories = 0
    for path, document in zip(args.shards, documents):
        require(document.get("result_type") == "partial_shard",
                f"{path} is not a partial shard")
        for field in stable_fields:
            require(canonical(document[field]) == canonical(first[field]),
                    f"scientific field {field} differs in {path}")
        for field in transport_fields:
            require(canonical(document["transport"][field]) ==
                    canonical(first["transport"][field]),
                    f"transport field {field} differs in {path}")
        parallel = document["parallel"]
        require(int(parallel["shard_count"]) == expected_count and
                int(parallel["global_requested_histories"]) == global_requested and
                int(parallel["base_random_seed"]) == base_seed,
                f"parallel definition differs in {path}")
        index = int(parallel["shard_index"])
        require(index not in indexes, f"duplicate shard index {index}")
        indexes.add(index)
        seeds.append(int(parallel["effective_random_seed"]))
        local = int(parallel["local_histories"])
        require(local == int(document["transport"]["histories"]),
                f"local history accounting differs in {path}")
        histories += local
    require(indexes == set(range(expected_count)), "shard index set is incomplete")
    require(histories == global_requested,
            f"merged histories {histories} != requested {global_requested}")
    require(len(set(seeds)) == expected_count, "shards did not use distinct random seeds")

    bin_count = len(first["bins"])
    require(all(len(document["bins"]) == bin_count for document in documents),
            "shard histogram lengths differ")
    merged_bins: list[dict] = []
    ideal_integral = 0.0
    cut_integral = 0.0
    observed_integral = 0.0
    assessed = 0
    meeting = 0
    log_residuals: list[float] = []
    target = float(first["accuracy"]["requested_relative_statistical_uncertainty"])
    for index in range(bin_count):
        rows = [document["bins"][index] for document in documents]
        for field in ("lower_keV_per_um", "upper_keV_per_um", "center_keV_per_um",
                      "observed_per_cm2_sr_day_keV_per_um"):
            require(all(row[field] == rows[0][field] for row in rows[1:]),
                    f"bin metadata {field} differs at row {index}")
        lower = float(rows[0]["lower_keV_per_um"])
        upper = float(rows[0]["upper_keV_per_um"])
        center = float(rows[0]["center_keV_per_um"])
        observed = float(rows[0]["observed_per_cm2_sr_day_keV_per_um"])
        width = upper - lower
        require(width > 0.0, f"nonpositive bin width at row {index}")
        ideal_sum = math.fsum(float(row["ideal_raw_weight_sum"]) for row in rows)
        ideal_squares = math.fsum(
            float(row["ideal_raw_weight_sum_squares"]) for row in rows
        )
        cut_sum = math.fsum(float(row["cuts_raw_weight_sum"]) for row in rows)
        cut_squares = math.fsum(
            float(row["cuts_raw_weight_sum_squares"]) for row in rows
        )
        conversion = SECONDS_PER_DAY / (CM2_PER_M2 * width)
        ideal = ideal_sum / histories * conversion
        cut = cut_sum / histories * conversion
        ideal_error = standard_error(ideal_sum, ideal_squares, histories) * conversion
        cut_error = standard_error(cut_sum, cut_squares, histories) * conversion
        ideal_integral += ideal * width
        cut_integral += cut * width
        observed_integral += observed * width
        if 0.5 <= center <= 10.0 and cut > 0.0:
            assessed += 1
            if cut_error / cut <= target:
                meeting += 1
        if cut > 0.0 and observed > 0.0:
            log_residuals.append(math.log(cut / observed))
        merged_bins.append({
            "lower_keV_per_um": lower,
            "upper_keV_per_um": upper,
            "center_keV_per_um": center,
            "ideal_foreground_per_cm2_sr_day_keV_per_um": ideal,
            "ideal_mc_standard_error": ideal_error,
            "foreground_after_cuts_per_cm2_sr_day_keV_per_um": cut,
            "cuts_mc_standard_error": cut_error,
            "observed_per_cm2_sr_day_keV_per_um": observed,
            "ideal_raw_weight_sum": ideal_sum,
            "ideal_raw_weight_sum_squares": ideal_squares,
            "cuts_raw_weight_sum": cut_sum,
            "cuts_raw_weight_sum_squares": cut_squares,
        })

    result = copy.deepcopy(first)
    result["result_type"] = "complete_parallel_merge"
    result["validation_status"] = (
        "diagnostic_parallel_benchmark_not_stage8_acceptance_or_full_instrument_validation"
    )
    result["parallel"] = {
        "shard_count": expected_count,
        "merged_shard_indexes": sorted(indexes),
        "local_histories_by_shard": [
            int(document["parallel"]["local_histories"])
            for document in sorted(documents, key=lambda item: int(item["parallel"]["shard_index"]))
        ],
        "global_requested_histories": global_requested,
        "merged_histories": histories,
        "base_random_seed": base_seed,
        "effective_random_seeds": sorted(seeds),
        "seed_strategy": "base_random_seed_plus_shard_index_distinct_Geant4_streams",
        "merge_status": "complete_all_shards_present",
    }
    result["transport"]["histories"] = histories
    result["transport"]["random_seed"] = base_seed
    result["transport"]["parallel_seed_count"] = expected_count
    result["event_counts"] = {
        key: sum(int(document["event_counts"][key]) for document in documents)
        for key in first["event_counts"]
    }
    result["accuracy"] = {
        "requested_relative_statistical_uncertainty": target,
        "assessment_range_keV_per_um": "0.5_to_10",
        "nonzero_bins_assessed": assessed,
        "bins_meeting_target": meeting,
        "status": (
            "requested_MC_precision_met_in_all_nonzero_assessed_bins"
            if assessed > 0 and meeting == assessed
            else "requested_MC_precision_not_met_in_all_nonzero_assessed_bins"
        ),
    }
    geometry_factor = float(result["geometry"]["D1_D6_geometry_factor_cm2_sr"])
    result["diagnostic_comparison"] = {
        "ideal_integral_per_cm2_sr_day": ideal_integral,
        "after_cuts_integral_per_cm2_sr_day": cut_integral,
        "observed_integral_per_cm2_sr_day": observed_integral,
        "ideal_telescope_rate_per_day": ideal_integral * geometry_factor,
        "after_cuts_telescope_rate_per_day": cut_integral * geometry_factor,
        "published_geometry_factor_m2_sr": geometry_factor / CM2_PER_M2,
        "log_space_RMSE_after_cuts_vs_observation": (
            math.sqrt(math.fsum(value * value for value in log_residuals) / len(log_residuals))
            if log_residuals else None
        ),
        "acceptance_decision": (
            "not_defined_no_observational_covariance_or_preregistered_threshold"
        ),
    }
    result["bins"] = merged_bins

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n",
                           encoding="utf-8")
    csv_path = pathlib.Path(str(args.output).replace("_result.json", "_spectrum.csv"))
    require(csv_path != args.output, "output filename must end in _result.json")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "lower_keV_per_um", "upper_keV_per_um", "center_keV_per_um",
            "ideal_foreground_per_cm2_sr_day_keV_per_um", "ideal_mc_standard_error",
            "foreground_after_cuts_per_cm2_sr_day_keV_per_um",
            "cuts_mc_standard_error", "observed_per_cm2_sr_day_keV_per_um",
        ])
        for row in merged_bins:
            writer.writerow([row[field] for field in (
                "lower_keV_per_um", "upper_keV_per_um", "center_keV_per_um",
                "ideal_foreground_per_cm2_sr_day_keV_per_um", "ideal_mc_standard_error",
                "foreground_after_cuts_per_cm2_sr_day_keV_per_um",
                "cuts_mc_standard_error", "observed_per_cm2_sr_day_keV_per_um",
            )])
    print(args.output)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

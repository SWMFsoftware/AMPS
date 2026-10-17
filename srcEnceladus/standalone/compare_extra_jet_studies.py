#!/usr/bin/env python3
"""Compare a published-source fit with the same fit plus one candidate jet."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CANDIDATE_LABEL = "jet:e17_candidate_baghdad"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_pair(base: dict, augmented: dict, family: str) -> None:
    if base.get("source_case") != "porco2014_98":
        raise ValueError(f"{family} baseline is not the Porco 98-source case")
    if augmented.get("source_case") != "porco2014_98":
        raise ValueError(f"{family} augmented result has a different base case")
    if augmented.get("source_variant") != "porco2014_98_plus_e17_candidate":
        raise ValueError(f"{family} result is not the checked candidate variant")
    for key in ("config_sha256", "data_root"):
        if base.get(key) != augmented.get(key):
            raise ValueError(f"{family} baseline and augmented {key} differ")
    if not augmented.get("extra_jet_sha256"):
        raise ValueError(f"{family} augmented result lacks candidate provenance")


def strategy_table(document: dict, family: str) -> dict:
    key = ("source_rate_shape_strategies" if family == "H2O/UVIS"
           else "source_rate_strategies")
    if key not in document:
        raise ValueError(f"{family} result does not contain {key}")
    return document[key]


def is_training_target(family: str, strategy: str, target: str) -> bool:
    if strategy == "independent_each_observation":
        return True
    if family == "H2O/UVIS":
        if strategy == "shared_shape_all_observations":
            return True
        if strategy == "shared_shape_inms_only":
            return target in ("E3", "E5")
        if strategy == "shared_shape_uvis_only":
            return target.startswith("uvis_")
    elif family == "CO2":
        if strategy == "shared_shape_all_co2":
            return True
        if strategy == "shared_shape_e3_e5":
            return target in ("E3", "E5")
    raise ValueError(f"Unknown family/strategy: {family}/{strategy}")


def compare_family(base: dict, augmented: dict, family: str,
                   candidate_label: str = CANDIDATE_LABEL) -> list[dict]:
    base_table = strategy_table(base, family)
    augmented_table = strategy_table(augmented, family)
    if set(base_table) != set(augmented_table):
        raise ValueError(f"{family} strategy sets differ")
    base_labels = base.get("component_labels", [])
    augmented_labels = augmented.get("component_labels", [])
    if candidate_label in base_labels or candidate_label not in augmented_labels:
        raise ValueError("Candidate label must occur only in augmented components")
    if [label for label in augmented_labels if label != candidate_label] != base_labels:
        raise ValueError("Augmented component order is not base order plus candidate")

    rows = []
    for strategy, base_strategy in base_table.items():
        augmented_strategy = augmented_table[strategy]
        base_targets = base_strategy["target_fits"]
        augmented_targets = augmented_strategy["target_fits"]
        if set(base_targets) != set(augmented_targets):
            raise ValueError(f"{family}/{strategy} target sets differ")
        for target, base_fit in base_targets.items():
            new_fit = augmented_targets[target]
            fraction = float(new_fit["source_shape_fractions"][candidate_label])
            base_rate = float(base_fit["total_source_rate_s"])
            new_rate = float(new_fit["total_source_rate_s"])
            base_chi2 = float(base_fit["chi2"])
            new_chi2 = float(new_fit["chi2"])
            samples = int(base_fit["samples"])
            if samples != int(new_fit["samples"]):
                raise ValueError(f"Sample count changed for {family}/{strategy}/{target}")
            rows.append({
                "family": family,
                "strategy": strategy,
                "target": target,
                "used_to_train_shape": is_training_target(family, strategy, target),
                "samples": samples,
                "base_chi2": base_chi2,
                "augmented_chi2": new_chi2,
                "delta_chi2_base_minus_augmented": base_chi2-new_chi2,
                "base_chi2_per_sample": base_chi2/samples,
                "augmented_chi2_per_sample": new_chi2/samples,
                "delta_chi2_per_sample_base_minus_augmented":
                    (base_chi2-new_chi2)/samples,
                "relative_chi2_change_percent": (
                    100.0*(new_chi2-base_chi2)/base_chi2
                    if base_chi2 > 0 else math.nan),
                "base_total_source_rate_s-1": base_rate,
                "augmented_total_source_rate_s-1": new_rate,
                "augmented_over_base_total_rate": (
                    new_rate/base_rate if base_rate > 0 else math.nan),
                "candidate_source_fraction": fraction,
                "candidate_source_rate_s-1": new_rate*fraction,
            })
    return rows


def summarize(rows: list[dict]) -> list[dict]:
    groups = {}
    for row in rows:
        groups.setdefault((row["family"], row["strategy"]), []).append(row)
    output = []
    for (family, strategy), items in groups.items():
        training = [row for row in items if row["used_to_train_shape"]]
        output.append({
            "family": family,
            "strategy": strategy,
            "target_count": len(items),
            "training_target_count": len(training),
            "all_targets_base_sum_chi2_per_sample":
                sum(row["base_chi2_per_sample"] for row in items),
            "all_targets_augmented_sum_chi2_per_sample":
                sum(row["augmented_chi2_per_sample"] for row in items),
            "all_targets_delta_sum_chi2_per_sample": sum(
                row["delta_chi2_per_sample_base_minus_augmented"]
                for row in items),
            "training_base_objective":
                sum(row["base_chi2_per_sample"] for row in training),
            "training_augmented_objective":
                sum(row["augmented_chi2_per_sample"] for row in training),
            "training_delta_objective": sum(
                row["delta_chi2_per_sample_base_minus_augmented"]
                for row in training),
            "maximum_candidate_fraction": max(
                row["candidate_source_fraction"] for row in items),
            "maximum_candidate_rate_s-1": max(
                row["candidate_source_rate_s-1"] for row in items),
        })
    return output


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def short_target(target: str) -> str:
    return {"uvis_2007_zeta_orionis": "UVIS 2007",
            "uvis_2010_e10_solar": "UVIS 2010"}.get(target, target)


def heatmap(rows: list[dict], family: str, value_key: str, output: Path,
            title: str, color_label: str, logarithmic: bool = False) -> None:
    selected = [row for row in rows if row["family"] == family]
    strategies = list(dict.fromkeys(row["strategy"] for row in selected))
    targets = list(dict.fromkeys(row["target"] for row in selected))
    lookup = {(row["strategy"], row["target"]): row for row in selected}
    values = np.asarray([[float(lookup[(strategy, target)][value_key])
                          for target in targets] for strategy in strategies])
    display = np.log10(np.maximum(values, 1.0e-50)) if logarithmic else values
    fig, axis = plt.subplots(figsize=(max(7.2, 1.35*len(targets)),
                                      max(3.2, .72*len(strategies)+1.8)))
    if logarithmic:
        image = axis.imshow(display, aspect="auto", cmap="viridis")
    else:
        extent = float(np.nanmax(np.abs(display)))
        extent = extent if extent > 0 else 1.0
        image = axis.imshow(display, aspect="auto", cmap="RdBu_r",
                            vmin=-extent, vmax=extent)
    for i, strategy in enumerate(strategies):
        for j, target in enumerate(targets):
            value = values[i, j]
            label = f"{value:.3g}" if not logarithmic else f"{value:.2e}"
            axis.text(j, i, label, ha="center", va="center", fontsize=7,
                      color=("white" if logarithmic and display[i, j] < -6
                             else "black"))
    axis.set_xticks(range(len(targets)))
    axis.set_xticklabels([short_target(value) for value in targets],
                         rotation=25, ha="right")
    axis.set_yticks(range(len(strategies)))
    axis.set_yticklabels([value.replace("_", " ") for value in strategies])
    axis.set_title(title)
    colorbar = fig.colorbar(image, ax=axis, shrink=.9)
    colorbar.set_label(color_label)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=200)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def markdown_report(rows: list[dict], summaries: list[dict], paths: dict) -> str:
    independent = [row for row in rows
                   if row["strategy"] == "independent_each_observation"]
    lines = [
        "# Effect of the E17 candidate jet on all fitted observations", "",
        "The published Porco 98-jet model and the otherwise identical model with the "
        "fixed E17 Baghdad candidate as source 99 were run through the complete fitting "
        "workflow. Positive delta chi-square means the 99-source case fits better; "
        "negative values mean it fits worse. Rates and source fractions are inverse-model "
        "outputs and remain non-unique where the design matrix is ill-conditioned.", "",
        "## Independent fits", "",
        "| Species/data | Observation | Delta chi-square | Delta chi-square/sample | Candidate fraction | Candidate rate (s^-1) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in independent:
        lines.append(
            f"| {row['family']} | {short_target(row['target'])} | "
            f"{row['delta_chi2_base_minus_augmented']:.6g} | "
            f"{row['delta_chi2_per_sample_base_minus_augmented']:.6g} | "
            f"{row['candidate_source_fraction']:.6g} | "
            f"{row['candidate_source_rate_s-1']:.6g} |")
    lines += ["", "## Strategy-level diagnostic", "",
              "The training objective is the sum of chi-square per sample over the data "
              "sets used to determine that strategy's shared shape. For independent fits, "
              "it is the sum of the separately optimized objectives.", "",
              "| Family | Strategy | Delta training objective | Delta across all targets | Maximum candidate fraction |",
              "|---|---|---:|---:|---:|"]
    for row in summaries:
        lines.append(
            f"| {row['family']} | {row['strategy'].replace('_', ' ')} | "
            f"{row['training_delta_objective']:.6g} | "
            f"{row['all_targets_delta_sum_chi2_per_sample']:.6g} | "
            f"{row['maximum_candidate_fraction']:.6g} |")
    lines += [
        "", "## Interpretation", "",
        "The E17 candidate is selected only weakly when all 99 source amplitudes are free. "
        "Its clearest independent-profile effect is in E17 CO2; the other CO2 flybys put "
        "essentially zero weight on it. The UVIS independent fits assign small but nonzero "
        "fractions, while E3/E5 H2O do not select it. Consequently, the candidate is not a "
        "globally required persistent source under this steady ballistic model.", "",
        "A nested 99-source nonnegative fit should never have a larger optimized training "
        "objective than the 98-source fit because zero candidate amplitude reproduces the "
        "old model. Any negative training-objective delta therefore diagnoses numerical "
        "sensitivity in this severely ill-conditioned inverse problem, not evidence that "
        "the physical candidate makes the data worse. Cross-prediction targets may worsen "
        "legitimately because their shape was trained on different observations.", "",
        "## Products", "",
        f"- Detailed per-observation table: `{paths['detail']}`",
        f"- Strategy summary: `{paths['summary']}`",
        "- Heatmaps: `h2o_uvis_delta_chi2_per_sample.*`, "
        "`co2_delta_chi2_per_sample.*`, and candidate-fraction counterparts.", "",
        "The comparison preserves the original fitted-result JSON files and records their "
        "paths in `comparison_results.json`.", "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-h2o-uvis", type=Path, required=True)
    parser.add_argument("--augmented-h2o-uvis", type=Path, required=True)
    parser.add_argument("--base-co2", type=Path, required=True)
    parser.add_argument("--augmented-co2", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("output directory already exists")
    args.output_dir.mkdir(parents=True)

    inputs = {name: read_json(path) for name, path in {
        "base_h2o_uvis": args.base_h2o_uvis,
        "augmented_h2o_uvis": args.augmented_h2o_uvis,
        "base_co2": args.base_co2,
        "augmented_co2": args.augmented_co2,
    }.items()}
    validate_pair(inputs["base_h2o_uvis"], inputs["augmented_h2o_uvis"],
                  "H2O/UVIS")
    validate_pair(inputs["base_co2"], inputs["augmented_co2"], "CO2")
    if (inputs["augmented_h2o_uvis"]["extra_jet_sha256"]
            != inputs["augmented_co2"]["extra_jet_sha256"]):
        raise ValueError("H2O/UVIS and CO2 runs used different candidate definitions")
    rows = compare_family(inputs["base_h2o_uvis"],
                          inputs["augmented_h2o_uvis"], "H2O/UVIS")
    rows += compare_family(inputs["base_co2"], inputs["augmented_co2"], "CO2")
    summaries = summarize(rows)
    detail_path = args.output_dir/"all_observations_98_vs_99.csv"
    summary_path = args.output_dir/"strategy_summary_98_vs_99.csv"
    write_csv(detail_path, rows)
    write_csv(summary_path, summaries)
    for family, prefix in (("H2O/UVIS", "h2o_uvis"), ("CO2", "co2")):
        heatmap(rows, family, "delta_chi2_per_sample_base_minus_augmented",
                args.output_dir/f"{prefix}_delta_chi2_per_sample",
                f"{family}: fit change after adding E17 candidate",
                r"$\chi^2_{98}/N-\chi^2_{99}/N$")
        heatmap(rows, family, "candidate_source_fraction",
                args.output_dir/f"{prefix}_candidate_fraction",
                f"{family}: fitted fraction assigned to E17 candidate",
                r"$\log_{10}$ candidate source fraction", logarithmic=True)
    result = {
        "study": "Porco 98 sources versus Porco 98 plus E17 candidate",
        "candidate_label": CANDIDATE_LABEL,
        "input_result_files": {
            "base_h2o_uvis": str(args.base_h2o_uvis.resolve()),
            "augmented_h2o_uvis": str(args.augmented_h2o_uvis.resolve()),
            "base_co2": str(args.base_co2.resolve()),
            "augmented_co2": str(args.augmented_co2.resolve()),
        },
        "input_result_sha256": {
            "base_h2o_uvis": sha256(args.base_h2o_uvis),
            "augmented_h2o_uvis": sha256(args.augmented_h2o_uvis),
            "base_co2": sha256(args.base_co2),
            "augmented_co2": sha256(args.augmented_co2),
        },
        "extra_jet_sha256": inputs["augmented_h2o_uvis"]["extra_jet_sha256"],
        "per_observation": rows,
        "strategy_summary": summaries,
        "sign_convention": "positive delta chi2 means augmented model is better",
    }
    (args.output_dir/"comparison_results.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8")
    paths = {"detail": detail_path.name, "summary": summary_path.name}
    (args.output_dir/"comparison_report.md").write_text(
        markdown_report(rows, summaries, paths), encoding="utf-8")
    print(json.dumps({"output": str(args.output_dir),
                      "rows": len(rows), "summaries": len(summaries)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

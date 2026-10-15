#!/usr/bin/env python3
"""Create a checked cross-case synthesis from three completed plume studies."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CASES = (
    "tenishev2014_hybrid_8",
    "porco2014_98",
    "portyankina2022_solar_41",
)
SHORT = {
    "tenishev2014_hybrid_8": "Hybrid 8+stripes",
    "porco2014_98": "Porco 98",
    "portyankina2022_solar_41": "Portyankina 41",
}
DATASETS = (
    "E3",
    "E5",
    "uvis_2007_zeta_orionis",
    "uvis_2010_e10_solar",
)
DISPLAY = {
    "E3": "E3 INMS",
    "E5": "E5 INMS",
    "uvis_2007_zeta_orionis": "2007 UVIS",
    "uvis_2010_e10_solar": "2010 UVIS",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_counts(result):
    labels = result["component_labels"]
    return (
        sum(label.startswith("jet:") for label in labels),
        sum(label.startswith("stripe:") for label in labels),
    )


def independent_rows(case, result):
    jets, stripes = source_counts(result)
    rows = []
    for flyby, dataset, total_key in (
            ("e3", "E3", "E3_independent"),
            ("e5", "E5", "E5_independent")):
        fit = result["inms_fits"][flyby]
        samples = len(fit["predicted"])
        dof = samples-fit["rank"]
        rows.append({
            "source_case": case,
            "dataset": dataset,
            "jets": jets,
            "stripe_components": stripes,
            "samples": samples,
            "fitted_columns": len(fit["coefficients"]),
            "rank": fit["rank"],
            "effective_rank_1e-3": fit[
                "effective_rank_by_relative_singular_cutoff"]["1e-3"],
            "residual_dof": dof,
            "chi2": fit["chi2"],
            "chi2_per_sample": fit["chi2"]/samples,
            "reduced_chi2": fit["chi2"]/dof if dof > 0 else None,
            "condition_number_scaled": fit["condition_number_scaled"],
            "total_source_rate_s": result["total_source_rate_s"][total_key],
            "quadrature_max_relative_fitted_column_difference": None,
        })
    for event in DATASETS[2:]:
        fit = result["uvis_fits"][event]
        samples = len(fit["tau_prediction"])
        rows.append({
            "source_case": case,
            "dataset": event,
            "jets": jets,
            "stripe_components": stripes,
            "samples": samples,
            "fitted_columns": len(fit["labels"])+1,
            "rank": fit["rank"],
            "effective_rank_1e-3": fit["effective_rank"]["1e-3"],
            "residual_dof": fit["shape_dof"],
            "chi2": fit["chi2_shape"],
            "chi2_per_sample": fit["chi2_shape"]/samples,
            "reduced_chi2": fit["shape_reduced_chi2"],
            "condition_number_scaled": fit["condition_number_scaled"],
            "total_source_rate_s": result["total_source_rate_s"][
                event+"_independent"],
            "quadrature_max_relative_fitted_column_difference": result[
                "quadrature_convergence"][event][
                    "max_relative_difference_fitted_total_column"],
        })
    return rows


def format_reduced(row):
    if row["reduced_chi2"] is None:
        return "undefined (0 dof)"
    return f"{row['reduced_chi2']:.3f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.study_root.resolve()
    outputs = {
        "csv": root/"source_case_comparison.csv",
        "json": root/"source_case_comparison.json",
        "md": root/"source_case_comparison.md",
        "png": root/"source_case_comparison.png",
        "pdf": root/"source_case_comparison.pdf",
    }
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing:
        parser.error("refusing to overwrite: "+", ".join(existing))

    results = {}
    input_hashes = {}
    rows = []
    for case in CASES:
        path = root/case/"study_results.json"
        if not path.is_file():
            parser.error("missing completed study: "+str(path))
        result = json.loads(path.read_text(encoding="utf-8"))
        if result.get("source_case") != case:
            parser.error("source-case mismatch in "+str(path))
        results[case] = result
        input_hashes[str(path)] = sha256(path)
        rows.extend(independent_rows(case, result))

    fields = list(rows[0])
    with outputs["csv"].open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    synthesis = {
        "study": "comparison of three Enceladus neutral-plume source cases",
        "study_root": str(root),
        "source_cases": list(CASES),
        "input_study_sha256": input_hashes,
        "software_sha256": results[CASES[0]]["software_sha256"],
        "quadrature_order": results[CASES[0]]["quadrature_order"],
        "independent_fit_summary": rows,
        "combined_fit_total_source_rate_s": {
            case: {
                model: results[case]["total_source_rate_s"][model]
                for model in ("INMS_combined", "UVIS_combined",
                              "steady_balanced_diagnostic")}
            for case in CASES},
        "combined_fit_chi2_per_sample": {
            case: {
                model: {
                    dataset: results[case]["cross_prediction"][model][dataset][
                        "chi2_per_sample"]
                    for dataset in DATASETS}
                for model in ("INMS_combined", "UVIS_combined",
                              "steady_balanced_diagnostic")}
            for case in CASES},
        "combined_objective_chi2_per_sample": {
            case: {
                "INMS_combined": {
                    dataset: item["chi2_per_sample"]
                    for dataset, item in results[case][
                        "inms_combined_fit"]["per_dataset"].items()},
                "UVIS_combined": {
                    dataset: item["chi2_per_sample"]
                    for dataset, item in results[case][
                        "uvis_combined_fit"]["per_dataset"].items()},
                "steady_balanced_diagnostic": {
                    dataset: item["chi2_per_sample"]
                    for dataset, item in results[case][
                        "steady_balanced_diagnostic"]["per_dataset"].items()},
            }
            for case in CASES},
        "shared_source_chi2_per_sample": {
            case: {
                dataset: results[case]["cross_prediction"][
                    "steady_balanced_diagnostic"][dataset]["chi2_per_sample"]
                for dataset in DATASETS}
            for case in CASES},
        "interpretation": [
            "No formal model-selection probability is assigned because the data operators, active nonnegative bounds, effective ranks, and systematic errors do not define comparable regular likelihoods.",
            "The 98-jet 2007 UVIS fit is saturated at rank 97 for 97 samples and has zero residual degrees of freedom.",
            "Individual source amplitudes are not uniquely identifiable in any case; total production remains conditional on source geometry and gas kinematics.",
        ],
    }
    with outputs["json"].open("x", encoding="utf-8") as stream:
        json.dump(synthesis, stream, indent=2, allow_nan=False)
        stream.write("\n")

    by_key = {(row["source_case"], row["dataset"]): row for row in rows}
    report = [
        "# Comparison of the three Enceladus source cases", "",
        "All cases use the same qualified E3/E5 density profiles, 2007 stellar "
        "and 2010 solar UVIS products, collisionless physics, background model, "
        "quadrature order 192, nonnegative source-rate constraints, and one "
        "independent fit per observation. Geometry and source kinematics differ "
        "as documented in the standalone README.", "",
        "## Independent fit quality", "",
        "The first number is chi-square per sample, which remains defined for a "
        "saturated model. Parentheses give reduced chi-square using formal matrix "
        "rank; it is undefined when no residual degrees of freedom remain.", "",
        "| Source case | E3 INMS | E5 INMS | 2007 UVIS shape | 2010 UVIS shape |",
        "|---|---:|---:|---:|---:|",
    ]
    for case in CASES:
        cells = []
        for dataset in DATASETS:
            row = by_key[case, dataset]
            cells.append(
                f"{row['chi2_per_sample']:.3f} ({format_reduced(row)})")
        report.append(f"| {SHORT[case]} | "+" | ".join(cells)+" |")

    report += ["", "## Independently inferred total source rate", "",
               "| Source case | E3 (s^-1) | E5 (s^-1) | 2007 UVIS (s^-1) | 2010 UVIS (s^-1) |",
               "|---|---:|---:|---:|---:|"]
    for case in CASES:
        rates = [by_key[case, dataset]["total_source_rate_s"]
                 for dataset in DATASETS]
        report.append(f"| {SHORT[case]} | "+" | ".join(
            f"{value:.6e}" for value in rates)+" |")

    report += ["", "## Identifiability", "",
               "| Source case | E3 rank / columns (effective rank at 1e-3) | E3 condition | E5 rank / columns (effective rank at 1e-3) | E5 condition |",
               "|---|---:|---:|---:|---:|"]
    for case in CASES:
        e3, e5 = by_key[case, "E3"], by_key[case, "E5"]
        report.append(
            f"| {SHORT[case]} | {e3['rank']} / {e3['fitted_columns']} "
            f"({e3['effective_rank_1e-3']}) | {e3['condition_number_scaled']:.3e} | "
            f"{e5['rank']} / {e5['fitted_columns']} "
            f"({e5['effective_rank_1e-3']}) | {e5['condition_number_scaled']:.3e} |")

    report += ["", "## Combined fitting by instrument", "",
               "One vector is fitted jointly to E3+E5 and a different vector is "
               "fitted jointly to the two UVIS occultations. Values are objective "
               "chi-square per sample in local-density space for INMS and in the "
               "documented physical-column-proxy space for UVIS.", "",
               "| Source case | INMS-combined rate (s^-1) | E3 | E5 | UVIS-combined rate (s^-1) | 2007 UVIS | 2010 UVIS |",
               "|---|---:|---:|---:|---:|---:|---:|"]
    for case in CASES:
        rate = synthesis["combined_fit_total_source_rate_s"][case]
        metric = synthesis["combined_objective_chi2_per_sample"][case]
        report.append(
            f"| {SHORT[case]} | {rate['INMS_combined']:.6e} | "
            f"{metric['INMS_combined']['E3']:.3f} | "
            f"{metric['INMS_combined']['E5']:.3f} | "
            f"{rate['UVIS_combined']:.6e} | "
            f"{metric['UVIS_combined']['uvis_2007_zeta_orionis']:.3f} | "
            f"{metric['UVIS_combined']['uvis_2010_e10_solar']:.3f} |")

    report += ["", "## Shared-source cross-epoch diagnostic", "",
               "A smaller chi-square per sample means that one balanced frozen "
               "source vector transfers less poorly across epochs. This is not a "
               "formal joint likelihood.", "",
               "| Source case | E3 | E5 | 2007 UVIS | 2010 UVIS |",
               "|---|---:|---:|---:|---:|"]
    for case in CASES:
        metric = synthesis["shared_source_chi2_per_sample"][case]
        report.append(f"| {SHORT[case]} | "+" | ".join(
            f"{metric[dataset]:.3f}" for dataset in DATASETS)+" |")

    report += [
        "", "## Interpretation", "",
        "- The hybrid 8+stripes case gives the best rank-adjusted E3 and E5 fits "
        "and the best nonsaturated 2007 UVIS fit. Its E3/E5 source totals remain "
        "in the physically familiar few-times-10^28 s^-1 range.",
        "- The 98-jet case gives the lowest raw E5 and UVIS shape residuals and "
        "the least-poor balanced cross-epoch diagnostic. However, the 2007 UVIS "
        "fit is saturated, its INMS matrices have condition numbers above 10^16, "
        "and its E3/E5 totals inflate to 10^30--10^31 s^-1. Those totals are a "
        "degeneracy/kinematics warning, not credible plume production estimates.",
        "- The 41-jet case avoids the extreme INMS source-rate inflation and has "
        "excellent fitted-column quadrature convergence, but it does not improve "
        "the INMS or 2007 UVIS fit over the hybrid case. It improves the 2010 UVIS "
        "shape relative to the hybrid case, consistent with its 2010-selection "
        "provenance.",
        "- No source case supplies a statistically adequate time-invariant map "
        "for every epoch. The data support continuing with source-family ensembles "
        "and explicit temporal/regularization physics rather than selecting 98 "
        "independent jet rates from these observations.", "",
        "## Numerical checks", "",
    ]
    for case in CASES:
        values = [by_key[case, event][
            "quadrature_max_relative_fitted_column_difference"]
            for event in DATASETS[2:]]
        report.append(
            f"- {SHORT[case]}: maximum order-192 versus order-224 fitted-column "
            f"difference is {values[0]:.3e} (2007) and {values[1]:.3e} (2010).")
    report += [
        "", "The CSV and JSON beside this report retain exact values, input-study "
        "hashes, software hashes, and machine-readable interpretation limits.", ""]
    outputs["md"].write_text("\n".join(report), encoding="utf-8")

    colors = ("#1b9e77", "#d95f02", "#7570b3")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    index = np.arange(len(DATASETS))
    width = 0.25
    for offset, (case, color) in enumerate(zip(CASES, colors)):
        axes[0, 0].bar(index+(offset-1)*width,
                       [by_key[case, d]["chi2_per_sample"] for d in DATASETS],
                       width, label=SHORT[case], color=color)
        axes[0, 1].plot(index,
                        [by_key[case, d]["total_source_rate_s"] for d in DATASETS],
                        "o-", label=SHORT[case], color=color)
        axes[1, 0].plot(index,
                        [by_key[case, d]["effective_rank_1e-3"]
                         / by_key[case, d]["fitted_columns"] for d in DATASETS],
                        "o-", label=SHORT[case], color=color)
        combined = synthesis["combined_objective_chi2_per_sample"][case]
        matched = {
            "E3": combined["INMS_combined"]["E3"],
            "E5": combined["INMS_combined"]["E5"],
            "uvis_2007_zeta_orionis": combined["UVIS_combined"][
                "uvis_2007_zeta_orionis"],
            "uvis_2010_e10_solar": combined["UVIS_combined"][
                "uvis_2010_e10_solar"],
        }
        axes[1, 1].bar(index+(offset-1)*width,
                       [matched[d] for d in DATASETS], width,
                       label=SHORT[case], color=color)
    for ax in axes.flat:
        ax.set_xticks(index)
        ax.set_xticklabels([DISPLAY[d] for d in DATASETS], rotation=18,
                           ha="right")
        ax.grid(True, axis="y", alpha=0.25)
    axes[0, 0].set_ylabel("Independent-fit chi-square / sample")
    axes[0, 0].set_yscale("log")
    axes[0, 1].set_ylabel("Total fitted source rate (s$^{-1}$)")
    axes[0, 1].set_yscale("log")
    axes[1, 0].set_ylabel("Effective rank / fitted columns (1e-3)")
    axes[1, 0].set_ylim(0, 1.05)
    axes[1, 1].set_ylabel("Within-instrument combined chi-square / sample")
    axes[1, 1].set_yscale("log")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Enceladus source-case fitting comparison")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outputs["png"], dpi=180)
    fig.savefig(outputs["pdf"])
    plt.close(fig)
    print(json.dumps({name: str(path) for name, path in outputs.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

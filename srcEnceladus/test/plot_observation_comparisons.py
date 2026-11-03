#!/usr/bin/env python3
"""Render fail-closed AMPS comparisons for requested INMS and UVIS data."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

import plot_flyby_density


def weighted_line(x: list[float], y: list[float], sigma: list[float]
                  ) -> tuple[float, float, float]:
    """Fit y=alpha*x+baseline with alpha>=0 and a signed baseline."""
    weight = [1.0/(value*value) for value in sigma]
    total = sum(weight)
    mean_x = sum(w*a for w, a in zip(weight, x))/total
    mean_y = sum(w*a for w, a in zip(weight, y))/total
    denominator = sum(w*(a-mean_x)**2 for w, a in zip(weight, x))
    alpha = (sum(w*(a-mean_x)*(b-mean_y)
                 for w, a, b in zip(weight, x, y))/denominator
             if denominator > 0.0 else 0.0)
    alpha = max(0.0, alpha)
    baseline = sum(w*(b-alpha*a) for w, a, b in zip(weight, x, y))/total
    chi2 = sum(((alpha*a+baseline-b)/s)**2
               for a, b, s in zip(x, y, sigma))
    return alpha, baseline, chi2


def validate_zones(path: Path, zones: list[plot_flyby_density.Zone],
                   expected: int | None) -> None:
    if expected is None:
        return
    outputs = [zone.output for zone in zones]
    row_counts = [len(zone.rows) for zone in zones]
    if (len(zones) != expected or outputs != list(range(expected)) or
            len(set(row_counts)) != 1):
        raise ValueError(
            f"{path} zone validation failed: outputs={outputs}, rows={row_counts}, "
            f"expected_zones={expected}")


def plot_inms(path: Path, output: Path) -> dict[str, object]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    variables, zones = plot_flyby_density.read_tecplot(path)
    zone = zones[-1]
    time_i = variables.index("time_s")
    domain_i = variables.index("in_domain")
    model_columns = [name for name in variables
                     if name.startswith("model_n_") and name.endswith("_m-3")]
    observed_columns = [name for name in variables
                        if name.startswith("observed_n_")]
    context_present = "observed_context_value" in variables
    text = path.read_text(encoding="utf-8")
    context_match = re.search(r"^# context_observable=([^ ]+)", text,
                              flags=re.MULTILINE)
    if context_present != (context_match is not None):
        raise ValueError(f"{path} has inconsistent native-context schema/receipt")
    if not model_columns:
        raise ValueError(f"{path} has no modeled INMS density columns")
    rows = [row for row in zone.rows if row[domain_i] == 1.0]
    if not rows:
        raise ValueError(f"{path} has no in-domain observation points")
    panel_count = len(model_columns) + int(context_present)
    figure, axes = plt.subplots(panel_count, 1,
                                figsize=(10.5, 3.7*panel_count),
                                sharex=True, constrained_layout=True)
    if panel_count == 1:
        axes = [axes]
    metrics: dict[str, object] = {}
    for axis, model_name in zip(axes, model_columns):
        species = model_name[len("model_n_"):-len("_m-3")]
        model_i = variables.index(model_name)
        times = [row[time_i] for row in rows]
        model = [row[model_i] for row in rows]
        axis.plot(times, [value if value > 0.0 else math.nan for value in model],
                  color="#2166ac", lw=2.0,
                  label=f"AMPS {species}, final sampling window")
        observed_name = f"observed_n_{species}_m-3"
        if observed_name in observed_columns:
            observed_i = variables.index(observed_name)
            sigma_i = variables.index(f"observed_sigma_{species}_m-3")
            observed = [row[observed_i] for row in rows]
            sigma = [row[sigma_i] for row in rows]
            positive = [(t, y, s) for t, y, s in zip(times, observed, sigma)
                        if y > 0.0]
            if positive:
                axis.errorbar([item[0] for item in positive],
                              [item[1] for item in positive],
                              yerr=[item[2] for item in positive], fmt="o",
                              color="black", ms=3.2, elinewidth=.6,
                              label=f"qualified INMS {species}")
            zero = [(t, s) for t, y, s in zip(times, observed, sigma) if y == 0.0]
            if zero:
                axis.scatter([item[0] for item in zero],
                             [item[1] for item in zero], marker="v", s=22,
                             facecolors="none", edgecolors="black",
                             label="zero count at one-count sensitivity")
            chi2 = sum(((a-b)/s)**2
                       for a, b, s in zip(model, observed, sigma))
            ratios = [a/b for a, b in zip(model, observed)
                      if a > 0.0 and b > 0.0]
            ratios.sort()
            metrics[species] = {
                "comparison": "response-qualified local number density",
                "samples_in_domain": len(rows),
                "chi2": chi2,
                "chi2_per_sample": chi2/len(rows),
                "median_model_over_observed_positive": (
                    ratios[len(ratios)//2] if ratios else None),
            }
        else:
            metrics[species] = {
                "comparison": "model_only_no_response_qualified_density",
                "samples_in_domain": len(rows),
            }
        if any(value > 0.0 for value in model):
            axis.set_yscale("log")
        axis.set_ylabel(f"n({species}) [m$^{{-3}}$]")
        axis.grid(True, which="both", alpha=.25)
        axis.legend(fontsize=8)
    if context_present:
        axis = axes[-1]
        context_i = variables.index("observed_context_value")
        context_sigma_i = variables.index("observed_context_sigma")
        all_times = [row[time_i] for row in zone.rows]
        context = [row[context_i] for row in zone.rows]
        context_sigma = [row[context_sigma_i] for row in zone.rows]
        axis.errorbar(all_times, context, yerr=context_sigma, fmt=".-",
                      color="#444444", ms=2.2, lw=.55, elinewidth=.25,
                      alpha=.75, label="archived native INMS observable")
        covered = [(row[time_i], row[context_i]) for row in zone.rows
                   if row[domain_i] == 1.0]
        if covered:
            axis.scatter([item[0] for item in covered],
                         [item[1] for item in covered], s=8,
                         color="#b2182b", zorder=3,
                         label="samples inside current AMPS domain")
        context_name = context_match.group(1)
        axis.set_ylabel(context_name.replace("_", " "), fontsize=8)
        axis.grid(True, alpha=.25)
        axis.legend(fontsize=8)
        metrics["native_context"] = {
            "observable": context_name,
            "comparison": "display_only_no_density_or_count_response",
            "samples_total": len(zone.rows),
            "samples_in_domain": len(covered),
        }
    axes[-1].set_xlabel("Time from closest approach [s]")
    figure.suptitle(
        f"{path.stem}: current E3 source distribution transferred without refit\n"
        "Density comparisons and native-count context retain distinct units")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    return {"instrument": "INMS", "file": str(path), "png": str(output),
            "final_output": zone.output, "zones": len(zones),
            "rows_per_zone": len(zone.rows), "metrics": metrics}


def plot_uvis(path: Path, output: Path) -> dict[str, object]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    variables, zones = plot_flyby_density.read_tecplot(path)
    zone = zones[-1]
    required = {"time_from_min_altitude_s", "model_column_H2O_m-2",
                "observed_optical_depth", "observed_sigma_optical_depth"}
    if not required.issubset(variables):
        raise ValueError(f"{path} is missing UVIS columns: {sorted(required-set(variables))}")
    time_i = variables.index("time_from_min_altitude_s")
    column_i = variables.index("model_column_H2O_m-2")
    tau_i = variables.index("observed_optical_depth")
    sigma_i = variables.index("observed_sigma_optical_depth")
    times = [row[time_i] for row in zone.rows]
    column = [row[column_i] for row in zone.rows]
    tau = [row[tau_i] for row in zone.rows]
    sigma = [row[sigma_i] for row in zone.rows]
    if not any(value > 0.0 for value in column):
        raise ValueError(f"{path} has no positive modeled H2O column")
    alpha, baseline, chi2 = weighted_line(column, tau, sigma)
    prediction = [alpha*value+baseline for value in column]
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r"published_h2o_column_m-2=([0-9.eE+-]+) "
        r"published_h2o_column_sigma_m-2=([0-9.eE+-]+)", text)
    if match is None:
        raise ValueError(f"{path} lacks its published-column receipt")
    published = float(match.group(1))
    published_sigma = float(match.group(2))

    figure, axes = plt.subplots(2, 1, figsize=(11.0, 7.3), sharex=True,
                                constrained_layout=True)
    axes[0].errorbar(times, tau, yerr=sigma, fmt=".", color="black", ms=3,
                     elinewidth=.35, label="UVIS optical-depth proxy")
    axes[0].plot(times, prediction, color="#762a83", lw=1.8,
                 label="weighted AMPS-column shape fit")
    axes[0].set_ylabel("Optical-depth proxy")
    axes[0].legend(fontsize=8)
    axes[1].plot(times, column, color="#1b7837", lw=1.8,
                 label="AMPS H$_2$O column in finite domain")
    axes[1].axhspan(max(0.0, published-published_sigma),
                    published+published_sigma, color="#d95f02", alpha=.18,
                    label="published representative/peak column ±1σ")
    axes[1].axhline(published, color="#d95f02", lw=1.0)
    axes[1].set_ylabel("H$_2$O column [m$^{-2}$]")
    axes[1].set_xlabel("Time from minimum tangent altitude [s]")
    axes[1].set_yscale("log")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(True, which="both", alpha=.25)
    figure.suptitle(
        f"{path.stem}: UVIS comparison for the current E3 source distribution\n"
        "Optical depth constrains shape; the published column constrains absolute scale")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    peak = max(column)
    return {
        "instrument": "UVIS", "file": str(path), "png": str(output),
        "final_output": zone.output, "zones": len(zones),
        "samples": len(zone.rows),
        "shape_fit": {"column_to_optical_depth": alpha,
                      "signed_baseline": baseline, "chi2": chi2,
                      "chi2_per_sample": chi2/len(zone.rows)},
        "absolute_column": {"model_peak_m-2": peak,
                            "published_m-2": published,
                            "published_sigma_m-2": published_sigma,
                            "model_peak_over_published": peak/published},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--expected-zones", type=int)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    results = []
    for item in args.inputs:
        path = item.resolve()
        output = output_dir/f"{path.stem}.comparison.png"
        variables, zones = plot_flyby_density.read_tecplot(path)
        validate_zones(path, zones, args.expected_zones)
        if "observed_optical_depth" in variables:
            results.append(plot_uvis(path, output))
        else:
            results.append(plot_inms(path, output))
        print(f"PASS {path.name} -> {output}")
    summary = (args.summary.resolve() if args.summary else
               output_dir/"observation_comparison_summary.json")
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(json.dumps({"status": "PASS",
                                  "expected_zones": args.expected_zones,
                                  "comparisons": results}, indent=2,
                                  sort_keys=True)+"\n", encoding="utf-8")
    print(f"SUMMARY comparisons={len(results)} receipt={summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

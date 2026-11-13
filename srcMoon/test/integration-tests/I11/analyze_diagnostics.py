#!/usr/bin/env python3
"""Reconstruct I11 block statistics and plot production budget diagnostics.

This processor reads actual linked output only. It reproduces the declared
statistics independently in Python; it neither selects a new startup cutoff
nor produces accepted convergence means or observational validation evidence.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

INTEGRATED = {"global_particle_inventory", "global_surface_return_rate_s"}
STUDENT_T_FIVE = 2.7764451051977987  # t(0.975, df=4), frozen production contract.


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative_change(a, b):
    scale = max(abs(a), abs(b))
    return abs(a-b) / scale if scale else 0.0


def statistics(values, batch_size):
    """Return independently summed lag/equilibrium/CI diagnostics in gate order."""
    count = len(values) // batch_size
    first = len(values) - count * batch_size
    means = [math.fsum(values[first+i*batch_size:first+(i+1)*batch_size]) / batch_size
             for i in range(count)]
    average = math.fsum(means) / count
    deltas = [x-average for x in means]
    variance_sum = math.fsum(x*x for x in deltas)
    lag = (math.fsum(a*b for a, b in zip(deltas, deltas[1:])) / variance_sum
           if variance_sum else 0.0)
    bound = 1.96 / math.sqrt(count)
    result = {"batch_size_base_blocks": batch_size, "batch_count": count,
              "mean_si": average, "lag1": lag, "lag1_bound": bound,
              "independence_pass": abs(lag) <= bound,
              "eligible": count >= 9}
    if count < 9:
        return dict(result, first_failed_gate="insufficient_blocks")
    equilibrium = means[-9:-5]
    changes = [relative_change(a, b) for a, b in zip(equilibrium, equilibrium[1:])]
    science = means[-5:]
    science_mean = math.fsum(science) / 5
    variance = math.fsum((x-science_mean)**2 for x in science) / 4
    half = STUDENT_T_FIVE * math.sqrt(variance / 5)
    # Preserve the production zero-variance convention as a diagnostic, while
    # explicitly flagging zero support for the separately documented audit.
    ci = half / abs(science_mean) if science_mean else (0.0 if half == 0 else math.inf)
    result.update(equilibrium_failed_checks=sum(x > 0.02 for x in changes),
                  equilibrium_max_relative_change=max(changes),
                  science_mean_si=science_mean, relative_ci95_half_width=ci,
                  zero_science_support=not any(science))
    return result


def read_blocks(path):
    """Reject malformed/incomplete blocks rather than silently discard data."""
    blocks = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            block = int(row["base_block"])
            metric = row["metric"]
            value = float(row["value_si"])
            if not math.isfinite(value) or metric in blocks.setdefault(block, {}):
                raise ValueError(f"duplicate/nonfinite metric at block {block}: {path}")
            blocks[block][metric] = value
    if not blocks or sorted(blocks) != list(range(1, len(blocks)+1)):
        raise ValueError(f"missing/nonconsecutive completed blocks: {path}")
    names = set(blocks[1])
    if any(set(item) != names for item in blocks.values()):
        raise ValueError(f"incomplete scoring-vector block: {path}")
    return {name: [blocks[i][name] for i in sorted(blocks)] for name in sorted(names)}


def produce(level, output):
    series = read_blocks(level / "block_metrics.csv")
    count = len(next(iter(series.values())))
    sizes, size = [], 1
    while size <= max(1, count // 9):
        sizes.append(size)
        size *= 2
    records, reconstructed = [], {}
    for size in sizes:
        group = []
        for name, values in series.items():
            row = dict(statistics(values, size), metric=name,
                       classification="integrated" if name in INTEGRATED else "pointwise")
            row["ci_relative_gate"] = 0.0167 if name in INTEGRATED else 0.0333
            if row["eligible"]:
                row["ci_pass"] = row["relative_ci95_half_width"] <= row["ci_relative_gate"]
                row["first_failed_gate"] = (
                    "independence" if not row["independence_pass"] else
                    "equilibrium" if row["equilibrium_failed_checks"] else
                    "confidence_interval" if not row["ci_pass"] else "pass")
            records.append(row)
            group.append(row)
        if count // size >= 9:
            independence = sum(not x["independence_pass"] for x in group)
            equilibrium = sum(x["equilibrium_failed_checks"] for x in group)
            ci = sum(not x["ci_pass"] for x in group)
            reconstructed[size] = {
                "first_failed_gate": "independence" if independence else
                    "equilibrium" if equilibrium else "confidence_interval" if ci else "pass",
                "independence_failed_metrics": independence,
                "equilibrium_failed_checks": equilibrium, "ci_failed_metrics": ci}
    # Compare only the gates actually reached by the production scorer. Later
    # computed diagnostics remain explanatory and never override that score.
    comparisons = []
    with (level / "analysis_diagnostics.csv").open(newline="") as stream:
        for actual in csv.DictReader(stream):
            if int(actual["base_blocks"]) != count:
                continue
            expected = reconstructed[int(actual["batch_size_base_blocks"])]
            matches = actual["first_failed_gate"] == expected["first_failed_gate"]
            for key in ("independence_failed_metrics", "equilibrium_failed_checks", "ci_failed_metrics"):
                if actual[key] != "":
                    matches = matches and int(actual[key]) == expected[key]
            comparisons.append(dict(actual, independent_reconstruction_matches=matches))
    if comparisons and not all(x["independent_reconstruction_matches"] for x in comparisons):
        raise ValueError("independent reconstruction disagrees with production gate diagnostics")
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite diagnostics: {output}")
    output.mkdir(parents=True, exist_ok=True)
    columns = sorted(set().union(*(row.keys() for row in records)))
    with (output / "metric_statistics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader();writer.writerows(records)
    render(series, records, level, output)
    inputs = [level / "block_metrics.csv", level / "analysis_diagnostics.csv"]
    if (level / "particle_budget.csv").exists():
        inputs.append(level / "particle_budget.csv")
    report = {"schema": "srcMoon-I11-independent-diagnostics-v1",
              "formal_i11_status": "NOT EVALUATED BY THIS PROCESSOR",
              "accepted_convergence_input": False, "observational_validation": False,
              "metric_count": len(series), "base_blocks": count,
              "gate_reconstruction": "PASS" if comparisons else "SKIPPED",
              "production_comparisons": comparisons,
              "source_files": {str(p): digest(p) for p in inputs},
              "processor_sha256": digest(Path(__file__)),
              "artifacts": {p.name: digest(p) for p in output.iterdir() if p.is_file()}}
    (output / "diagnostics.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def render(series, records, level, output):
    """Figures display actual numerical output; CSV retains every scored bin."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    first = [r for r in records if r["batch_size_base_blocks"] == 1]
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    for name in sorted(INTEGRATED.intersection(series)):
        values = series[name]
        scale = max(abs(x) for x in values) or 1
        axes[0].plot(range(1, len(values)+1), [x/scale for x in values], marker="o", label=name)
    axes[0].set(xlabel="500-step sample index", ylabel="value / maximum absolute value",
                title="Integrated observables: numerical histories")
    axes[0].legend()
    axes[1].plot([r["lag1"] for r in first], ".", markersize=2)
    bound = first[0]["lag1_bound"]
    axes[1].axhline(bound, color="red", linestyle="--")
    axes[1].axhline(-bound, color="red", linestyle="--")
    axes[1].set(xlabel="metric index (alphabetical; see CSV)", ylabel="lag-one correlation",
                title="All metrics and the unchanged independence bound")
    fig.suptitle("I11 DIAGNOSTICS ONLY — NO CONVERGENCE OR OBSERVATIONAL VALIDATION")
    for suffix in ("png", "eps"):
        fig.savefig(output / f"statistics.{suffix}", dpi=160)
    plt.close(fig)
    budget = level / "particle_budget.csv"
    if not budget.exists():
        return
    with budget.open(newline="") as stream:
        data = list(csv.DictReader(stream))
    times = [float(r["elapsed_particle_time_s"])/3600 for r in data]
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), constrained_layout=True)
    axes[0].plot(times, [float(r["gas_atoms"]) for r in data], "o-")
    axes[0].set(ylabel="instantaneous Na gas inventory [atoms]")
    # A cooperative stop usually leaves a shorter final checkpoint interval.
    # Plot actual rates [atoms/s], using each interval's elapsed particle time,
    # so that a short closing row cannot look like a sudden source/sink drop.
    # Keep every interval and its raw atom totals in the input ledger; this
    # conversion does not resample, smooth, or change any scoring quantity.
    seconds = [float(r["elapsed_particle_time_s"]) for r in data]
    intervals = [b-a for a, b in zip(seconds, seconds[1:])]
    if any(not math.isfinite(dt) or dt <= 0 for dt in intervals):
        raise ValueError("budget rate plot requires strictly increasing checkpoint times [s]")
    for name in ("emission", "adsorption", "photolytic_loss", "escape", "numerical_removal"):
        rates = [float(r[name+"_atoms"])/dt for r, dt in zip(data[1:], intervals)]
        axes[1].plot(times[1:], rates, marker=".", label=name)
    axes[1].set(ylabel="actual weighted event rate [atoms/s]")
    axes[1].legend()
    axes[2].plot(times, [float(r["weighted_residual_atoms"]) for r in data], "o-", label="weighted residual")
    axes[2].plot(times, [float(r["weighted_roundoff_bound_atoms"]) for r in data], "--", label="roundoff bound")
    axes[2].set(xlabel="elapsed particle time [hours]; external epoch frozen", ylabel="atoms", yscale="symlog")
    axes[2].legend()
    fig.suptitle("ACTUAL-EVENT Na GAS BUDGET — DIAGNOSTIC, NOT I19 QUALIFICATION")
    for suffix in ("png", "eps"):
        fig.savefig(output / f"particle_budget.{suffix}", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    report = produce(args.level_dir.resolve(), (args.output_dir or args.level_dir / "diagnostics").resolve())
    print(json.dumps(report, indent=2, sort_keys=True))

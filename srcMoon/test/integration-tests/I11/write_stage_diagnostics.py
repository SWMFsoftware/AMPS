#!/usr/bin/env python3
"""Render raw I11 stage diagnostics without creating validation evidence.

The accepted I11 convergence calculation uses only ``qualified_metrics.csv``
after the production independence, equilibrium, and confidence-interval gates
pass.  This renderer deliberately reads ``block_metrics.csv`` instead.  Its
CSV and figures help an investigator inspect a completed resolution stage,
including a capped FAIL, but they are never inputs to the formal I11 score.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


INTEGRATED_METRICS = {
    "global_surface_return_rate_s",
    "global_particle_inventory",
}
NOTICE = "DIAGNOSTIC ONLY - NOT OBSERVATIONAL VALIDATION"


def sha256(path: Path) -> str:
    """Return a streaming digest suitable for the stage provenance record."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_change(left: float, right: float) -> float:
    """Use I11's symmetric scale so zero-valued metrics remain well defined."""
    scale = max(abs(left), abs(right))
    return 0.0 if scale == 0.0 else abs(left - right) / scale


def signed_relative_deviation(value: float, mean: float) -> float:
    """Normalize a raw sample without dividing by a possibly zero mean."""
    scale = max(abs(value), abs(mean))
    return 0.0 if scale == 0.0 else (value - mean) / scale


def read_blocks(path: Path) -> tuple[list[int], dict[str, list[float]]]:
    """Read a complete rectangular production block table strictly.

    A truncated or duplicated row is infrastructure corruption, not an
    invitation to silently drop a metric from a visually convenient plot.
    """
    by_block: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected = ["base_block", "production_output_number", "metric",
                    "value_si"]
        if reader.fieldnames != expected:
            raise ValueError(f"unexpected block-metric header in {path}")
        for line_number, row in enumerate(reader, 2):
            try:
                block = int(row["base_block"])
                value = float(row["value_si"])
            except ValueError as exception:
                raise ValueError(
                    f"line {line_number} has a non-numeric block/value") from exception
            metric = row["metric"]
            if block <= 0 or not metric or not math.isfinite(value):
                raise ValueError(f"line {line_number} is not finite/identified")
            values = by_block.setdefault(block, {})
            if metric in values:
                raise ValueError(
                    f"line {line_number} duplicates block {block} metric {metric}")
            values[metric] = value
    if not by_block:
        raise ValueError(f"no completed production samples in {path}")

    blocks = sorted(by_block)
    if blocks != list(range(1, len(blocks) + 1)):
        raise ValueError(f"non-contiguous base-block sequence in {path}")
    names = set(by_block[blocks[0]])
    if not names:
        raise ValueError(f"empty scoring vector in {path}")
    for block in blocks[1:]:
        if set(by_block[block]) != names:
            raise ValueError(f"scoring-vector shape changed at block {block}")
    return blocks, {
        name: [by_block[block][name] for block in blocks]
        for name in sorted(names)
    }


def summarize_level(resolution_km: int, path: Path) -> dict[str, object]:
    """Summarize raw blocks; no equilibrium or CI qualification is inferred."""
    blocks, values = read_blocks(path)
    metrics: dict[str, dict[str, object]] = {}
    for name, samples in values.items():
        mean = statistics.fmean(samples)
        standard_deviation = statistics.stdev(samples) if len(samples) > 1 else 0.0
        metrics[name] = {
            "classification": (
                "integrated" if name in INTEGRATED_METRICS else "pointwise"),
            "raw_mean_si": mean,
            "raw_standard_deviation_si": standard_deviation,
            "relative_standard_deviation": (
                0.0 if mean == 0.0 and standard_deviation == 0.0 else
                math.inf if mean == 0.0 else standard_deviation / abs(mean)),
            "samples": samples,
        }
    return {
        "resolution_km": resolution_km,
        "block_metrics": str(path),
        "base_blocks": blocks,
        "metrics": metrics,
    }


def write_csv(levels: list[dict[str, object]], destination: Path) -> None:
    """Preserve every completed level mean and its raw predecessor change."""
    fields = [
        "resolution_km", "metric", "classification", "base_sample_count",
        "raw_mean_si", "raw_standard_deviation_si",
        "relative_standard_deviation", "previous_resolution_km",
        "previous_raw_mean_si", "relative_change_from_previous",
        "evidence_class",
    ]
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        previous: dict[str, object] | None = None
        for level in levels:
            metrics = level["metrics"]
            assert isinstance(metrics, dict)
            previous_metrics = None if previous is None else previous["metrics"]
            for name, metric_value in metrics.items():
                metric = metric_value
                assert isinstance(metric, dict)
                old = (None if previous_metrics is None else
                       previous_metrics.get(name))  # type: ignore[union-attr]
                writer.writerow({
                    "resolution_km": level["resolution_km"],
                    "metric": name,
                    "classification": metric["classification"],
                    "base_sample_count": len(level["base_blocks"]),
                    "raw_mean_si": f"{float(metric['raw_mean_si']):.17g}",
                    "raw_standard_deviation_si":
                        f"{float(metric['raw_standard_deviation_si']):.17g}",
                    "relative_standard_deviation":
                        f"{float(metric['relative_standard_deviation']):.17g}",
                    "previous_resolution_km": (
                        "" if previous is None else previous["resolution_km"]),
                    "previous_raw_mean_si": (
                        "" if old is None else
                        f"{float(old['raw_mean_si']):.17g}"),
                    "relative_change_from_previous": (
                        "" if old is None else f"{relative_change(float(old['raw_mean_si']), float(metric['raw_mean_si'])):.17g}"),
                    "evidence_class": "unqualified-raw-diagnostic",
                })
            previous = level


def quantile(values: list[float], fraction: float) -> float:
    """Return a deterministic linearly interpolated quantile for plotting."""
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def render(levels: list[dict[str, object]], current_status: str,
           png: Path, eps: Path) -> None:
    """Render within-stage stability and cumulative raw grid comparisons."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    current = levels[-1]
    current_metrics = current["metrics"]
    assert isinstance(current_metrics, dict)
    blocks = current["base_blocks"]
    assert isinstance(blocks, list)
    pointwise = [metric for metric in current_metrics.values()
                 if metric["classification"] == "pointwise"]
    integrated = [(name, metric) for name, metric in current_metrics.items()
                  if metric["classification"] == "integrated"]

    figure, axes = plt.subplots(3, 1, figsize=(9.0, 11.0),
                                constrained_layout=True)
    # Hundreds of profile/LOS metrics cannot be read as individual curves.
    # Quantile envelopes retain the distribution of their normalized sample
    # motion without selecting only those metrics that look well behaved.
    medians: list[float] = []
    p10: list[float] = []
    p90: list[float] = []
    for sample_index in range(len(blocks)):
        deviations = [abs(signed_relative_deviation(
            float(metric["samples"][sample_index]),
            float(metric["raw_mean_si"]))) for metric in pointwise]
        medians.append(quantile(deviations, 0.5))
        p10.append(quantile(deviations, 0.1))
        p90.append(quantile(deviations, 0.9))
    axes[0].plot(blocks, medians, marker="o", label="median")
    axes[0].plot(blocks, p10, marker=".", linestyle="--", label="10th percentile")
    axes[0].plot(blocks, p90, marker=".", linestyle="--", label="90th percentile")
    axes[0].set_ylabel("absolute relative deviation\nfrom raw stage mean")
    axes[0].set_xlabel("500-step production sample")
    axes[0].set_title("Pointwise metric stability (all declared metrics)")
    axes[0].grid(True, linewidth=0.4, color="#b0b0b0")
    axes[0].legend(framealpha=1.0)

    for name, metric in integrated:
        deviations = [signed_relative_deviation(float(value),
                      float(metric["raw_mean_si"]))
                      for value in metric["samples"]]
        axes[1].plot(blocks, deviations, marker="o", label=name)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_ylabel("signed relative deviation\nfrom raw stage mean")
    axes[1].set_xlabel("500-step production sample")
    axes[1].set_title("Integrated metric stability")
    axes[1].grid(True, linewidth=0.4, color="#b0b0b0")
    axes[1].legend(framealpha=1.0)

    if len(levels) == 1:
        axes[2].text(0.5, 0.5,
                     "First resolution stage: no preceding grid exists.\n"
                     "Cross-resolution comparison begins after stage 2.",
                     ha="center", va="center", transform=axes[2].transAxes)
        axes[2].set_axis_off()
    else:
        previous = levels[-2]["metrics"]
        changes = {"pointwise": [], "integrated": []}
        for name, metric in current_metrics.items():
            old = previous[name]
            changes[str(metric["classification"])].append(relative_change(
                float(old["raw_mean_si"]), float(metric["raw_mean_si"])))
        for classification, color in (("pointwise", "tab:blue"),
                                      ("integrated", "tab:orange")):
            ordered = sorted(changes[classification])
            axes[2].plot(range(1, len(ordered) + 1), ordered, marker=".",
                         linestyle="none", color=color,
                         label=f"{classification} raw means")
        axes[2].axhline(0.10, color="tab:blue", linestyle="--",
                       label="formal pointwise tolerance (context only)")
        axes[2].axhline(0.05, color="tab:orange", linestyle="--",
                       label="formal integrated tolerance (context only)")
        axes[2].set_ylabel("relative change from preceding grid")
        axes[2].set_xlabel("metrics sorted independently within each class")
        axes[2].set_title(
            f"Raw means: {levels[-2]['resolution_km']} km vs "
            f"{current['resolution_km']} km")
        axes[2].grid(True, linewidth=0.4, color="#b0b0b0")
        axes[2].legend(framealpha=1.0)

    figure.suptitle(
        f"I11 {current['resolution_km']} km stage; formal level status: "
        f"{current_status}\n{NOTICE}", fontsize=12)
    figure.savefig(png, dpi=180, transparent=False)
    figure.savefig(eps, format="eps", transparent=False)
    plt.close(figure)


def produce(campaign: Path, resolutions: list[int], current_status: str) -> dict[str, object]:
    """Create one immutable diagnostic set after the latest completed stage."""
    if not resolutions:
        raise ValueError("at least one completed resolution is required")
    levels = [summarize_level(
        resolution, campaign / f"{resolution}km" / "block_metrics.csv")
        for resolution in resolutions]
    destination = campaign / f"{resolutions[-1]}km"
    csv_path = destination / "stage_diagnostic.csv"
    png_path = destination / "stage_diagnostic.png"
    eps_path = destination / "stage_diagnostic.eps"
    json_path = destination / "stage_diagnostic.json"
    write_csv(levels, csv_path)
    render(levels, current_status, png_path, eps_path)
    metadata = {
        "schema": "srcMoon-I11-stage-diagnostic-v1",
        "notice": NOTICE,
        "formal_level_status": current_status,
        "current_resolution_km": resolutions[-1],
        "completed_resolutions_km": resolutions,
        "sample_length_particle_steps": 500,
        "accepted_convergence_input": False,
        "observational_validation": False,
        "source_block_metrics": [
            {"path": str(level["block_metrics"]),
             "sha256": sha256(Path(str(level["block_metrics"])))}
            for level in levels
        ],
        "artifacts": {
            "csv": str(csv_path), "csv_sha256": sha256(csv_path),
            "png": str(png_path), "png_sha256": sha256(png_path),
            "eps": str(eps_path), "eps_sha256": sha256(eps_path),
        },
    }
    json_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    metadata["artifacts"]["metadata"] = str(json_path)  # type: ignore[index]
    metadata["artifacts"]["metadata_sha256"] = sha256(json_path)  # type: ignore[index]
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render unqualified raw diagnostics after one I11 stage")
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--completed-resolutions-km", type=int, nargs="+",
                        required=True)
    parser.add_argument("--formal-level-status",
                        choices=("PASS", "FAIL", "SKIPPED", "ERROR"),
                        required=True)
    arguments = parser.parse_args()
    try:
        metadata = produce(arguments.campaign_dir.resolve(),
                           arguments.completed_resolutions_km,
                           arguments.formal_level_status)
    except (OSError, ValueError, ImportError) as exception:
        print(f"ERROR: {exception}")
        return 2
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

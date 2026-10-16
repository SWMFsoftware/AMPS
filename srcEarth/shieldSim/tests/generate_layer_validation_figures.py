#!/usr/bin/env python3
"""Generate figures for analytic/controlled comparisons in legacy layer tests.

The shell runners perform the numerical acceptance checks before invoking this
script.  Figures are review aids and are never compared pixel-by-pixel.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import sys
import textwrap

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as error:  # pragma: no cover - environment dependent
    raise SystemExit(f"layer validation figures require matplotlib: {error}")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric_rows(path: pathlib.Path) -> list[list[str]]:
    require(path.is_file(), f"missing layer-test artifact: {path}")
    rows = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if stripped and not stripped.startswith("#"):
            rows.append(stripped.split())
    require(bool(rows), f"no numeric rows in {path}")
    return rows


def save(fig, output: pathlib.Path, stem: str, caption: str, entries: list[dict]) -> None:
    fig.text(0.01, 0.008, textwrap.fill(caption, width=145), fontsize=8, color="#333333")
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.92))
    for extension in ("png", "svg"):
        path = output / f"{stem}.{extension}"
        kwargs = {"bbox_inches": "tight"}
        if extension == "png":
            kwargs["dpi"] = 170
        fig.savefig(path, **kwargs)
        require(path.stat().st_size > 1000, f"figure is unexpectedly small: {path}")
        entries.append(
            {
                "file": path.name,
                "caption": caption,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    plt.close(fig)


def layer2(run_dir: pathlib.Path, output: pathlib.Path, entries: list[dict]) -> None:
    beam = numeric_rows(run_dir / "beam_source" / "source.dat")
    beam_energy = [float(row[2]) for row in beam]
    beam_x = [float(row[3]) for row in beam]
    beam_y = [float(row[4]) for row in beam]
    beam_uz = [float(row[8]) for row in beam]
    events = list(range(len(beam)))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    axes[0].scatter(beam_x, beam_y, s=8, alpha=0.45, label="generated")
    axes[0].scatter([0.0], [0.0], marker="+", s=120, color="#d6604d", label="analytic beam")
    axes[0].set_xlabel("source x [mm]")
    axes[0].set_ylabel("source y [mm]")
    axes[0].legend(fontsize=8)
    axes[1].plot(events, beam_energy, ".", markersize=2, label="generated")
    axes[1].axhline(100.0, color="#d6604d", linestyle="--", label="analytic 100 MeV")
    axes[1].set_xlabel("event")
    axes[1].set_ylabel("kinetic energy [MeV]")
    axes[1].legend(fontsize=8)
    axes[2].plot(events, beam_uz, ".", markersize=2, label="generated")
    axes[2].axhline(1.0, color="#d6604d", linestyle="--", label="analytic +z")
    axes[2].set_xlabel("event")
    axes[2].set_ylabel("direction cosine uz")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.suptitle("Layer 2 pencil-beam analytic source contract", fontsize=14, fontweight="bold")
    save(
        fig, output, "layer2_beam_source_analytic",
        "The generated pencil-beam position, energy, and direction are compared with the exact "
        "monoenergetic on-axis source contract.", entries,
    )

    isotropic = numeric_rows(run_dir / "isotropic_source" / "source.dat")
    iso_x = [float(row[3]) for row in isotropic]
    iso_y = [float(row[4]) for row in isotropic]
    mu = [float(row[8]) for row in isotropic]
    mean_mu = sum(mu) / len(mu)
    mean_mu2 = sum(value * value for value in mu) / len(mu)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    axes[0].hist(mu, bins=40, range=(0.0, 1.0), density=True, alpha=0.65,
                 color="#2166ac", label="generated")
    grid = [index / 200.0 for index in range(201)]
    axes[0].plot(grid, [2.0 * value for value in grid], color="#d6604d", linewidth=2,
                 label="analytic p(µ)=2µ")
    axes[0].set_xlabel("µ = uz")
    axes[0].set_ylabel("probability density")
    axes[0].legend(fontsize=8)
    axes[0].text(0.04, 0.72, f"<µ>={mean_mu:.5f} (2/3)\n<µ²>={mean_mu2:.5f} (1/2)",
                 transform=axes[0].transAxes, fontsize=8)
    axes[1].hist(iso_x, bins=40, range=(-25.0, 25.0), density=True, alpha=0.65,
                 color="#2166ac", label="generated x")
    axes[1].axhline(1.0 / 50.0, color="#d6604d", linestyle="--",
                    label="uniform plane density")
    axes[1].set_xlabel("source x [mm]")
    axes[1].set_ylabel("probability density [1/mm]")
    axes[1].legend(fontsize=8)
    stride = max(1, len(iso_x) // 3000)
    axes[2].scatter(iso_x[::stride], iso_y[::stride], s=2, alpha=0.25)
    axes[2].set_xlim(-26, 26)
    axes[2].set_ylim(-26, 26)
    axes[2].set_aspect("equal", adjustable="box")
    axes[2].set_xlabel("source x [mm]")
    axes[2].set_ylabel("source y [mm]")
    for axis in axes:
        axis.grid(alpha=0.22)
    fig.suptitle("Layer 2 isotropic plane-crossing source", fontsize=14, fontweight="bold")
    save(
        fig, output, "layer2_isotropic_source_analytic",
        "The sampled plane-crossing direction distribution is compared with p(µ)=2µ and its "
        "analytic moments; the source-plane coordinates are checked against the uniform finite plane.",
        entries,
    )

    vacuum = numeric_rows(run_dir / "vacuum_transmission" / "exit.dat")
    exit_energy = [float(row[2]) for row in vacuum]
    local_z = [float(row[8]) for row in vacuum]
    local_uz = [float(row[11]) for row in vacuum]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    axes[0].hist(exit_energy, bins=35, color="#2166ac", alpha=0.7)
    axes[0].axvline(100.0, color="#d6604d", linestyle="--", label="incident 100 MeV")
    axes[0].set_xlabel("rear-face proton energy [MeV]")
    axes[0].set_ylabel("count")
    axes[0].legend(fontsize=8)
    axes[1].hist(local_z, bins=35, color="#2166ac", alpha=0.7)
    axes[1].axvline(0.0005, color="#d6604d", linestyle="--", label="analytic rear face")
    axes[1].set_xlabel("local rear-face z [mm]")
    axes[1].set_ylabel("count")
    axes[1].legend(fontsize=8)
    axes[2].hist(local_uz, bins=35, color="#2166ac", alpha=0.7)
    axes[2].axvline(1.0, color="#d6604d", linestyle="--", label="incident +z")
    axes[2].set_xlabel("local exit direction uz")
    axes[2].set_ylabel("count")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.22)
    fig.suptitle("Layer 2 near-vacuum transmission check", fontsize=14, fontweight="bold")
    save(
        fig, output, "layer2_vacuum_transmission_analytic",
        "Rear-face energy, local crossing position, and downstream direction are shown against the "
        "near-vacuum limiting expectations. This is a geometry/scoring check, not material validation.",
        entries,
    )


def read_summary(path: pathlib.Path) -> tuple[dict[str, float], list[list[str]]]:
    scalars: dict[str, float] = {}
    targets: list[list[str]] = []
    for fields in numeric_rows(path):
        if fields[0] == "scalar":
            scalars[fields[1]] = float(fields[2])
        elif fields[0] == "target":
            targets.append(fields)
    return scalars, targets


def mean_let(path: pathlib.Path) -> tuple[list[float], list[float], float]:
    let: list[float] = []
    charged: list[float] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith(("#", "TITLE", "VARIABLES", "ZONE")):
            continue
        fields = stripped.split()
        if len(fields) >= 4:
            let.append(float(fields[0]))
            charged.append(float(fields[3]))
    denominator = sum(charged)
    require(denominator > 0.0, f"no charged LET fluence in {path}")
    return let, charged, sum(x * weight for x, weight in zip(let, charged)) / denominator


def layer3(run_dir: pathlib.Path, output: pathlib.Path, entries: list[dict]) -> None:
    scalar1, targets1 = read_summary(run_dir / "rate_scale_1" / "summary.dat")
    scalar10, targets10 = read_summary(run_dir / "rate_scale_10" / "summary.dat")
    require(targets1 and targets10, "source scaling summaries lack targets")
    tid1, rate1 = float(targets1[0][4]), float(targets1[0][5])
    tid10, rate10 = float(targets10[0][4]), float(targets10[0][5])
    expected = [1.0, 10.0, 10.0]
    actual = [tid10 / tid1, rate10 / rate1, scalar10["source_norm"] / scalar1["source_norm"]]
    labels = ["per-primary TID", "TID rate", "source normalization"]
    positions = list(range(len(labels)))
    fig, axis = plt.subplots(figsize=(9, 5.6))
    axis.bar([p - 0.18 for p in positions], expected, width=0.36, color="#bdbdbd", label="analytic")
    axis.bar([p + 0.18 for p in positions], actual, width=0.36, color="#2166ac", label="computed")
    axis.set_xticks(positions)
    axis.set_xticklabels(labels)
    axis.set_ylabel("rate-10 / rate-1 ratio")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    for position, value in zip(positions, actual):
        axis.text(position + 0.18, value, f"{value:.5g}", ha="center", va="bottom", fontsize=8)
    fig.suptitle("Layer 3 source-normalization analytic scaling", fontsize=14, fontweight="bold")
    save(
        fig, output, "layer3_source_normalization_analytic",
        "A tenfold source intensity must leave per-primary TID unchanged while source-normalized "
        "rate and recorded normalization increase tenfold.", entries,
    )

    low, _ = read_summary(run_dir / "hardness_50MeV" / "summary.dat")
    high, _ = read_summary(run_dir / "hardness_150MeV" / "summary.dat")
    hardness_actual = [low["H100_10"], high["H100_10"]]
    fig, axis = plt.subplots(figsize=(8, 5.6))
    axis.bar([p - 0.18 for p in (0, 1)], [0.0, 1.0], width=0.36,
             color="#bdbdbd", label="ideal threshold limit")
    axis.bar([p + 0.18 for p in (0, 1)], hardness_actual, width=0.36,
             color="#2166ac", label="computed")
    axis.set_xticks([0, 1])
    axis.set_xticklabels(["50 MeV proton", "150 MeV proton"])
    axis.set_ylim(0.0, 1.08)
    axis.set_ylabel("H100/10")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    fig.suptitle("Layer 3 H100/10 analytic threshold limits", fontsize=14, fontweight="bold")
    save(
        fig, output, "layer3_hardness_threshold_analytic",
        "Near-vacuum monoenergetic protons below and above 100 MeV are compared with the ideal "
        "primary-proton H100/10 limits of zero and one.", entries,
    )

    proton_x, proton_y, proton_mean = mean_let(
        run_dir / "let_proton_100MeV" / "letp_let_spectrum.dat"
    )
    alpha_x, alpha_y, alpha_mean = mean_let(
        run_dir / "let_alpha_100MeV" / "leta_let_spectrum.dat"
    )
    fig, (spectrum, means) = plt.subplots(1, 2, figsize=(12, 5.6))
    spectrum.plot(proton_x, proton_y, label="100 MeV proton", color="#2166ac")
    spectrum.plot(alpha_x, alpha_y, label="100 MeV total-energy alpha", color="#d6604d")
    spectrum.set_xscale("log")
    spectrum.set_yscale("log")
    spectrum.set_xlabel("LET [MeV cm²/mg]")
    spectrum.set_ylabel("charged fluence per LET bin [test output units]")
    spectrum.grid(which="both", alpha=0.22)
    spectrum.legend(fontsize=8)
    means.bar([0, 1], [proton_mean, alpha_mean], color=["#2166ac", "#d6604d"])
    means.set_xticks([0, 1])
    means.set_xticklabels(["proton", "alpha"])
    means.set_ylabel("fluence-weighted mean LET [MeV cm²/mg]")
    means.grid(axis="y", alpha=0.22)
    means.text(0.5, 0.92, f"alpha/proton = {alpha_mean/proton_mean:.3g}",
               transform=means.transAxes, ha="center")
    fig.suptitle("Layer 3 controlled LET comparison", fontsize=14, fontweight="bold")
    save(
        fig, output, "layer3_let_controlled_comparison",
        "The controlled proton/alpha LET sanity test is visualized. It is not an analytic solution "
        "or observational validation and remains a prototype diagnostic on the legacy path.", entries,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("layer", choices=("layer2", "layer3"))
    parser.add_argument("--run-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    entries: list[dict] = []
    if args.layer == "layer2":
        layer2(run_dir, output, entries)
        expected_subjects = 3
    else:
        layer3(run_dir, output, entries)
        expected_subjects = 3
    require(len(entries) == 2 * expected_subjects, "layer figure set is incomplete")
    manifest = {
        "schema_version": "shieldSim.validation_figure_manifest.v1",
        "authority": "review_aid_only_numeric_tests_are_authoritative",
        "physical_validation_claim": False,
        "layer": args.layer,
        "figures": sorted(entries, key=lambda item: item["file"]),
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Generated {expected_subjects} {args.layer} validation figures in PNG and SVG: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError) as error:
        print(f"FAIL layer validation figure generation: {error}", file=sys.stderr)
        raise SystemExit(1)

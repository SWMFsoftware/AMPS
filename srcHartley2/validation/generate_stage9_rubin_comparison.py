#!/usr/bin/env python3
"""Publish authenticated Stage-9 ensemble profiles beside Rubin et al. (2011).

The input is the checksum-protected ``RUN_MANIFEST.json`` produced by the
Stage-9 campaign driver.  This utility deliberately reuses the strict science
gate before reading any profile: incomplete seed matrices, altered summaries,
different physical hashes, nonconverged runs, or failed Stage-9 acceptance
cannot be turned into a presentation plot by this program.

Number densities in the Stage-9 summaries are SI (m^-3).  Rubin's digitized
neutral curves are cm^-3, so every model value and confidence limit is divided
by 10^6.  The CSV and JSON retain both estimator arms and the published curve
variants; the graphics show the independent-seed mean and confidence bounds,
not one visually favorable seed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
TESTS_DIR = REPO_ROOT / "srcHartley2" / "tests"
sys.path.insert(0, str(TESTS_DIR))

import test_stage9_science_gate as science_gate  # noqa: E402


SCHEMA = "amps.comet.stage9-rubin-comparison.v1"
DEFAULT_PAPER_DATA = (
    REPO_ROOT
    / "srcHartley2"
    / "halley2011"
    / "reference"
    / "rubin2011_neutral_profiles.csv"
)
DEFAULT_PAPER_METADATA = DEFAULT_PAPER_DATA.with_suffix(".metadata.json")
PAPER_COLUMNS = {
    "figure",
    "pdf_page",
    "species",
    "variant",
    "quantity",
    "radius_km",
    "value",
    "unit",
}


class ComparisonError(RuntimeError):
    """Raised when an archive or paper reference is not publishable evidence."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_object(path: Path, label: str) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ComparisonError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ComparisonError(f"{label} root is not an object: {path}")
    return value


def finite(raw: str, label: str, *, positive: bool = False) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ComparisonError(f"{label} is not numeric") from exc
    if not math.isfinite(value) or (positive and value <= 0.0):
        raise ComparisonError(f"{label} is not finite and positive")
    return value


PaperKey = Tuple[int, str, str]


def load_paper_profiles(
    csv_path: Path, metadata_path: Path
) -> Tuple[Dict[PaperKey, List[Dict[str, float]]], Dict[str, Any]]:
    """Load and cross-check the committed vector-path digitization."""

    metadata = load_object(metadata_path, "Rubin metadata")
    if not isinstance(metadata.get("curves"), list) or not metadata["curves"]:
        raise ComparisonError("Rubin metadata has no curve inventory")
    source_hash = metadata.get("source_pdf_sha256")
    if (
        not isinstance(source_hash, str)
        or len(source_hash) != 64
        or any(character not in "0123456789abcdef" for character in source_hash)
    ):
        raise ComparisonError("Rubin metadata has no valid source-PDF SHA-256")

    grouped: Dict[PaperKey, List[Dict[str, float]]] = defaultdict(list)
    try:
        with csv_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if set(reader.fieldnames or ()) != PAPER_COLUMNS:
                raise ComparisonError("Rubin CSV columns differ from the committed schema")
            for line_number, row in enumerate(reader, start=2):
                if row["quantity"] != "number_density":
                    continue
                if row["unit"] != "cm^-3":
                    raise ComparisonError(
                        f"Rubin density row {line_number} has unit {row['unit']!r}"
                    )
                try:
                    figure = int(row["figure"])
                except ValueError as exc:
                    raise ComparisonError(
                        f"Rubin row {line_number} has an invalid figure"
                    ) from exc
                species, variant = row["species"], row["variant"]
                if not species or not variant:
                    raise ComparisonError(
                        f"Rubin row {line_number} lacks species or variant"
                    )
                grouped[(figure, species, variant)].append(
                    {
                        "radius_km": finite(
                            row["radius_km"], f"Rubin radius row {line_number}", positive=True
                        ),
                        "density_cm-3": finite(
                            row["value"], f"Rubin density row {line_number}", positive=True
                        ),
                    }
                )
    except OSError as exc:
        raise ComparisonError(f"cannot read Rubin CSV {csv_path}: {exc}") from exc
    if not grouped:
        raise ComparisonError("Rubin CSV contains no neutral-density curves")

    inventory: Dict[PaperKey, Mapping[str, Any]] = {}
    for row in metadata["curves"]:
        if not isinstance(row, dict) or row.get("quantity") != "number_density":
            continue
        try:
            key = (int(row["figure"]), str(row["species"]), str(row["variant"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise ComparisonError("Rubin metadata contains a malformed curve") from exc
        if key in inventory:
            raise ComparisonError(f"duplicate Rubin metadata curve {key}")
        inventory[key] = row
    if set(grouped) != set(inventory):
        raise ComparisonError("Rubin CSV and metadata curve inventories differ")
    for key, rows in grouped.items():
        last = -math.inf
        for row in rows:
            if row["radius_km"] <= last:
                raise ComparisonError(f"Rubin curve {key} is not strictly radial")
            last = row["radius_km"]
        expected_samples = inventory[key].get("samples")
        if expected_samples != len(rows):
            raise ComparisonError(
                f"Rubin curve {key} has {len(rows)} samples, metadata says {expected_samples}"
            )
    return dict(grouped), metadata


def authenticated_runs(manifest_path: Path) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """Return raw summaries only after the science gate authenticated them."""

    manifest = load_object(manifest_path, "Stage-9 manifest")
    runs: List[Dict[str, Any]] = []
    for record in manifest.get("runs", []):
        # The gate has already checked path containment and the digest.  We
        # repeat neither rule differently here; this second read supplies
        # only convergence/exposure metadata omitted from the compact result.
        path = (manifest_path.parent / record["summary"]).resolve()
        summary = load_object(path, "Stage-9 summary")
        runs.append(
            {
                "method": summary["method"],
                "seed": summary["seed"],
                "mpi_ranks": summary["mpi_ranks"],
                "omp_threads": summary["omp_threads"],
                "complete_steady_windows": summary["convergence"][
                    "complete_steady_windows"
                ],
                "discarded_burn_in_iterations": summary["convergence"][
                    "discarded_burn_in_iterations"
                ],
                "convergence_status": summary["convergence"]["status"],
                "hard_limit_reached": summary["convergence"]["hard_limit_reached"],
                "wall_seconds": summary["operations"]["wall_seconds"],
                "summary": record["summary"],
                "summary_sha256": record["sha256"],
            }
        )
    return manifest, runs


def build_product(
    manifest_path: Path,
    paper_csv: Path,
    paper_metadata: Path,
) -> Dict[str, Any]:
    """Authenticate inputs and form the unit-converted comparison object."""

    manifest_path = manifest_path.resolve()
    try:
        science = science_gate.evaluate(manifest_path)
    except (science_gate.EvidenceError, OSError) as exc:
        raise ComparisonError(f"Stage-9 evidence is invalid: {exc}") from exc
    if not science.get("accepted"):
        failures = science.get("failures", [])
        detail = "; ".join(str(item) for item in failures[:5])
        raise ComparisonError(
            "Stage-9 science gate did not pass; comparison publication is refused"
            + (f": {detail}" if detail else "")
        )
    manifest, runs = authenticated_runs(manifest_path)
    paper, paper_info = load_paper_profiles(paper_csv.resolve(), paper_metadata.resolve())

    profiles: Dict[str, Any] = {}
    for species, methods in science["profiles"].items():
        profile: Dict[str, Any] = {"ANALOG": [], "FORCED": [], "RUBIN2011": []}
        for method in science_gate.METHODS:
            for row in methods[method]:
                profile[method].append(
                    {
                        "radius_km": row["radius_m"] / 1000.0,
                        "mean_density_cm-3": row["mean"] * 1.0e-6,
                        "standard_error_cm-3": row["standard_error"] * 1.0e-6,
                        "confidence_low_cm-3": row["confidence_low"] * 1.0e-6,
                        "confidence_high_cm-3": row["confidence_high"] * 1.0e-6,
                    }
                )
        for (figure, paper_species, variant), rows in sorted(paper.items()):
            if paper_species == species:
                profile["RUBIN2011"].append(
                    {"figure": figure, "variant": variant, "samples": rows}
                )
        profiles[species] = profile

    policy = manifest["policy"]
    return {
        "schema": SCHEMA,
        "stage9_outcome": science["outcome"],
        "statistical_method": (
            "independent-seed arithmetic mean and policy-supplied symmetric "
            "confidence multiplier times the sample standard error"
        ),
        "confidence_multiplier": policy["confidence_multiplier"],
        "seeds": science["seeds"],
        "run_count": science["run_count"],
        "density_units": "cm^-3",
        "radius_units": "km cometocentric",
        "effective_exposure": {
            "unit": "complete post-burn-in sampling windows",
            "runs": runs,
        },
        "convergence": {
            "all_runs_converged_before_hard_limit": all(
                row["convergence_status"] == "converged"
                and not row["hard_limit_reached"]
                for row in runs
            ),
            "minimum_complete_steady_windows": policy["minimum_steady_windows"],
            "minimum_discarded_burn_in_iterations": policy["burn_in_iterations"],
        },
        "provenance": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256_file(manifest_path),
            "run_summary_sha256": {
                f"{row['method']}/seed={row['seed']}": row["summary_sha256"]
                for row in runs
            },
            "physical_input_sha256": sorted(
                {
                    load_object(
                        (manifest_path.parent / row["summary"]).resolve(),
                        "Stage-9 summary",
                    )["physical_input_sha256"]
                    for row in runs
                }
            )[0],
            "paper_csv": str(paper_csv.resolve()),
            "paper_csv_sha256": sha256_file(paper_csv.resolve()),
            "paper_metadata": str(paper_metadata.resolve()),
            "paper_metadata_sha256": sha256_file(paper_metadata.resolve()),
            "paper_source_pdf_sha256": paper_info["source_pdf_sha256"],
            "paper_digitization_method": paper_info.get("method", ""),
        },
        "cost": science["cost"],
        "profiles": profiles,
    }


def write_csv(path: Path, product: Mapping[str, Any]) -> None:
    columns = (
        "species",
        "source",
        "paper_figure",
        "paper_variant",
        "radius_km",
        "mean_density_cm-3",
        "standard_error_cm-3",
        "confidence_low_cm-3",
        "confidence_high_cm-3",
        "seed_count",
        "stage9_outcome",
    )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for species, profile in product["profiles"].items():
            for method in science_gate.METHODS:
                for row in profile[method]:
                    writer.writerow(
                        {
                            "species": species,
                            "source": method,
                            "paper_figure": "",
                            "paper_variant": "",
                            **row,
                            "seed_count": len(product["seeds"]),
                            "stage9_outcome": product["stage9_outcome"],
                        }
                    )
            for curve in profile["RUBIN2011"]:
                for row in curve["samples"]:
                    writer.writerow(
                        {
                            "species": species,
                            "source": "RUBIN2011",
                            "paper_figure": curve["figure"],
                            "paper_variant": curve["variant"],
                            "radius_km": row["radius_km"],
                            "mean_density_cm-3": row["density_cm-3"],
                            "standard_error_cm-3": "",
                            "confidence_low_cm-3": "",
                            "confidence_high_cm-3": "",
                            "seed_count": "",
                            "stage9_outcome": product["stage9_outcome"],
                        }
                    )


def positive(values: Iterable[float]) -> List[float]:
    """Represent nonpositive log-axis values as gaps rather than clipping."""

    return [value if value > 0.0 else math.nan for value in values]


def write_figure(base: Path, product: Mapping[str, Any], dpi: int) -> None:
    species_list = list(product["profiles"])
    columns = 2 if len(species_list) > 1 else 1
    rows = (len(species_list) + columns - 1) // columns
    figure, axes = plt.subplots(rows, columns, figsize=(7.2 * columns, 4.8 * rows), squeeze=False)
    colors = {"ANALOG": "#1f77b4", "FORCED": "#d62728"}
    for index, species in enumerate(species_list):
        axis = axes[index // columns][index % columns]
        profile = product["profiles"][species]
        for method in science_gate.METHODS:
            rows_for_method = profile[method]
            radius = [row["radius_km"] for row in rows_for_method]
            mean = positive(row["mean_density_cm-3"] for row in rows_for_method)
            low = positive(row["confidence_low_cm-3"] for row in rows_for_method)
            high = positive(row["confidence_high_cm-3"] for row in rows_for_method)
            axis.plot(radius, mean, color=colors[method], linewidth=1.8, label=f"{method} mean")
            axis.plot(radius, low, color=colors[method], linewidth=0.8, linestyle="--", label=f"{method} CI")
            axis.plot(radius, high, color=colors[method], linewidth=0.8, linestyle="--")
        for paper_index, curve in enumerate(profile["RUBIN2011"]):
            label = f"Rubin 2011 Fig. {curve['figure']} ({curve['variant']})"
            axis.plot(
                [row["radius_km"] for row in curve["samples"]],
                [row["density_cm-3"] for row in curve["samples"]],
                color="black",
                linewidth=1.3,
                linestyle=("-" if paper_index == 0 else ":"),
                label=label,
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Cometocentric radius [km]")
        axis.set_ylabel("Number density [cm$^{-3}$]")
        axis.set_title(species)
        axis.grid(True, which="both", linewidth=0.35, color="0.85")
        # A fully opaque legend is important for EPS: PostScript has no alpha
        # channel, and a backend warning would otherwise obscure a real plot
        # generation failure in an automated campaign.
        axis.legend(fontsize="small", framealpha=1.0)
    for index in range(len(species_list), rows * columns):
        axes[index // columns][index % columns].set_visible(False)
    figure.suptitle(
        f"Stage 9 {product['stage9_outcome']}: ANALOG/FORCED ensemble vs. Rubin et al. (2011)\n"
        f"seeds={','.join(str(seed) for seed in product['seeds'])}; dashed curves are confidence bounds"
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    figure.savefig(base.with_suffix(".png"), dpi=dpi)
    figure.savefig(base.with_suffix(".eps"))
    figure.savefig(base.with_suffix(".pdf"))
    plt.close(figure)


def publish(
    manifest: Path,
    output_directory: Path,
    paper_csv: Path = DEFAULT_PAPER_DATA,
    paper_metadata: Path = DEFAULT_PAPER_METADATA,
    dpi: int = 180,
) -> Dict[str, Any]:
    if dpi < 72:
        raise ComparisonError("DPI must be at least 72")
    output_directory.mkdir(parents=True, exist_ok=True)
    base = output_directory / "stage9_rubin2011_ensemble_comparison"
    destinations = [base.with_suffix(suffix) for suffix in (".json", ".csv", ".png", ".eps", ".pdf")]
    if any(path.exists() for path in destinations):
        existing = ", ".join(str(path) for path in destinations if path.exists())
        raise ComparisonError(f"refusing to overwrite comparison evidence: {existing}")
    product = build_product(manifest, paper_csv, paper_metadata)
    # Render every member in an isolated directory first.  A missing plotting
    # backend or full filesystem therefore cannot leave a JSON/CSV pair that
    # looks like a successfully published five-format comparison set.
    with tempfile.TemporaryDirectory(prefix=".stage9-comparison-", dir=output_directory) as temporary:
        temporary_base = Path(temporary) / base.name
        temporary_base.with_suffix(".json").write_text(
            json.dumps(product, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        write_csv(temporary_base.with_suffix(".csv"), product)
        write_figure(temporary_base, product, dpi)
        for suffix in (".json", ".csv", ".png", ".eps", ".pdf"):
            os.replace(temporary_base.with_suffix(suffix), base.with_suffix(suffix))
    return product


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Publish a passed, checksum-authenticated Stage-9 ANALOG/FORCED "
            "ensemble beside digitized Rubin et al. (2011) neutral profiles."
        )
    )
    parser.add_argument("manifest", type=Path, help="completed Stage-9 RUN_MANIFEST.json")
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--paper-data", type=Path, default=DEFAULT_PAPER_DATA)
    parser.add_argument("--paper-metadata", type=Path, default=DEFAULT_PAPER_METADATA)
    parser.add_argument("--dpi", type=int, default=180)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        product = publish(
            arguments.manifest,
            arguments.output_directory,
            arguments.paper_data,
            arguments.paper_metadata,
            arguments.dpi,
        )
    except ComparisonError as exc:
        print(f"generate_stage9_rubin_comparison.py: {exc}", file=sys.stderr)
        return 2
    print(
        f"PASS Stage-9 Rubin ensemble comparison species={len(product['profiles'])} "
        f"output={arguments.output_directory.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

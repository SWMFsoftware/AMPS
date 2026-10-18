#!/usr/bin/env python3
"""Render the diagnostic CRaTER foreground-only comparison.

The figure is intentionally watermarked with the benchmark scope.  It is a
review aid, not a pixel-tested acceptance artifact, and it must not be relabeled
as a full CRaTER instrument validation because the benchmark excludes the
housing, spacecraft, side-entry background, and lunar albedo.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=pathlib.Path)
    parser.add_argument("output_prefix", type=pathlib.Path)
    args = parser.parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    require(
        payload.get("schema_version") == "shieldSim.crater_foreground_axial_result.v1",
        "unexpected CRaTER result schema",
    )
    require(payload.get("benchmark_scope") == "foreground_only_axial", "scope label missing")
    bins = payload["bins"]
    x = [float(row["center_keV_per_um"]) for row in bins]
    observed = [float(row["observed_per_cm2_sr_day_keV_per_um"]) for row in bins]
    ideal = [float(row["ideal_foreground_per_cm2_sr_day_keV_per_um"]) for row in bins]
    cut = [float(row["foreground_after_cuts_per_cm2_sr_day_keV_per_um"]) for row in bins]
    ideal_error = [float(row["ideal_mc_standard_error"]) for row in bins]
    cut_error = [float(row["cuts_mc_standard_error"]) for row in bins]

    fig, (spectrum, ratio) = plt.subplots(
        2, 1, figsize=(10.5, 8.4), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.0]},
    )
    spectrum.loglog(x, observed, color="#111111", linewidth=1.5,
                    label="CRaTER v1.3 observation (2009259–2010065)")
    positive_ideal = [value if value > 0.0 else math.nan for value in ideal]
    positive_cut = [value if value > 0.0 else math.nan for value in cut]
    spectrum.loglog(x, positive_ideal, color="#d95f02", linewidth=1.2,
                    label="Geant4 ideal axial foreground")
    spectrum.loglog(x, positive_cut, color="#1b9e77", linewidth=1.2,
                    label="Geant4 axial foreground after published cuts")
    # Error envelopes are shown only where their lower limits remain positive
    # on the logarithmic axis; omitted bands are zero-statistics bins, not zero
    # physical uncertainty.
    ideal_low = [max(v-e, v*1.0e-6) if v > 0 else math.nan
                 for v, e in zip(ideal, ideal_error)]
    ideal_high = [v+e if v > 0 else math.nan for v, e in zip(ideal, ideal_error)]
    cut_low = [max(v-e, v*1.0e-6) if v > 0 else math.nan
               for v, e in zip(cut, cut_error)]
    cut_high = [v+e if v > 0 else math.nan for v, e in zip(cut, cut_error)]
    spectrum.fill_between(x, ideal_low, ideal_high, color="#d95f02", alpha=0.16)
    spectrum.fill_between(x, cut_low, cut_high, color="#1b9e77", alpha=0.16)
    spectrum.axvline(20.0, color="#666666", linestyle=":", linewidth=1.0)
    spectrum.text(20.0, spectrum.get_ylim()[1], " D6 → D5 splice",
                  rotation=90, va="top", ha="left", fontsize=8, color="#555555")
    spectrum.set_ylabel("Differential normalized energy-deposit spectrum\n"
                        r"particles (cm$^2$ sr day keV/$\mu$m)$^{-1}$")
    spectrum.grid(which="both", alpha=0.22)
    spectrum.legend(fontsize=8, loc="best")

    ratio_ideal = [a/b if a > 0.0 and b > 0.0 else math.nan
                   for a, b in zip(ideal, observed)]
    ratio_cut = [a/b if a > 0.0 and b > 0.0 else math.nan
                 for a, b in zip(cut, observed)]
    ratio.semilogx(x, ratio_ideal, color="#d95f02", linewidth=1.0,
                   label="ideal / observed")
    ratio.semilogx(x, ratio_cut, color="#1b9e77", linewidth=1.0,
                   label="after cuts / observed")
    ratio.axhline(1.0, color="#222222", linewidth=0.9)
    ratio.axvline(20.0, color="#666666", linestyle=":", linewidth=1.0)
    ratio.set_ylim(0.0, 3.0)
    ratio.set_xlabel(r"D5/D6 normalized energy deposit (keV/$\mu$m)")
    ratio.set_ylabel("model / observed")
    ratio.grid(which="both", alpha=0.22)
    ratio.legend(fontsize=8, loc="best")

    fig.suptitle(
        "FOREGROUND-ONLY AXIAL CRaTER BENCHMARK\n"
        "Published axial stack + BON2020; not a full-instrument validation",
        fontsize=13, fontweight="bold",
    )
    fig.text(
        0.01, 0.008,
        "Diagnostic comparison only. Excludes housing, boards, dead layers, thermal blanket, "
        "spacecraft, lunar albedo, and all out-of-aperture/side-entry background. Observation "
        "covariance is unavailable; no pass/fail goodness-of-fit threshold is applied.",
        fontsize=7.5, color="#333333",
    )
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 0.91))
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        target = args.output_prefix.with_suffix(f".{extension}")
        fig.savefig(target, dpi=180 if extension == "png" else None,
                    bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

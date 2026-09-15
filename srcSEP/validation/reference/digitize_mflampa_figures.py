#!/usr/bin/env python3
"""Reproduce the coordinate calibration used by XM02/XM03 references.

This utility intentionally separates rendering from test execution. Published
figures are immutable inputs, while the small reviewed CSVs are committed test
references. Regeneration is explicit and writes to a caller-selected directory
so a new trace can be diffed and reviewed before replacing any baseline.

The trace coordinates below are human-reviewed curve centers at 400 dpi. This
is more auditable than a color-only automatic trace because dashed curves,
legends, and observational overlays share pixels. The script verifies source
PDF checksums, can render the cited pages with ``pdftoppm``, and applies the
documented linear/logarithmic axis transforms without calling srcSEP.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, subprocess
from pathlib import Path

SOURCES = {
    "zhao": "945361cc4c27e6481eac042eaf9a0e3f6b097ed277e1b7fab793ab9305591bbe",
    "liu": "a5a613e9ad3b127c5c412366b0c4a2029339f9ac068fd9508325ab682a6dc357",
}

# Figure 7 samples: elapsed hour -> centerline pixel y. Axis calibration is
# log10(pfu)=2-(y-497)/145. Samples intersecting the in-panel legend are
# deliberately traced on the physical curve, not selected by color alone.
ZHAO_FIG7 = {
    "mfp_0.05au_integral_gt10mev": [587.5,539,511.5,498,490,486.5,486,487,491.5,495,502.5,508,516,522.5,530.5,538,546,553.5,561.5,568],
    "mfp_0.3au_integral_gt10mev": [465,446.5,460.5,479,498,519.5,536.5,551.5,567,581,596,608.5,622,633,645,655,667,674.5,683.5,690],
    "mfp_1.0au_integral_gt10mev": [448.5,498.5,540,578,613,641,667,692,708.5,722.5,732,739.5,745,750.5,753.5,758,760.5,764,767,769.5],
}

# Figure 15 green-curve samples after applying each panel's 0--44 h and
# 10^-2--10^3 intensity calibration. Figure 14 values use 5-Rs sampling and
# its 10^-3--1 au log axis. Keeping the calibrated trace here lets reviewers
# diff deterministic CSV regeneration even when Poppler anti-aliasing differs.
LIU_LOW_ENERGY = [0.0244,1.5159,14.7254,36.5616,48.8737,59.3074,61.6476,
    60.4662,57.0562,51.7947,47.9372,44.3669,40.2756,37.2759,34.4997,
    31.9302,30.1295,28.4303,26.8270,25.3140,24.8289]
LIU_HIGH_ENERGY = [1.4443,2.4829,2.1684,1.8218,1.5605,1.3111,1.0805,
    0.9078,0.7928,0.6791,0.5706,0.4887,0.4352,0.3875,0.3450,0.3013,
    0.2735,0.2531,0.2343,0.2168,0.2046]
LIU_MFP = [0.001648,0.001770,0.002405,0.003880,0.013240,0.002529,
    0.003200,0.003486,0.003798,0.004108,0.004412,0.004738,0.005053,
    0.005465,0.005828,0.006083,0.006440,0.006818,0.007168,0.007535,0.007864]

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            value.update(block)
    return value.hexdigest()

def check(path: Path, expected: str) -> None:
    observed = digest(path)
    if observed != expected:
        raise RuntimeError(f"unexpected publication checksum for {path}: {observed}")

def render(pdf: Path, page: int, destination: Path) -> None:
    """Render one 1-based PDF page to PNG using Poppler at exactly 400 dpi."""
    prefix = destination.with_suffix("")
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "400",
                    "-png", "-singlefile", str(pdf), str(prefix)], check=True)

def write_xm02(path: Path) -> None:
    temporary = path.with_name(path.name+".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n"); writer.writerow(["elapsed_hours", "series", "intensity"])
        for series, pixels in ZHAO_FIG7.items():
            for hour, pixel_y in zip(range(4, 44, 2), pixels):
                writer.writerow([hour, series, f"{10**(2-(pixel_y-497)/145):.8f}"])
    temporary.replace(path)

def write_xm03(path: Path) -> None:
    """Write the reviewed Liu Figure 14/15 trace and reported fit."""
    temporary = path.with_name(path.name+".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["observable", "coordinate", "value"])
        for name, values in (("earth_2.0_2.5mev_intensity_lambda0.3au", LIU_LOW_ENERGY),
                             ("earth_20.0_25.0mev_intensity_lambda0.3au", LIU_HIGH_ENERGY)):
            for hour, value in zip(range(4, 46, 2), values):
                writer.writerow([name, hour, f"{value:.4f}"])
        for radius, value in zip(range(5, 110, 5), LIU_MFP):
            writer.writerow(["earth_mean_free_path_au", radius, f"{value:.6f}"])
        # -1.78 is the colored lambda0=0.3 au Earth fit printed in Figure
        # 15(c), not a slope refitted from pixel samples.
        writer.writerow(["earth_fluence_spectral_index", 0, "-1.78"])
    temporary.replace(path)

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zhao-pdf", required=True, type=Path)
    parser.add_argument("--liu-pdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--render-pages", action="store_true")
    args = parser.parse_args()
    check(args.zhao_pdf, SOURCES["zhao"]); check(args.liu_pdf, SOURCES["liu"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_xm02(args.output_dir/"mflampa_2013_apr11_mfp_sensitivity.csv")
    write_xm03(args.output_dir/"mflampa_2013_apr11_event.csv")
    if args.render_pages:
        render(args.zhao_pdf, 30, args.output_dir/"zhao_figure7.png")
        render(args.liu_pdf, 23, args.output_dir/"liu_figure14.png")
        render(args.liu_pdf, 24, args.output_dir/"liu_figure15.png")
    (args.output_dir/"digitization_run.json").write_text(json.dumps({
        "schema": "srcsep-digitization-run-v1", "source_sha256": SOURCES,
        "render_dpi": 400, "note": "Outputs reproduce the reviewed XM02 and XM03 trace tables; compare before baseline replacement."}, indent=2)+"\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

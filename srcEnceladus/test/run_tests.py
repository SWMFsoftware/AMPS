#!/usr/bin/env python3
"""Global Enceladus forward-model test runner with reference comparisons."""

from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
TEST = ROOT / "test"
REFERENCE = TEST / "reference" / "forward_model_reference.json"


def run(command: list[str], *, cwd: Path = ROOT) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True)


def parse_receipt(text: str) -> dict[str, object]:
    result: dict[str, object] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            result[key] = float(value)
        except ValueError:
            result[key] = value
    return result


def write_input(path: Path, source_name: str) -> None:
    path.write_text(f"""#begin application: enceladus
schema_version=2
source_distribution_file={source_name}
simulated_species=H2O,CO2
injected_species=H2O,CO2
chemistry_species=H2O
selected_reactions=burger2007.photo.h2o_to_oh_h
reaction_products=discard
sample_length_steps=200
total_iterations=800
flyby_density_output=true
flyby_id=E3
flyby_trajectory_file=trajectory.csv
flyby_mesh_refinement=true
flyby_refinement_target_cell_size_m=20000
flyby_refinement_radius_m=150000
flyby_refinement_time_min_s=-120
flyby_refinement_time_max_s=100
random_seed=19
particle_splitting_merging=velocity_shift
min_particles_per_cell=10
max_particles_per_cell=20
split_position_shift_cell_diagonal_fraction=0.2
split_velocity_shift_thermal_speed_fraction=0.01
nonuniform_split_position_shift=true
solar_activity=quiet
plasma_background_mode=none
total_source_rate.H2O_s-1=1.0e28
total_source_rate.CO2_s-1=2.0e26
primary_macroparticles_per_step.H2O=100
primary_macroparticles_per_step.CO2=50
characteristic_speed.H2O_m_s-1=2000
characteristic_speed.CO2_m_s-1=1000
courant.H2O=0.5
courant.CO2=0.5
#end application
""", encoding="utf-8")


def write_transport_only_input(path: Path, source_name: str) -> None:
    path.write_text(f"""#begin application: enceladus
schema_version=2
source_distribution_file={source_name}
simulated_species=H2O,CO2
injected_species=H2O,CO2
chemistry_species=none
selected_reactions=none
reaction_products=discard
sample_length_steps=200
total_iterations=800
random_seed=19
particle_splitting_merging=off
solar_activity=quiet
plasma_background_mode=cache_only
plasma_cache_directory=.
total_source_rate.H2O_s-1=1.0e28
total_source_rate.CO2_s-1=2.0e26
primary_macroparticles_per_step.H2O=100
primary_macroparticles_per_step.CO2=50
characteristic_speed.H2O_m_s-1=2000
characteristic_speed.CO2_m_s-1=1000
courant.H2O=0.5
courant.CO2=0.5
#end application
    """, encoding="utf-8")


def write_wrong_geometry_source(path: Path, source_text: str) -> None:
    """Swap a canonical jet/stripe ID while keeping both blocks well formed."""
    source_text = source_text.replace(
        "#begin source: I_baghdad", "#begin source: __swap__", 1)
    source_text = source_text.replace(
        "#begin source: alexandria", "#begin source: I_baghdad", 1)
    source_text = source_text.replace(
        "#begin source: __swap__", "#begin source: alexandria", 1)
    path.write_text(source_text, encoding="utf-8")


def svg_comparison(path: Path, rows: list[tuple[str, float, float]]) -> None:
    width, height = 980, 80 + 34 * len(rows)
    plot_x, plot_width = 410, 500
    pieces = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="20" y="30" font-family="sans-serif" font-size="18">'
        'Enceladus forward-model reference comparison</text>',
        '<text x="410" y="55" font-family="sans-serif" font-size="12">'
        'relative error (log10 scale; green threshold = 1e-12)</text>',
    ]
    for index, (name, actual, expected) in enumerate(rows):
        y = 82 + 34 * index
        relative = abs(actual - expected) / max(abs(expected), 1e-300)
        magnitude = -16 if relative == 0 else max(-16, min(0, math.log10(relative)))
        bar = (magnitude + 16) / 16 * plot_width
        color = "#2e7d32" if relative <= 1e-12 else "#c62828"
        pieces.extend([
            f'<text x="20" y="{y}" font-family="monospace" font-size="12">'
            f'{html.escape(name)}</text>',
            f'<line x1="{plot_x}" y1="{y-5}" x2="{plot_x+plot_width}" y2="{y-5}" '
            'stroke="#dddddd"/>',
            f'<rect x="{plot_x}" y="{y-15}" width="{max(1, bar):.3f}" height="14" '
            f'fill="{color}"/>',
            f'<text x="{plot_x+plot_width+10}" y="{y-3}" font-family="monospace" '
            f'font-size="11">{relative:.3e}</text>',
        ])
    pieces.append("</svg>")
    path.write_text("\n".join(pieces) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "test_output" / "srcEnceladus")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    outcomes = []
    values: dict[str, object] = {}

    standalone = run([sys.executable, "-m", "unittest",
                      "standalone.tests.test_export_forward_sources"])
    outcomes.append(("U00_exporter", standalone.returncode == 0,
                     standalone.stdout + standalone.stderr))

    with tempfile.TemporaryDirectory(prefix="enceladus_forward_") as temporary:
        temp = Path(temporary)
        source = temp / "sources.in"
        exported = run([
            sys.executable, str(ROOT / "standalone" / "export_forward_sources.py"),
            "--manifest", str(TEST / "data" / "shape_manifest.json"),
            "--output", str(source),
        ])
        outcomes.append(("U01_source_export", exported.returncode == 0,
                         exported.stdout + exported.stderr))
        binary = temp / "forward_model_tests"
        compiled = run([
            "g++", "-std=c++17", "-O2", "-I..", "-I../src",
            "ForwardModel.cpp", "test/forward_model_tests.cpp", "-o", str(binary),
        ])
        outcomes.append(("U02_core_compile", compiled.returncode == 0,
                         compiled.stdout + compiled.stderr))
        if compiled.returncode == 0 and exported.returncode == 0:
            # Swap one canonical jet ID with one canonical stripe ID while
            # leaving both blocks otherwise well formed. The parser must bind
            # the hybrid-8 IDs to their geometry classes, not merely count 15
            # valid source blocks.
            wrong_geometry = temp / "wrong_geometry.in"
            write_wrong_geometry_source(
                wrong_geometry, source.read_text(encoding="utf-8"))
            wrong_geometry_result = run(
                [str(binary), "source", str(wrong_geometry)])
            wrong_geometry_text = (wrong_geometry_result.stdout +
                                   wrong_geometry_result.stderr)
            outcomes.append((
                "U03b_canonical_source_geometry",
                wrong_geometry_result.returncode != 0 and
                "must use geometry=" in wrong_geometry_text,
                wrong_geometry_text))
            runtime_input = temp / "amps.in"
            write_input(runtime_input, source.name)
            transport_input = temp / "transport_only.in"
            write_transport_only_input(transport_input, source.name)
            for name, command in (
                ("U03_source_parser", [str(binary), "source", str(source)]),
                ("U04_runtime_parser", [str(binary), "input", str(runtime_input), str(source)]),
                ("U04b_transport_only_parser", [str(binary), "input", str(transport_input), str(source)]),
                ("U05_chemistry", [str(binary), "chemistry"]),
                ("U06_source_law", [str(binary), "source_law"]),
                ("U07_numerics", [str(binary), "numerics", str(source)]),
                ("U08_validation", [str(binary), "validation", str(source)]),
            ):
                result = run(command)
                outcomes.append((name, result.returncode == 0,
                                 result.stdout + result.stderr))
                if result.returncode == 0:
                    values.update(parse_receipt(result.stdout))

    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    comparison = []
    for key, expected in reference.items():
        if not isinstance(expected, (int, float)):
            continue
        actual = values.get(key)
        passed = isinstance(actual, float) and math.isclose(
            actual, float(expected), rel_tol=1e-12, abs_tol=1e-300)
        outcomes.append((f"R_{key}", passed,
                         f"actual={actual!r}, expected={expected!r}\n"))
        if isinstance(actual, float):
            comparison.append((key, actual, float(expected)))
    svg_comparison(output / "reference_comparison.svg", comparison)

    receipt = {
        "runner": str(Path(__file__).resolve()),
        "elapsed_s": time.time() - started,
        "reference": str(REFERENCE.resolve()),
        "results": [
            {"test": name, "status": "PASS" if passed else "FAIL",
             "detail": detail}
            for name, passed, detail in outcomes
        ],
        "calculated_values": values,
        "visualization": str((output / "reference_comparison.svg").resolve()),
    }
    (output / "test_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for name, passed, detail in outcomes:
        print(f"{'PASS' if passed else 'FAIL'} {name}")
        if not passed and detail:
            print(detail.rstrip())
    failed = sum(not item[1] for item in outcomes)
    print(f"SUMMARY pass={len(outcomes)-failed} fail={failed} "
          f"receipt={output/'test_receipt.json'}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

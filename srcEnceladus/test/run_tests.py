#!/usr/bin/env python3
"""Global Enceladus forward-model test runner with reference comparisons."""

from __future__ import annotations

import argparse
import html
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
TEST = ROOT / "test"
REFERENCE = TEST / "reference" / "forward_model_reference.json"
STANDALONE_TESTS = ROOT / "standalone" / "tests"
# Every test of standalone/tests is reported as its own row S_<module>.<Class>.<method> (standalone_test_receipt.py).
# These four exporter/publisher modules run together, in the caller's environment, and their tests also make up the
# registered aggregate case U00_exporter (PASS when none of them fails).  test_prepare_forward_cases (the
# measurement--basis bundle publisher) belongs to U00 as documented in test/README.md.  Every other module runs in
# its own interpreter with the validation-data and SPICE roots below.
EXPORTER_MODULES = ("test_export_forward_sources", "test_fit_source_export", "test_prepare_forward_cases",
                    "test_export_forward_observations")
STANDALONE_RECEIPT = TEST / "standalone_test_receipt.py"
# Inputs of the standalone tests that are machine dependent, passed on in these environment variables
DATA_ROOT_ENV = "ENCELADUS_VALIDATION_DATA_ROOT"          # Cassini validation data (standalone/validation_paths.py)
SPICE_TOOLKIT_ENV = "ENCELADUS_SPICE_TOOLKIT_ROOT"        # CSPICE (standalone/tests/test_spice_ground_tracks.py)
SPICE_KERNELS_ENV = "ENCELADUS_SPICE_KERNEL_ROOT"         # its NAIF/ and CASSINI/ kernels
SKIP = None                                               # outcome of a case whose tests were all skipped


def run(command: list[str], *, cwd: Path = ROOT, env: dict | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True, env=env)


def standalone_tests(modules: list[str], env: dict | None, record: Path) -> list[tuple]:
    """Run standalone/tests modules in one fresh interpreter; one (name, passed, detail, duration_s) per test.

    passed is True, False or SKIP.  The per-test JSON written by standalone_test_receipt.py is kept at ``record``.
    If the interpreter ends without writing it (a crash, not a test failure), every module becomes one failing row
    S_<module> carrying the interpreter's output, so no test disappears silently.
    """
    record.parent.mkdir(parents=True, exist_ok=True)
    if record.exists():
        record.unlink()
    result = run([sys.executable, "-B", str(STANDALONE_RECEIPT), "--output", str(record)] + list(modules), env=env)
    transcript = result.stdout + result.stderr
    try:
        rows = json.loads(record.read_text(encoding="utf-8"))["tests"]
    except (OSError, ValueError, KeyError, TypeError) as error:
        return [("S_" + module[len("test_"):], False,
                 "ERROR: standalone test recorder wrote no outcomes ({}; exit {})\n{}".format(
                     error, result.returncode, transcript), None) for module in modules]
    if result.returncode != (1 if any(row["status"] == "FAIL" for row in rows) else 0):
        rows.append({"name": "S_recorder_" + modules[0][len("test_"):], "status": "FAIL", "duration_s": None,
                     "detail": "ERROR: recorder exit {} disagrees with its outcomes\n{}".format(
                         result.returncode, transcript)})
    status = {"PASS": True, "FAIL": False, "SKIP": SKIP}
    return [(row["name"], status[row["status"]], row["detail"] + "\n", row.get("duration_s")) for row in rows]


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
volume_field_output=false
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
observation_data_root=.
#begin observation: inms_test
instrument=INMS
event_id=E3
data_file=inms.csv
time_column=time_ca_s
x_column=x_m
y_column=y_m
z_column=z_m
comparison_species=H2O
value_column=density_m3
uncertainty_column=sigma_m3
#end observation
#begin observation: inms_context_test
instrument=INMS
event_id=E21
data_file=input:context.csv
time_column=time_ca_s
x_column=x_m
y_column=y_m
z_column=z_m
comparison_species=none
context_observable=INMS_OSNB_H2O_counts_per_integration
context_value_column=counts
context_uncertainty_column=sigma_counts
#end observation
#begin observation: uvis_test
instrument=UVIS
event_id=uvis_test_event
data_file=uvis.csv
ray_geometry=finite_observer_to_sun
observer_geometry_file=observer.csv
published_h2o_column_m-2=1.4e20
published_h2o_column_sigma_m-2=1.4e19
#end observation
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
    parser.add_argument(
        "--validation-data-root", type=Path, default=os.environ.get(DATA_ROOT_ENV) or None,
        help=("Cassini validation-data root for the standalone tests (default: ${}); without it the tests "
              "that need it are skipped".format(DATA_ROOT_ENV)))
    parser.add_argument(
        "--spice-path", type=Path, default=os.environ.get(SPICE_TOOLKIT_ENV) or None,
        help="CSPICE toolkit for the SPICE ground-track test (default: ${})".format(SPICE_TOOLKIT_ENV))
    parser.add_argument(
        "--spice-kernels", type=Path, default=os.environ.get(SPICE_KERNELS_ENV) or None,
        help=("SPICE kernel tree (NAIF/, CASSINI/) for the SPICE ground-track test (default: ${})"
              .format(SPICE_KERNELS_ENV)))
    parser.add_argument(
        "--cxx", nargs="+", default=["g++"], metavar="ARG",
        help=("C++ compiler argv used for the dependency-light forward-model "
              "binary (default: g++). The global test adapter passes the "
              "machine-profile compiler vector here without invoking a shell."),
    )
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    outcomes = []
    durations: dict[str, float] = {}
    values: dict[str, object] = {}
    records = output / "standalone_tests"

    def add_standalone(rows: list[tuple]) -> None:
        for name, passed, detail, duration in rows:
            outcomes.append((name, passed, detail))
            if isinstance(duration, (int, float)):
                durations[name] = duration

    # U00_exporter: the exporter/publisher modules in the caller's environment, as before; the aggregate row
    # comes first (registered order), then one row per test.
    exporter = standalone_tests(list(EXPORTER_MODULES), None, records / "exporter.json")
    exporter_failed = [name for name, passed, _, _ in exporter if passed is False]
    exporter_skipped = [name for name, passed, _, _ in exporter if passed is SKIP]
    summary = "{} exporter tests: {} PASS, {} SKIP, {} FAIL (their S_ rows follow)\n".format(
        len(exporter), len(exporter) - len(exporter_failed) - len(exporter_skipped), len(exporter_skipped),
        len(exporter_failed))
    outcomes.append(("U00_exporter", bool(exporter) and not exporter_failed,
                     summary + "".join("FAIL: {}\n".format(name) for name in exporter_failed)))
    add_standalone(exporter)

    # The other tests of the standalone plume model (standalone/README.md): each module in its own interpreter
    standalone_env = dict(os.environ)
    for variable, value in ((DATA_ROOT_ENV, args.validation_data_root), (SPICE_TOOLKIT_ENV, args.spice_path),
                            (SPICE_KERNELS_ENV, args.spice_kernels)):
        if value is not None:
            standalone_env[variable] = str(Path(value).expanduser().resolve())
    for path in sorted(STANDALONE_TESTS.glob("test_*.py")):
        if path.stem not in EXPORTER_MODULES:
            add_standalone(standalone_tests([path.stem], standalone_env, records / (path.stem + ".json")))

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
        # The compiler is configuration, not model physics.  Accepting an argv
        # vector makes the same test portable to compiler wrappers and cluster
        # modules while preserving exact flag ordering and avoiding shell
        # expansion.  The source/include contract and optimization level stay
        # owned by this application test.
        compiled = run(list(args.cxx) + [
            "-std=c++17", "-O2", "-I..", "-I../src",
            "ForwardModel.cpp", "MeshControls.cpp", "CaseSelection.cpp", "SourceCatalog.cpp", "test/forward_model_tests.cpp", "-o", str(binary),
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
            catalog = run([sys.executable, str(TEST/"run_source_export_roundtrip.py"),
                           "--executable", str(binary)])
            outcomes.append(("U03c_catalog_source_geometry", catalog.returncode == 0,
                             catalog.stdout+catalog.stderr))
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
                ("U09_case_runtime", [sys.executable,str(TEST/"run_case_runtime_tests.py"),"--executable",str(binary)]),
                ("U10_runtime_mesh", [str(binary),"runtime_mesh"]),
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
            dict({"test": name, "status": "SKIP" if passed is SKIP else "PASS" if passed else "FAIL",
                  "detail": detail}, **({"duration_s": durations[name]} if name in durations else {}))
            for name, passed, detail in outcomes
        ],
        "calculated_values": values,
        "visualization": str((output / "reference_comparison.svg").resolve()),
    }
    (output / "test_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for name, passed, detail in outcomes:
        print(f"{'SKIP' if passed is SKIP else 'PASS' if passed else 'FAIL'} {name}")
        if passed is SKIP:
            print(detail.rstrip().splitlines()[-1])
        elif not passed and detail:
            print(detail.rstrip())
    failed = sum(item[1] is False for item in outcomes)
    skipped = sum(item[1] is SKIP for item in outcomes)
    print(f"SUMMARY pass={len(outcomes)-failed-skipped} fail={failed} skip={skipped} "
          f"receipt={output/'test_receipt.json'}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

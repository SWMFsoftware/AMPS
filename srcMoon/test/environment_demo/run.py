#!/usr/bin/env python3
"""Export meeting-oriented production illumination and thermal-state products.

The runner is deliberately outside U01--U26/I01--I46.  It requires the local
U08/U09/U11 kernels, invokes a thin C++ adapter linked to libAMPS.a, and keeps
formal I15/I17/I20 status at SKIPPED/NOT VALIDATED.  Generated tables are model
output and are never substituted for D02 or D03 observations.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys


TEST_ROOT = Path(__file__).resolve().parents[1]
REPO = TEST_ROOT.parents[1]
sys.path.insert(0, str(TEST_ROOT / "common"))
from moon_testlib import compile_probe, run_one  # noqa: E402


NOTICE = "MEETING DEMONSTRATION - FORMAL I15/I17/I20 NOT RUN"
D01_ROOT = Path("/data/vtenishe/moon_validation_data/lola")


def sha256(path: Path) -> str:
    """Hash every executable, receipt, data dependency, and output artifact."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def replace_setting(text: str, key: str, value: str) -> str:
    """Replace one explicit setting; duplicates/missing keys are errors."""
    pattern = rf"(?m)^({re.escape(key)}\s*=\s*)([^!\n]+?)(\s*)$"
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} in Moon input")
    return re.sub(pattern, lambda match: match.group(1) + value + match.group(3),
                  text)


def setting(text: str, key: str) -> str:
    """Read one flat Moon assignment without reimplementing parser policy."""
    pattern = rf"(?m)^{re.escape(key)}\s*=\s*([^!\n]+?)\s*$"
    matches = re.findall(pattern, text)
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} in Moon input")
    return matches[0].strip().strip('"\'')


def prepare_input(template: Path, destination: Path, resolution_m: float,
                  cea: Path, tecplot: Path) -> str:
    """Create an immutable run receipt with output paths owned by this run."""
    text = template.read_text(encoding="utf-8")
    start_utc = setting(text, "simulation_start_utc")
    replacements = {
        "surface_mesh_resolution_m": f"{resolution_m:.17g}",
        "lola_image_file": str((D01_ROOT / "raw/LDEM_4.IMG").resolve()),
        "lola_label_file": str((D01_ROOT / "raw/LDEM_4.LBL").resolve()),
        "surface_cea_file": str(cea.resolve()),
        "surface_tecplot_file": str(tecplot.resolve()),
        # A read-only demonstration must not create restart/checkpoint state.
        "surface_inventory_checkpoint_file": "none",
    }
    for key, value in replacements.items():
        text = replace_setting(text, key, value)
    destination.write_text(text, encoding="utf-8")
    return start_utc


def utc_after(start: str, seconds: float) -> str:
    """Derive a requested UTC using civil-time arithmetic only for CLI staging."""
    normalized = start[:-1] + "+00:00" if start.endswith("Z") else start
    parsed = dt.datetime.fromisoformat(normalized)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    result = parsed.astimezone(dt.timezone.utc) + dt.timedelta(seconds=seconds)
    return result.isoformat(timespec="milliseconds").replace("+00:00", "Z")


def parse_summary(text: str) -> dict[str, object]:
    """Parse the adapter's deliberately flat key/value diagnostics strictly."""
    strings = {"thermal_initialization"}
    summary: dict[str, object] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in strings:
            summary[key] = value
        else:
            try:
                summary[key] = float(value)
            except ValueError as exception:
                raise ValueError(f"invalid adapter summary line: {line}") from exception
    required = {
        "surface_element_count", "thermal_profile_row_count",
        "terrain_occluded_count", "terrain_visible_count",
        "minimum_surface_temperature_k", "maximum_surface_temperature_k",
        "thermal_initialization", "snapshot_et_s", "advanced_duration_s",
        "spinup_completed_cycles", "spinup_maximum_node_change_k",
        "spinup_maximum_enthalpy_change_j_m2",
    }
    if set(summary) != required:
        raise ValueError("adapter summary differs: missing=" +
                         repr(sorted(required - set(summary))) + " extra=" +
                         repr(sorted(set(summary) - required)))
    return summary


def inspect_surface_csv(path: Path) -> dict[str, object]:
    """Independently check row count/ranges; this is integrity, not validation."""
    expected = [
        "surface_element", "zenith_index", "azimuth_index", "latitude_deg",
        "longitude_east_deg", "local_solar_time_hour", "x_body_fixed_m",
        "y_body_fixed_m", "z_body_fixed_m", "radius_m", "elevation_m",
        "normal_x_body_fixed", "normal_y_body_fixed", "normal_z_body_fixed",
        "solar_incidence_cosine", "sun_above_local_horizon",
        "terrain_sun_visible", "terrain_occluded",
        "solar_disc_visible_fraction", "earth_eclipse_fraction",
        "direct_beam_fraction", "absorbed_solar_flux_w_m2",
        "surface_temperature_k", "last_energy_residual_j_m2",
        "last_nonlinear_iterations",
    ]
    count = 0
    temperature: list[float] = []
    terrain = 0
    eclipse = 0
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected:
            raise ValueError("surface_environment.csv schema differs")
        for line_number, row in enumerate(reader, 2):
            try:
                numeric = {key: float(row[key]) for key in expected}
            except ValueError as exception:
                raise ValueError(f"non-numeric surface row {line_number}") from exception
            if not all(math.isfinite(value) for value in numeric.values()):
                raise ValueError(f"non-finite surface row {line_number}")
            if int(numeric["surface_element"]) != count:
                raise ValueError("surface element order is not contiguous")
            temperature.append(numeric["surface_temperature_k"])
            terrain += int(numeric["terrain_occluded"])
            eclipse += numeric["earth_eclipse_fraction"] > 0.0
            count += 1
    if count == 0 or min(temperature) <= 0.0:
        raise ValueError("surface environment table is empty or non-physical")
    return {"row_count": count, "minimum_temperature_k": min(temperature),
            "maximum_temperature_k": max(temperature),
            "terrain_occluded_count": terrain,
            "earth_eclipse_affected_count": eclipse}


def run(arguments: argparse.Namespace) -> dict[str, object]:
    """Run local gates, production adapter, integrity checks, and rendering."""
    output = arguments.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if not math.isfinite(arguments.surface_resolution_m) or not (
            10_000.0 <= arguments.surface_resolution_m <= 1_000_000.0):
        raise ValueError("surface resolution must lie within [10000,1000000] m")
    if not math.isfinite(arguments.advance_seconds) or arguments.advance_seconds < 0.0:
        raise ValueError("advance-seconds must be finite and non-negative")

    verification_root = output / "verification"
    verification: dict[str, dict[str, object]] = {}
    for test_id in ("U08", "U09", "U11"):
        result = run_one(test_id, arguments.build_dir.resolve(), verification_root)
        verification[test_id] = result
        if result.get("status") != "PASS":
            raise RuntimeError(
                f"environment demonstration requires {test_id} PASS; observed "
                f"{result.get('status')}: {result.get('reason', '')}")

    executable, compilation = compile_probe(
        arguments.build_dir.resolve(), output / "adapter")
    if executable is None:
        raise RuntimeError("cannot compile production adapter: " + repr(compilation))
    template = arguments.input.resolve()
    application = output / "application.in"
    cea = output / "lola_surface.cea"
    tecplot = output / "lola_surface.dat"
    surface_csv = output / "surface_environment.csv"
    profiles_csv = output / "thermal_profiles.csv"
    start_utc = prepare_input(template, application,
                              arguments.surface_resolution_m, cea, tecplot)
    snapshot_utc = arguments.snapshot_utc or utc_after(
        start_utc, arguments.advance_seconds)
    if arguments.snapshot_utc and arguments.advance_seconds != 0.0:
        raise ValueError("use either --snapshot-utc or --advance-seconds, not both")
    command = [str(executable), "surface-environment-demo", str(application),
               str(surface_csv), str(profiles_csv), snapshot_utc,
               arguments.initialization]
    completed = subprocess.run(command, cwd=REPO, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               check=False)
    (output / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (output / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"production environment adapter exited {completed.returncode}; "
            f"see {output / 'stderr.log'}")
    summary = parse_summary(completed.stdout)
    integrity = inspect_surface_csv(surface_csv)
    if int(summary["surface_element_count"]) != integrity["row_count"]:
        raise RuntimeError("adapter and CSV surface-element counts differ")

    visualizer_command = [sys.executable, "-B",
        str(REPO / "srcMoon/tools/visualize_environment.py"),
        str(surface_csv), str(profiles_csv), "--output-dir", str(output),
        "--prefix", "moon_environment", "--epoch-label", snapshot_utc]
    visualizer = subprocess.run(visualizer_command, cwd=REPO, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                check=False)
    (output / "visualizer_stdout.log").write_text(
        visualizer.stdout, encoding="utf-8")
    (output / "visualizer_stderr.log").write_text(
        visualizer.stderr, encoding="utf-8")
    if visualizer.returncode != 0:
        raise RuntimeError(
            f"environment visualizer exited {visualizer.returncode}; see "
            f"{output / 'visualizer_stderr.log'}")

    formal = {
        test_id: {"execution": "NOT_RUN", "status": "SKIPPED",
                  "validation_status": "NOT_VALIDATED", "reason": reason}
        for test_id, reason in {
            "I15": "No independent D02 comparison or same-campaign I13 gate was run.",
            "I17": "No time-step convergence, READY D03 comparison, or hold-out score was run.",
            "I20": "No independently qualified PSR/cold-trapping comparison was run.",
        }.items()
    }
    artifact_paths = {
        "application_input": application, "cea_surface": cea,
        "tecplot_surface": tecplot, "surface_environment_csv": surface_csv,
        "thermal_profiles_csv": profiles_csv,
        "illumination_png": output / "moon_environment_illumination_shadowing.png",
        "illumination_eps": output / "moon_environment_illumination_shadowing.eps",
        "thermal_png": output / "moon_environment_thermal_state.png",
        "thermal_eps": output / "moon_environment_thermal_state.eps",
        "visualization_receipt": output / "moon_environment_visualization.json",
        "stdout": output / "stdout.log", "stderr": output / "stderr.log",
        "visualizer_stdout": output / "visualizer_stdout.log",
        "visualizer_stderr": output / "visualizer_stderr.log",
    }
    for test_id in verification:
        artifact_paths[f"{test_id.lower()}_result"] = (
            verification_root / test_id / "result.json")
    result: dict[str, object] = {
        "schema": "srcMoon-environment-demonstration-v1",
        "execution_status": "PASS", "notice": NOTICE,
        "capability_label": "EXPERIMENTAL",
        "artifact_class": "MEETING_DEMONSTRATION",
        "snapshot_utc": snapshot_utc,
        "thermal_initialization": arguments.initialization,
        "initial_temperature_source": (
            "periodic production spin-up using the input cycle/tolerances"
            if arguments.initialization == "periodic-spinup" else
            "thermal_initial_temperature_k at simulation_start_utc, then production advance"),
        "adapter_summary": summary, "integrity_metrics": integrity,
        "local_prerequisites": {key: value["status"]
                                for key, value in verification.items()},
        "formal_vv": formal, "artifacts": {},
    }

    build_library = arguments.build_dir.resolve() / "libAMPS.a"
    compiler = subprocess.run(["mpicxx", "--version"], cwd=REPO, text=True,
                              capture_output=True, check=False)
    data_files = [D01_ROOT / "README.md", D01_ROOT / "provenance.json",
                  D01_ROOT / "qa_report.json", D01_ROOT / "raw/LDEM_4.IMG",
                  D01_ROOT / "raw/LDEM_4.LBL"]
    manifest = {
        "schema": "srcMoon-environment-demonstration-manifest-v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_revision": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
            capture_output=True, check=False).stdout.strip(),
        "source_dirty_paths": subprocess.run(
            ["git", "status", "--short"], cwd=REPO, text=True,
            capture_output=True, check=False).stdout.splitlines(),
        "command": command, "visualizer_command": visualizer_command,
        "snapshot_utc": snapshot_utc,
        "adapter": {"path": str(executable), "sha256": sha256(executable)},
        "build_library": {"path": str(build_library),
                          "sha256": sha256(build_library)},
        "compiler": {"command": ["mpicxx", "--version"],
                     "returncode": compiler.returncode,
                     "first_line": compiler.stdout.splitlines()[0]
                     if compiler.stdout.splitlines() else ""},
        "template": {"path": str(template), "sha256": sha256(template)},
        "roadmaps": {
            "source": {"path": str(REPO / "srcMoon/development_and_validation.md"),
                       "sha256": sha256(REPO / "srcMoon/development_and_validation.md")},
            "data": {"path": str(Path(
                "/data/vtenishe/moon_validation_data/development_and_validation.md")),
                     "sha256": sha256(Path(
                "/data/vtenishe/moon_validation_data/development_and_validation.md"))},
        },
        "external_data": [{"path": str(path), "sha256": sha256(path)}
                          for path in data_files],
        "formal_vv_status": formal,
    }
    manifest_path = output / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    artifact_paths["run_manifest"] = manifest_path
    result_artifacts = result["artifacts"]
    assert isinstance(result_artifacts, dict)
    for name, path in artifact_paths.items():
        if not path.is_file():
            raise RuntimeError(f"required demonstration artifact is missing: {path}")
        result_artifacts[name] = {"path": str(path), "sha256": sha256(path)}
    (output / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser = argparse.ArgumentParser(
        description="Export production LOLA illumination/shadow/thermal model data",
        epilog=("Default input-initial mode advances from the explicit input "
                "thermal_initial_temperature_k for --advance-seconds. It is an "
                "initial-condition transient. Use --initialization periodic-spinup "
                "only when the potentially long configured convergence run is intended."))
    parser.add_argument("--build-dir", type=Path, default=REPO / "build")
    parser.add_argument("--input", type=Path,
                        default=REPO / "srcMoon/examples/lola_surface.in")
    parser.add_argument("--surface-resolution-m", type=float, default=100_000.0)
    parser.add_argument("--advance-seconds", type=float, default=0.0,
                        help="advance from simulation_start_utc; default 0 preserves the explicit initial state")
    parser.add_argument("--snapshot-utc",
                        help="explicit UTC snapshot; mutually exclusive with nonzero --advance-seconds")
    parser.add_argument("--initialization",
                        choices=("input-initial", "periodic-spinup"),
                        default="input-initial")
    parser.add_argument("--output-dir", type=Path,
                        default=REPO / "test_output/srcMoon/environment-demo" /
                        timestamp)
    arguments = parser.parse_args()
    try:
        result = run(arguments)
    except (OSError, ValueError, RuntimeError, ImportError) as exception:
        print(f"ERROR: {exception}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

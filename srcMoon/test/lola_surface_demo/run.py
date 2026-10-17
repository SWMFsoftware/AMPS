#!/usr/bin/env python3
"""Produce a meeting-ready LOLA surface demonstration from production kernels.

This runner is intentionally outside the frozen U01--U26/I01--I46 registry.
It first requires U08 to pass, then calls a thin adapter linked against the
same srcMoon parser, LOLA decoder, triangulation, altitude-adjustment, and
surface writers used by production.  Its figures are demonstration artifacts;
they never satisfy the linked I14 mover/BVH/predecessor contract.
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
sys.path.insert(0, str(TEST_ROOT))
from moon_testlib import compile_probe, run_one  # noqa: E402
from integration_common import ldem4_elevation_m  # noqa: E402


NOTICE = "MEETING DEMONSTRATION - FORMAL I14 NOT RUN"
D01_ROOT = Path("/data/vtenishe/moon_validation_data/lola")


def sha256(path: Path) -> str:
    """Hash every source/data/artifact named by the demonstration manifest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def replace_setting(text: str, key: str, value: str) -> str:
    """Replace one explicit Moon assignment and reject ambiguous templates."""
    pattern = rf"(?m)^({re.escape(key)}\s*=\s*)([^!\n]+?)(\s*)$"
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} in LOLA template")
    return re.sub(pattern, lambda match: match.group(1) + value + match.group(3),
                  text)


def prepare_input(template: Path, destination: Path, resolution_m: float,
                  cea: Path, tecplot: Path) -> None:
    """Derive a preserved production receipt without editing the template."""
    text = template.read_text(encoding="utf-8")
    text = replace_setting(text, "surface_mesh_resolution_m",
                           f"{resolution_m:.17g}")
    text = replace_setting(text, "lola_image_file",
                           str((D01_ROOT / "raw/LDEM_4.IMG").resolve()))
    text = replace_setting(text, "lola_label_file",
                           str((D01_ROOT / "raw/LDEM_4.LBL").resolve()))
    text = replace_setting(text, "surface_cea_file", str(cea.resolve()))
    text = replace_setting(text, "surface_tecplot_file",
                           str(tecplot.resolve()))
    # The adapter stops after surface construction. Disabling checkpoint output
    # in the archived receipt makes that limited scope explicit and prevents a
    # reader from mistaking it for an interrupted particle/inventory run.
    text = replace_setting(text, "surface_inventory_checkpoint_file", "none")
    destination.write_text(text, encoding="utf-8")


def parse_summary(text: str) -> dict[str, float]:
    """Parse the deliberately simple production-adapter summary strictly."""
    result: dict[str, float] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            result[key] = float(value)
        except ValueError as exception:
            raise ValueError(f"non-numeric adapter summary value: {line}") from exception
    required = {
        "surface_resolution_m", "reference_radius_m", "subdivision_level",
        "vertex_count", "face_count", "maximum_edge_length_m",
        "minimum_radius_m", "maximum_radius_m", "degenerate_face_count",
        "inward_face_count",
    }
    if set(result) != required:
        raise ValueError(
            "adapter summary fields differ: missing=" +
            repr(sorted(required - set(result))) + " extra=" +
            repr(sorted(set(result) - required)))
    return result


def compare_vertices(vertex_path: Path, image_path: Path,
                     comparison_path: Path) -> tuple[list[dict[str, float]], dict[str, float]]:
    """Compare every production vertex with an independent native-D01 decode."""
    image = image_path.read_bytes()
    if len(image) != 720 * 1440 * 2:
        raise ValueError("LDEM_4 byte count violates the detached-label contract")
    rows: list[dict[str, float]] = []
    with vertex_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected = ["vertex_index", "x_m", "y_m", "z_m", "latitude_deg",
                    "longitude_east_deg", "radius_m", "model_elevation_m"]
        if reader.fieldnames != expected:
            raise ValueError("production vertex table has an unexpected schema")
        for line_number, row in enumerate(reader, 2):
            try:
                values = {key: float(row[key]) for key in expected}
            except ValueError as exception:
                raise ValueError(
                    f"non-numeric production vertex on line {line_number}") from exception
            if not all(math.isfinite(value) for value in values.values()):
                raise ValueError(f"non-finite production vertex on line {line_number}")
            reference = ldem4_elevation_m(
                image, values["latitude_deg"], values["longitude_east_deg"])
            rows.append({**values, "d01_elevation_m": reference,
                         "model_minus_d01_m":
                             values["model_elevation_m"] - reference})
    if not rows:
        raise ValueError("production vertex table is empty")
    fields = list(rows[0])
    with comparison_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{value:.17g}" for key, value in row.items()})
    errors = [row["model_minus_d01_m"] for row in rows]
    absolute = [abs(value) for value in errors]
    metrics = {
        "maximum_absolute_vertex_error_m": max(absolute),
        "rms_vertex_error_m": math.sqrt(
            sum(value * value for value in errors) / len(errors)),
        "mean_vertex_error_m": sum(errors) / len(errors),
    }
    return rows, metrics


def render(rows: list[dict[str, float]], metrics: dict[str, float],
           png: Path, eps: Path) -> None:
    """Render unchanged comparison rows in presentation and vector formats."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    longitude = [row["longitude_east_deg"] for row in rows]
    latitude = [row["latitude_deg"] for row in rows]
    model = [row["model_elevation_m"] for row in rows]
    reference = [row["d01_elevation_m"] for row in rows]
    residual = [row["model_minus_d01_m"] for row in rows]
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.5),
                                constrained_layout=True)
    surface = axes[0, 0].scatter(longitude, latitude, c=model, s=5,
                                 cmap="terrain", linewidths=0)
    axes[0, 0].set_title("Production triangulated surface elevation")
    axes[0, 0].set_xlabel("east longitude [deg]")
    axes[0, 0].set_ylabel("planetocentric latitude [deg]")
    figure.colorbar(surface, ax=axes[0, 0], label="elevation [m]")

    axes[0, 1].scatter(reference, model, s=5, color="tab:blue", linewidths=0)
    minimum = min(min(reference), min(model))
    maximum = max(max(reference), max(model))
    axes[0, 1].plot((minimum, maximum), (minimum, maximum), color="black",
                    linewidth=1.0, label="identity")
    axes[0, 1].set_title("Every mesh vertex versus native LDEM_4")
    axes[0, 1].set_xlabel("independently decoded D01 elevation [m]")
    axes[0, 1].set_ylabel("production surface elevation [m]")
    axes[0, 1].legend(framealpha=1.0)

    error_map = axes[1, 0].scatter(longitude, latitude, c=residual, s=5,
                                   cmap="coolwarm", linewidths=0)
    axes[1, 0].set_title("Production minus independently decoded D01")
    axes[1, 0].set_xlabel("east longitude [deg]")
    axes[1, 0].set_ylabel("planetocentric latitude [deg]")
    figure.colorbar(error_map, ax=axes[1, 0], label="residual [m]")

    axes[1, 1].hist(residual, bins=50, color="tab:blue", edgecolor="black")
    axes[1, 1].set_title(
        "Vertex residual distribution\n"
        f"max |error|={metrics['maximum_absolute_vertex_error_m']:.3e} m, "
        f"RMS={metrics['rms_vertex_error_m']:.3e} m")
    axes[1, 1].set_xlabel("model minus D01 [m]")
    axes[1, 1].set_ylabel("vertex count")
    for axis in axes.flat:
        axis.grid(True, linewidth=0.35, color="#b0b0b0")
    figure.suptitle(
        "LOLA LDEM_4 production-surface demonstration\n" + NOTICE,
        fontsize=13)
    figure.savefig(png, dpi=180, transparent=False)
    figure.savefig(eps, format="eps", transparent=False)
    plt.close(figure)


def run(arguments: argparse.Namespace) -> dict[str, object]:
    """Execute verified production kernels and preserve a separate demo receipt."""
    output: Path = arguments.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if not math.isfinite(arguments.surface_resolution_m) or not (
            10_000.0 <= arguments.surface_resolution_m <= 1_000_000.0):
        raise ValueError("surface resolution must be within [10000,1000000] m")

    # U08 is the local correctness gate for all kernels this demonstration
    # calls. Its PASS remains U08 only and is nested as prerequisite evidence;
    # it is never relabeled as an I14 PASS.
    verification_root = output / "verification"
    u08 = run_one("U08", arguments.build_dir.resolve(), verification_root)
    if u08.get("status") != "PASS":
        raise RuntimeError(
            "LOLA demonstration requires U08 PASS; observed " +
            str(u08.get("status")) + ": " + str(u08.get("reason", "")))

    executable, compilation = compile_probe(
        arguments.build_dir.resolve(), output / "adapter")
    if executable is None:
        raise RuntimeError("cannot compile production adapter: " + repr(compilation))
    template = (REPO / "srcMoon/examples/lola_surface.in").resolve()
    application = output / "application.in"
    cea = output / "lola_surface.cea"
    tecplot = output / "lola_surface.dat"
    vertices = output / "surface_vertices.csv"
    comparison = output / "lola_dem_comparison.csv"
    png = output / "lola_surface_demonstration.png"
    eps = output / "lola_surface_demonstration.eps"
    prepare_input(template, application, arguments.surface_resolution_m,
                  cea, tecplot)
    command = [str(executable), "lola-surface-demo", str(application),
               str(vertices)]
    completed = subprocess.run(command, cwd=REPO, text=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=False)
    (output / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (output / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"production surface adapter exited {completed.returncode}; "
            f"see {output / 'stderr.log'}")
    summary = parse_summary(completed.stdout)
    if int(summary["degenerate_face_count"]) != 0 or int(
            summary["inward_face_count"]) != 0:
        raise RuntimeError("generated surface has invalid face topology")
    if summary["maximum_edge_length_m"] > arguments.surface_resolution_m:
        raise RuntimeError("generated surface violates requested edge resolution")
    rows, metrics = compare_vertices(
        vertices, D01_ROOT / "raw/LDEM_4.IMG", comparison)
    render(rows, metrics, png, eps)

    data_files = [D01_ROOT / "README.md", D01_ROOT / "provenance.json",
                  D01_ROOT / "qa_report.json", D01_ROOT / "raw/LDEM_4.IMG",
                  D01_ROOT / "raw/LDEM_4.LBL"]
    result: dict[str, object] = {
        "schema": "srcMoon-lola-surface-demonstration-v1",
        "execution_status": "PASS",
        "notice": NOTICE,
        # Use the roadmap's frozen capability vocabulary.  The independent
        # artifact class identifies presentation intent without inventing a
        # fifth capability state that could be mistaken for VERIFIED.
        "capability_label": "EXPERIMENTAL",
        "artifact_class": "MEETING_DEMONSTRATION",
        "formal_vv": {
            "test_id": "I14", "execution": "NOT_RUN",
            "status": "SKIPPED", "validation_status": "NOT_VALIDATED",
            "reason": (
                "This surface-only adapter does not execute the same-campaign "
                "I11/I13 predecessors, production BVH radial queries, or "
                "second-order particle/surface collision fixture required by I14."),
        },
        "u08_prerequisite_status": u08["status"],
        "surface": summary,
        "comparison_metrics": {
            **metrics,
            "vertex_count": len(rows),
            "within_i14_vertex_threshold_50_m":
                metrics["maximum_absolute_vertex_error_m"] <= 50.0,
        },
        "artifacts": {},
    }
    artifacts = result["artifacts"]
    assert isinstance(artifacts, dict)
    for name, path in {
        "application_input": application, "cea_surface": cea,
        "tecplot_surface": tecplot, "production_vertices": vertices,
        "comparison_csv": comparison, "comparison_png": png,
        "comparison_eps": eps, "stdout": output / "stdout.log",
        "stderr": output / "stderr.log",
        "u08_result": verification_root / "U08/result.json",
    }.items():
        artifacts[name] = {"path": str(path), "sha256": sha256(path)}

    # A commit hash alone cannot reproduce a demonstration made from a dirty
    # work tree.  Preserve the exact dirty-path inventory, linked-library
    # digest, and compiler identity without treating them as formal V&V gates.
    dirty_paths = subprocess.run(
        ["git", "status", "--short"], cwd=REPO, text=True,
        capture_output=True, check=False).stdout.splitlines()
    compiler = subprocess.run(
        ["mpicxx", "--version"], cwd=REPO, text=True,
        capture_output=True, check=False)
    build_library = arguments.build_dir.resolve() / "libAMPS.a"
    if not build_library.is_file():
        raise RuntimeError(f"linked AMPS library is missing: {build_library}")
    manifest = {
        "schema": "srcMoon-lola-surface-demonstration-manifest-v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_revision": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
            capture_output=True, check=False).stdout.strip(),
        "source_dirty_paths": dirty_paths,
        "command": command,
        "adapter": {"path": str(executable), "sha256": sha256(executable)},
        "build_library": {
            "path": str(build_library), "sha256": sha256(build_library)},
        "compiler": {
            "command": ["mpicxx", "--version"],
            "returncode": compiler.returncode,
            "first_line": compiler.stdout.splitlines()[0]
            if compiler.stdout.splitlines() else "",
        },
        "template": {"path": str(template), "sha256": sha256(template)},
        "roadmaps": {
            "source": {"path": str(REPO / "srcMoon/development_and_validation.md"),
                       "sha256": sha256(REPO / "srcMoon/development_and_validation.md")},
            "data": {"path": str(Path("/data/vtenishe/moon_validation_data/development_and_validation.md")),
                     "sha256": sha256(Path("/data/vtenishe/moon_validation_data/development_and_validation.md"))},
        },
        "external_data": [
            {"path": str(path), "sha256": sha256(path)} for path in data_files],
        "formal_vv_status": result["formal_vv"],
    }
    manifest_path = output / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    artifacts["run_manifest"] = {
        "path": str(manifest_path), "sha256": sha256(manifest_path)}
    result_path = output / "result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    return result


def main() -> int:
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser = argparse.ArgumentParser(
        description="Produce a production-kernel LOLA surface demonstration")
    parser.add_argument("--build-dir", type=Path, default=REPO / "build")
    parser.add_argument("--surface-resolution-m", type=float, default=100_000.0)
    parser.add_argument("--output-dir", type=Path,
                        default=REPO / "test_output/srcMoon/lola-surface-demo" /
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

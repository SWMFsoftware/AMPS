#!/usr/bin/env python3
"""Visualize srcMoon production illumination and thermal-state CSV output.

This program contains no lunar physics.  It reads the complete per-element and
per-layer tables emitted by the production-kernel demonstration adapter,
validates their schemas and finite ranges, and renders those values unchanged.
Every figure is labelled as experimental model output rather than V&V evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Iterable


NOTICE = "EXPERIMENTAL MODEL OUTPUT - NOT V&V EVIDENCE"
SURFACE_FIELDS = [
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
PROFILE_FIELDS = [
    "surface_element", "latitude_deg", "longitude_east_deg", "layer_index",
    "layer_center_depth_m", "layer_temperature_k",
]


def sha256(path: Path) -> str:
    """Return the stable identity recorded beside every input/output."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_numeric_csv(path: Path, expected: list[str]) -> list[dict[str, float]]:
    """Read a strict numeric CSV without repairing, filling, or clipping rows."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected:
            raise ValueError(
                f"{path} schema differs; expected {expected}, got {reader.fieldnames}")
        rows: list[dict[str, float]] = []
        for line_number, row in enumerate(reader, 2):
            try:
                parsed = {name: float(row[name]) for name in expected}
            except (TypeError, ValueError) as exception:
                raise ValueError(
                    f"{path}:{line_number} contains a non-numeric field") from exception
            if not all(math.isfinite(value) for value in parsed.values()):
                raise ValueError(f"{path}:{line_number} contains a non-finite value")
            rows.append(parsed)
    if not rows:
        raise ValueError(f"{path} contains no model rows")
    return rows


def validate_surface(rows: list[dict[str, float]]) -> None:
    """Reject malformed classifications/ranges rather than hiding them in plots."""
    identifiers = [int(row["surface_element"]) for row in rows]
    if identifiers != list(range(len(rows))):
        raise ValueError("surface_element must be contiguous, unique, and zero-based")
    binary = ("sun_above_local_horizon", "terrain_sun_visible",
              "terrain_occluded")
    for row in rows:
        if not (-90.0 <= row["latitude_deg"] <= 90.0 and
                0.0 <= row["longitude_east_deg"] < 360.0 and
                0.0 <= row["local_solar_time_hour"] < 24.0):
            raise ValueError("surface coordinate/local-time range is invalid")
        if any(row[name] not in (0.0, 1.0) for name in binary):
            raise ValueError("surface visibility classifications must be 0 or 1")
        for name in ("solar_disc_visible_fraction", "earth_eclipse_fraction",
                     "direct_beam_fraction"):
            if not 0.0 <= row[name] <= 1.0:
                raise ValueError(f"{name} lies outside [0,1]")
        if abs(row["solar_disc_visible_fraction"] +
               row["earth_eclipse_fraction"] - 1.0) > 2.0e-14:
            raise ValueError("solar-disc and eclipse fractions do not close")
        if row["surface_temperature_k"] <= 0.0 or row["radius_m"] <= 0.0:
            raise ValueError("temperature and radius must be positive")


def validate_profiles(surface: list[dict[str, float]],
                      profiles: list[dict[str, float]]) -> None:
    """Require one complete, monotonically deepening layer set per element."""
    by_element: dict[int, list[dict[str, float]]] = {}
    for row in profiles:
        element = int(row["surface_element"])
        if element < 0 or element >= len(surface):
            raise ValueError("thermal profile references an unknown surface element")
        by_element.setdefault(element, []).append(row)
    if set(by_element) != set(range(len(surface))):
        raise ValueError("thermal profiles do not cover every surface element")
    layer_count: int | None = None
    for element, rows in by_element.items():
        rows.sort(key=lambda row: int(row["layer_index"]))
        if [int(row["layer_index"]) for row in rows] != list(range(len(rows))):
            raise ValueError(f"surface element {element} has a layer-index gap")
        depths = [row["layer_center_depth_m"] for row in rows]
        if any(depth <= 0.0 for depth in depths) or any(
                right <= left for left, right in zip(depths, depths[1:])):
            raise ValueError(f"surface element {element} has invalid layer depths")
        if any(row["layer_temperature_k"] <= 0.0 for row in rows):
            raise ValueError(f"surface element {element} has invalid temperature")
        if layer_count is None:
            layer_count = len(rows)
        elif len(rows) != layer_count:
            raise ValueError("thermal columns do not share one configured layer count")


def values(rows: Iterable[dict[str, float]], name: str):
    """Convert one complete CSV column to NumPy only after strict validation."""
    import numpy as np

    return np.asarray([row[name] for row in rows], dtype=float)


def configure_map(axis, title: str) -> None:
    axis.set_title(title)
    axis.set_xlabel("east longitude [deg]")
    axis.set_ylabel("planetocentric latitude [deg]")
    axis.set_xlim(0.0, 360.0)
    axis.set_ylim(-90.0, 90.0)
    axis.grid(True, linewidth=0.3, color="#aaaaaa")


def render_static(surface: list[dict[str, float]],
                  profiles: list[dict[str, float]], output_dir: Path,
                  prefix: str, epoch_label: str) -> list[Path]:
    """Render complete model tables into illumination and thermal summaries."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    # Older Matplotlib releases register the ``3d`` projection only when the
    # toolkit module is imported explicitly.  The symbol itself is unused.
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    lon = values(surface, "longitude_east_deg")
    lat = values(surface, "latitude_deg")
    x = values(surface, "x_body_fixed_m")
    y = values(surface, "y_body_fixed_m")
    z = values(surface, "z_body_fixed_m")
    absorbed = values(surface, "absorbed_solar_flux_w_m2")
    eclipse = values(surface, "earth_eclipse_fraction")
    terrain = values(surface, "terrain_occluded")
    above = values(surface, "sun_above_local_horizon")
    direct = values(surface, "direct_beam_fraction")
    temperature = values(surface, "surface_temperature_k")
    local_time = values(surface, "local_solar_time_hour")

    illumination = plt.figure(figsize=(13.0, 9.0), constrained_layout=True)
    axes = illumination.subplots(2, 2)
    plot = axes[0, 0].scatter(lon, lat, c=absorbed, s=8, cmap="inferno",
                              linewidths=0)
    configure_map(axes[0, 0], "Absorbed direct-solar flux")
    illumination.colorbar(plot, ax=axes[0, 0], label="W m$^{-2}$")

    # Classification is purely a display encoding of three exported flags:
    # 0=nightside/back-facing, 1=terrain occlusion, 2=visible solar centre.
    classification = np.where(above == 0.0, 0.0,
                              np.where(terrain == 1.0, 1.0, 2.0))
    plot = axes[0, 1].scatter(lon, lat, c=classification, s=8,
                              cmap=matplotlib.colors.ListedColormap(
                                  ["#202040", "#7f3c8d", "#f9d057"]),
                              vmin=-0.5, vmax=2.5, linewidths=0)
    configure_map(axes[0, 1],
                  "Geometric nightside / terrain shadow / visible centre")
    colorbar = illumination.colorbar(plot, ax=axes[0, 1], ticks=[0, 1, 2])
    colorbar.ax.set_yticklabels(["nightside", "terrain", "visible"])

    plot = axes[1, 0].scatter(lon, lat, c=eclipse, s=8, cmap="magma",
                              vmin=0.0, vmax=1.0, linewidths=0)
    configure_map(axes[1, 0], "Finite-disc Earth occultation")
    illumination.colorbar(plot, ax=axes[1, 0], label="occulted solar-disc fraction")

    axis3d = illumination.add_subplot(2, 2, 4, projection="3d")
    axes[1, 1].remove()
    plot = axis3d.scatter(x / 1000.0, y / 1000.0, z / 1000.0, c=direct,
                          s=3, cmap="viridis", vmin=0.0, vmax=1.0,
                          linewidths=0, depthshade=False)
    axis3d.set_title("Direct-beam fraction on production LOLA surface")
    axis3d.set_xlabel("X [km]")
    axis3d.set_ylabel("Y [km]")
    axis3d.set_zlabel("Z [km]")
    extent = max(float(abs(x).max()), float(abs(y).max()), float(abs(z).max())) / 1000.0
    axis3d.set_xlim(-extent, extent)
    axis3d.set_ylim(-extent, extent)
    axis3d.set_zlim(-extent, extent)
    if hasattr(axis3d, "set_box_aspect"):
        axis3d.set_box_aspect((1, 1, 1))
    for pane_axis in (axis3d.xaxis, axis3d.yaxis, axis3d.zaxis):
        pane_axis.pane.set_alpha(1.0)
    illumination.colorbar(plot, ax=axis3d, shrink=0.65,
                          label="visible-disc fraction after terrain")
    illumination.suptitle(
        f"srcMoon illumination and terrain shadowing — {epoch_label}\n{NOTICE}")

    thermal = plt.figure(figsize=(13.0, 9.0), constrained_layout=True)
    thermal_axes = thermal.subplots(2, 2)
    plot = thermal_axes[0, 0].scatter(lon, lat, c=temperature, s=8,
                                      cmap="plasma", linewidths=0)
    configure_map(thermal_axes[0, 0], "Production surface temperature")
    thermal.colorbar(plot, ax=thermal_axes[0, 0], label="K")

    plot = thermal_axes[0, 1].scatter(local_time, temperature, c=lat, s=7,
                                      cmap="coolwarm", linewidths=0)
    thermal_axes[0, 1].set_title("Surface temperature versus local solar time")
    thermal_axes[0, 1].set_xlabel("local apparent solar time [h]")
    thermal_axes[0, 1].set_ylabel("temperature [K]")
    thermal_axes[0, 1].grid(True, linewidth=0.3, color="#aaaaaa")
    thermal.colorbar(plot, ax=thermal_axes[0, 1], label="latitude [deg]")

    thermal_axes[1, 0].hist(temperature, bins=50, color="tab:red",
                            edgecolor="black", linewidth=0.35)
    thermal_axes[1, 0].set_title("Complete surface-temperature distribution")
    thermal_axes[1, 0].set_xlabel("temperature [K]")
    thermal_axes[1, 0].set_ylabel("logical surface elements")
    thermal_axes[1, 0].grid(True, linewidth=0.3, color="#aaaaaa")

    depth_groups: dict[float, list[float]] = {}
    for row in profiles:
        depth_groups.setdefault(row["layer_center_depth_m"], []).append(
            row["layer_temperature_k"])
    depths = np.asarray(sorted(depth_groups), dtype=float)
    minimum = np.asarray([min(depth_groups[depth]) for depth in depths])
    mean = np.asarray([sum(depth_groups[depth]) / len(depth_groups[depth])
                       for depth in depths])
    maximum = np.asarray([max(depth_groups[depth]) for depth in depths])
    # Keep the envelope opaque.  Apart from making the plot unambiguous on
    # projectors, this gives the vector EPS exactly the same rendering intent
    # as the PNG; PostScript cannot represent partial transparency.
    thermal_axes[1, 1].fill_betweenx(depths, minimum, maximum,
                                     color="#c6dbef", alpha=1.0,
                                     label="global min–max")
    thermal_axes[1, 1].plot(mean, depths, color="tab:blue", label="global mean")
    thermal_axes[1, 1].invert_yaxis()
    thermal_axes[1, 1].set_title("All-column subsurface thermal state")
    thermal_axes[1, 1].set_xlabel("temperature [K]")
    thermal_axes[1, 1].set_ylabel("layer-centre depth [m]")
    thermal_axes[1, 1].grid(True, linewidth=0.3, color="#aaaaaa")
    # Matplotlib defaults to an 80%-opaque legend frame.  Use an opaque frame
    # so the PNG and transparency-free EPS representations remain identical.
    thermal_axes[1, 1].legend(framealpha=1.0)
    thermal.suptitle(
        f"srcMoon depth-resolved lunar thermal state — {epoch_label}\n{NOTICE}")

    outputs: list[Path] = []
    for figure, stem in ((illumination, "illumination_shadowing"),
                         (thermal, "thermal_state")):
        for suffix in ("png", "eps"):
            path = output_dir / f"{prefix}_{stem}.{suffix}"
            figure.savefig(path, dpi=180 if suffix == "png" else None,
                           format=suffix, transparent=False)
            outputs.append(path)
        plt.close(figure)
    return outputs


def interactive_globe(surface: list[dict[str, float]], color_by: str,
                      selected_output: Path | None, force: bool) -> dict[str, object]:
    """Open a rotatable point globe and optionally save its selected camera."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    x = values(surface, "x_body_fixed_m") / 1000.0
    y = values(surface, "y_body_fixed_m") / 1000.0
    z = values(surface, "z_body_fixed_m") / 1000.0
    color = values(surface, color_by)
    figure = plt.figure(figsize=(9.0, 8.0))
    axis = figure.add_subplot(111, projection="3d")
    # Disable Matplotlib's artificial distance-dependent alpha modulation.
    # The displayed opacity must not be mistaken for a physical visibility
    # or attenuation quantity, and the selected view may be saved as EPS.
    plot = axis.scatter(x, y, z, c=color, s=4, cmap="viridis", linewidths=0,
                        depthshade=False)
    extent = max(float(abs(x).max()), float(abs(y).max()), float(abs(z).max()))
    axis.set_xlim(-extent, extent)
    axis.set_ylim(-extent, extent)
    axis.set_zlim(-extent, extent)
    if hasattr(axis, "set_box_aspect"):
        axis.set_box_aspect((1, 1, 1))
    axis.set_xlabel("X [km]")
    axis.set_ylabel("Y [km]")
    axis.set_zlabel("Z [km]")
    axis.set_title(f"{color_by}\n{NOTICE}\nDrag to rotate; press S to save, Q to close")
    figure.colorbar(plot, ax=axis, shrink=0.7, label=color_by)
    saved = False

    def on_key(event) -> None:
        nonlocal saved
        if event.key and event.key.lower() == "q":
            plt.close(figure)
        elif event.key and event.key.lower() == "s" and selected_output:
            if selected_output.exists() and not force and not saved:
                print(f"Refusing to replace {selected_output}; use --force",
                      file=sys.stderr)
                return
            figure.savefig(selected_output, dpi=180, transparent=False)
            saved = True
            print(f"Saved selected view: {selected_output}")

    figure.canvas.mpl_connect("key_press_event", on_key)
    plt.show()
    return {"color_by": color_by, "selected_view_saved": saved,
            "selected_view": str(selected_output) if selected_output else None,
            "elevation_deg": float(axis.elev), "azimuth_deg": float(axis.azim)}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Visualize production srcMoon illumination/shadow/thermal CSVs",
        epilog=("Static PNG/EPS figures are always produced. With --interactive-3d, "
                "drag to rotate, press S to save --selected-view-output, and Q "
                "to close. Figures are model demonstrations, not V&V evidence."))
    parser.add_argument("surface_csv", type=Path,
                        help="surface_environment.csv from environment_demo/run.py")
    parser.add_argument("thermal_profiles_csv", type=Path,
                        help="thermal_profiles.csv from the same run")
    parser.add_argument("--output-dir", type=Path,
                        help="default: directory containing surface_csv")
    parser.add_argument("--prefix", default="moon_environment")
    parser.add_argument("--epoch-label", default="epoch recorded in run manifest")
    parser.add_argument("--force", action="store_true",
                        help="replace existing figure/receipt files")
    parser.add_argument("--interactive-3d", action="store_true",
                        help="open an additional mouse-rotatable 3-D point globe")
    parser.add_argument("--color-by", choices=[
        "absorbed_solar_flux_w_m2", "direct_beam_fraction",
        "terrain_occluded", "earth_eclipse_fraction", "surface_temperature_k"],
        default="surface_temperature_k",
        help="quantity shown on the interactive globe")
    parser.add_argument("--selected-view-output", type=Path,
                        help="path written by S in interactive mode")
    parser.add_argument("--receipt", type=Path,
                        help="default: <output-dir>/<prefix>_visualization.json")
    arguments = parser.parse_args()
    try:
        surface_path = arguments.surface_csv.resolve()
        profiles_path = arguments.thermal_profiles_csv.resolve()
        output_dir = (arguments.output_dir or surface_path.parent).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        receipt = (arguments.receipt or
                   output_dir / f"{arguments.prefix}_visualization.json").resolve()
        planned = [
            output_dir / f"{arguments.prefix}_illumination_shadowing.png",
            output_dir / f"{arguments.prefix}_illumination_shadowing.eps",
            output_dir / f"{arguments.prefix}_thermal_state.png",
            output_dir / f"{arguments.prefix}_thermal_state.eps", receipt,
        ]
        existing = [path for path in planned if path.exists()]
        if existing and not arguments.force:
            raise ValueError("refusing to replace existing outputs: " +
                             ", ".join(str(path) for path in existing))
        surface = read_numeric_csv(surface_path, SURFACE_FIELDS)
        profiles = read_numeric_csv(profiles_path, PROFILE_FIELDS)
        validate_surface(surface)
        validate_profiles(surface, profiles)
        outputs = render_static(surface, profiles, output_dir, arguments.prefix,
                                arguments.epoch_label)
        interaction = None
        if arguments.interactive_3d:
            interaction = interactive_globe(
                surface, arguments.color_by,
                arguments.selected_view_output.resolve()
                if arguments.selected_view_output else None,
                arguments.force)
        result = {
            "schema": "srcMoon-environment-visualization-v1",
            "notice": NOTICE,
            "inputs": {
                "surface": {"path": str(surface_path), "sha256": sha256(surface_path)},
                "thermal_profiles": {"path": str(profiles_path),
                                     "sha256": sha256(profiles_path)},
            },
            "counts": {"surface_elements": len(surface),
                       "thermal_profile_rows": len(profiles)},
            "outputs": [{"path": str(path), "sha256": sha256(path)}
                        for path in outputs],
            "interactive": interaction,
        }
        receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, ImportError) as exception:
        print(f"ERROR: {exception}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

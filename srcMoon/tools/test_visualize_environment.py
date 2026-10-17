#!/usr/bin/env python3
"""Infrastructure test for visualize_environment.py; not lunar V&V."""

from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="moon-environment-visualizer-") as raw:
        root = Path(raw)
        surface = root / "surface.csv"
        profiles = root / "profiles.csv"
        surface_fields = [
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
        rows = [
            [0, 0, 0, 0, 0, 12, 1737400, 0, 0, 1737400, 0,
             1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1000, 350, 0, 4],
            [1, 1, 0, 0, 180, 0, -1737400, 0, 0, 1737400, 0,
             -1, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 100, 0, 4],
        ]
        with surface.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(surface_fields)
            writer.writerows(rows)
        with profiles.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["surface_element", "latitude_deg",
                             "longitude_east_deg", "layer_index",
                             "layer_center_depth_m", "layer_temperature_k"])
            writer.writerows([[0, 0, 0, 0, 0.5, 300], [0, 0, 0, 1, 1.5, 250],
                              [1, 0, 180, 0, 0.5, 120],
                              [1, 0, 180, 1, 1.5, 180]])
        command = [sys.executable, "-B", str(HERE / "visualize_environment.py"),
                   str(surface), str(profiles), "--output-dir", str(root),
                   "--epoch-label", "analytical infrastructure fixture"]
        completed = subprocess.run(command, text=True, capture_output=True,
                                   check=False)
        if completed.returncode != 0:
            print(completed.stdout)
            print(completed.stderr, file=sys.stderr)
            return 1
        required = [root / "moon_environment_illumination_shadowing.png",
                    root / "moon_environment_illumination_shadowing.eps",
                    root / "moon_environment_thermal_state.png",
                    root / "moon_environment_thermal_state.eps",
                    root / "moon_environment_visualization.json"]
        if not all(path.is_file() and path.stat().st_size > 0 for path in required):
            print("FAIL: visualizer did not create every required artifact")
            return 1
    print("PASS: environment visualizer infrastructure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Plot the Figure 2 map and Figure 3 curves requested for Stage 0."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates


CURRENT_IAU_WEST_SHIFT_DEG = 3.5
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"]


def polar_xy(latitude_deg, west_longitude_deg, shift=True):
    latitude = np.asarray(latitude_deg, dtype=float)
    west = np.asarray(west_longitude_deg, dtype=float)
    if shift:
        west = west + CURRENT_IAU_WEST_SHIFT_DEG
    radius = 90.0 + latitude
    angle = np.deg2rad(west)
    return radius * np.sin(angle), -radius * np.cos(angle)


def load_south_polar_mosaic(path, size=1400, latitude_limit=-65.0):
    Image.MAX_IMAGE_PIXELS = None
    image = Image.open(path)
    width, height = image.size
    radial_limit = 90.0 + latitude_limit
    y0 = max(0, int((90.0 - latitude_limit) / 180.0 * height) - 3)
    south = np.asarray(image.crop((0, y0, width, height)).convert("L"), dtype=float)
    south = np.pad(south, ((0, 0), (2, 2)), mode="wrap")

    horizontal = np.linspace(-radial_limit, radial_limit, size)
    vertical = np.linspace(radial_limit, -radial_limit, size)
    x, y = np.meshgrid(horizontal, vertical)
    radius = np.hypot(x, y)
    latitude = -90.0 + radius
    west = np.mod(np.degrees(np.arctan2(x, -y)), 360.0)
    east = np.mod(-west, 360.0)
    column = east / 360.0 * width + 2.0
    row = (90.0 - latitude) / 180.0 * height - y0
    sampled = map_coordinates(south, [row, column], order=1, mode="nearest")
    sampled[radius > radial_limit] = np.nan
    valid = sampled[np.isfinite(sampled)]
    low, high = np.percentile(valid, [1.0, 99.5])
    sampled = np.clip((sampled - low) / (high - low), 0.0, 1.0) ** 0.72
    return sampled, radial_limit


def load_track(path):
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    latitude = np.asarray([float(r["lat_planetocentric_deg"]) for r in rows])
    longitude = np.asarray([float(r["west_lon_planetocentric_deg"]) for r in rows])
    time = np.asarray([float(r["time_from_reference_s"]) for r in rows])
    return latitude, longitude, time


def stripe_lat_lon(points):
    points = np.asarray(points, dtype=float)
    radius = np.linalg.norm(points, axis=1)
    latitude = np.degrees(np.arcsin(points[:, 2] / radius))
    east = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    return latitude, np.mod(-east, 360.0)


def decorate_polar_grid(ax, limit):
    theta = np.linspace(0, 2 * np.pi, 721)
    for latitude in (-85, -80, -75, -70, -65):
        radius = 90.0 + latitude
        ax.plot(radius * np.sin(theta), -radius * np.cos(theta),
                color="white", lw=0.45, alpha=0.45, zorder=2)
        ax.text(0.35, -radius + 0.15, "%d°" % latitude, color="white",
                fontsize=7, ha="left", va="bottom", zorder=8)
    for west in range(0, 360, 30):
        angle = np.deg2rad(west)
        ax.plot([0, limit * np.sin(angle)], [0, -limit * np.cos(angle)],
                color="white", lw=0.35, alpha=0.34, zorder=2)
        x = (limit + 1.0) * np.sin(angle)
        y = -(limit + 1.0) * np.cos(angle)
        ax.text(x, y, "%d°W" % west, color="0.15", fontsize=7,
                ha="center", va="center", zorder=8)


def plot_figure2(args):
    background, limit = load_south_polar_mosaic(args.mosaic)
    with open(args.config) as stream:
        config = json.load(stream)
    with open(args.stripes) as stream:
        stripes = json.load(stream)["stripes"]

    fig, ax = plt.subplots(figsize=(9.2, 9.2), constrained_layout=True)
    ax.imshow(background, extent=(-limit, limit, -limit, limit), origin="upper",
              cmap="gray", vmin=0, vmax=1, zorder=0)
    decorate_polar_grid(ax, limit)

    stripe_families = {item["geometry_id"]: item["family"]
                       for item in config["tiger_stripes"]}
    stripe_handles = {}
    for name, points in stripes.items():
        latitude, longitude = stripe_lat_lon(points)
        x, y = polar_xy(latitude, longitude)
        family = stripe_families[name]
        style = dict(color="#00e5ff", lw=2.0, ls="-") if family == "normal" else dict(
            color="#ff4fd8", lw=1.6, ls="--")
        handle, = ax.plot(x, y, zorder=5, **style)
        stripe_handles.setdefault(family, handle)

    for index, source in enumerate(config["discrete_jets"]):
        x, y = polar_xy(source["latitude_deg"], source["west_longitude_deg"])
        ax.scatter(x, y, s=62, marker="o", facecolor="#ffe45e", edgecolor="black",
                   linewidth=0.9, zorder=7)
        radial = np.hypot(x, y)
        dx, dy = (0.55 * x / radial, 0.55 * y / radial) if radial else (0.5, 0.5)
        ax.text(x + dx, y + dy, ROMAN[index], fontsize=8, fontweight="bold",
                color="white", ha="center", va="center", zorder=8,
                bbox=dict(boxstyle="round,pad=0.12", facecolor="black",
                          edgecolor="none", alpha=0.62))

    track_handles = []
    for flyby, path, color in (("E3", args.e3_track, "#ff2b2b"),
                               ("E5", args.e5_track, "#ff9f1c")):
        latitude, longitude, time = load_track(path)
        keep = latitude <= -65.0
        x, y = polar_xy(latitude[keep], longitude[keep])
        handle, = ax.plot(x, y, color=color, lw=2.3, zorder=6, label=flyby)
        track_handles.append(handle)
        if len(x):
            middle = len(x) // 2
            ax.annotate(flyby, (x[middle], y[middle]), xytext=(6, 6),
                        textcoords="offset points", color="white", fontsize=9,
                        fontweight="bold", zorder=9,
                        bbox=dict(boxstyle="round,pad=0.18", facecolor=color,
                                  edgecolor="white", linewidth=0.6, alpha=0.9))

    vent = ax.scatter([], [], s=62, facecolor="#ffe45e", edgecolor="black",
                      label="Discrete source")
    handles = track_handles + [stripe_handles["normal"], stripe_handles["extra"], vent]
    labels = ["E3 SPICE ground track", "E5 SPICE ground track",
              "Tiger stripes", "Extra model fissures", "Discrete sources I–VIII"]
    ax.legend(handles, labels, loc="upper left", fontsize=8, framealpha=0.86)
    ax.set_xlim(-limit - 2.2, limit + 2.2)
    ax.set_ylim(-limit - 2.2, limit + 2.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Enceladus south pole: SPICE flyby tracks and 2014 plume sources\n"
                 "USGS Cassini ISS controlled mosaic; current-IAU basemap",
                 fontsize=13, pad=10)
    fig.text(0.5, 0.012,
             "Event-PCK coordinates shifted +3.5° west only for alignment with the current-IAU mosaic.",
             ha="center", fontsize=8)
    save_pair(fig, args.output_dir / "figure2_ground_tracks_sources")
    plt.close(fig)


def load_figure3(path):
    grouped = {}
    with open(path, newline="") as stream:
        for row in csv.DictReader(stream):
            key = (row["flyby"], row["series"])
            grouped.setdefault(key, [[], []])
            grouped[key][0].append(float(row["time_ca_s"]))
            grouped[key][1].append(float(row["density_m3"]))
    return {key: (np.asarray(value[0]), np.asarray(value[1]))
            for key, value in grouped.items()}


def plot_figure3(args):
    curves = load_figure3(args.figure3_curves)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True,
                             sharex=True)
    styles = {
        "INMS": dict(color="black", lw=1.0, label="INMS observations"),
        "with_tiger_stripes": dict(color="#e31a1c", lw=2.0,
                                   label="Model with tiger stripes"),
        "without_tiger_stripes": dict(color="#1769e0", lw=1.8,
                                      label="Model without tiger stripes"),
    }
    for ax, flyby in zip(axes, ("E3", "E5")):
        for series in ("INMS", "without_tiger_stripes", "with_tiger_stripes"):
            time, density = curves[(flyby, series)]
            order = np.argsort(time)
            ax.semilogy(time[order], density[order], **styles[series])
        ax.set_xlim(-400, 400)
        ax.set_ylim((1e9, 3e14) if flyby == "E3" else (1e10, 5e14))
        ax.set_title(flyby, fontweight="bold")
        ax.set_xlabel("Time from flyby reference (s)")
        ax.grid(True, which="both", color="0.88", lw=0.55)
    axes[0].set_ylabel(r"H$_2$O number density (m$^{-3}$)")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.05))
    fig.suptitle("Reconstruction of Tenishev et al. (2014), Figure 3",
                 y=1.12, fontsize=13)
    fig.text(0.5, -0.035,
             "Curves recovered from the published PDF vector paths; this is not a rerun of the 2014 simulation.",
             ha="center", fontsize=8)
    save_pair(fig, args.output_dir / "figure3_published_density_comparison")
    plt.close(fig)


def save_pair(figure, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [stem.with_suffix(".png"), stem.with_suffix(".pdf")]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise FileExistsError("refusing to overwrite: " + ", ".join(existing))
    figure.savefig(outputs[0], dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(outputs[1], bbox_inches="tight", facecolor="white")


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--mosaic", type=Path,
                        default=root / "assets/Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif")
    parser.add_argument("--config", type=Path,
                        default=root / "config/tenishev2014_stage0.json")
    parser.add_argument("--stripes", type=Path,
                        default=root / "data/published_stripes_2014.json")
    parser.add_argument("--e3-track", type=Path,
                        default=root / "figure_data/e3_spice_ground_track.csv")
    parser.add_argument("--e5-track", type=Path,
                        default=root / "figure_data/e5_spice_ground_track.csv")
    parser.add_argument("--figure3-curves", type=Path,
                        default=root / "figure_data/figure3_published_curves.csv")
    parser.add_argument("--output-dir", type=Path, default=root / "figures")
    args = parser.parse_args()
    plot_figure2(args)
    plot_figure3(args)


if __name__ == "__main__":
    main()

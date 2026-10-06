#!/usr/bin/env python3
"""Interactive 3-D viewer for srcSEP3D reduced-front Tecplot surfaces.

The maintained writer uses an ASCII FETRIANGLE zone with nodal geometry and
cell-centred shock state.  The reader also accepts historical POINT-packed
FEQUADRILATERAL evidence so changing the production topology does not make old
qualification products unreadable. Some VisIt installations do not recognize
either standalone surface as a
standalone database, so this tool reads the file directly and displays it with
Matplotlib.  It is a visualization utility: it does not reconstruct a sheath,
ejecta, or any spatial downstream CME volume.

New triangle values are already face-centred and are never averaged across a
classification boundary. For a legacy quad, values are arithmetic means of
the four nodes only when their ``shock_accepted`` classification agrees; a
mixed legacy quad remains neutral gray and is removed by accepted-only mode.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import glob
import math
import os
from pathlib import Path
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


SOLAR_RADIUS_M = 6.957e8
ASTRONOMICAL_UNIT_M = 149_597_870_700.0


class FrontFormatError(RuntimeError):
    """The input is not the finite reduced-front Tecplot subset we support."""


@dataclass(frozen=True)
class FrontSurface:
    """Parsed immutable surface in the HCI/SI convention of the front writer."""

    path: Path
    title: str
    variables: Tuple[str, ...]
    data: np.ndarray
    faces: np.ndarray
    cell_data: np.ndarray
    cell_centered: Tuple[bool, ...]
    zone_type: str
    auxiliary: Dict[str, str]

    def column(self, name: str) -> np.ndarray:
        try:
            index = self.variables.index(name)
        except ValueError as error:
            raise FrontFormatError(
                f"unknown variable {name!r}; use --list-variables"
            ) from error
        return self.cell_data[:, index] if self.cell_centered[index] \
            else self.data[:, index]

    @property
    def quads(self) -> np.ndarray:
        """Compatibility alias for callers that predate triangular output."""

        return self.faces


def _parse_variables(line: str) -> Tuple[str, ...]:
    if not line.startswith("VARIABLES="):
        raise FrontFormatError("line 2 must begin with VARIABLES=")
    try:
        names = next(csv.reader([line.split("=", 1)[1]], skipinitialspace=True))
    except (csv.Error, StopIteration) as error:
        raise FrontFormatError(f"cannot parse VARIABLES record: {error}") from error
    names = tuple(name.strip() for name in names)
    if not names or any(not name for name in names) or len(set(names)) != len(names):
        raise FrontFormatError("VARIABLES names must be nonempty and unique")
    return names


def load_front(path: Path) -> FrontSurface:
    """Parse one maintained triangle zone or historical quadrilateral zone.

    The geometry is not guessed from row ordering.  Connectivity is read from
    the file and converted from Tecplot's one-based node indices to NumPy's
    zero-based convention.  Strict counts and finite geometry keep a truncated
    output from looking like a legitimate physical surface.
    """

    path = path.expanduser().resolve()
    try:
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    except OSError as error:
        raise FrontFormatError(f"cannot open {path}: {error}") from error
    if len(lines) < 4:
        raise FrontFormatError("surface file is truncated before its data zone")
    title_line = lines[0]
    variables = _parse_variables(lines[1])
    zone_line = lines[2].strip()
    if not title_line.startswith("TITLE="):
        raise FrontFormatError("line 1 must begin with TITLE=")
    is_triangle = "ZONETYPE=FETRIANGLE" in zone_line
    is_legacy_quad = "ZONETYPE=FEQUADRILATERAL" in zone_line
    block = "DATAPACKING=BLOCK" in zone_line
    point = "DATAPACKING=POINT" in zone_line
    if not ((is_triangle and block) or (is_legacy_quad and point)):
        raise FrontFormatError(
            "ZONE must be BLOCK/FETRIANGLE or legacy POINT/FEQUADRILATERAL"
        )
    node_match = re.search(r"(?:^|[, ])N=(\d+)(?:,| |$)", zone_line)
    element_match = re.search(r"(?:^|[, ])E=(\d+)(?:,| |$)", zone_line)
    if not node_match or not element_match:
        raise FrontFormatError("ZONE does not declare integer N and E")
    node_count = int(node_match.group(1))
    element_count = int(element_match.group(1))
    if node_count <= 0 or element_count <= 0:
        raise FrontFormatError("ZONE N and E must both be positive")

    auxiliary: Dict[str, str] = {}
    line_index = 3
    while line_index < len(lines):
        stripped = lines[line_index].strip()
        if not stripped:
            line_index += 1
            continue
        if not stripped.startswith("AUXDATA"):
            break
        match = re.fullmatch(r'AUXDATA\s+([^=\s]+)="(.*)"', stripped)
        if not match:
            raise FrontFormatError(f"malformed AUXDATA record: {stripped}")
        auxiliary[match.group(1)] = match.group(2)
        line_index += 1
    tokens = " ".join(lines[line_index:]).split()
    if not tokens:
        raise FrontFormatError("ZONE contains no numerical data")

    variable_count = len(variables)
    data = np.full((node_count, variable_count), np.nan, dtype=float)
    cell_data = np.full((element_count, variable_count), np.nan, dtype=float)
    cell_centered = [False] * variable_count
    cursor = 0
    if is_triangle:
        location = re.search(
            r"VARLOCATION=\(\[4-(\d+)\]=CELLCENTERED\)", zone_line
        )
        if not location or int(location.group(1)) != variable_count:
            raise FrontFormatError(
                "triangle zone must declare variables 4..N as CELLCENTERED"
            )
        cell_centered[3:] = [True] * (variable_count - 3)
        for variable in range(variable_count):
            count = element_count if cell_centered[variable] else node_count
            try:
                values = [float(word) for word in tokens[cursor:cursor + count]]
            except ValueError as error:
                raise FrontFormatError(
                    f"variable block {variables[variable]!r} is non-numeric"
                ) from error
            if len(values) != count:
                raise FrontFormatError(
                    f"variable block {variables[variable]!r} is truncated"
                )
            if cell_centered[variable]:
                cell_data[:, variable] = values
            else:
                data[:, variable] = values
            cursor += count
        nodes_per_face = 3
        zone_type = "FETRIANGLE"
    else:
        value_count = node_count * variable_count
        try:
            values = [float(word) for word in tokens[:value_count]]
        except ValueError as error:
            raise FrontFormatError("legacy point records contain non-numeric data") from error
        if len(values) != value_count:
            raise FrontFormatError("legacy point records are truncated")
        data[:, :] = np.asarray(values).reshape(node_count, variable_count)
        cursor = value_count
        nodes_per_face = 4
        zone_type = "FEQUADRILATERAL"

    connectivity = []
    for element in range(element_count):
        words = tokens[cursor:cursor + nodes_per_face]
        if len(words) != nodes_per_face:
            raise FrontFormatError(f"element {element + 1} connectivity is truncated")
        try:
            face = [int(word) - 1 for word in words]
        except ValueError as error:
            raise FrontFormatError(
                f"element {element + 1} has a non-integer node index"
            ) from error
        if any(index < 0 or index >= node_count for index in face):
            raise FrontFormatError(
                f"element {element + 1} references a node outside 1..{node_count}"
            )
        connectivity.append(face)
        cursor += nodes_per_face
    if cursor != len(tokens):
        raise FrontFormatError("unexpected values follow surface connectivity")

    faces = np.asarray(connectivity, dtype=np.int64)
    for coordinate in ("x_m", "y_m", "z_m", "shock_accepted"):
        if coordinate not in variables:
            raise FrontFormatError(f"required variable {coordinate!r} is absent")
    geometry = data[:, [variables.index(name) for name in ("x_m", "y_m", "z_m")]]
    if not np.all(np.isfinite(geometry)):
        raise FrontFormatError("surface coordinates contain NaN or infinity")
    return FrontSurface(path, title_line, variables, data, faces, cell_data,
                        tuple(cell_centered), zone_type, auxiliary)


def classify_faces(surface: FrontSurface) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return accepted, non-shock, and legacy mixed-boundary face masks.

    ``shock_accepted`` is authoritative.  Status codes are deliberately not
    decoded here because their integer values are an implementation detail and
    because finite no-shock display values must never turn a patch into a shock.
    Requiring unanimous vertices is conservative at a classification boundary.
    """

    node_accepted = surface.column("shock_accepted") >= 0.5
    if node_accepted.shape[0] == surface.faces.shape[0]:
        accepted = node_accepted
        return accepted, ~accepted, np.zeros_like(accepted, dtype=bool)
    face_nodes = node_accepted[surface.faces]
    accepted = np.all(face_nodes, axis=1)
    nonshock = np.all(~face_nodes, axis=1)
    boundary = ~(accepted | nonshock)
    return accepted, nonshock, boundary


def face_values(surface: FrontSurface, variable: str) -> np.ndarray:
    """Return face data directly, or average a historical nodal field."""

    values = surface.column(variable)
    if values.shape[0] == surface.faces.shape[0]:
        return values
    return np.mean(values[surface.faces], axis=1)


def _length_scale(unit: str) -> Tuple[float, str]:
    if unit == "m":
        return 1.0, "m"
    if unit == "rs":
        return SOLAR_RADIUS_M, "R_sun"
    if unit == "au":
        return ASTRONOMICAL_UNIT_M, "AU"
    raise ValueError(unit)


def expand_front_paths(arguments: Sequence[str]) -> Tuple[Path, ...]:
    """Expand explicit files, directories, and quoted shell-style patterns.

    Directory input is intentionally limited to ``*-front.dat``.  This avoids
    accidentally parsing the much larger native ambient volume outputs as
    front surfaces.  Paths are de-duplicated after resolution so overlapping
    patterns do not duplicate movie frames.
    """

    paths: List[Path] = []
    for argument in arguments:
        candidate = Path(argument).expanduser()
        if candidate.is_dir():
            paths.extend(sorted(candidate.glob("*-front.dat")))
        elif glob.has_magic(argument):
            paths.extend(Path(item) for item in sorted(glob.glob(argument)))
        else:
            paths.append(candidate)
    unique = []
    seen = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            unique.append(resolved)
            seen.add(resolved)
    if not unique:
        raise FrontFormatError("no *-front.dat files matched the input")
    return tuple(unique)


def surface_time(surface: FrontSurface) -> float:
    """Return the unique physical epoch stored on every surface node."""

    values = surface.column("time_s")
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise FrontFormatError(f"{surface.path.name}: time_s is absent or non-finite")
    scale = max(1.0, abs(float(values[0])))
    if np.max(np.abs(values - values[0])) > 1.0e-12 * scale:
        raise FrontFormatError(f"{surface.path.name}: node epochs are inconsistent")
    return float(values[0])


def surface_generation(surface: FrontSurface) -> int:
    values = surface.column("generation")
    if values.size == 0 or not np.all(np.isfinite(values)) or \
            np.any(values != values[0]):
        raise FrontFormatError(f"{surface.path.name}: node generations are inconsistent")
    return int(values[0])


def sort_movie_surfaces(surfaces: Sequence[FrontSurface]) -> Tuple[FrontSurface, ...]:
    """Sort by embedded time, rejecting ambiguous duplicate movie epochs."""

    ordered = tuple(sorted(surfaces, key=lambda item: (
        surface_time(item), surface_generation(item), item.path.name)))
    times = [surface_time(item) for item in ordered]
    for left, right in zip(times, times[1:]):
        if right <= left:
            raise FrontFormatError(
                "movie inputs contain duplicate/non-increasing time_s values; "
                "select one MPI/cadence product series"
            )
    return ordered


def _coordinates(surface: FrontSurface, scale: float) -> np.ndarray:
    indices = [surface.variables.index(name) for name in ("x_m", "y_m", "z_m")]
    return surface.data[:, indices] / scale


def _shock_axis_direction(surface: FrontSurface) -> np.ndarray:
    """Infer the Sun-to-front symmetry direction from the sampled surface.

    The reduced output contains the complete finite surface but not a redundant
    event-axis column.  For the selected SSE geometry, the area-weighted surface
    centroid lies exactly on the symmetry axis in the continuous quadrature
    limit.  Equal-area sampling makes this a stable discrete estimator.  If a
    future output omits usable quadrature weights, equal node weights preserve
    the same symmetry identity.  This is a plotting construction and never
    feeds back into provider geometry or shock classification.
    """

    xyz = _coordinates(surface, 1.0)
    if "quadrature_area_m2" in surface.variables:
        weights = surface.column("quadrature_area_m2")
        if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
            raise FrontFormatError(
                f"{surface.path.name}: quadrature areas must be finite and positive"
            )
        # Production triangle areas are cell-centred, so weight chord
        # centroids. Historical quadrilateral products stored one area per
        # node and retain their original node-weighted estimate.
        points = np.mean(xyz[surface.faces], axis=1) \
            if weights.shape[0] == surface.faces.shape[0] else xyz
    else:
        points = xyz
        weights = np.ones(points.shape[0])
    centroid = np.average(points, axis=0, weights=weights)
    magnitude = float(np.linalg.norm(centroid))
    if not np.isfinite(magnitude) or magnitude == 0.0:
        raise FrontFormatError(
            f"{surface.path.name}: cannot infer a Sun-to-front symmetry axis"
        )
    return centroid / magnitude


def _nice_radial_step(maximum_rs: float) -> float:
    """Choose about six legible 1/2/2.5/5-decade radial tick intervals."""

    raw = max(maximum_rs / 6.0, 1.0e-12)
    decade = 10.0 ** math.floor(math.log10(raw))
    fraction = raw / decade
    for candidate in (1.0, 2.0, 2.5, 5.0, 10.0):
        if fraction <= candidate:
            return candidate * decade
    return 10.0 * decade


def _draw_sun(axis, scale: float) -> None:
    """Draw the one-solar-radius reference sphere at the HCI origin."""

    longitude = np.linspace(0.0, 2.0 * np.pi, 48)
    colatitude = np.linspace(0.0, np.pi, 24)
    radius = SOLAR_RADIUS_M / scale
    x = radius * np.outer(np.cos(longitude), np.sin(colatitude))
    y = radius * np.outer(np.sin(longitude), np.sin(colatitude))
    z = radius * np.outer(np.ones_like(longitude), np.cos(colatitude))
    axis.plot_surface(x, y, z, color="#f5b642", edgecolor="#d98c10",
                      linewidth=0.12, alpha=0.95, shade=True)


def _draw_distance_axis(axis, surface: FrontSurface, scale: float,
                        maximum_rs_override: Optional[float] = None) -> None:
    """Draw a Sun-centered radial ruler through the middle of the front.

    Locations are transformed to the selected plot unit, but labels are always
    heliocentric solar radii as requested.  The line ends at the farthest node
    projected onto the inferred symmetry axis, which approximates the apex.
    A movie supplies the series-wide maximum so ruler extent and ticks remain
    fixed while its direction continues to pass through each frame's middle.
    """

    direction = _shock_axis_direction(surface)
    xyz_m = _coordinates(surface, 1.0)
    local_maximum_rs = float(np.max(xyz_m @ direction) / SOLAR_RADIUS_M)
    maximum_rs = local_maximum_rs if maximum_rs_override is None else \
        maximum_rs_override
    if maximum_rs <= 0.0:
        raise FrontFormatError(f"{surface.path.name}: front axis has no positive extent")
    end = direction * (1.035 * maximum_rs * SOLAR_RADIUS_M / scale)
    axis.plot([0.0, end[0]], [0.0, end[1]], [0.0, end[2]],
              color="black", linewidth=1.35, zorder=20)

    reference = np.array([0.0, 0.0, 1.0])
    if abs(float(direction @ reference)) > 0.92:
        reference = np.array([0.0, 1.0, 0.0])
    perpendicular = np.cross(direction, reference)
    perpendicular /= np.linalg.norm(perpendicular)
    tick_half_length = max(0.006 * maximum_rs * SOLAR_RADIUS_M / scale,
                           0.04 * SOLAR_RADIUS_M / scale)
    label_offset = 2.2 * tick_half_length
    step = _nice_radial_step(maximum_rs)
    ticks = np.arange(0.0, maximum_rs + 0.25 * step, step)
    for tick_rs in ticks:
        point = direction * (tick_rs * SOLAR_RADIUS_M / scale)
        low = point - tick_half_length * perpendicular
        high = point + tick_half_length * perpendicular
        axis.plot([low[0], high[0]], [low[1], high[1]], [low[2], high[2]],
                  color="black", linewidth=0.8, zorder=21)
        label = point + label_offset * perpendicular
        axis.text(label[0], label[1], label[2], f"{tick_rs:g} R_sun",
                  color="black", fontsize=7, zorder=22)


def _global_bounds(surfaces: Sequence[FrontSurface], scale: float,
                   include_origin: bool) -> Tuple[np.ndarray, float]:
    """Return one fixed equal-aspect HCI cube for every requested frame."""

    coordinates = np.concatenate([_coordinates(item, scale) for item in surfaces],
                                 axis=0)
    if include_origin:
        solar_radius = SOLAR_RADIUS_M / scale
        references = np.vstack((np.zeros(3), np.eye(3) * solar_radius,
                                -np.eye(3) * solar_radius))
        coordinates = np.vstack((coordinates, references))
    minima = np.min(coordinates, axis=0)
    maxima = np.max(coordinates, axis=0)
    center = 0.5 * (minima + maxima)
    half_span = 0.525 * float(np.max(maxima - minima))
    if not np.isfinite(half_span) or half_span <= 0.0:
        raise FrontFormatError("surface series has zero or invalid spatial extent")
    return center, half_span


def _normalization_values(surfaces: Sequence[FrontSurface], variable: str) -> np.ndarray:
    values = []
    for surface in surfaces:
        if variable not in surface.variables:
            raise FrontFormatError(
                f"{surface.path.name}: unknown variable {variable!r}; "
                "use --list-variables"
            )
        accepted, nonshock, _boundary = classify_faces(surface)
        local = face_values(surface, variable)
        selected = local[(accepted | nonshock) & np.isfinite(local)]
        if selected.size:
            values.append(selected)
    if not values:
        raise FrontFormatError(f"{variable!r} has no finite homogeneous-face values")
    return np.concatenate(values)


def _color_limits(values: np.ndarray, vmin: Optional[float],
                  vmax: Optional[float]) -> Tuple[float, float]:
    lower = float(np.min(values)) if vmin is None else vmin
    upper = float(np.max(values)) if vmax is None else vmax
    if not np.isfinite(lower) or not np.isfinite(upper) or upper < lower:
        raise FrontFormatError("--vmin/--vmax must define a finite ordered interval")
    if upper == lower:
        padding = max(1.0, abs(lower)) * 1.0e-12
        lower -= padding
        upper += padding
    return lower, upper


def _draw_frame(axis, surface: FrontSurface, variable: str,
                accepted_only: bool, scale: float, unit_label: str,
                cmap, norm, center: np.ndarray, half_span: float,
                elevation: float, azimuth: float, show_sun: bool,
                show_distance_axis: bool,
                distance_axis_maximum_rs: Optional[float] = None):
    """Draw one frame with frozen camera/bounds and return face statistics."""

    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    axis.cla()
    coordinates = _coordinates(surface, scale)
    polygons = coordinates[surface.faces]
    values = face_values(surface, variable)
    accepted, nonshock, boundary = classify_faces(surface)
    homogeneous = accepted | nonshock
    finite = np.isfinite(values)
    visible = (accepted if accepted_only else homogeneous) & finite
    if np.any(visible):
        colored = Poly3DCollection(
            polygons[visible], facecolors=cmap(norm(values[visible])),
            edgecolors=(0.12, 0.12, 0.12, 0.22), linewidths=0.18,
            antialiased=True,
        )
        axis.add_collection3d(colored)
    # Only historical nodal quads can be mixed. New triangles carry one
    # authoritative face classification and therefore never enter this path.
    if not accepted_only and np.any(boundary):
        transition = Poly3DCollection(
            polygons[boundary], facecolors=(0.55, 0.55, 0.55, 0.40),
            edgecolors=(0.2, 0.2, 0.2, 0.25), linewidths=0.18,
        )
        axis.add_collection3d(transition)
    if show_sun:
        _draw_sun(axis, scale)
    if show_distance_axis:
        _draw_distance_axis(axis, surface, scale, distance_axis_maximum_rs)

    axis.set_xlim(center[0] - half_span, center[0] + half_span)
    axis.set_ylim(center[1] - half_span, center[1] + half_span)
    axis.set_zlim(center[2] - half_span, center[2] + half_span)
    axis.set_xlabel(f"HCI x [{unit_label}]")
    axis.set_ylabel(f"HCI y [{unit_label}]")
    axis.set_zlabel(f"HCI z [{unit_label}]")
    axis.view_init(elev=elevation, azim=azimuth)
    mode = "accepted shocks only" if accepted_only else \
        "complete front (mixed boundary gray)"
    axis.set_title(
        f"{surface.path.name}\n"
        f"t={surface_time(surface):g} s; generation={surface_generation(surface)}; "
        f"{variable}; {mode}; shown faces={int(np.count_nonzero(visible))}"
    )
    return (int(np.count_nonzero(accepted)), int(np.count_nonzero(nonshock)),
            int(np.count_nonzero(boundary)), int(np.count_nonzero(visible)))


def render(surface: FrontSurface, variable: str, accepted_only: bool,
           length_unit: str, cmap_name: str, vmin: Optional[float],
           vmax: Optional[float], elevation: float, azimuth: float,
           save: Optional[Path], dpi: int, show: bool, show_sun: bool,
           show_distance_axis: bool) -> None:
    """Render the HCI surface and install an interactive acceptance filter."""

    # Select a noninteractive backend before importing pyplot.  This supports
    # batch rendering on compute nodes without changing interactive behavior on
    # workstations that have a display server.
    import matplotlib
    if not show or not os.environ.get("DISPLAY"):
        matplotlib.use("Agg")
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    from matplotlib.widgets import CheckButtons
    # Matplotlib 3.1 registers the ``3d`` projection as an import side effect;
    # newer releases do it eagerly.  Retain this explicit compatibility import
    # because AMPS sites can carry the older module.
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    scale, unit_label = _length_scale(length_unit)
    normalizing_values = _normalization_values((surface,), variable)
    lower, upper = _color_limits(normalizing_values, vmin, vmax)
    norm = colors.Normalize(vmin=lower, vmax=upper)
    cmap = cm.get_cmap(cmap_name)
    center, half_span = _global_bounds(
        (surface,), scale, show_sun or show_distance_axis
    )

    figure = plt.figure(figsize=(11.5, 8.0))
    axis = figure.add_subplot(111, projection="3d")
    figure.subplots_adjust(left=0.04, right=0.86, bottom=0.08, top=0.92)
    scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_map.set_array(normalizing_values)
    colorbar = figure.colorbar(scalar_map, ax=axis, fraction=0.035, pad=0.08)
    colorbar.set_label(variable)
    state = {"accepted_only": bool(accepted_only),
             "elevation": float(elevation), "azimuth": float(azimuth)}

    def draw_surface() -> None:
        _draw_frame(
            axis, surface, variable, state["accepted_only"], scale, unit_label,
            cmap, norm, center, half_span, state["elevation"], state["azimuth"],
            show_sun, show_distance_axis, None,
        )
        figure.canvas.draw_idle()

    control_axis = figure.add_axes([0.865, 0.82, 0.13, 0.08])
    checkbox = CheckButtons(
        control_axis, ["accepted\nshock only"], [state["accepted_only"]]
    )

    def toggle(_label: str) -> None:
        state["elevation"] = float(axis.elev)
        state["azimuth"] = float(axis.azim)
        state["accepted_only"] = not state["accepted_only"]
        draw_surface()

    checkbox.on_clicked(toggle)
    draw_surface()

    def report_view(event) -> None:
        if event.key and event.key.lower() == "v":
            print(f"selected_view=--elev {axis.elev:.8g} --azim {axis.azim:.8g}")

    figure.canvas.mpl_connect("key_press_event", report_view)

    accepted, nonshock, boundary = classify_faces(surface)
    accepted_count = int(np.count_nonzero(accepted))
    nonshock_count = int(np.count_nonzero(nonshock))
    boundary_count = int(np.count_nonzero(boundary))
    print(
        f"nodes={surface.data.shape[0]} faces={surface.faces.shape[0]} "
        f"zone_type={surface.zone_type} accepted_faces={accepted_count} "
        f"nonshock_faces={nonshock_count} mixed_boundary_faces={boundary_count} "
        f"variable={variable} "
        f"range=[{lower:.17g},{upper:.17g}]"
    )
    if save is not None:
        save = save.expanduser().resolve()
        save.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(save, dpi=dpi, bbox_inches="tight")
        print(f"saved={save}")
    if show and os.environ.get("DISPLAY"):
        plt.show()
    elif show and not os.environ.get("DISPLAY"):
        print("DISPLAY is unset; rendered without opening a window", file=sys.stderr)
    plt.close(figure)


def render_movie(surfaces: Sequence[FrontSurface], variable: str,
                 accepted_only: bool, length_unit: str, cmap_name: str,
                 vmin: Optional[float], vmax: Optional[float],
                 elevation: float, azimuth: float, destination: Path,
                 fps: float, dpi: int, show_sun: bool,
                 show_distance_axis: bool) -> None:
    """Encode a time-ordered series with one immutable view and color scale.

    The global bounding cube is computed from every frame (and the origin when
    reference geometry is requested).  Thus neither camera orientation nor
    zoom changes as the front propagates.  A single color normalization is also
    frozen before encoding so apparent temporal changes are not colorbar
    rescaling artifacts.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.animation as animation
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    ordered = sort_movie_surfaces(surfaces)
    if len(ordered) < 2:
        raise FrontFormatError("--movie requires at least two distinct epochs")
    if not np.isfinite(fps) or fps <= 0.0:
        raise FrontFormatError("--fps must be finite and positive")
    destination = destination.expanduser().resolve()
    suffix = destination.suffix.lower()
    if suffix == ".gif":
        if not animation.writers.is_available("pillow"):
            raise FrontFormatError("GIF output requires the Matplotlib Pillow writer")
    elif suffix == ".mp4":
        if not animation.writers.is_available("ffmpeg"):
            raise FrontFormatError(
                "MP4 output requires ffmpeg; use a .gif destination on this system"
            )
    else:
        raise FrontFormatError("--movie destination must end in .gif or .mp4")
    scale, unit_label = _length_scale(length_unit)
    normalizing_values = _normalization_values(ordered, variable)
    lower, upper = _color_limits(normalizing_values, vmin, vmax)
    norm = colors.Normalize(vmin=lower, vmax=upper)
    cmap = cm.get_cmap(cmap_name)
    center, half_span = _global_bounds(
        ordered, scale, show_sun or show_distance_axis
    )
    # A movie uses one ruler extent and tick spacing determined from the whole
    # series.  Its direction still follows each frame's inferred front axis, so
    # it passes through that frame's middle even for a time-dependent event
    # direction.  The fixed extent prevents low-coronal frames from squeezing
    # many sub-R_sun labels into pixels set by the final 1-AU bounding cube.
    distance_axis_maximum_rs = None
    if show_distance_axis:
        distance_axis_maximum_rs = max(
            float(np.max(_coordinates(item, 1.0) @ _shock_axis_direction(item)) /
                  SOLAR_RADIUS_M) for item in ordered
        )

    figure = plt.figure(figsize=(11.5, 8.0))
    axis = figure.add_subplot(111, projection="3d")
    figure.subplots_adjust(left=0.04, right=0.88, bottom=0.07, top=0.93)
    scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_map.set_array(normalizing_values)
    colorbar = figure.colorbar(scalar_map, ax=axis, fraction=0.035, pad=0.08)
    colorbar.set_label(variable)

    def draw(index: int):
        _draw_frame(
            axis, ordered[index], variable, accepted_only, scale, unit_label,
            cmap, norm, center, half_span, elevation, azimuth, show_sun,
            show_distance_axis, distance_axis_maximum_rs,
        )
        return tuple(axis.collections) + tuple(axis.lines) + tuple(axis.texts)

    movie = animation.FuncAnimation(
        figure, draw, frames=len(ordered), interval=1000.0 / fps,
        repeat=True, blit=False,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if suffix == ".gif":
        writer = animation.PillowWriter(fps=fps)
    else:
        writer = animation.FFMpegWriter(fps=fps, metadata={
            "title": "srcSEP3D reduced shock front",
            "comment": "fixed HCI camera, bounds, and color normalization",
        })
    movie.save(str(destination), writer=writer, dpi=dpi)
    plt.close(figure)
    times = [surface_time(item) for item in ordered]
    print(
        f"movie_frames={len(ordered)} time_range_s=[{times[0]:.17g},"
        f"{times[-1]:.17g}] fixed_view_elev={elevation:.17g} "
        f"fixed_view_azim={azimuth:.17g} fixed_half_span={half_span:.17g} "
        f"color_range=[{lower:.17g},{upper:.17g}] saved={destination}"
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "View or animate srcSEP3D reduced-front FETRIANGLE files (and "
            "historical FEQUADRILATERAL files).\n\nThe production surface uses "
            "nodal geometry and authoritative face-centred shock variables. "
            "Legacy nodal values are averaged only over quads whose nodes "
            "agree on shock acceptance; mixed legacy boundaries are gray "
            "or hidden."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
Find the available front snapshots in a native-output directory:
  find OUTPUT_DIRECTORY -maxdepth 1 -name '*-front.dat' -print | sort

Inspect the variables stored in a front file:
  python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat --list-variables

Open an interactive fast-Mach view:
  python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat --variable fast_mach

Plot density compression, hide non-shock areas, and show the Sun/ruler:
  python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat \\
    --variable density_compression --accepted-only \\
    --show-sun --show-distance-axis

Choose a colormap and a fixed physical color interval:
  python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat \\
    --variable theta_Bn_rad --cmap plasma --vmin 0 --vmax 1.5707963268

List every colormap installed in the current Matplotlib environment:
  python3 -c "import matplotlib.pyplot as p; print(*p.colormaps(), sep=chr(10))"

View a signed IMF component from a fixed HCI camera direction, with axes in AU:
  python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat \\
    --variable b1z_T --cmap seismic --elev 20 --azim 135 --length-unit au

Save a headless still image on a compute node:
  python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat \\
    --variable magnetic_compression --accepted-only \\
    --save front-magnetic-compression.png --no-show

Make a GIF from every *-front.dat file in one output directory:
  python3 srcSEP3D/examples/shock-front/view_front.py OUTPUT_DIRECTORY \\
    --variable fast_mach --accepted-only \\
    --show-sun --show-distance-axis \\
    --elev 24 --azim -58 --fps 3 --movie front.gif

Make an MP4 with a reproducible view, physical range, and rendering quality:
  python3 srcSEP3D/examples/shock-front/view_front.py OUTPUT_DIRECTORY \\
    --variable density_compression --accepted-only --cmap plasma \\
    --vmin 1 --vmax 4 --elev 20 --azim 135 \\
    --show-sun --show-distance-axis --fps 12 --dpi 170 --movie front.mp4

The movie input may instead be a quoted pattern or an explicit file list:
  python3 srcSEP3D/examples/shock-front/view_front.py \\
    'products/*-front.dat' --movie front.gif
  python3 srcSEP3D/examples/shock-front/view_front.py \\
    front-0001.dat front-0002.dat --movie front.gif

Selecting a movie viewpoint
----------------------------
First open one representative file interactively, rotate/zoom the 3-D view,
and press the `v` key. The terminal prints reusable camera arguments:
  selected_view=--elev 24.5 --azim -61
Copy those values into the movie command. The selected elevation, azimuth,
global HCI bounds, ruler extent, and color normalization are then identical in
every frame. Use --vmin/--vmax when movies must share a prescribed color scale.
Elevation is measured above the HCI x-y plane; azimuth rotates about HCI +z.
These camera options change only the view, never the stored geometry or normals.

Variables and colormaps
-----------------------
--variable must exactly match a name printed by --list-variables, for example
fast_mach, density_compression, magnetic_compression, theta_Bn_rad, rho1_kg_m3,
p1_Pa, u1x_m_s, or b1z_T. Vector quantities are stored as separate x/y/z
components. Use shock_accepted/downstream_valid to interpret RH-limit fields.

--cmap accepts an installed Matplotlib colormap name. Portable examples include
viridis (default), plasma, inferno, magma, cividis, coolwarm, and seismic.
Sequential maps suit positive magnitudes; diverging maps such as coolwarm or
seismic suit signed vector components. --vmin and --vmax retain the same
physical color meaning across snapshots; omitting them uses the range over all
input frames. Invalid colormap names and inverted color limits fail explicitly.

Movie format and physical meaning
---------------------------------
A .gif destination uses Pillow. A .mp4 destination requires ffmpeg. Inputs are
ordered by embedded time_s, not filename. The viewer displays a prescribed
front and front-local RH limits; it does not create sheath/ejecta volume data.
""",
    )
    parser.add_argument(
        "front", nargs="+",
        help="front file(s), a directory, or a quoted *-front.dat pattern",
    )
    parser.add_argument("--variable", default="fast_mach",
                        help=("exact node-variable name used for color; run "
                              "--list-variables to discover valid names"))
    parser.add_argument("--list-variables", action="store_true",
                        help="print variables in the file and exit")
    parser.add_argument("--accepted-only", action="store_true",
                        help="initially hide every face not classified as an accepted shock")
    parser.add_argument("--show-sun", action="store_true",
                        help="draw the one-solar-radius Sun at the HCI origin")
    parser.add_argument(
        "--show-distance-axis", action="store_true",
        help="draw a Sun-centered ruler through the front middle, labeled in R_sun",
    )
    parser.add_argument("--length-unit", choices=("m", "rs", "au"), default="rs",
                        help="HCI axis unit (default: rs, solar radii)")
    parser.add_argument(
        "--cmap", default="viridis",
        help="Matplotlib colormap (e.g. viridis, plasma, coolwarm; default: viridis)",
    )
    parser.add_argument("--vmin", type=float,
                        help=("fixed physical color-scale minimum shared by all "
                              "movie frames (default: minimum over all inputs)"))
    parser.add_argument("--vmax", type=float,
                        help=("fixed physical color-scale maximum shared by all "
                              "movie frames (default: maximum over all inputs)"))
    parser.add_argument("--elev", type=float, default=24.0,
                        help=("camera elevation above the HCI x-y plane in degrees; "
                              "frozen for every movie frame"))
    parser.add_argument("--azim", type=float, default=-58.0,
                        help=("camera azimuth about HCI +z in degrees; frozen for "
                              "every movie frame"))
    parser.add_argument("--save", type=Path, help="also save a PNG/PDF/SVG image")
    parser.add_argument(
        "--movie", type=Path,
        help="encode all input epochs to .gif, or to .mp4 when ffmpeg is installed",
    )
    parser.add_argument("--fps", type=float, default=5.0,
                        help="movie frames per second (default: 5)")
    parser.add_argument("--dpi", type=int, default=170, help="saved-image DPI")
    parser.add_argument("--no-show", action="store_true",
                        help="render/save without opening an interactive window")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        paths = expand_front_paths(args.front)
        surfaces = tuple(load_front(path) for path in paths)
        if args.list_variables:
            print("\n".join(surfaces[0].variables))
            return 0
        if args.movie is not None:
            if args.save is not None:
                raise FrontFormatError("--save and --movie are mutually exclusive")
            render_movie(
                surfaces=surfaces,
                variable=args.variable,
                accepted_only=args.accepted_only,
                length_unit=args.length_unit,
                cmap_name=args.cmap,
                vmin=args.vmin,
                vmax=args.vmax,
                elevation=args.elev,
                azimuth=args.azim,
                destination=args.movie,
                fps=args.fps,
                dpi=args.dpi,
                show_sun=args.show_sun,
                show_distance_axis=args.show_distance_axis,
            )
            return 0
        if len(surfaces) != 1:
            raise FrontFormatError("multiple front files require --movie")
        render(
            surface=surfaces[0],
            variable=args.variable,
            accepted_only=args.accepted_only,
            length_unit=args.length_unit,
            cmap_name=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            elevation=args.elev,
            azimuth=args.azim,
            save=args.save,
            dpi=args.dpi,
            show=not args.no_show,
            show_sun=args.show_sun,
            show_distance_axis=args.show_distance_axis,
        )
    except (FrontFormatError, OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())

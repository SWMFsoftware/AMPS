#!/usr/bin/env python3
"""Visualize an AMPS lunar surface without changing scientific artifacts.

The utility reads the production CEA or Tecplot FETRIANGLE representation
written by ``Moon::Surface::WriteCeaSurface`` and
``Moon::Surface::WriteTecplotSurface``.  It does not reconstruct, smooth, or
otherwise alter the geometry.  Figures are visualization products, not
verification or validation evidence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
from typing import TextIO

import numpy as np


NOTICE = "VISUALIZATION ONLY - NOT V&V EVIDENCE"


@dataclass(frozen=True)
class SurfaceMesh:
    """A validated, centre-relative triangular surface in SI units."""

    vertices_m: np.ndarray
    faces: np.ndarray
    source_format: str
    stored_radius_m: np.ndarray | None = None


def sha256(path: Path) -> str:
    """Return a streaming source hash for the printed visualization receipt."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def next_nonempty(stream: TextIO, context: str) -> str:
    """Read one required non-empty record and report a precise truncation."""
    while True:
        line = stream.readline()
        if line == "":
            raise ValueError(f"unexpected end of file while reading {context}")
        stripped = line.strip()
        if stripped:
            return stripped


def parse_floats(line: str, count: int, context: str) -> list[float]:
    """Parse a fixed-width numeric record; extra columns are not discarded."""
    fields = line.split()
    if len(fields) != count:
        raise ValueError(f"{context} requires {count} columns, found {len(fields)}")
    try:
        values = [float(field) for field in fields]
    except ValueError as exception:
        raise ValueError(f"{context} contains a non-numeric value") from exception
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{context} contains a non-finite value")
    return values


def parse_indices(line: str, context: str) -> list[int]:
    """Read exactly three one-based integer node identifiers."""
    fields = line.split()
    if len(fields) != 3:
        raise ValueError(f"{context} requires three node indices")
    try:
        values = [int(field) for field in fields]
    except ValueError as exception:
        raise ValueError(f"{context} contains a non-integer node index") from exception
    return values


def require_no_extra_records(stream: TextIO) -> None:
    """Reject appended zones or partial records instead of silently hiding them."""
    for line in stream:
        if line.strip():
            raise ValueError("surface file contains records after declared connectivity")


def read_cea(path: Path) -> SurfaceMesh:
    """Read the exact one-based CEA long format emitted by srcMoon."""
    with path.open("r", encoding="utf-8") as stream:
        counts = next_nonempty(stream, "CEA node/face counts").split()
        if len(counts) != 2:
            raise ValueError("CEA first record must contain node and face counts")
        try:
            vertex_count, face_count = (int(value) for value in counts)
        except ValueError as exception:
            raise ValueError("CEA counts must be integers") from exception
        if vertex_count < 3 or face_count < 1:
            raise ValueError("CEA surface must contain at least 3 nodes and 1 face")

        vertices = np.empty((vertex_count, 3), dtype=np.float64)
        for index in range(vertex_count):
            vertices[index] = parse_floats(
                next_nonempty(stream, f"CEA vertex {index + 1}"), 3,
                f"CEA vertex {index + 1}")
        faces = np.empty((face_count, 3), dtype=np.int64)
        for index in range(face_count):
            one_based = parse_indices(
                next_nonempty(stream, f"CEA face {index + 1}"),
                f"CEA face {index + 1}")
            faces[index] = np.asarray(one_based, dtype=np.int64) - 1
        require_no_extra_records(stream)
    return SurfaceMesh(vertices, faces, "CEA")


def read_tecplot(path: Path) -> SurfaceMesh:
    """Read srcMoon's point-packed Tecplot FETRIANGLE inspection format."""
    with path.open("r", encoding="utf-8") as stream:
        variables = next_nonempty(stream, "Tecplot VARIABLES header")
        if not variables.upper().startswith("VARIABLES="):
            raise ValueError("Tecplot surface is missing its VARIABLES header")
        for required in ('"X [m]"', '"Y [m]"', '"Z [m]"', '"Radius [m]"'):
            if required.lower() not in variables.lower():
                raise ValueError(f"Tecplot VARIABLES header is missing {required}")

        zone = next_nonempty(stream, "Tecplot ZONE header")
        match_n = re.search(r"\bN\s*=\s*(\d+)", zone, re.IGNORECASE)
        match_e = re.search(r"\bE\s*=\s*(\d+)", zone, re.IGNORECASE)
        if match_n is None or match_e is None:
            raise ValueError("Tecplot ZONE must declare N and E")
        if "DATAPACKING=POINT" not in zone.upper().replace(" ", ""):
            raise ValueError("Tecplot surface must use DATAPACKING=POINT")
        if "ZONETYPE=FETRIANGLE" not in zone.upper().replace(" ", ""):
            raise ValueError("Tecplot surface must use ZONETYPE=FETRIANGLE")
        vertex_count = int(match_n.group(1))
        face_count = int(match_e.group(1))
        if vertex_count < 3 or face_count < 1:
            raise ValueError("Tecplot surface must contain at least 3 nodes and 1 face")

        vertices = np.empty((vertex_count, 3), dtype=np.float64)
        stored_radius = np.empty(vertex_count, dtype=np.float64)
        for index in range(vertex_count):
            values = parse_floats(
                next_nonempty(stream, f"Tecplot vertex {index + 1}"), 4,
                f"Tecplot vertex {index + 1}")
            vertices[index] = values[:3]
            stored_radius[index] = values[3]
        faces = np.empty((face_count, 3), dtype=np.int64)
        for index in range(face_count):
            one_based = parse_indices(
                next_nonempty(stream, f"Tecplot face {index + 1}"),
                f"Tecplot face {index + 1}")
            faces[index] = np.asarray(one_based, dtype=np.int64) - 1
        require_no_extra_records(stream)
    return SurfaceMesh(vertices, faces, "TECPLOT_FETRIANGLE", stored_radius)


def read_surface(path: Path, requested_format: str) -> SurfaceMesh:
    """Select a parser by explicit option, suffix, or the first file record."""
    if requested_format == "cea":
        return read_cea(path)
    if requested_format == "tecplot":
        return read_tecplot(path)
    suffix = path.suffix.lower()
    if suffix == ".cea":
        return read_cea(path)
    if suffix in {".dat", ".plt", ".tec"}:
        return read_tecplot(path)
    with path.open("r", encoding="utf-8") as stream:
        first = next_nonempty(stream, "surface header")
    return read_tecplot(path) if first.upper().startswith("VARIABLES=") else read_cea(path)


def validate_surface(mesh: SurfaceMesh) -> dict[str, float | int]:
    """Validate production geometry before allowing a potentially misleading plot."""
    vertices = mesh.vertices_m
    faces = mesh.faces
    if vertices.ndim != 2 or vertices.shape[1] != 3 or not np.isfinite(vertices).all():
        raise ValueError("surface vertices are not a finite N-by-3 array")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("surface faces are not an E-by-3 array")
    if np.any(faces < 0) or np.any(faces >= len(vertices)):
        raise ValueError("surface connectivity references a node outside [1,N]")
    if np.any(faces[:, 0] == faces[:, 1]) or np.any(faces[:, 1] == faces[:, 2]) or np.any(faces[:, 0] == faces[:, 2]):
        raise ValueError("surface contains a face with repeated node indices")

    radius = np.linalg.norm(vertices, axis=1)
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("surface contains a non-positive or non-finite radius")
    if mesh.stored_radius_m is not None:
        # The tolerance only covers decimal roundoff in the text writer.  A
        # larger disagreement means the displayed scalar is not the Cartesian
        # geometry and must not be plotted as though it were.
        tolerance_m = np.maximum(1.0e-6, 1.0e-10 * radius)
        if np.any(np.abs(mesh.stored_radius_m - radius) > tolerance_m):
            raise ValueError("Tecplot Radius column disagrees with Cartesian nodes")

    degenerate = 0
    inward = 0
    # Process topology in bounded chunks: fine LOLA surfaces can contain
    # millions of triangles, and validation should not triple their memory
    # footprint merely to create a picture.
    for start in range(0, len(faces), 200_000):
        triangle = vertices[faces[start:start + 200_000]]
        normal = np.cross(triangle[:, 1] - triangle[:, 0],
                          triangle[:, 2] - triangle[:, 0])
        area2 = np.linalg.norm(normal, axis=1)
        centre = triangle.sum(axis=1)
        degenerate += int(np.count_nonzero(~np.isfinite(area2) | (area2 <= 0.0)))
        inward += int(np.count_nonzero(np.einsum("ij,ij->i", normal, centre) <= 0.0))
    if degenerate:
        raise ValueError(f"surface contains {degenerate} degenerate faces")
    if inward:
        raise ValueError(f"surface contains {inward} inward-facing faces")
    return {
        "vertex_count": int(len(vertices)),
        "face_count": int(len(faces)),
        "minimum_radius_m": float(radius.min()),
        "maximum_radius_m": float(radius.max()),
        "degenerate_face_count": degenerate,
        "inward_face_count": inward,
    }


def selected_face_indices(face_count: int, limit: int) -> np.ndarray:
    """Return all faces or a deterministic, explicitly requested subset."""
    if limit == 0 or limit >= face_count:
        return np.arange(face_count, dtype=np.int64)
    if limit < 1:
        raise ValueError("--max-render-faces must be zero or positive")
    # Linspace samples the stored connectivity ordering reproducibly.  The
    # figure subtitle discloses the subset; all validation statistics above
    # still use the complete mesh.
    return np.linspace(0, face_count - 1, limit, dtype=np.int64)


def surface_scalar(mesh: SurfaceMesh, color_by: str,
                   reference_radius_m: float | None) -> tuple[np.ndarray, str]:
    """Return the requested vertex color without inventing a lunar datum."""
    radius_m = np.linalg.norm(mesh.vertices_m, axis=1)
    if color_by == "elevation-m":
        if (reference_radius_m is None or
                not math.isfinite(reference_radius_m) or
                reference_radius_m <= 0.0):
            raise ValueError(
                "--color-by elevation-m requires a positive "
                "--reference-radius-m")
        return (radius_m - reference_radius_m,
                f"elevation above {reference_radius_m:.6g} m sphere [m]")
    return radius_m / 1000.0, "centre-relative radius [km]"


def display_vertices(mesh: SurfaceMesh, reference_radius_m: float | None,
                     vertical_exaggeration: float) -> np.ndarray:
    """Return display-only vertices with explicitly exaggerated radial relief.

    LOLA elevation is a radial offset from its reference sphere.  Multiplying
    that offset, while retaining every vertex direction, makes kilometre-scale
    relief visible on a 1737-km globe.  The source ``mesh.vertices_m`` array is
    never modified and all validation/statistical values remain unexaggerated.
    """
    if (not math.isfinite(vertical_exaggeration) or
            vertical_exaggeration <= 0.0):
        raise ValueError("--vertical-exaggeration must be finite and positive")
    if vertical_exaggeration == 1.0:
        return mesh.vertices_m
    if (reference_radius_m is None or
            not math.isfinite(reference_radius_m) or
            reference_radius_m <= 0.0):
        raise ValueError(
            "--vertical-exaggeration other than 1 requires a positive "
            "--reference-radius-m")
    radius_m = np.linalg.norm(mesh.vertices_m, axis=1)
    displayed_radius_m = (reference_radius_m + vertical_exaggeration *
                          (radius_m - reference_radius_m))
    if np.any(displayed_radius_m <= 0.0):
        raise ValueError(
            "vertical exaggeration produces a non-positive display radius")
    return mesh.vertices_m * (displayed_radius_m / radius_m)[:, np.newaxis]


def exaggeration_notice(vertical_exaggeration: float) -> str:
    """Return an unavoidable figure label for display-only relief scaling."""
    if vertical_exaggeration == 1.0:
        return "true radial geometry (vertical exaggeration x1)"
    return (f"VERTICAL EXAGGERATION x{vertical_exaggeration:g} - "
            "DISPLAY GEOMETRY ONLY")


def configure_surface_axis(axis: object, triangle_km: np.ndarray,
                           face_scalar: np.ndarray, normalizer: object,
                           colormap: object, extent_km: float,
                           show_edges: bool, elevation_deg: float,
                           azimuth_deg: float, axis_title: str) -> None:
    """Populate one 3-D axis consistently for static and interactive views."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    edge_color = (0.08, 0.08, 0.08, 1.0) if show_edges else "none"
    collection = Poly3DCollection(
        triangle_km, linewidths=0.12 if show_edges else 0.0,
        edgecolors=edge_color)
    collection.set_facecolor(colormap(normalizer(face_scalar)))
    collection.set_rasterized(True)
    axis.add_collection3d(collection)
    axis.set_xlim(-extent_km, extent_km)
    axis.set_ylim(-extent_km, extent_km)
    axis.set_zlim(-extent_km, extent_km)
    # Matplotlib before 3.3 lacks set_box_aspect(). Equal numerical limits
    # still prevent axis-dependent data scaling there; newer versions also
    # receive an explicitly cubic display box.
    if hasattr(axis, "set_box_aspect"):
        axis.set_box_aspect((1.0, 1.0, 1.0))
    # Old PostScript backends warn and discard translucent 3-D panes. Opaque
    # panes render consistently in PNG, EPS, PDF, and SVG.
    for pane_axis in (axis.xaxis, axis.yaxis, axis.zaxis):
        pane_axis.pane.set_alpha(1.0)
    axis.set_xlabel("X [km]")
    axis.set_ylabel("Y [km]")
    axis.set_zlabel("Z [km]")
    axis.view_init(elev=elevation_deg, azim=azimuth_deg)
    axis.set_title(axis_title)


def render_surface(mesh: SurfaceMesh, source: Path, outputs: list[Path],
                   color_by: str, reference_radius_m: float | None,
                   frame_label: str, title: str | None, view_latitude_deg: float,
                   view_longitude_deg: float, max_render_faces: int,
                   show_edges: bool, dpi: int, show: bool,
                   gui_backend: bool = False,
                   vertical_exaggeration: float = 1.0) -> dict[str, object]:
    """Render two globe views, a body-fixed map, and the radial distribution."""
    if not show and not gui_backend:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    # Matplotlib 3.1 registers projection="3d" only when mplot3d is imported.
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    scalar, scalar_label = surface_scalar(mesh, color_by, reference_radius_m)
    displayed_vertices_m = display_vertices(
        mesh, reference_radius_m, vertical_exaggeration)

    indices = selected_face_indices(len(mesh.faces), max_render_faces)
    faces = mesh.faces[indices]
    actual_triangle_km = mesh.vertices_m[faces] / 1000.0
    triangle_km = displayed_vertices_m[faces] / 1000.0
    face_scalar = scalar[faces].mean(axis=1)
    scalar_min = float(scalar.min())
    scalar_max = float(scalar.max())
    if scalar_min == scalar_max:
        padding = max(1.0, abs(scalar_min) * 1.0e-12)
        normalizer = colors.Normalize(scalar_min - padding, scalar_max + padding)
    else:
        normalizer = colors.Normalize(scalar_min, scalar_max)
    colormap = cm.get_cmap("terrain")

    # A near-square cell is important with Matplotlib 3.1, whose 3-D axes do
    # not yet expose set_box_aspect().  An explicit GridSpec also reserves
    # enough inter-column space for the map colorbar and histogram labels.
    figure = plt.figure(figsize=(14.0, 12.0))
    grid = figure.add_gridspec(
        2, 2, left=0.055, right=0.965, bottom=0.065, top=0.885,
        wspace=0.34, hspace=0.24)
    axes_3d = [figure.add_subplot(grid[0, 0], projection="3d"),
               figure.add_subplot(grid[0, 1], projection="3d")]
    extent_km = float(np.max(np.abs(displayed_vertices_m)) / 1000.0)
    for index, axis in enumerate(axes_3d):
        configure_surface_axis(
            axis, triangle_km, face_scalar, normalizer, colormap, extent_km,
            show_edges, view_latitude_deg,
            view_longitude_deg + 180.0 * index,
            "requested view" if index == 0 else "antipodal view")

    # Use triangle-centre directions in the body-fixed Cartesian frame.  This
    # is a display projection only; vertices and topology remain unchanged.
    # The map uses the true production coordinates. Vertical exaggeration is
    # confined to 3-D display geometry and must not shift mapped directions.
    centres = actual_triangle_km.mean(axis=1)
    centre_radius = np.linalg.norm(centres, axis=1)
    latitude = np.degrees(np.arcsin(centres[:, 2] / centre_radius))
    longitude = np.degrees(np.arctan2(centres[:, 1], centres[:, 0]))
    map_axis = figure.add_subplot(grid[1, 0])
    map_plot = map_axis.scatter(longitude, latitude, c=face_scalar, s=4.0,
                                cmap=colormap, norm=normalizer,
                                linewidths=0, rasterized=True)
    map_axis.set_xlim(-180.0, 180.0)
    map_axis.set_ylim(-90.0, 90.0)
    map_axis.set_xlabel("east longitude [deg]")
    map_axis.set_ylabel("planetocentric latitude [deg]")
    map_axis.set_title("body-fixed triangle-centre map")
    map_axis.grid(True, linewidth=0.35, color="#b0b0b0")
    figure.colorbar(map_plot, ax=map_axis, fraction=0.046, pad=0.025,
                    label=scalar_label)

    histogram = figure.add_subplot(grid[1, 1])
    histogram.hist(scalar, bins=80, color="tab:blue", edgecolor="none")
    histogram.set_xlabel(scalar_label)
    histogram.set_ylabel("vertex count")
    histogram.set_title("complete vertex distribution")
    histogram.grid(True, linewidth=0.35, color="#b0b0b0")

    rendered = len(indices)
    subset_notice = ("all faces rendered" if rendered == len(mesh.faces) else
                     f"{rendered:,} of {len(mesh.faces):,} faces rendered; deterministic display subset")
    figure.suptitle(
        (title or f"AMPS lunar surface: {source.name}") + "\n" +
        f"{frame_label}; {subset_notice}\n" +
        f"{exaggeration_notice(vertical_exaggeration)}\n{NOTICE}", fontsize=13)
    saved: list[dict[str, str]] = []
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.suffix.lower() not in {".png", ".eps", ".pdf", ".svg"}:
            raise ValueError(f"unsupported figure suffix for {output}")
        # Backend metadata dictionaries differ among PNG, EPS, PDF, and SVG.
        # Keep provenance in the uniform JSON receipt rather than silently
        # dropping fields or failing only for one requested output format.
        figure.savefig(output, dpi=dpi)
        saved.append({"path": str(output.resolve()), "sha256": sha256(output)})
    if show:
        plt.show()
    plt.close(figure)
    return {
        "color_by": color_by,
        "reference_radius_m": reference_radius_m,
        "frame_label": frame_label,
        "rendered_face_count": rendered,
        "complete_face_count": int(len(mesh.faces)),
        "deterministic_subset": rendered != len(mesh.faces),
        "vertical_exaggeration": vertical_exaggeration,
        "display_geometry_is_exaggerated": vertical_exaggeration != 1.0,
        "display_radius_formula": (
            "reference_radius_m + vertical_exaggeration * "
            "(true_radius_m - reference_radius_m)"),
        "outputs": saved,
    }


def save_selected_view(mesh: SurfaceMesh, source: Path, outputs: list[Path],
                       color_by: str, reference_radius_m: float | None,
                       frame_label: str, title: str | None,
                       elevation_deg: float, azimuth_deg: float,
                       max_render_faces: int, show_edges: bool,
                       dpi: int,
                       vertical_exaggeration: float = 1.0) -> list[dict[str, str]]:
    """Save a clean single-panel image at the camera selected by the user."""
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    scalar, scalar_label = surface_scalar(mesh, color_by, reference_radius_m)
    displayed_vertices_m = display_vertices(
        mesh, reference_radius_m, vertical_exaggeration)
    indices = selected_face_indices(len(mesh.faces), max_render_faces)
    faces = mesh.faces[indices]
    triangle_km = displayed_vertices_m[faces] / 1000.0
    face_scalar = scalar[faces].mean(axis=1)
    scalar_min = float(scalar.min())
    scalar_max = float(scalar.max())
    if scalar_min == scalar_max:
        padding = max(1.0, abs(scalar_min) * 1.0e-12)
        normalizer = colors.Normalize(scalar_min - padding, scalar_max + padding)
    else:
        normalizer = colors.Normalize(scalar_min, scalar_max)
    colormap = cm.get_cmap("terrain")
    extent_km = float(np.max(np.abs(displayed_vertices_m)) / 1000.0)

    figure = plt.figure(figsize=(9.0, 9.0))
    axis = figure.add_subplot(111, projection="3d")
    configure_surface_axis(
        axis, triangle_km, face_scalar, normalizer, colormap, extent_km,
        show_edges, elevation_deg, azimuth_deg,
        f"elevation={elevation_deg:.2f} deg, azimuth={azimuth_deg:.2f} deg")
    scalar_map = cm.ScalarMappable(norm=normalizer, cmap=colormap)
    scalar_map.set_array(np.asarray([], dtype=np.float64))
    figure.colorbar(scalar_map, ax=axis, fraction=0.035, pad=0.08,
                    label=scalar_label)
    subset_notice = (
        "all faces rendered" if len(indices) == len(mesh.faces) else
        f"{len(indices):,} of {len(mesh.faces):,} faces rendered; "
        "deterministic display subset")
    figure.suptitle(
        (title or f"AMPS lunar surface: {source.name}") + "\n" +
        f"{frame_label}; {subset_notice}\n" +
        f"{exaggeration_notice(vertical_exaggeration)}\n{NOTICE}", fontsize=12)
    figure.subplots_adjust(left=0.02, right=0.91, bottom=0.03, top=0.88)
    saved: list[dict[str, str]] = []
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=dpi)
        saved.append({"path": str(output.resolve()), "sha256": sha256(output)})
    plt.close(figure)
    return saved


def run_interactive_3d(mesh: SurfaceMesh, source: Path,
                       selected_outputs: list[Path], color_by: str,
                       reference_radius_m: float | None, frame_label: str,
                       title: str | None, initial_elevation_deg: float,
                       initial_azimuth_deg: float, max_render_faces: int,
                       show_edges: bool, dpi: int,
                       vertical_exaggeration: float = 1.0) -> dict[str, object]:
    """Open a rotatable globe and save the investigator-selected camera."""
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from matplotlib.widgets import Button
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    scalar, scalar_label = surface_scalar(mesh, color_by, reference_radius_m)
    displayed_vertices_m = display_vertices(
        mesh, reference_radius_m, vertical_exaggeration)
    indices = selected_face_indices(len(mesh.faces), max_render_faces)
    faces = mesh.faces[indices]
    triangle_km = displayed_vertices_m[faces] / 1000.0
    face_scalar = scalar[faces].mean(axis=1)
    scalar_min = float(scalar.min())
    scalar_max = float(scalar.max())
    if scalar_min == scalar_max:
        padding = max(1.0, abs(scalar_min) * 1.0e-12)
        normalizer = colors.Normalize(scalar_min - padding, scalar_max + padding)
    else:
        normalizer = colors.Normalize(scalar_min, scalar_max)
    colormap = cm.get_cmap("terrain")
    extent_km = float(np.max(np.abs(displayed_vertices_m)) / 1000.0)

    try:
        figure = plt.figure(figsize=(10.0, 9.5))
    except Exception as exception:
        # GUI backend failures (for example an inaccessible forwarded DISPLAY)
        # use backend-specific exception classes. Normalize them into one
        # actionable CLI error without hiding failures after initialization.
        raise RuntimeError(
            "cannot open the interactive 3-D window; verify DISPLAY and X "
            "forwarding, or use static camera options") from exception
    try:
        figure.canvas.manager.set_window_title("srcMoon rotatable surface")
    except AttributeError:
        pass
    axis = figure.add_subplot(111, projection="3d")
    configure_surface_axis(
        axis, triangle_km, face_scalar, normalizer, colormap, extent_km,
        show_edges, initial_elevation_deg, initial_azimuth_deg,
        "Drag with the mouse to rotate the surface")
    scalar_map = cm.ScalarMappable(norm=normalizer, cmap=colormap)
    scalar_map.set_array(np.asarray([], dtype=np.float64))
    figure.colorbar(scalar_map, ax=axis, fraction=0.035, pad=0.08,
                    label=scalar_label)
    figure.suptitle(
        (title or f"AMPS lunar surface: {source.name}") + "\n" +
        f"{frame_label}\n{exaggeration_notice(vertical_exaggeration)}\n"
        f"{NOTICE}", fontsize=12)
    figure.subplots_adjust(left=0.02, right=0.91, bottom=0.14, top=0.87)
    status = figure.text(
        0.5, 0.075,
        "Rotate with the mouse. Press S or click the button to save; Q closes.",
        ha="center", va="center")
    button_axis = figure.add_axes((0.34, 0.015, 0.32, 0.045))
    save_button = Button(button_axis, "Save selected view")
    saved: list[dict[str, str]] = []

    def update_camera_text(event: object) -> None:
        """Expose the selected angles after each mouse rotation."""
        if getattr(event, "inaxes", None) is axis:
            status.set_text(
                f"Selected elevation={axis.elev:.2f} deg, "
                f"azimuth={axis.azim:.2f} deg. Press S to save; Q closes.")
            figure.canvas.draw_idle()

    def save_camera(event: object) -> None:
        """Write clean figures without toolbar/button decorations."""
        nonlocal saved
        try:
            saved = save_selected_view(
                mesh, source, selected_outputs, color_by, reference_radius_m,
                frame_label, title, float(axis.elev), float(axis.azim),
                max_render_faces, show_edges, dpi, vertical_exaggeration)
            status.set_text(
                f"Saved elevation={axis.elev:.2f} deg, "
                f"azimuth={axis.azim:.2f} deg to " +
                ", ".join(str(path) for path in selected_outputs))
        except (OSError, ValueError) as exception:
            status.set_text(f"SAVE ERROR: {exception}")
        figure.canvas.draw_idle()

    def key_press(event: object) -> None:
        """Provide keyboard controls in addition to the visible button."""
        key = str(getattr(event, "key", "")).lower()
        if key == "s":
            save_camera(event)
        elif key == "q":
            plt.close(figure)

    # Keep the callback/button objects alive until pyplot's blocking show()
    # returns; otherwise some backends may garbage-collect the controls.
    callbacks = [
        figure.canvas.mpl_connect("button_release_event", update_camera_text),
        figure.canvas.mpl_connect("key_press_event", key_press),
    ]
    save_button.on_clicked(save_camera)
    print(
        "Interactive 3-D surface opened: drag to rotate; press S or click "
        "'Save selected view'; press Q to close.", file=sys.stderr)
    try:
        plt.show()
    except Exception as exception:
        raise RuntimeError(
            "interactive 3-D window failed; verify DISPLAY and X forwarding") \
            from exception
    del callbacks, save_button
    interaction = {
        "opened": True,
        "controls": "mouse rotate; S/button save selected view; Q close",
        "selected_view_outputs_requested": [str(path) for path in selected_outputs],
        "selected_view_outputs_saved": saved,
        "final_elevation_deg": float(axis.elev),
        "final_azimuth_deg": float(axis.azim),
        "vertical_exaggeration": vertical_exaggeration,
    }
    plt.close(figure)
    return interaction


def build_argument_parser() -> argparse.ArgumentParser:
    """Construct a self-contained CLI guide as well as the option parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Visualize a production srcMoon CEA or Tecplot FETRIANGLE "
            "surface. By default the script writes a static four-panel PNG "
            "and does not open a window."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
OUTPUT BEHAVIOR
  With no --output option, the static figure is written beside the input as
  <input-stem>_visualization.png. The JSON receipt printed to stdout gives the
  absolute output path and SHA-256. Existing figures are protected; use
  --force only when replacement is intended.

STATIC AND ROTATABLE WINDOWS
  --show opens the normal four-panel figure. Its two 3-D panels can be rotated,
  but toolbar saving writes the complete four-panel layout.

  --interactive-3d additionally opens a dedicated rotatable globe after the
  normal static figure is written. In that window:
    * drag with the mouse to rotate the lunar surface;
    * press S or click "Save selected view" to write a clean single-panel plot;
    * press Q to close the interactive window.

  The selected-view filename defaults to <input-stem>_selected_view.png.
  Repeat --selected-view-output to save the selected camera in several formats,
  for example PNG and EPS. The JSON receipt is printed after the window closes
  and records final camera angles and files saved with S/the button.

DISPLAY REQUIREMENTS
  Interactive mode needs a GUI Matplotlib backend and a working X display.
  For a remote session, connect with X forwarding. If MPLBACKEND=Agg is set,
  unset it before using --interactive-3d. Without an X display, create a static
  camera with --view-latitude-deg and --view-longitude-deg instead.

COORDINATES AND COLOR
  Mesh files contain centre-relative Cartesian metres but do not encode a
  coordinate-frame name or reference datum. Radius coloring is therefore the
  safe default. Elevation coloring requires both an authoritative radius and
  an explicit frame label from the source metadata; for LDEM_4 these are
  1737400 m and MEAN EARTH/POLAR AXIS OF DE421.

VERTICAL EXAGGERATION
  --vertical-exaggeration F multiplies only radial relief relative to the
  explicitly supplied --reference-radius-m:

    displayed_radius = reference_radius + F * (true_radius - reference_radius)

  F=1 is the physically true geometry and is the default. Values other than 1
  require --reference-radius-m. Colors, input coordinates, geometry checks,
  maps, histograms, and receipts retain true unexaggerated values. Every
  affected figure is prominently labeled VERTICAL EXAGGERATION and the factor
  is recorded in JSON. Exaggeration is presentation only, never model or V&V
  evidence.

LARGE MESHES
  Every face is rendered by default. --max-render-faces N requests a
  deterministic display subset for memory-limited plotting. Geometry checks
  and the radial histogram still use the complete mesh, and the subset is
  disclosed in the figure and JSON receipt.

EXAMPLES
  Static CEA figure in the current directory:
    visualize_surface.py lola_surface.cea

  Static elevation-colored PNG and EPS:
    visualize_surface.py lola_surface.cea \\
      --color-by elevation-m --reference-radius-m 1737400 \\
      --frame-label "MEAN EARTH/POLAR AXIS OF DE421" \\
      --output lola_surface.png --output lola_surface.eps

  Rotatable globe and saved selected camera:
    visualize_surface.py lola_surface.cea --interactive-3d --force \\
      --color-by elevation-m --reference-radius-m 1737400 \\
      --vertical-exaggeration 20 \\
      --frame-label "MEAN EARTH/POLAR AXIS OF DE421" \\
      --selected-view-output selected_view.png \\
      --selected-view-output selected_view.eps

All products are labeled VISUALIZATION ONLY - NOT V&V EVIDENCE.
""")

    input_group = parser.add_argument_group("input surface")
    input_group.add_argument(
        "surface", type=Path,
        help="production .cea or point-packed Tecplot FETRIANGLE file")
    input_group.add_argument(
        "--format", choices=("auto", "cea", "tecplot"), default="auto",
        help="input parser; auto uses the suffix/header (default: auto)")

    static_group = parser.add_argument_group("static figure")
    static_group.add_argument(
        "--output", type=Path, action="append", default=[],
        help="PNG/EPS/PDF/SVG output; repeat to save multiple formats")
    static_group.add_argument(
        "--view-latitude-deg", type=float, default=25.0,
        help="initial/static camera elevation in degrees (default: 25)")
    static_group.add_argument(
        "--view-longitude-deg", type=float, default=-60.0,
        help="initial/static camera azimuth in degrees (default: -60)")
    static_group.add_argument(
        "--show", action="store_true",
        help="open the four-panel figure; mutually exclusive with --interactive-3d")

    appearance_group = parser.add_argument_group("appearance and metadata")
    appearance_group.add_argument(
        "--color-by", choices=("radius-km", "elevation-m"),
        default="radius-km",
        help="surface color quantity (default: radius-km)")
    appearance_group.add_argument(
        "--reference-radius-m", type=float,
        help="authoritative datum required by --color-by elevation-m")
    appearance_group.add_argument(
        "--vertical-exaggeration", type=float, default=1.0,
        help="display-only radial relief factor; !=1 requires reference radius")
    appearance_group.add_argument(
        "--frame-label",
        default="input Cartesian frame (not encoded in mesh file)",
        help="coordinate-frame text printed on figures; never inferred")
    appearance_group.add_argument("--title", help="custom figure title")
    appearance_group.add_argument(
        "--max-render-faces", type=int, default=0,
        help="0 renders all; N>0 selects a disclosed deterministic subset")
    appearance_group.add_argument(
        "--show-edges", action="store_true",
        help="draw triangle edges (slow/dense for fine meshes)")
    appearance_group.add_argument(
        "--dpi", type=int, default=180,
        help="raster resolution in [50,1200] (default: 180)")
    appearance_group.add_argument(
        "--force", action="store_true",
        help="replace existing requested static/selected-view figures")

    interactive_group = parser.add_argument_group("rotatable 3-D globe")
    interactive_group.add_argument(
        "--interactive-3d", action="store_true",
        help="open mouse-rotatable globe with S/button save and Q close")
    interactive_group.add_argument(
        "--selected-view-output", type=Path, action="append", default=[],
        help="clean selected-camera PNG/EPS/PDF/SVG; repeat for multiple formats")
    return parser


def main() -> int:
    parser = build_argument_parser()
    arguments = parser.parse_args()
    try:
        source = arguments.surface.resolve()
        if not source.is_file():
            raise ValueError(f"surface file does not exist: {source}")
        if not math.isfinite(arguments.view_latitude_deg) or not (-90.0 <= arguments.view_latitude_deg <= 90.0):
            raise ValueError("--view-latitude-deg must be within [-90,90]")
        if not math.isfinite(arguments.view_longitude_deg):
            raise ValueError("--view-longitude-deg must be finite")
        if arguments.dpi < 50 or arguments.dpi > 1200:
            raise ValueError("--dpi must be within [50,1200]")
        if arguments.show and arguments.interactive_3d:
            raise ValueError("use either --show or --interactive-3d, not both")
        if arguments.selected_view_output and not arguments.interactive_3d:
            raise ValueError(
                "--selected-view-output requires --interactive-3d")
        if (arguments.interactive_3d and sys.platform.startswith("linux") and
                not os.environ.get("DISPLAY")):
            raise ValueError(
                "--interactive-3d requires an X display; use SSH X forwarding "
                "or create a static view with --view-latitude-deg and "
                "--view-longitude-deg")
        outputs = arguments.output or [
            source.with_name(source.stem + "_visualization.png")]
        resolved_outputs = [output.resolve() for output in outputs]
        selected_outputs = ([output.resolve()
                             for output in arguments.selected_view_output] or
                            [source.with_name(
                                source.stem + "_selected_view.png").resolve()])
        all_outputs = (resolved_outputs + selected_outputs
                       if arguments.interactive_3d else resolved_outputs)
        if len(set(all_outputs)) != len(all_outputs):
            raise ValueError("each figure output path must be unique")
        for output in all_outputs:
            if output == source:
                raise ValueError("a figure output cannot replace the input surface")
            if output.suffix.lower() not in {".png", ".eps", ".pdf", ".svg"}:
                raise ValueError(f"unsupported figure suffix for {output}")
            if output.exists() and not arguments.force:
                raise ValueError(
                    f"refusing to overwrite existing figure without --force: {output}")
        mesh = read_surface(source, arguments.format)
        geometry = validate_surface(mesh)
        rendering = render_surface(
            mesh, source, resolved_outputs, arguments.color_by,
            arguments.reference_radius_m, arguments.frame_label,
            arguments.title, arguments.view_latitude_deg,
            arguments.view_longitude_deg, arguments.max_render_faces,
            arguments.show_edges, arguments.dpi, arguments.show,
            arguments.interactive_3d, arguments.vertical_exaggeration)
        if arguments.interactive_3d:
            import matplotlib
            if matplotlib.get_backend().lower() == "agg":
                raise RuntimeError(
                    "--interactive-3d selected a non-GUI Matplotlib backend; "
                    "unset MPLBACKEND=Agg and verify DISPLAY/X forwarding")
            rendering["interactive_3d"] = run_interactive_3d(
                mesh, source, selected_outputs, arguments.color_by,
                arguments.reference_radius_m, arguments.frame_label,
                arguments.title, arguments.view_latitude_deg,
                arguments.view_longitude_deg, arguments.max_render_faces,
                arguments.show_edges, arguments.dpi,
                arguments.vertical_exaggeration)
        receipt = {
            "schema": "srcMoon-surface-visualization-v1",
            "notice": NOTICE,
            "source": {"path": str(source), "sha256": sha256(source),
                       "format": mesh.source_format},
            "geometry": geometry,
            "rendering": rendering,
        }
        print(json.dumps(receipt, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, ImportError, RuntimeError) as exception:
        print(f"ERROR: {exception}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

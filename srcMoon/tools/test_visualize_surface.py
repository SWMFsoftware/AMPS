#!/usr/bin/env python3
"""Infrastructure tests for the surface visualizer; no science claim is made."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import matplotlib
import numpy as np

# Infrastructure tests never open GUI windows. Select the headless backend
# before importing the tool so unittest order cannot accidentally select Tk.
matplotlib.use("Agg")


SCRIPT = Path(__file__).with_name("visualize_surface.py")
SPEC = importlib.util.spec_from_file_location("visualize_surface", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
visualizer = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = visualizer
SPEC.loader.exec_module(visualizer)


def tetrahedron() -> tuple[np.ndarray, np.ndarray]:
    """Return a small closed fixture with algorithmically outward winding."""
    vertices = np.asarray(((1.0, 1.0, 1.0), (-1.0, -1.0, 1.0),
                           (-1.0, 1.0, -1.0), (1.0, -1.0, -1.0)))
    faces = []
    for face in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)):
        indices = list(face)
        triangle = vertices[indices]
        normal = np.cross(triangle[1] - triangle[0],
                          triangle[2] - triangle[0])
        if float(np.dot(normal, triangle.sum(axis=0))) <= 0.0:
            indices[1], indices[2] = indices[2], indices[1]
        faces.append(indices)
    return vertices, np.asarray(faces, dtype=np.int64)


def write_cea(path: Path, vertices: np.ndarray, faces: np.ndarray) -> None:
    """Write the documented production CEA layout into a temporary fixture."""
    lines = [f"{len(vertices)} {len(faces)}"]
    lines.extend(" ".join(str(value) for value in vertex)
                 for vertex in vertices)
    lines.extend(" ".join(str(int(value) + 1) for value in face)
                 for face in faces)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_tecplot(path: Path, vertices: np.ndarray, faces: np.ndarray,
                  radius_offset: float = 0.0) -> None:
    """Write the exact point-packed FETRIANGLE subset consumed by the tool."""
    lines = [
        'VARIABLES="X [m]","Y [m]","Z [m]","Radius [m]"',
        f"ZONE N={len(vertices)}, E={len(faces)}, "
        "DATAPACKING=POINT, ZONETYPE=FETRIANGLE",
    ]
    for vertex in vertices:
        radius = math.sqrt(float(np.dot(vertex, vertex))) + radius_offset
        lines.append(" ".join(str(value) for value in (*vertex, radius)))
    lines.extend(" ".join(str(int(value) + 1) for value in face)
                 for face in faces)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class SurfaceVisualizationInfrastructure(unittest.TestCase):
    def test_help_documents_interactive_workflow(self):
        """The printed help is sufficient to discover and operate 3-D mode."""
        help_text = visualizer.build_argument_parser().format_help()
        for required in (
                "drag with the mouse", 'press S or click "Save selected view"',
                "press Q to close", "--selected-view-output",
                "working X display", "MPLBACKEND=Agg", "--force",
                "--reference-radius-m 1737400", "PNG and EPS",
                "--vertical-exaggeration", "displayed_radius =",
                "VISUALIZATION ONLY - NOT V&V EVIDENCE"):
            self.assertIn(required, help_text)

    def test_vertical_exaggeration_is_radial_and_non_mutating(self):
        """Display scaling changes radius only and preserves source vertices."""
        vertices, faces = tetrahedron()
        mesh = visualizer.SurfaceMesh(vertices.copy(), faces, "CEA")
        original = mesh.vertices_m.copy()
        shown = visualizer.display_vertices(mesh, 1.0, 3.0)
        true_radius = np.linalg.norm(original, axis=1)
        shown_radius = np.linalg.norm(shown, axis=1)
        self.assertTrue(np.allclose(shown_radius,
                                    1.0 + 3.0 * (true_radius - 1.0)))
        self.assertTrue(np.allclose(shown / shown_radius[:, np.newaxis],
                                    original / true_radius[:, np.newaxis]))
        self.assertTrue(np.array_equal(mesh.vertices_m, original))
        with self.assertRaisesRegex(ValueError, "reference-radius-m"):
            visualizer.display_vertices(mesh, None, 3.0)

    def test_production_formats_have_identical_validated_geometry(self):
        """CEA and Tecplot parsers retain all fixture vertices and faces."""
        vertices, faces = tetrahedron()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cea = root / "surface.cea"
            tecplot = root / "surface.dat"
            write_cea(cea, vertices, faces)
            write_tecplot(tecplot, vertices, faces)
            parsed_cea = visualizer.read_surface(cea, "auto")
            parsed_tecplot = visualizer.read_surface(tecplot, "auto")
            self.assertTrue(np.array_equal(parsed_cea.vertices_m,
                                           parsed_tecplot.vertices_m))
            self.assertTrue(np.array_equal(parsed_cea.faces,
                                           parsed_tecplot.faces))
            self.assertEqual(visualizer.validate_surface(parsed_cea)["face_count"], 4)
            self.assertEqual(visualizer.validate_surface(parsed_tecplot)["vertex_count"], 4)

    def test_inconsistent_stored_radius_is_rejected(self):
        """A misleading Tecplot scalar must fail before rendering."""
        vertices, faces = tetrahedron()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.dat"
            write_tecplot(path, vertices, faces, radius_offset=1.0)
            with self.assertRaisesRegex(ValueError, "Radius column"):
                visualizer.validate_surface(
                    visualizer.read_surface(path, "tecplot"))

    def test_static_render_and_explicit_elevation_reference(self):
        """Rendering writes a figure and discloses any face subset."""
        vertices, faces = tetrahedron()
        mesh = visualizer.SurfaceMesh(vertices, faces, "CEA")
        visualizer.validate_surface(mesh)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "view.png"
            receipt = visualizer.render_surface(
                mesh, root / "fixture.cea", [output], "elevation-m",
                math.sqrt(3.0), "fixture Cartesian frame", "fixture", 25.0,
                -60.0, 2, True, 80, False)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)
            self.assertEqual(receipt["rendered_face_count"], 2)
            self.assertTrue(receipt["deterministic_subset"])

    def test_selected_camera_writer_uses_requested_angles(self):
        """The interactive save path creates a clean standalone 3-D image."""
        vertices, faces = tetrahedron()
        mesh = visualizer.SurfaceMesh(vertices, faces, "CEA")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "selected.png"
            artifacts = visualizer.save_selected_view(
                mesh, root / "fixture.cea", [output], "radius-km", None,
                "fixture Cartesian frame", "selected fixture", 37.5, 122.0,
                0, False, 80)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)
            self.assertEqual(artifacts[0]["path"], str(output.resolve()))
            self.assertEqual(artifacts[0]["sha256"],
                             visualizer.sha256(output))

    def test_interactive_controller_is_headless_testable(self):
        """The controller wires controls without saving before user action."""
        vertices, faces = tetrahedron()
        mesh = visualizer.SurfaceMesh(vertices, faces, "CEA")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            selected = root / "not_saved_without_user_action.png"
            with mock.patch("matplotlib.pyplot.show"):
                receipt = visualizer.run_interactive_3d(
                    mesh, root / "fixture.cea", [selected], "radius-km", None,
                    "fixture Cartesian frame", "interactive fixture", 25.0,
                    -60.0, 0, False, 80)
            self.assertTrue(receipt["opened"])
            self.assertEqual(receipt["selected_view_outputs_saved"], [])
            self.assertFalse(selected.exists())


if __name__ == "__main__":
    unittest.main()

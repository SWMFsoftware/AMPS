#!/usr/bin/env python3
"""Independent parser, classification, and headless-render tests.

``--help`` intentionally prints the production viewer's complete usage guide
after the test-specific commands. Sourcing that text from the viewer parser
keeps variable, camera, color-range, and movie examples synchronized with the
CLI that users actually run instead of maintaining a second stale help page.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
VIEWER = ROOT / "srcSEP3D" / "examples" / "shock-front" / "view_front.py"
SPEC = importlib.util.spec_from_file_location("shock_front_viewer", VIEWER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def print_usage_guide() -> None:
    """Print test invocation examples followed by authoritative viewer help."""

    print(
        "Shock-front viewer test and usage guide\n"
        "========================================\n\n"
        "View a front file through this convenience entry point:\n"
        "  python3 srcSEP3D/test/test_shock_front_viewer.py FRONT.dat "
        "--variable fast_mach\n\n"
        "List the variables in a front file:\n"
        "  python3 srcSEP3D/test/test_shock_front_viewer.py FRONT.dat "
        "--list-variables\n\n"
        "Run the regression suite:\n"
        "  python3 srcSEP3D/test/test_shock_front_viewer.py\n\n"
        "Run one named regression:\n"
        "  python3 srcSEP3D/test/test_shock_front_viewer.py "
        "ShockFrontViewerTests.test_help_is_a_self_contained_usage_guide\n\n"
        "The production-viewer guide follows. Its examples show how to find "
        "snapshots and variables, plot a field, select a colormap and physical "
        "range, fix the viewing direction, and create GIF or MP4 movies.\n"
    )
    # Format the real parser rather than copying its options here. Besides
    # preventing documentation drift, this makes the test runner's help an
    # executable check of the exact CLI installed beside the example.
    viewer_parser = MODULE._parser()
    viewer_parser.prog = "srcSEP3D/examples/shock-front/view_front.py"
    print(viewer_parser.format_help())


def is_viewer_invocation(arguments: list[str]) -> bool:
    """Distinguish front inputs from optional ``unittest`` selectors.

    ``unittest.main`` treats every positional argument as a Python test name,
    which produces a misleading ``_FailedTest`` error when a user naturally
    passes a Tecplot ``*-front.dat`` file. Front products have a stable ``.dat``
    suffix; directory and quoted-glob inputs are also part of the production
    viewer contract. Those inputs are therefore forwarded to ``MODULE.main``.
    Named tests such as ``ShockFrontViewerTests.test_parser...`` do not match
    any of these forms and retain the standard unittest behavior.
    """

    for argument in arguments:
        if argument.startswith("-"):
            continue
        if argument.lower().endswith(".dat"):
            return True
        if "*" in argument or "?" in argument or "[" in argument:
            return True
        if Path(argument).is_dir():
            return True
    return False


def fixture_text(time_s: float = 10.0, generation: int = 2) -> str:
    variables = (
        "x_m", "y_m", "z_m", "time_s", "generation",
        "shock_accepted", "fast_mach", "density_compression",
    )
    # Quads 1 and 2 are homogeneous accepted/non-shock controls.  Quad 3 is
    # deliberately mixed and proves the viewer does not assign it an averaged
    # physical classification.
    rows = (
        (0, 0, 0, time_s, generation, 1, 2.0, 2.5),
        (1, 0, 0, time_s, generation, 1, 2.2, 2.6),
        (1, 1, 0, time_s, generation, 1, 2.4, 2.7),
        (0, 1, 0, time_s, generation, 1, 2.6, 2.8),
        (0, 0, 1, time_s, generation, 0, 0.0, 1.0),
        (1, 0, 1, time_s, generation, 0, 0.0, 1.0),
        (1, 1, 1, time_s, generation, 0, 0.0, 1.0),
        (0, 1, 1, time_s, generation, 0, 0.0, 1.0),
    )
    text = [
        'TITLE="viewer fixture"',
        "VARIABLES=" + ",".join(f'\"{name}\"' for name in variables),
        'ZONE T="fixture", N=8, E=3, DATAPACKING=POINT, '
        'ZONETYPE=FEQUADRILATERAL',
        'AUXDATA volume_role="ambient-reference-only"',
    ]
    text.extend(" ".join(str(value) for value in row) for row in rows)
    text.extend(("1 2 3 4", "5 6 7 8", "3 4 8 7"))
    return "\n".join(text) + "\n"


def triangle_fixture_text(time_s: float = 10.0, generation: int = 2) -> str:
    """Small mixed-location BLOCK fixture matching the production contract."""

    variables = (
        "x_m", "y_m", "z_m", "time_s", "generation",
        "shock_accepted", "fast_mach", "density_compression",
    )
    nodal = (
        (0.0, 1.0, 1.0, 0.0, 0.5),
        (0.0, 0.0, 1.0, 1.0, 0.5),
        (0.0, 0.0, 0.0, 0.0, 1.0),
    )
    cell = (
        (time_s, time_s, time_s),
        (generation, generation, generation),
        (1.0, 0.0, 1.0),
        (2.0, 0.0, 2.4),
        (2.5, 1.0, 2.8),
    )
    text = [
        'TITLE="triangle viewer fixture"',
        "VARIABLES=" + ",".join(f'\"{name}\"' for name in variables),
        'ZONE T="fixture", N=5, E=3, DATAPACKING=BLOCK, '
        'ZONETYPE=FETRIANGLE, VARLOCATION=([4-8]=CELLCENTERED)',
        'AUXDATA surface_topology="triangular-sse-cap-v1"',
    ]
    text.extend(" ".join(str(value) for value in block) for block in nodal)
    text.extend(" ".join(str(value) for value in block) for block in cell)
    text.extend(("1 2 5", "2 3 5", "3 4 5"))
    return "\n".join(text) + "\n"


class ShockFrontViewerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.front = self.root / "front.dat"
        self.front.write_text(fixture_text(), encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_help_is_a_self_contained_usage_guide(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(VIEWER), "--help"],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        for required in (
                "Examples", "Find the available front snapshots",
                "--list-variables", "--variable fast_mach",
                "--cmap plasma", "--accepted-only", "--show-sun",
                "print(*p.colormaps(), sep=chr(10))",
                "--show-distance-axis", "Selecting a movie viewpoint",
                "press the `v` key", "--elev 24 --azim -58",
                "--elev 20 --azim 135 --length-unit au",
                "--movie front.gif", "--movie front.mp4",
                "A .mp4 destination requires ffmpeg",
                "does not create sheath/ejecta volume data"):
            self.assertIn(required, completed.stdout)

    def test_test_runner_help_exposes_viewer_workflows(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--help"],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        for required in (
                "Shock-front viewer test and usage guide",
                "Run the regression suite", "--list-variables",
                "--variable fast_mach", "--cmap plasma --vmin 0 --vmax",
                "--elev 20 --azim 135", "--movie front.gif",
                "--movie front.mp4", "Selecting a movie viewpoint"):
            self.assertIn(required, completed.stdout)

    def test_test_runner_forwards_front_file_to_production_viewer(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), str(self.front),
             "--list-variables"],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stdout.splitlines(), [
            "x_m", "y_m", "z_m", "time_s", "generation",
            "shock_accepted", "fast_mach", "density_compression",
        ])

    def test_parser_and_unanimous_acceptance_masks(self) -> None:
        surface = MODULE.load_front(self.front)
        accepted, nonshock, boundary = MODULE.classify_faces(surface)
        self.assertEqual(surface.data.shape, (8, 8))
        self.assertEqual(surface.quads.shape, (3, 4))
        self.assertEqual(accepted.tolist(), [True, False, False])
        self.assertEqual(nonshock.tolist(), [False, True, False])
        self.assertEqual(boundary.tolist(), [False, False, True])
        self.assertEqual(MODULE.face_values(surface, "fast_mach").tolist(),
                         [2.3, 0.0, 1.25])

    def test_triangle_block_parser_uses_face_centered_classification(self) -> None:
        triangle = self.root / "triangle-front.dat"
        triangle.write_text(triangle_fixture_text(), encoding="utf-8")
        surface = MODULE.load_front(triangle)
        accepted, nonshock, boundary = MODULE.classify_faces(surface)
        self.assertEqual(surface.zone_type, "FETRIANGLE")
        self.assertEqual(surface.faces.shape, (3, 3))
        self.assertEqual(accepted.tolist(), [True, False, True])
        self.assertEqual(nonshock.tolist(), [False, True, False])
        self.assertEqual(boundary.tolist(), [False, False, False])
        self.assertEqual(MODULE.face_values(surface, "fast_mach").tolist(),
                         [2.0, 0.0, 2.4])

    def test_connectivity_outside_node_range_is_rejected(self) -> None:
        self.front.write_text(fixture_text().replace("1 2 3 4\n", "1 2 3 9\n"),
                              encoding="utf-8")
        with self.assertRaisesRegex(MODULE.FrontFormatError,
                                    "outside 1..8"):
            MODULE.load_front(self.front)

    def test_headless_accepted_only_render_and_variable_listing(self) -> None:
        image = self.root / "accepted.png"
        completed = subprocess.run(
            [sys.executable, str(VIEWER), str(self.front),
             "--variable", "fast_mach", "--accepted-only",
             "--show-sun", "--show-distance-axis",
             "--save", str(image), "--no-show"],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("accepted_faces=1", completed.stdout)
        self.assertIn("nonshock_faces=1", completed.stdout)
        self.assertIn("mixed_boundary_faces=1", completed.stdout)
        self.assertTrue(image.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

        listed = subprocess.run(
            [sys.executable, str(VIEWER), str(self.front), "--list-variables"],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False,
        )
        self.assertEqual(listed.returncode, 0, listed.stderr)
        self.assertEqual(listed.stdout.splitlines(), [
            "x_m", "y_m", "z_m", "time_s", "generation",
            "shock_accepted", "fast_mach", "density_compression",
        ])

    def test_movie_sorts_epochs_and_freezes_requested_view(self) -> None:
        later = self.root / "front-later.dat"
        later.write_text(fixture_text(time_s=20.0, generation=3), encoding="utf-8")
        movie = self.root / "front.gif"
        # Pass the later file first.  The movie summary must still report the
        # embedded physical time range in ascending order; filenames do not
        # define temporal order.
        completed = subprocess.run(
            [sys.executable, str(VIEWER), str(later), str(self.front),
             "--movie", str(movie), "--fps", "2",
             "--accepted-only", "--show-sun", "--show-distance-axis",
             "--elev", "31", "--azim", "-42", "--no-show"],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("movie_frames=2", completed.stdout)
        self.assertIn("time_range_s=[10,20]", completed.stdout)
        self.assertIn("fixed_view_elev=31", completed.stdout)
        self.assertIn("fixed_view_azim=-42", completed.stdout)
        self.assertIn(movie.read_bytes()[:6], (b"GIF87a", b"GIF89a"))


if __name__ == "__main__":
    command_arguments = sys.argv[1:]
    if any(argument in ("-h", "--help") for argument in command_arguments):
        print_usage_guide()
        raise SystemExit(0)
    if is_viewer_invocation(command_arguments):
        raise SystemExit(MODULE.main(command_arguments))
    unittest.main(verbosity=2)

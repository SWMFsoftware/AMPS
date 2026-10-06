#!/usr/bin/env python3
"""Independent parser, classification, and headless-render tests."""

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
                "Examples", "--list-variables", "--variable fast_mach",
                "--cmap plasma", "--accepted-only", "--show-sun",
                "--show-distance-axis", "Selecting a movie viewpoint",
                "press the `v` key", "--elev 24 --azim -58",
                "--movie front.gif", "A .mp4 destination requires ffmpeg",
                "does not create sheath/ejecta volume data"):
            self.assertIn(required, completed.stdout)

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
        self.assertIn("accepted_quads=1", completed.stdout)
        self.assertIn("nonshock_quads=1", completed.stdout)
        self.assertIn("mixed_boundary_quads=1", completed.stdout)
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
    unittest.main(verbosity=2)

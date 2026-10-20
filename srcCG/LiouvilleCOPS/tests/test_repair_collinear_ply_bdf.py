"""Exact edge-flip invariants; near-degenerate geometry must not be moved."""
import sys
from pathlib import Path
import unittest
import tempfile
import struct
import json
import os
import subprocess

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from repair_collinear_ply_bdf import plan_repairs, repair
from convert_ply_to_bdf import convert
from gravity_topology import require_outward_closed_body


class RepairTests(unittest.TestCase):
    def setUp(self):
        self.v = np.array([[0., 0, 0], [2., 0, 0], [1., 0, 0], [0., 1, 0]])
        self.f = np.array([[0, 1, 2], [1, 0, 3]])

    def test_exact_partition(self):
        changes, records = plan_repairs(self.v, self.f)
        self.assertEqual(changes, {0: (2, 0, 3), 1: (2, 3, 1)})
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["interior_fraction"], .5)

    def test_existing_diagonal_rejected(self):
        with self.assertRaisesRegex(ValueError, "diagonal"):
            plan_repairs(self.v, np.vstack([self.f, [2, 3, 0]]))

    def test_open_edge_rejected(self):
        with self.assertRaisesRegex(ValueError, "two-sided"):
            plan_repairs(self.v, self.f[:1])

    def test_near_collinear_is_not_repaired(self):
        self.v[2, 2] = 1e-15
        self.assertEqual(plan_repairs(self.v, self.f), ({}, []))

    def test_wrong_winding_rejected(self):
        self.f[1] = [0, 1, 3]
        with self.assertRaisesRegex(ValueError, "orientation"):
            plan_repairs(self.v, self.f)

    def test_closed_surface_round_trip(self):
        # Insert a midpoint into one tetrahedron edge, with an extra collinear
        # face joining the unsplit opposite face. The flip must restore a
        # closed manifold without changing coordinates, area, or volume.
        vertices = np.array([[0., 0, 0], [2., 0, 0], [0., 2, 0],
                             [0., 0, 2], [1., 0, 0]])
        faces = np.array([[0, 2, 1], [0, 1, 4], [0, 4, 3],
                          [4, 1, 3], [1, 2, 3], [2, 0, 3]])
        changes, _ = plan_repairs(vertices, faces)
        fixed = faces.copy()
        for i, face in changes.items():
            fixed[i] = face
        audit = require_outward_closed_body(vertices, fixed)
        self.assertAlmostEqual(audit["signed_volume_m3"], 4./3.)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, original, output = (root/"a.ply", root/"a.bdf", root/"b.bdf")
            header = ("ply\nformat binary_little_endian 1.0\n"
                      "element vertex 5\nproperty float x\nproperty float y\n"
                      "property float z\nelement face 6\n"
                      "property list int int vertex_index\nend_header\n")
            with source.open("wb") as stream:
                stream.write(header.encode())
                stream.write(vertices.astype("<f4").tobytes())
                for face in faces:
                    stream.write(struct.pack("<4i", 3, *face))
            convert(source, original, "m", "test", "synthetic fixture", True)
            repair(source, original, output)
            old = original.read_text().splitlines()
            new = output.read_text().splitlines()
            self.assertEqual(sum(a != b for a, b in zip(old, new)), 2)
            self.assertEqual(len(old), len(new))
            report = json.loads(output.with_suffix(".bdf.repair.json").read_text())
            self.assertEqual(report["vertices_moved"], 0)
            if os.environ.get("COPS_MESH_PROBE"):
                subprocess.run([os.environ["COPS_MESH_PROBE"], str(output)],
                               check=True, capture_output=True, text=True)
            with self.assertRaisesRegex(ValueError, "already exists"):
                repair(source, original, output)


if __name__ == "__main__":
    unittest.main()

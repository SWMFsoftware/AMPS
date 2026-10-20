"""Synthetic endian/unit/order tests independent of the downloaded shape."""
from pathlib import Path
import os
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from convert_ply_to_bdf import convert, read_ply
from build_physical_activity_templates import read_nastran_surface
from convert_dsk_to_bdf import write_grid, write_triangle


class PlyConversionTests(unittest.TestCase):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])

    def test_eight_digit_adjacent_node_ids(self):
        # Four sparse IDs exercise the same boundary as the 44M mesh without
        # allocating tens of millions of synthetic vertices.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "large_ids.bdf"
            with path.open("w") as stream:
                for index, vertex in enumerate(self.vertices, 10000000):
                    write_grid(stream, index, vertex / 1000.)
                for index, face in enumerate(self.faces, 10000000):
                    write_triangle(stream, index, tuple(int(i) + 10000000 for i in face))
            vertices, faces = read_nastran_surface(path)
            np.testing.assert_array_equal(vertices, self.vertices)
            np.testing.assert_array_equal(faces, self.faces)
            # Optional production C++ reader check. The executable must call
            # ReadNastran and FinalizeGeometry(requireClosed=true); keeping it
            # opt-in avoids requiring CSPICE in lightweight Python test runs.
            if os.environ.get("COPS_MESH_PROBE"):
                result = subprocess.run([os.environ["COPS_MESH_PROBE"], str(path)],
                                        text=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
                self.assertEqual(result.returncode, 0, result.stdout)
                self.assertIn("4 CTRIA3", result.stdout)

    def fixture(self, path, endian=">", bad_count=False):
        order = "big" if endian == ">" else "little"
        header = ("ply\nformat binary_%s_endian 1.0\n"
                  "element vertex 4\nproperty float x\nproperty float y\nproperty float z\n"
                  "element face 4\nproperty list int int vertex_index\nend_header\n") % order
        with path.open("wb") as stream:
            stream.write(header.encode("ascii"))
            for v in self.vertices:
                stream.write(struct.pack(endian + "3f", *v))
            for f in self.faces:
                stream.write(struct.pack(endian + "4i", 4 if bad_count else 3, *f))

    def test_both_endians_and_metre_round_trip(self):
        for endian in (">", "<"):
            with tempfile.TemporaryDirectory() as directory:
                source, output = Path(directory) / "a.ply", Path(directory) / "a.bdf"
                self.fixture(source, endian)
                result = convert(source, output, "km", "synthetic", "test fixture")
                vertices, faces = read_nastran_surface(output)
                np.testing.assert_array_equal(vertices, self.vertices * 1000)
                np.testing.assert_array_equal(faces, self.faces)
                self.assertTrue(result["gravity_topology_audit"]["passed"])
                with self.assertRaisesRegex(ValueError, "existing"):
                    convert(source, output, "km", "synthetic", "test fixture")

    def test_refuse_nontriangular_and_truncated_payload(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "a.ply"
            self.fixture(source, bad_count=True)
            with self.assertRaisesRegex(ValueError, "non-triangular"):
                read_ply(source)
            self.fixture(source)
            with source.open("ab") as stream:
                stream.write(b"x")
            with self.assertRaisesRegex(ValueError, "size"):
                read_ply(source)

    def test_inward_surface_not_silently_repaired(self):
        original = self.faces
        try:
            self.faces = original[:, ::-1]
            with tempfile.TemporaryDirectory() as directory:
                source, output = Path(directory) / "a.ply", Path(directory) / "a.bdf"
                self.fixture(source)
                with self.assertRaisesRegex(ValueError, "inward"):
                    convert(source, output, "m", "synthetic", "test fixture")
                self.assertFalse(output.exists())
                report = convert(source, output, "m", "synthetic", "test fixture", True)
                self.assertFalse(report["gravity_topology_audit"]["passed"])
                _, faces = read_nastran_surface(output)
                np.testing.assert_array_equal(faces, self.faces)
        finally:
            self.faces = original


if __name__ == "__main__":
    unittest.main()

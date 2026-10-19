"""Deterministic tests for the DSK-export to NASTRAN conversion layer."""

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

import convert_dsk_to_bdf as converter


class DskToBdfTests(unittest.TestCase):
    def test_vertex_facet_export_becomes_meter_bdf(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            export = root / "tetra.obj"
            export.write_text(
                "v 1.0 0.0 0.0\n"
                "v 0.0 1.0 0.0\n"
                "v 0.0 0.0 1.0\n"
                "v 0.0 0.0 0.0\n"
                "f 1 3 2\n"
                "f 1 2 4\n"
                "f 2 3 4\n"
                "f 3 1 4\n")
            output = root / "tetra.bdf"
            temporary, vertices, input_plates, output_plates, removed = \
                converter.convert_export(
                    export, output, 4, 4, "tetra.bds", "abc123")
            text = temporary.read_text()
            self.assertEqual((vertices, input_plates, output_plates), (4, 4, 4))
            self.assertEqual(removed, [])
            self.assertEqual(text.count("GRID*"), 4)
            self.assertEqual(text.count("CTRIA3"), 4)
            # Mesh.cpp reads x from columns 40:56 and z from continuation
            # columns 8:24. Verify both the exact field placement and km->m.
            lines = text.splitlines()
            first = next(i for i, line in enumerate(lines)
                         if line.startswith("GRID*"))
            self.assertAlmostEqual(float(lines[first][40:56]), 1000.0)
            self.assertAlmostEqual(float(lines[first + 1][8:24]), 0.0)
            triangle = next(line for line in lines if line.startswith("CTRIA3"))
            self.assertEqual([int(triangle[offset:offset + 8])
                              for offset in range(8, 48, 8)],
                             [1, 1, 1, 3, 2])
            self.assertTrue(text.endswith("ENDDATA\n"))
            temporary.unlink()

    def test_bad_counts_and_faces_are_rejected_without_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            export = root / "bad.obj"
            export.write_text("v 0 0 0\nv 1 0 0\nf 1 1 2\n")
            output = root / "bad.bdf"
            with self.assertRaisesRegex(ValueError, "degenerate"):
                converter.convert_export(export, output, 2, 1, "bad.bds", "x")
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob("bad.bdf.*.tmp")), [])

    def test_continuation_is_fixed_width_and_bounded(self):
        self.assertEqual(converter.continuation(1572866), "*0180002")
        with self.assertRaises(ValueError):
            converter.continuation(0)

    def test_reversed_duplicate_fin_is_removed(self):
        faces = np.array([
            [1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4],
            [1, 2, 5], [1, 5, 2]], dtype=np.int32)
        retained, pairs = converter.sanitize_faces(faces)
        self.assertEqual(pairs, [(5, 6)])
        self.assertEqual(retained.tolist(), faces[:4].tolist())


if __name__ == "__main__":
    unittest.main()

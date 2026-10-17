#!/usr/bin/env python3
"""Test the demo receipt and independent comparison helpers without D01 writes."""

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("lola_surface_demo", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
demo = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(demo)


class DemoInfrastructure(unittest.TestCase):
    def test_receipt_derivation_changes_only_declared_surface_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / "application.in"
            demo.prepare_input(
                demo.REPO / "srcMoon/examples/lola_surface.in", destination,
                125000.0, root / "surface.cea", root / "surface.dat")
            text = destination.read_text(encoding="utf-8")
            self.assertIn("surface_mesh_resolution_m = 125000", text)
            self.assertIn(f"surface_cea_file = {(root / 'surface.cea').resolve()}",
                          text)
            self.assertIn("surface_inventory_checkpoint_file = none", text)
            self.assertEqual(text.count("#begin application: moon"), 1)
            self.assertEqual(text.count("#end application"), 1)

    def test_adapter_summary_requires_complete_topology_receipt(self):
        complete = "\n".join((
            "surface_resolution_m=100000", "reference_radius_m=1737400",
            "subdivision_level=6", "vertex_count=40962", "face_count=81920",
            "maximum_edge_length_m=95000", "minimum_radius_m=1730000",
            "maximum_radius_m=1747000", "degenerate_face_count=0",
            "inward_face_count=0"))
        parsed = demo.parse_summary(complete)
        self.assertEqual(parsed["vertex_count"], 40962.0)
        with self.assertRaises(ValueError):
            demo.parse_summary(complete.replace("inward_face_count=0", ""))

    def test_status_vocabulary_is_not_promoted_by_demo_success(self):
        """Freeze the formal status fields documented by the demo runner."""
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn('"capability_label": "EXPERIMENTAL"', source)
        self.assertIn('"artifact_class": "MEETING_DEMONSTRATION"', source)
        self.assertIn('"execution": "NOT_RUN"', source)
        self.assertIn('"status": "SKIPPED"', source)
        self.assertIn('"validation_status": "NOT_VALIDATED"', source)


if __name__ == "__main__":
    unittest.main()

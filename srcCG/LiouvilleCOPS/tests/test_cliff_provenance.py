"""A relocated HEAL bundle must never regain eligibility for cliff fitting."""
import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from activity_templates import select_templates
from gravity_cache import ALGORITHM
from gravity_topology import ORIENTATION_POLICY


class CliffProvenanceTests(unittest.TestCase):
    def bundle(self):
        return {"path": "relocated.npz", "names": ["late_ng_not_rg", "cliff_ge_60deg"],
                "values": np.eye(2), "metadata": {}}

    def test_old_cliff_rejected_but_geometry_column_remains_usable(self):
        b = self.bundle()
        with self.assertRaisesRegex(ValueError, "HEAL-era"):
            select_templates(b, ["cliff_ge_60deg"])
        np.testing.assert_array_equal(select_templates(b, ["late_ng_not_rg"]), [[1.], [0.]])

    def test_repaired_policy_is_required_even_when_algorithm_matches(self):
        b = self.bundle()
        b["metadata"]["gravity"] = {
            "integrity_check": "DISABLE after topology audit",
            "cache": {"definition": {"algorithm": ALGORITHM,
                                      "orientation_policy": ORIENTATION_POLICY}}}
        np.testing.assert_array_equal(select_templates(b, ["cliff_ge_60deg"]), [[0.], [1.]])
        b["metadata"]["gravity"]["cache"]["definition"].pop("orientation_policy")
        with self.assertRaisesRegex(ValueError, "orientation provenance"):
            select_templates(b, ["cliff_ge_60deg"])


if __name__ == "__main__":
    unittest.main()

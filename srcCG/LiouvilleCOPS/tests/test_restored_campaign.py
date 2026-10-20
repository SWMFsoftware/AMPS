"""Preparation must preserve free baseline parameters and reject HEAL cliffs."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from activity_templates import write_bundle
from prepare_restored_campaign import compare_selections, fit_argv, settings_audit, template_precheck
from gravity_cache import ALGORITHM
from gravity_topology import ORIENTATION_POLICY


class RestoredCampaignTests(unittest.TestCase):
    def test_common_selection_does_not_use_obsolete_exclusion_flags(self):
        # Real profile expects 1150 records, with 216 NG gaps. The final four
        # points become outside on the new mesh. All old useNG flags below
        # are zero: reconstruction must use the shared all_available policy.
        with tempfile.TemporaryDirectory() as directory:
            old, new = (Path(directory) / name for name in ("old.txt", "new.txt"))
            old_lines = ["#iPoint UTC pNG_used_Pa pRG_used_Pa"]
            new_lines = ["#iPoint UTC inside refNG refRG useNG useRG"]
            for index in range(1150):
                ng = 1 if index >= 216 else -1
                old_lines.append("%d 2016-09-30T0:00:04 %d %d" %
                                 (index, ng if index < 1146 else -1, 2 if index < 1146 else -1))
                new_lines.append("%d 2016-09-30T00:00:04 0 %d 2 0 1" % (index, ng))
            old.write_text("\n".join(old_lines))
            new.write_text("\n".join(new_lines))
            report = compare_selections(old, new)
            self.assertEqual(report["NG"]["common_selected"], 930)
            self.assertEqual(report["RG"]["common_selected"], 1146)
            self.assertEqual(report["NG"]["current_only_ids"], [1146, 1147, 1148, 1149])
            self.assertEqual(report["NG"]["changed_common_pressure_ids"], [])

    def test_full_legacy_baseline_is_free(self):
        command = fit_argv("base", "out")
        self.assertEqual(command[command.index("--lambda") + 1], "0")
        self.assertEqual(command[command.index("--jobs") + 1], "1")
        self.assertEqual(command[command.index("--degree") + 1], "4")
        self.assertFalse(any(flag.startswith("--fix-") for flag in command))
        self.assertNotIn("--ordered-temperature", command)

    def test_path_move_does_not_hide_changed_physics(self):
        a = {"effective_settings": {"ShapeModel": ["old"], "CastShadows": ["on"]},
             "dependencies": [{"role": "observations", "file": {"sha256": "abc"}}]}
        b = {"effective_settings": {"ShapeModel": ["new"], "CastShadows": ["off"]},
             "dependencies": [{"role": "observations", "file": {"sha256": "def"}}]}
        report = settings_audit(a, b)
        self.assertEqual(report["unexpected_setting_changes"], ["CastShadows"])
        self.assertFalse(report["dependency_roles"]["observations"]["same_recorded_content"])

    def test_heal_cliff_refused_but_independent_ng_allowed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bundle.npz"
            write_bundle(path, ["late_ng_not_rg", "cliff_ge_60deg"],
                         np.ones((1, 2)), np.ones((1, 3)), np.ones(1), {})
            self.assertEqual(template_precheck(path, ["late_ng_not_rg"]), [])
            self.assertTrue(template_precheck(path, ["cliff_ge_60deg"]))
            self.assertTrue(template_precheck(path, ["missing"]))

    def test_corrected_policy_passes_metadata_precheck(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bundle.npz"
            metadata = {"gravity": {"integrity_check": "DISABLE after topology audit", "cache": {"definition": {
                "algorithm": ALGORITHM, "orientation_policy": ORIENTATION_POLICY}}}}
            write_bundle(path, ["cliff_ge_60deg"], np.ones((1, 1)),
                         np.ones((1, 3)), np.ones(1), metadata)
            self.assertEqual(template_precheck(path, ["cliff_ge_60deg"]), [])


if __name__ == "__main__":
    unittest.main()

"""Fast Stage-A1 selection-profile tests; no mission kernels are required."""

import contextlib
import io
from pathlib import Path
import sys
import unittest


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

from cops_selection import (PROFILE_SCHEMA, SelectionError,
                            apply_selection_profile, load_profile)
import run_fit_profile


class SelectionProfileTests(unittest.TestCase):
    columns = {name: index for index, name in enumerate(
        ("iPoint", "UTC", "inside", "refNG", "refRG", "useNG", "useRG"))}

    def rows(self):
        # Source use flags are deliberately zero: named profiles derive the
        # selection from pressure/geometry/time policy rather than stale flags.
        return [
            ["0", "2016-09-30T10:00:00", "0", "1", "2", "0", "0"],
            ["1", "2016-09-30T10:00:10", "0", "1", "2", "0", "0"],
            ["2", "2016-09-30T10:00:20", "0", "1", "2", "0", "0"],
            ["3", "2016-09-30T10:00:30", "1", "1", "2", "0", "0"],
            ["4", "2016-09-30T10:00:40", "0", "-1", "-1", "0", "0"],
        ]

    def profile(self, **changes):
        value = {
            "schema": PROFILE_SCHEMA,
            "schema_version": 1,
            "name": "fixture",
            "description": "test",
            "candidate_policy": "positive_outside_shape",
            "apply_fit_interval": True,
            "apply_run_exclusions": False,
            "exclude_intervals": {"NG": [], "RG": []},
            "saturated_ng": "include",
            "expected_counts": {"rows": 5, "used_ng": 3, "used_rg": 3},
        }
        value.update(changes)
        return value

    def test_explicit_interval_is_closed_and_reasons_are_audited(self):
        profile = self.profile(
            exclude_intervals={
                "NG": [["2016-09-30T10:00:10", "2016-09-30T10:00:20"]],
                "RG": [],
            },
            expected_counts={"rows": 5, "used_ng": 1, "used_rg": 3})
        selected, report = apply_selection_profile(
            self.rows(), self.columns, profile,
            {"FitInterval": [], "ExcludeNG": [], "ExcludeRG": []},
            [False, True, False, False, False])
        self.assertEqual([row[self.columns["useNG"]] for row in selected],
                         ["1", "0", "0", "0", "0"])
        self.assertEqual(report["reasons"]["NG"],
                         ["selected", "profile_exclusion", "profile_exclusion",
                          "inside_shape", "nonpositive_or_missing"])
        self.assertEqual(report["counts"]["saturated_ng_selected"], 0)
        for original, changed in zip(self.rows(), selected):
            self.assertEqual(original[:5], changed[:5])

    def test_run_exclusions_and_saturation_policy(self):
        profile = self.profile(
            apply_run_exclusions=True, saturated_ng="exclude",
            expected_counts={"rows": 5, "used_ng": 1, "used_rg": 3})
        intervals = {
            "FitInterval": [], "ExcludeRG": [],
            "ExcludeNG": [(0, 0, "unused", "unused")],
        }
        # Replace the placeholder with seconds parsed by the implementation.
        from cops_selection import parse_utc
        intervals["ExcludeNG"] = [(parse_utc("2016-09-30T10:00:20"),
                                    parse_utc("2016-09-30T10:00:20"), "a", "b")]
        selected, report = apply_selection_profile(
            self.rows(), self.columns, profile, intervals,
            [False, True, False, False, False])
        self.assertEqual([row[self.columns["useNG"]] for row in selected],
                         ["1", "0", "0", "0", "0"])
        self.assertEqual(report["reasons"]["NG"][1], "saturated_excluded")
        self.assertEqual(report["reasons"]["NG"][2], "run_exclusion")

    def test_count_change_is_fatal(self):
        profile = self.profile(expected_counts={"rows": 5, "used_ng": 99,
                                                "used_rg": 3})
        with self.assertRaisesRegex(SelectionError, "does not match the audited dataset"):
            apply_selection_profile(
                self.rows(), self.columns, profile,
                {"FitInterval": [], "ExcludeNG": [], "ExcludeRG": []},
                [False] * 5)

    def test_shape_dependent_inside_count_preserves_data_audit(self):
        profile = self.profile(
            expected_counts_before_inside={"used_ng": 4, "used_rg": 4})
        selected, report = apply_selection_profile(
            self.rows(), self.columns, profile,
            {"FitInterval": [], "ExcludeNG": [], "ExcludeRG": []},
            [False] * 5)
        self.assertEqual(report["counts"]["used_ng"], 3)
        self.assertEqual(report["counts"]["eligible_ng_before_inside"], 4)
        self.assertEqual(report["counts"]["eligible_ng_inside"], 1)

        # A refined shape can move the same valid point outside. The final
        # count changes, while the audited pressure/time/profile candidates do
        # not, so the profile must remain valid.
        refined = self.rows()
        refined[3][self.columns["inside"]] = "0"
        _, refined_report = apply_selection_profile(
            refined, self.columns, profile,
            {"FitInterval": [], "ExcludeNG": [], "ExcludeRG": []},
            [False] * 5)
        self.assertEqual(refined_report["counts"]["used_ng"], 4)
        self.assertEqual(refined_report["counts"]["eligible_ng_inside"], 0)

        # Altering an observation remains fatal; the dynamic rule is not a
        # blanket exemption for count changes.
        changed_data = self.rows()
        changed_data[2][self.columns["refNG"]] = "0"
        with self.assertRaisesRegex(SelectionError,
                                    "eligible_ng_before_inside=3"):
            apply_selection_profile(
                changed_data, self.columns, profile,
                {"FitInterval": [], "ExcludeNG": [], "ExcludeRG": []},
                [False] * 5)

    def test_checked_in_profiles_are_valid(self):
        expected = {
            "legacy_base_mask": ((839, 1146), (839, 1150)),
            "all_available": ((930, 1146), (934, 1150)),
            "transition_cut_100m": ((869, 1146), (873, 1150)),
        }
        for name, counts in expected.items():
            profile, identity = load_profile(name, APP / "fit_profiles")
            self.assertEqual((profile["expected_counts"]["used_ng"],
                              profile["expected_counts"]["used_rg"]), counts[0])
            self.assertEqual(
                (profile["expected_counts_before_inside"]["used_ng"],
                 profile["expected_counts_before_inside"]["used_rg"]),
                counts[1])
            self.assertEqual(identity["role"], "selection_profile")

    def test_runner_dry_run_uses_fourteen_jobs_for_every_profile(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = run_fit_profile.main([
                "--all-profiles", "--data-only", "--dry-run",
                "--out-prefix", "a1_dry_run"])
        self.assertEqual(status, 0)
        text = output.getvalue()
        self.assertEqual(text.count("--jobs 14"), 3)
        self.assertIn("data-only mode runs no optimizer", text)
        for name in ("legacy_base_mask", "all_available", "transition_cut_100m"):
            self.assertIn("--profile " + name, text)

    def test_runner_forwards_controlled_enhancement_contract(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = run_fit_profile.main([
                "--profile", "all_available", "--dry-run",
                "--out", "a4_dry_run", "--enhancement-degree", "8",
                "--fix-legacy-activity", "--fix-temperature", "--fix-kappa",
                "--lambda", "1e-5", "1e-4"])
        self.assertEqual(status, 0)
        text = output.getvalue()
        self.assertIn("--enhancement-degree 8", text)
        for option in ("--fix-legacy-activity", "--fix-temperature", "--fix-kappa"):
            self.assertIn(option, text)
        self.assertIn("at most 2 of the requested 14 fit process", text)

    def test_runner_rejects_uncontrolled_enhancement(self):
        error = io.StringIO()
        with contextlib.redirect_stderr(error), self.assertRaises(SystemExit):
            run_fit_profile.main([
                "--profile", "all_available", "--dry-run",
                "--enhancement-degree", "5"])
        self.assertIn("currently requires", error.getvalue())


if __name__ == "__main__":
    unittest.main()

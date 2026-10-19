"""Unit tests for Stage-A7 planning logic.

Wall-clock thresholds would make CI host-dependent, so these tests validate
the deterministic parts of the performance workflow: conservative resource
selection, campaign arithmetic, bounded names, and report rendering.  The
real-cache measurements are intentionally an operator-run benchmark documented
in A7_PERFORMANCE.md, not a pass/fail speed requirement.
"""

from pathlib import Path
import sys
import unittest


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

import benchmark_fit_performance as a7


class PerformancePlanningTests(unittest.TestCase):
    def scaling_rows(self):
        """A deterministic curve with a known memory/efficiency boundary."""
        return [
            {"jobs_effective": 1, "evaluations_per_second": 2.0,
             "parallel_efficiency": 1.0},
            {"jobs_effective": 2, "evaluations_per_second": 3.6,
             "parallel_efficiency": 0.9},
            {"jobs_effective": 4, "evaluations_per_second": 6.4,
             "parallel_efficiency": 0.8},
            {"jobs_effective": 8, "evaluations_per_second": 8.0,
             "parallel_efficiency": 0.5},
        ]

    def test_recommendation_enforces_efficiency_cpu_and_memory_bounds(self):
        # 70% of 10,000 MiB with a conservative 2,000 MiB per worker permits
        # only three jobs; degree 4 is therefore rejected despite good scaling.
        result = a7.choose_jobs(
            self.scaling_rows(), cores=16, available_mib=10000.0,
            worker_peak_mib=2000.0)
        self.assertEqual(result["jobs"], 2)
        self.assertEqual(result["threads_per_process"], 1)
        self.assertEqual(result["conservative_memory_capacity_jobs"], 3)

    def test_campaign_estimate_uses_declared_upper_scenario(self):
        degrees = [{
            "degree": 5,
            "timing_seconds": {"ready_to_evaluate": 10.0,
                               "residual_median": 0.2,
                               "jacobian_median": 0.8},
            "output_probe": {"total_bytes": 1000},
        }]
        recommendation = {"jobs": 2, "measured_parallel_efficiency": 0.5}
        result = a7.campaign_estimate(
            degrees, recommendation, profiles=3, lambdas=4, max_nfev=5)
        row = result["rows"][0]
        # startup = 3*10; optimizer capacity = 3*4*5*(.2+.8) = 60.
        self.assertEqual(row["serial_cpu_capacity_seconds"], 90.0)
        self.assertEqual(row["estimated_wall_seconds_at_recommended_jobs"], 90.0)
        self.assertEqual(row["conservative_output_bytes"], 12000)
        self.assertTrue(result["not_a_convergence_forecast"])

    def test_bounded_names_are_unique_and_reject_unsafe_ids(self):
        names = a7.bounded_output_names(
            ("legacy_base_mask", "all_available"), (5, 6), "trial_01")
        self.assertEqual(len(names), len(set(names)))
        self.assertTrue(all(name.startswith("a8_") for name in names))
        with self.assertRaises(ValueError):
            a7.bounded_output_names(("all_available",), (5,), "../escape")

    def test_markdown_contains_resource_and_scaling_tables(self):
        degree = {
            "degree": 5, "coefficient_count": 11,
            "timing_seconds": {"ready_to_evaluate": 2.0,
                               "basis_median": 0.1,
                               "residual_median": 0.2,
                               "jacobian_median": 0.3},
            "memory": {"peak_rss_mib": 400.0},
            "output_probe": {"total_bytes": 1024},
        }
        scaling = [{"jobs_effective": 1, "evaluations_per_second": 2.0,
                    "throughput_speedup": 1.0, "parallel_efficiency": 1.0,
                    "maximum_worker_rss_mib": 400.0}]
        report = {
            "created_utc": "2026-10-18T00:00:00+00:00",
            "base": {"path": "/data/base"},
            "configuration": {"profile": "all_available",
                              "threads_per_process": 1},
            "degrees": [degree], "parallel_scaling": scaling,
            "recommendation": {"jobs": 1, "threads_per_process": 1,
                               "basis": "test recommendation"},
            "campaign_estimate": {
                "assumption": "declared test formula",
                "total_serial_cpu_capacity_seconds": 3600.0,
                "total_estimated_wall_seconds_at_recommended_jobs": 3600.0,
                "total_conservative_output_bytes": 1024,
            },
            "bounded_output_names": ["a8_all_available_l5_trial"],
        }
        markdown = a7.render_markdown(report)
        self.assertIn("Measured resource table", markdown)
        self.assertIn("Parallel scaling", markdown)
        self.assertIn("not a convergence or science forecast", markdown)
        self.assertIn("a8_all_available_l5_trial", markdown)


if __name__ == "__main__":
    unittest.main()

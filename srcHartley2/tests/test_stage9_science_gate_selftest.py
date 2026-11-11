#!/usr/bin/env python3
"""Deterministic mutation tests for the Stage-9 independent-seed gate.

The fixture values are dimensionally plausible but have no cometary meaning.
They test parser and decision boundaries only.  A synthetic PASS must never be
presented as transported AMPS science evidence.
"""

from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable, Dict

import test_stage9_science_gate as stage9


SEEDS = (104729, 1299709, 15485863, 32452843)
ANALOG_FACTORS = (0.70, 0.90, 1.10, 1.30)
FORCED_FACTORS = (0.97, 0.99, 1.01, 1.03)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Archive:
    """Small, fully authenticated ANALOG/FORCED archive used by every test."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.manifest_path = root / "RUN_MANIFEST.json"
        self.manifest: Dict[str, Any] = {
            "schema": stage9.MANIFEST_SCHEMA,
            "policy": {
                "minimum_seeds": 4,
                "mpi_ranks": 12,
                "omp_threads": 1,
                "burn_in_iterations": 10000,
                "minimum_steady_windows": 3,
                "maximum_buffer_utilization": 0.75,
                # This is the two-sided 95% Student-t multiplier for four
                # seeds (three degrees of freedom), supplied by policy rather
                # than hidden in the evaluator.
                "confidence_multiplier": 3.182446305,
                "minimum_forced_daughter_collisions": 1,
                "flux_closure": {
                    "method": "independent_seed_ci_contains_zero",
                    "relative_tolerance": None,
                },
                "primary_species": {
                    "H2O": {
                        "relative_tolerance": 0.05,
                        "combined_se_multiplier": 2.0,
                    }
                },
                "daughter_species": {
                    "C2": {
                        "relative_tolerance": 0.10,
                        "combined_se_multiplier": 2.0,
                        "minimum_ess_fraction": 0.20,
                        "maximum_zero_cell_fraction": 0.05,
                        "minimum_variance_cost_gain": 2.0,
                        "minimum_gain_bin_fraction": 1.0,
                    }
                },
            },
            "runs": [],
        }
        for method in stage9.METHODS:
            factors = ANALOG_FACTORS if method == "ANALOG" else FORCED_FACTORS
            for seed, factor in zip(SEEDS, factors):
                summary = self.summary(method, seed, factor)
                relative = Path("runs") / method.lower() / str(seed) / "summary.json"
                path = root / relative
                write_json(path, summary)
                self.manifest["runs"].append(
                    {
                        "method": method,
                        "seed": seed,
                        "summary": str(relative),
                        "sha256": digest(path),
                    }
                )
        self.save_manifest()

    @staticmethod
    def profile(base: float, factor: float) -> list[dict[str, float]]:
        return [
            {"radius_m": 7000.0, "value": base * factor},
            {"radius_m": 17000.0, "value": 0.5 * base * factor},
            {"radius_m": 107000.0, "value": 0.1 * base * factor},
        ]

    def summary(self, method: str, seed: int, daughter_factor: float) -> Dict[str, Any]:
        # Parent seed scatter is deliberately identical across methods so the
        # invariance check is independent of the daughter variance test.
        parent_factor = 1.0 + 0.001 * (SEEDS.index(seed) - 1.5)
        return {
            "schema": stage9.SUMMARY_SCHEMA,
            "method": method,
            "seed": seed,
            "mpi_ranks": 12,
            "omp_threads": 1,
            "physical_input_sha256": "1" * 64,
            "mesh_sha256": "2" * 64,
            "reaction_network_sha256": "3" * 64,
            "started_from_iteration_zero": True,
            "restart_used": False,
            "convergence": {
                "status": "converged",
                "hard_limit_reached": False,
                "complete_steady_windows": 3,
                "discarded_burn_in_iterations": 10000,
            },
            "operations": {
                "wall_seconds": 100.0 if method == "ANALOG" else 150.0,
                "peak_buffer_utilization": 0.50,
                "buffer_exhaustions": 0,
                "emergency_events": 0,
                "unfinished_trajectories": 0,
            },
            "collision_audit": {
                "transported_forced_daughter_events": 0 if method == "ANALOG" else 10,
                "correction_mismatch_events": 0,
            },
            "profiles": {
                "H2O": self.profile(1.0e15, parent_factor),
                "C2": self.profile(1.0e9, daughter_factor),
            },
            "daughters": {
                "C2": {
                    "zero_cell_fraction": 0.0,
                    "model_carriers": 1000.0,
                    "sum_weights": 1000.0,
                    "sum_squared_weights": 1000.0,
                    "surface_source_s-1": 0.0,
                    "chemical_production_s-1": 100.0,
                    "inner_reemitted_s-1": 0.0,
                    "chemical_destruction_s-1": 40.0,
                    "inner_incident_s-1": 0.0,
                    "outer_escape_s-1": 60.0,
                    "net_flux_s-1": 0.0,
                }
            },
        }

    def save_manifest(self) -> None:
        write_json(self.manifest_path, self.manifest)

    def mutate_summary(
        self, method: str, seed: int, mutation: Callable[[Dict[str, Any]], None], *, rehash: bool = True
    ) -> None:
        record = next(
            row for row in self.manifest["runs"] if row["method"] == method and row["seed"] == seed
        )
        path = self.root / record["summary"]
        summary = json.loads(path.read_text(encoding="utf-8"))
        mutation(summary)
        write_json(path, summary)
        if rehash:
            record["sha256"] = digest(path)
            self.save_manifest()


class Stage9ScienceGateTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.archive = Archive(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_complete_independent_seed_matrix_passes(self) -> None:
        result = stage9.evaluate(self.archive.manifest_path)
        self.assertEqual(result["outcome"], "PASS")
        self.assertEqual(result["run_count"], 8)
        self.assertEqual(result["seeds"], list(SEEDS))
        self.assertGreater(result["efficiency"]["C2"]["accepted_bin_fraction"], 0.99)
        for method in stage9.METHODS:
            self.assertEqual(len(result["profiles"]["C2"][method]), 3)
            self.assertIn("confidence_low", result["profiles"]["C2"][method][0])

    def test_missing_metric_is_invalid_evidence(self) -> None:
        self.archive.mutate_summary(
            "FORCED", SEEDS[0], lambda value: value["operations"].pop("wall_seconds")
        )
        with self.assertRaisesRegex(stage9.EvidenceError, "keys differ"):
            stage9.evaluate(self.archive.manifest_path)

    def test_summary_hash_protects_post_run_evidence(self) -> None:
        self.archive.mutate_summary(
            "FORCED", SEEDS[0], lambda value: value["operations"].__setitem__("wall_seconds", 1.0), rehash=False
        )
        with self.assertRaisesRegex(stage9.EvidenceError, "hash mismatch"):
            stage9.evaluate(self.archive.manifest_path)

    def test_seed_matrix_and_physical_provenance_are_strict(self) -> None:
        self.archive.manifest["runs"].pop()
        self.archive.save_manifest()
        with self.assertRaisesRegex(stage9.EvidenceError, "at least 8"):
            stage9.evaluate(self.archive.manifest_path)

        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory()
        self.archive = Archive(Path(self.temporary.name))
        self.archive.mutate_summary(
            "FORCED", SEEDS[0], lambda value: value.__setitem__("mesh_sha256", "4" * 64)
        )
        with self.assertRaisesRegex(stage9.EvidenceError, "one mesh_sha256"):
            stage9.evaluate(self.archive.manifest_path)

    def test_nonconvergence_and_hard_limit_fail_science(self) -> None:
        def mutation(value: Dict[str, Any]) -> None:
            value["convergence"]["status"] = "nonconverged"
            value["convergence"]["hard_limit_reached"] = True

        self.archive.mutate_summary("ANALOG", SEEDS[0], mutation)
        result = stage9.evaluate(self.archive.manifest_path)
        self.assertEqual(result["outcome"], "FAIL_SCIENCE")
        self.assertIn("hard limit", "\n".join(result["failures"]))

    def test_restart_or_nonzero_start_is_invalid_evidence(self) -> None:
        self.archive.mutate_summary(
            "ANALOG", SEEDS[0], lambda value: value.__setitem__("restart_used", True)
        )
        with self.assertRaisesRegex(stage9.EvidenceError, "restart"):
            stage9.evaluate(self.archive.manifest_path)

    def test_buffer_emergency_and_unfinished_trajectory_fail(self) -> None:
        def mutation(value: Dict[str, Any]) -> None:
            value["operations"].update(
                peak_buffer_utilization=0.80,
                buffer_exhaustions=1,
                emergency_events=1,
                unfinished_trajectories=1,
            )

        self.archive.mutate_summary("FORCED", SEEDS[0], mutation)
        failures = "\n".join(stage9.evaluate(self.archive.manifest_path)["failures"])
        self.assertIn("buffer utilization", failures)
        self.assertIn("buffer exhaustion", failures)
        self.assertIn("emergency adaptation", failures)
        self.assertIn("unfinished", failures)

    def test_primary_profile_change_fails_invariance(self) -> None:
        for seed in SEEDS:
            self.archive.mutate_summary(
                "FORCED",
                seed,
                lambda value: [row.__setitem__("value", 1.25 * row["value"]) for row in value["profiles"]["H2O"]],
            )
        failures = "\n".join(stage9.evaluate(self.archive.manifest_path)["failures"])
        self.assertIn("H2O profile differs", failures)

    def test_flux_identity_corruption_is_invalid_evidence(self) -> None:
        self.archive.mutate_summary(
            "FORCED",
            SEEDS[0],
            lambda value: value["daughters"]["C2"].__setitem__("net_flux_s-1", 5.0),
        )
        with self.assertRaisesRegex(stage9.EvidenceError, "flux identity"):
            stage9.evaluate(self.archive.manifest_path)

    def test_flux_residual_ci_ess_and_zero_fraction_fail_science(self) -> None:
        def quality_mutation(value: Dict[str, Any]) -> None:
            row = value["daughters"]["C2"]
            row["sum_squared_weights"] = 10000.0
            row["zero_cell_fraction"] = 0.10

        self.archive.mutate_summary("FORCED", SEEDS[0], quality_mutation)
        for seed in SEEDS:
            def biased_flux(value: Dict[str, Any]) -> None:
                row = value["daughters"]["C2"]
                row["outer_escape_s-1"] = 40.0
                row["net_flux_s-1"] = 20.0

            self.archive.mutate_summary("FORCED", seed, biased_flux)
        failures = "\n".join(stage9.evaluate(self.archive.manifest_path)["failures"])
        self.assertIn("flux-residual CI excludes zero", failures)
        self.assertIn("ESS fraction", failures)
        self.assertIn("zero-cell fraction", failures)

    def test_relative_flux_tolerance_method_passes_and_fails_per_run(self) -> None:
        self.archive.manifest["policy"]["flux_closure"] = {
            "method": "relative_tolerance",
            "relative_tolerance": 0.05,
        }
        self.archive.save_manifest()
        result = stage9.evaluate(self.archive.manifest_path)
        self.assertEqual(result["outcome"], "PASS")
        self.assertEqual(
            result["flux_closure"]["C2"]["FORCED"]["policy_method"],
            "relative_tolerance",
        )

        def biased_flux(value: Dict[str, Any]) -> None:
            row = value["daughters"]["C2"]
            row["outer_escape_s-1"] = 40.0
            row["net_flux_s-1"] = 20.0

        self.archive.mutate_summary("FORCED", SEEDS[0], biased_flux)
        failures = "\n".join(stage9.evaluate(self.archive.manifest_path)["failures"])
        self.assertIn("exceeds the flux relative-tolerance policy", failures)

    def test_flux_policy_is_explicit_and_strict(self) -> None:
        policy = self.archive.manifest["policy"]
        policy["flux_closure"] = {
            "method": "DECISION_REQUIRED",
            "relative_tolerance": "DECISION_REQUIRED",
        }
        self.archive.save_manifest()
        with self.assertRaisesRegex(stage9.EvidenceError, "unknown flux_closure method"):
            stage9.evaluate(self.archive.manifest_path)

        policy["flux_closure"] = {
            "method": "independent_seed_ci_contains_zero",
            "relative_tolerance": 0.05,
        }
        self.archive.save_manifest()
        with self.assertRaisesRegex(stage9.EvidenceError, "must be null"):
            stage9.evaluate(self.archive.manifest_path)

        policy["flux_closure"] = {
            "method": "relative_tolerance",
            "relative_tolerance": 0.0,
        }
        self.archive.save_manifest()
        with self.assertRaisesRegex(stage9.EvidenceError, "must be positive"):
            stage9.evaluate(self.archive.manifest_path)

    def test_transported_collision_evidence_is_required(self) -> None:
        self.archive.mutate_summary(
            "FORCED",
            SEEDS[0],
            lambda value: value["collision_audit"].__setitem__(
                "transported_forced_daughter_events", 0
            ),
        )
        failures = "\n".join(stage9.evaluate(self.archive.manifest_path)["failures"])
        self.assertIn("transported forced-daughter collision", failures)

    def test_collision_correction_mismatch_is_rejected(self) -> None:
        self.archive.mutate_summary(
            "FORCED",
            SEEDS[0],
            lambda value: value["collision_audit"].__setitem__("correction_mismatch_events", 1),
        )
        failures = "\n".join(stage9.evaluate(self.archive.manifest_path)["failures"])
        self.assertIn("changed a carrier correction", failures)

    def test_no_variance_cost_improvement_fails(self) -> None:
        for seed, factor in zip(SEEDS, ANALOG_FACTORS):
            self.archive.mutate_summary(
                "FORCED",
                seed,
                lambda value, factor=factor: value["profiles"].__setitem__(
                    "C2", self.archive.profile(1.0e9, factor)
                ),
            )
        result = stage9.evaluate(self.archive.manifest_path)
        self.assertEqual(result["outcome"], "FAIL_SCIENCE")
        self.assertIn("variance-per-cost gain", "\n".join(result["failures"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)

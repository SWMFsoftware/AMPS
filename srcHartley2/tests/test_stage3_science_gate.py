#!/usr/bin/env python3
"""Synthetic regression tests for the permanent Stage-3 science gate.

These tests deliberately construct a complete eight-run archive rather than
mocking the evaluator's file readers.  That makes the test enforce the same
manifest, sampling-report, exact-cell, block-policy, and Rubin-baseline
contracts that a production archive must satisfy, while keeping the fixture
small enough to run without AMPS or MPI.

The synthetic densities have no physical meaning.  They are positive,
well-conditioned inputs used only to exercise decision boundaries and outcome
classification.  Scientific evidence must always come from the archived
12-rank Halley runs, never from this fixture.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable, Dict, List


ROOT = Path(__file__).resolve().parents[2]
VALIDATION = ROOT / "srcHartley2" / "validation"
VISUALIZATION = ROOT / "srcHartley2" / "visualization"
sys.path.insert(0, str(VALIDATION))
sys.path.insert(0, str(VISUALIZATION))

from adaptive_radial_bins import Cell  # noqa: E402
from evaluate_stage3_science_gate import (  # noqa: E402
    _population_coverage,
    _support_statistics,
    derive_common_support,
    evaluate,
    main as evaluator_main,
)
from stage3_gate_common import (  # noqa: E402
    GateError,
    MANIFEST_SCHEMA,
    POLICY_SCHEMA,
    load_policy,
    sha256_file,
    sha256_json,
    validate_policy,
)


POLICY_SOURCE = VALIDATION / "stage3_science_policy.json"


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


class SyntheticArchive:
    """Create and mutate a minimal, structurally complete science archive."""

    cell_count = 20
    block_cell_count = 2

    def __init__(self, root: Path) -> None:
        self.root = root
        self.archive = root / "archive"
        self.policy = json.loads(POLICY_SOURCE.read_text(encoding="utf-8"))
        self.policy_path = root / "policy.json"
        write_json(self.policy_path, self.policy)

        inner = self.policy["nucleus_radius_m"] + 1000.0 * self.policy[
            "comparison_altitude_km"
        ][0]
        outer = self.policy["nucleus_radius_m"] + 1000.0 * self.policy[
            "comparison_altitude_km"
        ][1]
        logarithmic_step = math.log(outer / inner) / self.cell_count
        self.edges = [inner * math.exp(index * logarithmic_step) for index in range(self.cell_count + 1)]
        # Pin the endpoints exactly.  The gate intentionally selects only
        # fully contained cells, so a rounded endpoint must not accidentally
        # turn the last synthetic cell into an excluded partial cell.
        self.edges[0] = inner
        self.edges[-1] = outer
        self.run_roots: Dict[tuple[str, int], Path] = {}
        self._write_archive()

    def _write_archive(self) -> None:
        runs: List[Dict[str, Any]] = []
        frozen = self.archive / "frozen"
        frozen.mkdir(parents=True, exist_ok=True)
        executable = frozen / "amps.stage3-science"
        input_path = frozen / "halley1d.input"
        pic_global = frozen / "picGlobal.dfn"
        executable.write_text("synthetic executable\n", encoding="utf-8")
        input_path.write_text("synthetic input\n", encoding="utf-8")
        pic_global.write_text("synthetic picGlobal\n", encoding="utf-8")
        write_json(self.archive / "policy.json", self.policy)
        executable_hash = sha256_file(executable)
        input_hash = sha256_file(input_path)
        pic_global_hash = sha256_file(pic_global)
        for mode in self.policy["modes"]:
            for seed in self.policy["seeds"]:
                relative = Path("runs") / mode.lower() / str(seed)
                run_root = self.archive / relative
                self.run_roots[(mode, seed)] = run_root
                diagnostics = run_root / "PT" / "diagnostics"
                diagnostics.mkdir(parents=True, exist_ok=True)
                # The evaluator audits the text log for a final sampling-
                # population request that can be raised after the last JSON
                # report.  A quiet synthetic run therefore needs an explicit
                # empty log rather than an absent optional artifact.
                (run_root / "run.log").write_text("", encoding="utf-8")
                (run_root / "error.log").write_text("", encoding="utf-8")
                self._write_measures(diagnostics / "cell_measures.csv")
                for output_index in range(4):
                    self._write_report(mode, seed, output_index)
                    (diagnostics / f"comet_stage0.out={output_index}.csv").write_text(
                        "record_type\nmetadata\n", encoding="utf-8"
                    )
                    for species_index, species in enumerate(self.policy["important_species"]):
                        self._write_cells(
                            mode,
                            seed,
                            species,
                            species_index,
                            output_index,
                        )
                artifact_hash = self._write_artifact_manifest(run_root)
                runs.append(
                    {
                        "mode": mode,
                        "seed": seed,
                        "status": "completed",
                        "return_code": 0,
                        "mpi_ranks": self.policy["mpi_ranks"],
                        "executable_sha256": executable_hash,
                        "input_sha256": input_hash,
                        "pic_global_sha256": pic_global_hash,
                        "artifact_manifest_sha256": artifact_hash,
                        "attempt_path": str(relative),
                    }
                )
        write_json(
            self.archive / "RUN_MANIFEST.json",
            {
                "schema": MANIFEST_SCHEMA,
                "status": "completed",
                "policy_sha256": sha256_json(self.policy),
                "policy_path": "policy.json",
                "frozen_executable": "frozen/amps.stage3-science",
                "frozen_executable_sha256": executable_hash,
                "frozen_input": "frozen/halley1d.input",
                "frozen_input_sha256": input_hash,
                "frozen_pic_global": "frozen/picGlobal.dfn",
                "frozen_pic_global_sha256": pic_global_hash,
                "runs": runs,
            },
        )

    def _write_artifact_manifest(self, run_root: Path) -> str:
        diagnostics = run_root / "PT" / "diagnostics"
        artifacts = [
            run_root / "run.log",
            run_root / "error.log",
            diagnostics / "cell_measures.csv",
        ]
        for output_index in range(4):
            artifacts.extend(
                [
                    diagnostics / f"comet_stage0.out={output_index}.json",
                    diagnostics / f"comet_stage0.out={output_index}.csv",
                ]
            )
            artifacts.extend(
                diagnostics / f"stage3_cells.{species}.out={output_index}.csv"
                for species in self.policy["important_species"]
            )
        path = run_root / "artifacts.sha256"
        path.write_text(
            "".join(
                f"{sha256_file(artifact)}  {artifact.relative_to(run_root)}\n"
                for artifact in sorted(artifacts)
            ),
            encoding="utf-8",
        )
        return sha256_file(path)

    def _refresh_artifacts(self, mode: str, seed: int) -> None:
        """Refresh synthetic provenance after an intentional semantic edit."""
        artifact_hash = self._write_artifact_manifest(self.run_root(mode, seed))
        manifest_path = self.archive / "RUN_MANIFEST.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for run in manifest["runs"]:
            if run["mode"] == mode and run["seed"] == seed:
                run["artifact_manifest_sha256"] = artifact_hash
                break
        else:
            raise AssertionError("synthetic run record is missing")
        write_json(manifest_path, manifest)

    def _write_measures(self, path: Path) -> None:
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                ["cell_index", "radius_m", "r_inner_m", "r_outer_m", "measure_m3"]
            )
            for index, (inner, outer) in enumerate(zip(self.edges, self.edges[1:])):
                writer.writerow([index, math.sqrt(inner * outer), inner, outer, 1.0])

    def _write_report(self, mode: str, seed: int, output_index: int) -> None:
        flux = []
        for species in self.policy["important_species"]:
            flux.append(
                {
                    "species": species,
                    "configured_surface_source_s-1": 0.0,
                    "sampled_surface_source_s-1": 0.0,
                    "production_agreement": True,
                    "relative_l1_change": None if output_index == 0 else 0.01,
                }
            )
        report = {
            "schema": "amps.comet.stage0-diagnostics.v1",
            "iteration": (output_index + 1) * self.policy["sampling_length"],
            "output_index": output_index,
            "sampling_length": self.policy["sampling_length"],
            "random_seed": seed,
            "random_seed_configured": True,
            "mpi_ranks": self.policy["mpi_ranks"],
            "residual_sampling": mode,
            "emergency_requests": 0,
            "emergency_events": {
                "total": 0,
                "buffer_utilization": 0,
                "sampling_population": 0,
                "allocation_burst": 0,
            },
            "reactions": [],
            "species_flux": flux,
        }
        write_json(self.report_path(mode, seed, output_index), report)

    def _write_cells(
        self,
        mode: str,
        seed: int,
        species: str,
        species_index: int,
        output_index: int,
    ) -> None:
        path = self.table_path(mode, seed, species, output_index)
        fieldnames = [
            "cell_index",
            "block_index",
            "local_cell_index",
            "block_cell_count",
            "radius_m",
            "r_inner_m",
            "r_outer_m",
            "measure_m3",
            "number_density_m-3",
            "mean_model_particles",
            "population_target_min_per_cell",
            "population_target_max_per_cell",
            "population_locked",
            "population_class",
            "population_policy_source",
            "sampling_length",
            "output_index",
            "random_seed",
            "random_seed_configured",
            "residual_sampling",
            "mpi_ranks",
        ]
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            for index, (inner, outer) in enumerate(zip(self.edges, self.edges[1:])):
                writer.writerow(
                    {
                        "cell_index": index,
                        "block_index": index // self.block_cell_count,
                        "local_cell_index": index % self.block_cell_count,
                        "block_cell_count": self.block_cell_count,
                        "radius_m": math.sqrt(inner * outer),
                        "r_inner_m": inner,
                        "r_outer_m": outer,
                        "measure_m3": 1.0,
                        "number_density_m-3": 1.0e12 * (species_index + 1),
                        "mean_model_particles": 100.0,
                        "population_target_min_per_cell": 50,
                        "population_target_max_per_cell": 200,
                        "population_locked": 0,
                        # These integer values are the authoritative IDs emitted
                        # by the current C++ schema for CHEMICAL_ONLY policy.
                        "population_class": 2,
                        "population_policy_source": 2,
                        "sampling_length": self.policy["sampling_length"],
                        "output_index": output_index,
                        "random_seed": seed,
                        "random_seed_configured": 1,
                        "residual_sampling": mode,
                        "mpi_ranks": self.policy["mpi_ranks"],
                    }
                )

    def run_root(self, mode: str, seed: int) -> Path:
        return self.run_roots[(mode, seed)]

    def report_path(self, mode: str, seed: int, output_index: int) -> Path:
        return (
            self.run_root(mode, seed)
            / "PT"
            / "diagnostics"
            / f"comet_stage0.out={output_index}.json"
        )

    def table_path(self, mode: str, seed: int, species: str, output_index: int) -> Path:
        return (
            self.run_root(mode, seed)
            / "PT"
            / "diagnostics"
            / f"stage3_cells.{species}.out={output_index}.csv"
        )

    def mutate_report(
        self,
        mode: str,
        seed: int,
        output_index: int,
        mutation: Callable[[Dict[str, Any]], None],
    ) -> None:
        path = self.report_path(mode, seed, output_index)
        report = json.loads(path.read_text(encoding="utf-8"))
        mutation(report)
        write_json(path, report)
        self._refresh_artifacts(mode, seed)

    def mutate_table(
        self,
        mode: str,
        seed: int,
        species: str,
        output_index: int,
        mutation: Callable[[int, Dict[str, str]], None],
    ) -> None:
        path = self.table_path(mode, seed, species, output_index)
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
        for index, row in enumerate(rows):
            mutation(index, row)
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        self._refresh_artifacts(mode, seed)

    def accepted_baseline(self) -> Path:
        """Promote the fixture's generated candidate, mirroring user review."""
        first = evaluate(self.archive, self.policy_path, None)
        self._assert_baseline_required(first)
        candidate = json.loads(
            (self.archive / "evaluation" / "rubin_candidate_baseline.json").read_text(
                encoding="utf-8"
            )
        )
        candidate["accepted"] = True
        path = self.root / "accepted-rubin-baseline.json"
        write_json(path, candidate)
        return path

    @staticmethod
    def _assert_baseline_required(summary: Dict[str, Any]) -> None:
        if summary["outcome"] != "INCONCLUSIVE_RERUN" or summary["rubin_status"] != (
            "BASELINE_REQUIRED"
        ):
            raise AssertionError("synthetic candidate did not request baseline acceptance")


def simple_cell(index: int, density: float, occupancy: float) -> Cell:
    """Construct one exact cell for pure support/statistics boundary tests."""
    return Cell(
        index=index,
        radius_m=index + 1.5,
        inner_m=index + 1.0,
        outer_m=index + 2.0,
        volume_m3=1.0,
        density_m3=density,
        occupancy=occupancy,
        block_index=index // 2,
        local_cell_index=index % 2,
        block_cell_count=2,
        population_target_min=50,
        population_target_max=200,
        population_locked=False,
        population_class=2,
        population_policy_source=2,
    )


class PolicyAndStatisticBoundaryTest(unittest.TestCase):
    def test_policy_accepts_only_the_frozen_approved_values(self) -> None:
        policy = json.loads(POLICY_SOURCE.read_text(encoding="utf-8"))
        self.assertEqual(validate_policy(policy)["schema"], POLICY_SCHEMA)

        mutations = {
            "minimum exposure": ("minimum_particle_iteration_exposure", 99999.0),
            "CV target": ("systematic_cv_target", 0.26),
            "CV fraction": ("systematic_cv_required_bin_fraction", 0.89),
            "CV maximum": ("systematic_cv_hard_maximum", 0.51),
            "zero fraction": ("systematic_raw_zero_cell_fraction_maximum", 0.06),
            "population fraction": (
                "systematic_population_required_block_fraction",
                0.89,
            ),
            "mode tolerance": ("mode_mean_relative_tolerance", 0.06),
            "mode SE": ("mode_mean_combined_standard_error_multiplier", 3.0),
        }
        for label, (key, value) in mutations.items():
            with self.subTest(label=label):
                changed = dict(policy)
                changed[key] = value
                with self.assertRaises(GateError):
                    validate_policy(changed)

    def test_common_support_uses_the_minimum_across_all_eight_runs(self) -> None:
        fast = [simple_cell(index, 1.0, 10.0) for index in range(3)]
        slow = [simple_cell(index, 1.0, 5.0) for index in range(3)]
        support = derive_common_support(
            [fast for _ in range(7)] + [slow],
            sampling_length=10000,
            minimum_exposure=100000.0,
        )
        self.assertEqual([(row["first_cell"], row["last_cell"]) for row in support], [(0, 1), (2, 2)])
        self.assertEqual(support[0]["minimum_run_exposure"], 100000.0)
        self.assertTrue(support[0]["qualified"])
        self.assertFalse(support[1]["qualified"])

    def test_cv_and_population_thresholds_are_inclusive(self) -> None:
        mean = 100.0
        half_span = 0.25 * mean * math.sqrt(3.0) / 2.0
        runs = [
            [simple_cell(0, value, 100.0)]
            for value in (mean - half_span, mean - half_span, mean + half_span, mean + half_span)
        ]
        support = derive_common_support(runs, 10000, 100000.0)
        row = _support_statistics(runs, support, 10000)[0]
        self.assertTrue(math.isclose(row["coefficient_of_variation"], 0.25, rel_tol=1e-14))

        # Nine of ten populated blocks at the exact 50--200 boundary must pass
        # the inclusive 90% rule.  The tenth block is deliberately outside.
        cells: List[Cell] = []
        for block in range(10):
            occupancy = 50.0 if block < 9 else 40.0
            cells.extend(
                [simple_cell(2 * block, 1.0, occupancy), simple_cell(2 * block + 1, 1.0, occupancy)]
            )
        coverage = _population_coverage(cells)
        self.assertEqual(coverage["populated_blocks"], 10)
        self.assertEqual(coverage["blocks_within_target"], 9)
        self.assertEqual(coverage["within_target_fraction"], 0.9)

        malformed = [simple_cell(0, 1.0, 100.0), simple_cell(1, 1.0, 100.0)]
        malformed[1] = Cell(
            **{**malformed[1].__dict__, "local_cell_index": 0}
        )
        with self.assertRaises(GateError):
            _population_coverage(malformed)


class CompleteArchiveGateTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = SyntheticArchive(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_pass_with_explicitly_accepted_rubin_baseline(self) -> None:
        baseline = self.fixture.accepted_baseline()
        summary = evaluate(
            self.fixture.archive, self.fixture.policy_path, baseline
        )
        self.assertEqual(summary["outcome"], "PASS_SKIP_STAGE4")
        self.assertEqual(summary["rubin_status"], "PASS")
        self.assertFalse(summary["inconclusive_reasons"])
        self.assertFalse(summary["quality_failures"])

    def test_missing_baseline_is_inconclusive_and_emits_candidate(self) -> None:
        summary = evaluate(self.fixture.archive, self.fixture.policy_path, None)
        self.assertEqual(summary["outcome"], "INCONCLUSIVE_RERUN")
        self.assertEqual(summary["rubin_status"], "BASELINE_REQUIRED")
        self.assertTrue(
            (self.fixture.archive / "evaluation" / "rubin_candidate_baseline.json").is_file()
        )
        self.assertTrue(
            any("baseline" in reason for reason in summary["inconclusive_reasons"])
        )

    def test_exact_zero_and_mode_mean_thresholds_pass(self) -> None:
        # The Bernoulli mean is exactly 5% below the systematic mean.  Zero
        # within-mode variance makes the relative branch of the agreement rule
        # authoritative, and equality must pass.
        for seed in self.fixture.policy["seeds"]:
            for species in self.fixture.policy["important_species"]:
                self.fixture.mutate_table(
                    "BERNOULLI",
                    seed,
                    species,
                    3,
                    lambda _index, row: row.__setitem__(
                        "number_density_m-3",
                        str(0.95 * float(row["number_density_m-3"])),
                    ),
                )

        # One zero among 20 raw systematic cells is exactly the approved 5%.
        # Merge that low-exposure cell with its neighbor so the adaptive-bin CV
        # remains below the 50% hard ceiling.  The two-cell block mean is
        # exactly 50, testing the inclusive population lower bound as well.
        target_seed = self.fixture.policy["seeds"][0]
        for mode in self.fixture.policy["modes"]:
            for seed in self.fixture.policy["seeds"]:
                self.fixture.mutate_table(
                    mode,
                    seed,
                    "C2",
                    3,
                    lambda index, row: row.__setitem__(
                        "mean_model_particles", "5" if index == 0 else ("95" if index == 1 else row["mean_model_particles"])
                    ),
                )
        self.fixture.mutate_table(
            "SYSTEMATIC",
            target_seed,
            "C2",
            3,
            lambda index, row: row.__setitem__(
                "number_density_m-3", "0" if index == 0 else row["number_density_m-3"]
            ),
        )

        baseline = self.fixture.accepted_baseline()
        summary = evaluate(self.fixture.archive, self.fixture.policy_path, baseline)
        self.assertEqual(summary["outcome"], "PASS_SKIP_STAGE4")
        zero_values = summary["species"]["C2"]["raw_zero_cell_fraction"]["SYSTEMATIC"]
        self.assertIn(0.05, zero_values)
        self.assertTrue(
            all(item["passed"] for result in summary["species"].values() for item in result["mode_agreement"])
        )

    def test_cv_fraction_and_hard_maximum_equalities_pass(self) -> None:
        # Two of twenty bins are deliberately above the 25% target, leaving
        # exactly 90% at or below it.  One exceptional bin is exactly at the
        # 50% hard maximum.  The symmetric four-seed values preserve the mean,
        # so this test isolates inclusive CV comparisons from Rubin and
        # Bernoulli/systematic mean agreement.
        seeds = self.fixture.policy["seeds"]
        for cell_index, requested_cv in ((0, 0.50), (1, 0.30)):
            half_span = requested_cv * math.sqrt(3.0) / 2.0
            factors = (1.0 - half_span, 1.0 - half_span, 1.0 + half_span, 1.0 + half_span)
            for seed, factor in zip(seeds, factors):
                self.fixture.mutate_table(
                    "SYSTEMATIC",
                    seed,
                    "C2",
                    3,
                    lambda index, row, factor=factor, target=cell_index: row.__setitem__(
                        "number_density_m-3",
                        str(factor * float(row["number_density_m-3"]))
                        if index == target
                        else row["number_density_m-3"],
                    ),
                )

        baseline = self.fixture.accepted_baseline()
        summary = evaluate(self.fixture.archive, self.fixture.policy_path, baseline)
        self.assertEqual(summary["outcome"], "PASS_SKIP_STAGE4")
        rows = summary["species"]["C2"]["mode_statistics"]["SYSTEMATIC"]
        self.assertEqual(sum(row["coefficient_of_variation"] <= 0.25 for row in rows), 18)
        self.assertTrue(math.isclose(max(row["coefficient_of_variation"] for row in rows), 0.50))

    def test_cv_zero_and_population_violations_are_quality_failures(self) -> None:
        seeds = self.fixture.policy["seeds"]
        factors = (0.2, 0.6, 1.4, 1.8)
        for seed, factor in zip(seeds, factors):
            self.fixture.mutate_table(
                "SYSTEMATIC",
                seed,
                "C2",
                3,
                lambda _index, row, factor=factor: row.__setitem__(
                    "number_density_m-3",
                    str(factor * float(row["number_density_m-3"])),
                ),
            )

        # Two raw zero cells are 10%, strictly beyond the approved 5%.
        self.fixture.mutate_table(
            "SYSTEMATIC",
            seeds[0],
            "CN",
            3,
            lambda index, row: row.__setitem__(
                "number_density_m-3", "0" if index < 2 else row["number_density_m-3"]
            ),
        )

        # Two of ten populated blocks below 50 gives only 80% coverage.
        for seed in seeds:
            for output_index in (1, 2, 3):
                self.fixture.mutate_table(
                    "SYSTEMATIC",
                    seed,
                    "HCO",
                    output_index,
                    lambda index, row: row.__setitem__(
                        "mean_model_particles", "40" if index < 4 else row["mean_model_particles"]
                    ),
                )

        baseline = self.fixture.accepted_baseline()
        summary = evaluate(self.fixture.archive, self.fixture.policy_path, baseline)
        self.assertEqual(summary["outcome"], "FAIL_CONSIDER_STAGE4")
        failures = "\n".join(summary["quality_failures"])
        self.assertIn("C2 systematic CV", failures)
        self.assertIn("CN seed", failures)
        self.assertIn("raw zero fraction", failures)
        self.assertIn("HCO seed", failures)
        self.assertIn("population coverage", failures)

    def test_emergency_or_nonconvergence_is_inconclusive(self) -> None:
        seed = self.fixture.policy["seeds"][0]
        self.fixture.mutate_report(
            "SYSTEMATIC",
            seed,
            2,
            lambda report: report.__setitem__(
                "emergency_events",
                {
                    "total": 1,
                    "buffer_utilization": 1,
                    "sampling_population": 0,
                    "allocation_burst": 0,
                },
            ),
        )

        def break_convergence(report: Dict[str, Any]) -> None:
            report["species_flux"][0]["relative_l1_change"] = 0.0500001

        self.fixture.mutate_report("BERNOULLI", seed, 3, break_convergence)
        baseline = self.fixture.accepted_baseline()
        summary = evaluate(self.fixture.archive, self.fixture.policy_path, baseline)
        self.assertEqual(summary["outcome"], "INCONCLUSIVE_RERUN")
        reasons = "\n".join(summary["inconclusive_reasons"])
        self.assertIn("emergency adaptation", reasons)
        self.assertIn("maximum profile change", reasons)

    def test_invalid_block_policy_is_invalid_evidence(self) -> None:
        seed = self.fixture.policy["seeds"][0]
        self.fixture.mutate_table(
            "SYSTEMATIC",
            seed,
            "C2",
            3,
            lambda _index, row: row.__setitem__(
                "population_target_min_per_cell", "100"
            ),
        )
        exit_code = evaluator_main(
            [str(self.fixture.archive), "--policy", str(self.fixture.policy_path)]
        )
        self.assertEqual(exit_code, 12)
        report = json.loads(
            (self.fixture.archive / "evaluation" / "stage3_science_gate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(report["outcome"], "INVALID_EVIDENCE")
        self.assertIn("population policy", report["error"])

    def test_invalid_manifest_is_invalid_evidence(self) -> None:
        manifest_path = self.fixture.archive / "RUN_MANIFEST.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["runs"][0]["mpi_ranks"] = 4
        write_json(manifest_path, manifest)
        self.assertEqual(
            evaluator_main(
                [str(self.fixture.archive), "--policy", str(self.fixture.policy_path)]
            ),
            12,
        )
        report = json.loads(
            (self.fixture.archive / "evaluation" / "stage3_science_gate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn("MPI rank", report["error"])

    def test_post_run_artifact_mutation_is_invalid_evidence(self) -> None:
        seed = self.fixture.policy["seeds"][0]
        path = self.fixture.table_path("SYSTEMATIC", seed, "C2", 3)
        # Do not refresh the synthetic checksum: this represents damage or an
        # edit after the operational runner declared the case complete.
        path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        self.assertEqual(
            evaluator_main(
                [str(self.fixture.archive), "--policy", str(self.fixture.policy_path)]
            ),
            12,
        )
        report = json.loads(
            (self.fixture.archive / "evaluation" / "stage3_science_gate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn("artifact changed", report["error"])

    def test_malformed_cell_schema_is_reported_not_raised(self) -> None:
        seed = self.fixture.policy["seeds"][0]
        path = self.fixture.table_path("SYSTEMATIC", seed, "C2", 3)
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            rows = list(reader)
            fieldnames = [name for name in (reader.fieldnames or []) if name != "block_index"]
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        self.assertEqual(
            evaluator_main(
                [str(self.fixture.archive), "--policy", str(self.fixture.policy_path)]
            ),
            12,
        )

    def test_unaccepted_or_wrong_policy_baseline_is_invalid(self) -> None:
        baseline = self.fixture.accepted_baseline()
        record = json.loads(baseline.read_text(encoding="utf-8"))
        record["accepted"] = False
        write_json(baseline, record)
        with self.assertRaises(GateError):
            evaluate(self.fixture.archive, self.fixture.policy_path, baseline)

        record["accepted"] = True
        record["policy_sha256"] = "0" * 64
        write_json(baseline, record)
        with self.assertRaises(GateError):
            evaluate(self.fixture.archive, self.fixture.policy_path, baseline)


if __name__ == "__main__":
    unittest.main(verbosity=2)

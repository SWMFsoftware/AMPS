#!/usr/bin/env python3
"""Dependency-free synthetic tests for Stage-9 campaign/summary tooling.

The fixture is operational-contract evidence only.  Its profiles and fluxes
are synthetic and must never be cited as comet science results.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict

import stage9_run_summary as builder
import run_stage9_science_campaign as campaign

TEST_DIR = Path(__file__).resolve().parents[1] / "tests"
sys.path.insert(0, str(TEST_DIR))
import test_stage9_science_gate as science_gate  # noqa: E402


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def resolved_policy() -> Dict[str, Any]:
    value = json.loads(
        (Path(__file__).parent / "stage9_science_policy.template.json").read_text(
            encoding="utf-8"
        )
    )
    daughter = value["gate_policy"]["daughter_species"]["C2"]
    daughter.update(
        minimum_ess_fraction=0.20,
        minimum_variance_cost_gain=2.0,
        minimum_gain_bin_fraction=0.90,
    )
    value["gate_policy"]["flux_closure"] = {
        "method": "independent_seed_ci_contains_zero",
        "relative_tolerance": None,
    }
    return value


class Fixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.archive = root / "archive"
        self.policy_path = root / "policy.json"
        self.executable = root / "amps"
        self.input = root / "halley1d.input"
        self.pic_global = root / "picGlobal.dfn"
        self.header_dir = root / "headers"
        self.mpiexec = root / "mpiexec"
        write_json(self.policy_path, resolved_policy())
        self.executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        self.mpiexec.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        os.chmod(self.executable, 0o755)
        os.chmod(self.mpiexec, 0o755)
        self.input.write_text(
            "\n".join(
                (
                    "Stage9ScienceDiagnostics=on",
                    "RadialSamplingLength=10000",
                    "RadialMaxIterations=40000",
                    "RadialConvergenceWindows=3",
                    "CometDiagnostics=on",
                    "ChemistryMode=on",
                    "ChemistryResidualSampling=bernoulli",
                    "ParticleWeightCorrectionWindowMode=on",
                    "ChemistryDaughterPacketMode=on",
                    "ChemistryChannelImportanceMode=on",
                    "ChemistryTotalReactionForcingMode=on",
                    "",
                )
            ),
            encoding="utf-8",
        )
        self.pic_global.write_text("synthetic picGlobal\n", encoding="utf-8")
        self.header_dir.mkdir()
        for name in campaign.HEADER_NAMES:
            (self.header_dir / name).write_text(f"synthetic {name}\n", encoding="utf-8")
        self.state = campaign.initialize_archive(
            self.archive,
            self.policy_path,
            self.executable,
            self.input,
            self.pic_global,
            self.header_dir,
            str(self.mpiexec),
        )
        self.policy = campaign.load_policy(self.archive / "frozen" / "stage9_policy.json")

    def make_attempt(self, arm: str = "enhanced", seed: int = campaign.SEEDS[0]) -> Path:
        attempt, status = campaign.create_attempt(self.archive, self.state, seed, arm)
        status.update(
            state="summarizing",
            return_code=0,
            wall_seconds=120.0,
            fresh_start_verified=True,
            restart_used=False,
        )
        write_json(attempt / "status.json", status)
        (attempt / "error.log").write_text("", encoding="utf-8")
        (attempt / "run.log").write_text(
            "AMPS: radial simulation converged after 40000 iterations and 4 complete sampling windows\n",
            encoding="utf-8",
        )
        (attempt / "comet.radial.mesh.msh").write_text("synthetic mesh\n", encoding="utf-8")
        diagnostics = attempt / "PT" / "diagnostics"
        plots = attempt / "PT" / "plots"
        diagnostics.mkdir(parents=True)
        plots.mkdir(parents=True)
        method_name = "ENHANCED" if arm == "enhanced" else "ANALOG"
        for output in range(4):
            iteration = (output + 1) * 10000
            report = {
                "schema": "amps.comet.stage0-diagnostics.v1",
                "iteration": iteration,
                "output_index": output,
                "sampling_length": 10000,
                "random_seed": seed,
                "random_seed_configured": True,
                "mpi_ranks": 12,
                "emergency_events": {
                    "total": 0,
                    "buffer_utilization": 0,
                    "sampling_population": 0,
                    "allocation_burst": 0,
                },
                "species_flux": [
                    {
                        "species": "H2O",
                        "sampled_surface_source_s-1": 1000.0,
                        "chemical_represented_production_s-1": 0.0,
                        "inner_reemitted_s-1": 0.0,
                        "chemical_destruction_s-1": 0.0,
                        "inner_incident_s-1": 0.0,
                        "outer_escape_s-1": 1000.0,
                        "net_represented_flux_s-1": 0.0,
                    },
                    {
                        "species": "C2",
                        "sampled_surface_source_s-1": 0.0,
                        "chemical_represented_production_s-1": 100.0,
                        "inner_reemitted_s-1": 0.0,
                        "chemical_destruction_s-1": 40.0,
                        "inner_incident_s-1": 0.0,
                        "outer_escape_s-1": 60.0,
                        "net_represented_flux_s-1": 0.0,
                    },
                ],
            }
            write_json(diagnostics / f"comet_stage0.out={output}.json", report)
            (diagnostics / f"comet_stage0.out={output}.csv").write_text(
                "record_type\nmetadata\n", encoding="utf-8"
            )
            checkpoint = {
                "schema": builder.CHECKPOINT_SCHEMA,
                "iteration": iteration,
                "output_index": output,
                "sampling_length": 10000,
                "estimator_mode": method_name,
                "random_seed": seed,
                "mpi_ranks": 12,
                "window_wall_seconds": 30.0,
                "peak_rank_live_particles": 20,
                "total_live_particles": 40,
                "rank_capacity_min": 100,
                "rank_capacity_max": 100,
                "peak_buffer_utilization": 0.40,
                "safe_boundary_audits_min": 1,
                "safe_boundary_audits_max": 1,
                "pending_remainders": 0,
                "pending_reaction_events": 0,
                "species": [
                    {
                        "symbol": species,
                        "model_carriers": 20,
                        "physical_inventory": 20.0,
                        "sum_squared_weights": 20.0,
                        "carrier_ess": 20.0,
                        "occupied_cells": 2,
                        "minimum_positive_cell_ess": 10.0,
                        "correction_min": 1.0,
                        "correction_max": 1.0,
                    }
                    for species in ("H2O", "C2")
                ],
            }
            write_json(diagnostics / f"comet_stage9_checkpoint.out={output}.json", checkpoint)
            (diagnostics / f"comet_stage9_checkpoint.out={output}.csv").write_text(
                "synthetic\n", encoding="utf-8"
            )
            for species, density in (("H2O", 1.0e15), ("C2", 1.0e9)):
                carrier = diagnostics / f"stage9_carriers.{species}.out={output}.csv"
                with carrier.open("w", newline="", encoding="utf-8") as stream:
                    fields = [
                        "cell_index",
                        "model_carriers",
                        "sum_weights",
                        "sum_squared_weights",
                        "carrier_ess",
                        "correction_min",
                        "correction_max",
                        "output_index",
                        "sampling_length",
                        "estimator_mode",
                        "random_seed",
                        "mpi_ranks",
                    ]
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    for cell in range(2):
                        writer.writerow(
                            {
                                "cell_index": cell,
                                "model_carriers": 10,
                                "sum_weights": 10.0,
                                "sum_squared_weights": 10.0,
                                "carrier_ess": 10.0,
                                "correction_min": 1.0,
                                "correction_max": 1.0,
                                "output_index": output,
                                "sampling_length": 10000,
                                "estimator_mode": method_name,
                                "random_seed": seed,
                                "mpi_ranks": 12,
                            }
                        )
                (plots / f"pic.{species}.s=0.out={output}.dat").write_text(
                    'VARIABLES="X", "Number Density[1/m^3]"\n'
                    f"7000 {density}\n17000 {density * 0.5}\n",
                    encoding="utf-8",
                )
        if arm == "enhanced":
            collision = attempt / "collision-preflight"
            collision.mkdir()
            (collision / "error.log").write_text("", encoding="utf-8")
            (collision / "run.log").write_text(
                f"The total number of the MPI processes=12\n"
                f"AMPS: comet random seed={seed} rank-stream=seed+rank residual-phase-stream=keyed\n"
                "PASS transported-daughter-collision ranks=12 reaction=C2_2C "
                "product=1 partner=0 owner=0 daughters=16 partners=64 "
                "daughter-velocity-changes=16 partner-velocity-changes=64 "
                "unchanged-records=80 unfinished-trajectories=0 wrong-owner=0 "
                "cycles=1 trials=64 accepted=64 updates=80 probability-checks=128 "
                "invalid-probabilities=0 roundoff-clamps=0 probability-min=1.0 "
                "probability-max=1.0 correction-error=0.0 represented-rate-error=0.0 "
                "momentum-error=0.0\n",
                encoding="utf-8",
            )
        return attempt


class Stage9ToolingTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = Fixture(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_template_refuses_unresolved_science_decisions(self) -> None:
        template = Path(__file__).parent / "stage9_science_policy.template.json"
        with self.assertRaisesRegex(campaign.PolicyError, "DECISION_REQUIRED"):
            campaign.load_policy(template)

    def test_campaign_flux_policy_accepts_both_methods_and_rejects_mismatch(self) -> None:
        ci_policy = resolved_policy()
        self.assertEqual(
            campaign.validate_policy(ci_policy)["gate_policy"]["flux_closure"]["method"],
            "independent_seed_ci_contains_zero",
        )

        relative_policy = json.loads(json.dumps(ci_policy))
        relative_policy["gate_policy"]["flux_closure"] = {
            "method": "relative_tolerance",
            "relative_tolerance": 0.05,
        }
        self.assertEqual(
            campaign.validate_policy(relative_policy)["gate_policy"]["flux_closure"]["method"],
            "relative_tolerance",
        )

        malformed = json.loads(json.dumps(ci_policy))
        malformed["gate_policy"]["flux_closure"]["relative_tolerance"] = 0.05
        with self.assertRaisesRegex(campaign.PolicyError, "must be null"):
            campaign.validate_policy(malformed)

        malformed = json.loads(json.dumps(relative_policy))
        malformed["gate_policy"]["flux_closure"]["relative_tolerance"] = 0.0
        with self.assertRaisesRegex(campaign.PolicyError, "must be positive"):
            campaign.validate_policy(malformed)

    def test_plan_freezes_complete_provenance_and_balanced_order(self) -> None:
        state = self.fixture.state
        self.assertEqual(state["status"], "planned")
        self.assertEqual(len(state["execution_plan"]), 8)
        self.assertEqual(set(state["frozen"]["headers"]), set(campaign.HEADER_NAMES))
        for path in campaign.frozen_paths(self.fixture.archive)["headers"].values():
            self.assertTrue(path.is_file())
        self.assertEqual(campaign.status_payload(self.fixture.archive)["completed_cases"], 0)

    def test_enhanced_summary_closes_all_raw_evidence(self) -> None:
        attempt = self.fixture.make_attempt("enhanced")
        summary, artifacts = builder.build_summary(self.fixture.archive, attempt)
        self.assertEqual(summary["method"], "FORCED")
        self.assertEqual(summary["convergence"]["complete_steady_windows"], 3)
        self.assertEqual(summary["collision_audit"]["transported_forced_daughter_events"], 64)
        self.assertEqual(summary["daughters"]["C2"]["model_carriers"], 60.0)
        self.assertEqual(summary["daughters"]["C2"]["zero_cell_fraction"], 0.0)
        self.assertTrue(any(path.name.startswith("stage9_carriers.C2") for path in artifacts))

    def test_analog_summary_has_exact_zero_forced_collision_count(self) -> None:
        attempt = self.fixture.make_attempt("analog")
        summary, _artifacts = builder.build_summary(self.fixture.archive, attempt)
        self.assertEqual(summary["method"], "ANALOG")
        self.assertEqual(summary["collision_audit"]["transported_forced_daughter_events"], 0)

    def test_missing_checkpoint_or_collision_log_fails_closed(self) -> None:
        attempt = self.fixture.make_attempt("enhanced")
        (attempt / "PT" / "diagnostics" / "comet_stage9_checkpoint.out=3.json").unlink()
        with self.assertRaisesRegex(builder.SummaryError, "index sets differ"):
            builder.build_summary(self.fixture.archive, attempt)

        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = Fixture(Path(self.temporary.name))
        attempt = self.fixture.make_attempt("enhanced")
        (attempt / "collision-preflight" / "run.log").unlink()
        with self.assertRaisesRegex(builder.SummaryError, "collision-preflight"):
            builder.build_summary(self.fixture.archive, attempt)

    def test_carrier_checkpoint_disagreement_fails_closed(self) -> None:
        attempt = self.fixture.make_attempt("enhanced")
        path = attempt / "PT" / "diagnostics" / "stage9_carriers.C2.out=3.csv"
        text = path.read_text(encoding="utf-8").replace(",10.0,10.0,10.0,", ",11.0,10.0,10.0,", 1)
        path.write_text(text, encoding="utf-8")
        with self.assertRaises(builder.SummaryError):
            builder.build_summary(self.fixture.archive, attempt)

    def test_reconcile_marks_dead_running_attempt_interrupted(self) -> None:
        attempt, status = campaign.create_attempt(
            self.fixture.archive, self.fixture.state, campaign.SEEDS[0], "analog"
        )
        status.update(state="running", child_pid=99999999, child_start_ticks=1)
        write_json(attempt / "status.json", status)
        campaign.reconcile(self.fixture.archive, self.fixture.policy, self.fixture.state)
        updated = json.loads((attempt / "status.json").read_text(encoding="utf-8"))
        self.assertEqual(updated["state"], "interrupted")

    def test_published_campaign_is_valid_gate_evidence(self) -> None:
        for seed in campaign.SEEDS:
            for arm in campaign.ARMS:
                attempt = self.fixture.make_attempt(arm, seed)
                summary, _artifacts = builder.build_summary(self.fixture.archive, attempt)
                path = attempt / "stage9_summary.json"
                write_json(path, summary)
                self.fixture.state["completed"].append(
                    {
                        "method": "ANALOG" if arm == "analog" else "FORCED",
                        "seed": seed,
                        "summary": str(path.relative_to(self.fixture.archive)),
                        "sha256": campaign.sha256_file(path),
                    }
                )
        campaign.publish_manifest(
            self.fixture.archive, self.fixture.policy, self.fixture.state
        )
        result = science_gate.evaluate(self.fixture.archive / "RUN_MANIFEST.json")
        # Identical synthetic seed profiles intentionally demonstrate no
        # variance improvement, so the science decision is FAIL rather than
        # PASS.  The important integration assertion is that the manifest and
        # all eight production summaries are valid, complete evidence.
        self.assertEqual(result["outcome"], "FAIL_SCIENCE")
        self.assertEqual(result["run_count"], 8)


if __name__ == "__main__":
    unittest.main(verbosity=2)

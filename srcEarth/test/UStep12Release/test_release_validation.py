#!/usr/bin/env python3
"""Numerical/contract tests for the Step-12 fail-closed release gate.

The fixtures are deliberately synthetic: they test the evaluator, not the physical
validity of AMPS.  A synthetic PASS never replaces a U/I/C/F/O production result.  The
tests inject mismatches, relaxed tolerances, stale hashes, missing observation output,
retuning, and restart mutations to prove those conditions remain release failures.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, Mapping, MutableMapping, Tuple


HERE = Path(__file__).resolve().parent
SRC_EARTH = HERE.parents[1]
sys.path.insert(0, str(SRC_EARTH))

from release_validation.compare_cross_path import compare_campaigns  # noqa: E402
from release_validation.contract import (  # noqa: E402
    BUILD_SCHEMA,
    CAPABILITY_SCHEMA,
    COMMON_PHYSICS_TAG,
    EVIDENCE_SCHEMA,
    NORMALIZATION_POLICY,
    PHASE_1_SCOPE,
    REQUIRED_GATE_IDS,
    REQUIRED_PARITY_ROLES,
    REQUIRED_RESOURCE_CASES,
    RESOURCE_SCHEMA,
    SCHEMA,
    ReleaseValidationError,
    ResourceResolver,
    sha256_file,
)
from release_validation.run_release import evaluate_release  # noqa: E402


REVISION = "0123456789abcdef0123456789abcdef01234567"
SOURCE_DIGEST = "1" * 64
EXECUTABLE_DIGEST = "2" * 64
INPUT_DIGEST = "3" * 64
BOUNDARY_DIGEST = "4" * 64
RESPONSE_DIGEST = "5" * 64
PREREGISTRATION_DIGEST = "6" * 64


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def reference(path: Path, base: Path) -> Dict[str, str]:
    return {"path": str(path.relative_to(base)), "sha256": sha256_file(path)}


def write_parity_csv(path: Path, outcome: str = "ALLOWED", value: float = 1.0) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["id", "outcome", "value"])
        writer.writerow(["sample-0001", outcome, "%.17g" % value])


class Fixture:
    """Construct a complete release bundle with real SHA-256 relationships."""

    def __init__(self, root: Path):
        self.root = root
        self.manifest_path = root / "release.json"
        self.output_path = root / "summary.json"
        self.artifact_path = root / "evidence" / "result.dat"
        self.artifact_path.parent.mkdir(parents=True)
        self.artifact_path.write_text("closed numerical artifact\n", encoding="utf-8")
        self.manifest = self._build_manifest()
        self.save_manifest()

    def _build_file(self, kind: str) -> Path:
        path = self.root / (kind.lower() + "_build.json")
        write_json(path, {
            "schema": BUILD_SCHEMA,
            "build_kind": kind,
            "common_physics_tag": COMMON_PHYSICS_TAG,
            "common_physics_source_sha256": SOURCE_DIGEST,
            "source_revision": REVISION,
            "dirty": False,
            "compiler": "g++ 13.2.0 -std=c++17",
            "dependencies": {"AMPS": REVISION, "MPI": "OpenMPI-4.1"},
            "executable_sha256": EXECUTABLE_DIGEST,
        })
        return path

    def _evidence_file(self, gate_id: str) -> Path:
        path = self.root / "evidence" / (gate_id.replace("-", "_") + ".json")
        observational = gate_id.startswith("O") or gate_id in {"F8", "F9", "F10", "F17"}
        value: Dict[str, Any] = {
            "schema": EVIDENCE_SCHEMA,
            "gate_id": gate_id,
            "status": "PASS",
            "result_line": "RESULT: PASS",
            "exit_code": 0,
            "dry_run": False,
            "command": "python3 run_%s.py --production" % gate_id.lower(),
            "source_revision": REVISION,
            "dirty": False,
            "compiler": "g++ 13.2.0 -std=c++17",
            "dependencies": {"AMPS": REVISION, "MPI": "OpenMPI-4.1"},
            "common_physics_tag": COMMON_PHYSICS_TAG,
            "executable_sha256": EXECUTABLE_DIGEST,
            "input_sha256": INPUT_DIGEST,
            "field_id": "DIPOLE_REFERENCE",
            "snapshot_id": "snapshot-001",
            "boundary_sha256": BOUNDARY_DIGEST,
            "response_sha256": RESPONSE_DIGEST,
            "species": "PROTON",
            "epoch_utc": "2012-05-17T01:56:00Z",
            "frame": "GSM",
            "mover": "BORIS",
            "integrator_tolerances": {"gyro_angle_rad": 0.15, "boundary_m": 10.0},
            "grid_history": [{"level": 0, "energy_nodes": 80, "angular_nodes": 144}],
            "mpi_ranks": 4,
            "threads": 16,
            "scheduler": "STATIC",
            "random_seeds": [1729],
            "wall_time_s": 1.25,
            "termination_counts": {"sampled": 144, "allowed": 120, "forbidden": 24, "unresolved": 0},
            "artifacts": [reference(self.artifact_path, self.root)],
        }
        if observational:
            value["observation_comparison_count"] = 8
            value["normalization_policy"] = NORMALIZATION_POLICY
        if gate_id == "O4":
            value.update({"frozen": True, "retuned": False, "reported_outcome": "PASS"})
        write_json(path, value)
        return path

    def _capabilities_file(self) -> Path:
        path = self.root / "capabilities.json"
        supported = (
            "STANDALONE_CUTOFF", "STANDALONE_FLUX_SPECTRUM", "SWMF_CUTOFF",
            "SWMF_FLUX_SPECTRUM", "OFFLINE_SWMF_REPLAY",
        )
        unsupported = (
            "DYNAMIC_EB_CHARACTERISTICS", "LONG_DURATION_TRAPPING",
            "LOCAL_ACCELERATION", "LOSS_PHYSICS", "FORWARD_TRAPPED_FLUX",
        )
        write_json(path, {"schema": CAPABILITY_SCHEMA, "capabilities": [
            {"id": name, "status": "SUPPORTED", "basis": "validated Phase-1 path"}
            for name in supported
        ] + [
            {"id": name, "status": "UNSUPPORTED", "basis": "outside frozen-snapshot Phase-1 scope"}
            for name in unsupported
        ]})
        return path

    def _resources_file(self) -> Path:
        path = self.root / "resources.json"
        write_json(path, {"schema": RESOURCE_SCHEMA, "estimates": [
            {"case": case, "mpi_ranks": 4, "threads": 16, "memory_gib": 32,
             "wall_minutes": 30, "disk_gib": 4,
             "basis": "measured on the documented validation host"}
            for case in sorted(REQUIRED_RESOURCE_CASES)
        ]})
        return path

    def _campaign(self, kind: str) -> Dict[str, Any]:
        pairs = []
        for role in REQUIRED_PARITY_ROLES:
            stem = "%s_%s" % (kind.lower(), role.lower())
            left = self.root / "parity" / (stem + "_reference.csv")
            right = self.root / "parity" / (stem + "_candidate.csv")
            left.parent.mkdir(parents=True, exist_ok=True)
            write_parity_csv(left)
            write_parity_csv(right)
            pairs.append({
                "role": role,
                "reference": reference(left, self.root),
                "candidate": reference(right, self.root),
                "key_columns": ["id"],
                "exact_columns": ["outcome"],
                "numeric_columns": ["value"],
                "tolerance_class": "MESH" if kind == "SAMPLED_ANALYTIC_MESH" else "KERNEL",
                "rtol": 0.05 if kind == "SAMPLED_ANALYTIC_MESH" else 1.0e-10,
                "atol": 1.0e-12,
            })
        return {
            "campaign_id": "fixture-" + kind.lower(),
            "kind": kind,
            "common_physics_tag": COMMON_PHYSICS_TAG,
            "phase_1_scope": PHASE_1_SCOPE,
            "reference_field_source": (
                "ANALYTIC_OR_PHENOMENOLOGICAL" if kind == "SAMPLED_ANALYTIC_MESH" else "LIVE_SWMF"
            ),
            "candidate_field_source": (
                "SAMPLED_MESH" if kind == "SAMPLED_ANALYTIC_MESH" else "OFFLINE_SNAPSHOT_REPLAY"
            ),
            "invariants": {
                "field_state_id": "field-state-001", "species": "PROTON",
                "epoch_utc": "2012-05-17T01:56:00Z", "frame": "GSM",
                "boundary_sha256": BOUNDARY_DIGEST,
                "integration_sha256": "7" * 64,
                "output_schema_sha256": "8" * 64,
            },
            "pairs": pairs,
        }

    def _build_manifest(self) -> MutableMapping[str, Any]:
        standalone = self._build_file("STANDALONE")
        coupled = self._build_file("SWMF_COUPLED")
        capabilities = self._capabilities_file()
        resources = self._resources_file()
        gate_entries = []
        for gate_id in REQUIRED_GATE_IDS:
            evidence = self._evidence_file(gate_id)
            gate_entries.append({"gate_id": gate_id, "evidence": reference(evidence, self.root)})
        return {
            "schema": SCHEMA,
            "release_id": "synthetic-step12-contract-test",
            "phase_1_scope": PHASE_1_SCOPE,
            "common_physics": {
                "tag": COMMON_PHYSICS_TAG,
                "source_sha256": SOURCE_DIGEST,
                "standalone_build": reference(standalone, self.root),
                "coupled_build": reference(coupled, self.root),
            },
            "gate_evidence": gate_entries,
            "cross_path_campaigns": [
                self._campaign("SAMPLED_ANALYTIC_MESH"),
                self._campaign("SWMF_OFFLINE_REPLAY"),
            ],
            "holdout": {
                "gate_id": "O4", "frozen": True, "retuned": False,
                "selected_before_release": True, "outcome_reported": True,
                "normalization_policy": NORMALIZATION_POLICY,
                "preregistration_sha256": PREREGISTRATION_DIGEST,
            },
            "capabilities": reference(capabilities, self.root),
            "resource_estimates": reference(resources, self.root),
            "hooks": {
                "ccmc_archive": "python3 tools/archive_ccmc.py --manifest run.json",
                "swmf_snapshot_export": "./SWMF.exe -export-sep-snapshot snapshot.bin",
                "swmf_offline_replay": "./amps -mode 3d -i replay.in",
            },
            "commands": {
                "standalone_cutoff_only": "mpirun -np 4 ./amps -mode gridless -i cutoff.in",
                "standalone_combined": "mpirun -np 4 ./amps -mode 3d -i products.in",
                "swmf_cutoff_only": "./SWMF.exe -i coupled_cutoff.in",
                "swmf_combined": "./SWMF.exe -i coupled_products.in",
                "offline_replay": "mpirun -np 4 ./amps -mode 3d -i replay.in",
            },
        }

    def save_manifest(self) -> None:
        write_json(self.manifest_path, self.manifest)

    def gate_evidence(self, gate_id: str) -> Tuple[Path, MutableMapping[str, Any]]:
        entry = next(item for item in self.manifest["gate_evidence"] if item["gate_id"] == gate_id)
        path = self.root / entry["evidence"]["path"]
        return path, json.loads(path.read_text(encoding="utf-8"))

    def replace_gate_evidence(self, gate_id: str, value: Mapping[str, Any]) -> None:
        path, _ = self.gate_evidence(gate_id)
        write_json(path, value)
        entry = next(item for item in self.manifest["gate_evidence"] if item["gate_id"] == gate_id)
        entry["evidence"]["sha256"] = sha256_file(path)
        self.save_manifest()


class Step12ReleaseTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="step12-release-test-")
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_complete_matrix_and_both_parity_campaigns_pass(self) -> None:
        report, reused = evaluate_release(self.fixture.manifest_path, self.fixture.output_path)
        self.assertFalse(reused)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(len(report["prepared"]["gates"]), len(REQUIRED_GATE_IDS))
        self.assertEqual(
            {item["kind"] for item in report["cross_path"]["campaigns"]},
            {"SAMPLED_ANALYTIC_MESH", "SWMF_OFFLINE_REPLAY"},
        )
        for campaign in report["cross_path"]["campaigns"]:
            self.assertEqual({item["role"] for item in campaign["pairs"]}, set(REQUIRED_PARITY_ROLES))

    def test_numeric_and_exact_trajectory_mismatches_fail(self) -> None:
        sampled = self.fixture.manifest["cross_path_campaigns"][0]
        field_pair = next(pair for pair in sampled["pairs"] if pair["role"] == "FIELD_SAMPLES")
        field_path = self.root / field_pair["candidate"]["path"]
        write_parity_csv(field_path, value=1.2)
        field_pair["candidate"]["sha256"] = sha256_file(field_path)
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "cross-path comparisons failed"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

        # Rebuild the clean fixture, then prove a discrete terminal outcome is exact
        # even where the numeric mesh envelope would otherwise permit a difference.
        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory(prefix="step12-release-test-")
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)
        sampled = self.fixture.manifest["cross_path_campaigns"][0]
        trajectory_pair = next(pair for pair in sampled["pairs"] if pair["role"] == "TRAJECTORIES")
        trajectory_path = self.root / trajectory_pair["candidate"]["path"]
        write_parity_csv(trajectory_path, outcome="FORBIDDEN")
        trajectory_pair["candidate"]["sha256"] = sha256_file(trajectory_path)
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "cross-path comparisons failed"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

    def test_relaxed_mesh_or_replay_tolerance_is_rejected(self) -> None:
        sampled_pair = self.fixture.manifest["cross_path_campaigns"][0]["pairs"][0]
        sampled_pair["rtol"] = 0.0500001
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "relaxes the MESH ceiling"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory(prefix="step12-release-test-")
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)
        replay_pair = self.fixture.manifest["cross_path_campaigns"][1]["pairs"][0]
        replay_pair["tolerance_class"] = "MESH"
        replay_pair["rtol"] = 1.0e-9
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "replay comparison exceeds kernel precision"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

    def test_gate_registry_failed_dry_run_and_missing_comparison_fail(self) -> None:
        self.fixture.manifest["gate_evidence"].pop()
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "gate registry differs"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory(prefix="step12-release-test-")
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)
        _, evidence = self.fixture.gate_evidence("U-F01")
        evidence["dry_run"] = True
        self.fixture.replace_gate_evidence("U-F01", evidence)
        with self.assertRaisesRegex(ReleaseValidationError, "dry run"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory(prefix="step12-release-test-")
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)
        _, evidence = self.fixture.gate_evidence("F8")
        evidence["observation_comparison_count"] = 0
        self.fixture.replace_gate_evidence("F8", evidence)
        with self.assertRaisesRegex(ReleaseValidationError, "produced no observation comparison"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

    def test_common_build_identity_and_o4_freeze_are_mandatory(self) -> None:
        coupled_ref = self.fixture.manifest["common_physics"]["coupled_build"]
        coupled_path = self.root / coupled_ref["path"]
        coupled = json.loads(coupled_path.read_text(encoding="utf-8"))
        coupled["common_physics_tag"] = "incompatible-physics"
        write_json(coupled_path, coupled)
        coupled_ref["sha256"] = sha256_file(coupled_path)
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "different common physics tag"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory(prefix="step12-release-test-")
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)
        self.fixture.manifest["holdout"]["retuned"] = True
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "must not be retuned"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

    def test_restart_reuses_only_immutable_verified_pass(self) -> None:
        first, reused = evaluate_release(self.fixture.manifest_path, self.fixture.output_path)
        self.assertEqual(first["status"], "PASS")
        self.assertFalse(reused)
        second, reused = evaluate_release(
            self.fixture.manifest_path, self.fixture.output_path, restart=True
        )
        self.assertEqual(second["status"], "PASS")
        self.assertTrue(reused)

        # Modify a referenced result without changing its declared digest.  Restart
        # must hash inputs first and fail rather than returning the cached PASS.
        self.fixture.artifact_path.write_text("mutated after PASS\n", encoding="utf-8")
        with self.assertRaisesRegex(ReleaseValidationError, "SHA-256 mismatch"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path, restart=True)

    def test_placeholders_and_missing_capabilities_fail_fast(self) -> None:
        self.fixture.manifest["commands"]["swmf_combined"] = "<insert command>"
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "placeholder"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory(prefix="step12-release-test-")
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)
        capabilities_ref = self.fixture.manifest["capabilities"]
        capabilities_path = self.root / capabilities_ref["path"]
        capabilities = json.loads(capabilities_path.read_text(encoding="utf-8"))
        capabilities["capabilities"] = [
            row for row in capabilities["capabilities"]
            if row["id"] != "SWMF_FLUX_SPECTRUM"
        ]
        write_json(capabilities_path, capabilities)
        capabilities_ref["sha256"] = sha256_file(capabilities_path)
        self.fixture.save_manifest()
        with self.assertRaisesRegex(ReleaseValidationError, "required supported capability"):
            evaluate_release(self.fixture.manifest_path, self.fixture.output_path)

    def test_cli_result_lines_and_validate_only_semantics(self) -> None:
        command = [sys.executable, str(SRC_EARTH / "release_validation" / "run_release.py"),
                   "--manifest", str(self.fixture.manifest_path),
                   "--output", str(self.fixture.output_path)]
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("RESULT: PASS", completed.stdout)

        validate_output = self.root / "validate_only.json"
        completed = subprocess.run(command[:-1] + [str(validate_output), "--validate-only"],
                                   text=True, capture_output=True, check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("RESULT: VALID_NOT_RELEASE", completed.stdout)
        self.assertNotIn("RESULT: PASS", completed.stdout)

        # The dedicated comparator is independently runnable and follows the same
        # explicit RESULT/non-zero runner contract.
        compare_output = self.root / "comparison.json"
        compare_command = [sys.executable,
                           str(SRC_EARTH / "release_validation" / "compare_cross_path.py"),
                           "--manifest", str(self.fixture.manifest_path),
                           "--output", str(compare_output)]
        completed = subprocess.run(compare_command, text=True, capture_output=True, check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("RESULT: PASS", completed.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)

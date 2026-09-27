#!/usr/bin/env python3
"""Audit Step-12 production wiring and prove the old test list is unchanged.

This is intentionally a source contract, not a substring-only physics test.  Numerical
failure injection lives in test_release_validation.py; here we ensure all production
surfaces publish the common identity and that registering Step 12 did not edit any
pre-existing command, expected state, or last-pass hash.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SRC_EARTH = HERE.parents[1]
PRE_STEP12_LIST_SHA256 = "24997e40773679968370ff3c927fe882327e6e160567e254597ef17f5a85a685"
STEP12_LIST_BLOCK = "P srcEarth/test/UStep12Release/run_test.sh\nlast pass:\n"


def text(relative: str) -> str:
    return (SRC_EARTH / relative).read_text(encoding="utf-8")


class Step12SourceContract(unittest.TestCase):
    def test_common_identity_reaches_all_manifests_and_numeric_products(self) -> None:
        common = text("util/CommonPhysicsRelease.h")
        self.assertIn("sep-in-geospace-phase1-static-characteristics-v1", common)
        self.assertIn("INSTANTANEOUS_QUASI_STATIC_MAGNETIC", common)
        for relative in (
            "util/StandaloneProductContract.h",
            "util/SWMFCoupledAccessContract.h",
            "util/SWMFCoupledProductsContract.h",
        ):
            source = text(relative)
            self.assertIn("common_physics_tag", source, relative)
            self.assertIn("CommonPhysicsRelease::kTag", source, relative)
            self.assertIn("common_physics_scope", source, relative)
        for relative in ("3d/CutoffRigidityMode3D.cpp", "3d/DensityMode3D.cpp"):
            source = text(relative)
            self.assertIn("AUXDATA COMMON_PHYSICS_TAG", source, relative)
            self.assertIn("AUXDATA COMMON_PHYSICS_SCOPE", source, relative)

    def test_release_tooling_keeps_fixed_matrix_roles_and_ceiling_values(self) -> None:
        source = text("release_validation/contract.py")
        comparator = text("release_validation/compare_cross_path.py")
        runner = text("release_validation/run_release.py")
        for token in ("U-F%02d", "I-F%02d", "C%d", "F%d", '"O1", "O2", "O3", "O4"'):
            self.assertIn(token, source)
        for role in (
            "FIELD_SAMPLES", "TRAJECTORIES", "DIRECTIONAL_ACCESS", "CUTOFF",
            "SPECTRUM", "DENSITY", "DETECTOR",
        ):
            self.assertIn(role, source)
        self.assertIn('"KERNEL": (1.0e-10, 1.0e-12)', source)
        self.assertIn('"INTEGRATED": (0.02, 1.0e-12)', source)
        self.assertIn('"MESH": (0.05, 1.0e-12)', source)
        self.assertIn("must classify every CSV column exactly once", comparator)
        self.assertIn("NaN/Inf cannot pass parity", comparator)
        self.assertIn("RESTART: SKIPPED_VERIFIED_PASS", runner)
        self.assertIn("RESULT: VALID_NOT_RELEASE", runner)

    def test_help_and_readmes_document_commands_scope_and_no_retuning(self) -> None:
        help_source = text("util/cutoff_cli.cpp")
        self.assertIn("Phase-1 release and reproducible command forms", help_source)
        self.assertIn("standalone_cutoff.in", help_source)
        self.assertIn("standalone_products.in", help_source)
        self.assertIn("run_release.py --manifest", help_source)
        self.assertIn("long-duration trapping", help_source)
        required_docs = (
            "README.md", "READ.md", "util/README.md", "3d/README.md",
            "3d_forward_swmf/README.md", "examples/README.md",
            "examples/step12_release/README.md", "release_validation/README.md",
            "test/README.md", "test/UStep12Release/README.md",
        )
        for relative in required_docs:
            body = text(relative)
            self.assertTrue("Step 12" in body or "Step-12" in body, relative)
        release_readme = text("release_validation/README.md")
        for phrase in ("no-platform-scale", "retuned=false", "--restart",
                       "cutoff only", "cutoff + flux/spectrum"):
            self.assertIn(phrase, release_readme)

    def test_manifest_generator_has_every_gate_and_parity_pair(self) -> None:
        with tempfile.TemporaryDirectory(prefix="step12-skeleton-") as directory:
            output = Path(directory) / "release.json"
            completed = subprocess.run(
                [sys.executable,
                 str(SRC_EARTH / "release_validation" / "make_manifest_skeleton.py"),
                 "--output", str(output)],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("RESULT: SKELETON_WRITTEN", completed.stdout)
            manifest = json.loads(output.read_text(encoding="utf-8"))
            # 12 U + 8 I + 19 C + 17 F + 4 O = 60 fixed Phase-1 gates.
            self.assertEqual(len(manifest["gate_evidence"]), 60)
            self.assertEqual(len(manifest["cross_path_campaigns"]), 2)
            self.assertTrue(all(len(campaign["pairs"]) == 7
                                for campaign in manifest["cross_path_campaigns"]))

    def test_test_list_addition_is_the_only_step12_list_change(self) -> None:
        current = text("test/list")
        self.assertEqual(current.count(STEP12_LIST_BLOCK), 1)
        self.assertIn("Step 2 through Step 12", current)
        projected = current.replace(STEP12_LIST_BLOCK, "", 1)
        projected = projected.replace("Step 2 through Step 12", "Step 2 through Step 11", 1)
        digest = hashlib.sha256(projected.encode("utf-8")).hexdigest()
        self.assertEqual(digest, PRE_STEP12_LIST_SHA256)


if __name__ == "__main__":
    unittest.main(verbosity=2)

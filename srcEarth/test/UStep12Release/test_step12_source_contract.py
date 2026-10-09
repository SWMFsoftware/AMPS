#!/usr/bin/env python3
"""Audit Step-12 production wiring and the approved test-list structure.

This is intentionally a source contract, not a substring-only physics test.
Numerical failure injection lives in test_release_validation.py; here we ensure
all production surfaces publish the common identity and that every reviewed test
command, expected state, ordering rule, and comment remains frozen.

``last pass:`` values are different: they are runner-owned mutable provenance.
The production runner legitimately updates those hashes after a pass, so folding
their current values into an immutable source digest makes a successful run break
the next run.  This contract validates their syntax and placement, then normalizes
only the optional hash before calculating the structural digest.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SRC_EARTH = HERE.parents[1]
APPROVED_NORMALIZED_LIST_SHA256 = "e68727bf79510f78feac73c22e562e72354595b6bdc9a46d7d4b8e5bc52dcc13"
STEP12_COMMAND = "P srcEarth/test/UStep12Release/run_test.sh"
COMMIT_RE = r"[0-9a-f]{40}"
# Looped list entries are rewritten by test_runner.py as a brace-delimited
# value per expansion. Empty items preserve a variant with no passing commit.
# Accept that production syntax while still rejecting every partial SHA.
LOOP_COMMITS_RE = r"\{(?:%s)?(?:,(?:%s)?)*\}" % (COMMIT_RE, COMMIT_RE)
LAST_PASS_RE = re.compile(
    r"^(!?last pass:)(?:\s+(%s|%s))?\s*$" % (COMMIT_RE, LOOP_COMMITS_RE)
)


def text(relative: str) -> str:
    return (SRC_EARTH / relative).read_text(encoding="utf-8")


def normalized_test_list(body: str) -> str:
    """Remove only runner-written commit values from a validated list.

    The active/commented marker remains part of the digest.  Consequently this
    does not hide enabling, disabling, adding, removing, reordering, or editing
    a test; it makes only a valid 40-hex provenance value non-structural.
    """

    normalized = []
    for line in body.splitlines(keepends=True):
        content = line.rstrip("\r\n")
        ending = line[len(content):]
        match = LAST_PASS_RE.fullmatch(content)
        if match:
            normalized.append(match.group(1) + ending)
        else:
            normalized.append(line)
    return "".join(normalized)


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

    def test_test_list_provenance_and_approved_structure(self) -> None:
        current = text("test/list")

        # Every provenance line is either empty or a complete Git SHA-1. This
        # fails explicitly for truncated/malformed values instead of relying on
        # a less-informative whole-file digest mismatch.
        lines = current.splitlines()
        provenance_lines = [
            (index, line) for index, line in enumerate(lines)
            if line.startswith("last pass:") or line.startswith("!last pass:")
        ]
        self.assertTrue(provenance_lines)
        for index, line in provenance_lines:
            self.assertRegex(line, LAST_PASS_RE, "test/list line %d" % (index + 1))

        # Prove the normalization has exactly the intended mutability: changing
        # a valid scalar SHA or using the runner's loop-value syntax cannot
        # perturb normalized structure, while a truncated SHA is never valid.
        advanced = re.sub(
            r"(?m)^last pass:\s+[0-9a-f]{40}$",
            "last pass: " + "0" * 40,
            current,
            count=1,
        )
        self.assertEqual(normalized_test_list(current), normalized_test_list(advanced))
        loop_provenance = "last pass: {%s,,%s}\n" % ("1" * 40, "2" * 40)
        self.assertEqual(normalized_test_list(loop_provenance), "last pass:\n")
        self.assertIsNone(LAST_PASS_RE.fullmatch("last pass: " + "a" * 39))

        # Step 12 stays an independent active gate and owns exactly one adjacent
        # active provenance record. A hash may be present after a successful run.
        step12_indexes = [i for i, line in enumerate(lines) if line == STEP12_COMMAND]
        self.assertEqual(len(step12_indexes), 1)
        step12_index = step12_indexes[0]
        self.assertLess(step12_index + 1, len(lines))
        self.assertRegex(lines[step12_index + 1], r"^last pass:(?:\s+[0-9a-f]{40})?$")
        self.assertIn("Step 2 through Step 12", current)

        # Active tests must be checkout-portable. In particular, C9/C10 must
        # launch the Earth binary built by this checkout and write beneath this
        # checkout, rather than a developer's historical home directory.
        active_commands = [
            line for line in lines if line.startswith("P ") or line.startswith("F ")
        ]
        for command in active_commands:
            self.assertNotIn("/home/", command)
            self.assertNotIn("~/", command)
        c9_commands = [line for line in active_commands if "test/C9/run_C9.py" in line]
        c10_commands = [line for line in active_commands if "test/C10/run_C10.py" in line]
        self.assertTrue(c9_commands)
        self.assertTrue(c10_commands)
        self.assertTrue(all("--amps ./amps" in line for line in c9_commands + c10_commands))
        self.assertTrue(all("--output-root test_output/" in line for line in c10_commands))

        digest = hashlib.sha256(
            normalized_test_list(current).encode("utf-8")
        ).hexdigest()
        self.assertEqual(digest, APPROVED_NORMALIZED_LIST_SHA256)


if __name__ == "__main__":
    unittest.main(verbosity=2)

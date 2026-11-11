#!/usr/bin/env python3
"""Dependency-free mutation tests for the all-stage JSON receipt writer."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence


THIS_DIR = Path(__file__).resolve().parent
AMPS_ROOT = THIS_DIR.parents[1]
WRITER = THIS_DIR / "write_all_stage_result.py"
OBJECT_A = "1" * 40
OBJECT_B = "2" * 40


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def base_command(root: Path, log_root: Path, output: Path) -> list[str]:
    command = [
        sys.executable,
        str(WRITER),
        "--repository-root",
        str(root),
        "--head-commit",
        OBJECT_A,
        "--head-tree",
        OBJECT_B,
        "--branch",
        'branch/with-"quote"',
        "--log-root",
        str(log_root),
        "--summary",
        str(log_root / "summary.txt"),
        "--output",
        str(output),
        "--expected-overall",
        "PASS",
    ]
    for stage in range(4, 10):
        command.extend(("--mode", f"stage{stage}", "archived"))
    return command


def run(command: Sequence[str], expected: int = 0) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != expected:
        raise AssertionError(
            f"command returned {completed.returncode}, expected {expected}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def main() -> int:
    cases = 0
    # Receipt fixtures follow the same AMPS-local scratch policy as the
    # synthetic Git-history test. An explicit parent avoids any dependency on
    # /tmp or platform fallback; the context manager removes only its unique
    # child directory, including when an assertion raises an exception.
    with tempfile.TemporaryDirectory(
            prefix="stage10-all-stage-result-", dir=str(AMPS_ROOT)) as temporary:
        root = Path(temporary)

        # Quotation, backslashes, whitespace, Unicode, and a newline in a log
        # name all reach json.dumps through argv without shell interpolation.
        log_root = root / "logs focused"
        log_root.mkdir()
        summary = log_root / "summary.txt"
        summary.write_text("focused summary\n", encoding="utf-8")
        odd_log = log_root / 'test "quoted" \\ snowman-☃\nline.log'
        odd_log.write_text("PASS with arbitrary JSON characters\n", encoding="utf-8")
        output = log_root / "all_stage_result.json"
        command = base_command(root, log_root, output)
        command.extend(("--test", 'test."quoted"\\id', "stage0", "PASS", "0", str(odd_log)))
        run(command)
        value = json.loads(output.read_text(encoding="utf-8"))
        assert value["schema"] == "amps.comet.all-stage-test-result.v1"
        assert value["overall_result"] == "PASS"
        assert value["overall_scope"] == "focused_implementation_only"
        assert value["repository"]["branch"] == 'branch/with-"quote"'
        assert value["tests"][0]["id"] == 'test."quoted"\\id'
        assert value["tests"][0]["log_sha256"] == digest(odd_log)
        assert value["summary"]["sha256"] == digest(summary)
        assert value["counts"]["tests"] == {
            "run": 1,
            "passed": 1,
            "skipped": 0,
            "failed": 0,
        }
        assert value["interpretation"]["stage9_science"]["result"] == "NOT_EVALUATED"
        cases += 1

        # A supplied manifest has a distinct test record and interpretation;
        # its PASS does not alter the focused-only counts.
        science_root = root / "science"
        science_root.mkdir()
        (science_root / "summary.txt").write_text("science summary\n", encoding="utf-8")
        focused_log = science_root / "focused.log"
        focused_log.write_text("PASS focused\n", encoding="utf-8")
        science_log = science_root / "science.log"
        science_log.write_text("PASS science\n", encoding="utf-8")
        manifest = root / "RUN_MANIFEST.json"
        manifest.write_text('{"schema":"fixture"}\n', encoding="utf-8")
        science_output = science_root / "all_stage_result.json"
        command = base_command(root, science_root, science_output)
        command.extend(("--science-manifest", str(manifest)))
        command.extend(("--test", "focused", "stage9", "PASS", "0", str(focused_log)))
        command.extend(
            (
                "--test",
                "stage9.science-acceptance",
                "stage9",
                "PASS",
                "0",
                str(science_log),
            )
        )
        run(command)
        value = json.loads(science_output.read_text(encoding="utf-8"))
        assert value["overall_scope"] == "focused_implementation_and_stage9_science"
        assert value["interpretation"]["focused_implementation"]["counts"]["run"] == 1
        assert value["interpretation"]["stage9_science"]["accepted"] is True
        assert value["interpretation"]["stage9_science"]["manifest_sha256"] == digest(manifest)
        cases += 1

        # Mutation cases prove fail-closed consistency, containment, and
        # immutability rather than only exercising the successful serializer.
        mutations = []
        duplicate = base_command(root, log_root, root / "duplicate.json")
        duplicate.extend(("--test", "same", "stage0", "PASS", "0", str(odd_log)))
        duplicate.extend(("--test", "same", "stage1", "PASS", "0", str(odd_log)))
        mutations.append(duplicate)

        outside = root / "outside.log"
        outside.write_text("outside\n", encoding="utf-8")
        escaped = base_command(root, log_root, root / "escaped.json")
        escaped.extend(("--test", "escaped", "stage0", "PASS", "0", str(outside)))
        mutations.append(escaped)

        mismatch = base_command(root, log_root, root / "mismatch.json")
        mismatch.extend(("--test", "failed", "stage0", "FAIL", "7", str(odd_log)))
        mutations.append(mismatch)

        missing_science_test = base_command(root, log_root, root / "missing-science.json")
        missing_science_test.extend(("--science-manifest", str(manifest)))
        missing_science_test.extend(("--test", "focused", "stage9", "PASS", "0", str(odd_log)))
        mutations.append(missing_science_test)

        overwrite = base_command(root, log_root, output)
        overwrite.extend(("--test", "new", "stage0", "PASS", "0", str(odd_log)))
        mutations.append(overwrite)

        for mutation in mutations:
            run(mutation, expected=2)
            cases += 1

        # A skipped required test makes the overall/focused result fail while
        # retaining the shell runner's more specific SKIP group disposition.
        skipped_root = root / "skipped"
        skipped_root.mkdir()
        (skipped_root / "summary.txt").write_text("skipped summary\n", encoding="utf-8")
        skipped_log = skipped_root / "skipped.log"
        skipped_log.write_text("SKIP fixture\n", encoding="utf-8")
        skipped_output = skipped_root / "all_stage_result.json"
        command = base_command(root, skipped_root, skipped_output)
        command[command.index("--expected-overall") + 1] = "FAIL"
        command.extend(("--test", "skipped", "stage0", "SKIP", "0", str(skipped_log)))
        run(command)
        value = json.loads(skipped_output.read_text(encoding="utf-8"))
        assert value["overall_result"] == "FAIL"
        assert value["groups"][0]["result"] == "SKIP"
        assert value["counts"]["groups"]["skipped"] == 1
        cases += 1

    print(f"PASS Stage-10 all-stage result writer cases={cases}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

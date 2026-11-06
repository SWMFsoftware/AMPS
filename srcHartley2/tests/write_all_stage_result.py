#!/usr/bin/env python3
"""Write the checksum-bound machine-readable all-stage test receipt.

The shell orchestrator remains responsible for executing tests and retaining
their human-readable logs.  This helper owns JSON encoding, path containment,
checksum calculation, and consistency checks that are awkward to implement
correctly in shell.  In particular, it keeps focused implementation evidence
separate from the optional Stage-9 long-run science decision.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence


SCHEMA = "amps.comet.all-stage-test-result.v1"
SCIENCE_TEST_ID = "stage9.science-acceptance"
RESULTS = {"PASS", "SKIP", "FAIL"}
REQUIRED_MODES = {f"stage{stage}" for stage in range(4, 10)}
OBJECT_ID = re.compile(r"^[0-9a-f]{40}$")


class ResultError(ValueError):
    """Raised when the shell-provided evidence is incomplete or inconsistent."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _contained_file(path: Path, root: Path, label: str) -> tuple[Path, str]:
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        raise ResultError(f"cannot resolve {label}: {path}: {exc}") from exc
    if not resolved.is_file():
        raise ResultError(f"{label} is not a regular file: {path}")
    try:
        relative = resolved.relative_to(root)
    except ValueError as exc:
        raise ResultError(f"{label} escapes the log root: {path}") from exc
    return resolved, relative.as_posix()


def _counts(results: Iterable[str]) -> Dict[str, int]:
    values = list(results)
    return {
        "run": len(values),
        "passed": values.count("PASS"),
        "skipped": values.count("SKIP"),
        "failed": values.count("FAIL"),
    }


def _aggregate_result(counts: Mapping[str, int]) -> str:
    # The all-stage contract is deliberately stricter than many test runners:
    # a skipped required test is not an implementation PASS.
    return "PASS" if counts["failed"] == 0 and counts["skipped"] == 0 else "FAIL"


def _group_result(counts: Mapping[str, int]) -> str:
    """Mirror the shell summary's PASS/SKIP/FAIL group classification."""

    if counts["failed"]:
        return "FAIL"
    if counts["skipped"]:
        return "SKIP"
    return "PASS"


def build_result(
    *,
    repository_root: Path,
    head_commit: str,
    head_tree: str,
    branch: str | None,
    modes: Mapping[str, str],
    raw_tests: Sequence[Sequence[str]],
    log_root: Path,
    summary: Path,
    expected_overall: str,
    science_manifest: Path | None,
) -> Dict[str, Any]:
    """Validate shell records and return one schema-v1 receipt object."""

    repository_root = repository_root.resolve(strict=True)
    if not repository_root.is_dir():
        raise ResultError(f"repository root is not a directory: {repository_root}")
    log_root = log_root.resolve(strict=True)
    if not log_root.is_dir():
        raise ResultError(f"log root is not a directory: {log_root}")
    if not OBJECT_ID.fullmatch(head_commit):
        raise ResultError("HEAD commit is not a lowercase 40-character object ID")
    if not OBJECT_ID.fullmatch(head_tree):
        raise ResultError("HEAD tree is not a lowercase 40-character object ID")
    if branch is not None and not branch:
        raise ResultError("branch must be nonempty or null for detached HEAD")

    if set(modes) != REQUIRED_MODES:
        missing = sorted(REQUIRED_MODES - set(modes))
        extra = sorted(set(modes) - REQUIRED_MODES)
        raise ResultError(f"mode keys differ: missing={missing} extra={extra}")
    for stage, mode in modes.items():
        if mode not in {"archived", "fresh"}:
            raise ResultError(f"invalid mode for {stage}: {mode}")

    tests: list[Dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, raw in enumerate(raw_tests, start=1):
        if len(raw) != 5:
            raise ResultError(f"test record {index} does not contain five fields")
        test_id, stage, result, raw_exit_code, raw_log = raw
        if not test_id or not stage:
            raise ResultError(f"test record {index} has an empty ID or stage")
        if test_id in seen_ids:
            raise ResultError(f"duplicate test ID: {test_id}")
        seen_ids.add(test_id)
        if result not in RESULTS:
            raise ResultError(f"invalid result for {test_id}: {result}")
        try:
            exit_code = int(raw_exit_code, 10)
        except ValueError as exc:
            raise ResultError(f"invalid exit code for {test_id}: {raw_exit_code}") from exc
        if exit_code < 0:
            raise ResultError(f"negative exit code for {test_id}: {exit_code}")
        if result in {"PASS", "SKIP"} and exit_code != 0:
            raise ResultError(f"{result} test {test_id} has nonzero exit code")
        if result == "FAIL" and exit_code == 0:
            raise ResultError(f"failed test {test_id} has zero exit code")
        resolved_log, relative_log = _contained_file(
            Path(raw_log), log_root, f"log for {test_id}"
        )
        tests.append(
            {
                "id": test_id,
                "stage": stage,
                "result": result,
                "exit_code": exit_code,
                "log": relative_log,
                "log_sha256": sha256_file(resolved_log),
            }
        )
    if not tests:
        raise ResultError("at least one test record is required")

    summary_path, relative_summary = _contained_file(summary, log_root, "summary")

    # Ordered insertion follows test execution order, matching summary.txt
    # without imposing lexical assumptions on future stage identifiers.
    grouped: "OrderedDict[str, list[Dict[str, Any]]]" = OrderedDict()
    for test in tests:
        grouped.setdefault(test["stage"], []).append(test)
    groups: list[Dict[str, Any]] = []
    for stage, stage_tests in grouped.items():
        counts = _counts(test["result"] for test in stage_tests)
        groups.append(
            {
                "id": stage,
                "result": _group_result(counts),
                "counts": counts,
            }
        )

    test_counts = _counts(test["result"] for test in tests)
    group_counts = _counts(group["result"] for group in groups)
    overall = _aggregate_result(test_counts)
    if expected_overall not in {"PASS", "FAIL"}:
        raise ResultError(f"invalid expected overall result: {expected_overall}")
    if expected_overall != overall:
        raise ResultError(
            f"shell/Python overall result mismatch: shell={expected_overall} computed={overall}"
        )

    science_tests = [test for test in tests if test["id"] == SCIENCE_TEST_ID]
    if science_manifest is None:
        if science_tests:
            raise ResultError("science-acceptance test exists without a science manifest")
        science = {
            "result": "NOT_EVALUATED",
            "accepted": False,
            "manifest": None,
            "manifest_sha256": None,
            "meaning": (
                "Stage-9 long-run science was not evaluated; focused and archived "
                "implementation validation is not scientific acceptance."
            ),
        }
    else:
        if len(science_tests) != 1:
            raise ResultError(
                "a science manifest requires exactly one stage9.science-acceptance test"
            )
        manifest = science_manifest.resolve(strict=True)
        if not manifest.is_file():
            raise ResultError(f"science manifest is not a regular file: {manifest}")
        accepted = science_tests[0]["result"] == "PASS"
        science = {
            "result": "PASS" if accepted else "FAIL_OR_INVALID",
            "accepted": accepted,
            "manifest": str(manifest),
            "manifest_sha256": sha256_file(manifest),
            "meaning": (
                "The supplied checksum-protected Stage-9 manifest passed the "
                "independent-seed science gate."
                if accepted
                else "The supplied Stage-9 science evidence failed or was invalid."
            ),
        }

    focused_tests = [test for test in tests if test["id"] != SCIENCE_TEST_ID]
    focused_counts = _counts(test["result"] for test in focused_tests)
    focused_result = _aggregate_result(focused_counts)

    return {
        "schema": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z"),
        "repository": {
            "root": str(repository_root),
            "head_commit": head_commit,
            "head_tree": head_tree,
            "branch": branch,
            "detached_head": branch is None,
        },
        "modes": dict(sorted(modes.items())),
        "tests": tests,
        "groups": groups,
        "counts": {"tests": test_counts, "groups": group_counts},
        "summary": {
            "path": relative_summary,
            "sha256": sha256_file(summary_path),
        },
        "interpretation": {
            "focused_implementation": {
                "result": focused_result,
                "accepted": focused_result == "PASS",
                "counts": focused_counts,
                "meaning": (
                    "Focused, parser, synthetic, and archived/fresh implementation "
                    "checks only; this result is not Stage-9 science acceptance."
                ),
            },
            "stage9_science": science,
        },
        "overall_result": overall,
        "overall_scope": (
            "focused_implementation_and_stage9_science"
            if science_manifest is not None
            else "focused_implementation_only"
        ),
    }


def atomic_write(path: Path, value: Mapping[str, Any]) -> None:
    """Publish one complete receipt without overwriting prior evidence."""

    if path.exists():
        raise ResultError(f"refusing to overwrite all-stage result: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(rendered)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", required=True, type=Path)
    parser.add_argument("--head-commit", required=True)
    parser.add_argument("--head-tree", required=True)
    identity = parser.add_mutually_exclusive_group(required=True)
    identity.add_argument("--branch")
    identity.add_argument("--detached-head", action="store_true")
    parser.add_argument("--log-root", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-overall", required=True, choices=("PASS", "FAIL"))
    parser.add_argument(
        "--mode",
        action="append",
        nargs=2,
        metavar=("STAGE", "MODE"),
        default=[],
    )
    parser.add_argument(
        "--test",
        action="append",
        nargs=5,
        metavar=("ID", "STAGE", "RESULT", "EXIT", "LOG"),
        default=[],
    )
    parser.add_argument("--science-manifest", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    modes: Dict[str, str] = {}
    for stage, mode in arguments.mode:
        if stage in modes:
            print(f"write_all_stage_result.py: duplicate mode: {stage}", file=os.sys.stderr)
            return 2
        modes[stage] = mode
    try:
        value = build_result(
            repository_root=arguments.repository_root,
            head_commit=arguments.head_commit,
            head_tree=arguments.head_tree,
            branch=None if arguments.detached_head else arguments.branch,
            modes=modes,
            raw_tests=arguments.test,
            log_root=arguments.log_root,
            summary=arguments.summary,
            expected_overall=arguments.expected_overall,
            science_manifest=arguments.science_manifest,
        )
        output = arguments.output.resolve()
        log_root = arguments.log_root.resolve(strict=True)
        if output.parent != log_root:
            raise ResultError("all-stage result output must be directly inside the log root")
        atomic_write(output, value)
    except (OSError, ResultError) as exc:
        print(f"write_all_stage_result.py: {exc}", file=os.sys.stderr)
        return 2
    print(
        f"PASS all-stage JSON result output={output} "
        f"tests={value['counts']['tests']['run']} groups={value['counts']['groups']['run']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

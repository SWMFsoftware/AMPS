#!/usr/bin/env python3
"""Collect or verify a fail-closed Stage-10 merge-readiness bundle.

This program prepares evidence for review; it never grants permission to
merge and never changes a ref, the index, or the worktree.  ``collect``
requires the caller to spell out the exact 40-hex target, merge-base, and
development commits.  The current branch and HEAD, the target refs, the true
merge base, the tracked worktree state, and the focused-test evidence are all
checked independently before a new output directory is created.

Only the read-only git invocations enumerated by :class:`ReadOnlyGit` are
available to collection.  Keeping construction and validation in one small
wrapper makes a future accidental ``checkout``, ``merge``, ``commit``, or
index-writing command fail before git is executed.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple


SCHEMA = "amps.comet.stage10-merge-readiness.v1"
POLICY_SCHEMA = "amps.comet.stage10-release-policy.v1"
STATUS_FOCUSED = "FOCUSED_ONLY_NOT_SCIENCE_READY"
STATUS_READY = "READY_FOR_USER_REVIEW"
AUTHORIZATION = "NOT_GRANTED"
FULL_HASH = re.compile(r"[0-9a-f]{40}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
THIS_DIR = Path(__file__).resolve().parent
DEFAULT_POLICY = THIS_DIR / "stage10_release_policy.json"
BUNDLE_FILES = (
    "readiness.json",
    "readiness.md",
    "release_policy.json",
    "SHA256SUMS",
)


class ReadinessError(RuntimeError):
    """The requested report is ambiguous, incomplete, or unauthenticated."""


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path, label: str) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ReadinessError(f"cannot read {label} JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ReadinessError(f"{label} JSON root is not an object: {path}")
    return value


def exact_keys(value: Mapping[str, Any], expected: Iterable[str], label: str) -> None:
    actual, wanted = set(value), set(expected)
    if actual != wanted:
        raise ReadinessError(
            f"{label} keys differ: missing={sorted(wanted-actual)} "
            f"extra={sorted(actual-wanted)}"
        )


def full_hash(value: str, label: str) -> str:
    if FULL_HASH.fullmatch(value) is None:
        raise ReadinessError(f"{label} must be an explicit lowercase 40-hex commit")
    return value


def expected_sha(value: str, label: str) -> str:
    if SHA256.fullmatch(value) is None:
        raise ReadinessError(f"{label} must be an explicit lowercase SHA-256")
    return value


def authenticated_file(path: Path, digest: str, label: str) -> Dict[str, Any]:
    expected_sha(digest, f"{label} digest")
    if path.is_symlink():
        raise ReadinessError(f"{label} must not be a symlink: {path}")
    resolved = path.resolve(strict=True)
    if not resolved.is_file() or resolved.is_symlink():
        raise ReadinessError(f"{label} must be a regular non-symlink file: {resolved}")
    actual = sha256_file(resolved)
    if actual != digest:
        raise ReadinessError(
            f"{label} SHA-256 mismatch: expected={digest} actual={actual}"
        )
    return {"path": str(resolved), "sha256": actual, "size_bytes": resolved.stat().st_size}


def validate_policy(value: Dict[str, Any]) -> Dict[str, Any]:
    exact_keys(
        value,
        {
            "schema",
            "development_branch",
            "allowed_target_refs",
            "allowed_tracked_dirty_exclusions",
            "allowed_untracked_artifact_prefixes",
            "allowed_untracked_artifact_exact_paths",
            "focused_all_stage",
            "science_evidence",
            "path_classification",
            "authorization",
        },
        "Stage-10 policy",
    )
    if value["schema"] != POLICY_SCHEMA:
        raise ReadinessError("unknown Stage-10 policy schema")
    if value["development_branch"] != "1d_comet_dev":
        raise ReadinessError("Stage-10 policy must name 1d_comet_dev")
    if value["authorization"] != AUTHORIZATION:
        raise ReadinessError("policy is not fail-closed for merge authorization")
    for key in (
        "allowed_target_refs",
        "allowed_tracked_dirty_exclusions",
        "allowed_untracked_artifact_prefixes",
        "allowed_untracked_artifact_exact_paths",
    ):
        rows = value[key]
        if not isinstance(rows, list) or not rows or not all(
            isinstance(row, str) and row for row in rows
        ):
            raise ReadinessError(f"policy {key} is not a nonempty string list")
    focused = value["focused_all_stage"]
    exact_keys(
        focused,
        {"expected_groups", "expected_tests", "required_stages"},
        "focused-test policy",
    )
    if focused["expected_tests"] < 1 or focused["expected_groups"] < 1:
        raise ReadinessError("focused-test expected counts must be positive")
    if focused["required_stages"] != [f"stage{i}" for i in range(11)]:
        raise ReadinessError("focused-test policy must require stages 0 through 10")
    science = value["science_evidence"]
    exact_keys(
        science,
        {
            "manifest_schema",
            "minimum_balanced_runs",
            "minimum_independent_seeds",
            "result_schema",
            "rubin_schema",
            "required_rubin_companions",
        },
        "science-evidence policy",
    )
    if science["minimum_balanced_runs"] < 8 or science["minimum_independent_seeds"] < 4:
        raise ReadinessError("science policy weakens the Stage-9 balanced ensemble")
    classification = value["path_classification"]
    required_classification = {
        "core_source_prefixes",
        "input_policy_prefixes",
        "test_prefixes",
        "validation_prefixes",
        "reference_prefixes",
        "visualization_prefixes",
        "application_prefixes",
        "build_configuration_paths",
        "documentation_suffixes",
    }
    exact_keys(classification, required_classification, "path-classification policy")
    for key, rows in classification.items():
        if not isinstance(rows, list) or not all(isinstance(row, str) for row in rows):
            raise ReadinessError(f"classification field {key} is not a string list")
    return value


class ReadOnlyGit:
    """Execute only the exact read-only git command forms used by Stage 10."""

    # This list is also a human-reviewable assertion of what collection can do.
    ALLOWED_SUBCOMMANDS = {
        "rev-parse",
        "merge-base",
        "status",
        "diff",
        "rev-list",
        "show",
    }
    FORBIDDEN_SUBCOMMANDS = {
        "add",
        "am",
        "apply",
        "branch",
        "checkout",
        "cherry-pick",
        "clean",
        "commit",
        "fetch",
        "merge",
        "mv",
        "pull",
        "push",
        "rebase",
        "reset",
        "restore",
        "rm",
        "stash",
        "switch",
        "tag",
        "update-index",
        "update-ref",
    }

    def __init__(self, repository: Path):
        self.repository = repository.resolve(strict=True)

    @staticmethod
    def validate_arguments(arguments: Sequence[str]) -> None:
        if not arguments or arguments[0] not in ReadOnlyGit.ALLOWED_SUBCOMMANDS:
            raise ReadinessError(f"git command is not allowlisted: {list(arguments)!r}")
        command = arguments[0]
        rest = list(arguments[1:])
        valid = False
        if command == "rev-parse":
            valid = rest in (["--show-toplevel"], ["--abbrev-ref", "HEAD"]) or (
                len(rest) == 2
                and rest[0] == "--verify"
                and (
                    rest[1].endswith("^{commit}")
                    or re.fullmatch(r"[0-9a-f]{40}\^\{tree\}", rest[1]) is not None
                )
            )
        elif command == "merge-base":
            valid = len(rest) == 2 and all(FULL_HASH.fullmatch(item) for item in rest)
        elif command == "status":
            valid = rest == ["--porcelain=v1", "-z", "--untracked-files=all"]
        elif command == "diff":
            valid = (
                len(rest) == 5
                and rest[:3] in (
                    ["--name-status", "--no-renames", "-z"],
                    ["--numstat", "--no-renames", "-z"],
                )
                and all(FULL_HASH.fullmatch(item) for item in rest[3:])
            ) or (
                len(rest) == 3
                and rest[0] == "--shortstat"
                and all(FULL_HASH.fullmatch(item) for item in rest[1:])
            )
        elif command == "rev-list":
            valid = len(rest) == 2 and rest[0] == "--reverse" and re.fullmatch(
                r"[0-9a-f]{40}\.\.[0-9a-f]{40}", rest[1]
            ) is not None
        elif command == "show":
            metadata_form = (
                len(rest) == 4
                and rest[:2] == ["-s", "--format=%H%x00%P%x00%an%x00%aI%x00%s"]
                and FULL_HASH.fullmatch(rest[2]) is not None
                and rest[3] == "--"
            )
            blob_form = False
            if len(rest) == 1:
                match = re.fullmatch(r"([0-9a-f]{40}):(.+)", rest[0])
                if match is not None:
                    path = PurePosixPath(match.group(2))
                    blob_form = not path.is_absolute() and ".." not in path.parts
            valid = metadata_form or blob_form
        if not valid:
            raise ReadinessError(f"git command form is not allowlisted: {list(arguments)!r}")

    def run_bytes(self, *arguments: str, check: bool = True) -> bytes:
        self.validate_arguments(arguments)
        environment = dict(os.environ)
        environment.update({"LC_ALL": "C", "GIT_OPTIONAL_LOCKS": "0"})
        completed = subprocess.run(
            ["git", *arguments],
            cwd=self.repository,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if check and completed.returncode != 0:
            detail = completed.stderr.decode("utf-8", errors="replace").strip()
            raise ReadinessError(f"read-only git {' '.join(arguments)} failed: {detail}")
        if not check and completed.returncode != 0:
            return b""
        return completed.stdout

    def run_text(self, *arguments: str, check: bool = True) -> str:
        return self.run_bytes(*arguments, check=check).decode("utf-8", errors="strict").strip()


def resolve_commit(git: ReadOnlyGit, commit: str, label: str) -> str:
    requested = full_hash(commit, label)
    resolved = git.run_text("rev-parse", "--verify", f"{requested}^{{commit}}")
    if resolved != requested:
        raise ReadinessError(f"{label} does not resolve exactly to itself")
    return requested


def parse_status(data: bytes) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    for raw in data.split(b"\0"):
        if not raw:
            continue
        if len(raw) < 4 or raw[2:3] != b" ":
            raise ReadinessError("cannot parse tracked git status record")
        code = raw[:2].decode("ascii", errors="strict")
        path = raw[3:].decode("utf-8", errors="strict")
        if "R" in code or "C" in code or not path:
            raise ReadinessError(f"renamed/copied tracked dirt is not excludable: {path!r}")
        records.append((code, path))
    return records


def parse_dirty_exclusions(values: Sequence[str]) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for value in values:
        path, separator, digest = value.rpartition("=")
        if not separator or not path:
            raise ReadinessError("dirty exclusion must use PATH=SHA256")
        expected_sha(digest, f"dirty exclusion {path}")
        if path in result:
            raise ReadinessError(f"duplicate dirty exclusion: {path}")
        result[path] = digest
    return result


def inspect_worktree(
    git: ReadOnlyGit, repository: Path, dev: str, policy: Mapping[str, Any], supplied: Mapping[str, str]
) -> Dict[str, Any]:
    status = parse_status(
        git.run_bytes("status", "--porcelain=v1", "-z", "--untracked-files=all")
    )
    dirty = [(code, path) for code, path in status if code != "??"]
    untracked = sorted(path for code, path in status if code == "??")
    dirty_paths = {path for _, path in dirty}
    if set(supplied) != dirty_paths:
        raise ReadinessError(
            "tracked dirty exclusions must name every and only tracked dirty path: "
            f"dirty={sorted(dirty_paths)} supplied={sorted(supplied)}"
        )
    allowed = set(policy["allowed_tracked_dirty_exclusions"])
    unexpected = dirty_paths - allowed
    if unexpected:
        raise ReadinessError(f"tracked dirty paths are not excludable: {sorted(unexpected)}")
    exclusions: List[Dict[str, Any]] = []
    for code, path_text in sorted(dirty):
        path = repository / path_text
        if not path.is_file() or path.is_symlink():
            raise ReadinessError(f"dirty exclusion is not a regular file: {path_text}")
        actual = sha256_file(path)
        if actual != supplied[path_text]:
            raise ReadinessError(
                f"dirty exclusion digest mismatch for {path_text}: "
                f"expected={supplied[path_text]} actual={actual}"
            )
        baseline = git.run_bytes("show", f"{dev}:{path_text}")
        exclusions.append(
            {
                "path": path_text,
                "status": code,
                "worktree_sha256": actual,
                "development_commit_blob_sha256": sha256_bytes(baseline),
                "merge_scope": "EXCLUDED_MACHINE_LOCAL_CHANGE",
            }
        )

    # Untracked source is more dangerous than ordinary tracked dirt: it can be
    # compiled or executed by a validation run while remaining absent from the
    # proposed commit range.  Reviewable trees therefore block unconditionally.
    # Only explicit generated/runtime artifact families in the release policy
    # are excluded, and any new top-level family fails closed.
    reviewable = [
        path
        for path in untracked
        if path.startswith(("src/", "input/", "srcHartley2/"))
        or ("/" not in path and path.lower().endswith(".md"))
    ]
    if reviewable:
        raise ReadinessError(
            "untracked reviewable source/test/input/documentation paths block readiness: "
            f"{reviewable}"
        )
    exact_artifacts = set(policy["allowed_untracked_artifact_exact_paths"])
    prefix_artifacts = sorted(
        policy["allowed_untracked_artifact_prefixes"], key=lambda value: (-len(value), value)
    )
    artifact_groups: Dict[str, Dict[str, Any]] = {}
    rejected: List[str] = []
    for path_text in untracked:
        if path_text in exact_artifacts:
            rule = f"exact:{path_text}"
        else:
            prefix = next((value for value in prefix_artifacts if path_text.startswith(value)), None)
            if prefix is None:
                rejected.append(path_text)
                continue
            rule = f"prefix:{prefix}"
        path = repository / path_text
        try:
            size = path.lstat().st_size
        except OSError as exc:
            raise ReadinessError(f"cannot inspect allowed untracked artifact {path_text}: {exc}") from exc
        group = artifact_groups.setdefault(
            rule,
            {"rule": rule, "path_count": 0, "total_lstat_bytes": 0, "sample_paths": []},
        )
        group["path_count"] += 1
        group["total_lstat_bytes"] += size
        if len(group["sample_paths"]) < 10:
            group["sample_paths"].append(path_text)
    if rejected:
        raise ReadinessError(
            "untracked paths are outside explicit generated-artifact exclusions: "
            f"{rejected}"
        )
    return {
        "tracked_dirty_count": len(exclusions),
        "tracked_dirty_exclusions": exclusions,
        "untracked_artifact_exclusions": {
            "path_count": len(untracked),
            "groups": [artifact_groups[key] for key in sorted(artifact_groups)],
            "content_hashing": "OMITTED_FOR_RAW_GENERATED_ARTIFACTS",
            "merge_scope": "EXCLUDED_GENERATED_OR_RUNTIME_ARTIFACTS",
        },
    }


ALL_STAGE_SCHEMA = "amps.comet.all-stage-test-result.v1"
SCIENCE_TEST_ID = "stage9.science-acceptance"


def _receipt_file(log_root: Path, relative: object, digest: object, label: str) -> Dict[str, Any]:
    """Authenticate a receipt member without allowing path escape or aliases."""
    if not isinstance(relative, str) or not relative:
        raise ReadinessError(f"{label} path is not nonempty text")
    if not isinstance(digest, str):
        raise ReadinessError(f"{label} SHA-256 is not text")
    expected_sha(digest, f"{label} SHA-256")
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ReadinessError(f"{label} path is absolute or escapes the log root")
    unresolved = log_root / candidate
    if unresolved.is_symlink():
        raise ReadinessError(f"{label} must not be a symlink")
    resolved = unresolved.resolve(strict=True)
    try:
        resolved.relative_to(log_root)
    except ValueError as exc:
        raise ReadinessError(f"{label} resolves outside the log root") from exc
    if not resolved.is_file() or resolved.is_symlink():
        raise ReadinessError(f"{label} is not a regular non-symlink file")
    actual = sha256_file(resolved)
    if actual != digest:
        raise ReadinessError(f"{label} SHA-256 mismatch: expected={digest} actual={actual}")
    return {"path": relative, "sha256": actual, "size_bytes": resolved.stat().st_size}


def _pass_counts(value: object, expected: int, label: str) -> None:
    if not isinstance(value, dict):
        raise ReadinessError(f"{label} counts are not an object")
    exact_keys(value, {"run", "passed", "skipped", "failed"}, f"{label} counts")
    wanted = {"run": expected, "passed": expected, "skipped": 0, "failed": 0}
    if value != wanted:
        raise ReadinessError(f"{label} is not an exact PASS: {value}")


def parse_all_stage_result(
    path: Path,
    digest: str,
    repository: Path,
    dev: str,
    dev_tree: str,
    policy: Mapping[str, Any],
) -> Dict[str, Any]:
    """Strictly authenticate the receipt and every summary/log it names."""
    receipt_record = authenticated_file(path, digest, "all-stage result")
    receipt = read_json(path, "all-stage result")
    exact_keys(
        receipt,
        {
            "schema", "generated_at_utc", "repository", "modes", "tests",
            "groups", "counts", "summary", "interpretation",
            "overall_result", "overall_scope",
        },
        "all-stage result",
    )
    if receipt["schema"] != ALL_STAGE_SCHEMA:
        raise ReadinessError("unknown all-stage result schema")
    provenance = receipt["repository"]
    if not isinstance(provenance, dict):
        raise ReadinessError("all-stage repository provenance is not an object")
    exact_keys(
        provenance,
        {"root", "head_commit", "head_tree", "branch", "detached_head"},
        "all-stage repository provenance",
    )
    if provenance["head_commit"] != dev:
        raise ReadinessError("all-stage receipt HEAD does not equal --dev-commit")
    if provenance["branch"] != policy["development_branch"] or provenance["detached_head"] is not False:
        raise ReadinessError("all-stage receipt was not produced on the development branch")
    try:
        receipt_repository = Path(provenance["root"]).resolve(strict=True)
    except (TypeError, OSError) as exc:
        raise ReadinessError("all-stage repository root is invalid") from exc
    if receipt_repository != repository:
        raise ReadinessError("all-stage receipt belongs to a different repository")
    if provenance["head_tree"] != dev_tree:
        raise ReadinessError("all-stage receipt HEAD tree does not equal the development commit tree")
    if receipt["overall_result"] != "PASS":
        raise ReadinessError("all-stage overall result did not PASS")
    modes = receipt["modes"]
    if not isinstance(modes, dict) or set(modes) != {f"stage{i}" for i in range(4, 10)}:
        raise ReadinessError("all-stage modes do not cover exactly Stages 4--9")
    if any(mode not in {"archived", "fresh"} for mode in modes.values()):
        raise ReadinessError("all-stage receipt contains an invalid evidence mode")

    focused_policy = policy["focused_all_stage"]
    interpretation = receipt["interpretation"]
    if not isinstance(interpretation, dict):
        raise ReadinessError("all-stage interpretation is not an object")
    exact_keys(interpretation, {"focused_implementation", "stage9_science"}, "all-stage interpretation")
    focused = interpretation["focused_implementation"]
    if not isinstance(focused, dict) or focused.get("result") != "PASS" or focused.get("accepted") is not True:
        raise ReadinessError("focused implementation interpretation did not PASS")
    exact_keys(
        focused,
        {"result", "accepted", "counts", "meaning"},
        "focused implementation interpretation",
    )
    _pass_counts(focused.get("counts"), focused_policy["expected_tests"], "focused implementation")
    tests = receipt["tests"]
    if not isinstance(tests, list) or not tests:
        raise ReadinessError("all-stage receipt has no tests")
    log_root = path.resolve().parent
    authenticated_logs: List[Dict[str, Any]] = []
    seen_ids: set[str] = set()
    seen_paths: set[str] = set()
    stage_test_counts: Dict[str, int] = {}
    focused_test_count = 0
    science_test_count = 0
    for row in tests:
        if not isinstance(row, dict):
            raise ReadinessError("all-stage test record is not an object")
        exact_keys(row, {"id", "stage", "result", "exit_code", "log", "log_sha256"}, "all-stage test record")
        test_id, stage = row["id"], row["stage"]
        if not isinstance(test_id, str) or not test_id or test_id in seen_ids:
            raise ReadinessError(f"all-stage test ID is empty or duplicated: {test_id!r}")
        if stage not in focused_policy["required_stages"]:
            raise ReadinessError(f"all-stage test has unexpected stage: {stage!r}")
        if row["result"] != "PASS" or row["exit_code"] != 0:
            raise ReadinessError(f"all-stage test did not PASS exactly: {test_id}")
        if row["log"] in seen_paths:
            raise ReadinessError(f"multiple tests reference one log path: {row['log']}")
        seen_ids.add(test_id)
        seen_paths.add(row["log"])
        stage_test_counts[stage] = stage_test_counts.get(stage, 0) + 1
        authenticated = _receipt_file(log_root, row["log"], row["log_sha256"], f"log for {test_id}")
        authenticated.update({"id": test_id, "stage": stage})
        authenticated_logs.append(authenticated)
        if test_id == SCIENCE_TEST_ID:
            science_test_count += 1
        else:
            focused_test_count += 1
    if focused_test_count != focused_policy["expected_tests"]:
        raise ReadinessError(
            f"receipt has {focused_test_count} focused tests; expected {focused_policy['expected_tests']}"
        )
    if set(stage_test_counts) != set(focused_policy["required_stages"]):
        raise ReadinessError("receipt tests do not cover exactly Stages 0--10")
    if science_test_count not in {0, 1}:
        raise ReadinessError("receipt has multiple Stage-9 science acceptance tests")

    groups = receipt["groups"]
    if not isinstance(groups, list) or len(groups) != focused_policy["expected_groups"]:
        raise ReadinessError("all-stage group count differs from policy")
    group_ids: List[str] = []
    for group in groups:
        if not isinstance(group, dict):
            raise ReadinessError("all-stage group record is not an object")
        exact_keys(group, {"id", "result", "counts"}, "all-stage group")
        if group["id"] in group_ids or group["result"] != "PASS":
            raise ReadinessError("all-stage group is duplicate or did not PASS")
        group_ids.append(group["id"])
        _pass_counts(group["counts"], stage_test_counts.get(group["id"], -1), f"group {group['id']}")
    if group_ids != focused_policy["required_stages"]:
        raise ReadinessError("all-stage groups are not ordered exactly as Stages 0--10")
    counts = receipt["counts"]
    if not isinstance(counts, dict):
        raise ReadinessError("all-stage aggregate counts are not an object")
    exact_keys(counts, {"tests", "groups"}, "all-stage aggregate counts")
    _pass_counts(counts["tests"], len(tests), "all-stage aggregate tests")
    _pass_counts(counts["groups"], focused_policy["expected_groups"], "all-stage aggregate groups")
    summary = receipt["summary"]
    if not isinstance(summary, dict):
        raise ReadinessError("all-stage summary record is not an object")
    exact_keys(summary, {"path", "sha256"}, "all-stage summary record")
    summary_record = _receipt_file(log_root, summary["path"], summary["sha256"], "all-stage summary")
    receipt_science = interpretation["stage9_science"]
    if not isinstance(receipt_science, dict):
        raise ReadinessError("all-stage Stage-9 science interpretation is not an object")
    exact_keys(
        receipt_science,
        {"result", "accepted", "manifest", "manifest_sha256", "meaning"},
        "all-stage Stage-9 science interpretation",
    )
    if science_test_count == 0:
        if receipt_science.get("result") != "NOT_EVALUATED" or receipt_science.get("accepted") is not False:
            raise ReadinessError("receipt conflates focused PASS with Stage-9 science acceptance")
        if receipt["overall_scope"] != "focused_implementation_only":
            raise ReadinessError("receipt overall scope falsely claims Stage-9 science")
    elif receipt_science.get("result") != "PASS" or receipt_science.get("accepted") is not True:
        raise ReadinessError("science acceptance test and interpretation disagree")
    elif receipt["overall_scope"] != "focused_implementation_and_stage9_science":
        raise ReadinessError("receipt overall scope omits its Stage-9 science test")
    receipt_record.update(
        {
            "outcome": "PASS",
            "tests": focused_test_count,
            "groups": len(groups),
            "modes": modes,
            "group_results": groups,
            "test_logs": authenticated_logs,
            "summary": summary_record,
            "receipt_stage9_science": receipt_science,
        }
    )
    return receipt_record


def parse_diff_name_status(data: bytes) -> List[Tuple[str, str]]:
    fields = [field for field in data.split(b"\0") if field]
    if len(fields) % 2:
        raise ReadinessError("cannot parse git name-status output")
    result = []
    for index in range(0, len(fields), 2):
        status = fields[index].decode("ascii", errors="strict")
        path = fields[index + 1].decode("utf-8", errors="strict")
        if status not in {"A", "M", "D", "T", "U"} or not path:
            raise ReadinessError(f"unsupported diff entry: {status!r} {path!r}")
        result.append((status, path))
    return result


def parse_numstat(data: bytes) -> Dict[str, Dict[str, int | str]]:
    result: Dict[str, Dict[str, int | str]] = {}
    for raw in data.split(b"\0"):
        if not raw:
            continue
        text = raw.decode("utf-8", errors="strict")
        fields = text.split("\t", 2)
        if len(fields) != 3:
            raise ReadinessError("cannot parse git numstat output")
        added, deleted, path = fields
        if path in result:
            raise ReadinessError(f"duplicate numstat path: {path}")
        result[path] = {
            "added_lines": int(added) if added.isdigit() else "BINARY",
            "deleted_lines": int(deleted) if deleted.isdigit() else "BINARY",
        }
    return result


def classify_path(path: str, policy: Mapping[str, Any]) -> str:
    rules = policy["path_classification"]
    if path in rules["build_configuration_paths"] or path.endswith("/makefile"):
        return "build_configuration"
    if any(path.startswith(prefix) for prefix in rules["test_prefixes"]):
        return "tests"
    if any(path.startswith(prefix) for prefix in rules["validation_prefixes"]):
        return "validation"
    if any(path.startswith(prefix) for prefix in rules["reference_prefixes"]):
        return "reference_data"
    if any(path.startswith(prefix) for prefix in rules["visualization_prefixes"]):
        return "visualization"
    if any(path.startswith(prefix) for prefix in rules["core_source_prefixes"]):
        return "core_source"
    if any(path.startswith(prefix) for prefix in rules["input_policy_prefixes"]):
        return "input_policy"
    if any(path.endswith(suffix) for suffix in rules["documentation_suffixes"]):
        return "documentation"
    if any(path.startswith(prefix) for prefix in rules["application_prefixes"]):
        return "application_source"
    return "other_review_required"


def collect_diff(
    git: ReadOnlyGit, left: str, right: str, policy: Mapping[str, Any]
) -> Dict[str, Any]:
    names = parse_diff_name_status(
        git.run_bytes("diff", "--name-status", "--no-renames", "-z", left, right)
    )
    numbers = parse_numstat(
        git.run_bytes("diff", "--numstat", "--no-renames", "-z", left, right)
    )
    if {path for _, path in names} != set(numbers):
        raise ReadinessError("git name-status and numstat path sets differ")
    rows = []
    counts: Dict[str, int] = {}
    for status, path in names:
        category = classify_path(path, policy)
        counts[category] = counts.get(category, 0) + 1
        rows.append(
            {
                "path": path,
                "status": status,
                "classification": category,
                **numbers[path],
            }
        )
    unknown = [row["path"] for row in rows if row["classification"] == "other_review_required"]
    if unknown:
        raise ReadinessError(
            "diff contains paths without an approved semantic classification: "
            f"{unknown}"
        )
    return {
        "left_commit": left,
        "right_commit": right,
        "shortstat": git.run_text("diff", "--shortstat", left, right),
        "path_count": len(rows),
        "classification_counts": dict(sorted(counts.items())),
        "paths": rows,
    }


def commit_record(git: ReadOnlyGit, commit: str) -> Dict[str, Any]:
    raw = git.run_bytes(
        "show", "-s", "--format=%H%x00%P%x00%an%x00%aI%x00%s", commit, "--"
    ).rstrip(b"\n")
    fields = raw.decode("utf-8", errors="strict").split("\0")
    if len(fields) != 5 or fields[0] != commit:
        raise ReadinessError(f"cannot parse commit metadata for {commit}")
    return {
        "commit": fields[0],
        "parents": fields[1].split() if fields[1] else [],
        "author": fields[2],
        "authored_at": fields[3],
        "subject": fields[4],
    }


def commit_range(git: ReadOnlyGit, left: str, right: str) -> List[Dict[str, Any]]:
    text = git.run_text("rev-list", "--reverse", f"{left}..{right}")
    commits = [] if not text else text.splitlines()
    if not all(FULL_HASH.fullmatch(commit) for commit in commits):
        raise ReadinessError("git rev-list returned a non-full commit id")
    return [commit_record(git, commit) for commit in commits]


def validate_science_evidence(
    arguments: argparse.Namespace, policy: Mapping[str, Any]
) -> Tuple[Dict[str, Any] | None, List[str]]:
    science_paths = (arguments.stage9_science_manifest, arguments.stage9_science_result)
    science_hashes = (
        arguments.stage9_science_manifest_sha256,
        arguments.stage9_science_result_sha256,
    )
    any_science = any(value is not None for value in (*science_paths, *science_hashes))
    if any_science and not all(value is not None for value in (*science_paths, *science_hashes)):
        raise ReadinessError("Stage-9 science evidence requires manifest/result paths and hashes")
    any_rubin = arguments.rubin_comparison is not None or arguments.rubin_comparison_sha256 is not None
    if any_rubin and (
        arguments.rubin_comparison is None or arguments.rubin_comparison_sha256 is None
    ):
        raise ReadinessError("Rubin evidence requires both comparison path and SHA-256")
    if any_rubin and not any_science:
        raise ReadinessError("Rubin evidence cannot be accepted without Stage-9 science evidence")
    blockers: List[str] = []
    if not any_science:
        blockers.extend(
            ["STAGE9_SCIENCE_EVIDENCE_NOT_SUPPLIED", "RUBIN_COMPARISON_EVIDENCE_NOT_SUPPLIED"]
        )
        return None, blockers

    rules = policy["science_evidence"]
    manifest_record = authenticated_file(
        arguments.stage9_science_manifest,
        arguments.stage9_science_manifest_sha256,
        "Stage-9 science manifest",
    )
    result_record = authenticated_file(
        arguments.stage9_science_result,
        arguments.stage9_science_result_sha256,
        "Stage-9 science result",
    )
    manifest = read_json(arguments.stage9_science_manifest, "Stage-9 science manifest")
    result = read_json(arguments.stage9_science_result, "Stage-9 science result")
    if manifest.get("schema") != rules["manifest_schema"]:
        raise ReadinessError("Stage-9 manifest schema is not accepted")
    runs = manifest.get("runs")
    if not isinstance(runs, list) or len(runs) < rules["minimum_balanced_runs"]:
        raise ReadinessError("Stage-9 manifest lacks the minimum balanced runs")
    if result.get("schema") != rules["result_schema"]:
        raise ReadinessError("Stage-9 result schema is not accepted")
    # A caller-provided checksum authenticates bytes, not scientific truth.
    # Re-run the existing strict evaluator against the checksum-protected
    # manifest, then require byte-decoded semantic equality with its canonical
    # result.  A hand-written PASS-shaped JSON object therefore cannot promote
    # a deficient or incomplete campaign to merge review.
    evaluator_path = THIS_DIR.parent / "tests" / "test_stage9_science_gate.py"
    specification = importlib.util.spec_from_file_location(
        "stage10_strict_stage9_science_gate", evaluator_path
    )
    if specification is None or specification.loader is None:
        raise ReadinessError("cannot load the strict Stage-9 science evaluator")
    evaluator = importlib.util.module_from_spec(specification)
    try:
        specification.loader.exec_module(evaluator)
        computed_result = evaluator.evaluate(arguments.stage9_science_manifest.resolve())
    except Exception as exc:
        raise ReadinessError(f"strict Stage-9 science re-evaluation failed: {exc}") from exc
    if result != computed_result:
        raise ReadinessError("supplied Stage-9 result differs from strict canonical evaluation")
    if computed_result.get("outcome") != "PASS" or computed_result.get("accepted") is not True:
        raise ReadinessError("strict Stage-9 science evaluation is not an accepted PASS")
    if computed_result.get("failures") != []:
        raise ReadinessError("strict Stage-9 PASS result contains failures")
    seeds = result.get("seeds")
    if not isinstance(seeds, list) or len(set(seeds)) < rules["minimum_independent_seeds"]:
        raise ReadinessError("Stage-9 result lacks independent seeds")
    if result.get("run_count") != len(runs) or len(runs) < 2 * len(set(seeds)):
        raise ReadinessError("Stage-9 result and manifest run counts are inconsistent")
    record: Dict[str, Any] = {
        "manifest": manifest_record,
        "result": result_record,
        "outcome": "PASS",
        "run_count": len(runs),
        "independent_seeds": sorted(set(seeds)),
    }
    if not any_rubin:
        blockers.append("RUBIN_COMPARISON_EVIDENCE_NOT_SUPPLIED")
        record["rubin_comparison"] = None
        return record, blockers

    rubin_record = authenticated_file(
        arguments.rubin_comparison,
        arguments.rubin_comparison_sha256,
        "Rubin comparison JSON",
    )
    rubin = read_json(arguments.rubin_comparison, "Rubin comparison")
    if rubin.get("schema") != rules["rubin_schema"] or rubin.get("stage9_outcome") != "PASS":
        raise ReadinessError("Rubin comparison is not a Stage-9 PASS product")
    convergence = rubin.get("convergence")
    if not isinstance(convergence, dict) or convergence.get(
        "all_runs_converged_before_hard_limit"
    ) is not True:
        raise ReadinessError("Rubin comparison does not prove convergence before the hard limit")
    if rubin.get("run_count") != result["run_count"] or rubin.get("seeds") != result["seeds"]:
        raise ReadinessError("Rubin comparison and Stage-9 science result populations differ")
    provenance = rubin.get("provenance")
    if not isinstance(provenance, dict) or provenance.get("manifest_sha256") != manifest_record["sha256"]:
        raise ReadinessError("Rubin comparison is not bound to the supplied science manifest")
    profiles = rubin.get("profiles")
    if not isinstance(profiles, dict) or not profiles:
        raise ReadinessError("Rubin comparison has no model/paper profiles")
    companions = []
    for suffix in rules["required_rubin_companions"]:
        companion = arguments.rubin_comparison.with_suffix(suffix).resolve(strict=True)
        if not companion.is_file() or companion.is_symlink():
            raise ReadinessError(f"Rubin comparison companion is invalid: {companion}")
        companions.append(
            {"path": str(companion), "sha256": sha256_file(companion), "size_bytes": companion.stat().st_size}
        )
    rubin_record["companions"] = companions
    record["rubin_comparison"] = rubin_record
    return record, blockers


def render_markdown(report: Mapping[str, Any]) -> str:
    lines = [
        "# Hartley2 Stage-10 merge-readiness report",
        "",
        f"- Status: `{report['status']}`",
        f"- Merge authorization: `{report['authorization']['state']}`",
        f"- Development branch: `{report['repository']['branch']}`",
        f"- Requested target ref: `{report['commits']['requested_target_ref']}`",
        f"- Target commit: `{report['commits']['target']}`",
        f"- Merge base: `{report['commits']['base']}`",
        f"- Development commit: `{report['commits']['development']}`",
        f"- Focused tests: {report['focused_tests']['tests']}/{report['focused_tests']['tests']} "
        f"in {report['focused_tests']['groups']} groups",
        "",
        "## Blocking conditions",
        "",
    ]
    blockers = report["blockers"]
    lines.extend([f"- `{blocker}`" for blocker in blockers] if blockers else ["- None in the supplied evidence."])
    lines.extend(["", "## Exact integration diff", ""])
    integration = report["diffs"]["target_to_development"]
    lines.append(f"- Changed paths: {integration['path_count']}")
    lines.append(f"- Diffstat: {integration['shortstat'] or 'no changes'}")
    for category, count in integration["classification_counts"].items():
        lines.append(f"- `{category}`: {count}")
    lines.extend(["", "## Tracked worktree exclusions", ""])
    exclusions = report["worktree"]["tracked_dirty_exclusions"]
    if exclusions:
        for row in exclusions:
            lines.append(
                f"- `{row['path']}` status `{row['status']}`, worktree SHA-256 "
                f"`{row['worktree_sha256']}`; excluded from merge scope"
            )
    else:
        lines.append("- None")
    artifacts = report["worktree"]["untracked_artifact_exclusions"]
    lines.extend(["", "## Untracked generated-artifact exclusions", ""])
    lines.append(f"- Excluded paths: {artifacts['path_count']}")
    for group in artifacts["groups"]:
        lines.append(
            f"- `{group['rule']}`: {group['path_count']} paths, "
            f"{group['total_lstat_bytes']} lstat bytes"
        )
    lines.extend(
        [
            "",
            "## Safety boundary",
            "",
            "This bundle is read-only preparation. Even `READY_FOR_USER_REVIEW` does not grant "
            "permission to merge, commit, switch, reset, or modify `master`. A separate explicit "
            "user authorization is required.",
            "",
        ]
    )
    return "\n".join(lines)


def write_bundle(output: Path, policy_path: Path, report: Dict[str, Any]) -> None:
    output = output.resolve()
    if output.exists():
        raise ReadinessError(f"output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    try:
        policy_bytes = policy_path.read_bytes()
        (output / "release_policy.json").write_bytes(policy_bytes)
        json_text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        (output / "readiness.json").write_text(json_text, encoding="utf-8")
        (output / "readiness.md").write_text(render_markdown(report), encoding="utf-8")
        checksum_names = ("readiness.json", "readiness.md", "release_policy.json")
        checksum_text = "".join(
            f"{sha256_file(output / name)}  {name}\n" for name in checksum_names
        )
        (output / "SHA256SUMS").write_text(checksum_text, encoding="ascii")
    except BaseException:
        shutil.rmtree(output)
        raise


def collect(arguments: argparse.Namespace) -> int:
    repository = arguments.repo.resolve(strict=True)
    policy_path = arguments.policy.resolve(strict=True)
    policy_record = authenticated_file(policy_path, sha256_file(policy_path), "Stage-10 policy")
    policy = validate_policy(read_json(policy_path, "Stage-10 policy"))
    canonical_policy = validate_policy(read_json(DEFAULT_POLICY, "canonical Stage-10 policy"))
    if policy != canonical_policy:
        raise ReadinessError("supplied Stage-10 policy differs from the canonical release policy")
    git = ReadOnlyGit(repository)
    root = Path(git.run_text("rev-parse", "--show-toplevel")).resolve(strict=True)
    if root != repository:
        raise ReadinessError(f"--repo must be the repository root: actual={root}")
    target = resolve_commit(git, arguments.target_commit, "target commit")
    base = resolve_commit(git, arguments.base_commit, "base commit")
    dev = resolve_commit(git, arguments.dev_commit, "development commit")
    branch = git.run_text("rev-parse", "--abbrev-ref", "HEAD")
    if branch != policy["development_branch"]:
        raise ReadinessError(
            f"collection requires branch {policy['development_branch']}; current={branch}"
        )
    head = git.run_text("rev-parse", "--verify", "HEAD^{commit}")
    if head != dev:
        raise ReadinessError(f"development commit must equal current HEAD: head={head} supplied={dev}")
    dev_tree = git.run_text("rev-parse", "--verify", f"{dev}^{{tree}}")
    if FULL_HASH.fullmatch(dev_tree) is None:
        raise ReadinessError("development commit did not resolve to a full tree object ID")
    if arguments.target_ref not in policy["allowed_target_refs"]:
        raise ReadinessError("--target-ref is not an exact policy-allowed target ref")
    matched_refs = []
    target_ref_state = []
    for ref in policy["allowed_target_refs"]:
        resolved = git.run_text("rev-parse", "--verify", f"{ref}^{{commit}}", check=False)
        target_ref_state.append({"ref": ref, "commit": resolved or None})
        if resolved == target:
            matched_refs.append(ref)
    if not matched_refs:
        raise ReadinessError("target commit does not equal an allowed master target ref")
    requested_ref_commit = next(
        row["commit"] for row in target_ref_state if row["ref"] == arguments.target_ref
    )
    if requested_ref_commit != target:
        raise ReadinessError("--target-ref does not resolve to --target-commit")
    actual_base = git.run_text("merge-base", target, dev)
    if actual_base != base:
        raise ReadinessError(f"supplied base is not the exact merge base: actual={actual_base}")

    exclusions = parse_dirty_exclusions(arguments.dirty_exclusion)
    worktree = inspect_worktree(git, repository, dev, policy, exclusions)
    focused = parse_all_stage_result(
        arguments.all_stage_result,
        arguments.all_stage_result_sha256,
        repository,
        dev,
        dev_tree,
        policy,
    )
    science, blockers = validate_science_evidence(arguments, policy)
    status = STATUS_READY if science is not None and science.get("rubin_comparison") else STATUS_FOCUSED
    if status == STATUS_READY and blockers:
        raise ReadinessError("internal status error: ready report still has science blockers")
    report: Dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": status,
        "authorization": {
            "state": AUTHORIZATION,
            "required_next_action": "OBTAIN_SEPARATE_EXPLICIT_USER_AUTHORIZATION_BEFORE_ANY_MASTER_OPERATION",
        },
        "blockers": blockers,
        "policy": policy_record,
        "repository": {"root": str(repository), "branch": branch},
        "commits": {
            "target": target,
            "requested_target_ref": arguments.target_ref,
            "target_matching_refs": matched_refs,
            "target_ref_state": target_ref_state,
            "base": base,
            "development": dev,
            "development_commits_after_base": commit_range(git, base, dev),
            "integration_commits_not_in_target": commit_range(git, target, dev),
            "target_commits_not_in_development": commit_range(git, dev, target),
        },
        "diffs": {
            "base_to_development": collect_diff(git, base, dev, policy),
            "target_to_development": collect_diff(git, target, dev, policy),
        },
        "worktree": worktree,
        "focused_tests": focused,
        "stage9_science": science,
    }
    write_bundle(arguments.output_dir, policy_path, report)
    print(
        f"PASS Stage-10 readiness status={status} authorization={AUTHORIZATION} "
        f"output={arguments.output_dir.resolve()}"
    )
    return 0


def verify(arguments: argparse.Namespace) -> int:
    bundle = arguments.bundle.resolve(strict=True)
    if not bundle.is_dir() or bundle.is_symlink():
        raise ReadinessError(f"bundle is not a regular directory: {bundle}")
    actual_names = sorted(path.name for path in bundle.iterdir())
    if actual_names != sorted(BUNDLE_FILES):
        raise ReadinessError(
            f"bundle members differ: expected={sorted(BUNDLE_FILES)} actual={actual_names}"
        )
    checksum_lines = (bundle / "SHA256SUMS").read_text(encoding="ascii").splitlines()
    expected_names = ["readiness.json", "readiness.md", "release_policy.json"]
    if len(checksum_lines) != len(expected_names):
        raise ReadinessError("SHA256SUMS has the wrong number of records")
    for line, name in zip(checksum_lines, expected_names):
        match = re.fullmatch(r"([0-9a-f]{64})  ([A-Za-z0-9_.-]+)", line)
        if match is None or match.group(2) != name:
            raise ReadinessError(f"invalid or out-of-order checksum record: {line!r}")
        actual = sha256_file(bundle / name)
        if actual != match.group(1):
            raise ReadinessError(f"bundle member checksum mismatch: {name}")
    policy = validate_policy(read_json(bundle / "release_policy.json", "bundled policy"))
    canonical_policy = validate_policy(read_json(DEFAULT_POLICY, "canonical Stage-10 policy"))
    if policy != canonical_policy:
        raise ReadinessError("bundled policy differs from the canonical Stage-10 policy")
    report = read_json(bundle / "readiness.json", "readiness report")
    if report.get("schema") != SCHEMA:
        raise ReadinessError("unknown readiness-report schema")
    if report.get("status") not in {STATUS_FOCUSED, STATUS_READY}:
        raise ReadinessError("readiness report has an invalid status")
    authorization = report.get("authorization")
    if not isinstance(authorization, dict) or authorization.get("state") != AUTHORIZATION:
        raise ReadinessError("readiness report does not preserve NOT_GRANTED authorization")
    policy_record = report.get("policy")
    if not isinstance(policy_record, dict) or policy_record.get("sha256") != sha256_file(
        bundle / "release_policy.json"
    ):
        raise ReadinessError("readiness report is not bound to the bundled policy")
    focused = report.get("focused_tests")
    if not isinstance(focused, dict) or focused.get("outcome") != "PASS":
        raise ReadinessError("readiness report lacks a focused-test PASS")
    if focused.get("tests") != policy["focused_all_stage"]["expected_tests"]:
        raise ReadinessError("readiness report focused-test count differs from policy")
    science = report.get("stage9_science")
    blockers = report.get("blockers")
    if report["status"] == STATUS_READY:
        if not isinstance(science, dict) or not isinstance(science.get("rubin_comparison"), dict):
            raise ReadinessError("ready report lacks Stage-9/Rubin evidence")
        if blockers != []:
            raise ReadinessError("ready report contains blockers")
    else:
        if not isinstance(blockers, list) or not blockers:
            raise ReadinessError("focused-only report lacks explicit science blockers")
    expected_markdown = render_markdown(report)
    if (bundle / "readiness.md").read_text(encoding="utf-8") != expected_markdown:
        raise ReadinessError("readiness Markdown is not the canonical rendering of JSON")
    print(
        f"PASS verified Stage-10 readiness status={report['status']} "
        f"authorization={AUTHORIZATION} bundle={bundle}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Collect or verify a read-only Stage-10 merge-readiness bundle."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    collect_parser = subparsers.add_parser("collect", help="collect a new readiness bundle")
    collect_parser.add_argument("--repo", required=True, type=Path)
    collect_parser.add_argument("--target-commit", required=True)
    collect_parser.add_argument("--target-ref", required=True)
    collect_parser.add_argument("--base-commit", required=True)
    collect_parser.add_argument("--dev-commit", required=True)
    collect_parser.add_argument("--all-stage-result", required=True, type=Path)
    collect_parser.add_argument("--all-stage-result-sha256", required=True)
    collect_parser.add_argument("--output-dir", required=True, type=Path)
    collect_parser.add_argument(
        "--dirty-exclusion",
        action="append",
        default=[],
        metavar="PATH=SHA256",
        help="explicit hash for each policy-allowed tracked dirty file",
    )
    collect_parser.add_argument("--stage9-science-manifest", type=Path)
    collect_parser.add_argument("--stage9-science-manifest-sha256")
    collect_parser.add_argument("--stage9-science-result", type=Path)
    collect_parser.add_argument("--stage9-science-result-sha256")
    collect_parser.add_argument("--rubin-comparison", type=Path)
    collect_parser.add_argument("--rubin-comparison-sha256")
    collect_parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    collect_parser.set_defaults(handler=collect)
    verify_parser = subparsers.add_parser("verify", help="verify a preserved bundle")
    verify_parser.add_argument("--bundle", required=True, type=Path)
    verify_parser.set_defaults(handler=verify)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        return arguments.handler(arguments)
    except (ReadinessError, OSError, UnicodeError, ValueError) as exc:
        print(f"stage10_merge_readiness.py: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

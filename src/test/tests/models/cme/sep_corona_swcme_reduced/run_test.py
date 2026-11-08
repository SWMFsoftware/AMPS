#!/usr/bin/env python3
"""Import all 97 reduced-front rows with scope-qualified identities.

Routine T11 execution passes ``--skip-native``.  That still executes the
complete shared/architecture/application registry and produces 35 explicit
native SKIPs; it never drops native IDs.  A prior native receipt is not reused
as a fresh PASS, and reduced-front evidence can never qualify BG3D-4.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
TEST_ROOT = HERE.parents[4]
sys.path.insert(0, str(TEST_ROOT))
sys.path.insert(0, str(HERE.parents[1]))

from lib.runner_support import case, error_result, load_context, run_command, write_result  # noqa: E402
from common import load_json, prepare_model_repo, selected_rows, suite_status  # noqa: E402


SCOPE_COUNTS = {
    "shared-rsh": 57, "architecture": 1, "portable-srcsep3d": 4,
    "native-smoke-1": 4, "native-smoke-4": 4,
    "native-restart-reference-4": 4, "native-restart-checkpoint-4": 4,
    "native-restart-resume-4": 5, "native-restart-checkpoint-1": 4,
    "native-restart-repartition-4": 5, "native-1au-4": 5,
}


def normalized(scope: str, local_id: str) -> str:
    scope_token = re.sub(r"[^A-Z0-9]+", "-", scope.upper()).strip("-")
    if "#" in local_id:
        base, occurrence = local_id.split("#", 1)
        local_token = "{}-{:02d}".format(base, int(occurrence))
    else:
        local_token = re.sub(r"[^A-Z0-9]+", "-", local_id.upper()).strip("-")
    return "CSWC-REDUCED-{}-{}".format(scope_token, local_token)


def main() -> int:
    context = load_context("Import complete reduced-front plus ambient registry")
    try:
        repo = Path(context["repo_root"])
        work = Path(context["work_dir"])
        artifacts = Path(context["artifact_dir"])
        isolated = prepare_model_repo(repo, work, include_sep3d=True)
        local_output = artifacts / "local-report"
        print("[CASE] RUN reduced-front complete 97-row registry (native rows typed SKIP)", flush=True)
        completed = run_command(
            [sys.executable, "srcSEP3D/test/run_reduced_shock_front.py",
             "--skip-native", "--output-dir", str(local_output)],
            isolated, timeout=3600)
        payload = load_json(local_output / "summary.json")
        if (payload.get("schema") != "srcsep3d-reduced-shock-runner-v1" or
                payload.get("qualification") !=
                "reduced-front-plus-ambient-only; BG3D-4 remains unqualified"):
            raise RuntimeError("reduced report lost schema/qualification authority")
        local_rows = payload.get("results")
        if not isinstance(local_rows, list) or len(local_rows) != 97:
            raise RuntimeError("reduced registry must contain exactly 97 rows")
        counts = {scope: 0 for scope in SCOPE_COUNTS}
        imported = []
        seen = set()
        for row in local_rows:
            if not isinstance(row, dict):
                raise RuntimeError("malformed reduced result row")
            scope = row.get("scope")
            local_id = row.get("test_id")
            status = str(row.get("status", "")).upper()
            if scope not in SCOPE_COUNTS or not isinstance(local_id, str) or status not in {
                    "PASS", "FAIL", "SKIP", "ERROR"}:
                raise RuntimeError("invalid reduced row {!r}".format(row))
            counts[scope] += 1
            global_case = normalized(scope, local_id)
            if global_case in seen:
                raise RuntimeError("duplicate reduced global identity " + global_case)
            seen.add(global_case)
            imported.append(case(
                global_case, status,
                "scope={} local_id={} qualification=reduced-front-plus-ambient-only; {}".format(
                    scope, local_id, row.get("message", ""))))
        if counts != SCOPE_COUNTS:
            raise RuntimeError("reduced scoped cardinality drift: {}".format(counts))
        has_failure = any(row["status"] in ("FAIL", "ERROR") for row in imported)
        if completed.returncode != (1 if has_failure else 0):
            raise RuntimeError("reduced process exit disagrees with parsed statuses")
        imported = selected_rows(context, imported)
        return write_result(
            context, suite_status(imported),
            "complete reduced-front registry imported; 35 native cases remain explicit SKIP in routine mode",
            imported,
            evidence=[str(local_output / "summary.json"),
                      str(local_output / "summary.txt")],
            metrics={"local_registry_cases": 97, "native_cases": 35,
                     "qualification": "reduced-front-plus-ambient-only"})
    except BaseException as error:
        return error_result(context, error)


if __name__ == "__main__":
    raise SystemExit(main())


#!/usr/bin/env python3
"""Import all 222 authoritative staged coronal-CME software checks."""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve()
TEST_ROOT = HERE.parents[4]
sys.path.insert(0, str(TEST_ROOT))
sys.path.insert(0, str(HERE.parents[1]))

from lib.runner_support import case, error_result, load_context, run_command, write_result  # noqa: E402
from common import (exact_rows, global_id, load_json, prepare_model_repo,  # noqa: E402
                    selected_rows, suite_status)


def main() -> int:
    context = load_context("Import the complete sep_coronal_cme registry")
    try:
        repo = Path(context["repo_root"])
        work = Path(context["work_dir"])
        artifacts = Path(context["artifact_dir"])
        isolated = prepare_model_repo(repo, work)
        model = isolated / "src/models/sep_coronal_cme"
        report_dir = artifacts / "local-report"

        # --list is the local registry authority.  It is deliberately queried
        # in the same isolated source copy that is executed; a checkout/list
        # mismatch cannot be hidden by a hard-coded global expected-value set.
        listing = run_command(
            [sys.executable, "test/run_tests.py", "--list"], model, timeout=120)
        if listing.returncode != 0:
            raise RuntimeError("authoritative --list failed with exit {}".format(
                listing.returncode))
        local_ids = [line.split("\t", 1)[0] for line in listing.stdout.splitlines()
                     if "\tstage=" in line]
        if len(local_ids) != 222 or len(set(local_ids)) != 222:
            raise RuntimeError("coronal-CME registry is not the exact 222-case contract")

        print("[CASE] RUN SCCM complete 222-case software-verification registry", flush=True)
        completed = run_command(
            [sys.executable, "test/run_tests.py", "--all", "--expect-count", "222",
             "--output-dir", str(report_dir)], model, timeout=7200)
        payload = load_json(report_dir / "results.json")
        if (payload.get("suite") != "sep_coronal_cme" or
                payload.get("evidence_kind") != "software-verification" or
                payload.get("production_qualified") is not False or
                payload.get("observational_campaign_qualified") is not False):
            raise RuntimeError("local report lost its bounded qualification metadata")
        rows = exact_rows(payload.get("tests", []), local_ids)
        observed_exit = 1 if any(status == "FAIL" for _, status, _ in rows) else 0
        if completed.returncode != observed_exit:
            raise RuntimeError("local process exit disagrees with parsed statuses")

        imported = [case(global_id("SCCM", local_id), status,
                         "local_id={} scope=software-verification; production_qualified=false; {}".format(
                             local_id, message))
                    for local_id, status, message in rows]
        imported = selected_rows(context, imported)
        status = suite_status(imported)
        return write_result(
            context, status,
            "imported exact staged registry; this is software verification, not production or observational qualification",
            imported,
            evidence=[str(report_dir / "results.json"), str(report_dir / "junit.xml")],
            metrics={"local_registry_cases": 222, "qualification": "software-verification"})
    except BaseException as error:
        return error_result(context, error)


if __name__ == "__main__":
    raise SystemExit(main())


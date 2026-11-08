#!/usr/bin/env python3
"""Audit disposition of every legacy application named by the T11 roadmap.

Legacy programs are starting points, not automatically trustworthy tests.
This gate proves that every named directory still exists and has an explicit
global owner and disposition.  It does not execute diagnostic applications or
rename their output as a PASS.  Numerical behavior is exercised by the mapped
production-linked global suites; this audit prevents a legacy directory from
silently disappearing during that migration.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
TEST_ROOT = HERE.parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.manifest import load_manifest  # noqa: E402
from lib.runner_support import case, error_result, load_context, write_result  # noqa: E402


# ``diagnostic-retained`` is intentional for the four wave/solver programs:
# T8 replaced their visual output with independent residual/dispersion
# assertions, but the applications remain useful for manual investigations.
# ``mapped-not-qualified`` means the named global owner contains the current
# production-linked replacement family, but this disposition audit does not
# prove that every diagnostic path in the legacy application was migrated.
# The roadmap still requires per-application independent-oracle, failure-exit,
# deterministic-state and hard-coded-path review before legacy migration can
# be called complete.  That wording prevents directory existence from being
# misread as behavioral qualification.
MIGRATIONS = (
    ("LEGACY-UNIT-TESTS-001", "srcUnitTests", "CORE-PIC-NATIVE-001", "mapped-not-qualified"),
    ("LEGACY-INTERPOLATE-AMRTEST-001", "srcInterpolateAMRTest", "CORE-MESH-NATIVE-001", "mapped-not-qualified"),
    ("LEGACY-PERIODIC-BCTEST-001", "srcPeriodicBCTest", "CORE-PIC-PERIODIC-001", "mapped-not-qualified"),
    ("LEGACY-STENCIL-TEST-001", "srcStencilTest", "CORE-ECSIM-OPERATORS-001", "mapped-not-qualified"),
    ("LEGACY-LINEAR-SYSTEM-SOLVER-001", "srcLinearSystemSolver", "CORE-PIC-FIELD-SOLVER-001", "diagnostic-retained"),
    ("LEGACY-GMRESPERIODIC-BCTEST-001", "srcGMRESPeriodicBCTest", "CORE-PIC-FIELD-SOLVER-001", "diagnostic-retained"),
    ("LEGACY-FAST-WAVE-001", "srcFastWave", "CORE-PIC-FIELD-SOLVER-001", "diagnostic-retained"),
    ("LEGACY-LIGHT-WAVE-TEST-001", "srcLightWaveTest", "CORE-PIC-FIELD-SOLVER-001", "diagnostic-retained"),
    ("LEGACY-LOGGER-001", "srcLogger", "CORE-GENERAL-LOGGER-001", "mapped-not-qualified"),
    ("LEGACY-DATA-READER-TEST-001", "srcDataReaderTest", "CORE-PIC-DATAFILE-001", "diagnostic-retained"),
    ("LEGACY-MODEL-CHEMISTRY-TEST-001", "srcModelChemistryTest", "CORE-CHEMICAL-SOLVER-001", "mapped-not-qualified"),
    ("LEGACY-MODEL-COLLISION-TEST-001", "srcModelCollisionTest", "CORE-PIC-NATIVE-001", "mapped-not-qualified"),
    ("LEGACY-MODEL-DUST-TEST-001", "srcModelDustTest", "CORE-DUST-NATIVE-001", "mapped-not-qualified"),
)


def main() -> int:
    context = load_context("Audit T11 legacy-application migration")
    try:
        repo = Path(context["repo_root"])
        complete_ids = [item[0] for item in MIGRATIONS]
        configured = list(context["test"].get("expected_cases", ()))
        # The audit is useful only if every roadmap application remains
        # globally visible.  Require exact central membership/order before a
        # one-case selection so removing a manifest row cannot silently turn
        # an unreviewed legacy directory into an apparent smaller PASS.
        if configured != complete_ids:
            raise RuntimeError(
                "legacy migration registry drift: configured {} observed {}".format(
                    configured, complete_ids))
        manifest = load_manifest(TEST_ROOT / "manifest.json", repo, TEST_ROOT)
        known = {item["id"] for item in manifest}
        for item in manifest:
            known.update(item.get("aliases", []))
        rows = []
        evidence = []
        for case_id, directory, owner, disposition in MIGRATIONS:
            print("[CASE] RUN {} ({})".format(case_id, directory), flush=True)
            source = repo / "test" / directory
            if not source.is_dir():
                rows.append(case(case_id, "FAIL", "legacy directory is absent: {}".format(source)))
            elif owner not in known:
                rows.append(case(case_id, "ERROR", "mapped global owner is unregistered: " + owner))
            else:
                rows.append(case(
                    case_id, "PASS",
                    "directory=test/{} owner={} disposition={}; no legacy diagnostic was relabeled".format(
                        directory, owner, disposition)))
            evidence.append({"directory": "test/" + directory,
                             "owner": owner, "disposition": disposition})
        selected = list(context["test"].get(
            "selected_cases", configured))
        by_id = {row["id"]: row for row in rows}
        rows = [by_id[item] for item in selected]
        status = "ERROR" if any(row["status"] == "ERROR" for row in rows) else (
            "FAIL" if any(row["status"] == "FAIL" for row in rows) else "PASS")
        receipt = Path(context["artifact_dir"]) / "legacy-migration.json"
        receipt.write_text(json.dumps({"schema_version": 1, "migrations": evidence},
                                      indent=2) + "\n", encoding="utf-8")
        return write_result(
            context, status,
            "audited {} roadmap legacy applications and their bounded global owners".format(
                len(rows)), rows, evidence=[str(receipt)],
            metrics={"roadmap_legacy_applications": len(MIGRATIONS)})
    except BaseException as error:
        return error_result(context, error)


if __name__ == "__main__":
    raise SystemExit(main())


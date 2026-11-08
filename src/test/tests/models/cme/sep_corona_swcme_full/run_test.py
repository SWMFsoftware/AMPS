#!/usr/bin/env python3
"""Import dependency-light full-volume Corona--SWCME code checks.

These tests exercise bounded kernels used by the deferred BG3D-4 program.
Their marker PASS is never promoted into full-CME, G1, production, or
observational qualification.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
TEST_ROOT = HERE.parents[4]
sys.path.insert(0, str(TEST_ROOT))
sys.path.insert(0, str(HERE.parents[1]))

from lib.runner_support import (RunnerError, case, error_result, load_context,  # noqa: E402
                                run_command, write_result)
from common import global_id, prepare_model_repo, selected_rows, suite_status  # noqa: E402


LOCAL_IDS = [
    "CMBGU01", "CMBGU02", "CMBGU03", "CMBGU04", "CMBGU05",
    "CSWC0620", "CSWC0621", "CSWC0622", "CSWC0623", "CSWC0624",
    "CSWC0625", "CSWC0626", "CSWC0627", "CSWC0628-MASS-TIME",
    "ARCHCSWC01",
]

# This is the exact executable registry declared by the model Makefile.  The
# Make ``test`` recipe stops at the first failing executable, which would hide
# later case identities from a global report.  T11 therefore builds those same
# targets and invokes each one in Makefile order.  A zero exit must emit every
# expected PASS marker; a nonzero exit becomes FAIL for that executable's
# owned IDs even when the legacy binary reports only a prose diagnostic.
BINARY_CASES = [
    ("test_bg3d1", ["CMBGU01", "CMBGU02"]),
    ("test_bg3d2", ["CMBGU03"]),
    ("test_bg3d3", ["CMBGU04"]),
    ("test_bg3d4", ["CMBGU05"]),
    ("test_bg3d4_piston_contact", ["CSWC0620"]),
    ("test_bg3d4_planar_piston", ["CSWC0621"]),
    ("test_bg3d4_planar_formation", ["CSWC0622"]),
    ("test_bg3d4_planar_mhd", ["CSWC0623"]),
    ("test_bg3d4_spherical_piston", ["CSWC0624"]),
    ("test_bg3d4_ambient_balance", ["CSWC0625"]),
    ("test_bg3d4_append", ["CSWC0626"]),
    ("test_bg3d4_piston_sheath", ["CSWC0627"]),
    ("test_bg3d4_piston_refinement", ["CSWC0628-MASS-TIME"]),
]


def main() -> int:
    context = load_context("Import bounded full-volume Corona-SWCME checks")
    try:
        repo = Path(context["repo_root"])
        work = Path(context["work_dir"])
        artifacts = Path(context["artifact_dir"])
        isolated = prepare_model_repo(repo, work)
        model = isolated / "src/models/sep_corona_swcme"
        print("[CASE] RUN CSWC full-volume bounded software checks", flush=True)
        targets = ["build/" + binary for binary, _ in BINARY_CASES]
        completed = run_command(
            list(context["tools"]["make"]) + ["-j{}".format(
                context["execution"]["build_jobs"])] + targets,
            model, timeout=3600)
        transcript = artifacts / "full-volume-local.txt"
        transcript_text = completed.stdout
        if completed.returncode != 0:
            raise RuntimeError("full-volume exact-registry build failed")
        imported_by_local = {}
        for binary, owned_ids in BINARY_CASES:
            # The refinement study is intentionally much larger than the
            # other bounded binaries and can exceed ten minutes on this host.
            # A routine import gives it a documented three-minute ceiling and
            # publishes ERROR for CSWC0628 rather than hanging a nightly or
            # pretending the partially printed mesh levels passed.  A future
            # machine profile may raise the suite timeout only after supplying
            # a qualified performance budget.
            case_timeout = 180 if binary == "test_bg3d4_piston_refinement" else 600
            try:
                result = run_command(
                    [str(model / "build" / binary)], model, timeout=case_timeout)
            except RunnerError as error:
                transcript_text += "\n$ build/{}\nERROR: {}\n".format(binary, error)
                for local_id in owned_ids:
                    imported_by_local[local_id] = case(
                        global_id("CSWC-FULL", local_id), "ERROR",
                        "local_id={} executable={} exceeded routine {} s ceiling; "
                        "no PASS is inferred from partial refinement levels; "
                        "BG3D-4-unqualified/deferred".format(
                            local_id, binary, case_timeout))
                continue
            transcript_text += "\n$ build/{}\n{}".format(binary, result.stdout)
            markers = re.findall(
                r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL|SKIP|ERROR)\s*(.*)$",
                result.stdout, flags=re.MULTILINE)
            if result.returncode == 0:
                if ([row[0] for row in markers] != owned_ids or
                        any(row[1] != "PASS" for row in markers)):
                    raise RuntimeError("{} returned zero without exact PASS markers".format(
                        binary))
                for local_id, status, message in markers:
                    imported_by_local[local_id] = case(
                        global_id("CSWC-FULL", local_id), status,
                        "local_id={} BG3D-4-unqualified/deferred; {}".format(
                            local_id, message))
            elif markers:
                # Preserve a maintained executable's typed nonzero result.
                # Collapsing ERROR into FAIL would misclassify a crash/setup
                # problem as a rejected physical assertion, while collapsing
                # SKIP would invent execution.  Exact ownership and the
                # conventional 0/1/2 exit mapping are checked before import.
                marker_ids = [row[0] for row in markers]
                if marker_ids != owned_ids:
                    raise RuntimeError(
                        "{} nonzero result omitted/reordered owned markers".format(
                            binary))
                exit_by_status = {"PASS": 0, "SKIP": 0, "FAIL": 1, "ERROR": 2}
                expected_exit = max(exit_by_status[row[1]] for row in markers)
                if result.returncode != expected_exit:
                    raise RuntimeError(
                        "{} typed marker/exit mismatch: exit {} expected {}".format(
                            binary, result.returncode, expected_exit))
                for local_id, status, message in markers:
                    imported_by_local[local_id] = case(
                        global_id("CSWC-FULL", local_id), status,
                        "local_id={} executable={} exit={}; "
                        "BG3D-4-unqualified/deferred; {}".format(
                            local_id, binary, result.returncode, message))
            else:
                # Several legacy binaries print "BG3D-4 FAIL: ..." rather
                # than a bracketed case marker.  Accept that narrow historical
                # form only for exit 1 plus an explicit FAIL diagnostic.  A
                # signal, exit 2, or untyped nonzero result is infrastructure
                # ERROR: it must not be relabeled as a scientific rejection.
                detail = " ".join(result.stdout.split())[-500:]
                legacy_fail = (
                    result.returncode == 1 and
                    re.search(r"(?:BG3D-4|CSWC[0-9-]*)\s+FAIL\s*:",
                              result.stdout) is not None)
                for local_id in owned_ids:
                    imported_by_local[local_id] = case(
                        global_id("CSWC-FULL", local_id),
                        "FAIL" if legacy_fail else "ERROR",
                        "local_id={} executable={} exit={}; BG3D-4-unqualified/deferred; {}".format(
                            local_id, binary, result.returncode, detail))

        architecture = run_command(
            list(context["tools"]["python"]) + ["test/check_architecture.py"],
            model, timeout=300)
        transcript_text += "\n$ python test/check_architecture.py\n" + architecture.stdout
        arch_markers = re.findall(
            r"^\[(ARCHCSWC01)\]\s+(PASS|FAIL|SKIP|ERROR)\s*(.*)$",
            architecture.stdout, flags=re.MULTILINE)
        if len(arch_markers) != 1:
            raise RuntimeError("architecture check omitted its unique marker")
        local_id, arch_status, message = arch_markers[0]
        expected_arch_exit = {"PASS": 0, "SKIP": 0, "FAIL": 1, "ERROR": 2}[arch_status]
        if architecture.returncode != expected_arch_exit:
            raise RuntimeError("architecture exit/status mismatch")
        imported_by_local[local_id] = case(
            global_id("CSWC-FULL", local_id), arch_status,
            "local_id={} BG3D-4-unqualified/deferred; {}".format(local_id, message))
        transcript.write_text(transcript_text, encoding="utf-8")
        if list(imported_by_local) != LOCAL_IDS:
            raise RuntimeError("full-volume exact registry order/membership drift")
        imported = [imported_by_local[local_id] for local_id in LOCAL_IDS]
        imported = selected_rows(context, imported)
        return write_result(
            context, suite_status(imported),
            "bounded code checks imported; BG3D-4 and G1 remain unqualified and deferred",
            imported, evidence=[str(transcript)],
            metrics={"local_registry_cases": 15,
                     "qualification": "bounded-software-checks-only"})
    except BaseException as error:
        return error_result(context, error)


if __name__ == "__main__":
    raise SystemExit(main())


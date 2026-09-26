#!/usr/bin/env python3
"""Production-wiring and unchanged-gate audit for Roadmap Step 10."""

from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
EARTH = HERE.parents[1]


def require(text: str, needle: str, label: str) -> None:
    if needle not in text:
        raise AssertionError(f"missing {label}: {needle}")


def forbid(text: str, needle: str, label: str) -> None:
    if needle in text:
        raise AssertionError(f"obsolete/unsafe {label} remains: {needle}")


def ordered(text: str, first: str, second: str, label: str) -> None:
    if text.find(first) < 0 or text.find(second) < 0 or text.index(first) >= text.index(second):
        raise AssertionError(f"invalid order for {label}: {first!r} before {second!r}")


def main() -> int:
    main_lib = (EARTH / "main_lib.cpp").read_text(encoding="utf-8")
    swmf = (EARTH / "3d_forward_swmf" / "Mode3DForwardSWMF.cpp").read_text(
        encoding="utf-8"
    )
    swmf_header = (EARTH / "3d_forward_swmf" / "Mode3DForwardSWMF.h").read_text(
        encoding="utf-8"
    )
    cutoff = (EARTH / "3d" / "CutoffRigidityMode3D.cpp").read_text(encoding="utf-8")
    cutoff_header = (EARTH / "3d" / "CutoffRigidityMode3D.h").read_text(
        encoding="utf-8"
    )
    parser = (EARTH / "util" / "amps_param_parser.cpp").read_text(encoding="utf-8")
    test_list = (EARTH / "test" / "list").read_text(encoding="utf-8")
    comparator = (HERE / "compare_cutoff_access.py").read_text(encoding="utf-8")

    # Cadence must be a bridge-owned collective decision.  The old rank-local statics
    # could diverge before the first solver MPI call and deadlock the component.
    require(main_lib, "ShouldRunBackwardProductCalculation(true)",
            "collective cadence decision")
    require(main_lib, "MarkBackwardProductCalculationComplete()",
            "post-success cadence commit")
    forbid(main_lib, "IsFirstCoupledCutoffCalculation",
           "rank-local first-callback cadence state")
    forbid(main_lib, "LastCoupledCutoffCalculationTime_s",
           "rank-local last-completed cadence state")
    require(swmf_header, "ShouldRunBackwardProductCalculation",
            "public collective scheduling contract")
    require(swmf, "minAction!=maxAction", "cross-rank cadence-action gate")
    require(swmf, "RejectStaleTime", "stale-time fail-closed gate")
    ordered(main_lib, "ShouldRunBackwardProductCalculation(true)",
            "amps_cutoff_time_step()", "cadence before product collectives")
    ordered(main_lib, "amps_cutoff_time_step()",
            "MarkBackwardProductCalculationComplete()", "commit only after success")

    # Product identity and completeness are derived from the frozen snapshot, not from
    # callback order.  The root artifact check is broadcast before later collectives.
    require(swmf, "BuildProductSuffix", "time/snapshot output identity")
    require(swmf, "VerifyAndWriteCoupledAccessManifest_",
            "cutoff/access artifact transaction")
    require(swmf, "GetLastCutoffArtifactFiles", "root artifact inventory")
    require(swmf, "BuildAccessManifestJson", "machine-readable manifest")
    require(swmf, "RemoveSupersededAttemptStatus_",
            "provisional-to-content-identified status transaction")
    step_body = swmf[swmf.index("void amps_cutoff_time_step()") :]
    ordered(step_body, "PrepareGlobalSWMFCoupledMagneticFieldForCutoff(true)",
            "FormatCutoffOutputSuffix_", "content ID before final product suffix")
    ordered(step_body, "statusFile=CoupledSnapshotStatusFile_(suffix)",
            "RemoveSupersededAttemptStatus_",
            "final status published before provisional cleanup")
    ordered(step_body, "RunCutoffRigidity", "VerifyAndWriteCoupledAccessManifest_",
            "artifact verification after production writer")
    require(cutoff_header, "GetLastCutoffArtifactFiles",
            "cutoff writer transaction API")
    require(cutoff, "CloseAndRecordCutoffArtifact_", "close-before-record gate")
    require(cutoff, "AUXDATA SNAPSHOT_ID", "snapshot provenance in cutoff artifacts")
    require(cutoff, "AUXDATA SNAPSHOT_MESH_REVISION",
            "mesh provenance in cutoff artifacts")

    # BOX retains the historical event code; SHUE is active and missing data remains
    # INVALID_FIELD.  This is a scientific gate, not just an input-parser feature.
    require(parser, 'boundary=="SHUE"', "active SHUE validation")
    require(parser, "ResolveShueParameters", "AUTO/numeric Shue resolution")
    require(cutoff, "BuildOuterBoundaryPolicy3D_", "production boundary dispatch")
    require(cutoff, "SegmentDisposition::MeshUnavailable",
            "mesh-limit versus physical-escape distinction")
    require(cutoff, "TrajectoryTermination::InvalidField",
            "unavailable field termination")
    require(cutoff, "Preserve the historical BOX classifier exactly",
            "unchanged default BOX path")

    # The live/replay comparator is strict by default and identity-aware.  It cannot
    # pass two unrelated states merely because their numeric rows are close.
    require(comparator, "default=0.0", "exact default comparison gate")
    require(comparator, "REQUIRED_PROVENANCE", "required snapshot provenance")
    require(comparator, 'print("RESULT: PASS")', "explicit passing result line")
    require(comparator, 'print("RESULT: FAIL"', "explicit failing result line")

    # Representative existing scientific gates remain byte-for-byte present.  Step 10
    # adds an independent entry; it does not replace, disable, or relax these tests.
    require(test_list, "P srcEarth/test/F4/run_F4.py -np 4 -nt 16",
            "unchanged F4 production gate")
    require(test_list, "--min-access-state-agreement 0.999 --max-access-unresolved-fraction 0.01",
            "unchanged C9 access/unresolved gate")
    require(test_list, "--unresolved-extension-passes 2 --unresolved-extension-factor 2.0",
            "unchanged C19 extension gate")
    require(test_list, "P srcEarth/test/USWMFCoupledAccess/run_test.sh",
            "independent Step-10 test registration")

    print("PASS S10-SOURCE coupled-access wiring and unchanged validation gates")
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        print("RESULT: FAIL", file=sys.stderr)
        raise SystemExit(1)

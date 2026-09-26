#!/usr/bin/env python3
"""Production-wiring and unchanged-gate checks for Roadmap Step 9.

The compiled C++ test validates the portable numerical contract.  These checks ensure
that the live and replay production paths actually consume that contract and that Step
9 did not weaken representative pre-existing validation commands.
"""

from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
EARTH = HERE.parents[1]


def require(text: str, needle: str, label: str) -> None:
    if needle not in text:
        raise AssertionError(f"missing {label}: {needle}")


def ordered(text: str, first: str, second: str, label: str) -> None:
    if text.find(first) < 0 or text.find(second) < 0 or text.index(first) >= text.index(second):
        raise AssertionError(f"invalid order for {label}: {first!r} before {second!r}")


def main() -> int:
    global_field = (EARTH / "3d" / "GlobalMagneticField.cpp").read_text(encoding="utf-8")
    mode3d = (EARTH / "3d" / "Mode3D.cpp").read_text(encoding="utf-8")
    electric_field = (EARTH / "3d" / "ElectricField.cpp").read_text(encoding="utf-8")
    swmf = (EARTH / "3d_forward_swmf" / "Mode3DForwardSWMF.cpp").read_text(encoding="utf-8")
    parser = (EARTH / "util" / "amps_param_parser.cpp").read_text(encoding="utf-8")
    cli = (EARTH / "util" / "cutoff_cli.cpp").read_text(encoding="utf-8")
    contract = (EARTH / "util" / "StandaloneProductContract.h").read_text(encoding="utf-8")
    replay_template = (
        EARTH / "examples" / "standalone_step9_swmf_replay.in.template"
    ).read_text(encoding="utf-8")
    test_list = (EARTH / "test" / "list").read_text(encoding="utf-8")

    require(global_field, "BuildPortableSWMFSnapshot_", "content-derived live snapshot")
    require(global_field, "ExportCurrentSWMFSnapshot", "live snapshot export")
    require(global_field, "ImportSWMFSnapshot", "standalone snapshot replay")
    require(global_field, "DeriveElectricField", "shared E=-u x B convention")
    require(global_field, "minFingerprint!=maxFingerprint", "cross-rank identity gate")
    require(global_field, "ValidateOwnerCellParity_", "direct owner/compact parity gate")
    require(global_field, "ComputeMeshRevision", "mesh revision gate")
    require(global_field, "RequireCollectiveRootWriteSuccess_",
            "collective snapshot-export I/O gate")
    require(global_field, "SWMF replay content fingerprint differs across MPI ranks",
            "cross-rank replay identity gate")
    require(global_field, "allRanksReplayValidated",
            "collective replay topology/content gate")
    require(global_field, "BeginFrozenFieldBatch", "immutable batch lease")
    require(global_field, "EndFrozenFieldBatch", "batch lease release")
    require(swmf, "ValidateLiveSWMFStateCoherence_", "live epoch/offset coherence gate")
    require(swmf, "RequireCollectiveStatusWriteSuccess_",
            "collective product-status I/O gate")
    require(swmf, "ExportCurrentSWMFSnapshot", "coupled export dispatch")
    require(swmf, 'statusFile,"FAILED"', "fail-closed initial status")
    require(swmf, 'statusFile,"PASS"', "terminal passing status")
    ordered(swmf, 'statusFile,"FAILED"', "PrepareGlobalSWMFCoupledMagneticFieldForCutoff(true)",
            "FAILED status before field gather")
    ordered(swmf, "BeginFrozenFieldBatch", "ExportCurrentSWMFSnapshot",
            "field freeze before export")
    ordered(swmf, "ExportCurrentSWMFSnapshot", "RunCutoffRigidity", "export before products")
    ordered(swmf, "RunDensityAndFlux", 'statusFile,"PASS"', "all products before PASS")
    ordered(swmf, "EndFrozenFieldBatch", 'statusFile,"PASS"',
            "field release before terminal PASS")

    require(parser, 'uKey=="SWMF_SNAPSHOT_FILE"', "snapshot-file input")
    require(parser, 'uKey=="SWMF_SNAPSHOT_EXPORT"', "coupled export switch")
    require(parser, 'uKey=="SWMF_DERIVED_ELECTRIC_FIELD"',
            "explicit experimental electric-field selector")
    require(parser, 'mode=="OFF"', "magnetic-only default selector")
    require(parser, 'mode=="EXPERIMENTAL"', "experimental derived-E label")
    # The live-SWMF selector below is compiled from AMPS-generated configuration
    # macros.  Requiring pic.h, and requiring it before the preprocessor branch,
    # prevents standalone compiler invocations from seeing undefined identifiers
    # and silently treating them as zero in #if expressions.
    require(parser, '#include "pic.h"', "AMPS coupler configuration include")
    ordered(parser, '#include "pic.h"',
            '#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_',
            "coupler definitions before SWMF preprocessor selection")
    require(mode3d, 'model=="SWMF_SNAPSHOT"', "Mode3D replay selector")
    require(mode3d, "ImportSWMFSnapshot", "Mode3D replay import")
    require(mode3d, "SWMF_SNAPSHOT_FILE domain differs", "replay domain coherence")
    require(mode3d, "CurrentMeshRevision", "replay mesh-revision publication")
    require(electric_field, "SWMF_SNAPSHOT cannot be evaluated analytically",
            "fail-closed non-mesh replay evaluator")
    require(contract, 'fieldModel=="SWMF_SNAPSHOT" && representation!=FieldRepresentation::Mesh',
            "gridless replay rejection")
    require(replay_template, "SWMF_DERIVED_ELECTRIC_FIELD OFF",
            "released magnetic-only replay mode")
    require(cli, "SWMF_DERIVED_ELECTRIC_FIELD OFF | EXPERIMENTAL",
            "Step-9 command help")

    # Import must invalidate the previous generation before opening/parsing the new
    # artifact. Otherwise a caught corrupt-file exception could fall back to old B/u.
    import_body = global_field[global_field.index("MaterializationStats ImportSWMFSnapshot"):]
    ordered(import_body, "GlobalFieldsReady_=false", "SWMFSnapshot::Read(fileName)",
            "fail-closed replay replacement")

    # The documented replay deck must use real parser keywords.  This catches a subtle
    # but costly integration regression where an otherwise correct implementation is
    # shipped with a template that the strict parser rejects before field import.
    template_parser_keys = (
        "MODE3D_MESH_RESOLUTION_EARTH_KM",
        "MODE3D_MESH_RESOLUTION_BOUNDARY_KM",
        "MODE3D_MESH_OUTER_RADIUS_KM",
        "MODE3D_MESH_COARSENING",
        "MODE3D_MESH_EXPONENT",
    )
    for key in template_parser_keys:
        require(replay_template, key, f"replay template key {key}")
        require(parser, f'uKey=="{key}"', f"parser support for replay key {key}")

    # Existing scientific gates remain exact.  These strings intentionally include the
    # strict unresolved/access settings most vulnerable to accidental relaxation.
    require(test_list, "P srcEarth/test/F4/run_F4.py -np 4 -nt 16",
            "unchanged F4 production gate")
    require(test_list, "--min-access-state-agreement 0.999 --max-access-unresolved-fraction 0.01",
            "unchanged C9 access/unresolved gate")
    require(test_list, "--unresolved-extension-passes 2 --unresolved-extension-factor 2.0",
            "unchanged C19 extension gate")
    require(test_list, "P srcEarth/test/USWMFSnapshot/run_test.sh",
            "independent Step-9 test registration")

    print("PASS: Step-9 live/replay wiring and unchanged validation gates")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)

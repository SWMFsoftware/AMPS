#!/usr/bin/env python3
"""Structural wiring checks for the Step-7 production paths.

These checks do not replace linked C/F numerical validation.  They make sure the
dependency-free contract tested by test_standalone_products.cpp is actually consumed by
both production dispatchers, and that T01/TA15 empirical calls remain excluded from an
SWMF build by the existing compile-time guard.
"""

from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
EARTH = HERE.parents[1]


def require(text: str, needle: str, label: str) -> None:
    if needle not in text:
        raise AssertionError(f"missing {label}: {needle}")


def main() -> int:
    main_cpp = (EARTH / "main.cpp").read_text(encoding="utf-8")
    mode3d = (EARTH / "3d" / "Mode3D.cpp").read_text(encoding="utf-8")
    electric = (EARTH / "3d" / "ElectricField.cpp").read_text(encoding="utf-8")
    parser = (EARTH / "util" / "amps_param_parser.cpp").read_text(encoding="utf-8")

    for source, label in ((main_cpp, "gridless"), (mode3d, "Mode3D")):
        require(source, "namespace SP=Earth::StandaloneProducts;", f"{label} shared contract")
        require(source, "BuildManifestJson", f"{label} manifest")
        require(source, "RequireSameSnapshot", f"{label} immutable snapshot gate")

    require(main_cpp, "RunCutoffRigidity(p)", "gridless cutoff dispatch")
    require(main_cpp, "RunDensityAndSpectrum(p)", "gridless flux dispatch")
    require(mode3d, 'model=="T01"', "Mode3D T01 initialization")
    require(mode3d, 'model=="TA15N" || model=="TA15B"', "Mode3D TA15 initialization")
    require(electric, "::T01::GetMagneticField", "Mode3D T01 evaluation")
    require(electric, "::TA15::GetMagneticField", "Mode3D TA15 evaluation")
    require(parser, "ValidateTsDriverUnits", "driver unit gate")
    require(parser, "ParseFiniteDriverValue", "driver row-value gate")
    require(parser, "SetValidationProvenance", "driver validation provenance")

    # Empirical includes and calls must remain inside the non-SWMF branch.  Checking
    # ordering is intentionally strict enough to catch accidental movement above the
    # guard without encoding formatting or line numbers.
    guard = mode3d.index("#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_")
    endif = mode3d.index("#endif", guard)
    for include in ('#include "T01Interface.h"', '#include "TA15Interface.h"'):
        pos = mode3d.index(include)
        if not guard < pos < endif:
            raise AssertionError(f"{include} escaped the non-SWMF include guard")

    # The include guard alone is not enough: calls must also remain on the standalone
    # side of the function-level branch.  Check the complete configuration/evaluation
    # regions so a future edit cannot add a T01/TA15 call to the SWMF branch while
    # leaving the includes correctly guarded.
    configure = mode3d[
        mode3d.index("void ConfigureBackgroundFieldModel"):
        mode3d.index("// Traverse the full AMR tree")
    ]
    configure_else = configure.index("#else")
    configure_end = configure.rindex("#endif")
    for call in ("::T01::Init", "::TA15::Init"):
        pos = configure.index(call)
        if not configure_else < pos < configure_end:
            raise AssertionError(f"{call} escaped standalone Mode3D configuration")

    evaluate = electric[
        electric.index("void EvaluateBackgroundMagneticFieldSI"):
        electric.index("static void EvalCorotationSI")
    ]
    evaluate_else = evaluate.index("#else")
    evaluate_end = evaluate.rindex("#endif")
    for call in ("::T01::GetMagneticField", "::TA15::GetMagneticField"):
        pos = evaluate.index(call)
        if not evaluate_else < pos < evaluate_end:
            raise AssertionError(f"{call} escaped standalone field evaluation")

    print("PASS: Step-7 production wiring and SWMF isolation")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)

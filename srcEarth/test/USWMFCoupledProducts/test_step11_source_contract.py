#!/usr/bin/env python3
"""Production-wiring audit for Roadmap Step 11; no threshold is rewritten here."""

from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
EARTH = HERE.parents[1]


def require(text: str, needle: str, label: str) -> None:
    if needle not in text:
        raise AssertionError(f"missing {label}: {needle}")


def ordered(text: str, first: str, second: str, label: str) -> None:
    first_at, second_at = text.find(first), text.find(second)
    if first_at < 0 or second_at < 0 or first_at >= second_at:
        raise AssertionError(f"invalid order for {label}: {first!r} before {second!r}")


def main() -> int:
    bridge = (EARTH / "3d_forward_swmf" / "Mode3DForwardSWMF.cpp").read_text(
        encoding="utf-8"
    )
    density = (EARTH / "3d" / "DensityMode3D.cpp").read_text(encoding="utf-8")
    density_header = (EARTH / "3d" / "DensityMode3D.h").read_text(encoding="utf-8")
    contract = (EARTH / "util" / "SWMFCoupledProductsContract.h").read_text(
        encoding="utf-8"
    )
    comparator = (HERE / "compare_flux_products.py").read_text(encoding="utf-8")
    test_list = (EARTH / "test" / "list").read_text(encoding="utf-8")
    callback = bridge[bridge.index("void amps_cutoff_time_step()") :]

    # The coupled callback must invoke the shared backward-characteristic integrator,
    # not the forward residence-time density/sphere samplers.
    require(bridge, "Earth::Mode3D::RunDensityAndFlux(batchProductPrm)",
            "direct Step-6 integrator call")
    require(density, "EvaluateIsotropicProducts(",
            "shared product kernel")
    require(density, "This is NOT the same as 3d_forward/Density3D.cpp",
            "forward/backward distinction")
    ordered(callback, "RunDensityAndFlux(batchProductPrm)",
            "VerifyAndWriteCoupledProductsManifest_",
            "physics completes before coupled product manifest")
    ordered(callback, "VerifyAndWriteCoupledProductsManifest_",
            'statusFile,"PASS"', "manifest completes before PASS status")

    # Synchronization covers more than the field: controls include spectrum, channels,
    # response functions, and observation/ephemeris state and are compared collectively.
    require(bridge, "DescribeDensityFluxProductControl(s_prm)",
            "collective product-control description")
    require(bridge, "BuildSynchronizedCoupledProductParam_",
            "field-epoch ephemeris selection")
    require(bridge, "CoupledEpochsMatch_", "one-millisecond epoch matching")
    require(bridge, "RequireCollectiveCoupledProductControl_(batchProductPrm)",
            "post-selection collective control check")
    require(contract, "ObservationStateFingerprint", "spacecraft-state identity")
    require(contract, "BoundarySpectrumFingerprint", "selected spectrum-table identity")
    require(contract, "DetectorResponseFingerprint", "response identity")
    require(contract, "ChannelSchemaFingerprint", "channel identity")
    require(contract, "spectrumEvaluationEpochUTC!=epochUTC",
            "field/spectrum epoch equality gate")

    # Product publication is a close-verified transaction with explicit termination
    # closure and the retained unresolved threshold.
    require(density, "CloseAndRecordDensityArtifact_", "close-verified inventory")
    require(density, "out.close()", "explicit artifact close")
    require(density, "MPI_Bcast(&outputSucceeded", "collective root-output failure")
    require(density_header, "GetLastDensityFluxRunSummary",
            "coupled transaction summary API")
    require(contract, "terminationTotal!=summary.sampled",
            "termination-count closure")
    require(contract, "maximumUnresolvedFraction>summary.unresolvedTolerance",
            "unresolved release gate")
    require(bridge, "saveTerminationSummary=true",
            "mandatory coupled termination evidence")
    require(density, "AUXDATA BOUNDARY_SPECTRUM_EVALUATION_EPOCH_UTC",
            "boundary-spectrum epoch provenance")
    require(density, "AUXDATA PRODUCT_CONTROL_FINGERPRINT",
            "physics-control provenance")
    require(density, "AUXDATA BOUNDARY_SPECTRUM_FINGERPRINT",
            "selected boundary-distribution provenance")
    require(density, "AUXDATA RESPONSE_WEIGHTED_UNRESOLVED_UPPER_BOUND",
            "conservative response-support unresolved bound")

    # The complete-set comparator must remain exact by default and enforce all product
    # metadata.  A user can pass a predeclared tolerance, but there is no automatic
    # fallback or tolerance growth.
    require(comparator, "default=0.0", "exact default parity gate")
    require(comparator, "REQUIRED_PRODUCT_AUX", "required Step-11 provenance")
    require(comparator, 'print("RESULT: PASS")', "explicit comparator PASS")
    require(comparator, 'print(f"RESULT: FAIL:', "explicit comparator FAIL")

    # Representative retained validation gates are asserted verbatim.  The Step-11
    # suite is additive and cannot replace, disable, or weaken these entries.
    require(test_list, "P srcEarth/test/F4/run_F4.py -np 4 -nt 16",
            "unchanged F4 gate")
    require(test_list,
            "--min-access-state-agreement 0.999 --max-access-unresolved-fraction 0.01",
            "unchanged C9/C10 access gate")
    require(test_list,
            "--unresolved-extension-passes 2 --unresolved-extension-factor 2.0",
            "unchanged C19 unresolved-extension gate")
    require(test_list, "P srcEarth/test/USWMFCoupledProducts/run_test.sh",
            "independent Step-11 registration")

    print("PASS S11-SOURCE coupled product wiring and unchanged validation gates")
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"RESULT: FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)

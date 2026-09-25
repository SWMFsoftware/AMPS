#!/usr/bin/env python3
"""Strict schema tests for the observation-facing Step-5/6 C19 contract."""

from __future__ import annotations

import importlib.util
import re
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
C19_RUNNER = HERE.parent / "C19" / "run_C19.py"


def load_c19_module():
    spec = importlib.util.spec_from_file_location("c19_step5_reader_test", C19_RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load %s" % C19_RUNNER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


VARIABLES = [
    "lon_deg", "lat_deg", "direction_weight_sr", "rigidity_GV", "energy_MeV",
    "access_state", "allowed", "weighted_access_sr", "unresolved",
    "termination_code", "exit_state_valid", "x_exit_m", "y_exit_m", "z_exit_m",
    "px_exit_SI", "py_exit_SI", "pz_exit_SI", "vx_exit_unit", "vy_exit_unit",
    "vz_exit_unit", "cos_alpha_exit", "trace_time_at_exit_s",
    "rigidity_at_exit_GV", "adaptive_refined_intervals",
    "adaptive_estimated_error_GV", "adaptive_max_ambiguous_width_GV",
    "adaptive_target_reached", "adaptive_max_samples_reached",
    "response_weighted_unresolved_support",
]

# DIRECT_ACCESS is consumed by validation tools older than Step 5.  Those readers
# intentionally support extra trailing columns, but the first 24 values are a stable
# positional ABI.  Keeping this inventory in the test makes an accidental insertion
# fail before an expensive C8/C19 production run is attempted.
LEGACY_DIRECT_ACCESS_PREFIX = [
    "lon_deg", "lat_deg", "rigidity_GV", "energy_MeV", "access_state",
    "allowed", "unresolved", "termination_code", "trace_time_s",
    "trace_distance_Re", "trace_steps", "retry_count",
    "primary_termination_code", "primary_trace_time_s",
    "trace_extension_count", "initial_trace_limit_s", "final_trace_limit_s",
    "mirror_points", "bounce_cycles", "drift_revolutions", "drift_angle_deg",
    "drift_mean_radius_change_Re", "trap_mechanism", "momentum_relative_spread",
]

STEP5_DIRECT_ACCESS_SUFFIX = [
    "direction_weight_sr", "weighted_access_sr", "exit_state_valid", "x_exit_m",
    "y_exit_m", "z_exit_m", "px_exit_SI", "py_exit_SI", "pz_exit_SI",
    "vx_exit_unit", "vy_exit_unit", "vz_exit_unit", "cos_alpha_exit",
    "trace_time_at_exit_s", "rigidity_at_exit_GV", "adaptive_refined_intervals",
    "adaptive_estimated_error_GV", "adaptive_max_ambiguous_width_GV",
    "adaptive_target_reached", "adaptive_max_samples_reached",
    "response_weighted_unresolved_support",
]

# Step 6 is also append-only.  These fields turn each access row into a physical
# directional-intensity sample while retaining lower/upper uncertainty bounds for
# unresolved characteristics.  Keeping this as a separate exact inventory makes the
# compatibility rule explicit: neither the legacy nor Step-5 blocks may move, and the
# Step-6 block may not silently lose, rename, or reorder a physical quantity.
STEP6_DIRECT_ACCESS_SUFFIX = [
    "J_boundary_perMeV", "J_boundary_lower_perMeV", "J_boundary_upper_perMeV",
    "boundary_factor", "boundary_factor_lower", "boundary_factor_upper",
    "J_directional_local_perMeV", "J_directional_local_lower_perMeV",
    "J_directional_local_upper_perMeV",
]


def producer_direct_access_columns(source: str):
    """Extract the Step-5 DIRECT_ACCESS VARIABLES declaration from C++ source."""
    anchor = source.index('\\"direction_weight_sr\\"')
    start = source.rfind('"VARIABLES=', 0, anchor)
    end = source.find('\\n");', anchor)
    if start < 0 or end < 0:
        raise AssertionError("cannot locate DIRECT_ACCESS VARIABLES declaration")
    return re.findall(r'\\"([^"\\]+)\\"', source[start:end])


def record(*, allowed: bool, overrides=None):
    values = {
        "lon_deg": 0.0,
        "lat_deg": 0.0,
        "direction_weight_sr": 0.25,
        "rigidity_GV": 2.0 if allowed else 1.0,
        "energy_MeV": 200.0 if allowed else 100.0,
        "access_state": 1 if allowed else 0,
        "allowed": 1 if allowed else 0,
        "weighted_access_sr": 0.25 if allowed else 0.0,
        "unresolved": 0,
        "termination_code": 0 if allowed else 1,
        "exit_state_valid": 1 if allowed else 0,
        "x_exit_m": 12.0 if allowed else 0.0,
        "y_exit_m": -3.0 if allowed else 0.0,
        "z_exit_m": 8.0 if allowed else 0.0,
        "px_exit_SI": 2.0e-19 if allowed else 0.0,
        "py_exit_SI": 0.0,
        "pz_exit_SI": 0.0,
        "vx_exit_unit": 1.0 if allowed else 0.0,
        "vy_exit_unit": 0.0,
        "vz_exit_unit": 0.0,
        "cos_alpha_exit": 0.25 if allowed else 0.0,
        "trace_time_at_exit_s": 2.5 if allowed else 0.0,
        "rigidity_at_exit_GV": 2.0 if allowed else 0.0,
        "adaptive_refined_intervals": 5,
        "adaptive_estimated_error_GV": 0.01,
        "adaptive_max_ambiguous_width_GV": 0.01,
        "adaptive_target_reached": 1,
        "adaptive_max_samples_reached": 0,
        "response_weighted_unresolved_support": 0.0,
        "J_boundary_perMeV": 2.0,
        "J_boundary_lower_perMeV": 1.8,
        "J_boundary_upper_perMeV": 2.2,
        "boundary_factor": 1.0 if allowed else 0.0,
        "boundary_factor_lower": 1.0 if allowed else 0.0,
        "boundary_factor_upper": 1.0 if allowed else 0.0,
        "J_directional_local_perMeV": 2.0 if allowed else 0.0,
        "J_directional_local_lower_perMeV": 1.8 if allowed else 0.0,
        "J_directional_local_upper_perMeV": 2.2 if allowed else 0.0,
    }
    if overrides:
        values.update(overrides)
    return values


def write_cube(path: Path, rows, variables=VARIABLES):
    header = "VARIABLES=" + ",".join('"%s"' % name for name in variables)
    zone = 'ZONE T="loc=0 x_km=1 y_km=2 z_km=3 frame=SM" I=%d F=POINT' % len(rows)
    text_rows = [" ".join(str(row[name]) for name in variables) for row in rows]
    path.write_text("\n".join([header, zone] + text_rows) + "\n")


def expect_value_error(module, path: Path, message: str):
    try:
        module.parse_directional_access(path)
    except ValueError:
        return
    raise AssertionError(message)


def main() -> int:
    c19 = load_c19_module()

    # Source-integration gate: both producers must expose the same required public
    # columns and call the same detailed adaptive/error-control path.  Numerical tests
    # above validate the kernels; this check prevents one backend from silently keeping
    # the pre-Step-5 writer or omitting the captured exit state.
    src_earth = HERE.parent.parent
    for relative in ("gridless/CutoffRigidityGridless.cpp",
                     "3d/CutoffRigidityMode3D.cpp"):
        source_path = src_earth / relative
        source = source_path.read_text(errors="replace")
        step5_columns = LEGACY_DIRECT_ACCESS_PREFIX + STEP5_DIRECT_ACCESS_SUFFIX
        actual_columns = producer_direct_access_columns(source)
        if actual_columns[:len(step5_columns)] != step5_columns:
            raise AssertionError(
                "%s changed the legacy/Step-5 DIRECT_ACCESS prefix:\nexpected %s\nactual   %s" %
                (relative, step5_columns, actual_columns[:len(step5_columns)]))
        if actual_columns[len(step5_columns):] != STEP6_DIRECT_ACCESS_SUFFIX:
            raise AssertionError(
                "%s changed the exact Step-6 DIRECT_ACCESS suffix:\nexpected %s\nactual   %s" %
                (relative, STEP6_DIRECT_ACCESS_SUFFIX,
                 actual_columns[len(step5_columns):]))
        for name in VARIABLES:
            if '\\"%s\\"' % name not in source:
                raise AssertionError("%s writer is missing Step-5 column %s" %
                                     (relative, name))
        for marker in (
                "EvaluateAdaptiveDirectAccessDirectionDetailed",
                "controls.absoluteTolerance_GV",
                "controls.relativeTolerance",
                "controls.maximumSamples",
                "RegularLonLatCellWeightSr",
                "out.exitState=tr.exitState"):
            if marker not in source:
                raise AssertionError("%s is missing Step-5 integration marker %s" %
                                     (relative, marker))

        # Header order alone is insufficient: a reordered fprintf argument list would
        # produce a syntactically valid but physically corrupt file.  Verify the two
        # deliberately separated value groups mirror the legacy prefix and Step-5
        # suffix exactly in both backends.
        compact = re.sub(r"\s+", "", source)
        legacy_values = (
            "lon_deg,lat_deg,rigidity,energy,state,allowed,unresolved,"
            "d.terminationCode,d.traceTime_s,d.traceDistance_Re,d.steps,d.retryCount,"
            "d.primaryTerminationCode,d.primaryTraceTime_s,d.traceExtensionCount,"
            "d.initialTraceLimit_s,d.finalTraceLimit_s,d.mirrorPoints,d.bounceCycles,"
            "d.driftRevolutions,d.driftAngle_deg,d.driftMeanRadiusChange_Re,"
            "d.trapMechanism,d.momentumRelativeSpread")
        step5_values = (
            "directionWeight_sr,weightedAccess_sr,d.exitStateValid,d.xExit_m[0],"
            "d.xExit_m[1],d.xExit_m[2],d.pExit_SI[0],d.pExit_SI[1],d.pExit_SI[2],"
            "d.vExitUnit[0],d.vExitUnit[1],d.vExitUnit[2],d.cosAlphaExit,"
            "d.traceTimeAtExit_s,d.rigidityAtExit_GV,d.adaptiveRefinedIntervals,"
            "d.adaptiveEstimatedError_GV,d.adaptiveMaxAmbiguousWidth_GV,"
            "d.adaptiveTargetReached,d.adaptiveMaxSamplesReached,"
            "d.responseWeightedUnresolvedSupport")
        step6_values = (
            "boundaryIntensity.nominal,boundaryIntensity.lower,boundaryIntensity.upper,"
            "boundaryFactor.nominal,boundaryFactor.lower,boundaryFactor.upper,"
            "directionalIntensity.nominal,directionalIntensity.lower,"
            "directionalIntensity.upper")
        for label, values in (("legacy prefix", legacy_values),
                              ("Step-5 suffix", step5_values),
                              ("Step-6 suffix", step6_values)):
            if values not in compact:
                raise AssertionError("%s writer has mismatched %s value order" %
                                     (relative, label))

    # Capturing an allowed characteristic's boundary phase space is necessary for a
    # saved A(R,Omega) row, but it must not alter legacy scalar/PENUMBRA/RIGIDITY_LIST
    # trajectories.  These source gates require an explicit Boolean at every shared
    # classifier call site: false for legacy products and true for persisted direct
    # access diagnostics.
    gridless_source = re.sub(
        r"\s+", "", (src_earth / "gridless/CutoffRigidityGridless.cpp").read_text())
    for marker in (
            "doubleR_GV,boolcaptureExitState)->CutoffSampleDiagnosticGridless_",
            "R_GV,-1.0,captureExitState,GetDefaultMoverType()",
            "ClassifyCutoffSampleDetailed(taskField,x0_m,v0,R_GV,false)",
            "ClassifyCutoffSampleDetailed(taskField,x0_m,v0,rigidity_GV,true)",
            "prm.cutoff.rigidityList_GV[(std::size_t)iRigidity],true)"):
        if marker not in gridless_source:
            raise AssertionError(
                "gridless exit-state capture scope is missing marker %s" % marker)

    mode3d_source = re.sub(
        r"\s+", "", (src_earth / "3d/CutoffRigidityMode3D.cpp").read_text())
    for marker in (
            "constDomainBox3D&box,boolcaptureExitState)",
            "box,-1.0,captureExitState,GetDefaultMoverType()",
            "R_GV,q_C,m0_kg,box,false)",
            "rigidity_GV,q_C,m0,box,true)",
            "prm.cutoff.rigidityList_GV[(std::size_t)iRigidity],q_C,m0,box,true)"):
        if marker not in mode3d_source:
            raise AssertionError(
                "Mode3D exit-state capture scope is missing marker %s" % marker)
    if mode3d_source.count("R_GV,q_C,m0_kg,box,false)") < 2:
        raise AssertionError(
            "Mode3D scalar and PENUMBRA classifiers must both disable exit capture")

    with tempfile.TemporaryDirectory(prefix="c19_step5_reader_") as tmp:
        root = Path(tmp)

        valid = root / "valid.dat"
        write_cube(valid, [record(allowed=False), record(allowed=True)])
        cube = c19.parse_directional_access(valid)
        curve = cube.samples[(0.0, 0.0)]
        assert len(curve) == 2
        assert curve[1].direction_weight_sr == 0.25
        assert curve[1].weighted_access_sr == 0.25
        assert curve[1].exit_state_valid == 1
        assert curve[1].adaptive_target_reached == 1

        # The observation reader must consume the current 54-column producer output,
        # not merely tolerate its source declaration. This guards the prior C19
        # failure mode where a schema mismatch yielded no comparison at all.
        valid_step6 = root / "valid_step6.dat"
        write_cube(valid_step6,[record(allowed=False),record(allowed=True)],
                   variables=VARIABLES+STEP6_DIRECT_ACCESS_SUFFIX)
        cube_step6=c19.parse_directional_access(valid_step6)
        if cube_step6.samples != cube.samples:
            raise AssertionError(
                "Step-6 suffix changed C19 access samples or suppressed comparison input")

        bad_weight = root / "bad_weight.dat"
        write_cube(bad_weight, [
            record(allowed=False),
            record(allowed=True, overrides={"weighted_access_sr": 0.0}),
        ])
        expect_value_error(c19, bad_weight,
                           "reader accepted inconsistent weighted access")

        missing_exit = root / "missing_exit.dat"
        write_cube(missing_exit, [
            record(allowed=False),
            record(allowed=True, overrides={"exit_state_valid": 0}),
        ])
        expect_value_error(c19, missing_exit,
                           "reader accepted allowed row without exit state")

        nonunit_velocity = root / "nonunit_velocity.dat"
        write_cube(nonunit_velocity, [
            record(allowed=False),
            record(allowed=True, overrides={"vx_exit_unit": 0.5}),
        ])
        expect_value_error(c19, nonunit_velocity,
                           "reader accepted a non-unit boundary velocity")

        inconsistent_report = root / "inconsistent_report.dat"
        write_cube(inconsistent_report, [
            record(allowed=False),
            record(allowed=True, overrides={"adaptive_refined_intervals": 6}),
        ])
        expect_value_error(c19, inconsistent_report,
                           "reader accepted direction-varying convergence metadata")

        incomplete_schema = root / "incomplete_schema.dat"
        incomplete_variables = VARIABLES[:-1]
        write_cube(incomplete_schema,
                   [record(allowed=False), record(allowed=True)],
                   variables=incomplete_variables)
        expect_value_error(c19, incomplete_schema,
                           "reader accepted a partially populated Step-5 schema")

        # Archived seven-column files remain readable.  This is compatibility, not a
        # relaxation of the new schema: presence of direction_weight_sr selects the
        # strict Step-5 contract above.
        archive = root / "archive.dat"
        archive_variables = [
            "lon_deg", "lat_deg", "rigidity_GV", "energy_MeV",
            "access_state", "allowed", "unresolved",
        ]
        write_cube(archive,
                   [record(allowed=False), record(allowed=True)],
                   variables=archive_variables)
        archived_cube = c19.parse_directional_access(archive)
        assert archived_cube.samples[(0.0, 0.0)][1].direction_weight_sr is None

    print("UDirectionalAccess C19 reader: PASS (U-F16)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""C16 — trajectory termination and access-state closure validation.

C16 validates the *classification contract* at the boundary between the AMPS
trajectory integrator and the cutoff/access calculator.  It is intentionally
different from a cutoff-value comparison: every requested trajectory must carry
an explicit termination reason, and that reason must agree with the three-state
access value written to the DIRECT_ACCESS cube.

The routine profile performs two complementary checks:

1. Ordinary DIPOLE and T05 calculations are repeated with increasing physical
   trace-time budgets.  Once a sample resolves to ALLOWED or FORBIDDEN, a longer
   budget may not change that physical state or turn it back into UNRESOLVED.
2. Small DIPOLE probe runs deliberately exhaust time, step, and distance limits.
   At least 90 percent of each probe must report its intended termination code,
   and every such outcome must remain access_state=2 (UNRESOLVED).

The runner uses only Python's standard library.  ``--self-test`` exercises the
parser and every acceptance guard without an AMPS executable.
"""

from __future__ import print_function

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from collections import Counter, defaultdict
from pathlib import Path


TEST_ID = "C16"
TEST_NAME = "Trajectory termination and access-state closure"
RUNNER_SCHEMA_VERSION = 1
RUNNER_RELEASE = "2026-08-28"

RE_KM = 6371.2
DEFAULT_EPOCH = "2012-05-17T06:00:00"
DEFAULT_RIGIDITIES_GV = (0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)

# This dictionary is the executable scientific reference for C16.  The checked-
# in CSV contains the same rows in a human-readable form and is verified by the
# self-test, so documentation drift cannot silently alter the acceptance logic.
TERMINATION_CONTRACT = {
    0: ("OUTER_BOUNDARY_ALLOWED", 1, "physical", "accept"),
    1: ("INNER_BOUNDARY_FORBIDDEN", 0, "physical", "accept"),
    2: ("MAGNETICALLY_TRAPPED_FORBIDDEN", 0, "physical", "accept"),
    3: ("TIME_LIMIT", 2, "resource_limit", "accept_as_unresolved"),
    4: ("STEP_LIMIT", 2, "resource_limit", "accept_as_unresolved"),
    5: ("DISTANCE_LIMIT", 2, "resource_limit", "accept_as_unresolved"),
    6: ("INVALID_TIME_STEP", 2, "numerical_fatal", "fail_run"),
    7: ("INVALID_FIELD", 2, "numerical_fatal", "fail_run"),
    8: ("NUMERICAL_FAILURE", 2, "numerical_fatal", "fail_run"),
    9: ("DRIFT_TRAPPED_FORBIDDEN", 0, "physical", "accept"),
}

STATE_NAMES = {0: "FORBIDDEN", 1: "ALLOWED", 2: "UNRESOLVED"}
FATAL_CODES = {6, 7, 8}

# Probe values are deliberately tiny relative to an ordinary magnetospheric
# trace.  They make the requested resource limit dominate before a physical
# boundary can normally be reached, giving a deterministic integration test of
# each resource-exhaustion branch without adding test-only parser keywords.
DEFAULT_PROBES = {
    "time_limit": {
        "max_trace_time": 0.01,
        "max_steps": 2000000,
        "max_trace_distance": 0.0,
        "expected_code": 3,
    },
    "step_limit": {
        "max_trace_time": 600.0,
        "max_steps": 1,
        "max_trace_distance": 0.0,
        "expected_code": 4,
    },
    "distance_limit": {
        "max_trace_time": 600.0,
        "max_steps": 2000000,
        "max_trace_distance": 0.01,
        "expected_code": 5,
    },
}

# Profiles alter resolution and cost only; they never relax state/termination
# closure or allow fatal numerical reasons.  ROUTINE is the recommended model-
# validation profile and is the command documented in README.md.
PROFILE_DEFAULTS = {
    "SMOKE": {
        "models": "DIPOLE",
        "budgets": "60,600",
        "lats": "-60,0,60",
        "dir_lon_res": 60.0,
        "dir_lat_res": 60.0,
        "rigidities": "0.2,1,5",
        "max_steps": 1000000,
        "max_final_unresolved_fraction": 0.70,
    },
    "ROUTINE": {
        "models": "DIPOLE,T05",
        "budgets": "60,300,600",
        "lats": "-60,-30,0,30,60",
        "dir_lon_res": 30.0,
        "dir_lat_res": 30.0,
        "rigidities": ",".join("%.10g" % value for value in DEFAULT_RIGIDITIES_GV),
        "max_steps": 2000000,
        "max_final_unresolved_fraction": 0.50,
    },
    "THOROUGH": {
        "models": "DIPOLE,T05",
        "budgets": "60,300,600,2400",
        "lats": "-60,-30,0,30,60",
        "dir_lon_res": 15.0,
        "dir_lat_res": 15.0,
        "rigidities": ",".join("%.10g" % value for value in DEFAULT_RIGIDITIES_GV),
        "max_steps": 4000000,
        "max_final_unresolved_fraction": 0.25,
    },
}


def parse_float_list(text, label):
    """Parse a comma/semicolon list and reject empty or non-finite values."""
    values = []
    for token in str(text).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = float(token)
        except ValueError:
            raise SystemExit("Invalid %s value: %r" % (label, token))
        if not math.isfinite(value):
            raise SystemExit("Non-finite %s value: %r" % (label, token))
        values.append(value)
    if not values:
        raise SystemExit("No values supplied for %s" % label)
    return values


def parse_name_list(text, label, allowed):
    """Parse a case-insensitive name list while preserving request order."""
    result = []
    for token in str(text).replace(";", ",").split(","):
        name = token.strip().upper()
        if not name:
            continue
        if name not in allowed:
            raise SystemExit("Unsupported %s %r; choose from %s" %
                             (label, name, ",".join(sorted(allowed))))
        if name not in result:
            result.append(name)
    if not result:
        raise SystemExit("No values supplied for %s" % label)
    return result


def bool_token(value):
    """Normalize common boolean spellings to the AMPS T/F syntax."""
    token = str(value).strip().upper()
    if token in ("T", "TRUE", "1", "YES", "Y"):
        return "T"
    if token in ("F", "FALSE", "0", "NO", "N"):
        return "F"
    raise SystemExit("Expected a T/F boolean, got %r" % value)


def fmt_number(value):
    """Produce compact, deterministic text for directory names and inputs."""
    return "%.12g" % float(value)


def budget_label(seconds):
    """Return a filesystem-safe label without losing fractional probe limits."""
    text = fmt_number(seconds).replace("-", "m").replace(".", "p")
    return "t%s" % text


def norm_lon(lon_deg):
    """Canonicalize a sky longitude so 0 and 360 share one output key."""
    value = float(lon_deg) % 360.0
    return 0.0 if math.isclose(value, 0.0, abs_tol=1.0e-8) else value


def norm_lat(lat_deg):
    """Remove negative zero while retaining the requested polar labels."""
    value = float(lat_deg)
    return 0.0 if math.isclose(value, 0.0, abs_tol=1.0e-8) else value


def direction_key(lon_deg, lat_deg):
    """Quantize textual output coordinates at well below AMPS map resolution."""
    return (round(norm_lon(lon_deg), 7), round(norm_lat(lat_deg), 7))


def rigidity_key(rigidity_gv):
    """Quantize a rigidity only for dictionary lookup, never for comparison."""
    return round(float(rigidity_gv), 10)


def point_xyz_km(lon_deg, lat_deg, alt_km):
    """Convert a spherical GSM point definition into Cartesian kilometers."""
    radius = RE_KM + float(alt_km)
    lon = math.radians(float(lon_deg))
    lat = math.radians(float(lat_deg))
    cos_lat = math.cos(lat)
    return (
        radius * cos_lat * math.cos(lon),
        radius * cos_lat * math.sin(lon),
        radius * math.sin(lat),
    )


def build_points(lons_deg, lats_deg, alt_km):
    """Build the fixed observation-point table used by every baseline budget."""
    points = []
    for lon_deg in lons_deg:
        for lat_deg in lats_deg:
            xyz = point_xyz_km(lon_deg, lat_deg, alt_km)
            points.append({
                "point_id": len(points),
                "obs_lon_deg": float(lon_deg),
                "obs_lat_deg": float(lat_deg),
                "obs_alt_km": float(alt_km),
                "x_km": xyz[0],
                "y_km": xyz[1],
                "z_km": xyz[2],
            })
    return points


def points_block(points):
    """Render the parser-supported POINTS block with audit-friendly comments."""
    lines = []
    for point in points:
        lines.append(
            "POINT %.12e %.12e %.12e ! km point_id=%d lon=%g lat=%g alt=%g" %
            (point["x_km"], point["y_km"], point["z_km"],
             point["point_id"], point["obs_lon_deg"], point["obs_lat_deg"],
             point["obs_alt_km"]))
    return "\n".join(lines)


def expected_direction_keys(lon_res_deg, lat_res_deg):
    """Derive the complete regular sky grid independently of AMPS output.

    Unlike C17's reflection comparison, C16 retains polar rows.  At a pole the
    repeated longitude labels correspond to separate requested work items and
    must not disappear from the termination-accounting closure.
    """
    nlon = int(round(360.0 / float(lon_res_deg)))
    nlat_intervals = int(round(180.0 / float(lat_res_deg)))
    return {
        direction_key(i * float(lon_res_deg),
                      -90.0 + j * float(lat_res_deg))
        for j in range(nlat_intervals + 1)
        for i in range(nlon)
    }


def field_block(model):
    """Render a complete field-model block without leaving stale directives."""
    if model == "DIPOLE":
        return "\n".join([
            "FIELD_MODEL                    DIPOLE",
            "DIPOLE_MOMENT                  1.0",
            "DIPOLE_TILT                    0.0",
        ])
    if model == "T05":
        # Static values are parser-safe fallbacks.  The driver row at EPOCH is
        # authoritative and exercises the same T05 input path as C9/C19.
        return "\n".join([
            "FIELD_MODEL                    T05",
            "EPOCH                          %s" % DEFAULT_EPOCH,
            "DRIVER_FILE                    ts05_driver_C16.txt",
            "DST                            -22.0",
            "PDYN                           1.03",
            "IMF_BX                         -0.34",
            "IMF_BY                         -8.71",
            "IMF_BZ                         1.52",
            "SW_VX                          -359.6",
            "SW_N                           4.09",
            "G1                             0.0",
            "G2                             0.0",
            "G3                             0.0",
            "T05_W1                         0.21",
            "T05_W2                         0.08",
            "T05_W3                         0.26",
            "T05_W4                         0.09",
            "T05_W5                         0.02",
            "T05_W6                         0.03",
        ])
    raise ValueError("Unsupported field model %s" % model)


def render_input(template_path, output_path, run_id, model, points,
                 rigidities_gv, lon_res_deg, lat_res_deg, max_trace_time_s,
                 max_steps, max_trace_distance_re, args):
    """Render one immutable C16 run deck and reject forgotten placeholders."""
    text = template_path.read_text(errors="replace")
    replacements = {
        "__RUN_ID__": run_id,
        "__DIRMAP_LON_RES__": fmt_number(lon_res_deg),
        "__DIRMAP_LAT_RES__": fmt_number(lat_res_deg),
        "__CUTOFF_RIGIDITY_LIST_GV__": ",".join(
            fmt_number(value) for value in rigidities_gv),
        "__FIELD_BLOCK__": field_block(model),
        "__POINTS_BLOCK__": points_block(points),
        "__DT_TRACE__": fmt_number(args.dt_trace),
        "__ADAPTIVE_DT__": bool_token(args.adaptive_dt),
        "__MAX_STEPS__": str(int(max_steps)),
        "__MAX_TRACE_TIME__": fmt_number(max_trace_time_s),
        "__MAX_TRACE_DISTANCE_RE__": fmt_number(max_trace_distance_re),
        "__GRIDLESS_MPI_SCHEDULER__": args.scheduler,
        "__GRIDLESS_MPI_DYNAMIC_CHUNK__": str(args.dynamic_chunk),
        "__GRIDLESS_THREADS__": str(args.nt),
    }
    for placeholder, replacement in replacements.items():
        text = text.replace(placeholder, replacement)
    remaining = sorted(set(re.findall(r"__[A-Z0-9_]+__", text)))
    if remaining:
        raise RuntimeError("Unresolved template placeholders: %s" %
                           ", ".join(remaining))
    output_path.write_text(text)


def run_command(command, workdir, log_path):
    """Run AMPS while mirroring combined stdout/stderr into a durable log."""
    with log_path.open("w") as log:
        log.write("Command:\n  %s\n\n" % " ".join(command))
        log.flush()
        process = subprocess.Popen(
            command, cwd=str(workdir), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True)
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            log.write(line)
        return process.wait()


def normalize_variable_name(name):
    """Normalize Tecplot labels without weakening the required-column policy."""
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def parse_directional_access(path):
    """Parse one AMPS DIRECT_ACCESS file and detect malformed/duplicate rows.

    The parser accepts a Tecplot VARIABLES declaration split over multiple lines,
    but it requires the production C16 diagnostics by their normalized names.
    This avoids silently falling back to guessed column offsets when an older
    executable does not implement the termination audit output.
    """
    variables = []
    reading_variables = False
    samples_by_direction = defaultdict(list)
    seen_samples = set()

    with path.open("r", errors="replace") as stream:
        for line_number, raw in enumerate(stream, start=1):
            line = raw.strip()
            if not line:
                continue
            upper = line.upper()
            if upper.startswith("VARIABLES"):
                reading_variables = True
                variables.extend(normalize_variable_name(name)
                                 for name in re.findall(r'"([^"]+)"', line))
                continue
            if upper.startswith("ZONE"):
                reading_variables = False
                continue
            if reading_variables:
                quoted = re.findall(r'"([^"]+)"', line)
                if quoted:
                    variables.extend(normalize_variable_name(name)
                                     for name in quoted)
                    continue
                reading_variables = False
            if upper.startswith(("TITLE", "#", "!")):
                continue
            parts = line.replace(",", " ").split()
            if not variables or len(parts) < len(variables):
                continue
            record = dict(zip(variables, parts))

            required = ("lon_deg", "lat_deg", "rigidity_gv", "access_state",
                        "termination_code", "trace_time_s",
                        "trace_distance_re", "trace_steps")
            missing = [name for name in required if name not in record]
            if missing:
                raise RuntimeError(
                    "%s lacks required C16 columns: %s" %
                    (path, ", ".join(missing)))
            try:
                lon_deg = float(record["lon_deg"])
                lat_deg = float(record["lat_deg"])
                rigidity_gv = float(record["rigidity_gv"])
                state = int(float(record["access_state"]))
                termination_code = int(float(record["termination_code"]))
                trace_time_s = float(record["trace_time_s"])
                trace_distance_re = float(record["trace_distance_re"])
                trace_steps = int(float(record["trace_steps"]))
            except (ValueError, OverflowError) as exc:
                raise RuntimeError("Malformed C16 row %d in %s: %s" %
                                   (line_number, path, exc))

            numeric = (lon_deg, lat_deg, rigidity_gv, trace_time_s,
                       trace_distance_re, float(trace_steps))
            if not all(math.isfinite(value) for value in numeric):
                raise RuntimeError("Non-finite C16 diagnostic at row %d in %s" %
                                   (line_number, path))
            if rigidity_gv <= 0.0 or trace_time_s < 0.0 or \
                    trace_distance_re < 0.0 or trace_steps < 0:
                raise RuntimeError("Out-of-range C16 diagnostic at row %d in %s" %
                                   (line_number, path))
            if state not in STATE_NAMES:
                raise RuntimeError("Unknown access_state=%d at row %d in %s" %
                                   (state, line_number, path))

            dkey = direction_key(lon_deg, lat_deg)
            skey = (dkey, rigidity_key(rigidity_gv))
            if skey in seen_samples:
                raise RuntimeError(
                    "Duplicate direction/rigidity sample lon=%g lat=%g R=%g in %s" %
                    (lon_deg, lat_deg, rigidity_gv, path))
            seen_samples.add(skey)

            sample = {
                "lon_deg": lon_deg,
                "lat_deg": lat_deg,
                "rigidity_GV": rigidity_gv,
                "energy_MeV": _optional_float(record, "energy_mev"),
                "access_state": state,
                "termination_code": termination_code,
                "trace_time_s": trace_time_s,
                "trace_distance_Re": trace_distance_re,
                "trace_steps": trace_steps,
                "retry_count": _optional_int(record, "retry_count"),
                "reported_allowed": _optional_int(record, "allowed"),
                "reported_unresolved": _optional_int(record, "unresolved"),
            }
            samples_by_direction[dkey].append(sample)

    if not samples_by_direction:
        raise RuntimeError("No DIRECT_ACCESS rows parsed from %s" % path)
    for rows in samples_by_direction.values():
        rows.sort(key=lambda item: item["rigidity_GV"])
    return {"columns": tuple(variables), "directions": dict(samples_by_direction)}


def _optional_float(record, name):
    """Parse an optional finite diagnostic, returning None when it is absent."""
    if name not in record:
        return None
    try:
        value = float(record[name])
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _optional_int(record, name):
    """Parse an optional integer-like Tecplot diagnostic."""
    value = _optional_float(record, name)
    return None if value is None else int(round(value))


def find_access_file(workdir, point_id):
    """Locate the documented GRIDLESS name, then tolerate compatible variants."""
    preferred = workdir / ("cutoff_gridless_dir_access_point_%04d.dat" % point_id)
    if preferred.exists():
        return preferred
    matches = sorted(workdir.glob("*dir*access*point*%04d*.dat" % point_id))
    return matches[0] if matches else preferred


def validate_sample_contract(sample, limits, diagnostic_slack):
    """Check one state/reason pair and its limit-specific diagnostic evidence.

    The function returns error labels instead of raising so the final audit CSV
    remains complete even when many samples expose the same model defect.
    """
    errors = []
    code = sample["termination_code"]
    state = sample["access_state"]
    contract = TERMINATION_CONTRACT.get(code)
    if contract is None:
        errors.append("unknown_termination_code")
        return errors
    expected_state = contract[1]
    if state != expected_state:
        errors.append("state_termination_mismatch")
    if code in FATAL_CODES:
        errors.append("fatal_numerical_termination")

    # If AMPS emits redundant binary columns, C16 checks them too.  They are
    # optional because access_state is the authoritative production quantity.
    if sample["reported_allowed"] is not None:
        expected_allowed = 1 if state == 1 else 0
        if sample["reported_allowed"] != expected_allowed:
            errors.append("reported_allowed_mismatch")
    if sample["reported_unresolved"] is not None:
        expected_unresolved = 1 if state == 2 else 0
        if sample["reported_unresolved"] != expected_unresolved:
            errors.append("reported_unresolved_mismatch")

    # Resource-limit rows must carry evidence that the corresponding configured
    # budget was actually reached.  A generous default slack permits one adaptive
    # integration-step overshoot/roundoff but catches premature misclassification.
    lower_factor = max(0.0, 1.0 - float(diagnostic_slack))
    if code == 3 and sample["trace_time_s"] < \
            limits["max_trace_time"] * lower_factor:
        errors.append("time_limit_not_reached")
    if code == 4 and sample["trace_steps"] < \
            max(1, int(math.ceil(limits["max_steps"] * lower_factor))):
        errors.append("step_limit_not_reached")
    if code == 5:
        if limits["max_trace_distance"] <= 0.0:
            errors.append("distance_limit_with_disabled_cap")
        elif sample["trace_distance_Re"] < \
                limits["max_trace_distance"] * lower_factor:
            errors.append("distance_limit_not_reached")
    return errors


def match_requested_rigidity(value, requested):
    """Find a requested node with a tolerance tight enough to expose grid drift."""
    for target in requested:
        if math.isclose(value, target, rel_tol=5.0e-10, abs_tol=1.0e-12):
            return target
    return None


def audit_access_cube(parsed, point, expected_directions, rigidities_gv,
                      run_meta, limits, diagnostic_slack):
    """Audit completeness and state/reason closure for one observation point."""
    got_directions = set(parsed["directions"])
    missing_directions = expected_directions - got_directions
    extra_directions = got_directions - expected_directions
    sample_rows = []
    counter = Counter()
    sample_map = {}
    missing_rigidities = 0
    extra_rigidities = 0

    for dkey in sorted(expected_directions & got_directions):
        matched = set()
        for sample in parsed["directions"][dkey]:
            requested = match_requested_rigidity(sample["rigidity_GV"],
                                                 rigidities_gv)
            if requested is None:
                extra_rigidities += 1
                continue
            requested_key = rigidity_key(requested)
            if requested_key in matched:
                # parse_directional_access already catches exact duplicates;
                # this branch catches two distinct values within match tolerance.
                extra_rigidities += 1
                continue
            matched.add(requested_key)
            errors = validate_sample_contract(sample, limits, diagnostic_slack)
            code = sample["termination_code"]
            contract = TERMINATION_CONTRACT.get(code)
            counter["samples"] += 1
            counter["state_%d" % sample["access_state"]] += 1
            counter["code_%d" % code] += 1
            counter["contract_errors"] += len(errors)
            if code in FATAL_CODES:
                counter["fatal"] += 1

            key = (point["point_id"], dkey[0], dkey[1], requested_key)
            sample_map[key] = sample
            row = {
                "case_type": run_meta["case_type"],
                "probe": run_meta.get("probe", ""),
                "mover": run_meta["mover"],
                "field_model": run_meta["model"],
                "trace_budget_s": limits["max_trace_time"],
                "point_id": point["point_id"],
                "obs_lon_deg": point["obs_lon_deg"],
                "obs_lat_deg": point["obs_lat_deg"],
                "sky_lon_deg": sample["lon_deg"],
                "sky_lat_deg": sample["lat_deg"],
                "rigidity_GV": sample["rigidity_GV"],
                "access_state": sample["access_state"],
                "access_state_name": STATE_NAMES[sample["access_state"]],
                "termination_code": code,
                "termination_name": contract[0] if contract else "UNKNOWN",
                "expected_access_state": contract[1] if contract else "",
                "trace_time_s": sample["trace_time_s"],
                "trace_distance_Re": sample["trace_distance_Re"],
                "trace_steps": sample["trace_steps"],
                "retry_count": ("" if sample["retry_count"] is None
                                  else sample["retry_count"]),
                "passed": int(not errors),
                "errors": ";".join(errors),
            }
            sample_rows.append(row)
        missing_rigidities += len(rigidities_gv) - len(matched)

    expected_samples = len(expected_directions) * len(rigidities_gv)
    summary = {
        "case_type": run_meta["case_type"],
        "probe": run_meta.get("probe", ""),
        "mover": run_meta["mover"],
        "field_model": run_meta["model"],
        "trace_budget_s": limits["max_trace_time"],
        "max_steps": limits["max_steps"],
        "max_trace_distance_Re": limits["max_trace_distance"],
        "point_id": point["point_id"],
        "obs_lon_deg": point["obs_lon_deg"],
        "obs_lat_deg": point["obs_lat_deg"],
        "expected_directions": len(expected_directions),
        "actual_directions": len(got_directions),
        "missing_directions": len(missing_directions),
        "extra_directions": len(extra_directions),
        "expected_samples": expected_samples,
        "audited_samples": counter["samples"],
        "missing_rigidities": missing_rigidities,
        "extra_rigidities": extra_rigidities,
        "allowed": counter["state_1"],
        "forbidden": counter["state_0"],
        "unresolved": counter["state_2"],
        "unresolved_fraction": (float(counter["state_2"]) / counter["samples"]
                                if counter["samples"] else 1.0),
        "contract_errors": counter["contract_errors"],
        "fatal_terminations": counter["fatal"],
    }
    coverage_errors = (len(missing_directions) + len(extra_directions) +
                       missing_rigidities + extra_rigidities)
    summary["coverage_errors"] = coverage_errors
    summary["passed"] = int(coverage_errors == 0 and
                            counter["contract_errors"] == 0)
    return summary, sample_rows, sample_map, counter


def audit_run(run_spec, points, expected_directions, rigidities_gv,
              diagnostic_slack):
    """Audit all point files in one AMPS invocation and retain every failure."""
    summaries = []
    sample_rows = []
    sample_map = {}
    counts = Counter()
    messages = []
    limits = run_spec["limits"]

    for point in points:
        access_path = find_access_file(run_spec["workdir"], point["point_id"])
        if not access_path.exists():
            messages.append("missing DIRECT_ACCESS output for point %d: %s" %
                            (point["point_id"], access_path))
            continue
        try:
            parsed = parse_directional_access(access_path)
            summary, rows, point_map, point_counts = audit_access_cube(
                parsed, point, expected_directions, rigidities_gv, run_spec,
                limits, diagnostic_slack)
        except Exception as exc:
            messages.append("point %d output audit failed: %s" %
                            (point["point_id"], exc))
            continue
        summaries.append(summary)
        sample_rows.extend(rows)
        sample_map.update(point_map)
        counts.update(point_counts)
        if not summary["passed"]:
            messages.append(
                "point %d coverage=%d contract_errors=%d fatal=%d" %
                (point["point_id"], summary["coverage_errors"],
                 summary["contract_errors"], summary["fatal_terminations"]))

    if len(summaries) != len(points):
        messages.append("audited %d/%d expected point files" %
                        (len(summaries), len(points)))
    if run_spec["case_type"] == "baseline":
        # A matrix containing only one physical state would exercise bookkeeping
        # but not both sides of the access classifier, so it is not a useful
        # model-validation run.  This gate is applied to the full run, not to each
        # latitude separately, because real field topology can be locally one-sided.
        if counts["state_1"] == 0:
            messages.append("baseline contains no physical ALLOWED samples")
        if counts["state_0"] == 0:
            messages.append("baseline contains no physical FORBIDDEN samples")

    passed = (not messages and all(row["passed"] for row in summaries))
    return {
        "passed": passed,
        "summaries": summaries,
        "sample_rows": sample_rows,
        "sample_map": sample_map,
        "counts": counts,
        "messages": messages,
    }


def compare_budget_pair(short_spec, long_spec, short_audit, long_audit):
    """Apply the monotone-information rule to two increasing trace budgets."""
    short_map = short_audit["sample_map"]
    long_map = long_audit["sample_map"]
    all_keys = sorted(set(short_map) | set(long_map))
    rows = []
    summary = Counter()

    for key in all_keys:
        left = short_map.get(key)
        right = long_map.get(key)
        error = ""
        if left is None:
            error = "missing_short_sample"
        elif right is None:
            error = "missing_long_sample"
        else:
            old_state = left["access_state"]
            new_state = right["access_state"]
            if old_state in (0, 1) and new_state != old_state:
                error = ("resolved_state_flip" if new_state in (0, 1)
                         else "resolved_became_unresolved")
            elif old_state == 2 and new_state in (0, 1):
                summary["newly_resolved"] += 1
            elif old_state == 2 and new_state == 2:
                summary["still_unresolved"] += 1
            else:
                summary["stable_resolved"] += 1
        if error:
            summary["regressions"] += 1
        summary["samples"] += 1
        point_id, sky_lon, sky_lat, rigidity = key
        rows.append({
            "mover": short_spec["mover"],
            "field_model": short_spec["model"],
            "short_budget_s": short_spec["limits"]["max_trace_time"],
            "long_budget_s": long_spec["limits"]["max_trace_time"],
            "point_id": point_id,
            "sky_lon_deg": sky_lon,
            "sky_lat_deg": sky_lat,
            "rigidity_GV": rigidity,
            "short_state": "" if left is None else left["access_state"],
            "long_state": "" if right is None else right["access_state"],
            "short_code": "" if left is None else left["termination_code"],
            "long_code": "" if right is None else right["termination_code"],
            "passed": int(not error),
            "error": error,
        })

    summary_row = {
        "mover": short_spec["mover"],
        "field_model": short_spec["model"],
        "short_budget_s": short_spec["limits"]["max_trace_time"],
        "long_budget_s": long_spec["limits"]["max_trace_time"],
        "samples": summary["samples"],
        "stable_resolved": summary["stable_resolved"],
        "newly_resolved": summary["newly_resolved"],
        "still_unresolved": summary["still_unresolved"],
        "regressions": summary["regressions"],
        "passed": int(summary["regressions"] == 0),
    }
    return summary_row, rows


def probe_result(run_spec, audit, minimum_fraction):
    """Evaluate whether a deliberate probe reached its intended code often enough."""
    expected_code = run_spec["expected_code"]
    total = audit["counts"]["samples"]
    target = audit["counts"]["code_%d" % expected_code]
    fraction = float(target) / total if total else 0.0
    passed = audit["passed"] and fraction >= minimum_fraction
    row = {
        "mover": run_spec["mover"],
        "probe": run_spec["probe"],
        "expected_termination_code": expected_code,
        "expected_termination_name": TERMINATION_CONTRACT[expected_code][0],
        "expected_access_state": TERMINATION_CONTRACT[expected_code][1],
        "samples": total,
        "target_code_samples": target,
        "target_fraction": fraction,
        "minimum_target_fraction": minimum_fraction,
        "passed": int(passed),
    }
    return row


def termination_count_rows(run_spec, audit):
    """Expand a run's counter into one stable CSV row for every code 0–9."""
    rows = []
    total = audit["counts"]["samples"]
    for code in sorted(TERMINATION_CONTRACT):
        contract = TERMINATION_CONTRACT[code]
        count = audit["counts"]["code_%d" % code]
        rows.append({
            "case_type": run_spec["case_type"],
            "probe": run_spec.get("probe", ""),
            "mover": run_spec["mover"],
            "field_model": run_spec["model"],
            "trace_budget_s": run_spec["limits"]["max_trace_time"],
            "termination_code": code,
            "termination_name": contract[0],
            "expected_access_state": contract[1],
            "category": contract[2],
            "count": count,
            "fraction": float(count) / total if total else 0.0,
        })
    return rows


def write_dict_csv(rows, path, fieldnames=None):
    """Write unioned dictionary columns while preserving first-seen order."""
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path):
    """Calculate a streaming SHA-256 for input/reference provenance."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_t05_driver(path, required_epoch=DEFAULT_EPOCH):
    """Validate the compact T05 driver before any empirical-field run.

    Each data row must contain one ISO timestamp plus the 19 numeric quantities
    used by the existing AMPS TS05 reader.  C16 also requires strict time order
    and an exact row at its fixed epoch; otherwise a copied or truncated driver
    could change the validation field while leaving the test name unchanged.
    """
    rows = []
    with path.open("r", errors="replace") as stream:
        for line_number, raw in enumerate(stream, start=1):
            line = raw.strip()
            if not line or line.startswith(("#", "!")):
                continue
            parts = line.split()
            if len(parts) != 20:
                raise RuntimeError(
                    "T05 driver row %d must have timestamp + 19 values; got %d" %
                    (line_number, len(parts)))
            try:
                timestamp = dt.datetime.fromisoformat(parts[0].replace("Z", "+00:00"))
                if timestamp.tzinfo is not None:
                    timestamp = timestamp.astimezone(dt.timezone.utc).replace(tzinfo=None)
                values = [float(token) for token in parts[1:]]
            except ValueError as exc:
                raise RuntimeError("invalid T05 driver row %d: %s" %
                                   (line_number, exc))
            if not all(math.isfinite(value) for value in values):
                raise RuntimeError("non-finite T05 driver value at row %d" % line_number)
            rows.append((timestamp, parts[0]))
    if len(rows) < 3:
        raise RuntimeError("T05 driver must bracket the epoch with at least 3 rows")
    for left, right in zip(rows[:-1], rows[1:]):
        if not right[0] > left[0]:
            raise RuntimeError("T05 driver timestamps are not strictly increasing")
    normalized_epoch = required_epoch.rstrip("Z")
    if normalized_epoch not in {text.rstrip("Z") for _, text in rows}:
        raise RuntimeError("T05 driver has no exact row at %s" % required_epoch)
    required_time = dt.datetime.fromisoformat(normalized_epoch)
    if required_time.tzinfo is not None:
        required_time = required_time.astimezone(dt.timezone.utc).replace(tzinfo=None)
    if not rows[0][0] <= required_time <= rows[-1][0]:
        raise RuntimeError("T05 driver does not bracket %s" % required_epoch)
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "row_count": len(rows),
        "first_epoch": rows[0][1],
        "last_epoch": rows[-1][1],
        "selected_epoch": required_epoch,
    }


def validate_reference_files(script_dir):
    """Prove the checked-in CSV reference solutions match executable constants."""
    termination_path = script_dir / "reference_C16_termination_contract.csv"
    with termination_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(TERMINATION_CONTRACT):
        raise RuntimeError("termination reference must contain exactly codes 0-9")
    for row in rows:
        code = int(row["termination_code"])
        if code not in TERMINATION_CONTRACT:
            raise RuntimeError("termination reference contains unknown code %d" % code)
        expected = TERMINATION_CONTRACT[code]
        observed = (row["termination_name"], int(row["expected_access_state"]),
                    row["category"], row["validation_policy"])
        if observed != expected:
            raise RuntimeError("termination reference disagrees at code %d" % code)

    probe_path = script_dir / "reference_C16_probe_expectations.csv"
    with probe_path.open(newline="") as stream:
        probe_rows = list(csv.DictReader(stream))
    if {row["probe"] for row in probe_rows} != set(DEFAULT_PROBES):
        raise RuntimeError("probe reference does not contain the three C16 probes")
    controlled_names = {
        "time_limit": "MAX_TRACE_TIME",
        "step_limit": "MAX_STEPS",
        "distance_limit": "MAX_TRACE_DISTANCE",
    }
    for row in probe_rows:
        probe = DEFAULT_PROBES[row["probe"]]
        controlled_name = controlled_names[row["probe"]]
        configured_value = {
            "MAX_TRACE_TIME": probe["max_trace_time"],
            "MAX_STEPS": probe["max_steps"],
            "MAX_TRACE_DISTANCE": probe["max_trace_distance"],
        }[controlled_name]
        if (int(row["expected_termination_code"]) != probe["expected_code"] or
                row["expected_termination_name"] !=
                TERMINATION_CONTRACT[probe["expected_code"]][0] or
                int(row["expected_access_state"]) != 2 or
                row["controlled_limit"] != controlled_name or
                not math.isclose(float(row["configured_value"]),
                                 float(configured_value), rel_tol=0.0,
                                 abs_tol=1.0e-15) or
                not math.isclose(float(row["minimum_target_fraction"]), 0.90,
                                 rel_tol=0.0, abs_tol=1.0e-15)):
            raise RuntimeError("probe reference disagrees for %s" % row["probe"])
    return True


def validate_template_contract(template_path):
    """Reject parser-incompatible regressions before any expensive AMPS launch."""
    text = template_path.read_text(errors="replace")
    required_patterns = {
        "DIRECT_ACCESS algorithm":
            r"^CUTOFF_SEARCH_ALGORITHM\s+DIRECT_ACCESS\s*$",
        "VERTICAL parser contract":
            r"^CUTOFF_SAMPLING\s+VERTICAL\s*$",
        "explicit rigidity list":
            r"^CUTOFF_RIGIDITY_LIST_GV\s+__CUTOFF_RIGIDITY_LIST_GV__\s*$",
        "fixed access grid":
            r"^CUTOFF_DIRECT_ACCESS_ADAPTIVE\s+F\s*$",
        "three-state limit policy":
            r"^CUTOFF_TRACE_LIMIT_POLICY\s+UNRESOLVED\s*$",
        "step budget": r"^MAX_STEPS\s+__MAX_STEPS__\s*$",
        "time budget": r"^MAX_TRACE_TIME\s+__MAX_TRACE_TIME__\s*$",
        "distance budget":
            r"^MAX_TRACE_DISTANCE\s+__MAX_TRACE_DISTANCE_RE__\s*$",
    }
    for label, pattern in required_patterns.items():
        if not re.search(pattern, text, re.MULTILINE):
            raise RuntimeError("C16 template lacks %s" % label)

    # These optional controls caused real compatibility failures in neighboring
    # tests when used against the baseline parser.  C16 must stay deployable on
    # that parser because its validation requires no such extensions.
    unsupported = (
        "CUTOFF_UNRESOLVED_EXTENSION_PASSES",
        "CUTOFF_UNRESOLVED_EXTENSION_FACTOR",
        "CUTOFF_DEBUG_EXIT_TIME",
        "CUTOFF_DEBUG_EXIT_STEP",
        "CUTOFF_DEBUG_EXIT_DISTANCE",
        "TRAP_DETECTION",
        "TRAP_DRIFT_DETECTION",
    )
    for keyword in unsupported:
        if re.search(r"^\s*%s\b" % re.escape(keyword), text, re.MULTILINE):
            raise RuntimeError("C16 template contains unsupported keyword %s" %
                               keyword)
    return True


def _write_synthetic_access(path, directions, rigidities, code_for_sample,
                            corrupt_state=False, drop_last=False,
                            duplicate_first=False):
    """Create a production-shaped Tecplot cube for the package self-test."""
    lines = [
        'VARIABLES="lon_deg" "lat_deg" "rigidity_GV" "energy_MeV" '
        '"access_state" "termination_code" "trace_time_s" '
        '"trace_distance_Re" "trace_steps" "retry_count"',
        'ZONE T="synthetic C16"',
    ]
    data_lines = []
    index = 0
    for lon_deg, lat_deg in sorted(directions):
        for rigidity in rigidities:
            code = int(code_for_sample(index))
            state = TERMINATION_CONTRACT[code][1]
            if corrupt_state and index == 0:
                state = (state + 1) % 3
            # Synthetic limit diagnostics are chosen above all self-test limits.
            data_lines.append("%g %g %.12g %.12g %d %d 600 10 2000000 0" %
                              (lon_deg, lat_deg, rigidity, rigidity * 100.0,
                               state, code))
            index += 1
    if drop_last:
        data_lines.pop()
    if duplicate_first:
        data_lines.append(data_lines[0])
    path.write_text("\n".join(lines + data_lines) + "\n")


def self_test():
    """Exercise all critical C16 guards with small synthetic access products."""
    script_dir = Path(__file__).resolve().parent
    try:
        validate_template_contract(script_dir / "AMPS_PARAM_C16_gridless.in")
        validate_reference_files(script_dir)
        validate_t05_driver(script_dir / "data" / "ts05_driver_C16.txt")
    except Exception as exc:
        print("C16 self-test failed: %s" % exc, file=sys.stderr)
        return 1

    # Check all ten reference pairs, including the fatal-code distinction.
    generous_limits = {"max_trace_time": 1.0, "max_steps": 1,
                       "max_trace_distance": 0.01}
    for code, contract in TERMINATION_CONTRACT.items():
        sample = {
            "termination_code": code,
            "access_state": contract[1],
            "trace_time_s": 10.0,
            "trace_distance_Re": 10.0,
            "trace_steps": 10,
            "reported_allowed": None,
            "reported_unresolved": None,
        }
        errors = validate_sample_contract(sample, generous_limits, 0.10)
        expected_errors = (["fatal_numerical_termination"]
                           if code in FATAL_CODES else [])
        if errors != expected_errors:
            print("C16 self-test failed: code %d contract result %r" %
                  (code, errors), file=sys.stderr)
            return 1
        sample["access_state"] = (contract[1] + 1) % 3
        if "state_termination_mismatch" not in validate_sample_contract(
                sample, generous_limits, 0.10):
            print("C16 self-test failed: wrong state accepted for code %d" % code,
                  file=sys.stderr)
            return 1

    with tempfile.TemporaryDirectory(prefix="c16_selftest_") as raw:
        root = Path(raw)
        directions = expected_direction_keys(90.0, 90.0)
        rigidities = [1.0, 2.0]
        point = build_points([0.0], [0.0], 9000.0)[0]
        access_path = root / "cutoff_gridless_dir_access_point_0000.dat"
        run_spec = {
            "case_type": "baseline", "mover": "BORIS", "model": "DIPOLE",
            "limits": {"max_trace_time": 600.0, "max_steps": 2000000,
                       "max_trace_distance": 0.0},
        }

        # Alternate physical allowed/forbidden reasons so the cube exercises
        # completeness, closure, and the physical-diversity gate simultaneously.
        _write_synthetic_access(
            access_path, directions, rigidities,
            lambda index: 0 if index % 2 == 0 else 1)
        parsed = parse_directional_access(access_path)
        summary, _, sample_map, _ = audit_access_cube(
            parsed, point, directions, rigidities, run_spec,
            run_spec["limits"], 0.10)
        if not summary["passed"] or len(sample_map) != 24:
            print("C16 self-test failed: valid complete cube was rejected",
                  file=sys.stderr)
            return 1

        _write_synthetic_access(
            access_path, directions, rigidities, lambda index: 0,
            drop_last=True)
        parsed = parse_directional_access(access_path)
        truncated, _, _, _ = audit_access_cube(
            parsed, point, directions, rigidities, run_spec,
            run_spec["limits"], 0.10)
        if truncated["passed"] or truncated["missing_rigidities"] != 1:
            print("C16 self-test failed: truncated cube was not rejected",
                  file=sys.stderr)
            return 1

        _write_synthetic_access(
            access_path, directions, rigidities, lambda index: 0,
            duplicate_first=True)
        try:
            parse_directional_access(access_path)
        except RuntimeError:
            pass
        else:
            print("C16 self-test failed: duplicate cube row was not rejected",
                  file=sys.stderr)
            return 1

        # Directly exercise the convergence partial order.  UNRESOLVED may gain
        # information, but an established physical classification is immutable.
        key = (0, 0.0, 0.0, rigidity_key(1.0))
        short_spec = {"mover": "BORIS", "model": "DIPOLE",
                      "limits": {"max_trace_time": 60.0}}
        long_spec = {"mover": "BORIS", "model": "DIPOLE",
                     "limits": {"max_trace_time": 600.0}}
        unresolved = {"access_state": 2, "termination_code": 3}
        allowed = {"access_state": 1, "termination_code": 0}
        forbidden = {"access_state": 0, "termination_code": 1}
        good_summary, _ = compare_budget_pair(
            short_spec, long_spec, {"sample_map": {key: unresolved}},
            {"sample_map": {key: allowed}})
        if not good_summary["passed"] or good_summary["newly_resolved"] != 1:
            print("C16 self-test failed: unresolved-to-physical was rejected",
                  file=sys.stderr)
            return 1
        bad_summary, _ = compare_budget_pair(
            short_spec, long_spec, {"sample_map": {key: allowed}},
            {"sample_map": {key: forbidden}})
        if bad_summary["passed"] or bad_summary["regressions"] != 1:
            print("C16 self-test failed: physical state flip was accepted",
                  file=sys.stderr)
            return 1

        # A deliberate TIME_LIMIT probe must close to UNRESOLVED and satisfy
        # the target-fraction gate.  Corrupting one state must make it fail.
        probe_spec = {
            "case_type": "probe", "probe": "time_limit", "mover": "BORIS",
            "model": "DIPOLE", "expected_code": 3,
            "limits": {"max_trace_time": 0.01, "max_steps": 2000000,
                       "max_trace_distance": 0.0},
        }
        _write_synthetic_access(access_path, directions, rigidities,
                                lambda index: 3)
        probe_dir = root / "probe"
        probe_dir.mkdir()
        shutil.copy2(access_path,
                     probe_dir / "cutoff_gridless_dir_access_point_0000.dat")
        probe_spec["workdir"] = probe_dir
        audit = audit_run(probe_spec, [point], directions, rigidities, 0.10)
        if not probe_result(probe_spec, audit, 0.90)["passed"]:
            print("C16 self-test failed: valid time probe was rejected",
                  file=sys.stderr)
            return 1
        _write_synthetic_access(
            probe_dir / "cutoff_gridless_dir_access_point_0000.dat",
            directions, rigidities, lambda index: 3, corrupt_state=True)
        corrupted = audit_run(probe_spec, [point], directions, rigidities, 0.10)
        if probe_result(probe_spec, corrupted, 0.90)["passed"]:
            print("C16 self-test failed: corrupt time probe was accepted",
                  file=sys.stderr)
            return 1

    print("C16 self-test: PASS")
    return 0


def apply_profile(args):
    """Fill only unspecified CLI controls from the selected cost profile."""
    defaults = PROFILE_DEFAULTS[args.profile]
    for name, value in defaults.items():
        if getattr(args, name) is None:
            setattr(args, name, value)


def parse_args():
    """Define the complete command-line interface in one discoverable place."""
    parser = argparse.ArgumentParser(
        description="C16 trajectory termination/access-state closure validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Recommended model-validation run:
          python3 srcEarth/test/C16/run_C16.py --profile ROUTINE --amps ./amps -np 4 -nt 16

        Package-only verification and command inspection:
          python3 srcEarth/test/C16/run_C16.py --self-test
          python3 srcEarth/test/C16/run_C16.py --profile SMOKE --dry-run

        Reanalyze an existing output matrix without launching AMPS:
          python3 srcEarth/test/C16/run_C16.py --profile ROUTINE --skip-run
        """))
    parser.add_argument("--profile", choices=sorted(PROFILE_DEFAULTS),
                        default="ROUTINE", help="cost/resolution profile; default: ROUTINE")
    parser.add_argument("-np", type=int, default=4,
                        help="MPI ranks passed to the launcher; default: 4")
    parser.add_argument("-nt", type=int, default=16,
                        help="GRIDLESS worker threads; default: 16")
    parser.add_argument("--models", default=None,
                        help="comma-separated DIPOLE,T05 override")
    parser.add_argument("--movers", default="BORIS",
                        help="comma-separated mover list; default: BORIS")
    parser.add_argument("--budgets", default=None,
                        help="increasing baseline trace-time budgets in seconds")
    parser.add_argument("--lons", default="0",
                        help="observation longitudes, degrees; default: 0")
    parser.add_argument("--lats", default=None,
                        help="observation latitudes, degrees; profile default")
    parser.add_argument("--alt", type=float, default=9000.0,
                        help="observation altitude, km; default: 9000")
    parser.add_argument("--rigidities", default=None,
                        help="exact DIRECT_ACCESS rigidity list in GV")
    parser.add_argument("--dir-lon-res", type=float, default=None,
                        help="baseline sky longitude resolution, degrees")
    parser.add_argument("--dir-lat-res", type=float, default=None,
                        help="baseline sky latitude resolution, degrees")
    parser.add_argument("--max-steps", type=int, default=None,
                        help="ordinary baseline integration-step cap")
    parser.add_argument("--dt-trace", type=float, default=0.25,
                        help="maximum base trace step, seconds; default: 0.25")
    parser.add_argument("--adaptive-dt", default="T",
                        help="AMPS ADAPTIVE_DT setting T/F; default: T")
    parser.add_argument("--max-final-unresolved-fraction", type=float, default=None,
                        help="per-model final-budget unresolved ceiling")
    parser.add_argument("--limit-diagnostic-slack", type=float, default=0.10,
                        help="fractional lower slack for limit evidence; default: 0.10")
    parser.add_argument("--probe-target-fraction", type=float, default=0.90,
                        help="minimum intended-code fraction per probe; default: 0.90")
    parser.add_argument("--probe-dir-lon-res", type=float, default=90.0,
                        help="probe sky longitude resolution; default: 90")
    parser.add_argument("--probe-dir-lat-res", type=float, default=90.0,
                        help="probe sky latitude resolution; default: 90")
    parser.add_argument("--probe-rigidities", default="1,10",
                        help="probe rigidity list in GV; default: 1,10")
    parser.add_argument("--skip-probes", action="store_true",
                        help="omit deliberate limit probes (not recommended for validation)")
    parser.add_argument("--scheduler", default="DYNAMIC",
                        choices=("DYNAMIC", "BLOCK_CYCLIC", "STATIC"),
                        help="GRIDLESS MPI scheduler; default: DYNAMIC")
    parser.add_argument("--dynamic-chunk", type=int, default=0,
                        help="dynamic task chunk; 0 selects AMPS heuristic")
    parser.add_argument("--amps", default="./amps",
                        help="AMPS executable path relative to launch directory")
    parser.add_argument("--mpirun", default="mpirun",
                        help="MPI launcher executable; default: mpirun")
    parser.add_argument("--workdir", default="test_output/C16_exit_closure",
                        help="base output directory")
    parser.add_argument("--skip-run", action="store_true",
                        help="analyze existing outputs without launching AMPS")
    parser.add_argument("--keep", action="store_true",
                        help="keep existing workdir before a new run")
    parser.add_argument("--dry-run", action="store_true",
                        help="render inputs and commands without launching/analyzing")
    parser.add_argument("--self-test", action="store_true",
                        help="exercise package logic without launching AMPS")
    parser.add_argument("--version", action="version",
                        version="C16 runner schema %d (%s)" %
                                (RUNNER_SCHEMA_VERSION, RUNNER_RELEASE))
    return parser.parse_args()


def validate_args(args):
    """Validate numeric and grid invariants before deleting or creating output."""
    if args.np < 1 or args.nt < 1:
        raise SystemExit("-np and -nt must be >= 1")
    if args.dynamic_chunk < 0:
        raise SystemExit("--dynamic-chunk must be >= 0")
    if args.alt <= 0.0 or args.dt_trace <= 0.0:
        raise SystemExit("--alt and --dt-trace must be positive")
    if args.max_steps < 1:
        raise SystemExit("--max-steps must be >= 1")
    for name in ("dir_lon_res", "dir_lat_res", "probe_dir_lon_res",
                 "probe_dir_lat_res"):
        value = getattr(args, name)
        if value <= 0.0:
            raise SystemExit("--%s must be positive" % name.replace("_", "-"))
    for name, total in (("dir_lon_res", 360.0), ("dir_lat_res", 180.0),
                        ("probe_dir_lon_res", 360.0),
                        ("probe_dir_lat_res", 180.0)):
        value = getattr(args, name)
        if not math.isclose(total / value, round(total / value),
                            rel_tol=0.0, abs_tol=1.0e-9):
            raise SystemExit("--%s must divide %g exactly" %
                             (name.replace("_", "-"), total))
    for name in ("max_final_unresolved_fraction", "limit_diagnostic_slack",
                 "probe_target_fraction"):
        value = getattr(args, name)
        if not 0.0 <= value <= 1.0:
            raise SystemExit("--%s must be in [0,1]" % name.replace("_", "-"))
    bool_token(args.adaptive_dt)


def make_run_specs(base_workdir, movers, models, budgets, points, probe_points,
                   args):
    """Create the complete run matrix before launching any external process."""
    specs = []
    for mover in movers:
        for model in models:
            for budget in budgets:
                workdir = (base_workdir / mover.lower() / model.lower() /
                           budget_label(budget))
                specs.append({
                    "case_type": "baseline",
                    "probe": "",
                    "mover": mover,
                    "model": model,
                    "points": points,
                    "rigidities": args.rigidity_list,
                    "lon_res": args.dir_lon_res,
                    "lat_res": args.dir_lat_res,
                    "limits": {
                        "max_trace_time": budget,
                        "max_steps": args.max_steps,
                        "max_trace_distance": 0.0,
                    },
                    "workdir": workdir,
                })
        if not args.skip_probes:
            for probe_name, probe in DEFAULT_PROBES.items():
                specs.append({
                    "case_type": "probe",
                    "probe": probe_name,
                    "mover": mover,
                    "model": "DIPOLE",
                    "points": probe_points,
                    "rigidities": args.probe_rigidity_list,
                    "lon_res": args.probe_dir_lon_res,
                    "lat_res": args.probe_dir_lat_res,
                    "limits": {
                        "max_trace_time": probe["max_trace_time"],
                        "max_steps": probe["max_steps"],
                        "max_trace_distance": probe["max_trace_distance"],
                    },
                    "expected_code": probe["expected_code"],
                    "workdir": (base_workdir / mover.lower() / "probes" /
                                probe_name),
                })
    return specs


def prepare_run(run_spec, template_path, driver_path, amps_path, args):
    """Create one work directory, rendered deck, local driver, and command."""
    workdir = run_spec["workdir"]
    workdir.mkdir(parents=True, exist_ok=True)
    input_path = workdir / "AMPS_PARAM_C16.in"
    run_id = "C16_%s_%s_%s" % (
        run_spec["mover"].lower(), run_spec["model"].lower(),
        run_spec.get("probe") or budget_label(
            run_spec["limits"]["max_trace_time"]))
    render_input(
        template_path, input_path, run_id, run_spec["model"],
        run_spec["points"], run_spec["rigidities"], run_spec["lon_res"],
        run_spec["lat_res"], run_spec["limits"]["max_trace_time"],
        run_spec["limits"]["max_steps"],
        run_spec["limits"]["max_trace_distance"], args)
    if run_spec["model"] == "T05":
        shutil.copy2(driver_path, workdir / "ts05_driver_C16.txt")

    command = [
        args.mpirun, "-np", str(args.np), str(amps_path),
        "-mode", "gridless", "-i", input_path.name,
        "-mover", run_spec["mover"],
        "-gridless-mpi-scheduler", args.scheduler,
        "-gridless-mpi-dynamic-chunk", str(args.dynamic_chunk),
        "-gridless-parallel", "THREADS",
        "-gridless-threads", str(args.nt),
        "-cutoff-search", "DIRECT_ACCESS",
    ]
    log_path = workdir / "C16_amps.log"
    run_spec["input_file"] = input_path
    run_spec["log_file"] = log_path
    run_spec["command"] = command


def serializable_run_record(run_spec):
    """Return command provenance without leaking internal Path objects to JSON."""
    return {
        "case_type": run_spec["case_type"],
        "probe": run_spec.get("probe", ""),
        "mover": run_spec["mover"],
        "field_model": run_spec["model"],
        "points": run_spec["points"],
        "rigidities_GV": run_spec["rigidities"],
        "dir_lon_res_deg": run_spec["lon_res"],
        "dir_lat_res_deg": run_spec["lat_res"],
        "limits": run_spec["limits"],
        "expected_termination_code": run_spec.get("expected_code"),
        "workdir": str(run_spec["workdir"]),
        "input_file": str(run_spec.get("input_file", "")),
        "log_file": str(run_spec.get("log_file", "")),
        "command": run_spec.get("command", []),
        "return_code": run_spec.get("return_code"),
    }


def main():
    """Build, run, audit, summarize, and return a CI-friendly exit status."""
    args = parse_args()
    if args.self_test:
        return self_test()
    apply_profile(args)
    validate_args(args)

    models = parse_name_list(args.models, "field model", {"DIPOLE", "T05"})
    movers = parse_name_list(
        args.movers, "mover", {"BORIS", "RK2", "RK4", "RK6", "GC2", "GC4", "GC6"})
    budgets = sorted(set(parse_float_list(args.budgets, "--budgets")))
    if any(value <= 0.0 for value in budgets):
        raise SystemExit("all --budgets values must be positive")
    if len(budgets) < 2:
        raise SystemExit("C16 requires at least two increasing baseline budgets")
    args.rigidity_list = parse_float_list(args.rigidities, "--rigidities")
    args.probe_rigidity_list = parse_float_list(
        args.probe_rigidities, "--probe-rigidities")
    if any(value <= 0.0 for value in args.rigidity_list + args.probe_rigidity_list):
        raise SystemExit("all rigidity values must be positive")
    lons = parse_float_list(args.lons, "--lons")
    lats = parse_float_list(args.lats, "--lats")
    if any(abs(lat) > 90.0 for lat in lats):
        raise SystemExit("observation latitudes must be within [-90,90]")

    points = build_points(lons, lats, args.alt)
    # A single equatorial probe point makes the intended resource limit dominate
    # and keeps the three diagnostic runs inexpensive.
    probe_points = build_points([0.0], [0.0], args.alt)

    launch_dir = Path.cwd().resolve()
    script_dir = Path(__file__).resolve().parent
    template_path = script_dir / "AMPS_PARAM_C16_gridless.in"
    driver_path = script_dir / "data" / "ts05_driver_C16.txt"
    base_workdir = (launch_dir / args.workdir).resolve()
    amps_path = Path(args.amps)
    if not amps_path.is_absolute():
        amps_path = (launch_dir / amps_path).resolve()

    validate_template_contract(template_path)
    validate_reference_files(script_dir)
    driver_metadata = (validate_t05_driver(driver_path)
                       if "T05" in models else None)

    if args.skip_run:
        if not base_workdir.exists():
            raise SystemExit("--skip-run workdir does not exist: %s" % base_workdir)
    else:
        if base_workdir.exists() and not args.keep:
            shutil.rmtree(base_workdir)
        base_workdir.mkdir(parents=True, exist_ok=True)

    run_specs = make_run_specs(base_workdir, movers, models, budgets, points,
                               probe_points, args)
    if not args.skip_run:
        for run_spec in run_specs:
            prepare_run(run_spec, template_path, driver_path, amps_path, args)
    else:
        # Recreate the deterministic command/input paths for provenance without
        # overwriting an existing analyzed run.
        for run_spec in run_specs:
            run_spec["input_file"] = run_spec["workdir"] / "AMPS_PARAM_C16.in"
            run_spec["log_file"] = run_spec["workdir"] / "C16_amps.log"
            run_spec["command"] = []

    # Preserve the immutable reference solution alongside every result set.
    base_workdir.mkdir(parents=True, exist_ok=True)
    if not args.skip_run:
        shutil.copy2(script_dir / "reference_C16_termination_contract.csv",
                     base_workdir / "reference_C16_termination_contract.csv")
        shutil.copy2(script_dir / "reference_C16_probe_expectations.csv",
                     base_workdir / "reference_C16_probe_expectations.csv")

    if not args.skip_run:
        for run_spec in run_specs:
            print("[%s %s %s] %s" %
                  (run_spec["case_type"], run_spec["mover"],
                   run_spec.get("probe") or run_spec["model"],
                   " ".join(run_spec["command"])))
            if args.dry_run:
                continue
            print("  workdir: %s" % run_spec["workdir"])
            run_spec["return_code"] = run_command(
                run_spec["command"], run_spec["workdir"], run_spec["log_file"])

    # Dry runs are command-generation checks, so they deliberately stop before
    # requiring numerical products from AMPS.
    if args.dry_run:
        result = {
            "test_id": TEST_ID,
            "test_name": TEST_NAME,
            "runner_schema_version": RUNNER_SCHEMA_VERSION,
            "runner_release": RUNNER_RELEASE,
            "profile": args.profile,
            "dry_run": True,
            "passed": None,
            "t05_driver": driver_metadata,
            "run_records": [serializable_run_record(spec) for spec in run_specs],
        }
        (base_workdir / "C16_commands.json").write_text(
            json.dumps(result["run_records"], indent=2))
        (base_workdir / "C16_result.json").write_text(json.dumps(result, indent=2))
        print("\nC16 dry run complete: %d inputs prepared under %s" %
              (len(run_specs), base_workdir))
        return 0

    overall_passed = True
    messages = []
    audits = {}
    summary_rows = []
    sample_rows = []
    count_rows = []
    probe_rows = []

    for run_spec in run_specs:
        if run_spec.get("return_code") not in (None, 0):
            overall_passed = False
            messages.append("AMPS failed (%s/%s/%s), exit code %d; see %s" %
                            (run_spec["mover"], run_spec["model"],
                             run_spec.get("probe") or budget_label(
                                 run_spec["limits"]["max_trace_time"]),
                             run_spec["return_code"], run_spec["log_file"]))
            continue
        expected_directions = expected_direction_keys(
            run_spec["lon_res"], run_spec["lat_res"])
        audit = audit_run(
            run_spec, run_spec["points"], expected_directions,
            run_spec["rigidities"], args.limit_diagnostic_slack)
        audits[id(run_spec)] = audit
        summary_rows.extend(audit["summaries"])
        sample_rows.extend(audit["sample_rows"])
        count_rows.extend(termination_count_rows(run_spec, audit))
        if not audit["passed"]:
            overall_passed = False
            for message in audit["messages"]:
                messages.append("%s/%s/%s: %s" %
                                (run_spec["mover"], run_spec["model"],
                                 run_spec.get("probe") or budget_label(
                                     run_spec["limits"]["max_trace_time"]),
                                 message))
        if run_spec["case_type"] == "probe":
            row = probe_result(run_spec, audit, args.probe_target_fraction)
            probe_rows.append(row)
            if not row["passed"]:
                overall_passed = False
                messages.append(
                    "%s probe %s produced code %d for %d/%d samples (%.3f; required %.3f)" %
                    (run_spec["mover"], run_spec["probe"],
                     run_spec["expected_code"], row["target_code_samples"],
                     row["samples"], row["target_fraction"],
                     row["minimum_target_fraction"]))

    # Compare only successful baseline audits, in adjacent budget pairs.  Exact
    # key equality is checked inside compare_budget_pair and cannot be hidden by
    # aggregate unresolved counts.
    convergence_rows = []
    sample_convergence_rows = []
    for mover in movers:
        for model in models:
            series = [spec for spec in run_specs
                      if spec["case_type"] == "baseline" and
                      spec["mover"] == mover and spec["model"] == model]
            series.sort(key=lambda spec: spec["limits"]["max_trace_time"])
            for short_spec, long_spec in zip(series[:-1], series[1:]):
                short_audit = audits.get(id(short_spec))
                long_audit = audits.get(id(long_spec))
                if short_audit is None or long_audit is None:
                    overall_passed = False
                    messages.append("cannot compare %s/%s budgets %g -> %g" %
                                    (mover, model,
                                     short_spec["limits"]["max_trace_time"],
                                     long_spec["limits"]["max_trace_time"]))
                    continue
                conv, rows = compare_budget_pair(
                    short_spec, long_spec, short_audit, long_audit)
                convergence_rows.append(conv)
                sample_convergence_rows.extend(rows)
                if not conv["passed"]:
                    overall_passed = False
                    messages.append(
                        "%s/%s budget %g -> %g has %d resolved-state regressions" %
                        (mover, model, conv["short_budget_s"],
                         conv["long_budget_s"], conv["regressions"]))

            # The final budget must leave a bounded unresolved population.  This
            # is the only quantitative model-level gate; all earlier budgets are
            # convergence instruments and may legitimately be less complete.
            if series:
                final_spec = series[-1]
                final_audit = audits.get(id(final_spec))
                if final_audit is not None:
                    total = final_audit["counts"]["samples"]
                    unresolved = final_audit["counts"]["state_2"]
                    fraction = float(unresolved) / total if total else 1.0
                    if fraction > args.max_final_unresolved_fraction:
                        overall_passed = False
                        messages.append(
                            "%s/%s final budget %g s unresolved=%d/%d (%.3f), ceiling %.3f" %
                            (mover, model,
                             final_spec["limits"]["max_trace_time"],
                             unresolved, total, fraction,
                             args.max_final_unresolved_fraction))

    write_dict_csv(summary_rows, base_workdir / "C16_summary.csv")
    write_dict_csv(count_rows, base_workdir / "C16_termination_counts.csv")
    write_dict_csv(sample_rows, base_workdir / "C16_sample_audit.csv")
    write_dict_csv(convergence_rows, base_workdir / "C16_convergence.csv")
    write_dict_csv(sample_convergence_rows,
                   base_workdir / "C16_sample_convergence.csv")
    write_dict_csv(probe_rows, base_workdir / "C16_probe_summary.csv")
    commands = [serializable_run_record(spec) for spec in run_specs]
    (base_workdir / "C16_commands.json").write_text(json.dumps(commands, indent=2))

    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "runner_release": RUNNER_RELEASE,
        "profile": args.profile,
        "passed": overall_passed,
        "dry_run": False,
        "models": models,
        "movers": movers,
        "trace_budgets_s": budgets,
        "points": points,
        "rigidities_GV": args.rigidity_list,
        "direction_resolution_deg": {
            "longitude": args.dir_lon_res, "latitude": args.dir_lat_res},
        "max_final_unresolved_fraction": args.max_final_unresolved_fraction,
        "probe_target_fraction": args.probe_target_fraction,
        "t05_driver": driver_metadata,
        "termination_contract": {
            str(code): {
                "name": value[0], "expected_access_state": value[1],
                "category": value[2], "policy": value[3]}
            for code, value in sorted(TERMINATION_CONTRACT.items())
        },
        "run_records": commands,
        "summary": summary_rows,
        "convergence": convergence_rows,
        "probe_summary": probe_rows,
        "messages": messages,
    }
    (base_workdir / "C16_result.json").write_text(json.dumps(result, indent=2))

    print("\nC16 results written to:")
    for name in ("C16_summary.csv", "C16_termination_counts.csv",
                 "C16_sample_audit.csv", "C16_convergence.csv",
                 "C16_sample_convergence.csv", "C16_probe_summary.csv",
                 "C16_commands.json", "C16_result.json"):
        print("  %s" % (base_workdir / name))

    if overall_passed:
        print("\nC16 PASS: termination closure and budget monotonicity satisfied.")
        return 0
    print("\nC16 FAIL:")
    for message in messages:
        print("  - %s" % message)
    return 1


if __name__ == "__main__":
    sys.exit(main())

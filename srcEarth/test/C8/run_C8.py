#!/usr/bin/env python3
"""C8 — realistic-field directional access and East-West validation.

The runner executes the same GRIDLESS/DIRECT_ACCESS calculation for T96 and
T05, with an equal-mass positive and negative particle in each field.  It then
audits the complete three-state access cubes before reducing them to an
AMS-like local zenith field of view.

C8 intentionally distinguishes observables that the current AMPS file format
contains from observables that require a future producer extension.  Access
state, termination reason, high-rigidity transparency, model agreement, and the
charge-dependent East-West response are production acceptance gates.  The
published T96/T05 asymptotic-direction benchmark becomes a gate only when the
output contains explicit asymptotic/exit-direction columns and the user selects
``--asymptotic-policy REQUIRE``.

Only the Python standard library is required.  ``--self-test`` validates input
rendering, reference schemas, geometry, parser strictness, and the principal
acceptance reductions without launching AMPS.

Two compatibility rules are deliberately enforced here rather than left to
downstream CSV reductions:

* the DIRECT_ACCESS producer must retain its 24-column Step-4 prefix (enforced
  by UDirectionalAccess), while this recovery reader addresses values by their
  VARIABLES names.  It can therefore audit both the corrected append-only
  product and already-generated files from the brief historical Step-5 layout
  that inserted two named weight fields into the prefix;
* AMPS emits every longitude at both geographic poles.  Those rows are all
  required, checked for azimuthal state degeneracy (C8-G11), and collapsed to
  one exactly weighted polar-cap sample only after that check succeeds.

These are validation safeguards, not convenience fallbacks: malformed schemas,
missing rows, and disagreeing polar duplicates remain hard failures.
"""

from __future__ import print_function

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


TEST_ID = "C8"
TEST_NAME = "Realistic-field directional access and East-West validation"
RUNNER_SCHEMA_VERSION = 2
RUNNER_RELEASE = "2026-09-23"

RE_KM = 6371.2
DEFAULT_EPOCH = "2012-05-17T06:00:00"
DEFAULT_RIGIDITIES_GV = (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0)
DEFAULT_LONS_DEG = (0.0, 180.0)
DEFAULT_LATS_DEG = (-51.6, 0.0, 51.6)

# AMPS's DIRECT_ACCESS producer uses this closed three-state contract.  Codes
# 6--8 are kept as UNRESOLVED in the data product, but C8 additionally treats
# them as fatal numerical outcomes so a broken mover cannot pass by reporting
# a large population of nominally "unknown" samples.
TERMINATION_CONTRACT = {
    0: ("OUTER_BOUNDARY_ALLOWED", 1, False),
    1: ("INNER_BOUNDARY_FORBIDDEN", 0, False),
    2: ("MAGNETICALLY_TRAPPED_FORBIDDEN", 0, False),
    3: ("TIME_LIMIT", 2, False),
    4: ("STEP_LIMIT", 2, False),
    5: ("DISTANCE_LIMIT", 2, False),
    6: ("INVALID_TIME_STEP", 2, True),
    7: ("INVALID_FIELD", 2, True),
    8: ("NUMERICAL_FAILURE", 2, True),
    9: ("DRIFT_TRAPPED_FORBIDDEN", 0, False),
}
STATE_NAMES = {0: "FORBIDDEN", 1: "ALLOWED", 2: "UNRESOLVED"}

# Stable Step-4 DIRECT_ACCESS prefix.  The producer intentionally keeps these
# names and this order unchanged, while Step 5 appends another 21 columns.  C8
# maps row values through the VARIABLES declaration, so an extension cannot
# shift ``access_state`` (the failure mode that previously left all comparison
# products empty).  Requiring the complete prefix also makes an old/incomplete
# executable fail explicitly rather than silently validating a partial cube.
DIRECT_ACCESS_CORE_COLUMNS = (
    "lon_deg", "lat_deg", "rigidity_gv", "energy_mev",
    "access_state", "allowed", "unresolved", "termination_code",
    "trace_time_s", "trace_distance_re", "trace_steps", "retry_count",
    "primary_termination_code", "primary_trace_time_s",
    "trace_extension_count", "initial_trace_limit_s", "final_trace_limit_s",
    "mirror_points", "bounce_cycles", "drift_revolutions", "drift_angle_deg",
    "drift_mean_radius_change_re", "trap_mechanism",
    "momentum_relative_spread",
)

# Current Step-5 suffix, used by the self-test to reproduce the real 45-column
# product.  Runtime parsing intentionally permits additional uniquely named
# suffix fields so future append-only diagnostics do not break C8.
DIRECT_ACCESS_STEP5_COLUMNS = (
    "direction_weight_sr", "weighted_access_sr", "exit_state_valid",
    "x_exit_m", "y_exit_m", "z_exit_m", "px_exit_si", "py_exit_si",
    "pz_exit_si", "vx_exit_unit", "vy_exit_unit", "vz_exit_unit",
    "cos_alpha_exit", "trace_time_at_exit_s", "rigidity_at_exit_gv",
    "adaptive_refined_intervals", "adaptive_estimated_error_gv",
    "adaptive_max_ambiguous_width_gv", "adaptive_target_reached",
    "adaptive_max_samples_reached", "response_weighted_unresolved_support",
)

# Initial Step-5 builds emitted these same 45 named values in a shifted order.
# The producer has since restored the stable prefix, but C8 keeps read-only
# support so completed multi-hour runs can be re-audited with ``--skip-run``.
DIRECT_ACCESS_HISTORICAL_STEP5_COLUMNS = (
    "lon_deg", "lat_deg", "direction_weight_sr", "rigidity_gv",
    "energy_mev", "access_state", "allowed", "weighted_access_sr",
    "unresolved", "termination_code", "trace_time_s", "trace_distance_re",
    "trace_steps", "retry_count", "primary_termination_code",
    "primary_trace_time_s", "trace_extension_count", "initial_trace_limit_s",
    "final_trace_limit_s", "mirror_points", "bounce_cycles",
    "drift_revolutions", "drift_angle_deg", "drift_mean_radius_change_re",
    "trap_mechanism", "momentum_relative_spread", "exit_state_valid",
    "x_exit_m", "y_exit_m", "z_exit_m", "px_exit_si", "py_exit_si",
    "pz_exit_si", "vx_exit_unit", "vy_exit_unit", "vz_exit_unit",
    "cos_alpha_exit", "trace_time_at_exit_s", "rigidity_at_exit_gv",
    "adaptive_refined_intervals", "adaptive_estimated_error_gv",
    "adaptive_max_ambiguous_width_gv", "adaptive_target_reached",
    "adaptive_max_samples_reached", "response_weighted_unresolved_support",
)

# Profiles change numerical coverage/cost, never the meaning of a pass.  SMOKE
# is useful for parser and installation checks; ROUTINE is the validation run;
# THOROUGH adds angular resolution and a longer trajectory budget.
PROFILE_DEFAULTS = {
    "SMOKE": {
        "dir_res_deg": 30.0,
        "max_trace_time_s": 300.0,
        "max_steps": 1000000,
    },
    "ROUTINE": {
        "dir_res_deg": 15.0,
        "max_trace_time_s": 600.0,
        "max_steps": 2000000,
    },
    "THOROUGH": {
        "dir_res_deg": 10.0,
        "max_trace_time_s": 1200.0,
        "max_steps": 4000000,
    },
}


def parse_float_list(text, label):
    """Parse a comma/semicolon separated finite-number list."""
    values = []
    for token in str(text).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = float(token)
        except ValueError:
            raise SystemExit("Invalid %s value %r" % (label, token))
        if not math.isfinite(value):
            raise SystemExit("Non-finite %s value %r" % (label, token))
        values.append(value)
    if not values:
        raise SystemExit("No values supplied for %s" % label)
    return values


def fmt_number(value):
    """Render a compact decimal without locale-dependent formatting."""
    return "%.12g" % float(value)


def bool_token(value):
    """Return the AMPS input syntax for a Boolean."""
    return "T" if value else "F"


def direction_key(lon_deg, lat_deg):
    """Canonicalize a sky-grid coordinate while preserving polar work items."""
    lon = float(lon_deg) % 360.0
    if abs(lon - 360.0) < 1.0e-8:
        lon = 0.0
    return (round(lon, 8), round(float(lat_deg), 8))


def rigidity_key(value):
    """Canonical key tight enough to expose a changed requested grid."""
    return round(float(value), 10)


def vector_from_lon_lat(lon_deg, lat_deg):
    """Return the unit Cartesian vector for a GSM longitude/latitude."""
    lon = math.radians(float(lon_deg))
    lat = math.radians(float(lat_deg))
    c = math.cos(lat)
    return (c * math.cos(lon), c * math.sin(lon), math.sin(lat))


def dot(a, b):
    """Small dependency-free three-vector dot product."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a):
    """Euclidean norm used by local-aperture geometry."""
    return math.sqrt(max(0.0, dot(a, a)))


def unit(a):
    """Normalize a vector and reject degenerate geometry early."""
    magnitude = norm(a)
    if magnitude <= 0.0:
        raise ValueError("Cannot normalize a zero vector")
    return tuple(value / magnitude for value in a)


def point_xyz_km(lon_deg, lat_deg, alt_km):
    """Construct the spherical 400-km validation point in GSM coordinates."""
    radius = RE_KM + float(alt_km)
    return tuple(radius * value for value in vector_from_lon_lat(lon_deg, lat_deg))


def build_points(lons_deg, lats_deg, alt_km):
    """Build the representative two-local-time, three-latitude point set."""
    points = []
    for lon_deg in lons_deg:
        for lat_deg in lats_deg:
            xyz = point_xyz_km(lon_deg, lat_deg, alt_km)
            points.append({
                "point_id": len(points),
                "obs_lon_deg": float(lon_deg),
                "obs_lat_deg": float(lat_deg),
                "obs_alt_km": float(alt_km),
                "x_km": xyz[0], "y_km": xyz[1], "z_km": xyz[2],
            })
    return points


def points_block(points):
    """Render POINTS with comments that make point-to-file mapping auditable."""
    lines = []
    for point in points:
        lines.append(
            "POINT %.12e %.12e %.12e ! km point_id=%d lon=%g lat=%g alt=%g" %
            (point["x_km"], point["y_km"], point["z_km"], point["point_id"],
             point["obs_lon_deg"], point["obs_lat_deg"], point["obs_alt_km"]))
    return "\n".join(lines)


def expected_direction_keys(lon_res_deg, lat_res_deg):
    """Derive the non-polar part of the regular directional grid."""
    nlon_float = 360.0 / float(lon_res_deg)
    nlat_float = 180.0 / float(lat_res_deg)
    nlon = int(round(nlon_float))
    nlat_intervals = int(round(nlat_float))
    if abs(nlon - nlon_float) > 1.0e-9 or abs(nlat_intervals - nlat_float) > 1.0e-9:
        raise ValueError("Directional resolution must divide 360 and 180 degrees")
    return {
        direction_key(i * lon_res_deg, -90.0 + j * lat_res_deg)
        for j in range(1, nlat_intervals)
        for i in range(nlon)
    }


def expected_pole_direction_keys(lon_res_deg, lat_res_deg):
    """Return every longitude-tagged polar work item emitted by AMPS.

    Longitude is geometrically degenerate at +/-90 degrees, but the producer
    schedules and writes one row for every longitude.  C8 must audit all of
    those rows before selecting a canonical sample for solid-angle folds.
    """
    # Reuse the divisibility validation in ``expected_direction_keys``.
    expected_direction_keys(lon_res_deg, lat_res_deg)
    nlon = int(round(360.0 / float(lon_res_deg)))
    return {
        direction_key(i * lon_res_deg, pole_lat)
        for pole_lat in (-90.0, 90.0)
        for i in range(nlon)
    }


def expected_full_direction_keys(lon_res_deg, lat_res_deg):
    """Return the exact DIRECT_ACCESS output grid, including polar duplicates."""
    return (expected_direction_keys(lon_res_deg, lat_res_deg) |
            expected_pole_direction_keys(lon_res_deg, lat_res_deg))


def sky_cell_weight_sr(lon_res_deg, lat_res_deg, sky_lat_deg):
    """Exact solid angle of one non-polar longitude/latitude cell."""
    lower = max(-90.0, float(sky_lat_deg) - 0.5 * lat_res_deg)
    upper = min(90.0, float(sky_lat_deg) + 0.5 * lat_res_deg)
    return math.radians(lon_res_deg) * (
        math.sin(math.radians(upper)) - math.sin(math.radians(lower)))


def pole_cap_weight_sr(lat_res_deg):
    """Exact solid angle assigned to one canonical polar sample."""
    half_angle = math.radians(0.5 * float(lat_res_deg))
    return 2.0 * math.pi * (1.0 - math.cos(half_angle))


def fov_direction_terms(expected_directions, lon_res_deg, lat_res_deg):
    """Return canonical direction/weight pairs for all-sky reductions.

    ``expected_directions`` contains only regular, non-polar cells.  Each pole
    is added once at longitude zero after ``audit_case`` has confirmed that all
    longitude-tagged producer rows at that pole agree in access state.  This
    preserves exactly 4*pi steradians without copying a neighboring latitude.
    """
    terms = [
        (dkey, sky_cell_weight_sr(lon_res_deg, lat_res_deg, dkey[1]))
        for dkey in sorted(expected_directions)
    ]
    cap_weight = pole_cap_weight_sr(lat_res_deg)
    terms.extend(((direction_key(0.0, -90.0), cap_weight),
                  (direction_key(0.0, 90.0), cap_weight)))
    return terms


def local_aperture(point, sky_lon_deg, sky_lat_deg, fov_half_angle_deg,
                   east_west_deadband):
    """Classify one arrival direction in an AMS-like local look aperture.

    AMPS directional coordinates describe particle arrival velocity.  A
    detector's physical look vector is therefore the negative of that vector.
    Local zenith is the outward radial unit vector at the observation point.
    EAST/WEST are assigned from the horizontal look-vector projection; a small
    deadband removes directions too close to the north/south meridian to carry
    a robust East-West label.
    """
    arrival = vector_from_lon_lat(sky_lon_deg, sky_lat_deg)
    look = tuple(-value for value in arrival)
    radial = vector_from_lon_lat(point["obs_lon_deg"], point["obs_lat_deg"])
    cos_zenith = max(-1.0, min(1.0, dot(look, radial)))
    if cos_zenith < math.cos(math.radians(fov_half_angle_deg)) - 1.0e-12:
        return "OUTSIDE", math.degrees(math.acos(cos_zenith)), 0.0

    lon = math.radians(point["obs_lon_deg"])
    east = (-math.sin(lon), math.cos(lon), 0.0)
    horizontal = tuple(look[i] - cos_zenith * radial[i] for i in range(3))
    horizontal_norm = norm(horizontal)
    if horizontal_norm <= 1.0e-12:
        return "CENTER", 0.0, 0.0
    east_score = dot(tuple(value / horizontal_norm for value in horizontal), east)
    if east_score > east_west_deadband:
        sector = "EAST"
    elif east_score < -east_west_deadband:
        sector = "WEST"
    else:
        sector = "CENTER"
    return sector, math.degrees(math.acos(cos_zenith)), east_score


def field_block(model):
    """Render T96/T05 with identical epoch and solar-wind driver provenance."""
    if model not in ("T96", "T05"):
        raise ValueError("C8 supports only T96 and T05, not %s" % model)
    # Static values are parser-safe fallbacks; the driver row at EPOCH is the
    # authoritative value.  W1--W6 are used only by T05.
    return "\n".join([
        "FIELD_MODEL                    %s" % model,
        "EPOCH                          %s" % DEFAULT_EPOCH,
        "DRIVER_FILE                    ts05_driver_C8.txt",
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


def render_input(template_path, output_path, run_id, model, charge, points,
                 rigidities_gv, args):
    """Render one case deck and fail if any template token was forgotten."""
    text = template_path.read_text(errors="replace")
    replacements = {
        "__RUN_ID__": run_id,
        "__CUTOFF_NENERGY__": str(len(rigidities_gv)),
        "__DIRMAP_LON_RES__": fmt_number(args.dir_res),
        "__DIRMAP_LAT_RES__": fmt_number(args.dir_res),
        "__CUTOFF_RIGIDITY_LIST_GV__": ",".join(
            fmt_number(value) for value in rigidities_gv),
        "__MAX_TRACE_TIME__": fmt_number(args.max_trace_time),
        "__SPECIES_NAME__": "PROTON" if charge > 0 else "NEGATIVE_PROTON",
        "__CHARGE__": str(int(charge)),
        "__FIELD_BLOCK__": field_block(model),
        "__POINTS_BLOCK__": points_block(points),
        "__DT_TRACE__": fmt_number(args.dt_trace),
        "__MAX_STEPS__": str(int(args.max_steps)),
        "__GRIDLESS_MPI_SCHEDULER__": args.scheduler,
        "__GRIDLESS_MPI_DYNAMIC_CHUNK__": str(int(args.dynamic_chunk)),
        "__GRIDLESS_THREADS__": str(int(args.nt)),
    }
    for token, replacement in replacements.items():
        text = text.replace(token, replacement)
    remaining = sorted(set(re.findall(r"__[A-Z0-9_]+__", text)))
    if remaining:
        raise RuntimeError("Unresolved input tokens: %s" % ", ".join(remaining))
    output_path.write_text(text)


def run_command(command, workdir, log_path):
    """Run AMPS and mirror combined output to the console and a durable log."""
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
    """Normalize a quoted Tecplot label while retaining semantic names."""
    normalized = name.strip().lower().replace(" ", "_").replace("-", "_")
    return re.sub(r"[^a-z0-9_]+", "", normalized)


def validate_directional_access_schema(variables, path, line_number=None):
    """Require every uniquely named compatibility field in a Tecplot header.

    The compiled producer's *current* ABI is the exact 24-field prefix followed
    by append-only diagnostics, and a separate source-level regression test
    enforces it.  C8 is also a recovery postprocessor: it must be able to audit
    the named 45-column files already emitted by the initial Step-5 build, where
    two weight fields were inserted near the front.  Name-based lookup safely
    supports both layouts without weakening any scientific gate.

    Duplicate normalized names are rejected because converting such a header
    to a dictionary would otherwise discard one value silently.
    """
    where = "%s%s" % (path, (":%d" % line_number) if line_number else "")
    if not variables:
        raise RuntimeError("%s has an empty VARIABLES declaration" % where)
    duplicates = sorted(name for name, count in Counter(variables).items()
                        if count > 1)
    if duplicates:
        raise RuntimeError("%s has duplicate VARIABLES names: %s" %
                           (where, ", ".join(duplicates)))
    missing = [name for name in DIRECT_ACCESS_CORE_COLUMNS
               if name not in variables]
    if missing:
        raise RuntimeError("%s lacks required DIRECT_ACCESS columns: %s" %
                           (where, ", ".join(missing)))
    core_count = len(DIRECT_ACCESS_CORE_COLUMNS)
    current_layout = tuple(variables[:core_count]) == DIRECT_ACCESS_CORE_COLUMNS
    historical_count = len(DIRECT_ACCESS_HISTORICAL_STEP5_COLUMNS)
    historical_layout = (
        tuple(variables[:historical_count]) ==
        DIRECT_ACCESS_HISTORICAL_STEP5_COLUMNS)
    if not current_layout and not historical_layout:
        raise RuntimeError(
            "%s has a noncanonical DIRECT_ACCESS column order; expected the "
            "24-column compatibility prefix or the documented historical "
            "Step-5 45-column layout" % where)


def directional_access_schema_layout(variables):
    """Return provenance for the declared core-column ordering."""
    core_count = len(DIRECT_ACCESS_CORE_COLUMNS)
    if tuple(variables[:core_count]) == DIRECT_ACCESS_CORE_COLUMNS:
        return "legacy_prefix_append_only"
    if tuple(variables[:len(DIRECT_ACCESS_HISTORICAL_STEP5_COLUMNS)]) == \
            DIRECT_ACCESS_HISTORICAL_STEP5_COLUMNS:
        return "historical_step5_shifted"
    # validate_directional_access_schema has already rejected all other orders.
    raise RuntimeError("Unvalidated DIRECT_ACCESS schema layout")


def exact_int(record, name):
    """Parse an integer-valued numeric field without silently truncating it."""
    try:
        value = float(record[name])
    except (KeyError, ValueError, OverflowError):
        raise RuntimeError("Column %s contains a non-number %r" %
                           (name, record.get(name, "<missing>")))
    if not math.isfinite(value) or not value.is_integer():
        raise RuntimeError("Column %s is not a finite integer: %r" %
                           (name, record[name]))
    return int(value)


def optional_float(record, name):
    """Parse an optional finite column; malformed present values are errors."""
    if name not in record:
        return None
    try:
        value = float(record[name])
    except ValueError:
        raise RuntimeError("Column %s contains a non-number %r" %
                           (name, record[name]))
    if not math.isfinite(value):
        raise RuntimeError("Column %s contains a non-finite value" % name)
    return value


def asymptotic_lon_lat(record):
    """Read one of the documented candidate exit-direction schemas.

    The current baseline AMPS producer has no such columns.  Supporting both a
    spherical pair and a Cartesian unit vector lets C8 become a real published-
    benchmark test as soon as the producer exports either unambiguous schema.
    Partial schemas are rejected rather than silently treated as unavailable.
    """
    spherical_pairs = (
        ("asymptotic_lon_deg", "asymptotic_lat_deg"),
        ("exit_lon_deg", "exit_lat_deg"),
    )
    for lon_name, lat_name in spherical_pairs:
        present = (lon_name in record, lat_name in record)
        if any(present) and not all(present):
            raise RuntimeError("Incomplete asymptotic schema: %s/%s" %
                               (lon_name, lat_name))
        if all(present):
            lon = optional_float(record, lon_name) % 360.0
            lat = optional_float(record, lat_name)
            if lat < -90.0 or lat > 90.0:
                raise RuntimeError("Asymptotic latitude outside [-90,90]")
            return (lon, lat, lon_name.replace("_lon_deg", ""))

    vector_triplets = (
        ("asymptotic_x", "asymptotic_y", "asymptotic_z"),
        ("exit_direction_x", "exit_direction_y", "exit_direction_z"),
    )
    for x_name, y_name, z_name in vector_triplets:
        present = tuple(name in record for name in (x_name, y_name, z_name))
        if any(present) and not all(present):
            raise RuntimeError("Incomplete asymptotic vector schema")
        if all(present):
            direction = unit(tuple(optional_float(record, name)
                                   for name in (x_name, y_name, z_name)))
            lon = math.degrees(math.atan2(direction[1], direction[0])) % 360.0
            lat = math.degrees(math.asin(max(-1.0, min(1.0, direction[2]))))
            return (lon, lat, x_name.rsplit("_", 1)[0])
    return None


def parse_directional_access(path):
    """Parse a DIRECT_ACCESS cube with strict, append-tolerant schema checks.

    Modern files declare 45 or more VARIABLES: normally the stable 24-column
    Step-4 prefix followed by append-only diagnostics.  The header controls all
    name lookups (including recovery of historical named layouts), and every
    row must have exactly the declared width.  A strict
    24-value headerless fallback remains for archived Step-4 output only; a
    headerless extended row is ambiguous and is rejected.
    """
    variables = []
    reading_variables = False
    saw_variables = False
    schema_validated = False
    samples = {}

    with path.open("r", errors="replace") as stream:
        for line_number, raw in enumerate(stream, start=1):
            line = raw.strip()
            if not line:
                continue
            upper = line.upper()
            if upper.startswith("VARIABLES"):
                if saw_variables:
                    raise RuntimeError("%s:%d contains multiple VARIABLES blocks" %
                                       (path, line_number))
                saw_variables = True
                reading_variables = True
                variables.extend(normalize_variable_name(name)
                                 for name in re.findall(r'"([^"]+)"', line))
                continue
            if upper.startswith("ZONE"):
                if saw_variables and not schema_validated:
                    validate_directional_access_schema(
                        variables, path, line_number)
                    schema_validated = True
                reading_variables = False
                continue
            if reading_variables:
                names = re.findall(r'"([^"]+)"', line)
                if names:
                    variables.extend(normalize_variable_name(name) for name in names)
                    continue
                reading_variables = False
                validate_directional_access_schema(variables, path, line_number)
                schema_validated = True
            if upper.startswith(("TITLE", "AUXDATA", "#", "!")):
                continue
            parts = line.replace(",", " ").split()
            if not variables:
                if saw_variables:
                    # A declared-but-empty/malformed header must never fall back
                    # to positional parsing.
                    validate_directional_access_schema(
                        variables, path, line_number)
                if len(parts) != len(DIRECT_ACCESS_CORE_COLUMNS):
                    raise RuntimeError(
                        "%s:%d headerless row has width %d; only the exact "
                        "24-column legacy schema is accepted" %
                        (path, line_number, len(parts)))
                variables = list(DIRECT_ACCESS_CORE_COLUMNS)
                schema_validated = True
            elif not schema_validated:
                validate_directional_access_schema(variables, path, line_number)
                schema_validated = True
            if len(parts) != len(variables):
                raise RuntimeError("%s:%d row width %d does not match %d VARIABLES" %
                                   (path, line_number, len(parts), len(variables)))
            if any(token.strip().strip('"').strip() == "" for token in parts):
                raise RuntimeError("%s:%d contains a blank VARIABLES value" %
                                   (path, line_number))
            record = dict(zip(variables, parts))
            try:
                # Parse every core value, not just the fields used in a gate.
                # That guarantees a blank or non-numeric compatibility-prefix
                # field cannot be ignored merely because a later reduction does
                # not currently consume it.
                core_values = {
                    name: float(record[name])
                    for name in DIRECT_ACCESS_CORE_COLUMNS
                }
                lon = core_values["lon_deg"]
                lat = core_values["lat_deg"]
                rigidity = core_values["rigidity_gv"]
                state = exact_int(record, "access_state")
                code = exact_int(record, "termination_code")
                trace_time = core_values["trace_time_s"]
                trace_distance = core_values["trace_distance_re"]
                trace_steps = exact_int(record, "trace_steps")
                retry_count = exact_int(record, "retry_count")
            except (ValueError, OverflowError) as exc:
                raise RuntimeError("Malformed row %d in %s: %s" %
                                   (line_number, path, exc))
            except RuntimeError as exc:
                raise RuntimeError("Malformed row %d in %s: %s" %
                                   (line_number, path, exc))
            if not all(math.isfinite(value) for value in core_values.values()):
                raise RuntimeError("Non-finite C8 value at %s:%d" %
                                   (path, line_number))
            if not (-90.0 <= lat <= 90.0) or rigidity <= 0.0 or \
                    trace_time < 0.0 or trace_distance < 0.0 or trace_steps < 0:
                raise RuntimeError("Out-of-range C8 value at %s:%d" %
                                   (path, line_number))
            if state not in STATE_NAMES:
                raise RuntimeError("Unknown access_state=%d at %s:%d" %
                                   (state, path, line_number))
            if code not in TERMINATION_CONTRACT:
                raise RuntimeError("Unknown termination_code=%d at %s:%d" %
                                   (code, path, line_number))

            key = (direction_key(lon, lat), rigidity_key(rigidity))
            if key in samples:
                raise RuntimeError("Duplicate direction/rigidity row in %s: %s" %
                                   (path, key))
            expected_state = TERMINATION_CONTRACT[code][1]
            contract_error = state != expected_state
            fatal = TERMINATION_CONTRACT[code][2]
            samples[key] = {
                "lon_deg": lon % 360.0,
                "lat_deg": lat,
                "rigidity_GV": rigidity,
                "access_state": state,
                "termination_code": code,
                "termination_name": TERMINATION_CONTRACT[code][0],
                "trace_time_s": trace_time,
                "trace_distance_Re": trace_distance,
                "trace_steps": trace_steps,
                "retry_count": retry_count,
                "contract_error": contract_error,
                "fatal": fatal,
                "asymptotic": asymptotic_lon_lat(record),
            }
    if saw_variables and not schema_validated:
        validate_directional_access_schema(variables, path)
    if not samples:
        raise RuntimeError("No DIRECT_ACCESS samples parsed from %s" % path)
    return {"columns": tuple(variables),
            "schema_layout": directional_access_schema_layout(variables),
            "samples": samples}


def find_access_file(workdir, point_id):
    """Locate the production GRIDLESS name and compatible numbered variants."""
    preferred = workdir / ("cutoff_gridless_dir_access_point_%04d.dat" % point_id)
    if preferred.exists():
        return preferred
    matches = sorted(workdir.glob("*dir*access*point*%04d*.dat" % point_id))
    return matches[0] if matches else preferred


def match_requested_rigidity(value, requested):
    """Match only numerical roundoff, not a materially changed rigidity node."""
    for target in requested:
        if math.isclose(float(value), float(target), rel_tol=5.0e-10,
                        abs_tol=1.0e-12):
            return float(target)
    return None


def audit_case(case, points, expected_directions, expected_full_directions,
               rigidities_gv):
    """Audit the complete cube, including duplicate polar work items.

    ``expected_directions`` is the non-polar grid used by later canonical
    reductions; ``expected_full_directions`` is the producer contract audited
    by C8-G01.  All polar rows remain in ``sample_map``.  Reductions select the
    longitude-zero row only after this function proves that every longitude at
    the same pole/rigidity has an identical three-state classification (G11).
    """
    summaries = []
    sample_map = {}
    errors = []
    total = Counter()

    for point in points:
        path = find_access_file(case["workdir"], point["point_id"])
        if not path.exists():
            errors.append("%s: missing point %d output %s" %
                          (case["case_id"], point["point_id"], path))
            continue
        try:
            parsed = parse_directional_access(path)
        except Exception as exc:
            errors.append("%s point %d parse failed: %s" %
                          (case["case_id"], point["point_id"], exc))
            continue

        got_directions = {key[0] for key in parsed["samples"]}
        missing_directions = expected_full_directions - got_directions
        extra_directions = got_directions - expected_full_directions
        matched_keys = set()
        point_counts = Counter()
        extra_rigidities = 0

        for (dkey, _), sample in parsed["samples"].items():
            requested = match_requested_rigidity(sample["rigidity_GV"], rigidities_gv)
            if dkey not in expected_full_directions or requested is None:
                if requested is None:
                    extra_rigidities += 1
                continue
            out_key = (point["point_id"], dkey[0], dkey[1],
                       rigidity_key(requested))
            if out_key in sample_map:
                errors.append("%s duplicate matched sample %s" %
                              (case["case_id"], out_key))
                continue
            sample_map[out_key] = sample
            matched_keys.add((dkey, rigidity_key(requested)))
            point_counts["samples"] += 1
            point_counts["state_%d" % sample["access_state"]] += 1
            point_counts["contract_errors"] += int(sample["contract_error"])
            point_counts["fatal"] += int(sample["fatal"])

        expected_samples = len(expected_full_directions) * len(rigidities_gv)
        missing_samples = expected_samples - len(matched_keys)
        coverage_errors = (len(missing_directions) + len(extra_directions) +
                           missing_samples + extra_rigidities)

        # At either pole all longitude labels map to the same Cartesian arrival
        # direction.  Count a violation once per pole/rigidity group, not once
        # per pair, so diagnostics remain stable as angular resolution changes.
        pole_degeneracy_errors = 0
        pole_error_groups = []
        for pole_lat in (-90.0, 90.0):
            pole_keys = sorted(dkey for dkey in expected_full_directions
                               if dkey[1] == pole_lat)
            for rigidity in rigidities_gv:
                states = set()
                for dkey in pole_keys:
                    key = (point["point_id"], dkey[0], dkey[1],
                           rigidity_key(rigidity))
                    sample = sample_map.get(key)
                    if sample is not None:
                        states.add(sample["access_state"])
                if len(states) > 1:
                    pole_degeneracy_errors += 1
                    pole_error_groups.append("lat=%g R=%g states=%s" %
                                             (pole_lat, rigidity,
                                              sorted(states)))
        summary = {
            "case_id": case["case_id"], "field_model": case["model"],
            "charge": case["charge"], "point_id": point["point_id"],
            "schema_layout": parsed["schema_layout"],
            "obs_lon_deg": point["obs_lon_deg"],
            "obs_lat_deg": point["obs_lat_deg"],
            "expected_directions": len(expected_full_directions),
            "actual_directions": len(got_directions),
            "missing_directions": len(missing_directions),
            "extra_directions": len(extra_directions),
            "expected_samples": expected_samples,
            "actual_samples": point_counts["samples"],
            "missing_samples": missing_samples,
            "extra_rigidity_rows": extra_rigidities,
            "allowed": point_counts["state_1"],
            "forbidden": point_counts["state_0"],
            "unresolved": point_counts["state_2"],
            "unresolved_fraction": (float(point_counts["state_2"]) /
                                    point_counts["samples"]
                                    if point_counts["samples"] else 1.0),
            "contract_errors": point_counts["contract_errors"],
            "fatal_terminations": point_counts["fatal"],
            "coverage_errors": coverage_errors,
            "pole_degeneracy_errors": pole_degeneracy_errors,
            "passed": int(coverage_errors == 0 and
                          point_counts["contract_errors"] == 0 and
                          point_counts["fatal"] == 0 and
                          pole_degeneracy_errors == 0),
        }
        summaries.append(summary)
        total.update(point_counts)
        total["coverage_errors"] += coverage_errors
        total["pole_degeneracy_errors"] += pole_degeneracy_errors
        if coverage_errors:
            errors.append(
                "%s point %d coverage failed: missing_directions=%d "
                "extra_directions=%d missing_samples=%d "
                "extra_rigidity_rows=%d" %
                (case["case_id"], point["point_id"],
                 len(missing_directions), len(extra_directions),
                 missing_samples, extra_rigidities))
        if pole_degeneracy_errors:
            errors.append(
                "%s point %d pole azimuthal degeneracy failed in %d group(s): %s" %
                (case["case_id"], point["point_id"],
                 pole_degeneracy_errors, "; ".join(pole_error_groups)))

    if len(summaries) != len(points):
        total["coverage_errors"] += len(points) - len(summaries)
    return summaries, sample_map, total, errors


def case_is_structurally_auditable(totals, audit_messages, sample_map,
                                   expected_sample_count):
    """Gate every reduction on a complete, internally consistent cube.

    This predicate intentionally excludes the unresolved-population science
    gate, which is evaluated separately.  Its job is to ensure reducers never
    receive partial data and consequently never manufacture blank fractions.
    """
    return (
        not audit_messages and totals["coverage_errors"] == 0 and
        totals["contract_errors"] == 0 and totals["fatal"] == 0 and
        totals["pole_degeneracy_errors"] == 0 and
        totals["samples"] == expected_sample_count and
        len(sample_map) == expected_sample_count)


def reduce_fov(case, points, sample_map, expected_directions, rigidities_gv,
               args):
    """Reduce each cube to full, EAST, WEST, and meridian FOV solid angles."""
    rows = []
    errors = []
    direction_terms = fov_direction_terms(
        expected_directions, args.dir_res, args.dir_res)
    for point in points:
        geometry = {}
        for dkey, direction_weight in direction_terms:
            sector, zenith_deg, east_score = local_aperture(
                point, dkey[0], dkey[1], args.fov_half_angle,
                args.east_west_deadband)
            if sector == "OUTSIDE":
                continue
            geometry[dkey] = {
                "sector": sector, "zenith_deg": zenith_deg,
                "east_score": east_score,
                "weight_sr": direction_weight,
            }
        lobe_geometry = Counter()
        for item in geometry.values():
            lobe_geometry["FOV"] += item["weight_sr"]
            lobe_geometry[item["sector"]] += item["weight_sr"]
        if lobe_geometry["FOV"] <= 0.0 or lobe_geometry["EAST"] <= 0.0 or \
                lobe_geometry["WEST"] <= 0.0:
            errors.append("%s point %d has empty FOV/EAST/WEST geometry" %
                          (case["case_id"], point["point_id"]))
            continue

        for rigidity in rigidities_gv:
            accum = defaultdict(Counter)
            for dkey, item in geometry.items():
                key = (point["point_id"], dkey[0], dkey[1],
                       rigidity_key(rigidity))
                sample = sample_map.get(key)
                if sample is None:
                    continue
                for lobe in ("FOV", item["sector"]):
                    accum[lobe]["total"] += item["weight_sr"]
                    accum[lobe]["state_%d" % sample["access_state"]] += \
                        item["weight_sr"]

            row = {
                "case_id": case["case_id"], "field_model": case["model"],
                "charge": case["charge"], "point_id": point["point_id"],
                "obs_lon_deg": point["obs_lon_deg"],
                "obs_lat_deg": point["obs_lat_deg"],
                "rigidity_GV": rigidity,
            }
            for lobe in ("FOV", "EAST", "WEST", "CENTER"):
                total_weight = accum[lobe]["total"]
                row[lobe.lower() + "_solid_angle_sr"] = total_weight
                for state, name in ((1, "allowed"), (0, "forbidden"),
                                    (2, "unresolved")):
                    value = accum[lobe]["state_%d" % state]
                    row[lobe.lower() + "_" + name + "_sr"] = value
                    row[lobe.lower() + "_" + name + "_fraction"] = (
                        value / total_weight if total_weight > 0.0 else "")
            rows.append(row)
    return rows, errors


def model_comparison(cases, points, expected_directions, rigidities_gv, args):
    """Compare T96/T05 access states on identical resolved FOV samples."""
    rows = []
    messages = []
    passed = True
    direction_terms = fov_direction_terms(
        expected_directions, args.dir_res, args.dir_res)
    for charge in (+1, -1):
        left = cases[("T96", charge)]["sample_map"]
        right = cases[("T05", charge)]["sample_map"]
        for point in points:
            for rigidity in rigidities_gv:
                counts = Counter()
                for dkey, direction_weight in direction_terms:
                    sector, _, _ = local_aperture(
                        point, dkey[0], dkey[1], args.fov_half_angle,
                        args.east_west_deadband)
                    if sector == "OUTSIDE":
                        continue
                    key = (point["point_id"], dkey[0], dkey[1],
                           rigidity_key(rigidity))
                    a = left.get(key)
                    b = right.get(key)
                    if a is None or b is None:
                        counts["missing"] += 1
                        counts["missing_weight_sr"] += direction_weight
                    elif a["access_state"] == 2 and b["access_state"] == 2:
                        counts["both_unresolved"] += 1
                        counts["both_unresolved_weight_sr"] += direction_weight
                    elif a["access_state"] == 2 or b["access_state"] == 2:
                        counts["one_unresolved"] += 1
                        counts["one_unresolved_weight_sr"] += direction_weight
                    else:
                        counts["resolved"] += 1
                        counts["resolved_weight_sr"] += direction_weight
                        mismatch = int(a["access_state"] != b["access_state"])
                        counts["mismatch"] += mismatch
                        counts["mismatch_weight_sr"] += mismatch * direction_weight
                mismatch_fraction = (
                    counts["mismatch_weight_sr"] /
                    counts["resolved_weight_sr"]
                    if counts["resolved_weight_sr"] else 1.0)
                row_passed = (counts["missing"] == 0 and
                              (rigidity < args.high_rigidity_min or
                               mismatch_fraction <= args.max_model_mismatch))
                rows.append({
                    "charge": charge, "point_id": point["point_id"],
                    "obs_lon_deg": point["obs_lon_deg"],
                    "obs_lat_deg": point["obs_lat_deg"],
                    "rigidity_GV": rigidity, "resolved_pairs": counts["resolved"],
                    "state_mismatches": counts["mismatch"],
                    "resolved_weight_sr": counts["resolved_weight_sr"],
                    "mismatch_weight_sr": counts["mismatch_weight_sr"],
                    "state_mismatch_fraction": mismatch_fraction,
                    "one_sided_unresolved": counts["one_unresolved"],
                    "both_unresolved": counts["both_unresolved"],
                    "missing_samples": counts["missing"],
                    "high_rigidity_gate_passed": int(row_passed),
                })

        # The published 9.5% value came from a much larger orbit/FOV sample.  C8
        # uses 15% as a conservative aggregate bound for its representative grid.
        aggregate = Counter()
        for row in rows:
            if row["charge"] != charge:
                continue
            aggregate["resolved"] += row["resolved_pairs"]
            aggregate["mismatch"] += row["state_mismatches"]
            aggregate["resolved_weight_sr"] += row["resolved_weight_sr"]
            aggregate["mismatch_weight_sr"] += row["mismatch_weight_sr"]
            aggregate["missing"] += row["missing_samples"]
        fraction = (aggregate["mismatch_weight_sr"] /
                    aggregate["resolved_weight_sr"]
                    if aggregate["resolved_weight_sr"] else 1.0)
        if fraction > args.max_model_mismatch or aggregate["missing"]:
            passed = False
            messages.append(
                "T96/T05 q=%+d aggregate state mismatch %.3f exceeds %.3f "
                "or has %d missing samples" %
                (charge, fraction, args.max_model_mismatch, aggregate["missing"]))
    return rows, passed, messages


def east_west_comparison(fov_rows, args):
    """Test the charge-odd component of the integrated East-West response.

    The statistic D(q)=T_WEST(q)-T_EAST(q) is computed after averaging the six
    observation points.  A true geomagnetic East-West signature is charge odd:
    D(+)>0 and D(-)<0.  Only transition nodes with partial total transmission
    and sufficiently resolved EAST/WEST lobes are eligible, preventing fully
    blocked or fully transparent samples from manufacturing a zero-effect pass.
    """
    grouped = defaultdict(list)
    for row in fov_rows:
        grouped[(row["field_model"], row["charge"],
                 rigidity_key(row["rigidity_GV"]))].append(row)
    output = []
    model_pass = {"T96": False, "T05": False}

    for model in ("T96", "T05"):
        for rigidity in sorted({key[2] for key in grouped if key[0] == model}):
            charge_values = {}
            charge_informative = {}
            for charge in (+1, -1):
                rows = grouped.get((model, charge, rigidity), [])
                east = []
                west = []
                overall = []
                max_lobe_unresolved = 1.0
                if rows:
                    max_lobe_unresolved = 0.0
                for row in rows:
                    east.append(float(row["east_allowed_fraction"]))
                    west.append(float(row["west_allowed_fraction"]))
                    overall.append(float(row["fov_allowed_fraction"]))
                    max_lobe_unresolved = max(
                        max_lobe_unresolved,
                        float(row["east_unresolved_fraction"]),
                        float(row["west_unresolved_fraction"]))
                east_mean = sum(east) / len(east) if east else float("nan")
                west_mean = sum(west) / len(west) if west else float("nan")
                total_mean = sum(overall) / len(overall) if overall else float("nan")
                charge_values[charge] = west_mean - east_mean
                charge_informative[charge] = (
                    math.isfinite(total_mean) and
                    args.transition_min <= total_mean <= args.transition_max and
                    max_lobe_unresolved <= args.max_lobe_unresolved)

            dplus = charge_values.get(+1, float("nan"))
            dminus = charge_values.get(-1, float("nan"))
            informative = (charge_informative.get(+1, False) and
                           charge_informative.get(-1, False))
            sign_ok = (math.isfinite(dplus) and math.isfinite(dminus) and
                       dplus >= args.min_east_west_effect and
                       dminus <= -args.min_east_west_effect)
            charge_odd_amplitude = 0.5 * (dplus - dminus)
            row_passed = informative and sign_ok
            model_pass[model] = model_pass[model] or row_passed
            output.append({
                "field_model": model, "rigidity_GV": rigidity,
                "west_minus_east_q_plus": dplus,
                "west_minus_east_q_minus": dminus,
                "charge_odd_amplitude": charge_odd_amplitude,
                "informative_transition_node": int(informative),
                "opposite_sign_effect_passed": int(row_passed),
            })

    messages = []
    for model, value in model_pass.items():
        if not value:
            messages.append(
                "%s has no informative rigidity with D(+)>=%g and D(-)<=-%g" %
                (model, args.min_east_west_effect,
                 args.min_east_west_effect))
    return output, all(model_pass.values()), messages


def asymptotic_comparison(cases, points, expected_directions, rigidities_gv,
                          args):
    """Compare optional T96/T05 asymptotic coordinates and compute RMS values."""
    rows = []
    bins = defaultdict(list)
    available_samples = 0
    direction_terms = fov_direction_terms(
        expected_directions, args.dir_res, args.dir_res)
    for charge in (+1, -1):
        left = cases[("T96", charge)]["sample_map"]
        right = cases[("T05", charge)]["sample_map"]
        for point in points:
            for dkey, direction_weight in direction_terms:
                sector, _, _ = local_aperture(
                    point, dkey[0], dkey[1], args.fov_half_angle,
                    args.east_west_deadband)
                if sector == "OUTSIDE":
                    continue
                for rigidity in rigidities_gv:
                    key = (point["point_id"], dkey[0], dkey[1],
                           rigidity_key(rigidity))
                    a = left.get(key)
                    b = right.get(key)
                    if not a or not b or a["access_state"] != 1 or \
                            b["access_state"] != 1 or not a["asymptotic"] or \
                            not b["asymptotic"]:
                        continue
                    available_samples += 1
                    lon_delta = ((a["asymptotic"][0] - b["asymptotic"][0] +
                                  180.0) % 360.0) - 180.0
                    lat_delta = a["asymptotic"][1] - b["asymptotic"][1]
                    lon_mrad = math.radians(lon_delta) * 1000.0
                    lat_mrad = math.radians(lat_delta) * 1000.0
                    bin_name = ("gt_50_GV" if rigidity >= 50.0 else
                                "20_to_30_GV" if 20.0 <= rigidity <= 30.0 else
                                "other")
                    bins[bin_name].append(
                        (lat_mrad, lon_mrad, direction_weight))
                    rows.append({
                        "charge": charge, "point_id": point["point_id"],
                        "sky_lon_deg": dkey[0], "sky_lat_deg": dkey[1],
                        "rigidity_GV": rigidity, "rigidity_bin": bin_name,
                        "latitude_difference_mrad": lat_mrad,
                        "longitude_difference_mrad": lon_mrad,
                        "direction_weight_sr": direction_weight,
                        "t96_schema": a["asymptotic"][2],
                        "t05_schema": b["asymptotic"][2],
                    })

    summary = []
    for bin_name in ("20_to_30_GV", "gt_50_GV"):
        values = bins.get(bin_name, [])
        total_weight = sum(item[2] for item in values)
        lat_rms = (math.sqrt(sum(item[2] * item[0] ** 2 for item in values) /
                             total_weight) if total_weight else None)
        lon_rms = (math.sqrt(sum(item[2] * item[1] ** 2 for item in values) /
                             total_weight) if total_weight else None)
        summary.append({
            "rigidity_bin": bin_name, "samples": len(values),
            "solid_angle_weight_sr": total_weight,
            "latitude_rms_mrad": lat_rms,
            "longitude_rms_mrad": lon_rms,
        })

    available = available_samples > 0
    passed = True
    messages = []
    high = next(item for item in summary if item["rigidity_bin"] == "gt_50_GV")
    if args.asymptotic_policy == "REQUIRE":
        if not available or not high["samples"]:
            passed = False
            messages.append(
                "asymptotic-policy REQUIRE but no common T96/T05 exit directions "
                "were exported")
        elif (high["latitude_rms_mrad"] > args.max_asymptotic_rms_mrad or
              high["longitude_rms_mrad"] > args.max_asymptotic_rms_mrad):
            passed = False
            messages.append(
                "high-rigidity asymptotic RMS exceeds %.3g mrad: lat=%.3g lon=%.3g" %
                (args.max_asymptotic_rms_mrad,
                 high["latitude_rms_mrad"], high["longitude_rms_mrad"]))
    return rows, summary, available, passed, messages


def write_csv(path, rows, fieldnames=None):
    """Write a stable CSV even when a conditional diagnostic has no rows."""
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as stream:
        if not fieldnames:
            stream.write("status\nNO_ROWS\n")
            return
        writer = csv.DictWriter(stream, fieldnames=fieldnames,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_reference_files(script_dir):
    """Validate checked-in scientific and acceptance references as data."""
    published = script_dir / "reference_C8_boschini2013.csv"
    contract = script_dir / "reference_C8_acceptance_contract.csv"
    expected = script_dir / "reference_C8_expected_physics.csv"
    if not published.exists() or not contract.exists() or not expected.exists():
        raise RuntimeError("C8 reference files are missing")
    with published.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 9:
        raise RuntimeError("Published reference must contain exactly 9 rows")
    high_all = [row for row in rows
                if row["comparison"] == "T96_vs_T05" and
                row["rigidity_bin"] == "gt_50_GV" and
                row["pressure_subset"] == "all"]
    if len(high_all) != 1 or float(high_all[0]["latitude_rms_mrad"]) != 1.5 or \
            float(high_all[0]["longitude_rms_mrad"]) != 2.3:
        raise RuntimeError("Published high-rigidity benchmark drifted")
    classification = [row for row in rows
                      if row["comparison"] == "T96_vs_T05" and
                      row["rigidity_bin"] == "0.3_to_200_GV"]
    if len(classification) != 1 or \
            float(classification[0]["classification_difference_percent"]) != 9.5:
        raise RuntimeError("Published classification benchmark drifted")
    with contract.open(newline="") as stream:
        gates = list(csv.DictReader(stream))
    if len(gates) != 11 or {row["gate_id"] for row in gates} != \
            {"C8-G%02d" % index for index in range(1, 12)}:
        raise RuntimeError("Acceptance contract must contain gates C8-G01..G11")
    with expected.open(newline="") as stream:
        expected_rows = list(csv.DictReader(stream))
    if len(expected_rows) != 6 or {row["reference_id"] for row in expected_rows} != \
            {"C8-R%02d" % index for index in range(1, 7)}:
        raise RuntimeError("Expected-physics reference must contain C8-R01..R06")
    return rows, gates, expected_rows


def synthetic_access_record(dkey, rigidity, state):
    """Return one fully numeric 45-field DIRECT_ACCESS self-test record."""
    code = {0: 1, 1: 0, 2: 3}[int(state)]
    weight = 0.25
    record = {
        "lon_deg": dkey[0], "lat_deg": dkey[1],
        "rigidity_gv": rigidity, "energy_mev": rigidity * 1000.0,
        "access_state": int(state), "allowed": int(state == 1),
        "unresolved": int(state == 2), "termination_code": code,
        "trace_time_s": 1.0, "trace_distance_re": 2.0,
        "trace_steps": 10, "retry_count": 0,
        "primary_termination_code": code, "primary_trace_time_s": 1.0,
        "trace_extension_count": 0, "initial_trace_limit_s": 300.0,
        "final_trace_limit_s": 300.0, "mirror_points": 0,
        "bounce_cycles": 0, "drift_revolutions": 0,
        "drift_angle_deg": 0.0, "drift_mean_radius_change_re": 0.0,
        "trap_mechanism": 0, "momentum_relative_spread": 0.0,
        "direction_weight_sr": weight,
        "weighted_access_sr": weight if state == 1 else 0.0,
        "exit_state_valid": 0, "x_exit_m": 0.0, "y_exit_m": 0.0,
        "z_exit_m": 0.0, "px_exit_si": 0.0, "py_exit_si": 0.0,
        "pz_exit_si": 0.0, "vx_exit_unit": 0.0, "vy_exit_unit": 0.0,
        "vz_exit_unit": 0.0, "cos_alpha_exit": 0.0,
        "trace_time_at_exit_s": 0.0, "rigidity_at_exit_gv": rigidity,
        "adaptive_refined_intervals": 0,
        "adaptive_estimated_error_gv": 0.0,
        "adaptive_max_ambiguous_width_gv": 0.0,
        "adaptive_target_reached": 1, "adaptive_max_samples_reached": 0,
        "response_weighted_unresolved_support": 0.0,
    }
    return record


def write_synthetic_access(path, directions, rigidities, state_function,
                           columns=None, with_asymptotic=False,
                           include_header=True):
    """Create a valid named-schema cube used only by ``--self-test``."""
    if columns is None:
        columns = DIRECT_ACCESS_CORE_COLUMNS
    variables = list(columns)
    if with_asymptotic:
        variables += ["asymptotic_lon_deg", "asymptotic_lat_deg"]
    lines = []
    if include_header:
        lines.extend([
            "VARIABLES=" + " ".join('"%s"' % name for name in variables),
            "ZONE T=\"synthetic\"",
            "AUXDATA CUTOFF_RECONSTRUCTION=\"self-test\"",
        ])
    for dkey in sorted(directions):
        for rigidity in rigidities:
            state = int(state_function(dkey, rigidity))
            record = synthetic_access_record(dkey, rigidity, state)
            if with_asymptotic:
                record["asymptotic_lon_deg"] = dkey[0] + 0.01
                record["asymptotic_lat_deg"] = max(
                    -90.0, min(90.0, dkey[1] + 0.01))
            values = [record[name] for name in variables]
            lines.append(" ".join(fmt_number(value) for value in values))
    path.write_text("\n".join(lines) + "\n")


def self_test(script_dir):
    """Exercise schemas, full-grid auditing, physics signs, and hard gates."""
    validate_reference_files(script_dir)
    count_contracts = ((30.0, 60, 24, 84), (15.0, 264, 48, 312),
                       (10.0, 612, 72, 684))
    for resolution, nregular, npolar, nfull in count_contracts:
        if len(expected_direction_keys(resolution, resolution)) != nregular or \
                len(expected_pole_direction_keys(resolution,
                                                 resolution)) != npolar or \
                len(expected_full_direction_keys(resolution,
                                                 resolution)) != nfull:
            raise AssertionError("%g-degree direction-grid contract changed" %
                                 resolution)
    weights = sum(weight for _, weight in fov_direction_terms(
        expected_direction_keys(15.0, 15.0), 15.0, 15.0))
    if not math.isclose(weights, 4.0 * math.pi, rel_tol=0.0, abs_tol=1.0e-12):
        raise AssertionError("Directional cell weights do not close to 4*pi")
    independent_cap = 2.0 * math.pi * (
        1.0 - math.cos(math.radians(7.5)))
    if not math.isclose(pole_cap_weight_sr(15.0), independent_cap,
                        rel_tol=0.0, abs_tol=1.0e-15):
        raise AssertionError("Polar-cap solid-angle formula changed")

    point = build_points([0.0], [0.0], 400.0)[0]
    # At lon=0/lat=0, a detector look toward +Y is local EAST.  Arrival is -Y.
    sector, _, _ = local_aperture(point, 270.0, 0.0, 100.0, 0.2)
    if sector != "EAST":
        raise AssertionError("Local EAST convention changed")
    sector, _, _ = local_aperture(point, 90.0, 0.0, 100.0, 0.2)
    if sector != "WEST":
        raise AssertionError("Local WEST convention changed")

    class Args(object):
        pass
    args = Args()
    args.dir_res = 30.0
    args.max_trace_time = 300.0
    args.dt_trace = 0.25
    args.max_steps = 1000
    args.scheduler = "DYNAMIC"
    args.dynamic_chunk = 0
    args.nt = 2
    args.fov_half_angle = 100.0
    args.east_west_deadband = 0.2
    args.transition_min = 0.05
    args.transition_max = 0.95
    args.max_lobe_unresolved = 0.20
    args.min_east_west_effect = 0.01
    with tempfile.TemporaryDirectory(prefix="C8_selftest_") as temp_name:
        temp = Path(temp_name)
        rendered = temp / "AMPS_PARAM_C8.in"
        render_input(script_dir / "AMPS_PARAM_C8_gridless.in", rendered,
                     "C8_selftest", "T05", -1, [point], [1.0, 50.0], args)
        text = rendered.read_text()
        required_patterns = (
            r"^CUTOFF_SAMPLING\s+VERTICAL\s*$",
            r"^CUTOFF_SEARCH_ALGORITHM\s+DIRECT_ACCESS\s*$",
            r"^CUTOFF_BACKTRACE_CHARGE\s+REVERSED\s*$",
            r"^FIELD_MODEL\s+T05\s*$",
            r"^CHARGE\s+-1\s*$",
        )
        if any(not re.search(pattern, text, re.MULTILINE)
               for pattern in required_patterns) or \
                re.search(r"__[A-Z0-9_]+__", text) or \
                re.search(r"^CUTOFF_BACKTRACE_CHARGE\s+SAME\s*$",
                          text, re.MULTILINE):
            raise AssertionError("Rendered C8 input is incomplete")

        regular = expected_direction_keys(90.0, 90.0)
        full = expected_full_direction_keys(90.0, 90.0)
        rigidities = [1.0, 50.0]
        state_function = lambda dkey, rigidity: int(rigidity >= 50.0)

        # A Step-4 file, the corrected 45-field Step-5 file, and the historical
        # shifted 45-field file must recover identical core samples.  This is
        # the regression for the real C8 failure: positional indexing either
        # discarded these rows or shifted state fields into unrelated columns.
        schema_paths = []
        schema_columns = (
            DIRECT_ACCESS_CORE_COLUMNS,
            DIRECT_ACCESS_CORE_COLUMNS + DIRECT_ACCESS_STEP5_COLUMNS,
            DIRECT_ACCESS_HISTORICAL_STEP5_COLUMNS,
        )
        for index, columns in enumerate(schema_columns):
            path = temp / ("schema_%d.dat" % index)
            write_synthetic_access(path, full, rigidities, state_function,
                                   columns=columns)
            schema_paths.append(path)
        parsed_schemas = [parse_directional_access(path)
                          for path in schema_paths]
        if not all(item["samples"] == parsed_schemas[0]["samples"]
                   for item in parsed_schemas[1:]):
            raise AssertionError("24/45-column named schemas changed core samples")
        expected_layouts = ("legacy_prefix_append_only",
                            "legacy_prefix_append_only",
                            "historical_step5_shifted")
        if tuple(item["schema_layout"] for item in parsed_schemas) != \
                expected_layouts:
            raise AssertionError("DIRECT_ACCESS schema provenance changed")

        # The three schemas must also produce identical physical reductions,
        # not merely the same row count.
        reduced_rows = []
        for parsed in parsed_schemas:
            sample_map = {
                (point["point_id"], dkey[0], dkey[1], rkey): sample
                for (dkey, rkey), sample in parsed["samples"].items()
            }
            rows, errors = reduce_fov(
                {"case_id": "schema", "model": "T96", "charge": 1},
                [point], sample_map, regular, rigidities, args)
            if errors:
                raise AssertionError("Valid schema reduction failed: %s" % errors)
            reduced_rows.append(rows)
        if reduced_rows[1:] != [reduced_rows[0], reduced_rows[0]]:
            raise AssertionError("24/45-column schemas reduced differently")

        def assert_parse_failure(path, expected_fragment):
            try:
                parse_directional_access(path)
            except RuntimeError as exc:
                if expected_fragment.lower() not in str(exc).lower():
                    raise AssertionError(
                        "%s failed for the wrong reason: %s" % (path, exc))
                return
            raise AssertionError("Malformed schema was accepted: %s" % path)

        access_path = schema_paths[0]
        # Duplicate data rows and state/reason disagreement are distinct defects.
        duplicate_path = temp / "duplicate.dat"
        duplicate_path.write_text(access_path.read_text() +
                                  access_path.read_text().splitlines()[-1] + "\n")
        assert_parse_failure(duplicate_path, "duplicate direction")

        broken_lines = access_path.read_text().splitlines()
        first_values = broken_lines[3].split()
        first_values[DIRECT_ACCESS_CORE_COLUMNS.index("termination_code")] = "0"
        broken_lines[3] = " ".join(first_values)
        broken_path = temp / "broken.dat"
        broken_path.write_text("\n".join(broken_lines) + "\n")
        broken = parse_directional_access(broken_path)
        if not any(item["contract_error"] for item in broken["samples"].values()):
            raise AssertionError("State/reason mismatch was not detected")

        # Schema failures must be explicit: missing and duplicate names, row
        # widths both below and above the header, and blank numeric fields.
        missing_path = temp / "missing_core.dat"
        missing_columns = tuple(name for name in DIRECT_ACCESS_CORE_COLUMNS
                                if name != "energy_mev")
        write_synthetic_access(missing_path, {direction_key(0.0, 0.0)}, [1.0],
                               state_function, columns=missing_columns)
        assert_parse_failure(missing_path, "lacks required")

        duplicate_header_path = temp / "duplicate_header.dat"
        write_synthetic_access(
            duplicate_header_path, {direction_key(0.0, 0.0)}, [1.0],
            state_function, columns=DIRECT_ACCESS_CORE_COLUMNS + ("lon_deg",))
        assert_parse_failure(duplicate_header_path, "duplicate variables")

        current_lines = schema_paths[1].read_text().splitlines()
        short_path = temp / "short_row.dat"
        short_lines = list(current_lines)
        short_lines[3] = " ".join(short_lines[3].split()[:-1])
        short_path.write_text("\n".join(short_lines) + "\n")
        assert_parse_failure(short_path, "row width")
        long_path = temp / "long_row.dat"
        long_lines = list(current_lines)
        long_lines[3] += " 0"
        long_path.write_text("\n".join(long_lines) + "\n")
        assert_parse_failure(long_path, "row width")
        blank_path = temp / "blank_value.dat"
        blank_lines = access_path.read_text().splitlines()
        blank_values = blank_lines[3].split()
        blank_values[DIRECT_ACCESS_CORE_COLUMNS.index("energy_mev")] = '""'
        blank_lines[3] = " ".join(blank_values)
        blank_path.write_text("\n".join(blank_lines) + "\n")
        assert_parse_failure(blank_path, "blank variables value")

        # The only positional fallback is an exact 24-value headerless row.
        headerless_path = temp / "headerless_legacy.dat"
        write_synthetic_access(
            headerless_path, {direction_key(0.0, 0.0)}, [1.0],
            state_function, columns=DIRECT_ACCESS_CORE_COLUMNS,
            include_header=False)
        if len(parse_directional_access(headerless_path)["samples"]) != 1:
            raise AssertionError("Exact headerless legacy row did not parse")

        # C8-G01/G11: a complete regular+polar cube passes; missing polar rows
        # and one inconsistent polar duplicate each fail their own hard gate.
        def audit_synthetic(dirname, directions, state_fn):
            workdir = temp / dirname
            workdir.mkdir()
            write_synthetic_access(
                workdir / "cutoff_gridless_dir_access_point_0000.dat",
                directions, rigidities, state_fn,
                columns=DIRECT_ACCESS_CORE_COLUMNS +
                DIRECT_ACCESS_STEP5_COLUMNS)
            case = {"case_id": dirname, "model": "T96", "charge": 1,
                    "workdir": workdir}
            return audit_case(case, [point], regular, full, rigidities)

        summaries, sample_map, totals, errors = audit_synthetic(
            "complete", full, state_function)
        if errors or totals["coverage_errors"] or \
                totals["pole_degeneracy_errors"] or not summaries[0]["passed"]:
            raise AssertionError("Complete polar cube failed audit: %s" % errors)
        if not case_is_structurally_auditable(
                totals, errors, sample_map, len(full) * len(rigidities)):
            raise AssertionError("Complete cube was blocked from reductions")
        _, _, missing_totals, missing_errors = audit_synthetic(
            "missing_poles", regular, state_function)
        if not missing_totals["coverage_errors"] or not missing_errors:
            raise AssertionError("Missing polar rows passed C8-G01")
        if case_is_structurally_auditable(
                missing_totals, missing_errors, {},
                len(full) * len(rigidities)):
            raise AssertionError("Partial cube could reach a blank-fraction reduction")

        def corrupt_one_pole(dkey, rigidity):
            base = state_function(dkey, rigidity)
            if dkey == direction_key(90.0, -90.0) and rigidity == 1.0:
                return 1 - base
            return base
        corrupt_summaries, _, corrupt_totals, corrupt_errors = audit_synthetic(
            "corrupt_pole", full, corrupt_one_pole)
        if corrupt_totals["pole_degeneracy_errors"] != 1 or \
                corrupt_summaries[0]["pole_degeneracy_errors"] != 1 or \
                not any("pole azimuthal degeneracy" in item
                        for item in corrupt_errors):
            raise AssertionError("Corrupt polar duplicate did not fail C8-G11")

        # C8-G07 must accept the physical charge-odd sign and reject its exact
        # mirror.  No threshold is changed by the REVERSED-charge fix.
        def synthetic_ew_rows(mirrored=False):
            rows = []
            for model in ("T96", "T05"):
                for charge in (+1, -1):
                    sign = charge * (-1 if mirrored else 1)
                    rows.append({
                        "field_model": model, "charge": charge,
                        "rigidity_GV": 10.0,
                        "east_allowed_fraction": 0.5 - 0.1 * sign,
                        "west_allowed_fraction": 0.5 + 0.1 * sign,
                        "fov_allowed_fraction": 0.5,
                        "east_unresolved_fraction": 0.0,
                        "west_unresolved_fraction": 0.0,
                    })
            return rows
        if not east_west_comparison(synthetic_ew_rows(False), args)[1]:
            raise AssertionError("Physical East-West sign failed C8-G07")
        if east_west_comparison(synthetic_ew_rows(True), args)[1]:
            raise AssertionError("Mirrored East-West sign passed C8-G07")
    print("C8 self-test: PASS")


def apply_profile_defaults(args):
    """Fill cost controls from the selected profile unless explicitly set."""
    profile = PROFILE_DEFAULTS[args.profile]
    if args.dir_res is None:
        args.dir_res = profile["dir_res_deg"]
    if args.max_trace_time is None:
        args.max_trace_time = profile["max_trace_time_s"]
    if args.max_steps is None:
        args.max_steps = profile["max_steps"]


def build_parser():
    """Construct the documented command-line interface in one auditable place."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=tuple(PROFILE_DEFAULTS),
                        default="ROUTINE")
    parser.add_argument("--amps", default="./amps",
                        help="AMPS executable, relative to the launch directory")
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("-np", type=int, default=4, help="MPI ranks per AMPS case")
    parser.add_argument("-nt", type=int, default=16,
                        help="GRIDLESS worker threads per MPI rank")
    parser.add_argument("--mover", choices=("BORIS", "RK2", "RK4", "RK6"),
                        default="BORIS")
    parser.add_argument("--scheduler", choices=("DYNAMIC", "STATIC", "BLOCK_CYCLIC"),
                        default="DYNAMIC")
    parser.add_argument("--dynamic-chunk", type=int, default=0)
    parser.add_argument("--dir-res", type=float, default=None,
                        help="common longitude/latitude resolution in degrees")
    parser.add_argument("--rigidities", default=",".join(
        fmt_number(value) for value in DEFAULT_RIGIDITIES_GV))
    parser.add_argument("--lons", default=",".join(
        fmt_number(value) for value in DEFAULT_LONS_DEG))
    parser.add_argument("--lats", default=",".join(
        fmt_number(value) for value in DEFAULT_LATS_DEG))
    parser.add_argument("--alt-km", type=float, default=400.0)
    parser.add_argument("--dt-trace", type=float, default=0.25)
    parser.add_argument("--max-trace-time", type=float, default=None)
    parser.add_argument("--max-steps", type=int, default=None)
    parser.add_argument("--fov-half-angle", type=float, default=45.0)
    parser.add_argument("--east-west-deadband", type=float, default=0.25)
    parser.add_argument("--max-unresolved-fraction", type=float, default=0.25)
    parser.add_argument("--max-lobe-unresolved", type=float, default=0.20)
    parser.add_argument("--high-rigidity-min", type=float, default=50.0)
    parser.add_argument("--min-high-rigidity-transparency", type=float,
                        default=0.98)
    parser.add_argument("--max-high-rigidity-unresolved", type=float,
                        default=0.02)
    parser.add_argument("--max-model-mismatch", type=float, default=0.15)
    parser.add_argument("--transition-min", type=float, default=0.05)
    parser.add_argument("--transition-max", type=float, default=0.95)
    parser.add_argument("--min-east-west-effect", type=float, default=0.01)
    parser.add_argument("--asymptotic-policy", choices=("DIAGNOSTIC", "REQUIRE"),
                        default="DIAGNOSTIC")
    parser.add_argument("--max-asymptotic-rms-mrad", type=float, default=5.0)
    parser.add_argument("--workdir", default="test_output/C8_directional_access")
    parser.add_argument("--keep", action="store_true",
                        help="retain an existing output directory before running")
    parser.add_argument("--skip-run", action="store_true",
                        help="audit existing case directories without launching AMPS")
    parser.add_argument("--dry-run", action="store_true",
                        help="render inputs and print commands only")
    parser.add_argument("--validate-references", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser


def validate_arguments(args, rigidities, points):
    """Reject settings that would weaken or geometrically invalidate C8."""
    if args.np < 1 or args.nt < 1 or args.dynamic_chunk < 0:
        raise SystemExit("-np/-nt must be positive and dynamic chunk non-negative")
    expected_direction_keys(args.dir_res, args.dir_res)
    if any(value <= 0.0 for value in rigidities) or \
            sorted(set(rigidities)) != sorted(rigidities):
        raise SystemExit("Rigidities must be unique and positive")
    if not any(value >= args.high_rigidity_min for value in rigidities):
        raise SystemExit("Rigidity list needs a high-rigidity transparency node")
    if not any(args.transition_min < 0.5 < args.transition_max for _ in [0]):
        raise SystemExit("Transition bounds must contain 0.5")
    if not 0.0 < args.fov_half_angle < 90.0:
        raise SystemExit("FOV half-angle must be between 0 and 90 degrees")
    if not 0.0 <= args.east_west_deadband < 1.0:
        raise SystemExit("East-West deadband must be in [0,1)")
    if args.alt_km <= 0.0 or not points:
        raise SystemExit("C8 requires positive altitude and observation points")


def main(argv=None):
    """Render, run, audit, reduce, and report the complete C8 validation."""
    parser = build_parser()
    args = parser.parse_args(argv)
    apply_profile_defaults(args)
    script_dir = Path(__file__).resolve().parent

    if args.self_test:
        self_test(script_dir)
        return 0
    try:
        published_reference, acceptance_reference, expected_reference = \
            validate_reference_files(script_dir)
    except Exception as exc:
        raise SystemExit("C8 reference validation failed: %s" % exc)
    if args.validate_references:
        print("C8 references: PASS (%d published rows, %d acceptance gates, "
              "%d expected-physics rows)" %
              (len(published_reference), len(acceptance_reference),
               len(expected_reference)))
        return 0

    rigidities = parse_float_list(args.rigidities, "--rigidities")
    points = build_points(parse_float_list(args.lons, "--lons"),
                          parse_float_list(args.lats, "--lats"), args.alt_km)
    validate_arguments(args, rigidities, points)
    expected_directions = expected_direction_keys(args.dir_res, args.dir_res)
    expected_pole_directions = expected_pole_direction_keys(
        args.dir_res, args.dir_res)
    expected_full_directions = expected_full_direction_keys(
        args.dir_res, args.dir_res)

    launch_dir = Path.cwd().resolve()
    base_workdir = (launch_dir / args.workdir).resolve()
    if not args.skip_run:
        if base_workdir.exists() and not args.keep:
            shutil.rmtree(base_workdir)
        base_workdir.mkdir(parents=True, exist_ok=True)
    elif not base_workdir.exists():
        raise SystemExit("--skip-run workdir does not exist: %s" % base_workdir)

    amps_path = Path(args.amps)
    if not amps_path.is_absolute():
        amps_path = (launch_dir / amps_path).resolve()
    template = script_dir / "AMPS_PARAM_C8_gridless.in"
    driver = script_dir / "data" / "ts05_driver_C8.txt"
    cases = {}
    run_records = []
    overall_passed = True
    messages = []

    for model in ("T96", "T05"):
        for charge in (+1, -1):
            charge_label = "charge_plus" if charge > 0 else "charge_minus"
            case_id = "%s_%s" % (model.lower(), charge_label)
            workdir = base_workdir / model.lower() / charge_label
            case = {"case_id": case_id, "model": model, "charge": charge,
                    "workdir": workdir}
            cases[(model, charge)] = case
            if not args.skip_run:
                workdir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(str(driver), str(workdir / "ts05_driver_C8.txt"))
                render_input(template, workdir / "AMPS_PARAM_C8.in", case_id,
                             model, charge, points, rigidities, args)
            command = [
                args.mpirun, "-np", str(args.np), str(amps_path),
                "-mode", "gridless", "-i", "AMPS_PARAM_C8.in",
                "-mover", args.mover,
                "-gridless-parallel", "THREADS",
                "-gridless-threads", str(args.nt),
                "-gridless-mpi-scheduler", args.scheduler,
                "-gridless-mpi-dynamic-chunk", str(args.dynamic_chunk),
                "-cutoff-search", "DIRECT_ACCESS",
                "-cutoff-rigidity-list-gv", ",".join(
                    fmt_number(value) for value in rigidities),
                "-cutoff-trace-policy", "ACCURATE",
            ]
            log_path = workdir / "C8_AMPS.log"
            case["command"] = command
            run_records.append({
                "case_id": case_id, "workdir": str(workdir),
                "command": command, "log_file": str(log_path),
            })
            if args.dry_run:
                print("[%s] %s" % (case_id, " ".join(command)))
            elif not args.skip_run:
                print("\nC8 %s command:\n  %s" %
                      (case_id, " ".join(command)))
                rc = run_command(command, workdir, log_path)
                case["return_code"] = rc
                if rc != 0:
                    overall_passed = False
                    messages.append("%s: AMPS exited with code %d" %
                                    (case_id, rc))

    if args.dry_run:
        return 0

    coverage_rows = []
    fov_rows = []
    for key, case in cases.items():
        if case.get("return_code", 0) != 0:
            case["sample_map"] = {}
            case["auditable"] = False
            case["reduction_ok"] = False
            continue
        summaries, sample_map, totals, audit_messages = audit_case(
            case, points, expected_directions, expected_full_directions,
            rigidities)
        case["audit_totals"] = dict(totals)
        case["schema_layouts"] = sorted({row["schema_layout"]
                                         for row in summaries})
        coverage_rows.extend(summaries)
        messages.extend(audit_messages)
        expected_case_samples = (len(points) * len(expected_full_directions) *
                                 len(rigidities))
        structurally_auditable = case_is_structurally_auditable(
            totals, audit_messages, sample_map, expected_case_samples)
        case["auditable"] = structurally_auditable
        # Never pass a partial cube into a reduction.  In the original failure,
        # a parse/schema problem left empty lobe totals; float("") then aborted
        # before C8 could write any comparison artifact.  Structural gates now
        # fail first and the reporting phase still writes explicit NO_ROWS CSVs.
        case["sample_map"] = sample_map if structurally_auditable else {}
        if not structurally_auditable:
            overall_passed = False
            if totals["contract_errors"]:
                messages.append("%s has %d state/termination contract errors" %
                                (case["case_id"], totals["contract_errors"]))
            if totals["fatal"]:
                messages.append("%s has %d fatal numerical terminations" %
                                (case["case_id"], totals["fatal"]))
            if totals["pole_degeneracy_errors"] and not any(
                    "pole azimuthal degeneracy" in item
                    for item in audit_messages):
                messages.append("%s has %d polar degeneracy errors" %
                                (case["case_id"],
                                 totals["pole_degeneracy_errors"]))
            if totals["samples"] != expected_case_samples:
                messages.append("%s has %d/%d auditable samples" %
                                (case["case_id"], totals["samples"],
                                 expected_case_samples))
        unresolved_fraction = (float(totals["state_2"]) / totals["samples"]
                               if totals["samples"] else 1.0)
        case["unresolved_fraction"] = unresolved_fraction
        if unresolved_fraction > args.max_unresolved_fraction:
            overall_passed = False
            messages.append("%s unresolved fraction %.3f exceeds %.3f" %
                            (case["case_id"], unresolved_fraction,
                             args.max_unresolved_fraction))
        if not structurally_auditable:
            case["reduction_ok"] = False
            continue
        reduced, reduction_errors = reduce_fov(
            case, points, sample_map, expected_directions, rigidities, args)
        fov_rows.extend(reduced)
        messages.extend(reduction_errors)
        case["reduction_ok"] = not reduction_errors
        if reduction_errors:
            overall_passed = False

    # High-rigidity transparency is checked per point/case so an average cannot
    # conceal one broken local-time or latitude configuration.
    for row in fov_rows:
        if row["rigidity_GV"] < args.high_rigidity_min:
            continue
        allowed = float(row["fov_allowed_fraction"])
        unresolved = float(row["fov_unresolved_fraction"])
        if allowed < args.min_high_rigidity_transparency or \
                unresolved > args.max_high_rigidity_unresolved:
            overall_passed = False
            messages.append(
                "%s point %d R=%g high-rigidity transparency failed: "
                "allowed=%.3f unresolved=%.3f" %
                (row["case_id"], row["point_id"], row["rigidity_GV"],
                 allowed, unresolved))

    if all(cases[key].get("auditable", False) and
           cases[key].get("reduction_ok", False) for key in cases):
        model_rows, model_passed, model_messages = model_comparison(
            cases, points, expected_directions, rigidities, args)
        ew_rows, ew_passed, ew_messages = east_west_comparison(fov_rows, args)
        asym_rows, asym_summary, asym_available, asym_passed, asym_messages = \
            asymptotic_comparison(cases, points, expected_directions,
                                  rigidities, args)
        overall_passed = (overall_passed and model_passed and ew_passed and
                          asym_passed)
        messages.extend(model_messages + ew_messages + asym_messages)
    else:
        model_rows, ew_rows, asym_rows, asym_summary = [], [], [], []
        asym_available = False
        overall_passed = False
        messages.append("One or more C8 cases lack an auditable sample cube")

    base_workdir.mkdir(parents=True, exist_ok=True)
    write_csv(base_workdir / "C8_coverage_summary.csv", coverage_rows)
    write_csv(base_workdir / "C8_directional_access_summary.csv", fov_rows)
    write_csv(base_workdir / "C8_model_comparison.csv", model_rows)
    write_csv(base_workdir / "C8_east_west_comparison.csv", ew_rows)
    write_csv(base_workdir / "C8_asymptotic_comparison.csv", asym_rows,
              fieldnames=("charge", "point_id", "sky_lon_deg", "sky_lat_deg",
                          "rigidity_GV", "rigidity_bin",
                          "latitude_difference_mrad",
                          "longitude_difference_mrad", "direction_weight_sr",
                          "t96_schema", "t05_schema"))
    shutil.copy2(str(script_dir / "reference_C8_boschini2013.csv"),
                 str(base_workdir / "reference_C8_boschini2013.csv"))
    shutil.copy2(str(script_dir / "reference_C8_acceptance_contract.csv"),
                 str(base_workdir / "reference_C8_acceptance_contract.csv"))
    shutil.copy2(str(script_dir / "reference_C8_expected_physics.csv"),
                 str(base_workdir / "reference_C8_expected_physics.csv"))

    result = {
        "test_id": TEST_ID, "test_name": TEST_NAME,
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "runner_release": RUNNER_RELEASE, "passed": overall_passed,
        "profile": args.profile, "mover": args.mover, "np": args.np,
        "nt": args.nt, "epoch": DEFAULT_EPOCH, "points": points,
        "rigidities_GV": rigidities, "direction_resolution_deg": args.dir_res,
        "requested_direction_tasks_per_point": len(expected_directions),
        "requested_pole_direction_tasks_per_point":
            len(expected_pole_directions),
        "requested_full_direction_tasks_per_point":
            len(expected_full_directions),
        "fov_half_angle_deg": args.fov_half_angle,
        "asymptotic_policy": args.asymptotic_policy,
        "asymptotic_columns_available": asym_available,
        "asymptotic_summary": asym_summary,
        "cases": [{
            "case_id": case["case_id"], "model": case["model"],
            "charge": case["charge"], "return_code": case.get("return_code", 0),
            "auditable": case.get("auditable", False),
            "schema_layouts": case.get("schema_layouts", []),
            "unresolved_fraction": case.get("unresolved_fraction"),
        } for case in cases.values()],
        "messages": messages, "runs": run_records,
    }
    (base_workdir / "C8_result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n")
    summary_lines = [
        "C8 — %s" % TEST_NAME,
        "result: %s" % ("PASS" if overall_passed else "FAIL"),
        "profile: %s; cases: 4; points: %d; directions/point: %d "
        "(%d regular + %d polar); "
        "rigidities: %d" %
        (args.profile, len(points), len(expected_full_directions),
         len(expected_directions), len(expected_pole_directions),
         len(rigidities)),
        "MPI ranks / threads: %d / %d" % (args.np, args.nt),
        "asymptotic direction columns: %s (policy=%s)" %
        ("AVAILABLE" if asym_available else "NOT_AVAILABLE",
         args.asymptotic_policy),
        "",
    ]
    for case in cases.values():
        summary_lines.append("%-24s unresolved=%s" %
                             (case["case_id"],
                              ("n/a" if case.get("unresolved_fraction") is None
                               else "%.3f" % case["unresolved_fraction"])))
    if messages:
        summary_lines += ["", "Failure messages:"] + ["- " + item for item in messages]
    summary_lines += ["", "C8 results: %s" % base_workdir]
    summary_text = "\n".join(summary_lines) + "\n"
    (base_workdir / "C8_summary.txt").write_text(summary_text)
    print(summary_text)
    return 0 if overall_passed else 1


if __name__ == "__main__":
    sys.exit(main())

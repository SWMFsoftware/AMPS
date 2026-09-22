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

-------------------------------------------------------------------------------
ERRATA / FIX HISTORY (read this before changing acceptance thresholds)
-------------------------------------------------------------------------------
Two independent defects were found by auditing a ROUTINE C8 run that failed
with "T96/T05 has no informative rigidity with D(+)>=0.01 and D(-)<=-0.01".
Both are fixed in this revision.  Neither fix relaxes an acceptance threshold;
both fixes correct the test's own bookkeeping/configuration so that the
existing, unchanged thresholds are evaluated against physically and
structurally correct inputs.

FIX 1 — backtrace-charge convention (AMPS_PARAM_C8_gridless.in).
  ``CUTOFF_BACKTRACE_CHARGE`` was set to ``SAME``.  AMPS backtraces a particle
  by launching a forward-in-time numerical trajectory from the observation
  point with velocity reversed relative to the requested arrival direction.
  Time-reversing q(v x B) under that velocity flip additionally requires
  flipping the charge sign used in the integration (the standard cosmic-ray
  "antiparticle" cutoff construction).  ``REVERSED`` performs that flip;
  ``SAME`` does not and is documented in the AMPS source as a legacy mode.
  Under ``SAME`` the case labeled ``charge_plus`` numerically reproduces the
  physics of a real NEGATIVE particle and vice versa, which flips the sign of
  the computed East-West asymmetry for both charges -- exactly the failure
  observed.  Fixed by selecting ``REVERSED``, matching every other
  cutoff-physics test in this suite (C6, C7, C9, C10, C13, C15, C16, C18,
  C19).  See the comment block in AMPS_PARAM_C8_gridless.in for the full
  derivation.

FIX 2 — directional-grid completeness contract (this file).
  ``expected_direction_keys()`` assumed AMPS's DIRECT_ACCESS producer omits
  the longitude-degenerate +/-90-degree rows.  It does not: AMPS's directional
  grid unconditionally includes both poles, one row per longitude value (see
  ``nLatMap = ... + 1  // include poles`` in CutoffRigidityGridless.cpp), so a
  ROUTINE (15-degree) run emits 24 extra rows at each pole (48 total) beyond
  the 264-direction non-polar grid this file expected.  Because Section 5.1 /
  gate C8-G01 is a HARD "zero coverage errors" gate, those 48 legitimate rows
  were flagged as "unexpected directions" and silently failed the run via the
  ``coverage_errors`` counter -- silently, because that counter fed
  ``overall_passed`` without ever appending a message, so C8_summary.txt gave
  no indication this gate had failed at all.
  Fixed with three changes, all of which TIGHTEN the test rather than loosen
  it:
    (a) ``expected_full_direction_keys()`` now defines the true output
        contract (264 regular cells + 48 polar cells for ROUTINE) used by the
        Section-5.1 completeness/uniqueness audit, so genuine AMPS rows are
        no longer misclassified as errors.
    (b) A NEW hard check (gate C8-G11, "pole azimuthal degeneracy") verifies
        that the 24 longitude-tagged duplicates AMPS writes at each pole
        report an IDENTICAL access state per rigidity, since they all
        represent the same physical direction.  A disagreement among these
        duplicates is a genuine trace non-determinism/numerical defect that
        the previous version of C8 could never have detected, because it
        discarded pole rows outright instead of auditing them.
    (c) The FOV/East-West/model-comparison solid-angle reduction previously
        *approximated* the two polar caps by silently stretching the
        adjacent 15-degree latitude band all the way to the pole.  Now that
        genuine (and verified self-consistent) polar samples exist, each
        pole is folded into the reduction as its own direction with the
        exact polar-cap solid angle (``pole_cap_weight_sr``), replacing an
        approximation with the real traced value.  Total solid angle is
        unchanged (still exactly 4*pi sr, reapportioned rather than altered),
        so this does not change what "resolved"/"unresolved" or the 0.98 /
        0.15 / 0.01 thresholds mean -- it only makes the number that is
        compared against those thresholds slightly more accurate for the two
        observation points (|lat|=51.6 deg) whose 45-degree zenith cone
        reaches a pole.

Neither fix changes any CLI default or acceptance threshold in build_parser()
or reference_C8_expected_physics.csv / reference_C8_acceptance_contract.csv
(gate limits 0.98, 0.02, 0.25, 0.15, 0.01, 5 mrad are untouched); gate
C8-G11 was added to reference_C8_acceptance_contract.csv as a strictly
additional hard gate.
-------------------------------------------------------------------------------
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
# Schema version 2: adds the full (regular + polar) direction contract, gate
# C8-G11 (pole azimuthal degeneracy), and the
# requested_pole_direction_tasks_per_point / requested_full_direction_tasks_per_point
# result fields.  See the "ERRATA / FIX HISTORY" note at the top of this file.
RUNNER_SCHEMA_VERSION = 2
RUNNER_RELEASE = "2026-09-02"

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


def _direction_grid_counts(lon_res_deg, lat_res_deg):
    """Validate and return (nlon, nlat_intervals) shared by every grid helper.

    Factored out so the regular-grid, pole-grid, and full-grid builders below
    can never disagree about how many longitude/latitude steps a resolution
    implies.  ``nlat_intervals`` is the number of *non-polar* latitude bands
    (11 for 15 degrees); the two poles are handled separately because AMPS
    handles them separately (see ``expected_pole_direction_keys``).
    """
    nlon_float = 360.0 / float(lon_res_deg)
    nlat_float = 180.0 / float(lat_res_deg)
    nlon = int(round(nlon_float))
    nlat_intervals = int(round(nlat_float))
    if abs(nlon - nlon_float) > 1.0e-9 or abs(nlat_intervals - nlat_float) > 1.0e-9:
        raise ValueError("Directional resolution must divide 360 and 180 degrees")
    return nlon, nlat_intervals


def expected_direction_keys(lon_res_deg, lat_res_deg):
    """Derive the complete REGULAR (non-polar) sky grid requested from AMPS.

    This is the grid used for solid-angle-weighted reductions (FOV/EAST/WEST
    integration, T96/T05 model comparison, asymptotic-direction comparison):
    see ``sky_cell_weight_sr`` and ``fov_direction_terms``.  It intentionally
    excludes the two poles.

    NOTE (fix history): this function's docstring used to claim that AMPS's
    DIRECT_ACCESS writer omits the longitude-degenerate +/-90-degree rows
    entirely.  That is false -- AMPS's directional grid always includes both
    poles, with one row per longitude value (see
    ``nLatMap = ... + 1  // include poles`` in CutoffRigidityGridless.cpp).
    Excluding poles here is still the *right* choice for the solid-angle
    reduction (a regular lon/lat cell is a poor way to represent 24 identical
    duplicate directions), but the OUTPUT-COMPLETENESS audit must not use
    this function alone any more -- use ``expected_full_direction_keys``
    for that.  See ``expected_pole_direction_keys`` and
    ``pole_cap_weight_sr`` for how the poles are folded back in correctly.
    """
    nlon, nlat_intervals = _direction_grid_counts(lon_res_deg, lat_res_deg)
    return {
        direction_key(i * lon_res_deg, -90.0 + j * lat_res_deg)
        for j in range(1, nlat_intervals)
        for i in range(nlon)
    }


def expected_pole_direction_keys(lon_res_deg, lat_res_deg):
    """Derive the pole rows AMPS actually writes for a given resolution.

    AMPS's directional-map cell count is ``nLonMap * nLatMap`` with
    ``nLatMap = round(180/latRes) + 1`` (poles included), and it writes one
    row per (longitude, pole-latitude) pair even though every one of those
    rows is the same physical direction (cos(+/-90 deg) == 0, so the
    longitude has no effect on the traced Cartesian direction).  This
    function reproduces that same per-longitude duplication so the
    completeness audit expects exactly what AMPS produces, and so the new
    pole-degeneracy check (see ``audit_case``) has the full set of duplicate
    keys to compare against each other.
    """
    nlon, _ = _direction_grid_counts(lon_res_deg, lat_res_deg)
    return {
        direction_key(i * lon_res_deg, pole_lat)
        for pole_lat in (-90.0, 90.0)
        for i in range(nlon)
    }


def expected_full_direction_keys(lon_res_deg, lat_res_deg):
    """Return the true DIRECT_ACCESS output contract: regular cells + poles.

    This is the set used by the Section-5.1 / gate C8-G01 "no missing or
    unexpected directions" audit.  Using the regular-only set there (the
    pre-fix behavior) misclassified every genuine pole row AMPS writes as an
    "unexpected direction," silently failing gate C8-G01 on every ROUTINE
    run without ever explaining why.
    """
    return expected_direction_keys(lon_res_deg, lat_res_deg) | \
        expected_pole_direction_keys(lon_res_deg, lat_res_deg)


def sky_cell_weight_sr(lon_res_deg, lat_res_deg, sky_lat_deg):
    """Exact solid angle of one regular (non-polar) longitude/latitude cell.

    NOTE (fix history): earlier versions of this function stretched the
    first and last non-polar latitude bands all the way to the corresponding
    pole, because at that time no genuine polar sample existed to carry that
    solid angle.  AMPS *does* write genuine (duplicate-but-verifiable) polar
    samples -- see ``expected_pole_direction_keys`` and the new pole
    degeneracy check in ``audit_case`` -- so every regular band, including
    the first/last ones, is now a plain, unstretched cell, and the two polar
    caps are folded into reductions separately with their own exact weight
    via ``pole_cap_weight_sr``.  ``sky_cell_weight_sr`` plus two calls to
    ``pole_cap_weight_sr`` still sum to exactly 4*pi sr (verified in
    ``self_test``); the total solid angle represented is unchanged, only its
    apportionment between "regular band" and "explicit pole sample" is more
    faithful to what AMPS actually traced.
    """
    lower = max(-90.0, float(sky_lat_deg) - 0.5 * lat_res_deg)
    upper = min(90.0, float(sky_lat_deg) + 0.5 * lat_res_deg)
    return math.radians(lon_res_deg) * (
        math.sin(math.radians(upper)) - math.sin(math.radians(lower)))


def pole_cap_weight_sr(lat_res_deg):
    """Exact solid angle of the polar cap represented by ONE pole direction.

    The cap has angular radius ``lat_res_deg/2`` around the pole (it fills
    exactly the gap left by the now-unstretched regular bands, see
    ``sky_cell_weight_sr``).  Its solid angle is the standard spherical-cap
    formula ``2*pi*(1 - cos(half_angle))``, written here in the same
    ``sin(upper) - sin(lower)`` form used throughout this file for exact
    consistency with ``sky_cell_weight_sr``:

        upper = 90 deg (the pole itself)
        lower = 90 deg - lat_res_deg/2 (edge of the last regular band)

    This weight is independent of longitude resolution and is assigned ONCE
    per pole (to a single canonical direction), not once per AMPS longitude
    duplicate -- the 24-or-so duplicate rows AMPS writes at a pole are all
    the same physical direction, not 24 physically distinct sub-cells, so
    summing this weight once per duplicate would double- (or 24x-) count the
    cap.  See ``audit_case`` for how a single canonical pole sample is
    selected out of the verified-identical duplicates.
    """
    half_angle_deg = 0.5 * float(lat_res_deg)
    return 2.0 * math.pi * (1.0 - math.sin(math.radians(90.0 - half_angle_deg)))


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
    """Parse a DIRECT_ACCESS cube with strict schema and duplicate checks."""
    variables = []
    reading_variables = False
    samples = {}

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
                names = re.findall(r'"([^"]+)"', line)
                if names:
                    variables.extend(normalize_variable_name(name) for name in names)
                    continue
                reading_variables = False
            if upper.startswith(("TITLE", "#", "!")):
                continue
            parts = line.replace(",", " ").split()
            if not variables:
                continue
            if len(parts) != len(variables):
                raise RuntimeError("%s:%d row width %d does not match %d VARIABLES" %
                                   (path, line_number, len(parts), len(variables)))
            record = dict(zip(variables, parts))
            required = ("lon_deg", "lat_deg", "rigidity_gv", "access_state",
                        "termination_code", "trace_time_s",
                        "trace_distance_re", "trace_steps")
            missing = [name for name in required if name not in record]
            if missing:
                raise RuntimeError("%s lacks required C8 columns: %s" %
                                   (path, ", ".join(missing)))
            try:
                lon = float(record["lon_deg"])
                lat = float(record["lat_deg"])
                rigidity = float(record["rigidity_gv"])
                state = int(float(record["access_state"]))
                code = int(float(record["termination_code"]))
                trace_time = float(record["trace_time_s"])
                trace_distance = float(record["trace_distance_re"])
                trace_steps = int(float(record["trace_steps"]))
            except (ValueError, OverflowError) as exc:
                raise RuntimeError("Malformed row %d in %s: %s" %
                                   (line_number, path, exc))
            numeric = (lon, lat, rigidity, trace_time, trace_distance,
                       float(trace_steps))
            if not all(math.isfinite(value) for value in numeric):
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
                "retry_count": optional_float(record, "retry_count"),
                "contract_error": contract_error,
                "fatal": fatal,
                "asymptotic": asymptotic_lon_lat(record),
            }
    if not samples:
        raise RuntimeError("No DIRECT_ACCESS samples parsed from %s" % path)
    return {"columns": tuple(variables), "samples": samples}


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


def audit_case(case, points, expected_directions, pole_directions, rigidities_gv):
    """Audit complete output coverage, verify pole self-consistency, and
    build a common sample map.

    ``expected_directions`` is the regular (non-polar) grid.  ``pole_directions``
    is the set of per-longitude pole duplicates AMPS actually writes (see
    ``expected_pole_direction_keys``).  Their union is the true DIRECT_ACCESS
    output contract audited here for gate C8-G01 ("no missing or unexpected
    directions").  Using the regular-only grid for that audit (the pre-fix
    behavior) misclassified every genuine pole row as "unexpected" and
    silently failed every ROUTINE run -- see the module docstring's "FIX 2"
    note.

    A second, independent hard check implements the new gate C8-G11: every
    longitude-tagged duplicate AMPS writes at a given pole must report the
    SAME three-state access classification for a given rigidity, since
    cos(+/-90 deg) makes the traced Cartesian direction identical regardless
    of the longitude label.  A disagreement is a genuine trace
    non-determinism/numerical defect that the pre-fix runner could never
    observe because it discarded pole rows before auditing them.  When all
    duplicates for a pole/rigidity agree, ONE canonical sample (longitude
    0.0) is published into ``sample_map`` so the FOV/model/asymptotic
    reductions (see ``fov_direction_terms``) have a single, verified
    representative direction per pole instead of 24 redundant ones.

    Any coverage discrepancy (missing direction, unexpected direction,
    missing sample, or unrecognized extra rigidity) now also appends a
    specific message to the returned ``errors`` list.  Previously these
    discrepancies only incremented a counter that fed ``overall_passed``
    without any accompanying message, so a hard gate could fail completely
    silently; that reporting gap is fixed here as well.
    """
    full_directions = expected_directions | pole_directions
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
        missing_directions = full_directions - got_directions
        extra_directions = got_directions - full_directions
        matched_keys = set()
        point_counts = Counter()
        extra_rigidities = 0

        # Per (pole latitude, matched rigidity): {access_state -> {longitudes}}.
        # Used only to check gate C8-G11 below; regular (non-polar) directions
        # are not azimuthally degenerate and are not tracked here.
        pole_states = defaultdict(lambda: defaultdict(set))

        for (dkey, _), sample in parsed["samples"].items():
            requested = match_requested_rigidity(sample["rigidity_GV"], rigidities_gv)
            if dkey not in full_directions or requested is None:
                if requested is None:
                    extra_rigidities += 1
                continue
            rkey = rigidity_key(requested)
            out_key = (point["point_id"], dkey[0], dkey[1], rkey)
            if out_key in sample_map:
                errors.append("%s duplicate matched sample %s" %
                              (case["case_id"], out_key))
                continue
            sample_map[out_key] = sample
            matched_keys.add((dkey, rkey))
            point_counts["samples"] += 1
            point_counts["state_%d" % sample["access_state"]] += 1
            point_counts["contract_errors"] += int(sample["contract_error"])
            point_counts["fatal"] += int(sample["fatal"])

            if dkey in pole_directions:
                pole_states[(dkey[1], rkey)][sample["access_state"]].add(dkey[0])

        # --- Gate C8-G11: pole azimuthal degeneracy (see docstring above). ---
        pole_degeneracy_errors = 0
        for (pole_lat, rkey), states in pole_states.items():
            if len(states) > 1:
                pole_degeneracy_errors += 1
                detail = "; ".join(
                    "state=%d at lon in %s" % (state, sorted(lons))
                    for state, lons in sorted(states.items()))
                errors.append(
                    "%s point %d pole lat=%g R=%g GV: azimuthal degeneracy "
                    "violated -- duplicate longitudes disagree (%s)" %
                    (case["case_id"], point["point_id"], pole_lat, rkey, detail))
                continue
            # All present longitude duplicates agree: publish one canonical
            # sample (lowest longitude, deterministically 0.0 whenever the
            # full duplicate set is present) for downstream FOV/model/
            # asymptotic reductions.
            lon_set = next(iter(states.values()))
            canonical_lon = min(lon_set)
            canonical = sample_map[(point["point_id"], canonical_lon,
                                    pole_lat, rkey)]
            sample_map[(point["point_id"], 0.0, pole_lat, rkey)] = canonical

        expected_samples = len(full_directions) * len(rigidities_gv)
        missing_samples = expected_samples - len(matched_keys)
        coverage_errors = (len(missing_directions) + len(extra_directions) +
                           missing_samples + extra_rigidities)
        if coverage_errors:
            errors.append(
                "%s point %d coverage mismatch: missing_directions=%d "
                "extra_directions=%d missing_samples=%d "
                "extra_rigidity_rows=%d (see C8_coverage_summary.csv)" %
                (case["case_id"], point["point_id"], len(missing_directions),
                 len(extra_directions), missing_samples, extra_rigidities))
        summary = {
            "case_id": case["case_id"], "field_model": case["model"],
            "charge": case["charge"], "point_id": point["point_id"],
            "obs_lon_deg": point["obs_lon_deg"],
            "obs_lat_deg": point["obs_lat_deg"],
            "expected_directions": len(full_directions),
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
            "passed": int(coverage_errors == 0 and pole_degeneracy_errors == 0 and
                          point_counts["contract_errors"] == 0 and
                          point_counts["fatal"] == 0),
        }
        summaries.append(summary)
        total.update(point_counts)
        total["coverage_errors"] += coverage_errors
        total["pole_degeneracy_errors"] += pole_degeneracy_errors

    if len(summaries) != len(points):
        total["coverage_errors"] += len(points) - len(summaries)
    return summaries, sample_map, total, errors


def fov_direction_terms(point, expected_directions, dir_res_deg, args):
    """Yield (direction_key, sector, weight_sr) for every direction that
    contributes to one observation point's local-aperture reduction.

    This walks the regular (non-polar) sky grid exactly as before, PLUS the
    two poles -- each represented ONCE by its canonical (longitude 0.0)
    direction and its exact polar-cap solid angle (``pole_cap_weight_sr``).

    NOTE (fix history): before this fix, the two poles were never evaluated
    here at all; their solid angle was folded, as an approximation, into the
    adjacent regular band by ``sky_cell_weight_sr`` stretching that band's
    edge to 90 degrees.  For the C8 points at |latitude| = 51.6 degrees, one
    pole actually falls inside the 45-degree zenith aperture (its angular
    distance from zenith there is 90 - 51.6 = 38.4 degrees), so it is not a
    negligible edge case: the previous approximation replaced a real traced
    access state at that direction with a copy of the neighboring band's
    state.  This function is the SINGLE shared definition of "which
    directions, sectors, and weights make up a point's FOV/EAST/WEST/CENTER
    reduction," used identically by ``reduce_fov``, ``model_comparison``, and
    ``asymptotic_comparison`` so the three reductions can never disagree
    about the geometry.  Callers still look up the actual sample via
    ``sample_map``/case sample maps keyed by the yielded direction key; for
    poles that key is always ``direction_key(0.0, +/-90.0)``, the canonical
    sample published by ``audit_case`` once gate C8-G11 (pole azimuthal
    degeneracy) has verified every AMPS longitude duplicate agrees.
    """
    for dkey in expected_directions:
        sector, zenith_deg, east_score = local_aperture(
            point, dkey[0], dkey[1], args.fov_half_angle,
            args.east_west_deadband)
        if sector == "OUTSIDE":
            continue
        yield dkey, sector, sky_cell_weight_sr(dir_res_deg, dir_res_deg, dkey[1])
    for pole_lat in (-90.0, 90.0):
        canonical = direction_key(0.0, pole_lat)
        sector, zenith_deg, east_score = local_aperture(
            point, canonical[0], canonical[1], args.fov_half_angle,
            args.east_west_deadband)
        if sector == "OUTSIDE":
            continue
        yield canonical, sector, pole_cap_weight_sr(dir_res_deg)


def reduce_fov(case, points, sample_map, expected_directions, rigidities_gv,
               args):
    """Reduce each cube to full, EAST, WEST, and meridian FOV solid angles."""
    rows = []
    errors = []
    for point in points:
        geometry = {}
        for dkey, sector, weight_sr in fov_direction_terms(
                point, expected_directions, args.dir_res, args):
            geometry[dkey] = {
                "sector": sector, "weight_sr": weight_sr,
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
    """Compare T96/T05 access states on identical resolved FOV samples.

    Uses ``fov_direction_terms`` (shared with ``reduce_fov``) so the set of
    directions compared here -- regular cells plus the two canonical pole
    directions -- always matches what the FOV/East-West reduction actually
    integrated.  Before this fix the two functions could silently disagree,
    since poles were not evaluated by either one.
    """
    rows = []
    messages = []
    passed = True
    for charge in (+1, -1):
        left = cases[("T96", charge)]["sample_map"]
        right = cases[("T05", charge)]["sample_map"]
        for point in points:
            for rigidity in rigidities_gv:
                counts = Counter()
                for dkey, _sector, _weight_sr in fov_direction_terms(
                        point, expected_directions, args.dir_res, args):
                    key = (point["point_id"], dkey[0], dkey[1],
                           rigidity_key(rigidity))
                    a = left.get(key)
                    b = right.get(key)
                    if a is None or b is None:
                        counts["missing"] += 1
                    elif a["access_state"] == 2 and b["access_state"] == 2:
                        counts["both_unresolved"] += 1
                    elif a["access_state"] == 2 or b["access_state"] == 2:
                        counts["one_unresolved"] += 1
                    else:
                        counts["resolved"] += 1
                        counts["mismatch"] += int(
                            a["access_state"] != b["access_state"])
                mismatch_fraction = (float(counts["mismatch"]) /
                                     counts["resolved"]
                                     if counts["resolved"] else 1.0)
                row_passed = (counts["missing"] == 0 and
                              (rigidity < args.high_rigidity_min or
                               mismatch_fraction <= args.max_model_mismatch))
                rows.append({
                    "charge": charge, "point_id": point["point_id"],
                    "obs_lon_deg": point["obs_lon_deg"],
                    "obs_lat_deg": point["obs_lat_deg"],
                    "rigidity_GV": rigidity, "resolved_pairs": counts["resolved"],
                    "state_mismatches": counts["mismatch"],
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
            aggregate["missing"] += row["missing_samples"]
        fraction = (float(aggregate["mismatch"]) / aggregate["resolved"]
                    if aggregate["resolved"] else 1.0)
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
    """Compare optional T96/T05 asymptotic coordinates and compute RMS values.

    Uses ``fov_direction_terms`` (shared with ``reduce_fov`` and
    ``model_comparison``) so this comparison covers exactly the same
    direction set -- regular cells plus the two canonical pole directions --
    as the other two FOV-based reductions.
    """
    rows = []
    bins = defaultdict(list)
    available_samples = 0
    for charge in (+1, -1):
        left = cases[("T96", charge)]["sample_map"]
        right = cases[("T05", charge)]["sample_map"]
        for point in points:
            for dkey, _sector, _weight_sr in fov_direction_terms(
                    point, expected_directions, args.dir_res, args):
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
                    bins[bin_name].append((lat_mrad, lon_mrad))
                    rows.append({
                        "charge": charge, "point_id": point["point_id"],
                        "sky_lon_deg": dkey[0], "sky_lat_deg": dkey[1],
                        "rigidity_GV": rigidity, "rigidity_bin": bin_name,
                        "latitude_difference_mrad": lat_mrad,
                        "longitude_difference_mrad": lon_mrad,
                        "t96_schema": a["asymptotic"][2],
                        "t05_schema": b["asymptotic"][2],
                    })

    summary = []
    for bin_name in ("20_to_30_GV", "gt_50_GV"):
        values = bins.get(bin_name, [])
        lat_rms = (math.sqrt(sum(item[0] ** 2 for item in values) / len(values))
                   if values else None)
        lon_rms = (math.sqrt(sum(item[1] ** 2 for item in values) / len(values))
                   if values else None)
        summary.append({
            "rigidity_bin": bin_name, "samples": len(values),
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
    # 11 gates as of the FIX 2 revision: C8-G11 ("pole azimuthal degeneracy")
    # was added as a strictly ADDITIONAL hard gate -- see the module
    # docstring's "FIX 2" note.  No existing gate's policy or limit changed.
    if len(gates) != 11 or {row["gate_id"] for row in gates} != \
            {"C8-G%02d" % index for index in range(1, 12)}:
        raise RuntimeError("Acceptance contract must contain gates C8-G01..G11")
    with expected.open(newline="") as stream:
        expected_rows = list(csv.DictReader(stream))
    if len(expected_rows) != 6 or {row["reference_id"] for row in expected_rows} != \
            {"C8-R%02d" % index for index in range(1, 7)}:
        raise RuntimeError("Expected-physics reference must contain C8-R01..R06")
    return rows, gates, expected_rows


def write_synthetic_access(path, directions, rigidities, state_function,
                           with_asymptotic=False):
    """Create a tiny valid cube used only by ``--self-test``."""
    variables = ["lon_deg", "lat_deg", "rigidity_GV", "energy_MeV",
                 "access_state", "termination_code", "trace_time_s",
                 "trace_distance_Re", "trace_steps", "retry_count"]
    if with_asymptotic:
        variables += ["asymptotic_lon_deg", "asymptotic_lat_deg"]
    lines = ["VARIABLES=" + " ".join('"%s"' % name for name in variables),
             "ZONE T=\"synthetic\""]
    for dkey in sorted(directions):
        for rigidity in rigidities:
            state = int(state_function(dkey, rigidity))
            code = {0: 1, 1: 0, 2: 3}[state]
            values = [dkey[0], dkey[1], rigidity, rigidity * 1000.0,
                      state, code, 1.0, 2.0, 10, 0]
            if with_asymptotic:
                values += [dkey[0] + 0.01,
                           max(-90.0, min(90.0, dkey[1] + 0.01))]
            lines.append(" ".join(fmt_number(value) for value in values))
    path.write_text("\n".join(lines) + "\n")


def self_test(script_dir):
    """Exercise reference, geometry, rendering, parser, and failure guards."""
    validate_reference_files(script_dir)

    # --- Regular (non-polar) grid: unchanged behavior/counts. -----------
    if len(expected_direction_keys(30.0, 30.0)) != 60:
        raise AssertionError("30-degree regular (non-polar) task count changed")
    if len(expected_direction_keys(15.0, 15.0)) != 264:
        raise AssertionError("15-degree regular (non-polar) task count changed")

    # --- FIX 2 checks: pole grid, full contract, and 4*pi conservation. ---
    if len(expected_pole_direction_keys(15.0, 15.0)) != 48:
        raise AssertionError(
            "15-degree pole duplicate count changed (expected 24 longitudes "
            "x 2 poles = 48, matching AMPS's nLatMap = ... + 1 // include "
            "poles behavior)")
    if len(expected_full_direction_keys(15.0, 15.0)) != 312:
        raise AssertionError(
            "15-degree full (regular+polar) DIRECT_ACCESS task count "
            "changed; this must track AMPS's true producer contract")

    regular_weight = sum(sky_cell_weight_sr(15.0, 15.0, key[1])
                         for key in expected_direction_keys(15.0, 15.0))
    full_weight = regular_weight + 2.0 * pole_cap_weight_sr(15.0)
    if not math.isclose(full_weight, 4.0 * math.pi, rel_tol=0.0, abs_tol=1.0e-12):
        raise AssertionError(
            "Regular-band weights plus the two exact polar-cap weights do "
            "not close to 4*pi sr; the FOV/East-West solid-angle "
            "denominator would be wrong")
    # Independent re-derivation (standard spherical-cap formula
    # 2*pi*(1-cos(half_angle))) guards against a sign/formula slip inside
    # pole_cap_weight_sr itself rather than merely checking it against its
    # own algebraic rearrangement.
    half_angle_rad = math.radians(0.5 * 15.0)
    cross_check_cap = 2.0 * math.pi * (1.0 - math.cos(half_angle_rad))
    if not math.isclose(pole_cap_weight_sr(15.0), cross_check_cap,
                        rel_tol=1.0e-12, abs_tol=1.0e-15):
        raise AssertionError(
            "pole_cap_weight_sr disagrees with the independent spherical-"
            "cap cross-check formula")

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
    with tempfile.TemporaryDirectory(prefix="C8_selftest_") as temp_name:
        temp = Path(temp_name)
        rendered = temp / "AMPS_PARAM_C8.in"
        render_input(script_dir / "AMPS_PARAM_C8_gridless.in", rendered,
                     "C8_selftest", "T05", -1, [point], [1.0, 50.0], args)
        text = rendered.read_text()
        required_patterns = (
            r"^CUTOFF_SAMPLING\s+VERTICAL\s*$",
            r"^CUTOFF_SEARCH_ALGORITHM\s+DIRECT_ACCESS\s*$",
            # FIX 1: C8 must request the physically correct antiparticle
            # backtrace convention.  SAME (the pre-fix setting) numerically
            # swaps which real charge sign each case reports; see the
            # module docstring's "FIX 1" note and AMPS_PARAM_C8_gridless.in.
            r"^CUTOFF_BACKTRACE_CHARGE\s+REVERSED\s*$",
            r"^FIELD_MODEL\s+T05\s*$",
            r"^CHARGE\s+-1\s*$",
        )
        if any(not re.search(pattern, text, re.MULTILINE)
               for pattern in required_patterns) or \
                re.search(r"__[A-Z0-9_]+__", text):
            raise AssertionError("Rendered C8 input is incomplete")
        if re.search(r"^CUTOFF_BACKTRACE_CHARGE\s+SAME\s*$", text, re.MULTILINE):
            raise AssertionError(
                "Rendered C8 input still requests the SAME backtrace-charge "
                "convention; this is the exact defect FIX 1 corrects")

        directions = expected_direction_keys(90.0, 90.0)
        access_path = temp / "cube.dat"
        write_synthetic_access(
            access_path, directions, [1.0, 50.0],
            lambda dkey, rigidity: 1 if rigidity >= 50.0 else 0,
            with_asymptotic=True)
        parsed = parse_directional_access(access_path)
        if len(parsed["samples"]) != len(directions) * 2 or \
                not all(item["asymptotic"] for item in parsed["samples"].values()):
            raise AssertionError("Valid synthetic cube did not parse")

        # Duplicate rows and state/reason disagreement are distinct defects.
        duplicate_path = temp / "duplicate.dat"
        duplicate_path.write_text(access_path.read_text() +
                                  access_path.read_text().splitlines()[-1] + "\n")
        try:
            parse_directional_access(duplicate_path)
            raise AssertionError("Duplicate sample was accepted")
        except RuntimeError as exc:
            if "Duplicate" not in str(exc):
                raise
        broken_text = access_path.read_text().replace(" 1 0 1 2 10 0 ",
                                                      " 1 1 1 2 10 0 ", 1)
        broken_path = temp / "broken.dat"
        broken_path.write_text(broken_text)
        broken = parse_directional_access(broken_path)
        if not any(item["contract_error"] for item in broken["samples"].values()):
            raise AssertionError("State/reason mismatch was not detected")

        # ------------------------------------------------------------------
        # FIX 2, part (a)+(b): full-grid completeness and gate C8-G11 (pole
        # azimuthal degeneracy), exercised end-to-end through audit_case.
        # A coarse 90-degree grid keeps this fast: 4 regular directions (all
        # at latitude 0) plus 4 longitude-tagged duplicates at each pole.
        # ------------------------------------------------------------------
        pole_res = 90.0
        pole_regular = expected_direction_keys(pole_res, pole_res)
        pole_poles = expected_pole_direction_keys(pole_res, pole_res)
        pole_full = pole_regular | pole_poles
        if len(pole_regular) != 4 or len(pole_poles) != 8 or len(pole_full) != 12:
            raise AssertionError("90-degree pole self-test grid sizes changed")

        clean_dir = temp / "pole_clean"
        clean_dir.mkdir()
        clean_access = clean_dir / "cutoff_gridless_dir_access_point_0000.dat"
        write_synthetic_access(clean_access, pole_full, [10.0],
                               lambda dkey, rigidity: 1)  # everything ALLOWED
        clean_case = {"case_id": "selftest_pole_clean", "model": "T05",
                     "charge": -1, "workdir": clean_dir}
        _, clean_sample_map, clean_totals, clean_errors = audit_case(
            clean_case, [point], pole_regular, pole_poles, [10.0])
        if clean_totals["coverage_errors"] or clean_totals["pole_degeneracy_errors"] \
                or clean_errors:
            raise AssertionError(
                "Self-consistent full (regular+polar) cube unexpectedly "
                "failed the coverage/pole-degeneracy audit: %r" % clean_errors)
        rkey10 = rigidity_key(10.0)
        if (0, 0.0, -90.0, rkey10) not in clean_sample_map or \
                (0, 0.0, 90.0, rkey10) not in clean_sample_map:
            raise AssertionError(
                "audit_case did not publish canonical (lon=0) pole samples "
                "for downstream FOV/model/asymptotic reductions")

        # Corrupt exactly one south-pole longitude duplicate's access state
        # so the duplicates at that pole/rigidity no longer agree.
        lines = clean_access.read_text().splitlines()
        target = next((i for i, line in enumerate(lines)
                       if line.split()[:2] == [fmt_number(0.0), fmt_number(-90.0)]),
                      None)
        if target is None:
            raise AssertionError(
                "Self-test could not locate a south-pole row to corrupt")
        parts = lines[target].split()
        parts[4], parts[5] = "0", "1"  # ALLOWED(1)/code0 -> FORBIDDEN(0)/code1
        lines[target] = " ".join(parts)
        broken_pole_dir = temp / "pole_broken"
        broken_pole_dir.mkdir()
        broken_pole_access = (broken_pole_dir /
                              "cutoff_gridless_dir_access_point_0000.dat")
        broken_pole_access.write_text("\n".join(lines) + "\n")
        broken_pole_case = {"case_id": "selftest_pole_broken", "model": "T05",
                           "charge": -1, "workdir": broken_pole_dir}
        _, _, broken_pole_totals, broken_pole_errors = audit_case(
            broken_pole_case, [point], pole_regular, pole_poles, [10.0])
        if broken_pole_totals["pole_degeneracy_errors"] != 1:
            raise AssertionError(
                "Gate C8-G11 did not detect a single corrupted "
                "pole-duplicate row (found %d)" %
                broken_pole_totals["pole_degeneracy_errors"])
        if not any("azimuthal degeneracy" in item for item in broken_pole_errors):
            raise AssertionError(
                "Pole-degeneracy violation (gate C8-G11) produced no "
                "explanatory message")

        # A cube that omits the poles entirely (the pre-fix producer
        # assumption) must still be rejected: the completeness gate became
        # more ACCURATE, not more permissive.
        truncated_dir = temp / "pole_truncated"
        truncated_dir.mkdir()
        truncated_access = (truncated_dir /
                            "cutoff_gridless_dir_access_point_0000.dat")
        write_synthetic_access(truncated_access, pole_regular, [10.0],
                               lambda dkey, rigidity: 1)
        truncated_case = {"case_id": "selftest_pole_truncated", "model": "T05",
                         "charge": -1, "workdir": truncated_dir}
        _, _, truncated_totals, truncated_errors = audit_case(
            truncated_case, [point], pole_regular, pole_poles, [10.0])
        if truncated_totals["coverage_errors"] == 0:
            raise AssertionError(
                "A cube missing the pole rows AMPS actually writes was "
                "incorrectly accepted as complete -- gate C8-G01 was "
                "weakened, not fixed")
        if not any("coverage mismatch" in item for item in truncated_errors):
            raise AssertionError(
                "Missing-pole coverage failure (gate C8-G01) produced no "
                "explanatory message -- the silent-failure bug is back")

    # ------------------------------------------------------------------
    # FIX 1 check: gate C8-G07 (East-West charge-odd response) must accept
    # a physically correct sign pattern and reject the exact sign-inverted
    # pattern that CUTOFF_BACKTRACE_CHARGE SAME used to produce.  This does
    # not touch AMPS at all -- it proves the acceptance ARITHMETIC was
    # always correct and unchanged; only the upstream input convention
    # (fixed in AMPS_PARAM_C8_gridless.in) needed to change.
    # ------------------------------------------------------------------
    class EwArgs(object):
        pass
    ew_args = EwArgs()
    ew_args.transition_min = 0.05
    ew_args.transition_max = 0.95
    ew_args.max_lobe_unresolved = 0.20
    ew_args.min_east_west_effect = 0.01

    def make_ew_row(model, charge, rigidity, west_allowed, east_allowed):
        return {
            "field_model": model, "charge": charge, "rigidity_GV": rigidity,
            "west_allowed_fraction": west_allowed,
            "east_allowed_fraction": east_allowed,
            "fov_allowed_fraction": 0.5 * (west_allowed + east_allowed),
            "east_unresolved_fraction": 0.0, "west_unresolved_fraction": 0.0,
        }

    # Textbook East-West effect: positive charges see MORE access from the
    # west (T_WEST > T_EAST); negative charges see the opposite.
    correct_rows = [make_ew_row(model, charge, 10.0, west, east)
                    for model in ("T96", "T05")
                    for charge, west, east in ((+1, 0.60, 0.40),
                                               (-1, 0.40, 0.60))]
    _, correct_passed, correct_messages = east_west_comparison(correct_rows, ew_args)
    if not correct_passed or correct_messages:
        raise AssertionError(
            "Gate C8-G07 rejected a physically correct East-West "
            "asymmetry: %r" % correct_messages)

    # The exact defect CUTOFF_BACKTRACE_CHARGE SAME used to produce: both
    # charges' D(q) = T_WEST - T_EAST signs inverted from the textbook
    # pattern above.
    reversed_rows = [make_ew_row(model, charge, 10.0, west, east)
                     for model in ("T96", "T05")
                     for charge, west, east in ((+1, 0.40, 0.60),
                                                (-1, 0.60, 0.40))]
    _, reversed_passed, reversed_messages = east_west_comparison(
        reversed_rows, ew_args)
    if reversed_passed or not reversed_messages:
        raise AssertionError(
            "Gate C8-G07 failed to reject a sign-inverted charge-odd "
            "East-West response (the SAME-vs-REVERSED defect fixed in "
            "AMPS_PARAM_C8_gridless.in); the acceptance criterion may have "
            "been relaxed")

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
    # ``expected_directions`` remains the regular (non-polar) grid used for
    # solid-angle-weighted reductions.  ``pole_directions`` is the set of
    # per-longitude pole duplicates AMPS actually writes, and their union
    # (``full_directions``) is the true output contract audited for gate
    # C8-G01.  See the module docstring's "FIX 2" note.
    expected_directions = expected_direction_keys(args.dir_res, args.dir_res)
    pole_directions = expected_pole_direction_keys(args.dir_res, args.dir_res)
    full_directions = expected_directions | pole_directions

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
            continue
        summaries, sample_map, totals, audit_messages = audit_case(
            case, points, expected_directions, pole_directions, rigidities)
        case["sample_map"] = sample_map
        case["audit_totals"] = dict(totals)
        coverage_rows.extend(summaries)
        messages.extend(audit_messages)
        # NOTE (fix history): every one of these hard-gate counters now has
        # an accompanying entry in ``audit_messages`` (appended inside
        # ``audit_case``), so a nonzero counter can no longer flip
        # ``overall_passed`` to False without a corresponding, specific
        # message reaching C8_summary.txt.  Previously ``coverage_errors``
        # alone could fail the run silently.
        if audit_messages or totals["coverage_errors"] or \
                totals["pole_degeneracy_errors"] or \
                totals["contract_errors"] or totals["fatal"]:
            overall_passed = False
        unresolved_fraction = (float(totals["state_2"]) / totals["samples"]
                               if totals["samples"] else 1.0)
        case["unresolved_fraction"] = unresolved_fraction
        if unresolved_fraction > args.max_unresolved_fraction:
            overall_passed = False
            messages.append("%s unresolved fraction %.3f exceeds %.3f" %
                            (case["case_id"], unresolved_fraction,
                             args.max_unresolved_fraction))
        reduced, reduction_errors = reduce_fov(
            case, points, sample_map, expected_directions, rigidities, args)
        fov_rows.extend(reduced)
        messages.extend(reduction_errors)
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

    if all(cases[key].get("sample_map") for key in cases):
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
                          "longitude_difference_mrad", "t96_schema", "t05_schema"))
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
        # requested_direction_tasks_per_point: regular (non-polar) grid size,
        # unchanged in meaning from schema version 1 -- this is the grid used
        # for solid-angle-weighted FOV/East-West reductions.
        "requested_direction_tasks_per_point": len(expected_directions),
        # New in schema version 2: the pole duplicates AMPS actually writes,
        # and the true total DIRECT_ACCESS output contract (their union).
        # See the module docstring's "FIX 2" note and gate C8-G01/C8-G11.
        "requested_pole_direction_tasks_per_point": len(pole_directions),
        "requested_full_direction_tasks_per_point": len(full_directions),
        "fov_half_angle_deg": args.fov_half_angle,
        "asymptotic_policy": args.asymptotic_policy,
        "asymptotic_columns_available": asym_available,
        "asymptotic_summary": asym_summary,
        "cases": [{
            "case_id": case["case_id"], "model": case["model"],
            "charge": case["charge"], "return_code": case.get("return_code", 0),
            "unresolved_fraction": case.get("unresolved_fraction"),
        } for case in cases.values()],
        "messages": messages, "runs": run_records,
    }
    (base_workdir / "C8_result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n")
    summary_lines = [
        "C8 — %s" % TEST_NAME,
        "result: %s" % ("PASS" if overall_passed else "FAIL"),
        "profile: %s; cases: 4; points: %d; directions/point: %d regular + "
        "%d polar = %d total; rigidities: %d" %
        (args.profile, len(points), len(expected_directions),
         len(pole_directions), len(full_directions), len(rigidities)),
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

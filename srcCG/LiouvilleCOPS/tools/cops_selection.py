#!/usr/bin/env python3
"""Versioned observation-selection profiles for the LiouvilleCOPS fitter.

This module contains selection policy only.  It does not evaluate the gas
model, alter pressures, or choose fitting parameters.  Both ``fit_cops.py``
and ``run_fit_profile.py`` use it so named profiles cannot acquire subtly
different timestamp or validity rules in different entry points.

Profile intervals are closed at both ends.  A datum is fundamentally usable
only when its reference pressure is positive, the gauge is outside the shape,
and its timestamp is inside any FitInterval declared by the base run.  A
profile can then apply the base run's gauge exclusions, explicit intervals,
and the independently supplied saturated-NG mask.
"""

import datetime as dt
import json
from collections import Counter
from pathlib import Path

from cops_io import canonical_sha256, file_identity


PROFILE_SCHEMA = "liouville_cops.selection_profile"
PROFILE_SCHEMA_VERSION = 1
GAUGES = ("NG", "RG")
SATURATION_POLICIES = ("include", "exclude")


class SelectionError(ValueError):
    """A selection profile or its application is invalid."""


def parse_utc(value):
    """Return seconds after 2016-01-01 for ``YYYY-MM-DDTh:mm:ss[.sss]``.

    LiouvilleCOPS tables historically contain hours without a leading zero, so
    this parser deliberately handles the time fields numerically instead of
    requiring ISO zero padding.
    """
    try:
        date, clock = value.strip().rstrip("Z").split("T")
        hour, minute, second = clock.split(":")
        midnight = dt.datetime.strptime(date, "%Y-%m-%d")
        return ((midnight - dt.datetime(2016, 1, 1)).total_seconds()
                + 3600 * int(hour) + 60 * int(minute) + float(second))
    except (AttributeError, TypeError, ValueError):
        raise SelectionError("invalid UTC timestamp: %r" % value) from None


def profile_directory(app_root):
    """Return the application's checked-in selection-profile directory."""
    return Path(app_root) / "fit_profiles"


def available_profiles(directory):
    """Return sorted profile names advertised by valid-looking JSON files."""
    directory = Path(directory)
    return sorted(path.stem for path in directory.glob("*.json") if path.is_file())


def _validate_interval(value, label):
    if not isinstance(value, list) or len(value) != 2:
        raise SelectionError("%s must be a two-element UTC array" % label)
    start, stop = parse_utc(value[0]), parse_utc(value[1])
    if stop < start:
        raise SelectionError("%s ends before it starts" % label)
    return (start, stop, str(value[0]), str(value[1]))


def validate_profile(profile, source="profile"):
    """Validate and normalize a decoded profile without changing its meaning."""
    if not isinstance(profile, dict):
        raise SelectionError("%s is not a JSON object" % source)
    if profile.get("schema") != PROFILE_SCHEMA:
        raise SelectionError("%s has unexpected schema %r" %
                             (source, profile.get("schema")))
    if profile.get("schema_version") != PROFILE_SCHEMA_VERSION:
        raise SelectionError("%s has unsupported schema version %r" %
                             (source, profile.get("schema_version")))
    name = profile.get("name")
    if not isinstance(name, str) or not name or Path(name).name != name:
        raise SelectionError("%s has an invalid profile name" % source)
    if profile.get("candidate_policy") != "positive_outside_shape":
        raise SelectionError("%s has unsupported candidate_policy %r" %
                             (source, profile.get("candidate_policy")))
    if not isinstance(profile.get("apply_fit_interval"), bool):
        raise SelectionError("%s: apply_fit_interval must be boolean" % source)
    if not isinstance(profile.get("apply_run_exclusions"), bool):
        raise SelectionError("%s: apply_run_exclusions must be boolean" % source)
    if profile.get("saturated_ng") not in SATURATION_POLICIES:
        raise SelectionError("%s: saturated_ng must be include or exclude" % source)

    explicit = profile.get("exclude_intervals", {})
    if not isinstance(explicit, dict) or set(explicit) - set(GAUGES):
        raise SelectionError("%s: exclude_intervals may contain only NG and RG" % source)
    normalized = {}
    for gauge in GAUGES:
        values = explicit.get(gauge, [])
        if not isinstance(values, list):
            raise SelectionError("%s: exclude_intervals.%s must be an array" %
                                 (source, gauge))
        normalized[gauge] = [
            _validate_interval(value, "%s exclude_intervals.%s[%d]" %
                               (source, gauge, index))
            for index, value in enumerate(values)
        ]

    expected = profile.get("expected_counts")
    if not isinstance(expected, dict):
        raise SelectionError("%s: expected_counts is required" % source)
    for key in ("rows", "used_ng", "used_rg"):
        value = expected.get(key)
        if not isinstance(value, int) or value < 0:
            raise SelectionError("%s: expected_counts.%s must be a nonnegative integer" %
                                 (source, key))
    result = dict(profile)
    result["_parsed_exclude_intervals"] = normalized
    return result


def load_profile(name, directory):
    """Load ``<directory>/<name>.json`` and return profile plus file identity."""
    if not isinstance(name, str) or not name or Path(name).name != name:
        raise SelectionError("invalid profile name: %r" % name)
    path = Path(directory) / (name + ".json")
    try:
        profile = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        names = available_profiles(directory)
        raise SelectionError("cannot read selection profile %s: %s; available: %s" %
                             (path, exc, ", ".join(names) or "none")) from None
    except ValueError as exc:
        raise SelectionError("invalid JSON in selection profile %s: %s" %
                             (path, exc)) from None
    profile = validate_profile(profile, str(path))
    if profile["name"] != name:
        raise SelectionError("profile filename %s does not match its name %s" %
                             (name, profile["name"]))
    return profile, file_identity(path, "selection_profile")


def _in_intervals(timestamp, intervals):
    """Return true for closed intervals; interval tuples may include labels."""
    return any(interval[0] <= timestamp <= interval[1] for interval in intervals)


def apply_selection_profile(rows, columns, profile, run_intervals, saturated_ng):
    """Apply one validated profile and return copied rows plus an audit report.

    ``run_intervals`` is the dictionary returned by ``fit_cops.run_intervals``.
    ``saturated_ng`` must have one boolean per row.  The original rows are not
    mutated.  Count mismatches are fatal because the checked-in profiles define
    the currently audited 1,150-row mission dataset.
    """
    profile = validate_profile(profile, profile.get("name", "profile"))
    required = ("UTC", "inside", "refNG", "refRG", "useNG", "useRG")
    missing = [name for name in required if name not in columns]
    if missing:
        raise SelectionError("model table is missing selection columns: %s" %
                             ", ".join(missing))
    if len(saturated_ng) != len(rows):
        raise SelectionError("saturated-NG mask has %d entries for %d rows" %
                             (len(saturated_ng), len(rows)))

    parsed_time = [parse_utc(row[columns["UTC"]]) for row in rows]
    fit_intervals = run_intervals.get("FitInterval", [])
    result = [list(row) for row in rows]
    reasons = {gauge: [] for gauge in GAUGES}

    for index, (row, timestamp) in enumerate(zip(result, parsed_time)):
        inside = bool(int(row[columns["inside"]]))
        in_fit = (not profile["apply_fit_interval"] or not fit_intervals
                  or _in_intervals(timestamp, fit_intervals))
        for gauge in GAUGES:
            pressure = float(row[columns["ref" + gauge]])
            reason = "selected"
            if not pressure > 0.0:
                reason = "nonpositive_or_missing"
            elif inside:
                reason = "inside_shape"
            elif not in_fit:
                reason = "outside_fit_interval"
            elif (profile["apply_run_exclusions"] and
                  _in_intervals(timestamp, run_intervals.get("Exclude" + gauge, []))):
                reason = "run_exclusion"
            elif _in_intervals(timestamp,
                               profile["_parsed_exclude_intervals"][gauge]):
                reason = "profile_exclusion"
            elif (gauge == "NG" and saturated_ng[index] and
                  profile["saturated_ng"] == "exclude"):
                reason = "saturated_excluded"
            row[columns["use" + gauge]] = "1" if reason == "selected" else "0"
            reasons[gauge].append(reason)

    counts = {
        "rows": len(result),
        "used_ng": sum(int(row[columns["useNG"]]) for row in result),
        "used_rg": sum(int(row[columns["useRG"]]) for row in result),
        "saturated_ng_total": sum(bool(value) for value in saturated_ng),
        "saturated_ng_selected": sum(
            bool(value) and result[index][columns["useNG"]] == "1"
            for index, value in enumerate(saturated_ng)),
    }
    expected = profile["expected_counts"]
    mismatch = [key for key in ("rows", "used_ng", "used_rg")
                if counts[key] != expected[key]]
    if mismatch:
        detail = ", ".join("%s=%d (expected %d)" %
                           (key, counts[key], expected[key]) for key in mismatch)
        raise SelectionError(
            "selection profile %s does not match the audited dataset: %s; "
            "review the changed inputs before updating expected_counts" %
            (profile["name"], detail))

    reason_counts = {
        gauge: dict(sorted(Counter(reasons[gauge]).items())) for gauge in GAUGES
    }
    reason_rows = [[str(row[columns["UTC"]]), reasons["NG"][index],
                    reasons["RG"][index]]
                   for index, row in enumerate(result)]
    report = {
        "profile": profile["name"],
        "description": profile.get("description", ""),
        "counts": counts,
        "reason_counts": reason_counts,
        "reason_sha256": canonical_sha256(reason_rows),
        "reasons": reasons,
        "interval_semantics": "closed_start_and_stop",
    }
    return result, report


def public_profile(profile):
    """Remove parsed implementation fields before writing effective JSON."""
    return {key: value for key, value in profile.items() if not key.startswith("_")}


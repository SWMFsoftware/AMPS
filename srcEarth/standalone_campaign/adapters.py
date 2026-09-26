#!/usr/bin/env python3
"""Observation adapters, exposure averaging, and comparison metrics.

Each adapter converts a released instrument product into one canonical record
shape.  The conversion is intentionally explicit and conservative: missing or
bad-quality records are excluded with counts, units are recorded, and no
instrument- or platform-specific multiplicative normalization is permitted.

The cadence routine integrates a piecewise-linear model prediction over the
actual observation accumulation window.  It never substitutes the nearest
instantaneous snapshot for a finite exposure and never extrapolates beyond the
available model samples.
"""

from __future__ import annotations

import csv
import gzip
import math
import statistics
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .campaign import CampaignError, format_utc, parse_utc


def _truth(value: object) -> bool:
    return str(value).strip().upper() in ("1", "T", "TRUE", "Y", "YES")


def _finite(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise CampaignError("%s is not numeric: %s" % (label, value)) from exc
    if not math.isfinite(result):
        raise CampaignError("%s is not finite" % label)
    return result


def _read_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """Read plain or gzip CSV while retaining provenance comments."""

    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(str(path), "rt", encoding="utf-8", newline="") as stream:
        lines = [line for line in stream if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    fields = list(reader.fieldnames or [])
    rows = list(reader)
    if not fields:
        raise CampaignError("observation file has no CSV header: %s" % path)
    return fields, rows


def _require_columns(path: Path, fields: Sequence[str], required: Iterable[str]) -> None:
    missing = sorted(set(required).difference(fields))
    if missing:
        raise CampaignError("%s lacks columns: %s" % (path, ", ".join(missing)))


def _window_from_midpoint(utc: str, cadence_seconds: int) -> Tuple[str, str, str]:
    if cadence_seconds <= 0:
        raise CampaignError("instrument cadence_seconds must be positive")
    midpoint = parse_utc(utc)
    half = timedelta(seconds=0.5 * cadence_seconds)
    return format_utc(midpoint - half), format_utc(midpoint + half), format_utc(midpoint)


def _record(
    instrument_id: str,
    start: str,
    end: str,
    midpoint: str,
    platform: str,
    instrument: str,
    channel: str,
    direction: str,
    quantity: str,
    units: str,
    value: float,
    lower: Optional[float],
    upper: Optional[float],
    source_row: int,
) -> Dict[str, object]:
    start_time = parse_utc(start)
    end_time = parse_utc(end)
    midpoint_time = parse_utc(midpoint)
    if end_time < start_time or not start_time <= midpoint_time <= end_time:
        raise CampaignError("observation window/midpoint ordering is invalid")
    if not all(
        str(value).strip()
        for value in (
            instrument_id, platform, instrument, channel, direction, quantity, units
        )
    ):
        raise CampaignError("observation identity, quantity, and units must be explicit")
    if not math.isfinite(value):
        raise CampaignError("observation value is not finite")
    if lower is not None and (not math.isfinite(lower) or lower > value):
        raise CampaignError("observation lower bound is invalid")
    if upper is not None and (not math.isfinite(upper) or upper < value):
        raise CampaignError("observation upper bound is invalid")
    if lower is not None and upper is not None and lower > upper:
        raise CampaignError("observation bounds are reversed")
    return {
        "instrument_id": instrument_id,
        "start_utc": start,
        "end_utc": end,
        "midpoint_utc": midpoint,
        "platform": platform,
        "instrument": instrument,
        "channel": channel,
        "direction": direction,
        "quantity": quantity,
        "units": units,
        "value": value,
        "lower": lower,
        "upper": upper,
        "quality": "VALID",
        "source_row": source_row,
    }


def _declared_units(config: Mapping[str, object], expected: str) -> str:
    """Require the manifest and adapter to agree on one exact unit string."""

    declared = str(config.get("units", "")).strip()
    if declared != expected:
        raise CampaignError(
            "adapter units mismatch: manifest declares %s, adapter requires %s"
            % (declared or "<empty>", expected)
        )
    return declared


def adapt_pamela(path: Path, config: Mapping[str, object]) -> List[Dict[str, object]]:
    """Adapt the published PAMELA cutoff-latitude interval table used by C9."""

    fields, source = _read_rows(path)
    required = {
        "interval_midpoint_utc", "interval_start_utc", "interval_end_utc",
        "rigidity_min_gv", "rigidity_max_gv", "pamela_cutoff_aacgm_deg",
        "sigma_plus_deg", "sigma_minus_deg", "missing",
    }
    _require_columns(path, fields, required)
    result: List[Dict[str, object]] = []
    identifier = str(config["id"])
    units = _declared_units(config, "deg_AACGM")
    for line_number, row in enumerate(source, start=2):
        if _truth(row["missing"]):
            continue
        value = _finite(row["pamela_cutoff_aacgm_deg"], "PAMELA cutoff")
        sigma_minus = _finite(row["sigma_minus_deg"], "PAMELA sigma_minus")
        sigma_plus = _finite(row["sigma_plus_deg"], "PAMELA sigma_plus")
        rmin = _finite(row["rigidity_min_gv"], "PAMELA rigidity_min")
        rmax = _finite(row["rigidity_max_gv"], "PAMELA rigidity_max")
        result.append(_record(
            identifier, format_utc(parse_utc(row["interval_start_utc"])),
            format_utc(parse_utc(row["interval_end_utc"])),
            format_utc(parse_utc(row["interval_midpoint_utc"])), "PAMELA", "PAMELA",
            "%.12g-%.12g_GV" % (rmin, rmax), "OMNIDIRECTIONAL",
            "cutoff_latitude", units, value, value - sigma_minus,
            value + sigma_plus, line_number,
        ))
    return result


def adapt_poes_metop(path: Path, config: Mapping[str, object]) -> List[Dict[str, object]]:
    """Adapt the archive-derived POES/MetOp MEPED T50 boundary used by C10."""

    fields, source = _read_rows(path)
    required = {
        "interval_midpoint_utc", "interval_start_utc", "interval_end_utc",
        "rigidity_gv", "channel", "hemisphere", "boundary_aacgm_lat_deg",
        "sigma_deg", "acceptance_eligible", "missing", "quality_status",
    }
    _require_columns(path, fields, required)
    result: List[Dict[str, object]] = []
    identifier = str(config["id"])
    units = _declared_units(config, "deg_AACGM")
    for line_number, row in enumerate(source, start=2):
        if _truth(row["missing"]) or not _truth(row["acceptance_eligible"]):
            continue
        # C10 names the archive-derived, independent-window rows
        # PRIMARY_ACCEPTANCE.  PRIMARY_PLOT_ONLY and diagnostic rows are
        # deliberately excluded from the release comparison.
        if row["quality_status"].strip().upper() not in (
            "VALID", "PASS", "ELIGIBLE", "PRIMARY_ACCEPTANCE"
        ):
            continue
        value = _finite(row["boundary_aacgm_lat_deg"], "MEPED boundary")
        sigma = _finite(row["sigma_deg"], "MEPED sigma")
        rigidity = _finite(row["rigidity_gv"], "MEPED rigidity")
        platform = row.get("satellites", "POES_METOP").strip() or "POES_METOP"
        result.append(_record(
            identifier, format_utc(parse_utc(row["interval_start_utc"])),
            format_utc(parse_utc(row["interval_end_utc"])),
            format_utc(parse_utc(row["interval_midpoint_utc"])), platform,
            "MEPED", "%s_%.12g_GV" % (row["channel"].strip(), rigidity),
            row["hemisphere"].strip().upper(), "cutoff_latitude", units,
            value, value - sigma, value + sigma, line_number,
        ))
    return result


def adapt_goes_directional(path: Path, config: Mapping[str, object]) -> List[Dict[str, object]]:
    """Adapt GOES EPEAD physical east/west ratios from the C19 reference.

    C19 has already applied the spacecraft-specific telemetry-head mapping.  The
    adapter consumes the physical east/west columns and therefore never guesses
    a look direction from an ``E`` or ``W`` telemetry suffix.
    """

    fields, source = _read_rows(path)
    required = {
        "utc", "spacecraft", "channel", "east_west_ratio",
        "east_quality_flag", "west_quality_flag", "quality_status",
    }
    _require_columns(path, fields, required)
    cadence = int(config.get("cadence_seconds", 300))
    result: List[Dict[str, object]] = []
    identifier = str(config["id"])
    units = _declared_units(config, "1")
    for line_number, row in enumerate(source, start=2):
        if row["quality_status"].strip().upper() != "VALID":
            continue
        try:
            if int(row["east_quality_flag"]) != 0 or int(row["west_quality_flag"]) != 0:
                continue
        except ValueError as exc:
            raise CampaignError(
                "GOES quality flags must be integers at row %d" % line_number
            ) from exc
        value = _finite(row["east_west_ratio"], "GOES east/west ratio")
        if value <= 0.0:
            continue
        start, end, midpoint = _window_from_midpoint(row["utc"], cadence)
        result.append(_record(
            identifier, start, end, midpoint, row["spacecraft"].strip(), "EPEAD",
            row["channel"].strip(), "EAST_OVER_WEST", "directional_ratio", units,
            value, None, None, line_number,
        ))
    return result


def _first_present(row: Mapping[str, str], names: Sequence[str], label: str) -> str:
    for name in names:
        if name in row and row[name].strip():
            return row[name].strip()
    raise CampaignError("REPT row lacks %s (%s)" % (label, ", ".join(names)))


def adapt_rept(path: Path, config: Mapping[str, object]) -> List[Dict[str, object]]:
    """Adapt a released REPT proton differential-flux export.

    Public REPT exports use several harmless header variants.  The accepted
    aliases below are explicit; an unrecognized schema fails instead of shifting
    a neighboring electron or energy column into the proton comparison.
    """

    fields, source = _read_rows(path)
    if not any(name in fields for name in ("utc", "time_utc", "midpoint_utc")):
        raise CampaignError("%s lacks a REPT UTC column" % path)
    if not any(name in fields for name in ("proton_flux", "flux", "differential_flux")):
        raise CampaignError("%s lacks a REPT proton-flux column" % path)
    cadence = int(config.get("cadence_seconds", 300))
    result: List[Dict[str, object]] = []
    identifier = str(config["id"])
    declared_units = str(config.get("units", "")).strip()
    if not declared_units:
        raise CampaignError("REPT adapter requires explicit manifest units")
    for line_number, row in enumerate(source, start=2):
        quality = row.get("quality_status", row.get("quality", "VALID")).strip().upper()
        if quality not in ("VALID", "GOOD", "0", "PASS"):
            continue
        utc = _first_present(row, ("utc", "time_utc", "midpoint_utc"), "UTC")
        start, end, midpoint = _window_from_midpoint(utc, cadence)
        value = _finite(
            _first_present(row, ("proton_flux", "flux", "differential_flux"), "flux"),
            "REPT proton flux",
        )
        if value <= 0.0:
            continue
        energy = _first_present(
            row, ("energy_mev", "energy_center_mev", "channel"), "energy/channel"
        )
        platform = row.get("spacecraft", row.get("platform", "VAN_ALLEN_PROBES")).strip()
        row_units = row.get("units", declared_units).strip()
        if row_units != declared_units:
            raise CampaignError(
                "REPT units mismatch at row %d: %s != %s"
                % (line_number, row_units, declared_units)
            )
        result.append(_record(
            identifier, start, end, midpoint, platform, "REPT", energy,
            "OMNIDIRECTIONAL", "differential_flux",
            declared_units, value,
            _finite(row["lower"], "REPT lower") if row.get("lower", "").strip() else None,
            _finite(row["upper"], "REPT upper") if row.get("upper", "").strip() else None,
            line_number,
        ))
    return result


ADAPTERS = {
    "PAMELA_CUTOFF": adapt_pamela,
    "POES_METOP_MEPED_CUTOFF": adapt_poes_metop,
    "GOES_EPEAD_DIRECTIONAL": adapt_goes_directional,
    "REPT_PROTON_SPECTRUM": adapt_rept,
}


def validate_instrument_response(
    path: Path, instrument: Mapping[str, object]
) -> Dict[str, object]:
    """Validate the released response/bin definition used by an adapter.

    Cutoff instruments expose rigidity/threshold bins rather than a count-rate
    matrix; energetic-particle telescopes expose energy intervals plus a response
    or effective area.  Both are detector adapters in the campaign sense, and
    both are checked for finite, ordered, positive coordinates before AMPS runs.
    """

    fields, rows = _read_rows(path)
    adapter = str(instrument["adapter"]).upper()
    if adapter == "PAMELA_CUTOFF":
        alternatives = (
            ("rigidity_min_gv", "rigidity_max_gv"),
            ("energy_min_mev", "energy_max_mev"),
        )
    elif adapter == "POES_METOP_MEPED_CUTOFF":
        alternatives = (
            ("energy_threshold_mev", "rigidity_gv"),
            ("energy_min_mev", "energy_max_mev"),
        )
    else:
        alternatives = (("energy_min_mev", "energy_max_mev"),)
    coordinate_pair = next(
        (pair for pair in alternatives if pair[0] in fields and pair[1] in fields), None
    )
    if coordinate_pair is None:
        raise CampaignError(
            "%s lacks a supported %s response/bin coordinate pair"
            % (path, adapter)
        )

    response_fields = (
        "relative_response", "response", "effective_area_cm2_sr",
        "geometric_factor_m2_sr", "rigidity_gv", "rigidity_geometric_center_gv",
    )
    response_field = next((name for name in response_fields if name in fields), None)
    if response_field is None:
        raise CampaignError("%s has no response, area, or rigidity column" % path)
    if not rows:
        raise CampaignError("instrument response is empty: %s" % path)

    for line_number, row in enumerate(rows, start=2):
        left = _finite(row[coordinate_pair[0]], "%s coordinate" % adapter)
        right = _finite(row[coordinate_pair[1]], "%s coordinate" % adapter)
        # POES threshold/rigidity pairs are two different coordinates, not an
        # interval, so only positivity is applicable to that schema.
        if coordinate_pair == ("energy_threshold_mev", "rigidity_gv"):
            if left <= 0.0 or right <= 0.0:
                raise CampaignError("%s:%d nonpositive response coordinate" % (path, line_number))
        elif left < 0.0 or right <= left:
            raise CampaignError("%s:%d invalid response interval" % (path, line_number))
        response = _finite(row[response_field], "%s response" % adapter)
        if response <= 0.0:
            raise CampaignError("%s:%d response must be positive" % (path, line_number))
    return {
        "adapter": adapter,
        "row_count": len(rows),
        "coordinate_columns": list(coordinate_pair),
        "response_column": response_field,
        "status": "PASS",
    }


def adapt_observations(
    path: Path, instrument: Mapping[str, object]
) -> List[Dict[str, object]]:
    """Dispatch one manifest instrument through its named strict adapter."""

    adapter = str(instrument["adapter"]).upper()
    try:
        rows = ADAPTERS[adapter](path, instrument)
    except KeyError as exc:
        raise CampaignError("unsupported observation adapter: %s" % adapter) from exc
    if not rows:
        raise CampaignError("adapter %s produced no quality-controlled rows" % adapter)
    return rows


def write_canonical_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    """Write canonical rows with a stable schema for downstream review."""

    fields = [
        "instrument_id", "start_utc", "end_utc", "midpoint_utc", "platform",
        "instrument", "channel", "direction", "quantity", "units", "value",
        "lower", "upper", "quality", "source_row",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _interpolate(samples: Sequence[Tuple[datetime, float]], when: datetime) -> float:
    for time, value in samples:
        if time == when:
            return value
    for (left_time, left_value), (right_time, right_value) in zip(samples, samples[1:]):
        if left_time < when < right_time:
            fraction = (when - left_time).total_seconds() / (
                right_time - left_time
            ).total_seconds()
            return left_value + fraction * (right_value - left_value)
    raise CampaignError("model series does not bracket %s" % format_utc(when))


def cadence_average(
    samples: Sequence[Tuple[datetime, float]], start: datetime, end: datetime
) -> float:
    """Average a piecewise-linear series over one observation exposure.

    For a zero-duration window the exact interpolated value is returned.  For a
    finite window both endpoints must be bracketed; extrapolation is forbidden.
    Trapezoidal integration is exact for the declared piecewise-linear model.
    """

    ordered = sorted(samples, key=lambda item: item[0])
    if not ordered:
        raise CampaignError("cannot average an empty model series")
    for _, value in ordered:
        if not math.isfinite(value):
            raise CampaignError("model series contains a nonfinite value")
    for (left_time, _), (right_time, _) in zip(ordered, ordered[1:]):
        if right_time <= left_time:
            raise CampaignError("model-series epochs are not strictly increasing")
    if end < start:
        raise CampaignError("observation window ends before it starts")
    if start < ordered[0][0] or end > ordered[-1][0]:
        raise CampaignError("model series does not cover the observation window")
    if start == end:
        return _interpolate(ordered, start)

    knots = [(start, _interpolate(ordered, start))]
    knots.extend((time, value) for time, value in ordered if start < time < end)
    knots.append((end, _interpolate(ordered, end)))
    integral = 0.0
    for (left_time, left_value), (right_time, right_value) in zip(knots, knots[1:]):
        seconds = (right_time - left_time).total_seconds()
        integral += 0.5 * (left_value + right_value) * seconds
    return integral / (end - start).total_seconds()


def _key(record: Mapping[str, object]) -> Tuple[str, str, str, str, str, str]:
    return (
        str(record["instrument_id"]), str(record["platform"]),
        str(record["channel"]), str(record["direction"]), str(record["quantity"]),
        str(record["units"]),
    )


def _pearson(pairs: Sequence[Tuple[float, float]]) -> Optional[float]:
    """Return the ordinary Pearson correlation or ``None`` if undefined."""

    if len(pairs) < 2:
        return None
    observed = [item[0] for item in pairs]
    modeled = [item[1] for item in pairs]
    observed_mean = statistics.fmean(observed)
    modeled_mean = statistics.fmean(modeled)
    observed_ss = sum((value - observed_mean) ** 2 for value in observed)
    modeled_ss = sum((value - modeled_mean) ** 2 for value in modeled)
    if observed_ss == 0.0 or modeled_ss == 0.0:
        return None
    covariance = sum(
        (left - observed_mean) * (right - modeled_mean)
        for left, right in pairs
    )
    return covariance / math.sqrt(observed_ss * modeled_ss)


def compare_with_observations(
    predictions: Sequence[Mapping[str, object]],
    observations: Sequence[Mapping[str, object]],
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    """Exposure-average predictions and compute applicable Step-8 metrics."""

    series: Dict[
        Tuple[str, str, str, str, str, str],
        List[Tuple[datetime, float, Optional[float], Optional[float]]],
    ] = {}
    for row in predictions:
        value = _finite(row["value"], "model prediction")
        lower = (
            _finite(row["lower"], "model lower bound")
            if row.get("lower") is not None else None
        )
        upper = (
            _finite(row["upper"], "model upper bound")
            if row.get("upper") is not None else None
        )
        if lower is not None and lower > value:
            raise CampaignError("model lower bound exceeds its prediction")
        if upper is not None and upper < value:
            raise CampaignError("model upper bound is below its prediction")
        series.setdefault(_key(row), []).append((
            parse_utc(str(row["utc"])), value, lower, upper,
        ))

    comparisons: List[Dict[str, object]] = []
    missing_coverage = 0
    for observation in observations:
        values = sorted(series.get(_key(observation), []), key=lambda item: item[0])
        if not values:
            missing_coverage += 1
            continue
        start = parse_utc(str(observation["start_utc"]))
        end = parse_utc(str(observation["end_utc"]))
        try:
            model = cadence_average([(item[0], item[1]) for item in values], start, end)
            lower = None
            upper = None
            if all(item[2] is not None for item in values):
                lower = cadence_average([(item[0], float(item[2])) for item in values], start, end)
            if all(item[3] is not None for item in values):
                upper = cadence_average([(item[0], float(item[3])) for item in values], start, end)
        except CampaignError:
            missing_coverage += 1
            continue
        observed = float(observation["value"])
        log_ratio = math.log10(model / observed) if model > 0.0 and observed > 0.0 else None
        comparisons.append({
            "instrument_id": observation["instrument_id"],
            "platform": observation["platform"],
            "channel": observation["channel"],
            "direction": observation["direction"],
            "quantity": observation["quantity"],
            "units": observation["units"],
            "start_utc": observation["start_utc"],
            "end_utc": observation["end_utc"],
            "observed": observed,
            "modeled": model,
            "modeled_lower": lower,
            "modeled_upper": upper,
            "signed_error": model - observed,
            "absolute_error": abs(model - observed),
            "log10_ratio": log_ratio,
            "inside_model_interval": (
                lower <= observed <= upper if lower is not None and upper is not None else None
            ),
        })

    logs = [float(row["log10_ratio"]) for row in comparisons if row["log10_ratio"] is not None]
    pairs = [
        (float(row["observed"]), float(row["modeled"]))
        for row in comparisons
    ]
    cutoff = [row for row in comparisons if row["quantity"] == "cutoff_latitude"]
    directional = [
        row for row in comparisons
        if row["quantity"] == "directional_ratio" and row["log10_ratio"] is not None
    ]
    # A ratio above/below unity encodes the modeled east-west asymmetry sign.
    # Exact unity carries no sign and is counted as correct only when both model
    # and observation are numerically unity.
    directional_sign = []
    for row in directional:
        observed_offset = float(row["observed"]) - 1.0
        modeled_offset = float(row["modeled"]) - 1.0
        if observed_offset == 0.0:
            directional_sign.append(modeled_offset == 0.0)
        else:
            directional_sign.append(observed_offset * modeled_offset > 0.0)
    covered = [row for row in comparisons if row["inside_model_interval"] is not None]
    metrics: Dict[str, object] = {
        "observation_count": len(observations),
        "comparison_count": len(comparisons),
        "missing_model_coverage_count": missing_coverage,
        "comparison_fraction": (
            len(comparisons) / len(observations) if observations else None
        ),
        "mean_bias": (
            statistics.fmean(float(row["signed_error"]) for row in comparisons)
            if comparisons else None
        ),
        "mean_absolute_error": (
            statistics.fmean(float(row["absolute_error"]) for row in comparisons)
            if comparisons else None
        ),
        "rmse": (
            math.sqrt(
                statistics.fmean(
                    float(row["signed_error"]) ** 2 for row in comparisons
                )
            ) if comparisons else None
        ),
        "pearson_correlation": _pearson(pairs),
        "mean_log10_bias": statistics.fmean(logs) if logs else None,
        "log10_rmse": (
            math.sqrt(statistics.fmean(value ** 2 for value in logs))
            if logs else None
        ),
        "median_abs_log10_ratio": statistics.median(abs(value) for value in logs) if logs else None,
        "factor_two_fraction": (
            sum(abs(value) <= math.log10(2.0) for value in logs) / len(logs) if logs else None
        ),
        "cutoff_mae_deg": (
            statistics.fmean(float(row["absolute_error"]) for row in cutoff) if cutoff else None
        ),
        "cutoff_bias_deg": (
            statistics.fmean(float(row["signed_error"]) for row in cutoff) if cutoff else None
        ),
        "directional_log10_rmse": (
            math.sqrt(statistics.fmean(float(row["log10_ratio"]) ** 2 for row in directional))
            if directional else None
        ),
        "directional_sign_fraction": (
            sum(directional_sign) / len(directional_sign) if directional_sign else None
        ),
        "interval_coverage_fraction": (
            sum(bool(row["inside_model_interval"]) for row in covered) / len(covered)
            if covered else None
        ),
    }
    return comparisons, metrics


def write_comparison_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    fields = [
        "instrument_id", "platform", "channel", "direction", "quantity", "units",
        "start_utc", "end_utc", "observed", "modeled", "modeled_lower",
        "modeled_upper", "signed_error", "absolute_error", "log10_ratio",
        "inside_model_interval",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

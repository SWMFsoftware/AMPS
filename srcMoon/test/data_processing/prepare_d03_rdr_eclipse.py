#!/usr/bin/env python3
"""Create the frozen I17 Diviner eclipse product from immutable PDS4 RDRs.

This processor never reads AMPS output.  It applies the reviewed selection in
``I17/reference/selection_contract.json`` identically to the June training day
and December hold-out day.  The hold-out observations are processed only after
the selection contract exists on disk and its hash is recorded in the QA
receipt; no model residual can alter membership, binning, flags, or uncertainty.

The RDR supplies both calibrated radiance and its native calibrated brightness
temperature for every detector record.  Because the archive does not publish a
machine-readable inverse calibration function, this script reconstructs that
native monotone relation from the paired RDR fields and verifies its row-level
residual.  It does *not* substitute a rectangular detector bandpass for native
calibration.  The rectangular approximation is retained only by production as
a clearly labelled C6--C9 diagnostic, in accordance with the calibration
report; it is not the observational T_bol conversion or an acceptance oracle.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import dataclasses
import datetime as dt
import functools
import hashlib
import io
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile


REPO = Path(__file__).resolve().parents[3]
CONTRACT_PATH = (REPO / "srcMoon/test/integration-tests/I17/reference/"
                 "selection_contract.json")
DEFAULT_DATA_ROOT = Path("/data/vtenishe/moon_validation_data")
CHANNELS = tuple(range(3, 10))
MISSING = {-9999.0, -9998.0}
SIGMA = 5.670374419e-8
PLANCK = 6.62607015e-34
LIGHT = 299792458.0
BOLTZMANN = 1.380649e-23
PI = math.pi

# Diviner Radiometric Test Report Tables 4 and 10.  Table 4 exitances are
# converted from W cm^-2 to radiance W m^-2 sr^-1 by 1e4/pi.  These are
# channel-dependent calibration quantities, not empirical AMPS tolerances.
SIGNAL_300_COUNTS = {
    3: 1240.0, 4: 640.0, 5: 640.0, 6: 3730.0,
    7: 2850.0, 8: 740.0, 9: 245.0,
}
EXITANCE_300_W_CM2 = {
    3: 1.33e-3, 4: 8.80e-4, 5: 9.02e-4, 6: 1.65e-2,
    7: 5.44e-3, 8: 1.21e-3, 9: 1.87e-4,
}
INFLIGHT_FRACTIONAL_GAIN = {
    3: 1.22e-2, 4: 1.16e-2, 5: 1.13e-2, 6: 6.23e-3,
    7: 4.10e-3, 8: 3.94e-3, 9: 7.57e-3,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_utc(date_text: str, time_text: str) -> dt.datetime:
    # Native RDR time is SCET UTC with the non-ISO HH:MM:SS.sss field.  The
    # date has an English textual month.  Decode that month explicitly rather
    # than relying on the host's LC_TIME setting: a validation product must be
    # identical on an English or non-English machine.
    months = {name: number for number, name in enumerate(
        ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
         "Sep", "Oct", "Nov", "Dec"), 1)}
    date_fields = date_text.strip().split("-")
    if len(date_fields) != 3 or date_fields[1] not in months:
        raise ValueError(f"invalid native RDR date {date_text!r}")
    time_fields = time_text.strip().split(":")
    if len(time_fields) != 3:
        raise ValueError(f"invalid native RDR time {time_text!r}")
    seconds = float(time_fields[2])
    whole_seconds = int(seconds)
    microseconds = round((seconds - whole_seconds) * 1_000_000)
    base = dt.datetime(
        int(date_fields[2]), months[date_fields[1]], int(date_fields[0]),
        int(time_fields[0]), int(time_fields[1]), whole_seconds,
        tzinfo=dt.timezone.utc)
    # timedelta also handles the theoretically possible round-up from
    # 999999.5 microseconds without constructing an invalid datetime field.
    return base + dt.timedelta(microseconds=microseconds)


def planck_radiance_per_m(wavelength_m: float, temperature_k: float) -> float:
    exponent = PLANCK * LIGHT / (wavelength_m * BOLTZMANN * temperature_k)
    if exponent > 700.0:
        return 0.0
    return 2.0 * PLANCK * LIGHT**2 / wavelength_m**5 / math.expm1(exponent)


def planck_dimensionless_tail(x: float) -> float:
    """Return integral_x^infinity t^3/(exp(t)-1) dt reproducibly.

    Transforming wavelength to ``t=hc/(lambda*k*T)`` removes the very broad
    C9 wavelength interval and permits an exponentially convergent exact
    series.  Near zero the direct tail series converges slowly, so the
    Bernoulli expansion integrates 0..x and subtracts it from pi^4/15.  The
    x<0.5 expansion's first omitted term is below double-precision relevance
    for this application.  This is the same Planck integral as wavelength
    Simpson quadrature, not a rectangular/band-centre approximation.
    """
    if x < 0.0 or math.isnan(x):
        raise ValueError("dimensionless Planck limit must be non-negative")
    if math.isinf(x):
        return 0.0
    total = PI**4 / 15.0
    if x < 0.5:
        integral_zero_to_x = (
            x**3 / 3.0 - x**4 / 8.0 + x**5 / 60.0 -
            x**7 / 5040.0 + x**9 / 272160.0 -
            x**11 / 13305600.0 + x**13 / 622702080.0 -
            691.0 * x**15 / (15.0 * 1307674368000.0))
        return total - integral_zero_to_x

    result = 0.0
    for integer in range(1, 100000):
        inverse = 1.0 / integer
        exponential = math.exp(-integer * x)
        term = exponential * (
            x**3 * inverse + 3.0 * x**2 * inverse**2 +
            6.0 * x * inverse**3 + 6.0 * inverse**4)
        result += term
        if term <= 1.0e-15 * max(result, 1.0):
            return result
    raise RuntimeError("dimensionless Planck tail series did not converge")


@functools.lru_cache(maxsize=200000)
def integrate_planck_um(low_um: float, high_um: float,
                        temperature_k: float) -> float:
    # Integrated spectral radiance is sigma*T^4/pi times the blackbody energy
    # fraction between the two dimensionless limits.  Cache exact repeated
    # calls made by uncertainty perturbations; this changes neither values nor
    # reduction order and keeps memory bounded explicitly.
    if not (0.0 <= low_um < high_um and temperature_k > 0.0 and
            math.isfinite(temperature_k)):
        raise ValueError("invalid Planck integration interval/temperature")
    second_radiation_constant_m_k = PLANCK * LIGHT / BOLTZMANN
    lower_x = second_radiation_constant_m_k / (high_um * 1.0e-6 * temperature_k)
    upper_x = (math.inf if low_um == 0.0 else
               second_radiation_constant_m_k /
               (low_um * 1.0e-6 * temperature_k))
    fraction = 15.0 / PI**4 * (
        planck_dimensionless_tail(lower_x) -
        planck_dimensionless_tail(upper_x))
    if not (0.0 <= fraction <= 1.0 + 1.0e-14):
        raise RuntimeError("dimensionless Planck band fraction is non-physical")
    return fraction * SIGMA * temperature_k**4 / PI


PARTITIONS = ((0.0, 8.075), (8.075, 8.40), (8.40, 13.0),
              (13.0, 25.0), (25.0, 50.0), (50.0, 100.0),
              (100.0, 1000.0))
TMIN = (190.0, 190.0, 180.0, 95.0, 60.0, 40.0, None)


def tbol(channel_temperature_k: list[float]) -> tuple[float, list[float]]:
    if len(channel_temperature_k) != 7 or not math.isfinite(
            channel_temperature_k[-1]) or channel_temperature_k[-1] <= 0.0:
        raise ValueError("C9 is required for the Paige T_bol operator")
    used = list(channel_temperature_k)
    for index in range(5, -1, -1):
        threshold = TMIN[index]
        if (not math.isfinite(channel_temperature_k[index]) or
                channel_temperature_k[index] < float(threshold)):
            used[index] = used[index + 1]
    fourth_power = 0.0
    for temperature, limits in zip(used, PARTITIONS):
        fraction = integrate_planck_um(*limits, temperature) / \
            (SIGMA * temperature**4 / PI)
        fourth_power += temperature**4 * fraction
    return fourth_power**0.25, used


def tbol_uncertainty(channel_temperature_k: list[float],
                     channel_uncertainty_k: list[float]) -> float:
    """Propagate independent channel uncertainties through Paige T_bol.

    This is the first-order (Jacobian/delta-method) propagation appropriate to
    the one-sigma channel uncertainties produced above.  The derivative step
    is deliberately local and remains on the same side of every published
    sensitivity threshold.  Perturbing by a complete one-sigma interval would
    instead change the discrete channel-substitution branch and is not a
    derivative; for cold C9 it can also request a non-physical temperature at
    or below zero.  A below-threshold C3--C8 channel has exactly zero local
    derivative because the published operator substitutes the next usable
    longer-wavelength channel.

    Diviner does not publish the cross-channel covariance needed for a
    correlated propagation.  We consequently add the independent derivative
    contributions in quadrature and record that limitation in the QA receipt
    and I17 contract rather than inventing correlations.
    """
    if len(channel_temperature_k) != 7 or len(channel_uncertainty_k) != 7:
        raise ValueError("T_bol uncertainty requires C3-C9 values")
    if any(not math.isfinite(value) or value < 0.0
           for value in channel_uncertainty_k):
        raise ValueError("channel uncertainties must be finite and non-negative")

    baseline, _ = tbol(channel_temperature_k)
    variance = 0.0
    for index, (temperature, uncertainty) in enumerate(zip(
            channel_temperature_k, channel_uncertainty_k)):
        if uncertainty == 0.0:
            continue
        # Non-finite or sub-threshold C3--C8 inputs are absent from the Paige
        # sum after its prescribed substitution and therefore have no local
        # uncertainty contribution.  C9 validity was already enforced by
        # tbol() above.
        lower_branch_limit = (0.0 if index == 6 else float(TMIN[index]))
        if not math.isfinite(temperature) or temperature < lower_branch_limit:
            continue

        # A relative 1e-5 step is comfortably above floating-point roundoff
        # while remaining local compared with RDR calibration uncertainties.
        # Restrict a centred step to one quarter of the distance from the
        # positivity/sensitivity boundary so it cannot switch operator branch.
        step = max(abs(temperature) * 1.0e-5, 1.0e-5)
        distance = temperature - lower_branch_limit
        if distance > 0.0:
            step = min(step, 0.25 * distance)
            plus = list(channel_temperature_k)
            minus = list(channel_temperature_k)
            plus[index] += step
            minus[index] -= step
            derivative = (tbol(plus)[0] - tbol(minus)[0]) / (2.0 * step)
        else:
            # Exactly at a threshold the operator is non-differentiable.  Use
            # the valid-branch (right-hand) derivative and preserve this
            # explicit convention in the machine-readable contract.
            plus = list(channel_temperature_k)
            plus[index] += step
            derivative = (tbol(plus)[0] - baseline) / step
        variance += (derivative * uncertainty)**2

    return math.sqrt(variance)


@dataclasses.dataclass
class Sum:
    """Streaming statistics for one channel within one site/time bin."""

    n: int = 0
    radiance_sum: float = 0.0
    radiance_sum2: float = 0.0
    latitude_sum: float = 0.0
    longitude_x_sum: float = 0.0
    longitude_y_sum: float = 0.0
    epoch_sum: float = 0.0
    emission_sum: float = 0.0
    detector_mask: int = 0
    eclipse_records: int = 0
    qge_value_mask: int = 0
    qmi_value_mask: int = 0

    def add(self, radiance: float, latitude: float, longitude: float,
            epoch_s: float, emission: float, detector: int,
            qge: int, qmi: int) -> None:
        self.n += 1
        self.radiance_sum += radiance
        self.radiance_sum2 += radiance * radiance
        self.latitude_sum += latitude
        angle = math.radians(longitude)
        self.longitude_x_sum += math.cos(angle)
        self.longitude_y_sum += math.sin(angle)
        self.epoch_sum += epoch_s
        self.emission_sum += emission
        self.detector_mask |= 1 << (detector - 1)
        self.eclipse_records += int(bool(qmi & 2))
        # Values, not merely record counts, are retained.  A bit at index v
        # means that native decimal flag value v occurred in this aggregate.
        # The accepted masks restrict qge to 0/4/8/12 and qmi to 0/2, so these
        # small integers do not create a material per-bin memory cost.
        self.qge_value_mask |= 1 << qge
        self.qmi_value_mask |= 1 << qmi

    @property
    def mean_radiance(self) -> float:
        return self.radiance_sum / self.n

    @property
    def radiance_standard_error(self) -> float:
        if self.n < 2:
            return math.nan
        variance = max(0.0, (self.radiance_sum2 -
                            self.radiance_sum**2 / self.n) / (self.n - 1))
        return math.sqrt(variance / self.n)


class NativeCalibration:
    """Monotone inverse of the paired native RDR radiance/Tb calibration."""

    def __init__(self, pairs: dict[float, list[float]]) -> None:
        # RDR radiance is rounded more coarsely than brightness temperature.
        # Consequently many native Tb values can share one archived radiance;
        # grouping in the *radiance* coordinate and averaging those paired Tb
        # values is the correct inverse relation.  Grouping by Tb and deleting
        # duplicate radiances biases the long-wave channels. ``values`` is
        # [temperature_sum, record_count], so memory is proportional to the
        # number of archived radiance levels rather than tens of millions of
        # records.
        points = sorted((radiance, values[0] / values[1])
                        for radiance, values in pairs.items())
        self.radiance: list[float] = []
        self.temperature: list[float] = []
        for radiance, temperature in points:
            if not self.temperature or temperature > self.temperature[-1]:
                self.radiance.append(radiance)
                self.temperature.append(temperature)
        if len(self.radiance) < 100:
            raise ValueError("native RDR calibration relation has too few points")

    def brightness_temperature(self, radiance: float) -> float:
        index = bisect.bisect_left(self.radiance, radiance)
        if index == 0:
            return self.temperature[0]
        if index == len(self.radiance):
            return self.temperature[-1]
        x0, x1 = self.radiance[index - 1:index + 1]
        y0, y1 = self.temperature[index - 1:index + 1]
        return y0 + (radiance - x0) * (y1 - y0) / (x1 - x0)

    def reconstruction_statistics(
            self, pairs: dict[float, list[float]], minimum_temperature_k: float
    ) -> dict[str, float | int]:
        """Report archive-pair reconstruction above the usable sensitivity.

        Mineral channels contain physically uninformative, sometimes
        negative Tb values at low signal.  Paige's published sensitivity rule
        replaces them, so judging the monotone inverse on those values would
        reject the archive for precisely the samples the publication says not
        to use.  C9 has no sensitivity floor and is checked for positive Tb.
        """
        residuals = []
        for radiance, values in pairs.items():
            native_temperature = values[0] / values[1]
            if native_temperature >= minimum_temperature_k:
                residuals.append(abs(
                    self.brightness_temperature(radiance) - native_temperature))
        residuals.sort()
        if not residuals:
            raise ValueError(
                "native calibration has no nodes above its sensitivity floor")
        return {
            "usable_nodes": len(residuals),
            "maximum_absolute_residual_K": residuals[-1],
            "p99_absolute_residual_K": residuals[
                int(0.99 * (len(residuals) - 1))],
        }


def mask_values(mask: int, value_offset: int = 0) -> str:
    """Serialize the distinct native integer flag values without CSV commas."""
    return "|".join(str(value + value_offset)
                    for value in range(mask.bit_length())
                    if mask & (1 << value))


def population_count(mask: int) -> int:
    """Count set bits on the repository's supported Python 3.8 runtime.

    ``int.bit_count`` was added only in Python 3.8's successors.  The masks
    here contain at most 21 detector bits, so the portable binary-string count
    is exact and its negligible cost is paid once per completed channel bin.
    """
    if mask < 0:
        raise ValueError("population_count requires a non-negative mask")
    return bin(mask).count("1")


def record_fields(row: list[str]) -> tuple[object, ...]:
    """Parse only fields whose units/meaning are fixed by the native SIS."""
    if len(row) != 33:
        raise ValueError(f"RDR record has {len(row)} fields rather than 33")
    return (parse_utc(row[0].strip(' "'), row[1].strip(' "')),
            int(row[14]), int(row[17]), int(row[18]), float(row[22]),
            float(row[23]), float(row[24]), float(row[25]), float(row[26]),
            int(row[30]), int(row[31]), int(row[32]))


def validate_label(label: Path, archive_name: str) -> str:
    root = ET.parse(label).getroot()
    ns = {"p": "http://pds.nasa.gov/pds4/pds/v1"}
    lid = root.findtext("p:Identification_Area/p:logical_identifier", namespaces=ns)
    file_name = root.findtext(".//p:File/p:file_name", namespaces=ns)
    if not lid or ":lro_diviner_calibrated:data_calibrated1:" not in lid:
        raise ValueError(f"{label} is not a calibrated1 PDS4 product")
    if not file_name or Path(file_name).stem.lower() != Path(archive_name).stem.lower():
        raise ValueError(f"{label} does not identify {archive_name}")
    return lid


def calibration_sigma_radiance(channel: int, radiance: float) -> float:
    l300 = EXITANCE_300_W_CM2[channel] * 1.0e4 / PI
    offset = l300 * 3.6 / SIGNAL_300_COUNTS[channel]
    # The report states that in-flight uncertainty is also an upper limit for
    # ground calibration.  They are independent contributions in its roll-up.
    gain = math.sqrt(2.0) * INFLIGHT_FRACTIONAL_GAIN[channel] * abs(radiance)
    return math.hypot(offset, gain)


def record_is_selected(
        epoch: dt.datetime, begin: dt.datetime, end: dt.datetime,
        activity: int, latitude: float, longitude: float, emission: float,
        qca: int, qge: int, qmi: int, quality: dict[str, object]) -> bool:
    """Apply only the predeclared native-record selection.

    Keeping the predicate separate makes its mask semantics directly
    testable.  ``allowed_mask`` means that only those bits may be present; it
    does not require every allowed bit to be set.  Thus qge 0/4/8/12 and qmi
    0/2 are retained, while any predictive/missing geometry or instrument
    anomaly bit rejects the record.
    """
    allowed_activity = {int(value) for value in
                        quality["activity_flag_allowed"]}  # type: ignore[index]
    maximum_emission = float(
        quality["maximum_emission_angle_degrees_exclusive"])
    qca_mask = int(quality["qca_allowed_mask"])
    qge_mask = int(quality["qge_allowed_mask"])
    qmi_mask = int(quality["qmi_allowed_mask"])
    return (
        begin <= epoch <= end and activity in allowed_activity and
        math.isfinite(latitude) and latitude not in MISSING and
        -90.0 <= latitude <= 90.0 and
        math.isfinite(longitude) and longitude not in MISSING and
        math.isfinite(emission) and 0.0 <= emission < maximum_emission and
        (qca & ~qca_mask) == 0 and (qge & ~qge_mask) == 0 and
        (qmi & ~qmi_mask) == 0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=CONTRACT_PATH)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    bins: dict[tuple[str, int, int, int], dict[int, Sum]] = {}
    calibration_pairs = {channel: {} for channel in CHANNELS}
    counters = {"records_total": 0, "records_selected": 0,
                "records_rejected": 0, "archives": 0}
    selected_by_role = {"training": 0, "holdout": 0}
    raw_files: list[dict[str, object]] = []
    bin_seconds = int(contract["binning"]["time_bin_seconds"])
    spatial = float(contract["binning"]["latitude_bin_degrees"])
    longitude_spatial = float(contract["binning"]["longitude_bin_degrees"])
    quality = contract["quality"]

    for role in ("training", "holdout"):
        selection = contract["selections"][role]
        begin = dt.datetime.fromisoformat(selection["begin_utc"].replace("Z", "+00:00"))
        end = dt.datetime.fromisoformat(selection["end_utc"].replace("Z", "+00:00"))
        day = begin.date()
        directory = (args.data_root / "diviner/raw/rdr" /
                     f"{day.year:04d}" / f"{day.year:04d}{day.month:02d}" /
                     f"{day.year:04d}{day.month:02d}{day.day:02d}")
        archives = sorted(directory.glob("*_rdr.zip"))
        if len(archives) != 144:
            raise FileNotFoundError(f"{directory} has {len(archives)} RDR archives; expected 144")
        for archive in archives:
            # Filename minute provides a cheap prefilter only; each retained
            # record is still checked against its native UTC because products
            # can cross nominal file boundaries by one acquisition frame.
            stamp = dt.datetime.strptime(archive.name[:12], "%Y%m%d%H%M").replace(
                tzinfo=dt.timezone.utc)
            if stamp + dt.timedelta(minutes=11) < begin or stamp > end:
                continue
            label = archive.with_suffix(".xml")
            product_lid = validate_label(label, archive.name)
            raw_files.extend([
                {"path": str(archive), "bytes": archive.stat().st_size,
                 "sha256": sha256(archive), "product_lid": product_lid},
                {"path": str(label), "bytes": label.stat().st_size,
                 "sha256": sha256(label), "product_lid": product_lid},
            ])
            counters["archives"] += 1
            print(f"processing {role} archive {archive.name}", flush=True)
            with zipfile.ZipFile(archive) as compressed:
                if compressed.testzip() is not None or len(compressed.namelist()) != 1:
                    raise ValueError(f"ZIP CRC/member validation failed for {archive}")
                with compressed.open(compressed.namelist()[0]) as binary:
                    stream = io.TextIOWrapper(binary, encoding="ascii", newline="")
                    for _ in range(4):
                        next(stream)
                    for row in csv.reader(stream, skipinitialspace=True):
                        counters["records_total"] += 1
                        (epoch, activity, channel, detector, radiance, tb,
                         latitude, longitude, emission, qca, qge, qmi) = record_fields(row)
                        if (channel not in CHANNELS or radiance in MISSING or
                                tb in MISSING or not math.isfinite(radiance) or
                                not math.isfinite(tb)):
                            counters["records_rejected"] += 1
                            continue
                        accepted = record_is_selected(
                            epoch, begin, end, activity, latitude, longitude,
                            emission, qca, qge, qmi, quality)
                        if not accepted:
                            counters["records_rejected"] += 1
                            continue
                        # Construct the native inverse only from records that
                        # meet the exact science selection.  This prevents a
                        # poor-calibration record rejected by qca from
                        # influencing a retained bin through the conversion
                        # table.  The paired Tb is archive calibration—not an
                        # AMPS-derived or hand-tuned value.
                        accumulator = calibration_pairs[channel].setdefault(
                            radiance, [0.0, 0.0])
                        accumulator[0] += tb
                        accumulator[1] += 1.0
                        seconds = (epoch - dt.datetime.combine(
                            day, dt.time(), tzinfo=dt.timezone.utc)).total_seconds()
                        key = (role, int(seconds // bin_seconds),
                               math.floor((latitude + 90.0) / spatial),
                               math.floor((longitude % 360.0) /
                                          longitude_spatial))
                        channel_sum = bins.setdefault(key, {}).setdefault(channel, Sum())
                        channel_sum.add(radiance, latitude, longitude % 360.0,
                                        epoch.timestamp(), emission, detector,
                                        qge, qmi)
                        counters["records_selected"] += 1
                        selected_by_role[role] += 1

    calibration = {channel: NativeCalibration(calibration_pairs[channel])
                   for channel in CHANNELS}
    output_path = args.output_dir / "eclipse_thermal_history.csv"
    temporary_output_path = output_path.with_suffix(output_path.suffix + ".tmp")
    fieldnames = [
        "selection", "site_time_bin_id", "utc_mean", "time_bin_begin_utc",
        "latitude_bin_center_deg", "longitude_east_bin_center_deg",
        "bolometric_brightness_temperature_K",
        "bolometric_brightness_temperature_uncertainty_K",
        "eclipse_record_fraction",
    ]
    for channel in CHANNELS:
        fieldnames.extend([
            f"c{channel}_count", f"c{channel}_detector_count",
            f"c{channel}_detector_ids", f"c{channel}_qca_values",
            f"c{channel}_qge_values", f"c{channel}_qmi_values",
            f"c{channel}_latitude_mean_deg", f"c{channel}_longitude_east_mean_deg",
            f"c{channel}_utc_mean", f"c{channel}_emission_angle_mean_deg",
            f"c{channel}_radiance_W_m-2_sr-1",
            f"c{channel}_radiance_uncertainty_W_m-2_sr-1",
            f"c{channel}_brightness_temperature_K",
            f"c{channel}_brightness_temperature_uncertainty_K",
            f"c{channel}_temperature_used_for_tbol_K",
        ])

    produced = 0
    produced_by_role = {"training": 0, "holdout": 0}
    incomplete = 0
    nonpositive_c9 = 0
    # Write the potentially large science table under an unmistakable
    # temporary name.  The final CSV is published by an atomic rename only
    # after every bin and every QA statistic has been computed successfully.
    # A killed processor may leave ``*.tmp`` work, but never a new file whose
    # final name could be mistaken for qualified evidence.
    with temporary_output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for key in sorted(bins):
            role, time_index, latitude_index, longitude_index = key
            channel_sums = bins[key]
            minimum = int(contract["binning"]["minimum_records_per_channel"])
            if any(channel not in channel_sums or channel_sums[channel].n < minimum
                   for channel in CHANNELS):
                incomplete += 1
                continue
            temperatures: list[float] = []
            uncertainties: list[float] = []
            row: dict[str, object] = {"selection": role}
            epoch_weighted_sum = 0.0
            eclipse_count = 0
            total_count = 0
            for channel in CHANNELS:
                values = channel_sums[channel]
                mean_radiance = values.mean_radiance
                within = values.radiance_standard_error
                sigma_radiance = math.hypot(
                    calibration_sigma_radiance(channel, mean_radiance), within)
                temperature = calibration[channel].brightness_temperature(mean_radiance)
                lower = calibration[channel].brightness_temperature(
                    mean_radiance - sigma_radiance)
                upper = calibration[channel].brightness_temperature(
                    mean_radiance + sigma_radiance)
                sigma_temperature = 0.5 * abs(upper - lower)
                temperatures.append(temperature)
                uncertainties.append(sigma_temperature)
                mean_epoch = values.epoch_sum / values.n
                mean_longitude = math.degrees(math.atan2(
                    values.longitude_y_sum, values.longitude_x_sum)) % 360.0
                epoch_weighted_sum += values.epoch_sum
                eclipse_count += values.eclipse_records
                total_count += values.n
                prefix = f"c{channel}_"
                row.update({
                    prefix + "count": values.n,
                    prefix + "detector_count": population_count(
                        values.detector_mask),
                    prefix + "detector_ids": mask_values(
                        values.detector_mask, value_offset=1),
                    # qca is exactly zero by the frozen selection; preserving
                    # it explicitly keeps all three native flag families in
                    # every processed aggregate.
                    prefix + "qca_values": "0",
                    prefix + "qge_values": mask_values(values.qge_value_mask),
                    prefix + "qmi_values": mask_values(values.qmi_value_mask),
                    prefix + "latitude_mean_deg": values.latitude_sum / values.n,
                    prefix + "longitude_east_mean_deg": mean_longitude,
                    prefix + "utc_mean": dt.datetime.fromtimestamp(
                        mean_epoch, tz=dt.timezone.utc).isoformat().replace("+00:00", "Z"),
                    prefix + "emission_angle_mean_deg": values.emission_sum / values.n,
                    prefix + "radiance_W_m-2_sr-1": mean_radiance,
                    prefix + "radiance_uncertainty_W_m-2_sr-1": sigma_radiance,
                    prefix + "brightness_temperature_K": temperature,
                    prefix + "brightness_temperature_uncertainty_K": sigma_temperature,
                })
            # C9 is the required long-wave anchor.  A non-positive binned C9
            # value is a real low-signal gap, not corrupt infrastructure and
            # not permission to clamp or extrapolate it.  Preserve the count
            # in QA and omit the bin from observational evidence.
            if not math.isfinite(temperatures[-1]) or temperatures[-1] <= 0.0:
                nonpositive_c9 += 1
                incomplete += 1
                continue
            bolometric, used = tbol(temperatures)
            # Propagate the native calibration/dispersion uncertainty through
            # the exact nonlinear operator with its documented local Jacobian.
            # The helper keeps finite-difference samples on the same published
            # sensitivity branch and remains valid when cold-channel one-sigma
            # intervals extend below zero.
            bolometric_uncertainty = tbol_uncertainty(
                temperatures, uncertainties)
            for index in range(len(CHANNELS)):
                row[f"c{index + 3}_temperature_used_for_tbol_K"] = used[index]
            date = contract["selections"][role]["date_utc"]
            bin_begin = dt.datetime.fromisoformat(date).replace(
                tzinfo=dt.timezone.utc) + dt.timedelta(seconds=time_index * bin_seconds)
            row.update({
                "site_time_bin_id": f"{role}-t{time_index:03d}-lat{latitude_index:03d}-lon{longitude_index:03d}",
                "utc_mean": dt.datetime.fromtimestamp(
                    epoch_weighted_sum / total_count,
                    tz=dt.timezone.utc).isoformat().replace("+00:00", "Z"),
                "time_bin_begin_utc": bin_begin.isoformat().replace("+00:00", "Z"),
                "latitude_bin_center_deg": -90.0 + (latitude_index + 0.5) * spatial,
                "longitude_east_bin_center_deg":
                    (longitude_index + 0.5) * longitude_spatial,
                "bolometric_brightness_temperature_K": bolometric,
                "bolometric_brightness_temperature_uncertainty_K": bolometric_uncertainty,
                "eclipse_record_fraction": eclipse_count / total_count,
            })
            writer.writerow(row)
            produced += 1
            produced_by_role[role] += 1

    calibration_qa = {}
    for index, channel in enumerate(CHANNELS):
        # For C9 the published operator requires positive temperature.  For
        # C3--C8 use the published sensitivity floors; colder values are
        # deliberately replaced by a longer-wave channel.
        floor = TMIN[index] if TMIN[index] is not None else 0.0
        calibration_qa[str(channel)] = {
            "monotone_mapping_nodes": len(calibration[channel].radiance),
            **calibration[channel].reconstruction_statistics(
                calibration_pairs[channel], float(floor)),
        }

    contract_snapshot = args.output_dir / "selection_contract.json"
    temporary_contract = contract_snapshot.with_suffix(".json.tmp")
    temporary_contract.write_bytes(args.contract.read_bytes())
    temporary_contract.replace(contract_snapshot)

    qa = {
        "schema": "srcMoon-D03-RDR-eclipse-QA-v1",
        "status": ("READY" if all(produced_by_role[role] > 0
                                    for role in produced_by_role)
                   else "INCOMPLETE"),
        "processor": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "contract": {
            "source_path": str(args.contract.resolve()),
            "snapshot_path": str(contract_snapshot.resolve()),
            "sha256": sha256(contract_snapshot),
        },
        "raw_files": raw_files,
        "counts": {**counters,
                   "records_selected_by_role": selected_by_role,
                   "complete_bins": produced,
                   "complete_bins_by_role": produced_by_role,
                   "incomplete_bins_rejected": incomplete,
                   "nonpositive_c9_bins_rejected": nonpositive_c9},
        "processed": {"path": str(output_path.resolve()),
                      "bytes": temporary_output_path.stat().st_size,
                      "sha256": sha256(temporary_output_path)},
        "calibration": {
            "source": "paired native calibrated radiance and brightness-temperature RDR fields",
            "native_mapping_QA": calibration_qa,
            "uncertainty_source": "Diviner Radiometric Test Report equations 1-3, Tables 4 and 10",
            "unavailable_covariance": "The archive supplies no cross-channel calibration covariance; independent propagation is reported explicitly.",
        },
        "holdout_policy": "The processor does not read model output; training and holdout use the same pre-hashed selection contract.",
    }
    qa_path = args.output_dir / "d03_rdr_qa.json"
    temporary_qa_path = qa_path.with_suffix(qa_path.suffix + ".tmp")
    temporary_qa_path.write_text(
        json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # Publish the data first and its READY/INCOMPLETE receipt last.  If the
    # second rename is interrupted, an old QA hash cannot match the new CSV,
    # so the linked runner fails closed instead of accepting mixed versions.
    temporary_output_path.replace(output_path)
    temporary_qa_path.replace(qa_path)
    print(json.dumps(qa["counts"], sort_keys=True))
    return 0 if qa["status"] == "READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())

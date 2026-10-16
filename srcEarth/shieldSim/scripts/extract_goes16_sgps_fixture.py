#!/usr/bin/env python3
"""Extract auditable GOES-16 SGPS five-minute records from NOAA netCDF.

This is a data-preparation utility, not part of the shieldSim runtime.  It
requires netCDF4 only when the committed CSV fixture is regenerated.  The
result intentionally retains both directional sensors, channel bounds,
uncertainties, sample counts, quality flags, and the yaw-flip flag. It can
extract either one record or a half-open event interval. It does not infer an
omnidirectional or isotropic spectrum.

The NOAA coordinate is explicitly the *start* of each averaging interval. The
CSV retains that coordinate and adds the derived midpoint and end. This is
important: calling the raw coordinate a midpoint would shift every fluence
sample by 150 seconds and would mishandle event boundaries.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import math
import pathlib
import sys


CHANNEL_NAMES = (
    "P1",
    "P2A",
    "P2B",
    "P3",
    "P4",
    "P5",
    "P6",
    "P7",
    "P8A",
    "P8B",
    "P8C",
    "P9",
    "P10",
)

SOURCE_URL = (
    "https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/"
    "satellite-systems/goesr/solar_proton_events/sgps_sep2017_event_data/"
    "se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc"
)
SOURCE_DOI = "https://doi.org/10.1029/2021SW002750"


def parse_utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise argparse.ArgumentTypeError("timestamp must include UTC offset or Z")
    return parsed.astimezone(dt.timezone.utc)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=pathlib.Path, help="NOAA SGPS L2 netCDF file")
    parser.add_argument("output", type=pathlib.Path, help="output CSV path")
    parser.add_argument(
        "--time",
        type=parse_utc,
        default=None,
        help="extract one exact five-minute record start time",
    )
    parser.add_argument("--start", type=parse_utc, help="inclusive event-record start")
    parser.add_argument("--end", type=parse_utc, help="exclusive event-record end")
    args = parser.parse_args()

    if args.time is not None and (args.start is not None or args.end is not None):
        parser.error("--time cannot be combined with --start/--end")
    if (args.start is None) != (args.end is None):
        parser.error("--start and --end must be supplied together")
    if args.time is None and args.start is None:
        args.time = parse_utc("2017-09-11T07:30:00Z")
    if args.start is not None and args.start >= args.end:
        parser.error("--start must be earlier than --end")

    try:
        from netCDF4 import Dataset, num2date
        from numpy import ma
    except ImportError as exc:
        parser.error(f"regeneration requires the netCDF4 Python package: {exc}")

    required = (
        "L2_SciData_TimeStamp",
        "L1bRecordsInAvg",
        "YawFlipFlag",
        "AvgDiffProtonFlux",
        "AvgDiffProtonFluxUncert",
        "DiffValidL1bSamplesInAvg",
        "DiffDQFdtcSum",
        "DiffDQFoobSum",
        "DiffDQFerrSum",
        "DiffProtonLowerEnergy",
        "DiffProtonUpperEnergy",
    )

    with Dataset(args.input) as dataset:
        missing = [name for name in required if name not in dataset.variables]
        if missing:
            parser.error("input lacks required variables: " + ", ".join(missing))

        timestamp_var = dataset.variables["L2_SciData_TimeStamp"]
        expected_timestamp_description = (
            "Time stamp at the start of the averaging period, in seconds since the J2000 epoch"
        )
        if getattr(timestamp_var, "long_name", None) != expected_timestamp_description:
            parser.error(
                "L2_SciData_TimeStamp no longer declares the averaging-interval start; "
                "do not infer timestamp semantics"
            )
        timestamps = num2date(
            timestamp_var[:],
            timestamp_var.units,
            only_use_cftime_datetimes=False,
        )
        if args.time is not None:
            target_naive = args.time.replace(tzinfo=None)
            indices = [index for index, value in enumerate(timestamps) if value == target_naive]
            if len(indices) != 1:
                parser.error(
                    f"expected one record at {args.time.isoformat()}, found {len(indices)}"
                )
        else:
            start_naive = args.start.replace(tzinfo=None)
            end_naive = args.end.replace(tzinfo=None)
            indices = [
                index for index, value in enumerate(timestamps)
                if start_naive <= value < end_naive
            ]
            expected = int((args.end - args.start).total_seconds() / 300)
            if len(indices) != expected:
                parser.error(
                    f"event interval requires {expected} contiguous five-minute records, "
                    f"found {len(indices)}"
                )
        for previous, current in zip(indices, indices[1:]):
            if (timestamps[current] - timestamps[previous]).total_seconds() != 300:
                parser.error("selected SGPS source records are not contiguous five-minute intervals")

        flux = dataset.variables["AvgDiffProtonFlux"]
        uncertainty = dataset.variables["AvgDiffProtonFluxUncert"]
        valid = dataset.variables["DiffValidL1bSamplesInAvg"]
        dead_time = dataset.variables["DiffDQFdtcSum"]
        out_of_band = dataset.variables["DiffDQFoobSum"]
        dynamic_error = dataset.variables["DiffDQFerrSum"]
        energy_low = dataset.variables["DiffProtonLowerEnergy"]
        energy_high = dataset.variables["DiffProtonUpperEnergy"]

        if flux.shape[1:] != (2, 13):
            parser.error(f"unexpected AvgDiffProtonFlux shape {flux.shape}")
        for variable in (uncertainty, valid, dead_time, out_of_band, dynamic_error):
            if variable.shape != flux.shape:
                parser.error(
                    f"{variable.name} shape {variable.shape} differs from flux shape {flux.shape}"
                )
        if energy_low.shape != (2, 13) or energy_high.shape != (2, 13):
            parser.error("SGPS energy-bound arrays must have shape (2,13)")
        if dataset.variables["L1bRecordsInAvg"].shape != timestamp_var.shape or \
                dataset.variables["YawFlipFlag"].shape != timestamp_var.shape:
            parser.error("SGPS record-count/yaw arrays do not match the time coordinate")
        for variable in (flux, uncertainty):
            if getattr(variable, "units", None) != "protons/(cm^2 sr keV s)":
                parser.error(f"{variable.name} has unexpected units; refusing implicit conversion")
        for variable in (energy_low, energy_high):
            if getattr(variable, "units", None) != "keV":
                parser.error(f"{variable.name} has unexpected units; refusing implicit conversion")

        # Read each selected netCDF variable in one operation. Scalar HDF5
        # indexing inside the CSV loop is orders of magnitude slower and can
        # leave a truncated fixture if a job is interrupted mid-extraction.
        flux_values = flux[indices, :, :]
        uncertainty_values = uncertainty[indices, :, :]
        valid_values = valid[indices, :, :]
        dead_time_values = dead_time[indices, :, :]
        out_of_band_values = out_of_band[indices, :, :]
        dynamic_error_values = dynamic_error[indices, :, :]
        yaw_flip_values = dataset.variables["YawFlipFlag"][indices]
        records_values = dataset.variables["L1bRecordsInAvg"][indices]
        energy_low_values = energy_low[:, :]
        energy_high_values = energy_high[:, :]

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as output:
            output.write("# GOES-16 SGPS L2 five-minute differential proton spectrum\n")
            output.write(f"# source_url: {SOURCE_URL}\n")
            output.write(f"# instrument_paper: {SOURCE_DOI}\n")
            output.write(f"# source_file_sha256: {sha256(args.input)}\n")
            if args.time is not None:
                output.write(
                    f"# extraction_start_utc: "
                    f"{args.time.isoformat().replace('+00:00', 'Z')}\n"
                )
                output.write(
                    f"# extraction_end_utc: "
                    f"{(args.time + dt.timedelta(seconds=300)).isoformat().replace('+00:00', 'Z')}\n"
                )
            else:
                output.write(
                    f"# extraction_start_utc: "
                    f"{args.start.isoformat().replace('+00:00', 'Z')}\n"
                )
                output.write(
                    f"# extraction_end_utc: "
                    f"{args.end.isoformat().replace('+00:00', 'Z')}\n"
                )
            output.write("# source_timestamp_semantics: averaging_interval_start\n")
            output.write("# canonical_timestamp_semantics: derived_midpoint=start+150s\n")
            output.write("# sample_support: [start,start+300s)\n")
            output.write(
                "# sensor semantics: unit 0 is -X and unit 1 is +X; the NOAA "
                "file describes those units as westward and eastward only when "
                "the spacecraft is upright. No angular averaging was applied.\n"
            )
            output.write(
                "# conversion: energies keV->MeV (/1000); differential flux "
                "per keV->per MeV (x1000).\n"
            )
            output.write(
                "# scientific-use guard: this directional geosynchronous "
                "snapshot is a provenance/ingestion fixture, not an isotropic "
                "SEP mission environment.\n"
            )

            writer = csv.writer(output, lineterminator="\n")
            writer.writerow(
                (
                    "sample_start_utc",
                    "sample_midpoint_utc",
                    "sample_end_utc",
                    "l1b_records_in_average",
                    "sensor_unit",
                    "spacecraft_axis",
                    "upright_nominal_look",
                    "yaw_flip_flag",
                    "channel",
                    "energy_low_MeV",
                    "energy_high_MeV",
                    "proton_flux_per_cm2_sr_s_MeV",
                    "flux_uncertainty_per_cm2_sr_s_MeV",
                    "valid_1s_samples",
                    "dead_time_flagged_samples",
                    "out_of_band_flagged_samples",
                    "dynamic_error_flagged_samples",
                )
            )

            def timestamp_text(value: dt.datetime) -> str:
                return value.replace(tzinfo=dt.timezone.utc).isoformat().replace("+00:00", "Z")

            def converted_or_na(value, multiplier: float) -> str:
                if value is None or ma.is_masked(value):
                    return "NA"
                converted = float(value) * multiplier
                return f"{converted:.9g}" if math.isfinite(converted) else "NA"

            for selected_index, index in enumerate(indices):
                start = timestamps[index]
                midpoint = start + dt.timedelta(seconds=150)
                end = start + dt.timedelta(seconds=300)
                yaw_flip = int(yaw_flip_values[selected_index])
                records = int(records_values[selected_index])
                for sensor in range(2):
                    axis = "-X" if sensor == 0 else "+X"
                    nominal_look = "westward" if sensor == 0 else "eastward"
                    for channel, channel_name in enumerate(CHANNEL_NAMES):
                        writer.writerow(
                            (
                                timestamp_text(start),
                                timestamp_text(midpoint),
                                timestamp_text(end),
                                records,
                                sensor,
                                axis,
                                nominal_look,
                                yaw_flip,
                                channel_name,
                                f"{float(energy_low_values[sensor, channel]) / 1000.0:.9g}",
                                f"{float(energy_high_values[sensor, channel]) / 1000.0:.9g}",
                                converted_or_na(flux_values[selected_index, sensor, channel], 1000.0),
                                converted_or_na(uncertainty_values[selected_index, sensor, channel], 1000.0),
                                int(valid_values[selected_index, sensor, channel]),
                                int(dead_time_values[selected_index, sensor, channel]),
                                int(out_of_band_values[selected_index, sensor, channel]),
                                int(dynamic_error_values[selected_index, sensor, channel]),
                            )
                        )

    return 0


if __name__ == "__main__":
    sys.exit(main())

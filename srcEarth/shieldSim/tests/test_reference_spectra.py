#!/usr/bin/env python3
"""Offline integrity and semantics tests for published radiation fixtures."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import math
import pathlib
import re
import sys


ROOT = pathlib.Path(__file__).resolve().parent / "reference_spectra"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def finite_nonnegative(values: list[str], label: str) -> list[float]:
    converted = [float(value) for value in values]
    require(all(math.isfinite(value) and value >= 0.0 for value in converted), label)
    return converted


def test_manifest() -> None:
    entries = []
    for line in (ROOT / "manifest.sha256").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(None, 1)
        path = ROOT / relative.strip()
        require(path.is_file(), f"manifest file missing: {relative}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        require(actual == expected, f"SHA-256 mismatch for {relative}")
        entries.append(relative)
    require(len(entries) == 9, "manifest must pin exactly nine fixture and source files")


def parameter(metadata: dict, name: str) -> dict:
    matches = [item for item in metadata["parameters"] if item["name"] == name]
    require(len(matches) == 1, f"expected one HAPI parameter named {name}")
    return matches[0]


def test_ams02() -> None:
    info_path = ROOT / "gcr" / "ams02_20131118_hapi_info.json"
    data_path = ROOT / "gcr" / "ams02_20131118_hapi_data.csv"
    metadata = json.loads(info_path.read_text(encoding="utf-8"))
    require(metadata["status"]["code"] == 1200, "AMS-02 HAPI status is not OK")

    proton = parameter(metadata, "ProFlux")
    proton_error = parameter(metadata, "ProFlux_TotErr")
    helium = parameter(metadata, "HeFlux")
    helium_error = parameter(metadata, "HeFlux_TotErr")
    expected_units = "1/[m**2-sr-s-GeV/n]"
    for item in (proton, proton_error, helium, helium_error):
        require(item["units"] == expected_units, f"unexpected units for {item['name']}")
    require(proton["size"] == proton_error["size"] == [45], "proton width changed")
    require(helium["size"] == helium_error["size"] == [40], "helium width changed")
    require(len(proton["bins"][0]["centers"]) == 45, "proton coordinate count changed")
    require(len(helium["bins"][0]["centers"]) == 40, "helium coordinate count changed")

    with data_path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    require(len(rows) == 1, "AMS-02 fixture must contain exactly one time record")
    row = rows[0]
    require(len(row) == 172, "AMS-02 selected-field width changed")
    require(row[0] == "2013-11-18T00:00:00.000000000Z", "AMS-02 epoch changed")
    require(row[1] == "", "AMS-02 fixture unexpectedly populated the HAPI BR scalar")
    # Field order from the recorded HAPI request: Time, BR, p, p error, He, He error.
    p_flux = finite_nonnegative(row[2:47], "invalid AMS-02 proton flux")
    p_error = finite_nonnegative(row[47:92], "invalid AMS-02 proton uncertainty")
    he_flux = finite_nonnegative(row[92:132], "invalid AMS-02 helium flux")
    he_error = finite_nonnegative(row[132:172], "invalid AMS-02 helium uncertainty")
    require(all(value > 0.0 for value in p_flux + he_flux), "AMS-02 flux has nonpositive data")
    require(all(value > 0.0 for value in p_error + he_error), "AMS-02 errors are not retained")

    interpretation = json.loads(
        (ROOT / "gcr" / "ams02_20131118_interpretation.v1.json").read_text(encoding="utf-8")
    )
    require(
        interpretation["schema_version"] == "shieldSim.ams02_hapi_interpretation.v1",
        "AMS interpretation schema changed",
    )
    require(interpretation["bartels_rotation"] == 2460, "AMS supplement BR changed")
    proton_edges = interpretation["proton_rigidity_edges_GV"]
    helium_edges = interpretation["helium_rigidity_edges_GV"]
    require(len(proton_edges) == 46 and len(helium_edges) == 41, "AMS edge count changed")
    require(
        (ROOT / "gcr" / "raw" / "SM-p-He-vs-BR-REV.pdf").read_bytes().startswith(b"%PDF"),
        "AMS primary supplement is not a PDF",
    )

    # Primary Table SM XXXV reports flux per GV. CDAWeb reports flux per
    # GeV/n, so the bin-averaged differential value transforms by dR/dT. This
    # check independently proves that the HAPI arrays called "centers" are the
    # transformed lower edges named BinProLow/BinHeLow, not bin centers.
    mass = interpretation["nucleon_rest_mass_GeV"]

    def kinetic_per_nucleon(rigidity: float, charge_over_mass: float) -> float:
        momentum = rigidity * charge_over_mass
        return math.sqrt(momentum * momentum + mass * mass) - mass

    p_width = kinetic_per_nucleon(1.16, 1.0) - kinetic_per_nucleon(1.00, 1.0)
    he_width = kinetic_per_nucleon(2.15, 0.5) - kinetic_per_nucleon(1.92, 0.5)
    require(math.isclose(p_flux[0], 395.0 * 0.16 / p_width, rel_tol=1e-3), "AMS proton Jacobian")
    require(math.isclose(he_flux[0], 46.58 * 0.23 / he_width, rel_tol=2e-3), "AMS helium Jacobian")


def test_goes16() -> None:
    path = ROOT / "sep" / "goes16_sgps_20170911T0730Z.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(line for line in stream if not line.startswith("#")))
    require(len(rows) == 26, "GOES-16 fixture must contain 2 sensors x 13 channels")
    require({row["sensor_unit"] for row in rows} == {"0", "1"}, "GOES directions lost")
    require({row["spacecraft_axis"] for row in rows} == {"-X", "+X"}, "GOES axes lost")
    require(
        {row["upright_nominal_look"] for row in rows} == {"westward", "eastward"},
        "GOES nominal look directions lost",
    )
    require(all(row["sample_start_utc"] == "2017-09-11T07:30:00Z" for row in rows),
            "GOES start timestamp changed")
    require(all(row["sample_midpoint_utc"] == "2017-09-11T07:32:30Z" for row in rows),
            "GOES midpoint timestamp changed")
    require(all(row["sample_end_utc"] == "2017-09-11T07:35:00Z" for row in rows),
            "GOES end timestamp changed")
    require(all(row["l1b_records_in_average"] == "300" for row in rows),
            "GOES record-count completeness changed")
    require(all(row["yaw_flip_flag"] == "0" for row in rows), "GOES yaw flag changed")
    for row in rows:
        low = float(row["energy_low_MeV"])
        high = float(row["energy_high_MeV"])
        flux = float(row["proton_flux_per_cm2_sr_s_MeV"])
        uncertainty = float(row["flux_uncertainty_per_cm2_sr_s_MeV"])
        require(0.0 < low < high, f"invalid GOES energy interval in {row['channel']}")
        require(math.isfinite(flux) and flux > 0.0, f"invalid GOES flux in {row['channel']}")
        require(
            math.isfinite(uncertainty) and uncertainty > 0.0,
            f"invalid GOES uncertainty in {row['channel']}",
        )
        require(row["valid_1s_samples"] == "300", "GOES sample completeness changed")
        for name in (
            "dead_time_flagged_samples",
            "out_of_band_flagged_samples",
            "dynamic_error_flagged_samples",
        ):
            require(row[name] == "0", f"GOES quality flag changed: {name}")


def parse_utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    require(parsed.utcoffset() == dt.timedelta(0), f"timestamp is not UTC: {value}")
    return parsed


def test_goes16_event() -> None:
    definition_path = ROOT / "sep" / "goes16_gle72_event_definition.v1.json"
    csv_path = ROOT / "sep" / "goes16_sgps_gle72_timeseries.csv"
    raw_path = (
        ROOT
        / "sep"
        / "raw"
        / "se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc"
    )
    definition = json.loads(definition_path.read_text(encoding="utf-8"))
    require(
        definition["schema_version"] == "shieldSim.goes_sgps_event_definition.v1",
        "GOES event schema changed",
    )
    require(raw_path.read_bytes()[:8] == b"\x89HDF\r\n\x1a\n", "GOES raw source is not HDF5/netCDF4")
    require(
        "sha256:" + hashlib.sha256(raw_path.read_bytes()).hexdigest()
        == definition["source_file"]["sha256"],
        "GOES raw source hash differs from event definition",
    )

    windows = {item["id"]: item for item in definition["windows"]}
    require(set(windows) == {"gle72_gt100MeV_primary", "spe_gt10MeV_sensitivity"}, "GOES windows")
    require(
        windows["gle72_gt100MeV_primary"]["start_utc"] == "2017-09-10T16:25:00Z"
        and windows["gle72_gt100MeV_primary"]["end_utc"] == "2017-09-12T22:40:00Z",
        "GOES primary window changed",
    )
    require(
        windows["spe_gt10MeV_sensitivity"]["start_utc"] == "2017-09-10T16:45:00Z"
        and windows["spe_gt10MeV_sensitivity"]["end_utc"] == "2017-09-14T18:50:00Z",
        "GOES sensitivity window changed",
    )
    require(
        definition["markers"] == [
            {
                "id": "NMDB_GLE72_onset",
                "time_utc": "2017-09-10T16:08:00Z",
                "role": "independent_onset_marker_not_an_integration_boundary",
                "source_url": "https://www.nmdb.eu/nest/gle_list.php",
            }
        ],
        "NMDB onset marker changed",
    )

    with csv_path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(line for line in stream if not line.startswith("#")))
    require(len(rows) == 1181 * 2 * 13, "GOES event row count changed")
    keys: set[tuple[str, str, str]] = set()
    starts: set[str] = set()
    channel_bounds: dict[tuple[str, str], tuple[str, str]] = {}
    for row in rows:
        start = parse_utc(row["sample_start_utc"])
        midpoint = parse_utc(row["sample_midpoint_utc"])
        end = parse_utc(row["sample_end_utc"])
        require(midpoint - start == dt.timedelta(seconds=150), "GOES midpoint shifted")
        require(end - start == dt.timedelta(seconds=300), "GOES sample duration changed")
        key = (row["sample_start_utc"], row["sensor_unit"], row["channel"])
        require(key not in keys, f"duplicate GOES event row {key}")
        keys.add(key)
        starts.add(row["sample_start_utc"])
        component = (row["sensor_unit"], row["channel"])
        bounds = (row["energy_low_MeV"], row["energy_high_MeV"])
        require(component not in channel_bounds or channel_bounds[component] == bounds, "GOES bounds changed")
        channel_bounds[component] = bounds
    require(len(starts) == 1181 and len(channel_bounds) == 26, "GOES event grid is incomplete")
    require({row["yaw_flip_flag"] for row in rows} == {"0"},
            "GOES event no longer has a single upright yaw state")

    # Independent oracle for one published channel. It uses the same approved
    # quality predicate but no production C++ code.
    primary = windows["gle72_gt100MeV_primary"]
    window_start = parse_utc(primary["start_utc"])
    window_end = parse_utc(primary["end_utc"])
    fluence = valid_seconds = rejected_seconds = 0.0
    for row in rows:
        if row["sensor_unit"] != "0" or row["channel"] != "P7":
            continue
        overlap = max(
            0.0,
            (min(parse_utc(row["sample_end_utc"]), window_end)
             - max(parse_utc(row["sample_start_utc"]), window_start)).total_seconds(),
        )
        if not overlap:
            continue
        usable = (
            row["l1b_records_in_average"] == "300"
            and row["valid_1s_samples"] == "300"
            and row["dead_time_flagged_samples"] == "0"
            and row["out_of_band_flagged_samples"] == "0"
            and row["dynamic_error_flagged_samples"] == "0"
            and row["proton_flux_per_cm2_sr_s_MeV"] != "NA"
            and row["flux_uncertainty_per_cm2_sr_s_MeV"] != "NA"
        )
        if usable:
            valid_seconds += overlap
            fluence += float(row["proton_flux_per_cm2_sr_s_MeV"]) * overlap
        else:
            rejected_seconds += overlap
    require(math.isclose(fluence, 402238.9858300105, rel_tol=1e-13), "GOES fluence oracle changed")
    require(valid_seconds == 193800.0 and rejected_seconds == 1500.0, "GOES time ledger changed")


def test_crater() -> None:
    path = ROOT / "observations" / "crater_gcr_let_spectra_v1_3.csv"
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.reader(stream))
    require(len(rows) == 292, "CRaTER row count changed")
    require(len(rows[0]) == 35, "CRaTER must contain LET plus 34 time intervals")
    require(rows[0][0] == "LET (keV / micron)", "CRaTER LET unit changed")
    interval_pattern = re.compile(r"^20\d{5}-20\d{5}$")
    require(
        all(interval_pattern.match(value) for value in rows[0][1:]),
        "CRaTER interval labels changed",
    )
    require(all(len(row) == 35 for row in rows[1:]), "CRaTER rectangular shape changed")
    let_values = finite_nonnegative([row[0] for row in rows[1:]], "invalid CRaTER LET")
    require(all(a < b for a, b in zip(let_values, let_values[1:])), "CRaTER LET is not ordered")
    spectra = [float(value) for row in rows[1:] for value in row[1:]]
    require(all(math.isfinite(value) and value >= 0.0 for value in spectra), "invalid CRaTER flux")
    require(any(value > 0.0 for value in spectra), "CRaTER spectra are empty")


def main() -> int:
    tests = (test_manifest, test_ams02, test_goes16, test_goes16_event, test_crater)
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"All {len(tests)} published-reference tests passed.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"FAIL {exc}", file=sys.stderr)
        sys.exit(1)

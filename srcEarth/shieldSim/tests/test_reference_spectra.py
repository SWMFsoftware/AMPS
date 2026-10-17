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
    require(len(entries) == len(set(entries)), "manifest contains a duplicate path")
    expected_files = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.name not in {"README.md", "manifest.sha256"}
    }
    require(set(entries) == expected_files, "manifest must pin every fixture and source file exactly once")


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


def test_ams02_solar_cycle() -> None:
    directory = ROOT / "gcr" / "ams02_2025"
    interpretation = json.loads((directory / "interpretation.v1.json").read_text(encoding="utf-8"))
    require(
        interpretation["schema_version"] == "shieldSim.ams02_solar_cycle_interpretation.v1",
        "AMS solar-cycle interpretation schema changed",
    )
    require(len(interpretation["files"]) == 7, "AMS solar-cycle element count changed")
    common_rotations: list[int] | None = None
    common_grid: list[tuple[float, float]] | None = None
    for definition in interpretation["files"]:
        path = directory / definition["path"]
        require("sha256:" + hashlib.sha256(path.read_bytes()).hexdigest() == definition["sha256"],
                f"AMS solar-cycle hash mismatch for {definition['symbol']}")
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        require(len(rows) == 5880, f"AMS {definition['symbol']} row count changed")
        rotations = sorted({int(row["bartels_rotation_number"]) for row in rows})
        require(len(rotations) == 147, f"AMS {definition['symbol']} rotation count changed")
        counts = {rotation: 0 for rotation in rotations}
        grids: dict[int, list[tuple[float, float]]] = {rotation: [] for rotation in rotations}
        prefix = definition["published_name"]
        uncertainty_fields = (
            f"{prefix}_flux_error_statistical m^-2sr^-1s^-1GV^-1",
            f"{prefix}_flux_error_time_dependent m^-2sr^-1s^-1GV^-1",
            f"{prefix}_flux_error_time_independent m^-2sr^-1s^-1GV^-1",
        )
        for row in rows:
            rotation = int(row["bartels_rotation_number"])
            counts[rotation] += 1
            low, high = float(row["rigidity_min_GV"]), float(row["rigidity_max_GV"])
            flux = float(row[f"{prefix}_flux m^-2sr^-1s^-1GV^-1"])
            require(0.0 < low < high and math.isfinite(flux) and flux > 0.0,
                    f"invalid AMS {definition['symbol']} bin")
            require(all(math.isfinite(float(row[name])) and float(row[name]) >= 0.0
                        for name in uncertainty_fields),
                    f"invalid AMS {definition['symbol']} uncertainty")
            grids[rotation].append((low, high))
        require(set(counts.values()) == {40}, f"AMS {definition['symbol']} bin count changed")
        require(all(grid == grids[rotations[0]] for grid in grids.values()),
                f"AMS {definition['symbol']} grid varies with time")
        if common_rotations is None:
            common_rotations, common_grid = rotations, grids[rotations[0]]
        else:
            require(rotations == common_rotations and grids[rotations[0]] == common_grid,
                    f"AMS {definition['symbol']} time/grid coverage differs")

    require(common_rotations is not None and common_grid is not None, "AMS solar-cycle data absent")
    require(common_grid[0] == (1.92, 2.15) and common_grid[-1] == (56.1, 60.3),
            "AMS solar-cycle rigidity endpoints changed")
    partitions = json.loads(
        (ROOT / "gcr" / "validation_partitions.v1.json").read_text(encoding="utf-8")
    )["ams02_2025"]
    development = [r for r in common_rotations
                   if partitions["development"]["first"] <= r <= partitions["development"]["last"]]
    holdout = [r for r in common_rotations
               if partitions["holdout"]["first"] <= r <= partitions["holdout"]["last"]]
    require(len(development) == partitions["development"]["present_rotations"] == 118,
            "AMS development partition changed")
    require(len(holdout) == partitions["holdout"]["present_rotations"] == 29,
            "AMS holdout partition changed")
    require(development + holdout == common_rotations, "AMS partitions do not cover the data exactly")


def test_ace_cris() -> None:
    directory = ROOT / "gcr" / "ace_cris"
    interpretation = json.loads((directory / "interpretation.v1.json").read_text(encoding="utf-8"))
    require(
        interpretation["schema_version"] == "shieldSim.ace_cris_hapi_interpretation.v1",
        "ACE/CRIS interpretation schema changed",
    )
    info = json.loads((directory / "ac_h3_cris_hapi_info.json").read_text(encoding="utf-8"))
    require(info["status"]["code"] == 1200, "ACE/CRIS HAPI status is not OK")
    symbols = ["B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
               "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
    require(list(interpretation["species"]) == symbols,
            "ACE/CRIS species order must remain B through Ni")
    for symbol in symbols:
        flux, count = parameter(info, f"flux_{symbol}"), parameter(info, f"cnt_{symbol}")
        require(flux["size"] == count["size"] == [7], f"ACE/CRIS {symbol} width changed")
        require(flux["units"] == "1/(cm2 Sr sec MeV/nucleon)",
                f"ACE/CRIS {symbol} flux unit changed")
        require(count["units"] == "Counts/day", f"ACE/CRIS {symbol} count unit changed")
        bins = interpretation["species"][symbol]
        require(len(bins["bin_lower_MeV_per_nucleon"]) == 7
                and len(bins["bin_upper_MeV_per_nucleon"]) == 7,
                f"ACE/CRIS {symbol} bounds changed")
        require(all(low < high for low, high in zip(
            bins["bin_lower_MeV_per_nucleon"], bins["bin_upper_MeV_per_nucleon"])),
            f"ACE/CRIS {symbol} has invalid bounds")

    with (directory / "ac_h3_cris_20191201_all.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    require(len(rows) == 1 and len(rows[0]) == 337, "ACE/CRIS HAPI fixture shape changed")
    row = rows[0]
    require(row[0] == interpretation["record_utc"], "ACE/CRIS record time changed")
    fluxes = finite_nonnegative(row[1:169], "invalid ACE/CRIS flux")
    counts = finite_nonnegative(row[169:337], "invalid ACE/CRIS count")
    require(all(value.is_integer() for value in counts), "ACE/CRIS counts are not integral")
    require(all(flux == 0.0 for flux, count in zip(fluxes, counts) if count == 0.0),
            "ACE/CRIS zero-count bin has nonzero flux")
    require(sum(count == 0.0 for count in counts) == 19, "ACE/CRIS zero-count mask changed")
    oxygen_offset = 3 * 7
    require(fluxes[oxygen_offset] == 7.375e-7 and counts[oxygen_offset] == 139.0,
            "ACE/CRIS oxygen spot value changed")
    require(math.isclose(fluxes[oxygen_offset] / math.sqrt(counts[oxygen_offset]),
                         6.255393356389787e-8, rel_tol=1e-14),
            "ACE/CRIS count-derived statistical uncertainty changed")

    cdf = directory / interpretation["cdf"]["path"]
    require(cdf.read_bytes()[:4] == b"\xcd\xf3\x00\x01", "ACE/CRIS source is not a CDF3 file")
    hdf = directory / "raw" / "cris_data_Bartels.hdf"
    require(hdf.read_bytes()[:4] == b"\x0e\x03\x13\x01", "ACE/CRIS source is not HDF4")
    oxygen = interpretation["species"]["O"]
    require(oxygen["bin_lower_MeV_per_nucleon"][0] == 66.3
            and oxygen["bin_upper_MeV_per_nucleon"][-1] == 235.8,
            "ACE/CRIS CDF oxygen labels changed")
    factors = (directory / "raw" / "cris_factors.txt").read_text(encoding="utf-8")
    require(" 8      2      69.4    89.0" in factors,
            "ACE current factors documentation changed")
    require("mixing it with the V04 HAPI/CDF values is prohibited"
            in interpretation["energy_bounds_policy"],
            "ACE/CRIS energy-version mismatch is no longer guarded")
    require(interpretation["systematic_uncertainty_policy"].startswith("not_applied:"),
            "ACE/CRIS undocumented systematic uncertainty was enabled")


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
    tests = (
        test_manifest,
        test_ams02,
        test_ams02_solar_cycle,
        test_ace_cris,
        test_goes16,
        test_goes16_event,
        test_crater,
    )
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

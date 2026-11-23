#!/usr/bin/env python3
"""Reproduce the 1997 CDS/Ulysses observation-side connectivity screen."""

import argparse
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
RSUN_M = 695700000.0  # IAU nominal photospheric radius
AU_M = 149597870700.0
SHELL_RSUN = 1.05
RSS_RSUN = 2.5  # Source surface in Suess et al. Figure 4.
RUNS = [(12, 0.02), (12, 0.01), (8, 0.01)]
ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi
NATIVE_TYPES = {"I": ">i2", "D": ">f8", "E": ">f4", "J": ">i4", "B": "u1", "L": "S1"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def utc(text):
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def iso(epoch):
    return epoch.astimezone(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_inputs(folder):
    manifest_path = HERE / "retrieval.json"
    manifest = json.loads(manifest_path.read_bytes())
    for name, asset in manifest["products"].items():
        path = folder / name
        require(digest(path) == asset["sha256"],
                f"{path} differs from the pinned archive product")
    return manifest


def fits_header(raw, offset=0):
    # This reads the fixed FITS headers of the pinned NSO image and the two
    # historical CDS tables. It accepts their explicit metadata, including
    # the preliminary CDS coordinate cards. It is not a FITS/WCS calibrator.
    cards = {}
    for i in range(offset, len(raw), 80):
        card = raw[i:i + 80].decode("ascii")
        key = card[:8].strip()
        if key == "END":
            return cards, (i + 80 + 2879) // 2880 * 2880
        if card[8:10] == "= ":
            cards[key] = card[10:].split(" /", maxsplit=1)[0].strip().strip("'").strip()
    raise ValueError("missing FITS END card")


def radial_map(path):
    raw = path.read_bytes()
    header, offset = fits_header(raw)
    require((header["BITPIX"], header["NAXIS1"], header["NAXIS2"]) == ("16", "360", "180")
            and header["IRAFNAME"] == "m1923f", "unexpected KPVT map layout")
    stored = np.frombuffer(raw, dtype=">i2", count=360 * 180, offset=offset)
    require("BLANK" not in header or not np.any(stored == int(header["BLANK"])),
            "magnetogram includes missing values")
    return float(header["BSCALE"]) * stored + float(header["BZERO"])


def parse_tuple(value):
    return tuple(float(item) for item in value.strip("()").split(","))


def raster(path):
    raw = path.read_bytes()
    primary, offset = fits_header(raw)
    table, offset = fits_header(raw, offset)
    require(primary["NAXIS"] == "0" and primary["SEQVALID"] == "T"
            and primary["DETECTOR"] == "NIS" and table["EXTNAME"] == "DATA"
            and table["NAXIS2"] == "1", "unexpected CDS raster format")
    row_end = offset + int(table["NAXIS1"])
    cols, windows, geometry = {}, [], None
    for index in range(1, int(table["TFIELDS"]) + 1):
        number, kind = re.fullmatch(r"(\d+)([IDEJBL])", table[f"TFORM{index}"]).groups()
        array = np.frombuffer(raw, dtype=NATIVE_TYPES[kind], count=int(number), offset=offset)
        offset += array.nbytes
        name, units = table[f"TTYPE{index}"], table[f"TUNIT{index}"]
        if units == "ADC":
            # ADC is detector signal. Counts at saturation and the FITS
            # missing sentinel remain metadata rather than line intensities.
            windows.append({"name": name, "principal_angstrom": float(table[f"TWAVE{index}"]),
                            "samples": int(number), "missing": int(np.count_nonzero(
                                array == int(table[f"TNULL{index}"]))),
                            "samples_at_12_bit_ceiling": int(np.count_nonzero(array == 4095))})
        else:
            cols[name] = array
        if name == "MG_9_368_07":
            geometry = tuple(parse_tuple(table[f"{key}{index}"])
                             for key in ("TDIM", "TRPIX", "TRVAL", "TDELT", "TROTA"))
    require(offset == row_end and geometry is not None, "CDS table row or geometry differs")
    require(sum(win["name"] == "BACKGROUND" for win in windows) == 4,
            "expected all four CDS detector bias windows")
    dims, refpixel, refvalue, increment, rotation = geometry
    nx, ny = int(dims[1]), int(dims[2])
    # Wavelength is the first FITS axis. Raster columns run right to left
    # in time even though helioprojective X increases with column index.
    # Using the raster midtime for every column would blur the footpoint
    # longitude. The preliminary legacy geometry remains a diagnostic.
    xs = (np.arange(nx) + 1 - refpixel[1]) * increment[1]
    ys = (np.arange(ny) + 1 - refpixel[2]) * increment[2]
    angle = math.radians(rotation[1])
    require(rotation[1] == rotation[2], "unexpected rotation between CDS sky axes")
    start, stop = utc(primary["DATE_OBS"]), utc(primary["DATE_END"])
    exposure_sec = float(primary["EXPTIME"])
    epochs = [start + timedelta(seconds=float(dt) + exposure_sec / 2) for dt in cols["DEL_TIME"]]
    require(len(epochs) == nx and min(epochs) >= start and max(epochs) <= stop,
            "CDS exposure midpoints outside native raster interval")
    xs, ys = np.meshgrid(xs, ys, indexing="ij")
    x_arcsec = refvalue[1] + math.cos(angle) * xs - math.sin(angle) * ys
    y_arcsec = refvalue[2] + math.sin(angle) * xs + math.cos(angle) * ys
    return epochs, x_arcsec, y_arcsec, windows


def time_rows(path, columns):
    with path.open(newline="") as f:
        rows = [(utc(time), *(float(v) for v in values))
                for time, *values in csv.reader(f)]
    require(rows and all(len(row) == columns + 1 for row in rows)
            and all(a[0] < b[0] for a, b in zip(rows, rows[1:])),
            f"{path}: unexpected column count or time order")
    return rows


def angle_delta(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def soho_at(epoch, rows):
    # Interpolate spacecraft position only, with native one-hour bracketing.
    # Longitude is cyclic. No plasma or magnetic data gaps are filled.
    for a, b in zip(rows, rows[1:]):
        if a[0] <= epoch <= b[0]:
            dt = (b[0] - a[0]).total_seconds()
            require(dt == 3600.0 and min(a[1], b[1]) > 0.0, "invalid SOHO ephemeris bracket")
            weight = (epoch - a[0]).total_seconds() / dt
            au = a[1] + weight * (b[1] - a[1])
            lat = a[2] + weight * (b[2] - a[2])
            carr = a[3] + weight * angle_delta(b[3], a[3])
            hgi = a[5] + weight * angle_delta(b[5], a[5])
            return au, lat, carr % 360.0, hgi % 360.0
    return None


def ray_to_shell(x_arcsec, y_arcsec, au, latitude_deg, longitude_deg):
    # The helioprojective ray originates at SOHO. Solar west is positive X
    # and projects along increasing Carrington longitude at disk center.
    # Use the true minimum-distance impact parameter, followed by both
    # intersections with the chosen sphere. Neither side is assigned an
    # emissivity weight without an independently constrained LOS model.
    x, y = x_arcsec / ARCSEC_PER_RAD, y_arcsec / ARCSEC_PER_RAD
    lat, lon = math.radians(latitude_deg), math.radians(longitude_deg)
    toward_observer = np.array([math.cos(lat) * math.cos(lon),
                                math.cos(lat) * math.sin(lon), math.sin(lat)])
    west = np.array([-math.sin(lon), math.cos(lon), 0.0])
    north = np.cross(toward_observer, west)
    direction = (-math.cos(x) * math.cos(y) * toward_observer
                 + math.sin(x) * math.cos(y) * west + math.sin(y) * north)
    origin = (au * AU_M / RSUN_M) * toward_observer
    nearest = origin - np.dot(origin, direction) * direction
    impact = float(np.linalg.norm(nearest))
    if not 1.0 <= impact < SHELL_RSUN:
        return impact, ()
    half_chord = math.sqrt(SHELL_RSUN ** 2 - impact ** 2)
    return impact, ((nearest - half_chord * direction) * RSUN_M,
                    (nearest + half_chord * direction) * RSUN_M)


def calibrated_positions(folder, path, shape, receipt):
    # The optional calibration run ties both its products to the identical
    # pinned raw exposure. A WCS grid without its QLDS (or a package copied
    # from another raster) cannot be substituted for preliminary positions.
    record = receipt["rasters"][path.name]
    table = folder / (path.stem + ".wcs.csv")
    radiance = folder / (path.stem + ".sav")
    require(record["source_sha256"] == digest(path) and record["wcs_sha256"] == digest(table)
            and record["radiance_sha256"] == digest(radiance), "calibration provenance mismatch")
    with table.open(newline="") as f:
        rows = list(csv.DictReader(f))
    require(len(rows) == math.prod(shape), "corrected WCS must cover each native pixel")
    positions = np.empty((*shape, 2))
    occupied = set()
    for row in rows:
        index = int(row["x_index_zero"]), int(row["y_index_zero"])
        require(all(0 <= i < length for i, length in zip(index, shape)) and index not in occupied,
                "invalid or duplicate index in corrected WCS")
        positions[index] = float(row["x_arcsec"]), float(row["y_arcsec"])
        occupied.add(index)
    require(np.isfinite(positions).all(), "corrected WCS has invalid coordinates")
    return positions[:, :, 0], positions[:, :, 1]


def seeds_from_rasters(paths, observer, calibration_dir, calibration_receipt):
    rays, seeds, raster_metadata = [], [], {}
    for path in paths:
        epochs, xs, ys, windows = raster(path)
        if calibration_dir is not None:
            xs, ys = calibrated_positions(calibration_dir, path, xs.shape, calibration_receipt)
        # Freeze the observing interval before selecting any open or closed
        # rays. A pointing correction can change which pixels cross the
        # shell, but it must not change the Ulysses arrival cohort searched.
        raster_metadata[path.name] = {
            "first_exposure_mid_utc": iso(min(epochs)), "last_exposure_mid_utc": iso(max(epochs)),
            "window_metadata": windows, "off_limb_shell_pixels": 0,
        }
        for ix, epoch in enumerate(epochs):
            ephemeris = soho_at(epoch, observer)
            require(ephemeris is not None, "CDS exposure lies outside SOHO ephemeris")
            au, lat, lon, _ = ephemeris
            for iy, (x, y) in enumerate(zip(xs[ix], ys[ix])):
                impact, points = ray_to_shell(x, y, au, lat, lon)
                if not points:
                    continue
                raster_metadata[path.name]["off_limb_shell_pixels"] += 1
                for side, p in zip(("near", "far"), points):
                    rays.append((path.name, ix, iy, iso(epoch), x, y, impact, side))
                    seeds.append(p)
    require(seeds, "no off-limb pixels intersect the requested physical shell")
    return rays, np.asarray(seeds), raster_metadata


def ballistic_sources(rows, observer, first_epoch, last_epoch):
    candidates = []
    for arrival, range_au, lat, inertial_lon, br, speed, elev, azimuth in rows:
        # HAPI numeric fill is -1e31. Missing moments and signed sector
        # remain missing. Use radial velocity, not the RTN vector magnitude.
        if min(range_au, speed) <= 0 or any(abs(v) >= 1.0e30 or not math.isfinite(v)
                for v in (range_au, speed, lat, inertial_lon, br, elev, azimuth)) or br == 0:
            continue
        vr = speed * math.cos(math.radians(elev)) * math.cos(math.radians(azimuth))
        if vr <= 0.0:
            continue
        lag = (range_au * AU_M - RSS_RSUN * RSUN_M) / (vr * 1000.0)
        departure = arrival - timedelta(seconds=lag)
        # The ±12 h envelope is a coverage diagnostic around the raster
        # interval, not a same-parcel gate. Unknown travel time below the
        # source surface and departures near discontinuities remain visible.
        if not first_epoch - timedelta(hours=12) <= departure <= last_epoch + timedelta(hours=12):
            continue
        ephemeris = soho_at(departure, observer)
        require(ephemeris is not None, "ballistic departure outside SOHO ephemeris")
        _, _, carr_soho, hgi_soho = ephemeris
        lon = (carr_soho + angle_delta(inertial_lon, hgi_soho)) % 360.0
        candidates.append((iso(arrival), iso(departure), vr, br, lat, lon, lag / 86400.0))
    return candidates


def compile_trace(executable):
    models = ROOT / "src/models"
    core = models / "sep_coronal_cme"
    subprocess.run([os.environ.get("CXX", "c++"), "-std=c++17", "-O3", "-Wall",
                    "-Wextra", "-Werror", "-fopenmp", "-I", str(core / "include"),
                    "-I", str(models / "sep_common"), str(HERE / "trace_pfss.cpp"),
                    str(core / "src/pfss_harmonics.cpp"), "-o", str(executable)], check=True)


def trace(executable, map_g, seeds, degree, step):
    with tempfile.TemporaryFile(mode="w+") as source:
        source.write(f"{degree} {step} {RSS_RSUN} {len(seeds)}\n")
        np.savetxt(source, map_g.reshape(180, 360), fmt="%.10g")
        np.savetxt(source, seeds, fmt="%.17g")
        source.seek(0)
        proc = subprocess.run([str(executable)], stdin=source, stdout=subprocess.PIPE,
                              text=True, check=True)
    first, *lines = proc.stdout.splitlines()
    require(len(lines) == len(seeds), "native trace did not finish all requested rays")
    return float(first), [(parts[0], *(float(v) for v in parts[1:]))
                          for parts in map(str.split, lines)]


def coordinates(x, y, z):
    return math.degrees(math.atan2(z, math.hypot(x, y))), math.degrees(math.atan2(y, x)) % 360.0


def angular_distance(lat, lon, lat2, lon2):
    # Haversine avoids inverse-cosine cancellation for nearby footpoints.
    a = math.sin(math.radians(lat - lat2) / 2.0) ** 2
    a += math.cos(math.radians(lat)) * math.cos(math.radians(lat2)) * (
        math.sin(math.radians(angle_delta(lon, lon2)) / 2.0) ** 2)
    return math.degrees(2.0 * math.asin(min(1.0, math.sqrt(a))))


def write_trace(path, rays, traces, candidates):
    metrics = {"open": 0, "closed": 0, "unresolved": 0, "closest_same_sign_deg": None}
    with path.open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(("raster", "x_index_zero", "y_index_zero", "exposure_mid_utc",
                         "x_arcsec", "y_arcsec", "impact_rsun", "shell_intersection",
                         "topology", "source_lat_deg", "source_carrington_lon_deg",
                         "source_br_T", "closest_same_sign_deg", "ulysses_arrival_utc"))
        for ray, (topology, x, y, z, br) in zip(rays, traces):
            metrics[topology] += 1
            if topology != "open":
                writer.writerow((*ray, topology, "", "", "", "", ""))
                continue
            lat, lon = coordinates(x, y, z)
            matches = [(angular_distance(lat, lon, candidate[4], candidate[5]), candidate[0])
                       for candidate in candidates if br * candidate[3] > 0.0]
            nearest = min(matches, default=("", ""))
            writer.writerow((*ray, topology, lat, lon, br, *nearest))
            gap = metrics["closest_same_sign_deg"]
            if matches and (gap is None or nearest[0] < gap):
                metrics["closest_same_sign_deg"] = nearest[0]
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True, help="directory of pinned source files")
    parser.add_argument("--output-dir", type=Path, required=True, help="fresh directory for derived evidence")
    parser.add_argument("--calibration-dir", type=Path,
                        help="optional SolarSoft radiance and corrected WCS package with manifest")
    args = parser.parse_args()

    metadata = read_inputs(args.input_dir)
    raster_paths = sorted(args.input_dir / name for name in metadata["products"] if name.startswith("s7997"))
    soho = time_rows(args.input_dir / "soho_pos.csv", 5)
    ulysses = time_rows(args.input_dir / "uy_ctx.csv", 7)
    calibration_receipt = None
    if args.calibration_dir is not None:
        calibration_receipt = json.loads((args.calibration_dir / "calibration_receipt.json").read_bytes())
        require(calibration_receipt["software"] and calibration_receipt["database_sha256"],
                "calibration package requires SolarSoft version and database hashes")
    rays, seeds, rasters = seeds_from_rasters(raster_paths, soho, args.calibration_dir, calibration_receipt)
    first = min(utc(item["first_exposure_mid_utc"]) for item in rasters.values())
    last = max(utc(item["last_exposure_mid_utc"]) for item in rasters.values())
    sources = ballistic_sources(ulysses, soho, first, last)
    map_g = radial_map(args.input_dir / "m1923f.fits")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    with (output / "ulysses_ballistic_candidates.csv").open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(("arrival_utc", "departure_utc_at_source_surface", "vr_km_s",
                         "br_rtn_nT", "source_lat_deg", "source_carrington_lon_deg", "travel_days"))
        writer.writerows(sources)

    results = []
    with tempfile.TemporaryDirectory() as scratch:
        executable = Path(scratch) / "trace_pfss"
        compile_trace(executable)
        previous = None
        for degree, step in RUNS:
            mean, traces = trace(executable, map_g, seeds, degree, step)
            path = output / f"rays_l{degree}_step{step:g}.csv"
            entry = write_trace(path, rays, traces, sources)
            entry.update(degree=degree, step_rsun=step, removed_monopole_gauss=mean * 1.0e4,
                         table=path.name)
            if previous is not None:
                entry["topology_changes_from_previous"] = sum(
                    a[0] != b[0] for a, b in zip(previous, traces))
                entry["max_open_endpoint_delta_deg_from_previous"] = max((
                    angular_distance(*coordinates(*a[1:4]), *coordinates(*b[1:4]))
                    for a, b in zip(previous, traces) if a[0] == b[0] == "open"), default=None)
            previous = traces
            results.append(entry)
    receipt = {"schema": "srcsep-stage13-cds1997-screen-v1", "geometry_status": "PASS_DIAGNOSTIC",
               "wind_qualification": "SKIPPED",
               "radiometry": ("REQUIRES_SOLARSOFT_CALIBRATION" if calibration_receipt is None
                              else "SOLARSOFT_PACKAGE_REGISTERED_LINE_FITS_PENDING"),
               "emission_attribution": ("REQUIRES_CALIBRATED_POINTING_AND_LOS_MODEL"
                                        if calibration_receipt is None else "REQUIRES_LOS_MODEL"),
               "pointing": "legacy preliminary WCS" if calibration_receipt is None else "SolarSoft corrected WCS",
               "calibration": calibration_receipt,
               "departure_interpretation": "constant radial speed from source surface only",
               "raster_interval": "Suess et al. I2, variable streamer",
               "source_surface_rsun": RSS_RSUN, "shell_rsun": SHELL_RSUN,
               "rasters": rasters, "ballistic_rows": len(sources), "trace_runs": results,
               "products": metadata["products"]}
    tmp = output / "result.json.partial"
    tmp.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    tmp.replace(output / "result.json")
    print(f"Geometry PASS_DIAGNOSTIC. Mapped {len(seeds)} ray intersections "
          f"and {len(sources)} ballistic rows. Wind qualification SKIPPED.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

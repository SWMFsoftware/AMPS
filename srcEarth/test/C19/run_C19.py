#!/usr/bin/env python3
"""C19A — public GOES EPEAD east-west directional-access validation.

C19A compares AMPS directional cutoff maps with the background-subtracted
physical EAST/WEST P4 and P5 proton-flux ratios measured simultaneously by
GOES-13 and GOES-15 during the 17 May 2012 SEP event decay.

The observational reference is created by ``build_goes_reference.py``.  AMPS is
run once for each selected (UTC, spacecraft, solver, field-model) combination.
A global SM directional cutoff map is folded through a documented top-hat
approximation to the EPEAD aperture and an incident power-law proton spectrum.
The primary comparison quantity is log10(EAST/WEST).

Examples
--------
Routine public-data comparison with T96 and T05::

    python3 srcEarth/test/C19/run_C19.py --profile ROUTINE \
      --solver GRIDDED --models T96,T05 \
      --reference srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz \
      --driver /path/to/may2012_driver.txt \
      --amps ./amps -np 4 -nt 16

Quick command/input preview::

    python3 srcEarth/test/C19/run_C19.py --profile SMOKE --dry-run \
      --driver /path/to/may2012_driver.txt

Exercise parsing, response folding, metrics, and plot generation without AMPS::

    python3 srcEarth/test/C19/run_C19.py --self-test

C19A is deliberately a broad-aperture observational test.  It does not claim a
full detector response simulation: the routine model uses nominal P4/P5 energy
bounds, a uniform elliptical aperture, and a common isotropic power-law source.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MANIFEST = SCRIPT_DIR / "event_C19_may2012.json"
DEFAULT_REFERENCE = SCRIPT_DIR / "data" / "reference_C19_goes_epead_ew.csv.gz"
DEFAULT_TEMPLATE_GRIDLESS = SCRIPT_DIR / "AMPS_PARAM_C19_gridless.in"
DEFAULT_TEMPLATE_MODE3D = SCRIPT_DIR / "AMPS_PARAM_C19_mode3d.in"

SOLVERS = ("GRIDLESS", "GRIDDED", "BOTH")
FIELD_MODELS = ("T96", "T05")
PROFILE_STEP_MINUTES = {"SMOKE": None, "ROUTINE": 60, "FULL": 0}


@dataclass(frozen=True)
class ReferenceRow:
    utc: datetime
    spacecraft: str
    channel: str
    energy_min_mev: float
    energy_max_mev: float
    east_west_ratio: float
    log10_east_west_ratio: float
    longitude_deg_east: float
    latitude_deg: float
    altitude_km: float
    position_source: str


@dataclass(frozen=True)
class DirectionCell:
    lon_deg: float
    lat_deg: float
    rc_gv: float
    cutoff_energy_mev: float


@dataclass(frozen=True)
class DirectionMap:
    path: str
    frame: str
    x_km: float
    y_km: float
    z_km: float
    cells: Tuple[DirectionCell, ...]


@dataclass(frozen=True)
class ModelRow:
    utc: str
    spacecraft: str
    channel: str
    solver: str
    field_model: str
    observed_east_west_ratio: float
    observed_log10_east_west_ratio: float
    modeled_east_west_ratio: Optional[float]
    modeled_log10_east_west_ratio: Optional[float]
    residual_log10: Optional[float]
    east_transmission: Optional[float]
    west_transmission: Optional[float]
    n_east_cells: int
    n_west_cells: int
    unresolved_direction_fraction: float
    spectral_index: float
    map_frame: str
    map_path: str
    status: str


@dataclass(frozen=True)
class Metrics:
    solver: str
    field_model: str
    channel: str
    n_reference: int
    n_valid_model: int
    valid_fraction: float
    sign_agreement_fraction: float
    mean_bias_log10: Optional[float]
    mean_absolute_error_log10: Optional[float]
    rmse_log10: Optional[float]
    correlation: Optional[float]
    passed: bool


def parse_utc(text: str) -> datetime:
    value = str(text).strip().replace("Z", "+00:00")
    dt = datetime.fromisoformat(value)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


def format_utc(dt: datetime, suffix_z: bool = True) -> str:
    value = dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    return value + ("Z" if suffix_z else "")


def timestamp_token(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_csv_list(text: str, allowed: Optional[Sequence[str]] = None) -> List[str]:
    values = []
    for token in str(text).replace(";", ",").split(","):
        value = token.strip().upper()
        if value and value not in values:
            values.append(value)
    if not values:
        raise ValueError("empty comma-separated selection")
    if allowed is not None:
        bad = [value for value in values if value not in allowed]
        if bad:
            raise ValueError("unsupported value(s) %s; allowed: %s" %
                             (",".join(bad), ",".join(allowed)))
    return values


def load_reference(path: Path) -> List[ReferenceRow]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    rows: List[ReferenceRow] = []
    with opener(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {
            "utc", "spacecraft", "channel", "energy_min_mev", "energy_max_mev",
            "east_west_ratio", "log10_east_west_ratio", "longitude_deg_east",
            "latitude_deg", "altitude_km", "position_source",
        }
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError("reference is missing columns: %s" % ", ".join(sorted(missing)))
        for record in reader:
            ratio = float(record["east_west_ratio"])
            log_ratio = float(record["log10_east_west_ratio"])
            if ratio <= 0.0 or not math.isfinite(ratio):
                continue
            if not math.isclose(log_ratio, math.log10(ratio), rel_tol=0.0, abs_tol=5.0e-9):
                raise ValueError("reference ratio/log ratio disagree: %s" % record)
            rows.append(ReferenceRow(
                utc=parse_utc(record["utc"]),
                spacecraft=record["spacecraft"].strip().upper(),
                channel=record["channel"].strip().upper(),
                energy_min_mev=float(record["energy_min_mev"]),
                energy_max_mev=float(record["energy_max_mev"]),
                east_west_ratio=ratio,
                log10_east_west_ratio=log_ratio,
                longitude_deg_east=float(record["longitude_deg_east"]),
                latitude_deg=float(record["latitude_deg"]),
                altitude_km=float(record["altitude_km"]),
                position_source=record["position_source"].strip(),
            ))
    if not rows:
        raise ValueError("reference contains no valid rows: %s" % path)
    rows.sort(key=lambda row: (row.utc, row.spacecraft, row.channel))
    return rows


def select_reference_rows(rows: Sequence[ReferenceRow], args: argparse.Namespace) -> List[ReferenceRow]:
    selected = [row for row in rows if row.spacecraft in args.spacecraft_list
                and row.channel in args.channel_list]
    if args.start:
        start = parse_utc(args.start)
        selected = [row for row in selected if row.utc >= start]
    if args.end:
        end = parse_utc(args.end)
        selected = [row for row in selected if row.utc <= end]
    if not selected:
        raise ValueError("no reference rows remain after spacecraft/channel/time selection")

    key_times = sorted({(row.utc, row.spacecraft) for row in selected})
    if args.profile == "SMOKE":
        keep = set()
        for spacecraft in args.spacecraft_list:
            times = sorted({epoch for epoch, sc in key_times if sc == spacecraft})
            if not times:
                continue
            keep.update((times[0], spacecraft) for _ in (0,))
            keep.update((times[len(times) // 2], spacecraft) for _ in (0,))
            keep.update((times[-1], spacecraft) for _ in (0,))
        selected = [row for row in selected if (row.utc, row.spacecraft) in keep]
    else:
        step_minutes = PROFILE_STEP_MINUTES[args.profile]
        if args.time_step_minutes is not None:
            step_minutes = args.time_step_minutes
        if step_minutes and step_minutes > 0:
            keep = set()
            for spacecraft in args.spacecraft_list:
                times = sorted({epoch for epoch, sc in key_times if sc == spacecraft})
                next_allowed: Optional[datetime] = None
                for epoch in times:
                    if next_allowed is None or epoch >= next_allowed:
                        keep.add((epoch, spacecraft))
                        next_allowed = epoch + timedelta(minutes=step_minutes)
            selected = [row for row in selected if (row.utc, row.spacecraft) in keep]
    if not selected:
        raise ValueError("profile selection removed all reference rows")
    return selected


def group_reference(rows: Sequence[ReferenceRow]) -> Dict[Tuple[datetime, str], List[ReferenceRow]]:
    grouped: Dict[Tuple[datetime, str], List[ReferenceRow]] = {}
    for row in rows:
        grouped.setdefault((row.utc, row.spacecraft), []).append(row)
    for value in grouped.values():
        value.sort(key=lambda row: row.channel)
    return dict(sorted(grouped.items()))


def kinetic_energy_mev_from_rigidity_gv(rigidity_gv: float) -> float:
    rest_mev = 938.27208816
    momentum = 1000.0 * rigidity_gv
    return math.sqrt(momentum * momentum + rest_mev * rest_mev) - rest_mev


def rigidity_gv_from_kinetic_energy_mev(energy_mev: float) -> float:
    rest_mev = 938.27208816
    momentum = math.sqrt(max(0.0, energy_mev * (energy_mev + 2.0 * rest_mev)))
    return momentum / 1000.0


def render_template(template: Path, destination: Path, replacements: Mapping[str, str]) -> None:
    """Render one solver-specific input template.

    ``replacements`` contains the union of GRIDLESS and GRIDDED placeholders.
    A solver-specific template is allowed to omit placeholders used only by the
    other solver; only placeholders that remain in the rendered file are an
    error.
    """
    text = template.read_text()
    for key, value in replacements.items():
        text = text.replace(key, value)
    leftovers = sorted(set(re.findall(r"__[A-Z0-9_]+__", text)))
    if leftovers:
        raise ValueError(
            "template contains unresolved placeholder(s): %s" % ", ".join(leftovers))
    destination.write_text(text)


def write_trajectory(path: Path, row: ReferenceRow) -> None:
    path.write_text("%s %.12g %.12g %.12g\n" % (
        format_utc(row.utc, suffix_z=False), row.latitude_deg,
        row.longitude_deg_east, row.altitude_km))


def resolved_dynamic_chunk(args: argparse.Namespace, solver: str) -> int:
    if args.dynamic_chunk > 0:
        return args.dynamic_chunk
    return 1 if solver == "GRIDLESS" else max(1, args.nt)


def command_for(args: argparse.Namespace, amps: Path, solver: str) -> List[str]:
    command = [
        args.mpirun, "-np", str(args.np), str(amps),
        "-mode", "gridless" if solver == "GRIDLESS" else "3d",
        "-i", "AMPS_PARAM_C19.in",
        "-mover", args.mover,
        "-cutoff-search", "UPPER_SCAN",
        "-cutoff-upper-scan-n", str(args.cutoff_scan_n),
    ]
    chunk = resolved_dynamic_chunk(args, solver)
    if solver == "GRIDLESS":
        command += [
            "-gridless-mpi-scheduler", args.scheduler,
            "-gridless-mpi-dynamic-chunk", str(chunk),
            "-density-parallel", "THREADS",
            "-density-threads", str(args.nt),
        ]
    else:
        command += [
            "-mode3d-field-eval", "MESH",
            "-mode3d-parallel", "THREADS",
            "-mode3d-threads", str(args.nt),
            "-mode3d-mpi-scheduler", args.scheduler,
            "-mode3d-mpi-dynamic-chunk", str(chunk),
            "-mode3d-mesh-res-earth-re", str(args.mode3d_mesh_res_earth_re),
            "-mode3d-mesh-res-boundary-re", str(args.mode3d_mesh_res_boundary_re),
            "-mode3d-mesh-coarsening", args.mode3d_mesh_coarsening,
            "-mode3d-mesh-exponent", str(args.mode3d_mesh_exponent),
        ]
        if args.mode3d_parallel_field_init:
            command.append("-mode3d-parallel-field-init")
    return command


def run_process(command: Sequence[str], cwd: Path, log_path: Path) -> int:
    with log_path.open("w") as log:
        log.write("Command:\n  %s\n\n" % " ".join(shlex.quote(value) for value in command))
        log.flush()
        process = subprocess.Popen(
            list(command), cwd=str(cwd), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
            log.flush()
        return process.wait()


def parse_directional_map(path: Path) -> DirectionMap:
    variables: List[str] = []
    x_km = y_km = z_km = float("nan")
    frame = "SM"
    cells: List[DirectionCell] = []
    with path.open(errors="replace") as stream:
        for line_number, raw in enumerate(stream, start=1):
            text = raw.strip()
            if not text:
                continue
            upper = text.upper()
            if upper.startswith("VARIABLES"):
                variables = [value.strip().lower().replace("-", "_").replace(" ", "_")
                             for value in re.findall(r'"([^"]+)"', text)]
                continue
            if upper.startswith("ZONE"):
                for name in ("x_km", "y_km", "z_km"):
                    match = re.search(r"\b%s=([-+0-9.eE]+)" % name, text, re.IGNORECASE)
                    if match:
                        if name == "x_km": x_km = float(match.group(1))
                        elif name == "y_km": y_km = float(match.group(1))
                        else: z_km = float(match.group(1))
                match = re.search(r"\bframe=([^\s\"]+)", text, re.IGNORECASE)
                if match:
                    frame = match.group(1)
                continue
            if upper.startswith("TITLE"):
                continue
            parts = text.split()
            try:
                values = [float(value) for value in parts]
            except ValueError:
                continue
            if not variables:
                if len(values) < 4:
                    continue
                lon, lat, rc, energy = values[:4]
            else:
                index = {name: variables.index(name) for name in variables}
                required = ("lon_deg", "lat_deg", "rc_gv", "emin_mev")
                if not all(name in index for name in required):
                    raise ValueError("%s lacks required directional-map variables" % path)
                if len(values) != len(variables):
                    raise ValueError("%s line %d has %d values for %d variables" %
                                     (path, line_number, len(values), len(variables)))
                lon = values[index["lon_deg"]]
                lat = values[index["lat_deg"]]
                rc = values[index["rc_gv"]]
                energy = values[index["emin_mev"]]
            cells.append(DirectionCell(lon, lat, rc, energy))
    if not cells:
        raise ValueError("no directional-map cells parsed from %s" % path)
    if not all(math.isfinite(value) for value in (x_km, y_km, z_km)):
        raise ValueError("directional map does not identify the observation position: %s" % path)
    return DirectionMap(str(path.resolve()), frame, x_km, y_km, z_km, tuple(cells))


def locate_directional_map(run_dir: Path, solver: str) -> Path:
    exact = (run_dir / "cutoff_gridless_dir_map_point_0000.dat"
             if solver == "GRIDLESS"
             else run_dir / "cutoff_3d_dir_map_loc_000000.dat")
    if exact.exists():
        return exact
    patterns = (["cutoff_gridless_dir_map_point_0000*.dat"] if solver == "GRIDLESS"
                else ["cutoff_3d_dir_map_loc_000000*.dat"])
    matches: List[Path] = []
    for pattern in patterns:
        matches.extend(run_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError("directional-map output is missing: %s" % exact)
    return sorted(matches)[0]


def load_driver_tilts(
        path: Path, required_times: Sequence[datetime] = (),
        ) -> Tuple[List[Tuple[datetime, float]], Dict[str, object]]:
    """Validate the standard AMPS 5-minute T96/T05/TS05 driver and read tilt.

    Expected numerical columns are::

      UTC Bx By Bz Vx Vy Vz Np Temp SYM-H IMFflag SWflag Tilt Pdyn W1..W6

    The runner only needs ``Tilt`` for the SM/GSM aperture transform, while AMPS
    consumes the complete file through ``DRIVER_FILE``.
    """
    rows: List[Tuple[datetime, float]] = []
    header_seen = False
    with path.open() as stream:
        for line_number, raw in enumerate(stream, start=1):
            text = raw.strip()
            if not text:
                continue
            if text.startswith(("#", "!")):
                if "YYYY-MM-DDTHH:MM:SS" in text and "Tilt" in text and "W6" in text:
                    header_seen = True
                continue
            fields = text.split()
            if len(fields) != 20:
                raise ValueError(
                    "driver line %d has %d fields; expected timestamp + 19 values" %
                    (line_number, len(fields)))
            epoch = parse_utc(fields[0])
            try:
                values = [float(value) for value in fields[1:]]
            except ValueError as exc:
                raise ValueError("driver line %d contains a nonnumeric value" % line_number) from exc
            if rows and epoch <= rows[-1][0]:
                raise ValueError("driver timestamps are not strictly increasing at line %d" % line_number)
            rows.append((epoch, values[11]))
    if not rows:
        raise ValueError("driver contains no numerical records: %s" % path)

    gaps = [(second[0] - first[0]).total_seconds()
            for first, second in zip(rows, rows[1:])]
    median_cadence = statistics.median(gaps) if gaps else float("nan")
    maximum_gap = max(gaps) if gaps else float("nan")
    if gaps and not 299.0 <= median_cadence <= 301.0:
        raise ValueError("driver median cadence is not five minutes: %.1f s" % median_cadence)
    if gaps and maximum_gap > 600.0:
        raise ValueError("driver contains a gap larger than ten minutes: %.1f s" % maximum_gap)
    if required_times:
        first_required = min(required_times)
        last_required = max(required_times)
        if first_required < rows[0][0] or last_required > rows[-1][0]:
            raise ValueError(
                "driver coverage %s .. %s does not contain selected C19 epochs %s .. %s" %
                (format_utc(rows[0][0]), format_utc(rows[-1][0]),
                 format_utc(first_required), format_utc(last_required)))

    info: Dict[str, object] = {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_records": len(rows),
        "first_epoch_utc": format_utc(rows[0][0]),
        "last_epoch_utc": format_utc(rows[-1][0]),
        "median_cadence_seconds": median_cadence,
        "maximum_gap_seconds": maximum_gap,
        "standard_header_seen": header_seen,
    }
    return rows, info


def interpolate_tilt(rows: Sequence[Tuple[datetime, float]], epoch: datetime) -> float:
    if epoch <= rows[0][0]:
        return rows[0][1]
    if epoch >= rows[-1][0]:
        return rows[-1][1]
    for first, second in zip(rows, rows[1:]):
        if first[0] <= epoch <= second[0]:
            width = (second[0] - first[0]).total_seconds()
            fraction = (epoch - first[0]).total_seconds() / width if width else 0.0
            return first[1] + fraction * (second[1] - first[1])
    raise RuntimeError("could not interpolate driver tilt")


def norm(vector: Tuple[float, float, float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def unit(vector: Tuple[float, float, float]) -> Tuple[float, float, float]:
    magnitude = norm(vector)
    if magnitude <= 0.0:
        raise ValueError("cannot normalize zero vector")
    return tuple(value / magnitude for value in vector)  # type: ignore[return-value]


def dot(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def cross(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def add_scaled(a: Tuple[float, float, float], b: Tuple[float, float, float], scale: float
               ) -> Tuple[float, float, float]:
    return (a[0] + scale * b[0], a[1] + scale * b[1], a[2] + scale * b[2])


def scale(a: Tuple[float, float, float], factor: float) -> Tuple[float, float, float]:
    return (factor * a[0], factor * a[1], factor * a[2])


def spherical_direction(lon_deg: float, lat_deg: float) -> Tuple[float, float, float]:
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    clat = math.cos(lat)
    return (clat * math.cos(lon), clat * math.sin(lon), math.sin(lat))


def gsm_to_sm(vector: Tuple[float, float, float], tilt_rad: float) -> Tuple[float, float, float]:
    cosine = math.cos(tilt_rad)
    sine = math.sin(tilt_rad)
    x, y, z = vector
    return (cosine * x - sine * z, y, sine * x + cosine * z)


def detector_basis(position_sm: Tuple[float, float, float], direction: str
                   ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]:
    radial = unit(position_sm)
    sm_north = (0.0, 0.0, 1.0)
    east = cross(sm_north, radial)
    if norm(east) < 1.0e-12:
        east = (0.0, 1.0, 0.0)
    east = unit(east)
    boresight = east if direction == "EAST" else scale(east, -1.0)
    # Aperture horizontal axis follows the local equatorial/radial direction;
    # vertical axis follows SM north projected perpendicular to the boresight.
    horizontal = radial
    vertical = unit(cross(boresight, horizontal))
    if dot(vertical, sm_north) < 0.0:
        vertical = scale(vertical, -1.0)
    return boresight, horizontal, vertical


def aperture_coordinates(
        direction: Tuple[float, float, float],
        boresight: Tuple[float, float, float],
        horizontal: Tuple[float, float, float],
        vertical: Tuple[float, float, float],
        ) -> Optional[Tuple[float, float]]:
    forward = dot(direction, boresight)
    if forward <= 0.0:
        return None
    # Tangent-plane angular coordinates preserve the nominal half-angle boundary.
    alpha_h = math.degrees(math.atan2(dot(direction, horizontal), forward))
    alpha_v = math.degrees(math.atan2(dot(direction, vertical), forward))
    return alpha_h, alpha_v


def integrated_power_law(lower: float, upper: float, gamma: float) -> float:
    if upper <= lower or lower <= 0.0:
        return 0.0
    if abs(gamma - 1.0) < 1.0e-12:
        return math.log(upper / lower)
    exponent = 1.0 - gamma
    return (upper ** exponent - lower ** exponent) / exponent


def channel_transmission(cutoff_energy_mev: float, energy_min: float, energy_max: float,
                         gamma: float) -> Optional[float]:
    if not math.isfinite(cutoff_energy_mev) or cutoff_energy_mev < 0.0:
        return None
    denominator = integrated_power_law(energy_min, energy_max, gamma)
    if denominator <= 0.0:
        return None
    lower = max(energy_min, cutoff_energy_mev)
    if lower >= energy_max:
        return 0.0
    return max(0.0, min(1.0, integrated_power_law(lower, energy_max, gamma) / denominator))


def fold_aperture(
        direction_map: DirectionMap,
        position_sm: Tuple[float, float, float],
        detector_direction: str,
        energy_min: float,
        energy_max: float,
        equatorial_half_angle: float,
        north_south_half_angle: float,
        gamma: float,
        ) -> Tuple[Optional[float], int, int, List[Dict[str, object]]]:
    boresight, horizontal, vertical = detector_basis(position_sm, detector_direction)
    weighted_sum = 0.0
    weight_sum = 0.0
    n_cells = 0
    n_unresolved = 0
    diagnostic: List[Dict[str, object]] = []
    for cell in direction_map.cells:
        direction = spherical_direction(cell.lon_deg, cell.lat_deg)
        coordinates = aperture_coordinates(direction, boresight, horizontal, vertical)
        if coordinates is None:
            continue
        alpha_h, alpha_v = coordinates
        ellipse = ((alpha_h / equatorial_half_angle) ** 2
                   + (alpha_v / north_south_half_angle) ** 2)
        if ellipse > 1.0 + 1.0e-12:
            continue
        n_cells += 1
        transmission = channel_transmission(
            cell.cutoff_energy_mev, energy_min, energy_max, gamma)
        if transmission is None:
            n_unresolved += 1
            diagnostic.append({
                "lon_deg": cell.lon_deg, "lat_deg": cell.lat_deg,
                "detector_direction": detector_direction,
                "inside_aperture": True, "transmission": None,
                "cutoff_energy_mev": cell.cutoff_energy_mev,
            })
            continue
        # Spherical cell area is proportional to cos(latitude).  A top-hat
        # detector response then reduces to a solid-angle-weighted average.
        weight = max(0.0, math.cos(math.radians(cell.lat_deg)))
        weighted_sum += weight * transmission
        weight_sum += weight
        diagnostic.append({
            "lon_deg": cell.lon_deg, "lat_deg": cell.lat_deg,
            "detector_direction": detector_direction,
            "inside_aperture": True, "transmission": transmission,
            "cutoff_energy_mev": cell.cutoff_energy_mev,
        })
    value = weighted_sum / weight_sum if weight_sum > 0.0 else None
    return value, n_cells, n_unresolved, diagnostic


def evaluate_reference_row(
        reference: ReferenceRow,
        direction_map: DirectionMap,
        manifest: Mapping[str, object],
        solver: str,
        field_model: str,
        spectral_index: float,
        tilt_rad: float,
        ) -> Tuple[ModelRow, List[Dict[str, object]]]:
    channel = manifest["channels"][reference.channel]
    position_gsm = (direction_map.x_km, direction_map.y_km, direction_map.z_km)
    if direction_map.frame.upper().startswith("SM"):
        position_sm = gsm_to_sm(position_gsm, tilt_rad)
    else:
        # With a GSM fallback map, use the GSM position and axes consistently.
        position_sm = position_gsm

    common = dict(
        direction_map=direction_map,
        position_sm=position_sm,
        energy_min=reference.energy_min_mev,
        energy_max=reference.energy_max_mev,
        equatorial_half_angle=float(channel["equatorial_half_angle_deg"]),
        north_south_half_angle=float(channel["north_south_half_angle_deg"]),
        gamma=spectral_index,
    )
    east, n_east, unresolved_east, east_diag = fold_aperture(
        detector_direction="EAST", **common)
    west, n_west, unresolved_west, west_diag = fold_aperture(
        detector_direction="WEST", **common)
    total_cells = n_east + n_west
    unresolved_fraction = ((unresolved_east + unresolved_west) / float(total_cells)
                           if total_cells else 1.0)
    if east is None or west is None or west <= 0.0 or n_east == 0 or n_west == 0:
        ratio = log_ratio = residual = None
        status = "INSUFFICIENT_APERTURE_COVERAGE"
    else:
        ratio = east / west
        if ratio <= 0.0:
            log_ratio = residual = None
            status = "NONPOSITIVE_MODELED_RATIO"
        else:
            log_ratio = math.log10(ratio)
            residual = log_ratio - reference.log10_east_west_ratio
            status = "VALID"
    row = ModelRow(
        utc=format_utc(reference.utc), spacecraft=reference.spacecraft,
        channel=reference.channel, solver=solver, field_model=field_model,
        observed_east_west_ratio=reference.east_west_ratio,
        observed_log10_east_west_ratio=reference.log10_east_west_ratio,
        modeled_east_west_ratio=ratio,
        modeled_log10_east_west_ratio=log_ratio,
        residual_log10=residual,
        east_transmission=east, west_transmission=west,
        n_east_cells=n_east, n_west_cells=n_west,
        unresolved_direction_fraction=unresolved_fraction,
        spectral_index=spectral_index,
        map_frame=direction_map.frame, map_path=direction_map.path,
        status=status,
    )
    diagnostics = east_diag + west_diag
    for item in diagnostics:
        item.update({
            "utc": format_utc(reference.utc), "spacecraft": reference.spacecraft,
            "channel": reference.channel, "solver": solver, "field_model": field_model,
        })
    return row, diagnostics


def pearson(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    if len(x) < 2 or len(y) != len(x):
        return None
    mean_x = statistics.fmean(x)
    mean_y = statistics.fmean(y)
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    denominator = math.sqrt(sum(value * value for value in dx)
                            * sum(value * value for value in dy))
    if denominator <= 0.0:
        return None
    return sum(a * b for a, b in zip(dx, dy)) / denominator


def calculate_metrics(rows: Sequence[ModelRow], args: argparse.Namespace) -> List[Metrics]:
    groups = sorted({(row.solver, row.field_model, row.channel) for row in rows})
    result: List[Metrics] = []
    for solver, model, channel in groups:
        group = [row for row in rows if (row.solver, row.field_model, row.channel)
                 == (solver, model, channel)]
        valid = [row for row in group if row.status == "VALID"
                 and row.modeled_log10_east_west_ratio is not None]
        observed = [row.observed_log10_east_west_ratio for row in valid]
        modeled = [float(row.modeled_log10_east_west_ratio) for row in valid]
        residuals = [mod - obs for obs, mod in zip(observed, modeled)]
        valid_fraction = len(valid) / float(len(group)) if group else 0.0
        sign_agreement = (
            sum((obs < 0.0) == (mod < 0.0) for obs, mod in zip(observed, modeled))
            / float(len(valid)) if valid else 0.0)
        bias = statistics.fmean(residuals) if residuals else None
        mae = statistics.fmean(abs(value) for value in residuals) if residuals else None
        rmse = math.sqrt(statistics.fmean(value * value for value in residuals)) if residuals else None
        correlation = pearson(observed, modeled)
        passed = (
            valid_fraction + 1.0e-14 >= args.min_valid_fraction
            and sign_agreement + 1.0e-14 >= args.min_sign_agreement
            and mae is not None and mae <= args.max_mae_log10 + 1.0e-14
            and rmse is not None and rmse <= args.max_rmse_log10 + 1.0e-14
            and correlation is not None and correlation + 1.0e-14 >= args.min_correlation
        )
        result.append(Metrics(
            solver, model, channel, len(group), len(valid), valid_fraction,
            sign_agreement, bias, mae, rmse, correlation, passed))
    return result


def write_dict_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    fields: List[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_comparison_plots(rows: Sequence[ModelRow], output_root: Path) -> List[str]:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("C19 plot generation skipped: %s" % exc, file=sys.stderr)
        return []
    outputs: List[str] = []
    if not rows:
        return outputs

    for solver, model in sorted({(row.solver, row.field_model) for row in rows}):
        subset = [row for row in rows if row.solver == solver and row.field_model == model]
        panels = sorted({(row.spacecraft, row.channel) for row in subset})
        fig, axes = plt.subplots(len(panels), 1, figsize=(10.5, max(3.0, 2.5 * len(panels))),
                                 sharex=True, squeeze=False)
        for axis, (spacecraft, channel) in zip(axes[:, 0], panels):
            panel = sorted([row for row in subset
                            if row.spacecraft == spacecraft and row.channel == channel],
                           key=lambda row: row.utc)
            times = [parse_utc(row.utc) for row in panel]
            observed = [row.observed_log10_east_west_ratio for row in panel]
            modeled = [float("nan") if row.modeled_log10_east_west_ratio is None
                       else row.modeled_log10_east_west_ratio for row in panel]
            axis.plot(times, observed, marker="o", markersize=3, linewidth=1.2,
                      label="GOES observed")
            axis.plot(times, modeled, marker="x", markersize=3, linewidth=1.2,
                      label="AMPS modeled")
            axis.axhline(0.0, linewidth=0.8)
            axis.set_ylabel("log10(E/W)")
            axis.set_title("%s %s" % (spacecraft, channel))
            axis.grid(True, alpha=0.3)
            axis.legend(loc="best")
        axes[-1, 0].set_xlabel("UTC")
        fig.suptitle("C19A %s %s: GOES vs AMPS east/west ratio" % (solver, model))
        fig.tight_layout()
        path = output_root / ("C19_comparison_%s_%s.png" % (solver.lower(), model.lower()))
        fig.savefig(path, dpi=160)
        plt.close(fig)
        outputs.append(str(path))

        valid = [row for row in subset if row.status == "VALID"
                 and row.modeled_log10_east_west_ratio is not None]
        if valid:
            fig, ax = plt.subplots(figsize=(6.4, 6.0))
            for channel in sorted({row.channel for row in valid}):
                group = [row for row in valid if row.channel == channel]
                ax.scatter([row.observed_log10_east_west_ratio for row in group],
                           [row.modeled_log10_east_west_ratio for row in group],
                           label=channel, alpha=0.8)
            values = ([row.observed_log10_east_west_ratio for row in valid]
                      + [float(row.modeled_log10_east_west_ratio) for row in valid])
            lower, upper = min(values), max(values)
            margin = max(0.05, 0.05 * (upper - lower if upper > lower else 1.0))
            ax.plot([lower - margin, upper + margin], [lower - margin, upper + margin],
                    linestyle="--", linewidth=1.0)
            ax.set_xlabel("Observed log10(E/W)")
            ax.set_ylabel("Modeled log10(E/W)")
            ax.set_title("C19A %s %s comparison" % (solver, model))
            ax.grid(True, alpha=0.3)
            ax.legend()
            fig.tight_layout()
            path = output_root / ("C19_scatter_%s_%s.png" % (solver.lower(), model.lower()))
            fig.savefig(path, dpi=160)
            plt.close(fig)
            outputs.append(str(path))

        fig, axes = plt.subplots(len(panels), 1, figsize=(10.5, max(3.0, 2.5 * len(panels))),
                                 sharex=True, squeeze=False)
        for axis, (spacecraft, channel) in zip(axes[:, 0], panels):
            panel = sorted([row for row in subset
                            if row.spacecraft == spacecraft and row.channel == channel],
                           key=lambda row: row.utc)
            times = [parse_utc(row.utc) for row in panel]
            east = [float("nan") if row.east_transmission is None else row.east_transmission
                    for row in panel]
            west = [float("nan") if row.west_transmission is None else row.west_transmission
                    for row in panel]
            axis.plot(times, east, linewidth=1.2, label="East aperture")
            axis.plot(times, west, linewidth=1.2, label="West aperture")
            axis.set_ylabel("Transmission")
            axis.set_ylim(-0.03, 1.03)
            axis.set_title("%s %s" % (spacecraft, channel))
            axis.grid(True, alpha=0.3)
            axis.legend(loc="best")
        axes[-1, 0].set_xlabel("UTC")
        fig.suptitle("C19A %s %s modeled broad-aperture transmission" % (solver, model))
        fig.tight_layout()
        path = output_root / ("C19_transmission_%s_%s.png" % (solver.lower(), model.lower()))
        fig.savefig(path, dpi=160)
        plt.close(fig)
        outputs.append(str(path))
    return outputs


def make_aperture_plot(diagnostics: Sequence[Mapping[str, object]], output_path: Path) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None
    if not diagnostics:
        return None
    first_key = tuple(diagnostics[0].get(name) for name in
                      ("utc", "spacecraft", "channel", "solver", "field_model"))
    rows = [row for row in diagnostics if tuple(row.get(name) for name in
            ("utc", "spacecraft", "channel", "solver", "field_model")) == first_key]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    for direction, marker in (("EAST", "o"), ("WEST", "s")):
        group = [row for row in rows if row["detector_direction"] == direction
                 and row.get("transmission") is not None]
        if group:
            scatter = ax.scatter([float(row["lon_deg"]) for row in group],
                                 [float(row["lat_deg"]) for row in group],
                                 c=[float(row["transmission"]) for row in group],
                                 marker=marker, label=direction, vmin=0.0, vmax=1.0)
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_xlabel("SM direction longitude (deg)")
    ax.set_ylabel("SM direction latitude (deg)")
    ax.set_title("C19A aperture sampling: %s %s %s %s %s" % first_key)
    ax.grid(True, alpha=0.3)
    ax.legend()
    if 'scatter' in locals():
        fig.colorbar(scatter, ax=ax, label="Channel transmission")
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return str(output_path)


def render_case_input(
        args: argparse.Namespace, template: Path, run_dir: Path,
        reference: ReferenceRow, solver: str, field_model: str, driver: Path,
        ) -> None:
    replacements = {
        "__RUN_ID__": "C19_%s_%s_%s_%s" % (
            solver.lower(), field_model.lower(), reference.spacecraft.lower(),
            timestamp_token(reference.utc)),
        "__CUTOFF_EMIN_MEV__": "%.12g" % args.cutoff_emin_mev,
        "__CUTOFF_EMAX_MEV__": "%.12g" % args.cutoff_emax_mev,
        "__CUTOFF_SCAN_N__": str(args.cutoff_scan_n),
        "__MAX_TRACE_TIME__": "%.12g" % args.max_trace_time,
        "__DIRMAP_LON_RES__": "%.12g" % args.dir_lon_res_deg,
        "__DIRMAP_LAT_RES__": "%.12g" % args.dir_lat_res_deg,
        "__FIELD_MODEL__": field_model,
        "__EPOCH__": format_utc(reference.utc, suffix_z=False),
        "__DRIVER_FILE__": str(driver.resolve()),
        "__DT_TRACE__": "%.12g" % args.dt_trace,
        "__MAX_TRACE_DISTANCE_RE__": "%.12g" % args.max_trace_distance_re,
        "__MPI_SCHEDULER__": args.scheduler,
        "__DYNAMIC_CHUNK__": str(resolved_dynamic_chunk(args, solver)),
        "__THREADS__": str(args.nt),
        "__MESH_RES_EARTH_RE__": "%.12g" % args.mode3d_mesh_res_earth_re,
        "__MESH_RES_BOUNDARY_RE__": "%.12g" % args.mode3d_mesh_res_boundary_re,
        "__MESH_COARSENING__": args.mode3d_mesh_coarsening,
        "__MESH_EXPONENT__": "%.12g" % args.mode3d_mesh_exponent,
    }
    render_template(template, run_dir / "AMPS_PARAM_C19.in", replacements)
    write_trajectory(run_dir / "C19_trajectory.txt", reference)


def self_test() -> int:
    manifest = json.loads(DEFAULT_MANIFEST.read_text())
    with tempfile.TemporaryDirectory(prefix="C19_runner_selftest_") as temporary:
        root = Path(temporary)
        map_path = root / "cutoff_gridless_dir_map_point_0000.dat"
        lines = [
            'TITLE="synthetic C19 directional map"',
            'VARIABLES="lon_deg","lat_deg","Rc_GV","Emin_MeV"',
            'ZONE T="point=0 x_km=42157 y_km=0 z_km=0" I=24 J=13 F=POINT',
        ]
        for lat in range(-90, 91, 15):
            for lon in range(0, 360, 15):
                # At x>0, physical east is near +SM-y (lon 90). Give east a
                # larger cutoff than west to create E/W < 1.
                delta_e = abs(((lon - 90 + 180) % 360) - 180)
                delta_w = abs(((lon - 270 + 180) % 360) - 180)
                cutoff = 35.0 if delta_e < 45.0 and abs(lat) < 60.0 else 2.0
                if delta_w < 45.0 and abs(lat) < 60.0:
                    cutoff = 2.0
                lines.append("%g %g %g %g" %
                             (lon, lat, rigidity_gv_from_kinetic_energy_mev(cutoff), cutoff))
        map_path.write_text("\n".join(lines) + "\n")
        direction_map = parse_directional_map(map_path)
        reference = ReferenceRow(
            parse_utc("2012-05-17T06:00:00Z"), "GOES13", "P4", 15.0, 40.0,
            0.5, math.log10(0.5), -75.0, 0.0, 35786.0, "SYNTHETIC")
        # Exercise both solver templates and the supplied-driver contract.
        driver_path = root / "synthetic_driver.txt"
        driver_lines = [
            "# YYYY-MM-DDTHH:MM:SS Bx By Bz Vx Vy Vz Np Temp SYM-H IMFflag SWflag Tilt Pdyn W1 W2 W3 W4 W5 W6",
        ]
        for minute, tilt in ((55, -0.1), (0, 0.0), (5, 0.1)):
            hour = 5 if minute == 55 else 6
            epoch = datetime(2012, 5, 17, hour, minute, tzinfo=timezone.utc)
            values = [1.0, 2.0, -3.0, -450.0, 0.0, 0.0, 5.0, 100000.0,
                      -20.0, 1.0, 1.0, tilt, 2.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
            driver_lines.append("%s %s" % (format_utc(epoch, suffix_z=False),
                                            " ".join(str(value) for value in values)))
        driver_path.write_text("\n".join(driver_lines) + "\n")
        _, driver_info = load_driver_tilts(driver_path, [reference.utc])
        if driver_info["n_records"] != 3:
            raise AssertionError("synthetic driver validation failed")

        template_args = argparse.Namespace(
            cutoff_emin_mev=0.5, cutoff_emax_mev=500.0, cutoff_scan_n=20,
            max_trace_time=300.0, dir_lon_res_deg=15.0, dir_lat_res_deg=15.0,
            dt_trace=0.25, max_trace_distance_re=400.0, scheduler="STATIC",
            dynamic_chunk=1, nt=2, mode3d_mesh_res_earth_re=0.1,
            mode3d_mesh_res_boundary_re=2.0, mode3d_mesh_coarsening="LINEAR",
            mode3d_mesh_exponent=1.0)
        for solver, template in (("GRIDLESS", DEFAULT_TEMPLATE_GRIDLESS),
                                 ("GRIDDED", DEFAULT_TEMPLATE_MODE3D)):
            run_dir = root / solver.lower()
            run_dir.mkdir()
            render_case_input(template_args, template, run_dir, reference, solver,
                              "T05", driver_path)
            rendered = (run_dir / "AMPS_PARAM_C19.in").read_text()
            if re.search(r"__[A-Z0-9_]+__", rendered):
                raise AssertionError("%s template retained a placeholder" % solver)

        model, diagnostics = evaluate_reference_row(
            reference, direction_map, manifest, "GRIDLESS", "T05", 3.0, 0.0)
        if model.status != "VALID" or model.modeled_east_west_ratio is None:
            raise AssertionError("synthetic map did not produce a valid model row")
        if not (model.modeled_east_west_ratio < 1.0):
            raise AssertionError("synthetic east cutoff did not produce E/W < 1")
        rows = [model]
        plots = make_comparison_plots(rows, root)
        aperture = make_aperture_plot(diagnostics, root / "C19_aperture_diagnostic.png")
        if not plots or aperture is None:
            raise AssertionError("self-test did not generate plots")
        csv_path = root / "C19_model.csv"
        write_dict_rows(csv_path, [asdict(model)])
        if not csv_path.exists():
            raise AssertionError("self-test did not write model CSV")
    print("C19A runner self-test: PASS")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="C19A GOES EPEAD east-west directional-access validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 srcEarth/test/C19/run_C19.py --self-test
  python3 srcEarth/test/C19/run_C19.py --profile SMOKE --solver GRIDDED \\
    --models T96,T05 --driver /path/to/may2012_driver.txt \\
    --amps ./amps -np 4 -nt 16
  python3 srcEarth/test/C19/run_C19.py --profile ROUTINE --solver BOTH \\
    --mode3d-parallel-field-init --driver /path/to/may2012_driver.txt --amps ./amps
""",
    )
    parser.add_argument("--profile", choices=sorted(PROFILE_STEP_MINUTES), default="ROUTINE")
    parser.add_argument("--time-step-minutes", type=int,
                        help="override profile cadence; 0 keeps every reference epoch")
    parser.add_argument("--start", help="optional inclusive UTC start")
    parser.add_argument("--end", help="optional inclusive UTC end")
    parser.add_argument("--spacecraft", default="GOES13,GOES15")
    parser.add_argument("--channels", default="P4,P5")
    parser.add_argument("--solver", choices=SOLVERS, default="GRIDDED")
    parser.add_argument("--models", default="T96,T05")
    parser.add_argument("--event-manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--reference", default=str(DEFAULT_REFERENCE))
    parser.add_argument("--driver", help="AMPS-format five-minute T96/T05/TS05 event driver")
    parser.add_argument("--spectral-index", type=float, default=3.0,
                        help="incident J(E) proportional to E^-gamma; default gamma=3")
    parser.add_argument("--dir-lon-res-deg", type=float, default=10.0)
    parser.add_argument("--dir-lat-res-deg", type=float, default=10.0)
    parser.add_argument("--cutoff-emin-mev", type=float, default=0.5)
    parser.add_argument("--cutoff-emax-mev", type=float, default=500.0)
    parser.add_argument("--cutoff-scan-n", type=int, default=120)
    parser.add_argument("--dt-trace", type=float, default=0.25)
    parser.add_argument("--max-trace-time", type=float, default=300.0)
    parser.add_argument("--max-trace-distance-re", type=float, default=400.0)
    parser.add_argument("--mover", default="BORIS")
    parser.add_argument("--scheduler", choices=("STATIC", "BLOCK_CYCLIC", "DYNAMIC"),
                        default="DYNAMIC")
    parser.add_argument("--dynamic-chunk", type=int, default=0)
    parser.add_argument("--mode3d-mesh-res-earth-re", type=float, default=0.05)
    parser.add_argument("--mode3d-mesh-res-boundary-re", type=float, default=2.0)
    parser.add_argument("--mode3d-mesh-coarsening", choices=("LINEAR", "LOG", "POWER"),
                        default="LINEAR")
    parser.add_argument("--mode3d-mesh-exponent", type=float, default=1.0)
    parser.add_argument("--mode3d-parallel-field-init", action="store_true")
    parser.add_argument("-np", type=int, default=4)
    parser.add_argument("-nt", type=int, default=16)
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--amps", default="./amps")
    parser.add_argument("--output-root", default="test_output/C19_goes_epead_ew")
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--enforce-acceptance", action="store_true",
                        help="return exit status 1 when provisional observational gates fail")
    parser.add_argument("--min-valid-fraction", type=float, default=0.85)
    parser.add_argument("--min-sign-agreement", type=float, default=0.90)
    parser.add_argument("--min-correlation", type=float, default=0.60)
    parser.add_argument("--max-mae-log10", type=float, default=0.20)
    parser.add_argument("--max-rmse-log10", type=float, default=0.30)
    args = parser.parse_args(argv)
    try:
        args.spacecraft_list = parse_csv_list(args.spacecraft, ("GOES13", "GOES15"))
        args.channel_list = parse_csv_list(args.channels, ("P4", "P5"))
        args.model_list = parse_csv_list(args.models, FIELD_MODELS)
    except ValueError as exc:
        parser.error(str(exc))
    if args.np < 1 or args.nt < 1:
        parser.error("-np and -nt must be >= 1")
    if args.dynamic_chunk < 0:
        parser.error("--dynamic-chunk must be >= 0")
    if args.time_step_minutes is not None and args.time_step_minutes < 0:
        parser.error("--time-step-minutes must be >= 0")
    if args.dir_lon_res_deg <= 0.0 or args.dir_lat_res_deg <= 0.0:
        parser.error("directional-map resolutions must be positive")
    if args.cutoff_scan_n < 2:
        parser.error("--cutoff-scan-n must be >= 2")
    if args.spectral_index <= 0.0:
        parser.error("--spectral-index must be positive")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.self_test:
        return self_test()

    launch_dir = Path.cwd().resolve()
    reference_path = Path(args.reference).expanduser().resolve()
    manifest_path = Path(args.event_manifest).expanduser().resolve()
    driver_path = Path(args.driver).expanduser().resolve() if args.driver else None
    output_root = Path(args.output_root).expanduser()
    if not output_root.is_absolute():
        output_root = (launch_dir / output_root).resolve()
    amps = Path(args.amps).expanduser()
    if not amps.is_absolute():
        amps = (launch_dir / amps).resolve()

    if not reference_path.exists():
        print("C19A reference is missing: %s" % reference_path, file=sys.stderr)
        print("Generate it with build_goes_reference.py --download.", file=sys.stderr)
        return 2
    if not manifest_path.exists():
        print("C19A event manifest is missing: %s" % manifest_path, file=sys.stderr)
        return 2
    if driver_path is None or not driver_path.exists():
        print("C19A requires --driver with a May-2012 AMPS T96/T05/TS05 driver.",
              file=sys.stderr)
        return 2
    if not args.skip_run and not args.dry_run and not amps.exists():
        print("AMPS executable is missing: %s" % amps, file=sys.stderr)
        return 2

    try:
        reference_all = load_reference(reference_path)
        reference = select_reference_rows(reference_all, args)
        driver_tilts, driver_info = load_driver_tilts(
            driver_path, [row.utc for row in reference])
    except Exception as exc:
        print("C19A input validation failed: %s" % exc, file=sys.stderr)
        return 2
    manifest = json.loads(manifest_path.read_text())
    grouped = group_reference(reference)
    solvers = ("GRIDLESS", "GRIDDED") if args.solver == "BOTH" else (args.solver,)

    if not args.skip_run:
        if output_root.exists() and not args.keep:
            shutil.rmtree(output_root)
        output_root.mkdir(parents=True, exist_ok=True)
    elif not output_root.exists():
        print("--skip-run requested but output root does not exist: %s" % output_root,
              file=sys.stderr)
        return 2

    commands: List[Dict[str, object]] = []
    model_rows: List[ModelRow] = []
    aperture_diagnostics: List[Dict[str, object]] = []
    run_failures: List[Dict[str, object]] = []

    print("C19A selected %d reference rows at %d spacecraft epochs" %
          (len(reference), len(grouped)))
    print("C19A launches: %d" % (len(grouped) * len(solvers) * len(args.model_list)))

    for solver in solvers:
        template = DEFAULT_TEMPLATE_GRIDLESS if solver == "GRIDLESS" else DEFAULT_TEMPLATE_MODE3D
        for field_model in args.model_list:
            for (epoch, spacecraft), reference_group in grouped.items():
                representative = reference_group[0]
                run_dir = (output_root / solver.lower() / field_model.lower()
                           / spacecraft.lower() / timestamp_token(epoch))
                if not args.skip_run:
                    run_dir.mkdir(parents=True, exist_ok=True)
                    render_case_input(args, template, run_dir, representative,
                                      solver, field_model, driver_path)
                command = command_for(args, amps, solver)
                command_record = {
                    "solver": solver, "field_model": field_model,
                    "spacecraft": spacecraft, "utc": format_utc(epoch),
                    "cwd": str(run_dir), "command": command,
                }
                commands.append(command_record)
                if args.dry_run:
                    print("[%s %s %s %s]\n  %s" %
                          (solver, field_model, spacecraft, format_utc(epoch),
                           " ".join(shlex.quote(value) for value in command)))
                    continue
                if not args.skip_run:
                    return_code = run_process(command, run_dir, run_dir / "C19_amps.log")
                    if return_code != 0:
                        run_failures.append(dict(command_record, return_code=return_code))
                        continue
                try:
                    map_path = locate_directional_map(run_dir, solver)
                    direction_map = parse_directional_map(map_path)
                    tilt = interpolate_tilt(driver_tilts, epoch)
                    # Official TS05 files and AMPS use radians for the Tilt column.
                    for reference_row in reference_group:
                        model, diagnostics = evaluate_reference_row(
                            reference_row, direction_map, manifest, solver,
                            field_model, args.spectral_index, tilt)
                        model_rows.append(model)
                        if not aperture_diagnostics:
                            aperture_diagnostics.extend(diagnostics)
                except Exception as exc:
                    run_failures.append(dict(command_record, return_code=None,
                                             analysis_error=str(exc)))

    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "C19_commands.json").write_text(
        json.dumps(commands, indent=2, sort_keys=True) + "\n")
    if args.dry_run:
        print("C19A dry run complete; inputs and commands were generated in %s" % output_root)
        return 0

    write_dict_rows(output_root / "C19_reference_used.csv", [
        {
            "utc": format_utc(row.utc), "spacecraft": row.spacecraft,
            "channel": row.channel, "energy_min_mev": row.energy_min_mev,
            "energy_max_mev": row.energy_max_mev,
            "east_west_ratio": row.east_west_ratio,
            "log10_east_west_ratio": row.log10_east_west_ratio,
            "longitude_deg_east": row.longitude_deg_east,
            "latitude_deg": row.latitude_deg, "altitude_km": row.altitude_km,
            "position_source": row.position_source,
        } for row in reference])
    write_dict_rows(output_root / "C19_model.csv", [asdict(row) for row in model_rows])
    write_dict_rows(output_root / "C19_comparison.csv", [asdict(row) for row in model_rows])
    write_dict_rows(output_root / "C19_aperture_samples.csv", aperture_diagnostics)
    metrics = calculate_metrics(model_rows, args)
    write_dict_rows(output_root / "C19_metrics.csv", [asdict(row) for row in metrics])
    plot_paths = make_comparison_plots(model_rows, output_root)
    aperture_plot = make_aperture_plot(
        aperture_diagnostics, output_root / "C19_aperture_diagnostic.png")
    if aperture_plot:
        plot_paths.append(aperture_plot)

    numerical_complete = not run_failures and bool(model_rows)
    observational_passed = bool(metrics) and all(row.passed for row in metrics)
    result = {
        "test_id": "C19A",
        "test_name": "GOES EPEAD east-west directional-access validation",
        "profile": args.profile,
        "solver": args.solver,
        "field_models": args.model_list,
        "spacecraft": args.spacecraft_list,
        "channels": args.channel_list,
        "reference_path": str(reference_path),
        "reference_sha256": sha256(reference_path),
        "event_manifest_path": str(manifest_path),
        "event_manifest_sha256": sha256(manifest_path),
        "driver_path": str(driver_path),
        "driver_sha256": sha256(driver_path),
        "driver_validation": driver_info,
        "spectral_index": args.spectral_index,
        "instrument_response": "uniform elliptical top-hat inside nominal P4/P5 FOV",
        "observable": "log10(background-subtracted physical EAST/WEST flux ratio)",
        "n_reference_rows": len(reference),
        "n_model_rows": len(model_rows),
        "n_run_failures": len(run_failures),
        "run_failures": run_failures,
        "metrics": [asdict(row) for row in metrics],
        "acceptance_thresholds_provisional": {
            "min_valid_fraction": args.min_valid_fraction,
            "min_sign_agreement": args.min_sign_agreement,
            "min_correlation": args.min_correlation,
            "max_mae_log10": args.max_mae_log10,
            "max_rmse_log10": args.max_rmse_log10,
        },
        "numerical_complete": numerical_complete,
        "observational_passed": observational_passed,
        "acceptance_enforced": args.enforce_acceptance,
        "passed": numerical_complete and (observational_passed or not args.enforce_acceptance),
        "plot_files": plot_paths,
        "limitations": [
            "C19A uses nominal broad top-hat P4/P5 apertures, not a complete energy-angle response matrix.",
            "The incident spectrum is a common isotropic power law; prompt-onset interplanetary anisotropy is not modeled.",
            "The telemetry-head-to-physical-direction mapping is event-specific and fixed in the manifest.",
            "Publication runs should supply NOAA one-minute ephemeris rather than nominal GEO-slot positions.",
        ],
    }
    (output_root / "C19_result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n")

    summary_lines = [
        "C19A GOES EPEAD east-west directional-access validation",
        "reference rows: %d" % len(reference),
        "model rows: %d" % len(model_rows),
        "run failures: %d" % len(run_failures),
    ]
    for metric in metrics:
        summary_lines.append(
            "%s %s %s: valid=%.3f sign=%.3f MAE=%s RMSE=%s corr=%s -> %s" % (
                metric.solver, metric.field_model, metric.channel,
                metric.valid_fraction, metric.sign_agreement_fraction,
                "NA" if metric.mean_absolute_error_log10 is None else "%.4f" % metric.mean_absolute_error_log10,
                "NA" if metric.rmse_log10 is None else "%.4f" % metric.rmse_log10,
                "NA" if metric.correlation is None else "%.4f" % metric.correlation,
                "PASS" if metric.passed else "FAIL"))
    summary_lines.append("overall: %s" % ("PASS" if result["passed"] else "FAIL"))
    (output_root / "C19_summary.txt").write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))

    if not numerical_complete:
        return 2
    if args.enforce_acceptance and not observational_passed:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

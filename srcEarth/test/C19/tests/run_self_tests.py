#!/usr/bin/env python3
"""Run the C19A package checks without observational downloads or AMPS."""

from __future__ import annotations

import csv
import gzip
import math
import re
import py_compile
import subprocess
import sys
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(command):
    print("+", " ".join(str(value) for value in command), flush=True)
    completed = subprocess.run(command, cwd=str(ROOT), check=False)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


def write_reference(path: Path) -> None:
    fields = [
        "utc", "spacecraft", "channel", "energy_min_mev", "energy_max_mev",
        "telemetry_head_east", "telemetry_head_west",
        "east_west_ratio", "log10_east_west_ratio",
        "east_flux_background_subtracted", "west_flux_background_subtracted",
        "flux_product_policy", "east_flux_variable", "west_flux_variable",
        "east_flux_correction_state", "west_flux_correction_state",
        "longitude_deg_east", "latitude_deg", "altitude_km", "position_source",
    ]
    rows = []
    for hour in (6, 7):
        utc = "2012-05-17T%02d:00:00Z" % hour
        for spacecraft, longitude in (("GOES13", -75.0), ("GOES15", -135.0)):
            for channel, lo, hi, ratio in (("P4", 15.0, 40.0, 0.50),
                                           ("P5", 38.0, 82.0, 0.70)):
                rows.append({
                    "utc": utc,
                    "spacecraft": spacecraft,
                    "channel": channel,
                    "energy_min_mev": lo,
                    "energy_max_mev": hi,
                    # Mimic the May-2012 yaw mapping: GOES13 physical-east stream is
                    # telemetry W, whereas GOES15 physical-east stream is telemetry E.
                    # This proves that attitude lookup is by actual detector ID rather
                    # than by the compatibility words EAST/WEST.
                    "telemetry_head_east": ("W" if spacecraft == "GOES13" else "E"),
                    "telemetry_head_west": ("E" if spacecraft == "GOES13" else "W"),
                    "east_west_ratio": ratio,
                    "log10_east_west_ratio": math.log10(ratio),
                    # P1.3 needs two positive physical-WEST channel intensities per
                    # epoch.  Values are synthetic but decrease with energy so the
                    # default OBSERVED_WEST fit is well-defined.
                    "east_flux_background_subtracted": (50.0 if channel == "P4" else 10.0) * ratio,
                    "west_flux_background_subtracted": 50.0 if channel == "P4" else 10.0,
                    "flux_product_policy": "UNCORRECTED",
                    "east_flux_variable": "%s_SYNTH_UNCOR_FLUX" % channel,
                    "west_flux_variable": "%s_SYNTH_UNCOR_FLUX" % channel,
                    "east_flux_correction_state": "UNCORRECTED",
                    "west_flux_correction_state": "UNCORRECTED",
                    "longitude_deg_east": longitude,
                    "latitude_deg": 0.0,
                    "altitude_km": 35786.0,
                    "position_source": "SYNTHETIC",
                })
    with gzip.open(path, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_driver(path: Path) -> None:
    start = datetime(2012, 5, 17, 5, 55, tzinfo=timezone.utc)
    lines = [
        "# YYYY-MM-DDTHH:MM:SS Bx By Bz Vx Vy Vz Np Temp SYM-H "
        "IMFflag SWflag Tilt Pdyn W1 W2 W3 W4 W5 W6"
    ]
    for index in range(15):
        epoch = start + timedelta(minutes=5 * index)
        values = [1.0, 2.0, -3.0, -450.0, 0.0, 0.0, 5.0, 100000.0,
                  -20.0, 1.0, 1.0, 0.01 * index, 2.0,
                  0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        lines.append("%s %s" % (
            epoch.strftime("%Y-%m-%dT%H:%M:%S"),
            " ".join(str(value) for value in values)))
    path.write_text("\n".join(lines) + "\n")


def write_orientation(path: Path) -> None:
    """Write non-antipodal per-head attitude records for vector-aperture validation."""
    fields = [
        "utc", "spacecraft", "detector", "frame",
        "boresight_x", "boresight_y", "boresight_z",
        "aperture_north_x", "aperture_north_y", "aperture_north_z", "source",
    ]
    rows = []
    for hour in (6, 7):
        for spacecraft in ("GOES13", "GOES15"):
            # Deliberately make raw telemetry head W non-antipodal to head E to
            # prove that neither the runner nor AMPS derives one detector direction
            # by negating the other. The reference decides which head is numerator.
            for detector, bore in (("E", (0.0, 1.0, 0.0)),
                                   ("W", (0.15, -0.97, 0.18))):
                rows.append({
                    "utc": "2012-05-17T%02d:00:00Z" % hour,
                    "spacecraft": spacecraft, "detector": detector, "frame": "SM",
                    "boresight_x": bore[0], "boresight_y": bore[1], "boresight_z": bore[2],
                    "aperture_north_x": 0.0, "aperture_north_y": 0.0,
                    "aperture_north_z": 1.0, "source": "SYNTHETIC_C19_TEST",
                })
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def integration_dry_run() -> None:
    with tempfile.TemporaryDirectory(prefix="C19A_integration_") as temporary:
        root = Path(temporary)
        reference = root / "reference.csv.gz"
        driver = root / "driver.txt"
        output = root / "output"
        orientation = root / "orientation.csv"
        write_reference(reference)
        write_driver(driver)
        write_orientation(orientation)
        run([
            sys.executable, str(ROOT / "run_C19.py"),
            "--profile", "FULL",
            "--solver", "BOTH",
            "--models", "T05",
            "--reference", str(reference),
            "--driver", str(driver),
            "--output-root", str(output),
            "--amps", "./amps-not-required-for-dry-run",
            "-np", "2", "-nt", "2",
            "--dry-run",
        ])
        rendered = list(output.glob("**/AMPS_PARAM_C19.in"))
        trajectories = list(output.glob("**/C19_trajectory.txt"))
        if len(rendered) != 8 or len(trajectories) != 8:
            raise SystemExit(
                "C19A integration dry-run wrote %d inputs and %d trajectories; expected 8" %
                (len(rendered), len(trajectories)))
        aperture_files = list(output.glob("**/C19_directional_apertures.dat"))
        if len(aperture_files) != 8:
            raise SystemExit("C19A integration dry-run wrote %d aperture files; expected 8" %
                             len(aperture_files))
        for aperture_file in aperture_files:
            text = aperture_file.read_text()
            if "EAST LOCAL_SM" not in text or "WEST LOCAL_SM" not in text:
                raise SystemExit("SM_PROXY did not render generic LOCAL_SM aperture vectors: %s" %
                                 aperture_file)
        for path in rendered:
            text = path.read_text()
            if re.search(r"__[A-Z0-9_]+__", text):
                raise SystemExit("unresolved input placeholder in %s" % path)
            if "DRIVER_FILE" not in text or str(driver.resolve()) not in text:
                raise SystemExit("generated input does not contain absolute driver path: %s" % path)
            # P0.1 production contract: a normal C19 run must use the structured
            # three-state penumbra path and preserve trace-limit outcomes as
            # UNRESOLVED.  This integration assertion catches a regression in the
            # Python template renderer even when AMPS itself is not available.
            if not re.search(r"^CUTOFF_SEARCH_ALGORITHM\s+PENUMBRA_SCAN\s*$",
                             text, re.MULTILINE):
                raise SystemExit("normal C19 dry-run did not select PENUMBRA_SCAN: %s" % path)
            if not re.search(r"^CUTOFF_TRACE_LIMIT_POLICY\s+UNRESOLVED\s*$",
                             text, re.MULTILINE):
                raise SystemExit("normal C19 dry-run did not select UNRESOLVED policy: %s" % path)
            if "gridded" in str(path):
                if not re.search(r"^CUTOFF_RIGIDITY_LIST_GV\s+[0-9]", text, re.MULTILINE):
                    raise SystemExit("P1 GRIDDED input did not contain direct-access rigidity list: %s" % path)
            else:
                if re.search(r"^CUTOFF_RIGIDITY_LIST_GV\s+", text, re.MULTILINE):
                    raise SystemExit("GRIDLESS proxy input unexpectedly retained P1 rigidity list: %s" % path)

        # Single-workflow regression.  P0/P1/P2 are implementation history, not
        # mutually exclusive runner modes.  The ordinary dry run above must therefore
        # use the current production settings directly.
        import json
        commands_path = output / "C19_commands.json"
        if not commands_path.exists():
            raise SystemExit("single-workflow dry-run did not write C19_commands.json")
        commands = json.loads(commands_path.read_text())
        if not commands:
            raise SystemExit("single-workflow dry-run wrote an empty command list")
        for record in commands:
            command = record["command"]
            joined = " ".join(str(value) for value in command)
            if "-cutoff-search PENUMBRA_SCAN" not in joined:
                raise SystemExit("current workflow did not force PENUMBRA_SCAN")
            if "-cutoff-dirmap-coverage VECTOR_APERTURES" not in joined:
                raise SystemExit("current workflow command did not select VECTOR_APERTURES coverage")
            if "-cutoff-dirmap-aperture-file C19_directional_apertures.dat" not in joined:
                raise SystemExit("current workflow command lacks generic aperture-vector file")
            if record["solver"] == "GRIDDED":
                if "-mode3d-mesh-res-earth-re 0.025" not in joined:
                    raise SystemExit("current GRIDDED workflow did not use 0.025 Re near-Earth mesh")
                if "-mode3d-mesh-res-boundary-re 1.0" not in joined:
                    raise SystemExit("current GRIDDED workflow did not use 1.0 Re boundary mesh")
                if int(record.get("n_direct_access_rigidities", 0)) <= 0:
                    raise SystemExit("current GRIDDED workflow did not request direct A(E,Omega)")

        # The generated parameter decks must use the finest P2.1 production sky grid.
        for path in rendered:
            text = path.read_text()
            if not re.search(r"^DIRMAP_LON_RES\s+2.5\s*$", text, re.MULTILINE):
                raise SystemExit("current workflow did not use 2.5-degree longitude grid: %s" % path)
            if not re.search(r"^DIRMAP_LAT_RES\s+2.5\s*$", text, re.MULTILINE):
                raise SystemExit("current workflow did not use 2.5-degree latitude grid: %s" % path)
            if not re.search(r"^DIRMAP_COVERAGE\s+VECTOR_APERTURES\s*$", text, re.MULTILINE):
                raise SystemExit("current workflow did not select generic vector-aperture coverage: %s" % path)
            if not re.search(r"^DIRMAP_APERTURE_FILE\s+C19_directional_apertures\.dat\s*$", text, re.MULTILINE):
                raise SystemExit("current workflow did not reference its aperture-vector file: %s" % path)

        # Confirm that the public CLI no longer exposes historical alternate-runner
        # selectors.  Advanced scientific inputs remain configurable, but there is
        # only one execution path.
        help_text = subprocess.check_output(
            [sys.executable, str(ROOT / "run_C19.py"), "--help"],
            cwd=str(ROOT), text=True)
        for obsolete in ("--p0-diagnostic", "--p2-diagnostic",
                         "--cutoff-search", "--trace-limit-policy", "--response-fold"):
            if obsolete in help_text:
                raise SystemExit("obsolete alternate-mode option still exposed: %s" % obsolete)
        if "--direction-coverage" not in help_text:
            raise SystemExit("directional coverage selector is missing from runner CLI")
        if "--max-discrete-transition-fraction" not in help_text:
            raise SystemExit("finite-rigidity-grid uncertainty control is missing from runner CLI")
        if "epead_response_C19_uncorrected_extended.csv" not in help_text:
            # argparse normally exposes the default only when help includes it, so test
            # the source constant as a second layer below if formatting changes.
            runner_text = (ROOT / "run_C19.py").read_text()
            if 'DEFAULT_RESPONSE = SCRIPT_DIR / "data" / "epead_response_C19_uncorrected_extended.csv"' not in runner_text:
                raise SystemExit("runner no longer defaults to the extended uncorrected response")

        # Exact attitude regression: the preferred orientation schema supplies one
        # arbitrary vector per *actual telemetry head*. Head W is deliberately non-
        # antipodal to E. The reference maps physical numerator/denominator streams to
        # those IDs, and the generated AMPS aperture file must preserve both vectors.
        attitude_output = root / "integration_attitude_vectors"
        run([
            sys.executable, str(ROOT / "run_C19.py"),
            "--profile", "SMOKE", "--solver", "GRIDDED", "--models", "T05",
            "--reference", str(reference), "--driver", str(driver),
            "--output-root", str(attitude_output),
            "--amps", "./amps-not-required-for-dry-run", "-np", "1", "-nt", "1",
            "--detector-orientation-source", "FILE",
            "--detector-orientation-file", str(orientation), "--dry-run",
        ])
        vector_files = list(attitude_output.glob("**/C19_directional_apertures.dat"))
        if not vector_files:
            raise SystemExit("FILE attitude dry-run did not write a generic aperture file")
        vector_text = vector_files[0].read_text()
        if "E SM 0 1 0" not in vector_text:
            raise SystemExit("FILE attitude telemetry-head E vector was not propagated")
        west_lines = [line for line in vector_text.splitlines() if line.startswith("W SM ")]
        if len(west_lines) != 1:
            raise SystemExit("FILE attitude telemetry-head W vector was not propagated")
        parts = west_lines[0].split()
        west_bore = tuple(float(value) for value in parts[2:5])
        # Head E is exactly (0,1,0); a derived antipode would be (0,-1,0).
        if all(abs(a-b) < 1.0e-12 for a, b in zip(west_bore, (0.0, -1.0, 0.0))):
            raise SystemExit("FILE attitude head W vector was incorrectly derived as -head E")

        # FULL_SPHERE remains an explicit alternative to the optimized default.
        full_output = root / "integration_full_sphere"
        run([
            sys.executable, str(ROOT / "run_C19.py"),
            "--profile", "SMOKE", "--solver", "GRIDDED", "--models", "T05",
            "--reference", str(reference), "--driver", str(driver),
            "--output-root", str(full_output),
            "--amps", "./amps-not-required-for-dry-run", "-np", "1", "-nt", "1",
            "--direction-coverage", "FULL_SPHERE", "--dry-run",
        ])
        full_inputs = list(full_output.glob("**/AMPS_PARAM_C19.in"))
        if not full_inputs:
            raise SystemExit("FULL_SPHERE dry-run did not render an AMPS input")
        if not all(re.search(r"^DIRMAP_COVERAGE\s+FULL_SPHERE\s*$",
                             path.read_text(), re.MULTILINE) for path in full_inputs):
            raise SystemExit("FULL_SPHERE runner option was not propagated into input decks")


def validate_committed_inputs() -> None:
    required_directives = {
        "CALC_TARGET", "FIELD_EVAL_METHOD", "CUTOFF_EMIN", "CUTOFF_EMAX",
        "CUTOFF_NENERGY", "CUTOFF_SEARCH_ALGORITHM",
        "CUTOFF_TRACE_LIMIT_POLICY", "CUTOFF_UPPER_SCAN_N", "CUTOFF_RIGIDITY_LIST_GV",
        "CUTOFF_MAX_TRAJ_TIME", "DIRECTIONAL_MAP",
        "DIRMAP_LON_RES", "DIRMAP_LAT_RES", "DIRMAP_COVERAGE",
        "DIRMAP_APERTURE_FILE",
        "SPECIES", "CHARGE",
        "MASS_AMU", "FIELD_MODEL", "EPOCH", "DRIVER_FILE",
        "DOMAIN_X_MIN", "DOMAIN_X_MAX", "DOMAIN_Y_MIN", "DOMAIN_Y_MAX",
        "DOMAIN_Z_MIN", "DOMAIN_Z_MAX", "R_INNER", "SPECTRUM_TYPE",
        "SPEC_GAMMA", "OUTPUT_MODE", "TRAJ_FRAME", "TRAJ_FILE",
        "OUTPUT_COORDS", "DT_TRACE", "MAX_STEPS", "MAX_TRACE_TIME",
        "MAX_TRACE_DISTANCE", "TRAP_DETECTION",
    }
    for name in ("AMPS_PARAM_C19_gridless.in", "AMPS_PARAM_C19_mode3d.in"):
        path = ROOT / name
        text = path.read_text()
        if re.search(r"__[A-Z0-9_]+__", text):
            raise SystemExit("committed input contains a macro placeholder: %s" % path)
        directives = set()
        for raw in text.splitlines():
            stripped = raw.strip()
            if stripped and not stripped.startswith(("!", "#")):
                directives.add(stripped.split(None, 1)[0].upper())
        missing = sorted(required_directives - directives)
        if missing:
            raise SystemExit("%s lacks explicit directive(s): %s" %
                             (path, ", ".join(missing)))
    trajectory = ROOT / "C19_trajectory.txt"
    if not trajectory.exists() or not trajectory.read_text().strip():
        raise SystemExit("committed default C19_trajectory.txt is missing or empty")

    # The default C19 observational reference is UNCORRECTED, so its synthetic detector
    # model must not silently stop at the nominal P5 82-MeV boundary.  Protect the
    # committed response contract here without importing run_C19.py: both documented
    # high-energy secondary components and the 190-MeV P5 upper support must remain.
    response_path = ROOT / "data" / "epead_response_C19_uncorrected_extended.csv"
    if not response_path.exists():
        raise SystemExit("extended uncorrected EPEAD response is missing: %s" % response_path)
    with response_path.open(newline="") as stream:
        response_rows = list(csv.DictReader(stream))
    secondary = [row for row in response_rows
                 if row.get("response_component", "").upper().startswith("SECONDARY")]
    if not secondary:
        raise SystemExit("extended EPEAD response contains no secondary proton components")
    p5_hi = max(float(row["energy_max_mev"]) for row in response_rows
                if row.get("channel", "").upper() == "P5"
                and float(row.get("relative_response", 0.0)) > 0.0)
    if not math.isclose(p5_hi, 190.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise SystemExit("extended P5 response no longer reaches 190 MeV: %.12g" % p5_hi)



def validate_directional_coverage_source_contract() -> None:
    """Static regression checks for the AMPS-side aperture selector.

    The package self-tests run without building AMPS, so they cannot execute the C++
    trajectory scheduler.  These checks still protect the cross-layer contract that is
    easiest to break during future refactors: parser fields, generic cutoff CLI override,
    and both Mode3D/gridless compact full-grid-id mappings must remain present together.
    """
    src_earth = ROOT.parents[1]
    checks = {
        src_earth / "util" / "amps_param_parser.cpp": (
            "DIRMAP_COVERAGE", "DIRMAP_APERTURE_FILE",
            "DIRMAP_APERTURE", "VECTOR_APERTURES"),
        src_earth / "util" / "cutoff_cli.cpp": (
            "-cutoff-dirmap-coverage",
            "-cutoff-dirmap-aperture-file",
            "-cutoff-dirmap-aperture"),
        src_earth / "3d" / "CutoffRigidityMode3D.cpp": (
            "ApplyDirectionalMapCoverage3D", "VECTOR_APERTURES",
            "LOCAL_SM", "fullGridCellIds"),
        src_earth / "gridless" / "CutoffRigidityGridless.cpp": (
            "VECTOR_APERTURES", "LOCAL_SM", "dirMapFullCellIds",
            "DIRMAP coverage"),
    }
    for path, needles in checks.items():
        text = path.read_text(errors="replace")
        for needle in needles:
            if needle not in text:
                raise SystemExit("directional-coverage source contract missing %r in %s" %
                                 (needle, path))

    # Geometry sanity check for the documented C19 default.  This mirrors the
    # local SM_PROXY aperture test at a representative equatorial GEO position.
    # The exact retained count varies slightly with spacecraft longitude/grid phase,
    # but it must be far smaller than the 10,512-cell 2.5-degree full sphere while
    # remaining of order 1.7k cells.
    lon_res = lat_res = 2.5
    n_lon = round(360.0 / lon_res)
    n_lat = round(180.0 / lat_res) + 1
    full = n_lon * n_lat
    radial = (1.0, 0.0, 0.0)
    east = (0.0, 1.0, 0.0)
    horizontal = radial
    vertical = (0.0, 0.0, 1.0)
    selected = 0
    for j in range(n_lat):
        lat = -90.0 + lat_res * j
        clat = math.cos(math.radians(lat))
        for i in range(n_lon):
            lon = lon_res * i
            arrival = (clat * math.cos(math.radians(lon)),
                       clat * math.sin(math.radians(lon)),
                       math.sin(math.radians(lat)))
            look = tuple(-value for value in arrival)
            inside = False
            for sign in (1.0, -1.0):
                boresight = tuple(sign * value for value in east)
                forward = sum(a*b for a, b in zip(look, boresight))
                if forward <= 0.0:
                    continue
                ah = math.degrees(math.atan2(
                    sum(a*b for a, b in zip(look, horizontal)), forward))
                av = math.degrees(math.atan2(
                    sum(a*b for a, b in zip(look, vertical)), forward))
                if (ah/30.0)**2 + (av/60.0)**2 <= 1.0 + 1.0e-12:
                    inside = True
                    break
            if inside:
                selected += 1
    if full != 10512 or not (1600 <= selected <= 1900):
        raise SystemExit("unexpected generic local-vector aperture geometry: selected=%d full=%d" %
                         (selected, full))


def main() -> int:
    validate_committed_inputs()
    validate_directional_coverage_source_contract()
    scripts = (ROOT / "run_C19.py", ROOT / "build_goes_reference.py")
    for script in scripts:
        print("Compiling", script.name, flush=True)
        py_compile.compile(str(script), doraise=True)

    run([sys.executable, str(ROOT / "build_goes_reference.py"), "--self-test"])
    run([sys.executable, str(ROOT / "run_C19.py"), "--self-test"])
    integration_dry_run()
    run([sys.executable, str(ROOT / "build_goes_reference.py"), "--help"])
    run([sys.executable, str(ROOT / "run_C19.py"), "--help"])

    print("C19A package self-tests: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

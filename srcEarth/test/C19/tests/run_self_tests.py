#!/usr/bin/env python3
"""Run the C19A package checks without observational downloads or AMPS."""

from __future__ import annotations

import csv
import gzip
import math
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
        "east_west_ratio", "log10_east_west_ratio", "longitude_deg_east",
        "latitude_deg", "altitude_km", "position_source",
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
                    "east_west_ratio": ratio,
                    "log10_east_west_ratio": math.log10(ratio),
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


def integration_dry_run() -> None:
    with tempfile.TemporaryDirectory(prefix="C19A_integration_") as temporary:
        root = Path(temporary)
        reference = root / "reference.csv.gz"
        driver = root / "driver.txt"
        output = root / "output"
        write_reference(reference)
        write_driver(driver)
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
        for path in rendered:
            if "__" in path.read_text():
                raise SystemExit("unresolved input placeholder in %s" % path)


def main() -> int:
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

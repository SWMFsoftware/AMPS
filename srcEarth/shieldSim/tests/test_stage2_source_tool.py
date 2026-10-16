#!/usr/bin/env python3
"""End-to-end checks for the Stage-2 published-source command-line boundary."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def run(command: list[str], expect_success: bool = True) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if expect_success and completed.returncode != 0:
        raise AssertionError(f"command failed ({completed.returncode}): {' '.join(command)}\n{completed.stderr}")
    if not expect_success and completed.returncode == 0:
        raise AssertionError(f"command unexpectedly succeeded: {' '.join(command)}")
    return completed


def main() -> int:
    if len(sys.argv) != 3:
        raise AssertionError("usage: test_stage2_source_tool.py TOOL SOURCE_ROOT")
    tool = pathlib.Path(sys.argv[1])
    root = pathlib.Path(sys.argv[2])
    fixtures = root / "tests" / "reference_spectra"

    with tempfile.TemporaryDirectory(prefix="shieldsim-stage2-") as temporary:
        directory = pathlib.Path(temporary)
        ams_a = directory / "ams-a.json"
        ams_b = directory / "ams-b.json"
        ams_command = [
            str(tool),
            "ams02",
            f"--info={fixtures / 'gcr' / 'ams02_20131118_hapi_info.json'}",
            f"--data={fixtures / 'gcr' / 'ams02_20131118_hapi_data.csv'}",
            f"--interpretation={fixtures / 'gcr' / 'ams02_20131118_interpretation.v1.json'}",
        ]
        run([*ams_command, f"--output={ams_a}"])
        run([*ams_command, f"--output={ams_b}"])
        require(ams_a.read_bytes() == ams_b.read_bytes(), "AMS canonical output is not deterministic")
        ams = json.loads(ams_a.read_text(encoding="utf-8"))
        require(ams["schema_version"] == "shieldSim.spectrum_dataset.v2", "AMS output schema")
        require([item["species_id"] for item in ams["components"]] ==
                ["proton", "helium_mixed_isotopes"], "AMS species order/identity")

        goes_a = directory / "goes-a.json"
        goes_b = directory / "goes-b.json"
        goes_command = [
            str(tool),
            "goes16",
            f"--csv={fixtures / 'sep' / 'goes16_sgps_gle72_timeseries.csv'}",
            f"--definition={fixtures / 'sep' / 'goes16_gle72_event_definition.v1.json'}",
        ]
        run([*goes_command, f"--output={goes_a}"])
        run([*goes_command, f"--output={goes_b}"])
        require(goes_a.read_bytes() == goes_b.read_bytes(), "GOES canonical output is not deterministic")
        goes = json.loads(goes_a.read_text(encoding="utf-8"))
        require(goes["schema_version"] == "shieldSim.goes_event_fluence.v1", "GOES output schema")
        require(len(goes["windows"]) == 2, "GOES dual-window output")
        require(all(len(window["channels"]) == 26 for window in goes["windows"]),
                "GOES directions/channels were lost")
        require(any(not channel["complete"] for window in goes["windows"]
                    for channel in window["channels"]), "GOES quality gaps were hidden")

        failed = run([*ams_command, "--output=-", "--typo=ignored"], expect_success=False)
        require("unknown option" in failed.stderr, "unknown source-tool option was not diagnosed")

        bad_ams_contract = directory / "bad-ams-interpretation.json"
        bad_ams = json.loads(
            (fixtures / "gcr" / "ams02_20131118_interpretation.v1.json").read_text(encoding="utf-8")
        )
        bad_ams["data_sha256"] = "sha256:" + "0" * 64
        bad_ams_contract.write_text(json.dumps(bad_ams), encoding="utf-8")
        failed = run(
            [
                str(tool),
                "ams02",
                f"--info={fixtures / 'gcr' / 'ams02_20131118_hapi_info.json'}",
                f"--data={fixtures / 'gcr' / 'ams02_20131118_hapi_data.csv'}",
                f"--interpretation={bad_ams_contract}",
                "--output=-",
            ],
            expect_success=False,
        )
        require("SHA-256" in failed.stderr, "AMS source-hash mutation was not rejected")

        bad_goes_definition = directory / "bad-goes-definition.json"
        bad_goes = json.loads(
            (fixtures / "sep" / "goes16_gle72_event_definition.v1.json").read_text(encoding="utf-8")
        )
        bad_goes["timestamp_semantics"]["source_coordinate"] = "midpoint"
        bad_goes_definition.write_text(json.dumps(bad_goes), encoding="utf-8")
        failed = run(
            [
                str(tool),
                "goes16",
                f"--csv={fixtures / 'sep' / 'goes16_sgps_gle72_timeseries.csv'}",
                f"--definition={bad_goes_definition}",
                "--output=-",
            ],
            expect_success=False,
        )
        require("timestamp semantics" in failed.stderr,
                "GOES timestamp-contract mutation was not rejected")

    print("PASS Stage-2 source tool determinism, schemas, directions, and mutation rejection")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        sys.exit(1)

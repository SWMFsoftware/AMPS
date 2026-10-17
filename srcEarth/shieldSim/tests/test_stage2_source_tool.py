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

        cycle_a = directory / "ams-cycle-a.json"
        cycle_b = directory / "ams-cycle-b.json"
        cycle_command = [
            str(tool),
            "ams02-cycle",
            f"--interpretation={fixtures / 'gcr' / 'ams02_2025' / 'interpretation.v1.json'}",
            "--bartels-rotation=2581",
        ]
        run([*cycle_command, f"--output={cycle_a}"])
        run([*cycle_command, f"--output={cycle_b}"])
        require(cycle_a.read_bytes() == cycle_b.read_bytes(),
                "AMS solar-cycle canonical output is not deterministic")
        cycle = json.loads(cycle_a.read_text(encoding="utf-8"))
        require(cycle["bartels_rotation"] == 2581 and len(cycle["components"]) == 7,
                "AMS solar-cycle selection/species coverage")
        require(all(len(item["value"]) == 40
                    and len(item["statistical_uncertainty"]) == 40
                    and len(item["time_dependent_systematic_uncertainty"]) == 40
                    and len(item["time_independent_systematic_uncertainty"]) == 40
                    for item in cycle["components"]),
                "AMS solar-cycle uncertainty components were lost")
        require(all(item["mass_number"] == 0 and item["geant4_particle"] == ""
                    for item in cycle["components"]),
                "AMS inclusive elements acquired an assumed isotope")

        cris_a = directory / "cris-a.json"
        cris_b = directory / "cris-b.json"
        cris_command = [
            str(tool),
            "ace-cris",
            f"--info={fixtures / 'gcr' / 'ace_cris' / 'ac_h3_cris_hapi_info.json'}",
            f"--data={fixtures / 'gcr' / 'ace_cris' / 'ac_h3_cris_20191201_all.csv'}",
            f"--interpretation={fixtures / 'gcr' / 'ace_cris' / 'interpretation.v1.json'}",
        ]
        run([*cris_command, f"--output={cris_a}"])
        run([*cris_command, f"--output={cris_b}"])
        require(cris_a.read_bytes() == cris_b.read_bytes(),
                "ACE/CRIS canonical output is not deterministic")
        cris = json.loads(cris_a.read_text(encoding="utf-8"))
        require(len(cris["components"]) == 24
                and cris["components"][0]["atomic_number"] == 5
                and cris["components"][-1]["atomic_number"] == 28,
                "ACE/CRIS B--Ni coverage")
        require(sum(not valid for component in cris["components"] for valid in component["valid"]) == 19,
                "ACE/CRIS zero-count validity masks were lost")
        require(all("statistical_uncertainty" in component
                    and "time_dependent_systematic_uncertainty" not in component
                    for component in cris["components"]),
                "ACE/CRIS uncertainty boundary changed")

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

        bad_cris_contract = directory / "bad-cris-interpretation.json"
        bad_cris = json.loads(
            (fixtures / "gcr" / "ace_cris" / "interpretation.v1.json").read_text(encoding="utf-8")
        )
        bad_cris["hapi_data_sha256"] = "0" * 64
        bad_cris_contract.write_text(json.dumps(bad_cris), encoding="utf-8")
        failed = run(
            [
                str(tool),
                "ace-cris",
                f"--info={fixtures / 'gcr' / 'ace_cris' / 'ac_h3_cris_hapi_info.json'}",
                f"--data={fixtures / 'gcr' / 'ace_cris' / 'ac_h3_cris_20191201_all.csv'}",
                f"--interpretation={bad_cris_contract}",
                "--output=-",
            ],
            expect_success=False,
        )
        require("SHA-256" in failed.stderr, "ACE/CRIS source-hash mutation was not rejected")

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

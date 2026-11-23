#!/usr/bin/env python3
"""Exercise the named .in startup path with a linked srcSEP executable."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


SEP_ROOT = Path(__file__).resolve().parents[1]
DECK = SEP_ROOT.parent / "input/sep_stage13.in"
CASES = ["eclipse-2020-12-14", "gle-2012-05-17", "sep-2013-04-11",
         "sep-2001-09-24", "sep-1978-09-23", "sep-1979-03-01",
         "ulysses-2003-oct-nov"]


def expect(command, directory, code, snippet):
    result = subprocess.run(command, cwd=directory, text=True, check=False,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=30)
    if result.returncode != code or snippet not in result.stdout:
        raise RuntimeError(f"{command!r} returned {result.returncode}, expected {code}:\n{result.stdout}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    args = parser.parse_args()
    amps = str(args.amps.resolve())
    with tempfile.TemporaryDirectory(prefix="srcsep-stage13-cli-") as temporary:
        work = Path(temporary)
        command = [amps, "-i", str(DECK)]
        found = expect([*command, "--list-validation-cases"], work, 0, CASES[0]).splitlines()
        if found != CASES:
            raise RuntimeError(f"unexpected case order or inventory: {found}")
        expect(command, work, 1, "requires an explicit -case NAME")
        expect([*command, "-case", "missing"], work, 2, "available: " + ", ".join(CASES))
        expect([amps, "--case", CASES[0]], work, 1, "use -i FILE")
        expect([*command, "-case", CASES[0], "--particle-mover", "parker"], work, 1,
               "cannot override a validation case")
        report = expect([amps, f"--input={DECK}", f"--case={CASES[-1]}",
                         "--validate-runtime-input"], work, 0, "stage13.runtime_input=PASS")
        if "stage13.physical_outer_radius_m=747989353500\n" not in report or \
           "stage13.control_outer_radius_m=1196782965600\n" not in report:
            raise RuntimeError(f"the Ulysses control changed the physical domain:\n{report}")
        expect([*command, "-case", CASES[0]], work, 77, "stage13.comparison=SKIP")
        invalid = work / "bad.in"
        invalid.write_text(DECK.read_text(encoding="utf-8").replace(
            "domain.outer_radius_au = 5", "domain.outer_radius_au = 8"), encoding="utf-8")
        expect([amps, "-i", str(invalid), "-case", CASES[0], "--validate-runtime-input"],
               work, 2, "outer particle boundary must be 5 AU")
        # A completed integrity campaign proves each native overlay resolves
        # and the reviewed observation bytes match. Its schema explicitly
        # records that no spherical model prediction has been scored.
        audit = work / "integrity"
        runner = SEP_ROOT / "validation/run_stage13.py"
        expect([sys.executable, "-B", str(runner), "--amps", amps, "--all",
                "--validate-runtime-input", "--output-dir", str(audit)], work, 0,
               f"PASS {CASES[-1]}")
        report = json.loads((audit / "stage13-validation.json").read_text(encoding="utf-8"))
        if report["mode"] != "input-integrity" or report["totals"]["passed"] != len(CASES):
            raise RuntimeError(f"incorrect integrity report: {report['totals']}")
        if any(case["spherical_observation_validation"] != "NOT_RUN" for case in report["results"]):
            raise RuntimeError("integrity mode claimed that a spherical simulation was scored")

        # Preserve a red path for corrupted local data through the same
        # runtime deck and report layer used by the real campaign.
        case_path = SEP_ROOT / "validation/cases/stage13" / CASES[0] / "case.json"
        corrupt_package = json.loads(case_path.read_text(encoding="utf-8"))
        corrupt_package["artifacts"] = [{"path": "altered.txt", "bytes": 1, "sha256": "0" * 64}]
        (work / "altered.txt").write_text("x", encoding="utf-8")
        manifest = work / "altered.json"
        manifest.write_text(json.dumps(corrupt_package), encoding="utf-8")
        override = work / "altered.in"
        override.write_text(DECK.read_text(encoding="utf-8").replace(
            f"../srcSEP/validation/cases/stage13/{CASES[0]}/case.json", str(manifest)), encoding="utf-8")
        expect([sys.executable, "-B", str(runner), "--amps", amps, "-i", str(override),
                "-case", CASES[0], "--validate-runtime-input", "--output-dir", str(work / "negative")],
               work, 2, "observation artifact differs from its recorded bytes")
        if (work / "PT").exists():
            raise RuntimeError("startup/listing allocated a simulation output directory")
    print("PASS Stage 13 named input selection, physical/control boundaries, and observation integrity")


if __name__ == "__main__":
    main()

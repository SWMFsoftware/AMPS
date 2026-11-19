#!/usr/bin/env python3
"""Acquire/extract real isotope and constant evidence without changing physics.

Source records are immutable within this packet. The narrowly defined parser
requires exact native field names and unique isotope/constant records. The
result is a reference-mass inventory, not a physical state/kinetic policy.
No numerical values are inferred for unavailable molecular/atomic states.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
from html import unescape
import json
from pathlib import Path
import re
import sys
from urllib.error import URLError
from urllib.request import Request, urlopen


SOURCES = {
    "isotopes.html": "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all",
    "codata2022.txt": "https://physics.nist.gov/cuu/Constants/Table/allascii.txt",
    "hydrogen-ground.html": "https://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=1",
    "sodium-ground.html": "https://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=11",
}
ISOTOPES = (("H", 1), ("He", 4), ("O", 16), ("Ne", 20), ("Na", 23), ("Ar", 40))
CONSTANTS = ("atomic mass constant", "electron mass", "elementary charge", "speed of light in vacuum")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plain_html(value):
    """Preserve source field text while removing HTML presentation only."""
    return unescape(re.sub(r"<[^>]+>", "", value))


def atomic_ground_record(raw, filename, element):
    """Extract a native labelled energy row, preserving theory parentheses.

    Elemental NIST index records supply ground-state context and eV units.
    They are a threshold source, not an isotope-shift or uncertainty table;
    missing uncertainty remains null and is not inferred from decimal digits.
    """
    text = plain_html((raw / filename).read_text())
    match = re.search(r"Ionization Energy:\s*\(?([0-9]+\.[0-9]+)\)?\s*eV", text)
    configuration = re.search(r"Ground-state Configuration:\s*([^\n]+)", text)
    level = re.search(r"Ground-state Level:\s*([^\n]+)", text)
    if not match or not configuration or not level:
        raise ValueError("native atomic ground-state fields unavailable: " + filename)
    return {"element": element, "ionization_energy_eV": match.group(1),
        "native_energy_record": match.group(0),
        "ground_configuration_native_text": configuration.group(1).strip(),
        "ground_level_native_text": level.group(1).strip(),
        "uncertainty_eV": None, "isotope_shift_mapping_qualified": False,
        "source": filename, "source_sha256": digest(raw / filename)}


def extract(raw):
    """Return uniquely identified source records; never infer missing fields."""
    text = plain_html((raw / "isotopes.html").read_text())
    isotopes = []
    for symbol, number in ISOTOPES:
        matches = []
        for block in text.split("Atomic Number =")[1:]:
            fields = dict(re.findall(r"^\s*([^=\n]+?)\s*=\s*([^\n]*)", block, re.M))
            if fields.get("Atomic Symbol", "").strip() == symbol and fields.get("Mass Number", "").strip() == str(number):
                native = fields.get("Relative Atomic Mass", "").strip()
                match = re.fullmatch(r"(\d+(?:\.\d+)?)\((\d+)\)", native)
                if not match:
                    raise ValueError("unrecognized or estimated isotope mass: " + symbol + str(number))
                matches.append({"identity": str(number) + symbol, "element": symbol,
                    "mass_number": number, "neutral_relative_mass_u": match.group(1),
                    "native_mass_with_uncertainty": native,
                    "source": "isotopes.html", "source_sha256": digest(raw / "isotopes.html")})
        if len(matches) != 1:
            raise ValueError("expected one source isotope record: " + symbol + str(number))
        isotopes.extend(matches)
    table = (raw / "codata2022.txt").read_text()
    if "2022 CODATA adjustment" not in table:
        raise ValueError("the requested CODATA2022 source changed version")
    constants = {}
    for name in CONSTANTS:
        rows = [line for line in table.splitlines() if line[:60].strip() == name]
        if len(rows) != 1:
            raise ValueError("expected unique native CODATA constant: " + name)
        # NIST fixed-width columns preserve the distinction between numerical
        # uncertainty and exact definitions. Do not strip a unit into a value.
        row = rows[0]
        value = row[60:85].strip()
        uncertainty = row[85:110].strip()
        unit = row[110:].strip()
        number = Decimal(value.replace(" ", ""))
        if not number.is_finite() or number <= 0:
            raise ValueError("invalid native positive constant: " + name)
        constants[name] = {"value": str(number), "native_value": value,
            "native_uncertainty": uncertainty, "unit": unit,
            "source": "codata2022.txt", "source_sha256": digest(raw / "codata2022.txt")}
    if constants["atomic mass constant"]["unit"] != "kg" or constants["electron mass"]["unit"] != "kg":
        raise ValueError("native mass-constant units changed")
    if constants["elementary charge"]["unit"] != "C" or constants["speed of light in vacuum"]["unit"] != "m s^-1":
        raise ValueError("native electrical/light units changed")
    mass_unit = Decimal(constants["atomic mass constant"]["value"])
    for item in isotopes:
        item["reference_neutral_mass_kg"] = str(Decimal(item["neutral_relative_mass_u"]) * mass_unit)
    ground = [atomic_ground_record(raw, "hydrogen-ground.html", "H"),
              atomic_ground_record(raw, "sodium-ground.html", "Na")]
    return {"schema": "srcMoon-isotope-reference-mass-evidence-v1", "isotopes": isotopes,
        "atomic_ground_thresholds": ground,
        "constants": constants, "production_activation": False,
        "scope": "acquired reference-mass evidence; no complete formation/state/kinetic policy",
        "physical_state_transport_ready": False,
        "remaining_requirements": ["isotope/state formation and excitation crosswalk",
            "channel-specific parent/product threshold interpretation", "independently reviewed error bounds"]}


def run(output, raw_directory=None, fetch=False):
    """Preserve all attempts and an honest extraction result in a fresh packet."""
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "raw" if fetch else Path(raw_directory or output / "raw").resolve()
    if fetch:
        raw.mkdir()
    manifest = {"schema": "srcMoon-reference-mass-source-manifest-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(), "sources": [],
        "processing": "native-field extraction and documented isotope-u to kg conversion only",
        "canonical_data_modified": False}
    result = {"status": "ERROR", "science_tests_executed": 0,
        "formal_M3_status": "SKIPPED / NOT QUALIFIED"}
    try:
        for filename, url in SOURCES.items():
            record = {"url": url, "path": str(raw / filename)}
            manifest["sources"].append(record)
            if fetch:
                try:
                    request = Request(url, headers={"User-Agent": "AMPS-source-provenance/1.0"})
                    with urlopen(request, timeout=45) as response:
                        payload = response.read()
                        record["http_status"] = response.status
                    (raw / filename).write_bytes(payload)
                except (OSError, URLError) as error:
                    record["acquisition_error"] = str(error)
                    raise FileNotFoundError("public source acquisition unavailable: " + url) from error
            if not (raw / filename).is_file():
                raise FileNotFoundError("required native source missing: " + str(raw / filename))
            record.update(sha256=digest(raw / filename), byte_count=(raw / filename).stat().st_size)
        value = extract(raw)
        (output / "reference_mass_inventory.json").write_text(json.dumps(value, indent=2) + "\n")
        result.update(status="PASS", scope="source extraction only, not science qualification",
            isotope_records=len(value["isotopes"]), constant_records=len(value["constants"]),
            rejected_or_filtered_science_records=0)
    except FileNotFoundError as error:
        result.update(status="SKIPPED", reason=str(error))
    except (OSError, ValueError, ArithmeticError) as error:
        result.update(status="ERROR", reason=str(error))
    finally:
        manifest["script_sha256"] = digest(__file__)
        (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    if args.fetch and args.raw_dir:
        parser.error("--fetch and --raw-dir are mutually exclusive")
    try:
        result = run(args.output_dir, args.raw_dir, args.fetch)
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return {"PASS": 0, "ERROR": 2, "SKIPPED": 77}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())

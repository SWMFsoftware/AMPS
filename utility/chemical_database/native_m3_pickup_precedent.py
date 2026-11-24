#!/usr/bin/env python3
"""Build a source-access packet from the Poppe et al. lunar pickup ion deposit."""

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


RECORD_URL = "https://zenodo.org/api/records/6994648"
DATA_DOI = "10.5281/zenodo.6994648"
ARTICLE_DOI = "10.1029/2022JE007422"
NATIVE_NAME = "exospheric_ionization_rates.txt"
URLS = {
    "zenodo.json": RECORD_URL,
    NATIVE_NAME: RECORD_URL + "/files/" + NATIVE_NAME + "/content",
    "poppe.pdf": "https://research.ssl.berkeley.edu/~poppe/papers/Poppe_et_al_2022.pdf",
}
NATIVE_SHA = "a208a94f318bbe2c63c59e6e28cfd15bf78901b7f7560480c8298bb235b7826e"
ARTICLE_SHA = "2c8facf7d53366bb65c9f4285ed0e18e335536527926f9c898a1718747ea1ec0"
DATASET = "poppe2022-lunar-pickup-rates"
NATIVE_ID = "poppe2022:ionization-rates"
TABLE = "poppe2022:atomic-1au-rates"
PROCESSES = ("photoionization", "solar_wind_proton_charge_exchange", "solar_wind_electron_impact")
LEVELS = ("minimum", "median", "maximum")
SPECIES = ("He", "Ne", "Na", "Ar36", "Ar")
COLUMNS = ("native_file_line", "native_species") + tuple(
    f"{process}_{level}_per_s" for process in PROCESSES for level in LEVELS)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(directory):
    """Keep each successful HTTP response and its actual retrieval identity together."""
    directory.mkdir(parents=True, exist_ok=False)
    receipts = {}
    for name, url in URLS.items():
        path = directory / name
        status, effective_url = subprocess.check_output(
            ["curl", "--fail", "--silent", "--show-error", "--location", "--proto", "=https",
             "--proto-redir", "=https", "--max-time", "90", "--output", str(path),
             "--write-out", "%{response_code}\n%{url_effective}", url], text=True).splitlines()
        receipts[name] = dict(http_status=int(status), source_uri=url, final_url=effective_url,
                              sha256=digest(path), retrieved_utc=datetime.now(timezone.utc).isoformat())
    (directory / "retrieval.json").write_text(json.dumps(receipts, indent=2) + "\n")


def atomic_rows(path):
    """Project only the reviewed native atomic labels, preserving row and token identity."""
    if digest(path) != NATIVE_SHA:
        raise ValueError("the Zenodo 6994648 ionization rate file has changed")
    lines = path.read_text(encoding="ascii").splitlines()
    if ("s^(-1)" not in lines[1] or "minimum, median, and maximum" not in lines[2]
            or "at 1 au" not in lines[2] or lines[5].split() !=
            ["Species", "Photoionization", "SW", "Charge", "Exchange", "Electron", "Impact"]):
        raise ValueError("unexpected units, conditions, or ordering in the native header")
    rows = []
    for number, line in enumerate(lines[7:], 8):
        name, *tokens = line.split()
        if name in SPECIES:
            if len(tokens) != 9:
                raise ValueError(f"expected three rate triplets at native line {number}")
            for start in (0, 3, 6):
                low, mid, high = map(Decimal, tokens[start:start + 3])
                if not 0 < low <= mid <= high:
                    raise ValueError(f"invalid frequency ordering at native line {number}")
            rows.append((number, name, *tokens))
    if tuple(row[1] for row in rows) != SPECIES:
        raise ValueError("missing, duplicated, or reordered selected species")
    return rows


def write_table(native, output):
    """Rebuild the published CSV using only the licensed native numerical parent."""
    with output.open("x", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(atomic_rows(native))


def prepare(acquisition, output):
    receipts = json.loads((acquisition / "retrieval.json").read_text())
    for name, url in URLS.items():
        receipt = receipts[name]
        if (receipt["http_status"] != 200 or receipt["source_uri"] != url or
                receipt["final_url"] != url or receipt["sha256"] != digest(acquisition / name)):
            raise ValueError(f"missing or mismatched native retrieval for {name}")

    # Zenodo's API can refresh record statistics. Bind the grant to this version
    # DOI, the creator list, and the exact deposited member, then save the hash
    # of the actual rights record used. A licence from a sibling deposit or the
    # numerical table alone cannot establish the required permission.
    record = json.loads((acquisition / "zenodo.json").read_text())
    authors = [item["name"] for item in record["metadata"]["creators"]]
    member, = [item for item in record["files"] if item["key"] == NATIVE_NAME]
    native = acquisition / NATIVE_NAME
    if (record["id"] != 6994648 or record["doi"] != DATA_DOI or
            record["metadata"]["license"]["id"] != "cc-by-4.0" or
            authors != ["Poppe, Andrew", "Halekas, Jasper", "Harada, Yuki"] or
            member["links"]["self"] != URLS[NATIVE_NAME] or
            member["size"] != native.stat().st_size or
            member["checksum"] != "md5:149afa8c4f7ee6d8d276f097ffc59ec1" or
            digest(native) != NATIVE_SHA or digest(acquisition / "poppe.pdf") != ARTICLE_SHA):
        raise ValueError("deposit, grant, numerical member, or cited publication identity changed")
    output.mkdir(parents=True, exist_ok=False)
    for name in URLS:
        shutil.copyfile(acquisition / name, output / name)
    recipe = output / Path(__file__).name
    shutil.copyfile(Path(__file__), recipe)
    rights = [dict(path="zenodo.json", sha256=receipts["zenodo.json"]["sha256"], source_uri=RECORD_URL)]
    metadata = dict(data_doi=DATA_DOI, associated_article_doi=ARTICLE_DOI,
                    attribution="A. R. Poppe, J. S. Halekas and Y. Harada (2022), Zenodo 6994648.",
                    licence="CC BY 4.0", licence_url="https://creativecommons.org/licenses/by/4.0/",
                    unit="s^-1 per neutral", heliocentric_distance_au=1,
                    quantity="Model ionization frequencies for three solar and solar-wind conditions.",
                    change_to_native="None. The acquired numerical file is byte-identical to the deposit.",
                    figure_2_ranges="Photoionization varies with solar activity. Charge exchange and electron impact span the authors' solar-wind proton-flux and electron-density ranges.",
                    measurement_uncertainty=None, daughter_state_probabilities=None,
                    limitations="The endpoints describe condition variability. The source supplies no state-resolved cross sections or event-level product kinematics.")
    source = dict(receipts[NATIVE_NAME], id=NATIVE_ID, path=NATIVE_NAME, logical_name=NATIVE_NAME,
                  role="source-native", policy="zenodo-cc-by-4.0", rights_evidence=rights,
                  quantity_metadata=metadata, physical_activation=False)
    # The publication explicitly reserves copyright on its first page.
    # Its physics description is cited separately and contributes no source
    # bytes to the numerical projection. Its PDF stays in private acquisition.
    article = dict(receipts["poppe.pdf"], id="poppe2022:article", path="poppe.pdf",
                   logical_name="poppe.pdf", role="source-metadata", policy="metadata-only",
                   rights_evidence=[], physical_activation=False,
                   quantity_metadata=dict(doi=ARTICLE_DOI, locator="Table 1, Figure 2 and section 2.2",
                       photoionization_source="Literature compilation, primarily Huebner and Mukherjee (2015).",
                       kinetics="Singly charged test ions in uniform, steady perpendicular solar-wind velocity and IMF. Initial parent-neutral motion and lunar gravity are neglected.",
                       observable="High-energy pickup flux compared to ARTEMIS observations without mass composition.",
                       crosswalk="The published rows are He, Ne, Na, 36Ar and 40Ar. The deposit records the two argon rows as Ar36 and Ar, respectively.",
                       scope="No electronic-state or post-collision energy/angle probabilities."))
    csv_path = output / "poppe2022-atomic-1au-rates.csv"
    write_table(native, csv_path)
    table = dict(id=TABLE + ":csv", logical_name=csv_path.name, path=csv_path.name, sha256=digest(csv_path),
                 source_uri=URLS[NATIVE_NAME], role="derived-table", policy=source["policy"],
                 rights_evidence=rights, derived_from=[NATIVE_ID],
                 recipe=dict(path=recipe.name, sha256=digest(recipe)), physical_activation=False,
                 quantity_metadata=dict(metadata, change_to_native="Five atomic rows, native numerical strings, one-based source line and named columns. All conditions and individual frequency units retained."),
                 table=dict(id=TABLE, name=csv_path.stem, format="csv", columns=COLUMNS,
                            row_count=len(SPECIES), schema_id="poppe-atomic-rates-v1"))
    packet = dict(schema="amps-chemical-source-packet-v1", id=DATASET, version="20261124",
                  title="Lunar pickup ion model frequencies at 1 au", provider=metadata["attribution"],
                  description="Original model frequencies and a selected native-row projection.",
                  limitations=[metadata["limitations"]], physical_activation=False,
                  sources=[source, table, article])
    manifest = output / "source_packet.json"
    manifest.write_text(json.dumps(packet, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("fetch", help="retrieve and receipt three reviewed originals").add_argument("directory", type=Path)
    build = commands.add_parser("prepare", help="prepare a packet from retained originals offline")
    build.add_argument("acquisition", type=Path)
    build.add_argument("output", type=Path)
    projection = commands.add_parser("project", help="replay the exact numerical CSV from the licensed native file")
    projection.add_argument("native", type=Path)
    projection.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "fetch":
        acquire(args.directory)
    elif args.command == "prepare":
        print(prepare(args.acquisition, args.output))
    else:
        write_table(args.native, args.output)

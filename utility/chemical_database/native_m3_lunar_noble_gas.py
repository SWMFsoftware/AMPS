#!/usr/bin/env python3
"""Reproduce Peschel et al. (2026) Table 1 neutral losses from the source PDF."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


DOI = "10.3847/PSJ/ae4901"
URL = "https://iopscience.iop.org/article/" + DOI + "/pdf"
PDF = "peschel-et-al-2026.pdf"
ARTICLE_SHA = "2aa3fc35bf93027e5879c904250c840e0a3db2ee37fd4330f5bdd5feffd978ae"
PACKET = "peschel2026-lunar-noble-gas-loss"
SOURCE = "peschel2026:article"
TABLE = "peschel2026:quiet-sun-neutral-loss"
POLICY = "aas-peschel-2026-cc-by-4.0"
COLUMNS = ("species", "rate_coefficient_printed", "power_of_ten_printed",
           "photoionization_frequency_per_s")
SPECIES = ("He", "Ne", "40Ar")
D04_SHA = "3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(directory):
    """Retain the original PDF with the time and final URL of this transfer."""
    directory.mkdir(parents=True, exist_ok=False)
    path = directory / PDF
    status, final_url = subprocess.check_output(
        ["curl", "-fsSL", "--proto", "=https", "--proto-redir", "=https", "--max-time", "60",
         "-o", str(path), "-w", "%{http_code}\n%{url_effective}", URL], text=True).splitlines()
    if status != "200" or digest(path) != ARTICLE_SHA:
        raise ValueError("the response is not the reviewed Peschel et al. original")
    receipt = dict(source_uri=URL, http_status=int(status), final_url=final_url,
                   retrieved_utc=datetime.now(timezone.utc).isoformat(), sha256=ARTICLE_SHA)
    (directory / "retrieval.json").write_text(json.dumps(receipt, indent=2) + "\n")


def project(pdf, destination):
    """Preserve the precision and signed powers printed in the neutral loss row."""
    if digest(pdf) != ARTICLE_SHA:
        raise ValueError("unreviewed article edition")
    page = subprocess.check_output(
        ["pdftotext", "-layout", "-f", "3", "-l", "3", str(pdf), "-"], encoding="utf-8")
    # The printed isotope superscript precedes its symbol in the PDF text
    # layer. Check the complete header and table footnote before binding the
    # row to He, Ne and radiogenic Ar-40, rather than arbitrary element labels.
    header = r"(?m)^\s*40\s*\nSpecies\s+Helium\s+Neon\s+Ar\s*$"
    if (not re.search(header, page) or "Photoionization rates are adopted from quiet-Sun" not in page or
            "Illuminated volume of the exosphere" not in page):
        raise ValueError("Table 1 species or source conditions unavailable")
    lines = [line for line in page.splitlines() if line.startswith("Photoionization rate (s−1)")]
    if len(lines) != 1:
        raise ValueError("expected a unique species-resolved loss row")
    fields = re.findall(r"(\d+\.\d+) × 10(−\d+)", lines[0])
    if len(fields) != len(SPECIES):
        raise ValueError("incomplete neutral loss row")
    with destination.open("x", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(COLUMNS)
        for species, (coefficient, power) in zip(SPECIES, fields):
            exponent = power.replace("−", "-")
            writer.writerow([species, coefficient, exponent, coefficient + "e" + exponent])


def prepare(acquisition, output):
    """Attach physical scope, provenance and source-local licence to the projection."""
    source = acquisition / PDF
    receipt = json.loads((acquisition / "retrieval.json").read_text())
    if (receipt["source_uri"] != URL or receipt["http_status"] != 200 or
            receipt["sha256"] != ARTICLE_SHA or not receipt["final_url"].startswith("https://") or
            digest(source) != ARTICLE_SHA):
        raise ValueError("altered article or invalid retrieval")
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(source, output / PDF)
    recipe = output / Path(__file__).name
    shutil.copyfile(__file__, recipe)
    csv_path = output / "peschel2026-table1-photoionization.csv"
    project(source, csv_path)
    metadata = dict(doi=DOI, citation="A. Peschel, C. Grava, J. Breternitz and P. Reiss (2026), PSJ 7, 69.",
                    licence="CC BY 4.0", licence_url="https://creativecommons.org/licenses/by/4.0/",
                    locator="Printed page 3, Table 1 and equation (4).",
                    quantity="First-order neutral removal through photoionization, per second per neutral.",
                    spectrum_label="quiet-Sun, as printed in Table 1", spatial_domain="illuminated exosphere",
                    lifetime_law="Step survival exp(-k dt), gravitational escape is a separate loss.",
                    measurement_uncertainty=None, daughter_state_probabilities=None,
                    limitation="The model removes the neutral. The neon rate differs from the D04 quiet-Sun record. The corresponding solar conditions have not been reconciled.")
    # The article itself carries its authors' CC BY 4.0 notice on page 1.
    # Its permission is limited to this specific edition and projection.
    rights = [dict(path=PDF, sha256=ARTICLE_SHA, source_uri=URL)]
    article = dict(receipt, id=SOURCE, role="source-native", path=PDF, logical_name=PDF, policy=POLICY,
                   rights_evidence=rights, physical_activation=False,
                   quantity_metadata=dict(metadata, source_flux_note="Page 3 gives 1.12e7 in Table 1 and 1.12e-7 in its prose for global Ar-40 supply (m^-2 s^-1). No source flux is projected."))
    table = dict(id=TABLE + ":csv", path=csv_path.name, logical_name=csv_path.name, role="derived-table",
                 source_uri=URL, policy=POLICY, sha256=digest(csv_path), rights_evidence=rights,
                 derived_from=[SOURCE], recipe=dict(path=recipe.name, sha256=digest(recipe)),
                 physical_activation=False, quantity_metadata=metadata,
                 table=dict(id=TABLE, name=csv_path.stem, format="csv", columns=COLUMNS, row_count=3,
                            schema_id="peschel2026-neutral-frequency-v1"))
    packet = dict(schema="amps-chemical-source-packet-v1", id=PACKET, version="1",
                  title="Lunar neutral noble-gas losses in the Peschel et al. model", provider=metadata["citation"],
                  description="Published Table 1 quiet-Sun labels and three loss frequencies.",
                  limitations=[metadata["limitation"]], physical_activation=False, sources=[article, table])
    manifest = output / "source_packet.json"
    manifest.write_text(json.dumps(packet, indent=2) + "\n")
    return manifest


def compare_d04(path):
    """Return exact source-labelled endpoints; retain their own decimal tokens."""
    if digest(path) != D04_SHA:
        raise ValueError("unreviewed D04 snapshot")
    with path.open(newline="") as stream:
        selected = [{"line":line, "species":row["species"], "solar_activity":row["solar_activity"],
                     "frequency_per_s":row["rate_1au_s-1"]}
                    for line, row in enumerate(csv.DictReader(stream), 2) if
                    row["species"] in ("He", "Ne", "Ar") and row["reaction_type"] == "photoionization"
                    and row["distance_au"] == "1.0" and row["optical_depth"] == "0.0"]
    if len(selected) != 6:
        raise ValueError("expected a quiet and active endpoint for each D04 element")
    print(json.dumps(selected, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("fetch").add_argument("acquisition", type=Path)
    build = commands.add_parser("prepare")
    build.add_argument("acquisition", type=Path)
    build.add_argument("output", type=Path)
    replay = commands.add_parser("project")
    replay.add_argument("pdf", type=Path)
    replay.add_argument("output", type=Path)
    commands.add_parser("compare-d04").add_argument("csv", type=Path)
    args = parser.parse_args()
    if args.command == "fetch":
        fetch(args.acquisition)
    elif args.command == "prepare":
        print(prepare(args.acquisition, args.output))
    elif args.command == "project":
        project(args.pdf, args.output)
    else:
        compare_d04(args.csv)

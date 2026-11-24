#!/usr/bin/env python3
"""Project the three rows in the 2010 IUPAC evaluation of OH uptake on ice."""

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
from urllib.request import urlopen

from native_m3_literature import exactly_one, pdf_text, sha256


URL = "https://acp.copernicus.org/articles/10/9059/2010/acp-10-9059-2010.pdf"
FILENAME = "acp-10-9059-2010.pdf"
DIGEST = "2b6dbc6195f42e196ba1e8de83abbe34c0334d8dfc669032f26927c5760397f2"
SOURCE = "iupac2010-oh-ice:evaluation"
TABLE = "iupac2010-oh-ice:uptake"
COLUMNS = ("evidence", "surface_condition", "temperature_min_K", "temperature_max_K",
           "uptake_gamma", "absolute_gamma_plus_minus", "delta_log_gamma_plus_minus")


def fetch(directory):
    """Capture the bounded response and its actual acquisition identity once."""
    directory.mkdir(parents=True, exist_ok=False)
    path = directory / FILENAME
    with urlopen(URL, timeout=35) as response, path.open("xb") as destination:
        # The reviewed publication is under 1 MiB. Reject HTTP challenge pages,
        # larger replacement documents, and newer revisions before use.
        payload = response.read(2**20)
        if response.status != 200 or response.read(1):
            raise ValueError("expected the complete reviewed publisher PDF")
        destination.write(payload)
        receipt = dict(requested_url=URL, final_url=response.geturl(), http_status=response.status,
                       retrieved_utc=datetime.now(timezone.utc).isoformat(), bytes=len(payload))
    receipt["sha256"] = sha256(path)
    (directory / "retrieval.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if receipt["sha256"] != DIGEST:
        raise ValueError("newly acquired publication does not match the reviewed PDF")


def uptake_rows(page):
    """Preserve the different conditions, temperature domains and error conventions."""
    markers = ("V.A1.3", "HO + ice", "Experimental data", "Preferred values",
               "generally 1 Torr", "Comments on preferred values")
    if any(page.count(marker) != 1 for marker in markers):
        raise ValueError("failed to identify the HO + ice data sheet")
    initial, low, high = exactly_one(
        r"^\s*γ0\s*=\s*(\d+\.\d+)\s+(\d+)–(\d+)\s+Cooper and Abbatt, 1996\s+CWFT-RF \(a\)\s*$", page)
    steady, error = exactly_one(r"^\s*γss\s*=\s*(\d+\.\d+)±(\d+\.\d+)\s*$", page)
    preferred, preferred_low, preferred_high = exactly_one(r"^\s*γss\s+(\d+\.\d+)\s+(\d+)–(\d+)\s*$", page)
    # Poppler extracts the printed Greek delta in this PDF font as "1".
    # Keep the published log reliability in its own column. Its confidence
    # level is unspecified and no linear error is derived from it.
    log_error = exactly_one(r"^\s*(?:1|Δ|∆)log\s*\(γ\s*\)\s+±(\d+\.\d+)\s*$", page)
    return (("experimental", "fresh", low, high, initial, "", ""),
            ("experimental", "conditioned", low, high, steady, error, ""),
            ("evaluated", "steady_state_clean_ice", preferred_low, preferred_high,
             preferred, "", log_error))


def prepare(directory, output):
    """Create an offline-replayable packet. Every published file is explicitly identified."""
    path = directory / FILENAME
    receipt = json.loads((directory / "retrieval.json").read_text())
    if (receipt["requested_url"] != URL or receipt["final_url"] != URL or receipt["http_status"] != 200
            or receipt["bytes"] != path.stat().st_size or receipt["sha256"] != DIGEST or sha256(path) != DIGEST):
        raise ValueError("missing or altered OH/ice source or receipt")
    if "CC Attribution 3.0 License." not in pdf_text(path, 1, 1):
        raise ValueError("publisher licence notice is missing from cover")
    rows = uptake_rows(pdf_text(path, 17, 17))
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(path, output / FILENAME)
    helpers = []
    for name in (Path(__file__).name, "native_m3_literature.py"):
        target = output / name
        shutil.copyfile(Path(__file__).with_name(name), target)
        helpers.append(dict(path=name, sha256=sha256(target)))

    table_file = output / "oh-ice-uptake.csv"
    with table_file.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    evidence = dict(path=FILENAME, sha256=DIGEST, source_uri=URL)
    metadata = dict(
        doi="10.5194/acp-10-9059-2010", primary_experiment_doi="10.1021/jp952142z",
        licence="CC BY 3.0", licence_url="https://creativecommons.org/licenses/by/3.0/",
        attribution="J. N. Crowley et al., Atmospheric Chemistry and Physics 10 (2010) 9059–9223. Experimental data from P. L. Cooper and J. P. D. Abbatt (1996).",
        locator="Printed page 9075, section V.A1.3, experimental and preferred values",
        method="Coated-wall flow tube and resonance fluorescence. Typical total pressure 1 Torr. Table estimates approach 0.1 on fresh ice and 0.03 after conditioning.",
        units="Dimensionless uptake and log reliability. Temperature bounds in kelvin.",
        limitation="Irreversible gas-phase loss includes surface chemistry. Incident sticking, thermal accommodation, and regolith transfer require other evidence.",
    )
    common = dict(source_uri=URL, rights_evidence=[evidence], policy="iupac-acp-2010-oh-ice-cc-by-3.0",
                  physical_activation=False)
    article = dict(common, **{k: receipt[k] for k in ("retrieved_utc", "final_url", "http_status")},
                   id=SOURCE, logical_name=FILENAME, role="source-native", path=FILENAME,
                   sha256=DIGEST, quantity_metadata=metadata)
    csv_record = dict(common, id=TABLE + ":csv", logical_name=table_file.name, path=table_file.name,
                      sha256=sha256(table_file), role="derived-table", derived_from=[SOURCE], recipe=helpers[0],
                      table=dict(id=TABLE, name="OH uptake on ice", format="csv", columns=list(COLUMNS),
                                 row_count=len(rows), schema_id="oh-ice-uptake-v1"),
                      quantity_metadata=dict(metadata, recipe_dependencies=helpers[1:],
                                             modification="Reformatted the published dimensionless values; preserved both uncertainty scales and the absence of a fresh-ice uncertainty."))
    packet = dict(schema="amps-chemical-source-packet-v1", id="iupac2010-oh-ice",
                  title="Uptake of OH radicals on ice", provider="IUPAC Task Group, Crowley et al.",
                  description="2010 laboratory summary and evaluated steady-state OH uptake for water ice.",
                  version="2010-V.A1.3", limitations=[metadata["limitation"]],
                  sources=[article, csv_record], physical_activation=False)
    (output / "source_packet.json").write_text(json.dumps(packet, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True, help="Source directory and retrieval.json")
    parser.add_argument("--download", action="store_true", help="Acquire the publication in a new directory first")
    parser.add_argument("--output-dir", type=Path, required=True, help="New output packet directory")
    args = parser.parse_args()
    if args.download:
        fetch(args.acquisition)
    prepare(args.acquisition, args.output_dir)

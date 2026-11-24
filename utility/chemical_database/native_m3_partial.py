#!/usr/bin/env python3
"""Prepare source-preserving projections of reviewed atomic data offline."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

from native_m3_literature import exactly_one, pdf_text, sha256
from native_workbooks import decode_xlsx


SHEETS = {
    "2_He": ("He", "B"),
    "10_Ne": ("Ne", "BCD"),
    "11_Na": ("Na", "BCDE"),
    "18_Ar": ("Ar", "BCDE"),
}
CC_BY_4 = "https://creativecommons.org/licenses/by/4.0/"


def subshell_rows(path):
    """Keep each source numeric token and its position, including blank σ."""
    cells, _ = decode_xlsx(path.read_bytes())
    if not any(row[1] == "Explanation of Tables" and row[3] == "A15" and
               row[6].startswith("The atomic subshell photoionisation cross sections are given in Mb")
               for row in cells):
        raise ValueError("missing native megabarn convention")
    selected = {(row[1], row[3]): row for row in cells if row[1] in SHEETS}
    for sheet, (atom, columns) in SHEETS.items():
        if selected[sheet, "A1"][6] != "Photon Energy / eV":
            raise ValueError("unexpected photon-energy label in " + sheet)
        # Workbook explanation specifies Mb. Only the native photon rows
        # are selected. Row 19 gives occupancies, not cross sections.
        for index in range(2, 18):
            energy_cell = f"A{index}"
            energy = selected[sheet, energy_cell]
            for column in columns:
                sigma_cell = f"{column}{index}"
                sigma = selected[sheet, sigma_cell]
                orbital = selected[sheet, column + "1"][6]
                if (not re.fullmatch(r"[1-9][spdf]", orbital) or not energy[5] or
                        energy[4] != "n" or sigma[4] != "n" or energy[7] or sigma[7]):
                    raise ValueError(f"unexpected workbook encoding in {sheet}/{sigma_cell}")
                yield atom, orbital, energy[5], sigma[5], sheet, energy_cell, sigma_cell


def sutton_ratios(path):
    """Read printed p. 31 and retain both final ratios and width sensitivity."""
    energy = exactly_one(r"by (\d+)-eV electrons", pdf_text(path, 2, 4))
    page = pdf_text(path, 36, 36)
    table, final = page.split("giving the final results", 1)
    # The PDF text layer inserts an internal space in the left 0.81
    # cell and misreads 1σ as "l a". Compare the resulting numeric
    # strings and the orbital order with the printed page.
    widths = exactly_one(r"T = (2\.0) ns\s+T = (2\.5) ns\s+T = (3\.0) ns", table)
    width_sensitivity = re.findall(
        r"^ *(\d+(?:\.\d ?\d)?) {2,}(\d+(?:\.\d+)?) {2,}(\d+(?:\.\d+)?) *$", table, re.M)
    estimates = re.findall(r"=\)?\s*(\d+(?:\.\d+)?) \(with an estimated l a error of (\d+) percent\)", final)
    orbitals = ("4p", "4d", "4f")
    if widths != ("2.0", "2.5", "3.0") or len(width_sensitivity) != 3 or len(estimates) != 3:
        raise ValueError("unexpected Sutton table layout")
    for orbital, width_row, (ratio, error_percent) in zip(orbitals, width_sensitivity, estimates):
        sensitivity = tuple(value.replace(" ", "") for value in width_row)
        if ratio != sensitivity[1]:
            raise ValueError("Sutton final ratio differs from the reported 2.5 ns result")
        yield energy, orbital, "4s", ratio, error_percent, *sensitivity


def bhatia_hybrid_rows(path):
    page = pdf_text(path, 6, 6)
    table = page.split("Table 2. Photoionization cross sections (Mb)", 1)[1]
    table = table.split("It should be pointed out", 1)[0]
    # Requiring the other three printed columns prevents a row-index
    # shift. Their values are separate cited calculations or experiments.
    rows = re.findall(r"^\s*(\d\.\d)\s+(\d+\.\d+)(?:\s+(?:\d+\.\d+|-)){3}\s*$", table, re.M)
    if tuple(row[0] for row in rows) != (
            "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8",
            "0.9", "1.0", "1.1", "1.3", "1.4", "1.5", "1.6"):
        raise ValueError("unexpected Bhatia Table 2 momentum grid")
    return rows


def bhatia_excited_ratios(path):
    page = pdf_text(path, 12, 12)
    table = page.split("Table 8. Ratios", 1)[1]
    helium = table.split("\n", 1)[1].split(" He\n", 1)[1].split(" Li", 1)[0]
    # Keep the printed mantissas and exponent signs. This is a source
    # record of the exchange approximation, including two unconfirmed R2
    # exponents that differ by a decade from adjacent entries.
    rows = re.findall(r"^ *(\d\.\d)\s+(\d+\.\d+) \(([−-]\d+)\)\s+(\d+\.\d+) \(([−-]\d+)\) *$",
                      helium, re.M)
    if tuple(row[0] for row in rows) != ("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"):
        raise ValueError("unexpected helium momentum grid in Bhatia Table 8")
    return rows


def prepare(acquisition, output):
    """Build an explicit packet from successful, hash-verified acquisitions."""
    native_dir = output / "native"
    native_dir.mkdir(parents=True, exist_ok=False)
    helpers = ("native_m3_partial.py", "native_m3_literature.py", "native_workbooks.py")
    for name in helpers:
        shutil.copyfile(Path(__file__).with_name(name), output / name)
    recipe, pdf_helper, xlsx_helper = [
        {"path": name, "sha256": sha256(output / name)} for name in helpers
    ]
    receipts = {directory: json.loads((directory / "retrieval.json").read_text())
                for directory in (acquisition, acquisition / "followup")}
    sources = []

    def acquire(name, followup=False):
        directory = acquisition / "followup" if followup else acquisition
        receipt = receipts[directory][name]
        path = directory / receipt["file"]
        digest = sha256(path)
        if (receipt["status"] != "RETRIEVED" or receipt["http_code"] != "200" or
                digest != receipt["sha256"] or path.stat().st_size != receipt["bytes"]):
            raise ValueError("missing or altered native source: " + name)
        target = native_dir / name
        if not target.exists():
            shutil.copyfile(path, target)
        return target, dict(path=str(target.relative_to(output)), sha256=digest,
                            source_uri=receipt["requested_url"]), receipt

    def record(native, identifier, metadata, *, policy="metadata-only", rights=()):
        path, evidence, receipt = native
        row = dict(evidence, id=identifier, logical_name=path.name, policy=policy,
                   role="source-native" if rights else "source-metadata",
                   rights_evidence=list(rights), retrieved_utc=receipt["retrieved_utc"],
                   final_url=receipt["destination"], http_status=200,
                   quantity_metadata=metadata, physical_activation=False)
        sources.append(row)
        return row

    def project(parent, identifier, columns, rows, metadata):
        # A queryable view always keeps both the licensed original and a
        # pinned recipe. Only this delimiter/column projection is new.
        name = identifier.replace(":", "-")
        path = output / (name + ".csv")
        values = list(rows)
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(values)
        sources.append(dict(id=identifier + ":csv", logical_name=path.name,
                            path=path.name, sha256=sha256(path), source_uri=parent["source_uri"],
                            role="derived-table", policy=parent["policy"],
                            derived_from=[parent["id"]], rights_evidence=parent["rights_evidence"],
                            recipe=recipe, physical_activation=False,
                            quantity_metadata=dict(parent["quantity_metadata"], **metadata,
                                                   recipe_dependencies=[pdf_helper, xlsx_helper]),
                            table=dict(id=identifier, name=name, format="csv", columns=list(columns),
                                       row_count=len(values), schema_id=name + "-v1")))

    figshare_record, figshare_rights, _ = acquire("kalha-figshare.json")
    workbook = acquire("kalha-y-l-workbook.xlsx", followup=True)
    version = json.loads(figshare_record.read_text())
    file, = version["files"]
    # Match both the archive's version and its producer-supplied file identity.
    if (version["id"] != 12389750 or version["version"] != 2 or version["license"]["url"] != CC_BY_4 or
            file["id"] != 22860461 or file["is_link_only"] or file["download_url"] != workbook[1]["source_uri"] or
            file["size"] != workbook[0].stat().st_size or
            file["supplied_md5"] != hashlib.md5(workbook[0].read_bytes(), usedforsecurity=False).hexdigest()):
        raise ValueError("unexpected Figshare v2 producer record or workbook")
    parent = record(workbook, "kalha-yeh-lindau-v2:workbook",
                    dict(doi="10.6084/m9.figshare.12389750.v2", theory_doi="10.1016/0092-640X(85)90016-6",
                         attribution="C. Kalha et al., from theoretical tables by J. J. Yeh and I. Lindau.",
                         licence="CC BY 4.0", licence_url=CC_BY_4, native_md5=file["supplied_md5"],
                         method="Producer transcription of theoretical atomic subshell photoionization cross sections"),
                    policy="figshare-12389750-v2-cc-by-4.0", rights=[figshare_rights])
    project(parent, "kalha-yeh-lindau-v2:he-ne-na-ar",
            ("neutral_atom", "subshell", "photon_energy_eV_native", "photoionization_sigma_Mb_native",
             "sheet", "energy_cell", "cross_section_cell"), subshell_rows(workbook[0]),
            dict(locator="Four element sheets, rows 2–17, and units from Explanation of Tables",
                 limitation="Blank is native. Subshell cross sections do not give state-resolved daughter probabilities."))

    nasa_record, nasa_rights, _ = acquire("sutton-ntrs.json", followup=True)
    report = acquire("sutton-ntrs.pdf", followup=True)
    ntrs = json.loads(nasa_record.read_text())
    if (ntrs["id"] != 19720018903 or ntrs["distribution"] != "PUBLIC" or
            ntrs["disseminated"] != "DOCUMENT_AND_METADATA" or ntrs["copyright"]["containsThirdPartyMaterial"] or
            ntrs["copyright"]["determinationType"] != "GOV_PUBLIC_USE_PERMITTED"):
        raise ValueError("unexpected NASA source-specific public-use determination")
    parent = record(report, "sutton1972-he:report",
                    dict(ntrs_id=19720018903, attribution="J. F. Sutton, NASA Goddard Space Flight Center.",
                         method="Fast delayed-coincidence measurement of simultaneous He ionization and excitation"),
                    policy="ntrs-19720018903-public-use", rights=[nasa_rights])
    project(parent, "sutton1972-he:excitation-ratios",
            ("incident_electron_energy_eV_native", "numerator_He_plus_orbital", "denominator_He_plus_orbital",
             "cross_section_ratio_native", "estimated_relative_one_sigma_percent_native",
             "ratio_T_2_0_ns_native", "ratio_T_2_5_ns_native", "ratio_T_3_0_ns_native"),
            sutton_ratios(report[0]),
            dict(locator="Printed pp. 31–32, Table 2 and final ratios",
                 limitation="The author warns of a possible factor 2–3 bias in 4p. No absolute excited fraction or recoil distribution."))

    article = acquire("bhatia2022.pdf")
    parent = record(article, "bhatia2022-he:article",
                    dict(doi="10.3390/atoms10010009", attribution="A. K. Bhatia.", licence="CC BY 4.0",
                         licence_url=CC_BY_4, method="Hybrid-theory photon ionization of ground-state He"),
                    policy="mdpi-atoms10010009-cc-by-4.0", rights=[article[1]])
    project(parent, "bhatia2022-he:hybrid-table2", ("outgoing_electron_k_native", "photoionization_sigma_Mb_native"),
            bhatia_hybrid_rows(article[0]),
            dict(locator="Table 2, Hybrid Theory [8] column",
                 k_convention="Outgoing electron momentum k in the article's Rydberg convention, not photon energy",
                 limitation="No assigned model uncertainty or residual excited-state branching."))
    project(parent, "bhatia2022-he:excited-table8",
            ("outgoing_electron_k_native", "R1_mantissa_native", "R1_power10_native",
             "R2_mantissa_native", "R2_power10_native"), bhatia_excited_ratios(article[0]),
            dict(locator="Section 3, Table 8, neutral He rows",
                 method="Author's exchange approximation for photon ionization with excitation",
                 k_convention="Native outgoing electron momentum convention",
                 ratios="R1 = sigma(2s)/sigma(1s); R2 = [sigma(2s)+sigma(2p)]/sigma(1s)",
                 limitation="Preserves R2 exponents at k=0.3 and 0.4 as printed. Verify independently before fitting, interpolating or weighting. Section 3 reports poor agreement with more detailed calculations."))

    # The observational/physical applicability of these leads requires more
    # than a descriptive record. Their bytes cannot be used as source tables.
    archive = acquire("sciencedb-v1-en.html", followup=True)
    if '"conditionsOfAccess":"EMBARGO"' not in archive[0].read_text():
        raise ValueError("review current Science Data Bank access state before replaying")
    record(archive, "wang2026-he:sciencedb-metadata",
           dict(doi="10.57760/sciencedb.30583", article_doi="10.1088/1674-1056/ae64d6",
                archive_conditions_of_access="EMBARGO",
                limitation="Workbook files unavailable. Laboratory angular convention described in archive metadata."))
    record(acquire("suarez-pubmed.xml", followup=True), "suarez2018-water:publication-metadata",
           dict(doi="10.1039/C8CP03725F", pmid=30417901,
                limitation="Prepared H2O+ wave packet surfaces do not give the first-electron and fragment joint distributions."))

    packet = dict(schema="amps-chemical-source-packet-v1", id="m3-atomic-partial-20261124",
                  title="Native atomic subshell photoionization and residual helium-ion measurements",
                  description="Four reviewed numerical source projections and two metadata-only leads",
                  provider="Figshare, NASA NTRS, MDPI, Science Data Bank and Europe PMC", version="20261124",
                  limitations=["Source access only. No broadband state mapping, recoil kernel or application profile."],
                  sources=sources, physical_activation=False)
    (output / "source_packet.json").write_text(json.dumps(packet, indent=2) + "\n")
    return {row["table"]["id"]: row["table"]["row_count"] for row in sources if "table" in row}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.acquisition, args.output_dir), indent=2))


if __name__ == "__main__":
    main()

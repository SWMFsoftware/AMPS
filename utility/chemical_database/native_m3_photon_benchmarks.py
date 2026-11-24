#!/usr/bin/env python3
"""Reproduce the Moore–Mapaye OH/H2O fits and Laurell Ar fitted peak centers."""

import argparse
import csv
import json
from pathlib import Path
import re
import shutil

from native_m3_literature import pdf_text, sha256
from native_workbooks import decode_xlsx


def solar_rate_rows(paper):
    """Retain the printed coefficient tokens, including the shared A0 exponent."""
    page = pdf_text(paper, 6, 6)
    if "Coefficients of Power-law Fits" not in page or "2004–2014" not in page:
        raise ValueError("missing solar cycle Table 1")
    # The base species is printed only on the first line of each group.
    # Read all 11 rows in the OH and H2O groups. Keep the product states
    # exactly as printed, even where the author omits a state assignment.
    observed = []
    parent = None
    for line in page.splitlines():
        match = re.match(r"^\s*(OH|H 2O)?\s*\[(\d+)\]\s+(.+?)\s+"
                         r"\(([\d.]+) ± ([\d.]+)\) × 10−(\d+)\s+([\d.]+)\s+±\s+([\d.]+)\s*$", line)
        if match and 6 <= int(match[2]) <= 16:
            if match[1]:
                parent = match[1].replace(" ", "")
            observed.append((match[2], parent, match[3], match[4], match[5],
                             "-" + match[6], match[7], match[8]))
    if [(int(r[0]), r[1]) for r in observed] != [
            (i, "OH" if i < 10 else "H2O") for i in range(6, 17)]:
        raise ValueError("incomplete OH or H2O reaction group in Table 1")
    return observed


def argon_centers(workbook):
    """One pair per posterior sample. Use displayed shared strings, never their indices."""
    cells, book = decode_xlsx(workbook.read_bytes(), workbook.name)
    if ([sheet["dimension"] for sheet in book["worksheets"]] != ["A1:B801"] or
            len(cells) != 1602 or [row[6] for row in cells[:2]] != ["Center x_0", "Center x_1"]):
        raise ValueError("unexpected Extended Data Fig. 4 source layout")
    rows = []
    for excel_row in range(2, 802):
        pair = cells[2 * (excel_row - 1):2 * excel_row]
        if ([cell[3] for cell in pair] != [f"A{excel_row}", f"B{excel_row}"] or
                any(cell[4] != "s" or cell[7] or not re.fullmatch(r"\d+\.\d+", cell[6])
                    for cell in pair)):
            raise ValueError(f"missing or ambiguous argon fit sample at row {excel_row}")
        rows.append((excel_row, pair[0][6], pair[1][6]))
    return rows


def prepare(root, output):
    """Package only successful HTTP 200 objects with matching byte counts and digests."""
    output.mkdir(parents=True, exist_ok=False)
    sources = []
    recipe_files = []
    for name in (Path(__file__).name, "native_m3_literature.py", "native_workbooks.py"):
        dest = output / name
        shutil.copyfile(Path(__file__).with_name(name), dest)
        recipe_files.append(dict(path=name, sha256=sha256(dest)))

    def fetch(group, name):
        directory = root / group
        receipt = json.loads((directory / "retrieval.json").read_text())[name]
        original = directory / receipt["file"]
        if (receipt["http_code"] != "200" or receipt["status"] != "RETRIEVED" or
                sha256(original) != receipt["sha256"] or original.stat().st_size != receipt["bytes"]):
            raise ValueError("unauthenticated source " + name)
        shutil.copyfile(original, output / name)
        return dict(path=name, sha256=receipt["sha256"], source_uri=receipt["requested_url"]), receipt

    def add_source(group, name, identifier, policy, metadata, rights=()):
        file, receipt = fetch(group, name)
        source = dict(file, id=identifier, logical_name=name, policy=policy,
                      role="source-native" if rights else "source-metadata",
                      rights_evidence=list(rights), http_status=200, final_url=receipt["destination"],
                      retrieved_utc=receipt["retrieved_utc"], quantity_metadata=metadata, physical_activation=False)
        sources.append(source)
        return source

    def add_table(parent, identifier, columns, rows):
        filename = identifier.replace(":", "-") + ".csv"
        with (output / filename).open("x", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
        sources.append(dict(id=identifier + ":csv", path=filename, logical_name=filename,
                            sha256=sha256(output / filename), role="derived-table",
                            derived_from=[parent["id"]], policy=parent["policy"], source_uri=parent["source_uri"],
                            rights_evidence=parent["rights_evidence"], recipe=recipe_files[0],
                            quantity_metadata=dict(parent["quantity_metadata"], recipe_dependencies=recipe_files[1:]),
                            table=dict(id=identifier, name=filename[:-4], format="csv", columns=columns,
                                       row_count=len(rows), schema_id=filename[:-4] + "-v1"),
                            physical_activation=False))

    notice, _ = fetch("retrieved", "moore.json")
    moore = add_source("retrieved", "moore.pdf", "moore-mapaye2025:article",
                       "aas-moore-mapaye-2025-cc-by-4.0", dict(
                           attribution="Luke Moore and Russell Mapaye, Planet. Sci. J. 6, 131 (2025), doi:10.3847/PSJ/add6a2. CC BY 4.0.",
                           locator="Table 1, rows 6–16, page 6",
                           formula="j = A0 * (F10.7P in sfu)^A1, j in s^-1. A0 and its uncertainty share the printed power of ten.",
                           radiation_field="1 au optically thin, TIMED/SEE level 3 v12, 1 nm, 2004–2014. F10.7P is the mean of daily F10.7 and its 81-day running average.",
                           scope="Fits to frequencies computed with PHIDRATES cross sections. No independent daughter kinematics or fit covariance supplied."), [notice])
    licence = json.loads((output / notice["path"]).read_text())["message"]
    if (licence["DOI"].casefold() != "10.3847/psj/add6a2" or not any(item["content-version"] == "vor" and
            item["URL"] == "https://creativecommons.org/licenses/by/4.0/" for item in licence["license"])):
        raise ValueError("missing matching photorate publication licence")
    add_table(moore, "moore-mapaye2025:oh-water-fits",
              ["source_reaction_number", "parent", "products_native", "A0_mantissa", "A0_uncertainty_mantissa",
               "A0_scale_power_of_ten", "A1", "A1_uncertainty"], solar_rate_rows(output / moore["path"]))

    notice, _ = fetch("retrieved", "laurell-nature.html")
    argon = add_source("supplement", "laurell-ext4.xlsx", "laurell2025:ext4-workbook",
                       "nature-laurell-ext4-cc-by-4.0", dict(
                           attribution="H. Laurell et al., Nat. Photon. 19, 352–357 (2025), doi:10.1038/s41566-024-01607-8. CC BY 4.0.",
                           source="Extended Data Fig. 4 source workbook. Source experiment DOI 10.5878/m4fc-9a58.",
                           axes="Laboratory photoelectron kinetic energies in eV on the KRAKEN probe scale.",
                           sample="Paired centers from 800 Bayesian fit iterations for one experiment. The mean difference constrains the Ar+ fine-structure splitting.",
                           limitation="Neither column is an absolute cross section or a channel probability. The paper reports 173 +/- 9 meV."), [notice])
    html = (output / notice["path"]).read_text()
    if argon["source_uri"] not in html or "creativecommons.org/licenses/by/4.0" not in html:
        raise ValueError("missing publisher source-data linkage or CC BY notice")
    add_table(argon, "laurell2025:argon-peak-centers",
              ["source_excel_row", "center_x0_eV", "center_x1_eV"], argon_centers(output / argon["path"]))
    add_source("retrieved", "laurell-archive.json", "laurell2024:snd-record", "metadata-only",
               dict(doi="10.5878/m4fc-9a58", scope="Source archive of He and Ar tomography and instrument response. This catalogue entry identifies the underlying experiment."))
    packet = dict(schema="amps-chemical-source-packet-v1", id="photon-benchmarks-20261124",
                  title="Partial water and OH photorates and argon fine-structure photoelectron energies",
                  description="Published solar-rate fits and narrow-band paired laboratory fit samples.",
                  provider="Moore and Mapaye; Laurell et al.", version="20261124",
                  limitations=["Separate native solar and laboratory conditions. No solar-resolved joint product distribution."],
                  physical_activation=False, sources=sources)
    (output / "source_packet.json").write_text(json.dumps(packet, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.acquisition, args.output_dir)

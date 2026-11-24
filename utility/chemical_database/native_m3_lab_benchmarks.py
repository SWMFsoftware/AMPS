#!/usr/bin/env python3
"""Publish reviewed OH loss and He excitation table tokens from exact PDFs."""

import argparse
import csv
import json
from pathlib import Path
import re
import shutil

from native_m3_literature import exactly_one, pdf_text, sha256


def _quantity_cells(text, expected_count):
    """Keep the published estimate and uncertainty; a dash stays unavailable."""
    cells = re.split(r" {2,}", text.strip())
    if len(cells) != expected_count:
        raise ValueError("unexpected number of source table cells")
    for cell in cells:
        if cell == "—":
            yield "", ""
        else:
            yield exactly_one(r"^(\d+(?:\.\d+)?) *± *(\d+(?:\.\d+)?)$", cell)


def wang_tables(path):
    """Project source Tables 2 and 3, with counts and units pinned by scope."""
    page = pdf_text(path, 11, 11)
    condition_table, resolution_table = page.split("Table 3. Impact of the spatial resolution", 1)
    condition_table = condition_table.split("Table 2. The OH diffusive flux", 1)[1]
    if "Dry helium" not in condition_table or "Humid helium" not in condition_table:
        raise ValueError("unexpected Wang Table 2 gas columns")
    conditions = []
    positions = re.findall(r"^ *(\d+\.\d+) mm +(.*)$", condition_table, re.M)
    if tuple(position for position, _ in positions) != ("1.0", "3.5", "5.2"):
        raise ValueError("unexpected Wang Table 2 positions")
    for distance, cells in positions:
        estimates = list(_quantity_cells(cells, 6))
        # The PDF pairs the three columns for each gas. One row per gas
        # keeps flux scales and the source's missing humid entries explicit.
        for medium, offset in (("dry helium", 0), ("humid helium", 3)):
            values = [token for pair in estimates[offset:offset + 3] for token in pair]
            conditions.append((distance, medium, *values))
    resolutions = []
    rows = re.findall(r"^ *(\d+) µm +(.*)$", resolution_table, re.M)
    if tuple(resolution for resolution, _ in rows) != ("10", "50"):
        raise ValueError("unexpected Wang Table 3 resolutions")
    for resolution, cells in rows:
        values = [token for pair in _quantity_cells(cells, 4) for token in pair]
        resolutions.append((resolution, *values))
    return conditions, resolutions


def sutton_kay_ratios(path):
    energy = exactly_one(r"excitation of helium by (\d+) eV electrons", pdf_text(path, 3, 3))
    page = pdf_text(path, 27, 27)
    if "Table III" not in page or "Experiment" not in page:
        raise ValueError("missing experimental column in Sutton–Kay Table III")
    ratios = re.findall(
        r"^Q \(He *\+ *4 2 ([PD])\) +(\+\d+\.\d+) *\n"
        r"Q \(He *\+ *4 2 S\) +(\d+\.\d+) +(-\d+\.\d+) +\d+\.\d+ *$", page, re.M)
    if tuple(orbital for orbital, *_ in ratios) != ("P", "D"):
        raise ValueError("unexpected Sutton–Kay asymmetric uncertainty rows")
    # The native text layer of the scanned last row drops the lower
    # stroke of ±. The printed final column confirms a single symmetric
    # error. Preserve its magnitude and repeat it as both signed bounds.
    f_ratio, f_uncertainty = exactly_one(
        r"^Q \(He *\+ *4 F\) +(\d+\.\d+) [±+](\d+\.\d+) +\d+\.\d+ *$", page)
    return [(energy, "4" + orbital.lower(), "4s", ratio, upper, lower)
            for orbital, upper, ratio, lower in ratios] + [
                (energy, "4f", "4s", f_ratio, "+" + f_uncertainty, "-" + f_uncertainty)]


def prepare(acquisition, output):
    """Authenticate native acquisitions, rights evidence and decoder dependencies."""
    output.mkdir(parents=True, exist_ok=False)
    native = output / "native"
    native.mkdir()
    helpers = []
    for name in (Path(__file__).name, "native_m3_literature.py"):
        target = output / name
        shutil.copyfile(Path(__file__).with_name(name), target)
        helpers.append(dict(path=name, sha256=sha256(target)))
    recipe, dependency = helpers
    receipts = json.loads((acquisition / "retrieval.json").read_text())
    sources = []

    def acquire(filename):
        receipt = receipts[filename]
        path = acquisition / receipt["file"]
        if (receipt["status"] != "RETRIEVED" or receipt["http_code"] != "200" or
                sha256(path) != receipt["sha256"] or path.stat().st_size != receipt["bytes"]):
            raise ValueError("missing or altered primary source: " + filename)
        target = native / filename
        shutil.copyfile(path, target)
        evidence = dict(path=str(target.relative_to(output)), sha256=receipt["sha256"],
                        source_uri=receipt["requested_url"])
        return target, evidence, receipt

    def record(filename, identifier, metadata, policy="metadata-only", rights=()):
        target, evidence, receipt = acquire(filename)
        row = dict(evidence, id=identifier, logical_name=filename, policy=policy,
                   role="source-native" if rights else "source-metadata",
                   rights_evidence=list(rights), retrieved_utc=receipt["retrieved_utc"],
                   final_url=receipt["destination"], http_status=200,
                   quantity_metadata=metadata, physical_activation=False)
        sources.append(row)
        return row, target

    def project(parent, identifier, columns, rows, metadata):
        name = identifier.replace(":", "-")
        path = output / (name + ".csv")
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
        sources.append(dict(id=identifier + ":csv", logical_name=path.name, path=path.name,
                            sha256=sha256(path), source_uri=parent["source_uri"],
                            policy=parent["policy"], rights_evidence=parent["rights_evidence"],
                            role="derived-table", derived_from=[parent["id"]], recipe=recipe,
                            quantity_metadata=dict(parent["quantity_metadata"], **metadata,
                                                   recipe_dependencies=[dependency]),
                            table=dict(id=identifier, name=name, format="csv", columns=list(columns),
                                       row_count=len(rows), schema_id=name + "-v1"),
                            physical_activation=False))

    crossref, crossref_rights, _ = acquire("wang-crossref.json")
    publisher = json.loads(crossref.read_text())["message"]
    if publisher["DOI"] != "10.1088/1361-6595/aceff9" or not any(
            item.get("content-version") == "vor" and
            item["URL"].rstrip("/") == "http://creativecommons.org/licenses/by/4.0"
            for item in publisher["license"]):
        raise ValueError("the reviewed CC BY 4.0 publisher record is required")
    wang, paper = record(
        "wang-quartz.pdf", "wang2023-quartz-oh:article",
        dict(doi=publisher["DOI"], attribution="Jianan Wang and Peter J. Bruggeman. OH density, flux and loss probability measurements in a room temperature atmospheric pressure surface discharge by microscopic laser induced fluorescence. Plasma Sources Sci. Technol. 32 (2023) 085016.",
             licence="CC BY 4.0", licence_url="https://creativecommons.org/licenses/by/4.0/",
             substrate="fused silica (quartz)", conditions="room-temperature impinging atmospheric-pressure helium plasma in ambient air",
             limitation="Effective OH surface loss derived from near-surface flux. Neither vacuum sticking nor kinetic-energy accommodation."),
        "iop-aceff9-cc-by-4.0", [crossref_rights])
    if "Creative Commons Attribution 4.0 licence" not in pdf_text(paper, 2, 2):
        raise ValueError("missing CC BY 4.0 notice in the matching PDF")
    conditions, resolutions = wang_tables(paper)
    flux_columns = ("diffusive_flux_in_1e14_cm_minus2_s_minus1_native",
                    "diffusive_flux_uncertainty_same_units_native",
                    "thermal_flux_in_1e16_cm_minus2_s_minus1_native",
                    "thermal_flux_uncertainty_same_units_native",
                    "loss_probability_beta_native", "loss_probability_uncertainty_native")
    project(wang, "wang2023-quartz-oh:table2", ("Y_mm_native", "plasma_feed", *flux_columns),
            conditions, dict(locator="Section 4, Table 2, printed page 10",
                             missing="The authors did not report diffusive flux or loss at humid Y=1.0 mm."))
    project(wang, "wang2023-quartz-oh:table3",
            ("measurement_resolution_um_native", "OH_density_in_1e12_cm_minus3_native",
             "OH_density_uncertainty_same_units_native", *flux_columns),
            resolutions, dict(locator="Section 4, Table 3, printed page 10",
                              conditions="Dry He at Y=3.5 mm. Paired spatial-resolution analysis of one experiment."))

    nasa_record, nasa_rights, _ = acquire("sutton-kay-ntrs.json")
    ntrs = json.loads(nasa_record.read_text())
    if (ntrs["id"] != 19730006011 or ntrs["distribution"] != "PUBLIC" or
            ntrs["disseminated"] != "DOCUMENT_AND_METADATA" or ntrs["copyright"]["containsThirdPartyMaterial"] or
            ntrs["copyright"]["determinationType"] != "GOV_PUBLIC_USE_PERMITTED"):
        raise ValueError("the reviewed NTRS public-use determination is required")
    nasa, report = record(
        "sutton-kay-ntrs.pdf", "sutton-kay1972-he:preprint",
        dict(ntrs_id=19730006011, report="NASA-TM-X-66137",
             attribution="J. F. Sutton and R. B. Kay, NASA Goddard Space Flight Center.",
             method="Fast delayed coincidence measurement of He+ n=4 relative excitation-ionization cross sections"),
        "ntrs-19730006011-public-use", [nasa_rights])
    project(nasa, "sutton-kay1972-he:table3",
            ("incident_electron_energy_eV_native", "numerator_He_plus_orbital", "denominator_He_plus_orbital",
             "cross_section_ratio_native", "upper_absolute_error_native", "lower_absolute_error_native"),
            sutton_kay_ratios(report),
            dict(locator="Table III, experiment column, printed page 22",
                 limitation="Ratio relative to He+ 4s, no total excited yield or recoil. Same experimental programme as NASA 19720018903."))

    # Descriptive records retain the useful primary leads and their limits.
    # Their numerical payloads are absent, so none can become an active table.
    for filename, identifier, doi_or_id, scope in (
        ("water-pepico-ntrs.json", "water1991-pepico:ntrs", "19920033386",
         "PEPICO of H2O, D2O and water dimer. Parent survival over 10 microseconds or more; catalogue only."),
        ("bellm-crossref.json", "bellm2007-he-ion-states:crossref", "10.1103/PhysRevA.75.042704",
         "Electron differential excitation-ionization ratios for residual He+ n=1–4 at selected settings."),
        ("dogan-crossref.json", "dogan1998-he-ion-line:crossref", "10.1088/0953-4075/31/7/023",
         "Helium excitation-ionization photon-emission measurement from threshold to 300 eV. No reviewed full native table."),
        ("durr-crossref.json", "durr2006-he-recoil:crossref", "10.1103/PhysRevLett.96.243202",
         "Coincidence benchmark at 102 eV and 1 keV. No retrieved differential arrays for Maxwellian integration."),
    ):
        record(filename, identifier, dict(publication_id=doi_or_id, scope=scope))
    packet = dict(schema="amps-chemical-source-packet-v1", id="m3-lab-benchmarks-20261124",
                  title="OH surface loss and helium excitation-ionization laboratory benchmarks",
                  description="Exact evidence-bound table tokens for OH surface loss and relative He+ excitation.",
                  provider="Wang and Bruggeman (IOP), NASA NTRS and publisher metadata", version="20261124",
                  limitations=["These are laboratory references for their stated conditions and metadata access."],
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

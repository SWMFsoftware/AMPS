#!/usr/bin/env python3
"""Register the archival locator for finite-energy H2O fragment-ion yields."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

from native_m3_literature import exactly_one, pdf_text


DOI = "10.1016/0020-7381(76)80088-5"
URL = "https://api.crossref.org/works/" + DOI
PROCEEDINGS = "https://www.asms.org/docs/default-source/proceedings-archive/1975_asms_23rd_conference.pdf"
PROCEEDINGS_SHA = "40d1c49380c06395bf22ce0a228d8ec1efe21fb110955db4f9ce5a1687352678"
RECORD = "kronebusch-crossref.json"


def prepare(acquisition, output):
    """Attach documentary facts to the reviewed original, retaining held-source access."""
    receipts = json.loads((acquisition / "retrieval-with-asms.json").read_text())
    output.mkdir(parents=True, exist_ok=False)

    def held(name, url, identifier, expected_hash):
        source = acquisition / "native" / name
        receipt = receipts[name]
        payload = source.read_bytes()
        if (receipt["status"] != "RETRIEVED" or receipt["http_code"] != "200"
                or receipt["requested_url"] != url or receipt["destination"] != url
                or receipt["bytes"] != len(payload) or receipt["sha256"] != expected_hash
                or hashlib.sha256(payload).hexdigest() != expected_hash):
            raise ValueError("incorrect or changed native document: " + name)
        # The private packet needs the authenticated local file for admission;
        # the shared publisher withholds the payload under metadata-only policy.
        shutil.copyfile(source, output / name)
        return dict(id=identifier, logical_name=name, path=name, sha256=expected_hash,
                    source_uri=url, final_url=receipt["destination"], http_status=200,
                    retrieved_utc=receipt["retrieved_utc"], policy="metadata-only", role="source-metadata",
                    rights_evidence=[], physical_activation=False)

    citation = held(RECORD, URL, "kronebusch-water:crossref",
                    "1b1428ddbb66426121c9a17bbb68619ead35d3be7439e54a4be320fdf68b7f5a")
    article = json.loads((output / RECORD).read_text())["message"]
    if article["DOI"] != DOI or [a["family"] for a in article["author"]] != ["Kronebusch", "Berkowitz"]:
        raise ValueError("unexpected DOI or authors in journal metadata")
    citation["quantity_metadata"] = dict(journal_doi=DOI, related_conference_uri=PROCEEDINGS,
                                         limitation="The journal metadata list a text-and-data-mining licence.")

    proceedings = held("asms-1975.pdf", PROCEEDINGS, "kronebusch-water:asms1975", PROCEEDINGS_SHA)
    # The text layer misreads charge superscripts and spaces. Restrict all
    # extraction to the water section on the visually reviewed scan rather
    # than letting the neighbouring CO2 measurements enter the assessment.
    page = pdf_text(output / "asms-1975.pdf", 117, 117, "-raw")
    before_water, water = page.split("\nH2O\n", 1)
    if "kinetic energy discrimination" not in before_water or "regarded as lower limits" not in before_water:
        raise ValueError("missing kinetic-energy discrimination notice")
    wavelength, energy = exactly_one(r"The (\d+) / \(([\d .]+) eV\).*?photoelectron", water)
    measured = re.findall(r"\((\d+\.\d+) total ionization\)", water)
    if len(measured) != 3 or "competitive channels" not in water or "is very weak" not in water:
        raise ValueError("failed to identify the three water ion yields")
    collection = exactly_one(r"g r e a t e r than (\d+) eV", before_water)
    proceedings["quantity_metadata"] = dict(
        journal_doi=DOI, attribution="P. L. Kronebusch and J. Berkowitz, 23rd ASMS conference (1975), D-1",
        locator="Printed page 114, PDF page 117 of 712, water paragraph and apparatus description",
        parent="isolated H2O", wavelength_angstrom=wavelength, photon_energy_eV=energy.replace(" ", ""),
        assessed_fraction_by_detected_ion=dict(zip(("OH+", "H+", "O+"), measured)),
        denominator="total ionization", fraction_uncertainty=None,
        collection_discriminates_above_ion_kinetic_energy_eV=collection,
        interpretation="The authors identify affected yields as lower limits. The weak O+ signal is distinguished from O+ + 2H at higher incident photon energy.",
        limitation="These are documentary facts from a single photon energy. A joint distribution of electron, ion, neutral and daughter states is unavailable.")
    packet = dict(schema="amps-chemical-source-packet-v1", id="kronebusch-water-photofragmentation",
                  title="Water fragment-ion yields in the extreme ultraviolet",
                  provider="Kronebusch and Berkowitz, ASMS proceedings and Crossref",
                  description="Metadata and archive locators for single-photon water fragment yields.",
                  version="20261124", limitations=[proceedings["quantity_metadata"]["limitation"]],
                  sources=[citation, proceedings], physical_activation=False)
    (output / "source_packet.json").write_text(json.dumps(packet, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.acquisition, args.output_dir)

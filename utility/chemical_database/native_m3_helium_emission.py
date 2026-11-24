#!/usr/bin/env python3
"""Describe reviewed measurements of helium ion line emission in a source packet.

The inputs are locally held publisher originals. Only individually cited
scientific facts and source identities enter shared catalogue metadata.
Reading the local PDFs here grants no permission to distribute them.
"""

import argparse
import json
from pathlib import Path
import re

from native_m3_literature import exactly_one, pdf_text, sha256


# A publisher can generate a fresh dated cover for each download. Releasing
# with these whole-file pins makes the retained copies the offline oracle.
# Changed bytes require a new review even if the bibliographic DOI persists.
ORIGINALS = (
    ("bloemen", "10.1088/0022-3700/14/4/019",
     "f368830ffbb532f1ca4dbe05876e25516aae691ae20823f033100cf1d8e769cd"),
    ("forand", "10.1088/0022-3700/18/7/020",
     "4f3d9e3465cef11a7ec49fe532472222c6c5de36a3f90c941884f03273502c75"),
    ("zatsarinny", "10.1088/0953-4075/47/6/061001",
     "c9c5a3efa5bc7e627585431e64457d5bdd9a76f5b96331c6481e68923d25acd6"),
)


def measurement(energy, wavelength, wavelength_unit, value, uncertainty):
    # The old PDF text layers lose superscript glyphs. The common scale and
    # reported error convention were read from the actual printed pages
    # (Bloemen p. 723 and Forand p. 1416) of the pinned originals.
    return dict(target="neutral helium", projectile="electron",
                observable="He II 2p to 1s photon-emission cross section",
                incident_energy_eV=energy, wavelength=wavelength,
                wavelength_unit=wavelength_unit, cross_section=value,
                reported_absolute_error=uncertainty, cross_section_unit="1e-19 cm2",
                error_confidence_level=None, direct_population_fraction=None,
                limitation="Line emission. No all-state product probability or ion recoil law.")


def describe(name, path):
    if name == "bloemen":
        # Select only the electron-beam experiment in section 4. Section 5
        # reports a different projectile, beam energy, and observed species.
        page = pdf_text(path, 8, 8)
        wavelength, energy, value, uncertainty = exactly_one(
            r"He II (\d+\.\d+) nm radiation at (\d+) eV electron\s+"
            r"energy we thus find \((\d+\.\d+) f (\d+\.\d+)\)", page)
        return dict(measurement(energy, wavelength, "nm", value, uncertainty),
                    attribution="E. W. P. Bloemen and coauthors, 1981, J. Phys. B 14, 717.",
                    method="Modified absolute photoionization chamber.",
                    locator="Section 4, p. 723", polarization_correction="unapplied")
    if name == "forand":
        page = pdf_text(path, 9, 9)
        table = pdf_text(path, 10, 10)
        # In the native scanned text the ± glyph in the p. 1416 central
        # result appears as the digit 1. Read the central value from the
        # 200 eV row on p. 1417, and separately check the verified error
        # token on p. 1416. Do not assign that error to other energies.
        energy, value = exactly_one(r"^ *(200) +(\d+\.\d+) *$", table)
        uncertainty = exactly_one(r"\( 6 \. 2 1 (\d+\.\d+)\) x", page)
        wavelength = exactly_one(r"emission of the He\+ .* line at (\d+) A is", page)
        energy_nodes = re.findall(r"^ *(\d+) +[\dI][ .\d]*$", table, flags=re.M)
        if len(energy_nodes) != 20 or (energy_nodes[0], energy_nodes[-1]) != ("80", "1000"):
            raise ValueError("unexpected energy coverage in Forand Table 2")
        return dict(measurement(energy, wavelength, "angstrom", value, uncertainty),
                    attribution="J. L. Forand, K. Becker and J. W. McConkey, 1985, J. Phys. B 18, 1409.",
                    method="Neutral-He reference radiation with filter-transmission and detector-yield calibration.",
                    locator="Section 3.1.2, p. 1416 and Table 2, p. 1417",
                    cascade_correction="unapplied", polarization_correction="unapplied",
                    published_curve=dict(incident_energy_range_eV=[energy_nodes[0], energy_nodes[-1]],
                                         number_of_points=len(energy_nodes), normalization_energy_eV=energy))
    # This third source compares two published absolute calibrations against
    # integrated calculations. Cite it as theory rather than counting it as
    # a third independent 200 eV measurement.
    page = pdf_text(path, 4, 4)
    if "He+ (2p)" not in page or "He+ (2s)" not in page or "Figure 2." not in page:
        raise ValueError("missing final-state distinction in the theoretical source")
    return dict(attribution="O. Zatsarinny and K. Bartschat, 2014, J. Phys. B 47, 061001.",
                method="B-spline R-matrix with pseudostates",
                locator="Figures 1 and 2, p. 3",
                observable="Angle- and energy-integrated ionization-excitation from He(1s2) into He+ 2p and 2s.",
                result_scope="The calculation supports the original Bloemen and Forand 2p absolute scale.",
                limitation="2s is calculated. Higher ionic states and an ion recoil distribution are outside these figures.")


def prepare(acquisition):
    receipts = json.loads((acquisition / "retrieval.json").read_text(encoding="utf-8"))
    sources = []
    for name, doi, expected in ORIGINALS:
        filename = name + ".pdf"
        relative = "native/" + filename
        path = acquisition / relative
        row = receipts[filename]
        url = "https://iopscience.iop.org/article/" + doi + "/pdf"
        if (row["requested_url"] != url or row["destination"] != url or row["file"] != relative or
                row["status"] != "RETRIEVED" or row["http_code"] != "200" or
                row["sha256"] != expected or not path.is_file() or
                row["bytes"] != path.stat().st_size or sha256(path) != expected):
            raise ValueError("missing or altered reviewed primary source: " + filename)
        metadata = describe(name, path)
        sources.append(dict(id="he-ion-emission:" + name, logical_name=filename, path=relative,
                            sha256=expected, source_uri=url, final_url=url, http_status=200,
                            retrieved_utc=row["retrieved_utc"], policy="metadata-only",
                            role="source-metadata", rights_evidence=[],
                            quantity_metadata=dict(metadata, doi=doi), physical_activation=False))
    return dict(schema="amps-chemical-source-packet-v1", id="he-electron-ion-emission-20261124",
                title="Helium ion photon emission following electron impact on neutral helium",
                provider="Bloemen et al., Forand et al., Zatsarinny and Bartschat (IOP)",
                version="20261124", physical_activation=False, sources=sources,
                description="Two 200 eV experimental anchors and a comparison with integrated 2s/2p theory.",
                limitations=["Original PDFs and full data tables are held pending redistribution permission.",
                             "The experimental observable is line emission at the specified beam energy."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True,
                        help="retained native/ and retrieval.json; the new packet is written beside them")
    parser.add_argument("--check", action="store_true", help="verify exact offline replay of the existing packet")
    args = parser.parse_args()
    packet = args.acquisition / "source_packet.json"
    content = (json.dumps(prepare(args.acquisition), ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    if args.check:
        if packet.read_bytes() != content:
            raise ValueError("existing helium emission packet differs from the reviewed reduction")
    else:
        with packet.open("xb") as out:
            out.write(content)
    print(packet, sha256(packet))


if __name__ == "__main__":
    main()

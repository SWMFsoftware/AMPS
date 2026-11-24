#!/usr/bin/env python3
"""Package exact numeric text from the three reviewed M3 articles."""

import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


NIST_POLICY = "nist-jres-water-1966-kim-1992"
CUI_POLICY = "mdpi-atoms13090074-cc-by-4.0"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pdf_text(path, first, last, mode="-layout"):
    """Read the existing text layer rather than digitizing figures or using OCR."""
    result = subprocess.run(
        ["pdftotext", mode, "-f", str(first), "-l", str(last), str(path), "-"],
        capture_output=True, encoding="utf-8", check=True, timeout=30,
    )
    return result.stdout


def exactly_one(pattern, text):
    matches = re.findall(pattern, text, re.M)
    if len(matches) != 1:
        raise ValueError("missing or ambiguous published text: " + pattern)
    return matches[0]


def water_views(path):
    text = pdf_text(path, 1, 3, "-raw")
    if "Dibeler" not in text or "National Bureau of Standards" not in text:
        raise ValueError("unexpected source for water photoionization")
    # Use the author's observed onsets, with no conversion or inferred eV
    # uncertainty from the threshold location's wavelength uncertainty.
    water_wl, water_error, water_ev = exactly_one(
        r"ascent is (\d+\.\d+) ± (\d+) A \((\d+\.\d+) eV\)", text)
    oh_wl, oh_ev = exactly_one(r"occurs at (\d+) A\((\d+\.\d+) eV\)", text)
    onsets = [("H2O+", water_wl, water_error, water_ev),
              ("OH+", oh_wl, "", oh_ev)]
    # The ratio is an ion signal normalized to the parent, with no correction
    # for detector response versus mass. Preserve it without renormalization.
    wl, ev, h2o, oh, oxygen = exactly_one(
        r"ions of water observed at (\d+)A\s+\((\d+\.\d+) eV\)[\s\S]*?"
        r"namely, (\d+\.\d+), (\d+\.\d+),\s+and (\d+\.\d+)\.", text)
    signals = [("H2O+", wl, ev, h2o), ("OH+", wl, ev, oh), ("O+", wl, ev, oxygen)]
    return onsets, signals


def kim_neutral_fits(path):
    page = pdf_text(path, 2, 2)
    if "Table 1." not in page or "National Institute of Standards and Technology" not in page:
        raise ValueError("unexpected source for Kim Table 1")
    # In the two-column layout, anchor the entire numeric row to its exact
    # terminal columns so it cannot match a prose occurrence of H or He.
    rows = []
    for atom in ("H", "He"):
        cells = exactly_one(r"\b" + atom + r"\s+(\d+\.\d+)\s+(\d+\.\d+)"
                            r"\s+(-\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$", page)
        rows.append((atom, *cells))
    return rows


def cui_capture_fractions(path):
    page = pdf_text(path, 4, 4)
    if "Table 1. Relative cross sections (%)" not in page or "E (keV/u)" not in page:
        raise ValueError("unexpected source for Cui Table 1")
    # Fixed-width text output preserves the genuinely blank excited-state
    # measurement of the first row. An empty cell never becomes a measured 0.
    rows = []
    for line in page.split("Table 1. Relative cross sections (%)", 1)[1].splitlines():
        if not re.match(r"^\s*\d", line):
            continue
        cells = re.split(r"\s{2,}", line.strip())
        if cells[0] == "0.5" and len(cells) == 4 and cells[1] == "100":
            energy, ground, computed_ground, computed_excited = cells
            ground_error = excited = excited_error = ""
        elif len(cells) == 5:
            energy, ground_data, computed_ground, excited_data, computed_excited = cells
            ground, ground_error = ground_data.split(" ± ")
            excited, excited_error = excited_data.split(" ± ")
        else:
            raise ValueError("unexpected native branching-table layout")
        if Decimal(computed_ground) + Decimal(computed_excited) != 100:
            raise ValueError("native calculated branch percentages inconsistent")
        if excited and Decimal(ground) + Decimal(excited) != 100:
            raise ValueError("native measured branch percentages inconsistent")
        rows.append((energy, ground, ground_error, computed_ground, excited, excited_error, computed_excited))
    if tuple(row[0] for row in rows) != ("0.5", "1.25", "2.5", "3.75", "5"):
        raise ValueError("unexpected native projectile energy support")
    return rows


def prepare(acquisition, output):
    """Publish an offline packet after binding every object to its receipt."""
    output.mkdir(parents=True, exist_ok=False)
    (output / "native").mkdir()
    recipe = output / Path(__file__).name
    shutil.copyfile(__file__, recipe)
    sources = []
    catalogues = {
        section: json.loads((acquisition / section / "retrieval.json").read_text())
        for section in ("", "followup")
    }

    def original(section, filename, expected_sha=None):
        receipt = catalogues[section][filename]
        path = acquisition / section / receipt["file"]
        digest = sha256(path)
        if (digest != receipt["sha256"] or (expected_sha and digest != expected_sha) or
                receipt["status"] != "RETRIEVED" or receipt["http_code"] != "200"):
            raise ValueError("missing or changed research source: " + filename)
        target = output / "native" / filename
        if not target.exists():
            shutil.copyfile(path, target)
        return dict(path=str(target.relative_to(output)), sha256=digest,
                    source_uri=receipt["requested_url"]), target

    def add_source(section, filename, identifier, metadata, *, policy="metadata-only", sha=None, rights=()):
        receipt = dict(catalogues[section][filename])
        fields, target = original(section, filename, sha)
        row = dict(fields, role="source-native", id=identifier, policy=policy,
                   logical_name=filename, retrieved_utc=receipt["retrieved_utc"], http_status=200,
                   quantity_metadata=metadata, rights_evidence=list(rights), physical_activation=False)
        sources.append(row)
        return row, target

    def add_table(parent, identifier, name, columns, rows, metadata):
        path = output / (name + ".csv")
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
        sources.append(dict(
            id=identifier + ":csv", logical_name=path.name, path=path.name, sha256=sha256(path),
            source_uri=parent["source_uri"], role="derived-table", policy=parent["policy"],
            derived_from=[parent["id"]], rights_evidence=parent["rights_evidence"],
            recipe=dict(path=recipe.name, sha256=sha256(recipe)),
            quantity_metadata=dict(parent["quantity_metadata"], **metadata),
            table=dict(id=identifier, name=name, format="csv", columns=list(columns),
                       row_count=len(rows), schema_id=name + "-v1"),
            physical_activation=False,
        ))

    nist_rights, _ = original("", "nist-licence.html",
                             "5127d122b9b00973a0a92d98f462595d2dd39e831105d7ffc8d7ad40cd7e4b57")
    water, water_path = add_source(
        "", "water-nbs.pdf", "dibeler1966-water:pdf",
        dict(doi="10.6028/jres.070A.038", attribution="V. H. Dibeler, J. A. Walker, H. M. Rosenstock. Republished courtesy of the National Institute of Standards and Technology.",
             method="gas-phase photon ionization with mass-selected ion signals"),
        policy=NIST_POLICY, sha="caab1a72ff76e4d02fe59c559efd2eff4f2a5542730505eb7fc2776dcc9eff04",
        rights=[nist_rights])
    onsets, signals = water_views(water_path)
    add_table(water, "dibeler1966-water:observed-onsets", "water-photoionization-onsets",
              ("detected_ion", "wavelength_angstrom_native", "wavelength_uncertainty_angstrom_native",
               "onset_eV_native"), onsets,
              dict(locator="Section 3 H2O, printed pp. 460–461",
                   limitation="Measured historical appearance energies. Blank uncertainty is not zero."))
    add_table(water, "dibeler1966-water:relative-ion-signals", "water-relative-ion-signals",
              ("detected_ion", "wavelength_angstrom_native", "photon_energy_eV_native",
               "observed_signal_over_parent_native"), signals,
              dict(locator="Section 3 H2O, printed pp. 460–461",
                   limitation="No mass/structure detector correction. Signals are not branching fractions or neutral identities."))
    kim, kim_path = add_source(
        "followup", "he-kim-correct.pdf", "kim1992-neutral-ionization:pdf",
        dict(doi="10.6028/jres.097.032", attribution="Yong-Ki Kim. Republished courtesy of the National Institute of Standards and Technology.",
             method="integral ionization cross-section fit"),
        policy=NIST_POLICY, sha="b777277069851715e5810c17f36ad8fcbaf005fd6c0071c2b9baba28082be961",
        rights=[nist_rights])
    add_table(kim, "kim1992-neutral-ionization:table1", "kim-neutral-ionization-fit",
              ("neutral_target", "I_eV_native", "A_native", "B_native", "C_native", "D_native"),
              kim_neutral_fits(kim_path),
              dict(locator="Table 1 and Eq. (2), printed p. 690",
                   energy_convention="incident electron energy, eV; use native Eq. (2) domain E < 10 keV",
                   limitation="No residual-state split, slope uncertainty, recoil, or application selection."))
    cui_sha = "7ed7cbf026341d5c1c7be67c25a48000f0b8243a366ba9c5d259040badb922d0"
    cui_rights, _ = original("followup", "cui-mdpi.pdf", cui_sha)
    cui, cui_path = add_source(
        "followup", "cui-mdpi.pdf", "cui2025-he-capture:pdf",
        dict(doi="10.3390/atoms13090074", attribution="Shucheng Cui et al.", licence="CC BY 4.0",
             licence_url="https://creativecommons.org/licenses/by/4.0/",
             reaction="He+ + He, single capture; recoil He+ in the ground state"),
        policy=CUI_POLICY, sha=cui_sha, rights=[cui_rights])
    add_table(cui, "cui2025-he-capture:table1", "he-capture-relative-branches",
              ("projectile_energy_keV_per_u_native", "ground_measured_percent_native",
               "ground_measured_uncertainty_percent_native", "ground_calculated_percent_native",
               "excited_measured_percent_native", "excited_measured_uncertainty_percent_native",
               "excited_calculated_percent_native"), cui_capture_fractions(cui_path),
              dict(locator="Table 1, section 3.1",
                   limitation="Percent of single capture. No assigned uncertainty confidence level or absent measurement inferred. No angular data digitized."))

    # Preserve source identity and applicability without including held pages
    # or mixed third-party data in the shared managed-object store.
    for section, filename, identifier, limitation in (
        ("", "he-differential-adelaide.pdf", "bellm2007-he-ion-states",
         "He+(n=1–4) differential ratios at selected geometries and different incident energies. APS article permission is limited."),
        ("", "superfluid-native.zip", "signalhelium-14173434",
         "Superfluid simulation repository includes embedded third-party cross sections and object arrays. Separate permissions and gas-phase provenance required."),
        ("followup", "photon-he-osti.pdf", "samson-starace-he-photon",
         "Fluorescence identifies He+(2s,2p) partials with density-dependent quenching. No numerical solar crosswalk or source-specific managed-data grant established."),
        ("followup", "groningen-knoop-catalog.html", "groningen-knoop-repository",
         "Repository lists thesis chapters without an original ground-capture angular data deposit. Numerical bins and tails remain unresolved."),
    ):
        add_source(section, filename, "m3-literature-20261124:" + identifier, dict(limitation=limitation))

    packet = dict(schema="amps-chemical-source-packet-v1", id="m3-literature-20261124",
                  title="Source-specific laboratory and fit data for gas-phase reaction products",
                  description="Historical water signals/onsets, neutral electron-ionization fits, and helium single-capture branches",
                  provider="NBS/NIST and Cui et al. (MDPI)", version="20261124",
                  limitations=["Source access only. Joint water breakup data, gas-phase residual ion recoil envelope, and selected sodium angles remain unresolved."],
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

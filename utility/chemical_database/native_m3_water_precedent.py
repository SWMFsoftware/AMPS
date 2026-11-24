#!/usr/bin/env python3
"""Prepare an optional catalogue packet for published water-photofragment model choices."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import subprocess

from native_m3_literature import exactly_one, pdf_text, sha256


PACKET = "water-photofragment-model-precedents-20261124"
REV = "62e670f9d0344ced1c1414a86b40fbba43b373ef"
RAW = f"https://raw.githubusercontent.com/charligolea/EXOSPHID/{REV}/src/"
# Pin the exact Julia program inspected, including each primitive called by
# the top-level dispatcher. The classification below relies on their combined
# behavior. A future correction upstream needs a separately reviewed record.
SOURCES = {
    "H2O.jl": (RAW + "database/species/H2O.jl", "a023116aa657d95d166c0c5891534562d0ca122d957cc72e6be567781026bfed"),
    "photodestruction_logic.jl": (RAW + "photoreactions/photodestruction_logic.jl", "87fe278ddb68274c2315504152c75c78e3474b171187b8d45b466357a1533e9e"),
    "SimplePhotoionisation.jl": (RAW + "photoreactions/SimplePhotoionisation.jl", "3b69cc326d0de8226d25c9378933b08604442373bf2d50bc0322f48c908ae135"),
    "SimplePhotodissociation.jl": (RAW + "photoreactions/SimplePhotodissociation.jl", "3a254895f3fe053aeb285838dbfffeb6fa43df134894aa1cf880ad07f3505f90"),
    "thesis.pdf": ("https://mediatum.ub.tum.de/doc/1784196/1784196.pdf", "fde2941220fb1007c76a7a3e717f12cde2cf38d5cf2775b00271bf0b859943ee"),
    "shinnaka2017.pdf": ("https://iopscience.iop.org/article/10.3847/1538-3881/153/2/76/pdf", "58cbd3d4061ca18265642ffedc8e31d62c6ced55fb40f7dd6d82862b74441aa8"),
}


def acquire(directory):
    """Perform the optional online step, leaving native responses for offline review."""
    directory.mkdir(parents=True, exist_ok=False)
    receipts = {}
    for name, (url, expected) in SOURCES.items():
        path = directory / name
        # Capture the final URL and status from the same transfer as the payload.
        response = subprocess.check_output([
            "curl", "-fsSL", "--max-time", "120", "-o", str(path), "-w",
            "%{response_code}\n%{url_effective}", url], text=True)
        status, final = response.splitlines()
        if status != "200" or sha256(path) != expected:
            raise ValueError(f"unexpected edition or response from {url}")
        receipts[name] = dict(source_uri=url, final_url=final, http_status=int(status),
                              retrieved_utc=datetime.now(timezone.utc).isoformat())
    (directory / "retrieval.json").write_text(json.dumps(receipts, indent=2) + "\n")


def source_metadata(directory):
    """Separate software parameters and a coma-model comparison from laboratory data."""
    source = (directory / "H2O.jl").read_text()
    # Read the last value in each ionisation row: this is the source's
    # shortest wavelength band. These are conditional *model fractions*.
    # Match the producer's comments to preserve species/charge assignments.
    ions = ["H2O(+) + e-", "H + OH(+) + e-", "H2 + O(+) + e-", "H(+) + OH + e-"]
    fractions = []
    thresholds = []
    probabilities = source.split("const rp_H2O = (", 1)[1].split("\n)", 1)[0]
    appearances = source.split("const tsh_H2O = (", 1)[1].split("\n)", 1)[0]
    for products in ions:
        comment = r'#\s*"\s*H2O \+ γ -> ' + re.escape(products) + '"'
        fractions.append(exactly_one(r"^ +\([^\n]*, ([\d.]+)f0\),? +" + comment, probabilities))
        # Dissociative rows start with the intact-intermediate threshold,
        # followed by the total appearance threshold for both heavy products.
        pattern = r"^ +(?:\(12\.60f0, )?([\d.]+)f0\)?,? +" + comment
        thresholds.append(exactly_one(pattern, appearances))
    band = exactly_one(r"\((\d+\.\d+)f0, (\d+\.\d+)f0\)\)\s*\n", source)
    model = dict(commit=REV, evidence_type="published model parameters",
                 locator="src/database/species/H2O.jl, tsh_H2O, rp_H2O and wv_H2O",
                 photon_wavelength_bin_angstrom=list(band),
                 branches=[dict(products=name, fraction=f, threshold_eV=t) for name, f, t in
                           zip(ions, fractions, thresholds)], uncertainty=None)

    # Compare the published comet model in the context of its own observable.
    # Page 3 is the original Table 2. Two missing-speed glyphs are read as L
    # in this PDF's text layer; the printed table shows ellipses in both cells.
    table = pdf_text(directory / "shinnaka2017.pdf", 3, 3)
    a, b, f = exactly_one(r"H2O \+ hν → H \+ OH +1357–1860 +(\d+\.\d+) \(H\) +(\d+\.\d+) \(OH\) +(\d+\.\d+)", table)
    ion_fraction = exactly_one(r"H2O \+ hν → ionization product +<984 +L +L +(\d+\.\d+)", table)
    comet = dict(evidence_type="published comet model inputs", doi="10.3847/1538-3881/153/2/76",
                 locator="Table 2, 67P hydrogen-coma model", observable="neutral H Lyman-alpha images",
                 neutral_channel=dict(products="H + OH", wavelength_angstrom=[1357, 1860],
                                      H_km_s=a, OH_km_s=b, fraction=f),
                 ionization=dict(wavelength_angstrom_less_than=984, fraction=ion_fraction,
                                 product_speeds_km_s=None, state_resolved=False))

    return {"H2O.jl": model, "shinnaka2017.pdf": comet,
            "thesis.pdf": dict(evidence_type="author's model description",
                               locator="2025 TUM thesis, sections 3.2.1 and 4.2, PDF pages 21 and 41",
                               assumptions="Simplified wavelength ranges and branch probabilities. "
                                           "The thesis describes isotropic emission and translational excess energy."),
            "photodestruction_logic.jl": dict(commit=REV, locator="DiPI_logic, lines 229–249",
                                             evidence_type="published model algorithm",
                                             method="Sequential ionization and heavy-pair breakup use the same input photon energy.",
                                             returned_heavy_products=2, returned_electron=False,
                                             single_photon_electron_heavy_budget_demonstrated=False),
            "SimplePhotoionisation.jl": dict(commit=REV, locator="allocate_velocity_ionisation; simulate_photoionisation",
                                             method="Adds photon momentum to parent momentum, samples trial directions in the forward hemisphere, and returns only the ion velocity."),
            "SimplePhotodissociation.jl": dict(commit=REV, locator="allocate_velocity_dissociation; simulate_photodissociation",
                                               method="Adds photon momentum to the intermediate parent and solves for the two heavy product velocities.")}


def prepare(acquisition, output):
    """Replay offline; require every cited native file to match the reviewed bytes."""
    receipts = json.loads((acquisition / "retrieval.json").read_text())
    for name, (url, expected) in SOURCES.items():
        receipt = receipts[name]
        if (receipt["http_status"] != 200 or receipt["source_uri"] != url or
                sha256(acquisition / name) != expected):
            raise ValueError("changed source or acquisition identity: " + name)
    quantities = source_metadata(acquisition)
    output.mkdir(parents=True, exist_ok=False)
    entries = []
    for name, (_, digest) in SOURCES.items():
        # The source packet is a private input to publication. Its catalogue
        # entries carry model facts; they grant no downstream redistribution
        # permissions or activation of a physical event generator.
        shutil.copyfile(acquisition / name, output / name)
        entries.append(dict(id="water-model:" + name, logical_name=name, path=name,
                            sha256=digest, policy="metadata-only", role="source-metadata",
                            rights_evidence=[], physical_activation=False,
                            quantity_metadata=quantities[name], **receipts[name]))
    packet = dict(schema="amps-chemical-source-packet-v1", id=PACKET,
                  title="Water photoionization and fragment treatment in published models",
                  provider="EXOSPHID, Technical University of Munich and Shinnaka et al.",
                  version="20261124", description="Published model assumptions, code behavior and neutral-coma comparator.",
                  limitations=["No measured joint charged-fragment, neutral and electron recoil distribution or its uncertainty."],
                  physical_activation=False, sources=entries)
    manifest = output / "source_packet.json"
    manifest.write_text(json.dumps(packet, indent=2, ensure_ascii=False) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Download and authenticate sources into a new acquisition directory before replay.")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.fetch:
        acquire(args.acquisition)
    print(prepare(args.acquisition, args.output_dir))

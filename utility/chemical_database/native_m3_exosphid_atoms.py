#!/usr/bin/env python3
"""Read the single-ionization inputs in the pinned EXOSPHID noble-gas species records."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess

from native_m3_literature import exactly_one, sha256


DATASET = "exosphid-noble-gas-photoionization"
COMMIT = "62e670f9d0344ced1c1414a86b40fbba43b373ef"
RAW = f"https://raw.githubusercontent.com/charligolea/EXOSPHID/{COMMIT}/src/database/species/"
NIST_AR = "https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=argon&units=1"
DIGESTS = {
    "He": "33266e99326a251edc1cd4f54d8be2972e5f3e5e01f123f557ae7b871cbe35f3",
    "Ne": "f9780f69abf03784b2d122fff1d48817220bc7f78202cf269a78279c900a4b3b",
    "Ar": "f4cd3611cc0eb67ea2f4e25fd3bcb9adcc76007bbc86853d2a5e4fd4c4aca1de",
}


def acquire(directory):
    """Retain immutable native modules and the status of each HTTPS transfer for offline replay."""
    directory.mkdir(parents=True, exist_ok=False)
    receipts = {}
    for atom, expected in DIGESTS.items():
        name = atom + ".jl"
        url = RAW + name
        result = subprocess.check_output([
            "curl", "-fsSL", "--proto", "=https", "--proto-redir", "=https", "--max-time", "45",
            "--output", str(directory / name), "--write-out", "%{http_code}\n%{url_effective}", url,
        ], text=True)
        status, final = result.splitlines()
        if status != "200" or final != url or sha256(directory / name) != expected:
            raise ValueError(f"EXOSPHID source transfer changed for {name}")
        receipts[name] = dict(source_uri=url, final_url=final, http_status=int(status),
                              retrieved_utc=datetime.now(timezone.utc).isoformat())
    (directory / "retrieval.json").write_text(json.dumps(receipts, indent=2) + "\n")


def read_record(path, atom):
    """Read producer literals without loading or executing the Julia module."""
    source = path.read_text(encoding="utf-8")
    # Float32 literals use 'f' for the decimal exponent. Preserve the native
    # spelling/precision and export ordinary decimal strings for science users.
    number = r"(\d+(?:\.\d+)?)f0"
    quiet = exactly_one(rf"^const qr_{atom} = (\d+\.\d+)f(-\d+)\s*$", source)
    active = exactly_one(rf"^const ar_{atom} = (\d+\.\d+)f(-\d+)\s*$", source)
    wavelength = exactly_one(rf"^const wvl_thr_{atom} = {number}\s*$", source)
    low, high = exactly_one(rf"^const wv_{atom} = \(\({number}, {number}\),\)\s*$", source)
    threshold = exactly_one(rf"^const tsh_{atom} = \({number},\)\s*#", source)
    probability = exactly_one(rf"^const rp_{atom} = \(\({number}\),\)\s*$", source)
    products = exactly_one(rf'^const sn_{atom} = \(\(("{atom}", "{atom}\(\+\)", "e\(-\)")\),\)\s*$', source)
    reaction_type = exactly_one(rf'^const rt_{atom} = \("(SPI)",\)\s*$', source)
    reaction_id = exactly_one(rf'^const rn_{atom} = \("({atom}-PI1)",\)\s*$', source)
    if high != wavelength:
        raise ValueError(f"{atom} band upper boundary disagrees with its declared wavelength threshold")
    # A probability of unity selects the only reaction listed for the atom.
    # It gives no measured population weights among the ion's fine-structure
    # or excited electronic levels. Keep this absence explicit for all species.
    result = dict(commit=COMMIT, atom=atom, reaction_id=reaction_id, reaction_type=reaction_type,
                  species_tokens=json.loads("[" + products + "]"), reaction_probability=probability,
                  quiet_rate_s_inv="e".join(quiet), active_rate_s_inv="e".join(active),
                  threshold_eV=threshold, wavelength_cutoff_angstrom=wavelength,
                  wavelength_range_angstrom=[low, high], residual_ion_state_probabilities=None,
                  evidence_type="lunar model input", measurement_uncertainty=None,
                  locator=f"src/database/species/{atom}.jl; wvl_thr, qr, ar, tsh, sn, rp, rn, rt and wv constants")
    if atom == "Ar":
        # NIST supplies the independent ground-state comparator. The EXOSPHID
        # source leaves the product state unspecified. Preserve both values.
        result["ground_threshold_comparison"] = dict(species="Ar I", source_uri=NIST_AR,
            ground_ionization_energy_eV="15.7596119",
            applicability="17.56 eV in this model exceeds the Ar I ground-state ionization threshold. The source gives no excited-state assignment for this difference.")
    return result


def prepare(directory, output):
    """The publication admits model facts as metadata and keeps the native modules private."""
    receipts = json.loads((directory / "retrieval.json").read_text())
    entries = []
    for atom, expected in DIGESTS.items():
        name = atom + ".jl"
        receipt = receipts[name]
        path = directory / name
        if (sha256(path) != expected or receipt["http_status"] != 200 or
                receipt["source_uri"] != RAW + name or receipt["final_url"] != RAW + name):
            raise ValueError(f"original or retrieval identity differs for {name}")
        entries.append(dict(id="exosphid-atoms:" + atom, logical_name=name, path=name, sha256=expected,
            policy="metadata-only", role="source-metadata", rights_evidence=[], physical_activation=False,
            quantity_metadata=read_record(path, atom), **receipt))
    output.mkdir(parents=True, exist_ok=False)
    for entry in entries:
        shutil.copyfile(directory / entry["path"], output / entry["path"])
    packet = dict(schema="amps-chemical-source-packet-v1", id=DATASET, version="20261124",
                  title="Noble gas photoionization in the EXOSPHID lunar exosphere model",
                  provider="EXOSPHID", description="Input constants of three single-channel lunar model records.",
                  limitations=["Residual ion states are unspecified. The Ar model threshold differs from the ground-state value."],
                  physical_activation=False, sources=entries)
    manifest = output / "source_packet.json"
    manifest.write_text(json.dumps(packet, indent=2, ensure_ascii=False) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Download the three reviewed native modules before preparing the packet.")
    args = parser.parse_args()
    if args.fetch:
        acquire(args.acquisition)
    print(prepare(args.acquisition, args.output_dir))

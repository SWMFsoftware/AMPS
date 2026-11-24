#!/usr/bin/env python3
"""Record three documented lunar OH surface-model choices for sensitivity studies."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
from urllib.parse import urlsplit, urlunsplit

from native_m3_literature import exactly_one, pdf_text, sha256


PACKET = "lunar-oh-surface-model-precedents"
JONES = "jones-crossref.json"
SMOLKA = "smolka-lpsc2023-2089.pdf"
LAFERRIERE = "lunar_processes.py"
SOURCES = {
    JONES: ("https://api.crossref.org/works/10.1029/2018GL080008",
            "6997e0a85b328df58ac8cc480a04119c9dbcc496e74d7d678e26b865d60005a0"),
    SMOLKA: ("https://www.hou.usra.edu/meetings/lpsc2023/pdf/2089.pdf",
             "e6e71c615fe7ace61b855655791dfdb51e6b232cfcec60a705bc1f936908e6cf"),
    LAFERRIERE: ("https://ndownloader.figshare.com/files/52841993",
                 "9128965e06ed2fa21b39f2df3f97210ef719fbcbf7b627481ac50cc859c3c10c"),
}


def fetch(directory):
    """Acquire only the three reviewed items. Retain private inputs for offline replay."""
    directory.mkdir(parents=True, exist_ok=False)
    receipts = {}
    for name, (url, expected) in SOURCES.items():
        path = directory / name
        request = subprocess.run(
            ["curl", "-fsSL", "--max-time", "45", "--max-filesize", "1048576",
             "--proto", "=https", "--proto-redir", "=https", "--output", str(path),
             "--write-out", "%{http_code}\n%{url_effective}", url],
            check=True, capture_output=True, text=True)
        code, redirect = request.stdout.splitlines()
        if code != "200" or sha256(path) != expected:
            raise ValueError(f"unexpected HTTP response or changed bytes for {name}")
        # The Figshare redirect contains expiring signatures. Its host and
        # object path identify the destination without storing credentials.
        location = urlsplit(redirect)
        final_url = urlunsplit((location.scheme, location.netloc, location.path, "", ""))
        receipts[name] = dict(source_uri=url, final_url=final_url, http_status=int(code),
                              retrieved_utc=datetime.now(timezone.utc).isoformat(),
                              bytes=path.stat().st_size, sha256=expected)
    (directory / "retrieval.json").write_text(json.dumps(receipts, indent=2) + "\n")


def jones_metadata(path):
    """Identify the DOI deposit and record the separate manual reading of k6/k7."""
    article = json.loads(path.read_text())["message"]
    if article["DOI"].lower() != "10.1029/2018gl080008" or article["author"][0]["family"] != "Jones":
        raise ValueError("unexpected bibliographic deposit for Jones et al.")
    # Crossref supplies the identity and title only. The model values were
    # read in the primary article itself, in the three paragraphs on k5–k9.
    return dict(doi=article["DOI"], title=article["title"][0],
                primary_url="https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GL080008",
                locator="Model description immediately before section 3, terms k5 through k9.",
                provenance="Curated from the primary publisher text. The held Crossref bytes identify the publication; they do not contain the k7 parameters.",
                event="Gas-phase OH photofragments returning to the regolith form bound M–OH by chemisorption (k7).",
                impact_sticking_probability=1.0, activation_barrier_assumed=False,
                sticking_temperature_dependence="constant", energy_accommodation_coefficient=None,
                modeled_radical_escape_fraction=0.5,
                interpretation="The escape fraction acts before surface impact. Surviving radicals return to random sites on the sunlit surface.",
                limitation="The source models second-order surface recombination to molecular water. It provides no intact-OH thermal re-emission or measured incident-OH energy-transfer distribution.")


def smolka_metadata(path):
    """Read the numerical-model paragraph in the author's primary conference report."""
    text = " ".join(pdf_text(path, 1, 1, "-raw").split())
    if "Maxwell-Boltzmann flux distributed particle energies, and full thermal accommodation." not in text:
        raise ValueError("missing Smolka thermal launch or accommodation assumption")
    return dict(locator="LPSC 54 (2023), abstract 2089, page 1, The Numerical Model",
                author="Alexander Smolka", related_journal_doi="10.1016/j.icarus.2023.115508",
                species=["H", "H2", "OH", "H2O"],
                energy_accommodation_coefficient=1.0, impact_sticking_probability=None,
                launch_speed_weight="v^3 exp(-m v^2/(2 k_B T))",
                interpretation="The speed weighting is the three-dimensional thermal effusion flux written per unit speed; the source specifies full accommodation and a Maxwell–Boltzmann flux.",
                geometry="Spherical Moon, subsolar-stationary temperatures.",
                limitation="A surface launch assumption for thermally accommodated particles. The paper supplies no incident OH sticking measurement or recoil spectrum for a freshly produced gas-phase photofragment.")


def laferriere_metadata(path):
    """Inspect the archived author source as text, without running downloaded code."""
    code = path.read_text()
    residence = code.split("def surftime(", 1)[1].split("#%%", 1)[0]
    rate = code.split("def desorption_R21(", 1)[1].split("def desorp_prob(", 1)[0]
    speed = code.split("def maxwell_boltz_dist(", 1)[1].split("def ballistic_distance(", 1)[0]
    # Restrict each extraction to the called function. The module also
    # declares a distinct water-on-water binding energy and lattice frequency.
    barrier = float(exactly_one(r"^\s*Ed = (\d+\.\d+) +# eV", residence))
    coverage = float(exactly_one(r"^\s*f = (\d+\.\d+) +#", residence))
    exponent = int(exactly_one(r"^\s*A = 10\*\*(\d+) +# 1/s", rate))
    if ("tau_surf = 1/ (f * desorp)" not in residence or "desorption_R21(Temperature, Ed)" not in residence
            or "R_des = A*np.exp(-Ed/(kb_eV*T))" not in rate or "(vel**2)*(np.exp(" not in speed
            or "s1 = maxwell_boltz_dist(vel_dist, pMass, temp)" not in code
            or "tau_surf = surftime(results[3, i+1], mole)" not in code or 'elif molecule == "OH":' not in code):
        raise ValueError("author code no longer supports the recorded residence and speed choices")
    return dict(doi="10.1029/2024JE008628", software_doi="10.6084/m9.figshare.28548857.v1",
                locator="Figshare v1 file 52841993, lunar_processes.py: Model_MonteCarlo, surftime, desorption_R21, ballistic, maxwell_boltz_dist.",
                species=["OH", "H2O"], desorption_barrier_eV=barrier,
                attempt_frequency_s_inv=10**exponent, coverage_fraction=coverage,
                residence_time_s="1 / (f A exp(-E_d/(k_B T)))",
                launch_speed_weight="v^2 exp(-m v^2/(2 k_B T))",
                interpretation="In the smooth-surface calculation the particle resides at the landing site, with T in K, then launches with a newly sampled site-temperature speed. f is a dimensionless multiplier in the archived implementation.",
                impact_sticking_probability=None, energy_accommodation_coefficient=None,
                limitation="These values characterize a ballistic hopping simulation. It supplies no separate incident-energy-dependent radical sticking or kinetic-energy accommodation data.")


def prepare(directory, output):
    """Build the metadata packet, retaining checksums and provenance for each interpretation."""
    receipts = json.loads((directory / "retrieval.json").read_text())
    records = []
    for name, key, describe in ((JONES, "jones-k7", jones_metadata),
                                (SMOLKA, "smolka-launch", smolka_metadata),
                                (LAFERRIERE, "laferriere-hops", laferriere_metadata)):
        path = directory / name
        url, expected_hash = SOURCES[name]
        receipt = receipts[name]
        if (receipt["source_uri"] != url or receipt["http_status"] != 200 or receipt["bytes"] != path.stat().st_size
                or receipt["sha256"] != expected_hash or sha256(path) != expected_hash
                or urlsplit(receipt["final_url"]).scheme != "https"):
            raise ValueError("invalid original or receipt for " + name)
        metadata = describe(path)
        metadata.update(evidence_class="published model assumption", measured_lunar_oh_scattering=False)
        records.append(dict(id=PACKET + ":" + key, logical_name=name, path=name, sha256=expected_hash,
                            source_uri=url, final_url=receipt["final_url"], http_status=200,
                            retrieved_utc=receipt["retrieved_utc"], policy="metadata-only",
                            role="source-metadata" if name == JONES else "source-native",
                            rights_evidence=[], quantity_metadata=metadata, physical_activation=False))
    output.mkdir(parents=True, exist_ok=False)
    for name in SOURCES:
        shutil.copyfile(directory / name, output / name)
    packet = dict(schema="amps-chemical-source-packet-v1", id=PACKET, version="1",
                  title="Lunar models of OH capture, thermal emission and residence",
                  provider="Jones et al., Smolka, and Laferriere et al.",
                  description="Three individually attributed surface model choices for gas-phase hydroxyl.",
                  limitations=["No measured paired sticking and translational accommodation law for OH on lunar regolith.",
                               "Each entry describes a separate model. Combining parameters requires an independent physical justification."],
                  sources=records, physical_activation=False)
    target = output / "source_packet.json"
    target.write_text(json.dumps(packet, indent=2, ensure_ascii=False) + "\n")
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Download the reviewed originals into a new acquisition directory")
    args = parser.parse_args()
    if args.fetch:
        fetch(args.acquisition)
    print(prepare(args.acquisition, args.output_dir))

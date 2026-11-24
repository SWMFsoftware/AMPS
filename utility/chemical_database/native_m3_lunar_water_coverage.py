#!/usr/bin/env python3
"""Preserve model inputs used in two distinct lunar H2O calculations by Schörghofer."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess

from native_m3_literature import exactly_one, sha256


DATASET = "schorghofer-lunar-water-models-v123"
REVISION = "ab8baaecb0c8f4414f651b7d31cad610339be0cb"
RELEASE = "10.5281/zenodo.15161597"
FILES = {
    "desorptionrates.f90": ("Lunar", "0c15a23cc629a096913085e7df12f5039bb9e90003696a0efbb0f54b29e33e9f"),
    "montecarlo.f90": ("Exospheres", "ffe2012fae505c7cd57241215f1155fdba22c53f3c310555ba53d38a73eddef0"),
}


def source_uri(name):
    folder, _ = FILES[name]
    return f"https://raw.githubusercontent.com/nschorgh/Planetary-Code-Collection/{REVISION}/{folder}/{name}"


def fetch(destination):
    """Acquire the exact two source-code objects. Preserve the retrieval result for offline replay."""
    destination.mkdir(parents=True, exist_ok=False)
    receipts = {}
    for name, (_, digest) in FILES.items():
        path = destination / name
        url = source_uri(name)
        response = subprocess.run(["curl", "-fsSL", "--max-time", "45", "--max-filesize", "262144",
                                   "--proto", "=https", "--proto-redir", "=https", "--output", str(path),
                                   "--write-out", "%{http_code}\\n%{url_effective}", url],
                                  check=True, capture_output=True, text=True)
        status, final = response.stdout.splitlines()
        if status != "200" or final != url or sha256(path) != digest:
            raise ValueError(f"changed source or unsuccessful retrieval for {name}")
        receipts[name] = dict(source_uri=url, final_url=final, http_status=200, sha256=digest,
                              bytes=path.stat().st_size, retrieved_utc=datetime.now(timezone.utc).isoformat())
    (destination / "retrieval.json").write_text(json.dumps(receipts, indent=2) + "\n")


def active_lines(text):
    """Remove Fortran comments before checking the statements selected in this tag."""
    return [line.partition("!")[0].strip() for line in text.splitlines() if line.partition("!")[0].strip()]


def coverage(text):
    # Split at the routine boundaries. Adjacent spectral and pure-ice functions
    # use some of the same names. They have different domains and equations.
    params = "\n".join(active_lines(text.split("end module desorp", 1)[0]))
    body = text.split("elemental function desorptionrate(T,theta)", 1)[1].split("end function desorptionrate", 1)[0]
    lines = active_lines(body)
    ice, cutoff, width = exactly_one(r"Eice=([0-9.]+), Ec=([0-9.]+), W=([0-9.]+)$", params)
    monolayer = exactly_one(r"thetam=([\de+.]+)$", params)
    nu = exactly_one(r"^nu = ([\de+.]+)/sqrt\(T\)$", "\n".join(lines))
    required = ("v = theta/thetam", "b = 1./(kBeV*T)", "if (v<=0.) then", "S = 0.",
                "elseif (v<=1) then", "S = nu * theta * exp(-Ec*b) *v**(W*b) / (1.+W*b)",
                "gamma = 1./v**2", "S1 = nu*thetam * exp(-Eice*b) * &",
                "& exp(-gamma*(Ec-Eice)*b) / (1 + gamma* W*b)", "S = S1")
    if not all(statement in lines for statement in required):
        raise ValueError("water surface-flux branches differ from the reviewed tagged implementation")
    return dict(locator="Lunar/desorptionrates.f90, module desorp and elemental function desorptionrate",
                study_doi="10.3847/PSJ/acf19b", subsequent_study_doi="10.3847/PSJ/ade5b2",
                mechanism="thermal desorption of adsorbed H2O from heterogeneous lunar-grain sites",
                assumption="Efficient diffusion over the grain surface prioritizes the most strongly binding sites.",
                inputs="T is surface temperature in K, theta is surface coverage in molecules m^-2",
                output="desorbing H2O flux in molecules m^-2 s^-1",
                monolayer_molecules_m2=float(monolayer), Eice_eV=float(ice), Ec_eV=float(cutoff), W_eV=float(width),
                attempt_frequency_prefactor_s_inv=float(nu), attempt_frequency_s_inv="prefactor / sqrt(T/K)",
                coverage_ratio="v = theta / monolayer_molecules_m2", inverse_energy="b = 1/(k_B T), eV^-1",
                zero_coverage_flux=0,
                submonolayer_formula="nu theta exp(-Ec b) v^(W b) / (1+W b)", submonolayer_domain="0 < v <= 1",
                multilayer_formula="nu theta_m exp(-Eice b) exp[-gamma (Ec-Eice) b] / (1+gamma W b)",
                multilayer_domain="v > 1", multilayer_gamma="1/v^2",
                experimental_uncertainty=None, incident_sticking_probability=None,
                kinetic_energy_accommodation=None)


def ballistic(text):
    # Only one of the many commented species and lifetimes is active. The
    # small-hazard branch uses a linear probability, which must stay explicit
    # when comparing this model with a continuous-time survival calculation.
    lines = active_lines(text)
    nu = exactly_one(r"^real\(8\), parameter :: taudissoc = 1/([\d.e+-]+)$", "\n".join(lines))
    required = ("real(8), parameter :: mmass = 18.015", "v(1) = gasdev(idum)*sigma",
                "v(2) = gasdev(idum)*sigma", "v(3) = abs(gasdev(idum))*sigma",
                "if (Q>0. .and. taudissoc>0.) then",
                "destr_rate = flighttime/(taudissoc*semia**2)",
                "if (destr_rate>0.2) destr_rate = 1-exp(-destr_rate)", "p_s = -1")
    if not all(statement in lines for statement in required):
        raise ValueError("species, launch velocity or loss calculation differs from the tagged water model")
    return dict(locator="Exospheres/montecarlo.f90, module exo_species and subroutine hop1",
                study_doi="10.1002/2014GL060820", mechanism="ballistic hopping and in-flight destruction of H2O",
                loss_frequency_s_inv_1au=float(nu), rate_origin="Crovisier (1989) as cited in native code",
                illumination="in-flight loss enabled for Q > 0", distance_scaling="inverse square, distance in au",
                launch_frame="local surface at emission",
                launch_law="two tangential Gaussian components and one outward half-normal component at surface T",
                loss_probability="for integrated hazard h <= 0.2 use h, otherwise use 1-exp(-h)",
                outcome="absorbing parent loss", daughter_states=None, daughter_velocities=None)


def prepare(acquisition, output):
    """Validate private originals and publish queryable parameter metadata from their active routines."""
    receipts = json.loads((acquisition / "retrieval.json").read_text())
    entries = []
    for name, (directory, digest) in FILES.items():
        path = acquisition / name
        receipt = receipts[name]
        if (receipt["source_uri"] != source_uri(name) or receipt["final_url"] != source_uri(name)
                or receipt["http_status"] != 200 or receipt["bytes"] != path.stat().st_size
                or receipt["sha256"] != digest or sha256(path) != digest):
            raise ValueError(f"original or acquisition receipt mismatch for {name}")
        metadata = coverage(path.read_text()) if directory == "Lunar" else ballistic(path.read_text())
        metadata.update(species="H2O", software_doi=RELEASE, revision=REVISION, tag="v1.2.3",
                        evidence_class="published model prescription", qualified_for_incident_OH=False)
        entries.append(dict(id=DATASET + ":" + name, logical_name=name, path=name, role="source-native",
                            policy="metadata-only", rights_evidence=[], physical_activation=False,
                            quantity_metadata=metadata, **receipt))
    # Metadata-only preserves a primary reference even when store policy has
    # no source-specific grant for redistributing upstream Fortran code.
    output.mkdir(parents=True, exist_ok=False)
    for name in FILES:
        shutil.copyfile(acquisition / name, output / name)
    record = dict(schema="amps-chemical-source-packet-v1", id=DATASET, version="1.2.3",
                  title="Lunar H2O coverage-dependent desorption and model photodestruction",
                  provider="Norbert Schörghofer, Planetary Code Collection",
                  description="Two separately published calculations for retained and gas-phase molecular water.",
                  limitations=["The flux presumes efficient grain-surface diffusion. No prediction of incident OH radical sticking, recoil, or energy accommodation.",
                               "Ballistic photodestruction removes the parent. The routine gives no charged-product kinematics."],
                  sources=entries, physical_activation=False)
    target = output / "source_packet.json"
    target.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="retrieve the immutable tagged sources before preparing")
    parser.add_argument("--acquisition", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.fetch:
        fetch(args.acquisition)
    print(prepare(args.acquisition, args.output_dir))

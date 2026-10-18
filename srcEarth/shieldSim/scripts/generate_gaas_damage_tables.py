#!/usr/bin/env python3
"""Canonicalize the downloaded SR-NIEL/ASTM GaAs response pages.

This is a fixture-preparation tool, not a physics calculator.  The upstream
SR-NIEL web forms are session-backed and return transient result URLs, so the
download step and the deterministic canonicalization step are deliberately
separate.  An operator must save result pages produced with the parameters
checked below and then pass those pages to this script.

The charged-particle tables use the SR-NIEL page's published Ga=21.5 eV and
As=21.5 eV displacement thresholds, hadronic plus Coulomb contributions, and
no hadronic threshold scaling.  The neutron table is the ASTM E722-19 GaAs
damage function with its 10 eV threshold.  The 1-MeV-equivalent denominator is
the E722 GaAs value of 70 MeV-mb.  It is kept explicitly distinct from the
95 MeV-mb natural-silicon convention.

Generated files are candidate reference data.  This script does not itself
enable GaAs production outputs or claim accelerator/device validation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pathlib
import shutil
import sys


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from generate_stage5_niel_tables import parse_response  # noqa: E402


AVOGADRO = 6.02214076e23
GA_ATOMIC_WEIGHT = 69.723
AS_ATOMIC_WEIGHT = 74.921595
GAAS_AVERAGE_LATTICE_ATOM_WEIGHT = (GA_ATOMIC_WEIGHT + AS_ATOMIC_WEIGHT) / 2.0
GAAS_REFERENCE_DAMAGE_MEV_MB = 70.0
GAAS_REFERENCE_NIEL = (
    GAAS_REFERENCE_DAMAGE_MEV_MB
    * 1.0e-27
    * AVOGADRO
    / GAAS_AVERAGE_LATTICE_ATOM_WEIGHT
)


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_page(path: pathlib.Path, required_text: list[str]) -> None:
    text = path.read_text(encoding="latin-1")
    missing = [token for token in required_text if token not in text]
    if missing:
        raise RuntimeError(f"{path} does not document required parameters: {missing}")


def write_charged(path: pathlib.Path, rows: list[tuple[float, float]],
                  expected_count: int) -> tuple[float, float]:
    positive = [(energy, niel) for energy, niel in rows if niel > 0.0]
    if len(positive) != expected_count:
        raise RuntimeError(
            f"{path.name}: expected {expected_count} positive nodes, found {len(positive)}"
        )
    lines = ["energy_MeV,niel_MeV_cm2_g"]
    lines.extend(f"{energy:.12e},{niel:.12e}" for energy, niel in positive)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return positive[0][0], positive[-1][0]


def read_edges(si_groups: pathlib.Path) -> list[tuple[float, float]]:
    with si_groups.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    edges = [(float(row["lower_MeV"]), float(row["upper_MeV"])) for row in rows]
    if len(edges) != 640 or edges[0][0] != 1.0e-10 or edges[-1][1] != 20.0:
        raise RuntimeError("the pinned E192 edge source is not the expected 640-group grid")
    if any(edges[index][1] != edges[index + 1][0] for index in range(639)):
        raise RuntimeError("the pinned E192 groups are not exactly contiguous")
    return edges


def write_neutron(path: pathlib.Path, rows: list[tuple[float, float]],
                  edges: list[tuple[float, float]]) -> None:
    if len(rows) != 640:
        raise RuntimeError(f"expected 640 GaAs neutron rows, found {len(rows)}")
    lines = [
        "lower_MeV,upper_MeV,representative_MeV,niel_MeV_cm2_g,"
        "hardness_factor_to_1MeV_neutron"
    ]
    for index, ((representative, niel), (lower, upper)) in enumerate(zip(rows, edges)):
        midpoint = 0.5 * (lower + upper)
        if not math.isclose(representative, midpoint, rel_tol=5.0e-12,
                            abs_tol=1.0e-300):
            raise RuntimeError(
                f"GaAs neutron row {index}: {representative} != E192 midpoint {midpoint}"
            )
        lines.append(
            f"{lower:.12e},{upper:.12e},{representative:.12e},{niel:.12e},"
            f"{niel / GAAS_REFERENCE_NIEL:.12e}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proton-html", type=pathlib.Path, required=True)
    parser.add_argument("--alpha-html", type=pathlib.Path, required=True)
    parser.add_argument("--neutron-html", type=pathlib.Path, required=True)
    parser.add_argument("--si-neutron-groups", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    require_page(args.proton_html, [
        "<SUP>1</SUP>H in Ga", "21.50; 21.50", "Hadron+Coulomb", "No Scaling"
    ])
    require_page(args.alpha_html, [
        "<SUP>4</SUP>He in Ga", "21.50; 21.50", "Hadron+Coulomb", "No Scaling"
    ])
    require_page(args.neutron_html, [
        "Neutrons in GaAs", "ASTM Standard E722-19",
        "Displacement Threshold Energies [eV]: 10"
    ])

    output = args.output.resolve()
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    copied_raw: dict[str, pathlib.Path] = {}
    for name, source in (
        ("sr_niel_v11_1_gaas_ed21p5_proton.html", args.proton_html),
        ("sr_niel_v11_1_gaas_ed21p5_alpha.html", args.alpha_html),
        ("astm_e722_19_gaas_neutron.html", args.neutron_html),
    ):
        destination = raw / name
        shutil.copyfile(source, destination)
        copied_raw[name] = destination

    proton = output / "sr_niel_v11_1_gaas_ed21p5_proton.csv"
    alpha = output / "sr_niel_v11_1_gaas_ed21p5_alpha.csv"
    neutron = output / "astm_e722_19_gaas_neutron_groups.csv"
    proton_domain = write_charged(proton, parse_response(args.proton_html), 121)
    alpha_domain = write_charged(alpha, parse_response(args.alpha_html), 132)
    write_neutron(neutron, parse_response(args.neutron_html),
                  read_edges(args.si_neutron_groups))

    manifest = {
        "schema_version": "shieldSim.damage_response_candidate.v1",
        "preparation_status": "canonical_data_prepared_not_production_enabled",
        "target": {
            "material_id": "GaAs_G4_11_4",
            "description": "stoichiometric gallium arsenide",
            "average_lattice_atom_weight_g_mol": GAAS_AVERAGE_LATTICE_ATOM_WEIGHT,
            "atomic_weights": {"Ga": GA_ATOMIC_WEIGHT, "As": AS_ATOMIC_WEIGHT},
        },
        "charged_particle_definition": {
            "displacement_threshold_eV": {"Ga": 21.5, "As": 21.5},
            "model": "SR-NIEL-7 v11.1 hadronic plus Coulomb",
            "hadronic_threshold_scaling": "none",
            "interpolation": "log_log_between_positive_nodes",
        },
        "n_eq": {
            "convention": "ASTM_E722_19_gallium_arsenide_neutron_contribution_only",
            "reference_damage_function_MeV_mb": GAAS_REFERENCE_DAMAGE_MEV_MB,
            "reference_niel_MeV_cm2_g": GAAS_REFERENCE_NIEL,
            "avogadro_constant_mol_inverse": AVOGADRO,
            "supported_species": ["neutron"],
        },
        "tables": [
            {
                "species": "proton", "path": proton.name,
                "sha256": f"sha256:{sha256(proton)}", "nodes": 121,
                "domain_min_MeV": proton_domain[0], "domain_max_MeV": proton_domain[1],
                "raw_response_sha256": f"sha256:{sha256(copied_raw['sr_niel_v11_1_gaas_ed21p5_proton.html'])}",
            },
            {
                "species": "alpha", "path": alpha.name,
                "sha256": f"sha256:{sha256(alpha)}", "nodes": 132,
                "domain_min_MeV": alpha_domain[0], "domain_max_MeV": alpha_domain[1],
                "raw_response_sha256": f"sha256:{sha256(copied_raw['sr_niel_v11_1_gaas_ed21p5_alpha.html'])}",
            },
            {
                "species": "neutron", "path": neutron.name,
                "sha256": f"sha256:{sha256(neutron)}", "groups": 640,
                "domain_min_MeV": 1.0e-10, "domain_max_MeV": 20.0,
                "raw_response_sha256": f"sha256:{sha256(copied_raw['astm_e722_19_gaas_neutron.html'])}",
            },
        ],
        "sources": {
            "charged_calculator": "https://www.sr-niel.asi.it/index.php/sr-niel-web-calculators/niel-calculator-for-electrons-protons-and-ions/protons-ions-niel-calculator",
            "neutron_calculator": "https://www.sr-niel.asi.it/index.php/sr-niel-web-calculators/niel-calculator-for-neutrons-using-astm-standards/neutrons-niel-calculator",
            "gaas_reference_value": "SAND2022-15938, E722-94 GaAs response, 70 MeV-mb",
            "gaas_reference_url": "https://www.osti.gov/servlets/purl/1898256",
            "gaas_reference_sha256": "sha256:fa116f8ebf4f47b187327c5d384fff6f30c6c6dc84bf70b339c266e0ab11814d",
            "edge_source": "pinned E192 640-group edges from the existing E722-19 silicon package",
            "edge_source_sha256": f"sha256:{sha256(args.si_neutron_groups)}",
        },
        "limitations": [
            "candidate data are not enabled by the production loader",
            "no accelerator/device-degradation validation is claimed",
            "21.5 eV per Ga/As sublattice is distinct from the optional 24 eV effective triple-junction solar-cell convention",
            "the E722-19 GaAs neutron response and 70 MeV-mb GaAs denominator must never be mixed with the silicon response or 95 MeV-mb denominator",
        ],
    }
    manifest_path = output / "candidate_manifest.v1.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    print(json.dumps({
        "manifest": str(manifest_path),
        "reference_niel_MeV_cm2_g": GAAS_REFERENCE_NIEL,
        "table_hashes": {path.name: sha256(path) for path in (proton, alpha, neutron)},
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

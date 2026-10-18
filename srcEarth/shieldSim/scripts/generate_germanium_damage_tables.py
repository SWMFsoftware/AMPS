#!/usr/bin/env python3
"""Canonicalize candidate SR-NIEL proton/He-4 tables for natural germanium.

The SR-NIEL calculator documents an experimentally determined 40.5 eV
displacement threshold for Ge.  No ASTM E722 neutron damage function is
published for Ge, so this tool deliberately emits charged-particle DDD data
only.  It does not invent a neutron table or a germanium ``n_eq`` convention.
The current production damage-library schema requires all three species, so
the resulting package remains candidate data pending a schema that can express
quantity-specific species availability.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import shutil
import sys


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from generate_stage5_niel_tables import parse_response  # noqa: E402


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proton-html", type=pathlib.Path, required=True)
    parser.add_argument("--alpha-html", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    require_page(args.proton_html, [
        "<SUP>1</SUP>H in Ge", "40.50", "Hadron+Coulomb", "No Scaling"
    ])
    require_page(args.alpha_html, [
        "<SUP>4</SUP>He in Ge", "40.50", "Hadron+Coulomb", "No Scaling"
    ])

    output = args.output.resolve()
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    raw_proton = raw / "sr_niel_v11_1_ge_ed40p5_proton.html"
    raw_alpha = raw / "sr_niel_v11_1_ge_ed40p5_alpha.html"
    shutil.copyfile(args.proton_html, raw_proton)
    shutil.copyfile(args.alpha_html, raw_alpha)

    proton = output / "sr_niel_v11_1_ge_ed40p5_proton.csv"
    alpha = output / "sr_niel_v11_1_ge_ed40p5_alpha.csv"
    proton_domain = write_charged(proton, parse_response(args.proton_html), 113)
    alpha_domain = write_charged(alpha, parse_response(args.alpha_html), 130)

    manifest = {
        "schema_version": "shieldSim.damage_response_candidate.v1",
        "preparation_status": "charged_DDD_data_prepared_not_production_enabled",
        "target": {
            "material_id": "Ge_natural_G4_11_4",
            "description": "natural germanium",
        },
        "charged_particle_definition": {
            "displacement_threshold_eV": {"Ge": 40.5},
            "model": "SR-NIEL-7 v11.1 hadronic plus Coulomb",
            "hadronic_threshold_scaling": "none",
            "interpolation": "log_log_between_positive_nodes",
        },
        "tables": [
            {
                "species": "proton", "path": proton.name,
                "sha256": f"sha256:{sha256(proton)}", "nodes": 113,
                "domain_min_MeV": proton_domain[0], "domain_max_MeV": proton_domain[1],
                "raw_response_sha256": f"sha256:{sha256(raw_proton)}",
            },
            {
                "species": "alpha", "path": alpha.name,
                "sha256": f"sha256:{sha256(alpha)}", "nodes": 130,
                "domain_min_MeV": alpha_domain[0], "domain_max_MeV": alpha_domain[1],
                "raw_response_sha256": f"sha256:{sha256(raw_alpha)}",
            },
        ],
        "source": {
            "calculator": "https://www.sr-niel.asi.it/index.php/sr-niel-web-calculators/niel-calculator-for-electrons-protons-and-ions/protons-ions-niel-calculator",
            "version": "SR-NIEL-7 v11.1",
            "parameter_provenance": "calculator page cites Campesato et al. (2019) for Ge Ed=40.5 eV",
        },
        "unavailable": {
            "neutron_DDD": "no approved germanium neutron damage function was found",
            "n_eq": "no ASTM E722 germanium 1-MeV-equivalent convention was found",
        },
        "limitations": [
            "candidate data are not enabled by the production loader",
            "no extrapolation beyond the positive-node domains is permitted",
            "no accelerator/device-degradation validation is claimed",
        ],
    }
    manifest_path = output / "candidate_manifest.v1.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    print(json.dumps({
        "manifest": str(manifest_path),
        "table_hashes": {path.name: sha256(path) for path in (proton, alpha)},
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

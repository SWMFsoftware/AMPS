#!/usr/bin/env python3
"""Authenticate the cross-material reference-data audit and candidate tables.

This test is intentionally deterministic and offline.  It proves provenance
bytes, canonical-table structure, physical-domain bookkeeping, and the
material/quantity status matrix.  It does *not* promote candidate tables to a
production response or replace transport/experimental validation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import pathlib
import sys


if len(sys.argv) != 3:
    raise SystemExit(
        "usage: test_material_function_reference_data.py SOURCE_ROOT EVIDENCE_DIR"
    )

ROOT = pathlib.Path(sys.argv[1]).resolve()
EVIDENCE = pathlib.Path(sys.argv[2]).resolve()
STATUS_PATH = (
    ROOT / "tests/reference_material_functions/material_quantity_data_status.v1.json"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unprefix(value: str) -> str:
    return value[7:] if value.startswith("sha256:") else value


def csv_rows(path: pathlib.Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def authenticate_nist() -> int:
    base = ROOT / "tests/reference_tid/nist"
    manifest = json.loads((base / "manifest.sha256.json").read_text(encoding="utf-8"))
    for relative, record in manifest["files"].items():
        path = base / relative
        require(path.is_file(), f"missing NIST reference file: {relative}")
        require(path.stat().st_size == record["size_bytes"],
                f"NIST byte count changed: {relative}")
        require(sha256(path) == record["sha256"],
                f"NIST SHA-256 changed: {relative}")
    for material in ("germanium", "silicon_dioxide"):
        for particle in ("pstar", "astar"):
            path = base / f"{particle}_{material}.csv"
            header, rows = csv_rows(path)
            require(header == [
                "energy_MeV", "electronic_stopping_MeV_cm2_g",
                "nuclear_stopping_MeV_cm2_g", "total_stopping_MeV_cm2_g",
                "csda_range_g_cm2", "projected_range_g_cm2", "detour_factor",
            ], f"unexpected STAR header: {path.name}")
            energies = [float(row["energy_MeV"]) for row in rows]
            require(len(rows) >= 100 and all(a < b for a, b in zip(energies, energies[1:])),
                    f"invalid STAR energy grid: {path.name}")
            require(all(float(row["electronic_stopping_MeV_cm2_g"]) > 0.0
                        for row in rows), f"nonpositive STAR stopping power: {path.name}")
    return len(manifest["files"])


def authenticate_candidate(relative: str) -> dict:
    manifest_path = ROOT / relative
    base = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(manifest["schema_version"] == "shieldSim.damage_response_candidate.v1",
            f"unexpected candidate schema: {relative}")
    require("not_production_enabled" in manifest["preparation_status"],
            f"candidate accidentally claims production status: {relative}")

    raw_hashes = {sha256(path): path for path in (base / "raw").iterdir() if path.is_file()}
    species: set[str] = set()
    for table in manifest["tables"]:
        require(table["species"] not in species, f"duplicate species in {relative}")
        species.add(table["species"])
        path = base / table["path"]
        require(path.is_file(), f"missing candidate table: {path}")
        require(sha256(path) == unprefix(table["sha256"]),
                f"candidate SHA-256 changed: {path.name}")
        require(unprefix(table["raw_response_sha256"]) in raw_hashes,
                f"raw response is missing or changed for {path.name}")
        header, rows = csv_rows(path)
        count_key = "groups" if table["species"] == "neutron" else "nodes"
        require(len(rows) == int(table[count_key]), f"row count changed: {path.name}")
        if table["species"] == "neutron":
            require(header == ["lower_MeV", "upper_MeV", "representative_MeV",
                               "niel_MeV_cm2_g", "hardness_factor_to_1MeV_neutron"],
                    f"unexpected neutron header: {path.name}")
            reference = float(manifest["n_eq"]["reference_niel_MeV_cm2_g"])
            previous_upper: float | None = None
            for index, row in enumerate(rows):
                lower = float(row["lower_MeV"])
                upper = float(row["upper_MeV"])
                representative = float(row["representative_MeV"])
                niel = float(row["niel_MeV_cm2_g"])
                hardness = float(row["hardness_factor_to_1MeV_neutron"])
                require(lower > 0.0 and upper > lower and lower < representative < upper,
                        f"invalid neutron group {index}: {path.name}")
                if previous_upper is not None:
                    require(lower == previous_upper,
                            f"noncontiguous neutron groups at {index}: {path.name}")
                previous_upper = upper
                require(math.isclose(hardness, niel / reference, rel_tol=5.0e-12),
                        f"wrong hardness factor at {index}: {path.name}")
            first, last = float(rows[0]["lower_MeV"]), float(rows[-1]["upper_MeV"])
        else:
            require(header == ["energy_MeV", "niel_MeV_cm2_g"],
                    f"unexpected charged-particle header: {path.name}")
            energies = [float(row["energy_MeV"]) for row in rows]
            require(all(a < b for a, b in zip(energies, energies[1:])),
                    f"nonmonotonic charged-particle nodes: {path.name}")
            require(all(float(row["niel_MeV_cm2_g"]) > 0.0 for row in rows),
                    f"nonpositive charged-particle NIEL: {path.name}")
            first, last = energies[0], energies[-1]
        require(first == float(table["domain_min_MeV"]) and
                last == float(table["domain_max_MeV"]),
                f"declared domain does not match table endpoints: {path.name}")
    return {"target": manifest["target"]["material_id"],
            "species": sorted(species), "tables": len(manifest["tables"])}


def validate_status() -> dict:
    status = json.loads(STATUS_PATH.read_text(encoding="utf-8"))
    require(status["schema_version"] == "shieldSim.material_quantity_data_status.v1",
            "material/quantity audit has an unexpected schema")
    expected = {
        "Skin": "Skin_ICRP_G4_11_4",
        "EyeLens": "EyeLens_ICRP_G4_11_4",
        "BFO": "BFO_soft_tissue_proxy",
        "CNS": "Brain_ICRP_G4_11_4",
        "SoftTissue": "SoftTissue_ICRU4_G4_11_4",
        "Si": "Si_natural_G4_11_4",
        "SiO2": "SiO2_G4_11_4",
        "SiC": "SiC_6H_20C",
        "GaAs": "GaAs_G4_11_4",
        "InGaAs": "Ga0p47In0p53As_300K",
        "Ge": "Ge_natural_G4_11_4",
        "Water": "Water_G4_11_4",
    }
    quantities = {"TID", "DDD", "n_eq", "LET", "H100/10"}
    records = {record["label"]: record for record in status["materials"]}
    require(set(records) == set(expected), "status matrix does not cover exactly 12 absorbers")
    valid_statuses = set(status["status_definitions"])
    for label, material_id in expected.items():
        record = records[label]
        require(record["material_id"] == material_id,
                f"material identity mismatch for {label}")
        require(set(record["quantities"]) == quantities,
                f"quantity matrix is incomplete for {label}")
        require(set(record["quantities"].values()) <= valid_statuses,
                f"undefined status used for {label}")
        require(record["quantities"]["H100/10"] == "absorber_independent",
                f"H100/10 must not depend on absorber data: {label}")
    require(records["Si"]["quantities"]["n_eq"] == "production_enabled",
            "silicon n_eq production boundary changed")
    require(records["GaAs"]["quantities"]["n_eq"] ==
            "reference_data_prepared_pending_transport_validation",
            "GaAs target-specific n_eq candidate is not recorded")
    for name, relative in status["prepared_packages"].items():
        require((ROOT / relative).is_file(), f"prepared package is missing: {name}")
    return {"materials": len(records), "quantities_per_material": len(quantities)}


def main() -> int:
    nist_files = authenticate_nist()
    packages = {}
    status = json.loads(STATUS_PATH.read_text(encoding="utf-8"))
    for name, relative in status["prepared_packages"].items():
        if relative.endswith("candidate_manifest.v1.json"):
            packages[name] = authenticate_candidate(relative)
    matrix = validate_status()
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_id": "material_reference_data",
        "status": "passed",
        "scope": "offline integrity and coverage audit; no production promotion",
        "status_matrix_sha256": sha256(STATUS_PATH),
        "nist_authenticated_files": nist_files,
        "candidate_packages": packages,
        "matrix": matrix,
    }
    (EVIDENCE / "material_reference_data.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"Material reference-data checks passed: {nist_files} NIST files, "
        f"{len(packages)} candidate packages, {matrix['materials']} materials"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

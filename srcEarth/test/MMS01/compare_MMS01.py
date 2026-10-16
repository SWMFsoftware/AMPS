#!/usr/bin/env python3
"""Compare native AMPS numerical products with archived HPCA Level-2 flux."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from mms01_common import (
    ENERGIES, MAX_RELATIVE_ERROR, PROFILES, atomic_json, read_csv, sha256_file,
    utc_now, write_csv,
)


def _number(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{context}: invalid number {value!r}") from error
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{context}: expected a finite nonnegative value")
    return result


def _relative(actual: float, reference: float) -> float:
    """Compare physical zero exactly; never add an arbitrary logarithmic floor."""
    if reference == 0.0:
        return 0.0 if actual == 0.0 else math.inf
    return abs(actual - reference) / reference


def _key(row: dict[str, str]) -> tuple[int, int, int, int, int]:
    return tuple(int(row[name]) for name in
                 ("distribution_index", "sweep_index", "start_azimuth", "anode", "energy_index"))


def compare(output_root: Path) -> dict[str, Any]:
    """Validate provenance/completeness before computing unchanged numerical gates."""
    prepared_path = output_root / "prepared_manifest.json"
    native_path = output_root / "native_run_manifest.json"
    prepared = json.loads(prepared_path.read_text(encoding="utf-8"))
    native = json.loads(native_path.read_text(encoding="utf-8"))
    if prepared.get("schema") != "amps-mms01-prepared-v1":
        raise ValueError("unrecognized MMS01 preparation manifest")
    if native.get("schema") != "amps-mms01-native-run-v2":
        raise ValueError("unrecognized MMS01 native reconstruction manifest")
    if prepared.get("profile") not in PROFILES or native.get("profile") != prepared.get("profile"):
        raise ValueError("MMS01 preparation/native profile mismatch")
    if native.get("return_code") != 0:
        raise ValueError("native MMS01 calculation did not complete successfully")
    for artifact in prepared["artifacts"].values():
        if sha256_file(Path(artifact["path"])) != artifact["sha256"]:
            raise ValueError(f"prepared artifact changed: {artifact['path']}")
    for artifact in native["artifacts"].values():
        if sha256_file(Path(artifact["path"])) != artifact["sha256"]:
            raise ValueError(f"native artifact changed: {artifact['path']}")
    for role in ("production_header", "adapter_source"):
        item = native.get(role, {})
        path = Path(item.get("path", ""))
        if not path.is_file() or sha256_file(path) != item.get("sha256"):
            raise ValueError(f"native calculation code changed since execution: {role}")
    prepared_input = output_root / "input" / "hpca_psd_native_input.dat"
    if sha256_file(prepared_input) != native.get("prepared_input_sha256"):
        raise ValueError("native run does not match the current prepared PSD input")

    reference_directional = read_csv(output_root / "reference" / "directional_reference.csv")
    native_directional = read_csv(output_root / "native" / "directional_products.csv")
    if len(reference_directional) != prepared["directional_rows"] or len(native_directional) != len(reference_directional):
        raise ValueError("directional output is incomplete")
    directional_rows: list[dict[str, Any]] = []
    directional_error: list[float] = []
    zero_cells = 0
    for reference, reconstructed in zip(reference_directional, native_directional):
        if _key(reference) != _key(reconstructed):
            raise ValueError("native directional row identity/order differs from reference")
        ref_energy = _number(reference["energy_eV"], "reference directional energy")
        reconstructed_energy = _number(reconstructed["energy_eV"], "native directional energy")
        if ref_energy != reconstructed_energy:
            raise ValueError("native directional energy coordinate changed")
        ref = _number(reference["reference_differential_intensity"], "reference intensity")
        calculated = _number(
            reconstructed["kernel_reconstructed_differential_intensity"],
            "kernel-reconstructed intensity",
        )
        error = _relative(calculated, ref)
        directional_error.append(error)
        if ref == 0.0:
            zero_cells += 1
        directional_rows.append({
            **reference,
            "kernel_reconstructed_differential_intensity": f"{calculated:.17g}",
            "relative_error": f"{error:.17g}",
        })

    reference_spectrum = read_csv(output_root / "reference" / "spectrum_reference.csv")
    native_spectrum = read_csv(output_root / "native" / "spectrum_products.csv")
    if len(reference_spectrum) != prepared["spectrum_rows"] or len(native_spectrum) != len(reference_spectrum):
        raise ValueError("spectrum output is incomplete")
    spectrum_rows: list[dict[str, Any]] = []
    spectrum_error: list[float] = []
    for reference, reconstructed in zip(reference_spectrum, native_spectrum):
        ref_key = (int(reference["distribution_index"]), int(reference["energy_index"]))
        reconstructed_key = (int(reconstructed["distribution_index"]),
                             int(reconstructed["energy_index"]))
        if ref_key != reconstructed_key:
            raise ValueError("native spectrum row identity/order differs from reference")
        ref_energy = _number(reference["energy_eV"], "reference spectrum energy")
        reconstructed_energy = _number(reconstructed["energy_eV"], "native spectrum energy")
        if ref_energy != reconstructed_energy:
            raise ValueError("native spectrum energy coordinate changed")
        ref = _number(reference["reference_angular_integral"], "reference spectrum")
        calculated = _number(
            reconstructed["kernel_reconstructed_angular_integral"],
            "kernel-reconstructed spectrum",
        )
        error = _relative(calculated, ref)
        spectrum_error.append(error)
        spectrum_rows.append({
            **reference,
            "kernel_reconstructed_angular_integral": f"{calculated:.17g}",
            "relative_error": f"{error:.17g}",
        })

    comparison_dir = output_root / "comparison"
    write_csv(comparison_dir / "directional_comparison.csv", directional_rows[0].keys(), directional_rows)
    write_csv(comparison_dir / "spectrum_comparison.csv", spectrum_rows[0].keys(), spectrum_rows)
    max_directional = max(directional_error)
    max_spectrum = max(spectrum_error)
    gates = {
        "directional_complete": len(directional_rows) == prepared["directional_rows"],
        "spectrum_complete": len(spectrum_rows) == prepared["selected_distributions"] * ENERGIES,
        "directional_max_relative_error": max_directional <= MAX_RELATIVE_ERROR,
        "spectrum_max_relative_error": max_spectrum <= MAX_RELATIVE_ERROR,
    }
    summary = {
        "schema": "amps-mms01-comparison-v2",
        "created_utc": utc_now(),
        "profile": prepared["profile"],
        "status": "PASS" if all(gates.values()) else "FAIL",
        "gate_limit_relative_error": MAX_RELATIVE_ERROR,
        "gates": gates,
        "directional": {
            "rows": len(directional_rows), "zero_reference_cells": zero_cells,
            "maximum_relative_error": max_directional,
        },
        "spectrum": {"rows": len(spectrum_rows), "maximum_relative_error": max_spectrum},
        "interpretation": (
            "PASS validates the production PSD-to-differential-intensity and calibrated "
            "angular-integration kernels by reconstructing a same-instrument HPCA Level-2 product. "
            "It is not transport, density, energy-integrated-flux, or independent absolute-"
            "calibration validation, and it is not a simulation-versus-observation comparison."
        ),
        "artifacts": {
            path.name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
            for path in (comparison_dir / "directional_comparison.csv",
                         comparison_dir / "spectrum_comparison.csv")
        },
    }
    atomic_json(output_root / "comparison_summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    summary = compare(args.output_root.resolve())
    print(f"MMS01 comparison: {summary['status']}")
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare authenticated AMPS spectra with qualified GOES-16/SGPS P8A data."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from goes01_common import (
    GATES, M2_TO_CM2, atomic_json, power_law_bin_average, ratio_metrics, read_csv,
    read_spectrum_zones, sha256_file, utc_now, write_csv,
)


PATTERNS = {
    "spectrum": "mode3d_points_spectrum_snapshot_*.dat",
    "density": "mode3d_points_density_snapshot_*.dat",
    "flux": "mode3d_points_flux_snapshot_*.dat",
    "termination": "mode3d_termination_summary_snapshot_*.dat",
    "manifest": "standalone_mode3d_manifest_snapshot_*.json",
}


def _validated_products(branch: Path, prepared: dict[str, Any]) -> dict[str, list[Path]]:
    """Require a complete, unchanged native artifact set for every prepared epoch."""
    status_path = branch / "run_status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"missing physical-run evidence {status_path}")
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "PASS" or status.get("return_code") != 0:
        raise ValueError("GOES01 physical execution did not pass")
    for key in ("input_sha256", "driver_sha256", "spectrum_sha256",
                "snapshot_list_sha256", "observations_sha256"):
        if status.get(key) != prepared.get(key):
            raise ValueError(f"run/preparation provenance mismatch for {key}")
    count = int(prepared["trajectory"]["kept"])
    fingerprints = status.get("output_sha256", {})
    result: dict[str, list[Path]] = {}
    for role, pattern in PATTERNS.items():
        paths = sorted(branch.glob(pattern))
        if len(paths) != count:
            raise FileNotFoundError(
                f"{branch}: expected {count} {role} snapshots, found {len(paths)}"
            )
        for path in paths:
            if fingerprints.get(path.name) != sha256_file(path):
                raise ValueError(f"native artifact changed or is unverified: {path}")
        result[role] = paths
    return result


def _metric_gates(prefix: str, metrics: dict[str, Any], factor_min: float) -> list[dict[str, Any]]:
    median = metrics["median_abs_log10_ratio"]
    return [
        {"name": f"{prefix}.median_abs_log10_ratio", "operator": "<=",
         "threshold": GATES["median_abs_log10_ratio_max"], "value": median,
         "passed": median is not None and median <= GATES["median_abs_log10_ratio_max"]},
        {"name": f"{prefix}.factor_two_fraction", "operator": ">=",
         "threshold": factor_min, "value": metrics["factor_two_fraction"],
         "passed": metrics["factor_two_fraction"] >= factor_min},
    ]


def compare(output_root: Path) -> dict[str, Any]:
    """Fold native spectra into each documented head passband and score all rows."""
    branch = output_root / "GOES16"
    prepared_path = branch / "prepared_run.json"
    prepared = json.loads(prepared_path.read_text(encoding="utf-8"))
    for path_key, hash_key in (
        ("input_file", "input_sha256"), ("driver_file", "driver_sha256"),
        ("spectrum_file", "spectrum_sha256"),
        ("snapshot_list_file", "snapshot_list_sha256"),
        ("observations_file", "observations_sha256"),
        ("selection_ledger_file", "selection_ledger_sha256"),
    ):
        if sha256_file(Path(prepared[path_key])) != prepared[hash_key]:
            raise ValueError(f"prepared input changed after selection: {prepared[path_key]}")
    products = _validated_products(branch, prepared)
    zones = []
    for path in products["spectrum"]:
        # Unit conversion below is valid only for this explicit native contract.
        # Do not infer semantics from a familiar column name after an output-format
        # change; require the producer's per-file metadata at every snapshot.
        content = path.read_text(encoding="utf-8")
        for key, expected in (
            ("SPECTRUM_INTENSITY_UNIT", "m^-2 s^-1 sr^-1 MeV^-1"),
            ("SPECTRUM_ENERGY_BASIS", "PER_PARTICLE"),
            ("STEP6_CHARACTERISTIC_MAPPING", "STATIC_MAGNETIC"),
        ):
            values = re.findall(r'AUXDATA\s+' + key + r'="([^"]+)"', content)
            if not values or any(value != expected for value in values):
                raise ValueError(f"{path}: missing or inconsistent native {key}")
        zones.extend(read_spectrum_zones(path))
    zone_by_epoch = {zone.utc.rstrip("Z") + "Z": zone for zone in zones}
    if len(zone_by_epoch) != len(zones):
        raise ValueError("duplicate AMPS spectrum epoch")

    observations = read_csv(Path(prepared["observations_file"]))
    expected_epochs = {row["comparison_midpoint_utc"] for row in observations}
    if set(zone_by_epoch) != expected_epochs:
        raise ValueError("AMPS spectrum epochs differ from the prepared observation population")
    rows: list[dict[str, Any]] = []
    maximum_unresolved = 0.0
    for observed in observations:
        epoch = observed["comparison_midpoint_utc"]
        zone = zone_by_epoch[epoch]
        columns = {name: index for index, name in enumerate(zone.columns)}
        required = ("E_MeV", "J_local_perMeV", "unresolved_fraction", "N_sampled",
                    "N_resolved")
        missing = [name for name in required if name not in columns]
        if missing:
            raise ValueError(f"{epoch}: native spectrum lacks columns {missing}")
        energy = zone.rows[:, columns["E_MeV"]]
        intensity_si = zone.rows[:, columns["J_local_perMeV"]]
        unresolved = zone.rows[:, columns["unresolved_fraction"]]
        sampled = zone.rows[:, columns["N_sampled"]]
        resolved = zone.rows[:, columns["N_resolved"]]
        requested = prepared["profile_parameters"]["requested_directions_per_energy"]
        expected_nodes = prepared["profile_parameters"]["spectrum_nodes"]
        if len(energy) != expected_nodes or np.any(sampled != requested):
            raise ValueError(f"{epoch}: effective quadrature differs from preparation")
        if np.any(resolved < 0) or np.any(resolved > sampled) or not np.allclose(
                unresolved, (sampled - resolved) / sampled, atol=1.0e-12, rtol=0.0):
            raise ValueError(f"{epoch}: inconsistent trajectory-resolution accounting")
        epoch_unresolved = float(np.max(unresolved))
        maximum_unresolved = max(maximum_unresolved, epoch_unresolved)
        low, high = float(observed["energy_min_mev"]), float(observed["energy_max_mev"])
        model = power_law_bin_average(energy, intensity_si, low, high) * M2_TO_CM2
        reference = float(observed["observed_intensity_cm2_s_sr_mev"])
        uncertainty = float(observed["observed_uncertainty_cm2_s_sr_mev"])
        if reference <= 0.0 or uncertainty <= 0.0 or model < 0.0:
            raise ValueError(f"{epoch}/{observed['head']}: invalid comparison value")
        rows.append({
            "comparison_midpoint_utc": epoch,
            "record_start_utc": observed["record_start_utc"],
            "head": observed["head"], "energy_min_mev": low, "energy_max_mev": high,
            "model": model, "reference": reference,
            "reference_uncertainty": uncertainty,
            "model_observed_ratio": model / reference,
            # This is only a scale ratio to the delivered uncertainty magnitude;
            # no confidence level or Gaussian interpretation is assigned to it.
            "uncertainty_scaled_residual": (model - reference) / uncertainty,
            "two_head_envelope_separation": float(observed["two_head_envelope_separation"]),
            "gse_local_time_hours": float(observed["gse_local_time_hours"]),
            "max_unresolved_fraction_at_epoch": epoch_unresolved,
            "physical_zero_model": int(model == 0.0),
        })
    comparison_path = output_root / "comparison" / "p8a_intensity_comparison.csv"
    write_csv(comparison_path, rows)

    heads = sorted({row["head"] for row in rows})
    per_head = {head: ratio_metrics([row for row in rows if row["head"] == head])
                for head in heads}
    pooled = ratio_metrics(rows)
    epoch_count = len(expected_epochs)
    gates: list[dict[str, Any]] = [{
        "name": "coverage.matched_epochs", "operator": ">=",
        "threshold": GATES["minimum_matched_epochs"], "value": epoch_count,
        "passed": epoch_count >= GATES["minimum_matched_epochs"],
    }, {
        "name": "numerics.maximum_unresolved_fraction", "operator": "<=",
        "threshold": GATES["maximum_unresolved_fraction"], "value": maximum_unresolved,
        "passed": maximum_unresolved <= GATES["maximum_unresolved_fraction"],
    }]
    for head in heads:
        gates.extend(_metric_gates(head, per_head[head],
                                   GATES["factor_two_fraction_each_head_min"]))
    gates.extend(_metric_gates("pooled", pooled,
                               GATES["factor_two_fraction_pooled_min"]))
    result = {
        "schema": "amps-goes01-comparison-v1", "created_utc": utc_now(),
        "status": "PASS" if all(gate["passed"] for gate in gates) else "FAIL",
        "profile": prepared["profile"], "matched_epochs": epoch_count,
        "comparison_rows": len(rows), "maximum_unresolved_fraction": maximum_unresolved,
        "metrics": {"pooled": pooled, "by_head": per_head}, "gates": gates,
        "comparison_table": str(comparison_path.resolve()),
        "comparison_table_sha256": sha256_file(comparison_path),
        "interpretation": (
            "PASS supports near-Earth high-energy P8A intensity transport and absolute "
            "normalization for the prequalified population only; it is not a density, "
            "directional-response, SEP-source, heliospheric-propagation, or onset test."
        ),
    }
    atomic_json(output_root / "comparison_summary.json", result)
    return result


def compare_convergence(coarse_root: Path, fine_root: Path) -> dict[str, Any]:
    """Compare exact common native-model rows from QUICK and COMPLETE campaigns.

    Observation values play no role in this check.  Re-running ``compare`` on both
    roots first revalidates their immutable inputs and hashed native artifacts. The
    QUICK population must be a complete subset of COMPLETE; otherwise cadence loss
    could masquerade as convergence by silently dropping a difficult epoch.
    """
    coarse_root, fine_root = coarse_root.resolve(), fine_root.resolve()
    coarse_summary, fine_summary = compare(coarse_root), compare(fine_root)
    coarse_prepared = json.loads(
        (coarse_root / "GOES16" / "prepared_run.json").read_text(encoding="utf-8"))
    fine_prepared = json.loads(
        (fine_root / "GOES16" / "prepared_run.json").read_text(encoding="utf-8"))
    if coarse_prepared["profile"] != "QUICK" or fine_prepared["profile"] != "COMPLETE":
        raise ValueError("convergence requires QUICK coarse evidence and COMPLETE fine evidence")
    coarse_profile, fine_profile = (coarse_prepared["profile_parameters"],
                                    fine_prepared["profile_parameters"])
    if not (fine_profile["cadence_minutes"] < coarse_profile["cadence_minutes"] and
            fine_profile["ds_intervals"] > coarse_profile["ds_intervals"] and
            fine_profile["requested_directions_per_energy"] >
            coarse_profile["requested_directions_per_energy"]):
        raise ValueError("COMPLETE is not a strict cadence/energy/angular refinement of QUICK")
    for key in ("driver_sha256", "spectrum_sha256"):
        if coarse_prepared[key] != fine_prepared[key]:
            raise ValueError(f"convergence source mismatch: {key}")
    for key in ("source", "science_contract", "numerical_controls", "mode3d_mesh"):
        if coarse_prepared[key] != fine_prepared[key]:
            raise ValueError(f"convergence campaigns differ in non-refinement contract: {key}")
    coarse_rows = read_csv(Path(coarse_summary["comparison_table"]))
    fine_rows = read_csv(Path(fine_summary["comparison_table"]))

    def index(rows: list[dict[str, str]], label: str) -> dict[tuple[str, str], dict[str, str]]:
        result: dict[tuple[str, str], dict[str, str]] = {}
        for row in rows:
            identity = (row["comparison_midpoint_utc"], row["head"])
            if identity in result:
                raise ValueError(f"duplicate {label} convergence identity {identity}")
            result[identity] = row
        return result

    coarse, fine = index(coarse_rows, "coarse"), index(fine_rows, "fine")
    if not set(coarse).issubset(fine) or len(coarse) != 2 * coarse_summary["matched_epochs"]:
        raise ValueError("COMPLETE does not contain the full QUICK per-head population")
    evidence: list[dict[str, Any]] = []
    for identity in sorted(coarse):
        coarse_model, fine_model = float(coarse[identity]["model"]), float(fine[identity]["model"])
        scale = max(abs(coarse_model), abs(fine_model))
        change = (0.0 if scale <= 1.0e-30 and abs(coarse_model - fine_model) <= 1.0e-30
                  else abs(coarse_model - fine_model) / max(scale, 1.0e-30))
        evidence.append({
            "comparison_midpoint_utc": identity[0], "head": identity[1],
            "quick_model": coarse_model, "complete_model": fine_model,
            "symmetric_fractional_change": change,
        })
    table = fine_root / "comparison" / "quick_complete_convergence.csv"
    write_csv(table, evidence)
    maximum = max(row["symmetric_fractional_change"] for row in evidence)
    threshold = 0.02
    result = {
        "schema": "amps-goes01-convergence-v1", "created_utc": utc_now(),
        "status": "PASS" if maximum <= threshold else "FAIL",
        "coarse_root": str(coarse_root), "fine_root": str(fine_root),
        "common_rows": len(evidence), "maximum_symmetric_fractional_change": maximum,
        "threshold": threshold, "comparison_table": str(table.resolve()),
        "comparison_table_sha256": sha256_file(table),
        "coarse_run_status_sha256": sha256_file(coarse_root / "GOES16" / "run_status.json"),
        "fine_run_status_sha256": sha256_file(fine_root / "GOES16" / "run_status.json"),
    }
    atomic_json(fine_root / "convergence_summary.json", result)
    return result

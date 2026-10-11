#!/usr/bin/env python3
"""Compare saved SA01 AMPS products with immutable RBSP/REPT observations."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from sa01_common import (
    CHANNELS, GATES, M2_TO_CM2, SPACECRAFT, atomic_json, bin_average,
    fractional_change, ratio_metrics, read_csv, read_point_table,
    read_spectrum_zones, resolve_data_root, spectrum_integral, utc_now,
    validate_data_package, validate_spacecraft, write_csv, sha256_file,
)


PRODUCTS = ("spectrum", "channel_flux", "total_flux", "density")


def _prepared_solver(branch: Path) -> tuple[str, int]:
    prepared_path = branch / "prepared_run.json"
    if not prepared_path.is_file():
        raise FileNotFoundError(f"missing preparation provenance {prepared_path}")
    prepared = json.loads(prepared_path.read_text(encoding="utf-8"))
    solver = str(prepared.get("solver", "GRIDLESS")).upper()
    return solver, int(prepared["trajectory"]["kept"])


def _branch_product_paths(branch: Path, product: str) -> list[Path]:
    """Resolve one logical product without hiding incomplete Mode3D snapshots."""

    solver, expected_count = _prepared_solver(branch)
    if solver == "GRIDLESS":
        names = {
            "spectrum": "gridless_points_spectrum.dat",
            "flux": "gridless_points_flux.dat",
            "density": "gridless_points_density.dat",
            "termination": "gridless_termination_summary.dat",
            "manifest": "standalone_run_manifest.json",
        }
        path = branch / names[product]
        if not path.is_file():
            raise FileNotFoundError(f"missing required physical output {path}")
        return [path]
    if solver == "GRIDDED":
        patterns = {
            "spectrum": "mode3d_points_spectrum_snapshot_*.dat",
            "flux": "mode3d_points_flux_snapshot_*.dat",
            "density": "mode3d_points_density_snapshot_*.dat",
            "termination": "mode3d_termination_summary_snapshot_*.dat",
            "manifest": "standalone_mode3d_manifest_snapshot_*.json",
        }
        paths = sorted(branch.glob(patterns[product]))
        if len(paths) != expected_count:
            raise FileNotFoundError(
                f"{branch}: expected {expected_count} GRIDDED {product} snapshots, "
                f"found {len(paths)} matching {patterns[product]}"
            )
        return paths
    raise ValueError(f"{branch}: unsupported prepared solver {solver!r}")


def _read_all_spectrum_zones(paths: list[Path]):
    zones = []
    for path in paths:
        zones.extend(read_spectrum_zones(path))
    return zones


def _read_all_point_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        rows.extend(read_point_table(path))
    return rows


def validate_branch_evidence(branch: Path) -> set[str]:
    """Bind comparison coverage to the prepared ephemeris and completed run.

    Matching the three product tables to each other is insufficient: all three
    could omit the same failed epoch. Hash checks also prevent a previous PASS
    record from authorizing subsequently replaced or truncated physical output.
    """
    prepared = json.loads((branch / "prepared_run.json").read_text())
    status = json.loads((branch / "run_status.json").read_text())
    if status.get("status") != "PASS" or status.get("return_code") != 0:
        raise ValueError(f"{branch}: physical run did not complete successfully")
    for key in ("solver", "trajectory", "profile_parameters", "numerical_controls",
                "mode3d_mesh", "trajectory_closure", "science_contract"):
        if status.get(key) != prepared.get(key):
            raise ValueError(f"{branch}: run/preparation {key} mismatch")
    for path_key, hash_key in (("input_file", "input_sha256"),
                               ("driver_file", "driver_sha256"),
                               ("spectrum_file", "spectrum_sha256")):
        if sha256_file(Path(prepared[path_key])) != prepared[hash_key]:
            raise ValueError(f"{branch}: changed prepared {path_key}")
        if status.get(hash_key) != prepared[hash_key]:
            raise ValueError(f"{branch}: run/preparation {hash_key} mismatch")
    trajectory = prepared["trajectory"]
    path = Path(trajectory["output"])
    if sha256_file(path) != trajectory["sha256"]:
        raise ValueError(f"{branch}: prepared trajectory changed")
    epochs = [line.split()[0].rstrip("Z") + "Z" for line in path.read_text().splitlines()
              if line.strip() and not line.lstrip().startswith(("#", "!"))]
    if len(epochs) != len(set(epochs)) or len(epochs) != trajectory["kept"]:
        raise ValueError(f"{branch}: duplicate or inconsistent prepared epochs")
    fingerprints = status.get("output_sha256", {})
    for role in ("spectrum", "flux", "density", "termination", "manifest"):
        for artifact in _branch_product_paths(branch, role):
            if fingerprints.get(artifact.name) != sha256_file(artifact):
                raise ValueError(f"{branch}: unverified or changed output {artifact.name}")
    return set(epochs)


def _number(value: str, context: str, *, nonnegative: bool = True) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{context}: invalid number {value!r}") from error
    if not math.isfinite(result) or (nonnegative and result < 0.0):
        raise ValueError(f"{context}: expected finite {'nonnegative ' if nonnegative else ''}value")
    return result


def _index_unique(rows: list[dict[str, str]], key: str, path: Path) -> dict[str, dict[str, str]]:
    """Reject duplicate identities before dictionary construction could hide rows."""
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row[key]
        if value in result:
            raise ValueError(f"{path}: duplicate {key}={value}")
        result[value] = row
    return result


def compare_branch(validation_root: Path, output_root: Path,
                   spacecraft: str) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    """Validate complete native evidence, then fold spectra and compare all products.

    Native integrated products remain the modeled observables. Independent spectral
    reconstruction is a diagnostic, never a replacement for a failing native value.
    """
    branch = output_root / spacecraft
    expected_epochs = validate_branch_evidence(branch)
    spectrum_paths = _branch_product_paths(branch, "spectrum")
    flux_paths = _branch_product_paths(branch, "flux")
    density_paths = _branch_product_paths(branch, "density")
    _branch_product_paths(branch, "termination")
    _branch_product_paths(branch, "manifest")

    zones = _read_all_spectrum_zones(spectrum_paths)
    zone_by_time = {zone.utc: zone for zone in zones}
    if len(zone_by_time) != len(zones):
        raise ValueError(f"{spacecraft}: duplicate spectrum epoch across physical products")
    flux_label = flux_paths[0]
    density_label = density_paths[0]
    flux_by_time = _index_unique(_read_all_point_rows(flux_paths), "TimeUTC", flux_label)
    density_by_time = _index_unique(_read_all_point_rows(density_paths), "TimeUTC", density_label)
    if set(zone_by_time) != set(flux_by_time) or set(zone_by_time) != set(density_by_time):
        raise ValueError(f"{spacecraft}: spectrum/flux/density epochs differ")
    if set(zone_by_time) != expected_epochs:
        raise ValueError(f"{spacecraft}: physical products do not cover exactly the prepared epochs")
    for path in spectrum_paths:
        content = path.read_text()
        for key, expected in (("SPECTRUM_INTENSITY_UNIT", "m^-2 s^-1 sr^-1 MeV^-1"),
                              ("SPECTRUM_ENERGY_BASIS", "PER_PARTICLE"),
                              ("STEP6_CHARACTERISTIC_MAPPING", "STATIC_MAGNETIC")):
            values = re.findall(r'AUXDATA\s+' + key + r'="([^"]+)"', content)
            if not values or any(value != expected for value in values):
                raise ValueError(f"{path}: missing or inconsistent declared {key}")
    profile = json.loads((branch / "prepared_run.json").read_text())["profile_parameters"]

    processed = validation_root / "processed" / "SA01"
    observed_spectrum = read_csv(processed / f"{spacecraft.lower()}_rept_spectrum_5min.csv")
    observed_density_path = processed / f"{spacecraft.lower()}_rept_density_5min.csv"
    observed_density = _index_unique(read_csv(observed_density_path), "utc", observed_density_path)
    obs_by_time: dict[str, dict[int, dict[str, str]]] = defaultdict(dict)
    for row in observed_spectrum:
        bounds = (float(row["energy_min_mev"]), float(row["energy_max_mev"]))
        matches = [index for index, (_, low, high) in enumerate(CHANNELS)
                   if abs(bounds[0] - low) < 1e-8 and abs(bounds[1] - high) < 1e-8]
        if len(matches) != 1:
            raise ValueError(f"{spacecraft}: unknown REPT interval {bounds}")
        index = matches[0]
        if index in obs_by_time[row["utc"]]:
            raise ValueError(f"{spacecraft}: duplicate observation time/channel")
        obs_by_time[row["utc"]][index] = row

    products: dict[str, list[dict[str, Any]]] = {name: [] for name in PRODUCTS}
    reconstruction: list[dict[str, Any]] = []
    maximum_unresolved = 0.0
    for utc, zone in sorted(zone_by_time.items()):
        if utc not in obs_by_time or len(obs_by_time[utc]) != len(CHANNELS):
            raise ValueError(f"{spacecraft}/{utc}: modeled epoch lacks all seven observations")
        columns = {name: index for index, name in enumerate(zone.columns)}
        required = ("E_MeV", "J_local_perMeV", "unresolved_fraction", "N_sampled", "N_resolved")
        missing = [name for name in required if name not in columns]
        if missing:
            raise ValueError(f"{spectrum_paths[0]}: missing columns {missing}")
        energy = zone.rows[:, columns["E_MeV"]]
        local_si = zone.rows[:, columns["J_local_perMeV"]]
        unresolved = zone.rows[:, columns["unresolved_fraction"]]
        sampled = zone.rows[:, columns["N_sampled"]]
        resolved = zone.rows[:, columns["N_resolved"]]
        if len(energy) != profile["spectrum_nodes"] or np.any(sampled != profile["requested_directions_per_energy"]):
            raise ValueError(f"{spacecraft}/{utc}: effective energy/direction grid differs from preparation")
        if np.any(resolved < 0) or np.any(resolved > sampled) or np.any(resolved != np.floor(resolved)) or not np.allclose(
                unresolved, (sampled - resolved) / sampled, atol=1e-12, rtol=0):
            raise ValueError(f"{spacecraft}/{utc}: inconsistent trajectory-resolution accounting")
        if np.any(unresolved < 0.0) or np.any(unresolved > 1.0):
            raise ValueError(f"{spacecraft}/{utc}: invalid unresolved fraction")
        epoch_unresolved = float(np.max(unresolved))
        maximum_unresolved = max(maximum_unresolved, epoch_unresolved)

        reconstructed_channels_si: list[float] = []
        observed_channels_cm: list[float] = []
        context = obs_by_time[utc][0]
        for index, (name, low, high) in enumerate(CHANNELS):
            observed = obs_by_time[utc][index]
            reference = _number(observed["differential_flux_cm2_s_sr_mev"], "REPT intensity")
            if reference <= 0.0:
                raise ValueError(f"{spacecraft}/{utc}/{name}: nonpositive reference")
            folded_si = bin_average(energy, local_si, low, high)
            folded_cm = folded_si * M2_TO_CM2
            common = {
                "spacecraft": spacecraft, "utc": utc, "channel": name,
                "energy_min_mev": low, "energy_max_mev": high,
                "l_shell": float(observed["l_shell"]), "mlt_hour": float(observed["mlt_hour"]),
                "model": folded_cm, "reference": reference,
                "model_observed_ratio": folded_cm / reference,
                "physical_zero_model": int(folded_cm == 0.0),
                "max_unresolved_fraction_at_epoch": epoch_unresolved,
            }
            products["spectrum"].append(common)
            observed_integral = 4.0 * math.pi * reference * (high - low)
            native_key = f"F_{name}_m2s1"
            if native_key not in flux_by_time[utc]:
                raise ValueError(f"{flux_label}: missing {native_key}")
            native_integral = _number(flux_by_time[utc][native_key], native_key) * M2_TO_CM2
            products["channel_flux"].append({
                **common, "model": native_integral, "reference": observed_integral,
                "model_observed_ratio": native_integral / observed_integral,
                "physical_zero_model": int(native_integral == 0.0),
            })
            reconstructed_channels_si.append(4.0 * math.pi * folded_si * (high - low))
            observed_channels_cm.append(reference)

        reconstructed_total_si = 4.0 * math.pi * spectrum_integral(energy, local_si, 18.5, 116.0)
        reconstructed_density_cm3 = (
            4.0 * math.pi * spectrum_integral(
                energy, local_si * M2_TO_CM2, 18.5, 116.0, speed_weighted=True
            )
        )
        native_total_si = _number(flux_by_time[utc]["F_tot_m2s1"], "F_tot_m2s1")
        native_density_cm3 = _number(density_by_time[utc]["N_cm^-3"], "N_cm^-3")
        if utc not in observed_density:
            raise ValueError(f"{spacecraft}/{utc}: missing complete observed density row")
        density_obs = observed_density[utc]
        observed_total = _number(
            density_obs["isotropic_omnidirectional_flux_cm2_s"], "observed total flux"
        )
        observed_n = _number(
            density_obs["isotropic_equivalent_density_cm3"], "observed density"
        )
        base = {
            "spacecraft": spacecraft, "utc": utc,
            "l_shell": float(context["l_shell"]), "mlt_hour": float(context["mlt_hour"]),
            "max_unresolved_fraction_at_epoch": epoch_unresolved,
        }
        products["total_flux"].append({
            **base, "model": native_total_si * M2_TO_CM2, "reference": observed_total,
            "model_observed_ratio": native_total_si * M2_TO_CM2 / observed_total,
            "physical_zero_model": int(native_total_si == 0.0),
        })
        products["density"].append({
            **base, "model": native_density_cm3, "reference": observed_n,
            "model_observed_ratio": native_density_cm3 / observed_n,
            "physical_zero_model": int(native_density_cm3 == 0.0),
            "reference_interpretation": "isotropic-equivalent spectrum-derived density",
        })
        reconstruction.append({
            "spacecraft": spacecraft, "utc": utc,
            "native_total_flux_m2_s": native_total_si,
            "reconstructed_total_flux_m2_s": reconstructed_total_si,
            "total_flux_fractional_residual": fractional_change(native_total_si, reconstructed_total_si),
            "native_density_cm3": native_density_cm3,
            "reconstructed_density_cm3": reconstructed_density_cm3,
            "density_fractional_residual": fractional_change(native_density_cm3, reconstructed_density_cm3),
            "maximum_channel_fractional_residual": max(
                fractional_change(
                    _number(flux_by_time[utc][f"F_{name}_m2s1"], name), reconstructed
                ) for (name, _, _), reconstructed in zip(CHANNELS, reconstructed_channels_si)
            ),
        })

    comparison_dir = branch / "comparison"
    for name, rows in products.items():
        write_csv(comparison_dir / f"{name}_comparison.csv", rows)
    write_csv(comparison_dir / "native_reconstruction.csv", reconstruction)

    metrics: dict[str, Any] = {
        "spacecraft": spacecraft, "modeled_epochs": len(zone_by_time),
        "matched_epochs_by_product": {
            name: len({row["utc"] for row in rows}) for name, rows in products.items()
        },
        "max_unresolved_fraction": maximum_unresolved,
        "products": {name: ratio_metrics(rows) for name, rows in products.items()},
        "by_channel": {
            product: {
                name: ratio_metrics([row for row in products[product] if row["channel"] == name])
                for name, _, _ in CHANNELS
            } for product in ("spectrum", "channel_flux")
        },
        "reconstruction": {
            "max_total_flux_fractional_residual": max(row["total_flux_fractional_residual"] for row in reconstruction),
            "max_density_fractional_residual": max(row["density_fractional_residual"] for row in reconstruction),
            "max_channel_fractional_residual": max(row["maximum_channel_fractional_residual"] for row in reconstruction),
            "note": "Residuals compare native AMPS quadrature with independent zero-preserving spectral reconstruction.",
        },
    }
    return metrics, products


def _metric_gates(prefix: str, metric: dict[str, Any], fraction_min: float) -> list[dict[str, Any]]:
    median = metric["median_abs_log10_ratio"]
    return [
        {"name": prefix + ".median_abs_log10_ratio", "value": median, "operator": "<=",
         "threshold": GATES["median_abs_log10_ratio_max"],
         "passed": median is not None and median <= GATES["median_abs_log10_ratio_max"]},
        {"name": prefix + ".factor_two_fraction", "value": metric["factor_two_fraction"],
         "operator": ">=", "threshold": fraction_min,
         "passed": metric["factor_two_fraction"] >= fraction_min},
    ]


def compare_all(data_path: Path, output_root: Path,
                spacecraft: tuple[str, ...]) -> dict[str, Any]:
    _, validation_root = resolve_data_root(data_path)
    validate_data_package(validation_root)
    output = output_root.expanduser().resolve()
    spacecraft_metrics: dict[str, Any] = {}
    all_products: dict[str, list[dict[str, Any]]] = {name: [] for name in PRODUCTS}
    gates: list[dict[str, Any]] = []
    for item in spacecraft:
        metrics, products = compare_branch(validation_root, output, item)
        spacecraft_metrics[item] = metrics
        for name in PRODUCTS:
            all_products[name].extend(products[name])
            gates.extend(_metric_gates(
                f"{item}.{name}", metrics["products"][name],
                GATES["factor_two_fraction_each_min"],
            ))
            count = metrics["matched_epochs_by_product"][name]
            gates.append({
                "name": f"{item}.{name}.matched_epochs", "value": count,
                "operator": ">=", "threshold": GATES["min_matched_epochs"],
                "passed": count >= GATES["min_matched_epochs"],
            })
        gates.append({
            "name": f"{item}.max_unresolved_fraction",
            "value": metrics["max_unresolved_fraction"], "operator": "<=",
            "threshold": GATES["max_unresolved_fraction"],
            "passed": metrics["max_unresolved_fraction"] <= GATES["max_unresolved_fraction"],
        })

    pooled = {name: ratio_metrics(rows) for name, rows in all_products.items()}
    for name in PRODUCTS:
        gates.extend(_metric_gates(
            f"pooled.{name}", pooled[name], GATES["factor_two_fraction_pooled_min"]
        ))
    result = {
        "schema": "amps-sa01-observation-comparison-v2", "created_utc": utc_now(),
        # Comparison-only mode may run after source changes or on another checkout.
        # Record its own implementation identity, not merely the launcher's old one.
        "analysis_source_sha256": {name: sha256_file(Path(__file__).with_name(name))
                                   for name in ("compare_SA01.py", "sa01_common.py")},
        "status": "PASS" if all(gate["passed"] for gate in gates) else "FAIL",
        "spacecraft": spacecraft_metrics, "pooled": pooled, "gates": gates,
        "gate_policy": GATES,
        "notes": [
            "Physical-zero model predictions remain in every denominator and fail factor-two agreement.",
            "REPT density is isotropic-equivalent and spectrum-derived, not an independent direct density measurement.",
            "REPT rectangular channel folding is an approximation, not a calibrated angular/energy instrument simulator.",
            "Each modeled epoch is instantaneous at the selected trajectory time; observations are five-minute means.",
        ],
    }
    atomic_json(output / "comparison_summary.json", result)
    return result


def compare_convergence(coarse_root: Path, fine_root: Path,
                        spacecraft: tuple[str, ...], destination: Path) -> dict[str, Any]:
    """Evaluate all nine integrated products on a provenance-bound common population.

    Coverage, physical inputs and independent numerical settings are prerequisites.
    A cadence-only change or incomplete product set cannot establish convergence.
    """
    rows: list[dict[str, Any]] = []
    evidence: dict[str, Any] = {}
    for item in spacecraft:
        coarse_epochs = validate_branch_evidence(coarse_root / item)
        fine_epochs = validate_branch_evidence(fine_root / item)
        coarse_prepared = json.loads((coarse_root / item / "prepared_run.json").read_text())
        fine_prepared = json.loads((fine_root / item / "prepared_run.json").read_text())
        for key in ("driver_sha256", "spectrum_sha256", "science_contract"):
            if coarse_prepared.get(key) != fine_prepared.get(key):
                raise ValueError(f"{item}: refinement changed physical input {key}")
        def numerical_controls(prepared: dict[str, Any]) -> dict[str, Any]:
            profile = prepared["profile_parameters"]
            return {"energy_nodes": profile["spectrum_nodes"],
                    "angular_samples": profile["requested_directions_per_energy"],
                    "angular_grid": [profile.get("direction_nmu", 24), profile.get("direction_nphi", 48)],
                    "solver": prepared.get("solver"),
                    "mesh": prepared.get("mode3d_mesh"),
                    "trajectory_closure": prepared.get("trajectory_closure"),
                    "numerics": {key: value for key, value in prepared.get("numerical_controls", {}).items()
                                 if key != "cadence_minutes"}}
        if numerical_controls(coarse_prepared) == numerical_controls(fine_prepared):
            raise ValueError(f"{item}: identical numerical controls; cadence alone is not refinement")
        evidence[item] = {
            "coarse_controls": numerical_controls(coarse_prepared),
            "fine_controls": numerical_controls(fine_prepared),
            "common_epochs": sorted(coarse_epochs & fine_epochs),
            "coarse_only_epochs": sorted(coarse_epochs - fine_epochs),
            "fine_only_epochs": sorted(fine_epochs - coarse_epochs),
            "coarse_status_sha256": sha256_file(coarse_root / item / "run_status.json"),
            "fine_status_sha256": sha256_file(fine_root / item / "run_status.json"),
        }
        for product, filename, columns in (
            ("total_flux", "gridless_points_flux.dat", ["F_tot_m2s1"]),
            ("channel_flux", "gridless_points_flux.dat", [f"F_{name}_m2s1" for name, _, _ in CHANNELS]),
            ("density", "gridless_points_density.dat", ["N_cm^-3"]),
        ):
            logical = "density" if product == "density" else "flux"
            coarse_paths = _branch_product_paths(coarse_root / item, logical)
            fine_paths = _branch_product_paths(fine_root / item, logical)
            coarse = _index_unique(
                _read_all_point_rows(coarse_paths), "TimeUTC", coarse_paths[0])
            fine = _index_unique(
                _read_all_point_rows(fine_paths), "TimeUTC", fine_paths[0])
            common = sorted(set(coarse) & set(fine))
            if set(coarse) != coarse_epochs or set(fine) != fine_epochs:
                raise ValueError(f"{item}/{product}: refinement products omit prepared epochs")
            if len(common) < GATES["min_matched_epochs"]:
                raise ValueError(f"{item}/{product}: only {len(common)} common refinement epochs")
            for utc in common:
                for column in columns:
                    a, b = _number(coarse[utc][column], column), _number(fine[utc][column], column)
                    rows.append({"spacecraft": item, "product": product, "utc": utc,
                                 "column": column, "coarse": a, "fine": b,
                                 "fractional_change": fractional_change(a, b)})
    maximum = max(row["fractional_change"] for row in rows)
    result = {
        "schema": "amps-sa01-convergence-v1", "created_utc": utc_now(),
        "analysis_source_sha256": {name: sha256_file(Path(__file__).with_name(name))
                                   for name in ("compare_SA01.py", "sa01_common.py")},
        "coarse_root": str(coarse_root.resolve()), "fine_root": str(fine_root.resolve()),
        "common_population_rows": len(rows), "maximum_fractional_change": maximum,
        "evidence": evidence,
        "threshold": GATES["convergence_fractional_change_max"],
        "status": "PASS" if maximum <= GATES["convergence_fractional_change_max"] else "FAIL",
        "absolute_near_zero_criterion": 1e-30,
    }
    write_csv(destination.with_suffix(".csv"), rows)
    atomic_json(destination.with_suffix(".json"), result)
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    source = result.add_mutually_exclusive_group(required=True)
    source.add_argument("--data-root", type=Path)
    source.add_argument("--data-archive", type=Path,
                        help="Alias for an extracted package directory; not a tar archive.")
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--spacecraft", nargs="+", default=list(SPACECRAFT))
    result.add_argument("--convergence-with", type=Path,
                        help="Fine-resolution output root; compare on the common epoch population.")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    spacecraft = validate_spacecraft(args.spacecraft)
    root = args.output_dir.expanduser().resolve()
    if args.convergence_with:
        result = compare_convergence(root, args.convergence_with.expanduser().resolve(), spacecraft,
                                     root / "convergence_summary")
    else:
        result = compare_all(args.data_root or args.data_archive, root, spacecraft)
    print(f"[{result['status']}] {root}")
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, KeyError, ValueError) as error:
        print(f"ERROR: {error}")
        raise SystemExit(2)

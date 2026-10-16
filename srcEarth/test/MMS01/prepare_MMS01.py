#!/usr/bin/env python3
"""Authenticate and extract the frozen MMS1/HPCA calculation population."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Optional

import numpy as np

from mms01_common import (
    ANODE_CORRECTION, ANODES, DELTA_ALPHA_RAD, DELTA_BETA_RAD, ENERGIES,
    CDFLIB_VERSION, EVENT_START, EVENT_STOP, EXPECTED_DISTRIBUTIONS, EXPECTED_SWEEPS, ION_FILE,
    MEC_FILE, PROFILES, REQUIRED_QUALITY, REQUIRED_SCIENCE_MODE,
    SWEEPS_PER_DISTRIBUTION, atomic_json, ensure_cdflib, resolve_event_root, sha256_file,
    utc_now, verify_sources, write_csv,
)


def _cdf_modules():
    """Import the exact reader selected by the shared dependency contract."""
    cdflib = ensure_cdflib()
    return cdflib.CDF, cdflib.cdfepoch


def _utc_strings(values: np.ndarray) -> list[str]:
    """Format UTC timestamps without host-timezone conversion or midpoint truncation."""
    return [np.datetime_as_string(value, unit="us") + "Z" for value in values]


def _angular_weights(elevation_deg: np.ndarray) -> np.ndarray:
    """Evaluate CMAD equation 5-15's calibrated per-look solid-angle factors."""
    if elevation_deg.shape != (ANODES,) or not np.all(np.isfinite(elevation_deg)):
        raise ValueError("HPCA elevation-centroid array is invalid")
    correction = np.asarray(ANODE_CORRECTION, dtype=float)
    weights = (2.0 * DELTA_ALPHA_RAD * math.sin(DELTA_BETA_RAD / 2.0)
               * np.sin(np.deg2rad(elevation_deg)) / correction)
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("CMAD angular response produced an invalid weight")
    return weights


def _interpolate_positions(mec_path: Path, centers: np.ndarray) -> np.ndarray:
    """Linearly interpolate authoritative 30 s MEC GSE positions at acquisition midpoints.

    MEC position is a smooth Cartesian ephemeris, so component-wise linear
    interpolation is the documented local matching rule.  Extrapolation is rejected.
    Attitude is deliberately not read because this test integrates scalar intensity
    over the calibrated local sky and makes no directional frame comparison.
    """
    CDF, cdfepoch = _cdf_modules()
    mec = CDF(str(mec_path))
    if mec.varattsget("mms1_mec_r_gse").get("UNITS") != "km":
        raise ValueError("MEC GSE position unit is not km")
    times = np.asarray(cdfepoch.to_datetime(mec.varget("Epoch"))).astype("datetime64[ns]")
    positions = np.asarray(mec.varget("mms1_mec_r_gse"), dtype=float)
    targets = centers.astype("datetime64[ns]")
    if targets.min() < times.min() or targets.max() > times.max():
        raise ValueError("MEC ephemeris does not bracket every comparison midpoint")
    x = times.astype(np.int64).astype(float)
    tx = targets.astype(np.int64).astype(float)
    result = np.column_stack([np.interp(tx, x, positions[:, axis]) for axis in range(3)])
    if not np.all(np.isfinite(result)):
        raise ValueError("non-finite interpolated MMS1 position")
    return result


def prepare(data_root: Path, output_root: Path, profile_name: str) -> dict[str, Any]:
    """Create immutable native input and independent Level-2 reference tables."""
    profile_name = profile_name.upper()
    if profile_name not in PROFILES:
        raise ValueError(f"profile must be one of {sorted(PROFILES)}")
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"{output_root} is nonempty; prior evidence is immutable")
    output_root.mkdir(parents=True, exist_ok=True)
    input_dir = output_root / "input"
    reference_dir = output_root / "reference"
    input_dir.mkdir()
    reference_dir.mkdir()

    event_root = resolve_event_root(data_root)
    sources = verify_sources(event_root)
    CDF, cdfepoch = _cdf_modules()
    ion = CDF(str(event_root / ION_FILE))
    epochs_all = np.asarray(cdfepoch.to_datetime(ion.varget("Epoch"))).astype("datetime64[ms]")
    selected = ((epochs_all >= np.datetime64(EVENT_START)) &
                (epochs_all < np.datetime64(EVENT_STOP)))
    indices = np.flatnonzero(selected)
    if indices.size != EXPECTED_SWEEPS or not np.all(np.diff(indices) == 1):
        raise ValueError("frozen HPCA interval is not 624 contiguous sweep records")
    first, last = int(indices[0]), int(indices[-1])
    epochs = epochs_all[first:last + 1]

    def records(name: str, dtype: Optional[Any] = None) -> np.ndarray:
        value = np.asarray(ion.varget(name, startrec=first, endrec=last))
        return value.astype(dtype) if dtype is not None else value

    quality = records("mms1_hpca_hplus_data_quality")
    mode = records("mms1_hpca_science_mode")
    if not np.all(quality == REQUIRED_QUALITY) or not np.all(mode == REQUIRED_SCIENCE_MODE):
        raise ValueError("HPCA quality/science-mode gate failed")
    decimation = {
        "azimuth": records("mms1_hpca_azimuth_decimation_factor"),
        "polar": records("mms1_hpca_polar_decimation_factor"),
        "energy": records("mms1_hpca_energy_decimation_factor"),
    }
    for name, expected in (("azimuth", 2), ("polar", 2), ("energy", 4)):
        if not np.all(decimation[name] == expected):
            raise ValueError(f"HPCA {name} decimation differs from frozen value {expected}")
    if not np.all(records("mms1_hpca_sweep_table_number") == 3):
        raise ValueError("HPCA sweep-table number differs from frozen value 3")

    start_azimuth = records("mms1_hpca_start_azimuth").astype(int)
    for block in range(EXPECTED_DISTRIBUTIONS):
        values = start_azimuth[block * 16:(block + 1) * 16]
        if sorted(values.tolist()) != list(range(16)):
            raise ValueError(f"distribution {block} does not contain every start azimuth")
    cadence_ms = np.diff(epochs).astype("timedelta64[ms]").astype(int)
    if not np.all((cadence_ms == 625) | (cadence_ms == 624)):
        raise ValueError("HPCA sweep cadence is not the expected 625 ms within UTC rounding")

    energy = np.asarray(ion.varget("mms1_hpca_ion_energy"), dtype=float)
    if energy.shape != (ENERGIES,) or not np.all(np.isfinite(energy)) or not np.all(np.diff(energy) > 0):
        raise ValueError("HPCA energy-centroid support is invalid")
    flux = records("mms1_hpca_hplus_flux", float)
    psd = records("mms1_hpca_hplus_phase_space_density", float)
    expected_shape = (EXPECTED_SWEEPS, ANODES, ENERGIES)
    if flux.shape != expected_shape or psd.shape != expected_shape:
        raise ValueError("HPCA directional product shape changed")
    if (not np.all(np.isfinite(flux)) or not np.all(np.isfinite(psd)) or
            np.any(flux < 0.0) or np.any(psd < 0.0)):
        raise ValueError("HPCA directional products contain invalid values")
    if not np.array_equal(flux == 0.0, psd == 0.0):
        raise ValueError("HPCA flux/PSD zero masks differ")
    if ion.varattsget("mms1_hpca_hplus_flux").get("UNITS") != "1/(cm^2 s sr eV)":
        raise ValueError("unexpected HPCA H+ flux unit")
    if ion.varattsget("mms1_hpca_hplus_phase_space_density").get("UNITS") != "cm^-6 s^3":
        raise ValueError("unexpected HPCA H+ PSD unit")
    if ion.varattsget("mms1_hpca_ion_energy").get("UNITS") != "eV/q":
        raise ValueError("unexpected HPCA energy unit")

    elevations = np.asarray(ion.varget("mms1_hpca_centroid_elevation_angle"), dtype=float)
    anode_weights = _angular_weights(elevations)
    chosen = PROFILES[profile_name]
    # The CDF declares a 9.999 s acquisition support.  Work in microseconds so the
    # midpoint retains the final 0.5 ms instead of being truncated to an integer ms.
    centers = np.asarray([
        epochs[block * 16].astype("datetime64[us]") + np.timedelta64(4_999_500, "us")
        for block in chosen
    ])
    positions = _interpolate_positions(event_root / MEC_FILE, centers)
    starts = epochs[np.asarray(chosen) * 16]
    start_strings = _utc_strings(starts)
    center_strings = _utc_strings(centers)

    native_path = input_dir / "hpca_psd_native_input.dat"
    directional_path = reference_dir / "directional_reference.csv"
    spectrum_path = reference_dir / "spectrum_reference.csv"
    position_path = reference_dir / "comparison_positions.csv"
    directional_rows: list[dict[str, Any]] = []
    spectrum_rows: list[dict[str, Any]] = []
    with native_path.open("w", encoding="ascii") as native:
        native.write("# distribution sweep start_azimuth anode energy_index energy_eV psd_cm^-6_s^3 weight_sr\n")
        for output_index, source_block in enumerate(chosen):
            block_start = source_block * SWEEPS_PER_DISTRIBUTION
            for local_sweep in range(SWEEPS_PER_DISTRIBUTION):
                sweep = block_start + local_sweep
                for anode in range(ANODES):
                    weight = float(anode_weights[anode])
                    for energy_index, energy_eV in enumerate(energy):
                        native.write(
                            f"{output_index} {local_sweep} {start_azimuth[sweep]} {anode} "
                            f"{energy_index} {energy_eV:.17g} {psd[sweep, anode, energy_index]:.17g} "
                            f"{weight:.17g}\n"
                        )
                        directional_rows.append({
                            "distribution_index": output_index,
                            "source_distribution_index": source_block,
                            "sweep_index": local_sweep,
                            "start_azimuth": int(start_azimuth[sweep]),
                            "anode": anode,
                            "energy_index": energy_index,
                            "energy_eV": f"{energy_eV:.17g}",
                            "reference_differential_intensity":
                                f"{flux[sweep, anode, energy_index]:.17g}",
                        })
            for energy_index, energy_eV in enumerate(energy):
                terms = []
                for local_sweep in range(SWEEPS_PER_DISTRIBUTION):
                    sweep = block_start + local_sweep
                    terms.extend(float(flux[sweep, anode, energy_index]) *
                                 float(anode_weights[anode]) for anode in range(ANODES))
                spectrum_rows.append({
                    "distribution_index": output_index,
                    "source_distribution_index": source_block,
                    "utc_start": start_strings[output_index],
                    "utc_center": center_strings[output_index],
                    "energy_index": energy_index,
                    "energy_eV": f"{energy_eV:.17g}",
                    "reference_angular_integral": f"{math.fsum(terms):.17g}",
                })

    write_csv(directional_path, directional_rows[0].keys(), directional_rows)
    write_csv(spectrum_path, spectrum_rows[0].keys(), spectrum_rows)
    position_rows = [{
        "distribution_index": index,
        "source_distribution_index": source_block,
        "utc_center": center_strings[index],
        "x_gse_km": f"{positions[index, 0]:.17g}",
        "y_gse_km": f"{positions[index, 1]:.17g}",
        "z_gse_km": f"{positions[index, 2]:.17g}",
    } for index, source_block in enumerate(chosen)]
    write_csv(position_path, position_rows[0].keys(), position_rows)

    artifacts = (native_path, directional_path, spectrum_path, position_path)
    manifest = {
        "schema": "amps-mms01-prepared-v1",
        "created_utc": utc_now(),
        "case_id": "MMS01",
        "profile": profile_name,
        "scientific_scope": [
            "local HPCA H+ phase-space-density to differential-intensity conversion",
            "calibrated angular integration at the 63 archived energy centroids",
        ],
        "explicitly_excluded": [
            "trajectory or boundary-to-MMS transport validation",
            "number density and any energy-integrated flux",
            "absolute instrument calibration independence",
        ],
        "source_files": sources,
        "cdf_reader": {"name": "cdflib", "version": CDFLIB_VERSION},
        "source_event_sweeps": EXPECTED_SWEEPS,
        "source_event_distributions": EXPECTED_DISTRIBUTIONS,
        "selected_source_distribution_indices": list(chosen),
        "selected_distributions": len(chosen),
        "directional_rows": len(directional_rows),
        "spectrum_rows": len(spectrum_rows),
        "energy_centroids": ENERGIES,
        "angular_weight_sum_sr": float(anode_weights.sum() * 16.0),
        "artifacts": {path.name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
                      for path in artifacts},
    }
    atomic_json(output_root / "prepared_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--profile", choices=sorted(PROFILES), required=True)
    args = parser.parse_args()
    manifest = prepare(args.data_root, args.output_root.resolve(), args.profile)
    print(f"Prepared MMS01 {manifest['profile']}: {manifest['selected_distributions']} distributions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

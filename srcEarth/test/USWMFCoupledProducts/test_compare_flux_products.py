#!/usr/bin/env python3
"""Manufactured references and fail-closed cases for the Step-11 set comparator."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from compare_flux_products import ComparisonFailure, compare_product_set


SNAPSHOT = "field-v1-reference"
EPOCH = "2024-05-10T12:00:00.000000000Z"
SUFFIX = ".swmf_t0000003600.000000000s_sidfield-v1-reference"


def require_failure(function, label: str) -> None:
    try:
        function()
    except ComparisonFailure:
        return
    raise AssertionError(f"comparison unexpectedly passed: {label}")


def header(variables: list[str], zone: str = 'ZONE T="reference" I=1 F=POINT') -> str:
    aux = {
        "PHASE_1_INTERPRETATION": "INSTANTANEOUS_QUASI_STATIC",
        "STEP6_CHARACTERISTIC_MAPPING": "STATIC_MAGNETIC",
        "SNAPSHOT_ID": SNAPSHOT,
        "SNAPSHOT_EPOCH_UTC": EPOCH,
        "SNAPSHOT_MESH_REVISION": "swmf-mesh-v1-reference",
        "SNAPSHOT_CONTENT_FINGERPRINT": "swmf-state-v1-reference",
        "OUTER_BOUNDARY_POLICY": "BOX",
        "BOUNDARY_SPECTRUM_EVALUATION_EPOCH_UTC": EPOCH,
        "ACTIVE_SPECTRUM_TABLE_EPOCH_UTC": EPOCH,
        "PRODUCT_CONTROL_FINGERPRINT": "products-v1-reference",
        "BOUNDARY_SPECTRUM_FINGERPRINT": "spectrum-v1-reference",
        "CHANNEL_SCHEMA_FINGERPRINT": "channel-v1-reference",
        "DETECTOR_RESPONSE_FINGERPRINT": "response-v1-reference",
        "OBSERVATION_STATE_FINGERPRINT": "observation-v1-reference",
        "SPECTRUM_ENERGY_BASIS": "PER_PARTICLE",
        "SPECTRUM_MASS_NUMBER": "1",
        "SPECTRUM_INTENSITY_UNIT": "m^-2 s^-1 sr^-1 MeV^-1",
        "SPECTRUM_RELATIVE_UNCERTAINTY": "0",
        "SPECTRUM_TEMPORAL_STATUS": "EXACT",
        "SPECTRUM_TEMPORAL_GAP": "0",
        "SPECTRUM_TEMPORAL_FRACTION": "0",
        "PLANAR_FLUX_CONVENTION": "ISOTROPIC_EQUIVALENT_PI_J",
        "MAXIMUM_UNRESOLVED_FRACTION": "0.01",
        "RESPONSE_WEIGHTED_UNRESOLVED_UPPER_BOUND": "0.01",
        "UNRESOLVED_TOLERANCE": "0.01",
    }
    lines = ['TITLE="manufactured Step-11 reference"']
    lines.extend(f'AUXDATA {key}="{value}"' for key, value in aux.items())
    lines.append("VARIABLES=" + " ".join(f'\"{value}\"' for value in variables))
    lines.append(zone)
    return "\n".join(lines) + "\n"


def write_product_set(live: Path, replay: Path) -> tuple[Path, dict]:
    products = [
        (
            "DENSITY",
            "mode3d_points_density",
            ["X_km", "N_m^-3", "N_lower_m^-3", "N_upper_m^-3"],
            "7000 1.5 1.4 1.6\n",
        ),
        (
            "SPECTRUM",
            "mode3d_points_spectrum",
            [
                "E_MeV", "T", "T_lower", "T_upper", "unresolved_fraction",
                "J_boundary_perMeV", "J_local_perMeV", "J_local_lower_perMeV",
                "J_local_upper_perMeV", "J_omni_perMeV",
                "J_omni_lower_perMeV", "J_omni_upper_perMeV",
                "J_planar_perMeV", "J_planar_lower_perMeV",
                "J_planar_upper_perMeV",
            ],
            "10 0.5 0.49 0.51 0.01 4 2 1.96 2.04 25.13 24.63 25.64 6.28 6.16 6.41\n"
            "20 0.75 0.74 0.76 0.01 4 3 2.96 3.04 37.70 37.20 38.20 9.42 9.30 9.55\n",
        ),
        (
            "FLUX",
            "mode3d_points_flux",
            [
                "X_km", "F_tot_m2s1", "F_tot_lower_m2s1", "F_tot_upper_m2s1",
                "F_P15_25_m2s1", "F_P15_25_lower_m2s1",
                "F_P15_25_upper_m2s1", "F_planar_m2s1",
                "F_planar_lower_m2s1", "F_planar_upper_m2s1", "R_TOPHAT_s1",
                "R_TOPHAT_lower_s1", "R_TOPHAT_upper_s1",
            ],
            "7000 12 11.8 12.2 8 7.8 8.2 3 2.95 3.05 0.2 0.19 0.21\n",
        ),
        (
            "TERMINATION",
            "mode3d_termination_summary",
            [
                "location_index", "E_MeV", "N_sampled", "N_retried",
                "N_resolved", "N_allowed", "T", "T_lower", "T_upper",
                "unresolved_fraction",
            ],
            "0 10 100 1 99 80 0.8080808 0.8 0.81 0.01\n",
        ),
    ]
    artifacts = []
    for role, stem, variables, rows in products:
        live_path = live / f"{stem}{SUFFIX}.dat"
        replay_path = replay / f"{stem}.dat"
        zone = f'ZONE T="{role.lower()}" I={len(rows.splitlines())} F=POINT'
        text = header(variables, zone) + rows
        live_path.write_text(text, encoding="utf-8")
        replay_path.write_text(text, encoding="utf-8")
        artifacts.append({"role": role, "path": live_path.name})

    manifest = {
        "schema": "sep-in-geospace/swmf-coupled-products/v1",
        "RESULT": "PASS",
        "phase_1_interpretation": "INSTANTANEOUS_QUASI_STATIC",
        "characteristic_mapping": "STATIC_MAGNETIC",
        "snapshot_id": SNAPSHOT,
        "content_fingerprint": "swmf-state-v1-reference",
        "mesh_revision": "swmf-mesh-v1-reference",
        "epoch_utc": EPOCH,
        "boundary_spectrum_evaluation_epoch_utc": EPOCH,
        "active_spectrum_table_epoch_utc": EPOCH,
        "simulation_time_s": 3600.0,
        "output_mode": "POINTS",
        "outer_boundary_policy": "BOX",
        "output_suffix": SUFFIX,
        "product_control_fingerprint": "products-v1-reference",
        "boundary_spectrum_fingerprint": "spectrum-v1-reference",
        "channel_schema_fingerprint": "channel-v1-reference",
        "detector_response_fingerprint": "response-v1-reference",
        "observation_state_fingerprint": "observation-v1-reference",
        "spectrum_energy_basis": "PER_PARTICLE",
        "spectrum_mass_number": 1.0,
        "spectrum_intensity_unit": "m^-2 s^-1 sr^-1 MeV^-1",
        "spectrum_relative_uncertainty": 0.0,
        "spectrum_temporal_status": "EXACT",
        "spectrum_temporal_gap": False,
        "spectrum_temporal_fraction": 0.0,
        "location_count": 1,
        "energy_count": 2,
        "direction_count": 48,
        "termination": {
            "sampled": 100,
            "retried": 1,
            "resolved": 99,
            "allowed": 80,
            "maximum_unresolved_fraction": 0.01,
            "response_weighted_unresolved_upper_bound": 0.01,
            "unresolved_tolerance": 0.01,
            "counts": [80, 19, 1],
        },
        "channels": ["P15_25"],
        "detector_responses": ["TOPHAT"],
        "artifacts": artifacts,
        "message": "manufactured reference",
    }
    manifest_path = live / f"swmf_flux_spectrum_manifest{SUFFIX}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest_path, manifest


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="step11_compare_") as directory:
        root = Path(directory)
        live, replay = root / "live", root / "replay"
        live.mkdir()
        replay.mkdir()
        manifest_path, manifest = write_product_set(live, replay)

        result = compare_product_set(manifest_path, replay)
        assert result["RESULT"] == "PASS"
        assert result["artifacts_compared"] == 4
        assert result["maximum_absolute_difference"] == 0.0
        print("PASS S11-U05 exact complete-set live/replay reference")

        flux = replay / "mode3d_points_flux.dat"
        original_flux = flux.read_text(encoding="utf-8")
        flux.write_text(original_flux.replace("7000 12 11.8", "7000 12.1 11.8"),
                        encoding="utf-8")
        require_failure(lambda: compare_product_set(manifest_path, replay),
                        "numeric change under exact default")
        # Tolerance is honored only when explicitly supplied; this call demonstrates
        # the mechanism without changing the default or any registered campaign gate.
        compare_product_set(manifest_path, replay, rtol=0.01)
        flux.write_text(original_flux, encoding="utf-8")
        print("PASS S11-U06 exact default and explicit-only tolerance")

        spectrum = replay / "mode3d_points_spectrum.dat"
        original_spectrum = spectrum.read_text(encoding="utf-8")
        spectrum.write_text(
            original_spectrum.replace(
                'AUXDATA BOUNDARY_SPECTRUM_EVALUATION_EPOCH_UTC="' + EPOCH + '"\n',
                "",
            ),
            encoding="utf-8",
        )
        require_failure(lambda: compare_product_set(manifest_path, replay),
                        "missing boundary-spectrum epoch")
        spectrum.write_text(original_spectrum, encoding="utf-8")

        flux.write_text(
            original_flux.replace(' "R_TOPHAT_s1"', "").replace(
                "3.05 0.2 0.19 0.21", "3.05 0.19 0.21"
            ),
            encoding="utf-8",
        )
        require_failure(lambda: compare_product_set(manifest_path, replay),
                        "missing configured detector output")
        flux.write_text(original_flux, encoding="utf-8")
        print("PASS S11-U07 missing time/schema/provenance fails closed")

        damaged = dict(manifest)
        damaged["termination"] = dict(manifest["termination"])
        damaged["termination"]["maximum_unresolved_fraction"] = 0.010001
        manifest_path.write_text(json.dumps(damaged), encoding="utf-8")
        require_failure(lambda: compare_product_set(manifest_path, replay),
                        "relaxed unresolved acceptance")
        print("PASS S11-U08 unresolved release gate retained")

    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

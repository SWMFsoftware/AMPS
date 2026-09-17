#!/usr/bin/env python3
"""Independent analytical background and characteristic-speed oracle."""

from __future__ import annotations

import math

AU_M = 149_597_870_700.0
SOLAR_RADIUS_M = 695_700_000.0
PROTON_MASS_KG = 1.672_621_925_95e-27
BOLTZMANN_J_K = 1.380_649e-23
MU0_H_M = 1.256_637_061_27e-6
OMEGA_RAD_S = 2.86533e-6
GAMMA_AD = 5.0 / 3.0


def leblanc_density_cm3(radius_au: float, n1au_cm3: float = 5.0) -> float:
    """Evaluate the published three-term profile normalized exactly at 1 AU."""

    radius_rs = radius_au * AU_M / SOLAR_RADIUS_M
    one_au_rs = AU_M / SOLAR_RADIUS_M
    raw = 3.3e5 / radius_rs**2 + 4.1e6 / radius_rs**4 + 8.0e7 / radius_rs**6
    raw_1au = 3.3e5 / one_au_rs**2 + 4.1e6 / one_au_rs**4 + 8.0e7 / one_au_rs**6
    return n1au_cm3 * raw / raw_1au


def background_state(radius_au: float, latitude_deg: float) -> dict[str, float]:
    """Return the fixed default SWCME background from independent equations."""

    speed_km_s, temperature_k, b1au_nt = 400.0, 1.2e5, 5.0
    density_cm3 = leblanc_density_cm3(radius_au)
    pressure_pa = density_cm3 * 1.0e6 * BOLTZMANN_J_K * temperature_k
    sound_km_s = math.sqrt(GAMMA_AD * BOLTZMANN_J_K * temperature_k / PROTON_MASS_KG) / 1000.0
    # The 1-D API interprets B1AU on its configured fixed-latitude ray, so the
    # same local sin(colatitude) enters both 1-AU normalization and evaluation.
    reference_pitch = (OMEGA_RAD_S * AU_M * math.cos(math.radians(latitude_deg)) /
                       (speed_km_s * 1000.0))
    br_1au_nt = b1au_nt / math.sqrt(1.0 + reference_pitch**2)
    br_nt = br_1au_nt / radius_au**2
    local_pitch = OMEGA_RAD_S * radius_au * AU_M * math.cos(math.radians(latitude_deg)) / (speed_km_s * 1000.0)
    bmag_nt = abs(br_nt) * math.sqrt(1.0 + local_pitch**2)
    mass_density = density_cm3 * 1.0e6 * PROTON_MASS_KG
    alfven_km_s = bmag_nt * 1.0e-9 / math.sqrt(MU0_H_M * mass_density) / 1000.0
    return {
        "density_cm3": density_cm3,
        "speed_km_s": speed_km_s,
        "pressure_pa": pressure_pa,
        "sound_speed_km_s": sound_km_s,
        "alfven_speed_km_s": alfven_km_s,
        # This scalar is the perpendicular fast-mode upper branch used by the
        # existing shock preparation when only magnitudes are summarized.
        "fast_speed_km_s": math.sqrt(sound_km_s**2 + alfven_km_s**2),
    }


def observed_state(density_cm3: float, temperature_k: float, bmag_nt: float) -> dict[str, float]:
    """Close measured proton-core moments using the same stated proton-only physics."""

    pressure = density_cm3 * 1.0e6 * BOLTZMANN_J_K * temperature_k
    sound = math.sqrt(GAMMA_AD * BOLTZMANN_J_K * temperature_k / PROTON_MASS_KG) / 1000.0
    alfven = bmag_nt * 1.0e-9 / math.sqrt(MU0_H_M * density_cm3 * 1.0e6 * PROTON_MASS_KG) / 1000.0
    return {"pressure_pa": pressure, "sound_speed_km_s": sound,
            "alfven_speed_km_s": alfven, "fast_speed_km_s": math.sqrt(sound**2 + alfven**2)}

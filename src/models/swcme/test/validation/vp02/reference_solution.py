#!/usr/bin/env python3
"""Independent Parker-spiral equations for VP02."""

from __future__ import annotations

import math

AU_M = 149_597_870_700.0
SOLAR_RADIUS_M = 695_700_000.0
# This is the model's documented rotation convention, copied as a physical
# input rather than importing any production header or generated output.
OMEGA_RAD_S = 2.86533e-6
# SWCME's released Parker background winds from the coordinate origin.  The
# independent oracle records that convention explicitly instead of importing
# the production default, so a future default change is visible as a mismatch.
SOURCE_RADIUS_RS = 0.0


def parker_pitch_ratio(radius_au: float, latitude_deg: float, speed_km_s: float) -> float:
    """Return the positive winding ratio ``-Bphi/Br`` for an outward sector."""

    if not (radius_au > 0.0 and speed_km_s > 0.0):
        raise ValueError("radius and speed must be positive")
    winding_radius_m = radius_au * AU_M - SOURCE_RADIUS_RS * SOLAR_RADIUS_M
    sin_colatitude = math.cos(math.radians(latitude_deg))
    return OMEGA_RAD_S * winding_radius_m * sin_colatitude / (speed_km_s * 1000.0)


def parker_angle_deg(radius_au: float, latitude_deg: float, speed_km_s: float) -> float:
    """Return the signed outward-sector Parker angle in degrees."""

    return math.degrees(math.atan(parker_pitch_ratio(radius_au, latitude_deg, speed_km_s)))


def fold_observed_angle_deg(br_nt: float, bt_nt: float) -> float:
    """Remove sector sign while retaining departures from Parker winding."""

    if br_nt == 0.0:
        raise ValueError("radial field must be nonzero")
    polarity = 1.0 if br_nt > 0.0 else -1.0
    return math.degrees(math.atan2(-bt_nt * polarity, abs(br_nt)))

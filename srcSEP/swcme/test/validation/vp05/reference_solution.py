#!/usr/bin/env python3
"""Independent sign-aware drag-based propagation reference for VP05."""

from __future__ import annotations

import math

AU_KM = 149_597_870.7
VSW_KM_S = 400.0
GAMMA_KM_INV = 1.0e-7


def dbm_state(time_s: float, radius0_km: float, speed0_km_s: float) -> tuple[float, float]:
    """Evaluate the exact constant-wind DBM for fast or slow disturbances."""

    delta = speed0_km_s - VSW_KM_S
    if delta == 0.0 or GAMMA_KM_INV == 0.0:
        return radius0_km + speed0_km_s * time_s, speed0_km_s
    x = GAMMA_KM_INV * abs(delta) * time_s
    if 1.0 + x <= 0.0:
        raise ValueError("time lies beyond DBM pole")
    sign = 1.0 if delta > 0.0 else -1.0
    radius = radius0_km + VSW_KM_S * time_s + sign * math.log1p(x) / GAMMA_KM_INV
    speed = VSW_KM_S + delta / (1.0 + x)
    return radius, speed


def arrival(radius0_au: float, speed0_km_s: float, target_au: float) -> tuple[float, float]:
    """Solve the monotone outward DBM trajectory by bracketed bisection."""

    if not (0.0 < radius0_au < target_au and speed0_km_s > 0.0):
        raise ValueError("arrival requires an outward target and positive speed")
    radius0_km, target_km = radius0_au * AU_KM, target_au * AU_KM
    low, high = 0.0, 4.0 * (target_km - radius0_km) / min(speed0_km_s, VSW_KM_S)
    while dbm_state(high, radius0_km, speed0_km_s)[0] < target_km:
        high *= 2.0
    for _ in range(120):
        middle = 0.5 * (low + high)
        if dbm_state(middle, radius0_km, speed0_km_s)[0] < target_km:
            low = middle
        else:
            high = middle
    time_s = 0.5 * (low + high)
    return time_s, dbm_state(time_s, radius0_km, speed0_km_s)[1]

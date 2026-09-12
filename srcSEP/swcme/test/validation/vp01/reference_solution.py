#!/usr/bin/env python3
"""Independent published-equation reference for the VP01 Leblanc profile.

This module intentionally does not import or translate any SWCME production
routine.  Repeating the published coefficients here allows VP01 to detect a
mistake in model normalization, units, or radial-power evaluation.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np


# IAU 2012 nominal astronomical unit and IAU 2015 nominal solar radius.  These
# literal values make the reference executable independent of SWCME headers.
AU_M = 149_597_870_700.0
SOLAR_RADIUS_M = 695_700_000.0
PROTON_MASS_KG = 1.672_621_923_69e-27

# Equation (1) of Leblanc, Dulk, and Bougeret (1998), in cm^-3 when R is solar
# radii.  A single scale factor is fitted through n(1 AU), as prescribed in the
# SWCME model description.
LEBLANC_COEFFICIENTS_CM3 = (3.3e5, 4.1e6, 8.0e7)


def unscaled_leblanc_cm3(radius_au: np.ndarray | float) -> np.ndarray:
    """Evaluate the unnormalised published density at heliocentric radius."""

    radius = np.asarray(radius_au, dtype=float)
    radius_rs = radius * AU_M / SOLAR_RADIUS_M
    a2, a4, a6 = LEBLANC_COEFFICIENTS_CM3
    return a2 / radius_rs**2 + a4 / radius_rs**4 + a6 / radius_rs**6


def normalized_leblanc_cm3(
    radius_au: np.ndarray | float, density_at_1au_cm3: float
) -> np.ndarray:
    """Evaluate the profile after exact normalization at one astronomical unit."""

    if not math.isfinite(density_at_1au_cm3) or density_at_1au_cm3 <= 0.0:
        raise ValueError("density_at_1au_cm3 must be finite and positive")
    raw = unscaled_leblanc_cm3(radius_au)
    return density_at_1au_cm3 * raw / float(unscaled_leblanc_cm3(1.0))


def fit_outer_normalization_cm3(
    radii_au: Sequence[float],
    observed_density_cm3: Sequence[float],
    calibration_radius_min_au: float,
) -> float:
    """Robustly fit only the free amplitude using outer-heliosphere bins.

    The median is taken in logarithmic space because solar-wind density is
    positive and approximately lognormally distributed.  Restricting the fit
    to the outer bins leaves the inner radial behavior as held-out validation.
    """

    radii = np.asarray(radii_au, dtype=float)
    density = np.asarray(observed_density_cm3, dtype=float)
    mask = radii >= calibration_radius_min_au
    if mask.sum() < 2:
        raise ValueError("at least two outer radial bins are required")
    unit_shape = normalized_leblanc_cm3(radii[mask], 1.0)
    log_amplitudes = np.log(density[mask] / unit_shape)
    return float(np.exp(np.median(log_amplitudes)))


def log_log_slope(radii_au: Iterable[float], density_cm3: Iterable[float]) -> float:
    """Return the ordinary least-squares power-law exponent of binned medians."""

    radii = np.asarray(tuple(radii_au), dtype=float)
    density = np.asarray(tuple(density_cm3), dtype=float)
    if len(radii) < 2 or np.any(radii <= 0.0) or np.any(density <= 0.0):
        raise ValueError("positive values in at least two bins are required")
    return float(np.polyfit(np.log(radii), np.log(density), 1)[0])

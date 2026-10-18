#!/usr/bin/env python3
"""Normalized real spherical harmonics for localized H2O activity.

The normative definition is ``HARMONIC_BASIS.md`` and its identifier is
``liouville_cops.real_orthonormal_cs_v1``.  This module intentionally has no
dependency on ``fit_cops.py``. Stage A4 imports this evaluator to construct its
cached facet design matrix, while Stage A3 uses the matching C++ evaluator
through an explicit, validated forward-model coefficient file.

The implementation is vectorized over Cartesian points but uses the explicit
associated-Legendre recurrence from the specification.  It does not delegate
normalization or phase choices to SciPy, avoiding a hidden convention mismatch
with the C++ implementation.
"""

import math

import numpy as np


BASIS_ID = "liouville_cops.real_orthonormal_cs_v1"
MIN_SUPPORTED_DEGREE = 0
MAX_SUPPORTED_DEGREE = 8


class BasisError(ValueError):
    """The requested degree range or coordinate array is invalid."""


def _validate_degrees(min_degree, max_degree):
    """Return integer degree bounds after strict (non-boolean) validation."""
    if (isinstance(min_degree, bool) or isinstance(max_degree, bool)
            or not isinstance(min_degree, (int, np.integer))
            or not isinstance(max_degree, (int, np.integer))):
        raise BasisError("degree bounds must be integers")
    min_degree, max_degree = int(min_degree), int(max_degree)
    if not (MIN_SUPPORTED_DEGREE <= min_degree <= max_degree <= MAX_SUPPORTED_DEGREE):
        raise BasisError("degree range must satisfy 0 <= min_degree <= max_degree <= 8")
    return min_degree, max_degree


def coefficient_count(min_degree=5, max_degree=8):
    """Number of `(l,m)` coefficients in the inclusive degree range."""
    min_degree, max_degree = _validate_degrees(min_degree, max_degree)
    return (max_degree + 1) ** 2 - min_degree ** 2


def modes(min_degree=5, max_degree=8):
    """Return the normative increasing-`l`, increasing-`m` coefficient order."""
    min_degree, max_degree = _validate_degrees(min_degree, max_degree)
    return [(degree, order)
            for degree in range(min_degree, max_degree + 1)
            for order in range(-degree, degree + 1)]


def coefficient_index(degree, order, min_degree=5, max_degree=8):
    """Return the zero-based index of `(degree,order)` in ``modes()``."""
    min_degree, max_degree = _validate_degrees(min_degree, max_degree)
    if (isinstance(degree, bool) or isinstance(order, bool)
            or not isinstance(degree, (int, np.integer))
            or not isinstance(order, (int, np.integer))):
        raise BasisError("degree and order must be integers")
    degree, order = int(degree), int(order)
    if not min_degree <= degree <= max_degree or not -degree <= order <= degree:
        raise BasisError("mode (%d,%d) is outside the requested range" %
                         (degree, order))
    return degree * degree + degree + order - min_degree * min_degree


def _normalization(degree, abs_order):
    """Fully normalized real-harmonic factor, including sqrt(2) for |m|>0."""
    factorial_ratio = math.exp(math.lgamma(degree - abs_order + 1.0)
                               - math.lgamma(degree + abs_order + 1.0))
    value = math.sqrt((2.0 * degree + 1.0) / (4.0 * math.pi)
                      * factorial_ratio)
    return value if abs_order == 0 else math.sqrt(2.0) * value


def evaluate_cartesian(points, min_degree=5, max_degree=8):
    """Evaluate the basis at one or more nonzero Cartesian vectors.

    Parameters
    ----------
    points : array-like, shape (3,) or (n,3)
        Mesh/body-frame vectors.  Only direction matters.
    min_degree, max_degree : int
        Inclusive degree range, restricted to 0..8.

    Returns
    -------
    numpy.ndarray
        Shape ``(n,nmodes)`` for a point table, or ``(nmodes,)`` for one input
        vector.  Columns follow :func:`modes`.

    Notes
    -----
    Each associated-Legendre order keeps only two recurrence arrays.  Apart
    from the returned dense table, temporary storage is O(n), rather than one
    array for every `(l,m)` mode.
    """
    min_degree, max_degree = _validate_degrees(min_degree, max_degree)
    xyz = np.asarray(points, dtype=np.float64)
    scalar = xyz.ndim == 1
    if scalar:
        if xyz.shape != (3,):
            raise BasisError("one Cartesian point must have shape (3,)")
        xyz = xyz.reshape(1, 3)
    elif xyz.ndim != 2 or xyz.shape[1] != 3:
        raise BasisError("Cartesian point table must have shape (n,3)")
    if not np.all(np.isfinite(xyz)):
        raise BasisError("Cartesian coordinates must be finite")

    # Nested hypot avoids avoidable overflow/underflow in x*x+y*y+z*z.
    rho = np.hypot(xyz[:, 0], xyz[:, 1])
    radius = np.hypot(rho, xyz[:, 2])
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise BasisError("each Cartesian point must have finite positive radius")
    mu = np.clip(xyz[:, 2] / radius, -1.0, 1.0)
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    phi = np.where(rho == 0.0, 0.0, phi)  # explicit, deterministic pole convention
    root = np.sqrt(np.maximum(0.0, 1.0 - mu * mu))

    result = np.empty((len(xyz), coefficient_count(min_degree, max_degree)),
                      dtype=np.float64)

    def store(degree, abs_order, legendre, sine, cosine):
        """Store both signed real branches for one P_l^|m| array."""
        if degree < min_degree:
            return
        factor = _normalization(degree, abs_order)
        if abs_order == 0:
            result[:, coefficient_index(degree, 0, min_degree, max_degree)] = factor * legendre
        else:
            result[:, coefficient_index(degree, -abs_order, min_degree, max_degree)] = factor * legendre * sine
            result[:, coefficient_index(degree, +abs_order, min_degree, max_degree)] = factor * legendre * cosine

    # P_m^m is advanced diagonally between orders.  For each fixed order, the
    # l recurrence advances with two arrays and writes requested degrees.
    p_mm = np.ones(len(xyz), dtype=np.float64)
    for abs_order in range(0, max_degree + 1):
        if abs_order > 0:
            p_mm = -(2 * abs_order - 1) * root * p_mm
        sine = np.sin(abs_order * phi) if abs_order else None
        cosine = np.cos(abs_order * phi) if abs_order else None
        store(abs_order, abs_order, p_mm, sine, cosine)
        if abs_order == max_degree:
            continue

        p_lm2 = p_mm
        p_lm1 = (2 * abs_order + 1) * mu * p_mm
        store(abs_order + 1, abs_order, p_lm1, sine, cosine)
        for degree in range(abs_order + 2, max_degree + 1):
            p_lm = (((2 * degree - 1) * mu * p_lm1
                     - (degree + abs_order - 1) * p_lm2)
                    / (degree - abs_order))
            store(degree, abs_order, p_lm, sine, cosine)
            p_lm2, p_lm1 = p_lm1, p_lm

    return result[0] if scalar else result


def evaluate_angles(theta, phi, min_degree=5, max_degree=8):
    """Evaluate from colatitude/longitude using the same Cartesian evaluator."""
    theta, phi = np.broadcast_arrays(np.asarray(theta, dtype=np.float64),
                                     np.asarray(phi, dtype=np.float64))
    if (not np.all(np.isfinite(theta)) or not np.all(np.isfinite(phi))
            or np.any(theta < 0.0) or np.any(theta > math.pi)):
        raise BasisError("angles must be finite and colatitude must lie in [0,pi]")
    points = np.column_stack((np.sin(theta).ravel() * np.cos(phi).ravel(),
                              np.sin(theta).ravel() * np.sin(phi).ravel(),
                              np.cos(theta).ravel()))
    values = evaluate_cartesian(points, min_degree, max_degree)
    shape = theta.shape + (values.shape[-1],)
    return values.reshape(shape)

#!/usr/bin/env python3
"""Paired NLGCE equation reference; not a production library.

Run beside the delivered data files. The two NLGCE backends are checked
against their own fixtures; their disagreement is reported, not corrected.
"""

import argparse
import csv
import hashlib
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy.special import gamma


ROOT = Path(__file__).resolve().parent
NU = 5.0 / 6.0
C = float(gamma(NU) / (2.0 * math.sqrt(math.pi) * gamma(NU - 0.5)))


def coefficients(kind):
    """Parse the printed j,k,l rows into explicit d[i,j,k,l] indices."""
    result = {}
    with (ROOT / f"NLGCE_F_2014_{kind}.csv").open(newline="") as handle:
        reader = csv.DictReader(handle)
        expected_header = ["j", "k", "l"] + [f"d_i{i}" for i in range(6)]
        if reader.fieldnames != expected_header:
            raise ValueError("Coefficient header does not match the data contract")
        for row in reader:
            key = tuple(int(row[name]) for name in ("j", "k", "l"))
            if key in result:
                raise ValueError(f"Duplicate coefficient row {key}")
            result[key] = tuple(float(row[f"d_i{i}"]) for i in range(6))
    expected_keys = {(j, k, l) for j in range(4) for k in range(4) for l in range(3)}
    if set(result) != expected_keys:
        raise ValueError("Coefficient tuples are missing or out of range")
    return result


COEFFICIENTS = {kind: coefficients(kind) for kind in ("parallel", "perpendicular")}


def fitted_decimal(kind, r, fs, e2, ratio):
    """Independent 50-digit arithmetic on literal CSV coefficients and input decimals."""
    with localcontext() as context:
        context.prec = 50
        inputs = [Decimal(str(value)).ln() for value in (r, fs, e2, ratio)]
        powers = []
        for value, length in zip(inputs, (6, 4, 4, 3)):
            row = [Decimal(1)]
            for _ in range(1, length):
                row.append(row[-1] * value)
            powers.append(row)
        total = Decimal(0)
        with (ROOT / f"NLGCE_F_2014_{kind}.csv").open(newline="") as handle:
            for row in csv.DictReader(handle):
                j, k, l = (int(row[name]) for name in ("j", "k", "l"))
                for i in range(6):
                    total += Decimal(row[f"d_i{i}"]) * powers[0][i] * powers[1][j] * powers[2][k] * powers[3][l]
        return float(total.exp())


def fitted(kind, r, fs, e2, ratio):
    """Equation (53), direct compensated sum independent of production Horner code."""
    x = np.log([r, fs, e2, ratio])
    value = math.fsum(
        row[i] * x[0] ** i * x[1] ** j * x[2] ** k * x[3] ** l
        for (j, k, l), row in COEFFICIENTS[kind].items()
        for i in range(6)
    )
    return math.exp(value)


def solve(r, fs, e2, ratio, seed_scale=1.0, cutoff=40.0, tol=1e-10):
    """Equations (48),(50),(51), v=ell_s=B0=1; outputs lambda/ell_s.

    Each component uses t=ln(k*ell_a). These finite integration limits
    are quadrature controls, not a physical truncation of the spectra.
    """
    if not (r > 0 and 0 < fs < 1 and e2 > 0 and ratio > 0 and cutoff > 15):
        raise ValueError("Inadmissible mixed-component reference inputs")
    epsilon = math.sqrt(e2)
    xi = r / (C * epsilon)
    ax = 0.5 * math.sqrt(
        fs / ((xi / (1.0 + xi)) / epsilon + epsilon / (2.0 * xi))
    )
    ap2 = 1.0 / (1.0 / (math.sqrt(ratio) * fs) + 4.0 / (3.0 * (1.0 - fs)))
    omega2 = 1.0 / r**2
    edges = [-cutoff, -15.0, -5.0, 0.0, 5.0, 15.0, cutoff]

    def integrate(function):
        return math.fsum(
            quad(function, a, b, epsabs=tol * 1e-3, epsrel=tol, limit=300)[0]
            for a, b in zip(edges[:-1], edges[1:])
        )

    def predictions(loglambda):
        lp, lt = np.exp(loglambda)

        def weight(t):
            return math.exp(t - NU * np.logaddexp(0.0, 2.0 * t))

        def rates(t):
            return (
                1.0 / lp + math.exp(2.0 * t) * lp / 3.0,
                1.0 / lp + math.exp(2.0 * t) * ratio**2 * lt / 3.0,
            )

        js = integrate(lambda t: weight(t) * rates(t)[0] / (omega2 + rates(t)[0] ** 2))
        j2 = integrate(lambda t: weight(t) * rates(t)[1] / (omega2 + rates(t)[1] ** 2))
        jperp = integrate(lambda t: weight(t) / rates(t)[1])
        return [
            1.0 / (4.0 * ax * omega2 * C * e2 * (fs * js + (1.0 - fs) * j2)),
            2.0 * ap2 * C * e2 * (1.0 - fs) * jperp,
        ]

    def residual(loglambda):
        return loglambda - np.log(predictions(loglambda))

    seed = np.log(
        [fitted("parallel", r, fs, e2, ratio), fitted("perpendicular", r, fs, e2, ratio)]
    ) + math.log(seed_scale)
    result = least_squares(residual, seed, xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=300)
    error = float(np.max(np.abs(residual(result.x))))
    if not result.success or error >= 1e-8:
        raise RuntimeError(f"Nonlinear solve failed: {result.message}; residual={error}")
    lp, lt = np.exp(result.x)
    return {
        "nonlinear_parallel": float(lp),
        "nonlinear_perpendicular": float(lt),
        "ax": ax,
        "a_prime_squared": ap2,
        "maximum_log_residual": error,
    }


def verify():
    fixture = json.loads((ROOT / "paired_benchmark_points.json").read_text())
    results = []
    for point in fixture["nlgce"]:
        inputs = [point[k] for k in ("r", "slab_fraction", "epsilon_squared", "slab_to_2d_length_ratio")]
        for kind in ("parallel", "perpendicular"):
            for name, value in (("sum", fitted(kind, *inputs)), ("Decimal", fitted_decimal(kind, *inputs))):
                if not math.isclose(value, point[f"fit_{kind}"], rel_tol=1e-11):
                    raise AssertionError(f"{point['label']} {kind} {name}")
        direct = solve(*inputs)
        refined = solve(*inputs, cutoff=45.0, tol=3e-11)
        other = solve(*inputs, seed_scale=3.0)
        for name in ("nonlinear_parallel", "nonlinear_perpendicular", "ax", "a_prime_squared"):
            if not math.isclose(direct[name], point[name], rel_tol=5e-8):
                raise AssertionError(f"{point['label']} {name}")
        changes = [abs(direct[k]/candidate[k]-1) for candidate in (refined, other)
                   for k in ("nonlinear_parallel", "nonlinear_perpendicular")]
        if max(changes) > 5e-8:
            raise AssertionError(f"{point['label']} refinement/seed agreement")
        results.append({"label": point["label"], "maximum_log_residual": direct["maximum_log_residual"],
                        "maximum_refinement_or_seed_change": max(changes)})
    return {"checked_paired_states": len(results), "checked_polynomial_coefficients": 576,
            "states": results, "physical_validation_claimed": False}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))

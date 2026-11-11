#!/usr/bin/env python3
"""Independent analytic checks for Stage-8 forced reaction segments.

This file intentionally has no AMPS or MPI dependency.  It defines the
mathematics that the mover-integrated deterministic two-branch executor must
match:

* an incoming carrier is partitioned into exp(-tau) survivor weight and
  -expm1(-tau) reacted weight;
* a reaction point is sampled from the physical exponential distribution
  conditioned on a reaction inside the traversed segment; and
* a piecewise-constant spatial hazard is inverted through cumulative optical
  depth, not by interpolating geometric distance uniformly.

Keeping this oracle outside the production C++ prevents a common error where
the validator repeats the same incorrect inverse or branch formula as the
implementation it is meant to test.
"""

from __future__ import annotations

import math


def reaction_probability(optical_depth: float) -> float:
    """Return 1-exp(-tau) accurately for small positive optical depth."""

    if not math.isfinite(optical_depth) or optical_depth < 0.0:
        raise ValueError("optical depth must be finite and nonnegative")
    return -math.expm1(-optical_depth)


def conditional_optical_depth(unit_uniform: float, optical_depth: float) -> float:
    """Invert the reacted branch's conditional optical-depth CDF.

    For 0 <= theta <= tau,

        F(theta | reaction in segment) = (1-exp(-theta))/(1-exp(-tau)).

    The endpoint u=1 is excluded exactly like the AMPS random-number contract.
    A zero-depth segment has no reacted branch and is therefore rejected.
    """

    if not math.isfinite(unit_uniform) or not 0.0 <= unit_uniform < 1.0:
        raise ValueError("uniform variate must be finite in [0,1)")
    probability = reaction_probability(optical_depth)
    if probability <= 0.0:
        raise ValueError("a zero-depth segment has no conditional reaction point")
    result = -math.log1p(-unit_uniform * probability)
    if not 0.0 <= result <= optical_depth:
        raise ArithmeticError("conditional optical depth escaped its segment")
    return result


def constant_hazard_expected_fraction(optical_depth: float) -> float:
    """Mean geometric segment fraction for a constant conditional hazard."""

    if not math.isfinite(optical_depth) or optical_depth <= 0.0:
        raise ValueError("positive finite optical depth is required")
    # E[t/T | t<T] = 1/tau - 1/(exp(tau)-1).  expm1 retains the cancellation
    # limit near tau=0; the series branch avoids subtracting two large terms.
    if optical_depth < 1.0e-4:
        tau = optical_depth
        return 0.5 - tau / 12.0 + tau**3 / 720.0
    return 1.0 / optical_depth - 1.0 / math.expm1(optical_depth)


def piecewise_position(
    unit_uniform: float,
    first_fraction: float,
    first_optical_depth: float,
    second_optical_depth: float,
) -> float:
    """Map conditional cumulative optical depth into a two-zone path."""

    if not 0.0 < first_fraction < 1.0:
        raise ValueError("piecewise split must be strictly inside the segment")
    if first_optical_depth <= 0.0 or second_optical_depth <= 0.0:
        raise ValueError("both piecewise hazards must be positive")
    total = first_optical_depth + second_optical_depth
    theta = conditional_optical_depth(unit_uniform, total)
    if theta <= first_optical_depth:
        return first_fraction * theta / first_optical_depth
    return first_fraction + (1.0 - first_fraction) * (
        theta - first_optical_depth
    ) / second_optical_depth


def piecewise_expected_fraction(
    first_fraction: float,
    first_optical_depth: float,
    second_optical_depth: float,
) -> float:
    """Analytic mean position for the two-zone conditional distribution."""

    if not 0.0 < first_fraction < 1.0:
        raise ValueError("piecewise split must be strictly inside the segment")
    if first_optical_depth <= 0.0 or second_optical_depth <= 0.0:
        raise ValueError("both piecewise hazards must be positive")

    a = first_fraction
    b = 1.0 - a
    rate1 = first_optical_depth / a
    rate2 = second_optical_depth / b
    e1 = math.exp(-first_optical_depth)
    e2 = math.exp(-second_optical_depth)
    probability = reaction_probability(first_optical_depth + second_optical_depth)

    first_moment = 1.0 / rate1 - (a + 1.0 / rate1) * e1
    second_moment = e1 * (
        a * (1.0 - e2) + 1.0 / rate2 - (b + 1.0 / rate2) * e2
    )
    return (first_moment + second_moment) / probability


def close(left: float, right: float, relative: float = 2.0e-12) -> bool:
    return math.isclose(left, right, rel_tol=relative, abs_tol=1.0e-15)


def check_homogeneous() -> None:
    tau = 0.73
    incoming = 3.75
    probability = reaction_probability(tau)
    survivor = incoming * math.exp(-tau)
    reacted = incoming * probability
    if not close(survivor + reacted, incoming):
        raise AssertionError("deterministic branch carriers do not close")

    samples = 100_000
    positions = [
        conditional_optical_depth((index + 0.5) / samples, tau) / tau
        for index in range(samples)
    ]
    mean = math.fsum(positions) / samples
    expected = constant_hazard_expected_fraction(tau)
    if abs(mean - expected) > 2.0e-10:
        raise AssertionError("homogeneous conditional-position mean is wrong")
    if positions[0] < 0.0 or positions[-1] > 1.0:
        raise AssertionError("homogeneous conditional position left the segment")


def check_piecewise() -> None:
    first_fraction = 0.35
    tau1 = 0.20
    tau2 = 1.40
    samples = 100_000
    positions = [
        piecewise_position(
            (index + 0.5) / samples, first_fraction, tau1, tau2
        )
        for index in range(samples)
    ]
    mean = math.fsum(positions) / samples
    expected = piecewise_expected_fraction(first_fraction, tau1, tau2)
    if abs(mean - expected) > 2.0e-10:
        raise AssertionError("piecewise conditional-position mean is wrong")

    first_probability = reaction_probability(tau1) / reaction_probability(tau1 + tau2)
    first_count = sum(position <= first_fraction for position in positions)
    if abs(first_count - samples * first_probability) > 1.0:
        raise AssertionError("piecewise first-zone probability is wrong")
    if positions[0] < 0.0 or positions[-1] > 1.0:
        raise AssertionError("piecewise conditional position left the segment")


def check_multistep_decay() -> None:
    optical_depths = (0.03, 0.20, 0.75, 0.01, 1.30)
    survivor = 1.0
    reacted = 0.0
    for tau in optical_depths:
        incoming = survivor
        step_reacted = incoming * reaction_probability(tau)
        survivor = incoming * math.exp(-tau)
        reacted += step_reacted
        if not close(survivor + reacted, 1.0):
            raise AssertionError("multistep carrier closure failed")
    expected_survivor = math.exp(-math.fsum(optical_depths))
    if not close(survivor, expected_survivor):
        raise AssertionError("multistep survival is not exponential")
    if not close(reacted, 1.0 - expected_survivor):
        raise AssertionError("multistep reacted carrier is wrong")


def check_rejections() -> None:
    invalid_calls = (
        lambda: reaction_probability(-1.0),
        lambda: reaction_probability(float("nan")),
        lambda: conditional_optical_depth(0.0, 0.0),
        lambda: conditional_optical_depth(1.0, 0.5),
        lambda: piecewise_position(0.5, 0.0, 0.2, 0.4),
        lambda: piecewise_position(0.5, 0.5, 0.0, 0.4),
    )
    for call in invalid_calls:
        try:
            call()
        except ValueError:
            continue
        raise AssertionError("invalid forced-segment input was accepted")


def main() -> int:
    check_homogeneous()
    check_piecewise()
    check_multistep_decay()
    check_rejections()
    print(
        "PASS Stage-8 forced-segment analytic reference "
        "homogeneous=1 piecewise=1 multistep=1 errors=1"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Independent monotone cubic Hermite reference for VP04."""

from __future__ import annotations

from bisect import bisect_right


def endpoint_slope(h0: float, h1: float, d0: float, d1: float) -> float:
    slope = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
    if slope * d0 <= 0.0:
        return 0.0
    if d0 * d1 < 0.0 and abs(slope) > abs(3.0 * d0):
        return 3.0 * d0
    return slope


def slopes(times: list[float], radii: list[float]) -> list[float]:
    """Construct Fritsch-Butland weighted harmonic slopes independently."""

    if len(times) != len(radii) or len(times) < 2:
        raise ValueError("matching arrays with at least two knots are required")
    h = [times[i + 1] - times[i] for i in range(len(times) - 1)]
    if any(value <= 0.0 for value in h) or any(radii[i + 1] < radii[i] for i in range(len(radii) - 1)):
        raise ValueError("times must increase and radii must not decrease")
    delta = [(radii[i + 1] - radii[i]) / h[i] for i in range(len(h))]
    if len(times) == 2:
        return [delta[0], delta[0]]
    result = [0.0] * len(times)
    result[0] = endpoint_slope(h[0], h[1], delta[0], delta[1])
    result[-1] = endpoint_slope(h[-1], h[-2], delta[-1], delta[-2])
    for index in range(1, len(times) - 1):
        if delta[index - 1] == 0.0 or delta[index] == 0.0 or delta[index - 1] * delta[index] <= 0.0:
            result[index] = 0.0
        else:
            w1 = 2.0 * h[index] + h[index - 1]
            w2 = h[index] + 2.0 * h[index - 1]
            result[index] = (w1 + w2) / (w1 / delta[index - 1] + w2 / delta[index])
    return result


def evaluate(times: list[float], radii: list[float], query: float) -> tuple[float, float]:
    """Return radius [Rsun] and derivative [Rsun/s] inside the knot interval."""

    if query < times[0] or query > times[-1]:
        raise ValueError("query outside knot interval")
    tangent = slopes(times, radii)
    index = min(len(times) - 2, max(0, bisect_right(times, query) - 1))
    h = times[index + 1] - times[index]
    s = (query - times[index]) / h
    y0, y1, m0, m1 = radii[index], radii[index + 1], tangent[index], tangent[index + 1]
    radius = ((2*s**3 - 3*s**2 + 1)*y0 + (s**3 - 2*s**2 + s)*h*m0 +
              (-2*s**3 + 3*s**2)*y1 + (s**3 - s**2)*h*m1)
    derivative = (((6*s**2 - 6*s)*y0 + (3*s**2 - 4*s + 1)*h*m0 +
                   (-6*s**2 + 6*s)*y1 + (3*s**2 - 2*s)*h*m1) / h)
    return radius, derivative

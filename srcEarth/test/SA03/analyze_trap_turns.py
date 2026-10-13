#!/usr/bin/env python3
"""Analyze opt-in ``[trap-turn]`` Poincare records from an SA03 AMPS log.

This is a forensic tool, not a classifier.  It deliberately reports measured return
distances and never changes a trajectory termination.  Its purpose is to determine
whether a proposed turn-boundary recurrence gate is supported by the unresolved real
orbit before such a gate is implemented in production physics.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


EARTH_RADIUS_M = 6_371_200.0


def vector(value: str) -> tuple[float, float, float]:
    parts = tuple(float(item) for item in value.split(","))
    if len(parts) != 3:
        raise ValueError(f"expected a three-component vector, received {value!r}")
    return parts


def parse_records(path: Path, direction: int, energy_mev: float) -> list[dict]:
    """Return the selected final-attempt records, rejecting malformed evidence."""
    records: list[dict] = []
    for raw in path.read_text(encoding="utf-8", errors="strict").splitlines():
        if not raw.startswith("[trap-turn] "):
            continue
        fields = {}
        for token in raw.split()[1:]:
            if "=" not in token:
                raise ValueError(f"malformed trap-turn token: {token!r}")
            key, value = token.split("=", 1)
            fields[key] = value
        if int(fields["direction"]) != direction:
            continue
        if not math.isclose(float(fields["particle_energy_MeV"]), energy_mev,
                            rel_tol=0.0, abs_tol=1.0e-9):
            continue
        records.append({
            "revolution": int(fields["revolution"]),
            "trace_s": float(fields["trace_s"]),
            "radius_m": float(fields["radius_m"]),
            "latitude": float(fields["latitude"]),
            "pitch_cos2": float(fields["pitch_cos2"]),
            "x": vector(fields["x_m"]),
            "p": vector(fields["p_unit"]),
            "b": vector(fields["b_unit"]),
        })
    if not records:
        raise ValueError("no trap-turn records matched the requested trajectory")
    expected = list(range(1, len(records) + 1))
    actual = [record["revolution"] for record in records]
    if actual != expected:
        raise ValueError(f"non-contiguous revolution series: {actual[:3]}...{actual[-3:]}")
    return records


def angle_degrees(a: tuple[float, float, float],
                  b: tuple[float, float, float]) -> float:
    # Clamp the dot product because printed unit vectors can exceed one by roundoff.
    dot = sum(x * y for x, y in zip(a, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot))))


def metrics(first: dict, second: dict) -> dict:
    radial_tolerance = max(
        EARTH_RADIUS_M,
        0.20 * max(abs(first["radius_m"]), abs(second["radius_m"])),
    )
    return {
        "radial_norm": abs(second["radius_m"] - first["radius_m"]) / radial_tolerance,
        "latitude_norm": abs(second["latitude"] - first["latitude"]) / 0.20,
        "pitch_norm": abs(second["pitch_cos2"] - first["pitch_cos2"]) / 0.25,
        "position_m": math.sqrt(sum((x - y) ** 2 for x, y in zip(first["x"], second["x"]))),
        "momentum_angle_deg": angle_degrees(first["p"], second["p"]),
        "field_angle_deg": angle_degrees(first["b"], second["b"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    # SA03 has no pre-established problematic seed.  Require the caller to identify
    # the exact direction and energy from preserved diagnostic output rather than
    # silently analyzing the unrelated SA02 control trajectory copied by architecture.
    parser.add_argument("--direction", type=int, required=True)
    parser.add_argument("--energy-mev", type=float, required=True)
    parser.add_argument("--max-lag", type=int, default=64)
    args = parser.parse_args()

    records = parse_records(args.log, args.direction, args.energy_mev)
    candidates = []
    print(f"records={len(records)} trace_s={records[-1]['trace_s']:.17g}")
    print("lag pairs scalar_matches consecutive_scalar_matches best_worst_norm "
          "best_pair momentum_angle_deg position_m")
    for lag in range(1, min(args.max_lag, len(records) - 1) + 1):
        comparisons = [metrics(records[i], records[i + lag])
                       for i in range(len(records) - lag)]
        scalar_pass = [max(item["radial_norm"], item["latitude_norm"],
                           item["pitch_norm"]) <= 1.0 for item in comparisons]
        longest = current = 0
        for passed in scalar_pass:
            current = current + 1 if passed else 0
            longest = max(longest, current)
        scored = [(max(item["radial_norm"], item["latitude_norm"], item["pitch_norm"]),
                   index, item) for index, item in enumerate(comparisons)]
        score, index, best = min(scored)
        candidates.append((score, lag, index, best, sum(scalar_pass), longest))
        print(f"{lag} {len(comparisons)} {sum(scalar_pass)} {longest} {score:.9g} "
              f"{index + 1}:{index + lag + 1} {best['momentum_angle_deg']:.9g} "
              f"{best['position_m']:.9g}")

    print("best_returns")
    for score, lag, index, item, matches, longest in sorted(candidates)[:10]:
        print(f"lag={lag} pair={index + 1}:{index + lag + 1} worst_norm={score:.9g} "
              f"scalar_matches={matches} consecutive={longest} "
              f"radial_norm={item['radial_norm']:.9g} "
              f"latitude_norm={item['latitude_norm']:.9g} "
              f"pitch_norm={item['pitch_norm']:.9g} "
              f"momentum_angle_deg={item['momentum_angle_deg']:.9g} "
              f"field_angle_deg={item['field_angle_deg']:.9g} "
              f"position_m={item['position_m']:.9g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

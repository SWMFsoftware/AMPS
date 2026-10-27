#!/usr/bin/env python3
"""Build explicit Stage-3 matrix weights without extrapolation.

The input source is a probability distribution over exact species, energy,
mu=cos(theta), and phi nodes.  For ``exact_nodes`` every node must exist in the
matrix.  The interpolation policy performs multilinear interpolation of the
*response* in log(E), mu, and periodic phi; the resulting corner weights are
written explicitly for review and consumed by shieldSimResponseTool.

This program does not convert intensity to plane-crossing probability.  Source
adapters must perform and document that physical normalization first.  Keeping
that boundary explicit prevents an instrument flux, omnidirectional flux, or
directional intensity from being silently relabeled as a response weight.
"""

from __future__ import annotations

import argparse
import bisect
import json
import math
import pathlib


def bracket(values: list[float], value: float, *, periodic: bool = False) -> list[tuple[float, float]]:
    """Return ``[(coordinate, interpolation_weight), ...]`` for one axis."""
    values = sorted(set(values))
    if not values:
        raise ValueError("matrix interpolation axis is empty")
    tolerance = 1.0e-12
    for candidate in values:
        if math.isclose(value, candidate, rel_tol=tolerance, abs_tol=tolerance):
            return [(candidate, 1.0)]
    if len(values) == 1:
        raise ValueError(f"query coordinate {value} is absent from a one-node axis")
    if periodic:
        period = 2.0 * math.pi
        value %= period
        extended = values + [values[0] + period]
        shifted = value if value >= values[0] else value + period
        index = bisect.bisect_right(extended, shifted) - 1
        lower, upper = extended[index], extended[index + 1]
        fraction = (shifted - lower) / (upper - lower)
        return [(lower % period, 1.0 - fraction), (upper % period, fraction)]
    if value < values[0] or value > values[-1]:
        raise ValueError(
            f"query coordinate {value} is outside [{values[0]}, {values[-1]}]; extrapolation is prohibited"
        )
    index = bisect.bisect_right(values, value) - 1
    lower, upper = values[index], values[index + 1]
    fraction = (value - lower) / (upper - lower)
    return [(lower, 1.0 - fraction), (upper, fraction)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", type=pathlib.Path)
    parser.add_argument("source_nodes", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument(
        "--policy",
        choices=("exact_nodes", "linear_response_in_log_energy_bilinear_mu_phi"),
        default="linear_response_in_log_energy_bilinear_mu_phi",
    )
    args = parser.parse_args()
    matrix = json.loads(args.matrix.read_text(encoding="utf-8"))
    source = json.loads(args.source_nodes.read_text(encoding="utf-8"))
    if matrix.get("schema_version") != "shieldSim.response_matrix.v1":
        raise ValueError("unsupported matrix schema")
    if source.get("schema_version") != "shieldSim.response_source_nodes.v1":
        raise ValueError("unsupported source-node schema")

    cells = matrix["cells"]
    lookup: dict[tuple[str, int, int, float, float, float], str] = {}
    groups: dict[tuple[str, int, int], list[dict]] = {}
    for cell in cells:
        incident = cell["incident"]
        group = (incident["particle"], int(incident["Z"]), int(incident["A"]))
        groups.setdefault(group, []).append(incident)
        key = group + (
            float(incident["total_kinetic_energy_MeV"]),
            float(incident["mu"]),
            float(incident["phi_rad"]) % (2.0 * math.pi),
        )
        if key in lookup:
            raise ValueError("matrix has duplicate species/energy/direction coordinates")
        lookup[key] = incident["incident_id"]

    accumulated: dict[str, float] = {}
    source_total = 0.0
    for node in source["nodes"]:
        probability = float(node["probability"])
        if not math.isfinite(probability) or probability < 0.0:
            raise ValueError("source-node probabilities must be finite and nonnegative")
        source_total += probability
        group = (node["particle"], int(node["Z"]), int(node["A"]))
        if group not in groups:
            raise ValueError(f"source species {group} is absent from the response matrix")
        energy = float(node["total_kinetic_energy_MeV"])
        mu = float(node["mu"])
        phi = float(node["phi_rad"]) % (2.0 * math.pi)
        if energy <= 0.0 or not (0.0 < mu <= 1.0):
            raise ValueError("source energy must be positive and mu must be in (0,1]")
        available = groups[group]
        energies = sorted({float(item["total_kinetic_energy_MeV"]) for item in available})
        mus = sorted({float(item["mu"]) for item in available})
        phis = sorted({float(item["phi_rad"]) % (2.0 * math.pi) for item in available})
        if args.policy == "exact_nodes":
            key = group + (energy, mu, phi)
            if key not in lookup:
                raise ValueError("exact source node is absent; extrapolation/interpolation was not requested")
            accumulated[lookup[key]] = accumulated.get(lookup[key], 0.0) + probability
            continue
        e_bracket = [
            (math.exp(value), weight)
            for value, weight in bracket([math.log(value) for value in energies], math.log(energy))
        ]
        mu_bracket = bracket(mus, mu)
        phi_bracket = bracket(phis, phi, periodic=len(phis) > 1)
        for e_value, e_weight in e_bracket:
            for mu_value, mu_weight in mu_bracket:
                for phi_value, phi_weight in phi_bracket:
                    candidates = [
                        (key, incident_id)
                        for key, incident_id in lookup.items()
                        if key[:3] == group
                        and math.isclose(key[3], e_value, rel_tol=1e-12)
                        and math.isclose(key[4], mu_value, rel_tol=1e-12)
                        and math.isclose(key[5], phi_value, rel_tol=1e-12, abs_tol=1e-12)
                    ]
                    if len(candidates) != 1:
                        raise ValueError("matrix is not a complete tensor grid around a requested node")
                    incident_id = candidates[0][1]
                    weight = probability * e_weight * mu_weight * phi_weight
                    accumulated[incident_id] = accumulated.get(incident_id, 0.0) + weight

    if not math.isclose(source_total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(f"source probabilities sum to {source_total}, not one")
    total = sum(accumulated.values())
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(f"interpolated matrix weights sum to {total}, not one")
    output = {
        "schema_version": "shieldSim.response_weights.v1",
        "matrix_hash": matrix["matrix_hash"],
        "normalization": source["normalization"],
        "quadrature_policy": source["quadrature_policy"],
        "interpolation_policy": args.policy,
        "weights": [
            {"incident_id": incident_id, "weight": weight}
            for incident_id, weight in sorted(accumulated.items())
            if weight > 0.0
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

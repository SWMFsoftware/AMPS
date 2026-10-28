#!/usr/bin/env python3
"""Report geometry invariants for an ASCII NASTRAN triangular surface mesh."""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path
from typing import Dict, Iterable, TextIO, Tuple

Vector = Tuple[float, float, float]
Triangle = Tuple[int, int, int]


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("rt", encoding="ascii")


def fields(line: str, width: int = 8) -> Iterable[str]:
    return (line[index : index + width].strip() for index in range(0, len(line), width))


def read_nastran(path: Path) -> Tuple[Dict[int, Vector], list[Triangle]]:
    vertices: Dict[int, Vector] = {}
    triangles: list[Triangle] = []

    with open_text(path) as source:
        iterator = iter(source)
        for line in iterator:
            card = line[:8].strip().upper()
            if card == "GRID*":
                # Long-field GRID stores ID, CP, X and Y in the first record,
                # followed by Z in the continuation record.
                vertex_id = int(line[8:24])
                x = float(line[40:56])
                y = float(line[56:72])
                continuation = next(iterator)
                z = float(continuation[8:24])
                vertices[vertex_id] = (x, y, z)
            elif card == "GRID":
                entry = list(fields(line))
                vertex_id = int(entry[1])
                vertices[vertex_id] = tuple(float(value) for value in entry[3:6])  # type: ignore[assignment]
            elif card == "CTRIA3":
                entry = list(fields(line))
                triangles.append(tuple(int(value) for value in entry[3:6]))  # type: ignore[arg-type]

    if not vertices:
        raise ValueError(f"{path}: no GRID vertices found")
    if not triangles:
        raise ValueError(f"{path}: no CTRIA3 faces found")

    missing = sorted({vertex for face in triangles for vertex in face if vertex not in vertices})
    if missing:
        raise ValueError(f"{path}: faces reference {len(missing)} missing vertices")
    return vertices, triangles


def subtract(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a: Vector, b: Vector) -> Vector:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a: Vector, b: Vector) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def properties(vertices: Dict[int, Vector], triangles: list[Triangle]) -> dict[str, object]:
    area = 0.0
    signed_volume = 0.0
    edge_counts: Dict[Tuple[int, int], int] = {}

    for i, j, k in triangles:
        a, b, c = vertices[i], vertices[j], vertices[k]
        area_vector = cross(subtract(b, a), subtract(c, a))
        area += 0.5 * math.sqrt(dot(area_vector, area_vector))
        signed_volume += dot(a, cross(b, c)) / 6.0
        for edge in ((i, j), (j, k), (k, i)):
            ordered = tuple(sorted(edge))
            edge_counts[ordered] = edge_counts.get(ordered, 0) + 1

    bounds = tuple(
        (min(point[axis] for point in vertices.values()),
         max(point[axis] for point in vertices.values()))
        for axis in range(3)
    )
    volume = abs(signed_volume)
    return {
        "vertices": len(vertices),
        "triangles": len(triangles),
        "edges": len(edge_counts),
        "nonmanifold_edges": sum(count != 2 for count in edge_counts.values()),
        "bounds": bounds,
        "area": area,
        "volume": volume,
        "area_radius": math.sqrt(area / (4.0 * math.pi)),
        "volume_radius": (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument(
        "--input-unit",
        choices=("m", "km"),
        default="km",
        help="coordinate unit in the mesh (default: km)",
    )
    args = parser.parse_args()

    result = properties(*read_nastran(args.mesh))
    scale = 1000.0 if args.input_unit == "km" else 1.0
    bounds = result["bounds"]

    print(f"vertices={result['vertices']}")
    print(f"triangles={result['triangles']}")
    print(f"edges={result['edges']}")
    print(f"nonmanifold_edges={result['nonmanifold_edges']}")
    for axis, name in enumerate("xyz"):
        lower, upper = bounds[axis]  # type: ignore[index]
        print(f"{name}_bounds_m={lower * scale:.12g},{upper * scale:.12g}")
    print(f"surface_area_m2={float(result['area']) * scale**2:.12f}")
    print(f"enclosed_volume_m3={float(result['volume']) * scale**3:.12f}")
    print(f"area_equivalent_radius_m={float(result['area_radius']) * scale:.12f}")
    print(f"volume_equivalent_radius_m={float(result['volume_radius']) * scale:.12f}")


if __name__ == "__main__":
    main()

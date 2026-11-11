#!/usr/bin/env python3
"""Outer-boundary (edge) correction for spherical AMPS comet brightness.

A spherical AMPS run removes every particle that crosses ``mesh.outer_radius_m``.
Species whose lifetime times speed exceeds that radius (H in particular) are
therefore missing from the far coma, and the column on every line of sight is
truncated at the domain edge. This module estimates, for every impact
parameter b, the factor

    C(b) = N_full(b) / N_trunc(b)

from a collisionless steady-state test-particle reference model that uses the
run's own inputs. ``N_full`` is the column of an unbounded coma.
``N_trunc`` counts only atoms whose entire history (parent decay point,
intermediate decay points, and present position) lies inside the AMPS domain,
which is what AMPS can contain. Multiply the AMPS brightness by ``C(b)`` to
estimate the brightness of a domain without an edge.

The reference model reads everything it needs from the run directory:

* the runtime input (``i3.in`` or another ``#begin application: comet`` file):
  the selected case's heliocentric distance, source rates, 1-au reaction rates,
  nucleus radius, and ``mesh.outer_radius_m``;
* the compile-time input (``i3.conf``): reaction topology and the recoil
  (excess) energy of each channel;
* the AMPS radial profile of every source species at the same output, whose
  radial bulk speed carries the parents to their decay radii.

Products receive the AMPS photochemical recoil law: a Gaussian draw in the
3N-dimensional product velocity space, centre-of-mass motion removed, and
rescaling to the full excess energy. Collisions are ignored. See
``srcHartley2/EDGE_CORRECTION.md`` for the derivation, validity, and examples.

Command-line use prints the factor at every impact parameter of an AMPS
brightness table and can write it to CSV::

    python3 srcHartley2/visualization/comet_edge_correction.py \
        RUN/PT/plots/pic.H.brightness.s=0.out=6.dat --case 4
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

AMU = 1.66053906660e-27
EV = 1.602176634e-19
MASS_AMU = {
    "H": 1.008,
    "H2": 2.016,
    "O": 15.999,
    "OH": 17.007,
    "H2O": 18.015,
    "CO": 28.010,
    "CO2": 44.009,
    "C": 12.011,
    "N2": 28.014,
    "CH4": 16.043,
}
DEFAULT_PARENTS = 3_000_000
DEFAULT_SEED = 20261109
R_BIG_M = 1.0e12
AGE_MIN_S, AGE_MAX_S = 1.0, 1.0e9
INNER_SHELLS, OUTER_SHELLS = 700, 300


# --------------------------------------------------------------------------
# Input parsing
# --------------------------------------------------------------------------
@dataclass
class Channel:
    reaction_id: str
    parent: str
    products: List[str]
    energy_j: float
    rate_s: float = 0.0


@dataclass
class RunInputs:
    case: str
    heliocentric_distance_au: float
    nucleus_radius_m: float
    outer_radius_m: float
    sources: Dict[str, float]
    channels: Dict[str, List[Channel]] = field(default_factory=dict)

    def loss_rate(self, species: str) -> float:
        return sum(c.rate_s for c in self.channels.get(species, []))


def parse_products(text: str) -> List[str]:
    """``'2*H + O'`` -> ``['H', 'H', 'O']``; ``'NONE'`` -> ``[]``."""

    text = text.strip()
    if text.upper() == "NONE":
        return []
    products: List[str] = []
    for term in text.split("+"):
        term = term.strip()
        match = re.fullmatch(r"(?:(\d+)\s*\*\s*)?([A-Za-z][A-Za-z0-9]*)", term)
        if not match:
            raise ValueError(f"cannot parse reaction product term {term!r}")
        products += [match.group(2)] * int(match.group(1) or 1)
    return products


def parse_runtime_input(path: Path, case: str | int) -> Dict[str, str]:
    """Common ``key = value`` settings overlaid with the selected case section."""

    common: Dict[str, str] = {}
    cases: Dict[str, Dict[str, str]] = {}
    current: Dict[str, str] | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("!"):
            continue
        begin = re.fullmatch(r"#begin case:\s*(\S+)", line)
        if begin:
            current = cases.setdefault(begin.group(1), {})
            continue
        if line.startswith("#end case"):
            current = None
            continue
        if line.startswith("#"):
            continue
        key, separator, value = line.partition("=")
        if separator:
            (current if current is not None else common)[key.strip()] = value.strip()
    if cases:
        if str(case) not in cases:
            raise ValueError(f"{path}: no '#begin case: {case}' section")
        common.update(cases[str(case)])
    return common


def parse_reactions(path: Path) -> Dict[str, Tuple[str, List[str], float]]:
    """``ChemistryReaction=ID | parent | products | ... | energy_eV``."""

    reactions: Dict[str, Tuple[str, List[str], float]] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line.startswith("ChemistryReaction="):
            continue
        fields = [f.strip() for f in line.partition("=")[2].split("|")]
        if len(fields) < 4:
            raise ValueError(f"{path}: malformed reaction line {raw!r}")
        reactions[fields[0]] = (fields[1], parse_products(fields[2]), float(fields[-1]))
    if not reactions:
        raise ValueError(f"{path}: no ChemistryReaction lines")
    return reactions


def load_run_inputs(runtime_input: Path, build_input: Path, case: str | int) -> RunInputs:
    values = parse_runtime_input(runtime_input, case)
    rh = float(values["heliocentric_distance_au"])
    channels: Dict[str, List[Channel]] = {}
    for rid, (parent, products, energy_ev) in parse_reactions(build_input).items():
        key = f"chemistry.rate.{rid}_1au_s-1"
        if key not in values:
            raise ValueError(f"{runtime_input}: case {case} has no {key}")
        channel = Channel(rid, parent, products, energy_ev * EV, float(values[key]) / rh**2)
        channels.setdefault(parent, []).append(channel)
    sources = {
        m.group(1): float(v)
        for k, v in values.items()
        if (m := re.fullmatch(r"source\.(\w+)\.rate_s-1", k)) and float(v) > 0.0
    }
    return RunInputs(
        case=str(case),
        heliocentric_distance_au=rh,
        nucleus_radius_m=float(values["nucleus.radius_m"]),
        outer_radius_m=float(values["mesh.outer_radius_m"]),
        sources=sources,
        channels=channels,
    )


def read_radial_profile(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return unique radii, number density, and radial bulk speed of an AMPS
    spherical Tecplot profile (block-boundary duplicates averaged)."""

    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if line.startswith("VARIABLES="))
    names = [n.strip().strip('"').strip() for n in header.partition("=")[2].split(",")]
    rows = [line.split() for line in lines if line.strip() and line.strip()[0] in "+-.0123456789"]
    data = np.array([[float(x) for x in row] for row in rows if len(row) == len(names)])
    radius, inverse = np.unique(data[:, 0], return_inverse=True)
    count = np.bincount(inverse)

    def mean(name: str) -> np.ndarray:
        return np.bincount(inverse, weights=data[:, names.index(name)]) / count

    return radius, mean("Number Density[1/m^3]"), mean("Vx [m/s]")


def find_profile(plots_dir: Path, species: str, output: int) -> Path:
    matches = sorted(plots_dir.glob(f"pic.{species}.s=*.out={output}.dat"))
    if len(matches) != 1:
        raise ValueError(f"{plots_dir}: expected one pic.{species}.s=*.out={output}.dat, found {len(matches)}")
    return matches[0]


# --------------------------------------------------------------------------
# Collisionless reference model
# --------------------------------------------------------------------------
def isotropic(n: int, rng: np.random.Generator) -> np.ndarray:
    mu = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - mu**2)
    return np.column_stack([s * np.cos(phi), s * np.sin(phi), mu])


def recoil(masses_kg: Sequence[float], energy_j: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """AMPS recoil law: Gaussian in 3N, CM motion removed, scaled to energy_j.
    Returns velocities (n, len(masses), 3) in the parent rest frame."""

    m = np.asarray(masses_kg)[None, :, None]
    if energy_j <= 0.0 or len(masses_kg) < 2:
        return np.zeros((n, len(masses_kg), 3))
    u = rng.standard_normal((n, len(masses_kg), 3))
    u -= (m * u).sum(axis=1, keepdims=True) / m.sum()
    e = 0.5 * (m * u**2).sum(axis=(1, 2))
    return u * np.sqrt(energy_j / e)[:, None, None]


class ShellAccumulator:
    """Steady-state number per spherical shell, full and edge-truncated."""

    def __init__(self, r_inner: float, r_edge: float, r_big: float = R_BIG_M):
        inner = np.geomspace(r_inner, r_edge, INNER_SHELLS + 1)
        outer = np.geomspace(r_edge, r_big, OUTER_SHELLS + 1)[1:]
        self.r_edge = r_edge
        self.edges = np.concatenate([inner, outer])
        self.full = np.zeros(len(self.edges) - 1)
        self.trunc = np.zeros(len(self.edges) - 1)

    def add(self, x: np.ndarray, weight: np.ndarray, inside: np.ndarray) -> None:
        r = np.linalg.norm(x, axis=1)
        k = np.searchsorted(self.edges, r) - 1
        ok = (k >= 0) & (k < len(self.full))
        np.add.at(self.full, k[ok], weight[ok])
        t = ok & inside & (r <= self.r_edge)
        np.add.at(self.trunc, k[t], weight[t])

    def densities(self) -> Tuple[np.ndarray, np.ndarray]:
        volume = 4.0 / 3.0 * np.pi * np.diff(self.edges**3)
        return self.full / volume, self.trunc / volume

    def columns(self, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Exact projection of the piecewise-constant shell densities."""

        return tuple(project_shells(self.edges, n, b) for n in self.densities())


def project_shells(edges: np.ndarray, density: np.ndarray, b: Iterable[float]) -> np.ndarray:
    r1, r2 = edges[:-1], edges[1:]
    out = []
    for bb in np.asarray(b, dtype=float):
        m = r2 > bb
        seg = np.sqrt(r2[m] ** 2 - bb**2) - np.sqrt(np.maximum(r1[m], bb) ** 2 - bb**2)
        out.append(2.0 * float(np.sum(density[m] * seg)))
    return np.asarray(out)


def reference_columns(
    inputs: RunInputs,
    parent_profiles: Dict[str, Tuple[np.ndarray, np.ndarray]],
    b: np.ndarray,
    observed: Sequence[str] = ("H", "OH"),
    parents: int = DEFAULT_PARENTS,
    seed: int = DEFAULT_SEED,
    r_big: float = R_BIG_M,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Full and truncated columns of each observed species at impact
    parameters ``b``. ``parent_profiles[S] = (radius, radial_speed)``."""

    rng = np.random.default_rng(seed)
    r_edge = inputs.outer_radius_m
    observed = list(observed)

    # Species that are observed or can decay, directly or through a chain,
    # into an observed species.
    relevant = set(observed)
    changed = True
    while changed:
        changed = False
        for parent, chans in inputs.channels.items():
            if parent not in relevant and any(p in relevant for c in chans for p in c.products):
                relevant.add(parent)
                changed = True
    for species in relevant:
        if species not in MASS_AMU:
            raise ValueError(f"no mass for species {species}; add it to MASS_AMU")

    acc = {s: ShellAccumulator(inputs.nucleus_radius_m, r_edge, r_big) for s in observed}

    def snapshot(pos, vel, beta, weight_scale):
        age = np.exp(rng.uniform(np.log(AGE_MIN_S), np.log(AGE_MAX_S), len(pos)))
        weight = weight_scale * age * np.log(AGE_MAX_S / AGE_MIN_S) * np.exp(-beta * age)
        return pos + vel * age[:, None], weight

    def decay_products(species, pos, vel, inside, weight_scale, at_decay=False):
        chans = inputs.channels.get(species, [])
        beta = inputs.loss_rate(species)
        if not chans or beta <= 0.0 or len(pos) == 0:
            return
        if at_decay:  # source parents arrive at their sampled decay points
            q, inq = pos, inside
        else:
            age = rng.exponential(1.0 / beta, len(pos))
            q = pos + vel * age[:, None]
            inq = inside & (np.linalg.norm(q, axis=1) <= r_edge)
        pick = rng.choice(len(chans), size=len(pos), p=np.array([c.rate_s for c in chans]) / beta)
        for j, chan in enumerate(chans):
            if not any(p in relevant for p in chan.products):
                continue
            s = pick == j
            masses = [MASS_AMU[p] * AMU for p in chan.products]
            u = recoil(masses, chan.energy_j, int(s.sum()), rng)
            for k, product in enumerate(chan.products):
                if product in relevant:
                    follow(product, q[s], vel[s] + u[:, k], inq[s], weight_scale)

    def follow(species, pos, vel, inside, weight_scale):
        if species in acc:
            xs, w = snapshot(pos, vel, inputs.loss_rate(species), weight_scale)
            acc[species].add(xs, w, inside)
        decay_products(species, pos, vel, inside, weight_scale)

    for source, rate in inputs.sources.items():
        if source not in relevant:
            continue
        if source in acc:
            raise ValueError(f"source species {source} is also observed; not supported")
        radius, speed = parent_profiles[source]
        beta = inputs.loss_rate(source)
        good = speed > 0.0
        grid = np.geomspace(inputs.nucleus_radius_m, r_big, 20000)
        vg = np.interp(grid, radius[good], speed[good])
        optical = np.concatenate([[0.0], np.cumsum(np.diff(grid) * beta / (0.5 * (vg[1:] + vg[:-1])))])
        tau = rng.exponential(size=parents)
        tau = tau[tau < optical[-1]]
        r_decay = np.interp(tau, optical, grid)
        direction = isotropic(len(r_decay), rng)
        position = direction * r_decay[:, None]
        velocity = direction * np.interp(r_decay, grid, vg)[:, None]
        decay_products(source, position, velocity, r_decay <= r_edge, rate / parents, at_decay=True)

    return {s: acc[s].columns(b) for s in observed}


def correction_factor(full: np.ndarray, trunc: np.ndarray, b: np.ndarray, r_edge: float) -> np.ndarray:
    """``N_full / N_trunc``; NaN where the truncated column vanishes or b is at
    or beyond the domain edge."""

    factor = np.full(len(b), np.nan)
    ok = (trunc > 0.0) & (np.asarray(b) < r_edge * (1.0 - 1.0e-9))
    factor[ok] = full[ok] / trunc[ok]
    return factor


# --------------------------------------------------------------------------
# Run-directory front end
# --------------------------------------------------------------------------
def output_index(path: Path) -> int:
    match = re.search(r"\.out=(\d+)\.dat$", path.name)
    if not match:
        raise ValueError(f"cannot infer the output index from {path.name}")
    return int(match.group(1))


def run_directory(brightness_file: Path) -> Path:
    """``RUN/PT/plots/pic.*.dat`` -> ``RUN``."""

    return brightness_file.resolve().parents[2]


def default_build_input(run_dir: Path) -> Path:
    for candidate in (run_dir / "i3.conf", Path(__file__).resolve().parents[2] / "input/i3.conf"):
        if candidate.is_file():
            return candidate
    raise ValueError(f"no i3.conf in {run_dir} or input/; supply --build-input")


def correction_for_run(
    run_dir: Path,
    case: str | int,
    output: int,
    b: np.ndarray,
    observed: Sequence[str] = ("H", "OH"),
    parents: int = DEFAULT_PARENTS,
    seed: int = DEFAULT_SEED,
    runtime_input: Path | None = None,
    build_input: Path | None = None,
) -> Tuple[Dict[str, np.ndarray], RunInputs]:
    """Edge-correction factor for each observed species at impact parameters b."""

    runtime_input = runtime_input or run_dir / "i3.in"
    build_input = build_input or default_build_input(run_dir)
    inputs = load_run_inputs(runtime_input, build_input, case)
    profiles = {}
    for source in inputs.sources:
        path = find_profile(run_dir / "PT/plots", source, output)
        radius, _, speed = read_radial_profile(path)
        profiles[source] = (radius, speed)
    columns = reference_columns(inputs, profiles, b, observed, parents, seed)
    factors = {s: correction_factor(full, trunc, b, inputs.outer_radius_m) for s, (full, trunc) in columns.items()}
    return factors, inputs


def infer_case(path: Path) -> str:
    match = re.search(r"(?:^|[^a-z0-9])(?:i3[_-]?)?case[_-]?(\d+)(?:[^0-9]|$)", str(path.resolve()), re.I)
    if not match:
        raise ValueError("case could not be inferred from the path; supply --case")
    return match.group(1)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("brightness_file", type=Path, help="AMPS pic.<species>.brightness.s=<i>.out=<K>.dat")
    parser.add_argument("--case", help="case section name (otherwise inferred from the path)")
    parser.add_argument("--run-dir", type=Path, help="run directory (default: three levels above the file)")
    parser.add_argument("--runtime-input", type=Path, help="runtime input (default: RUN/i3.in)")
    parser.add_argument("--build-input", type=Path, help="compile-time input (default: RUN/i3.conf, else input/i3.conf)")
    parser.add_argument("--parents", type=int, default=DEFAULT_PARENTS, help="parent test particles per source")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--csv", type=Path, help="write b, AMPS kR, factor, corrected kR")
    args = parser.parse_args(argv)

    species = "OH" if re.search(r"pic\.OH\.brightness", args.brightness_file.name) else "H"
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from plot_i3_brightness_comparison import find_column, read_brightness

    variables, data, _ = read_brightness(args.brightness_file)
    b = data[:, find_column(variables, "Impact Parameter [m]", ("impact parameter", "[m]"))]
    kr = data[:, find_column(variables, "Brightness [kR]", ("brightness", "[kr]"))]
    run_dir = args.run_dir or run_directory(args.brightness_file)
    case = args.case or infer_case(args.brightness_file)
    factors, inputs = correction_for_run(
        run_dir, case, output_index(args.brightness_file), b, (species,), args.parents, args.seed,
        args.runtime_input, args.build_input)
    factor = factors[species]
    print(f"run {run_dir}  case {case}  species {species}  edge {inputs.outer_radius_m:.3e} m  "
          f"parents {args.parents}")
    print(f"{'b [m]':>12} {'AMPS [kR]':>12} {'factor':>8} {'corrected [kR]':>15}")
    for bb, k, f in zip(b, kr, factor):
        print(f"{bb:12.4e} {k:12.4e} {f:8.4f} {k * f:15.4e}" if math.isfinite(f) else f"{bb:12.4e} {k:12.4e} {'-':>8}")
    if args.csv:
        with args.csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["impact_parameter_m", "amps_kR", "edge_factor", "corrected_kR"])
            for bb, k, f in zip(b, kr, factor):
                writer.writerow([f"{bb:.6e}", f"{k:.6e}", f"{f:.6f}" if math.isfinite(f) else "",
                                 f"{k * f:.6e}" if math.isfinite(f) else ""])
        print(f"wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

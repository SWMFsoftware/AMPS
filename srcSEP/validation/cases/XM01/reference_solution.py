#!/usr/bin/env python3
"""Independent conservative finite-volume reference for XM01.

The solver advances cell-integrated probability on a periodic arclength grid
and a zero-flux pitch-angle grid. It does not import or call any srcSEP mover.
Periodic streaming is translated with a Fourier phase (spectrally accurate for
this smooth packet); deterministic and diffusive pitch fluxes use a
finite-volume discretization. Strang splitting and a small reviewed reference
timestep keep splitting and explicit diffusion errors below sampling error.
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np

PROTON_MASS_KG = 1.67262192369e-27
LIGHT_SPEED_M_PER_S = 299792458.0

def momentum_from_speed(speed: float) -> float:
    """Return relativistic proton momentum in kg m s^-1."""
    gamma = 1.0 / math.sqrt(1.0 - (speed / LIGHT_SPEED_M_PER_S) ** 2)
    return gamma * PROTON_MASS_KG * speed

def initial_probability(ns: int, nm: int, length: float) -> np.ndarray:
    """Integrate the reviewed Gaussian-by-linear initial density per cell."""
    sigma, center = 0.06 * length, 0.35 * length
    s_edges = np.linspace(0.0, length, ns + 1)
    spatial = np.array([
        0.5 * (math.erf((b-center)/(math.sqrt(2)*sigma)) -
               math.erf((a-center)/(math.sqrt(2)*sigma)))
        for a, b in zip(s_edges[:-1], s_edges[1:])])
    spatial /= spatial.sum()
    m_edges = np.linspace(-1.0, 1.0, nm + 1)
    angular = np.array([
        0.5 * ((b-a) + 0.3*(b*b-a*a))
        for a, b in zip(m_edges[:-1], m_edges[1:])])
    return spatial[:, None] * angular[None, :]

def advance(probability: np.ndarray, *, length: float, speed: float,
            d0: float, dlnb: float, duration: float, dt_requested: float) -> np.ndarray:
    """Advance the conservative Fokker-Planck flux form to ``duration``."""
    ns, nm = probability.shape
    ds, dm = length/ns, 2.0/nm
    mu = -1.0 + (np.arange(nm) + 0.5) * dm
    elapsed = 0.0
    value = probability.copy()
    wave_number = 2.0*math.pi*np.fft.fftfreq(ns, d=ds)
    def stream(state: np.ndarray, interval: float) -> np.ndarray:
        # Each pitch ordinate has constant characteristic speed v*mu. A Fourier
        # phase performs its periodic translation without the artificial
        # diffusion of a low-order upwind reference.
        transformed = np.fft.fft(state, axis=0)
        phase = np.exp(-1j*wave_number[:, None]*velocity[None, :]*interval)
        return np.fft.ifft(transformed*phase, axis=0).real
    velocity = speed * mu
    while elapsed < duration - 1.0e-14:
        dt = min(dt_requested, duration-elapsed)
        value = stream(value, 0.5*dt)

        # Pitch flux F=A_focus*f-D*df/dmu at interior faces. D vanishes at the
        # endpoints and both boundary fluxes are set exactly to zero.
        face = -1.0 + np.arange(nm + 1) * dm
        focus = -0.5 * (1.0-face*face) * speed * dlnb
        diffusion = d0 * np.maximum(0.0, 1.0-face*face)
        flux_m = np.zeros((ns, nm + 1))
        left, right = value[:, :-1], value[:, 1:]
        upwind = np.where(focus[1:-1][None, :] >= 0.0, left, right)
        flux_m[:, 1:-1] = focus[1:-1][None, :] * upwind / dm - \
            diffusion[1:-1][None, :] * (right-left) / (dm*dm)
        value -= dt * (flux_m[:, 1:] - flux_m[:, :-1])
        value = stream(value, 0.5*dt)
        # Fourier interpolation of a cell-averaged Gaussian can create a
        # roundoff-scale negative undershoot in cells whose exact mass is
        # effectively zero. Values below -1e-8 remain a hard failure; smaller
        # undershoots are projected to zero and the conservative total is
        # restored explicitly.
        if np.min(value) < -1.0e-8 or not np.all(np.isfinite(value)):
            raise RuntimeError("XM01 finite-volume reference lost positivity")
        value = np.maximum(value, 0.0)
        value /= value.sum()
        elapsed += dt
    return value

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    case = json.loads(args.input.read_text(encoding="utf-8"))
    p, n = case["physics"], case["numerics"]
    scenarios = {
        "streaming": (0.0, 0.0, 0.0),
        "scattering": (p["d0_per_s"], 0.0, 0.0),
        "focusing": (0.0, p["dlnb_ds_per_m"], 0.0),
        "adiabatic": (0.0, 0.0, p["divergence_per_s"]),
        "combined": (p["d0_per_s"], p["dlnb_ds_per_m"], p["divergence_per_s"]),
    }
    initial = initial_probability(n["s_bins"], n["mu_bins"], p["length_m"])
    p0 = momentum_from_speed(p["speed_m_per_s"])
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["scenario", "s_bin", "mu_bin", "s_left_m", "s_right_m",
                         "mu_left", "mu_right", "probability", "intensity",
                         "anisotropy", "mean_log_p"])
        for name, (d0, dlnb, divergence) in scenarios.items():
            solution = advance(initial, length=p["length_m"], speed=p["speed_m_per_s"],
                               d0=d0, dlnb=dlnb, duration=n["duration_s"],
                               dt_requested=n["reference_time_step_s"])
            mu = -1.0 + (np.arange(n["mu_bins"]) + 0.5) * 2.0/n["mu_bins"]
            anisotropy = float(3.0 * np.sum(solution * mu[None, :]))
            mean_log_p = math.log(p0) - divergence*n["duration_s"]/3.0
            intensity = solution.sum(axis=1)
            for i in range(n["s_bins"]):
                for j in range(n["mu_bins"]):
                    writer.writerow([name, i, j, p["length_m"]*i/n["s_bins"],
                        p["length_m"]*(i+1)/n["s_bins"], -1+2*j/n["mu_bins"],
                        -1+2*(j+1)/n["mu_bins"], solution[i, j], intensity[i],
                        anisotropy, mean_log_p])
    temporary.replace(args.output)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

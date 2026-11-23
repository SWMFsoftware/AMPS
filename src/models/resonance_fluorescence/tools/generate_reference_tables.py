#!/usr/bin/env python3
"""Generate immutable C++ fluorescence tables from reviewed reference data.

The source text/CSV files remain the review and provenance surface.  AMPS
does not yet have a stable installed reference-data root, so production code
uses this generated include and never opens a data file in a particle or
sampling loop.  ``--check`` verifies that the checked-in include is current
without rewriting it.

The atomic-oxygen O I 130.4-nm tables are derived here, deterministically and
with the standard library only, from three reviewed sources (see
``reference/g_factors/manifest.json`` and the module README):

* Killen, Vervack & Burger (2022) data behind their Figure 14: g values of
  the three O I 2p4 3P_J - 3s 3S1 lines, computed at 0.352 au from Gladstone
  (1992) solar line shapes scaled to SORCE SOLSTICE on 2009-03-01;
* the SOHO/SUMER quiet-Sun disk-centre atlas (Curdt et al. 2001), retrieved
  from BASS2000, around the O I triplet and around H Lyman-beta;
* NIST ASD oscillator strengths and transition probabilities.

Every O I model assumes that the absorbing atoms are in the 2p4 3P2 ground
level, i.e. that the fine-structure levels have relaxed radiatively (the
3P1 and 3P0 lifetimes are 3.1 h and 16 h), and returns the rate of 130.4-nm
triplet photons (130.2 + 130.5 + 130.6 nm) per atom.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import math
from pathlib import Path
import tempfile


MODULE = Path(__file__).resolve().parents[1]
REFERENCE = MODULE / "reference/g_factors"
H_SOURCE = REFERENCE / "h_lya_ikl2020/table1v4.txt"
OH_SOURCE = REFERENCE / (
    "oh_ax00_schleicher_ahearn1988/lowell_oh00_velocity_gfactor_1au.csv"
)
O_KVB_ARCHIVE = REFERENCE / "o_i_1304_kvb2022/apjsac9eabf14.tar.gz"
O_KVB_TABLES = {
    "1302": REFERENCE / "o_i_1304_kvb2022/oI_1302_gvalues.txt",
    "1304": REFERENCE / "o_i_1304_kvb2022/oI_1304_gvalues.txt",
    "1306": REFERENCE / "o_i_1304_kvb2022/oI_1306_gvalues.txt",
}
O_SUMER_TRIPLET = REFERENCE / "o_i_sumer2001_bass2000/sumer_disk_center_1297_1311A.txt"
O_SUMER_LYB = REFERENCE / "o_i_sumer2001_bass2000/sumer_disk_center_1019_1032A.txt"
O_NIST_UV = REFERENCE / "o_i_nist_asd/o_i_lines_1020_1310A_vacuum.csv"
O_NIST_IR = REFERENCE / "o_i_nist_asd/o_i_lines_11270_11310A_air.csv"
OUTPUT = MODULE / "generated/reference_tables.inc"
PARAMETERS = ("AK", "muK", "sigmaK", "kappa", "AR", "dmu", "sigmaR", "bbkg", "abkg")

# CODATA 2018 (c, h exact).  The classical integrated line strength
# pi e^2/(m_e c) is evaluated in CGS because the published g-value formula
# (Killen et al. 2022, Eq. 1) is written in CGS.
LIGHT_SPEED_KM_S = 299792.458
LIGHT_SPEED_CM_S = 2.99792458e10
LIGHT_SPEED_M_S = 299792458.0
PLANCK_J_S = 6.62607015e-34
ELECTRON_CHARGE_ESU = 4.803204712570263e-10
ELECTRON_MASS_G = 9.1093837015e-28
PI_E2_MC_CM2_S = math.pi * ELECTRON_CHARGE_ESU**2 / (ELECTRON_MASS_G * LIGHT_SPEED_CM_S)

# Killen et al. (2022) evaluate all g values at Mercury's 0.352 au.  The
# generated Gladstone profiles end near +-48.3 km/s (piecewise-linear
# artefacts and exact zeros follow), so the reviewed domain stops at 48 km/s.
KVB_DISTANCE_AU = 0.352
KVB_LIMIT_KM_S = 48.0
SUMER_LIMIT_KM_S = 100.0
VELOCITY_STEP_KM_S = 0.1

# Atlas normalisation windows.  The triplet window holds the complete
# 130.217-nm line; the continuum is the line-free 1303.55-1304.05 A interval.
# The Lyman-beta window is identical to the four FISM2 0.1-nm bins
# 102.4-102.8 nm that the I3 case generator integrates.
TRIPLET_HALF_WINDOW_A = 0.45
TRIPLET_CONTINUUM_A = (1303.55, 1304.05)
LYB_WINDOW_A = (1024.0, 1028.0)

GROUND = ("2s2.2p4", "3P")
UPPER_3S = ("2s2.2p3.(4S*).3s", "3S*", "1")
UPPER_3D = ("2s2.2p3.(4S*).3d", "3D*")
MIDDLE_3P = ("2s2.2p3.(4S*).3p", "3P")

# Every float summation below uses math.fsum, the correctly rounded sum.
# Python 3.12 changed the built-in sum() of floats from left-to-right
# addition to compensated summation.  With sum(), the emitted %.17g values
# (and so the --check comparison) depended on the interpreter: the O I
# tables differed by 1-3 units in the last place between Python 3.8 and
# 3.14.  math.fsum gives the same bits in every Python version; on
# 2026-11-23 it reproduced the checked-in include byte for byte under
# Python 3.8.10 and 3.14.10.


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_h() -> list[tuple[str, float, float]]:
    rows: dict[str, tuple[float, float]] = {}
    for line in H_SOURCE.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) == 3 and fields[0] in PARAMETERS:
            rows[fields[0]] = (float(fields[1]), float(fields[2]))
    missing = set(PARAMETERS) - rows.keys()
    if missing:
        raise RuntimeError(f"missing IKL parameters: {sorted(missing)}")
    return [(name, *rows[name]) for name in PARAMETERS]


def read_oh() -> tuple[list[float], list[float]]:
    velocity: list[float] = []
    energy_g: list[float] = []
    with OH_SOURCE.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["api_extrapolated"] != "false":
                raise RuntimeError("an extrapolated Lowell value cannot be embedded")
            velocity.append(float(row["heliocentric_radial_velocity_km_s"]))
            energy_g.append(float(row["oh00_energy_g_1au_erg_s-1_molecule-1"]))
    if len(velocity) != 1201 or velocity[0] != -60.0 or velocity[-1] != 60.0:
        raise RuntimeError("Lowell OH table must contain -60..60 km/s at 0.1 km/s")
    for left, right in zip(velocity, velocity[1:]):
        if abs((right - left) - 0.1) > 1.0e-10:
            raise RuntimeError("Lowell OH velocity grid is not uniformly spaced")
    if any(value <= 0.0 for value in energy_g):
        raise RuntimeError("Lowell OH energy g factors must be positive")
    return velocity, energy_g


# ---------------------------------------------------------------- O I data
def velocity_grid(limit_km_s: float) -> list[float]:
    count = int(round(2.0 * limit_km_s / VELOCITY_STEP_KM_S)) + 1
    return [-limit_km_s + index * VELOCITY_STEP_KM_S for index in range(count)]


def interpolate(xs: list[float], ys: list[float], x: float) -> float:
    if not xs[0] <= x <= xs[-1]:
        raise RuntimeError(f"interpolation abscissa {x} outside [{xs[0]}, {xs[-1]}]")
    upper = bisect.bisect_right(xs, x)
    if upper >= len(xs):
        return ys[-1]
    lower = upper - 1
    fraction = (x - xs[lower]) / (xs[upper] - xs[lower])
    return ys[lower] + fraction * (ys[upper] - ys[lower])


def integrate_linear(xs: list[float], ys: list[float], left: float, right: float) -> float:
    """Exact integral of the piecewise-linear interpolant on [left, right]."""
    knots = [left] + [x for x in xs if left < x < right] + [right]
    values = [interpolate(xs, ys, x) for x in knots]
    return math.fsum(
        0.5 * (values[i] + values[i + 1]) * (knots[i + 1] - knots[i])
        for i in range(len(knots) - 1)
    )


def line_coefficient(wavelength_a: float, oscillator_strength: float) -> float:
    """g = coefficient * F_lambda, F_lambda in photons cm^-2 s^-1 A^-1."""
    return (
        PI_E2_MC_CM2_S
        * oscillator_strength
        * wavelength_a
        * (wavelength_a * 1.0e-8)
        / LIGHT_SPEED_CM_S
    )


def read_kvb(name: str) -> tuple[list[float], list[float]]:
    lines = O_KVB_TABLES[name].read_text(encoding="ascii").splitlines()
    if "g_value" not in lines[0] or not lines[1].startswith("Velocity (km/s)"):
        raise RuntimeError(f"unexpected Killen et al. header in {O_KVB_TABLES[name]}")
    velocity: list[float] = []
    g_value: list[float] = []
    for line in lines[2:]:
        if line.strip():
            left, right = line.split()
            velocity.append(float(left))
            g_value.append(float(right))
    if len(velocity) != 1001 or velocity[0] != -50.0 or velocity[-1] != 50.0:
        raise RuntimeError("Killen et al. tables must contain -50..50 km/s")
    for left, right in zip(velocity, velocity[1:]):
        if abs((right - left) - 0.1) > 1.0e-9:
            raise RuntimeError("Killen et al. velocity grid is not uniform")
    if any(value < 0.0 for value in g_value):
        raise RuntimeError("Killen et al. g values must be non-negative")
    return velocity, g_value


def nist_rows(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        for raw in reader:
            if raw:
                rows.append(
                    {key: value.lstrip("=").strip('"') for key, value in zip(header, raw)}
                )
    return rows


def nist_line(rows, lower, lower_j, upper, upper_j) -> dict[str, str]:
    found = [
        row
        for row in rows
        if (row["conf_i"], row["term_i"], row["J_i"]) == (*lower, lower_j)
        and (row["conf_k"], row["term_k"], row["J_k"]) == (*upper, upper_j)
    ]
    if len(found) != 1:
        raise RuntimeError(
            f"expected one NIST line {lower} J={lower_j} - {upper} J={upper_j}, found {len(found)}"
        )
    return found[0]


def nist_ritz_vacuum_a(row: dict[str, str]) -> float:
    return float(row["ritz_wl_vac(A)"])


def oxygen_atomic_data() -> dict[str, object]:
    uv = nist_rows(O_NIST_UV)
    ir = nist_rows(O_NIST_IR)
    triplet = []
    for lower_j in ("2", "1", "0"):
        row = nist_line(uv, GROUND, lower_j, UPPER_3S[:2], UPPER_3S[2])
        triplet.append(
            (lower_j, nist_ritz_vacuum_a(row), float(row["fik"]), float(row["Aki(s^-1)"]))
        )
    # 3d 3D*_J' populated from 3P2; branching into 3p 3P (11287 A), which
    # cascades through 8446 A into the 3s 3S1 level and so yields exactly
    # one 130.4-nm triplet photon.  The other channel returns directly to
    # 2p4 3P_J at 1025.8-1028.2 nm.
    bowen = []
    for upper_j in ("3", "2", "1"):
        pump = nist_line(uv, GROUND, "2", UPPER_3D, upper_j)
        direct = [
            float(row["Aki(s^-1)"])
            for row in uv
            if (row["conf_i"], row["term_i"]) == GROUND
            and (row["conf_k"], row["term_k"], row["J_k"]) == (*UPPER_3D, upper_j)
        ]
        cascade = [
            float(row["Aki(s^-1)"])
            for row in ir
            if (row["conf_i"], row["term_i"]) == MIDDLE_3P
            and (row["conf_k"], row["term_k"], row["J_k"]) == (*UPPER_3D, upper_j)
        ]
        if not direct or not cascade:
            raise RuntimeError(f"incomplete NIST decay data for 3d 3D J={upper_j}")
        branching = math.fsum(cascade) / (math.fsum(cascade) + math.fsum(direct))
        bowen.append(
            (upper_j, nist_ritz_vacuum_a(pump), float(pump["fik"]), branching)
        )
    return {"triplet": triplet, "bowen": bowen}


def read_atlas(path: Path) -> tuple[list[float], list[float]]:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines[0].startswith("Wavelength (Angstroms), Normalized intensity"):
        raise RuntimeError(f"unexpected BASS2000 header in {path}")
    pairs = []
    for line in lines[1:]:
        if line.strip():
            fields = line.split(",")
            pairs.append((float(fields[0]), float(fields[1])))
    pairs.sort()
    wavelength = [pair[0] for pair in pairs]
    intensity = [pair[1] for pair in pairs]
    if any(right <= left for left, right in zip(wavelength, wavelength[1:])):
        raise RuntimeError(f"duplicate atlas wavelengths in {path}")
    if any(value <= 0.0 for value in intensity):
        raise RuntimeError(f"atlas intensities must be positive in {path}")
    return wavelength, intensity


def require_uniform_atlas(wavelength: list[float], left: float, right: float) -> None:
    # BASS2000 stitches atlas sections; a seam inside a used window would
    # silently join two exposures.  The reviewed windows contain none.
    inside = [x for x in wavelength if left <= x <= right]
    spacing = [b - a for a, b in zip(inside, inside[1:])]
    if not inside or min(spacing) < 0.03 or max(spacing) > 0.05:
        raise RuntimeError(f"atlas window {left}-{right} A is not a single uniform section")


def oxygen_tables() -> dict[str, object]:
    atomic = oxygen_atomic_data()
    triplet = atomic["triplet"]
    bowen = atomic["bowen"]
    j2_wavelength, j2_f = triplet[0][1], triplet[0][2]

    # Killen et al.: recover each line's implied solar photon flux at 1 AU
    # from g = coefficient * F_lambda, F_lambda sampled at lambda0/(1+v/c).
    kvb_velocity, kvb_1302 = read_kvb("1302")
    line_flux = {}
    for (name, row) in zip(("1302", "1304", "1306"), triplet):
        velocity, g_value = read_kvb(name)
        coefficient = line_coefficient(row[1], row[2])
        flux_lambda = [value / coefficient for value in g_value]
        dlambda_dv = row[1] / LIGHT_SPEED_KM_S
        line_flux[name] = (
            integrate_linear(velocity, flux_lambda, velocity[0], velocity[-1])
            * dlambda_dv
            * KVB_DISTANCE_AU**2
        )
    reference_flux_cm2 = math.fsum(line_flux.values())
    fraction_1302 = line_flux["1302"] / reference_flux_cm2
    kvb_grid = velocity_grid(KVB_LIMIT_KM_S)
    kvb_rate = [
        interpolate(kvb_velocity, kvb_1302, v) * KVB_DISTANCE_AU**2 for v in kvb_grid
    ]
    if any(value <= 0.0 for value in kvb_rate):
        raise RuntimeError("Killen et al. 130.2-nm g values vanish inside the domain")

    # SUMER direct pumping of 3P2 at 130.217 nm, normalised so that the
    # continuum-subtracted 130.217-nm line carries fraction_1302 of the
    # 1-AU triplet photon flux supplied at run time.
    wavelength, intensity = read_atlas(O_SUMER_TRIPLET)
    window = (j2_wavelength - TRIPLET_HALF_WINDOW_A, j2_wavelength + TRIPLET_HALF_WINDOW_A)
    sumer_limit = j2_wavelength * SUMER_LIMIT_KM_S / LIGHT_SPEED_KM_S
    require_uniform_atlas(
        wavelength,
        min(window[0], j2_wavelength - 1.01 * sumer_limit),
        max(window[1], j2_wavelength + 1.01 * sumer_limit),
    )
    continuum_samples = sorted(
        value
        for x, value in zip(wavelength, intensity)
        if TRIPLET_CONTINUUM_A[0] <= x <= TRIPLET_CONTINUUM_A[1]
    )
    middle = len(continuum_samples) // 2
    continuum = (
        continuum_samples[middle]
        if len(continuum_samples) % 2
        else 0.5 * (continuum_samples[middle - 1] + continuum_samples[middle])
    )
    line_integral = integrate_linear(
        wavelength, [value - continuum for value in intensity], *window
    )
    sumer_grid = velocity_grid(SUMER_LIMIT_KM_S)
    # Per unit 1-AU triplet photon flux in photons m^-2 s^-1 (factor 1e-4
    # converts the SI flux into the CGS flux of the line coefficient).
    direct_scale = line_coefficient(j2_wavelength, j2_f) * fraction_1302 / line_integral * 1.0e-4
    sumer_direct = [
        direct_scale
        * interpolate(wavelength, intensity, j2_wavelength / (1.0 + v / LIGHT_SPEED_KM_S))
        for v in sumer_grid
    ]

    # Lyman-beta (Bowen) pumping of 3P2 -> 3d 3D*, followed by the
    # 11287 A -> 8446 A cascade into the 130.4-nm triplet.
    lyb_wavelength, lyb_intensity = read_atlas(O_SUMER_LYB)
    require_uniform_atlas(lyb_wavelength, *LYB_WINDOW_A)
    lyb_integral = integrate_linear(lyb_wavelength, lyb_intensity, *LYB_WINDOW_A)
    bowen_rate = []
    for v in sumer_grid:
        total = 0.0
        for _, pump_wavelength, pump_f, branching in bowen:
            total += (
                line_coefficient(pump_wavelength, pump_f)
                * branching
                * interpolate(
                    lyb_wavelength,
                    lyb_intensity,
                    pump_wavelength / (1.0 + v / LIGHT_SPEED_KM_S),
                )
            )
        bowen_rate.append(total / lyb_integral * 1.0e-4)

    # Mean energy of an emitted triplet photon: the 3s 3S1 level decays to
    # 3P2, 3P1, 3P0 in proportion to the three Einstein A values.
    total_a = math.fsum(row[3] for row in triplet)
    mean_energy = math.fsum(
        row[3] / total_a * PLANCK_J_S * LIGHT_SPEED_M_S / (row[1] * 1.0e-10)
        for row in triplet
    )
    return {
        "kvb_grid": kvb_grid,
        "kvb_rate": kvb_rate,
        "reference_flux_m2": reference_flux_cm2 * 1.0e4,
        "line_flux_cm2": line_flux,
        "fraction_1302": fraction_1302,
        "sumer_grid": sumer_grid,
        "sumer_direct": sumer_direct,
        "bowen_rate": bowen_rate,
        "continuum": continuum,
        "line_integral": line_integral,
        "lyb_integral": lyb_integral,
        "mean_energy": mean_energy,
        "triplet": triplet,
        "bowen": bowen,
    }


def emit_array(lines: list[str], name: str, values: list[float]) -> None:
    lines.append(f"constexpr double {name}[] = {{")
    for start in range(0, len(values), 6):
        chunk = ",".join(f"{value:.17g}" for value in values[start : start + 6])
        lines.append(f"  {chunk},")
    lines.append("};")


def render() -> str:
    h_rows = read_h()
    velocity, energy_g = read_oh()
    oxygen = oxygen_tables()
    lines = [
        "// Generated by tools/generate_reference_tables.py; do not edit.",
        f'constexpr char HReferenceSha256[] = "{sha256(H_SOURCE)}";',
        f'constexpr char OHReferenceSha256[] = "{sha256(OH_SOURCE)}";',
        "constexpr IKLCoefficient HIKL2020Coefficients[] = {",
    ]
    for name, beta, alpha in h_rows:
        lines.append(f'  {{"{name}",{beta:.17g},{alpha:.17g}}},')
    lines.extend(
        [
            "};",
            f"constexpr std::size_t OHTableLength = {len(energy_g)};",
            f"constexpr double OHVelocityMinimumMPS = {velocity[0] * 1.0e3:.17g};",
            # The reviewed grid is exactly 0.1 km/s.  Emit the exact SI value
            # after validating every source abscissa above; subtracting two
            # binary floats would otherwise encode 100.00000000000142 m/s.
            "constexpr double OHVelocityStepMPS = 100.0;",
            "constexpr double OHEnergyGFactor1AUErgPerSecond[] = {",
        ]
    )
    for start in range(0, len(energy_g), 6):
        values = ",".join(f"{value:.17g}" for value in energy_g[start : start + 6])
        lines.append(f"  {values},")
    lines.append("};")

    line_flux = oxygen["line_flux_cm2"]
    lines.extend(
        [
            "",
            "// O I 130.4-nm triplet models; absorbers in 2p4 3P2.",
            f'constexpr char OKVB2022ReferenceSha256[] = "{sha256(O_KVB_ARCHIVE)}";',
            f'constexpr char OSUMER2001TripletReferenceSha256[] = "{sha256(O_SUMER_TRIPLET)}";',
            f'constexpr char OSUMER2001LymanBetaReferenceSha256[] = "{sha256(O_SUMER_LYB)}";',
            f'constexpr char ONISTUltravioletReferenceSha256[] = "{sha256(O_NIST_UV)}";',
            f'constexpr char ONISTInfraredReferenceSha256[] = "{sha256(O_NIST_IR)}";',
            "// Killen et al. (2022) 1-AU solar line photon fluxes implied by their",
            "// g values and the NIST oscillator strengths [photons cm^-2 s^-1]:",
            f"//   130.217 nm {line_flux['1302']:.6e}, 130.486 nm {line_flux['1304']:.6e},"
            f" 130.603 nm {line_flux['1306']:.6e}; 130.217-nm fraction"
            f" {oxygen['fraction_1302']:.6f}.",
            "// SUMER atlas: 130.217-nm continuum-subtracted line integral"
            f" {oxygen['line_integral']:.6f} atlas*A over +-{TRIPLET_HALF_WINDOW_A} A,"
            f" continuum {oxygen['continuum']:.6f} atlas units;",
            f"// Lyman-beta window {LYB_WINDOW_A[0]}-{LYB_WINDOW_A[1]} A integral"
            f" {oxygen['lyb_integral']:.6f} atlas*A.",
            "// Bowen pump lines (3P2 -> 3d 3D*_J', vacuum A, f, 11287-A branching):",
        ]
    )
    for upper_j, wavelength, oscillator, branching in oxygen["bowen"]:
        lines.append(
            f"//   J'={upper_j} {wavelength:.5f} A f={oscillator:.3e} branching={branching:.6f}"
        )
    lines.extend(
        [
            "constexpr double OKVB2022ReferenceTripletFlux1AUPhotonM2S = "
            f"{oxygen['reference_flux_m2']:.17g};",
            f"constexpr double OTripletMeanPhotonEnergyJ = {oxygen['mean_energy']:.17g};",
            f"constexpr std::size_t OKVB2022TableLength = {len(oxygen['kvb_rate'])};",
            f"constexpr double OKVB2022VelocityMinimumMPS = {-KVB_LIMIT_KM_S * 1.0e3:.17g};",
            "constexpr double OKVB2022VelocityStepMPS = 100.0;",
        ]
    )
    emit_array(lines, "OKVB2022DirectPhotonRate1AUS1", oxygen["kvb_rate"])
    lines.extend(
        [
            f"constexpr std::size_t OSUMER2001TableLength = {len(oxygen['sumer_direct'])};",
            f"constexpr double OSUMER2001VelocityMinimumMPS = {-SUMER_LIMIT_KM_S * 1.0e3:.17g};",
            "constexpr double OSUMER2001VelocityStepMPS = 100.0;",
            "// Photon rates per unit 1-AU solar photon flux [s^-1 / (photons m^-2 s^-1)].",
        ]
    )
    emit_array(lines, "OSUMER2001DirectRatePerTripletFlux1AU", oxygen["sumer_direct"])
    emit_array(lines, "OSUMER2001BowenRatePerLymanBetaFlux1AU", oxygen["bowen_rate"])
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = render()
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != content:
            raise SystemExit(f"generated table is stale: {OUTPUT}")
        return 0
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    # Atomic replacement prevents an interrupted generator from leaving a
    # compilable-looking partial array.
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=OUTPUT.parent, delete=False
    ) as stream:
        stream.write(content)
        temporary = Path(stream.name)
    temporary.replace(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

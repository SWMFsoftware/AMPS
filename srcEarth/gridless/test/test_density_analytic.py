#!/usr/bin/env python3
"""
test_density_analytic.py
========================
Compare AMPS gridless density/spectrum output against analytic predictions
for four test cases in a pure DIPOLE magnetic field.

USAGE
-----
  # Run after AMPS has produced its output files in the current directory:
  python3 test_density_analytic.py [--spectrum-dir <path>] [--tol <fraction>]

  Default: looks for gridless_points_density.dat and gridless_points_spectrum.dat
  in the current directory.

  Options:
    --spectrum-dir PATH   directory containing the AMPS output files
    --tol FRACTION        relative tolerance for pass/fail (default 0.05 = 5 %)
    --no-plots            skip matplotlib figures even if matplotlib is available

PHYSICS SUMMARY
---------------
The AMPS gridless solver computes, at each observation point x0:

  J_loc(E; x0) = T(E; x0) * J_b(E)

  n_tot(x0) = 4pi * integral[Emin,Emax] T(E; x0) * J_b(E) / v(E) dE

where
  J_b(E)  = boundary differential intensity  [m^-2 s^-1 sr^-1 MeV^-1]
  T(E; x0) = transmissivity = fraction of isotropic directions that
             backtrack to the outer domain box (= ALLOWED trajectories)
  v(E)    = relativistic proton speed at kinetic energy E

For a power-law boundary spectrum  J_b(E) = J0 * (E/E0)^{-gamma}
and observation points where ALL energies are above the maximum Størmer
cutoff, the transmissivity equals the geometric transmissivity:

  T(E; x0) = T_geo(r)  for all E in [Emin, Emax]

  T_geo(r) = (1 + sqrt(1 - (r_inner/r)^2)) / 2

This is the EXACT analytic prediction — no free parameters.

DIPOLE STØRMER CUTOFFS (for protons, Earth dipole B_eq = 3.12e-5 T)
--------------------------------------------------------------------
The maximum Størmer cutoff at equatorial radius r (in Earth radii):

  Rc_max(r) = (mu0/4pi * M_E / Re^2) * (c/1e9) / r^2   [GV]
            = 59.6 / r^2  GV                             [Re units]

At the geomagnetic pole: Rc = 0 for all r (cos(lambda)=0 exactly).

TEST CASES
----------
Point index  Location (GSM km)       Energy range      Expected T       Test type
-----------  ---------------------   ---------------   ---------------  ---------
0  TC1        (50969.6, 0, 0)        2–20 GeV          T_geo = 0.99600  EXACT
1  TC2        (0, 0, 31856)          50–500 MeV        T_geo = 0.98969  EXACT
2  TC3        (38227.2, 0, 0)        3–30 GeV          T_geo = 0.99287  EXACT
3  TC4        (12742.4, 0, 0)        10–500 MeV        T << T_geo       UPPER BOUND

For TC1–TC3, the pass criterion is:
  |n_model - n_analytic| / n_analytic < tolerance

For TC4, the pass criterion is:
  n_model < 0.5 * n_analytic_Tgeo   (strong magnetic shielding)

EFFECT OF EARTH'S PRESENCE ON THE ANALYTIC ESTIMATE
----------------------------------------------------
Earth's inner boundary (loss sphere, r < r_inner = 1.01 Re) affects the
density in TWO ways:

1. GEOMETRIC SHADOW (always present):
   From any point at radius r, Earth subtends a half-angle
     alpha = arcsin(r_inner / r)
   The fraction of 4pi solid angle blocked is
     Omega_blocked / 4pi = (1 - cos alpha) / 2
   So T_geo = (1 + cos alpha) / 2 = (1 + sqrt(1 - (r_inner/r)^2)) / 2

   This is the factor encoded in T_geo and IS captured by the solver
   automatically (trajectories hitting the inner sphere are FORBIDDEN).

   At r = 2 Re: T_geo = 0.932  (7% blocked by Earth disk)
   At r = 5 Re: T_geo = 0.990  (1% blocked)
   At r = 8 Re: T_geo = 0.996  (0.4% blocked)

2. GEOMAGNETIC SHIELDING (additional, energy-dependent):
   For energies below the Størmer cutoff, the geomagnetic field additionally
   forbids trajectories that would otherwise escape. This is an extra factor
   T_Stormer(E; x0) <= 1 on top of the geometric shadow.

   For TC1–TC3 (energies above all cutoffs): T_Stormer = 1.
   For TC4 (energies below cutoffs): T_Stormer << 1  →  n << n_free * T_geo.

The total transmissivity is:
  T_total = T_geo * T_Stormer   (schematically; the solver computes this jointly)

TOLERANCE JUSTIFICATION
-----------------------
The solver uses a finite direction sample (typically 24 azimuths × 24 elevations
= 576 directions per energy per point). The statistical uncertainty in T from
sampling is approximately 1/sqrt(N_dirs) ~ 4%. We therefore use 5% tolerance
for TC1–TC3. TC4 uses a factor-of-2 margin (50% of the upper bound) which is
very conservative given the deep shielding expected.
"""

import sys
import os
import math
import argparse
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Physical constants (must match DipoleInterface.h and DensityGridless.cpp)
# ---------------------------------------------------------------------------
C_LIGHT  = 299792458.0           # m/s
QE       = 1.602176634e-19       # C
AMU      = 1.66053906660e-27     # kg
MP_KG    = 1.007276466621 * AMU  # proton mass kg
MP_MEV   = MP_KG * C_LIGHT**2 / (1e6 * QE)  # ~938.272 MeV
RE_M     = 6371.2e3              # m   (must match kEarthRadiusKm in amps_param_parser)
RE_KM    = 6371.2                # km
MU0_4PI  = 1e-7                  # T·m/A
B_EQ_RE  = 3.12e-5               # T   (from DipoleInterface.h)
M_E      = B_EQ_RE * RE_M**3 / MU0_4PI  # A·m^2
CS_GV    = MU0_4PI * M_E / RE_M**2 * C_LIGHT / 1e9  # ~59.6 GV·Re²

# ---------------------------------------------------------------------------
# Spectrum parameters (MUST match the #SPECTRUM section in the .in file)
# ---------------------------------------------------------------------------
J0_SPEC  = 1.0e4    # m^-2 s^-1 sr^-1 MeV^-1
E0_SPEC  = 100.0    # MeV
GAMMA    = 2.0      # spectral index

# Inner loss sphere radius (must match R_INNER in the .in file)
R_INNER_KM = 6434.9   # km = 1.01 Re

# ---------------------------------------------------------------------------
# Test case definitions
# ---------------------------------------------------------------------------
# Each entry: (label, x_km, y_km, z_km, Emin_test, Emax_test, r_Re_approx)
#
#  Emin_test, Emax_test — the energy sub-range used for the analytic comparison
#  (the full AMPS energy grid is wider; we integrate only over the sub-range
#   where T = T_geo is guaranteed, i.e., above Rc_max for that point)
TEST_CASES = [
    {
        "label":      "TC1 — 8Re equator (above-cutoff, free-space)",
        "x_km":       8.0 * RE_KM,
        "y_km":       0.0,
        "z_km":       0.0,
        "Emin_test":  2000.0,   # MeV  (well above Rc_max(8Re)=384 MeV)
        "Emax_test":  20000.0,  # MeV
        "r_Re":       8.0,
        "lam_deg":    0.0,      # geomagnetic latitude (equator)
        "test_type":  "exact",
        "description": (
            "At r=8 Re equator, Rc_max = 0.93 GV → 384 MeV. "
            "All energies in 2–20 GeV are above this, so T = T_geo exactly. "
            "Tests the free-space (no-field) density integral."
        ),
    },
    {
        "label":      "TC2 — 5Re pole (geomagnetic pole, Rc=0)",
        "x_km":       0.0,
        "y_km":       0.0,
        "z_km":       5.0 * RE_KM,
        "Emin_test":  50.0,    # MeV
        "Emax_test":  500.0,   # MeV
        "r_Re":       5.0,
        "lam_deg":    90.0,    # geomagnetic pole: cos(lambda) = 0
        "test_type":  "exact",
        "description": (
            "At the geomagnetic pole, cos(lambda)=0 makes Størmer Rc=0 for ALL "
            "energies and ALL directions. So T = T_geo regardless of energy. "
            "Tests that the solver does NOT over-shield at high geomagnetic latitudes."
        ),
    },
    {
        "label":      "TC3 — 6Re equator (above-cutoff, mid-distance)",
        "x_km":       6.0 * RE_KM,
        "y_km":       0.0,
        "z_km":       0.0,
        "Emin_test":  3000.0,   # MeV  (Rc_max(6Re)=1.66 GV → 965 MeV)
        "Emax_test":  30000.0,  # MeV
        "r_Re":       6.0,
        "lam_deg":    0.0,
        "test_type":  "exact",
        "description": (
            "At r=6 Re equator, Rc_max = 1.66 GV → 965 MeV. "
            "All energies in 3–30 GeV are above this → T = T_geo. "
            "Closer-in than TC1; tests that the geometric shadow scales correctly."
        ),
    },
    {
        "label":      "TC4 — 2Re equator (deep Størmer shielding)",
        "x_km":       2.0 * RE_KM,
        "y_km":       0.0,
        "z_km":       0.0,
        "Emin_test":  10.0,    # MeV
        "Emax_test":  500.0,   # MeV
        "r_Re":       2.0,
        "lam_deg":    0.0,
        "test_type":  "upper_bound",
        "description": (
            "At r=2 Re equator, Rc_vert = 3.73 GV → 2.9 GeV. "
            "All test energies (up to 500 MeV, R=1.09 GV) are far below Rc_vert. "
            "Strong magnetic shielding expected: n << T_geo * n_free. "
            "Upper-bound test only — exact T_Stormer not analytically tractable."
        ),
    },
]


# ---------------------------------------------------------------------------
# Analytic helpers
# ---------------------------------------------------------------------------

def v_from_E_MeV(E_MeV: float) -> float:
    """Relativistic proton speed [m/s] from kinetic energy [MeV]."""
    gamma = 1.0 + E_MeV / MP_MEV
    beta2 = max(0.0, 1.0 - 1.0 / gamma**2)
    return C_LIGHT * math.sqrt(beta2)


def Jb_from_E_MeV(E_MeV: float) -> float:
    """Power-law boundary spectrum [m^-2 s^-1 sr^-1 MeV^-1]."""
    return J0_SPEC * (E_MeV / E0_SPEC) ** (-GAMMA)


def T_geometric(r_Re: float, r_inner_km: float = R_INNER_KM) -> float:
    """
    Geometric transmissivity: fraction of 4pi NOT blocked by the inner sphere.

    T_geo = (1 + cos(alpha)) / 2
    where alpha = arcsin(r_inner / r) is the half-angle of the blocked cone.
    """
    r_km = r_Re * RE_KM
    sin_alpha = r_inner_km / r_km
    if sin_alpha >= 1.0:
        return 0.0
    cos_alpha = math.sqrt(1.0 - sin_alpha**2)
    return (1.0 + cos_alpha) / 2.0


def Rc_max_GV(r_Re: float, lam_deg: float) -> float:
    """Maximum Størmer cutoff [GV] at (r_Re, lam_deg)."""
    cos_lam = math.cos(math.radians(lam_deg))
    return CS_GV * cos_lam**4 / r_Re**2


def E_kin_from_R_GV(R_GV: float) -> float:
    """Kinetic energy [MeV] corresponding to rigidity R_GV [GV] for a proton."""
    pc_MeV = R_GV * 1e3   # R [GV] × 1000 = pc [MeV] for Z=1
    return math.sqrt(pc_MeV**2 + MP_MEV**2) - MP_MEV


def n_integral_analytic(Emin_MeV: float, Emax_MeV: float, T: float,
                        N: int = 4000) -> float:
    """
    Compute analytic density [m^-3] for a given T (constant over energy):

      n = 4pi * T * integral[Emin,Emax] J_b(E) / v(E) dE

    Uses log-spaced quadrature with N points.
    """
    Es = np.logspace(np.log10(Emin_MeV), np.log10(Emax_MeV), N)
    Jb = J0_SPEC * (Es / E0_SPEC) ** (-GAMMA)   # m^-2 s^-1 sr^-1 MeV^-1
    vs = np.array([v_from_E_MeV(e) for e in Es], dtype=float)  # m/s
    integrand = Jb / vs  # m^-3 sr^-1 MeV^-1 × (MeV/MeV) = m^-3 sr^-1
    return 4.0 * math.pi * T * float(np.trapezoid(integrand, Es))


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_density_dat(filepath: str) -> list[dict]:
    """
    Parse gridless_points_density.dat (Tecplot ASCII, POINT data format).

    Returns list of dicts with keys: X_km, Y_km, Z_km, N_m3, N_cm3
    (order matches the order of points in the input file).
    """
    records = []
    in_zone = False
    var_names = []

    with open(filepath) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('!') or s.startswith('#'):
                continue
            if s.upper().startswith('VARIABLES'):
                # VARIABLES = "X_km" "Y_km" ...
                parts = s.split('=', 1)[1]
                var_names = [v.strip().strip('"') for v in parts.split()]
                continue
            if s.upper().startswith('ZONE'):
                in_zone = True
                continue
            if in_zone:
                try:
                    vals = [float(v) for v in s.split()]
                    if len(vals) == len(var_names):
                        records.append(dict(zip(var_names, vals)))
                except ValueError:
                    pass

    return records


def read_spectrum_dat(filepath: str) -> list[dict]:
    """
    Parse gridless_points_spectrum.dat (Tecplot ASCII, one ZONE per point).

    Returns a list of per-point dicts:
      { 'E_MeV': [...], 'T': [...], 'J_boundary': [...], 'J_local': [...] }
    Column names expected: E_MeV  T  J_boundary_perMeV  J_local_perMeV

    One entry per ZONE (= one entry per observation point).
    """
    points = []
    current_zone = None
    var_names = []

    with open(filepath) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('!') or s.startswith('#'):
                continue
            if s.upper().startswith('VARIABLES'):
                parts = s.split('=', 1)[1]
                var_names = [v.strip().strip('"') for v in parts.split()]
                continue
            if s.upper().startswith('ZONE'):
                if current_zone:
                    points.append(_finalise_zone(current_zone, var_names))
                current_zone = {v: [] for v in var_names}
                continue
            if current_zone is not None:
                try:
                    vals = [float(v) for v in s.split()]
                    if len(vals) == len(var_names):
                        for v, x in zip(var_names, vals):
                            current_zone[v].append(x)
                except ValueError:
                    pass

    if current_zone:
        points.append(_finalise_zone(current_zone, var_names))

    return points


def _finalise_zone(zone_dict: dict, var_names: list) -> dict:
    """Convert lists to numpy arrays inside a zone dict."""
    return {k: np.array(v, dtype=float) for k, v in zone_dict.items()}


# ---------------------------------------------------------------------------
# Per-point analytic comparison
# ---------------------------------------------------------------------------

def compare_point(tc: dict, density_rec: dict, spectrum_zone: dict,
                  tol: float = 0.05, verbose: bool = True) -> dict:
    """
    Compare a single test case against its analytic prediction.

    Parameters
    ----------
    tc           : test case descriptor dict (from TEST_CASES)
    density_rec  : row from gridless_points_density.dat
    spectrum_zone: zone from gridless_points_spectrum.dat
    tol          : relative tolerance for pass/fail
    verbose      : print a detailed report

    Returns
    -------
    dict with 'passed' (bool), 'ratio' (model/analytic), 'message' (str)
    """
    label     = tc["label"]
    r_Re      = tc["r_Re"]
    lam_deg   = tc["lam_deg"]
    E1        = tc["Emin_test"]
    E2        = tc["Emax_test"]
    test_type = tc["test_type"]

    # --- Analytic prediction ---
    Tg          = T_geometric(r_Re)
    n_analytic  = n_integral_analytic(E1, E2, Tg)
    n_free      = n_integral_analytic(E1, E2, 1.0)
    Rc_mx       = Rc_max_GV(r_Re, lam_deg)
    E_cutoff_MeV = E_kin_from_R_GV(Rc_mx) if Rc_mx > 0 else 0.0

    # --- Model: integrate spectrum output over the test energy sub-range ---
    E_arr = spectrum_zone.get("E_MeV", np.array([]))
    T_arr = spectrum_zone.get("T",     np.array([]))
    J_arr = spectrum_zone.get("J_local_perMeV", np.array([]))

    if E_arr.size == 0:
        return {"passed": False, "ratio": float("nan"),
                "message": f"{label}: spectrum zone is empty"}

    # Select energy bins within [E1, E2]
    mask = (E_arr >= E1) & (E_arr <= E2)
    if mask.sum() < 2:
        return {"passed": False, "ratio": float("nan"),
                "message": f"{label}: fewer than 2 energy bins in [{E1},{E2}] MeV"}

    E_sub = E_arr[mask]
    T_sub = T_arr[mask]

    # Compute model density in the sub-range:
    #   n_sub = 4pi * trapezoid( T(E)*J_b(E)/v(E) dE )
    Jb_sub = J0_SPEC * (E_sub / E0_SPEC) ** (-GAMMA)
    vs_sub = np.array([v_from_E_MeV(e) for e in E_sub])
    n_model = 4.0 * math.pi * float(np.trapezoid(T_sub * Jb_sub / vs_sub, E_sub))

    # Mean transmissivity in the sub-range
    T_mean_model = float(np.mean(T_sub))

    # --- Pass/fail ---
    if test_type == "exact":
        rel_err = abs(n_model - n_analytic) / n_analytic
        passed  = rel_err < tol
        ratio   = n_model / n_analytic if n_analytic > 0 else float("nan")
        verdict = "PASS" if passed else "FAIL"
        msg = (f"{label}\n"
               f"  Test type  : {test_type}  (|n_model-n_analytic|/n_analytic < {tol:.0%})\n"
               f"  T_geo      = {Tg:.5f}  (geometric shadow)\n"
               f"  Rc_max     = {Rc_mx:.3f} GV  → {E_cutoff_MeV:.0f} MeV\n"
               f"  E test range: {E1:.0f}–{E2:.0f} MeV  ({mask.sum()} bins used)\n"
               f"  n_free     = {n_free:.4e} m^-3  (T=1 reference)\n"
               f"  n_analytic = {n_analytic:.4e} m^-3  (T=T_geo={Tg:.5f})\n"
               f"  n_model    = {n_model:.4e} m^-3  (from spectrum output)\n"
               f"  T_mean(model) = {T_mean_model:.5f}  (should ≈ {Tg:.5f})\n"
               f"  rel_err    = {rel_err:.2%}  (tol={tol:.0%})\n"
               f"  >>> {verdict} <<<")

    else:  # upper_bound
        upper = 0.5 * n_analytic   # n_model must be < 50% of T_geo reference
        passed = (n_model < upper)
        ratio  = n_model / n_analytic if n_analytic > 0 else float("nan")
        verdict = "PASS" if passed else "FAIL"
        msg = (f"{label}\n"
               f"  Test type  : {test_type}  (n_model < 50% of T_geo reference)\n"
               f"  T_geo      = {Tg:.5f}  (geometric shadow only)\n"
               f"  Rc_vert    = {Rc_max_GV(r_Re,lam_deg)/4:.3f} GV → {E_kin_from_R_GV(Rc_max_GV(r_Re,lam_deg)/4):.0f} MeV  (vertical)\n"
               f"  Rc_max     = {Rc_mx:.3f} GV → {E_cutoff_MeV:.0f} MeV  (horizontal-east)\n"
               f"  E test range: {E1:.0f}–{E2:.0f} MeV  ({mask.sum()} bins used)\n"
               f"  n_free     = {n_free:.4e} m^-3  (T=1 reference)\n"
               f"  n_Tgeo     = {n_analytic:.4e} m^-3  (T=T_geo upper bound)\n"
               f"  50% upper  = {upper:.4e} m^-3  (pass threshold)\n"
               f"  n_model    = {n_model:.4e} m^-3  (from spectrum output)\n"
               f"  T_mean(model) = {T_mean_model:.5f}  (should << {Tg:.5f})\n"
               f"  ratio n_model/n_Tgeo = {ratio:.3f}  (should << 0.5)\n"
               f"  >>> {verdict} <<<")

    if verbose:
        print(msg)
        print()

    return {"passed": passed, "ratio": ratio, "message": msg,
            "n_model": n_model, "n_analytic": n_analytic, "T_geo": Tg,
            "T_mean_model": T_mean_model, "E_sub": E_sub, "T_sub": T_sub}


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary_header():
    print("=" * 70)
    print("AMPS DENSITY/SPECTRUM ANALYTIC TEST SUITE — DIPOLE FIELD")
    print("=" * 70)
    print()
    print(f"Proton rest mass     : {MP_MEV:.3f} MeV")
    print(f"Earth dipole B_eq(Re): {B_EQ_RE:.3e} T")
    print(f"Størmer const Cs     : {CS_GV:.3f} GV·Re²  (= mu0/4pi*M_E/Re^2 * c)")
    print(f"Inner sphere r_inner : {R_INNER_KM:.1f} km = {R_INNER_KM/RE_KM:.3f} Re")
    print()
    print(f"Spectrum: J_b(E) = {J0_SPEC:.1e} * (E/{E0_SPEC:.0f} MeV)^(-{GAMMA:.1f})")
    print(f"          [m^-2 s^-1 sr^-1 MeV^-1]")
    print()


def print_summary_table(results: list[dict]):
    print()
    print("-" * 70)
    print("SUMMARY")
    print("-" * 70)
    print(f"{'Case':<35} {'n_model':>12} {'n_analytic':>12} {'ratio':>8} {'Result':>8}")
    print("-" * 70)
    n_pass = 0
    for i, r in enumerate(results):
        tc = TEST_CASES[i]
        status = "PASS" if r["passed"] else "FAIL"
        n_pass += int(r["passed"])
        nm = f"{r.get('n_model', float('nan')):.3e}"
        na = f"{r.get('n_analytic', float('nan')):.3e}"
        ratio = f"{r.get('ratio', float('nan')):.3f}"
        print(f"  {tc['label'][:33]:<33} {nm:>12} {na:>12} {ratio:>8} {status:>8}")
    print("-" * 70)
    print(f"  Passed: {n_pass}/{len(results)}")
    print("=" * 70)
    if n_pass == len(results):
        print("ALL TESTS PASSED")
    else:
        print(f"WARNING: {len(results) - n_pass} TEST(S) FAILED")
    print("=" * 70)


# ---------------------------------------------------------------------------
# Optional matplotlib T(E) plots
# ---------------------------------------------------------------------------

def make_plots(results: list[dict], out_dir: str = "."):
    """
    Plot T(E) from AMPS output vs analytic expectation for each test case.
    Requires matplotlib.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[plots] matplotlib not available — skipping figures")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.ravel()

    colors = {"TC1": "#1f77b4", "TC2": "#2ca02c", "TC3": "#ff7f0e", "TC4": "#d62728"}

    for i, (r, ax) in enumerate(zip(results, axes)):
        tc = TEST_CASES[i]
        key = f"TC{i+1}"
        color = colors.get(key, "black")

        E_sub = r.get("E_sub", np.array([]))
        T_sub = r.get("T_sub", np.array([]))
        Tg    = r.get("T_geo", float("nan"))

        if E_sub.size > 0:
            ax.semilogx(E_sub, T_sub, "o-", color=color, ms=4, lw=1.5,
                        label="AMPS model T(E)")

        # Analytic expectation
        E_dense = np.logspace(np.log10(tc["Emin_test"]),
                              np.log10(tc["Emax_test"]), 200)
        if tc["test_type"] == "exact":
            ax.axhline(Tg, color="k", ls="--", lw=1.5,
                       label=f"T_geo = {Tg:.4f} (analytic)")
        else:
            ax.axhline(Tg, color="k", ls="--", lw=1.5,
                       label=f"T_geo = {Tg:.4f} (upper bound)")
            ax.axhline(0.5 * Tg, color="gray", ls=":", lw=1.5,
                       label=f"50% × T_geo (pass threshold)")

        # Mark the Størmer cutoff
        Rc_max_here = Rc_max_GV(tc["r_Re"], tc["lam_deg"])
        if Rc_max_here > 0:
            E_cm = E_kin_from_R_GV(Rc_max_here)
            if tc["Emin_test"] < E_cm < tc["Emax_test"]:
                ax.axvline(E_cm, color="purple", ls="-.", lw=1.0, alpha=0.6,
                           label=f"Rc_max = {Rc_max_here:.2f} GV ({E_cm:.0f} MeV)")

        status = "PASS ✓" if r["passed"] else "FAIL ✗"
        ax.set_title(f"{key}: {status}\n{tc['r_Re']:.0f} Re, lam={tc['lam_deg']:.0f}°",
                     fontsize=9)
        ax.set_xlabel("E [MeV]")
        ax.set_ylabel("Transmissivity T")
        ax.set_ylim(-0.05, 1.15)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    fig.suptitle("AMPS gridless: T(E) — model vs analytic (DIPOLE field)", fontsize=11)
    plt.tight_layout()
    out_path = os.path.join(out_dir, "test_density_T_comparison.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[plots] Saved: {out_path}")
    plt.close(fig)

    # Second figure: density bar chart
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    labels  = [f"TC{i+1}" for i in range(len(results))]
    n_model = [r.get("n_model", 0) for r in results]
    n_analy = [r.get("n_analytic", 0) for r in results]

    x = np.arange(len(labels))
    w = 0.35
    bars_m = ax2.bar(x - w/2, n_model, w, label="n_model (AMPS)", color="#1f77b4")
    bars_a = ax2.bar(x + w/2, n_analy, w, label="n_analytic (T_geo)", color="#ff7f0e",
                     alpha=0.7)

    ax2.set_yscale("log")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("Number density  n  [m⁻³]")
    ax2.set_title("AMPS model vs analytic density (DIPOLE field)")
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis="y")

    out_path2 = os.path.join(out_dir, "test_density_n_comparison.png")
    plt.savefig(out_path2, dpi=150, bbox_inches="tight")
    print(f"[plots] Saved: {out_path2}")
    plt.close(fig2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spectrum-dir", default=".",
                        help="Directory containing AMPS output files (default: '.')")
    parser.add_argument("--tol", type=float, default=0.05,
                        help="Relative tolerance for TC1–TC3 pass/fail (default: 0.05)")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip matplotlib figure generation")
    args = parser.parse_args()

    density_file  = os.path.join(args.spectrum_dir, "gridless_points_density.dat")
    spectrum_file = os.path.join(args.spectrum_dir, "gridless_points_spectrum.dat")

    print_summary_header()

    # --- Load AMPS output ---
    if not os.path.isfile(density_file):
        print(f"ERROR: density file not found: {density_file}")
        print("       Run AMPS with test_density_dipole_4pts.in first.")
        sys.exit(1)
    if not os.path.isfile(spectrum_file):
        print(f"ERROR: spectrum file not found: {spectrum_file}")
        sys.exit(1)

    density_rows = read_density_dat(density_file)
    spectrum_zones = read_spectrum_dat(spectrum_file)

    if len(density_rows) < len(TEST_CASES):
        print(f"WARNING: expected {len(TEST_CASES)} points in density file, "
              f"found {len(density_rows)}. Will test what is available.")

    if len(spectrum_zones) < len(TEST_CASES):
        print(f"WARNING: expected {len(TEST_CASES)} zones in spectrum file, "
              f"found {len(spectrum_zones)}. Will test what is available.")

    # --- Print analytic cutoff table (informational) ---
    print("Størmer cutoff table (DIPOLE, protons):")
    print(f"  {'r [Re]':>8}  {'lam':>5}  {'Rc_vert [GV]':>14}  {'E_vert [MeV]':>14}  "
          f"{'Rc_max [GV]':>12}  {'E_max [MeV]':>12}  {'T_geo':>8}")
    for tc in TEST_CASES:
        Rv = Rc_max_GV(tc["r_Re"], tc["lam_deg"]) / 4
        Rm = Rc_max_GV(tc["r_Re"], tc["lam_deg"])
        Ev = E_kin_from_R_GV(Rv) if Rv > 1e-6 else 0.0
        Em = E_kin_from_R_GV(Rm) if Rm > 1e-6 else 0.0
        Tg = T_geometric(tc["r_Re"])
        print(f"  {tc['r_Re']:>8.1f}  {tc['lam_deg']:>5.0f}  "
              f"{Rv:>14.3f}  {Ev:>14.0f}  "
              f"{Rm:>12.3f}  {Em:>12.0f}  {Tg:>8.5f}")
    print()

    # --- Run comparisons ---
    print("DETAILED TEST RESULTS")
    print("-" * 70)
    results = []
    for i, tc in enumerate(TEST_CASES):
        if i >= len(spectrum_zones):
            results.append({"passed": False, "ratio": float("nan"),
                             "message": f"Point {i} missing from output",
                             "n_model": float("nan"), "n_analytic": float("nan")})
            continue
        dr = density_rows[i] if i < len(density_rows) else {}
        sz = spectrum_zones[i]
        result = compare_point(tc, dr, sz, tol=args.tol, verbose=True)
        results.append(result)

    # --- Summary table ---
    print_summary_table(results)

    # --- Optional plots ---
    if not args.no_plots:
        make_plots(results, out_dir=args.spectrum_dir)

    # Exit code: 0 if all passed, 1 if any failed
    n_failed = sum(1 for r in results if not r["passed"])
    sys.exit(n_failed)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Collection-radius convergence report of the srcSEP3D 2013-04-11 test.

Usage:
  convergence_report.py [options] RADIUS_AU=RUN_DIR RADIUS_AU=RUN_DIR ...

Each RUN_DIR is a finished run of event3d.convergence-r<radius>au.<mover>.in
(run_convergence.sh). The script reads every publication's spacecraft.csv
under RUN_DIR/PT/plots/sep3d. For each observer (default: tube-earth,
tube-soho, tube-ace) it builds the time profile of the differential
intensity in broad energy bands:

  j = Y v / (4 pi V),  V = (4/3) pi rho^3,

where Y = represented_particles_J-1 (sum of weights in the collection sphere
per joule) and v is the proton speed at the geometric centre of each energy
bin. The uncertainty is the same transform of standard_uncertainty_J-1.

Consecutive radii are compared pairwise (rho_small < rho_large). The two
runs build different tubes (the tube width is rho + s), so they sample the
flux with different macroparticle densities, and the importance weights vary
by more than an order of magnitude inside one band. Cells are therefore
never chosen by their own value, which would select upward fluctuations:
  - an (observer, band) enters when both runs hold at least --min-count
    effective macroparticles in it over the whole run,
    N_eff = (sum j)^2 / sum sigma^2 over its publications;
  - in an entered band every publication where either run is nonzero is a
    cell (a cell where both are zero carries no information), with

      z = (j_large - j_small) / sqrt(sigma_large^2 + sigma_small^2).

An empty sphere has sigma = 0, so a cell where only one run holds particles
is scored by that run's uncertainty alone. That overstates |z| when the
other run should have held a few particles, which makes the test stricter,
never more lenient.

A pair agrees when chi^2 = sum z^2 over nu cells does not exceed the
upper --sigma quantile of the chi^2_nu distribution (Wilson-Hilferty).
Successive publications of one run share particles, so the cells are
positively correlated in time. chi^2 then keeps its mean nu but has a larger
spread, which makes the test conservative (it rejects too often), never
lenient. The selected radius is the largest one that agrees with the next
smaller radius (README.md, section 6).

Exit status: 0 when a radius is selected, 1 when none is, 2 on bad input.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

C = 299792458.0                 # speed of light [m/s]
PROTON_KG = 1.67262192369e-27   # proton mass [kg] (CODATA 2018)
EV_J = 1.602176634e-19          # 1 eV [J]
AU_M = 1.495978707e11           # astronomical unit [m]
# Kinetic-energy bands [MeV]: one decade each over the 20 keV - 200 MeV band
# of the source and the tube observers (README.md section 6).
DEFAULT_BANDS_MEV = ((0.02, 0.2), (0.2, 2.0), (2.0, 20.0), (20.0, 200.0))


def speed(kinetic_j):
    """Proton speed [m/s] at kinetic energy kinetic_j [J] (relativistic)."""
    gamma = 1.0 + kinetic_j / (PROTON_KG * C * C)
    return C * math.sqrt(1.0 - 1.0 / (gamma * gamma))


def publications(run_dir):
    """[(time_s, spacecraft.csv rows)] of one run, in time order."""
    root = Path(run_dir) / 'PT' / 'plots' / 'sep3d'
    found = []
    # Publication directories are <prefix>.<zero-padded sequence>; the digit
    # test skips files such as <prefix>.radial-energy.csv.  The list is put in
    # time order by the manifest's simulation time below, not by the name.
    for directory in sorted(root.glob('*.[0-9]*')):
        if not directory.is_dir() or not directory.name.rsplit('.', 1)[-1].isdigit():
            continue
        manifest = directory / 'manifest.txt'
        table = directory / 'spacecraft.csv'
        if not manifest.is_file() or not table.is_file():
            continue
        time_s = None
        for line in manifest.read_text().splitlines():
            if line.startswith('simulation_time_s='):
                time_s = float(line.split('=', 1)[1])
        if time_s is None:
            raise ValueError(f'{manifest} has no simulation_time_s')
        with table.open(newline='') as stream:
            found.append((time_s, list(csv.DictReader(stream))))
    if not found:
        raise ValueError(f'no publication under {root}')
    return found


def band_profiles(run_dir, rho_m, observers, bands_j):
    """{(observer, band index, time_s): (j, sigma_j, n_eff)} of one run.

    rho_m is the collection radius of the run [m]; observers is the set of
    observer names to read; bands_j the band edges [J].  For each publication
    and (observer, band):
      j        band-averaged differential intensity [1/(m^2 s sr J)],
               sum over the band's bins of Y dE v/(4 pi V), divided by the
               band width;
      sigma_j  its standard uncertainty (bins added in quadrature);
      n_eff    (sum w)^2 / sum w^2 of the band's represented particles, the
               effective number of macroparticles (not used by compare(),
               which forms its own whole-run count from j and sigma_j).
    Only (observer, band, time) cells with at least one row are present.
    """
    volume = 4.0 / 3.0 * math.pi * rho_m ** 3
    out = {}
    for time_s, rows in publications(run_dir):
        acc = {}
        for row in rows:
            name = row['spacecraft']
            if name not in observers:
                continue
            lo, hi = float(row['energy_min_J']), float(row['energy_max_J'])
            # A bin belongs to the band that holds its geometric centre (the
            # observer bins are logarithmic).  A bin across a band edge is
            # counted whole in one band; both runs use the same bins, so the
            # comparison is unaffected.
            centre = math.sqrt(lo * hi)
            band = next((b for b, (blo, bhi) in enumerate(bands_j)
                         if centre >= blo and centre < bhi), None)
            if band is None:
                continue
            width = hi - lo
            # Y = represented_particles_J-1 is sum(w)/dE inside the sphere, so
            # Y dE is the represented particle number of the bin; the
            # uncertainty sqrt(sum w^2)/dE transforms the same way.
            weight = float(row['represented_particles_J-1']) * width
            sigma = float(row['standard_uncertainty_J-1']) * width
            # Number in the sphere -> intensity: density N/V times v/(4 pi),
            # with v at the bin centre, averaged over the band width.
            scale = speed(centre) / (4.0 * math.pi * volume * (bands_j[band][1] - bands_j[band][0]))
            a = acc.setdefault((name, band), [0.0, 0.0, 0.0, 0.0])
            a[0] += weight * scale           # band-averaged intensity
            a[1] += (sigma * scale) ** 2     # its variance
            a[2] += weight                   # sum w
            a[3] += sigma * sigma            # sum w^2
        for (name, band), (j, var, sw, sw2) in acc.items():
            n_eff = sw * sw / sw2 if sw2 > 0.0 else 0.0
            out[(name, band, time_s)] = (j, math.sqrt(var), n_eff)
    return out


def chi2_upper_quantile(nu, sigma):
    """Upper quantile of chi^2_nu at a one-sided Gaussian z = sigma (Wilson-Hilferty)."""
    a = 2.0 / (9.0 * nu)
    return nu * (1.0 - a + sigma * math.sqrt(a)) ** 3


def compare(small, large, min_count):
    """chi^2 comparison of two band_profiles() results (rho_small < rho_large).

    Returns (chi2, nu, worst, bands):
      chi2, nu  sum of z^2 over the nu scored cells;
      worst     (z, (observer, band, time_s)) of the cell with the largest
                |z|, or None when nu = 0;
      bands     {(observer, band): [sum j_small, sum j_large, cells]} over the
                scored cells, for the per-band ratio printed by main().
    The entry rule and the cell rule are those of the module docstring.
    """
    # Whole-run totals per (observer, band) and run: [sum j, sum sigma^2] of
    # the small run, then of the large run.
    totals = {}
    for which, profile in ((0, small), (1, large)):
        for (name, band, _), (j, sigma, _) in profile.items():
            t = totals.setdefault((name, band), [0.0, 0.0, 0.0, 0.0])
            t[2 * which] += j
            t[2 * which + 1] += sigma * sigma
    # An (observer, band) enters when both runs hold at least min_count
    # effective macroparticles over the run, N_eff = (sum j)^2 / sum sigma^2.
    # The selection uses whole-run totals of both runs, never the value of a
    # single cell.
    entered = {key for key, (js, vs, jl, vl) in totals.items()
               if vs > 0.0 and vl > 0.0 and js * js / vs >= min_count and jl * jl / vl >= min_count}
    chi2, nu, worst = 0.0, 0, None
    bands = {}
    for key in sorted(set(small) | set(large)):
        if (key[0], key[1]) not in entered:
            continue
        # A cell missing from one run (no row of that observer and band at that
        # time) counts as an empty sphere of that run: j = sigma = 0.
        j_s, s_s, _ = small.get(key, (0.0, 0.0, 0.0))
        j_l, s_l, _ = large.get(key, (0.0, 0.0, 0.0))
        if j_s == 0.0 and j_l == 0.0:
            continue
        # standard_uncertainty_J-1 is sqrt(sum w^2)/dE (output/sampling.cpp), so
        # a run with a nonzero j here has sigma > 0: the denominator is positive.
        z = (j_l - j_s) / math.hypot(s_l, s_s)
        chi2 += z * z
        nu += 1
        b = bands.setdefault((key[0], key[1]), [0.0, 0.0, 0])
        b[0] += j_s
        b[1] += j_l
        b[2] += 1
        if worst is None or abs(z) > abs(worst[0]):
            worst = (z, key)
    return chi2, nu, worst, bands


def main():
    """Compare consecutive radii pairwise, print the result and select a radius."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('runs', nargs='+', metavar='RADIUS_AU=RUN_DIR')
    parser.add_argument('--observers', default='tube-earth,tube-soho,tube-ace')
    parser.add_argument('--min-count', type=float, default=10.0,
                        help='minimum effective macroparticles of an (observer, band) over a run, in both runs (default 10)')
    parser.add_argument('--sigma', type=float, default=3.0,
                        help='one-sided Gaussian level of the chi^2 acceptance (default 3)')
    args = parser.parse_args()
    try:
        runs = sorted((float(item.split('=', 1)[0]), item.split('=', 1)[1]) for item in args.runs)
    except (ValueError, IndexError):
        parser.error('each run is RADIUS_AU=RUN_DIR')
    if len(runs) < 2:
        parser.error('at least two radii are needed')
    observers = set(args.observers.split(','))
    # Band edges in joules, for comparison with the energy_min_J/energy_max_J
    # columns of spacecraft.csv.
    bands_j = [(lo * 1e6 * EV_J, hi * 1e6 * EV_J) for lo, hi in DEFAULT_BANDS_MEV]
    try:
        profiles = [(r, band_profiles(d, r * AU_M, observers, bands_j)) for r, d in runs]
    except (OSError, ValueError, KeyError) as error:
        print(f'convergence_report: {error}', file=sys.stderr)
        return 2
    print(f'observers {sorted(observers)}; bands (MeV) {list(DEFAULT_BANDS_MEV)}; '
          f'min band count {args.min_count:g}; acceptance chi2 <= {args.sigma:g}-sigma quantile')
    # Radii are sorted ascending, so the loop visits each consecutive pair
    # from the smallest up and keeps the larger radius of the last pair that
    # agrees: the largest radius that agrees with the next smaller one.  A
    # pair without statistics is reported and skipped.
    selected = None
    for (r_s, small), (r_l, large) in zip(profiles, profiles[1:]):
        chi2, nu, worst, bands = compare(small, large, args.min_count)
        if nu == 0:
            print(f'{r_s:g} au vs {r_l:g} au: no cell has enough statistics; not decided')
            continue
        limit = chi2_upper_quantile(nu, args.sigma)
        agree = chi2 <= limit
        print(f'{r_s:g} au vs {r_l:g} au: chi2/nu = {chi2 / nu:.3f} over nu = {nu} cells '
              f'(limit {limit / nu:.3f}) -> {"agree" if agree else "DIFFER"}; '
              f'largest |z| = {abs(worst[0]):.2f} at {worst[1][0]} band {worst[1][1]} t = {worst[1][2]:.0f} s')
        for (name, band), (sum_s, sum_l, count) in sorted(bands.items()):
            print(f'    {name:11s} {DEFAULT_BANDS_MEV[band][0]:g}-{DEFAULT_BANDS_MEV[band][1]:g} MeV: '
                  f'sum j({r_l:g} au)/sum j({r_s:g} au) = '
                  f'{sum_l / sum_s if sum_s > 0 else float("inf"):.3f} over {count} publications')
        if agree:
            selected = r_l
    if selected is None:
        print('no radius agrees with the next smaller one: extend the study or the statistics')
        return 1
    print(f'selected collection radius: {selected:g} au')
    return 0


if __name__ == '__main__':
    sys.exit(main())

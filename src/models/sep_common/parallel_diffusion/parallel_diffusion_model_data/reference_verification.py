#!/usr/bin/env python3
"""Standalone verification of the parallel-diffusion companion bundle (bundle revision 1.4).

Usage, from the extracted parallel_diffusion_model_data/ directory:

    python3 reference_verification.py              # digests, schema, deterministic fixtures
    python3 reference_verification.py --audit      # also re-solve the 300-state NLGCE-F/NLGCE-N audit
    python3 reference_verification.py --broadened  # also recompute the broadened-slab fixtures (slow)

Dependencies: Python 3.9+, NumPy, SciPy. Every reference value is recomputed from the
equations of PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md and compared with the stored value at
the tolerance declared in benchmark_points.json. Exit status is nonzero if any check fails.
"""
import csv
import hashlib
import json
import math
import os
import sys
from decimal import Decimal, getcontext

import numpy as np
from scipy import integrate, optimize
from scipy.special import gamma as G

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC_DIGESTS = {
    'NLGCE_F_2014_parallel.csv': '7cdc5ab9cddda0c7295a25bfcaff4ba3422002da12ac1b4660a98dcc64eaa762',
    'NLGCE_F_2014_perpendicular.csv': '7371b6557f40c5efa7a48460a5f1b2971374a2c43d0eb3e209f0a8e62670356a',
}
NU = 5.0 / 6.0
C_NU = G(NU) / (2.0 * math.sqrt(math.pi) * G(NU - 0.5))
BOX = {'r': (1e-5, 6.3), 'fs': (1e-3, 0.85), 'eps2': (1e-4, 1e2), 'rho': (1.0, 1e3)}

RESULTS = []


def record(name, ok, detail=''):
    RESULTS.append((name, ok))
    print(('PASS ' if ok else 'FAIL ') + name + (('  ' + detail) if detail else ''))
    sys.stdout.flush()


def path(f):
    return os.path.join(HERE, f)


def rel(a, b):
    return abs(a / b - 1.0) if b != 0 else abs(a)


# ------------------------------------------------------------------------- data integrity
def check_digests():
    with open(path('SHA256SUMS')) as fh:
        lines = [ln.split() for ln in fh if ln.strip()]
    for digest, name in lines:
        got = hashlib.sha256(open(path(name), 'rb').read()).hexdigest()
        record('sha256 ' + name, got == digest)
    for name, digest in SPEC_DIGESTS.items():
        got = hashlib.sha256(open(path(name), 'rb').read()).hexdigest()
        record('specification digest ' + name, got == digest)


def load_coeffs(name):
    with open(path(name), newline='') as fh:
        rd = csv.DictReader(fh)
        ok = rd.fieldnames == ['j', 'k', 'l', 'd_i0', 'd_i1', 'd_i2', 'd_i3', 'd_i4', 'd_i5']
        d = {}
        dup = False
        for row in rd:
            j, k, l = int(row['j']), int(row['k']), int(row['l'])
            for i in range(6):
                key = (i, j, k, l)
                dup |= key in d
                d[key] = row['d_i%d' % i]
    full = sorted(d) == [(i, j, k, l) for i in range(6) for j in range(4) for k in range(4) for l in range(3)]
    record('schema and index coverage ' + name, ok and full and not dup, '%d coefficients' % len(d))
    return d


# ------------------------------------------------------------------------- polynomial
def poly_F(d, x, use_decimal=False):
    if use_decimal:
        getcontext().prec = 50
        x = [t if isinstance(t, Decimal) else Decimal(repr(float(t))) for t in x]
        conv = Decimal
    else:
        conv = float
    a_coef = []
    for i in range(6):
        b_coef = []
        for j in range(4):
            c_coef = []
            for k in range(4):
                acc = conv(0)
                for l in (2, 1, 0):
                    acc = acc * x[3] + conv(d[(i, j, k, l)])
                c_coef.append(acc)
            acc = conv(0)
            for k in (3, 2, 1, 0):
                acc = acc * x[2] + c_coef[k]
            b_coef.append(acc)
        acc = conv(0)
        for j in (3, 2, 1, 0):
            acc = acc * x[1] + b_coef[j]
        a_coef.append(acc)
    acc = conv(0)
    for i in (5, 4, 3, 2, 1, 0):
        acc = acc * x[0] + a_coef[i]
    return acc


def poly_grad_decimal(d, x, h='1e-12'):
    getcontext().prec = 50
    xd = [Decimal(repr(float(t))) for t in x]
    hd = Decimal(h)
    g = []
    for a in range(4):
        xp, xm = list(xd), list(xd)
        xp[a] += hd
        xm[a] -= hd
        g.append(float((poly_F(d, xp, True) - poly_F(d, xm, True)) / (2 * hd)))
    return g


# ------------------------------------------------------------------------- closures
def a_x(r, fs, eps2):
    eps = math.sqrt(eps2)
    xi = r / (C_NU * eps)
    return 0.5 * math.sqrt(fs / ((xi / (1.0 + xi)) / eps + eps / (2.0 * xi)))


def a_prime2(fs, rho):
    return 1.0 / (math.sqrt(1.0 / rho) / fs + (4.0 / 3.0) / (1.0 - fs))


def _quad_u(f, pts, U=40.0, epsrel=1e-13):
    pts = sorted(p for p in set(pts) if -U < p < U)
    edges = [-U] + pts + [U]
    return sum(integrate.quad(f, a, b, epsabs=0.0, epsrel=epsrel, limit=800)[0] for a, b in zip(edges[:-1], edges[1:]))


def _bp(*pairs):
    out = [0.0]
    for t, b in pairs:
        if b > 0 and t > 0 and math.isfinite(b) and math.isfinite(t) and math.isfinite(t / b) and t / b > 0:
            out.append(0.5 * math.log(t / b))
    return out


def resonant(alpha, beta, Om):
    def f(u):
        x = math.exp(u)
        A = alpha + beta * x * x
        return 4.0 * C_NU * (1.0 + x * x) ** (-NU) * A / (Om * Om + A * A) * x
    return _quad_u(f, _bp((alpha, beta), (Om, beta)))


def inverse(alpha, beta):
    def f(u):
        x = math.exp(u)
        return 4.0 * C_NU * (1.0 + x * x) ** (-NU) / (alpha + beta * x * x) * x
    return _quad_u(f, _bp((alpha, beta)))


def residuals(y, st, model, kx_given=None):
    if any((not math.isfinite(t)) or abs(t) > 300.0 for t in y):
        return [1e3, 1e3]
    r, fs, eps2, rho = st
    kz = math.exp(y[0])
    kx = kx_given if kx_given is not None else math.exp(y[1])
    Om, l2 = 1.0 / r, 1.0 / rho
    dBs2, dB22 = fs * eps2, (1.0 - fs) * eps2
    alpha = 1.0 / (3.0 * kz)
    R1 = y[0] + math.log(3.0 * a_x(r, fs, eps2) * Om * Om *
                         (dBs2 * resonant(alpha, kz, Om) + dB22 * resonant(alpha, kx / l2 ** 2, Om)))
    if model == 'nlpa_given_perp':
        return [R1, 0.0]
    J2 = inverse(alpha, kx / l2 ** 2)
    if model == 'nlgce_n':
        rhs = a_prime2(fs, rho) / 6.0 * dB22 * J2
    else:  # nlgc_e, a^2 = 1/3
        rhs = (1.0 / 3.0) / 6.0 * (dBs2 * inverse(alpha, kz) + dB22 * J2)
    return [R1, y[1] - math.log(rhs)]


def solve(st, model, y0, kx_given=None):
    if model == 'nlpa_given_perp':
        sol = optimize.root(lambda v: [residuals([v[0], 0.0], st, model, kx_given)[0]], [y0[0]], method='hybr',
                            options={'xtol': 1e-15})
        y = [sol.x[0], math.log(kx_given)]
    else:
        sol = optimize.root(lambda v: residuals(v, st, model), y0, method='hybr', options={'xtol': 1e-15})
        y = list(sol.x)
    res = max(abs(t) for t in residuals(y, st, model, kx_given))
    return 3.0 * math.exp(y[0]), 3.0 * math.exp(y[1]), res, y


# ------------------------------------------------------------------------- deterministic fixtures
def check_particles(B):
    getcontext().prec = 50
    c = Decimal(299792458)
    e = Decimal('1.602176634e-19')
    AU = Decimal(149597870700)
    tol = B['tolerances']['particles']['rel']
    for name, P in B['particles'].items():
        m = Decimal(P['mass_kg'])
        T = Decimal(P['kinetic_energy_total_eV']) * e
        if 'kinetic_energy_per_nucleon_eV' in P:
            ok_n = Decimal(P['kinetic_energy_per_nucleon_eV']) * P['nucleon_count'] == Decimal(P['kinetic_energy_total_eV'])
            record('energy-per-nucleon conversion ' + name, ok_n)
        mc2 = m * c * c
        pc = (T * (T + 2 * mc2)).sqrt()
        v = pc * c / (T + mc2)
        lam = AU / 10
        calc = {'gamma': 1 + T / mc2, 'beta': v / c, 'speed_m_s': v, 'momentum_kg_m_s': pc / c,
                'rigidity_V': pc / (abs(P['charge_number']) * e), 'lambda_parallel_m': lam,
                'kappa_parallel_m2_s': v * lam / 3}
        worst = max(rel(float(calc[k]), P[k]) for k in calc)
        record('particle kinematics ' + name, worst <= tol, 'max rel %.1e' % worst)


def check_spectrum(B):
    S = B['spectrum_and_pitch_angle']
    tol = B['tolerances']['spectrum_and_pitch_angle']['rel']
    s = 2 * NU
    norm = integrate.quad(lambda x: 4 * C_NU * (1 + x * x) ** (-NU), 0, np.inf, epsabs=0, epsrel=1e-13)[0]
    record('Eq. (27) spectrum normalisation', abs(norm - 1) < 1e-12, '%.3e' % (norm - 1))
    calc = {'C_5_6': C_NU, 'Lc_over_ell_s_5_6': 2 * math.pi * C_NU,
            'inertial_prefactor_s_5_3': 3 / (2 * math.pi * C_NU * (2 - s) * (4 - s))}
    h = 0.01
    I = 2 * integrate.quad(lambda z: 3 * (1 - z ** 6) * z * z / (z * z + h), 0, 1, epsabs=0, epsrel=1e-13)[0]
    calc['I_q5_3_h0_01'] = I
    calc['D0_lambda_over_v_q5_3_h0_01'] = 3 * I / 8
    calc['I_q5_3_h0_exact_Eq24'] = 4 / ((2 - 5 / 3) * (4 - 5 / 3))
    I0 = 2 * integrate.quad(lambda z: 3 * (1 - z ** 6), 0, 1, epsabs=0, epsrel=1e-13)[0]
    record('Eq. (24) closed form against quadrature', rel(I0, calc['I_q5_3_h0_exact_Eq24']) < 1e-12)
    iso = 0.375 * integrate.quad(lambda m: (1 - m * m), -1, 1)[0]
    record('Eq. (21) isotropic normalisation', abs(iso - S['isotropic_check_Eq21']['lambda_times_D0_over_v']) < 1e-14)
    for k, v in calc.items():
        record('spectrum/pitch-angle ' + k, rel(v, S[k]) <= tol, 'rel %.1e' % rel(v, S[k]))
    tq = B['tolerances']['qlt_eq30']['rel']
    for key, Q in B['qlt_eq30'].items():
        r = Q['r_L_over_ell_s']
        J = integrate.quad(lambda z: 3 * (1 - z ** 6) * (1 + r * r * z ** 6) ** (s / 2), 0, 1, epsabs=0, epsrel=1e-13)[0]
        lam = 3 / (4 * math.pi * C_NU * Q['eps_s2']) * r ** (2 - s) * J
        record('QLT Eq. (30) r*=%s' % key, rel(lam, Q['lambda_par_over_ell_s_Eq30']) <= tq,
               'rel %.1e' % rel(lam, Q['lambda_par_over_ell_s_Eq30']))


def check_nlgce(B, dpar, dperp):
    tv = B['tolerances']['nlgce_f_2014_values']['rel']
    td = B['tolerances']['nlgce_f_2014_derivatives']['abs']
    tp = B['tolerances']['closure_parameters_a_x_a_prime2']['rel']
    tn = B['tolerances']['nlgce_n']['rel']
    tg = B['tolerances']['nlgc_e']['rel']
    tfd = B['tolerances']['nlgce_n_log_derivatives_finite_difference']['abs']
    for p, P in B['nlgce_f_2014'].items():
        i = P['inputs']
        st = (i['r_L_over_ell_s'], i['f_s'], i['eps2_total'], i['ell_s_over_ell_2'])
        inbox = all(BOX[k][0] <= v <= BOX[k][1] for k, v in zip(('r', 'fs', 'eps2', 'rho'), st))
        x = [math.log(t) for t in st]
        Fp, Fq = float(poly_F(dpar, x, True)), float(poly_F(dperp, x, True))
        w = max(rel(math.exp(Fp), P['lambda_par_over_ell_s']), rel(math.exp(Fq), P['lambda_perp_over_ell_s']))
        record('NLGCE-F values point %s' % p, inbox and w <= tv, 'max rel %.1e' % w)
        g = poly_grad_decimal(dpar, x) + poly_grad_decimal(dperp, x)
        wd = max(abs(a - b) for a, b in zip(g, P['dF_par_dx'] + P['dF_perp_dx']))
        record('NLGCE-F derivatives point %s' % p, wd <= td, 'max abs %.1e' % wd)
        N = B['nlgce_n'][p]
        wp = max(rel(a_x(*st[:3]), N['a_x']), rel(a_prime2(st[1], st[3]), N['a_prime2']))
        record('a_x and a_prime2 point %s' % p, wp <= tp, 'max rel %.1e' % wp)
        y0 = [math.log(math.exp(Fp) / 3) + 1.0, math.log(math.exp(Fq) / 3) - 1.0]  # deliberately offset start
        lp, lq, res, y = solve(st, 'nlgce_n', y0)
        w = max(rel(lp, N['lambda_par_over_ell_s']), rel(lq, N['lambda_perp_over_ell_s']))
        record('NLGCE-N point %s' % p, res < 1e-8 and w <= tn, 'max rel %.1e, residual %.1e' % (w, res))
        h = 1e-3
        names = (0, 1, 2, 3)
        worst = 0.0
        for a in names:
            vals = []
            for m in (-2, -1, 1, 2):
                st2 = list(st)
                st2[a] = st[a] * math.exp(m * h)
                l1, l2, rs, _ = solve(tuple(st2), 'nlgce_n', y)
                vals.append((math.log(l1), math.log(l2)))
            for kk, key in ((0, 'dlnlambda_par_dx'), (1, 'dlnlambda_perp_dx')):
                fm2, fm1, fp1, fp2 = (v[kk] for v in vals)
                worst = max(worst, abs((fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h) - N[key][a]))
        record('NLGCE-N log-derivatives point %s' % p, worst <= tfd, 'max abs %.1e' % worst)
        E = B['nlgc_e'][p]
        lp, lq, res, _ = solve(st, 'nlgc_e', y)
        w = max(rel(lp, E['lambda_par_over_ell_s']), rel(lq, E['lambda_perp_over_ell_s']))
        record('NLGC-E point %s' % p, res < 1e-8 and w <= tg, 'max rel %.1e, residual %.1e' % (w, res))
    ta = B['tolerances']['nlpa_given_perp']['rel']
    for key, P in B['nlpa_given_perp'].items():
        i = B['nlgce_f_2014'][P['state']]['inputs']
        st = (i['r_L_over_ell_s'], i['f_s'], i['eps2_total'], i['ell_s_over_ell_2'])
        lp, _, res, _ = solve(st, 'nlpa_given_perp', [math.log(P['lambda_par_over_ell_s'] / 3) - 2.0, 0.0],
                              kx_given=P['kappa_perp_supplied_over_v_ell_s'])
        w = rel(lp, P['lambda_par_over_ell_s'])
        record('NLPA given perp %s' % key, res < 1e-8 and w <= ta, 'rel %.1e, residual %.1e' % (w, res))
    # a deliberately wrong logarithm base must be detected
    i = B['nlgce_f_2014']['A']['inputs']
    x10 = [math.log10(t) for t in (i['r_L_over_ell_s'], i['f_s'], i['eps2_total'], i['ell_s_over_ell_2'])]
    record('log10 substitution is detected', rel(math.exp(poly_F(dpar, x10)), B['nlgce_f_2014']['A']['lambda_par_over_ell_s']) > 0.1)
    # the superseded exponentiated a_x must disagree at eps != 1 (point F) and agree at eps = 1 (point A)
    def ax_exp(r, fs, e2):
        eps = math.sqrt(e2)
        xi = r / (C_NU * eps)
        return 0.5 * math.sqrt(fs / ((xi / (1 + xi)) ** (1 / eps) + eps / (2 * xi)))
    iF = B['nlgce_f_2014']['F']['inputs']
    sF = (iF['r_L_over_ell_s'], iF['f_s'], iF['eps2_total'])
    record('exponentiated a_x is detected at point F', rel(ax_exp(*sF), B['nlgce_n']['F']['a_x']) > 1e-3)


# ------------------------------------------------------------------------- broadened slab (optional)
def check_broadened(B):
    tol = B['tolerances']['broadened_slab']['rel']

    def D(mu, rstar, eps2, kern, width):
        Om = 1.0 / rstar
        if kern == 'lorentzian':
            R = lambda w: width / (w * w + width * width)
        else:
            R = lambda w: math.sqrt(math.pi) / width * math.exp(-(w / width) ** 2)

        def f(u):
            k = math.exp(u)
            return 4.0 * C_NU * eps2 * (1.0 + k * k) ** (-NU) * (R(k * mu - Om) + R(k * mu + Om)) * k
        pts = [0.0]
        for n in [0.0] + [sg * 2.0 ** m for m in range(60) for sg in (-1.0, 1.0)]:
            kk = (Om + n * width) / mu
            if kk > 0:
                pts.append(math.log(kk))
        lo, hi = -50.0, 150.0
        pts = sorted(p for p in set(pts) if lo < p < hi)
        edges = [lo] + pts + [hi]
        tot = sum(integrate.quad(f, a, b, epsabs=0.0, epsrel=1e-12, limit=500)[0] for a, b in zip(edges[:-1], edges[1:]))
        return Om * Om * (1 - mu * mu) / 4.0 * tot

    xg, wg = np.polynomial.legendre.leggauss(30)
    edges = np.concatenate([[0.0], np.geomspace(1e-9, 1.0, 90)])
    for key, P in B['broadened_slab']['cases'].items():
        rs, e2, kern = P['r_L_over_ell_s'], P['eps_s2'], P['kernel']
        width = P['width_over_Omega'] / rs
        tot = 0.0
        for a, b in zip(edges[:-1], edges[1:]):
            z = 0.5 * (b - a) * xg + 0.5 * (b + a)
            tot += 0.5 * (b - a) * float(np.dot(wg, [3 * t * t * (1 - t ** 6) ** 2 / D(t ** 3, rs, e2, kern, width) for t in z]))
        lam = 0.75 * tot
        record('broadened slab ' + key, rel(lam, P['lambda_par_over_ell_s']) <= tol,
               'rel %.1e' % rel(lam, P['lambda_par_over_ell_s']))


# ------------------------------------------------------------------------- audit
def check_audit(dpar, dperp):
    S = json.load(open(path('fit_error_summary.json')))
    g = S['generator']
    lo = np.log([g['box'][k][0] for k in g['column_order']])
    hi = np.log([g['box'][k][1] for k in g['column_order']])
    X = np.random.default_rng(g['seed']).uniform(lo, hi, size=(300, 4))
    rows = list(csv.DictReader(open(path('fit_error_samples.csv'), newline='')))
    record('audit sample size', len(rows) == 300)
    draw_ok = all(abs(float(r_['ln_r_L_over_ell_s']) - X[n, 0]) == 0 and abs(float(r_['ln_f_s']) - X[n, 1]) == 0 and
                  abs(float(r_['ln_eps2_total']) - X[n, 2]) == 0 and abs(float(r_['ln_ell_s_over_ell_2']) - X[n, 3]) == 0
                  for n, r_ in enumerate(rows))
    record('audit draw reproduces fit_error_samples.csv inputs exactly', draw_ok)
    wF = wN = wres = 0.0
    ap, aq, eps2s = [], [], []
    for n, r_ in enumerate(rows):
        x = [float(t) for t in X[n]]
        st = tuple(math.exp(t) for t in x)
        Fp, Fq = math.exp(poly_F(dpar, x)), math.exp(poly_F(dperp, x))
        lp, lq, res, _ = solve(st, 'nlgce_n', [math.log(Fp / 3), math.log(Fq / 3)])
        wF = max(wF, rel(Fp, float(r_['F_lambda_par_over_ell_s'])), rel(Fq, float(r_['F_lambda_perp_over_ell_s'])))
        wN = max(wN, rel(lp, float(r_['N_lambda_par_over_ell_s'])), rel(lq, float(r_['N_lambda_perp_over_ell_s'])))
        wres = max(wres, res)
        ap.append(abs(Fp / lp - 1))
        aq.append(abs(Fq / lq - 1))
        eps2s.append(st[2])
    record('audit polynomial values', wF <= 1e-11, 'max rel %.1e' % wF)
    record('audit NLGCE-N values', wN <= 5e-8 and wres < 1e-8, 'max rel %.1e, max residual %.1e' % (wN, wres))
    ap, aq, eps2s = np.array(ap), np.array(aq), np.array(eps2s)
    print('\n eps2 bin          states  median par  >25% par  median perp  >25% perp')
    allok = True
    for b in S['results']['bins']:
        lo_, hi_ = b['eps2_lower'], b['eps2_upper']
        sel = (eps2s >= lo_) & ((eps2s < hi_) | (b['upper_inclusive'] & (eps2s == hi_)))
        mp_, sp_ = np.median(ap[sel]), np.mean(ap[sel] > 0.25)
        mq_, sq_ = np.median(aq[sel]), np.mean(aq[sel] > 0.25)
        print(' %-7g to %-7g %5d   %8.1f%%  %8.1f%%   %8.1f%%   %8.1f%%' % (lo_, hi_, sel.sum(), 100 * mp_, 100 * sp_, 100 * mq_, 100 * sq_))
        allok &= int(sel.sum()) == b['states'] and abs(mp_ - b['median_parallel_discrepancy']) < 1e-6 and \
            abs(mq_ - b['median_perpendicular_discrepancy']) < 1e-6 and abs(sp_ - b['parallel_share_above_0_25']) < 1e-12 \
            and abs(sq_ - b['perpendicular_share_above_0_25']) < 1e-12
    print()
    record('audit six-bin summary', allok)


def main():
    args = set(sys.argv[1:])
    check_digests()
    dpar = load_coeffs('NLGCE_F_2014_parallel.csv')
    dperp = load_coeffs('NLGCE_F_2014_perpendicular.csv')
    B = json.load(open(path('benchmark_points.json')))
    check_particles(B)
    check_spectrum(B)
    check_nlgce(B, dpar, dperp)
    if '--broadened' in args:
        check_broadened(B)
    if '--audit' in args:
        check_audit(dpar, dperp)
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print('\n%d checks, %d failed' % (len(RESULTS), nfail))
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    import warnings
    from scipy.integrate import IntegrationWarning
    warnings.simplefilter('ignore', IntegrationWarning)
    main()

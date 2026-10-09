#!/usr/bin/env python3
"""reference_verification.py -- MEAN_FREE_PATH_MODEL_DATA, specification version 1.2.

1. Recomputes every numerical value in benchmark_points.json independently of the
   fixture generator: closed forms are used where they exist (exact slab-QLT integrals
   as Gauss hypergeometric and piecewise-polynomial forms, the closed form of I(5/3,H),
   hypergeometric forms of the Droge integral, analytic logarithmic slopes), and the
   GCR inputs are read from the data files of this bundle.
2. Compares with relative tolerances (Section 13 of the specification).  Reference
   values that are 40-digit differences (keys 'relative_difference', 'rel_diff') are
   checked against absolute bounds instead.  Values printed with fewer digits
   (numerical slopes, extremum locations) carry the tolerance of their digits.
3. Checks the structure and provenance fields of the data files, the runtime states,
   that every equation label is resolved, and the SHA256SUMS digests.

Requires Python 3.9+ and mpmath.  Exit status 0 when every check passes.
This script does not verify that a transcription matches its source.
"""
import csv, hashlib, json, os, re, sys
import mpmath as mp

mp.mp.dps = 30
BASE = os.path.dirname(os.path.abspath(__file__))
def P(*a): return os.path.join(BASE, *a)
M = mp.mpf

# ---------------------------------------------------------------- constants (Section 2)
c = M(299792458); e = M('1.602176634e-19'); AU = M(149597870700); RSUN = M('6.957e8')
MASS = {'p': M('1.67262192369e-27'), 'e': M('9.1093837015e-31'), 'alpha': M('6.6446573357e-27')}
def rest(sp): return MASS[sp] * c**2 / (e * 10**6)                 # MeV
def pc(T, sp): T = M(T); return mp.sqrt(T * (T + 2 * rest(sp)))       # MeV (= rigidity in MV for |Z| = 1)
def beta(T, sp): T = M(T); return pc(T, sp) / (T + rest(sp))
def T_of_pc(p, sp): p = M(p); return mp.sqrt(p**2 + rest(sp)**2) - rest(sp)

# ---------------------------------------------------------------- data files used as inputs
def rows(fn): return list(csv.DictReader(open(P(fn), encoding='utf-8')))
NWU = rows('parameters/gcr_nwu_family_parameter_sets.csv')
def nwu(src, loc, label, par):
    hit = [r['value_as_printed'] for r in NWU if r['source'] == src and r['location'] == loc and r['label'] == label and r['parameter'] == par]
    if len(hit) != 1: raise SystemExit('data lookup failed: %s %s %s %s' % (src, loc, label, par))
    return M(hit[0])
HELMOD = json.load(open(P('parameters/gcr_helmod.json'), encoding='utf-8'))
OTHER = json.load(open(P('parameters/gcr_other_models.json'), encoding='utf-8'))
def central(s): return M(re.match(r'\s*([-+0-9.eE]+)', s).group(1))

# ---------------------------------------------------------------- closed forms
def I53(H):
    """I(5/3,H) = int_0^1 (1-mu^2)/(mu^(2/3)+H) dmu; mu = t^3 makes the integrand rational."""
    H = M(H)
    if H == 0: return M(18) / 7
    return 3 * (M(6) / 7 + H / 5 - H**2 / 3 + H**3 - (H + H**4) * mp.atan(1 / mp.sqrt(H)) / mp.sqrt(H))

def phi_eps(eps):
    """phi = (3/2) int_0^1 (1-mu)(1+mu)^2/(a mu + eps) dmu, a = 1 + eps, by polynomial division."""
    eps = M(eps); a = 1 + eps; m0 = -eps / a
    N = [M(-1), M(-1), M(1), M(1)]            # -mu^3 - mu^2 + mu + 1, highest power first
    q = [N[0]]
    for coef in N[1:]: q.append(coef + q[-1] * m0)
    rem = q.pop()                               # N(m0)
    deg = len(q) - 1
    intQ = sum(qk / (deg - k + 1) for k, qk in enumerate(q))
    return M(3) / 2 * (intQ / a + rem / a * mp.log((a + eps) / eps))

def droge_ratio(a):
    """int_0^1 (1-mu^2)(mu^2+a^2)^(-1/3) dmu / (18/7), via 2F1."""
    a = M(a); z = -1 / a**2; p = M(1) / 3
    J = a ** (-2 * p) * (mp.hyp2f1(p, M(1) / 2, M(3) / 2, z) - mp.hyp2f1(p, M(3) / 2, M(5) / 2, z) / 3)
    return J / (M(18) / 7)

def hewan(x):
    x = M(x)
    with mp.workdps(80):
        return 3 * (x - mp.tanh(x)) / x**3

NU = M(5) / 6
def zank_exact_ratio(x): return mp.hyp2f1(-NU, 1 - NU, 3 - NU, -M(x)**2)
def zank_aux(sz):
    s2 = M(sz)**2
    A = mp.expm1(M(5) / 6 * mp.log1p(s2))
    q = M(2) if s2 == 0 else (M(5) / 3 * s2) / (s2 - mp.expm1(mp.log1p(s2) / 6))
    return A, q, 1 + (M(7) / 9 * A) / ((q + M(1) / 3) * (q + M(7) / 3))

def ts_J(R, s):
    R = M(R); s = M(s)
    if R <= 1: return 2 * R ** (-s) / ((2 - s) * (4 - s))
    return M(1) / 4 + (1 / (2 - s) - M(1) / 2) / R**2 + (M(1) / 4 - 1 / (4 - s)) / R**4
def ts_ratio(R, s):
    R = M(R); s = M(s)
    return ts_J(R, s) / (M(1) / 4 + 2 * R ** (-s) / ((2 - s) * (4 - s)))

def nwu_slope(Pv, a, b, cc, Pk):          # d ln G / d ln P for Eq. (33)
    Pv, a, b, cc, Pk = map(M, (Pv, a, b, cc, Pk)); return a + (b - a) * Pv**cc / (Pv**cc + Pk**cc)
def G(Pv, a, b, cc, Pk):
    Pv, a, b, cc, Pk = map(M, (Pv, a, b, cc, Pk)); return Pv**a * ((Pv**cc + Pk**cc) / (1 + Pk**cc)) ** ((b - a) / cc)
def duan_slope(R, a, b, cc, Rk):          # Eq. (47)
    R, a, b, cc, Rk = map(M, (R, a, b, cc, Rk)); y = (R / Rk) ** ((b - a) / cc); return a + (b - a) * y / (1 + y)
def corti_slope(R, a, b, s, Rk):
    R, a, b, s, Rk = map(M, (R, a, b, s, Rk)); x = (R / Rk) ** s; return a + (b - a) * x / (1 + x)

# ---------------------------------------------------------------- TS / Lang evaluations (Eqs. (21), (22), (23))
def lam_ts(T, sp, pr, branch, variant=None):
    B0 = M(pr['B0_nT']) * M('1e-9'); dB2 = M(pr['dB2_nT2']); s = M(pr['s']); p = M(pr['p'])
    kmin = M(pr['kmin_km']) / 1000; kd = M(pr['kd_km']) / 1000; VA = M(pr['VA_kms']) * 1000; aD = M(pr['aD'])
    RL = pc(T, sp) * 10**6 / (c * B0); R = RL * kmin; Q = RL * kd
    base = 3 * s * RL**2 * kmin / (4 * mp.pi * (s - 1)) * M(pr['B0_nT'])**2 / dB2
    br = 1 + 8 / ((2 - s) * (4 - s) * R**s)
    if branch == 'none': return base * br / AU
    a = beta(T, sp) * c / (aD * VA)
    if branch == 'RS':
        K = (mp.sqrt(mp.pi) / mp.gamma(p / 2) + 1 / (p - 2)) * (a / 2) ** (p - 2); Qp = Q ** (p - s)
    else:
        f1 = 2 * (p - s) / (mp.pi * (p - 2) * (2 - s))
        arg = -a / (f1 * Q ** (p - 2)) if variant == 'EB13' else -a / (f1 * Q)
        K = mp.hyp2f1(1, 1 / (p - 1), p / (p - 1), arg) * a / f1
        Qp = Q ** (3 - s) if variant == 'TS-Q3MS' else Q ** (p - s)
    return base * (br + 4 * K / (Qp * R**s)) / AU

# ---------------------------------------------------------------- expected outputs per fixture
def expected(fx):
    i = fx['id']; o = fx['outputs']; inp = fx['inputs']
    if i == 'F-KIN-01':
        return {'proton_1MeV_R_MV': pc(1, 'p'), 'proton_1MeV_beta': beta(1, 'p'), 'proton_pc_1GeV_T_MeV': T_of_pc(1000, 'p'),
                'electron_0.094MeV_R_MV': pc('0.094', 'e'), 'electron_0.94MeV_R_MV': pc('0.94', 'e'),
                'electron_R_0.324MV_T_MeV': T_of_pc('0.324', 'e'), 'electron_2.0MeV_R_MV': pc(2, 'e'),
                'electron_2.5MeV_R_MV': pc('2.5', 'e'), 'alpha_Z2_pc_1GeV_R_GV': M('0.5'), 'alpha_pc_1GeV_T_total_MeV': T_of_pc(1000, 'alpha')}
    if i == 'F-SEP-01':
        return {'SEP-PATH09_lambda_au': M('0.8'), 'SEP-EPREM13_lambda_par_au': M('0.05'), 'SEP-MFLAMPA25_lambda_par_au': M('0.3'),
                'SEP-ZHANG23_lambda_r_au': 200 * RSUN / AU, 'SEP-ZHANG23_200Rsun_m': 200 * RSUN}
    if i == 'F-SEP-02':
        imp = {}
        for k in o['implied_cos2psi']:
            a, b = re.search(r'\(([0-9.]+), ([0-9.]+)\)', k).groups(); imp[k] = M(a) / M(b)
        t = M('2.66e-6') * (AU - M('0.005') * AU) / M(400000)
        return {'lambda_par_over_lambda_r_at_psi45': M(2), 'implied_cos2psi': imp, 'Parker_example_tan_psi': t,
                'Parker_example_psi_deg': mp.degrees(mp.atan(t)), 'Parker_example_cos2psi': 1 / (1 + t**2)}
    if i == 'F-SEP-03':
        out = []
        for r in o['rows']:
            T = M(repr(r['T_MeV'])); E0 = rest('p')
            lam = M('0.3') * (pc(T, 'p') / 1000) ** (M(1) / 3)
            D1 = beta(T, 'p') * c * lam * AU / 3
            D2 = c * M('0.3') * AU / 3 * (T * (T + 2 * E0) / M(10)**6) ** (M(1) / 6) * mp.sqrt(T * (T + 2 * E0) / (T + E0)**2)
            out.append({'T_MeV': T, 'lambda_par_au': lam, 'D_par_m2_s_from_v_lambda_over_3': D1, 'D_par_m2_s_Eq16': D2, 'relative_difference': abs(D1 / D2 - 1)})
        return {'rows': out}
    if i == 'F-SEP-05':
        R1 = pc(1, 'p')
        out = [{'T_MeV': M(repr(r['T_MeV'])), 'R_MV': pc(repr(r['T_MeV']), 'p'),
                'lambda0_au_Rref43': M('0.1') * (pc(repr(r['T_MeV']), 'p') / 43) ** (M(1) / 3),
                'lambda0_au_Rref_exact': M('0.1') * (pc(repr(r['T_MeV']), 'p') / R1) ** (M(1) / 3)} for r in o['rows']]
        return {'R_1MeV_proton_MV': R1, 'ratio_(R1/43)^(1/3)': (R1 / 43) ** (M(1) / 3), 'rows': out}
    if i == 'F-SEP-06':
        out = []
        for r in o['rows']:
            rr = M(repr(r['r_au'])); E = M(repr(r['E_keV'])); T = E / 1000
            kap = M('5.16e18') * rr ** M('1.17') * E ** M('0.71')
            out.append({'r_au': rr, 'E_keV': E, 'kappa_par_cm2_s': kap,
                        'derived_lambda_par_au_proton': 3 * kap / (beta(T, 'p') * c * 100) / (AU * 100),
                        'MFLAMPA25_lambda_par_au': M('0.3') * rr * (pc(T, 'p') / 1000) ** (M(1) / 3)})
        return {'rows': out}
    if i == 'F-PA-01':
        return {'I_exact_H0_18_over_7': M(18) / 7, 'values': {H: {'I': I53(H), 'D0_lambda_over_v': 3 * I53(H) / 4} for H in o['values']}}
    if i == 'F-PA-02':
        return {'values': {k: {'phi_exact': phi_eps(k), 'asymptote_1.5_ln(1/eps)': M(3) / 2 * mp.log(1 / M(k)),
                               'ratio_exact_over_asymptote': phi_eps(k) / (M(3) / 2 * mp.log(1 / M(k)))} for k in o['values']}}
    if i == 'F-PA-03':
        return {'ratio': {k: droge_ratio(k) for k in o['ratio']}}
    if i == 'F-PA-04':
        return {'lambda_effective_over_lambda': M(2)}
    if i == 'F-PA-05':
        tp = (2 * mp.pi) ** (-M(2) / 3)
        return {'lambda_xx_over_lambda_mumu_numeric': M(27) / 14, 'exact_27_over_14': M(27) / 14,
                'printed_ratio_(54/7pi)/(4/pi)': (54 / (7 * mp.pi)) / (4 / mp.pi), '81/(7pi)': 81 / (7 * mp.pi),
                'lambda_xx_coeff_with_k0=2pi/Lmax_and_int_I=w': 81 / (7 * mp.pi) * tp,
                'lambda_mumu_coeff_same_assumptions': 6 / mp.pi * tp, 'printed_approximate_coefficients': [M('0.9'), M('0.5')]}
    if i == 'F-PA-06':
        return {'lambda_over_lambda0': {x: hewan(x) for x in o['lambda_over_lambda0']}}
    if i == 'F-QLT-01':
        B0 = M('4.12e-9'); kmin = M('1e-10'); s = M(5) / 3
        RL = M(10)**6 / (c * B0); R = RL * kmin
        mid02 = 6 * R ** (2 - s) / (mp.pi * (s - 1) * (2 - s) * (4 - s) * kmin) / AU
        return {'TS2002_mid_coeff_AU': mid02, 'printed_TS2002': M('0.0106'), 'TS2003_mid_coeff_AU': s * mid02, 'printed_TS2003': M('0.018'),
                'TS2003_high_coeff_AU': 3 * s * kmin * RL**2 / (4 * mp.pi * (s - 1)) / AU, 'printed_TS2003_high': M('2.62e-10'),
                'P_at_R_eq_1_MV': c * B0 / kmin / M(10)**6, 'printed_limit_MV': M('1.23e4')}
    if i == 'F-QLT-02':
        C = mp.gamma(NU) / (2 * mp.sqrt(mp.pi) * mp.gamma(NU - M(1) / 2)); s = 2 * NU
        cl = 3 / (2 * mp.pi * C * (2 - s) * (4 - s)); ct = cl * (2 * mp.pi * C) ** (-M(2) / 3)
        acc = {}
        for x in o['accuracy_vs_RL_over_l']:
            if x == 'max': continue
            F = zank_exact_ratio(x); br = zank_aux(x)[2]
            acc[x] = {'exact_over_inertial': F, 'zank_bracket_at_s_eq_x': br, 'bracket_over_exact': br / F}
        g = lambda y: zank_aux(y)[2] / zank_exact_ratio(y)
        xm = mp.findroot(lambda y: mp.diff(g, y), M('5.8'))
        acc['max'] = {'R_L_over_l': xm, 'bracket_over_exact': g(xm)}
        return {'C(5/6)': C, '2 pi C': 2 * mp.pi * C, 'inertial_coeff_vs_l_total_variance': cl, 'coeff_vs_lambda_s_total_variance': ct,
                'coeff_vs_lambda_s_per_component': ct / 2, 'printed': {'Zank1998_per_component': M('3.1371'), 'Chhiber2017_total': M('6.2742'), 's_constant': M('0.746834')},
                'accuracy_vs_RL_over_l': acc}
    if i == 'F-QLT-03':
        return {'rows': [{'T_MeV': M(repr(r['T_MeV'])), 'lambda_par_au': lam_ts(repr(r['T_MeV']), 'p', inp, 'none')} for r in o['rows']]}
    if i == 'F-QLT-04':
        out = []
        for r in o['lambda_par_au_rows']:
            T = repr(r['T_MeV'])
            out.append({'T_MeV': M(T), 'DT_LANG24': lam_ts(T, 'e', inp['DT'], 'DT', 'LANG24'), 'DT_TS-Q3MS': lam_ts(T, 'e', inp['DT'], 'DT', 'TS-Q3MS'),
                        'DT_EB13': lam_ts(T, 'e', inp['DT'], 'DT', 'EB13'), 'RS': lam_ts(T, 'e', inp['RS'], 'RS'),
                        'no_dissipation_term_DTparams': lam_ts(T, 'e', inp['DT'], 'none')})
        return {'lambda_par_au_rows': out}
    if i == 'F-QLT-05':
        km = M(repr(inp['k_min_per_km'])); kd = M(repr(inp['k_d_per_km'])); s = M(repr(inp['s'])); p = M(repr(inp['p']))
        return {'lambda_slab_km': mp.pi * (s - 1) / (2 * km) / (s + (s - p) / (p - 1) * (km / kd) ** (s - 1)),
                'printed_Table5_lambda_slab_km': M('2.81e6'), 'k_min_to_k_d_term_limit_km': mp.pi * (s - 1) / (2 * s * km)}
    if i == 'F-QLT-06':
        s = M('1.65'); v = {}
        for k in o['values']:
            p = M(k)
            v[k] = {'2^(p-1) G((p+1)/2)/(sqrt(pi) G(p))': 2 ** (p - 1) * mp.gamma((p + 1) / 2) / (mp.sqrt(mp.pi) * mp.gamma(p)),
                    '1/G(p/2)': 1 / mp.gamma(p / 2), 'f1_TS/pi': (2 / (p - 2) + 2 / (2 - s)) / mp.pi,
                    'f1_Lang(s=1.65)': 2 * (p - s) / (mp.pi * (p - 2) * (2 - s))}
        return {'values': v}
    if i == 'F-QLT-07':
        s = M(5) / 3
        Rm = mp.findroot(lambda y: mp.diff(lambda z: ts_ratio(z, s), y), M(3))
        return {'exact_over_form_vs_R': {R: ts_ratio(R, s) for R in o['exact_over_form_vs_R']},
                'minimum': {'R': Rm, 'exact_over_form': ts_ratio(Rm, s)}}
    if i == 'F-SH-01':
        out = []
        for r in o['rows']:
            T = repr(r['T_MeV']); rg = pc(T, 'p') * 10**6 / (c * M('5e-9'))
            out.append({'T_MeV': M(T), 'B_nT': M(5), 'r_g_m': rg, 'r_g_au': rg / AU, 'kappa_Bohm_m2_s': beta(T, 'p') * c * rg / 3})
        return {'rows': out}
    if i == 'F-SH-02':
        return {'lambda_from_nu_over_printed_lambda': M(1)}
    if i == 'F-SH-03':
        out = []
        for r in o['rows']:
            krr = M(286000) * M(556000) * M(repr(r['Delta_t_min'])) * 60
            out.append({'channel': r['channel'], 'Delta_t_min': M(repr(r['Delta_t_min'])), 'kappa_rr_cm2_s': krr * 10**4,
                        'kappa_par_cm2_s': krr / mp.cos(mp.radians(66))**2 * 10**4})
        return {'rows': out}
    if i == 'F-SH-04':
        return {'D_min_m2_s': M('0.1') * RSUN * 10**5, 'D_min_cm2_s': M('0.1') * RSUN * 10**9}
    if i == 'F-GCR-01':
        return {'G(1 GV)': G(1, '0.56', '1.95', '3.0', '4.0'),
                'logarithmic_slope_at_P_GV': {k: nwu_slope(k, '0.56', '1.95', '3.0', '4.0') for k in o['logarithmic_slope_at_P_GV']}}
    if i == 'F-GCR-02':
        by = {}
        for y in o['by_year']:
            lam = nwu('Potgieter2014', 'Table 1', y, 'lambda_par_Earth_100MV'); a = nwu('Potgieter2014', 'Table 1', y, 'a')
            Pk = nwu('Potgieter2014', 'Table 1', y, 'P_k'); b = nwu('Potgieter2014', 'Sec. 3 (Eq. 5)', y, 'b'); cc = nwu('Potgieter2014', 'Sec. 3 (Eq. 5)', y, 'c')
            g = G('0.1', a, b, cc, Pk); by[y] = {'G(0.1 GV)': g, 'implied_lambda_par_1GV_au': lam / g}
        return {'by_year': by, 'printed_text_Sec5': 'arXiv v3 (final): increased by a factor of ~2.3, from ~0.13 AU in 2006, to ~0.3 AU in 2009; v1-v2 printed ~23 and ~30 AU'}
    if i == 'F-GCR-03':
        a, b, s, Rk = M('0.8'), M('1.7'), M('2.2'), M('4.3'); k0 = Rk**b / (1 + Rk**s) ** ((b - a) / s); out = []
        for r in o['rows']:
            R = M(repr(r['R_GV'])); n = G(R, a, b, s, Rk); co = k0 * (R / Rk) ** a * (1 + (R / Rk) ** s) ** ((b - a) / s)
            out.append({'R_GV': R, 'NWU_shape': n, 'Corti_shape_with_converted_k0': co, 'rel_diff': abs(n / co - 1)})
        return {'k0_Corti_over_K0_NWU': k0, 'rows': out}
    if i == 'F-GCR-04':
        by = {}
        for y in o['by_period']:
            lam = nwu('VosPotgieter2015', 'Table 2', y, 'lambda_par_Earth_1GV'); B = nwu('VosPotgieter2015', 'Table 1', y, 'B_e')
            K = lam * AU * 100 * c * 100 * B / 3; by[y] = {'implied_(K_par)_0_cm2_s': K, 'in_units_of_1e22': K / M(10)**22}
        return {'by_period': by}
    if i == 'F-GCR-05':
        AU2 = AU**2 * 10**4; T1 = HELMOD['Boschini2018ASR']['K0_SSN']['rows']
        rowmap = {'A<0 asc': 'A<0 ascending', 'A<0 desc': 'A<0 descending', 'A>0 asc': 'A>0 ascending', 'A>0 desc': 'A>0 descending'}
        ssn = {}
        for k, src in rowmap.items():
            cs = [M(T1[src][n]) if T1[src][n] is not None else M(0) for n in ('c0', 'c1', 'c2', 'c3')]
            ssn[k] = {S: cs[0] + cs[1] * M(S) + cs[2] * M(S)**2 + cs[3] * M(S)**3 for S in o['K0_SSN_by_row_and_SSN'][k]}
        mc = HELMOD['Boschini2018ASR']['K0_NMCR']['rows']['MCMU']
        nm = {N: M(mc['p0']) * mp.exp(M(mc['p1']) * M(N) + M(mc['p2']) * M(N)**2) for N in o['K0_NMCR_MCMU']}
        c0 = M(T1['A<0 ascending']['c0']); kp = []
        for r in o['K_par_examples']:
            Pv = M(repr(r['P_GV'])); bt = Pv * 1000 / mp.sqrt((Pv * 1000)**2 + rest('p')**2); K = bt / 3 * c0 * (Pv + M('0.3')) * 2
            kp.append({'P_GV': Pv, 'beta': bt, 'K_par_AU2_s': K, 'K_par_cm2_s': K * AU2})
        return {'1 AU^2/s in cm^2/s': AU2, 'c0_A<0asc_cm2_s_GV': c0 * AU2, 'K0_SSN_by_row_and_SSN': ssn, 'K0_NMCR_MCMU': nm, 'K_par_examples': kp}
    if i == 'F-GCR-06':
        tab = OTHER['GCR-TOMASSETTI17']['table_II']
        return {'kappa0': {st: {ph: central(tab[st]['a_MV']) / M(ph) + central(tab[st]['b']) for ph in o['kappa0'][st]} for st in o['kappa0']}}
    if i == 'F-GCR-07':
        st = [{'P_GV': M(repr(r['P_GV'])), 'r_au': M(repr(r['r_au'])), 'lambda_par_au': M('0.15') * max(M(repr(r['P_GV'])), M(1)) * (1 + M(repr(r['r_au'])))} for r in o['Strauss2011']]
        return {'Strauss2011': st, 'Wang2019_continuity_beta_B_factor_1': {'R=0.1 GV low branch (k units)': M(1) / 30, 'R=0.1 GV high branch (k units)': M('0.1') / 3}}
    if i in ('F-GCR-08', 'F-GCR-09'):
        a, b = (M('0.8'), M('1.7')) if i == 'F-GCR-08' else (M('1.7'), M('0.8')); cc = M('2.2'); Rk = M('4.3')
        if i == 'F-GCR-08':
            sl = {R: {'Duan': duan_slope(R, a, b, cc, Rk), 'Corti': corti_slope(R, a, b, cc, Rk)} for R in o['log_slopes']}
        else:
            sl = {R: {'Duan_numerical': duan_slope(R, a, b, cc, Rk), 'Duan_slope_formula': duan_slope(R, a, b, cc, Rk),
                      'Corti_numerical': corti_slope(R, a, b, cc, Rk)} for R in o['log_slopes']}
        res = {'log_slopes': sl, 'value_at_Rk': {'Duan': 2**cc, 'Corti': 2 ** ((b - a) / cc)}}
        if i == 'F-GCR-09': res['limits'] = {'Duan_low': b, 'Duan_high': a, 'Corti_low': a, 'Corti_high': b}
        return res
    if i == 'F-NUM-01':
        hw = {}
        for x in o['hewan_series']:
            X = M(x)
            hw[x] = {'lambda_over_lambda0': hewan(x), 'series_through_x8': 1 - 2 * X**2 / 5 + 17 * X**4 / 105 - 62 * X**6 / 945 + 1382 * X**8 / M(51975),
                     'first_omitted_term': M(21844) / 2027025 * X**10}
        za = {sv: dict(zip(('A', 'q', 'bracket'), zank_aux(sv))) for sv in o['zank_aux_small_s']}
        zx = {x: {'2F1(-nu,1-nu;3-nu;-x^2)': zank_exact_ratio(x), 'quadrature': zank_exact_ratio(x)} for x in o['zank_exact_over_inertial']}
        s = M(5) / 3
        tx = {R: {'J_tilde': ts_J(R, s), 'exact_over_form_closed': ts_ratio(R, s), 'exact_over_form_quadrature': ts_ratio(R, s)} for R in o['ts_exact_closed_form']}
        Rm = mp.findroot(lambda y: mp.diff(lambda z: ts_ratio(z, s), y), M(3))
        return {'hewan_series': hw, 'zank_aux_small_s': za, 'zank_q_limit_s_to_0': M(2), 'zank_exact_over_inertial': zx,
                'ts_exact_closed_form': tx, 'ts_minimum_closed_form': {'R': Rm, 'exact_over_form': ts_ratio(Rm, s)}}
    if i == 'F-PERP-01':
        return {'<sqrt(1-mu^2)>': mp.pi / 4, 'pi/4': mp.pi / 4}
    return None

# ---------------------------------------------------------------- comparison
fails = []; n_ok = 0
count = {'derived': 0, 'coordinate': 0, 'printed': 0}
COORD = {'T_MeV', 'r_au', 'E_keV', 'R_GV', 'P_GV', 'Delta_t_min', 'B_nT'}
PRINTED = ('printed', 'SEP-PATH09_lambda_au', 'SEP-EPREM13_lambda_par_au', 'SEP-MFLAMPA25_lambda_par_au', '/limits/')
def ok(cond, msg):
    global n_ok
    if cond: n_ok += 1
    else: fails.append(msg)

def rtol_for(path):
    if '/log_slopes/' in path or '/logarithmic_slope_at_P_GV/' in path or 'first_omitted_term' in path: return M('1e-4')
    if path.endswith('/minimum/R') or path.endswith('/max/R_L_over_l') or path.endswith('ts_minimum_closed_form/R'): return M('1e-7')
    return M('1e-12')

def compare(stored, exp, path):
    if isinstance(stored, dict):
        ok(isinstance(exp, dict) and set(stored) == set(exp), 'keys differ at %s: %s' % (path, sorted(set(stored) ^ set(exp or {}))))
        if isinstance(exp, dict):
            for k in stored:
                if k in exp: compare(stored[k], exp[k], path + '/' + k)
    elif isinstance(stored, list):
        ok(isinstance(exp, list) and len(stored) == len(exp), 'list length differs at ' + path)
        if isinstance(exp, list):
            for j, (a, b) in enumerate(zip(stored, exp)): compare(a, b, path + '[%d]' % j)
    elif isinstance(stored, str):
        ok(stored == exp, 'string differs at %s' % path)
    else:
        key = path.rsplit('/', 1)[-1].split('[')[0]
        count['coordinate' if key in COORD else 'printed' if any(t in path for t in PRINTED) else 'derived'] += 1
        if 'relative_difference' in path or 'rel_diff' in path:
            ok(abs(stored) <= 1e-30 and abs(exp) <= M('1e-20'), 'difference not negligible at %s: stored %r, recomputed %s' % (path, stored, mp.nstr(exp, 5)))
        else:
            ref = M(repr(stored)); tol = rtol_for(path)
            good = (abs(M(exp) - ref) <= tol * abs(ref)) if ref != 0 else (abs(M(exp)) <= M('1e-30'))
            ok(good, 'value differs at %s: stored %r, recomputed %s' % (path, stored, mp.nstr(M(exp), 17)))

bp = json.load(open(P('benchmark_points.json')))
EXPECTED_IDS = ['F-KIN-01', 'F-SEP-01', 'F-SEP-02', 'F-SEP-03', 'F-SEP-05', 'F-SEP-06', 'F-PA-01', 'F-PA-02', 'F-PA-03', 'F-PA-04', 'F-PA-05',
                'F-PA-06', 'F-QLT-01', 'F-QLT-02', 'F-QLT-03', 'F-QLT-04', 'F-QLT-05', 'F-QLT-06', 'F-QLT-07', 'F-SH-01', 'F-SH-02', 'F-SH-03',
                'F-SH-04', 'F-GCR-01', 'F-GCR-02', 'F-GCR-03', 'F-GCR-04', 'F-GCR-05', 'F-GCR-06', 'F-GCR-07', 'F-GCR-08', 'F-GCR-09',
                'F-NUM-01', 'F-PERP-01']
ok([x['id'] for x in bp['fixtures']] == EXPECTED_IDS, 'fixture identifiers differ from the 34 groups of Section 13')
for fx in bp['fixtures']:
    ok(all(k in fx for k in ('id', 'title', 'kind', 'inputs', 'outputs', 'reference', 'note')), fx['id'] + ': missing fixture field')
    exp = expected(fx)
    ok(exp is not None, fx['id'] + ': no independent recomputation')
    if exp is not None: compare(fx['outputs'], exp, fx['id'])
print('fixture output values compared: %d (34 groups): %d derived values recomputed without the generator, '
      '%d row coordinates, %d printed or illustrative constants' % (sum(count.values()), count['derived'], count['coordinate'], count['printed']))

# ---------------------------------------------------------------- data files
srcmap = json.load(open(P('source_key_map.json')))
refs = srcmap['references']
eqnums = srcmap['equation_numbers']
ok(sorted(int(v) for v in eqnums.values()) == list(range(1, len(eqnums) + 1)), 'equation numbers are not 1..N')
for fn, n, need in [('observations/event_fitted_mfp.csv', 114, ['source', 'location', 'quantity', 'value_au', 'verification']),
                    ('parameters/gcr_nwu_family_parameter_sets.csv', 461, ['source', 'location', 'parameter', 'value_as_printed', 'verification']),
                    ('turbulence/published_turbulence_values.csv', 139, ['source', 'location', 'quantity', 'value_as_printed', 'verification'])]:
    R = rows(fn)
    ok(len(R) == n, '%s: %d rows, expected %d' % (fn, len(R), n))
    for j, r in enumerate(R):
        for col in need: ok(r.get(col, '') != '', '%s row %d: empty %s' % (fn, j + 2, col))
        ok(r['source'] in refs, '%s row %d: source %s not in source_key_map' % (fn, j + 2, r['source']))
STATES = {'READY_EXPLICIT_INPUTS', 'READY_PUBLISHED_VARIANT', 'READY_PUBLISHED_KAPPA_ONLY', 'REQUIRES_EXTERNAL_INPUT',
          'REQUIRES_USER_DECISION', 'REQUIRES_SOURCE_OR_CODE_AUDIT', 'REFERENCE_DATA_ONLY'}
J = {}
for fn in ['parameters/sep_code_presets.json', 'parameters/gcr_helmod.json', 'parameters/gcr_other_models.json',
           'parameters/qlt_parameter_sets.json', 'parameters/model_registry.json']:
    try:
        J[fn] = json.load(open(P(fn), encoding='utf-8')); ok(True, '')
    except Exception as ex:
        ok(False, fn + ': ' + str(ex))
pre = J['parameters/sep_code_presets.json']['presets']
ok(len(pre) == 22 and len({p['id'] for p in pre}) == 22, 'sep_code_presets: expected 22 unique presets')
for p in pre:
    ok(p['source'] in refs, 'preset source not in key map: ' + p['source'])
    ok(p.get('runtime_state') in STATES, 'preset runtime_state invalid: ' + p['id'])
    for k in ('location', 'form', 'quantity', 'verification'): ok(p.get(k, '') != '', 'preset %s: empty %s' % (p['id'], k))
reg = J['parameters/model_registry.json']['models']
ok(len(reg) == 57 and len({m['id'] for m in reg}) == 57, 'model_registry: expected 57 unique IDs')
for m in reg:
    ok(m.get('runtime_state') in STATES, 'registry runtime_state invalid: ' + m['id'])
    ok(m.get('status', '') != '' and m.get('section', '') != '', 'registry entry incomplete: ' + m['id'])
    for s_ in re.findall(r'\b[A-Z][A-Za-z]+[0-9]{4}[A-Za-z]*\b', m['source']): ok(s_ in refs, 'registry source not in key map: %s (%s)' % (s_, m['id']))
ok({p['id'] for p in pre} <= {m['id'] for m in reg}, 'presets missing from registry')
KEYTOK = re.compile(r'\b[A-Z][A-Za-z]+[0-9]{4}[A-Za-z]*\b')
LOCWORD = re.compile(r'(Eq|Table|Sec|Fig|text|Appendix)')
def has_status(v): return bool(v.get('status')) or any(isinstance(x, dict) and x.get('status') for x in v.values())
for fn in ['parameters/gcr_other_models.json', 'parameters/qlt_parameter_sets.json', 'parameters/gcr_helmod.json']:
    for k, v in J[fn].items():
        if k.startswith('_'): continue
        src = k if fn.endswith('gcr_helmod.json') else v.get('source', '')
        toks = KEYTOK.findall(src)
        ok(len(toks) > 0 and all(t in refs for t in toks), '%s %s: source key missing or not in key map' % (fn, k))
        ok(has_status(v), '%s %s: no verification status' % (fn, k))
        ok(bool(v.get('location') or v.get('equation') or v.get('table') or LOCWORD.search(v.get('source', ''))), '%s %s: no location' % (fn, k))

# every equation label resolved; no stray label syntax in data or text files
for root, _, fns in os.walk(BASE):
    for f in fns:
        if f.endswith(('.json', '.csv', '.md')):
            txt = open(os.path.join(root, f), encoding='utf-8').read()
            ok(re.search(r'\bE:[A-Za-z]', txt) is None, 'unresolved equation label in ' + os.path.relpath(os.path.join(root, f), BASE))

# ---------------------------------------------------------------- digests
listed = set()
for line in open(P('SHA256SUMS')):
    h, fn = line.split(None, 1); fn = fn.strip(); listed.add(fn)
    ok(hashlib.sha256(open(P(fn), 'rb').read()).hexdigest() == h, 'digest mismatch: ' + fn)
present = {os.path.relpath(os.path.join(r, f), BASE) for r, _, fs in os.walk(BASE) for f in fs if f != 'SHA256SUMS' and '__pycache__' not in r}
ok(listed == present, 'SHA256SUMS does not list exactly the files present: %s' % sorted(listed ^ present))

print('checks passed: %d, failed: %d' % (n_ok, len(fails)))
for msg in fails[:60]: print('FAIL:', msg)
sys.exit(1 if fails else 0)

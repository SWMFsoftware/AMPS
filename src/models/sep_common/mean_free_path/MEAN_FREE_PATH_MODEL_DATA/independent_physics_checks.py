#!/usr/bin/env python3
"""independent_physics_checks.py -- MEAN_FREE_PATH_MODEL_DATA, specification version 1.2.

Checks the mathematical identities that the specification uses (Sections 6-10 and 21)
against numerical quadrature, numerical differentiation and a seeded stochastic
simulation.  Each check prints PASS or FAIL with the measured error and the tolerance.
These checks test mathematics only; they do not verify that any transcription matches
its source.  Requires Python 3.9+ and mpmath.  Exit status 0 when every check passes.
"""
import math, random, sys
import mpmath as mp

mp.mp.dps = 30
M = mp.mpf
results = []
def check(name, err, tol, detail=''):
    passed = bool(err <= tol)
    results.append(passed)
    print('%s  %-62s error %-10s tolerance %-8s %s' % ('PASS' if passed else 'FAIL', name, mp.nstr(M(err), 3), mp.nstr(M(tol), 2), detail))
def rel(a, b): return abs(M(a) / M(b) - 1)

def quad_sym(f):        # integral over [-1, 1] with break points at -1, 0, 1
    return mp.quad(f, [-1, -M('0.5'), 0, M('0.5'), 1])
def quad_pow0(g, alpha, pts=None):
    """int_0^1 mu^alpha g(mu) dmu for alpha > -1, with mu = t^(1/(1+alpha)) to remove the endpoint singularity"""
    k = 1 / (1 + alpha)
    return k * mp.quad(lambda t: g(t ** k), pts or [0, M('0.5'), 1])

# 1. q-form normalization, Eq. (13), inserted in Eq. (11)
def I_closed(q, H):
    q = M(q); H = M(H)
    if H == 0: return 2 / ((2 - q) * (4 - q))
    assert q == M(5) / 3
    return 3 * (M(6) / 7 + H / 5 - H**2 / 3 + H**3 - (H + H**4) * mp.atan(1 / mp.sqrt(H)) / mp.sqrt(H))
v = M(1); lam = M('0.3')
for q, H in [(M(5) / 3, 0), (M(5) / 3, '0.01'), (M(5) / 3, '0.05'), (M(5) / 3, '0.2'), (M('1.3'), 0), (M('1.5'), 0), (M('1.9'), 0)]:
    H = M(H); D0 = 3 * v * I_closed(q, H) / (4 * lam)
    if H == 0:   # even integrand; mu^(1-q) singular at 0
        lam_num = 3 * v / 4 * quad_pow0(lambda m: (1 - m**2) / D0, 1 - q)
    else:
        lam_num = 3 * v / 8 * quad_sym(lambda m: (1 - m**2)**2 / (D0 * (1 - m**2) * (abs(m) ** (q - 1) + H)))
    check('q-form: lambda from D0 (q=%s, H=%s)' % (mp.nstr(q, 4), mp.nstr(H, 3)), rel(lam_num, lam), M('1e-10'))

# 2. EPREM operator factor (Section 6.6): D = (1-mu^2) v/(2 lambda) in d_mu[(D/2) d_mu f]
lam_eff = 3 * v / 8 * quad_sym(lambda m: (1 - m**2)**2 / ((1 - m**2) * v / (2 * lam) / 2))
check('EPREM half-D operator: lambda_eff / lambda = 2', abs(lam_eff / lam - 2), M('1e-12'))

# 3. Duan slope, Eq. (47), both orderings, against numerical differentiation
def lnk2(lnR, a, b, c, Rk):
    x = mp.exp(lnR) / Rk
    return a * mp.log(x) + c * mp.log(1 + x ** ((b - a) / c))
def slope_eq(R, a, b, c, Rk):
    y = (R / Rk) ** ((b - a) / c); return a + (b - a) * y / (1 + y)
worst = M(0)
for a, b in [(M('0.8'), M('1.7')), (M('1.7'), M('0.8'))]:
    for c in [M(1), M('2.2')]:
        for R in ['0.01', '1', '4.3', '30', '1000']:
            R = M(R); d = mp.diff(lambda t: lnk2(t, a, b, c, M('4.3')), mp.log(R))
            worst = max(worst, abs(d - slope_eq(R, a, b, c, M('4.3'))))
check('Duan: Eq. (47) versus numerical derivative', worst, M('1e-20'))
for a, b, lo, hi in [(M('0.8'), M('1.7'), M('0.8'), M('1.7')), (M('1.7'), M('0.8'), M('0.8'), M('1.7'))]:
    e_lo = abs(slope_eq(M('1e-40'), a, b, M('2.2'), M('4.3')) - lo); e_hi = abs(slope_eq(M('1e40'), a, b, M('2.2'), M('4.3')) - hi)
    check('Duan limits a=%s b=%s: low slope %s, high slope %s' % (a, b, lo, hi), max(e_lo, e_hi), M('1e-6'))
def corti_slope(R, a, b, s, Rk):
    x = (R / Rk) ** s; return a + (b - a) * x / (1 + x)
e_c = max(abs(corti_slope(M('1e-40'), M('1.7'), M('0.8'), M('2.2'), M('4.3')) - M('1.7')), abs(corti_slope(M('1e40'), M('1.7'), M('0.8'), M('2.2'), M('4.3')) - M('0.8')))
check('Corti form keeps a below and b above for b < a', e_c, M('1e-6'))

# 4. Zank exact integral, Eq. (44), against quadrature, and the Pfaff form
worst = M(0); worst_pf = M(0)
for nu in [M('0.7'), M(5) / 6, M('0.9')]:
    s = 2 * nu
    for x in ['0.01', '1', '5.787', '100']:
        x = M(x)
        # mu (1 - mu^2)(1 + 1/(mu x)^2)^nu = mu^(1-2nu) (1 - mu^2)(mu^2 + 1/x^2)^nu
        J = quad_pow0(lambda m: (1 - m**2) * (m**2 + 1 / x**2) ** nu, 1 - 2 * nu, [0, M('1e-6'), M('1e-3'), M('0.1'), 1])
        ratio_quad = (2 - s) * (4 - s) / 2 * x**s * J
        F = mp.hyp2f1(-nu, 1 - nu, 3 - nu, -x**2)
        Fp = (1 + x**2) ** nu * mp.hyp2f1(-nu, 2, 3 - nu, x**2 / (1 + x**2))
        worst = max(worst, rel(ratio_quad, F)); worst_pf = max(worst_pf, rel(Fp, F))
check('Zank: 2F1 form versus quadrature (nu = 0.7, 5/6, 0.9)', worst, M('1e-12'))
check('Zank: Pfaff-transformed 2F1 equals 2F1', worst_pf, M('1e-20'))

# 5. TS2003 exact integral, Eq. (45), against quadrature
def Jt(R, s):
    if R <= 1: return 2 * R ** (-s) / ((2 - s) * (4 - s))
    return M(1) / 4 + (1 / (2 - s) - M(1) / 2) / R**2 + (M(1) / 4 - 1 / (4 - s)) / R**4
worst = M(0)
for s in [M('1.5'), M('1.65'), M(5) / 3, M('1.69'), M('1.9')]:
    for R in ['0.1', '1', '2', '3.03', '10', '100']:
        R = M(R)
        # inertial part mu^(1-s)(1-mu^2) R^-s on [0, min(1, 1/R)], flat part mu(1-mu^2) on [1/R, 1]
        u = min(M(1), 1 / R)
        Jq = R ** (-s) * u ** (2 - s) * quad_pow0(lambda t: 1 - (u * t)**2, 1 - s)
        if R > 1: Jq += mp.quad(lambda m: m * (1 - m**2), [1 / R, 1])
        worst = max(worst, rel(Jq, Jt(R, s)))
check('TS2003: closed form versus quadrature (5 values of s)', worst, M('1e-12'))
s = M(5) / 3
form = lambda R: M(1) / 4 + 2 * R ** (-s) / ((2 - s) * (4 - s))
check('TS2003: exact/form = 72/79 at R = 1 (s = 5/3)', rel(Jt(M(1), s) / form(M(1)), M(72) / 79), M('1e-25'))
Rm = mp.findroot(lambda y: mp.diff(lambda z: Jt(z, s) / form(z), y), M(3))
check('TS2003: minimum exact/form = 0.793821 (6 digits)', abs(Jt(Rm, s) / form(Rm) - M('0.793821')), M('5e-7'))
check('TS2003: minimum at R = 3.0284 (5 digits)', abs(Rm - M('3.0284')), M('5e-5'))

# 5b. Full Eqs. (44) and (45) with their prefactors: lambda from Eq. (11) with the QLT
#     D_mumu = 2 pi^2 Omega^2 (1 - mu^2) g(k_res) / (B^2 v |mu|), k_res = 1/(R_L |mu|), and the normalized spectra
#     (8 pi int_0^inf g dk = dB^2), which gives lambda = (3 R_L^2 B^2 / (8 pi^2)) int_0^1 mu (1 - mu^2) / g dmu.
Bf, dB2, ell = M('5e-9'), M('4e-18'), M('2e9')
worst = M(0)
for nu in [M('0.7'), M(5) / 6]:
    Cn = mp.gamma(nu) / (2 * mp.sqrt(mp.pi) * mp.gamma(nu - M(1) / 2)); s = 2 * nu
    norm = 8 * mp.pi * mp.quad(lambda k: Cn * ell * dB2 * (1 + k**2 * ell**2) ** (-nu) / (2 * mp.pi), [0, 1 / ell, 10 / ell, mp.inf])
    worst = max(worst, rel(norm, dB2))
    for x in [M('0.01'), M('0.5'), M(5)]:
        RL = x * ell
        # mu (1 - mu^2)/g = (2 pi/(C l dB^2)) mu^(1-2nu) (1 - mu^2) (mu^2 + (l/RL)^2)^nu
        J = 2 * mp.pi / (Cn * ell * dB2) * quad_pow0(lambda m: (1 - m**2) * (m**2 + (ell / RL)**2) ** nu, 1 - 2 * nu, [0, M('1e-3'), M('0.1'), 1])
        lam_q = 3 * RL**2 * Bf**2 / (8 * mp.pi**2) * J
        lam_eq = 3 / (2 * mp.pi * Cn * (2 - s) * (4 - s)) * Bf**2 / dB2 * RL ** (2 - s) * ell ** (s - 1) * mp.hyp2f1(-nu, 1 - nu, 3 - nu, -x**2)
        worst = max(worst, rel(lam_q, lam_eq))
check('Zank: spectrum normalization and full Eq. (44) from Eq. (11)', worst, M('1e-12'))
worst = M(0)
kmin = M('1e-9')
for s in [M('1.5'), M(5) / 3]:
    g0 = dB2 * kmin ** (s - 1) * (s - 1) / (8 * mp.pi * s)
    gfun = lambda k: g0 * kmin ** (-s) if k <= kmin else g0 * k ** (-s)
    norm = 8 * mp.pi * (mp.quad(gfun, [0, kmin]) + mp.quad(gfun, [kmin, mp.inf]))
    worst = max(worst, rel(norm, dB2))
    for R in [M('0.1'), M(1), M(3), M(10)]:
        RL = R / kmin; u = min(M(1), 1 / R)
        Jin = (RL * u) ** (-s) * u**2 / g0 * quad_pow0(lambda t: 1 - (u * t)**2, 1 - s)   # mu < 1/R: 1/g = (R_L mu)^(-s)/g0
        Jfl = mp.quad(lambda m: m * (1 - m**2), [1 / R, 1]) * kmin**s / g0 if R > 1 else M(0)
        lam_q = 3 * RL**2 * Bf**2 / (8 * mp.pi**2) * (Jin + Jfl)
        lam_eq = 3 * s * RL**2 * kmin / (mp.pi * (s - 1)) * Bf**2 / dB2 * Jt(R, s)
        worst = max(worst, rel(lam_q, lam_eq))
check('TS2003: spectrum normalization and full Eq. (45) from Eq. (11)', worst, M('1e-12'))

# 6. He & Wan series, Eq. (43), and a float64 evaluation policy
def hw_exact(x):
    with mp.workdps(80):
        x = M(x); return 3 * (x - mp.tanh(x)) / x**3
def hw_series(x):
    x2 = x * x
    return 1 - 2 * x2 / 5 + 17 * x2**2 / 105 - 62 * x2**3 / 945 + 1382 * x2**4 / 51975
worst = M(0)
for x in ['1e-6', '0.001', '0.01', '0.03', '0.05']:
    worst = max(worst, rel(hw_series(M(x)), hw_exact(x)))
check('He-Wan: series through x^8 for |x| <= 0.05', worst, M('1.1e-15'))
def hw_float(x):
    return hw_series(x) if abs(x) <= 0.05 else 3.0 * (x - math.tanh(x)) / x**3
worst = M(0)
for k in range(0, 121):
    x = 10 ** (-10 + k * 11.3 / 120)          # 1e-10 .. 20
    worst = max(worst, rel(hw_float(x), hw_exact(repr(x))))
check('He-Wan: float64 series/direct policy on 1e-10..20', worst, M('1e-12'))
naive = 3.0 * (1e-6 - math.tanh(1e-6)) / 1e-6**3
print('INFO  He-Wan: the direct float64 formula at x = 1e-6 has relative error %s (cancellation; not counted)' % mp.nstr(rel(naive, hw_exact('1e-6')), 3))

# 7. Zank auxiliary functions in float64 with expm1/log1p
def aux_float(sv):
    s2 = sv * sv
    A = math.expm1(5.0 / 6.0 * math.log1p(s2))
    q = (5.0 / 3.0 * s2) / (s2 - math.expm1(math.log1p(s2) / 6.0))
    return A, q
def aux_mp(sv):
    with mp.workdps(60):
        s2 = M(sv)**2
        return (1 + s2) ** (M(5) / 6) - 1, (M(5) / 3 * s2) / (1 + s2 - (1 + s2) ** (M(1) / 6))
worst = M(0)
for sv in ['1e-8', '1e-4', '0.01', '1', '10']:
    Af, qf = aux_float(float(sv)); Am, qm = aux_mp(sv)
    worst = max(worst, rel(Af, Am), rel(qf, qm))
check('Zank auxiliary A, q in float64 (expm1/log1p) at s = 1e-8..10', worst, M('1e-12'))
s2 = M('1e-40')
check('Zank q(s) -> 2 as s -> 0 (stable form at s = 1e-20)', abs((M(5) / 3 * s2) / (s2 - mp.expm1(mp.log1p(s2) / 6)) - 2), M('1e-30'))

# 8. Ito SDE, Eq. (46): forward operator equals d_mu(D d_mu f)
def D(m): return M('0.7') * (1 - m**2) * (abs(m) ** (M(2) / 3) + M('0.05'))
def f(m): return mp.exp(m) * mp.cos(2 * m)
worst = M(0)
for m0 in ['-0.8', '-0.3', '0.2', '0.6', '0.9']:
    m0 = M(m0)
    L1 = mp.diff(lambda m: D(m) * mp.diff(f, m), m0)
    A = lambda m: mp.diff(D, m)
    L2 = -mp.diff(lambda m: A(m) * f(m), m0) + mp.diff(lambda m: D(m) * f(m), m0, 2)
    Lh = mp.diff(lambda m: D(m) / 2 * mp.diff(f, m), m0)
    Lh2 = -mp.diff(lambda m: A(m) / 2 * f(m), m0) + mp.diff(lambda m: D(m) / 2 * f(m), m0, 2)
    worst = max(worst, abs(L1 - L2), abs(Lh - Lh2))
check('Ito SDE forward operator = d_mu(D d_mu f), also with D/2', worst, M('1e-15'))

# 9. Moment decay and lambda for D = (nu0/2)(1 - mu^2): quadrature and a seeded simulation
nu0 = M(1)
lam_iso = 3 * v / 8 * quad_sym(lambda m: (1 - m**2)**2 / (nu0 / 2 * (1 - m**2)))
check('Isotropic D: lambda = v/nu0', rel(lam_iso, v / nu0), M('1e-20'))
rng = random.Random(20261009)
N, dt, T = 4000, 0.001, 0.3        # t = 0.3/nu0: the transient of <mu^2> (0.19) dominates the sampling error
mus = [0.9] * N
sq = math.sqrt(dt)
for _ in range(int(round(T / dt))):
    for j in range(N):
        m = mus[j]
        m += -m * dt + math.sqrt(max(0.0, 1.0 - m * m)) * sq * rng.gauss(0.0, 1.0)
        if m > 1.0: m = 2.0 - m
        elif m < -1.0: m = -2.0 - m
        mus[j] = m
m1 = sum(mus) / N; m2 = sum(x * x for x in mus) / N
se1 = math.sqrt(max(m2 - m1 * m1, 0.0) / N); se2 = math.sqrt(sum((x * x - m2)**2 for x in mus) / N / N)
th1 = 0.9 * math.exp(-T); th2 = 1 / 3 + (0.81 - 1 / 3) * math.exp(-3.0 * T)
bias = 5 * dt          # allowance for the O(dt) weak error of the Euler-Maruyama scheme with reflection
check('SDE <mu(t)> = mu0 exp(-nu0 t) (N=4000, dt=0.001, t=0.3)', abs(m1 - th1), 4 * se1 + bias, '(sim %.4f, theory %.4f)' % (m1, th1))
check('SDE <mu^2(t)> = 1/3 + (mu0^2 - 1/3) exp(-3 nu0 t), t=0.3', abs(m2 - th2), 4 * se2 + bias, '(sim %.4f, theory %.4f, equilibrium 1/3)' % (m2, th2))

# 10. Tensor projection, Eqs. (3) and (32), by rotating diag(K_par, K_perp)
worst = M(0)
for psi_deg in ['10', '44.7078', '80']:
    psi = mp.radians(M(psi_deg)); Kp, Kq = M(3), M('0.2')
    b = [mp.cos(psi), -mp.sin(psi)]                     # unit field vector in (r, phi), B_phi < 0
    K = [[Kp * b[i] * b[j] + Kq * ((1 if i == j else 0) - b[i] * b[j]) for j in range(2)] for i in range(2)]
    worst = max(worst, abs(K[0][0] - (Kp * mp.cos(psi)**2 + Kq * mp.sin(psi)**2)), abs(K[1][1] - (Kq * mp.cos(psi)**2 + Kp * mp.sin(psi)**2)),
                abs(K[0][1] - (Kq - Kp) * mp.cos(psi) * mp.sin(psi)))
check('Tensor rotation reproduces K_rr, K_phiphi, K_rphi', worst, M('1e-25'))

# 11. Parker field, Eq. (42): |B| / |B_r| = 1/|cos psi|; value at 1 AU in the ecliptic
AU = M(149597870700)
tpsi = M('2.66e-6') * (AU - M('0.005') * AU) / M(400000)
ratio = mp.sqrt(1 + tpsi**2)
check('Parker |B|/|B_r| = 1/cos psi at 1 AU (equals 1.4071 for 400 km/s)', max(abs(ratio - 1 / mp.cos(mp.atan(tpsi))), abs(ratio - M('1.4071')) / M('1.4071')), M('5e-5'), '(value %s)' % mp.nstr(ratio, 8))

# 12. Afanasiev et al. (2015): steady-state I_w gives lambda = 3(u1 - VA)(x + x0)/v
rngm = random.Random(7)
worst = M(0)
for _ in range(5):
    Om, B, u, VA, x, x0, vv = [M(rngm.uniform(0.5, 5)) for _ in range(7)]; u += VA
    k = Om / vv; Iw = Om * B**2 * abs(k) ** -3 / (3 * mp.pi * (u - VA) * (x + x0))
    nu = mp.pi * Om**2 * Iw / (vv * B**2)
    worst = max(worst, rel(vv / nu, 3 * (u - VA) * (x + x0) / vv))
check('Afanasiev: lambda = v/nu equals the printed lambda', worst, M('1e-25'))

# 13. Droge / PARADISE perpendicular identity (Section 10)
avg = mp.quad(lambda m: mp.sqrt(1 - m**2), [-1, 1]) / 2
check('<sqrt(1 - mu^2)> = pi/4', abs(avg - mp.pi / 4), M('1e-25'))

# 14. Printed (1 - mu)^2 variant (Section 6.2): D(-1) = 4 D0 (1 + H); normalization finite
D0, H = M('0.37'), M('0.05')
Dv = lambda m: D0 * (abs(m) ** (M(2) / 3) + H) * (1 - m)**2
check('Lang printed variant: D(-1) = 4 D0 (1 + H)', abs(Dv(M(-1)) - 4 * D0 * (1 + H)), M('1e-25'))
In = quad_sym(lambda m: (1 - m**2)**2 / Dv(m)); In2 = quad_sym(lambda m: (1 + m)**2 / (D0 * (abs(m) ** (M(2) / 3) + H)))
check('Lang printed variant: normalization integrand reduces to (1+mu)^2/(...)', rel(In, In2), M('1e-12'))

# 15. Corti <-> NWU normalization identity, Eq. (35), random parameters
worst = M(0)
for _ in range(5):
    a, b, s_, Rk = M(rngm.uniform(0.2, 1)), M(rngm.uniform(1.2, 2.2)), M(rngm.uniform(1, 4)), M(rngm.uniform(1, 6))
    k0 = Rk**b / (1 + Rk**s_) ** ((b - a) / s_)
    for R in [M('0.3'), M(2), M(50)]:
        nwu = R**a * ((R**s_ + Rk**s_) / (1 + Rk**s_)) ** ((b - a) / s_)
        cor = k0 * (R / Rk) ** a * (1 + (R / Rk) ** s_) ** ((b - a) / s_)
        worst = max(worst, rel(nwu, cor))
check('Corti <-> NWU identity for random a, b, s, R_k', worst, M('1e-25'))

# 16. Liu et al. (2025) Eq. (16) equals v lambda/3, random energies
c = M(299792458); E0 = M('1.67262192369e-27') * c**2 / (M('1.602176634e-19') * 10**6)
worst = M(0)
for T in [M('0.05'), M('3.7'), M(250), M(9000)]:
    p = mp.sqrt(T * (T + 2 * E0)); bt = p / (T + E0)
    D1 = bt * c * (p / 1000) ** (M(1) / 3) / 3
    D2 = c / 3 * (T * (T + 2 * E0) / M(10)**6) ** (M(1) / 6) * mp.sqrt(T * (T + 2 * E0) / (T + E0)**2)
    worst = max(worst, rel(D1, D2))
check('Liu Eq. (16) equals v lambda/3', worst, M('1e-25'))

# 17. Closed forms used by reference_verification.py: phi(eps), Droge ratio
def phi_q(eps): return M(3) / 2 * mp.quad(lambda m: (1 - m**2) * (1 + m) / (m + eps * (1 + m)), [0, eps, 10 * eps, 1])
def phi_c(eps):
    a = 1 + eps; m0 = -eps / a; N = [M(-1), M(-1), M(1), M(1)]; q = [N[0]]
    for coef in N[1:]: q.append(coef + q[-1] * m0)
    rem = q.pop(); deg = len(q) - 1
    return M(3) / 2 * (sum(qk / (deg - k + 1) for k, qk in enumerate(q)) / a + rem / a * mp.log((a + eps) / eps))
check('phi(eps) closed form versus quadrature', max(rel(phi_c(M(e)), phi_q(M(e))) for e in ['0.048', '0.01', '0.001']), M('1e-15'))
def dr_q(a): return mp.quad(lambda m: (1 - m**2) / (m**2 + a**2) ** (M(1) / 3), [0, a, 1])
def dr_c(a):
    z = -1 / a**2; p = M(1) / 3
    return a ** (-2 * p) * (mp.hyp2f1(p, M(1) / 2, M(3) / 2, z) - mp.hyp2f1(p, M(3) / 2, M(5) / 2, z) / 3)
check('Droge integral hypergeometric form versus quadrature', max(rel(dr_c(M(a)), dr_q(M(a))) for a in ['0.001', '0.01', '0.1', '0.3']), M('1e-15'))

n_pass = sum(results)
print('\nindependent physics checks passed: %d of %d' % (n_pass, len(results)))
sys.exit(0 if n_pass == len(results) else 1)

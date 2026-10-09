#!/usr/bin/env python3
"""Independent equation calculations for perpendicular specification 2.1.

This is a verification artifact, not AMPS code or a physical validation.
All unit-scale states and arbitrary fit inputs below are synthetic tests.
Run `python3 reference_verification.py` to check the delivered fixtures.
Only a deliberate specification update should use --generate-fixtures.
"""

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erfc, erfcx, gamma, hyp2f1

import paired_reference_verification as paired

ROOT = Path(__file__).resolve().parent
S = 5.0 / 3.0
C = float(gamma(S / 2) / (2 * math.sqrt(math.pi) * gamma((S - 1) / 2)))
CHECKS = []


def check(actual, expected, label, rtol=2e-9, atol=1e-13):
    if not (math.isfinite(actual) and math.isfinite(expected)):
        raise AssertionError(f"Nonfinite {label}")
    error = abs(actual - expected)
    budget = atol + rtol * abs(expected)
    if error > budget:
        raise AssertionError(f"{label}: {actual:.17g} versus {expected:.17g}; budget {budget:g}")
    CHECKS.append({"check": label, "absolute_error": error,
                   "relative_error": error / abs(expected) if expected else None,
                   "relative_tolerance": rtol, "absolute_tolerance": atol})


def dnorm(q, s=S):
    return float(gamma((s + q) / 2) / (2 * gamma((s - 1) / 2) * gamma((q + 1) / 2)))


def log_integral(function, bound=48.0, tol=2e-12):
    edges = [-bound, -15.0, -5.0, 0.0, 5.0, 15.0, bound]
    return math.fsum(quad(function, a, b, epsabs=tol, epsrel=tol, limit=300)[0]
                     for a, b in zip(edges[:-1], edges[1:]))


def log_shape(y, q, s=S):
    return q * y - (s + q) / 2 * float(np.logaddexp(0.0, 2 * y))


def spectral_moment(q, moment, variance=1.0, ell=1.0, s=S):
    if q + moment <= -1 or s <= moment + 1:
        raise ValueError("The uncut spectrum has a divergent moment")
    return (dnorm(q, s) / math.pi * variance * ell ** (-moment)
            * gamma((q + moment + 1) / 2) * gamma((s - moment - 1) / 2)
            / gamma((s + q) / 2))


def lengths(q):
    lu = math.sqrt((S - 1) / (q - 1))
    lp = 2 * gamma(q / 2) * gamma(S / 2) / (gamma((q + 1) / 2) * gamma((S - 1) / 2))
    return {"q": q, "L_U_over_ell_2": lu, "L_perp_over_ell_2": float(lp),
            "CLRR_shape_prefactor": float(lp * lp / (4 * lu * lu))}


def fieldlines(q, fs, variance):
    ks = math.pi * C * fs * variance
    k2 = math.sqrt((S - 1) / (2 * (q - 1)) * (1 - fs) * variance)
    kf = (ks + math.sqrt(ks * ks + 4 * k2 * k2)) / 2
    return ks, k2, kf


def flpd(q, fs, variance, lp, bound=42.0, tol=1e-11):
    """Dimensionless D2; unknown is ln(kappa_perp/kappa_parallel)."""
    _, _, kf = fieldlines(q, fs, variance)
    bx2 = (1 - fs) * variance / 2

    def residual(log_eta):
        log_alpha = math.log(kf) - log_eta / 2
        log_a = 2 * math.log(lp) + log_eta - math.log(3)

        def integrand(y):
            log_u = float(np.logaddexp(0.0, log_a + 2 * y))
            log_den = float(np.logaddexp(log_u, log_alpha + y + log_u / 2))
            return math.exp(log_shape(y, q) + y - log_den)

        integral = log_integral(integrand, bound, tol)
        return math.log(4 * dnorm(q) * bx2 * integral) - log_eta

    result = brentq(residual, -45.0, math.log(bx2), xtol=2e-13)
    check(residual(result), 0.0, f"FLPD log residual q={q},fs={fs},var={variance},lp={lp}", atol=1e-10)
    return math.exp(result)


def dimensional_flpd(q, fs, variance, lp):
    """Direct D1, with nonunit v, B0, ell; ordinary k quadrature."""
    v, b0, ell = 2.3, 7.0, 4.1  # Chosen unit-restoration test inputs.
    kp = v * lp * ell / 3
    rate = v * v / (3 * kp)
    _, _, kf_unit = fieldlines(q, fs, variance)
    kf = kf_unit * ell
    v2 = (1 - fs) * variance * b0 * b0

    def rhs(kt):
        def integrand(k):
            x = k * ell
            g = 2 * dnorm(q) / math.pi * v2 * ell * x**q / (1 + x*x)**((S + q) / 2)
            u = rate + kt * k*k
            return g / (u + v * kf / math.sqrt(3 * kt) * k * math.sqrt(u))
        return math.pi * v*v / (3 * b0*b0) * quad(
            integrand, 0, np.inf, epsabs=1e-12, epsrel=2e-12, limit=400)[0]

    def residual(log_kt):
        kt = math.exp(log_kt)
        return math.log(rhs(kt)) - log_kt

    upper = kp * v2 / (2 * b0 * b0)
    kt = math.exp(brentq(residual, math.log(upper) - 35, math.log(upper), xtol=2e-13))
    return kt / kp


def enlgc(lp, bound=48.0, tol=2e-12):
    def residual(log_lt):
        log_a = math.log(lp) + log_lt - math.log(3)
        integral = log_integral(lambda y: math.exp(
            log_shape(y, 0) + y - float(np.logaddexp(0.0, log_a + 2*y))), bound, tol)
        return math.log(2 * C * 0.8 * lp * integral) - log_lt
    result = brentq(residual, -45, math.log(0.4 * lp), xtol=2e-13)
    check(residual(result), 0.0, f"ENLGC log residual lp={lp}", atol=1e-10)
    return math.exp(result)


def kernel_positive(xi):
    if xi < 20:
        return quad(lambda u: 2*u*math.exp(-u*u-2*xi*u), 0, np.inf,
                    epsabs=1e-14, epsrel=2e-12)[0]
    return quad(lambda t: t*math.exp(-t-t*t/(4*xi*xi)), 0, np.inf,
                epsabs=1e-14, epsrel=2e-12)[0] / (2*xi*xi)


def rbd(variance, lp, nu=S, p=2.0):
    a2 = 1 / 3
    vx = a2 * variance / 6
    c2 = (nu-1)*(p+2) / (2*math.pi*(p+nu+1))
    qrate = 0 if math.isinf(lp) else 1 / lp

    def integrand(x, exponent):
        if x == 0:
            return 0.0
        return x**exponent * erfc(qrate/(x*math.sqrt(2*vx)))
    integral = (quad(lambda x: integrand(x,p), 0,1, epsabs=1e-13, epsrel=2e-12)[0]
                + quad(lambda x: integrand(x,-nu-1),1,np.inf,epsabs=1e-13,epsrel=2e-12)[0])
    return a2 / 6 * math.sqrt(math.pi/2) * 2*math.pi*c2*variance/math.sqrt(vx)*integral


def closed(lp, d=0.1):
    return 4*d*d*lp / (math.sqrt(1+8*d*lp/3)+1)**2


def calculate_fixtures():
    cases = [(3,0,1),(2,0.2,1),(3,0.2,1),(3,0.5,1),(3,0.2,0.25),(3,0.2,2)]
    lps = [0.01,0.1,1,10,100,1000]
    rows = []
    for i,(q,fs,var) in enumerate(cases,1):
        ks,k2,kf = fieldlines(q,fs,var)
        values = [flpd(q,fs,var,lp) for lp in lps]
        refined = [flpd(q,fs,var,lp,bound=48,tol=2e-12) for lp in lps]
        for lp,a,b in zip(lps,values,refined):
            check(a,b,f"FLPD refinement state={i},lp={lp}",rtol=2e-9)
        rows.append({"state":i,"q":q,"slab_fraction":fs,"total_variance_over_B0_squared":var,
                     "kappa_s_over_ell":ks,"kappa_2_over_ell":k2,"kappa_FL_over_ell":kf,
                     "parallel_lengths":lps,"ratio":refined,
                     "maximum_relative_refinement_change":max(abs(a/b-1) for a,b in zip(values,refined)),
                     "long_limit_lambda_perp_over_ell":math.sqrt(3)/2*(math.sqrt(kf*kf+4*k2*k2)-kf)})
    for q,fs,var,lp in [(3,0.2,0.25,1),(2,0.2,1,10),(3,0,1,100)]:
        check(dimensional_flpd(q,fs,var,lp),flpd(q,fs,var,lp),
              f"FLPD dimensional versus dimensionless q={q},fs={fs},lp={lp}",rtol=2e-9)

    kr = []
    for xi in [0.01,0.1,0.5,1,2,5,10]:
        exact = float(1-math.sqrt(math.pi)*xi*erfcx(xi))
        positive = kernel_positive(xi)
        check(exact,positive,f"Implicit kernel positive integral xi={xi}",rtol=3e-11)
        kr.append({"xi":xi,"K_exact":positive,
                   "rational_relative_error":1/(1+2*xi*xi)/positive-1})
    check(kernel_positive(1e4),1/(2*1e8)-3/(4*1e16),
          "Implicit large-xi asymptotic",rtol=1e-12,atol=1e-20)
    rbd_rows = [{"variance":var,"lambda_parallel":lp,"kappa_perp_over_v_lambda_2":rbd(var,lp)}
                for var,lp in [(0.1,100),(0.1,1000),(0.4,100),(0.4,1000),(0.8,100)]]
    c2 = (S-1)*4 / (2*math.pi*(2+S+1))
    limit = math.sqrt(math.pi)/6 * 2*math.pi*c2 * (1/3+1/S)
    for var in [0.1,0.4,0.8,1.0]:
        check(rbd(var,math.inf),limit*math.sqrt(var),f"RBD analytic long limit var={var}",rtol=2e-11)
    closed_rows = [{"lambda_parallel_over_ell":lp,"lambda_perp_over_ell":closed(lp)}
                   for lp in [0.1,1,10,100,10000]]
    for row in closed_rows:
        lp,lt = row["lambda_parallel_over_ell"],row["lambda_perp_over_ell"]
        check(2*lt/3+math.sqrt(lt/lp),0.1,f"Composite B5 substitution lp={lp}",rtol=1e-12)
    alpha = float(gamma(7/6)/math.sqrt(math.pi)*(18*C*math.sqrt(math.pi/2))**(2/3))
    fld_msd_bracket = (9*C*math.sqrt(math.pi/2))**(2/3)
    normal_moment = 2 / math.sqrt(2*math.pi) * quad(
        lambda z:z**(4/3)*math.exp(-z*z/2),0,np.inf,epsabs=1e-13,epsrel=2e-12)[0]
    check(fld_msd_bracket*normal_moment,alpha,"GCD alpha direct Gaussian moment",rtol=1e-11)
    au,omega,wind = 149597870700.0,2.87e-6,400000.0
    psi = math.atan(omega*au/wind)
    return {"schema_version":"2.0","normalization":"Dimensionless states; B0=v=ell=1 unless explicitly stated.",
            "numerical_controls":{
                "classification":"Chosen engineering controls, not physical parameters or cutoffs.",
                "FLPD_log_k_bounds_initial":[-42.0,42.0],"FLPD_initial_quadrature_tolerance":1e-11,
                "FLPD_log_k_bounds_refined":[-48.0,48.0],"FLPD_refined_quadrature_tolerance":2e-12,
                "FLPD_root_log_absolute_tolerance":2e-13,"FLPD_log_residual_acceptance":1e-10,
                "dimensional_FLPD_state":{"v":2.3,"B0":7.0,"ell":4.1},
                "ENLGC_log_k_bounds":[-48.0,48.0],"ENLGC_quadrature_tolerance":2e-12,
                "NLGCE_log_k_bounds_initial":[-40.0,40.0],"NLGCE_initial_quadrature_relative_tolerance":1e-10,
                "NLGCE_log_k_bounds_refined":[-45.0,45.0],"NLGCE_refined_quadrature_relative_tolerance":3e-11,
                "NLGCE_other_seed_scale":3.0,
                "NLGCE_nonlinear_fixture_relative_tolerance":5e-8,
                "NLGCE_polynomial_fixture_relative_tolerance":1e-11},
            "physical_validation_claimed":False,
            "constants":{"s":S,"C_s":C,"pi_C_s":math.pi*C,"L_s_over_ell_s":2*math.pi*C},
            "lengths":[lengths(q) for q in [1.5,2,3]],"flpd_states":rows,
            "flpd_limits":{"pure2d_long_factor":(math.sqrt(5)-1)/math.sqrt(3),
                           "pure2d_ratio_at_parallel_100000":flpd(3,0,1,1e5,bound=48,tol=2e-12),
                           "pure2d_CLRR_ratio":lengths(3)["CLRR_shape_prefactor"]/2},
            "enlgc":[{"lambda_parallel":lp,"lambda_perp":enlgc(lp)} for lp in [0.01,0.1,1,10,1000]],
            "rbd":rbd_rows,"rbd_long_limit_coefficient":limit,"closed_form":closed_rows,
            "implicit_slab_kernel":kr,"gcd":{"alpha_s_5_over_3":alpha,
                                             "fieldline_msd_shape_coefficient":fld_msd_bracket,
                                             "secant_average_lambda_prefactor_variance_0_8":(3/2)**(4/3)*alpha*0.8**(2/3)},
            "chosen_Parker_example":{"rotation_rate_s_inverse":omega,"wind_m_s":wind,"r_m":au,
                                     "psi_degrees":math.degrees(psi),"cos_squared":math.cos(psi)**2,
                                     "cos":math.cos(psi),"alpha_D":0.3,
                                     "averaged_ratio":0.3*math.cos(psi)*math.pi/4}}


def verify_spectral_equations():
    slab = log_integral(lambda y:math.exp(y-S/2*float(np.logaddexp(0,2*y))))
    check(2*C*slab,0.5,"Slab per-axis S7 variance",rtol=3e-11)
    check(4*C*slab,1.0,"NLGCE E1 total variance",rtol=3e-11)
    for q in [0,1.5,2,3]:
        integral = log_integral(lambda y:math.exp(log_shape(y,q)+y))
        check(2*dnorm(q)*integral,0.5,f"2D S7 per-axis q={q}",rtol=3e-11)
        # S2=g2/k independently inserted into its total-variance integral.
        check(2*math.pi*(2*dnorm(q)/math.pi)*integral,1.0,
              f"Area adapter S11 total variance q={q}",rtol=3e-11)
        for moment in [-2,-1,0]:
            if q+moment <= -1:
                try:
                    spectral_moment(q,moment)
                except ValueError:
                    continue
                raise AssertionError("Divergent moment was accepted")
            numerical = 2*dnorm(q)/math.pi*log_integral(
                lambda y:math.exp(log_shape(y,q)+(moment+1)*y),bound=70)
            check(numerical,spectral_moment(q,moment),f"Moment S8 q={q},m={moment}",rtol=3e-10)
        for a in [0,0.01,1,100]:
            integral_a = log_integral(lambda y:math.exp(log_shape(y,q)+y)/
                                     (1+a*math.exp(2*y)))
            direct = 4*dnorm(q)*integral_a
            special = (S-1)/(S+q)*hyp2f1(1,(q+1)/2,(S+q)/2+1,1-a)
            check(float(special),direct,f"Hypergeometric N3 q={q},A={a}",rtol=2e-9)
    try:
        spectral_moment(3,1)
    except ValueError:
        pass
    else:
        raise AssertionError("Divergent I1 accepted")
    nu,p = S,2
    c2 = (nu-1)*(p+2)/(2*math.pi*(p+nu+1))
    area_norm = 2*math.pi*c2*(1/(p+2)+1/(nu-1))
    check(area_norm,1.0,"Chhiber area S12 total variance",rtol=1e-12)
    outer,ell,q = 8.0,1.3,2.0
    c0 = 1/(1/(q+1)+math.log(outer/ell)+1/(nu-1))
    check(c0*(1/(q+1)+math.log(outer/ell)+1/(nu-1)),1.0,"Supplied S13 normalization",rtol=1e-12)


def verify_reduced_kernels():
    # Nonunit test state, with distinct bend-over lengths.
    v,b0,ls,l2,kp,kt,vs,v2,a2,q = 1.7,2.2,0.8,1.3,2.1,0.17,0.3,0.8,1/3,3.0
    rate = v*v/(3*kp)
    gs = lambda k:C/(2*math.pi)*vs*ls/(1+(k*ls)**2)**(S/2)
    g2 = lambda k:2*dnorm(q)/math.pi*v2*l2*(k*l2)**q/(1+(k*l2)**2)**((S+q)/2)
    direct = a2*v*v/(3*b0*b0)*(4*math.pi*quad(lambda k:gs(k)/(rate+kp*k*k),0,np.inf)[0]
             +math.pi*quad(lambda k:g2(k)/(rate+kt*k*k),0,np.inf)[0])
    h = lambda q,a:(S-1)/(S+q)*hyp2f1(1,(q+1)/2,(S+q)/2+1,1-a)
    hyper = kp*a2/(2*b0*b0)*(vs*h(0,kp/(rate*ls*ls))+v2*h(q,kt/(rate*l2*l2)))
    check(direct,float(hyper),"NLGC N2 versus N4 RHS",rtol=2e-9)
    unlt = a2*v*v*math.pi/(3*b0*b0)*quad(lambda k:g2(k)/(rate+4*kt*k*k/3),0,np.inf)[0]
    hyper_unlt = kp*a2*v2/(2*b0*b0)*h(q,4*kt/(3*rate*l2*l2))
    check(unlt,float(hyper_unlt),"UNLT U2 versus N3 RHS",rtol=2e-9)
    # NRMHD normalization and the independent k_parallel integration.
    ell,kc,variance = 0.9,0.7,0.6
    a0 = 8/9*ell**4*variance
    ak = lambda k:a0/(1+(k*ell)**2)**(7/3)
    check(0.5*quad(lambda k:k**3*ak(k),0,np.inf,epsabs=1e-12)[0],
          variance/2,"NRMHD U3 per-axis variance",rtol=2e-9)
    for label in ["NLGC","UNLT"]:
        def uv(k):
            return ((rate+kt*k*k,kp) if label=="NLGC" else
                    (rate+4*kt*k*k/3,v*v/(3*kt*k*k)))
        def analytic_k(k):
            u,w = uv(k)
            return k**3*ak(k)/(2*kc)*math.atan(kc*math.sqrt(w/u))/math.sqrt(u*w)
        def nested_k(k):
            u,w = uv(k)
            return k**3*ak(k)/(4*kc)*quad(lambda z:1/(u+w*z*z),-kc,kc,epsabs=2e-11)[0]
        analytic = quad(analytic_k,0,np.inf,epsabs=2e-11,epsrel=2e-10,limit=300)[0]
        nested = quad(nested_k,0,np.inf,epsabs=2e-11,epsrel=2e-10,limit=300)[0]
        check(analytic,nested,f"NRMHD U4 analytic versus nested integral {label}",rtol=2e-8)
    # Rational implicit-slab U9 follows algebraically from U5/U6/U8.
    slab_fl,lp,lt,l2,q = 0.13,4.0,0.25,1.1,3.0
    aa = lp*lt/(3*l2*l2)
    gg = 2*lp*lp*slab_fl*slab_fl/(3*math.pi*l2**4)
    for x in [0.01,0.3,1,7,100]:
        xi = slab_fl*lp*(x/l2)**2/(math.sqrt(3*math.pi)*math.sqrt(1+aa*x*x))
        check(1/((1+aa*x*x)*(1+2*xi*xi)),1/(1+aa*x*x+gg*x**4),
              f"Rational implicit-slab reduction x={x}",rtol=2e-12)


def verify_running_and_classical():
    kp,v = 2.7,1.3
    tau = 3*kp/(v*v)
    mz = lambda t:2*kp*t*t/(math.sqrt(t*t+tau*tau)+tau)
    dp = lambda t:kp*t/math.sqrt(t*t+tau*tau)
    for factor in [1e-4,0.1,1,10,1000]:
        t = factor*tau
        numerical = 2*quad(dp,0,t,epsabs=1e-12,epsrel=2e-12)[0]
        check(mz(t),numerical,f"Kuhlen I10 versus I9 integral t/tau={factor}",rtol=2e-10)
    # Chosen artificial parameters: these are not a source perpendicular preset.
    ck,gk,z1,z2,bx = 0.8,-0.3,3.0,16.0,0.7
    dfl = lambda z:ck*z*bx*(1+(z/z1)**((1-gk)/1.5))**(-1.5)*(1+(z/z2)**(-gk/0.2))**0.2
    for factor in [0.1,1,10]:
        t = factor*tau
        h = t*1e-5
        mx = lambda tt:2*quad(dfl,0,math.sqrt(mz(tt)),epsabs=1e-12,epsrel=2e-12)[0]
        running = dfl(math.sqrt(mz(t)))/math.sqrt(mz(t))*dp(t)
        derivative = (mx(t+h)-mx(t-h))/(4*h)
        check(derivative,running,f"Kuhlen I12 chain rule t/tau={factor}",rtol=2e-8)
    alpha = float(gamma(7/6)/math.sqrt(math.pi)*(18*C*math.sqrt(math.pi/2))**(2/3))
    ell,var,lp = 0.9,0.8,2.5
    tc = lp/v
    ak = alpha*var**(2/3)*ell**(2/3)*(2*v*lp/3)**(2/3)
    average = 3/v/tc*quad(lambda t:ak/2*t**(-1/3),0,tc,epsabs=1e-12,epsrel=2e-12)[0]
    closed_avg = (3/2)**(4/3)*alpha*var**(2/3)*ell**(2/3)*lp**(1/3)
    check(average,closed_avg,"GCD C6 secant scattering-time average",rtol=2e-10)
    a1,tl,t1,t2,am,bm = 0.23,2.0,20.0,60.0,1.5,3.0
    median = lambda t:a1*(t/tl)**am*(1+(t/t1)**(bm-am))/(1+(t/t2)**(bm-1))
    a2 = a1*t2**(bm-1)/(tl**(am-1)*t1**(bm-am))
    check(median(1e10),a2*1e10/tl,"Conditional median C10/C11 late asymptotic",rtol=2e-12,atol=0)
    check(median(1e-8),a1*(1e-8/tl)**am,"Conditional median C10 early asymptotic",rtol=2e-12,atol=0)
    for x in [0.01,0.2,1,5,100]:
        kpar = x/3  # v=rL=Omega=1; tau=x.
        kper = kpar/(1+x*x)
        ka = (1/3)*x*x/(1+x*x)
        check(kper*kper+ka*ka,kpar*kper,f"Classical Hall identity x={x}",rtol=1e-12)
        check(ka/kper,x,f"Classical Hall/perpendicular ratio x={x}",rtol=1e-12)


def verify_tensors():
    x = np.array([0.4,-0.7,0.2])
    jb = np.array([[0.3,0,0],[0,-0.3,0],[0.1,0,0]])

    def field(position):
        raw = np.array([1+0.3*position[0],-0.3*position[1],0.7+0.1*position[0]])
        return raw/np.linalg.norm(raw),np.linalg.norm(raw)

    def scalars(position):
        xx,yy,zz = position
        return 0.7+0.2*xx*xx+0.3*yy*yy,0.12+0.04*xx+0.02*zz*zz

    def tensor(position):
        b,_ = field(position)
        kp,kt = scalars(position)
        return kt*np.eye(3)+(kp-kt)*np.outer(b,b)

    b,norm = field(x)
    kp,kt = scalars(x)
    db = (np.eye(3)-np.outer(b,b))@jb/norm
    gpar = np.array([0.4*x[0],0.6*x[1],0])
    gper = np.array([0.04,0,0.04*x[2]])
    exact = gper+b*np.dot(b,gpar-gper)+(kp-kt)*(db@b+b*np.trace(db))
    errors = []
    for h in [1e-3,1e-4,1e-5]:
        finite = np.zeros(3)
        for j in range(3):
            offset = np.eye(3)[j]*h
            finite += (tensor(x+offset)[:,j]-tensor(x-offset)[:,j])/(2*h)
        errors.append(float(np.linalg.norm(finite-exact)))
        for i in range(3):
            check(float(finite[i]),float(exact[i]),f"Tensor divergence SDE5 versus SDE11 h={h},i={i}",
                  rtol=3e-6 if h==1e-3 else 5e-8,atol=2e-10)
    if errors[1] >= errors[0]/30:
        raise AssertionError("Tensor finite-difference refinement did not reduce the error")
    e1 = np.array([0,0,1.0])-b*b[2]
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(b,e1)
    basis = np.column_stack([b,e1,e2])
    for eig in ([kp,kt,kt],[kp,0.13,0.27],[kp,0,0]):
        k = basis@np.diag(eig)@basis.T
        noise = basis@np.diag(np.sqrt(2*np.asarray(eig)))
        for i in range(3):
            for j in range(3):
                check(float(k[i,j]),float(k[j,i]),f"Tensor symmetry eig={eig},i={i},j={j}",rtol=1e-12)
                check(float((noise@noise.T)[i,j]),float(2*k[i,j]),
                      f"SDE8 covariance eig={eig},i={i},j={j}",rtol=1e-12)
        for a,z in zip(np.linalg.eigvalsh(k),sorted(eig)):
            check(float(a),float(z),f"Tensor eigenvalue eig={eig},target={z}",rtol=1e-12)
    for form in [lambda mu:1.0,lambda mu:2*abs(mu),lambda mu:4/math.pi*math.sqrt(1-mu*mu)]:
        check(quad(form,-1,1,epsabs=1e-12)[0]/2,1.0,"P3 isotropic normalization",rtol=2e-10)


def compare_tree(actual,expected,path="fixtures"):
    if isinstance(expected,dict):
        if set(actual) != set(expected):
            raise AssertionError(f"Different keys at {path}")
        for key in expected:
            compare_tree(actual[key],expected[key],f"{path}.{key}")
    elif isinstance(expected,list):
        if len(actual) != len(expected):
            raise AssertionError(f"Different array length at {path}")
        for i,(a,b) in enumerate(zip(actual,expected)):
            compare_tree(a,b,f"{path}[{i}]")
    elif isinstance(expected,(int,float)) and not isinstance(expected,bool):
        if "maximum_relative_refinement_change" not in path:
            check(float(actual),float(expected),path,rtol=3e-9,atol=2e-13)
    elif actual != expected:
        raise AssertionError(f"Different metadata at {path}")


def verify_checksums():
    manifest = ROOT/"SHA256SUMS"
    if not manifest.exists():
        raise FileNotFoundError("Missing SHA256SUMS; use --skip-checksums only during an intentional asset revision")
    count = 0
    seen = set()
    for line in manifest.read_text().splitlines():
        digest,name = line.split("  ",1)
        if Path(name).name != name or name in seen:
            raise ValueError("Invalid checksum manifest path")
        actual = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        if actual != digest:
            raise AssertionError(f"Checksum mismatch: {name}")
        seen.add(name)
        count += 1
    delivered = {p.name for p in ROOT.iterdir() if p.is_file() and p.name!="SHA256SUMS"}
    if seen!=delivered:
        raise AssertionError("Checksum manifest does not cover every delivered regular file")
    return {"state":"verified","files":count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate-fixtures",action="store_true")
    parser.add_argument("--skip-checksums",action="store_true")
    parser.add_argument("--write-results",action="store_true",help="Write the current check report after a successful verification")
    args = parser.parse_args()
    if not args.skip_checksums:
        verify_checksums()
    values = calculate_fixtures()
    if args.generate_fixtures:
        (ROOT/"benchmark_points.json").write_text(json.dumps(values,indent=2,allow_nan=False)+"\n")
    else:
        compare_tree(values,json.loads((ROOT/"benchmark_points.json").read_text()))
    verify_spectral_equations()
    verify_reduced_kernels()
    verify_running_and_classical()
    verify_tensors()
    pair = paired.verify()
    result = {"state":"passed","scope":"Equation, numeric fixture and asset checks; no AMPS or particle-simulation validation.",
              "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
              "equation_comparisons":len(CHECKS),"checks":CHECKS,"paired_NLGCE":pair,
              "physical_validation_claimed":False}
    if args.generate_fixtures or args.write_results:
        (ROOT/"verification_results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:val for k,val in result.items() if k!="checks"},indent=2),flush=True)


if __name__ == "__main__":
    main()

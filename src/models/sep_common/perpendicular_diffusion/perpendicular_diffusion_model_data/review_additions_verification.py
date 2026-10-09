#!/usr/bin/env python3
"""Recompute revision 2.1 review additions; synthetic equation checks only.

Requires NumPy, SciPy and mpmath. Default execution is read-only.
The original reference script separately checks the manifest and old fixtures.
"""
import argparse
import json
import math
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import erfc, erfcx

import reference_verification as old

ROOT = Path(__file__).resolve().parent
CHECKS = []


def check(actual, expected, label, rtol=3e-10, atol=1e-13):
    actual, expected = float(actual), float(expected)
    error = abs(actual-expected)
    assert math.isfinite(actual) and math.isfinite(expected), label
    assert error <= atol+rtol*abs(expected), (label,actual,expected)
    CHECKS.append({"check":label,"absolute_error":error,
                   "relative_error":error/abs(expected) if expected else None,
                   "relative_tolerance":rtol,"absolute_tolerance":atol})


def require(condition, label):
    assert condition, label
    CHECKS.append({"check":label,"state":"passed"})


def calculate():
    mp.mp.dps = 50
    s = mp.mpf(5)/3
    c = mp.gamma(s/2)/(2*mp.sqrt(mp.pi)*mp.gamma((s-1)/2))
    alpha_g = mp.gamma(mp.mpf(7)/6)/mp.sqrt(mp.pi)*(18*c*mp.sqrt(mp.pi/2))**(mp.mpf(2)/3)
    out = {
        "schema_version":"1.0","specification_version":"2.1",
        "scope":"Algebraic review additions; synthetic inputs, no physical calibration or AMPS tests.",
        "precision_decimal_digits":50,
        "constants":{
            "C_s":float(c),"pi_C_s":float(mp.pi*c),
            "L_s_over_ell_s":float(2*mp.pi*c),
            "GCD_alpha":float(alpha_g),
            "GCD_secant_average_shape":float((mp.mpf(3)/2)**(mp.mpf(4)/3)*alpha_g),
            "RBD_long_shape":float(4*mp.sqrt(mp.pi)/45)},
        "short_flpd":[],"nonanalytic_tail":[],"rbd_plain_vs_bc":[],
        "closed_bounds":[]}
    for q in map(mp.mpf,["1.5","2","3"]):
        d = mp.gamma((s+q)/2)/(2*mp.gamma((s-1)/2)*mp.gamma((q+1)/2))
        h = lambda x:x**q/(1+x*x)**((s+q)/2)
        lu = mp.sqrt((s-1)/(q-1))
        lper = 2*mp.gamma(q/2)*mp.gamma(s/2)/(mp.gamma((q+1)/2)*mp.gamma((s-1)/2))
        integral = lambda z:mp.quad(lambda x:h(x)/(1+lu*x/mp.sqrt(z)),[0,1,mp.inf])
        zeta = mp.findroot(lambda z:z-4*d*integral(z),(.1,.2))
        alpha_root = lu/mp.sqrt(zeta)
        check(zeta,4*d*integral(zeta),f"D10 root q={q}",rtol=1e-14)
        i0 = mp.quad(h,[0,1,mp.inf])
        im1 = mp.quad(lambda x:h(x)/x,[0,1,mp.inf])
        check(im1/i0,lper/2,f"B8 integral length q={q}",rtol=1e-14)
        alpha_d4 = 2*lu**2/lper
        gamma_form = mp.gamma((q-1)/2)*mp.gamma((s+1)/2)/(mp.gamma(q/2)*mp.gamma(s/2))
        check(alpha_d4,gamma_form,f"D9 gamma identity q={q}",rtol=1e-14)
        # Derivative of the positive implicit map is < 1/2 at its root.
        derivative = 4*d*mp.quad(
            lambda x:h(x)*(alpha_root*x)/(2*zeta*(1+alpha_root*x)**2),[0,1,mp.inf])
        require(0<derivative<mp.mpf(1)/2,f"D3 fixed-point derivative q={q}")
        # Cross-check the independent double-precision short-limit solve.
        fp = old.flpd(float(q),0,1,1e-6)
        check(fp,float(zeta/2),f"D2 tends to D10 q={q}",rtol=1e-8)
        out["short_flpd"].append({
            "q":float(q),"alpha_fluid":float(lu),"alpha_D4":float(alpha_d4),
            "zeta":float(zeta),"alpha_root":float(alpha_root),
            "integral_length_over_ell_2":float(lper/2),
            "fixed_point_derivative":float(derivative)})
    q = mp.mpf(3)
    alpha = alpha_root  # Preserve the high-precision q=3 root in the integral.
    h = lambda x:x**q/(1+x*x)**((s+q)/2)
    tail = mp.quad(lambda y:y**(-s-1)*(y*y/(mp.sqrt(1+y*y)*(mp.sqrt(1+y*y)+1))),
                   [0,1,mp.inf])
    tail_beta = mp.gamma(1-s/2)*mp.gamma((1+s)/2)/(s*mp.sqrt(mp.pi))
    check(tail,tail_beta,"D13 convergent tail integral",rtol=1e-12)
    previous = None
    target = tail/alpha
    for aval in ["1e-8","1e-10","1e-12"]:
        a = mp.mpf(aval)
        def delta(x):
            u=a*x*x
            root=mp.sqrt(1+u)
            difference=u+alpha*x*u/(root+1)
            return h(x)*difference/((1+alpha*x)*(1+u+alpha*x*root))
        delta_j=mp.quad(delta,[0,1,1/mp.sqrt(a),mp.inf])
        scaled=delta_j/a**(s/2)
        require(0<scaled<target,f"D13 scaled difference A={aval}")
        exponent=None
        if previous:
            exponent=mp.log(delta_j/previous[1])/mp.log(a/previous[0])
            require(previous[2]<scaled,f"D13 approach to tail limit A={aval}")
            require(0<exponent<s/2,f"D13 local exponent A={aval}")
        out["nonanalytic_tail"].append({"A":float(a),"fixed_alpha":float(alpha),
            "J_zero_minus_J_A":float(delta_j),"scaled_difference":float(scaled),
            "asymptotic_scaled_limit":float(target),
            "local_log_exponent":float(exponent) if exponent else None})
        previous=(a,delta_j,scaled)
    for Q,V in [(0.2,0.4),(1.3,0.7),(0.5,0.02)]:
        beta=Q/math.sqrt(2*V)
        plain=math.sqrt(math.pi/(2*V))*erfcx(beta)
        bc=math.sqrt(math.pi/(2*V))*erfc(beta)
        numerical=quad(lambda t:math.exp(-Q*t-V*t*t/2),0,np.inf,epsabs=1e-13)[0]
        check(plain,numerical,f"B9 direct Gaussian transform Q={Q},V={V}")
        check(bc,math.exp(-beta*beta)*plain,f"B9 backtracking factor Q={Q},V={V}")
        require(0<=bc<=plain<=1/Q,f"B9 finite-Q bound Q={Q},V={V}")
    # Source S12 indices; chosen unit v=lambda_2=B0=1, a^2=1/3.
    for variance in [0.1,0.4,0.8]:
        for lp in [100.,1000.]:
            vx=variance/18
            plain=old.log_integral(lambda y:
                math.exp((3 if y<0 else -old.S)*y)
                *erfcx(math.exp(-y)/(lp*math.sqrt(2*vx))))
            factor=(1/18)*math.sqrt(math.pi/2)*(4/7)*variance/math.sqrt(vx)
            plain*=factor
            bc=old.rbd(variance,lp)
            require(0<=bc<=plain,f"B3 corrected <= plain var={variance},lp={lp}")
            out["rbd_plain_vs_bc"].append({"variance":variance,"lambda_parallel":lp,
                "corrected":bc,"plain":plain})
    for lp in [0,1e-4,0.01,1,100,10000]:
        kp=old.closed(lp)
        upper=min(.1**2*lp,1.5*.1)
        require(0<=kp<=upper*(1+1e-14),f"B6 bound lp={lp}")
        if lp:
            check((2/3)*kp+math.sqrt(kp/lp),.1,f"B5 residual lp={lp}")
        out["closed_bounds"].append({"lambda_parallel":lp,"lambda_perp":kp,"upper_bound":upper})
    # C12 is an estimator identity. These amplitudes/time are arbitrary tests.
    ag,tc=2.,3.
    derivative_mean=ag/2*tc**(-1/3)
    secant_mean=3*ag/4*tc**(-1/3)
    check(derivative_mean/secant_mean,2/3,"C12 derivative vs secant average")
    # G4/G6 identity for chosen h_C only; no historical angular default implied.
    for theta in [0.,.4,math.pi/2,2.4,math.pi]:
        hc=1.7
        u=1.5+.5*math.tanh(hc*(abs(theta-math.pi/2)-35*math.pi/180))
        g=1.5-.5*math.tanh(hc*(min(theta,math.pi-theta)-math.pi/2+35*math.pi/180))
        check(u,g,f"G4/G6 folded identity theta={theta}")
    paired=json.loads((ROOT/"paired_benchmark_points.json").read_text())
    for r in paired["nlgce"]:
        bound=r["a_prime_squared"]*r["nonlinear_parallel"]*(1-r["slab_fraction"])*r["epsilon_squared"]/2
        require(0<=r["nonlinear_perpendicular"]<=bound,f"E4 paired bound {r['label']}")
    out["synthetic_input_notes"]=[
        "RBD transform Q,V, finite-path states, closed-form kappa_FL/ell=0.1, GCD amplitudes, and angular width h_C=1.7 rad^-1 are chosen equation tests.",
        "The angular width test does not resolve the Corti source's multiplier 8.",
        "Nonanalytic-tail points fix alpha at the q=3 short-path root; they are not full finite-A self-consistent roots.",
        "All tolerances are engineering check budgets, not physical errors."]
    return out


def compare(actual,expected,path="review_fixtures"):
    if isinstance(actual,dict):
        require(actual.keys()==expected.keys(),path+" keys")
        for k in actual:compare(actual[k],expected[k],path+"."+k)
    elif isinstance(actual,list):
        require(len(actual)==len(expected),path+" length")
        for i,(a,e) in enumerate(zip(actual,expected)):compare(a,e,path+f"[{i}]")
    elif isinstance(actual,(float,int)):
        check(actual,expected,path,rtol=3e-10,atol=1e-14)
    else:
        require(actual==expected,path)


def verify_printed_constants(document,values):
    c=values["constants"]
    expected={
        "C(s)=0.118862354635443":c["C_s"],
        r"\pi C(s)=0.373417100111093":c["pi_C_s"],
        r"L_s/\ell_s=0.746834200222187":c["L_s_over_ell_s"]}
    for printed,val in expected.items():
        require(document.count(printed)==1,"Printed constant "+printed)
        require(printed.rsplit("=",1)[1]==format(val,".15f"),"Correct rounding "+printed)
    for val in [c["GCD_alpha"],c["GCD_secant_average_shape"]]:
        require(format(val,".15g") in document,"Printed GCD coefficient "+format(val,".15g"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate-fixtures",action="store_true")
    parser.add_argument("--write-results",action="store_true")
    parser.add_argument("--document",type=Path,default=ROOT.parent/"PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md")
    args=parser.parse_args()
    values=calculate()
    path=ROOT/"review_additions_fixtures.json"
    if args.generate_fixtures:
        path.write_text(json.dumps(values,indent=2,allow_nan=False)+"\n")
    compare(values,json.loads(path.read_text()))
    verify_printed_constants(args.document.read_text(),values)
    result={"state":"passed","specification_version":"2.1","checks_passed":len(CHECKS),
        "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
        "mpmath":mp.__version__,"precision_decimal_digits":mp.mp.dps,
        "physical_validation_claimed":False,"AMPS_tests_claimed":False,"checks":CHECKS}
    if args.generate_fixtures or args.write_results:
        (ROOT/"review_additions_results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="checks"},indent=2))


if __name__=="__main__":
    main()

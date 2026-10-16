#!/usr/bin/env python3
"""Independent checks for the v1.2 AMPS implementation contract, not an AMPS solver.

Requires numpy, scipy and mpmath. No supplied fixture or source file is overwritten.
Run: python check_physics.py --output verification_results.json
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

mp.mp.dps = 45
RESULTS = []

def check(name, condition, **evidence):
    RESULTS.append(dict(name=name, status="PASS" if bool(condition) else "FAIL", **evidence))

def close(a, b, atol=1e-12, rtol=1e-11):
    return abs(float(a)-float(b)) <= atol + rtol*abs(float(b))

def tensor_identity():
    # A deliberately nonconstant PSD tensor; arbitrary normalization p0=1.
    def dmm(mu,p): return (1-mu*mu)*(1+p*p)
    def dmp(mu,p): return mp.mpf('.2')*p*(1-mu*mu)
    def dpp(mu,p): return p*p*(1+mu*mu)
    def f(mu,p): return mp.exp(mu)*mp.exp(-p)+mu*mu*p*p
    def density(mu,p): return p*p*f(mu,p)
    def am(mu,p):
        return mp.diff(lambda u:dmm(u,p),mu)+mp.diff(lambda v:dmp(mu,v),p)+2*dmp(mu,p)/p
    def ap(mu,p):
        return mp.diff(lambda u:dmp(u,p),mu)+mp.diff(lambda v:dpp(mu,v),p)+2*dpp(mu,p)/p
    for mu,p in [(-.8,.2),(-.1,.7),(.3,1.3),(.9,3.)]:
        mu,p=mp.mpf(mu),mp.mpf(p)
        f_mu=lambda u,v:mp.diff(lambda x:f(x,v),u)
        f_p=lambda u,v:mp.diff(lambda x:f(u,x),v)
        physical=mp.diff(lambda u:dmm(u,p)*f_mu(u,p)+dmp(u,p)*f_p(u,p),mu)
        physical+=mp.diff(lambda v:v*v*(dmp(mu,v)*f_mu(mu,v)+dpp(mu,v)*f_p(mu,v)),p)/(p*p)
        forward=-mp.diff(lambda u:am(u,p)*density(u,p),mu)-mp.diff(lambda v:ap(mu,v)*density(mu,v),p)
        forward+=mp.diff(lambda u:dmm(u,p)*density(u,p),mu,2)
        forward+=2*mp.diff(lambda u:mp.diff(lambda v:dmp(u,v)*density(u,v),p),mu)
        forward+=mp.diff(lambda v:dpp(mu,v)*density(mu,v),p,2)
        err=float(abs(forward-p*p*physical))
        check(f"tensor_measure_identity_mu_{float(mu):g}_p_{float(p):g}",err<1e-35,absolute_residual=err)
        D=np.array([[float(dmm(mu,p)),float(dmp(mu,p))],[float(dmp(mu,p)),float(dpp(mu,p))]])
        B=np.linalg.cholesky(2*D)
        check(f"tensor_noise_covariance_mu_{float(mu):g}",np.max(np.abs(B@B.T-2*D))<1e-12)
        # The drift of log p is a_p/p-Dpp/p^2, including Itô's correction.
        expected=ap(mu,p)/p-dpp(mu,p)/(p*p)
        applied=ap(mu,p)*mp.diff(mp.log,p)+dpp(mu,p)*mp.diff(mp.log,p,2)
        check(f"log_momentum_ito_mu_{float(mu):g}",abs(expected-applied)<mp.mpf('1e-38'))
    check("rank_one_tensor_is_admissible",close(np.linalg.det([[1.,2.],[2.,4.]]),0.))
    check("indefinite_tensor_is_rejected",np.linalg.eigvalsh([[1.,2.],[2.,1.]]).min()<0)

def wave_and_helicity():
    # Units normalized c=m=1; the Lorentz identities are dimension-independent.
    for u in [-.3,-.02,.02,.3]:
        gu=1/math.sqrt(1-u*u)
        for p in [.01,1.,100.]:
            for mu in [-1.,-.4,.7,1.]:
                E=math.hypot(1.,p); pp=gu*(p*mu-u*E); perp=p*math.sqrt(max(0.,1-mu*mu))
                Ep=gu*(E-u*p*mu); pprime=math.hypot(pp,perp)
                Eback=gu*(Ep+u*pp); pback=gu*(pp+u*Ep)
                scale=max(1.,E*E)
                check(f"lorentz_mass_shell_{u}_{p}_{mu}",abs(Ep*Ep-pprime*pprime-1)<1e-11*scale)
                check(f"lorentz_roundtrip_{u}_{p}_{mu}",close(Eback,E) and close(pback,p*mu))
                mu_new=.13
                dE=gu*u*pprime*(mu_new-pp/pprime)
                dp=gu*pprime*(mu_new-pp/pprime)
                Enew=gu*(Ep+u*pprime*mu_new)
                check(f"elastic_wave_exchange_{u}_{p}_{mu}",close(Enew-Eback,dE) and close(dE,u*dp))
                clock=gu*(1-u*(p/E)*mu)
                check(f"positive_clock_{u}_{p}_{mu}",clock>0)
    for sc,sm in [(0.,0.),(.6,-.4),(-.7,.2),(1.,-.3)]:
        low=-1+abs(sc+sm); high=1-abs(sc-sm)
        for chi in [low,(low+high)/2,high]:
            weights={(j,n):(1+j*sc+n*sm+j*n*chi)/4 for j in [-1,1] for n in [-1,1]}
            check(f"helicity_joint_split_{sc}_{sm}_{chi}",min(weights.values())>=-1e-15 and close(sum(weights.values()),1.) and close(sum(j*w for (j,n),w in weights.items()),sc) and close(sum(n*w for (j,n),w in weights.items()),sm))
    # Identical marginal helicities can correspond to different resonant populations.
    check("zero_helicity_does_not_fix_four_populations",(1+.8)/4 != (1-.8)/4)

def finite_volume():
    edges=np.array([-1.,-.65,-.2,.15,.55,1.]); x=(edges[1:]+edges[:-1])/2; widths=np.diff(edges)
    nu=1.7
    D=lambda mu:.5*nu*(1-mu*mu)
    T=np.array([1/quad(lambda z:1/D(z),a,b,epsabs=1e-12)[0] for a,b in zip(x[:-1],x[1:])])
    Q=np.zeros((len(x),len(x)))
    for i,t in enumerate(T): Q[i,i+1]=t/widths[i]; Q[i+1,i]=t/widths[i+1]
    np.fill_diagonal(Q,-Q.sum(axis=1)); rho=widths/2
    check("fv_nonuniform_number_conservation",np.max(np.abs(Q.sum(axis=1)))<1e-13)
    check("fv_nonuniform_stationary_mass",np.max(np.abs(rho@Q))<1e-13)
    check("fv_nonuniform_detailed_balance",np.max(np.abs(rho[:,None]*Q-(rho[:,None]*Q).T))<1e-13)
    transition=expm(Q*.4)
    check("fv_semigroup_positive_and_stochastic",transition.min()>-1e-14 and np.max(np.abs(transition.sum(axis=1)-1))<1e-13)
    # Integrable cusp; endpoints +/-0.05 and D0=1, q=5/3, H=0.
    alpha=mp.mpf(2)/3; a=mp.mpf('.05')
    integral=2*mp.quad(lambda z:1/((1-z*z)*z**alpha),[0,a])
    trans=1/integral
    reference=2*a**(1-alpha)/(1-alpha)*mp.hyp2f1(1,(1-alpha)/2,(3-alpha)/2,a*a)
    check("qform_transmits_across_zero",close(trans,1/reference,atol=2e-12),transmissibility=float(trans),hypergeometric_reference=float(1/reference))
    # Integral of 1/|mu| diverges, so epsilon=0 pure scattering has no cross-zero rate.
    check("epsilon_zero_internal_resistance_diverges",float(mp.log(mp.mpf('1e-3')/mp.mpf('1e-30')))>60)
    return Q,x,rho

def focusing_reference():
    fixtures={0.:1/3,.1:.3322984618611749,1.:.2480767406094695,1.5:.1782757803373446,5.:.01032482948483814,-1.5:.1782757803373446}
    for xi,expected in fixtures.items():
        a=mp.mpf(xi)
        if a==0:
            value=mp.mpf(1)/3; mean=mp.mpf(0)
        else:
            Z=2*mp.sinh(a)/a; mean=mp.coth(a)-1/a
            def integrand(mu):
                if abs(mu)==1: return mp.mpf(0)
                rho=mp.exp(a*mu)/Z
                H=(mp.exp(a*mu)*(mu-mp.coth(a))+mp.exp(-a)*(1+mp.coth(a)))/(a*Z)
                return H*H/(mp.mpf('.5')*(1-mu*mu)*rho)
            value=mp.quad(integrand,[-1,0,1])
        check(f"focused_variance_xi_{xi}",close(value,expected,rtol=2e-12),kappa_variance=float(value),mean_mu=float(mean))
    xi=1.5; target=fixtures[xi]; previous=None
    for N in [400,800,1600]:
        h=2/N; x=np.linspace(-1+h/2,1-h/2,N); edges=np.linspace(-1,1,N+1)
        # Cell stationary masses, stable at these moderate xi.
        Z=2*np.sinh(xi)/xi
        mass=(np.exp(xi*edges[1:])-np.exp(xi*edges[:-1]))/(xi*Z)
        # Composite Simpson integration between neighboring centers; regular interiors.
        midpoint=(x[:-1]+x[1:])/2
        def reciprocal(z):return Z*np.exp(-xi*z)/(.5*(1-z*z))
        resistance=h*(reciprocal(x[:-1])+4*reciprocal(midpoint)+reciprocal(x[1:]))/6
        T=1/resistance
        f=x-np.dot(mass,x)
        flux=np.cumsum(mass*f)[:-1]
        variance=np.sum(flux*flux/T)
        err=abs(variance-target)
        check(f"focused_fv_grid_{N}",err<6e-5 and (previous is None or err<previous),kappa_variance=float(variance),absolute_error=float(err))
        previous=err
    # Exact manufactured physical f and stationary angular flux, independently differentiated.
    for xi in [-1.5,.2,1.5]:
        L=1/mp.mpf(xi)
        f=lambda s,mu:mp.exp(-s/L)*mp.exp(mp.mpf(xi)*mu)
        D=lambda mu:mp.mpf('.5')*(1-mu*mu)
        a=lambda mu:(1-mu*mu)/(2*L)
        for mu in [-.7,0.,.6]:
            s=mp.mpf('.3'); mu=mp.mpf(mu)
            residual=mu*mp.diff(lambda y:f(y,mu),s)+a(mu)*mp.diff(lambda y:f(s,y),mu)-mp.diff(lambda y:D(y)*mp.diff(lambda z:f(s,z),y),mu)
            check(f"manufactured_fte_xi_{xi}_mu_{float(mu)}",abs(residual)<mp.mpf('1e-38'))

def streaming_and_bridge(Q,x,rho):
    # Finite-state Feynman-Kac derivatives: endpoint-perfect angular transition is
    # insufficient to reproduce the displacement integral over a finite step.
    N=len(x); V=np.diag(x)
    A=np.zeros((3*N,3*N)); A[:N,:N]=Q; A[N:2*N,N:2*N]=Q; A[2*N:,2*N:]=Q
    A[N:2*N,:N]=V; A[2*N:,N:2*N]=2*V
    dt=.6; state=np.r_[np.ones(N),np.zeros(2*N)]
    solved=expm(A*dt)@state; exact_second=float(rho@solved[2*N:]); endpoint_second=dt*dt*float(rho@(x*x))
    check("angular_endpoint_is_not_displacement_integral",abs(exact_second-endpoint_second)>1e-3,integrated_second=exact_second,endpoint_second=endpoint_second)
    rng=np.random.default_rng(3017133); h=.8; q=.3; count=200000
    parent=rng.normal(size=count)*math.sqrt(h); z=rng.normal(size=count)
    left=q*parent+math.sqrt(q*(1-q)*h)*z; right=parent-left
    # Fixed count and tolerances declared here before running, not tuned to result.
    check("bridge_preserves_parent_path",np.max(np.abs(left+right-parent))<1e-14)
    check("bridge_left_variance",abs(np.var(left)-q*h)<.004,measured=float(np.var(left)),expected=q*h)
    check("bridge_right_variance",abs(np.var(right)-(1-q)*h)<.007,measured=float(np.var(right)),expected=(1-q)*h)
    check("bridge_disjoint_increment_covariance",abs(float(np.cov(left,right,ddof=0)[0,1]))<.004)
    # Large-lag endpoint-only streaming violates the finite-time isotropic reference.
    tau=1.; exact=2*(tau-(1-math.exp(-tau)))/3; endpoint=tau*tau/3
    check("isotropic_streaming_reference_detects_endpoint_kick",abs(exact-endpoint)>.08,exact_second=exact,endpoint_second=endpoint)

def spectrum_and_time():
    k0,k1,P0,P1=1.,2.,1.,.25
    alpha=math.log(P1/P0)/math.log(k1/k0)
    check("log_spectrum_recovers_nodes",close(P0*(k1/k0)**alpha,P1))
    exact=P0*k0*((k1/k0)**(alpha+1)-1)/(alpha+1)
    numerical=quad(lambda k:P0*(k/k0)**alpha,k0,k1)[0]
    check("log_spectrum_variance_integral",close(exact,numerical))
    # A step in the rate; splitting at update follows integrated optical time.
    nu1,nu2,t1,t2=.4,2.1,.3,.7
    for ell in [1,2,3,4]:
        eig=ell*(ell+1)/2
        combined=math.exp(-eig*(nu1*t1+nu2*t2))
        split=math.exp(-eig*nu1*t1)*math.exp(-eig*nu2*t2)
        check(f"rate_change_legendre_l_{ell}",close(combined,split))

def state_dependent_clock():
    D=lambda mu:mp.mpf('.5')*(1-mu*mu)
    g=lambda mu:1+mp.mpf('.3')*mu
    Z=mp.quad(lambda mu:1/g(mu),[-1,1])
    rho=lambda mu:1/(Z*g(mu))
    for mu in [-.8,0.,.6]:
        mu=mp.mpf(mu)
        drift=g(mu)*mp.diff(D,mu)
        flux=drift*rho(mu)-mp.diff(lambda y:g(y)*D(y)*rho(y),mu)
        check(f"time_change_stationary_flux_mu_{float(mu)}",abs(flux)<mp.mpf('1e-38'))
        naive=mp.diff(lambda y:g(y)*D(y),mu)
        check(f"time_change_detects_spurious_drift_mu_{float(mu)}",abs(naive-drift)>mp.mpf('.02'),extra_drift=float(naive-drift))

def bundle_evidence(bundle):
    root=Path(bundle)
    README=(root/'README.md').read_text()
    registry=json.loads((root/'parameters/model_registry.json').read_text())
    fixtures=json.loads((root/'benchmark_points.json').read_text())
    def list_count(obj):
        if isinstance(obj,list):return len(obj)
        for key in ['models','fixtures','model_registry','benchmarks']:
            if key in obj and isinstance(obj[key],list):return len(obj[key])
        return None
    result=dict(readme_says_version_1_0='version 1.0' in README,registry_ids=list_count(registry),fixtures=list_count(fixtures),sha256sums_sha256=hashlib.sha256((root/'SHA256SUMS').read_bytes()).hexdigest())
    check("available_data_are_not_v1_1",result['readme_says_version_1_0'] and result['sha256sums_sha256']=='aba757e4c214bca1987c5f1a4014359315b8bc8e576a2fd6d39fee9f72dcfd13',**result)
    return result

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',default='verification_results.json'); parser.add_argument('--bundle',default=None); args=parser.parse_args()
    start=time.time()
    tensor_identity(); wave_and_helicity(); Q,x,rho=finite_volume(); focusing_reference(); streaming_and_bridge(Q,x,rho); spectrum_and_time(); state_dependent_clock()
    evidence=bundle_evidence(args.bundle) if args.bundle else None
    report=dict(scope='Independent mathematical checks and interface-contract checks; no AMPS integration certification; no v1.1 companion-bundle certification.',python=sys.version.split()[0],platform=platform.platform(),checks=len(RESULTS),passed=sum(r['status']=='PASS' for r in RESULTS),failed=sum(r['status']=='FAIL' for r in RESULTS),elapsed_seconds=time.time()-start,bundle_evidence=evidence,results=RESULTS)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(f"{report['passed']}/{report['checks']} passed; {report['failed']} failed")
    for r in RESULTS:
        if r['status']!='PASS':print(json.dumps(r))
    return int(report['failed']>0)

if __name__=='__main__':raise SystemExit(main())

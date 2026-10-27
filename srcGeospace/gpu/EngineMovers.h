#ifndef _SRC_GEOSPACE_GPU_ENGINE_MOVERS_H_
#define _SRC_GEOSPACE_GPU_ENGINE_MOVERS_H_

//======================================================================================
// EngineMovers.h -- device port of gridless/GridlessParticleMovers.cpp.
//
// Ported functions (same arithmetic, same operand order):
//   BorisStep, HigueraCaryMagneticMidpointStep, HC4StepUnchecked, RKStepGeneric<N>,
//   GCStepGeneric<N>, EvaluateFieldGeometry, ComputeGuidingCenterAdiabaticity,
//   GuidingCenterRhs, ReconstructGuidingCenterMomentum, SelectPerpendicularUnitVector,
//   the HYBRID branch selector/commit, and StepParticleChecked().
//
// Differences that do not change results:
//   * The HYBRID context is a per-trajectory struct instead of thread_local storage.
//   * Butcher tableaus are returned by switch functions instead of static std::arrays;
//     the explicit stage loops still include every zero coefficient, as the legacy movers do.
//   * Length constants (1 m, 100 m) and the speed of light go through the unit policy;
//     they are compile-time identities in the SI build.
//======================================================================================

#include "EngineMath.h"
#include "EngineTypes.h"

namespace Geospace {
namespace GPU {

template <class R>
struct EgHybrid {
  int initialized;
  Vec3<R> xLaunch;
  R rInner;
  int totalOuterSteps;
  int consecutiveEligibleSteps;
  int currentlyUsingGC;
  int pendingDecisionValid;
  int pendingEligible;
  int pendingUseGC;
};

template <class R>
EG_HD void ResetHybrid(EgHybrid<R>& h,const Vec3<R>& xLaunch,R rInner) {
  h.initialized=1;
  h.xLaunch=xLaunch;
  h.rInner=rInner;
  h.totalOuterSteps=0;
  h.consecutiveEligibleSteps=0;
  h.currentlyUsingGC=0;
  h.pendingDecisionValid=0;
  h.pendingEligible=0;
  h.pendingUseGC=0;
}

template <class R,class U>
struct Movers {
  typedef Vec3<R> V;
  struct State { V x; V p; };
  struct GCState { V x; R pPar; };
  struct FieldGeometry { V B; V bHat; R Bmag; V gradBmag; V curvature; };

  static EG_HD R C() { return R(U::kSpeedOfLight); }
  static EG_HD R Len(double meters) { return R(ToLength<U>(meters)); }

  static EG_HD R GammaFromMomentum(const V& p,R m0) {
    const R p2=vdot(p,p);
    const R mc=m0*C();
    return std::sqrt(R(1.0)+p2/(mc*mc));
  }
  static EG_HD V VelocityFromMomentum(const V& p,R m0) {
    const R gamma=GammaFromMomentum(p,m0);
    return vmul(R(1.0)/(gamma*m0),p);
  }

  static EG_HD State add_state(const State& a,const State& b) {
    return { vadd(a.x,b.x), vadd(a.p,b.p) };
  }
  static EG_HD State mul_state(R s,const State& a) {
    return { vmul(s,a.x), vmul(s,a.p) };
  }
  static EG_HD GCState add_gc_state(const GCState& a,const GCState& b) {
    return { vadd(a.x,b.x), a.pPar+b.pPar };
  }
  static EG_HD GCState mul_gc_state(R s,const GCState& a) {
    return { vmul(s,a.x), s*a.pPar };
  }

  static EG_HD bool InsideInnerSphereM(const V& x,R rInner) {
    return (rInner>R(0.0)) && (vdot(x,x)<=rInner*rInner);
  }

  static EG_HD bool SegmentHitsInnerSphere(const V& a,const V& b,R rInner) {
    if (rInner<=R(0.0)) return false;
    if (InsideInnerSphereM(a,rInner) || InsideInnerSphereM(b,rInner)) return true;
    const V d=vsub(b,a);
    const R dd=vdot(d,d);
    if (dd<=R(0.0)) return false;
    R t=-vdot(a,d)/dd;
    if (t<R(0.0)) t=R(0.0);
    if (t>R(1.0)) t=R(1.0);
    const V q=vadd(a,vmul(t,d));
    return vdot(q,q)<=rInner*rInner;
  }

  static EG_HD bool CheckIntermediateSphereHit(const V& from,const V& to,R rInner,V& x) {
    if (SegmentHitsInnerSphere(from,to,rInner)) {
      x=to;
      return true;
    }
    return false;
  }

  template <class F>
  static EG_HD State LorentzRhs(const State& s,R q,R m0,F& field) {
    V B; field.GetB(s.x,B);
    const V v=VelocityFromMomentum(s.p,m0);
    return { v, vmul(q,vcross(v,B)) };
  }

  static EG_HD R Clamp(R x,R lo,R hi) { return (x<lo) ? lo : ((x>hi) ? hi : x); }

  static EG_HD V SelectPerpendicularUnitVector(const V& bHat) {
    V ref=(std::fabs(bHat.z)<R(0.9)) ? V{R(0.0),R(0.0),R(1.0)} : V{R(0.0),R(1.0),R(0.0)};
    V e1=vsub(ref,vmul(vdot(ref,bHat),bHat));
    const R n=vnorm(e1);
    if (n>R(0.0)) return vmul(R(1.0)/n,e1);
    ref=V{R(1.0),R(0.0),R(0.0)};
    e1=vsub(ref,vmul(vdot(ref,bHat),bHat));
    const R n2=vnorm(e1);
    return (n2>R(0.0)) ? vmul(R(1.0)/n2,e1) : V{R(1.0),R(0.0),R(0.0)};
  }

  template <class F>
  static EG_HD FieldGeometry EvaluateFieldGeometry(const V& x,F& field) {
    FieldGeometry g;
    field.GetB(x,g.B);
    g.Bmag=vnorm(g.B);
    if (g.Bmag<=R(0.0)) {
      g.bHat=V{R(0.0),R(0.0),R(1.0)};
      g.gradBmag=V{R(0.0),R(0.0),R(0.0)};
      g.curvature=V{R(0.0),R(0.0),R(0.0)};
      return g;
    }
    g.bHat=vmul(R(1.0)/g.Bmag,g.B);

    const R r=emax(Len(1.0),vnorm(x));
    const R h=emax(Len(100.0),R(1.0e-4)*r);

    V dbdx[3];
    g.gradBmag=V{R(0.0),R(0.0),R(0.0)};
    for (int idir=0;idir<3;++idir) {
      V dx{R(0.0),R(0.0),R(0.0)};
      if (idir==0) dx.x=h;
      if (idir==1) dx.y=h;
      if (idir==2) dx.z=h;

      V Bp,Bm;
      field.GetB(vadd(x,dx),Bp);
      field.GetB(vsub(x,dx),Bm);

      const R BpMag=vnorm(Bp);
      const R BmMag=vnorm(Bm);
      const R dB=(BpMag-BmMag)/(R(2.0)*h);
      if (idir==0) g.gradBmag.x=dB;
      if (idir==1) g.gradBmag.y=dB;
      if (idir==2) g.gradBmag.z=dB;

      const V bpHat=(BpMag>R(0.0)) ? vmul(R(1.0)/BpMag,Bp) : g.bHat;
      const V bmHat=(BmMag>R(0.0)) ? vmul(R(1.0)/BmMag,Bm) : g.bHat;
      dbdx[idir]=vmul(R(1.0)/(R(2.0)*h),vsub(bpHat,bmHat));
    }
    g.curvature=vadd(vadd(vmul(g.bHat.x,dbdx[0]),vmul(g.bHat.y,dbdx[1])),
                     vmul(g.bHat.z,dbdx[2]));
    return g;
  }

  template <class F>
  static EG_HD R ComputeGuidingCenterAdiabaticity(const V& x,const V& p,R q,F& field) {
    const R tiny=std::numeric_limits<R>::min();
    const R inf=std::numeric_limits<R>::infinity();
    const R qAbs=std::fabs(q);
    if (qAbs<=tiny) return inf;

    const FieldGeometry g=EvaluateFieldGeometry(x,field);
    if (g.Bmag<=tiny) return inf;

    const R pMag=vnorm(p);
    const R pPar=vdot(p,g.bHat);
    const R pPerp2=emax(R(0.0),pMag*pMag-pPar*pPar);
    const R pPerp=std::sqrt(pPerp2);
    const R rho=pPerp/(qAbs*g.Bmag);

    const R gradBmagNorm=vnorm(g.gradBmag);
    R LB=inf;
    if (gradBmagNorm>tiny) LB=g.Bmag/gradBmagNorm;

    const R curvatureNorm=vnorm(g.curvature);
    R Rc=inf;
    if (curvatureNorm>tiny) Rc=R(1.0)/curvatureNorm;

    const R Leff=emin(LB,Rc);
    if (!(Leff>tiny) || !efinite(Leff)) return R(0.0);
    return rho/Leff;
  }

  template <class F>
  static EG_HD bool EvaluateHybridInstantaneousEligibility(EgHybrid<R>& h,const V& x,
                                                           const V& p,R q,F& field) {
    if (!h.initialized) return false;
    const R tiny=std::numeric_limits<R>::min();
    const R qAbs=std::fabs(q);
    if (qAbs<=tiny) return false;

    V B;
    field.GetB(x,B);
    const R Bmag=vnorm(B);
    if (Bmag<=tiny) return false;

    const V bHat=vmul(R(1.0)/Bmag,B);
    const R pMag=vnorm(p);
    if (pMag<=tiny) return false;

    const R pPar=vdot(p,bHat);
    const R pParAbsFrac=std::fabs(pPar)/pMag;
    const R pPerp2=emax(R(0.0),pMag*pMag-pPar*pPar);
    const R pPerp=std::sqrt(pPerp2);
    const R rho=pPerp/(qAbs*Bmag);

    const R dLaunch=vnorm(vsub(x,h.xLaunch));
    R dInner=std::numeric_limits<R>::infinity();
    if (h.rInner>R(0.0)) dInner=std::fabs(vnorm(x)-h.rInner);

    const R epsilon=ComputeGuidingCenterAdiabaticity(x,p,q,field);

    const R epsilonEnter=R(0.03);
    const R epsilonStay=R(0.05);
    const R pParAbsFracMin=R(0.35);
    const R minDistanceFactorLaunch=R(8.0);
    const R minDistanceFactorInner=R(8.0);
    const int minWarmupStepsBeforeAnyGC=4;

    const bool epsilonOK=efinite(epsilon) &&
                         (epsilon<=(h.currentlyUsingGC ? epsilonStay : epsilonEnter));
    const bool pitchOK=(pParAbsFrac>=pParAbsFracMin);
    const bool launchOK=(dLaunch>=minDistanceFactorLaunch*emax(rho,Len(1.0)));
    const bool innerOK=(!efinite(dInner)) ||
                       (dInner>=minDistanceFactorInner*emax(rho,Len(1.0)));
    const bool warmupOK=(h.totalOuterSteps>=minWarmupStepsBeforeAnyGC);
    return epsilonOK && pitchOK && launchOK && innerOK && warmupOK;
  }

  template <class F>
  static EG_HD bool HybridPrepare(EgHybrid<R>& h,const V& x,const V& p,R q,F& field) {
    const bool eligibleNow=EvaluateHybridInstantaneousEligibility(h,x,p,q,field);
    const int minConsecutiveEligibleToEnterGC=2;
    const int prospectiveEligibleCount=eligibleNow ? (h.consecutiveEligibleSteps+1) : 0;
    bool useGC=false;
    if (eligibleNow) {
      if (h.currentlyUsingGC) useGC=true;
      else if (prospectiveEligibleCount>=minConsecutiveEligibleToEnterGC) useGC=true;
    }
    h.pendingDecisionValid=1;
    h.pendingEligible=eligibleNow ? 1 : 0;
    h.pendingUseGC=useGC ? 1 : 0;
    return useGC;
  }

  static EG_HD void HybridCommit(EgHybrid<R>& h) {
    const bool eligibleNow=h.pendingDecisionValid ? (h.pendingEligible!=0) : false;
    const bool useGC=h.pendingDecisionValid ? (h.pendingUseGC!=0) : false;
    ++h.totalOuterSteps;
    if (eligibleNow) {
      ++h.consecutiveEligibleSteps;
      h.currentlyUsingGC=useGC ? 1 : 0;
    }
    else {
      h.consecutiveEligibleSteps=0;
      h.currentlyUsingGC=0;
    }
    h.pendingDecisionValid=0;
    h.pendingEligible=0;
    h.pendingUseGC=0;
  }

  template <class F>
  static EG_HD bool HybridConsume(EgHybrid<R>& h,const V& x,const V& p,R q,F& field) {
    if (!h.pendingDecisionValid) HybridPrepare(h,x,p,q,field);
    return h.pendingUseGC!=0;
  }

  template <class F>
  static EG_HD V ReconstructGuidingCenterMomentum(const V& x,R pPar,R pMag,F& field) {
    V B;
    field.GetB(x,B);
    const R Bmag=vnorm(B);
    const V bHat=(Bmag>R(0.0)) ? vmul(R(1.0)/Bmag,B) : V{R(0.0),R(0.0),R(1.0)};
    const R pParClamped=Clamp(pPar,-pMag,pMag);
    const R pPerp2=emax(R(0.0),pMag*pMag-pParClamped*pParClamped);
    const R pPerp=std::sqrt(pPerp2);
    const V e1=SelectPerpendicularUnitVector(bHat);
    return vadd(vmul(pParClamped,bHat),vmul(pPerp,e1));
  }

  template <class F>
  static EG_HD GCState GuidingCenterRhs(const GCState& s,R q,R m0,R pMag,F& field) {
    const FieldGeometry g=EvaluateFieldGeometry(s.x,field);
    if (g.Bmag<=R(0.0) || std::fabs(q)<=std::numeric_limits<R>::min()) {
      const V pApprox=ReconstructGuidingCenterMomentum(s.x,s.pPar,pMag,field);
      return { VelocityFromMomentum(pApprox,m0), R(0.0) };
    }
    const R mc=m0*C();
    const R gamma=std::sqrt(R(1.0)+(pMag*pMag)/(mc*mc));
    const R pPar=Clamp(s.pPar,-pMag,pMag);
    const R pPerp2=emax(R(0.0),pMag*pMag-pPar*pPar);
    const R mu=pPerp2/(R(2.0)*gamma*m0*g.Bmag);
    const R vPar=pPar/(gamma*m0);

    const V vParallel=vmul(vPar,g.bHat);
    const V vGradB=vmul(mu/(q*g.Bmag*g.Bmag),vcross(g.B,g.gradBmag));
    const V vCurv=vmul((pPar*pPar)/(q*gamma*m0*g.Bmag),vcross(g.bHat,g.curvature));
    const V xDot=vadd(vadd(vParallel,vGradB),vCurv);
    const R dpParDt=-mu*vdot(g.bHat,g.gradBmag);
    return { xDot, dpParDt };
  }

  // Butcher tableaus: RK2 (midpoint), RK4 (classical), RK6 (7-stage).
  static EG_HD R TabA(int n,int i,int j) {
    if (n==2) return (i==1 && j==0) ? R(0.5) : R(0.0);
    if (n==4) {
      if (i==1 && j==0) return R(0.5);
      if (i==2 && j==1) return R(0.5);
      if (i==3 && j==2) return R(1.0);
      return R(0.0);
    }
    switch (i*7+j) {
      case 1*7+0: return R(1.0/3.0);
      case 2*7+1: return R(2.0/3.0);
      case 3*7+0: return R(1.0/12.0);
      case 3*7+1: return R(1.0/3.0);
      case 3*7+2: return R(-1.0/12.0);
      case 4*7+0: return R(-1.0/16.0);
      case 4*7+1: return R(9.0/8.0);
      case 4*7+2: return R(-3.0/16.0);
      case 4*7+3: return R(-3.0/8.0);
      case 5*7+1: return R(9.0/8.0);
      case 5*7+2: return R(-3.0/8.0);
      case 5*7+3: return R(-3.0/4.0);
      case 5*7+4: return R(0.5);
      case 6*7+0: return R(9.0/44.0);
      case 6*7+1: return R(-9.0/11.0);
      case 6*7+2: return R(63.0/44.0);
      case 6*7+3: return R(18.0/11.0);
      case 6*7+4: return R(0.0);
      case 6*7+5: return R(-16.0/11.0);
      default: return R(0.0);
    }
  }
  static EG_HD R TabB(int n,int i) {
    if (n==2) return (i==1) ? R(1.0) : R(0.0);
    if (n==4) {
      switch (i) {
        case 0: return R(1.0/6.0);
        case 1: return R(1.0/3.0);
        case 2: return R(1.0/3.0);
        default: return R(1.0/6.0);
      }
    }
    switch (i) {
      case 0: return R(11.0/120.0);
      case 1: return R(0.0);
      case 2: return R(27.0/40.0);
      case 3: return R(27.0/40.0);
      case 4: return R(-4.0/15.0);
      case 5: return R(-4.0/15.0);
      default: return R(11.0/120.0);
    }
  }

  template <int N,class F>
  static EG_HD bool RKStepGeneric(V& x,V& p,R q,R m0,R dt,F& field,R rInner) {
    const State y0{x,p};
    State k[N];
    V stageX[N+1];
    stageX[0]=x;
    for (int i=0;i<N;i++) {
      State yi=y0;
      for (int j=0;j<i;j++) yi=add_state(yi,mul_state(dt*TabA(N,i,j),k[j]));
      stageX[i+1]=yi.x;
      if (CheckIntermediateSphereHit(stageX[i],stageX[i+1],rInner,x)) return false;
      if (InsideInnerSphereM(yi.x,rInner)) { x=yi.x; p=yi.p; return false; }
      k[i]=LorentzRhs(yi,q,m0,field);
    }
    State y1=y0;
    for (int i=0;i<N;i++) y1=add_state(y1,mul_state(dt*TabB(N,i),k[i]));
    if (CheckIntermediateSphereHit(stageX[N],y1.x,rInner,x)) { p=y1.p; return false; }
    if (InsideInnerSphereM(y1.x,rInner)) { x=y1.x; p=y1.p; return false; }
    x=y1.x;
    p=y1.p;
    return true;
  }

  template <int N,class F>
  static EG_HD bool GCStepGeneric(V& x,V& p,R q,R m0,R dt,F& field,R rInner) {
    const R pMag0=vnorm(p);
    const FieldGeometry g0=EvaluateFieldGeometry(x,field);
    const R pPar0=(g0.Bmag>R(0.0)) ? vdot(p,g0.bHat) : p.z;

    const GCState y0{x,pPar0};
    GCState k[N];
    V stageX[N+1];
    stageX[0]=x;
    for (int i=0;i<N;++i) {
      GCState yi=y0;
      for (int j=0;j<i;++j) yi=add_gc_state(yi,mul_gc_state(dt*TabA(N,i,j),k[j]));
      stageX[i+1]=yi.x;
      if (CheckIntermediateSphereHit(stageX[i],stageX[i+1],rInner,x)) {
        p=ReconstructGuidingCenterMomentum(stageX[i+1],yi.pPar,pMag0,field);
        return false;
      }
      if (InsideInnerSphereM(yi.x,rInner)) {
        x=yi.x;
        p=ReconstructGuidingCenterMomentum(yi.x,yi.pPar,pMag0,field);
        return false;
      }
      k[i]=GuidingCenterRhs(yi,q,m0,pMag0,field);
    }
    GCState y1=y0;
    for (int i=0;i<N;++i) y1=add_gc_state(y1,mul_gc_state(dt*TabB(N,i),k[i]));
    if (CheckIntermediateSphereHit(stageX[N],y1.x,rInner,x)) {
      p=ReconstructGuidingCenterMomentum(y1.x,y1.pPar,pMag0,field);
      return false;
    }
    if (InsideInnerSphereM(y1.x,rInner)) {
      x=y1.x;
      p=ReconstructGuidingCenterMomentum(y1.x,y1.pPar,pMag0,field);
      return false;
    }
    x=y1.x;
    p=ReconstructGuidingCenterMomentum(y1.x,y1.pPar,pMag0,field);
    return true;
  }

  template <class F>
  static EG_HD void HigueraCaryMagneticMidpointStep(V& x,V& p,R q,R m0,R dt,F& field) {
    const V v0=VelocityFromMomentum(p,m0);
    const V x_mid=vadd(x,vmul(R(0.5)*dt,v0));
    V B;
    field.GetB(x_mid,B);
    const R gamma0=GammaFromMomentum(p,m0);
    const V t=vmul((q*dt)/(R(2.0)*gamma0*m0),B);
    const R t2=vdot(t,t);
    const V s=vmul(R(2.0)/(R(1.0)+t2),t);
    const V p_minus=p;
    const V p_prime=vadd(p_minus,vcross(p_minus,t));
    const V p_plus=vadd(p_minus,vcross(p_prime,s));
    p=p_plus;
    const V v1=VelocityFromMomentum(p,m0);
    x=vadd(x_mid,vmul(R(0.5)*dt,v1));
  }

  template <class F>
  static EG_HD void HC4StepUnchecked(V& x,V& p,R q,R m0,R dt,F& field) {
    const double cbrt2=1.25992104989487316477;
    const R w1=R(1.0/(2.0-cbrt2));
    const R w0=R(-cbrt2/(2.0-cbrt2));
    HigueraCaryMagneticMidpointStep(x,p,q,m0,w1*dt,field);
    HigueraCaryMagneticMidpointStep(x,p,q,m0,w0*dt,field);
    HigueraCaryMagneticMidpointStep(x,p,q,m0,w1*dt,field);
  }

  template <class F>
  static EG_HD void BorisStep(V& x,V& p,R q,R m0,R dt,F& field) {
    V B; field.GetB(x,B);
    const R gam=GammaFromMomentum(p,m0);
    const V v=vmul(R(1.0)/(gam*m0),p);
    x=vadd(x,vmul(R(0.5)*dt,v));
    const V t=vmul((q*dt)/(R(2.0)*gam*m0),B);
    const R t2=vdot(t,t);
    const V s=vmul(R(2.0)/(R(1.0)+t2),t);
    const V p_minus=p;
    const V p_prime=vadd(p_minus,vcross(p_minus,t));
    const V p_plus=vadd(p_minus,vcross(p_prime,s));
    p=p_plus;
    const R gam2=GammaFromMomentum(p,m0);
    const V v2=vmul(R(1.0)/(gam2*m0),p);
    x=vadd(x,vmul(R(0.5)*dt,v2));
  }

  // StepParticleChecked(mover,...)
  template <class F>
  static EG_HD bool StepChecked(int mover,EgHybrid<R>& hyb,V& x,V& p,R q,R m0,R dt,
                                F& field,R rInner) {
    switch (mover) {
      case kMoverHC4: {
        const int nCheck=4;
        const R dtSub=dt/R(nCheck);
        for (int i=0;i<nCheck;++i) {
          const V x0=x;
          HC4StepUnchecked(x,p,q,m0,dtSub,field);
          if (SegmentHitsInnerSphere(x0,x,rInner) || InsideInnerSphereM(x,rInner)) return false;
        }
        return true;
      }
      case kMoverRK2: return RKStepGeneric<2>(x,p,q,m0,dt,field,rInner);
      case kMoverRK4: return RKStepGeneric<4>(x,p,q,m0,dt,field,rInner);
      case kMoverRK6: return RKStepGeneric<7>(x,p,q,m0,dt,field,rInner);
      case kMoverGC2: return GCStepGeneric<2>(x,p,q,m0,dt,field,rInner);
      case kMoverGC4: return GCStepGeneric<4>(x,p,q,m0,dt,field,rInner);
      case kMoverGC6: return GCStepGeneric<7>(x,p,q,m0,dt,field,rInner);
      case kMoverHybrid: {
        const bool useGC=HybridConsume(hyb,x,p,q,field);
        bool ok=false;
        if (useGC) ok=GCStepGeneric<4>(x,p,q,m0,dt,field,rInner);
        else ok=RKStepGeneric<4>(x,p,q,m0,dt,field,rInner);
        HybridCommit(hyb);
        return ok;
      }
      case kMoverBoris:
      default: {
        const V x0=x;
        BorisStep(x,p,q,m0,dt,field);
        return !SegmentHitsInnerSphere(x0,x,rInner) && !InsideInnerSphereM(x,rInner);
      }
    }
  }
};

} // namespace GPU
} // namespace Geospace

#endif

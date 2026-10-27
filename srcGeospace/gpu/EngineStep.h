#ifndef _SRC_GEOSPACE_GPU_ENGINE_STEP_H_
#define _SRC_GEOSPACE_GPU_ENGINE_STEP_H_

//======================================================================================
// EngineStep.h -- resumable per-trajectory tracer.
//
// The legacy tracer integrates one trajectory to completion inside one function call
// (gridless TraceTrajectoryImpl, Mode3D TraceTrajectory3D).  The engine splits the same
// algorithm into
//
//   TraceBegin()     everything before the legacy `while` loop,
//   TraceAdvance(K)  at most K iterations of the loop body (plus the loop exit tests),
//
// keeping the complete loop state in a slot so that a kernel can advance thousands of
// trajectories by K steps each, return to the host to refill finished slots, and
// continue.  The loop body is a literal port of both legacy variants:
//
//   kVariantGridless  CutoffRigidityGridless.cpp  TraceTrajectoryImpl()
//   kVariantMode3D    CutoffRigidityMode3D.cpp    TraceTrajectory3D()
//
// including SelectLegacyCutoffTraceDt[3D]/SelectStructuredTraceDt[3D],
// TrajectoryBoundary::FindFirstEvent, PopulateExitKinematics/CompleteExitPitchAngle,
// the Shue outer boundary (SWMFCoupledAccess::ClassifyPosition/ClassifySegment), the
// Mode3D field-availability checks and the trap detector.  Time step control is per
// trajectory and fully adaptive, exactly as in the legacy tracer; there is no global time step.
//======================================================================================

#include "EngineField.h"
#include "EngineMovers.h"
#include "EngineTrap.h"

namespace Geospace {
namespace GPU {

template <class R>
struct EgSlot {
  Vec3<R> x,p;
  double t;                // trace time, internal units
  double s;                // path length, internal units
  int n;
  int done;
  int lastNode;            // mesh search hint
  EgHybrid<R> hyb;
  EgTrapState trap;
  EgResult res;
};

template <class R>
struct EgBox {
  R min[3],max[3],innerRadius;
};

template <class R>
struct EgEvent {
  int type;                // 0 none, 1 inner sphere, 2 outer box
  R fraction;
  R position[3];
};

enum : int { kEventNone=0, kEventInner=1, kEventOuter=2 };

// ---- TrajectoryBoundary.h port ----------------------------------------------------
template <class R> EG_HD bool IsFinitePointR(const R x[3]) {
  return efinite(x[0]) && efinite(x[1]) && efinite(x[2]);
}
template <class R> EG_HD bool InsideBoxR(const R x[3],const EgBox<R>& b,R tol) {
  return x[0]>=b.min[0]-tol && x[0]<=b.max[0]+tol &&
         x[1]>=b.min[1]-tol && x[1]<=b.max[1]+tol &&
         x[2]>=b.min[2]-tol && x[2]<=b.max[2]+tol;
}
template <class R> EG_HD bool InsideInnerSphereR(const R x[3],const EgBox<R>& b,R tol) {
  if (!(b.innerRadius>R(0.0))) return false;
  const R r=emax(R(0.0),b.innerRadius+tol);
  return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<=r*r;
}
template <class R> EG_HD void FillPositionR(EgEvent<R>& e,const R a[3],const R d[3]) {
  for (int i=0;i<3;++i) e.position[i]=a[i]+e.fraction*d[i];
}
template <class R>
EG_HD void ConsiderInnerSphereR(const R a[3],const R d[3],const EgBox<R>& b,R tol,
                                EgEvent<R>& best) {
  if (!(b.innerRadius>R(0.0))) return;
  if (InsideInnerSphereR(a,b,R(0.0))) {
    best.type=kEventInner;
    best.fraction=R(0.0);
    best.position[0]=a[0]; best.position[1]=a[1]; best.position[2]=a[2];
    return;
  }
  const R aa=d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
  if (!(aa>R(0.0))) return;
  const R segmentLength=std::sqrt(aa);
  const R fractionTol=(segmentLength>R(0.0)) ? emax(R(0.0),tol)/segmentLength : R(0.0);
  const R bb=R(2.0)*(a[0]*d[0]+a[1]*d[1]+a[2]*d[2]);
  const R cc=a[0]*a[0]+a[1]*a[1]+a[2]*a[2]-b.innerRadius*b.innerRadius;
  R disc=bb*bb-R(4.0)*aa*cc;
  if (disc<R(0.0)) return;
  disc=emax(R(0.0),disc);
  const R root=std::sqrt(disc);
  const R candidates[2]={(-bb-root)/(R(2.0)*aa),(-bb+root)/(R(2.0)*aa)};
  for (int ic=0;ic<2;++ic) {
    const R candidate=candidates[ic];
    if (candidate< -fractionTol || candidate>R(1.0)+fractionTol) continue;
    const R t=emax(R(0.0),emin(R(1.0),candidate));
    if (t>=best.fraction) continue;
    best.type=kEventInner;
    best.fraction=t;
    FillPositionR(best,a,d);
  }
}
template <class R>
EG_HD void ConsiderOuterBoxR(const R a[3],const R z[3],const R d[3],const EgBox<R>& b,
                             R tol,EgEvent<R>& best) {
  if (!InsideBoxR(a,b,R(0.0))) {
    if (R(0.0)<best.fraction) {
      best.type=kEventOuter;
      best.fraction=R(0.0);
      best.position[0]=a[0]; best.position[1]=a[1]; best.position[2]=a[2];
    }
    return;
  }
  if (InsideBoxR(z,b,R(0.0))) return;
  const R segmentLength=std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
  const R fractionTol=(segmentLength>R(0.0)) ? emax(R(0.0),tol)/segmentLength : R(0.0);
  for (int axis=0;axis<3;++axis) {
    if (d[axis]==R(0.0)) continue;
    const int side=(d[axis]>R(0.0)) ? 1 : 0;
    const R plane=side ? b.max[axis] : b.min[axis];
    const R candidate=(plane-a[axis])/d[axis];
    if (candidate< -fractionTol || candidate>R(1.0)+fractionTol) continue;
    const R t=emax(R(0.0),emin(R(1.0),candidate));
    if (t>=best.fraction) continue;
    R q[3]={a[0]+t*d[0],a[1]+t*d[1],a[2]+t*d[2]};
    bool onFace=true;
    for (int j=0;j<3;++j) {
      if (j==axis) continue;
      if (q[j]<b.min[j]-tol || q[j]>b.max[j]+tol) { onFace=false; break; }
    }
    if (!onFace) continue;
    best.type=kEventOuter;
    best.fraction=t;
    best.position[0]=q[0]; best.position[1]=q[1]; best.position[2]=q[2];
  }
}
template <class R>
EG_HD EgEvent<R> FindFirstEventR(const R start[3],const R end[3],const EgBox<R>& b,R tol) {
  EgEvent<R> best;
  best.type=kEventNone;
  best.fraction=std::numeric_limits<R>::infinity();
  best.position[0]=best.position[1]=best.position[2]=R(0.0);
  if (!IsFinitePointR(start) || !IsFinitePointR(end)) return best;
  const R d[3]={end[0]-start[0],end[1]-start[1],end[2]-start[2]};
  ConsiderInnerSphereR(start,d,b,tol,best);
  ConsiderOuterBoxR(start,end,d,b,tol,best);
  return best;
}

// ---- SWMFCoupledAccessContract.h Shue port (SI double) ----------------------------
enum : int { kPosInterior=0, kPosPhysicalEscape=1, kPosMeshUnavailable=2, kPosInvalid=3 };
enum : int { kSegNoExit=0, kSegPhysicalEscape=1, kSegMeshUnavailable=2, kSegInvalid=3 };

EG_HD bool ShueFinite(const double x[3]) {
  return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}
EG_HD bool ShueInsideBox(const EgShue& s,const double x[3],double tol) {
  if (!ShueFinite(x)) return false;
  for (int d=0;d<3;++d)
    if (x[d]<s.boxMin_m[d]-tol || x[d]>s.boxMax_m[d]+tol) return false;
  return true;
}
EG_HD double ShueRadiusMeters(const EgShue& s,double cosTheta) {
  cosTheta=emax(-1.0,emin(1.0,cosTheta));
  const double denominator=1.0+cosTheta;
  if (denominator<=64.0*std::numeric_limits<double>::epsilon())
    return std::numeric_limits<double>::infinity();
  return s.r0_Re*s.earthRadius_m*std::pow(2.0/denominator,s.alpha);
}
EG_HD double ShueSignedMarginMeters(const EgShue& s,const double x[3]) {
  if (!ShueFinite(x)) return -std::numeric_limits<double>::infinity();
  const double r=std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  if (r==0.0) return s.r0_Re*s.earthRadius_m;
  return ShueRadiusMeters(s,x[0]/r)-r;
}
EG_HD int ShueClassifyPosition(const EgShue& s,const double x[3],bool meshSampleAvailable,
                               double tol) {
  if (!ShueFinite(x)) return kPosInvalid;
  const bool insideBox=ShueInsideBox(s,x,tol);
  if (x[0]<=s.tailCapX_m+tol) return kPosPhysicalEscape;
  if (ShueSignedMarginMeters(s,x)<=tol) return kPosPhysicalEscape;
  if (!insideBox || !meshSampleAvailable) return kPosMeshUnavailable;
  return kPosInterior;
}
EG_HD void ShueInterpolate(const double a[3],const double b[3],double f,double out[3]) {
  for (int d=0;d<3;++d) out[d]=a[d]+f*(b[d]-a[d]);
}
// face: 0 none, 1 XMIN, 2 other face
EG_HD bool ShueFirstBoxExit(const EgShue& s,const double a[3],const double b[3],
                            double& fraction,int& face) {
  if (!ShueInsideBox(s,a,0.0)) return false;
  double best=2.0;
  int bestFace=0;
  for (int d=0;d<3;++d) {
    const double delta=b[d]-a[d];
    if (delta>0.0) {
      const double t=(s.boxMax_m[d]-a[d])/delta;
      if (t>=0.0 && t<=1.0 && t<best) { best=t; bestFace=2; }
    }
    else if (delta<0.0) {
      const double t=(s.boxMin_m[d]-a[d])/delta;
      if (t>=0.0 && t<=1.0 && t<best) { best=t; bestFace=(d==0) ? 1 : 2; }
    }
  }
  if (best>1.0 || ShueInsideBox(s,b,0.0)) return false;
  fraction=best;
  face=bestFace;
  return true;
}
EG_HD bool ShueFirstExit(const EgShue& s,const double a[3],const double b[3],double& fraction) {
  double previous=ShueSignedMarginMeters(s,a);
  if (!(previous>0.0)) { fraction=0.0; return true; }
  const int brackets=64;
  double lo=0.0,hi=0.0;
  bool found=false;
  double pt[3];
  for (int n=1;n<=brackets;++n) {
    const double t=static_cast<double>(n)/static_cast<double>(brackets);
    ShueInterpolate(a,b,t,pt);
    const double current=ShueSignedMarginMeters(s,pt);
    if (current<=0.0 && previous>0.0) {
      lo=static_cast<double>(n-1)/static_cast<double>(brackets);
      hi=t;
      found=true;
      break;
    }
    previous=current;
  }
  if (!found) return false;
  for (int iteration=0;iteration<64;++iteration) {
    const double mid=0.5*(lo+hi);
    ShueInterpolate(a,b,mid,pt);
    if (ShueSignedMarginMeters(s,pt)>0.0) lo=mid;
    else hi=mid;
  }
  fraction=hi;
  return true;
}
EG_HD int ShueClassifySegment(const EgShue& s,const double a[3],const double b[3],
                              double& fraction,double position[3]) {
  fraction=1.0;
  position[0]=position[1]=position[2]=0.0;
  if (!ShueFinite(a) || !ShueFinite(b)) return kSegInvalid;
  double boxFraction=2.0;
  int boxFace=0;
  const bool crossesBox=ShueFirstBoxExit(s,a,b,boxFraction,boxFace);
  double physicalFraction=2.0;
  double shueFraction=2.0;
  if (ShueFirstExit(s,a,b,shueFraction)) physicalFraction=shueFraction;
  if (crossesBox && boxFace==1 && boxFraction<physicalFraction) physicalFraction=boxFraction;
  if (physicalFraction<=1.0 && (!crossesBox || physicalFraction<=boxFraction+1.0e-13)) {
    fraction=physicalFraction;
    ShueInterpolate(a,b,physicalFraction,position);
    return kSegPhysicalEscape;
  }
  if (crossesBox) {
    fraction=boxFraction;
    ShueInterpolate(a,b,boxFraction,position);
    return kSegMeshUnavailable;
  }
  return kSegNoExit;
}

// ---- TrajectoryContract.h exit-state port (SI double) -----------------------------
EG_HD void ClearExitState(EgExitState& s) {
  for (int d=0;d<3;++d) { s.x_exit_m[d]=0.0; s.p_exit_SI[d]=0.0; s.v_exit_unit[d]=0.0; }
  s.cosAlpha=0.0; s.traceTimeAtExit_s=0.0; s.rigidityAtExit_GV=0.0; s.valid=0;
}
EG_HD bool PopulateExitKinematicsSI(EgExitState& state,const double eventPosition_m[3],
                                    const double pBefore[3],const double pAfter[3],
                                    double eventFraction,double traceTimeAtStepEnd_s,
                                    double acceptedStep_s,double absoluteCharge_C) {
  ClearExitState(state);
  if (!std::isfinite(eventFraction) || eventFraction<0.0 || eventFraction>1.0 ||
      !std::isfinite(traceTimeAtStepEnd_s) || traceTimeAtStepEnd_s<0.0 ||
      !std::isfinite(acceptedStep_s) || acceptedStep_s<0.0 ||
      !std::isfinite(absoluteCharge_C) || !(absoluteCharge_C>0.0)) return false;
  const double timeTolerance=16.0*std::numeric_limits<double>::epsilon()*
      emax(1.0,traceTimeAtStepEnd_s);
  if (acceptedStep_s>traceTimeAtStepEnd_s+timeTolerance) return false;
  double p2=0.0;
  for (int d=0;d<3;++d) {
    if (!std::isfinite(eventPosition_m[d]) || !std::isfinite(pBefore[d]) ||
        !std::isfinite(pAfter[d])) return false;
    state.x_exit_m[d]=eventPosition_m[d];
    state.p_exit_SI[d]=pBefore[d]+eventFraction*(pAfter[d]-pBefore[d]);
    p2+=state.p_exit_SI[d]*state.p_exit_SI[d];
  }
  if (!(p2>0.0) || !std::isfinite(p2)) return false;
  const double p=std::sqrt(p2);
  for (int d=0;d<3;++d) state.v_exit_unit[d]=state.p_exit_SI[d]/p;
  const double stepStart=traceTimeAtStepEnd_s-acceptedStep_s;
  state.traceTimeAtExit_s=emax(0.0,stepStart+eventFraction*acceptedStep_s);
  state.rigidityAtExit_GV=(p*299792458.0/absoluteCharge_C)*1.0e-9;
  return std::isfinite(state.traceTimeAtExit_s) && std::isfinite(state.rigidityAtExit_GV);
}
EG_HD bool CompleteExitPitchAngleSI(EgExitState& state,const double B[3]) {
  state.valid=0;
  double b2=0.0,v2=0.0;
  for (int d=0;d<3;++d) {
    if (!std::isfinite(B[d]) || !std::isfinite(state.v_exit_unit[d])) return false;
    b2+=B[d]*B[d];
    v2+=state.v_exit_unit[d]*state.v_exit_unit[d];
  }
  if (!(b2>0.0) || !std::isfinite(b2) || !std::isfinite(v2) ||
      std::fabs(v2-1.0)>1.0e-10) return false;
  const double invB=1.0/std::sqrt(b2);
  double cosine=0.0;
  for (int d=0;d<3;++d) cosine+=state.v_exit_unit[d]*B[d]*invB;
  if (!std::isfinite(cosine)) return false;
  state.cosAlpha=emax(-1.0,emin(1.0,cosine));
  state.valid=1;
  return true;
}

//======================================================================================
// Tracer
//======================================================================================
template <class R,class U>
struct Tracer {
  typedef Vec3<R> V;
  typedef Movers<R,U> M;

  // Per-call constants derived from the request/context.
  struct Consts {
    R q;              // signed trace charge (internal)
    R m0;             // rest mass (internal)
    double qabs_C;    // |q| SI, for the exit rigidity
    double maxTime;   // internal time
    double maxDist;   // internal length, <=0 disables
    double dtTrace;   // internal time
    double dtTrace_s; // SI, used by the legacy floor exactly as the legacy tracer
    int maxSteps;
    int mover;
    int policy;
    int variant;
    int captureExit;
    EgBox<R> box;
    R boundaryTol;    // internal length
    double boundaryTol_m;
    R boxRe[7];       // gridless legacy selector: domain in Re (internal scaled)
    R reInt;          // _EARTH__RADIUS_ in internal length units
  };

  static EG_HD void MakeConsts(const EgRequest& rq,const EgContext& cx,Consts& c) {
    c.q=R(ToCharge<U>(rq.qTrace_C));
    c.m0=R(ToMass<U>(rq.m0_kg));
    c.qabs_C=rq.qabs_C;
    c.maxTime=ToTime<U>(rq.maxTime_s);
    c.maxDist=(rq.maxDist_m>0.0) ? ToLength<U>(rq.maxDist_m) : -1.0;
    c.dtTrace=ToTime<U>(rq.dtTrace_s);
    c.dtTrace_s=rq.dtTrace_s;
    c.maxSteps=rq.maxSteps;
    c.mover=rq.mover;
    c.policy=rq.policy;
    c.variant=cx.variant;
    c.captureExit=rq.captureExit;
    for (int d=0;d<3;++d) {
      c.box.min[d]=R(ToLength<U>(cx.boxMin_m[d]));
      c.box.max[d]=R(ToLength<U>(cx.boxMax_m[d]));
    }
    c.box.innerRadius=R(ToLength<U>(cx.rInner_m));
    c.boundaryTol=R(ToLength<U>(cx.boundaryEventTolerance_m));
    c.boundaryTol_m=cx.boundaryEventTolerance_m;
    // DomainBoxRe values are already dimensionless (Re).  The legacy selector converts
    // them with *_EARTH__RADIUS_, so the internal length of 1 Re is used.
    for (int i=0;i<7;++i) c.boxRe[i]=R(cx.boxRe[i]);
    c.reInt=R(ToLength<U>(6371.0E3));
  }

  // ---- SelectLegacyCutoffTraceDt (gridless) -----------------------------------------
  template <class F>
  static EG_HD double LegacyDtGridless(const Consts& c,const EgContext& cx,F& field,
                                       const V& x,const V& p,double timeRemaining,
                                       bool useGC) {
    double dt=c.dtTrace;
    if (timeRemaining<dt) dt=timeRemaining;
    if (!cx.adaptiveDt) return dt;

    const R p2=vdot(p,p);
    const R mc=c.m0*M::C();
    const R gamma=std::sqrt(R(1.0)+p2/(mc*mc));
    const R pMag=std::sqrt(emax(R(0.0),p2));
    const R vMag=(gamma>R(0.0) && c.m0>R(0.0)) ? pMag/(gamma*c.m0) : R(0.0);

    V B; field.GetB(x,B);
    const R Bmag=vnorm(B);
    if (!useGC) {
      const R dphiMax=R(0.15);
      if (Bmag>R(0.0)) {
        const R omega=std::fabs(c.q)*Bmag/(gamma*c.m0);
        if (omega>R(0.0)) dt=emin(dt,double(dphiMax/omega));
      }
    }

    const R Re=c.reInt;
    const V xRe{x.x/Re,x.y/Re,x.z/Re};
    const R rRe=std::sqrt(xRe.x*xRe.x+xRe.y*xRe.y+xRe.z*xRe.z);
    R dInner=Sentinel<R>::kBig;
    if (rRe>c.boxRe[6]) dInner=(rRe-c.boxRe[6])*Re;

    R dBox=Sentinel<R>::kBig;
    const bool inside=(xRe.x>=c.boxRe[0] && xRe.x<=c.boxRe[1] &&
                       xRe.y>=c.boxRe[2] && xRe.y<=c.boxRe[3] &&
                       xRe.z>=c.boxRe[4] && xRe.z<=c.boxRe[5]);
    if (inside) {
      const R dxp=(c.boxRe[1]-xRe.x)*Re;
      const R dxm=(xRe.x-c.boxRe[0])*Re;
      const R dyp=(c.boxRe[3]-xRe.y)*Re;
      const R dym=(xRe.y-c.boxRe[2])*Re;
      const R dzp=(c.boxRe[5]-xRe.z)*Re;
      const R dzm=(xRe.z-c.boxRe[4])*Re;
      dBox=emin(emin(emin(dxp,dxm),emin(dyp,dym)),emin(dzp,dzm));
    }

    const R dNear=emin(dInner,dBox);
    R vFor=vMag;
    if (useGC && Bmag>R(0.0)) {
      const V bHat=vmul(R(1.0)/Bmag,B);
      const R pPar=vdot(p,bHat);
      const R vParAbs=std::fabs(pPar)/(gamma*c.m0);
      vFor=emax(vParAbs,R(ToSpeed<U>(1.0e-12)));
    }
    if (vFor>R(0.0) && dNear<Sentinel<R>::kBigLimit)
      dt=emin(dt,double(R(0.20)*dNear/vFor));

    const R safeSpeed=emax(vFor,R(ToSpeed<U>(1.0e-12)));
    const double dtFloor=emax(
        ToTime<U>(emax(1.0e-12,1.0e-9*emax(c.dtTrace_s,1.0))),
        double(R(ToLength<U>(100.0e3))/safeSpeed));
    dt=emax(dtFloor,dt);
    if (timeRemaining>0.0) dt=emin(dt,timeRemaining);
    return dt;
  }

  // ---- SelectLegacyCutoffTraceDt3D (Mode3D) -----------------------------------------
  template <class F>
  static EG_HD double LegacyDtMode3D(const Consts& c,const EgContext& cx,F& field,
                                     const V& x,const V& p,double timeRemaining,
                                     bool useGC) {
    double dt=c.dtTrace;
    if (timeRemaining<dt) dt=timeRemaining;
    if (!cx.adaptiveDt) return dt;

    const R p2=vdot(p,p);
    const R mc=c.m0*M::C();
    const R gamma=std::sqrt(R(1.0)+p2/(mc*mc));
    const R pMag=std::sqrt(emax(R(0.0),p2));
    const R vMag=(gamma>R(0.0) && c.m0>R(0.0)) ? pMag/(gamma*c.m0) : R(0.0);

    V B; field.GetB(x,B);
    const R Bmag=vnorm(B);
    if (!useGC) {
      const R gyroAngleLimit=R(0.15);
      if (Bmag>R(0.0)) {
        const R omega=std::fabs(c.q)*Bmag/(gamma*c.m0);
        if (omega>R(0.0)) dt=emin(dt,double(gyroAngleLimit/omega));
      }
    }

    const R dxp=c.box.max[0]-x.x;
    const R dxm=x.x-c.box.min[0];
    const R dyp=c.box.max[1]-x.y;
    const R dym=x.y-c.box.min[1];
    const R dzp=c.box.max[2]-x.z;
    const R dzm=x.z-c.box.min[2];
    // std::min({dxp,dxm,dyp,dym,dzp,dzm}) == left fold of (b<a)?b:a
    R dBox=dxp;
    dBox=emin(dBox,dxm); dBox=emin(dBox,dyp); dBox=emin(dBox,dym);
    dBox=emin(dBox,dzp); dBox=emin(dBox,dzm);
    const R dBoxSafe=(dBox<R(0.0)) ? Sentinel<R>::kBig : dBox;
    const R dInner=emax(R(0.0),vnorm(x)-c.box.innerRadius);
    const R dNear=emin(dBoxSafe,dInner);

    R vFor=vMag;
    if (useGC && Bmag>R(0.0)) {
      const V bHat=vmul(R(1.0)/Bmag,B);
      const R pPar=vdot(p,bHat);
      const R vParAbs=std::fabs(pPar)/(gamma*c.m0);
      vFor=emax(vParAbs,R(ToSpeed<U>(1.0e-12)));
    }
    if (vFor>R(0.0) && dNear<Sentinel<R>::kBigLimit)
      dt=emin(dt,double(R(0.20)*dNear/vFor));

    const R safeSpeed=emax(vFor,R(ToSpeed<U>(1.0e-12)));
    const double dtFloor=emax(
        ToTime<U>(emax(1.0e-12,1.0e-9*emax(c.dtTrace_s,1.0))),
        double(R(ToLength<U>(100.0e3))/safeSpeed));
    dt=emax(dtFloor,dt);
    if (timeRemaining>0.0) dt=emin(dt,timeRemaining);
    return dt;
  }

  static EG_HD double FinalizeUpperBoundedStep(double candidate,double timeRemaining) {
    if (!(candidate>0.0) || !std::isfinite(candidate))
      return std::numeric_limits<double>::quiet_NaN();
    if (timeRemaining>0.0) candidate=emin(candidate,timeRemaining);
    if (!(candidate>0.0) || !std::isfinite(candidate))
      return std::numeric_limits<double>::quiet_NaN();
    return candidate;
  }

  // ---- SelectStructuredTraceDt[3D] (identical in both variants) ---------------------
  template <class F>
  static EG_HD double StructuredDt(const Consts& c,const EgContext& cx,F& field,
                                   const V& x,const V& p,double timeRemaining,bool useGC) {
    double dt=c.dtTrace;
    if (!cx.adaptiveDt) return FinalizeUpperBoundedStep(dt,timeRemaining);
    const R p2=vdot(p,p);
    const R mc=c.m0*M::C();
    const R gamma=std::sqrt(R(1.0)+p2/(mc*mc));
    V B; field.GetB(x,B);
    const R Bmag=vnorm(B);
    if (!useGC) {
      const R dphiMax=R(0.15);
      if (Bmag>R(0.0)) {
        const R omega=std::fabs(c.q)*Bmag/(gamma*c.m0);
        if (omega>R(0.0)) dt=emin(dt,double(dphiMax/omega));
      }
    }
    return FinalizeUpperBoundedStep(dt,timeRemaining);
  }

  template <class F>
  static EG_HD double SelectDt(const Consts& c,const EgContext& cx,F& field,const V& x,
                               const V& p,double timeRemaining,bool useGC) {
    if (c.policy==kPolicyLegacy) {
      return (c.variant==kVariantMode3D)
          ? LegacyDtMode3D(c,cx,field,x,p,timeRemaining,useGC)
          : LegacyDtGridless(c,cx,field,x,p,timeRemaining,useGC);
    }
    return StructuredDt(c,cx,field,x,p,timeRemaining,useGC);
  }

  // PopulateOuterExit / PopulateExit.  eventPos is in internal units.
  template <class F>
  static EG_HD bool PopulateExit(const Consts& c,EgSlot<R>& s,F& field,
                                 const R eventPosInt[3],R fraction,
                                 const V& xBefore,const V& xAfter,
                                 const V& pBefore,const V& pAfter,double acceptedStep) {
    if (!c.captureExit) return true;
    double pos[3],pb[3],pa[3];
    pos[0]=FromLength<U>(eventPosInt[0]); pos[1]=FromLength<U>(eventPosInt[1]);
    pos[2]=FromLength<U>(eventPosInt[2]);
    pb[0]=FromMomentum<U>(pBefore.x); pb[1]=FromMomentum<U>(pBefore.y);
    pb[2]=FromMomentum<U>(pBefore.z);
    pa[0]=FromMomentum<U>(pAfter.x); pa[1]=FromMomentum<U>(pAfter.y);
    pa[2]=FromMomentum<U>(pAfter.z);
    if (!PopulateExitKinematicsSI(s.res.exit,pos,pb,pa,double(fraction),
                                  FromTime<U>(s.t),FromTime<U>(acceptedStep),c.qabs_C))
      return false;
    const R a=fraction;
    const V chord=vsub(xAfter,xBefore);
    const R chordLength=vnorm(chord);
    R insetLength=R(ToLength<U>(emax(1.0,c.boundaryTol_m)));
    if constexpr (!std::is_same<R,double>::value) {
      // Float: the legacy tracer samples the exit field max(1 m, boundary tolerance) inside the
      // event point.  That inset is below the float spacing of positions near the outer
      // box (about 24 m at 40 Re in Re units), so xEval can round onto or past the face,
      // where the mesh has no leaf, B is zero and the pitch angle is undefined
      // (INVALID_FIELD, reason 10 in C19).  Use at least 64 ulps of the event position;
      // the double build keeps the legacy inset.
      const R scale=emax(std::fabs(eventPosInt[0]),
                         emax(std::fabs(eventPosInt[1]),std::fabs(eventPosInt[2])));
      insetLength=emax(insetLength,R(64)*std::numeric_limits<R>::epsilon()*scale);
    }
    R insetFraction=R(0.0);
    if (chordLength>R(0.0))
      insetFraction=emin(a,insetLength/chordLength);
    const V xEval=vadd(xBefore,vmul(emax(R(0.0),a-insetFraction),chord));
    if (c.variant==kVariantMode3D) field.ResetAvailability();
    V Bexit; field.GetB(xEval,Bexit);
    if (c.variant==kVariantMode3D && !field.Available()) return false;
    const double b[3]={FromField<U>(Bexit.x),FromField<U>(Bexit.y),FromField<U>(Bexit.z)};
    return CompleteExitPitchAngleSI(s.res.exit,b);
  }

  static EG_HD void ToSI(const V& v,double scale,double out[3]) {
    out[0]=double(v.x)*scale; out[1]=double(v.y)*scale; out[2]=double(v.z)*scale;
  }

  template <class F>
  static EG_HD bool TrapUpdate(EgTrapDetector& trap,F& field,const V& x,const V& p,
                               const V& B,double t) {
    double xs[3],ps[3],bs[3];
    if constexpr (U::kSI) {
      xs[0]=x.x; xs[1]=x.y; xs[2]=x.z;
      ps[0]=p.x; ps[1]=p.y; ps[2]=p.z;
      bs[0]=B.x; bs[1]=B.y; bs[2]=B.z;
    }
    else {
      ToSI(x,U::kLength_m,xs); ToSI(p,U::kMomentum_SI,ps); ToSI(B,U::kField_T,bs);
    }
    (void)field;
    return trap.Update(xs,ps,bs,FromTime<U>(t));
  }

  static EG_HD void Finish(EgSlot<R>& s,EgTrapDetector& trap,int termination) {
    EgResult& r=s.res;
    const EgTrapState& ts=s.trap;
    r.termination=termination;
    r.traceTime_s=FromTime<U>(s.t);
    r.traceDistance_m=FromLength<U>(s.s);
    r.steps=s.n;
    r.mirrorPoints=ts.mirrorPoints;
    r.bounceCycles=ts.bounceCycles;
    r.driftRevolutions=ts.driftRevolutions;
    r.driftAngle_rad=trap.driftAngleRadians();
    r.driftMeanRadiusChange_m=ts.lastDriftMeanRadiusChange_m;
    r.driftMeanRadiusSpan_m=ts.lastDriftMeanRadiusSpan_m;
    r.driftStableComparisons=ts.stableDriftComparisons;
    r.driftRecurrenceLagRevolutions=ts.driftRecurrenceLagRevolutions;
    r.minimumOuterMargin_m=ts.minimumOuterMargin;
    r.lastBounceEnvelopeChange_m=ts.lastBounceEnvelopeChange_m;
    const EgDriftDiag& d=ts.lastDriftComparison;
    r.driftComparisonAvailable=d.available;
    r.driftPreviousCoveragePassed=d.previousCoveragePassed;
    r.driftCurrentCoveragePassed=d.currentCoveragePassed;
    r.driftSecularDriftPassed=d.secularDriftPassed;
    r.driftProfilesRecurred=d.profilesRecurred;
    r.driftPreviousCoverageFraction=d.previousCoverageFraction;
    r.driftCurrentCoverageFraction=d.currentCoverageFraction;
    r.driftCommonBins=d.commonBins;
    r.driftMatchedBins=d.matchedBins;
    r.driftRadialMismatchBins=d.radialMismatchBins;
    r.driftLatitudeMismatchBins=d.latitudeMismatchBins;
    r.driftPitchMismatchBins=d.pitchMismatchBins;
    r.driftMatchedBinFraction=d.matchedBinFraction;
    r.driftMaximumNormalizedRadialResidual=d.maximumNormalizedRadialResidual;
    r.driftMaximumNormalizedLatitudeResidual=d.maximumNormalizedLatitudeResidual;
    r.driftMaximumNormalizedPitchResidual=d.maximumNormalizedPitchResidual;
    r.trapMechanism=ts.mechanism;
    r.momentumRelativeSpread=trap.momentumRelativeSpread();
    r.turnStateCount=ts.turnCount;
    r.finalX_m[0]=FromLength<U>(s.x.x); r.finalX_m[1]=FromLength<U>(s.x.y);
    r.finalX_m[2]=FromLength<U>(s.x.z);
    if (ts.turnOverflow && r.error==kErrNone) r.error=kErrTurnStateOverflow;
    s.done=1;
  }

  static EG_HD void MergeStats(EgSlot<R>& s,const EgDipoleErrorStats& st) {
    if (st.count==0) return;
    EgResult& r=s.res;
    r.dipoleErrorCount+=st.count;
    r.dipoleErrorSum+=st.sum;
    if (st.maxValue>r.dipoleErrorMax) {
      r.dipoleErrorMax=st.maxValue;
      r.dipoleErrorMaxLocation_m[0]=st.maxLocation[0];
      r.dipoleErrorMaxLocation_m[1]=st.maxLocation[1];
      r.dipoleErrorMaxLocation_m[2]=st.maxLocation[2];
    }
  }

  // Everything that precedes the legacy integration loop.  Returns true when the
  // trajectory already finished (Mode3D InvalidField at the trap initialization).
  template <class F>
  static EG_HD bool Begin(EgSlot<R>& s,const EgRequest& rq,const EgContext& cx,F& field,
                          EgTrapDetector& trap) {
    Consts c;
    MakeConsts(rq,cx,c);
    EgResult& r=s.res;
    r.tag=rq.tag;
    r.termination=kTermNumericalFailure;
    r.steps=0; r.traceTime_s=0.0; r.traceDistance_m=0.0;
    ClearExitState(r.exit);
    r.error=kErrNone;
    r.invalidReason=0;
    r.unavailableKind=0;
    r.unavailableX_m[0]=r.unavailableX_m[1]=r.unavailableX_m[2]=0.0;
    r.turnStateCount=0;
    r.dipoleErrorCount=0; r.dipoleErrorSum=0.0; r.dipoleErrorMax=-1.0;
    r.dipoleErrorMaxLocation_m[0]=r.dipoleErrorMaxLocation_m[1]=
        r.dipoleErrorMaxLocation_m[2]=0.0;

    // pMag = MomentumFromRigidity_GV(R_GV,qabs) evaluated in SI, then p=pMag*v0_unit.
    const double pMag_SI=(rq.R_GV*1.0e9*rq.qabs_C)/299792458.0;
    if constexpr (U::kSI) {
      s.p.x=pMag_SI*rq.v0_unit[0]; s.p.y=pMag_SI*rq.v0_unit[1]; s.p.z=pMag_SI*rq.v0_unit[2];
      s.x.x=rq.x0_m[0]; s.x.y=rq.x0_m[1]; s.x.z=rq.x0_m[2];
    }
    else {
      const double pm=ToMomentum<U>(pMag_SI);
      s.p.x=R(pm*rq.v0_unit[0]); s.p.y=R(pm*rq.v0_unit[1]); s.p.z=R(pm*rq.v0_unit[2]);
      s.x.x=R(ToLength<U>(rq.x0_m[0])); s.x.y=R(ToLength<U>(rq.x0_m[1]));
      s.x.z=R(ToLength<U>(rq.x0_m[2]));
    }
    s.t=0.0; s.s=0.0; s.n=0; s.done=0;
    trap.Reset();
    ResetHybrid(s.hyb,s.x,c.box.innerRadius);

    if (cx.trap.enabled) {
      if (c.variant==kVariantMode3D) field.ResetAvailability();
      V B0; field.GetB(s.x,B0);
      if (c.variant==kVariantMode3D && !field.Available()) {
        s.res.invalidReason=1; Finish(s,trap,kTermInvalidField);
        return true;
      }
      TrapUpdate(trap,field,s.x,s.p,B0,0.0);
    }
    return false;
  }

  // Up to maxIterations iterations of the legacy loop.  Returns true when finished.
  template <class F>
  static EG_HD bool Advance(EgSlot<R>& s,const EgRequest& rq,const EgContext& cx,F& field,
                            EgTrapDetector& trap,int maxIterations) {
    Consts c;
    MakeConsts(rq,cx,c);
    const bool mode3d=(c.variant==kVariantMode3D);
    const bool shue=mode3d && cx.shue.enabled;
    const bool structured=(c.policy==kPolicyStructured);
    V x=s.x,p=s.p;

#define EG_FINISH(code) do { s.x=x; s.p=p; Finish(s,trap,(code)); return true; } while (0)

    for (int it=0;it<maxIterations;++it) {
      if (!(s.n<c.maxSteps && s.t<c.maxTime && (c.maxDist<=0.0 || s.s<c.maxDist))) {
        s.x=x; s.p=p;
        if (s.n>=c.maxSteps) { Finish(s,trap,kTermStepLimit); return true; }
        if (c.maxDist>0.0 && s.s>=c.maxDist) { Finish(s,trap,kTermDistanceLimit); return true; }
        if (s.t>=c.maxTime) { Finish(s,trap,kTermTimeLimit); return true; }
        Finish(s,trap,kTermNumericalFailure);
        return true;
      }

      const R xArr[3]={x.x,x.y,x.z};
      if (!IsFinitePointR(xArr) || !efinite(p.x) || !efinite(p.y) || !efinite(p.z))
        EG_FINISH(kTermNumericalFailure);
      if (InsideInnerSphereR(xArr,c.box,R(0.0))) EG_FINISH(kTermInnerBoundaryForbidden);

      if (shue) {
        double xs[3];
        ToSI(x,U::kLength_m,xs);
        const int disposition=ShueClassifyPosition(cx.shue,xs,true,c.boundaryTol_m);
        if (disposition==kPosInvalid || disposition==kPosMeshUnavailable)
          { s.res.invalidReason=2; EG_FINISH(kTermInvalidField); }
        if (disposition==kPosPhysicalEscape) {
          if (!PopulateExit(c,s,field,xArr,R(0.0),x,x,p,p,0.0))
            { s.res.invalidReason=3; EG_FINISH(kTermInvalidField); }
          EG_FINISH(kTermOuterBoundaryAllowed);
        }
      }
      else if (!InsideBoxR(xArr,c.box,R(0.0))) {
        if (!PopulateExit(c,s,field,xArr,R(0.0),x,x,p,p,0.0))
          { s.res.invalidReason=4; EG_FINISH(kTermInvalidField); }
        EG_FINISH(kTermOuterBoundaryAllowed);
      }

      const double timeRemaining=c.maxTime-s.t;
      if (mode3d) field.ResetAvailability();
      bool useGC=false;
      if (c.mover==kMoverGC2 || c.mover==kMoverGC4 || c.mover==kMoverGC6) useGC=true;
      else if (c.mover==kMoverHybrid) useGC=M::HybridPrepare(s.hyb,x,p,c.q,field);
      if (mode3d && !field.Available()) { s.res.invalidReason=5; EG_FINISH(kTermInvalidField); }

      const double dt=SelectDt(c,cx,field,x,p,timeRemaining,useGC);
      if (mode3d && !field.Available()) { s.res.invalidReason=6; EG_FINISH(kTermInvalidField); }
      if (!(dt>0.0) || !std::isfinite(dt)) EG_FINISH(kTermInvalidTimeStep);

      const V xBefore=x;
      const V pBefore=p;
      if (mode3d) field.ResetAvailability();
      const bool stepAccepted=M::StepChecked(c.mover,s.hyb,x,p,c.q,c.m0,R(dt),field,
                                             c.box.innerRadius);
      if (mode3d && !field.Available()) { s.res.invalidReason=7; EG_FINISH(kTermInvalidField); }
      if (!stepAccepted) EG_FINISH(kTermInnerBoundaryForbidden);

      const R xBeforeArr[3]={xBefore.x,xBefore.y,xBefore.z};
      const R xAfterArr[3]={x.x,x.y,x.z};
      if (!IsFinitePointR(xAfterArr) || !efinite(p.x) || !efinite(p.y) || !efinite(p.z))
        EG_FINISH(kTermNumericalFailure);

      EgEvent<R> event;
      event.type=kEventNone;
      event.fraction=std::numeric_limits<R>::infinity();
      event.position[0]=event.position[1]=event.position[2]=R(0.0);
      if (structured || shue)
        event=FindFirstEventR(xBeforeArr,xAfterArr,c.box,c.boundaryTol);

      int outerDisposition=kSegNoExit;
      double outerFraction=1.0;
      double outerPosition[3]={0.0,0.0,0.0};
      if (shue) {
        double a[3],b[3];
        ToSI(xBefore,U::kLength_m,a);
        ToSI(x,U::kLength_m,b);
        outerDisposition=ShueClassifySegment(cx.shue,a,b,outerFraction,outerPosition);
      }

      const R segmentLength=vnorm(vsub(x,xBefore));
      if (mode3d) { if (efinite(segmentLength)) s.s+=segmentLength; }
      else s.s+=efinite(segmentLength) ? segmentLength : R(0.0);
      s.t+=dt;
      ++s.n;

      if (shue) {
        if (event.type==kEventInner &&
            (outerDisposition==kSegNoExit || double(event.fraction)<=outerFraction))
          EG_FINISH(kTermInnerBoundaryForbidden);
        if (outerDisposition==kSegPhysicalEscape) {
          bool ok;
          if constexpr (U::kSI) {
            // The legacy tracer copies outerEvent.position_m into the event: no unit round trip.
            ok=PopulateExitSIPosition(c,s,field,outerPosition,R(outerFraction),
                                      xBefore,x,pBefore,p,dt);
          }
          else {
            R pos[3];
            for (int d=0;d<3;++d) pos[d]=R(ToLength<U>(outerPosition[d]));
            ok=PopulateExit(c,s,field,pos,R(outerFraction),xBefore,x,pBefore,p,dt);
          }
          if (!ok) { s.res.invalidReason=8; EG_FINISH(kTermInvalidField); }
          EG_FINISH(kTermOuterBoundaryAllowed);
        }
        if (outerDisposition==kSegMeshUnavailable || outerDisposition==kSegInvalid)
          { s.res.invalidReason=9; EG_FINISH(kTermInvalidField); }
      }
      else if (structured) {
        if (event.type==kEventInner) EG_FINISH(kTermInnerBoundaryForbidden);
        if (event.type==kEventOuter) {
          if (!PopulateExit(c,s,field,event.position,event.fraction,xBefore,x,pBefore,p,dt))
            { s.res.invalidReason=10; EG_FINISH(kTermInvalidField); }
          EG_FINISH(kTermOuterBoundaryAllowed);
        }
      }

      if (cx.trap.enabled) {
        if (mode3d) field.ResetAvailability();
        V Btrap; field.GetB(x,Btrap);
        if (mode3d && !field.Available()) { s.res.invalidReason=11; EG_FINISH(kTermInvalidField); }
        if (TrapUpdate(trap,field,x,p,Btrap,s.t)) {
          EG_FINISH(s.trap.mechanism==2 ? kTermDriftTrappedForbidden
                                        : kTermMagneticallyTrappedForbidden);
        }
      }
    }
#undef EG_FINISH
    s.x=x; s.p=p;
    return false;
  }

  // SI-only variant used for the Shue escape so the event position is the exact double
  // value computed by ClassifySegment (the legacy tracer copies it without conversion).
  template <class F>
  static EG_HD bool PopulateExitSIPosition(const Consts& c,EgSlot<R>& s,F& field,
                                           const double pos[3],R fraction,
                                           const V& xBefore,const V& xAfter,
                                           const V& pBefore,const V& pAfter,
                                           double acceptedStep) {
    if (!c.captureExit) return true;
    const double pb[3]={double(pBefore.x),double(pBefore.y),double(pBefore.z)};
    const double pa[3]={double(pAfter.x),double(pAfter.y),double(pAfter.z)};
    if (!PopulateExitKinematicsSI(s.res.exit,pos,pb,pa,double(fraction),s.t,acceptedStep,
                                  c.qabs_C))
      return false;
    const R a=fraction;
    const V chord=vsub(xAfter,xBefore);
    const R chordLength=vnorm(chord);
    R insetFraction=R(0.0);
    if (chordLength>R(0.0))
      insetFraction=emin(a,R(emax(1.0,c.boundaryTol_m))/chordLength);
    const V xEval=vadd(xBefore,vmul(emax(R(0.0),a-insetFraction),chord));
    field.ResetAvailability();
    V Bexit; field.GetB(xEval,Bexit);
    if (!field.Available()) return false;
    const double b[3]={double(Bexit.x),double(Bexit.y),double(Bexit.z)};
    return CompleteExitPitchAngleSI(s.res.exit,b);
  }
};

} // namespace GPU
} // namespace Geospace

#endif

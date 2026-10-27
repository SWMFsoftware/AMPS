#ifndef _SRC_GEOSPACE_GPU_ENGINE_TRAP_H_
#define _SRC_GEOSPACE_GPU_ENGINE_TRAP_H_

//======================================================================================
// EngineTrap.h -- fixed-capacity port of Geospace::TrajectoryTrap::Detector.
//
// The legacy detector keeps the drift-recurrence profiles in std::vector members that grow with
// the trajectory.  A device kernel cannot allocate, so the state is split into
//
//   EgTrapState   scalars (one per trajectory slot), and
//   an arena      per-slot doubles/ints whose size follows from the run configuration:
//                 nb = max(8,TRAP_DRIFT_PROFILE_BINS),
//                 window = max(1,TRAP_DRIFT_WINDOW_REVOLUTIONS),
//                 maxLag = max(max(1,TRAP_DRIFT_MAX_RECURRENCE_LAG),window),
//                 retained = maxLag+window completed profiles (the legacy detector erases older ones).
//
// The legacy detector's vector of completed profiles (push_back followed by erase(begin) above
// `retained`) is represented by a ring buffer with the same logical ordering, and
// CombineProfiles() is evaluated bin by bin with the identical summation order, so the
// detector is bitwise identical to the legacy tracer in double precision.  The detector always
// works in SI double precision regardless of the engine Real type.
//======================================================================================

#include "EngineTypes.h"

namespace Geospace {
namespace GPU {

struct EgDriftDiag {
  int available;
  int previousCoveragePassed;
  int currentCoveragePassed;
  int secularDriftPassed;
  int profilesRecurred;
  double previousCoverageFraction;
  double currentCoverageFraction;
  int commonBins;
  int matchedBins;
  int radialMismatchBins;
  int latitudeMismatchBins;
  int pitchMismatchBins;
  double matchedBinFraction;
  double maximumNormalizedRadialResidual;
  double maximumNormalizedLatitudeResidual;
  double maximumNormalizedPitchResidual;
  double windowMeanRadiusChange_m;
  double windowMeanRadiusSpan_m;
  int recurrenceLagRevolutions;
};

EG_HD void ClearDriftDiag(EgDriftDiag& d) {
  d.available=0; d.previousCoveragePassed=0; d.currentCoveragePassed=0;
  d.secularDriftPassed=0; d.profilesRecurred=0;
  d.previousCoverageFraction=0.0; d.currentCoverageFraction=0.0;
  d.commonBins=0; d.matchedBins=0; d.radialMismatchBins=0;
  d.latitudeMismatchBins=0; d.pitchMismatchBins=0;
  d.matchedBinFraction=0.0;
  d.maximumNormalizedRadialResidual=0.0;
  d.maximumNormalizedLatitudeResidual=0.0;
  d.maximumNormalizedPitchResidual=0.0;
  d.windowMeanRadiusChange_m=0.0; d.windowMeanRadiusSpan_m=0.0;
  d.recurrenceLagRevolutions=0;
}

struct EgTrapState {
  int initialized;
  int trapped;
  int mechanism;            // 0 none, 1 bounce, 2 drift

  int lastParallelSign;
  int mirrorPoints;
  int bounceCycles;
  int stableCycleComparisons;
  double lastBounceEnvelopeChange_m;
  int havePreviousCycle;
  double currentRMin;
  double currentRMax;
  double previousRMin;
  double previousRMax;

  double pReference;
  double pMin;
  double pMax;
  double minimumOuterMargin;

  int driftPhiInitialized;
  double lastDriftPhi;
  double driftAccumulatedAngle;
  int driftRevolutions;
  int stableDriftComparisons;
  int driftRecurrenceLagRevolutions;
  double lastDriftMeanRadiusChange_m;
  double lastDriftMeanRadiusSpan_m;
  EgDriftDiag lastDriftComparison;

  int ringHead;             // physical index of the oldest completed profile
  int ringCount;            // number of retained completed profiles
  int turnCount;            // recorded drift-turn states
  int turnOverflow;
};

// Geometry of the per-slot arena.
struct EgTrapLayout {
  int nb;
  int window;
  int maxLag;
  int retained;
  long long doublesPerSlot;   // (1+retained)*9*nb + retained
  long long intsPerSlot;      // (1+retained)*nb + maxLag+1
  int turnCapacity;
};

inline EgTrapLayout MakeTrapLayout(const EgTrapConfig& cfg,int turnCapacity) {
  EgTrapLayout l;
  l.nb=emax(8,cfg.driftProfileBins);
  l.window=emax(1,cfg.driftWindowRevolutions);
  l.maxLag=emax(emax(1,cfg.driftMaxRecurrenceLagRevolutions),l.window);
  l.retained=l.maxLag+l.window;
  const bool drift=(cfg.enabled!=0 && cfg.driftEnabled!=0);
  l.doublesPerSlot=drift ? (long long)(1+l.retained)*9LL*l.nb+l.retained : 0;
  l.intsPerSlot=drift ? (long long)(1+l.retained)*l.nb+(l.maxLag+1) : 0;
  l.turnCapacity=(drift && cfg.recordDriftTurnStates) ? turnCapacity : 0;
  return l;
}

// Profile quantity index inside a 9*nb double block.
enum : int { kPrRSum=0, kPrLatSum, kPrPitchSum, kPrRMin, kPrRMax, kPrLatMin, kPrLatMax,
             kPrPitchMin, kPrPitchMax, kPrCount };

class EgTrapDetector {
 public:
  EG_HD EgTrapDetector(const EgTrapConfig& cfg,const double* boxMin,const double* boxMax,
                       const EgTrapLayout& layout,EgTrapState& state,
                       double* dArena,int* iArena,EgTurnState* turns)
      : cfg_(cfg), boxMin_(boxMin), boxMax_(boxMax), l_(layout), s_(state),
        d_(dArena), i_(iArena), turns_(turns) {}

  // Equivalent of the legacy constructor (member initializers + currentProfile_.Reset).
  EG_HD void Reset() {
    const double inf=std::numeric_limits<double>::infinity();
    s_.initialized=0; s_.trapped=0; s_.mechanism=0;
    s_.lastParallelSign=0; s_.mirrorPoints=0; s_.bounceCycles=0;
    s_.stableCycleComparisons=0; s_.lastBounceEnvelopeChange_m=inf;
    s_.havePreviousCycle=0;
    s_.currentRMin=inf; s_.currentRMax=0.0; s_.previousRMin=0.0; s_.previousRMax=0.0;
    s_.pReference=0.0; s_.pMin=inf; s_.pMax=0.0; s_.minimumOuterMargin=inf;
    s_.driftPhiInitialized=0; s_.lastDriftPhi=0.0; s_.driftAccumulatedAngle=0.0;
    s_.driftRevolutions=0; s_.stableDriftComparisons=0;
    s_.driftRecurrenceLagRevolutions=0;
    s_.lastDriftMeanRadiusChange_m=0.0; s_.lastDriftMeanRadiusSpan_m=0.0;
    ClearDriftDiag(s_.lastDriftComparison);
    s_.ringHead=0; s_.ringCount=0; s_.turnCount=0; s_.turnOverflow=0;
    if (d_!=nullptr) {
      ResetProfile(CurD(),CurI());
      for (int lag=0;lag<=l_.maxLag;++lag) StableByLag()[lag]=0;
    }
  }

  EG_HD bool Update(const double x[3],const double p[3],const double B[3],
                    double traceTime_s) {
    if (!cfg_.enabled || s_.trapped) return s_.trapped!=0;
    if (!Finite3(x) || !Finite3(p) || !Finite3(B)) return false;

    const double pMag=Norm3(p);
    const double bMag=Norm3(B);
    const double r=Norm3(x);
    if (!(pMag>0.0) || !(bMag>0.0) || !(r>0.0) || !std::isfinite(r)) return false;

    const double pitchCos=(p[0]*B[0]+p[1]*B[1]+p[2]*B[2])/(pMag*bMag);
    const double pitchCos2=emax(0.0,emin(1.0,pitchCos*pitchCos));
    const double latitudeProxy=emax(-1.0,emin(1.0,x[2]/r));

    if (!s_.initialized) {
      s_.initialized=1;
      s_.pReference=pMag;
      s_.pMin=pMag;
      s_.pMax=pMag;
      ResetCurrentCycle(r,OuterMargin(x));
      InitializeDriftAzimuth(x);
    }
    else {
      s_.pMin=emin(s_.pMin,pMag);
      s_.pMax=emax(s_.pMax,pMag);
      s_.currentRMin=emin(s_.currentRMin,r);
      s_.currentRMax=emax(s_.currentRMax,r);
      s_.minimumOuterMargin=emin(s_.minimumOuterMargin,OuterMargin(x));
    }

    if (cfg_.driftEnabled) {
      UpdateDriftRecurrence(x,p,B,pMag,bMag,r,latitudeProxy,pitchCos2,traceTime_s);
      if (DriftCriteriaSatisfied()) {
        s_.trapped=1;
        s_.mechanism=2;
        return true;
      }
    }

    const int sign=(pitchCos>cfg_.parallelDeadband) ? 1 :
                   ((pitchCos< -cfg_.parallelDeadband) ? -1 : 0);
    if (sign==0) return s_.trapped!=0;
    if (s_.lastParallelSign==0) {
      s_.lastParallelSign=sign;
      return s_.trapped!=0;
    }
    if (sign==s_.lastParallelSign) return s_.trapped!=0;

    s_.lastParallelSign=sign;
    ++s_.mirrorPoints;

    if ((s_.mirrorPoints%2)!=0) return s_.trapped!=0;
    ++s_.bounceCycles;

    if (s_.havePreviousCycle) {
      const double envelopeChange=emax(
          std::fabs(s_.currentRMin-s_.previousRMin),
          std::fabs(s_.currentRMax-s_.previousRMax));
      s_.lastBounceEnvelopeChange_m=envelopeChange;
      if (envelopeChange<=cfg_.radialEnvelopeTolerance_m)
        ++s_.stableCycleComparisons;
      else
        s_.stableCycleComparisons=0;
    }

    s_.previousRMin=s_.currentRMin;
    s_.previousRMax=s_.currentRMax;
    s_.havePreviousCycle=1;

    const int requiredStableComparisons=2;
    s_.trapped=(s_.mirrorPoints>=cfg_.minMirrorPoints &&
                s_.bounceCycles>=cfg_.minBounceCycles &&
                s_.stableCycleComparisons>=requiredStableComparisons &&
                s_.minimumOuterMargin>=cfg_.outerMargin_m &&
                MomentumWithinTolerance()) ? 1 : 0;
    if (s_.trapped) s_.mechanism=1;

    ResetCurrentCycle(r,OuterMargin(x));
    return s_.trapped!=0;
  }

  EG_HD double momentumRelativeSpread() const {
    return (s_.pReference>0.0) ? (s_.pMax-s_.pMin)/s_.pReference : 0.0;
  }
  EG_HD double driftAngleRadians() const { return std::fabs(s_.driftAccumulatedAngle); }

 private:
  // ---- arena addressing -------------------------------------------------------------
  EG_HD double* CurD() const { return d_; }
  EG_HD int* CurI() const { return i_; }
  EG_HD double* RingD(int physical) const {
    return d_+(long long)(1+physical)*9LL*l_.nb;
  }
  EG_HD int* RingI(int physical) const { return i_+(long long)(1+physical)*l_.nb; }
  EG_HD double* MeanRadii() const { return d_+(long long)(1+l_.retained)*9LL*l_.nb; }
  EG_HD int* StableByLag() const { return i_+(long long)(1+l_.retained)*l_.nb; }
  EG_HD int Physical(int logical) const { return (s_.ringHead+logical)%l_.retained; }

  EG_HD void ResetProfile(double* pd,int* pi) const {
    const double inf=std::numeric_limits<double>::infinity();
    const int n=l_.nb;
    for (int b=0;b<n;++b) {
      pd[kPrRSum*n+b]=0.0; pd[kPrLatSum*n+b]=0.0; pd[kPrPitchSum*n+b]=0.0;
      pd[kPrRMin*n+b]=inf; pd[kPrRMax*n+b]=-inf;
      pd[kPrLatMin*n+b]=inf; pd[kPrLatMax*n+b]=-inf;
      pd[kPrPitchMin*n+b]=inf; pd[kPrPitchMax*n+b]=-inf;
      pi[b]=0;
    }
  }

  static EG_HD bool Finite3(const double a[3]) {
    return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]);
  }
  static EG_HD double Norm3(const double a[3]) {
    return std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  }

  EG_HD double OuterMargin(const double x[3]) const {
    return emin(
        emin(emin(x[0]-boxMin_[0],boxMax_[0]-x[0]),
             emin(x[1]-boxMin_[1],boxMax_[1]-x[1])),
        emin(x[2]-boxMin_[2],boxMax_[2]-x[2]));
  }

  EG_HD bool MomentumWithinTolerance() const {
    return momentumRelativeSpread()<=cfg_.energyRelativeTolerance;
  }

  EG_HD void ResetCurrentCycle(double r,double outerMargin) {
    s_.currentRMin=r;
    s_.currentRMax=r;
    s_.minimumOuterMargin=emin(s_.minimumOuterMargin,outerMargin);
  }

  EG_HD void InitializeDriftAzimuth(const double x[3]) {
    const double rho=std::hypot(x[0],x[1]);
    if (!(rho>0.0) || !std::isfinite(rho)) return;
    s_.lastDriftPhi=std::atan2(x[1],x[0]);
    s_.driftPhiInitialized=1;
  }

  EG_HD int CurrentProfileBin() const {
    const double kPi=3.141592653589793238462643383279502884;
    const double kTwoPi=2.0*kPi;
    const int n=l_.nb;
    if (n<=0) return 0;
    double phase=efmod(std::fabs(s_.driftAccumulatedAngle),kTwoPi)/kTwoPi;
    if (!std::isfinite(phase)) phase=0.0;
    int bin=static_cast<int>(std::floor(phase*n));
    if (bin<0) bin=0;
    if (bin>=n) bin=n-1;
    return bin;
  }

  EG_HD void AddCurrentDriftProfileSample(double r,double latitudeProxy,double pitchCos2) {
    const int i=CurrentProfileBin();
    const int n=l_.nb;
    double* pd=CurD();
    pd[kPrRSum*n+i]+=r;
    pd[kPrLatSum*n+i]+=latitudeProxy;
    pd[kPrPitchSum*n+i]+=pitchCos2;
    pd[kPrRMin*n+i]=emin(pd[kPrRMin*n+i],r);
    pd[kPrRMax*n+i]=emax(pd[kPrRMax*n+i],r);
    pd[kPrLatMin*n+i]=emin(pd[kPrLatMin*n+i],latitudeProxy);
    pd[kPrLatMax*n+i]=emax(pd[kPrLatMax*n+i],latitudeProxy);
    pd[kPrPitchMin*n+i]=emin(pd[kPrPitchMin*n+i],pitchCos2);
    pd[kPrPitchMax*n+i]=emax(pd[kPrPitchMax*n+i],pitchCos2);
    ++CurI()[i];
  }

  // Combined values of one bin over completed profiles [first,first+count) in the
  // order used by the legacy CombineProfiles(): start from Reset() values and fold.
  struct BinValues { double v[9]; int count; };

  EG_HD void CombinedBin(int first,int count,int bin,BinValues& out) const {
    const double inf=std::numeric_limits<double>::infinity();
    out.v[kPrRSum]=0.0; out.v[kPrLatSum]=0.0; out.v[kPrPitchSum]=0.0;
    out.v[kPrRMin]=inf; out.v[kPrRMax]=-inf;
    out.v[kPrLatMin]=inf; out.v[kPrLatMax]=-inf;
    out.v[kPrPitchMin]=inf; out.v[kPrPitchMax]=-inf;
    out.count=0;
    const int end=emin(s_.ringCount,first+count);
    const int n=l_.nb;
    for (int ip=first;ip<end;++ip) {
      const int ph=Physical(ip);
      const double* sd=RingD(ph);
      const int* si=RingI(ph);
      if (si[bin]<=0) continue;
      out.v[kPrRSum]+=sd[kPrRSum*n+bin];
      out.v[kPrLatSum]+=sd[kPrLatSum*n+bin];
      out.v[kPrPitchSum]+=sd[kPrPitchSum*n+bin];
      out.v[kPrRMin]=emin(out.v[kPrRMin],sd[kPrRMin*n+bin]);
      out.v[kPrRMax]=emax(out.v[kPrRMax],sd[kPrRMax*n+bin]);
      out.v[kPrLatMin]=emin(out.v[kPrLatMin],sd[kPrLatMin*n+bin]);
      out.v[kPrLatMax]=emax(out.v[kPrLatMax],sd[kPrLatMax*n+bin]);
      out.v[kPrPitchMin]=emin(out.v[kPrPitchMin],sd[kPrPitchMin*n+bin]);
      out.v[kPrPitchMax]=emax(out.v[kPrPitchMax],sd[kPrPitchMax*n+bin]);
      out.count+=si[bin];
    }
  }

  EG_HD int CombinedPopulatedBins(int first,int count) const {
    int populated=0;
    BinValues v;
    for (int b=0;b<l_.nb;++b) {
      CombinedBin(first,count,b,v);
      if (v.count>0) ++populated;
    }
    return populated;
  }

  EG_HD double CombinedCoverageFraction(int first,int count) const {
    const int n=l_.nb;
    return n>0 ? static_cast<double>(CombinedPopulatedBins(first,count))/
                 static_cast<double>(n) : 0.0;
  }

  EG_HD bool CombinedHasCoverage(int first,int count) const {
    if (l_.nb<=0) return false;
    return CombinedCoverageFraction(first,count)>=cfg_.driftMinProfileCoverage;
  }

  EG_HD double CombinedMeanRadius(int first,int count) const {
    double sum=0.0;
    int used=0;
    BinValues v;
    for (int b=0;b<l_.nb;++b) {
      CombinedBin(first,count,b,v);
      const int n=v.count;
      if (n<=0) continue;
      sum+=v.v[kPrRSum]/static_cast<double>(n);
      ++used;
    }
    return used>0 ? sum/static_cast<double>(used) :
                    std::numeric_limits<double>::quiet_NaN();
  }

  // ProfileMeanRadius(currentProfile_) for the profile being closed.
  EG_HD double CurrentMeanRadius() const {
    double sum=0.0;
    int used=0;
    const int n=l_.nb;
    const double* pd=CurD();
    const int* pi=CurI();
    for (int b=0;b<n;++b) {
      const int c=pi[b];
      if (c<=0) continue;
      sum+=pd[kPrRSum*n+b]/static_cast<double>(c);
      ++used;
    }
    return used>0 ? sum/static_cast<double>(used) :
                    std::numeric_limits<double>::quiet_NaN();
  }

  static EG_HD double NormalizedResidual(double residual,double tolerance) {
    if (tolerance>0.0) return residual/tolerance;
    return residual==0.0 ? 0.0 : std::numeric_limits<double>::infinity();
  }

  EG_HD bool ProfilesRecur(int previousFirst,int currentFirst,int window,
                           EgDriftDiag& diagnostics) const {
    diagnostics.available=1;
    diagnostics.previousCoverageFraction=CombinedCoverageFraction(previousFirst,window);
    diagnostics.currentCoverageFraction=CombinedCoverageFraction(currentFirst,window);
    diagnostics.previousCoveragePassed=CombinedHasCoverage(previousFirst,window) ? 1 : 0;
    diagnostics.currentCoveragePassed=CombinedHasCoverage(currentFirst,window) ? 1 : 0;
    if (!diagnostics.previousCoveragePassed || !diagnostics.currentCoveragePassed)
      return false;
    const int n=l_.nb;
    int common=0;
    int matched=0;
    BinValues prv,cur;
    for (int i=0;i<n;++i) {
      CombinedBin(previousFirst,window,i,prv);
      CombinedBin(currentFirst,window,i,cur);
      const int n0=prv.count;
      const int n1=cur.count;
      if (n0<=0 || n1<=0) continue;
      ++common;

      const double r0=prv.v[kPrRSum]/n0;
      const double r1=cur.v[kPrRSum]/n1;
      const double z0=prv.v[kPrLatSum]/n0;
      const double z1=cur.v[kPrLatSum]/n1;
      const double a0=prv.v[kPrPitchSum]/n0;
      const double a1=cur.v[kPrPitchSum]/n1;

      const double rTol=emax(cfg_.driftRadialAbsoluteTolerance_m,
                             cfg_.driftRadialRelativeTolerance*
                             emax(std::fabs(r0),std::fabs(r1)));
      const double radialMeanResidual=std::fabs(r1-r0);
      const double radialEnvelopeResidual=emax(
          std::fabs(cur.v[kPrRMin]-prv.v[kPrRMin]),
          std::fabs(cur.v[kPrRMax]-prv.v[kPrRMax]));
      const double latitudeMeanResidual=std::fabs(z1-z0);
      const double latitudeEnvelopeResidual=emax(
          std::fabs(cur.v[kPrLatMin]-prv.v[kPrLatMin]),
          std::fabs(cur.v[kPrLatMax]-prv.v[kPrLatMax]));
      const double pitchMeanResidual=std::fabs(a1-a0);
      const double pitchEnvelopeResidual=emax(
          std::fabs(cur.v[kPrPitchMin]-prv.v[kPrPitchMin]),
          std::fabs(cur.v[kPrPitchMax]-prv.v[kPrPitchMax]));

      const bool requireEnvelope=cfg_.driftWindowRevolutions>1;
      const double radialDecisionResidual=requireEnvelope
          ? emax(radialMeanResidual,radialEnvelopeResidual) : radialMeanResidual;
      const double latitudeDecisionResidual=requireEnvelope
          ? emax(latitudeMeanResidual,latitudeEnvelopeResidual) : latitudeMeanResidual;
      const double pitchDecisionResidual=requireEnvelope
          ? emax(pitchMeanResidual,pitchEnvelopeResidual) : pitchMeanResidual;
      const bool radialOk=radialDecisionResidual<=rTol;
      const bool latitudeOk=latitudeDecisionResidual<=cfg_.driftLatitudeTolerance;
      const bool pitchOk=pitchDecisionResidual<=cfg_.driftPitchCos2Tolerance;
      diagnostics.maximumNormalizedRadialResidual=emax(
          diagnostics.maximumNormalizedRadialResidual,
          NormalizedResidual(radialDecisionResidual,rTol));
      diagnostics.maximumNormalizedLatitudeResidual=emax(
          diagnostics.maximumNormalizedLatitudeResidual,
          NormalizedResidual(latitudeDecisionResidual,cfg_.driftLatitudeTolerance));
      diagnostics.maximumNormalizedPitchResidual=emax(
          diagnostics.maximumNormalizedPitchResidual,
          NormalizedResidual(pitchDecisionResidual,cfg_.driftPitchCos2Tolerance));
      if (!radialOk) ++diagnostics.radialMismatchBins;
      if (!latitudeOk) ++diagnostics.latitudeMismatchBins;
      if (!pitchOk) ++diagnostics.pitchMismatchBins;
      if (radialOk && latitudeOk && pitchOk) ++matched;
    }

    diagnostics.commonBins=common;
    diagnostics.matchedBins=matched;
    diagnostics.matchedBinFraction=common>0
        ? static_cast<double>(matched)/static_cast<double>(common) : 0.0;

    const int requiredCommon=static_cast<int>(std::ceil(
        cfg_.driftMinProfileCoverage*static_cast<double>(n)));
    if (common<requiredCommon) return false;
    return diagnostics.matchedBinFraction>=cfg_.driftMinMatchedBinFraction;
  }

  EG_HD void PushCompletedProfile(double meanRadius) {
    // legacy: push_back(currentProfile_) then erase(begin()) while size>retained.
    int slot;
    if (s_.ringCount<l_.retained) {
      slot=Physical(s_.ringCount);
      ++s_.ringCount;
    }
    else {
      slot=s_.ringHead;
      s_.ringHead=(s_.ringHead+1)%l_.retained;
    }
    const int n=l_.nb;
    double* dst=RingD(slot);
    int* dstI=RingI(slot);
    const double* src=CurD();
    const int* srcI=CurI();
    for (int k=0;k<9*n;++k) dst[k]=src[k];
    for (int b=0;b<n;++b) dstI[b]=srcI[b];
    MeanRadii()[slot]=meanRadius;
  }

  EG_HD double MeanRadiusLogical(int logical) const { return MeanRadii()[Physical(logical)]; }

  EG_HD void CompleteDriftRevolution(const double x[3],const double p[3],const double B[3],
                                     double pMag,double bMag,double r,
                                     double latitudeProxy,double pitchCos2,
                                     double traceTime_s) {
    ++s_.driftRevolutions;
    PushCompletedProfile(CurrentMeanRadius());

    if (cfg_.recordDriftTurnStates) {
      if (turns_!=nullptr && s_.turnCount<l_.turnCapacity) {
        EgTurnState& state=turns_[s_.turnCount];
        state.revolution=s_.driftRevolutions;
        state.traceTime_s=traceTime_s;
        state.accumulatedAngle_rad=s_.driftAccumulatedAngle;
        state.radius_m=r;
        state.latitude=latitudeProxy;
        state.pitchCos2=pitchCos2;
        state.outerMargin_m=OuterMargin(x);
        state.momentumRelativeSpread=momentumRelativeSpread();
        for (int d=0;d<3;++d) {
          state.x_m[d]=x[d];
          state.p_unit[d]=p[d]/pMag;
          state.b_unit[d]=B[d]/bMag;
        }
        ++s_.turnCount;
      }
      else s_.turnOverflow=1;
    }

    const int window=l_.window;
    const int maximumLag=l_.maxLag;

    ClearDriftDiag(s_.lastDriftComparison);
    s_.stableDriftComparisons=0;
    s_.driftRecurrenceLagRevolutions=0;

    int* stableByLag=StableByLag();
    bool haveDiagnostic=false;
    int bestStable=-1;
    double bestMatched=-1.0;
    for (int lag=window;lag<=maximumLag;++lag) {
      if (s_.ringCount<lag+window) {
        stableByLag[lag]=0;
        continue;
      }

      const int currentFirst=s_.ringCount-window;
      const int previousFirst=currentFirst-lag;
      const double previousMean=CombinedMeanRadius(previousFirst,window);
      const double currentMean=CombinedMeanRadius(currentFirst,window);
      double meanChange=std::numeric_limits<double>::infinity();
      if (std::isfinite(previousMean) && std::isfinite(currentMean))
        meanChange=std::fabs(currentMean-previousMean);

      double minimumMean=std::numeric_limits<double>::infinity();
      double maximumMean=-std::numeric_limits<double>::infinity();
      bool finiteMeans=true;
      for (int i=previousFirst;i<s_.ringCount;++i) {
        const double value=MeanRadiusLogical(i);
        if (!std::isfinite(value)) {
          finiteMeans=false;
          break;
        }
        minimumMean=emin(minimumMean,value);
        maximumMean=emax(maximumMean,value);
      }
      const double meanSpan=finiteMeans ? maximumMean-minimumMean :
          std::numeric_limits<double>::infinity();

      const bool secularDriftOk=(cfg_.driftMaxMeanRadiusChange_m<=0.0) ||
          (meanChange<=cfg_.driftMaxMeanRadiusChange_m &&
           meanSpan<=cfg_.driftMaxMeanRadiusChange_m);

      EgDriftDiag candidate;
      ClearDriftDiag(candidate);
      candidate.recurrenceLagRevolutions=lag;
      candidate.secularDriftPassed=secularDriftOk ? 1 : 0;
      candidate.windowMeanRadiusChange_m=meanChange;
      candidate.windowMeanRadiusSpan_m=meanSpan;
      const bool recurrent=ProfilesRecur(previousFirst,currentFirst,window,candidate);
      candidate.profilesRecurred=recurrent ? 1 : 0;
      if (recurrent && secularDriftOk)
        ++stableByLag[lag];
      else
        stableByLag[lag]=0;

      const int stable=stableByLag[lag];
      s_.stableDriftComparisons=emax(s_.stableDriftComparisons,stable);

      if (!haveDiagnostic || stable>bestStable ||
          (stable==bestStable && candidate.matchedBinFraction>bestMatched)) {
        haveDiagnostic=true;
        bestStable=stable;
        bestMatched=candidate.matchedBinFraction;
        s_.lastDriftComparison=candidate;
        s_.lastDriftMeanRadiusChange_m=meanChange;
        s_.lastDriftMeanRadiusSpan_m=meanSpan;
        s_.driftRecurrenceLagRevolutions=lag;
      }
    }
    ResetProfile(CurD(),CurI());
  }

  EG_HD void UpdateDriftRecurrence(const double x[3],const double p[3],const double B[3],
                                   double pMag,double bMag,double r,
                                   double latitudeProxy,double pitchCos2,
                                   double traceTime_s) {
    const double rho=std::hypot(x[0],x[1]);
    if (!(rho>0.0) || !std::isfinite(rho)) return;

    const double phi=std::atan2(x[1],x[0]);
    if (!s_.driftPhiInitialized) {
      s_.lastDriftPhi=phi;
      s_.driftPhiInitialized=1;
      AddCurrentDriftProfileSample(r,latitudeProxy,pitchCos2);
      return;
    }

    double dphi=phi-s_.lastDriftPhi;
    const double kPi=3.141592653589793238462643383279502884;
    const double kTwoPi=2.0*kPi;
    while (dphi>kPi) dphi-=kTwoPi;
    while (dphi< -kPi) dphi+=kTwoPi;
    if (std::isfinite(dphi)) s_.driftAccumulatedAngle+=dphi;
    s_.lastDriftPhi=phi;

    const int completed=static_cast<int>(
        std::floor(std::fabs(s_.driftAccumulatedAngle)/kTwoPi));
    while (s_.driftRevolutions<completed)
      CompleteDriftRevolution(x,p,B,pMag,bMag,r,latitudeProxy,pitchCos2,traceTime_s);

    AddCurrentDriftProfileSample(r,latitudeProxy,pitchCos2);
  }

  EG_HD bool DriftCriteriaSatisfied() const {
    const int window=emax(1,cfg_.driftWindowRevolutions);
    const int requiredComparisons=(window==1)
        ? emax(1,cfg_.minDriftRevolutions-1)
        : emax(1,cfg_.driftMinStableWindowComparisons);
    const int selectedLag=emax(window,s_.driftRecurrenceLagRevolutions);
    const int requiredRevolutions=emax(
        cfg_.minDriftRevolutions,
        selectedLag+window+requiredComparisons-1);
    if (!cfg_.driftEnabled || s_.driftRevolutions<requiredRevolutions)
      return false;
    return s_.stableDriftComparisons>=requiredComparisons &&
           s_.minimumOuterMargin>=cfg_.outerMargin_m &&
           MomentumWithinTolerance();
  }

  const EgTrapConfig& cfg_;
  const double* boxMin_;
  const double* boxMax_;
  const EgTrapLayout& l_;
  EgTrapState& s_;
  double* d_;
  int* i_;
  EgTurnState* turns_;
};

} // namespace GPU
} // namespace Geospace

#endif

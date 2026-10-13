#ifndef _SRC_EARTH_UTIL_TRAJECTORY_TRAP_DETECTOR_H_
#define _SRC_EARTH_UTIL_TRAJECTORY_TRAP_DETECTOR_H_

#include "TrajectoryBoundary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace Earth {
namespace TrajectoryTrap {

// Conservative trapped-orbit detector for one FROZEN magnetic-field snapshot.
//
// C19 needs a positive physical FORBIDDEN verdict for long-lived GEO trajectories;
// simply changing TIME_LIMIT/DISTANCE_LIMIT into FORBIDDEN would turn a numerical
// budget into physics.  This class therefore contains two independent recurrence
// tests which may terminate a trace only after the orbit itself supplies evidence of
// bounded motion:
//
//   (1) BOUNCE recurrence (legacy/F3 path): repeated parallel-velocity reversals,
//       several complete bounce cycles, and a stable radial envelope.
//
//   (2) DRIFT recurrence (C19 path): repeated complete azimuthal circuits plus a
//       repeatable *azimuth-resolved phase-space profile*.  The profile compares
//       gyro-averaged geocentric radius, normalized z/r, and cos^2(pitch angle) in
//       fixed drift-phase bins from one complete revolution to the next.  This is
//       deliberately stronger than the first C19 implementation, which compared only
//       global r_min/r_max envelopes and could confuse two geometrically different
//       shells that happened to have similar extrema.
//
// The drift branch continues to integrate the full Lorentz orbit.  It does NOT use a
// guiding-centre approximation to classify 15--190 MeV C19 protons: their gyroradii
// at GEO can be an appreciable fraction of the magnetospheric scale length, precisely
// where the penumbra is non-adiabatic.  Averaging only the recurrence observables in
// azimuth bins removes most gyrophase sensitivity while preserving the authoritative
// full-orbit trajectory.
//
// Neither branch ever infers trapping from elapsed time/path length.  If recurrence is
// not demonstrated before a configured numerical budget is exhausted, the caller must
// keep the trajectory UNRESOLVED.  The detector is valid only when the field is frozen
// throughout one backtrace.
struct Config {
  bool enabled{false};
  int minMirrorPoints{8};
  int minBounceCycles{4};
  double outerMargin_m{0.0};
  double radialEnvelopeTolerance_m{0.0};
  double energyRelativeTolerance{1.0e-4};
  double parallelDeadband{1.0e-6};

  // Optional C19 drift-shell recurrence branch.
  bool driftEnabled{false};
  int minDriftRevolutions{3};

  // Absolute + relative radius tolerances are applied per azimuth-profile bin:
  // |r_n-r_{n-1}| <= max(absTol, relTol*max(r_n,r_{n-1})).
  // The relative term is essential in non-axisymmetric T05 snapshots, where a closed
  // shell can have a large day/night excursion yet repeat from one drift turn to the
  // next.  The absolute term protects small-r bins from an unrealistically tiny bound.
  double driftRadialAbsoluteTolerance_m{0.0};
  double driftRadialRelativeTolerance{0.20};

  // Additional full-orbit recurrence observables.  z/r is dimensionless and insensitive
  // to a uniform radial breathing of the shell.  cos^2(alpha) removes the sign flip at
  // mirror points while retaining pitch-angle structure.
  double driftLatitudeTolerance{0.20};
  double driftPitchCos2Tolerance{0.25};

  // Optional secular-drift veto.  A recurrent shell should not only look similar in
  // each azimuth bin; its mean radial profile should also stop migrating systematically
  // from one completed revolution to the next.  A value <=0 disables this extra gate
  // for backward compatibility.  C19 enables a finite value as an additional safety
  // check against classifying a slowly leaking/quasi-trapped orbit as permanently trapped.
  double driftMaxMeanRadiusChange_m{0.0};

  // A revolution must sample enough azimuth bins to be meaningful, and a configurable
  // fraction of the bins present in both turns must satisfy all recurrence tolerances.
  // These gates prevent a sparse or accidental single-point return from being called a
  // closed drift shell.
  int driftProfileBins{24};
  double driftMinProfileCoverage{0.70};
  double driftMinMatchedBinFraction{0.75};

  // Number of complete drift revolutions combined into one phase-mixed profile.
  // The historical detector compares one revolution with the next, so the default
  // remains one and is bit-for-bit compatible with existing input decks.  Values
  // greater than one are intended for bounded quasi-periodic full orbits whose gyro,
  // bounce, and drift periods are incommensurate: adjacent single turns can start at
  // different bounce/gyro phases even though their multi-turn phase-space support is
  // recurrent.  Two compared windows are always non-overlapping.
  int driftWindowRevolutions{1};

  // Window mode must succeed repeatedly before it may emit a physical trapped state.
  // A value of two means that, for a three-turn window, at least seven completed turns
  // are required: turns [1..3] versus [4..6], followed by [2..4] versus [5..7].
  // This is deliberately stronger than accepting one fortuitous multi-turn match.
  int driftMinStableWindowComparisons{2};

  // Largest separation, in complete drift revolutions, at which the detector may
  // look for a recurring profile.  A value of one preserves the historical
  // adjacent-turn comparison.  Values greater than one do not loosen any profile
  // tolerance: they allow a quasi-periodic full orbit to demonstrate its actual
  // bounce/gyro--drift super-period.  The same lag must pass on consecutive turns,
  // and the secular-radius veto covers the entire interval between the profiles.
  int driftMaxRecurrenceLagRevolutions{1};

  // Forensic recording is deliberately independent of the classification policy.
  // Production tests leave this false.  A focused diagnostic may retain one compact
  // state at each completed drift turn to determine whether a true Poincare return
  // exists before any new recurrence observable is admitted into the detector.
  bool recordDriftTurnStates{false};
};

struct DriftTurnState {
  int revolution{0};
  double traceTime_s{0.0};
  double accumulatedAngle_rad{0.0};
  double x_m[3]{0.0,0.0,0.0};
  double p_unit[3]{0.0,0.0,0.0};
  double b_unit[3]{0.0,0.0,0.0};
  double radius_m{0.0};
  double latitude{0.0};
  double pitchCos2{0.0};
  double outerMargin_m{0.0};
  double momentumRelativeSpread{0.0};
};

// Complete audit of the most recent drift-profile comparison.  These values are
// propagated into trajectory diagnostics so a TIME_LIMIT can be attributed to a
// specific failed physical gate instead of being explained by guesswork.  Fractions
// are in [0,1]; normalized residuals are dimensionless ratios to their configured
// tolerance (zero when both residual and tolerance are zero).
struct DriftComparisonDiagnostics {
  bool available{false};
  bool previousCoveragePassed{false};
  bool currentCoveragePassed{false};
  bool secularDriftPassed{false};
  bool profilesRecurred{false};
  double previousCoverageFraction{0.0};
  double currentCoverageFraction{0.0};
  int commonBins{0};
  int matchedBins{0};
  int radialMismatchBins{0};
  int latitudeMismatchBins{0};
  int pitchMismatchBins{0};
  double matchedBinFraction{0.0};
  double maximumNormalizedRadialResidual{0.0};
  double maximumNormalizedLatitudeResidual{0.0};
  double maximumNormalizedPitchResidual{0.0};
  double windowMeanRadiusChange_m{0.0};
  double windowMeanRadiusSpan_m{0.0};
  int recurrenceLagRevolutions{0};
};

enum class Mechanism {
  None=0,
  Bounce=1,
  Drift=2
};

class Detector {
 public:
  Detector(const Config& config, const Earth::TrajectoryBoundary::Box& box)
      : cfg_(config), box_(box) {
    const int n=std::max(8,cfg_.driftProfileBins);
    currentProfile_.Reset(n);
  }

  bool Update(const double x[3], const double p[3], const double B[3],
              double traceTime_s=0.0) {
    if (!cfg_.enabled || trapped_) return trapped_;
    if (!Finite3(x) || !Finite3(p) || !Finite3(B)) return false;

    const double pMag=Norm3(p);
    const double bMag=Norm3(B);
    const double r=Norm3(x);
    if (!(pMag>0.0) || !(bMag>0.0) || !(r>0.0) || !std::isfinite(r)) return false;

    const double pitchCos=(p[0]*B[0]+p[1]*B[1]+p[2]*B[2])/(pMag*bMag);
    const double pitchCos2=std::max(0.0,std::min(1.0,pitchCos*pitchCos));
    const double latitudeProxy=std::max(-1.0,std::min(1.0,x[2]/r));

    if (!initialized_) {
      initialized_=true;
      pReference_=pMag;
      pMin_=pMag;
      pMax_=pMag;
      ResetCurrentCycle(r,OuterMargin(x));
      InitializeDriftAzimuth(x);
    }
    else {
      pMin_=std::min(pMin_,pMag);
      pMax_=std::max(pMax_,pMag);
      currentRMin_=std::min(currentRMin_,r);
      currentRMax_=std::max(currentRMax_,r);
      minimumOuterMargin_=std::min(minimumOuterMargin_,OuterMargin(x));
    }

    // Update the drift recurrence before the bounce branch.  Near 90-degree pitch
    // angle, v_parallel can hover around zero and the legacy mirror-point detector may
    // see no clean sign reversals even though the full orbit repeatedly drifts around
    // Earth.  The drift branch therefore needs every spatial sample.
    if (cfg_.driftEnabled) {
      UpdateDriftRecurrence(x,p,B,pMag,bMag,r,latitudeProxy,pitchCos2,traceTime_s);
      if (DriftCriteriaSatisfied()) {
        trapped_=true;
        mechanism_=Mechanism::Drift;
        return true;
      }
    }

    const int sign=(pitchCos>cfg_.parallelDeadband) ? 1 :
                   ((pitchCos<-cfg_.parallelDeadband) ? -1 : 0);
    if (sign==0) return trapped_;
    if (lastParallelSign_==0) {
      lastParallelSign_=sign;
      return trapped_;
    }
    if (sign==lastParallelSign_) return trapped_;

    lastParallelSign_=sign;
    ++mirrorPoints_;

    // Two successive mirror-point reversals form one complete bounce cycle.
    if ((mirrorPoints_%2)!=0) return trapped_;
    ++bounceCycles_;

    if (havePreviousCycle_) {
      const double envelopeChange=std::max(
          std::fabs(currentRMin_-previousRMin_),
          std::fabs(currentRMax_-previousRMax_));
      lastBounceEnvelopeChange_m_=envelopeChange;
      if (envelopeChange<=cfg_.radialEnvelopeTolerance_m)
        ++stableCycleComparisons_;
      else
        stableCycleComparisons_=0;
    }

    previousRMin_=currentRMin_;
    previousRMax_=currentRMax_;
    havePreviousCycle_=true;

    // Require stability across the last two cycle-to-cycle comparisons.  Together
    // with the default four-cycle minimum, this prevents a single fortuitously similar
    // pair of bounce envelopes from being classified as trapped.
    const int requiredStableComparisons=2;
    trapped_=(mirrorPoints_>=cfg_.minMirrorPoints &&
              bounceCycles_>=cfg_.minBounceCycles &&
              stableCycleComparisons_>=requiredStableComparisons &&
              minimumOuterMargin_>=cfg_.outerMargin_m &&
              MomentumWithinTolerance());
    if (trapped_) mechanism_=Mechanism::Bounce;

    ResetCurrentCycle(r,OuterMargin(x));
    return trapped_;
  }

  bool trapped() const { return trapped_; }
  Mechanism mechanism() const { return mechanism_; }
  int mirrorPoints() const { return mirrorPoints_; }
  int bounceCycles() const { return bounceCycles_; }
  int driftRevolutions() const { return driftRevolutions_; }
  int stableDriftComparisons() const { return stableDriftComparisons_; }
  int driftRecurrenceLagRevolutions() const { return driftRecurrenceLagRevolutions_; }
  double driftAngleRadians() const { return std::fabs(driftAccumulatedAngle_); }
  double driftMeanRadiusChangeMeters() const { return lastDriftMeanRadiusChange_m_; }
  double driftMeanRadiusSpanMeters() const { return lastDriftMeanRadiusSpan_m_; }
  double minimumOuterMarginMeters() const { return minimumOuterMargin_; }
  double lastBounceEnvelopeChangeMeters() const { return lastBounceEnvelopeChange_m_; }
  const DriftComparisonDiagnostics& driftComparisonDiagnostics() const {
    return lastDriftComparison_;
  }
  const std::vector<DriftTurnState>& driftTurnStates() const {
    return driftTurnStates_;
  }
  double momentumRelativeSpread() const {
    return (pReference_>0.0) ? (pMax_-pMin_)/pReference_ : 0.0;
  }

 private:
  struct Profile {
    std::vector<double> rSum;
    std::vector<double> latitudeSum;
    std::vector<double> pitchCos2Sum;
    // Extrema prevent a multi-turn average from concealing a large phase-space
    // excursion.  A bin matches only when both its means and its support envelopes
    // recur within the same configured physical tolerances.
    std::vector<double> rMin,rMax;
    std::vector<double> latitudeMin,latitudeMax;
    std::vector<double> pitchCos2Min,pitchCos2Max;
    std::vector<int> count;

    void Reset(int n) {
      rSum.assign((std::size_t)n,0.0);
      latitudeSum.assign((std::size_t)n,0.0);
      pitchCos2Sum.assign((std::size_t)n,0.0);
      rMin.assign((std::size_t)n,std::numeric_limits<double>::infinity());
      rMax.assign((std::size_t)n,-std::numeric_limits<double>::infinity());
      latitudeMin.assign((std::size_t)n,std::numeric_limits<double>::infinity());
      latitudeMax.assign((std::size_t)n,-std::numeric_limits<double>::infinity());
      pitchCos2Min.assign((std::size_t)n,std::numeric_limits<double>::infinity());
      pitchCos2Max.assign((std::size_t)n,-std::numeric_limits<double>::infinity());
      count.assign((std::size_t)n,0);
    }
    int bins() const { return static_cast<int>(count.size()); }
    int populatedBins() const {
      return static_cast<int>(std::count_if(count.begin(),count.end(),
                                            [](int n){return n>0;}));
    }
  };

  static bool Finite3(const double a[3]) {
    return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]);
  }

  static double Norm3(const double a[3]) {
    return std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  }

  double OuterMargin(const double x[3]) const {
    return std::min(
        std::min(std::min(x[0]-box_.min[0],box_.max[0]-x[0]),
                 std::min(x[1]-box_.min[1],box_.max[1]-x[1])),
        std::min(x[2]-box_.min[2],box_.max[2]-x[2]));
  }

  bool MomentumWithinTolerance() const {
    return momentumRelativeSpread()<=cfg_.energyRelativeTolerance;
  }

  void ResetCurrentCycle(double r, double outerMargin) {
    currentRMin_=r;
    currentRMax_=r;
    minimumOuterMargin_=std::min(minimumOuterMargin_,outerMargin);
  }

  void InitializeDriftAzimuth(const double x[3]) {
    const double rho=std::hypot(x[0],x[1]);
    if (!(rho>0.0) || !std::isfinite(rho)) return;
    lastDriftPhi_=std::atan2(x[1],x[0]);
    driftPhiInitialized_=true;
  }

  int CurrentProfileBin() const {
    constexpr double kPi=3.141592653589793238462643383279502884;
    constexpr double kTwoPi=2.0*kPi;
    const int n=currentProfile_.bins();
    if (n<=0) return 0;
    double phase=std::fmod(std::fabs(driftAccumulatedAngle_),kTwoPi)/kTwoPi;
    if (!std::isfinite(phase)) phase=0.0;
    int bin=static_cast<int>(std::floor(phase*n));
    if (bin<0) bin=0;
    if (bin>=n) bin=n-1;
    return bin;
  }

  void AddCurrentDriftProfileSample(double r,double latitudeProxy,double pitchCos2) {
    const int bin=CurrentProfileBin();
    const std::size_t i=static_cast<std::size_t>(bin);
    currentProfile_.rSum[i]+=r;
    currentProfile_.latitudeSum[i]+=latitudeProxy;
    currentProfile_.pitchCos2Sum[i]+=pitchCos2;
    currentProfile_.rMin[i]=std::min(currentProfile_.rMin[i],r);
    currentProfile_.rMax[i]=std::max(currentProfile_.rMax[i],r);
    currentProfile_.latitudeMin[i]=std::min(currentProfile_.latitudeMin[i],latitudeProxy);
    currentProfile_.latitudeMax[i]=std::max(currentProfile_.latitudeMax[i],latitudeProxy);
    currentProfile_.pitchCos2Min[i]=std::min(currentProfile_.pitchCos2Min[i],pitchCos2);
    currentProfile_.pitchCos2Max[i]=std::max(currentProfile_.pitchCos2Max[i],pitchCos2);
    ++currentProfile_.count[i];
  }

  double ProfileCoverageFraction(const Profile& profile) const {
    const int n=profile.bins();
    return n>0 ? static_cast<double>(profile.populatedBins())/static_cast<double>(n) : 0.0;
  }

  bool ProfileHasCoverage(const Profile& profile) const {
    const int n=profile.bins();
    if (n<=0) return false;
    return ProfileCoverageFraction(profile)>=cfg_.driftMinProfileCoverage;
  }

  double ProfileMeanRadius(const Profile& profile) const {
    // Average the per-bin mean radii rather than all raw samples.  Each azimuth bin
    // therefore contributes equally even when adaptive timesteps place more full-orbit
    // samples in one sector.  This is a diagnostic/gating observable only; the mover
    // itself still integrates the complete Lorentz orbit without guiding-centre averaging.
    double sum=0.0;
    int used=0;
    for (int i=0;i<profile.bins();++i) {
      const int n=profile.count[(std::size_t)i];
      if (n<=0) continue;
      sum += profile.rSum[(std::size_t)i]/static_cast<double>(n);
      ++used;
    }
    return used>0 ? sum/static_cast<double>(used) :
                    std::numeric_limits<double>::quiet_NaN();
  }

  static double NormalizedResidual(double residual,double tolerance) {
    if (tolerance>0.0) return residual/tolerance;
    return residual==0.0 ? 0.0 : std::numeric_limits<double>::infinity();
  }

  // Sum several complete revolution profiles without reweighting individual samples.
  // The resulting per-bin means are phase-mixed over the requested number of turns,
  // while the extrema retain the complete support envelope across those turns.
  Profile CombineProfiles(std::size_t first,std::size_t count) const {
    Profile combined;
    combined.Reset(std::max(8,cfg_.driftProfileBins));
    const std::size_t end=std::min(completedProfiles_.size(),first+count);
    for (std::size_t ip=first;ip<end;++ip) {
      const Profile& source=completedProfiles_[ip];
      for (int bin=0;bin<combined.bins();++bin) {
        const std::size_t i=static_cast<std::size_t>(bin);
        if (source.count[i]<=0) continue;
        combined.rSum[i]+=source.rSum[i];
        combined.latitudeSum[i]+=source.latitudeSum[i];
        combined.pitchCos2Sum[i]+=source.pitchCos2Sum[i];
        combined.rMin[i]=std::min(combined.rMin[i],source.rMin[i]);
        combined.rMax[i]=std::max(combined.rMax[i],source.rMax[i]);
        combined.latitudeMin[i]=std::min(combined.latitudeMin[i],source.latitudeMin[i]);
        combined.latitudeMax[i]=std::max(combined.latitudeMax[i],source.latitudeMax[i]);
        combined.pitchCos2Min[i]=std::min(combined.pitchCos2Min[i],source.pitchCos2Min[i]);
        combined.pitchCos2Max[i]=std::max(combined.pitchCos2Max[i],source.pitchCos2Max[i]);
        combined.count[i]+=source.count[i];
      }
    }
    return combined;
  }

  bool ProfilesRecur(const Profile& previous,const Profile& current,
                      DriftComparisonDiagnostics& diagnostics) const {
    diagnostics.available=true;
    diagnostics.previousCoverageFraction=ProfileCoverageFraction(previous);
    diagnostics.currentCoverageFraction=ProfileCoverageFraction(current);
    diagnostics.previousCoveragePassed=ProfileHasCoverage(previous);
    diagnostics.currentCoveragePassed=ProfileHasCoverage(current);
    if (!diagnostics.previousCoveragePassed || !diagnostics.currentCoveragePassed)
      return false;
    const int n=std::min(previous.bins(),current.bins());
    int common=0;
    int matched=0;
    for (int i=0;i<n;++i) {
      const int n0=previous.count[(std::size_t)i];
      const int n1=current.count[(std::size_t)i];
      if (n0<=0 || n1<=0) continue;
      ++common;

      const double r0=previous.rSum[(std::size_t)i]/n0;
      const double r1=current.rSum[(std::size_t)i]/n1;
      const double z0=previous.latitudeSum[(std::size_t)i]/n0;
      const double z1=current.latitudeSum[(std::size_t)i]/n1;
      const double a0=previous.pitchCos2Sum[(std::size_t)i]/n0;
      const double a1=current.pitchCos2Sum[(std::size_t)i]/n1;

      const double rTol=std::max(cfg_.driftRadialAbsoluteTolerance_m,
                                 cfg_.driftRadialRelativeTolerance*
                                 std::max(std::fabs(r0),std::fabs(r1)));
      const double radialMeanResidual=std::fabs(r1-r0);
      const double radialEnvelopeResidual=std::max(
          std::fabs(current.rMin[(std::size_t)i]-previous.rMin[(std::size_t)i]),
          std::fabs(current.rMax[(std::size_t)i]-previous.rMax[(std::size_t)i]));
      const double latitudeMeanResidual=std::fabs(z1-z0);
      const double latitudeEnvelopeResidual=std::max(
          std::fabs(current.latitudeMin[(std::size_t)i]-previous.latitudeMin[(std::size_t)i]),
          std::fabs(current.latitudeMax[(std::size_t)i]-previous.latitudeMax[(std::size_t)i]));
      const double pitchMeanResidual=std::fabs(a1-a0);
      const double pitchEnvelopeResidual=std::max(
          std::fabs(current.pitchCos2Min[(std::size_t)i]-previous.pitchCos2Min[(std::size_t)i]),
          std::fabs(current.pitchCos2Max[(std::size_t)i]-previous.pitchCos2Max[(std::size_t)i]));

      // Means establish recurrence of the central profile.  Multi-turn mode also
      // requires the complete within-bin envelopes to recur, which prevents phase
      // mixing from hiding a large excursion.  The legacy one-turn mode intentionally
      // retains its historical mean-only decision so existing decks do not change.
      const bool requireEnvelope=cfg_.driftWindowRevolutions>1;
      const double radialDecisionResidual=requireEnvelope
          ? std::max(radialMeanResidual,radialEnvelopeResidual) : radialMeanResidual;
      const double latitudeDecisionResidual=requireEnvelope
          ? std::max(latitudeMeanResidual,latitudeEnvelopeResidual) : latitudeMeanResidual;
      const double pitchDecisionResidual=requireEnvelope
          ? std::max(pitchMeanResidual,pitchEnvelopeResidual) : pitchMeanResidual;
      const bool radialOk=radialDecisionResidual<=rTol;
      const bool latitudeOk=latitudeDecisionResidual<=cfg_.driftLatitudeTolerance;
      const bool pitchOk=pitchDecisionResidual<=cfg_.driftPitchCos2Tolerance;
      diagnostics.maximumNormalizedRadialResidual=std::max(
          diagnostics.maximumNormalizedRadialResidual,
          NormalizedResidual(radialDecisionResidual,rTol));
      diagnostics.maximumNormalizedLatitudeResidual=std::max(
          diagnostics.maximumNormalizedLatitudeResidual,
          NormalizedResidual(latitudeDecisionResidual,
                             cfg_.driftLatitudeTolerance));
      diagnostics.maximumNormalizedPitchResidual=std::max(
          diagnostics.maximumNormalizedPitchResidual,
          NormalizedResidual(pitchDecisionResidual,
                             cfg_.driftPitchCos2Tolerance));
      if (!radialOk) ++diagnostics.radialMismatchBins;
      if (!latitudeOk) ++diagnostics.latitudeMismatchBins;
      if (!pitchOk) ++diagnostics.pitchMismatchBins;
      if (radialOk && latitudeOk && pitchOk) ++matched;
    }

    diagnostics.commonBins=common;
    diagnostics.matchedBins=matched;
    diagnostics.matchedBinFraction=common>0
        ? static_cast<double>(matched)/static_cast<double>(common) : 0.0;

    // Require enough common bins independently of each turn's own coverage.  This
    // protects against two well-sampled revolutions whose missing bins are disjoint.
    const int requiredCommon=static_cast<int>(std::ceil(
        cfg_.driftMinProfileCoverage*static_cast<double>(n)));
    if (common<requiredCommon) return false;
    return diagnostics.matchedBinFraction>=cfg_.driftMinMatchedBinFraction;
  }

  void CompleteDriftRevolution(const double x[3],const double p[3],const double B[3],
                               double pMag,double bMag,double r,
                               double latitudeProxy,double pitchCos2,
                               double traceTime_s) {
    ++driftRevolutions_;
    completedProfiles_.push_back(currentProfile_);
    completedMeanRadii_.push_back(ProfileMeanRadius(currentProfile_));

    if (cfg_.recordDriftTurnStates) {
      DriftTurnState state;
      state.revolution=driftRevolutions_;
      state.traceTime_s=traceTime_s;
      state.accumulatedAngle_rad=driftAccumulatedAngle_;
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
      driftTurnStates_.push_back(state);
    }

    const std::size_t window=static_cast<std::size_t>(
        std::max(1,cfg_.driftWindowRevolutions));
    // Retain enough complete turns to compare the newest window with every allowed
    // earlier non-overlapping window.  This is bounded by the parser (64 turns), so
    // long production traces cannot grow detector memory without limit.
    const std::size_t maximumLag=static_cast<std::size_t>(
        std::max(std::max(1,cfg_.driftMaxRecurrenceLagRevolutions),
                 static_cast<int>(window)));
    const std::size_t retained=maximumLag+window;
    if (completedProfiles_.size()>retained) {
      completedProfiles_.erase(completedProfiles_.begin());
      completedMeanRadii_.erase(completedMeanRadii_.begin());
    }

    lastDriftComparison_=DriftComparisonDiagnostics{};
    stableDriftComparisons_=0;
    driftRecurrenceLagRevolutions_=0;

    // Each candidate lag has an independent consecutive-success counter.  Mixing
    // successful comparisons from different lags would manufacture recurrence that
    // the trajectory never exhibited, so a counter is reset whenever that exact lag
    // fails either the profile or secular gate.
    if (stableComparisonsByLag_.size()<=maximumLag)
      stableComparisonsByLag_.resize(maximumLag+1,0);

    bool haveDiagnostic=false;
    int bestStable=-1;
    double bestMatched=-1.0;
    for (std::size_t lag=window;lag<=maximumLag;++lag) {
      if (completedProfiles_.size()<lag+window) {
        stableComparisonsByLag_[lag]=0;
        continue;
      }

      const std::size_t currentFirst=completedProfiles_.size()-window;
      const std::size_t previousFirst=currentFirst-lag;
      const Profile previous=CombineProfiles(previousFirst,window);
      const Profile current=CombineProfiles(currentFirst,window);
      const double previousMean=ProfileMeanRadius(previous);
      const double currentMean=ProfileMeanRadius(current);
      double meanChange=std::numeric_limits<double>::infinity();
      if (std::isfinite(previousMean) && std::isfinite(currentMean))
        meanChange=std::fabs(currentMean-previousMean);

      // The span is evaluated across every individual turn from the beginning of
      // the earlier window through the end of the current window.  Thus a long-lag
      // return cannot hide a large intermediate radial excursion or secular walk.
      double minimumMean=std::numeric_limits<double>::infinity();
      double maximumMean=-std::numeric_limits<double>::infinity();
      bool finiteMeans=true;
      for (std::size_t i=previousFirst;i<completedMeanRadii_.size();++i) {
        const double value=completedMeanRadii_[i];
        if (!std::isfinite(value)) {
          finiteMeans=false;
          break;
        }
        minimumMean=std::min(minimumMean,value);
        maximumMean=std::max(maximumMean,value);
      }
      const double meanSpan=finiteMeans ? maximumMean-minimumMean :
          std::numeric_limits<double>::infinity();

      // Window means must agree and every individual turn mean in the two windows
      // must remain inside the same bound.  The second condition is the secular-leak
      // veto: phase mixing may remove oscillatory bounce/gyro phase, but it may never
      // average away a systematic radial walk.
      const bool secularDriftOk=(cfg_.driftMaxMeanRadiusChange_m<=0.0) ||
          (meanChange<=cfg_.driftMaxMeanRadiusChange_m &&
           meanSpan<=cfg_.driftMaxMeanRadiusChange_m);

      // Profile recurrence is positive physical evidence; the secular check is only
      // a veto.  Neither a long elapsed time nor a large revolution count contributes
      // toward a trapped verdict on its own.
      DriftComparisonDiagnostics candidate;
      candidate.recurrenceLagRevolutions=static_cast<int>(lag);
      candidate.secularDriftPassed=secularDriftOk;
      candidate.windowMeanRadiusChange_m=meanChange;
      candidate.windowMeanRadiusSpan_m=meanSpan;
      const bool recurrent=ProfilesRecur(previous,current,candidate);
      candidate.profilesRecurred=recurrent;
      if (recurrent && secularDriftOk)
        ++stableComparisonsByLag_[lag];
      else
        stableComparisonsByLag_[lag]=0;

      const int stable=stableComparisonsByLag_[lag];
      stableDriftComparisons_=std::max(stableDriftComparisons_,stable);

      // Prefer a candidate with demonstrated consecutive recurrence.  Before any
      // lag passes, retain the closest failed comparison by matched-bin fraction;
      // this makes a TIME_LIMIT diagnostic identify the strongest candidate instead
      // of merely reporting the numerically last lag examined.
      if (!haveDiagnostic || stable>bestStable ||
          (stable==bestStable && candidate.matchedBinFraction>bestMatched)) {
        haveDiagnostic=true;
        bestStable=stable;
        bestMatched=candidate.matchedBinFraction;
        lastDriftComparison_=candidate;
        lastDriftMeanRadiusChange_m_=meanChange;
        lastDriftMeanRadiusSpan_m_=meanSpan;
        driftRecurrenceLagRevolutions_=static_cast<int>(lag);
      }
    }
    currentProfile_.Reset(std::max(8,cfg_.driftProfileBins));
  }

  void UpdateDriftRecurrence(const double x[3],const double p[3],const double B[3],
                             double pMag,double bMag,double r,
                             double latitudeProxy,double pitchCos2,
                             double traceTime_s) {
    const double rho=std::hypot(x[0],x[1]);
    if (!(rho>0.0) || !std::isfinite(rho)) return;

    const double phi=std::atan2(x[1],x[0]);
    if (!driftPhiInitialized_) {
      lastDriftPhi_=phi;
      driftPhiInitialized_=true;
      AddCurrentDriftProfileSample(r,latitudeProxy,pitchCos2);
      return;
    }

    // Sum the signed, locally unwrapped azimuth increment.  Fast gyromotion can move
    // phi back and forth, but those signed excursions cancel; a systematic gradient /
    // curvature drift accumulates.  A complete revolution is recognized only from the
    // net accumulated angle, not from absolute per-step motion.
    double dphi=phi-lastDriftPhi_;
    constexpr double kPi=3.141592653589793238462643383279502884;
    constexpr double kTwoPi=2.0*kPi;
    while (dphi>kPi) dphi-=kTwoPi;
    while (dphi<-kPi) dphi+=kTwoPi;
    if (std::isfinite(dphi)) driftAccumulatedAngle_+=dphi;
    lastDriftPhi_=phi;

    const int completed=static_cast<int>(
        std::floor(std::fabs(driftAccumulatedAngle_)/kTwoPi));
    while (driftRevolutions_<completed)
      CompleteDriftRevolution(x,p,B,pMag,bMag,r,latitudeProxy,pitchCos2,traceTime_s);

    // The current sample belongs to the current (possibly just-started) revolution.
    AddCurrentDriftProfileSample(r,latitudeProxy,pitchCos2);
  }

  bool DriftCriteriaSatisfied() const {
    const int window=std::max(1,cfg_.driftWindowRevolutions);
    const int requiredComparisons=(window==1)
        ? std::max(1,cfg_.minDriftRevolutions-1)
        : std::max(1,cfg_.driftMinStableWindowComparisons);
    const int selectedLag=std::max(window,driftRecurrenceLagRevolutions_);
    const int requiredRevolutions=std::max(
        cfg_.minDriftRevolutions,
        selectedLag+window+requiredComparisons-1);
    if (!cfg_.driftEnabled || driftRevolutions_<requiredRevolutions)
      return false;

    // Window comparisons advance by one complete revolution.  Requiring consecutive
    // successes rejects isolated returns in chaotic penumbral motion.  With the legacy
    // one-turn window, this reduces exactly to N-1 adjacent-turn comparisons.
    return stableDriftComparisons_>=requiredComparisons &&
           minimumOuterMargin_>=cfg_.outerMargin_m &&
           MomentumWithinTolerance();
  }

  Config cfg_;
  Earth::TrajectoryBoundary::Box box_;
  bool initialized_{false};
  bool trapped_{false};
  Mechanism mechanism_{Mechanism::None};

  // Legacy bounce-recurrence state.
  int lastParallelSign_{0};
  int mirrorPoints_{0};
  int bounceCycles_{0};
  int stableCycleComparisons_{0};
  double lastBounceEnvelopeChange_m_{std::numeric_limits<double>::infinity()};
  bool havePreviousCycle_{false};
  double currentRMin_{std::numeric_limits<double>::infinity()};
  double currentRMax_{0.0};
  double previousRMin_{0.0};
  double previousRMax_{0.0};

  // Shared numerical/physical safety diagnostics.
  double pReference_{0.0};
  double pMin_{std::numeric_limits<double>::infinity()};
  double pMax_{0.0};
  double minimumOuterMargin_{std::numeric_limits<double>::infinity()};

  // Full-orbit drift recurrence state.
  bool driftPhiInitialized_{false};
  double lastDriftPhi_{0.0};
  double driftAccumulatedAngle_{0.0};
  int driftRevolutions_{0};
  int stableDriftComparisons_{0};
  int driftRecurrenceLagRevolutions_{0};
  double lastDriftMeanRadiusChange_m_{0.0};
  double lastDriftMeanRadiusSpan_m_{0.0};
  DriftComparisonDiagnostics lastDriftComparison_{};
  std::vector<Profile> completedProfiles_;
  std::vector<double> completedMeanRadii_;
  std::vector<int> stableComparisonsByLag_;
  std::vector<DriftTurnState> driftTurnStates_;
  Profile currentProfile_;
};

} // namespace TrajectoryTrap
} // namespace Earth

#endif

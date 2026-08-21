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

  // Secular-growth guard.  A genuinely closed or quasi-closed drift shell may breathe
  // strongly with magnetic local time, especially in T05, but its *turn-averaged* radius
  // should not march systematically outward from one complete revolution to the next.
  // This guard is deliberately much tighter than the per-bin recurrence tolerance: the
  // latter admits day/night shell splitting, whereas this term rejects a slowly escaping
  // orbit that could otherwise satisfy a loose relative per-bin match for several turns.
  double driftMaxSecularGrowthAbsolute_m{0.25*6371200.0};
  double driftMaxSecularGrowthRelative{0.03};

  // A revolution must sample enough azimuth bins to be meaningful, and a configurable
  // fraction of the bins present in both turns must satisfy all recurrence tolerances.
  // These gates prevent a sparse or accidental single-point return from being called a
  // closed drift shell.
  int driftProfileBins{24};
  double driftMinProfileCoverage{0.70};
  double driftMinMatchedBinFraction{0.75};
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
    previousProfile_.Reset(n);
  }

  bool Update(const double x[3], const double p[3], const double B[3]) {
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
      UpdateDriftRecurrence(x,r,latitudeProxy,pitchCos2);
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
  double driftAngleRadians() const { return std::fabs(driftAccumulatedAngle_); }
  double lastDriftMatchedBinFraction() const { return lastDriftMatchedBinFraction_; }
  double lastDriftRmsRadiusError_m() const { return lastDriftRmsRadiusError_m_; }
  double lastDriftSecularGrowth_m() const { return lastDriftSecularGrowth_m_; }
  double minimumOuterMargin_m() const {
    return std::isfinite(minimumOuterMargin_) ? minimumOuterMargin_ : 0.0;
  }
  double momentumRelativeSpread() const {
    return (pReference_>0.0) ? (pMax_-pMin_)/pReference_ : 0.0;
  }

 private:
  struct Profile {
    std::vector<double> rSum;
    std::vector<double> latitudeSum;
    std::vector<double> pitchCos2Sum;
    std::vector<int> count;

    void Reset(int n) {
      rSum.assign((std::size_t)n,0.0);
      latitudeSum.assign((std::size_t)n,0.0);
      pitchCos2Sum.assign((std::size_t)n,0.0);
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
    currentProfile_.rSum[(std::size_t)bin]+=r;
    currentProfile_.latitudeSum[(std::size_t)bin]+=latitudeProxy;
    currentProfile_.pitchCos2Sum[(std::size_t)bin]+=pitchCos2;
    ++currentProfile_.count[(std::size_t)bin];
  }

  bool ProfileHasCoverage(const Profile& profile) const {
    const int n=profile.bins();
    if (n<=0) return false;
    return static_cast<double>(profile.populatedBins())/static_cast<double>(n)
           >= cfg_.driftMinProfileCoverage;
  }

  bool ProfilesRecur(const Profile& previous,const Profile& current) {
    lastDriftMatchedBinFraction_=0.0;
    lastDriftRmsRadiusError_m_=std::numeric_limits<double>::infinity();
    lastDriftSecularGrowth_m_=std::numeric_limits<double>::infinity();

    if (!ProfileHasCoverage(previous) || !ProfileHasCoverage(current)) return false;
    const int n=std::min(previous.bins(),current.bins());
    int common=0;
    int matched=0;
    double sumRadiusError2=0.0;
    double sumPreviousRadius=0.0;
    double sumCurrentRadius=0.0;
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

      const double dr=r1-r0;
      sumRadiusError2+=dr*dr;
      sumPreviousRadius+=r0;
      sumCurrentRadius+=r1;

      const double rTol=std::max(cfg_.driftRadialAbsoluteTolerance_m,
                                 cfg_.driftRadialRelativeTolerance*
                                 std::max(std::fabs(r0),std::fabs(r1)));
      const bool radialOk=std::fabs(dr)<=rTol;
      const bool latitudeOk=std::fabs(z1-z0)<=cfg_.driftLatitudeTolerance;
      const bool pitchOk=std::fabs(a1-a0)<=cfg_.driftPitchCos2Tolerance;
      if (radialOk && latitudeOk && pitchOk) ++matched;
    }

    // Require enough common bins independently of each turn's own coverage.  This
    // protects against two well-sampled revolutions whose missing bins are disjoint.
    const int requiredCommon=static_cast<int>(std::ceil(
        cfg_.driftMinProfileCoverage*static_cast<double>(n)));
    if (common<requiredCommon) return false;

    lastDriftMatchedBinFraction_=static_cast<double>(matched)/static_cast<double>(common);
    lastDriftRmsRadiusError_m_=std::sqrt(sumRadiusError2/static_cast<double>(common));

    // A turn-averaged radius is intentionally used only as a *secular* guard.  The
    // detailed per-azimuth recurrence above remains the actual shell-shape test.  This
    // distinction allows a strongly asymmetric T05 shell to breathe by several Re over
    // one orbit while rejecting a trajectory whose entire shell drifts outward each
    // turn.  Inward secular motion is not rejected here because it will ultimately
    // encounter the physical inner boundary; only outward motion can masquerade as a
    // stable trapped shell while actually escaping toward the outer box.
    const double meanPrevious=sumPreviousRadius/static_cast<double>(common);
    const double meanCurrent=sumCurrentRadius/static_cast<double>(common);
    lastDriftSecularGrowth_m_=meanCurrent-meanPrevious;
    const double secularTol=std::max(cfg_.driftMaxSecularGrowthAbsolute_m,
                                     cfg_.driftMaxSecularGrowthRelative*
                                     std::max(std::fabs(meanPrevious),std::fabs(meanCurrent)));
    const bool secularGrowthOk=(lastDriftSecularGrowth_m_<=secularTol);

    return lastDriftMatchedBinFraction_>=cfg_.driftMinMatchedBinFraction &&
           secularGrowthOk;
  }

  void CompleteDriftRevolution() {
    ++driftRevolutions_;
    if (havePreviousDriftProfile_) {
      if (ProfilesRecur(previousProfile_,currentProfile_))
        ++stableDriftComparisons_;
      else
        stableDriftComparisons_=0;
    }
    previousProfile_=currentProfile_;
    havePreviousDriftProfile_=true;
    currentProfile_.Reset(std::max(8,cfg_.driftProfileBins));
  }

  void UpdateDriftRecurrence(const double x[3],double r,
                             double latitudeProxy,double pitchCos2) {
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
      CompleteDriftRevolution();

    // The current sample belongs to the current (possibly just-started) revolution.
    AddCurrentDriftProfileSample(r,latitudeProxy,pitchCos2);
  }

  bool DriftCriteriaSatisfied() const {
    if (!cfg_.driftEnabled || driftRevolutions_<cfg_.minDriftRevolutions)
      return false;

    // N completed revolutions provide N-1 independent turn-to-turn comparisons.  All
    // comparisons demanded by minDriftRevolutions must be consecutive.  Any failed
    // recurrence resets stableDriftComparisons_ to zero, so a chaotic penumbral orbit
    // cannot accumulate isolated matches over a long trace.
    const int requiredStableComparisons=std::max(1,cfg_.minDriftRevolutions-1);
    return stableDriftComparisons_>=requiredStableComparisons &&
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
  bool havePreviousDriftProfile_{false};
  // Metrics from the most recent complete turn-to-turn profile comparison.  They are
  // serialized by DIRECT_ACCESS so C19 can distinguish "nearly recurring but not yet
  // proven trapped" from a completely non-recurring long-lived penumbral trajectory.
  double lastDriftMatchedBinFraction_{0.0};
  double lastDriftRmsRadiusError_m_{std::numeric_limits<double>::infinity()};
  double lastDriftSecularGrowth_m_{std::numeric_limits<double>::infinity()};
  Profile currentProfile_;
  Profile previousProfile_;
};

} // namespace TrajectoryTrap
} // namespace Earth

#endif

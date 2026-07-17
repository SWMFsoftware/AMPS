#ifndef _SRC_EARTH_UTIL_TRAJECTORY_TRAP_DETECTOR_H_
#define _SRC_EARTH_UTIL_TRAJECTORY_TRAP_DETECTOR_H_

#include "TrajectoryBoundary.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Earth {
namespace TrajectoryTrap {

// Conservative trapped-orbit detector for static magnetic fields.
//
// The detector never infers trapping from elapsed time alone.  It requires repeated
// parallel-velocity reversals (mirror points), several complete bounce cycles, a
// stable radial envelope, adequate clearance from the outer Cartesian boundary, and
// bounded momentum-magnitude variation.  The caller must enable it deliberately only
// for a field configuration where a trapped-orbit classification is physically
// meaningful (for example the static centered dipole used by F3).
struct Config {
  bool enabled{false};
  int minMirrorPoints{8};
  int minBounceCycles{4};
  double outerMargin_m{0.0};
  double radialEnvelopeTolerance_m{0.0};
  double energyRelativeTolerance{1.0e-4};
  double parallelDeadband{1.0e-6};
};

class Detector {
 public:
  Detector(const Config& config, const Earth::TrajectoryBoundary::Box& box)
      : cfg_(config), box_(box) {}

  bool Update(const double x[3], const double p[3], const double B[3]) {
    if (!cfg_.enabled || trapped_) return trapped_;
    if (!Finite3(x) || !Finite3(p) || !Finite3(B)) return false;

    const double pMag=Norm3(p);
    const double bMag=Norm3(B);
    const double r=Norm3(x);
    if (!(pMag>0.0) || !(bMag>0.0) || !std::isfinite(r)) return false;

    if (!initialized_) {
      initialized_=true;
      pReference_=pMag;
      pMin_=pMag;
      pMax_=pMag;
      ResetCurrentCycle(r,OuterMargin(x));
    }
    else {
      pMin_=std::min(pMin_,pMag);
      pMax_=std::max(pMax_,pMag);
      currentRMin_=std::min(currentRMin_,r);
      currentRMax_=std::max(currentRMax_,r);
      minimumOuterMargin_=std::min(minimumOuterMargin_,OuterMargin(x));
    }

    const double mu=(p[0]*B[0]+p[1]*B[1]+p[2]*B[2])/(pMag*bMag);
    const int sign=(mu>cfg_.parallelDeadband) ? 1 :
                   ((mu<-cfg_.parallelDeadband) ? -1 : 0);
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
    // with the default four-cycle minimum, this prevents a single fortuitously
    // similar pair of bounce envelopes from being classified as trapped.
    const int requiredStableComparisons=2;
    const double pSpread=(pReference_>0.0) ? (pMax_-pMin_)/pReference_
                                           : std::numeric_limits<double>::infinity();
    trapped_=(mirrorPoints_>=cfg_.minMirrorPoints &&
              bounceCycles_>=cfg_.minBounceCycles &&
              stableCycleComparisons_>=requiredStableComparisons &&
              minimumOuterMargin_>=cfg_.outerMargin_m &&
              pSpread<=cfg_.energyRelativeTolerance);

    ResetCurrentCycle(r,OuterMargin(x));
    return trapped_;
  }

  bool trapped() const { return trapped_; }
  int mirrorPoints() const { return mirrorPoints_; }
  int bounceCycles() const { return bounceCycles_; }
  double momentumRelativeSpread() const {
    return (pReference_>0.0) ? (pMax_-pMin_)/pReference_ : 0.0;
  }

 private:
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

  void ResetCurrentCycle(double r, double outerMargin) {
    currentRMin_=r;
    currentRMax_=r;
    minimumOuterMargin_=std::min(minimumOuterMargin_,outerMargin);
  }

  Config cfg_;
  Earth::TrajectoryBoundary::Box box_;
  bool initialized_{false};
  bool trapped_{false};
  int lastParallelSign_{0};
  int mirrorPoints_{0};
  int bounceCycles_{0};
  int stableCycleComparisons_{0};
  bool havePreviousCycle_{false};
  double pReference_{0.0};
  double pMin_{std::numeric_limits<double>::infinity()};
  double pMax_{0.0};
  double currentRMin_{std::numeric_limits<double>::infinity()};
  double currentRMax_{0.0};
  double previousRMin_{0.0};
  double previousRMax_{0.0};
  double minimumOuterMargin_{std::numeric_limits<double>::infinity()};
};

} // namespace TrajectoryTrap
} // namespace Earth

#endif

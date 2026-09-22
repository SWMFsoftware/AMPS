#ifndef _SRC_EARTH_UTIL_DIRECTIONAL_ACCESS_H_
#define _SRC_EARTH_UTIL_DIRECTIONAL_ACCESS_H_

//======================================================================================
// DirectionalAccess.h
//======================================================================================
// Production data contract and field-independent reconstruction kernels for A(E,Omega).
//
// Step 5 promotes directional access from a temporary scalar-transmissivity input to a
// saved physical product.  One sample therefore retains the requested energy/rigidity,
// direction and quadrature weight, exact three-state classification, detailed
// termination, retry provenance, and the boundary exit state returned by Step 4.
//
// The routines below reconstruct cutoff diagnostics from the saved samples alone.  No
// field evaluator or particle mover is called.  This is an intentional release
// invariant: if a cutoff in an output file cannot be reproduced from its A(E,Omega)
// rows, the producer has discarded information and the result is rejected.
//======================================================================================

#include "CutoffBandSearch.h"
#include "TrajectoryContract.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Earth {
namespace DirectionalAccess {

inline double RegularLonLatCellWeightSr(double longitudeWidth_deg,
                                        double latitudeCenter_deg,
                                        double latitudeHeight_deg) {
  if (!(longitudeWidth_deg>0.0) || !(latitudeHeight_deg>0.0) ||
      !std::isfinite(longitudeWidth_deg) ||
      !std::isfinite(latitudeCenter_deg) ||
      !std::isfinite(latitudeHeight_deg) ||
      latitudeCenter_deg<-90.0 || latitudeCenter_deg>90.0)
    throw std::invalid_argument("invalid regular lon/lat angular cell");
  const double pi=std::acos(-1.0);
  const double latitude=latitudeCenter_deg*pi/180.0;
  const double halfHeight=0.5*latitudeHeight_deg*pi/180.0;
  const double lower=std::max(-0.5*pi,latitude-halfHeight);
  const double upper=std::min(0.5*pi,latitude+halfHeight);
  return longitudeWidth_deg*pi/180.0*(std::sin(upper)-std::sin(lower));
}

struct Sample {
  double energy_MeV{0.0};
  double rigidity_GV{0.0};
  double direction_unit[3]{0.0,0.0,1.0};
  double directionWeight_sr{0.0};
  double responseWeight{1.0};
  EarthUtil::CutoffSampleState state{EarthUtil::CutoffSampleState::Unresolved};
  GridlessMode::TrajectoryTermination termination{
      GridlessMode::TrajectoryTermination::NumericalFailure};
  Trajectory::ExitState exitState{};
  int retryCount{0};
  int traceExtensionCount{0};
};

struct CutoffDiagnostics {
  double lower_GV{std::numeric_limits<double>::quiet_NaN()};
  double effective_GV{std::numeric_limits<double>::quiet_NaN()};
  double upper_GV{std::numeric_limits<double>::quiet_NaN()};
  double penumbraWidth_GV{std::numeric_limits<double>::quiet_NaN()};

  // Conservative bounds obtained by treating every unresolved interval first as
  // forbidden and then as allowed.  They are product uncertainty fields, not fit
  // errors or statistical confidence intervals.
  double effectiveLower_GV{std::numeric_limits<double>::quiet_NaN()};
  double effectiveUpper_GV{std::numeric_limits<double>::quiet_NaN()};
  double responseWeightedUnresolvedSupport{0.0};

  int transitions{0};
  int allowedIntervals{0};
  int unresolvedSamples{0};
  bool reconstructable{false};
};

inline void ValidateCurve(const std::vector<Sample>& samples) {
  if (samples.size()<2)
    throw std::invalid_argument("directional access reconstruction needs at least two samples");
  for (std::size_t i=0;i<samples.size();++i) {
    const Sample& sample=samples[i];
    if (!(sample.rigidity_GV>0.0) || !std::isfinite(sample.rigidity_GV) ||
        !(sample.energy_MeV>=0.0) || !std::isfinite(sample.energy_MeV) ||
        !(sample.directionWeight_sr>=0.0) ||
        !std::isfinite(sample.directionWeight_sr) ||
        !(sample.responseWeight>=0.0) || !std::isfinite(sample.responseWeight))
      throw std::invalid_argument("directional access sample has invalid coordinates or weights");
    if (i>0 && !(sample.rigidity_GV>samples[i-1].rigidity_GV))
      throw std::invalid_argument("directional access rigidities must be strictly increasing");
    if (sample.state==EarthUtil::CutoffSampleState::Allowed &&
        (!GridlessMode::IsAllowedTermination(sample.termination) ||
         !sample.exitState.valid))
      throw std::invalid_argument(
          "ALLOWED directional access sample must retain a valid outer-boundary state");
    if (sample.state==EarthUtil::CutoffSampleState::PhysicalForbidden &&
        !GridlessMode::IsPhysicalForbiddenTermination(sample.termination))
      throw std::invalid_argument(
          "PHYSICAL_FORBIDDEN sample has a non-physical termination reason");
  }
}

inline CutoffDiagnostics ReconstructCutoff(const std::vector<Sample>& samples) {
  ValidateCurve(samples);
  CutoffDiagnostics result;

  std::vector<EarthUtil::CutoffSampleState> states;
  states.reserve(samples.size());
  for (const Sample& sample:samples) states.push_back(sample.state);
  const EarthUtil::CutoffBandTopology topology=
      EarthUtil::AnalyzeCutoffBandSamples(states);
  result.transitions=topology.nTransitions;
  result.allowedIntervals=topology.nAllowedIntervals;
  result.unresolvedSamples=topology.nUnresolved;

  if (topology.lowerBelowRange) result.lower_GV=samples.front().rigidity_GV;
  else if (topology.lowerAllowedIndex>=0)
    result.lower_GV=samples[static_cast<std::size_t>(topology.lowerAllowedIndex)].rigidity_GV;

  if (topology.upperBelowRange) result.upper_GV=samples.front().rigidity_GV;
  else if (topology.upperAllowedIndex>=0)
    result.upper_GV=samples[static_cast<std::size_t>(topology.upperAllowedIndex)].rigidity_GV;

  if (std::isfinite(result.lower_GV) && std::isfinite(result.upper_GV) &&
      result.upper_GV>=result.lower_GV)
    result.penumbraWidth_GV=result.upper_GV-result.lower_GV;

  double allowedNominal=0.0;
  double allowedLower=0.0;
  double allowedUpper=0.0;
  double responseTotal=0.0;
  double responseUnresolved=0.0;

  for (std::size_t i=0;i+1<samples.size();++i) {
    const Sample& left=samples[i];
    const Sample& right=samples[i+1];
    const double fullWidth=right.rigidity_GV-left.rigidity_GV;
    const double response=0.5*(left.responseWeight+right.responseWeight)*fullWidth;
    responseTotal+=response;

    const bool leftAllowed=left.state==EarthUtil::CutoffSampleState::Allowed;
    const bool rightAllowed=right.state==EarthUtil::CutoffSampleState::Allowed;
    const bool unresolved=left.state==EarthUtil::CutoffSampleState::Unresolved ||
                          right.state==EarthUtil::CutoffSampleState::Unresolved;
    if (unresolved) {
      responseUnresolved+=response;
    }

    // Rc_effective integrates only across [Rc_lower,Rc_upper].  Samples below the
    // first access boundary and above the continuously allowed branch remain in the
    // saved curve (and in the unresolved-support diagnostic), but cannot contribute
    // to penumbral allowed width.
    if (!std::isfinite(result.lower_GV) || !std::isfinite(result.upper_GV)) continue;
    const double a=std::max(left.rigidity_GV,result.lower_GV);
    const double b=std::min(right.rigidity_GV,result.upper_GV);
    if (!(b>a)) continue;
    const double width=b-a;
    if (unresolved) {
      allowedUpper+=width;
      continue;
    }

    // Piecewise-linear reconstruction of the binary access indicator.  Equal endpoint
    // states integrate exactly; a resolved transition contributes half the bracket
    // width and its remaining width is the declared sampling uncertainty.
    const double contribution=(leftAllowed?0.5:0.0)+(rightAllowed?0.5:0.0);
    allowedNominal+=contribution*width;
    allowedLower+=contribution*width;
    allowedUpper+=contribution*width;
  }

  if (responseTotal>0.0)
    result.responseWeightedUnresolvedSupport=responseUnresolved/responseTotal;

  if (std::isfinite(result.upper_GV)) {
    result.effective_GV=result.upper_GV-allowedNominal;
    // More allowed support lowers effective cutoff.  Hence the naming below is
    // deliberately reversed relative to allowed-width bounds.
    result.effectiveLower_GV=result.upper_GV-allowedUpper;
    result.effectiveUpper_GV=result.upper_GV-allowedLower;
  }

  result.reconstructable=std::isfinite(result.lower_GV) &&
      std::isfinite(result.effective_GV) && std::isfinite(result.upper_GV) &&
      !topology.lowerBracketUnresolved && !topology.upperBracketUnresolved;
  return result;
}

} // namespace DirectionalAccess
} // namespace Earth

#endif

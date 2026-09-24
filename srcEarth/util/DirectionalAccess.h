#ifndef _SRC_EARTH_UTIL_DIRECTIONAL_ACCESS_H_
#define _SRC_EARTH_UTIL_DIRECTIONAL_ACCESS_H_

//======================================================================================
// DirectionalAccess.h
//======================================================================================
// Backend-independent data contract and reconstruction kernels for A(E,Omega).
//
// Roadmap Step 5 promotes directional access from an internal trajectory result to a
// saved physical product.  A sample therefore carries its energy/rigidity coordinate,
// look direction and exact solid-angle weight, three-state access classification,
// termination/retry provenance, and (for ALLOWED characteristics) the complete Step-4
// outer-boundary phase-space state.
//
// ReconstructCutoff() intentionally consumes only saved samples.  It never calls a
// field model or mover.  This is an important reproducibility invariant: lower,
// effective, and upper cutoff diagnostics written by a producer must be derivable from
// its A(E,Omega) rows without rerunning trajectories or assuming monotonic access.
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

// Exact area of a regular longitude/latitude cell on the unit sphere.  Latitude rows
// are centered at latitudeCenter_deg; the polar rows are clipped at +/-90 degrees.
// This formula, rather than cos(latitude)*dlat*dlon, makes a complete regular grid sum
// to 4*pi including the half-height polar caps.
inline double RegularLonLatCellWeightSr(double longitudeWidth_deg,
                                        double latitudeCenter_deg,
                                        double latitudeHeight_deg) {
  if (!(longitudeWidth_deg>0.0) || longitudeWidth_deg>360.0 ||
      !(latitudeHeight_deg>0.0) || latitudeHeight_deg>180.0 ||
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

  // Optional non-negative response/spectrum factor used only for unresolved-support
  // integration.  It does not change the binary/three-state access classification.
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

  // Conservative effective-cutoff bounds obtained by treating intervals touching an
  // UNRESOLVED sample first as maximally allowed and then as forbidden.  These are
  // numerical-classification bounds, not statistical confidence intervals.
  double effectiveLower_GV{std::numeric_limits<double>::quiet_NaN()};
  double effectiveUpper_GV{std::numeric_limits<double>::quiet_NaN()};
  double responseWeightedUnresolvedSupport{0.0};

  int transitions{0};
  int allowedIntervals{0};
  int unresolvedSamples{0};
  bool reconstructable{false};
};

inline bool IsFiniteCompleteExitState(const Trajectory::ExitState& state) {
  if (!state.valid || !std::isfinite(state.cosAlpha) ||
      state.cosAlpha<-1.0 || state.cosAlpha>1.0 ||
      !std::isfinite(state.traceTimeAtExit_s) || state.traceTimeAtExit_s<0.0 ||
      !std::isfinite(state.rigidityAtExit_GV) || !(state.rigidityAtExit_GV>0.0))
    return false;

  double p2=0.0;
  double v2=0.0;
  double pDotV=0.0;
  for (int d=0;d<3;++d) {
    if (!std::isfinite(state.x_exit_m[d]) ||
        !std::isfinite(state.p_exit_SI[d]) ||
        !std::isfinite(state.v_exit_unit[d])) return false;
    p2+=state.p_exit_SI[d]*state.p_exit_SI[d];
    v2+=state.v_exit_unit[d]*state.v_exit_unit[d];
    pDotV+=state.p_exit_SI[d]*state.v_exit_unit[d];
  }
  if (!(p2>0.0) || !std::isfinite(p2) ||
      std::fabs(v2-1.0)>1.0e-10) return false;

  // Momentum and velocity direction are produced at the same Step-4 crossing event.
  // A negative or non-parallel pair would indicate that fields from different steps
  // were mixed while serializing the product.
  const double p=std::sqrt(p2);
  return pDotV>0.0 && std::fabs(pDotV/p-1.0)<=1.0e-10;
}

inline void ValidateCurve(const std::vector<Sample>& samples) {
  if (samples.size()<2)
    throw std::invalid_argument(
        "directional access reconstruction needs at least two samples");

  for (std::size_t i=0;i<samples.size();++i) {
    const Sample& sample=samples[i];
    double directionNorm2=0.0;
    for (int d=0;d<3;++d) {
      if (!std::isfinite(sample.direction_unit[d]))
        throw std::invalid_argument(
            "directional access sample has a non-finite direction");
      directionNorm2+=sample.direction_unit[d]*sample.direction_unit[d];
    }
    if (!(sample.rigidity_GV>0.0) || !std::isfinite(sample.rigidity_GV) ||
        !(sample.energy_MeV>=0.0) || !std::isfinite(sample.energy_MeV) ||
        !(sample.directionWeight_sr>0.0) ||
        !std::isfinite(sample.directionWeight_sr) ||
        !(sample.responseWeight>=0.0) || !std::isfinite(sample.responseWeight) ||
        std::fabs(directionNorm2-1.0)>1.0e-10)
      throw std::invalid_argument(
          "directional access sample has invalid coordinates or weights");
    if (i>0 && !(sample.rigidity_GV>samples[i-1].rigidity_GV))
      throw std::invalid_argument(
          "directional access rigidities must be strictly increasing");

    switch (sample.state) {
      case EarthUtil::CutoffSampleState::Allowed:
        if (!GridlessMode::IsAllowedTermination(sample.termination) ||
            !IsFiniteCompleteExitState(sample.exitState))
          throw std::invalid_argument(
              "ALLOWED directional access sample lacks a complete outer-boundary state");
        break;
      case EarthUtil::CutoffSampleState::PhysicalForbidden:
        if (!GridlessMode::IsPhysicalForbiddenTermination(sample.termination))
          throw std::invalid_argument(
              "PHYSICAL_FORBIDDEN sample has a non-physical termination reason");
        if (sample.exitState.valid)
          throw std::invalid_argument(
              "PHYSICAL_FORBIDDEN sample must not retain an allowed exit state");
        break;
      case EarthUtil::CutoffSampleState::Unresolved:
        if (GridlessMode::IsResolvedTermination(sample.termination) ||
            sample.exitState.valid)
          throw std::invalid_argument(
              "UNRESOLVED sample has resolved termination or a valid exit state");
        break;
      default:
        throw std::invalid_argument("directional access sample has an invalid state");
    }
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
    result.lower_GV=
        samples[static_cast<std::size_t>(topology.lowerAllowedIndex)].rigidity_GV;

  if (topology.upperBelowRange) result.upper_GV=samples.front().rigidity_GV;
  else if (topology.upperAllowedIndex>=0)
    result.upper_GV=
        samples[static_cast<std::size_t>(topology.upperAllowedIndex)].rigidity_GV;

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
    const double response=
        0.5*(left.responseWeight+right.responseWeight)*fullWidth;
    responseTotal+=response;

    const bool leftAllowed=
        left.state==EarthUtil::CutoffSampleState::Allowed;
    const bool rightAllowed=
        right.state==EarthUtil::CutoffSampleState::Allowed;
    const bool unresolved=
        left.state==EarthUtil::CutoffSampleState::Unresolved ||
        right.state==EarthUtil::CutoffSampleState::Unresolved;
    if (unresolved) responseUnresolved+=response;

    // Effective cutoff is defined only across the reconstructed penumbra.  Support
    // outside [lower,upper] is retained in the file and unresolved metric, but cannot
    // be counted as penumbral allowed width.
    if (!std::isfinite(result.lower_GV) ||
        !std::isfinite(result.upper_GV)) continue;
    const double a=std::max(left.rigidity_GV,result.lower_GV);
    const double b=std::min(right.rigidity_GV,result.upper_GV);
    if (!(b>a)) continue;
    const double width=b-a;

    if (unresolved) {
      // Nominal/lower assume no resolved allowed support in this interval; the upper
      // access bound assumes the full interval is allowed.
      allowedUpper+=width;
      continue;
    }

    // Piecewise-linear integration of the binary endpoint indicators.  Equal states
    // are exact on the sampled representation; a resolved transition contributes half
    // of its finite bracket.  AdaptiveDirectAccess separately reports that bracket's
    // full width as discretization support.
    const double contribution=
        (leftAllowed ? 0.5 : 0.0)+(rightAllowed ? 0.5 : 0.0);
    allowedNominal+=contribution*width;
    allowedLower+=contribution*width;
    allowedUpper+=contribution*width;
  }

  if (responseTotal>0.0)
    result.responseWeightedUnresolvedSupport=
        responseUnresolved/responseTotal;

  if (std::isfinite(result.upper_GV)) {
    result.effective_GV=result.upper_GV-allowedNominal;
    // More allowed support lowers the effective cutoff; therefore cutoff-bound names
    // are reversed relative to the corresponding allowed-width bounds.
    result.effectiveLower_GV=result.upper_GV-allowedUpper;
    result.effectiveUpper_GV=result.upper_GV-allowedLower;
  }

  result.reconstructable=
      std::isfinite(result.lower_GV) &&
      std::isfinite(result.effective_GV) &&
      std::isfinite(result.upper_GV) &&
      !topology.lowerBracketUnresolved &&
      !topology.upperBracketUnresolved;
  return result;
}

} // namespace DirectionalAccess
} // namespace Earth

#endif // _SRC_EARTH_UTIL_DIRECTIONAL_ACCESS_H_

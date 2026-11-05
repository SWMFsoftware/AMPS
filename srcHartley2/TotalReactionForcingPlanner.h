#ifndef _HARTLEY2_TOTAL_REACTION_FORCING_PLANNER_H_
#define _HARTLEY2_TOTAL_REACTION_FORCING_PLANNER_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

// Pure Stage-8 deterministic total-reaction-forcing planner.  It deliberately
// has no AMPS, MPI, mesh, particle-buffer, or RNG dependency.  Production can
// therefore preflight the complete two-branch transaction before mutating a
// parent, while focused tests prove the exponential and conservation formulas.
namespace Comet {
namespace TotalReactionForcingPlanner {

struct cSplitPlan {
  double TotalRate=0.0;
  double TimeStep=0.0;
  double OpticalDepth=0.0;

  // For a parent carrier C, the two physical carriers are C*(1-P) and C*P.
  // SplitCarrier below computes the second as an exact arithmetic remainder
  // so the stored binary64 carriers add to the original C.
  double SurvivorFraction=1.0;
  double ReactedFraction=0.0;
};

struct cCarrierSplit {
  double Survivor=0.0;
  double Reacted=0.0;
};

inline cSplitPlan BuildSplitPlan(double totalRate,double timeStep) {
  if (!std::isfinite(totalRate) || totalRate<0.0)
    throw std::invalid_argument(
        "total-reaction forcing requires a finite nonnegative rate");
  if (!std::isfinite(timeStep) || !(timeStep>0.0))
    throw std::invalid_argument(
        "total-reaction forcing requires a finite positive timestep");

  const long double opticalDepthLong=
      static_cast<long double>(totalRate)*timeStep;
  if (!std::isfinite(opticalDepthLong) ||
      opticalDepthLong>
          static_cast<long double>(std::numeric_limits<double>::max()))
    throw std::overflow_error(
        "total-reaction optical depth exceeds binary64 range");
  const double opticalDepth=static_cast<double>(opticalDepthLong);

  cSplitPlan result;
  result.TotalRate=totalRate;
  result.TimeStep=timeStep;
  result.OpticalDepth=opticalDepth;
  if (opticalDepth==0.0) return result;

  // expm1 retains small reaction probabilities that 1-exp(-tau) would lose.
  // The survivor is the final arithmetic remainder.  This deliberately makes
  // SurvivorFraction+ReactedFraction exactly one in production precision.
  result.ReactedFraction=-std::expm1(-opticalDepth);
  if (!std::isfinite(result.ReactedFraction) ||
      result.ReactedFraction<0.0 || result.ReactedFraction>1.0)
    throw std::runtime_error(
        "total-reaction forcing produced an invalid reaction fraction");
  result.SurvivorFraction=1.0-result.ReactedFraction;
  if (result.SurvivorFraction<0.0 ||
      result.SurvivorFraction+result.ReactedFraction!=1.0)
    throw std::runtime_error(
        "total-reaction forcing fractions do not close exactly");
  return result;
}

inline cCarrierSplit SplitCarrier(const cSplitPlan &plan,
    double parentCarrier) {
  if (!std::isfinite(parentCarrier) || !(parentCarrier>0.0))
    throw std::invalid_argument(
        "total-reaction forcing requires a positive finite parent carrier");
  if (!std::isfinite(plan.ReactedFraction) ||
      plan.ReactedFraction<0.0 || plan.ReactedFraction>1.0 ||
      plan.SurvivorFraction<0.0 ||
      plan.SurvivorFraction+plan.ReactedFraction!=1.0)
    throw std::invalid_argument(
        "total-reaction forcing split plan is invalid");

  cCarrierSplit result;
  const double requestedReacted=parentCarrier*plan.ReactedFraction;
  double lowerReacted=requestedReacted,upperReacted=requestedReacted;
  bool foundExact=false;
  for (int adjustment=0;adjustment<9 && !foundExact;adjustment++) {
    if (adjustment!=0) {
      if (adjustment%2==1) {
        lowerReacted=std::nextafter(lowerReacted,0.0);
        result.Reacted=lowerReacted;
      }
      else {
        upperReacted=std::nextafter(upperReacted,
            std::numeric_limits<double>::infinity());
        result.Reacted=upperReacted;
      }
    }
    else result.Reacted=requestedReacted;
    result.Survivor=parentCarrier-result.Reacted;
    foundExact=result.Survivor+result.Reacted==parentCarrier;
  }
  // For some exponents, the nearest product and its subtraction remainder
  // straddle rather than hit the parent when re-added.  Searching adjacent
  // reacted values retains the physical split to a few ulps while giving the
  // particle commit an exact binary64 conservation identity.  Adjusting the
  // reacted value as well as recomputing the remainder is essential; changing
  // the remainder alone cannot bridge every rounding interval.
  if (!std::isfinite(result.Reacted) || !std::isfinite(result.Survivor))
    throw std::runtime_error(
        "total-reaction forcing carrier split is nonfinite");
  if (result.Reacted<0.0 || result.Survivor<0.0)
    throw std::runtime_error(
        "total-reaction forcing carrier split is negative");
  if (!foundExact)
    throw std::runtime_error(
        "total-reaction forcing carrier split does not close exactly");
  return result;
}

inline void ValidateUnitUniform(double unitUniform) {
  if (!std::isfinite(unitUniform) || unitUniform<0.0 || unitUniform>=1.0)
    throw std::invalid_argument(
        "conditional-reaction uniform variate must be finite in [0,1)");
}

// Sample the physical exponential reaction time conditional on a reaction
// occurring inside this timestep.  For P=1-exp(-lambda*dt), inversion of
// F(t|t<dt)=(1-exp(-lambda*t))/P gives the expression below.  Stage 8 changes
// occurrence representation only; it must not bias this conditional position.
inline double SampleConditionalReactionTime(const cSplitPlan &plan,
    double unitUniform) {
  ValidateUnitUniform(unitUniform);
  if (!(plan.TotalRate>0.0) || !(plan.ReactedFraction>0.0) ||
      !(plan.TimeStep>0.0) || !std::isfinite(plan.TotalRate) ||
      !std::isfinite(plan.TimeStep))
    throw std::invalid_argument(
        "conditional reaction time requires an active split plan");
  const double opticalPosition=
      -std::log1p(-unitUniform*plan.ReactedFraction);
  double reactionTime=opticalPosition/plan.TotalRate;
  const double roundoff=64.0*std::numeric_limits<double>::epsilon()*
      std::max(1.0,plan.TimeStep);
  if (reactionTime<0.0 || reactionTime>plan.TimeStep+roundoff ||
      !std::isfinite(reactionTime))
    throw std::runtime_error(
        "conditional reaction time is outside its physical segment");
  if (reactionTime>plan.TimeStep) reactionTime=plan.TimeStep;
  return reactionTime;
}

struct cHazardSegment {
  double Rate=0.0;
  double Duration=0.0;
};

struct cPiecewisePlan {
  std::vector<cHazardSegment> Segments;
  std::vector<double> CumulativeStartTime;
  std::vector<double> CumulativeStartOpticalDepth;
  double TotalTime=0.0;
  double TotalOpticalDepth=0.0;
  double SurvivorFraction=1.0;
  double ReactedFraction=0.0;
};

inline cPiecewisePlan BuildPiecewisePlan(
    const std::vector<cHazardSegment> &segments) {
  if (segments.empty())
    throw std::invalid_argument("piecewise hazard plan is empty");
  cPiecewisePlan result;
  result.Segments=segments;
  result.CumulativeStartTime.reserve(segments.size());
  result.CumulativeStartOpticalDepth.reserve(segments.size());
  long double elapsed=0.0L,opticalDepth=0.0L;
  for (const cHazardSegment &segment:segments) {
    if (!std::isfinite(segment.Rate) || segment.Rate<0.0 ||
        !std::isfinite(segment.Duration) || !(segment.Duration>0.0))
      throw std::invalid_argument(
          "piecewise hazard segments require nonnegative rates and positive durations");
    if (elapsed>std::numeric_limits<double>::max() ||
        opticalDepth>std::numeric_limits<double>::max())
      throw std::overflow_error("piecewise hazard cumulative value overflow");
    result.CumulativeStartTime.push_back(static_cast<double>(elapsed));
    result.CumulativeStartOpticalDepth.push_back(
        static_cast<double>(opticalDepth));
    elapsed+=segment.Duration;
    opticalDepth+=static_cast<long double>(segment.Rate)*segment.Duration;
  }
  if (!std::isfinite(elapsed) || !std::isfinite(opticalDepth) ||
      elapsed>std::numeric_limits<double>::max() ||
      opticalDepth>std::numeric_limits<double>::max())
    throw std::overflow_error("piecewise hazard total value overflow");
  result.TotalTime=static_cast<double>(elapsed);
  result.TotalOpticalDepth=static_cast<double>(opticalDepth);
  if (result.TotalOpticalDepth>0.0)
    result.ReactedFraction=-std::expm1(-result.TotalOpticalDepth);
  result.SurvivorFraction=1.0-result.ReactedFraction;
  if (result.SurvivorFraction<0.0 || result.ReactedFraction<0.0 ||
      result.SurvivorFraction+result.ReactedFraction!=1.0)
    throw std::runtime_error("piecewise hazard fractions do not close exactly");
  return result;
}

struct cPiecewiseReactionPosition {
  std::size_t Segment=0;
  double TimeInSegment=0.0;
  double ElapsedTime=0.0;
  double OpticalDepth=0.0;
};

inline cPiecewiseReactionPosition SampleConditionalReactionPosition(
    const cPiecewisePlan &plan,double unitUniform) {
  ValidateUnitUniform(unitUniform);
  if (!(plan.TotalOpticalDepth>0.0) || !(plan.ReactedFraction>0.0) ||
      plan.Segments.empty() ||
      plan.CumulativeStartTime.size()!=plan.Segments.size() ||
      plan.CumulativeStartOpticalDepth.size()!=plan.Segments.size())
    throw std::invalid_argument(
        "conditional piecewise position requires an active valid plan");
  const double target=-std::log1p(-unitUniform*plan.ReactedFraction);
  std::size_t lastPositive=plan.Segments.size();
  for (std::size_t i=0;i<plan.Segments.size();i++) {
    const cHazardSegment &segment=plan.Segments[i];
    if (!(segment.Rate>0.0)) continue;
    lastPositive=i;
    const double start=plan.CumulativeStartOpticalDepth[i];
    const double end=start+segment.Rate*segment.Duration;
    if (target<=end) {
      cPiecewiseReactionPosition result;
      result.Segment=i;
      result.OpticalDepth=target;
      result.TimeInSegment=(target-start)/segment.Rate;
      if (result.TimeInSegment<0.0) result.TimeInSegment=0.0;
      if (result.TimeInSegment>segment.Duration)
        result.TimeInSegment=segment.Duration;
      result.ElapsedTime=plan.CumulativeStartTime[i]+result.TimeInSegment;
      return result;
    }
  }
  if (lastPositive==plan.Segments.size())
    throw std::runtime_error("piecewise plan lost every positive hazard");

  // Only a final binary64 boundary residue can reach this fallback.  Return
  // the end of the final positive segment, never a later zero-hazard region.
  cPiecewiseReactionPosition result;
  result.Segment=lastPositive;
  result.TimeInSegment=plan.Segments[lastPositive].Duration;
  result.ElapsedTime=plan.CumulativeStartTime[lastPositive]+
      result.TimeInSegment;
  result.OpticalDepth=plan.TotalOpticalDepth;
  return result;
}

struct cCapacityPlan {
  std::size_t ConcreteProductsPerPacket=0;
  std::size_t PacketCount=0;
  std::size_t ReactionDaughters=0;
  std::size_t OutputParticles=0;
  std::size_t RequiredNewSlots=0;
};

// A normal forced transaction has one survivor plus K*N reaction daughters.
// One output reuses the parent slot, so it needs exactly K*N new slots.  At
// numerical P=0 or P=1 the zero-carrier branch is suppressed and the exact
// bound is correspondingly smaller.  This check must run before recoil RNG,
// diagnostics, correction state, or the parent list is mutated.
inline cCapacityPlan BuildCapacityPlan(const cSplitPlan &split,
    std::size_t concreteProductsPerPacket,std::size_t packetCount,
    std::size_t maximumPacketCount,std::size_t maximumReactionDaughters) {
  if (packetCount==0 || maximumPacketCount==0 ||
      packetCount>maximumPacketCount)
    throw std::invalid_argument("forced-reaction packet count is invalid");
  if (concreteProductsPerPacket!=0 &&
      packetCount>std::numeric_limits<std::size_t>::max()/
          concreteProductsPerPacket)
    throw std::overflow_error("forced-reaction daughter count overflow");
  const std::size_t reactionDaughters=
      concreteProductsPerPacket*packetCount;
  if (reactionDaughters>maximumReactionDaughters)
    throw std::length_error(
        "forced-reaction daughters exceed the transaction limit");

  const bool survivorActive=split.SurvivorFraction>0.0;
  const bool reactionActive=split.ReactedFraction>0.0;
  const std::size_t activeDaughters=reactionActive ? reactionDaughters : 0;
  if (activeDaughters==std::numeric_limits<std::size_t>::max() &&
      survivorActive)
    throw std::overflow_error("forced-reaction output count overflow");

  cCapacityPlan result;
  result.ConcreteProductsPerPacket=concreteProductsPerPacket;
  result.PacketCount=packetCount;
  result.ReactionDaughters=activeDaughters;
  result.OutputParticles=activeDaughters+(survivorActive ? 1 : 0);
  result.RequiredNewSlots=result.OutputParticles>0 ?
      result.OutputParticles-1 : 0;
  return result;
}

inline bool FitsAvailableSlots(const cCapacityPlan &plan,
    std::size_t availableNewSlots) {
  return plan.RequiredNewSlots<=availableNewSlots;
}

} // namespace TotalReactionForcingPlanner
} // namespace Comet

#endif

#ifndef _HARTLEY2_CORRECTION_WINDOW_PLANNER_H_
#define _HARTLEY2_CORRECTION_WINDOW_PLANNER_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

// This header contains only the pure planning half of the individual
// correction-window controller.  It deliberately has no AMPS, MPI, mesh, or
// particle-buffer dependency.  The runtime can therefore construct and
// capacity-check a complete batch before it changes a particle list.  Keeping
// random roulette decisions and particle cloning outside the planner also
// makes a rejected plan exactly retryable.
namespace Comet {
namespace CorrectionWindowPlanner {

enum class eAction {
  Keep,
  Roulette,
  Split
};

struct cPolicy {
  bool Enabled;
  double LowerCorrection,TargetCorrection,UpperCorrection;
  std::size_t MaximumSplitParticleNumber;

  cPolicy():Enabled(false),LowerCorrection(1.0),TargetCorrection(1.0),
      UpperCorrection(1.0),MaximumSplitParticleNumber(1) {}

  cPolicy(double lower,double target,double upper,
      std::size_t maximumSplitParticleNumber):Enabled(true),
      LowerCorrection(lower),TargetCorrection(target),
      UpperCorrection(upper),
      MaximumSplitParticleNumber(maximumSplitParticleNumber) {}
};

struct cRequest {
  int Species;
  double Correction;

  cRequest(int species,double correction):Species(species),
      Correction(correction) {}
};

struct cParticlePlan {
  int Species;
  eAction Action;
  double InputCorrection;

  // Keep and Split have probability one.  Roulette stores the conditional
  // survivor correction in OutputCorrections[0] and this survival
  // probability.  The expectation is therefore exactly the input correction.
  double SurvivalProbability;
  std::vector<double> OutputCorrections;

  cParticlePlan():Species(-1),Action(eAction::Keep),InputCorrection(0.0),
      SurvivalProbability(1.0) {}
};

struct cBatchPlan {
  std::vector<cParticlePlan> ParticlePlans;
  std::size_t RequiredAdditionalParticles;

  cBatchPlan():RequiredAdditionalParticles(0) {}
};

inline void ValidatePolicy(const cPolicy &policy) {
  if (!policy.Enabled) return;

  // Strict ordering gives both sides of the window a real hysteresis band.
  // Degenerate lower==target or target==upper policies cause particles on a
  // numerical boundary to alternate between roulette and splitting.
  if (!(policy.LowerCorrection>0.0) ||
      !(policy.LowerCorrection<policy.TargetCorrection) ||
      !(policy.TargetCorrection<policy.UpperCorrection) ||
      !std::isfinite(policy.LowerCorrection) ||
      !std::isfinite(policy.TargetCorrection) ||
      !std::isfinite(policy.UpperCorrection))
    throw std::invalid_argument(
        "correction window requires finite 0<lower<target<upper");

  if (policy.MaximumSplitParticleNumber<2)
    throw std::invalid_argument(
        "enabled correction window requires a split limit of at least two");

  // The required Stage-5 split rule is n=ceil(c/c_target).  Immediately above
  // the inclusive upper boundary it produces the smallest possible child
  // correction.  Reject a window whose lower edge would make that mandatory
  // result illegal; otherwise a finite, valid input could fail only because
  // of a contradictory policy rather than the explicit split/capacity bound.
  const double firstSplitCount=
      std::floor(policy.UpperCorrection/policy.TargetCorrection)+1.0;
  if (!std::isfinite(firstSplitCount) || firstSplitCount<2.0 ||
      policy.LowerCorrection>
      policy.UpperCorrection/firstSplitCount)
    throw std::invalid_argument(
        "correction window is incompatible with ceil(c/target) splitting");
}

inline std::size_t CheckedCount(double value,const char *message) {
  if (!(value>=0.0) || !std::isfinite(value) ||
      static_cast<long double>(value)>
          static_cast<long double>(std::numeric_limits<std::size_t>::max()))
    throw std::length_error(message);
  return static_cast<std::size_t>(value);
}

inline cParticlePlan BuildParticlePlan(const cRequest &request,
    const cPolicy &policy) {
  ValidatePolicy(policy);
  if (request.Species<0)
    throw std::invalid_argument("negative species in correction-window request");
  if (!(request.Correction>0.0) || !std::isfinite(request.Correction))
    throw std::invalid_argument(
        "particle correction must be positive and finite");

  cParticlePlan result;
  result.Species=request.Species;
  result.InputCorrection=request.Correction;

  if (!policy.Enabled ||
      (request.Correction>=policy.LowerCorrection &&
       request.Correction<=policy.UpperCorrection)) {
    result.Action=eAction::Keep;
    result.OutputCorrections.push_back(request.Correction);
    return result;
  }

  if (request.Correction<policy.LowerCorrection) {
    result.Action=eAction::Roulette;
    result.SurvivalProbability=request.Correction/policy.TargetCorrection;
    if (!(result.SurvivalProbability>0.0) ||
        !(result.SurvivalProbability<1.0) ||
        !std::isfinite(result.SurvivalProbability))
      throw std::invalid_argument(
          "invalid correction-window roulette probability");
    result.OutputCorrections.push_back(policy.TargetCorrection);
    return result;
  }

  result.Action=eAction::Split;

  // Stage 5 deliberately fixes the split rule rather than letting capacity or
  // implementation details choose a statistical representation.  Ceil puts
  // every equal child at or below the target; the policy compatibility check
  // above guarantees the lower window edge as well.
  std::size_t splitCount=CheckedCount(
      std::ceil(request.Correction/policy.TargetCorrection),
      "correction-window split count is not representable");
  splitCount=std::max<std::size_t>(2,splitCount);
  if (splitCount>policy.MaximumSplitParticleNumber)
    throw std::length_error(
        "split limit cannot place correction inside the configured window");

  const double equalCorrection=request.Correction/splitCount;
  if (!(equalCorrection>=policy.LowerCorrection) ||
      !(equalCorrection<=policy.UpperCorrection) ||
      !std::isfinite(equalCorrection))
    throw std::runtime_error(
        "internal correction-window split did not land inside its window");

  result.OutputCorrections.reserve(splitCount);
  double assigned=0.0;
  for (std::size_t child=1;child<splitCount;child++) {
    result.OutputCorrections.push_back(equalCorrection);
    assigned+=equalCorrection;
  }

  // Put the floating-point remainder in the final child instead of dropping
  // it.  The production executor assigns children in this same order.  The
  // explicit equality check below makes carrier reconstruction a hard
  // contract rather than accepting a merely small relative error.
  const double finalCorrection=request.Correction-assigned;
  if (!(finalCorrection>0.0) || !std::isfinite(finalCorrection) ||
      finalCorrection<policy.LowerCorrection ||
      finalCorrection>policy.UpperCorrection)
    throw std::runtime_error("invalid final correction in split plan");
  if (assigned+finalCorrection!=request.Correction)
    throw std::runtime_error(
        "split corrections do not exactly reconstruct the input correction");
  result.OutputCorrections.push_back(finalCorrection);
  return result;
}

inline cBatchPlan BuildPlan(const std::vector<cRequest> &requests,
    const std::vector<cPolicy> &speciesPolicies,
    std::size_t availableAdditionalParticles=
        std::numeric_limits<std::size_t>::max()) {
  cBatchPlan result;
  result.ParticlePlans.reserve(requests.size());

  // Build the whole plan without touching a particle.  Only after every
  // request and per-particle limit is valid do we compare the aggregate clone
  // requirement with capacity.  A failure therefore cannot partially mutate
  // an earlier species or particle.
  for (const cRequest &request:requests) {
    if (request.Species<0 ||
        static_cast<std::size_t>(request.Species)>=speciesPolicies.size())
      throw std::out_of_range(
          "correction-window request has no species policy");
    cParticlePlan particlePlan=BuildParticlePlan(
        request,speciesPolicies[request.Species]);
    const std::size_t additional=
        particlePlan.OutputCorrections.size()-1;
    if (additional>std::numeric_limits<std::size_t>::max()-
        result.RequiredAdditionalParticles)
      throw std::length_error(
          "correction-window batch clone count overflow");
    result.RequiredAdditionalParticles+=additional;
    result.ParticlePlans.push_back(particlePlan);
  }

  if (result.RequiredAdditionalParticles>availableAdditionalParticles)
    throw std::length_error(
        "insufficient particle capacity for correction-window batch");
  return result;
}

} // namespace CorrectionWindowPlanner
} // namespace Comet

#endif

#ifndef _HARTLEY2_CHANNEL_IMPORTANCE_PLANNER_H_
#define _HARTLEY2_CHANNEL_IMPORTANCE_PLANNER_H_

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

// Pure Stage-7 reaction-channel importance planner.  This header deliberately
// has no AMPS, MPI, mesh, particle, or random-number dependency.  The input
// layer can validate a requested q* family, the production executor can select
// a channel from one supplied uniform variate, and focused tests can prove the
// likelihood/variance identities without constructing a simulation.
namespace Comet {
namespace ChannelImportancePlanner {

struct cPlan {
  double TotalPhysicalRate=0.0;
  bool ImportanceSamplingEnabled=false;

  // All arrays retain the caller's stable channel order.  A disabled physical
  // channel has zero in every array and is never selected.  Positive channels
  // use normalized probabilities whose production-order sum is exactly one.
  std::vector<double> PhysicalProbability;
  std::vector<double> BiasedProbability;
  std::vector<double> LikelihoodRatio;

  // Per-event variance of L_j I(selected=j), whose expectation is q_j.
  // This is q_j^2(1/q*_j-1), and is zero for disabled channels.
  std::vector<double> WeightedChannelVariance;

  // E_q*[L^2].  For M independent selected channels, the analytic family
  // effective sample size is M/ExpectedLikelihoodSquare.
  double ExpectedLikelihoodSquare=0.0;
  double AsymptoticEffectiveSampleFraction=0.0;
};

inline void ValidateUnitUniform(double unitUniform) {
  if (!std::isfinite(unitUniform) || unitUniform<0.0 || unitUniform>=1.0)
    throw std::invalid_argument(
        "channel-selection uniform variate must be finite in [0,1)");
}

// Preserve the pre-Stage-7 analog selector as an explicit helper.  Production
// compatibility depends on multiplying the one existing uniform draw by the
// raw total rate and subtracting raw rates left-to-right.  Do not replace this
// path with normalized cumulative probabilities when importance mode is off:
// even mathematically equivalent boundaries can change binary64 tie behavior.
inline std::size_t SelectAnalogRateSubtract(
    const std::vector<double> &physicalRates,double totalPhysicalRate,
    double unitUniform) {
  ValidateUnitUniform(unitUniform);
  if (physicalRates.empty())
    throw std::invalid_argument("analog channel family is empty");
  if (!(totalPhysicalRate>0.0) || !std::isfinite(totalPhysicalRate))
    throw std::invalid_argument("analog channel total rate is invalid");

  long double verifiedTotal=0.0L;
  std::size_t lastPositive=physicalRates.size();
  for (std::size_t channel=0;channel<physicalRates.size();channel++) {
    const double rate=physicalRates[channel];
    if (!std::isfinite(rate) || rate<0.0)
      throw std::invalid_argument("analog channel rate is invalid");
    verifiedTotal+=static_cast<long double>(rate);
    if (rate>0.0) lastPositive=channel;
  }
  if (lastPositive==physicalRates.size())
    throw std::invalid_argument("analog channel family has no positive rate");
  const long double scale=std::max(
      std::fabs(verifiedTotal),
      std::fabs(static_cast<long double>(totalPhysicalRate)));
  const long double tolerance=64.0L*
      std::numeric_limits<double>::epsilon()*
      static_cast<long double>(physicalRates.size())*scale;
  if (std::fabs(verifiedTotal-
      static_cast<long double>(totalPhysicalRate))>tolerance)
    throw std::invalid_argument(
        "analog channel total does not equal the supplied rates");

  double selector=unitUniform*totalPhysicalRate;
  for (std::size_t channel=0;channel<physicalRates.size();channel++) {
    const double rate=physicalRates[channel];
    if (!(rate>0.0)) continue;
    selector-=rate;
    if (selector<=0.0) return channel;
  }

  // Rounding in unitUniform*totalPhysicalRate or in the subtraction sequence
  // can leave a tiny positive residue.  The historical production code used
  // the last active channel as the final interval; retain that exact contract.
  return lastPositive;
}

namespace Detail {

inline std::vector<double> CanonicalNormalizedProbabilities(
    const std::vector<double> &requested,
    const std::vector<double> &physicalRates,const char *name) {
  if (requested.size()!=physicalRates.size())
    throw std::invalid_argument(
        "channel probability and physical-rate sizes differ");

  long double requestedSum=0.0L;
  std::size_t lastPositive=physicalRates.size();
  for (std::size_t channel=0;channel<physicalRates.size();channel++) {
    const double probability=requested[channel];
    if (!std::isfinite(probability) || probability<0.0 || probability>1.0)
      throw std::invalid_argument(
          "channel probability is not finite in [0,1]");
    if (physicalRates[channel]>0.0) {
      if (!(probability>0.0))
        throw std::invalid_argument(
            "a positive physical channel has zero sampling probability");
      lastPositive=channel;
    }
    else if (probability!=0.0) {
      throw std::invalid_argument(
          "a disabled physical channel has nonzero sampling probability");
    }
    requestedSum+=static_cast<long double>(probability);
  }
  if (lastPositive==physicalRates.size())
    throw std::invalid_argument("channel family has no positive rate");

  const long double tolerance=64.0L*
      std::numeric_limits<double>::epsilon()*
      static_cast<long double>(requested.size());
  if (std::fabs(requestedSum-1.0L)>tolerance)
    throw std::invalid_argument(name);

  // Copy every requested interval except the final active one.  Assign that
  // interval the arithmetic remainder in production order.  This provides an
  // exact terminal boundary at one without silently renormalizing a malformed
  // policy or allowing a floating-point selection hole.
  std::vector<double> result=requested;
  double assigned=0.0;
  for (std::size_t channel=0;channel<result.size();channel++) {
    if (physicalRates[channel]<=0.0 || channel==lastPositive) continue;
    assigned+=result[channel];
  }
  result[lastPositive]=1.0-assigned;
  if (!(result[lastPositive]>0.0) ||
      !std::isfinite(result[lastPositive]) ||
      assigned+result[lastPositive]!=1.0)
    throw std::runtime_error(
        "final channel probability remainder is not positive and exact");
  if (std::fabs(static_cast<long double>(result[lastPositive])-
      static_cast<long double>(requested[lastPositive]))>tolerance)
    throw std::invalid_argument(name);
  return result;
}

} // namespace Detail

// An empty requestedBiasedProbability vector means strict analog mode.  In
// that mode the caller must still select with SelectAnalogRateSubtract(); the
// normalized arrays are diagnostic/analytic values and do not replace the raw
// rate path.  A nonempty vector is the exact q* policy for the whole family.
inline cPlan BuildPlan(const std::vector<double> &physicalRates,
    const std::vector<double> &requestedBiasedProbability={}) {
  if (physicalRates.empty())
    throw std::invalid_argument("channel family is empty");

  long double total=0.0L;
  for (double rate:physicalRates) {
    if (!std::isfinite(rate) || rate<0.0)
      throw std::invalid_argument(
          "physical channel rate must be finite and nonnegative");
    total+=static_cast<long double>(rate);
  }
  if (!(total>0.0L) ||
      total>static_cast<long double>(std::numeric_limits<double>::max()))
    throw std::invalid_argument("channel family total rate is invalid");

  cPlan result;
  result.TotalPhysicalRate=static_cast<double>(total);
  std::vector<double> rawPhysical(physicalRates.size(),0.0);
  for (std::size_t channel=0;channel<physicalRates.size();channel++)
    if (physicalRates[channel]>0.0)
      rawPhysical[channel]=static_cast<double>(
          static_cast<long double>(physicalRates[channel])/total);
  result.PhysicalProbability=Detail::CanonicalNormalizedProbabilities(
      rawPhysical,physicalRates,
      "physical channel probabilities do not normalize to one");

  if (requestedBiasedProbability.empty()) {
    result.BiasedProbability=result.PhysicalProbability;
    result.ImportanceSamplingEnabled=false;
  }
  else {
    result.BiasedProbability=Detail::CanonicalNormalizedProbabilities(
        requestedBiasedProbability,physicalRates,
        "biased channel probabilities do not normalize to one");
    for (std::size_t channel=0;channel<physicalRates.size();channel++)
      if (result.BiasedProbability[channel]!=
          result.PhysicalProbability[channel]) {
        result.ImportanceSamplingEnabled=true;
        break;
      }
  }

  result.LikelihoodRatio.assign(physicalRates.size(),0.0);
  result.WeightedChannelVariance.assign(physicalRates.size(),0.0);
  long double expectedSquare=0.0L;
  for (std::size_t channel=0;channel<physicalRates.size();channel++) {
    const double physical=result.PhysicalProbability[channel];
    const double biased=result.BiasedProbability[channel];
    if (!(physical>0.0)) continue;
    const double likelihood=physical/biased;
    if (!(likelihood>0.0) || !std::isfinite(likelihood))
      throw std::runtime_error("channel likelihood ratio is invalid");
    result.LikelihoodRatio[channel]=likelihood;
    const long double q=physical,p=biased;
    long double variance=q*q/p-q*q;
    const long double roundoff=64.0L*
        std::numeric_limits<double>::epsilon()*q*q/p;
    if (variance<0.0L && -variance<=roundoff) variance=0.0L;
    if (variance<0.0L ||
        variance>static_cast<long double>(
            std::numeric_limits<double>::max()))
      throw std::runtime_error(
          "weighted channel variance is invalid");
    result.WeightedChannelVariance[channel]=
        static_cast<double>(variance);
    expectedSquare+=q*q/p;
  }
  if (!(expectedSquare>0.0L) ||
      expectedSquare>static_cast<long double>(
          std::numeric_limits<double>::max()))
    throw std::runtime_error(
        "channel likelihood second moment is invalid");
  result.ExpectedLikelihoodSquare=static_cast<double>(expectedSquare);
  result.AsymptoticEffectiveSampleFraction=
      1.0/result.ExpectedLikelihoodSquare;
  if (!(result.AsymptoticEffectiveSampleFraction>0.0) ||
      result.AsymptoticEffectiveSampleFraction>1.0+
          128.0*std::numeric_limits<double>::epsilon())
    throw std::runtime_error(
        "channel effective-sample fraction is invalid");
  if (result.AsymptoticEffectiveSampleFraction>1.0)
    result.AsymptoticEffectiveSampleFraction=1.0;
  return result;
}

// Select with q* using the same left-to-right subtraction/tie convention as
// the analog path.  The last positive interval is an exact remainder, and the
// fallback protects the endpoint from one final subtraction roundoff.
inline std::size_t SelectBiasedProbabilitySubtract(const cPlan &plan,
    double unitUniform) {
  ValidateUnitUniform(unitUniform);
  if (plan.BiasedProbability.empty())
    throw std::invalid_argument("channel-importance plan is empty");
  double selector=unitUniform;
  std::size_t lastPositive=plan.BiasedProbability.size();
  for (std::size_t channel=0;
      channel<plan.BiasedProbability.size();channel++) {
    const double probability=plan.BiasedProbability[channel];
    if (!(probability>0.0)) continue;
    lastPositive=channel;
    selector-=probability;
    if (selector<=0.0) return channel;
  }
  if (lastPositive==plan.BiasedProbability.size())
    throw std::invalid_argument("channel-importance plan has no active channel");
  return lastPositive;
}

inline double ExpectedEffectiveSampleSize(const cPlan &plan,
    std::size_t sampleCount) {
  if (sampleCount==0 || !(plan.ExpectedLikelihoodSquare>0.0) ||
      !std::isfinite(plan.ExpectedLikelihoodSquare))
    throw std::invalid_argument(
        "expected ESS requires a positive sample count and valid plan");
  const double result=static_cast<double>(sampleCount)/
      plan.ExpectedLikelihoodSquare;
  if (!(result>0.0) || !std::isfinite(result))
    throw std::runtime_error("expected channel ESS is invalid");
  return result;
}

inline double SampleEffectiveSampleSize(double weightSum,
    double weightSquareSum) {
  if (!(weightSum>0.0) || !std::isfinite(weightSum) ||
      !(weightSquareSum>0.0) || !std::isfinite(weightSquareSum))
    throw std::invalid_argument(
        "sample ESS requires positive finite weight moments");
  const double result=weightSum*weightSum/weightSquareSum;
  if (!(result>0.0) || !std::isfinite(result))
    throw std::runtime_error("sample channel ESS is invalid");
  return result;
}

} // namespace ChannelImportancePlanner
} // namespace Comet

#endif

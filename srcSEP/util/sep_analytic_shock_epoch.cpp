#include "sep_analytic_shock_epoch.h"
#include <cmath>
#include <utility>

namespace SEP { namespace Shock {
Transport::Status EvaluateSourceAtEpoch(const TrajectoryConfiguration& configuration,
    double departureEpochS,SourceEpochView* output) {
  if(!output)return Transport::Status::Error(Transport::StatusCode::InvalidArgument,"analytic source epoch view requires an output");
  SourceEpochView candidate;
  // StateAtEpoch uses the existing exact piecewise dr/dt=v(r) law. It writes
  // only this private candidate: even a rejected trajectory cannot partially
  // replace the caller's previously valid immutable source view.
  const auto state=StateAtEpoch(configuration,departureEpochS,&candidate.trajectory);
  if(!state.ok())return state;
  const auto speed=SpeedAtRadius(configuration,candidate.trajectory.radiusM);
  if(!speed.status.ok())return speed.status;
  if(!std::isfinite(candidate.trajectory.radiusM)||candidate.trajectory.radiusM<=0||
      !std::isfinite(speed.value)||speed.value<=0)
    return Transport::Status::Error(Transport::StatusCode::UnresolvedCoefficient,
        "analytic source radius/speed at its declared departure is not finite and positive");
  candidate.speedMPerS=speed.value;
  *output=std::move(candidate);return Transport::Status::Ok();
}
Transport::Status ValidateSourceEpochView(const SourceEpochView& view) {
  SourceEpochView expected;
  const auto status=EvaluateSourceAtEpoch(view.trajectory.configuration,view.trajectory.epochS,&expected);
  if(!status.ok())return status;
  // This is an identity/provenance check on a view made by the deterministic
  // evaluator, not an approximate fit or permission to adjust source geometry.
  if(view.trajectory.radiusM!=expected.trajectory.radiusM||view.speedMPerS!=expected.speedMPerS)
    return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
        "analytic source view mixes an epoch/configuration with a different radius or speed");
  return Transport::Status::Ok();
}
Transport::Status ValidateSourceDepartureView(const SourceEpochView& view,
    const TrajectoryState* preparedOwner,double departureEpochS) {
  if(!preparedOwner)return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
      "analytic source consumption requires an explicitly prepared trajectory owner");
  if(!std::isfinite(departureEpochS))return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
      "analytic source departure must be finite [s]");
  const auto valid=ValidateSourceEpochView(view);if(!valid.ok())return valid;
  SourceEpochView owner;
  const auto evaluated=EvaluateSourceAtEpoch(preparedOwner->configuration,preparedOwner->epochS,&owner);
  if(!evaluated.ok())return evaluated;
  if(owner.trajectory.radiusM!=preparedOwner->radiusM||
      preparedOwner->epochS!=departureEpochS||view.trajectory.epochS!=departureEpochS)
    return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
        "analytic source view/owner do not describe the same prepared departure epoch");
  const auto& expected=preparedOwner->configuration;const auto& actual=view.trajectory.configuration;
  bool same=expected.launchEpochS==actual.launchEpochS&&expected.launchRadiusM==actual.launchRadiusM&&
      expected.knots.size()==actual.knots.size();
  for(std::size_t index=0;same&&index<expected.knots.size();++index)
    same=expected.knots[index].radiusM==actual.knots[index].radiusM&&
        expected.knots[index].speedMPerS==actual.knots[index].speedMPerS;
  if(!same)return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
      "analytic source view belongs to a different governing trajectory configuration");
  return Transport::Status::Ok();
}
} }

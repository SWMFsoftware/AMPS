#include "sep_corona_swcme/surface_shock.h"

#include <cmath>

namespace SEP { namespace CoronaSwcme {

Core::Result<CoronalCME::EllipsoidKinematics> ContactKinematics(
    const EventConfiguration& event,
    const CoronalCME::EllipsoidKinematics& front) {
  using Return=Core::Result<CoronalCME::EllipsoidKinematics>;
  const double fraction=event.regional.contactApexFraction;
  if(!(std::isfinite(fraction)&&fraction>0&&fraction<1))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"contact span fraction is invalid");
  const auto combine=[&](const CoronalCME::KinematicValue& center,
      const CoronalCME::KinematicValue& radial) {
    // Rear + fraction*radial = center - (1-fraction)*radial.
    return CoronalCME::KinematicValue{
        center.value-(1-fraction)*radial.value,
        center.firstDerivative-(1-fraction)*radial.firstDerivative,
        center.secondDerivative-(1-fraction)*radial.secondDerivative};
  };
  const auto scale=[&](const CoronalCME::KinematicValue& value) {
    return CoronalCME::KinematicValue{fraction*value.value,
        fraction*value.firstDerivative,fraction*value.secondDerivative};
  };
  CoronalCME::EllipsoidKinematics contact={
      combine(front.centerDistanceM,front.radialSemiAxisM),
      scale(front.radialSemiAxisM),scale(front.firstLateralSemiAxisM),
      scale(front.secondLateralSemiAxisM)};
  const double dimensions[]={contact.centerDistanceM.value,
      contact.radialSemiAxisM.value,contact.firstLateralSemiAxisM.value,
      contact.secondLateralSemiAxisM.value};
  for(double value:dimensions)if(!(std::isfinite(value)&&value>0))
    return Return::Failure(Core::StatusCode::InvalidState,
        "nested contact has nonpositive geometry");
  return Return::Success(contact);
}

Core::Result<std::shared_ptr<SurfaceShockModel>> SurfaceShockModel::Create(
    std::shared_ptr<const EventConfiguration> event,
    std::shared_ptr<const AmbientModel> ambient) {
  using Return=Core::Result<std::shared_ptr<SurfaceShockModel>>;
  if(!event||!ambient)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "surface/shock model requires event and ambient authorities");
  if(event->physicsFingerprint!=ambient->Event().physicsFingerprint)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "surface and ambient authorities have different event identities");
  std::shared_ptr<SurfaceShockModel> model(new SurfaceShockModel);
  model->event_=std::move(event);
  model->ambient_=std::move(ambient);
  return Return::Success(std::move(model));
}

Core::Result<std::shared_ptr<const SurfaceShockEpoch>> SurfaceShockModel::Prepare(
    double epochS,std::uint64_t backgroundGeneration,
    int polarCells,int azimuthCells) {
  using Return=Core::Result<std::shared_ptr<const SurfaceShockEpoch>>;
  if(backgroundGeneration==0)return Return::Failure(Core::StatusCode::InvalidState,
      "surface/shock background generation zero is reserved");
  const auto evolution=event_->At(epochS);
  if(!evolution.ok())return Return::Failure(evolution.status.code,evolution.status.message);
  const auto contactState=ContactKinematics(*event_,evolution.value.ellipsoid);
  if(!contactState.ok())return Return::Failure(
      contactState.status.code,contactState.status.message);
  const auto front=CoronalCME::FixedOrientationEllipsoid::FromCenter(
      event_->basis,evolution.value.ellipsoid,event_->support.solarRadiusM);
  const auto contact=CoronalCME::FixedOrientationEllipsoid::FromCenter(
      event_->basis,contactState.value,event_->support.solarRadiusM);
  if(!front.ok()||!contact.ok())return Return::Failure(
      Core::StatusCode::InvalidState,!front.ok()?front.status.message:contact.status.message);
  const auto frontPatches=front.value.Tessellate(polarCells,azimuthCells);
  const auto contactPatches=contact.value.Tessellate(polarCells,azimuthCells);
  if(!frontPatches.ok()||!contactPatches.ok())return Return::Failure(
      Core::StatusCode::InvalidState,
      !frontPatches.ok()?frontPatches.status.message:contactPatches.status.message);

  // Same-parameter rear-aligned scaling is analytically nested. Retain an
  // independent direct implicit check so a later geometry change cannot turn
  // the contact into a crossing surface without failing preparation.
  for(const auto& patch:contactPatches.value) {
    const auto inside=front.value.Evaluate(patch.centerM);
    if(!inside.ok()||inside.value.implicitValue>1e-10)return Return::Failure(
        Core::StatusCode::InvalidState,"contact surface is not nested in front");
  }

  std::shared_ptr<SurfaceShockEpoch> candidate(new SurfaceShockEpoch);
  candidate->event=evolution.value;
  candidate->front=front.value;
  candidate->contact=contact.value;
  candidate->contactPatches=contactPatches.value;
  candidate->backgroundGeneration=backgroundGeneration;
  candidate->eventIdentity=event_->physicsFingerprint;
  candidate->frontPatches.reserve(frontPatches.value.size());
  std::vector<CoronalCME::ShockPatchInput> shockInputs;
  shockInputs.reserve(frontPatches.value.size());
  for(const auto& patch:frontPatches.value) {
    const auto surface=front.value.Evaluate(patch.centerM);
    if(!surface.ok())return Return::Failure(surface.status.code,surface.status.message);
    const auto upstream=ambient_->Evaluate(patch.centerM,epochS);
    if(!upstream.ok())return Return::Failure(upstream.status.code,
        "front patch "+std::to_string(patch.physicalId)+": "+upstream.status.message);
    candidate->frontPatches.push_back({patch,surface.value.normalSpeedMPerS,
        upstream.value});
    CoronalCME::ShockPatchInput input;
    input.stableId=patch.physicalId;
    input.parentId=patch.physicalId;
    input.areaM2=patch.areaM2;
    input.outwardNormal=patch.outwardNormal;
    input.upstream={upstream.value.plasma.massDensityKgM3,
        upstream.value.plasma.pressurePa,upstream.value.velocityMPerS,
        upstream.value.magneticFieldT};
    // The scalar speed is positive along the outward (CME-to-ambient) normal.
    // RH therefore sees positive upstream inflow as Vn-U1.n.  The complete
    // parameterized surface velocity is retained separately for the sheath
    // material map; substituting apex speed here would misclassify flanks.
    input.shockNormalSpeedMPerS=surface.value.normalSpeedMPerS;
    // G1/G2 own diagnostic plasma shocks, never particle sources.  Disabling
    // eligibility as well as leaving the rates at zero keeps that separation
    // explicit in every exported patch and aggregate measure.
    input.sourceEnabled=false;
    input.interface={CoronalCME::InterfacePolicy::DiagnosticKinematic,
        backgroundGeneration,true,true,true,true,event_->physicsFingerprint};
    shockInputs.push_back(input);
  }
  CoronalCME::ShockPreparationOptions options;
  options.gammaAdiabatic=event_->composition.gammaAdiabatic;
  options.productionIntent=false;
  options.initialGate=CoronalCME::InitialShockGate::None;
  const auto shocks=shockProvider_.Prepare(epochS,backgroundGeneration,
      shockInputs,options);
  if(!shocks.ok())return Return::Failure(shocks.status.code,shocks.status.message);
  candidate->shocks=shocks.value;
  current_=candidate;
  return Return::Success(std::move(candidate));
}

} } // namespace SEP::CoronaSwcme

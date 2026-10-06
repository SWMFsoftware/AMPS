#include "provider.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace ShockFront { namespace {

constexpr double kPi=3.141592653589793238462643383279502884;

double SmoothStep5(double x) {
  return x*x*x*(10+x*(-15+6*x));
}

double SmoothStep5Derivative(double x) {
  return 30*x*x*(1-x)*(1-x);
}

// Construct a unit vector at polar angle alpha about the event direction.
// The auxiliary basis is selected by its smallest direction component so the
// chart remains well conditioned near every Cartesian pole.
CoronalCME::Vec3 About(CoronalCME::Vec3 direction,double cosineAlpha,
    double azimuth) {
  CoronalCME::Vec3 seed=std::abs(direction.z)<0.8?CoronalCME::Vec3{0,0,1}:
      CoronalCME::Vec3{0,1,0};
  const auto first=CoronalCME::Unit(CoronalCME::Cross(seed,direction));
  const auto second=CoronalCME::Cross(direction,first);
  const double sineAlpha=std::sqrt(std::max(0.0,1-cosineAlpha*cosineAlpha));
  return cosineAlpha*direction+sineAlpha*(std::cos(azimuth)*first+
      std::sin(azimuth)*second);
}

CoronalCME::Vec3 Tangential(CoronalCME::Vec3 value,
    CoronalCME::Vec3 normal) {
  return value-CoronalCME::Dot(value,normal)*normal;
}

double RelativeCrossResidual(CoronalCME::Vec3 velocity,
    CoronalCME::Vec3 field) {
  const double scale=CoronalCME::Norm(velocity)*CoronalCME::Norm(field);
  if(scale==0)return 0;
  return CoronalCME::Norm(CoronalCME::Cross(velocity,field))/scale;
}

} // namespace

Core::Result<LocalJumpEvaluation> EvaluateLocalJump(
    const CoronalCME::MhdPrimitiveState& upstream,
    CoronalCME::Vec3 normal,double normalSpeed,double gamma,
    double weakTolerance,double residualTolerance,
    double magneticDirectionThreshold) {
  using Return=Core::Result<LocalJumpEvaluation>;
  normal=CoronalCME::Unit(normal);
  if(CoronalCME::Norm(normal)==0||!std::isfinite(normalSpeed)||
      !std::isfinite(weakTolerance)||weakTolerance<=0||
      !std::isfinite(residualTolerance)||residualTolerance<=0||
      !std::isfinite(magneticDirectionThreshold)||magneticDirectionThreshold<0)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "local jump evaluation requires a normal, finite speed and positive tolerances");

  LocalJumpEvaluation out;
  const auto characteristics=CoronalCME::EvaluateMhdCharacteristics(
      upstream,normal,gamma);
  if(!characteristics.ok())return Return::Failure(
      characteristics.status.code,characteristics.status.message);
  out.characteristics=characteristics.value;
  out.inflowMPerS=normalSpeed-CoronalCME::Dot(upstream.velocityMPerS,normal);
  out.fastSpeedMarginMPerS=out.inflowMPerS-characteristics.value.fastSpeedMPerS;
  out.fastMach=out.inflowMPerS/characteristics.value.fastSpeedMPerS;
  const double magnetic=CoronalCME::Norm(upstream.magneticFieldT);
  if(magnetic>magneticDirectionThreshold) {
    out.signedMagneticNormalCosine=CoronalCME::Dot(
        upstream.magneticFieldT,normal)/magnetic;
    out.magneticDirectionValid=true;
  }
  if(out.inflowMPerS<=0) {
    out.status=FrontStatus::NonForwardInflow;
    out.reason="canonical ambient overtakes or equals the front normal speed";
    return Return::Success(std::move(out));
  }
  if(out.fastMach<=1) {
    out.status=FrontStatus::SubfastFront;
    out.reason="positive normal inflow is not super-fast";
    return Return::Success(std::move(out));
  }

  const auto jump=CoronalCME::SolveObliqueFastShock(upstream,normal,
      normalSpeed,gamma,residualTolerance);
  if(!jump.ok()) {
    out.status=out.fastMach-1<=weakTolerance?
        FrontStatus::NumericallyUnresolvedWeakShock:FrontStatus::InvalidJump;
    out.reason=jump.status.message;
    return Return::Success(std::move(out));
  }
  out.status=FrontStatus::SolvedFastShock;
  out.downstreamValid=true;
  out.jump=jump.value;

  // Magnetic compression is a norm ratio, not density compression.  It is
  // undefined at a true upstream null even though the hydrodynamic RH limit
  // is valid.  The threshold is used only for diagnostic direction validity;
  // it never modifies B or the jump solve.
  if(magnetic>magneticDirectionThreshold) {
    out.diagnostics.magneticCompressionValid=true;
    out.diagnostics.magneticCompression=
        CoronalCME::Norm(jump.value.downstream.magneticFieldT)/magnetic;
  }

  const double bn=CoronalCME::Dot(upstream.magneticFieldT,normal);
  if(magnetic<=magneticDirectionThreshold) {
    out.diagnostics.htStatus=HtDiagnosticStatus::UpstreamMagneticNull;
    out.diagnostics.reason="HT direction is absent at the upstream magnetic null";
  } else if(bn==0) {
    out.diagnostics.absoluteNormalFieldCosine=0;
    out.diagnostics.htStatus=HtDiagnosticStatus::PerpendicularNoFiniteBoost;
    out.diagnostics.reason=
        "a perpendicular shock with nonzero inflow has no finite conventional HT boost";
  } else {
    // In the shock frame u=U-Vn*n.  The tangential boost v_HT,t is chosen so
    // u-v_HT,t is parallel to B on the upstream side.  Ideal RH tangential-E
    // continuity then cancels E on the downstream side with the same boost.
    const CoronalCME::Vec3 u1=upstream.velocityMPerS-normalSpeed*normal;
    const CoronalCME::Vec3 u2=jump.value.downstream.velocityMPerS-
        normalSpeed*normal;
    const double u1n=CoronalCME::Dot(u1,normal);
    out.diagnostics.htTangentialBoostMPerS=Tangential(u1,normal)-
        (u1n/bn)*Tangential(upstream.magneticFieldT,normal);
    out.diagnostics.upstreamVelocityHtMPerS=
        u1-out.diagnostics.htTangentialBoostMPerS;
    out.diagnostics.downstreamVelocityHtMPerS=
        u2-out.diagnostics.htTangentialBoostMPerS;
    out.diagnostics.incidentHtSpeedMPerS=
        CoronalCME::Norm(out.diagnostics.upstreamVelocityHtMPerS);
    out.diagnostics.absoluteNormalFieldCosine=std::abs(bn)/magnetic;
    out.diagnostics.upstreamElectricCancellation=RelativeCrossResidual(
        out.diagnostics.upstreamVelocityHtMPerS,upstream.magneticFieldT);
    out.diagnostics.downstreamElectricCancellation=RelativeCrossResidual(
        out.diagnostics.downstreamVelocityHtMPerS,
        jump.value.downstream.magneticFieldT);
    const bool finite=std::isfinite(out.diagnostics.incidentHtSpeedMPerS)&&
        std::isfinite(out.diagnostics.upstreamElectricCancellation)&&
        std::isfinite(out.diagnostics.downstreamElectricCancellation);
    out.diagnostics.htStatus=finite?HtDiagnosticStatus::Valid:
        HtDiagnosticStatus::NumericalFailure;
    if(!finite)out.diagnostics.reason=
        "finite nonzero Bn produced a nonfinite HT reconstruction";
  }
  return Return::Success(std::move(out));
}

Core::Result<GeometrySample> EvaluateSseRay(const SseKinematicState& state,
    CoronalCME::Vec3 ray,std::uint64_t stableId) {
  using Return=Core::Result<GeometrySample>;
  const double directionNorm=CoronalCME::Norm(state.direction);
  const double rayNorm=CoronalCME::Norm(ray);
  if(!std::isfinite(directionNorm)||!std::isfinite(rayNorm)||
      std::abs(directionNorm-1)>1e-12||rayNorm<=0||
      !std::isfinite(state.apexRadiusM)||state.apexRadiusM<=0||
      !std::isfinite(state.apexSpeedMPerS)||
      !std::isfinite(state.halfWidthRad)||state.halfWidthRad<=0||
      state.halfWidthRad>=0.5*kPi||
      !std::isfinite(state.halfWidthRateRadPerS)||
      !std::isfinite(state.directionRatePerS.x)||
      !std::isfinite(state.directionRatePerS.y)||
      !std::isfinite(state.directionRatePerS.z)||
      std::abs(CoronalCME::Dot(state.direction,state.directionRatePerS))>
        1e-12*std::max(1.0,CoronalCME::Norm(state.directionRatePerS)))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "finite SSE requires positive regular geometry, a unit direction and d dot d_dot=0");
  ray=ray/rayNorm;
  const double cosine=CoronalCME::Dot(ray,state.direction);
  const double minimum=std::cos(state.halfWidthRad);
  if(cosine<minimum-64*std::numeric_limits<double>::epsilon())
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "ray lies outside finite SSE support");
  const double sine=std::sin(state.halfWidthRad);
  const double cosineWidth=std::cos(state.halfWidthRad);
  const double center=state.apexRadiusM/(1+sine);
  const double radius=center*sine;
  const double transverse2=std::max(0.0,1-cosine*cosine);
  double radicand=radius*radius-center*center*transverse2;
  const double scale=std::max(radius*radius,center*center);
  if(radicand<-128*std::numeric_limits<double>::epsilon()*scale)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "ray has no finite SSE leading intersection");
  // The support edge is a double ray/sphere root.  Independently rounded
  // sin/cos calls can leave a tiny *positive* discriminant as well as a tiny
  // negative one; retaining its square root produces an O(sqrt(epsilon))
  // spurious edge normal speed.  Collapse only the roundoff-sized symmetric
  // interval.  A genuinely positive near-edge root is left untouched.
  if(std::abs(radicand)<=128*std::numeric_limits<double>::epsilon()*scale)
    radicand=0.0;
  else radicand=std::max(0.0,radicand);
  GeometrySample out;out.stableId=stableId;
  out.positionM=(center*cosine+std::sqrt(radicand))*ray;
  out.outwardNormal=CoronalCME::Unit(out.positionM-center*state.direction);
  const double directionNormal=CoronalCME::Dot(state.direction,out.outwardNormal);
  out.normalSpeedMPerS=state.apexSpeedMPerS*(directionNormal+sine)/(1+sine)+
      state.apexRadiusM*cosineWidth*state.halfWidthRateRadPerS*
        (1-directionNormal)/((1+sine)*(1+sine))+
      center*CoronalCME::Dot(state.directionRatePerS,out.outwardNormal);
  out.supportEdge=std::abs(cosine-minimum)<=1e-10;
  return Return::Success(out);
}

Core::Result<TrajectoryState> EvaluateQuadraticDrag(double epoch,
    double handoffTime,double handoffRadius,double handoffSpeed,double wind,
    double gamma) {
  using Return=Core::Result<TrajectoryState>;
  if(!std::isfinite(epoch)||!std::isfinite(handoffTime)||epoch<handoffTime||
      !std::isfinite(handoffRadius)||handoffRadius<=0||
      !std::isfinite(handoffSpeed)||!std::isfinite(wind)||
      !std::isfinite(gamma)||gamma<0)return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
          "quadratic drag requires a finite post-handoff state and gamma>=0");
  TrajectoryState out;out.timeS=epoch;out.phase=Phase::SwcmeOuter;
  const double dt=epoch-handoffTime;const double z0=handoffSpeed-wind;
  if(gamma==0||z0==0) {
    out.apexSpeedMPerS=handoffSpeed;
    out.apexRadiusM=handoffRadius+handoffSpeed*dt;
    out.apexAccelerationMPerS2=0;
  } else {
    const double denominator=1+gamma*std::abs(z0)*dt;
    out.apexSpeedMPerS=wind+z0/denominator;
    out.apexRadiusM=handoffRadius+wind*dt+
        std::copysign(std::log1p(gamma*std::abs(z0)*dt),z0)/gamma;
    out.apexAccelerationMPerS2=-gamma*(out.apexSpeedMPerS-wind)*
        std::abs(out.apexSpeedMPerS-wind);
  }
  return Return::Success(out);
}

AreaLedger SummarizeAreas(const std::vector<ShockRecord>& records) noexcept {
  AreaLedger out;
  for(const auto& record:records) {
    const double area=std::isfinite(record.geometry.areaM2)&&
        record.geometry.areaM2>0?record.geometry.areaM2:0;
    out.geometricSupportM2+=area;
    if(record.status==FrontStatus::BelowPhysicalInnerBoundary)
      out.belowInnerBoundaryM2+=area;
    if(record.fastMach>1&&(record.status==FrontStatus::SolvedFastShock||
        record.status==FrontStatus::NumericallyUnresolvedWeakShock||
        record.status==FrontStatus::WrongBranch||
        record.status==FrontStatus::InvalidJump))out.superfastCandidateM2+=area;
    if(record.status==FrontStatus::SolvedFastShock)out.acceptedShockM2+=area;
    if(record.status==FrontStatus::AmbientUnavailable||
        record.status==FrontStatus::NumericallyUnresolvedWeakShock||
        record.status==FrontStatus::WrongBranch||
        record.status==FrontStatus::InvalidJump)out.numericalFailureM2+=area;
  }
  return out;
}

Core::Result<std::shared_ptr<Provider>> Provider::Create(
    std::shared_ptr<const Configuration> configuration) {
  using Return=Core::Result<std::shared_ptr<Provider>>;
  if(!configuration)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "shock-front provider requires a resolved event");
  const auto ambient=AmbientModel::Create(configuration->ambient);
  if(!ambient.ok())return Return::Failure(ambient.status.code,ambient.status.message);
  std::shared_ptr<Provider> provider(new Provider);
  provider->configuration_=std::move(configuration);
  provider->ambient_=ambient.value;
  const auto handoff=provider->HandoffTimeS();
  if(!handoff.ok())return Return::Failure(handoff.status.code,handoff.status.message);
  const auto endpoint=provider->EndpointTimeS();
  if(!endpoint.ok())return Return::Failure(endpoint.status.code,endpoint.status.message);
  if(endpoint.value>provider->configuration_->ambient.support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "resolved endpoint lies after declared event coverage");
  return Return::Success(std::move(provider));
}

Core::Result<double> Provider::HandoffTimeS() const {
  using Return=Core::Result<double>;
  const auto& c=*configuration_;
  const auto radius=[&](double t) {
    if(c.historyModel==HistoryModel::ConstantSpeed)
      return c.initialApexRadiusM+c.initialApexSpeedMPerS*t;
    if(t>=c.accelerationDurationS)return c.initialApexRadiusM+
        0.5*c.accelerationDurationS*(c.initialApexSpeedMPerS+
        c.finalPulseSpeedMPerS)+c.finalPulseSpeedMPerS*
        (t-c.accelerationDurationS);
    const double s=t/c.accelerationDurationS;
    return c.initialApexRadiusM+c.initialApexSpeedMPerS*t+
        (c.finalPulseSpeedMPerS-c.initialApexSpeedMPerS)*
        c.accelerationDurationS*(2.5*std::pow(s,4)-3*std::pow(s,5)+std::pow(s,6));
  };
  if(radius(c.ambient.support.startS)>c.handoffApexRadiusM||
      radius(c.ambient.support.endS)<c.handoffApexRadiusM)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "front history does not bracket the handoff apex radius");
  double low=c.ambient.support.startS,high=c.ambient.support.endS;
  for(int i=0;i<100&&high-low>1e-9;++i) {
    const double mid=0.5*(low+high);
    if(radius(mid)<c.handoffApexRadiusM)low=mid;else high=mid;
  }
  return Return::Success(0.5*(low+high));
}

Core::Result<HandoffReceipt> Provider::Handoff() const {
  using Return=Core::Result<HandoffReceipt>;
  const auto root=HandoffTimeS();if(!root.ok())return Return::Failure(
      root.status.code,root.status.message);
  const auto matched=Trajectory(root.value);if(!matched.ok())return Return::Failure(
      matched.status.code,matched.status.message);
  const auto outer=EvaluateQuadraticDrag(root.value,root.value,
      configuration_->handoffApexRadiusM,configuration_->finalPulseSpeedMPerS,
      configuration_->effectiveTrajectoryWindMPerS,configuration_->dragGammaPerM);
  if(!outer.ok())return Return::Failure(outer.status.code,outer.status.message);
  HandoffReceipt receipt;receipt.timeS=root.value;
  receipt.apexRadiusM=configuration_->handoffApexRadiusM;
  receipt.apexSpeedMPerS=configuration_->finalPulseSpeedMPerS;
  receipt.coronalAccelerationMPerS2=matched.value.apexAccelerationMPerS2;
  receipt.outerAccelerationMPerS2=-configuration_->dragGammaPerM*
      (configuration_->finalPulseSpeedMPerS-
       configuration_->effectiveTrajectoryWindMPerS)*
      std::abs(configuration_->finalPulseSpeedMPerS-
       configuration_->effectiveTrajectoryWindMPerS);
  receipt.maximumSurfacePositionMismatchM=
      std::abs(matched.value.apexRadiusM-outer.value.apexRadiusM);
  receipt.maximumNormalSpeedMismatchMPerS=
      std::abs(matched.value.apexSpeedMPerS-outer.value.apexSpeedMPerS);
  receipt.c1Matched=receipt.maximumSurfacePositionMismatchM<=1e-6&&
      receipt.maximumNormalSpeedMismatchMPerS<=1e-9;
  receipt.c2Matched=std::abs(receipt.coronalAccelerationMPerS2-
      receipt.outerAccelerationMPerS2)<=1e-9;
  receipt.authorityReset=false;
  receipt.effectiveTrajectoryWindMPerS=configuration_->effectiveTrajectoryWindMPerS;
  const auto upstream=ambient_->Evaluate(
      configuration_->handoffApexRadiusM*configuration_->direction,root.value);
  if(upstream.ok()) {
    receipt.canonicalWindValid=true;
    receipt.canonicalApexWindMPerS=CoronalCME::Dot(
        upstream.value.velocityMPerS,configuration_->direction);
    receipt.effectiveMinusCanonicalWindMPerS=
        receipt.effectiveTrajectoryWindMPerS-receipt.canonicalApexWindMPerS;
  }
  return Return::Success(receipt);
}

Core::Result<TrajectoryState> Provider::Trajectory(double t) const {
  using Return=Core::Result<TrajectoryState>;
  const auto& c=*configuration_;
  if(!std::isfinite(t)||t<c.ambient.support.startS||t>c.ambient.support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "front epoch is outside declared event coverage");
  const auto handoff=HandoffTimeS();
  if(!handoff.ok())return Return::Failure(handoff.status.code,handoff.status.message);
  TrajectoryState out;out.timeS=t;
  if(t<=handoff.value) {
    out.phase=Phase::CoronalHistory;
    if(c.historyModel==HistoryModel::ConstantSpeed) {
      out.apexRadiusM=c.initialApexRadiusM+c.initialApexSpeedMPerS*t;
      out.apexSpeedMPerS=c.initialApexSpeedMPerS;out.apexAccelerationMPerS2=0;
    } else if(t<c.accelerationDurationS) {
      const double s=t/c.accelerationDurationS;
      const double delta=c.finalPulseSpeedMPerS-c.initialApexSpeedMPerS;
      out.apexSpeedMPerS=c.initialApexSpeedMPerS+delta*SmoothStep5(s);
      out.apexRadiusM=c.initialApexRadiusM+c.initialApexSpeedMPerS*t+
          delta*c.accelerationDurationS*(2.5*std::pow(s,4)-3*std::pow(s,5)+
          std::pow(s,6));
      out.apexAccelerationMPerS2=delta/c.accelerationDurationS*
          SmoothStep5Derivative(s);
    } else {
      const double pulseRadius=c.initialApexRadiusM+0.5*c.accelerationDurationS*
          (c.initialApexSpeedMPerS+c.finalPulseSpeedMPerS);
      out.apexRadiusM=pulseRadius+c.finalPulseSpeedMPerS*
          (t-c.accelerationDurationS);
      out.apexSpeedMPerS=c.finalPulseSpeedMPerS;out.apexAccelerationMPerS2=0;
    }
    // Root finding can return the crossing a few ulps to either side.  The
    // matched state is authoritative, so report its exact configured radius.
    if(std::abs(t-handoff.value)<=2e-9)out.apexRadiusM=c.handoffApexRadiusM;
    return Return::Success(out);
  }
  return EvaluateQuadraticDrag(t,handoff.value,c.handoffApexRadiusM,
      c.finalPulseSpeedMPerS,c.effectiveTrajectoryWindMPerS,c.dragGammaPerM);
}

Core::Result<double> Provider::EndpointTimeS() const {
  using Return=Core::Result<double>;
  const auto& c=*configuration_;
  double low=c.ambient.support.startS,high=c.ambient.support.endS;
  const auto lower=Trajectory(low);const auto upper=Trajectory(high);
  if(!lower.ok()||!upper.ok()||lower.value.apexRadiusM>c.endpointRadiusM||
      upper.value.apexRadiusM<c.endpointRadiusM)return Return::Failure(
          Core::StatusCode::OutOfDomain,
          "front trajectory does not bracket the requested endpoint");
  for(int i=0;i<120&&high-low>1e-8;++i) {
    const double mid=0.5*(low+high);const auto state=Trajectory(mid);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    if(state.value.apexRadiusM<c.endpointRadiusM)low=mid;else high=mid;
  }
  return Return::Success(0.5*(low+high));
}

Core::Result<GeometrySample> Provider::EvaluateRay(CoronalCME::Vec3 ray,
    double epoch,std::uint64_t stableId) const {
  using Return=Core::Result<GeometrySample>;
  const auto state=Trajectory(epoch);
  if(!state.ok())return Return::Failure(state.status.code,state.status.message);
  return EvaluateSseRay({state.value.apexRadiusM,state.value.apexSpeedMPerS,
      configuration_->direction,{},configuration_->halfWidthRad,0},ray,stableId);
}

Core::Result<ShockRecord> Provider::EvaluateFrontPoint(
    CoronalCME::Vec3 position,double epoch,std::uint64_t stableId) const {
  using Return=Core::Result<ShockRecord>;
  const double radius=CoronalCME::Norm(position);
  if(!std::isfinite(radius)||radius<=0)return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "front-point query requires a finite nonzero HCI position");
  const auto geometry=EvaluateRay(position/radius,epoch,stableId);
  if(!geometry.ok())return Return::Failure(geometry.status.code,geometry.status.message);
  const double mismatch=CoronalCME::Norm(geometry.value.positionM-position);
  const double tolerance=256*std::numeric_limits<double>::epsilon()*
      std::max(1.0,radius);
  if(mismatch>tolerance)return Return::Failure(Core::StatusCode::OutOfDomain,
      "point is not on the supported leading front");
  ShockRecord out;out.geometry=geometry.value;
  if(radius<configuration_->physicalInnerRadiusM) {
    out.status=FrontStatus::BelowPhysicalInnerBoundary;
    out.reason="supported mathematical front is below the physical ambient boundary";
    return Return::Success(std::move(out));
  }
  const auto upstream=ambient_->Evaluate(position,epoch);
  if(!upstream.ok())return Return::Failure(upstream.status.code,upstream.status.message);
  out.upstream=upstream.value;
  const CoronalCME::MhdPrimitiveState primitive{
      upstream.value.plasma.massDensityKgM3,upstream.value.plasma.pressurePa,
      upstream.value.velocityMPerS,upstream.value.magneticFieldT};
  const auto local=EvaluateLocalJump(primitive,geometry.value.outwardNormal,
      geometry.value.normalSpeedMPerS,
      configuration_->ambient.composition.gammaAdiabatic,
      configuration_->weakMachTolerance,configuration_->rhResidualTolerance,
      configuration_->ambient.ambient.minimumMagneticFieldT);
  if(!local.ok())return Return::Failure(local.status.code,local.status.message);
  out.characteristics=local.value.characteristics;out.inflowMPerS=local.value.inflowMPerS;
  out.fastSpeedMarginMPerS=local.value.fastSpeedMarginMPerS;
  out.fastMach=local.value.fastMach;
  out.signedMagneticNormalCosine=local.value.signedMagneticNormalCosine;
  out.magneticDirectionValid=local.value.magneticDirectionValid;
  out.downstreamValid=local.value.downstreamValid;out.jump=local.value.jump;
  out.diagnostics=local.value.diagnostics;out.status=local.value.status;
  out.reason=local.value.reason;return Return::Success(std::move(out));
}

Core::Result<AmbientState> Provider::QueryAmbient(CoronalCME::Vec3 position,
    double epoch,std::uint64_t generation) const {
  return ambient_->EvaluateWithDerivatives(position,epoch,generation);
}

Core::Status Provider::RequireCapability(VolumeCapability capability) const {
  if(capability==VolumeCapability::AmbientReference)return Core::Status::Success();
  return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
      "reduced shock-front mode has no physical downstream/sheath/ejecta volume");
}

Core::Result<std::shared_ptr<const Epoch>> Provider::Prepare(
    double epoch,std::uint64_t generation) {
  using Return=Core::Result<std::shared_ptr<const Epoch>>;
  if(generation==0)return Return::Failure(Core::StatusCode::InvalidState,
      "shock-front generation zero is reserved");
  const auto trajectory=Trajectory(epoch);
  if(!trajectory.ok())return Return::Failure(trajectory.status.code,trajectory.status.message);
  std::shared_ptr<Epoch> candidate(new Epoch);candidate->trajectory=trajectory.value;
  candidate->generation=generation;candidate->ambientGeneration=generation;
  candidate->eventIdentity=configuration_->physicsFingerprint;
  const auto handoff=Handoff();if(!handoff.ok())return Return::Failure(
      handoff.status.code,handoff.status.message);
  candidate->handoff=handoff.value;

  const double sine=std::sin(configuration_->halfWidthRad);
  const double center=trajectory.value.apexRadiusM/(1+sine);
  const double sphereRadius=center*sine;
  // Uniform bins in mu=n.d give equal exact generating-sphere area.  Samples
  // are cell centroids; area is not inferred from heliocentric radial angles,
  // whose metric is singular at the tangent support edge.
  const double muMinimum=-sine;
  const double dmu=(1-muMinimum)/configuration_->polarCells;
  const double dphi=2*kPi/configuration_->azimuthCells;
  const double area=sphereRadius*sphereRadius*dmu*dphi;
  candidate->records.reserve(static_cast<std::size_t>(configuration_->polarCells)*
      configuration_->azimuthCells);
  for(int i=0;i<configuration_->polarCells;++i)for(int j=0;
      j<configuration_->azimuthCells;++j) {
    const double mu=muMinimum+(i+0.5)*dmu;
    const double phi=(j+0.5)*dphi;
    const auto normal=About(configuration_->direction,mu,phi);
    ShockRecord record;record.geometry.stableId=
        static_cast<std::uint64_t>(i)*configuration_->azimuthCells+j+1;
    record.geometry.outwardNormal=normal;
    record.geometry.positionM=center*configuration_->direction+sphereRadius*normal;
    record.geometry.normalSpeedMPerS=trajectory.value.apexSpeedMPerS*(mu+sine)/(1+sine);
    record.geometry.areaM2=area;
    record.geometry.supportEdge=i==0;
    candidate->area.geometricSupportM2+=area;
    if(CoronalCME::Norm(record.geometry.positionM)<configuration_->physicalInnerRadiusM) {
      record.status=FrontStatus::BelowPhysicalInnerBoundary;
      record.reason="supported mathematical front is below the physical ambient boundary";
      candidate->area.belowInnerBoundaryM2+=area;
      candidate->records.push_back(std::move(record));continue;
    }
    const auto upstream=ambient_->Evaluate(record.geometry.positionM,epoch);
    if(!upstream.ok()) {
      record.status=FrontStatus::AmbientUnavailable;record.reason=upstream.status.message;
      candidate->area.numericalFailureM2+=area;
      candidate->records.push_back(std::move(record));continue;
    }
    record.upstream=upstream.value;
    CoronalCME::MhdPrimitiveState primitive{upstream.value.plasma.massDensityKgM3,
        upstream.value.plasma.pressurePa,upstream.value.velocityMPerS,
        upstream.value.magneticFieldT};
    const auto local=EvaluateLocalJump(primitive,normal,
        record.geometry.normalSpeedMPerS,
        configuration_->ambient.composition.gammaAdiabatic,
        configuration_->weakMachTolerance,configuration_->rhResidualTolerance,
        configuration_->ambient.ambient.minimumMagneticFieldT);
    if(!local.ok()) {
      record.status=FrontStatus::InvalidJump;record.reason=local.status.message;
      candidate->area.numericalFailureM2+=area;
      candidate->records.push_back(std::move(record));continue;
    }
    record.characteristics=local.value.characteristics;
    record.inflowMPerS=local.value.inflowMPerS;
    record.fastSpeedMarginMPerS=local.value.fastSpeedMarginMPerS;
    record.fastMach=local.value.fastMach;
    record.signedMagneticNormalCosine=local.value.signedMagneticNormalCosine;
    record.magneticDirectionValid=local.value.magneticDirectionValid;
    record.downstreamValid=local.value.downstreamValid;
    record.jump=local.value.jump;record.diagnostics=local.value.diagnostics;
    record.status=local.value.status;record.reason=local.value.reason;
    candidate->records.push_back(std::move(record));
  }
  // Recompute from immutable records rather than relying on branch-local
  // increments.  This is the single accounting authority used by output and
  // by disconnected/unknown-area tests; physical no-shock area remains in the
  // geometric denominator and is never confused with numerical unknown.
  candidate->area=SummarizeAreas(candidate->records);
  const auto endpoint=EndpointTimeS();
  candidate->geometricEndpointReached=endpoint.ok()&&epoch>=endpoint.value-1e-7;
  // The apex is stable sample 1 only in a separate exact query; use the
  // closest normal-centroid ring here and keep this aggregate explicitly an
  // approximation for receipts. Observer qualification evaluates its root.
  const auto best=std::max_element(candidate->records.begin(),candidate->records.end(),
      [&](const ShockRecord& a,const ShockRecord& b){return
        CoronalCME::Dot(a.geometry.outwardNormal,configuration_->direction)<
        CoronalCME::Dot(b.geometry.outwardNormal,configuration_->direction);});
  candidate->apexShockAccepted=best!=candidate->records.end()&&
      best->status==FrontStatus::SolvedFastShock;
  if(configuration_->validityPolicy==TrajectoryValidityPolicy::RequireDeclaredShockCoverage&&
      (candidate->area.numericalFailureM2>0||
       (configuration_->requireFastShockAtObserver&&
        candidate->geometricEndpointReached&&!candidate->apexShockAccepted)))
    return Return::Failure(Core::StatusCode::InvalidState,
        "candidate violates declared numerical/shock coverage; committed epoch retained");
  current_=candidate;
  return Return::Success(std::move(candidate));
}

} } } // namespace SEP::CoronaSwcme::ShockFront

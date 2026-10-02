#include "transport_adapter.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Adapters {
namespace {

bool Finite(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

double ResolveKappaPerpendicular(const MoverInput& input) {
  switch (input.perpendicularDiffusion) {
    case RuntimeModel::PerpendicularDiffusionMode::None: return 0.0;
    case RuntimeModel::PerpendicularDiffusionMode::Constant:
      return input.constantKappaPerpendicularM2PerS;
    case RuntimeModel::PerpendicularDiffusionMode::ConstantRatio:
      return input.kappaPerpendicularToParallelRatio *
          input.local.kappaParallelM2PerS;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

Core::Vec3 GradientOfMagnitude(const Background::BackgroundSample& background) {
  // gradB(i,j)=d B_i/d x_j. Therefore d|B|/d x_j=sum_i b_i gradB(i,j).
  return {
      background.bHat.x*background.gradB(0,0) +
          background.bHat.y*background.gradB(1,0) +
          background.bHat.z*background.gradB(2,0),
      background.bHat.x*background.gradB(0,1) +
          background.bHat.y*background.gradB(1,1) +
          background.bHat.z*background.gradB(2,1),
      background.bHat.x*background.gradB(0,2) +
          background.bHat.y*background.gradB(1,2) +
          background.bHat.z*background.gradB(2,2)};
}

MoverResult Failed(const MoverInput& input, const Core::Status& status) {
  MoverResult result;
  result.particle = input.particle;
  result.disposition = ParticleDisposition::Failed;
  result.status = status;
  return result;
}

}  // namespace

const char* Name(ParticleDisposition disposition) {
  switch (disposition) {
    case ParticleDisposition::Active: return "active";
    case ParticleDisposition::Escaped: return "escaped";
    case ParticleDisposition::Absorbed: return "absorbed";
    case ParticleDisposition::Failed: return "failed";
  }
  return "unknown";
}

ShockIntersection FirstShockIntersection(
    const Core::Vec3& initialPositionM,
    const Core::Vec3& finalPositionM,
    double dtS, const ExpandingShock& shock,
    std::uint64_t lastShockGeneration) {
  ShockIntersection result;
  result.generation=shock.generation;
  if (!shock.active||shock.generation==0||shock.generation==lastShockGeneration) {
    result.status=Core::Status::OK(); return result;
  }
  result.status=ValidateShockGeometry(shock);
  if (!result.status.ok()) return result;
  if (!Finite(initialPositionM)||!Finite(finalPositionM)||
      !std::isfinite(dtS)||dtS<=0||
      !std::isfinite(shock.radiusAtStepStartM+shock.radialSpeedMPerS*dtS)||
      shock.radiusAtStepStartM+shock.radialSpeedMPerS*dtS<=0) {
    result.status=Invalid("expanding-shock intersection input is invalid"); return result;
  }
  const auto sphere=ShockGeneratingSphere(shock);
  const auto segment=finalPositionM-initialPositionM;
  const auto relative=initialPositionM-sphere.centerM;
  const auto relativeSegment=segment-sphere.centerVelocityMPerS*dtS;
  const double radiusChange=sphere.radiusSpeedMPerS*dtS;
  // |p0-C0+f*(dp-dC)|^2=(a0+f*da)^2. Scale lengths before forming
  // coefficients so heliospheric distances do not overflow the discriminant,
  // and use the cancellation-resistant q formula for a root near an endpoint.
  const double length=std::max({1.0,relative.Norm(),relativeSegment.Norm(),
      sphere.radiusM,std::fabs(radiusChange)});
  if (!std::isfinite(length)) {
    result.status=Invalid("shock-relative segment is not representable");return result;
  }
  const auto r=relative/length, d=relativeSegment/length;
  const double radius=sphere.radiusM/length, change=radiusChange/length;
  const double a=d.Dot(d)-change*change;
  const double b=2*(r.Dot(d)-radius*change);
  const double c=r.Dot(r)-radius*radius;
  const double eps=64*std::numeric_limits<double>::epsilon();
  double roots[2]={std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
  if (std::fabs(a)<=eps*std::max({std::fabs(b),std::fabs(c),1e-300})) {
    if (b!=0) roots[0]=-c/b;
  } else {
    double discriminant=b*b-4*a*c;
    const double allowance=eps*(b*b+std::fabs(4*a*c));
    if (discriminant>=-allowance) {
      discriminant=std::max(0.0,discriminant);
      const double q=-0.5*(b+std::copysign(std::sqrt(discriminant),b));
      if (q==0) roots[0]=-b/(2*a);
      else { roots[0]=q/a; roots[1]=c/q; }
    }
  }
  std::sort(roots,roots+2);
  for (double f:roots) {
    if (!std::isfinite(f)||f<-eps||f>1+eps) continue;
    f=std::max(0.0,std::min(1.0,f));
    auto at=sphere;
    at.centerM+=sphere.centerVelocityMPerS*(dtS*f);
    at.radiusM+=radiusChange*f;
    const auto point=initialPositionM+segment*f;
    // Checking only the cone would also accept the rear sphere intersection.
    // Reject it, then continue to the second quadratic root if it is physical.
    if (!OnOutwardShockCap(shock,point,at)) continue;
    result.crossed=true; result.stepFraction=f; result.positionM=point;
    result.outwardNormal=(point-at.centerM).Normalized();
    result.normalSpeedMPerS=sphere.centerVelocityMPerS.Dot(result.outwardNormal)+
        sphere.radiusSpeedMPerS;
    break;
  }
  result.status=Core::Status::OK(); return result;
}

const std::vector<std::string>& ProductionMoverRegistry::CanonicalNames() {
  static const std::vector<std::string> names = {
      "parker", "focused-diffusion", "focused-scattering"};
  return names;
}

Core::Status ProductionMoverRegistry::Validate(
    RuntimeModel::TransportModel selected) {
  switch (selected) {
    case RuntimeModel::TransportModel::Parker3D:
    case RuntimeModel::TransportModel::FocusedDiffusion3D:
    case RuntimeModel::TransportModel::FocusedScattering3D:
      return Core::Status::OK();
  }
  return Invalid("transport selection is not in the production mover registry");
}

MoverResult AdvanceParticle(const MoverInput& input) {
  const Core::Status registered = ProductionMoverRegistry::Validate(input.model);
  if (!registered.ok()) return Failed(input, registered);
  const ParticleRecord& particle = input.particle;
  const Background::BackgroundSample& background = input.local.background;
  if (particle.stableId == 0 || particle.species < 0 ||
      !Finite(particle.positionM) ||
      !std::isfinite(particle.momentumKgMPerS) ||
      particle.momentumKgMPerS < 0.0 || !std::isfinite(particle.mu) ||
      particle.mu < -1.0 || particle.mu > 1.0 ||
      !std::isfinite(particle.statisticalWeight) ||
      particle.statisticalWeight <= 0.0 ||
      !std::isfinite(input.speciesMassKg) || input.speciesMassKg <= 0.0 ||
      !std::isfinite(input.speciesChargeC) ||
      (input.drift != RuntimeModel::DriftMode::None &&
       input.speciesChargeC == 0.0) ||
      !std::isfinite(input.requestedDtS) || input.requestedDtS <= 0.0 ||
      !std::isfinite(input.innerRadiusM) || input.innerRadiusM <= 0.0 ||
      !std::isfinite(input.outerRadiusM) ||
      input.outerRadiusM <= input.innerRadiusM || input.campaignSeed == 0) {
    return Failed(input, Invalid("particle mover input record is invalid"));
  }
  if (!background.status.ok() || !background.valid ||
      !Finite(background.B) || !Finite(background.U) ||
      !std::isfinite(background.absB) || background.absB <= 0.0 ||
      std::fabs(background.bHat.Norm() - 1.0) > 1.0e-12) {
    return Failed(input, Core::Status(
        Core::StatusCode::BackgroundInvalid,
        "particle mover received an invalid background sample"));
  }

  const double speed = Transport::RelativisticSpeed(
      particle.momentumKgMPerS, input.speciesMassKg);
  const double kappaPerpendicularM2PerS = ResolveKappaPerpendicular(input);
  if (!std::isfinite(kappaPerpendicularM2PerS) ||
      kappaPerpendicularM2PerS < 0.0)
    return Failed(input, Invalid("resolved perpendicular diffusion is invalid"));
  Core::Vec3 guidingCenterDrift;
  Transport::GuidingCenterInput driftInput;
  driftInput.bHat = background.bHat;
  driftInput.gradAbsBTPerM = GradientOfMagnitude(background);
  driftInput.curvaturePerM = background.curvature;
  driftInput.absBT = background.absB;
  driftInput.momentumKgMPerS = particle.momentumKgMPerS;
  driftInput.speedMPerS = speed;
  driftInput.pitchCosine = particle.mu;
  driftInput.chargeC = input.speciesChargeC;
  driftInput.pitchAveraged =
      input.model == RuntimeModel::TransportModel::Parker3D;
  driftInput.includeGradientB = input.drift == RuntimeModel::DriftMode::GradientB ||
      input.drift == RuntimeModel::DriftMode::GradientAndCurvature;
  driftInput.includeCurvature = input.drift == RuntimeModel::DriftMode::Curvature ||
      input.drift == RuntimeModel::DriftMode::GradientAndCurvature;
  Core::Status driftStatus =
      Transport::EvaluateGuidingCenterDrift(driftInput, &guidingCenterDrift);
  if (!driftStatus.ok()) return Failed(input, driftStatus);
  Transport::TimeStepPhysics physics;
  physics.requestedS = input.requestedDtS;
  physics.cellSizeM = input.local.cellSizeM;
  physics.characteristicSpeedMPerS = background.U.Norm() + speed;
  physics.kappaParallelM2PerS = input.model == RuntimeModel::TransportModel::Parker3D
      ? input.local.kappaParallelM2PerS : 0.0;
  physics.kappaPerpendicularM2PerS = kappaPerpendicularM2PerS;
  physics.focusingRatePerS = input.model != RuntimeModel::TransportModel::Parker3D
      ? std::fabs(Transport::FocusedPitchDriftPerS(
            particle.mu, speed, [&]() {
              Transport::FocusedLocalState local;
              local.bHat = background.bHat;
              local.bulkVelocityMPerS = background.U;
              local.divBhatPerM = background.divBhat;
              local.divUPerS = background.divU;
              local.fieldAlignedStrainPerS = background.fieldAlignedStrain;
              return local;
            }())) : 0.0;
  physics.coolingRatePerS = std::fabs(background.divU) / 3.0;
  physics.fractionalFieldVariationPerS =
      input.local.fractionalFieldVariationPerS;
  physics.timeToSnapshotBoundaryS = input.local.timeToSnapshotBoundaryS;
  // Distance to the full generating sphere is a lower bound on distance to
  // its finite outward cap. Using it may subcycle early outside the SSE cone,
  // but never uses a fabricated Sun-centered flank. Center translation plus
  // radius growth is bounded by |V_apex| for SSE, as for the spherical case.
  if (input.shock.active) {
    const auto valid=ValidateShockGeometry(input.shock);
    if (!valid.ok()) return Failed(input,valid);
    const auto sphere=ShockGeneratingSphere(input.shock);
    const double gapM=std::fabs((particle.positionM-sphere.centerM).Norm()-sphere.radiusM);
    const double closingSpeedMPerS=speed+background.U.Norm()+
        std::fabs(input.shock.radialSpeedMPerS);
    if (gapM>0&&closingSpeedMPerS>0) {
      const double arrivalS=gapM/closingSpeedMPerS;
      // A particle can start on, or asymptotically approach, the exact surface.
      // Do not halve a purely geometric gap below the numerical time floor:
      // that would make an injected particle fail before its first move.
      // The analytic post-step intersection resolves this final approach;
      // all cell/diffusion/cooling/field/snapshot bounds remain in force.
      if (input.timeStepControls.shockCrossingFraction*arrivalS>=
          input.timeStepControls.minimumSubstepS)
        physics.timeToShockCrossingS=arrivalS;
    }
  }
  MoverResult result;
  result.particle = particle;
  result.selectedStep = Transport::SelectTimeStep(input.timeStepControls, physics);
  if (!result.selectedStep.status.ok())
    return Failed(input, result.selectedStep.status);

  Transport::RandomKey key;
  key.campaignSeed = input.campaignSeed;
  key.particleId = particle.stableId;
  key.step = particle.completedStep;
  key.substep = particle.substep;
  const double dtS = result.selectedStep.valueS;
  result.consumedTimeS = dtS;
  result.acceptedSubsteps = 1;
  result.finalLocal = input.local;
  if (input.model == RuntimeModel::TransportModel::Parker3D) {
    key.purpose = Transport::RandomPurpose::ParkerParallel;
    Transport::KeyedRandomStream random(key);
    key.purpose = Transport::RandomPurpose::PerpendicularFirst;
    Transport::KeyedRandomStream perpendicularFirst(key);
    key.purpose = Transport::RandomPurpose::PerpendicularSecond;
    Transport::KeyedRandomStream perpendicularSecond(key);
    Transport::ParkerParticleState state;
    state.positionM = particle.positionM;
    state.momentumKgMPerS = particle.momentumKgMPerS;
    Transport::ParkerLocalState local;
    local.bulkVelocityMPerS = background.U;
    local.bHat = background.bHat;
    local.curvaturePerM = background.curvature;
    local.divBhatPerM = background.divBhat;
    local.divUPerS = background.divU;
    local.kappaParallelM2PerS = input.local.kappaParallelM2PerS;
    local.dKappaParallelDsMPerS = input.local.dKappaParallelDsMPerS;
    local.kappaPerpendicularM2PerS = kappaPerpendicularM2PerS;
    local.dKappaPerpendicularDsMPerS =
        input.perpendicularDiffusion ==
                RuntimeModel::PerpendicularDiffusionMode::ConstantRatio
            ? input.kappaPerpendicularToParallelRatio *
                input.local.dKappaParallelDsMPerS
            : 0.0;
    local.driftVelocityMPerS = guidingCenterDrift;
    Transport::ParkerRandomStreams streams;
    streams.parallel = &random;
    streams.perpendicularFirst = kappaPerpendicularM2PerS > 0.0
        ? &perpendicularFirst : nullptr;
    streams.perpendicularSecond = kappaPerpendicularM2PerS > 0.0
        ? &perpendicularSecond : nullptr;
    const auto moved = Transport::AdvanceParker(state, local, dtS, streams);
    if (!moved.status.ok()) return Failed(input, moved.status);
    result.particle.positionM = moved.state.positionM;
    result.particle.momentumKgMPerS = moved.state.momentumKgMPerS;
  } else if (input.model ==
             RuntimeModel::TransportModel::FocusedDiffusion3D) {
    key.purpose = Transport::RandomPurpose::FocusedPitch;
    Transport::KeyedRandomStream random(key);
    key.purpose = Transport::RandomPurpose::PerpendicularFirst;
    Transport::KeyedRandomStream perpendicularFirst(key);
    key.purpose = Transport::RandomPurpose::PerpendicularSecond;
    Transport::KeyedRandomStream perpendicularSecond(key);
    Transport::FocusedParticleState state;
    state.positionM = particle.positionM;
    state.momentumKgMPerS = particle.momentumKgMPerS;
    state.mu = particle.mu;
    Transport::FocusedLocalState local;
    local.bulkVelocityMPerS = background.U;
    local.bHat = background.bHat;
    local.divBhatPerM = background.divBhat;
    local.divUPerS = background.divU;
    local.fieldAlignedStrainPerS = background.fieldAlignedStrain;
    local.dMuMuPerS = input.local.dMuMuPerS;
    local.dDmuMuDmuPerS = input.local.dDmuMuDmuPerS;
    local.kappaPerpendicularM2PerS = kappaPerpendicularM2PerS;
    local.driftVelocityMPerS = guidingCenterDrift;
    local.scheme = input.pitchScheme ==
            RuntimeModel::PitchAngleSchemeMode::ReflectingMilstein
        ? Transport::PitchAngleScheme::ReflectingMilstein
        : Transport::PitchAngleScheme::ReflectingEulerMaruyama;
    Transport::FocusedRandomStreams streams;
    streams.pitch = &random;
    streams.perpendicularFirst = kappaPerpendicularM2PerS > 0.0
        ? &perpendicularFirst : nullptr;
    streams.perpendicularSecond = kappaPerpendicularM2PerS > 0.0
        ? &perpendicularSecond : nullptr;
    const auto moved = Transport::AdvanceFocused(
        state, local, input.speciesMassKg, dtS, streams);
    if (!moved.status.ok()) return Failed(input, moved.status);
    result.particle.positionM = moved.state.positionM;
    result.particle.momentumKgMPerS = moved.state.momentumKgMPerS;
    result.particle.mu = moved.state.mu;
  } else {
    Transport::FocusedScatteringState state;
    state.particle.positionM = particle.positionM;
    state.particle.momentumKgMPerS = particle.momentumKgMPerS;
    state.particle.mu = particle.mu;
    state.remainingOpticalDepth =
        particle.remainingScatteringOpticalDepth;
    state.nextEventIndex = particle.nextScatteringEvent;
    Transport::FocusedScatteringLocalState local;
    local.deterministic.bulkVelocityMPerS = background.U;
    local.deterministic.bHat = background.bHat;
    local.deterministic.divBhatPerM = background.divBhat;
    local.deterministic.divUPerS = background.divU;
    local.deterministic.fieldAlignedStrainPerS =
        background.fieldAlignedStrain;
    local.deterministic.driftVelocityMPerS = guidingCenterDrift;
    local.meanFreePathM = input.local.meanFreePathM;
    local.alfvenSpeedMPerS = background.alfvenSpeedMpS;
    local.plusWaveFraction = input.local.plusWaveFraction;
    local.minusWaveFraction = input.local.minusWaveFraction;
    local.frame = input.focusedScatteringFrame ==
            RuntimeModel::FocusedScatteringFrame::PlasmaFrameIsotropic
        ? Transport::FocusedScatteringFrame::PlasmaFrameIsotropic
        : Transport::FocusedScatteringFrame::AlfvenWaveFrameIsotropic;
    const Transport::FocusedScatteringStepResult moved =
        Transport::AdvanceFocusedScattering(
            state, local, input.speciesMassKg, dtS, input.campaignSeed,
            particle.stableId, particle.completedStep,
            input.maximumScatteringEventsPerSubstep);
    if (!moved.status.ok()) return Failed(input, moved.status);
    result.particle.positionM = moved.state.particle.positionM;
    result.particle.momentumKgMPerS =
        moved.state.particle.momentumKgMPerS;
    result.particle.mu = moved.state.particle.mu;
    result.particle.remainingScatteringOpticalDepth =
        moved.state.remainingOpticalDepth;
    result.particle.nextScatteringEvent = moved.state.nextEventIndex;
  }

  ++result.particle.substep;
  const double radius = result.particle.positionM.Norm();
  if (radius < input.innerRadiusM) {
    result.disposition = ParticleDisposition::Absorbed;
    result.status = Core::Status(Core::StatusCode::InnerBoundary,
                                 "particle crossed the Parker/CME transport "
                                 "source shell");
    return result;
  }
  if (radius > input.outerRadiusM) {
    result.disposition = ParticleDisposition::Escaped;
    result.status = Core::Status(Core::StatusCode::DomainExit,
                                 "particle crossed the outer boundary");
    return result;
  }

  result.shockIntersection = FirstShockIntersection(
      particle.positionM, result.particle.positionM, dtS,
      input.shock, particle.lastShockGeneration);
  if (!result.shockIntersection.status.ok())
    return Failed(input, result.shockIntersection.status);
  if (result.shockIntersection.crossed)
    result.particle.lastShockGeneration = input.shock.generation;
  result.disposition = ParticleDisposition::Active;
  result.status = Core::Status::OK();
  return result;
}

MoverResult AdvanceParticleRequestedTime(const RequestedTimeAdvance& request) {
  const MoverInput& original = request.input;
  if (request.resolveLocal == nullptr || request.maximumSubsteps == 0) {
    return Failed(original, Invalid(
        "complete mover request requires a resolver and positive substep cap"));
  }
  if (!std::isfinite(original.requestedDtS) || original.requestedDtS <= 0.0) {
    return Failed(original, Invalid(
        "complete mover requested time must be finite and positive"));
  }

  // The tolerance is relative to the requested interval, but has an absolute
  // floor near one second.  It is used only to absorb subtraction roundoff;
  // it is never used to skip a physically meaningful minimum substep.
  const double toleranceS = 64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, original.requestedDtS);
  double consumedS = 0.0;
  ParticleRecord current = original.particle;
  MoverResult aggregate;
  aggregate.particle = current;
  aggregate.disposition = ParticleDisposition::Active;
  aggregate.status = Core::Status::OK();

  for (std::uint64_t accepted = 0; accepted < request.maximumSubsteps;
       ++accepted) {
    double remainingS = original.requestedDtS - consumedS;
    if (remainingS <= toleranceS) {
      // Make the public accounting exact after proving that only roundoff is
      // left.  This avoids restart drift caused by repeatedly carrying a
      // 1-ulp residue into the next global step.
      aggregate.consumedTimeS = original.requestedDtS;
      aggregate.particle = current;
      aggregate.disposition = ParticleDisposition::Active;
      aggregate.status = Core::Status::OK();
      return aggregate;
    }

    MoverInput substep = original;
    substep.particle = current;
    substep.requestedDtS = remainingS;
    if (substep.shock.active) {
      substep.shock.radiusAtStepStartM =
          original.shock.radiusAtStepStartM +
          original.shock.radialSpeedMPerS * consumedS;
    }
    Core::Status resolved = request.resolveLocal(
        current, consumedS, request.resolverContext, &substep.local);
    if (!resolved.ok()) {
      aggregate.status = resolved;
      aggregate.disposition = ParticleDisposition::Failed;
      aggregate.particle = current;
      aggregate.consumedTimeS = consumedS;
      return aggregate;
    }

    MoverResult moved = AdvanceParticle(substep);
    if (!std::isfinite(moved.consumedTimeS) || moved.consumedTimeS <= 0.0 ||
        moved.consumedTimeS > remainingS + toleranceS) {
      aggregate.status = Core::Status(
          Core::StatusCode::StepUnderflow,
          "accepted transport substep made no valid progress");
      aggregate.disposition = ParticleDisposition::Failed;
      aggregate.particle = current;
      aggregate.consumedTimeS = consumedS;
      return aggregate;
    }

    consumedS += moved.consumedTimeS;
    current = moved.particle;
    aggregate.selectedStep = moved.selectedStep;
    aggregate.finalLocal = moved.finalLocal;
    ++aggregate.acceptedSubsteps;
    if (moved.shockIntersection.crossed &&
        !aggregate.shockIntersection.crossed) {
      aggregate.shockIntersection = moved.shockIntersection;
    }

    if (!moved.status.ok() ||
        moved.disposition != ParticleDisposition::Active) {
      aggregate.status = moved.status;
      aggregate.disposition = moved.disposition;
      aggregate.particle = current;
      aggregate.consumedTimeS = consumedS;
      return aggregate;
    }
  }

  aggregate.status = Core::Status(
      Core::StatusCode::StepUnderflow,
      "transport substep cap reached before requested time was consumed");
  aggregate.disposition = ParticleDisposition::Failed;
  aggregate.particle = current;
  aggregate.consumedTimeS = consumedS;
  return aggregate;
}

}  // namespace Adapters
}  // namespace SEP3D

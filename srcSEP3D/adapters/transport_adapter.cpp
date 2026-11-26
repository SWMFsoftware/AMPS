#include "transport_adapter.h"
#include "shock_return.h"   // absorbing-shock test of one substep

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

// Rotating-frame inertial terms of the focused movers.  a_frame is the
// acceleration of the plasma frame due to the rotation (2 Omega x U' and,
// optionally, Omega x (Omega x x)); A = b.a_frame / v enters the pitch and
// momentum equations (focused_transport.h).  The same inertial force on the
// particle, F = -gamma m [2 Omega x (v mu b) + Omega x (Omega x x)], has the
// guiding-centre drift (F x B)/(q B^2), added to the spatial drift.
//
// Option application.sep.transport.frame_inertial_forces (CONFIGURATION.md;
// srcSEP3D/test/event_2013_04_11/README.md section 2.6).  It exists for uses
// of the focused movers in a rotating frame where these terms matter (for
// example electrons in the corona); the 2013-04-11 decks leave it off so the
// 3-D equations are those of the 1-D movers, which have no inertial terms.
//
// Derivation.  In a frame rotating with constant Omega, the plasma velocity
// is U' = U - Omega x x and a fluid element's inertial acceleration, written
// in rotating coordinates, is
//   a_inertial = dU'/dt + 2 Omega x U' + Omega x (Omega x x).
// The focused transport equation is written in the local plasma frame,
// which accelerates with a_inertial; its plasma-frame-acceleration terms
// therefore acquire the rotation part
//   a_frame = 2 Omega x U'  [+ Omega x (Omega x x)].
// Only this rotation part is added: the focused coefficients of
// transport/focused_transport.h carry no dU'/dt (b . DU/Dt) term.
// Only the component along b enters the gyrotropic equations, through
// A = (b . a_frame) / v [1/s]:
//   dmu/dt += -(1 - mu^2) A,   d ln p/dt += -mu A.
// When U' is exactly parallel to b (the corotating Parker/PFSS ambient),
// b . (Omega x U') = 0, so only the centrifugal part survives in A; the
// Coriolis force then acts only through the drift below.
//
// Drift.  The fictitious force on the particle's guiding centre in the
// rotating frame is F = -gamma m [2 Omega x v_gc + Omega x (Omega x x)];
// the code uses v_gc = v mu b, the plasma-frame parallel velocity.  Its
// guiding-centre drift is
//   V_F = (F x B) / (q B^2) = (F x b) / (q |B|)   [m/s],
// added to the existing V01 guiding-centre drift.  The drift is skipped for
// neutral species and where |B| = 0.
//
// Inputs: the background sample at the substep start (frame U' [m/s], b,
// |B| [T]), the frame position x [m], the plasma-frame speed v [m/s] and mu.
// The terms are evaluated once per substep and frozen over it, like every
// other local coefficient of the focused movers.  Returns zero terms when
// the option is off or v <= 0.
struct FrameInertialTerms {
  double accelerationOverSpeedPerS = 0.0;   // A = (b . a_frame) / v [1/s]
  Core::Vec3 driftVelocityMPerS;            // V_F [m/s, frame]
};

static FrameInertialTerms EvaluateFrameInertialTerms(
    const MoverInput& input, const Background::BackgroundSample& background,
    const Core::Vec3& positionM, double speedMPerS, double mu) {
  FrameInertialTerms terms;
  if (input.inertialForces == RuntimeModel::FrameInertialForces::None ||
      !(speedMPerS > 0.0)) return terms;
  const Core::Vec3 omega = input.frameAngularVelocityRadPerS;
  const bool centrifugal = input.inertialForces ==
      RuntimeModel::FrameInertialForces::CoriolisAndCentrifugal;
  // Omega x (Omega x x) = -Omega^2 x_perp (x_perp: x normal to the axis),
  // the centripetal acceleration of a point at rest in the rotating frame.
  const Core::Vec3 centrifugalAcceleration =
      centrifugal ? omega.Cross(omega.Cross(positionM)) : Core::Vec3();
  // a_frame of the plasma, from the frame bulk velocity U'.
  const Core::Vec3 plasma = omega.Cross(background.U) * 2.0 + centrifugalAcceleration;
  terms.accelerationOverSpeedPerS = background.bHat.Dot(plasma) / speedMPerS;
  if (input.speciesChargeC != 0.0 && background.absB > 0.0) {
    // Lorentz factor of the plasma-frame speed (the floor only protects the
    // square root against v rounding to c).
    const double gamma = 1.0 / std::sqrt(std::max(1e-300,
        1.0 - speedMPerS * speedMPerS / (Core::Const::c * Core::Const::c)));
    // Acceleration the rotation imposes on the guiding centre moving with
    // v mu b: Coriolis plus the optional centrifugal part.
    const Core::Vec3 particle =
        omega.Cross(background.bHat * (speedMPerS * mu)) * 2.0 + centrifugalAcceleration;
    // F = -gamma m a  (relativistic mass gamma m), and V_F = F x b / (q |B|).
    const Core::Vec3 force = particle * (-gamma * input.speciesMassKg);
    terms.driftVelocityMPerS = force.Cross(background.bHat) *
        (1.0 / (input.speciesChargeC * background.absB));
  }
  return terms;
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
  // Only the focused-diffusion mover realizes D_mumu continuously; the
  // library's scattering limit is meaningless for the other movers.
  physics.pitchAngleScatteringLimitS =
      input.model == RuntimeModel::TransportModel::FocusedDiffusion3D
          ? input.local.pitchAngleScatteringLimitS
          : std::numeric_limits<double>::infinity();
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
    {
      // Optional rotating-frame terms (zero unless frame_inertial_forces is
      // set): A for the pitch and momentum equations, V_F added to the drift.
      const FrameInertialTerms inertial = EvaluateFrameInertialTerms(
          input, background, particle.positionM, speed, particle.mu);
      local.frameAccelerationOverSpeedPerS = inertial.accelerationOverSpeedPerS;
      local.driftVelocityMPerS += inertial.driftVelocityMPerS;
    }
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
    {
      // Optional rotating-frame terms, as for focused diffusion: they enter
      // the deterministic (focusing/flow/cooling) part between scatterings.
      const FrameInertialTerms inertial = EvaluateFrameInertialTerms(
          input, background, particle.positionM, speed, particle.mu);
      local.deterministic.frameAccelerationOverSpeedPerS =
          inertial.accelerationOverSpeedPerS;
      local.deterministic.driftVelocityMPerS += inertial.driftVelocityMPerS;
    }
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
  // Absorbing shock (source.shock.returning_particles = absorb; README of
  // test/event_2013_04_11 section 2.5).  The tube-line source injects the
  // 1-D NET escaping spectrum, which already counts only escaping particles
  // (1-D README section 2.2), so a particle that returns into the downstream
  // side of an accepted shock leaves the population, as in 1-D.  The test runs after the
  // substep's move from x0 = particle.positionM to x1 =
  // result.particle.positionM, against the front moving linearly over the
  // substep (input.shock holds its geometry at this substep's start).
  if (input.absorbAtShock && input.shock.active) {
    // The diffusion normal to the front uses the parallel coefficient and the
    // field direction at the substep start (Parker only).
    //   kappa_n dt = kappa_par (b . e_n)^2 dt   [m^2],
    // e_n the front normal at x0 (radial for the heliocentric sphere).  A
    // displacement ds along b moves the particle ds (b . e_n) normal to the
    // front, so on a locally straight line the 1-D bridge along the line,
    // exp(-d0 d1 / (kappa_par dt)) with arc distances d = g / (b . e_n),
    // equals exp(-g0 g1 / (kappa_n dt)) with the normal distances g of
    // shock_return.cpp.  Perpendicular diffusion, when configured, is not
    // included in kappa_n.
    // The focused movers keep kappa_n = 0: their substep path is
    // represented by the straight chord x0 -> x1, tested exactly for a dip
    // below the moving front.
    double kappaNormalDtM2 = 0.0;
    // Separate keyed stream for the bridge draw: the mover's own draws do not
    // depend on whether the test runs.
    Transport::RandomKey bridgeKey = key;
    bridgeKey.purpose = Transport::RandomPurpose::ShockReturnBridge;
    Transport::KeyedRandomStream bridge(bridgeKey);
    if (input.model == RuntimeModel::TransportModel::Parker3D) {
      const GeneratingSphere sphere = ShockGeneratingSphere(input.shock);
      const Core::Vec3 normal = (particle.positionM - sphere.centerM).Normalized();
      const double cosine = background.bHat.Dot(normal);
      kappaNormalDtM2 = input.local.kappaParallelM2PerS * cosine * cosine * dtS;
    }
    const ShockReturnTest test = TestShockReturn(
        particle.positionM, result.particle.positionM, dtS, input.shock,
        kappaNormalDtM2, &bridge);
    if (test.hit) {
      // test.fraction tau in [0,1] locates the crossing on the chord and in
      // time (t = moverEpochS + tau dt).  Rebuild the front at that time and
      // project the crossing onto it (from the generating-sphere centre) to
      // get the front point whose acceptance decides the removal.
      ExpandingShock atHit = input.shock;
      atHit.radiusAtStepStartM += input.shock.radialSpeedMPerS * dtS * test.fraction;
      const GeneratingSphere sphere = ShockGeneratingSphere(atHit);
      const Core::Vec3 crossing = particle.positionM +
          (result.particle.positionM - particle.positionM) * test.fraction;
      const Core::Vec3 frontPoint = sphere.centerM +
          (crossing - sphere.centerM).Normalized() * sphere.radiusM;
      // Geometric gate: the outward cap of a finite SSE front (always true
      // for the heliocentric sphere).  Physical gate: the provider's front
      // state at the hit epoch must be a supported, accepted shock.  A
      // crossing through a non-shock part of the front (sub-critical or
      // unsupported) leaves the particle active.
      bool accepted = OnOutwardShockCap(atHit, frontPoint, sphere);
      if (accepted && input.shockAccepted != nullptr) {
        const Core::Status queried = input.shockAccepted(
            frontPoint, input.moverEpochS + test.fraction * dtS, &accepted);
        if (!queried.ok()) return Failed(input, queried);
      }
      // Removal: the particle is placed at the crossing and reported as
      // Absorbed with status InnerBoundary; the AMPS adapter deletes it and
      // the ledger records the absorbed disposition.
      if (accepted) {
        result.particle.positionM = crossing;
        result.disposition = ParticleDisposition::Absorbed;
        result.status = Core::Status(Core::StatusCode::InnerBoundary,
                                     "particle returned to the accepted shock");
        return result;
      }
    }
  }
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
    // Model time at this substep's start: the absorbing shock queries the
    // front's acceptance at moverEpochS + tau dt.
    substep.moverEpochS = original.moverEpochS + consumedS;
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

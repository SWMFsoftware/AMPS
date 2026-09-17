// ============================================================================
// Phase-A host-neutral particle mover adapter.
//
// This layer is deliberately free of pic.h.  It validates the complete record
// produced by an AMPS (or test) buffer translator, selects exactly one Phase-P
// core, classifies physical boundaries, and records a possible expanding-shock
// intersection.  It contains no duplicate Parker or focused equations.
// ============================================================================

#ifndef SEP3D_ADAPTERS_TRANSPORT_ADAPTER_H
#define SEP3D_ADAPTERS_TRANSPORT_ADAPTER_H

#include "../background/bg_provider.h"
#include "../runtime/run_configuration.h"
#include "../transport/focused_transport.h"
#include "../transport/parker_transport.h"
#include "../transport/time_step.h"

#include <cstdint>
#include <string>
#include <vector>

namespace SEP3D {
namespace Adapters {

struct ParticleRecord {
  std::uint64_t stableId = 0;
  int species = -1;
  Core::Vec3 positionM;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
  double gyrophaseRad = 0.0;
  double statisticalWeight = 0.0;
  std::uint64_t completedStep = 0;
  std::uint64_t substep = 0;
  std::uint64_t lastShockGeneration = 0;
};

struct LocalTransportRecord {
  Background::BackgroundSample background;
  double cellSizeM = 0.0;
  double kappaParallelM2PerS = 0.0;
  double dKappaParallelDsMPerS = 0.0;
  double dMuMuPerS = 0.0;
  double dDmuMuDmuPerS = 0.0;
  double fractionalFieldVariationPerS = 0.0;
  double timeToSnapshotBoundaryS = 0.0;
};

struct ExpandingSphericalShock {
  Core::Vec3 centerM;
  double radiusAtStepStartM = 0.0;
  double radialSpeedMPerS = 0.0;
  std::uint64_t generation = 0;
  bool active = false;
};

struct ShockIntersection {
  Core::Status status;
  bool crossed = false;
  double stepFraction = 0.0;
  Core::Vec3 positionM;
  std::uint64_t generation = 0;
};

enum class ParticleDisposition { Active, Escaped, Absorbed, Failed };

const char* Name(ParticleDisposition disposition);

struct MoverInput {
  RuntimeModel::TransportModel model = RuntimeModel::TransportModel::Parker3D;
  ParticleRecord particle;
  LocalTransportRecord local;
  ExpandingSphericalShock shock;
  double speciesMassKg = 0.0;
  double requestedDtS = 0.0;
  double innerRadiusM = 0.0;
  double outerRadiusM = 0.0;
  std::uint64_t campaignSeed = 0;
  Transport::TimeStepControls timeStepControls;
  RuntimeModel::PitchAngleSchemeMode pitchScheme =
      RuntimeModel::PitchAngleSchemeMode::ReflectingMilstein;
};

struct MoverResult {
  Core::Status status;
  ParticleRecord particle;
  ParticleDisposition disposition = ParticleDisposition::Failed;
  Transport::TimeStepSelection selectedStep;
  ShockIntersection shockIntersection;
};

// Solve the first intersection between a straight particle substep and a
// sphere whose radius changes linearly over the same interval. The result is
// the smallest root in [0,1]. A generation already recorded by the particle
// is suppressed so a trajectory cannot inject twice on one shock surface.
ShockIntersection FirstShockIntersection(
    const Core::Vec3& initialPositionM,
    const Core::Vec3& finalPositionM,
    double dtS,
    const ExpandingSphericalShock& shock,
    std::uint64_t lastShockGeneration);

class ProductionMoverRegistry final {
 public:
  static const std::vector<std::string>& CanonicalNames();
  static Core::Status Validate(RuntimeModel::TransportModel selected);
};

MoverResult AdvanceParticle(const MoverInput& input);

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_TRANSPORT_ADAPTER_H

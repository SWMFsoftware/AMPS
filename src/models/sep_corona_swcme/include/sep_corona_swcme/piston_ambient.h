#ifndef SEP_CORONA_SWCME_PISTON_AMBIENT_H
#define SEP_CORONA_SWCME_PISTON_AMBIENT_H

#include "sep_corona_swcme/ambient_state.h"
#include "sep_corona_swcme/piston_solver.h"

#include <array>
#include <memory>
#include <vector>

namespace SEP { namespace CoronaSwcme {

// One frozen Eulerian projection of the maintained 3-D ambient authority onto
// a fixed HCI radial tube.  Freezing the projection epoch is intentional: the
// Level-B equations prescribe time-independent f_amb, h_amb and S_b along a
// ray.  Solar-rotation/time dependence remains available through the ambient
// authority outside the disturbed tube; changing this contract would require
// an additional explicit source and energy term.
struct PistonAmbientRayState {
  double radiusM = 0.0;
  double densityKgM3 = 0.0;
  double pressurePa = 0.0;
  double specificInternalEnergyJPerKg = 0.0;
  double radialVelocityMPerS = 0.0;
  double radialVelocityGradientPerS = 0.0;
  double materialRadialAccelerationMPerS2 = 0.0;
  std::array<double,2> transverseVelocityMPerS{};
  double radialMagneticFieldT = 0.0;
  std::array<double,2> transverseMagneticFieldT{};
  double reducedFastSpeedMPerS = 0.0;
  double canonicalFastSpeedMPerS = 0.0;
  PistonVolumeSource source;
  // Kept separately for the validity ledger even though the solver consumes
  // their well-conditioned sum source.radialAccelerationMPerS2.
  double gravityAccelerationMPerS2 = 0.0;
  double ambientMaintainingAccelerationMPerS2 = 0.0;
  AmbientRegion region = AmbientRegion::PfssClosed;
  int magneticSector = 0;
};

class PistonAmbientProjection final {
 public:
  static Core::Result<std::shared_ptr<const PistonAmbientProjection>> Create(
      std::shared_ptr<const AmbientModel> ambient,
      CoronalCME::Vec3 rayDirection,double frozenEpochS,
      double derivativeRelativeStep=1e-5);

  Core::Result<PistonAmbientRayState> Evaluate(double radiusM) const;
  // Primitive-only sampling avoids radial derivative/topology work when
  // initializing many material cells.  Source fields are zero in this
  // result; callers needing f_amb/h_amb/S_b must use Evaluate().
  Core::Result<PistonAmbientRayState> EvaluatePrimitive(double radiusM) const;

  const CoronalCME::Vec3& RayDirection() const noexcept { return ray_; }
  const CoronalCME::Vec3& FirstTangent() const noexcept { return tangent1_; }
  const CoronalCME::Vec3& SecondTangent() const noexcept { return tangent2_; }
  double FrozenEpochS() const noexcept { return epochS_; }

 private:
  std::shared_ptr<const AmbientModel> ambient_;
  CoronalCME::Vec3 ray_,tangent1_,tangent2_;
  double epochS_ = 0.0;
  double relativeStep_ = 0.0;
};

// Frozen source samples are tabulated once per ray.  Calling the 3-D ambient
// topology tracer and its branch-safe derivative stencil at every RK stage
// would make a production propagation prohibitively expensive and could let
// MPI scheduling change which adaptive derivative attempt was used.  Linear
// interpolation is deliberately simple; table-spacing convergence is an
// acceptance dimension and no extrapolation is permitted.
class PistonAmbientSourceTable final {
 public:
  static Core::Result<std::shared_ptr<const PistonAmbientSourceTable>> Create(
      std::shared_ptr<const PistonAmbientProjection> projection,
      double minimumRadiusM,double maximumRadiusM,int points,
      bool logarithmicSpacing=false);

  Core::Result<PistonVolumeSource> Evaluate(double radiusM) const;
  double MaximumAmbientForceOverGravity() const noexcept {
    return maximumAmbientForceOverGravity_;
  }
  double MaximumAbsoluteHeatingWPerM3() const noexcept {
    return maximumAbsoluteHeatingWPerM3_;
  }

 private:
  std::shared_ptr<const PistonAmbientProjection> projection_;
  std::vector<double> radius_;
  std::vector<PistonVolumeSource> source_;
  std::vector<AmbientRegion> region_;
  std::vector<int> sector_;
  double maximumAmbientForceOverGravity_ = 0.0;
  double maximumAbsoluteHeatingWPerM3_ = 0.0;
};

// Material path dr/dt=u_a(r) through the frozen projected ambient.  It is used
// for a comoving piston/outer boundary in well-balance tests and later for
// ambient-buffer insertion.  Each interval is evaluated with quintic Hermite
// data (r,u,u du/dr), so position, velocity and acceleration are continuous at
// table nodes; no large-position finite difference defines boundary speed.
class PistonAmbientTrajectory final {
 public:
  static Core::Result<std::shared_ptr<const PistonAmbientTrajectory>> Create(
      std::shared_ptr<const PistonAmbientProjection> projection,
      double initialRadiusM,double startS,double endS,double maximumStepS);

  Core::Result<CoronalCME::KinematicValue> Evaluate(double timeS) const;

 private:
  std::vector<double> time_,radius_,velocity_,acceleration_;
};

} } // namespace SEP::CoronaSwcme

#endif

#ifndef SEP_CORONA_SWCME_SPHERICAL_SHOCK_ONLY_H
#define SEP_CORONA_SWCME_SPHERICAL_SHOCK_ONLY_H

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "sep_corona_swcme/ambient_state.h"
#include "sep_corona_swcme/open_tube_isopoly.h"
#include "sep_coronal_cme/mhd_jump_solver.h"
#include "swcme_kinematics.hpp"

namespace SEP {
namespace CoronaSwcme {
namespace SphericalShockOnly {

// This is an explicitly prescribed ISOTHERMAL KINEMATIC candidate, not an
// observationally calibrated wind or a qualified downstream CME volume.
// Physical particle spheres and provider/derivative/material guards differ:
// guards are actual caller-supplied AmbientDefinition support, never 8/10-AU
// defaults or an extension inferred from a requested query.
struct Definition {
  AmbientDefinition ambient;
  swcme::kinematics::Config front;
  double physicalInnerRadiusM = 0, physicalOuterRadiusM = 0;
  std::string inputIdentity;
};
struct Numerics {
  double positionToleranceM = 0, deformationTolerance = 0;
  double maximumTimeStepS = 0, maximumArcStepM = 0;
  std::size_t maximumAttempts = 0;
};
enum class PhysicalSupport { InnerGuard, ParticleDomain, OuterGuard };
enum class FrontSide { Upstream, OnFront, BehindFrontReferenceOnly };
struct Sample {
  AmbientState ambient;
  PhysicalSupport support = PhysicalSupport::ParticleDomain;
  FrontSide side = FrontSide::BehindFrontReferenceOnly;
  double continuityResidualKgM3S = 0;
  CoronalCME::Vec3 inductionResidualTPerS;
  // f_ext=rho DU/Dt+grad(p)-rho*g-(curl B)xB/mu0 [N/m^3]. Its
  // field-aligned value diagnoses the radial-Parker/nonradial-tube momentum
  // mismatch; no unapproved coronal heating/force closure is inferred.
  CoronalCME::Vec3 requiredExternalForceNPerM3;
  double fieldAlignedForcePerMassMPerS2 = 0, requiredMechanicalPowerWPerM3 = 0;
  double requiredIsothermalHeatWPerM3 = 0;
};
struct Motion {
  CoronalCME::Vec3 positionM, velocityMPerS;
  std::array<double, 9> deformation{};
  double jacobian = 0;
};
struct Shock {
  double epochS = 0, radiusM = 0, normalSpeedMPerS = 0, fastMach = 0;
  CoronalCME::Vec3 positionM, outwardNormal;
  Sample ambient;
  bool physicalFast = false;
  // A sub-fast geometric front has no fabricated jump. The canonical RH
  // residual/entropy/branch gates are unchanged; solver failures propagate.
  std::optional<CoronalCME::MhdShockSolution> jump;
};
// A resolved tube provides the actual local upstream primitive and thermal
// budget at the spherical front. It has no transverse gradient tensor or
// downstream volume. A sub-fast front retains its geometry with no jump.
struct TubeShock {
  double epochS=0, radiusM=0, normalSpeedMPerS=0, fastMach=0;
  CoronalCME::Vec3 outwardNormal;
  OpenTubeIsoPolySample upstream;
  PhysicalSupport support=PhysicalSupport::ParticleDomain;
  std::optional<CoronalCME::MhdShockSolution> jump;
};
struct Connection {
  std::size_t segment = 0;
  double fraction = 0, arcLengthM = 0, surfaceResidualM = 0;
  bool tangent = false;
  Shock shock;
};
class Event final {
 public:
  static Core::Result<std::shared_ptr<const Event>> Create(const Definition&);
  const Definition& Configuration() const {
    return definition_;
  }
  const std::string& Identity() const {
    return identity_;
  }
  const std::shared_ptr<const AmbientModel>& Ambient() const {
    return ambient_;
  }
  Core::Result<Sample> Evaluate(CoronalCME::Vec3, double epochS) const;
  Core::Result<Motion> Move(CoronalCME::Vec3, double beginS, double endS,
                            const Numerics&) const;
  // Exactly a complete caller-owned AMPS step. Numerical RK stages are
  // quadrature only; this API never shortens dt or touches sampling clocks.
  Core::Result<Motion> AdvanceFixed(CoronalCME::Vec3, double beginS, double endS,
                                    double fixedDtS, const Numerics&) const;
  Core::Result<Shock> Surface(CoronalCME::Vec3 ray, double epochS) const;
  Core::Result<TubeShock> SurfaceOnTube(const OpenTubeIsoPolyWind& wind,
                                      double inertialAzimuthRad,double epochS) const;
  Core::Result<std::vector<Connection>> Connect(const std::vector<CoronalCME::Vec3>&,
                                                double epochS,
                                                double spatialToleranceM) const;

 private:
  Event() = default;
  Core::Result<Motion> Integrate(CoronalCME::Vec3, double beginS, double durationS,
                                 const Numerics&) const;
  Definition definition_;
  std::shared_ptr<const AmbientModel> ambient_;
  std::string identity_;
};

}  // namespace SphericalShockOnly
}  // namespace CoronaSwcme
}  // namespace SEP
#endif

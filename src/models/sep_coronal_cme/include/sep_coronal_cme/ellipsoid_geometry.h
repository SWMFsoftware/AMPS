#ifndef SEP_CORONAL_CME_ELLIPSOID_GEOMETRY_H
#define SEP_CORONAL_CME_ELLIPSOID_GEOMETRY_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <cstdint>
#include <utility>
#include <vector>

namespace SEP { namespace CoronalCME {

struct KinematicValue {
  double value = 0.0;
  double firstDerivative = 0.0;
  double secondDerivative = 0.0;
};

// Independent quintic rate transition from Section 8.4.  qReference is the
// value at referenceTime; the pre-transition linear history is continued
// exactly to transitionStart before the C2 rate blend begins.
Core::Result<KinematicValue> EvaluateSmoothKinematics(
    double timeS, double referenceTimeS, double qReference,
    double transitionStartS, double transitionDurationS,
    double initialRatePerS, double finalRatePerS);

struct HermiteKnot {
  double timeS = 0.0;
  double value = 0.0;
  double ratePerS = 0.0;
};
Core::Result<KinematicValue> EvaluateCubicHermiteHistory(
    const std::vector<HermiteKnot>& knots, double timeS);

struct RadialPrincipalBasis {
  Vec3 radial;
  Vec3 firstLateral;
  Vec3 secondLateral;
};
Core::Result<RadialPrincipalBasis> BuildRadialPrincipalBasis(
    double latitudeRad, double longitudeRad, double lateralTiltRad);

struct EllipsoidKinematics {
  KinematicValue centerDistanceM;
  KinematicValue radialSemiAxisM;
  KinematicValue firstLateralSemiAxisM;
  KinematicValue secondLateralSemiAxisM;
};

struct SurfaceEvaluation {
  double implicitValue = 0.0;
  Vec3 outwardNormal;
  double normalSpeedMPerS = 0.0;
  double meanCurvaturePerM = 0.0;
  double gaussianCurvaturePerM2 = 0.0;
};

struct SurfacePatch {
  std::uint64_t physicalId = 0;
  Vec3 centerM;
  Vec3 outwardNormal;
  double areaM2 = 0.0;
  double polarParameterRad = 0.0;
  double azimuthParameterRad = 0.0;
};

class FixedOrientationEllipsoid {
 public:
  static Core::Result<FixedOrientationEllipsoid> FromCenter(
      const RadialPrincipalBasis& basis,
      const EllipsoidKinematics& kinematics, double solarRadiusM = 0.0);
  static Core::Result<FixedOrientationEllipsoid> FromApex(
      const RadialPrincipalBasis& basis, KinematicValue apexRadiusM,
      KinematicValue radialSemiAxisM,
      KinematicValue firstLateralSemiAxisM,
      KinematicValue secondLateralSemiAxisM,
      double solarRadiusM = 0.0);

  Vec3 Point(double polarParameterRad,
             double azimuthParameterRad) const;
  Vec3 SurfaceVelocity(double polarParameterRad,
                       double azimuthParameterRad) const;
  Core::Result<SurfaceEvaluation> Evaluate(Vec3 positionM) const;

  // Midpoint surface quadrature keeps physical IDs tied to parameter cells,
  // not MPI rank or load-balancing order.  A clipped-out cell leaves an ID
  // gap rather than renumbering every later physical patch.
  Core::Result<std::vector<SurfacePatch>> Tessellate(
      int polarCells, int azimuthCells) const;

  double CenterDistanceM() const noexcept {
    return kinematics_.centerDistanceM.value;
  }
  double ApexRadiusM() const noexcept {
    return kinematics_.centerDistanceM.value +
        kinematics_.radialSemiAxisM.value;
  }
  const EllipsoidKinematics& Kinematics() const noexcept {
    return kinematics_;
  }
  const RadialPrincipalBasis& Basis() const noexcept { return basis_; }

 private:
  RadialPrincipalBasis basis_;
  EllipsoidKinematics kinematics_;
  double solarRadiusM_ = 0.0;
};

// A dependency-light deterministic index.  It stores immutable physical
// patches and returns stable IDs within a query ball; the AMPS adapter may
// replace its linear scan with a mesh-aware acceleration structure without
// changing surface ownership or IDs.
class PatchSpatialIndex {
 public:
  explicit PatchSpatialIndex(std::vector<SurfacePatch> patches)
      : patches_(std::move(patches)) {}
  std::vector<std::uint64_t> Query(Vec3 pointM, double radiusM) const;
 private:
  std::vector<SurfacePatch> patches_;
};

struct PistonNestingDiagnostics {
  double minimumRadialSeparationM = 0.0;
  std::uint64_t limitingPatchId = 0;
};
Core::Result<PistonNestingDiagnostics> CheckPistonNesting(
    const FixedOrientationEllipsoid& front,
    const FixedOrientationEllipsoid& piston,
    double requiredMinimumSeparationM, int polarSamples,
    int azimuthSamples);

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_ELLIPSOID_GEOMETRY_H

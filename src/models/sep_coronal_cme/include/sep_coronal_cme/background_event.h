#ifndef SEP_CORONAL_CME_BACKGROUND_EVENT_H
#define SEP_CORONAL_CME_BACKGROUND_EVENT_H

#include <array>
#include <functional>
#include <memory>

#include "sep_coronal_cme/configuration_parser.h"
#include "sep_coronal_cme/ellipsoid_geometry.h"
#include "sep_coronal_cme/pfss_harmonics.h"
#include "sep_coronal_cme/plasma_eos.h"
#include "sep_coronal_cme/shock_provider.h"

namespace SEP {
namespace CoronalCME {

// This is background-only assembly, independent of source plans, particle
// records, random streams and wave feedback. "Materialize" means resolving
// immutable typed input, not claiming a downstream fluid/material-map closure.
enum class HistoryKind { SmoothRateTransition, CubicHermiteC1 };
struct ScalarGeometryHistory {
  HistoryKind kind = HistoryKind::SmoothRateTransition;
  double referenceTimeS = 0, referenceValueM = 0;
  double initialRateMPerS = 0, finalRateMPerS = 0;
  double transitionStartS = 0, transitionDurationS = 0;
  std::vector<HermiteKnot> knots;
  Core::Result<KinematicValue> At(double timeS) const;
  // Exact scalar extrema (analytic derivative roots), not checks at knots.
  Core::Result<std::pair<double, double>> Range(double beginS, double endS) const;
};
struct GeometryHistory {
  RadialPrincipalBasis basis;
  // center distance, radial semiaxis, first/second lateral semiaxis, in SI.
  std::array<ScalarGeometryHistory, 4> components;
  Core::Result<EllipsoidKinematics> At(double timeS) const;
};
enum class AttachmentPolicy { RequireSolarIntersection, AllowDetachedFront };
enum class PlasmaCoveragePolicy {
  RequireEveryExposedPatch,
  GeometryOnlyOutsideCoverage
};
struct EventSupport {
  double beginS = 0, endS = 0, solarRadiusM = 0, qualifiedSourceInnerRadiusM = 0;
  double plasmaMinimumRadiusM = 0, plasmaMaximumRadiusM = 0;
  std::string frame, epochUtc;
  AttachmentPolicy attachment = AttachmentPolicy::RequireSolarIntersection;
  PlasmaCoveragePolicy plasmaCoverage = PlasmaCoveragePolicy::RequireEveryExposedPatch;
};
struct EventNumerics {
  // Explicit certificate work bounds. An unresolved interval rejects; it is
  // never accepted because a finite time sample happened to miss an extremum.
  unsigned maximumCertificateDepth = 16;
  double minimumCertificateIntervalS = 1e-8;
  int polarCells = 0, azimuthCells = 0;
};
struct BackgroundEventDefinition {
  GeometryHistory front, piston;
  bool pistonEnabled = false;
  double minimumPistonSeparationM = 0, gammaAdiabatic = 0;
  EventSupport support;
  EventNumerics numerics;
  // Input provenance is immutable authority metadata. Create derives the
  // final physicsIdentity from ALL typed physical fields and this provenance;
  // a caller cannot carry an old ID across changed geometry/EOS/support.
  std::string inputIdentity, physicsIdentity, ambientIdentity;
  bool includeElectronMass = false;
};

struct GeometryAsset {
  std::string bytes, expectedSha256;
};
// The application/asset layer owns acquisition and the expected content hash.
// No library I/O, path guessing or implicit research-format conversion occurs.
using GeometryAssetReader =
    std::function<Core::Result<GeometryAsset>(const std::string&)>;
Core::Result<BackgroundEventDefinition> MaterializeBackgroundEvent(
    const ModelConfiguration&, const EventSupport&, const EventNumerics&,
    const std::string& ambientIdentity, const GeometryAssetReader& = {});

// A retained immutable ambient epoch is supplied by its actual authority.
// Existing PFSS/open-wind/closed-plasma/Parker implementations can implement
// this contract without making the coronal archive depend on a composite or
// application library. Every query is bound to the same frame, EOS and epoch.
struct AmbientEpochIdentity {
  std::string physicsIdentity, frame, epochUtc;
  double timeS = 0, validFromS = 0, validUntilS = 0, gammaAdiabatic = 0;
  std::uint64_t generation = 0;
  bool includeElectronMass = false;
};
enum class AmbientRegion { PfssOpen, PfssClosed, FiniteShellScs, ParkerExterior };
enum class QuerySide { Ordinary, Minus, Plus };
enum class EventDifferenceStencil {
  CentralSecondOrder,
  ForwardSecondOrder,
  BackwardSecondOrder
};
struct EventAmbientSample {
  MhdPrimitiveState primitive;
  double electronDensityM3 = 0, electronTemperatureK = 0;
  std::vector<IonSpecies> ions;
  AmbientRegion region = AmbientRegion::PfssClosed;
  FieldLineTopology topology = FieldLineTopology::Invalid;
  int sector = 0;
  bool oneSided = false;
  QuerySide side = QuerySide::Ordinary;
  std::array<double, 9> gradientBPerM{}, gradientUPerS{};
  // Numerical Cartesian stencil directions are not physical radial traces.
  std::array<EventDifferenceStencil, 3> derivativeStencils{};
  bool derivativesValid = false;
};
class AmbientEpoch {
 public:
  virtual ~AmbientEpoch() = default;
  virtual const AmbientEpochIdentity& Identity() const = 0;
  virtual Core::Result<EventAmbientSample> At(Vec3 positionM, QuerySide) const = 0;
};

enum class PatchSupport { PlasmaAvailable, BelowPlasmaCoverage, AbovePlasmaCoverage };
// The center sample is a quadrature point, not a whole surface cell. A
// conservative parameter-cell enclosure separately records whether complete
// cell coverage is proved or an interface crossing remains unresolved.
enum class ParameterCellCoverage {
  InsidePlasmaSupport,
  OutsidePlasmaSupport,
  UnresolvedBoundaryCrossing
};
enum class SourceRadialSupport { QualifiedRadius, BelowQualifiedSourceRadius };
struct BackgroundEventPatch {
  SurfacePatch geometry;
  double normalSpeedMPerS = 0;
  PatchSupport plasmaSupport = PatchSupport::PlasmaAvailable;
  ParameterCellCoverage parameterCellCoverage =
      ParameterCellCoverage::UnresolvedBoundaryCrossing;
  double parameterCellMinimumRadiusM = 0, parameterCellMaximumRadiusM = 0;
  SourceRadialSupport sourceRadialSupport =
      SourceRadialSupport::BelowQualifiedSourceRadius;
  EventAmbientSample ambient;
};
struct BackgroundEventEpoch {
  std::uint64_t generation = 0;
  double timeS = 0;
  std::string eventIdentity;
  FixedOrientationEllipsoid front, piston;
  EllipsoidRadialExtent frontRadialExtent;
  bool pistonEnabled = false;
  std::vector<BackgroundEventPatch> frontPatches;
  std::vector<SurfacePatch> pistonPatches;
  std::shared_ptr<const AmbientEpoch> ambient;
  std::shared_ptr<const ShockSurfaceSnapshot> shocks;
};
class BackgroundEvent final {
 public:
  static Core::Result<std::shared_ptr<BackgroundEvent>> Create(
      BackgroundEventDefinition);
  // All geometry, ambient queries, EOS checks and local RH states are prepared
  // privately. Failure leaves the prior complete epoch AND its local-shock
  // generation unchanged. Source eligibility/rates are explicitly disabled.
  Core::Result<std::shared_ptr<const BackgroundEventEpoch>> Prepare(
      double timeS, std::shared_ptr<const AmbientEpoch>,
      QuerySide = QuerySide::Ordinary);
  std::shared_ptr<const BackgroundEventEpoch> Current() const {
    return current_;
  }
  const BackgroundEventDefinition& Definition() const {
    return definition_;
  }

 private:
  BackgroundEventDefinition definition_;
  TransactionalShockProvider shockProvider_;
  std::shared_ptr<const BackgroundEventEpoch> current_;
};

}  // namespace CoronalCME
}  // namespace SEP
#endif

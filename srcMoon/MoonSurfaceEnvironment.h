#ifndef SRCMOON_MOON_SURFACE_ENVIRONMENT_H
#define SRCMOON_MOON_SURFACE_ENVIRONMENT_H

#include "LunarSurface.h"
#include "species/IonCollector.h"

#include <array>
#include <cstdint>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace Moon {
namespace SurfaceEnvironment {

struct RayHit {
  bool hit = false;
  double distanceM = std::numeric_limits<double>::infinity();
  std::size_t face = std::numeric_limits<std::size_t>::max();
  std::array<double, 3> positionM{{0.0, 0.0, 0.0}};
  std::array<double, 3> outwardNormal{{0.0, 0.0, 0.0}};
};

// Reusable body-fixed terrain query service.  Its BVH is derived only from the
// production triangulation; D02 illumination products never enter this class,
// preserving them as independent I15/I20 evidence.
class TerrainRayTracer {
 public:
  bool Initialize(const Surface::Triangulation& triangulation,
                  std::string* error);
  bool IntersectRay(const std::array<double, 3>& originM,
                    const std::array<double, 3>& direction,
                    double minimumDistanceM, double maximumDistanceM,
                    std::size_t excludedFace, RayHit* hit,
                    std::string* error) const;
  bool IsSunVisible(const std::array<double, 3>& surfacePositionM,
                    const std::array<double, 3>& outwardNormal,
                    const std::array<double, 3>& unitSunDirection,
                    double rayOriginOffsetM, std::size_t sourceFace,
                    bool* visible, std::string* error) const;
  const Surface::Triangulation& triangulation() const { return mesh_; }

 private:
  struct Node {
    std::array<double, 3> minimum{{0.0, 0.0, 0.0}};
    std::array<double, 3> maximum{{0.0, 0.0, 0.0}};
    std::size_t begin = 0;
    std::size_t end = 0;
    int left = -1;
    int right = -1;
  };

  int BuildNode(std::size_t begin, std::size_t end);
  Surface::Triangulation mesh_;
  std::vector<std::size_t> faceOrder_;
  std::vector<Node> nodes_;
};

// Deterministic quadrature of a uniformly radiating apparent solar disc.
// Positions/distances are metres in one Cartesian frame; the caller supplies
// an actual point-to-Sun vector, not a Moon-centred direction. A zero opaque
// radius explicitly disables Earth occultation. Radial samples are uniform
// in projected disc area, with a symmetric azimuth grid. Resolution is always
// supplied by the input/fixture: no scientific angular resolution is hidden.
struct SolarDiscRequest {
  std::array<double,3> surfacePositionM{{0.0,0.0,0.0}};
  std::array<double,3> outwardNormal{{0.0,0.0,1.0}};
  std::array<double,3> pointToSunM{{0.0,0.0,0.0}};
  std::array<double,3> pointToEarthM{{0.0,0.0,0.0}};
  double solarRadiusM=0.0;
  double opaqueEarthRadiusM=0.0;
  double rayOriginOffsetM=0.0;
  std::size_t sourceFace=std::numeric_limits<std::size_t>::max();
  int radialPoints=0;
  int azimuthPoints=0;
  // Surface irradiance uses the facet's self-horizon. Elevated-observer sky
  // visibility instead lets the actual terrain intersections define the
  // horizon; it must not discard a Sun below the local horizontal plane.
  bool applyLocalSelfHorizon=true;
  // Gas absorption starts at an actual volume point rather than an offset
  // facet. This opt-in mode requires self-horizon OFF, offset exactly zero
  // and the no-face sentinel: a real blocking facet must never be skipped;
  // surface callers retain their positive-offset/self-intersection contract.
  // The dummy normal is used only by surface-irradiance diagnostics, never by
  // the gas solid-angle visibility returned below.
  bool volumePointOrigin=false;
};

struct SolarDiscSample {
  std::size_t rayCount=0;
  double terrainVisibleAreaFraction=0.0;
  double earthVisibleAreaFraction=0.0;
  double jointlyVisibleAreaFraction=0.0;
  // Scalar absorption of uniform radiance integrates dOmega, whereas the
  // area fractions above integrate projected disc area. The equal-area ray
  // grid therefore needs a 1/cos(rho) Jacobian and a full-disc normalization.
  // These fractions have no facet incidence cosine or surface absorptivity.
  // They are quadrature estimates; finite grids can miss a thin visible sliver.
  double terrainVisibleSolidAngleFraction=0.0;
  double earthVisibleSolidAngleFraction=0.0;
  double jointlyVisibleSolidAngleFraction=0.0;
  double unobscuredSolidAngleSr=0.0; // Exact full apparent-disc solid angle [sr].
  // Incident irradiance divided by irradiance on a perpendicular plane.
  // This includes the positive local incidence cosine for EACH visible ray.
  double projectedIncidence=0.0;
  double absorbedProjectedIncidence=0.0;
  // A finite grid cannot establish that an arbitrarily small visible sliver
  // is absent. These names deliberately retain that numerical limitation.
  bool sampledAnyTerrainVisible=false;
  bool sampledAnyJointlyVisible=false;
};

// Actual fixed-site observer geometry, not a body-centre aberration followed
// by geometric parallax. et is the reception epoch [TDB s past J2000]; pointM
// is constant in MOON_ME_DE421 relative to MOON. Returned Sun/Earth vectors
// are observer-to-target metres on that frame's axes at the OBSERVER epoch.
// Explicit reception corrections: NONE, LT, LT+S, CN or CN+S. Kernels must
// already be loaded. This belongs to the application's serialized SPICE
// phase (CSPICE global error/kernel state is not a thread-local service).
// Invalid/missing-SPICE/runtime input leaves BOTH distinct outputs unchanged.
bool ApparentSolarSystemGeometryAtPoint(double et,
    const std::array<double,3>& pointM,const std::string& abcorr,
    std::array<double,3>* sunRelativeM,std::array<double,3>* earthRelativeM,
    std::string* error);

// An empty absorptivity callback means unity. Otherwise it is evaluated at
// each positive incidence cosine and must return a finite value in [0,1].
// Production passes the SAME directional Bond-albedo kernel used by Column,
// so a partially visible disc cannot acquire a central-ray albedo surrogate.
// A null tracer supplies the planar self-horizon limit; normal lunar terrain
// runs pass their initialized production BVH. Failure leaves output unchanged.
bool EvaluateSolarDisc(const TerrainRayTracer* terrain,
    const SolarDiscRequest& request, SolarDiscSample* output,
    std::string* error,
    const std::function<double(double)>& absorptivity={});

// Equal-epoch averaging only. D02 occurrence, mean visible area and incident
// flux are different quantities and are never substituted for one another.
// Epoch generation and endpoint/time-system interpretation belong to the
// explicit external contract, not this purely numerical accumulator.
struct TemporalVisibility {
  std::size_t epochCount=0;
  std::size_t sampledAnyVisibleCount=0;
  double meanVisibleAreaFraction=0.0;
  double meanProjectedIncidence=0.0;
  bool Add(const SolarDiscSample& sample,std::string* error);
  double SampledAnyDiscOccurrenceFraction() const;
};

// Process-wide production terrain state.  PrepareProductionSurface installs
// the exact triangulation generated from D01 on every MPI rank.  These hooks
// match PIC::Mover's optional moving-surface ABI and therefore exercise the
// normal particle mover rather than a test-only propagator.
bool InitializeProductionTerrain(const Surface::Triangulation& triangulation,
                                 const Runtime::Configuration& configuration,
                                 std::string* error);
bool ProductionTerrainInitialized();
const TerrainRayTracer* ProductionTerrain();
// Optional P08 artificial material collector. node is the actual PIC AMR node
// as an opaque pointer so geometry-only header users need no PIC/generated
// headers. The explicit resolver identifies EACH incident catalog state; no
// compiled-species-to-ground-state fallback exists. Installation validates
// compiled species/mass/charge and is allowed only for the opt-in selector.
using IonCollectorStateResolver=bool(*)(long int particle,void* node,
    void* user,std::size_t* state,std::string* error);
bool InstallProductionIonCollector(
    const std::vector<Species::IonCollector::StateDefinition>& states,
    IonCollectorStateResolver resolver,void* user,
    Species::IonCollector::Ledger* ledger,std::string* error);
// Physical source owners supply a source-bound unit explicitly. The overload
// above forwards the unchanged compiled legacy constant; no SI replacement is
// hidden in a geometry-only header or written into shared core constants.
bool InstallProductionIonCollector(
    const std::vector<Species::IonCollector::StateDefinition>& states,
    IonCollectorStateResolver resolver,void* user,
    Species::IonCollector::Ledger* ledger,std::string* error,
    double elementaryChargeC,const std::string& elementaryChargeSourceIdentity);
bool ProductionIonCollectorInstalled();
// Call at the ACTUAL accepted sphere/facet collision BEFORE querying neutral
// temperature/sticking. false means a neutral and preserves its old path.
// Unsupported/incomplete ion physics terminates with a detailed diagnostic.
// true means the transfer is recorded: caller returns the mover's Delete code,
// leaving the mover as the SOLE owner of actual particle deletion/list repair.
bool HandleProductionIonImpact(long int particle,const double* positionLsoM,
    const double* velocityLsoMS,void* node,
    const Species::IonCollector::ImpactFrame& frame);
void BeginProductionParticleStep(double ephemerisTimeBeginS);
// Fixed source-CDF/hazard epoch [TDB s] for the current operator-split
// particle interval. Unlike ProductionInjectionEpochS, it never includes an
// individual particle's departure offset.
double ProductionSourceEpochS();
// Actual source-departure epoch [TDB s] when the scoped PIC injection context
// is active. Explicit common-time uses the advancing TRANSPORT origin/step;
// independent-history retains the legacy external-step epoch plus offset.
// This must not be substituted for the frozen coefficient/geometry epoch:
// those providers use ProductionSourceEpochS. A common birth outside an
// active physical step is fatal rather than reusing a stale/frozen birth.
double ProductionInjectionEpochS();
// Shared explicit new-tracer/frozen-body policy. Preserve both rotation
// blocks; remove only dR/dt*x material velocity. Legacy matrices are untouched.
bool ApplyFrozenTracingBodyTransform(const Runtime::Configuration& configuration,
    double* stateTransform36,std::string* error);
bool ProductionInjectionBodyToLsoTransform(double* stateTransform36,
    std::string* error);
bool FindProductionMovingSurfaceIntersection(long int particle,
    const double* segmentPositionLsoM,const double* segmentVelocityLsoMPerS,
    const double* segmentAccelerationLsoMPerS2,double elapsedFlightTimeS,
    double maximumFlightTimeS,double& hitFlightTimeS,
    double* hitPositionLsoM,double* outwardNormalLso,
    double* surfaceVelocityLsoMPerS);
bool LastProductionCollisionFrame(double* outwardNormalBodyFixed,
                                  double* lsoToBodyStateTransform36,
                                  double* bodyToLsoStateTransform36);
// NEW tracer accepted-event transaction. The finder is a pure preview; only
// this callback commits the actual selected normal/frame and one collision.
bool CommitProductionMovingSurfaceCollision(long int particle,
    const double* selectedNormalLso,double* normalBody,
    double* lsoToBody36,double* bodyToLso36,std::string* error);

// Monotonic process-local diagnostic count of geometrically resolved moving
// terrain intersections. Legacy mode retains finder-report counting; new
// tracers count accepted interaction commits, never adaptive previews. Linked
// I14 snapshots it before/after a mover fixture to prove that the normal mover
// traversed the production LOLA callback rather than a test-only ray tracer.
std::uint64_t ProductionCollisionCount();

}  // namespace SurfaceEnvironment
}  // namespace Moon

#endif  // SRCMOON_MOON_SURFACE_ENVIRONMENT_H

#ifndef SHIELDSIM_CASE_DEFINITION_HH
#define SHIELDSIM_CASE_DEFINITION_HH

/*
 * Strict, versioned case-v1/core and case-v2/Stage-3/4/5/6 run definition.
 *
 * All values below are stored in canonical engineering units (mm, MeV,
 * g/cm2).  Unit conversion happens once while parsing.  Geometry builders and
 * transport actions therefore never need to guess the unit or apply an
 * undocumented default. Case v2 additionally carries explicit local-to-world
 * transforms and exact response/direct-validation incident nodes; it permits
 * raw transport tallies plus the bounded Stage-4 TID, Stage-5 silicon damage,
 * and Stage-6 LET/hardness observables. Each remains inside its own explicit
 * species/material/model validation gate.
 */

#include "MaterialLibrary.hh"
#include "SpectrumDataset.hh"
#include "StrictJson.hh"
#include "GeometryTransform.hh"
#include "DamageResponseLibrary.hh"

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace shield_sim {

struct ShieldLayerDefinition {
  std::string id;
  std::string materialId;
  std::string materialHash;
  double densityGPerCm3 = 0.0;
  std::string inputThicknessType;
  double thicknessMm = 0.0;
  double arealDensityGPerCm2 = 0.0;
};

// Stage-10's first non-slab provider deliberately uses a small, auditable CSG
// vocabulary.  Nodes are immutable after parsing; every dimensional member is
// already converted to millimetres and every transform is active/local-to-
// parent.  Primitive members irrelevant to a node type remain zero.
struct AnalyticCsgNodeDefinition {
  std::string nodeId;
  std::string type; // box, sphere, cylinder, union, intersection, subtraction
  double sizeXmm = 0.0;
  double sizeYmm = 0.0;
  double sizeZmm = 0.0;
  double innerRadiusMm = 0.0;
  double outerRadiusMm = 0.0;
  double lengthMm = 0.0;
  std::string relationship;
  std::shared_ptr<AnalyticCsgNodeDefinition> left;
  std::shared_ptr<AnalyticCsgNodeDefinition> right;
  RigidTransform rightTransform;
  Vector3d boundsMinimumMm;
  Vector3d boundsMaximumMm;
  Json canonicalJson;
};

struct AnalyticCsgRegionDefinition {
  std::string id;
  std::string role; // exactly one "sensitive"; remaining regions are passive
  std::string materialId;
  std::string materialHash;
  double densityGPerCm3 = 0.0;
  RigidTransform transform;
  std::shared_ptr<AnalyticCsgNodeDefinition> solid;
  double expectedVolumeCm3 = 0.0;
  Vector3d expectedCentroidMm;
  double massRelativeTolerance = 0.0;
  std::string massPropertiesProvenance;
  std::string entrySurface;
  Json canonicalJson;
};

// A bounded Stage-10 surface-mesh contract. Provider v1 accepts only closed,
// consistently outward-oriented convex polyhedra. That restriction makes
// watertightness, orientation, self-intersection freedom, volume, centroid,
// and surface area deterministically checkable without an implicit repair or
// an unqualified general-purpose CAD kernel.
struct TessellatedMeshFacetDefinition {
  std::string id;
  std::array<std::size_t,3> vertexIndices{{0,0,0}};
};

struct TessellatedMeshRegionDefinition {
  std::string id;
  std::string role;
  std::string materialId;
  std::string materialHash;
  double densityGPerCm3 = 0.0;
  RigidTransform transform;
  std::vector<Vector3d> verticesMm;
  std::vector<TessellatedMeshFacetDefinition> facets;
  double validationToleranceMm = 0.0;
  double exactVolumeCm3 = 0.0;
  double exactSurfaceAreaCm2 = 0.0;
  Vector3d exactCentroidMm;
  double expectedVolumeCm3 = 0.0;
  Vector3d expectedCentroidMm;
  double massRelativeTolerance = 0.0;
  std::string massPropertiesProvenance;
  std::string entrySurfaceType;
  std::vector<std::string> entryFacetIds;
  Json canonicalJson;
};

struct SlabAbsorberDefinition {
  std::string id;
  std::string materialId;
  std::string materialHash;
  double densityGPerCm3 = 0.0;
  double sizeXmm = 0.0;
  double sizeYmm = 0.0;
  double thicknessMm = 0.0;
  // Reporting-only local-z extent of the one sensitive region. It equals
  // thicknessMm for slab.v1 but is distinct from an analytic provider's
  // overall conservative z bound when passive overburden is present.
  double scoringThicknessMm = 0.0;
  double frontZmm = 0.0;
  double centerXmm = 0.0;
  double centerYmm = 0.0;
  // In case-v2 the local origin is the center of the upstream slab face and
  // local +z points through the slab.  The legacy scalar fields above remain
  // populated so old reporting code can read v1 cases during migration.
  RigidTransform transform;
  // `slab.v1` preserves all Stage 1--9 behaviour. `analytic_csg.v1` and
  // `tessellated_mesh.v1` store one or more non-overlapping material regions
  // below. The legacy scalar
  // material members describe the sole sensitive region, while size/thickness
  // are conservative full-provider bounds used only for containment. The
  // current one-sensitive-region limitation is explicit at the parser boundary.
  std::string provider = "slab.v1";
  std::string providerVersion = "1.0.0";
  std::string geometrySourceId;
  std::string geometrySourceVersion;
  std::string geometrySourceProvenance;
  std::string geometryContentHash;
  std::vector<AnalyticCsgRegionDefinition> csgRegions;
  std::vector<TessellatedMeshRegionDefinition> meshRegions;
  std::string sensitiveRegionId;
  Vector3d localBoundsMinimumMm;
  Vector3d localBoundsMaximumMm;
  Json canonicalGeometry;
  std::string geometryHash;
};

struct ResponseIncidentDefinition {
  bool enabled = false;
  std::string incidentId;
  std::string particleName;
  int atomicNumber = 0;
  int massNumber = 0;
  std::string energyCoordinate;
  double energyCoordinateValue = 0.0;
  std::string energyCoordinateUnit;
  double totalKineticEnergyMeV = 0.0;
  double mu = 1.0;
  double phiRad = 0.0;
};

struct WeightedResponseIncident {
  ResponseIncidentDefinition incident;
  double weight = 0.0;
};

// Stage-6 LET scoring is deliberately configured in the case instead of
// hidden in plotting code.  The resolved edges are retained so every result
// can be interpreted without knowing a parser default. Values use linear LET
// in keV/um; the output also reports the density-dependent mass stopping-power
// representation in MeV cm2/g.
struct LETResponseDefinition {
  bool enabled = false;
  std::string definitionId;
  double minimumKeVPerUm = 0.0;
  double maximumKeVPerUm = 0.0;
  int binsPerDecade = 0;
  std::vector<double> edgesKeVPerUm;
  std::string outOfRangePolicy;
};

// One standard-uncertainty term in a combined Stage-4 budget. Keeping the
// category and provenance with the number prevents an aggregate percentage
// from concealing whether source, material/mass, numerical, or physics-
// reference uncertainty was omitted.
struct SystematicUncertaintyComponent {
  std::string id;
  std::string category;
  double relativeStandardUncertainty = 0.0;
  std::string provenance;
};

// Correlation is declared separately because the same component values can
// combine differently when, for example, material density and absorber mass
// share a calibration. A complete Stage-4 budget requires every unordered
// component pair to be explicit, including pairs with zero correlation.
struct SystematicCorrelation {
  std::string firstComponent;
  std::string secondComponent;
  double coefficient = 0.0;
  std::string rationale;
};

struct AccuracyRequirement {
  double target = 0.10;
  std::string scope = "combined";
  double standardK = 1.0;
  double standardProbability = 0.68;
  double intervalProbability = 0.95;
  std::string intervalMethod = "student_t_or_exact_bounded";
  // Case-v1 retains the aggregate field for compatibility. A case-v2 TID
  // budget is complete only when the component and all-pairs correlation
  // records below have been checked and combined by the parser.
  bool systematicBudgetComplete = false;
  double systematicRelativeUncertainty = 0.0;
  std::vector<SystematicUncertaintyComponent> systematicComponents;
  std::vector<SystematicCorrelation> systematicCorrelations;
  std::string systematicPairPolicy;
  bool hasKnownSystematicFloor = false;
  double knownSystematicFloor = 0.0;
};

struct AccuracyPlan {
  AccuracyRequirement global;
  std::map<std::string, AccuracyRequirement> quantityOverrides;
  long minHistories = 1000;
  long maxHistories = 10000;
  long batchSize = 1000;
  bool targetWasDefaulted = false;
};

class CaseDefinition {
public:
  static std::shared_ptr<CaseDefinition> Load(const std::string& casePath);

  std::string casePath;
  std::string schemaVersion;
  std::string caseId;
  std::string purpose;
  std::shared_ptr<MaterialLibrary> materials;
  std::vector<ShieldLayerDefinition> shieldLayers;
  std::vector<SlabAbsorberDefinition> absorbers;
  SpectrumDataset spectrum;
  // The resolved source format and selected component are archived separately
  // from the adapter's full canonical metadata.  This makes it unambiguous
  // which one of a multi-component published dataset entered Geant4.
  std::string sourceFormat;
  std::string sourceComponentId;

  std::string coordinateFrame;
  double worldSizeXmm = 0.0, worldSizeYmm = 0.0, worldSizeZmm = 0.0;
  double worldCenterXmm = 0.0, worldCenterYmm = 0.0, worldCenterZmm = 0.0;
  double sourceSizeXmm = 0.0, sourceSizeYmm = 0.0, sourceZmm = 0.0;
  double sourceCenterXmm = 0.0, sourceCenterYmm = 0.0;
  RigidTransform sourceTransform;
  std::string sourceMode;
  double shieldSizeXmm = 0.0, shieldSizeYmm = 0.0, shieldFrontZmm = 0.0;
  double shieldCenterXmm = 0.0, shieldCenterYmm = 0.0;
  RigidTransform shieldTransform;

  double protonMinMeV = 0.0, protonMaxMeV = 0.0;
  double alphaMinMeV = 0.0, alphaMaxMeV = 0.0;
  std::vector<std::string> quantities;
  LETResponseDefinition letResponse;
  // Non-null exactly when a Stage-5 DDD or n_eq quantity was requested. The
  // manifest and all referenced tables are authenticated during case parsing.
  std::shared_ptr<DamageResponseLibrary> damageResponse;
  ResponseIncidentDefinition responseIncident;
  std::vector<WeightedResponseIncident> directIncidents;
  AccuracyPlan accuracy;

  std::string physicsList;
  std::string ionPhysics = "reference_list_default";
  // The reference-list default is the production baseline. Option 4 is an
  // explicitly archived electromagnetic-model sensitivity, never an implicit
  // replacement selected by the executable.
  std::string emPhysics = "reference_list_default";
  double productionCutMm = -1.0;
  double maxStepMm = -1.0;
  long randomSeed = 0;
  std::string outputPrefix;

  Json resolvedJson;
  std::string geometryHash;
  std::string caseHash;

  double ShieldTotalThicknessMm() const;
  double ShieldTotalArealDensityGPerCm2() const;
  const AccuracyRequirement& AccuracyFor(const std::string& quantity) const;
};

} // namespace shield_sim

#endif

#ifndef SHIELDSIM_CASE_DEFINITION_HH
#define SHIELDSIM_CASE_DEFINITION_HH

/*
 * Strict, versioned case-v1/core and case-v2/Stage-3 run definition.
 *
 * All values below are stored in canonical engineering units (mm, MeV,
 * g/cm2).  Unit conversion happens once while parsing.  Geometry builders and
 * transport actions therefore never need to guess the unit or apply an
 * undocumented default. Case v2 additionally carries explicit local-to-world
 * transforms and exact response/direct-validation incident nodes; it permits
 * only raw transport tallies until later quantity gates are completed.
 */

#include "MaterialLibrary.hh"
#include "SpectrumDataset.hh"
#include "StrictJson.hh"
#include "GeometryTransform.hh"

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

struct SlabAbsorberDefinition {
  std::string id;
  std::string materialId;
  std::string materialHash;
  double densityGPerCm3 = 0.0;
  double sizeXmm = 0.0;
  double sizeYmm = 0.0;
  double thicknessMm = 0.0;
  double frontZmm = 0.0;
  double centerXmm = 0.0;
  double centerYmm = 0.0;
  // In case-v2 the local origin is the center of the upstream slab face and
  // local +z points through the slab.  The legacy scalar fields above remain
  // populated so old reporting code can read v1 cases during migration.
  RigidTransform transform;
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

struct AccuracyRequirement {
  double target = 0.10;
  std::string scope = "combined";
  double standardK = 1.0;
  double standardProbability = 0.68;
  double intervalProbability = 0.95;
  std::string intervalMethod = "student_t_or_exact_bounded";
  // A complete budget supplies the combined standard relative uncertainty.
  // An incomplete budget may still carry a known lower floor; if that floor
  // already exceeds the target, additional Monte Carlo histories cannot help.
  bool systematicBudgetComplete = false;
  double systematicRelativeUncertainty = 0.0;
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
  ResponseIncidentDefinition responseIncident;
  std::vector<WeightedResponseIncident> directIncidents;
  AccuracyPlan accuracy;

  std::string physicsList;
  std::string ionPhysics = "reference_list_default";
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

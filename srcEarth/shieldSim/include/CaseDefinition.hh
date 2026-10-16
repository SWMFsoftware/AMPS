#ifndef SHIELDSIM_CASE_DEFINITION_HH
#define SHIELDSIM_CASE_DEFINITION_HH

/*
 * Strict, versioned Stage-1 run definition.
 *
 * All values below are stored in canonical engineering units (mm, MeV,
 * g/cm2).  Unit conversion happens once while parsing.  Geometry builders and
 * transport actions therefore never need to guess the unit or apply an
 * undocumented default.
 */

#include "MaterialLibrary.hh"
#include "SpectrumDataset.hh"
#include "StrictJson.hh"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace shield_sim {

struct ShieldLayerDefinition {
  std::string id;
  std::string materialId;
  std::string materialHash;
  std::string inputThicknessType;
  double thicknessMm = 0.0;
  double arealDensityGPerCm2 = 0.0;
};

struct SlabAbsorberDefinition {
  std::string id;
  std::string materialId;
  std::string materialHash;
  double sizeXmm = 0.0;
  double sizeYmm = 0.0;
  double thicknessMm = 0.0;
  double frontZmm = 0.0;
  double centerXmm = 0.0;
  double centerYmm = 0.0;
  std::string geometryHash;
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
  std::string caseId;
  std::string purpose;
  std::shared_ptr<MaterialLibrary> materials;
  std::vector<ShieldLayerDefinition> shieldLayers;
  std::vector<SlabAbsorberDefinition> absorbers;
  SpectrumDataset spectrum;

  std::string coordinateFrame;
  double worldSizeXmm = 0.0, worldSizeYmm = 0.0, worldSizeZmm = 0.0;
  double worldCenterXmm = 0.0, worldCenterYmm = 0.0, worldCenterZmm = 0.0;
  double sourceSizeXmm = 0.0, sourceSizeYmm = 0.0, sourceZmm = 0.0;
  double sourceCenterXmm = 0.0, sourceCenterYmm = 0.0;
  std::string sourceMode;
  double shieldSizeXmm = 0.0, shieldSizeYmm = 0.0, shieldFrontZmm = 0.0;
  double shieldCenterXmm = 0.0, shieldCenterYmm = 0.0;

  double protonMinMeV = 0.0, protonMaxMeV = 0.0;
  double alphaMinMeV = 0.0, alphaMaxMeV = 0.0;
  std::vector<std::string> quantities;
  AccuracyPlan accuracy;

  std::string physicsList;
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

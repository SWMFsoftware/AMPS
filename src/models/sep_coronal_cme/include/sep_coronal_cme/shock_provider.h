#ifndef SEP_CORONAL_CME_SHOCK_PROVIDER_H
#define SEP_CORONAL_CME_SHOCK_PROVIDER_H

#include "sep_coronal_cme/mhd_jump_solver.h"
#include "sep_coronal_cme/model_configuration.h"
#include "sep_coronal_cme/source_surface_coupling.h"
#include "sep_status.h"

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

// The initial gate controls only whether a prepared zero-source surface may be
// published.  It does not turn a sub-fast patch into a shock and therefore
// cannot manufacture injection at the beginning of a run.
enum class InitialShockGate {
  None,
  AnyFastPatch,
  MinimumFastAreaFraction
};

// Every one-sided patch carries the exact interface evidence used to prepare
// its upstream state.  A policy label alone is insufficient: production must
// also retain the generation, provenance, topology, kinematic, and balance
// results which qualified that policy.
struct PatchInterfaceEvidence {
  InterfacePolicy policy = InterfacePolicy::DiagnosticKinematic;
  std::uint64_t backgroundGeneration = 0;
  bool provenanceComplete = false;
  bool topologyPassed = false;
  bool kinematicPassed = false;
  bool balancePassed = false;
  std::string interfaceIdentity;
};

struct ShockPatchInput {
  std::uint64_t stableId = 0;
  std::uint64_t parentId = 0;
  double areaM2 = 0.0;
  Vec3 outwardNormal;
  MhdPrimitiveState upstream;
  double shockNormalSpeedMPerS = 0.0;
  double incidentNumberRatePerS = 0.0;
  double incidentKineticEnergyRateW = 0.0;
  bool intersectsTransitionClearance = false;
  bool sourceTerminated = false;
  PatchInterfaceEvidence interface;
};

struct ShockPatchSnapshot {
  std::uint64_t stableId = 0;
  std::uint64_t parentId = 0;
  double areaM2 = 0.0;
  double fastMach = 0.0;
  bool geometric = true;
  bool fast = false;
  bool supercritical = false;
  bool sourceEligibleBeforeClearance = false;
  bool sourceActive = false;
  bool sourceTerminated = false;
  bool excludedByTransitionClearance = false;
  double physicalNumberRatePerS = 0.0;
  double physicalKineticEnergyRateW = 0.0;
  std::string rejectionReason;
  PatchInterfaceEvidence interface;
  MhdShockSolution jump;
};

struct ShockSurfaceMeasures {
  double geometricAreaM2 = 0.0;
  double fastAreaM2 = 0.0;
  double supercriticalAreaM2 = 0.0;
  double sourceActiveAreaM2 = 0.0;
  double sourceTerminatedAreaM2 = 0.0;
  double counterfactualAreaM2 = 0.0;
  double excludedAreaM2 = 0.0;
  double counterfactualNumberRatePerS = 0.0;
  double excludedNumberRatePerS = 0.0;
  double counterfactualKineticEnergyRateW = 0.0;
  double excludedKineticEnergyRateW = 0.0;
  MeasureValidity areaFractionValidity = MeasureValidity::Valid;
  MeasureValidity numberFractionValidity = MeasureValidity::Valid;
  MeasureValidity energyFractionValidity = MeasureValidity::Valid;
};

struct ShockSurfaceSnapshot {
  std::uint64_t generation = 0;
  std::uint64_t backgroundGeneration = 0;
  double timeS = 0.0;
  bool eventGradeable = true;
  std::vector<ShockPatchSnapshot> patches;
  ShockSurfaceMeasures measures;
};

struct ShockPreparationOptions {
  bool productionIntent = false;
  bool requireSupercritical = false;
  double criticalFastMach = 0.0;
  InitialShockGate initialGate = InitialShockGate::None;
  double minimumFastAreaFraction = 0.0;
  double maximumTransitionAreaFraction = 1.0;
  double maximumTransitionNumberFraction = 1.0;
  double maximumTransitionEnergyFraction = 1.0;
};

// Split a parent patch into deterministic one-sided children.  Fractions are
// physical area fractions and must close to one; rates are divided with the
// same fractions so geometry and both source measures remain conservative.
Core::Result<std::vector<ShockPatchInput>> SplitShockPatch(
    const ShockPatchInput& parent, const std::vector<double>& areaFractions,
    const std::vector<PatchInterfaceEvidence>& childEvidence);

// Locate the first M_f=1 crossing.  The supplied function must be continuous
// on the bracket; no update-cadence sample is used as the event time.
Core::Result<double> LocateFirstFastCrossing(
    double beginTimeS, double endTimeS,
    const std::function<double(double)>& fastMach, double timeToleranceS,
    int maximumIterations = 100);

class ShockProvider {
 public:
  virtual ~ShockProvider() = default;
  virtual std::shared_ptr<const ShockSurfaceSnapshot> PreparedSurface() const = 0;
};

// Preparation is a two-phase transaction: all patches and integral ledgers
// are built in a private candidate, then one immutable shared snapshot is
// published.  Returning failure leaves CurrentGeneration() and every owning
// handle obtained earlier unchanged.
class TransactionalShockProvider final : public ShockProvider {
 public:
  Core::Result<std::shared_ptr<const ShockSurfaceSnapshot>> Prepare(
      double timeS, std::uint64_t backgroundGeneration,
      const std::vector<ShockPatchInput>& patches,
      const ShockPreparationOptions& options);

  std::shared_ptr<const ShockSurfaceSnapshot> PreparedSurface() const override {
    return current_;
  }
  std::uint64_t CurrentGeneration() const noexcept {
    return current_ ? current_->generation : 0;
  }

 private:
  std::shared_ptr<const ShockSurfaceSnapshot> current_;
};

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_SHOCK_PROVIDER_H

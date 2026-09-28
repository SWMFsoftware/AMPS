#ifndef SEP_CORONAL_CME_MHD_JUMP_SOLVER_H
#define SEP_CORONAL_CME_MHD_JUMP_SOLVER_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

struct MhdPrimitiveState {
  double massDensityKgM3 = 0.0;
  double pressurePa = 0.0;
  Vec3 velocityMPerS;
  Vec3 magneticFieldT;
};

struct MhdCharacteristicSpeeds {
  double soundSpeedMPerS = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double normalAlfvenSpeedMPerS = 0.0;
  double slowSpeedMPerS = 0.0;
  double fastSpeedMPerS = 0.0;
  double obliquityRad = 0.0;
  double plasmaBeta = 0.0;
};
Core::Result<MhdCharacteristicSpeeds> EvaluateMhdCharacteristics(
    const MhdPrimitiveState& state, Vec3 normal, double gammaAdiabatic);

struct RankineHugoniotResiduals {
  double mass = 0.0;
  double normalMagnetic = 0.0;
  double tangentialElectric = 0.0;
  double normalMomentum = 0.0;
  double tangentialMomentum = 0.0;
  double totalEnergy = 0.0;
  double maximum = 0.0;
};

struct MhdShockSolution {
  MhdPrimitiveState downstream;
  MhdCharacteristicSpeeds upstreamCharacteristics;
  MhdCharacteristicSpeeds downstreamCharacteristics;
  double compressionRatio = 0.0;
  double upstreamFastMach = 0.0;
  double upstreamAlfvenMach = 0.0;
  double entropyLogIncrement = 0.0;
  RankineHugoniotResiduals residuals;
  std::string branch;
};

// Solves the regular compressive fast branch.  Failure returns no downstream
// state: a sub-fast front, inadmissible root, or residual failure cannot leak a
// partially constructed primitive state to a source provider.
Core::Result<MhdShockSolution> SolveObliqueFastShock(
    const MhdPrimitiveState& upstream, Vec3 outwardNormal,
    double shockNormalSpeedMPerS, double gammaAdiabatic,
    double residualTolerance = 1.0e-9);

enum class MachConvention { Fast, TotalAlfven, NormalAlfven };

struct CriticalMachTable {
  std::string version;
  std::string checksum;
  double gammaAdiabatic = 0.0;
  MachConvention convention = MachConvention::Fast;
  std::vector<double> betaGrid;
  std::vector<double> obliquityGridRad;
  // Row-major beta x obliquity values.
  std::vector<double> values;
};

enum class CriticalQueryValidity { Valid, BetaOutsideCoverage,
                                   NormalAlfvenInapplicable };
struct CriticalMachQuery {
  CriticalQueryValidity validity = CriticalQueryValidity::Valid;
  double value = 0.0;
  std::string reason;
};
Core::Result<CriticalMachQuery> QueryCriticalMach(
    const CriticalMachTable& table, double beta, double obliquityRad,
    double gammaAdiabatic, MachConvention requestedConvention,
    bool exactNormalFieldIsZero = false);

enum class CriticalityMissPolicy { DiagnosticOnlyFast,
                                   FailPreflight,
                                   ExcludeSourceBudgeted };
struct CriticalityDecision {
  bool sourceEligible = false;
  bool criticalityAvailable = false;
  bool supercritical = false;
  bool excludedWithoutRenormalization = false;
};
Core::Result<CriticalityDecision> ApplyCriticalityPolicy(
    double matchingMachNumber, const CriticalMachQuery& query,
    bool sourceRequiresSupercritical, CriticalityMissPolicy policy);

struct CriticalCoverageSample {
  double beta = 0.0;
  double obliquityRad = 0.0;
  double areaM2 = 0.0;
  double incidentNumberRatePerS = 0.0;
  double incidentKineticEnergyRateW = 0.0;
  bool exactNormalFieldIsZero = false;
};
struct CriticalCoverageLedger {
  double candidateAreaM2 = 0.0;
  double unavailableAreaM2 = 0.0;
  double candidateNumberRatePerS = 0.0;
  double unavailableNumberRatePerS = 0.0;
  double candidateEnergyRateW = 0.0;
  double unavailableEnergyRateW = 0.0;
  double areaFraction = 0.0;
  double numberFraction = 0.0;
  double energyFraction = 0.0;
  bool noCandidateSupport = false;
};
Core::Result<CriticalCoverageLedger> PreflightCriticalCoverage(
    const CriticalMachTable& table,
    const std::vector<CriticalCoverageSample>& samples,
    double gammaAdiabatic, MachConvention convention,
    double maximumAreaFraction, double maximumNumberFraction,
    double maximumEnergyFraction);

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_MHD_JUMP_SOLVER_H

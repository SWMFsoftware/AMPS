#ifndef SEP_COMMON_MEAN_FREE_PATH_MEAN_FREE_PATH_H
#define SEP_COMMON_MEAN_FREE_PATH_MEAN_FREE_PATH_H

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace SEP {
namespace MeanFreePath {

// Public constants are exact definitions used by the specification.  The
// library accepts and returns SI values; these constants are supplied only so
// adapters and tests can make explicit, one-time conversions at the boundary.
constexpr double SpeedOfLightMPerS = 299792458.0;
constexpr double AstronomicalUnitM = 149597870700.0;
constexpr double ElementaryChargeC = 1.602176634e-19;
constexpr double ElectronVoltJ = ElementaryChargeC;

enum class StatusCode {
  Success,
  InvalidParticle,
  MissingInput,
  InvalidInput,
  InvalidConfiguration,
  OutsideModelDomain,
  RequiresExternalInput,
  RequiresUserDecision,
  RequiresSourceOrCodeAudit,
  ReferenceDataOnly,
  BlockedDependency,
  UnsupportedModel,
  NumericalFailure
};

struct Status {
  StatusCode code = StatusCode::Success;
  std::string detail;
  bool ok() const { return code == StatusCode::Success; }
  static Status Success();
  static Status Error(StatusCode code, const std::string& detail);
};

// LambdaKind prevents the common but physically incorrect interchange of a
// field-parallel MFP, the SEP radial convention, a radial tensor projection,
// and a full-orbit isotropic-scattering length (specification Eqs. 1--3).
enum class LambdaKind {
  Parallel,
  RadialSEP,
  RadialTensor,
  IsotropicScattering,
  Unspecified
};

enum class MomentumVariable {
  RigidityV,
  MomentumPcEV,
  KineticTotalEV,
  KineticPerNucleonEV
};

enum class RuntimeState {
  ReadyExplicitInputs,
  ReadyPublishedVariant,
  ReadyPublishedKappaOnly,
  RequiresExternalInput,
  RequiresUserDecision,
  RequiresSourceOrCodeAudit,
  ReferenceDataOnly,
  BlockedDependency
};

enum class OperatorConvention {
  Standard,
  HalfD
};

enum class OutOfDomainPolicy {
  Error,
  WarnAndEvaluate
};

enum Diagnostic : std::uint32_t {
  DiagnosticNone = 0,
  Extrapolated = 1u << 0,
  DerivedQuantity = 1u << 1,
  PublishedVariant = 1u << 2,
  NumericalFloorApplied = 1u << 3,
  NominalAndIntegralLambdaDiffer = 1u << 4
};

struct ParticleState {
  // momentumKgMPerS is the magnitude of total relativistic momentum.  A
  // nucleon count is mandatory only for KineticPerNucleonEV; it is never
  // inferred from mass or charge.
  double massKg = 0.0;
  double chargeC = 0.0;
  double momentumKgMPerS = 0.0;
  std::optional<double> nucleonCount;
  std::string speciesId;
};

struct Kinematics {
  double gamma = 0.0;
  double beta = 0.0;
  double speedMPerS = 0.0;
  double rigidityV = 0.0;
  double momentumPcEV = 0.0;
  double kineticTotalEV = 0.0;
  std::optional<double> kineticPerNucleonEV;
};

Status ComputeKinematics(const ParticleState& particle, Kinematics* output);

struct LocalState {
  // All optional values are unavailable until supplied.  Zero is never used
  // as a synonym for a missing turbulence, field, shock, or activity input.
  std::optional<double> radiusM;
  std::optional<double> meanFieldMagnitudeT;
  std::optional<double> parkerSpiralAngleRad;
  std::optional<double> kappaPerpendicularM2PerS;
  std::optional<double> pitchAngleCosine;

  // Canonical slab inputs.  slabVarianceT2 is the total two-component slab
  // magnetic variance unless a model explicitly selects PerComponent.
  std::optional<double> slabVarianceT2;
  // Total Alfvén-wave magnetic variance [T^2] used by prescriptions such as
  // M-FLAMPA Eq. (30).  It is intentionally separate from slabVarianceT2:
  // an application may populate one from a provider only when that provider
  // establishes the corresponding physical decomposition/convention.
  std::optional<double> waveVarianceT2;
  std::optional<double> slabCorrelationLengthM;
  std::optional<double> kMinPerM;
  std::optional<double> kDPerM;
  std::optional<double> alfvenSpeedMPerS;

  // Signed-distance shock contract: upstreamDistanceM and downstreamDistanceM
  // are non-negative coordinates on their named side, never unsigned distance
  // to a surface.  A model validates only the side it consumes.
  enum class ShockSide { Unspecified, Upstream, Downstream };
  ShockSide shockSide = ShockSide::Unspecified;
  std::optional<double> upstreamDistanceM;
  std::optional<double> downstreamDistanceM;
  std::optional<double> upstreamFlowShockFrameMPerS;

  std::uint64_t backgroundRevision = 0;
  std::uint64_t turbulenceRevision = 0;
  std::string sampleIdentity;
};

struct LambdaValue {
  double metres = 0.0;
  LambdaKind kind = LambdaKind::Unspecified;
  std::optional<double> psiRad;
};

struct Provenance {
  std::string requestedModelId;
  std::string evaluatedModelId;
  std::string sourceKey;
  std::string sourceLocation;
  std::string formula;
  std::string configurationFingerprint;
  std::string sampleIdentity;
  std::uint64_t backgroundRevision = 0;
  std::uint64_t turbulenceRevision = 0;
  RuntimeState runtimeState = RuntimeState::RequiresSourceOrCodeAudit;
  bool derived = false;
};

struct Result {
  Status status;
  std::optional<LambdaValue> lambda;
  std::optional<double> kappaParallelM2PerS;
  std::optional<double> dMuMuPerS;
  std::optional<double> pitchAmplitudePerS;
  std::uint32_t diagnostics = DiagnosticNone;
  Provenance provenance;
};

struct ModelDescriptor {
  const char* stableId;
  RuntimeState declaredState;
  bool executable;
  const char* sourceKey;
  const char* sourceLocation;
  const char* directOutput;
  const char* note;
};

const std::vector<ModelDescriptor>& ModelRegistry();
const ModelDescriptor* FindModel(const std::string& stableId);

// Parser-neutral configuration.  The public maps are an immutable value once
// installed by SetActiveConfiguration.  They retain exact input text for
// restart/provenance, while numbers/choices contain the validated typed view.
struct Configuration {
  std::string modelId;
  std::map<std::string, std::string> raw;
  std::map<std::string, double> numbers;
  std::map<std::string, std::string> choices;
  LambdaKind outputKind = LambdaKind::Unspecified;
  MomentumVariable momentumVariable = MomentumVariable::RigidityV;
  OperatorConvention operatorConvention = OperatorConvention::Standard;
  OutOfDomainPolicy domainPolicy = OutOfDomainPolicy::Error;
  RuntimeState runtimeState = RuntimeState::RequiresSourceOrCodeAudit;
  std::string fingerprint;
};

struct InputParameter {
  std::string name;
  std::string value;
};

Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          Configuration* output);
Status ValidateConfiguration(const Configuration& configuration);
std::string ConfigurationFingerprint(const Configuration& configuration);

using ModelFunction = Result (*)(const ParticleState&, const LocalState&,
                                 const Configuration&);
extern ModelFunction ActiveModelFunction;
ModelFunction FunctionForModel(const std::string& modelId);
Result Evaluate(const ParticleState& particle, const LocalState& local,
                const Configuration& configuration);
Result EvaluateActive(const ParticleState& particle, const LocalState& local);
Status SetActiveConfiguration(const Configuration& configuration);
Status ConfigureActiveModel(const std::string& modelId,
                            const std::vector<InputParameter>& parameters);
Configuration GetActiveConfiguration();

Status EvaluateBatch(const std::vector<ParticleState>& particles,
                     const std::vector<LocalState>& states,
                     const Configuration& configuration,
                     std::vector<Result>* results);

// Geometry conversions are explicit operations.  No evaluator silently turns
// lambda_r into lambda_parallel, or drops the perpendicular term in Eq. (3).
Status ParallelToRadialSEP(const LambdaValue& parallel, double psiRad,
                           LambdaValue* radial);
Status RadialSEPToParallel(const LambdaValue& radial, double psiRad,
                           LambdaValue* parallel);
Status ParallelToRadialTensor(const LambdaValue& parallel,
                              double kappaPerpendicularM2PerS,
                              double speedMPerS, double psiRad,
                              LambdaValue* radialTensor);

// Pitch-angle helpers implement Eq. (11).  amplitude is D0 for q-form and
// nu0 for epsilon/isotropic forms.  Printed variants require their distinct
// model ID and are never selected by these helpers implicitly.
Status PitchAmplitudeFromLambda(const std::string& pitchModelId,
                                double q, double gapParameter,
                                double alfvenToParticleSpeed,
                                OperatorConvention op,
                                const LambdaValue& parallel,
                                double speedMPerS, double* amplitudePerS);
Status LambdaFromPitchAmplitude(const std::string& pitchModelId,
                                double q, double gapParameter,
                                double alfvenToParticleSpeed,
                                OperatorConvention op,
                                double amplitudePerS, double speedMPerS,
                                LambdaValue* parallel);
Status EvaluatePitchAngleDiffusion(double mu, const ParticleState& particle,
                                   const LocalState& local,
                                   const Configuration& configuration,
                                   double* dMuMuPerS);

// Stable He--Wan Eq. (43).  The input is x=lambda0/L.  Returning the ratio
// rather than reconstructing L avoids unnecessary overflow and is directly
// testable against F-NUM-01.
Status HeWanFocusingRatio(double x, double* ratio);

// Companion-data parsing helpers.  PublishedScalar never guesses at a range,
// uncertainty, inequality, or unit.  Exact rational syntax is accepted only
// when dimensionless=true.  Raw text is always retained for provenance.
enum class PublishedScalarKind {
  PlainNumber,
  ExactRational,
  NumberWithUnit,
  NonScalar
};
struct PublishedScalar {
  PublishedScalarKind kind = PublishedScalarKind::NonScalar;
  std::string raw;
  std::optional<double> value;
  std::string unit;
};
Status ParsePublishedScalar(const std::string& text, bool dimensionless,
                            PublishedScalar* output);

struct SourceRecord {
  std::string relativeFile;
  std::string recordId;
  std::string rawJson;
  std::string manifestSha256;
};
Status LoadSourceRecord(const std::string& bundleDirectory,
                        const std::string& relativeFile,
                        const std::string& recordId,
                        SourceRecord* output);

}  // namespace MeanFreePath
}  // namespace SEP

#endif

#ifndef SEP_COMMON_PARALLEL_DIFFUSION_PARALLEL_DIFFUSION_H
#define SEP_COMMON_PARALLEL_DIFFUSION_PARALLEL_DIFFUSION_H

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace SEP {
namespace ParallelDiffusion {

// The scientific contract is PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md v1.3.
// This API returns the parallel eigenvalue of the *symmetric spatial*
// diffusion tensor used by a Parker transport operator.  It does not return a
// radial or shock-normal projection, add particle drifts, construct the full
// tensor, or infer a momentum-diffusion coefficient.  Those operations need
// geometry and/or physics that are deliberately owned by the transport
// consumer (specification Equations (2), (3), (5), and (6)).
//
// All values crossing this boundary are SI.  In particular, momentum is the
// magnitude of the total particle momentum [kg m s^-1], charge is signed [C],
// rigidity is positive [V], magnetic field is [T], lambda_parallel is [m],
// and kappa_parallel is [m^2 s^-1].  Input-file adapters must perform any AU,
// GV, nT, MeV, or cm^2/s conversion exactly once before calling this API.
constexpr double SpeedOfLightMPerS = 299792458.0;

// Status reports whether a requested numerical value exists.  It is separate
// from Diagnostic: a finite coefficient can be numerically successful while
// a diagnostic warns that a physical approximation needs review.
enum class StatusCode {
  Success,
  // Invalid mass, charge, momentum, or a failed relativistic conversion.
  InvalidParticle,
  // A required local provider quantity is absent from its physical domain,
  // for example a non-positive magnitude for a selected magnetic field.
  InvalidBackground,
  // A model-selected optional input was not supplied at all.
  MissingInput,
  // Inputs are finite and well formed but outside the model's stated domain,
  // or the positive result cannot be represented in double precision.
  OutsideModelDomain,
  // Reserved for scattering models whose mathematical lambda diverges.
  InfiniteMeanFreePath,
  // Reserved for quadrature-based models in later roadmap stages.
  IntegrationFailed,
  // Reserved for coupled nonlinear closures in later roadmap stages.
  NonlinearSolverFailed,
  // Reserved for spectrum providers that violate their normalization.
  InconsistentSpectrum,
  // The selected model schema or parser-neutral key/value set is invalid.
  InvalidConfiguration,
  // The identifier is stable and known, but its backend is not implemented.
  UnsupportedModel
};

struct Status {
  StatusCode code = StatusCode::Success;
  std::string detail;

  bool ok() const { return code == StatusCode::Success; }
  static Status Success();
  static Status Error(StatusCode code, const std::string& detail);
};

// Bit mask values intentionally match the normative v1.3 interface sketch.
// PD01/PD02 currently sets DerivativeUnavailable for spatially varying laws;
// the other bits are reserved for the later physical backends that can
// actually establish those conditions.  Reserved bits must not be guessed
// from a finite scalar result.
enum Diagnostic : std::uint32_t {
  DiagnosticNone = 0,
  DiffusionLimitConcern = 1u << 0,
  StrongTurbulence = 1u << 1,
  QltWeakPerturbationConcern = 1u << 2,
  SurrogateErrorUnbounded = 1u << 3,
  FallbackUsed = 1u << 5,
  UserBoundApplied = 1u << 6,
  DerivativeUnavailable = 1u << 7
};

enum class ModelId {
  ConstantLambda,
  ConstantKappa,
  PowerLawLambda,
  BrokenRigidityKappa,
  QltSlabSpectrum,
  QltSlabInertial,
  PrescribedLambdaMuShape,
  BroadenedSlab,
  NlpaGivenPerp,
  NlgcE,
  NlgceN,
  NlgceF2014,
  TurbulenceAdapter,
  WaveSpectrumAdapter,
  Bohm,
  TabulatedParallel
};

struct ModelDescriptor {
  ModelId id;
  // Stable, case-insensitively parsed input spelling.  Never use an enum
  // ordinal as an input/restart identity.
  const char* stableId;
  // False means registered for an explicit UnsupportedModel result, not that
  // another model will be substituted.
  bool implemented;
  const char* directOutput;
  const char* implementationStage;
};

const std::vector<ModelDescriptor>& ModelRegistry();
const char* ModelName(ModelId model);
bool ParseModelId(const std::string& text, ModelId* model);

struct ParticleState {
  // Total rest mass [kg], signed charge [C], and positive momentum magnitude
  // [kg m s^-1].  Present scalar models depend on |chargeC|, while retaining
  // its sign keeps this state suitable for future directional-wave adapters.
  double massKg = 0.0;
  double chargeC = 0.0;
  double momentumKgMPerS = 0.0;

  // Required only by the energy-per-nucleon power-law variant.  The library
  // never infers mass number from charge or rest mass.
  std::optional<double> nucleonCount;
};

struct ParticleKinematics {
  // Centralized exact-relativistic values derived from ParticleState.  gamma
  // and beta are dimensionless; rigidity is pc/|q| [V], and kinetic energy is
  // the total-particle energy above rest energy [J].
  double gamma = 0.0;
  double beta = 0.0;
  double speedMPerS = 0.0;
  double rigidityV = 0.0;
  double kineticEnergyJ = 0.0;
};

Status ComputeParticleKinematics(const ParticleState& particle,
                                 ParticleKinematics* kinematics);

struct LocalState {
  // Coordinate time [s].  The current explicit models do not interpret this
  // number directly; a provider may use it to construct timeFactor.
  double timeS = 0.0;
  // Heliocentric Cartesian position [m].  The explicit power-law model uses
  // its Euclidean norm as r; adapters for translated meshes must subtract the
  // configured solar origin before constructing this value.
  std::array<double, 3> positionM{{0.0, 0.0, 0.0}};
  // The resolved mean field vector [T] defines B0 for empirical field factors
  // and mean-field Bohm.  effectiveFieldMagnitudeT [T] is a separately named
  // convention for the explicitly selected effective-field Bohm option; the
  // library never relabels it as the mean field.
  std::optional<std::array<double, 3> > meanFieldT;
  std::optional<double> effectiveFieldMagnitudeT;

  // These positive dimensionless multipliers are not physical defaults.  A
  // model that enables one requires the application/background adapter to
  // supply it.  Their construction, regional discontinuities, and time
  // interpretation (simultaneous, convected, or retarded) remain provider
  // responsibilities; the library only multiplies the declared value.
  std::optional<double> timeFactor;
  std::optional<double> regionFactor;
  std::optional<double> radialFactor;

  // Immutable-provider generation identities copied to provenance.  They
  // will also be required cache keys once caching is introduced at PD10.
  std::uint64_t backgroundRevision = 0;
  std::uint64_t turbulenceRevision = 0;
};

enum class IndependentVariable {
  Rigidity,
  TotalKineticEnergy,
  EnergyPerNucleon,
  Speed
};

enum class BohmFieldDefinition { MeanField, EffectiveField };

struct ConstantLambdaParameters {
  // Equation (11): lambda_parallel=lambdaParallelM [m].  Kappa still varies
  // with particle speed through kappa=v*lambda/3.
  double lambdaParallelM = 0.0;
};

struct ConstantKappaParameters {
  // Equation (12): kappa_parallel=kappaParallelM2PerS [m^2 s^-1].  The
  // corresponding mean free path is therefore species/speed dependent.
  double kappaParallelM2PerS = 0.0;
};

struct PowerLawLambdaParameters {
  // Equations (13)--(15).  lambda0M [m] is reached when every enabled ratio
  // and provider factor equals one.  independentReferenceSI has units set by
  // independentVariable: V, J, J/nucleon, or m/s, respectively.
  double lambda0M = 0.0;
  IndependentVariable independentVariable = IndependentVariable::Rigidity;
  double independentReferenceSI = 0.0;
  double independentExponent = 0.0;

  // If enabled, multiply by (|positionM|/radius0M)^radialExponent.  positionM
  // must already be relative to the solar origin, so this is heliocentric
  // radius rather than an arbitrary mesh-coordinate norm.
  bool useRadialFactor = false;
  double radius0M = 0.0;
  double radialExponent = 0.0;

  // If enabled with a nonzero exponent, multiply by
  // (|meanFieldT|/fieldReferenceT)^(-fieldExponent).  A zero exponent is the
  // exact unity identity and imposes no mean-field input requirement.
  bool useFieldFactor = false;
  double fieldReferenceT = 0.0;
  double fieldExponent = 0.0;

  bool useTimeFactor = false;
  bool useRegionFactor = false;
};

struct BrokenRigidityKappaParameters {
  // Equation (17) uses the speed-factored K_star convention.  This is not the
  // actual kappa at R0: at the reference state kappa=K_star*beta(R0).
  double kStarM2PerS = 0.0;
  double rigidity0V = 0.0;
  double breakRigidityV = 0.0;
  double lowSlope = 0.0;
  double highSlope = 0.0;
  double smoothness = 0.0;

  // Equation (17) field multiplier is
  // (fieldReferenceT/|meanFieldT|)^fieldExponent.  The radial and regional
  // forms are intentionally not invented here: when selected, their already
  // evaluated positive multipliers must be supplied in LocalState.
  bool useFieldFactor = false;
  double fieldReferenceT = 0.0;
  double fieldExponent = 0.0;
  bool useRadialFactor = false;
  bool useRegionFactor = false;
};

struct BohmParameters {
  // Equation (61): lambda_parallel=etaB*p/(|q|B).  etaB=1 is the conventional
  // comparison value but is not installed as a default and is not a bound on
  // any other backend.
  double etaB = 0.0;
  BohmFieldDefinition fieldDefinition = BohmFieldDefinition::MeanField;
};

struct ModelConfiguration {
  // Only the parameter member corresponding to model participates in
  // validation, evaluation, and the canonical fingerprint.  Zero member
  // initializers are invalid sentinels, not physical defaults.
  ModelId model = ModelId::ConstantLambda;
  ConstantLambdaParameters constantLambda;
  ConstantKappaParameters constantKappa;
  PowerLawLambdaParameters powerLawLambda;
  BrokenRigidityKappaParameters brokenRigidityKappa;
  BohmParameters bohm;
};

struct Provenance {
  // requested and evaluated IDs are separately retained so later explicit
  // fallback policies can report backend changes.  They are equal in PD02,
  // which has no fallback.
  std::string requestedModelId;
  std::string evaluatedModelId;
  std::string specificationVersion;
  std::string softwareVersion;
  // Deterministic configuration identity for restart/output comparison.  It
  // is not a cryptographic digest and is not a substitute for the SHA-256
  // identities required by future coefficient data sets.
  std::string configurationFingerprint;
  std::uint64_t backgroundRevision = 0;
  std::uint64_t turbulenceRevision = 0;
};

struct ParallelResult {
  // On success, both scalar optionals are present and obey Equation (4),
  // kappa_parallel=v*lambda_parallel/3.  On failure they remain absent; a
  // caller must inspect status before dereferencing either optional.
  Status status;
  std::optional<double> kappaParallelM2PerS;
  std::optional<double> lambdaParallelM;

  // These dimensionless logarithmic derivatives are analytical for the
  // explicit PD02 laws.  They mean partial ln(quantity)/partial ln(rigidity)
  // at fixed species and fixed LocalState, including the relativistic speed
  // contribution d ln(v)/d ln(rigidity)=1/gamma^2 where applicable.
  //
  // The Cartesian spatial gradient has units
  // [kappa]/[position]=m s^-1 and is defined at fixed particle momentum.
  // It is exactly zero only for declared constant models.  It remains absent
  // for other models until PD09 can chain coherent provider gradients;
  // absence must never be interpreted as a zero vector.
  std::optional<double> dLnLambdaDLnRigidity;
  std::optional<double> dLnKappaDLnRigidity;
  std::optional<std::array<double, 3> > gradKappaParallelMPerS;
  std::uint32_t diagnosticMask = DiagnosticNone;
  Provenance provenance;
};

Status ValidateConfiguration(const ModelConfiguration& configuration);
std::string ConfigurationFingerprint(const ModelConfiguration& configuration);

using ModelFunction = ParallelResult (*)(
    const ParticleState&, const LocalState&, const ModelConfiguration&);

// Public dispatch pointer requested by the application contract.  Input-file
// parsers must not assign it directly: call SetActiveConfiguration or
// ConfigureActiveModel so schema validation and the active parameters change
// transactionally.  Configure once during serial initialization, before any
// mover thread calls EvaluateActive.  It is a non-owning pointer to a
// static-lifetime evaluator; callers must not delete or replace it.
extern ModelFunction ActiveModelFunction;

ModelFunction FunctionForModel(ModelId model);
ParallelResult Evaluate(const ParticleState& particle,
                        const LocalState& local,
                        const ModelConfiguration& configuration);
ParallelResult EvaluateActive(const ParticleState& particle,
                              const LocalState& local);
Status SetActiveConfiguration(const ModelConfiguration& configuration);
ModelConfiguration GetActiveConfiguration();

// Parser-neutral input bridge.  Application parsers retain ownership of file
// syntax and source-location diagnostics, then pass the model identifier and
// exact key/value text here.  Keys name SI units; unknown and duplicate keys
// fail closed.  BuildConfiguration writes its output only after the complete
// candidate validates.  ConfigureActiveModel likewise leaves the previous
// active configuration and dispatch pointer unchanged after any failure.
struct InputParameter {
  std::string name;
  std::string value;
};

Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          ModelConfiguration* configuration);
Status ConfigureActiveModel(const std::string& modelId,
                            const std::vector<InputParameter>& parameters);

// Explicit converter for callers whose source gives the actual reference
// kappa instead of Equation (17)'s K_star.  It evaluates
// K_star=kappa_reference/beta(referenceParticle).  The reference particle is
// required because beta at a stated rigidity is species dependent.
Status KStarFromReferenceKappa(double referenceKappaM2PerS,
                               const ParticleState& referenceParticle,
                               double* kStarM2PerS);

}  // namespace ParallelDiffusion
}  // namespace SEP

#endif

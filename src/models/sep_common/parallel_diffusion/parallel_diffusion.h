#ifndef SEP_COMMON_PARALLEL_DIFFUSION_PARALLEL_DIFFUSION_H
#define SEP_COMMON_PARALLEL_DIFFUSION_PARALLEL_DIFFUSION_H

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace SEP {
namespace ParallelDiffusion {

// The scientific contract is PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md v1.4.
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
  // A selected quadrature failed its declared convergence controls.
  IntegrationFailed,
  // A coupled nonlinear closure failed its residual/iteration controls.
  NonlinearSolverFailed,
  // A spectrum provider violates its declared canonical normalization.
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

// Bit mask values intentionally match the normative revision-1.4 interface.
// A bit is set only when its evaluator establishes that condition; diagnostics
// must not be guessed from a finite scalar result.
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
  // Matrix convention: dMeanFieldDxTPerM[i][j] = partial B_i/partial x_j.
  // This is sufficient to differentiate |B| while the transport consumer
  // retains ownership of field-direction derivatives in Equation (62b).
  std::optional<std::array<std::array<double, 3>, 3> > dMeanFieldDxTPerM;
  std::optional<std::array<double, 3> >
      gradEffectiveFieldMagnitudeTPerM;

  // These positive dimensionless multipliers are not physical defaults.  A
  // model that enables one requires the application/background adapter to
  // supply it.  Their construction, regional discontinuities, and time
  // interpretation (simultaneous, convected, or retarded) remain provider
  // responsibilities; the library only multiplies the declared value.
  std::optional<double> timeFactor;
  std::optional<double> regionFactor;
  std::optional<double> radialFactor;
  std::optional<std::array<double, 3> > gradTimeFactorPerM;
  std::optional<std::array<double, 3> > gradRegionFactorPerM;
  std::optional<std::array<double, 3> > gradRadialFactorPerM;

  // Canonical turbulence quantities used by the spectral and nonlinear
  // closures.  They follow Equations (25), (27), and (46): both variances are
  // total two-component magnetic variances [T^2], and both lengths are
  // spectral bend-over lengths [m], not silently substituted integral
  // correlation lengths.  A model validates only the members it consumes.
  struct TurbulenceState {
    std::optional<double> densityKgPerM3;
    std::optional<double> slabVarianceT2;
    std::optional<double> twoDVarianceT2;
    std::optional<double> slabBendoverLengthM;
    std::optional<double> twoDBendoverLengthM;
    std::optional<double> inertialIndex;
    // Raw moment fields are consumed only by turbulence_adapter after its
    // named convention has been configured. They are not canonical magnetic
    // variances and therefore remain separately named.
    std::optional<double> providerMomentM2PerS2;
    std::optional<double> residualEnergy;
    std::optional<double> slabFraction;

    // Optional Cartesian gradients use component j = partial/partial x_j.
    // They are absent unless the provider supplies one coherent snapshot;
    // the library never interprets absence as a spatially uniform field.
    std::optional<std::array<double, 3> > gradSlabVarianceT2PerM;
    std::optional<std::array<double, 3> > gradTwoDVarianceT2PerM;
    std::optional<std::array<double, 3> > gradSlabLengthMPerM;
    std::optional<std::array<double, 3> > gradTwoDLengthMPerM;
  };
  std::optional<TurbulenceState> turbulence;

  // Some closures deliberately consume an independently calculated
  // perpendicular coefficient.  Its identity and revision are required so
  // nlpa_given_perp cannot be mistaken for a self-contained coupled model.
  std::optional<double> suppliedKappaPerpendicularM2PerS;
  std::string suppliedPerpendicularModelId;
  std::uint64_t suppliedPerpendicularRevision = 0;

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

enum class SpectrumForm { SmoothBendover, Multirange, SuppliedLogLog };
enum class SpectrumTailPolicy { OutOfDomain, PowerLaw, Zero };
enum class PitchAngleAmplitudeMode { TargetLambda, FixedAmplitude };
enum class BroadeningKernel { LorentzianConstant, LorentzianLinear, Gaussian };
enum class AdapterClosure { QltSlabSpectrum, QltSlabInertial, BroadenedSlab };
enum class StoredCoefficient { LambdaParallel, KappaParallel };
enum class TimeInterpolation { Linear, StepPrevious };
enum class TableAxis {
  Rigidity,
  TotalKineticEnergy,
  EnergyPerNucleon,
  Speed,
  HeliocentricRadius,
  Time,
  MeanFieldMagnitude
};
enum class MomentEnergyConvention {
  KineticPlusMagneticVariance,
  HalfElsasserSum,
  ElsasserSum,
  SpecificTotalFluctuationEnergy
};
enum class ResidualEnergyConvention {
  KineticMinusMagnetic,
  MagneticMinusKinetic
};
enum class MomentSource { InputParameters, LocalState };

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

// Numerical controls are algorithmic accuracy requirements, not physical
// calibration values.  The defaults are the initial float64 controls stated
// in Sections 10.5 and 14.4; callers may request tighter values but cannot
// relax the nonlinear residual beyond the specification's 1e-8 gate.
struct NumericalParameters {
  double relativeTolerance = 1.0e-8;
  int maximumRefinements = 18;
  int maximumIterations = 80;
};

struct SpectrumParameters {
  SpectrumForm form = SpectrumForm::SmoothBendover;
  // Equation (35) multirange parameters.  The runtime turbulence snapshot
  // supplies delta B_s^2, ell_s, and the inertial index s.  These three
  // configuration values supply q_E, s_d, and k_d [rad m^-1].  Selection is
  // valid only when q_E>-1, s_d>1, and the runtime product k_d*ell_s>1.
  double energyRangeIndex = 0.0;
  double dissipationIndex = 0.0;
  double dissipationWavenumberRadPerM = 0.0;

  // A supplied spectrum is canonical Equation (25) data: strictly increasing
  // k [rad m^-1] and nonnegative P_s [T^2 m].  Log-log interpolation is used
  // only between positive ordinates.  No extrapolation or artificial tail is
  // installed; a resonance outside coverage is MissingInput.
  std::vector<double> wavenumberRadPerM;
  std::vector<double> powerT2M;
  double declaredVarianceT2 = 0.0;
  std::string sourceIdentity;
  SpectrumTailPolicy lowKPolicy = SpectrumTailPolicy::OutOfDomain;
  SpectrumTailPolicy highKPolicy = SpectrumTailPolicy::OutOfDomain;
  double lowKPowerIndex = 0.0;
  double highKPowerIndex = 0.0;
};

struct QltSlabParameters {
  SpectrumParameters spectrum;
  NumericalParameters numerical;
};

struct QltInertialParameters {
  // The state supplies B0, slab variance, bend-over length, and s.  This
  // explicit selector keeps Equation (31) distinct from full Equation (30).
};

struct PrescribedLambdaMuParameters {
  PitchAngleAmplitudeMode amplitudeMode =
      PitchAngleAmplitudeMode::TargetLambda;
  double qMu = 0.0;
  double hMu = 0.0;
  double targetLambdaM = 0.0;
  double fixedD0PerS = 0.0;
  NumericalParameters numerical;
};

struct BroadenedSlabParameters {
  SpectrumParameters spectrum;
  BroadeningKernel kernel = BroadeningKernel::LorentzianConstant;
  // Gamma0 or Delta0 [s^-1]; uDec [m s^-1] appears only in Equation (41).
  // A zero width explicitly selects the exact QLT branch, never a narrow
  // fixed-grid approximation.
  double width0PerS = 0.0;
  double decorrelationSpeedMPerS = 0.0;
  NumericalParameters numerical;
};

struct NonlinearParameters {
  NumericalParameters numerical;
};

struct NlgceFParameters {
  // The published arrays and natural-log convention are fixed by the model
  // identity.  This field exists to make the required coefficient identity
  // explicit in fingerprints and restart metadata.
  std::string coefficientSet = "Qin_Zhang_2014_Tables_3_4";
};

struct TurbulenceAdapterParameters {
  MomentEnergyConvention energyConvention =
      MomentEnergyConvention::KineticPlusMagneticVariance;
  ResidualEnergyConvention residualConvention =
      ResidualEnergyConvention::KineticMinusMagnetic;
  MomentSource momentSource = MomentSource::InputParameters;
  AdapterClosure closure = AdapterClosure::QltSlabSpectrum;
  // Vacuum permeability [H m^-1] is explicit because the post-2019 SI value
  // is measured rather than exact and revision 1.4 does not freeze a CODATA
  // release for this adapter conversion.
  double vacuumPermeabilityHPerM = 0.0;
  std::optional<double> providerMomentM2PerS2;
  std::optional<double> residualEnergy;
  std::optional<double> slabFraction;
  QltSlabParameters qlt;
  BroadenedSlabParameters broadened;
};

struct WaveSpectrumAdapterParameters {
  // Revision 1.4 only permits the balanced, magnetostatic reduction.  These
  // named declarations are required rather than inferred from a handle.
  std::string propagation = "";
  std::string polarization = "";
  std::string frame = "";
  AdapterClosure closure = AdapterClosure::QltSlabSpectrum;
  QltSlabParameters qlt;
  BroadenedSlabParameters broadened;
};

struct TabulatedParallelParameters {
  StoredCoefficient storedCoefficient = StoredCoefficient::LambdaParallel;
  // Axes and flattened values use row-major order with the final axis varying
  // fastest. Positive axes are interpolated in log space; time alone is
  // linear. Every boundary is out-of-domain, as required by the undecided
  // D11 default—there is no extrapolation selector in this release.
  std::vector<TableAxis> axes;
  std::vector<std::vector<double> > axisSI;
  std::vector<double> coefficientSI;
  std::string generationIdentity;
  std::optional<TimeInterpolation> timeInterpolation;
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
  QltSlabParameters qltSlab;
  QltInertialParameters qltInertial;
  PrescribedLambdaMuParameters prescribedLambdaMu;
  BroadenedSlabParameters broadenedSlab;
  NonlinearParameters nonlinear;
  NlgceFParameters nlgceF;
  TurbulenceAdapterParameters turbulenceAdapter;
  WaveSpectrumAdapterParameters waveSpectrumAdapter;
  TabulatedParallelParameters table;
};

struct Provenance {
  // requested and evaluated IDs are separately retained so later explicit
  // fallback/adapter policies can report backend changes. They remain equal
  // whenever the requested backend itself produced the coefficient.
  std::string requestedModelId;
  std::string evaluatedModelId;
  std::string specificationVersion;
  std::string softwareVersion;
  // Populated by data-backed/adapted models.  The NLGCE-F string contains the
  // two audited CSV SHA-256 values; adapter strings preserve the provider
  // convention rather than only the converted canonical moment.
  std::string coefficientSetIdentity;
  std::string inputMomentConvention;
  // SHA-256 of the model's canonical active-schema serialization for
  // restart/output comparison. External coefficient/data identities remain
  // separate because their bytes are not replaced by this configuration hash.
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

  struct PerpendicularPair {
    double kappaPerpendicularM2PerS = 0.0;
    double lambdaPerpendicularM = 0.0;
  };
  std::optional<PerpendicularPair> perpendicular;

  // D_mu_mu is evaluated at a caller-selected pitch angle through the
  // separate EvaluatePitchAngleDiffusion API below.  The coefficient result
  // records nonlinear convergence information when applicable.
  std::optional<int> nonlinearIterations;
  std::optional<double> nonlinearMaxLogResidual;

  // These dimensionless logarithmic derivatives are analytical for the
  // supported explicit/spectral laws. They mean partial ln(quantity)/partial ln(rigidity)
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

// Evaluate the pitch-angle coefficient [s^-1] for models that explicitly
// define D_mu_mu.  Unsupported eigenvalue-only models fail with
// UnsupportedModel; endpoints use their mathematical limits and are never
// regularized by an implicit mu cutoff or scattering floor.
Status EvaluatePitchAngleDiffusion(double mu,
                                   const ParticleState& particle,
                                   const LocalState& local,
                                   const ModelConfiguration& configuration,
                                   double* dMuMuPerS);

// Equal-length batch evaluation required by Section 14.1.  Shape validation
// occurs before writing results.  Individual physics failures remain local to
// their element; an empty batch is valid.
Status EvaluateBatch(const std::vector<ParticleState>& particles,
                     const std::vector<LocalState>& states,
                     const ModelConfiguration& configuration,
                     std::vector<ParallelResult>* results);

// Named Section 8.7 conversions into the canonical Equation (25) spectrum:
// k>0 in radians per metre, one-sided, total transverse slab power P_s(k) in
// T^2 m.  These scalar adapters are intentionally explicit; a provider loops
// over its immutable samples and records the named conversion in provenance.
// Nonnegative powers are allowed because a genuine spectral zero must remain
// distinguishable from missing coverage.
Status ConvertOneSidedComponentsToCanonical(
    double pXxT2M, double pYyT2M, double* canonicalPowerT2M);
Status ConvertTwoSidedComponentsToCanonical(
    double sXxPositiveT2M, double sXxNegativeT2M,
    double sYyPositiveT2M, double sYyNegativeT2M,
    double* canonicalPowerT2M);
Status ConvertEvenTwoSidedTotalToCanonical(
    double signedTotalPowerT2M, double* canonicalPowerT2M);
Status ConvertOneSidedCyclesPerMToCanonical(
    double wavenumberCyclesPerM, double powerT2PerCyclePerM,
    double* wavenumberRadPerM, double* canonicalPowerT2M);

struct FrozenFlowMapping {
  // Positive magnitude of the sampling-velocity projection [m s^-1] along
  // the direction used to identify k.  The nonempty identity records the
  // caller's frozen-flow and geometry assumptions; the library cannot infer
  // them from a temporal PSD.
  double samplingVelocityProjectionMPerS = 0.0;
  std::string assumptionIdentity;
};

Status ConvertFrozenFlowFrequencyToCanonical(
    double frequencyHz, double oneSidedFrequencyPsdT2PerHz,
    const FrozenFlowMapping& mapping, double* wavenumberRadPerM,
    double* canonicalPowerT2M);
Status ConvertQinZhangSlabComponentToCanonical(
    double sourcePowerT2M, double* canonicalPowerT2M);
Status ConvertQinZhangTwoDComponentToCanonical(
    double sourcePowerT2M, double* canonicalPowerT2M);

// Equation (57b) conversion used by turbulence_adapter.  The provider's
// convention and residual-energy sign are explicit inputs, eliminating the
// otherwise ambiguous factors of two and sign.
Status ConvertTurbulenceMoment(double providerMomentM2PerS2,
                               double densityKgPerM3,
                               double residualEnergy,
                               double vacuumPermeabilityHPerM,
                               MomentEnergyConvention energyConvention,
                               ResidualEnergyConvention residualConvention,
                               double* totalMagneticVarianceT2);

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

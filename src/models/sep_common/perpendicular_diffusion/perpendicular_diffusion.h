#ifndef SEP_COMMON_PERPENDICULAR_DIFFUSION_PERPENDICULAR_DIFFUSION_H
#define SEP_COMMON_PERPENDICULAR_DIFFUSION_PERPENDICULAR_DIFFUSION_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace SEP {
namespace PerpendicularDiffusion {

// Scientific contract: PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md v2.1.
// This dependency-free C++17 interface returns coefficient eigenvalues and
// explicitly typed diagnostics.  It does not own a mesh, stochastic mover,
// field-line realization, MPI exchange, turbulence evolution, or injection.
// All dimensional inputs and outputs are SI after one host-adapter conversion.
constexpr double SpeedOfLightMPerS = 299792458.0;

enum class StatusCode {
  Success,
  InvalidInput,
  MissingInput,
  MissingCalibration,
  OutsideModelDomain,
  IncompatibleGeometry,
  InconsistentPair,
  DivergentMoment,
  IntegrationFailed,
  NonlinearSolverFailed,
  DerivativeUnavailable,
  SourceGate,
  InvalidConfiguration,
  UnsupportedModel
};

struct Status {
  StatusCode code = StatusCode::Success;
  std::string detail;
  bool ok() const { return code == StatusCode::Success; }
  static Status Success();
  static Status Error(StatusCode code, const std::string& detail);
};

enum class Observable {
  SymmetricCoefficient,
  PairedCoefficients,
  PitchAngleCoefficient,
  SignedHallCoefficient,
  ParticleMsd,
  FieldLineCoefficient,
  ConditionalReturnMedian,
  PerturbativeDiagnostic
};
enum class Estimator {
  Asymptotic,
  InstantaneousDerivative,
  Secant,
  FiniteWindowAverage,
  MomentValue,
  NotApplicable
};
enum class FrameKind { AxisymmetricOrdered, OrientedUnequal, Isotropic };
enum class DomainState { InsideDeclaredDomain, Extrapolated, NotApplicable };
enum class Quality {
  SourceFit,
  ParameterizedClosure,
  ScalingEstimate,
  Diagnostic
};
enum class DiffusionRegime {
  NormalDiffusion,
  NoNormalDiffusion,
  NotEstablished,
  NotApplicable
};
enum class DependencyOwner { SuppliedParallel, PairedBackend, None };
enum class GeometryKind { CompositeSlab2D, Pure2D, PureSlab, Isotropic3D };

enum class ModelId {
  ConstantKappaPerp,
  ConstantLambdaPerp,
  RatioKappa,
  PowerLawPerp,
  PitchAnglePerp,
  DrogeLambdaPerp,
  ParadiseAlpha,
  FieldLineSlab,
  FieldLine2D,
  FieldLineComposite,
  FlrwParticle,
  Nlgc,
  Enlgc2D,
  NlgceN,
  NlgceF2014,
  Unlt,
  ImplicitSlabExact2016,
  ImplicitSlabRational2016,
  FlpdComplete,
  RbdBc,
  CompositeClosed2019,
  CompoundDiffusiveLines,
  GcdCompound,
  FrozenFieldline,
  PrediffusiveFit,
  IsoFitCandiaRoulet2004,
  IsoFitSnodin2016,
  IsoKappa0Snodin2016,
  IsoFitKuhlen2025,
  ClassicalScattering,
  YanLazarianMa4,
  NwuRatioPolar,
  CortiAms02,
  HelmodRatio,
  DriftWeakScattering,
  DriftRigidityReduction,
  DriftClassical,
  DriftCandiaRoulet2004,
  TabulatedPerp,
  QltPerp,
  UnltFgrPerturbativeDiagnostic,
  UnltFgr,
  EnlgcSlabSecantDiagnostic,
  NlgcSlabKernelDiagnostic,
  IsoFitCasse2002,
  RestrictedScattering
};

struct ModelDescriptor {
  ModelId id;
  const char* stableId;
  bool executable;
  bool productionMoverEligible;
  Observable observable;
  const char* equation;
};

const std::vector<ModelDescriptor>& ModelRegistry();
const char* ModelName(ModelId model);
bool ParseModelId(const std::string& text, ModelId* model);

struct ParticleState {
  // Total rest mass [kg], signed charge [C], total momentum magnitude
  // [kg m s^-1], and optional pitch-angle cosine.  Charge sign is retained
  // for Hall coefficients; symmetric coefficients use |charge|.
  double massKg = 0.0;
  double chargeC = 0.0;
  double momentumKgMPerS = 0.0;
  std::optional<double> mu;
};

struct ParticleKinematics {
  double gamma = 0.0;
  double beta = 0.0;
  double speedMPerS = 0.0;
  double rigidityV = 0.0;
};

Status ComputeParticleKinematics(const ParticleState& particle,
                                 ParticleKinematics* kinematics);

struct ParallelInput {
  // A supplied dependency is a coefficient, not merely a number: its model,
  // equation version, and coherent sample identity are mandatory provenance.
  double kappaM2PerS = 0.0;
  std::string modelId;
  std::string equationVersion;
  std::string sampleFingerprint;
};

struct TurbulenceSample {
  GeometryKind geometry = GeometryKind::CompositeSlab2D;
  std::string geometryId;
  std::string spectrumId;
  std::string energyConvention;
  std::string providerRevision;
  std::string sampleFingerprint;

  // Composite quantities follow (S2): variances are total two-component
  // transverse magnetic variances [T^2], and ell_s/ell_2 are spectral
  // bend-over scales [m].  Isotropic totalVarianceT2 instead includes all
  // three random-field components and must not be repartitioned implicitly.
  std::optional<double> slabVarianceT2;
  std::optional<double> twoDVarianceT2;
  std::optional<double> totalVarianceT2;
  std::optional<double> slabBendoverLengthM;
  std::optional<double> twoDBendoverLengthM;
  std::optional<double> inertialIndex;
  std::optional<double> energyRangeIndex;

  // Optional source-specific scales. Their meanings stay explicit; no
  // bend-over/integral/outer-scale conversion occurs without a named model.
  std::optional<double> outerScaleM;
  std::optional<double> correlationLengthM;
  std::optional<double> transverseCorrelationLengthM;
};

struct LocalState {
  std::array<double, 3> positionM{{0.0, 0.0, 0.0}};
  double timeS = 0.0;
  std::optional<std::array<double, 3> > meanFieldT;
  TurbulenceSample turbulence;
  std::optional<ParallelInput> parallelDependency;

  // Geometry-dependent inputs are supplied rather than inferred.  The two
  // perpendicular axes are required when a model returns unequal eigenvalues.
  std::optional<double> heliocentricColatitudeRad;
  std::optional<double> parkerSpiralCosine;
  std::optional<std::array<double, 3> > perpendicularAxis1;
  std::optional<std::array<double, 3> > perpendicularAxis2;

  // A pre-evaluated field-line closure can be consumed by FLRW, FLPD, and
  // composite models. Its identity prevents an unlabeled length from crossing
  // the API boundary.
  std::optional<double> fieldLineCoefficientM;
  std::string fieldLineModelId;
  std::string fieldLineSampleFingerprint;
};

struct InputParameter {
  std::string name;
  std::string value;
};

struct ModelConfiguration {
  ModelId model = ModelId::ConstantKappaPerp;
  // Model-specific readers populate only their declared schema.  SI-valued
  // numbers and enumerated text remain separated so canonical fingerprints do
  // not depend on parser order.  Callers should use BuildConfiguration rather
  // than populate these maps directly.
  std::map<std::string, double> numeric;
  std::map<std::string, std::string> text;
  std::vector<double> tableAxisSI;
  std::vector<double> tableValuesSI;
};

struct PhysicalFrame {
  FrameKind kind = FrameKind::AxisymmetricOrdered;
  std::optional<std::array<double, 3> > b;
  std::optional<std::array<double, 3> > e1;
  std::optional<std::array<double, 3> > e2;
  std::string perpendicularAxisDefinition;
};

struct SpatialCoefficients {
  PhysicalFrame frame;
  std::optional<double> parallelM2PerS;
  std::optional<double> perpendicular1M2PerS;
  std::optional<double> perpendicular2M2PerS;
  std::optional<double> isotropicM2PerS;
  DependencyOwner parallelOwner = DependencyOwner::None;
};

struct ParticleMoments {
  double ageS = 0.0;
  std::array<double, 3> rawMsdM2{{0.0, 0.0, 0.0}};
  std::optional<std::array<double, 3> > derivativeM2PerS;
  std::optional<std::array<double, 3> > secantM2PerS;
};

struct NamedSensitivity {
  std::string quantity;
  std::string independentInput;
  std::string unit;
  double value = 0.0;
  std::optional<double> estimatedAbsoluteError;
};

struct Provenance {
  std::string requestedModel;
  std::string usedModel;
  std::string equationVersion;
  std::string sourceVersion;
  std::string configurationFingerprint;
  std::string sampleFingerprint;
  std::optional<std::string> calibrationId;
  std::optional<std::string> tableChecksum;
  std::map<std::string, std::string> conventions;
  DependencyOwner dependencyOwner = DependencyOwner::None;
  std::optional<std::string> parallelModelId;
};

struct NumericalDiagnostics {
  std::string backend;
  std::string integrationMethod;
  std::string rootMethod;
  std::optional<double> absoluteErrorEstimate;
  std::optional<double> relativeErrorEstimate;
  std::optional<double> residualNorm;
  std::vector<std::string> approximations;
  std::vector<std::string> unavailableDerivatives;
  std::optional<std::size_t> iterations;
  std::map<std::string, double> controls;
};

struct ModelResult {
  Status status;
  Observable observable = Observable::SymmetricCoefficient;
  Estimator estimator = Estimator::Asymptotic;
  DomainState domain = DomainState::NotApplicable;
  Quality quality = Quality::ParameterizedClosure;
  DiffusionRegime diffusionRegime = DiffusionRegime::NotApplicable;
  std::optional<SpatialCoefficients> coefficients;
  std::optional<double> pitchAngleM2PerS;
  std::optional<double> fieldLineM;
  std::optional<ParticleMoments> particleMoments;
  std::optional<double> conditionalMedianM2;
  std::optional<double> signedHallM2PerS;
  std::vector<NamedSensitivity> sensitivities;
  Provenance provenance;
  NumericalDiagnostics numerical;
};

Status ValidateConfiguration(const ModelConfiguration& configuration);
Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          ModelConfiguration* configuration);
std::string ConfigurationFingerprint(const ModelConfiguration& configuration);

using ModelFunction = ModelResult (*)(const ParticleState&, const LocalState&,
                                      const ModelConfiguration&);
// ConfigureActiveModel installs a distinct static-lifetime function for the
// selected registry identity together with its validated immutable parameter
// value. Configure only during serial initialization; changing this pair while
// worker threads evaluate it is outside the API contract.
extern ModelFunction ActiveModelFunction;
ModelFunction FunctionForModel(ModelId model);
ModelResult Evaluate(const ParticleState& particle, const LocalState& local,
                     const ModelConfiguration& configuration);
ModelResult EvaluateActive(const ParticleState& particle,
                           const LocalState& local);
Status SetActiveConfiguration(const ModelConfiguration& configuration);
Status ConfigureActiveModel(const std::string& modelId,
                            const std::vector<InputParameter>& parameters);
ModelConfiguration GetActiveConfiguration();

Status EvaluateBatch(const std::vector<ParticleState>& particles,
                     const std::vector<LocalState>& states,
                     const ModelConfiguration& configuration,
                     std::vector<ModelResult>* results);

// Construct the symmetric Cartesian tensor from returned eigenvalues.  The
// matrix uses physical Cartesian components [m^2 s^-1], not spherical
// coordinate components.  Unequal transverse eigenvalues require a complete
// oriented frame; isotropic results require no axis.
Status AssembleSymmetricTensor(
    const SpatialCoefficients& coefficients,
    std::array<std::array<double, 3>, 3>* tensorM2PerS);

// Public equation helpers allow focused tests without exposing numerical
// implementation details. They retain the exact source normalizations.
double SpectrumC(double inertialIndex);
double SpectrumD(double inertialIndex, double energyRangeIndex);
Status SmoothTwoDSpectralMoment(double inertialIndex,
                                double energyRangeIndex,
                                int moment,
                                double varianceT2,
                                double bendoverLengthM,
                                double* value);
Status SmoothSpectrumLengths(double inertialIndex,
                             double energyRangeIndex,
                             double bendoverLengthM,
                             double* ultraScaleM,
                             double* perpendicularIntegralScaleM);
Status EvaluateImplicitSlabKernel(double xi, bool rational,
                                  double* kernel);

}  // namespace PerpendicularDiffusion
}  // namespace SEP

#endif

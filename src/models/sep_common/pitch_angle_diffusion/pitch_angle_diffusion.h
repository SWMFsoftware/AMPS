#ifndef SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_DIFFUSION_H
#define SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_DIFFUSION_H

// ============================================================================
// Pitch-angle diffusion coefficient library (Release A).
//
// Scientific contract: PITCH_ANGLE_SCATTERING_MODEL.md, specification v1.2
// (13 October 2026), with the implementation supplement in
// pitch_angle_diffusion_model_data/ (Sections 16.12 and 17, roadmap Section
// 26).  This header is the coefficient layer of the six-layer design of
// Section 17: model registry and provenance, validated SI coefficients, the
// lambda_parallel normalization, scattering step limits and the coefficient-
// kernel readiness matrix.  The stochastic kernels and the independent
// finite-volume reference operator are in pitch_angle_kernels.h.
//
// What the library realizes (Eq. (N1), Section 16.12).  In the local frame in
// which scattering is elastic (here: the local plasma/solar-wind frame in
// which the host stores p and mu), the pitch-angle operator is
//
//     (df/dt)_scat = d/dmu ( D_mumu(mu) df/dmu ),      -1 <= mu <= 1,     (21)
//
// with D_mumu >= 0 [s^-1] the coefficient of the *standard* operator.  Every
// value returned here is that standard coefficient: a source written with
// another operator factor (EPREM: d/dmu (D^E/2 df/dmu)) is converted exactly
// once, inside the model (Section 3.3).  The forward Ito SDE realizing (21)
// together with the host's deterministic pitch-angle rate is
//
//     dmu = (dmu/dt + dD_mumu/dmu) dt + sqrt(2 D_mumu) dW,                (25)
//
// so a kernel needs D and its derivative; both are returned together.  The
// library never adds streaming, focusing, flow or adiabatic terms: those are
// owned by the host's transport composer (Section 16.12, "ownership").
//
// Release A scope (user decision 2026-10-13; roadmap Section 26.1):
//   PA-ISO    D = D0 (1 - mu^2)                               (24), (5.1)
//   PA-QFORM  D = D0 (1 - mu^2)(|mu|^(q-1) + H), 1 < q < 2, H >= 0     (26)
//   PA-EPS    D = (nu0/2)(|mu|/(1+|mu|) + eps)(1 - mu^2),   eps > 0   (29)
//   PA-EPREM  printed D^E = (1-mu^2) v/(2 lambda) with operator D^E/2:
//             standard D = (1-mu^2) v/(4 lambda), lambda(Eq. 23) = 2 lambda
//             (34), (5.4), MFP v1.2 D-2.
// The amplitude of every shape is normalized per particle to a parallel mean
// free path given by the configured law (user decision 2026-10-13)
//
//     lambda_law = lambda_ref (R/R_ref)^a (r/r_ref)^b,    R = p c/|q| [V],
//
// with r the heliocentric radius [m], through the exact relation (23)
//
//     lambda_par = (3 v / 8) Int_{-1}^{1} (1-mu^2)^2 / D_mumu dmu,        (23)
//
// i.e. D0 = N_shape v / lambda_law with the shape constant N_shape of
// Section 5.6 (1/2 for PA-ISO, (3/4) I(q,H) for PA-QFORM, phi(eps)/2 for
// PA-EPS as D0 = nu0/2, and 1/4 for PA-EPREM by its printed convention).
// Every other inventory ID of Section 4.2 is registered (Registry()) and is
// refused with a typed status naming its section and runtime state; no other
// model is ever substituted.
//
// Units: everything crossing this boundary is SI.  Momentum is the total
// relativistic momentum magnitude [kg m s^-1] in the frame in which mu is
// defined, charge is signed [C], rigidity is positive [V], lengths [m],
// D_mumu and rates [s^-1], times [s].  Input-file adapters convert AU, GV or
// MeV exactly once before entering this API; the library never does.
//
// Frame and clock contract (user decision 2026-10-13).  p and mu are taken in
// the local plasma frame and D_mumu is applied per unit host time.  The exact
// clock factor of Eq. (173), dt_plasma/dt_host = gamma_u (1 - u v mu/c^2),
// differs from one by O(U/c) (about 1e-3 for U = 400 km/s, v ~ c) and is not
// applied: the coefficient is declared as a mixed-frame coefficient retained
// to order 0 in U/c (FrameContract below).  A host that needs the exact clock
// must multiply the whole generator per Eq. (N7), never rescale D and then
// differentiate it.
//
// Thread safety: configuration (Set/Configure*) happens once during serial
// initialization; evaluation functions are pure and reentrant afterwards.
// ============================================================================

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace SEP {
namespace PitchAngleDiffusion {

// Exact SI constants (2019 SI) used by the kinematics.  No other constant is
// needed by Release A; masses are supplied by the host per species.
constexpr double SpeedOfLightMPerS = 299792458.0;

// Identity of the scientific contract, copied to provenance.
constexpr const char* SpecificationVersion =
    "PITCH_ANGLE_SCATTERING_MODEL.md v1.2 (2026-10-13)";
constexpr const char* SoftwareVersion = "pitch_angle_diffusion Release A";

// ---------------------------------------------------------------------------
// Typed status (Section 17, Result<T>/Diagnostic; roadmap Section 26.7).
// Unavailable data, physical refusal and nonconverged numerics have distinct
// codes; a failed call never returns a coefficient (no NaN-as-success).
enum class StatusCode {
  Success,
  // Invalid particle state: mass, charge, momentum or mu outside its domain.
  InvalidInput,
  // Invalid local background (for example a non-positive heliocentric radius
  // where the lambda law needs r).
  InvalidField,
  // Finite, well-formed inputs outside a model's stated domain.
  OutsideDomain,
  // A source formula is not established for execution (open audit item).
  UnresolvedSource,
  // A printed form that is not a diffusion coefficient as printed.
  Inadmissible,
  // The coefficient cannot be paired with the requested kernel (for example
  // an unbounded drift with a derivative-based kernel, Section 16.11).
  IncompatibleKernel,
  // A quadrature or closed-form cross-check failed its declared tolerance.
  Nonconverged,
  // A step rule would return a non-positive interval.
  Nonprogress,
  // A configuration value is malformed, duplicated, unknown or inactive.
  InvalidConfiguration,
  // A required configuration value is absent.
  MissingInput,
  // A registered inventory ID whose backend is outside Release A.
  UnsupportedModel
};

struct Status {
  StatusCode code = StatusCode::Success;
  std::string detail;
  bool ok() const { return code == StatusCode::Success; }
  static Status Success() { return Status(); }
  static Status Error(StatusCode code, const std::string& detail);
};

const char* Name(StatusCode code);

// ---------------------------------------------------------------------------
// Model registry (Section 4).  The rows are generated from the specification
// inventory tables of Section 4.2 by generate_registry.py; the canonical
// companion archive (with admissibility/capability/domain per ID) is not
// available, so only the fields the specification prints are compiled.
enum class RuntimeState {
  ReadyExplicitInputs,
  ReadyPublishedVariant,
  RequiresExternalInput,
  RequiresUserDecision,
  RequiresSourceOrCodeAudit,
  ReferenceDataOnly,
  // Assigned by an implementation when a companion text is unavailable
  // (Section 4.1); no Release A model is in this state.
  BlockedDependency
};
const char* Name(RuntimeState state);

struct RegistryEntry {
  const char* id;         // stable ID, e.g. "PA-QFORM"
  const char* aliases;    // comma-separated alias IDs (Section 4.2), or ""
  const char* section;    // defining subsection, e.g. "5.2"
  const char* equations;  // equation labels of this document
  RuntimeState state;     // registry runtime state (Section 4.1)
};

// All 172 inventory IDs in specification order.
const std::vector<RegistryEntry>& Registry();
// Case-insensitive lookup by ID or alias; null when the ID is not in the
// inventory (then it is not implementable under this specification).
const RegistryEntry* FindRegistryEntry(const std::string& id);
// SHA-256 of the specification file the registry was generated from.
const char* SpecificationSha256();

// The Release A backends.  Stable input spellings are the inventory IDs.
enum class ModelId { Isotropic, QForm, EpsilonForm, Eprem };
const char* ModelName(ModelId model);  // "PA-ISO", "PA-QFORM", ...
// Typed parse: Success for an implemented ID; UnsupportedModel (naming the
// section and runtime state) for a registered ID outside Release A;
// InvalidConfiguration for an ID absent from the inventory.
Status ParseModelId(const std::string& text, ModelId* model);

// Operator convention of the source (Section 3.3).  The returned coefficient
// is always the Standard one; this records which conversion was applied.
enum class OperatorConvention { Standard, EpremHalfD };
const char* Name(OperatorConvention convention);

// Behaviour of the coefficient at mu = 0 (Sections 16.5, 16.11).
//   Regular          D smooth through 0 (PA-ISO, PA-EPREM);
//   DriftJump        D' bounded but discontinuous: -D0 -> +D0 (PA-EPS, Eq. 164);
//   UnboundedDrift   D' ~ |mu|^(q-2) diverges for 1<q<2 and every H (Eq. 163);
//   InternalBarrier  Int dmu/D diverges across 0: the hemispheres do not
//                    communicate by scattering (PA-EPS with eps = 0).
enum class ZeroClass { Regular, DriftJump, UnboundedDrift, InternalBarrier };
const char* Name(ZeroClass zero);

// Pusher kernels known to the readiness matrix (Section 16.11).  The library
// implements all four in pitch_angle_kernels.h; the finite-volume chain is a
// reference operator for verification, not a production pusher kernel.
enum class Kernel { EulerMaruyama, Milstein, RandomRotation, FiniteVolumeChain };
const char* Name(Kernel kernel);

// Scattering step rules (Section 16.6; user decisions 2026-10-13):
//   RateFraction    dt = f / nu_max,
//   VanDenBerg2026  dt = max( min{ ell sqrt(1-mu^2)/|A|, ell^2 (1-mu^2)/(2D) },
//                             f / nu_max )       (Eq. 165, floored everywhere),
//   Sandroos2014    dt = min( dmu_max / a_max, dmu_max^2 / b_max^2 ) (Eq. 166),
// where nu_max = sup_mu 2 D/(1-mu^2) is the largest local scattering
// frequency, A = dD/dmu + dmu/dt the forward Ito drift (the source prints the
// backward drift dD/dmu - dmu/dt), a_max = sup_mu |dD/dmu| and
// b_max^2 = sup_mu 2 D, all at the particle's momentum and position.
enum class StepRule { RateFraction, VanDenBerg2026, Sandroos2014 };
const char* Name(StepRule rule);

// ---------------------------------------------------------------------------
// Particle and local state.
struct ParticleState {
  double massKg = 0.0;           // rest mass of the particle [kg]
  double chargeC = 0.0;          // signed charge [C]; nonzero
  double momentumKgMPerS = 0.0;  // total momentum magnitude [kg m s^-1], > 0
};

struct Kinematics {
  // Exact relativistic quantities derived once from ParticleState
  // (Section 1.6): gamma = sqrt(1 + (p/mc)^2), v = p c^2/E, R = p c/|q|,
  // T = p^2 c^2/(E + m c^2) (stable at p << m c).
  double gamma = 0.0;
  double speedMPerS = 0.0;
  double rigidityV = 0.0;
  double totalEnergyJ = 0.0;
  double kineticEnergyJ = 0.0;
};
Status ComputeKinematics(const ParticleState& particle, Kinematics* out);

struct LocalState {
  // Heliocentric Cartesian position [m] (the host subtracts its solar
  // origin).  Its norm is the r of the lambda law.
  std::array<double, 3> heliocentricPositionM{{0.0, 0.0, 0.0}};
  // Coordinate time [s] and immutable provider revisions; copied to
  // provenance only (the Release A laws are time independent).
  double timeS = 0.0;
  std::uint64_t backgroundRevision = 0;
  std::uint64_t turbulenceRevision = 0;
};

// ---------------------------------------------------------------------------
// Configuration.  Zero initializers are invalid sentinels, never defaults: a
// parameter-loading utility must not choose H, eps, a step rule or a
// lambda-law value (Section 17).
struct LambdaLawParameters {
  double lambdaReferenceM = 0.0;        // lambda_ref > 0 [m]
  double referenceRigidityV = 0.0;      // R_ref > 0 [V]
  double rigidityExponent = 0.0;        // a (finite)
  double referenceRadiusM = 0.0;        // r_ref > 0 [m]
  double radialExponent = 0.0;          // b (finite)
};

struct ShapeParameters {
  // PA-QFORM: spectral index q (1<q<2) and gap constant H >= 0 (the source's
  // symbol H; h/h0 of other sources are not converted, Section 5.2).
  double q = 0.0;
  double gapH = -1.0;
  // PA-EPS: eps > 0.  eps = 0 makes lambda_par logarithmically divergent
  // (Section 3.4), so it cannot be normalized to a finite lambda law.
  double epsilon = -1.0;
};

struct StepRuleParameters {
  StepRule rule = StepRule::RateFraction;
  double rateFraction = 0.0;  // f > 0: RateFraction and VanDenBerg2026
  double ell = 0.0;           // ell > 0: VanDenBerg2026 (source: 2^-6)
  double deltaMuMax = 0.0;    // dmu_max > 0: Sandroos2014 (value not printed)
};

struct ModelConfiguration {
  ModelId model = ModelId::Isotropic;
  ShapeParameters shape;
  LambdaLawParameters lambdaLaw;
  StepRuleParameters step;
};

// Constants derived from a validated configuration (never accepted from
// input).  D = D0 * g(mu) with g = (1-mu^2) h(mu); every supremum below is a
// property of the shape alone, so the per-particle limits are D0 times it.
struct ShapeConstants {
  // N_shape with D0 = N_shape v / lambda_law (Section 5.6).
  double normalization = 0.0;
  // The closed-form value and the independent quadrature value of the shape
  // integral behind N_shape (I(q,H) or phi(eps)), and their relative
  // discrepancy.  PA-ISO/PA-EPREM are exact rationals (both fields equal).
  double integralClosedForm = 0.0;
  double integralQuadrature = 0.0;
  double integralRelativeDiscrepancy = 0.0;
  std::string closedFormMethod;  // which closed form was used
  // Ratio lambda(Eq. 23)/lambda_law: 1, except 2 for PA-EPREM.
  double lambdaEquation23OverLaw = 1.0;
  double supH = 0.0;          // sup h(mu): nu_max = 2 D0 supH
  double infH = 0.0;          // inf h(mu) (eps for PA-EPS, H for PA-QFORM)
  double supAbsGPrime = 0.0;  // sup |g'(mu)|; +inf for PA-QFORM
  double supG = 0.0;          // sup g(mu) = sup (1-mu^2) h(mu)
  ZeroClass zero = ZeroClass::Regular;
  OperatorConvention convention = OperatorConvention::Standard;
};

// Validate every field of the selected model (only the selected model's
// members participate) and compute its ShapeConstants.  The closed-form
// normalization and an independent adaptive quadrature must agree to 1e-8
// relative (roadmap Section 26.5, smooth-quadrature budget) or the call
// fails with Nonconverged.
Status ValidateConfiguration(const ModelConfiguration& configuration,
                             ShapeConstants* constants = nullptr);
// SHA-256 of the canonical serialization of the active inputs (restart and
// output identity).
std::string ConfigurationFingerprint(const ModelConfiguration& configuration);

// ---------------------------------------------------------------------------
// Evaluation result.
enum class LambdaStatus { Finite, Divergent, Unavailable, Nonconverged };
const char* Name(LambdaStatus status);

// Provenance carried by every result.  The text fields point to static
// strings and the fingerprint is a fixed buffer, so filling a result never
// allocates (the coefficient is evaluated once per particle substep).
struct Provenance {
  const char* modelId = "";
  const char* specificationVersion = "";
  const char* specificationSha256 = "";
  const char* softwareVersion = "";
  // SHA-256 hex digest of the active configuration (empty for EvaluateShape).
  char configurationFingerprint[65] = {0};
  std::uint64_t backgroundRevision = 0;
  std::uint64_t turbulenceRevision = 0;
};

// Frame and clock of the returned coefficient (Section 17, FrameContract).
struct FrameContract {
  const char* frame = "local plasma (solar-wind) frame";
  const char* muReference = "mean magnetic field B";
  const char* clock = "host time; mixed-frame, order 0 in U/c (Eq. 173 "
                      "factor not applied)";
  int retainedOrderUOverC = 0;
};

struct CoefficientResult {
  Status status;
  // D_mumu at mu [s^-1], standard operator (21).
  double dMuMuPerS = 0.0;
  // dD/dmu used by an Ito kernel [s^-1].  At mu = 0 exactly PA-EPS uses
  // Eq. (164) as printed with sgn(0) = 0 (value 0); PA-QFORM leaves it
  // absent there (unbounded).  Absent means "no finite value", never zero.
  std::optional<double> dDmuMuDmuPerS;
  // Analytic one-sided derivatives where finite (equal away from mu = 0).
  std::optional<double> derivativeMinusPerS;
  std::optional<double> derivativePlusPerS;
  double amplitudeD0PerS = 0.0;              // D0 [s^-1]
  double localScatteringFrequencyPerS = 0.0;  // 2 D/(1-mu^2) = 2 D0 h(mu)
  double maximumScatteringFrequencyPerS = 0.0;  // nu_max = 2 D0 sup h
  // Tagged mean free paths: the configured law value and the Eq. (23)
  // unfocused value of the returned coefficient (lambda_parallel_unfocused,
  // Section 16.12).  They differ only for PA-EPREM (factor 2).
  LambdaStatus lambdaStatus = LambdaStatus::Unavailable;
  double lambdaLawM = 0.0;
  double lambdaParallelUnfocusedM = 0.0;
  ZeroClass zero = ZeroClass::Regular;
  Kinematics kinematics;
  Provenance provenance;
};

// Direct evaluation with an explicit configuration (validates it first).
CoefficientResult Evaluate(double mu, const ParticleState& particle,
                           const LocalState& local,
                           const ModelConfiguration& configuration);

// Shape-level evaluation with a given amplitude D0 [s^-1] (no lambda law,
// no particle).  Used by the kernels' verification and by the finite-volume
// reference; eps = 0 and H = 0 are admitted here because no normalization is
// involved.  The returned struct carries D, D' and one-sided derivatives.
CoefficientResult EvaluateShape(double mu, double amplitudeD0PerS,
                                ModelId model, const ShapeParameters& shape);

// Allocation-free core of every evaluator above: D = D0 g(mu) and, when a
// finite value exists, D' = D0 g'(mu) (PA-EPS at mu = 0: Eq. 164 with
// sgn(0) = 0).  Returns false only when D' has no finite value (PA-QFORM at
// mu = 0); *dPrime is then left unchanged.  No validation: the caller passes
// a validated shape, |mu| <= 1 and D0 >= 0.  Intended for inner loops
// (Monte Carlo verification); hosts use the active pointer below.
bool ShapeDmumu(double mu, double amplitudeD0PerS, ModelId model,
                const ShapeParameters& shape, double* dMuMuPerS,
                double* dDmuMuDmuPerS);

// ---------------------------------------------------------------------------
// Scattering step limit for an Ito kernel (Section 16.6).
struct StepLimitResult {
  Status status;
  double maximumStepS = 0.0;
  // Which entry limited the step: "rate-fraction", "drift" or "noise"
  // (van den Berg), "rate-fraction-floor", "sandroos-drift" or
  // "sandroos-noise".
  std::string limitingEntry;
  // The individual candidates [s] (+inf where an entry does not apply):
  // the drift and noise entries of Eq. (165) or (166) and the f/nu_max
  // rate-fraction value (the floor of the van den Berg rule).
  double driftEntryS = 0.0;
  double noiseEntryS = 0.0;
  double rateFractionS = 0.0;
};

// deterministicMuRatePerS is the host's dmu/dt at mu (focusing, flow); it
// enters only the van den Berg drift entry.  A non-positive result is
// returned as Nonprogress, never as a zero-length success.
StepLimitResult ScatteringStepLimit(double mu, double deterministicMuRatePerS,
                                    const ParticleState& particle,
                                    const LocalState& local,
                                    const ModelConfiguration& configuration);

// Worst-case dimensionless step nu_max * dt that a configured step rule can
// produce for this shape (rule-independent certification measure):
//   RateFraction:   f;
//   VanDenBerg2026: max(f, ell^2 supH/infH)   (noise entry bounds the step);
//   Sandroos2014:   min(dmu_max * 2 supH / sup|g'|, dmu_max^2 supH / supG).
// PA-QFORM has sup|g'| = inf, and infH = 0 when H = 0: the bound is then
// infinite for the rules that need it.
double WorstCaseScaledStep(const ModelConfiguration& configuration,
                           const ShapeConstants& constants);

// ---------------------------------------------------------------------------
// Coefficient-kernel readiness (Sections 16.11, 26.3 PS05).  Release is a
// property of a (model, kernel) pair and of the step it is run with.  The
// matrix is the evidence of the PS05 acceptance campaign recorded in
// IMPLEMENTATION_STATUS.md; outside it the call fails with
// IncompatibleKernel and names the reason.
struct PairCertification {
  ModelId model;
  Kernel kernel;
  // Largest certified dimensionless step nu_max * dt.
  double maximumScaledStep;
  // Certified PA-EPS range (ignored for the other models).
  double epsilonMinimum;
  double epsilonMaximum;
  const char* evidence;  // test IDs that certify the entry
};
const std::vector<PairCertification>& CertifiedPairs();
Status CheckKernelPair(const ModelConfiguration& configuration, Kernel kernel);

// ---------------------------------------------------------------------------
// Parser-neutral input bridge (same contract as the parallel-diffusion
// library).  Keys are case-sensitive and name SI units; values are SI text
// without unit suffixes.  A value may be an exact dimensionless rational
// "a/b" (for example q = 5/3; Section 26.3 PS01); every other value is a
// decimal number.  Unknown, duplicate, missing and inactive keys fail.
struct InputParameter {
  std::string name;
  std::string value;
};
Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          ModelConfiguration* configuration);

// Library-owned section grammar for the host block
// '#begin physics: pitch_angle_diffusion'.  The host passes the body lines
// with their physical line numbers; '#' or '!' starts a comment; every other
// nonblank line is exactly one 'key = value'; 'model' selects the inventory
// ID (case-insensitive) and must appear once.  Diagnostics begin with
// "line N:".  Keys of the block:
//
//   model                        PA-ISO | PA-QFORM | PA-EPS | PA-EPREM
//   q                            PA-QFORM: 1 < q < 2 (rational allowed)
//   gap_h                        PA-QFORM: H >= 0
//   epsilon                      PA-EPS:   eps > 0
//   lambda_reference_m           lambda_ref > 0 [m]
//   lambda_reference_rigidity_v  R_ref > 0 [V]
//   lambda_rigidity_exponent     a
//   lambda_reference_radius_m    r_ref > 0 [m]
//   lambda_radial_exponent       b
//   step_rule                    rate-fraction | van-den-berg-2026 |
//                                sandroos-2014
//   rate_fraction                f > 0 (rate-fraction, van-den-berg-2026)
//   ell                          ell > 0 (van-den-berg-2026)
//   delta_mu_max                 0 < dmu_max <= 2 (sandroos-2014)
struct SectionLine {
  std::size_t lineNumber = 0;  // 0 = unknown (no "line N:" prefix)
  std::string text;
};
struct SectionParameter {
  std::string name;
  std::string value;
  std::size_t lineNumber = 0;
};
struct ParsedSection {
  std::string modelId;
  std::size_t modelLineNumber = 0;
  std::vector<SectionParameter> parameters;
  ModelConfiguration configuration;
  ShapeConstants constants;
  std::string configurationFingerprint;
};
// Transactional: *result is written only after the whole section validates.
Status ParseSection(const std::vector<SectionLine>& lines,
                    ParsedSection* result);

// ---------------------------------------------------------------------------
// Active model and mover-facing dispatch pointers.  Configure once during
// serial initialization (before any mover thread).  Before a successful
// configuration both pointers target evaluators that return
// InvalidConfiguration, never a coefficient.  A focused-transport pusher
// needs no model switch:
//
//   CoefficientResult r = ActivePitchAngleDiffusion(mu, particle, local);
//   if (!r.status.ok()) { /* typed failure; never substitute a value */ }
//   StepLimitResult s = ActiveScatteringStepLimit(mu, dmu_dt, particle, local);
Status SetActiveConfiguration(const ModelConfiguration& configuration);
Status ConfigureActiveModel(const std::string& modelId,
                            const std::vector<InputParameter>& parameters);
Status ConfigureActiveModelFromSection(const std::vector<SectionLine>& lines,
                                       ParsedSection* result);
ModelConfiguration GetActiveConfiguration();
ShapeConstants GetActiveShapeConstants();
bool HasActiveConfiguration();

using ActiveCoefficientEvaluator = CoefficientResult (*)(
    double mu, const ParticleState&, const LocalState&);
using ActiveStepLimitEvaluator = StepLimitResult (*)(
    double mu, double deterministicMuRatePerS, const ParticleState&,
    const LocalState&);
extern ActiveCoefficientEvaluator ActivePitchAngleDiffusion;
extern ActiveStepLimitEvaluator ActiveScatteringStepLimit;
// The evaluator ActivePitchAngleDiffusion takes when `model` is active.
ActiveCoefficientEvaluator BoundFunctionForModel(ModelId model);

// ---------------------------------------------------------------------------
// Numerical building blocks exposed for verification (pure functions).
//
// I(q,H) = Int_0^1 (1-mu^2)/(mu^(q-1)+H) dmu (Section 5).  Closed forms:
// H = 0 exactly 2/((2-q)(4-q)); q = 5/3 Eq. (38) for H <= 10; otherwise the
// convergent second form of Eq. (37), 2F1(1,1;c;z) with z = 1/(1+H), summed
// with a rigorous geometric tail bound.  Nonconverged if the series would
// need more than 2e8 terms (H < ~1e-7 for q != 5/3).
Status QFormIntegralClosedForm(double q, double gapH, double* value,
                               std::string* method);
// Independent adaptive Gauss-Kronrod quadrature of the same integral after
// the substitution mu = t^(1/(2-q)), which removes the |mu|^(1-q) endpoint
// behaviour.
Status QFormIntegralQuadrature(double q, double gapH, double* value);
// phi(eps) = nu0 lambda/v of PA-EPS, Eq. (39) (closed form) and an
// independent quadrature after mu = (eps/a)(e^u - 1).
Status EpsilonPhiClosedForm(double epsilon, double* value);
Status EpsilonPhiQuadrature(double epsilon, double* value);

}  // namespace PitchAngleDiffusion
}  // namespace SEP

#endif  // SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_DIFFUSION_H

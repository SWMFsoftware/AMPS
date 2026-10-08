#include "parallel_diffusion.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <sstream>

namespace SEP {
namespace ParallelDiffusion {
namespace {

constexpr const char* SpecificationVersion = "1.3";
constexpr const char* SoftwareVersion = "parallel-diffusion-pd02-v1";

// The active configuration is process-local state selected during serial
// startup.  It is intentionally not protected by a mutex: the public contract
// forbids reconfiguration after mover threads begin, while direct Evaluate()
// remains reentrant because its configuration is passed by const reference.
ModelConfiguration gActiveConfiguration;

bool FinitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

double Magnitude(const std::array<double, 3>& value) {
  // Nested hypot avoids the unnecessary overflow/underflow risk of forming
  // x*x+y*y+z*z.  Callers separately decide whether zero is physically valid.
  return std::hypot(value[0], std::hypot(value[1], value[2]));
}

double LogSumExp(double left, double right) {
  // max(a,b)+log(exp(a-max)+exp(b-max)) is algebraically log(exp(a)+exp(b))
  // but keeps both exponent arguments non-positive.  This is the stability
  // primitive required by the smooth broken-rigidity law.
  const double largest = std::max(left, right);
  return largest + std::log(std::exp(left - largest) +
                            std::exp(right - largest));
}

ParallelResult Failure(StatusCode code, const std::string& detail,
                       const ModelConfiguration& configuration,
                       const LocalState& local) {
  // A failed result deliberately carries no coefficient optionals.  Provenance
  // is still populated so a host can diagnose which configuration and
  // provider generations rejected the state without mistaking failure for a
  // zero coefficient.
  ParallelResult result;
  result.status = Status::Error(code, detail);
  result.provenance.requestedModelId = ModelName(configuration.model);
  result.provenance.evaluatedModelId = ModelName(configuration.model);
  result.provenance.specificationVersion = SpecificationVersion;
  result.provenance.softwareVersion = SoftwareVersion;
  result.provenance.configurationFingerprint =
      ConfigurationFingerprint(configuration);
  result.provenance.backgroundRevision = local.backgroundRevision;
  result.provenance.turbulenceRevision = local.turbulenceRevision;
  return result;
}

ParallelResult Success(const ModelConfiguration& configuration,
                       const LocalState& local, double lambdaM,
                       double kappaM2PerS) {
  // Every implemented PD02 law is a finite-scattering model with strictly
  // positive lambda and kappa.  Infinity and underflow-to-zero therefore mean
  // that a mathematically valid-looking input left the representable domain;
  // neither value is clamped because that would add unconfigured physics.
  if (!FinitePositive(lambdaM) || !FinitePositive(kappaM2PerS)) {
    return Failure(StatusCode::OutsideModelDomain,
        "parallel coefficient evaluation overflowed or was non-positive",
        configuration, local);
  }
  ParallelResult result;
  result.status = Status::Success();
  result.lambdaParallelM = lambdaM;
  result.kappaParallelM2PerS = kappaM2PerS;
  result.provenance.requestedModelId = ModelName(configuration.model);
  result.provenance.evaluatedModelId = ModelName(configuration.model);
  result.provenance.specificationVersion = SpecificationVersion;
  result.provenance.softwareVersion = SoftwareVersion;
  result.provenance.configurationFingerprint =
      ConfigurationFingerprint(configuration);
  result.provenance.backgroundRevision = local.backgroundRevision;
  result.provenance.turbulenceRevision = local.turbulenceRevision;
  return result;
}

Status MeanFieldMagnitude(const LocalState& local, double* magnitudeT) {
  if (!magnitudeT) {
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null mean-field magnitude output");
  }
  if (!local.meanFieldT.has_value()) {
    return Status::Error(StatusCode::MissingInput,
                         "selected model requires mean_B_T");
  }
  // meanFieldT is the resolved field defining the local field-aligned basis,
  // not a fluctuation RMS or an automatically constructed total field.
  const std::array<double, 3>& field = *local.meanFieldT;
  if (!std::isfinite(field[0]) || !std::isfinite(field[1]) ||
      !std::isfinite(field[2])) {
    return Status::Error(StatusCode::InvalidBackground,
                         "mean_B_T contains a non-finite component");
  }
  *magnitudeT = Magnitude(field);
  return FinitePositive(*magnitudeT)
      ? Status::Success()
      : Status::Error(StatusCode::InvalidBackground,
                      "field-aligned diffusion requires |mean_B_T|>0");
}

Status RequiredFactor(const std::optional<double>& input,
                      const char* name, double* value) {
  // Provider factors enter multiplicatively and are evaluated in logarithmic
  // form, hence the strictly positive domain.  Zero is not silently treated
  // as absence, and absence is not silently replaced by unity once enabled.
  if (!input.has_value()) {
    return Status::Error(StatusCode::MissingInput,
                         std::string("selected model requires ") + name);
  }
  if (!FinitePositive(*input)) {
    return Status::Error(StatusCode::InvalidBackground,
                         std::string(name) + " must be finite and positive");
  }
  *value = *input;
  return Status::Success();
}

ParallelResult EvaluateUnconfigured(const ParticleState&, const LocalState& local,
                                    const ModelConfiguration& configuration) {
  // This sentinel makes an accidental pre-initialization EvaluateActive call
  // fail deterministically instead of invoking a null function pointer or a
  // physically invalid zero-initialized ConstantLambda configuration.
  return Failure(StatusCode::InvalidConfiguration,
                 "no validated active parallel-diffusion model is configured",
                 configuration, local);
}

ParallelResult EvaluateConstantLambda(
    const ParticleState& particle, const LocalState& local,
    const ModelConfiguration& configuration) {
  ParticleKinematics kinematics;
  const Status particleStatus =
      ComputeParticleKinematics(particle, &kinematics);
  if (!particleStatus.ok())
    return Failure(particleStatus.code, particleStatus.detail,
                   configuration, local);

  // Equation (11): lambda is independent of position and momentum, while
  // kappa=v*lambda/3 retains the exact relativistic particle speed.  The
  // spatial derivative is therefore exactly zero at fixed momentum.  For a
  // fixed species, rigidity is proportional to p and
  // d ln(v)/d ln(p)=1/gamma^2, which supplies the kappa rigidity slope.
  ParallelResult result = Success(configuration, local,
      configuration.constantLambda.lambdaParallelM,
      kinematics.speedMPerS *
          configuration.constantLambda.lambdaParallelM / 3.0);
  if (result.status.ok()) {
    result.dLnLambdaDLnRigidity = 0.0;
    result.dLnKappaDLnRigidity =
        1.0 / (kinematics.gamma * kinematics.gamma);
    result.gradKappaParallelMPerS = std::array<double, 3>{{0.0, 0.0, 0.0}};
  }
  return result;
}

ParallelResult EvaluateConstantKappa(
    const ParticleState& particle, const LocalState& local,
    const ModelConfiguration& configuration) {
  ParticleKinematics kinematics;
  const Status particleStatus =
      ComputeParticleKinematics(particle, &kinematics);
  if (!particleStatus.ok())
    return Failure(particleStatus.code, particleStatus.detail,
                   configuration, local);

  // Equation (12) is deliberately not implemented through Equation (11): a
  // constant kappa has lambda=3*kappa/v and hence a different speed law.
  // Thus d ln(lambda)/d ln(R)=-1/gamma^2 at fixed species, while the declared
  // kappa and its spatial gradient are exactly constant.
  ParallelResult result = Success(configuration, local,
      3.0 * configuration.constantKappa.kappaParallelM2PerS /
          kinematics.speedMPerS,
      configuration.constantKappa.kappaParallelM2PerS);
  if (result.status.ok()) {
    result.dLnLambdaDLnRigidity =
        -1.0 / (kinematics.gamma * kinematics.gamma);
    result.dLnKappaDLnRigidity = 0.0;
    result.gradKappaParallelMPerS = std::array<double, 3>{{0.0, 0.0, 0.0}};
  }
  return result;
}

ParallelResult EvaluatePowerLawLambda(
    const ParticleState& particle, const LocalState& local,
    const ModelConfiguration& configuration) {
  ParticleKinematics kinematics;
  const Status particleStatus =
      ComputeParticleKinematics(particle, &kinematics);
  if (!particleStatus.ok())
    return Failure(particleStatus.code, particleStatus.detail,
                   configuration, local);

  const PowerLawLambdaParameters& p = configuration.powerLawLambda;
  double independent = 0.0;
  double dLnIndependentDLnRigidity = 0.0;
  switch (p.independentVariable) {
    case IndependentVariable::Rigidity:
      independent = kinematics.rigidityV;
      dLnIndependentDLnRigidity = 1.0;
      break;
    case IndependentVariable::TotalKineticEnergy:
      independent = kinematics.kineticEnergyJ;
      // Since dT/dp=v and R is proportional to p for a fixed species,
      // d ln(T)/d ln(R)=p*v/T.  This exact expression spans both the
      // nonrelativistic limit (2) and ultrarelativistic limit (1).
      dLnIndependentDLnRigidity =
          particle.momentumKgMPerS * kinematics.speedMPerS /
          kinematics.kineticEnergyJ;
      break;
    case IndependentVariable::EnergyPerNucleon:
      if (!particle.nucleonCount.has_value() ||
          !FinitePositive(*particle.nucleonCount)) {
        return Failure(StatusCode::MissingInput,
            "energy_per_nucleon requires an explicit positive nucleon_count",
            configuration, local);
      }
      independent = kinematics.kineticEnergyJ / *particle.nucleonCount;
      // Nucleon count is constant under the rigidity derivative, so division
      // by A changes the normalization but not the logarithmic slope.
      dLnIndependentDLnRigidity =
          particle.momentumKgMPerS * kinematics.speedMPerS /
          kinematics.kineticEnergyJ;
      break;
    case IndependentVariable::Speed:
      independent = kinematics.speedMPerS;
      // Exact fixed-species identity from v=cp/sqrt((mc)^2+p^2).
      dLnIndependentDLnRigidity =
          1.0 / (kinematics.gamma * kinematics.gamma);
      break;
  }

  // Equations (13)--(15) are accumulated in logarithmic form:
  //   ln(lambda)=ln(lambda0)+a ln(X/X0)+alpha ln(r/r0)
  //              -eta ln(B/Bref)+ln(g_t)+ln(g_region).
  // This preserves the stated separable prescription without overflowing
  // intermediate powers.  It does not remove the final double-precision
  // range check in Success().  Optional factors introduce no state
  // requirement when disabled.
  double logLambda = std::log(p.lambda0M) +
      p.independentExponent *
          std::log(independent / p.independentReferenceSI);
  if (p.useRadialFactor) {
    // LocalState::positionM is contractually heliocentric.  The evaluator has
    // no mesh origin and therefore cannot correct a translated host position.
    const double radiusM = Magnitude(local.positionM);
    if (!FinitePositive(radiusM))
      return Failure(StatusCode::InvalidBackground,
          "radial power law requires a finite positive heliocentric radius",
          configuration, local);
    logLambda += p.radialExponent * std::log(radiusM / p.radius0M);
  }
  // Equation (13) defines eta=0 as a disabled field factor.  Honor that
  // identity before looking for B so a mathematically inactive term creates
  // no background-provider requirement.
  if (p.useFieldFactor && p.fieldExponent != 0.0) {
    double fieldT = 0.0;
    const Status field = MeanFieldMagnitude(local, &fieldT);
    if (!field.ok())
      return Failure(field.code, field.detail, configuration, local);
    logLambda -= p.fieldExponent *
        std::log(fieldT / p.fieldReferenceT);
  }
  if (p.useTimeFactor) {
    double factor = 0.0;
    const Status status = RequiredFactor(local.timeFactor, "time_factor", &factor);
    if (!status.ok())
      return Failure(status.code, status.detail, configuration, local);
    logLambda += std::log(factor);
  }
  if (p.useRegionFactor) {
    double factor = 0.0;
    const Status status =
        RequiredFactor(local.regionFactor, "region_factor", &factor);
    if (!status.ok())
      return Failure(status.code, status.detail, configuration, local);
    logLambda += std::log(factor);
  }

  const double lambdaM = std::exp(logLambda);
  ParallelResult result = Success(configuration, local, lambdaM,
                                  kinematics.speedMPerS * lambdaM / 3.0);
  if (result.status.ok()) {
    // Enabled spatial factors can make grad(kappa) nonzero, but the current
    // LocalState does not yet carry their coherent gradients.  Preserve the
    // scalar and analytic rigidity slopes while explicitly withholding the
    // spatial derivative needed by a complete stochastic drift.
    const double lambdaSlope =
        p.independentExponent * dLnIndependentDLnRigidity;
    result.dLnLambdaDLnRigidity = lambdaSlope;
    result.dLnKappaDLnRigidity = lambdaSlope +
        1.0 / (kinematics.gamma * kinematics.gamma);
    result.diagnosticMask |= DerivativeUnavailable;
  }
  return result;
}

ParallelResult EvaluateBrokenRigidityKappa(
    const ParticleState& particle, const LocalState& local,
    const ModelConfiguration& configuration) {
  ParticleKinematics kinematics;
  const Status particleStatus =
      ComputeParticleKinematics(particle, &kinematics);
  if (!particleStatus.ok())
    return Failure(particleStatus.code, particleStatus.detail,
                   configuration, local);

  const BrokenRigidityKappaParameters& p =
      configuration.brokenRigidityKappa;
  const double logR = std::log(kinematics.rigidityV / p.rigidity0V);
  const double logRb = std::log(p.breakRigidityV / p.rigidity0V);

  // Equation (17) in log-sum-exp form.  With x=ln(R/R0) and xb=ln(Rb/R0),
  // ln(H)=a*x+(b-a)/h * {log(exp(h*x)+exp(h*xb))
  //                       -log(1+exp(h*xb))}.
  // Directly raising R and R_b to h can overflow even when this normalized
  // transition and H are representable.
  const double hLogR = p.smoothness * logR;
  const double hLogRb = p.smoothness * logRb;
  const double logTransition =
      LogSumExp(hLogR, hLogRb) - LogSumExp(0.0, hLogRb);
  double logH = p.lowSlope * logR +
      (p.highSlope - p.lowSlope) / p.smoothness * logTransition;

  double logFactors = 0.0;
  if (p.useFieldFactor && p.fieldExponent != 0.0) {
    double fieldT = 0.0;
    const Status field = MeanFieldMagnitude(local, &fieldT);
    if (!field.ok())
      return Failure(field.code, field.detail, configuration, local);
    // Equation (17) uses (Bref/B)^eta, hence the positive eta multiplier on
    // ln(Bref/B).  This is equivalent to the minus sign in Equation (13).
    logFactors += p.fieldExponent *
        std::log(p.fieldReferenceT / fieldT);
  }
  if (p.useRadialFactor) {
    double factor = 0.0;
    const Status status =
        RequiredFactor(local.radialFactor, "radial_factor", &factor);
    if (!status.ok())
      return Failure(status.code, status.detail, configuration, local);
    logFactors += std::log(factor);
  }
  if (p.useRegionFactor) {
    double factor = 0.0;
    const Status status =
        RequiredFactor(local.regionFactor, "region_factor", &factor);
    if (!status.ok())
      return Failure(status.code, status.detail, configuration, local);
    logFactors += std::log(factor);
  }

  // Equations (17),(18): beta appears in kappa=K_star*beta*H*factors.
  // Substitution in lambda=3*kappa/v with v=beta*c cancels beta exactly,
  // making equal-rigidity species share lambda while their kappas can differ.
  const double hAndFactors = std::exp(logH + logFactors);
  const double kappaM2PerS =
      p.kStarM2PerS * kinematics.beta * hAndFactors;
  const double lambdaM =
      3.0 * p.kStarM2PerS / SpeedOfLightMPerS * hAndFactors;
  ParallelResult result = Success(configuration, local, lambdaM, kappaM2PerS);
  if (result.status.ok()) {
    // Equation (19) uses w=[1+(Rb/R)^h]^-1.  Evaluate the logistic with a
    // bounded exponential argument: values beyond +/-700 are already at the
    // double-precision asymptote, and bounding avoids an otherwise irrelevant
    // overflow while preserving the limiting slopes a and b.
    const double transitionWeight =
        1.0 / (1.0 + std::exp(std::clamp(
            p.smoothness * std::log(p.breakRigidityV /
                                    kinematics.rigidityV),
            -700.0, 700.0)));
    const double lambdaSlope = p.lowSlope +
        (p.highSlope - p.lowSlope) * transitionWeight;
    result.dLnLambdaDLnRigidity = lambdaSlope;
    result.dLnKappaDLnRigidity = lambdaSlope +
        1.0 / (kinematics.gamma * kinematics.gamma);
    result.diagnosticMask |= DerivativeUnavailable;
  }
  return result;
}

ParallelResult EvaluateBohm(
    const ParticleState& particle, const LocalState& local,
    const ModelConfiguration& configuration) {
  ParticleKinematics kinematics;
  const Status particleStatus =
      ComputeParticleKinematics(particle, &kinematics);
  if (!particleStatus.ok())
    return Failure(particleStatus.code, particleStatus.detail,
                   configuration, local);

  double fieldT = 0.0;
  if (configuration.bohm.fieldDefinition == BohmFieldDefinition::MeanField) {
    const Status field = MeanFieldMagnitude(local, &fieldT);
    if (!field.ok())
      return Failure(field.code, field.detail, configuration, local);
  } else {
    if (!local.effectiveFieldMagnitudeT.has_value())
      return Failure(StatusCode::MissingInput,
          "effective-field Bohm model requires effective_field_magnitude_T",
          configuration, local);
    fieldT = *local.effectiveFieldMagnitudeT;
    if (!FinitePositive(fieldT))
      return Failure(StatusCode::InvalidBackground,
          "effective_field_magnitude_T must be finite and positive",
          configuration, local);
  }

  // Equation (61), using the maximum (90-degree pitch-angle) gyroradius
  // r_L=p/(|q|B).  Charge sign is intentionally removed for this scalar
  // balanced prescription.  eta_B is an explicitly supplied comparison
  // parameter; this evaluator never clamps a different backend to the Bohm
  // value and never claims Bohm is a universal lower bound.
  const double larmorRadiusM = particle.momentumKgMPerS /
      (std::fabs(particle.chargeC) * fieldT);
  const double lambdaM = configuration.bohm.etaB * larmorRadiusM;
  ParallelResult result = Success(configuration, local, lambdaM,
                                  kinematics.speedMPerS * lambdaM / 3.0);
  if (result.status.ok()) {
    result.dLnLambdaDLnRigidity = 1.0;
    result.dLnKappaDLnRigidity = 1.0 +
        1.0 / (kinematics.gamma * kinematics.gamma);
    result.diagnosticMask |= DerivativeUnavailable;
  }
  return result;
}

bool ParseDouble(const std::string& text, double* value) {
  // strtod is used only as a syntax conversion.  Full consumption, ERANGE,
  // and finiteness checks make inputs such as "1 m", NaN, infinity, and
  // under/overflow failures explicit instead of partially accepting them.
  if (!value || text.empty()) return false;
  char* end = nullptr;
  errno = 0;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ParseBool(const std::string& text, bool* value) {
  // The bridge accepts only explicit true/false (case-insensitive); numeric
  // booleans and implicit presence flags are intentionally not guessed.
  if (!value) return false;
  const std::string normalized = Lower(text);
  if (normalized == "true") *value = true;
  else if (normalized == "false") *value = false;
  else return false;
  return true;
}

Status InputMap(const std::vector<InputParameter>& parameters,
                std::map<std::string, std::string>* values) {
  if (!values) return Status::Error(StatusCode::InvalidConfiguration,
                                    "null parsed-parameter output");
  // emplace, rather than assignment, preserves duplicate-key detection.  A
  // file parser must not let the last spelling silently replace an earlier
  // physical parameter.
  for (const InputParameter& parameter : parameters) {
    if (parameter.name.empty())
      return Status::Error(StatusCode::InvalidConfiguration,
                           "parallel-diffusion parameter name is empty");
    if (!values->emplace(parameter.name, parameter.value).second)
      return Status::Error(StatusCode::InvalidConfiguration,
          "duplicate parallel-diffusion parameter '" + parameter.name + "'");
  }
  return Status::Success();
}

Status TakeDouble(std::map<std::string, std::string>* values,
                  const char* key, bool required, double* output) {
  // Consumed entries are erased.  RejectUnknown can then prove that every
  // supplied key belongs to the selected model's exact schema.
  const auto found = values->find(key);
  if (found == values->end()) {
    return required
        ? Status::Error(StatusCode::MissingInput,
                        std::string("missing parameter '") + key + "'")
        : Status::Success();
  }
  if (!ParseDouble(found->second, output))
    return Status::Error(StatusCode::InvalidConfiguration,
        std::string("parameter '") + key + "' is not a finite number");
  values->erase(found);
  return Status::Success();
}

Status TakeBool(std::map<std::string, std::string>* values,
                const char* key, bool required, bool* output) {
  const auto found = values->find(key);
  if (found == values->end()) {
    return required
        ? Status::Error(StatusCode::MissingInput,
                        std::string("missing parameter '") + key + "'")
        : Status::Success();
  }
  if (!ParseBool(found->second, output))
    return Status::Error(StatusCode::InvalidConfiguration,
        std::string("parameter '") + key + "' must be true or false");
  values->erase(found);
  return Status::Success();
}

Status RejectUnknown(const std::map<std::string, std::string>& values) {
  if (values.empty()) return Status::Success();
  return Status::Error(StatusCode::InvalidConfiguration,
      "unknown parameter '" + values.begin()->first + "'");
}

void FingerprintDouble(std::ostringstream* output, const char* name,
                       double value) {
  *output << ';' << name << '=' << std::setprecision(17) << value;
}

std::uint64_t Fnv1a(const std::string& text) {
  // This deterministic hash is a compact configuration identity, not a
  // cryptographic checksum.  Future external tables retain independent
  // SHA-256 provenance as required by the specification.
  std::uint64_t hash = UINT64_C(1469598103934665603);
  for (unsigned char value : text) {
    hash ^= value;
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

}  // namespace

Status Status::Success() { return Status(); }

Status Status::Error(StatusCode code, const std::string& detail) {
  Status status;
  status.code = code;
  status.detail = detail;
  return status;
}

const std::vector<ModelDescriptor>& ModelRegistry() {
  // Keep the complete v1.3 identifier inventory visible even when a backend
  // is deferred.  This lets configuration handling distinguish a known but
  // unavailable physical model from a misspelling and prevents an implicit
  // fallback to the nearest implemented law.
  static const std::vector<ModelDescriptor> models = {
      {ModelId::ConstantLambda, "constant_lambda", true, "lambda,kappa", "PD02"},
      {ModelId::ConstantKappa, "constant_kappa", true, "kappa,lambda", "PD02"},
      {ModelId::PowerLawLambda, "power_law_lambda", true, "lambda,kappa", "PD02"},
      {ModelId::BrokenRigidityKappa, "broken_rigidity_kappa", true, "kappa,lambda", "PD02"},
      {ModelId::QltSlabSpectrum, "qlt_slab_spectrum", false, "D_mumu,lambda,kappa", "PD04"},
      {ModelId::QltSlabInertial, "qlt_slab_inertial", false, "lambda,kappa", "PD04"},
      {ModelId::PrescribedLambdaMuShape, "prescribed_lambda_mu_shape", false, "D_mumu,lambda,kappa", "PD04"},
      {ModelId::BroadenedSlab, "broadened_slab", false, "D_mumu,lambda,kappa", "PD07"},
      {ModelId::NlpaGivenPerp, "nlpa_given_perp", false, "kappa_parallel,lambda_parallel", "PD06"},
      {ModelId::NlgcE, "nlgc_e", false, "parallel,perpendicular", "PD06"},
      {ModelId::NlgceN, "nlgce_n", false, "parallel,perpendicular", "PD06"},
      {ModelId::NlgceF2014, "nlgce_f_2014", false, "parallel,perpendicular", "PD05"},
      {ModelId::TurbulenceAdapter, "turbulence_adapter", false, "selected closure", "PD08"},
      {ModelId::WaveSpectrumAdapter, "wave_spectrum_adapter", false, "selected closure", "PD08"},
      {ModelId::Bohm, "bohm", true, "lambda,kappa", "PD02"},
      {ModelId::TabulatedParallel, "tabulated_parallel", false, "lambda or kappa", "PD08"}};
  return models;
}

const char* ModelName(ModelId model) {
  for (const ModelDescriptor& descriptor : ModelRegistry())
    if (descriptor.id == model) return descriptor.stableId;
  return "unknown";
}

bool ParseModelId(const std::string& text, ModelId* model) {
  if (!model) return false;
  // Stable IDs are case-insensitive but otherwise exact: no punctuation,
  // whitespace, alias, or legacy-name normalization occurs in the physics
  // library.  Compatibility aliases, if required, belong in a host adapter.
  const std::string normalized = Lower(text);
  for (const ModelDescriptor& descriptor : ModelRegistry()) {
    if (normalized == descriptor.stableId) {
      *model = descriptor.id;
      return true;
    }
  }
  return false;
}

Status ComputeParticleKinematics(const ParticleState& particle,
                                 ParticleKinematics* kinematics) {
  if (!kinematics)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null particle-kinematics output");
  if (!FinitePositive(particle.massKg) || !std::isfinite(particle.chargeC) ||
      particle.chargeC == 0.0 || !FinitePositive(particle.momentumKgMPerS)) {
    return Status::Error(StatusCode::InvalidParticle,
        "particle requires positive mass and momentum and nonzero finite charge");
  }
  // Momentum form of Equations (7),(8).  hypot(mc,p) evaluates
  // sqrt((mc)^2+p^2) without forming squares that can overflow.  Then
  // gamma=hypot(mc,p)/(mc), beta=p/hypot(mc,p), and v=c*beta.
  const double mc = particle.massKg * SpeedOfLightMPerS;
  const double denominator = std::hypot(mc, particle.momentumKgMPerS);
  const double mc2 = mc * SpeedOfLightMPerS;
  const double pc = particle.momentumKgMPerS * SpeedOfLightMPerS;
  const double totalEnergyJ = std::hypot(mc2, pc);
  kinematics->gamma = denominator / mc;
  kinematics->beta = particle.momentumKgMPerS / denominator;
  kinematics->speedMPerS = SpeedOfLightMPerS * kinematics->beta;
  // Rigidity is positive pc/|q| in volts.  The signed charge remains in the
  // source ParticleState for future models in which polarity matters.
  kinematics->rigidityV = pc / std::fabs(particle.chargeC);
  // T=(pc)^2/(E+mc^2) is algebraically identical to E-mc^2 and avoids
  // catastrophic cancellation for non-relativistic particles.
  kinematics->kineticEnergyJ = pc * pc / (totalEnergyJ + mc2);
  if (!FinitePositive(kinematics->gamma) ||
      !FinitePositive(kinematics->beta) ||
      !(kinematics->beta < 1.0) ||
      !FinitePositive(kinematics->speedMPerS) ||
      !FinitePositive(kinematics->rigidityV) ||
      !FinitePositive(kinematics->kineticEnergyJ)) {
    return Status::Error(StatusCode::InvalidParticle,
                         "relativistic particle conversion overflowed");
  }
  return Status::Success();
}

Status ValidateConfiguration(const ModelConfiguration& configuration) {
  // Validation is deliberately selective: only fields used by the chosen
  // model are inspected, so a constant law does not acquire a magnetic-field
  // or heliocentric-position dependency.  Local/provider inputs are checked
  // later, also selectively, at evaluation time.
  switch (configuration.model) {
    case ModelId::ConstantLambda:
      return FinitePositive(configuration.constantLambda.lambdaParallelM)
          ? Status::Success()
          : Status::Error(StatusCode::InvalidConfiguration,
                          "constant_lambda requires lambda_parallel_m>0");
    case ModelId::ConstantKappa:
      return FinitePositive(configuration.constantKappa.kappaParallelM2PerS)
          ? Status::Success()
          : Status::Error(StatusCode::InvalidConfiguration,
                          "constant_kappa requires kappa_parallel_m2_per_s>0");
    case ModelId::PowerLawLambda: {
      const PowerLawLambdaParameters& p = configuration.powerLawLambda;
      if (!FinitePositive(p.lambda0M) ||
          !FinitePositive(p.independentReferenceSI) ||
          !std::isfinite(p.independentExponent))
        return Status::Error(StatusCode::InvalidConfiguration,
            "power_law_lambda requires positive lambda/reference and a finite exponent");
      if (p.useRadialFactor &&
          (!FinitePositive(p.radius0M) || !std::isfinite(p.radialExponent)))
        return Status::Error(StatusCode::InvalidConfiguration,
            "enabled radial factor requires radius0_m>0 and finite radial_exponent");
      if (p.useFieldFactor && p.fieldExponent != 0.0 &&
          (!FinitePositive(p.fieldReferenceT) || !std::isfinite(p.fieldExponent)))
        return Status::Error(StatusCode::InvalidConfiguration,
            "enabled field factor requires field_reference_T>0 and finite field_exponent");
      return Status::Success();
    }
    case ModelId::BrokenRigidityKappa: {
      const BrokenRigidityKappaParameters& p =
          configuration.brokenRigidityKappa;
      if (!FinitePositive(p.kStarM2PerS) || !FinitePositive(p.rigidity0V) ||
          !FinitePositive(p.breakRigidityV) ||
          !std::isfinite(p.lowSlope) || !std::isfinite(p.highSlope) ||
          !FinitePositive(p.smoothness))
        return Status::Error(StatusCode::InvalidConfiguration,
            "broken_rigidity_kappa requires positive K_star, rigidities and smoothness, and finite slopes");
      if (p.useFieldFactor && p.fieldExponent != 0.0 &&
          (!FinitePositive(p.fieldReferenceT) || !std::isfinite(p.fieldExponent)))
        return Status::Error(StatusCode::InvalidConfiguration,
            "enabled field factor requires field_reference_T>0 and finite field_exponent");
      return Status::Success();
    }
    case ModelId::Bohm:
      return FinitePositive(configuration.bohm.etaB)
          ? Status::Success()
          : Status::Error(StatusCode::InvalidConfiguration,
                          "bohm requires eta_B>0");
    default:
      return Status::Error(StatusCode::UnsupportedModel,
          std::string("model '") + ModelName(configuration.model) +
          "' is specified but not implemented in the PD02 release slice");
  }
}

std::string ConfigurationFingerprint(const ModelConfiguration& configuration) {
  // Serialize only the active model schema in a fixed order with enough
  // decimal precision to round-trip a binary64 value.  Inactive union-like
  // parameter members do not change physical behavior and therefore do not
  // change this identity.  This routine can also describe an invalid
  // configuration in a failure record; validity remains a separate status.
  std::ostringstream canonical;
  canonical << "spec=" << SpecificationVersion << ";model="
            << ModelName(configuration.model);
  switch (configuration.model) {
    case ModelId::ConstantLambda:
      FingerprintDouble(&canonical, "lambda_parallel_m",
                        configuration.constantLambda.lambdaParallelM);
      break;
    case ModelId::ConstantKappa:
      FingerprintDouble(&canonical, "kappa_parallel_m2_per_s",
                        configuration.constantKappa.kappaParallelM2PerS);
      break;
    case ModelId::PowerLawLambda: {
      const PowerLawLambdaParameters& p = configuration.powerLawLambda;
      FingerprintDouble(&canonical, "lambda0_m", p.lambda0M);
      canonical << ";independent_variable="
                << static_cast<int>(p.independentVariable);
      FingerprintDouble(&canonical, "independent_reference_si",
                        p.independentReferenceSI);
      FingerprintDouble(&canonical, "independent_exponent",
                        p.independentExponent);
      canonical << ";use_radial_factor=" << p.useRadialFactor
                << ";use_field_factor=" << p.useFieldFactor
                << ";use_time_factor=" << p.useTimeFactor
                << ";use_region_factor=" << p.useRegionFactor;
      if (p.useRadialFactor) {
        FingerprintDouble(&canonical, "radius0_m", p.radius0M);
        FingerprintDouble(&canonical, "radial_exponent", p.radialExponent);
      }
      if (p.useFieldFactor && p.fieldExponent != 0.0) {
        FingerprintDouble(&canonical, "field_reference_T", p.fieldReferenceT);
        FingerprintDouble(&canonical, "field_exponent", p.fieldExponent);
      }
      break;
    }
    case ModelId::BrokenRigidityKappa: {
      const BrokenRigidityKappaParameters& p =
          configuration.brokenRigidityKappa;
      FingerprintDouble(&canonical, "K_star_m2_per_s", p.kStarM2PerS);
      FingerprintDouble(&canonical, "rigidity0_V", p.rigidity0V);
      FingerprintDouble(&canonical, "break_rigidity_V", p.breakRigidityV);
      FingerprintDouble(&canonical, "low_slope", p.lowSlope);
      FingerprintDouble(&canonical, "high_slope", p.highSlope);
      FingerprintDouble(&canonical, "smoothness", p.smoothness);
      canonical << ";use_field_factor=" << p.useFieldFactor
                << ";use_radial_factor=" << p.useRadialFactor
                << ";use_region_factor=" << p.useRegionFactor;
      if (p.useFieldFactor && p.fieldExponent != 0.0) {
        FingerprintDouble(&canonical, "field_reference_T", p.fieldReferenceT);
        FingerprintDouble(&canonical, "field_exponent", p.fieldExponent);
      }
      break;
    }
    case ModelId::Bohm:
      FingerprintDouble(&canonical, "eta_B", configuration.bohm.etaB);
      canonical << ";field_definition="
                << static_cast<int>(configuration.bohm.fieldDefinition);
      break;
    default:
      break;
  }
  const std::uint64_t hash = Fnv1a(canonical.str());
  std::ostringstream result;
  result << std::hex << std::setw(16) << std::setfill('0') << hash;
  return result.str();
}

ModelFunction FunctionForModel(ModelId model) {
  // Returning nullptr for a registered-but-unimplemented ID is intentional.
  // Evaluate and SetActiveConfiguration convert it to UnsupportedModel rather
  // than selecting a numerically convenient substitute.
  switch (model) {
    case ModelId::ConstantLambda: return &EvaluateConstantLambda;
    case ModelId::ConstantKappa: return &EvaluateConstantKappa;
    case ModelId::PowerLawLambda: return &EvaluatePowerLawLambda;
    case ModelId::BrokenRigidityKappa: return &EvaluateBrokenRigidityKappa;
    case ModelId::Bohm: return &EvaluateBohm;
    default: return nullptr;
  }
}

ModelFunction ActiveModelFunction = &EvaluateUnconfigured;

ParallelResult Evaluate(const ParticleState& particle,
                        const LocalState& local,
                        const ModelConfiguration& configuration) {
  // Direct callers receive the same schema checks as parser-configured calls;
  // constructing ModelConfiguration manually cannot bypass model validation.
  const Status valid = ValidateConfiguration(configuration);
  if (!valid.ok())
    return Failure(valid.code, valid.detail, configuration, local);
  ModelFunction function = FunctionForModel(configuration.model);
  if (!function)
    return Failure(StatusCode::UnsupportedModel,
                   "selected parallel-diffusion model has no evaluator",
                   configuration, local);
  return function(particle, local, configuration);
}

ParallelResult EvaluateActive(const ParticleState& particle,
                              const LocalState& local) {
  // Under the startup-only mutation contract, this reads a stable matching
  // pair: ActiveModelFunction and gActiveConfiguration.  Concurrent calls are
  // safe only after configuration has been frozen.
  return ActiveModelFunction(particle, local, gActiveConfiguration);
}

Status SetActiveConfiguration(const ModelConfiguration& configuration) {
  const Status valid = ValidateConfiguration(configuration);
  if (!valid.ok()) return valid;
  ModelFunction function = FunctionForModel(configuration.model);
  if (!function)
    return Status::Error(StatusCode::UnsupportedModel,
                         "selected model has no implementation function");
  // Assignment order is intentional for the documented serial-initialization
  // contract: the complete validated parameter value is installed before the
  // pointer makes the backend reachable by EvaluateActive.  This ordering is
  // not a synchronization mechanism; concurrent reconfiguration remains
  // outside the API contract.
  gActiveConfiguration = configuration;
  ActiveModelFunction = function;
  return Status::Success();
}

ModelConfiguration GetActiveConfiguration() { return gActiveConfiguration; }

Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          ModelConfiguration* configuration) {
  if (!configuration)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null model-configuration output");
  ModelId model;
  if (!ParseModelId(modelId, &model))
    return Status::Error(StatusCode::UnsupportedModel,
                         "unknown parallel-diffusion model '" + modelId + "'");
  // Parse into a local candidate and publish to *configuration only after all
  // keys, paired-option rules, numerical domains, and backend availability
  // validate.  Thus even a caller that reuses its output object sees
  // transactional behavior on failure.
  ModelConfiguration candidate;
  candidate.model = model;
  std::map<std::string, std::string> values;
  Status status = InputMap(parameters, &values);
  if (!status.ok()) return status;

  if (model == ModelId::ConstantLambda) {
    status = TakeDouble(&values, "lambda_parallel_m", true,
                        &candidate.constantLambda.lambdaParallelM);
  } else if (model == ModelId::ConstantKappa) {
    status = TakeDouble(&values, "kappa_parallel_m2_per_s", true,
                        &candidate.constantKappa.kappaParallelM2PerS);
  } else if (model == ModelId::PowerLawLambda) {
    PowerLawLambdaParameters& p = candidate.powerLawLambda;
    status = TakeDouble(&values, "lambda0_m", true, &p.lambda0M);
    if (!status.ok()) return status;
    const auto variable = values.find("independent_variable");
    if (variable == values.end())
      return Status::Error(StatusCode::MissingInput,
                           "missing parameter 'independent_variable'");
    const std::string name = Lower(variable->second);
    values.erase(variable);
    const char* referenceKey = nullptr;
    if (name == "rigidity") {
      p.independentVariable = IndependentVariable::Rigidity;
      referenceKey = "rigidity0_V";
    } else if (name == "total_kinetic_energy") {
      p.independentVariable = IndependentVariable::TotalKineticEnergy;
      referenceKey = "kinetic_energy0_J";
    } else if (name == "energy_per_nucleon") {
      p.independentVariable = IndependentVariable::EnergyPerNucleon;
      referenceKey = "energy_per_nucleon0_J";
    } else if (name == "speed") {
      p.independentVariable = IndependentVariable::Speed;
      referenceKey = "speed0_m_per_s";
    } else {
      return Status::Error(StatusCode::InvalidConfiguration,
          "independent_variable must be rigidity, total_kinetic_energy, energy_per_nucleon, or speed");
    }
    status = TakeDouble(&values, referenceKey, true,
                        &p.independentReferenceSI);
    if (!status.ok()) return status;
    status = TakeDouble(&values, "independent_exponent", true,
                        &p.independentExponent);
    if (!status.ok()) return status;

    const bool hasRadius0 = values.count("radius0_m") != 0;
    const bool hasRadialExponent = values.count("radial_exponent") != 0;
    // Reference and exponent form one physical factor and must be present or
    // absent together.  Partial specification is more likely a typo than an
    // intended unity factor, so it fails closed.
    if (hasRadius0 != hasRadialExponent)
      return Status::Error(StatusCode::MissingInput,
          "radial factor requires both radius0_m and radial_exponent");
    p.useRadialFactor = hasRadius0;
    if (hasRadius0) {
      status = TakeDouble(&values, "radius0_m", true, &p.radius0M);
      if (!status.ok()) return status;
      status = TakeDouble(&values, "radial_exponent", true,
                          &p.radialExponent);
      if (!status.ok()) return status;
    }
    const bool hasFieldReference = values.count("field_reference_T") != 0;
    const bool hasFieldExponent = values.count("field_exponent") != 0;
    if (hasFieldReference != hasFieldExponent)
      return Status::Error(StatusCode::MissingInput,
          "field factor requires both field_reference_T and field_exponent");
    p.useFieldFactor = hasFieldReference;
    if (hasFieldReference) {
      status = TakeDouble(&values, "field_reference_T", true,
                          &p.fieldReferenceT);
      if (!status.ok()) return status;
      status = TakeDouble(&values, "field_exponent", true, &p.fieldExponent);
      if (!status.ok()) return status;
      // Canonicalize the exact eta=0 identity to a disabled factor.  The
      // supplied B reference then has no effect and no field provider is
      // required at runtime.
      if (p.fieldExponent == 0.0) p.useFieldFactor = false;
    }
    status = TakeBool(&values, "use_time_factor", false, &p.useTimeFactor);
    if (!status.ok()) return status;
    status = TakeBool(&values, "use_region_factor", false,
                      &p.useRegionFactor);
  } else if (model == ModelId::BrokenRigidityKappa) {
    BrokenRigidityKappaParameters& p = candidate.brokenRigidityKappa;
    status = TakeDouble(&values, "K_star_m2_per_s", true, &p.kStarM2PerS);
    if (!status.ok()) return status;
    status = TakeDouble(&values, "rigidity0_V", true, &p.rigidity0V);
    if (!status.ok()) return status;
    status = TakeDouble(&values, "break_rigidity_V", true,
                        &p.breakRigidityV);
    if (!status.ok()) return status;
    status = TakeDouble(&values, "low_slope", true, &p.lowSlope);
    if (!status.ok()) return status;
    status = TakeDouble(&values, "high_slope", true, &p.highSlope);
    if (!status.ok()) return status;
    status = TakeDouble(&values, "smoothness", true, &p.smoothness);
    if (!status.ok()) return status;
    const bool hasFieldReference = values.count("field_reference_T") != 0;
    const bool hasFieldExponent = values.count("field_exponent") != 0;
    if (hasFieldReference != hasFieldExponent)
      return Status::Error(StatusCode::MissingInput,
          "field factor requires both field_reference_T and field_exponent");
    p.useFieldFactor = hasFieldReference;
    if (hasFieldReference) {
      status = TakeDouble(&values, "field_reference_T", true,
                          &p.fieldReferenceT);
      if (!status.ok()) return status;
      status = TakeDouble(&values, "field_exponent", true, &p.fieldExponent);
      if (!status.ok()) return status;
      if (p.fieldExponent == 0.0) p.useFieldFactor = false;
    }
    status = TakeBool(&values, "use_radial_factor", false,
                      &p.useRadialFactor);
    if (!status.ok()) return status;
    status = TakeBool(&values, "use_region_factor", false,
                      &p.useRegionFactor);
  } else if (model == ModelId::Bohm) {
    status = TakeDouble(&values, "eta_B", true, &candidate.bohm.etaB);
    if (!status.ok()) return status;
    const auto definition = values.find("field_definition");
    if (definition == values.end())
      return Status::Error(StatusCode::MissingInput,
                           "missing parameter 'field_definition'");
    const std::string name = Lower(definition->second);
    values.erase(definition);
    if (name == "mean_field")
      candidate.bohm.fieldDefinition = BohmFieldDefinition::MeanField;
    else if (name == "effective_field")
      candidate.bohm.fieldDefinition = BohmFieldDefinition::EffectiveField;
    else
      return Status::Error(StatusCode::InvalidConfiguration,
          "field_definition must be mean_field or effective_field");
  } else {
    return Status::Error(StatusCode::UnsupportedModel,
        std::string("model '") + ModelName(model) +
        "' is reserved by the specification but not implemented");
  }
  if (!status.ok()) return status;
  status = RejectUnknown(values);
  if (!status.ok()) return status;
  status = ValidateConfiguration(candidate);
  if (!status.ok()) return status;
  *configuration = candidate;
  return Status::Success();
}

Status ConfigureActiveModel(const std::string& modelId,
                            const std::vector<InputParameter>& parameters) {
  // Do not call SetActiveConfiguration unless the complete parser-neutral
  // schema succeeds.  This preserves both the old parameter value and old
  // dispatch pointer after any rejected input block.
  ModelConfiguration candidate;
  const Status built = BuildConfiguration(modelId, parameters, &candidate);
  return built.ok() ? SetActiveConfiguration(candidate) : built;
}

Status KStarFromReferenceKappa(double referenceKappaM2PerS,
                               const ParticleState& referenceParticle,
                               double* kStarM2PerS) {
  if (!kStarM2PerS || !FinitePositive(referenceKappaM2PerS))
    return Status::Error(StatusCode::InvalidConfiguration,
        "K_star conversion requires positive reference kappa and non-null output");
  ParticleKinematics kinematics;
  const Status status =
      ComputeParticleKinematics(referenceParticle, &kinematics);
  if (!status.ok()) return status;
  // At R0 and unit field/radial/region factors, Equation (17) is
  // kappa_reference=K_star*beta(referenceParticle).  No species-independent
  // beta can be assumed at a fixed rigidity, so the particle is mandatory.
  *kStarM2PerS = referenceKappaM2PerS / kinematics.beta;
  return FinitePositive(*kStarM2PerS)
      ? Status::Success()
      : Status::Error(StatusCode::OutsideModelDomain,
                      "K_star conversion overflowed");
}

}  // namespace ParallelDiffusion
}  // namespace SEP

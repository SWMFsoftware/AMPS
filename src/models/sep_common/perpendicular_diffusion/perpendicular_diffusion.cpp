#include "perpendicular_diffusion.h"
#include "perpendicular_diffusion_internal.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

namespace SEP {
namespace PerpendicularDiffusion {
namespace {

ModelConfiguration ActiveConfiguration;

bool Finite(double value) { return std::isfinite(value); }
bool Positive(double value) { return Finite(value) && value > 0.0; }

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
      [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
      });
  return value;
}

bool ParseDouble(const std::string& text, double* value) {
  if (!value || text.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || !end || *end != '\0' || !Finite(parsed)) return false;
  *value = parsed;
  return true;
}

Status ParseList(const std::string& text, std::vector<double>* values) {
  if (!values)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null list output");
  std::vector<double> parsed;
  std::istringstream stream(text);
  std::string item;
  while (std::getline(stream, item, ',')) {
    double value = 0.0;
    if (!ParseDouble(item, &value))
      return Status::Error(StatusCode::InvalidConfiguration,
                           "invalid comma-separated numeric list");
    parsed.push_back(value);
  }
  if (parsed.empty())
    return Status::Error(StatusCode::InvalidConfiguration,
                         "numeric list must not be empty");
  *values = std::move(parsed);
  return Status::Success();
}

struct Schema {
  std::vector<std::string> requiredNumeric;
  std::vector<std::string> optionalNumeric;
  std::vector<std::string> requiredText;
  std::vector<std::string> optionalText;
};

Schema NumericalSchema(std::vector<std::string> required = {},
                       std::vector<std::string> optional = {}) {
  optional.push_back("relative_tolerance");
  optional.push_back("maximum_refinements");
  optional.push_back("maximum_iterations");
  return Schema{std::move(required), std::move(optional), {}, {}};
}

Schema SchemaFor(ModelId model) {
  switch (model) {
    case ModelId::ConstantKappaPerp:
      return Schema{{"kappa_perp_m2_per_s"}, {}, {}, {}};
    case ModelId::ConstantLambdaPerp:
      return Schema{{"lambda_perp_m"}, {}, {}, {}};
    case ModelId::RatioKappa:
      return Schema{{"eta_kappa"}, {}, {}, {}};
    case ModelId::PowerLawPerp:
      return Schema{{"kappa0_m2_per_s", "beta_exponent", "rigidity0_V",
                     "rigidity_exponent", "field_reference_T",
                     "field_exponent", "radius0_m", "radius_exponent"},
                    {}, {}, {}};
    case ModelId::PitchAnglePerp:
      return Schema{{"D0_m2_per_s"}, {}, {"shape"}, {}};
    case ModelId::DrogeLambdaPerp:
      return Schema{{"alpha_D", "radius0_m"}, {}, {}, {}};
    case ModelId::ParadiseAlpha:
      return Schema{{"alpha_P", "field_reference_T"}, {}, {}, {}};
    case ModelId::FieldLineSlab:
    case ModelId::FieldLine2D:
    case ModelId::FieldLineComposite:
    case ModelId::NlgceN:
    case ModelId::NlgceF2014:
    case ModelId::ClassicalScattering:
    case ModelId::DriftWeakScattering:
    case ModelId::DriftClassical:
      return Schema{};
    case ModelId::FlrwParticle:
      return Schema{{"a_FL"}, {}, {"pitch_angle_mode"}, {}};
    case ModelId::Nlgc:
    case ModelId::NlgcSlabKernelDiagnostic:
      return NumericalSchema({"a_squared"});
    case ModelId::Enlgc2D:
    case ModelId::ImplicitSlabExact2016:
    case ModelId::ImplicitSlabRational2016:
    case ModelId::FlpdComplete:
      return NumericalSchema();
    case ModelId::Unlt:
      return NumericalSchema({"a_squared"});
    case ModelId::RbdBc:
      return NumericalSchema({"a_squared", "area_energy_index",
                              "area_inertial_index"});
    case ModelId::CompositeClosed2019:
      return Schema{{"perpendicular_length_m"}, {}, {"length_profile"}, {}};
    case ModelId::CompoundDiffusiveLines:
    case ModelId::GcdCompound:
    case ModelId::EnlgcSlabSecantDiagnostic:
      return Schema{{"age_s"}, {}, {}, {}};
    case ModelId::PrediffusiveFit:
      return Schema{{"age_s", "A1_m2", "gyroperiod_s", "t1_s", "t2_s",
                     "alpha_m", "beta_m"}, {}, {"calibration_id"}, {}};
    case ModelId::IsoFitCandiaRoulet2004:
      return Schema{{"outer_scale_m"}, {}, {"spectrum", "domain_policy"}, {}};
    case ModelId::IsoFitSnodin2016:
      return Schema{{"outer_scale_m"}, {}, {"domain_policy"}, {}};
    case ModelId::IsoKappa0Snodin2016:
      return Schema{{"outer_scale_m"}, {}, {"spectrum", "domain_policy"}, {}};
    case ModelId::IsoFitKuhlen2025:
      return Schema{{"A", "rho_star", "s_kappa", "C_K", "z1_m", "z2_m",
                     "gamma_K", "transverse_correlation_length_m"},
                    {"age_s", "relative_tolerance", "maximum_refinements",
                     "maximum_iterations"},
                    {"calibration_id", "root_selection"}, {}};
    case ModelId::YanLazarianMa4:
      return Schema{{"c_M4", "alfven_mach", "injection_scale_m"}, {}, {}, {}};
    case ModelId::NwuRatioPolar:
      return Schema{{"eta_r", "eta_theta", "polar_enhancement", "theta_F_rad",
                     "width_per_rad"}, {}, {"angular_convention"}, {}};
    case ModelId::CortiAms02:
      return Schema{{"K0_m2_per_s", "field_reference_T", "break_rigidity_V",
                     "smoothness", "parallel_low_slope", "parallel_high_slope",
                     "perp_low_slope", "perp_high_slope", "width_per_rad"},
                    {}, {"angular_convention", "calibration_id"}, {}};
    case ModelId::HelmodRatio:
      return Schema{{"rho_H"}, {}, {"polar_theta_rad", "polar_factor",
                                    "polar_function_identity",
                                    "polar_interpolation"}, {}};
    case ModelId::DriftRigidityReduction:
      return Schema{{"K_A0", "rigidity_A_V"}, {}, {}, {}};
    case ModelId::DriftCandiaRoulet2004:
      return Schema{{"outer_scale_m", "rho_min", "rho_max", "sigma2_min",
                     "sigma2_max"}, {}, {"spectrum", "domain_policy",
                                         "domain_identity"}, {}};
    case ModelId::TabulatedPerp:
      return Schema{{}, {}, {"axis", "axis_values_SI", "values_m2_per_s",
                             "generation_identity", "table_checksum",
                             "boundary_policy", "interpolation",
                             "zero_policy"}, {}};
    case ModelId::UnltFgrPerturbativeDiagnostic:
      return Schema{{"maximum_relative_correction"}, {}, {}, {}};
    case ModelId::FrozenFieldline:
    case ModelId::QltPerp:
    case ModelId::UnltFgr:
    case ModelId::IsoFitCasse2002:
    case ModelId::RestrictedScattering:
      return Schema{};
  }
  return Schema{};
}

bool Contains(const std::vector<std::string>& values, const std::string& key) {
  return std::find(values.begin(), values.end(), key) != values.end();
}

bool IsSourceGated(ModelId model) {
  return model == ModelId::FrozenFieldline || model == ModelId::QltPerp ||
      model == ModelId::UnltFgr || model == ModelId::IsoFitCasse2002 ||
      model == ModelId::RestrictedScattering;
}

bool Nonnegative(const ModelConfiguration& configuration,
                 const std::string& key) {
  const auto found = configuration.numeric.find(key);
  return found != configuration.numeric.end() && Finite(found->second) &&
         found->second >= 0.0;
}

bool PositiveParameter(const ModelConfiguration& configuration,
                       const std::string& key) {
  const auto found = configuration.numeric.find(key);
  return found != configuration.numeric.end() && Positive(found->second);
}

std::string CanonicalConfiguration(const ModelConfiguration& configuration) {
  std::ostringstream output;
  output << "spec=2.1;model=" << ModelName(configuration.model);
  output << std::setprecision(17);
  for (const auto& value : configuration.numeric)
    output << ';' << value.first << '=' << value.second;
  for (const auto& value : configuration.text)
    output << ';' << value.first << '=' << value.second;
  for (double value : configuration.tableAxisSI) output << ";axis=" << value;
  for (double value : configuration.tableValuesSI) output << ";value=" << value;
  return output.str();
}

ModelResult ConfigurationFailure(const ModelConfiguration& configuration,
                                 StatusCode code,
                                 const std::string& detail) {
  ModelResult result;
  result.status=Status::Error(code,detail);
  result.provenance.requestedModel=ModelName(configuration.model);
  result.provenance.usedModel=ModelName(configuration.model);
  result.provenance.equationVersion="perpendicular-spec-2.1";
  result.provenance.configurationFingerprint=
      ConfigurationFingerprint(configuration);
  return result;
}

template<ModelId Selected>
ModelResult EvaluateSelected(const ParticleState& particle,
                             const LocalState& local,
                             const ModelConfiguration& configuration) {
  // The active pointer is public, so defend against a caller pairing a pointer
  // obtained for one model with another model's configuration. Normal manager
  // use installs an already validated matching pair during serial startup.
  if(configuration.model!=Selected)
    return ConfigurationFailure(configuration,StatusCode::InvalidConfiguration,
                                "model function/configuration identity mismatch");
  const Status valid=ValidateConfiguration(configuration);
  if(!valid.ok()) return ConfigurationFailure(configuration,valid.code,valid.detail);
  return Internal::EvaluateModel(particle,local,configuration);
}

ModelResult EvaluateUnconfigured(const ParticleState&,
                                 const LocalState&,
                                 const ModelConfiguration& configuration) {
  return ConfigurationFailure(configuration,StatusCode::InvalidConfiguration,
                              "no active perpendicular model is configured");
}

}  // namespace

Status Status::Success() { return Status(); }

Status Status::Error(StatusCode code, const std::string& detail) {
  Status result;
  result.code = code;
  result.detail = detail;
  return result;
}

const std::vector<ModelDescriptor>& ModelRegistry() {
  static const std::vector<ModelDescriptor> models = {
      {ModelId::ConstantKappaPerp,"constant_kappa_perp",true,true,Observable::SymmetricCoefficient,"P1"},
      {ModelId::ConstantLambdaPerp,"constant_lambda_perp",true,true,Observable::SymmetricCoefficient,"P1"},
      {ModelId::RatioKappa,"ratio_kappa",true,true,Observable::SymmetricCoefficient,"P1"},
      {ModelId::PowerLawPerp,"power_law_perp",true,true,Observable::SymmetricCoefficient,"P2"},
      {ModelId::PitchAnglePerp,"pitch_angle_perp",true,true,Observable::PitchAngleCoefficient,"P3"},
      {ModelId::DrogeLambdaPerp,"droge_lambda_perp",true,true,Observable::PitchAngleCoefficient,"P5-P6"},
      {ModelId::ParadiseAlpha,"paradise_alpha",true,true,Observable::SymmetricCoefficient,"P7"},
      {ModelId::FieldLineSlab,"fl_slab",true,false,Observable::FieldLineCoefficient,"F2"},
      {ModelId::FieldLine2D,"fl_2d",true,false,Observable::FieldLineCoefficient,"F3"},
      {ModelId::FieldLineComposite,"fl_composite",true,false,Observable::FieldLineCoefficient,"F4"},
      {ModelId::FlrwParticle,"flrw_particle",true,true,Observable::PitchAngleCoefficient,"F5/P4"},
      {ModelId::Nlgc,"nlgc",true,true,Observable::SymmetricCoefficient,"N1-N4"},
      {ModelId::Enlgc2D,"enlgc_2d",true,true,Observable::SymmetricCoefficient,"N5"},
      {ModelId::NlgceN,"nlgce_n",true,true,Observable::PairedCoefficients,"E1-E6"},
      {ModelId::NlgceF2014,"nlgce_f_2014",true,true,Observable::PairedCoefficients,"E7-E9"},
      {ModelId::Unlt,"unlt",true,true,Observable::SymmetricCoefficient,"U1-U2"},
      {ModelId::ImplicitSlabExact2016,"implicit_slab_exact_2016",true,true,Observable::SymmetricCoefficient,"U5-U7"},
      {ModelId::ImplicitSlabRational2016,"implicit_slab_rational_2016",true,true,Observable::SymmetricCoefficient,"U8-U9"},
      {ModelId::FlpdComplete,"flpd_complete",true,true,Observable::SymmetricCoefficient,"D1-D2"},
      {ModelId::RbdBc,"rbd_bc",true,true,Observable::SymmetricCoefficient,"B1-B4/B9"},
      {ModelId::CompositeClosed2019,"composite_closed_2019",true,true,Observable::SymmetricCoefficient,"B5-B8"},
      {ModelId::CompoundDiffusiveLines,"compound_diffusive_lines",true,false,Observable::ParticleMsd,"F8"},
      {ModelId::GcdCompound,"gcd_compound",true,false,Observable::ParticleMsd,"C3-C6"},
      {ModelId::FrozenFieldline,"frozen_fieldline",false,false,Observable::FieldLineCoefficient,"C7-C9"},
      {ModelId::PrediffusiveFit,"prediffusive_fit",true,false,Observable::ConditionalReturnMedian,"C10-C11"},
      {ModelId::IsoFitCandiaRoulet2004,"iso_fit_candia_roulet_2004",true,true,Observable::PairedCoefficients,"I4-I5"},
      {ModelId::IsoFitSnodin2016,"iso_fit_snodin_2016",true,true,Observable::PairedCoefficients,"I6-I7"},
      {ModelId::IsoKappa0Snodin2016,"iso_kappa0_snodin_2016",true,true,Observable::SymmetricCoefficient,"I6"},
      {ModelId::IsoFitKuhlen2025,"iso_fit_kuhlen_2025",true,true,Observable::PairedCoefficients,"I8-I13"},
      {ModelId::ClassicalScattering,"classical_scattering",true,true,Observable::PairedCoefficients,"I14"},
      {ModelId::YanLazarianMa4,"yan_lazarian_ma4",true,true,Observable::SymmetricCoefficient,"I16"},
      {ModelId::NwuRatioPolar,"nwu_ratio_polar",true,true,Observable::PairedCoefficients,"G3-G4"},
      {ModelId::CortiAms02,"corti_ams02",true,true,Observable::PairedCoefficients,"G2/G5-G6"},
      {ModelId::HelmodRatio,"helmod_ratio",true,true,Observable::PairedCoefficients,"G7"},
      {ModelId::DriftWeakScattering,"drift_weak_scattering",true,true,Observable::SignedHallCoefficient,"G8"},
      {ModelId::DriftRigidityReduction,"drift_rigidity_reduction",true,true,Observable::SignedHallCoefficient,"G9"},
      {ModelId::DriftClassical,"drift_classical",true,true,Observable::SignedHallCoefficient,"I14"},
      {ModelId::DriftCandiaRoulet2004,"drift_candia_roulet_2004",true,true,Observable::SignedHallCoefficient,"G10"},
      {ModelId::TabulatedPerp,"tabulated_perp",true,true,Observable::SymmetricCoefficient,"Section 16.3"},
      {ModelId::QltPerp,"qlt_perp",false,false,Observable::PerturbativeDiagnostic,"source gate"},
      {ModelId::UnltFgrPerturbativeDiagnostic,"unlt_fgr_perturbative_diagnostic",true,false,Observable::PerturbativeDiagnostic,"U10"},
      {ModelId::UnltFgr,"unlt_fgr",false,false,Observable::PerturbativeDiagnostic,"source gate"},
      {ModelId::EnlgcSlabSecantDiagnostic,"enlgc_slab_secant_diagnostic",true,false,Observable::ParticleMsd,"N6"},
      {ModelId::NlgcSlabKernelDiagnostic,"nlgc_slab_kernel_diagnostic",true,false,Observable::PerturbativeDiagnostic,"N2/N8"},
      {ModelId::IsoFitCasse2002,"iso_fit_casse_2002",false,false,Observable::SymmetricCoefficient,"source gate"},
      {ModelId::RestrictedScattering,"restricted_scattering",false,false,Observable::SymmetricCoefficient,"source gate"}};
  return models;
}

const char* ModelName(ModelId model) {
  for (const ModelDescriptor& descriptor : ModelRegistry())
    if (descriptor.id == model) return descriptor.stableId;
  return "unknown";
}

bool ParseModelId(const std::string& text, ModelId* model) {
  if (!model) return false;
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
  if (!kinematics || !Positive(particle.massKg) ||
      !Finite(particle.chargeC) || particle.chargeC == 0.0 ||
      !Positive(particle.momentumKgMPerS))
    return Status::Error(StatusCode::InvalidInput,
                         "particle requires positive mass/momentum and nonzero finite charge");
  const double mc = particle.massKg * SpeedOfLightMPerS;
  const double ratio = particle.momentumKgMPerS / mc;
  ParticleKinematics result;
  result.gamma = std::hypot(1.0, ratio);
  result.beta = ratio / result.gamma;
  result.speedMPerS = result.beta * SpeedOfLightMPerS;
  result.rigidityV = particle.momentumKgMPerS * SpeedOfLightMPerS /
      std::fabs(particle.chargeC);
  if (!Positive(result.speedMPerS) || !Positive(result.rigidityV))
    return Status::Error(StatusCode::OutsideModelDomain,
                         "particle kinematics are not representable");
  *kinematics = result;
  return Status::Success();
}

Status BuildConfiguration(const std::string& modelId,
                          const std::vector<InputParameter>& parameters,
                          ModelConfiguration* configuration) {
  if (!configuration)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null configuration output");
  ModelId model;
  if (!ParseModelId(modelId, &model))
    return Status::Error(StatusCode::UnsupportedModel,
                         "unknown perpendicular model '" + modelId + "'");
  if (IsSourceGated(model))
    return Status::Error(StatusCode::SourceGate,
        "model '" + std::string(ModelName(model)) +
        "' remains unavailable under specification Section 20.4");

  const Schema schema = SchemaFor(model);
  ModelConfiguration candidate;
  candidate.model = model;
  std::set<std::string> seen;
  for (const InputParameter& parameter : parameters) {
    if (!seen.insert(parameter.name).second)
      return Status::Error(StatusCode::InvalidConfiguration,
                           "duplicate parameter '" + parameter.name + "'");
    if (Contains(schema.requiredNumeric, parameter.name) ||
        Contains(schema.optionalNumeric, parameter.name)) {
      double value = 0.0;
      if (!ParseDouble(parameter.value, &value))
        return Status::Error(StatusCode::InvalidConfiguration,
                             "parameter '" + parameter.name + "' must be a finite SI number");
      candidate.numeric.emplace(parameter.name, value);
    } else if (Contains(schema.requiredText, parameter.name) ||
               Contains(schema.optionalText, parameter.name)) {
      candidate.text.emplace(parameter.name, parameter.value);
    } else {
      return Status::Error(StatusCode::InvalidConfiguration,
                           "unknown parameter '" + parameter.name + "'");
    }
  }
  for (const std::string& key : schema.requiredNumeric)
    if (!candidate.numeric.count(key))
      return Status::Error(model==ModelId::IsoFitKuhlen2025 &&
                               (key=="gamma_K" ||
                                key=="transverse_correlation_length_m") ?
                           StatusCode::MissingCalibration : StatusCode::MissingInput,
                           "missing parameter '" + key + "'");
  for (const std::string& key : schema.requiredText)
    if (!candidate.text.count(key))
      return Status::Error(model==ModelId::IsoFitKuhlen2025 ?
                           StatusCode::MissingCalibration : StatusCode::MissingInput,
                           "missing parameter '" + key + "'");

  // Tables and polar functions are parsed into immutable numeric arrays; the
  // original text identity/checksum stays in the configuration provenance.
  if (model == ModelId::TabulatedPerp) {
    Status status = ParseList(candidate.text["axis_values_SI"],
                              &candidate.tableAxisSI);
    if (!status.ok()) return status;
    status = ParseList(candidate.text["values_m2_per_s"],
                       &candidate.tableValuesSI);
    if (!status.ok()) return status;
  } else if (model == ModelId::HelmodRatio) {
    Status status = ParseList(candidate.text["polar_theta_rad"],
                              &candidate.tableAxisSI);
    if (!status.ok()) return status;
    status = ParseList(candidate.text["polar_factor"],
                       &candidate.tableValuesSI);
    if (!status.ok()) return status;
  }
  const Status valid = ValidateConfiguration(candidate);
  if (!valid.ok()) return valid;
  *configuration = std::move(candidate);
  return Status::Success();
}

Status ValidateConfiguration(const ModelConfiguration& c) {
  if (IsSourceGated(c.model))
    return Status::Error(StatusCode::SourceGate,
                         "selected model remains source gated");
  // Direct callers may construct ModelConfiguration without using the parser
  // manager. Reapply the complete selected-model schema here before any
  // map::at access: malformed direct values must return a typed status rather
  // than throw, and inactive keys must not survive a model change.
  const Schema schema=SchemaFor(c.model);
  for(const std::string& key:schema.requiredNumeric)
    if(!c.numeric.count(key))
      return Status::Error(c.model==ModelId::IsoFitKuhlen2025 &&
                               (key=="gamma_K" ||
                                key=="transverse_correlation_length_m") ?
                           StatusCode::MissingCalibration : StatusCode::MissingInput,
                           "missing parameter '"+key+"'");
  for(const std::string& key:schema.requiredText)
    if(!c.text.count(key))
      return Status::Error(c.model==ModelId::IsoFitKuhlen2025 ?
                           StatusCode::MissingCalibration : StatusCode::MissingInput,
                           "missing parameter '"+key+"'");
  for(const auto& value:c.numeric)
    if(!Contains(schema.requiredNumeric,value.first) &&
       !Contains(schema.optionalNumeric,value.first))
      return Status::Error(StatusCode::InvalidConfiguration,
                           "inactive numeric parameter '"+value.first+"'");
  for(const auto& value:c.text)
    if(!Contains(schema.requiredText,value.first) &&
       !Contains(schema.optionalText,value.first))
      return Status::Error(StatusCode::InvalidConfiguration,
                           "inactive text parameter '"+value.first+"'");
  auto finiteAll = [&]() {
    for (const auto& value : c.numeric)
      if (!Finite(value.second)) return false;
    return true;
  };
  if (!finiteAll())
    return Status::Error(StatusCode::InvalidConfiguration,
                         "configuration contains a nonfinite number");
  if (c.numeric.count("relative_tolerance") &&
      (!PositiveParameter(c, "relative_tolerance") ||
       c.numeric.at("relative_tolerance") > 1.0e-2))
    return Status::Error(StatusCode::InvalidConfiguration,
                         "relative_tolerance must be in (0,1e-2]");
  for (const char* key : {"maximum_refinements", "maximum_iterations"})
    if (c.numeric.count(key) &&
        (!(c.numeric.at(key) >= 1.0) ||
         c.numeric.at(key) != std::floor(c.numeric.at(key))))
      return Status::Error(StatusCode::InvalidConfiguration,
                           std::string(key) + " must be a positive integer");

  switch (c.model) {
    case ModelId::ConstantKappaPerp:
      if (Nonnegative(c,"kappa_perp_m2_per_s")) return Status::Success();
      break;
    case ModelId::ConstantLambdaPerp:
      if (Nonnegative(c,"lambda_perp_m")) return Status::Success();
      break;
    case ModelId::RatioKappa:
      if (Nonnegative(c,"eta_kappa")) return Status::Success();
      break;
    case ModelId::PowerLawPerp:
      if (Nonnegative(c,"kappa0_m2_per_s") && PositiveParameter(c,"rigidity0_V") &&
          PositiveParameter(c,"field_reference_T") && PositiveParameter(c,"radius0_m"))
        return Status::Success();
      break;
    case ModelId::PitchAnglePerp: {
      const std::string& shape = c.text.at("shape");
      if (Nonnegative(c,"D0_m2_per_s") &&
          (shape=="isotropic" || shape=="abs_mu" || shape=="sqrt_one_minus_mu2"))
        return Status::Success();
      break;
    }
    case ModelId::FlrwParticle:
      if (Nonnegative(c,"a_FL") &&
          (c.text.at("pitch_angle_mode")=="isotropic_average" ||
           c.text.at("pitch_angle_mode")=="local_mu"))
        return Status::Success();
      break;
    case ModelId::DrogeLambdaPerp:
      if (Nonnegative(c,"alpha_D") && PositiveParameter(c,"radius0_m")) return Status::Success();
      break;
    case ModelId::ParadiseAlpha:
      if (Nonnegative(c,"alpha_P") && PositiveParameter(c,"field_reference_T")) return Status::Success();
      break;
    case ModelId::Nlgc:
    case ModelId::NlgcSlabKernelDiagnostic:
    case ModelId::Unlt:
      if (Nonnegative(c,"a_squared")) return Status::Success();
      break;
    case ModelId::RbdBc:
      if (Nonnegative(c,"a_squared") &&
          c.numeric.at("area_energy_index") > -2.0 &&
          c.numeric.at("area_inertial_index") > 1.0)
        return Status::Success();
      break;
    case ModelId::CompositeClosed2019:
      if (PositiveParameter(c,"perpendicular_length_m") &&
          (c.text.at("length_profile")=="shalchi_2019_bend_over" ||
           c.text.at("length_profile")=="snodin_2022_integral" ||
           c.text.at("length_profile")=="parameterized_tagged_length"))
        return Status::Success();
      break;
    case ModelId::CompoundDiffusiveLines:
    case ModelId::GcdCompound:
    case ModelId::EnlgcSlabSecantDiagnostic:
      if (PositiveParameter(c,"age_s")) return Status::Success();
      break;
    case ModelId::PrediffusiveFit:
      if (Nonnegative(c,"A1_m2") && PositiveParameter(c,"age_s") &&
          PositiveParameter(c,"gyroperiod_s") && PositiveParameter(c,"t1_s") &&
          PositiveParameter(c,"t2_s") && c.numeric.at("t1_s")<c.numeric.at("t2_s") &&
          c.numeric.at("beta_m")>c.numeric.at("alpha_m") &&
          c.numeric.at("beta_m")>1.0 && !c.text.at("calibration_id").empty())
        return Status::Success();
      break;
    case ModelId::YanLazarianMa4:
      if (Nonnegative(c,"c_M4") && PositiveParameter(c,"alfven_mach") &&
          c.numeric.at("alfven_mach")<1.0 && PositiveParameter(c,"injection_scale_m"))
        return Status::Success();
      break;
    case ModelId::IsoFitCandiaRoulet2004:
      if (PositiveParameter(c,"outer_scale_m") &&
          (c.text.at("spectrum")=="kraichnan" ||
           c.text.at("spectrum")=="kolmogorov" ||
           c.text.at("spectrum")=="bykov_toptygin") &&
          (c.text.at("domain_policy")=="reject" ||
           c.text.at("domain_policy")=="tag_extrapolated"))
        return Status::Success();
      break;
    case ModelId::IsoFitSnodin2016:
      if (PositiveParameter(c,"outer_scale_m") &&
          (c.text.at("domain_policy")=="reject" ||
           c.text.at("domain_policy")=="tag_extrapolated"))
        return Status::Success();
      break;
    case ModelId::IsoKappa0Snodin2016:
      if (PositiveParameter(c,"outer_scale_m") &&
          (c.text.at("spectrum")=="sharp_kolmogorov" ||
           c.text.at("spectrum")=="sharp_kraichnan" ||
           c.text.at("spectrum")=="peaked_kolmogorov") &&
          (c.text.at("domain_policy")=="reject" ||
           c.text.at("domain_policy")=="tag_extrapolated"))
        return Status::Success();
      break;
    case ModelId::IsoFitKuhlen2025:
      if (PositiveParameter(c,"A") && PositiveParameter(c,"rho_star") &&
          PositiveParameter(c,"s_kappa") && PositiveParameter(c,"C_K") &&
          PositiveParameter(c,"z1_m") && PositiveParameter(c,"z2_m") &&
          c.numeric.at("gamma_K") < 0.0 &&
          PositiveParameter(c,"transverse_correlation_length_m") &&
          !c.text.at("calibration_id").empty() &&
          c.text.at("root_selection") == "first_upward_crossing" &&
          (!c.numeric.count("age_s") || Nonnegative(c,"age_s")))
        return Status::Success();
      break;
    case ModelId::DriftRigidityReduction:
      if (Nonnegative(c,"K_A0") && c.numeric.at("K_A0")<=1.0 &&
          PositiveParameter(c,"rigidity_A_V")) return Status::Success();
      break;
    case ModelId::NwuRatioPolar:
      if (Nonnegative(c,"eta_r") && Nonnegative(c,"eta_theta") &&
          PositiveParameter(c,"polar_enhancement") &&
          c.numeric.at("theta_F_rad")>=0.0 &&
          PositiveParameter(c,"width_per_rad") &&
          c.text.at("angular_convention")=="folded_radian")
        return Status::Success();
      break;
    case ModelId::CortiAms02:
      if (Nonnegative(c,"K0_m2_per_s") && PositiveParameter(c,"field_reference_T") &&
          PositiveParameter(c,"break_rigidity_V") && PositiveParameter(c,"smoothness") &&
          PositiveParameter(c,"width_per_rad") &&
          c.text.at("angular_convention")=="parameterized_folded_radian" &&
          !c.text.at("calibration_id").empty())
        return Status::Success();
      break;
    case ModelId::DriftCandiaRoulet2004:
      if (PositiveParameter(c,"outer_scale_m") &&
          PositiveParameter(c,"rho_min") &&
          c.numeric.at("rho_max")>=c.numeric.at("rho_min") &&
          Nonnegative(c,"sigma2_min") &&
          c.numeric.at("sigma2_max")>=c.numeric.at("sigma2_min") &&
          (c.text.at("spectrum")=="kraichnan" ||
           c.text.at("spectrum")=="kolmogorov" ||
           c.text.at("spectrum")=="bykov_toptygin") &&
          (c.text.at("domain_policy")=="reject" ||
           c.text.at("domain_policy")=="tag_extrapolated") &&
          !c.text.at("domain_identity").empty())
        return Status::Success();
      break;
    case ModelId::UnltFgrPerturbativeDiagnostic:
      if (PositiveParameter(c,"maximum_relative_correction") &&
          c.numeric.at("maximum_relative_correction")<1.0) return Status::Success();
      break;
    case ModelId::TabulatedPerp:
    case ModelId::HelmodRatio:
      if (c.tableAxisSI.size()>=2 && c.tableAxisSI.size()==c.tableValuesSI.size()) {
        for (std::size_t i=0;i<c.tableAxisSI.size();++i)
          if (!Finite(c.tableAxisSI[i]) || !Finite(c.tableValuesSI[i]) ||
              c.tableValuesSI[i]<0.0 || (i && !(c.tableAxisSI[i]>c.tableAxisSI[i-1])))
            return Status::Error(StatusCode::InvalidConfiguration,
                                 "table axes must increase and values be nonnegative");
        if (c.model==ModelId::HelmodRatio &&
            (!Nonnegative(c,"rho_H") ||
             c.text.at("polar_function_identity").empty() ||
             c.text.at("polar_interpolation")!="linear")) break;
        if (c.model==ModelId::TabulatedPerp &&
            !((c.text.at("axis")=="rigidity" || c.text.at("axis")=="time" ||
               c.text.at("axis")=="heliocentric_radius" ||
               c.text.at("axis")=="mean_field_magnitude") &&
              c.text.at("boundary_policy")=="reject" &&
              ((c.text.at("interpolation")=="linear" &&
                c.text.at("zero_policy")=="linear_explicit_zero") ||
               (c.text.at("interpolation")=="log_log" &&
                c.text.at("zero_policy")=="strictly_positive")) &&
              !c.text.at("generation_identity").empty() &&
              !c.text.at("table_checksum").empty())) break;
        if(c.model==ModelId::TabulatedPerp &&
           c.text.at("interpolation")=="log_log") {
          for(std::size_t i=0;i<c.tableAxisSI.size();++i)
            if(!Positive(c.tableAxisSI[i])||!Positive(c.tableValuesSI[i]))
              return Status::Error(StatusCode::InvalidConfiguration,
                  "log_log table interpolation requires strictly positive axes and values");
        }
        return Status::Success();
      }
      break;
    default:
      return Status::Success();
  }
  return Status::Error(StatusCode::InvalidConfiguration,
                       "invalid configuration for model '"+
                       std::string(ModelName(c.model))+"'");
}

std::string ConfigurationFingerprint(const ModelConfiguration& configuration) {
  // FNV-1a is used only as a reproducible configuration identity, not as a
  // security primitive or an asset checksum. Scientific data retain their
  // separately verified SHA-256 identities.
  std::uint64_t hash = 1469598103934665603ull;
  for (unsigned char byte : CanonicalConfiguration(configuration)) {
    hash ^= byte;
    hash *= 1099511628211ull;
  }
  std::ostringstream output;
  output << "fnv1a64:" << std::hex << std::setfill('0') << std::setw(16) << hash;
  return output.str();
}

ModelFunction FunctionForModel(ModelId model) {
#define PERP_FUNCTION_CASE(name) case ModelId::name: return &EvaluateSelected<ModelId::name>
  switch(model) {
    PERP_FUNCTION_CASE(ConstantKappaPerp);
    PERP_FUNCTION_CASE(ConstantLambdaPerp);
    PERP_FUNCTION_CASE(RatioKappa);
    PERP_FUNCTION_CASE(PowerLawPerp);
    PERP_FUNCTION_CASE(PitchAnglePerp);
    PERP_FUNCTION_CASE(DrogeLambdaPerp);
    PERP_FUNCTION_CASE(ParadiseAlpha);
    PERP_FUNCTION_CASE(FieldLineSlab);
    PERP_FUNCTION_CASE(FieldLine2D);
    PERP_FUNCTION_CASE(FieldLineComposite);
    PERP_FUNCTION_CASE(FlrwParticle);
    PERP_FUNCTION_CASE(Nlgc);
    PERP_FUNCTION_CASE(Enlgc2D);
    PERP_FUNCTION_CASE(NlgceN);
    PERP_FUNCTION_CASE(NlgceF2014);
    PERP_FUNCTION_CASE(Unlt);
    PERP_FUNCTION_CASE(ImplicitSlabExact2016);
    PERP_FUNCTION_CASE(ImplicitSlabRational2016);
    PERP_FUNCTION_CASE(FlpdComplete);
    PERP_FUNCTION_CASE(RbdBc);
    PERP_FUNCTION_CASE(CompositeClosed2019);
    PERP_FUNCTION_CASE(CompoundDiffusiveLines);
    PERP_FUNCTION_CASE(GcdCompound);
    PERP_FUNCTION_CASE(FrozenFieldline);
    PERP_FUNCTION_CASE(PrediffusiveFit);
    PERP_FUNCTION_CASE(IsoFitCandiaRoulet2004);
    PERP_FUNCTION_CASE(IsoFitSnodin2016);
    PERP_FUNCTION_CASE(IsoKappa0Snodin2016);
    PERP_FUNCTION_CASE(IsoFitKuhlen2025);
    PERP_FUNCTION_CASE(ClassicalScattering);
    PERP_FUNCTION_CASE(YanLazarianMa4);
    PERP_FUNCTION_CASE(NwuRatioPolar);
    PERP_FUNCTION_CASE(CortiAms02);
    PERP_FUNCTION_CASE(HelmodRatio);
    PERP_FUNCTION_CASE(DriftWeakScattering);
    PERP_FUNCTION_CASE(DriftRigidityReduction);
    PERP_FUNCTION_CASE(DriftClassical);
    PERP_FUNCTION_CASE(DriftCandiaRoulet2004);
    PERP_FUNCTION_CASE(TabulatedPerp);
    PERP_FUNCTION_CASE(QltPerp);
    PERP_FUNCTION_CASE(UnltFgrPerturbativeDiagnostic);
    PERP_FUNCTION_CASE(UnltFgr);
    PERP_FUNCTION_CASE(EnlgcSlabSecantDiagnostic);
    PERP_FUNCTION_CASE(NlgcSlabKernelDiagnostic);
    PERP_FUNCTION_CASE(IsoFitCasse2002);
    PERP_FUNCTION_CASE(RestrictedScattering);
  }
#undef PERP_FUNCTION_CASE
  return nullptr;
}

ModelResult Evaluate(const ParticleState& particle, const LocalState& local,
                     const ModelConfiguration& configuration) {
  ModelFunction function=FunctionForModel(configuration.model);
  return function ? function(particle,local,configuration) :
      ConfigurationFailure(configuration,StatusCode::UnsupportedModel,
                           "invalid perpendicular model identifier");
}

ModelFunction ActiveModelFunction = &EvaluateUnconfigured;

Status SetActiveConfiguration(const ModelConfiguration& configuration) {
  const Status valid = ValidateConfiguration(configuration);
  if (!valid.ok()) return valid;
  ActiveConfiguration = configuration;
  ActiveModelFunction = FunctionForModel(configuration.model);
  return Status::Success();
}

Status ConfigureActiveModel(const std::string& modelId,
                            const std::vector<InputParameter>& parameters) {
  ModelConfiguration candidate;
  const Status built = BuildConfiguration(modelId, parameters, &candidate);
  return built.ok() ? SetActiveConfiguration(candidate) : built;
}

ModelConfiguration GetActiveConfiguration() { return ActiveConfiguration; }

ModelResult EvaluateActive(const ParticleState& particle,
                           const LocalState& local) {
  return ActiveModelFunction(particle, local, ActiveConfiguration);
}

Status EvaluateBatch(const std::vector<ParticleState>& particles,
                     const std::vector<LocalState>& states,
                     const ModelConfiguration& configuration,
                     std::vector<ModelResult>* results) {
  if (!results)
    return Status::Error(StatusCode::InvalidConfiguration,
                         "null batch output");
  if (particles.size()!=states.size())
    return Status::Error(StatusCode::InvalidConfiguration,
                         "particle and local-state batch shapes differ");
  const Status valid=ValidateConfiguration(configuration);
  if (!valid.ok()) return valid;
  std::vector<ModelResult> candidate;
  candidate.reserve(particles.size());
  for (std::size_t i=0;i<particles.size();++i)
    candidate.push_back(Evaluate(particles[i],states[i],configuration));
  *results=std::move(candidate);
  return Status::Success();
}

Status AssembleSymmetricTensor(
    const SpatialCoefficients& c,
    std::array<std::array<double,3>,3>* tensor) {
  if (!tensor)
    return Status::Error(StatusCode::InvalidInput,"null tensor output");
  std::array<std::array<double,3>,3> value{};
  if (c.frame.kind==FrameKind::Isotropic) {
    if (!c.isotropicM2PerS.has_value() ||
        !Finite(*c.isotropicM2PerS) || *c.isotropicM2PerS < 0.0)
      return Status::Error(StatusCode::InvalidInput,
                           "isotropic frame requires nonnegative isotropic coefficient");
    for (int i=0;i<3;++i) value[i][i]=*c.isotropicM2PerS;
  } else {
    if (!c.frame.b.has_value() || !c.frame.e1.has_value() ||
        !c.frame.e2.has_value() || !c.parallelM2PerS.has_value() ||
        !c.perpendicular1M2PerS.has_value() ||
        !c.perpendicular2M2PerS.has_value())
      return Status::Error(StatusCode::MissingInput,
                           "ordered tensor requires three axes and three eigenvalues");
    const std::array<std::array<double,3>,3> axes{{*c.frame.b,*c.frame.e1,*c.frame.e2}};
    const std::array<double,3> eigen{{*c.parallelM2PerS,
        *c.perpendicular1M2PerS,*c.perpendicular2M2PerS}};
    for (int a=0;a<3;++a) {
      double norm=0.0;
      for (double x:axes[a]) norm+=x*x;
      if (std::fabs(norm-1.0)>1.0e-10 || !Finite(eigen[a]) ||
          eigen[a] < 0.0)
        return Status::Error(StatusCode::InvalidInput,
                             "tensor frame/eigenvalue is invalid");
      for (int b=a+1;b<3;++b) {
        double dot=0.0; for(int i=0;i<3;++i) dot+=axes[a][i]*axes[b][i];
        if (std::fabs(dot)>1.0e-10)
          return Status::Error(StatusCode::InvalidInput,
                               "tensor axes are not orthogonal");
      }
      for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        value[i][j]+=eigen[a]*axes[a][i]*axes[a][j];
    }
  }
  *tensor=value;
  return Status::Success();
}

}  // namespace PerpendicularDiffusion
}  // namespace SEP

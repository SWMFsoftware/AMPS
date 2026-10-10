#include "parallel_diffusion_adapter.h"

#include "parallel_diffusion/parallel_diffusion.h"

#include <cmath>
#include <vector>

namespace SEP {
namespace ParallelDiffusionBinding {
namespace {

namespace PD = SEP::ParallelDiffusion;
using Transport::Status;
using Transport::StatusCode;

// Set only by a successful Configure() during serial startup.
bool gInstalled = false;

Status Unsupported(const std::string& message) {
  return Status::Error(StatusCode::UnsupportedConfiguration, message);
}

Status Invalid(const std::string& message) {
  return Status::Error(StatusCode::InvalidCoefficient, message);
}

bool FinitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

}  // namespace

Status Configure(bool librarySelected, const std::string& moverName,
                 Initialization::Configuration* initialization) {
  const bool selected = librarySelected;
  const bool sectionPresent =
      initialization != NULL && initialization->hasParallelDiffusionSection;

  // Selection is two-sided, as in srcSEP3D schema 5: the CLI provider names
  // the coefficient source and the input file supplies its model.  Either
  // half alone is an operator error, not an implicit default.
  if (!selected && sectionPresent)
    return Unsupported("[parallel_diffusion] is present in the --input file "
                       "but --spatial-diffusion-provider is not "
                       "parallel-diffusion-library");
  if (!selected) return Status::Ok();
  if (initialization == NULL)
    return Unsupported("--spatial-diffusion-provider parallel-diffusion-library "
                       "requires --input with a schema-4 [parallel_diffusion] "
                       "section");
  if (!sectionPresent)
    return Unsupported("--spatial-diffusion-provider parallel-diffusion-library "
                       "requires a [parallel_diffusion] section in the --input "
                       "file (run.schema_version=4)");
  // The registry's ValidateMoverCompatibility already enforces this for CLI
  // input; repeat it at the installation boundary for non-CLI callers.
  if (moverName != "parker")
    return Unsupported("the parallel-diffusion library supplies the Parker "
                       "spatial coefficient and is available only to the "
                       "parker mover");

  std::vector<PD::SectionLine> lines;
  lines.reserve(initialization->parallelDiffusionLines.size());
  for (const Initialization::ParallelDiffusionLine& line :
       initialization->parallelDiffusionLines) {
    PD::SectionLine record;
    record.lineNumber = line.line;
    record.text = line.text;
    lines.push_back(record);
  }
  PD::ParsedSection parsed;
  const PD::Status parsedStatus = PD::ParseSection(lines, &parsed);
  if (!parsedStatus.ok())
    return Invalid("invalid [parallel_diffusion] section: " +
                   parsedStatus.detail);

  // srcSEP constructs ParticleState from the compiled species' mass and
  // charge and LocalState from the snapshot position, B vector, and time.  It
  // has no authoritative nucleon count (PIC exposes mass only; a rounded
  // mass/m_p is not a nucleon count), no slab/2-D variance or bend-over
  // length (its turbulence state is E+/E- Alfven-wave energy, and no mapping
  // to the library's canonical spectra is specified), no external
  // time/region/radial factors, and no effective-field definition.  All
  // HostInputAvailability flags therefore stay false.
  const PD::HostInputAvailability srcSepInputs;
  const PD::Status hostStatus =
      PD::CheckHostInputAvailability(parsed.configuration, srcSepInputs);
  if (!hostStatus.ok())
    return Unsupported("srcSEP cannot evaluate the selected parallel-"
                       "diffusion configuration: " + hostStatus.detail);

  const PD::Status installed = PD::SetActiveConfiguration(parsed.configuration);
  if (!installed.ok())
    return Invalid("cannot install parallel-diffusion model: " +
                   installed.detail);
  if (PD::ActiveParallelDiffusion !=
      PD::BoundFunctionForModel(parsed.configuration.model))
    return Invalid("mover-facing parallel-diffusion pointer does not match "
                   "the installed model");

  initialization->parallelDiffusionModelId =
      PD::ModelName(parsed.configuration.model);
  initialization->parallelDiffusionConfigurationFingerprint =
      parsed.configurationFingerprint;
  gInstalled = true;
  return Status::Ok();
}

bool IsInstalled() { return gInstalled; }

KappaSample EvaluateKappa(double massKg, double signedChargeC,
                          double momentumKgMPerS, const double positionM[3],
                          const double meanFieldT[3], double timeS,
                          std::uint64_t backgroundRevision) {
  KappaSample sample;
  if (!gInstalled) {
    sample.status = Invalid("parallel-diffusion library selected but no "
                            "model was installed before particle transport");
    return sample;
  }
  PD::ParticleState particle;
  particle.massKg = massKg;
  particle.chargeC = signedChargeC;
  particle.momentumKgMPerS = momentumKgMPerS;
  // nucleonCount deliberately remains absent (see Configure()).

  PD::LocalState local;
  local.timeS = timeS;
  local.positionM = {{positionM[0], positionM[1], positionM[2]}};
  local.meanFieldT = std::array<double, 3>{{meanFieldT[0], meanFieldT[1],
                                            meanFieldT[2]}};
  local.backgroundRevision = backgroundRevision;
  // turbulence, effective field, external factors, and a supplied
  // perpendicular coefficient are not populated: srcSEP has no approved
  // source for those exact quantities.

  const PD::ParallelResult result = PD::ActiveParallelDiffusion(particle, local);
  sample.modelId = result.provenance.evaluatedModelId;
  sample.configurationFingerprint = result.provenance.configurationFingerprint;
  if (!result.status.ok()) {
    sample.status = Status::Error(StatusCode::UnresolvedCoefficient,
        "parallel-diffusion model '" + sample.modelId + "' failed: " +
        result.status.detail);
    return sample;
  }
  if (!result.kappaParallelM2PerS.has_value() ||
      !result.lambdaParallelM.has_value() ||
      !FinitePositive(*result.kappaParallelM2PerS) ||
      !FinitePositive(*result.lambdaParallelM)) {
    sample.status = Status::Error(StatusCode::UnresolvedCoefficient,
        "parallel-diffusion library returned an incomplete or non-positive "
        "kappa/lambda pair");
    return sample;
  }
  sample.kappaParallelM2PerS = *result.kappaParallelM2PerS;
  sample.lambdaParallelM = *result.lambdaParallelM;
  sample.status = Status::Ok();
  return sample;
}

}  // namespace ParallelDiffusionBinding
}  // namespace SEP

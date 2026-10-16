// ============================================================================
// Coefficient-kernel certification matrix and readiness check (roadmap stage
// PS05; specification Section 16.11).
//
// Release is a property of a (model, kernel) pair and of the step it runs
// with.  Each row of pitch_angle_certified_pairs.inc is the evidence of the
// acceptance campaign test_pitch_angle_montecarlo.cpp recorded in
// IMPLEMENTATION_STATUS.md: it is added only after every case of that pair
// passes its predeclared budget at the stated largest nu_max*dt.  Outside
// the matrix the check fails with IncompatibleKernel and names the reason;
// no pair is ever released by default.
//
// This file is deliberately separate from pitch_angle_diffusion.cpp: the
// acceptance binary does not link it, so recording evidence here cannot
// change the binary that produced that evidence.
// ============================================================================

#include "pitch_angle_diffusion.h"

#include <sstream>

namespace SEP {
namespace PitchAngleDiffusion {

const std::vector<PairCertification>& CertifiedPairs() {
  static const std::vector<PairCertification> pairs = {
#include "pitch_angle_certified_pairs.inc"
  };
  return pairs;
}

Status CheckKernelPair(const ModelConfiguration& configuration, Kernel kernel) {
  ShapeConstants constants;
  const Status valid = ValidateConfiguration(configuration, &constants);
  if (!valid.ok()) return valid;
  const std::string pair = std::string(ModelName(configuration.model)) +
                           " with " + Name(kernel);
  if (configuration.model == ModelId::QForm &&
      (kernel == Kernel::EulerMaruyama || kernel == Kernel::Milstein))
    return Status::Error(StatusCode::IncompatibleKernel,
        pair + ": the Ito drift of Eq. (163) is unbounded at mu = 0 for every "
        "H >= 0; the pair is not released (Sections 16.5, 16.11, user "
        "decision 2026-10-13)");
  const double scaled = WorstCaseScaledStep(configuration, constants);
  for (const PairCertification& entry : CertifiedPairs()) {
    if (entry.model != configuration.model || entry.kernel != kernel) continue;
    if (configuration.model == ModelId::EpsilonForm &&
        (configuration.shape.epsilon < entry.epsilonMinimum ||
         configuration.shape.epsilon > entry.epsilonMaximum)) {
      std::ostringstream message;
      message << pair << ": epsilon = " << configuration.shape.epsilon
              << " is outside the certified range [" << entry.epsilonMinimum
              << ", " << entry.epsilonMaximum << "] (" << entry.evidence << ")";
      return Status::Error(StatusCode::IncompatibleKernel, message.str());
    }
    if (!(scaled <= entry.maximumScaledStep)) {
      std::ostringstream message;
      message << pair << ": the step rule " << Name(configuration.step.rule)
              << " allows nu_max*dt up to " << scaled
              << ", above the certified " << entry.maximumScaledStep << " ("
              << entry.evidence << ")";
      return Status::Error(StatusCode::IncompatibleKernel, message.str());
    }
    return Status::Success();
  }
  return Status::Error(StatusCode::IncompatibleKernel,
      pair + " has no certification entry (see IMPLEMENTATION_STATUS.md)");
}

}  // namespace PitchAngleDiffusion
}  // namespace SEP

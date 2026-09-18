#ifndef SEP_UTIL_SEP_CONFIGURATION_MATRIX_H
#define SEP_UTIL_SEP_CONFIGURATION_MATRIX_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_registry.h)
#include "sep_production_mover.h"
#include "sep_turbulence_core.h"

#include <string>
#include <vector>

namespace SEP {
namespace ConfigurationMatrix {

enum class SupportStatus { Supported, Conditional, Unsupported };

struct Combination {
  Mover::ProductionMover mover = Mover::ProductionMover::Parker;
  Transport::Coefficient::SourceMode coefficientSource =
      Transport::Coefficient::SourceMode::Prescribed;
  Turbulence::Source turbulenceSource = Turbulence::Source::Prescribed;
  Turbulence::CouplingPolicy coupling = Turbulence::CouplingPolicy::Disabled;
};

struct Classification {
  SupportStatus support = SupportStatus::Unsupported;
  std::string diagnosticCode;
  std::string reason;
  std::string suggestedAlternative;
};

const char* SupportStatusName(SupportStatus status);
Classification Classify(const Combination& combination);

// Conditional SWMF handoff combinations resolve against an explicit phase.
// The function is suitable for preflight because it never reads mutable global
// state and returns a stable diagnostic before particles or wave arrays exist.
Transport::Status Preflight(const Combination& combination,
                            bool swmfHandoffCompleted);

// The complete Cartesian product is generated from production registries and
// enum domains.  Documentation and tests consume the same records, preventing
// a hand-maintained support table from omitting a reachable combination.
std::vector<std::pair<Combination, Classification> > Enumerate();
std::string RenderMarkdown();

}  // namespace ConfigurationMatrix
}  // namespace SEP

#endif  // SEP_UTIL_SEP_CONFIGURATION_MATRIX_H

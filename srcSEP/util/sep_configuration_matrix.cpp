#include "sep_configuration_matrix.h"

#include <sstream>

namespace SEP {
namespace ConfigurationMatrix {
namespace {

bool MatrixEvolvesLocally(Turbulence::Source source) {
  return source == Turbulence::Source::SelfConsistentIntegrated ||
         source == Turbulence::Source::SelfConsistentSpectral ||
         source == Turbulence::Source::SwmfInitialThenEvolveLocal;
}

bool IsSwmf(Turbulence::Source source) {
  return source == Turbulence::Source::SwmfReadOnly ||
         source == Turbulence::Source::SwmfInitialThenEvolveLocal;
}

}  // namespace

const char* SupportStatusName(SupportStatus status) {
  switch (status) {
    case SupportStatus::Supported: return "supported";
    case SupportStatus::Conditional: return "conditional";
    case SupportStatus::Unsupported: return "unsupported";
  }
  return "unclassified";
}

Classification Classify(const Combination& c) {
  Classification result;
  if (c.coupling == Turbulence::CouplingPolicy::StreamingEnergyExchange &&
      !MatrixEvolvesLocally(c.turbulenceSource)) {
    result.diagnosticCode = "CFG-WAVE-IMMUTABLE";
    result.reason = "streaming exchange requires locally owned turbulence";
    result.suggestedAlternative = "disable coupling or select a self-consistent source";
    return result;
  }
  if (c.coefficientSource == Transport::Coefficient::SourceMode::SelfConsistent &&
      !MatrixEvolvesLocally(c.turbulenceSource)) {
    result.diagnosticCode = "CFG-COEFF-NO-LOCAL-WAVES";
    result.reason = "self-consistent coefficients require locally evolved waves";
    result.suggestedAlternative = "use prescribed coefficients or local turbulence";
    return result;
  }
  if (c.coefficientSource == Transport::Coefficient::SourceMode::Swmf &&
      !IsSwmf(c.turbulenceSource)) {
    result.diagnosticCode = "CFG-COEFF-NO-SWMF";
    result.reason = "SWMF coefficients require an SWMF-owned source generation";
    result.suggestedAlternative = "use prescribed coefficients or an SWMF source";
    return result;
  }
  if (c.turbulenceSource == Turbulence::Source::SwmfInitialThenEvolveLocal &&
      (c.coefficientSource == Transport::Coefficient::SourceMode::Swmf ||
       c.coefficientSource == Transport::Coefficient::SourceMode::SelfConsistent)) {
    result.support = SupportStatus::Conditional;
    result.diagnosticCode = "CFG-SWMF-HANDOFF-PHASE";
    result.reason = "coefficient authority changes exactly once at SWMF handoff";
    result.suggestedAlternative =
        "use SWMF coefficients before handoff and self-consistent coefficients after";
    return result;
  }
  result.support = SupportStatus::Supported;
  result.diagnosticCode = "CFG-SUPPORTED";
  result.reason = "mover, coefficient source, turbulence ownership, and coupling agree";
  return result;
}

Transport::Status Preflight(const Combination& c, bool handoffCompleted) {
  const Classification classification = Classify(c);
  if (classification.support == SupportStatus::Unsupported)
    return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
        classification.diagnosticCode + ": " + classification.reason +
        "; alternative: " + classification.suggestedAlternative);
  if (classification.support == SupportStatus::Conditional) {
    const bool correctSource = handoffCompleted
        ? c.coefficientSource == Transport::Coefficient::SourceMode::SelfConsistent
        : c.coefficientSource == Transport::Coefficient::SourceMode::Swmf;
    if (!correctSource)
      return Transport::Status::Error(
          Transport::StatusCode::UnsupportedConfiguration,
          classification.diagnosticCode + ": coefficient source does not match handoff phase");
  }
  return Transport::Status::Ok();
}

std::vector<std::pair<Combination, Classification> > Enumerate() {
  std::vector<std::pair<Combination, Classification> > result;
  const std::vector<Mover::Descriptor>& movers = Mover::Registry();
  const Transport::Coefficient::SourceMode coefficientSources[] = {
      Transport::Coefficient::SourceMode::Prescribed,
      Transport::Coefficient::SourceMode::SelfConsistent,
      Transport::Coefficient::SourceMode::Swmf};
  const Turbulence::Source turbulenceSources[] = {
      Turbulence::Source::Prescribed,
      Turbulence::Source::SelfConsistentIntegrated,
      Turbulence::Source::SelfConsistentSpectral,
      Turbulence::Source::SwmfReadOnly,
      Turbulence::Source::SwmfInitialThenEvolveLocal};
  const Turbulence::CouplingPolicy couplingPolicies[] = {
      Turbulence::CouplingPolicy::Disabled,
      Turbulence::CouplingPolicy::StreamingEnergyExchange};
  for (std::size_t m = 0; m < movers.size(); ++m)
    for (std::size_t c = 0; c < 3; ++c)
      for (std::size_t t = 0; t < 5; ++t)
        for (std::size_t p = 0; p < 2; ++p) {
          Combination combination;
          combination.mover = movers[m].mover;
          combination.coefficientSource = coefficientSources[c];
          combination.turbulenceSource = turbulenceSources[t];
          combination.coupling = couplingPolicies[p];
          result.push_back(std::make_pair(combination, Classify(combination)));
        }
  return result;
}

std::string RenderMarkdown() {
  std::ostringstream out;
  out << "| Mover | Coefficients | Turbulence | Coupling | Support | Code |\n"
      << "| --- | --- | --- | --- | --- | --- |\n";
  const std::vector<std::pair<Combination, Classification> > rows = Enumerate();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const Combination& c = rows[i].first;
    const Classification& s = rows[i].second;
    out << "| " << Mover::Describe(c.mover).canonicalName << " | "
        << Transport::Coefficient::SourceName(c.coefficientSource) << " | "
        << Turbulence::SourceName(c.turbulenceSource) << " | "
        << Turbulence::CouplingPolicyName(c.coupling) << " | "
        << SupportStatusName(s.support) << " | " << s.diagnosticCode << " |\n";
  }
  return out.str();
}

}  // namespace ConfigurationMatrix
}  // namespace SEP

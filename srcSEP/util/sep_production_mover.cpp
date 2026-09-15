#include "sep_production_mover.h"

#include <algorithm>
#include <cctype>
#include <ostream>
#include <stdexcept>

namespace {

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  return value;
}

}  // namespace

const std::vector<SEP::Mover::Descriptor>& SEP::Mover::Registry() {
  // Function-local construction avoids static initialization order coupling to
  // PIC while retaining a single immutable registry for the process lifetime.
  static const std::vector<Descriptor> registry = {
      {ProductionMover::Parker,
       "parker",
       "pitch-angle-averaged Parker transport along a field line",
       {true, false, true, false, CoefficientContract::SpatialDiffusion}},
      {ProductionMover::FocusedTransportDiffusion,
       "fte-dmumu",
       "focused transport with an explicit pitch-angle diffusion coefficient",
       {true, true, true, false, CoefficientContract::PitchAngleDiffusion}},
      {ProductionMover::FocusedTransportMeanFreePath,
       "fte-mfp",
       "focused transport with event scattering derived from a mean free path",
       {true, true, true, false, CoefficientContract::MeanFreePath}},
  };
  return registry;
}

const SEP::Mover::Descriptor& SEP::Mover::Describe(ProductionMover mover) {
  const std::vector<Descriptor>& registry = Registry();
  for (std::size_t i = 0; i < registry.size(); ++i) {
    if (registry[i].mover == mover) return registry[i];
  }
  throw std::invalid_argument("unknown ProductionMover value");
}

const char* SEP::Mover::CoefficientContractName(CoefficientContract contract) {
  switch (contract) {
    case CoefficientContract::SpatialDiffusion: return "spatial-diffusion Dxx";
    case CoefficientContract::PitchAngleDiffusion: return "pitch-angle diffusion Dmumu";
    case CoefficientContract::MeanFreePath: return "parallel mean free path lambda";
  }
  return "unknown";
}

bool SEP::Mover::ParseProductionMover(const std::string& rawName,
                                      ProductionMover& mover,
                                      std::string& warning) {
  const std::string name = Lower(rawName);
  warning.clear();

  const std::vector<Descriptor>& registry = Registry();
  for (std::size_t i = 0; i < registry.size(); ++i) {
    if (name == registry[i].canonicalName) {
      mover = registry[i].mover;
      return true;
    }
  }

  // Step 14 closes the one-release migration window.  Accepting an alias here
  // would keep an undocumented fourth naming surface alive indefinitely and
  // could make two input decks appear different while selecting identical
  // physics.  The migration manifest records replacements; runtime selection
  // now accepts only the three names printed by --list-movers.
  return false;
}

void SEP::Mover::PrintProductionMovers(std::ostream& out) {
  out << "srcSEP production movers:\n";
  const std::vector<Descriptor>& registry = Registry();
  for (std::size_t i = 0; i < registry.size(); ++i) {
    const Descriptor& descriptor = registry[i];
    out << "  " << descriptor.canonicalName << "\n"
        << "    " << descriptor.description << "\n"
        << "    coefficient: "
        << CoefficientContractName(descriptor.capabilities.coefficientContract)
        << "\n"
        << "    field-line attachment: "
        << (descriptor.capabilities.requiresFieldLineAttachment ? "required" : "not required")
        << "; pitch angle: "
        << (descriptor.capabilities.usesPitchAngleState ? "explicit" : "averaged")
        << "; wave streaming: "
        << (descriptor.capabilities.accumulatesWaveStreaming ? "accumulated" : "not accumulated")
        << "\n";
  }
}

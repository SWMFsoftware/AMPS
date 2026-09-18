#ifndef SEP_UTIL_SEP_PRODUCTION_MOVER_H
#define SEP_UTIL_SEP_PRODUCTION_MOVER_H

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace SEP {
namespace Mover {

// This is the complete public production-mover set.  Historical experimental
// and fully three-dimensional movers remain implementation sources for now, but
// cannot be selected through the production registry or standalone CLI.
enum class ProductionMover {
  Parker,
  FocusedTransportDiffusion,
  FocusedTransportMeanFreePath
};

enum class CoefficientContract {
  SpatialDiffusion,
  PitchAngleDiffusion,
  MeanFreePath
};

struct MoverCapabilities {
  bool requiresFieldLineAttachment;
  bool usesPitchAngleState;
  bool accumulatesWaveStreaming;
  bool evolvesWaveStateDirectly;
  CoefficientContract coefficientContract;
};

struct Descriptor {
  ProductionMover mover;
  const char* canonicalName;
  const char* description;
  MoverCapabilities capabilities;
};

// Registry() contains exactly the three stable public choices in deterministic
// order.  CLI help, --list-movers, metadata, and runtime selection all consume
// these descriptors so their names and capabilities cannot drift apart.
const std::vector<Descriptor>& Registry();
const Descriptor& Describe(ProductionMover mover);
const char* CoefficientContractName(CoefficientContract contract);

// Parse canonical names and only aliases with an unambiguous physical mapping.
// A nonempty warning denotes a deprecated alias accepted for this transition
// release.  Ambiguous, direct-wave, drift, Boris, and 3-D names return false.
bool ParseProductionMover(const std::string& rawName,
                          ProductionMover& mover,
                          std::string& warning);
void PrintProductionMovers(std::ostream& out);

// Runtime implementations live in production_mover_runtime.cpp so the registry
// and parser remain dependency-light.  Selection and capability queries are
// declared here because neither exposes PIC implementation types.
void SelectProductionMover(ProductionMover mover);
ProductionMover CurrentProductionMover();
const MoverCapabilities& CurrentCapabilities();

// Return the number of particle-dispatch calls that reached and returned from
// the selected production implementation on this process.  This is an
// integration-evidence counter, not a physical particle population: a mover
// may legitimately delete a particle at a boundary after advancing it, and
// that completed dispatch is still counted.  Selection resets the counter so
// a linked test run cannot inherit evidence from an earlier mover choice.
std::uint64_t CompletedDispatchCount();
void PrintRuntimeConfiguration(std::ostream& out);

}  // namespace Mover
}  // namespace SEP

#endif

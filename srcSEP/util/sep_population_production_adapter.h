#ifndef SRCSEP_UTIL_SEP_POPULATION_PRODUCTION_ADAPTER_H
#define SRCSEP_UTIL_SEP_POPULATION_PRODUCTION_ADAPTER_H
#include "sep_conservative_resampling.h"
#include <string>
#include <vector>
namespace SEP { namespace PopulationControl {
enum class ConservativeMergeStatus { Reduced,SparseNoReduction,NumericalRejection };
struct ConservativeMergeReport {
  ConservativeMergeStatus status=ConservativeMergeStatus::SparseNoReduction;
  int fieldLine=-1,species=-1,segment=-1,before=0,after=0,target=0;
  unsigned sparseBins=0;
  unsigned numericalStops=0;
  double maximumRelativeResidual=0.0;
  // Diagnostics below use actual represented physical weight, not only the
  // correction multiplier stored by the core. Current has units C m/s; stress
  // moments have J before division by a physical sampling volume.
  Conservative::UnpreservedMoments beforeUnpreserved,afterUnpreserved;
  std::string detail;
};
void RegisterConservativeMerger();
// Collective on every rank, once before the core wrapper. Persisted source
// scheduler events are reserved for every initialized physical line/species;
// segment index keys a separate population stream, never a mover RNG stream.
void BeginConservativePopulationPass();
const std::vector<ConservativeMergeReport>& ConservativeMergeReports();
}} // namespace SEP::PopulationControl
#endif

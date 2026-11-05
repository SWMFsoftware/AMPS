#ifndef _GENERIC_COMET_DIAGNOSTICS_H_
#define _GENERIC_COMET_DIAGNOSTICS_H_

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Comet {
namespace Diagnostics {

// Compile-time input selects this runtime constant.  OFF is intentionally the
// source default: calls from hot paths return before touching tables, files,
// MPI, random numbers, or model state.
extern bool Enabled;

void Initialize();

void RecordEligibleParentVisit(std::size_t reaction);
void RecordReactionSelection(std::size_t reaction,bool neutralSink);
void RecordAllocationDeferral(std::size_t reaction);
void RecordEmergencyRequest();

// Record one executed common-factor emergency and the trigger(s) that caused
// it. More than one trigger can be true for one execution; the total counter
// advances once while each applicable trigger counter advances independently.
// The radial callback is collective, so the implementation records on rank
// zero and the sampling writer reduces exactly once.
void RecordEmergencyEvent(bool bufferUtilization,
    bool samplingPopulation,bool allocationBurst);

// Add one fixed-carrier-epoch contribution to the conservative systematic
// discrepancy bound. The caller invokes this once per active concrete
// product accumulator, not once per chemistry event.
void RecordSystematicResidualBound(std::size_t reaction,
    std::size_t productTerm,int species,double carrierRate);

void RecordChemistryProduct(std::size_t reaction,std::size_t productTerm,
    std::uint64_t modelMultiplicity,std::uint64_t allocatedSlots,
    std::uint64_t reusedSlots,double expectedNumberRate,
    double expectedMassRate,const double expectedMomentumRate[3],
    double expectedEnergyRate,double representedNumberRate,
    double representedMassRate,const double representedMomentumRate[3],
    double representedEnergyRate,double representedNumberRateVariance);

void RecordSurfaceInjection(int species,double representedNumberRate);
void RecordChemicalProduction(int species,double expectedNumberRate,
    double representedNumberRate,double representedNumberRateVariance);
void RecordChemicalDestruction(int species,double representedNumberRate);
void RecordInnerBoundary(int species,double representedNumberRate);
void RecordOuterBoundary(int species,double representedNumberRate);

// Collective sampling-window writer.  relativeL1 is indexed by compiled
// species and contains the exact cell-measure-weighted convergence change;
// it is omitted for the first reference window when hasPreviousProfile=false.
void WriteSamplingReport(long int completedIteration,long int outputIndex,
    long int sampleLength,const std::vector<double> &relativeL1,
    bool hasPreviousProfile);

}
}

#endif

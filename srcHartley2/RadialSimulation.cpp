#include "RadialSimulation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "ComponentCarrierPlanner.h"
#include "CorrectionWindowPlanner.h"
#include "pic.h"
#include "CometChemistry.h"
#include "CometDiagnostics.h"
#include "CometModel.h"
#include "CometRadiativeCooling.h"

bool Comet::RadialSimulation::Enabled=false;
double Comet::RadialSimulation::CourantNumber=-1.0;
double Comet::RadialSimulation::CharacteristicSpeed[PIC::nTotalSpecies]={0.0};
double Comet::RadialSimulation::InnerBoundaryTemperature=-1.0;
long int Comet::RadialSimulation::MaxIterations=-1;
long int Comet::RadialSimulation::SamplingLength=-1;
long int Comet::RadialSimulation::LoadBalanceInterval=-1;
bool Comet::RadialSimulation::ParticleWeightAdaptationEnabled=false;
long int Comet::RadialSimulation::ParticleWeightAdaptationInterval=-1;
long int Comet::RadialSimulation::MinimumParticlesPerCell=-1;
long int Comet::RadialSimulation::MaximumParticlesPerCell=-1;
long int Comet::RadialSimulation::SpeciesMinimumParticlesPerCell[PIC::nTotalSpecies]={0};
long int Comet::RadialSimulation::SpeciesMaximumParticlesPerCell[PIC::nTotalSpecies]={0};
int Comet::RadialSimulation::ParticleWeightLocked[PIC::nTotalSpecies]={0};
int Comet::RadialSimulation::PopulationPolicySource[PIC::nTotalSpecies]={0};
int Comet::RadialSimulation::PopulationClass[PIC::nTotalSpecies]={0};
double Comet::RadialSimulation::EmergencyParticleBufferUtilization=-1.0;
double Comet::RadialSimulation::AllocationLowWatermark=-1.0;
double Comet::RadialSimulation::AllocationReserveFraction=-1.0;
long int Comet::RadialSimulation::MaximumFaceTransferMultiplicity=-1;
long int Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction=-1;
bool Comet::RadialSimulation::IndividualCorrectionWindowsEnabled=false;
double Comet::RadialSimulation::IndividualCorrectionWindowLower[PIC::nTotalSpecies]={0.0};
double Comet::RadialSimulation::IndividualCorrectionWindowTarget[PIC::nTotalSpecies]={0.0};
double Comet::RadialSimulation::IndividualCorrectionWindowUpper[PIC::nTotalSpecies]={0.0};
int Comet::RadialSimulation::IndividualCorrectionWindowPolicySource[PIC::nTotalSpecies]={0};
long int Comet::RadialSimulation::MaximumIndividualCorrectionWindowSplits=-1;
double Comet::RadialSimulation::ConvergenceTolerance=-1.0;
int Comet::RadialSimulation::RequiredConvergedWindows=-1;
bool Comet::RadialSimulation::Stage9ScienceDiagnosticsEnabled=false;

namespace {
// Stage-9 monitoring is deliberately kept outside the production diagnostic
// accumulators.  The values below describe one wall-clock sampling window and
// are reset only after its immutable checkpoint has been written.  Keeping
// them process-local avoids a new collective on every particle iteration;
// the checkpoint reduction reports the slowest rank and the largest observed
// per-rank allocation pressure.
long int Stage9WindowPeakLiveParticles=0;
double Stage9WindowPeakBufferUtilization=0.0;
double Stage9WindowStartWallTime=0.0;
unsigned long long Stage9SafeBoundaryAudits=0;

bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
}

void UpdateStage9BufferPeak() {
  if (!Comet::RadialSimulation::Stage9ScienceDiagnosticsEnabled) return;
  const long int live=PIC::ParticleBuffer::NAllPart;
  const long int capacity=PIC::ParticleBuffer::MaxNPart;
  if (live<0 || capacity<=0 || live>capacity) {
    std::fprintf(stderr,
        "invalid particle-buffer state in Stage-9 monitoring: live=%ld capacity=%ld\n",
        live,capacity);
    std::abort();
  }
  Stage9WindowPeakLiveParticles=std::max(
      Stage9WindowPeakLiveParticles,live);
  Stage9WindowPeakBufferUtilization=std::max(
      Stage9WindowPeakBufferUtilization,
      static_cast<double>(live)/static_cast<double>(capacity));
}

// Local-timestep continuation particles and Stage-8 event markers are valid
// only inside the production mover's drained execution track.  Population
// adaptation, output, and AMR redistribution must never observe either state.
// Scan at those comparatively infrequent boundaries instead of every
// iteration so the science diagnostic cannot dominate a 40,000-step run.
void AuditStage9SafeBoundary(const char *boundary) {
  if (!Comet::RadialSimulation::Stage9ScienceDiagnosticsEnabled) return;
  long long live=0,pendingRemainders=0,pendingReactionEvents=0;
  const int cellCount=
      _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL) continue;
    for (int cell=0;cell<cellCount;cell++)
      for (long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;
           ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        live++;
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        const double remainder=
            PIC::Mover::GetLocalTimeStepRemainderFraction(data);
        if (!std::isfinite(remainder) || remainder!=0.0)
          pendingRemainders++;
        if (Comet::Chemistry::HasPendingTotalReactionEvent(data))
          pendingReactionEvents++;
      }
  }
  if (live!=PIC::ParticleBuffer::NAllPart || pendingRemainders!=0 ||
      pendingReactionEvents!=0) {
    std::fprintf(stderr,
        "Stage-9 unsafe %s boundary on rank %d: list-live=%lld buffer-live=%ld pending-remainders=%lld pending-reaction-events=%lld\n",
        boundary,PIC::ThisThread,live,PIC::ParticleBuffer::NAllPart,
        pendingRemainders,pendingReactionEvents);
    std::abort();
  }
  Stage9SafeBoundaryAudits++;
}

bool HasConfiguredDirectSource(int species) {
  return Comet::Model::SurfaceProductionRate[species]>0.0 ||
      Comet::Model::ExtendedProductionRate[species]>0.0;
}

double ParticlePhysicalRate(long int ptr,int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL || node->block==NULL) {
    std::fprintf(stderr,"missing radial block in particle-rate diagnostic\n");
    std::abort();
  }
  PIC::ParticleBuffer::byte *data=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  double weight=node->block->GetLocalParticleWeight(species);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  weight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
#endif
  const double timeStep=node->block->GetLocalTimeStep(species);
  if (!(weight>0.0) || !std::isfinite(weight) ||
      !(timeStep>0.0) || !std::isfinite(timeStep)) {
    std::fprintf(stderr,
        "invalid particle weight/timestep in radial flux diagnostic\n");
    std::abort();
  }
  return weight/timeStep;
}

bool HasInnerFace(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  bool externalFaces[2*DIM];
  if (PIC::Mesh::mesh->ExternalBoundaryBlock(node,externalFaces)!=
      _EXTERNAL_BOUNDARY_BLOCK_) return false;
  return externalFaces[0];
}

cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *FindInnerBoundaryNode(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return NULL;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_)
    return node->block!=NULL && HasInnerFace(node) ? node : NULL;

  for (int child=0;child<(1<<DIM);child++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *result=
        FindInnerBoundaryNode(node->downNode[child]);
    if (result!=NULL) return result;
  }
  return NULL;
}

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
void SetWeightOverTimeStep(int species,double weightOverTimeStep,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->block!=NULL) {
      // W/dt is held constant with radius.  Because dt is proportional to the
      // local cell width, a steady spherical flux then uses approximately the
      // same number of model particles in every logarithmic radial cell.  The
      // face mover converts between any unequal neighboring W/dt values, but
      // this radial policy deliberately makes them equal for a given species
      // so ordinary radial transport needs neither systematic multiplication
      // nor systematic roulette.
      const double localTimeStep=node->block->GetLocalTimeStep(species);
      const double localWeight=weightOverTimeStep*localTimeStep;
      if (!(localTimeStep>0.0) || !std::isfinite(localTimeStep) ||
          !(localWeight>0.0) || !std::isfinite(localWeight)) {
        std::fprintf(stderr,
            "invalid radial local weight/timestep for species %d\n",species);
        std::abort();
      }
      node->block->SetLocalParticleWeight(localWeight,species);
    }
    return;
  }

  for (int child=0;child<(1<<DIM);child++)
    SetWeightOverTimeStep(species,weightOverTimeStep,node->downNode[child]);
}

const char *PopulationPolicySourceName(int source) {
  switch (source) {
  case Comet::RadialSimulation::PolicySource::Global: return "global";
  case Comet::RadialSimulation::PolicySource::DirectSource:
    return "DIRECT_SOURCE";
  case Comet::RadialSimulation::PolicySource::ChemicalOnly:
    return "CHEMICAL_ONLY";
  case Comet::RadialSimulation::PolicySource::ExplicitSpecies:
    return "explicit-species";
  }
  return "invalid";
}

const char *PopulationClassName(int populationClass) {
  switch (populationClass) {
  case Comet::RadialSimulation::PolicyClass::Other: return "other";
  case Comet::RadialSimulation::PolicyClass::DirectSource:
    return "direct-source";
  case Comet::RadialSimulation::PolicyClass::ChemicalOnly:
    return "chemical-only";
  }
  return "invalid";
}

void ReportResolvedPopulationPolicies(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *innerNode) {
  // The generated policy table is immutable, but the surface carrier is a
  // runtime result of AMPS source-weight initialization.  Report both in one
  // startup table so validation can prove that changing a daughter target did
  // not change a locked direct-source injection carrier or physical rate.
  std::vector<double> localCarrier(PIC::nTotalSpecies,0.0);
  std::vector<double> surfaceCarrier(PIC::nTotalSpecies,0.0);
  if (innerNode!=NULL && innerNode->block!=NULL)
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const double weight=innerNode->block->GetLocalParticleWeight(species);
      const double timeStep=innerNode->block->GetLocalTimeStep(species);
      if (!(weight>0.0) || !std::isfinite(weight) ||
          !(timeStep>0.0) || !std::isfinite(timeStep)) {
        std::fprintf(stderr,
            "invalid surface carrier while reporting species %d policy\n",
            species);
        std::abort();
      }
      localCarrier[species]=weight/timeStep;
    }
  MPI_Allreduce(localCarrier.data(),surfaceCarrier.data(),PIC::nTotalSpecies,
      MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread!=0) return;
  std::fprintf(PIC::DiagnospticMessageStream,
      "$PREFIX: resolved radial population policies adaptation-cadence=%s",
      Comet::RadialSimulation::ParticleWeightAdaptationInterval>0 ?
          "fixed-iteration" : "legacy-sampling-break");
  if (Comet::RadialSimulation::ParticleWeightAdaptationInterval>0)
    std::fprintf(PIC::DiagnospticMessageStream," interval=%ld",
        Comet::RadialSimulation::ParticleWeightAdaptationInterval);
  std::fprintf(PIC::DiagnospticMessageStream,"\n");
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    char symbol[_MAX_STRING_LENGTH_PIC_];
    PIC::MolecularData::GetChemSymbol(symbol,species);
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: population-policy species=%s class=%s selector=%s locked=%s range=%ld:%ld surface-source-s-1=%.17e surface-carrier-W-dt-s-1=%.17e\n",
        symbol,PopulationClassName(
            Comet::RadialSimulation::PopulationClass[species]),
        PopulationPolicySourceName(
            Comet::RadialSimulation::PopulationPolicySource[species]),
        Comet::RadialSimulation::ParticleWeightLocked[species] ? "yes" : "no",
        Comet::RadialSimulation::SpeciesMinimumParticlesPerCell[species],
        Comet::RadialSimulation::SpeciesMaximumParticlesPerCell[species],
        Comet::Model::SurfaceProductionRate[species],surfaceCarrier[species]);
  }
  std::fflush(PIC::DiagnospticMessageStream);
}

void ReportResolvedCorrectionWindowPolicies() {
  if (PIC::ThisThread!=0) return;
  std::fprintf(PIC::DiagnospticMessageStream,
      "$PREFIX: individual correction-window controller mode=%s max-split-particles=%ld cadence=completed-time-step\n",
      Comet::RadialSimulation::IndividualCorrectionWindowsEnabled ?
          "on" : "off",
      Comet::RadialSimulation::MaximumIndividualCorrectionWindowSplits);
  if (Comet::RadialSimulation::IndividualCorrectionWindowsEnabled)
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      char symbol[_MAX_STRING_LENGTH_PIC_];
      PIC::MolecularData::GetChemSymbol(symbol,species);
      std::fprintf(PIC::DiagnospticMessageStream,
          "$PREFIX: correction-window-policy species=%s selector=%s lower=%.17e target=%.17e upper=%.17e direct-source=%s\n",
          symbol,PopulationPolicySourceName(
              Comet::RadialSimulation::
                  IndividualCorrectionWindowPolicySource[species]),
          Comet::RadialSimulation::IndividualCorrectionWindowLower[species],
          Comet::RadialSimulation::IndividualCorrectionWindowTarget[species],
          Comet::RadialSimulation::IndividualCorrectionWindowUpper[species],
          HasConfiguredDirectSource(species) ? "yes" : "no");
    }
  std::fflush(PIC::DiagnospticMessageStream);
}

bool RegularAdaptationDue(long int completedIterations,
    bool completeSamplingBreak) {
  if (!Comet::RadialSimulation::ParticleWeightAdaptationEnabled) return false;
  if (Comet::RadialSimulation::ParticleWeightAdaptationInterval>0)
    return completedIterations%
        Comet::RadialSimulation::ParticleWeightAdaptationInterval==0;
  // An input without the Stage-1 cadence record follows the exact historical
  // trigger: one regular adaptation after every complete sampling window.
  return completeSamplingBreak;
}

constexpr long int BlockCellCount=
    _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;

// A completed sampling window already pays the cost of counting every
// species in every owner block.  Use that same snapshot to detect a local
// population excursion that is severe even when the owning MPI rank still
// has ample aggregate buffer capacity.  The threshold requested for this
// second emergency trigger is
//
//   10 * nSpecies * maximum requested particles per species per cell.
//
// Counts are accumulated for a whole block, so the comparison below uses the
// algebraically equivalent threshold multiplied by BlockCellCount.  The
// pending request is set identically on every MPI rank from a global maximum.
// It is deliberately consumed only by the existing pre-particle-exchange
// callback: sampling occurs after the preceding TimeStep has completed, when
// calling particle exchange recursively would be unsafe.
constexpr long long SamplingPopulationEmergencyMultiplier=10;
bool SamplingPopulationEmergencyPending=false;
bool AllocationBurstEmergencyPending=false;
long long SamplingPopulationEmergencyMaximumBlockParticles=0;
long long SamplingPopulationEmergencyThresholdBlockParticles=0;
long int SamplingPopulationEmergencyIteration=-1;

long int MaximumRequestedParticlesPerCell() {
  long int result=0;
  for (int species=0;species<PIC::nTotalSpecies;species++)
    result=std::max(result,
        Comet::RadialSimulation::SpeciesMaximumParticlesPerCell[species]);
  return result;
}

long long GetSamplingPopulationEmergencyBlockThreshold() {
  const long int maximumRequested=MaximumRequestedParticlesPerCell();
  if (maximumRequested<=0 ||
      PIC::nTotalSpecies<=0 || BlockCellCount<=0) {
    std::fprintf(stderr,
        "invalid sampling-population emergency threshold inputs\n");
    std::abort();
  }

  // Form the integer threshold with explicit overflow checks.  Comparing
  // integer block populations avoids a floating-point rounding ambiguity at
  // the exact requested boundary; only a strict exceedance raises an alarm.
  long long threshold=SamplingPopulationEmergencyMultiplier;
  if (PIC::nTotalSpecies>LLONG_MAX/threshold) {
    std::fprintf(stderr,
        "sampling-population emergency species threshold overflows long long\n");
    std::abort();
  }
  threshold*=PIC::nTotalSpecies;
  if (maximumRequested>LLONG_MAX/threshold) {
    std::fprintf(stderr,
        "sampling-population emergency per-cell threshold overflows long long\n");
    std::abort();
  }
  threshold*=maximumRequested;
  if (BlockCellCount>LLONG_MAX/threshold) {
    std::fprintf(stderr,
        "sampling-population emergency block threshold overflows long long\n");
    std::abort();
  }
  return threshold*BlockCellCount;
}

long long GetTotalBlockParticleCount(
    const std::array<long int,PIC::nTotalSpecies> &counts) {
  long long total=0;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (counts[species]<0 || counts[species]>LLONG_MAX-total) {
      std::fprintf(stderr,
          "sampling-population emergency count overflows long long\n");
      std::abort();
    }
    total+=counts[species];
  }
  return total;
}

struct cWeightAdaptationPlan {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node=NULL;
  int Species=-1;
  long int ParticleCount=0;
  // Positive values double the statistical weight and roulette particles;
  // negative values halve the weight and split every represented particle.
  int WeightPower=0;
  // A zero-population photochemical daughter has no particles to roulette.
  // Its weight can nevertheless be raised before its first reaction event,
  // preventing residual resampling from creating an unbounded clone burst.
  // This flag is deliberately limited to positive powers.  An empty,
  // unlocked parent can instead receive a negative chemistry power when its
  // locked daughters cannot be coarsened; that safe, allocation-free weight
  // decrease follows the ordinary split path (whose particle list is empty).
  bool EmptyPhotochemicalDaughter=false;
  // True when the final power was raised above the population controller's
  // choice to impose the hard all-product reaction allocation bound.  The
  // species may be present; in that case normal factor-two roulette is
  // required rather than a weight-only edit.
  bool PhotochemicalCarrierConstraint=false;
};

struct cLockedWeightSnapshot {
  PIC::Mesh::cDataBlockAMR *Block=NULL;
  int Species=-1;
  double Weight=0.0;
  std::uint64_t WeightBits=0;
};

std::uint64_t DoubleBits(double value) {
  static_assert(sizeof(std::uint64_t)==sizeof(double),
      "double bit snapshot requires a 64-bit representation");
  std::uint64_t result=0;
  std::memcpy(&result,&value,sizeof(result));
  return result;
}

long int CountSpeciesInBlock(PIC::Mesh::cDataBlockAMR *block,int species) {
  long int count=0;
  for (long int cell=0;cell<BlockCellCount;cell++) {
    for (long int ptr=block->FirstCellParticleTable[cell];ptr!=-1;
         ptr=PIC::ParticleBuffer::GetNext(ptr)) {
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
      // Adaptation runs only after TimeStep() has completed all continuation
      // sweeps.  A nonzero marker would mean the particle is still in flight
      // and cannot safely be cloned or removed from its ordinary cell list.
      if (PIC::Mover::GetLocalTimeStepRemainderFraction(data)!=0.0) {
        std::fprintf(stderr,
            "particle-weight adaptation found an unfinished trajectory\n");
        std::abort();
      }
#endif
      if (PIC::ParticleBuffer::GetI(data)==species) count++;
    }
  }
  return count;
}

std::array<long int,PIC::nTotalSpecies> CountAllSpeciesInBlock(
    PIC::Mesh::cDataBlockAMR *block) {
  std::array<long int,PIC::nTotalSpecies> counts{};
  for (long int cell=0;cell<BlockCellCount;cell++) {
    for (long int ptr=block->FirstCellParticleTable[cell];ptr!=-1;
         ptr=PIC::ParticleBuffer::GetNext(ptr)) {
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
      if (PIC::Mover::GetLocalTimeStepRemainderFraction(data)!=0.0) {
        std::fprintf(stderr,
            "particle-weight adaptation found an unfinished trajectory\n");
        std::abort();
      }
#endif
      const int species=PIC::ParticleBuffer::GetI(data);
      if (species<0 || species>=PIC::nTotalSpecies) {
        std::fprintf(stderr,
            "particle-weight adaptation found invalid species %d\n",species);
        std::abort();
      }
      counts[species]++;
    }
  }
  return counts;
}

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
// Compact counters returned by the correction-window transaction.  The
// planner itself is AMPS-independent (CorrectionWindowPlanner.h); these
// fields describe only the particle-buffer commit performed here.
struct cCorrectionWindowStatistics {
  long long Examined=0;
  long long Kept=0;
  long long RouletteCandidates=0;
  long long RouletteSurvivors=0;
  long long SplitParents=0;
  long long SplitChildren=0;
  long long RequiredAdditionalParticles=0;
};

std::vector<Comet::CorrectionWindowPlanner::cPolicy>
BuildCorrectionWindowPolicies() {
  std::vector<Comet::CorrectionWindowPlanner::cPolicy> policies;
  policies.reserve(PIC::nTotalSpecies);
  for (int species=0;species<PIC::nTotalSpecies;species++)
    policies.push_back(Comet::CorrectionWindowPlanner::cPolicy(
        Comet::RadialSimulation::IndividualCorrectionWindowLower[species],
        Comet::RadialSimulation::IndividualCorrectionWindowTarget[species],
        Comet::RadialSimulation::IndividualCorrectionWindowUpper[species],
        static_cast<std::size_t>(
            Comet::RadialSimulation::MaximumIndividualCorrectionWindowSplits)));
  return policies;
}

// Apply one collective, mutation-free-preflighted correction-window batch.
// The capacityOverride argument exists only for the focused rollback test;
// production passes -1 and reserves the same allocation headroom used by
// chemistry.  Two mesh traversals are intentional.  Storing a plan object for
// every live particle would itself cause an avoidable memory burst in a
// multi-million-particle comet run.  The first traversal validates every
// correction and counts the exact additional slots; only after all ranks
// accept does the second traversal rebuild the deterministic per-particle
// plan and commit it.
bool ApplyIndividualCorrectionWindowsInternal(long int completedIterations,
    long long capacityOverride,cCorrectionWindowStatistics *statistics,
    bool emitDiagnostics) {
  namespace Window=Comet::CorrectionWindowPlanner;
  cCorrectionWindowStatistics local;
  if (!Comet::RadialSimulation::IndividualCorrectionWindowsEnabled) {
    if (statistics!=NULL) *statistics=local;
    return true;
  }

  const std::vector<Window::cPolicy> policies=
      BuildCorrectionWindowPolicies();
  bool localPlanAccepted=true;
  std::string localFailure;

  // Preflight only authoritative owner blocks.  At a completed AMPS time
  // step no live particle belongs in a boundary-copy list; processing copies
  // would double-count and could publish a correction on the wrong rank.
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL) continue;
    for (long int cell=0;cell<BlockCellCount;cell++)
      for (long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;
           ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
        // A nonzero remainder denotes an unfinished face continuation.  The
        // controller is deliberately forbidden from cloning/deleting that
        // state, even if a future caller invokes it at the wrong location.
        if (PIC::Mover::GetLocalTimeStepRemainderFraction(data)!=0.0) {
          localPlanAccepted=false;
          localFailure="correction-window controller found an unfinished trajectory";
          continue;
        }
#endif
        const int species=PIC::ParticleBuffer::GetI(data);
        const double correction=
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
        try {
          const Window::cParticlePlan plan=Window::BuildParticlePlan(
              Window::cRequest(species,correction),policies.at(species));
          const std::size_t additional=plan.OutputCorrections.size()-1;
          if (additional>static_cast<std::size_t>(LLONG_MAX) ||
              local.RequiredAdditionalParticles>
                  LLONG_MAX-static_cast<long long>(additional))
            throw std::length_error(
                "correction-window additional-particle count overflows");
          local.RequiredAdditionalParticles+=
              static_cast<long long>(additional);
          local.Examined++;
        }
        catch (const std::exception &exception) {
          localPlanAccepted=false;
          if (localFailure.empty()) localFailure=exception.what();
        }
      }
  }

  // Do not count stochastic roulette deletions as available capacity.  A
  // split transaction must be safe even if every low-correction particle
  // survives.  The normal ceiling also preserves the input chemistry/MPI
  // reserve, matching the allocation contract used by ProcessParticle().
  long long available=0;
  if (capacityOverride>=0) available=capacityOverride;
  else {
    const long double ceilingValue=std::floor(
        static_cast<long double>(PIC::ParticleBuffer::MaxNPart)*
        (1.0L-Comet::RadialSimulation::AllocationReserveFraction));
    if (ceilingValue<0.0L || ceilingValue>LLONG_MAX) {
      localPlanAccepted=false;
      if (localFailure.empty())
        localFailure="correction-window capacity ceiling is invalid";
    }
    else {
      const long long ceiling=static_cast<long long>(ceilingValue);
      available=ceiling>PIC::ParticleBuffer::NAllPart ?
          ceiling-PIC::ParticleBuffer::NAllPart : 0;
    }
  }
  if (local.RequiredAdditionalParticles>available) {
    localPlanAccepted=false;
    if (localFailure.empty())
      localFailure="insufficient particle capacity for correction-window batch";
  }

  const int localAccepted=localPlanAccepted ? 1 : 0;
  int globallyAccepted=0;
  MPI_Allreduce(&localAccepted,&globallyAccepted,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (globallyAccepted==0) {
    if (emitDiagnostics && !localPlanAccepted)
      std::fprintf(stderr,
          "rank %d rejected correction-window transaction: %s (required=%lld available=%lld)\n",
          PIC::ThisThread,localFailure.c_str(),
          local.RequiredAdditionalParticles,available);
    if (statistics!=NULL) *statistics=local;
    return false;
  }

  // The complete batch is now known to fit on every rank.  Snapshot `next`
  // before each mutation: split children are inserted at the list head and
  // must not be processed again in this batch; a roulette deletion invalidates
  // only the current link.  CloneParticle preserves position, velocity,
  // internal modes, chemistry stage, and all application fields.  Tracker
  // identity is the sole field intentionally regenerated for each sibling.
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL) continue;
    for (long int cell=0;cell<BlockCellCount;cell++) {
      long int &first=node->block->FirstCellParticleTable[cell];
      for (long int ptr=first;ptr!=-1;) {
        const long int next=PIC::ParticleBuffer::GetNext(ptr);
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        const int species=PIC::ParticleBuffer::GetI(data);
        const double correction=
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
        const Window::cParticlePlan plan=Window::BuildParticlePlan(
            Window::cRequest(species,correction),policies[species]);

        if (plan.Action==Window::eAction::Keep) local.Kept++;
        else if (plan.Action==Window::eAction::Roulette) {
          local.RouletteCandidates++;
          if (rnd()<plan.SurvivalProbability) {
            PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
                plan.OutputCorrections[0],data);
            local.RouletteSurvivors++;
          }
          else PIC::ParticleBuffer::DeleteParticle(ptr,first);
        }
        else {
          local.SplitParents++;
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
              plan.OutputCorrections[0],data);
          for (std::size_t child=1;
               child<plan.OutputCorrections.size();child++) {
            const long int clone=PIC::ParticleBuffer::GetNewParticle(first);
            PIC::ParticleBuffer::byte *copy=
                PIC::ParticleBuffer::GetParticleDataPointer(clone);
            PIC::ParticleBuffer::CloneParticle(copy,data);
            PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
                plan.OutputCorrections[child],copy);
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
            PIC::ParticleTracker::InitParticleID(copy);
            double trackedPosition[3],trackedVelocity[3];
            PIC::ParticleBuffer::GetX(trackedPosition,copy);
            PIC::ParticleBuffer::GetV(trackedVelocity,copy);
            PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(
                trackedPosition,trackedVelocity,species,copy,
                static_cast<void*>(node));
#endif
            local.SplitChildren++;
          }
        }
        ptr=next;
      }
    }
  }

  if (statistics!=NULL) *statistics=local;
  if (emitDiagnostics) {
    long long localValues[8]={local.Examined,local.Kept,
        local.RouletteCandidates,local.RouletteSurvivors,
        local.SplitParents,local.SplitChildren,
        local.RequiredAdditionalParticles,PIC::ParticleBuffer::NAllPart};
    long long globalValues[8]={0,0,0,0,0,0,0,0};
    MPI_Allreduce(localValues,globalValues,8,MPI_LONG_LONG,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread==0 &&
        (globalValues[2]!=0 || globalValues[4]!=0)) {
      std::fprintf(PIC::DiagnospticMessageStream,
          "$PREFIX: individual correction-window transaction iteration=%ld examined=%lld kept=%lld roulette=%lld survivors=%lld split-parents=%lld split-children=%lld required-additional=%lld particles-after=%lld\n",
          completedIterations,globalValues[0],globalValues[1],
          globalValues[2],globalValues[3],globalValues[4],globalValues[5],
          globalValues[6],globalValues[7]);
      std::fflush(PIC::DiagnospticMessageStream);
    }
  }
  return true;
}

void ApplyIndividualCorrectionWindowsOrAbort(long int completedIterations) {
  if (!ApplyIndividualCorrectionWindowsInternal(completedIterations,-1,NULL,
      true)) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "individual correction-window transaction rejected before iteration %ld could complete\n",
          completedIterations);
    std::abort();
  }
}
#endif

double GetCarrierRate(PIC::Mesh::cDataBlockAMR *block,int species,
    int pendingWeightPower=0) {
  const double oldWeight=block->GetLocalParticleWeight(species);
  const double timeStep=block->GetLocalTimeStep(species);
  const double plannedWeight=std::ldexp(oldWeight,pendingWeightPower);
  const double carrierRate=plannedWeight/timeStep;
  if (!(oldWeight>0.0) || !std::isfinite(oldWeight) ||
      !(timeStep>0.0) || !std::isfinite(timeStep) ||
      !(plannedWeight>0.0) || !std::isfinite(plannedWeight) ||
      !(carrierRate>0.0) || !std::isfinite(carrierRate)) {
    std::fprintf(stderr,
        "invalid particle carrier rate while adapting species %d\n",species);
    std::abort();
  }
  return carrierRate;
}

long long MaximumReactionDaughterAllocation(
    PIC::Mesh::cDataBlockAMR *block,
    const std::array<long int,PIC::nTotalSpecies> &counts) {
  if (!Comet::Chemistry::Enabled) return 0;

  // ProcessParticle first selects one mutually exclusive reaction channel.
  // In a fixed-correction build it then residual-resamples every concrete
  // product with expectation
  // (W_parent/dt_parent)/(W_product/dt_product); floor(mu)+Bernoulli(frac)
  // can allocate ceil(mu).  With individual corrections enabled, chemistry
  // instead creates exactly N model daughters per integer-stoichiometric
  // product, where N is the resolved Stage-6 reaction packet count, and
  // stores an exact partition of the carrier conversion in those daughters.
  // The hard allocation bound must mirror those two executors exactly.
  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  std::array<bool,PIC::nTotalSpecies> reachable{};
  for (int species=0;species<PIC::nTotalSpecies;species++)
    // Direct surface parents can react during iteration one, before a model
    // particle of that species happens to be present in this block or before
    // the first sampling break.  Treat the configured source itself as
    // reachability evidence so initialization can prove the hard bound.
    reachable[species]=counts[species]>0 ||
        HasConfiguredDirectSource(species);

  // A daughter absent at the sampling break may be created and react again
  // before the next break.  Close the reaction graph transitively from the
  // species currently present so the diagnostic includes those later
  // generations rather than only direct daughters.
  bool reachabilityChanged;
  do {
    reachabilityChanged=false;
    for (std::size_t reactionIndex=0;
         reactionIndex<network.GetReactionCount();reactionIndex++) {
      if (!network.IsUnitUnimolecular(reactionIndex) ||
          !(network.GetRateCoefficient(reactionIndex,0.0,
              Comet::Chemistry::HeliocentricDistanceAU)>0.0)) continue;
      const int parent=
          network.GetReactantTerm(reactionIndex,0).Species;
      if (!reachable[parent]) continue;
      const ChemicalSolver::cReaction &reaction=
          network.GetReaction(reactionIndex);
      for (std::size_t term=0;term<reaction.ProductCount;term++) {
        const int product=
            network.GetProductTerm(reactionIndex,term).Species;
        if (!reachable[product]) {
          reachable[product]=true;
          reachabilityChanged=true;
        }
      }
    }
  } while (reachabilityChanged);

  long long maximum=0;
  for (std::size_t reactionIndex=0;
       reactionIndex<network.GetReactionCount();reactionIndex++) {
    if (!network.IsUnitUnimolecular(reactionIndex) ||
        !(network.GetRateCoefficient(reactionIndex,0.0,
            Comet::Chemistry::HeliocentricDistanceAU)>0.0)) continue;

    const int parent=
        network.GetReactantTerm(reactionIndex,0).Species;
    if (!reachable[parent]) continue;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
    const double parentCarrier=GetCarrierRate(block,parent);
#endif
    const ChemicalSolver::cReaction &reaction=
        network.GetReaction(reactionIndex);
    long long maximumProducts=0;
    for (std::size_t term=0;term<reaction.ProductCount;term++) {
      const ChemicalSolver::cStoichiometricTerm &product=
          network.GetProductTerm(reactionIndex,term);
      const long long coefficient=static_cast<long long>(
          std::floor(product.Coefficient+0.5));
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      // Stage-6 creates one instance of every concrete product in every
      // complete packet. The reaction-resolved N replaces Stage-4's unit
      // multiplicity while the correction factor carries 1/N of the event.
      const long long oneProductMaximum=
          Comet::Chemistry::GetDaughterPacketCount(reactionIndex);
#else
      const double mean=parentCarrier/
          GetCarrierRate(block,product.Species);
      const long long oneProductMaximum=
          static_cast<long long>(std::ceil(mean));
#endif
      if (oneProductMaximum<=0 || coefficient<=0 ||
          oneProductMaximum>LLONG_MAX/coefficient ||
          maximumProducts>LLONG_MAX-oneProductMaximum*coefficient) {
        std::fprintf(stderr,
            "reaction %s hard daughter allocation overflows long long\n",
            reaction.Id);
        std::abort();
      }
      maximumProducts+=oneProductMaximum*coefficient;
    }
    maximum=std::max(maximum,maximumProducts);
  }
  return maximum;
}

long int RequiredCloneCount(const cWeightAdaptationPlan &plan) {
  if (plan.WeightPower>=0) return 0;
  long int multiplier=1;
  for (int step=0;step < -plan.WeightPower;step++) {
    if (multiplier>LONG_MAX/2) {
      std::fprintf(stderr,
          "particle-weight adaptation clone multiplier overflows long int\n");
      std::abort();
    }
    multiplier*=2;
  }
  if (plan.ParticleCount>0 &&
      multiplier-1>LONG_MAX/plan.ParticleCount) {
    std::fprintf(stderr,
        "particle-weight adaptation clone count overflows long int\n");
    std::abort();
  }
  return plan.ParticleCount*(multiplier-1);
}

void SplitSpeciesByTwo(const cWeightAdaptationPlan &plan) {
  PIC::Mesh::cDataBlockAMR *block=plan.Node->block;
  const double oldWeight=block->GetLocalParticleWeight(plan.Species);
  const double newWeight=0.5*oldWeight;
  if (!(oldWeight>0.0) || !std::isfinite(oldWeight) ||
      !(newWeight>0.0) || !std::isfinite(newWeight)) {
    std::fprintf(stderr,
        "invalid particle weight while splitting species %d\n",plan.Species);
    std::abort();
  }

  // Snapshot the next pointer before inserting the clone at the cell-list
  // head.  Thus every particle present at the beginning of this factor-two
  // step is cloned exactly once, while the new copies are not visited until a
  // possible later factor-two step.  Position, velocity, internal energy,
  // and all other physical state are copied unchanged.
  for (long int cell=0;cell<BlockCellCount;cell++) {
    long int &first=block->FirstCellParticleTable[cell];
    for (long int ptr=first;ptr!=-1;) {
      const long int next=PIC::ParticleBuffer::GetNext(ptr);
      PIC::ParticleBuffer::byte *source=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      if (PIC::ParticleBuffer::GetI(source)==plan.Species) {
        const long int clone=PIC::ParticleBuffer::GetNewParticle(first);
        PIC::ParticleBuffer::byte *copy=
            PIC::ParticleBuffer::GetParticleDataPointer(clone);
        PIC::ParticleBuffer::CloneParticle(copy,source);
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
        PIC::ParticleTracker::InitParticleID(copy);
#endif
      }
      ptr=next;
    }
  }
  block->SetLocalParticleWeight(newWeight,plan.Species);
}

void RouletteSpeciesByTwo(const cWeightAdaptationPlan &plan) {
  PIC::Mesh::cDataBlockAMR *block=plan.Node->block;
  const double oldWeight=block->GetLocalParticleWeight(plan.Species);
  const double newWeight=2.0*oldWeight;
  if (!(oldWeight>0.0) || !std::isfinite(oldWeight) ||
      !(newWeight>0.0) || !std::isfinite(newWeight)) {
    std::fprintf(stderr,
        "invalid particle weight while reducing species %d\n",plan.Species);
    std::abort();
  }

  // Each survivor represents the new, larger weight.  The requested
  // probability is W_old/W_new (exactly one half for this step), so the
  // represented number, mass, momentum, energy, and chemical population are
  // conserved in expectation without an individual weight correction.
  const double keepProbability=oldWeight/newWeight;
  for (long int cell=0;cell<BlockCellCount;cell++) {
    long int &first=block->FirstCellParticleTable[cell];
    for (long int ptr=first;ptr!=-1;) {
      const long int next=PIC::ParticleBuffer::GetNext(ptr);
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      if (PIC::ParticleBuffer::GetI(data)==plan.Species &&
          rnd()>=keepProbability)
        PIC::ParticleBuffer::DeleteParticle(ptr,first);
      ptr=next;
    }
  }
  block->SetLocalParticleWeight(newWeight,plan.Species);
}

void RefreshGlobalMinimumParticleWeights() {
  // SetLocalParticleWeight() only lowers this legacy global-minimum cache.
  // Recompute it explicitly whenever an adaptation raises block weights.
  std::array<double,PIC::nTotalSpecies> localMinimum;
  std::array<double,PIC::nTotalSpecies> globalMinimum;
  localMinimum.fill(std::numeric_limits<double>::max());
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    PIC::Mesh::cDataBlockAMR *block=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex]->block;
    if (block==NULL) continue;
    for (int species=0;species<PIC::nTotalSpecies;species++)
      localMinimum[species]=std::min(localMinimum[species],
          block->GetLocalParticleWeight(species));
  }
  MPI_Allreduce(localMinimum.data(),globalMinimum.data(),PIC::nTotalSpecies,
      MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!(globalMinimum[species]>0.0) ||
        !std::isfinite(globalMinimum[species])) {
      std::fprintf(stderr,
          "invalid global minimum after adapting species %d\n",species);
      std::abort();
    }
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]=
        globalMinimum[species];
  }
}

void SynchronizeAdaptedBlockWeights() {
  // The owner block is the only copy whose live particle list is resampled.
  // AMPS boundary-layer blocks also store local weights because the face
  // mover needs both sides of an MPI crossing.  The standard block-data
  // exchange propagates the owner's updated associated data to those copies.
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  PIC::Mesh::mesh->ParallelBlockDataExchange();
  RefreshGlobalMinimumParticleWeights();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

struct cFaceCarrierPropagationStatistics {
  long long LocalChangedBlockSpecies=0;
  long long GlobalChangedBlockSpecies=0;
  long long ActiveComponents=0;
  long long MaximumPropagationHops=0;
  double MaximumCloneMeanBefore=0.0;
  double MaximumCloneMeanAfter=0.0;
};

cFaceCarrierPropagationStatistics LastFaceCarrierPropagationStatistics;

void CollectCalculationLeafNodes(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,
    std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> &nodes) {
  if (node==NULL) return;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->IsUsedInCalculationFlag) nodes.push_back(node);
    return;
  }

  for (int child=0;child<(1<<DIM);child++)
    CollectCalculationLeafNodes(node->downNode[child],nodes);
}

// Stage 2 plans the complete component transaction before it changes either
// a particle list or a block weight.  The old implementation first committed
// independent block splits and only then raised carriers across faces; one
// empty/coarse block could consequently roulette away the daughter split.
// This record owns every piece of state required to preflight, commit, and
// verify one collective transaction without re-planning from mutated data.
struct cComponentCarrierTransaction {
  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> Nodes;
  std::vector<std::vector<std::size_t> > Neighbors;
  std::vector<long long> GlobalCounts;
  std::vector<double> GlobalCarriers;
  std::vector<double> TargetCarriers;
  std::vector<std::size_t> ComponentIndex;
  std::vector<std::vector<std::size_t> > ComponentNodes;
  std::vector<long long> ComponentParticleCount;
  std::vector<long long> ComponentCellCount;
  std::vector<double> ComponentCarrier;
  std::vector<int> ComponentPowerAdjustment;
  std::vector<unsigned char> ChemistryConstrained;
  std::vector<unsigned char> Reachable;
  std::vector<cWeightAdaptationPlan> LocalPlans;
  std::vector<cLockedWeightSnapshot> LockedWeights;
  cFaceCarrierPropagationStatistics PropagationStatistics;
  long long LocalRequiredClones=0;
  long long LocalCapacityAvailable=0;
  long long GlobalMaximumTotalBlockParticles=0;
  long long GlobalMaximumProductsBefore=0;
  long long GlobalMaximumProductsPlanned=0;
  bool CapacityAccepted=false;
};

long long CheckedReactionCoefficient(double coefficient,const char *reaction) {
  if (!(coefficient>0.0) || !std::isfinite(coefficient)) {
    std::fprintf(stderr,"reaction %s has an invalid product coefficient\n",
        reaction);
    std::abort();
  }
  const long long rounded=static_cast<long long>(std::llround(coefficient));
  if (rounded<=0 || std::fabs(coefficient-rounded)>
      64.0*std::numeric_limits<double>::epsilon()*
          std::max(std::fabs(coefficient),1.0)) {
    std::fprintf(stderr,
        "reaction %s has a non-integral product coefficient\n",reaction);
    std::abort();
  }
  return rounded;
}

cComponentCarrierTransaction BuildComponentCarrierTransaction(
    long long capacityOverride=-1) {
  namespace Planner=
      Comet::RadialSimulation::ComponentCarrierPlanner;
  cComponentCarrierTransaction transaction;

  // The replicated AMR tree provides the same stable node order on every
  // rank.  Only node->Thread owns authoritative particles/weights; collective
  // sums below turn those single-owner values into an identical planning
  // table on all ranks.
  CollectCalculationLeafNodes(PIC::Mesh::mesh->rootTree,transaction.Nodes);
  if (transaction.Nodes.empty()) {
    std::fprintf(stderr,
        "component carrier planner found no calculation blocks\n");
    std::abort();
  }
  if (transaction.Nodes.size()*PIC::nTotalSpecies>
      static_cast<std::size_t>(INT_MAX)) {
    std::fprintf(stderr,
        "component carrier planner exceeds the MPI count limit\n");
    std::abort();
  }

  std::unordered_map<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,std::size_t>
      nodeIndex;
  nodeIndex.reserve(2*transaction.Nodes.size());
  for (std::size_t block=0;block<transaction.Nodes.size();block++)
    nodeIndex[transaction.Nodes[block]]=block;

  constexpr int nFaces=2*DIM;
  constexpr int nFaceSlotsI=(DIM>=2) ? 2 : 1;
  constexpr int nFaceSlotsJ=(DIM>=3) ? 2 : 1;
  std::vector<Planner::cEdge> edges;
  transaction.Neighbors.resize(transaction.Nodes.size());
  for (std::size_t block=0;block<transaction.Nodes.size();block++) {
    for (int face=0;face<nFaces;face++)
      for (int iFace=0;iFace<nFaceSlotsI;iFace++)
        for (int jFace=0;jFace<nFaceSlotsJ;jFace++) {
          cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neighbor=
              transaction.Nodes[block]->GetNeibFace(face,iFace,jFace,
                  PIC::Mesh::mesh);
          if (neighbor==NULL) continue;
          const auto location=nodeIndex.find(neighbor);
          if (location==nodeIndex.end()) {
            std::fprintf(stderr,
                "component carrier planner found a face neighbor outside the leaf table\n");
            std::abort();
          }
          const std::size_t other=location->second;
          if (other==block) continue;
          edges.push_back(Planner::cEdge(block,other));
          transaction.Neighbors[block].push_back(other);
          transaction.Neighbors[other].push_back(block);
        }
  }
  for (std::vector<std::size_t> &list:transaction.Neighbors) {
    std::sort(list.begin(),list.end());
    list.erase(std::unique(list.begin(),list.end()),list.end());
  }

  const std::size_t valueCount=
      transaction.Nodes.size()*PIC::nTotalSpecies;
  std::vector<long long> localCounts(valueCount,0);
  transaction.GlobalCounts.assign(valueCount,0);
  std::vector<double> localCarriers(valueCount,0.0);
  transaction.GlobalCarriers.assign(valueCount,0.0);
  std::vector<int> localOwners(transaction.Nodes.size(),0),
      globalOwners(transaction.Nodes.size(),0);
  long long localMaximumProductsBefore=0;
  for (std::size_t block=0;block<transaction.Nodes.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        transaction.Nodes[block];
    if (node->Thread!=PIC::ThisThread) continue;
    if (node->block==NULL) {
      std::fprintf(stderr,
          "component carrier owner block is not allocated\n");
      std::abort();
    }
    localOwners[block]=1;
    const std::array<long int,PIC::nTotalSpecies> counts=
        CountAllSpeciesInBlock(node->block);
    localMaximumProductsBefore=std::max(localMaximumProductsBefore,
        MaximumReactionDaughterAllocation(node->block,counts));
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      localCounts[offset]=counts[species];
      localCarriers[offset]=GetCarrierRate(node->block,species);
    }
  }
  MPI_Allreduce(localOwners.data(),globalOwners.data(),
      static_cast<int>(transaction.Nodes.size()),MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localCounts.data(),transaction.GlobalCounts.data(),
      static_cast<int>(valueCount),MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localCarriers.data(),transaction.GlobalCarriers.data(),
      static_cast<int>(valueCount),MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localMaximumProductsBefore,
      &transaction.GlobalMaximumProductsBefore,1,MPI_LONG_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);

  for (std::size_t block=0;block<transaction.Nodes.size();block++) {
    if (globalOwners[block]!=1) {
      std::fprintf(stderr,
          "component carrier block has %d authoritative owners\n",
          globalOwners[block]);
      std::abort();
    }
    long long total=0;
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      if (transaction.GlobalCounts[offset]<0 ||
          !(transaction.GlobalCarriers[offset]>0.0) ||
          !std::isfinite(transaction.GlobalCarriers[offset]) ||
          total>LLONG_MAX-transaction.GlobalCounts[offset]) {
        std::fprintf(stderr,
            "invalid authoritative component state for block %zu species %d\n",
            block,species);
        std::abort();
      }
      total+=transaction.GlobalCounts[offset];
    }
    transaction.GlobalMaximumTotalBlockParticles=std::max(
        transaction.GlobalMaximumTotalBlockParticles,total);
  }

  std::vector<Planner::cPlan> speciesPlans(PIC::nTotalSpecies);
  transaction.TargetCarriers=transaction.GlobalCarriers;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    std::vector<Planner::cNodeState> state;
    state.reserve(transaction.Nodes.size());
    for (std::size_t block=0;block<transaction.Nodes.size();block++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      state.push_back(Planner::cNodeState(transaction.GlobalCounts[offset],
          BlockCellCount,transaction.GlobalCarriers[offset]));
    }

    try {
      speciesPlans[species]=Planner::BuildPlan(state,edges,Planner::cPolicy(
          Comet::RadialSimulation::SpeciesMinimumParticlesPerCell[species],
          Comet::RadialSimulation::SpeciesMaximumParticlesPerCell[species],
          Comet::RadialSimulation::ParticleWeightLocked[species]!=0));
    }
    catch (const std::exception &exception) {
      std::fprintf(stderr,
          "component carrier planning failed for species %d: %s\n",species,
          exception.what());
      std::abort();
    }

    const Planner::cPlan &plan=speciesPlans[species];
    if (species==0) {
      transaction.ComponentIndex=plan.ComponentIndex;
      transaction.ComponentNodes.reserve(plan.Components.size());
      for (const Planner::cComponentPlan &component:plan.Components)
        transaction.ComponentNodes.push_back(component.Nodes);
    }
    else if (plan.ComponentIndex!=transaction.ComponentIndex ||
        plan.Components.size()!=transaction.ComponentNodes.size()) {
      std::fprintf(stderr,
          "component topology changed between species plans\n");
      std::abort();
    }

    for (std::size_t component=0;component<plan.Components.size();component++) {
      if (plan.Components[component].Nodes!=
          transaction.ComponentNodes[component]) {
        std::fprintf(stderr,
            "component node ordering changed between species plans\n");
        std::abort();
      }
      if (plan.Components[component].Nodes.size()>1)
        transaction.PropagationStatistics.ActiveComponents++;
    }
    for (std::size_t block=0;block<transaction.Nodes.size();block++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      const int power=plan.WeightPower[block];
      const double target=std::ldexp(transaction.GlobalCarriers[offset],power);
      if (!(target>0.0) || !std::isfinite(target)) {
        std::fprintf(stderr,
            "component carrier target is invalid for block %zu species %d\n",
            block,species);
        std::abort();
      }
      transaction.TargetCarriers[offset]=target;
    }
  }

  const std::size_t componentCount=transaction.ComponentNodes.size();
  const std::size_t componentValueCount=componentCount*PIC::nTotalSpecies;
  transaction.ComponentParticleCount.assign(componentValueCount,0);
  transaction.ComponentCellCount.assign(componentCount,0);
  transaction.ComponentCarrier.assign(componentValueCount,0.0);
  transaction.ComponentPowerAdjustment.assign(componentValueCount,0);
  transaction.ChemistryConstrained.assign(componentValueCount,0);
  transaction.Reachable.assign(componentValueCount,0);
  for (std::size_t component=0;component<componentCount;component++) {
    if (transaction.ComponentNodes[component].size()>
        static_cast<std::size_t>(LLONG_MAX/BlockCellCount)) {
      std::fprintf(stderr,"component cell count overflows long long\n");
      std::abort();
    }
    transaction.ComponentCellCount[component]=
        static_cast<long long>(transaction.ComponentNodes[component].size())*
            BlockCellCount;
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t componentOffset=
          component*PIC::nTotalSpecies+species;
      long long count=0;
      for (std::size_t block:transaction.ComponentNodes[component]) {
        const long long addition=transaction.GlobalCounts[
            block*PIC::nTotalSpecies+species];
        if (count>LLONG_MAX-addition) {
          std::fprintf(stderr,
              "component particle count overflows long long\n");
          std::abort();
        }
        count+=addition;
      }
      transaction.ComponentParticleCount[componentOffset]=count;
      transaction.ComponentCarrier[componentOffset]=
          speciesPlans[species].Components[component].TargetCarrier;
      transaction.Reachable[componentOffset]=count>0 ||
          HasConfiguredDirectSource(species);
    }
  }

  // Close chemistry reachability on the same component scalars selected by
  // the population planner and enforce the irreducible stoichiometric bound.
  // In correction-off builds, also solve the residual-allocation cap by
  // changing base carriers.  In correction-on builds each concrete product is
  // represented by one daughter whose c stores the carrier conversion, so a
  // base-carrier constraint here would be both unnecessary and incorrect.
  if (Comet::Chemistry::Enabled) {
    const ChemicalSolver::cReactionNetwork &network=
        Comet::Chemistry::GetNetwork();
    for (std::size_t component=0;component<componentCount;component++) {
      bool changed=true;
      std::size_t iterations=0;
      const std::size_t maximumIterations=
          4096*PIC::nTotalSpecies*(network.GetReactionCount()+1)+1;
      while (changed) {
        changed=false;
        if (++iterations>maximumIterations) {
          std::fprintf(stderr,
              "component chemistry carrier constraints did not converge\n");
          std::abort();
        }
        for (std::size_t reactionIndex=0;
             reactionIndex<network.GetReactionCount();reactionIndex++) {
          if (!network.IsUnitUnimolecular(reactionIndex) ||
              !(network.GetRateCoefficient(reactionIndex,0.0,
                  Comet::Chemistry::HeliocentricDistanceAU)>0.0)) continue;
          const ChemicalSolver::cReaction &reaction=
              network.GetReaction(reactionIndex);
          const int parent=
              network.GetReactantTerm(reactionIndex,0).Species;
          const std::size_t parentOffset=
              component*PIC::nTotalSpecies+parent;
          if (!transaction.Reachable[parentOffset]) continue;

          std::array<long long,PIC::nTotalSpecies> coefficient{};
          long long irreducibleFloor=0;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
          const long long reactionPacketCount=
              Comet::Chemistry::GetDaughterPacketCount(reactionIndex);
#endif
          for (std::size_t term=0;term<reaction.ProductCount;term++) {
            const ChemicalSolver::cStoichiometricTerm &product=
                network.GetProductTerm(reactionIndex,term);
            const long long value=
                CheckedReactionCoefficient(product.Coefficient,reaction.Id);
            long long irreducibleContribution=value;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
            if (reactionPacketCount<=0 ||
                value>LLONG_MAX/reactionPacketCount) {
              std::fprintf(stderr,
                  "reaction %s packet product count overflows\n",
                  reaction.Id);
              std::abort();
            }
            irreducibleContribution*=reactionPacketCount;
#endif
            if (coefficient[product.Species]>LLONG_MAX-value ||
                irreducibleFloor>LLONG_MAX-irreducibleContribution) {
              std::fprintf(stderr,
                  "reaction %s product coefficient sum overflows\n",
                  reaction.Id);
              std::abort();
            }
            coefficient[product.Species]+=value;
            irreducibleFloor+=irreducibleContribution;
            const std::size_t productOffset=
                component*PIC::nTotalSpecies+product.Species;
            if (!transaction.Reachable[productOffset]) {
              transaction.Reachable[productOffset]=1;
              changed=true;
            }
          }
          if (irreducibleFloor>
              Comet::RadialSimulation::
                  MaximumPhotochemicalProductsPerReaction) {
            std::fprintf(stderr,
                "reaction %s has irreducible product count %lld above limit %ld\n",
                reaction.Id,irreducibleFloor,
                Comet::RadialSimulation::
                    MaximumPhotochemicalProductsPerReaction);
            std::abort();
          }

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
          // Reachability changes still drive another closure pass, but no
          // component carrier is changed for chemistry in this mode.
          continue;
#else
          auto hardAllocation = [&]() -> long long {
            long long total=0;
            for (int product=0;product<PIC::nTotalSpecies;product++) {
              if (coefficient[product]==0) continue;
              const double ratio=
                  transaction.ComponentCarrier[parentOffset]/
                  transaction.ComponentCarrier[
                      component*PIC::nTotalSpecies+product];
              const long long multiplicity=
                  static_cast<long long>(std::ceil(ratio));
              if (multiplicity<=0 ||
                  multiplicity>LLONG_MAX/coefficient[product] ||
                  total>LLONG_MAX-multiplicity*coefficient[product]) {
                std::fprintf(stderr,
                    "reaction %s component allocation overflows\n",
                    reaction.Id);
                std::abort();
              }
              total+=multiplicity*coefficient[product];
            }
            return total;
          };

          while (hardAllocation()>
              Comet::RadialSimulation::
                  MaximumPhotochemicalProductsPerReaction) {
            int selectedProduct=-1;
            long long selectedContribution=-1;
            for (int product=0;product<PIC::nTotalSpecies;product++) {
              if (coefficient[product]==0 || product==parent ||
                  Comet::RadialSimulation::ParticleWeightLocked[product])
                continue;
              const std::size_t productOffset=
                  component*PIC::nTotalSpecies+product;
              const long long multiplicity=static_cast<long long>(std::ceil(
                  transaction.ComponentCarrier[parentOffset]/
                  transaction.ComponentCarrier[productOffset]));
              if (multiplicity<=1) continue;
              const long long contribution=
                  coefficient[product]*multiplicity;
              if (contribution>selectedContribution) {
                selectedContribution=contribution;
                selectedProduct=product;
              }
            }

            if (selectedProduct>=0) {
              const std::size_t productOffset=
                  component*PIC::nTotalSpecies+selectedProduct;
              if (transaction.ComponentPowerAdjustment[productOffset]>=1024 ||
                  !std::isfinite(
                      2.0*transaction.ComponentCarrier[productOffset])) {
                std::fprintf(stderr,
                    "reaction %s requires an invalid daughter carrier\n",
                    reaction.Id);
                std::abort();
              }
              transaction.ComponentPowerAdjustment[productOffset]++;
              transaction.ComponentCarrier[productOffset]*=2.0;
              transaction.ChemistryConstrained[productOffset]=1;
            }
            else if (!Comet::RadialSimulation::
                ParticleWeightLocked[parent]) {
              if (transaction.ComponentPowerAdjustment[parentOffset]<=-1024) {
                std::fprintf(stderr,
                    "reaction %s requires an invalid parent carrier\n",
                    reaction.Id);
                std::abort();
              }
              transaction.ComponentPowerAdjustment[parentOffset]--;
              transaction.ComponentCarrier[parentOffset]*=0.5;
              if (!(transaction.ComponentCarrier[parentOffset]>0.0) ||
                  !std::isfinite(
                      transaction.ComponentCarrier[parentOffset])) {
                std::fprintf(stderr,
                    "reaction %s underflowed its parent carrier\n",
                    reaction.Id);
                std::abort();
              }
              transaction.ChemistryConstrained[parentOffset]=1;
            }
            else {
              std::fprintf(stderr,
                  "reaction %s cannot meet the component allocation limit without changing locked carriers\n",
                  reaction.Id);
              std::abort();
            }
            changed=true;
          }
#endif
        }
      }
    }
  }

  // Convert the final component scalars back into per-owner dyadic powers.
  // No state has been mutated yet.  Exact clone requirements therefore cover
  // every population and chemistry split in the transaction, not merely the
  // initial per-block request as in the old two-phase controller.
  for (std::size_t block=0;block<transaction.Nodes.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=transaction.Nodes[block];
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      const std::size_t component=transaction.ComponentIndex[block];
      const std::size_t componentOffset=
          component*PIC::nTotalSpecies+species;
      const long long combinedPower=static_cast<long long>(
          speciesPlans[species].WeightPower[block])+
          transaction.ComponentPowerAdjustment[componentOffset];
      if (combinedPower<-1024 || combinedPower>1024) {
        std::fprintf(stderr,
            "final component power is out of range for species %d\n",species);
        std::abort();
      }
      const int power=static_cast<int>(combinedPower);
      const double target=
          std::ldexp(transaction.GlobalCarriers[offset],power);
      if (!(target>0.0) || !std::isfinite(target)) {
        std::fprintf(stderr,
            "final component carrier is invalid for species %d\n",species);
        std::abort();
      }
      const double expected=transaction.ComponentCarrier[componentOffset];
      const double scale=std::max(std::fabs(target),std::fabs(expected));
      if (std::fabs(target-expected)>
          1024.0*std::numeric_limits<double>::epsilon()*scale) {
        std::fprintf(stderr,
            "final carrier left the component dyadic lattice for species %d\n",
            species);
        std::abort();
      }
      transaction.TargetCarriers[offset]=target;

      if (Comet::RadialSimulation::ParticleWeightLocked[species]) {
        if (power!=0) {
          std::fprintf(stderr,
              "component plan would change locked species %d\n",species);
          std::abort();
        }
        if (node->Thread==PIC::ThisThread) {
          cLockedWeightSnapshot snapshot;
          snapshot.Block=node->block;
          snapshot.Species=species;
          snapshot.Weight=node->block->GetLocalParticleWeight(species);
          snapshot.WeightBits=DoubleBits(snapshot.Weight);
          transaction.LockedWeights.push_back(snapshot);
        }
        continue;
      }
      if (power==0) continue;
      transaction.PropagationStatistics.GlobalChangedBlockSpecies++;
      if (node->Thread!=PIC::ThisThread) continue;
      if (transaction.GlobalCounts[offset]>LONG_MAX) {
        std::fprintf(stderr,
            "local component population exceeds long int\n");
        std::abort();
      }
      cWeightAdaptationPlan plan;
      plan.Node=node;
      plan.Species=species;
      plan.ParticleCount=
          static_cast<long int>(transaction.GlobalCounts[offset]);
      plan.WeightPower=power;
      plan.EmptyPhotochemicalDaughter=plan.ParticleCount==0 &&
          transaction.ChemistryConstrained[componentOffset] && power>0;
      plan.PhotochemicalCarrierConstraint=
          transaction.ChemistryConstrained[componentOffset];
      const long int clones=RequiredCloneCount(plan);
      if (transaction.LocalRequiredClones>LLONG_MAX-clones) {
        std::fprintf(stderr,
            "component transaction clone count overflows long long\n");
        std::abort();
      }
      transaction.LocalRequiredClones+=clones;
      transaction.LocalPlans.push_back(plan);
    }
  }

  // Quantify both directed senses of every face before and after the plan.
  // Equality is a plan invariant; the mover sees no systematic cloning even
  // when a daughter carrier was lowered through hundreds of empty blocks.
  for (std::size_t source=0;source<transaction.Nodes.size();source++)
    for (std::size_t destination:transaction.Neighbors[source])
      for (int species=0;species<PIC::nTotalSpecies;species++) {
        const std::size_t sourceOffset=
            source*PIC::nTotalSpecies+species;
        const std::size_t destinationOffset=
            destination*PIC::nTotalSpecies+species;
        transaction.PropagationStatistics.MaximumCloneMeanBefore=std::max(
            transaction.PropagationStatistics.MaximumCloneMeanBefore,
            transaction.GlobalCarriers[sourceOffset]/
                transaction.GlobalCarriers[destinationOffset]);
        transaction.PropagationStatistics.MaximumCloneMeanAfter=std::max(
            transaction.PropagationStatistics.MaximumCloneMeanAfter,
            transaction.TargetCarriers[sourceOffset]/
                transaction.TargetCarriers[destinationOffset]);
      }
  const double faceTolerance=
      Comet::RadialSimulation::MaximumFaceTransferMultiplicity*
      (1.0+256.0*std::numeric_limits<double>::epsilon());
  if (transaction.PropagationStatistics.MaximumCloneMeanAfter>
      faceTolerance) {
    std::fprintf(stderr,
        "component plan leaves face multiplicity %.16e above limit %ld\n",
        transaction.PropagationStatistics.MaximumCloneMeanAfter,
        Comet::RadialSimulation::MaximumFaceTransferMultiplicity);
    std::abort();
  }

  // Record the maximum graph distance once per spatial component.  This is a
  // validation diagnostic proving that a plan reached beyond immediate
  // neighbors; it does not affect the carrier decision.
  for (const std::vector<std::size_t> &component:transaction.ComponentNodes) {
    if (component.empty()) continue;
    std::vector<long long> distance(transaction.Nodes.size(),-1);
    std::vector<std::size_t> work(1,component.front());
    distance[component.front()]=0;
    for (std::size_t cursor=0;cursor<work.size();cursor++) {
      const std::size_t block=work[cursor];
      transaction.PropagationStatistics.MaximumPropagationHops=std::max(
          transaction.PropagationStatistics.MaximumPropagationHops,
          distance[block]);
      for (std::size_t neighbor:transaction.Neighbors[block])
        if (distance[neighbor]<0) {
          distance[neighbor]=distance[block]+1;
          work.push_back(neighbor);
        }
    }
  }

  // The capacity ceiling preserves the input reserve for chemistry and MPI
  // bookkeeping.  Roulette releases are deliberately ignored: accepting a
  // split may not depend on stochastic deletions that have not happened yet.
  const long double ceilingValue=std::floor(
      static_cast<long double>(PIC::ParticleBuffer::MaxNPart)*
      (1.0L-Comet::RadialSimulation::AllocationReserveFraction));
  if (ceilingValue<0.0L || ceilingValue>LLONG_MAX) {
    std::fprintf(stderr,"component transaction capacity is invalid\n");
    std::abort();
  }
  const long long ceiling=static_cast<long long>(ceilingValue);
  transaction.LocalCapacityAvailable=ceiling>
      PIC::ParticleBuffer::NAllPart ?
      ceiling-PIC::ParticleBuffer::NAllPart : 0;
  if (capacityOverride>=0)
    transaction.LocalCapacityAvailable=capacityOverride;
  const int localCapacityFailure=
      transaction.LocalRequiredClones>
          transaction.LocalCapacityAvailable ? 1 : 0;
  int globalCapacityFailure=0;
  MPI_Allreduce(&localCapacityFailure,&globalCapacityFailure,1,MPI_INT,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  transaction.CapacityAccepted=globalCapacityFailure==0;

  // Audit every final reachable reaction against the same component values
  // that will be committed.  This is separate from the mutation-free solver
  // loop so a future network change cannot accidentally skip a channel.
  if (Comet::Chemistry::Enabled) {
    const ChemicalSolver::cReactionNetwork &network=
        Comet::Chemistry::GetNetwork();
    for (std::size_t component=0;component<componentCount;component++)
      for (std::size_t reactionIndex=0;
           reactionIndex<network.GetReactionCount();reactionIndex++) {
        if (!network.IsUnitUnimolecular(reactionIndex) ||
            !(network.GetRateCoefficient(reactionIndex,0.0,
                Comet::Chemistry::HeliocentricDistanceAU)>0.0)) continue;
        const int parent=
            network.GetReactantTerm(reactionIndex,0).Species;
        if (!transaction.Reachable[
            component*PIC::nTotalSpecies+parent]) continue;
        const ChemicalSolver::cReaction &reaction=
            network.GetReaction(reactionIndex);
        long long total=0;
        for (std::size_t term=0;term<reaction.ProductCount;term++) {
          const ChemicalSolver::cStoichiometricTerm &product=
              network.GetProductTerm(reactionIndex,term);
          const long long coefficient=
              CheckedReactionCoefficient(product.Coefficient,reaction.Id);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
          // Packet splitting still avoids the legacy base-carrier ratio: it
          // creates exactly N corrected instances of each concrete product.
          // Keep this audit identical to MaximumReactionDaughterAllocation()
          // and ProcessParticle().
          const long long multiplicity=
              Comet::Chemistry::GetDaughterPacketCount(reactionIndex);
#else
          const long long multiplicity=static_cast<long long>(std::ceil(
              transaction.ComponentCarrier[
                  component*PIC::nTotalSpecies+parent]/
              transaction.ComponentCarrier[
                  component*PIC::nTotalSpecies+product.Species]));
#endif
          if (multiplicity<=0 || multiplicity>LLONG_MAX/coefficient ||
              total>LLONG_MAX-multiplicity*coefficient) {
            std::fprintf(stderr,
                "reaction %s final allocation overflows\n",reaction.Id);
            std::abort();
          }
          total+=multiplicity*coefficient;
        }
        transaction.GlobalMaximumProductsPlanned=std::max(
            transaction.GlobalMaximumProductsPlanned,total);
      }
  }
  if (transaction.GlobalMaximumProductsPlanned>
      Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction) {
    std::fprintf(stderr,
        "component chemistry plan leaves %lld products above limit %ld\n",
        transaction.GlobalMaximumProductsPlanned,
        Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction);
    std::abort();
  }
  return transaction;
}

void SetEmptySpeciesWeight(const cWeightAdaptationPlan &plan) {
  if (plan.ParticleCount!=0 ||
      CountSpeciesInBlock(plan.Node->block,plan.Species)!=0) {
    std::fprintf(stderr,
        "component empty-weight plan found a represented particle\n");
    std::abort();
  }
  const double oldWeight=
      plan.Node->block->GetLocalParticleWeight(plan.Species);
  const double newWeight=std::ldexp(oldWeight,plan.WeightPower);
  if (!(oldWeight>0.0) || !std::isfinite(oldWeight) ||
      !(newWeight>0.0) || !std::isfinite(newWeight)) {
    std::fprintf(stderr,
        "invalid empty component weight for species %d\n",plan.Species);
    std::abort();
  }
  // No list is created.  The empty block merely receives the same W/dt
  // lattice value as every populated member of its face component.
  plan.Node->block->SetLocalParticleWeight(newWeight,plan.Species);
}

void EmergencyRouletteAllAllocatedBlockCopies(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return;
  if (node->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    for (int child=0;child<(1<<DIM);child++)
      EmergencyRouletteAllAllocatedBlockCopies(node->downNode[child]);
    return;
  }
  if (node->block==NULL) return;

  PIC::Mesh::cDataBlockAMR *block=node->block;
  for (long int cell=0;cell<BlockCellCount;cell++) {
    long int &first=block->FirstCellParticleTable[cell];
    for (long int ptr=first;ptr!=-1;) {
      const long int next=PIC::ParticleBuffer::GetNext(ptr);

      // This emergency operation may run after the mover has staged a
      // particle in a non-owner boundary-block copy and before MPI migration.
      // Deletion is safe for both completed and continuation particles; every
      // survivor retains its complete state, including any remaining local-
      // timestep travel fraction.
      if (rnd()>=0.5) PIC::ParticleBuffer::DeleteParticle(ptr,first);
      ptr=next;
    }
  }

  // Every allocated owner or boundary copy is changed by the same factor.
  // Therefore a particle staged for MPI transfer sees the same new weight on
  // the sending and receiving copies, and no per-particle correction is
  // needed.  Applying the factor to every species also preserves all W/dt
  // ratios used by chemistry and face-crossing flux conversion.
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    const double oldWeight=block->GetLocalParticleWeight(species);
    const double newWeight=2.0*oldWeight;
    if (!(oldWeight>0.0) || !std::isfinite(oldWeight) ||
        !(newWeight>0.0) || !std::isfinite(newWeight)) {
      std::fprintf(stderr,
          "invalid emergency particle weight for species %d\n",species);
      std::abort();
    }
    block->SetLocalParticleWeight(newWeight,species);
  }
}

void AdaptParticleWeights(long int completedIterations) {
  if (!Comet::RadialSimulation::ParticleWeightAdaptationEnabled) return;

  cComponentCarrierTransaction transaction=
      BuildComponentCarrierTransaction();

  // Preserve the existing independent all-species emergency trigger.  It is
  // evaluated from the same authoritative pre-transaction counts used by the
  // component planner and still executes only at the next safe pre-exchange
  // callback.
  const long long samplingEmergencyBlockThreshold=
      GetSamplingPopulationEmergencyBlockThreshold();
  if (transaction.GlobalMaximumTotalBlockParticles>
      samplingEmergencyBlockThreshold) {
    if (!SamplingPopulationEmergencyPending)
      Comet::Diagnostics::RecordEmergencyRequest();
    SamplingPopulationEmergencyPending=true;
    SamplingPopulationEmergencyMaximumBlockParticles=
        transaction.GlobalMaximumTotalBlockParticles;
    SamplingPopulationEmergencyThresholdBlockParticles=
        samplingEmergencyBlockThreshold;
    SamplingPopulationEmergencyIteration=completedIterations;
    if (PIC::ThisThread==0) {
      std::fprintf(PIC::DiagnospticMessageStream,
          "$PREFIX: radial sampling-population emergency requested iteration=%ld max-mean-total-particles-per-cell=%.6e threshold=%.6e multiplier=%lld species=%d execution=before-next-regular-mpi-particle-exchange\n",
          completedIterations,
          static_cast<double>(
              transaction.GlobalMaximumTotalBlockParticles)/BlockCellCount,
          static_cast<double>(samplingEmergencyBlockThreshold)/
              BlockCellCount,
          SamplingPopulationEmergencyMultiplier,PIC::nTotalSpecies);
      std::fflush(PIC::DiagnospticMessageStream);
    }
  }

  if (!transaction.CapacityAccepted) {
    if (transaction.LocalRequiredClones>
        transaction.LocalCapacityAvailable)
      std::fprintf(stderr,
          "rank %d component transaction needs %lld clones but has %lld reserved particle slots; no state was changed\n",
          PIC::ThisThread,transaction.LocalRequiredClones,
          transaction.LocalCapacityAvailable);
    std::abort();
  }

  const long long localBefore=PIC::ParticleBuffer::NAllPart;
  long long localSplitPlans=0,localRoulettePlans=0;
  long long localEmptyDaughterPlans=0;
  long long localChemistryConstrainedPlans=0;
  const long long localPlanCount=static_cast<long long>(
      transaction.LocalPlans.size());
  long long globalPlanCount=0;
  MPI_Allreduce(&localPlanCount,&globalPlanCount,1,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  if (globalPlanCount>0)
    Comet::Chemistry::PrepareForRegularCarrierChange();

  for (const cWeightAdaptationPlan &plan:transaction.LocalPlans) {
    if (plan.PhotochemicalCarrierConstraint)
      localChemistryConstrainedPlans++;
    if (plan.ParticleCount==0) {
      if (plan.EmptyPhotochemicalDaughter)
        localEmptyDaughterPlans++;
      SetEmptySpeciesWeight(plan);
    }
    else if (plan.WeightPower<0) {
      localSplitPlans++;
      for (int step=0;step < -plan.WeightPower;step++)
        SplitSpeciesByTwo(plan);
    }
    else {
      localRoulettePlans++;
      for (int step=0;step<plan.WeightPower;step++)
        RouletteSpeciesByTwo(plan);
    }
  }

  long long localStatistics[7]={localBefore,PIC::ParticleBuffer::NAllPart,
      static_cast<long long>(transaction.LocalPlans.size()),localSplitPlans,
      localRoulettePlans,localEmptyDaughterPlans,
      localChemistryConstrainedPlans};
  long long globalStatistics[7]={0,0,0,0,0,0,0};
  MPI_Allreduce(localStatistics,globalStatistics,7,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  transaction.PropagationStatistics.LocalChangedBlockSpecies=
      static_cast<long long>(transaction.LocalPlans.size());

  if (globalStatistics[2]>0) SynchronizeAdaptedBlockWeights();

  // A regular transaction never writes a locked owner value.  Compare its
  // exact IEEE representation after synchronization; a tolerance would hide
  // a one-bit source carrier drift and invalidate the bitwise source-lock
  // contract used by the Stage-2 pilot.
  for (const cLockedWeightSnapshot &snapshot:transaction.LockedWeights) {
    const double current=
        snapshot.Block->GetLocalParticleWeight(snapshot.Species);
    if (DoubleBits(current)!=snapshot.WeightBits) {
      std::fprintf(stderr,
          "component transaction changed locked species %d weight bits\n",
          snapshot.Species);
      std::abort();
    }
  }

  long long localMaximumProductsAfter=0;
  for (std::size_t block=0;block<transaction.Nodes.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=transaction.Nodes[block];
    if (node->Thread!=PIC::ThisThread) continue;
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      const double actual=GetCarrierRate(node->block,species);
      const double expected=transaction.TargetCarriers[offset];
      const double scale=std::max(std::fabs(actual),std::fabs(expected));
      if (std::fabs(actual-expected)>
          1024.0*std::numeric_limits<double>::epsilon()*scale) {
        std::fprintf(stderr,
            "committed component carrier differs from its plan for block %zu species %d\n",
            block,species);
        std::abort();
      }
    }
    const std::array<long int,PIC::nTotalSpecies> counts=
        CountAllSpeciesInBlock(node->block);
    localMaximumProductsAfter=std::max(localMaximumProductsAfter,
        MaximumReactionDaughterAllocation(node->block,counts));
  }
  long long globalMaximumProductsAfter=0;
  MPI_Allreduce(&localMaximumProductsAfter,&globalMaximumProductsAfter,1,
      MPI_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if (Comet::Chemistry::Enabled &&
      globalMaximumProductsAfter>
          Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction) {
    std::fprintf(stderr,
        "component transaction left %lld products per event above limit %ld\n",
        globalMaximumProductsAfter,
        Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction);
    std::abort();
  }

  LastFaceCarrierPropagationStatistics=
      transaction.PropagationStatistics;
  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: radial particle-weight adaptation iteration=%ld range=%ld:%ld particles=%lld->%lld changed-block-species=%lld split=%lld roulette=%lld empty-photochemical-daughter=%lld chemistry-constrained=%lld face-carrier-changes=%lld face-components=%lld max-face-hops=%lld max-face-clone-mean=%.6e->%.6e max-products-per-reaction=%lld->%lld limit=%ld component-planner=on required-clones=%lld\n",
        completedIterations,
        Comet::RadialSimulation::MinimumParticlesPerCell,
        Comet::RadialSimulation::MaximumParticlesPerCell,
        globalStatistics[0],globalStatistics[1],globalStatistics[2],
        globalStatistics[3],globalStatistics[4],globalStatistics[5],
        globalStatistics[6],
        transaction.PropagationStatistics.GlobalChangedBlockSpecies,
        transaction.PropagationStatistics.ActiveComponents,
        transaction.PropagationStatistics.MaximumPropagationHops,
        transaction.PropagationStatistics.MaximumCloneMeanBefore,
        transaction.PropagationStatistics.MaximumCloneMeanAfter,
        transaction.GlobalMaximumProductsBefore,
        globalMaximumProductsAfter,
        Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction,
        transaction.LocalRequiredClones);
    std::fflush(PIC::DiagnospticMessageStream);
  }
}

constexpr double ValidationCorrectionValues[]={0.0625,0.25,1.0,4.0,16.0};
constexpr std::size_t ValidationCorrectionValueCount=
    sizeof(ValidationCorrectionValues)/sizeof(ValidationCorrectionValues[0]);

void SeedValidationParticles(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int species,
    long int particleCount);

struct cValidationCorrectionAudit {
  long long ParticleCount=0;
  long long InvalidCount=0;
  double MaximumRelativeError=0.0;
  double CorrectionSum=0.0;
  double CorrectionSquareSum=0.0;
};

cValidationCorrectionAudit AuditValidationCorrections(int species) {
  cValidationCorrectionAudit local,global;
  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> nodes;
  CollectCalculationLeafNodes(PIC::Mesh::mesh->rootTree,nodes);
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node:nodes) {
    if (node->Thread!=PIC::ThisThread || node->block==NULL) continue;
    for (long int cell=0;cell<BlockCellCount;cell++)
      for (long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;
           ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        if (PIC::ParticleBuffer::GetI(data)!=species) continue;
        local.ParticleCount++;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        const double correction=
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
        double minimumError=std::numeric_limits<double>::infinity();
        for (double expected:ValidationCorrectionValues)
          minimumError=std::min(minimumError,
              std::fabs(correction-expected)/expected);
        if (!(correction>0.0) || !std::isfinite(correction) ||
            minimumError>256.0*std::numeric_limits<double>::epsilon())
          local.InvalidCount++;
        else {
          local.CorrectionSum+=correction;
          local.CorrectionSquareSum+=correction*correction;
        }
        if (std::isfinite(minimumError))
          local.MaximumRelativeError=std::max(
              local.MaximumRelativeError,minimumError);
#endif
      }
  }
  MPI_Allreduce(&local.ParticleCount,&global.ParticleCount,1,MPI_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.InvalidCount,&global.InvalidCount,1,MPI_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.MaximumRelativeError,&global.MaximumRelativeError,1,
      MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.CorrectionSum,&global.CorrectionSum,1,MPI_DOUBLE,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.CorrectionSquareSum,&global.CorrectionSquareSum,1,
      MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  return global;
}

bool NearlyEqualValidation(double left,double right,double tolerance=2.0e-12) {
  const double scale=std::max({std::fabs(left),std::fabs(right),
      std::numeric_limits<double>::min()});
  return std::fabs(left-right)<=tolerance*scale;
}

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
// Exercise the real core sampler with correction correlated to velocity and
// rotational energy.  Reading the collecting buffer directly avoids a
// sampling-buffer switch that would perturb the focused validator's global
// state, while still invoking exactly the production ProcessCell() path.
bool ValidateIndividualCorrectionSampling(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int &localChecks) {
  if (node==NULL || node->block==NULL) return true;
  int species=-1;
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  for (int candidate=0;candidate<PIC::nTotalSpecies;candidate++)
    if (PIC::IDF::nTotalRotationalModes[candidate]>0) {
      species=candidate;
      break;
    }
#endif
  if (species<0) return false;

  PIC::ParticleBuffer::DeleteAllParticles();
  SeedValidationParticles(node,species,
      static_cast<long int>(ValidationCorrectionValueCount));
  PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(
      _getCenterNodeLocalNumber(0,0,0));
  if (cell==NULL) return false;
  PIC::Mesh::flushCollectingSamplingBuffer(cell);

  const double baseWeight=node->block->GetLocalParticleWeight(species);
  double expectedWeight=0.0,expectedDensity=0.0;
  double expectedVelocity[3]={0.0,0.0,0.0};
  double expectedRotationalEnergy=0.0;
  int particleIndex=0;
  for (long int ptr=node->block->FirstCellParticleTable[0];ptr!=-1;
       ptr=PIC::ParticleBuffer::GetNext(ptr),particleIndex++) {
    PIC::ParticleBuffer::byte *data=
        PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    const double correction=
        PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
    const double representedWeight=baseWeight*correction;
    double velocity[3]={100.0+37.0*particleIndex,
        -21.0+5.0*particleIndex,9.0-3.0*particleIndex};
    PIC::ParticleBuffer::SetV(velocity,data);
    const double rotationalEnergy=
        ValidationCorrectionValues[particleIndex]*1.0e-21;
    PIC::IDF::LB::SetRotE(rotationalEnergy,data);
    expectedWeight+=representedWeight;
    expectedDensity+=representedWeight/cell->Measure;
    for (int dimension=0;dimension<3;dimension++)
      expectedVelocity[dimension]+=representedWeight*velocity[dimension];
    expectedRotationalEnergy+=representedWeight*rotationalEnergy;
  }

  std::vector<int> sampledParticleCounts(PIC::nTotalSpecies,0);
  int *sampledParticleCountTable[1]={sampledParticleCounts.data()};
  PIC::Sampling::ProcessCell(0,0,0,sampledParticleCountTable,node,0);

  const double *sampledWeight=
      cell->GetSampleDatum_ptr(&PIC::Mesh::DatumParticleWeight,species);
  const double *sampledDensity=
      cell->GetSampleDatum_ptr(&PIC::Mesh::DatumNumberDensity,species);
  const double *sampledVelocity=
      cell->GetSampleDatum_ptr(&PIC::Mesh::DatumParticleVelocity,species);
  if (sampledWeight==NULL || sampledDensity==NULL || sampledVelocity==NULL)
    return false;
  char *samplingData=cell->GetAssociatedDataBufferPointer()+
      PIC::Mesh::collectingCellSampleDataPointerOffset;
  const double sampledInternalWeight=*(species+reinterpret_cast<double*>(
      samplingData+PIC::IDF::_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_));
  const double sampledRotational=*(species+reinterpret_cast<double*>(
      samplingData+PIC::IDF::_ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_));

  bool pass=sampledParticleCounts[species]==
      static_cast<int>(ValidationCorrectionValueCount) &&
      NearlyEqualValidation(sampledWeight[0],expectedWeight) &&
      NearlyEqualValidation(sampledDensity[0],expectedDensity) &&
      NearlyEqualValidation(sampledInternalWeight,expectedWeight) &&
      NearlyEqualValidation(sampledRotational,expectedRotationalEnergy);
  for (int dimension=0;dimension<3;dimension++)
    pass=pass && NearlyEqualValidation(sampledVelocity[dimension],
        expectedVelocity[dimension]);
  if (pass) localChecks++;
  PIC::Mesh::flushCollectingSamplingBuffer(cell);
  return pass;
}

struct cCollisionSpeciesMoment {
  double Weight=0.0;
  double Velocity[3]={0.0,0.0,0.0};
};

cCollisionSpeciesMoment CollectCollisionSpeciesMoment(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int species) {
  cCollisionSpeciesMoment result;
  const double baseWeight=node->block->GetLocalParticleWeight(species);
  for (long int ptr=node->block->FirstCellParticleTable[0];ptr!=-1;
       ptr=PIC::ParticleBuffer::GetNext(ptr)) {
    PIC::ParticleBuffer::byte *data=
        PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    if (PIC::ParticleBuffer::GetI(data)!=species) continue;
    const double weight=baseWeight*
        PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
    double velocity[3];
    PIC::ParticleBuffer::GetV(velocity,data);
    result.Weight+=weight;
    for (int dimension=0;dimension<3;dimension++)
      result.Velocity[dimension]+=weight*velocity[dimension];
  }
  if (result.Weight>0.0)
    for (int dimension=0;dimension<3;dimension++)
      result.Velocity[dimension]/=result.Weight;
  return result;
}

int ValidationNTCSpecies0=-1;
int ValidationNTCSpecies1=-1;

bool SelectValidationNTCCrossSpeciesPair(
    PIC::Mesh::cDataCenterNode*,int species0,int species1) {
  // The Stage-4 collision test seeds two deliberately separated streams.  A
  // pair filter is required because updates from same-species collisions add
  // work to the core diagnostic but cannot relax the separation that the
  // test is intended to measure.  This hook is the normal general-NTC pair
  // selector; production state is restored immediately after the focused
  // kernel loop.
  return (species0==ValidationNTCSpecies0 &&
              species1==ValidationNTCSpecies1) ||
      (species0==ValidationNTCSpecies1 &&
              species1==ValidationNTCSpecies0);
}

// Run the selected production general-NTC cell kernel.  This test is aimed at
// compatibility/invariants rather than a new collision-physics calibration:
// it requires nonzero relaxation work, correction retention, and instrumented
// finite update probabilities in [0,1].
bool ValidateIndividualCorrectionNTC(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int &localChecks,
    unsigned long long &localProbabilityChecks,double &localMinimumProbability,
    double &localMaximumProbability) {
  if (node==NULL || node->block==NULL) return true;
  int species0=-1,species1=-1;
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  for (int species=0;species<PIC::nTotalSpecies;species++)
    if (PIC::MolecularData::GetSpecieType(species)==
            _PIC_SPECIE_TYPE__GAS_ &&
        (PIC::IDF::nTotalRotationalModes[species]>0 ||
         PIC::IDF::nTotalVibtationalModes[species]>0)) {
      species0=species;
      break;
    }
#endif
  for (int species=0;species<PIC::nTotalSpecies;species++)
    if (PIC::MolecularData::GetSpecieType(species)==
        _PIC_SPECIE_TYPE__GAS_ && species!=species0) {
      if (species0<0) species0=species;
      else {species1=species;break;}
    }
  if (species0<0 || species1<0) return false;

  PIC::ParticleBuffer::DeleteAllParticles();
  constexpr long int ParticlesPerSpecies=96;
  SeedValidationParticles(node,species0,ParticlesPerSpecies);
  SeedValidationParticles(node,species1,ParticlesPerSpecies);
  int index0=0,index1=0;
  for (long int ptr=node->block->FirstCellParticleTable[0];ptr!=-1;
       ptr=PIC::ParticleBuffer::GetNext(ptr)) {
    PIC::ParticleBuffer::byte *data=
        PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    const int species=PIC::ParticleBuffer::GetI(data);
    const double correction=
        PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
    double velocity[3];
    if (species==species0) {
      velocity[0]=700.0+8.0*index0+3.0*correction;
      velocity[1]=30.0; velocity[2]=-10.0; index0++;
    }
    else {
      velocity[0]=-500.0-6.0*index1-2.0*correction;
      velocity[1]=-20.0; velocity[2]=15.0; index1++;
    }
    PIC::ParticleBuffer::SetV(velocity,data);
  }
  const cCollisionSpeciesMoment before0=
      CollectCollisionSpeciesMoment(node,species0);
  const cCollisionSpeciesMoment before1=
      CollectCollisionSpeciesMoment(node,species1);
  const double initialSeparation=
      std::fabs(before0.Velocity[0]-before1.Velocity[0]);

  PIC::MolecularCollisions::ParticleCollisionModel::
      ResetNTCValidationDiagnostic(true);
  // A fixed number of kernel calls is not a meaningful amount of collision
  // work when c spans 1/16...16: the exact unequal-carrier update probability
  // can be small, and 32 calls produced only a few updates in the first
  // Stage-4 smoke test.  Stop on completed update work instead, retaining a
  // hard cycle ceiling so a broken/zero-probability operator fails promptly.
  // Two species with 96 particles each are seeded.  Requiring 256 updates
  // therefore exercises more than one cross-species update per model
  // particle on average; the independent drift-relaxation predicate below
  // still prevents an update counter alone from declaring success.
  constexpr unsigned long long TargetParticleUpdates=256;
  constexpr int MaximumCollisionCycles=65536;
  const PIC::MolecularCollisions::ParticleCollisionModel::
      fDoSimulateCellCollisions savedCollisionSelector=
      PIC::MolecularCollisions::ParticleCollisionModel::
          DoSimulateCellCollisions;
  ValidationNTCSpecies0=species0;
  ValidationNTCSpecies1=species1;
  PIC::MolecularCollisions::ParticleCollisionModel::
      DoSimulateCellCollisions=SelectValidationNTCCrossSpeciesPair;
  int collisionCycles=0;
  for (;collisionCycles<MaximumCollisionCycles;collisionCycles++) {
    PIC::MolecularCollisions::ParticleCollisionModel::
        ModelCellCollisions_ntc(0,0,0,node);
    if ((collisionCycles+1)%32==0 &&
        PIC::MolecularCollisions::ParticleCollisionModel::
            GetNTCValidationDiagnostic().ParticleUpdateCount>=
                TargetParticleUpdates) {
      collisionCycles++;
      break;
    }
  }
  PIC::MolecularCollisions::ParticleCollisionModel::
      DoSimulateCellCollisions=savedCollisionSelector;
  ValidationNTCSpecies0=-1;
  ValidationNTCSpecies1=-1;
  const PIC::MolecularCollisions::ParticleCollisionModel::
      cNTCValidationDiagnostic diagnostic=
      PIC::MolecularCollisions::ParticleCollisionModel::
          GetNTCValidationDiagnostic();
  PIC::MolecularCollisions::ParticleCollisionModel::
      ResetNTCValidationDiagnostic(false);

  const cCollisionSpeciesMoment after0=
      CollectCollisionSpeciesMoment(node,species0);
  const cCollisionSpeciesMoment after1=
      CollectCollisionSpeciesMoment(node,species1);
  const double finalSeparation=
      std::fabs(after0.Velocity[0]-after1.Velocity[0]);
  const cValidationCorrectionAudit audit0=
      AuditValidationCorrections(species0);
  const cValidationCorrectionAudit audit1=
      AuditValidationCorrections(species1);
  const bool pass=diagnostic.TrialCount>0 &&
      diagnostic.AcceptedCollisionCount>0 &&
      diagnostic.ParticleUpdateCount>=TargetParticleUpdates &&
      diagnostic.UpdateProbabilityCount>0 &&
      diagnostic.InvalidUpdateProbabilityCount==0 &&
      std::isfinite(diagnostic.MinimumUpdateProbability) &&
      std::isfinite(diagnostic.MaximumUpdateProbability) &&
      diagnostic.MinimumUpdateProbability>=0.0 &&
      diagnostic.MaximumUpdateProbability<=1.0 &&
      audit0.InvalidCount==0 && audit1.InvalidCount==0 &&
      audit0.ParticleCount==ParticlesPerSpecies*PIC::nTotalThreads &&
      audit1.ParticleCount==ParticlesPerSpecies*PIC::nTotalThreads &&
      finalSeparation<initialSeparation;
  if (pass) localChecks++;
  localProbabilityChecks=diagnostic.UpdateProbabilityCount;
  localMinimumProbability=diagnostic.MinimumUpdateProbability;
  localMaximumProbability=diagnostic.MaximumUpdateProbability;
  if (PIC::ThisThread==0) {
    // Keep a compact, machine-readable failure record in every archive.  In
    // particular, a valid probability range alone does not prove that the
    // unequal-collision operator performed work or relaxed the two streams.
    // These fields distinguish a selection/update failure from correction
    // corruption or an over-strict relaxation assertion without requiring a
    // debugger rerun of the MPI validation case.
    std::fprintf(PIC::DiagnospticMessageStream,
        "%s individual-weight NTC validation cycles=%d target-updates=%llu trials=%llu accepted=%llu updates=%llu probability-checks=%llu invalid-probabilities=%llu probability-min=%.17e probability-max=%.17e species0-count=%lld species1-count=%lld species0-invalid-c=%lld species1-invalid-c=%lld initial-separation=%.17e final-separation=%.17e\n",
        pass ? "PASS" : "FAIL",collisionCycles,TargetParticleUpdates,
        diagnostic.TrialCount,
        diagnostic.AcceptedCollisionCount,diagnostic.ParticleUpdateCount,
        diagnostic.UpdateProbabilityCount,
        diagnostic.InvalidUpdateProbabilityCount,
        diagnostic.MinimumUpdateProbability,
        diagnostic.MaximumUpdateProbability,audit0.ParticleCount,
        audit1.ParticleCount,audit0.InvalidCount,audit1.InvalidCount,
        initialSeparation,finalSeparation);
    std::fflush(PIC::DiagnospticMessageStream);
  }
  return pass;
}

bool ValidateIndividualCorrectionBoundaryAndRestart(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int species,
    int &localBoundaryChecks,int &localRestartChecks) {
  if (node==NULL || node->block==NULL) return true;
  PIC::ParticleBuffer::DeleteAllParticles();
  SeedValidationParticles(node,species,1);
  long int ptr=node->block->FirstCellParticleTable[0];
  if (ptr<0) return false;
  PIC::ParticleBuffer::byte *data=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(4.0,data);
  double position[3],velocity[3];
  PIC::ParticleBuffer::GetX(position,data);
  PIC::ParticleBuffer::GetV(velocity,data);
  const int innerResult=Comet::RadialSimulation::ProcessDomainBoundary(ptr,
      position,velocity,0,node);
  const bool innerPass=innerResult==_PARTICLE_REJECTED_ON_THE_FACE_ &&
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data)==4.0;
  const int outerResult=Comet::RadialSimulation::ProcessDomainBoundary(ptr,
      position,velocity,1,node);
  const bool outerPass=outerResult==_PARTICLE_DELETED_ON_THE_FACE_ &&
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data)==4.0;
  if (innerPass && outerPass) localBoundaryChecks++;

  const unsigned long checksumBefore=PIC::Restart::GetParticleDataCheckSum();
  const double correctionBefore=
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
  char restartName[256];
  std::snprintf(restartName,sizeof(restartName),
      "/tmp/amps-hartley2-stage4.ranks=%d.restart",PIC::nTotalThreads);
  PIC::Restart::SaveParticleData(restartName);
  PIC::ParticleBuffer::DeleteAllParticles();
  PIC::Restart::ReadParticleData(restartName);
  const unsigned long checksumAfter=PIC::Restart::GetParticleDataCheckSum();
  const cValidationCorrectionAudit restartAudit=
      AuditValidationCorrections(species);
  const bool restartPass=checksumAfter==checksumBefore &&
      restartAudit.ParticleCount==PIC::nTotalThreads &&
      restartAudit.InvalidCount==0 && correctionBefore==4.0;
  if (restartPass) localRestartChecks++;
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) std::remove(restartName);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  return innerPass && outerPass && restartPass;
}

double Stage4MigrationLoadMeasure(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  // A stable, deliberately irregular AMR-ID hash changes multi-rank
  // ownership boundaries without depending on their old owners.  It is used
  // only by the focused validator, after all other initialized-mesh tests.
  std::uint64_t key=1469598103934665603ULL^
      static_cast<std::uint64_t>(node->AMRnodeID.ResolutionLevel);
  for (unsigned int digit=0;digit<node->AMRnodeID.Length();digit++) {
    key^=static_cast<std::uint64_t>(node->AMRnodeID.id[digit])+1ULL;
    key*=1099511628211ULL;
  }
  return 1.0+static_cast<double>(key%97ULL);
}

bool ValidateIndividualCorrectionMigration(int species,int &localChecks) {
  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> nodes;
  CollectCalculationLeafNodes(PIC::Mesh::mesh->rootTree,nodes);
  std::vector<int> oldOwner(nodes.size(),-1);
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *seedNode=NULL;
  for (std::size_t block=0;block<nodes.size();block++) {
    oldOwner[block]=nodes[block]->Thread;
    if (seedNode==NULL && nodes[block]->Thread==PIC::ThisThread &&
        nodes[block]->block!=NULL) seedNode=nodes[block];
  }
  if (seedNode==NULL) return false;

  PIC::ParticleBuffer::DeleteAllParticles();
  SeedValidationParticles(seedNode,species,1);
  PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(
      seedNode->block->FirstCellParticleTable[0]);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
      ValidationCorrectionValues[PIC::ThisThread%
          ValidationCorrectionValueCount],data);
  const cValidationCorrectionAudit before=
      AuditValidationCorrections(species);

  PIC::Mesh::mesh->SetParallelLoadMeasure(Stage4MigrationLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists(
      _PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);
  PIC::DomainBlockDecomposition::UpdateBlockTable();

  long long localMigratedBlocks=0;
  for (std::size_t block=0;block<nodes.size();block++)
    if (nodes[block]->Thread==PIC::ThisThread &&
        oldOwner[block]!=nodes[block]->Thread) localMigratedBlocks++;
  long long migratedBlocks=0;
  MPI_Allreduce(&localMigratedBlocks,&migratedBlocks,1,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  const cValidationCorrectionAudit after=
      AuditValidationCorrections(species);
  const bool pass=after.InvalidCount==0 &&
      after.ParticleCount==before.ParticleCount &&
      after.CorrectionSum==before.CorrectionSum &&
      after.CorrectionSquareSum==before.CorrectionSquareSum &&
      (PIC::nTotalThreads==1 || migratedBlocks>0);
  if (pass) localChecks++;
  return pass;
}
#endif

void SeedValidationParticles(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int species,
    long int particleCount) {
  if (node==NULL || node->block==NULL) return;
  long int &first=node->block->FirstCellParticleTable[0];
  const double cellWidth=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
  const double position[3]={node->xmin[0]+0.5*cellWidth,0.0,0.0};
  // Use a non-zero, identical state for every synthetic particle.  Exact
  // splitting must duplicate every stored degree of freedom, while roulette
  // must retain an unbiased subset.  With a common state the represented
  // number, mass, momentum, kinetic energy, and internal energy all have the
  // same analytically known split/roulette scaling and can be checked without
  // confusing adaptation error with thermal sampling noise.
  const double velocity[3]={417.0,-31.0,12.0};
  // A wide deterministic correction cycle makes the adaptation validator
  // sensitive to CloneParticle() accidentally resetting c and to roulette
  // moments that use W instead of the represented carrier W*c.  The cycle is
  // compile-time inert in the default correction-off production build.
  for (long int particle=0;particle<particleCount;particle++) {
    const long int ptr=PIC::ParticleBuffer::GetNewParticle(first);
    PIC::ParticleBuffer::byte *data=
        PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    PIC::ParticleBuffer::SetX(const_cast<double*>(position),data);
    PIC::ParticleBuffer::SetV(const_cast<double*>(velocity),data);
    PIC::ParticleBuffer::SetI(species,data);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
        ValidationCorrectionValues[(particle+species)%
            ValidationCorrectionValueCount],data);
#endif
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
    PIC::IDF::InitRotTemp(
        Comet::RadialSimulation::InnerBoundaryTemperature,data);
    PIC::IDF::InitVibTemp(
        Comet::RadialSimulation::InnerBoundaryTemperature,data);
    if (PIC::IDF::nTotalRotationalModes[species]>0)
      PIC::IDF::LB::SetRotE(1.0e-21*(species+1),data);
    for (int mode=0;mode<PIC::IDF::nTotalVibtationalModes[species];mode++)
      PIC::IDF::LB::SetVibE(2.0e-22*(mode+1),mode,data);
#endif
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
    PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
#endif
  }
}

// Integrated represented moments used by the Stage-2 transaction tests.
// WeightSquared is not a physical moment; it is the exact variance coefficient
// for independent 1/2 roulette with unequal represented weights W*c.
struct cValidationMoments {
  double Number=0.0;
  double Mass=0.0;
  double Momentum[3]={0.0,0.0,0.0};
  double KineticEnergy=0.0;
  double InternalEnergy=0.0;
  double WeightSquared=0.0;
};

cValidationMoments CollectValidationMoments(int species) {
  cValidationMoments local,global;
  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> nodes;
  CollectCalculationLeafNodes(PIC::Mesh::mesh->rootTree,nodes);
  const double mass=PIC::MolecularData::GetMass(species);
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node:nodes) {
    if (node->Thread!=PIC::ThisThread || node->block==NULL) continue;
    const double baseWeight=node->block->GetLocalParticleWeight(species);
    for (long int cell=0;cell<BlockCellCount;cell++)
      for (long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;
           ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        if (PIC::ParticleBuffer::GetI(data)!=species) continue;
        double weight=baseWeight;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        const double correction=
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
        if (!(correction>0.0) || !std::isfinite(correction)) {
          std::fprintf(stderr,
              "invalid individual correction in adaptation validation\n");
          std::abort();
        }
        weight*=correction;
#endif
        double velocity[3];
        PIC::ParticleBuffer::GetV(velocity,data);
        local.Number+=weight;
        local.Mass+=weight*mass;
        double speedSquared=0.0;
        for (int dimension=0;dimension<3;dimension++) {
          local.Momentum[dimension]+=weight*mass*velocity[dimension];
          speedSquared+=velocity[dimension]*velocity[dimension];
        }
        local.KineticEnergy+=0.5*weight*mass*speedSquared;
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        double internal=0.0;
        if (PIC::IDF::nTotalRotationalModes[species]>0)
          internal+=PIC::IDF::GetRotE(data);
        for (int mode=0;mode<PIC::IDF::nTotalVibtationalModes[species];mode++)
          internal+=PIC::IDF::GetVibE(mode,data);
        local.InternalEnergy+=weight*internal;
#endif
        local.WeightSquared+=weight*weight;
      }
  }

  // Keep this as a contiguous explicit pack instead of relying on struct
  // padding in MPI.  Every entry is a sum and therefore decomposition
  // independent to ordinary floating-point reduction tolerance.
  double localValues[8]={local.Number,local.Mass,local.Momentum[0],
      local.Momentum[1],local.Momentum[2],local.KineticEnergy,
      local.InternalEnergy,local.WeightSquared};
  double globalValues[8]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  MPI_Allreduce(localValues,globalValues,8,MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  global.Number=globalValues[0];
  global.Mass=globalValues[1];
  for (int dimension=0;dimension<3;dimension++)
    global.Momentum[dimension]=globalValues[2+dimension];
  global.KineticEnergy=globalValues[5];
  global.InternalEnergy=globalValues[6];
  global.WeightSquared=globalValues[7];
  return global;
}

bool ValidationMomentsEqual(const cValidationMoments &left,
    const cValidationMoments &right,double relativeTolerance) {
  const double leftValues[7]={left.Number,left.Mass,left.Momentum[0],
      left.Momentum[1],left.Momentum[2],left.KineticEnergy,
      left.InternalEnergy};
  const double rightValues[7]={right.Number,right.Mass,right.Momentum[0],
      right.Momentum[1],right.Momentum[2],right.KineticEnergy,
      right.InternalEnergy};
  for (int value=0;value<7;value++) {
    const double scale=std::max(std::fabs(leftValues[value]),
        std::fabs(rightValues[value]));
    if (scale==0.0) continue;
    if (std::fabs(leftValues[value]-rightValues[value])>
        relativeTolerance*scale) return false;
  }
  return true;
}

bool OtherValidationCellsAreEmpty(PIC::Mesh::cDataBlockAMR *block) {
  for (long int cell=1;cell<BlockCellCount;cell++)
    if (block->FirstCellParticleTable[cell]!=-1) return false;
  return true;
}

cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *SelectValidationOwnerBlock() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *fallback=NULL;
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL) continue;
    if (fallback==NULL) fallback=node;

    // Prefer an owner block touching another rank.  Its non-owner boundary
    // copy lets the MPI validator measure, rather than merely execute, the
    // block-weight synchronization performed after adaptation.
    for (int face=0;face<2;face++) {
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neighbor=
          node->GetNeibFace(face,0,0,PIC::Mesh::mesh);
      if (neighbor!=NULL && neighbor->Thread!=PIC::ThisThread) return node;
    }
  }
  return fallback;
}

bool ValidateSynchronizedWeightCopies(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *ownerNode,int species,
    double expectedOwnerWeight,int &localGhostChecks) {
  // All ranks publish the center and expected weight of their selected owner
  // block.  A rank that has a boundary-layer allocation for another owner's
  // block must see exactly the synchronized value in that non-owner copy.
  const double localDescriptor[3]={
      ownerNode!=NULL ? 1.0 : 0.0,
      ownerNode!=NULL ? 0.5*(ownerNode->xmin[0]+ownerNode->xmax[0]) : 0.0,
      ownerNode!=NULL ? expectedOwnerWeight : 0.0};
  std::vector<double> descriptors(3*PIC::nTotalThreads,0.0);
  MPI_Allgather(localDescriptor,3,MPI_DOUBLE,descriptors.data(),3,MPI_DOUBLE,
      MPI_GLOBAL_COMMUNICATOR);

  bool pass=true;
  for (int owner=0;owner<PIC::nTotalThreads;owner++) {
    if (descriptors[3*owner]==0.0 || owner==PIC::ThisThread) continue;
    double position[3]={descriptors[3*owner+1],0.0,0.0};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *copy=
        PIC::Mesh::mesh->findTreeNode(position);
    if (copy==NULL || copy->block==NULL || copy->Thread==PIC::ThisThread)
      continue;

    localGhostChecks++;
    const double expected=descriptors[3*owner+2];
    const double actual=copy->block->GetLocalParticleWeight(species);
    const double scale=std::max(std::fabs(expected),1.0);
    if (std::fabs(actual-expected)>1.0e-13*scale) pass=false;
  }
  return pass;
}
#endif

long int CountProfileCells(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return 0;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_)
    return _BLOCK_CELLS_X_;

  long int result=0;
  for (int child=0;child<(1<<DIM);child++)
    result+=CountProfileCells(node->downNode[child]);
  return result;
}

void FillLocalIntegratedDensityProfile(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,long int &cellIndex,
    std::vector<double> &profile) {
  if (node==NULL) return;

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    // Every MPI rank traverses the complete AMR tree in the same child order,
    // so cellIndex is a stable global radial-cell index even when ownership is
    // redistributed at a sampling break. Only the current owner contributes;
    // MPI_Allreduce below assembles the complete profile on every rank.
    if (node->Thread==PIC::ThisThread && node->block!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) {
        PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i,0,0));
        if (cell==NULL || !(cell->Measure>0.0) ||
            !std::isfinite(cell->Measure)) {
          std::fprintf(stderr,
              "invalid cell measure while evaluating radial convergence\n");
          std::abort();
        }

        for (int species=0;species<PIC::nTotalSpecies;species++) {
          const double density=cell->GetNumberDensity(species);
          if (!(density>=0.0) || !std::isfinite(density)) {
            std::fprintf(stderr,
                "invalid number density for species %d while evaluating radial convergence\n",
                species);
            std::abort();
          }

          // Storing n*V makes the subsequent discrete L1 norm exactly the
          // shell-volume-weighted norm. Measure is the core's 4*pi/3 times
          // (r_outer^3-r_inner^3) value in DIM1 spherical symmetry.
          profile[(cellIndex+i)*PIC::nTotalSpecies+species]=
              density*cell->Measure;
        }
      }
    }

    cellIndex+=_BLOCK_CELLS_X_;
    return;
  }

  for (int child=0;child<(1<<DIM);child++)
    FillLocalIntegratedDensityProfile(node->downNode[child],cellIndex,profile);
}

std::vector<double> CollectIntegratedDensityProfile() {
  const long int cellCount=CountProfileCells(PIC::Mesh::mesh->rootTree);
  if (cellCount<=0) {
    std::fprintf(stderr,"radial convergence found no mesh cells\n");
    std::abort();
  }

  const size_t valueCount=
      static_cast<size_t>(cellCount)*PIC::nTotalSpecies;
  if (valueCount>static_cast<size_t>(INT_MAX)) {
    std::fprintf(stderr,
        "radial convergence profile exceeds the MPI count limit\n");
    std::abort();
  }

  std::vector<double> localProfile(valueCount,0.0);
  std::vector<double> globalProfile(valueCount,0.0);
  long int visitedCells=0;
  FillLocalIntegratedDensityProfile(PIC::Mesh::mesh->rootTree,visitedCells,
      localProfile);
  if (visitedCells!=cellCount) {
    std::fprintf(stderr,"inconsistent radial convergence cell traversal\n");
    std::abort();
  }

  MPI_Allreduce(localProfile.data(),globalProfile.data(),
      static_cast<int>(valueCount),MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  return globalProfile;
}

void FillLocalStage3CellSamples(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,long int &cellIndex,
    std::vector<double> &samples,std::vector<double> &geometry) {
  if (node==NULL) return;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->Thread==PIC::ThisThread && node->block!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) {
        PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i,0,0));
        if (cell==NULL) {
          std::fprintf(stderr,"missing center node in Stage-3 cell output\n");
          std::abort();
        }
        const double width=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
        const std::size_t geometryOffset=4*static_cast<std::size_t>(cellIndex+i);
        geometry[geometryOffset]=cell->GetX()[0];
        geometry[geometryOffset+1]=node->xmin[0]+i*width;
        geometry[geometryOffset+2]=node->xmin[0]+(i+1)*width;
        geometry[geometryOffset+3]=cell->Measure;
        for (int species=0;species<PIC::nTotalSpecies;species++) {
          const double density=cell->GetNumberDensity(species);
          const double occupancy=cell->GetDatumAverage(
              PIC::Mesh::DatumParticleNumber,species);
          if (!(density>=0.0) || !std::isfinite(density) ||
              !(occupancy>=0.0) || !std::isfinite(occupancy)) {
            std::fprintf(stderr,
                "invalid Stage-3 cell sample for species %d\n",species);
            std::abort();
          }
          const std::size_t offset=2*(
              static_cast<std::size_t>(cellIndex+i)*PIC::nTotalSpecies+
              species);
          samples[offset]=density;
          samples[offset+1]=occupancy;
        }
      }
    }
    cellIndex+=_BLOCK_CELLS_X_;
    return;
  }
  for (int child=0;child<(1<<DIM);child++)
    FillLocalStage3CellSamples(node->downNode[child],cellIndex,samples,
        geometry);
}

void WriteStage3CellSamplesInternal(long int outputIndex,long int sampleLength) {
  if (!Comet::Diagnostics::Enabled) return;
  const long int cellCount=CountProfileCells(PIC::Mesh::mesh->rootTree);
  const std::size_t valueCount=2*static_cast<std::size_t>(cellCount)*
      PIC::nTotalSpecies;
  const std::size_t geometryValueCount=4*static_cast<std::size_t>(cellCount);
  if (cellCount<=0 || valueCount>static_cast<std::size_t>(INT_MAX) ||
      geometryValueCount>static_cast<std::size_t>(INT_MAX)) {
    std::fprintf(stderr,"Stage-3 cell output exceeds the MPI count limit\n");
    std::abort();
  }
  std::vector<double> local(valueCount,0.0),global(valueCount,0.0);
  std::vector<double> localGeometry(geometryValueCount,0.0);
  std::vector<double> globalGeometry(localGeometry.size(),0.0);
  long int visited=0;
  FillLocalStage3CellSamples(PIC::Mesh::mesh->rootTree,visited,local,
      localGeometry);
  if (visited!=cellCount) {
    std::fprintf(stderr,"inconsistent Stage-3 cell-output traversal\n");
    std::abort();
  }
  MPI_Reduce(local.data(),global.data(),static_cast<int>(valueCount),
      MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localGeometry.data(),globalGeometry.data(),
      static_cast<int>(localGeometry.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread!=0) return;

  // Tecplot tables are node-interpolated and duplicate AMR interfaces.  They
  // cannot be joined exactly to shell measures for an exposure decision.
  // Emit a separate center-cell table, preserving the original plot files.
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    char symbol[_MAX_STRING_LENGTH_PIC_];
    PIC::MolecularData::GetChemSymbol(symbol,species);
    // Build the diagnostic name dynamically.  Chemical symbols normally are
    // short, but avoiding a fixed path buffer makes the output contract safe
    // for user-defined species names as the generic comet application grows.
    const std::string path=std::string("PT/diagnostics/stage3_cells.")+symbol+
        ".out="+std::to_string(outputIndex)+".csv";
    FILE *file=std::fopen(path.c_str(),"w");
    if (file==NULL) {
      std::fprintf(stderr,"cannot open Stage-3 cell diagnostic %s\n",
          path.c_str());
      std::abort();
    }
    std::fprintf(file,
        "cell_index,block_index,local_cell_index,block_cell_count,radius_m,r_inner_m,r_outer_m,measure_m3,number_density_m-3,mean_model_particles,population_target_min_per_cell,population_target_max_per_cell,population_locked,population_class,population_policy_source,sampling_length,output_index,random_seed,random_seed_configured,residual_sampling,mpi_ranks\n");
    const char *residualMode=Comet::Chemistry::ResidualSamplingMode==
        Comet::Chemistry::ResidualSampling::RandomPhaseSystematic ?
        "SYSTEMATIC" : "BERNOULLI";
    for (long int cell=0;cell<cellCount;cell++) {
      const std::size_t offset=2*(
          static_cast<std::size_t>(cell)*PIC::nTotalSpecies+species);
      const std::size_t geometryOffset=4*static_cast<std::size_t>(cell);
      if (!(globalGeometry[geometryOffset+3]>0.0) ||
          !(globalGeometry[geometryOffset+1]<globalGeometry[geometryOffset]) ||
          !(globalGeometry[geometryOffset]<globalGeometry[geometryOffset+2])) {
        std::fprintf(stderr,"invalid reduced Stage-3 cell geometry\n");
        std::abort();
      }
      std::fprintf(file,
          "%ld,%ld,%ld,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%ld,%ld,%d,%d,%d,%ld,%ld,%llu,%d,%s,%d\n",
          cell,cell/_BLOCK_CELLS_X_,cell%_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,
          globalGeometry[geometryOffset],globalGeometry[geometryOffset+1],
          globalGeometry[geometryOffset+2],globalGeometry[geometryOffset+3],
          global[offset],global[offset+1],
          Comet::RadialSimulation::SpeciesMinimumParticlesPerCell[species],
          Comet::RadialSimulation::SpeciesMaximumParticlesPerCell[species],
          Comet::RadialSimulation::ParticleWeightLocked[species],
          Comet::RadialSimulation::PopulationClass[species],
          Comet::RadialSimulation::PopulationPolicySource[species],
          sampleLength,outputIndex,
          static_cast<unsigned long long>(
              Comet::Model::RandomSeedConfigured ?
              Comet::Model::RandomSeed : 0),
          Comet::Model::RandomSeedConfigured ? 1 : 0,
          residualMode,PIC::nTotalThreads);
    }
    std::fclose(file);
  }
}

// Accumulate instantaneous, individually weighted carriers in the same
// authoritative radial-cell ordering used by the Stage-3 density tables.
// q=W_block,species*c_particle is the physical carrier represented by one
// model particle.  The Kish carrier ESS in a cell is therefore
// (sum q)^2/sum(q^2); using corrections alone would be wrong after local
// block-weight adaptation changes W between radial blocks.
void FillLocalStage9CarrierSamples(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,long int &cellIndex,
    std::vector<double> &weightSum,std::vector<double> &weightSquareSum,
    std::vector<double> &minimumCorrection,
    std::vector<double> &maximumCorrection,
    std::vector<unsigned long long> &carrierCount) {
  if (node==NULL) return;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->Thread==PIC::ThisThread && node->block!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++)
        for (long int ptr=node->block->FirstCellParticleTable[i];ptr!=-1;
             ptr=PIC::ParticleBuffer::GetNext(ptr)) {
          PIC::ParticleBuffer::byte *data=
              PIC::ParticleBuffer::GetParticleDataPointer(ptr);
          const int species=PIC::ParticleBuffer::GetI(data);
          if (species<0 || species>=PIC::nTotalSpecies) {
            std::fprintf(stderr,
                "invalid species in Stage-9 carrier diagnostic\n");
            std::abort();
          }
          const double correction=
              PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
          const double carrier=node->block->GetLocalParticleWeight(species)*
              correction;
          if (!(correction>0.0) || !std::isfinite(correction) ||
              !(carrier>0.0) || !std::isfinite(carrier) ||
              !std::isfinite(carrier*carrier)) {
            std::fprintf(stderr,
                "invalid Stage-9 carrier for species %d\n",species);
            std::abort();
          }
          const std::size_t offset=
              static_cast<std::size_t>(cellIndex+i)*PIC::nTotalSpecies+
              species;
          weightSum[offset]+=carrier;
          weightSquareSum[offset]+=carrier*carrier;
          minimumCorrection[offset]=std::min(minimumCorrection[offset],
              correction);
          maximumCorrection[offset]=std::max(maximumCorrection[offset],
              correction);
          carrierCount[offset]++;
        }
    }
    cellIndex+=_BLOCK_CELLS_X_;
    return;
  }
  for (int child=0;child<(1<<DIM);child++)
    FillLocalStage9CarrierSamples(node->downNode[child],cellIndex,weightSum,
        weightSquareSum,minimumCorrection,maximumCorrection,carrierCount);
}

void WriteStage9CheckpointInternal(long int completedIteration,
    long int outputIndex,long int sampleLength) {
  if (!Comet::RadialSimulation::Stage9ScienceDiagnosticsEnabled) return;
  AuditStage9SafeBoundary("science-checkpoint");
  UpdateStage9BufferPeak();

  const long int cellCount=CountProfileCells(PIC::Mesh::mesh->rootTree);
  const std::size_t valueCount=
      static_cast<std::size_t>(cellCount)*PIC::nTotalSpecies;
  if (cellCount<=0 || valueCount>static_cast<std::size_t>(INT_MAX)) {
    std::fprintf(stderr,"Stage-9 carrier output exceeds the MPI count limit\n");
    std::abort();
  }
  std::vector<double> localWeight(valueCount,0.0),
      localWeight2(valueCount,0.0),globalWeight(valueCount,0.0),
      globalWeight2(valueCount,0.0);
  std::vector<double> localMinimum(valueCount,
      std::numeric_limits<double>::infinity());
  std::vector<double> localMaximum(valueCount,0.0),globalMinimum(valueCount),
      globalMaximum(valueCount);
  std::vector<unsigned long long> localCount(valueCount,0),
      globalCount(valueCount,0);
  long int visited=0;
  FillLocalStage9CarrierSamples(PIC::Mesh::mesh->rootTree,visited,
      localWeight,localWeight2,localMinimum,localMaximum,localCount);
  if (visited!=cellCount) {
    std::fprintf(stderr,"inconsistent Stage-9 carrier traversal\n");
    std::abort();
  }
  MPI_Reduce(localWeight.data(),globalWeight.data(),static_cast<int>(valueCount),
      MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localWeight2.data(),globalWeight2.data(),
      static_cast<int>(valueCount),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localMinimum.data(),globalMinimum.data(),
      static_cast<int>(valueCount),MPI_DOUBLE,MPI_MIN,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localMaximum.data(),globalMaximum.data(),
      static_cast<int>(valueCount),MPI_DOUBLE,MPI_MAX,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localCount.data(),globalCount.data(),
      static_cast<int>(valueCount),MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);

  const long int localPeak=Stage9WindowPeakLiveParticles;
  const long int localCurrent=PIC::ParticleBuffer::NAllPart;
  const long int localCapacity=PIC::ParticleBuffer::MaxNPart;
  long int globalPeak=0,globalCurrent=0,minimumCapacity=0,maximumCapacity=0;
  MPI_Reduce(&localPeak,&globalPeak,1,MPI_LONG,MPI_MAX,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&localCurrent,&globalCurrent,1,MPI_LONG,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&localCapacity,&minimumCapacity,1,MPI_LONG,MPI_MIN,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&localCapacity,&maximumCapacity,1,MPI_LONG,MPI_MAX,0,
      MPI_GLOBAL_COMMUNICATOR);
  double globalPeakUtilization=0.0;
  MPI_Reduce(&Stage9WindowPeakBufferUtilization,&globalPeakUtilization,1,
      MPI_DOUBLE,MPI_MAX,0,MPI_GLOBAL_COMMUNICATOR);
  const double localWallSeconds=std::max(0.0,
      MPI_Wtime()-Stage9WindowStartWallTime);
  double globalWallSeconds=0.0;
  MPI_Reduce(&localWallSeconds,&globalWallSeconds,1,MPI_DOUBLE,MPI_MAX,0,
      MPI_GLOBAL_COMMUNICATOR);
  const unsigned long long localAudits=Stage9SafeBoundaryAudits;
  unsigned long long minimumAudits=0,maximumAudits=0;
  MPI_Reduce(&localAudits,&minimumAudits,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&localAudits,&maximumAudits,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,0,
      MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    const char *mode=Comet::Chemistry::GetEstimatorModeName();
    const std::string base="PT/diagnostics/comet_stage9_checkpoint.out="+
        std::to_string(outputIndex);
    FILE *csv=std::fopen((base+".csv").c_str(),"w");
    FILE *json=std::fopen((base+".json").c_str(),"w");
    if (csv==NULL || json==NULL) {
      if (csv!=NULL) std::fclose(csv);
      if (json!=NULL) std::fclose(json);
      std::fprintf(stderr,"cannot open Stage-9 checkpoint output\n");
      std::abort();
    }
    std::fprintf(csv,
        "iteration,output_index,sampling_length,estimator_mode,random_seed,mpi_ranks,species,model_carriers,physical_inventory,sum_squared_weights,carrier_ess,occupied_cells,minimum_positive_cell_ess,correction_min,correction_max,window_wall_seconds,peak_rank_live_particles,total_live_particles,rank_capacity_min,rank_capacity_max,peak_buffer_utilization,safe_boundary_audits,pending_remainders,pending_reaction_events\n");
    std::fprintf(json,
        "{\n  \"schema\": \"amps.comet.stage9-checkpoint.v1\",\n  \"iteration\": %ld,\n  \"output_index\": %ld,\n  \"sampling_length\": %ld,\n  \"estimator_mode\": \"%s\",\n  \"random_seed\": %llu,\n  \"mpi_ranks\": %d,\n  \"window_wall_seconds\": %.17e,\n  \"peak_rank_live_particles\": %ld,\n  \"total_live_particles\": %ld,\n  \"rank_capacity_min\": %ld,\n  \"rank_capacity_max\": %ld,\n  \"peak_buffer_utilization\": %.17e,\n  \"safe_boundary_audits_min\": %llu,\n  \"safe_boundary_audits_max\": %llu,\n  \"pending_remainders\": 0,\n  \"pending_reaction_events\": 0,\n  \"species\": [\n",
        completedIteration,outputIndex,sampleLength,mode,
        static_cast<unsigned long long>(Comet::Model::RandomSeed),
        PIC::nTotalThreads,globalWallSeconds,globalPeak,globalCurrent,
        minimumCapacity,maximumCapacity,globalPeakUtilization,minimumAudits,
        maximumAudits);

    for (int species=0;species<PIC::nTotalSpecies;species++) {
      char symbol[_MAX_STRING_LENGTH_PIC_];
      PIC::MolecularData::GetChemSymbol(symbol,species);
      double weight=0.0,weight2=0.0,minimum=INFINITY,maximum=0.0;
      unsigned long long count=0,occupiedCells=0;
      double minimumCellEss=INFINITY;
      const std::string cellPath=std::string("PT/diagnostics/stage9_carriers.")+
          symbol+".out="+std::to_string(outputIndex)+".csv";
      FILE *cellFile=std::fopen(cellPath.c_str(),"w");
      if (cellFile==NULL) {
        std::fprintf(stderr,"cannot open Stage-9 carrier cell output %s\n",
            cellPath.c_str());
        std::abort();
      }
      std::fprintf(cellFile,
          "cell_index,model_carriers,sum_weights,sum_squared_weights,carrier_ess,correction_min,correction_max,output_index,sampling_length,estimator_mode,random_seed,mpi_ranks\n");
      for (long int cell=0;cell<cellCount;cell++) {
        const std::size_t offset=
            static_cast<std::size_t>(cell)*PIC::nTotalSpecies+species;
        const unsigned long long cellCountValue=globalCount[offset];
        const double cellWeight=globalWeight[offset];
        const double cellWeight2=globalWeight2[offset];
        const double cellEss=cellWeight2>0.0 ?
            cellWeight*cellWeight/cellWeight2 : 0.0;
        if (cellCountValue>0) {
          occupiedCells++;
          minimumCellEss=std::min(minimumCellEss,cellEss);
          minimum=std::min(minimum,globalMinimum[offset]);
          maximum=std::max(maximum,globalMaximum[offset]);
        }
        count+=cellCountValue;
        weight+=cellWeight;
        weight2+=cellWeight2;
        std::fprintf(cellFile,
            "%ld,%llu,%.17e,%.17e,%.17e,",
            cell,cellCountValue,cellWeight,cellWeight2,cellEss);
        if (cellCountValue>0)
          std::fprintf(cellFile,"%.17e,%.17e,",globalMinimum[offset],
              globalMaximum[offset]);
        else std::fprintf(cellFile,",,");
        std::fprintf(cellFile,"%ld,%ld,%s,%llu,%d\n",outputIndex,
            sampleLength,mode,
            static_cast<unsigned long long>(Comet::Model::RandomSeed),
            PIC::nTotalThreads);
      }
      std::fclose(cellFile);
      const double ess=weight2>0.0 ? weight*weight/weight2 : 0.0;
      std::fprintf(csv,"%ld,%ld,%ld,%s,%llu,%d,%s,%llu,%.17e,%.17e,%.17e,%llu,",
          completedIteration,outputIndex,sampleLength,mode,
          static_cast<unsigned long long>(Comet::Model::RandomSeed),
          PIC::nTotalThreads,symbol,count,weight,weight2,ess,occupiedCells);
      if (occupiedCells>0)
        std::fprintf(csv,"%.17e,%.17e,%.17e,",minimumCellEss,minimum,
            maximum);
      else std::fprintf(csv,",,,");
      std::fprintf(csv,"%.17e,%ld,%ld,%ld,%ld,%.17e,%llu,0,0\n",
          globalWallSeconds,globalPeak,globalCurrent,minimumCapacity,
          maximumCapacity,globalPeakUtilization,minimumAudits);
      std::fprintf(json,
          "    {\"symbol\": \"%s\", \"model_carriers\": %llu, \"physical_inventory\": %.17e, \"sum_squared_weights\": %.17e, \"carrier_ess\": %.17e, \"occupied_cells\": %llu, \"minimum_positive_cell_ess\": ",
          symbol,count,weight,weight2,ess,occupiedCells);
      if (occupiedCells>0)
        std::fprintf(json,"%.17e, \"correction_min\": %.17e, \"correction_max\": %.17e}",
            minimumCellEss,minimum,maximum);
      else
        std::fprintf(json,"null, \"correction_min\": null, \"correction_max\": null}");
      std::fprintf(json,species+1==PIC::nTotalSpecies ? "\n" : ",\n");
    }
    std::fprintf(json,"  ]\n}\n");
    std::fclose(csv);
    std::fclose(json);
  }

  // Begin the next window from the current population; this preserves a
  // post-checkpoint baseline even when adaptation immediately changes it.
  Stage9WindowPeakLiveParticles=PIC::ParticleBuffer::NAllPart;
  Stage9WindowPeakBufferUtilization=
      static_cast<double>(PIC::ParticleBuffer::NAllPart)/
      static_cast<double>(PIC::ParticleBuffer::MaxNPart);
  Stage9WindowStartWallTime=MPI_Wtime();
  Stage9SafeBoundaryAudits=0;
}

void FillLocalCellMeasureAudit(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,long int &cellIndex,
    std::vector<double> &radius,std::vector<double> &innerRadius,
    std::vector<double> &outerRadius,std::vector<double> &measure) {
  if (node==NULL) return;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->Thread==PIC::ThisThread && node->block!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) {
        PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i,0,0));
        if (cell==NULL || !(cell->Measure>0.0) ||
            !std::isfinite(cell->Measure) ||
            !std::isfinite(cell->GetX()[0])) {
          std::fprintf(stderr,"invalid radial cell-measure audit data\n");
          std::abort();
        }
        radius[cellIndex+i]=cell->GetX()[0];
        const double width=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
        innerRadius[cellIndex+i]=node->xmin[0]+i*width;
        outerRadius[cellIndex+i]=node->xmin[0]+(i+1)*width;
        measure[cellIndex+i]=cell->Measure;
      }
    }
    cellIndex+=_BLOCK_CELLS_X_;
    return;
  }
  for (int child=0;child<(1<<DIM);child++)
    FillLocalCellMeasureAudit(node->downNode[child],cellIndex,radius,
        innerRadius,outerRadius,measure);
}

void WriteCellMeasureAudit() {
  if (!Comet::Diagnostics::Enabled) return;
  const long int cellCount=CountProfileCells(PIC::Mesh::mesh->rootTree);
  std::vector<double> localRadius(cellCount,0.0),localInner(cellCount,0.0),
      localOuter(cellCount,0.0),localMeasure(cellCount,0.0);
  std::vector<double> globalRadius(cellCount,0.0),globalInner(cellCount,0.0),
      globalOuter(cellCount,0.0),globalMeasure(cellCount,0.0);
  long int visited=0;
  FillLocalCellMeasureAudit(PIC::Mesh::mesh->rootTree,visited,localRadius,
      localInner,localOuter,localMeasure);
  if (visited!=cellCount || cellCount>INT_MAX) {
    std::fprintf(stderr,"invalid radial cell-measure audit traversal\n");
    std::abort();
  }
  MPI_Reduce(localRadius.data(),globalRadius.data(),static_cast<int>(cellCount),
      MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localMeasure.data(),globalMeasure.data(),
      static_cast<int>(cellCount),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localInner.data(),globalInner.data(),static_cast<int>(cellCount),
      MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localOuter.data(),globalOuter.data(),static_cast<int>(cellCount),
      MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    FILE *file=std::fopen("PT/diagnostics/cell_measures.csv","w");
    if (file==NULL) {
      std::fprintf(stderr,"cannot open radial cell-measure audit file\n");
      std::abort();
    }
    std::fprintf(file,
        "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3\n");
    for (long int cell=0;cell<cellCount;cell++) {
      if (!(globalMeasure[cell]>0.0) || !std::isfinite(globalRadius[cell])) {
        std::fprintf(stderr,"incomplete reduced radial cell-measure audit\n");
        std::abort();
      }
      std::fprintf(file,"%ld,%.17e,%.17e,%.17e,%.17e\n",cell,
          globalRadius[cell],globalInner[cell],globalOuter[cell],
          globalMeasure[cell]);
    }
    std::fclose(file);
  }
}

std::vector<double> EvaluateRelativeL1Changes(
    const std::vector<double> &current,const std::vector<double> &previous) {
  if (current.size()!=previous.size() ||
      current.size()%PIC::nTotalSpecies!=0) {
    std::fprintf(stderr,"inconsistent radial convergence profile sizes\n");
    std::abort();
  }
  std::vector<double> changes(PIC::nTotalSpecies,0.0);
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    double differenceNorm=0.0,currentNorm=0.0,previousNorm=0.0;
    for (size_t offset=species;offset<current.size();
         offset+=PIC::nTotalSpecies) {
      differenceNorm+=std::abs(current[offset]-previous[offset]);
      currentNorm+=std::abs(current[offset]);
      previousNorm+=std::abs(previous[offset]);
    }
    const double referenceNorm=std::max(currentNorm,previousNorm);
    changes[species]=referenceNorm>0.0 ? differenceNorm/referenceNorm : 0.0;
  }
  return changes;
}

double EvaluateRelativeL1Change(const std::vector<double> &current,
    const std::vector<double> &previous,int &limitingSpecies) {
  const std::vector<double> changes=
      EvaluateRelativeL1Changes(current,previous);
  double maximumChange=-1.0;
  limitingSpecies=0;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (changes[species]>maximumChange) {
      maximumChange=changes[species];
      limitingSpecies=species;
    }
  }

  return maximumChange;
}

void RebalanceRadialDomain(long int completedIterations) {
  // Use the core's normal block-redistribution path only after a complete
  // outer iteration.  At this point the local-timestep continuation sweeps,
  // particle exchange, collisions, chemistry, sampling, and radiative
  // cooling have all finished, so every particle belongs to an ordinary cell
  // list and the entire block (including accumulated samples and IDF state)
  // can be migrated without changing the physical solution.
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  const double start=MPI_Wtime();
  PIC::Mesh::mesh->CreateNewParallelDistributionLists(
      _PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);

  // Match the accounting performed by the core's emergency and sampling-break
  // rebalancers.  UpdateBlockTable is essential because movers, collisions,
  // and chemistry traverse this rank-local table on the next iteration.
  PIC::Parallel::IterationNumberAfterRebalancing=0;
  PIC::Parallel::CumulativeLatency=0.0;
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  PIC::Parallel::RebalancingTime=MPI_Wtime()-start;

  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: radial scheduled load balance iteration=%ld interval=%ld elapsed=%e sec\n",
        completedIterations,Comet::RadialSimulation::LoadBalanceInterval,
        PIC::Parallel::RebalancingTime);
    std::fflush(PIC::DiagnospticMessageStream);
  }
}
}

void Comet::RadialSimulation::WriteStage3CellDiagnostics(
    long int outputIndex,long int sampleLength) {
  WriteStage3CellSamplesInternal(outputIndex,sampleLength);
}

void Comet::RadialSimulation::EmergencyParticleWeightAdaptation() {
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  return;
#else
  if (!ParticleWeightAdaptationEnabled) return;
  if (!(EmergencyParticleBufferUtilization>0.0) ||
      !(EmergencyParticleBufferUtilization<1.0) ||
      !std::isfinite(EmergencyParticleBufferUtilization) ||
      PIC::ParticleBuffer::MaxNPart<=0) {
    std::fprintf(stderr,
        "invalid radial emergency particle-buffer utilization threshold\n");
    std::abort();
  }

  const long long localBefore=PIC::ParticleBuffer::NAllPart;
  long long maximumRankBefore=0;
  MPI_Allreduce(&localBefore,&maximumRankBefore,1,MPI_LONG_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  const double maximumUtilizationBefore=
      static_cast<double>(maximumRankBefore)/PIC::ParticleBuffer::MaxNPart;
  const bool bufferUtilizationTrigger=
      maximumUtilizationBefore>=EmergencyParticleBufferUtilization;
  const bool samplingPopulationTrigger=
      SamplingPopulationEmergencyPending;
  const bool allocationBurstTrigger=AllocationBurstEmergencyPending;
  if (!bufferUtilizationTrigger && !samplingPopulationTrigger &&
      !allocationBurstTrigger) return;

  // Preserve the sampling evidence for the diagnostic before clearing the
  // one-shot request below.  The request is collective: AdaptParticleWeights
  // derives it from an MPI global maximum and assigns it on every rank.
  const long long sampledMaximumBlockParticles=
      SamplingPopulationEmergencyMaximumBlockParticles;
  const long long sampledThresholdBlockParticles=
      SamplingPopulationEmergencyThresholdBlockParticles;
  const long int sampledIteration=SamplingPopulationEmergencyIteration;

  // The maximum reduction above is the only added collective on the normal
  // no-action path.  The more expensive global population accounting is used
  // solely for diagnostics after a threshold crossing has been established.
  long long globalBefore=0;
  MPI_Allreduce(&localBefore,&globalBefore,1,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);

  // The ordinary callback is reached after the mover and immediately before
  // MPI migration.  Transactional chemistry may also call this routine at a
  // synchronized pass boundary, where every particle is again linked in an
  // ordinary cell list.  In either case a rank may contain particles in
  // both owner blocks and allocated non-owner boundary copies.  A selective
  // block adaptation would require exchanging those staged counts before the
  // exchange it is meant to protect.  Instead, one emergency step doubles
  // every local block/species weight and applies 1/2 roulette to every model
  // particle.  A sampling-population request guarantees at least this one
  // step even if aggregate rank utilization is below its independent input
  // threshold.  The common factor preserves every neighboring-block and
  // inter-species W/dt ratio, so flux conversion and chemistry remain
  // unchanged in expectation while approximately half of every rank's slots
  // are released before receive buffers are allocated.
  int weightDoublings=0;
  long long globalAfter=globalBefore,maximumRankAfter=maximumRankBefore;
  double maximumUtilizationAfter=maximumUtilizationBefore;
  const double targetUtilization=allocationBurstTrigger ?
      AllocationLowWatermark : EmergencyParticleBufferUtilization;
  if (!(targetUtilization>0.0) || !std::isfinite(targetUtilization) ||
      (allocationBurstTrigger &&
       !(targetUtilization<EmergencyParticleBufferUtilization))) {
    std::fprintf(stderr,
        "invalid radial allocation-burst low-water target\n");
    std::abort();
  }
  do {
    // Every species and every allocated block copy is multiplied by the same
    // factor. Chemistry source/product W/dt ratios and systematic residual
    // means are therefore bitwise unchanged; preserving phases is the
    // explicit Stage-3 transformation for this emergency common scaling.
    EmergencyRouletteAllAllocatedBlockCopies(PIC::Mesh::mesh->rootTree);
    RefreshGlobalMinimumParticleWeights();
    weightDoublings++;

    const long long localAfter=PIC::ParticleBuffer::NAllPart;
    MPI_Allreduce(&localAfter,&globalAfter,1,MPI_LONG_LONG,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&localAfter,&maximumRankAfter,1,MPI_LONG_LONG,MPI_MAX,
        MPI_GLOBAL_COMMUNICATOR);
    maximumUtilizationAfter=
        static_cast<double>(maximumRankAfter)/PIC::ParticleBuffer::MaxNPart;

    // A factor-two roulette has a finite probability of retaining every
    // particle in a very small test population.  Repeat collectively until
    // the requested headroom is actually present rather than assuming the
    // expected one-half reduction occurred.
    if (weightDoublings>=1024 &&
        maximumUtilizationAfter>targetUtilization) {
      std::fprintf(stderr,
          "radial emergency particle-weight adaptation did not free the requested buffer capacity\n");
      std::abort();
    }
  } while (maximumUtilizationAfter>targetUtilization);

  // Count the execution, rather than merely the allocation subsystem's
  // request. These counters are reset only after a complete diagnostic
  // window is written, which makes the science gate authoritative even when
  // an event was caused by pure buffer pressure or by the sampling-population
  // alarm. Multiple simultaneous triggers still describe one event.
  Comet::Diagnostics::RecordEmergencyEvent(bufferUtilizationTrigger,
      samplingPopulationTrigger,allocationBurstTrigger);

  // Consume the sampling request only after the roulette and cache rebuild
  // have completed successfully.  A later sampling window may raise a new
  // request if its pre-adaptation block mean again exceeds the threshold.
  SamplingPopulationEmergencyPending=false;
  SamplingPopulationEmergencyMaximumBlockParticles=0;
  SamplingPopulationEmergencyThresholdBlockParticles=0;
  SamplingPopulationEmergencyIteration=-1;
  AllocationBurstEmergencyPending=false;

  static long int emergencyEvent=0;
  emergencyEvent++;
  if (PIC::ThisThread==0) {
    const char *trigger=allocationBurstTrigger ?
        (bufferUtilizationTrigger ? "allocation-burst+buffer" :
            "allocation-burst") :
        (bufferUtilizationTrigger ?
            (samplingPopulationTrigger ? "buffer+sampling-population" :
                "buffer") : "sampling-population");
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: radial emergency particle-weight adaptation event=%ld trigger=%s high-watermark=%.6f target=%.6f max-buffer-utilization=%.6e->%.6e particles=%lld->%lld max-rank-particles=%lld->%lld weight-doublings=%d sampling-iteration=%ld sampling-max-mean-total-particles-per-cell=%.6e sampling-threshold=%.6e\n",
        emergencyEvent,trigger,EmergencyParticleBufferUtilization,
        targetUtilization,maximumUtilizationBefore,maximumUtilizationAfter,globalBefore,
        globalAfter,maximumRankBefore,maximumRankAfter,weightDoublings,
        samplingPopulationTrigger ? sampledIteration : -1,
        samplingPopulationTrigger ?
            static_cast<double>(sampledMaximumBlockParticles)/
                BlockCellCount : 0.0,
        samplingPopulationTrigger ?
            static_cast<double>(sampledThresholdBlockParticles)/
                BlockCellCount : 0.0);
    std::fflush(PIC::DiagnospticMessageStream);
  }
#endif
}

void Comet::RadialSimulation::RequestAllocationBurstEmergency() {
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  std::fprintf(stderr,
      "allocation-burst emergency requires species-local block weights\n");
  std::abort();
#else
  // This function is collective by contract.  Setting the request on every
  // rank before entering the existing emergency routine guarantees identical
  // control flow even when only one rank reported capacity pressure.
  AllocationBurstEmergencyPending=true;
  Comet::Diagnostics::RecordEmergencyRequest();
  EmergencyParticleWeightAdaptation();
#endif
}

bool Comet::RadialSimulation::Validate(std::string *error) {
  if (!Enabled) return true;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ != _PIC_MODE_ON_
  return Fail("individual corrections in the radial backend require local-timestep flux resampling",error);
#endif
#if _PIC__PARTICLE_COLLISION_MODEL__MODE_ == _PIC_MODE_ON_ && \
    _PIC__PARTICLE_COLLISION_MODEL_ != _PIC__PARTICLE_COLLISION_MODEL__NTC_
  // Stage 4 validates the correction-aware general NTC path. MF and
  // user-defined collision executors have different pair samplers and must
  // not silently inherit that claim until they receive their own audit.
  return Fail("individual corrections with collisions require the validated general NTC model",error);
#endif
#endif
  if (IndividualCorrectionWindowsEnabled) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    return Fail("individual correction windows require individual particle-weight corrections",error);
#else
#if _COMPILATION_MODE_ != _COMPILATION_MODE__MPI_
    return Fail("individual correction windows currently require the MPI-only host executor",error);
#endif
    if (MaximumIndividualCorrectionWindowSplits<2)
      return Fail("individual correction windows require a split limit of at least two",error);
    if (!(AllocationReserveFraction>0.0) ||
        !(AllocationReserveFraction<1.0) ||
        !std::isfinite(AllocationReserveFraction))
      return Fail("individual correction windows require a valid allocation reserve fraction",error);
    try {
      const std::vector<Comet::CorrectionWindowPlanner::cPolicy> policies=
          BuildCorrectionWindowPolicies();
      for (int species=0;species<PIC::nTotalSpecies;species++) {
        Comet::CorrectionWindowPlanner::ValidatePolicy(policies[species]);
        // Ordered bounds alone are insufficient for the fixed
        // ceil(c/target) split rule.  Prove that the first representable value
        // above the upper edge can be returned to the window; the input
        // parser applies the same compatibility gate.
        Comet::CorrectionWindowPlanner::BuildParticlePlan(
            Comet::CorrectionWindowPlanner::cRequest(species,
                std::nextafter(IndividualCorrectionWindowUpper[species],
                    std::numeric_limits<double>::infinity())),
            policies[species]);
        if (IndividualCorrectionWindowPolicySource[species]<
                PolicySource::Global ||
            IndividualCorrectionWindowPolicySource[species]>
                PolicySource::ExplicitSpecies)
          return Fail("an individual correction-window policy source is invalid",error);
        if (HasConfiguredDirectSource(species) &&
            !(IndividualCorrectionWindowLower[species]<=1.0 &&
              IndividualCorrectionWindowUpper[species]>=1.0))
          return Fail("a direct-source correction window must contain the injection correction one",error);
      }
    }
    catch (const std::exception &exception) {
      return Fail(std::string("invalid individual correction window: ")+
          exception.what(),error);
    }
#endif
  }
  if (!(CourantNumber>0.0) || CourantNumber>1.0 ||
      !std::isfinite(CourantNumber))
    return Fail("radial Courant number must be in (0,1]",error);
  if (!(InnerBoundaryTemperature>0.0) ||
      !std::isfinite(InnerBoundaryTemperature))
    return Fail("radial inner-boundary temperature must be positive",error);
  if (MaxIterations<=0)
    return Fail("radial maximum iteration count must be positive",error);
  if (SamplingLength<=0)
    return Fail("radial sampling length must be positive",error);
  if (SamplingLength>MaxIterations)
    return Fail("radial sampling length cannot exceed the maximum iteration count",error);
  if (LoadBalanceInterval<=0)
    return Fail("radial load-balance interval must be positive",error);
  if (LoadBalanceInterval>SamplingLength)
    return Fail("radial load-balance interval cannot exceed the sampling length",error);
  if (ParticleWeightAdaptationEnabled) {
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
    return Fail("particle-weight adaptation requires SpeciesLocalParticleWeight",error);
#else
    if (MinimumParticlesPerCell<=0 || MaximumParticlesPerCell<=0)
      return Fail("particle-weight adaptation range must be positive",error);
    if (MaximumParticlesPerCell<2*MinimumParticlesPerCell)
      return Fail("particle-weight adaptation maximum must be at least twice the minimum",error);
    if (ParticleWeightAdaptationInterval==0 ||
        ParticleWeightAdaptationInterval < -1)
      return Fail("particle-weight adaptation interval must be positive or use the legacy -1 sampling-break value",error);
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      if (SpeciesMinimumParticlesPerCell[species]<=0 ||
          SpeciesMaximumParticlesPerCell[species]<
              2*SpeciesMinimumParticlesPerCell[species])
        return Fail("a resolved species population range is invalid",error);
      if (ParticleWeightLocked[species]!=0 &&
          ParticleWeightLocked[species]!=1)
        return Fail("a resolved species weight-lock flag is invalid",error);
      if (PopulationPolicySource[species]<PolicySource::Global ||
          PopulationPolicySource[species]>PolicySource::ExplicitSpecies)
        return Fail("a resolved population-policy source is invalid",error);
      if (PopulationClass[species]<PolicyClass::Other ||
          PopulationClass[species]>PolicyClass::ChemicalOnly)
        return Fail("a resolved population class is invalid",error);
      const bool direct=HasConfiguredDirectSource(species);
      if (direct!=(PopulationClass[species]==PolicyClass::DirectSource))
        return Fail("resolved direct-source population class disagrees with configured sources",error);
      if (direct && PopulationClass[species]==PolicyClass::ChemicalOnly)
        return Fail("a direct source cannot be classified as chemical-only",error);
    }
    if (!(EmergencyParticleBufferUtilization>0.0) ||
        !(EmergencyParticleBufferUtilization<1.0) ||
        !std::isfinite(EmergencyParticleBufferUtilization))
      return Fail("emergency particle-buffer utilization must be in (0,1)",error);
    if (!(AllocationLowWatermark>0.0) ||
        !(AllocationLowWatermark<EmergencyParticleBufferUtilization) ||
        !std::isfinite(AllocationLowWatermark))
      return Fail("allocation low watermark must be positive and below the emergency high watermark",error);
    if (!(AllocationReserveFraction>0.0) ||
        !(AllocationReserveFraction<1.0) ||
        !std::isfinite(AllocationReserveFraction) ||
        EmergencyParticleBufferUtilization>1.0-AllocationReserveFraction)
      return Fail("allocation reserve must leave the emergency high watermark below the transaction ceiling",error);
    if (MaximumFaceTransferMultiplicity!=1)
      return Fail("the validated face-allocation controller currently requires maximum face multiplicity 1",error);
    if (Comet::Chemistry::Enabled &&
        MaximumPhotochemicalProductsPerReaction<=0)
      return Fail("particle-weight adaptation with chemistry requires a positive maximum photochemical product count per reaction",error);
#endif
  }
  else if (ParticleWeightAdaptationInterval>0) {
    return Fail("particle-weight adaptation interval requires adaptation enabled",error);
  }
  if (!(ConvergenceTolerance>0.0) || ConvergenceTolerance>1.0 ||
      !std::isfinite(ConvergenceTolerance))
    return Fail("radial convergence tolerance must be in (0,1]",error);
  if (RequiredConvergedWindows<=0)
    return Fail("required converged-window count must be positive",error);
  if (Stage9ScienceDiagnosticsEnabled) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    return Fail("Stage-9 science diagnostics require individual particle corrections",error);
#else
    if (!Comet::Diagnostics::Enabled)
      return Fail("Stage-9 science diagnostics require comet diagnostics",error);
    if (!Comet::Chemistry::Enabled)
      return Fail("Stage-9 science diagnostics require comet chemistry",error);
#endif
  }

  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!(CharacteristicSpeed[species]>0.0) ||
        !std::isfinite(CharacteristicSpeed[species])) {
      char message[160];
      std::snprintf(message,sizeof(message),
          "radial characteristic speed is missing for species %d",species);
      return Fail(message,error);
    }
  }
  return true;
}

double Comet::RadialSimulation::LocalTimeStep(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  // Restrict one move to the selected fraction of the local radial cell size.
  // Each species has its own explicitly supplied speed so light, energetic
  // photoproducts cannot inherit the much slower parent-gas timestep.
  return CourantNumber*node->GetCharacteristicCellSize()/
      CharacteristicSpeed[species];
}

double Comet::RadialSimulation::LocalBlockInjectionRate(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  // A uniform spherical source has one global rate Q. In DIM1 the inner box
  // face represents the entire 4*pi*R^2 shell, so Q is assigned once to the
  // block touching face 0; no area division is needed here.
  return HasInnerFace(node) ? Comet::Model::SurfaceProductionRate[species] : 0.0;
}

bool Comet::RadialSimulation::InjectionBlockIndicator(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (!HasInnerFace(node)) return false;
  for (int species=0;species<PIC::nTotalSpecies;species++)
    if (Comet::Model::SurfaceProductionRate[species]>0.0) return true;
  return false;
}

long int Comet::RadialSimulation::InjectBlock(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (!HasInnerFace(node)) return 0;

  long int injected=0;
  const double externalNormal[3]={-1.0,0.0,0.0};
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    const double sourceRate=Comet::Model::SurfaceProductionRate[species];
    if (!(sourceRate>0.0)) continue;

    const double particleWeight=node->block->GetLocalParticleWeight(species);
    const double timeStep=node->block->GetLocalTimeStep(species);
    const double modelRate=sourceRate/particleWeight;
    double eventTime=0.0;
    const double bulkVelocity[3]={
        Comet::Model::SurfaceBulkSpeed[species],0.0,0.0};

    // A Poisson process avoids systematic rounding when Q*dt/weight is not an
    // integer. The remaining fraction of the current step is passed to the
    // mover so injected particles are uniformly distributed in injection time.
    while ((eventTime+=-std::log(rnd())/modelRate)<timeStep) {
      const long int ptr=PIC::ParticleBuffer::GetNewParticle();
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      double position[3]={node->xmin[0]+PIC::Mesh::mesh->EPS,0.0,0.0};
      double velocity[3];
      PIC::Distribution::InjectMaxwellianDistribution(velocity,bulkVelocity,
          Comet::Model::SurfaceTemperature[species],
          const_cast<double*>(externalNormal),species,-1);

      PIC::ParticleBuffer::SetX(position,data);
      PIC::ParticleBuffer::SetV(velocity,data);
      PIC::ParticleBuffer::SetI(species,data);
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
      // A surface molecule leaves the thermalized Knudsen layer with its
      // configured translational and internal temperatures in equilibrium.
      // AMPS owns the Larsen-Borgnakke sampling of the internal energies.
      PIC::IDF::InitRotTemp(Comet::Model::SurfaceTemperature[species],data);
      PIC::IDF::InitVibTemp(Comet::Model::SurfaceTemperature[species],data);
#endif
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,data);
#endif
      // One injected model particle carries W/dt molecules per second through
      // one outer iteration.  Accumulating that carrier rate and dividing by
      // the completed window length gives the realized physical source rate.
      Comet::Diagnostics::RecordSurfaceInjection(species,
          particleWeight/timeStep);
      PIC::Mover::DispatchParticleMover(ptr,timeStep-eventTime,node);
      injected++;
    }
  }
  return injected;
}

int Comet::RadialSimulation::ProcessDomainBoundary(long int ptr,
    double *position,double *velocity,int intersectionFace,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  (void)position;
  const int species=PIC::ParticleBuffer::GetI(ptr);
  const double representedRate=ParticlePhysicalRate(ptr,species,node);
  if (intersectionFace==1) {
    // Rubin et al. remove every particle crossing the outer boundary.
    Comet::Diagnostics::RecordOuterBoundary(species,representedRate);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }
  if (intersectionFace!=0) {
    std::fprintf(stderr,"unexpected DIM1 comet boundary face %d\n",
        intersectionFace);
    std::abort();
  }

  // Gas returning to the top of the Knudsen layer is re-emitted. The paper
  // states that the gas is completely thermalized there, so discard the
  // incident velocity and draw a new outward half-Maxwellian at the explicit
  // inner-boundary temperature. Dust is outside this radial gas backend.
  Comet::Diagnostics::RecordInnerBoundary(species,representedRate);
  const double externalNormal[3]={-1.0,0.0,0.0};
  const double bulkVelocity[3]={0.0,0.0,0.0};
  PIC::Distribution::InjectMaxwellianDistribution(velocity,bulkVelocity,
      InnerBoundaryTemperature,const_cast<double*>(externalNormal),species,-1);
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  PIC::ParticleBuffer::byte *data=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::IDF::InitRotTemp(InnerBoundaryTemperature,data);
  PIC::IDF::InitVibTemp(InnerBoundaryTemperature,data);
#endif
  (void)node;
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

void Comet::RadialSimulation::Initialize() {
  // Sampling storage and cell blocks must exist before AMPS initializes
  // weights, timesteps, collision data, or its particle mover.
  PIC::Mesh::initCellSamplingDataBuffer();
  // The optional systematic chemistry accumulator is block-owned state.
  // Register it after the core's local W/dt fields but before block allocation
  // so ordinary AMR migration and restart include the complete state vector.
  Comet::Chemistry::RequestBlockData();
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  // PIC::TimeStep() normally builds this owner-block table at the beginning
  // of the first iteration.  Build it here as well because focused validators
  // run after Initialize() but before TimeStep(), and because every sampling-
  // break adaptation must traverse the same authoritative owner list used by
  // the mover and collision modules.
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  PIC::Mesh::mesh->InitCellMeasure();

  PIC::Mover::ProcessOutsideDomainParticles=ProcessDomainBoundary;
  // The transactional chemistry byte must be appended before
  // PIC::Init_AfterParser initializes the fixed-length particle buffer.
  // Other comet modes that disable chemistry reserve nothing.
  Comet::Chemistry::RequestParticleData();
  PIC::Init_AfterParser();
  // Stage 8 installs an all-or-none, application-owned trajectory-event
  // callback set. Installation must precede Mover::Init(), where the generic
  // core validates the contract even for a zero-particle validation run.
  Comet::Chemistry::ConfigureTotalReactionForcingMover();
  PIC::Mover::Init();
  Comet::Diagnostics::Initialize();
  WriteCellMeasureAudit();

  if (Comet::Chemistry::Enabled) {
    // Install the complete generic core batch interface as one unit.  The
    // core rejects a partial callback set; non-comet applications leave all
    // pointers NULL and retain the historical single photochemistry pass.
    PIC::ChemicalReactions::PhotolyticReactions::BeginBatchedModel=
        Comet::Chemistry::BeginBatchedModel;
    PIC::ChemicalReactions::PhotolyticReactions::IsBatchedParticleEligible=
        Comet::Chemistry::IsBatchedParticleEligible;
    PIC::ChemicalReactions::PhotolyticReactions::GetBatchedParticleLimit=
        Comet::Chemistry::GetBatchedParticleLimit;
    PIC::ChemicalReactions::PhotolyticReactions::BatchedReactionProcessor=
        Comet::Chemistry::ProcessParticle;
    PIC::ChemicalReactions::PhotolyticReactions::FinishBatchedPass=
        Comet::Chemistry::FinishBatchedPass;
  }

  // The core invokes this collective hook immediately before each ordinary
  // particle migration, including local-timestep continuation exchanges.
  // Installing it here keeps the generic core inactive for other AMPS cases.
  PIC::Parallel::BeforeParticleDataExchange=
      ParticleWeightAdaptationEnabled ? EmergencyParticleWeightAdaptation : NULL;

  // Use the case-controlled window length for AMPS macroscopic sampling.
  // Convergence comparisons are made only after a complete window; a partial
  // final window at the hard iteration ceiling is never declared converged.
  PIC::RequiredSampleLength=SamplingLength;

  PIC::ParticleWeightTimeStep::LocalTimeStep=LocalTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();
  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=
      LocalBlockInjectionRate;

  // Direct parent sources determine their surface statistical weights.
  // Daughter species have no injection rate, so their W/dt carrier rate is
  // tied to the largest initialized parent carrier rate.  This prevents a
  // single parent event from requiring a large product-clone burst while the
  // chemistry executor still conserves the physical production rate through
  // explicit residual split/roulette whenever species carrier rates differ.
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *innerNode=
      FindInnerBoundaryNode(PIC::Mesh::mesh->rootTree);

  // After domain decomposition, only the rank that owns the block adjacent to
  // the nucleus has that block allocated. Confirm that there is exactly one
  // owner globally, then broadcast surface values with an all-reduction. The
  // previous per-rank NULL check aborted every valid rank that did not own the
  // nucleus-adjacent block.
  const int localInnerOwner=(innerNode!=NULL) ? 1 : 0;
  int globalInnerOwners=0;
  MPI_Allreduce(&localInnerOwner,&globalInnerOwners,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  if (globalInnerOwners!=1) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "radial comet mesh has %d allocated inner-boundary blocks; expected exactly one\n",
          globalInnerOwners);
    std::abort();
  }

  // Start each direct parent with AMPS' source-derived constant weight.  Use
  // its value at the surface as the normalization, then scale W with the
  // local timestep.  Global weights would make the number of particles per
  // cell grow in proportion to the logarithmic cell width and render the
  // 10^6 km Halley domain impractical at sub-metre surface resolution.
  double daughterWeightOverTimeStep=0.0;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (Comet::Model::SurfaceProductionRate[species]>0.0) {
      PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(species);

      // Only the owning rank reads the local block. Summation acts as a
      // broadcast because the ownership check above guarantees one source.
      const double localSurfaceValues[2]={
          innerNode!=NULL ? innerNode->block->GetLocalParticleWeight(species) : 0.0,
          innerNode!=NULL ? innerNode->block->GetLocalTimeStep(species) : 0.0};
      double surfaceValues[2]={0.0,0.0};
      MPI_Allreduce(localSurfaceValues,surfaceValues,2,MPI_DOUBLE,MPI_SUM,
          MPI_GLOBAL_COMMUNICATOR);
      const double surfaceWeight=surfaceValues[0];
      const double surfaceTimeStep=surfaceValues[1];
      if (!(surfaceWeight>0.0) || !std::isfinite(surfaceWeight) ||
          !(surfaceTimeStep>0.0) || !std::isfinite(surfaceTimeStep)) {
        if (PIC::ThisThread==0)
          std::fprintf(stderr,
              "invalid surface weight or timestep for radial species %d\n",
              species);
        std::abort();
      }

      const double weightOverTimeStep=surfaceWeight/surfaceTimeStep;
      if (!(weightOverTimeStep>0.0) || !std::isfinite(weightOverTimeStep)) {
        if (PIC::ThisThread==0)
          std::fprintf(stderr,
              "invalid surface W/dt carrier rate for radial species %d\n",
              species);
        std::abort();
      }
      SetWeightOverTimeStep(species,weightOverTimeStep,
          PIC::Mesh::mesh->rootTree);
      if (weightOverTimeStep>daughterWeightOverTimeStep)
        daughterWeightOverTimeStep=weightOverTimeStep;
    }
  }
  if (!(daughterWeightOverTimeStep>0.0) ||
      !std::isfinite(daughterWeightOverTimeStep)) {
    std::fprintf(stderr,"radial comet run has no direct surface source\n");
    std::abort();
  }

  // Give every non-source daughter the same conservative W/dt carrier rate.
  // Its actual block weight remains species- and block-dependent because it
  // is this carrier rate multiplied by that species' local block timestep.
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!(Comet::Model::SurfaceProductionRate[species]>0.0)) {
      SetWeightOverTimeStep(species,daughterWeightOverTimeStep,
          PIC::Mesh::mesh->rootTree);
    }
  }
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  double daughterWeight=0.0;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (Comet::Model::SurfaceProductionRate[species]>0.0) {
      PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(species);
      if (PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]>daughterWeight)
        daughterWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[species];
    }
  }
  if (!(daughterWeight>0.0) || !std::isfinite(daughterWeight)) {
    std::fprintf(stderr,"radial comet run has no direct surface source\n");
    std::abort();
  }

  // A daughter with no direct surface source needs a statistical weight even
  // though the core injection-based initializer cannot derive one.  Use the
  // largest direct-parent weight.  Consequently a chemistry event never has
  // to create an unbounded burst merely because the target species inherited
  // the smallest source weight; products are retained by unbiased roulette
  // when their parent represents fewer molecules than this target weight.
  // This choice changes variance, not the expected number/mass/momentum/
  // energy source rates, which CometChemistry converts by
  // (W_parent/dt_parent)/(W_product/dt_product).
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!(PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]>0.0))
      PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(
          species,daughterWeight);
  }
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SINGLE_GLOBAL_PARTICLE_WEIGHT_
  // Compute the source-based candidate weight for every directly emitted
  // parent.  cDataBlockAMR::SetLocalParticleWeight() stores the minimum in
  // GlobalParticleWeight[0] in single-global mode.  Selecting the minimum is
  // intentional: the least represented source still receives the requested
  // ReferenceInjectionParticleNumber, while stronger sources receive more
  // particles and every chemical daughter retains exactly the same physical
  // weight without an individual correction field.
  bool foundDirectSource=false;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (Comet::Model::SurfaceProductionRate[species]>0.0) {
      PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(species);
      foundDirectSource=true;
    }
  }
  if (!foundDirectSource ||
      !(PIC::ParticleWeightTimeStep::GlobalParticleWeight[0]>0.0) ||
      !std::isfinite(PIC::ParticleWeightTimeStep::GlobalParticleWeight[0])) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "radial comet run could not initialize a common global particle weight\n");
    std::abort();
  }
#else
#error "The radial comet backend requires a supported global or local particle-weight mode"
#endif

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  if (ParticleWeightAdaptationEnabled) {
    // Establish every transport and chemistry allocation invariant before
    // the first source particle is injected.  Counts are still zero here,
    // but the planner treats configured surface species as graph roots and
    // follows all enabled daughter reactions transitively.  Thus iteration
    // one cannot encounter an unsafe source-to-source carrier ratio while
    // waiting for the first 500-step sampling break.  Empty blocks remain
    // empty because every resulting plan is a weight-only increase.
    AdaptParticleWeights(0);
  }
  ReportResolvedPopulationPolicies(innerNode);
  ReportResolvedCorrectionWindowPolicies();
#endif

  PIC::BC::BlockInjectionBCindicatior=InjectionBlockIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=InjectBlock;
  PIC::BC::InitBoundingBoxInjectionBlockList();
}

bool Comet::RadialSimulation::ValidateIndividualCorrectionWindows() {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  if (PIC::ThisThread==0)
    std::fprintf(stderr,
        "correction-window validation requires individual corrections\n");
  return false;
#else
  namespace Window=Comet::CorrectionWindowPlanner;
  if (!IndividualCorrectionWindowsEnabled) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "correction-window validation requires the input switch ON\n");
    return false;
  }

  // The focused validator owns all particle lists and is selected instead of
  // a science run.  One owner block per rank exercises the collective commit
  // and rank-local capacity calculation without manufacturing ghost copies.
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
      SelectValidationOwnerBlock();
  const int localHasBlock=(node!=NULL && node->block!=NULL) ? 1 : 0;
  int globalBlockCount=0;
  MPI_Allreduce(&localHasBlock,&globalBlockCount,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  if (globalBlockCount!=PIC::nTotalThreads) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "correction-window validation requires one owner block per MPI rank\n");
    return false;
  }

  int species=-1;
  for (int candidate=0;candidate<PIC::nTotalSpecies;candidate++)
    if (!HasConfiguredDirectSource(candidate)) {
      species=candidate;
      break;
    }
  if (species<0) species=0;
  const Window::cPolicy policy(
      IndividualCorrectionWindowLower[species],
      IndividualCorrectionWindowTarget[species],
      IndividualCorrectionWindowUpper[species],
      static_cast<std::size_t>(MaximumIndividualCorrectionWindowSplits));

  // The first value above the inclusive upper edge is the hardest place to
  // form two legal children.  Reject a mathematically ordered but unusable
  // policy before constructing any synthetic particle.
  const double highCorrection=std::nextafter(policy.UpperCorrection,
      std::numeric_limits<double>::infinity());
  Window::cParticlePlan highPlan;
  try {
    highPlan=Window::BuildParticlePlan(
        Window::cRequest(species,highCorrection),policy);
  }
  catch (const std::exception &exception) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "correction-window validation cannot split above the upper edge: %s\n",
          exception.what());
    return false;
  }
  if (highPlan.Action!=Window::eAction::Split ||
      highPlan.OutputCorrections.size()<2) return false;

  struct cAudit {
    long long Count=0;
    long long Outside=0;
    double Sum=0.0;
    double SquareSum=0.0;
  };
  auto auditSpecies = [&](int auditSpecies) -> cAudit {
    cAudit local,global;
    for (unsigned int blockIndex=0;
         blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;
         blockIndex++) {
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *auditNode=
          PIC::DomainBlockDecomposition::BlockTable[blockIndex];
      if (auditNode==NULL || auditNode->block==NULL) continue;
      for (long int cell=0;cell<BlockCellCount;cell++)
        for (long int ptr=auditNode->block->FirstCellParticleTable[cell];
             ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
          PIC::ParticleBuffer::byte *data=
              PIC::ParticleBuffer::GetParticleDataPointer(ptr);
          if (PIC::ParticleBuffer::GetI(data)!=auditSpecies) continue;
          const double correction=
              PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
          local.Count++;
          local.Sum+=correction;
          local.SquareSum+=correction*correction;
          if (!(correction>=IndividualCorrectionWindowLower[auditSpecies]) ||
              !(correction<=IndividualCorrectionWindowUpper[auditSpecies]) ||
              !std::isfinite(correction)) local.Outside++;
        }
    }
    long long localCounts[2]={local.Count,local.Outside};
    long long globalCounts[2]={0,0};
    MPI_Allreduce(localCounts,globalCounts,2,MPI_LONG_LONG,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    double localSums[2]={local.Sum,local.SquareSum};
    double globalSums[2]={0.0,0.0};
    MPI_Allreduce(localSums,globalSums,2,MPI_DOUBLE,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    global.Count=globalCounts[0];
    global.Outside=globalCounts[1];
    global.Sum=globalSums[0];
    global.SquareSum=globalSums[1];
    return global;
  };
  auto setCorrection = [&](double correction) {
    if (!localHasBlock) return;
    for (long int ptr=node->block->FirstCellParticleTable[0];ptr!=-1;
         ptr=PIC::ParticleBuffer::GetNext(ptr))
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(correction,
          PIC::ParticleBuffer::GetParticleDataPointer(ptr));
  };

  bool pass=true;
  int keepChecks=0,splitChecks=0,rouletteChecks=0;
  int capacityRollbackChecks=0,restartChecks=0,directSourceChecks=0;

  // Inclusive boundaries and the target must be bitwise identity operations.
  PIC::ParticleBuffer::DeleteAllParticles();
  SeedValidationParticles(node,species,3);
  if (localHasBlock) {
    const double values[3]={policy.LowerCorrection,
        policy.TargetCorrection,policy.UpperCorrection};
    int index=0;
    for (long int ptr=node->block->FirstCellParticleTable[0];ptr!=-1;
         ptr=PIC::ParticleBuffer::GetNext(ptr))
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(values[index++],
          PIC::ParticleBuffer::GetParticleDataPointer(ptr));
  }
  const cValidationMoments keepBefore=CollectValidationMoments(species);
  cCorrectionWindowStatistics keepStatistics;
  const bool keepAccepted=ApplyIndividualCorrectionWindowsInternal(0,-1,
      &keepStatistics,false);
  const cValidationMoments keepAfter=CollectValidationMoments(species);
  const cAudit keepAudit=auditSpecies(species);
  const bool keepPass=keepAccepted && keepAudit.Outside==0 &&
      keepAudit.Count==3LL*globalBlockCount &&
      ValidationMomentsEqual(keepBefore,keepAfter,0.0) &&
      keepStatistics.Kept==3 && keepStatistics.RouletteCandidates==0 &&
      keepStatistics.SplitParents==0;
  pass=pass && keepPass;
  if (keepPass) keepChecks++;

  // A zero-capacity preflight of a required split must be a true rollback:
  // particle count, corrections, and all represented physical moments remain
  // unchanged, and roulette RNG has not yet been consulted.
  PIC::ParticleBuffer::DeleteAllParticles();
  SeedValidationParticles(node,species,1);
  setCorrection(highCorrection);
  const long int capacityCountBefore=PIC::ParticleBuffer::NAllPart;
  const cValidationMoments capacityBefore=CollectValidationMoments(species);
  const cAudit capacityAuditBefore=auditSpecies(species);
  cCorrectionWindowStatistics rejectedStatistics;
  const bool capacityAccepted=ApplyIndividualCorrectionWindowsInternal(0,0,
      &rejectedStatistics,false);
  const cValidationMoments capacityAfter=CollectValidationMoments(species);
  const cAudit capacityAuditAfter=auditSpecies(species);
  const bool capacityPass=!capacityAccepted &&
      PIC::ParticleBuffer::NAllPart==capacityCountBefore &&
      capacityAuditAfter.Count==capacityAuditBefore.Count &&
      capacityAuditAfter.Sum==capacityAuditBefore.Sum &&
      capacityAuditAfter.SquareSum==capacityAuditBefore.SquareSum &&
      ValidationMomentsEqual(capacityBefore,capacityAfter,0.0);
  pass=pass && capacityPass;
  if (capacityPass) capacityRollbackChecks++;

  // Commit the same split with normal capacity.  Exact child corrections and
  // CloneParticle preserve represented number, mass, momentum, kinetic
  // energy, internal energy, chemistry stage, and zero mover remainder.
  const cValidationMoments splitBefore=CollectValidationMoments(species);
  cCorrectionWindowStatistics splitStatistics;
  const bool splitAccepted=ApplyIndividualCorrectionWindowsInternal(0,-1,
      &splitStatistics,false);
  const cValidationMoments splitAfter=CollectValidationMoments(species);
  const cAudit splitAudit=auditSpecies(species);
  const long long expectedSplitCount=
      static_cast<long long>(highPlan.OutputCorrections.size())*
      globalBlockCount;
  const bool splitPass=splitAccepted && splitAudit.Outside==0 &&
      splitAudit.Count==expectedSplitCount &&
      ValidationMomentsEqual(splitBefore,splitAfter,5.0e-13) &&
      splitStatistics.SplitParents==1 &&
      splitStatistics.SplitChildren==
          static_cast<long long>(highPlan.OutputCorrections.size()-1);
  pass=pass && splitPass;
  if (splitPass) splitChecks++;

  // Restart serializes the full particle record, including each child
  // correction and internal state.  Tracker IDs and RNG state are purposely
  // excluded from this assertion because the core reinitializes tracking on
  // reload and does not checkpoint the ordinary roulette stream.
  const cValidationMoments restartBefore=CollectValidationMoments(species);
  const cAudit restartAuditBefore=auditSpecies(species);
  char restartName[256];
  std::snprintf(restartName,sizeof(restartName),
      "/tmp/amps-hartley2-stage5.ranks=%d.restart",PIC::nTotalThreads);
  PIC::Restart::SaveParticleData(restartName);
  PIC::ParticleBuffer::DeleteAllParticles();
  PIC::Restart::ReadParticleData(restartName);
  const cValidationMoments restartAfter=CollectValidationMoments(species);
  const cAudit restartAuditAfter=auditSpecies(species);
  const bool restartPass=restartAuditAfter.Count==restartAuditBefore.Count &&
      restartAuditAfter.Outside==0 &&
      restartAuditAfter.Sum==restartAuditBefore.Sum &&
      restartAuditAfter.SquareSum==restartAuditBefore.SquareSum &&
      ValidationMomentsEqual(restartBefore,restartAfter,0.0);
  pass=pass && restartPass;
  if (restartPass) restartChecks++;
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) std::remove(restartName);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  // Low-side roulette is assessed against its analytic mean and variance.
  // All seeded particles have identical state, while base weights may differ
  // by radial block; WeightSquared retains the exact unequal-weight variance
  // coefficient needed for a decomposition-independent physical-moment test.
  PIC::ParticleBuffer::DeleteAllParticles();
  constexpr long int RouletteParticlesPerRank=65536;
  SeedValidationParticles(node,species,RouletteParticlesPerRank);
  const double lowCorrection=0.5*policy.LowerCorrection;
  setCorrection(lowCorrection);
  const cValidationMoments rouletteBefore=CollectValidationMoments(species);
  rnd_seed(49979687+PIC::ThisThread);
  cCorrectionWindowStatistics rouletteStatistics;
  const bool rouletteAccepted=ApplyIndividualCorrectionWindowsInternal(0,-1,
      &rouletteStatistics,false);
  const cValidationMoments rouletteAfter=CollectValidationMoments(species);
  const cAudit rouletteAudit=auditSpecies(species);
  long long localRoulette[2]={rouletteStatistics.RouletteCandidates,
      rouletteStatistics.RouletteSurvivors};
  long long globalRoulette[2]={0,0};
  MPI_Allreduce(localRoulette,globalRoulette,2,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  const double survivalProbability=
      lowCorrection/policy.TargetCorrection;
  const double expectedSurvivors=
      RouletteParticlesPerRank*globalBlockCount*survivalProbability;
  const double survivorSigma=std::sqrt(
      RouletteParticlesPerRank*globalBlockCount*survivalProbability*
      (1.0-survivalProbability));
  const double physicalSigma=std::sqrt(rouletteBefore.WeightSquared*
      (policy.TargetCorrection-lowCorrection)/lowCorrection);
  const double physicalTolerance=rouletteBefore.Number>0.0 ?
      7.0*physicalSigma/rouletteBefore.Number+5.0e-13 : 0.0;
  const bool roulettePass=rouletteAccepted && rouletteAudit.Outside==0 &&
      globalRoulette[0]==
          static_cast<long long>(RouletteParticlesPerRank)*globalBlockCount &&
      std::fabs(globalRoulette[1]-expectedSurvivors)<=7.0*survivorSigma &&
      std::fabs(rouletteAfter.Number-rouletteBefore.Number)<=
          7.0*physicalSigma &&
      ValidationMomentsEqual(rouletteAfter,rouletteBefore,
          physicalTolerance);
  pass=pass && roulettePass;
  if (roulettePass) rouletteChecks++;

  // Direct injection explicitly initializes c=1.  The input resolver rejects
  // a direct-source window that excludes one; exercise that identity through
  // the production controller rather than relying on parser evidence alone.
  int directSpecies=-1;
  for (int candidate=0;candidate<PIC::nTotalSpecies;candidate++)
    if (HasConfiguredDirectSource(candidate)) {
      directSpecies=candidate;
      break;
    }
  PIC::ParticleBuffer::DeleteAllParticles();
  if (directSpecies>=0) {
    SeedValidationParticles(node,directSpecies,1);
    if (localHasBlock) {
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(
              node->block->FirstCellParticleTable[0]);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,data);
    }
    cCorrectionWindowStatistics directStatistics;
    const bool directAccepted=ApplyIndividualCorrectionWindowsInternal(0,-1,
        &directStatistics,false);
    const cAudit directAudit=auditSpecies(directSpecies);
    const bool directPass=directAccepted &&
        directAudit.Count==globalBlockCount && directAudit.Outside==0 &&
        directAudit.Sum==static_cast<double>(globalBlockCount) &&
        directStatistics.Kept==1;
    pass=pass && directPass;
    if (directPass) directSourceChecks++;
  }
  else pass=false;

  // An empty species/list remains empty; windows never manufacture a
  // trajectory without a represented input particle.
  PIC::ParticleBuffer::DeleteAllParticles();
  cCorrectionWindowStatistics emptyStatistics;
  const bool emptyPass=ApplyIndividualCorrectionWindowsInternal(0,-1,
      &emptyStatistics,false) && PIC::ParticleBuffer::NAllPart==0 &&
      emptyStatistics.Examined==0;
  pass=pass && emptyPass;

  const int localResult=pass ? 1 : 0;
  int globalResult=0;
  MPI_Allreduce(&localResult,&globalResult,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "%s individual-correction window validation correction-mode=ON mpi-ranks=%d owner-blocks=%d lower=%.17e target=%.17e upper=%.17e max-split-particles=%ld keep-checks=%d split-checks=%d roulette-checks=%d capacity-rollback-checks=%d restart-checks=%d direct-source-checks=%d empty-checks=%d roulette-candidates=%lld roulette-survivors=%lld\n",
        globalResult!=0 ? "PASS" : "FAIL",PIC::nTotalThreads,
        globalBlockCount,policy.LowerCorrection,policy.TargetCorrection,
        policy.UpperCorrection,MaximumIndividualCorrectionWindowSplits,
        keepChecks,splitChecks,rouletteChecks,capacityRollbackChecks,
        restartChecks,directSourceChecks,emptyPass ? 1 : 0,
        globalRoulette[0],globalRoulette[1]);
    std::fflush(PIC::DiagnospticMessageStream);
  }
  return globalResult!=0;
#endif
}

bool Comet::RadialSimulation::ValidateParticleWeightAdaptation() {
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  if (PIC::ThisThread==0)
    std::fprintf(stderr,
        "particle-weight adaptation validation requires local weights\n");
  return false;
#else
  if (!ParticleWeightAdaptationEnabled) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "particle-weight adaptation validation requires the input switch ON\n");
    return false;
  }

  // One owner block per rank is sufficient to test simultaneous independent
  // adaptations and the following owner-to-boundary-copy synchronization.
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
      SelectValidationOwnerBlock();
  const int localHasBlock=(node!=NULL && node->block!=NULL) ? 1 : 0;
  int globalBlockCount=0;
  MPI_Allreduce(&localHasBlock,&globalBlockCount,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  if (globalBlockCount==0) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "particle-weight adaptation validation found no owner block\n");
    return false;
  }

  int species=-1;
  for (int candidate=0;candidate<PIC::nTotalSpecies;candidate++)
    if (!ParticleWeightLocked[candidate]) {
      species=candidate;
      break;
    }
  if (species<0) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "particle-weight adaptation validation requires one unlocked species\n");
    return false;
  }
  const double originalCarrier=localHasBlock ?
      GetCarrierRate(node->block,species) : 0.0;
  double globalOriginalCarrier=0.0;
  MPI_Allreduce(&originalCarrier,&globalOriginalCarrier,1,MPI_DOUBLE,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  bool localPass=true;
  int localGhostChecks=0;
  int localFaceCarrierPropagationChecks=0;
  int localComponentPopulationChecks=0;
  int localEmptyBlockChecks=0;
  int localTwoEndedChecks=0;
  int localLockedSourceChecks=0;
  int localSplitConservationChecks=0;
  int localRouletteConservationChecks=0;
  int localCapacityRollbackChecks=0;
  int localPartialMutations=0;
  int localSplitCorrectionChecks=0;
  int localRouletteCorrectionChecks=0;
  int localEmergencyCorrectionChecks=0;
  double localMaximumCorrectionError=0.0;
  int localSamplingCorrectionChecks=0;
  int localCollisionCorrectionChecks=0;
  int localBoundaryCorrectionChecks=0;
  int localRestartCorrectionChecks=0;
  int localMigrationCorrectionChecks=0;
  unsigned long long localCollisionProbabilityChecks=0;
  double localMinimumCollisionProbability=
      std::numeric_limits<double>::infinity();
  double localMaximumCollisionProbability=
      -std::numeric_limits<double>::infinity();
  double stage2MaximumFaceMultiplicity=0.0;
  long long stage2MaximumReactionAllocation=0;

  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> calculationNodes;
  CollectCalculationLeafNodes(PIC::Mesh::mesh->rootTree,calculationNodes);
  if (calculationNodes.size()>
      static_cast<std::size_t>(LLONG_MAX/BlockCellCount) ||
      globalBlockCount!=PIC::nTotalThreads) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "Stage-2 validation requires one selected owner per rank and a representable cell count\n");
    localPass=false;
  }
  const long long totalCalculationCells=
      static_cast<long long>(calculationNodes.size())*BlockCellCount;

  // Use a compact validation range.  Its 5:20 interval has exact dyadic test
  // points: a component mean of one must split three times to eight, and a
  // mean of forty must roulette once to twenty.  Only the selected species is
  // changed and the input policy is restored before the chemistry/emergency
  // regression tests below.
  const long int configuredSpeciesMinimum=
      SpeciesMinimumParticlesPerCell[species];
  const long int configuredSpeciesMaximum=
      SpeciesMaximumParticlesPerCell[species];
  SpeciesMinimumParticlesPerCell[species]=5;
  SpeciesMaximumParticlesPerCell[species]=20;
  const bool configuredChemistryEnabled=Comet::Chemistry::Enabled;
  Comet::Chemistry::Enabled=false;

  int lockedSpecies=-1;
  for (int candidate=0;candidate<PIC::nTotalSpecies;candidate++)
    if (ParticleWeightLocked[candidate]) {
      lockedSpecies=candidate;
      break;
    }
  std::vector<std::uint64_t> lockedWeightBits;
  if (lockedSpecies>=0) {
    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *testNode:calculationNodes)
      if (testNode->Thread==PIC::ThisThread && testNode->block!=NULL)
        lockedWeightBits.push_back(DoubleBits(
            testNode->block->GetLocalParticleWeight(lockedSpecies)));
  }
  else localPass=false;

  auto distributedCount = [&](long long total) -> long int {
    if (globalBlockCount<=0) return 0;
    const long long base=total/globalBlockCount;
    const long long remainder=total%globalBlockCount;
    const long long result=base+(PIC::ThisThread<remainder ? 1 : 0);
    if (result<0 || result>LONG_MAX) {
      std::fprintf(stderr,"Stage-2 validation population exceeds long int\n");
      std::abort();
    }
    return localHasBlock ? static_cast<long int>(result) : 0;
  };

  // Component split: particles occupy one owner block per rank while every
  // other cell contributes zero to the component mean.  The plan must lower
  // W/dt on the complete face component, clone only existing particles, and
  // leave every initially empty cell empty.
  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  SynchronizeAdaptedBlockWeights();
  const long int splitInitial=distributedCount(totalCalculationCells);
  const double splitInitialWeight=localHasBlock ?
      node->block->GetLocalParticleWeight(species) : 1.0;
  if (localHasBlock) SeedValidationParticles(node,species,splitInitial);
  const cValidationMoments splitBefore=CollectValidationMoments(species);
  AdaptParticleWeights(0);
  const cValidationMoments splitAfter=CollectValidationMoments(species);
  const cValidationCorrectionAudit splitCorrectionAudit=
      AuditValidationCorrections(species);
  localMaximumCorrectionError=std::max(localMaximumCorrectionError,
      splitCorrectionAudit.MaximumRelativeError);
  const double splitExpectedWeight=std::ldexp(splitInitialWeight,-3);
  bool splitLocalPass=true;
  if (localHasBlock) {
    splitLocalPass=
        CountSpeciesInBlock(node->block,species)==8*splitInitial &&
        DoubleBits(node->block->GetLocalParticleWeight(species))==
            DoubleBits(splitExpectedWeight) &&
        OtherValidationCellsAreEmpty(node->block);
  }
  const bool splitMomentPass=ValidationMomentsEqual(splitBefore,splitAfter,
      5.0e-13);
  cComponentCarrierTransaction splitState=
      BuildComponentCarrierTransaction();
  bool componentPopulationPass=true,hasPopulatedComponent=false;
  bool emptyBlocksPass=true;
  long long localEmptyBlocks=0;
  for (std::size_t component=0;
       component<splitState.ComponentNodes.size();component++) {
    const std::size_t componentOffset=
        component*PIC::nTotalSpecies+species;
    const long long count=
        splitState.ComponentParticleCount[componentOffset];
    if (count==0) continue;
    hasPopulatedComponent=true;
    const long double mean=static_cast<long double>(count)/
        splitState.ComponentCellCount[component];
    componentPopulationPass=componentPopulationPass &&
        mean>=5.0L && mean<=20.0L;
  }
  componentPopulationPass=componentPopulationPass && hasPopulatedComponent;
  for (std::size_t block=0;block<splitState.Nodes.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *testNode=splitState.Nodes[block];
    if (testNode->Thread!=PIC::ThisThread || testNode->block==NULL) continue;
    const std::size_t offset=block*PIC::nTotalSpecies+species;
    if (splitState.GlobalCounts[offset]!=0) continue;
    localEmptyBlocks++;
    emptyBlocksPass=emptyBlocksPass &&
        CountSpeciesInBlock(testNode->block,species)==0 &&
        std::fabs(GetCarrierRate(testNode->block,species)-
            std::ldexp(globalOriginalCarrier,-3))<=
            512.0*std::numeric_limits<double>::epsilon()*
                std::ldexp(globalOriginalCarrier,-3);
  }
  long long globalEmptyBlocks=0;
  MPI_Allreduce(&localEmptyBlocks,&globalEmptyBlocks,1,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  const bool componentSplitPass=splitLocalPass && splitMomentPass &&
      splitCorrectionAudit.InvalidCount==0 &&
      componentPopulationPass && emptyBlocksPass && globalEmptyBlocks>0 &&
      splitState.PropagationStatistics.MaximumCloneMeanAfter<=
          1.0+256.0*std::numeric_limits<double>::epsilon();
  localPass=localPass && componentSplitPass;
  if (componentSplitPass && localHasBlock) {
    localComponentPopulationChecks++;
    localEmptyBlockChecks++;
    localSplitConservationChecks++;
    localSplitCorrectionChecks++;
    localFaceCarrierPropagationChecks++;
  }
  stage2MaximumFaceMultiplicity=std::max(stage2MaximumFaceMultiplicity,
      splitState.PropagationStatistics.MaximumCloneMeanAfter);
  localPass=ValidateSynchronizedWeightCopies(node,species,
      splitExpectedWeight,localGhostChecks) && localPass;

  // A second exact split places the population only at the two ends of the
  // radial chain.  It closes the regression in which an empty interior block
  // used to raise the carrier again after local splitting.
  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  SynchronizeAdaptedBlockWeights();
  bool twoEndedPass=calculationNodes.size()>=2;
  if (twoEndedPass) {
    const long long firstCount=totalCalculationCells/2;
    const long long lastCount=totalCalculationCells-firstCount;
    if (calculationNodes.front()->Thread==PIC::ThisThread)
      SeedValidationParticles(calculationNodes.front(),species,
          static_cast<long int>(firstCount));
    if (calculationNodes.back()->Thread==PIC::ThisThread)
      SeedValidationParticles(calculationNodes.back(),species,
          static_cast<long int>(lastCount));
    const cValidationMoments before=CollectValidationMoments(species);
    AdaptParticleWeights(0);
    const cValidationMoments after=CollectValidationMoments(species);
    twoEndedPass=ValidationMomentsEqual(before,after,5.0e-13) &&
        LastFaceCarrierPropagationStatistics.MaximumCloneMeanAfter<=
            1.0+256.0*std::numeric_limits<double>::epsilon();
    if (calculationNodes.front()->Thread==PIC::ThisThread)
      twoEndedPass=twoEndedPass &&
          CountSpeciesInBlock(calculationNodes.front()->block,species)==
              8*firstCount;
    if (calculationNodes.back()->Thread==PIC::ThisThread)
      twoEndedPass=twoEndedPass &&
          CountSpeciesInBlock(calculationNodes.back()->block,species)==
              8*lastCount;
  }
  int twoEndedLocal=twoEndedPass ? 1 : 0,twoEndedGlobal=0;
  MPI_Allreduce(&twoEndedLocal,&twoEndedGlobal,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  localPass=localPass && twoEndedGlobal!=0;
  if (twoEndedGlobal!=0 && localHasBlock) localTwoEndedChecks++;

  // Component roulette: start from an exact mean of forty.  A single
  // doubling gives an expected mean of twenty.  The represented moments are
  // tested against the exact weighted Bernoulli variance, so differing local
  // timesteps/weights cannot make the statistical bound rank dependent.
  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  SynchronizeAdaptedBlockWeights();
  const long int rouletteInitial=distributedCount(40*totalCalculationCells);
  const double rouletteInitialWeight=localHasBlock ?
      node->block->GetLocalParticleWeight(species) : 1.0;
  if (localHasBlock) SeedValidationParticles(node,species,rouletteInitial);
  const cValidationMoments rouletteBefore=CollectValidationMoments(species);
  rnd_seed(99173+PIC::ThisThread);
  AdaptParticleWeights(0);
  const cValidationMoments rouletteAfter=CollectValidationMoments(species);
  const cValidationCorrectionAudit rouletteCorrectionAudit=
      AuditValidationCorrections(species);
  localMaximumCorrectionError=std::max(localMaximumCorrectionError,
      rouletteCorrectionAudit.MaximumRelativeError);
  // Survival probability is 1/2 but every survivor receives twice the old
  // weight.  Represented moments are therefore unbiased at their original
  // values.  For unequal old weights, Var[sum(2*w_i*Bernoulli(1/2))]
  // is exactly sum(w_i^2), which was accumulated before mutation.
  const cValidationMoments rouletteExpected=rouletteBefore;
  const double rouletteSigma=std::sqrt(rouletteBefore.WeightSquared);
  const double rouletteTolerance=rouletteExpected.Number>0.0 ?
      7.0*rouletteSigma/rouletteExpected.Number+5.0e-13 : 0.0;
  const bool rouletteMomentPass=
      std::fabs(rouletteAfter.Number-rouletteExpected.Number)<=
          7.0*rouletteSigma &&
      ValidationMomentsEqual(rouletteAfter,rouletteExpected,
          rouletteTolerance);
  bool rouletteLocalPass=true;
  if (localHasBlock) {
    const long int actualCount=CountSpeciesInBlock(node->block,species);
    const double expectation=0.5*rouletteInitial;
    const double sigma=std::sqrt(0.25*rouletteInitial);
    rouletteLocalPass=std::fabs(actualCount-expectation)<=7.0*sigma &&
        DoubleBits(node->block->GetLocalParticleWeight(species))==
            DoubleBits(2.0*rouletteInitialWeight) &&
        OtherValidationCellsAreEmpty(node->block);
  }
  localPass=localPass && rouletteMomentPass && rouletteLocalPass &&
      rouletteCorrectionAudit.InvalidCount==0;
  if (rouletteMomentPass && rouletteLocalPass && localHasBlock)
    localRouletteConservationChecks++,localRouletteCorrectionChecks++;
  localPass=ValidateSynchronizedWeightCopies(node,species,
      2.0*rouletteInitialWeight,localGhostChecks) && localPass;

  // Capacity rejection is a pure planning operation.  Force zero additional
  // capacity while a three-step split is needed, then compare exact particle
  // count, represented moments, and every local owner carrier bit before and
  // after the rejected plan.  No mutation routine is called on failure.
  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  SynchronizeAdaptedBlockWeights();
  if (localHasBlock) SeedValidationParticles(node,species,splitInitial);
  const long long capacityParticlesBefore=PIC::ParticleBuffer::NAllPart;
  const cValidationMoments capacityMomentsBefore=
      CollectValidationMoments(species);
  std::vector<std::uint64_t> capacityWeightBits;
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *testNode:calculationNodes)
    if (testNode->Thread==PIC::ThisThread && testNode->block!=NULL)
      capacityWeightBits.push_back(DoubleBits(
          testNode->block->GetLocalParticleWeight(species)));
  const cComponentCarrierTransaction rejectedTransaction=
      BuildComponentCarrierTransaction(0);
  const cValidationMoments capacityMomentsAfter=
      CollectValidationMoments(species);
  bool capacityRollbackPass=!rejectedTransaction.CapacityAccepted &&
      PIC::ParticleBuffer::NAllPart==capacityParticlesBefore &&
      ValidationMomentsEqual(capacityMomentsBefore,capacityMomentsAfter,0.0);
  std::size_t capacityBitIndex=0;
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *testNode:calculationNodes)
    if (testNode->Thread==PIC::ThisThread && testNode->block!=NULL) {
      capacityRollbackPass=capacityRollbackPass &&
          DoubleBits(testNode->block->GetLocalParticleWeight(species))==
              capacityWeightBits[capacityBitIndex++];
    }
  int capacityLocal=capacityRollbackPass ? 1 : 0,capacityGlobal=0;
  MPI_Allreduce(&capacityLocal,&capacityGlobal,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (capacityGlobal==0) localPartialMutations++;
  localPass=localPass && capacityGlobal!=0;
  if (capacityGlobal!=0 && localHasBlock) localCapacityRollbackChecks++;

  // Source locks cover every ordinary transaction above.  Compare exact
  // IEEE-754 bits on every owner, independently of the internal assertion in
  // AdaptParticleWeights, so the validation record proves the contract.
  bool lockedSourcePass=lockedSpecies>=0;
  std::size_t lockedBitIndex=0;
  if (lockedSpecies>=0)
    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *testNode:calculationNodes)
      if (testNode->Thread==PIC::ThisThread && testNode->block!=NULL) {
        lockedSourcePass=lockedSourcePass &&
            DoubleBits(testNode->block->GetLocalParticleWeight(
                lockedSpecies))==lockedWeightBits[lockedBitIndex++];
      }
  int lockedLocal=lockedSourcePass ? 1 : 0,lockedGlobal=0;
  MPI_Allreduce(&lockedLocal,&lockedGlobal,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  localPass=localPass && lockedGlobal!=0;
  if (lockedGlobal!=0 && localHasBlock) localLockedSourceChecks++;

  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  SynchronizeAdaptedBlockWeights();
  SpeciesMinimumParticlesPerCell[species]=configuredSpeciesMinimum;
  SpeciesMaximumParticlesPerCell[species]=configuredSpeciesMaximum;
  Comet::Chemistry::Enabled=configuredChemistryEnabled;

  int localPhotochemicalChecks=0;
  if (Comet::Chemistry::Enabled) {
    const ChemicalSolver::cReactionNetwork &network=
        Comet::Chemistry::GetNetwork();
    std::size_t validationReaction=network.GetReactionCount();
    for (std::size_t reactionIndex=0;
         reactionIndex<network.GetReactionCount();reactionIndex++) {
      if (!network.IsUnitUnimolecular(reactionIndex) ||
          !(network.GetRateCoefficient(reactionIndex,0.0,
              Comet::Chemistry::HeliocentricDistanceAU)>0.0)) continue;
      const int parent=
          network.GetReactantTerm(reactionIndex,0).Species;
      const ChemicalSolver::cReaction &reaction=
          network.GetReaction(reactionIndex);
      bool hasDistinctProduct=false,allDistinctProductsUnlocked=true;
      for (std::size_t term=0;term<reaction.ProductCount;term++) {
        const int product=
            network.GetProductTerm(reactionIndex,term).Species;
        if (product==parent) continue;
        hasDistinctProduct=true;
        if (ParticleWeightLocked[product])
          allDistinctProductsUnlocked=false;
      }
      if (hasDistinctProduct && allDistinctProductsUnlocked) {
        validationReaction=reactionIndex;
        break;
      }
    }

    if (validationReaction==network.GetReactionCount()) {
      if (PIC::ThisThread==0)
        std::fprintf(stderr,
            "particle-weight adaptation validation found no photochemical daughter channel\n");
      localPass=false;
    }
    else {
      const int parent=
          network.GetReactantTerm(validationReaction,0).Species;
      const ChemicalSolver::cReaction &reaction=
          network.GetReaction(validationReaction);
      std::array<double,PIC::nTotalSpecies> localSavedCarriers{};
      std::array<double,PIC::nTotalSpecies> savedCarriers{};
      std::array<bool,PIC::nTotalSpecies> testedDaughter{};

      PIC::ParticleBuffer::DeleteAllParticles();
      if (localHasBlock) {
        for (int testSpecies=0;testSpecies<PIC::nTotalSpecies;testSpecies++)
          localSavedCarriers[testSpecies]=
              GetCarrierRate(node->block,testSpecies);

        // Keep the parent exactly at the lower population bound so ordinary
        // population adaptation leaves its weight unchanged.  Every distinct
        // product begins absent from all cells in the selected block.
        const long int parentCount=
            SpeciesMinimumParticlesPerCell[parent]*BlockCellCount;
        SeedValidationParticles(node,parent,parentCount);
      }
      MPI_Allreduce(localSavedCarriers.data(),savedCarriers.data(),
          PIC::nTotalSpecies,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

      double absentCoefficient=0.0;
      for (std::size_t term=0;term<reaction.ProductCount;term++) {
        const ChemicalSolver::cStoichiometricTerm &product=
            network.GetProductTerm(validationReaction,term);
        if (product.Species!=parent) {
          absentCoefficient+=product.Coefficient;
          testedDaughter[product.Species]=true;
        }
      }

      // Correction-off chemistry must raise an unsafe residual daughter
      // carrier before allocation.  Correction-on chemistry always creates
      // one daughter per concrete stoichiometric product, so its hard count
      // is carrier-independent and this test instead proves that adaptation
      // leaves the carrier alone while retaining the irreducible bound.
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      const double unsafeCarrier=savedCarriers[parent]*absentCoefficient/
          (8.0*MaximumPhotochemicalProductsPerReaction);
      for (int product=0;product<PIC::nTotalSpecies;product++)
        if (testedDaughter[product])
          SetWeightOverTimeStep(product,unsafeCarrier,
              PIC::Mesh::mesh->rootTree);
      SynchronizeAdaptedBlockWeights();
#else
      (void)absentCoefficient;
#endif

      long long multiplicityBefore=0;
      if (localHasBlock) {
        const std::array<long int,PIC::nTotalSpecies> counts=
            CountAllSpeciesInBlock(node->block);
        multiplicityBefore=
            MaximumReactionDaughterAllocation(node->block,counts);
      }
      AdaptParticleWeights(0);
      if (localHasBlock) {
        const std::array<long int,PIC::nTotalSpecies> counts=
            CountAllSpeciesInBlock(node->block);
        const long long multiplicityAfter=
            MaximumReactionDaughterAllocation(node->block,counts);
        stage2MaximumReactionAllocation=std::max(
            stage2MaximumReactionAllocation,multiplicityAfter);
        bool daughtersRemainEmpty=true;
        for (int product=0;product<PIC::nTotalSpecies;product++) {
          if (!testedDaughter[product]) continue;
          daughtersRemainEmpty=daughtersRemainEmpty && counts[product]==0;
          localPass=ValidateSynchronizedWeightCopies(node,product,
              node->block->GetLocalParticleWeight(product),localGhostChecks) &&
              localPass;
        }
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        localPass=localPass && multiplicityBefore==multiplicityAfter &&
            multiplicityAfter<=MaximumPhotochemicalProductsPerReaction &&
            daughtersRemainEmpty;
#else
        localPass=localPass &&
            multiplicityBefore>MaximumPhotochemicalProductsPerReaction &&
            multiplicityAfter<=MaximumPhotochemicalProductsPerReaction &&
            daughtersRemainEmpty;
#endif
        localPhotochemicalChecks++;
      }
      for (int testSpecies=0;testSpecies<PIC::nTotalSpecies;testSpecies++)
        SetWeightOverTimeStep(testSpecies,savedCarriers[testSpecies],
            PIC::Mesh::mesh->rootTree);
      SynchronizeAdaptedBlockWeights();
    }
  }

  // Exercise the actual core hook, not only the roulette helper.  A modest
  // synthetic population is enough: temporarily set the threshold to 75% of
  // that population's measured maximum-rank utilization, then invoke the
  // regular particle-exchange API.  The callback must run before exchange,
  // double every block/species weight by a common power of two, and leave the
  // maximum rank below the temporary trigger.  The production input value is
  // restored immediately afterward.
  int localEmergencyChecks=0;
  PIC::ParticleBuffer::DeleteAllParticles();
  const long int emergencyInitialCount=4096;
  const double emergencyInitialWeight=localHasBlock ?
      node->block->GetLocalParticleWeight(species) : 1.0;
  if (localHasBlock)
    SeedValidationParticles(node,species,emergencyInitialCount);

  const long long localEmergencyInitial=PIC::ParticleBuffer::NAllPart;
  long long maximumEmergencyInitial=0;
  MPI_Allreduce(&localEmergencyInitial,&maximumEmergencyInitial,1,
      MPI_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  const double configuredEmergencyThreshold=
      EmergencyParticleBufferUtilization;
  const double validationEmergencyThreshold=
      0.75*static_cast<double>(maximumEmergencyInitial)/
      PIC::ParticleBuffer::MaxNPart;
  if (!(validationEmergencyThreshold>0.0) ||
      !(validationEmergencyThreshold<1.0)) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "cannot construct emergency particle-weight validation threshold\n");
    localPass=false;
  }
  else {
    EmergencyParticleBufferUtilization=validationEmergencyThreshold;
    rnd_seed(77123+PIC::ThisThread);
    PIC::Parallel::ExchangeParticleData();
    EmergencyParticleBufferUtilization=configuredEmergencyThreshold;
    const cValidationCorrectionAudit emergencyCorrectionAudit=
        AuditValidationCorrections(species);
    localMaximumCorrectionError=std::max(localMaximumCorrectionError,
        emergencyCorrectionAudit.MaximumRelativeError);
    if (emergencyCorrectionAudit.InvalidCount!=0) localPass=false;
    else if (localHasBlock) localEmergencyCorrectionChecks++;

    const long long localEmergencyAfter=PIC::ParticleBuffer::NAllPart;
    long long maximumEmergencyAfter=0;
    MPI_Allreduce(&localEmergencyAfter,&maximumEmergencyAfter,1,
        MPI_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    const double maximumUtilizationAfter=
        static_cast<double>(maximumEmergencyAfter)/
        PIC::ParticleBuffer::MaxNPart;

    if (localHasBlock) {
      const long int actualCount=CountSpeciesInBlock(node->block,species);
      const double actualWeight=
          node->block->GetLocalParticleWeight(species);
      double weightRatio=actualWeight/emergencyInitialWeight;
      int weightPower=0;
      while (weightRatio>1.5 && weightPower<1024) {
        weightRatio*=0.5;
        weightPower++;
      }
      const double survivalProbability=std::ldexp(1.0,-weightPower);
      const double expectation=emergencyInitialCount*survivalProbability;
      const double sigma=std::sqrt(emergencyInitialCount*
          survivalProbability*(1.0-survivalProbability));
      const bool countPass=std::fabs(actualCount-expectation)<=6.0*sigma;
      const bool weightPass=weightPower>0 &&
          std::fabs(weightRatio-1.0)<=1.0e-13;
      const bool capacityPass=
          maximumUtilizationAfter<validationEmergencyThreshold;
      localPass=localPass && countPass && weightPass && capacityPass &&
          OtherValidationCellsAreEmpty(node->block);
      localPass=ValidateSynchronizedWeightCopies(node,species,actualWeight,
          localGhostChecks) && localPass;
      localEmergencyChecks++;
    }
  }

  // Exercise the transactional allocation-pressure path independently of
  // the ordinary high-water callback.  A production 50% target would require
  // more than ten million validation particles in the configured 20-million
  // slot buffer, so scale the target from this compact synthetic population.
  // The control law is identical: the collective request must continue
  // factor-two roulette until the measured maximum rank is at or below the
  // temporary low-water mark, consume its one-shot flag, and keep owner/ghost
  // weights identical.
  int localAllocationBurstChecks=0;
  PIC::ParticleBuffer::DeleteAllParticles();
  const long int allocationInitialCount=4096;
  const double allocationInitialWeight=localHasBlock ?
      node->block->GetLocalParticleWeight(species) : 1.0;
  if (localHasBlock)
    SeedValidationParticles(node,species,allocationInitialCount);
  const long long localAllocationInitial=PIC::ParticleBuffer::NAllPart;
  long long maximumAllocationInitial=0;
  MPI_Allreduce(&localAllocationInitial,&maximumAllocationInitial,1,
      MPI_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  const double configuredAllocationLowWatermark=AllocationLowWatermark;
  const double validationAllocationLowWatermark=
      0.30*static_cast<double>(maximumAllocationInitial)/
      PIC::ParticleBuffer::MaxNPart;
  if (!(validationAllocationLowWatermark>0.0) ||
      !(validationAllocationLowWatermark<EmergencyParticleBufferUtilization)) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "cannot construct allocation-burst low-water validation target\n");
    localPass=false;
  }
  else {
    AllocationLowWatermark=validationAllocationLowWatermark;
    rnd_seed(83117+PIC::ThisThread);
    RequestAllocationBurstEmergency();
    AllocationLowWatermark=configuredAllocationLowWatermark;
    const cValidationCorrectionAudit allocationCorrectionAudit=
        AuditValidationCorrections(species);
    localMaximumCorrectionError=std::max(localMaximumCorrectionError,
        allocationCorrectionAudit.MaximumRelativeError);
    if (allocationCorrectionAudit.InvalidCount!=0) localPass=false;
    else if (localHasBlock) localEmergencyCorrectionChecks++;

    const long long localAllocationAfter=PIC::ParticleBuffer::NAllPart;
    long long maximumAllocationAfter=0;
    MPI_Allreduce(&localAllocationAfter,&maximumAllocationAfter,1,
        MPI_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    const double maximumAllocationUtilizationAfter=
        static_cast<double>(maximumAllocationAfter)/
        PIC::ParticleBuffer::MaxNPart;
    if (localHasBlock) {
      const long int actualCount=CountSpeciesInBlock(node->block,species);
      const double actualWeight=
          node->block->GetLocalParticleWeight(species);
      double weightRatio=actualWeight/allocationInitialWeight;
      int weightPower=0;
      while (weightRatio>1.5 && weightPower<1024) {
        weightRatio*=0.5;
        weightPower++;
      }
      const double survivalProbability=std::ldexp(1.0,-weightPower);
      const double expectation=allocationInitialCount*survivalProbability;
      const double sigma=std::sqrt(allocationInitialCount*
          survivalProbability*(1.0-survivalProbability));
      const bool countPass=std::fabs(actualCount-expectation)<=6.0*sigma;
      const bool weightPass=weightPower>=2 &&
          std::fabs(weightRatio-1.0)<=1.0e-13;
      const bool lowWaterPass=maximumAllocationUtilizationAfter<=
          validationAllocationLowWatermark;
      localPass=localPass && countPass && weightPass && lowWaterPass &&
          !AllocationBurstEmergencyPending &&
          OtherValidationCellsAreEmpty(node->block);
      localPass=ValidateSynchronizedWeightCopies(node,species,actualWeight,
          localGhostChecks) && localPass;
      localAllocationBurstChecks++;
    }
  }
  // Restore the configured value on every error path as well as the normal
  // path before continuing to the unrelated sampling-population test.
  AllocationLowWatermark=configuredAllocationLowWatermark;

  // Sampling-population trigger test: temporarily use a compact but valid
  // per-species range so the exact production rule can be crossed without
  // allocating 10*nSpecies*200 particles per cell on every validation rank.
  // AdaptParticleWeights must record the pre-adaptation block-mean excess but
  // must not execute the global emergency there.  The subsequent ordinary
  // particle exchange must consume the request in the existing callback and
  // add exactly one or more common factor-two roulette steps.
  int localSamplingPopulationEmergencyChecks=0;
  PIC::ParticleBuffer::DeleteAllParticles();
  const long int configuredMinimumParticlesPerCell=MinimumParticlesPerCell;
  const long int configuredMaximumParticlesPerCell=MaximumParticlesPerCell;
  std::array<long int,PIC::nTotalSpecies> configuredMinimumBySpecies{};
  std::array<long int,PIC::nTotalSpecies> configuredMaximumBySpecies{};
  for (int testSpecies=0;testSpecies<PIC::nTotalSpecies;testSpecies++) {
    configuredMinimumBySpecies[testSpecies]=
        SpeciesMinimumParticlesPerCell[testSpecies];
    configuredMaximumBySpecies[testSpecies]=
        SpeciesMaximumParticlesPerCell[testSpecies];
  }
  constexpr long int validationMinimumParticlesPerCell=5;
  constexpr long int validationMaximumParticlesPerCell=20;
  MinimumParticlesPerCell=validationMinimumParticlesPerCell;
  MaximumParticlesPerCell=validationMaximumParticlesPerCell;
  for (int testSpecies=0;testSpecies<PIC::nTotalSpecies;testSpecies++) {
    SpeciesMinimumParticlesPerCell[testSpecies]=
        validationMinimumParticlesPerCell;
    SpeciesMaximumParticlesPerCell[testSpecies]=
        validationMaximumParticlesPerCell;
  }
  const long long validationSamplingBlockThreshold=
      GetSamplingPopulationEmergencyBlockThreshold();
  if (validationSamplingBlockThreshold>LONG_MAX/2) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "cannot construct sampling-population emergency validation count\n");
    localPass=false;
  }
  else {
    const long int samplingInitialCount=
        static_cast<long int>(2*validationSamplingBlockThreshold);
    const double samplingInitialWeight=localHasBlock ?
        node->block->GetLocalParticleWeight(species) : 1.0;
    if (localHasBlock)
      SeedValidationParticles(node,species,samplingInitialCount);

    rnd_seed(88129+PIC::ThisThread);
    AdaptParticleWeights(0);
    const bool requestRecorded=SamplingPopulationEmergencyPending &&
        SamplingPopulationEmergencyMaximumBlockParticles>
            SamplingPopulationEmergencyThresholdBlockParticles &&
        SamplingPopulationEmergencyThresholdBlockParticles==
            validationSamplingBlockThreshold;

    // Restore the production input range before exercising the callback.  Its
    // pending record contains the threshold that actually caused the request,
    // so diagnostics remain traceable to the sampling decision.
    MinimumParticlesPerCell=configuredMinimumParticlesPerCell;
    MaximumParticlesPerCell=configuredMaximumParticlesPerCell;
    for (int testSpecies=0;testSpecies<PIC::nTotalSpecies;testSpecies++) {
      SpeciesMinimumParticlesPerCell[testSpecies]=
          configuredMinimumBySpecies[testSpecies];
      SpeciesMaximumParticlesPerCell[testSpecies]=
          configuredMaximumBySpecies[testSpecies];
    }
    rnd_seed(99241+PIC::ThisThread);
    PIC::Parallel::ExchangeParticleData();
    const bool requestConsumed=!SamplingPopulationEmergencyPending;
    const cValidationCorrectionAudit samplingCorrectionAudit=
        AuditValidationCorrections(species);
    localMaximumCorrectionError=std::max(localMaximumCorrectionError,
        samplingCorrectionAudit.MaximumRelativeError);
    if (samplingCorrectionAudit.InvalidCount!=0) localPass=false;
    else if (localHasBlock) localEmergencyCorrectionChecks++;

    if (localHasBlock) {
      const long int actualCount=CountSpeciesInBlock(node->block,species);
      const double actualWeight=
          node->block->GetLocalParticleWeight(species);
      double weightRatio=actualWeight/samplingInitialWeight;
      int weightPower=0;
      while (weightRatio>1.5 && weightPower<1024) {
        weightRatio*=0.5;
        weightPower++;
      }
      const double survivalProbability=std::ldexp(1.0,-weightPower);
      const double expectation=samplingInitialCount*survivalProbability;
      const double sigma=std::sqrt(samplingInitialCount*
          survivalProbability*(1.0-survivalProbability));
      const bool countPass=std::fabs(actualCount-expectation)<=6.0*sigma;
      const bool weightPass=weightPower>0 &&
          std::fabs(weightRatio-1.0)<=1.0e-13;
      localPass=localPass && requestRecorded && requestConsumed &&
          countPass && weightPass && OtherValidationCellsAreEmpty(node->block);
      localPass=ValidateSynchronizedWeightCopies(node,species,actualWeight,
          localGhostChecks) && localPass;
      localSamplingPopulationEmergencyChecks++;
    }
  }

  // The error branch above does not enter the callback, so restore the input
  // range unconditionally before any final cleanup or result reporting.
  MinimumParticlesPerCell=configuredMinimumParticlesPerCell;
  MaximumParticlesPerCell=configuredMaximumParticlesPerCell;
  for (int testSpecies=0;testSpecies<PIC::nTotalSpecies;testSpecies++) {
    SpeciesMinimumParticlesPerCell[testSpecies]=
        configuredMinimumBySpecies[testSpecies];
    SpeciesMaximumParticlesPerCell[testSpecies]=
        configuredMaximumBySpecies[testSpecies];
  }

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  // Stage-4 paths that are independent of the population controller itself:
  // the production sampler, general NTC kernel, radial boundary callback, and
  // full-record restart.  Keeping them in this initialized focused driver
  // reuses the exact mesh, local W/dt values, MPI decomposition, and particle
  // layout of the application rather than a reduced mock.
  localPass=ValidateIndividualCorrectionSampling(node,
      localSamplingCorrectionChecks) && localPass;
  localPass=ValidateIndividualCorrectionNTC(node,
      localCollisionCorrectionChecks,localCollisionProbabilityChecks,
      localMinimumCollisionProbability,localMaximumCollisionProbability) &&
      localPass;
  localPass=ValidateIndividualCorrectionBoundaryAndRestart(node,species,
      localBoundaryCorrectionChecks,localRestartCorrectionChecks) &&
      localPass;
  // Run ownership migration last because it deliberately changes the AMR
  // distribution.  All following cleanup uses the rebuilt owner table and no
  // longer relies on the pre-migration validation-node pointer.
  localPass=ValidateIndividualCorrectionMigration(species,
      localMigrationCorrectionChecks) && localPass;
#endif

  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  SynchronizeAdaptedBlockWeights();

  const int localResult=localPass ? 1 : 0;
  int globalResult=0;
  MPI_Allreduce(&localResult,&globalResult,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  int globalGhostChecks=0;
  MPI_Allreduce(&localGhostChecks,&globalGhostChecks,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  int globalPhotochemicalChecks=0;
  MPI_Allreduce(&localPhotochemicalChecks,&globalPhotochemicalChecks,1,
      MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  int globalEmergencyChecks=0;
  MPI_Allreduce(&localEmergencyChecks,&globalEmergencyChecks,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  int globalSamplingPopulationEmergencyChecks=0;
  MPI_Allreduce(&localSamplingPopulationEmergencyChecks,
      &globalSamplingPopulationEmergencyChecks,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  int globalAllocationBurstChecks=0;
  MPI_Allreduce(&localAllocationBurstChecks,&globalAllocationBurstChecks,1,
      MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  int globalFaceCarrierPropagationChecks=0;
  MPI_Allreduce(&localFaceCarrierPropagationChecks,
      &globalFaceCarrierPropagationChecks,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  int localStage2Checks[8]={localComponentPopulationChecks,
      localEmptyBlockChecks,localTwoEndedChecks,localLockedSourceChecks,
      localSplitConservationChecks,localRouletteConservationChecks,
      localCapacityRollbackChecks,localPartialMutations};
  int globalStage2Checks[8]={0,0,0,0,0,0,0,0};
  MPI_Allreduce(localStage2Checks,globalStage2Checks,8,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  double globalStage2MaximumFaceMultiplicity=0.0;
  MPI_Allreduce(&stage2MaximumFaceMultiplicity,
      &globalStage2MaximumFaceMultiplicity,1,MPI_DOUBLE,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  long long globalStage2MaximumReactionAllocation=0;
  MPI_Allreduce(&stage2MaximumReactionAllocation,
      &globalStage2MaximumReactionAllocation,1,MPI_LONG_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  int localCorrectionChecks[3]={localSplitCorrectionChecks,
      localRouletteCorrectionChecks,localEmergencyCorrectionChecks};
  int globalCorrectionChecks[3]={0,0,0};
  MPI_Allreduce(localCorrectionChecks,globalCorrectionChecks,3,MPI_INT,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  double globalMaximumCorrectionError=0.0;
  MPI_Allreduce(&localMaximumCorrectionError,
      &globalMaximumCorrectionError,1,MPI_DOUBLE,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  int localCoreCorrectionChecks[5]={localSamplingCorrectionChecks,
      localCollisionCorrectionChecks,localBoundaryCorrectionChecks,
      localRestartCorrectionChecks,localMigrationCorrectionChecks};
  int globalCoreCorrectionChecks[5]={0,0,0,0,0};
  MPI_Allreduce(localCoreCorrectionChecks,globalCoreCorrectionChecks,5,
      MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  unsigned long long globalCollisionProbabilityChecks=0;
  MPI_Allreduce(&localCollisionProbabilityChecks,
      &globalCollisionProbabilityChecks,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  double globalMinimumCollisionProbability=0.0;
  double globalMaximumCollisionProbability=0.0;
  MPI_Allreduce(&localMinimumCollisionProbability,
      &globalMinimumCollisionProbability,1,MPI_DOUBLE,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localMaximumCollisionProbability,
      &globalMaximumCollisionProbability,1,MPI_DOUBLE,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::nTotalThreads>1 && globalGhostChecks==0) globalResult=0;
  if (Comet::Chemistry::Enabled &&
      globalPhotochemicalChecks!=globalBlockCount) globalResult=0;
  if (globalEmergencyChecks!=globalBlockCount) globalResult=0;
  if (globalAllocationBurstChecks!=globalBlockCount) globalResult=0;
  if (globalSamplingPopulationEmergencyChecks!=globalBlockCount)
    globalResult=0;
  if (globalFaceCarrierPropagationChecks!=globalBlockCount)
    globalResult=0;
  for (int check=0;check<7;check++)
    if (globalStage2Checks[check]!=globalBlockCount) globalResult=0;
  if (globalStage2Checks[7]!=0 ||
      globalCorrectionChecks[0]!=globalBlockCount ||
      globalCorrectionChecks[1]!=globalBlockCount ||
      globalCorrectionChecks[2]!=3*globalBlockCount ||
      globalMaximumCorrectionError>
          256.0*std::numeric_limits<double>::epsilon() ||
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      globalCoreCorrectionChecks[0]!=globalBlockCount ||
      globalCoreCorrectionChecks[1]!=globalBlockCount ||
      globalCoreCorrectionChecks[2]!=globalBlockCount ||
      globalCoreCorrectionChecks[3]!=globalBlockCount ||
      globalCoreCorrectionChecks[4]!=globalBlockCount ||
      globalCollisionProbabilityChecks==0 ||
      !(globalMinimumCollisionProbability>=0.0) ||
      !(globalMaximumCollisionProbability<=1.0) ||
#endif
      globalStage2MaximumFaceMultiplicity>
          1.0+256.0*std::numeric_limits<double>::epsilon() ||
      globalStage2MaximumReactionAllocation>
          MaximumPhotochemicalProductsPerReaction)
    globalResult=0;
  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "%s particle-weight adaptation validation individual-correction=%s range=%ld:%ld emergency-buffer-utilization=%.6f allocation-low-watermark=%.6f allocation-reserve=%.6f max-face-transfer-multiplicity=%ld max-photochemical-products=%ld owner-blocks=%d emergency-checks=%d allocation-burst-checks=%d sampling-population-emergency-checks=%d face-carrier-propagation-checks=%d photochemical-checks=%d ghost-weight-checks=%d\n",
        globalResult!=0 ? "PASS" : "FAIL",
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        "on",
#else
        "off",
#endif
        MinimumParticlesPerCell,
        MaximumParticlesPerCell,EmergencyParticleBufferUtilization,
        AllocationLowWatermark,AllocationReserveFraction,
        MaximumFaceTransferMultiplicity,
        MaximumPhotochemicalProductsPerReaction,globalBlockCount,
        globalEmergencyChecks,globalAllocationBurstChecks,
        globalSamplingPopulationEmergencyChecks,
        globalFaceCarrierPropagationChecks,
        globalPhotochemicalChecks,globalGhostChecks);
    std::fprintf(PIC::DiagnospticMessageStream,
        "%s component-carrier validation mpi-ranks=%d owner-blocks=%d component-population-checks=%d empty-block-checks=%d two-ended-checks=%d locked-source-checks=%d reaction-allocation-checks=%d split-conservation-checks=%d roulette-conservation-checks=%d capacity-rollback-checks=%d ghost-weight-checks=%d max-face-multiplicity=%.17e max-reaction-allocation=%lld reaction-limit=%ld partial-mutations=%d\n",
        globalResult!=0 ? "PASS" : "FAIL",PIC::nTotalThreads,
        globalBlockCount,globalStage2Checks[0],globalStage2Checks[1],
        globalStage2Checks[2],globalStage2Checks[3],
        globalPhotochemicalChecks,globalStage2Checks[4],
        globalStage2Checks[5],globalStage2Checks[6],globalGhostChecks,
        globalStage2MaximumFaceMultiplicity,
        globalStage2MaximumReactionAllocation,
        MaximumPhotochemicalProductsPerReaction,globalStage2Checks[7]);
    std::fprintf(PIC::DiagnospticMessageStream,
        "%s individual-weight adaptation validation correction-mode=%s mpi-ranks=%d correction-min=%.17g correction-max=%.17g correction-values=%zu max-correction-error=%.17e partial-mutations=%d split-correction-checks=%d roulette-correction-checks=%d emergency-correction-checks=%d corrected-moment-checks=%d capacity-rollback-checks=%d ghost-weight-publication-checks=%d\n",
        globalResult!=0 ? "PASS" : "FAIL",
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        "ON",
#else
        "OFF",
#endif
        PIC::nTotalThreads,ValidationCorrectionValues[0],
        ValidationCorrectionValues[ValidationCorrectionValueCount-1],
        ValidationCorrectionValueCount,globalMaximumCorrectionError,
        globalStage2Checks[7],globalCorrectionChecks[0],
        globalCorrectionChecks[1],globalCorrectionChecks[2],
        globalStage2Checks[4]+globalStage2Checks[5],globalStage2Checks[6],
        globalGhostChecks);
    std::fprintf(PIC::DiagnospticMessageStream,
        "%s individual-weight core compatibility correction-mode=%s mpi-ranks=%d sampling-checks=%d collision-checks=%d boundary-checks=%d restart-checks=%d migration-checks=%d collision-probability-checks=%llu collision-probability-min=%.17e collision-probability-max=%.17e\n",
        globalResult!=0 ? "PASS" : "FAIL",
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        "ON",
#else
        "OFF",
#endif
        PIC::nTotalThreads,globalCoreCorrectionChecks[0],
        globalCoreCorrectionChecks[1],globalCoreCorrectionChecks[2],
        globalCoreCorrectionChecks[3],globalCoreCorrectionChecks[4],
        globalCollisionProbabilityChecks,
        globalMinimumCollisionProbability,globalMaximumCollisionProbability);
    std::fflush(PIC::DiagnospticMessageStream);
  }
  return globalResult!=0;
#endif
}

void Comet::RadialSimulation::Run() {
  // Compare only complete, non-overlapping AMPS sampling windows. The first
  // window establishes the reference profile; each later window contributes
  // one pass or failure to the consecutive-window convergence counter.
  std::vector<double> previousProfile;
  int convergedWindows=0;
  long int completedWindows=0;
  long int previousOutputNumber=PIC::DataOutputFileNumber;
  if (Stage9ScienceDiagnosticsEnabled) {
    Stage9WindowPeakLiveParticles=PIC::ParticleBuffer::NAllPart;
    Stage9WindowPeakBufferUtilization=
        static_cast<double>(PIC::ParticleBuffer::NAllPart)/
        static_cast<double>(PIC::ParticleBuffer::MaxNPart);
    Stage9WindowStartWallTime=MPI_Wtime();
    Stage9SafeBoundaryAudits=0;
  }

  for (long int iteration=0;iteration<MaxIterations;iteration++) {
    PIC::TimeStep();
    UpdateStage9BufferPeak();
    // Cooling changes particle internal energy after this iteration's
    // transport/collisions and before the next iteration.  It deliberately
    // uses only the previously completed macroscopic sample.
    Comet::RadiativeCooling::Apply();

    // PIC::PrintData already redistributes blocks at a sampling break.  Its
    // core path resets IterationNumberAfterRebalancing to zero, so this guard
    // avoids an immediate duplicate redistribution when a science-output
    // iteration is also a multiple of the shorter scheduled interval.
    if ((iteration+1)%LoadBalanceInterval==0 &&
        PIC::Parallel::IterationNumberAfterRebalancing!=0) {
      AuditStage9SafeBoundary("load-balance");
      RebalanceRadialDomain(iteration+1);
    }

    if (PIC::DataOutputFileNumber==previousOutputNumber) {
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
      if (RegularAdaptationDue(iteration+1,false)) {
        AuditStage9SafeBoundary("regular-weight-adaptation");
        AdaptParticleWeights(iteration+1);
      }
#endif
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      // PIC::TimeStep has drained every local-timestep continuation and
      // chemistry batch.  If a base-weight adaptation is due, it has already
      // committed above; individual correction windows are therefore the
      // final carrier normalization at this safe completed-stage boundary.
      ApplyIndividualCorrectionWindowsOrAbort(iteration+1);
#endif
      continue;
    }
    previousOutputNumber=PIC::DataOutputFileNumber;

    // Alarm-triggered partial output must never satisfy a convergence test
    // configured for a full SamplingLength window.
    if (PIC::LastSampleLength!=SamplingLength) {
      if (PIC::ThisThread==0) {
        std::fprintf(PIC::DiagnospticMessageStream,
            "$PREFIX: radial convergence skipped partial window of %ld iterations\n",
            PIC::LastSampleLength);
        std::fflush(PIC::DiagnospticMessageStream);
      }
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
      if (RegularAdaptationDue(iteration+1,false))
        AdaptParticleWeights(iteration+1);
#endif
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      ApplyIndividualCorrectionWindowsOrAbort(iteration+1);
#endif
      continue;
    }

    WriteStage3CellDiagnostics(previousOutputNumber-1,PIC::LastSampleLength);
    WriteStage9CheckpointInternal(iteration+1,previousOutputNumber-1,
        PIC::LastSampleLength);
    completedWindows++;
    std::vector<double> currentProfile=CollectIntegratedDensityProfile();
    if (previousProfile.empty()) {
      Comet::Diagnostics::WriteSamplingReport(iteration+1,
          previousOutputNumber-1,PIC::LastSampleLength,
          std::vector<double>(PIC::nTotalSpecies,0.0),false);
      previousProfile.swap(currentProfile);
      if (PIC::ThisThread==0) {
        std::fprintf(PIC::DiagnospticMessageStream,
            "$PREFIX: radial convergence reference window=%ld iteration=%ld\n",
            completedWindows,iteration+1);
        std::fflush(PIC::DiagnospticMessageStream);
      }
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
      if (RegularAdaptationDue(iteration+1,true))
        AdaptParticleWeights(iteration+1);
#endif
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      ApplyIndividualCorrectionWindowsOrAbort(iteration+1);
#endif
      continue;
    }

    int limitingSpecies=0;
    const std::vector<double> relativeChanges=Comet::Diagnostics::Enabled ?
        EvaluateRelativeL1Changes(currentProfile,previousProfile) :
        std::vector<double>();
    const double maximumChange=EvaluateRelativeL1Change(currentProfile,
        previousProfile,limitingSpecies);
    Comet::Diagnostics::WriteSamplingReport(iteration+1,
        previousOutputNumber-1,PIC::LastSampleLength,relativeChanges,true);
    previousProfile.swap(currentProfile);
    if (maximumChange<=ConvergenceTolerance)
      convergedWindows++;
    else
      convergedWindows=0;

    if (PIC::ThisThread==0) {
      char symbol[_MAX_STRING_LENGTH_PIC_];
      PIC::MolecularData::GetChemSymbol(symbol,limitingSpecies);
      std::fprintf(PIC::DiagnospticMessageStream,
          "$PREFIX: radial convergence window=%ld iteration=%ld max-relative-L1=%e limiting-species=%s consecutive=%d/%d tolerance=%e\n",
          completedWindows,iteration+1,maximumChange,symbol,convergedWindows,
          RequiredConvergedWindows,ConvergenceTolerance);
      std::fflush(PIC::DiagnospticMessageStream);
    }

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
    if (RegularAdaptationDue(iteration+1,true))
      AdaptParticleWeights(iteration+1);
#endif
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    ApplyIndividualCorrectionWindowsOrAbort(iteration+1);
#endif

    if (convergedWindows>=RequiredConvergedWindows) {
      if (PIC::ThisThread==0) {
        std::fprintf(PIC::DiagnospticMessageStream,
            "$PREFIX: radial simulation converged after %ld iterations and %ld complete sampling windows\n",
            iteration+1,completedWindows);
        std::fflush(PIC::DiagnospticMessageStream);
      }
      return;
    }
  }

  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: radial simulation reached the hard limit of %ld iterations without satisfying convergence (%d/%d consecutive windows)\n",
        MaxIterations,convergedWindows,RequiredConvergedWindows);
    std::fflush(PIC::DiagnospticMessageStream);
  }
}

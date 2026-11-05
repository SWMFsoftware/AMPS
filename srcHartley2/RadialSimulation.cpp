#include "RadialSimulation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <unordered_map>
#include <vector>

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
long int Comet::RadialSimulation::MinimumParticlesPerCell=-1;
long int Comet::RadialSimulation::MaximumParticlesPerCell=-1;
double Comet::RadialSimulation::EmergencyParticleBufferUtilization=-1.0;
double Comet::RadialSimulation::AllocationLowWatermark=-1.0;
double Comet::RadialSimulation::AllocationReserveFraction=-1.0;
long int Comet::RadialSimulation::MaximumFaceTransferMultiplicity=-1;
long int Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction=-1;
double Comet::RadialSimulation::ConvergenceTolerance=-1.0;
int Comet::RadialSimulation::RequiredConvergedWindows=-1;

namespace {
bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
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

long long GetSamplingPopulationEmergencyBlockThreshold() {
  if (Comet::RadialSimulation::MaximumParticlesPerCell<=0 ||
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
  if (Comet::RadialSimulation::MaximumParticlesPerCell>
      LLONG_MAX/threshold) {
    std::fprintf(stderr,
        "sampling-population emergency per-cell threshold overflows long long\n");
    std::abort();
  }
  threshold*=Comet::RadialSimulation::MaximumParticlesPerCell;
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
  bool EmptyPhotochemicalDaughter=false;
  // True when the final power was raised above the population controller's
  // choice to impose the hard all-product reaction allocation bound.  The
  // species may be present; in that case normal factor-two roulette is
  // required rather than a weight-only edit.
  bool PhotochemicalCarrierConstraint=false;
};

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

int SelectWeightPower(long int particleCount) {
  if (particleCount<=0) return 0;

  // The range is stated per cell, while the statistical weight belongs to
  // the whole block.  Empty cells therefore contribute zero to this mean;
  // they are not seeded by the resampling operations below.
  long double expectedMean=
      static_cast<long double>(particleCount)/BlockCellCount;
  int power=0;
  while (expectedMean<Comet::RadialSimulation::MinimumParticlesPerCell) {
    expectedMean*=2.0L;
    power--;
    if (power < -1024) {
      std::fprintf(stderr,
          "particle-weight adaptation requires an invalid split exponent\n");
      std::abort();
    }
  }
  while (expectedMean>Comet::RadialSimulation::MaximumParticlesPerCell) {
    expectedMean*=0.5L;
    power++;
    if (power > 1024) {
      std::fprintf(stderr,
          "particle-weight adaptation requires an invalid roulette exponent\n");
      std::abort();
    }
  }
  return power;
}

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

int SelectMinimumCarrierWeightPower(PIC::Mesh::cDataBlockAMR *block,
    int species,double minimumCarrierRate) {
  if (!(minimumCarrierRate>0.0) || !std::isfinite(minimumCarrierRate)) {
    std::fprintf(stderr,
        "invalid minimum photochemical carrier rate for species %d\n",species);
    std::abort();
  }

  long double carrierRate=GetCarrierRate(block,species);
  int power=0;
  while (carrierRate<minimumCarrierRate) {
    carrierRate*=2.0L;
    power++;
    if (power>1024 || !std::isfinite(carrierRate)) {
      std::fprintf(stderr,
          "photochemical daughter adaptation requires an invalid weight exponent\n");
      std::abort();
    }
  }
  return power;
}

long long MaximumReactionDaughterAllocation(
    PIC::Mesh::cDataBlockAMR *block,
    const std::array<long int,PIC::nTotalSpecies> &counts) {
  if (!Comet::Chemistry::Enabled) return 0;

  // ProcessParticle first selects one mutually exclusive reaction channel and
  // then residual-resamples every concrete product with expectation
  // (W_parent/dt_parent)/(W_product/dt_product).  floor(mu)+Bernoulli(frac)
  // can allocate ceil(mu), so the hard capacity diagnostic must use that
  // upper value for every stoichiometric product, present or absent.
  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  std::array<bool,PIC::nTotalSpecies> reachable{};
  for (int species=0;species<PIC::nTotalSpecies;species++)
    // Direct surface parents can react during iteration one, before a model
    // particle of that species happens to be present in this block or before
    // the first sampling break.  Treat the configured source itself as
    // reachability evidence so initialization can prove the hard bound.
    reachable[species]=counts[species]>0 ||
        Comet::Model::SurfaceProductionRate[species]>0.0;

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
    const double parentCarrier=GetCarrierRate(block,parent);
    const ChemicalSolver::cReaction &reaction=
        network.GetReaction(reactionIndex);
    long long maximumProducts=0;
    for (std::size_t term=0;term<reaction.ProductCount;term++) {
      const ChemicalSolver::cStoichiometricTerm &product=
          network.GetProductTerm(reactionIndex,term);
      const double mean=parentCarrier/
          GetCarrierRate(block,product.Species);
      const long long oneProductMaximum=
          static_cast<long long>(std::ceil(mean));
      const long long coefficient=static_cast<long long>(
          std::floor(product.Coefficient+0.5));
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

void ConstrainPhotochemicalDaughterPlans(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,
    const std::array<long int,PIC::nTotalSpecies> &counts,
    std::array<int,PIC::nTotalSpecies> &weightPower,
    std::array<bool,PIC::nTotalSpecies> &chemistryConstrained) {
  if (!Comet::Chemistry::Enabled) return;

  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  std::array<double,PIC::nTotalSpecies> plannedCarrierRate{};
  std::array<bool,PIC::nTotalSpecies> reachable{};
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    // Population adaptation of a present parent is applied in this same
    // transaction.  Start from its planned, not old, carrier rate; otherwise
    // a roulette step that raises W would make daughter protection obsolete
    // immediately. Zero-population species have no population plan.
    plannedCarrierRate[species]=GetCarrierRate(node->block,species,
        weightPower[species]);
    reachable[species]=counts[species]>0 ||
        Comet::Model::SurfaceProductionRate[species]>0.0;
  }

  // Propagate reachability and hard carrier requirements together.  The old
  // controller considered only absent daughters and constrained an expected
  // sum.  A present daughter can have an even smaller W/dt after population
  // adaptation, and residual Bernoulli rounding can exceed an expectation.
  // Here every reachable channel is constrained by the actual worst-case sum
  // of ceil(parentCarrier/productCarrier) over all concrete products.
  bool changed=true;
  std::size_t iteration=0;
  const std::size_t maximumIterations=
      4096*PIC::nTotalSpecies*(network.GetReactionCount()+1)+1;
  while (changed) {
    changed=false;
    if (++iteration>maximumIterations) {
      std::fprintf(stderr,
          "photochemical daughter adaptation did not reach a finite carrier-rate fixed point\n");
      std::abort();
    }

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
        const ChemicalSolver::cStoichiometricTerm &product=
            network.GetProductTerm(reactionIndex,term);
        if (!reachable[product.Species]) {
          reachable[product.Species]=true;
          changed=true;
        }
      }

      auto maximumChannelAllocation = [&]() -> long long {
        long long total=0;
        for (std::size_t term=0;term<reaction.ProductCount;term++) {
          const ChemicalSolver::cStoichiometricTerm &product=
              network.GetProductTerm(reactionIndex,term);
          const long long coefficient=static_cast<long long>(
              std::floor(product.Coefficient+0.5));
          const double mean=plannedCarrierRate[parent]/
              plannedCarrierRate[product.Species];
          const long long multiplicity=
              static_cast<long long>(std::ceil(mean));
          if (multiplicity<=0 || coefficient<=0 ||
              multiplicity>LLONG_MAX/coefficient ||
              total>LLONG_MAX-multiplicity*coefficient) {
            std::fprintf(stderr,
                "invalid hard allocation bound in reaction %s\n",
                reaction.Id);
            std::abort();
          }
          total+=multiplicity*coefficient;
        }
        return total;
      };

      while (maximumChannelAllocation()>
          Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction) {
        // Double the carrier of the reducible product currently responsible
        // for the largest hard contribution.  Products identical to the
        // parent cannot reduce their own ratio and are excluded.  Repeating
        // this dyadic step yields the least invasive fixed point compatible
        // with the application's power-of-two weight policy.
        int selectedProduct=-1;
        long long selectedContribution=-1;
        for (std::size_t term=0;term<reaction.ProductCount;term++) {
          const ChemicalSolver::cStoichiometricTerm &product=
              network.GetProductTerm(reactionIndex,term);
          if (product.Species==parent) continue;
          const long long coefficient=static_cast<long long>(
              std::floor(product.Coefficient+0.5));
          const long long contribution=coefficient*
              static_cast<long long>(std::ceil(
                  plannedCarrierRate[parent]/
                  plannedCarrierRate[product.Species]));
          if (contribution>selectedContribution) {
            selectedContribution=contribution;
            selectedProduct=product.Species;
          }
        }
        if (selectedProduct<0) {
          std::fprintf(stderr,
              "reaction %s cannot satisfy the hard daughter limit by product-weight adaptation\n",
              reaction.Id);
          std::abort();
        }
        if (weightPower[selectedProduct]>=1024) {
          std::fprintf(stderr,
              "reaction %s requires an invalid product-weight exponent\n",
              reaction.Id);
          std::abort();
        }
        weightPower[selectedProduct]++;
        plannedCarrierRate[selectedProduct]*=2.0;
        if (!std::isfinite(plannedCarrierRate[selectedProduct])) {
          std::fprintf(stderr,
              "reaction %s overflowed a product carrier rate\n",reaction.Id);
          std::abort();
        }
        chemistryConstrained[selectedProduct]=true;
        changed=true;
      }
    }
  }
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

void RaiseEmptyPhotochemicalDaughterWeight(
    const cWeightAdaptationPlan &plan) {
  PIC::Mesh::cDataBlockAMR *block=plan.Node->block;
  if (!plan.EmptyPhotochemicalDaughter || plan.ParticleCount!=0 ||
      plan.WeightPower<=0 || CountSpeciesInBlock(block,plan.Species)!=0) {
    std::fprintf(stderr,
        "invalid empty photochemical daughter adaptation plan\n");
    std::abort();
  }
  const double oldWeight=block->GetLocalParticleWeight(plan.Species);
  const double newWeight=std::ldexp(oldWeight,plan.WeightPower);
  if (!(oldWeight>0.0) || !std::isfinite(oldWeight) ||
      !(newWeight>0.0) || !std::isfinite(newWeight)) {
    std::fprintf(stderr,
        "invalid empty daughter weight for species %d\n",plan.Species);
    std::abort();
  }

  // There is no particle list to resample.  Changing only the block weight
  // leaves every cell empty; when the first physical daughter is eventually
  // produced, the chemistry executor's ordinary unbiased residual sampling
  // uses this larger carrier rate and respects the configured event cap.
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

double RaiseCarrierByPowersOfTwo(double carrier,double minimumCarrier,
    int *power=NULL) {
  if (!(carrier>0.0) || !std::isfinite(carrier) ||
      !(minimumCarrier>0.0) || !std::isfinite(minimumCarrier)) {
    std::fprintf(stderr,
        "invalid carrier while enforcing the face-transfer limit\n");
    std::abort();
  }

  // All radial carrier rates start from one species-wide value and every
  // subsequent adaptation changes them by an exact power of two.  Retain
  // that invariant here.  The tolerance prevents roundoff in W/dt from
  // manufacturing an unnecessary additional doubling at an equal-carrier
  // interface.
  constexpr long double tolerance=
      128.0L*std::numeric_limits<double>::epsilon();
  long double result=carrier;
  int resultPower=0;
  while (result<static_cast<long double>(minimumCarrier)*(1.0L-tolerance)) {
    result*=2.0L;
    resultPower++;
    if (resultPower>1024 || !std::isfinite(result)) {
      std::fprintf(stderr,
          "face-carrier propagation requires an invalid weight exponent\n");
      std::abort();
    }
  }
  if (power!=NULL) *power=resultPower;
  return static_cast<double>(result);
}

cFaceCarrierPropagationStatistics EnforceFaceCarrierMultiplicityLimit() {
  cFaceCarrierPropagationStatistics statistics;

  // Build one deterministic list of calculation leaves on every MPI rank.
  // The tree topology is replicated even though live particle lists and the
  // authoritative block weights exist only on each block's owner.  Stable
  // global indexing lets the following reductions communicate population and
  // carrier information without attempting recursive point-to-point traffic
  // while sampling-break adaptation is already collective.
  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> nodes;
  CollectCalculationLeafNodes(PIC::Mesh::mesh->rootTree,nodes);
  if (nodes.empty()) {
    std::fprintf(stderr,
        "face-carrier propagation found no calculation blocks\n");
    std::abort();
  }
  if (nodes.size()*PIC::nTotalSpecies>static_cast<std::size_t>(INT_MAX)) {
    std::fprintf(stderr,
        "face-carrier propagation exceeds the MPI count limit\n");
    std::abort();
  }

  std::unordered_map<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,std::size_t>
      nodeIndex;
  nodeIndex.reserve(2*nodes.size());
  for (std::size_t block=0;block<nodes.size();block++)
    nodeIndex[nodes[block]]=block;

  // Particles cross block faces, not edges or corners.  DIM1 has one slot per
  // face, DIM2 has two transverse subfaces, and DIM3 has four.  Iterating the
  // active 2*DIM faces and all subface slots is the part that must be retained
  // when this algorithm is validated in higher-dimensional comet geometries.
  // Coarse/fine AMR interfaces may report the same neighbor more than once;
  // duplicates are harmless but are removed to keep component traversals and
  // diagnostics deterministic.
  constexpr int nFaces=2*DIM;
  constexpr int nFaceSlotsI=(DIM>=2) ? 2 : 1;
  constexpr int nFaceSlotsJ=(DIM>=3) ? 2 : 1;
  std::vector<std::vector<std::size_t> > neighbors(nodes.size());
  for (std::size_t block=0;block<nodes.size();block++) {
    for (int face=0;face<nFaces;face++)
      for (int iFace=0;iFace<nFaceSlotsI;iFace++)
        for (int jFace=0;jFace<nFaceSlotsJ;jFace++) {
          cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neighbor=
              nodes[block]->GetNeibFace(face,iFace,jFace,PIC::Mesh::mesh);
          if (neighbor==NULL) continue;
          const auto location=nodeIndex.find(neighbor);
          if (location==nodeIndex.end()) {
            std::fprintf(stderr,
                "face-carrier propagation found a face neighbor outside the leaf table\n");
            std::abort();
          }
          if (location->second==block) continue;
          std::vector<std::size_t> &list=neighbors[block];
          if (std::find(list.begin(),list.end(),location->second)==list.end())
            list.push_back(location->second);
        }
  }

  const std::size_t valueCount=nodes.size()*PIC::nTotalSpecies;
  std::vector<long long> localCounts(valueCount,0),globalCounts(valueCount,0);
  std::vector<double> localCarriers(valueCount,0.0),
      globalCarriers(valueCount,0.0);
  std::vector<int> localOwners(nodes.size(),0),globalOwners(nodes.size(),0);
  for (std::size_t block=0;block<nodes.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=nodes[block];
    if (node->Thread!=PIC::ThisThread) continue;
    if (node->block==NULL) {
      std::fprintf(stderr,
          "owner block is not allocated during face-carrier propagation\n");
      std::abort();
    }
    localOwners[block]=1;
    const std::array<long int,PIC::nTotalSpecies> counts=
        CountAllSpeciesInBlock(node->block);
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      localCounts[offset]=counts[species];
      localCarriers[offset]=GetCarrierRate(node->block,species);
    }
  }
  MPI_Allreduce(localOwners.data(),globalOwners.data(),
      static_cast<int>(nodes.size()),MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localCounts.data(),globalCounts.data(),
      static_cast<int>(valueCount),MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localCarriers.data(),globalCarriers.data(),
      static_cast<int>(valueCount),MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  for (std::size_t block=0;block<nodes.size();block++) {
    if (globalOwners[block]!=1) {
      std::fprintf(stderr,
          "calculation block has %d owners during face-carrier propagation\n",
          globalOwners[block]);
      std::abort();
    }
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      if (globalCounts[offset]<0 || !(globalCarriers[offset]>0.0) ||
          !std::isfinite(globalCarriers[offset])) {
        std::fprintf(stderr,
            "invalid global block state while propagating species %d\n",
            species);
        std::abort();
      }
    }
  }

  std::vector<double> targetCarriers=globalCarriers;
  std::vector<unsigned char> activeComponent(valueCount,0);
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    std::vector<unsigned char> visited(nodes.size(),0);
    for (std::size_t start=0;start<nodes.size();start++) {
      if (visited[start]) continue;

      // A particle can cross any face in either direction.  Consequently a
      // maximum transfer multiplicity of one implies
      //
      //   carrier(source) / carrier(destination) <= 1
      //
      // in both directed senses of every face, which is possible only when
      // all blocks in that face-connected component use the same W/dt for
      // this species.  This is deliberately stronger than the former
      // empty-neighbor repair: a populated low-weight destination can create
      // the same allocation burst as an empty one.  The explicit work list
      // is dimension-independent and cannot overflow the call stack in a
      // long radial mesh.  In DIM2/DIM3 the face table constructed above
      // naturally includes all AMR subfaces; no edge/corner coupling is
      // introduced.
      std::vector<std::size_t> component,workList(1,start);
      visited[start]=1;
      double componentCarrier=
          globalCarriers[start*PIC::nTotalSpecies+species];
      for (std::size_t cursor=0;cursor<workList.size();cursor++) {
        const std::size_t block=workList[cursor];
        component.push_back(block);
        componentCarrier=std::max(componentCarrier,
            globalCarriers[block*PIC::nTotalSpecies+species]);
        for (std::size_t neighbor:neighbors[block]) {
          if (!visited[neighbor]) {
            visited[neighbor]=1;
            workList.push_back(neighbor);
          }
        }
      }
      if (component.size()>1) statistics.ActiveComponents++;

      // All initial carrier rates for one species share a common value and
      // every adaptive change is a power of two.  The fixed-point loop finds
      // a common member of those dyadic lattices without ever lowering a
      // weight (lowering would require immediate particle cloning).  Fail
      // loudly if a future initialization or higher-dimensional extension
      // violates this invariant instead of silently retaining a burst face.
      bool targetChanged=true;
      int targetIterations=0;
      while (targetChanged) {
        targetChanged=false;
        if (++targetIterations>1024) {
          std::fprintf(stderr,
              "face-carrier propagation did not reach a dyadic fixed point\n");
          std::abort();
        }
        for (std::size_t block:component) {
          const double raised=RaiseCarrierByPowersOfTwo(
              globalCarriers[block*PIC::nTotalSpecies+species],
              componentCarrier);
          if (raised>componentCarrier*(1.0+
              128.0*std::numeric_limits<double>::epsilon())) {
            componentCarrier=raised;
            targetChanged=true;
          }
        }
      }

      // Record the graph distance reached from one component member.  The
      // validator uses this to prove propagation crossed more than the first
      // face; it is not part of the physical decision.
      std::unordered_map<std::size_t,long long> distance;
      std::vector<std::size_t> distanceWork(1,start);
      distance[start]=0;
      for (std::size_t cursor=0;cursor<distanceWork.size();cursor++) {
        const std::size_t block=distanceWork[cursor];
        const long long nextDistance=distance[block]+1;
        statistics.MaximumPropagationHops=std::max(
            statistics.MaximumPropagationHops,distance[block]);
        for (std::size_t neighbor:neighbors[block]) {
          if (distance.count(neighbor)!=0) continue;
          distance[neighbor]=nextDistance;
          distanceWork.push_back(neighbor);
        }
      }

      for (std::size_t block:component) {
        const std::size_t offset=block*PIC::nTotalSpecies+species;
        activeComponent[offset]=1;
        targetCarriers[offset]=RaiseCarrierByPowersOfTwo(
            globalCarriers[offset],componentCarrier);
      }
    }
  }

  if (Comet::Chemistry::Enabled) {
    // Face propagation can raise a parent carrier because some remote block
    // already uses a larger value.  A daughter bound computed only from the
    // population currently present in one block would then be stale, and an
    // incoming parent could create a large reaction burst before the next
    // sampling break.  Close the globally present reaction graph and impose
    // the hard event bound on the *final* component carriers.  This is
    // conservative for a future disconnected mesh (it couples components by
    // their maxima) but is exact for the connected radial domain; a future
    // DIM2/DIM3 refinement can perform the identical closure per component.
    const ChemicalSolver::cReactionNetwork &network=
        Comet::Chemistry::GetNetwork();
    std::array<bool,PIC::nTotalSpecies> reachable{};
    std::array<double,PIC::nTotalSpecies> finalCarrier{};
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      for (std::size_t block=0;block<nodes.size();block++) {
        const std::size_t offset=block*PIC::nTotalSpecies+species;
        reachable[species]=reachable[species] || globalCounts[offset]>0 ||
            Comet::Model::SurfaceProductionRate[species]>0.0;
        finalCarrier[species]=std::max(finalCarrier[species],
            targetCarriers[offset]);
      }
    }

    bool changed=true;
    std::size_t iterations=0;
    const std::size_t maximumIterations=
        4096*PIC::nTotalSpecies*(network.GetReactionCount()+1)+1;
    while (changed) {
      changed=false;
      if (++iterations>maximumIterations) {
        std::fprintf(stderr,
            "final face/chemistry carrier constraints did not converge\n");
        std::abort();
      }

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
            changed=true;
          }
        }

        auto hardAllocation = [&]() -> long long {
          long long total=0;
          for (std::size_t term=0;term<reaction.ProductCount;term++) {
            const ChemicalSolver::cStoichiometricTerm &product=
                network.GetProductTerm(reactionIndex,term);
            const long long coefficient=static_cast<long long>(
                std::floor(product.Coefficient+0.5));
            const long long multiplicity=static_cast<long long>(std::ceil(
                finalCarrier[parent]/finalCarrier[product.Species]));
            if (multiplicity<=0 || coefficient<=0 ||
                multiplicity>LLONG_MAX/coefficient ||
                total>LLONG_MAX-multiplicity*coefficient) {
              std::fprintf(stderr,
                  "invalid final hard allocation bound in reaction %s\n",
                  reaction.Id);
              std::abort();
            }
            total+=multiplicity*coefficient;
          }
          return total;
        };

        while (hardAllocation()>
            Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction) {
          int selectedProduct=-1;
          long long selectedContribution=-1;
          for (std::size_t term=0;term<reaction.ProductCount;term++) {
            const ChemicalSolver::cStoichiometricTerm &product=
                network.GetProductTerm(reactionIndex,term);
            if (product.Species==parent) continue;
            const long long coefficient=static_cast<long long>(
                std::floor(product.Coefficient+0.5));
            const long long contribution=coefficient*
                static_cast<long long>(std::ceil(finalCarrier[parent]/
                    finalCarrier[product.Species]));
            if (contribution>selectedContribution) {
              selectedContribution=contribution;
              selectedProduct=product.Species;
            }
          }
          if (selectedProduct<0 ||
              !std::isfinite(2.0*finalCarrier[selectedProduct])) {
            std::fprintf(stderr,
                "reaction %s cannot satisfy the final hard daughter limit\n",
                reaction.Id);
            std::abort();
          }
          finalCarrier[selectedProduct]*=2.0;
          changed=true;
        }
      }
    }

    // Raise every spatial component to the chemistry-safe value using the
    // same dyadic operation as the face controller.  Never lower a carrier:
    // lowering here would allocate particles before capacity is protected.
    for (int species=0;species<PIC::nTotalSpecies;species++)
      for (std::size_t block=0;block<nodes.size();block++) {
        const std::size_t offset=block*PIC::nTotalSpecies+species;
        targetCarriers[offset]=RaiseCarrierByPowersOfTwo(
            targetCarriers[offset],finalCarrier[species]);
      }
  }

  // Quantify the exact mover multiplicity that motivated this propagation on
  // every directed face, including populated-populated and populated-empty
  // pairs.  The configured production value is currently required to be one;
  // retaining it in this check makes a future bounded-multiplicity extension
  // explicit rather than hidden in a numerical tolerance.
  for (std::size_t source=0;source<nodes.size();source++)
    for (std::size_t destination:neighbors[source])
      for (int species=0;species<PIC::nTotalSpecies;species++) {
        const std::size_t sourceOffset=
            source*PIC::nTotalSpecies+species;
        const std::size_t destinationOffset=
            destination*PIC::nTotalSpecies+species;
        if (!activeComponent[destinationOffset]) continue;
        statistics.MaximumCloneMeanBefore=std::max(
            statistics.MaximumCloneMeanBefore,
            globalCarriers[sourceOffset]/globalCarriers[destinationOffset]);
        statistics.MaximumCloneMeanAfter=std::max(
            statistics.MaximumCloneMeanAfter,
            targetCarriers[sourceOffset]/targetCarriers[destinationOffset]);
      }
  const double ratioTolerance=
      Comet::RadialSimulation::MaximumFaceTransferMultiplicity*
      (1.0+256.0*std::numeric_limits<double>::epsilon());
  if (statistics.MaximumCloneMeanAfter>ratioTolerance) {
    std::fprintf(stderr,
        "face-carrier propagation left an excessive transfer multiplicity: mean=%.16e limit=%ld\n",
        statistics.MaximumCloneMeanAfter,
        Comet::RadialSimulation::MaximumFaceTransferMultiplicity);
    std::abort();
  }

  for (std::size_t block=0;block<nodes.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=nodes[block];
    if (node->Thread!=PIC::ThisThread) continue;
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      const std::size_t offset=block*PIC::nTotalSpecies+species;
      if (!activeComponent[offset]) continue;
      int power=0;
      RaiseCarrierByPowersOfTwo(globalCarriers[offset],
          targetCarriers[offset],&power);
      if (power==0) continue;
      const long int localCount=CountSpeciesInBlock(node->block,species);
      if (localCount!=globalCounts[offset]) {
        std::fprintf(stderr,
            "owner population changed during face-carrier planning for species %d\n",
            species);
        std::abort();
      }

      if (localCount==0) {
        // No list exists to resample.  Raising only the block weight keeps
        // every cell empty and ensures the first arriving particle is never
        // cloned by this face conversion.
        const double oldWeight=node->block->GetLocalParticleWeight(species);
        const double newWeight=std::ldexp(oldWeight,power);
        if (!(newWeight>0.0) || !std::isfinite(newWeight)) {
          std::fprintf(stderr,
              "invalid propagated empty-block weight for species %d\n",
              species);
          std::abort();
        }
        node->block->SetLocalParticleWeight(newWeight,species);
      }
      else {
        // A populated block must conserve its represented flux while its
        // carrier is raised.  Apply the same unbiased factor-two roulette as
        // ordinary weight adaptation once per exponent step.  No allocation
        // occurs in this branch; it can only release particle-buffer slots.
        cWeightAdaptationPlan plan;
        plan.Node=node;
        plan.Species=species;
        plan.ParticleCount=localCount;
        plan.WeightPower=power;
        for (int step=0;step<power;step++) RouletteSpeciesByTwo(plan);
      }
      statistics.LocalChangedBlockSpecies++;
    }
  }
  MPI_Allreduce(&statistics.LocalChangedBlockSpecies,
      &statistics.GlobalChangedBlockSpecies,1,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  LastFaceCarrierPropagationStatistics=statistics;
  return statistics;
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

  std::vector<cWeightAdaptationPlan> plans;
  long int requiredClones=0;
  long long localMaximumProductsBefore=0;
  long long localMaximumTotalBlockParticles=0;
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL) continue;

    // Count all species in one list traversal.  The decision is based on the
    // instantaneous population at the break, not on the completed time-
    // averaged density sample that has already been written.
    const std::array<long int,PIC::nTotalSpecies> counts=
        CountAllSpeciesInBlock(node->block);
    localMaximumTotalBlockParticles=std::max(
        localMaximumTotalBlockParticles,GetTotalBlockParticleCount(counts));
    localMaximumProductsBefore=std::max(localMaximumProductsBefore,
        MaximumReactionDaughterAllocation(node->block,counts));
    std::array<int,PIC::nTotalSpecies> finalWeightPower{};
    std::array<bool,PIC::nTotalSpecies> chemistryConstrained{};

    for (int species=0;species<PIC::nTotalSpecies;species++) {
      if (counts[species]==0) continue;
      finalWeightPower[species]=SelectWeightPower(counts[species]);
    }

    // Apply the hard reaction constraint to every reachable product, not only
    // absent species.  This may override a population-driven split: keeping a
    // reaction event below its actual allocation limit takes precedence over
    // the requested sampling range because violating it can exhaust the rank
    // before any later emergency checkpoint is reachable.
    ConstrainPhotochemicalDaughterPlans(node,counts,finalWeightPower,
        chemistryConstrained);

    for (int species=0;species<PIC::nTotalSpecies;species++) {
      if (finalWeightPower[species]==0) continue;
      cWeightAdaptationPlan plan;
      plan.Node=node;
      plan.Species=species;
      plan.ParticleCount=counts[species];
      plan.WeightPower=finalWeightPower[species];
      plan.EmptyPhotochemicalDaughter=counts[species]==0 &&
          chemistryConstrained[species];
      plan.PhotochemicalCarrierConstraint=chemistryConstrained[species];

      const double oldWeight=node->block->GetLocalParticleWeight(species);
      const double newWeight=std::ldexp(oldWeight,plan.WeightPower);
      if (!(oldWeight>0.0) || !std::isfinite(oldWeight) ||
          !(newWeight>0.0) || !std::isfinite(newWeight)) {
        std::fprintf(stderr,
            "invalid planned particle weight for species %d\n",species);
        std::abort();
      }
      const long int clones=RequiredCloneCount(plan);
      if (requiredClones>LONG_MAX-clones) {
        std::fprintf(stderr,
            "particle-weight adaptation total clone count overflows long int\n");
        std::abort();
      }
      requiredClones+=clones;
      plans.push_back(plan);
    }
  }

  // The regular controller above treats each species independently.  A
  // pathological block can nevertheless contain a very large sum over many
  // species, or a single reaction-driven population can overshoot far beyond
  // the requested range between two sampling windows.  Detect that condition
  // from the pre-adaptation snapshot.  All ranks store the same global result
  // so the later collective callback cannot diverge across MPI processes.
  long long globalMaximumTotalBlockParticles=0;
  MPI_Allreduce(&localMaximumTotalBlockParticles,
      &globalMaximumTotalBlockParticles,1,MPI_LONG_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  const long long samplingEmergencyBlockThreshold=
      GetSamplingPopulationEmergencyBlockThreshold();
  if (globalMaximumTotalBlockParticles>samplingEmergencyBlockThreshold) {
    SamplingPopulationEmergencyPending=true;
    SamplingPopulationEmergencyMaximumBlockParticles=
        globalMaximumTotalBlockParticles;
    SamplingPopulationEmergencyThresholdBlockParticles=
        samplingEmergencyBlockThreshold;
    SamplingPopulationEmergencyIteration=completedIterations;

    if (PIC::ThisThread==0) {
      std::fprintf(PIC::DiagnospticMessageStream,
          "$PREFIX: radial sampling-population emergency requested iteration=%ld max-mean-total-particles-per-cell=%.6e threshold=%.6e multiplier=%lld species=%d execution=before-next-regular-mpi-particle-exchange\n",
          completedIterations,
          static_cast<double>(globalMaximumTotalBlockParticles)/
              BlockCellCount,
          static_cast<double>(samplingEmergencyBlockThreshold)/
              BlockCellCount,
          SamplingPopulationEmergencyMultiplier,PIC::nTotalSpecies);
      std::fflush(PIC::DiagnospticMessageStream);
    }
  }

  // Preflight against the per-rank buffer before changing any weight or list.
  // This deliberately ignores particles that roulette may free, making the
  // check conservative and preventing a half-applied split transaction.
  const long int available=
      PIC::ParticleBuffer::MaxNPart-PIC::ParticleBuffer::NAllPart;
  const int localCapacityFailure=requiredClones>available ? 1 : 0;
  int globalCapacityFailure=0;
  MPI_Allreduce(&localCapacityFailure,&globalCapacityFailure,1,MPI_INT,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  if (globalCapacityFailure!=0) {
    if (localCapacityFailure!=0)
      std::fprintf(stderr,
          "rank %d needs %ld adaptation clones but has only %ld particle slots\n",
          PIC::ThisThread,requiredClones,available);
    std::abort();
  }

  const long long localBefore=PIC::ParticleBuffer::NAllPart;
  long long localSplitPlans=0,localRoulettePlans=0;
  long long localEmptyDaughterPlans=0;
  long long localChemistryConstrainedPlans=0;
  for (const cWeightAdaptationPlan &plan:plans) {
    if (plan.PhotochemicalCarrierConstraint)
      localChemistryConstrainedPlans++;
    if (plan.EmptyPhotochemicalDaughter) {
      localEmptyDaughterPlans++;
      RaiseEmptyPhotochemicalDaughterWeight(plan);
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

  // Population and chemistry plans above modify only their owner block.
  // Before the next mover step, enforce the configured multiplicity on every
  // directed face.  At the validated value of one this equalizes W/dt over a
  // complete face-connected component.  Empty blocks receive only a weight;
  // populated blocks are rouletted, so the operation allocates no particles.
  // This closes both empty- and populated-destination transport bursts.
  const cFaceCarrierPropagationStatistics propagationStatistics=
      EnforceFaceCarrierMultiplicityLimit();

  long long localMaximumProductsAfter=0;
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL) continue;
    const std::array<long int,PIC::nTotalSpecies> counts=
        CountAllSpeciesInBlock(node->block);
    localMaximumProductsAfter=std::max(localMaximumProductsAfter,
        MaximumReactionDaughterAllocation(node->block,counts));
  }
  if (Comet::Chemistry::Enabled &&
      localMaximumProductsAfter>
          Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction) {
    std::fprintf(stderr,
        "rank %d photochemical adaptation left %lld possible products per event; limit=%ld\n",
        PIC::ThisThread,localMaximumProductsAfter,
        Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction);
    std::abort();
  }

  long long localStatistics[7]={localBefore,PIC::ParticleBuffer::NAllPart,
      static_cast<long long>(plans.size()),localSplitPlans,
      localRoulettePlans,localEmptyDaughterPlans,
      localChemistryConstrainedPlans};
  long long globalStatistics[7]={0,0,0,0,0,0,0};
  MPI_Allreduce(localStatistics,globalStatistics,7,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  long long localProductStatistics[2]={localMaximumProductsBefore,
      localMaximumProductsAfter};
  long long globalProductStatistics[2]={0,0};
  MPI_Allreduce(localProductStatistics,globalProductStatistics,2,MPI_LONG_LONG,
      MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  if (globalStatistics[2]>0 ||
      propagationStatistics.GlobalChangedBlockSpecies>0)
    SynchronizeAdaptedBlockWeights();
  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: radial particle-weight adaptation iteration=%ld range=%ld:%ld particles=%lld->%lld changed-block-species=%lld split=%lld roulette=%lld empty-photochemical-daughter=%lld chemistry-constrained=%lld face-carrier-changes=%lld face-components=%lld max-face-hops=%lld max-face-clone-mean=%.6e->%.6e max-products-per-reaction=%lld->%lld limit=%ld\n",
        completedIterations,
        Comet::RadialSimulation::MinimumParticlesPerCell,
        Comet::RadialSimulation::MaximumParticlesPerCell,
        globalStatistics[0],globalStatistics[1],globalStatistics[2],
        globalStatistics[3],globalStatistics[4],globalStatistics[5],
        globalStatistics[6],
        propagationStatistics.GlobalChangedBlockSpecies,
        propagationStatistics.ActiveComponents,
        propagationStatistics.MaximumPropagationHops,
        propagationStatistics.MaximumCloneMeanBefore,
        propagationStatistics.MaximumCloneMeanAfter,
        globalProductStatistics[0],globalProductStatistics[1],
        Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction);
    std::fflush(PIC::DiagnospticMessageStream);
  }
}

void SeedValidationParticles(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int species,
    long int particleCount) {
  if (node==NULL || node->block==NULL) return;
  long int &first=node->block->FirstCellParticleTable[0];
  const double cellWidth=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
  const double position[3]={node->xmin[0]+0.5*cellWidth,0.0,0.0};
  const double velocity[3]={0.0,0.0,0.0};
  for (long int particle=0;particle<particleCount;particle++) {
    const long int ptr=PIC::ParticleBuffer::GetNewParticle(first);
    PIC::ParticleBuffer::byte *data=
        PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    PIC::ParticleBuffer::SetX(const_cast<double*>(position),data);
    PIC::ParticleBuffer::SetV(const_cast<double*>(velocity),data);
    PIC::ParticleBuffer::SetI(species,data);
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
    PIC::IDF::InitRotTemp(
        Comet::RadialSimulation::InnerBoundaryTemperature,data);
    PIC::IDF::InitVibTemp(
        Comet::RadialSimulation::InnerBoundaryTemperature,data);
#endif
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
    PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
#endif
  }
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

void FillLocalCellMeasureAudit(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,long int &cellIndex,
    std::vector<double> &radius,std::vector<double> &measure) {
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
        measure[cellIndex+i]=cell->Measure;
      }
    }
    cellIndex+=_BLOCK_CELLS_X_;
    return;
  }
  for (int child=0;child<(1<<DIM);child++)
    FillLocalCellMeasureAudit(node->downNode[child],cellIndex,radius,measure);
}

void WriteCellMeasureAudit() {
  if (!Comet::Diagnostics::Enabled) return;
  const long int cellCount=CountProfileCells(PIC::Mesh::mesh->rootTree);
  std::vector<double> localRadius(cellCount,0.0),localMeasure(cellCount,0.0);
  std::vector<double> globalRadius(cellCount,0.0),globalMeasure(cellCount,0.0);
  long int visited=0;
  FillLocalCellMeasureAudit(PIC::Mesh::mesh->rootTree,visited,localRadius,
      localMeasure);
  if (visited!=cellCount || cellCount>INT_MAX) {
    std::fprintf(stderr,"invalid radial cell-measure audit traversal\n");
    std::abort();
  }
  MPI_Reduce(localRadius.data(),globalRadius.data(),static_cast<int>(cellCount),
      MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localMeasure.data(),globalMeasure.data(),
      static_cast<int>(cellCount),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    FILE *file=std::fopen("PT/diagnostics/cell_measures.csv","w");
    if (file==NULL) {
      std::fprintf(stderr,"cannot open radial cell-measure audit file\n");
      std::abort();
    }
    std::fprintf(file,"cell_index,radius_m,measure_m3\n");
    for (long int cell=0;cell<cellCount;cell++) {
      if (!(globalMeasure[cell]>0.0) || !std::isfinite(globalRadius[cell])) {
        std::fprintf(stderr,"incomplete reduced radial cell-measure audit\n");
        std::abort();
      }
      std::fprintf(file,"%ld,%.17e,%.17e\n",cell,globalRadius[cell],
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
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
    return Fail("particle-weight adaptation requires individual weight correction off",error);
#else
    if (MinimumParticlesPerCell<=0 || MaximumParticlesPerCell<=0)
      return Fail("particle-weight adaptation range must be positive",error);
    if (MaximumParticlesPerCell<2*MinimumParticlesPerCell)
      return Fail("particle-weight adaptation maximum must be at least twice the minimum",error);
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
      return Fail("the validated correction-free face-allocation controller currently requires maximum face multiplicity 1",error);
    if (Comet::Chemistry::Enabled &&
        MaximumPhotochemicalProductsPerReaction<=0)
      return Fail("particle-weight adaptation with chemistry requires a positive maximum photochemical product count per reaction",error);
#endif
  }
  if (!(ConvergenceTolerance>0.0) || ConvergenceTolerance>1.0 ||
      !std::isfinite(ConvergenceTolerance))
    return Fail("radial convergence tolerance must be in (0,1]",error);
  if (RequiredConvergedWindows<=0)
    return Fail("required converged-window count must be positive",error);

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
#endif

  PIC::BC::BlockInjectionBCindicatior=InjectionBlockIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=InjectBlock;
  PIC::BC::InitBoundingBoxInjectionBlockList();
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

  constexpr int species=0;
  const double originalWeight=localHasBlock ?
      node->block->GetLocalParticleWeight(species) : 1.0;
  const double originalCarrier=localHasBlock ?
      GetCarrierRate(node->block,species) : 0.0;
  double globalOriginalCarrier=0.0;
  MPI_Allreduce(&originalCarrier,&globalOriginalCarrier,1,MPI_DOUBLE,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  bool localPass=true;
  int localGhostChecks=0;
  int localFaceCarrierPropagationChecks=0;

  // Splitting test: choose a nonzero block mean well below the lower bound.
  // The selector decides the required power, so the test remains valid for
  // any mathematically supported input range.  All particles begin in cell
  // zero; the other cells must remain empty after every split.
  PIC::ParticleBuffer::DeleteAllParticles();
  const long int splitInitial=std::max(1L,
      (MinimumParticlesPerCell*BlockCellCount+7)/8);
  const int splitPower=SelectWeightPower(splitInitial);
  const double splitTargetCarrier=
      std::ldexp(globalOriginalCarrier,splitPower);
  // Put every other block at the carrier the selected blocks will reach.
  // Otherwise the strict face controller would correctly raise the selected
  // carrier after splitting and its roulette would obscure this focused test
  // of the split primitive.
  SetWeightOverTimeStep(species,splitTargetCarrier,
      PIC::Mesh::mesh->rootTree);
  if (localHasBlock) {
    node->block->SetLocalParticleWeight(originalWeight,species);
    SeedValidationParticles(node,species,splitInitial);
  }
  SynchronizeAdaptedBlockWeights();
  AdaptParticleWeights(0);
  if (localHasBlock) {
    long int multiplier=1;
    for (int step=0;step < -splitPower;step++) multiplier*=2;
    const long int expectedCount=splitInitial*multiplier;
    const long int actualCount=CountSpeciesInBlock(node->block,species);
    const double expectedWeight=std::ldexp(originalWeight,splitPower);
    const double actualWeight=node->block->GetLocalParticleWeight(species);
    const double weightScale=std::max(std::fabs(expectedWeight),1.0);
    localPass=localPass && splitPower<0 && actualCount==expectedCount &&
        std::fabs(actualWeight-expectedWeight)<=1.0e-13*weightScale &&
        OtherValidationCellsAreEmpty(node->block);
  }
  localPass=ValidateSynchronizedWeightCopies(node,species,
      std::ldexp(originalWeight,splitPower),localGhostChecks) && localPass;

  // Roulette test: twice the upper-bound mean requires exactly one weight
  // doubling.  The survivor count is stochastic, so use a six-sigma binomial
  // bound while also checking that no empty cell was populated.
  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,2.0*globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  if (localHasBlock) node->block->SetLocalParticleWeight(
      originalWeight,species);
  SynchronizeAdaptedBlockWeights();
  const long int rouletteInitial=
      2*MaximumParticlesPerCell*BlockCellCount;
  if (localHasBlock)
    SeedValidationParticles(node,species,rouletteInitial);
  rnd_seed(99173+PIC::ThisThread);
  AdaptParticleWeights(0);
  if (localHasBlock) {
    const long int actualCount=CountSpeciesInBlock(node->block,species);
    const double expectation=0.5*rouletteInitial;
    const double sigma=std::sqrt(0.25*rouletteInitial);
    const double actualWeight=node->block->GetLocalParticleWeight(species);
    const double expectedWeight=2.0*originalWeight;
    const double weightScale=std::max(std::fabs(expectedWeight),1.0);
    localPass=localPass &&
        std::fabs(actualCount-expectation)<=6.0*sigma &&
        std::fabs(actualWeight-expectedWeight)<=1.0e-13*weightScale &&
        OtherValidationCellsAreEmpty(node->block);
  }
  localPass=ValidateSynchronizedWeightCopies(node,species,
      2.0*originalWeight,localGhostChecks) && localPass;
  // Test recursive face propagation separately from population adaptation.
  // All cells are empty; selected blocks start at twice the common carrier.
  // A directed-face limit of one must propagate that carrier through each
  // complete connected component without allocating any particle.
  PIC::ParticleBuffer::DeleteAllParticles();
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  if (localHasBlock)
    node->block->SetLocalParticleWeight(2.0*originalWeight,species);
  SynchronizeAdaptedBlockWeights();
  AdaptParticleWeights(0);
  if (localHasBlock) {
    const cFaceCarrierPropagationStatistics &propagation=
        LastFaceCarrierPropagationStatistics;
    const bool propagationPass=
        propagation.GlobalChangedBlockSpecies>0 &&
        propagation.ActiveComponents>0 &&
        propagation.MaximumPropagationHops>=2 &&
        propagation.MaximumCloneMeanBefore>1.0 &&
        propagation.MaximumCloneMeanAfter<=
            1.0+256.0*std::numeric_limits<double>::epsilon();
    localPass=localPass && propagationPass &&
        CountSpeciesInBlock(node->block,species)==0;
    if (propagationPass) localFaceCarrierPropagationChecks++;
  }
  SetWeightOverTimeStep(species,globalOriginalCarrier,
      PIC::Mesh::mesh->rootTree);
  SynchronizeAdaptedBlockWeights();

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
      bool hasDistinctProduct=false;
      for (std::size_t term=0;term<reaction.ProductCount;term++)
        if (network.GetProductTerm(reactionIndex,term).Species!=parent)
          hasDistinctProduct=true;
      if (hasDistinctProduct) {
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
        const long int parentCount=MinimumParticlesPerCell*BlockCellCount;
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

      // Force the hard residual-allocation bound—not only its expectation—
      // above the configured cap.  Install the unsafe carrier on every mesh
      // copy so the subsequent face controller cannot mask the chemistry
      // preparation test with an older, larger carrier in another block.
      const double unsafeCarrier=savedCarriers[parent]*absentCoefficient/
          (8.0*MaximumPhotochemicalProductsPerReaction);
      for (int product=0;product<PIC::nTotalSpecies;product++)
        if (testedDaughter[product])
          SetWeightOverTimeStep(product,unsafeCarrier,
              PIC::Mesh::mesh->rootTree);
      SynchronizeAdaptedBlockWeights();

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
        bool daughtersRemainEmpty=true;
        for (int product=0;product<PIC::nTotalSpecies;product++) {
          if (!testedDaughter[product]) continue;
          daughtersRemainEmpty=daughtersRemainEmpty && counts[product]==0;
          localPass=ValidateSynchronizedWeightCopies(node,product,
              node->block->GetLocalParticleWeight(product),localGhostChecks) &&
              localPass;
        }
        localPass=localPass &&
            multiplicityBefore>MaximumPhotochemicalProductsPerReaction &&
            multiplicityAfter<=MaximumPhotochemicalProductsPerReaction &&
            daughtersRemainEmpty;
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
  constexpr long int validationMinimumParticlesPerCell=5;
  constexpr long int validationMaximumParticlesPerCell=20;
  MinimumParticlesPerCell=validationMinimumParticlesPerCell;
  MaximumParticlesPerCell=validationMaximumParticlesPerCell;
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
    rnd_seed(99241+PIC::ThisThread);
    PIC::Parallel::ExchangeParticleData();
    const bool requestConsumed=!SamplingPopulationEmergencyPending;

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
  if (PIC::nTotalThreads>1 && globalGhostChecks==0) globalResult=0;
  if (Comet::Chemistry::Enabled &&
      globalPhotochemicalChecks!=globalBlockCount) globalResult=0;
  if (globalEmergencyChecks!=globalBlockCount) globalResult=0;
  if (globalAllocationBurstChecks!=globalBlockCount) globalResult=0;
  if (globalSamplingPopulationEmergencyChecks!=globalBlockCount)
    globalResult=0;
  if (globalFaceCarrierPropagationChecks!=globalBlockCount)
    globalResult=0;
  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "%s particle-weight adaptation validation range=%ld:%ld emergency-buffer-utilization=%.6f allocation-low-watermark=%.6f allocation-reserve=%.6f max-face-transfer-multiplicity=%ld max-photochemical-products=%ld owner-blocks=%d emergency-checks=%d allocation-burst-checks=%d sampling-population-emergency-checks=%d face-carrier-propagation-checks=%d photochemical-checks=%d ghost-weight-checks=%d\n",
        globalResult!=0 ? "PASS" : "FAIL",MinimumParticlesPerCell,
        MaximumParticlesPerCell,EmergencyParticleBufferUtilization,
        AllocationLowWatermark,AllocationReserveFraction,
        MaximumFaceTransferMultiplicity,
        MaximumPhotochemicalProductsPerReaction,globalBlockCount,
        globalEmergencyChecks,globalAllocationBurstChecks,
        globalSamplingPopulationEmergencyChecks,
        globalFaceCarrierPropagationChecks,
        globalPhotochemicalChecks,globalGhostChecks);
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

  for (long int iteration=0;iteration<MaxIterations;iteration++) {
    PIC::TimeStep();
    // Cooling changes particle internal energy after this iteration's
    // transport/collisions and before the next iteration.  It deliberately
    // uses only the previously completed macroscopic sample.
    Comet::RadiativeCooling::Apply();

    // PIC::PrintData already redistributes blocks at a sampling break.  Its
    // core path resets IterationNumberAfterRebalancing to zero, so this guard
    // avoids an immediate duplicate redistribution when a science-output
    // iteration is also a multiple of the shorter scheduled interval.
    if ((iteration+1)%LoadBalanceInterval==0 &&
        PIC::Parallel::IterationNumberAfterRebalancing!=0)
      RebalanceRadialDomain(iteration+1);

    if (PIC::DataOutputFileNumber==previousOutputNumber) continue;
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
      continue;
    }

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
      AdaptParticleWeights(iteration+1);
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

    if (convergedWindows>=RequiredConvergedWindows) {
      if (PIC::ThisThread==0) {
        std::fprintf(PIC::DiagnospticMessageStream,
            "$PREFIX: radial simulation converged after %ld iterations and %ld complete sampling windows\n",
            iteration+1,completedWindows);
        std::fflush(PIC::DiagnospticMessageStream);
      }
      return;
    }
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
    AdaptParticleWeights(iteration+1);
#endif
  }

  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: radial simulation reached the hard limit of %ld iterations without satisfying convergence (%d/%d consecutive windows)\n",
        MaxIterations,convergedWindows,RequiredConvergedWindows);
    std::fflush(PIC::DiagnospticMessageStream);
  }
}

#include "MoonParticleBudget.h"
#include "Moon.h"
#include "MoonInput.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm>

namespace Moon { namespace ParticleBudget {
namespace {
using namespace PIC::ParticleEvents;
constexpr int Events=EventCount;
std::array<long long,Events> localCount{},globalCount{},previousCount{};
std::array<long double,Events> localAtoms{},globalAtoms{},previousAtoms{};
std::array<long long,1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_> sourceCount{};
std::array<long double,1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_> sourceAtoms{};
long long previousGasCount=0,lastGasCount=0,maximumResidual=0;
long double previousGasAtoms=0.0L;
bool active=false,valid=true;
std::ofstream ledger,sources;
const char* names[Events]={"emission","adsorption","reflection",
    "photolytic_loss","photolytic_product","escape","numerical_removal",
    "reservoir_supply","surface_deposit","surface_release"};

bool GasEvent(int event) {
  return event!=Reflection && event!=ReservoirSupply &&
      event!=SurfaceDeposit && event!=SurfaceRelease;
}
int Sign(int event) {
  return event==Emission || event==PhotolyticProduct ? 1 : -1;
}

void Observe(Event event,int species,double weight,int detail) {
  if (species!=_NA_SPEC_) return;
  if (!std::isfinite(weight) || (event!=ReservoirSupply && !(weight>0.0)))
    exit(__LINE__,__FILE__,"Na event ledger received an invalid weight [atoms]");
  // MPI ranks own disjoint event streams. A named critical section also
  // makes the observer safe if a supported build executes threaded callbacks.
#ifdef _OPENMP
#pragma omp critical (MoonParticleEventLedger)
#endif
  {
    ++localCount[event];
    localAtoms[event]+=static_cast<long double>(weight);
    if (event==Emission) {
      if (detail<0 || detail>_EXOSPHERE__SOURCE_MAX_ID_VALUE_)
        exit(__LINE__,__FILE__,"Na event ledger received an unknown source ID");
      ++sourceCount[detail];
      sourceAtoms[detail]+=static_cast<long double>(weight);
    }
  }
}

std::array<long double,2> Inventory() {
  long double local[2]={0.0L,0.0L},global[2];
  for (auto* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node!=NULL;node=node->nextNodeThisThread) {
    if (node->block==NULL) continue;
    for (int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for (long int particle=node->block->FirstCellParticleTable[cell];
           particle!=-1;particle=PIC::ParticleBuffer::GetNext(particle)) {
        if (PIC::ParticleBuffer::GetI(particle)!=_NA_SPEC_) continue;
        double weight=node->block->GetLocalParticleWeight(_NA_SPEC_);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        weight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle);
#endif
        local[0]+=1.0L;
        local[1]+=static_cast<long double>(weight);
      }
  }
  MPI_Allreduce(local,global,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  return {global[0],global[1]};
}
}

void Begin(const std::filesystem::path& directory) {
  if (PIC::ParticleEvents::Observer!=NULL)
    exit(__LINE__,__FILE__,"Na event ledger refuses to replace an existing event observer");
  if (PIC::nTotalSpecies!=1 || _NA_SPEC_<0 ||
      !Moon::Runtime::HasConfiguration() ||
      Moon::Runtime::GetConfiguration().surfaceGeometry!=Moon::Runtime::SurfaceGeometry::Sphere)
    exit(__LINE__,__FILE__,"I11 event ledger currently requires the neutral-Na sphere input");
  localCount.fill(0);globalCount.fill(0);previousCount.fill(0);
  localAtoms.fill(0);globalAtoms.fill(0);previousAtoms.fill(0);
  sourceCount.fill(0);sourceAtoms.fill(0);
  maximumResidual=0;valid=true;
  const auto inventory=Inventory();
  previousGasCount=lastGasCount=static_cast<long long>(inventory[0]);
  previousGasAtoms=inventory[1];
  if (PIC::ThisThread==0) {
    ledger.open(directory/"particle_budget.csv");
    sources.open(directory/"particle_sources.csv");
    if (!ledger || !sources)
      exit(__LINE__,__FILE__,"cannot open Na particle budget/source diagnostics");
    ledger << std::setprecision(std::numeric_limits<long double>::max_digits10)
           << "checkpoint,completed_particle_steps,elapsed_particle_time_s,gas_model_particles,gas_atoms,surface_atoms,model_count_residual,weighted_residual_atoms,weighted_roundoff_bound_atoms";
    for (const char* name:names) ledger << ',' << name << "_count," << name << "_atoms";
    ledger << '\n';
    sources << std::setprecision(std::numeric_limits<long double>::max_digits10)
            << "checkpoint,source_id,source_symbol,cumulative_emitted_model_particles,cumulative_emitted_atoms\n";
  }
  active=true;
  PIC::ParticleEvents::Observer=Observe;
  Snapshot(0,"initial");
}

void Snapshot(long long completedParticleSteps,const char* checkpoint) {
  if (!active) return;
  const auto inventory=Inventory();
  lastGasCount=static_cast<long long>(inventory[0]);
  MPI_Allreduce(localCount.data(),globalCount.data(),Events,MPI_LONG_LONG,
                MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localAtoms.data(),globalAtoms.data(),Events,MPI_LONG_DOUBLE,
                MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  long long expectedCount=0,operations=lastGasCount+previousGasCount;
  long double expectedAtoms=0.0L,absoluteTerms=
      std::fabs(inventory[1])+std::fabs(previousGasAtoms);
  for (int event=0;event<Events;++event) if (GasEvent(event)) {
    expectedCount+=Sign(event)*(globalCount[event]-previousCount[event]);
    expectedAtoms+=Sign(event)*(globalAtoms[event]-previousAtoms[event]);
    absoluteTerms+=std::fabs(globalAtoms[event])+std::fabs(previousAtoms[event]);
    operations+=globalCount[event];
  }
  const long long residual=lastGasCount-previousGasCount-expectedCount;
  const long double weightedResidual=inventory[1]-previousGasAtoms-expectedAtoms;
  // A conservative forward-error bound includes both cumulative sums, their
  // subtraction, two list sums and distributed reductions. Values are sums of
  // actual double weights promoted exactly to long double; no fitted physical
  // tolerance or expected Poisson rate enters this bookkeeping identity.
  const long double n=4.0L*operations+32.0L*PIC::nTotalThreads*(Events+1);
  const long double nu=n*std::numeric_limits<long double>::epsilon();
  if (!(nu<1.0L)) exit(__LINE__,__FILE__,"particle ledger roundoff bound exhausted");
  const long double bound=nu/(1.0L-nu)*absoluteTerms;
  maximumResidual=std::max(maximumResidual,std::llabs(residual));
  valid=valid && residual==0 && std::fabs(weightedResidual)<=bound &&
      globalCount[NumericalRemoval]==0;
  long long allSources[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];
  long double allSourceAtoms[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];
  MPI_Reduce(sourceCount.data(),allSources,sourceCount.size(),MPI_LONG_LONG,
             MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(sourceAtoms.data(),allSourceAtoms,sourceAtoms.size(),MPI_LONG_DOUBLE,
             MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    long double surfaceAtoms=0.0L;
    // Logical surface populations are already replicated by their production
    // MPI exchange. Summing them over ranks would multiply the reservoir.
    for (int element=0;element<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;++element)
      surfaceAtoms+=Moon::Planet->SurfaceElementPopulation[_NA_SPEC_][element];
    ledger << checkpoint << ',' << completedParticleSteps << ','
           << completedParticleSteps*PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_]
           << ',' << lastGasCount << ',' << inventory[1] << ',' << surfaceAtoms
           << ',' << residual << ',' << weightedResidual << ',' << bound;
    for (int event=0;event<Events;++event)
      ledger << ',' << globalCount[event]-previousCount[event] << ','
             << globalAtoms[event]-previousAtoms[event];
    ledger << '\n';ledger.flush();
    for (int source=0;source<=_EXOSPHERE__SOURCE_MAX_ID_VALUE_;++source)
      sources << checkpoint << ',' << source << ','
              << _EXOSPHERE__SOURCE_SYMBOLIC_ID_[source] << ','
              << allSources[source] << ',' << allSourceAtoms[source] << '\n';
    sources.flush();
    if (!ledger || !sources) exit(__LINE__,__FILE__,"Na ledger write failed; check output disk capacity");
  }
  previousGasCount=lastGasCount;previousGasAtoms=inventory[1];
  previousCount=globalCount;previousAtoms=globalAtoms;
}

void End() {
  if (active) PIC::ParticleEvents::Observer=NULL;
  active=false;ledger.close();sources.close();
}
bool Valid() { return valid; }
long long MaximumCountResidual() { return maximumResidual; }
long long Count(PIC::ParticleEvents::Event event) { return globalCount[event]; }
long long GasCount() { return lastGasCount; }
} }

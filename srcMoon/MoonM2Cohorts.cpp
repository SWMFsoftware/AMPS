#include "pic.h"
#include "Moon.h"
#include "MoonM2Cohorts.h"
#include "MoonSurfaceTemperature.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

namespace Moon { namespace M2Cohorts { namespace {
long int offset=-1;
int K=0,E=0;
std::vector<int> sites;
std::vector<long double> stock,pending,initial,first,restricted,exposure,loss;
std::vector<double> packetFractions;
std::vector<long long> localOperations,replicatedOperations;
std::ofstream output,hazards;
int currentSpecies=-1,currentElement=-1;
double elapsed=0.0;
std::size_t Index(int species,int element,int cohort) {
  return (static_cast<std::size_t>(species)*E+element)*K+cohort;
}
int Species(int s) {
  if (s==_NA_SPEC_) return 0;
  if (s==_AR_SPEC_) return 1;
  exit(__LINE__,__FILE__,"M2 cohort observer encountered unsupported species");
  return -1;
}
double* Fractions(long int p) {
  return reinterpret_cast<double*>(PIC::ParticleBuffer::GetParticleDataPointer(p)+offset);
}
void BeforeRelease(int spec,int element,double population,double remaining,
                   double hazard,double dt,double thermal,double photon) {
  // Identical on all ranks: the physical reservoir is replicated. This
  // updates only passive composition; the finite injector owns actual stock.
  const int s=Species(spec);
  if (!(population>0) || remaining<0 || remaining>population || hazard<0 || dt<=0)
    exit(__LINE__,__FILE__,"invalid actual finite-release cohort interval");
  currentSpecies=spec;currentElement=element;
  const long double survival=static_cast<long double>(remaining)/population;
  for (int c=0;c<K;++c) {
    const auto i=Index(s,element,c);replicatedOperations[s*K+c]+=4;
    packetFractions[c]=static_cast<double>(stock[i]/population);
    stock[i]*=survival;
    if (sites[c]!=element) continue;
    const int j=s*K+c;replicatedOperations[j]+=16;
    // Integral of never-released stock over this piecewise-constant hazard.
    // Returned atoms never re-enter this first-episode survival population.
    const long double integral=hazard==0 ? dt : -std::expm1(-hazard*dt)/hazard;
    restricted[j]+=first[j]*integral;
    first[j]*=survival;exposure[j]+=hazard*dt;
    if (PIC::ThisThread==0)
      hazards << std::setprecision(17) << elapsed << ',' << dt << ','
        << (s==0?"Na":"Ar") << ',' << c << ',' << element << ','
        << Moon::SurfaceTemperature::ProductionElementTemperatureK(element)
        << ',' << hazard << ',' << thermal << ',' << photon << '\n';
  }
}
void Packet(int spec,int element,long int particle,double) {
  if (spec!=currentSpecies || element!=currentElement)
    exit(__LINE__,__FILE__,"cohort packet does not match collective release element");
  std::memcpy(Fractions(particle),packetFractions.data(),K*sizeof(double));
}
void Deposit(int spec,int element,long int particle,double atoms) {
  const int s=Species(spec);const double* f=Fractions(particle);
  for (int c=0;c<K;++c) {pending[Index(s,element,c)]+=atoms*f[c];localOperations[s*K+c]+=2;}
}
void ParticleEvent(PIC::ParticleEvents::Event event,int spec,long int p,double atoms,int) {
  if (event!=PIC::ParticleEvents::Escape && event!=PIC::ParticleEvents::NumericalRemoval) return;
  const int s=Species(spec);const double* f=Fractions(p);
  const int kind=event==PIC::ParticleEvents::Escape ? 0:1;
  for (int c=0;c<K;++c) {loss[(kind*2+s)*K+c]+=atoms*f[c];localOperations[s*K+c]+=2;}
}
} // namespace

void RegisterStorage(int cohorts) {
  if (offset!=-1 || cohorts<=0 || cohorts>48)
    exit(__LINE__,__FILE__,"M2 cohort field registration must occur once for 1..48 cohorts");
  K=cohorts;
  PIC::ParticleBuffer::RequestDataStorage(offset,K*sizeof(double),alignof(double));
}
void Begin(const std::vector<int>& elements,const std::vector<double>& initialAtoms,
           const std::filesystem::path& directory) {
  using namespace Exosphere::SourceProcesses;
  if (offset<0 || K!=static_cast<int>(elements.size()) || initialAtoms.size()!=2u*K ||
      FiniteReleaseBegin || FiniteReleasePacket || FiniteReleaseDeposit ||
      PIC::ParticleEvents::ParticleObserver)
    exit(__LINE__,__FILE__,"M2 cohorts require registered storage and unclaimed passive observers");
  E=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;sites=elements;
  stock.assign(2u*E*K,0);pending=stock;initial.assign(initialAtoms.begin(),initialAtoms.end());
  first=initial;restricted.assign(2*K,0);exposure.assign(2*K,0);loss.assign(4*K,0);
  packetFractions.assign(K,0);elapsed=0;
  localOperations.assign(2*K,0);replicatedOperations.assign(2*K,0);
  for (int s=0;s<2;++s) for (int c=0;c<K;++c) {
    if (sites[c]<0 || sites[c]>=E || initial[s*K+c]<=0)
      exit(__LINE__,__FILE__,"M2 cohort site/physical initial atoms invalid");
    stock[Index(s,sites[c],c)]=initial[s*K+c];
  }
  if (PIC::ThisThread==0) {
    output.open(directory/"m2_cohorts.csv");hazards.open(directory/"m2_hazards.csv");
    output << "elapsed_s,species,cohort,element,initial_atoms,first_survival,restricted_first_residence_s,surface_fraction,gas_fraction,escape_fraction,numerical_loss_fraction,tag_residual_atoms,integrated_hazard,arithmetic_operation_bound,tag_roundoff_bound_atoms\n";
    hazards << "begin_elapsed_s,dt_s,species,cohort,element,temperature_K,total_release_hazard_s-1,thermal_release_hazard_s-1,photon_release_hazard_s-1\n";
    if (!output || !hazards) exit(__LINE__,__FILE__,"cannot create cohort diagnostics; inspect output directory/disk");
  }
  FiniteReleaseBegin=BeforeRelease;FiniteReleasePacket=Packet;FiniteReleaseDeposit=Deposit;
  PIC::ParticleEvents::ParticleObserver=ParticleEvent;
  Snapshot(0);
}
void EndStep(double elapsedS) {
  elapsed=elapsedS;
  // Only deposits are rank-local. Never MPI-sum the replicated logical stock.
  std::vector<long double> combined(pending.size());
  MPI_Allreduce(pending.data(),combined.data(),static_cast<int>(pending.size()),
                MPI_LONG_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  for (std::size_t i=0;i<stock.size();++i) stock[i]+=combined[i];
  std::fill(pending.begin(),pending.end(),0);
  // Bound logical pending-stock reductions once, not once per replicated
  // rank. Each reduction has at most ranks-1 additions plus the stock update.
  for (auto& operations:replicatedOperations)
    operations+=static_cast<long long>(E)*(PIC::nTotalThreads+1);

}
void Snapshot(double elapsedS) {
  elapsed=elapsedS;
  std::vector<long double> gas(2*K,0),totalGas(2*K),totalLoss(4*K);
  std::vector<long long> operations=localOperations,totalOperations(2*K);
  for (auto* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node;node=node->nextNodeThisThread) {
    if (!node->block) continue;
    for (int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for (long int p=node->block->FirstCellParticleTable[cell];p!=-1;p=PIC::ParticleBuffer::GetNext(p)) {
        const int species=PIC::ParticleBuffer::GetI(p),s=Species(species);
        const double atoms=node->block->GetLocalParticleWeight(species)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        const double* f=Fractions(p);
        for (int c=0;c<K;++c) {gas[s*K+c]+=atoms*f[c];operations[s*K+c]+=4;}
      }
  }
  MPI_Allreduce(gas.data(),totalGas.data(),2*K,MPI_LONG_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(loss.data(),totalLoss.data(),4*K,MPI_LONG_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(operations.data(),totalOperations.data(),2*K,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread!=0) return;
  output << std::setprecision(std::numeric_limits<long double>::max_digits10);
  for (int s=0;s<2;++s) for (int c=0;c<K;++c) {
    long double surface=0;for (int e=0;e<E;++e) surface+=stock[Index(s,e,c)];
    const int j=s*K+c;const long double n=initial[j],escaped=totalLoss[j],removed=totalLoss[2*K+j];
    const long double residual=surface+totalGas[j]+escaped+removed-n;
    const long long operationBound=totalOperations[j]+replicatedOperations[j]+E+32*PIC::nTotalThreads;
    const long double nu=operationBound*std::numeric_limits<double>::epsilon();
    if (!(nu<1.0L)) exit(__LINE__,__FILE__,"M2 cohort roundoff bound exhausted; improve summation/window accounting before rerun");
    const long double bound=nu/(1-nu)*(n+surface+totalGas[j]+escaped+removed);

    output << elapsedS << ',' << (s==0?"Na":"Ar") << ',' << c << ',' << sites[c] << ',' << n << ','
      << first[j]/n << ',' << restricted[j]/n << ',' << surface/n << ',' << totalGas[j]/n << ','
      << escaped/n << ',' << removed/n << ',' << residual << ',' << exposure[j] << ','
      << operationBound << ',' << bound << '\n';
  }
  output.flush();hazards.flush();
  if (!output || !hazards) exit(__LINE__,__FILE__,"cohort diagnostics write failed; preserve output and free disk");
}
void End() {
  Exosphere::SourceProcesses::FiniteReleaseBegin=nullptr;
  Exosphere::SourceProcesses::FiniteReleasePacket=nullptr;
  Exosphere::SourceProcesses::FiniteReleaseDeposit=nullptr;
  PIC::ParticleEvents::ParticleObserver=nullptr;
  output.close();hazards.close();
}
} }

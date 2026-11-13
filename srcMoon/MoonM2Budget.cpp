#include "pic.h"
#include "MoonM2Budget.h"
#include "Moon.h"
#include "MoonSurfaceTemperature.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

namespace Moon { namespace M2Budget { namespace {
using namespace PIC::ParticleEvents;
constexpr int E=EventCount;
const char* names[E]={"emission","adsorption","reflection","photolytic_loss",
  "photolytic_product","escape","numerical_removal","reservoir_supply",
  "surface_deposit","surface_release"};
std::vector<std::array<long long,E>> counts;
std::vector<std::array<long double,E>> atoms;
std::vector<std::array<long double,3>> initial;
std::ofstream output;
bool valid=true;

void Observe(Event event,int species,double weight,int) {
  if (species<0 || species>=PIC::nTotalSpecies || !std::isfinite(weight) ||
      (event!=ReservoirSupply && weight<=0.0))
    exit(__LINE__,__FILE__,"M2 ledger received invalid species/event weight [atoms]");
  // Injection is serial; retain thread safety for ordinary mover callbacks.
#ifdef _OPENMP
#pragma omp critical (MoonM2EventLedger)
#endif
  {++counts[species][event];atoms[species][event]+=weight;}
}

Event ClassifyExit(const double* x,const double* v) {
  // Reuse the already approved I09 definition, not ultimate three-body
  // escape: epsilon=|v_LSO+Omega x r|^2/2-GM_Moon/r. Bound finite-box exits
  // are numerical-domain losses and remain distinct from positive-energy
  // lunar escapes. Invalid state is fatal, never silently called escape.
  double inertial[3],r2=0.0,v2=0.0;
  const double* w=Moon::RotationVector_SO_FROZEN;
  inertial[0]=v[0]+w[1]*x[2]-w[2]*x[1];
  inertial[1]=v[1]+w[2]*x[0]-w[0]*x[2];
  inertial[2]=v[2]+w[0]*x[1]-w[1]*x[0];
  for (int i=0;i<3;++i) {r2+=x[i]*x[i];v2+=inertial[i]*inertial[i];}
  if (!(r2>0.0) || !std::isfinite(r2+v2))
    exit(__LINE__,__FILE__,"M2 domain-loss classification has invalid LSO state");
  return 0.5*v2-GravityConstant*_MASS_(_MOON_)/std::sqrt(r2)>0.0 ? Escape:NumericalRemoval;
}

std::array<long double,3> Inventory(int species) {
  // [model count, gas atoms, effective surface atoms]. The logical surface
  // stock is replicated and counted ONCE. Unexchanged local adsorption minus
  // release flux is added per rank, so partial-step checkpoints remain exact
  // without forcing an extra exchange or changing source/thermal scheduling.
  long double local[3]={},global[3];
  for (auto* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node!=nullptr;node=node->nextNodeThisThread) {
    if (node->block==nullptr) continue;
    for (int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for (long int p=node->block->FirstCellParticleTable[cell];p!=-1;
           p=PIC::ParticleBuffer::GetNext(p)) if (PIC::ParticleBuffer::GetI(p)==species) {
        double weight=node->block->GetLocalParticleWeight(species);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        weight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
#endif
        local[0]+=1;local[1]+=weight;
      }
  }
  for (int e=0;e<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;++e) {
    if (PIC::ThisThread==0) local[2]+=Moon::Planet->SurfaceElementPopulation[species][e];
    local[2]+=static_cast<long double>(Moon::Planet->SurfaceElementAdsorptionFluxDOWN[species][e])-
        Moon::Planet->SurfaceElementDesorptionFluxUP[species][e];
  }
  MPI_Allreduce(local,global,3,MPI_LONG_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  return {global[0],global[1],global[2]};
}
} // namespace

void Begin(const std::filesystem::path& csv) {
  if (PIC::ParticleEvents::Observer!=nullptr || PIC::ParticleEvents::DomainExitClassifier!=nullptr ||
      _NA_SPEC_<0 || _AR_SPEC_<0 || PIC::nTotalSpecies!=2 ||
      !Moon::SurfaceTemperature::UseProductionFluxBalancedInventory())
    exit(__LINE__,__FILE__,"M2 budget needs the dedicated finite-inventory Na/Ar configuration and unclaimed diagnostic observers");
  counts.assign(PIC::nTotalSpecies,{});atoms.assign(PIC::nTotalSpecies,{});
  initial.resize(PIC::nTotalSpecies);valid=true;
  for (int s=0;s<PIC::nTotalSpecies;++s) initial[s]=Inventory(s);
  if (PIC::ThisThread==0) {
    output.open(csv);
    output << std::setprecision(std::numeric_limits<long double>::max_digits10)
      << "completed_steps,species,mass_kg,initial_count,gas_count,initial_gas_atoms,gas_atoms,initial_surface_atoms,surface_atoms,count_residual,gas_residual_atoms,surface_residual_atoms,combined_residual_atoms,roundoff_bound_atoms,combined_residual_kg,arithmetic_epsilon";
    for (const char* name:names) output << ',' << name << "_count," << name << "_atoms";
    output << '\n';
    if (!output) exit(__LINE__,__FILE__,"cannot open M2 budget CSV");
  }
  PIC::ParticleEvents::Observer=Observe;
  PIC::ParticleEvents::DomainExitClassifier=ClassifyExit;
  Snapshot(0);
}

void Snapshot(long long steps) {
  for (int s=0;s<PIC::nTotalSpecies;++s) {
    auto state=Inventory(s);
    std::array<long long,E> c{};std::array<long double,E> a{};
    MPI_Allreduce(counts[s].data(),c.data(),E,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(atoms[s].data(),a.data(),E,MPI_LONG_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    const long long countResidual=static_cast<long long>(state[0]-initial[s][0])-
      c[Emission]+c[Adsorption]+c[PhotolyticLoss]-c[PhotolyticProduct]+c[Escape]+c[NumericalRemoval];
    const long double gasResidual=state[1]-initial[s][1]-a[Emission]+
      a[Adsorption]+a[PhotolyticLoss]-a[PhotolyticProduct]+a[Escape]+a[NumericalRemoval];
    const long double surfaceResidual=state[2]-initial[s][2]-a[SurfaceDeposit]+
      a[SurfaceRelease]-a[ReservoirSupply];
    long double magnitude=std::fabs(state[1])+std::fabs(state[2])+
      std::fabs(initial[s][1])+std::fabs(initial[s][2]);
    long double operations=state[0]+initial[s][0]+32*PIC::nTotalThreads+
      4.0L*(steps+1)*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
    for (int e=0;e<E;++e) {magnitude+=std::fabs(a[e]);operations+=4.0L*c[e];}
    // Surface arrays and local weight products use double, even though the
    // counters/list reductions use long double. Bound BOTH sources using a
    // conservative double forward-error gamma_n; do not fit a physical error
    // allowance to observed residuals. The independent Python oracle also
    // enforces the separate frozen <0.5% dominant-source criterion.
    const long double nu=operations*std::numeric_limits<double>::epsilon();
    if (!(nu<1.0L)) exit(__LINE__,__FILE__,"M2 ledger roundoff bound exhausted");
    const long double bound=nu/(1.0L-nu)*magnitude;
    valid=valid && countResidual==0 && std::fabs(gasResidual)<=bound &&
      std::fabs(surfaceResidual)<=bound && a[PhotolyticLoss]==0 && a[PhotolyticProduct]==0;
    if (PIC::ThisThread==0) {
      const double mass=PIC::MolecularData::GetMass(s);
      output << steps << ',' << (s==_NA_SPEC_?"Na":"Ar") << ',' << mass << ','
        << initial[s][0] << ',' << state[0] << ',' << initial[s][1] << ',' << state[1]
        << ',' << initial[s][2] << ',' << state[2] << ',' << countResidual << ','
        << gasResidual << ',' << surfaceResidual << ',' << gasResidual+surfaceResidual
        << ',' << bound << ',' << (gasResidual+surfaceResidual)*mass << ','
        << std::numeric_limits<long double>::epsilon();
      for (int e=0;e<E;++e) output << ',' << c[e] << ',' << a[e];
      output << '\n';output.flush();
      if (!output) exit(__LINE__,__FILE__,"M2 budget write failed; check disk space");
    }
  }
}
bool Valid() {return valid;}
void End() {PIC::ParticleEvents::Observer=nullptr;
  PIC::ParticleEvents::DomainExitClassifier=nullptr;output.close();}
} }

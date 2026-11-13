#include "pic.h"
#include "Exosphere.h"
#include "FiniteSurfaceRelease.h"
#include <algorithm>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

int Exosphere::SourceProcesses::FiniteReservoirElement=-1;

long int Exosphere::SourceProcesses::InjectFiniteSurfaceInventory(
    int spec,int boundaryType,void* boundary) {
  // This path is enabled only by an application's explicit finite-inventory
  // policy. Ordinary legacy/prescribed reservoirs never enter this routine.
  // Each MPI rank has the same synchronized logical reservoir; generating
  // independent candidates per rank would lose/duplicate physical atoms.
  // Root draws a candidate using the NORMAL source generator, broadcasts it,
  // and exactly one final terrain/AMR owner allocates the weighted particle.
  if (boundaryType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_)
    exit(__LINE__,__FILE__,"finite inventory requires the logical spherical reservoir");
#ifdef _OPENMP
  if (omp_get_max_threads()!=1)
    exit(__LINE__,__FILE__,"finite reservoir source dispatch currently requires OMP_NUM_THREADS=1; use MPI ranks for parallelism until hybrid source scheduling is qualified");
#endif
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  exit(__LINE__,__FILE__,"finite inventory needs individual statistical weights for exact partial packets");
#endif
#if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  exit(__LINE__,__FILE__,"finite inventory requires a synchronized global species timestep");
#endif
  auto* sphere=static_cast<cInternalSphericalData*>(boundary);
  const double dt=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
  const double weight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
  if (!(dt>0.0) || !(weight>0.0) || !std::isfinite(dt+weight))
    exit(__LINE__,__FILE__,"finite inventory has invalid seconds/particle weight");
  double *centre=nullptr,radius=0.0;
  sphere->GetSphereGeometricalParameters(centre,radius);
  long int injected=0;
  for (int element=0;element<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;++element) {
    const double population=sphere->SurfaceElementPopulation[spec][element];
    if (!std::isfinite(population) || population<0.0)
      exit(__LINE__,__FILE__,"finite reservoir population is invalid [atoms]");
    if (population==0.0) continue;
    std::vector<double> hazards(1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_,0.0);
#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
    hazards[_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]=
        ThermalDesorption::GetSurfaceElementProductionRate(spec,element,sphere)/population;
#endif
#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
    hazards[_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]=
        PhotonStimulatedDesorption::GetSurfaceElementProductionRate(spec,element,sphere)/population;
#endif
#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
    // The existing sputtering law is population independent. Treating it as
    // an exponential hazard would invent physics. Preserve the explicit gap.
    if (SolarWindSputtering::GetSurfaceElementProductionRate(spec,element,sphere)>0.0)
      exit(__LINE__,__FILE__,"finite-reservoir sputtering needs a separately specified population-independent depletion law; disable sputtering or qualify that coupling");
#endif
    FiniteRelease release;
    try {release=IntegrateFiniteRelease(population,hazards,dt);}
    catch (const std::exception& error) {exit(__LINE__,__FILE__,error.what());}
    for (int source=0;source<=_EXOSPHERE__SOURCE_MAX_ID_VALUE_;++source) {
      const double atoms=release.atoms[source];
      if (!(atoms>0.0)) continue;
      if (!Source_DeplitSurfaceSpeciesAbundance_Flag[source])
        exit(__LINE__,__FILE__,"finite release source must debit its surface inventory");
      const double packetCount=std::ceil(atoms/weight);
      if (!std::isfinite(packetCount) || packetCount>std::numeric_limits<long int>::max())
        exit(__LINE__,__FILE__,"finite release packet count exceeds addressable range; increase particle weight or reduce inventory/resolution and rerun");
      for (long int packet=0;packet<static_cast<long int>(packetCount);++packet) {
        const double packetAtoms=std::min(weight,atoms-packet*weight);
        if (!(packetAtoms>0.0)) continue;
        double state[13]={}; // SO x, body x, SO v, body v, departure seconds.
        std::vector<PIC::ParticleBuffer::byte> data(PIC::ParticleBuffer::ParticleDataLength,0);
        int generated=1;
        if (PIC::ThisThread==0) {
          FiniteReservoirElement=element;
          auto* node=static_cast<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*>(nullptr);
          bool ok=false;
#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
          if (source==_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_)
            ok=ThermalDesorption::GenerateParticleProperties(spec,data.data(),state,state+3,state+6,state+9,centre,radius,node,boundaryType,boundary);
#endif
#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
          if (source==_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_)
            ok=PhotonStimulatedDesorption::GenerateParticleProperties(spec,data.data(),state,state+3,state+6,state+9,centre,radius,node,boundaryType,boundary);
#endif
          FiniteReservoirElement=-1;
          generated=ok ? 1:0;
          // Conditional first-order departure-time CDF. In particular a very
          // short Ar residence time releases near the START, not uniformly
          // throughout the particle step. k*dt may overflow safely to +inf.
          state[12]=-std::log1p(-rnd()*release.probability)/release.hazardS;
        }
        MPI_Bcast(&generated,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
        if (!generated) exit(__LINE__,__FILE__,"normal source generator rejected a forced finite-reservoir candidate");
        MPI_Bcast(state,13,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Bcast(data.data(),data.size(),MPI_BYTE,0,MPI_GLOBAL_COMMUNICATOR);
        PIC::Mover::cMovingInjectionEpochScope injectionEpoch(state[12]);
        auto* node=PIC::Mesh::mesh->findTreeNode(state);
        bool owned=node!=nullptr && node->Thread==PIC::ThisThread && node->block!=nullptr;
        if (PostProcessSurfaceInjection!=nullptr)
          owned=PostProcessSurfaceInjection(spec,source,state,state+3,state+6,state+9,node);
        int owner=owned?1:0,owners=0;
        MPI_Allreduce(&owner,&owners,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        if (owners!=1) exit(__LINE__,__FILE__,"finite release needs exactly one final terrain/AMR owner; inspect surface geometry and domain");
        if (!owned) continue;
        const long int particle=PIC::ParticleBuffer::GetNewParticle();
        PIC::ParticleBuffer::CloneParticle(
            PIC::ParticleBuffer::GetParticleDataPointer(particle),data.data());
        PIC::ParticleBuffer::SetI(spec,particle);
        PIC::ParticleBuffer::SetX(state,particle);
        PIC::ParticleBuffer::SetV(state+6,particle);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(packetAtoms/weight,particle);
        Sampling::SetParticleSourceID(source,PIC::ParticleBuffer::GetParticleDataPointer(particle));
        Sampling::SetParicleOriginSurfaceElementNumber(element,PIC::ParticleBuffer::GetParticleDataPointer(particle));
        const double representedAtoms=weight*
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle);
        PIC::ParticleEvents::Record(PIC::ParticleEvents::Emission,spec,representedAtoms,source);
        PIC::ParticleEvents::Record(PIC::ParticleEvents::SurfaceRelease,spec,representedAtoms,element);
        Sampling::CalculatedSourceRate[spec][source]+=representedAtoms/dt;
        PIC::BC::nInjectedParticles[spec]++;
        PIC::BC::ParticleProductionRate[spec]+=representedAtoms/dt;
        if (sphere->SampleSpeciesSurfaceInjectionFlux!=nullptr)
          sphere->SampleSpeciesSurfaceInjectionFlux[spec][element]+=representedAtoms;
        PIC::ParticleBuffer::SetParticleAllocated(particle);
        ++injected;
        _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(particle,
            std::max(0.0,dt-state[12]),node,true);
      }
    }
    // Store the analytic continuous remainder directly on every replicated
    // rank. Do NOT additionally queue these releases in DesorptionFluxUP:
    // doing so would debit twice. Subtracting independently MPI-summed packet
    // totals from a fully depleted double stock can produce a tiny negative
    // population through cancellation. This representation needs no empirical
    // clamp or fictitious replenishment; the ledger independently bounds the
    // analytic stock versus actual emitted-weight summation error. Incoming
    // adsorption stays pending and is added by the ordinary end-step exchange.
    sphere->SurfaceElementPopulation[spec][element]=release.remainingAtoms;
  }
  return injected;
}

#include "CometChemistry.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "pic.h"
#include "CometChemistry.generated.h"

bool Comet::Chemistry::Enabled=false;
double Comet::Chemistry::HeliocentricDistanceAU=-1.0;

const ChemicalSolver::cReactionNetwork& Comet::Chemistry::GetNetwork() {
  static const ChemicalSolver::cReactionNetwork network(
      Generated::Terms,Generated::TermCount,
      Generated::Reactions,Generated::ReactionCount,
      PIC::nTotalSpecies);
  return network;
}

bool Comet::Chemistry::Validate(std::string *error) {
  if (Enabled && !(HeliocentricDistanceAU>0.0)) {
    if (error!=NULL)
      *error="enabled comet chemistry requires a positive heliocentric distance";
    return false;
  }

  if (!GetNetwork().Validate(error)) return false;

  // The particle executor below represents discrete decay events. General
  // mass-action reactions remain useful to continuum solvers, but enabling
  // them here would assign a density-dependent cell rate to one particle
  // without the required reactant-pair selection.
  if (Enabled) for (std::size_t reaction=0;
      reaction<GetNetwork().GetReactionCount();reaction++) {
    if (!GetNetwork().IsUnitUnimolecular(reaction)) {
      if (error!=NULL)
        *error="the PIC comet chemistry executor accepts only unit unimolecular reactions";
      return false;
    }

    const ChemicalSolver::cReaction &entry=GetNetwork().GetReaction(reaction);
    double reactantMass=0.0,productMass=0.0;
    for (std::size_t reactant=0;reactant<entry.ReactantCount;reactant++) {
      const ChemicalSolver::cStoichiometricTerm &term=
          GetNetwork().GetReactantTerm(reaction,reactant);
      reactantMass+=term.Coefficient*PIC::MolecularData::GetMass(term.Species);
    }
    for (std::size_t product=0;product<entry.ProductCount;product++) {
      const ChemicalSolver::cStoichiometricTerm &term=
          GetNetwork().GetProductTerm(reaction,product);
      const double coefficient=term.Coefficient;
      if (coefficient<1.0 || std::fabs(coefficient-std::floor(coefficient+0.5))>
          1.0e-12) {
        if (error!=NULL)
          *error="the PIC comet chemistry executor requires positive integer product coefficients";
        return false;
      }
      productMass+=coefficient*PIC::MolecularData::GetMass(term.Species);
    }

    // Product-free ionization entries are neutral sinks. Every neutral
    // dissociation channel must conserve the compiled species masses within
    // the precision of the shared species catalog. Some legacy entries use
    // rounded SI atomic masses while molecules use nominal integer amu, so an
    // exact floating-point equality would reject valid stoichiometry. The
    // 0.1% relative tolerance remains far below a missing-nucleon error and
    // still catches transcription errors such as C2+H2 versus C2H+H2.
    if (entry.ProductCount!=0 &&
        std::fabs(productMass-reactantMass)>
            5.0e-4*(productMass+reactantMass)) {
      if (error!=NULL)
        *error="a neutral comet chemistry channel does not conserve species mass";
      return false;
    }
  }

  return true;
}

namespace {

// Reinsert an unchanged particle into the cell list being rebuilt by the
// AMPS photochemistry driver. The driver has detached the old list before it
// invokes the application callback, so every surviving parent must be linked
// explicitly.
void KeepParticle(long int ptr,long int &firstParticleCell) {
  PIC::ParticleBuffer::SetNext(firstParticleCell,ptr);
  PIC::ParticleBuffer::SetPrev(-1,ptr);
  if (firstParticleCell!=-1)
    PIC::ParticleBuffer::SetPrev(ptr,firstParticleCell);
  firstParticleCell=ptr;
}

double StandardNormal() {
  // Box--Muller converts two open-interval uniform samples into an isotropic
  // normal component. Avoid log(0) without changing any representable
  // nonzero sample.
  double u1;
  do u1=rnd(); while (!(u1>0.0));
  const double u2=rnd();
  return std::sqrt(-2.0*std::log(u1))*std::cos(2.0*Pi*u2);
}

void GenerateProductVelocities(const std::vector<int> &species,
    const double parentVelocity[3],double excessEnergyJ,
    std::vector<double> &velocity) {
  const std::size_t count=species.size();
  velocity.assign(3*count,0.0);
  if (count==0) return;

  if (count==1) {
    // A single tracked product cannot carry internal relative kinetic energy.
    // The current Rubin network has no such neutral channel; retain this
    // branch for zero-excess generic reactions and reject unsupported input.
    if (excessEnergyJ>0.0) {
      std::fprintf(stderr,
          "one-product comet reaction has nonzero excess energy\n");
      std::abort();
    }
    for (int idim=0;idim<3;idim++) velocity[idim]=parentVelocity[idim];
    return;
  }

  // Draw an orientation in the 3N-dimensional product-velocity space. Then
  // remove its mass-weighted center-of-mass motion and rescale the remaining
  // velocities to the channel's exact excess energy. This produces no
  // preferred direction and enforces both sum(m_i u_i)=0 and
  // sum(m_i |u_i|^2/2)=E for two- and many-body channels.
  double totalMass=0.0,centerOfMass[3]={0.0,0.0,0.0};
  for (std::size_t product=0;product<count;product++) {
    const double mass=PIC::MolecularData::GetMass(species[product]);
    totalMass+=mass;
    for (int idim=0;idim<3;idim++) {
      const double value=StandardNormal();
      velocity[3*product+idim]=value;
      centerOfMass[idim]+=mass*value;
    }
  }
  for (int idim=0;idim<3;idim++) centerOfMass[idim]/=totalMass;

  double sampledEnergy=0.0;
  for (std::size_t product=0;product<count;product++) {
    const double mass=PIC::MolecularData::GetMass(species[product]);
    for (int idim=0;idim<3;idim++)
      velocity[3*product+idim]-=centerOfMass[idim];
    const double u2=velocity[3*product]*velocity[3*product]+
        velocity[3*product+1]*velocity[3*product+1]+
        velocity[3*product+2]*velocity[3*product+2];
    sampledEnergy+=0.5*mass*u2;
  }

  const double scale=(excessEnergyJ>0.0) ?
      std::sqrt(excessEnergyJ/sampledEnergy) : 0.0;
  for (std::size_t product=0;product<count;product++)
    for (int idim=0;idim<3;idim++)
      velocity[3*product+idim]=parentVelocity[idim]+
          scale*velocity[3*product+idim];
}

} // namespace

void Comet::Chemistry::ProcessParticle(long int ptr,long int &firstParticleCell,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (!Enabled) {
    KeepParticle(ptr,firstParticleCell);
    return;
  }

  PIC::ParticleBuffer::byte *parentData=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  const int parentSpecies=PIC::ParticleBuffer::GetI(parentData);
  const ChemicalSolver::cReactionNetwork &network=GetNetwork();

  // Build the total decay rate for this parent. Selection by a single random
  // channel below is essential: independently sampling expected daughter
  // yields could combine products from mutually exclusive Table 4 branches.
  std::vector<std::size_t> channels;
  std::vector<double> rates;
  double totalRate=0.0;
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
    if (!network.IsUnitUnimolecular(reaction) ||
        network.GetReactantTerm(reaction,0).Species!=parentSpecies) continue;
    const double rate=network.GetRateCoefficient(
        reaction,0.0,HeliocentricDistanceAU);
    if (rate>0.0) {
      channels.push_back(reaction);
      rates.push_back(rate);
      totalRate+=rate;
    }
  }

  if (!(totalRate>0.0)) {
    KeepParticle(ptr,firstParticleCell);
    return;
  }

  const double timeStep=node->block->GetLocalTimeStep(parentSpecies);
  if (!(timeStep>0.0)) {
    std::fprintf(stderr,"invalid comet chemistry timestep for species %d\n",
        parentSpecies);
    std::abort();
  }
  if (rnd()>=1.0-std::exp(-totalRate*timeStep)) {
    KeepParticle(ptr,firstParticleCell);
    return;
  }

  double selector=rnd()*totalRate;
  std::size_t selected=channels.back();
  for (std::size_t i=0;i<channels.size();i++) {
    selector-=rates[i];
    if (selector<=0.0) {
      selected=channels[i];
      break;
    }
  }
  const ChemicalSolver::cReaction &reaction=network.GetReaction(selected);

  // Expand integer stoichiometric coefficients into concrete daughter
  // particles. A zero-product channel is a photoionization loss from this
  // neutral-only network and therefore simply removes the parent.
  std::vector<int> productSpecies;
  for (std::size_t term=0;term<reaction.ProductCount;term++) {
    const ChemicalSolver::cStoichiometricTerm &product=
        network.GetProductTerm(selected,term);
    const int multiplicity=static_cast<int>(std::floor(product.Coefficient+0.5));
    for (int copy=0;copy<multiplicity;copy++)
      productSpecies.push_back(product.Species);
  }

  double xParent[3],vParent[3];
  PIC::ParticleBuffer::GetX(xParent,parentData);
  PIC::ParticleBuffer::GetV(vParent,parentData);
  std::vector<double> productVelocity;
  GenerateProductVelocities(productSpecies,vParent,reaction.ExcessEnergyJ,
      productVelocity);

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  std::fprintf(stderr,
      "comet chemistry requires ParticleWeightCorrectionMode=on\n");
  std::abort();
#else
  const double parentPhysicalWeight=
      node->block->GetLocalParticleWeight(parentSpecies)*
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(parentData);

  // Clone auxiliary fields from the parent, then replace the state and scale
  // each daughter's individual correction so one reaction event carries the
  // same physical statistical weight for every product species.
  char daughterTemplate[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated(
      reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate));
  PIC::ParticleBuffer::CloneParticle(
      reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate),parentData);

  for (std::size_t product=0;product<productSpecies.size();product++) {
    const int species=productSpecies[product];
    const double baseWeight=node->block->GetLocalParticleWeight(species);
    if (!(baseWeight>0.0)) {
      std::fprintf(stderr,
          "invalid base particle weight for comet chemistry species %d\n",
          species);
      std::abort();
    }

    PIC::ParticleBuffer::byte *data=
        reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate);
    PIC::ParticleBuffer::SetX(xParent,data);
    PIC::ParticleBuffer::SetV(productVelocity.data()+3*product,data);
    PIC::ParticleBuffer::SetI(species,data);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
        parentPhysicalWeight/baseWeight,data);

#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    PIC::ParticleTracker::InitParticleID(daughterTemplate);
    PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(
        xParent,productVelocity.data()+3*product,species,daughterTemplate,
        static_cast<void*>(node));
#endif

    const long int daughter=PIC::ParticleBuffer::GetNewParticle(firstParticleCell);
    PIC::ParticleBuffer::CloneParticle(
        PIC::ParticleBuffer::GetParticleDataPointer(daughter),data);
  }
#endif

  PIC::ParticleBuffer::DeleteParticle(ptr);
}

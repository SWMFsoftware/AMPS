#include "CometChemistry.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

#include "pic.h"
#include "CometChemistry.generated.h"
#include "CometDiagnostics.h"
#include "RadialSimulation.h"

bool Comet::Chemistry::Enabled=false;
double Comet::Chemistry::HeliocentricDistanceAU=-1.0;

namespace {
long int ChemistryProcessedFlagOffset=-1;
bool LocalChemistryAllocationPressure=false;

// Per-parent reaction dispatch is immutable during a chemistry stage.  The
// generic network supplies applicable reaction indices in generated-table
// order; this application cache adds the rates evaluated at the selected
// heliocentric distance.  Constructing it outside ProcessParticle removes a
// full-network scan and two vector allocations from every model particle.
struct cParentReactionRateCache {
  std::vector<std::size_t> Channels;
  std::vector<double> Rates;
  double TotalRate=0.0;
};

std::vector<cParentReactionRateCache> ParentReactionRates;
double CachedHeliocentricDistanceAU=
    std::numeric_limits<double>::quiet_NaN();
bool ParentReactionRatesReady=false;

bool RefreshParentReactionRates(std::string *error) {
  ParentReactionRates.clear();
  ParentReactionRates.resize(PIC::nTotalSpecies);
  ParentReactionRatesReady=false;

  if (!Comet::Chemistry::Enabled) {
    CachedHeliocentricDistanceAU=Comet::Chemistry::HeliocentricDistanceAU;
    ParentReactionRatesReady=true;
    return true;
  }
  if (!(Comet::Chemistry::HeliocentricDistanceAU>0.0) ||
      !std::isfinite(Comet::Chemistry::HeliocentricDistanceAU)) {
    if (error!=NULL)
      *error="enabled comet chemistry requires a positive finite heliocentric distance";
    return false;
  }

  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    cParentReactionRateCache &cache=ParentReactionRates[species];
    const std::size_t channelCount=
        network.GetUnitUnimolecularReactionCount(species);
    cache.Channels.reserve(channelCount);
    cache.Rates.reserve(channelCount);

    // Retain input order and the original left-to-right rate accumulation.
    // Besides preserving mutually exclusive channel semantics, this keeps the
    // floating-point selection boundaries and random-number stream unchanged.
    for (std::size_t channel=0;channel<channelCount;channel++) {
      const std::size_t reaction=
          network.GetUnitUnimolecularReactionIndex(species,channel);
      if (reaction>=network.GetReactionCount()) {
        if (error!=NULL)
          *error="invalid indexed comet chemistry reaction";
        return false;
      }
      const double rate=network.GetRateCoefficient(
          reaction,0.0,Comet::Chemistry::HeliocentricDistanceAU);
      if (!std::isfinite(rate) || rate<0.0) {
        if (error!=NULL)
          *error="comet chemistry produced an invalid cached reaction rate";
        return false;
      }
      if (rate>0.0) {
        cache.Channels.push_back(reaction);
        cache.Rates.push_back(rate);
        cache.TotalRate+=rate;
      }
    }
    if (!std::isfinite(cache.TotalRate) || cache.TotalRate<0.0) {
      if (error!=NULL)
        *error="comet chemistry produced an invalid cached total rate";
      return false;
    }
  }

  CachedHeliocentricDistanceAU=Comet::Chemistry::HeliocentricDistanceAU;
  ParentReactionRatesReady=true;
  return true;
}

unsigned char& ChemistryProcessedFlag(
    PIC::ParticleBuffer::byte *particleData) {
  if (ChemistryProcessedFlagOffset<0) {
    std::fprintf(stderr,
        "comet chemistry processing flag was not reserved before particle-buffer initialization\n");
    std::abort();
  }
  return *reinterpret_cast<unsigned char*>(
      particleData+ChemistryProcessedFlagOffset);
}

void MarkChemistryProcessed(PIC::ParticleBuffer::byte *particleData) {
  ChemistryProcessedFlag(particleData)=1;
}
}

void Comet::Chemistry::RequestParticleData() {
  if (!Enabled || ChemistryProcessedFlagOffset>=0) return;
  PIC::ParticleBuffer::RequestDataStorage(ChemistryProcessedFlagOffset,
      sizeof(unsigned char));
}

const ChemicalSolver::cReactionNetwork& Comet::Chemistry::GetNetwork() {
  static const ChemicalSolver::cReactionNetwork network(
      Generated::Terms,Generated::TermCount,
      Generated::Reactions,Generated::ReactionCount,
      PIC::nTotalSpecies);
  return network;
}

bool Comet::Chemistry::Validate(std::string *error) {
  if (Enabled && (!(HeliocentricDistanceAU>0.0) ||
      !std::isfinite(HeliocentricDistanceAU))) {
    if (error!=NULL)
      *error="enabled comet chemistry requires a positive finite heliocentric distance";
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

  // Build application runtime rates only after the generic network and the
  // particle executor's stronger mass/multiplicity contract have succeeded.
  // A subsequent input/state change is detected at BeginBatchedModel, which is
  // a serial stage boundary rather than a particle-parallel hot path.
  return RefreshParentReactionRates(error);
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

double GetResidualMultiplicityMean(double sourceWeight,double sourceTimeStep,
    double targetWeight,double targetTimeStep) {
  // Species-local pseudo-time requires conversion of physical rate, not just
  // physical weight.  One source particle carries W_source/dt_source through
  // an outer iteration, while a target particle carries W_target/dt_target.
  // Residual resampling with
  //
  //   mu=(W_source/dt_source)/(W_target/dt_target)
  //
  // is therefore the species-changing analogue of the mover's block-face
  // conversion.  It is Russian roulette for mu<1 and cloning for mu>1.
  if (!(sourceWeight>0.0) || !std::isfinite(sourceWeight) ||
      !(sourceTimeStep>0.0) || !std::isfinite(sourceTimeStep) ||
      !(targetWeight>0.0) || !std::isfinite(targetWeight)) {
    std::fprintf(stderr,"invalid comet chemistry source weight or timestep\n");
    std::abort();
  }
  if (!(targetTimeStep>0.0) || !std::isfinite(targetTimeStep)) {
    std::fprintf(stderr,"invalid comet chemistry target weight or timestep\n");
    std::abort();
  }

  const double mean=sourceWeight*targetTimeStep/
      (targetWeight*sourceTimeStep);
  if (!(mean>0.0) || !std::isfinite(mean) ||
      mean>static_cast<double>(std::numeric_limits<long int>::max()-1)) {
    std::fprintf(stderr,"invalid comet chemistry resampling multiplicity\n");
    std::abort();
  }

  return mean;
}

long int DrawResidualMultiplicity(double sourceWeight,double sourceTimeStep,
    double targetWeight,double targetTimeStep) {
  const double mean=GetResidualMultiplicityMean(sourceWeight,sourceTimeStep,
      targetWeight,targetTimeStep);

  const double integerPart=std::floor(mean);
  long int result=static_cast<long int>(integerPart);
  if (rnd()<mean-integerPart) result++;
  return result;
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

struct cChemistryEnsembleMoments {
  std::vector<double> ExpectedNumber,ActualNumber,NumberVariance;
  double ExpectedMass=0.0,ActualMass=0.0,MassVariance=0.0;
  double ExpectedMomentum[3]={0.0,0.0,0.0};
  double ActualMomentum[3]={0.0,0.0,0.0};
  double MomentumVariance[3]={0.0,0.0,0.0};
  double ExpectedEnergy=0.0,ActualEnergy=0.0,EnergyVariance=0.0;

  cChemistryEnsembleMoments() :
      ExpectedNumber(PIC::nTotalSpecies,0.0),
      ActualNumber(PIC::nTotalSpecies,0.0),
      NumberVariance(PIC::nTotalSpecies,0.0) {}
};

bool CheckEnsembleMoment(const char *reactionName,const char *momentName,
    double expected,double actual,double variance) {
  if (!(variance>=0.0) || !std::isfinite(variance) ||
      !std::isfinite(expected) || !std::isfinite(actual)) {
    if (PIC::ThisThread==0)
      std::printf("  FAIL chemistry reaction=%s moment=%s invalid state "
          "expected=%.9e actual=%.9e variance=%.9e\n",reactionName,
          momentName,expected,actual,variance);
    return false;
  }

  const double difference=std::fabs(actual-expected);
  const double roundoff=1024.0*std::numeric_limits<double>::epsilon()*
      std::max({std::fabs(expected),std::fabs(actual),
          std::numeric_limits<double>::min()});
  const double allowed=std::max({0.02*std::fabs(expected),
      5.0*std::sqrt(variance),roundoff});
  if (difference<=allowed) return true;

  if (PIC::ThisThread==0) {
    const double z=(variance>0.0) ? difference/std::sqrt(variance) : INFINITY;
    const double relative=(expected!=0.0) ? difference/std::fabs(expected) :
        INFINITY;
    std::printf("  FAIL chemistry reaction=%s moment=%s difference=%.9e "
        "allowed=%.9e relative=%.6g z=%.6g\n",reactionName,momentName,
        difference,allowed,relative,z);
  }
  return false;
}

} // namespace

bool Comet::Chemistry::ValidateEnsemble(long int eventsPerReaction) {
  if (!Enabled) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "chemistry ensemble validation requires ChemistryMode=on\n");
    return false;
  }
  if (eventsPerReaction<=0) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,"chemistry validation event count must be positive\n");
    return false;
  }

  // RadialSimulation::Initialize() allocates all local data blocks, but the
  // core builds DomainBlockDecomposition::BlockTable lazily at the start of a
  // normal time step.  The focused validator exits before that time-step path,
  // so explicitly synchronize the table with the already distributed mesh.
  // This uses the same core API as the local-time-step transport validator and
  // guarantees that the weight/timestep samples below come from an actual
  // locally owned production block on every MPI rank.
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  if (PIC::DomainBlockDecomposition::nLocalBlocks<=0 ||
      PIC::DomainBlockDecomposition::BlockTable[0]==NULL ||
      PIC::DomainBlockDecomposition::BlockTable[0]->block==NULL) {
    std::fprintf(stderr,
        "rank %d has no local block for chemistry ensemble validation\n",
        PIC::ThisThread);
    std::abort();
  }

  // Every rank uses one of its actual configured blocks.  Consequently the
  // ensemble samples the same species-dependent W/dt carrier rates and
  // block-local weights/timesteps used by the production Halley run at several
  // radial locations.  Each reaction's requested event count is divided
  // across ranks, not duplicated per rank.
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
      PIC::DomainBlockDecomposition::BlockTable[0];
  const long int localEvents=eventsPerReaction/PIC::nTotalThreads+
      ((PIC::ThisThread<eventsPerReaction%PIC::nTotalThreads) ? 1 : 0);
  const double parentVelocity[3]={431.0,-137.0,59.0};
  const double parentSpeed2=parentVelocity[0]*parentVelocity[0]+
      parentVelocity[1]*parentVelocity[1]+
      parentVelocity[2]*parentVelocity[2];
  const ChemicalSolver::cReactionNetwork &network=GetNetwork();
  bool allPass=true;
  int testedReactionCount=0,sinkReactionCount=0,unequalReactionCount=0;

  if (PIC::ThisThread==0)
    std::printf("Chemistry ensemble validation: events/reaction=%ld ranks=%d "
        "relative-limit=0.02 z-limit=5.0\n",eventsPerReaction,
        PIC::nTotalThreads);

  for (std::size_t reactionIndex=0;
       reactionIndex<network.GetReactionCount();reactionIndex++) {
    const ChemicalSolver::cReaction &reaction=
        network.GetReaction(reactionIndex);
    if (reaction.ProductCount==0) {
      sinkReactionCount++;
      continue;
    }

    const int parentSpecies=
        network.GetReactantTerm(reactionIndex,0).Species;
    std::vector<int> products;
    for (std::size_t term=0;term<reaction.ProductCount;term++) {
      const ChemicalSolver::cStoichiometricTerm &product=
          network.GetProductTerm(reactionIndex,term);
      const int coefficient=
          static_cast<int>(std::floor(product.Coefficient+0.5));
      for (int copy=0;copy<coefficient;copy++)
        products.push_back(product.Species);
    }

    const double parentWeight=
        node->block->GetLocalParticleWeight(parentSpecies);
    const double parentTimeStep=node->block->GetLocalTimeStep(parentSpecies);
    if (!(parentWeight>0.0) || !(parentTimeStep>0.0)) {
      std::fprintf(stderr,
          "invalid parent weight/timestep in chemistry validation\n");
      std::abort();
    }
    const double sourceRateWeight=parentWeight/parentTimeStep;
    const double parentMass=PIC::MolecularData::GetMass(parentSpecies);
    double productMass=0.0;
    for (int species:products)
      productMass+=PIC::MolecularData::GetMass(species);

    // Apply the same mass tolerance as startup validation.  This separately
    // exposes catalog rounding while the statistical checks below measure the
    // residual-resampling implementation itself.
    bool reactionPass=std::fabs(productMass-parentMass)<=
        5.0e-4*(productMass+parentMass);
    const double catalogMassRelative=std::fabs(productMass-parentMass)/
        std::max(parentMass,std::numeric_limits<double>::min());

    cChemistryEnsembleMoments local;
    double localMaximumMomentumError=0.0;
    double localMaximumEnergyError=0.0;
    double localMinimumMean=std::numeric_limits<double>::max();
    double localMaximumMean=0.0;
    int localUnequal=0;
    rnd_seed(830000+1009*PIC::ThisThread+
        7919*static_cast<int>(reactionIndex));

    for (long int event=0;event<localEvents;event++) {
      std::vector<double> velocity;
      GenerateProductVelocities(products,parentVelocity,
          reaction.ExcessEnergyJ,velocity);

      double exactMomentum[3]={0.0,0.0,0.0};
      double exactEnergy=0.0;
      for (std::size_t product=0;product<products.size();product++) {
        const double mass=PIC::MolecularData::GetMass(products[product]);
        const double *v=velocity.data()+3*product;
        for (int idim=0;idim<3;idim++)
          exactMomentum[idim]+=mass*v[idim];
        exactEnergy+=0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
      }

      double momentumError2=0.0;
      for (int idim=0;idim<3;idim++) {
        const double difference=
            exactMomentum[idim]-productMass*parentVelocity[idim];
        momentumError2+=difference*difference;
      }
      const double momentumScale=productMass*std::sqrt(parentSpeed2);
      localMaximumMomentumError=std::max(localMaximumMomentumError,
          std::sqrt(momentumError2)/
              std::max(momentumScale,std::numeric_limits<double>::min()));
      const double expectedExactEnergy=
          0.5*productMass*parentSpeed2+reaction.ExcessEnergyJ;
      localMaximumEnergyError=std::max(localMaximumEnergyError,
          std::fabs(exactEnergy-expectedExactEnergy)/
              std::max(expectedExactEnergy,
                  std::numeric_limits<double>::min()));

      for (std::size_t product=0;product<products.size();product++) {
        const int species=products[product];
        const double targetWeight=node->block->GetLocalParticleWeight(species);
        const double targetTimeStep=node->block->GetLocalTimeStep(species);
        const double targetRateWeight=targetWeight/targetTimeStep;
        const double mean=GetResidualMultiplicityMean(parentWeight,
            parentTimeStep,targetWeight,targetTimeStep);
        const double residual=mean-std::floor(mean);
        const double bernoulliVariance=residual*(1.0-residual);
        const long int multiplicity=DrawResidualMultiplicity(parentWeight,
            parentTimeStep,targetWeight,targetTimeStep);
        const double actualRate=multiplicity*targetRateWeight;
        const double mass=PIC::MolecularData::GetMass(species);
        const double *v=velocity.data()+3*product;
        const double kineticPerMolecule=
            0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

        localMinimumMean=std::min(localMinimumMean,mean);
        localMaximumMean=std::max(localMaximumMean,mean);
        if (std::fabs(mean-1.0)>1.0e-12) localUnequal=1;

        local.ExpectedNumber[species]+=sourceRateWeight;
        local.ActualNumber[species]+=actualRate;
        local.NumberVariance[species]+=bernoulliVariance*
            targetRateWeight*targetRateWeight;
        local.ExpectedMass+=sourceRateWeight*mass;
        local.ActualMass+=actualRate*mass;
        local.MassVariance+=bernoulliVariance*
            std::pow(targetRateWeight*mass,2);
        for (int idim=0;idim<3;idim++) {
          local.ExpectedMomentum[idim]+=sourceRateWeight*mass*v[idim];
          local.ActualMomentum[idim]+=actualRate*mass*v[idim];
          local.MomentumVariance[idim]+=bernoulliVariance*
              std::pow(targetRateWeight*mass*v[idim],2);
        }
        local.ExpectedEnergy+=sourceRateWeight*kineticPerMolecule;
        local.ActualEnergy+=actualRate*kineticPerMolecule;
        local.EnergyVariance+=bernoulliVariance*
            std::pow(targetRateWeight*kineticPerMolecule,2);
      }
    }

    cChemistryEnsembleMoments global;
    MPI_Allreduce(local.ExpectedNumber.data(),global.ExpectedNumber.data(),
        PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.ActualNumber.data(),global.ActualNumber.data(),
        PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.NumberVariance.data(),global.NumberVariance.data(),
        PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    double localScalars[15]={local.ExpectedMass,local.ActualMass,
        local.MassVariance,local.ExpectedMomentum[0],
        local.ExpectedMomentum[1],local.ExpectedMomentum[2],
        local.ActualMomentum[0],local.ActualMomentum[1],
        local.ActualMomentum[2],local.MomentumVariance[0],
        local.MomentumVariance[1],local.MomentumVariance[2],
        local.ExpectedEnergy,local.ActualEnergy,local.EnergyVariance};
    double globalScalars[15]={};
    MPI_Allreduce(localScalars,globalScalars,15,MPI_DOUBLE,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    global.ExpectedMass=globalScalars[0];
    global.ActualMass=globalScalars[1];
    global.MassVariance=globalScalars[2];
    for (int idim=0;idim<3;idim++) {
      global.ExpectedMomentum[idim]=globalScalars[3+idim];
      global.ActualMomentum[idim]=globalScalars[6+idim];
      global.MomentumVariance[idim]=globalScalars[9+idim];
    }
    global.ExpectedEnergy=globalScalars[12];
    global.ActualEnergy=globalScalars[13];
    global.EnergyVariance=globalScalars[14];

    double globalMaximumMomentumError=0.0,globalMaximumEnergyError=0.0;
    double globalMinimumMean=0.0,globalMaximumMean=0.0;
    MPI_Allreduce(&localMaximumMomentumError,&globalMaximumMomentumError,1,
        MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&localMaximumEnergyError,&globalMaximumEnergyError,1,
        MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&localMinimumMean,&globalMinimumMean,1,MPI_DOUBLE,MPI_MIN,
        MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&localMaximumMean,&globalMaximumMean,1,MPI_DOUBLE,MPI_MAX,
        MPI_GLOBAL_COMMUNICATOR);
    int globalUnequal=0;
    MPI_Allreduce(&localUnequal,&globalUnequal,1,MPI_INT,MPI_MAX,
        MPI_GLOBAL_COMMUNICATOR);

    reactionPass=(globalMaximumMomentumError<=1.0e-11 &&
        globalMaximumEnergyError<=1.0e-11) && reactionPass;
    for (int species=0;species<PIC::nTotalSpecies;species++)
      if (global.ExpectedNumber[species]>0.0) {
        char symbol[_MAX_STRING_LENGTH_PIC_];
        char name[_MAX_STRING_LENGTH_PIC_+16];
        PIC::MolecularData::GetChemSymbol(symbol,species);
        std::snprintf(name,sizeof(name),"number[%s]",symbol);
        reactionPass=CheckEnsembleMoment(reaction.Id,name,
            global.ExpectedNumber[species],global.ActualNumber[species],
            global.NumberVariance[species]) && reactionPass;
      }
    reactionPass=CheckEnsembleMoment(reaction.Id,"mass",
        global.ExpectedMass,global.ActualMass,global.MassVariance) &&
        reactionPass;
    for (int idim=0;idim<3;idim++) {
      char name[32];
      std::snprintf(name,sizeof(name),"momentum[%d]",idim);
      reactionPass=CheckEnsembleMoment(reaction.Id,name,
          global.ExpectedMomentum[idim],global.ActualMomentum[idim],
          global.MomentumVariance[idim]) && reactionPass;
    }
    reactionPass=CheckEnsembleMoment(reaction.Id,"kinetic+excess",
        global.ExpectedEnergy,global.ActualEnergy,global.EnergyVariance) &&
        reactionPass;

    testedReactionCount++;
    if (globalUnequal) unequalReactionCount++;
    allPass=reactionPass && allPass;
    if (PIC::ThisThread==0)
      std::printf("%s chemistry reaction=%s events=%ld products=%zu "
          "multiplicity-mean=[%.6g,%.6g] catalog-mass-relative=%.3e "
          "kinematic-momentum=%.3e kinematic-energy=%.3e\n",
          reactionPass ? "PASS" : "FAIL",reaction.Id,eventsPerReaction,
          products.size(),globalMinimumMean,globalMaximumMean,
          catalogMassRelative,globalMaximumMomentumError,
          globalMaximumEnergyError);
  }

  allPass=(testedReactionCount>0 && unequalReactionCount>0) && allPass;
  int localPass=allPass ? 1 : 0,globalPass=0;
  MPI_Allreduce(&localPass,&globalPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0)
    std::printf("Chemistry ensemble validation result: %s "
        "neutral-reactions=%d unequal-weight/timestep-reactions=%d "
        "sink-reactions=%d\n",globalPass ? "PASS" : "FAIL",
        testedReactionCount,unequalReactionCount,sinkReactionCount);
  return globalPass!=0;
}

void Comet::Chemistry::BeginBatchedModel() {
  LocalChemistryAllocationPressure=false;

  if (!ParentReactionRatesReady ||
      CachedHeliocentricDistanceAU!=HeliocentricDistanceAU) {
    std::string error;
    if (!RefreshParentReactionRates(&error)) {
      std::fprintf(stderr,"Cannot refresh comet chemistry rates: %s\n",
          error.c_str());
      std::abort();
    }
  }

  // Reset the per-stage flag on every owner particle, including particles
  // injected immediately before chemistry.  A byte is used instead of an
  // ever-growing epoch so a reused particle-buffer slot can never inherit a
  // value that accidentally marks a new physical particle as completed.
  constexpr long int blockCellCount=
      _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
  for (unsigned int blockIndex=0;
       blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;blockIndex++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node==NULL || node->block==NULL) continue;
    for (long int cell=0;cell<blockCellCount;cell++)
      for (long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;
           ptr=PIC::ParticleBuffer::GetNext(ptr))
        ChemistryProcessedFlag(
            PIC::ParticleBuffer::GetParticleDataPointer(ptr))=0;
  }
}

bool Comet::Chemistry::IsBatchedParticleEligible(long int ptr) {
  return ChemistryProcessedFlag(
      PIC::ParticleBuffer::GetParticleDataPointer(ptr))==0;
}

long int Comet::Chemistry::GetBatchedParticleLimit() {
  const long int capacity=PIC::ParticleBuffer::MaxNPart;
  if (capacity<=0 ||
      Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction<=0) {
    std::fprintf(stderr,
        "invalid particle capacity or chemistry event limit in batch controller\n");
    std::abort();
  }

  // The high-water mark, rather than the physical end of the buffer, defines
  // the maximum population after one pass.  Because every eligible parent is
  // charged the hard worst-case event allocation, the pass cannot cross this
  // point even if every tested parent reacts and every residual Bernoulli
  // draw rounds upward.  The separately configured reserve is checked again
  // by ProcessParticle before an individual transaction is committed.
  const long int highWaterSlots=static_cast<long int>(std::floor(
      Comet::RadialSimulation::EmergencyParticleBufferUtilization*capacity));
  const long int available=highWaterSlots-PIC::ParticleBuffer::NAllPart;
  if (available<=0) {
    LocalChemistryAllocationPressure=true;
    return 0;
  }
  const long int limit=available/
      Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction;
  if (limit<=0) {
    LocalChemistryAllocationPressure=true;
    return 0;
  }
  return limit;
}

void Comet::Chemistry::FinishBatchedPass(long int attemptedParents,
    long int remainingParents) {
  (void)attemptedParents;
  const long long localPopulation=PIC::ParticleBuffer::NAllPart;
  long long maximumPopulation=0;
  MPI_Allreduce(&localPopulation,&maximumPopulation,1,MPI_LONG_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  int localPressure=(LocalChemistryAllocationPressure ||
      (remainingParents>0 && attemptedParents==0)) ? 1 : 0;
  int globalPressure=0;
  MPI_Allreduce(&localPressure,&globalPressure,1,MPI_INT,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);

  const double maximumUtilization=static_cast<double>(maximumPopulation)/
      PIC::ParticleBuffer::MaxNPart;
  if (globalPressure!=0 || maximumUtilization>=
      Comet::RadialSimulation::EmergencyParticleBufferUtilization)
    Comet::RadialSimulation::RequestAllocationBurstEmergency();

  LocalChemistryAllocationPressure=false;
}

bool Comet::Chemistry::ProcessParticle(long int ptr,long int &firstParticleCell,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (!Enabled) {
    MarkChemistryProcessed(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
    KeepParticle(ptr,firstParticleCell);
    return true;
  }

  PIC::ParticleBuffer::byte *parentData=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  const int parentSpecies=PIC::ParticleBuffer::GetI(parentData);
  const ChemicalSolver::cReactionNetwork &network=GetNetwork();

  if (!ParentReactionRatesReady ||
      CachedHeliocentricDistanceAU!=HeliocentricDistanceAU ||
      parentSpecies<0 || parentSpecies>=PIC::nTotalSpecies) {
    std::fprintf(stderr,
        "comet chemistry rate cache is stale or particle species %d is invalid\n",
        parentSpecies);
    std::abort();
  }
  const cParentReactionRateCache &parentRates=
      ParentReactionRates[parentSpecies];
  const double totalRate=parentRates.TotalRate;

  // The approved Stage-0 definition counts one eligible visit for every
  // positive-rate channel applicable to this parent.  Deferred parents are
  // intentionally counted again when a later batch revisits them; this is an
  // executor-work metric, while committed selections below remain physical
  // event counts and are never double counted.
  if (Comet::Diagnostics::Enabled)
    for (std::size_t channel=0;channel<parentRates.Channels.size();channel++)
      Comet::Diagnostics::RecordEligibleParentVisit(
          parentRates.Channels[channel]);

  if (!(totalRate>0.0)) {
    MarkChemistryProcessed(parentData);
    KeepParticle(ptr,firstParticleCell);
    return true;
  }

  const double timeStep=node->block->GetLocalTimeStep(parentSpecies);
  if (!(timeStep>0.0)) {
    std::fprintf(stderr,"invalid comet chemistry timestep for species %d\n",
        parentSpecies);
    std::abort();
  }
  if (rnd()>=1.0-std::exp(-totalRate*timeStep)) {
    MarkChemistryProcessed(parentData);
    KeepParticle(ptr,firstParticleCell);
    return true;
  }

  double selector=rnd()*totalRate;
  std::size_t selected=parentRates.Channels.back();
  for (std::size_t i=0;i<parentRates.Channels.size();i++) {
    selector-=parentRates.Rates[i];
    if (selector<=0.0) {
      selected=parentRates.Channels[i];
      break;
    }
  }
  const ChemicalSolver::cReaction &reaction=network.GetReaction(selected);

  // Expand integer stoichiometric coefficients into concrete daughter
  // particles. A zero-product channel is a photoionization loss from this
  // neutral-only network and therefore simply removes the parent.
  std::vector<int> productSpecies;
  std::vector<std::size_t> productTermIndex;
  for (std::size_t term=0;term<reaction.ProductCount;term++) {
    const ChemicalSolver::cStoichiometricTerm &product=
        network.GetProductTerm(selected,term);
    const int multiplicity=static_cast<int>(std::floor(product.Coefficient+0.5));
    for (int copy=0;copy<multiplicity;copy++) {
      productSpecies.push_back(product.Species);
      productTermIndex.push_back(term);
    }
  }

  double xParent[3],vParent[3];
  PIC::ParticleBuffer::GetX(xParent,parentData);
  PIC::ParticleBuffer::GetV(vParent,parentData);
  std::vector<double> productVelocity;
  GenerateProductVelocities(productSpecies,vParent,reaction.ExcessEnergyJ,
      productVelocity);

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  const double parentPhysicalWeight=
      node->block->GetLocalParticleWeight(parentSpecies)*
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(parentData);
#elif (_SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_) || \
      (_SIMULATION_PARTICLE_WEIGHT_MODE_ == _SINGLE_GLOBAL_PARTICLE_WEIGHT_) || \
      (_SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_)
  const double parentPhysicalWeight=
      node->block->GetLocalParticleWeight(parentSpecies);
#else
#error "Comet chemistry without individual corrections requires a supported global or species-local particle weight"
#endif

  // Draw every residual multiplicity before changing the parent or allocating
  // a daughter.  This makes the reaction a transaction: either capacity for
  // the complete event is available, or the untouched parent is reinserted
  // and retried after the collective batch checkpoint creates headroom.
  std::vector<long int> daughterMultiplicity(productSpecies.size(),0);
  long int totalDaughters=0;
  for (std::size_t product=0;product<productSpecies.size();product++) {
    const int species=productSpecies[product];
    const double baseWeight=node->block->GetLocalParticleWeight(species);
    const double productTimeStep=node->block->GetLocalTimeStep(species);
    if (!(baseWeight>0.0) || !std::isfinite(baseWeight) ||
        !(productTimeStep>0.0) || !std::isfinite(productTimeStep)) {
      std::fprintf(stderr,
          "invalid particle weight or timestep for comet chemistry species %d\n",
          species);
      std::abort();
    }

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    daughterMultiplicity[product]=1;
#else
    daughterMultiplicity[product]=
        DrawResidualMultiplicity(parentPhysicalWeight,timeStep,baseWeight,
            productTimeStep);
#endif

    if (daughterMultiplicity[product]<0 ||
        totalDaughters>LONG_MAX-daughterMultiplicity[product]) {
      std::fprintf(stderr,
          "photochemical daughter allocation count overflow in reaction %s\n",
          reaction.Id);
      std::abort();
    }
    totalDaughters+=daughterMultiplicity[product];
  }

  // RadialMaxPhotochemicalProductsPerReaction is a hard actual allocation
  // limit, not merely an expectation.  Sampling-break carrier planning is
  // responsible for making this invariant true for present as well as absent
  // products.  Never truncate here: truncation would bias reaction rates.
  if (totalDaughters>
      Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction) {
    std::fprintf(stderr,
        "reaction %s requests %ld model daughters; hard limit=%ld parent-species=%d W/dt=%.16e\n",
        reaction.Id,totalDaughters,
        Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction,
        parentSpecies,parentPhysicalWeight/timeStep);
    std::abort();
  }

  // Commit diagnostics only after the transaction is known to fit.  A
  // capacity-deferred parent is retried by the core and must not appear as a
  // second physical reaction.  Expected product rates follow stoichiometry
  // from the selected parent carrier; represented rates follow the daughters
  // actually produced by residual split/roulette.
  auto recordCommittedDiagnostics=[&]() {
    if (!Comet::Diagnostics::Enabled) return;
    const double parentRate=parentPhysicalWeight/timeStep;
    Comet::Diagnostics::RecordReactionSelection(selected,
        reaction.ProductCount==0);
    Comet::Diagnostics::RecordChemicalDestruction(parentSpecies,parentRate);

    struct cTermDiagnostic {
      std::uint64_t Multiplicity=0,Allocated=0,Reused=0;
      double ExpectedNumber=0.0,ExpectedMass=0.0;
      double ExpectedMomentum[3]={0.0,0.0,0.0};
      double ExpectedEnergy=0.0;
      double RepresentedNumber=0.0,RepresentedMass=0.0;
      double RepresentedMomentum[3]={0.0,0.0,0.0};
      double RepresentedEnergy=0.0,NumberVariance=0.0;
    };
    std::vector<cTermDiagnostic> terms(reaction.ProductCount);
    bool reuseAvailable=totalDaughters>0;
    for (std::size_t product=0;product<productSpecies.size();product++) {
      const std::size_t termIndex=productTermIndex[product];
      cTermDiagnostic &term=terms[termIndex];
      const int species=productSpecies[product];
      const double mass=PIC::MolecularData::GetMass(species);
      const double *velocity=productVelocity.data()+3*product;
      const double speedSquared=velocity[0]*velocity[0]+
          velocity[1]*velocity[1]+velocity[2]*velocity[2];
      const std::uint64_t multiplicity=static_cast<std::uint64_t>(
          daughterMultiplicity[product]);
      term.Multiplicity+=multiplicity;
      if (multiplicity>0) {
        if (reuseAvailable) {
          term.Reused++;
          term.Allocated+=multiplicity-1;
          reuseAvailable=false;
        }
        else term.Allocated+=multiplicity;
      }

      term.ExpectedNumber+=parentRate;
      term.ExpectedMass+=parentRate*mass;
      for (int component=0;component<3;component++)
        term.ExpectedMomentum[component]+=
            parentRate*mass*velocity[component];
      term.ExpectedEnergy+=parentRate*0.5*mass*speedSquared;

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      const double representedRate=parentRate;
#else
      const double productCarrierRate=
          node->block->GetLocalParticleWeight(species)/
          node->block->GetLocalTimeStep(species);
      const double representedRate=multiplicity*productCarrierRate;
      const double mean=GetResidualMultiplicityMean(parentPhysicalWeight,
          timeStep,node->block->GetLocalParticleWeight(species),
          node->block->GetLocalTimeStep(species));
      const double fraction=mean-std::floor(mean);
      term.NumberVariance+=productCarrierRate*productCarrierRate*
          fraction*(1.0-fraction);
#endif
      term.RepresentedNumber+=representedRate;
      term.RepresentedMass+=representedRate*mass;
      for (int component=0;component<3;component++)
        term.RepresentedMomentum[component]+=
            representedRate*mass*velocity[component];
      term.RepresentedEnergy+=representedRate*0.5*mass*speedSquared;
    }

    for (std::size_t product=0;product<reaction.ProductCount;product++) {
      const ChemicalSolver::cStoichiometricTerm &networkTerm=
          network.GetProductTerm(selected,product);
      const cTermDiagnostic &term=terms[product];
      Comet::Diagnostics::RecordChemistryProduct(selected,product,
          term.Multiplicity,term.Allocated,term.Reused,term.ExpectedNumber,
          term.ExpectedMass,term.ExpectedMomentum,term.ExpectedEnergy,
          term.RepresentedNumber,term.RepresentedMass,
          term.RepresentedMomentum,term.RepresentedEnergy,
          term.NumberVariance);
      Comet::Diagnostics::RecordChemicalProduction(networkTerm.Species,
          term.ExpectedNumber,term.RepresentedNumber,term.NumberVariance);
    }
  };

  if (totalDaughters==0) {
    recordCommittedDiagnostics();
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return true;
  }

  // Reuse the parent slot for the first daughter.  Only K-1 new slots are
  // needed for a K-daughter event, which both lowers peak memory and makes the
  // capacity preflight exact.  The reserve is a hard ceiling independent of
  // the lower 75% checkpoint: even a future caller with a larger batch cannot
  // consume slots retained for MPI receipt and other core operations.
  const long int transactionCeiling=static_cast<long int>(std::floor(
      (1.0-Comet::RadialSimulation::AllocationReserveFraction)*
      PIC::ParticleBuffer::MaxNPart));
  const long int requiredNewSlots=totalDaughters-1;
  if (requiredNewSlots>transactionCeiling-PIC::ParticleBuffer::NAllPart) {
    LocalChemistryAllocationPressure=true;
    Comet::Diagnostics::RecordAllocationDeferral(selected);
    KeepParticle(ptr,firstParticleCell);
    return false;
  }

  recordCommittedDiagnostics();

  // Clone auxiliary fields from the parent, then replace position, velocity,
  // and species for each integer-stoichiometry product.  Fixed block-local
  // weights use residual split/roulette; number, mass, momentum, and excess-
  // energy rates are conserved in expectation without an individual weight
  // correction field.
  char daughterTemplate[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated(
      reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate));
  PIC::ParticleBuffer::SetPrev(-1,
      reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate));
  PIC::ParticleBuffer::SetNext(-1,
      reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate));
  PIC::ParticleBuffer::CloneParticle(
      reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate),parentData);

  bool parentSlotReused=false;
  for (std::size_t product=0;product<productSpecies.size();product++) {
    const int species=productSpecies[product];
    const double baseWeight=node->block->GetLocalParticleWeight(species);
    const double productTimeStep=node->block->GetLocalTimeStep(species);

    PIC::ParticleBuffer::byte *data=
        reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate);
    PIC::ParticleBuffer::SetX(xParent,data);
    PIC::ParticleBuffer::SetV(productVelocity.data()+3*product,data);
    PIC::ParticleBuffer::SetI(species,data);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    // A single corrected daughter carries the same physical production rate
    // as its parent event, including the species-local timestep conversion.
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
        parentPhysicalWeight*productTimeStep/(baseWeight*timeStep),data);
#endif
    MarkChemistryProcessed(data);

    // Each residual-resampling sibling represents the same physical product
    // cohort and therefore starts with the same event velocity.  Subsequent
    // motion and collisions decorrelate the model particles.  Assign a fresh
    // tracker identity to every sibling rather than cloning one identity.
    for (long int copy=0;copy<daughterMultiplicity[product];copy++) {
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      PIC::ParticleTracker::InitParticleID(daughterTemplate);
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(
          xParent,productVelocity.data()+3*product,species,daughterTemplate,
          static_cast<void*>(node));
#endif
      if (!parentSlotReused) {
        // Detach the parent slot's stale links before CloneParticle, which
        // intentionally preserves the destination links.  Reinsert it into
        // the newly rebuilt cell list exactly like an allocated daughter.
        PIC::ParticleBuffer::SetPrev(-1,parentData);
        PIC::ParticleBuffer::SetNext(-1,parentData);
        PIC::ParticleBuffer::CloneParticle(parentData,data);
        KeepParticle(ptr,firstParticleCell);
        parentSlotReused=true;
      }
      else {
        const long int daughter=
            PIC::ParticleBuffer::GetNewParticle(firstParticleCell);
        PIC::ParticleBuffer::CloneParticle(
            PIC::ParticleBuffer::GetParticleDataPointer(daughter),data);
      }
    }
  }
  return true;
}

#include "CometChemistry.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

#include "pic.h"
#include "src/models/chemical_solver/SystematicResidualAccumulator.h"
#include "CometChemistry.generated.h"
#include "CometDaughterPacketPolicy.generated.h"
#include "CometDiagnostics.h"
#include "DaughterPacketPlanner.h"
#include "RadialSimulation.h"

bool Comet::Chemistry::Enabled=false;
double Comet::Chemistry::HeliocentricDistanceAU=-1.0;
// Keep this definition on one physical line.  The AMPS input preprocessor
// replaces complete variable-definition lines; a wrapped initializer would
// leave its continuation behind when ChemistryResidualSampling is parsed.
int Comet::Chemistry::ResidualSamplingMode=Comet::Chemistry::ResidualSampling::IndependentBernoulli;
std::uint64_t Comet::Chemistry::ResidualPhaseSeed=1;
bool Comet::Chemistry::DaughterPacketSplittingEnabled=false;
long int Comet::Chemistry::MaximumDaughterPacketCount=1;

static_assert(
    Comet::Chemistry::DaughterPacketGenerated::ReactionCount==
        Comet::Chemistry::Generated::ReactionCount,
    "generated daughter-packet policy does not match reaction network");

namespace {
long int ChemistryProcessedFlagOffset=-1;
bool LocalChemistryAllocationPressure=false;

// Systematic residual state occupies only the block's application-owned data
// region, so ordinary AMR ownership migration carries it byte-for-byte. The
// header and each concrete-product record use doubles because the core's
// block clean routine deliberately treats all application storage as an
// array of doubles. Flags are consequently exact 0.0/1.0 values.
constexpr std::size_t ResidualHeaderInitialized=0;
constexpr std::size_t ResidualHeaderCarrierEpoch=1;
constexpr std::size_t ResidualHeaderLayoutSignature=2;
constexpr std::size_t ResidualHeaderLength=3;
constexpr std::size_t ResidualRecordPhase=0;
constexpr std::size_t ResidualRecordDiagnosticActive=1;
constexpr std::size_t ResidualRecordLength=2;
// Bump this whenever header/record byte semantics change even if the total
// number of doubles does not.  Network contents are mixed in separately.
constexpr std::uint64_t ResidualLayoutAbiVersion=1;
long int ResidualBlockDataOffset=-1;
std::vector<std::size_t> ResidualReactionOffset;
std::vector<int> ResidualConcreteSpecies;
std::vector<std::size_t> ResidualConcreteTerm;
double ResidualLayoutSignature=0.0;

void BuildResidualLayout() {
  if (!ResidualReactionOffset.empty()) return;

  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  ResidualReactionOffset.assign(network.GetReactionCount()+1,0);
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
    ResidualReactionOffset[reaction]=ResidualConcreteSpecies.size();
    const ChemicalSolver::cReaction &entry=network.GetReaction(reaction);
    for (std::size_t term=0;term<entry.ProductCount;term++) {
      const ChemicalSolver::cStoichiometricTerm &product=
          network.GetProductTerm(reaction,term);
      const long int coefficient=static_cast<long int>(
          std::floor(product.Coefficient+0.5));
      for (long int copy=0;copy<coefficient;copy++) {
        ResidualConcreteSpecies.push_back(product.Species);
        ResidualConcreteTerm.push_back(term);
      }
    }
  }
  ResidualReactionOffset[network.GetReactionCount()]=
      ResidualConcreteSpecies.size();

  // The restart core validates byte length, while the application knows the
  // semantic ordering of those bytes.  Store a 53-bit, exactly representable
  // fingerprint in every initialized block so a reordered reaction network
  // with the same product count cannot silently reinterpret saved phases.
  namespace Accumulator=ChemicalSolver::SystematicResidualAccumulator;
  std::uint64_t signature=Accumulator::Mix(ResidualLayoutAbiVersion);
  signature=Accumulator::Combine(signature,ResidualHeaderLength);
  signature=Accumulator::Combine(signature,ResidualRecordLength);
  signature=Accumulator::Combine(signature,ResidualHeaderInitialized);
  signature=Accumulator::Combine(signature,ResidualHeaderCarrierEpoch);
  signature=Accumulator::Combine(signature,ResidualHeaderLayoutSignature);
  signature=Accumulator::Combine(signature,ResidualRecordPhase);
  signature=Accumulator::Combine(signature,ResidualRecordDiagnosticActive);
  signature=Accumulator::Combine(signature,network.GetReactionCount());
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
    const ChemicalSolver::cReaction &entry=network.GetReaction(reaction);
    for (const unsigned char *p=
        reinterpret_cast<const unsigned char*>(entry.Id);*p!='\0';p++)
      signature=Accumulator::Combine(signature,*p);
    signature=Accumulator::Combine(signature,ResidualReactionOffset[reaction]);
    signature=Accumulator::Combine(signature,
        ResidualReactionOffset[reaction+1]-ResidualReactionOffset[reaction]);
  }
  for (std::size_t product=0;product<ResidualConcreteSpecies.size();product++) {
    signature=Accumulator::Combine(signature,
        static_cast<std::uint64_t>(ResidualConcreteSpecies[product]));
    signature=Accumulator::Combine(signature,ResidualConcreteTerm[product]);
  }
  signature&=(UINT64_C(1)<<53)-1;
  if (signature==0) signature=1;
  ResidualLayoutSignature=static_cast<double>(signature);
}

std::size_t ResidualRecordIndex(std::size_t concreteProduct,
    std::size_t field) {
  return ResidualHeaderLength+ResidualRecordLength*concreteProduct+field;
}

double* GetResidualBlockData(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (ResidualBlockDataOffset<0 || node==NULL || node->block==NULL ||
      node->Thread!=PIC::ThisThread) {
    std::fprintf(stderr,
        "systematic comet residual state was requested from a non-owner or uninitialized block\n");
    std::abort();
  }
  char *raw=node->block->GetAssociatedDataBufferPointer()+
      ResidualBlockDataOffset;
  if (reinterpret_cast<std::uintptr_t>(raw)%alignof(double)!=0) {
    std::fprintf(stderr,
        "systematic comet residual block state is not double aligned\n");
    std::abort();
  }
  return reinterpret_cast<double*>(raw);
}

std::uint64_t ResidualNodeKey(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  namespace Accumulator=ChemicalSolver::SystematicResidualAccumulator;
  std::uint64_t key=Accumulator::Combine(
      Comet::Chemistry::ResidualPhaseSeed,
      static_cast<std::uint64_t>(node->AMRnodeID.ResolutionLevel));
  for (unsigned int i=0;i<node->AMRnodeID.Length();i++)
    key=Accumulator::Combine(key,node->AMRnodeID.id[i]);
  return key;
}

double ResidualCarrierRate(PIC::Mesh::cDataBlockAMR *block,int species) {
  const double weight=block->GetLocalParticleWeight(species);
  const double timeStep=block->GetLocalTimeStep(species);
  if (!(weight>0.0) || !std::isfinite(weight) || !(timeStep>0.0) ||
      !std::isfinite(timeStep)) {
    std::fprintf(stderr,
        "invalid block carrier while accessing systematic comet residual state\n");
    std::abort();
  }
  return weight/timeStep;
}

void BuildInitialResidualPhases(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,
    std::vector<double> &phases) {
  namespace Accumulator=ChemicalSolver::SystematicResidualAccumulator;
  double *data=GetResidualBlockData(node);
  if (data[ResidualHeaderInitialized]!=0.0 ||
      !(data[ResidualHeaderCarrierEpoch]>=0.0) ||
      !std::isfinite(data[ResidualHeaderCarrierEpoch]) ||
      (data[ResidualHeaderLayoutSignature]!=0.0 &&
       data[ResidualHeaderLayoutSignature]!=ResidualLayoutSignature)) {
    std::fprintf(stderr,"invalid systematic residual block header\n");
    std::abort();
  }
  phases.assign(ResidualConcreteSpecies.size(),0.0);

  const std::uint64_t epoch=static_cast<std::uint64_t>(
      data[ResidualHeaderCarrierEpoch]);
  const std::uint64_t nodeKey=ResidualNodeKey(node);
  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();

  // Every concrete product retains its own state slot. Within a reaction,
  // exact-equal W/dt carriers receive one common initial packet phase. Since
  // they add the same mu, they remain in lockstep until an accepted regular
  // carrier transaction starts a new epoch. Different carriers intentionally
  // diverge and conserve production only in cumulative expectation.
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
    const std::size_t begin=ResidualReactionOffset[reaction];
    const std::size_t end=ResidualReactionOffset[reaction+1];
    for (std::size_t product=begin;product<end;product++) {
      const double carrier=ResidualCarrierRate(node->block,
          ResidualConcreteSpecies[product]);
      std::size_t leader=product;
      for (std::size_t candidate=begin;candidate<product;candidate++)
        if (ResidualCarrierRate(node->block,
            ResidualConcreteSpecies[candidate])==carrier) {
          leader=candidate;
          break;
        }

      double phase;
      if (leader!=product)
        phase=phases[leader];
      else {
        std::uint64_t key=Accumulator::Combine(nodeKey,epoch);
        key=Accumulator::Combine(key,reaction);
        key=Accumulator::Combine(key,product-begin);
        phase=Accumulator::UnitPhase(key);
      }
      phases[product]=phase;
    }
  }
}

void InitializeResidualBlock(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double *data=GetResidualBlockData(node);
  if (data[ResidualHeaderInitialized]==1.0) {
    if (data[ResidualHeaderLayoutSignature]!=ResidualLayoutSignature) {
      std::fprintf(stderr,
          "systematic residual restart layout fingerprint mismatch\n");
      std::abort();
    }
    return;
  }
  std::vector<double> phases;
  BuildInitialResidualPhases(node,phases);
  for (std::size_t product=0;product<phases.size();product++) {
    data[ResidualRecordIndex(product,ResidualRecordPhase)]=phases[product];
    data[ResidualRecordIndex(product,ResidualRecordDiagnosticActive)]=0.0;
  }
  data[ResidualHeaderLayoutSignature]=ResidualLayoutSignature;
  data[ResidualHeaderInitialized]=1.0;
}

void CloseResidualBlockDiagnosticEpoch(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,bool resetCarrierEpoch) {
  double *data=GetResidualBlockData(node);
  if (data[ResidualHeaderInitialized]==1.0) {
    if (data[ResidualHeaderLayoutSignature]!=ResidualLayoutSignature) {
      std::fprintf(stderr,
          "systematic residual diagnostic layout fingerprint mismatch\n");
      std::abort();
    }
    for (std::size_t product=0;product<ResidualConcreteSpecies.size();product++) {
      double &active=data[ResidualRecordIndex(product,
          ResidualRecordDiagnosticActive)];
      if (active==1.0 && Comet::Diagnostics::Enabled) {
        const auto upper=std::upper_bound(ResidualReactionOffset.begin(),
            ResidualReactionOffset.end(),product);
        const std::size_t reaction=static_cast<std::size_t>(
            upper-ResidualReactionOffset.begin()-1);
        const int species=ResidualConcreteSpecies[product];
        Comet::Diagnostics::RecordSystematicResidualBound(reaction,
            ResidualConcreteTerm[product],species,
            ResidualCarrierRate(node->block,species));
      }
      active=0.0;
    }
  }

  if (resetCarrierEpoch) {
    const double nextEpoch=data[ResidualHeaderCarrierEpoch]+1.0;
    if (!(nextEpoch>data[ResidualHeaderCarrierEpoch]) ||
        !std::isfinite(nextEpoch)) {
      std::fprintf(stderr,"systematic residual carrier epoch overflow\n");
      std::abort();
    }
    data[ResidualHeaderCarrierEpoch]=nextEpoch;
    data[ResidualHeaderInitialized]=0.0;
  }
}

void ForEachOwnerResidualBlock(void (*function)(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,bool),bool argument) {
  if (Comet::Chemistry::ResidualSamplingMode!=
      Comet::Chemistry::ResidualSampling::RandomPhaseSystematic) return;
  for (unsigned int block=0;
       block<PIC::DomainBlockDecomposition::nLocalBlocks;block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[block];
    if (node==NULL || node->block==NULL || node->Thread!=PIC::ThisThread)
      continue;
    function(node,argument);
  }
}

std::size_t ResidualStateDoubleCount() {
  return ResidualHeaderLength+
      ResidualRecordLength*ResidualConcreteSpecies.size();
}

void CollectResidualValidationLeaves(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,
    std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> &leaves) {
  if (node==NULL) return;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->IsUsedInCalculationFlag) leaves.push_back(node);
    return;
  }
  for (int child=0;child<(1<<DIM);child++)
    CollectResidualValidationLeaves(node->downNode[child],leaves);
}

double Stage3MigrationLoadMeasure(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  namespace Accumulator=ChemicalSolver::SystematicResidualAccumulator;
  std::uint64_t key=static_cast<std::uint64_t>(
      node->AMRnodeID.ResolutionLevel);
  for (unsigned int i=0;i<node->AMRnodeID.Length();i++)
    key=Accumulator::Combine(key,node->AMRnodeID.id[i]);
  // A deliberately irregular but deterministic measure changes multi-rank
  // ownership boundaries without depending on the previous MPI owner.
  return 1.0+static_cast<double>(Accumulator::Mix(key)%97);
}

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

void Comet::Chemistry::RequestBlockData() {
  if (!Enabled || ResidualSamplingMode!=
      ResidualSampling::RandomPhaseSystematic) return;
  if (ResidualBlockDataOffset>=0) return;

  BuildResidualLayout();
  const std::size_t doubleCount=ResidualHeaderLength+
      ResidualRecordLength*ResidualConcreteSpecies.size();
  if (doubleCount>static_cast<std::size_t>(INT_MAX)/sizeof(double)) {
    std::fprintf(stderr,
        "systematic comet residual block state exceeds the core offset range\n");
    std::abort();
  }
  // RequestInternalBlockData() is byte-granular and makes no alignment
  // promise.  Insert explicit padding before exposing the storage as doubles;
  // the runtime accessor also checks the actual allocated base address.
  const int currentLength=
      PIC::Mesh::cDataBlockAMR_static_data::totalAssociatedDataLength;
  const int alignment=static_cast<int>(alignof(double));
  const int padding=(alignment-currentLength%alignment)%alignment;
  if (padding!=0)
    PIC::Mesh::cDataBlockAMR::RequestInternalBlockData(padding);
  ResidualBlockDataOffset=PIC::Mesh::cDataBlockAMR::RequestInternalBlockData(
      static_cast<int>(doubleCount*sizeof(double)));
  if (ResidualBlockDataOffset%alignment!=0) {
    std::fprintf(stderr,
        "cannot align systematic comet residual block storage\n");
    std::abort();
  }
}

void Comet::Chemistry::PrepareForRegularCarrierChange() {
  // This hook is called only after the complete component transaction has
  // passed capacity preflight and immediately before any carrier is changed.
  // Closing the active error bound under the old carriers and incrementing
  // the block epoch guarantees a fresh, independently keyed phase after the
  // new W/dt ratios have been committed.
  ForEachOwnerResidualBlock(CloseResidualBlockDiagnosticEpoch,true);
}

void Comet::Chemistry::CloseResidualDiagnosticWindow() {
  // Keep the physical phases continuous across output windows. Only the
  // activity markers used to form a conservative per-window discrepancy
  // bound are cleared.
  ForEachOwnerResidualBlock(CloseResidualBlockDiagnosticEpoch,false);
}

const ChemicalSolver::cReactionNetwork& Comet::Chemistry::GetNetwork() {
  static const ChemicalSolver::cReactionNetwork network(
      Generated::Terms,Generated::TermCount,
      Generated::Reactions,Generated::ReactionCount,
      PIC::nTotalSpecies);
  return network;
}

long int Comet::Chemistry::GetDaughterPacketCount(std::size_t reaction) {
  if (reaction>=GetNetwork().GetReactionCount()) {
    std::fprintf(stderr,"daughter-packet reaction index is out of range\n");
    std::abort();
  }
  const long int count=DaughterPacketSplittingEnabled ?
      DaughterPacketGenerated::PacketCount[reaction] : 1;
  if (count<1 ||
      (DaughterPacketSplittingEnabled &&
       count>MaximumDaughterPacketCount)) {
    std::fprintf(stderr,
        "invalid generated daughter-packet count for reaction %zu\n",
        reaction);
    std::abort();
  }
  return count;
}

int Comet::Chemistry::GetDaughterPacketPolicySource(std::size_t reaction) {
  if (reaction>=GetNetwork().GetReactionCount()) {
    std::fprintf(stderr,"daughter-packet policy index is out of range\n");
    std::abort();
  }
  return DaughterPacketSplittingEnabled ?
      DaughterPacketGenerated::PolicySource[reaction] : 0;
}

bool Comet::Chemistry::Validate(std::string *error) {
  if (ResidualSamplingMode!=ResidualSampling::IndependentBernoulli &&
      ResidualSamplingMode!=ResidualSampling::RandomPhaseSystematic) {
    if (error!=NULL)
      *error="unknown comet chemistry residual-sampling mode";
    return false;
  }
  if (Enabled && !Comet::RadialSimulation::Enabled) {
    // The generated reaction database is geometry neutral, but its current
    // transactional particle executor still obtains allocation limits,
    // reserve policy, and emergency adaptation from the radial controller.
    // Reject a 3-D request explicitly instead of silently falling back to the
    // legacy photolysis path or aborting with unset (-1) capacity controls.
    if (error!=NULL)
      *error="generic comet chemistry execution currently requires the radial backend; define a geometry-neutral batch policy before enabling it in 2-D/3-D";
    return false;
  }
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  if (ResidualSamplingMode==ResidualSampling::RandomPhaseSystematic) {
    if (error!=NULL)
      *error="systematic residual sampling requires individual particle corrections off";
    return false;
  }
#endif
  if (DaughterPacketSplittingEnabled) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    if (error!=NULL)
      *error="daughter packet splitting requires individual particle corrections";
    return false;
#else
    if (!Enabled) {
      if (error!=NULL)
        *error="daughter packet splitting requires comet chemistry";
      return false;
    }
    if (ResidualSamplingMode!=ResidualSampling::IndependentBernoulli) {
      if (error!=NULL)
        *error="daughter packet splitting requires Bernoulli chemistry residual mode";
      return false;
    }
    if (MaximumDaughterPacketCount<2) {
      if (error!=NULL)
        *error="daughter packet splitting requires a maximum of at least two";
      return false;
    }
    if (!Comet::RadialSimulation::IndividualCorrectionWindowsEnabled) {
      if (error!=NULL)
        *error="daughter packet splitting requires Stage-5 correction windows";
      return false;
    }
#endif
  }
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
  bool effectivePacketSplit=false;
  if (Enabled) for (std::size_t reaction=0;
      reaction<GetNetwork().GetReactionCount();reaction++) {
    if (!GetNetwork().IsUnitUnimolecular(reaction)) {
      if (error!=NULL)
        *error="the PIC comet chemistry executor accepts only unit unimolecular reactions";
      return false;
    }

    const ChemicalSolver::cReaction &entry=GetNetwork().GetReaction(reaction);
    const long int packetCount=GetDaughterPacketCount(reaction);
    // Generated rules may cover a database reaction disabled by a zero rate.
    // Such a row is deterministic metadata but cannot make packet mode
    // operational, so it does not satisfy the effective-split startup gate.
    if (packetCount>1 && GetNetwork().GetRateCoefficient(reaction,0.0,
        HeliocentricDistanceAU)>0.0) effectivePacketSplit=true;
    if (entry.ProductCount==0 && packetCount!=1) {
      if (error!=NULL)
        *error="a product-free chemistry sink cannot use daughter packets";
      return false;
    }
    double reactantMass=0.0,productMass=0.0;
    long int concreteProductCount=0;
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
      const long int integerCoefficient=static_cast<long int>(
          std::floor(coefficient+0.5));
      if (concreteProductCount>LONG_MAX-integerCoefficient) {
        if (error!=NULL)
          *error="concrete chemistry product count overflows long int";
        return false;
      }
      concreteProductCount+=integerCoefficient;
      productMass+=coefficient*PIC::MolecularData::GetMass(term.Species);
    }

    // Repeat the input generator's K*N transaction proof at startup.  This
    // catches a stale generated packet header, a hand-edited source default,
    // or a hard-limit change made after preprocessing before any reaction can
    // partially mutate a particle list.
    if (concreteProductCount!=0 &&
        (packetCount>LONG_MAX/concreteProductCount ||
         packetCount*concreteProductCount>
             Comet::RadialSimulation::
                 MaximumPhotochemicalProductsPerReaction)) {
      if (error!=NULL)
        *error="daughter packet transaction exceeds the radial hard product limit";
      return false;
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

  if (DaughterPacketSplittingEnabled && !effectivePacketSplit) {
    if (error!=NULL)
      *error="daughter packet mode has no reaction with packet count above one";
    return false;
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
        "individual-correction=%s relative-limit=0.02 z-limit=5.0\n",
        eventsPerReaction,PIC::nTotalThreads,
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        "on"
#else
        "off"
#endif
        );

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
    const double baseSourceRateWeight=parentWeight/parentTimeStep;
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
      double parentCorrection=1.0;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      // Span the same deliberately wide Stage-4 correction range used by the
      // mover/adaptation validators.  This exercises corrected production
      // rates below and above the base W/dt carrier without enabling any
      // Stage-5 weight-window policy.
      static const double correctionCycle[]={0.0625,0.25,1.0,4.0,16.0};
      parentCorrection=correctionCycle[
          (event+static_cast<long int>(reactionIndex))%5];
#endif
      const double sourceRateWeight=
          baseSourceRateWeight*parentCorrection;
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
        double mean,bernoulliVariance,actualRate;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        // The production executor emits exactly one concrete product and
        // assigns c_d=C_parent/(W_d/dt_d).  Reconstruct that same carrier here
        // and require its represented rate to equal the corrected parent rate
        // exactly; there is no residual-allocation variance in this branch.
        const double daughterCorrection=
            sourceRateWeight/targetRateWeight;
        if (!(daughterCorrection>0.0) ||
            !std::isfinite(daughterCorrection)) {
          std::fprintf(stderr,
              "invalid daughter correction in chemistry validation\n");
          std::abort();
        }
        mean=daughterCorrection;
        bernoulliVariance=0.0;
        actualRate=targetRateWeight*daughterCorrection;
#else
        mean=GetResidualMultiplicityMean(parentWeight,
            parentTimeStep,targetWeight,targetTimeStep);
        const double residual=mean-std::floor(mean);
        bernoulliVariance=residual*(1.0-residual);
        const long int multiplicity=DrawResidualMultiplicity(parentWeight,
            parentTimeStep,targetWeight,targetTimeStep);
        actualRate=multiplicity*targetRateWeight;
#endif
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

bool Comet::Chemistry::ValidateDaughterPacketSplitting(
    long int eventsPerReaction) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  (void)eventsPerReaction;
  if (PIC::ThisThread==0)
    std::fprintf(stderr,
        "daughter-packet validation requires individual particle corrections\n");
  return false;
#else
  if (!Enabled || !DaughterPacketSplittingEnabled) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "daughter-packet validation requires enabled Stage-6 chemistry\n");
    return false;
  }
  if (eventsPerReaction<1000) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "daughter-packet validation requires at least 1000 events per reaction\n");
    return false;
  }

  // This validator deliberately exercises the same recoil generator and pure
  // correction planner used by ProcessParticle, but does not allocate real
  // particles.  It can therefore compare N=1/2/4/8 at identical physical
  // event counts without buffer pressure or transport obscuring a packet
  // error.  The separate transported science pilot remains the test of final
  // radial density statistics and Stage-5 post-window retention.
  const long int packetCounts[]={1,2,4,8};
  constexpr int packetCaseCount=4;
  const double parentVelocity[3]={431.0,-137.0,59.0};
  const double parentSpeed2=parentVelocity[0]*parentVelocity[0]+
      parentVelocity[1]*parentVelocity[1]+
      parentVelocity[2]*parentVelocity[2];
  const double launchRadius=5000.0;
  const double flightTime=1000.0;
  const double referenceRadius=launchRadius+
      std::sqrt(parentSpeed2)*flightTime;
  const double shellHalfWidth=5.0e5;
  const double shellNormalization=1.0/
      (4.0*Pi*referenceRadius*referenceRadius*Pi*shellHalfWidth);
  const ChemicalSolver::cReactionNetwork &network=GetNetwork();

  struct cEstimatorResult {
    double VelocityMean=0.0,VelocityVariance=0.0;
    double DensityMean=0.0,DensityVariance=0.0;
  };

  bool allPass=true;
  int testedReactionCount=0,sinkReactionCount=0;
  if (PIC::ThisThread==0)
    std::printf("Daughter-packet validation: events/reaction=%ld ranks=%d "
        "packets=1,2,4,8 deterministic-global-events=on\n",
        eventsPerReaction,PIC::nTotalThreads);

  for (std::size_t reactionIndex=0;
       reactionIndex<network.GetReactionCount();reactionIndex++) {
    const ChemicalSolver::cReaction &reaction=
        network.GetReaction(reactionIndex);
    if (reaction.ProductCount==0) {
      sinkReactionCount++;
      if (GetDaughterPacketCount(reactionIndex)!=1) allPass=false;
      continue;
    }

    std::vector<int> products;
    for (std::size_t term=0;term<reaction.ProductCount;term++) {
      const ChemicalSolver::cStoichiometricTerm &product=
          network.GetProductTerm(reactionIndex,term);
      const long int coefficient=static_cast<long int>(
          std::floor(product.Coefficient+0.5));
      for (long int copy=0;copy<coefficient;copy++)
        products.push_back(product.Species);
    }
    if (products.empty()) {
      allPass=false;
      continue;
    }

    std::vector<cEstimatorResult> result(packetCaseCount);
    bool reactionPass=true;
    for (int packetCase=0;packetCase<packetCaseCount;packetCase++) {
      const long int packetCount=packetCounts[packetCase];
      const std::size_t transactionCount=products.size()*packetCount;
      if (packetCount>MaximumDaughterPacketCount ||
          transactionCount>static_cast<std::size_t>(
              Comet::RadialSimulation::
                  MaximumPhotochemicalProductsPerReaction)) {
        if (PIC::ThisThread==0)
          std::printf("FAIL daughter-packet reaction=%s packets=%ld "
              "reason=validation-case-exceeds-configured-hard-limit\n",
              reaction.Id,packetCount);
        reactionPass=false;
        continue;
      }

      // The strict checker compares the same global event ensemble at 1, 4,
      // and 12 ranks.  Maintain centered online moments in extended precision
      // instead of subtracting two raw O(N*mean^2) moments.  The latter loses
      // significant digits for the ~1e-21 density proxy and makes harmless
      // MPI regrouping look like a different sample set.  These statistics
      // are validator bookkeeping only; production chemistry is unchanged.
      long double localVelocityMean=0.0L,localVelocityM2=0.0L;
      long double localDensityMean=0.0L,localDensityM2=0.0L;
      double localMaximumMomentumError=0.0;
      double localMaximumEnergyError=0.0;
      int localCorrectionFailure=0,localIndependentSample=0;
      long long localEventCount=0;

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      const double start=MPI_Wtime();
      // Strided global event IDs plus a seed derived from the ID make the
      // physical sample set invariant under 1/4/12-rank decomposition.  This
      // is a validator-only RNG discipline; production retains its configured
      // seed+rank event streams.
      for (long int event=PIC::ThisThread;event<eventsPerReaction;
           event+=PIC::nTotalThreads) {
        const unsigned long long seedValue=
            UINT64_C(104729)+UINT64_C(1000003)*reactionIndex+
            UINT64_C(10007)*packetCase+UINT64_C(7919)*event;
        rnd_seed(static_cast<int>(seedValue%UINT64_C(2000000000))+1);

        // Vary the unsplit carrier through a wide, deterministic range.  The
        // estimator below divides by it again, so this tests exact correction
        // partitioning without changing the physical packet distribution.
        static const double correctionCycle[]={0.0625,0.25,1.0,4.0,16.0};
        std::vector<double> unsplit(products.size());
        for (std::size_t product=0;product<products.size();product++)
          unsplit[product]=correctionCycle[
              (event+static_cast<long int>(product)+
               static_cast<long int>(reactionIndex))%5];

        Comet::DaughterPacketPlanner::cPlan plan;
        try {
          plan=Comet::DaughterPacketPlanner::BuildPlan(unsplit,
              static_cast<std::size_t>(packetCount),
              static_cast<std::size_t>(MaximumDaughterPacketCount),
              static_cast<std::size_t>(Comet::RadialSimulation::
                  MaximumPhotochemicalProductsPerReaction));
        }
        catch (const std::exception &exception) {
          std::fprintf(stderr,
              "daughter-packet validation planner failed for %s: %s\n",
              reaction.Id,exception.what());
          localCorrectionFailure=1;
          continue;
        }

        for (std::size_t product=0;product<products.size();product++) {
          double reconstructed=0.0;
          for (long int packet=0;packet<packetCount;packet++)
            reconstructed+=plan.Correction(packet,product);
          if (reconstructed!=unsplit[product])
            localCorrectionFailure=1;
        }

        double velocityEstimator=0.0,densityEstimator=0.0;
        double firstPacketVelocity=0.0;
        for (long int packet=0;packet<packetCount;packet++) {
          std::vector<double> velocity;
          GenerateProductVelocities(products,parentVelocity,
              reaction.ExcessEnergyJ,velocity);

          double momentum[3]={0.0,0.0,0.0};
          double energy=0.0,totalMass=0.0;
          for (std::size_t product=0;product<products.size();product++) {
            const double mass=PIC::MolecularData::GetMass(products[product]);
            const double *v=velocity.data()+3*product;
            totalMass+=mass;
            for (int idim=0;idim<3;idim++)
              momentum[idim]+=mass*v[idim];
            energy+=0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
          }
          double momentumError2=0.0;
          for (int idim=0;idim<3;idim++) {
            const double error=
                momentum[idim]-totalMass*parentVelocity[idim];
            momentumError2+=error*error;
          }
          localMaximumMomentumError=std::max(localMaximumMomentumError,
              std::sqrt(momentumError2)/std::max(
                  totalMass*std::sqrt(parentSpeed2),
                  std::numeric_limits<double>::min()));
          const double expectedEnergy=
              0.5*totalMass*parentSpeed2+reaction.ExcessEnergyJ;
          localMaximumEnergyError=std::max(localMaximumEnergyError,
              std::fabs(energy-expectedEnergy)/std::max(expectedEnergy,
                  std::numeric_limits<double>::min()));

          // Give every concrete product nonzero representation in the two
          // estimators. Velocity uses deterministic index weights 1..K; an
          // equal-weight average would collapse to the conserved center-of-
          // mass velocity for equal-mass reactions such as C2 -> C+C and
          // falsely report zero recoil variance. Density uses the same index
          // weights for a consistent all-column test. Each term also uses its
          // own exact correction share, so a lost/reordered planner column
          // changes the estimator. The radial-shell kernel remains a focused
          // density proxy, not a substitute for transported comet output.
          const double productWeightSum=0.5*products.size()*
              (products.size()+1.0);
          for (std::size_t product=0;product<products.size();product++) {
            const double *selectedVelocity=velocity.data()+3*product;
            const double share=(product+1.0)/productWeightSum*
                plan.Correction(packet,product)/unsplit[product];
            velocityEstimator+=share*selectedVelocity[0];
            const double finalX=launchRadius+
                selectedVelocity[0]*flightTime;
            const double finalY=selectedVelocity[1]*flightTime;
            const double finalZ=selectedVelocity[2]*flightTime;
            const double radius=std::sqrt(finalX*finalX+finalY*finalY+
                finalZ*finalZ);
            const double normalizedOffset=
                (radius-referenceRadius)/shellHalfWidth;
            const double shellKernel=shellNormalization/
                (1.0+normalizedOffset*normalizedOffset);
            densityEstimator+=share*shellKernel;
          }
          if (packet==0) firstPacketVelocity=velocity[0];
          else if (velocity[0]!=firstPacketVelocity)
            localIndependentSample=1;
        }

        localEventCount++;
        const long double sampleCount=localEventCount;
        const long double velocityValue=velocityEstimator;
        const long double velocityDelta=velocityValue-localVelocityMean;
        localVelocityMean+=velocityDelta/sampleCount;
        localVelocityM2+=velocityDelta*(velocityValue-localVelocityMean);
        const long double densityValue=densityEstimator;
        const long double densityDelta=densityValue-localDensityMean;
        localDensityMean+=densityDelta/sampleCount;
        localDensityM2+=densityDelta*(densityValue-localDensityMean);
      }
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      const double localWall=MPI_Wtime()-start;

      long long globalEventCount=0;
      MPI_Allreduce(&localEventCount,&globalEventCount,1,MPI_LONG_LONG,
          MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      // Gather and combine (count,mean,M2) tuples with the parallel centered-
      // moment identity.  A plain MPI_SUM cannot combine centered moments.
      // Rank-order merging on the root plus a broadcast also gives every rank
      // the identical finalized validator statistic.
      long double localMoments[4]={localVelocityMean,localVelocityM2,
          localDensityMean,localDensityM2};
      std::vector<long double> gatheredMoments;
      std::vector<long long> gatheredCounts;
      if (PIC::ThisThread==0) {
        gatheredMoments.resize(4*PIC::nTotalThreads);
        gatheredCounts.resize(PIC::nTotalThreads);
      }
      MPI_Gather(localMoments,4,MPI_LONG_DOUBLE,
          PIC::ThisThread==0 ? gatheredMoments.data() : NULL,4,
          MPI_LONG_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Gather(&localEventCount,1,MPI_LONG_LONG,
          PIC::ThisThread==0 ? gatheredCounts.data() : NULL,1,
          MPI_LONG_LONG,0,MPI_GLOBAL_COMMUNICATOR);
      long double globalMoments[4]={0.0L,0.0L,0.0L,0.0L};
      if (PIC::ThisThread==0) {
        long long accumulatedCount=0;
        for (int rank=0;rank<PIC::nTotalThreads;rank++) {
          const long long rankCount=gatheredCounts[rank];
          if (rankCount<=0) continue;
          if (accumulatedCount==0) {
            for (int field=0;field<4;field++)
              globalMoments[field]=gatheredMoments[4*rank+field];
            accumulatedCount=rankCount;
            continue;
          }
          const long double oldCount=accumulatedCount;
          const long double addedCount=rankCount;
          const long double mergedCount=oldCount+addedCount;
          for (int estimator=0;estimator<2;estimator++) {
            const int meanIndex=2*estimator;
            const int m2Index=meanIndex+1;
            const long double addedMean=
                gatheredMoments[4*rank+meanIndex];
            const long double delta=addedMean-globalMoments[meanIndex];
            globalMoments[m2Index]+=
                gatheredMoments[4*rank+m2Index]+delta*delta*
                oldCount*addedCount/mergedCount;
            globalMoments[meanIndex]+=delta*addedCount/mergedCount;
          }
          accumulatedCount+=rankCount;
        }
      }
      MPI_Bcast(globalMoments,4,MPI_LONG_DOUBLE,0,
          MPI_GLOBAL_COMMUNICATOR);
      double globalMaximumMomentumError=0.0;
      double globalMaximumEnergyError=0.0,globalWall=0.0;
      MPI_Allreduce(&localMaximumMomentumError,
          &globalMaximumMomentumError,1,MPI_DOUBLE,MPI_MAX,
          MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&localMaximumEnergyError,&globalMaximumEnergyError,1,
          MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&localWall,&globalWall,1,MPI_DOUBLE,MPI_MAX,
          MPI_GLOBAL_COMMUNICATOR);
      int globalCorrectionFailure=0,globalIndependentSample=0;
      MPI_Allreduce(&localCorrectionFailure,&globalCorrectionFailure,1,
          MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&localIndependentSample,&globalIndependentSample,1,
          MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

      if (globalEventCount!=eventsPerReaction || globalEventCount<=1) {
        if (PIC::ThisThread==0)
          std::printf("FAIL daughter-packet reaction=%s packets=%ld "
              "reason=incomplete-global-event-count actual=%lld expected=%ld\n",
              reaction.Id,packetCount,globalEventCount,eventsPerReaction);
        reactionPass=false;
        continue;
      }
      const long double denominator=static_cast<long double>(globalEventCount);
      result[packetCase].VelocityMean=static_cast<double>(globalMoments[0]);
      result[packetCase].DensityMean=static_cast<double>(globalMoments[2]);
      result[packetCase].VelocityVariance=static_cast<double>(std::max(0.0L,
          globalMoments[1]/(denominator-1.0L)));
      result[packetCase].DensityVariance=static_cast<double>(std::max(0.0L,
          globalMoments[3]/(denominator-1.0L)));

      bool packetPass=globalCorrectionFailure==0 &&
          globalMaximumMomentumError<=1.0e-11 &&
          globalMaximumEnergyError<=1.0e-11;
      if (packetCount>1 && reaction.ExcessEnergyJ>0.0)
        packetPass=globalIndependentSample!=0 && packetPass;

      if (packetCase>0) {
        const cEstimatorResult &baseline=result[0];
        const double velocityStandardError=std::sqrt(
            (baseline.VelocityVariance+result[packetCase].VelocityVariance)/
            denominator);
        const double densityStandardError=std::sqrt(
            (baseline.DensityVariance+result[packetCase].DensityVariance)/
            denominator);
        const double velocityAllowed=std::max(8.0*velocityStandardError,
            4096.0*std::numeric_limits<double>::epsilon()*
                std::max(1.0,std::fabs(baseline.VelocityMean)));
        const double densityAllowed=std::max(8.0*densityStandardError,
            4096.0*std::numeric_limits<double>::epsilon()*
                std::max(std::fabs(baseline.DensityMean),
                    std::numeric_limits<double>::min()));
        packetPass=std::fabs(result[packetCase].VelocityMean-
            baseline.VelocityMean)<=velocityAllowed && packetPass;
        packetPass=std::fabs(result[packetCase].DensityMean-
            baseline.DensityMean)<=densityAllowed && packetPass;

        // Independent packet averaging has theoretical event-estimator
        // variance Var(X)/N.  Allow generous finite-sample margins while
        // still detecting cloned velocities (scaled ratio near N) or a lost
        // correction factor (mean and variance both wrong).
        if (baseline.VelocityVariance>0.0) {
          const double scaled=result[packetCase].VelocityVariance*
              packetCount/baseline.VelocityVariance;
          packetPass=(scaled>=0.5 && scaled<=1.6) && packetPass;
        }
        if (baseline.DensityVariance>0.0) {
          const double scaled=result[packetCase].DensityVariance*
              packetCount/baseline.DensityVariance;
          packetPass=(scaled>=0.5 && scaled<=1.6) && packetPass;
        }
      }

      reactionPass=packetPass && reactionPass;
      if (PIC::ThisThread==0) {
        const double velocityScaled=(packetCase==0 ||
            result[0].VelocityVariance==0.0) ? 1.0 :
            result[packetCase].VelocityVariance*packetCount/
                result[0].VelocityVariance;
        const double densityScaled=(packetCase==0 ||
            result[0].DensityVariance==0.0) ? 1.0 :
            result[packetCase].DensityVariance*packetCount/
                result[0].DensityVariance;
        static const char *sourceNames[]={"implicit","global",
            "chemical-only","species","reaction"};
        const int policySource=GetDaughterPacketPolicySource(reactionIndex);
        const char *sourceName=(policySource>=0 && policySource<5) ?
            sourceNames[policySource] : "invalid";
        std::printf("%s daughter-packet reaction=%s packets=%ld products=%zu "
            "configured=%ld source=%d selector=%s events=%lld daughters=%llu "
            // Full binary64 round-trip precision is part of the validator
            // log interface.  The independent checker compares deterministic
            // global-event ensembles across MPI decompositions; nine printed
            // digits can straddle a decimal rounding boundary even when the
            // underlying reductions differ only at machine precision.
            "velocity-mean=%.17e velocity-variance=%.17e "
            "velocity-scaled=%.17e density-mean=%.17e "
            "density-variance=%.17e density-scaled=%.17e "
            "momentum-error=%.3e energy-error=%.3e wall=%.6g\n",
            packetPass ? "PASS" : "FAIL",reaction.Id,packetCount,
            products.size(),GetDaughterPacketCount(reactionIndex),
            policySource,sourceName,
            globalEventCount,
            static_cast<unsigned long long>(globalEventCount)*
                transactionCount,result[packetCase].VelocityMean,
            result[packetCase].VelocityVariance,velocityScaled,
            result[packetCase].DensityMean,
            result[packetCase].DensityVariance,densityScaled,
            globalMaximumMomentumError,globalMaximumEnergyError,globalWall);
      }
    }

    testedReactionCount++;
    allPass=reactionPass && allPass;
  }

  // Exercise the real particle transaction, not only its pure planner.  The
  // Rubin network's C2_2C channel is especially useful: it has a repeated
  // concrete product, the staged species policy resolves it to N=8, and its
  // only competing C2 channel is a sink.  A very large temporary parent
  // timestep makes one of those channels occur with probability effectively
  // one; deterministic seed search locates the neutral channel without
  // modifying physical rates or adding a validator-only executor hook.
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  bool productionPass=PIC::DomainBlockDecomposition::nLocalBlocks>0 &&
      PIC::DomainBlockDecomposition::BlockTable[0]!=NULL &&
      PIC::DomainBlockDecomposition::BlockTable[0]->block!=NULL;
  long int productionCount=0,parentReuseCount=0,processedCount=0;
  long int linkCheckCount=0,rollbackCount=0,seedAttempts=0;
  double maximumCorrectionError=0.0;
  double maximumRepresentedMomentumError=0.0;
  double maximumRepresentedEnergyError=0.0;
  if (productionPass) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
        PIC::DomainBlockDecomposition::BlockTable[0];
    int c2Species=-1,cSpecies=-1;
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      char symbol[_MAX_STRING_LENGTH_PIC_];
      PIC::MolecularData::GetChemSymbol(symbol,species);
      if (std::strcmp(symbol,"C2")==0) c2Species=species;
      if (std::strcmp(symbol,"C")==0) cSpecies=species;
    }
    std::size_t c2Reaction=network.GetReactionCount();
    for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++)
      if (std::strcmp(network.GetReaction(reaction).Id,"C2_2C")==0)
        c2Reaction=reaction;
    productionPass=c2Species>=0 && cSpecies>=0 &&
        c2Reaction<network.GetReactionCount() &&
        GetDaughterPacketCount(c2Reaction)==8;

    const double validationVelocity[3]={431.0,-137.0,59.0};
    const double cellWidth=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
    const double validationPosition[3]={
        node->xmin[0]+0.5*cellWidth,0.0,0.0};
    auto seedParent=[&]() -> long int {
      long int &first=node->block->FirstCellParticleTable[0];
      const long int ptr=PIC::ParticleBuffer::GetNewParticle(first);
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      PIC::ParticleBuffer::SetX(const_cast<double*>(validationPosition),data);
      PIC::ParticleBuffer::SetV(const_cast<double*>(validationVelocity),data);
      PIC::ParticleBuffer::SetI(c2Species,data);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,data);
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
      PIC::IDF::InitRotTemp(180.0,data);
      PIC::IDF::InitVibTemp(180.0,data);
#endif
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
      PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
#endif
      return ptr;
    };
    auto detachAndProcess=[&](long int ptr,int seed) -> bool {
      long int &first=node->block->FirstCellParticleTable[0];
      first=-1;
      PIC::ParticleBuffer::SetPrev(-1,ptr);
      PIC::ParticleBuffer::SetNext(-1,ptr);
      rnd_seed(seed);
      return ProcessParticle(ptr,first,node);
    };

    const double originalTimeStep=productionPass ?
        node->block->GetLocalTimeStep(c2Species) : 1.0;
    const double forcedTimeStep=1.0e12;
    if (!(originalTimeStep>0.0) || !std::isfinite(originalTimeStep))
      productionPass=false;
    int acceptedSeed=-1;
    long int acceptedParentPointer=-1;
    bool timeStepChanged=false;
    if (productionPass) {
      node->block->SetLocalTimeStep(forcedTimeStep,c2Species);
      timeStepChanged=true;
      for (int trial=0;trial<256 && acceptedSeed<0;trial++) {
        PIC::ParticleBuffer::DeleteAllParticles();
        const long int ptr=seedParent();
        BeginBatchedModel();
        const int seed=1900001+7919*trial+PIC::ThisThread;
        const bool accepted=detachAndProcess(ptr,seed);
        seedAttempts=trial+1;
        if (accepted && PIC::ParticleBuffer::NAllPart==16) {
          acceptedSeed=seed;
          acceptedParentPointer=ptr;
        }
      }
      productionPass=acceptedSeed>=0;
    }

    if (productionPass) {
      // GetNewParticle(first) head-inserts every allocated daughter.  Thus the
      // two products of each complete packet remain adjacent, in reverse
      // packet/product commit order; the reused parent is the final tail.
      std::vector<long int> particles;
      bool listLinksPass=true;
      long int previous=-1;
      for (long int ptr=node->block->FirstCellParticleTable[0];ptr!=-1;
           ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        particles.push_back(ptr);
        if (PIC::ParticleBuffer::GetPrev(ptr)!=previous)
          listLinksPass=false;
        previous=ptr;
      }
      productionCount=static_cast<long int>(particles.size());
      productionPass=productionCount==16 && productionPass;
      const double targetCarrier=node->block->GetLocalParticleWeight(cSpecies)/
          node->block->GetLocalTimeStep(cSpecies);
      const double parentCarrier=
          node->block->GetLocalParticleWeight(c2Species)/forcedTimeStep;
      const double unsplitCorrection=parentCarrier/targetCarrier;
      double correctionSum=0.0;
      const double carbonMass=PIC::MolecularData::GetMass(cSpecies);
      const double parentSpeed2=validationVelocity[0]*validationVelocity[0]+
          validationVelocity[1]*validationVelocity[1]+
          validationVelocity[2]*validationVelocity[2];
      const double packetMolecularEnergy=carbonMass*parentSpeed2+
          network.GetReaction(c2Reaction).ExcessEnergyJ;
      double firstPacketVx=0.0;
      bool independentPackets=false;
      for (std::size_t index=0;index<particles.size();index++) {
        const long int ptr=particles[index];
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        if (PIC::ParticleBuffer::GetI(data)!=cSpecies)
          productionPass=false;
        if (!IsBatchedParticleEligible(ptr)) processedCount++;
        if (ptr==acceptedParentPointer) parentReuseCount++;
        double position[3];
        PIC::ParticleBuffer::GetX(position,data);
        // GetX() only guarantees coordinates compiled into this application.
        // In the radial DIM=1 build, the inactive stack slots are deliberately
        // unspecified and therefore must not participate in the transaction
        // identity check.  The production mover and chemistry both operate on
        // exactly DIM active coordinates.
        for (int idim=0;idim<DIM;idim++)
          if (position[idim]!=validationPosition[idim]) productionPass=false;
        correctionSum+=
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
      }
      for (std::size_t packet=0;packet<particles.size()/2;packet++) {
        double representedRate[2],velocity[2][3];
        for (int product=0;product<2;product++) {
          PIC::ParticleBuffer::byte *data=
              PIC::ParticleBuffer::GetParticleDataPointer(
                  particles[2*packet+product]);
          representedRate[product]=targetCarrier*
              PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
          PIC::ParticleBuffer::GetV(velocity[product],data);
        }
        const double rateScale=std::max(representedRate[0],
            std::numeric_limits<double>::min());
        maximumCorrectionError=std::max(maximumCorrectionError,
            std::fabs(representedRate[0]-representedRate[1])/rateScale);
        double representedMomentumError2=0.0;
        for (int idim=0;idim<3;idim++) {
          const double represented=carbonMass*(
              representedRate[0]*velocity[0][idim]+
              representedRate[1]*velocity[1][idim]);
          const double expected=carbonMass*(representedRate[0]+
              representedRate[1])*validationVelocity[idim];
          const double error=represented-expected;
          representedMomentumError2+=error*error;
        }
        const double momentumScale=carbonMass*(representedRate[0]+
            representedRate[1])*std::sqrt(parentSpeed2);
        maximumRepresentedMomentumError=std::max(
            maximumRepresentedMomentumError,
            std::sqrt(representedMomentumError2)/std::max(momentumScale,
                std::numeric_limits<double>::min()));
        double representedEnergy=0.0;
        for (int product=0;product<2;product++)
          representedEnergy+=0.5*carbonMass*representedRate[product]*(
              velocity[product][0]*velocity[product][0]+
              velocity[product][1]*velocity[product][1]+
              velocity[product][2]*velocity[product][2]);
        const double representedEnergyExpected=0.5*
            (representedRate[0]+representedRate[1])*packetMolecularEnergy;
        maximumRepresentedEnergyError=std::max(
            maximumRepresentedEnergyError,
            std::fabs(representedEnergy-representedEnergyExpected)/
                std::max(representedEnergyExpected,
                    std::numeric_limits<double>::min()));
        if (packet==0) firstPacketVx=velocity[0][0];
        else if (velocity[0][0]!=firstPacketVx) independentPackets=true;
      }
      maximumCorrectionError=std::max(maximumCorrectionError,
          std::fabs(correctionSum-2.0*unsplitCorrection)/
              (2.0*unsplitCorrection));
      // Report the list invariant independently of the physical transaction
      // checks.  This keeps a future velocity/correction failure from being
      // mislabeled as a particle-list corruption in the strict log parser.
      linkCheckCount=(listLinksPass && particles.size()==16) ? 1 : 0;
      productionPass=processedCount==16 && parentReuseCount==1 &&
          linkCheckCount==1 &&
          independentPackets && maximumCorrectionError<=1.0e-12 &&
          maximumRepresentedMomentumError<=1.0e-11 &&
          maximumRepresentedEnergyError<=1.0e-11 && productionPass;

      // Repeat the identical selected event with a deliberately reduced
      // logical capacity. MaxNPart is restored before any allocation/deletion;
      // the production preflight must return false, retain the one parent,
      // and leave every particle byte unchanged. Controller-only eligible and
      // allocation-deferral diagnostics are intentionally allowed to change.
      PIC::ParticleBuffer::DeleteAllParticles();
      const long int rollbackParent=seedParent();
      BeginBatchedModel();
      std::vector<PIC::ParticleBuffer::byte> before(
          PIC::ParticleBuffer::ParticleDataLength);
      std::memcpy(before.data(),
          PIC::ParticleBuffer::GetParticleDataPointer(rollbackParent),
          PIC::ParticleBuffer::ParticleDataLength);
      const long int originalMaximum=PIC::ParticleBuffer::MaxNPart;
      PIC::ParticleBuffer::MaxNPart=16;
      const bool rollbackAccepted=
          detachAndProcess(rollbackParent,acceptedSeed);
      PIC::ParticleBuffer::MaxNPart=originalMaximum;
      const long int rollbackFirst=
          node->block->FirstCellParticleTable[0];
      const bool rollbackBytes=rollbackFirst==rollbackParent &&
          PIC::ParticleBuffer::NAllPart==1 &&
          std::memcmp(before.data(),
              PIC::ParticleBuffer::GetParticleDataPointer(rollbackParent),
              PIC::ParticleBuffer::ParticleDataLength)==0;
      if (!rollbackAccepted && rollbackBytes) rollbackCount=1;
      else productionPass=false;
    }

    if (timeStepChanged)
      node->block->SetLocalTimeStep(originalTimeStep,c2Species);
    PIC::ParticleBuffer::DeleteAllParticles();
  }

  int localProductionPass=productionPass ? 1 : 0;
  int globalProductionPass=0;
  MPI_Allreduce(&localProductionPass,&globalProductionPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  long int localProductionCounts[6]={productionCount,parentReuseCount,
      processedCount,linkCheckCount,rollbackCount,seedAttempts};
  long int globalProductionCounts[6]={};
  MPI_Allreduce(localProductionCounts,globalProductionCounts,5,MPI_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localProductionCounts+5,globalProductionCounts+5,1,MPI_LONG,
      MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  double localProductionErrors[3]={maximumCorrectionError,
      maximumRepresentedMomentumError,maximumRepresentedEnergyError};
  double globalProductionErrors[3]={};
  MPI_Allreduce(localProductionErrors,globalProductionErrors,3,MPI_DOUBLE,
      MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0)
    std::printf("%s daughter-packet production-transaction ranks=%d "
        "reaction=C2_2C packets=8 products=2 emitted=%ld reused=%ld "
        "processed=%ld link-checks=%ld capacity-rollbacks=%ld "
        "correction-error=%.3e represented-momentum=%.3e "
        "represented-energy=%.3e max-seed-attempts=%ld\n",
        globalProductionPass ? "PASS" : "FAIL",PIC::nTotalThreads,
        globalProductionCounts[0],globalProductionCounts[1],
        globalProductionCounts[2],globalProductionCounts[3],
        globalProductionCounts[4],globalProductionErrors[0],
        globalProductionErrors[1],globalProductionErrors[2],
        globalProductionCounts[5]);
  allPass=globalProductionPass!=0 && allPass;

  allPass=testedReactionCount>0 && allPass;
  int localPass=allPass ? 1 : 0,globalPass=0;
  MPI_Allreduce(&localPass,&globalPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0)
    std::printf("Daughter-packet validation result: %s "
        "neutral-reactions=%d sink-reactions=%d packet-cases=%d\n",
        globalPass ? "PASS" : "FAIL",testedReactionCount,
        sinkReactionCount,packetCaseCount);
  return globalPass!=0;
#endif
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
  const long int reactionPacketCount=GetDaughterPacketCount(selected);

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
  if (reactionPacketCount==1) {
    // Compatibility boundary: before Stage 6, one recoil solution was drawn
    // here, before correction-off residual allocation and before capacity
    // preflight. Preserve that call position exactly whenever packet mode is
    // disabled, including the historical random consumption of a transaction
    // that is later deferred. An explicit Stage-6 N=1 suppression uses this
    // same identity path. Only N>1 has a new random-stream contract and waits
    // until acceptance before drawing its independent solutions.
    GenerateProductVelocities(productSpecies,vParent,reaction.ExcessEnergyJ,
        productVelocity);
  }

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  const double parentPhysicalWeight=
      node->block->GetLocalParticleWeight(parentSpecies)*
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(parentData);
  if (!(parentPhysicalWeight>0.0) || !std::isfinite(parentPhysicalWeight)) {
    std::fprintf(stderr,
        "invalid corrected parent carrier in comet chemistry\n");
    std::abort();
  }
#elif (_SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_) || \
      (_SIMULATION_PARTICLE_WEIGHT_MODE_ == _SINGLE_GLOBAL_PARTICLE_WEIGHT_) || \
      (_SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_)
  const double parentPhysicalWeight=
      node->block->GetLocalParticleWeight(parentSpecies);
#else
#error "Comet chemistry without individual corrections requires a supported global or species-local particle weight"
#endif

  // Plan every residual multiplicity before changing the parent, accumulator,
  // diagnostics, or particle list. This makes the reaction a transaction:
  // either capacity for the complete event is available, or the untouched
  // parent and bitwise-identical phases are retried after the collective batch
  // checkpoint creates headroom.
  std::vector<long int> daughterMultiplicity(productSpecies.size(),0);
  std::vector<double> unsplitDaughterCorrection(productSpecies.size(),0.0);
  Comet::DaughterPacketPlanner::cPlan daughterPacketPlan;
  std::vector<double> daughterNextPhase(productSpecies.size(),0.0);
  std::vector<unsigned char> daughterPhaseActive(productSpecies.size(),0);
  double *residualBlockData=NULL;
  bool residualInitializationPending=false;
  std::vector<double> provisionalInitialPhases;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
  if (ResidualSamplingMode==ResidualSampling::RandomPhaseSystematic) {
    residualBlockData=GetResidualBlockData(node);
    if (residualBlockData[ResidualHeaderInitialized]==0.0) {
      // First use is planned from deterministic keyed phases without writing
      // the block.  A capacity deferral therefore leaves even a never-used
      // block bitwise unchanged; initialization is committed with the event.
      BuildInitialResidualPhases(node,provisionalInitialPhases);
      residualInitializationPending=true;
    }
    else if (residualBlockData[ResidualHeaderInitialized]!=1.0) {
      std::fprintf(stderr,"invalid systematic residual initialization flag\n");
      std::abort();
    }
  }
#endif
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
    // Every packet contains the complete concrete product list.  Store N as
    // the per-product multiplicity for diagnostics/capacity, while the pure
    // planner below divides each exact unsplit correction across packets.
    daughterMultiplicity[product]=reactionPacketCount;
    unsplitDaughterCorrection[product]=
        parentPhysicalWeight*productTimeStep/(baseWeight*timeStep);
#else
    if (ResidualSamplingMode==ResidualSampling::RandomPhaseSystematic) {
      namespace Accumulator=
          ChemicalSolver::SystematicResidualAccumulator;
      const std::size_t concrete=
          ResidualReactionOffset[selected]+product;
      if (concrete>=ResidualReactionOffset[selected+1]) {
        std::fprintf(stderr,
            "systematic residual concrete-product index is inconsistent with the generated network\n");
        std::abort();
      }
      const double mean=GetResidualMultiplicityMean(parentPhysicalWeight,
          timeStep,baseWeight,productTimeStep);
      const double phase=residualInitializationPending ?
          provisionalInitialPhases[concrete] :
          residualBlockData[ResidualRecordIndex(concrete,
              ResidualRecordPhase)];
      const Accumulator::cPlan plan=Accumulator::Plan(phase,mean);
      if (plan.Multiplicity<0 || !Accumulator::IsValidPhase(plan.NextPhase)) {
        std::fprintf(stderr,
            "invalid systematic residual plan for reaction %s product %zu\n",
            reaction.Id,product);
        std::abort();
      }
      daughterMultiplicity[product]=plan.Multiplicity;
      daughterNextPhase[product]=plan.NextPhase;
      daughterPhaseActive[product]=
          (mean!=std::floor(mean)) ? 1 : 0;
    }
    else {
      daughterMultiplicity[product]=
          DrawResidualMultiplicity(parentPhysicalWeight,timeStep,baseWeight,
              productTimeStep);
    }
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

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  try {
    daughterPacketPlan=Comet::DaughterPacketPlanner::BuildPlan(
        unsplitDaughterCorrection,
        static_cast<std::size_t>(reactionPacketCount),
        static_cast<std::size_t>(
            std::max<long int>(1,MaximumDaughterPacketCount)),
        static_cast<std::size_t>(Comet::RadialSimulation::
            MaximumPhotochemicalProductsPerReaction));
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,
        "invalid daughter-packet plan for reaction %s: %s\n",
        reaction.Id,exception.what());
    std::abort();
  }
  if (daughterPacketPlan.TotalDaughterCount!=
      static_cast<std::size_t>(totalDaughters)) {
    std::fprintf(stderr,
        "daughter-packet transaction count disagrees with chemistry allocation for reaction %s\n",
        reaction.Id);
    std::abort();
  }
#endif

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

  auto commitResidualPhases=[&]() {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
    if (ResidualSamplingMode==ResidualSampling::RandomPhaseSystematic) {
      if (residualInitializationPending) InitializeResidualBlock(node);
      for (std::size_t product=0;product<productSpecies.size();product++) {
        const std::size_t concrete=
            ResidualReactionOffset[selected]+product;
        residualBlockData[ResidualRecordIndex(concrete,
            ResidualRecordPhase)]=daughterNextPhase[product];
        if (daughterPhaseActive[product])
          residualBlockData[ResidualRecordIndex(concrete,
              ResidualRecordDiagnosticActive)]=1.0;
      }
    }
#endif
  };

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

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      // Each packet is an independently oriented reaction solution carrying
      // 1/N of the event.  Expected moments use the mathematical equal share;
      // represented moments use the exact stored correction, including the
      // final binary64 remainder.  Their totals therefore expose any carrier
      // reconstruction error instead of assuming it away.
      const double expectedPacketRate=
          parentRate/static_cast<double>(reactionPacketCount);
      const double productCarrierRate=
          node->block->GetLocalParticleWeight(species)/
          node->block->GetLocalTimeStep(species);
      for (long int packet=0;packet<reactionPacketCount;packet++) {
        const double *velocity=productVelocity.data()+3*(
            packet*productSpecies.size()+product);
        const double speedSquared=velocity[0]*velocity[0]+
            velocity[1]*velocity[1]+velocity[2]*velocity[2];
        const double representedRate=productCarrierRate*
            daughterPacketPlan.Correction(packet,product);
        term.ExpectedNumber+=expectedPacketRate;
        term.ExpectedMass+=expectedPacketRate*mass;
        term.RepresentedNumber+=representedRate;
        term.RepresentedMass+=representedRate*mass;
        for (int component=0;component<3;component++) {
          term.ExpectedMomentum[component]+=
              expectedPacketRate*mass*velocity[component];
          term.RepresentedMomentum[component]+=
              representedRate*mass*velocity[component];
        }
        term.ExpectedEnergy+=
            expectedPacketRate*0.5*mass*speedSquared;
        term.RepresentedEnergy+=
            representedRate*0.5*mass*speedSquared;
      }
#else
      const double *velocity=productVelocity.data()+3*product;
      const double speedSquared=velocity[0]*velocity[0]+
          velocity[1]*velocity[1]+velocity[2]*velocity[2];
      term.ExpectedNumber+=parentRate;
      term.ExpectedMass+=parentRate*mass;
      for (int component=0;component<3;component++)
        term.ExpectedMomentum[component]+=
            parentRate*mass*velocity[component];
      term.ExpectedEnergy+=parentRate*0.5*mass*speedSquared;
      const double productCarrierRate=
          node->block->GetLocalParticleWeight(species)/
          node->block->GetLocalTimeStep(species);
      const double representedRate=multiplicity*productCarrierRate;
      if (ResidualSamplingMode==ResidualSampling::IndependentBernoulli) {
        const double mean=GetResidualMultiplicityMean(parentPhysicalWeight,
            timeStep,node->block->GetLocalParticleWeight(species),
            node->block->GetLocalTimeStep(species));
        const double fraction=mean-std::floor(mean);
        term.NumberVariance+=productCarrierRate*productCarrierRate*
            fraction*(1.0-fraction);
      }
      term.RepresentedNumber+=representedRate;
      term.RepresentedMass+=representedRate*mass;
      for (int component=0;component<3;component++)
        term.RepresentedMomentum[component]+=
            representedRate*mass*velocity[component];
      term.RepresentedEnergy+=representedRate*0.5*mass*speedSquared;
#endif
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
    // A physical reaction with an all-zero systematic allocation still
    // consumes its cumulative residual. Commit before deleting the parent.
    commitResidualPhases();
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

  // For N>1, generate one independent, complete recoil
  // solution per accepted packet. Products inside a packet conserve momentum
  // and excess energy together; different packets use independent random
  // orientations. The N=1 identity path already generated its historical
  // solution above and must not draw it a second time.
  if (reactionPacketCount>1) {
    productVelocity.resize(3*productSpecies.size()*reactionPacketCount);
    for (long int packet=0;packet<reactionPacketCount;packet++) {
      std::vector<double> packetVelocity;
      GenerateProductVelocities(productSpecies,vParent,reaction.ExcessEnergyJ,
          packetVelocity);
      std::copy(packetVelocity.begin(),packetVelocity.end(),
          productVelocity.begin()+3*packet*productSpecies.size());
    }
  }

  commitResidualPhases();
  recordCommittedDiagnostics();

  // Clone auxiliary fields from the parent, then replace position, velocity,
  // species, and (when present) correction for each concrete daughter.  The
  // packet-major commit order keeps every independently sampled complete
  // recoil solution contiguous.  The parent slot is reused exactly once, so
  // the earlier K*N-1 capacity preflight remains exact.
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
  auto emitDaughter=[&](std::size_t product,long int packet,
      double correction) {
    const int species=productSpecies[product];
    double *velocity=productVelocity.data()+3*(
        packet*productSpecies.size()+product);
    PIC::ParticleBuffer::byte *data=
        reinterpret_cast<PIC::ParticleBuffer::byte*>(daughterTemplate);
    PIC::ParticleBuffer::SetX(xParent,data);
    PIC::ParticleBuffer::SetV(velocity,data);
    PIC::ParticleBuffer::SetI(species,data);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    if (!(correction>0.0) || !std::isfinite(correction)) {
      std::fprintf(stderr,
          "invalid daughter-packet correction in reaction %s species %d\n",
          reaction.Id,species);
      std::abort();
    }
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
        correction,data);
#else
    (void)correction;
#endif
    MarkChemistryProcessed(data);

#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    PIC::ParticleTracker::InitParticleID(daughterTemplate);
    PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(
        xParent,velocity,species,daughterTemplate,static_cast<void*>(node));
#endif
    if (!parentSlotReused) {
      // Detach the parent slot's stale links before CloneParticle, which
      // intentionally preserves destination links. Reinsert it into the
      // rebuilt cell list exactly like an allocated daughter.
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
  };

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  for (long int packet=0;packet<reactionPacketCount;packet++)
    for (std::size_t product=0;product<productSpecies.size();product++)
      emitDaughter(product,packet,
          daughterPacketPlan.Correction(packet,product));
#else
  // Correction-off residual siblings intentionally share the one historical
  // recoil solution generated before residual allocation. Stage-6 packet
  // mode is rejected in this build, so both RNG order and representation are
  // identical to the pre-Stage-6 executor.
  for (std::size_t product=0;product<productSpecies.size();product++) {
    for (long int copy=0;copy<daughterMultiplicity[product];copy++)
      emitDaughter(product,0,1.0);
  }
#endif
  return true;
}

bool Comet::Chemistry::ValidateSystematicResidualState() {
  namespace Accumulator=ChemicalSolver::SystematicResidualAccumulator;
  if (!Enabled || ResidualSamplingMode!=
      ResidualSampling::RandomPhaseSystematic || ResidualBlockDataOffset<0) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "systematic residual validation requires enabled SYSTEMATIC chemistry\n");
    return false;
  }

  BuildResidualLayout();
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> leaves;
  CollectResidualValidationLeaves(PIC::Mesh::mesh->rootTree,leaves);
  if (leaves.empty() || ResidualConcreteSpecies.empty()) return false;

  long long localOwnerBlocks=0,localAccumulators=0;
  long long localSameCarrierChecks=0,localRollbackChecks=0;
  long long localCarrierEpochChecks=0,localFirstUseRollbackChecks=0;
  long long localPartialMutations=0,localOwnerStateChecks=0;
  double localMaximumPrefixError=0.0;
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node:leaves) {
    if (node->Thread!=PIC::ThisThread) continue;
    if (node->block==NULL) return false;
    localOwnerBlocks++;
    InitializeResidualBlock(node);
    double *data=GetResidualBlockData(node);
    localOwnerStateChecks++;
    localAccumulators+=ResidualConcreteSpecies.size();

    const ChemicalSolver::cReactionNetwork &network=GetNetwork();
    for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
      const std::size_t begin=ResidualReactionOffset[reaction];
      const std::size_t end=ResidualReactionOffset[reaction+1];
      for (std::size_t product=begin;product<end;product++) {
        const double phase=data[ResidualRecordIndex(product,
            ResidualRecordPhase)];
        if (!Accumulator::IsValidPhase(phase)) return false;
        for (std::size_t earlier=begin;earlier<product;earlier++) {
          const double carrier=ResidualCarrierRate(node->block,
              ResidualConcreteSpecies[product]);
          const double earlierCarrier=ResidualCarrierRate(node->block,
              ResidualConcreteSpecies[earlier]);
          if (carrier==earlierCarrier) {
            if (phase!=data[ResidualRecordIndex(earlier,
                ResidualRecordPhase)]) return false;
            const Accumulator::cPlan first=Accumulator::Plan(phase,0.375);
            const Accumulator::cPlan second=Accumulator::Plan(
                data[ResidualRecordIndex(earlier,ResidualRecordPhase)],0.375);
            if (first.Multiplicity!=second.Multiplicity ||
                first.NextPhase!=second.NextPhase) return false;
            localSameCarrierChecks++;
            break;
          }
        }
      }
    }

    // Preview a transaction and deliberately discard it. The state must be
    // bitwise unchanged, and a retry must return the identical allocation.
    const std::size_t concrete=0;
    double &phase=data[ResidualRecordIndex(concrete,ResidualRecordPhase)];
    const double phaseBefore=phase;
    const Accumulator::cPlan preview=Accumulator::Plan(phaseBefore,0.625);
    const Accumulator::cPlan retry=Accumulator::Plan(phase,0.625);
    if (std::memcmp(&phase,&phaseBefore,sizeof(double))!=0 ||
        preview.Multiplicity!=retry.Multiplicity ||
        preview.NextPhase!=retry.NextPhase) return false;
    localRollbackChecks++;

    // Commit an all-zero allocation when possible; even K=0 must consume its
    // residual. Restore the phase afterward so migration compares the keyed
    // initialized state rather than a validation-specific perturbation.
    const double zeroMean=0.5*(1.0-phaseBefore);
    const Accumulator::cPlan zeroPlan=Accumulator::Plan(phaseBefore,zeroMean);
    if (zeroPlan.Multiplicity!=0) return false;
    phase=zeroPlan.NextPhase;
    localMaximumPrefixError=std::max(localMaximumPrefixError,
        std::fabs(zeroPlan.Multiplicity-zeroMean));
    phase=phaseBefore;

    const double epochBefore=data[ResidualHeaderCarrierEpoch];
    CloseResidualBlockDiagnosticEpoch(node,true);
    if (data[ResidualHeaderInitialized]!=0.0 ||
        data[ResidualHeaderCarrierEpoch]!=epochBefore+1.0) return false;

    // Exercise the exact first-use planning helper used by ProcessParticle.
    // It must derive common keyed phases without touching any block byte; a
    // capacity-rejected event returns at this point without a compensating
    // write or rollback operation.
    std::vector<double> uninitializedState(data,data+ResidualStateDoubleCount());
    std::vector<double> provisionalPhases;
    BuildInitialResidualPhases(node,provisionalPhases);
    if (std::memcmp(data,uninitializedState.data(),
        ResidualStateDoubleCount()*sizeof(double))!=0)
      localPartialMutations++;
    if (provisionalPhases.size()!=ResidualConcreteSpecies.size()) return false;
    localFirstUseRollbackChecks++;
    InitializeResidualBlock(node);
    if (data[ResidualHeaderInitialized]!=1.0) return false;
    localCarrierEpochChecks++;
  }

  const std::size_t stateLength=ResidualStateDoubleCount();
  if (leaves.size()*stateLength>static_cast<std::size_t>(INT_MAX))
    return false;
  std::vector<double> localExpected(leaves.size()*stateLength,0.0);
  std::vector<double> globalExpected(localExpected.size(),0.0);
  std::vector<int> oldOwner(leaves.size(),-1);
  for (std::size_t block=0;block<leaves.size();block++) {
    oldOwner[block]=leaves[block]->Thread;
    if (leaves[block]->Thread==PIC::ThisThread) {
      const double *state=GetResidualBlockData(leaves[block]);
      std::copy(state,state+stateLength,
          localExpected.begin()+block*stateLength);
    }
  }
  MPI_Allreduce(localExpected.data(),globalExpected.data(),
      static_cast<int>(globalExpected.size()),MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);

  PIC::Mesh::mesh->SetParallelLoadMeasure(Stage3MigrationLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists(
      _PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);
  PIC::DomainBlockDecomposition::UpdateBlockTable();

  long long localMigrated=0,localMigrationChecks=0;
  for (std::size_t block=0;block<leaves.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=leaves[block];
    if (node->Thread!=PIC::ThisThread) continue;
    if (node->block==NULL) return false;
    if (oldOwner[block]!=node->Thread) localMigrated++;
    const double *actual=GetResidualBlockData(node);
    const double *expected=globalExpected.data()+block*stateLength;
    if (std::memcmp(actual,expected,stateLength*sizeof(double))!=0)
      return false;
    localMigrationChecks++;
  }

  // The core restart now serializes every local-weight owner block, including
  // empty blocks. Put one controlled particle on each MPI rank so the same
  // round trip covers both kinds of record: most blocks remain empty and one
  // block per rank carries a real particle payload.  This is intentionally a
  // restart-format test, not a chemistry event; no sampling counter is changed.
  if (PIC::ParticleBuffer::NAllPart!=0) return false;
  long long localRestartParticles=0;
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node:leaves)
    if (node->Thread==PIC::ThisThread) {
      long int &first=node->block->FirstCellParticleTable[
          _getCenterNodeLocalNumber(0,0,0)];
      const long int particle=PIC::ParticleBuffer::GetNewParticle(first);
      PIC::ParticleBuffer::byte *particleData=
          PIC::ParticleBuffer::GetParticleDataPointer(particle);
      PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(
          _getCenterNodeLocalNumber(0,0,0));
      if (cell==NULL) return false;
      double position[3]={0.0,0.0,0.0};
      double velocity[3]={
          101.0+PIC::ThisThread,202.0+PIC::ThisThread,
          303.0+PIC::ThisThread};
      position[0]=cell->GetX()[0];
      PIC::ParticleBuffer::SetX(position,particleData);
      PIC::ParticleBuffer::SetV(velocity,particleData);
      PIC::ParticleBuffer::SetI(PIC::ThisThread%PIC::nTotalSpecies,
          particleData);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      // The particle checksum spans the complete particle record.  A
      // rank-distinct, non-unit correction therefore makes this restart test
      // explicitly sensitive to loss or reinitialization of c rather than
      // merely relying on the default value one.
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
          std::ldexp(1.0,(PIC::ThisThread%7)-3),particleData);
#endif
      localRestartParticles=1;
      break;
    }
  long long globalRestartParticles=0;
  MPI_Allreduce(&localRestartParticles,&globalRestartParticles,1,
      MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if (globalRestartParticles!=PIC::nTotalThreads) return false;
  const unsigned long particleChecksumBefore=
      PIC::Restart::GetParticleDataCheckSum();

  // Save nonzero phases and the controlled particles, clear both classes of
  // state, reload, and require bitwise phase and particle continuation.
  char restartName[256];
  std::snprintf(restartName,sizeof(restartName),
      "/tmp/amps-hartley2-stage3.seed=%llu.ranks=%d.restart",
      static_cast<unsigned long long>(ResidualPhaseSeed),PIC::nTotalThreads);
  PIC::Restart::SaveParticleData(restartName);
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node:leaves)
    if (node->block!=NULL) {
      // Owners are restored directly from the restart record.  Poison ghost
      // copies as well so the post-read verification proves that the core's
      // owner-to-ghost publication actually occurred.
      double *state=(node->Thread==PIC::ThisThread) ?
          GetResidualBlockData(node) : reinterpret_cast<double*>(
              node->block->GetAssociatedDataBufferPointer()+
              ResidualBlockDataOffset);
      std::fill(state,state+stateLength,-1.0);
    }
  PIC::ParticleBuffer::DeleteAllParticles();
  if (PIC::ParticleBuffer::NAllPart!=0) return false;
  PIC::Restart::ReadParticleData(restartName);
  const unsigned long particleChecksumAfter=
      PIC::Restart::GetParticleDataCheckSum();
  long long localReloadedParticles=PIC::ParticleBuffer::NAllPart;
  long long globalReloadedParticles=0;
  MPI_Allreduce(&localReloadedParticles,&globalReloadedParticles,1,
      MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  const bool particleRestartPassed=
      globalReloadedParticles==globalRestartParticles &&
      particleChecksumAfter==particleChecksumBefore;

  long long localRestartChecks=0,localRestartGhostChecks=0;
  for (std::size_t block=0;block<leaves.size();block++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=leaves[block];
    const double *expected=globalExpected.data()+block*stateLength;
    if (node->Thread==PIC::ThisThread) {
      double *state=GetResidualBlockData(node);
      if (std::memcmp(state,expected,stateLength*sizeof(double))!=0)
        return false;
      localRestartChecks++;
    }
    else if (node->block!=NULL) {
      const double *state=reinterpret_cast<const double*>(
          node->block->GetAssociatedDataBufferPointer()+
          ResidualBlockDataOffset);
      if (std::memcmp(state,expected,stateLength*sizeof(double))!=0)
        return false;
      localRestartGhostChecks++;
    }
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) std::remove(restartName);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  // Exercise the exact center-cell diagnostic collectively. The focused
  // validator has no completed physical sample, so values are intentionally
  // zero; the test proves mesh indexing, MPI reduction, schema, and one file
  // per species without pretending to provide scientific statistics.
  Comet::RadialSimulation::WriteStage3CellDiagnostics(0,1);
  // Exercise the complete systematic diagnostic serialization as well.  The
  // focused run has no sampled chemistry events, but it must still produce an
  // aligned CSV and valid JSON with null Bernoulli sigma, an explicit fixed-
  // epoch bound, and the selected uncertainty method.
  Comet::Diagnostics::WriteSamplingReport(1,0,1,
      std::vector<double>(PIC::nTotalSpecies,0.0),false);

  long long localCounts[12]={localOwnerBlocks,localAccumulators,
      localSameCarrierChecks,localRollbackChecks,localCarrierEpochChecks,
      localMigrated,localMigrationChecks,localRestartChecks,
      localOwnerStateChecks,localRestartGhostChecks,
      localFirstUseRollbackChecks,localPartialMutations};
  long long globalCounts[12]={0,0,0,0,0,0,0,0,0,0,0,0};
  MPI_Allreduce(localCounts,globalCounts,12,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  double maximumPrefixError=0.0;
  MPI_Allreduce(&localMaximumPrefixError,&maximumPrefixError,1,MPI_DOUBLE,
      MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  const bool pass=globalCounts[0]>0 && globalCounts[1]>0 &&
      globalCounts[2]>0 && globalCounts[3]==globalCounts[0] &&
      globalCounts[4]==globalCounts[0] &&
      globalCounts[6]==static_cast<long long>(leaves.size()) &&
      globalCounts[7]==static_cast<long long>(leaves.size()) &&
      globalCounts[8]==globalCounts[0] &&
      globalCounts[10]==globalCounts[0] && globalCounts[11]==0 &&
      maximumPrefixError<1.0 &&
      particleRestartPassed &&
      (PIC::nTotalThreads==1 ? globalCounts[9]==0 : globalCounts[9]>0) &&
      (PIC::nTotalThreads==1 || globalCounts[5]>0);
  if (PIC::ThisThread==0)
    std::printf(
        "%s systematic-residual validation mode=SYSTEMATIC individual-correction=%s mpi-ranks=%d owner-blocks=%lld accumulators=%lld same-carrier-packet-checks=%lld plan-retry-checks=%lld first-use-rollback-checks=%lld carrier-epoch-reset-checks=%lld migrated-nonzero-states=%lld migration-checks=%lld restart-continuation-checks=%lld restart-particles=%lld restart-particle-checksum=%d restart-ghost-sync-checks=%lld owner-only-state-checks=%lld max-prefix-error=%.17e center-cell-output=1 sampling-report-output=1 partial-mutations=%lld\n",
        pass ? "PASS" : "FAIL",
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        "on",
#else
        "off",
#endif
        PIC::nTotalThreads,globalCounts[0],
        globalCounts[1],globalCounts[2],globalCounts[3],globalCounts[10],
        globalCounts[4],
        globalCounts[5],globalCounts[6],globalCounts[7],
        globalReloadedParticles,particleRestartPassed ? 1 : 0,
        globalCounts[9],globalCounts[8],maximumPrefixError,globalCounts[11]);
  return pass;
}

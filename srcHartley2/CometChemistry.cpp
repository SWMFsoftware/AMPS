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
#include "CometChannelImportancePolicy.generated.h"
#include "CometDaughterPacketPolicy.generated.h"
#include "CometTotalReactionForcingPolicy.generated.h"
#include "CometDiagnostics.h"
#include "CometModel.h"
#include "ChannelImportancePlanner.h"
#include "DaughterPacketPlanner.h"
#include "RadialSimulation.h"
#include "TotalReactionForcingPlanner.h"

bool Comet::Chemistry::Enabled=false;
double Comet::Chemistry::HeliocentricDistanceAU=-1.0;
// Keep this definition on one physical line.  The AMPS input preprocessor
// replaces complete variable-definition lines; a wrapped initializer would
// leave its continuation behind when ChemistryResidualSampling is parsed.
int Comet::Chemistry::ResidualSamplingMode=Comet::Chemistry::ResidualSampling::IndependentBernoulli;
std::uint64_t Comet::Chemistry::ResidualPhaseSeed=1;
bool Comet::Chemistry::DaughterPacketSplittingEnabled=false;
long int Comet::Chemistry::MaximumDaughterPacketCount=1;
bool Comet::Chemistry::ChannelImportanceSamplingEnabled=false;
// Keep this definition on one physical line: the compile-time input parser
// replaces the complete line when Stage 8 is enabled for a radial case.
bool Comet::Chemistry::TotalReactionForcingEnabled=false;
int Comet::Chemistry::ActiveEstimatorMode=
    Comet::Chemistry::EstimatorMode::Configured;

static_assert(
    Comet::Chemistry::DaughterPacketGenerated::ReactionCount==
        Comet::Chemistry::Generated::ReactionCount,
    "generated daughter-packet policy does not match reaction network");
static_assert(
    Comet::Chemistry::ChannelImportanceGenerated::ReactionCount==
        Comet::Chemistry::Generated::ReactionCount,
    "generated channel-importance policy does not match reaction network");

namespace {
long int ChemistryProcessedFlagOffset=-1;
long int TotalReactionForcingStateOffset=-1;
bool LocalChemistryAllocationPressure=false;
bool ValidationSuppressTotalReactionForcing=false;

// A forced reacted carrier is still physically the parent while it travels
// to its conditional reaction point.  The mover therefore needs a compact,
// application-owned marker that survives block cloning and MPI exchange.
// Access is always through memcpy because RequestDataStorage() is byte
// granular and does not promise alignment for application records.
struct cTotalReactionForcingState {
  // Code==0 means an ordinary particle.  Otherwise the upper/lower 16 bits
  // store reaction+1 and concrete-product+1.  A zero concrete product denotes
  // a product-free neutral sink.  The +1 encoding reserves the all-zero state
  // for source particles and for cleared/reused particle-buffer slots.
  std::uint32_t Code=0;

  // Bit 31 marks the first packet's first concrete product (or the sink),
  // which owns once-per-reaction diagnostics.  The remaining 31 bits seed a
  // counter-based recoil generator; every product in one packet receives the
  // same seed and therefore reconstructs one jointly conservative solution.
  std::uint32_t KeyAndFlags=0;

  // Normalized time still to the conditional event.  It is decremented by
  // elapsed/old-local-dt and is copied unchanged across a block face, so the
  // event clock is automatically rescaled to the new block timestep.
  double RemainingFraction=0.0;
};
static_assert(sizeof(cTotalReactionForcingState)==16,
    "unexpected Stage-8 particle-state layout");

constexpr std::uint32_t TotalReactionLeaderFlag=UINT32_C(0x80000000);
constexpr std::uint32_t TotalReactionSeedMask=UINT32_C(0x7fffffff);

// Focused Stage-8 validation observes the application callback at the point
// where the *core mover* invokes it.  Keeping this audit state in the
// application, rather than adding test-only branches to pic_mover.cpp, leaves
// the generic trajectory-event contract unchanged.  Active is false in every
// science run.  The validator enables it only around one compact trajectory,
// after deleting all other particles, and reduces the rank-local counters
// immediately after the mover/exchange continuation loop has drained.
struct cTotalReactionMoverValidationEndpoint {
  long int Particle=-1;
  double Radius=0.0;
};

struct cTotalReactionMoverValidationAudit {
  bool Active=false;
  int ExpectedOwner=-1;
  double ExpectedProductEventRadius=0.0;
  double ExpectedSinkEventRadius=0.0;
  double ExpectedProductRemainingFraction=0.0;
  double ExpectedSinkRemainingFraction=0.0;
  long long ProductProcesses=0;
  long long ProductPositiveRemaining=0;
  long long SinkProcesses=0;
  long long WrongOwnerProcesses=0;
  double MaximumEventPositionError=0.0;
  double MaximumRemainingFractionError=0.0;
  std::vector<cTotalReactionMoverValidationEndpoint> ProductEndpoints;
};

cTotalReactionMoverValidationAudit TotalReactionMoverValidationAudit;

// The Stage-9 transported-daughter fixture temporarily restricts the normal
// general-NTC pair loop to the exact daughter/background pair under test.
// These selectors are never enabled by a science run: the validator saves
// and restores the production callback around one owner-local cell-kernel
// invocation.  Keeping the filter here avoids adding comet-specific state to
// the reusable collision core.
int TransportedCollisionSpecies0=-1;
int TransportedCollisionSpecies1=-1;

bool SelectTransportedCollisionPair(
    PIC::Mesh::cDataCenterNode*,int species0,int species1) {
  return (species0==TransportedCollisionSpecies0 &&
              species1==TransportedCollisionSpecies1) ||
      (species0==TransportedCollisionSpecies1 &&
              species1==TransportedCollisionSpecies0);
}

cTotalReactionForcingState GetTotalReactionState(
    const PIC::ParticleBuffer::byte *particleData) {
  if (TotalReactionForcingStateOffset<0) {
    std::fprintf(stderr,
        "total-reaction-forcing state was not reserved before particle-buffer initialization\n");
    std::abort();
  }
  cTotalReactionForcingState state;
  std::memcpy(&state,particleData+TotalReactionForcingStateOffset,
      sizeof(state));
  return state;
}

void SetTotalReactionState(PIC::ParticleBuffer::byte *particleData,
    const cTotalReactionForcingState &state) {
  if (TotalReactionForcingStateOffset<0) {
    std::fprintf(stderr,
        "total-reaction-forcing state was not reserved before particle-buffer initialization\n");
    std::abort();
  }
  std::memcpy(particleData+TotalReactionForcingStateOffset,&state,
      sizeof(state));
}

void ClearTotalReactionState(PIC::ParticleBuffer::byte *particleData) {
  SetTotalReactionState(particleData,cTotalReactionForcingState());
}

std::uint32_t EncodeTotalReactionState(std::size_t reaction,
    std::size_t concreteProduct) {
  if (reaction>=UINT16_MAX || concreteProduct>=UINT16_MAX) {
    std::fprintf(stderr,
        "Stage-8 reaction/product index exceeds the compact particle-state range\n");
    std::abort();
  }
  return (static_cast<std::uint32_t>(reaction+1)<<16) |
      static_cast<std::uint32_t>(concreteProduct+1);
}

std::uint32_t EncodeTotalReactionSink(std::size_t reaction) {
  if (reaction>=UINT16_MAX) {
    std::fprintf(stderr,
        "Stage-8 reaction index exceeds the compact particle-state range\n");
    std::abort();
  }
  return static_cast<std::uint32_t>(reaction+1)<<16;
}

std::size_t DecodeTotalReactionIndex(std::uint32_t code) {
  const std::uint32_t encoded=code>>16;
  if (encoded==0) {
    std::fprintf(stderr,"invalid empty Stage-8 reaction code\n");
    std::abort();
  }
  return static_cast<std::size_t>(encoded-1);
}

// Return -1 for a product-free sink, otherwise the zero-based concrete
// product index selected by the lower half of the compact code.
long int DecodeTotalReactionProduct(std::uint32_t code) {
  const std::uint32_t encoded=code&UINT32_C(0xffff);
  return encoded==0 ? -1 : static_cast<long int>(encoded-1);
}

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
  std::vector<double> PhysicalProbabilities;
  std::vector<double> SelectionProbabilities;
  std::vector<double> LikelihoodRatios;
  double TotalRate=0.0;
  bool ImportanceConfigured=false;
};

bool ProcessTotalReactionForcedParticle(long int ptr,
    long int &firstParticleCell,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int parentSpecies,
    const cParentReactionRateCache &parentRates,double timeStep);

std::vector<cParentReactionRateCache> ParentReactionRates;
std::vector<int> RuntimeChannelParentSpecies;
std::vector<double> RuntimePhysicalChannelProbability;
std::vector<double> RuntimeBiasedChannelProbability;
std::vector<double> RuntimeChannelLikelihoodRatio;
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
  RuntimeChannelParentSpecies.assign(network.GetReactionCount(),-1);
  RuntimePhysicalChannelProbability.assign(network.GetReactionCount(),0.0);
  RuntimeBiasedChannelProbability.assign(network.GetReactionCount(),0.0);
  RuntimeChannelLikelihoodRatio.assign(network.GetReactionCount(),1.0);
  bool foundConfiguredFamily=false;

  // Repeat the generated policy's structural checks at runtime.  This catches
  // stale/hand-edited headers before a configured row can be silently omitted
  // from the unit-unimolecular parent index.
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
    const bool configured=Comet::Chemistry::
        ChannelImportanceGenerated::Configured[reaction]!=0;
    if (network.IsUnitUnimolecular(reaction))
      RuntimeChannelParentSpecies[reaction]=
          network.GetReactantTerm(reaction,0).Species;
    if (configured && !Comet::Chemistry::ChannelImportanceSamplingEnabled) {
      if (error!=NULL)
        *error="generated channel-importance rows require enabled mode";
      return false;
    }
    if (configured && !network.IsUnitUnimolecular(reaction)) {
      if (error!=NULL)
        *error="channel importance targets a non-unit-unimolecular reaction";
      return false;
    }
  }

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
      if (!(rate>0.0) && Comet::Chemistry::
          ChannelImportanceGenerated::Configured[reaction]!=0) {
        if (error!=NULL)
          *error="channel importance targets a reaction with zero runtime rate";
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

    if (!(cache.TotalRate>0.0)) continue;
    cache.PhysicalProbabilities.resize(cache.Channels.size(),0.0);
    cache.SelectionProbabilities.resize(cache.Channels.size(),0.0);
    cache.LikelihoodRatios.resize(cache.Channels.size(),1.0);

    std::size_t configuredCount=0;
    double configuredSum=0.0;
    for (std::size_t channel=0;channel<cache.Channels.size();channel++) {
      const std::size_t reaction=cache.Channels[channel];
      RuntimeChannelParentSpecies[reaction]=species;
      const bool configured=Comet::Chemistry::
          ChannelImportanceGenerated::Configured[reaction]!=0;
      if (configured) {
        configuredCount++;
        const double probability=Comet::Chemistry::
            ChannelImportanceGenerated::BiasedProbability[reaction];
        if (!(probability>0.0) || probability>1.0 ||
            !std::isfinite(probability)) {
          if (error!=NULL)
            *error="generated channel-importance probability is invalid";
          return false;
        }
        configuredSum+=probability;
      }
    }
    if (configuredCount!=0 && configuredCount!=cache.Channels.size()) {
      if (error!=NULL)
        *error="channel importance must configure a complete positive-rate parent family";
      return false;
    }
    cache.ImportanceConfigured=configuredCount!=0;
    if (cache.ImportanceConfigured) {
      foundConfiguredFamily=true;
      if (!std::isfinite(configuredSum) ||
          std::fabs(configuredSum-1.0)>1.0e-12) {
        if (error!=NULL)
          *error="generated channel-importance family probabilities do not sum to one";
        return false;
      }
    }

    double biasedPrefix=0.0;
    for (std::size_t channel=0;channel<cache.Channels.size();channel++) {
      const std::size_t reaction=cache.Channels[channel];
      const double physical=cache.Rates[channel]/cache.TotalRate;
      // The final selection interval is the exact arithmetic remainder.  Its
      // likelihood uses that same actual interval, never a rounded input
      // value that differs from the selector's measure.
      double biased=physical;
      if (cache.ImportanceConfigured) {
        biased=(channel+1==cache.Channels.size()) ? 1.0-biasedPrefix :
            Comet::Chemistry::ChannelImportanceGenerated::
                BiasedProbability[reaction];
        biasedPrefix+=biased;
      }
      if (!(physical>0.0) || !(biased>0.0) || !std::isfinite(physical) ||
          !std::isfinite(biased)) {
        if (error!=NULL)
          *error="channel importance produced an invalid runtime probability";
        return false;
      }
      const double likelihood=physical/biased;
      if (!(likelihood>0.0) || !std::isfinite(likelihood)) {
        if (error!=NULL)
          *error="channel importance produced an invalid likelihood ratio";
        return false;
      }
      cache.PhysicalProbabilities[channel]=physical;
      cache.SelectionProbabilities[channel]=biased;
      cache.LikelihoodRatios[channel]=likelihood;
      RuntimePhysicalChannelProbability[reaction]=physical;
      RuntimeBiasedChannelProbability[reaction]=biased;
      RuntimeChannelLikelihoodRatio[reaction]=likelihood;
    }
  }

  if (Comet::Chemistry::ChannelImportanceSamplingEnabled &&
      !foundConfiguredFamily) {
    if (error!=NULL)
      *error="enabled channel importance has no configured positive-rate family";
    return false;
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
  if (TotalReactionForcingEnabled)
    PIC::ParticleBuffer::RequestDataStorage(TotalReactionForcingStateOffset,
        sizeof(cTotalReactionForcingState));
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

bool Comet::Chemistry::ConfigureEstimatorMode(int mode,std::string *error) {
  if (mode!=EstimatorMode::Configured && mode!=EstimatorMode::Analog &&
      mode!=EstimatorMode::Enhanced) {
    if (error!=NULL) *error="unknown comet chemistry estimator mode";
    return false;
  }

  if (ChemistryProcessedFlagOffset>=0 && mode!=ActiveEstimatorMode) {
    if (error!=NULL)
      *error="comet chemistry estimator mode must be selected before particle storage is reserved";
    return false;
  }

  if (mode!=EstimatorMode::Configured) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    if (error!=NULL)
      *error="Stage-9 estimator comparison requires individual particle corrections";
    return false;
#else
    // Both comparison arms intentionally require one fully enhanced binary.
    // Analog does not turn these compiled flags off: RequestParticleData and
    // ConfigureTotalReactionForcingMover therefore reserve/install the exact
    // same state layout and callback contract as Enhanced.  Only the runtime
    // estimator decisions below differ.
    if (!DaughterPacketSplittingEnabled ||
        !ChannelImportanceSamplingEnabled ||
        !TotalReactionForcingEnabled) {
      if (error!=NULL)
        *error="Stage-9 estimator comparison requires one Stage-6/7/8-enabled executable";
      return false;
    }
#endif
  }

  ActiveEstimatorMode=mode;
  return true;
}

const char* Comet::Chemistry::GetEstimatorModeName() {
  switch (ActiveEstimatorMode) {
  case EstimatorMode::Configured: return "configured";
  case EstimatorMode::Analog: return "analog";
  case EstimatorMode::Enhanced: return "enhanced";
  default: return "invalid";
  }
}

bool Comet::Chemistry::IsDaughterPacketSplittingActive() {
  return DaughterPacketSplittingEnabled &&
      ActiveEstimatorMode!=EstimatorMode::Analog;
}

bool Comet::Chemistry::IsChannelImportanceSamplingActive() {
  return ChannelImportanceSamplingEnabled &&
      ActiveEstimatorMode!=EstimatorMode::Analog;
}

bool Comet::Chemistry::IsTotalReactionForcingActive() {
  return TotalReactionForcingEnabled &&
      ActiveEstimatorMode!=EstimatorMode::Analog &&
      !ValidationSuppressTotalReactionForcing;
}

namespace {
long int GetConfiguredDaughterPacketCount(std::size_t reaction) {
  if (reaction>=Comet::Chemistry::GetNetwork().GetReactionCount()) {
    std::fprintf(stderr,"daughter-packet reaction index is out of range\n");
    std::abort();
  }
  const long int count=Comet::Chemistry::DaughterPacketSplittingEnabled ?
      Comet::Chemistry::DaughterPacketGenerated::PacketCount[reaction] : 1;
  if (count<1 ||
      (Comet::Chemistry::DaughterPacketSplittingEnabled &&
       count>Comet::Chemistry::MaximumDaughterPacketCount)) {
    std::fprintf(stderr,
        "invalid generated daughter-packet count for reaction %zu\n",
        reaction);
    std::abort();
  }
  return count;
}
}

long int Comet::Chemistry::GetDaughterPacketCount(std::size_t reaction) {
  if (reaction>=GetNetwork().GetReactionCount()) {
    std::fprintf(stderr,"daughter-packet reaction index is out of range\n");
    std::abort();
  }
  if (!IsDaughterPacketSplittingActive()) return 1;
  return GetConfiguredDaughterPacketCount(reaction);
}

int Comet::Chemistry::GetDaughterPacketPolicySource(std::size_t reaction) {
  if (reaction>=GetNetwork().GetReactionCount()) {
    std::fprintf(stderr,"daughter-packet policy index is out of range\n");
    std::abort();
  }
  return IsDaughterPacketSplittingActive() ?
      DaughterPacketGenerated::PolicySource[reaction] : 0;
}

bool Comet::Chemistry::IsChannelImportanceConfigured(std::size_t reaction) {
  if (reaction>=GetNetwork().GetReactionCount()) {
    std::fprintf(stderr,"channel-importance reaction index is out of range\n");
    std::abort();
  }
  return IsChannelImportanceSamplingActive() &&
      ChannelImportanceGenerated::Configured[reaction]!=0;
}

int Comet::Chemistry::GetChannelParentSpecies(std::size_t reaction) {
  if (reaction>=RuntimeChannelParentSpecies.size()) return -1;
  return RuntimeChannelParentSpecies[reaction];
}

double Comet::Chemistry::GetPhysicalChannelProbability(
    std::size_t reaction) {
  if (reaction>=RuntimePhysicalChannelProbability.size()) return 0.0;
  return RuntimePhysicalChannelProbability[reaction];
}

double Comet::Chemistry::GetBiasedChannelProbability(
    std::size_t reaction) {
  if (reaction>=RuntimeBiasedChannelProbability.size()) return 0.0;
  if (!IsChannelImportanceSamplingActive())
    return RuntimePhysicalChannelProbability[reaction];
  return RuntimeBiasedChannelProbability[reaction];
}

double Comet::Chemistry::GetChannelLikelihoodRatio(std::size_t reaction) {
  if (reaction>=RuntimeChannelLikelihoodRatio.size()) return 1.0;
  if (!IsChannelImportanceSamplingActive()) return 1.0;
  return RuntimeChannelLikelihoodRatio[reaction];
}

bool Comet::Chemistry::IsTotalReactionForcingConfigured(int species) {
  if (species<0 || species>=PIC::nTotalSpecies) {
    std::fprintf(stderr,
        "total-reaction-forcing species index is out of range: %d\n",species);
    std::abort();
  }
  if (!IsTotalReactionForcingActive()) return false;
  if (TotalReactionForcingGenerated::SpeciesCount!=
      static_cast<std::size_t>(PIC::nTotalSpecies)) {
    std::fprintf(stderr,
        "generated total-reaction-forcing policy has a stale species count\n");
    std::abort();
  }
  return TotalReactionForcingGenerated::Configured[species]!=0;
}

bool Comet::Chemistry::HasPendingTotalReactionEvent(
    const void *particleData) {
  // Check both the compiled mode and the reserved offset.  RequestParticleData
  // deliberately reserves no 16-byte marker in ordinary Stage-8-off builds,
  // so a checkpoint auditor must not read at an invalid negative offset.
  if (!TotalReactionForcingEnabled || TotalReactionForcingStateOffset<0 ||
      particleData==NULL)
    return false;
  const PIC::ParticleBuffer::byte *data=
      static_cast<const PIC::ParticleBuffer::byte*>(particleData);
  return GetTotalReactionState(data).Code!=0;
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
  if (ChannelImportanceSamplingEnabled) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    if (error!=NULL)
      *error="channel importance sampling requires individual particle corrections";
    return false;
#else
    if (!Enabled) {
      if (error!=NULL)
        *error="channel importance sampling requires comet chemistry";
      return false;
    }
    if (ResidualSamplingMode!=ResidualSampling::IndependentBernoulli) {
      if (error!=NULL)
        *error="channel importance sampling requires Bernoulli chemistry residual mode";
      return false;
    }
    if (!Comet::RadialSimulation::IndividualCorrectionWindowsEnabled) {
      if (error!=NULL)
        *error="channel importance sampling requires Stage-5 correction windows";
      return false;
    }
    if (!Comet::Diagnostics::Enabled) {
      if (error!=NULL)
        *error="channel importance sampling requires comet diagnostics for likelihood and ESS metadata";
      return false;
    }
#endif
  }
  bool generatedTotalReactionRow=false;
  if (TotalReactionForcingGenerated::SpeciesCount==
      static_cast<std::size_t>(PIC::nTotalSpecies))
    for (int species=0;species<PIC::nTotalSpecies;species++)
      generatedTotalReactionRow=generatedTotalReactionRow ||
          TotalReactionForcingGenerated::Configured[species]!=0;
  else if (TotalReactionForcingEnabled) {
    if (error!=NULL)
      *error="generated total-reaction-forcing policy does not match the species catalog";
    return false;
  }
  if (!TotalReactionForcingEnabled && generatedTotalReactionRow) {
    if (error!=NULL)
      *error="generated total-reaction-forcing rows require enabled mode";
    return false;
  }
  if (TotalReactionForcingEnabled) {
#if DIM != 1
    if (error!=NULL)
      *error="total-reaction forcing currently requires the DIM1 radial first-order mover";
    return false;
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    if (error!=NULL)
      *error="total-reaction forcing requires individual particle corrections";
    return false;
#else
    if (!Enabled || !Comet::RadialSimulation::Enabled) {
      if (error!=NULL)
        *error="total-reaction forcing requires radial comet chemistry";
      return false;
    }
    if (ResidualSamplingMode!=ResidualSampling::IndependentBernoulli) {
      if (error!=NULL)
        *error="total-reaction forcing requires Bernoulli residual mode";
      return false;
    }
    if (!Comet::RadialSimulation::IndividualCorrectionWindowsEnabled) {
      if (error!=NULL)
        *error="total-reaction forcing requires Stage-5 correction windows";
      return false;
    }
    if (!Comet::Diagnostics::Enabled) {
      if (error!=NULL)
        *error="total-reaction forcing requires comet diagnostics";
      return false;
    }
    if (!generatedTotalReactionRow) {
      if (error!=NULL)
        *error="enabled total-reaction forcing has no configured parent species";
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
    // Structural startup validation always checks the compiled policy, even
    // when Stage 9 selects the analog runtime arm.  Otherwise ANALOG would
    // hide a stale/unsafe N>1 row in the very executable used by ENHANCED.
    const long int packetCount=GetConfiguredDaughterPacketCount(reaction);
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
  if (!RefreshParentReactionRates(error)) return false;
  if (TotalReactionForcingEnabled) {
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      if (TotalReactionForcingGenerated::Configured[species]==0) continue;
      if (!(ParentReactionRates[species].TotalRate>0.0) ||
          ParentReactionRates[species].Channels.empty()) {
        if (error!=NULL)
          *error="total-reaction forcing targets a parent with no positive runtime reaction rate";
        return false;
      }
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

// SplitMix64 is used as a counter transform, not as mutable global state.
// Every product marker in a complete reaction packet stores the same 31-bit
// key, so it can independently reconstruct the same full K-body recoil vector
// after MPI migration.  This avoids storing 3*K doubles in every particle and
// makes the result independent of mover/list traversal order.
std::uint64_t TotalReactionCounterWord(std::uint32_t key,
    std::uint64_t counter) {
  std::uint64_t value=(static_cast<std::uint64_t>(key)<<32)^
      (counter+UINT64_C(0x9e3779b97f4a7c15));
  value=(value^(value>>30))*UINT64_C(0xbf58476d1ce4e5b9);
  value=(value^(value>>27))*UINT64_C(0x94d049bb133111eb);
  return value^(value>>31);
}

double TotalReactionUnitOpen(std::uint32_t key,std::uint64_t counter) {
  // Use the upper 53 bits and add one half ulp.  The result is strictly inside
  // (0,1), which makes Box--Muller finite without a rejection loop or mutable
  // RNG state.
  const std::uint64_t mantissa=TotalReactionCounterWord(key,counter)>>11;
  return (static_cast<double>(mantissa)+0.5)/9007199254740992.0;
}

double TotalReactionStandardNormal(std::uint32_t key,
    std::uint64_t component) {
  const std::uint64_t pair=component/2;
  const double u1=TotalReactionUnitOpen(key,2*pair);
  const double u2=TotalReactionUnitOpen(key,2*pair+1);
  const double radius=std::sqrt(-2.0*std::log(u1));
  const double angle=2.0*Pi*u2;
  return (component%2==0) ? radius*std::cos(angle) :
      radius*std::sin(angle);
}

void GenerateProductVelocitiesDeterministic(
    const std::vector<int> &species,const double parentVelocity[3],
    double excessEnergyJ,std::uint32_t key,std::vector<double> &velocity) {
  const std::size_t count=species.size();
  velocity.assign(3*count,0.0);
  if (count==0) return;
  if (count==1) {
    if (excessEnergyJ>0.0) {
      std::fprintf(stderr,
          "one-product Stage-8 reaction has nonzero excess energy\n");
      std::abort();
    }
    for (int dimension=0;dimension<3;dimension++)
      velocity[dimension]=parentVelocity[dimension];
    return;
  }

  double totalMass=0.0,centerOfMass[3]={0.0,0.0,0.0};
  for (std::size_t product=0;product<count;product++) {
    const double mass=PIC::MolecularData::GetMass(species[product]);
    totalMass+=mass;
    for (int dimension=0;dimension<3;dimension++) {
      const double value=TotalReactionStandardNormal(key,
          3*product+dimension);
      velocity[3*product+dimension]=value;
      centerOfMass[dimension]+=mass*value;
    }
  }
  for (int dimension=0;dimension<3;dimension++)
    centerOfMass[dimension]/=totalMass;

  double sampledEnergy=0.0;
  for (std::size_t product=0;product<count;product++) {
    const double mass=PIC::MolecularData::GetMass(species[product]);
    for (int dimension=0;dimension<3;dimension++)
      velocity[3*product+dimension]-=centerOfMass[dimension];
    const double *relative=velocity.data()+3*product;
    sampledEnergy+=0.5*mass*(relative[0]*relative[0]+
        relative[1]*relative[1]+relative[2]*relative[2]);
  }
  if (excessEnergyJ>0.0 && !(sampledEnergy>0.0)) {
    std::fprintf(stderr,"degenerate deterministic Stage-8 recoil sample\n");
    std::abort();
  }
  const double scale=excessEnergyJ>0.0 ?
      std::sqrt(excessEnergyJ/sampledEnergy) : 0.0;
  for (std::size_t product=0;product<count;product++)
    for (int dimension=0;dimension<3;dimension++)
      velocity[3*product+dimension]=parentVelocity[dimension]+
          scale*velocity[3*product+dimension];
}

void ExpandConcreteProducts(std::size_t reactionIndex,
    std::vector<int> &species,std::vector<std::size_t> *termIndices=NULL) {
  species.clear();
  if (termIndices!=NULL) termIndices->clear();
  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  const ChemicalSolver::cReaction &reaction=
      network.GetReaction(reactionIndex);
  for (std::size_t term=0;term<reaction.ProductCount;term++) {
    const ChemicalSolver::cStoichiometricTerm &entry=
        network.GetProductTerm(reactionIndex,term);
    const long int copies=static_cast<long int>(
        std::floor(entry.Coefficient+0.5));
    for (long int copy=0;copy<copies;copy++) {
      species.push_back(entry.Species);
      if (termIndices!=NULL) termIndices->push_back(term);
    }
  }
}

int QueryTotalReactionEvent(long int ptr,int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,double localTimeStep,
    double maximumSegment,double &eventTime) {
  (void)node;
  const PIC::ParticleBuffer::byte *data=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  const cTotalReactionForcingState state=GetTotalReactionState(data);
  if (state.Code==0) return PIC::Mover::TrajectoryEventNone;
  const std::size_t reaction=DecodeTotalReactionIndex(state.Code);
  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  if (reaction>=network.GetReactionCount() ||
      network.GetReactantTerm(reaction,0).Species!=species ||
      !std::isfinite(state.RemainingFraction) ||
      state.RemainingFraction<0.0 || state.RemainingFraction>1.0 ||
      !(localTimeStep>0.0) || !std::isfinite(localTimeStep) ||
      maximumSegment<0.0 || !std::isfinite(maximumSegment)) {
    std::fprintf(stderr,"invalid pending Stage-8 trajectory event\n");
    std::abort();
  }
  const double candidate=state.RemainingFraction*localTimeStep;
  const double tolerance=64.0*std::numeric_limits<double>::epsilon()*
      std::max(localTimeStep,maximumSegment);
  if (candidate>maximumSegment+tolerance)
    return PIC::Mover::TrajectoryEventNone;
  eventTime=std::min(candidate,maximumSegment);
  return PIC::Mover::TrajectoryEventPending;
}

void AdvanceTotalReactionEvent(long int ptr,int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,double localTimeStep,
    double elapsedTime) {
  (void)node;
  PIC::ParticleBuffer::byte *data=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  cTotalReactionForcingState state=GetTotalReactionState(data);
  if (state.Code==0) return;
  const std::size_t reaction=DecodeTotalReactionIndex(state.Code);
  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  if (reaction>=network.GetReactionCount() ||
      network.GetReactantTerm(reaction,0).Species!=species ||
      !(localTimeStep>0.0) || !std::isfinite(localTimeStep) ||
      elapsedTime<0.0 || !std::isfinite(elapsedTime)) {
    std::fprintf(stderr,"invalid Stage-8 trajectory-event advance\n");
    std::abort();
  }
  state.RemainingFraction-=elapsedTime/localTimeStep;
  const double tolerance=128.0*std::numeric_limits<double>::epsilon();
  if (std::fabs(state.RemainingFraction)<=tolerance)
    state.RemainingFraction=0.0;
  if (state.RemainingFraction<0.0) {
    if (state.RemainingFraction < -tolerance) {
      std::fprintf(stderr,"Stage-8 event clock advanced past zero\n");
      std::abort();
    }
    state.RemainingFraction=0.0;
  }
  SetTotalReactionState(data,state);
}

int ProcessTotalReactionEvent(long int ptr,double *position,double *velocity,
    int &species,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,
    double oldLocalTimeStep,double &remainingTime) {
  PIC::ParticleBuffer::byte *data=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  const cTotalReactionForcingState state=GetTotalReactionState(data);
  if (state.Code==0 || state.RemainingFraction!=0.0) {
    std::fprintf(stderr,"Stage-8 processor received an event that is not due\n");
    std::abort();
  }
  const std::size_t reactionIndex=DecodeTotalReactionIndex(state.Code);
  const long int concreteProduct=DecodeTotalReactionProduct(state.Code);
  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();
  if (reactionIndex>=network.GetReactionCount() ||
      network.GetReactantTerm(reactionIndex,0).Species!=species) {
    std::fprintf(stderr,"Stage-8 processor reaction/parent mismatch\n");
    std::abort();
  }
  const int parentSpecies=species;
  const double parentWeight=node->block->GetLocalParticleWeight(parentSpecies);
  const double parentTimeStep=node->block->GetLocalTimeStep(parentSpecies);
  const double parentCorrection=
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
  if (!(parentWeight>0.0) || !(parentTimeStep>0.0) ||
      !(parentCorrection>0.0) || !std::isfinite(parentCorrection)) {
    std::fprintf(stderr,"invalid Stage-8 pending parent carrier\n");
    std::abort();
  }

  std::vector<int> products;
  std::vector<std::size_t> terms;
  ExpandConcreteProducts(reactionIndex,products,&terms);
  const ChemicalSolver::cReaction &reaction=
      network.GetReaction(reactionIndex);
  const bool leader=(state.KeyAndFlags&TotalReactionLeaderFlag)!=0 &&
      concreteProduct<=0;
  const std::uint32_t recoilKey=state.KeyAndFlags&TotalReactionSeedMask;
  const long int packetCount=
      Comet::Chemistry::GetDaughterPacketCount(reactionIndex);
  const double pendingParentRate=parentWeight/parentTimeStep*parentCorrection;
  const double likelihood=
      Comet::Chemistry::GetChannelLikelihoodRatio(reactionIndex);

  // The validation-only audit is deliberately recorded inside the callback:
  // a test that called Query/Advance/Process directly could reproduce the
  // chemistry algebra but could not prove where the core trajectory mover
  // stopped.  In particular, the sink fixture is due exactly at a block face.
  // Seeing it at the destination radius and owner proves that the face won
  // the tie, Advance committed a zero clock before the MPI copy, and the
  // continuation mover consumed that due-zero event exactly once.
  if (TotalReactionMoverValidationAudit.Active) {
    const bool sink=concreteProduct<0;
    const double radius=std::sqrt(position[0]*position[0]+position[1]*position[1]+
        position[2]*position[2]);
    const double expectedRadius=sink ?
        TotalReactionMoverValidationAudit.ExpectedSinkEventRadius :
        TotalReactionMoverValidationAudit.ExpectedProductEventRadius;
    const double expectedRemaining=sink ?
        TotalReactionMoverValidationAudit.ExpectedSinkRemainingFraction :
        TotalReactionMoverValidationAudit.ExpectedProductRemainingFraction;
    TotalReactionMoverValidationAudit.MaximumEventPositionError=std::max(
        TotalReactionMoverValidationAudit.MaximumEventPositionError,
        std::fabs(radius-expectedRadius));
    TotalReactionMoverValidationAudit.MaximumRemainingFractionError=std::max(
        TotalReactionMoverValidationAudit.MaximumRemainingFractionError,
        std::fabs(remainingTime/oldLocalTimeStep-expectedRemaining));
    if (PIC::ThisThread!=TotalReactionMoverValidationAudit.ExpectedOwner ||
        node==NULL || node->Thread!=
            TotalReactionMoverValidationAudit.ExpectedOwner)
      TotalReactionMoverValidationAudit.WrongOwnerProcesses++;
    if (sink) TotalReactionMoverValidationAudit.SinkProcesses++;
    else TotalReactionMoverValidationAudit.ProductProcesses++;
  }

  // Clear before returning either result.  The core immediately performs a
  // zero-length re-query and rejects a marker that was not consumed exactly
  // once; clearing first also makes a product continuation ordinary chemistry
  // state for all later faces and MPI exchanges.
  ClearTotalReactionState(data);

  if (concreteProduct<0) {
    if (!products.empty() || packetCount!=1) {
      std::fprintf(stderr,"invalid Stage-8 sink marker\n");
      std::abort();
    }
    if (Comet::Diagnostics::Enabled && leader) {
      Comet::Diagnostics::RecordReactionSelection(reactionIndex,true,
          likelihood,pendingParentRate*likelihood);
      Comet::Diagnostics::RecordChemicalDestruction(parentSpecies,
          pendingParentRate);
    }
    return PIC::Mover::TrajectoryEventDelete;
  }
  if (concreteProduct>=static_cast<long int>(products.size()) ||
      products.empty()) {
    std::fprintf(stderr,"Stage-8 concrete product index is invalid\n");
    std::abort();
  }

  std::vector<double> productVelocities;
  GenerateProductVelocitiesDeterministic(products,velocity,
      reaction.ExcessEnergyJ,recoilKey,productVelocities);
  const int productSpecies=products[concreteProduct];
  const double productWeight=node->block->GetLocalParticleWeight(productSpecies);
  const double productTimeStep=node->block->GetLocalTimeStep(productSpecies);
  if (!(productWeight>0.0) || !(productTimeStep>0.0) ||
      !std::isfinite(productWeight) || !std::isfinite(productTimeStep) ||
      !(oldLocalTimeStep>0.0) || !std::isfinite(oldLocalTimeStep)) {
    std::fprintf(stderr,"invalid Stage-8 product carrier\n");
    std::abort();
  }

  // The K pending parent markers partition one packet's reacted carrier.
  // Multiplication by K restores the complete reaction yield for each
  // concrete product.  L=q/q* applies Stage-7 channel importance, and the
  // parent/product W/dt ratio converts physical rate without any boundary
  // cloning bias.  Stage-6 already divided the pending correction among N
  // packet keys, so no second division is applied here.
  const double productCorrection=parentCorrection*
      static_cast<double>(products.size())*likelihood*
      (parentWeight/parentTimeStep)/(productWeight/productTimeStep);
  if (!(productCorrection>0.0) || !std::isfinite(productCorrection)) {
    std::fprintf(stderr,"invalid Stage-8 converted product correction\n");
    std::abort();
  }
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(productCorrection,
      data);
  PIC::ParticleBuffer::SetI(productSpecies,data);
  species=productSpecies;
  for (int dimension=0;dimension<3;dimension++)
    velocity[dimension]=productVelocities[3*concreteProduct+dimension];

  // remainingTime is expressed in seconds of the old parent timestep.  The
  // logical fraction still to travel is invariant, so converting species
  // requires dt_new/dt_old before the mover recomputes geometry and continues.
  remainingTime*=productTimeStep/oldLocalTimeStep;
  if (remainingTime<0.0 || !std::isfinite(remainingTime)) {
    std::fprintf(stderr,"invalid Stage-8 product continuation interval\n");
    std::abort();
  }

  if (TotalReactionMoverValidationAudit.Active) {
    // The production mover will integrate this exact post-conversion velocity
    // for remainingTime.  The compact fixture chooses a product timestep that
    // keeps the endpoint inside the destination block, so the stored radius is
    // an independent, exact check that the event did not truncate the unused
    // part of the logical local timestep.
    double endpoint[3];
    for (int dimension=0;dimension<3;dimension++)
      endpoint[dimension]=position[dimension]+
          remainingTime*velocity[dimension];
    TotalReactionMoverValidationAudit.ProductEndpoints.push_back(
        {ptr,std::sqrt(endpoint[0]*endpoint[0]+endpoint[1]*endpoint[1]+
             endpoint[2]*endpoint[2])});
    if (remainingTime>0.0)
      TotalReactionMoverValidationAudit.ProductPositiveRemaining++;
  }

  if (Comet::Diagnostics::Enabled) {
    const double productRate=productWeight/productTimeStep*productCorrection;
    const double mass=PIC::MolecularData::GetMass(productSpecies);
    const double speedSquared=velocity[0]*velocity[0]+
        velocity[1]*velocity[1]+velocity[2]*velocity[2];
    double momentum[3];
    for (int dimension=0;dimension<3;dimension++)
      momentum[dimension]=productRate*mass*velocity[dimension];
    const double energy=productRate*0.5*mass*speedSquared;
    Comet::Diagnostics::RecordChemistryProduct(reactionIndex,
        terms[concreteProduct],1,1,0,productRate,productRate*mass,momentum,
        energy,productRate,productRate*mass,momentum,energy,0.0);
    Comet::Diagnostics::RecordChemicalProduction(productSpecies,productRate,
        productRate,0.0);
    if (leader) {
      const double fullReactedRate=pendingParentRate*
          static_cast<double>(products.size())*packetCount;
      Comet::Diagnostics::RecordReactionSelection(reactionIndex,false,
          likelihood,fullReactedRate*likelihood);
      Comet::Diagnostics::RecordChemicalDestruction(parentSpecies,
          fullReactedRate);
    }
  }
  return PIC::Mover::TrajectoryEventContinue;
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

void Comet::Chemistry::ConfigureTotalReactionForcingMover() {
  if (!TotalReactionForcingEnabled) return;
  if (TotalReactionForcingStateOffset<0) {
    std::fprintf(stderr,
        "cannot install Stage-8 mover before reserving particle state\n");
    std::abort();
  }
  // Install the contract atomically from the application's serial
  // initialization path. PIC::Mover::Init independently rejects a partial
  // set, which protects future callers that bypass this helper.
  PIC::Mover::QueryTrajectoryEvent=QueryTotalReactionEvent;
  PIC::Mover::AdvanceTrajectoryEvent=AdvanceTotalReactionEvent;
  PIC::Mover::ProcessTrajectoryEvent=ProcessTotalReactionEvent;
}

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

  // A chemistry case is not invalid merely because every reacting species
  // has the same configured W/dt carrier rate.  This occurs naturally in the
  // I3 H2O-family network: H2O is the largest direct-source carrier and all
  // of its descendants inherit that rate, so every residual multiplicity is
  // exactly one.  Individual reactions above still exercise stoichiometry,
  // recoil, mass, momentum, and energy.  Keep reporting the unequal count so
  // dedicated importance/roulette fixtures (including the Halley Stage-4
  // archive checker) can and do require a positive value, but do not make
  // that fixture-specific coverage condition a universal chemistry-validity
  // requirement.
  allPass=(testedReactionCount>0) && allPass;
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
      // This is explicitly the accepted Stage-6 endpoint-transaction fixture.
      // A Stage-8-enabled executable must replay it unchanged before the
      // separate Stage-8 validator tests in-flight conversion. Temporarily
      // suppress only the occurrence-forcing dispatch around this one direct
      // call; generated packet and Stage-7 channel policies remain active.
      const bool savedTotalReactionSuppression=
          ValidationSuppressTotalReactionForcing;
      ValidationSuppressTotalReactionForcing=true;
      const bool result=ProcessParticle(ptr,first,node);
      ValidationSuppressTotalReactionForcing=
          savedTotalReactionSuppression;
      return result;
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
      // When Stage 7 is active the production executor must apply L before
      // Stage 6 partitions the product across packets.  Reusing the real
      // transaction check here proves the composed correction L/N, including
      // its final arithmetic remainder, rather than testing either feature in
      // isolation.
      const double unsplitCorrection=parentCarrier/targetCarrier*
          GetChannelLikelihoodRatio(c2Reaction);
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

bool Comet::Chemistry::ValidateChannelImportanceSampling(
    long int eventsPerFamily) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  (void)eventsPerFamily;
  if (PIC::ThisThread==0)
    std::fprintf(stderr,
        "channel-importance validation requires individual corrections\n");
  return false;
#else
  namespace Accumulator=ChemicalSolver::SystematicResidualAccumulator;
  if (!Enabled || !ChannelImportanceSamplingEnabled) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "channel-importance validation requires enabled Stage-7 chemistry\n");
    return false;
  }
  if (eventsPerFamily<1000 || !ParentReactionRatesReady) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "channel-importance validation requires at least 1000 events and initialized rates\n");
    return false;
  }

  const ChemicalSolver::cReactionNetwork &network=GetNetwork();
  auto uniformFor=[](std::uint64_t event,std::uint64_t stream,
      std::uint64_t parent) {
    std::uint64_t key=Accumulator::Combine(event,stream);
    key=Accumulator::Combine(key,parent);
    const std::uint64_t bits=Accumulator::Mix(key);
    return static_cast<double>(bits>>11)*0x1.0p-53;
  };
  auto select=[](const std::vector<double> &probability,double uniform) {
    std::size_t selected=probability.size()-1;
    double remainder=uniform;
    for (std::size_t channel=0;channel<probability.size();channel++) {
      remainder-=probability[channel];
      if (remainder<=0.0) {
        selected=channel;
        break;
      }
    }
    return selected;
  };
  auto normalizedError=[](long double actual,long double expected,
      long double variance,long long samples) {
    const long double sigma=std::sqrt(std::max(0.0L,
        static_cast<long double>(samples)*variance));
    const long double floor=4096.0L*std::numeric_limits<double>::epsilon()*
        std::max(1.0L,std::fabs(expected));
    return static_cast<double>(std::fabs(actual-expected)/
        std::max(8.0L*sigma,floor));
  };

  if (PIC::ThisThread==0)
    std::printf("Channel-importance validation: events/family=%ld ranks=%d "
        "deterministic-global-events=on\n",eventsPerFamily,
        PIC::nTotalThreads);

  bool allPass=true;
  int familyCount=0,biasedFamilyCount=0,reactionCount=0,sinkCount=0;
  bool coveredMany=false,coveredRare=false,coveredReducedRare=false;
  bool coveredPacket=false;
  for (int parent=0;parent<PIC::nTotalSpecies;parent++) {
    const cParentReactionRateCache &cache=ParentReactionRates[parent];
    if (!(cache.TotalRate>0.0) || cache.Channels.empty()) continue;
    familyCount++;
    if (cache.ImportanceConfigured) biasedFamilyCount++;
    if (cache.Channels.size()>=3 && cache.ImportanceConfigured)
      coveredMany=true;

    const std::size_t channels=cache.Channels.size();
    std::vector<unsigned long long> localCount(2*channels,0),
        globalCount(2*channels,0);
    long long localReacted=0,globalReacted=0;
    // Derive a nontrivial survival test from the actual total rate.  Both arms
    // reuse the same occurrence uniform and therefore must have exactly the
    // same reacted/surviving parent set; only their conditional channel draw
    // is allowed to differ.
    const double validationTimeStep=-std::log(0.375)/cache.TotalRate;
    const double reactionProbability=
        1.0-std::exp(-cache.TotalRate*validationTimeStep);
    for (long int event=PIC::ThisThread;event<eventsPerFamily;
         event+=PIC::nTotalThreads) {
      if (uniformFor(static_cast<std::uint64_t>(event),0,
          static_cast<std::uint64_t>(parent))>=reactionProbability) continue;
      localReacted++;
      const double channelUniform=uniformFor(
          static_cast<std::uint64_t>(event),1,
          static_cast<std::uint64_t>(parent));
      const std::size_t analog=select(cache.PhysicalProbabilities,
          channelUniform);
      const std::size_t biased=select(cache.SelectionProbabilities,
          channelUniform);
      localCount[analog]++;
      localCount[channels+biased]++;

      // Drive the real production diagnostic accumulator with the same
      // deterministic local selections used by this validator.  Every MPI
      // rank records only its event-ID partition; WriteSamplingReport below
      // performs the ordinary collective reduction.  This closes the durable
      // CSV/JSON metadata path in addition to checking the estimator printed
      // in this log.  cache.TotalRate is used as a finite synthetic parent
      // carrier rate, so the stored weighted-channel rate obeys the exact
      // production definition parentRate*L without requiring a model particle.
      const std::size_t selectedReaction=cache.Channels[biased];
      Comet::Diagnostics::RecordReactionSelection(selectedReaction,
          network.GetReaction(selectedReaction).ProductCount==0,
          cache.LikelihoodRatios[biased],
          cache.TotalRate*cache.LikelihoodRatios[biased]);
    }
    MPI_Allreduce(localCount.data(),globalCount.data(),
        static_cast<int>(globalCount.size()),MPI_UNSIGNED_LONG_LONG,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&localReacted,&globalReacted,1,MPI_LONG_LONG,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);

    bool familyPass=globalReacted>0;
    long double familyWeight=0.0L,familyWeight2=0.0L;
    long double exactMoment[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    long double analogMoment[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    long double biasedMoment[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    long double analogFirst[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    long double analogSecond[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    long double biasedFirst[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    long double biasedSecond[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    const double parentVelocity[3]={400.0,-120.0,55.0};
    const double parentMass=PIC::MolecularData::GetMass(parent);
    const double parentSpeed2=parentVelocity[0]*parentVelocity[0]+
        parentVelocity[1]*parentVelocity[1]+
        parentVelocity[2]*parentVelocity[2];

    for (std::size_t channel=0;channel<channels;channel++) {
      const std::size_t reactionIndex=cache.Channels[channel];
      const ChemicalSolver::cReaction &reaction=
          network.GetReaction(reactionIndex);
      const double physical=cache.PhysicalProbabilities[channel];
      const double biased=cache.SelectionProbabilities[channel];
      const double likelihood=cache.LikelihoodRatios[channel];
      const unsigned long long analogSelections=globalCount[channel];
      const unsigned long long biasedSelections=globalCount[channels+channel];
      const double weighted=biasedSelections*likelihood;
      const double expected=globalReacted*physical;
      const double selectionSigma=std::sqrt(globalReacted*biased*
          (1.0-biased));
      const double weightedVariance=globalReacted*physical*physical*
          (1.0/biased-1.0);
      const double weightedSigma=std::sqrt(std::max(0.0,weightedVariance));
      const double selectionAllowed=std::max(8.0*selectionSigma,
          4096.0*std::numeric_limits<double>::epsilon()*
              std::max(1.0,globalReacted*biased));
      const double weightedAllowed=std::max(8.0*weightedSigma,
          4096.0*std::numeric_limits<double>::epsilon()*
              std::max(1.0,expected));
      bool reactionPass=std::fabs(static_cast<double>(biasedSelections)-
          globalReacted*biased)<=selectionAllowed &&
          std::fabs(weighted-expected)<=weightedAllowed &&
          likelihood==physical/biased && physical>0.0 && biased>0.0;
      const double analogVariance=physical*(1.0-physical);
      const double biasedVariance=physical*physical*(1.0/biased-1.0);
      const double varianceRatio=analogVariance>0.0 ?
          biasedVariance/analogVariance : 1.0;
      if (physical<0.05 && cache.ImportanceConfigured) {
        coveredRare=true;
        if (biased>physical && varianceRatio<1.0)
          coveredReducedRare=true;
      }
      if (GetDaughterPacketCount(reactionIndex)>1 &&
          cache.ImportanceConfigured) coveredPacket=true;
      if (reaction.ProductCount==0) sinkCount++;
      reactionCount++;
      familyWeight+=biasedSelections*likelihood;
      familyWeight2+=biasedSelections*likelihood*likelihood;

      long double observable[6]={0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
      if (reaction.ProductCount!=0) {
        long double concrete=0.0L,trackedMass=0.0L;
        for (std::size_t term=0;term<reaction.ProductCount;term++) {
          const ChemicalSolver::cStoichiometricTerm &product=
              network.GetProductTerm(reactionIndex,term);
          concrete+=product.Coefficient;
          trackedMass+=product.Coefficient*
              PIC::MolecularData::GetMass(product.Species);
        }
        observable[0]=concrete;
        observable[1]=trackedMass;
        observable[2]=parentMass*parentVelocity[0];
        observable[3]=parentMass*parentVelocity[1];
        observable[4]=parentMass*parentVelocity[2];
        observable[5]=0.5L*parentMass*parentSpeed2+
            reaction.ExcessEnergyJ;
      }
      for (int moment=0;moment<6;moment++) {
        const long double value=observable[moment];
        exactMoment[moment]+=globalReacted*physical*value;
        analogMoment[moment]+=analogSelections*value;
        biasedMoment[moment]+=biasedSelections*likelihood*value;
        analogFirst[moment]+=physical*value;
        analogSecond[moment]+=physical*value*value;
        biasedFirst[moment]+=physical*value;
        biasedSecond[moment]+=biased*likelihood*likelihood*value*value;
      }

      if (PIC::ThisThread==0)
        std::printf("%s channel-importance reaction=%s parent=%d "
            "physical=%.17e biased=%.17e likelihood=%.17e analog=%llu "
            "selected=%llu weighted=%.17e expected=%.17e sigma=%.17e "
            "variance-ratio=%.17e sink=%d packets=%ld\n",
            reactionPass ? "PASS" : "FAIL",reaction.Id,parent,physical,
            biased,likelihood,analogSelections,biasedSelections,weighted,
            expected,weightedSigma,varianceRatio,
            reaction.ProductCount==0 ? 1 : 0,
            GetDaughterPacketCount(reactionIndex));
      familyPass=reactionPass && familyPass;
    }

    const double familyEss=familyWeight2>0.0L ?
        static_cast<double>(familyWeight*familyWeight/familyWeight2) : 0.0;
    const double familyEssFraction=globalReacted>0 ?
        familyEss/globalReacted : 0.0;
    familyPass=familyEss>=1.0 &&
        familyEss<=globalReacted*(1.0+1.0e-12) && familyPass;

    double maximumAnalogZ=0.0,maximumBiasedZ=0.0;
    for (int moment=0;moment<6;moment++) {
      const long double analogVariance=std::max(0.0L,
          analogSecond[moment]-analogFirst[moment]*analogFirst[moment]);
      const long double biasedVariance=std::max(0.0L,
          biasedSecond[moment]-biasedFirst[moment]*biasedFirst[moment]);
      maximumAnalogZ=std::max(maximumAnalogZ,
          normalizedError(analogMoment[moment],exactMoment[moment],
              analogVariance,globalReacted));
      maximumBiasedZ=std::max(maximumBiasedZ,
          normalizedError(biasedMoment[moment],exactMoment[moment],
              biasedVariance,globalReacted));
    }
    familyPass=maximumAnalogZ<=1.0 && maximumBiasedZ<=1.0 && familyPass;
    if (PIC::ThisThread==0) {
      std::printf("%s channel-importance family parent=%d channels=%zu "
          "total-rate=%.17e trials=%ld reacted=%lld analog-loss=%lld "
          "biased-loss=%lld family-ess=%.17e ess-fraction=%.17e "
          "configured=%d\n",familyPass ? "PASS" : "FAIL",parent,
          channels,cache.TotalRate,eventsPerFamily,globalReacted,
          globalReacted,globalReacted,familyEss,familyEssFraction,
          cache.ImportanceConfigured ? 1 : 0);
      std::printf("%s channel-importance moments parent=%d "
          "analog-max-z=%.17e biased-max-z=%.17e\n",
          familyPass ? "PASS" : "FAIL",parent,maximumAnalogZ,
          maximumBiasedZ);
    }
    allPass=familyPass && allPass;
  }

  // Require the staged validation policy to cover a multi-channel family, a
  // rare branch whose analytic variance is reduced, a neutral/sink family,
  // and a packet-split biased reaction.  Passing with only identity rows
  // would not close the Stage-7 implementation gate.
  allPass=familyCount>0 && biasedFamilyCount>0 && reactionCount>0 &&
      sinkCount>0 && coveredMany && coveredRare && coveredReducedRare &&
      coveredPacket && allPass;
  int localPass=allPass ? 1 : 0,globalPass=0;
  MPI_Allreduce(&localPass,&globalPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);

  // Exercise the actual Stage-7 serializer in every accepted decomposition.
  // A high, reserved output index makes the focused artifact unambiguous and
  // avoids collision with ordinary science output indices.  The integrated
  // runner launches each decomposition in its own immutable case directory.
  const long int diagnosticOutputIndex=7000000+PIC::nTotalThreads;
  const std::vector<double> noPreviousProfile(PIC::nTotalSpecies,0.0);
  Comet::Diagnostics::WriteSamplingReport(eventsPerFamily,
      diagnosticOutputIndex,eventsPerFamily,noPreviousProfile,false);
  if (PIC::ThisThread==0) {
    std::printf("%s channel-importance diagnostics output=%ld "
        "csv=PT/diagnostics/comet_stage7_importance.out=%ld.csv "
        "json=PT/diagnostics/comet_stage7_importance.out=%ld.json\n",
        globalPass ? "PASS" : "FAIL",diagnosticOutputIndex,
        diagnosticOutputIndex,diagnosticOutputIndex);
    std::printf("Channel-importance validation result: %s families=%d "
        "biased-families=%d reactions=%d sinks=%d events=%ld\n",
        globalPass ? "PASS" : "FAIL",familyCount,biasedFamilyCount,
        reactionCount,sinkCount,eventsPerFamily);
  }
  return globalPass!=0;
#endif
}

bool Comet::Chemistry::ValidateTotalReactionForcing(
    long int eventsPerFixture,bool validateTransportedCollision) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  if (PIC::ThisThread==0)
    std::fprintf(stderr,
        "total-reaction-forcing validation requires individual corrections\n");
  return false;
#else
  namespace Forcing=Comet::TotalReactionForcingPlanner;
  if (!Enabled || !TotalReactionForcingEnabled ||
      eventsPerFixture<1000 || TotalReactionForcingStateOffset<0) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "total-reaction-forcing validation requires an enabled Stage-8 build and at least 1000 events\n");
    return false;
  }
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  if (PIC::DomainBlockDecomposition::nLocalBlocks<=0 ||
      PIC::DomainBlockDecomposition::BlockTable[0]==NULL ||
      PIC::DomainBlockDecomposition::BlockTable[0]->block==NULL) {
    std::fprintf(stderr,"rank %d has no block for Stage-8 validation\n",
        PIC::ThisThread);
    std::abort();
  }
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
      PIC::DomainBlockDecomposition::BlockTable[0];

  auto validationUniform=[](std::uint64_t event,std::uint32_t stream) {
    const std::uint32_t key=static_cast<std::uint32_t>(
        (event*UINT64_C(0x9e3779b9)+stream*UINT64_C(0x85ebca6b))&
        TotalReactionSeedMask);
    return TotalReactionUnitOpen(key,event+stream);
  };
  auto conditionalCdf=[](double opticalPosition,double totalOpticalDepth) {
    return -std::expm1(-opticalPosition)/
        -std::expm1(-totalOpticalDepth);
  };
  auto ksDistance=[](std::vector<double> values) {
    std::sort(values.begin(),values.end());
    double distance=0.0;
    const double count=static_cast<double>(values.size());
    for (std::size_t i=0;i<values.size();i++) {
      const double lower=static_cast<double>(i)/count;
      const double upper=static_cast<double>(i+1)/count;
      distance=std::max(distance,std::fabs(values[i]-lower));
      distance=std::max(distance,std::fabs(upper-values[i]));
    }
    return distance;
  };

  // Homogeneous split and multigeneration decay are deterministic carrier
  // identities.  Repeating the production split exposes an accidental linear
  // probability approximation that a one-step comparison could miss.
  const double homogeneousHazard=0.37;
  const double homogeneousDuration=0.4;
  const int homogeneousSteps=3;
  double homogeneousSurvival=1.0;
  for (int step=0;step<homogeneousSteps;step++)
    homogeneousSurvival=Forcing::SplitCarrier(
        Forcing::BuildSplitPlan(homogeneousHazard,homogeneousDuration),
        homogeneousSurvival).Survivor;
  const double homogeneousTau=homogeneousHazard*homogeneousDuration*
      homogeneousSteps;
  const double homogeneousAnalyticSurvival=std::exp(-homogeneousTau);
  const double homogeneousAnalyticReacted=-std::expm1(-homogeneousTau);
  const double homogeneousRepresentedReacted=1.0-homogeneousSurvival;
  const double homogeneousClosure=std::fabs(homogeneousSurvival+
      homogeneousRepresentedReacted-1.0);
  const bool homogeneousPass=std::fabs(homogeneousSurvival-
      homogeneousAnalyticSurvival)<=3.0e-12*homogeneousAnalyticSurvival &&
      homogeneousClosure<=2.0e-12;

  struct cPositionFixtureResult {
    long int FirstCount=0;
    double FirstExpected=0.0;
    double Mean=0.0;
    double ExpectedMean=0.0;
    double CdfDistance=0.0;
    long int EndpointViolations=0;
    bool Pass=false;
  };
  cPositionFixtureResult constantResult,piecewiseResult;
  const double constantTau=0.75;
  const Forcing::cSplitPlan constantPlan=
      Forcing::BuildSplitPlan(constantTau,1.0);
  std::vector<double> constantCdf;
  constantCdf.reserve(eventsPerFixture);
  for (long int event=0;event<eventsPerFixture;event++) {
    const double position=Forcing::SampleConditionalReactionTime(
        constantPlan,validationUniform(event,11));
    if (position<0.0 || position>1.0) constantResult.EndpointViolations++;
    constantResult.FirstCount++;
    constantResult.Mean+=position;
    constantCdf.push_back(conditionalCdf(constantTau*position,constantTau));
  }
  constantResult.Mean/=eventsPerFixture;
  constantResult.FirstExpected=1.0;
  constantResult.ExpectedMean=1.0/constantTau-
      1.0/std::expm1(constantTau);
  constantResult.CdfDistance=ksDistance(constantCdf);
  constantResult.Pass=constantResult.EndpointViolations==0 &&
      std::fabs(constantResult.Mean-constantResult.ExpectedMean)<=
          8.0/std::sqrt(static_cast<double>(eventsPerFixture)) &&
      constantResult.CdfDistance<=
          8.0/std::sqrt(static_cast<double>(eventsPerFixture));

  const double firstFraction=0.35,piecewiseTau1=0.25,piecewiseTau2=0.90;
  const Forcing::cPiecewisePlan piecewisePlan=Forcing::BuildPiecewisePlan({
      {piecewiseTau1/firstFraction,firstFraction},
      {piecewiseTau2/(1.0-firstFraction),1.0-firstFraction}});
  std::vector<double> piecewiseCdf;
  piecewiseCdf.reserve(eventsPerFixture);
  for (long int event=0;event<eventsPerFixture;event++) {
    const Forcing::cPiecewiseReactionPosition position=
        Forcing::SampleConditionalReactionPosition(piecewisePlan,
            validationUniform(event,29));
    if (position.ElapsedTime<0.0 || position.ElapsedTime>1.0)
      piecewiseResult.EndpointViolations++;
    if (position.Segment==0) piecewiseResult.FirstCount++;
    piecewiseResult.Mean+=position.ElapsedTime;
    piecewiseCdf.push_back(conditionalCdf(position.OpticalDepth,
        piecewiseTau1+piecewiseTau2));
  }
  piecewiseResult.Mean/=eventsPerFixture;
  piecewiseResult.FirstExpected=-std::expm1(-piecewiseTau1)/
      -std::expm1(-(piecewiseTau1+piecewiseTau2));
  {
    const double rate1=piecewiseTau1/firstFraction;
    const double rate2=piecewiseTau2/(1.0-firstFraction);
    const double e1=std::exp(-piecewiseTau1);
    const double e2=std::exp(-piecewiseTau2);
    const double moment1=1.0/rate1-
        (firstFraction+1.0/rate1)*e1;
    const double secondDuration=1.0-firstFraction;
    const double moment2=e1*(firstFraction*(1.0-e2)+1.0/rate2-
        (secondDuration+1.0/rate2)*e2);
    piecewiseResult.ExpectedMean=(moment1+moment2)/
        -std::expm1(-(piecewiseTau1+piecewiseTau2));
  }
  piecewiseResult.CdfDistance=ksDistance(piecewiseCdf);
  const double regionSigma=std::sqrt(eventsPerFixture*
      piecewiseResult.FirstExpected*(1.0-piecewiseResult.FirstExpected));
  piecewiseResult.Pass=piecewiseResult.EndpointViolations==0 &&
      std::fabs(piecewiseResult.FirstCount-
          eventsPerFixture*piecewiseResult.FirstExpected)<=8.0*regionSigma &&
      std::fabs(piecewiseResult.Mean-piecewiseResult.ExpectedMean)<=
          8.0/std::sqrt(static_cast<double>(eventsPerFixture)) &&
      piecewiseResult.CdfDistance<=
          8.0/std::sqrt(static_cast<double>(eventsPerFixture));

  const int decaySteps=7;
  const double decayTotalTau=0.91;
  double decaySurvival=1.0;
  for (int step=0;step<decaySteps;step++)
    decaySurvival=Forcing::SplitCarrier(
        Forcing::BuildSplitPlan(decayTotalTau/decaySteps,1.0),
        decaySurvival).Survivor;
  const double decayAnalyticSurvival=std::exp(-decayTotalTau);
  const double decayAnalyticReacted=-std::expm1(-decayTotalTau);
  const double decayRepresentedReacted=1.0-decaySurvival;
  const double decayClosure=std::fabs(decaySurvival+
      decayRepresentedReacted-1.0);
  const bool decayPass=std::fabs(decaySurvival-decayAnalyticSurvival)<=
      3.0e-12*decayAnalyticSurvival && decayClosure<=3.0e-12;

  // The configured C2 -> 2C fixture composes every accepted estimator layer:
  // Stage-8 occurrence P, Stage-7 likelihood L, and Stage-6 N=8 packet split.
  const ChemicalSolver::cReactionNetwork &network=GetNetwork();
  std::size_t compositionReaction=network.GetReactionCount();
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++)
    if (std::strcmp(network.GetReaction(reaction).Id,"C2_2C")==0) {
      compositionReaction=reaction;
      break;
    }
  bool compositionPass=compositionReaction<network.GetReactionCount();
  int compositionParent=-1;
  long int compositionPackets=0;
  std::vector<int> compositionProducts;
  double occurrenceFactor=0.0,compositionLikelihood=0.0;
  double expectedCorrection=0.0,packetCorrectionSum=0.0;
  double compositionMomentumError=INFINITY,compositionEnergyError=INFINITY;
  if (compositionPass) {
    compositionParent=network.GetReactantTerm(compositionReaction,0).Species;
    compositionPackets=GetDaughterPacketCount(compositionReaction);
    ExpandConcreteProducts(compositionReaction,compositionProducts);
    occurrenceFactor=Forcing::BuildSplitPlan(0.25,0.7).ReactedFraction;
    compositionLikelihood=GetChannelLikelihoodRatio(compositionReaction);
    expectedCorrection=occurrenceFactor*compositionLikelihood;
    std::vector<double> unsplit(compositionProducts.size(),
        occurrenceFactor/static_cast<double>(compositionProducts.size()));
    const Comet::DaughterPacketPlanner::cPlan packetPlan=
        Comet::DaughterPacketPlanner::BuildPlan(unsplit,
            static_cast<std::size_t>(compositionPackets),
            static_cast<std::size_t>(std::max<long int>(1,
                MaximumDaughterPacketCount)),
            static_cast<std::size_t>(Comet::RadialSimulation::
                MaximumPhotochemicalProductsPerReaction));
    for (long int packet=0;packet<compositionPackets;packet++)
      packetCorrectionSum+=packetPlan.Correction(packet,0)*
          compositionProducts.size()*compositionLikelihood;

    const double parentVelocity[3]={400.0,-120.0,55.0};
    std::vector<double> velocities;
    GenerateProductVelocitiesDeterministic(compositionProducts,
        parentVelocity,network.GetReaction(compositionReaction).ExcessEnergyJ,
        UINT32_C(731245),velocities);
    double productMass=0.0,momentum[3]={0.0,0.0,0.0},energy=0.0;
    for (std::size_t product=0;product<compositionProducts.size();product++) {
      const double mass=PIC::MolecularData::GetMass(
          compositionProducts[product]);
      productMass+=mass;
      const double *v=velocities.data()+3*product;
      for (int dimension=0;dimension<3;dimension++)
        momentum[dimension]+=mass*v[dimension];
      energy+=0.5*mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }
    double momentumDifference2=0.0,parentSpeed2=0.0;
    for (int dimension=0;dimension<3;dimension++) {
      const double difference=momentum[dimension]-
          productMass*parentVelocity[dimension];
      momentumDifference2+=difference*difference;
      parentSpeed2+=parentVelocity[dimension]*parentVelocity[dimension];
    }
    compositionMomentumError=std::sqrt(momentumDifference2)/
        std::max(productMass*std::sqrt(parentSpeed2),
            std::numeric_limits<double>::min());
    const double expectedEnergy=0.5*productMass*parentSpeed2+
        network.GetReaction(compositionReaction).ExcessEnergyJ;
    compositionEnergyError=std::fabs(energy-expectedEnergy)/
        std::max(expectedEnergy,std::numeric_limits<double>::min());
    compositionPass=compositionPackets==8 &&
        compositionProducts.size()==2 &&
        std::fabs(packetCorrectionSum-expectedCorrection)<=
            2.0e-12*expectedCorrection &&
        compositionMomentumError<=1.0e-11 &&
        compositionEnergyError<=1.0e-11;
  }

  long long localParents=0,localSurvivors=0,localReacted=0;
  long long localEmitted=0,localProcessed=0,localRemainingMoves=0;
  long long localCapacityRollbacks=0,localByteRollbacks=0;
  long long localRngRollbacks=0,localPathRollbacks=0;
  long long localDoubleChemistry=0;
  double localCorrectionError=0.0,localPositionError=0.0;
  bool transactionPass=compositionPass &&
      IsTotalReactionForcingConfigured(compositionParent);
  if (transactionPass) {
    const double savedReserve=
        Comet::RadialSimulation::AllocationReserveFraction;
    const bool savedDiagnostics=Comet::Diagnostics::Enabled;
    Comet::Diagnostics::Enabled=false;
    double x[3]={0.5*(node->xmin[0]+node->xmax[0]),0.0,0.0};
    double v[3]={400.0,-120.0,55.0};
    const double parentDt=node->block->GetLocalTimeStep(compositionParent);
    auto initializeParent=[&](long int particle) {
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(particle);
      PIC::ParticleBuffer::SetI(compositionParent,data);
      PIC::ParticleBuffer::SetX(x,data);
      PIC::ParticleBuffer::SetV(v,data);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,data);
      PIC::ParticleBuffer::SetPrev(-1,data);
      PIC::ParticleBuffer::SetNext(-1,data);
      ClearTotalReactionState(data);
      ChemistryProcessedFlag(data)=0;
    };

    // Capacity failure is forced before the first Stage-8 random draw. A
    // reseed/replay comparison proves that the failed transaction consumed no
    // RNG, while a whole-record memcmp covers correction, state, path, and all
    // application bytes.
    long int failed=PIC::ParticleBuffer::GetNewParticle(true);
    initializeParent(failed);
    long int failedList=-1;
    std::vector<unsigned char> beforeFailure(
        PIC::ParticleBuffer::ParticleDataLength);
    std::memcpy(beforeFailure.data(),
        PIC::ParticleBuffer::GetParticleDataPointer(failed),
        beforeFailure.size());
    const int rollbackSeed=913001+PIC::ThisThread;
    rnd_seed(rollbackSeed);
    Comet::RadialSimulation::AllocationReserveFraction=1.0;
    const bool incorrectlyAccepted=ProcessTotalReactionForcedParticle(failed,
        failedList,node,compositionParent,
        ParentReactionRates[compositionParent],parentDt);
    const double afterFailureDraw=rnd();
    rnd_seed(rollbackSeed);
    const double referenceDraw=rnd();
    const bool byteRollback=std::memcmp(beforeFailure.data(),
        PIC::ParticleBuffer::GetParticleDataPointer(failed),
        beforeFailure.size())==0;
    const bool pathRollback=failedList==failed &&
        PIC::ParticleBuffer::GetNext(failed)==-1 &&
        PIC::ParticleBuffer::GetPrev(failed)==-1;
    localCapacityRollbacks=incorrectlyAccepted ? 0 : 1;
    localByteRollbacks=byteRollback ? 1 : 0;
    localRngRollbacks=afterFailureDraw==referenceDraw ? 1 : 0;
    localPathRollbacks=pathRollback ? 1 : 0;
    while (failedList!=-1)
      PIC::ParticleBuffer::DeleteParticle(failedList,failedList);

    Comet::RadialSimulation::AllocationReserveFraction=savedReserve;
    LocalChemistryAllocationPressure=false;
    rnd_seed(914003+PIC::ThisThread);
    long int parent=PIC::ParticleBuffer::GetNewParticle(true);
    initializeParent(parent);
    long int transactionList=-1;
    const bool committed=ProcessTotalReactionForcedParticle(parent,
        transactionList,node,compositionParent,
        ParentReactionRates[compositionParent],parentDt);
    localParents=committed ? 1 : 0;
    // One incoming parent creates one logical reacted branch. Its K*N pending
    // markers are numerical product carriers, not additional reactions.
    localReacted=committed ? 1 : 0;
    for (long int particle=transactionList;particle!=-1;
         particle=PIC::ParticleBuffer::GetNext(particle)) {
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(particle);
      const cTotalReactionForcingState state=GetTotalReactionState(data);
      if (state.Code==0) {
        localSurvivors++;
        if (ChemistryProcessedFlag(data)==1) localProcessed++;
        continue;
      }
      const double oldCorrection=
          PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
      const double parentCarrier=node->block->GetLocalParticleWeight(
          compositionParent)/parentDt;
      double eventTime=-1.0;
      const int query=QueryTotalReactionEvent(particle,compositionParent,node,
          parentDt,parentDt,eventTime);
      if (query!=PIC::Mover::TrajectoryEventPending) transactionPass=false;
      AdvanceTotalReactionEvent(particle,compositionParent,node,parentDt,
          eventTime);
      localPositionError=std::max(localPositionError,
          std::fabs(GetTotalReactionState(data).RemainingFraction));
      double eventPosition[3],eventVelocity[3];
      PIC::ParticleBuffer::GetX(eventPosition,data);
      PIC::ParticleBuffer::GetV(eventVelocity,data);
      int eventSpecies=compositionParent;
      double remainingTime=parentDt-eventTime;
      const int result=ProcessTotalReactionEvent(particle,eventPosition,
          eventVelocity,eventSpecies,node,parentDt,remainingTime);
      if (result!=PIC::Mover::TrajectoryEventContinue ||
          eventSpecies<0 || eventSpecies>=PIC::nTotalSpecies)
        transactionPass=false;
      const double productCarrier=node->block->GetLocalParticleWeight(
          eventSpecies)/node->block->GetLocalTimeStep(eventSpecies);
      const double represented=productCarrier*
          PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
      const double expected=parentCarrier*oldCorrection*
          compositionProducts.size()*compositionLikelihood;
      localCorrectionError=std::max(localCorrectionError,
          std::fabs(represented-expected)/
              std::max(std::fabs(expected),
                  std::numeric_limits<double>::min()));
      if (remainingTime>0.0) localRemainingMoves++;
      if (GetTotalReactionState(data).Code!=0) localDoubleChemistry++;
      if (ChemistryProcessedFlag(data)==1) localProcessed++;
      localEmitted++;
    }
    transactionPass=committed && localSurvivors==1 && localReacted==1 &&
        localEmitted==16 && localProcessed==17 &&
        localRemainingMoves==16 && localCapacityRollbacks==1 &&
        localByteRollbacks==1 && localRngRollbacks==1 &&
        localPathRollbacks==1 && localDoubleChemistry==0 &&
        localCorrectionError<=2.0e-12 && localPositionError<=2.0e-12 &&
        transactionPass;
    while (transactionList!=-1)
      PIC::ParticleBuffer::DeleteParticle(transactionList,transactionList);
    Comet::Diagnostics::Enabled=savedDiagnostics;
    rnd_seed(static_cast<int>(Comet::Model::RandomSeed)+PIC::ThisThread);
  }

  // The transaction fixture above intentionally calls the three application
  // callbacks directly so it can isolate byte/RNG rollback.  It does not,
  // however, exercise the core's event-versus-geometry arbitration.  This
  // second, deliberately small fixture therefore sends one real C2 -> 2C
  // transaction through ExecutionTrackDefault(), which is the same first-
  // order mover, face-resampling, MPI exchange, and continuation loop used by
  // PIC::TimeStep().  Only one global source is used, so all physical counts
  // below are invariant at 1, 4, and 12 ranks.
  //
  // Product markers are scheduled strictly after an internal block face.  A
  // separate C2_ION sink marker is scheduled at the exact face time.  The
  // latter is the sharpest ownership check: geometry must win the tie,
  // Advance must commit RemainingFraction==0 before the particle is copied,
  // and the destination continuation must invoke Process at zero flight time;
  // only the mover may then release the sink's particle-buffer slot.
  bool moverPass=compositionPass;
  bool moverProductionPass=false;
  int moverMpiFace=0;
  long long moverGlobalCounts[12]={};
  long long moverGlobalFaceTransfers=0;
  long int moverGlobalInitialPending=0,moverGlobalContinuationMoves=0;
  int moverGlobalSweeps=0;
  long int moverGlobalFinalPending=0;
  double moverCorrectionError=INFINITY,moverEventPositionError=INFINITY;
  double moverRemainingFractionError=INFINITY,moverEndpointError=INFINITY;

  // Stage 9 extends the Stage-8 mover fixture only when explicitly selected.
  // Keep every result outside the nested mover scope so rank zero can emit a
  // single strict record after all reductions have completed.
  bool transportedCollisionPass=!validateTransportedCollision;
  int transportedCollisionPartner=-1,transportedCollisionOwner=-1;
  long long transportedDaughters=0,transportedPartners=0;
  long long transportedDaughtersChanged=0,transportedPartnersChanged=0;
  long long transportedUnchangedRecords=0,transportedUnfinished=0;
  long long transportedWrongOwner=0,transportedCollisionCycles=0;
  unsigned long long transportedTrials=0,transportedAccepted=0;
  unsigned long long transportedUpdates=0,transportedProbabilityChecks=0;
  unsigned long long transportedInvalidProbabilities=0;
  unsigned long long transportedRoundoffClamps=0;
  double transportedProbabilityMinimum=INFINITY;
  double transportedProbabilityMaximum=-INFINITY;
  double transportedCorrectionError=INFINITY;
  double transportedRateError=INFINITY;
  double transportedMomentumError=INFINITY;
  if (moverPass) {
    std::size_t sinkReaction=network.GetReactionCount();
    for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++)
      if (std::strcmp(network.GetReaction(reaction).Id,"C2_ION")==0) {
        sinkReaction=reaction;
        break;
      }
    moverPass=sinkReaction<network.GetReactionCount() &&
        network.GetReaction(sinkReaction).ProductCount==0 &&
        network.GetReactantTerm(sinkReaction,0).Species==compositionParent &&
        compositionProducts.size()==2 &&
        compositionProducts[0]==compositionProducts[1];

    using cNode=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
    std::vector<cNode*> leaves;
    for (cNode *leaf=PIC::Mesh::mesh->BranchBottomNodeList;leaf!=NULL;
         leaf=leaf->nextBranchBottomNode)
      if (leaf->IsUsedInCalculationFlag) leaves.push_back(leaf);
    std::sort(leaves.begin(),leaves.end(),[](cNode *left,cNode *right) {
      return left->xmin[0]<right->xmin[0];
    });

    cNode *source=NULL,*destination=NULL;
    for (std::size_t index=1;index<leaves.size();index++) {
      cNode *left=leaves[index-1],*right=leaves[index];
      const double leftWidth=left->xmax[0]-left->xmin[0];
      const double rightWidth=right->xmax[0]-right->xmin[0];
      const double scale=std::max({std::fabs(left->xmax[0]),leftWidth,
          rightWidth,1.0});
      if (std::fabs(left->xmax[0]-right->xmin[0])>1.0e-12*scale)
        continue;
      const bool mpiFace=left->Thread!=right->Thread;
      if (mpiFace!=(PIC::nTotalThreads>1)) continue;
      source=left;
      destination=right;
      break;
    }
    moverPass=source!=NULL && destination!=NULL && moverPass;

    if (moverPass) {
      moverMpiFace=source->Thread!=destination->Thread ? 1 : 0;
      const int productSpecies=compositionProducts[0];
      const double sourceWidth=source->xmax[0]-source->xmin[0];
      const double destinationWidth=
          destination->xmax[0]-destination->xmin[0];
      const double minimumWidth=std::min(sourceWidth,destinationWidth);
      const double face=source->xmax[0];

      // Unit W/dt on both sides removes stochastic face multiplicity from this
      // event-ordering fixture.  The Stage-4 mover suite independently covers
      // clone/roulette statistics over unequal local weights and timesteps.
      // Update every allocated owner/boundary copy so the source rank and the
      // destination rank observe identical face metadata.
      for (cNode *leaf:leaves) if (leaf->block!=NULL) {
        leaf->block->SetLocalTimeStep(1.0,compositionParent);
        leaf->block->SetLocalParticleWeight(1.0,compositionParent);
      }

      PIC::ParticleBuffer::DeleteAllParticles();
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      const double faceDistance=0.10*minimumWidth;
      const double parentVelocity=0.40*minimumWidth;
      double initialPosition[3]={face-faceDistance,0.0,0.0};
      double initialVelocity[3]={parentVelocity,0.0,0.0};

      // Reproduce the mover's stable q-form intersection calculation.  The
      // resulting double, with parent dt exactly one, is stored directly in
      // the sink marker.  Query therefore returns bitwise the same time as the
      // geometric face solver and exercises the documented equality rule
      // without a tolerance-dependent near-tie.
      auto faceIntersectionTime=[&]() {
        const double a=parentVelocity*parentVelocity;
        const double b=2.0*initialPosition[0]*parentVelocity;
        const double c=initialPosition[0]*initialPosition[0]-face*face;
        const double discriminant=b*b-4.0*a*c;
        const double root=std::sqrt(discriminant);
        const double q=-0.5*(b+std::copysign(root,b));
        double first=1.0;
        const double candidates[2]={q/a,c/q};
        for (double candidate:candidates)
          if (candidate>0.0 && candidate<first) first=candidate;
        return first;
      };
      const double faceTime=faceIntersectionTime();
      const double afterFaceEventDelay=0.25;
      const double productEventFraction=faceTime+afterFaceEventDelay;
      moverPass=faceTime>0.0 && faceTime<0.5 &&
          productEventFraction<1.0 && std::isfinite(faceTime) && moverPass;

      long int transportList=-1;
      long long localPendingProducts=0;
      double localMaximumProductSpeed=0.0;
      const bool savedDiagnostics=Comet::Diagnostics::Enabled;
      Comet::Diagnostics::Enabled=false;
      if (source->Thread==PIC::ThisThread && moverPass) {
        auto initializeMoverParent=[&](long int particle) {
          PIC::ParticleBuffer::byte *data=
              PIC::ParticleBuffer::GetParticleDataPointer(particle);
          PIC::ParticleBuffer::SetI(compositionParent,data);
          PIC::ParticleBuffer::SetX(initialPosition,data);
          PIC::ParticleBuffer::SetV(initialVelocity,data);
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,data);
          PIC::ParticleBuffer::SetPrev(-1,data);
          PIC::ParticleBuffer::SetNext(-1,data);
          PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
          ClearTotalReactionState(data);
          ChemistryProcessedFlag(data)=0;
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
          // Face diagnostics read the internal-energy fields.  A production
          // parent already owns initialized values; make the focused raw
          // buffer allocation obey the same invariant before it is cloned.
          PIC::IDF::InitRotTemp(180.0,data);
          PIC::IDF::InitVibTemp(180.0,data);
#endif
        };

        long int parent=PIC::ParticleBuffer::GetNewParticle(true);
        initializeMoverParent(parent);
        // Internal-energy initialization may consume random values.  Reseed
        // afterward so this uses the same known productive C2 channel as the
        // accepted transaction fixture on every possible source owner.
        rnd_seed(914003+PIC::ThisThread);
        const bool committed=ProcessTotalReactionForcedParticle(parent,
            transportList,source,compositionParent,
            ParentReactionRates[compositionParent],1.0);
        moverPass=committed;
        for (long int particle=transportList;particle!=-1;
             particle=PIC::ParticleBuffer::GetNext(particle)) {
          PIC::ParticleBuffer::byte *data=
              PIC::ParticleBuffer::GetParticleDataPointer(particle);
          PIC::ParticleBuffer::SetX(initialPosition,data);
          PIC::ParticleBuffer::SetV(initialVelocity,data);
          PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
          cTotalReactionForcingState state=GetTotalReactionState(data);
          if (state.Code==0) continue;
          if (DecodeTotalReactionIndex(state.Code)!=compositionReaction)
            moverPass=false;
          state.RemainingFraction=productEventFraction;
          SetTotalReactionState(data,state);
          localPendingProducts++;

          std::vector<double> velocities;
          GenerateProductVelocitiesDeterministic(compositionProducts,
              initialVelocity,
              network.GetReaction(compositionReaction).ExcessEnergyJ,
              state.KeyAndFlags&TotalReactionSeedMask,velocities);
          for (std::size_t product=0;product<compositionProducts.size();product++) {
            const double *velocity=velocities.data()+3*product;
            localMaximumProductSpeed=std::max(localMaximumProductSpeed,
                std::sqrt(velocity[0]*velocity[0]+velocity[1]*velocity[1]+
                    velocity[2]*velocity[2]));
          }
        }
      }

      long long globalPendingProducts=0;
      double globalMaximumProductSpeed=0.0;
      MPI_Allreduce(&localPendingProducts,&globalPendingProducts,1,
          MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&localMaximumProductSpeed,&globalMaximumProductSpeed,1,
          MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
      moverPass=globalPendingProducts==16 &&
          globalMaximumProductSpeed>0.0 &&
          std::isfinite(globalMaximumProductSpeed) && moverPass;

      // Choose the product dt from the actually generated recoil velocities.
      // Every product endpoint then remains at least forty continuation
      // displacements from the nearest destination face, so an endpoint
      // mismatch cannot be hidden by a second face transfer.
      // The accepted Stage-8 fixture retains its original 0.05-block
      // continuation.  Stage 9 shortens only its private copy to 0.01 block:
      // the conditional event is 0.10 block inside the destination and a
      // radial cell is 0.125 block wide, so every recoil endpoint then stays
      // in the same physical cell and can enter one real NTC cell call.
      const double productTimeStep=
          (validateTransportedCollision ? 0.01 : 0.05)*minimumWidth/
          globalMaximumProductSpeed;
      moverPass=productTimeStep>0.0 && std::isfinite(productTimeStep) &&
          moverPass;
      for (cNode *leaf:leaves) if (leaf->block!=NULL) {
        leaf->block->SetLocalTimeStep(productTimeStep,productSpecies);
        // W_product/dt_product == W_parent/dt_parent == 1.  Consequently the
        // post-event correction sum itself is the represented product rate.
        leaf->block->SetLocalParticleWeight(productTimeStep,productSpecies);
      }

      if (source->Thread==PIC::ThisThread && moverPass) {
        // Add one explicit product-free marker to the already-built list.  It
        // is not a second chemistry-selection test; C2_ION selection and
        // likelihood are covered by Stage 7.  Its sole purpose is to require
        // the mover-owned delete branch after an exact face/event tie.
        const long int sink=PIC::ParticleBuffer::GetNewParticle(true);
        PIC::ParticleBuffer::byte *sinkData=
            PIC::ParticleBuffer::GetParticleDataPointer(sink);
        PIC::ParticleBuffer::SetI(compositionParent,sinkData);
        PIC::ParticleBuffer::SetX(initialPosition,sinkData);
        PIC::ParticleBuffer::SetV(initialVelocity,sinkData);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(0.5,sinkData);
        PIC::ParticleBuffer::SetPrev(-1,sinkData);
        PIC::ParticleBuffer::SetNext(-1,sinkData);
        PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,sinkData);
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        PIC::IDF::InitRotTemp(180.0,sinkData);
        PIC::IDF::InitVibTemp(180.0,sinkData);
#endif
        cTotalReactionForcingState sinkState;
        sinkState.Code=EncodeTotalReactionSink(sinkReaction);
        sinkState.KeyAndFlags=TotalReactionLeaderFlag;
        sinkState.RemainingFraction=faceTime;
        SetTotalReactionState(sinkData,sinkState);
        ChemistryProcessedFlag(sinkData)=1;
        KeepParticle(sink,transportList);

        int i=0,j=0,k=0;
        const long int centerNode=PIC::Mesh::mesh->FindCellIndex(
            initialPosition,i,j,k,source,false);
        // FindCellIndex returns the center-node storage index, which includes
        // ghost-cell offsets.  Particle lists use the compact physical-cell
        // index instead.  Using centerNode here can address the adjacent
        // temporary-list array for a particle in the outermost cell, making
        // the first MoveParticles() call appear to skip the fixture entirely.
        const long int cell=i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
        if (centerNode<0 || cell<0 ||
            cell>=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_ ||
            source->block==NULL ||
            source->block->FirstCellParticleTable[cell]!=-1) {
          moverPass=false;
        }
        else source->block->FirstCellParticleTable[cell]=transportList;
      }

      int localSetupPass=moverPass ? 1 : 0,globalSetupPass=0;
      MPI_Allreduce(&localSetupPass,&globalSetupPass,1,MPI_INT,MPI_MIN,
          MPI_GLOBAL_COMMUNICATOR);
      moverPass=globalSetupPass!=0;
      if (moverPass) {
        TotalReactionMoverValidationAudit=
            cTotalReactionMoverValidationAudit();
        TotalReactionMoverValidationAudit.Active=true;
        TotalReactionMoverValidationAudit.ExpectedOwner=destination->Thread;
        TotalReactionMoverValidationAudit.ExpectedSinkEventRadius=
            face+PIC::Mesh::mesh->EPS;
        TotalReactionMoverValidationAudit.ExpectedProductEventRadius=
            face+PIC::Mesh::mesh->EPS+
            parentVelocity*afterFaceEventDelay;
        TotalReactionMoverValidationAudit.ExpectedSinkRemainingFraction=
            1.0-faceTime;
        TotalReactionMoverValidationAudit.ExpectedProductRemainingFraction=
            1.0-productEventFraction;

        PIC::Mover::ResetLocalTimeStepFluxDiagnostic(true);
        double movingTime=0.0,exchangeTime=0.0;
        PIC::TimeStepInternal::ExecutionTrackDefault(movingTime,exchangeTime);
        TotalReactionMoverValidationAudit.Active=false;

        const PIC::Mover::cLocalTimeStepFluxDiagnostic &diagnostic=
            PIC::Mover::LocalTimeStepFluxDiagnostic;
        MPI_Allreduce(&diagnostic.FaceTransferCount,
            &moverGlobalFaceTransfers,1,MPI_LONG_LONG,MPI_SUM,
            MPI_GLOBAL_COMMUNICATOR);
        MPI_Allreduce(&diagnostic.InitialPendingParticleCount,
            &moverGlobalInitialPending,1,MPI_LONG,MPI_MAX,
            MPI_GLOBAL_COMMUNICATOR);
        MPI_Allreduce(&diagnostic.ContinuationMovedParticleCount,
            &moverGlobalContinuationMoves,1,MPI_LONG,MPI_SUM,
            MPI_GLOBAL_COMMUNICATOR);
        MPI_Allreduce(&diagnostic.ContinuationSweepCount,
            &moverGlobalSweeps,1,MPI_INT,MPI_MAX,
            MPI_GLOBAL_COMMUNICATOR);
        MPI_Allreduce(&diagnostic.FinalPendingParticleCount,
            &moverGlobalFinalPending,1,MPI_LONG,MPI_MAX,
            MPI_GLOBAL_COMMUNICATOR);

        long long localCounts[12]={
            TotalReactionMoverValidationAudit.ProductProcesses,
            TotalReactionMoverValidationAudit.SinkProcesses,
            TotalReactionMoverValidationAudit.WrongOwnerProcesses,0,0,0,0,0,
            0,0,0,TotalReactionMoverValidationAudit.ProductPositiveRemaining};
        double localProductRate=0.0,localSurvivorRate=0.0;
        double localEndpointError=0.0;
        const int cellCount=
            _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
        for (unsigned int blockIndex=0;
             blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;
             blockIndex++) {
          cNode *scanNode=
              PIC::DomainBlockDecomposition::BlockTable[blockIndex];
          if (scanNode==NULL || scanNode->block==NULL) continue;
          for (int cell=0;cell<cellCount;cell++)
            for (long int particle=
                     scanNode->block->FirstCellParticleTable[cell];
                 particle!=-1;particle=PIC::ParticleBuffer::GetNext(particle)) {
              localCounts[3]++;
              PIC::ParticleBuffer::byte *data=
                  PIC::ParticleBuffer::GetParticleDataPointer(particle);
              if (GetTotalReactionState(data).Code!=0) localCounts[6]++;
              if (PIC::Mover::GetLocalTimeStepRemainderFraction(data)!=0.0)
                localCounts[7]++;
              if (ChemistryProcessedFlag(data)==1) localCounts[8]++;
              const int species=PIC::ParticleBuffer::GetI(data);
              const double correction=
                  PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
              if (species==compositionParent) {
                localCounts[4]++;
                localSurvivorRate+=scanNode->block->GetLocalParticleWeight(
                    species)/scanNode->block->GetLocalTimeStep(species)*
                    correction;
              }
              else if (species==productSpecies) {
                localCounts[5]++;
                localProductRate+=scanNode->block->GetLocalParticleWeight(
                    species)/scanNode->block->GetLocalTimeStep(species)*
                    correction;
                double position[3];
                PIC::ParticleBuffer::GetX(position,data);
                for (const cTotalReactionMoverValidationEndpoint &endpoint:
                     TotalReactionMoverValidationAudit.ProductEndpoints)
                  if (endpoint.Particle==particle) {
                    localCounts[9]++;
                    const double radius=std::sqrt(position[0]*position[0]+
                        position[1]*position[1]+position[2]*position[2]);
                    localEndpointError=std::max(localEndpointError,
                        std::fabs(radius-endpoint.Radius));
                    break;
                  }
              }
              else localCounts[10]++;
            }
        }
        MPI_Allreduce(localCounts,moverGlobalCounts,12,MPI_LONG_LONG,MPI_SUM,
            MPI_GLOBAL_COMMUNICATOR);
        double globalRates[2]={0.0,0.0};
        const double localRates[2]={localProductRate,localSurvivorRate};
        MPI_Allreduce(localRates,globalRates,2,MPI_DOUBLE,MPI_SUM,
            MPI_GLOBAL_COMMUNICATOR);
        double localErrors[3]={
            TotalReactionMoverValidationAudit.MaximumEventPositionError,
            TotalReactionMoverValidationAudit.MaximumRemainingFractionError,
            localEndpointError};
        double globalErrors[3]={0.0,0.0,0.0};
        MPI_Allreduce(localErrors,globalErrors,3,MPI_DOUBLE,MPI_MAX,
            MPI_GLOBAL_COMMUNICATOR);

        const Forcing::cSplitPlan moverSplit=Forcing::BuildSplitPlan(
            ParentReactionRates[compositionParent].TotalRate,1.0);
        const Forcing::cCarrierSplit moverCarriers=
            Forcing::SplitCarrier(moverSplit,1.0);
        const double expectedProductRate=moverCarriers.Reacted*
            compositionLikelihood*compositionProducts.size();
        moverCorrectionError=std::max(
            std::fabs(globalRates[0]-expectedProductRate)/
                std::max(expectedProductRate,
                    std::numeric_limits<double>::min()),
            std::fabs(globalRates[1]-moverCarriers.Survivor)/
                std::max(moverCarriers.Survivor,
                    std::numeric_limits<double>::min()));
        moverEventPositionError=globalErrors[0]/minimumWidth;
        moverRemainingFractionError=globalErrors[1];
        moverEndpointError=globalErrors[2]/minimumWidth;
        moverPass=moverGlobalFaceTransfers==18 &&
            moverGlobalInitialPending==18 &&
            moverGlobalContinuationMoves==18 && moverGlobalSweeps==1 &&
            moverGlobalFinalPending==0 &&
            moverGlobalCounts[0]==16 && moverGlobalCounts[1]==1 &&
            moverGlobalCounts[2]==0 && moverGlobalCounts[3]==17 &&
            moverGlobalCounts[4]==1 && moverGlobalCounts[5]==16 &&
            moverGlobalCounts[6]==0 && moverGlobalCounts[7]==0 &&
            moverGlobalCounts[8]==17 && moverGlobalCounts[9]==16 &&
            moverGlobalCounts[10]==0 && moverGlobalCounts[11]==16 &&
            moverCorrectionError<=2.0e-12 &&
            moverEventPositionError<=2.0e-12 &&
            moverRemainingFractionError<=2.0e-12 &&
            moverEndpointError<=2.0e-12;
        moverProductionPass=moverPass;

        if (validateTransportedCollision && moverPass) {
#if _PIC__PARTICLE_COLLISION_MODEL_ != _PIC__PARTICLE_COLLISION_MODEL__NTC_
          // Stage 9 is specifically a qualification of the correction-aware
          // general NTC executor used by the Halley science input.  Refuse to
          // reinterpret another collision algorithm as equivalent evidence.
          transportedCollisionPass=false;
#else
          struct cTransportedCollisionRecord {
            long int Particle=-1;
            int Species=-1;
            double Correction=0.0;
            std::uint64_t CorrectionBits=0;
            cTotalReactionForcingState State;
            double Remainder=0.0;
            double Velocity[3]={0.0,0.0,0.0};
          };
          static_assert(sizeof(double)==sizeof(std::uint64_t),
              "Stage-9 correction-bit audit requires binary64 storage");

          transportedCollisionOwner=destination->Thread;
          for (int species=0;species<PIC::nTotalSpecies;species++) {
            char symbol[100];
            PIC::MolecularData::GetChemSymbol(symbol,species);
            if (std::strcmp(symbol,"H2O")==0) {
              transportedCollisionPartner=species;
              break;
            }
          }

          std::vector<cTransportedCollisionRecord> daughters,partners;
          int collisionCell=-1,collisionI=0,collisionJ=0,collisionK=0;
          bool localCollisionPass=transportedCollisionPartner>=0 &&
              transportedCollisionPartner!=productSpecies;
          double savedProductWeight=0.0,savedProductTimeStep=0.0;
          double savedPartnerWeight=0.0,savedPartnerTimeStep=0.0;
          bool carrierChanged=false;

          // ProductEndpoints contains the pointer recorded by the application
          // callback before the core continuation mover integrates the unused
          // product interval.  Requiring those exact pointers here proves the
          // collision input was created by Stage 8 and transported by the
          // production mover; a freshly seeded C particle cannot satisfy this
          // identity test.
          if (destination->Thread==PIC::ThisThread) {
            if (destination->block==NULL) localCollisionPass=false;
            else {
              const int cellCount=
                  _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
              for (int cell=0;cell<cellCount;cell++)
                for (long int particle=
                         destination->block->FirstCellParticleTable[cell];
                     particle!=-1;
                     particle=PIC::ParticleBuffer::GetNext(particle)) {
                  bool isTransportedDaughter=false;
                  for (const cTotalReactionMoverValidationEndpoint &endpoint:
                       TotalReactionMoverValidationAudit.ProductEndpoints)
                    if (endpoint.Particle==particle) {
                      isTransportedDaughter=true;
                      break;
                    }
                  if (!isTransportedDaughter) continue;

                  PIC::ParticleBuffer::byte *data=
                      PIC::ParticleBuffer::GetParticleDataPointer(particle);
                  cTransportedCollisionRecord record;
                  record.Particle=particle;
                  record.Species=PIC::ParticleBuffer::GetI(data);
                  record.Correction=
                      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(
                          data);
                  std::memcpy(&record.CorrectionBits,&record.Correction,
                      sizeof(record.CorrectionBits));
                  record.State=GetTotalReactionState(data);
                  record.Remainder=
                      PIC::Mover::GetLocalTimeStepRemainderFraction(data);
                  PIC::ParticleBuffer::GetV(record.Velocity,data);
                  daughters.push_back(record);
                  if (collisionCell<0) collisionCell=cell;
                  else if (collisionCell!=cell) localCollisionPass=false;
                }
            }

            localCollisionPass=daughters.size()==16 && collisionCell>=0 &&
                localCollisionPass;
            if (localCollisionPass) {
              const double commonCorrection=daughters.front().Correction;
              localCollisionPass=commonCorrection>0.0 &&
                  std::isfinite(commonCorrection);
              for (const cTransportedCollisionRecord &record:daughters)
                localCollisionPass=record.Species==productSpecies &&
                    record.CorrectionBits==
                        daughters.front().CorrectionBits &&
                    record.State.Code==0 && record.Remainder==0.0 &&
                    localCollisionPass;

              collisionI=collisionCell%_BLOCK_CELLS_X_;
              const int yz=collisionCell/_BLOCK_CELLS_X_;
              collisionJ=yz%_BLOCK_CELLS_Y_;
              collisionK=yz/_BLOCK_CELLS_Y_;
              PIC::Mesh::cDataCenterNode *center=
                  destination->block->GetCenterNode(
                      _getCenterNodeLocalNumber(collisionI,collisionJ,
                          collisionK));
              localCollisionPass=center!=NULL && center->Measure>0.0 &&
                  std::isfinite(center->Measure) && localCollisionPass;

              // A counter-stream separated from every recoil product avoids
              // a zero-relative-speed fixture.  The real comet cross-section
              // adapter is evaluated for every transported daughter, so the
              // scale below remains valid if a future database makes C-H2O
              // speed dependent.
              double partnerVelocity[3]={-12000.0,1750.0,-625.0};
              double minimumSigmaCr=INFINITY;
              if (localCollisionPass)
                for (cTransportedCollisionRecord &record:daughters) {
                  const double dv0=record.Velocity[0]-partnerVelocity[0];
                  const double dv1=record.Velocity[1]-partnerVelocity[1];
                  const double dv2=record.Velocity[2]-partnerVelocity[2];
                  const double relativeSpeed=
                      std::sqrt(dv0*dv0+dv1*dv1+dv2*dv2);
                  const double crossSection=PIC::MolecularCollisions::
                      ParticleCollisionModel::UserDefined::
                          GetTotalCrossSection(record.Velocity,
                              partnerVelocity,productSpecies,
                              transportedCollisionPartner,
                              destination->block,center);
                  const double sigmaCr=crossSection*relativeSpeed;
                  if (!(sigmaCr>0.0) || !std::isfinite(sigmaCr))
                    localCollisionPass=false;
                  else minimumSigmaCr=std::min(minimumSigmaCr,sigmaCr);
                }

              constexpr int PartnerCount=64;
              constexpr double TargetTrials=256.0;
              double collisionCarrier=0.0;
              if (localCollisionPass) {
                // With W=dt=S for both species and one common correction c,
                // the unequal-species NTC prospective count is
                //
                //   N_C N_H S^2 c (sigma*g)_max / V.
                //
                // Solving with the minimum observed sigma*g guarantees at
                // least TargetTrials before collision limiting.  W/dt stays
                // exactly one, preserving the represented Stage-8 product
                // rate, while equal c, W, and dt make both update
                // probabilities analytically one.  Every accepted collision
                // then conserves represented momentum event by event.
                collisionCarrier=std::sqrt(TargetTrials*center->Measure/
                    (static_cast<double>(daughters.size())*PartnerCount*
                     commonCorrection*minimumSigmaCr));
                localCollisionPass=collisionCarrier>0.0 &&
                    std::isfinite(collisionCarrier);
              }

              if (localCollisionPass) {
                savedProductWeight=destination->block->
                    GetLocalParticleWeight(productSpecies);
                savedProductTimeStep=destination->block->
                    GetLocalTimeStep(productSpecies);
                savedPartnerWeight=destination->block->
                    GetLocalParticleWeight(transportedCollisionPartner);
                savedPartnerTimeStep=destination->block->
                    GetLocalTimeStep(transportedCollisionPartner);
                destination->block->SetLocalParticleWeight(collisionCarrier,
                    productSpecies);
                destination->block->SetLocalTimeStep(collisionCarrier,
                    productSpecies);
                destination->block->SetLocalParticleWeight(collisionCarrier,
                    transportedCollisionPartner);
                destination->block->SetLocalTimeStep(collisionCarrier,
                    transportedCollisionPartner);
                carrierChanged=true;

                double partnerPosition[3];
                PIC::ParticleBuffer::GetX(partnerPosition,
                    PIC::ParticleBuffer::GetParticleDataPointer(
                        daughters.front().Particle));
                for (int index=0;index<PartnerCount;index++) {
                  const long int particle=
                      PIC::ParticleBuffer::GetNewParticle(true);
                  PIC::ParticleBuffer::byte *data=
                      PIC::ParticleBuffer::GetParticleDataPointer(particle);
                  PIC::ParticleBuffer::SetI(transportedCollisionPartner,data);
                  PIC::ParticleBuffer::SetX(partnerPosition,data);
                  PIC::ParticleBuffer::SetV(partnerVelocity,data);
                  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
                      commonCorrection,data);
                  PIC::ParticleBuffer::SetPrev(-1,data);
                  PIC::ParticleBuffer::SetNext(-1,data);
                  PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
                  ClearTotalReactionState(data);
                  ChemistryProcessedFlag(data)=1;
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
                  PIC::IDF::InitRotTemp(180.0,data);
                  PIC::IDF::InitVibTemp(180.0,data);
#endif
                  KeepParticle(particle,destination->block->
                      FirstCellParticleTable[collisionCell]);

                  cTransportedCollisionRecord record;
                  record.Particle=particle;
                  record.Species=transportedCollisionPartner;
                  record.Correction=commonCorrection;
                  std::memcpy(&record.CorrectionBits,&record.Correction,
                      sizeof(record.CorrectionBits));
                  record.State=GetTotalReactionState(data);
                  record.Remainder=0.0;
                  PIC::ParticleBuffer::GetV(record.Velocity,data);
                  partners.push_back(record);
                }
              }
            }
          }

          int localSetupPass=localCollisionPass ? 1 : 0;
          int globalSetupPass=0;
          MPI_Allreduce(&localSetupPass,&globalSetupPass,1,MPI_INT,MPI_MIN,
              MPI_GLOBAL_COMMUNICATOR);
          localCollisionPass=globalSetupPass!=0;

          auto weightedMomentum=[](
              const std::vector<cTransportedCollisionRecord> &records,
              double result[3],double &scale) {
            for (const cTransportedCollisionRecord &record:records) {
              const double mass=PIC::MolecularData::GetMass(record.Species);
              double speedSquared=0.0;
              for (int dimension=0;dimension<3;dimension++) {
                result[dimension]+=mass*record.Correction*
                    record.Velocity[dimension];
                speedSquared+=record.Velocity[dimension]*
                    record.Velocity[dimension];
              }
              scale+=mass*record.Correction*std::sqrt(speedSquared);
            }
          };

          double localMomentumBefore[3]={0.0,0.0,0.0};
          double localMomentumScale=0.0;
          if (destination->Thread==PIC::ThisThread && localCollisionPass) {
            weightedMomentum(daughters,localMomentumBefore,
                localMomentumScale);
            weightedMomentum(partners,localMomentumBefore,
                localMomentumScale);
          }

          const auto savedCollisionSelector=PIC::MolecularCollisions::
              ParticleCollisionModel::DoSimulateCellCollisions;
          TransportedCollisionSpecies0=productSpecies;
          TransportedCollisionSpecies1=transportedCollisionPartner;
          PIC::MolecularCollisions::ParticleCollisionModel::
              DoSimulateCellCollisions=SelectTransportedCollisionPair;
          PIC::MolecularCollisions::ParticleCollisionModel::
              ResetNTCValidationDiagnostic(true);

          long long localCycles=0;
          if (destination->Thread==PIC::ThisThread && localCollisionPass) {
            // Re-seed after internal-state initialization so the collision
            // fixture, not the number of IDF random draws, owns this stream.
            rnd_seed(971927+PIC::ThisThread);
            constexpr int MaximumCollisionCycles=16;
            for (;localCycles<MaximumCollisionCycles;localCycles++) {
              PIC::MolecularCollisions::ParticleCollisionModel::
                  ModelCellCollisions_ntc(collisionI,collisionJ,collisionK,
                      destination);
              long long daughterChanges=0,partnerChanges=0;
              for (const cTransportedCollisionRecord &record:daughters) {
                double velocity[3];
                PIC::ParticleBuffer::GetV(velocity,
                    PIC::ParticleBuffer::GetParticleDataPointer(
                        record.Particle));
                if (std::memcmp(velocity,record.Velocity,
                    sizeof(velocity))!=0) daughterChanges++;
              }
              for (const cTransportedCollisionRecord &record:partners) {
                double velocity[3];
                PIC::ParticleBuffer::GetV(velocity,
                    PIC::ParticleBuffer::GetParticleDataPointer(
                        record.Particle));
                if (std::memcmp(velocity,record.Velocity,
                    sizeof(velocity))!=0) partnerChanges++;
              }
              const auto diagnostic=PIC::MolecularCollisions::
                  ParticleCollisionModel::GetNTCValidationDiagnostic();
              if (daughterChanges>=8 && partnerChanges>=8 &&
                  diagnostic.ParticleUpdateCount>=64) {
                localCycles++;
                break;
              }
            }
          }

          const auto localDiagnostic=PIC::MolecularCollisions::
              ParticleCollisionModel::GetNTCValidationDiagnostic();
          PIC::MolecularCollisions::ParticleCollisionModel::
              DoSimulateCellCollisions=savedCollisionSelector;
          TransportedCollisionSpecies0=-1;
          TransportedCollisionSpecies1=-1;
          PIC::MolecularCollisions::ParticleCollisionModel::
              ResetNTCValidationDiagnostic(false);

          long long localDaughters=0,localPartners=0;
          long long localDaughtersChanged=0,localPartnersChanged=0;
          long long localUnchangedRecords=0,localUnfinished=0;
          long long localWrongOwner=0;
          double localCorrectionError=0.0,localRateBefore=0.0;
          double localRateAfter=0.0;
          double localMomentumAfter[3]={0.0,0.0,0.0};
          double ignoredMomentumScale=0.0;
          if (destination->Thread==PIC::ThisThread && localCollisionPass) {
            localDaughters=daughters.size();
            localPartners=partners.size();
            for (const cTransportedCollisionRecord &record:daughters) {
              PIC::ParticleBuffer::byte *data=
                  PIC::ParticleBuffer::GetParticleDataPointer(
                      record.Particle);
              double velocity[3],correction;
              std::uint64_t correctionBits=0;
              PIC::ParticleBuffer::GetV(velocity,data);
              correction=PIC::ParticleBuffer::
                  GetIndividualStatWeightCorrection(data);
              std::memcpy(&correctionBits,&correction,
                  sizeof(correctionBits));
              const cTotalReactionForcingState state=
                  GetTotalReactionState(data);
              const double remainder=PIC::Mover::
                  GetLocalTimeStepRemainderFraction(data);
              if (std::memcmp(velocity,record.Velocity,
                  sizeof(velocity))!=0) localDaughtersChanged++;
              const bool unchanged=PIC::ParticleBuffer::GetI(data)==
                      record.Species &&
                  correctionBits==record.CorrectionBits &&
                  std::memcmp(&state,&record.State,sizeof(state))==0 &&
                  remainder==record.Remainder &&
                  ChemistryProcessedFlag(data)==1;
              if (unchanged) localUnchangedRecords++;
              if (state.Code!=0 || remainder!=0.0) localUnfinished++;
              localCorrectionError=std::max(localCorrectionError,
                  std::fabs(correction-record.Correction)/record.Correction);
              localRateBefore+=record.Correction;
              localRateAfter+=correction;

              cTransportedCollisionRecord after=record;
              after.Correction=correction;
              std::memcpy(after.Velocity,velocity,sizeof(velocity));
              std::vector<cTransportedCollisionRecord> one(1,after);
              weightedMomentum(one,localMomentumAfter,
                  ignoredMomentumScale);
            }
            for (const cTransportedCollisionRecord &record:partners) {
              PIC::ParticleBuffer::byte *data=
                  PIC::ParticleBuffer::GetParticleDataPointer(
                      record.Particle);
              double velocity[3],correction;
              std::uint64_t correctionBits=0;
              PIC::ParticleBuffer::GetV(velocity,data);
              correction=PIC::ParticleBuffer::
                  GetIndividualStatWeightCorrection(data);
              std::memcpy(&correctionBits,&correction,
                  sizeof(correctionBits));
              const cTotalReactionForcingState state=
                  GetTotalReactionState(data);
              const double remainder=PIC::Mover::
                  GetLocalTimeStepRemainderFraction(data);
              if (std::memcmp(velocity,record.Velocity,
                  sizeof(velocity))!=0) localPartnersChanged++;
              const bool unchanged=PIC::ParticleBuffer::GetI(data)==
                      record.Species &&
                  correctionBits==record.CorrectionBits &&
                  std::memcmp(&state,&record.State,sizeof(state))==0 &&
                  remainder==record.Remainder &&
                  ChemistryProcessedFlag(data)==1;
              if (unchanged) localUnchangedRecords++;
              if (state.Code!=0 || remainder!=0.0) localUnfinished++;
              localCorrectionError=std::max(localCorrectionError,
                  std::fabs(correction-record.Correction)/record.Correction);

              cTransportedCollisionRecord after=record;
              after.Correction=correction;
              std::memcpy(after.Velocity,velocity,sizeof(velocity));
              std::vector<cTransportedCollisionRecord> one(1,after);
              weightedMomentum(one,localMomentumAfter,
                  ignoredMomentumScale);
            }
            if (PIC::ThisThread!=destination->Thread) localWrongOwner++;
          }

          if (carrierChanged) {
            destination->block->SetLocalParticleWeight(savedProductWeight,
                productSpecies);
            destination->block->SetLocalTimeStep(savedProductTimeStep,
                productSpecies);
            destination->block->SetLocalParticleWeight(savedPartnerWeight,
                transportedCollisionPartner);
            destination->block->SetLocalTimeStep(savedPartnerTimeStep,
                transportedCollisionPartner);
          }

          long long localIntegerResults[8]={localDaughters,localPartners,
              localDaughtersChanged,localPartnersChanged,
              localUnchangedRecords,localUnfinished,localWrongOwner,
              localCycles};
          long long globalIntegerResults[8]={};
          MPI_Allreduce(localIntegerResults,globalIntegerResults,8,
              MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
          transportedDaughters=globalIntegerResults[0];
          transportedPartners=globalIntegerResults[1];
          transportedDaughtersChanged=globalIntegerResults[2];
          transportedPartnersChanged=globalIntegerResults[3];
          transportedUnchangedRecords=globalIntegerResults[4];
          transportedUnfinished=globalIntegerResults[5];
          transportedWrongOwner=globalIntegerResults[6];
          transportedCollisionCycles=globalIntegerResults[7];

          unsigned long long localDiagnosticCounts[6]={
              localDiagnostic.TrialCount,
              localDiagnostic.AcceptedCollisionCount,
              localDiagnostic.ParticleUpdateCount,
              localDiagnostic.UpdateProbabilityCount,
              localDiagnostic.InvalidUpdateProbabilityCount,
              localDiagnostic.RoundoffClampCount};
          unsigned long long globalDiagnosticCounts[6]={};
          MPI_Allreduce(localDiagnosticCounts,globalDiagnosticCounts,6,
              MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
          transportedTrials=globalDiagnosticCounts[0];
          transportedAccepted=globalDiagnosticCounts[1];
          transportedUpdates=globalDiagnosticCounts[2];
          transportedProbabilityChecks=globalDiagnosticCounts[3];
          transportedInvalidProbabilities=globalDiagnosticCounts[4];
          transportedRoundoffClamps=globalDiagnosticCounts[5];

          const double localProbabilityMinimum=
              localDiagnostic.UpdateProbabilityCount>0 ?
                  localDiagnostic.MinimumUpdateProbability : INFINITY;
          const double localProbabilityMaximum=
              localDiagnostic.UpdateProbabilityCount>0 ?
                  localDiagnostic.MaximumUpdateProbability : -INFINITY;
          MPI_Allreduce(&localProbabilityMinimum,
              &transportedProbabilityMinimum,1,MPI_DOUBLE,MPI_MIN,
              MPI_GLOBAL_COMMUNICATOR);
          MPI_Allreduce(&localProbabilityMaximum,
              &transportedProbabilityMaximum,1,MPI_DOUBLE,MPI_MAX,
              MPI_GLOBAL_COMMUNICATOR);

          double globalMomentumBefore[3]={0.0,0.0,0.0};
          double globalMomentumAfter[3]={0.0,0.0,0.0};
          double globalMomentumScale=0.0;
          MPI_Allreduce(localMomentumBefore,globalMomentumBefore,3,
              MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
          MPI_Allreduce(localMomentumAfter,globalMomentumAfter,3,
              MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
          MPI_Allreduce(&localMomentumScale,&globalMomentumScale,1,
              MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
          double momentumDifferenceSquared=0.0;
          for (int dimension=0;dimension<3;dimension++) {
            const double difference=globalMomentumAfter[dimension]-
                globalMomentumBefore[dimension];
            momentumDifferenceSquared+=difference*difference;
          }
          transportedMomentumError=std::sqrt(momentumDifferenceSquared)/
              std::max(globalMomentumScale,
                  std::numeric_limits<double>::min());
          MPI_Allreduce(&localCorrectionError,
              &transportedCorrectionError,1,MPI_DOUBLE,MPI_MAX,
              MPI_GLOBAL_COMMUNICATOR);
          double globalRatesBeforeAfter[2]={0.0,0.0};
          const double localRatesBeforeAfter[2]={localRateBefore,
              localRateAfter};
          MPI_Allreduce(localRatesBeforeAfter,globalRatesBeforeAfter,2,
              MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
          transportedRateError=std::fabs(globalRatesBeforeAfter[1]-
              globalRatesBeforeAfter[0])/
              std::max(globalRatesBeforeAfter[0],
                  std::numeric_limits<double>::min());

          localCollisionPass=localCollisionPass &&
              transportedDaughters==16 && transportedPartners==64 &&
              transportedDaughtersChanged>=8 &&
              transportedPartnersChanged>=8 &&
              transportedUnchangedRecords==80 &&
              transportedUnfinished==0 && transportedWrongOwner==0 &&
              transportedCollisionCycles>0 && transportedTrials>0 &&
              transportedAccepted>0 && transportedUpdates>=64 &&
              transportedProbabilityChecks==2*transportedAccepted &&
              transportedInvalidProbabilities==0 &&
              transportedRoundoffClamps==0 &&
              std::fabs(transportedProbabilityMinimum-1.0)<=2.0e-12 &&
              std::fabs(transportedProbabilityMaximum-1.0)<=2.0e-12 &&
              transportedCorrectionError==0.0 &&
              transportedRateError==0.0 &&
              transportedMomentumError<=2.0e-12;
          int localCollisionResult=localCollisionPass ? 1 : 0;
          int globalCollisionResult=0;
          MPI_Allreduce(&localCollisionResult,&globalCollisionResult,1,
              MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
          transportedCollisionPass=globalCollisionResult!=0;
          moverPass=moverPass && transportedCollisionPass;
#endif
        }
        else if (validateTransportedCollision) {
          transportedCollisionPass=false;
          moverPass=false;
        }
        PIC::Mover::ResetLocalTimeStepFluxDiagnostic(false);
      }
      Comet::Diagnostics::Enabled=savedDiagnostics;
      TotalReactionMoverValidationAudit.Active=false;
      PIC::ParticleBuffer::DeleteAllParticles();
      rnd_seed(static_cast<int>(Comet::Model::RandomSeed)+PIC::ThisThread);
    }
  }

  int localMoverPass=moverPass ? 1 : 0,globalMoverPass=0;
  MPI_Allreduce(&localMoverPass,&globalMoverPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  moverPass=globalMoverPass!=0;

  long long localCounts[10]={localParents,localSurvivors,localReacted,
      localEmitted,localProcessed,localRemainingMoves,localCapacityRollbacks,
      localByteRollbacks,localRngRollbacks,localPathRollbacks};
  long long globalCounts[10]={};
  MPI_Allreduce(localCounts,globalCounts,10,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  double globalCorrectionError=0.0,globalPositionError=0.0;
  MPI_Allreduce(&localCorrectionError,&globalCorrectionError,1,MPI_DOUBLE,
      MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localPositionError,&globalPositionError,1,MPI_DOUBLE,
      MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  long long globalDoubleChemistry=0;
  MPI_Allreduce(&localDoubleChemistry,&globalDoubleChemistry,1,MPI_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  int localTransactionPass=transactionPass ? 1 : 0;
  int globalTransactionPass=0;
  MPI_Allreduce(&localTransactionPass,&globalTransactionPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);

  const bool allPass=homogeneousPass && constantResult.Pass &&
      piecewiseResult.Pass && decayPass && compositionPass &&
      globalTransactionPass!=0 && moverPass;
  if (PIC::ThisThread==0) {
    std::printf("Total-reaction forcing validation: events/fixture=%ld ranks=%d deterministic-global-events=on branch-mode=deterministic-two-branch\n",
        eventsPerFixture,PIC::nTotalThreads);
    std::printf("%s total-reaction-forcing homogeneous parent=%d events=%ld steps=%d hazard=%.17e duration=%.17e optical-depth=%.17e analytic-survival=%.17e represented-survival=%.17e analytic-reacted=%.17e represented-reacted=%.17e closure-error=%.17e\n",
        homogeneousPass ? "PASS" : "FAIL",compositionParent,
        eventsPerFixture,homogeneousSteps,homogeneousHazard,
        homogeneousDuration,homogeneousTau,homogeneousAnalyticSurvival,
        homogeneousSurvival,homogeneousAnalyticReacted,
        homogeneousRepresentedReacted,homogeneousClosure);
    std::printf("%s total-reaction-forcing position fixture=constant parent=%d samples=%ld tau-first=%.17e tau-second=%.17e first-fraction=%.17e first-region-count=%ld first-region-expected=%.17e mean-fraction=%.17e expected-mean-fraction=%.17e cdf-distance=%.17e endpoint-violations=%ld\n",
        constantResult.Pass ? "PASS" : "FAIL",compositionParent,
        eventsPerFixture,constantTau,0.0,1.0,constantResult.FirstCount,
        constantResult.FirstExpected,constantResult.Mean,
        constantResult.ExpectedMean,constantResult.CdfDistance,
        constantResult.EndpointViolations);
    std::printf("%s total-reaction-forcing position fixture=piecewise parent=%d samples=%ld tau-first=%.17e tau-second=%.17e first-fraction=%.17e first-region-count=%ld first-region-expected=%.17e mean-fraction=%.17e expected-mean-fraction=%.17e cdf-distance=%.17e endpoint-violations=%ld\n",
        piecewiseResult.Pass ? "PASS" : "FAIL",compositionParent,
        eventsPerFixture,piecewiseTau1,piecewiseTau2,firstFraction,
        piecewiseResult.FirstCount,piecewiseResult.FirstExpected,
        piecewiseResult.Mean,piecewiseResult.ExpectedMean,
        piecewiseResult.CdfDistance,piecewiseResult.EndpointViolations);
    std::printf("%s total-reaction-forcing decay parent=%d events=%ld steps=%d total-optical-depth=%.17e analytic-survival=%.17e represented-survival=%.17e analytic-reacted=%.17e represented-reacted=%.17e closure-error=%.17e\n",
        decayPass ? "PASS" : "FAIL",compositionParent,eventsPerFixture,
        decaySteps,decayTotalTau,decayAnalyticSurvival,decaySurvival,
        decayAnalyticReacted,decayRepresentedReacted,decayClosure);
    std::printf("%s total-reaction-forcing composition reaction=C2_2C parent=%d packets=%ld products=%zu occurrence-factor=%.17e channel-likelihood=%.17e expected-correction=%.17e packet-correction-sum=%.17e emitted=%d momentum-error=%.17e energy-error=%.17e\n",
        compositionPass ? "PASS" : "FAIL",compositionParent,
        compositionPackets,compositionProducts.size(),occurrenceFactor,
        compositionLikelihood,expectedCorrection,packetCorrectionSum,
        16*PIC::nTotalThreads,compositionMomentumError,
        compositionEnergyError);
    std::printf("%s total-reaction-forcing production-transaction ranks=%d parents=%lld survivors=%lld reacted-branches=%lld emitted=%lld processed=%lld remaining-moves=%lld capacity-rollbacks=%lld byte-rollbacks=%lld rng-rollbacks=%lld path-rollbacks=%lld double-chemistry=%lld correction-error=%.17e position-error=%.17e\n",
        globalTransactionPass ? "PASS" : "FAIL",PIC::nTotalThreads,
        globalCounts[0],globalCounts[1],globalCounts[2],globalCounts[3],
        globalCounts[4],globalCounts[5],globalCounts[6],globalCounts[7],
        globalCounts[8],globalCounts[9],globalDoubleChemistry,
        globalCorrectionError,globalPositionError);
    std::printf("%s total-reaction-forcing mover-production ranks=%d mpi-face=%d face-transfers=%lld initial-pending=%ld continuation-moves=%ld sweeps=%d final-pending=%ld product-events=%lld sink-events=%lld wrong-owner-events=%lld final-particles=%lld survivors=%lld products=%lld uncleared-markers=%lld nonzero-remainders=%lld processed=%lld endpoint-checks=%lld moved-endpoints=%lld correction-error=%.17e event-position-error=%.17e remaining-fraction-error=%.17e endpoint-error=%.17e\n",
        moverProductionPass ? "PASS" : "FAIL",PIC::nTotalThreads,
        moverMpiFace,
        moverGlobalFaceTransfers,moverGlobalInitialPending,
        moverGlobalContinuationMoves,moverGlobalSweeps,
        moverGlobalFinalPending,moverGlobalCounts[0],moverGlobalCounts[1],
        moverGlobalCounts[2],moverGlobalCounts[3],moverGlobalCounts[4],
        moverGlobalCounts[5],moverGlobalCounts[6],moverGlobalCounts[7],
        moverGlobalCounts[8],moverGlobalCounts[9],moverGlobalCounts[11],
        moverCorrectionError,moverEventPositionError,
        moverRemainingFractionError,moverEndpointError);
    if (validateTransportedCollision)
      std::printf("%s transported-daughter-collision ranks=%d reaction=C2_2C product=%d partner=%d owner=%d daughters=%lld partners=%lld daughter-velocity-changes=%lld partner-velocity-changes=%lld unchanged-records=%lld unfinished-trajectories=%lld wrong-owner=%lld cycles=%lld trials=%llu accepted=%llu updates=%llu probability-checks=%llu invalid-probabilities=%llu roundoff-clamps=%llu probability-min=%.17e probability-max=%.17e correction-error=%.17e represented-rate-error=%.17e momentum-error=%.17e\n",
          transportedCollisionPass ? "PASS" : "FAIL",
          PIC::nTotalThreads,compositionProducts.empty() ? -1 :
              compositionProducts[0],transportedCollisionPartner,
          transportedCollisionOwner,transportedDaughters,
          transportedPartners,transportedDaughtersChanged,
          transportedPartnersChanged,transportedUnchangedRecords,
          transportedUnfinished,transportedWrongOwner,
          transportedCollisionCycles,transportedTrials,
          transportedAccepted,transportedUpdates,
          transportedProbabilityChecks,transportedInvalidProbabilities,
          transportedRoundoffClamps,
          transportedProbabilityMinimum,transportedProbabilityMaximum,
          transportedCorrectionError,transportedRateError,
          transportedMomentumError);
    std::printf("Total-reaction forcing validation result: %s homogeneous=%d piecewise=%d decay=%d composition=%d transactions=%d mover=%d events=%ld\n",
        allPass ? "PASS" : "FAIL",homogeneousPass ? 1 : 0,
        piecewiseResult.Pass ? 1 : 0,decayPass ? 1 : 0,
        compositionPass ? 1 : 0,globalTransactionPass ? PIC::nTotalThreads : 0,
        moverPass ? 1 : 0,eventsPerFixture);
  }
  return allPass;
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
           ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        ChemistryProcessedFlag(data)=0;
        if (TotalReactionForcingEnabled) {
          // Chemistry is followed by a complete mover plus all local-timestep
          // continuation sweeps. Consequently no legitimate pending event may
          // cross into the next chemistry stage. Clearing every slot here also
          // initializes application bytes in freshly injected particles and
          // in reused buffer slots without adding a generic core dependency.
          ClearTotalReactionState(data);
        }
      }
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

namespace {

bool ProcessTotalReactionForcedParticle(long int ptr,
    long int &firstParticleCell,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int parentSpecies,
    const cParentReactionRateCache &parentRates,double timeStep) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  (void)ptr;
  (void)firstParticleCell;
  (void)node;
  (void)parentSpecies;
  (void)parentRates;
  (void)timeStep;
  std::fprintf(stderr,
      "Stage-8 total-reaction forcing requires individual corrections\n");
  std::abort();
#else
  namespace Forcing=Comet::TotalReactionForcingPlanner;
  PIC::ParticleBuffer::byte *parentData=
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  const double parentCorrection=
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(parentData);
  if (!(parentCorrection>0.0) || !std::isfinite(parentCorrection)) {
    std::fprintf(stderr,"invalid Stage-8 incoming parent correction\n");
    std::abort();
  }

  Forcing::cSplitPlan split;
  Forcing::cCarrierSplit carriers;
  try {
    split=Forcing::BuildSplitPlan(parentRates.TotalRate,timeStep);
    carriers=Forcing::SplitCarrier(split,parentCorrection);
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,"cannot plan Stage-8 parent split: %s\n",
        exception.what());
    std::abort();
  }
  if (!(carriers.Reacted>0.0)) {
    // A positive rate can still round to a zero reacted fraction only below
    // the representable optical-depth range. Preserve the exact incoming
    // carrier as a survivor and consume no channel/event RNG in that limit.
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(carriers.Survivor,
        parentData);
    ClearTotalReactionState(parentData);
    MarkChemistryProcessed(parentData);
    KeepParticle(ptr,firstParticleCell);
    return true;
  }

  const ChemicalSolver::cReactionNetwork &network=
      Comet::Chemistry::GetNetwork();

  // Preflight the largest channel in this parent family before the first
  // channel, position, or recoil-key draw.  A capacity deferral therefore
  // leaves the parent bytes, linked-list path, physical selection/product
  // diagnostics, and global RNG stream exactly unchanged.  The outer driver
  // has already counted an eligible visit, and a rejected preflight records
  // one allocation deferral; those are work/pressure diagnostics rather than
  // a committed reaction. Product-free sinks still need one temporary parent
  // marker so their loss occurs at the physical in-flight position.
  std::size_t maximumMarkers=0;
  for (std::size_t channel:parentRates.Channels) {
    std::vector<int> products;
    ExpandConcreteProducts(channel,products);
    const std::size_t packetCount=static_cast<std::size_t>(
        Comet::Chemistry::GetDaughterPacketCount(channel));
    std::size_t markerCount=1;
    if (!products.empty()) {
      if (packetCount>std::numeric_limits<std::size_t>::max()/products.size()) {
        std::fprintf(stderr,"Stage-8 marker count overflows size_t\n");
        std::abort();
      }
      markerCount=packetCount*products.size();
    }
    maximumMarkers=std::max(maximumMarkers,markerCount);
  }
  if (maximumMarkers==0 || maximumMarkers>static_cast<std::size_t>(
      Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction)) {
    std::fprintf(stderr,
        "Stage-8 family exceeds the photochemical transaction limit\n");
    std::abort();
  }
  const long int transactionCeiling=static_cast<long int>(std::floor(
      (1.0-Comet::RadialSimulation::AllocationReserveFraction)*
      PIC::ParticleBuffer::MaxNPart));
  const std::size_t maximumNewSlots=maximumMarkers+
      (carriers.Survivor>0.0 ? 1U : 0U)-1U;
  const long int availableNewSlots=
      transactionCeiling-PIC::ParticleBuffer::NAllPart;
  if (availableNewSlots<0 || maximumNewSlots>
      static_cast<std::size_t>(availableNewSlots)) {
    LocalChemistryAllocationPressure=true;
    if (Comet::Diagnostics::Enabled)
      Comet::Diagnostics::RecordAllocationDeferral(
          parentRates.Channels.front());
    KeepParticle(ptr,firstParticleCell);
    return false;
  }

  std::size_t selectedSlot=parentRates.Channels.size()-1;
  if (!Comet::Chemistry::IsChannelImportanceSamplingActive() ||
      !parentRates.ImportanceConfigured) {
    double selector=rnd()*parentRates.TotalRate;
    for (std::size_t slot=0;slot<parentRates.Channels.size();slot++) {
      selector-=parentRates.Rates[slot];
      if (selector<=0.0) {
        selectedSlot=slot;
        break;
      }
    }
  }
  else {
    double selector=rnd();
    for (std::size_t slot=0;slot<parentRates.Channels.size();slot++) {
      selector-=parentRates.SelectionProbabilities[slot];
      if (selector<=0.0) {
        selectedSlot=slot;
        break;
      }
    }
  }
  const std::size_t selected=parentRates.Channels[selectedSlot];
  const double likelihood=
      Comet::Chemistry::IsChannelImportanceSamplingActive() ?
          parentRates.LikelihoodRatios[selectedSlot] : 1.0;
  if (!(likelihood>0.0) || !std::isfinite(likelihood)) {
    std::fprintf(stderr,"invalid Stage-8 channel likelihood\n");
    std::abort();
  }

  std::vector<int> products;
  ExpandConcreteProducts(selected,products);
  const long int packetCount=
      Comet::Chemistry::GetDaughterPacketCount(selected);
  if (products.empty() && packetCount!=1) {
    std::fprintf(stderr,"Stage-8 sink cannot use daughter packets\n");
    std::abort();
  }
  const std::size_t markerCount=products.empty() ? 1 :
      products.size()*static_cast<std::size_t>(packetCount);

  double unit=rnd();
  if (!(unit>=0.0) || !std::isfinite(unit)) {
    std::fprintf(stderr,"invalid Stage-8 conditional-position draw\n");
    std::abort();
  }
  if (unit>=1.0) unit=std::nextafter(1.0,0.0);
  double eventFraction;
  try {
    eventFraction=Forcing::SampleConditionalReactionTime(split,unit)/timeStep;
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,"cannot sample Stage-8 event position: %s\n",
        exception.what());
    std::abort();
  }
  if (eventFraction<0.0 || eventFraction>1.0 ||
      !std::isfinite(eventFraction)) {
    std::fprintf(stderr,"invalid normalized Stage-8 event position\n");
    std::abort();
  }

  // Build per-product packet corrections while the marker still has the
  // parent's W/dt carrier. For each concrete product, packet corrections sum
  // exactly to P*C/K. The event callback multiplies by K, L, and the carrier
  // ratio, yielding P*C*L/N per packet and P*C*L per concrete product.
  Comet::DaughterPacketPlanner::cPlan markerPlan;
  if (!products.empty()) {
    std::vector<double> unsplit(products.size(),
        carriers.Reacted/static_cast<double>(products.size()));
    try {
      markerPlan=Comet::DaughterPacketPlanner::BuildPlan(unsplit,
          static_cast<std::size_t>(packetCount),
          static_cast<std::size_t>(
              std::max<long int>(1,Comet::Chemistry::MaximumDaughterPacketCount)),
          static_cast<std::size_t>(Comet::RadialSimulation::
              MaximumPhotochemicalProductsPerReaction));
    }
    catch (const std::exception &exception) {
      std::fprintf(stderr,"cannot plan Stage-8 packet markers: %s\n",
          exception.what());
      std::abort();
    }
    if (markerPlan.TotalDaughterCount!=markerCount) {
      std::fprintf(stderr,"Stage-8 marker-plan count mismatch\n");
      std::abort();
    }
  }

  char markerTemplate[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated(
      reinterpret_cast<PIC::ParticleBuffer::byte*>(markerTemplate));
  PIC::ParticleBuffer::SetPrev(-1,
      reinterpret_cast<PIC::ParticleBuffer::byte*>(markerTemplate));
  PIC::ParticleBuffer::SetNext(-1,
      reinterpret_cast<PIC::ParticleBuffer::byte*>(markerTemplate));
  PIC::ParticleBuffer::CloneParticle(
      reinterpret_cast<PIC::ParticleBuffer::byte*>(markerTemplate),parentData);

  bool parentSlotUsed=false;
  if (carriers.Survivor>0.0) {
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(carriers.Survivor,
        parentData);
    ClearTotalReactionState(parentData);
    MarkChemistryProcessed(parentData);
    KeepParticle(ptr,firstParticleCell);
    parentSlotUsed=true;
  }

  auto emitMarker=[&](const cTotalReactionForcingState &state,
      double correction) {
    if (!(correction>0.0) || !std::isfinite(correction)) {
      std::fprintf(stderr,"invalid Stage-8 pending-marker correction\n");
      std::abort();
    }
    PIC::ParticleBuffer::byte *templateData=
        reinterpret_cast<PIC::ParticleBuffer::byte*>(markerTemplate);
    PIC::ParticleBuffer::SetI(parentSpecies,templateData);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(correction,
        templateData);
    SetTotalReactionState(templateData,state);
    MarkChemistryProcessed(templateData);

    if (!parentSlotUsed) {
      PIC::ParticleBuffer::SetPrev(-1,parentData);
      PIC::ParticleBuffer::SetNext(-1,parentData);
      PIC::ParticleBuffer::CloneParticle(parentData,templateData);
      KeepParticle(ptr,firstParticleCell);
      parentSlotUsed=true;
    }
    else {
      const long int marker=
          PIC::ParticleBuffer::GetNewParticle(firstParticleCell);
      PIC::ParticleBuffer::CloneParticle(
          PIC::ParticleBuffer::GetParticleDataPointer(marker),templateData);
    }
  };

  if (products.empty()) {
    cTotalReactionForcingState state;
    state.Code=EncodeTotalReactionSink(selected);
    state.KeyAndFlags=TotalReactionLeaderFlag;
    state.RemainingFraction=eventFraction;
    emitMarker(state,carriers.Reacted);
  }
  else {
    for (long int packet=0;packet<packetCount;packet++) {
      double keyDraw=rnd();
      if (!(keyDraw>=0.0) || !std::isfinite(keyDraw)) {
        std::fprintf(stderr,"invalid Stage-8 recoil-key draw\n");
        std::abort();
      }
      if (keyDraw>=1.0) keyDraw=std::nextafter(1.0,0.0);
      const std::uint32_t key=static_cast<std::uint32_t>(
          std::floor(keyDraw*2147483648.0))&TotalReactionSeedMask;
      for (std::size_t product=0;product<products.size();product++) {
        cTotalReactionForcingState state;
        state.Code=EncodeTotalReactionState(selected,product);
        state.KeyAndFlags=key;
        if (packet==0) state.KeyAndFlags|=TotalReactionLeaderFlag;
        state.RemainingFraction=eventFraction;
        emitMarker(state,markerPlan.Correction(packet,product));
      }
    }
  }
  return true;
#endif
}

} // namespace

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
  if (IsTotalReactionForcingConfigured(parentSpecies))
    return ProcessTotalReactionForcedParticle(ptr,firstParticleCell,node,
        parentSpecies,parentRates,timeStep);
  if (rnd()>=1.0-std::exp(-totalRate*timeStep)) {
    MarkChemistryProcessed(parentData);
    KeepParticle(ptr,firstParticleCell);
    return true;
  }

  std::size_t selectedSlot=parentRates.Channels.size()-1;
  std::size_t selected=parentRates.Channels.back();
  if (!IsChannelImportanceSamplingActive() ||
      !parentRates.ImportanceConfigured) {
    // Compatibility boundary: retain the exact historical multiply and
    // left-to-right rate subtraction for every analog family.  Enabling a
    // policy for another parent therefore does not perturb this family's
    // selection thresholds or random-number consumption.
    double selector=rnd()*totalRate;
    for (std::size_t i=0;i<parentRates.Channels.size();i++) {
      selector-=parentRates.Rates[i];
      if (selector<=0.0) {
        selectedSlot=i;
        selected=parentRates.Channels[i];
        break;
      }
    }
  }
  else {
    double selector=rnd();
    for (std::size_t i=0;i<parentRates.Channels.size();i++) {
      selector-=parentRates.SelectionProbabilities[i];
      if (selector<=0.0) {
        selectedSlot=i;
        selected=parentRates.Channels[i];
        break;
      }
    }
  }
  const double channelLikelihood=IsChannelImportanceSamplingActive() ?
      parentRates.LikelihoodRatios[selectedSlot] : 1.0;
  if (!(channelLikelihood>0.0) || !std::isfinite(channelLikelihood)) {
    std::fprintf(stderr,
        "invalid channel-importance likelihood for reaction %zu\n",selected);
    std::abort();
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
        parentPhysicalWeight*productTimeStep/(baseWeight*timeStep)*
        channelLikelihood;
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
        reaction.ProductCount==0,channelLikelihood,
        parentRate*channelLikelihood);
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
          parentRate*channelLikelihood/
          static_cast<double>(reactionPacketCount);
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

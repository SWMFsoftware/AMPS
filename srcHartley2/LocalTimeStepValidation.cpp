#include "LocalTimeStepValidation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "pic.h"

namespace {

using cNode=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;

constexpr double RequiredRelativeAccuracy=0.02;
constexpr double MaximumStandardDeviations=5.0;
constexpr double CrossingFraction=0.25;

// Exercise H2O whenever the configured species table provides it.  The Halley
// case gives water rotational degrees of freedom, so the same face tests also
// verify that cloning and MPI continuation preserve internal-energy state.
// The fallback keeps reduced, non-comet species tables buildable.
#ifdef _H2O_SPEC_
constexpr int TestSpecies=_H2O_SPEC_;
#else
constexpr int TestSpecies=0;
#endif

// Exact powers of two give a wide, bitwise reproducible correction cycle.  A
// copied or exchanged particle must therefore return one of these exact bins.
// Correction-off builds retain the historical c=1 behavior through the same
// code path instead of maintaining a separate, weaker validator.
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
const double ValidationCorrectionCycle[]={1.0/16.0,1.0/4.0,1.0,4.0,16.0};
#else
const double ValidationCorrectionCycle[]={1.0};
#endif
constexpr int ValidationCorrectionCount=
    sizeof(ValidationCorrectionCycle)/sizeof(ValidationCorrectionCycle[0]);
constexpr double ValidationInternalEnergyQuantum=1.0e-21;
bool UseUnitCorrections=false;

struct cFacePair {
  cNode *Left=NULL;
  cNode *Right=NULL;
  bool IsValid() const {return Left!=NULL && Right!=NULL;}
};

struct cParticleScan {
  long long Count=0;
  long long WrongNode=0;
  long long NonzeroRemainder=0;
  long long InvalidState=0;
  long long InvalidSpeed=0;
  long long InvalidRadialProbe=0;
  long long InvalidCorrection=0;
  long long InvalidInternalEnergy=0;
  long long CorrectionCounts[ValidationCorrectionCount]={};
  double CorrectionSum=0.0;
  double CorrectionSquareSum=0.0;
  double CorrectionWeightedInternalEnergy=0.0;
};

struct cDiagnostic {
  long long FaceTransferCount=0,DestinationParticleCount=0;
  long long RouletteDeletionCount=0,CloneCount=0;
  double IncomingNumberFlux=0.0,OutgoingNumberFlux=0.0;
  double NumberFluxVariance=0.0;
  double IncomingMassFlux=0.0,OutgoingMassFlux=0.0;
  double MassFluxVariance=0.0;
  double IncomingMomentumFlux[3]={0.0,0.0,0.0};
  double OutgoingMomentumFlux[3]={0.0,0.0,0.0};
  double MomentumFluxVariance[3]={0.0,0.0,0.0};
  double IncomingKineticEnergyFlux=0.0,OutgoingKineticEnergyFlux=0.0;
  double KineticEnergyFluxVariance=0.0;
  double IncomingInternalEnergyFlux=0.0,OutgoingInternalEnergyFlux=0.0;
  double InternalEnergyFluxVariance=0.0;
  double MaximumRemainderInvariantError=0.0;
  int ContinuationSweepCount=0;
  long int MaximumPendingParticleCount=0,FinalPendingParticleCount=0;
  long int InitialPendingParticleCount=0;
  long int FirstPostExchangePendingParticleCount=0;
  long int ContinuationMovedParticleCount=0;
  long int ContinuationDeletedParticleCount=0;
};

struct cValidationSummary {
  long long LocalFaceChecks=0;
  long long MpiFaceChecks=0;
  long long MultiFaceChecks=0;
  long long CorrectionPreservationChecks=0;
  long long CloneCorrectionChecks=0;
  long long InternalEnergyChecks=0;
  long long PartialMutations=0;
  double MaximumCorrectionError=0.0;
  double MaximumFluxChecksumError=0.0;
};

int UnitCorrectionIndex() {
  for (int index=0;index<ValidationCorrectionCount;index++)
    if (ValidationCorrectionCycle[index]==1.0) return index;
  return -1;
}

int ValidationCorrectionIndexForParticle(long int particle) {
  if (UseUnitCorrections) return UnitCorrectionIndex();
  return particle%ValidationCorrectionCount;
}

int FindValidationCorrection(double correction) {
  for (int index=0;index<ValidationCorrectionCount;index++)
    if (correction==ValidationCorrectionCycle[index]) return index;
  return -1;
}

double ParticleInternalEnergy(PIC::ParticleBuffer::byte *data,int species) {
  double energy=0.0;
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  if (PIC::IDF::nTotalRotationalModes[species]>0)
    energy+=PIC::IDF::GetRotE(data);
  for (int mode=0;mode<PIC::IDF::nTotalVibtationalModes[species];mode++)
    energy+=PIC::IDF::GetVibE(mode,data);
#else
  (void)data;
  (void)species;
#endif
  return energy;
}

bool HasDeterministicValidationInternalEnergy() {
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_ && \
    _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_ == \
        _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_
  return true;
#else
  return false;
#endif
}

bool HasValidationInternalEnergyModes() {
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  return PIC::IDF::nTotalRotationalModes[TestSpecies]>0 ||
      PIC::IDF::nTotalVibtationalModes[TestSpecies]>0;
#else
  return false;
#endif
}

double ExpectedValidationInternalEnergy(int correctionIndex,int species) {
  double energy=0.0;
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_ && \
    _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_ == \
        _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_
  if (PIC::IDF::nTotalRotationalModes[species]>0)
    energy+=ValidationInternalEnergyQuantum*(correctionIndex+1);
  for (int mode=0;mode<PIC::IDF::nTotalVibtationalModes[species];mode++)
    energy+=0.125*ValidationInternalEnergyQuantum*
        (correctionIndex+mode+1);
#else
  (void)correctionIndex;
  (void)species;
#endif
  return energy;
}

void InitializeValidationParticleState(long int particle,
    PIC::ParticleBuffer::byte *data) {
  const int correctionIndex=ValidationCorrectionIndexForParticle(particle);
  if (correctionIndex<0) {
    std::fprintf(stderr,"validation correction cycle does not contain unity\n");
    std::abort();
  }
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
      ValidationCorrectionCycle[correctionIndex],data);
#endif

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  // Initialize all model-owned fields first.  In the validated LB build the
  // random values are then replaced by an exact correction-indexed pattern,
  // making a state-copy error distinguishable from a merely plausible value.
  PIC::IDF::InitRotTemp(180.0,data);
  PIC::IDF::InitVibTemp(180.0,data);
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_ == \
    _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_
  if (PIC::IDF::nTotalRotationalModes[TestSpecies]>0)
    PIC::IDF::LB::SetRotE(
        ValidationInternalEnergyQuantum*(correctionIndex+1),data);
  for (int mode=0;mode<PIC::IDF::nTotalVibtationalModes[TestSpecies];mode++)
    PIC::IDF::LB::SetVibE(0.125*ValidationInternalEnergyQuantum*
        (correctionIndex+mode+1),mode,data);
#endif
#else
  (void)correctionIndex;
  (void)data;
#endif
}

double GetOwnerCarrier(cNode *node) {
  const double local=(node!=NULL && node->Thread==PIC::ThisThread &&
      node->block!=NULL) ?
      node->block->GetLocalParticleWeight(TestSpecies)/
          node->block->GetLocalTimeStep(TestSpecies) : 0.0;
  double global=0.0;
  MPI_Allreduce(&local,&global,1,MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  return global;
}

std::vector<cNode*> GetOrderedLeaves() {
  std::vector<cNode*> leaves;
  for (cNode *node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;
       node=node->nextBranchBottomNode)
    if (node->IsUsedInCalculationFlag) leaves.push_back(node);

  std::sort(leaves.begin(),leaves.end(),[](cNode *left,cNode *right) {
    return left->xmin[0]<right->xmin[0];
  });
  return leaves;
}

bool AdjacentEqualWidth(cNode *left,cNode *right) {
  const double leftWidth=left->xmax[0]-left->xmin[0];
  const double rightWidth=right->xmax[0]-right->xmin[0];
  const double scale=std::max({std::fabs(left->xmax[0]),leftWidth,rightWidth,1.0});
  return std::fabs(left->xmax[0]-right->xmin[0])<=1.0e-12*scale &&
      std::fabs(leftWidth-rightWidth)<=1.0e-12*scale;
}

cFacePair FindFacePair(const std::vector<cNode*> &leaves,bool differentOwner) {
  for (std::size_t index=1;index<leaves.size();index++) {
    cNode *left=leaves[index-1],*right=leaves[index];
    if (!AdjacentEqualWidth(left,right)) continue;
    if ((left->Thread!=right->Thread)!=differentOwner) continue;
    return {left,right};
  }
  return {};
}

void SetFaceTransportParameters(const std::vector<cNode*> &leaves,
    cNode *source,cNode *destination,double destinationTimeStepRatio,
    double destinationWeightRatio) {
  // Boundary-layer copies are present only on ranks that need them.  Setting
  // every allocated copy from global node identity keeps the source and
  // destination values identical on both sides of an MPI ownership face.  A
  // unit source value makes the configured destination values equal to the
  // ratios printed by the regression.
  for (cNode *node:leaves) if (node->block!=NULL) {
    double timeStep=1.0;
    if (node==source) timeStep=1.0;
    if (node==destination) timeStep=destinationTimeStepRatio;
    node->block->SetLocalTimeStep(timeStep,TestSpecies);

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
    double weight=1.0;
    if (node==source) weight=1.0;
    if (node==destination) weight=destinationWeightRatio;
    node->block->SetLocalParticleWeight(weight,TestSpecies);
#else
    // A spatially global build cannot represent a non-unit block-to-block
    // weight ratio.  Keep the older timestep-only validation usable, but fail
    // loudly if a caller accidentally requests a local-weight case.
    if (destinationWeightRatio!=1.0)
      exit(__LINE__,__FILE__,
          "local-weight face validation requires SpeciesLocalParticleWeight");
#endif
  }
}

std::array<std::array<double,3>,4> BuildVelocities(double width,int direction) {
  const double scale=width;
  return {{{direction*0.040*scale,0.0,0.0},
           {direction*0.045*scale,0.004*scale,0.0},
           {direction*0.050*scale,-0.003*scale,0.002*scale},
           {direction*0.055*scale,0.002*scale,-0.005*scale}}};
}

double Speed(const std::array<double,3> &velocity) {
  return std::sqrt(velocity[0]*velocity[0]+velocity[1]*velocity[1]+
      velocity[2]*velocity[2]);
}

void InjectSingleFaceParticles(cNode *source,double face,long int particleNumber,
    const std::array<std::array<double,3>,4> &velocities) {
  if (source->Thread!=PIC::ThisThread) return;
  if (source->block==NULL)
    exit(__LINE__,__FILE__,"validation source block is not allocated");

  for (long int particle=0;particle<particleNumber;particle++) {
    const std::array<double,3> &velocity=velocities[particle%velocities.size()];
    const double flightTime=CrossingFraction;
    const double transverse2=velocity[1]*velocity[1]+
        velocity[2]*velocity[2];
    const double radicand=face*face-transverse2*flightTime*flightTime;
    if (!(radicand>0.0))
      exit(__LINE__,__FILE__,"invalid validation face geometry");

    // At t=CrossingFraction the Cartesian trajectory lies exactly on the
    // spherical interface.  The sign of vx places the initial radius on the
    // requested source side for both outward and inward tests.
    double position[3]={
        std::sqrt(radicand)-velocity[0]*flightTime,0.0,0.0};
    double velocityBuffer[3]={velocity[0],velocity[1],velocity[2]};
    int species=TestSpecies;
    const long int ptr=PIC::ParticleBuffer::InitiateParticle(position,
        velocityBuffer,NULL,&species,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,
        source);
    PIC::ParticleBuffer::byte *data=
        PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    InitializeValidationParticleState(particle,data);
    PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
  }
}

cDiagnostic ReduceDiagnostic() {
  const PIC::Mover::cLocalTimeStepFluxDiagnostic &local=
      PIC::Mover::LocalTimeStepFluxDiagnostic;
  cDiagnostic result;
  long long localCounts[4]={local.FaceTransferCount,
      local.DestinationParticleCount,local.RouletteDeletionCount,
      local.CloneCount};
  long long globalCounts[4]={0,0,0,0};
  MPI_Allreduce(localCounts,globalCounts,4,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  result.FaceTransferCount=globalCounts[0];
  result.DestinationParticleCount=globalCounts[1];
  result.RouletteDeletionCount=globalCounts[2];
  result.CloneCount=globalCounts[3];

  const double localValues[21]={
      local.IncomingNumberFlux,local.OutgoingNumberFlux,
      local.NumberFluxVariance,local.IncomingMassFlux,
      local.OutgoingMassFlux,local.MassFluxVariance,
      local.IncomingMomentumFlux[0],local.IncomingMomentumFlux[1],
      local.IncomingMomentumFlux[2],local.OutgoingMomentumFlux[0],
      local.OutgoingMomentumFlux[1],local.OutgoingMomentumFlux[2],
      local.MomentumFluxVariance[0],local.MomentumFluxVariance[1],
      local.MomentumFluxVariance[2],local.IncomingKineticEnergyFlux,
      local.OutgoingKineticEnergyFlux,local.KineticEnergyFluxVariance,
      local.IncomingInternalEnergyFlux,local.OutgoingInternalEnergyFlux,
      local.InternalEnergyFluxVariance};
  double globalValues[21]={};
  MPI_Allreduce(localValues,globalValues,21,MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  result.IncomingNumberFlux=globalValues[0];
  result.OutgoingNumberFlux=globalValues[1];
  result.NumberFluxVariance=globalValues[2];
  result.IncomingMassFlux=globalValues[3];
  result.OutgoingMassFlux=globalValues[4];
  result.MassFluxVariance=globalValues[5];
  for (int dimension=0;dimension<3;dimension++) {
    result.IncomingMomentumFlux[dimension]=globalValues[6+dimension];
    result.OutgoingMomentumFlux[dimension]=globalValues[9+dimension];
    result.MomentumFluxVariance[dimension]=globalValues[12+dimension];
  }
  result.IncomingKineticEnergyFlux=globalValues[15];
  result.OutgoingKineticEnergyFlux=globalValues[16];
  result.KineticEnergyFluxVariance=globalValues[17];
  result.IncomingInternalEnergyFlux=globalValues[18];
  result.OutgoingInternalEnergyFlux=globalValues[19];
  result.InternalEnergyFluxVariance=globalValues[20];

  MPI_Allreduce(&local.MaximumRemainderInvariantError,
      &result.MaximumRemainderInvariantError,1,MPI_DOUBLE,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.ContinuationSweepCount,
      &result.ContinuationSweepCount,1,MPI_INT,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  const long int localGlobalLifecycle[2]={local.InitialPendingParticleCount,
      local.FirstPostExchangePendingParticleCount};
  long int globalLifecycle[2]={0,0};
  MPI_Allreduce(localGlobalLifecycle,globalLifecycle,2,MPI_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  result.InitialPendingParticleCount=globalLifecycle[0];
  result.FirstPostExchangePendingParticleCount=globalLifecycle[1];

  const long int localEventCounts[2]={local.ContinuationMovedParticleCount,
      local.ContinuationDeletedParticleCount};
  long int globalEventCounts[2]={0,0};
  MPI_Allreduce(localEventCounts,globalEventCounts,2,MPI_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  result.ContinuationMovedParticleCount=globalEventCounts[0];
  result.ContinuationDeletedParticleCount=globalEventCounts[1];
  MPI_Allreduce(&local.MaximumPendingParticleCount,
      &result.MaximumPendingParticleCount,1,MPI_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.FinalPendingParticleCount,
      &result.FinalPendingParticleCount,1,MPI_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  return result;
}

bool IsExpectedSpeed(double speed,
    const std::array<std::array<double,3>,4> &velocities,int *group=NULL) {
  for (std::size_t index=0;index<velocities.size();index++) {
    const double expected=Speed(velocities[index]);
    if (std::fabs(speed-expected)<=1.0e-11*std::max(expected,1.0)) {
      if (group!=NULL) *group=static_cast<int>(index);
      return true;
    }
  }
  return false;
}

cParticleScan ScanParticles(cNode *expectedNode,
    const std::array<std::array<double,3>,4> *velocities,
    double expectedRadialProbe) {
  cParticleScan local;
  const int cellCount=
      _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
  const double radialTolerance=std::max(256.0*PIC::Mesh::mesh->EPS,
      1.0e-10*(expectedNode->xmax[0]-expectedNode->xmin[0]));

  for (int blockIndex=0;
       blockIndex<static_cast<int>(PIC::DomainBlockDecomposition::nLocalBlocks);
       blockIndex++) {
    cNode *node=PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node->block==NULL) continue;
    for (int cell=0;cell<cellCount;cell++) {
      for (long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;
           ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        local.Count++;
        if (node!=expectedNode) local.WrongNode++;
        PIC::ParticleBuffer::byte *data=
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        const double remainder=
            PIC::Mover::GetLocalTimeStepRemainderFraction(data);
        if (remainder!=0.0) local.NonzeroRemainder++;

        double position[3],velocity[3];
        PIC::ParticleBuffer::GetX(position,data);
        PIC::ParticleBuffer::GetV(velocity,data);
        bool finite=true;
        for (int dimension=0;dimension<3;dimension++)
          finite=finite && std::isfinite(position[dimension]) &&
              std::isfinite(velocity[dimension]);
        if (!finite) {
          local.InvalidState++;
          continue;
        }

        // This compact histogram is a correction-preservation checksum.  Exact
        // bin membership catches a transformed c, while sum(c), sum(c^2), and
        // sum(c*E_internal) expose record corruption across MPI reductions.
        const double correction=
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
        const int correctionIndex=FindValidationCorrection(correction);
        const double internalEnergy=ParticleInternalEnergy(data,TestSpecies);
        if (!(correction>0.0) || !std::isfinite(correction) ||
            correctionIndex<0) {
          local.InvalidCorrection++;
        }
        else {
          local.CorrectionCounts[correctionIndex]++;
          local.CorrectionSum+=correction;
          local.CorrectionSquareSum+=correction*correction;
        }

        if (!(internalEnergy>=0.0) || !std::isfinite(internalEnergy)) {
          local.InvalidInternalEnergy++;
        }
        else {
          if (correctionIndex>=0)
            local.CorrectionWeightedInternalEnergy+=
                correction*internalEnergy;
          if (correctionIndex>=0 &&
              HasDeterministicValidationInternalEnergy()) {
            const double expected=ExpectedValidationInternalEnergy(
                correctionIndex,TestSpecies);
            const double tolerance=512.0*
                std::numeric_limits<double>::epsilon()*
                std::max(expected,std::numeric_limits<double>::min());
            if (std::fabs(internalEnergy-expected)>tolerance)
              local.InvalidInternalEnergy++;
          }
        }

        if (velocities!=NULL) {
          int group=-1;
          const double speed=std::sqrt(velocity[0]*velocity[0]+
              velocity[1]*velocity[1]+velocity[2]*velocity[2]);
          if (!IsExpectedSpeed(speed,*velocities,&group))
            local.InvalidSpeed++;
          if (group==0 && std::fabs(position[0]-expectedRadialProbe)>
              radialTolerance)
            local.InvalidRadialProbe++;
        }
      }
    }
  }

  cParticleScan global;
  long long localValues[8]={local.Count,local.WrongNode,
      local.NonzeroRemainder,local.InvalidState,local.InvalidSpeed,
      local.InvalidRadialProbe,local.InvalidCorrection,
      local.InvalidInternalEnergy};
  long long globalValues[8]={};
  MPI_Allreduce(localValues,globalValues,8,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  global.Count=globalValues[0];
  global.WrongNode=globalValues[1];
  global.NonzeroRemainder=globalValues[2];
  global.InvalidState=globalValues[3];
  global.InvalidSpeed=globalValues[4];
  global.InvalidRadialProbe=globalValues[5];
  global.InvalidCorrection=globalValues[6];
  global.InvalidInternalEnergy=globalValues[7];
  MPI_Allreduce(local.CorrectionCounts,global.CorrectionCounts,
      ValidationCorrectionCount,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  const double localChecksums[3]={local.CorrectionSum,
      local.CorrectionSquareSum,local.CorrectionWeightedInternalEnergy};
  double globalChecksums[3]={0.0,0.0,0.0};
  MPI_Allreduce(localChecksums,globalChecksums,3,MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  global.CorrectionSum=globalChecksums[0];
  global.CorrectionSquareSum=globalChecksums[1];
  global.CorrectionWeightedInternalEnergy=globalChecksums[2];
  return global;
}

bool CheckCorrectionChecksum(const cParticleScan &scan,long int sourceCount,
    double meanMultiplicity,double &maximumError,
    bool &cloneCheckPerformed) {
  bool pass=scan.InvalidCorrection==0;
  long long checksumCount=0;
  double checksumSum=0.0,checksumSquareSum=0.0;
  for (int index=0;index<ValidationCorrectionCount;index++) {
    checksumCount+=scan.CorrectionCounts[index];
    checksumSum+=scan.CorrectionCounts[index]*
        ValidationCorrectionCycle[index];
    checksumSquareSum+=scan.CorrectionCounts[index]*
        ValidationCorrectionCycle[index]*ValidationCorrectionCycle[index];
  }

  const double sumScale=std::max(std::fabs(scan.CorrectionSum),1.0);
  const double squareScale=
      std::max(std::fabs(scan.CorrectionSquareSum),1.0);
  const double sumError=std::fabs(checksumSum-scan.CorrectionSum)/sumScale;
  const double squareError=
      std::fabs(checksumSquareSum-scan.CorrectionSquareSum)/squareScale;
  maximumError=std::max(maximumError,std::max(sumError,squareError));
  pass=checksumCount==scan.Count &&
      sumError<=512.0*std::numeric_limits<double>::epsilon() &&
      squareError<=512.0*std::numeric_limits<double>::epsilon() && pass;

  // Integer multiplicities contain no roulette.  Their exact per-correction
  // histogram proves that each original c was retained by every local clone,
  // MPI pack/unpack, and continuation.  Fractional cases use the independent
  // analytic flux variance because their per-bin populations are stochastic.
  cloneCheckPerformed=false;
  const double roundedMultiplicity=std::round(meanMultiplicity);
  if (roundedMultiplicity>=1.0 &&
      std::fabs(meanMultiplicity-roundedMultiplicity)<=
          64.0*std::numeric_limits<double>::epsilon()*
              std::max(meanMultiplicity,1.0)) {
    const long long multiplier=static_cast<long long>(roundedMultiplicity);
    cloneCheckPerformed=multiplier>1;
    const int unitIndex=UnitCorrectionIndex();
    for (int index=0;index<ValidationCorrectionCount;index++) {
      const long long sourceBin=UseUnitCorrections ?
          (index==unitIndex ? sourceCount : 0) :
          sourceCount/ValidationCorrectionCount+
              (index<sourceCount%ValidationCorrectionCount ? 1 : 0);
      pass=scan.CorrectionCounts[index]==multiplier*sourceBin && pass;
    }
  }
  return pass;
}

bool CheckScanFluxChecksum(const cParticleScan &scan,
    const cDiagnostic &diagnostic,double destinationCarrier,
    long long repeatedFaces,double &maximumError) {
  if (!(destinationCarrier>0.0) || !std::isfinite(destinationCarrier) ||
      repeatedFaces<=0) return false;
  const double representedNumber=repeatedFaces*destinationCarrier*
      scan.CorrectionSum;
  const double representedInternalEnergy=repeatedFaces*destinationCarrier*
      scan.CorrectionWeightedInternalEnergy;
  const double numberScale=std::max({std::fabs(representedNumber),
      std::fabs(diagnostic.OutgoingNumberFlux),1.0});
  const double internalScale=std::max({
      std::fabs(representedInternalEnergy),
      std::fabs(diagnostic.OutgoingInternalEnergyFlux),
      std::numeric_limits<double>::min()});
  const double numberError=
      std::fabs(representedNumber-diagnostic.OutgoingNumberFlux)/numberScale;
  const double internalError=std::fabs(representedInternalEnergy-
      diagnostic.OutgoingInternalEnergyFlux)/internalScale;
  maximumError=std::max(maximumError,std::max(numberError,internalError));

  // The diagnostic and the post-move scan visit the same particles in
  // different linked-list orders.  Their long sums are therefore not
  // bitwise reproducible even though every summand is identical.  A fixed
  // ULP allowance becomes invalid once an integer face split creates tens of
  // thousands of clones.  Bound the comparison by the number of additions
  // in both accumulations.  The factor 16 covers multiplication and MPI
  // reduction rounding while retaining a sub-ppb bound in the largest Stage
  // 4 validation case.  Exact per-correction histograms and per-particle
  // internal-energy checks above still detect a missing or altered clone.
  const double accumulationTolerance=std::max(
      1024.0*std::numeric_limits<double>::epsilon(),
      16.0*std::numeric_limits<double>::epsilon()*
          static_cast<double>(scan.Count+diagnostic.FaceTransferCount));
  return numberError<=accumulationTolerance &&
      internalError<=accumulationTolerance;
}

bool CheckMoment(const char *name,double incoming,double outgoing,
    double variance) {
  const double difference=outgoing-incoming;
  const double scale=std::max(std::fabs(incoming),std::fabs(outgoing));
  if (!(variance>=0.0) || !std::isfinite(variance)) {
    if (PIC::ThisThread==0)
      std::printf("  FAIL %-12s invalid variance %.9e\n",name,variance);
    return false;
  }

  if (variance==0.0) {
    const double tolerance=512.0*std::numeric_limits<double>::epsilon()*
        std::max(scale,std::numeric_limits<double>::min());
    if (std::fabs(difference)>tolerance) {
      if (PIC::ThisThread==0)
        std::printf("  FAIL %-12s exact difference %.9e > %.9e\n",
            name,difference,tolerance);
      return false;
    }
    return true;
  }

  const double sigma=std::sqrt(variance);
  const double z=std::fabs(difference)/sigma;
  const double relative=(std::fabs(incoming)>0.0) ?
      std::fabs(difference/incoming) : 0.0;

  // A fixed relative tolerance is ill-conditioned for a moment whose expected
  // value is zero or small compared with its particle-to-particle variance.
  // Conversely, requiring only a sigma bound can make a well-resolved nonzero
  // flux less precise than the requested engineering tolerance.  Accept the
  // difference when it lies inside the larger of the 2% deterministic band and
  // the five-sigma Monte Carlo band.  This is equivalent to a relative OR a
  // statistical check, and prevents a legitimate low-signal transverse moment
  // from failing solely because its relative denominator is nearly zero.
  const double relativeTolerance=
      RequiredRelativeAccuracy*std::fabs(incoming);
  const double statisticalTolerance=MaximumStandardDeviations*sigma;
  const double allowedDifference=
      std::max(relativeTolerance,statisticalTolerance);
  const bool pass=std::fabs(difference)<=allowedDifference;
  if (!pass && PIC::ThisThread==0)
    std::printf("  FAIL %-12s difference=% .9e allowed=%.9e "
        "z=%.4f relative=%.4g\n",name,difference,allowedDifference,z,
        relative);
  return pass;
}

bool CheckFluxes(const cDiagnostic &diagnostic) {
  bool pass=true;
  pass=CheckMoment("number",diagnostic.IncomingNumberFlux,
      diagnostic.OutgoingNumberFlux,diagnostic.NumberFluxVariance) && pass;
  pass=CheckMoment("mass",diagnostic.IncomingMassFlux,
      diagnostic.OutgoingMassFlux,diagnostic.MassFluxVariance) && pass;
  for (int dimension=0;dimension<3;dimension++) {
    char name[20];
    std::snprintf(name,sizeof(name),"momentum[%d]",dimension);
    pass=CheckMoment(name,diagnostic.IncomingMomentumFlux[dimension],
        diagnostic.OutgoingMomentumFlux[dimension],
        diagnostic.MomentumFluxVariance[dimension]) && pass;
  }
  pass=CheckMoment("kinetic",diagnostic.IncomingKineticEnergyFlux,
      diagnostic.OutgoingKineticEnergyFlux,
      diagnostic.KineticEnergyFluxVariance) && pass;
  pass=CheckMoment("internal",diagnostic.IncomingInternalEnergyFlux,
      diagnostic.OutgoingInternalEnergyFlux,
      diagnostic.InternalEnergyFluxVariance) && pass;
  return pass;
}

bool RunSingleFaceCase(const std::vector<cNode*> &leaves,
    const cFacePair &pair,double timeStepRatio,double weightRatio,int direction,
    long int particleNumber,const char *moverName,bool mpiFace,int caseNumber,
    cValidationSummary &summary) {
  cNode *source=(direction>0) ? pair.Left : pair.Right;
  cNode *destination=(direction>0) ? pair.Right : pair.Left;
  const double face=pair.Left->xmax[0];
  const double width=std::min(source->xmax[0]-source->xmin[0],
      destination->xmax[0]-destination->xmin[0]);
  const auto velocities=BuildVelocities(width,direction);

  PIC::ParticleBuffer::DeleteAllParticles();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  SetFaceTransportParameters(leaves,source,destination,timeStepRatio,
      weightRatio);
  rnd_seed(71000+caseNumber);
  InjectSingleFaceParticles(source,face,particleNumber,velocities);

  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(true);
  double movingTime=0.0,exchangeTime=0.0;
  PIC::TimeStepInternal::ExecutionTrackDefault(movingTime,exchangeTime);
  const cDiagnostic diagnostic=ReduceDiagnostic();

  const double radialVelocity=velocities[0][0];
  const double expectedRadialProbe=face+direction*PIC::Mesh::mesh->EPS+
      radialVelocity*(1.0-CrossingFraction)*timeStepRatio;
  const cParticleScan scan=
      ScanParticles(destination,&velocities,expectedRadialProbe);

  bool pass=true;
  pass=(diagnostic.FaceTransferCount==particleNumber) && pass;
  pass=(scan.Count==diagnostic.DestinationParticleCount) && pass;
  pass=(scan.WrongNode==0 && scan.NonzeroRemainder==0 &&
      scan.InvalidState==0 && scan.InvalidSpeed==0 &&
      scan.InvalidRadialProbe==0 && scan.InvalidCorrection==0 &&
      scan.InvalidInternalEnergy==0) && pass;
  pass=(diagnostic.FinalPendingParticleCount==0) && pass;
  pass=(diagnostic.InitialPendingParticleCount==
      diagnostic.DestinationParticleCount) && pass;
  pass=(diagnostic.FirstPostExchangePendingParticleCount==
      diagnostic.InitialPendingParticleCount) && pass;
  pass=(diagnostic.ContinuationDeletedParticleCount==0) && pass;
  pass=(diagnostic.ContinuationSweepCount==1) && pass;
  pass=(diagnostic.MaximumRemainderInvariantError<=
      128.0*std::numeric_limits<double>::epsilon()) && pass;
  const double meanMultiplicity=timeStepRatio/weightRatio;
  bool cloneCheckPerformed=false;
  const bool correctionPass=CheckCorrectionChecksum(scan,particleNumber,
      meanMultiplicity,summary.MaximumCorrectionError,cloneCheckPerformed);
  const bool scanFluxPass=CheckScanFluxChecksum(scan,diagnostic,
      GetOwnerCarrier(destination),1,summary.MaximumFluxChecksumError);
  pass=correctionPass && scanFluxPass && pass;
  const bool fluxPass=CheckFluxes(diagnostic);
  pass=fluxPass && pass;

  if (mpiFace) summary.MpiFaceChecks++;
  else summary.LocalFaceChecks++;
  if (correctionPass && scanFluxPass)
    summary.CorrectionPreservationChecks++;
  else summary.PartialMutations++;
  if (cloneCheckPerformed && correctionPass)
    summary.CloneCorrectionChecks++;
  if (HasValidationInternalEnergyModes() &&
      scan.InvalidInternalEnergy==0 && scanFluxPass && fluxPass)
    summary.InternalEnergyChecks++;

  if (PIC::ThisThread==0) {
    std::printf("%s mover=%s face=%s direction=%s dt-ratio=%.3g "
        "weight-ratio=%.3g mu=%.3g input=%ld output=%lld sweeps=%d "
        "max-pending=%ld\n",
        pass ? "PASS" : "FAIL",moverName,mpiFace ? "MPI" : "local",
        direction>0 ? "out" : "in",timeStepRatio,weightRatio,
        timeStepRatio/weightRatio,particleNumber,scan.Count,
        diagnostic.ContinuationSweepCount,
        diagnostic.MaximumPendingParticleCount);
    if (!pass)
      std::printf("  transfers=%lld destination-count=%lld wrong-node=%lld "
          "remainders=%lld invalid-state=%lld invalid-speed=%lld "
          "invalid-position=%lld invalid-correction=%lld invalid-internal=%lld "
          "correction-checksum=%s scan-flux-checksum=%s "
          "initial-pending=%ld post-exchange=%ld "
          "moved=%ld deleted=%ld final-pending=%ld flux=%s "
          "remainder-error=%.9e\n",
          diagnostic.FaceTransferCount,diagnostic.DestinationParticleCount,
          scan.WrongNode,scan.NonzeroRemainder,scan.InvalidState,
          scan.InvalidSpeed,scan.InvalidRadialProbe,scan.InvalidCorrection,
          scan.InvalidInternalEnergy,correctionPass ? "pass" : "fail",
          scanFluxPass ? "pass" : "fail",
          diagnostic.InitialPendingParticleCount,
          diagnostic.FirstPostExchangePendingParticleCount,
          diagnostic.ContinuationMovedParticleCount,
          diagnostic.ContinuationDeletedParticleCount,
          diagnostic.FinalPendingParticleCount,fluxPass ? "pass" : "fail",
          diagnostic.MaximumRemainderInvariantError);
  }

  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(false);
  return pass;
}

struct cMultiFaceSpan {
  std::size_t Begin=0,End=0;
  int OwnershipTransitions=0;
  bool IsValid() const {return End>Begin;}
};

cMultiFaceSpan FindMultiFaceSpan(const std::vector<cNode*> &leaves,
    int requiredOwnershipTransitions) {
  // The continuation regression must cross several block faces even in a
  // single-rank run.  Ownership transitions add an MPI requirement; they do
  // not replace the minimum geometric path length.
  constexpr std::size_t MinimumFaces=3;
  cMultiFaceSpan best;
  std::size_t bestLength=std::numeric_limits<std::size_t>::max();
  for (std::size_t begin=1;begin+2<leaves.size();begin++) {
    int transitions=0;
    for (std::size_t end=begin+1;end+1<leaves.size();end++) {
      if (leaves[end-1]->Thread!=leaves[end]->Thread) transitions++;
      if (end-begin<MinimumFaces || transitions<requiredOwnershipTransitions)
        continue;
      if (transitions>=requiredOwnershipTransitions && end-begin<bestLength) {
        best={begin,end,transitions};
        bestLength=end-begin;
        break;
      }
    }
  }
  return best;
}

bool RunMultiFaceCase(const std::vector<cNode*> &leaves,int rankCount,
    cValidationSummary &summary) {
  const int requiredTransitions=(rankCount>=3) ? 2 : (rankCount==2 ? 1 : 0);
  const cMultiFaceSpan span=FindMultiFaceSpan(leaves,requiredTransitions);
  if (!span.IsValid()) {
    if (PIC::ThisThread==0)
      std::printf("FAIL mover=first multi-face: no suitable ownership span\n");
    return false;
  }

  PIC::ParticleBuffer::DeleteAllParticles();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  for (cNode *node:leaves) if (node->block!=NULL) {
    node->block->SetLocalTimeStep(1.0,TestSpecies);
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
    node->block->SetLocalParticleWeight(1.0,TestSpecies);
#endif
  }

  cNode *source=leaves[span.Begin],*destination=leaves[span.End];
  const double sourceWidth=source->xmax[0]-source->xmin[0];
  const double initialRadius=source->xmax[0]-0.10*sourceWidth;
  const double targetRadius=0.5*(destination->xmin[0]+destination->xmax[0]);
  double position[3]={initialRadius,0.0,0.0};
  double velocity[3]={targetRadius-initialRadius,0.0,0.0};
  const long int particleNumber=32;
  if (source->Thread==PIC::ThisThread) {
    for (long int particle=0;particle<particleNumber;particle++) {
      int species=TestSpecies;
      const long int ptr=PIC::ParticleBuffer::InitiateParticle(position,
          velocity,NULL,&species,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,source);
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      InitializeValidationParticleState(particle,data);
      PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,data);
    }
  }

  rnd_seed(99001);
  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(true);
  double movingTime=0.0,exchangeTime=0.0;
  PIC::TimeStepInternal::ExecutionTrackDefault(movingTime,exchangeTime);
  const cDiagnostic diagnostic=ReduceDiagnostic();
  const cParticleScan scan=ScanParticles(destination,NULL,0.0);
  const long long faceNumber=static_cast<long long>(span.End-span.Begin);

  bool cloneCheckPerformed=false;
  const bool correctionPass=CheckCorrectionChecksum(scan,particleNumber,1.0,
      summary.MaximumCorrectionError,cloneCheckPerformed);
  const bool scanFluxPass=CheckScanFluxChecksum(scan,diagnostic,
      GetOwnerCarrier(destination),faceNumber,summary.MaximumFluxChecksumError);
  bool pass=scan.Count==particleNumber && scan.WrongNode==0 &&
      scan.NonzeroRemainder==0 && scan.InvalidState==0 &&
      scan.InvalidCorrection==0 && scan.InvalidInternalEnergy==0 &&
      diagnostic.FaceTransferCount==faceNumber*particleNumber &&
      diagnostic.DestinationParticleCount==faceNumber*particleNumber &&
      diagnostic.ContinuationSweepCount==faceNumber &&
      diagnostic.FinalPendingParticleCount==0 && correctionPass &&
      scanFluxPass && CheckFluxes(diagnostic);

  summary.MultiFaceChecks++;
  if (correctionPass && scanFluxPass)
    summary.CorrectionPreservationChecks++;
  else summary.PartialMutations++;
  if (HasValidationInternalEnergyModes() &&
      scan.InvalidInternalEnergy==0 && scanFluxPass)
    summary.InternalEnergyChecks++;

  if (PIC::ThisThread==0)
    std::printf("%s mover=first multi-face faces=%lld ownership-crossings=%d "
        "particles=%lld sweeps=%d max-pending=%ld correction-checksum=%s "
        "scan-flux-checksum=%s\n",
        pass ? "PASS" : "FAIL",faceNumber,span.OwnershipTransitions,
        scan.Count,diagnostic.ContinuationSweepCount,
        diagnostic.MaximumPendingParticleCount,
        correctionPass ? "pass" : "fail",scanFluxPass ? "pass" : "fail");

  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(false);
  return pass;
}

bool RunMoverSuite(const std::vector<cNode*> &leaves,const cFacePair &localPair,
    const cFacePair &mpiPair,long int particleNumber,const char *name,
    PIC::Mover::fSpeciesDependentParticleMover mover,int &caseNumber,
    cValidationSummary &summary) {
  PIC::Mover::SetLocalTimeStepValidationMover(mover);
  struct cRatioCase {double TimeStep,Weight;};
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  // These cases independently exercise roulette and cloning caused by the
  // timestep ratio, the weight ratio, both ratios reinforcing one another,
  // and exact cancellation when W is proportional to dt.  The expected
  // destination multiplicity is mu=dt_ratio/weight_ratio because the source
  // timestep and weight are both one.
  const cRatioCase cases[]={
      {0.25,1.0},{0.5,2.0},{1.0,0.4},{2.0,0.5},
      {2.5,2.5},{0.5,0.25},{2.0,4.0}};
#else
  const cRatioCase cases[]={
      {0.25,1.0},{0.5,1.0},{1.0,1.0},{2.0,1.0},{2.5,1.0}};
#endif
  bool pass=true;
  for (const cRatioCase &ratio:cases) for (int direction: {-1,1})
    pass=RunSingleFaceCase(leaves,localPair,ratio.TimeStep,ratio.Weight,
        direction,particleNumber,name,false,caseNumber++,summary) && pass;

  if (mpiPair.IsValid()) {
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
    const cRatioCase mpiCases[]={{0.5,2.0},{1.0,0.4},{2.5,2.5}};
#else
    const cRatioCase mpiCases[]={{0.5,1.0},{1.0,1.0},{2.5,1.0}};
#endif
    for (const cRatioCase &ratio:mpiCases) for (int direction: {-1,1})
      pass=RunSingleFaceCase(leaves,mpiPair,ratio.TimeStep,ratio.Weight,
          direction,particleNumber,name,true,caseNumber++,summary) && pass;
  }
  return pass;
}

}

bool Comet::LocalTimeStepValidation::Run(const cOptions &options) {
#if DIM != 1
  if (PIC::ThisThread==0)
    std::fprintf(stderr,"local-timestep validation requires DIM1\n");
  return false;
#elif _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ != _PIC_MODE_ON_
  if (PIC::ThisThread==0)
    std::fprintf(stderr,"local-timestep flux resampling is not compiled in\n");
  return false;
#else
  if (options.ParticlesPerCase<=0) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,"validation particle count must be positive\n");
    return false;
  }

  PIC::DomainBlockDecomposition::UpdateBlockTable();
  const std::vector<cNode*> leaves=GetOrderedLeaves();
  const cFacePair localPair=FindFacePair(leaves,false);
  const cFacePair mpiPair=FindFacePair(leaves,true);
  if (!localPair.IsValid()) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,"no equal-width local block pair for validation\n");
    return false;
  }
  if (PIC::nTotalThreads>1 && !mpiPair.IsValid()) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,"no equal-width MPI ownership face for validation\n");
    return false;
  }

  UseUnitCorrections=options.UseUnitCorrections;
  const double reportedCorrectionMinimum=UseUnitCorrections ? 1.0 :
      ValidationCorrectionCycle[0];
  const double reportedCorrectionMaximum=UseUnitCorrections ? 1.0 :
      ValidationCorrectionCycle[ValidationCorrectionCount-1];
  const int reportedCorrectionCount=UseUnitCorrections ? 1 :
      ValidationCorrectionCount;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  const char *correctionPattern=UseUnitCorrections ? "unit" : "wide";
#else
  // A correction-off build has no per-particle field and the compile-time
  // cycle contains only the identity value, regardless of the runtime option.
  const char *correctionPattern="unit";
#endif

  if (PIC::ThisThread==0)
    std::printf("Local-timestep validation: particles/case=%ld ranks=%d "
        "species=%s correction-pattern=%s correction-bins=%d "
        "correction-range=%.5g:%.5g "
        "deterministic-internal-energy=%s relative-limit=%.3g z-limit=%.1f\n",
        options.ParticlesPerCase,PIC::nTotalThreads,
        PIC::MolecularData::GetChemSymbol(TestSpecies),
        correctionPattern,reportedCorrectionCount,reportedCorrectionMinimum,
        reportedCorrectionMaximum,
        HasDeterministicValidationInternalEnergy() ? "yes" : "no",
        RequiredRelativeAccuracy,MaximumStandardDeviations);

  bool pass=true;
  int caseNumber=0;
  cValidationSummary summary;
  if (options.TestFirstOrder)
    pass=RunMoverSuite(leaves,localPair,mpiPair,options.ParticlesPerCase,
        "first",PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory,
        caseNumber,summary) && pass;

  if (options.TestSecondOrder)
    pass=RunMoverSuite(leaves,localPair,mpiPair,options.ParticlesPerCase,
        "second",
        PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder,
        caseNumber,summary) && pass;

  // The production Halley case uses the first-order tracer.  Its long path is
  // the required check of repeated exchange/continuation sweeps and, with at
  // least three ranks, two distinct ownership transitions.
  if (options.TestFirstOrder && options.TestMultiFace) {
    PIC::Mover::SetLocalTimeStepValidationMover(
        PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory);
    pass=RunMultiFaceCase(leaves,PIC::nTotalThreads,summary) && pass;
  }

  PIC::ParticleBuffer::DeleteAllParticles();
  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(false);
  PIC::Mover::SetLocalTimeStepValidationMover(NULL);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  int localPass=pass ? 1 : 0,globalPass=0;
  MPI_Allreduce(&localPass,&globalPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    const char *correctionMode="ON";
#else
    const char *correctionMode="OFF";
#endif
    std::printf("%s individual-weight local-timestep validation "
        "correction-mode=%s correction-pattern=%s mpi-ranks=%d "
        "correction-min=%.17g "
        "correction-max=%.17g correction-values=%d "
        "max-correction-error=%.17e max-flux-checksum-error=%.17e "
        "partial-mutations=%lld "
        "local-face-checks=%lld mpi-face-checks=%lld multiface-checks=%lld "
        "correction-preservation-checks=%lld clone-correction-checks=%lld "
        "internal-energy-checks=%lld\n",
        globalPass ? "PASS" : "FAIL",correctionMode,correctionPattern,
        PIC::nTotalThreads,reportedCorrectionMinimum,
        reportedCorrectionMaximum,reportedCorrectionCount,
        summary.MaximumCorrectionError,summary.MaximumFluxChecksumError,
        summary.PartialMutations,summary.LocalFaceChecks,
        summary.MpiFaceChecks,summary.MultiFaceChecks,
        summary.CorrectionPreservationChecks,summary.CloneCorrectionChecks,
        summary.InternalEnergyChecks);
    // Preserve the historical focused-test sentinel for correction-off
    // regression tooling while the richer Stage-4 record above carries the
    // new compatibility evidence.
    std::printf("Local-timestep validation result: %s\n",
        globalPass ? "PASS" : "FAIL");
  }
  UseUnitCorrections=false;
  return globalPass!=0;
#endif
}

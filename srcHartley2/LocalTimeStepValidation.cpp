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
constexpr int TestSpecies=0;

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
  double MaximumRemainderInvariantError=0.0;
  int ContinuationSweepCount=0;
  long int MaximumPendingParticleCount=0,FinalPendingParticleCount=0;
  long int InitialPendingParticleCount=0;
  long int FirstPostExchangePendingParticleCount=0;
  long int ContinuationMovedParticleCount=0;
  long int ContinuationDeletedParticleCount=0;
};

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

void SetTimeSteps(const std::vector<cNode*> &leaves,cNode *source,
    cNode *destination,double destinationRatio) {
  // Boundary-layer copies are present only on ranks that need them.  Setting
  // every allocated copy from global node identity keeps the source and
  // destination values identical on both sides of an MPI ownership face.
  for (cNode *node:leaves) if (node->block!=NULL) {
    double timeStep=1.0;
    if (node==source) timeStep=1.0;
    if (node==destination) timeStep=destinationRatio;
    node->block->SetLocalTimeStep(timeStep,TestSpecies);
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
    PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,
        PIC::ParticleBuffer::GetParticleDataPointer(ptr));
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

  const double localValues[18]={
      local.IncomingNumberFlux,local.OutgoingNumberFlux,
      local.NumberFluxVariance,local.IncomingMassFlux,
      local.OutgoingMassFlux,local.MassFluxVariance,
      local.IncomingMomentumFlux[0],local.IncomingMomentumFlux[1],
      local.IncomingMomentumFlux[2],local.OutgoingMomentumFlux[0],
      local.OutgoingMomentumFlux[1],local.OutgoingMomentumFlux[2],
      local.MomentumFluxVariance[0],local.MomentumFluxVariance[1],
      local.MomentumFluxVariance[2],local.IncomingKineticEnergyFlux,
      local.OutgoingKineticEnergyFlux,local.KineticEnergyFluxVariance};
  double globalValues[18]={};
  MPI_Allreduce(localValues,globalValues,18,MPI_DOUBLE,MPI_SUM,
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
  long long localValues[6]={local.Count,local.WrongNode,
      local.NonzeroRemainder,local.InvalidState,local.InvalidSpeed,
      local.InvalidRadialProbe};
  long long globalValues[6]={};
  MPI_Allreduce(localValues,globalValues,6,MPI_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  global.Count=globalValues[0];
  global.WrongNode=globalValues[1];
  global.NonzeroRemainder=globalValues[2];
  global.InvalidState=globalValues[3];
  global.InvalidSpeed=globalValues[4];
  global.InvalidRadialProbe=globalValues[5];
  return global;
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
  return pass;
}

bool RunSingleFaceCase(const std::vector<cNode*> &leaves,
    const cFacePair &pair,double ratio,int direction,long int particleNumber,
    const char *moverName,bool mpiFace,int caseNumber) {
  cNode *source=(direction>0) ? pair.Left : pair.Right;
  cNode *destination=(direction>0) ? pair.Right : pair.Left;
  const double face=pair.Left->xmax[0];
  const double width=std::min(source->xmax[0]-source->xmin[0],
      destination->xmax[0]-destination->xmin[0]);
  const auto velocities=BuildVelocities(width,direction);

  PIC::ParticleBuffer::DeleteAllParticles();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  SetTimeSteps(leaves,source,destination,ratio);
  rnd_seed(71000+caseNumber);
  InjectSingleFaceParticles(source,face,particleNumber,velocities);

  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(true);
  double movingTime=0.0,exchangeTime=0.0;
  PIC::TimeStepInternal::ExecutionTrackDefault(movingTime,exchangeTime);
  const cDiagnostic diagnostic=ReduceDiagnostic();

  const double radialVelocity=velocities[0][0];
  const double expectedRadialProbe=face+direction*PIC::Mesh::mesh->EPS+
      radialVelocity*(1.0-CrossingFraction)*ratio;
  const cParticleScan scan=
      ScanParticles(destination,&velocities,expectedRadialProbe);

  bool pass=true;
  pass=(diagnostic.FaceTransferCount==particleNumber) && pass;
  pass=(scan.Count==diagnostic.DestinationParticleCount) && pass;
  pass=(scan.WrongNode==0 && scan.NonzeroRemainder==0 &&
      scan.InvalidState==0 && scan.InvalidSpeed==0 &&
      scan.InvalidRadialProbe==0) && pass;
  pass=(diagnostic.FinalPendingParticleCount==0) && pass;
  pass=(diagnostic.InitialPendingParticleCount==
      diagnostic.DestinationParticleCount) && pass;
  pass=(diagnostic.FirstPostExchangePendingParticleCount==
      diagnostic.InitialPendingParticleCount) && pass;
  pass=(diagnostic.ContinuationDeletedParticleCount==0) && pass;
  pass=(diagnostic.ContinuationSweepCount==1) && pass;
  pass=(diagnostic.MaximumRemainderInvariantError<=
      128.0*std::numeric_limits<double>::epsilon()) && pass;
  const bool fluxPass=CheckFluxes(diagnostic);
  pass=fluxPass && pass;

  if (PIC::ThisThread==0) {
    std::printf("%s mover=%s face=%s direction=%s dt-ratio=%.3g "
        "input=%ld output=%lld sweeps=%d max-pending=%ld\n",
        pass ? "PASS" : "FAIL",moverName,mpiFace ? "MPI" : "local",
        direction>0 ? "out" : "in",ratio,particleNumber,scan.Count,
        diagnostic.ContinuationSweepCount,
        diagnostic.MaximumPendingParticleCount);
    if (!pass)
      std::printf("  transfers=%lld destination-count=%lld wrong-node=%lld "
          "remainders=%lld invalid-state=%lld invalid-speed=%lld "
          "invalid-position=%lld initial-pending=%ld post-exchange=%ld "
          "moved=%ld deleted=%ld final-pending=%ld flux=%s "
          "remainder-error=%.9e\n",
          diagnostic.FaceTransferCount,diagnostic.DestinationParticleCount,
          scan.WrongNode,scan.NonzeroRemainder,scan.InvalidState,
          scan.InvalidSpeed,scan.InvalidRadialProbe,
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

bool RunMultiFaceCase(const std::vector<cNode*> &leaves,int rankCount) {
  const int requiredTransitions=(rankCount>=3) ? 2 : (rankCount==2 ? 1 : 0);
  const cMultiFaceSpan span=FindMultiFaceSpan(leaves,requiredTransitions);
  if (!span.IsValid()) {
    if (PIC::ThisThread==0)
      std::printf("FAIL mover=first multi-face: no suitable ownership span\n");
    return false;
  }

  PIC::ParticleBuffer::DeleteAllParticles();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  for (cNode *node:leaves) if (node->block!=NULL)
    node->block->SetLocalTimeStep(1.0,TestSpecies);

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
      PIC::Mover::SetLocalTimeStepRemainderFraction(0.0,
          PIC::ParticleBuffer::GetParticleDataPointer(ptr));
    }
  }

  rnd_seed(99001);
  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(true);
  double movingTime=0.0,exchangeTime=0.0;
  PIC::TimeStepInternal::ExecutionTrackDefault(movingTime,exchangeTime);
  const cDiagnostic diagnostic=ReduceDiagnostic();
  const cParticleScan scan=ScanParticles(destination,NULL,0.0);
  const long long faceNumber=static_cast<long long>(span.End-span.Begin);

  bool pass=scan.Count==particleNumber && scan.WrongNode==0 &&
      scan.NonzeroRemainder==0 && scan.InvalidState==0 &&
      diagnostic.FaceTransferCount==faceNumber*particleNumber &&
      diagnostic.DestinationParticleCount==faceNumber*particleNumber &&
      diagnostic.ContinuationSweepCount==faceNumber &&
      diagnostic.FinalPendingParticleCount==0 && CheckFluxes(diagnostic);

  if (PIC::ThisThread==0)
    std::printf("%s mover=first multi-face faces=%lld ownership-crossings=%d "
        "particles=%lld sweeps=%d max-pending=%ld\n",
        pass ? "PASS" : "FAIL",faceNumber,span.OwnershipTransitions,
        scan.Count,diagnostic.ContinuationSweepCount,
        diagnostic.MaximumPendingParticleCount);

  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(false);
  return pass;
}

bool RunMoverSuite(const std::vector<cNode*> &leaves,const cFacePair &localPair,
    const cFacePair &mpiPair,long int particleNumber,const char *name,
    PIC::Mover::fSpeciesDependentParticleMover mover,int &caseNumber) {
  PIC::Mover::SetLocalTimeStepValidationMover(mover);
  const double ratios[]={0.25,0.5,1.0,2.0,2.5};
  bool pass=true;
  for (double ratio:ratios) for (int direction: {-1,1})
    pass=RunSingleFaceCase(leaves,localPair,ratio,direction,particleNumber,
        name,false,caseNumber++) && pass;

  if (mpiPair.IsValid()) {
    const double mpiRatios[]={0.5,1.0,2.5};
    for (double ratio:mpiRatios) for (int direction: {-1,1})
      pass=RunSingleFaceCase(leaves,mpiPair,ratio,direction,particleNumber,
          name,true,caseNumber++) && pass;
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

  if (PIC::ThisThread==0)
    std::printf("Local-timestep validation: particles/case=%ld ranks=%d "
        "relative-limit=%.3g z-limit=%.1f\n",options.ParticlesPerCase,
        PIC::nTotalThreads,RequiredRelativeAccuracy,MaximumStandardDeviations);

  bool pass=true;
  int caseNumber=0;
  if (options.TestFirstOrder)
    pass=RunMoverSuite(leaves,localPair,mpiPair,options.ParticlesPerCase,
        "first",PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory,
        caseNumber) && pass;

  if (options.TestSecondOrder)
    pass=RunMoverSuite(leaves,localPair,mpiPair,options.ParticlesPerCase,
        "second",
        PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder,
        caseNumber) && pass;

  // The production Halley case uses the first-order tracer.  Its long path is
  // the required check of repeated exchange/continuation sweeps and, with at
  // least three ranks, two distinct ownership transitions.
  if (options.TestFirstOrder && options.TestMultiFace) {
    PIC::Mover::SetLocalTimeStepValidationMover(
        PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory);
    pass=RunMultiFaceCase(leaves,PIC::nTotalThreads) && pass;
  }

  PIC::ParticleBuffer::DeleteAllParticles();
  PIC::Mover::ResetLocalTimeStepFluxDiagnostic(false);
  PIC::Mover::SetLocalTimeStepValidationMover(NULL);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  int localPass=pass ? 1 : 0,globalPass=0;
  MPI_Allreduce(&localPass,&globalPass,1,MPI_INT,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0)
    std::printf("Local-timestep validation result: %s\n",
        globalPass ? "PASS" : "FAIL");
  return globalPass!=0;
#endif
}

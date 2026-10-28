#include "RadialSimulation.h"

#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

#include "pic.h"
#include "CometModel.h"

bool Comet::RadialSimulation::Enabled=false;
double Comet::RadialSimulation::CourantNumber=-1.0;
double Comet::RadialSimulation::CharacteristicSpeed[PIC::nTotalSpecies]={0.0};
double Comet::RadialSimulation::InnerBoundaryTemperature=-1.0;
long int Comet::RadialSimulation::MaxIterations=-1;
long int Comet::RadialSimulation::SamplingLength=-1;
double Comet::RadialSimulation::ConvergenceTolerance=-1.0;
int Comet::RadialSimulation::RequiredConvergedWindows=-1;

namespace {
bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
}

bool HasInnerFace(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  bool externalFaces[2*DIM];
  if (PIC::Mesh::mesh->ExternalBoundaryBlock(node,externalFaces)!=
      _EXTERNAL_BOUNDARY_BLOCK_) return false;
  return externalFaces[0];
}

cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *FindInnerBoundaryNode(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return NULL;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_)
    return node->block!=NULL && HasInnerFace(node) ? node : NULL;

  for (int child=0;child<(1<<DIM);child++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *result=
        FindInnerBoundaryNode(node->downNode[child]);
    if (result!=NULL) return result;
  }
  return NULL;
}

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
void SetWeightProportionalToTimeStep(int species,double surfaceWeight,
    double surfaceTimeStep,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->block!=NULL) {
      // W/dt is held constant with radius.  Because dt is proportional to the
      // local cell width, a steady spherical flux then uses approximately the
      // same number of model particles in every logarithmic radial cell.
      const double localTimeStep=node->block->GetLocalTimeStep(species);
      node->block->SetLocalParticleWeight(
          surfaceWeight*localTimeStep/surfaceTimeStep,species);
    }
    return;
  }

  for (int child=0;child<(1<<DIM);child++)
    SetWeightProportionalToTimeStep(species,surfaceWeight,surfaceTimeStep,
        node->downNode[child]);
}
#endif

long int CountProfileCells(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return 0;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_)
    return _BLOCK_CELLS_X_;

  long int result=0;
  for (int child=0;child<(1<<DIM);child++)
    result+=CountProfileCells(node->downNode[child]);
  return result;
}

void FillLocalIntegratedDensityProfile(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,long int &cellIndex,
    std::vector<double> &profile) {
  if (node==NULL) return;

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    // Every MPI rank traverses the complete AMR tree in the same child order,
    // so cellIndex is a stable global radial-cell index even when ownership is
    // redistributed at a sampling break. Only the current owner contributes;
    // MPI_Allreduce below assembles the complete profile on every rank.
    if (node->Thread==PIC::ThisThread && node->block!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) {
        PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(
            _getCenterNodeLocalNumber(i,0,0));
        if (cell==NULL || !(cell->Measure>0.0) ||
            !std::isfinite(cell->Measure)) {
          std::fprintf(stderr,
              "invalid cell measure while evaluating radial convergence\n");
          std::abort();
        }

        for (int species=0;species<PIC::nTotalSpecies;species++) {
          const double density=cell->GetNumberDensity(species);
          if (!(density>=0.0) || !std::isfinite(density)) {
            std::fprintf(stderr,
                "invalid number density for species %d while evaluating radial convergence\n",
                species);
            std::abort();
          }

          // Storing n*V makes the subsequent discrete L1 norm exactly the
          // shell-volume-weighted norm. Measure is the core's 4*pi/3 times
          // (r_outer^3-r_inner^3) value in DIM1 spherical symmetry.
          profile[(cellIndex+i)*PIC::nTotalSpecies+species]=
              density*cell->Measure;
        }
      }
    }

    cellIndex+=_BLOCK_CELLS_X_;
    return;
  }

  for (int child=0;child<(1<<DIM);child++)
    FillLocalIntegratedDensityProfile(node->downNode[child],cellIndex,profile);
}

std::vector<double> CollectIntegratedDensityProfile() {
  const long int cellCount=CountProfileCells(PIC::Mesh::mesh->rootTree);
  if (cellCount<=0) {
    std::fprintf(stderr,"radial convergence found no mesh cells\n");
    std::abort();
  }

  const size_t valueCount=
      static_cast<size_t>(cellCount)*PIC::nTotalSpecies;
  if (valueCount>static_cast<size_t>(INT_MAX)) {
    std::fprintf(stderr,
        "radial convergence profile exceeds the MPI count limit\n");
    std::abort();
  }

  std::vector<double> localProfile(valueCount,0.0);
  std::vector<double> globalProfile(valueCount,0.0);
  long int visitedCells=0;
  FillLocalIntegratedDensityProfile(PIC::Mesh::mesh->rootTree,visitedCells,
      localProfile);
  if (visitedCells!=cellCount) {
    std::fprintf(stderr,"inconsistent radial convergence cell traversal\n");
    std::abort();
  }

  MPI_Allreduce(localProfile.data(),globalProfile.data(),
      static_cast<int>(valueCount),MPI_DOUBLE,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  return globalProfile;
}

double EvaluateRelativeL1Change(const std::vector<double> &current,
    const std::vector<double> &previous,int &limitingSpecies) {
  if (current.size()!=previous.size() ||
      current.size()%PIC::nTotalSpecies!=0) {
    std::fprintf(stderr,"inconsistent radial convergence profile sizes\n");
    std::abort();
  }

  double maximumChange=-1.0;
  limitingSpecies=0;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    double differenceNorm=0.0,currentNorm=0.0,previousNorm=0.0;
    for (size_t offset=species;offset<current.size();
         offset+=PIC::nTotalSpecies) {
      differenceNorm+=std::abs(current[offset]-previous[offset]);
      currentNorm+=std::abs(current[offset]);
      previousNorm+=std::abs(previous[offset]);
    }

    // max(||n_k||_1,||n_{k-1}||_1) gives a symmetric reference scale and
    // avoids division by zero. A species absent from both windows has zero
    // change; appearance or disappearance produces a relative change of one.
    const double referenceNorm=std::max(currentNorm,previousNorm);
    const double relativeChange=
        referenceNorm>0.0 ? differenceNorm/referenceNorm : 0.0;
    if (relativeChange>maximumChange) {
      maximumChange=relativeChange;
      limitingSpecies=species;
    }
  }

  return maximumChange;
}
}

bool Comet::RadialSimulation::Validate(std::string *error) {
  if (!Enabled) return true;
  if (!(CourantNumber>0.0) || CourantNumber>1.0 ||
      !std::isfinite(CourantNumber))
    return Fail("radial Courant number must be in (0,1]",error);
  if (!(InnerBoundaryTemperature>0.0) ||
      !std::isfinite(InnerBoundaryTemperature))
    return Fail("radial inner-boundary temperature must be positive",error);
  if (MaxIterations<=0)
    return Fail("radial maximum iteration count must be positive",error);
  if (SamplingLength<=0)
    return Fail("radial sampling length must be positive",error);
  if (SamplingLength>MaxIterations)
    return Fail("radial sampling length cannot exceed the maximum iteration count",error);
  if (!(ConvergenceTolerance>0.0) || ConvergenceTolerance>1.0 ||
      !std::isfinite(ConvergenceTolerance))
    return Fail("radial convergence tolerance must be in (0,1]",error);
  if (RequiredConvergedWindows<=0)
    return Fail("required converged-window count must be positive",error);

  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!(CharacteristicSpeed[species]>0.0) ||
        !std::isfinite(CharacteristicSpeed[species])) {
      char message[160];
      std::snprintf(message,sizeof(message),
          "radial characteristic speed is missing for species %d",species);
      return Fail(message,error);
    }
  }
  return true;
}

double Comet::RadialSimulation::LocalTimeStep(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  // Restrict one move to the selected fraction of the local radial cell size.
  // Each species has its own explicitly supplied speed so light, energetic
  // photoproducts cannot inherit the much slower parent-gas timestep.
  return CourantNumber*node->GetCharacteristicCellSize()/
      CharacteristicSpeed[species];
}

double Comet::RadialSimulation::LocalBlockInjectionRate(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  // A uniform spherical source has one global rate Q. In DIM1 the inner box
  // face represents the entire 4*pi*R^2 shell, so Q is assigned once to the
  // block touching face 0; no area division is needed here.
  return HasInnerFace(node) ? Comet::Model::SurfaceProductionRate[species] : 0.0;
}

bool Comet::RadialSimulation::InjectionBlockIndicator(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (!HasInnerFace(node)) return false;
  for (int species=0;species<PIC::nTotalSpecies;species++)
    if (Comet::Model::SurfaceProductionRate[species]>0.0) return true;
  return false;
}

long int Comet::RadialSimulation::InjectBlock(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (!HasInnerFace(node)) return 0;

  long int injected=0;
  const double externalNormal[3]={-1.0,0.0,0.0};
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    const double sourceRate=Comet::Model::SurfaceProductionRate[species];
    if (!(sourceRate>0.0)) continue;

    const double particleWeight=node->block->GetLocalParticleWeight(species);
    const double timeStep=node->block->GetLocalTimeStep(species);
    const double modelRate=sourceRate/particleWeight;
    double eventTime=0.0;
    const double bulkVelocity[3]={
        Comet::Model::SurfaceBulkSpeed[species],0.0,0.0};

    // A Poisson process avoids systematic rounding when Q*dt/weight is not an
    // integer. The remaining fraction of the current step is passed to the
    // mover so injected particles are uniformly distributed in injection time.
    while ((eventTime+=-std::log(rnd())/modelRate)<timeStep) {
      const long int ptr=PIC::ParticleBuffer::GetNewParticle();
      PIC::ParticleBuffer::byte *data=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      double position[3]={node->xmin[0]+PIC::Mesh::mesh->EPS,0.0,0.0};
      double velocity[3];
      PIC::Distribution::InjectMaxwellianDistribution(velocity,bulkVelocity,
          Comet::Model::SurfaceTemperature[species],
          const_cast<double*>(externalNormal),species,-1);

      PIC::ParticleBuffer::SetX(position,data);
      PIC::ParticleBuffer::SetV(velocity,data);
      PIC::ParticleBuffer::SetI(species,data);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,data);
#endif
      PIC::Mover::DispatchParticleMover(ptr,timeStep-eventTime,node);
      injected++;
    }
  }
  return injected;
}

int Comet::RadialSimulation::ProcessDomainBoundary(long int ptr,
    double *position,double *velocity,int intersectionFace,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  (void)position;
  if (intersectionFace==1) {
    // Rubin et al. remove every particle crossing the outer boundary.
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }
  if (intersectionFace!=0) {
    std::fprintf(stderr,"unexpected DIM1 comet boundary face %d\n",
        intersectionFace);
    std::abort();
  }

  // Gas returning to the top of the Knudsen layer is re-emitted. The paper
  // states that the gas is completely thermalized there, so discard the
  // incident velocity and draw a new outward half-Maxwellian at the explicit
  // inner-boundary temperature. Dust is outside this radial gas backend.
  const int species=PIC::ParticleBuffer::GetI(ptr);
  const double externalNormal[3]={-1.0,0.0,0.0};
  const double bulkVelocity[3]={0.0,0.0,0.0};
  PIC::Distribution::InjectMaxwellianDistribution(velocity,bulkVelocity,
      InnerBoundaryTemperature,const_cast<double*>(externalNormal),species,-1);
  (void)node;
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

void Comet::RadialSimulation::Initialize() {
  // Sampling storage and cell blocks must exist before AMPS initializes
  // weights, timesteps, collision data, or its particle mover.
  PIC::Mesh::initCellSamplingDataBuffer();
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();

  PIC::Mover::ProcessOutsideDomainParticles=ProcessDomainBoundary;
  PIC::Init_AfterParser();
  PIC::Mover::Init();

  // Use the case-controlled window length for AMPS macroscopic sampling.
  // Convergence comparisons are made only after a complete window; a partial
  // final window at the hard iteration ceiling is never declared converged.
  PIC::RequiredSampleLength=SamplingLength;

  PIC::ParticleWeightTimeStep::LocalTimeStep=LocalTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();
  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=
      LocalBlockInjectionRate;

  // Direct parent sources determine their statistical weights. Daughter
  // species have no injection rate, so their base weight is tied to the
  // smallest initialized parent weight; event processing applies an exact
  // individual correction whenever a reaction changes species.
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *innerNode=
      FindInnerBoundaryNode(PIC::Mesh::mesh->rootTree);

  // After domain decomposition, only the rank that owns the block adjacent to
  // the nucleus has that block allocated. Confirm that there is exactly one
  // owner globally, then broadcast surface values with an all-reduction. The
  // previous per-rank NULL check aborted every valid rank that did not own the
  // nucleus-adjacent block.
  const int localInnerOwner=(innerNode!=NULL) ? 1 : 0;
  int globalInnerOwners=0;
  MPI_Allreduce(&localInnerOwner,&globalInnerOwners,1,MPI_INT,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  if (globalInnerOwners!=1) {
    if (PIC::ThisThread==0)
      std::fprintf(stderr,
          "radial comet mesh has %d allocated inner-boundary blocks; expected exactly one\n",
          globalInnerOwners);
    std::abort();
  }

  // Start each direct parent with AMPS' source-derived constant weight.  Use
  // its value at the surface as the normalization, then scale W with the
  // local timestep.  Global weights would make the number of particles per
  // cell grow in proportion to the logarithmic cell width and render the
  // 10^6 km Halley domain impractical at sub-metre surface resolution.
  double daughterWeight=std::numeric_limits<double>::max();
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (Comet::Model::SurfaceProductionRate[species]>0.0) {
      PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(species);

      // Only the owning rank reads the local block. Summation acts as a
      // broadcast because the ownership check above guarantees one source.
      const double localSurfaceValues[2]={
          innerNode!=NULL ? innerNode->block->GetLocalParticleWeight(species) : 0.0,
          innerNode!=NULL ? innerNode->block->GetLocalTimeStep(species) : 0.0};
      double surfaceValues[2]={0.0,0.0};
      MPI_Allreduce(localSurfaceValues,surfaceValues,2,MPI_DOUBLE,MPI_SUM,
          MPI_GLOBAL_COMMUNICATOR);
      const double surfaceWeight=surfaceValues[0];
      const double surfaceTimeStep=surfaceValues[1];
      if (!(surfaceWeight>0.0) || !std::isfinite(surfaceWeight) ||
          !(surfaceTimeStep>0.0) || !std::isfinite(surfaceTimeStep)) {
        if (PIC::ThisThread==0)
          std::fprintf(stderr,
              "invalid surface weight or timestep for radial species %d\n",
              species);
        std::abort();
      }

      SetWeightProportionalToTimeStep(species,surfaceWeight,surfaceTimeStep,
          PIC::Mesh::mesh->rootTree);
      if (surfaceWeight<daughterWeight) daughterWeight=surfaceWeight;
    }
  }
  if (!(daughterWeight<std::numeric_limits<double>::max())) {
    std::fprintf(stderr,"radial comet run has no direct surface source\n");
    std::abort();
  }

  // Give each daughter the same surface base weight and scale it with that
  // daughter's timestep.  A reaction's individual statistical correction
  // still preserves the exact physical weight inherited from its parent.
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!(Comet::Model::SurfaceProductionRate[species]>0.0)) {
      const double localSurfaceTimeStep=
          innerNode!=NULL ? innerNode->block->GetLocalTimeStep(species) : 0.0;
      double surfaceTimeStep=0.0;
      MPI_Allreduce(&localSurfaceTimeStep,&surfaceTimeStep,1,MPI_DOUBLE,MPI_SUM,
          MPI_GLOBAL_COMMUNICATOR);
      if (!(surfaceTimeStep>0.0) || !std::isfinite(surfaceTimeStep)) {
        if (PIC::ThisThread==0)
          std::fprintf(stderr,
              "invalid surface timestep for radial daughter species %d\n",
              species);
        std::abort();
      }

      SetWeightProportionalToTimeStep(species,daughterWeight,surfaceTimeStep,
          PIC::Mesh::mesh->rootTree);
    }
  }
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  double daughterWeight=std::numeric_limits<double>::max();
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (Comet::Model::SurfaceProductionRate[species]>0.0) {
      PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(species);
      if (PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]<daughterWeight)
        daughterWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[species];
    }
  }
  if (!(daughterWeight<std::numeric_limits<double>::max())) {
    std::fprintf(stderr,"radial comet run has no direct surface source\n");
    std::abort();
  }
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!(PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]>0.0))
      PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(
          species,daughterWeight);
  }
#else
#error "The radial comet backend requires species-dependent particle weights"
#endif

  PIC::BC::BlockInjectionBCindicatior=InjectionBlockIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=InjectBlock;
  PIC::BC::InitBoundingBoxInjectionBlockList();
}

void Comet::RadialSimulation::Run() {
  // Compare only complete, non-overlapping AMPS sampling windows. The first
  // window establishes the reference profile; each later window contributes
  // one pass or failure to the consecutive-window convergence counter.
  std::vector<double> previousProfile;
  int convergedWindows=0;
  long int completedWindows=0;
  long int previousOutputNumber=PIC::DataOutputFileNumber;

  for (long int iteration=0;iteration<MaxIterations;iteration++) {
    PIC::TimeStep();

    if (PIC::DataOutputFileNumber==previousOutputNumber) continue;
    previousOutputNumber=PIC::DataOutputFileNumber;

    // Alarm-triggered partial output must never satisfy a convergence test
    // configured for a full SamplingLength window.
    if (PIC::LastSampleLength!=SamplingLength) {
      if (PIC::ThisThread==0) {
        std::fprintf(PIC::DiagnospticMessageStream,
            "$PREFIX: radial convergence skipped partial window of %ld iterations\n",
            PIC::LastSampleLength);
        std::fflush(PIC::DiagnospticMessageStream);
      }
      continue;
    }

    completedWindows++;
    std::vector<double> currentProfile=CollectIntegratedDensityProfile();
    if (previousProfile.empty()) {
      previousProfile.swap(currentProfile);
      if (PIC::ThisThread==0) {
        std::fprintf(PIC::DiagnospticMessageStream,
            "$PREFIX: radial convergence reference window=%ld iteration=%ld\n",
            completedWindows,iteration+1);
        std::fflush(PIC::DiagnospticMessageStream);
      }
      continue;
    }

    int limitingSpecies=0;
    const double maximumChange=EvaluateRelativeL1Change(currentProfile,
        previousProfile,limitingSpecies);
    previousProfile.swap(currentProfile);
    if (maximumChange<=ConvergenceTolerance)
      convergedWindows++;
    else
      convergedWindows=0;

    if (PIC::ThisThread==0) {
      char symbol[_MAX_STRING_LENGTH_PIC_];
      PIC::MolecularData::GetChemSymbol(symbol,limitingSpecies);
      std::fprintf(PIC::DiagnospticMessageStream,
          "$PREFIX: radial convergence window=%ld iteration=%ld max-relative-L1=%e limiting-species=%s consecutive=%d/%d tolerance=%e\n",
          completedWindows,iteration+1,maximumChange,symbol,convergedWindows,
          RequiredConvergedWindows,ConvergenceTolerance);
      std::fflush(PIC::DiagnospticMessageStream);
    }

    if (convergedWindows>=RequiredConvergedWindows) {
      if (PIC::ThisThread==0) {
        std::fprintf(PIC::DiagnospticMessageStream,
            "$PREFIX: radial simulation converged after %ld iterations and %ld complete sampling windows\n",
            iteration+1,completedWindows);
        std::fflush(PIC::DiagnospticMessageStream);
      }
      return;
    }
  }

  if (PIC::ThisThread==0) {
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: radial simulation reached the hard limit of %ld iterations without satisfying convergence (%d/%d consecutive windows)\n",
        MaxIterations,convergedWindows,RequiredConvergedWindows);
    std::fflush(PIC::DiagnospticMessageStream);
  }
}

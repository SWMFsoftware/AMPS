// One physical trial of the approved I24 total-capture case (2026-11-25)
// through the normal Moon startup, force tracer, particle allocator, reaction
// controller, deferred reinsertion and MPI exchange. The separate runner
// reconstructs the expected event clock, shell, weights and budgets.
#include "total_capture_setup.h"
#include "MoonSodiumTotalCaptureCohort.h"
#include "MoonInput.h"
#include "MoonTransportClock.h"
#include "Moon.h"
#include "SpiceUsr.h"
#include "../../sodium_cohort/zero_force_recipe.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

void amps_init();
#define I24_STRING_INNER(...) #__VA_ARGS__
#define I24_STRING(...) I24_STRING_INNER(__VA_ARGS__)
namespace S=Moon::Species::SodiumTotalCapture;
namespace C=Moon::Chemistry::SodiumTotalCaptureCohort;
namespace M=Moon::Chemistry::SemanticRandom;
namespace T=PIC::Mover::ForceCutSurface;
namespace fs=std::filesystem;
namespace {
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
const std::array<double,3> weights={{7,56,.875}};
constexpr double represented=2.625;

void InstallTables(double duration,double scale){
  const double ratio=Moon::Runtime::GetConfiguration().hydrogenToSodiumReferenceWeightRatio;
  // Read the H weight installed by normal startup before replacing the entire
  // cohort's weight table. This catches an ignored .in key or an uninitialised
  // linked neutral without attributing a test override to production.
  for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){
    auto* n=PIC::DomainBlockDecomposition::BlockTable[b];
    if (!n || !n->block) continue;
    const double sodium=n->block->GetLocalParticleWeight(_NA_SPEC_);
    const double hydrogen=n->block->GetLocalParticleWeight(_H_SPEC_);
    if (!(sodium>0) || hydrogen!=ratio*sodium || !(n->block->GetLocalTimeStep(_H_SPEC_)>0))
      throw std::runtime_error("normal startup did not initialize the configured hydrogen timestep and weight");
  }
  // Global native timestep setters reduce to their minimum. Reinitialize the
  // reduction sentinel exactly as in the dedicated I21 cohort before the first
  // step. Local-mode tables instead accept direct assignments per block.
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  for(int s=0;s<PIC::nTotalSpecies;++s)PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1;
#elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  PIC::ParticleWeightTimeStep::GlobalTimeStep[0]=-1;
#elif _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  throw std::runtime_error("unsupported PIC timestep mode for prescribed I24 step");
#endif
  for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){
    auto* n=PIC::DomainBlockDecomposition::BlockTable[b];
    if(!n||!n->block)continue;
    for(int s=0;s<PIC::nTotalSpecies;++s){
      n->block->SetLocalTimeStep(duration,s);
      const double referenceWeight=scale*weights.at(s);
      n->block->SetLocalParticleWeight(referenceWeight,s);
      if(n->block->GetLocalTimeStep(s)!=duration||n->block->GetLocalParticleWeight(s)!=referenceWeight)
        throw std::runtime_error("readback disagrees with prescribed local dt or model-particle weight");
    }
  }
}

void InjectParent(double epoch){
  double x[]={1e7,1001,1003},v[]={0,0,0};
  Node* node=PIC::Mesh::mesh->findTreeNode(x);
  if(!node)throw std::runtime_error("I24 physical source falls outside the model mesh");
  if(node->Thread!=PIC::ThisThread)return;
  if(PIC::ParticleBuffer::GetMaxNPart()-PIC::ParticleBuffer::GetAllPartNum()<3)
    throw std::runtime_error("one parent and its two reserved products require three available particle slots");
  int i,j,k; PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false);
  long& head=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
  const long p=PIC::ParticleBuffer::GetNewParticle();
  auto* data=PIC::ParticleBuffer::GetParticleDataPointer(p);
  PIC::ParticleBuffer::SetI(0,data);PIC::ParticleBuffer::SetX(x,data);PIC::ParticleBuffer::SetV(v,data);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(represented/node->block->GetLocalParticleWeight(_NA_SPEC_),data);
  T::ResetParticleState(data);Moon::Chemistry::Controller::ResetForNewParticle(data);T::SetParticleEpochS(p,epoch);
  std::string error;
  if(!M::InitializeSource(p,0,&error))throw std::runtime_error(error);
  PIC::ParticleBuffer::SetPrev(-1,data);PIC::ParticleBuffer::SetNext(head,data);
  if(head!=-1)PIC::ParticleBuffer::SetPrev(p,head);
  head=p;
}

std::ofstream Open(const fs::path& dir,const std::string& prefix,const char* header){
  auto name=dir/(prefix+"-rank-"+std::to_string(PIC::ThisThread)+".csv");
  if(fs::exists(name))throw std::runtime_error("refusing to replace native receipt "+name.string());
  std::ofstream stream(name);stream.exceptions(std::ios::badbit|std::ios::failbit);
  stream<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<header<<'\n';
  return stream;
}
template<class V> void Vector(std::ostream& out,const V& values){
  for(int j=0;j<3;++j)out<<','<<values[j];
}

void Record(const fs::path& out,const C::Prepared& owner,double origin){
  if(PIC::ThisThread==0){
    // SPICE exposes the actual loaded paths. The external runner resolves
    // and hashes each one for the exact epoch and LSO frame used in this run.
    std::ofstream paths(out/"spice-kernels.txt");
    paths.exceptions(std::ios::badbit|std::ios::failbit);
    SpiceInt count=0;ktotal_c("ALL",&count);
    if(count==0)throw std::runtime_error("no actual SPICE kernels loaded");
    for(SpiceInt i=0;i<count;++i){
      SpiceChar file[4096],kind[64],source[4096];SpiceInt handle;SpiceBoolean found;
      kdata_c(i,"ALL",sizeof(file),sizeof(kind),sizeof(source),file,kind,source,&handle,&found);
      if(!found || failed_c())throw std::runtime_error("could not enumerate actual loaded kernels");
      paths<<fs::absolute(file).string()<<'\n';
    }
    paths.close();
  }
  auto particles=Open(out,"particles","parent,child,species,x,y,z,vx,vy,vz,reference_weight,correction,represented,epoch_et");
  auto events=Open(out,"events","parent,elapsed_s,epoch_et,origin_et,rate_s,shell,shell_uniform,radiated_J,Q_J,bath_number,bath_mass_kg,bath_px,bath_py,bath_pz,bath_energy_J,tracked_px,tracked_py,tracked_pz,tracked_energy_J,tracked_charge_e,bath_charge_e,integer_ok,kinematic_ok");
  auto products=Open(out,"products","parent,ordinal,species,x,y,z,vx,vy,vz,reference_weight,correction,represented,epoch_et,remaining_s");
  if(owner.receipts().size()!=owner.ledger().events.size())
    throw std::runtime_error("prospective sodium total captures disagree with atomic PIC publications");
  // Read the fresh product bytes stored by ParticleTransaction at publication.
  // Final per-particle rows below come from live native mesh lists after the
  // allocated packets have moved for the remainder of the original step.
  for(std::size_t i=0;i<owner.ledger().events.size();++i){
    const auto& event=owner.ledger().events[i];const auto& receipt=owner.receipts()[i];
    const auto& l=event.ledger;
    if(event.identity!=S::ModelId||event.physicalEpochS!=receipt.epochEt||
       l.bathTerms.size()!=1||!l.bathTerms[0].reactant||l.bathTerms[0].species!=1||
       event.initialProducts.size()!=2)throw std::runtime_error("source, epoch, bath or product identity mismatch");
    events<<receipt.key.originalParent<<','<<receipt.interval.elapsedS()<<','<<event.physicalEpochS<<','<<origin<<','
          <<receipt.composition.evaluation.rateS<<','<<receipt.composition.shell<<','<<receipt.shellUniform<<','
          <<receipt.composition.radiatedEnergyJ<<','<<receipt.composition.shellQJ<<','
          <<l.bathTerms[0].representedNumber<<','<<l.bathModelMassChangeKg;
    Vector(events,l.bathMomentumChangeKgMS);events<<','<<l.bathEnergyChangeJ;
    Vector(events,l.trackedMomentumChangeKgMS);
    events<<','<<l.trackedEnergyChangeJ<<','<<l.trackedChargeChangeElementary<<','<<l.bathChargeChangeElementary
          <<','<<l.integerConservationVerified<<','<<l.kinematicConservationVerified<<'\n';
    for(std::size_t j=0;j<event.initialProducts.size();++j){
      const auto& p=event.initialProducts[j];
      products<<receipt.key.originalParent<<','<<j<<','<<p.compiledSpecies;
      Vector(products,p.positionM);Vector(products,p.velocityMS);
      products<<','<<p.referenceWeight<<','<<p.individualCorrection<<','<<p.referenceWeight*p.individualCorrection
              <<','<<p.physicalEpochS<<','<<p.interval.remainingS()<<'\n';
    }
  }
  for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){
    auto* n=PIC::DomainBlockDecomposition::BlockTable[b];if(!n||!n->block)continue;
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell){
      for(long ptr=n->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)){
        M::State key;std::string error;if(!M::Read(ptr,&key,&error))throw std::runtime_error(error);
        if(key.genealogy.size()>1)throw std::runtime_error("unexpected multigeneration capture genealogy");
        double x[3],v[3];PIC::ParticleBuffer::GetX(x,ptr);PIC::ParticleBuffer::GetV(v,ptr);
        const auto s=PIC::ParticleBuffer::GetI(ptr);const double weight=n->block->GetLocalParticleWeight(s);
        const double correction=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
        particles<<key.originalParent<<','<<(key.genealogy.empty()?-1:static_cast<long>(key.genealogy.front()))
                 <<','<<static_cast<int>(s);Vector(particles,x);Vector(particles,v);
        particles<<','<<weight<<','<<correction<<','<<weight*correction<<','<<T::GetParticleEpochS(ptr)<<'\n';
      }
    }
  }
  particles.close();events.close();products.close();
}
}

// argv: 1-10 sources (total_capture_setup.h); 11 Moon .in; 12 output;
// 13 horizon [s]; 14 beam energy [eV/u]; 15-16 energies below and above the
// approved envelope [eV/u] for the native refusal controls; optional
// 17 steps and 18 reference-weight scale.
int main(int argc,char** argv){
  if(argc!=17 && argc!=19)return 2;
  try{
    PIC::InitMPI();
    Moon::Runtime::Configuration config;std::string error;
    if(!Moon::Runtime::ParseApplicationInput(argv[11],&config,&error)||
       !Moon::Runtime::InstallConfiguration(config,&error))throw std::runtime_error(error);
    S::Definition definition;C::Request request;
    I24TotalSetup::Load(17,argv,&definition,&request.provider,_NA_SPEC_,_NA_PLUS_SPEC_,_H_SPEC_,Moon::Frames::SolarOrbital);
    const auto beam=[&](long double energyEvPerU){
      const long double speed=std::sqrt(2*energyEvPerU*definition.electronVoltJ/definition.atomicMassUnitKg);
      return S::Beam{2.5e6L,{{.6L*speed,.8L*speed,0}}};};
    request.beam=beam(std::stold(argv[14]));
    request.random.seed=731993;request.random.cohort="I24-PIC-total-capture-v1";request.random.genealogyCapacity=2;
    request.maximumRootRefinements=256;request.elementaryChargeC=static_cast<double>(definition.electronVoltJ);
    request.elementaryChargeIdentity="2019 SI exact elementary charge";
    // Native refusal controls on the production owner, before any particle
    // buffer exists. Each invalid request must be refused without leaving a
    // registered owner, so the valid registration below must still succeed.
    const auto refused=[](const C::Request& trial){
      std::shared_ptr<C::Prepared> candidate;std::string reason;
      return !C::Prepared::RegisterForNormalStartup(trial,&candidate,&reason)&&!candidate&&!reason.empty();
    };
    S::Definition otherFrame;auto trial=request;
    I24TotalSetup::Load(17,argv,&otherFrame,&trial.provider,_NA_SPEC_,_NA_PLUS_SPEC_,_H_SPEC_,"source-target-rest-inertial");
    const bool frameRefused=refused(trial);
    trial=request;trial.beam=beam(std::stold(argv[15]));const bool lowRefused=refused(trial);
    trial=request;trial.beam=beam(std::stold(argv[16]));const bool highRefused=refused(trial);
    trial=request;trial.beam.protonDensityM3=-1;const bool densityRefused=refused(trial);
    trial=request;trial.maximumRootRefinements=0;const bool effortRefused=refused(trial);
    trial=request;trial.elementaryChargeC=1.6e-19;const bool chargeRefused=refused(trial);
    if(PIC::ThisThread==0){
      const fs::path controlFile=fs::path(argv[12])/"registration-controls.csv";
      std::ofstream controls(controlFile);
      controls.exceptions(std::ios::badbit|std::ios::failbit);
      controls<<"non_lso_frame_refused,below_envelope_refused,above_envelope_refused,"
                "negative_density_refused,zero_root_effort_refused,wrong_elementary_charge_refused\n"
              <<frameRefused<<','<<lowRefused<<','<<highRefused<<','<<densityRefused<<','
              <<effortRefused<<','<<chargeRefused<<'\n';
    }
    std::shared_ptr<C::Prepared> owner;
    if(!C::Prepared::RegisterForNormalStartup(request,&owner,&error))throw std::runtime_error(error);
    amps_init();PIC::DomainBlockDecomposition::UpdateBlockTable();
    if(PIC::nTotalSpecies!=3||!Moon::CommonTransportTimeSelected()||!Moon::GetCommonTransportClock())
      throw std::runtime_error("capture requires three generated transport species and the common physical clock");
    if(T::Acceleration||!Moon::SodiumCohortInput::ExactZeroAccelerationMacro(
        I24_STRING(_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(a,s,p,x,v,n))))
      throw std::runtime_error("I24 controlled first collision requires actual zero-force native acceleration");
    const double horizon=std::stod(argv[13]);
    const int steps=argc==19?std::stoi(argv[17]):1;
    const double scale=argc==19?std::stod(argv[18]):1.;
    if(!std::isfinite(horizon)||horizon<=0||(steps!=1&&steps!=2&&steps!=4)||
        (scale!=1.&&scale!=.5&&scale!=.25))throw std::runtime_error("invalid I24 horizon or refinement level");
    // Powers of two are exact in binary64. Across refinements the original
    // represented population and final physical horizon are held fixed.
    InstallTables(horizon/steps,scale);
    const double epoch=Moon::GetCommonTransportClock()->originEtS();
    InjectParent(epoch);
    auto history=Open(argv[12],"steps","step,relative_end_s,committed_events,neutral_parents,sodium_ions,hydrogen,represented_Na,represented_Na_ion,represented_H");
    for(int step=1;step<=steps;++step) {
      if(!Moon::BeginCommonTransportStep(epoch,&error))throw std::runtime_error(error);
      Moon::PrepareProductionParticleEnvironment();
      double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
      if(!Moon::CommitCommonTransportStep(&error))throw std::runtime_error(error);
      // Read owned PIC lists after native reinsertion and MPI exchange.
      std::array<long,3> count{};
      std::array<double,3> population{};
      for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){
        const auto* node=PIC::DomainBlockDecomposition::BlockTable[b];
        if(!node||!node->block)continue;
        for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
          for(long p=node->block->FirstCellParticleTable[cell];p!=-1;p=PIC::ParticleBuffer::GetNext(p)){
            const auto s=PIC::ParticleBuffer::GetI(p);
            ++count.at(s);
            population.at(s)+=node->block->GetLocalParticleWeight(s)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
          }
      }
      history<<step<<','<<Moon::GetCommonTransportClock()->elapsedS()<<','<<owner->ledger().events.size();
      Vector(history,count);Vector(history,population);history<<'\n';
    }
    history.close();
    Record(argv[12],*owner,epoch);
    MPI_Finalize();return 0;
  }catch(const std::exception& e){
    std::cerr<<"I24 PIC total-capture error: "<<e.what()<<'\n';
    int initialized=0;MPI_Initialized(&initialized);
    if(initialized)MPI_Abort(MPI_GLOBAL_COMMUNICATOR,2);
    return 2;
  }
}

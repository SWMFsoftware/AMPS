// One physical I23 trial through the normal Moon startup (including the new
// He/He+ timestep/weight keys), force tracer, particle allocator, reaction
// controller, deferred reinsertion and MPI exchange. The separate runner
// reconstructs the expected event clock, weights and signed electron bath.
#include "helium_impact_setup.h"
#include "MoonHeliumImpactCohort.h"
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
#define I23_STRING_INNER(...) #__VA_ARGS__
#define I23_STRING(...) I23_STRING_INNER(__VA_ARGS__)
namespace C=Moon::Chemistry::HeliumImpactCohort;
namespace M=Moon::Chemistry::SemanticRandom;
namespace T=PIC::Mover::ForceCutSurface;
namespace fs=std::filesystem;
using namespace I23Setup;
namespace {
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;

double Weight(const Arguments& a,int species){
  if(species==_NA_SPEC_)return static_cast<double>(a.Real("weight_na"));
  if(species==_HE_SPEC_)return static_cast<double>(a.Real("weight_he"));
  if(species==_HE_PLUS_SPEC_)return static_cast<double>(a.Real("weight_he_plus"));
  throw std::runtime_error("unexpected compiled species in the I23 build");
}

void InstallTables(const Arguments& a,double duration,double scale){
  const auto& configuration=Moon::Runtime::GetConfiguration();
  // Read the He/He+ weights and He+ timestep installed by NORMAL startup from
  // the new .in keys before replacing the whole cohort table. This detects an
  // ignored key or an uninitialised reaction-born species in production.
  for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){
    auto* n=PIC::DomainBlockDecomposition::BlockTable[b];
    if(!n||!n->block)continue;
    const double sodium=n->block->GetLocalParticleWeight(_NA_SPEC_);
    const double helium=n->block->GetLocalParticleWeight(_HE_SPEC_);
    const double ion=n->block->GetLocalParticleWeight(_HE_PLUS_SPEC_);
    const double ionStep=n->block->GetLocalTimeStep(_HE_PLUS_SPEC_);
    if(!(sodium>0)||helium!=configuration.heliumToSodiumReferenceWeightRatio*sodium||
       ion!=configuration.heliumIonToSodiumReferenceWeightRatio*sodium||!(ionStep>0)||!std::isfinite(ionStep))
      throw std::runtime_error("normal startup did not initialize the configured He/He+ weights and He+ timestep");
  }
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  for(int s=0;s<PIC::nTotalSpecies;++s)PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1;
#elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  PIC::ParticleWeightTimeStep::GlobalTimeStep[0]=-1;
#elif _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  throw std::runtime_error("unsupported PIC timestep mode for the prescribed I23 step");
#endif
  // The common transport clock requires one identical constant step for all
  // species and blocks; the exact binary horizon division is installed here.
  for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){
    auto* n=PIC::DomainBlockDecomposition::BlockTable[b];
    if(!n||!n->block)continue;
    for(int s=0;s<PIC::nTotalSpecies;++s){
      n->block->SetLocalTimeStep(duration,s);
      const double reference=scale*Weight(a,s);
      n->block->SetLocalParticleWeight(reference,s);
      if(n->block->GetLocalTimeStep(s)!=duration||n->block->GetLocalParticleWeight(s)!=reference)
        throw std::runtime_error("readback disagrees with prescribed local dt or model-particle weight");
    }
  }
}

void InjectParent(const Arguments& a,double epoch){
  double x[]={static_cast<double>(a.Real("x_m")),static_cast<double>(a.Real("y_m")),static_cast<double>(a.Real("z_m"))},v[]={0,0,0};
  const double represented=static_cast<double>(a.Real("represented"));
  Node* node=PIC::Mesh::mesh->findTreeNode(x);
  if(!node)throw std::runtime_error("I23 parent position falls outside the model mesh");
  if(node->Thread!=PIC::ThisThread)return;
  if(PIC::ParticleBuffer::GetMaxNPart()-PIC::ParticleBuffer::GetAllPartNum()<2)
    throw std::runtime_error("one parent and its reserved He+ product require two available particle slots");
  int i,j,k;PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false);
  long& head=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
  const long p=PIC::ParticleBuffer::GetNewParticle();
  auto* data=PIC::ParticleBuffer::GetParticleDataPointer(p);
  PIC::ParticleBuffer::SetI(_HE_SPEC_,data);PIC::ParticleBuffer::SetX(x,data);PIC::ParticleBuffer::SetV(v,data);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(represented/node->block->GetLocalParticleWeight(_HE_SPEC_),data);
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
template<class V> void Vector(std::ostream& out,const V& values){for(int j=0;j<3;++j)out<<','<<values[j];}

void Record(const fs::path& out,const C::Prepared& owner,double origin){
  if(PIC::ThisThread==0){
    std::ofstream paths(out/"spice-kernels.txt");
    paths.exceptions(std::ios::badbit|std::ios::failbit);
    SpiceInt count=0;ktotal_c("ALL",&count);
    if(count==0)throw std::runtime_error("no actual SPICE kernels loaded");
    for(SpiceInt i=0;i<count;++i){
      SpiceChar file[4096],kind[64],source[4096];SpiceInt handle;SpiceBoolean found;
      kdata_c(i,"ALL",sizeof(file),sizeof(kind),sizeof(source),file,kind,source,&handle,&found);
      if(!found||failed_c())throw std::runtime_error("could not enumerate actual loaded kernels");
      paths<<fs::absolute(file).string()<<'\n';
    }
    paths.close();
  }
  auto particles=Open(out,"particles","parent,child,species,x,y,z,vx,vy,vz,reference_weight,correction,represented,epoch_et");
  auto events=Open(out,"events","parent,elapsed_s,epoch_et,origin_et,rate_s,Q_J,bath_reactant_number,bath_product_number,"
      "bath_mass_kg,bath_px,bath_py,bath_pz,bath_energy_J,tracked_px,tracked_py,tracked_pz,tracked_energy_J,"
      "tracked_charge_e,bath_charge_e,integer_ok,kinematic_ok,in_kinetic_J,out_kinetic_J,bath_number_change,"
      "bath_charge_change_C,heavy_velocity_bound_m_s");
  auto products=Open(out,"products","parent,ordinal,species,x,y,z,vx,vy,vz,reference_weight,correction,represented,epoch_et,remaining_s");
  if(owner.receipts().size()!=owner.ledger().events.size())
    throw std::runtime_error("prospective helium events disagree with atomic PIC publications");
  for(std::size_t i=0;i<owner.ledger().events.size();++i){
    const auto& event=owner.ledger().events[i];const auto& receipt=owner.receipts()[i];
    const auto& l=event.ledger;
    if(event.physicalEpochS!=receipt.epochEt||l.bathTerms.size()!=2||event.initialProducts.size()!=1)
      throw std::runtime_error("epoch, bath-term or product layout mismatch at a committed I23 event");
    double reactant=0,product=0;
    for(const auto& term:l.bathTerms)(term.reactant?reactant:product)+=term.representedNumber;
    const auto& rr=receipt.composition.reservoir;
    events<<receipt.key.originalParent<<','<<receipt.interval.elapsedS()<<','<<event.physicalEpochS<<','<<origin<<','
          <<receipt.rateS<<','<<rr.stateQJ<<','<<reactant<<','<<product<<','<<l.bathModelMassChangeKg;
    Vector(events,l.bathMomentumChangeKgMS);events<<','<<l.bathEnergyChangeJ;
    Vector(events,l.trackedMomentumChangeKgMS);
    events<<','<<l.trackedEnergyChangeJ<<','<<l.trackedChargeChangeElementary<<','<<l.bathChargeChangeElementary
          <<','<<l.integerConservationVerified<<','<<l.kinematicConservationVerified<<','<<rr.incomingBath.kineticEnergyJ
          <<','<<rr.outgoingBath.kineticEnergyJ<<','<<rr.bathNumberChange<<','<<rr.bathChargeChangeC<<','
          <<rr.error.heavyVelocityUpperMS<<'\n';
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
        if(key.genealogy.size()>1)throw std::runtime_error("unexpected multigeneration helium genealogy");
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

int main(int argc,char** argv){
  try{
    PIC::InitMPI();
    const auto a=Parse(argc,argv);
    Moon::Runtime::Configuration config;std::string error;
    if(!Moon::Runtime::ParseApplicationInput(a["input"],&config,&error)||
       !Moon::Runtime::InstallConfiguration(config,&error))throw std::runtime_error(error);
    P::Registry profile;LoadProfile(a,&profile);
    C::Request request;
    Require(Prepare(a,profile,_HE_SPEC_,_HE_PLUS_SPEC_,&request.provider));
    request.bath=Bath(a,a.Real("density_m3"),Moon::Frames::SolarOrbital);
    request.admission=Admission(a,profile,Moon::Frames::SolarOrbital);
    request.random.seed=std::stoull(a["seed"]);request.random.cohort=a["cohort"];request.random.genealogyCapacity=2;
    request.maximumRootRefinements=256;
    P::Constant charge;
    if(!profile.GetConstant("elementary charge","C",&charge,&error))throw std::runtime_error(error);
    request.elementaryChargeC=static_cast<double>(charge.value);request.elementaryChargeIdentity="CODATA-2022 exact SI elementary charge (state profile)";
    // Native refusal controls on the production owner, before any particle
    // buffer exists. Each invalid request must be refused without leaving a
    // registered owner, so the valid registration below must still succeed.
    const auto refused=[](const C::Request& trial){
      std::shared_ptr<C::Prepared> candidate;std::string reason;
      return !C::Prepared::RegisterForNormalStartup(trial,&candidate,&reason)&&!candidate&&!reason.empty();
    };
    auto trial=request;trial.bath.frameIdentity="J2000";const bool frameRefused=refused(trial);
    trial=request;trial.bath.densityM3=1e7;trial.admission.crosswalkPath.clear();const bool admissionRefused=refused(trial);
    trial=request;trial.bath.densityM3=1e7;trial.admission.evidence.maximumHeavyVelocityErrorMS=1;const bool toleranceRefused=refused(trial);
    trial=request;trial.bath.densityM3=-1;const bool densityRefused=refused(trial);
    trial=request;trial.bath.thermalEnergyEv=11;const bool temperatureRefused=refused(trial);
    trial=request;trial.maximumRootRefinements=0;const bool effortRefused=refused(trial);
    if(PIC::ThisThread==0){
      const fs::path controlFile=fs::path(a["output"])/"registration-controls.csv";
      std::ofstream controls(controlFile);
      controls.exceptions(std::ios::badbit|std::ios::failbit);
      controls<<"non_lso_frame_refused,missing_admission_refused,tolerance_below_bound_refused,"
                "negative_density_refused,unapproved_temperature_refused,zero_root_effort_refused\n"
              <<frameRefused<<','<<admissionRefused<<','<<toleranceRefused<<','<<densityRefused<<','
              <<temperatureRefused<<','<<effortRefused<<'\n';
    }
    std::shared_ptr<C::Prepared> owner;
    if(!C::Prepared::RegisterForNormalStartup(request,&owner,&error))throw std::runtime_error(error);
    amps_init();PIC::DomainBlockDecomposition::UpdateBlockTable();
    if(PIC::nTotalSpecies!=3||!Moon::CommonTransportTimeSelected()||!Moon::GetCommonTransportClock())
      throw std::runtime_error("I23 requires three generated transport species and the common physical clock");
    if(T::Acceleration||!Moon::SodiumCohortInput::ExactZeroAccelerationMacro(
        I23_STRING(_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(a,s,p,x,v,n))))
      throw std::runtime_error("I23 stationary-target bath requires the actual zero-force native acceleration");
    const double horizon=static_cast<double>(a.Real("horizon_s"));
    const int steps=std::stoi(a["steps"]);
    const double scale=static_cast<double>(a.Real("weight_scale"));
    if(!std::isfinite(horizon)||horizon<=0||(steps!=1&&steps!=2&&steps!=4)||(scale!=1.&&scale!=.5&&scale!=.25))
      throw std::runtime_error("invalid I23 horizon or refinement level");
    InstallTables(a,horizon/steps,scale);
    const double epoch=Moon::GetCommonTransportClock()->originEtS();
    InjectParent(a,epoch);
    const fs::path out=a["output"];
    auto history=Open(out,"steps","step,relative_end_s,committed_events,helium,helium_ions,represented_He,represented_He_ion");
    for(int step=1;step<=steps;++step){
      if(!Moon::BeginCommonTransportStep(epoch,&error))throw std::runtime_error(error);
      Moon::PrepareProductionParticleEnvironment();
      double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
      if(!Moon::CommitCommonTransportStep(&error))throw std::runtime_error(error);
      long helium=0,ions=0;double representedHe=0,representedIon=0;
      for(unsigned b=0;b<PIC::DomainBlockDecomposition::nLocalBlocks;++b){
        const auto* node=PIC::DomainBlockDecomposition::BlockTable[b];
        if(!node||!node->block)continue;
        for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
          for(long p=node->block->FirstCellParticleTable[cell];p!=-1;p=PIC::ParticleBuffer::GetNext(p)){
            const auto s=PIC::ParticleBuffer::GetI(p);
            const double represented=node->block->GetLocalParticleWeight(s)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
            if(s==_HE_SPEC_){++helium;representedHe+=represented;}
            else if(s==_HE_PLUS_SPEC_){++ions;representedIon+=represented;}
            else throw std::runtime_error("unexpected species population in the I23 trial");
          }
      }
      history<<step<<','<<Moon::GetCommonTransportClock()->elapsedS()<<','<<owner->ledger().events.size()
             <<','<<helium<<','<<ions<<','<<representedHe<<','<<representedIon<<'\n';
    }
    history.close();
    Record(out,*owner,epoch);
    MPI_Finalize();return 0;
  }catch(const std::exception& e){
    std::cerr<<"I23 PIC helium error: "<<e.what()<<'\n';
    int initialized=0;MPI_Initialized(&initialized);
    if(initialized)MPI_Abort(MPI_GLOBAL_COMMUNICATOR,2);
    return 2;
  }
}

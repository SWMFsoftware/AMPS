// I24 total-capture MPI ensemble adapter. MPI partitions original parent IDs;
// rank plays no role in the physical draws. The adapter feeds authenticated
// inputs to the production provider and writes its outputs. Python scores
// first-hit times, the selected shells and both product velocities.
#include "total_capture_setup.h"
#include "srcMoon/species/QualificationRandom.h"
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace S=Moon::Species::SodiumTotalCapture;
namespace R=Moon::Species::ReactionExecution;
namespace G=Moon::Species::QualificationRandom;

namespace {
using I24TotalSetup::Require;
double Draw(const G::Key& key,G::Purpose purpose,unsigned index) {
  std::string error;
  double value=0;
  Require(G::OpenUniform(key,purpose,index,&value,&error),error);
  return value;
}
R::OpticalDepthClock Clock(double uniform) {
  R::OpticalDepthClock clock;
  std::string error;
  Require(R::OpticalDepthClock::FromOpenUniform(uniform,&clock,&error),error);
  return clock;
}
S::Beam Beam(long double density,long double speed) {
  S::Beam beam;beam.protonDensityM3=density;
  // A non-axis-aligned incident direction exercises the pair-frame geometry.
  beam.protonVelocityMS={{0.6L*speed,0.8L*speed,0}};
  return beam;
}

// argv: 1-10 sources (total_capture_setup.h); 11 density [m^-3]; 12 speed
// [m/s]; 13 horizon [s]; 14 seed; 15 cohort; 16 parents; 17 output;
// 18 parent weight; 19-20 Na+ and H reference weights; 21-22 speeds just
// below and above the approved envelope [m/s].
int Run(int argc,char** argv,int rank,int size) {
  if (argc!=23) throw std::runtime_error("expected sources, beam, cohort and refusal arguments");
  S::Definition definition;
  S::Provider provider;
  I24TotalSetup::Load(argc,argv,&definition,&provider);
  const auto beam=Beam(std::stold(argv[11]),std::stold(argv[12]));
  R::ParentState parent;
  parent.species=0;parent.referenceWeight=std::stod(argv[18]);parent.individualCorrection=1;
  const std::array<double,2> weights={{std::stod(argv[19]),std::stod(argv[20])}};
  S::Evaluation rate;
  Require(provider.Evaluate(beam,parent.velocityMS,&rate));
  const auto output=std::string(argv[17]);
  if (!rank) {
    // Each refusal must leave the previously published evaluation unchanged.
    auto refused=[&](const S::Beam& trial,const R::Vector3& velocity){
      auto prior=rate;return !provider.Evaluate(trial,velocity,&prior).ok() && prior.rateS==rate.rateS;};
    const bool low=refused(Beam(beam.protonDensityM3,std::stold(argv[21])),parent.velocityMS);
    const bool high=refused(Beam(beam.protonDensityM3,std::stold(argv[22])),parent.velocityMS);
    auto trial=beam;trial.protonDensityM3=-1;
    const bool negative=refused(trial,parent.velocityMS);
    trial.protonDensityM3=std::numeric_limits<long double>::quiet_NaN();
    const bool nonfinite=refused(trial,parent.velocityMS);
    trial.protonDensityM3=0;
    S::Evaluation zero=rate;
    const bool zeroRate=provider.Evaluate(trial,parent.velocityMS,&zero).ok() && zero.rateS==0;
    std::ofstream receipt(output+"/rate.csv");
    receipt<<std::setprecision(21)
           <<"rate_s,sigma_m2,speed_m_s,energy_ev_per_u,f1,f2,f3,excitation2_J,excitation3_J,"
             "low_energy_refused,high_energy_refused,negative_density_refused,nonfinite_density_refused,zero_density_zero_rate\n"
           <<rate.rateS<<','<<rate.totalCrossSectionM2<<','<<rate.relativeSpeedMS<<','<<rate.energyEvPerU;
    for (auto f:rate.shellFraction) receipt<<','<<f;
    receipt<<','<<provider.shellExcitationJ(2)<<','<<provider.shellExcitationJ(3)
           <<','<<low<<','<<high<<','<<negative<<','<<nonfinite<<','<<zeroRate<<'\n';
    receipt.close();Require(!receipt.fail(),"cannot write rate receipt");
  }
  const double horizon=std::stod(argv[13]);
  G::Key key;
  key.campaignSeed=std::stoull(argv[14]);key.cohort=argv[15];
  const auto parents=std::stoull(argv[16]);
  R::TraceInterval initial;
  std::string error;
  Require(R::TraceInterval::Create(0,horizon,0,&initial,&error),error);
  std::ofstream native(output+"/parents-rank-"+std::to_string(rank)+".csv");
  native<<std::setprecision(21)
        <<"parent,fate,time_s,shell,nax,nay,naz,hx,hy,hz,na_weight,h_weight,bath_weight,identity_ok,"
          "Q_J,radiated_J,bath_energy_change_J,impulse_bound_kg_m_s,na_recoil_bound_J,photon_recoil_m_s\n";
  const std::vector<R::Hazard> hazards={{0,static_cast<double>(rate.rateS)}};
  bool testedGuards=false;
  for (key.originalParent=rank;key.originalParent<parents;key.originalParent+=size) {
    const auto origin=Clock(Draw(key,G::Purpose::SourceClock,0));
    R::HazardAdvance advance;
    Require(R::AdvanceConstantHazards(hazards,initial,origin,Draw(key,G::Purpose::Channel,0),&advance,&error),error);
    native<<key.originalParent<<',';
    if (!advance.event) {
      native<<"survived,"<<horizon<<",,,,,,,,,,,,,,,,,\n";
      continue;
    }
    const std::array<R::OpticalDepthClock,2> fresh={{Clock(Draw(key,G::Purpose::ProductClock,0)),
                                                     Clock(Draw(key,G::Purpose::ProductClock,1))}};
    S::Composition event;
    const auto uniform=Draw(key,G::Purpose::Kinematics,0);
    Require(provider.Compose(beam,parent,advance.interval,uniform,weights,fresh,&event));
    const auto& products=event.plan.products();
    const auto& ledger=event.plan.ledger();
    Require(products.size()==2 && ledger.bathTerms.size()==1,"unreadable product/bath layout");
    const auto& bath=ledger.bathTerms.front();
    if (rank==0 && !testedGuards) {
      // Reuse the first accepted event for the endpoint and reweighting
      // controls. No parent draw or provider state is consumed.
      auto saved=event;
      bool endpoint=true;
      for (const long double edge:{0.0L,1.0L})
        endpoint=endpoint && !provider.Compose(beam,parent,advance.interval,edge,weights,fresh,&saved).ok() &&
                 saved.plan.valid() && saved.shell==event.shell &&
                 saved.plan.products()[0].velocityMS==products[0].velocityMS;
      auto scaled=weights;
      for (auto& w:scaled) w*=2;
      Require(provider.Compose(beam,parent,advance.interval,uniform,scaled,fresh,&saved));
      bool invariant=saved.plan.products().size()==2;
      for (std::size_t j=0; j<saved.plan.products().size() && invariant; ++j) {
        const auto& packet=saved.plan.products()[j];
        invariant=packet.referenceWeight==scaled[j] &&
            packet.individualCorrection==products[j].individualCorrection/2 &&
            packet.representedNumber==products[j].representedNumber && packet.velocityMS==products[j].velocityMS;
      }
      std::ofstream guards(output+"/event-controls.csv");
      guards<<"closed_endpoints_refused_prior_plan_preserved,rescaled_daughter_weights_invariant\n"
            <<endpoint<<','<<invariant<<'\n';
      guards.close();Require(!guards.fail(),"cannot write provider controls");
      testedGuards=true;
    }
    // Identities as committed to the immutable plan.
    const bool identities=event.plan.valid() && event.plan.scope()==R::PlanScope::Transport &&
       products[0].compiledSpecies==definition.sodiumIonSpecies &&
       products[1].compiledSpecies==definition.hydrogenSpecies &&
       products[1].stateIdentity=="1s" &&
       bath.species==1 && bath.reactant && ledger.integerConservationVerified && ledger.kinematicConservationVerified &&
       ledger.reactantChargeElementary==1 && ledger.productChargeElementary==1 &&
       ledger.trackedChargeChangeElementary==1 && ledger.bathChargeChangeElementary==-1;
    native<<"capture,"<<advance.eventOffsetS<<','<<event.shell;
    for (const auto& packet:products) for (const auto component:packet.velocityMS) native<<','<<component;
    native<<','<<products[0].representedNumber<<','<<products[1].representedNumber
          <<','<<bath.representedNumber<<','<<identities<<','<<event.shellQJ<<','<<event.radiatedEnergyJ
          <<','<<ledger.bathEnergyChangeJ<<','<<event.omittedImpulseBoundKgMS<<','<<event.sodiumIonRecoilBoundJ
          <<','<<event.photonRecoilSpeedBoundMS<<'\n';
  }
  native.close();Require(!native.fail(),"cannot write parent histories");
  return 0;
}
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank=0,size=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  int local=0;
  try { local=Run(argc,argv,rank,size); }
  catch (const std::exception& e) {std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';local=2;}
  int code=0;MPI_Allreduce(&local,&code,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Finalize();return code;
}

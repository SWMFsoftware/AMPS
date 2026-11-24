// MPI partitions original parent IDs; rank has no role in the physical draws.
// The adapter feeds source-bound inputs to normal production functions and
// writes their outputs. Python independently scores first-hit times and recoil.
#include "ground_capture_setup.h"
#include "srcMoon/species/QualificationRandom.h"
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace S=Moon::Species::SodiumGroundCapture;
namespace K=Moon::Species::ReducedKinematics;
namespace R=Moon::Species::ReactionExecution;
namespace G=Moon::Species::QualificationRandom;
namespace DB=AMPS::ChemicalDatabase;

namespace {
using I24Setup::Require;
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

int Run(int argc,char** argv,int rank,int size) {
  if (argc!=26) throw std::runtime_error("expected pinned bundle, source, beam and fixed cohort arguments");
  S::Definition definition;
  S::Provider provider;
  I24Setup::Load(argc,argv,&definition,&provider);
  S::Beam beam;
  beam.protonDensityM3=std::stold(argv[16]);
  const long double speed=std::stold(argv[17]);
  // A non-axis-aligned incident direction exercises the pair-frame rotation.
  beam.protonVelocityMS={{0.6L*speed,0.8L*speed,0}};
  R::ParentState parent;
  parent.species=0;parent.referenceWeight=std::stod(argv[23]);parent.individualCorrection=1;
  const std::array<double,2> weights={{std::stod(argv[24]),std::stod(argv[25])}};
  S::Evaluation rate;
  Require(provider.Evaluate(beam,parent.velocityMS,&rate));
  const auto output=std::string(argv[22]);
  if (!rank) {
    S::Beam trial=beam;
    auto prior=rate;
    trial.protonVelocityMS[0]*=1.01L;
    const bool rejected=!provider.Evaluate(trial,parent.velocityMS,&prior).ok() && prior.rateS==rate.rateS;
    trial=beam;trial.protonDensityM3=-1;
    const bool negative=!provider.Evaluate(trial,parent.velocityMS,&prior).ok() && prior.rateS==rate.rateS;
    trial.protonDensityM3=0;
    const bool zero=provider.Evaluate(trial,parent.velocityMS,&prior).ok() && prior.rateS==0;
    std::ofstream receipt(output+"/rate.csv");
    receipt<<std::setprecision(22)<<"rate_s,sigma_m2,speed_m_s,off_energy_refusal,negative_density_refusal,zero_density_zero_rate\n"
           <<rate.rateS<<','<<rate.crossSectionM2<<','<<rate.relativeSpeedMS<<','<<rejected<<','<<negative<<','<<zero<<'\n';
    receipt.close();Require(!receipt.fail(),"cannot write rate receipt");
  }
  const double horizon=std::stod(argv[18]);
  G::Key key;
  key.campaignSeed=std::stoull(argv[19]);key.cohort=argv[20];
  const auto parents=std::stoull(argv[21]);
  R::TraceInterval initial;
  std::string error;
  Require(R::TraceInterval::Create(0,horizon,0,&initial,&error),error);
  std::ofstream native(output+"/parents-rank-"+std::to_string(rank)+".csv");
  native<<std::setprecision(21)
      <<"parent,fate,time_s,nax,nay,naz,hx,hy,hz,na_weight,h_weight,bath_weight,identity_ok,Q_J\n";
  const std::vector<R::Hazard> hazards={{0,static_cast<double>(rate.rateS)}};
  bool testedGuards=false;
  for (key.originalParent=rank;key.originalParent<parents;key.originalParent+=size) {
    const auto origin=Clock(Draw(key,G::Purpose::SourceClock,0));
    R::HazardAdvance advance;
    Require(R::AdvanceConstantHazards(hazards,initial,origin,Draw(key,G::Purpose::Channel,0),&advance,&error),error);
    native<<key.originalParent<<',';
    if (!advance.event) {
      native<<"survived,"<<horizon<<",,,,,,,,,,,\n";
      continue;
    }
    const std::array<R::OpticalDepthClock,2> fresh={{Clock(Draw(key,G::Purpose::ProductClock,0)),
                                                     Clock(Draw(key,G::Purpose::ProductClock,1))}};
    S::Composition event;
    const auto polar=Draw(key,G::Purpose::Kinematics,0), azimuth=Draw(key,G::Purpose::Kinematics,1);
    Require(provider.Compose(beam,parent,advance.interval,polar,azimuth,weights,fresh,&event));
    const auto& products=event.plan.products();
    const auto& ledger=event.plan.ledger();
    Require(products.size()==2 && ledger.bathTerms.size()==1,"unreadable product/bath layout");
    const auto& bath=ledger.bathTerms.front();
    if (rank==0 && !testedGuards) {
      // Reuse the identical first accepted event for the rejection and
      // reweighting controls. No parent draw or provider state is consumed.
      auto saved=event;
      const bool endpoint=!provider.Compose(beam,parent,advance.interval,0,azimuth,weights,fresh,&saved).ok() &&
           saved.plan.valid() && saved.plan.reactionIdentity()==event.plan.reactionIdentity() &&
           saved.plan.products()[0].velocityMS==products[0].velocityMS;
      auto scaled=weights;
      for (auto& w:scaled) w*=2;
      Require(provider.Compose(beam,parent,advance.interval,polar,azimuth,scaled,fresh,&saved));
      bool invariant=saved.plan.products().size()==2;
      for(std::size_t j=0; j<saved.plan.products().size() && invariant; ++j) {
        const auto& packet=saved.plan.products()[j];
        invariant=packet.referenceWeight==scaled[j] &&
            packet.individualCorrection==products[j].individualCorrection/2 &&
            packet.representedNumber==products[j].representedNumber && packet.velocityMS==products[j].velocityMS;
      }
      std::ofstream guards(output+"/event-controls.csv");
      guards<<"closed_endpoint_refused_prior_plan_preserved,rescaled_daughter_weights_invariant\n"
            <<endpoint<<','<<invariant<<'\n';
      guards.close();Require(!guards.fail(),"cannot write provider controls");
      testedGuards=true;
    }
    // Expose identities as actually committed to the immutable plan. Neither
    // the native flags alone nor an independently sampled scattering histogram
    // proves that the planned velocities conserve the incoming four species.
    const bool identities=event.plan.valid() && event.plan.scope()==R::PlanScope::Transport &&
       products[0].compiledSpecies==definition.sodiumIonSpecies &&
       products[1].compiledSpecies==definition.hydrogenSpecies &&
       bath.species==1 && bath.reactant && ledger.integerConservationVerified && ledger.kinematicConservationVerified &&
       ledger.reactantChargeElementary==1 && ledger.productChargeElementary==1 &&
       ledger.trackedChargeChangeElementary==1 && ledger.bathChargeChangeElementary==-1;
    native<<"capture,"<<advance.eventOffsetS;
    for (const auto& packet:products) for (const auto component:packet.velocityMS) native<<','<<component;
    native<<','<<products[0].representedNumber<<','<<products[1].representedNumber
          <<','<<bath.representedNumber<<','<<identities<<','<<event.products.stateQJ<<'\n';
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

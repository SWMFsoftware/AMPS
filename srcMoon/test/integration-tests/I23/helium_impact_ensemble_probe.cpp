// I23 homogeneous electron-bath ensemble through the production provider.
//
// MPI partitions original parent IDs; rank has no role in the physical draws.
// For every parent this adapter calls the SAME production services used by
// the linked owner: HeliumElectronImpact::Evaluate for n_e*k(T_e), the
// production constant-hazard optical clock, and the guarded reviewed
// HeliumElectronImpact::Compose for the He+ product and signed electron
// reservoir. It writes their outputs; Python independently scores first-hit
// times, endpoint lifetime and every reservoir moment. No PIC mesh is used.
#include "helium_impact_setup.h"
#include "srcMoon/species/QualificationRandom.h"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

namespace G=Moon::Species::QualificationRandom;
using namespace I23Setup;

namespace {
double Draw(const G::Key& key,G::Purpose purpose,unsigned index) {
  std::string error;double value=0;
  Require(G::OpenUniform(key,purpose,index,&value,&error),error);
  return value;
}
R::OpticalDepthClock Clock(double uniform) {
  R::OpticalDepthClock clock;std::string error;
  Require(R::OpticalDepthClock::FromOpenUniform(uniform,&clock,&error),error);
  return clock;
}
template<class V> void Vector(std::ostream& out,const V& v){for(int i=0;i<3;++i)out<<','<<v[i];}

int Run(int argc,char** argv,int rank,int size) {
  const auto a=Parse(argc,argv);
  const std::string frame=a["frame"],output=a["output"];
  P::Registry profile;LoadProfile(a,&profile);
  H::Provider provider;
  const auto prepared=Prepare(a,profile,0,1,&provider);
  if(prepared.code==H::Code::Unavailable){std::cerr<<prepared.message<<'\n';return 77;}
  Require(prepared);
  const auto density=a.Real("density_m3"),horizon=a.Real("horizon_s");
  const auto bath=Bath(a,density,frame);
  const auto admission=Admission(a,profile,frame);
  H::Evaluation rate;
  Require(provider.Evaluate(bath,&rate));
  if(!rank){
    // Explicit refusal controls. Each must refuse and preserve the prior
    // output object; none may be accepted by substituting a default.
    auto guard=[&](H::Bath trial,H::Code expected){
      H::Evaluation prior=rate;const auto status=provider.Evaluate(trial,&prior);
      return status.code==expected&&prior.rateS==rate.rateS&&prior.coefficientM3S==rate.coefficientM3S;
    };
    auto trial=bath;trial.densityM3=-1;const bool negative=guard(trial,H::Code::InvalidInput);
    trial=bath;trial.densityM3=std::numeric_limits<double>::quiet_NaN();const bool nanDensity=guard(trial,H::Code::InvalidInput);
    trial=bath;trial.thermalEnergyEv=std::numeric_limits<double>::quiet_NaN();const bool nanTemperature=guard(trial,H::Code::InvalidInput);
    trial=bath;trial.thermalEnergyEv=-10;const bool negativeTemperature=guard(trial,H::Code::InvalidInput);
    trial=bath;trial.thermalEnergyEv=11;const bool unapprovedTemperature=guard(trial,H::Code::OutsideDomain);
    trial=bath;trial.bulkVelocityMS[0]=1;const bool drift=guard(trial,H::Code::OutsideDomain);
    trial=bath;trial.targetVelocityMS[2]=1;const bool movingTarget=guard(trial,H::Code::OutsideDomain);
    trial=bath;trial.frameIdentity.clear();const bool noFrame=guard(trial,H::Code::InvalidInput);
    std::ofstream receipt(output+"/rate.csv");
    receipt<<std::setprecision(std::numeric_limits<long double>::max_digits10)
      <<"density_m3,coefficient_m3_s,rate_s,log_slope,mean_incoming_energy_J,fit_threshold_eV,physical_threshold_J,"
        "negative_density_refused,nan_density_refused,nan_temperature_refused,negative_temperature_refused,"
        "unapproved_temperature_refused,drift_refused,moving_target_refused,missing_frame_refused\n"
      <<static_cast<double>(density)<<','<<rate.coefficientM3S<<','<<rate.rateS<<','<<rate.logSlope<<','
      <<rate.meanIncomingEnergyJ<<','<<rate.fitThresholdEv<<','<<rate.physicalThresholdJ<<','
      <<negative<<','<<nanDensity<<','<<nanTemperature<<','<<negativeTemperature<<','<<unapprovedTemperature<<','
      <<drift<<','<<movingTarget<<','<<noFrame<<'\n';
    receipt.close();Require(!receipt.fail(),"cannot write rate receipt");
  }
  G::Key key;key.campaignSeed=std::stoull(a["seed"]);key.cohort=a["cohort"];
  const auto parents=std::stoull(a["parents"]);
  const double parentWeight=static_cast<double>(a.Real("parent_weight"));
  const double ionWeight=static_cast<double>(a.Real("ion_weight"));
  R::TraceInterval initial;std::string error;
  Require(R::TraceInterval::Create(0,static_cast<double>(horizon),0,&initial,&error),error);
  std::ofstream native(output+"/parents-rank-"+std::to_string(rank)+".csv");
  const std::string header="parent,fate,time_s,ion_species,ion_vx,ion_vy,ion_vz,ion_reference_weight,ion_correction,ion_represented,"
      "bath_reactant_represented,bath_product_represented,integer_ok,kinematic_ok,tracked_charge_e,bath_charge_e,"
      "ledger_bath_mass_kg,ledger_bath_energy_J,ledger_tracked_energy_J,ledger_bath_px,ledger_bath_py,ledger_bath_pz,"
      "in_number,in_mass_kg,in_charge_C,in_kinetic_J,in_px,in_py,in_pz,out_number,out_mass_kg,out_charge_C,out_kinetic_J,"
      "out_px,out_py,out_pz,state_Q_J,bath_number_change,bath_mass_change_kg,bath_charge_change_C,bath_energy_change_J,"
      "incoming_mean_J,heavy_velocity_bound_m_s,heavy_kinetic_bound_J,relative_mass_error,represented_bath_charge_C,"
      "represented_bath_energy_J\n";
  native<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<header;
  // A survivor has no event fields; pad to the same column count.
  const auto columns=static_cast<std::size_t>(std::count(header.begin(),header.end(),','));
  const std::vector<R::Hazard> hazards={{0,rate.rateS}};
  bool testedGuards=false;
  for(key.originalParent=rank;key.originalParent<parents;key.originalParent+=size){
    const auto origin=Clock(Draw(key,G::Purpose::SourceClock,0));
    R::HazardAdvance advance;
    Require(R::AdvanceConstantHazards(hazards,initial,origin,Draw(key,G::Purpose::Channel,0),&advance,&error),error);
    native<<key.originalParent<<',';
    if(!advance.event){native<<"survived,"<<static_cast<double>(horizon)<<std::string(columns-2,',')<<'\n';continue;}
    R::ParentState parent;parent.species=0;parent.referenceWeight=parentWeight;parent.individualCorrection=1;
    const auto fresh=Clock(Draw(key,G::Purpose::ProductClock,0));
    const double epoch=advance.interval.startS()+advance.interval.elapsedS();
    H::Composition event;
    Require(provider.Compose(bath,admission,parent,advance.interval,epoch,ionWeight,fresh,&event));
    if(rank==0&&!testedGuards){
      // Reuse the first accepted event. No parent draw is consumed.
      auto saved=event;
      const auto refused=[&](const H::Admission& changed,const R::ParentState& p,H::Code code){
        auto copy=event;const auto status=provider.Compose(bath,changed,p,advance.interval,epoch,ionWeight,fresh,&copy);
        return status.code==code&&copy.plan.valid()&&copy.plan.reactionIdentity()==event.plan.reactionIdentity();
      };
      const bool missing=refused(H::Admission{},parent,H::Code::Unavailable);
      auto tight=admission;tight.evidence.maximumHeavyVelocityErrorMS=a.Real("central_heavy_velocity_bound_m_s");
      const bool centralRefused=refused(tight,parent,H::Code::OutsideDomain);
      auto corrupt=admission;corrupt.evidence.modelErrorEvidenceSha256=std::string(64,'0');
      const bool corruptRefused=refused(corrupt,parent,H::Code::CorruptSource);
      auto moving=parent;moving.velocityMS[0]=1;
      const bool movingRefused=refused(admission,moving,H::Code::OutsideDomain);
      Require(provider.Compose(bath,admission,parent,advance.interval,epoch,2*ionWeight,fresh,&saved));
      const auto& p0=event.plan.products()[0];const auto& p1=saved.plan.products()[0];
      const bool invariant=saved.plan.products().size()==1&&p1.referenceWeight==2*ionWeight&&
          p1.individualCorrection==p0.individualCorrection/2&&p1.representedNumber==p0.representedNumber&&p1.velocityMS==p0.velocityMS;
      std::ofstream guards(output+"/event-controls.csv");
      guards<<"missing_admission_refused,central_velocity_tolerance_refused,corrupt_evidence_refused,"
              "moving_parent_refused,rescaled_ion_weight_invariant\n"
            <<missing<<','<<centralRefused<<','<<corruptRefused<<','<<movingRefused<<','<<invariant<<'\n';
      guards.close();Require(!guards.fail(),"cannot write provider controls");
      testedGuards=true;
    }
    const auto& products=event.plan.products();
    const auto& l=event.plan.ledger();
    Require(products.size()==1&&l.bathTerms.size()==2,"unreadable He+/bath plan layout");
    double reactantBath=0,productBath=0;
    for(const auto& term:l.bathTerms)(term.reactant?reactantBath:productBath)+=term.representedNumber;
    const auto& ion=products[0];const auto& rr=event.reservoir;
    native<<"ionized,"<<advance.eventOffsetS<<','<<ion.compiledSpecies;Vector(native,ion.velocityMS);
    native<<','<<ion.referenceWeight<<','<<ion.individualCorrection<<','<<ion.representedNumber<<','
          <<reactantBath<<','<<productBath<<','<<l.integerConservationVerified<<','<<l.kinematicConservationVerified<<','
          <<l.trackedChargeChangeElementary<<','<<l.bathChargeChangeElementary<<','<<l.bathModelMassChangeKg<<','
          <<l.bathEnergyChangeJ<<','<<l.trackedEnergyChangeJ;Vector(native,l.bathMomentumChangeKgMS);
    for(const auto* m:{&rr.incomingBath,&rr.outgoingBath}){
      native<<','<<m->number<<','<<m->massKg<<','<<m->chargeC<<','<<m->kineticEnergyJ;Vector(native,m->momentumKgMS);
    }
    native<<','<<rr.stateQJ<<','<<rr.bathNumberChange<<','<<rr.bathMassChangeKg<<','<<rr.bathChargeChangeC<<','
          <<rr.bathEnergyChangeJ<<','<<rr.incomingConditionalMeanEnergyJ<<','<<rr.error.heavyVelocityUpperMS<<','
          <<rr.error.heavyKineticMeanAbsoluteUpperJ<<','<<rr.error.relativeMassError<<','
          <<rr.representedBathChargeChangeC<<','<<rr.representedBathEnergyChangeJ<<'\n';
  }
  native.close();Require(!native.fail(),"cannot write parent histories");
  return 0;
}
}
int main(int argc,char** argv){
  MPI_Init(&argc,&argv);
  int rank=0,size=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  int local=0;
  try{local=Run(argc,argv,rank,size);}
  catch(const std::exception& e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';local=2;}
  int code=0;MPI_Allreduce(&local,&code,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Finalize();return code;
}

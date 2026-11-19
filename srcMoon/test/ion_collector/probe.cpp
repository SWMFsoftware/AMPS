// Analytical numerical-boundary fixture only, NOT a lunar ion-surface model.
// All masses/energies/charges below are explicitly chosen dimensioned scales.
// The executable calls production IonCollector and the actual runtime parser.
#include "species/IonCollector.h"
#include "MoonInput.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
using namespace Moon::Species;
namespace C=IonCollector;
static int failures=0;
void Check(const char* name,bool ok) { std::cout<<"CHECK "<<name<<' '<<ok<<'\n';if(!ok)++failures; }
C::Incident Incident() {
  C::Incident x;
  x.state.material={"analytical-A2B-ion","declared-excited",0,{{"A:1",2},{"B:2",1}},2,4.0};
  x.state.formationEnergyJ=-7;x.state.excitationEnergyJ=3;
  x.state.energyReferenceIdentity="explicit-analytical-zero";
  x.state.sourceEvidenceIdentity="analytical-bookkeeping-fixture-not-lunar-data";
  x.representedNumber=2.5;x.elementaryChargeC=0.5;
  x.positionLsoM={{10,20,30}};x.velocityLsoMS={{3,-4,-5}};
  x.impactEpochEtS=101.25;x.frame.geometryEpochEtS=100;
  x.frame.geometryIdentity="analytical-moving-facet";
  x.frame.outwardNormalLso={{0,0,2}};x.frame.surfaceVelocityLsoMS={{1,-2,1}};
  x.frame.stateTransformsSupplied=true;
  for(int i=0;i<6;++i)x.frame.lsoToBody36[7*i]=x.frame.bodyToLso36[7*i]=1;
  return x;
}
int main(int argc,char** argv) {
  if(argc==3 && std::string(argv[1])=="parse") {
    Moon::Runtime::Configuration c;c.randomSeed=-777;std::string error;
    if(!Moon::Runtime::ParseApplicationInput(argv[2],&c,&error)) {
      std::cout<<"preserved="<<(c.randomSeed==-777)<<'\n';std::cerr<<error<<'\n';return 3;
    }
    std::cout<<"ion_surface_boundary="<<c.ionSurfaceBoundary<<'\n';return 0;
  }
  if(argc!=2)return 2;
  std::string error;auto incident=Incident();C::Ledger ledger;
  Check("collect-explicit",C::Collect(incident,&ledger,&error));
  Check("one-publication",ledger.records.size()==1);
  if(ledger.records.size()!=1)return 1;
  const auto record=ledger.records.front();
  Check("identity-state",record.incident.state.material.identity==incident.state.material.identity && record.incident.state.material.state==incident.state.material.state);
  Check("integer-nuclei",record.nuclei.size()==2 && record.nuclei[0].multiplicity==2 && record.nuclei[1].multiplicity==1);
  Check("actual-two-epochs",record.incident.impactEpochEtS==101.25 && record.incident.frame.geometryEpochEtS==100);
  Check("frame-normal-retained",record.incident.frame.outwardNormalLso==incident.frame.outwardNormalLso && record.incident.frame.lsoToBody36==incident.frame.lsoToBody36);
  Check("energy-provenance",record.incident.state.energyReferenceIdentity==incident.state.energyReferenceIdentity && record.incident.state.sourceEvidenceIdentity==incident.state.sourceEvidenceIdentity);
  // The file exposes maintained-kernel output; Python uses independent simple
  // algebra to compare these quantities and preserves every comparison row.
  std::ofstream csv(argv[1]);csv<<"quantity,actual\n"<<std::setprecision(24);
  csv<<"nucleus_A,"<<record.nuclei[0].representedNumber<<'\n'
     <<"nucleus_B,"<<record.nuclei[1].representedNumber<<'\n'
     <<"charge_elementary,"<<record.representedChargeElementary<<'\n'
     <<"charge_C,"<<record.chargeC<<'\n'
     <<"mass_kg,"<<record.modelMassKg<<'\n'
     <<"momentum_x,"<<record.momentumFacetKgMS[0]<<'\n'
     <<"momentum_y,"<<record.momentumFacetKgMS[1]<<'\n'
     <<"momentum_z,"<<record.momentumFacetKgMS[2]<<'\n'
     <<"kinetic_facet_J,"<<record.kineticEnergyFacetJ<<'\n'
     <<"kinetic_lso_J,"<<record.kineticEnergyLsoJ<<'\n'
     <<"formation_J,"<<record.formationEnergyJ<<'\n'
     <<"excitation_J,"<<record.excitationEnergyJ<<'\n'
     <<"total_facet_J,"<<record.totalEnergyFacetJ<<'\n';
  csv.close();if(!csv)return 4;
  auto boosted=incident;
  for(int d=0;d<3;++d) {boosted.velocityLsoMS[d]+=8;boosted.frame.surfaceVelocityLsoMS[d]+=8;}
  C::Record shifted;
  Check("moving-facet-galilean",C::Prepare(boosted,&shifted,&error) && shifted.kineticEnergyFacetJ==record.kineticEnergyFacetJ && shifted.momentumFacetKgMS==record.momentumFacetKgMS && shifted.totalEnergyFacetJ==record.totalEnergyFacetJ);
  auto reject=[&](const char* name,const C::Incident& bad) {
    const auto size=ledger.records.size();C::Record unchanged=record;
    const bool prepared=C::Prepare(bad,&unchanged,&error);
    const bool collected=C::Collect(bad,&ledger,&error);
    Check(name,!prepared && !collected && !error.empty() && ledger.records.size()==size && unchanged.incident.impactEpochEtS==record.incident.impactEpochEtS && unchanged.totalEnergyFacetJ==record.totalEnergyFacetJ);
  };
  auto bad=incident;bad.state.formationEnergyJ=NAN;reject("reject-missing-formation",bad);
  bad=incident;bad.state.excitationEnergyJ=NAN;reject("reject-missing-excitation",bad);
  bad=incident;bad.state.excitationEnergyJ=-1;reject("reject-negative-excitation",bad);
  bad=incident;bad.state.energyReferenceIdentity.clear();reject("reject-missing-energy-zero",bad);
  bad=incident;bad.state.sourceEvidenceIdentity.clear();reject("reject-missing-evidence",bad);
  bad=incident;bad.state.material.state.clear();reject("reject-missing-state",bad);
  bad=incident;bad.state.material.elementaryCharge=0;reject("reject-neutral",bad);
  bad=incident;bad.state.material.compiledSpecies=-1;reject("reject-uncompiled",bad);
  bad=incident;bad.state.material.modelMassKg=0;reject("reject-zero-mass",bad);
  bad=incident;bad.state.material.nuclei[0].count=0;reject("reject-invalid-nuclei",bad);
  bad=incident;bad.state.material.nuclei.clear();reject("reject-electron-without-nuclei",bad);
  bad=incident;bad.representedNumber=0;reject("reject-zero-weight",bad);
  bad=incident;bad.representedNumber=-1;reject("reject-negative-weight",bad);
  bad=incident;bad.representedNumber=INFINITY;reject("reject-nonfinite-weight",bad);
  bad=incident;bad.elementaryChargeC=0;reject("reject-missing-charge-unit",bad);
  bad=incident;bad.impactEpochEtS=NAN;reject("reject-missing-impact-epoch",bad);
  bad=incident;bad.frame.geometryEpochEtS=NAN;reject("reject-missing-geometry-epoch",bad);
  bad=incident;bad.positionLsoM[0]=NAN;reject("reject-nonfinite-position",bad);
  bad=incident;bad.velocityLsoMS[1]=INFINITY;reject("reject-nonfinite-velocity",bad);
  bad=incident;bad.frame.surfaceVelocityLsoMS[2]=NAN;reject("reject-nonfinite-surface-velocity",bad);
  bad=incident;bad.frame.outwardNormalLso={{0,0,0}};reject("reject-zero-normal",bad);
  bad=incident;bad.velocityLsoMS[2]=2;reject("reject-outgoing",bad);
  bad=incident;bad.frame.geometryIdentity.clear();reject("reject-missing-geometry",bad);
  bad=incident;bad.frame.stateTransformsSupplied=false;reject("reject-missing-transforms",bad);
  bad=incident;bad.frame.lsoToBody36[0]=NAN;reject("reject-nonfinite-transform",bad);
  bad=incident;bad.frame.bodyToLso36[0]=0;reject("reject-noninvertible-transform",bad);
  auto second=incident.state;second.material.state="another-state";second.excitationEnergyJ=0;
  Check("explicit-zero-excitation",C::ValidateStates({second},&error));
  Check("distinct-explicit-states",C::ValidateStates({incident.state,second},&error));
  Check("reject-duplicate-state",!C::ValidateStates({incident.state,incident.state},&error));
  second.energyReferenceIdentity="different-zero";
  Check("reject-mixed-energy-zero",!C::ValidateStates({incident.state,second},&error));
  auto mixed=incident;mixed.state.energyReferenceIdentity="different-zero";
  const auto count=ledger.records.size();
  Check("reject-ledger-energy-zero-change",!C::Collect(mixed,&ledger,&error) && ledger.records.size()==count);
  Check("reject-empty-inventory",!C::ValidateStates({},&error));
  Check("reject-null-ledger",!C::Collect(incident,nullptr,&error));
  auto negative=incident;negative.state.material.identity="negative-analytical-ion";negative.state.material.elementaryCharge=-1;C::Record negativeRecord;
  Check("negative-ion-signed-charge",C::Prepare(negative,&negativeRecord,&error) && negativeRecord.chargeC==-1.25L);
  auto zeroFormation=incident;zeroFormation.state.formationEnergyJ=0;
  Check("explicit-zero-formation",C::Prepare(zeroFormation,&shifted,&error) && shifted.formationEnergyJ==0);
  std::cout<<"SUMMARY failed="<<failures<<'\n';return failures?1:0;
}

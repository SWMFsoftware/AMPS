// Analytical nominal-mass controls ONLY; these are not Na/He measurements or
// an authorized physical source profile. The real shared product services are
// called with explicit identities so every conservation/refusal path is real.
#include "../../src/models/chemical_database/ReservoirProducts.h"
#include <iomanip>
#include <iostream>
#include <string>
namespace R=AMPS::ChemicalDatabase::ReservoirProducts;
namespace D=AMPS::ChemicalDatabase;
static R::State state(std::string id,long double mass,int charge,long double energy,std::string nucleus="A") {
 R::State s;s.identity=id;s.state="fixture-state";s.sourceIdentity="ANALYTICAL-NOT-OBSERVATIONS";
 s.energyReferenceIdentity="fixture-zero";s.kineticMassKg=mass;s.physicalRestMassKg=mass;
 s.chargeElementary=charge;s.stateEnergyJ=energy;if(!nucleus.empty())s.nuclei.push_back({nucleus,1});return s;
}
static R::ElectronImpactInput input() {
 R::ElectronImpactInput x;auto& e=x.evidence;
 e.modelIdentity="analytical-reservoir-control";e.sourceProfileIdentity="fixture-profile";e.sourceProfileSha256=std::string(64,'a');
 e.stateCrosswalkIdentity="fixture-states";e.stateCrosswalkSha256=std::string(64,'b');e.domainEvidenceIdentity="fixture-domain";
 e.modelErrorEvidenceIdentity="fixture-local-bound";e.modelErrorEvidenceSha256=std::string(64,'c');
 e.frameIdentity="fixture-inertial";e.timeScaleIdentity="fixture-seconds";e.energyReferenceIdentity="fixture-zero";e.observableIdentity="fixture-local-speed";
 e.maximumHeavyVelocityErrorMS=100;e.maximumHeavyEnergyErrorJ=1000;e.maximumRelativeMassError=0;e.maximumClassicalHeavySpeedMS=100;
 x.event.positionM={{3,4,5}};x.event.parentVelocityMS={{1,2,0}};x.event.physicalStartS=100;x.event.elapsedS=.375L;x.event.remainingS=.625L;x.event.physicalEventS=100.375L;x.event.representedNumber=7.5L;x.event.frameIdentity=e.frameIdentity;x.event.timeScaleIdentity=e.timeScaleIdentity;
 x.numerical.massAbsoluteKg=0;x.numerical.massRelative=1e-17L;x.numerical.timeProjectionAbsoluteS=0;
 x.numerical.conservation={1e-14L,1e-17L,1e-14L,1e-17L,1e-17L};
 x.parent=state("A",11,0,0);x.ion=state("A+",10,1,3);x.electron=state("electron",1,-1,0,"");
 x.elementaryChargeC=1.602176634e-19L;x.sourceThresholdCrosswalkMaximumDifferenceJ=0;
 x.incoming.derivativeSourceIdentity="constant-cross-section-with-explicit-threshold-analytic";x.incoming.temperatureDomainIdentity="fixture-range";
 x.incoming.bulkVelocityParentRestMS={{0,0,0}};x.incoming.thermalEnergyJ=2;x.incoming.minimumThermalEnergyJ=1;x.incoming.maximumThermalEnergyJ=3;x.incoming.logarithmicRateDerivative=1.4L;
 x.incoming.thresholdJ=3;x.incoming.isotropicUnshiftedParentRest=true;x.incoming.crossSectionTemperatureIndependent=true;x.incoming.collisionWeightRotationallyInvariant=true;return x;
}
int main(int argc,char** argv) {
 const std::string mode=argc>1?argv[1]:"ei";
 auto x=input();R::Result out;out.scope="SENTINEL";D::Status status;
 if(mode.rfind("cex",0)==0) {
  R::ForwardChargeExchangeInput c;c.evidence=x.evidence;c.numerical=x.numerical;c.event=x.event;
  c.parent=state("A",11,0,0);c.heavyIon=state("A+",10,1,3);
  c.plasmaIon=state("B+",2,1,5,"B");c.fastNeutral=state("B",3,0,0,"B");
  c.plasmaVelocityMS={{5,2,0}};c.elementaryChargeC=x.elementaryChargeC;
  c.sourceEnergyDomainMappingIdentity="fixture-COM-J";c.minimumIncomingRelativeKineticEnergyJ=1;c.maximumIncomingRelativeKineticEnergyJ=100;
  c.partialChannelIdentity="fixture-partial";c.angleSupportEvidenceIdentity="explicit-analytical-angle-domain";c.sourceSupportedMaximumCOMDeflectionRad=0;
  if(mode=="cex-no-nuclei"){c.parent.nuclei.clear();c.heavyIon.nuclei.clear();c.plasmaIon.nuclei.clear();c.fastNeutral.nuclei.clear();}
  if(mode=="cex-not-transfer"){c.plasmaIon.chargeElementary=0;c.heavyIon.chargeElementary=0;}
  if(mode=="cex-nuclei-overflow"){c.parent.nuclei[0].count=std::numeric_limits<std::int64_t>::max();c.heavyIon.nuclei=c.parent.nuclei;c.plasmaIon.nuclei[0].identity="A";c.fastNeutral.nuclei[0].identity="A";}
  if(mode=="cex-zero")c.plasmaVelocityMS=c.event.parentVelocityMS;
  if(mode=="cex-angle")c.sourceSupportedMaximumCOMDeflectionRad=3.141592653589793238462643383279502884L;
  if(mode=="cex-impossible")c.heavyIon.stateEnergyJ=1000;
  if(mode=="cex-domain")c.maximumIncomingRelativeKineticEnergyJ=1;
  if(mode=="cex-missing-angle")c.angleSupportEvidenceIdentity.clear();
  if(mode=="cex-missing-neutral")c.fastNeutral.state.clear();
  if(mode=="cex-bad-nuclei")c.fastNeutral.nuclei[0].identity="C";
  if(mode=="cex-bad-angle")c.sourceSupportedMaximumCOMDeflectionRad=4;
  if(mode=="cex-rotated"){c.event.parentVelocityMS={{0,1,2}};c.plasmaVelocityMS={{0,5,2}};}
  if(mode=="cex-boosted"){c.event.parentVelocityMS[0]+=3;c.plasmaVelocityMS[0]+=3;}
  status=R::ForwardPartialChargeExchange(c,&out);
 } else {
  if(mode=="infinite-derivative")x.incoming.logarithmicRateDerivative=std::numeric_limits<long double>::infinity();
  if(mode=="missing")x.evidence.stateCrosswalkSha256.clear();
  if(mode=="bad-hash")x.evidence.sourceProfileSha256=std::string(64,'z');
  if(mode=="bad-reference")x.ion.energyReferenceIdentity="OTHER";
  if(mode=="missing-Q")x.ion.stateEnergyJ=std::numeric_limits<long double>::quiet_NaN();
  if(mode=="anisotropic")x.incoming.collisionWeightRotationallyInvariant=false;
  if(mode=="rotated")x.event.parentVelocityMS={{0,1,2}};
  if(mode=="boosted")x.event.parentVelocityMS[0]+=3;
  if(mode=="N0")x.event.representedNumber=0;
  if(mode=="Noverflow")x.event.representedNumber=std::numeric_limits<long double>::max();
  if(mode=="Nunderflow")x.event.representedNumber=std::numeric_limits<long double>::denorm_min();
  if(mode=="source-mismatch")x.incoming.thresholdJ=2;
  if(mode=="charge")x.ion.chargeElementary=2;
  if(mode=="energy-budget")x.evidence.maximumHeavyEnergyErrorJ=0;
  if(mode=="mass-budget")x.ion.physicalRestMassKg=11;
  if(mode=="drift")x.incoming.isotropicUnshiftedParentRest=false;
  if(mode=="drift-vector")x.incoming.bulkVelocityParentRestMS={{1,0,0}};
  if(mode=="threshold")x.incoming.thresholdJ=6;
  if(mode=="temperature")x.incoming.thermalEnergyJ=4;
  if(mode=="budget")x.evidence.maximumHeavyVelocityErrorMS=0;
  if(mode=="nuclei")x.ion.nuclei[0].identity="B";
  if(mode=="mass")x.ion.kineticMassKg=9;
  if(mode=="clock")x.event.physicalEventS=101;
  if(mode=="mean"){x.incoming.thresholdJ=0;x.incoming.logarithmicRateDerivative=.5L;}
  if(mode=="mean"||mode=="threshold-mean") {long double v=-1;status=R::MaxwellianCollisionMeanEnergy(x.incoming,&v);std::cout<<std::setprecision(25)<<"mean "<<v<<"\n";}
  else status=R::ElectronImpactInheritedVelocity(x,&out);
 }
 std::cout<<std::setprecision(25)<<"code "<<D::CodeName(status.code)<<"\nmessage "<<status.message<<"\nsentinel "<<(out.scope=="SENTINEL")<<'\n';
 if(status.ok()&&mode!="mean"&&mode!="threshold-mean") {
  auto scalar=[](const char* key,long double v){std::cout<<key<<' '<<v<<'\n';};
  auto vector=[](const char* key,const R::Vector3& v){std::cout<<key<<' '<<v[0]<<' '<<v[1]<<' '<<v[2]<<'\n';};
  scalar("conditional-mean",out.incomingConditionalMeanEnergyJ);scalar("incoming-relative-K",out.incomingRelativeKineticEnergyJ);
  scalar("weighted-charge",out.representedBathChargeChangeC);scalar("mass-error",out.error.relativeMassError);
  std::cout<<"incoming-id "<<out.incomingBath.speciesIdentity<<"\noutgoing-id "<<out.outgoingBath.speciesIdentity<<'\n';
  scalar("Q",out.stateQJ);scalar("incoming-number",out.incomingBath.number);scalar("outgoing-number",out.outgoingBath.number);
  scalar("incoming-K",out.incomingBath.kineticEnergyJ);scalar("outgoing-K",out.outgoingBath.kineticEnergyJ);
  scalar("net-E",out.bathEnergyChangeJ);scalar("net-charge",out.bathChargeChangeC);scalar("weighted-E",out.representedBathEnergyChangeJ);
  scalar("weight",out.event.representedNumber);scalar("elapsed",out.event.elapsedS);scalar("remaining",out.event.remainingS);scalar("epoch",out.event.physicalEventS);
  scalar("speed-error",out.error.heavyVelocityUpperMS);scalar("energy-error",out.error.heavyKineticMeanAbsoluteUpperJ);scalar("energy-residual",out.energyResidualJ);
  vector("v",out.trackedVelocityMS);vector("incoming-P",out.incomingBath.momentumKgMS);vector("outgoing-P",out.outgoingBath.momentumKgMS);
  vector("net-P",out.bathMomentumChangeKgMS);vector("p-residual",out.momentumResidualKgMS);
 }
 return 0;
}

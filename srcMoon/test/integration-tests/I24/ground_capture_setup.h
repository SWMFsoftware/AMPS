// Shared by the provider and live PIC probes. The Python source reader supplies
// argv[1..15] from the verified r14 release and NIST inventory. Both adapters
// therefore use the same physical states and explicit SI arithmetic controls.
#ifndef AMPS_I24_GROUND_CAPTURE_SETUP_H
#define AMPS_I24_GROUND_CAPTURE_SETUP_H
#include "srcMoon/species/SodiumGroundCapture.h"
#include <stdexcept>
namespace I24Setup {
namespace S=Moon::Species::SodiumGroundCapture;
namespace K=Moon::Species::ReducedKinematics;
namespace DB=AMPS::ChemicalDatabase;
inline void Require(bool value,const std::string& message) {
  if (!value) throw std::runtime_error(message);
}
inline void Require(const DB::Status& status) {Require(status.ok(),status.message);}

inline K::ResolvedState State(const char* id,const char* state,const char* nucleus,int charge,
    long double energy,const std::vector<K::IsotopeMass>& isotopes,const K::MassValue& electron) {
  // These energies and masses arrive from the independently authenticated NIST
  // inventory. The two neutral ground states define the formation-energy zero.
  K::StateDefinition definition;
  definition.identity=id;definition.state=state;definition.nuclei={{nucleus,1}};
  definition.elementaryCharge=charge;definition.formation={energy,"NIST-elemental-ground-IE"};
  definition.excitation={0,"ground-state-definition"};
  K::ResolvedState resolved;
  std::string error;
  Require(K::ResolveState(definition,isotopes,electron,&resolved,&error),error);
  return resolved;
}

inline void Load(int argc,char** argv,S::Definition* definitionOut,S::Provider* providerOut,
                 int sodium=0,int ion=1,int hydrogen=2,
                 const std::string& frame="source-target-rest-inertial") {
  Require(argc>=16 && definitionOut && providerOut,"missing I24 source arguments");
  DB::BundleRequest request;
  request.root=argv[1];request.exportPath=argv[2];request.expectedExportSha256=argv[3];
  request.expectedReleaseId=argv[4];request.expectedReleaseSha256=argv[5];
  request.expectedProfileId=argv[6];request.expectedProfileSha256=argv[7];
  DB::Catalogue catalogue;
  Require(DB::OpenRuntimeBundle(request,&catalogue));
  S::Definition definition;
  definition.atomicMassUnitKg=std::stold(argv[12]);definition.electronVoltJ=std::stold(argv[13]);
  definition.constantsIdentity="CODATA-2022/NIST";definition.frame=frame;
  std::vector<K::IsotopeMass> isotope={{"1H",{std::stold(argv[9]),"NIST-1H"}},
                                      {"23Na",{std::stold(argv[10]),"NIST-23Na"}}};
  K::MassValue electron={std::stold(argv[11]),"CODATA-2022-electron"};
  definition.states={{State("Na","3s","23Na",0,0,isotope,electron),
                      State("H+","bare-proton","1H",1,std::stold(argv[14])*definition.electronVoltJ,isotope,electron),
                      State("Na+","2p6","23Na",1,std::stold(argv[15])*definition.electronVoltJ,isotope,electron),
                      State("H","1s","1H",0,0,isotope,electron)}};
  definition.sodiumSpecies=sodium;definition.sodiumIonSpecies=ion;definition.hydrogenSpecies=hydrogen;
  // Arithmetic error allowances cover double packet velocities and cancellation
  // in components perpendicular to the beam. They are independent of the
  // physical-model uncertainty and of the percent-level rate/lifetime gates.
  definition.numerical.maximumClassicalSpeedMS=3e6L;
  definition.numerical.unitVectorTolerance=5e-14L;
  definition.numerical.referenceMassKg={0,5e-14L};
  definition.numerical.momentumKgMS={1e-36L,5e-14L};
  definition.numerical.energyJ={1e-31L,5e-14L};
  definition.numerical.thresholdJ=definition.numerical.energyJ;
  definition.ledgerBounds={{1e-40,5e-13},{1e-34,5e-13},{1e-29,5e-13}};
  S::Provider provider;
  Require(S::Provider::Prepare(catalogue,argv[8],definition,&provider));
  *definitionOut=definition;
  *providerOut=provider;
}
 }
#endif

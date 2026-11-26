// Shared by the I24 total-capture ensemble and PIC probes. The Python runner
// supplies argv[1..10] from the authenticated case sources and the NIST mass
// inventory, so both adapters build the same physical states and controls:
//   1 Morgan publication directory   2 shingal_table4.csv   3 hydrogen_levels.csv
//   4 m(1H) [kg]   5 m(23Na) [kg]   6 m_e [kg]   7 u [kg]   8 e [J/eV]
//   9 IE(H) [eV]   10 IE(Na) [eV]
#ifndef AMPS_I24_TOTAL_CAPTURE_SETUP_H
#define AMPS_I24_TOTAL_CAPTURE_SETUP_H
#include "srcMoon/species/SodiumTotalCapture.h"
#include <stdexcept>
namespace I24TotalSetup {
namespace S=Moon::Species::SodiumTotalCapture;
namespace K=Moon::Species::ReducedKinematics;
namespace DB=AMPS::ChemicalDatabase;
constexpr int SourceArguments=10;
inline void Require(bool value,const std::string& message) {
  if (!value) throw std::runtime_error(message);
}
inline void Require(const DB::Status& status) {Require(status.ok(),status.message);}

inline K::ResolvedState State(const char* id,const char* state,const char* nucleus,int charge,
    long double energy,const std::vector<K::IsotopeMass>& isotopes,const K::MassValue& electron) {
  // Energies and masses come from the authenticated NIST inventory. The two
  // neutral ground states define the formation-energy zero.
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
  Require(argc>SourceArguments && definitionOut && providerOut,"missing I24 total-capture source arguments");
  S::Definition definition;
  definition.morganPublicationDirectory=argv[1];
  definition.shingalTablePath=argv[2];definition.hydrogenLevelsPath=argv[3];
  definition.atomicMassUnitKg=std::stold(argv[7]);definition.electronVoltJ=std::stold(argv[8]);
  definition.constantsIdentity="CODATA-2022/NIST";definition.frame=frame;
  std::vector<K::IsotopeMass> isotope={{"1H",{std::stold(argv[4]),"NIST-1H"}},
                                      {"23Na",{std::stold(argv[5]),"NIST-23Na"}}};
  K::MassValue electron={std::stold(argv[6]),"CODATA-2022-electron"};
  definition.states={{State("Na","3s","23Na",0,0,isotope,electron),
                      State("H+","bare-proton","1H",1,std::stold(argv[9])*definition.electronVoltJ,isotope,electron),
                      State("Na+","2p6","23Na",1,std::stold(argv[10])*definition.electronVoltJ,isotope,electron),
                      State("H","1s","1H",0,0,isotope,electron)}};
  definition.sodiumSpecies=sodium;definition.sodiumIonSpecies=ion;definition.hydrogenSpecies=hydrogen;
  // Arithmetic allowances (not physical uncertainty), identical to the
  // existing 4 keV/u ground-capture fixture.
  definition.numerical.maximumClassicalSpeedMS=3e6L;
  definition.numerical.unitVectorTolerance=5e-14L;
  definition.numerical.referenceMassKg={0,5e-14L};
  definition.numerical.momentumKgMS={1e-36L,5e-14L};
  definition.numerical.energyJ={1e-31L,5e-14L};
  definition.numerical.thresholdJ=definition.numerical.energyJ;
  definition.ledgerBounds={{1e-40,5e-13},{1e-34,5e-13},{1e-29,5e-13}};
  S::Provider provider;
  Require(S::Provider::Prepare(definition,&provider));
  *definitionOut=definition;
  *providerOut=provider;
}
}
#endif

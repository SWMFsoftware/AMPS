// Thin adapter for the production runtime owner. No rate law, parser or
// particle-transfer implementation belongs here: expected values are computed
// independently in Python, while this file calls the normal production APIs.
#include "MoonInput.h"

#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char** argv) {
  if (argc<2 || argc>3) return 2;
  Moon::Runtime::Configuration configuration;
  // A failed native parse must leave a caller's previous configuration intact.
  configuration.randomSeed=-12345;
  std::string error;
  if (!Moon::Runtime::ParseApplicationInput(argv[1],&configuration,&error)) {
    std::cerr<<error<<'\n';
    std::cout<<"parse_preserved="<<(configuration.randomSeed==-12345)<<'\n';
    return 3;
  }
  if (argc==3 && std::string(argv[2])=="rollback") {
    // Direct programmatic callers also need transactional installation.
    // Deliberately fail BEFORE a valid retry; no fake data/observation is made.
    auto invalid=configuration;
    invalid.electronTemperatureEv=0.0;
    const bool rejected=!Moon::Runtime::InstallConfiguration(invalid,&error);
    if (!rejected || Moon::Runtime::HasConfiguration() ||
        Moon::Runtime::GetElectronImpactTable()!=nullptr) return 5;
    std::cout<<"failed_install_preserved=1\n";
  }
  if (!Moon::Runtime::InstallConfiguration(configuration,&error)) {
    std::cerr<<error<<'\n';
    std::cout<<"install_preserved="<<(!Moon::Runtime::HasConfiguration())<<'\n';
    return 4;
  }
  const auto& actual=Moon::Runtime::GetConfiguration();
  const auto* table=Moon::Runtime::GetElectronImpactTable();
  std::cout<<std::setprecision(17)
      <<"model="<<actual.electronImpactRateModel<<'\n'
      <<"native_table="<<actual.electronImpactCoefficientsFile<<'\n'
      <<"density_m3="<<actual.electronDensityM3<<'\n'
      <<"temperature_ev="<<actual.electronTemperatureEv<<'\n'
      <<"has_table="<<(table!=nullptr)<<'\n';
  if (table) {
    using namespace Moon::Species::ElectronImpact;
    for (NeutralSpecies atom : {NeutralSpecies::He,NeutralSpecies::Na,
                               NeutralSpecies::Ne,NeutralSpecies::Ar}) {
      double rate=0.0;
      if (!table->IonizationRateS(atom,actual.electronTemperatureEv,
              actual.electronDensityM3,&rate,&error)) {
        std::cerr<<error<<'\n';return 6;
      }
      std::cout<<NeutralSpeciesName(atom)<<"_rate_s="<<rate<<'\n';
    }
  }
  if (argc==3 && std::string(argv[2])=="rollback") {
    // A second successful configuration may never replace an installed run.
    auto changed=configuration;changed.electronDensityM3*=2.0;
    if (Moon::Runtime::InstallConfiguration(changed,&error) ||
        Moon::Runtime::GetConfiguration().electronDensityM3!=
          configuration.electronDensityM3) return 7;
    std::cout<<"repeat_install_rejected=1\n";
  }
  return 0;
}

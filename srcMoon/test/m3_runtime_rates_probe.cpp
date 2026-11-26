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
  configuration.hydrogenCharacteristicSpeedMS=17.5;
  configuration.hydrogenToSodiumReferenceWeightRatio=1.25;
  std::string error;
  if (!Moon::Runtime::ParseApplicationInput(argv[1],&configuration,&error)) {
    std::cerr<<error<<'\n';
    std::cout<<"parse_preserved="<<(configuration.randomSeed==-12345 &&
        configuration.hydrogenCharacteristicSpeedMS==17.5 &&
        configuration.hydrogenToSodiumReferenceWeightRatio==1.25)<<'\n';
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
  if (argc==3 && std::string(argv[2])=="rollback-m3") {
    auto invalid=configuration;
    invalid.photochemistryDataDirectory+="/missing-authoritative-source";
    if (Moon::Runtime::InstallConfiguration(invalid,&error) ||
        Moon::Runtime::HasConfiguration() ||
        Moon::Runtime::GetElectronImpactTable()!=nullptr ||
        Moon::Runtime::GetPhotochemistryRegistry()!=nullptr ||
        Moon::Runtime::GetChargeExchangeRegistry()!=nullptr) return 8;
    std::cout<<"failed_m3_install_preserved=1\n";
  }
  if (argc==3 && std::string(argv[2])=="reject-dynamic-m2") {
    // A programmatic caller must not bypass the same local/global-clock
    // incompatibility refused by the .in parser. No thermal law is invented.
    auto invalid=configuration;invalid.m2EnvironmentEnabled=true;
    if (Moon::Runtime::InstallConfiguration(invalid,&error) ||
        Moon::Runtime::HasConfiguration() ||
        error.find("m2_environment=disabled")==std::string::npos) return 10;
    std::cout<<"incompatible_m2_install_preserved=1\n";
  }
  if (argc==3 && std::string(argv[2])=="reject-incomplete-hydrogen") {
    // A caller installing a Configuration directly can also lose one of
    // these inputs. Verify that failing this check leaves the singleton free
    // for the valid, complete retry below.
    auto incomplete=configuration;
    incomplete.hydrogenToSodiumReferenceWeightRatio=0;
    if (Moon::Runtime::InstallConfiguration(incomplete,&error) ||
        Moon::Runtime::HasConfiguration() ||
        error.find("hydrogen_to_sodium_reference_weight_ratio")==std::string::npos) return 11;
    std::cout<<"incomplete_hydrogen_install_preserved=1\n";
  }
  if (argc==3 && std::string(argv[2])=="reject-incomplete-helium-ion") {
    // The He+ pair follows the same direct-install refusal as tracked H.
    auto incomplete=configuration;
    incomplete.heliumIonToSodiumReferenceWeightRatio=0;
    if (Moon::Runtime::InstallConfiguration(incomplete,&error) ||
        Moon::Runtime::HasConfiguration() ||
        error.find("helium_ion_to_sodium_reference_weight_ratio")==std::string::npos) return 12;
    std::cout<<"incomplete_helium_ion_install_preserved=1\n";
  }
  if (!Moon::Runtime::InstallConfiguration(configuration,&error)) {
    std::cerr<<error<<'\n';
    std::cout<<"install_preserved="<<(!Moon::Runtime::HasConfiguration() &&
        !Moon::Runtime::GetElectronImpactTable() &&
        !Moon::Runtime::GetPhotochemistryRegistry() &&
        !Moon::Runtime::GetChargeExchangeRegistry())<<'\n';
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
  const auto* photons=Moon::Runtime::GetPhotochemistryRegistry();
  const auto* capture=Moon::Runtime::GetChargeExchangeRegistry();
  std::cout<<"chemistry_mode="<<actual.chemistryMode<<'\n'
      <<"has_photon_registry="<<(photons!=nullptr)<<'\n'
      <<"has_capture_registry="<<(capture!=nullptr)<<'\n';
  std::cout<<"tracing_mover="<<actual.particleTracingMover<<'\n'
      <<"tracing_clock_mode="<<actual.tracingClockMode<<'\n'
      <<"tracing_maximum_substep_s="<<actual.tracingMaximumSubstepS<<'\n'
      <<"tracing_field_frame="<<actual.tracingFieldFrame<<'\n'
      <<"hydrogen_characteristic_speed_m_s="<<actual.hydrogenCharacteristicSpeedMS<<'\n'
      <<"hydrogen_to_sodium_reference_weight_ratio="<<actual.hydrogenToSodiumReferenceWeightRatio<<'\n'
      <<"helium_ion_characteristic_speed_m_s="<<actual.heliumIonCharacteristicSpeedMS<<'\n'
      <<"helium_ion_to_sodium_reference_weight_ratio="<<actual.heliumIonToSodiumReferenceWeightRatio<<'\n'
      <<"helium_to_sodium_reference_weight_ratio="<<actual.heliumToSodiumReferenceWeightRatio<<'\n';
  if (photons) {
    using Moon::Species::Photochemistry::SolarCase;
    const auto solar=actual.photochemistrySolarCase=="quiet" ? SolarCase::Quiet : SolarCase::Active;
    std::cout<<"photon_catalogue_rows="<<photons->catalogueSize()<<'\n'
        <<"photon_reactive_rows="<<photons->reactiveSize()<<'\n';
    for (const char* atom : {"He","Na","Ne","Ar"}) {
      double oneAu=0,twoAu=0;
      if (!photons->TotalDestructionRateS(atom,solar,1,1,&oneAu,&error) ||
          !photons->TotalDestructionRateS(atom,solar,2,1,&twoAu,&error)) {
        std::cerr<<error<<'\n';return 9;
      }
      std::cout<<atom<<"_photon_rate_1au_s="<<oneAu<<'\n'
          <<atom<<"_photon_rate_2au_s="<<twoAu<<'\n';
    }
  }
  if (capture) {
    using namespace Moon::Species::ChargeExchange;
    for (ChannelId id : {ChannelId::ProtonHeOrnlA32,ChannelId::ProtonNaMorgan1985LogLog}) {
      ChannelMetadata metadata;
      if (!capture->GetMetadata(id,&metadata,&error)) return 10;
      std::cout<<metadata.id<<"_daughter_kinetics_available="<<metadata.daughterKineticsAvailable<<'\n'
          <<metadata.id<<"_particle_execution_enabled="<<metadata.productionParticleExecutionEnabled<<'\n';
    }
  }
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

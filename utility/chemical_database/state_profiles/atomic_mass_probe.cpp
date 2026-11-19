// Thin original-research preparation adapter: no duplicate mass formula.
// All inputs are explicit SI central values supplied from authenticated sources.
#include "src/models/chemical_database/ReducedModels.h"
#include <iomanip>
#include <iostream>
#include <limits>
int main() {
  long double neutral,electron,energy,c;
  if(!(std::cin>>neutral>>electron>>energy>>c)) return 2;
  AMPS::ChemicalDatabase::ReducedModels::IonMass value{};
  const auto status=AMPS::ChemicalDatabase::ReducedModels::IonPhysicalMass(neutral,electron,energy,c,&value);
  if(!status.ok()) {std::cerr<<status.message<<'\n';return 2;}
  std::cout<<std::setprecision(std::numeric_limits<long double>::max_digits10)
    <<"reference_mass_kg\t"<<value.referenceMassKg<<'\n'
    <<"binding_mass_kg\t"<<value.bindingMassKg<<'\n'
    <<"physical_mass_kg\t"<<value.physicalMassKg<<'\n'
    <<"machine_epsilon\t"<<std::numeric_limits<long double>::epsilon()<<'\n';
}

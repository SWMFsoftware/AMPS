// Standalone caller example for the shared explicit-input thermal helper.
// Compile with -I/path/to/sdk/src/models/chemical_database. There is no PIC/MPI,
// database-source selection, species constant or physical parameter default.
// Every argument is supplied by the caller in the units printed by Usage().
#include "ReducedModels.h"
#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {
void Usage() {
  std::cerr << "usage: reduced-models nu_s^-1 binding_J kB_J_per_K T_K population dt_s\n"
            << "First-order constant-temperature fixed-site release only.\n"
            << "Supply an admitted, source-consistent parameter set; this executable selects none.\n";
}
bool Number(const char* text,long double* output) {
  // A malformed/overflowed/underflowed command-line token is an input error,
  // never silently changed to zero. Physical-domain guards live in the helper.
  errno=0;
  char* end=nullptr;
  const long double value=std::strtold(text,&end);
  if(end==text || *end!='\0' || errno==ERANGE) return false;
  *output=value;
  return true;
}
}

int main(int argc,char** argv) {
  if(argc!=7) {Usage();return 2;}
  long double input[6];
  for(int i=0;i<6;++i) {
    if(!Number(argv[i+1],&input[i])) {Usage();return 2;}
  }
  namespace Model=AMPS::ChemicalDatabase::ReducedModels;
  Model::ThermalState result{};
  const auto status=Model::ThermalRelease(input[0],input[1],input[2],input[3],input[4],input[5],&result);
  if(!status.ok()) {std::cerr << status.message << '\n';return 2;}
  std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10)
            << "{\"scope\":\"explicit-parameter software example; not physical qualification\","
            << "\"rate_s_inverse\":" << result.rateS
            << ",\"remaining_population\":" << result.remainingPopulation
            << ",\"released_population\":" << result.releasedPopulation
            << ",\"rate_underflowed\":" << (result.rateUnderflowed?"true":"false")
            << ",\"exposure_underflowed\":" << (result.exposureUnderflowed?"true":"false")
            << ",\"release_underflowed\":" << (result.releasedUnderflowed?"true":"false") << "}\n";
  return 0;
}

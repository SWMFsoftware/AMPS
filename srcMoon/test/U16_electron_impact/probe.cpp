// Thin standalone adapter: every coefficient/rate/unit conversion comes from
// the production service. This file contains no electron-ionization formula.
#include "ElectronImpact.h"
#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
namespace EI=Moon::Species::ElectronImpact;

bool number(const char* text,double* value) {
  char* end=nullptr;errno=0;*value=std::strtod(text,&end);
  return end!=text && *end=='\0'; // production validates finite/domain values
}

int main(int argc,char** argv) {
  std::cout << std::setprecision(17);
  if (argc==3 && std::string(argv[1])=="temperature-k") {
    double input=0,output=-123;std::string error;
    if (!number(argv[2],&input)) return 2;
    const bool ok=EI::ElectronTemperatureEvFromK(input,&output,&error);
    std::cout << "temperature_ev=" << output << '\n';
    if (!ok) std::cerr << error << '\n';
    return ok?0:4;
  }
  if (argc<3) return 2;
  const std::string mode(argv[1]);EI::MaxwellianVoronovTable table;
  std::string error;
  if (mode!="uninitialized" && !table.LoadNativeTable(argv[2],&error)) {
    std::cout << "loaded=" << table.loaded() << "\nrecord_count=" << table.recordCount() << '\n';
    std::cerr << error << '\n';return 4;
  }
  if (mode=="load" && argc==3) {
    std::cout << "loaded=" << table.loaded() << "\nrecord_count=" << table.recordCount() << '\n';
    return 0;
  }
  if (mode=="reload" && argc==4) {
    const bool reloaded=table.LoadNativeTable(argv[3],&error);
    std::cout << "reload_ok=" << reloaded << "\nloaded=" << table.loaded() << '\n';
    if (!reloaded) std::cerr << error << '\n';
    double coefficient=-123;
    const bool ok=table.RateCoefficientM3S(EI::NeutralSpecies::He,10,&coefficient,&error);
    std::cout << "coefficient_m3_s=" << coefficient << '\n';return ok?0:4;
  }
  if ((mode!="evaluate" && mode!="uninitialized" && mode!="enum-invalid") || argc!=6) return 2;
  EI::NeutralSpecies species=EI::NeutralSpecies::He;
  if (!EI::NeutralSpeciesFromName(argv[3],&species,&error)) {
    std::cout << "species_after=" << EI::NeutralSpeciesName(species) << '\n';
    std::cerr << error << '\n';return 4;
  }
  if (mode=="enum-invalid") species=static_cast<EI::NeutralSpecies>(999);
  double temperature=0,density=0;
  if (!number(argv[4],&temperature) || !number(argv[5],&density)) return 2;
  double coefficient=-123,rate=-456;EI::Coefficient row;
  const bool coefficientOk=table.RateCoefficientM3S(species,temperature,&coefficient,&error);
  const std::string coefficientError=error;
  const bool rateOk=table.IonizationRateS(species,temperature,density,&rate,&error);
  std::cout << "coefficient_m3_s=" << coefficient << "\nrate_s=" << rate << '\n';
  if (!coefficientOk || !rateOk) {
    std::cerr << coefficientError << '\n' << error << '\n';return 4;
  }
  if (!table.GetCoefficient(species,&row,&error)) return 4;
  std::cout << "atomic_number=" << row.atomicNumber << "\nelectron_count=" << row.electronCount
            << "\ntemperature_min_ev=" << row.temperatureMinEv << "\ntemperature_max_ev="
            << row.temperatureMaxEv << "\nsource_line=" << row.sourceLine << '\n';
  return 0;
}

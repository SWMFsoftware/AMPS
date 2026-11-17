// Validation-only extraction of the actual double-valued owner grid.  The
// study runner compiles this against the same value API as the native owner;
// no independent Python reconstruction is labeled as the native grid.
#include "srcSEP/util/sep_polarized_wave_state.h"
#include "srcSEP/util/sep_turbulence_core.h"
#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
double Number(const char* value) {
  char* end=nullptr;errno=0;const double result=std::strtod(value,&end);
  if(errno||end==value||*end)throw std::invalid_argument("invalid explicit grid number");
  return result; // ResolveCenteredLogGrid validates finite positive bounds.
}
}
int main(int argc,char** argv) {
  try {
    SEP::Turbulence::Configuration configuration;
    if(argc!=1&&argc!=4)throw std::invalid_argument("supply either no grid overrides, or N first_center last_center");
    if(argc==4) {
      char* end=nullptr;errno=0;const auto count=std::strtoull(argv[1],&end,10);
      if(errno||argv[1][0]=='-'||end==argv[1]||*end||count>std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("invalid explicit bin count");
      configuration.spectralBins=static_cast<std::size_t>(count);
      configuration.spectralKMinPerM=Number(argv[2]);configuration.spectralKMaxPerM=Number(argv[3]);
    }
    const auto result=SEP::Turbulence::Polarized::ResolveCenteredLogGrid(configuration.spectralBins,
        configuration.spectralKMinPerM,configuration.spectralKMaxPerM);
    if(!result.status.ok())throw std::invalid_argument(result.status.message);
    const auto& grid=result.grid;
    std::cout<<"{\"schema\":\"srcsep-native-centered-grid-v1\",\"bins\":"<<grid.centersPerM.size()
        <<",\"configured_first_center_per_m\":\""<<std::setprecision(17)<<configuration.spectralKMinPerM
        <<"\",\"configured_last_center_per_m\":\""<<configuration.spectralKMaxPerM
        <<"\",\"log_width_hex\":\""<<std::hexfloat<<grid.logWidth<<"\",\"centers_hex\":[";
    for(std::size_t j=0;j<grid.centersPerM.size();++j) {
      if(j)std::cout<<',';
      std::cout<<'\"'<<grid.centersPerM[j]<<'\"';
    }
    std::cout<<"],\"edges_hex\":[";
    for(std::size_t j=0;j<grid.edgesPerM.size();++j) {
      if(j)std::cout<<',';
      std::cout<<'\"'<<grid.edgesPerM[j]<<'\"';
    }
    std::cout<<"]}\n";
    if(!std::cout)throw std::runtime_error("cannot write native grid image");
    return 0;
  } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}

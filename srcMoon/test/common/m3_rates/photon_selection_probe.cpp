#include "Photochemistry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

// This adapter supplies only reproducible uniform draws and counts the IDs
// returned by the real immutable production registry. Source probabilities,
// independent binomial bounds and scientific scoring belong to Python; no
// duplicate reaction probabilities or acceptance law is compiled here.
namespace Photo=Moon::Species::Photochemistry;

int main(int argc,char** argv) {
  if (argc!=2) {
    std::cerr<<"Usage: photon_selection_probe <frozen-D04-processed-directory>\n";
    return 2;
  }
  Photo::Registry registry;
  std::string error;
  if (!registry.LoadFrozenD04(Photo::FrozenV1Sources(argv[1]),&error)) {
    std::cerr<<"D04 production registry refused source: "<<error<<'\n';
    return 2;
  }
  constexpr std::uint32_t masterSeed=731993u;
  constexpr std::uint64_t draws=100000u;
  const std::array<std::string,6> parents={{"He","Na","Ne","Ar","H2O","OH_exp"}};
  std::uint32_t caseIndex=0;
  std::cout<<std::setprecision(17);
  for (const std::string& parent:parents) for (int solar=0;solar<2;++solar,++caseIndex) {
    // Each case starts from its own explicitly declared two-word seed vector.
    // This never shares a draw with a previous case or selects a favorable seed.
    std::seed_seq seed{masterSeed,caseIndex};
    std::mt19937_64 generator(seed);
    std::map<std::string,std::uint64_t> counts;
    std::vector<Photo::Channel> channels;
    if (!registry.ChannelsForParent(parent,&channels,&error)) {
      std::cerr<<"Production registry refused family metadata: "<<error<<'\n';
      return 2;
    }
    // An unobserved rare outcome still has an explicit measured count of0.
    // Omitting it would confuse a valid count with a missing native record.
    for (const auto& channel:channels) counts[channel.id]=0;
    double minimum=1.0,maximum=0.0;
    for (std::uint64_t draw=0;draw<draws;++draw) {
      const std::uint64_t bits=generator()>>12;
      // The integer has52bits. Adding1/2 requires at most53significand bits,
      // so even the largest midpoint is exactly below1 in binary64. A53bit
      // integer plus1/2 could instead round its maximum upward to1.
      const double uniform=std::ldexp(static_cast<double>(bits)+0.5,-52);
      if (!(uniform>0.0 && uniform<1.0)) {
        std::cerr<<"Open-uniform adapter violated its declared support\n";
        return 2;
      }
      minimum=std::min(minimum,uniform);maximum=std::max(maximum,uniform);
      std::string channel;
      if (!registry.SelectChannel(parent,static_cast<Photo::SolarCase>(solar),uniform,&channel,&error)) {
        std::cerr<<"Production selector refused case "<<caseIndex<<": "<<error<<'\n';
        return 2;
      }
      ++counts[channel];
    }
    // Source IDs are strict frozen ASCII identifiers; none contains JSON
    // quotes/escapes. Python independently checks the exact source ID set.
    std::cout<<"{\"case_index\":"<<caseIndex<<",\"parent\":\""<<parent
      <<"\",\"solar\":"<<solar<<",\"draws\":"<<draws
      <<",\"seed_words\":["<<masterSeed<<','<<caseIndex<<"]"
      <<",\"uniform_min\":"<<minimum<<",\"uniform_max\":"<<maximum<<",\"counts\":{";
    bool first=true;
    for (const auto& entry:counts) {
      if (!first) std::cout<<',';
      first=false;
      std::cout<<'"'<<entry.first<<"\":"<<entry.second;
    }
    std::cout<<"}}\n";
  }
  return 0;
}

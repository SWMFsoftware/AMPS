// Thin adapter: every digest/uniform comes from the maintained production RNG.
#include "../../species/QualificationRandom.h"
#include <iomanip>
#include <iostream>
#include <sstream>
namespace R=Moon::Species::QualificationRandom;
int main() {
  std::string line;
  while(std::getline(std::cin,line)) {
    std::istringstream in(line);R::Key key;std::uint64_t purpose,index,count;
    if(!(in>>key.campaignSeed>>key.originalParent>>purpose>>index>>std::quoted(key.cohort)>>count))return 2;
    for(std::uint64_t i=0;i<count;++i) {std::uint64_t value;if(!(in>>value))return 2;key.genealogy.push_back(value);}
    std::string digest="UNCHANGED",error;double uniform=17.125;
    const bool a=R::DigestHex(key,static_cast<R::Purpose>(purpose),index,&digest,&error);
    const bool b=R::OpenUniform(key,static_cast<R::Purpose>(purpose),index,&uniform,&error);
    std::cout<<a<<' '<<b<<' '<<std::quoted(digest)<<' '<<std::setprecision(17)<<uniform<<' '<<std::quoted(error)<<'\n';
  }
  return 0;
}

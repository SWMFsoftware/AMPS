#ifndef SRCMOON_TEST_EXACT_ZERO_FORCE_RECIPE_H
#define SRCMOON_TEST_EXACT_ZERO_FORCE_RECIPE_H
#include <algorithm>
#include <cctype>
#include <string>
namespace Moon { namespace SodiumCohortInput {
// Admit ONLY the complete existing core zero-acceleration call. A marker,
// bare type/function name, additional statement or suffix is not proof.
inline bool ExactZeroAccelerationMacro(std::string text){
 text.erase(std::remove_if(text.begin(),text.end(),[](unsigned char c){return std::isspace(c)!=0;}),text.end());
 if(!text.empty()&&text.back()==';')text.pop_back();
 return text=="TotalParticleAcceleration_default(a,s,p,x,v,n)"||text=="PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n)";
}
}}
#endif

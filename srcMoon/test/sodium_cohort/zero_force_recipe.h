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
// Admit ONLY the complete gravity-only reference call (I21 decision 1,
// 2026-11-26): the production lunar point-mass term selected by the
// lunar-point-mass-gravity compile recipe.
inline bool ExactLunarGravityReferenceMacro(std::string text){
 text.erase(std::remove_if(text.begin(),text.end(),[](unsigned char c){return std::isspace(c)!=0;}),text.end());
 if(!text.empty()&&text.back()==';')text.pop_back();
 return text=="Moon::LunarPointMassReferenceAcceleration(a,s,p,x,v,n)";
}
}}
#endif

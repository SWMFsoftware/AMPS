#ifndef SRCMOON_NATIVE_REACTION_REGISTRATION_H
#define SRCMOON_NATIVE_REACTION_REGISTRATION_H
#include <string>
namespace Moon { namespace Chemistry { namespace NativeRegistration {
// A caller owns the immutable prepared source/profile/controller context until
// shutdown. It registers before amps_init; the normal application invokes it
// after normal SPICE/mesh/species/weights setup, BEFORE particle allocation.
using Callback=bool(*)(void*,std::string*);
bool SetBeforeParticleBuffer(Callback,void*,const std::string& identity,std::string* error);
bool RunBeforeParticleBuffer(std::string* error);
const std::string& InstalledIdentity();
}}}
#endif

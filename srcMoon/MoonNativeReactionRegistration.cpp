#include "MoonNativeReactionRegistration.h"
#include "pic.h"
#include <exception>
namespace Moon { namespace Chemistry { namespace NativeRegistration {
namespace {Callback callback=nullptr;void* user=nullptr;std::string identity;bool called=false;}
bool SetBeforeParticleBuffer(Callback value,void* context,const std::string& id,std::string* error){
  if(!error)return false;
  if(!value||id.empty()||callback||called||PIC::ParticleBuffer::ParticleDataBuffer){
    *error="normal Moon chemistry registration requires one named callback, installed before amps_init/particle allocation";return false;}
  // Construct the fallible identity before publishing process state.
  std::string copied=id;identity.swap(copied);user=context;callback=value;error->clear();return true;
}
bool RunBeforeParticleBuffer(std::string* error){
  if(!error)return false;
  if(!callback){error->clear();return true;} // Unselected legacy path is exact.
  if(called||PIC::ParticleBuffer::ParticleDataBuffer){
    *error="normal Moon chemistry registration was called twice or after particle allocation";return false;}
  // A callback can register opaque storage/hooks; after failure these cannot be
  // unregistered safely. The caller MUST abort startup rather than retry it.
  called=true;
  try {if(!callback(user,error)){if(error->empty())*error="named Moon registration callback refused without a diagnostic";return false;}}
  catch(const std::bad_alloc&){*error="normal Moon source/profile registration exhausted memory before particle allocation. Increase available memory or reduce the explicit cohort/state inventory/genealogy capacity, then restart; no physical event was published";return false;}
  catch(const std::exception& e){*error="normal Moon pre-buffer registration exception: "+std::string(e.what());return false;}
  if(PIC::ParticleBuffer::ParticleDataBuffer){*error="normal Moon registration callback illegally allocated the particle buffer; abort startup before a second initialization";return false;}
  error->clear();return true;
}
const std::string& InstalledIdentity(){return identity;}
}}}

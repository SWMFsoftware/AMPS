#ifndef SRCSEP_UTIL_SEP_BIRTH_NATIVE_VALIDATION_H
#define SRCSEP_UTIL_SEP_BIRTH_NATIVE_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
namespace SEP { namespace Testing {
std::vector<Descriptor> BirthNativeDescriptors();
// A selected native lifecycle fixture supplies its scientific startup record
// explicitly.  The three operator roles are separate freshly initialized
// processes with the same initialized AMPS timestep. Direct controls query
// the existing explicit-duration diagnostic operator without physical sampling.  BIRTHN06 instead
// requests the existing, genuinely supported reduced-material checkpoint host.
// Missing input is reported by the dispatcher; malformed selected protocols
// throw before mesh allocation rather than guessing a physical background.
struct BirthNativeRequest {
  bool selected=false,material=false;
  std::uint32_t campaignReplicateIndex=0;
  std::string inputPath,eventPath,role,prefix;
};
BirthNativeRequest ReadBirthNativeRequest(
    const std::string& testId,const std::string& requestPath);
} }
#endif

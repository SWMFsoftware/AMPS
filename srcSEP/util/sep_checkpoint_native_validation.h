#ifndef SRCSEP_UTIL_SEP_CHECKPOINT_NATIVE_VALIDATION_H
#define SRCSEP_UTIL_SEP_CHECKPOINT_NATIVE_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
#include <vector>
namespace SEP { namespace Testing {
std::vector<Descriptor> CheckpointNativeDescriptors();
std::vector<Descriptor> FreshProcessCheckpointDescriptors();
// Dedicated frozen-operator executable; no particle motion on fixed geometry.
std::vector<Descriptor> OperatorCheckpointDescriptors();
std::vector<Descriptor> OperatorFreshCheckpointDescriptors();
// Read the explicit startup deck from a selected fresh-process role protocol.
// Unrelated tests return empty; malformed selected requests throw before model
// allocation. This consumes test input, not environment or production argv.
struct MaterialCheckpointRequest {
  bool selected=false;
  std::string inputPath,eventPath,role,prefix;
};
MaterialCheckpointRequest ReadMaterialCheckpointRequest(
    const std::string& testId,const std::string& requestPath);
} }
#endif

#ifndef SRCSEP_UTIL_SEP_MATERIAL_TUBE_KINEMATICS_VALIDATION_H
#define SRCSEP_UTIL_SEP_MATERIAL_TUBE_KINEMATICS_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
// Shared native/dependency-light callbacks use explicit SI material maps and
// independent analytic geometry references. They never publish a native grid.
namespace SEP { namespace Testing {std::vector<Descriptor> MaterialTubeKinematicsDescriptors();} }
#endif

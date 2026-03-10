#ifndef _SRC_EARTH_3D_MODE3D_H_
#define _SRC_EARTH_3D_MODE3D_H_

#include "../util/amps_param_parser.h"

namespace Earth {
namespace Mode3D {

extern bool ParsedDomainActive;
extern double ParsedDomainMin[3];
extern double ParsedDomainMax[3];

int Run(const EarthUtil::AmpsParam& prm);

} // namespace Mode3D
} // namespace Earth

#endif

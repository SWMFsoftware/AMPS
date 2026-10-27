#ifndef _SRC_GEOSPACE_3D_MODE3D_H_
#define _SRC_GEOSPACE_3D_MODE3D_H_

#include "../util/amps_param_parser.h"

namespace Geospace {
namespace Mode3D {

extern bool ParsedDomainActive;
extern double ParsedDomainMin[3];
extern double ParsedDomainMax[3];

// Optional user-defined AMR resolution profile for standalone Mode3D.
// If inactive, main_lib.cpp keeps using its historical hard-coded resolution
// function.  ConfigureMeshResolutionProfile() is called by Mode3D::Run() before
// amps_init_mesh(), and localResolution()/localSphericalSurfaceResolution() query
// ConfiguredMeshResolutionSI() while the AMR tree is being built.
extern bool MeshResolutionProfileActive;
extern double MeshResolutionEarth_m;
extern double MeshResolutionBoundary_m;
extern double MeshResolutionOuterRadius_Re;
extern double MeshResolutionExponent;
extern int MeshResolutionCoarseningCode;

void ConfigureMeshResolutionProfile(const GeospaceUtil::AmpsParam& prm);
double ConfiguredMeshResolutionSI(const double *x);

int Run(const GeospaceUtil::AmpsParam& prm);

} // namespace Mode3D
} // namespace Geospace

#endif

#ifndef SEP3D_CORE_DOMAIN_GEOMETRY_H
#define SEP3D_CORE_DOMAIN_GEOMETRY_H

#include "parker_geometry.h"

namespace SEP3D { namespace Core {

// Pure geometry inputs: no AMPS, MPI, parser or runtime ownership lives here.
struct DomainGeometryParameters {
  bool cornerCube = false;
  // x-y corner mode retains heliocentric coordinates and makes the two z
  // faces equidistant from the Sun. The cube may grow to contain a tilted
  // field line; no corridor or solar-neighbourhood volume is clipped.
  bool centerZ = false;
  Vec3 originM;
  Vec3 cornerDirection;  // -1/+1/0 automatic; z must be 0 when centerZ
  double endpointRadiusM = Const::AU;
  double solarSphereRadiusM = 0.0;
  double corridorPaddingM = 0.0;
  double cornerMarginM = 0.0;
  ParkerSpiralGeometry parker;
};

Status BuildDomainBoundsM(const DomainGeometryParameters& parameters,
                          Vec3* minimumM, Vec3* maximumM);

} }  // namespace SEP3D::Core
#endif

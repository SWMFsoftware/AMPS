#ifndef _GENERIC_COMET_MODEL_H_
#define _GENERIC_COMET_MODEL_H_

#include <string>

namespace Comet {
namespace Model {

namespace Geometry {
// Geometry selects the backend that interprets the same species/source data.
// Spherical is implemented by RadialMain, Triangulated by the legacy 3D main,
// and Axisymmetric is reserved for the planned DIM2 backend.
enum Type {
  Spherical=0,
  Axisymmetric=1,
  Triangulated=2
};
}

namespace SurfaceDistribution {
// These identifiers describe where direct nucleus gas is released.  Generic
// input records store the identifier per species; a geometry backend maps it
// to an area measure and an injection location.
enum Type {
  Disabled=0,
  Uniform=1,
  Cosine=2,
  Tabulated=3,
  LegacyBjorn=4,
  LegacyJet=5
};
}

namespace ExtendedDistribution {
// Extended sources are represented separately from surface production so a
// case can enable or omit grains, distributed parents, or a custom source
// without changing the common comet state.
enum Type {
  Disabled=0,
  RadialProfile=1,
  GrainSublimation=2,
  UserDefined=3
};
}

extern int GeometryMode;
extern double HeliocentricDistanceAU;
extern char SurfaceMeshFileName[];
extern double SurfaceMeshLengthUnit;

// Every array is indexed by the AMPS species number. Rates are molecules/s,
// temperatures are K, and speeds are m/s. A negative temperature or speed
// means the input did not specify that quantity.
extern double SurfaceProductionRate[];
extern int SurfaceDistributionMode[];
extern double SurfaceTemperature[];
extern double SurfaceBulkSpeed[];
extern double ExtendedProductionRate[];
extern int ExtendedDistributionMode[];

// Validate always checks geometry, rates, and distribution identifiers.  A
// backend passes requireRunnableSources=true when it will inject particles;
// that additionally requires explicit temperatures and bulk speeds.
bool Validate(std::string *error=NULL,bool requireRunnableSources=false);
const char *GeometryName();

}
}

#endif

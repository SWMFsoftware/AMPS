#ifndef SEP_CORONA_SWCME_SPHERICAL_SECTOR_H
#define SEP_CORONA_SWCME_SPHERICAL_SECTOR_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <array>

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

// Static curvilinear control-volume measures for the initial cross-ray mesh
// increment. mu=sin(latitude)=cos(colatitude); phi increases right-handedly
// about basis.z. Radii are Sun-centred HCI metres and all angles are radians.
// These measures do not yet provide a moving piston mesh or a CT update.
struct SectorBounds {
  double innerRadiusM=0,outerRadiusM=0;
  double muLower=0,muUpper=0;
  double phiLowerRad=0,phiUpperRad=0;
};

struct SectorBasis {
  // Explicit right-handed orthonormal HCI axes. No input frame is silently
  // normalised/rotated to make inconsistent geometry look valid.
  CoronalCME::Vec3 x{1,0,0},y{0,1,0},z{0,0,1};
};

enum class SectorFace : int { Inner=0,Outer=1,MuLower=2,MuUpper=3,
                              PhiLower=4,PhiUpper=5 };

struct SectorFaceMeasure {
  // Integral n_out dA in HCI m^2 versus integral dA in m^2. Their magnitudes
  // differ on curved faces. Zero polar faces are exact geometric degeneracies,
  // never finite-area patches dropped by a numerical tolerance.
  CoronalCME::Vec3 outwardAreaVectorM2;
  double physicalAreaM2=0;
};

struct SectorMeasure {
  double solidAngleSr=0;
  double volumeM3=0;
  // Exact first moment/volume. A broad nonconvex sector may have its centroid
  // outside the physical cell; this is not automatically a plasma sample site.
  CoronalCME::Vec3 volumeCentroidHciM;
  std::array<SectorFaceMeasure,6> faces;
  // Raw discrete geometric residual, not corrected by modifying a face.
  CoronalCME::Vec3 areaClosureResidualM2;
};

// Exact analytic angular/radial moments of a spherical sector. Requires
// 0<=rL<rR, -1<=muL<muR<=1 and 0<phiR-phiL<=2pi. Shared angular edges must
// be supplied with identical coordinates by the future mesh authority.
// Invalid configuration and nonrepresentable measures reject with typed
// status; no volume/area/angle floor or inferred legacy-ray boundary is used.
Core::Result<SectorMeasure> MeasureSphericalSector(
    const SectorBounds& bounds,const SectorBasis& basis=SectorBasis{});

} } } // namespace SEP::CoronaSwcme::CrossRay

#endif

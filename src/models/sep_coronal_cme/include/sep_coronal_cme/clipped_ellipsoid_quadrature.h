#ifndef SEP_CORONAL_CME_CLIPPED_ELLIPSOID_QUADRATURE_H
#define SEP_CORONAL_CME_CLIPPED_ELLIPSOID_QUADRATURE_H

#include "sep_coronal_cme/ellipsoid_geometry.h"

namespace SEP { namespace CoronalCME {

// Integrate only the part of each ORIGINAL parameter cell outside the physical
// solar sphere. Unlike the legacy midpoint tessellator, a cell whose midpoint
// is buried can retain a positive exposed cap/sliver. The returned area [m^2]
// is the integral on that exposed part; its representative position and normal
// are canonical ellipsoid evaluations at an actually exposed point in the SAME
// cell. It is a quadrature representative, not the area centroid or a claim that
// every point of a partly cut cell has the representative's plasma coverage.
//
// IDs retain 1 + polarIndex*azimuthCells + azimuthIndex; cells with zero exposed
// area leave holes. Solar clipping, source qualification, and plasma coverage
// remain different operations. The caller must still certify plasma coverage
// on the complete physical patch, not merely at this exposed representative.
//
// Numerical contract: exact quadratic cuts in cos(theta), explicit azimuthal
// topology splits, and positive Gauss8/16 area quadrature with relative targets
// 2e-13 (inner) and 2e-11 (outer). Recursion is bounded by 20 subdivisions and
// the complete call by 50 million integrand evaluations / 1,048,576 original
// cells. Exhaustion or unrepresentable positive samples returns an error; no
// unresolved sliver is silently discarded. The AMPS physical dt is unrelated.
Core::Result<std::vector<SurfacePatch>> TessellateSolarExposed(
    const FixedOrientationEllipsoid&,double solarRadiusM,
    int polarCells,int azimuthCells);

} }
#endif

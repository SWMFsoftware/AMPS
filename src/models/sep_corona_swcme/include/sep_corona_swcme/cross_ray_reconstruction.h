#ifndef SEP_CORONA_SWCME_CROSS_RAY_RECONSTRUCTION_H
#define SEP_CORONA_SWCME_CROSS_RAY_RECONSTRUCTION_H

#include "sep_corona_swcme/cross_ray_ct.h"

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

struct ReconstructedMagneticCell {
  CoronalCME::Vec3 fieldHciT;
  // Residual of the piecewise-constant fit, in T, on the actual curved faces.
  // A nonuniform field generally has a nonzero fit residual even with exact
  // div B=0. This is retained, not mistaken for or subtracted as a monopole.
  double maximumNormalFitResidualT=0;
  // Ratio of largest to smallest Cholesky pivot of the dimensionless normal
  // equations; a conditioning diagnostic, not a physical magnetic tolerance.
  double pivotRatio=0;
};

// First-order constant reconstruction from the authoritative CT face flux:
// minimize sum_f [(S_f/A_f).B - Phi_f/A_f]^2. S_f=integral n dA in HCI m^2,
// A_f=integral dA in m^2, and Phi_f in T m^2. Using the physical face area
// rather than |S_f| avoids pretending a curved face is a plane. The result
// is exact for constant Cartesian B on a resolvable closed cell, but is not
// a face-field projection, magnetic-energy correction or higher-order fit.
//
// The 3x3 symmetric normal system is formed and solved in long double. Its
// pivots must exceed an explicit floating-point conditioning bound; no
// regularization, arbitrary axis selection or field floor fills a deficient
// cell. The input face flux, its raw divergence and its identity stay intact.
Core::Result<std::vector<ReconstructedMagneticCell>> ReconstructCellField(
    const StaticMagneticFlux& field);

} } } // namespace SEP::CoronaSwcme::CrossRay

#endif

#ifndef SRCSEP_REDUCED_SOURCE_H
#define SRCSEP_REDUCED_SOURCE_H
#include "sep_reduced_geometry.h"
#include "sep_reduced_particle_options.h"
#include <vector>

namespace SEP { namespace ReducedParticles {
// A particle is located by its segment/fraction at its own birth epoch. The
// fraction labels the same material chord in the frozen departure mesh; xM
// is its physical position on that chord at birth, not an off-line normal
// displacement of the shock intersection. No native storage is mutated.
struct ReleasePoint {
  std::size_t segment=0;
  double fraction=0;
  ReducedShock::Position xM{},normal{};
  ReducedShock::FrontSample normalFoot;
};
Core::Result<std::vector<ReleasePoint>> ReferenceIntersections(
    const ReducedShock::IntervalView&,const std::vector<ReducedShock::Position>&
    verticesAtBirth,double birthEpochS,double normalDistanceM,double toleranceM);

// eta is an explicitly prescribed NET first-passage number fraction. Phi is
// the conserved tube flux [T m^2], so A=Phi/|B| and A/|b.n| is the shock patch
// intersecting that tube. This returns a rate [physical protons/s], never an
// entire-front rate, macro count, electron flux, or gross-to-net conversion.
Core::Result<double> NetProtonRate(const ReducedShock::FrontSample&,
    double conservedTubeFluxTm2,double netFraction);
} }
#endif

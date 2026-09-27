// ============================================================================
// AMPS-independent relativistic population-control kernel.
//
// AMPS' generic 3-to-2 merger conserves sum(w |v|^2), which is the
// nonrelativistic kinetic-energy moment.  SEP particles may have relativistic
// energies, so srcSEP3D uses this exact momentum-space reconstruction instead.
// The AMPS adapter owns linked-list mutations and metadata; this file owns only
// the reviewable conservation algebra and is directly unit tested.
// ============================================================================

#ifndef SEP3D_TRANSPORT_POPULATION_CONTROL_H
#define SEP3D_TRANSPORT_POPULATION_CONTROL_H

#include "../core/sep3d_types.h"

#include <array>

namespace SEP3D {
namespace Transport {

struct WeightedPhasePoint {
  double weight = 0.0;
  Core::Vec3 positionM;
  Core::Vec3 momentumKgMPerS;
};

struct RelativisticMergeResult {
  Core::Status status;
  double outputWeight = 0.0;
  Core::Vec3 outputPositionM;
  Core::Vec3 firstMomentumKgMPerS;
  Core::Vec3 secondMomentumKgMPerS;
  double relativeWeightResidual = 0.0;
  double relativeMomentumResidual = 0.0;
  double relativeEnergyResidual = 0.0;
};

double RelativisticTotalEnergyJ(const Core::Vec3& momentumKgMPerS,
                                double massKg);

// Merge three weighted phase points into two equal-weight points.  Both output
// positions are the input weighted centroid.  For a caller-provided unit
// direction n, p_A=p_bar+q*n and p_B=p_bar-q*n conserve momentum identically;
// q is the unique non-negative root that also conserves relativistic total
// energy (and therefore kinetic energy because represented rest mass is fixed).
RelativisticMergeResult MergeRelativisticThreeToTwo(
    const std::array<WeightedPhasePoint, 3>& input,
    double massKg, const Core::Vec3& direction);

}  // namespace Transport
}  // namespace SEP3D

#endif  // SEP3D_TRANSPORT_POPULATION_CONTROL_H

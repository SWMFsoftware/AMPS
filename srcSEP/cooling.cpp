#include <cmath>
#include <limits>

#include "cooling.h"

/*
 * Namespace: COOLING_FACTOR_PARKER
 * This namespace contains functions for calculating the adiabatic cooling factor
 * for Solar Energetic Particles (SEPs) in the heliosphere, based on the Parker
 * transport equation. The adiabatic cooling accounts for the energy loss of SEPs
 * due to the expansion of the solar wind.
 *
 * The radial distance from the Sun (r) is the heliocentric distance in meters.
 *
 * References:
 * - Parker, E. N. (1965). "The passage of energetic charged particles through interplanetary space."
 *   Planetary and Space Science, 13(1), 9-49.
 * - Schlickeiser, R. (2002). *Cosmic Ray Astrophysics*. Springer.
 * - Zhang, M. (2000). "Adiabatic cooling of solar energetic particles."
 *   The Astrophysical Journal, 541(2), 428-436.
 */

namespace COOLING_FACTOR_PARKER {

    SEP::Transport::ScalarResult CalculateAdiabaticMomentum(
        double p_old, double r, double dt, double V_sw) {
        SEP::Transport::ScalarResult result;
        if (!std::isfinite(r) || !std::isfinite(V_sw) || r <= 0.0 ||
            V_sw < 0.0) {
            result.status = SEP::Transport::Status::Error(
                SEP::Transport::StatusCode::InvalidArgument,
                "radial cooling requires r>0 and finite V_sw>=0");
            return result;
        }

        // No artificial inner-radius floor is applied.  The caller must supply
        // a valid physical domain, so a geometry error cannot silently change
        // the represented cooling rate.
        const double divergencePerS = 2.0 * V_sw / r;
        return SEP::Transport::ApplyAdiabaticMomentum(
            p_old, divergencePerS, dt);
    }

    double calculateAdiabaticCooling(double p_old, double r, double dt,
                                     double V_sw) {
        const SEP::Transport::ScalarResult result =
            CalculateAdiabaticMomentum(p_old, r, dt, V_sw);
        return result.status.ok()
            ? result.value
            : std::numeric_limits<double>::quiet_NaN();
    }

} // namespace COOLING_FACTOR_PARKER


#ifndef _COOLING_FACTOR_PARKER_
#define _COOLING_FACTOR_PARKER_

#include "util/sep_transport_common.h"

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

    // Physical constants
    const double c = 2.99792458e8;      // Speed of light in m/s
    const double e = 1.602176634e-19;   // Elementary charge in C
    const double mp = 1.67262192369e-27; // Proton mass in kg
    const double AU = 1.495978707e11;   // Astronomical Unit in meters

    // Status-returning radial-wind adapter for the common plasma-frame
    // adiabatic update.  The radial model supplies div(U)=2*V_sw/r; the common
    // kernel then integrates dp/dt=-(p/3)div(U) exactly over dt.
    SEP::Transport::ScalarResult CalculateAdiabaticMomentum(
        double p_old, double r, double dt, double V_sw);

    // Compatibility wrapper retained for external callers of the historical
    // API.  It returns quiet NaN for invalid input; new code must call the
    // status-returning function above and handle the reported cause explicitly.
    extern double calculateAdiabaticCooling(double p_old, double r, double dt, double V_sw);
} // namespace COOLING_FACTOR_PARKER

#endif

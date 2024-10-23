#include <iostream>
#include <cmath>

#ifndef _QLT1_
#define _QLT1_

/*
 * Namespace: QLTApproximation
 * This namespace contains functions based on the Quasi-Linear Theory (QLT) approximation
 * for calculating the mean free path of Solar Energetic Particle (SEP) protons in the heliosphere,
 * including the effect of heliocentric distance on the turbulence energy density.
 *
 * References:
 * - Schlickeiser, R. (2002). *Cosmic Ray Astrophysics*. Springer.
 * - Zank, G. P., & Matthaeus, W. H. (1992). "The radial evolution of interplanetary turbulence."
 *   *Journal of Geophysical Research*, 97(A11), 17189–17200.
 * - Bruno, R., & Carbone, V. (2013). "The solar wind as a turbulence laboratory."
 *   *Living Reviews in Solar Physics*, 10(1), 2.
 */
namespace QLT1 {

    // Physical constants
    const double c = 2.99792458e8;          // Speed of light in m/s
    const double e = 1.602176634e-19;       // Elementary charge in C
    const double mp = 1.67262192369e-27;    // Proton mass in kg
    const double AU = 1.495978707e11;       // Astronomical Unit in meters
    const double B1AU = 5e-9;               // Magnetic field strength at 1 AU in Tesla
    const double deltaB2_1AU = 1e-10;       // Turbulence energy density at 1 AU in (T^2)
    const double s = 3.5;                   // Turbulence radial decay exponent (3 ≤ s ≤ 4)
    const double q = 5.0 / 3.0;             // Spectral index for Kolmogorov turbulence

    /*
     * Function: calculateMeanFreePath
     * Calculates the mean free path (λₚₐᵣₐₗₗₑₗ) of SEP protons using the Quasi-Linear Theory approximation,
     * including the effect of heliocentric distance on the turbulence energy density.
     *
     * Inputs:
     * - v: Particle speed in m/s (must be less than the speed of light)
     * - r: Heliocentric distance in meters 
     *
     * Output:
     * - Mean free path in meters
     *
     * Assumptions and Approximations:
     * - The interplanetary magnetic field (IMF) follows the Parker spiral model.
     * - The mean magnetic field strength decreases with distance as B ∝ r⁻².
     * - The turbulence energy density decreases with distance as δB² ∝ r⁻ˢ.
     * - The turbulence level δB²/B₀² varies with distance accordingly.
     * - Magnetic turbulence follows a Kolmogorov spectrum with spectral index q = 5/3.
     * - Uses an average pitch-angle cosine μ to simplify integration over pitch angles.
     * - Neglects perpendicular diffusion and nonlinear effects.
     *
     * References:
     * - Schlickeiser, R. (2002). *Cosmic Ray Astrophysics*. Springer.
     * - Zank, G. P., & Matthaeus, W. H. (1992). "The radial evolution of interplanetary turbulence."
     *   *Journal of Geophysical Research*, 97(A11), 17189–17200.
     */
    extern double calculateMeanFreePath(double v, double r);
} // namespace QLTApproximation


#endif

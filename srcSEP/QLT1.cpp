#include <iostream>
#include <cmath>

#include "QLT1.h"

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
namespace QLTApproximation {
    /*
     * Function: calculateMeanFreePath
     * Calculates the mean free path (λₚₐᵣₐₗₗₑₗ) of SEP protons using the Quasi-Linear Theory approximation,
     * including the effect of heliocentric distance on the turbulence energy density.
     *
     * Inputs:
     * - v: Particle speed in m/s (must be less than the speed of light)
     * - r: Heliocentric distance in Astronomical Units (AU)
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
    double calculateMeanFreePath(double v, double r_meters) {
        // Validate inputs
        if (v <= 0 || v >= c) {
            std::cerr << "Error: Particle speed must be between 0 and the speed of light." << std::endl;
            return -1.0;
        }
        if (r_meters <= 0) {
            std::cerr << "Error: Heliocentric distance must be positive." << std::endl;
            return -1.0;
        }

        // Mean magnetic field strength at distance r (B ∝ r⁻²)
        double B = B1AU * pow(AU / r_meters, 2); // Tesla

        // Turbulence energy density at distance r (δB² ∝ r⁻ˢ)
        double deltaB2 = deltaB2_1AU * pow(AU / r_meters, s); // T^2

        // Turbulence level δB² / B₀² at distance r
        double deltaB2_over_B0_2 = deltaB2 / (B * B);

        // Calculate Lorentz factor γ
        double beta = v / c;
        double gamma = 1.0 / sqrt(1 - beta * beta);

        // Proton momentum p = γ * m_p * v
        double p = gamma * mp * v; // kg·m/s

        // Proton rigidity R = p / (Z * e), Z = 1 for protons
        double R = p / e; // V·s (Volt·seconds)

        // Proton gyrofrequency Ω = (Z * e * B) / (γ * m_p)
        double Omega = e * B / (gamma * mp); // rad/s

        // Average pitch-angle cosine μ (approximation)
        double mu = 0.7; // Dimensionless
        double mu2 = mu * mu;

        // Resonant wave number k_res = Ω / (v * μ)
        double k_res = Omega / (v * fabs(mu)); // 1/m

        // Reference wave number k0 corresponding to the turbulence outer scale (~1e-10 m⁻¹)
        double k0 = 1e-10; // 1/m

        // Turbulence power spectrum at k_res: [δB²(k_res) / B₀²] ∝ (k_res / k0)⁻q
        double deltaB2_over_B0_2_kres = deltaB2_over_B0_2 * pow(k_res / k0, -q);

        // Pitch-angle diffusion coefficient D_μμ at μ
        // D_μμ = (π / 2) * Ω² * (1 - μ²) * [δB²(k_res) / B₀²] / (v * |μ|)
        double D_mu_mu = (M_PI / 2.0) * Omega * Omega * (1.0 - mu2) * deltaB2_over_B0_2_kres / (v * fabs(mu)); // 1/s

        // Mean free path λ_parallel = (3 * v) / (8 * D_μμ)
        double lambda_parallel = (3.0 * v) / (8.0 * D_mu_mu); // meters

        return lambda_parallel;
    }

} // namespace QLTApproximation



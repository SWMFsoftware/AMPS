#include <iostream>
#include <cmath>
#include "QLT.h"

namespace QLT {
    // Function to calculate proton gyrofrequency Omega (rad/s)
    //  \[ \Omega = \frac{q B}{m_p} \]
    // Where \( q \) is the proton charge, \( B \) is the magnetic field, and \( m_p \) is the proton mass.
    // The gyrofrequency \( \Omega \) is used to calculate the resonant wavenumber for particle interactions with magnetic turbulence.

    double calculateOmega(double B) {
        return proton_charge * B / proton_mass;
    }

    // Resonant Wavenumber Calculation (calculateKParallel)
    // The resonant wavenumber is given by:
    // \[ k_{\parallel} = \frac{\Omega}{v |\mu|} \]
    // Where \( \Omega \) is the proton gyrofrequency, \( v \) is the particle velocity, and \( \mu \) is the pitch-angle cosine.
    double calculateKParallel(double B, double v, double mu) {
        double omega = calculateOmega(B);
        return omega / (v * std::fabs(mu));
    }

    // Heliocentric Dependence of the Power Spectrum (calculateDmuMu)
    // The Kolmogorov turbulence spectrum is given by \( P(k) \propto k^{-5/3} \).
    // The minimum wavenumber \( k_{\text{min}} \) scales as \( \frac{1}{r} \), and the maximum wavenumber \( k_{\text{max}} \) scales as \( \frac{1}{r^2} \).
    double calculateDmuMu(double B, double dB, double v, double mu, double r) {
        double k_parallel = calculateKParallel(B, v, mu);

        // Heliocentric dependence for k_min and k_max
        double k_min_r = k_min_1AU * (r0 / r);
        double k_max_r = k_max_1AU * std::pow(r0 / r, 2);

        // Normalization constant C based on magnetic field fluctuations dB
        double C = (dB * dB) / (3.0 / 2.0) * (std::pow(k_min_r, -2.0 / 3.0) - std::pow(k_max_r, -2.0 / 3.0));

        // Power spectrum follows Kolmogorov's scaling: \( P(k_{\parallel}) = C k_{\parallel}^{-5/3} \)
        double P_k = C * std::pow(k_parallel, -5.0 / 3.0);

        // Pitch-angle diffusion coefficient \( D_{\mu\mu} = \frac{k_{\parallel} P(k_{\parallel})}{B^2} \)
        return (k_parallel * P_k) / (B * B);
    }

    // Spatial Diffusion Coefficient Calculation (calculateDxx)
    // The spatial diffusion coefficient \( D_{xx} \) is related to \( D_{\mu\mu} \) by:
    // \[ D_{xx} = \frac{v^2}{5 D_{\mu\mu}} \]
    // Where the integral over pitch angles is simplified by using a small value for \( \mu \).
    double calculateDxx(double B, double dB, double v, double r) {
        double Dmu_mu = calculateDmuMu(B, dB, v, 0.01, r);
        return v * v / (5.0 * Dmu_mu);
    }

    // Finite Difference Derivative Calculation (calculateDxxDerivative)
    // The derivative \( \frac{d D_{xx}}{dx} \) is calculated using the finite difference method:
    // \[ \frac{d D_{xx}}{dr} \approx \frac{D_{xx}(r + \Delta r) - D_{xx}(r)}{\Delta r} \]
    // Where \( \Delta r \) is a small increment in heliocentric distance.
    double calculateDxxDerivative(double B, double dB, double v, double r, double delta_r) {
        double Dxx_r = calculateDxx(B, dB, v, r);

        double B_next = B * std::pow((r / (r + delta_r)), 2);
        double dB_next = dB * std::pow((r / (r + delta_r)), 2);

        double Dxx_r_next = calculateDxx(B_next, dB_next, v, r + delta_r);

        return (Dxx_r_next - Dxx_r) / delta_r;
    }

    // Wrapper Function (calculateAtHeliocentricDistance)
    // This function calculates the magnetic field \( B(r) \) and fluctuation \( dB(r) \) at a given heliocentric distance \( r \),
    // and then calls the functions to compute \( D_{xx} \) and its derivative.
    // The scaling of the magnetic field is given by \( B(r) = B_{1AU} \left( \frac{r_0}{r} \right)^2 \).
    void calculateAtHeliocentricDistance(double r, double v) {
        double B_1AU = 5e-9;      // Magnetic field at 1 AU (Tesla)
        double dB_1AU = 2e-9;     // Magnetic field fluctuation at 1 AU (Tesla)

        // Magnetic field scaling \( B(r) = B_{1AU} \left( \frac{r_0}{r} \right)^2 \)
        double B = B_1AU * std::pow((r0 / r), 2);

        // Magnetic field fluctuation scaling \( dB(r) = dB_{1AU} \left( \frac{r_0}{r} \right)^2 \)
        double dB = dB_1AU * std::pow((r0 / r), 2);

        std::cout << "At heliocentric distance r = " << r << " meters:\n";
        std::cout << "Magnetic field B = " << B << " Tesla\n";
        std::cout << "Magnetic field fluctuation dB = " << dB << " Tesla\n";

        double Dxx = calculateDxx(B, dB, v, r);
        std::cout << "Dxx = " << Dxx << " m^2/s\n";

        double delta_r = 1e6;  // Small increment for finite difference
        double dDxx_dx = calculateDxxDerivative(B, dB, v, r, delta_r);
        std::cout << "d(Dxx)/dx = " << dDxx_dx << " m^2/s per meter\n";
    }
}



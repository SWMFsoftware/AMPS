#include <cmath>
#include <iostream>
#include <limits>

// Constants
const double q = 1.602e-19; // Proton charge in C
const double m = 1.673e-27; // Proton mass in kg
const double pi = 3.141592653589793;

// Function to calculate the gyrofrequency
double calculateGyrofrequency(double B) {
    return (q * B) / m; // Gyrofrequency Omega = qB/m
}

// Function to calculate the resonant wavenumber k_parallel
double calculateKParallel(double Omega, double v_parallel, double mu) {
    // The historical signature also carries mu, but v_parallel already equals
    // v*mu and is the dimensional resonance denominator. Keep the parameter
    // for source compatibility while making the non-use explicit.
    (void)mu;
    return fabs(Omega / v_parallel );
}

// Function to calculate the power spectrum P(k) assuming a Kolmogorov spectrum
double calculatePowerSpectrum(double k_parallel, double dB, double B, double k_min, double k_max) {
    (void)B;
    if (!std::isfinite(k_parallel) || !std::isfinite(dB) ||
        !std::isfinite(k_min) || !std::isfinite(k_max) || dB < 0.0 ||
        !(k_parallel > 0.0) || !(k_max > k_min) || k_parallel < k_min ||
        k_parallel > k_max) return 0.0;

    // P=C k^-5/3 is normalized by integral P(k)dk=deltaB^2.  P therefore
    // carries T^2 m, which is required for the QLT expression below to return
    // a pitch-angle diffusion rate in s^-1.
    const double bandwidth = pow(k_min, -2.0 / 3.0) -
                             pow(k_max, -2.0 / 3.0);
    if (!(bandwidth > 0.0)) return 0.0;
    double C = (2.0 / 3.0) * dB * dB / bandwidth;
    return C * pow(k_parallel, -5.0 / 3.0);
}

// Function to calculate Dmu_mu
double calculateDmuMu(double dB, double B, double r, double mu, double v_parallel) {
    if (!std::isfinite(dB) || !std::isfinite(B) || !std::isfinite(r) ||
        !std::isfinite(mu) || !std::isfinite(v_parallel) || dB < 0.0 ||
        B <= 0.0 || r <= 0.0 || std::fabs(mu) > 1.0) return 0.0;
    // Step 1: Calculate gyrofrequency Omega
    double Omega = calculateGyrofrequency(B);

    double au=149598000.0E3;
    
    // Step 2: Calculate resonant wavenumber k_parallel
    const double resonantSpeed = std::fabs(v_parallel);
    if (!(resonantSpeed > 0.0)) return 0.0;
    double k_parallel = calculateKParallel(Omega, v_parallel, mu);
    
    // Step 3: Define wavenumber range (k_min and k_max) based on heliocentric distance r
    double k_min_1AU = 1e-6;  // Wavenumber at 1 AU for large structures
    double k_max_1AU = 1e-2;  // Wavenumber at 1 AU for dissipation scale
    double k_min = k_min_1AU * (au / r);      // k_min scales as 1/r
    double k_max = k_max_1AU * (au*au / (r * r)); // k_max scales as 1/r^2
    
    // Step 4: Calculate the power spectrum P(k_parallel)
    double P_k_parallel = calculatePowerSpectrum(k_parallel, dB, B, k_min, k_max);
    
    // Magnetostatic slab QLT in SI units [s^-1].  An out-of-band resonance
    // returns P=0 above rather than being silently clamped to an edge bin.
    double DmuMu = (pi / 2.0) * Omega * Omega * (1.0 - mu * mu) *
                   P_k_parallel / (B * B * resonantSpeed);
    
    return DmuMu;
}

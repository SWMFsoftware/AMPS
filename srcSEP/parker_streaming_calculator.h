/*
================================================================================
                    PARKER STREAMING: HEADER FILE
================================================================================

PURPOSE:
--------
Header file for Parker streaming growth/damping rate calculations using
Kolmogorov turbulence spectrum and isotropic SEP approximation.

NAMESPACE ORGANIZATION:
-----------------------
SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP
  └── Solar Energetic Particle transport with Alfvén wave turbulence
      └── Kolmogorov spectrum assumption  
          └── Isotropic pitch angle distribution approximation

MAIN FUNCTIONS:
---------------
1. CalculateParkerGrowthRatesFromScalar() - Per-energy growth rates
2. CalculateParkerGrowthRatesIntegrated() - Spatially integrated rates
3. CreateMomentumGridFromSEPParams() - Grid from SEP namespace parameters

DEPENDENCIES:
-------------
- AMPS PIC framework (pic.h)
- MPI parallelization (mpi.h)
- SEP namespace parameters
- Relativistic kinematics functions

USAGE:
------
#include "parker_streaming_calculator.h"

using namespace SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP;

// Basic usage with SEP parameters
CalculateParkerGrowthRatesFromScalar(S_scalar, Gamma_plus, Gamma_minus);

// Advanced usage with custom parameters
auto kGrid = CreateLogUniformKGrid(1e-6, 1e-3, 64);
auto pGrid = CreateMomentumGridFromSEPParams();
CalculateParkerGrowthRatesFromScalar(S_scalar, Gamma_plus, Gamma_minus,
                                     kGrid, pGrid, B0, rho);

================================================================================
*/

#ifndef PARKER_STREAMING_CALCULATOR_H
#define PARKER_STREAMING_CALCULATOR_H

// ============================================================================
// SYSTEM INCLUDES
// ============================================================================
#include <vector>             // STL containers
#include <iostream>           // Error output

// ============================================================================
// AMPS FRAMEWORK INCLUDES
// ============================================================================
#include "pic.h"              // AMPS PIC framework
#include <mpi.h>              // MPI parallelization

// ============================================================================
// NAMESPACE DECLARATION
// ============================================================================
namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {
namespace IsotropicSEP {

// ============================================================================
// MAIN CALCULATION FUNCTIONS
// ============================================================================

/**
 * @brief Calculate Parker growth rates from scalar SEP distribution
 * 
 * Computes growth and damping rates (Γ+ and Γ-) directly from scalar particle 
 * distribution S using Parker transport theory and Kolmogorov turbulence spectrum.
 * No intermediate S±(k) arrays are stored - everything computed in single pass.
 * 
 * @param S_scalar     [in]  Scalar SEP distribution S(p) [particles m⁻² s⁻¹ (Δp)⁻¹]
 * @param Gamma_plus   [out] Growth rates Γ+ for outward waves [s⁻¹]
 * @param Gamma_minus  [out] Damping rates Γ- for inward waves [s⁻¹]
 * @param kGrid        [in]  Wavenumber grid (log-uniform) [rad/m]
 * @param pGrid        [in]  Momentum grid [kg⋅m/s]
 * @param B0           [in]  Background magnetic field [T]
 * @param rho          [in]  Background plasma mass density [kg/m³]
 * 
 * @note Uses relativistic kinematics via Relativistic::Momentum2Speed()
 * @note Processes only segments assigned to current MPI rank
 * @note Growth rates calculated per energy bin for each segment
 */
void CalculateParkerGrowthRatesFromScalar(
    const PIC::Datum::cDatumStored& S_scalar,
    const PIC::Datum::cDatumStored& Gamma_plus,
    const PIC::Datum::cDatumStored& Gamma_minus,
    const std::vector<double>& kGrid,
    const std::vector<double>& pGrid,
    double B0,
    double rho
);

/**
 * @brief Calculate spatially integrated Parker growth rates
 * 
 * Computes single growth/damping rate values by integrating over all segments
 * assigned to current MPI process. Useful for global turbulence analysis.
 * 
 * @param S_scalar     [in]  Scalar SEP distribution S(p)
 * @param kGrid        [in]  Wavenumber grid [rad/m]
 * @param pGrid        [in]  Momentum grid [kg⋅m/s]
 * @param B0           [in]  Background magnetic field [T]
 * @param rho          [in]  Background plasma density [kg/m³]
 * @param gammaPlus    [out] Integrated growth rate Γ+ [s⁻¹]
 * @param gammaMinus   [out] Integrated damping rate Γ- [s⁻¹]
 */
void CalculateParkerGrowthRatesIntegrated(
    const PIC::Datum::cDatumStored& S_scalar,
    const std::vector<double>& kGrid,
    const std::vector<double>& pGrid,
    double B0,
    double rho,
    double& gammaPlus,
    double& gammaMinus
);

/**
 * @brief Convenience wrapper using default parameters and SEP momentum grid
 * 
 * Uses default k-grid and gets momentum grid from SEP namespace parameters:
 * - SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_min
 * - SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_max  
 * - SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::n_stream_intervals
 * 
 * @param S_scalar     [in]  Scalar SEP distribution
 * @param Gamma_plus   [out] Growth rates Γ+
 * @param Gamma_minus  [out] Damping rates Γ-
 */
void CalculateParkerGrowthRatesFromScalar(
    const PIC::Datum::cDatumStored& S_scalar,
    const PIC::Datum::cDatumStored& Gamma_plus,
    const PIC::Datum::cDatumStored& Gamma_minus
);

// ============================================================================
// GRID CREATION UTILITY FUNCTIONS
// ============================================================================

/**
 * @brief Create log-uniform wavenumber grid
 * 
 * @param kMin   Minimum wavenumber [rad/m]
 * @param kMax   Maximum wavenumber [rad/m]  
 * @param Nk     Number of grid points
 * @return       Vector of k values in log-uniform spacing
 */
std::vector<double> CreateLogUniformKGrid(double kMin, double kMax, size_t Nk);

/**
 * @brief Create linear momentum grid
 * 
 * @param pMin   Minimum momentum [kg⋅m/s]
 * @param pMax   Maximum momentum [kg⋅m/s]
 * @param Np     Number of grid points
 * @return       Vector of p values in linear spacing
 */
std::vector<double> CreateLinearPGrid(double pMin, double pMax, size_t Np);

/**
 * @brief Create log-uniform momentum grid from log values
 * 
 * @param log_pMin   log(minimum momentum)
 * @param log_pMax   log(maximum momentum)
 * @param Np         Number of grid points
 * @return           Vector of p values in log-uniform spacing
 */
std::vector<double> CreateLogPGrid(double log_pMin, double log_pMax, size_t Np);

/**
 * @brief Create momentum grid from SEP namespace parameters
 * 
 * Reads momentum grid specification from:
 * - SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_min
 * - SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_max
 * - SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::n_stream_intervals
 * 
 * @return  Log-uniform momentum grid [kg⋅m/s]
 * @note    Returns empty vector on error
 */
std::vector<double> CreateMomentumGridFromSEPParams();

// ============================================================================
// PLASMA PHYSICS UTILITY FUNCTIONS
// ============================================================================

/**
 * @brief Calculate plasma parameters from magnetic field and density
 * 
 * @param B0                     [in]  Magnetic field [T]
 * @param rho                    [in]  Mass density [kg/m³]
 * @param VA                     [out] Alfvén speed [m/s]
 * @param Omega                  [out] Cyclotron frequency [rad/s]
 * @param kolmogorov_prefactor   [out] Kolmogorov turbulence prefactor [s⁻¹]
 */
void CalculatePlasmaParameters(
    double B0,
    double rho,
    double& VA,
    double& Omega,
    double& kolmogorov_prefactor
);

// ============================================================================
// PHYSICAL CONSTANTS (for reference)
// ============================================================================

/**
 * @brief Physical constants used in calculations
 * 
 * These are defined in the implementation file:
 * - mu0     = 4π×10⁻⁷ H/m       (permeability of free space)
 * - eCharge = 1.602176634×10⁻¹⁹ C (elementary charge)
 * - mp      = 1.6726219×10⁻²⁷ kg  (proton mass)
 * - CLIGHT  = 2.99792458×10⁸ m/s  (speed of light)
 */

// ============================================================================
// PHYSICS DOCUMENTATION
// ============================================================================

/**
 * @section physics Physics Background
 * 
 * @subsection resonance Wave-Particle Resonance Conditions
 * For charged particles streaming along magnetic field lines, wave-particle
 * resonances occur when:
 * 
 *   k± = Ω / (v⋅μ ± VA)
 * 
 * Where:
 * - k±  = wavenumber for outward (+) or inward (-) propagating waves [rad/m]
 * - Ω   = cyclotron frequency = eB₀/mₚ [rad/s]
 * - v   = particle velocity from Relativistic::Momentum2Speed(p,mₚ) [m/s]
 * - μ   = pitch angle cosine = VA/v (isotropic approximation)
 * - VA  = Alfvén speed = B₀/√(μ₀ρ) [m/s]
 * 
 * @subsection kolmogorov Kolmogorov Turbulence
 * For Kolmogorov turbulence with spectrum P(k) ∝ k^(-5/3), the growth/damping
 * rates are given by:
 * 
 *   Γ± = C × Σ S±(k)/k
 * 
 * Where:
 * - C = (π²e²VA)/(cB₀²) = Kolmogorov prefactor from quasi-linear theory
 * - S±(k) = wave-sense particle flux from resonance conditions
 * 
 * @subsection algorithm Algorithm Summary
 * 1. For each momentum bin in scalar distribution S(p):
 *    - Calculate relativistic velocity v = Relativistic::Momentum2Speed(p,mₚ)
 *    - Determine resonant wavenumbers k± using Parker theory
 *    - Split scalar flux between ± modes based on validity
 *    - Directly accumulate contribution to Σ S±/k using CIC interpolation
 * 2. Apply Kolmogorov prefactor to get final growth/damping rates
 * 3. Store results in Γ+ and Γ- data arrays for each segment
 */

} // namespace IsotropicSEP
} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP

#endif // PARKER_STREAMING_CALCULATOR_H

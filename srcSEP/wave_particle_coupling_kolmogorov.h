/*
================================================================================
                    WAVE-PARTICLE ENERGY COUPLING HEADER
================================================================================

PURPOSE:
--------
Header file for wave-particle energy coupling functions in Alfvén wave turbulence 
for solar energetic particle (SEP) transport simulations. Provides interface for 
bidirectional energy exchange between particles and turbulent magnetic field 
fluctuations using Parker transport theory.

MAIN FUNCTIONS:
---------------

1. UpdateWaveEnergyWithParticleCoupling()
   - Updates integrated wave energies for a single field line segment
   - Calculates Parker growth rates from particle distribution
   - Evolves wave energies: E±^{n+1} = E±^n + Δt[2γ± E± - E±/τ_cas + Q_shock × V_cell]
   - Redistributes energy gain/loss among particles attached to segment

2. UpdateAllSegmentsWaveEnergyWithParticleCoupling()
   - Global function that loops through all field lines and segments
   - Calls single-segment function for each segment assigned to MPI process
   - Accesses wave data from SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity

3. Diagnostic Functions:
   - CalculateTotalWaveEnergyInSystem() - Total wave energy across simulation
   - CalculateTotalParticleEnergyInSystem() - Total particle energy across simulation  
   - CheckEnergyConservation() - Monitors energy conservation violations

USAGE EXAMPLES:
---------------

// Update single segment with full parameter control
PIC::FieldLine::cFieldLineSegment* segment = "get segment";
double E_plus_old = 1.5e12;  // 1.5 TJ outward wave energy [J]
double E_minus_old = 0.8e12; // 0.8 TJ inward wave energy [J]
double E_plus_new, E_minus_new;
double dt = 10.0; // 10 second time step
double tau_cas = 1000.0; // cascade time
double Q_shock = 0.01; // shock injection

UpdateWaveEnergyWithParticleCoupling(
    segment, E_plus_old, E_minus_old, particle_distribution, dt,
    tau_cas, Q_shock, E_plus_new, E_minus_new, B0, density
);

// Simplified global update with defaults
UpdateAllSegmentsWaveEnergyWithParticleCoupling(particle_distribution, dt);

// Energy conservation monitoring
CheckEnergyConservation(B0, true); // verbose output

================================================================================
*/

#ifndef WAVE_PARTICLE_COUPLING_H
#define WAVE_PARTICLE_COUPLING_H

#include "pic.h"              // AMPS PIC framework
#include <mpi.h>              // MPI parallelization
#include <vector>             // STL containers
#include <iostream>           // Error output
#include <cmath>              // Mathematical functions

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {
namespace IsotropicSEP {

// ============================================================================
// SINGLE SEGMENT WAVE-PARTICLE COUPLING FUNCTIONS
// ============================================================================

/**
 * @brief Update wave energy with particle coupling for a single segment
 * 
 * Main function that calculates Parker growth rates and evolves integrated wave energies
 * while redistributing energy among particles. Uses iterative energy redistribution
 * with energy floor protection (10% minimum particle energy).
 * 
 * @param segment           Pointer to field line segment to process
 * @param E_plus_initial    Initial outward integrated wave energy [J]
 * @param E_minus_initial   Initial inward integrated wave energy [J]
 * @param S_scalar          Scalar particle distribution
 * @param dt                Time step [s]
 * @param Q_shock           Shock injection rate [dimensionless]
 * @param E_plus_final      Output: final outward integrated wave energy [J]
 * @param E_minus_final     Output: final inward integrated wave energy [J]
 * @param B0                Background magnetic field [T]
 * @param rho               Mass density [kg/m³]
 */
void UpdateWaveEnergyWithParticleCoupling(
    int field_line_idx,
    PIC::FieldLine::cFieldLineSegment* segment,
    double& E_plus_initial,
    double& E_minus_initial,
    const PIC::Datum::cDatumStored& S_scalar,
    double dt,
    double& E_plus_final,
    double& E_minus_final,
    double B0,
    double rho
);

/**
 * @brief Convenience wrapper with default parameters
 * 
 * @param segment           Pointer to field line segment to process
 * @param E_plus_initial    Initial outward integrated wave energy [J]
 * @param E_minus_initial   Initial inward integrated wave energy [J]
 * @param S_scalar          Scalar particle distribution
 * @param dt                Time step [s]
 * @param E_plus_final      Output: final outward integrated wave energy [J]
 * @param E_minus_final     Output: final inward integrated wave energy [J]
 */
void UpdateWaveEnergyWithParticleCoupling(
    int field_line_idx,
    PIC::FieldLine::cFieldLineSegment* segment,
    double& E_plus_initial,
    double& E_minus_initial,
    const PIC::Datum::cDatumStored& S_scalar,
    double dt,
    double& E_plus_final,
    double& E_minus_final
);

// ============================================================================
// GLOBAL WAVE-PARTICLE COUPLING FUNCTIONS
// ============================================================================

/**
 * @brief Update wave energy for all field line segments (MPI parallel)
 * 
 * Global function that loops through all field lines and segments, calling
 * the single-segment update function for each segment assigned to the current
 * MPI process. Accesses wave data from WaveEnergyDensity datum.
 * 
 * @param S_scalar          Scalar particle distribution
 * @param dt                Time step [s]
 * @param Q_shock           Shock injection rate [dimensionless]
 * @param B0                Background magnetic field [T]
 * @param rho               Mass density [kg/m³]
 */
void UpdateAllSegmentsWaveEnergyWithParticleCoupling(
    PIC::Datum::cDatumStored& WaveEnergyDensity,
    PIC::Datum::cDatumStored& S_scalar,
    double dt 
);

// ============================================================================
// DIAGNOSTIC AND UTILITY FUNCTIONS
// ============================================================================

/**
 * @brief Calculate total particle energy in a single segment
 * 
 * @param segment           Pointer to field line segment
 * @return Total particle energy in segment [J]
 */
double CalculateTotalParticleEnergyInSegment(PIC::FieldLine::cFieldLineSegment* segment);

/**
 * @brief Calculate wave energy density from integrated wave energies
 * 
 * @param E_plus            Integrated outward wave energy [J]
 * @param E_minus           Integrated inward wave energy [J]
 * @param V_segment         Segment volume [m³]
 * @param B0                Background magnetic field [T]
 * @return Wave energy density [J/m³]
 */
double CalculateWaveEnergyDensity(double E_plus, double E_minus, double V_segment, double B0);

/**
 * @brief Convert total energy density to integrated wave energies
 * 
 * @param energy_density_total  Total wave energy density [J/m³]
 * @param energy_ratio          Ratio E_plus / E_minus
 * @param V_segment            Segment volume [m³]
 * @param E_plus               Output: integrated outward wave energy [J]
 * @param E_minus              Output: integrated inward wave energy [J]
 */
void ConvertEnergyDensityToIntegratedEnergies(
    double energy_density_total,
    double energy_ratio,
    double V_segment,
    double& E_plus,
    double& E_minus
);

/**
 * @brief Calculate total wave energy across entire simulation domain
 * 
 * @param B0                Background magnetic field [T]
 * @return Total wave energy in system [J]
 */
double CalculateTotalWaveEnergyInSystem(double B0);

/**
 * @brief Calculate total particle energy across entire simulation domain
 * 
 * @return Total particle energy in system [J]
 */
double CalculateTotalParticleEnergyInSystem();

/**
 * @brief Monitor energy conservation in wave-particle system
 * 
 * @param B0                Background magnetic field [T]
 * @param verbose          Enable detailed output
 */
void CheckEnergyConservation(double B0, bool verbose = false);

// ============================================================================
// PHYSICS CONSTANTS (for reference)
// ============================================================================

namespace PhysicsConstants {
    constexpr double MU0         = 4.0e-7 * M_PI;      // Permeability [H/m]
    constexpr double E_CHARGE    = 1.602176634e-19;    // Elementary charge [C]
    constexpr double C_LIGHT     = 2.99792458e8;       // Speed of light [m/s]
}

namespace TypicalSolarWind {
    constexpr double B0_TYPICAL     = 5.0e-9;          // 5 nT magnetic field
    constexpr double RHO_TYPICAL    = 5.0e-21;         // kg/m³ mass density
    constexpr double Q_SHOCK_DEFAULT = 0.0;            // No shock injection
}

} // namespace IsotropicSEP
} // namespace AlfvenTurbulence_Kolmogorov  
} // namespace SEP

#endif // WAVE_PARTICLE_COUPLING_H

/*
================================================================================
                    TURBULENCE ENERGY ADVECTION HEADER FILE
================================================================================

PURPOSE:
--------
Header file for turbulence energy advection functions in Alfvén wave turbulence 
for solar energetic particle (SEP) transport simulations. Provides interface for 
Alfvén wave energy transport along magnetic field lines using Lagrangian mesh 
with spatially varying wave speeds calculated from local plasma conditions.

MAIN FUNCTIONS:
---------------

1. AdvectTurbulenceEnergyAlongFieldLine()
   - Advects integrated wave energies along a single field line
   - Calculates local Alfvén velocities from vertex magnetic field and plasma density
   - Uses conservative upwind scheme for energy transport
   - Applies energy updates only to segments assigned to current MPI process

2. AdvectTurbulenceEnergyAllFieldLines()
   - Global function that processes all field lines in parallel
   - Calls single field line function for each field line with local segments
   - MPI-aware: processes only field lines with segments assigned to current process

3. CalculateTotalAdvectedEnergy()
   - Diagnostic function to monitor total wave energy in system
   - Sums energy across all local segments with MPI reduction
   - Used for energy conservation verification

PHYSICS:
--------
- Lagrangian mesh: Field line segments move with solar wind plasma
- Wave propagation: E± transport at local Alfvén speeds V_A = |B|/√(μ₀ρ)
- Energy flux: Flux = V_A × (E_total/V_segment) × A_boundary
- Conservative scheme: Total energy conserved during advection
- Boundary conditions: Inward waves exit at Sun, outward waves damped at outer boundary

MPI PARALLELIZATION:
--------------------
- Global data access: Wave energies read from all segments
- Local processing: Flux calculations only for local segments (segment->Thread == PIC::ThisThread)
- Local updates: Energy changes applied only to local segments
- Conservative fluxes: Energy transfer between segments properly handled across processes

USAGE EXAMPLES:
---------------

// Single field line advection
int field_line_idx = 5;
double dt = 10.0;  // seconds
SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAlongFieldLine(
    field_line_idx,
    SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,
    dt,
    true  // apply boundary conditions
);

// All field lines advection
SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(
    SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,
    dt
);

// Energy conservation check
double total_energy = SEP::AlfvenTurbulence_Kolmogorov::CalculateTotalAdvectedEnergy(
    WaveEnergyDensity
);

REQUIREMENTS:
-------------
- AMPS PIC framework with field line structure
- SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity datum
- SEP::FieldLine::GetSegmentVolume() and SEP::FieldLine::MagneticTubeRadius()
- FL::DatumAtVertexMagneticField and FL::DatumAtVertexPlasmaDensity
- MPI parallelization support

================================================================================
*/

#ifndef TURBULENCE_ADVECTION_H
#define TURBULENCE_ADVECTION_H

#include "pic.h"              // AMPS PIC framework
#include <mpi.h>              // MPI parallelization

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

// ============================================================================
// TURBULENCE ENERGY ADVECTION FUNCTIONS
// ============================================================================

/**
 * @brief Advect turbulence energy along a single field line with local Alfvén speeds
 * 
 * Main function that handles energy transport of integrated wave energies E± along
 * a magnetic field line. Calculates local Alfvén velocities at each segment boundary
 * from vertex magnetic field and plasma density data. Uses conservative upwind scheme
 * to maintain energy conservation. Only processes segments assigned to current MPI process.
 * 
 * Algorithm:
 * 1. Read wave energy data from ALL segments (global information)
 * 2. Calculate local Alfvén velocities V_A = |B|/√(μ₀ρ) at each boundary
 * 3. Calculate energy fluxes using upwind scheme (only for local segments)
 * 4. Apply conservative energy updates (only to local segments)
 * 5. Apply boundary conditions (only to local boundary segments)
 * 6. Write updated wave energies (only to local segments)
 * 
 * @param field_line_idx        Index of field line to process
 * @param WaveEnergyDensity     Wave energy datum containing E+ and E- data
 * @param dt                    Time step [s]
 * @param apply_boundary_conditions  Enable boundary conditions at inner/outer boundaries
 * 
 * @note Wave energy data is read globally but updates applied only to local segments
 * @note Requires valid magnetic field and plasma density data at all vertices
 * @note Energy flux calculation: Flux = V_A × (E_total/V_segment) × A_boundary
 */
void AdvectTurbulenceEnergyAlongFieldLine(
    int field_line_idx,
    const PIC::Datum::cDatumStored& WaveEnergyDensity,
    double dt,
    bool apply_boundary_conditions = true
);

/**
 * @brief Advect turbulence energy for all field lines in MPI parallel fashion
 * 
 * Global function that processes all field lines in the simulation. For each field line,
 * checks if it has segments assigned to the current MPI process and calls the single
 * field line advection function. Provides automatic load balancing across MPI processes.
 * 
 * @param WaveEnergyDensity     Wave energy datum containing E+ and E- data
 * @param dt                    Time step [s]
 * @param apply_boundary_conditions  Enable boundary conditions at inner/outer boundaries
 * 
 * @note Only processes field lines that have at least one segment assigned to current process
 * @note Automatically handles MPI load balancing
 * @note May require MPI_Barrier() depending on implementation needs
 */
void AdvectTurbulenceEnergyAllFieldLines(
    const PIC::Datum::cDatumStored& WaveEnergyDensity,
    double dt,
    bool apply_boundary_conditions = true
);

/**
 * @brief Calculate total wave energy across entire simulation domain
 * 
 * Diagnostic function that sums wave energies (E+ + E-) across all segments
 * assigned to the current MPI process, then performs MPI reduction to get
 * global total. Used for monitoring energy conservation during advection.
 * 
 * @param WaveEnergyDensity     Wave energy datum containing E+ and E- data
 * @return Total wave energy in simulation [J]
 * 
 * @note Performs MPI_Allreduce to sum energy across all processes
 * @note Only counts energy from segments assigned to current process
 * @note Returns same value on all MPI processes after reduction
 */
double CalculateTotalAdvectedEnergy(
    const PIC::Datum::cDatumStored& WaveEnergyDensity
);

// ============================================================================
// PHYSICS CONSTANTS AND TYPICAL VALUES
// ============================================================================

namespace PhysicsConstants {
    constexpr double MU0            = 4.0e-7 * M_PI;      // Permeability [H/m]
    constexpr double PROTON_MASS    = 1.67262192e-27;     // Proton mass [kg]
    constexpr double C_LIGHT        = 2.99792458e8;       // Speed of light [m/s]
}

namespace TypicalSolarWind {
    constexpr double B_FIELD_1AU    = 5.0e-9;             // 5 nT magnetic field at 1 AU
    constexpr double DENSITY_1AU    = 5.0e6;              // 5 cm⁻³ number density at 1 AU
    constexpr double V_SW_TYPICAL   = 400.0e3;            // 400 km/s solar wind speed
    constexpr double V_A_TYPICAL    = 50.0e3;             // 50 km/s Alfvén speed (typical)
    constexpr double ENERGY_SCALE   = 1.0e12;             // TJ scale for wave energies
}

// ============================================================================
// ADVANCED USAGE EXAMPLES AND INTEGRATION PATTERNS
// ============================================================================

/*
EXAMPLE 1: Basic integration in simulation timestep
-------------------------------------------------------
void SimulationTimestep(double dt) {
    // ... particle transport, magnetic field evolution, etc. ...
    
    // Advect turbulence energy
    SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(
        SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,
        dt,
        true  // apply boundary conditions
    );
    
    // ... continue with other physics modules ...
}

EXAMPLE 2: Energy conservation monitoring
------------------------------------------
double energy_before = SEP::AlfvenTurbulence_Kolmogorov::CalculateTotalAdvectedEnergy(
    WaveEnergyDensity
);

SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(
    WaveEnergyDensity, dt
);

double energy_after = SEP::AlfvenTurbulence_Kolmogorov::CalculateTotalAdvectedEnergy(
    WaveEnergyDensity
);

double energy_change = energy_after - energy_before;
double relative_change = energy_change / energy_before;

if (std::abs(relative_change) > 1.0e-6) {
    std::cout << "Warning: Energy conservation violation: " 
              << relative_change * 100.0 << "%" << std::endl;
}

EXAMPLE 3: Adaptive time stepping with CFL condition
-----------------------------------------------------
double CalculateAdaptiveTimestep(double dx_min, double V_A_max, double CFL_factor = 0.5) {
    return CFL_factor * dx_min / V_A_max;
}

void AdaptiveAdvection(const PIC::Datum::cDatumStored& WaveEnergyDensity, 
                      double dt_desired) {
    double dx_min =  calculate minimum segment length ;
    double V_A_max =  calculate maximum Alfvén speed ;
    double dt_stable = CalculateAdaptiveTimestep(dx_min, V_A_max);
    
    double dt_actual = std::min(dt_desired, dt_stable);
    int num_substeps = static_cast<int>(std::ceil(dt_desired / dt_actual));
    dt_actual = dt_desired / num_substeps;
    
    for (int step = 0; step < num_substeps; ++step) {
        SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(
            WaveEnergyDensity, dt_actual
        );
    }
}

EXAMPLE 4: Field line specific processing
------------------------------------------
void ProcessSpecificFieldLines(const std::vector<int>& field_line_indices,
                               const PIC::Datum::cDatumStored& WaveEnergyDensity,
                               double dt) {
    for (int field_line_idx : field_line_indices) {
        // Check if this process has segments on this field line
        bool has_local_segments = false;
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        
        for (int i = 0; i < num_segments; ++i) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
            if (segment && segment->Thread == PIC::ThisThread) {
                has_local_segments = true;
                break;
            }
        }
        
        if (has_local_segments) {
            SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAlongFieldLine(
                field_line_idx, WaveEnergyDensity, dt
            );
        }
    }
}
*/

} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP

#endif // TURBULENCE_ADVECTION_H

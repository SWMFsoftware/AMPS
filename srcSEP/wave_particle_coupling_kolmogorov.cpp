/* 
================================================================================
                    WAVE-PARTICLE ENERGY COUPLING FUNCTIONS
================================================================================

PURPOSE:
--------
This file implements wave-particle energy coupling for Alfvén wave turbulence 
in solar energetic particle (SEP) transport simulations. The functions handle 
the bidirectional energy exchange between particles and turbulent magnetic 
field fluctuations using Parker transport theory.

MAIN FUNCTIONS:
---------------

1. UpdateWaveEnergyWithParticleCoupling()
   - Updates integrated wave energies for a single field line segment
   - Calculates Parker growth rates from particle distribution
   - Evolves integrated wave energies: E±^{n+1} = E±^n + Δt[2γ± E± - E±/τ_cas + Q_shock × V_cell]
   - Redistributes energy gain/loss among particles attached to segment

2. UpdateAllSegmentsWaveEnergyWithParticleCoupling()
   - Global function that loops through all field lines and segments
   - Calls single-segment function for each segment assigned to MPI process
   - Accesses wave data from SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity

3. Diagnostic Functions:
   - CalculateTotalWaveEnergyInSystem() - Total wave energy across simulation
   - CalculateTotalParticleEnergyInSystem() - Total particle energy across simulation  
   - CheckEnergyConservation() - Monitors energy conservation violations

PHYSICS:
--------
- Uses Parker transport theory for wave-particle resonances: k± = Ω/(v‖ ± VA)
- Implements Kolmogorov turbulence spectrum scaling: Γ± = C × Σ S±(k)/k
- Ensures energy conservation: wave energy gain = -particle energy loss
- Handles relativistic particle dynamics with proper statistical weights
- Applies energy floors (10% minimum) to prevent particle energy depletion

USAGE EXAMPLES:
---------------

// Example 1: Update single segment with full parameter control
PIC::FieldLine::cFieldLineSegment* segment = /get segment/;
double A_plus_old = 0.1, A_minus_old = 0.05;
double A_plus_new, A_minus_new;
double dt = 10.0; // 10 second time step
double tau_cas = 1000.0; // cascade time
double Q_shock = 0.01; // shock injection

UpdateWaveEnergyWithParticleCoupling(
    segment, A_plus_old, A_minus_old, particle_distribution, dt,
    tau_cas, Q_shock, A_plus_new, A_minus_new, B0, density
);

// Example 2: Update single segment with default parameters
UpdateWaveEnergyWithParticleCoupling(
    segment, A_plus_old, A_minus_old, particle_distribution, dt,
    A_plus_new, A_minus_new
);

// Example 3: Update all segments in simulation (typical main loop usage)
UpdateAllSegmentsWaveEnergyWithParticleCoupling(
    particle_distribution, dt, tau_cas, Q_shock, B0, density
);

// Example 4: Simplified global update with defaults
UpdateAllSegmentsWaveEnergyWithParticleCoupling(particle_distribution, dt);

// Example 5: Energy conservation monitoring
double B0 = 5.0e-9; // Tesla
CheckEnergyConservation(B0, true); // verbose output

// Example 6: Complete simulation timestep
void SimulationTimestep(double dt) {
    // ... other physics updates ...
    
    // Update wave-particle coupling
    UpdateAllSegmentsWaveEnergyWithParticleCoupling(
        SEP_scalar_distribution, dt
    );
    
    // Check energy conservation every 10 steps
    static int step_counter = 0;
    if (++step_counter % 10 == 0) {
        CheckEnergyConservation(magnetic_field_strength);
    }
    
    // ... continue with other updates ...
}

REQUIREMENTS:
-------------
- AMPS PIC framework with field line structure
- SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity datum for wave amplitudes
- MPI parallelization support
- Relativistic particle dynamics (Relativistic::Vel2E function)
- Parker growth rate calculation functions (from parker_streaming_calculator.cpp)

PERFORMANCE NOTES:
------------------
- Memory efficient: No intermediate arrays stored during calculation
- MPI parallel: Each process handles assigned field line segments only
- Iterative energy redistribution: Handles energy floor constraints robustly
- Direct particle loop: Works with actual computational particles, not distributions

================================================================================
*/

#include "sep.h"

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {
namespace IsotropicSEP {

// ============================================================================
// WAVE ENERGY EVOLUTION WITH PARTICLE COUPLING
// ============================================================================

void UpdateWaveEnergyWithParticleCoupling(
    int field_line_idx,
    PIC::FieldLine::cFieldLineSegment* segment,
    double& E_plus_initial,                    // Initial integrated outward wave energy [J]
    double& E_minus_initial,                   // Initial integrated inward wave energy [J]
    const PIC::Datum::cDatumStored& S_scalar,  // Scalar particle distribution
    double dt,                                 // Time step [s]
    double& E_plus_final,                     // Output: final outward wave energy [J]
    double& E_minus_final,                    // Output: final inward wave energy [J]
    double B0,                                // Background magnetic field [T]
    double rho                                // Mass density [kg/m³]
) {
    // ========================================================================
    // INPUT VALIDATION
    // ========================================================================
    if (!segment) {
        std::cerr << "Error: Null segment pointer in UpdateWaveEnergyWithParticleCoupling" << std::endl;
        E_plus_final = E_plus_initial;
        E_minus_final = E_minus_initial;
        return;
    }

    if (segment->Thread != PIC::ThisThread) {
        E_plus_final = E_plus_initial;
        E_minus_final = E_minus_initial;
        return;
    }
    
    if (dt <= 0.0) {
        std::cerr << "Error: Invalid time parameters (dt=" << dt << ")" << std::endl;
        E_plus_final = E_plus_initial;
        E_minus_final = E_minus_initial;
        return;
    }
    
    // ========================================================================
    // CALCULATE SEGMENT VOLUME
    // ========================================================================
    double V_cell = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx); 
    
    if (V_cell <= 0.0) {
        std::cerr << "Warning: Invalid cell volume (" << V_cell << ") in segment" << std::endl;
        E_plus_final = E_plus_initial;
        E_minus_final = E_minus_initial;
        return;
    }
    
    // ========================================================================
    // STEP 1: CALCULATE GROWTH RATES FOR THIS SEGMENT
    // ========================================================================
    double gamma_plus = 0.0;   // Growth rate for outward waves [s⁻¹]
    double gamma_minus = 0.0;  // Growth rate for inward waves [s⁻¹]
    
    // Use single segment Parker growth rate calculation
    CalculateParkerGrowthRatesFromScalarSingleSegment(
        segment, S_scalar, gamma_plus, gamma_minus
    );
    
    // ========================================================================
    // STEP 2: INTEGRATED WAVE ENERGY EVOLUTION
    // ========================================================================
    // From A±^{n+1} = A±^n + Δt[2γ± A±] 
    // Since E± = A± * V_cell, we get:
    // E±^{n+1} = E±^n + Δt[2γ± E± ]
    
    double dE_plus_dt = 2.0 * gamma_plus * E_plus_initial; 
    double dE_minus_dt = 2.0 * gamma_minus * E_minus_initial; 
    
    // Update integrated wave energies directly
    E_plus_final = E_plus_initial + dt * dE_plus_dt;
    E_minus_final = E_minus_initial + dt * dE_minus_dt;
    
    // Ensure non-negative energies
    E_plus_final = std::max(0.0, E_plus_final);
    E_minus_final = std::max(0.0, E_minus_final);
    
    // ========================================================================
    // STEP 3: CALCULATE TOTAL WAVE ENERGY CHANGE
    // ========================================================================
    double E_wave_gain_plus = dE_plus_dt * dt;    // [J]
    double E_wave_gain_minus = dE_minus_dt * dt;  // [J]
    
    double total_wave_energy_change = E_wave_gain_plus + E_wave_gain_minus; // [J]
    
    // ========================================================================
    // STEP 5: DISTRIBUTE ENERGY CHANGE TO PARTICLES
    // ========================================================================
    // Energy conservation: particles lose energy equal to wave energy gain
    double particle_energy_change = -total_wave_energy_change;
    
    // FIRST LOOP: Calculate total particle energy and count particles
    double total_particle_energy = 0.0;
    int total_particle_count = 0;
    long int p = segment->FirstParticleIndex;
    
    while (p != -1) {
        double v[3];
        PIC::ParticleBuffer::GetV(v, p);
        
        // Calculate particle kinetic energy
        double v_magnitude = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        
        // Relativistic kinetic energy calculation
        double particle_mass = _H__MASS_; // [kg]
        double kinetic_energy = Relativistic::Vel2E(v, particle_mass); // [J]
        
        // Get particle statistical weight
        double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] * 
                            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        
        // Model particle energy (physical energy represented by this computational particle)
        double model_particle_energy = kinetic_energy * stat_weight;
        
        total_particle_energy += model_particle_energy;
        total_particle_count++;
        
        // Get next particle
        p = PIC::ParticleBuffer::GetNext(p);
    }
    
    if (total_particle_count == 0 || total_particle_energy <= 0.0) {
        // No particles to redistribute energy to
        return;
    }
    
    // Check if energy removal exceeds available energy
    if (particle_energy_change < 0 && std::abs(particle_energy_change) > total_particle_energy) {
        char error_msg[512];
        sprintf(error_msg, "Energy removal (%.6e J) exceeds total particle energy (%.6e J) in segment", 
                std::abs(particle_energy_change), total_particle_energy);
        exit( __LINE__, __FILE__,error_msg);
    }
    
    // ITERATIVE ENERGY REDISTRIBUTION LOOP
    double remaining_energy_to_distribute = particle_energy_change;
    
    while (std::abs(remaining_energy_to_distribute) > 1.0e-20) {
        // Split energy equally among all particles
        double energy_per_particle = remaining_energy_to_distribute / total_particle_count;
        double energy_actually_redistributed = 0.0;
        
        // SECOND LOOP: Update particle velocities
        p = segment->FirstParticleIndex;
        
        while (p != -1) {
            double v_current[3];
            PIC::ParticleBuffer::GetV(v_current, p);
            
            // Calculate current particle kinetic energy
            double particle_mass = _H__MASS_;
            double kinetic_energy = Relativistic::Vel2E(v_current, particle_mass);
            
            // Get particle statistical weight
            double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] * 
                                PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
            
            // Model particle energy (physical energy represented by this computational particle)
            double current_energy = kinetic_energy * stat_weight;
            
            // Determine energy change for this particle
            double energy_change_this_particle = energy_per_particle;
            
            // If removing energy, check energy floor (10% of current energy)
            if (energy_change_this_particle < 0) {
                double energy_floor = 0.1 * current_energy;
                double max_removable = current_energy - energy_floor;
                
                // Limit energy removal to not go below floor
                if (std::abs(energy_change_this_particle) > max_removable) {
                    energy_change_this_particle = -max_removable;
                }
            }
            
            // Calculate new kinetic energy (model particle level)
            double new_kinetic_energy = current_energy + energy_change_this_particle;
            
            // Ensure positive energy (additional safety check)
            if (new_kinetic_energy <= 0.0) {
                new_kinetic_energy = 0.1 * current_energy; // 10% floor
                energy_change_this_particle = new_kinetic_energy - current_energy;
            }
            
            // Track actual energy redistributed
            energy_actually_redistributed += energy_change_this_particle;
            
            // Convert back to single-particle kinetic energy for velocity calculation
            double new_single_particle_kinetic_energy = new_kinetic_energy / stat_weight;
            
            // Convert back to velocity using AMPS relativistic functions
            // Note: This requires implementing inverse of Vel2E or using iterative approach
            // For now, using the original approach but with proper mass
            double v_magnitude = std::sqrt(v_current[0]*v_current[0] + v_current[1]*v_current[1] + v_current[2]*v_current[2]);
            double gamma_new = new_single_particle_kinetic_energy / (particle_mass * PhysicsConstants::C_LIGHT * PhysicsConstants::C_LIGHT) + 1.0;
            
            // Calculate new velocity magnitude
            double v_new_magnitude = PhysicsConstants::C_LIGHT * std::sqrt(1.0 - 1.0/(gamma_new*gamma_new));
            
            if (v_magnitude > 1.0e-20) {
                // Scale velocity to new magnitude while preserving direction
                double scale_factor = v_new_magnitude / v_magnitude;
                double v_new[3] = {
                    v_current[0] * scale_factor,
                    v_current[1] * scale_factor,
                    v_current[2] * scale_factor
                };
                
                // Update particle velocity in the buffer
                PIC::ParticleBuffer::SetV(v_new, p);
            }
            
            // Get next particle
            p = PIC::ParticleBuffer::GetNext(p);
        }
        
        // Update remaining energy to distribute
        remaining_energy_to_distribute -= energy_actually_redistributed;
        
        // Safety check to prevent infinite loop
        if (std::abs(remaining_energy_to_distribute) < 1.0e-20) {
            break;
        }
        
        // If we can't redistribute any more energy, exit
        if (std::abs(energy_actually_redistributed) < 1.0e-20) {
            if (std::abs(remaining_energy_to_distribute) > 1.0e-15) {
                std::cerr << "Warning: Could not redistribute all energy. Remaining: " 
                          << remaining_energy_to_distribute << " J" << std::endl;
            }
            break;
        }
    }
}

// ============================================================================
// CONVENIENCE WRAPPER WITH DEFAULT PARAMETERS
// ============================================================================

void UpdateWaveEnergyWithParticleCoupling(
    int field_line_idx,
    PIC::FieldLine::cFieldLineSegment* segment,
    double& E_plus_initial,
    double& E_minus_initial,
    const PIC::Datum::cDatumStored& S_scalar,
    double dt,
    double& E_plus_final,
    double& E_minus_final
) {
    // Default parameters for typical solar wind conditions
    const double Q_shock_default = 0.0;       // No shock injection
    const double B0_default = 5.0e-9;         // 5 nT
    const double rho_default = 5.0e-21;       // kg/m³
    
    UpdateWaveEnergyWithParticleCoupling(
        field_line_idx,segment, E_plus_initial, E_minus_initial, S_scalar, dt,
        E_plus_final, E_minus_final,
        B0_default, rho_default
    );
}

// ============================================================================
// UTILITY FUNCTION: CALCULATE TOTAL PARTICLE ENERGY IN SEGMENT
// ============================================================================

double CalculateTotalParticleEnergyInSegment(PIC::FieldLine::cFieldLineSegment* segment) {
    if (!segment) return 0.0;
    
    double total_energy = 0.0;
    long int p = segment->FirstParticleIndex;
    
    while (p != -1) {
        double v[3];
        PIC::ParticleBuffer::GetV(v, p);
        
        // Calculate particle kinetic energy
        double v_magnitude = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        
        // Relativistic kinetic energy
        double particle_mass = _H__MASS_;
        double kinetic_energy = Relativistic::Vel2E(v, particle_mass);
        
        // Get particle statistical weight
        double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] * 
                            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        
        // Model particle energy
        double model_particle_energy = kinetic_energy * stat_weight;
        
        total_energy += model_particle_energy;
        
        // Get next particle
        p = PIC::ParticleBuffer::GetNext(p);
    }
    
    return total_energy; // [J]
}

double CalculateWaveEnergyDensity(double E_plus, double E_minus, double V_segment, double B0) {
    const double mu0 = 4.0e-7 * M_PI;
    
    // Convert integrated energies back to amplitudes
    double A_plus = E_plus / V_segment;
    double A_minus = E_minus / V_segment;
    
    // Wave energy density: ε = (B₀²/2μ₀) * (A₊² + A₋²)
    // where A± are normalized wave amplitudes
    double energy_density = (B0 * B0 / (2.0 * mu0)) * (A_plus * A_plus + A_minus * A_minus);
    
    return energy_density;  // [J/m³]
}

// ============================================================================
// UTILITY FUNCTION: CONVERT BETWEEN AMPLITUDE AND ENERGY DENSITY
// ============================================================================

void ConvertEnergyDensityToIntegratedEnergies(
    double energy_density_total,  // Total wave energy density [J/m³]
    double energy_ratio,          // E_plus / E_minus
    double V_segment,             // Segment volume [m³]
    double& E_plus,              // Output: integrated outward wave energy [J]
    double& E_minus              // Output: integrated inward wave energy [J]
) {
    if (energy_density_total <= 0.0 || V_segment <= 0.0) {
        E_plus = E_minus = 0.0;
        return;
    }
    
    // Total integrated energy in segment
    double total_integrated_energy = energy_density_total * V_segment;
    
    // From: E_plus = r * E_minus, where r = energy_ratio
    // And: total = E_plus + E_minus = r * E_minus + E_minus = E_minus * (r + 1)
    // Solve for E_minus: E_minus = total / (r + 1)
    
    E_minus = total_integrated_energy / (energy_ratio + 1.0);
    E_plus = energy_ratio * E_minus;
}

// ============================================================================
// GLOBAL UPDATE FUNCTION FOR ALL SEGMENTS
// ============================================================================

void UpdateAllSegmentsWaveEnergyWithParticleCoupling(
    PIC::Datum::cDatumStored& WaveEnergyDensity,
    PIC::Datum::cDatumStored& S_scalar,
    double dt
) {
    double B0,rho;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int processed_segments = 0;

    // Loop through all field lines
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();

        // Loop through all segments in this field line
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);

            // Only process segments assigned to this MPI process
            if (!segment || segment->Thread != PIC::ThisThread) {
                continue;
            }

	    //Get background plasma parameters 
            segment->GetPlasmaDensity(0.5,rho);  // Number density 
            rho*=_H__MASS_; //convert to mass density

            // Get current wave energy data
            double* wave_data = segment->GetDatum_ptr(WaveEnergyDensity);
            if (!wave_data) continue;

            double E_plus_initial = wave_data[0];   // Outward wave energy [J]
            double E_minus_initial = wave_data[1];  // Inward wave energy [J]

            double E_plus_final, E_minus_final;

            // Update wave energy with particle coupling
            UpdateWaveEnergyWithParticleCoupling(
                field_line_idx, segment,
                E_plus_initial, E_minus_initial,
                S_scalar, dt, 
                E_plus_final, E_minus_final,
                B0, rho
            );

            // Write updated wave energies back to segment
            wave_data[0] = E_plus_final;
            wave_data[1] = E_minus_final;

            processed_segments++;
        }
    }

    if (rank == 0) {
        std::cout << "Wave-particle coupling updated for " << processed_segments
                  << " segments" << std::endl;
    }
}

// ============================================================================
// DIAGNOSTIC FUNCTIONS
// ============================================================================

double CalculateTotalWaveEnergyInSystem(PIC::Datum::cDatumStored& WaveEnergyDensity) {
    double local_total_energy = 0.0;

    // Sum energy across all segments assigned to this process
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();

        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;

            double* wave_data = segment->GetDatum_ptr(WaveEnergyDensity);
            if (wave_data) {
                local_total_energy += wave_data[0] + wave_data[1];  // E+ + E-
            }
        }
    }

    // Sum across all MPI processes
    double total_energy = 0.0;
    MPI_Allreduce(&local_total_energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return total_energy;  // [J]
}

double CalculateTotalParticleEnergyInSystem() {
    double local_total_energy = 0.0;

    // Sum energy across all segments assigned to this process
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();

        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;

            local_total_energy += CalculateTotalParticleEnergyInSegment(segment);
        }
    }

    // Sum across all MPI processes
    double total_energy = 0.0;
    MPI_Allreduce(&local_total_energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return total_energy;  // [J]
}

void CheckEnergyConservation(PIC::Datum::cDatumStored& WaveEnergyDensity,
                           bool verbose = false) {
    static double previous_total_energy = -1.0;

    double wave_energy = CalculateTotalWaveEnergyInSystem(WaveEnergyDensity);
    double particle_energy = CalculateTotalParticleEnergyInSystem();
    double total_energy = wave_energy + particle_energy;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        if (verbose) {
            std::cout << "Energy Conservation Check:" << std::endl;
            std::cout << "  Wave Energy:     " << wave_energy << " J" << std::endl;
            std::cout << "  Particle Energy: " << particle_energy << " J" << std::endl;
            std::cout << "  Total Energy:    " << total_energy << " J" << std::endl;
        }

        if (previous_total_energy > 0.0) {
            double energy_change = total_energy - previous_total_energy;
            double relative_change = energy_change / previous_total_energy;

            if (std::abs(relative_change) > 1.0e-6) {
                std::cout << "Warning: Energy conservation violation: "
                          << relative_change * 100.0 << "%" << std::endl;
            } else if (verbose) {
                std::cout << "  Energy Conservation: PASS (" << relative_change * 100.0 << "%)" << std::endl;
            }
        }

        previous_total_energy = total_energy;
    }
}

} // namespace IsotropicSEP
} // namespace AlfvenTurbulence_Kolmogorov  
} // namespace SEP

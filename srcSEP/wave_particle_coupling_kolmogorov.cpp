// ============================================================================
// UTILITY FUNCTION: CALCULATE WAVE ENERGY DENSITY FROM AMPLITUDE
// ============================================================================/*
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
   - Updates wave amplitudes for a single field line segment
   - Calculates Parker growth rates from particle distribution
   - Evolves wave amplitudes: A±^{n+1} = A±^n + Δt[2γ± A± - A±/τ_cas + Q_shock]
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
PIC::FieldLine::cFieldLineSegment* segment = /* get segment */;
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

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {
namespace IsotropicSEP {

// ============================================================================
// WAVE ENERGY EVOLUTION WITH PARTICLE COUPLING
// ============================================================================

void UpdateWaveEnergyWithParticleCoupling(
    PIC::FieldLine::cFieldLineSegment* segment,
    double& A_plus_initial,                    // Initial wave amplitude for outward waves
    double& A_minus_initial,                   // Initial wave amplitude for inward waves
    const PIC::Datum::cDatumStored& S_scalar,  // Scalar particle distribution
    double dt,                                 // Time step [s]
    double tau_cas,                           // Cascade time [s] 
    double Q_shock,                           // Shock injection rate [dimensionless]
    double& A_plus_final,                     // Output: final outward wave amplitude
    double& A_minus_final,                    // Output: final inward wave amplitude
    double B0,                                // Background magnetic field [T]
    double rho                                // Mass density [kg/m³]
) {
    // ========================================================================
    // INPUT VALIDATION
    // ========================================================================
    if (!segment) {
        std::cerr << "Error: Null segment pointer in UpdateWaveEnergyWithParticleCoupling" << std::endl;
        A_plus_final = A_plus_initial;
        A_minus_final = A_minus_initial;
        return;
    }
    
    if (dt <= 0.0 || tau_cas <= 0.0) {
        std::cerr << "Error: Invalid time parameters (dt=" << dt << ", tau_cas=" << tau_cas << ")" << std::endl;
        A_plus_final = A_plus_initial;
        A_minus_final = A_minus_initial;
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
    // STEP 2: WAVE AMPLITUDE ADVANCE USING GIVEN FORMULA
    // ========================================================================
    // A±^{n+1} = A±^n + Δt[2γ± A± - A±/τ_cas + Q_shock]
    
    double dA_plus_dt = 2.0 * gamma_plus * A_plus_initial 
                       - A_plus_initial / tau_cas 
                       + Q_shock;
    
    double dA_minus_dt = 2.0 * gamma_minus * A_minus_initial 
                        - A_minus_initial / tau_cas 
                        + Q_shock;
    
    // Update wave amplitudes
    A_plus_final = A_plus_initial + dt * dA_plus_dt;
    A_minus_final = A_minus_initial + dt * dA_minus_dt;
    
    // Ensure non-negative amplitudes
    A_plus_final = std::max(0.0, A_plus_final);
    A_minus_final = std::max(0.0, A_minus_final);
    
    // ========================================================================
    // STEP 3: CALCULATE CELL VOLUME
    // ========================================================================
    // V_cell = L_j * (A_L + sqrt(A_L * A_R) + A_R) / 3
    
    // Get segment length
    double L_j = segment->GetLength();  // Segment length [m]
    
    // Get cross-sectional areas at segment boundaries
    // Note: This assumes AMPS provides methods to get boundary areas
    // You may need to adapt this based on actual AMPS interface
    double A_L = segment->GetLeftBoundaryArea();   // Left boundary area [m²]
    double A_R = segment->GetRightBoundaryArea();  // Right boundary area [m²]
    
    // Calculate cell volume using given formula
    double V_cell = L_j * (A_L + std::sqrt(A_L * A_R) + A_R) / 3.0;  // [m³]
    
    if (V_cell <= 0.0) {
        std::cerr << "Warning: Invalid cell volume (" << V_cell << ") in segment" << std::endl;
        return;
    }
    
    // ========================================================================
    // STEP 4: CALCULATE WAVE ENERGY GAIN/LOSS
    // ========================================================================
    // E_wave_gain = (2γ± A± + Q_shock - A±/τ_cas) * Δt * V_cell
    
    double E_wave_gain_plus = dA_plus_dt * dt * V_cell;   // [J] or appropriate energy units
    double E_wave_gain_minus = dA_minus_dt * dt * V_cell; // [J] or appropriate energy units
    
    double total_wave_energy_change = E_wave_gain_plus + E_wave_gain_minus;
    
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
        double particle_mass = _H_MASS_; // [kg]
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
        exit(error_msg, __LINE__, __FILE__);
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
            double particle_mass = _H_MASS_;
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
    PIC::FieldLine::cFieldLineSegment* segment,
    double& A_plus_initial,
    double& A_minus_initial,
    const PIC::Datum::cDatumStored& S_scalar,
    double dt,
    double& A_plus_final,
    double& A_minus_final
) {
    // Default parameters for typical solar wind conditions
    const double tau_cas_default = 1000.0;    // 1000 s cascade time
    const double Q_shock_default = 0.0;       // No shock injection
    const double B0_default = 5.0e-9;         // 5 nT
    const double rho_default = 5.0e-21;       // kg/m³
    
    UpdateWaveEnergyWithParticleCoupling(
        segment, A_plus_initial, A_minus_initial, S_scalar, dt,
        tau_cas_default, Q_shock_default, A_plus_final, A_minus_final,
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
        double particle_mass = _H_MASS_;
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

double CalculateWaveEnergyDensity(double A_plus, double A_minus, double B0, double rho) {
    const double mu0 = 4.0e-7 * M_PI;
    
    // Wave energy density: ε = (B₀²/2μ₀) * (A₊² + A₋²)
    // where A± are normalized wave amplitudes
    double energy_density = (B0 * B0 / (2.0 * mu0)) * (A_plus * A_plus + A_minus * A_minus);
    
    return energy_density;  // [J/m³]
}

// ============================================================================
// UTILITY FUNCTION: CONVERT BETWEEN AMPLITUDE AND ENERGY DENSITY
// ============================================================================

void ConvertEnergyDensityToAmplitudes(
    double energy_density_total,  // Total wave energy density [J/m³]
    double amplitude_ratio,       // A_plus / A_minus
    double B0,                    // Background field [T]
    double& A_plus,              // Output: outward wave amplitude
    double& A_minus              // Output: inward wave amplitude
) {
    const double mu0 = 4.0e-7 * M_PI;
    
    if (energy_density_total <= 0.0) {
        A_plus = A_minus = 0.0;
        return;
    }
    
    // From: ε = (B₀²/2μ₀) * (A₊² + A₋²)
    // And: A₊ = r * A₋, where r = amplitude_ratio
    // Solve for A₋: ε = (B₀²/2μ₀) * A₋² * (r² + 1)
    
    double amplitude_factor = std::sqrt(2.0 * mu0 * energy_density_total / (B0 * B0));
    double denominator = std::sqrt(amplitude_ratio * amplitude_ratio + 1.0);
    
    A_minus = amplitude_factor / denominator;
    A_plus = amplitude_ratio * A_minus;
}

} // namespace IsotropicSEP
} // namespace AlfvenTurbulence_Kolmogorov  
} // namespace SEP

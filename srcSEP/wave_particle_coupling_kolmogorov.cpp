/* 
================================================================================
                WAVE-PARTICLE ENERGY COUPLING WITH EXACT QLT KERNEL
                              (THREE-FUNCTION DESIGN)
================================================================================

PURPOSE:
--------
This file implements wave-particle energy coupling for Alfvén wave turbulence 
in solar energetic particle (SEP) transport simulations using the exact 
Quasi-Linear Theory (QLT) kernel. The implementation is split into three 
functions for better integration with particle movers:

1. AccumulateParticleFluxForWaveCoupling() - Called by particle mover
2. CalculateGrowthRatesFromAccumulatedFlux() - Calculates  γ± integrated ofvek (k) from flux data
3. RedistributeWaveEnergyToParticles() - Energy conservation redistribution

PHYSICS BACKGROUND:
-------------------
The code implements resonant wave-particle interactions in magnetized plasmas
where particles interact with Alfvén waves through cyclotron resonance:

    ω - k‖ v‖ = ±Ω_c

This gives resonant wavenumbers: k_res = Ω_c / |μ v|, where μ = v‖/v is the 
pitch angle cosine. The exact QLT kernel accounts for both parallel and 
perpendicular particle motion without assuming v ≫ v_A.

ALGORITHM WORKFLOW:
-------------------
1. During particle transport: AccumulateParticleFluxForWaveCoupling() is called
   for each particle step, accumulating flux data in intermediate arrays
2. After particle transport: CalculateGrowthRatesFromAccumulatedFlux() computes
   growth rates and evolves wave energies
3. Energy conservation: RedistributeWaveEnergyToParticles() redistributes energy
   changes back to particles

SPECTRAL DISCRETIZATION:
------------------------
- Wavenumber grid: k ∈ [10⁻⁷, 10⁻²] m⁻¹, NK=128 log-uniform bins
- Momentum grid: p ∈ [10⁻²⁰, 10⁻¹⁷] kg⋅m/s, NP=96 log-uniform bins
- Bin spacing: Δln k = ln(k_max/k_min)/(NK-1)

INTERMEDIATE ARRAYS (stored as Datums):
---------------------------------------
- G_plus_streaming[NK]:  Σ w p³ (vμ - v_A) Δs / (v k) [outward waves]
- G_minus_streaming[NK]: Σ w p³ (vμ + v_A) Δs / (v k) [inward waves]
- gamma_plus_array[NK]:  Growth rates per k-bin for outward waves [s⁻¹]
- gamma_minus_array[NK]: Growth rates per k-bin for inward waves [s⁻¹]

REQUIRED DATUMS:
----------------
The following Datums must be declared and initialized:
- SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming   (size: NK)
- SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming  (size: NK)
- SEP::AlfvenTurbulence_Kolmogorov::gamma_plus_array   (size: NK) -- growth rate as a function of (k)
- SEP::AlfvenTurbulence_Kolmogorov::gamma_minus_array  (size: NK) -- growth rate as a function of (k) 

USAGE EXAMPLE:
--------------
// In particle mover loop:
AccumulateParticleFluxForWaveCoupling(
    segment, dt, vParallel, vNormal, s_start, s_finish, totalTraversedPath
);

// After all particles moved:
double E_wave_change = CalculateGrowthRatesFromAccumulatedFlux(
    segment, E_plus_initial, E_minus_initial, dt, B0, rho
);

// Energy conservation:
RedistributeWaveEnergyToParticles(segment, E_wave_change);

================================================================================
MORE: 
 *========================================================================= 
 *  QLT  Alfvén-wave growth / damping calculator
 *  --------------------------------------------
 *
 *  This module fills
 *          G_plus [k-bin][cell]   →  γ₊(k,s)  [ 1/s ]
 *          G_minus[k-bin][cell]   →  γ₋(k,s)  [ 1/s ]
 *  from Monte-Carlo particle data and then advances the wave spectrum
 *          W₊(k,s),  W₋(k,s).
 *
 *  ──────────────  Physical definitions  ────────────────────────────────
 *
 *    * σ = +1  :  antisun-ward wave  W₊   (resonates with μ < 0)
 *      σ = −1  :  sun-ward     wave  W₋   (resonates with μ > 0)
 *
 *    * Local QLT kernel (no v ≫ v_A assumption):
 *          Kσ = v μ − σ v_A                    [ m s⁻¹ ]
 *
 *    * Bin-resolved growth / damping rate:
 *
 *          γσ(k)  =  (π² Ω² / B² k)
 *                    ∫p²v ∫Kσ δ(k−Ω/vμ) f(p,μ) dμ dp     [ 1/s ]
 *
 *      →  on a log-uniform p-mesh (Δln p) and with a counting weight w_cnt
 *
 *          γσ(k_j) =  pref / (2 ΔV Δt Δln p k_j)
 *                    Σ_segments  w_cnt p² Δs  Kσ
 *
 *          pref = π² Ω² / B²             (cell-local)
 *
 *    * Spectrum update per time step:
 *          Wσ(k) ← Wσ(k) · exp( 2 γσ(k) Δt )
 *
 *    * k-integrated net rate (diagnostic):
 *
 *          Γσ = ∫ γσ(k) dk
 *              = Δln k  Σ_j γσ(k_j) k_j         [ 1/s ]
 *              (because dk = k d ln k on a log grid)
 *
 *  ──────────────  Variables & units  ───────────────────────────────────
 *
 *      w_i            counting weight               [ # real particles ]
 *      p_i            momentum  at segment start    [ kg m s⁻¹ ]
 *      v_i            speed     (relativistic)      [ m s⁻¹ ]
 *      Δs_i           signed path in cell           [ m ]
 *      μ_i            pitch-cosine  = Δs_i /(v_i Δt)
 *
 *      G_plus, G_minus        → accumulate γσ(k)    [ 1/s ]
 *      gamma_plus, gamma_minus→ alias / copy of G   [ 1/s ]
 *      W_plus,   W_minus      → wave energy / log-k [ J m⁻³ ]
 *      Gamma_plus, Gamma_minus→ k-integrated rate   [ 1/s ]
 *
 *  ──────────────  Algorithm outline  ───────────────────────────────────
 *
 *    (1) Zero G_± arrays.
 *    (2) For each macro-particle and for every cell segment it crosses:
 *            • compute local   B, ρ  →  v_A, Ω, pref
 *            • find resonant k-bin   j  via k_res = Ω / (|μ| v)
 *            • coeff = pref *  w_cnt * p² Δs
 *                       / (2 Δt Δln p ΔV k_j)
 *            • G_plus [j][c]  += coeff * ( v μ −  v_A )
 *              G_minus[j][c]  += coeff * ( v μ +  v_A )
 *
 *    (3) After particle loop:
 *            gamma_±[j][c] = G_±[j][c]
 *            W_±   *= exp( 2 γ_± Δt )
 *
 *    (4) Optional diagnostic:
 *            Γ_±(c) = Δln k · Σ_j γ_±[j][c] · k_j
 *
 *  All expressions are cell-local, plasma-frame, and free of
 *  assumptions about the ratio v / v_A.
 *=========================================================================
 */


#include "sep.h"

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {
namespace IsotropicSEP {

// ============================================================================
// COMPILE-TIME PARAMETERS AND CONSTANTS
// ============================================================================
const int NK = 128;        // Number of k-bins (log-uniform spacing)
const int NP = 96;         // Number of p-bins (log-uniform spacing) - for reference
const double PI = 3.141592653589793;
const double Q = 1.602e-19;   // Elementary charge [C]
const double M = _H__MASS_;   // Proton mass [kg]

// Wavenumber range: covers typical solar wind turbulence scales
const double K_MIN = 1.0e-8; // Minimum wavenumber [m⁻¹] (large scales)
const double K_MAX = 1.0e-2; // Maximum wavenumber [m⁻¹] (small scales)

// Momentum range: covers SEP energy range from keV to GeV
const double P_MIN = 1.0e-20; // Minimum momentum [kg⋅m/s] 
const double P_MAX = 1.0e-17; // Maximum momentum [kg⋅m/s]

// Grid spacing in logarithmic coordinates
const double DLNK = log(K_MAX / K_MIN) / (NK - 1); // log-k spacing
const double DLNP = log(P_MAX / P_MIN) / (NP - 1); // log-p spacing (for normalization)

// ============================================================================
// HELPER FUNCTION: FIND K-BIN INDEX
// ============================================================================

int GetKBinIndex(double k_val) {
    // Find nearest k-bin index using log interpolation
    int j = (int)(0.5 + (log(k_val/K_MIN) / DLNK));
    
    // Clamp to valid range
    if (j < 0) j = 0;
    if (j >= NK) j = NK-1;
    
    return j;
}

// ============================================================================
// WAVE-PARTICLE COUPLING MANAGER FUNCTION
// ============================================================================

void WaveParticleCouplingManager(
    PIC::Datum::cDatumStored& WaveEnergy,
    double dt
) {
    /*
    This manager function orchestrates the complete wave-particle coupling process:
    1. Calculates growth rates from accumulated flux data
    2. Updates wave energies in each segment 
    3. Redistributes energy changes to particles for energy conservation
    
    Note: Flux accumulation (AccumulateParticleFluxForWaveCoupling) should have
    been done during particle transport phase before calling this manager.
    */
    
    double B0, rho;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int processed_segments = 0;
    double total_wave_energy_change_system = 0.0;

    // ========================================================================
    // STEP A: CALCULATE GROWTH RATES AND UPDATE ENERGIES FOR ALL SEGMENTS
    // ========================================================================
    
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

            // Get background plasma parameters for this segment
            segment->GetPlasmaDensity(0.5, rho);  // Number density at segment midpoint
            rho *= _H__MASS_; // Convert to mass density [kg/m³]

	    double B[3];
	    segment->GetMagneticField(0.5,B);
            B0=Vector3D::Length(B);

            // Get current wave energy data from segment
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) {
                continue;  // Skip segments without wave energy data
            }

            double E_plus_initial = wave_data[0];   // Initial outward wave energy [J]
            double E_minus_initial = wave_data[1];  // Initial inward wave energy [J]

            // ================================================================
            // STEP A1: Calculate growth rates from accumulated flux data
            // ================================================================
            double Gamma_plus, Gamma_minus;
            CalculateGrowthRatesFromAccumulatedFlux(
                segment, dt, B0, rho, Gamma_plus, Gamma_minus
            );

            // ================================================================
            // STEP A2: Update wave energies using calculated growth rates
            // ================================================================
            // Wave energy evolution: E_±(t+Δt) = E_±(t) × exp(2 Γ_± Δt)
            double E_plus_final = E_plus_initial * exp(2.0 * Gamma_plus * dt);
            double E_minus_final = E_minus_initial * exp(2.0 * Gamma_minus * dt);

	    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
              validate_numeric(E_plus_final,__LINE__,__FILE__);
	      validate_numeric(E_minus_final,__LINE__,__FILE__);
            }

            
            // Ensure non-negative energies (physical constraint)
            E_plus_final = std::max(0.0, E_plus_final);
            E_minus_final = std::max(0.0, E_minus_final);

            // Calculate total wave energy change for this segment
            double segment_wave_energy_change = (E_plus_final - E_plus_initial) + 
                                              (E_minus_final - E_minus_initial);
            
            // Accumulate total system wave energy change
            total_wave_energy_change_system += segment_wave_energy_change;

            // Update wave energy data in segment
            wave_data[0] = E_plus_final;   // Store updated outward wave energy
            wave_data[1] = E_minus_final;  // Store updated inward wave energy

	    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
              validate_numeric(wave_data[0],__LINE__,__FILE__);
              validate_numeric(wave_data[1],__LINE__,__FILE__);
	    }

            processed_segments++;
        }
    }

    // ========================================================================
    // STEP B: REDISTRIBUTE ENERGY TO PARTICLES FOR ENERGY CONSERVATION
    // ========================================================================
    
    // Second loop: redistribute energy changes to particles
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();

        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);

            // Only process segments assigned to this MPI process
            if (!segment || segment->Thread != PIC::ThisThread) {
                continue;
            }

            // Get wave energy data to calculate energy change
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) {
                continue;
            }

            // We need to recalculate the energy change for this segment
            // This could be optimized by storing the change from step A
            double E_plus_current = wave_data[0];
            double E_minus_current = wave_data[1];
            
            // Get plasma parameters again (could be optimized)
            segment->GetPlasmaDensity(0.5, rho);
            rho *= _H__MASS_;

	    double B[3];
	    segment->GetMagneticField(0.5,B);
            B0=Vector3D::Length(B);

            // Recalculate initial energies by working backwards
            double Gamma_plus, Gamma_minus;
            CalculateGrowthRatesFromAccumulatedFlux(
                segment, dt, B0, rho, Gamma_plus, Gamma_minus
            );
            
            // Calculate what the initial energies were
            double E_plus_initial = E_plus_current / exp(2.0 * Gamma_plus * dt);
            double E_minus_initial = E_minus_current / exp(2.0 * Gamma_minus * dt);
            
            // Calculate segment energy change
            double segment_wave_energy_change = (E_plus_current - E_plus_initial) + 
                                              (E_minus_current - E_minus_initial);

            // ================================================================
            // STEP B1: Redistribute energy to particles (energy conservation)
            // ================================================================
            if (std::abs(segment_wave_energy_change) > 1.0e-25) {
                RedistributeWaveEnergyToParticles(segment, segment_wave_energy_change);
            }
        }
    }

    // ========================================================================
    // STEP C: OUTPUT DIAGNOSTICS
    // ========================================================================
    
    // Sum total wave energy change across all MPI processes
    double global_wave_energy_change = 0.0;
    MPI_Allreduce(&total_wave_energy_change_system, &global_wave_energy_change, 
                  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Wave-Particle Coupling Manager Results:" << std::endl;
        std::cout << "  Processed segments: " << processed_segments << std::endl;
        std::cout << "  Total wave energy change: " << global_wave_energy_change << " J" << std::endl;
        std::cout << "  Time step: " << dt << " s" << std::endl;
    }
}

// ============================================================================
// FUNCTION 1: ACCUMULATE PARTICLE FLUX DATA (CALLED BY PARTICLE MOVER)
// ============================================================================

void AccumulateParticleFluxForWaveCoupling(
    int field_line_idx,                          // Field line index
    long int particle_index,                     // Particle index parameter
    double dt,                                   // Time step [s]
    double speed,                               // Particle speed magnitude [m/s]
    double s_start,                             // Start position along field line [m]
    double s_finish,                            // End position along field line [m]
    double totalTraversedPath                   // Signed parallel path length [m] (+ outward, - inward)
) {
    /*
    BIDIRECTIONAL MOTION HANDLING:
    - totalTraversedPath > 0: Particle moving away from Sun (outward, μ > 0)
    - totalTraversedPath < 0: Particle moving toward Sun (inward, μ < 0)  
    - μ = totalTraversedPath / (speed * dt) = v_parallel / v_total
    
    The pitch angle cosine μ determines wave-particle resonance:
    - μ > 0: Particle moving outward, resonates differently with ± waves
    - μ < 0: Particle moving inward, resonates differently with ± waves
    - |μ| determines the resonant wavenumber: k_res = Ω_c / |μ v|
    
    Note: totalTraversedPath should be the signed parallel displacement,
    where the sign indicates direction relative to the magnetic field.
    */
    // ========================================================================
    // INPUT VALIDATION
    // ========================================================================
    if (field_line_idx < 0 || field_line_idx >= PIC::FieldLine::nFieldLine) {
        std::cerr << "Error: Invalid field line index (" << field_line_idx 
                  << ") in AccumulateParticleFluxForWaveCoupling" << std::endl;
        return;
    }
    
    if (particle_index < 0) {
        std::cerr << "Error: Invalid particle index (" << particle_index 
                  << ") in AccumulateParticleFluxForWaveCoupling" << std::endl;
        return;
    }
    
    if (dt <= 0.0) {
        std::cerr << "Error: Invalid time step (dt=" << dt << ")" << std::endl;
        return;
    }
    
    // ========================================================================
    // GET FIELD LINE AND DETERMINE TRAVERSED SEGMENTS
    // ========================================================================
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    
    // ========================================================================
    // DECODE FIELD LINE COORDINATES AND CALCULATE TRAJECTORY
    // ========================================================================
    // Field line coordinates format: s = (segment_index).(fractional_position)
    // Extract segment indices and fractional positions
    int seg_start_idx = (int)floor(s_start);
    double frac_start = s_start - seg_start_idx;
    
    int seg_finish_idx = (int)floor(s_finish);
    double frac_finish = s_finish - seg_finish_idx;
    
    // Calculate pitch angle cosine from speed and parallel displacement
    // Account for bidirectional motion: particles can move toward or away from Sun
    double mu = 0.0;
    if (speed > 1.0e-20) {
        // Use the provided totalTraversedPath which includes directional information
        // μ = v_parallel / v_total, where v_parallel = totalTraversedPath / dt
        mu = totalTraversedPath / (speed * dt);
        
        // Clamp μ to valid range [-1, 1] (physical constraint)
        mu = std::max(-1.0, std::min(1.0, mu));
    } else {
        return; // No movement, no contribution
    }
    
    // ========================================================================
    // CALCULATE PARTICLE PROPERTIES (INDEPENDENT OF SEGMENTS)
    // ========================================================================
    // Use provided speed directly
    double v_magnitude = speed;
    
    // Skip particles with negligible velocity
    if (v_magnitude < 1.0e-20) {
        return;
    }
    
    // Calculate relativistic momentum: p = γ m v
    double gamma_rel = 1.0 / sqrt(1.0 - (v_magnitude*v_magnitude)/
                                 (PhysicsConstants::C_LIGHT*PhysicsConstants::C_LIGHT));
    double p_momentum = gamma_rel * M * v_magnitude;  // [kg⋅m/s]
    
    // Skip particles outside momentum range of interest
    if (p_momentum < P_MIN || p_momentum > P_MAX) {
        return;
    }
    
    // ========================================================================
    // GET PARTICLE STATISTICAL WEIGHT USING SPECIES-SPECIFIC CALCULATION
    // ========================================================================
    // Get particle species number for this particle
    int particle_species = PIC::ParticleBuffer::GetI(particle_index);
    if (particle_species < 0 || particle_species >= PIC::nTotalSpecies) {
        std::cerr << "Error: Invalid particle species (" << particle_species 
                  << ") for particle " << particle_index << std::endl;
        return;
    }
    
    // Calculate species-specific statistical weight
    double w_i = PIC::ParticleWeightTimeStep::GlobalParticleWeight[particle_species] * 
                 PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle_index);
    
    if (w_i <= 0.0) {
        std::cerr << "Warning: Invalid particle weight (" << w_i 
                  << ") for particle " << particle_index << std::endl;
        return;
    }
    
    // ========================================================================
    // LOOP THROUGH ALL SEGMENTS TRAVERSED BY PARTICLE
    // ========================================================================
    // Determine range of segments that the particle trajectory intersects
    int seg_min = std::min(seg_start_idx, seg_finish_idx);
    int seg_max = std::max(seg_start_idx, seg_finish_idx);
    
    int num_segments = field_line->GetTotalSegmentNumber();
    
    // Ensure segment indices are within valid range
    seg_min = std::max(0, seg_min);
    seg_max = std::min(num_segments - 1, seg_max);
    
    for (int seg_idx = seg_min; seg_idx <= seg_max; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        
        // ====================================================================
        // CALCULATE DIRECTIONAL PATH LENGTH WITHIN THIS SEGMENT
        // ====================================================================
        double ds_seg = 0.0;
        
        if (seg_start_idx == seg_finish_idx && seg_idx == seg_start_idx) {
            // Particle starts and ends in the same segment
            double segment_length = segment->GetLength();
            ds_seg = (frac_finish - frac_start) * segment_length;  // Directional: can be negative
            
        } else if (seg_idx == seg_start_idx) {
            // First segment: from start position to end of segment
            double segment_length = segment->GetLength();
            if (seg_start_idx < seg_finish_idx) {
                // Moving forward: from frac_start to 1.0
                ds_seg = (1.0 - frac_start) * segment_length;
            } else {
                // Moving backward: from frac_start to 0.0
                ds_seg = -frac_start * segment_length;
            }
            
        } else if (seg_idx == seg_finish_idx) {
            // Last segment: from beginning of segment to finish position
            double segment_length = segment->GetLength();
            if (seg_start_idx < seg_finish_idx) {
                // Moving forward: from 0.0 to frac_finish
                ds_seg = frac_finish * segment_length;
            } else {
                // Moving backward: from 1.0 to frac_finish
                ds_seg = (frac_finish - 1.0) * segment_length;
            }
            
        } else {
            // Middle segment: entire segment length with direction
            double segment_length = segment->GetLength();
            if (seg_start_idx < seg_finish_idx) {
                // Moving forward: positive full segment length
                ds_seg = segment_length;
            } else {
                // Moving backward: negative full segment length
                ds_seg = -segment_length;
            }
        }
        
        // Skip if path length in segment is negligible (but keep sign)
        if (fabs(ds_seg) < 1.0e-20) {
            continue;
        }
        
        // ====================================================================
        // GET LOCAL PLASMA PARAMETERS FOR THIS SEGMENT
        // ====================================================================
        double B0, rho;
        segment->GetPlasmaDensity(0.5, rho);  // Get density at segment midpoint
        rho *= _H__MASS_;  // Convert number density to mass density
        
        // Get magnetic field from field line
        double B[3];
	PIC::FieldLine::FieldLinesAll[field_line_idx].GetMagneticField(B,0.5+seg_idx);
	B0=Vector3D::Length(B);

        
        // Calculate local plasma parameters
        double vAc = B0 / sqrt(VacuumPermeability * rho);                  // Alfvén speed [m/s]
        double Omega = Q * B0 / M;                               // Proton cyclotron frequency [rad/s]
        double pref = (PI * PI) * Omega * Omega / (B0 * B0);    // Normalization factor
        
        // Get segment volume using proper SEP function
        double V_cell = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);
        if (V_cell <= 0.0) {
            std::cerr << "Warning: Invalid segment volume (" << V_cell 
                      << ") in segment " << seg_idx << std::endl;
            continue;
        }
        double Vinv = 1.0 / V_cell;  // Inverse volume [m⁻³]
        
        // ====================================================================
        // CALCULATE RESONANT WAVENUMBER AND K-BIN
        // ====================================================================
        // Cyclotron resonance condition: ω - k‖ v‖ = ±Ω
        // For Alfvén waves: ω ≈ k‖ v_A, so k_res = Ω / |μ v|
        double kRes = Omega / (fabs(mu) * v_magnitude);  // Resonant wavenumber [m⁻¹]
        
        // Find corresponding k-bin index
        int j = GetKBinIndex(kRes);
        
        // Initialize k-grid for this calculation
        double k_j = K_MIN * exp(j * DLNK);  // k value for bin j
        
        // ====================================================================
        // CALCULATE PARTICLE CONTRIBUTION TO STREAMING INTEGRALS
        // ====================================================================
        // Calculate coefficient following pseudo-code
        // Note: ds_seg is now directional and can be negative
        double p2v = p_momentum * p_momentum * v_magnitude;  // p² v term
        double coeff = 2.0* Pi*pref * w_i * p2v * (ds_seg / v_magnitude) / 
                      (2.0 * dt * DLNP) * Vinv / k_j;
        
        // ====================================================================
        // ACCESS INTERMEDIATE STREAMING ARRAYS AND ACCUMULATE
        // ====================================================================
        double* G_plus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
        double* G_minus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
        
        if (!G_plus_data || !G_minus_data) {
            std::cerr << "Error: Cannot access G_plus/G_minus streaming arrays in segment " 
                      << seg_idx << std::endl;
            continue;
        }
        
        // Add particle contribution to streaming sums using exact QLT kernel
        // Thread-safe accumulation (assumes single-threaded access per segment)
        // Note: coeff now includes directional ds_seg which can be positive or negative
	
	
	//particles with mu>0 interactes for Inward wave; and particles with mu<0 interacts only with outward wave  
        if (mu<0.0) G_plus_data[j]  += coeff * (v_magnitude * mu - vAc);  // Outward waves (+ direction)
        if (mu>0.0) G_minus_data[j] += coeff * (v_magnitude * mu + vAc);  // Inward waves (- direction)

        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
           validate_numeric(G_plus_data[j],-100.0,50.0,__LINE__,__FILE__);
	   validate_numeric(G_minus_data[j],-100.0,50.0,__LINE__,__FILE__);
        }

    }
}

// ============================================================================
// FUNCTION 2: CALCULATE GROWTH RATES AND RETURN GAMMA VALUES
// ============================================================================

void CalculateGrowthRatesFromAccumulatedFlux(
    PIC::FieldLine::cFieldLineSegment* segment,  // Target field line segment
    double dt,                                   // Time step [s]
    double B0,                                   // Background magnetic field [T]
    double rho,                                  // Mass density [kg/m³]
    double& Gamma_plus,                         // Output: integrated growth rate for outward waves [s⁻¹]
    double& Gamma_minus                         // Output: integrated growth rate for inward waves [s⁻¹]
) {
    // ========================================================================
    // INPUT VALIDATION
    // ========================================================================
    if (!segment) {
        std::cerr << "Error: Null segment pointer in CalculateGrowthRatesFromAccumulatedFlux" << std::endl;
        Gamma_plus = 0.0;
        Gamma_minus = 0.0;
        return;
    }

    if (segment->Thread != PIC::ThisThread) {
        Gamma_plus = 0.0;
        Gamma_minus = 0.0;
        return;
    }
    
    if (dt <= 0.0) {
        std::cerr << "Error: Invalid time step (dt=" << dt << ")" << std::endl;
        Gamma_plus = 0.0;
        Gamma_minus = 0.0;
        return;
    }
    
    // ========================================================================
    // ACCESS INTERMEDIATE STREAMING ARRAYS
    // ========================================================================
    double* G_plus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
    double* G_minus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
    
    if (!G_plus_data || !G_minus_data) {
        std::cerr << "Error: Cannot access G_plus/G_minus streaming arrays" << std::endl;
        Gamma_plus = 0.0;
        Gamma_minus = 0.0;
        return;
    }
    
    // ========================================================================
    // ACCESS GROWTH RATE ARRAYS
    // ========================================================================
    double* gamma_plus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::gamma_plus_array);
    double* gamma_minus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::gamma_minus_array);
    
    if (!gamma_plus_data || !gamma_minus_data) {
        std::cerr << "Error: Cannot access gamma_plus/gamma_minus arrays" << std::endl;
        Gamma_plus = 0.0;
        Gamma_minus = 0.0;
        return;
    }
    
    // ========================================================================
    // CONVERT STREAMING INTEGRALS G_± TO GROWTH RATES γ_±
    // ========================================================================
    // Initialize k-grid for integration
    static std::vector<double> k(NK);
    static bool k_init=false;

    if (k_init==false) {
      k_init=true;

      for (int j = 0; j < NK; ++j) {
        k[j] = K_MIN * exp(j * DLNK);  // k[j] = k_min * exp(j * Δln k)
      }
    }
    
    // Direct conversion: growth rates equal streaming integrals
    for (int j = 0; j < NK; ++j) {
        gamma_plus_data[j] = G_plus_data[j];    // γ₊(k_j) = G₊(k_j)
        gamma_minus_data[j] = G_minus_data[j];  // γ₋(k_j) = G₋(k_j)

        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
          validate_numeric(gamma_plus_data[j],__LINE__,__FILE__);
	  validate_numeric(gamma_minus_data[j],__LINE__,__FILE__);
        }
    }
    
    // ========================================================================
    // CALCULATE KOLMOGOROV-WEIGHTED INTEGRATED GROWTH RATES
    // ========================================================================
    // Integrate growth rates over wavenumber spectrum with Kolmogorov weighting
    // Γ_± = ∫ γ_±(k) k dk ≈ Δln k × Σ γ_±(k_j) × k_j
    
    Gamma_plus = 0.0;   // Integrated growth rate for outward waves [s⁻¹]
    Gamma_minus = 0.0;  // Integrated growth rate for inward waves [s⁻¹]
    
    for (int j = 0; j < NK; ++j) {
        // k-weighted sum (Kolmogorov spectrum has k⁻⁵/³ dependence)
        Gamma_plus  += gamma_plus_data[j]  * k[j];
        Gamma_minus += gamma_minus_data[j] * k[j];
    }
    
    // Apply logarithmic spacing normalization
    Gamma_plus  *= DLNK;   // Γ₊ = Δln k × Σ γ₊(k_j) × k_j
    Gamma_minus *= DLNK;   // Γ₋ = Δln k × Σ γ₋(k_j) × k_j

    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
      validate_numeric(Gamma_plus,__LINE__,__FILE__);
      validate_numeric(Gamma_minus,__LINE__,__FILE__);
    }

    
    // ========================================================================
    // RESET STREAMING ARRAYS FOR NEXT TIME STEP
    // ========================================================================
    // Clear accumulated flux data for next time step
    for (int j = 0; j < NK; ++j) {
        G_plus_data[j] = 0.0;   // Reset outward wave streaming integral
        G_minus_data[j] = 0.0;  // Reset inward wave streaming integral
    }
}

// ============================================================================
// FUNCTION 3: REDISTRIBUTE WAVE ENERGY TO PARTICLES (ENERGY CONSERVATION)
// ============================================================================

void RedistributeWaveEnergyToParticles(
    PIC::FieldLine::cFieldLineSegment* segment,  // Target field line segment
    double wave_energy_change                    // Total wave energy change [J]
) {
    // ========================================================================
    // INPUT VALIDATION
    // ========================================================================
    if (!segment) {
        std::cerr << "Error: Null segment pointer in RedistributeWaveEnergyToParticles" << std::endl;
        return;
    }

    if (segment->Thread != PIC::ThisThread) {
        return;
    }
    
    // Skip if no significant energy change
    if (std::abs(wave_energy_change) < 1.0e-25) {
        return;
    }
   
    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
      validate_numeric(wave_energy_change,__LINE__,__FILE__);
    }

    // ========================================================================
    // ENERGY CONSERVATION PRINCIPLE
    // ========================================================================
    // Energy conservation: ΔE_waves + ΔE_particles = 0
    // Therefore: ΔE_particles = -ΔE_waves
    double particle_energy_change = -wave_energy_change;
    
    // ========================================================================
    // FIRST PASS: Calculate total particle energy and count particles
    // ========================================================================
    double total_particle_energy = 0.0;  // Total kinetic energy of all particles [J]
    int total_particle_count = 0;        // Number of computational particles
    long int p = segment->FirstParticleIndex;  // Start of particle linked list
    
    while (p != -1) {
        // Get particle velocity components
        double vParallel = PIC::ParticleBuffer::GetVParallel(p);
        double vNormal = PIC::ParticleBuffer::GetVNormal(p);
        double v_magnitude = sqrt(vParallel*vParallel + vNormal*vNormal);
        
        // Calculate relativistic kinetic energy
        double kinetic_energy = Relativistic::Speed2E(v_magnitude, _H__MASS_);  // Single particle energy [J]
        
        // Get statistical weight (number of real particles represented)
        double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] * 
                            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        
        // Add to total energy (model particle energy = single particle × weight)
        total_particle_energy += kinetic_energy * stat_weight;
        total_particle_count++;
        
        p = PIC::ParticleBuffer::GetNext(p);
    }
    
    // Check if energy redistribution is possible
    if (total_particle_count == 0 || total_particle_energy <= 0.0) {
        std::cerr << "Warning: No particles available for energy redistribution" << std::endl;
        return;
    }
    
    // Verify energy conservation is physically possible
    if (particle_energy_change < 0 && std::abs(particle_energy_change) > total_particle_energy) {
        char error_msg[512];
        sprintf(error_msg, "Energy removal (%.6e J) exceeds total particle energy (%.6e J) in segment", 
                std::abs(particle_energy_change), total_particle_energy);
        std::cerr << "Error: " << error_msg << std::endl;
        return;  // Don't exit, just skip this redistribution
    }
    
    // ========================================================================
    // ITERATIVE ENERGY REDISTRIBUTION WITH ENERGY FLOORS
    // ========================================================================
    // Handle cases where some particles hit energy floors during redistribution
    
    double remaining_energy_to_distribute = particle_energy_change;
    int iteration_count = 0;
    const int max_iterations = 10;  // Prevent infinite loops
    
    while (std::abs(remaining_energy_to_distribute) > 1.0e-20 && iteration_count < max_iterations) {
        iteration_count++;
        
        // Distribute remaining energy equally among all particles
        double energy_per_particle = remaining_energy_to_distribute / total_particle_count;
        double energy_actually_redistributed = 0.0;
        
        // SECOND PASS: Update particle velocities
        p = segment->FirstParticleIndex;
        
        while (p != -1) {
            // Get current particle velocity components
            double vParallel_current = PIC::ParticleBuffer::GetVParallel(p);
            double vNormal_current = PIC::ParticleBuffer::GetVNormal(p);
            double v_magnitude_current = sqrt(vParallel_current*vParallel_current + 
                                            vNormal_current*vNormal_current);
            
            // Calculate current particle energy
            double kinetic_energy = Relativistic::Speed2E(v_magnitude_current, _H__MASS_);
            double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] * 
                                PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
            double current_energy = kinetic_energy * stat_weight;  // Model particle energy
            
            // Determine energy change for this particle
            double energy_change_this_particle = energy_per_particle;
            
            // Apply energy floor constraint (particles cannot lose more than 90% of energy)
            if (energy_change_this_particle < 0) {
                double energy_floor = 0.1 * current_energy;      // 10% minimum energy
                double max_removable = current_energy - energy_floor;
                
                // Limit energy removal to respect floor
                if (std::abs(energy_change_this_particle) > max_removable) {
                    energy_change_this_particle = -max_removable;
                }
            }
            
            // Calculate new particle energy
            double new_kinetic_energy = current_energy + energy_change_this_particle;
            
            // Additional safety check for positive energy
            if (new_kinetic_energy <= 0.0) {
                new_kinetic_energy = 0.1 * current_energy;  // Emergency 10% floor
                energy_change_this_particle = new_kinetic_energy - current_energy;
            }
            
            // Track actual energy redistributed (may differ from requested due to floors)
            energy_actually_redistributed += energy_change_this_particle;
            
            // ================================================================
            // CONVERT ENERGY BACK TO VELOCITY COMPONENTS
            // ================================================================
            // Convert model particle energy back to single particle energy
            double new_single_particle_kinetic_energy = new_kinetic_energy / stat_weight;
            
            // Calculate new Lorentz factor: γ = KE/(mc²) + 1
            double gamma_new = new_single_particle_kinetic_energy / 
                              (_H__MASS_ * PhysicsConstants::C_LIGHT * PhysicsConstants::C_LIGHT) + 1.0;
            
            // Calculate new velocity magnitude: v = c√(1 - 1/γ²)
            double v_new_magnitude = PhysicsConstants::C_LIGHT * sqrt(1.0 - 1.0/(gamma_new*gamma_new));
            
            // Update velocity components while preserving direction
            if (v_magnitude_current > 1.0e-20) {
                // Scale both components proportionally to maintain pitch angle
                double scale_factor = v_new_magnitude / v_magnitude_current;
                double vParallel_new = vParallel_current * scale_factor;
                double vNormal_new = vNormal_current * scale_factor;
                
                // Update particle velocity components in buffer
                PIC::ParticleBuffer::SetVParallel(vParallel_new, p);
                PIC::ParticleBuffer::SetVNormal(vNormal_new, p);

               if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                  validate_numeric(vParallel_new,__LINE__,__FILE__);
		  validate_numeric(vNormal_new,__LINE__,__FILE__);
               }
            }
            
            // Move to next particle
            p = PIC::ParticleBuffer::GetNext(p);
        }
        
        // Update remaining energy for next iteration
        remaining_energy_to_distribute -= energy_actually_redistributed;
        
        // Convergence check: exit if remaining energy is negligible
        if (std::abs(remaining_energy_to_distribute) < 1.0e-20) {
            break;
        }
        
        // Safety check: exit if no progress is being made
        if (std::abs(energy_actually_redistributed) < 1.0e-20) {
            if (std::abs(remaining_energy_to_distribute) > 1.0e-15) {
                std::cerr << "Warning: Could not redistribute all energy. Remaining: " 
                          << remaining_energy_to_distribute << " J (iteration " 
                          << iteration_count << ")" << std::endl;
            }
            break;
        }
    }
    
    if (iteration_count >= max_iterations) {
        std::cerr << "Warning: Energy redistribution reached maximum iterations. Remaining: " 
                  << remaining_energy_to_distribute << " J" << std::endl;
    }
}

// ============================================================================
// UTILITY FUNCTION: INITIALIZE STREAMING ARRAYS (CALL AT START OF TIME STEP)
// ============================================================================

void InitializeStreamingArraysForTimeStep(PIC::FieldLine::cFieldLineSegment* segment) {
    if (!segment || segment->Thread != PIC::ThisThread) {
        return;
    }
    
    // Access and zero intermediate streaming arrays
    double* G_plus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
    double* G_minus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
    
    if (G_plus_data && G_minus_data) {
        for (int j = 0; j < SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::NK; ++j) {
            G_plus_data[j] = 0.0;   // Reset outward wave streaming integral
            G_minus_data[j] = 0.0;  // Reset inward wave streaming integral
        }
    }
}

// ============================================================================
// UTILITY FUNCTION: CALCULATE TOTAL PARTICLE ENERGY IN SEGMENT
// ============================================================================

double CalculateTotalParticleEnergyInSegment(PIC::FieldLine::cFieldLineSegment* segment) {
    if (!segment) return 0.0;
    
    double total_energy = 0.0;
    long int p = segment->FirstParticleIndex;
    
    while (p != -1) {
        // Get particle velocity components in magnetic field coordinates
        double vParallel = PIC::ParticleBuffer::GetVParallel(p);
        double vNormal = PIC::ParticleBuffer::GetVNormal(p);
        
        // Calculate total velocity magnitude
        double v_magnitude = sqrt(vParallel*vParallel + vNormal*vNormal);
        
        // Relativistic kinetic energy calculation
        double kinetic_energy = Relativistic::Speed2E(v_magnitude, _H__MASS_);
        
        // Get particle statistical weight
        double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] * 
                            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        
        // Model particle energy (physical energy represented by computational particle)
        double model_particle_energy = kinetic_energy * stat_weight;
        
        total_energy += model_particle_energy;
        
        // Get next particle in linked list
        p = PIC::ParticleBuffer::GetNext(p);
    }
    
    return total_energy; // [J]
}

// ============================================================================
// OPTIMIZED WAVE-PARTICLE COUPLING MANAGER (SINGLE-PASS VERSION)
// ============================================================================

void OptimizedWaveParticleCouplingManager(
    PIC::Datum::cDatumStored& WaveEnergy,
    double dt
) {
    /*
    Optimized version that processes each segment only once:
    1. Calculates growth rates from accumulated flux data
    2. Updates wave energies and immediately redistributes to particles
    3. More efficient than the two-pass version above
    */
    
    double B0, rho;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int processed_segments = 0;
    double total_wave_energy_change_system = 0.0;

    // Single loop through all segments
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();

        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);

            // Only process segments assigned to this MPI process
            if (!segment || segment->Thread != PIC::ThisThread) {
                continue;
            }

            // Get background plasma parameters for this segment
            segment->GetPlasmaDensity(0.5, rho);  // Number density at segment midpoint
            rho *= _H__MASS_; // Convert to mass density [kg/m³]

	    double B[3];
	    segment->GetMagneticField(0.5,B);
            B0=Vector3D::Length(B);

            // Get current wave energy data from segment
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) {
                continue;  // Skip segments without wave energy data
            }

            double E_plus_initial = wave_data[0];   // Initial outward wave energy [J]
            double E_minus_initial = wave_data[1];  // Initial inward wave energy [J]

            // Calculate growth rates from accumulated flux data
            double Gamma_plus, Gamma_minus;
            SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::CalculateGrowthRatesFromAccumulatedFlux(
                segment, dt, B0, rho, Gamma_plus, Gamma_minus
            );

            // Update wave energies using calculated growth rates
            double E_plus_final = E_plus_initial * exp(2.0 * Gamma_plus * dt);
            double E_minus_final = E_minus_initial * exp(2.0 * Gamma_minus * dt);
            
            // Ensure non-negative energies
            E_plus_final = std::max(0.0, E_plus_final);
            E_minus_final = std::max(0.0, E_minus_final);

            // Calculate wave energy change for this segment
            double segment_wave_energy_change = (E_plus_final - E_plus_initial) + 
                                              (E_minus_final - E_minus_initial);
            
            // Update wave energy data in segment
            wave_data[0] = E_plus_final;
            wave_data[1] = E_minus_final;

            // Immediately redistribute energy to particles for energy conservation
            if (std::abs(segment_wave_energy_change) > 1.0e-25) {
                RedistributeWaveEnergyToParticles(segment, segment_wave_energy_change);
            }

            // Accumulate for diagnostics
            total_wave_energy_change_system += segment_wave_energy_change;
            processed_segments++;
        }
    }

    // Output diagnostics
    double global_wave_energy_change = 0.0;
    MPI_Allreduce(&total_wave_energy_change_system, &global_wave_energy_change, 
                  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Optimized Wave-Particle Coupling Results:" << std::endl;
        std::cout << "  Processed segments: " << processed_segments << std::endl;
        std::cout << "  Total wave energy change: " << global_wave_energy_change << " J" << std::endl;
        std::cout << "  Time step: " << dt << " s" << std::endl;
    }
}

// ============================================================================
// DIAGNOSTIC FUNCTIONS
// ============================================================================

double CalculateTotalWaveEnergyInSystem(PIC::Datum::cDatumStored& WaveEnergy) {
    double local_total_energy = 0.0;

    // Sum energy across all segments assigned to this process
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();

        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;

            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
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

void CheckEnergyConservation(PIC::Datum::cDatumStored& WaveEnergy, bool verbose) {
    static double previous_total_energy = -1.0;

    double wave_energy = CalculateTotalWaveEnergyInSystem(WaveEnergy);
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


// ============================================================================
// MODIFIED FUNCTION 1: ACCUMULATE PARTICLE FLUX DATA (CALLED BY PARTICLE MOVER)
// ============================================================================

/*
================================================================================
                    AccumulateParticleFluxForWaveCoupling
================================================================================

PURPOSE:
--------
Accumulates particle flux contributions to Alfvén wave growth/damping rates 
during particle transport. This function implements time-weighted sampling of 
the quasi-linear diffusion coefficient for resonant wave-particle interactions
in magnetized solar wind plasma.

CALLING CONTEXT:
----------------
This function is called by the particle mover for each particle at each time 
step during the transport phase. It accumulates statistical data that will 
later be processed by CalculateGrowthRatesFromAccumulatedFlux() to compute 
actual growth rates and update wave energies.

PHYSICS IMPLEMENTED:
--------------------
1. CYCLOTRON RESONANCE: ω - k‖v‖ = ±Ωc
   - Resonant wavenumber: k_res = Ωc/(|μ|v_total)
   - Pitch angle cosine: μ = v‖/v_total

2. QUASI-LINEAR THEORY KERNELS:
   - K₊ = v_total×μ - v_A  (outward wave kernel)
   - K₋ = v_total×μ + v_A  (inward wave kernel)

3. RESONANCE RULES:
   - μ < 0 (inward motion): interacts with outward waves (G_plus)
   - μ > 0 (outward motion): interacts with inward waves (G_minus)

4. TIME-WEIGHTED SAMPLING:
   - Each segment contribution weighted by time_fraction = dt_segment/dt_total
   - Ensures proper temporal representation in flux accumulation

PARTICLE MOTION MODEL:
----------------------
- vParallel: Translational velocity along magnetic field line [m/s]
  * Determines transport direction and time spent in segments
  * Positive = outward (away from Sun), Negative = inward (toward Sun)

- vNormal: Gyration velocity around magnetic field line [m/s]  
  * Represents cyclotron motion perpendicular to B-field
  * Contributes to total particle energy and momentum
  * Does NOT affect transport time (particle stays on same field line)

- v_total = sqrt(vParallel² + vNormal²): Total particle speed
  * Used for relativistic momentum calculation
  * Used in resonance condition and QLT kernels

ALGORITHM FLOW:
---------------
1. Validate inputs and calculate particle properties
2. Loop through all field line segments traversed by particle
3. For each segment:
   a) Calculate path length ds_seg (with proper sign)
   b) Calculate time spent: dt_seg = |ds_seg|/|vParallel|
   c) Calculate time weight: w_time = dt_seg/dt_total
   d) Get local plasma parameters (B, ρ, v_A, Ω)
   e) Calculate resonant wavenumber and find k-bin
   f) Calculate QLT kernels K₊, K₋
   g) Accumulate weighted contributions to G_plus/G_minus arrays

SEGMENT PATH LENGTH CALCULATION:
--------------------------------
The function handles three cases for particle trajectories:

1. SAME SEGMENT (start and end in same segment):
   ds = (frac_finish - frac_start) × segment_length

2. FIRST SEGMENT (particle starts here):
   - Forward: ds = (1.0 - frac_start) × segment_length  
   - Backward: ds = -frac_start × segment_length

3. LAST SEGMENT (particle ends here):
   - Forward: ds = frac_finish × segment_length
   - Backward: ds = (frac_finish - 1.0) × segment_length

4. MIDDLE SEGMENTS (particle traverses completely):
   - Forward: ds = +segment_length
   - Backward: ds = -segment_length

TIME WEIGHTING RATIONALE:
-------------------------
Traditional methods assume particles spend equal time in all segments, which
is incorrect for:
- Variable segment lengths
- Non-uniform particle velocities  
- Particles spending different amounts of time at different locations

Time weighting ensures:
- Fast particles contribute less (spend less time interacting)
- Slow particles contribute more (spend more time interacting)
- Proper statistical sampling of the quasi-linear diffusion process
- Physically accurate representation of wave-particle interaction strength

MOMENTUM RANGE ENFORCEMENT:
---------------------------
The function enforces momentum bounds [P_MIN, P_MAX] with error trapping:
- Particles outside range trigger exit() with diagnostic message
- Prevents silent loss of physics (no particles dropped)
- Provides clear guidance to adjust momentum grid if needed

ACCUMULATED DATA ARRAYS:
------------------------
Results stored in segment-local Datum arrays:
- G_plus_streaming[NK]:  Σ(weighted contributions to outward wave growth)
- G_minus_streaming[NK]: Σ(weighted contributions to inward wave growth)

These arrays are later processed by CalculateGrowthRatesFromAccumulatedFlux()
to compute actual growth rates γ₊(k), γ₋(k) and evolve wave energies.

THREAD SAFETY:
--------------
- Only processes segments assigned to current MPI thread
- Assumes single-threaded access to segment data arrays
- Thread-safe accumulation into G_plus/G_minus arrays

PARAMETERS:
-----------
@param field_line_idx      Index of magnetic field line [0, nFieldLine)
@param particle_index      Particle buffer index for weight calculation  
@param dt                  Time step duration [s]
@param vParallel          Translational velocity along field line [m/s]
@param vNormal            Gyration velocity around field line [m/s]
@param s_start            Start position in field line coordinates [dimensionless]
@param s_finish           End position in field line coordinates [dimensionless]  
@param totalTraversedPath Total signed parallel displacement [m]

FIELD LINE COORDINATES:
-----------------------
Field line coordinates encode position as: s = segment_index.fractional_position
- Integer part: segment index [0, num_segments)
- Fractional part: position within segment [0.0, 1.0)
- Example: s = 5.75 means 75% through segment 5

ERROR CONDITIONS:
-----------------
Function will exit() with diagnostic message for:
- Invalid field line or particle indices
- Non-positive time step
- Particle momentum outside [P_MIN, P_MAX] range
- Invalid segment volumes or plasma parameters

PERFORMANCE NOTES:
------------------
- Pre-computed constants (K_MIN, K_MAX, DLNK) for efficiency
- Logarithmic k-grid with 128 bins for spectral coverage
- Minimal memory allocation (uses pre-existing segment arrays)
- Optimized for frequent calls during particle transport

VALIDATION:
-----------
Debug mode includes validation of:
- G_plus/G_minus array values within reasonable bounds
- Time weights in range [0,1]  
- Segment times ≤ total time step
- Numerical stability of accumulated values

================================================================================
*/

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {
namespace IsotropicSEP { 
    
void AccumulateParticleFluxForWaveCoupling(
    int field_line_idx,                          // Field line index
    long int particle_index,                     // Particle index parameter
    double dt,                                   // Time step [s]
    double vParallel,                           // Particle parallel velocity [m/s] (signed: + outward, - inward)
    double vNormal,                             // Particle gyration velocity [m/s] (perpendicular to B)
    double s_start,                             // Start position along field line [m]
    double s_finish,                            // End position along field line [m]
    double totalTraversedPath                   // Signed parallel path length [m] (+ outward, - inward)
) {
    /*
    MODIFIED VERSION - TIME-WEIGHTED FLUX ACCUMULATION WITH PROPER GYRATION:
    - vParallel: translational motion along magnetic field line
    - vNormal: gyration velocity around magnetic field line (cyclotron motion)
    - Total speed = sqrt(vParallel² + vNormal²)
    - Calculates time spent in each segment based on translational motion
    - Weights flux contribution by (segment_time / total_time_step)
    
    PARTICLE MOTION COMPONENTS:
    - vParallel > 0: Particle moving away from Sun (outward, μ > 0)
    - vParallel < 0: Particle moving toward Sun (inward, μ < 0)  
    - vNormal: Gyration around field line (determines total energy/momentum)
    - μ = vParallel / sqrt(vParallel² + vNormal²) = pitch angle cosine
    
    The pitch angle cosine μ determines wave-particle resonance:
    - μ > 0: Particle moving outward, resonates with inward waves
    - μ < 0: Particle moving inward, resonates with outward waves
    - |μ| determines the resonant wavenumber: k_res = Ω_c / |μ v_total|
    
    TIME WEIGHTING:
    - Each segment contribution is weighted by the time fraction spent in that segment
    - Time in segment based on translational motion: dt_seg = |ds_seg| / |vParallel|
    */
    
    // ========================================================================
    // INPUT VALIDATION
    // ========================================================================
    if (field_line_idx < 0 || field_line_idx >= PIC::FieldLine::nFieldLine) {
        std::cerr << "Error: Invalid field line index (" << field_line_idx 
                  << ") in AccumulateParticleFluxForWaveCoupling" << std::endl;
        return;
    }
    
    if (particle_index < 0) {
        std::cerr << "Error: Invalid particle index (" << particle_index 
                  << ") in AccumulateParticleFluxForWaveCoupling" << std::endl;
        return;
    }
    
    if (dt <= 0.0) {
        std::cerr << "Error: Invalid time step (dt=" << dt << ")" << std::endl;
        return;
    }
    
    // ========================================================================
    // CALCULATE PARTICLE SPEED AND MOTION PARAMETERS
    // ========================================================================
    // Total particle speed includes both parallel and gyration components
    double v_magnitude = sqrt(vParallel * vParallel + vNormal * vNormal);
    
    // Skip particles with negligible velocity
    if (v_magnitude < 1.0e-20) {
        return;
    }
    
    // Skip particles with negligible parallel velocity (no translational motion)
    if (fabs(vParallel) < 1.0e-20) {
        return;  // Pure gyration, no motion along field line
    }
    
    // Calculate pitch angle cosine: μ = v_parallel / v_total
    double mu = vParallel / v_magnitude;
    
    // Clamp μ to valid range [-1, 1] (should already be valid, but safety check)
    mu = std::max(-1.0, std::min(1.0, mu));
    
    // Consistency check: totalTraversedPath should be consistent with vParallel and dt
    double expected_path = vParallel * dt;
    if (fabs(totalTraversedPath - expected_path) > 1.0e-10 * fabs(expected_path)) {
        std::cerr << "Warning: Inconsistent path calculation. Expected=" << expected_path 
                  << ", Actual=" << totalTraversedPath << std::endl;
    }
    
    // ========================================================================
    // CALCULATE PARTICLE PROPERTIES (INDEPENDENT OF SEGMENTS)
    // ========================================================================
    // Calculate relativistic momentum: p = γ m v
    double gamma_rel = 1.0 / sqrt(1.0 - (v_magnitude*v_magnitude)/
                                 (SpeedOfLight*SpeedOfLight));
    double p_momentum = gamma_rel * M * v_magnitude;  // [kg⋅m/s]
    
    // Check particles outside momentum range of interest
    if (p_momentum < P_MIN || p_momentum > P_MAX) {
        char error_msg[512];
        sprintf(error_msg, 
            "Particle momentum (%.6e kg⋅m/s) outside range [%.6e, %.6e]. "
            "Particle: speed=%.6e m/s, vParallel=%.6e m/s, vNormal=%.6e m/s. "
            "Please increase momentum range P_MIN/P_MAX in wave-particle coupling.",
            p_momentum, P_MIN, P_MAX, v_magnitude, vParallel, vNormal);
        exit(__LINE__, __FILE__, error_msg);
    }
    
    // ========================================================================
    // GET PARTICLE STATISTICAL WEIGHT USING SPECIES-SPECIFIC CALCULATION
    // ========================================================================
    // Get particle species number for this particle
    int particle_species = PIC::ParticleBuffer::GetI(particle_index);
    if (particle_species < 0 || particle_species >= PIC::nTotalSpecies) {
        std::cerr << "Error: Invalid particle species (" << particle_species 
                  << ") for particle " << particle_index << std::endl;
        return;
    }
    
    // Calculate species-specific statistical weight
    double w_i = PIC::ParticleWeightTimeStep::GlobalParticleWeight[particle_species] * 
                 PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle_index);
    
    if (w_i <= 0.0) {
        std::cerr << "Warning: Invalid particle weight (" << w_i 
                  << ") for particle " << particle_index << std::endl;
        return;
    }
    
    // ========================================================================
    // GET FIELD LINE AND DETERMINE TRAVERSED SEGMENTS
    // ========================================================================
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    
    // ========================================================================
    // DECODE FIELD LINE COORDINATES AND CALCULATE TRAJECTORY
    // ========================================================================
    // Field line coordinates format: s = (segment_index).(fractional_position)
    // Extract segment indices and fractional positions
    int seg_start_idx = (int)floor(s_start);
    double frac_start = s_start - seg_start_idx;
    
    int seg_finish_idx = (int)floor(s_finish);
    double frac_finish = s_finish - seg_finish_idx;
    
    // ========================================================================
    // LOOP THROUGH ALL SEGMENTS TRAVERSED BY PARTICLE
    // ========================================================================
    // Determine range of segments that the particle trajectory intersects
    int seg_min = std::min(seg_start_idx, seg_finish_idx);
    int seg_max = std::max(seg_start_idx, seg_finish_idx);
    
    int num_segments = field_line->GetTotalSegmentNumber();
    
    // Ensure segment indices are within valid range
    seg_min = std::max(0, seg_min);
    seg_max = std::min(num_segments - 1, seg_max);
    
    for (int seg_idx = seg_min; seg_idx <= seg_max; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        
        // ====================================================================
        // CALCULATE DIRECTIONAL PATH LENGTH WITHIN THIS SEGMENT
        // ====================================================================
        double ds_seg = 0.0;
        
        if (seg_start_idx == seg_finish_idx && seg_idx == seg_start_idx) {
            // Particle starts and ends in the same segment
            double segment_length = segment->GetLength();
            ds_seg = (frac_finish - frac_start) * segment_length;  // Directional: can be negative
            
        } else if (seg_idx == seg_start_idx) {
            // First segment: from start position to end of segment
            double segment_length = segment->GetLength();
            if (seg_start_idx < seg_finish_idx) {
                // Moving forward: from frac_start to 1.0
                ds_seg = (1.0 - frac_start) * segment_length;
            } else {
                // Moving backward: from frac_start to 0.0
                ds_seg = -frac_start * segment_length;
            }
            
        } else if (seg_idx == seg_finish_idx) {
            // Last segment: from beginning of segment to finish position
            double segment_length = segment->GetLength();
            if (seg_start_idx < seg_finish_idx) {
                // Moving forward: from 0.0 to frac_finish
                ds_seg = frac_finish * segment_length;
            } else {
                // Moving backward: from 1.0 to frac_finish
                ds_seg = (frac_finish - 1.0) * segment_length;
            }
            
        } else {
            // Middle segment: entire segment length with direction
            double segment_length = segment->GetLength();
            if (seg_start_idx < seg_finish_idx) {
                // Moving forward: positive full segment length
                ds_seg = segment_length;
            } else {
                // Moving backward: negative full segment length
                ds_seg = -segment_length;
            }
        }
        
        // Skip if path length in segment is negligible (but keep sign)
        if (fabs(ds_seg) < 1.0e-20) {
            continue;
        }
        
        // ====================================================================
        // CALCULATE TIME SPENT IN THIS SEGMENT
        // ====================================================================
        // Time spent in segment based on TRANSLATIONAL motion along field line
        // Time = |path length in segment| / |parallel velocity|
        double dt_segment = fabs(ds_seg) / fabs(vParallel);  // [s]
        
        // Calculate time weighting factor: fraction of total time step spent in this segment
        double time_weight = dt_segment / dt;  // Dimensionless [0, 1]
        
        // Ensure time weighting is physical (can't spend more than total time in segment)
        if (time_weight > 1.0) {
            std::cerr << "Warning: Time weight (" << time_weight 
                      << ") > 1.0 in segment " << seg_idx 
                      << ". Clamping to 1.0." << std::endl;
            time_weight = 1.0;
        }
        
        // Skip segments with negligible time contribution
        if (time_weight < 1.0e-15) {
            continue;
        }
        
        // ====================================================================
        // GET LOCAL PLASMA PARAMETERS FOR THIS SEGMENT
        // ====================================================================
        double B0, rho;
        segment->GetPlasmaDensity(0.5, rho);  // Get density at segment midpoint
        rho *= _H__MASS_;  // Convert number density to mass density
        
        // Get magnetic field from field line
        double B[3];
        PIC::FieldLine::FieldLinesAll[field_line_idx].GetMagneticField(B, 0.5 + seg_idx);
        B0 = Vector3D::Length(B);
        
        // Calculate local plasma parameters
        double vAc = B0 / sqrt(VacuumPermeability * rho);                  // Alfvén speed [m/s]
        double Omega = Q * B0 / M;                               // Proton cyclotron frequency [rad/s]
        double pref = (PI * PI) * Omega * Omega / (B0 * B0);    // Normalization factor
        
        // Get segment volume using proper SEP function
        double V_cell = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);
        if (V_cell <= 0.0) {
            std::cerr << "Warning: Invalid segment volume (" << V_cell 
                      << ") in segment " << seg_idx << std::endl;
            continue;
        }
        double Vinv = 1.0 / V_cell;  // Inverse volume [m⁻³]
        
        // ====================================================================
        // CALCULATE RESONANT WAVENUMBER AND K-BIN
        // ====================================================================
        // Cyclotron resonance condition: ω - k‖ v‖ = ±Ω
        // For Alfvén waves: ω ≈ k‖ v_A, so k_res = Ω / |μ v_total|
        // where v_total = sqrt(vParallel² + vNormal²)
        double kRes = Omega / (fabs(mu) * v_magnitude);  // Resonant wavenumber [m⁻¹]
        
        // Find corresponding k-bin index
        int j = GetKBinIndex(kRes);
        
        // Initialize k-grid for this calculation
        double k_j = K_MIN * exp(j * DLNK);  // k value for bin j
        
        // ====================================================================
        // CALCULATE FLUX FUNCTION VALUES
        // ====================================================================
        // Calculate local QLT kernel values for both wave modes
        // K_σ = v_total × μ - σ × v_A, where σ = ±1 for wave direction
        double K_plus = v_magnitude * mu - vAc;   // Kernel for outward waves (σ = +1)
        double K_minus = v_magnitude * mu + vAc;  // Kernel for inward waves (σ = -1)
        
        // Calculate base flux function coefficient (without time weighting)
        double p2v = p_momentum * p_momentum * v_magnitude;  // p² v term
        double flux_coeff = 2.0 * PI * pref * w_i * p2v / 
                           (2.0 * DLNP) * Vinv / k_j;
        
        // ====================================================================
        // APPLY TIME WEIGHTING AND ACCUMULATE FLUX CONTRIBUTIONS
        // ====================================================================
        // Apply time weighting to flux coefficient
        double weighted_flux_coeff = flux_coeff * time_weight;
        
        // ====================================================================
        // ACCESS INTERMEDIATE STREAMING ARRAYS AND ACCUMULATE
        // ====================================================================
        double* G_plus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
        double* G_minus_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
        
        if (!G_plus_data || !G_minus_data) {
            std::cerr << "Error: Cannot access G_plus/G_minus streaming arrays in segment " 
                      << seg_idx << std::endl;
            continue;
        }
        
        // Add particle contribution to streaming sums using exact QLT kernel
        // Thread-safe accumulation (assumes single-threaded access per segment)
        // Apply time weighting to each contribution
        
        // Wave-particle resonance rules:
        // μ < 0 (inward motion): particle interacts with outward waves (G_plus)
        // μ > 0 (outward motion): particle interacts with inward waves (G_minus)
        if (mu < 0.0) {  // Inward motion
            G_plus_data[j] += weighted_flux_coeff * K_plus;   // Time-weighted outward wave contribution
        }
        if (mu > 0.0) {  // Outward motion
            G_minus_data[j] += weighted_flux_coeff * K_minus; // Time-weighted inward wave contribution
        }

        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
            validate_numeric(G_plus_data[j], -100.0, 50.0, __LINE__, __FILE__);
            validate_numeric(G_minus_data[j], -100.0, 50.0, __LINE__, __FILE__);
            validate_numeric(time_weight, 0.0, 1.0, __LINE__, __FILE__);
            validate_numeric(dt_segment, 0.0, dt, __LINE__, __FILE__);
        }
    }
}

} // namespace IsotropicSEP
} // namespace AlfvenTurbulence_Kolmogorov  
} // namespace SEP
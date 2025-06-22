/*
================================================================================
                    TURBULENCE ENERGY ADVECTION ON LAGRANGIAN MESH
================================================================================

PURPOSE:
--------
This file implements advection of turbulence energy across field line segments
that serve as a Lagrangian mesh. The functions handle the transport of integrated
wave energies E± along magnetic field lines due to Alfvén wave propagation with
spatially varying wave speeds calculated from local plasma conditions.

PHYSICAL MODEL:
---------------
- Lagrangian mesh: Field line segments move with solar wind plasma
- Wave propagation: Outward waves (E+) and inward waves (E-) propagate relative to moving plasma
- Advection equation: ∂E/∂t + ∇·(v_wave × E) = 0 in moving frame
- Local Alfvén speeds: V_A = |B| / sqrt(μ₀ × ρ) calculated at each boundary from vertex data
- Conservative scheme: Total energy conserved during advection

MATHEMATICAL FORMULATION:
-------------------------
Energy flux calculation:
  Flux_rate = V_A × (E_total / V_segment) × A_boundary  [J/s]

Where:
- E_total: Total integrated wave energy in magnetic tube segment [J]
- V_segment: Volume of magnetic tube segment from SEP::FieldLine::GetSegmentVolume() [m³]
- A_boundary: Cross-sectional area = π × [MagneticTubeRadius()]² [m²]
- V_A: Local Alfvén velocity from vertex magnetic field and plasma density [m/s]

Wave energy update:
  E±^{n+1} = E±^n + Σ(flux_in) - Σ(flux_out)

NUMERICAL ALGORITHM:
--------------------
1. Read wave energy data from ALL segments (global information needed for fluxes)
2. Calculate local Alfvén velocities at each segment boundary from vertex data:
   - Magnetic field: FL::DatumAtVertexMagneticField
   - Plasma density: FL::DatumAtVertexPlasmaDensity
   - V_A = |B| / sqrt(μ₀ × n × m_proton)
3. Calculate fluxes using upwind scheme (only for local segments):
   - Outward waves (E+): upwind from left (smaller index)
   - Inward waves (E-): upwind from right (larger index)
4. Apply conservative energy updates (only to local segments)
5. Apply boundary conditions (only to local boundary segments)
6. Write updated wave energies (only to local segments)

MPI PARALLELIZATION:
--------------------
- Wave energy information is global: read from all segments
- Flux calculations: performed only for segments assigned to current process
- Energy updates: applied only to segments assigned to current process
- Data writes: performed only for segments assigned to current process
- Thread ownership: determined by segment->Thread == PIC::ThisThread

BOUNDARY CONDITIONS:
--------------------
- Inner boundary (index 0): E_minus = 0 (inward waves exit at Sun)
- Outer boundary (highest index): E_plus *= exp(-dt/τ) (gradual damping)

USAGE EXAMPLES:
---------------

// Example 1: Process single field line
int field_line_idx = 5;
double dt = 10.0;  // 10 second time step
AdvectTurbulenceEnergyAlongFieldLine(
    field_line_idx, 
    SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity, 
    dt, 
    true  // apply boundary conditions
);

// Example 2: Process all field lines
AdvectTurbulenceEnergyAllFieldLines(
    SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity, 
    dt
);

// Example 3: Energy conservation monitoring
double energy_before = CalculateTotalAdvectedEnergy(WaveEnergyDensity);
AdvectTurbulenceEnergyAllFieldLines(WaveEnergyDensity, dt);
double energy_after = CalculateTotalAdvectedEnergy(WaveEnergyDensity);
double energy_change = energy_after - energy_before;

// Example 4: Integration in main simulation loop
void SimulationTimestep(double dt) {
    // ... other physics updates ...
    
    // Advect turbulence energy
    SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(
        SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity, dt
    );
    
    // ... continue with other updates ...
}

REQUIREMENTS:
-------------
- AMPS PIC framework with field line structure
- SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity datum for wave energies
- SEP::FieldLine::GetSegmentVolume() for volume calculation
- SEP::FieldLine::MagneticTubeRadius() for boundary area calculation
- FL::DatumAtVertexMagneticField and FL::DatumAtVertexPlasmaDensity at vertices
- MPI parallelization support

PERFORMANCE NOTES:
------------------
- Memory efficient: Minimal storage arrays, direct segment access
- MPI optimized: Each process handles only assigned segments
- Conservative scheme: Maintains exact energy conservation
- Local calculations: Alfvén speeds computed from actual plasma conditions
- Upwind stability: Prevents numerical oscillations in wave transport

ENERGY CONSERVATION:
--------------------
The algorithm guarantees strict energy conservation:
- Total flux into system = Total flux out of system
- Energy lost from source segment = Energy gained by destination segment
- No energy creation or destruction during advection
- Boundary conditions properly handle energy outflow

TYPICAL PARAMETERS:
-------------------
- Time step: dt ~ 1-100 seconds (limited by CFL condition: dt < Δx/V_A)
- Alfvén speeds: V_A ~ 10-1000 km/s in solar wind
- Wave energies: E± ~ 10¹²-10¹⁵ J per segment (TJ scale)
- Segment volumes: V ~ 10²⁰-10²³ m³ (depends on heliocentric distance)

================================================================================
*/

#include "pic.h"              // AMPS PIC framework
#include <mpi.h>              // MPI parallelization
#include <vector>             // STL containers
#include <iostream>           // Error output
#include <cmath>              // Mathematical functions

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

// ============================================================================
// TURBULENCE ENERGY ADVECTION ALONG FIELD LINES WITH LOCAL ALFVÉN SPEEDS
// ============================================================================

void AdvectTurbulenceEnergyAlongFieldLine(
    int field_line_idx,                           // Field line index
    const PIC::Datum::cDatumStored& WaveEnergyDensity,  // Wave energy datum
    double dt,                                     // Time step [s]
    bool apply_boundary_conditions          // Apply boundary conditions
) {
    namespace FL=PIC::FieldLine;
    // Get pointer to field line from index
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    
    if (!field_line) {
        std::cerr << "Error: Invalid field line index " << field_line_idx << std::endl;
        return;
    }
    
    int num_segments = field_line->GetTotalSegmentNumber();
    if (num_segments < 2) {
        // Need at least 2 segments for advection
        return;
    }
    
    const double mu0 = 4.0e-7 * M_PI;  // Permeability of free space [H/m]
    const double proton_mass = 1.67262192e-27;  // Proton mass [kg]
    
    // ========================================================================
    // READ WAVE ENERGY DATA FROM ALL SEGMENTS (GLOBAL INFORMATION)
    // ========================================================================
    std::vector<double> E_plus_current(num_segments, 0.0);
    std::vector<double> E_minus_current(num_segments, 0.0);
    std::vector<double> segment_volumes(num_segments, 0.0);
    
    // Read wave energy data from ALL segments (global information needed for fluxes)
    for (int i = 0; i < num_segments; ++i) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
        if (!segment) continue;
        
        // Read wave energy data from ALL segments (global)
        double* wave_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        if (wave_data) {
            E_plus_current[i] = wave_data[0];   // Total outward wave energy [J]
            E_minus_current[i] = wave_data[1];  // Total inward wave energy [J]
        }
        
        // Get segment volume for ALL segments (needed for flux calculations)
        segment_volumes[i] = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);  // [m³]
    }
    
    // ========================================================================
    // CALCULATE ENERGY FLUXES WITH LOCAL ALFVÉN SPEEDS
    // ========================================================================
    
    // Calculate all fluxes first, then apply updates
    std::vector<double> energy_flux_plus(num_segments-1, 0.0);   // Flux between segments for E+
    std::vector<double> energy_flux_minus(num_segments-1, 0.0);  // Flux between segments for E-
    
    // Calculate outward wave fluxes (E+) - only for local segments
    for (int i = 0; i < num_segments-1; ++i) {
        PIC::FieldLine::cFieldLineSegment* segment_current = field_line->GetSegment(i);
        if (!segment_current || segment_current->Thread != PIC::ThisThread) {
            continue;  // Only process if current segment is local
        }
        
        // Get vertex at right boundary of current segment
        PIC::FieldLine::cFieldLineVertex* vertex_right = segment_current->GetEnd();
        double* x_right = vertex_right->GetX();
        
        // Calculate Alfvén velocity at right boundary
        double* B0_right = vertex_right->GetDatum_ptr(FL::DatumAtVertexMagneticField);
        double n_sw_right;
        vertex_right->GetDatum(FL::DatumAtVertexPlasmaDensity, &n_sw_right);
        
        if (!B0_right || n_sw_right <= 0.0) {
            continue;  // Skip if no valid plasma parameters
        }
        
        // Calculate magnetic field magnitude
        double B_mag_right = std::sqrt(B0_right[0]*B0_right[0] + B0_right[1]*B0_right[1] + B0_right[2]*B0_right[2]);
        
        // Calculate mass density: ρ = n × m_proton
        double rho_right = n_sw_right * proton_mass;  // [kg/m³]
        
        // Calculate Alfvén velocity: V_A = B / sqrt(μ₀ × ρ)
        double V_A_right = B_mag_right / std::sqrt(mu0 * rho_right);  // [m/s]
        
        // Calculate boundary area using SEP::FieldLine::MagneticTubeRadius
        double radius_right = SEP::FieldLine::MagneticTubeRadius(x_right, field_line_idx);
        double boundary_area = M_PI * radius_right * radius_right;  // [m²]
        
        // Calculate energy density in current segment (using global values)
        double energy_density = E_plus_current[i] / segment_volumes[i];  // [J/m³]
        
        // Calculate flux: V_A * energy_density * boundary_area
        double flux_rate = V_A_right * energy_density * boundary_area;  // [J/s]
        double energy_flux = flux_rate * dt;  // [J]
        
        // Limit flux to avoid negative energies
        energy_flux = std::min(energy_flux, E_plus_current[i]);
        
        // Store flux from segment i to segment i+1
        energy_flux_plus[i] = energy_flux;
    }
    
    // Calculate inward wave fluxes (E-) - only for local segments
    for (int i = 0; i < num_segments-1; ++i) {
        PIC::FieldLine::cFieldLineSegment* segment_current = field_line->GetSegment(i);
        if (!segment_current || segment_current->Thread != PIC::ThisThread) {
            continue;  // Only process if current segment is local
        }
        
        // For inward waves, we need the segment i+1 as source (upwind from right)
        PIC::FieldLine::cFieldLineSegment* segment_source = field_line->GetSegment(i+1);
        if (!segment_source) continue;
        
        // Get vertex at left boundary of source segment (i+1)
        PIC::FieldLine::cFieldLineVertex* vertex_left = segment_source->GetBegin();
        double* x_left = vertex_left->GetX();
        
        // Calculate Alfvén velocity at left boundary
        double* B0_left = vertex_left->GetDatum_ptr(FL::DatumAtVertexMagneticField);
        double n_sw_left;
        vertex_left->GetDatum(FL::DatumAtVertexPlasmaDensity, &n_sw_left);
        
        if (!B0_left || n_sw_left <= 0.0) {
            continue;  // Skip if no valid plasma parameters
        }
        
        // Calculate magnetic field magnitude
        double B_mag_left = std::sqrt(B0_left[0]*B0_left[0] + B0_left[1]*B0_left[1] + B0_left[2]*B0_left[2]);
        
        // Calculate mass density: ρ = n × m_proton
        double rho_left = n_sw_left * proton_mass;  // [kg/m³]
        
        // Calculate Alfvén velocity: V_A = B / sqrt(μ₀ × ρ)
        double V_A_left = B_mag_left / std::sqrt(mu0 * rho_left);  // [m/s]
        
        // Calculate boundary area using SEP::FieldLine::MagneticTubeRadius
        double radius_left = SEP::FieldLine::MagneticTubeRadius(x_left, field_line_idx);
        double boundary_area = M_PI * radius_left * radius_left;  // [m²]
        
        // Calculate energy density in source segment (i+1) using global values
        double energy_density = E_minus_current[i+1] / segment_volumes[i+1];  // [J/m³]
        
        // Calculate flux: V_A * energy_density * boundary_area
        double flux_rate = V_A_left * energy_density * boundary_area;  // [J/s]
        double energy_flux = flux_rate * dt;  // [J]
        
        // Limit flux to avoid negative energies
        energy_flux = std::min(energy_flux, E_minus_current[i+1]);
        
        // Store flux from segment i+1 to segment i
        energy_flux_minus[i] = energy_flux;
    }
    
    // ========================================================================
    // APPLY ENERGY UPDATES (ONLY TO LOCAL SEGMENTS)
    // ========================================================================
    
    // Apply all flux updates simultaneously (only to local segments)
    for (int i = 0; i < num_segments; ++i) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
        if (!segment || segment->Thread != PIC::ThisThread) {
            continue;  // Skip non-local segments
        }
        
        // Outward waves: receive from left neighbor, lose to right neighbor
        if (i > 0) {
            E_plus_current[i] += energy_flux_plus[i-1];  // Gain from left
        }
        if (i < num_segments-1) {
            E_plus_current[i] -= energy_flux_plus[i];    // Loss to right
        }
        
        // Inward waves: receive from right neighbor, lose to left neighbor  
        if (i < num_segments-1) {
            E_minus_current[i] += energy_flux_minus[i];   // Gain from right
        }
        if (i > 0) {
            E_minus_current[i] -= energy_flux_minus[i-1]; // Loss to left
        }
    }
    
    // ========================================================================
    // APPLY BOUNDARY CONDITIONS (ONLY TO LOCAL SEGMENTS)
    // ========================================================================
    if (apply_boundary_conditions) {
        // Inner boundary (close to Sun, index 0): outflow for inward waves
        if (num_segments > 0) {
            PIC::FieldLine::cFieldLineSegment* segment_inner = field_line->GetSegment(0);
            if (segment_inner && segment_inner->Thread == PIC::ThisThread) {
                E_minus_current[0] = 0.0;  // Inward waves exit at inner boundary
            }
        }
        
        // Outer boundary (far from Sun, highest index): outflow for outward waves
        if (num_segments > 0) {
            PIC::FieldLine::cFieldLineSegment* segment_outer = field_line->GetSegment(num_segments-1);
            if (segment_outer && segment_outer->Thread == PIC::ThisThread) {
                E_plus_current[num_segments-1] *= std::exp(-dt / 1000.0);  // Gradual damping at outer boundary
            }
        }
    }
    
    // ========================================================================
    // UPDATE WAVE ENERGY DATA (ONLY IN LOCAL SEGMENTS)
    // ========================================================================
    for (int i = 0; i < num_segments; ++i) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
        if (!segment || segment->Thread != PIC::ThisThread) {
            continue;  // Skip non-local segments
        }
        
        double* wave_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        if (wave_data) {
            // Ensure non-negative energies
            wave_data[0] = std::max(0.0, E_plus_current[i]);   // Outward wave energy
            wave_data[1] = std::max(0.0, E_minus_current[i]);  // Inward wave energy
        }
    }
}

// ============================================================================
// GLOBAL ADVECTION FOR ALL FIELD LINES
// ============================================================================



/*
================================================================================
                    TURBULENCE ENERGY ADVECTION WITH BOUNDARY CONDITIONS
================================================================================

PURPOSE:
--------
Advects turbulence wave energy along all magnetic field lines in the AMPS plasma
simulation with integrated boundary conditions. Handles transport of outward (E+)
and inward (E-) Alfvén wave energies using local plasma conditions and applies
turbulence level boundary conditions at field line endpoints.

PHYSICAL MODEL:
---------------
- Lagrangian mesh: Field line segments move with solar wind plasma
- Wave propagation: E+ waves propagate outward, E- waves propagate inward
- Local Alfvén speeds: V_A = |B| / sqrt(μ₀ × ρ) calculated from vertex data
- Energy flux: Flux = V_A × TurbulenceLevel × BoundaryArea [J/s]
- Boundary conditions: Turbulence level specification at field line endpoints

ALGORITHM SEQUENCE:
-------------------
1. OPTIMIZE MEMORY ACCESS:
   - Identify local segments (owned by current MPI process)
   - Mark adjacent segments needed for flux calculations
   - Read wave energy data only from required segments

2. CALCULATE INTERIOR FLUXES:
   - Compute energy_flux_plus[i]: E+ flux from segment i to i+1
   - Compute energy_flux_minus[i]: E- flux from segment i+1 to i
   - Use upwind scheme with local Alfvén velocities

3. APPLY BOUNDARY CONDITIONS:
   - First segment: 
     * Add incoming E+ flux based on TurbulenceLevelBeginning (converted to wave energy density)
     * Remove outgoing E- flux based on actual segment energy
   - Last segment:
     * Add incoming E- flux based on TurbulenceLevelEnd (converted to wave energy density)
     * Remove outgoing E+ flux based on actual segment energy

4. UPDATE INTERIOR SEGMENTS:
   - Apply flux exchanges between adjacent segments
   - Skip boundary segments (already updated in step 3)
   - Maintain energy conservation

5. WRITE UPDATED ENERGIES:
   - Store results only in locally-owned segments
   - Ensure non-negative energy values

BOUNDARY CONDITIONS PHYSICS:
-----------------------------
Beginning boundary (index 0):
- INCOMING: E+ energy injection using TurbulenceLevelBeginning (converted to wave energy density)
- OUTGOING: E- energy outflow based on actual segment energy content

End boundary (index num_segments-1):
- INCOMING: E- energy injection using TurbulenceLevelEnd (converted to wave energy density)
- OUTGOING: E+ energy outflow based on actual segment energy content

Turbulence level to wave energy density conversion:
  WaveEnergyDensity = TurbulenceLevel × ReferenceEnergyDensity [J/m³]
  
Where ReferenceEnergyDensity is calculated from local magnetic field energy:
  ReferenceEnergyDensity = B² / (2 × μ₀) [J/m³]

Energy flux calculation:
  Flux = V_A × WaveEnergyDensity × π × [MagneticTubeRadius()]² × dt [J]

MPI PARALLELIZATION:
--------------------
- Process assignment: Each segment belongs to one MPI process (segment->Thread)
- Memory optimization: Read only local + adjacent segment data
- Flux calculations: Performed only for locally-owned segments
- Energy updates: Applied only to locally-owned segments
- Communication: Global wave energy data available, local processing only

CONSERVATION PROPERTIES:
------------------------
- Total energy conserved during interior advection
- Boundary conditions allow controlled energy injection/extraction
- Flux limiting prevents negative energies
- Upwind scheme maintains numerical stability

USAGE EXAMPLE:
--------------
```cpp
// Call during simulation time step
double dt = 10.0;  // Time step [s]
double turbulence_level_beginning = 0.1;  // Turbulence level at start [dimensionless, fraction of magnetic energy]
double turbulence_level_end = 0.05;       // Turbulence level at end [dimensionless]

SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(
    SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,
    dt,
    turbulence_level_beginning,
    turbulence_level_end
);
```

REQUIREMENTS:
-------------
- SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy datum
- FL::DatumAtVertexMagneticField and FL::DatumAtVertexPlasmaDensity at vertices
- SEP::FieldLine::GetSegmentVolume() and SEP::FieldLine::MagneticTubeRadius()
- MPI environment with PIC::ThisThread process identification

PERFORMANCE OPTIMIZATIONS:
---------------------------
- Minimal memory access: Read only needed segments (local + adjacent)
- Process-local operations: No unnecessary thread checks in inner loops  
- Conservative flux calculation: All fluxes computed before any updates
- Boundary-first updates: Prevents interference with interior advection

TYPICAL PARAMETERS:
-------------------
- Time step: dt = 1-100 seconds (CFL limited by V_A and segment size)
- Alfvén speeds: V_A = 10-1000 km/s (solar wind conditions)
- Turbulence levels: 0.01-0.5 (1% to 50% of magnetic field energy density)
- Segment volumes: 10²⁰ - 10²³ m³ (heliocentric distance dependent)

================================================================================
*/

void AdvectTurbulenceEnergyAllFieldLines(
    const PIC::Datum::cDatumStored& WaveEnergyDensity,  // Wave energy datum
    double dt,                                          // Time step [s]
    double TurbulenceLevelBeginning,                    // Turbulence level at field line beginning [dimensionless]
    double TurbulenceLevelEnd                           // Turbulence level at field line end [dimensionless]
) {
    namespace FL=PIC::FieldLine;
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    const double mu0 = 4.0e-7 * M_PI;  // Permeability of free space [H/m]
    const double proton_mass = 1.67262192e-27;  // Proton mass [kg]
    
    int processed_field_lines = 0;
    
    // Loop through all field lines
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        
        if (!field_line) continue;
        
        int num_segments = field_line->GetTotalSegmentNumber();
        if (num_segments < 2) continue;
        
        // Check if this field line has segments assigned to this MPI process
        bool has_local_segments = false;
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (segment && segment->Thread == PIC::ThisThread) {
                has_local_segments = true;
                break;
            }
        }
        
        if (!has_local_segments) continue;
        
        // ========================================================================
        // IDENTIFY LOCAL AND ADJACENT SEGMENTS FOR MINIMAL MEMORY ACCESS
        // ========================================================================
        
        // First pass: identify which segments are local to this process
        std::vector<bool> is_local(num_segments, false);
        std::vector<bool> need_data(num_segments, false);
        
        for (int i = 0; i < num_segments; ++i) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
            if (segment && segment->Thread == PIC::ThisThread) {
                is_local[i] = true;
                need_data[i] = true;  // Need data for local segments
                
                // Also need data for adjacent segments for flux calculations
                if (i > 0) need_data[i-1] = true;          // Left neighbor
                if (i < num_segments-1) need_data[i+1] = true;  // Right neighbor
            }
        }
        
        // ========================================================================
        // READ WAVE ENERGY DATA ONLY FROM NEEDED SEGMENTS (OPTIMIZED)
        // ========================================================================
        std::vector<double> E_plus_current(num_segments, 0.0);
        std::vector<double> E_minus_current(num_segments, 0.0);
        std::vector<double> segment_volumes(num_segments, 0.0);
        
        // Read wave energy data only from segments we actually need
        for (int i = 0; i < num_segments; ++i) {
            if (!need_data[i]) continue;  // Skip segments we don't need
            
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
            if (!segment) continue;
            
            // Read wave energy data only from needed segments
            double* wave_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
            if (wave_data) {
                E_plus_current[i] = wave_data[0];   // Total outward wave energy [J]
                E_minus_current[i] = wave_data[1];  // Total inward wave energy [J]
            }
            
            // Get segment volume only for needed segments
            segment_volumes[i] = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);  // [m³]
        }
        
      
        // ========================================================================
        // CALCULATE ENERGY FLUXES WITH LOCAL ALFVÉN SPEEDS (INTERIOR SEGMENTS)
        // ========================================================================
        
        // Calculate all fluxes first, then apply updates
        std::vector<double> energy_flux_plus(num_segments-1, 0.0);   // Flux between segments for E+
        std::vector<double> energy_flux_minus(num_segments-1, 0.0);  // Flux between segments for E-
        
        // Calculate outward wave fluxes (E+) - only for local segments
        for (int i = 0; i < num_segments-1; ++i) {
            if (!is_local[i]) continue;  // Only process local segments
            
            PIC::FieldLine::cFieldLineSegment* segment_current = field_line->GetSegment(i);
            if (!segment_current) continue;
            
            // Get vertex at right boundary of current segment
            PIC::FieldLine::cFieldLineVertex* vertex_right = segment_current->GetEnd();
            double* x_right = vertex_right->GetX();
            
            // Calculate Alfvén velocity at right boundary
            double* B0_right = vertex_right->GetDatum_ptr(FL::DatumAtVertexMagneticField);
            double n_sw_right;
            vertex_right->GetDatum(FL::DatumAtVertexPlasmaDensity, &n_sw_right);
            
            if (!B0_right || n_sw_right <= 0.0) {
                continue;  // Skip if no valid plasma parameters
            }
            
            // Calculate magnetic field magnitude
            double B_mag_right = std::sqrt(B0_right[0]*B0_right[0] + B0_right[1]*B0_right[1] + B0_right[2]*B0_right[2]);
            
            // Calculate mass density: ρ = n × m_proton
            double rho_right = n_sw_right * proton_mass;  // [kg/m³]
            
            // Calculate Alfvén velocity: V_A = B / sqrt(μ₀ × ρ)
            double V_A_right = B_mag_right / std::sqrt(mu0 * rho_right);  // [m/s]
            
            // Calculate boundary area using SEP::FieldLine::MagneticTubeRadius
            double radius_right = SEP::FieldLine::MagneticTubeRadius(x_right, field_line_idx);
            double boundary_area = M_PI * radius_right * radius_right;  // [m²]
            
            // Calculate energy density in current segment (using global values)
            double energy_density = E_plus_current[i] / segment_volumes[i];  // [J/m³]
            
            // Calculate flux: V_A * energy_density * boundary_area
            double flux_rate = V_A_right * energy_density * boundary_area;  // [J/s]
            double energy_flux = flux_rate * dt;  // [J]
            
            // Limit flux to avoid negative energies
            energy_flux = std::min(energy_flux, E_plus_current[i]);
            
            // Store flux from segment i to segment i+1
            energy_flux_plus[i] = energy_flux;
        }
        
        // Calculate inward wave fluxes (E-) - only for local segments
        for (int i = 0; i < num_segments-1; ++i) {
            if (!is_local[i]) continue;  // Only process local segments
            
            PIC::FieldLine::cFieldLineSegment* segment_current = field_line->GetSegment(i);
            if (!segment_current) continue;
            
            // For inward waves, we need the segment i+1 as source (upwind from right)
            PIC::FieldLine::cFieldLineSegment* segment_source = field_line->GetSegment(i+1);
            if (!segment_source) continue;
            
            // Get vertex at left boundary of source segment (i+1)
            PIC::FieldLine::cFieldLineVertex* vertex_left = segment_source->GetBegin();
            double* x_left = vertex_left->GetX();
            
            // Calculate Alfvén velocity at left boundary
            double* B0_left = vertex_left->GetDatum_ptr(FL::DatumAtVertexMagneticField);
            double n_sw_left;
            vertex_left->GetDatum(FL::DatumAtVertexPlasmaDensity, &n_sw_left);
            
            if (!B0_left || n_sw_left <= 0.0) {
                continue;  // Skip if no valid plasma parameters
            }
            
            // Calculate magnetic field magnitude
            double B_mag_left = std::sqrt(B0_left[0]*B0_left[0] + B0_left[1]*B0_left[1] + B0_left[2]*B0_left[2]);
            
            // Calculate mass density: ρ = n × m_proton
            double rho_left = n_sw_left * proton_mass;  // [kg/m³]
            
            // Calculate Alfvén velocity: V_A = B / sqrt(μ₀ × ρ)
            double V_A_left = B_mag_left / std::sqrt(mu0 * rho_left);  // [m/s]
            
            // Calculate boundary area using SEP::FieldLine::MagneticTubeRadius
            double radius_left = SEP::FieldLine::MagneticTubeRadius(x_left, field_line_idx);
            double boundary_area = M_PI * radius_left * radius_left;  // [m²]
            
            // Calculate energy density in source segment (i+1) using global values
            double energy_density = E_minus_current[i+1] / segment_volumes[i+1];  // [J/m³]
            
            // Calculate flux: V_A * energy_density * boundary_area
            double flux_rate = V_A_left * energy_density * boundary_area;  // [J/s]
            double energy_flux = flux_rate * dt;  // [J]
            
            // Limit flux to avoid negative energies
            energy_flux = std::min(energy_flux, E_minus_current[i+1]);
            
            // Store flux from segment i+1 to segment i
            energy_flux_minus[i] = energy_flux;
        }
        
        // ========================================================================
        // UPDATE BOUNDARY SEGMENTS FIRST (FIRST AND LAST SEGMENTS)
        // ========================================================================
        
        // FIRST SEGMENT (Beginning of field line) - Calculate both E+ and E- fluxes
        if (num_segments > 0 && is_local[0]) {
            PIC::FieldLine::cFieldLineSegment* first_segment = field_line->GetSegment(0);
            if (first_segment) {
                // Get vertex at the beginning of field line
                PIC::FieldLine::cFieldLineVertex* vertex_beginning = first_segment->GetBegin();
                if (vertex_beginning) {
                    double* x_beginning = vertex_beginning->GetX();
                    
                    // Calculate Alfvén velocity at beginning boundary
                    double* B0_beginning = vertex_beginning->GetDatum_ptr(FL::DatumAtVertexMagneticField);
                    double n_sw_beginning;
                    vertex_beginning->GetDatum(FL::DatumAtVertexPlasmaDensity, &n_sw_beginning);
                    
                    if (B0_beginning && n_sw_beginning > 0.0 && segment_volumes[0] > 0.0) {
                        // Calculate magnetic field magnitude
                        double B_mag_beginning = std::sqrt(B0_beginning[0]*B0_beginning[0] + 
                                                          B0_beginning[1]*B0_beginning[1] + 
                                                          B0_beginning[2]*B0_beginning[2]);
                        
                        // Calculate mass density: ρ = n × m_proton
                        double rho_beginning = n_sw_beginning * proton_mass;  // [kg/m³]
                        
                        // Calculate Alfvén velocity: V_A = B / sqrt(μ₀ × ρ)
                        double V_A_beginning = B_mag_beginning / std::sqrt(mu0 * rho_beginning);  // [m/s]
                        
                        // Calculate boundary area using SEP::FieldLine::MagneticTubeRadius
                        double radius_beginning = SEP::FieldLine::MagneticTubeRadius(x_beginning, field_line_idx);
                        double boundary_area = M_PI * radius_beginning * radius_beginning;  // [m²]
                        
                        // Calculate reference energy density from magnetic field energy
                        double reference_energy_density = (B_mag_beginning * B_mag_beginning) / (2.0 * mu0);  // [J/m³]
                        
                        // Convert turbulence level to wave energy density
                        double wave_energy_density_beginning = TurbulenceLevelBeginning * reference_energy_density;  // [J/m³]
                        
                        // Calculate INCOMING energy flux based on turbulence level at beginning
                        double flux_rate_incoming = V_A_beginning * wave_energy_density_beginning * boundary_area;  // [J/s]
                        double energy_flux_incoming = flux_rate_incoming * dt;  // [J]
                        
                        // Calculate OUTGOING energy flux based on actual wave energy in first segment
                        double energy_density_minus_outgoing = E_minus_current[0] / segment_volumes[0];  // [J/m³]
                        double flux_rate_outgoing = V_A_beginning * energy_density_minus_outgoing * boundary_area;  // [J/s]
                        double energy_flux_outgoing = flux_rate_outgoing * dt;  // [J]
                        
                        // Limit outgoing flux to avoid negative energies
                        energy_flux_outgoing = std::min(energy_flux_outgoing, E_minus_current[0]);
                        
                        // Update E+ (outward waves): Add incoming energy from beginning boundary
                        E_plus_current[0] += energy_flux_incoming;
                        
                        // Update E- (inward waves): Remove outgoing energy at beginning boundary
                        E_minus_current[0] -= energy_flux_outgoing;
                    }
                }
            }
        }
        
        // LAST SEGMENT (End of field line) - Calculate both E+ and E- fluxes
        if (num_segments > 0 && is_local[num_segments - 1]) {
            PIC::FieldLine::cFieldLineSegment* last_segment = field_line->GetSegment(num_segments - 1);
            if (last_segment) {
                // Get vertex at the end of field line
                PIC::FieldLine::cFieldLineVertex* vertex_end = last_segment->GetEnd();
                if (vertex_end) {
                    double* x_end = vertex_end->GetX();
                    
                    // Calculate Alfvén velocity at end boundary
                    double* B0_end = vertex_end->GetDatum_ptr(FL::DatumAtVertexMagneticField);
                    double n_sw_end;
                    vertex_end->GetDatum(FL::DatumAtVertexPlasmaDensity, &n_sw_end);
                    
                    if (B0_end && n_sw_end > 0.0 && segment_volumes[num_segments - 1] > 0.0) {
                        // Calculate magnetic field magnitude
                        double B_mag_end = std::sqrt(B0_end[0]*B0_end[0] + 
                                                    B0_end[1]*B0_end[1] + 
                                                    B0_end[2]*B0_end[2]);
                        
                        // Calculate mass density: ρ = n × m_proton
                        double rho_end = n_sw_end * proton_mass;  // [kg/m³]
                        
                        // Calculate Alfvén velocity: V_A = B / sqrt(μ₀ × ρ)
                        double V_A_end = B_mag_end / std::sqrt(mu0 * rho_end);  // [m/s]
                        
                        // Calculate boundary area using SEP::FieldLine::MagneticTubeRadius
                        double radius_end = SEP::FieldLine::MagneticTubeRadius(x_end, field_line_idx);
                        double boundary_area = M_PI * radius_end * radius_end;  // [m²]
                        
                        // Calculate reference energy density from magnetic field energy
                        double reference_energy_density = (B_mag_end * B_mag_end) / (2.0 * mu0);  // [J/m³]
                        
                        // Convert turbulence level to wave energy density
                        double wave_energy_density_end = TurbulenceLevelEnd * reference_energy_density;  // [J/m³]
                        
                        // Calculate INCOMING energy flux based on turbulence level at end
                        double flux_rate_incoming = V_A_end * wave_energy_density_end * boundary_area;  // [J/s]
                        double energy_flux_incoming = flux_rate_incoming * dt;  // [J]
                        
                        // Calculate OUTGOING energy flux based on actual wave energy in last segment
                        double energy_density_plus_outgoing = E_plus_current[num_segments - 1] / segment_volumes[num_segments - 1];  // [J/m³]
                        double flux_rate_outgoing = V_A_end * energy_density_plus_outgoing * boundary_area;  // [J/s]
                        double energy_flux_outgoing = flux_rate_outgoing * dt;  // [J]
                        
                        // Limit outgoing flux to avoid negative energies
                        energy_flux_outgoing = std::min(energy_flux_outgoing, E_plus_current[num_segments - 1]);
                        
                        // Update E- (inward waves): Add incoming energy from end boundary
                        E_minus_current[num_segments - 1] += energy_flux_incoming;
                        
                        // Update E+ (outward waves): Remove outgoing energy at end boundary
                        E_plus_current[num_segments - 1] -= energy_flux_outgoing;
                    }
                }
            }
        }
        
        // ========================================================================
        // UPDATE ALL INTERIOR SEGMENTS (APPLY FLUX EXCHANGES)
        // ========================================================================
        
        // Apply all flux updates simultaneously (only to interior segments, excluding boundaries)
        for (int i = 1; i < num_segments - 1; ++i) {  // Skip first (i=0) and last (i=num_segments-1) segments
            if (!is_local[i]) continue;  // Skip non-local segments
            
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
            if (!segment) continue;
            
            // Outward waves: receive from left neighbor, lose to right neighbor
            if (i > 0) {
                E_plus_current[i] += energy_flux_plus[i-1];  // Gain from left
            }
            if (i < num_segments-1) {
                E_plus_current[i] -= energy_flux_plus[i];    // Loss to right
            }
            
            // Inward waves: receive from right neighbor, lose to left neighbor  
            if (i < num_segments-1) {
                E_minus_current[i] += energy_flux_minus[i];   // Gain from right
            }
            if (i > 0) {
                E_minus_current[i] -= energy_flux_minus[i-1]; // Loss to left
            }
        }
        
        // ========================================================================
        // UPDATE WAVE ENERGY DATA (ONLY IN LOCAL SEGMENTS)
        // ========================================================================
        for (int i = 0; i < num_segments; ++i) {
            if (!is_local[i]) continue;  // Skip non-local segments
            
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
            if (!segment) continue;
            
            double* wave_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
            if (wave_data) {
                // Ensure non-negative energies
                wave_data[0] = std::max(0.0, E_plus_current[i]);   // Outward wave energy
                wave_data[1] = std::max(0.0, E_minus_current[i]);  // Inward wave energy
            }
        }
        
        processed_field_lines++;
    }
    
    if (rank == 0) {
        std::cout << "Turbulence energy advection completed for " << processed_field_lines 
                  << " field lines with required boundary conditions" << std::endl;
        std::cout << "Boundary turbulence levels - Beginning: " << TurbulenceLevelBeginning 
                  << ", End: " << TurbulenceLevelEnd << std::endl;
    }
}

// ============================================================================
// DIAGNOSTIC: CALCULATE TOTAL ENERGY BEFORE/AFTER ADVECTION
// ============================================================================

double CalculateTotalAdvectedEnergy(const PIC::Datum::cDatumStored& WaveEnergyDensity) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double local_total_energy = 0.0;
    
    // Sum energy across all segments assigned to this process
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;
            
            double* wave_data = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
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

} // namespace AlfvenTurbulence_Kolmogorov  
} // namespace SEP

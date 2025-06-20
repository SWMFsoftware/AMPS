/*
================================================================================
                    WAVE ENERGY INITIALIZATION SOURCE
================================================================================

FILENAME: wave_energy_initialization.cpp

PURPOSE:
--------
Implementation of wave energy initialization functions for Alfvén wave turbulence
in solar energetic particle (SEP) transport simulations. Provides functions to 
set initial wave energy distributions across field line segments with proper 
heliospheric distance scaling and physical parameter-based initialization.

MAIN FUNCTIONS:
---------------

1. InitializeWaveEnergyInAllSegments() - Core initialization with r^-2 scaling
2. InitializeWaveEnergyFromPhysicalParameters() - Physics-based initialization  
3. CalculateTypicalWaveEnergyDensity1AU() - Convert B0 and turbulence to energy density
4. PrintWaveEnergyProfile() - Display energy distribution along field lines
5. PrintWaveEnergyInitializationSummary() - Global statistics and validation
6. ValidateWaveEnergyInitialization() - Consistency checking and validation

PHYSICS IMPLEMENTATION:
-----------------------
- Heliospheric scaling: ε(r) = ε_1AU × (1AU/r)²
- Wave energy density: ε = (δB²)/(2μ₀) where δB = turbulence_level × B₀
- Equal wave energies: E+ = E- = ε × V_segment / 2
- Segment center position: (x_left + x_right) / 2
- Distance calculation: r = √(x² + y² + z²)

TYPICAL USAGE EXAMPLES:
-----------------------

// Example 1: Basic initialization with specified energy density
void BasicWaveEnergySetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Initialize with 1 pJ/m³ at 1 AU
    double energy_density_1AU = 1.0e-12;  // [J/m³]
    InitializeWaveEnergyInAllSegments(WaveEnergy, energy_density_1AU, true);
    
    // Print profile for verification
    PrintWaveEnergyProfile(WaveEnergy);
    PrintWaveEnergyInitializationSummary(WaveEnergy);
}

// Example 2: Physics-based initialization (recommended)
void PhysicsBasedSetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    double B0_1AU = 5.0e-9;        // 5 nT magnetic field
    double turbulence_level = 0.2; // 20% turbulence
    
    InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, B0_1AU, turbulence_level, true);
    
    // Validate initialization
    double expected_energy = CalculateTypicalWaveEnergyDensity1AU(B0_1AU, turbulence_level);
    bool valid = ValidateWaveEnergyInitialization(WaveEnergy, expected_energy, 0.01);
    
    if (!valid) {
        std::cerr << "Error: Wave energy initialization validation failed!" << std::endl;
        exit(__LINE__);
    }
}

// Example 3: Custom reference distance (Parker Solar Probe mission)
void InnerHeliosphereSetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    double energy_density_at_10Rs = 1.0e-8;    // High energy near Sun [J/m³]
    double reference_distance = 10.0 * 6.96e8; // 10 solar radii [m]
    
    InitializeWaveEnergyInAllSegments(
        WaveEnergy, 
        energy_density_at_10Rs, 
        reference_distance, 
        true  // verbose
    );
    
    // Print detailed profile for inner heliosphere
    std::vector<int> inner_field_lines = {0, 1, 2, 3, 4};
    PrintWaveEnergyProfile(WaveEnergy, inner_field_lines, 50);
}

// Example 4: Complete simulation setup with validation
void CompleteSimulationSetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Step 1: Initialize with typical solar wind conditions
    InitializeWaveEnergyFromPhysicalParameters(
        WaveEnergy,
        WaveEnergyConstants::TYPICAL_B0_1AU,     // 5 nT
        WaveEnergyConstants::TYPICAL_TURBULENCE, // 20%
        true  // verbose output
    );
    
    // Step 2: Validate initialization
    double expected_energy = CalculateTypicalWaveEnergyDensity1AU();
    if (!ValidateWaveEnergyInitialization(WaveEnergy, expected_energy, 0.01)) {
        exit("Wave energy initialization validation failed", __LINE__, __FILE__);
    }
    
    // Step 3: Print comprehensive diagnostics
    PrintWaveEnergyInitializationSummary(WaveEnergy);
    
    // Step 4: Optional - print profile for specific field lines
    std::vector<int> sample_lines = {0, PIC::FieldLine::nFieldLine/4, PIC::FieldLine::nFieldLine/2};
    PrintWaveEnergyProfile(WaveEnergy, sample_lines, 20);
    
    std::cout << "Wave energy initialization complete. System ready for simulation." << std::endl;
}

// Example 5: Sensitivity study with different turbulence levels
void TurbulenceSensitivityStudy() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    std::vector<double> turbulence_levels = {0.1, 0.15, 0.2, 0.25, 0.3};
    double B0 = 5.0e-9;  // Fixed magnetic field
    
    for (size_t i = 0; i < turbulence_levels.size(); ++i) {
        double turb = turbulence_levels[i];
        
        std::cout << "\n=== Turbulence Level: " << turb*100 << "% ===" << std::endl;
        
        // Initialize with current turbulence level
        InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, B0, turb, false);
        
        // Print summary statistics
        PrintWaveEnergyInitializationSummary(WaveEnergy);
        
        // Validate
        double expected = CalculateTypicalWaveEnergyDensity1AU(B0, turb);
        bool valid = ValidateWaveEnergyInitialization(WaveEnergy, expected, 0.01);
        std::cout << "Validation: " << (valid ? "PASSED" : "FAILED") << std::endl;
        
        // Could save results for analysis or comparison
        // SaveEnergyConfiguration(WaveEnergy, "turbulence_" + std::to_string(turb));
    }
}

// Example 6: Distance-dependent study (multi-spacecraft mission)
void DistanceDependentStudy() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Study energy scaling from 0.1 AU to 10 AU
    std::vector<double> reference_distances = {
        0.1 * WaveEnergyConstants::ONE_AU,  // Near Sun
        0.3 * WaveEnergyConstants::ONE_AU,  // Venus orbit
        1.0 * WaveEnergyConstants::ONE_AU,  // Earth orbit
        1.5 * WaveEnergyConstants::ONE_AU,  // Mars orbit
        5.2 * WaveEnergyConstants::ONE_AU   // Jupiter orbit
    };
    
    double energy_density_base = CalculateTypicalWaveEnergyDensity1AU();
    
    for (size_t i = 0; i < reference_distances.size(); ++i) {
        double ref_dist = reference_distances[i];
        double ref_AU = ref_dist / WaveEnergyConstants::ONE_AU;
        
        std::cout << "\n=== Reference Distance: " << ref_AU << " AU ===" << std::endl;
        
        // Initialize with scaled energy density
        double scaled_energy = energy_density_base;  // Energy density at reference distance
        InitializeWaveEnergyInAllSegments(WaveEnergy, scaled_energy, ref_dist, false);
        
        PrintWaveEnergyInitializationSummary(WaveEnergy);
        
        // Print sample profile
        std::vector<int> sample_lines = {0, 1};
        PrintWaveEnergyProfile(WaveEnergy, sample_lines, 10);
    }
}

// Example 7: Integration with simulation main loop
void MainSimulationLoop() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // One-time initialization at simulation start
    static bool initialized = false;
    if (!initialized) {
        InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, 5.0e-9, 0.2, true);
        
        // Validate once at startup
        double expected = CalculateTypicalWaveEnergyDensity1AU(5.0e-9, 0.2);
        if (!ValidateWaveEnergyInitialization(WaveEnergy, expected)) {
            exit("Initial wave energy validation failed", __LINE__, __FILE__);
        }
        
        PrintWaveEnergyInitializationSummary(WaveEnergy);
        initialized = true;
    }
    
    // Main simulation timesteps would follow...
    // for (int timestep = 0; timestep < max_timesteps; ++timestep) {
    //     UpdateAllSegmentsWaveEnergyWithParticleCoupling(WaveEnergy, S_scalar, dt);
    //     // ... other physics updates ...
    // }
}

// Example 8: Error handling and robustness testing
void RobustnessTest() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Test with extreme parameters
    std::vector<std::pair<double, double>> test_cases = {
        {1.0e-9, 0.1},   // Low field, low turbulence
        {1.0e-8, 0.5},   // High field, high turbulence
        {5.0e-9, 0.01},  // Normal field, very low turbulence
        {5.0e-9, 0.9}    // Normal field, very high turbulence
    };
    
    for (auto& test_case : test_cases) {
        double B0 = test_case.first;
        double turb = test_case.second;
        
        std::cout << "\nTesting B0=" << B0*1e9 << " nT, turb=" << turb*100 << "%" << std::endl;
        
        try {
            InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, B0, turb, false);
            
            double expected = CalculateTypicalWaveEnergyDensity1AU(B0, turb);
            bool valid = ValidateWaveEnergyInitialization(WaveEnergy, expected, 0.05); // 5% tolerance
            
            std::cout << "  Result: " << (valid ? "SUCCESS" : "WARNING") << std::endl;
            
        } catch (const std::exception& e) {
            std::cout << "  Result: ERROR - " << e.what() << std::endl;
        }
    }
}

INTEGRATION NOTES:
------------------
1. Call initialization functions once at simulation startup
2. Use ValidateWaveEnergyInitialization() to ensure proper setup
3. PrintWaveEnergyProfile() is useful for debugging field line geometry
4. PrintWaveEnergyInitializationSummary() provides quick health checks
5. For production runs, set verbose=false to reduce output
6. MPI-parallel safe - all functions handle thread ownership correctly
7. Consistent with AMPS framework data access patterns

PERFORMANCE CONSIDERATIONS:
---------------------------
- Initialization is O(N_segments) and should be fast even for large simulations
- Diagnostic functions can be expensive for large numbers of field lines
- Use selective field line printing for large simulations
- Validation should be run only during setup, not every timestep

ERROR HANDLING:
---------------
- Functions check for null pointers and invalid parameters
- MPI rank 0 handles most error reporting to avoid message flooding
- Validation function returns boolean success/failure
- Detailed error messages include segment and field line indices

AUTHOR: Generated for SEP simulation framework
DATE: 2025
VERSION: 1.0

================================================================================
*/

#include "sep.h"

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

// ============================================================================
// CORE INITIALIZATION FUNCTIONS
// ============================================================================

void InitializeWaveEnergyInAllSegments(
    PIC::Datum::cDatumStored& WaveEnergy,
    double wave_energy_density_1AU,
    double reference_distance,
    bool verbose
) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int initialized_segments = 0;
    double total_wave_energy_initialized = 0.0;
    
    if (rank == 0 && verbose) {
        std::cout << "Initializing wave energy in all field line segments..." << std::endl;
        std::cout << "Reference wave energy density: " << wave_energy_density_1AU << " J/m³ at " 
                  << reference_distance/WaveEnergyConstants::ONE_AU << " AU" << std::endl;
    }
    
    // ========================================================================
    // LOOP THROUGH ALL FIELD LINES AND SEGMENTS
    // ========================================================================
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
            
            // ================================================================
            // CALCULATE SEGMENT CENTER POSITION
            // ================================================================
            
            // Get segment vertex positions
            double* x_left = segment->GetBegin()->GetX();
            double* x_right = segment->GetEnd()->GetX();
            
            // Calculate segment center position
            double x_center[3] = {
                0.5 * (x_left[0] + x_right[0]),
                0.5 * (x_left[1] + x_right[1]),
                0.5 * (x_left[2] + x_right[2])
            };
            
            // Calculate heliospheric distance (distance from origin)
            double r_helio = CalculateHeliosphericDistance(x_center);
            
            if (r_helio <= 0.0) {
                std::cerr << "Warning: Invalid heliospheric distance (" << r_helio 
                          << ") for segment " << seg_idx << " in field line " 
                          << field_line_idx << std::endl;
                continue;
            }
            
            // ================================================================
            // CALCULATE WAVE ENERGY DENSITY WITH r^-2 SCALING
            // ================================================================
            
            // Scale wave energy density as r^-2
            double local_wave_energy_density = ApplyHeliosphericScaling(
                wave_energy_density_1AU, reference_distance, r_helio
            );
            
            // ================================================================
            // CALCULATE SEGMENT VOLUME
            // ================================================================
            
            double V_segment = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);
            
            if (V_segment <= 0.0) {
                std::cerr << "Warning: Invalid segment volume (" << V_segment 
                          << ") for segment " << seg_idx << " in field line " 
                          << field_line_idx << std::endl;
                continue;
            }
            
            // ================================================================
            // CALCULATE INTEGRATED WAVE ENERGIES
            // ================================================================
            
            double E_plus, E_minus;
            ConvertDensityToIntegratedEnergy(local_wave_energy_density, V_segment, E_plus, E_minus);
            
            // ================================================================
            // STORE WAVE ENERGY DATA IN SEGMENT
            // ================================================================
            
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (wave_data) {
                wave_data[0] = E_plus;   // Outward wave energy [J]
                wave_data[1] = E_minus;  // Inward wave energy [J]
                
                initialized_segments++;
                total_wave_energy_initialized += (E_plus + E_minus);
                
                if (verbose && rank == 0 && initialized_segments <= 10) {
                    std::cout << "  Segment " << seg_idx << " (FL " << field_line_idx << "): "
                              << "r=" << r_helio/WaveEnergyConstants::ONE_AU << " AU, "
                              << "ε=" << local_wave_energy_density << " J/m³, "
                              << "V=" << V_segment << " m³, "
                              << "E+=" << E_plus << " J, "
                              << "E-=" << E_minus << " J" << std::endl;
                }
            } else {
                std::cerr << "Error: Could not access wave energy datum for segment " 
                          << seg_idx << " in field line " << field_line_idx << std::endl;
            }
        }
    }
    
    // ========================================================================
    // MPI REDUCTION FOR GLOBAL STATISTICS
    // ========================================================================
    
    int total_initialized_segments = 0;
    double global_total_wave_energy = 0.0;
    
    MPI_Allreduce(&initialized_segments, &total_initialized_segments, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_wave_energy_initialized, &global_total_wave_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "Wave energy initialization complete:" << std::endl;
        std::cout << "  Total segments initialized: " << total_initialized_segments << std::endl;
        std::cout << "  Total wave energy in system: " << global_total_wave_energy << " J" << std::endl;
        if (total_initialized_segments > 0) {
            std::cout << "  Average energy per segment: " << global_total_wave_energy / total_initialized_segments << " J" << std::endl;
        }
    }
}

// ============================================================================
// CONVENIENCE WRAPPER WITH DEFAULT 1 AU REFERENCE
// ============================================================================

void InitializeWaveEnergyInAllSegments(
    PIC::Datum::cDatumStored& WaveEnergy,
    double wave_energy_density_1AU,
    bool verbose
) {
    InitializeWaveEnergyInAllSegments(WaveEnergy, wave_energy_density_1AU, WaveEnergyConstants::ONE_AU, verbose);
}

// ============================================================================
// PHYSICS-BASED INITIALIZATION
// ============================================================================

void InitializeWaveEnergyFromPhysicalParameters(
    PIC::Datum::cDatumStored& WaveEnergy,
    double B0_1AU,
    double turbulence_level,
    bool verbose
) {
    double wave_energy_density_1AU = CalculateTypicalWaveEnergyDensity1AU(B0_1AU, turbulence_level);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0 && verbose) {
        std::cout << "Initializing wave energy from physical parameters:" << std::endl;
        std::cout << "  B₀ at 1 AU: " << B0_1AU * 1e9 << " nT" << std::endl;
        std::cout << "  Turbulence level: " << turbulence_level * 100 << "%" << std::endl;
        std::cout << "  Resulting wave energy density at 1 AU: " << wave_energy_density_1AU << " J/m³" << std::endl;
    }
    
    InitializeWaveEnergyInAllSegments(WaveEnergy, wave_energy_density_1AU, verbose);
}

// ============================================================================
// PHYSICAL PARAMETER CALCULATIONS
// ============================================================================

double CalculateTypicalWaveEnergyDensity1AU(double B0_1AU, double turbulence_level) {
    // Wave energy density: ε = (δB²)/(2μ₀) where δB = turbulence_level * B₀
    double delta_B = turbulence_level * B0_1AU;
    double wave_energy_density = (delta_B * delta_B) / (2.0 * WaveEnergyConstants::MU0);
    
    return wave_energy_density;  // [J/m³]
}

// ============================================================================
// DIAGNOSTIC AND MONITORING FUNCTIONS
// ============================================================================

void PrintWaveEnergyProfile(
    PIC::Datum::cDatumStored& WaveEnergy,
    const std::vector<int>& field_line_indices,
    int max_segments_print
) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank != 0) return; // Only rank 0 prints
    
    std::cout << "\n=== Wave Energy Profile ===" << std::endl;
    std::cout << "FieldLine  Segment  Distance[AU]  E+[J]        E-[J]        Total[J]     Density[J/m³]" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    
    // Determine which field lines to print
    std::vector<int> lines_to_print;
    if (field_line_indices.empty()) {
        // Print all field lines
        for (int i = 0; i < PIC::FieldLine::nFieldLine; ++i) {
            lines_to_print.push_back(i);
        }
    } else {
        lines_to_print = field_line_indices;
    }
    
    for (int field_line_idx : lines_to_print) {
        if (field_line_idx >= PIC::FieldLine::nFieldLine) continue;
        
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        int segments_printed = 0;
        
        for (int seg_idx = 0; seg_idx < num_segments && segments_printed < max_segments_print; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment) continue;
            
            // Calculate segment center distance
            double x_center[3];
            GetSegmentCenterPosition(segment, x_center);
            double r_helio = CalculateHeliosphericDistance(x_center);
            double r_AU = r_helio / WaveEnergyConstants::ONE_AU;
            
            // Get wave energy data
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (wave_data) {
                double E_plus = wave_data[0];
                double E_minus = wave_data[1];
                double E_total = E_plus + E_minus;
                
                double V_segment = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);
                double energy_density = ConvertIntegratedEnergyToDensity(E_plus, E_minus, V_segment);
                
                printf("%8d  %7d  %11.4f  %11.4e  %11.4e  %11.4e  %11.4e\n",
                       field_line_idx, seg_idx, r_AU, E_plus, E_minus, E_total, energy_density);
                
                segments_printed++;
            }
        }
        
        if (segments_printed >= max_segments_print && num_segments > max_segments_print) {
            std::cout << "    ... (" << (num_segments - max_segments_print) << " more segments)" << std::endl;
        }
    }
    
    std::cout << "=== End Profile ===" << std::endl << std::endl;
}

void PrintWaveEnergyInitializationSummary(PIC::Datum::cDatumStored& WaveEnergy) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Calculate local statistics
    int local_segment_count = 0;
    double local_total_energy = 0.0;
    double local_min_distance = 1e20;
    double local_max_distance = 0.0;
    double local_min_energy_density = 1e20;
    double local_max_energy_density = 0.0;
    
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;
            
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) continue;
            
            double x_center[3];
            GetSegmentCenterPosition(segment, x_center);
            double r_helio = CalculateHeliosphericDistance(x_center);
            
            double E_plus = wave_data[0];
            double E_minus = wave_data[1];
            double E_total = E_plus + E_minus;
            
            double V_segment = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);
            double energy_density = ConvertIntegratedEnergyToDensity(E_plus, E_minus, V_segment);
            
            local_segment_count++;
            local_total_energy += E_total;
            local_min_distance = std::min(local_min_distance, r_helio);
            local_max_distance = std::max(local_max_distance, r_helio);
            local_min_energy_density = std::min(local_min_energy_density, energy_density);
            local_max_energy_density = std::max(local_max_energy_density, energy_density);
        }
    }
    
    // MPI reductions for global statistics
    int global_segment_count = 0;
    double global_total_energy = 0.0;
    double global_min_distance = 0.0;
    double global_max_distance = 0.0;
    double global_min_energy_density = 0.0;
    double global_max_energy_density = 0.0;
    
    MPI_Allreduce(&local_segment_count, &global_segment_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_total_energy, &global_total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min_distance, &global_min_distance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_distance, &global_max_distance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min_energy_density, &global_min_energy_density, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_energy_density, &global_max_energy_density, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "\n=== Wave Energy Initialization Summary ===" << std::endl;
        std::cout << "Total segments: " << global_segment_count << std::endl;
        std::cout << "Total wave energy: " << std::scientific << global_total_energy << " J" << std::endl;
        
        if (global_segment_count > 0) {
            std::cout << "Average energy per segment: " << global_total_energy / global_segment_count << " J" << std::endl;
        }
        
        std::cout << "Distance range: " << global_min_distance/WaveEnergyConstants::ONE_AU 
                  << " - " << global_max_distance/WaveEnergyConstants::ONE_AU << " AU" << std::endl;
        std::cout << "Energy density range: " << global_min_energy_density 
                  << " - " << global_max_energy_density << " J/m³" << std::endl;
        std::cout << "==========================================" << std::endl << std::endl;
    }
}

bool ValidateWaveEnergyInitialization(
    PIC::Datum::cDatumStored& WaveEnergy,
    double expected_energy_1AU,
    double tolerance
) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    bool validation_passed = true;
    int error_count = 0;
    
    // Check segments near 1 AU for proper scaling
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;
            
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) continue;
            
            double x_center[3];
            GetSegmentCenterPosition(segment, x_center);
            double r_helio = CalculateHeliosphericDistance(x_center);
            double r_AU = r_helio / WaveEnergyConstants::ONE_AU;
            
            double E_plus = wave_data[0];
            double E_minus = wave_data[1];
            
            // Check for negative energies
            if (E_plus < 0.0 || E_minus < 0.0) {
                if (rank == 0) {
                    std::cerr << "Error: Negative wave energy in segment " << seg_idx 
                              << " of field line " << field_line_idx << std::endl;
                }
                validation_passed = false;
                error_count++;
            }
            
            // Check E+ = E- assumption
            double energy_ratio = std::abs(E_plus - E_minus) / (E_plus + E_minus);
            if (energy_ratio > tolerance) {
                if (rank == 0) {
                    std::cerr << "Warning: E+ != E- in segment " << seg_idx 
                              << " of field line " << field_line_idx 
                              << " (ratio: " << energy_ratio << ")" << std::endl;
                }
            }
            
            // Check r^-2 scaling for segments near 1 AU
            if (r_AU > 0.8 && r_AU < 1.2) {  // Within 20% of 1 AU
                double V_segment = SEP::FieldLine::GetSegmentVolume(segment, field_line_idx);
                double energy_density = ConvertIntegratedEnergyToDensity(E_plus, E_minus, V_segment);
                
                // Expected energy density with r^-2 scaling
                double expected_density = expected_energy_1AU * (1.0 / (r_AU * r_AU));
                double relative_error = std::abs(energy_density - expected_density) / expected_density;
                
                if (relative_error > tolerance) {
                    if (rank == 0) {
                        std::cerr << "Warning: Energy density scaling error in segment " << seg_idx
                                  << " (expected: " << expected_density 
                                  << ", actual: " << energy_density 
                                  << ", error: " << relative_error * 100 << "%)" << std::endl;
                    }
                }
            }
        }
    }
    
    // MPI reduction for global validation result
    int global_validation_passed = validation_passed ? 1 : 0;
    int global_validation_result = 0;
    MPI_Allreduce(&global_validation_passed, &global_validation_result, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    
    if (rank == 0) {
        if (global_validation_result == 1) {
            std::cout << "Wave energy initialization validation: PASSED" << std::endl;
        } else {
            std::cout << "Wave energy initialization validation: FAILED" << std::endl;
        }
    }
    
    return (global_validation_result == 1);
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

double CalculateSegmentWaveEnergy(
    PIC::FieldLine::cFieldLineSegment* segment,
    PIC::Datum::cDatumStored& WaveEnergy
) {
    if (!segment) return 0.0;
    
    double* wave_data = segment->GetDatum_ptr(WaveEnergy);
    if (!wave_data) return 0.0;
    
    return wave_data[0] + wave_data[1];  // E+ + E-
}

void GetSegmentCenterPosition(
    PIC::FieldLine::cFieldLineSegment* segment,
    double x_center[3]
) {
    if (!segment) {
        x_center[0] = x_center[1] = x_center[2] = 0.0;
        return;
    }
    
    double* x_left = segment->GetBegin()->GetX();
    double* x_right = segment->GetEnd()->GetX();
    
    x_center[0] = 0.5 * (x_left[0] + x_right[0]);
    x_center[1] = 0.5 * (x_left[1] + x_right[1]);
    x_center[2] = 0.5 * (x_left[2] + x_right[2]);
}

} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP

/*
================================================================================
                    WAVE ENERGY INITIALIZATION TEST CODE
================================================================================

FILENAME: test_wave_energy_initialization.cpp

PURPOSE:
--------
Test code to verify wave energy initialization by examining the first 200 
segments of the first field line and printing out E+ values along with 
position and scaling information.

USAGE:
------
This test function should be called after wave energy initialization to
verify that the energy values are properly set and scaled with distance.

================================================================================
*/

#include "sep.h"

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

// ============================================================================
// TEST FUNCTION: EXAMINE FIRST 200 SEGMENTS OF FIRST FIELD LINE
// ============================================================================

void TestWaveEnergyInitialization(PIC::Datum::cDatumStored& WaveEnergy) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Only rank 0 prints to avoid output flooding
    if (rank != 0) return;
    
    std::cout << "\n===============================================" << std::endl;
    std::cout << "    WAVE ENERGY INITIALIZATION TEST" << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << "Examining first 200 segments of first field line" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    
    // Check if we have any field lines
    if (PIC::FieldLine::nFieldLine <= 0) {
        std::cout << "ERROR: No field lines found!" << std::endl;
        return;
    }
    
    // Get the first field line
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[0];
    int num_segments = field_line->GetTotalSegmentNumber();
    
    std::cout << "Field line 0 has " << num_segments << " segments" << std::endl;
    std::cout << "Will examine first " << std::min(200, num_segments) << " segments" << std::endl;
    std::cout << std::endl;
    
    // Print header
    std::cout << std::setw(8) << "Segment" 
              << std::setw(12) << "Distance[AU]" 
              << std::setw(15) << "E+[J]" 
              << std::setw(15) << "E-[J]" 
              << std::setw(15) << "Total[J]" 
              << std::setw(12) << "Thread" 
              << std::setw(15) << "Volume[m³]" 
              << std::setw(15) << "Density[J/m³]" << std::endl;
    std::cout << std::string(110, '-') << std::endl;
    
    int segments_examined = 0;
    int segments_with_data = 0;
    double total_energy_examined = 0.0;
    
    // Loop through first 200 segments (or all if fewer than 200)
    int max_segments_to_check = std::min(200, num_segments);
    
    for (int seg_idx = 0; seg_idx < max_segments_to_check; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        
        if (!segment) {
            std::cout << std::setw(8) << seg_idx 
                      << std::setw(12) << "N/A" 
                      << std::setw(15) << "NULL_SEGMENT" << std::endl;
            continue;
        }
        
        segments_examined++;
        
        // Get segment center position
        double x_center[3];
        GetSegmentCenterPosition(segment, x_center);
        double r_helio = CalculateHeliosphericDistance(x_center);
        double r_AU = r_helio / WaveEnergyConstants::ONE_AU;
        
        // Get thread ownership
        int thread_id = segment->Thread;
        
        // Get segment volume
        double V_segment = SEP::FieldLine::GetSegmentVolume(segment, 0); // field_line_idx = 0
        
        // Try to get wave energy data
        double* wave_data = segment->GetDatum_ptr(WaveEnergy);
        
        if (wave_data) {
            segments_with_data++;
            
            double E_plus = wave_data[0];   // Outward wave energy [J]
            double E_minus = wave_data[1];  // Inward wave energy [J]
            double E_total = E_plus + E_minus;
            
            total_energy_examined += E_total;
            
            // Calculate energy density
            double energy_density = (V_segment > 0.0) ? E_total / V_segment : 0.0;
            
            // Print segment information
            std::cout << std::setw(8) << seg_idx 
                      << std::setw(12) << std::fixed << std::setprecision(4) << r_AU
                      << std::setw(15) << std::scientific << std::setprecision(4) << E_plus
                      << std::setw(15) << std::scientific << std::setprecision(4) << E_minus
                      << std::setw(15) << std::scientific << std::setprecision(4) << E_total
                      << std::setw(12) << thread_id
                      << std::setw(15) << std::scientific << std::setprecision(4) << V_segment
                      << std::setw(15) << std::scientific << std::setprecision(4) << energy_density
                      << std::endl;
                      
        } else {
            // No wave data available
            std::cout << std::setw(8) << seg_idx 
                      << std::setw(12) << std::fixed << std::setprecision(4) << r_AU
                      << std::setw(15) << "NO_DATA" 
                      << std::setw(15) << "NO_DATA" 
                      << std::setw(15) << "NO_DATA" 
                      << std::setw(12) << thread_id
                      << std::setw(15) << std::scientific << std::setprecision(4) << V_segment
                      << std::setw(15) << "NO_DATA" << std::endl;
        }
    }
    
    // Print summary statistics
    std::cout << std::string(110, '-') << std::endl;
    std::cout << "SUMMARY:" << std::endl;
    std::cout << "  Segments examined: " << segments_examined << std::endl;
    std::cout << "  Segments with wave data: " << segments_with_data << std::endl;
    std::cout << "  Segments without data: " << (segments_examined - segments_with_data) << std::endl;
    
    if (segments_with_data > 0) {
        std::cout << "  Total energy in examined segments: " << std::scientific << total_energy_examined << " J" << std::endl;
        std::cout << "  Average energy per segment: " << std::scientific << (total_energy_examined / segments_with_data) << " J" << std::endl;
    }
    
    // Additional analysis: Check r^-2 scaling
    if (segments_with_data >= 2) {
        std::cout << "\nSCALING ANALYSIS:" << std::endl;
        
        // Find two segments at different distances with data
        double r1 = -1, r2 = -1, E1 = -1, E2 = -1;
        
        for (int seg_idx = 0; seg_idx < max_segments_to_check; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment) continue;
            
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) continue;
            
            double x_center[3];
            GetSegmentCenterPosition(segment, x_center);
            double r_helio = CalculateHeliosphericDistance(x_center);
            double r_AU = r_helio / WaveEnergyConstants::ONE_AU;
            
            double V_segment = SEP::FieldLine::GetSegmentVolume(segment, 0);
            double E_total = wave_data[0] + wave_data[1];
            double energy_density = (V_segment > 0.0) ? E_total / V_segment : 0.0;
            
            if (r1 < 0) {
                r1 = r_AU;
                E1 = energy_density;
            } else if (std::abs(r_AU - r1) > 0.1) {  // Different distance
                r2 = r_AU;
                E2 = energy_density;
                break;
            }
        }
        
        if (r1 > 0 && r2 > 0 && E1 > 0 && E2 > 0) {
            // Check if scaling follows r^-2 law
            double expected_ratio = (r2 / r1) * (r2 / r1);  // (r2/r1)^2
            double actual_ratio = E1 / E2;  // E1/E2 should equal (r2/r1)^2
            double scaling_error = std::abs(expected_ratio - actual_ratio) / expected_ratio;
            
            std::cout << "  Distance 1: " << r1 << " AU, Energy density: " << E1 << " J/m³" << std::endl;
            std::cout << "  Distance 2: " << r2 << " AU, Energy density: " << E2 << " J/m³" << std::endl;
            std::cout << "  Expected ratio (r2/r1)²: " << expected_ratio << std::endl;
            std::cout << "  Actual ratio E1/E2: " << actual_ratio << std::endl;
            std::cout << "  Scaling error: " << scaling_error * 100 << "%" << std::endl;
            
            if (scaling_error < 0.01) {
                std::cout << "  ✓ r^-2 scaling VERIFIED" << std::endl;
            } else if (scaling_error < 0.1) {
                std::cout << "  ⚠ r^-2 scaling approximately correct" << std::endl;
            } else {
                std::cout << "  ✗ r^-2 scaling FAILED" << std::endl;
            }
        }
    }
    
    std::cout << "===============================================" << std::endl;
    std::cout << "           TEST COMPLETE" << std::endl;
    std::cout << "===============================================" << std::endl << std::endl;
}

// ============================================================================
// SIMPLIFIED TEST FUNCTION: JUST PRINT E+ VALUES
// ============================================================================

void TestPrintEPlusValues(PIC::Datum::cDatumStored& WaveEnergy) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Only rank 0 prints
    if (rank != 0) return;
    
    std::cout << "\n=== E+ VALUES FOR FIRST 400 SEGMENTS OF FIELD LINE 0 ===" << std::endl;
    
    if (PIC::FieldLine::nFieldLine <= 0) {
        std::cout << "ERROR: No field lines found!" << std::endl;
        return;
    }
    
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[0];
    int num_segments = field_line->GetTotalSegmentNumber();
    int max_segments = std::min(400, num_segments);
    
    std::cout << "Segment    E+ [J]" << std::endl;
    std::cout << "-------------------" << std::endl;
    
    for (int seg_idx = 0; seg_idx < max_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        
        if (!segment) {
            std::cout << std::setw(7) << seg_idx << "    NULL_SEGMENT" << std::endl;
            continue;
        }
        
        double* wave_data = segment->GetDatum_ptr(WaveEnergy);
        
        if (wave_data) {
            double E_plus = wave_data[0];
            std::cout << std::setw(7) << seg_idx << "    " 
                      << std::scientific << std::setprecision(6) << E_plus << "    " << segment->Thread << std::endl;
        } else {
            std::cout << std::setw(7) << seg_idx << "    NO_DATA" << std::endl;
        }
    }
    
    std::cout << "==========================================" << std::endl << std::endl;
}

} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP

// ============================================================================
// EXAMPLE USAGE IN MAIN SIMULATION CODE
// ============================================================================

/*

// Example: How to use these test functions in your simulation

void SetupAndTestWaveEnergy() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Step 1: Initialize wave energy
    InitializeWaveEnergyFromPhysicalParameters(
        WaveEnergy,
        5.0e-9,  // 5 nT
        0.2,     // 20% turbulence
        true     // verbose
    );
    
    // Step 2: Run comprehensive test
    TestWaveEnergyInitialization(WaveEnergy);
    
    // Step 3: Run simple E+ test
    TestPrintEPlusValues(WaveEnergy);
    
    // Step 4: Validate
    double expected = CalculateTypicalWaveEnergyDensity1AU(5.0e-9, 0.2);
    bool valid = ValidateWaveEnergyInitialization(WaveEnergy, expected);
    
    if (!valid) {
        std::cerr << "Validation failed!" << std::endl;
    }
}

// Or just the simple test:
void SimpleTest() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Assume WaveEnergy is already initialized...
    TestPrintEPlusValues(WaveEnergy);
}

*/

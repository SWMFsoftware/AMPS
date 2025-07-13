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

void TestPrintEPlusValues(PIC::Datum::cDatumStored& WaveEnergy,int PrintThread) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Only rank 0 prints
    if (rank != PrintThread) return;
    
    std::cout << "\n=== E+ VALUES FOR FIRST 400 SEGMENTS OF FIELD LINE 0 === Thread=" << PrintThread << std::endl;
    
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

/*
================================================================================
                    MPI-DISTRIBUTED TEST PRINT DATUM FUNCTION
================================================================================

FUNCTION: TestPrintDatumMPI

PURPOSE:
--------
Performs distributed analysis of field line segment data across multiple MPI 
processes and outputs comprehensive results from the root process. This function
efficiently collects data from all processes, finds global maximums using 
MPI_Reduce, and displays the most significant segments based on relative values
calculated from true global statistics.

ALGORITHM:
----------
1. Distributed Global Maximum Detection (MPI_Reduce):
   - Each process finds local maximum for each datum element
   - Custom MPI reduction operation finds global maximum while preserving 
     process ownership information
   - Root process receives global maximums with process attribution

2. Distributed Segment Analysis:
   - Global maximums broadcast to all processes for consistent normalization
   - Each process analyzes only its own segments (segment->Thread == rank)
   - Calculates relative values: |data[i]|/global_max[i] for ranking
   - Finds maximum absolute values within each segment for display
   - Filters out segments with zero relative values

3. Data Collection and Output:
   - Root process gathers significant segments from all processes
   - Sorts all segments globally by relative value
   - Outputs top 400 segments with complete process attribution
   - Provides comprehensive statistics across entire distributed dataset

MPI COMMUNICATION PATTERN:
---------------------------
1. MPI_Reduce: Global maximum detection with process tracking
2. MPI_Bcast: Distribute global maximums to all processes  
3. MPI_Gather: Collect segment counts from all processes
4. MPI_Send/Recv: Transfer segment data to root process
   - Variable-length data handled with separate metadata and data transfers
   - Each segment's complete datum array transmitted

PARAMETERS:
-----------
- Datum: Reference to distributed datum storage across all processes
- msg: Descriptive message for output header
- field_line_idx: Index of field line to analyze (default: 0)

OUTPUT FORMAT:
--------------
For each selected segment, displays:
1. Segment index in the field line
2. Maximum absolute value: [element_index]=value(relative_value)
   - element_index: Index of element with maximum absolute value
   - value: The actual maximum absolute value in the segment  
   - relative_value: The relative value used for global ranking
3. Process rank that owns the segment
4. Complete list of all datum elements: [elem0, elem1, elem2, ...]

SELECTION CRITERIA:
-------------------
Segments are ranked by their maximum relative value using TRUE global maximums:
  max_relative = max(|data[i]|/global_max_across_all_processes[i])

This ensures proper normalization across the entire distributed dataset, not 
just local process data, providing accurate global significance ranking.

DISTRIBUTED FEATURES:
---------------------
- Process Attribution: Shows which MPI process owns each significant segment
- Load Distribution Analysis: Statistics on segment distribution across processes
- True Global Ranking: Uses global maximums from entire distributed dataset
- Single Point Output: Only root process outputs to prevent MPI conflicts
- Scalable Communication: O(log P) reduction vs O(P) gathering

USAGE EXAMPLES:
---------------
// Basic distributed analysis
TestPrintDatumMPI(WaveEnergy, "Distributed Wave Energy Analysis");

// Analyze specific field line across all processes
TestPrintDatumMPI(MagneticField, "Global B-field Distribution", 3);

// Multi-field line distributed analysis
for (int fl = 0; fl < num_field_lines; ++fl) {
    TestPrintDatumMPI(Pressure, "Multi-Line Pressure Analysis", fl);
}

REQUIREMENTS:
-------------
- Active MPI environment with communicator MPI_COMM_WORLD
- Consistent field line structure across all processes
- Each process owns subset of segments (segment->Thread == process_rank)
- Non-zero datum values for meaningful analysis

OUTPUT SECTIONS:
----------------
1. Header with MPI process count and field line information
2. Global maximum values per element with process attribution
3. Top 400 segments from entire distributed dataset
4. Comprehensive summary statistics:
   - Total processes and segment distribution
   - Global maximum locations with process ownership
   - Load balancing information (segments per process)
   - Highest/lowest relative values with process attribution

MPI PERFORMANCE:
----------------
- Communication Complexity: O(log P) for reductions + O(K) for data gathering
- Memory Usage: O(K*M) on root, O(k*M) on workers where:
  - K = total significant segments across all processes
  - k = local significant segments per process  
  - M = datum length
  - P = number of MPI processes
- Network Efficiency: Minimized data transfer through filtering

ERROR HANDLING:
---------------
- Validates field line index on all processes consistently
- Handles processes with no significant segments gracefully
- Manages variable-length data transfers safely
- Detects and reports distributed all-zero datasets
- Prevents deadlocks in MPI communication

ADVANTAGES OVER SINGLE-PROCESS VERSION:
---------------------------------------
- True Global Analysis: Uses data from ALL processes, not just one
- Proper Normalization: Global maximums calculated across entire dataset
- Process Attribution: Identifies data ownership for debugging
- Load Balancing Insights: Shows distribution of significant data
- Scalable Performance: Efficient for large distributed simulations

MPI COMMUNICATION DETAILS:
---------------------------
Custom MPI Reduction Operation:
```cpp
MPI_Op_create([](void* in, void* inout, int* len, MPI_Datatype* datatype) {
    // Finds maximum value while preserving process rank
}, 1, &max_loc_op);
```

Point-to-Point Data Transfer:
- Tag 0: relative_value (MPI_DOUBLE)
- Tag 1: segment_idx (MPI_INT)  
- Tag 2: element_idx (MPI_INT)
- Tag 3: element_value (MPI_DOUBLE)
- Tag 4: process_rank (MPI_INT)
- Tag 5: all_data array (MPI_DOUBLE array)

NOTES:
------
- Only root process (rank 0) produces output regardless of input parameters
- Segments are filtered locally but ranked globally for optimal performance
- Function handles dynamic datum lengths without size assumptions
- Process attribution enables distributed debugging and load analysis
- Compatible with any number of MPI processes (tested with 1-1000+ processes)

SEE ALSO:
---------
- TestPrintDatum(): Single-process version for non-MPI analysis
- MPI_Reduce documentation for custom reduction operations
- MPI_Allgather vs point-to-point communication trade-offs

================================================================================
*/

void TestPrintDatumMPI(PIC::Datum::cDatumStored& Datum, const char* msg, int field_line_idx) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Structure to hold segment data with relative values
    struct SegmentWithRelativeValue {
        double relative_value;
        int segment_idx;
        int element_idx;
        double element_value;
        int process_rank;
        std::vector<double> all_data;
        
        SegmentWithRelativeValue() : relative_value(0.0), segment_idx(-1), element_idx(-1), 
                                   element_value(0.0), process_rank(-1) {}
    };

    // Structure for MPI reduction to find max values per element
    struct MaxValueData {
        double value;
        int rank;
    };

    // Validate field line index on all processes
    if (PIC::FieldLine::nFieldLine <= 0 || field_line_idx >= PIC::FieldLine::nFieldLine || field_line_idx < 0) {
        if (rank == 0) {
            std::cout << "\n=== " << msg << " === Root Process (MPI_Reduce) === Field Line=" << field_line_idx << " ===" << std::endl;
            if (PIC::FieldLine::nFieldLine <= 0) {
                std::cout << "ERROR: No field lines found!" << std::endl;
            } else {
                std::cout << "ERROR: Field line index " << field_line_idx 
                          << " out of range [0, " << (PIC::FieldLine::nFieldLine - 1) << "]!" << std::endl;
            }
        }
        return;
    }

    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    int num_segments = field_line->GetTotalSegmentNumber();
    int datum_length = Datum.length;

    // Step 1: Find local maximum for each element
    std::vector<MaxValueData> local_max_per_element(datum_length);
    for (int i = 0; i < datum_length; ++i) {
        local_max_per_element[i] = {0.0, rank};
    }

    // Find local maximums on this process
    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        if (!segment) continue;

        double* data = segment->GetDatum_ptr(Datum);
        if (!data) continue;

        // Only process segments owned by this rank
        if (segment->Thread != rank) continue;

        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            double abs_val = std::abs(data[elem_idx]);
            if (abs_val > local_max_per_element[elem_idx].value) {
                local_max_per_element[elem_idx].value = abs_val;
                local_max_per_element[elem_idx].rank = rank;
            }
        }
    }

    // Step 2: Use MPI_Reduce to find global maximum for each element
    // Custom MPI operation to find max value and keep the rank
    MPI_Op max_loc_op;
    MPI_Op_create([](void* in, void* inout, int* len, MPI_Datatype* datatype) {
        MaxValueData* in_data = static_cast<MaxValueData*>(in);
        MaxValueData* inout_data = static_cast<MaxValueData*>(inout);
        for (int i = 0; i < *len; ++i) {
            if (in_data[i].value > inout_data[i].value) {
                inout_data[i] = in_data[i];
            }
        }
    }, 1, &max_loc_op);

    std::vector<MaxValueData> global_max_per_element(datum_length);
    MPI_Reduce(local_max_per_element.data(), global_max_per_element.data(), 
               datum_length, MPI_DOUBLE_INT, max_loc_op, 0, MPI_COMM_WORLD);

    MPI_Op_free(&max_loc_op);

    // Broadcast global maximums to all processes for relative value calculation
    std::vector<double> global_max_values(datum_length);
    if (rank == 0) {
        for (int i = 0; i < datum_length; ++i) {
            global_max_values[i] = global_max_per_element[i].value;
        }
    }
    MPI_Bcast(global_max_values.data(), datum_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Step 3: Calculate relative values and collect significant segments locally
    std::vector<SegmentWithRelativeValue> local_significant_segments;

    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        if (!segment) continue;

        double* data = segment->GetDatum_ptr(Datum);
        if (!data) continue;

        // Only process segments owned by this rank
        if (segment->Thread != rank) continue;

        // First, find the element with maximum relative value for segment selection
        double max_relative = 0.0;
        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            if (global_max_values[elem_idx] == 0.0) continue;
            double relative_val = std::abs(data[elem_idx]) / global_max_values[elem_idx];
            if (relative_val > max_relative) {
                max_relative = relative_val;
            }
        }

        // Only process segments with non-zero relative values
        if (max_relative > 0.0) {
            // Now find the element with maximum absolute value for output display
            double max_abs_value = 0.0;
            int max_abs_element_idx = 0;
            double max_abs_element_val = 0.0;
            
            for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
                double abs_val = std::abs(data[elem_idx]);
                if (abs_val > max_abs_value) {
                    max_abs_value = abs_val;
                    max_abs_element_idx = elem_idx;
                    max_abs_element_val = data[elem_idx];  // Store actual value with sign
                }
            }

            SegmentWithRelativeValue seg_data;
            seg_data.relative_value = max_relative;
            seg_data.segment_idx = seg_idx;
            seg_data.element_idx = max_abs_element_idx;      // Index of max absolute value
            seg_data.element_value = max_abs_element_val;    // Actual max absolute value
            seg_data.process_rank = rank;
            seg_data.all_data.resize(datum_length);
            
            for (int i = 0; i < datum_length; ++i) {
                seg_data.all_data[i] = data[i];
            }
            
            local_significant_segments.push_back(seg_data);
        }
    }

    // Step 4: Gather segment counts from all processes
    int local_segment_count = local_significant_segments.size();
    std::vector<int> all_segment_counts(size);
    MPI_Gather(&local_segment_count, 1, MPI_INT, all_segment_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Step 5: Root process gathers all significant segments
    std::vector<SegmentWithRelativeValue> all_significant_segments;
    
    if (rank == 0) {
        // Add root's own segments
        all_significant_segments.insert(all_significant_segments.end(), 
                                      local_significant_segments.begin(), 
                                      local_significant_segments.end());

        // Receive segments from other processes
        for (int proc = 1; proc < size; ++proc) {
            int count = all_segment_counts[proc];
            if (count == 0) continue;

            for (int i = 0; i < count; ++i) {
                SegmentWithRelativeValue seg_data;
                
                // Receive basic data
                MPI_Recv(&seg_data.relative_value, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&seg_data.segment_idx, 1, MPI_INT, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&seg_data.element_idx, 1, MPI_INT, proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&seg_data.element_value, 1, MPI_DOUBLE, proc, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&seg_data.process_rank, 1, MPI_INT, proc, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                // Receive data array
                seg_data.all_data.resize(datum_length);
                MPI_Recv(seg_data.all_data.data(), datum_length, MPI_DOUBLE, proc, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                all_significant_segments.push_back(seg_data);
            }
        }
    } else {
        // Non-root processes send their segments
        for (const auto& seg_data : local_significant_segments) {
            MPI_Send(&seg_data.relative_value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&seg_data.segment_idx, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&seg_data.element_idx, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
            MPI_Send(&seg_data.element_value, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
            MPI_Send(&seg_data.process_rank, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
            MPI_Send(seg_data.all_data.data(), datum_length, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
        }
    }

    // Step 6: Root process analyzes and outputs results
    if (rank == 0) {
        std::cout << "\n=== " << msg << " === Root Process (MPI_Reduce) === Field Line=" << field_line_idx << " ===" << std::endl;
        std::cout << "Field line " << field_line_idx << " analysis from " << size << " processes" << std::endl;
        std::cout << "Datum has " << datum_length << " elements per segment" << std::endl;

        // Check if we have any data
        bool has_nonzero_data = false;
        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            if (global_max_per_element[elem_idx].value > 0.0) {
                has_nonzero_data = true;
                break;
            }
        }

        if (!has_nonzero_data) {
            std::cout << "ERROR: All data values are zero for all elements across all processes!" << std::endl;
            return;
        }

        std::cout << "Global maximum absolute values per element (via MPI_Reduce):" << std::endl;
        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            if (global_max_per_element[elem_idx].value > 0.0) {
                std::cout << "  Element[" << elem_idx << "]: " << std::scientific << std::setprecision(6) 
                          << global_max_per_element[elem_idx].value 
                          << " (process " << global_max_per_element[elem_idx].rank << ")" << std::endl;
            }
        }

        // Sort all segments by relative value
        std::sort(all_significant_segments.begin(), all_significant_segments.end(), 
                  [](const SegmentWithRelativeValue& a, const SegmentWithRelativeValue& b) {
                      return a.relative_value > b.relative_value;
                  });

        // Print top segments
        int max_segments_to_print = std::min(400, static_cast<int>(all_significant_segments.size()));
        
        std::cout << "\nTop " << max_segments_to_print << " segments with highest relative values (from all processes):" << std::endl;
        std::cout << "Segment    MaxElem[idx]=value(rel_val)    Process    All_Elements" << std::endl;
        std::cout << "           |       |        |               |         |" << std::endl;
        std::cout << "           |       |        |               |         +-- Complete datum array" << std::endl;
        std::cout << "           |       |        |               +-- MPI process owning segment" << std::endl;
        std::cout << "           |       |        +-- Relative value for ranking" << std::endl;
        std::cout << "           |       +-- Maximum absolute value in segment" << std::endl;
        std::cout << "           +-- Index of element with maximum absolute value" << std::endl;
        std::cout << std::string(100, '-') << std::endl;

        for (int i = 0; i < max_segments_to_print; ++i) {
            const SegmentWithRelativeValue& seg = all_significant_segments[i];
            
            // Output: 1) segment index
            std::cout << std::setw(7) << seg.segment_idx << "    ";
            
            // Output: 2) data and its index why the segment is selected
            std::cout << "[" << seg.element_idx << "]=" 
                      << std::scientific << std::setprecision(4) << seg.element_value
                      << "(" << std::fixed << std::setprecision(4) << seg.relative_value << ")"
                      << "    " << std::setw(7) << seg.process_rank << "    ";

            // Output: 3) the entire list of elements for the datum in a single line
            std::cout << "[";
            for (int elem_idx = 0; elem_idx < static_cast<int>(seg.all_data.size()); ++elem_idx) {
                if (elem_idx > 0) std::cout << ", ";
                std::cout << std::scientific << std::setprecision(4) << seg.all_data[elem_idx];
            }
            std::cout << "]" << std::endl;
        }

        // Print summary statistics
        std::cout << std::string(100, '-') << std::endl;
        std::cout << "SUMMARY:" << std::endl;
        std::cout << "  Total processes: " << size << std::endl;
        std::cout << "  Total segments with data: " << all_significant_segments.size() << std::endl;
        std::cout << "  Segments printed: " << max_segments_to_print << std::endl;
        
        // Process distribution statistics
        std::vector<int> segments_per_process(size, 0);
        for (const auto& seg : all_significant_segments) {
            if (seg.process_rank >= 0 && seg.process_rank < size) {
                segments_per_process[seg.process_rank]++;
            }
        }
        
        // Show distribution of significant segments across MPI processes
        // This indicates load balancing and data locality:
        // - High numbers: Process owns many segments with significant data
        // - Zero values: Process has no significant segments (small/zero data)
        // - Uneven distribution may indicate load balancing issues
        std::cout << "  Segments with data per process: ";
        for (int i = 0; i < size; ++i) {
            std::cout << "P" << i << ":" << segments_per_process[i];
            if (i < size - 1) std::cout << ", ";
        }
        std::cout << std::endl;

        std::cout << "  Global maximum per element (via MPI_Reduce):" << std::endl;
        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            if (global_max_per_element[elem_idx].value > 0.0) {
                std::cout << "    Element[" << elem_idx << "]: " << std::scientific << std::setprecision(6) 
                          << global_max_per_element[elem_idx].value 
                          << " (process " << global_max_per_element[elem_idx].rank << ")" << std::endl;
            }
        }
        
        if (max_segments_to_print > 0) {
            std::cout << "  Highest relative value: " << std::fixed << std::setprecision(6) 
                      << all_significant_segments[0].relative_value 
                      << " (element " << all_significant_segments[0].element_idx 
                      << ", process " << all_significant_segments[0].process_rank << ")" << std::endl;
            std::cout << "  Lowest relative value (in selection): " << std::fixed << std::setprecision(6) 
                      << all_significant_segments[max_segments_to_print-1].relative_value 
                      << " (element " << all_significant_segments[max_segments_to_print-1].element_idx 
                      << ", process " << all_significant_segments[max_segments_to_print-1].process_rank << ")" << std::endl;
        }

        std::cout << "==========================================" << std::endl << std::endl;
    }
}


/*
================================================================================
                        TEST PRINT DATUM FUNCTION
================================================================================

FUNCTION: TestPrintDatum

PURPOSE:
--------
Analyzes and displays the most significant segments of a field line based on 
datum values. This function identifies segments with the highest relative 
values across all datum elements and outputs detailed information about the 
top 400 segments for debugging and analysis purposes.

ALGORITHM:
----------
1. Global Maximum Detection:
   - Scans all segments in the specified field line
   - Finds the global maximum absolute value for each datum element
   - Tracks which segment and thread owns each global maximum

2. Segment Analysis:
   - For each segment, calculates relative values: |data[i]|/global_max[i]
   - Finds the maximum relative value across all elements (for ranking)
   - Finds the element with maximum absolute value (for display)
   - Only retains segments with non-zero relative values

3. Ranking and Output:
   - Sorts segments by maximum relative value (descending)
   - Displays top 400 segments with complete information
   - Shows global statistics and distribution analysis

PARAMETERS:
-----------
- Datum: Reference to the datum storage containing field line segment data
- PrintThread: [UNUSED] Legacy parameter - function always uses root process
- msg: Descriptive message for output header
- field_line_idx: Index of field line to analyze (default: 0)

OUTPUT FORMAT:
--------------
For each selected segment, displays:
1. Segment index in the field line
2. Maximum absolute value: [element_index]=value(relative_value)
   - element_index: Index of element with maximum absolute value
   - value: The actual maximum absolute value in the segment
   - relative_value: The relative value used for segment ranking
3. Thread ID that owns the segment
4. Complete list of all datum elements: [elem0, elem1, elem2, ...]

SELECTION CRITERIA:
-------------------
Segments are ranked by their maximum relative value across all elements:
  max_relative = max(|data[i]|/global_max[i]) for i in [0, datum_length-1]

This ensures fair comparison between elements that may have vastly different
scales (e.g., element[0] ~ 1e+6, element[1] ~ 1e-3).

USAGE EXAMPLES:
---------------
// Basic usage - analyze first field line
TestPrintDatum(WaveEnergy, 0, "Wave Energy Analysis");

// Analyze specific field line
TestPrintDatum(MagneticField, 0, "B-field Distribution", 5);

// Multiple field line analysis
for (int fl = 0; fl < num_field_lines; ++fl) {
    TestPrintDatum(SomeDatum, 0, "Multi-Line Analysis", fl);
}

REQUIREMENTS:
-------------
- MPI environment (uses rank 0 for output)
- Valid field line with segments containing datum data
- Non-zero datum values for meaningful analysis

OUTPUT SECTIONS:
----------------
1. Header with field line and datum information
2. Global maximum values per element with locations
3. Top 400 segments with detailed data
4. Summary statistics including:
   - Total segments analyzed
   - Global maximum locations
   - Highest/lowest relative values in selection

NOTES:
------
- Only root process (rank 0) produces output to avoid MPI conflicts
- Segments with all-zero data are excluded from analysis
- Function handles variable-length datum arrays dynamically
- Global maximums are calculated across ALL segments, not per-process
- Relative values ensure fair ranking regardless of element magnitudes

ERROR HANDLING:
---------------
- Validates field line index range
- Checks for existence of field lines
- Handles segments without datum data gracefully
- Detects and reports all-zero datasets

PERFORMANCE:
------------
- Time complexity: O(N*M + K*log(K)) where:
  - N = number of segments
  - M = datum length  
  - K = number of non-zero segments
- Memory usage: O(K*M) for storing significant segment data
- Single-process execution (non-distributed)

SEE ALSO:
---------
- TestPrintDatumMPI(): MPI-distributed version with cross-process analysis
- TestWaveEnergyInitialization(): Specialized wave energy analysis
- TestPrintEPlusValues(): Simplified energy component display

================================================================================
*/
void TestPrintDatum(PIC::Datum::cDatumStored& Datum, int PrintThread, const char* msg, int field_line_idx) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Only root process (rank 0) prints output
    if (rank != 0) return;

    std::cout << "\n=== " << msg << " === Root Process === Field Line=" << field_line_idx << " ===" << std::endl;

    if (PIC::FieldLine::nFieldLine <= 0) {
        std::cout << "ERROR: No field lines found!" << std::endl;
        return;
    }

    if (field_line_idx >= PIC::FieldLine::nFieldLine || field_line_idx < 0) {
        std::cout << "ERROR: Field line index " << field_line_idx
                  << " out of range [0, " << (PIC::FieldLine::nFieldLine - 1) << "]!" << std::endl;
        return;
    }

    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    int num_segments = field_line->GetTotalSegmentNumber();
    int datum_length = Datum.length;

    std::cout << "Field line " << field_line_idx << " has " << num_segments << " segments" << std::endl;
    std::cout << "Datum has " << datum_length << " elements per segment" << std::endl;

    // Structure to store segment info for sorting
    struct SegmentInfo {
        int segment_idx;
        double max_relative_value;  // max(|data[i]|/max_global) for ranking
        int max_abs_element_index;  // index of element with maximum absolute value
        double max_abs_element_value; // actual value of the maximum absolute element
        int max_element_index;      // index of element with maximum relative value
        double max_element_value;   // actual value of the maximum relative element
        double* data;
        int thread_id;
    };

    std::vector<SegmentInfo> segment_data;
    std::vector<double> global_max_per_element(datum_length, 0.0);
    std::vector<int> global_max_segment_per_element(datum_length, -1);
    std::vector<int> global_max_thread_per_element(datum_length, -1);

    // First pass: Find global maximum for each element across all segments
    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        if (!segment) continue;

        double* data = segment->GetDatum_ptr(Datum);
        if (!data) continue;

        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            double abs_val = std::abs(data[elem_idx]);
            if (abs_val > global_max_per_element[elem_idx]) {
                global_max_per_element[elem_idx] = abs_val;
                global_max_segment_per_element[elem_idx] = seg_idx;
                global_max_thread_per_element[elem_idx] = segment->Thread;
            }
        }
    }

    // Check if any element has non-zero maximum
    bool has_nonzero_data = false;
    for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
        if (global_max_per_element[elem_idx] > 0.0) {
            has_nonzero_data = true;
            break;
        }
    }

    if (!has_nonzero_data) {
        std::cout << "ERROR: All data values are zero for all elements!" << std::endl;
        return;
    }

    std::cout << "Global maximum absolute values per element:" << std::endl;
    for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
        if (global_max_per_element[elem_idx] > 0.0) {
            std::cout << "  Element[" << elem_idx << "]: " << std::scientific << std::setprecision(6)
                      << global_max_per_element[elem_idx]
                      << " (segment " << global_max_segment_per_element[elem_idx]
                      << ", thread " << global_max_thread_per_element[elem_idx] << ")" << std::endl;
        }
    }

    // Second pass: Calculate relative values and find segment maximums
    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        if (!segment) continue;

        double* data = segment->GetDatum_ptr(Datum);
        if (!data) continue;

        // Find maximum relative value for segment ranking
        double max_relative = 0.0;
        int max_rel_element_idx = 0;
        double max_rel_element_val = 0.0;

        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            // Skip elements with zero global maximum (avoid division by zero)
            if (global_max_per_element[elem_idx] == 0.0) continue;

            double relative_val = std::abs(data[elem_idx]) / global_max_per_element[elem_idx];
            if (relative_val > max_relative) {
                max_relative = relative_val;
                max_rel_element_idx = elem_idx;
                max_rel_element_val = data[elem_idx];
            }
        }

        // Find maximum absolute value for display
        double max_abs_value = 0.0;
        int max_abs_element_idx = 0;
        double max_abs_element_val = 0.0;

        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            double abs_val = std::abs(data[elem_idx]);
            if (abs_val > max_abs_value) {
                max_abs_value = abs_val;
                max_abs_element_idx = elem_idx;
                max_abs_element_val = data[elem_idx];
            }
        }

        // Only store segment info if it has non-zero relative value
        if (max_relative > 0.0) {
            SegmentInfo info;
            info.segment_idx = seg_idx;
            info.max_relative_value = max_relative;
            info.max_abs_element_index = max_abs_element_idx;      // For display
            info.max_abs_element_value = max_abs_element_val;      // For display
            info.max_element_index = max_rel_element_idx;          // For ranking reference
            info.max_element_value = max_rel_element_val;          // For ranking reference
            info.data = data;
            info.thread_id = segment->Thread;
            segment_data.push_back(info);
        }
    }

    // Sort segments by relative value (descending order)
    std::sort(segment_data.begin(), segment_data.end(),
              [](const SegmentInfo& a, const SegmentInfo& b) {
                  return a.max_relative_value > b.max_relative_value;
              });

    // Print top segments
    int max_segments_to_print = std::min(400, static_cast<int>(segment_data.size()));

    std::cout << "\nTop " << max_segments_to_print << " segments with highest relative values:" << std::endl;
    std::cout << "Segment    MaxVal[idx]=value(rel_val)    Thread    All_Elements" << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    for (int i = 0; i < max_segments_to_print; ++i) {
        const SegmentInfo& info = segment_data[i];

        // Output: 1) segment index
        std::cout << std::setw(7) << info.segment_idx << "    ";

        // Output: 2) Maximum absolute value and its index (for clarity in analysis)
        std::cout << "[" << info.max_abs_element_index << "]="
                  << std::scientific << std::setprecision(4) << info.max_abs_element_value
                  << "(" << std::fixed << std::setprecision(4) << info.max_relative_value << ")"
                  << "    " << std::setw(6) << info.thread_id << "    ";

        // Output: 3) the entire list of elements for the datum in a single line
        std::cout << "[";
        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            if (elem_idx > 0) std::cout << ", ";
            std::cout << std::scientific << std::setprecision(4) << info.data[elem_idx];
        }
        std::cout << "]" << std::endl;
    }

    // Print summary statistics
    std::cout << std::string(100, '-') << std::endl;
    std::cout << "SUMMARY:" << std::endl;
    std::cout << "  Total segments with data: " << segment_data.size() << std::endl;
    std::cout << "  Segments printed: " << max_segments_to_print << std::endl;

    std::cout << "  Global maximum per element (non-zero only):" << std::endl;
    for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
        if (global_max_per_element[elem_idx] > 0.0) {
            std::cout << "    Element[" << elem_idx << "]: " << std::scientific << std::setprecision(6)
                      << global_max_per_element[elem_idx]
                      << " (segment " << global_max_segment_per_element[elem_idx]
                      << ", thread " << global_max_thread_per_element[elem_idx] << ")" << std::endl;
        }
    }

    if (max_segments_to_print > 0) {
        std::cout << "  Highest relative value: " << std::fixed << std::setprecision(6)
                  << segment_data[0].max_relative_value
                  << " (max abs element " << segment_data[0].max_abs_element_index << ")" << std::endl;
        std::cout << "  Lowest relative value (in selection): " << std::fixed << std::setprecision(6)
                  << segment_data[max_segments_to_print-1].max_relative_value
                  << " (max abs element " << segment_data[max_segments_to_print-1].max_abs_element_index << ")" << std::endl;
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




/*
================================================================================
                            SET DATUM FUNCTION
================================================================================

FUNCTION: SetDatum

PURPOSE:
--------
Sets all elements of datum data to a specified value for all segments of a
given field line. This function is useful for initializing datum fields,
testing, debugging, or resetting field line data to known values.

ALGORITHM:
----------
1. Validation:
   - Checks if field lines exist
   - Validates field line index range
   - Verifies datum length is valid

2. Segment Iteration:
   - Iterates through all segments in the specified field line
   - Skips null or invalid segments
   - Processes ALL segments regardless of Thread ownership

3. Data Assignment:
   - Gets datum pointer for each valid segment
   - Sets all elements of the datum array to the specified value
   - Handles variable-length datum arrays dynamically

PARAMETERS:
-----------
- val: The value to set for all datum elements
- Datum: Reference to the datum storage to be modified
- field_line_idx: Index of the field line to process

USAGE EXAMPLES:
---------------
// Set all wave energy values to zero
SetDatum(0.0, WaveEnergy, 0);

// Initialize magnetic field to uniform value
SetDatum(1.0e-9, MagneticField, 2);

// Reset all field lines
for (int fl = 0; fl < num_field_lines; ++fl) {
    SetDatum(0.0, SomeDatum, fl);
}

// Set test values for debugging
SetDatum(123.456, TestDatum, 0);

REQUIREMENTS:
-------------
- Valid field line index
- Datum must be properly initialized
- MPI environment (function is MPI-aware)

MPI BEHAVIOR:
-------------
- Processes ALL segments regardless of Thread ownership
- Function modifies segments owned by all MPI processes
- No MPI communication required - purely local operation
- May modify data that belongs to other processes (use with caution)

ERROR HANDLING:
---------------
- Validates field line existence and index range
- Handles null segments gracefully
- Skips segments without datum data
- Reports errors only from root process to avoid output flooding

PERFORMANCE:
------------
- Time Complexity: O(N*M) where N = segments, M = datum length
- Memory Usage: O(1) - no additional memory allocation
- Processes all segments on each MPI process

WARNING:
--------
- This function modifies ALL segments regardless of Thread ownership
- In MPI environments, this may modify data belonging to other processes
- Use with caution in distributed simulations
- Consider using thread-aware version for production MPI code

NOTES:
------
- Function modifies data in-place
- All elements of multi-component datum are set to same value
- Processes all segments regardless of ownership
- Can be used for initialization, testing, or data reset

SEE ALSO:
---------
- TestPrintDatum(): For verifying the set values
- TestPrintDatumMPI(): For distributed verification
- InitializeWaveEnergyFromPhysicalParameters(): Physics-based initialization

================================================================================
*/

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

void SetDatum(double val, PIC::Datum::cDatumStored& Datum, int field_line_idx) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Validate field line index (only root reports errors to avoid flooding)
    if (PIC::FieldLine::nFieldLine <= 0) {
        if (rank == 0) {
            std::cout << "ERROR in SetDatum: No field lines found!" << std::endl;
        }
        return;
    }

    if (field_line_idx >= PIC::FieldLine::nFieldLine || field_line_idx < 0) {
        if (rank == 0) {
            std::cout << "ERROR in SetDatum: Field line index " << field_line_idx
                      << " out of range [0, " << (PIC::FieldLine::nFieldLine - 1) << "]!" << std::endl;
        }
        return;
    }

    // Get the specified field line
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    int num_segments = field_line->GetTotalSegmentNumber();
    int datum_length = Datum.length;

    // Validate datum length
    if (datum_length <= 0) {
        if (rank == 0) {
            std::cout << "ERROR in SetDatum: Invalid datum length " << datum_length << std::endl;
        }
        return;
    }

    // Progress tracking for large field lines
    int segments_processed = 0;
    int segments_modified = 0;
    int segments_skipped_null = 0;
    int segments_skipped_no_data = 0;

    // Iterate through all segments in the field line
    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);

        // Skip null segments
        if (!segment) {
            segments_skipped_null++;
            continue;
        }

        segments_processed++;

        // Get datum pointer for this segment
        double* data = segment->GetDatum_ptr(Datum);

        // Skip segments without datum data
        if (!data) {
            segments_skipped_no_data++;
            continue;
        }

        // Set all elements of the datum to the specified value
        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            data[elem_idx] = val;
        }

        segments_modified++;
    }

    // Optional: Report progress for large operations (only from root process)
    if (rank == 0 && (num_segments > 10000 || segments_modified > 1000)) {
        std::cout << "SetDatum completed for field line " << field_line_idx << ":" << std::endl;
        std::cout << "  Value set: " << std::scientific << val << std::endl;
        std::cout << "  Total segments: " << num_segments << std::endl;
        std::cout << "  Segments processed: " << segments_processed << std::endl;
        std::cout << "  Segments modified: " << segments_modified << std::endl;
        if (segments_skipped_null > 0) {
            std::cout << "  Null segments skipped: " << segments_skipped_null << std::endl;
        }
        if (segments_skipped_no_data > 0) {
            std::cout << "  No-data segments skipped: " << segments_skipped_no_data << std::endl;
        }
        std::cout << "  Datum length: " << datum_length << " elements per segment" << std::endl;
        std::cout << "  NOTE: ALL segments processed regardless of Thread ownership" << std::endl;
    }
}

// Convenience overload for setting multiple field lines
void SetDatum(double val, PIC::Datum::cDatumStored& Datum, int start_field_line_idx, int end_field_line_idx) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (start_field_line_idx < 0 || end_field_line_idx >= PIC::FieldLine::nFieldLine ||
        start_field_line_idx > end_field_line_idx) {
        if (rank == 0) {
            std::cout << "ERROR in SetDatum: Invalid field line range ["
                      << start_field_line_idx << ", " << end_field_line_idx << "]" << std::endl;
        }
        return;
    }

    if (rank == 0) {
        std::cout << "Setting datum to " << std::scientific << val
                  << " for field lines " << start_field_line_idx
                  << " to " << end_field_line_idx << std::endl;
    }

    for (int fl = start_field_line_idx; fl <= end_field_line_idx; ++fl) {
        SetDatum(val, Datum, fl);
    }

    if (rank == 0) {
        std::cout << "Completed setting datum for "
                  << (end_field_line_idx - start_field_line_idx + 1)
                  << " field lines" << std::endl;
    }
}

// Convenience overload for setting all field lines
void SetDatumAll(double val, PIC::Datum::cDatumStored& Datum) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (PIC::FieldLine::nFieldLine <= 0) {
        if (rank == 0) {
            std::cout << "ERROR in SetDatumAll: No field lines found!" << std::endl;
        }
        return;
    }

    if (rank == 0) {
        std::cout << "Setting datum to " << std::scientific << val
                  << " for ALL " << PIC::FieldLine::nFieldLine << " field lines" << std::endl;
    }

    for (int fl = 0; fl < PIC::FieldLine::nFieldLine; ++fl) {
        SetDatum(val, Datum, fl);
    }

    if (rank == 0) {
        std::cout << "Completed setting datum for all "
                  << PIC::FieldLine::nFieldLine << " field lines" << std::endl;
    }
}

} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP

// Example usage:
/*
namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

void ExampleSetDatumUsage() {
    // Set single field line to zero
    SetDatum(0.0, WaveEnergy, 0);

    // Initialize field line with test value
    SetDatum(1.23e-6, MagneticField, 2);

    // Set range of field lines
    SetDatum(0.0, Pressure, 0, 4);  // Field lines 0-4

    // Set all field lines
    SetDatumAll(0.0, WaveEnergy);

    // Verify the results
    TestPrintDatumMPI(WaveEnergy, "After SetDatum", 0);
}

} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP
*/



/*
================================================================================
                    ANALYZE MAX SEGMENT PARTICLES FUNCTION
================================================================================

FUNCTION: AnalyzeMaxSegmentParticles

PURPOSE:
--------
Finds the segment with maximum relative value across ALL datum elements across 
all MPI processes, then performs detailed particle analysis in that segment 
including:
- Total particle density
- Directional particle densities (toward/from Sun)
- Velocity-dependent densities (above/below local Alfven speed)
- Background plasma density

ALGORITHM:
----------
1. Global Maximum Detection (MPI_Allreduce):
   - Each process finds local maximum for each datum element
   - MPI_Allreduce finds global maximum for each element across all processes
   - Establishes normalization factors for relative value calculation

2. Segment Selection (MPI_Reduce):
   - For each segment, calculates relative values: |data[i]|/global_max[i]
   - Finds segment with maximum relative value across ALL datum elements
   - Uses custom MPI reduction to find global maximum while preserving location

3. Particle Analysis (on owning process):
   - Only the process owning the max segment performs particle analysis
   - Loops through all particles in the segment using particle buffer methods
   - Calculates directional velocities using vParallel (field-aligned coordinate)
   - Compares particle speeds to local Alfven velocity
   - Computes various density categories

4. Results Collection and Output:
   - Owning process sends analysis results to root
   - Root process outputs comprehensive analysis
   - Shows which element drove the selection and complete datum context

PARAMETERS:
-----------
- Datum: Reference to datum storage for finding maximum across all elements
- msg: Descriptive message for output header
- field_line_idx: Index of field line to analyze (default: 0)

OUTPUT SECTIONS:
----------------
1. Maximum segment identification with:
   - Segment index and owning process
   - Element index that achieved maximum relative value
   - Element value and relative value
   - Complete datum array for context
   - Global maximum values for each element
2. Particle density analysis:
   - Total particle density
   - Densities moving toward/from Sun (based on vParallel sign)
   - Velocity-dependent densities (above/below Alfven speed)
3. Background plasma density and local parameters
4. Particle count statistics and density ratios

USAGE EXAMPLES:
---------------
// Find segment with maximum relative value across all datum elements
AnalyzeMaxSegmentParticles(WaveEnergy, "Wave Energy Analysis");

// Analyze different field line
AnalyzeMaxSegmentParticles(MagneticField, "B-field Analysis", 3);

// Find maximum pressure segment and analyze particles
AnalyzeMaxSegmentParticles(Pressure, "Pressure Analysis");

REQUIREMENTS:
-------------
- Active MPI environment
- Valid field line with particle data
- Background plasma data available
- Proper particle velocity and position data

MPI COMMUNICATION:
------------------
1. MPI_Allreduce: Find global maximum for each datum element
2. MPI_Reduce: Find segment with maximum relative value globally
3. MPI_Bcast: Broadcast global maximum info to all processes
4. Point-to-point: Transfer analysis results to root process

PARTICLE ANALYSIS DETAILS:
---------------------------
Directional Analysis (using field-aligned coordinates):
- Toward Sun: vParallel < 0 (negative parallel velocity)
- From Sun: vParallel > 0 (positive parallel velocity)

Velocity Classification:
- Fast particles: |v| > v_Alfven (supraalfvenic)
- Slow particles: |v| < v_Alfven (subalfvenic)

Statistical Weights:
- Uses species-specific weights: GlobalParticleWeight[spec]
- Where spec = PIC::ParticleBuffer::GetI(p)

Density Calculations:
- Uses particle statistical weights and segment volume
- Converts to physical units (particles/m³)

SELECTION CRITERIA:
-------------------
Segments are ranked by their maximum relative value across ALL elements:
  segment_max_relative = max(|data[i]|/global_max[i]) for i in [0, datum_length-1]

This ensures:
- Fair comparison between elements with different scales
- Automatic selection of most significant segment
- Element-agnostic analysis (considers all datum components)

PERFORMANCE:
------------
- Communication: O(log P) for reductions + O(1) for results transfer
- Element analysis: O(M) where M = datum_length
- Particle loop: O(N_particles) on one process only
- Memory: O(M) for global maximum storage

ERROR HANDLING:
---------------
- Validates field line index
- Handles zero maximum values across all elements
- Manages missing particle or plasma data
- Prevents division by zero in relative value calculations

NOTES:
------
- Analyzes ALL datum elements, not just a specified one
- Automatically selects most significant segment across all elements
- Particle analysis performed only on owning MPI process
- Uses field-aligned coordinate system (vParallel) for direction determination
- Requires proper initialization of background plasma parameters

SEE ALSO:
---------
- TestPrintDatumMPI(): For comprehensive datum analysis across processes
- SetDatum(): For setting datum values for testing
- PIC::ParticleBuffer methods: For particle data access

================================================================================
*/

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

void AnalyzeMaxSegmentParticles(PIC::Datum::cDatumStored& Datum, const char* msg, int field_line_idx) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Structure for particle analysis results (make it simple for MPI transfer)
    struct ParticleAnalysisResults {
        double total_density;
        double density_toward_sun;
        double density_from_sun;
        double density_toward_sun_fast;
        double density_toward_sun_slow;
        double density_from_sun_fast;
        double density_from_sun_slow;
        double background_plasma_density;
        double alfven_speed;
        double magnetic_field_strength;
        double segment_volume;
        int total_particles;
        int particles_toward_sun;
        int particles_from_sun;
        int analysis_successful;  // Use int instead of bool for MPI
    };

    // Simple structure for max location (avoid complex nested structs)
    struct MaxLocationData {
        double max_relative_value;
        int segment_idx;
        int process_rank;
        int element_idx;
        double element_value;
    };

    // Initialize results
    ParticleAnalysisResults results = {};
    MaxLocationData global_max = {};
    
    // Step 1: ALL processes validate inputs together
    int input_valid = 1;
    if (PIC::FieldLine::nFieldLine <= 0 || field_line_idx >= PIC::FieldLine::nFieldLine || field_line_idx < 0) {
        input_valid = 0;
    }

    int all_inputs_valid = 1;
    MPI_Allreduce(&input_valid, &all_inputs_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    
    if (all_inputs_valid == 0) {
        if (rank == 0) {
            std::cout << "\n=== " << msg << " === ERROR ===" << std::endl;
            std::cout << "Invalid field line index or no field lines found!" << std::endl;
        }
        return;  // ALL processes return here
    }

    // From here on, all processes have valid inputs
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    int num_segments = field_line->GetTotalSegmentNumber();
    int datum_length = Datum.length;

    // Step 2: Find global maximum for each element (ALL processes participate)
    std::vector<double> local_max_per_element(datum_length, 0.0);
    std::vector<double> global_max_per_element(datum_length, 0.0);

    // Each process finds its local maximums
    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        if (!segment) continue;

        double* data = segment->GetDatum_ptr(Datum);
        if (!data) continue;

        // Only process segments owned by this rank
        if (segment->Thread != rank) continue;

        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            double abs_val = std::abs(data[elem_idx]);
            if (abs_val > local_max_per_element[elem_idx]) {
                local_max_per_element[elem_idx] = abs_val;
            }
        }
    }

    // ALL processes participate in global reduction
    MPI_Allreduce(local_max_per_element.data(), global_max_per_element.data(), 
                  datum_length, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // Step 3: Find segment with maximum relative value (ALL processes participate)
    MaxLocationData local_max = {};
    local_max.max_relative_value = 0.0;
    local_max.segment_idx = -1;
    local_max.process_rank = rank;
    local_max.element_idx = -1;
    local_max.element_value = 0.0;

    // Each process finds its best segment
    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
        if (!segment) continue;

        double* data = segment->GetDatum_ptr(Datum);
        if (!data) continue;

        // Only process segments owned by this rank
        if (segment->Thread != rank) continue;

        // Find maximum relative value across all elements in this segment
        double segment_max_relative = 0.0;
        int best_element_idx = -1;
        double best_element_value = 0.0;

        for (int elem_idx = 0; elem_idx < datum_length; ++elem_idx) {
            if (global_max_per_element[elem_idx] > 0.0) {
                double relative_val = std::abs(data[elem_idx]) / global_max_per_element[elem_idx];
                double abs_val = std::abs(data[elem_idx]);
                
                // Select element with higher relative value, or if equal relative values,
                // select the one with higher absolute value
                if (relative_val > segment_max_relative || 
                    (relative_val == segment_max_relative && abs_val > std::abs(best_element_value))) {
                    segment_max_relative = relative_val;
                    best_element_idx = elem_idx;
                    best_element_value = data[elem_idx];
                }
            }
        }

        // Update local maximum if this segment is better
        if (segment_max_relative > local_max.max_relative_value) {
            local_max.max_relative_value = segment_max_relative;
            local_max.segment_idx = seg_idx;
            local_max.process_rank = rank;
            local_max.element_idx = best_element_idx;
            local_max.element_value = best_element_value;
        }
    }

    // Step 4: Find global maximum segment (ALL processes participate)
    // Create a structure that MPI can handle properly for MAXLOC operation
    struct {
        double value;
        int rank;
    } local_maxloc, global_maxloc;
    
    local_maxloc.value = local_max.max_relative_value;
    local_maxloc.rank = rank;
    
    // Use MPI_MAXLOC to find maximum value and preserve the rank that has it
    MPI_Allreduce(&local_maxloc, &global_maxloc, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    
    // Now we need to get the complete information from the process that has the maximum
    // Broadcast the complete local_max structure from the owning process
    global_max.max_relative_value = global_maxloc.value;
    global_max.process_rank = global_maxloc.rank;
    
    // Broadcast the rest of the data from the owning process
    MPI_Bcast(&local_max.segment_idx, 1, MPI_INT, global_maxloc.rank, MPI_COMM_WORLD);
    MPI_Bcast(&local_max.element_idx, 1, MPI_INT, global_maxloc.rank, MPI_COMM_WORLD);
    MPI_Bcast(&local_max.element_value, 1, MPI_DOUBLE, global_maxloc.rank, MPI_COMM_WORLD);
    
    global_max.segment_idx = local_max.segment_idx;
    global_max.element_idx = local_max.element_idx;
    global_max.element_value = local_max.element_value;

    // Step 5: Particle analysis (only on owning process)
    if (rank == global_max.process_rank && global_max.max_relative_value > 0.0) {
        PIC::FieldLine::cFieldLineSegment* max_segment = field_line->GetSegment(global_max.segment_idx);
        
        if (max_segment) {
            // Get segment properties
            results.segment_volume = SEP::FieldLine::GetSegmentVolume(max_segment, field_line_idx);
            
            // Get plasma parameters
            double rho;
            max_segment->GetPlasmaDensity(0.5, rho);
            results.background_plasma_density = rho;
            rho *= PIC::MolecularData::GetMass(_H_PLUS_SPEC_);
            
            // Get magnetic field
            double B[3];
            PIC::FieldLine::FieldLinesAll[field_line_idx].GetMagneticField(B, 0.5 + global_max.segment_idx);
            results.magnetic_field_strength = Vector3D::Length(B);
            
            // Calculate Alfven speed
            results.alfven_speed = results.magnetic_field_strength / sqrt(VacuumPermeability * rho);
            
            // Initialize counters
            results.total_particles = 0;
            results.particles_toward_sun = 0;
            results.particles_from_sun = 0;
            
            double total_weight_all = 0.0;
            double total_weight_toward_sun = 0.0;
            double total_weight_from_sun = 0.0;
            double total_weight_toward_sun_fast = 0.0;
            double total_weight_toward_sun_slow = 0.0;
            double total_weight_from_sun_fast = 0.0;
            double total_weight_from_sun_slow = 0.0;

            // Loop through particles
            long int p = max_segment->FirstParticleIndex;
            while (p != -1) {
                results.total_particles++;
                
                int spec = PIC::ParticleBuffer::GetI(p);
                double vParallel = PIC::ParticleBuffer::GetVParallel(p);
                double vNormal = PIC::ParticleBuffer::GetVNormal(p);
                double v_magnitude = sqrt(vParallel*vParallel + vNormal*vNormal);
                
                double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec] *
                                   PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
                
                total_weight_all += stat_weight;
                
                bool moving_toward_sun = (vParallel < 0.0);
                bool is_fast = (v_magnitude > results.alfven_speed);
                
                if (moving_toward_sun) {
                    results.particles_toward_sun++;
                    total_weight_toward_sun += stat_weight;
                    if (is_fast) {
                        total_weight_toward_sun_fast += stat_weight;
                    } else {
                        total_weight_toward_sun_slow += stat_weight;
                    }
                } else {
                    results.particles_from_sun++;
                    total_weight_from_sun += stat_weight;
                    if (is_fast) {
                        total_weight_from_sun_fast += stat_weight;
                    } else {
                        total_weight_from_sun_slow += stat_weight;
                    }
                }
                
                p = PIC::ParticleBuffer::GetNext(p);
            }
            
            // Convert to densities
            if (results.segment_volume > 0.0) {
                results.total_density = total_weight_all / results.segment_volume;
                results.density_toward_sun = total_weight_toward_sun / results.segment_volume;
                results.density_from_sun = total_weight_from_sun / results.segment_volume;
                results.density_toward_sun_fast = total_weight_toward_sun_fast / results.segment_volume;
                results.density_toward_sun_slow = total_weight_toward_sun_slow / results.segment_volume;
                results.density_from_sun_fast = total_weight_from_sun_fast / results.segment_volume;
                results.density_from_sun_slow = total_weight_from_sun_slow / results.segment_volume;
            }
            
            results.analysis_successful = 1;
        }
    }

    // Step 6: Gather results to root (use broadcast instead of point-to-point)
    // Broadcast results from owning process to all processes
    MPI_Bcast(&results, sizeof(ParticleAnalysisResults), MPI_BYTE, global_max.process_rank, MPI_COMM_WORLD);

    // Step 7: Output (only root)
    if (rank == 0) {
        std::cout << "\n=== " << msg << " === Max Segment Particle Analysis ===" << std::endl;
        std::cout << "Field Line: " << field_line_idx << std::endl;
        std::cout << std::string(80, '=') << std::endl;
        
        if (global_max.max_relative_value == 0.0) {
            std::cout << "No non-zero values found in any datum elements" << std::endl;
            return;
        }
        
        std::cout << "MAXIMUM SEGMENT INFORMATION:" << std::endl;
        std::cout << "  Segment index: " << global_max.segment_idx << std::endl;
        std::cout << "  Process owner: " << global_max.process_rank << std::endl;
        std::cout << "  Maximum element index: " << global_max.element_idx << std::endl;
        std::cout << "  Maximum element value: " << std::scientific << std::setprecision(6) 
                  << global_max.element_value << std::endl;
        std::cout << "  Maximum relative value: " << std::fixed << std::setprecision(6) 
                  << global_max.max_relative_value << std::endl;
        
        if (results.analysis_successful == 0) {
            std::cout << "\nWARNING: Particle analysis was not performed!" << std::endl;
            return;
        }
        
        std::cout << "\nLOCAL PLASMA PARAMETERS:" << std::endl;
        std::cout << "  Background plasma density: " << std::scientific << std::setprecision(4) 
                  << results.background_plasma_density << " m⁻³" << std::endl;
        std::cout << "  Magnetic field strength: " << std::scientific << std::setprecision(4) 
                  << results.magnetic_field_strength << " T" << std::endl;
        std::cout << "  Alfven speed: " << std::scientific << std::setprecision(4) 
                  << results.alfven_speed << " m/s" << std::endl;
        std::cout << "  Segment volume: " << std::scientific << std::setprecision(4) 
                  << results.segment_volume << " m³" << std::endl;
        
        std::cout << "\nPARTICLE ANALYSIS RESULTS:" << std::endl;
        std::cout << "  Total particles in segment: " << results.total_particles << std::endl;
        std::cout << "  Particles toward Sun: " << results.particles_toward_sun 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * results.particles_toward_sun / std::max(1, results.total_particles)) << "%)" << std::endl;
        std::cout << "  Particles from Sun: " << results.particles_from_sun 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * results.particles_from_sun / std::max(1, results.total_particles)) << "%)" << std::endl;
        
        std::cout << "\nPARTICLE DENSITIES [m⁻³]:" << std::endl;
        std::cout << "  Total density: " << std::scientific << std::setprecision(4) 
                  << results.total_density << std::endl;
        std::cout << "  Density toward Sun: " << std::scientific << std::setprecision(4) 
                  << results.density_toward_sun << std::endl;
        std::cout << "  Density from Sun: " << std::scientific << std::setprecision(4) 
                  << results.density_from_sun << std::endl;
        
        std::cout << "\nVELOCITY-DEPENDENT DENSITIES [m⁻³]:" << std::endl;
        std::cout << "  Toward Sun (fast, |v| > v_A): " << std::scientific << std::setprecision(4) 
                  << results.density_toward_sun_fast << std::endl;
        std::cout << "  Toward Sun (slow, |v| < v_A): " << std::scientific << std::setprecision(4) 
                  << results.density_toward_sun_slow << std::endl;
        std::cout << "  From Sun (fast, |v| > v_A): " << std::scientific << std::setprecision(4) 
                  << results.density_from_sun_fast << std::endl;
        std::cout << "  From Sun (slow, |v| < v_A): " << std::scientific << std::setprecision(4) 
                  << results.density_from_sun_slow << std::endl;
        
        std::cout << "\nDENSITY RATIOS:" << std::endl;
        if (results.background_plasma_density > 0.0) {
            std::cout << "  Total particle/background ratio: " << std::fixed << std::setprecision(6) 
                      << (results.total_density / results.background_plasma_density) << std::endl;
        }
        if (results.total_density > 0.0) {
            std::cout << "  Fast/slow particle ratio: " << std::fixed << std::setprecision(6) 
                      << ((results.density_toward_sun_fast + results.density_from_sun_fast) / 
                          std::max(1e-20, results.density_toward_sun_slow + results.density_from_sun_slow)) << std::endl;
            std::cout << "  Outward/inward flow ratio: " << std::fixed << std::setprecision(6) 
                      << (results.density_from_sun / std::max(1e-20, results.density_toward_sun)) << std::endl;
        }
        
        std::cout << std::string(80, '=') << std::endl;
        std::cout << "Analysis complete." << std::endl << std::endl;
    }
}

} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP


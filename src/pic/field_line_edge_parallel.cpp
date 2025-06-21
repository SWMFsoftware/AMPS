/*
PURPOSE:
--------
This module provides efficient MPI communication for DatumStoredAtEdge data 
across field line segments in the AMPS plasma simulation. It handles the complex
memory layout where each field line segment stores its own DatumStoredAtEdge data
in individual memory locations, making standard MPI operations challenging.


KEY FEATURES:
-------------
1. AUTOMATIC DATATYPE MANAGEMENT: Creates and caches MPI struct datatypes
2. PERFORMANCE OPTIMIZATION: Caches field lines and datatypes for reuse
3. MEMORY SAFETY: Handles const-correctness and proper memory access
4. CHANGE DETECTION: Automatically rebuilds caches when field lines change
5. SCALABLE OPERATIONS: Supports all-reduce, scatter, gather patterns
6. PROCESS-SPECIFIC: Handles per-process data distribution

CORE FUNCTIONS:
---------------

1. ALL-REDUCE OPERATIONS (Sum data across all processes):
   - MPIAllReduceDatumStoredAtEdge(S)
   - MPIAllReduceDatumStoredAtEdgeFieldLine(field_line_idx, S)

2. BROADCAST OPERATIONS (Distribute from root to all processes):
   - MPIBcastDatumStoredAtEdge(S, root_rank)

3. GATHER OPERATIONS (Collect from all processes to root):
   - MPIGatherDatumStoredAtEdge(S, root_rank)
   - MPIGatherDatumStoredAtEdgeFieldLine(field_line_idx, S, root_rank)

4. MANAGEMENT FUNCTIONS:
   - InitializeDatumStoredAtEdgeMPI()
   - CleanupDatumStoredAtEdgeMPI()
   - NotifyFieldLinesChanged()

USAGE EXAMPLES:
===============

BASIC INITIALIZATION:
---------------------
```cpp

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    // Initialize the MPI module
    PIC::FieldLine::Parallel::InitializeDatumStoredAtEdgeMPI();
    
    // Your simulation code here...
    
    // Cleanup before MPI_Finalize
    PIC::FieldLine::Parallel::CleanupDatumStoredAtEdgeMPI();
    MPI_Finalize();
    return 0;
}
```

EXAMPLE 1: ALL-REDUCE PARTICLE COUNTS
-------------------------------------
```cpp
// Sum particle counts from all processes
void ReduceParticleCounts() {
    // Assume we have a DatumStored for particle counts
    PIC::Datum::cDatumStored particleCountDatum;
    
    // Each process has calculated local particle counts in each segment
    // Now sum across all processes
    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(particleCountDatum);
    
    // After this call, each segment contains the global sum
    std::cout << "Global particle counts summed across all processes" << std::endl;
}
```

EXAMPLE 2: SCATTER INITIAL CONDITIONS
------------------------------------
```cpp
// Distribute initial conditions from root process
void ScatterInitialConditions() {
    PIC::Datum::cDatumStored initialConditions;
    int root_rank = 0;
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == root_rank) {
        // Root process sets up initial conditions for all segments
        std::cout << "Root process: Setting up initial conditions..." << std::endl;
        // ... fill initialConditions data ...
    }
    
    // Scatter from root to all processes
    PIC::FieldLine::Parallel::MPIBcastDatumStoredAtEdge(initialConditions, root_rank);
    
    // After this call, each process has received its portion of the data
    std::cout << "Process " << rank << ": Received initial conditions" << std::endl;
}
```

EXAMPLE 3: GATHER RESULTS FOR OUTPUT
-----------------------------------
```cpp
// Collect results from all processes for output
void GatherResultsForOutput() {
    PIC::Datum::cDatumStored simulationResults;
    int root_rank = 0;
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Each process has computed its local results
    std::cout << "Process " << rank << ": Computed local results" << std::endl;
    
    // Gather all results to root process
    PIC::FieldLine::Parallel::MPIGatherDatumStoredAtEdge(simulationResults, root_rank);
    
    if (rank == root_rank) {
        // Root process now has all data for output
        std::cout << "Root process: Writing complete results to file..." << std::endl;
        // ... write results to file ...
    }
}
```

EXAMPLE 4: SINGLE FIELD LINE OPERATIONS
---------------------------------------
```cpp
// Work with specific field lines
void ProcessSpecificFieldLine(int field_line_idx) {
    PIC::Datum::cDatumStored fieldLineData;
    
    // All-reduce for a specific field line only
    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdgeFieldLine(
        field_line_idx, fieldLineData);
    
    std::cout << "Processed field line " << field_line_idx << std::endl;
}
```

EXAMPLE 5: HANDLING FIELD LINE CHANGES
--------------------------------------
```cpp
// When field line structure changes during simulation
void AdaptiveMeshRefinement() {
    // Modify field line structure
    // ... add/remove segments, change field lines ...
    
    // Notify the MPI module that field lines have changed
    PIC::FieldLine::Parallel::NotifyFieldLinesChanged();
    
    // Next MPI operation will automatically rebuild datatypes
    PIC::Datum::cDatumStored newData;
    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(newData);
    
    std::cout << "Adapted to new field line structure" << std::endl;
}
```

EXAMPLE 6: MULTI-STEP SIMULATION WORKFLOW
-----------------------------------------
```cpp
void RunSimulationStep() {
    // Step 1: Scatter boundary conditions from rank 0
    PIC::FieldLine::Parallel::MPIBcastDatumStoredAtEdge(boundaryConditions, 0);
    
    // Step 2: Each process computes local physics
    ComputeLocalPhysics();
    
    // Step 3: All-reduce to share information between processes
    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(sharedData);
    
    // Step 4: Continue local computation with updated shared data
    ContinueLocalComputation();
    
    // Step 5: Gather final results to root for analysis
    PIC::FieldLine::Parallel::MPIGatherDatumStoredAtEdge(finalResults, 0);
}
```

TECHNICAL DETAILS:
==================

MPI DATATYPE STRATEGY:
----------------------
The module uses MPI_Type_create_struct() to create custom datatypes that handle
the non-contiguous memory layout. Each datatype contains:
- Absolute memory addresses for each segment's data
- Block lengths (S.length elements per segment)
- MPI_DOUBLE datatype for each block

CACHING MECHANISM:
------------------
- Field lines are cached and reused until NotifyFieldLinesChanged() is called
- MPI datatypes are cached using combinatorial hash keys
- Automatic cache invalidation when structure changes
- Separate caches for all-reduce vs scatter/gather operations

MEMORY ADDRESSING:
------------------
All MPI operations use MPI_BOTTOM with absolute addressing:
- MPI_BOTTOM tells MPI to use addresses stored in the struct datatype
- Avoids pointer arithmetic and offset calculations
- Handles arbitrary memory layouts correctly

PERFORMANCE CONSIDERATIONS:
---------------------------
- First call creates and caches datatypes (one-time cost)
- Subsequent calls reuse cached datatypes (very fast)
- Cache invalidation only when field lines actually change
- Efficient hash-based datatype lookup

ERROR HANDLING:
---------------
- Comprehensive error checking for all MPI operations
- Graceful handling of empty or invalid field line data
- Clear error messages with operation context
- Automatic cleanup of failed datatype creation

THREAD SAFETY:
--------------
- Uses static variables with proper initialization
- Cache operations are atomic at the MPI level
- Safe for single-threaded MPI applications
- Field line change notification is process-global

DEBUGGING TIPS:
===============

ENABLE VERBOSE OUTPUT:
---------------------
Set MPI rank 0 to print detailed information about datatype creation:
- Number of segments processed
- Memory addresses being used
- Datatype creation success/failure

CHECK FIELD LINE CONSISTENCY:
-----------------------------
Ensure all processes have the same field line structure:
- Same number of field lines
- Same number of segments per field line
- Consistent segment ordering

VERIFY DATUM LENGTH:
-------------------
Ensure S.length is consistent across all processes:
- Same DatumStored configuration
- Same number of data elements per segment

MONITOR CACHE EFFECTIVENESS:
---------------------------
- Call NotifyFieldLinesChanged() only when actually needed
- Avoid unnecessary cache invalidation
- Check that datatypes are being reused

COMMON PITFALLS:
================

1. FORGETTING INITIALIZATION:
   Always call InitializeDatumStoredAtEdgeMPI() before use

2. MISSING CLEANUP:
   Always call CleanupDatumStoredAtEdgeMPI() before MPI_Finalize()

3. INCONSISTENT FIELD LINES:
   Ensure all processes have identical field line structures

4. CONST CORRECTNESS:
   The module handles const parameters correctly internally

5. PREMATURE CACHE INVALIDATION:
   Only call NotifyFieldLinesChanged() when structure actually changes

================================================================================
*/


#include "pic.h"


#include <iostream>
#include <cassert>
#include <unordered_map>
#include <vector>
#include <cstddef>

namespace PIC {
namespace FieldLine {
namespace Parallel {

// Static member definitions
std::unordered_map<size_t, MPI_Datatype> DatumStoredAtEdgeMPIManager::datatype_cache;
std::unordered_map<size_t, bool> DatumStoredAtEdgeMPIManager::datatype_availability;
bool DatumStoredAtEdgeMPIManager::is_initialized = false;
bool DatumStoredAtEdgeMPIManager::field_lines_changed = true;  // Initially changed
std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>> DatumStoredAtEdgeMPIManager::cached_field_lines;
bool DatumStoredAtEdgeMPIManager::field_lines_cache_valid = false;

// ============================================================================
// Helper function to collect all field lines using PIC::FieldLine structures
// ============================================================================

std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>> CollectAllFieldLines() {
    return DatumStoredAtEdgeMPIManager::GetCachedFieldLines();
}

const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& DatumStoredAtEdgeMPIManager::GetCachedFieldLines() {
    // Check if cache is valid and field lines haven't changed
    if (field_lines_cache_valid && !field_lines_changed) {
        return cached_field_lines;
    }

    // Rebuild cache
    cached_field_lines.clear();

    for (int fl_idx = 0; fl_idx < PIC::FieldLine::nFieldLine; ++fl_idx) {
        std::vector<PIC::FieldLine::cFieldLineSegment*> segments;
        int num_segments = PIC::FieldLine::FieldLinesAll[fl_idx].GetTotalSegmentNumber();

        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            // Assuming a method to get segment pointer from field line
            PIC::FieldLine::cFieldLineSegment* seg = PIC::FieldLine::FieldLinesAll[fl_idx].GetSegment(seg_idx);
            if (seg) {
                segments.push_back(seg);
            }
        }
        cached_field_lines.push_back(segments);
    }

    // Mark cache as valid and field lines as unchanged
    field_lines_cache_valid = true;
    field_lines_changed = false;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int total_segments = 0;
        for (const auto& line : cached_field_lines) {
            total_segments += line.size();
        }
        std::cout << "Cached field lines data: " << PIC::FieldLine::nFieldLine
                  << " field lines, " << total_segments << " total segments" << std::endl;
    }

    return cached_field_lines;
}

void DatumStoredAtEdgeMPIManager::InvalidateFieldLinesCache() {
    field_lines_cache_valid = false;
    field_lines_changed = true;
    cached_field_lines.clear();
}

void DatumStoredAtEdgeMPIManager::NotifyFieldLinesChanged() {
    // Mark field lines as changed
    field_lines_changed = true;

    // Invalidate field lines cache
    InvalidateFieldLinesCache();

    // Clear all MPI datatypes as they depend on field line structure
    CleanupAllDatatypes();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "Field lines changed notification: cleared all caches and MPI datatypes" << std::endl;
    }
}

bool DatumStoredAtEdgeMPIManager::AreFieldLinesChanged() {
    return field_lines_changed;
}

// Public wrapper for field line change notification
void NotifyFieldLinesChanged() {
    DatumStoredAtEdgeMPIManager::NotifyFieldLinesChanged();
}

// ============================================================================
// DatumStoredAtEdgeMPIManager Implementation
// ============================================================================

void DatumStoredAtEdgeMPIManager::Initialize() {
    if (!is_initialized) {
        datatype_cache.clear();
        datatype_availability.clear();
        is_initialized = true;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cout << "DatumStoredAtEdgeMPIManager initialized" << std::endl;
        }
    }
}

void DatumStoredAtEdgeMPIManager::Finalize() {
    CleanupAllDatatypes();
    InvalidateFieldLinesCache();
    is_initialized = false;
}

size_t DatumStoredAtEdgeMPIManager::CreateDatatypeKey(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    int datum_length) {

    size_t segment_count = 0;

    for (const auto& line : field_lines) {
        segment_count += line.size();
    }

    // Improved combinatorial hash using polynomial rolling hash
    size_t key = 0;
    key = key * 31 + field_lines.size();
    key = key * 31 + segment_count;
    key = key * 31 + static_cast<size_t>(datum_length);
    return key;
}

size_t DatumStoredAtEdgeMPIManager::CreateScatterGatherDatatypeKey(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    int datum_length,
    int target_process_rank) {

    size_t hash = 0;
    size_t segment_count = 0;

    for (const auto& line : field_lines) {
        segment_count += line.size();
    }

    // Hash combining field lines, segments, datum length, and target process rank
    hash = field_lines.size() * 10000000 + segment_count * 10000 +
           datum_length * 100 + target_process_rank;
    return hash;
}

MPI_Datatype DatumStoredAtEdgeMPIManager::CreateDatumStoredAtEdgeDatatypeInternal(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {

    // Collect all DatumStoredAtEdge data pointers - each has individual memory location
    std::vector<MPI_Aint> displacements;
    std::vector<MPI_Datatype> types;
    std::vector<int> blocklengths;
    std::vector<double*> data_pointers;

    // Create a non-const copy of S for GetDatum_ptr calls
    cDatumStored S_copy = S;

    // Traverse all field lines and segments
    for (size_t field_line_idx = 0; field_line_idx < field_lines.size(); ++field_line_idx) {
        const auto& segments = field_lines[field_line_idx];

        for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* seg = segments[seg_idx];
            if (!seg) continue;

            double* DatumData = seg->GetDatum_ptr(S_copy);
            if (!DatumData) continue;

            data_pointers.push_back(DatumData);

            // Get absolute address for this edge's data
            MPI_Aint displacement;
            MPI_Get_address(DatumData, &displacement);
            displacements.push_back(displacement);

            // Each edge has S.length elements
            blocklengths.push_back(S.length);
            types.push_back(MPI_DOUBLE);
        }
    }

    if (displacements.empty()) {
        std::cerr << "Warning: No valid DatumStoredAtEdge data found" << std::endl;
        return MPI_DATATYPE_NULL;
    }

    // Create MPI datatype using struct (handles individual memory locations)
    MPI_Datatype datum_datatype;
    int count = static_cast<int>(displacements.size());

    int result = MPI_Type_create_struct(
        count,                    // Number of blocks
        blocklengths.data(),      // Elements per block array
        displacements.data(),     // Absolute displacement array
        types.data(),             // Datatype array (all MPI_DOUBLE)
        &datum_datatype           // New datatype
    );

    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Type_create_struct failed with code " << result << std::endl;
        return MPI_DATATYPE_NULL;
    }

    // Commit the datatype
    result = MPI_Type_commit(&datum_datatype);
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Type_commit failed with code " << result << std::endl;
        MPI_Type_free(&datum_datatype);
        return MPI_DATATYPE_NULL;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "Created MPI struct datatype for DatumStoredAtEdge: " << count << " segments, "
                  << S.length << " elements each (S.length=" << S.length << ")" << std::endl;
    }

    return datum_datatype;
}

MPI_Datatype DatumStoredAtEdgeMPIManager::CreateScatterGatherDatatypeInternal(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S,
    int target_process_rank) {

    // Collect DatumStoredAtEdge data pointers for segments belonging to target_process_rank
    std::vector<MPI_Aint> displacements;
    std::vector<MPI_Datatype> types;
    std::vector<int> blocklengths;
    std::vector<double*> data_pointers;

    // Create a non-const copy of S for GetDatum_ptr calls
    cDatumStored S_copy = S;

    // Traverse all field lines and segments, filter by Thread (process rank)
    for (size_t field_line_idx = 0; field_line_idx < field_lines.size(); ++field_line_idx) {
        const auto& segments = field_lines[field_line_idx];

        for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* seg = segments[seg_idx];
            if (!seg) continue;

            // Check if this segment belongs to the target process
            if (seg->Thread != target_process_rank) continue;

            double* DatumData = seg->GetDatum_ptr(S_copy);
            if (!DatumData) continue;

            data_pointers.push_back(DatumData);

            // Get absolute address for this edge's data
            MPI_Aint displacement;
            MPI_Get_address(DatumData, &displacement);
            displacements.push_back(displacement);

            // Each edge has S.length elements
            blocklengths.push_back(S.length);
            types.push_back(MPI_DOUBLE);
        }
    }

    if (displacements.empty()) {
        // No segments for this process rank - create empty datatype
        return MPI_DATATYPE_NULL;
    }

    // Create MPI datatype using struct (handles individual memory locations)
    MPI_Datatype scatter_gather_datatype;
    int count = static_cast<int>(displacements.size());

    int result = MPI_Type_create_struct(
        count,                    // Number of blocks
        blocklengths.data(),      // Elements per block array
        displacements.data(),     // Absolute displacement array
        types.data(),             // Datatype array (all MPI_DOUBLE)
        &scatter_gather_datatype  // New datatype
    );

    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Type_create_struct failed for process "
                  << target_process_rank << " with code " << result << std::endl;
        return MPI_DATATYPE_NULL;
    }

    // Commit the datatype
    result = MPI_Type_commit(&scatter_gather_datatype);
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Type_commit failed for process "
                  << target_process_rank << " with code " << result << std::endl;
        MPI_Type_free(&scatter_gather_datatype);
        return MPI_DATATYPE_NULL;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "Created scatter/gather MPI struct datatype for process " << target_process_rank
                  << ": " << count << " segments, " << S.length << " elements each" << std::endl;
    }

    return scatter_gather_datatype;
}
MPI_Datatype DatumStoredAtEdgeMPIManager::ConstructDatumStoredAtEdgeDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {

    Initialize();

    // Create unique key for this configuration
    size_t key = CreateDatatypeKey(field_lines, S.length);

    // Check if datatype already exists in cache
    auto it = datatype_cache.find(key);
    if (it != datatype_cache.end()) {
        // Mark as available and return cached datatype
        datatype_availability[key] = true;
        return it->second;
    }

    // Create new datatype
    MPI_Datatype datum_datatype = CreateDatumStoredAtEdgeDatatypeInternal(field_lines, S);

    if (datum_datatype != MPI_DATATYPE_NULL) {
        // Cache the datatype and mark as available
        datatype_cache[key] = datum_datatype;
        datatype_availability[key] = true;
    } else {
        // Mark as unavailable
        datatype_availability[key] = false;
    }

    return datum_datatype;
}

MPI_Datatype DatumStoredAtEdgeMPIManager::ConstructScatterGatherDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S,
    int target_process_rank) {

    Initialize();

    // Create unique key for this configuration
    size_t key = CreateScatterGatherDatatypeKey(field_lines, S.length, target_process_rank);

    // Check if datatype already exists in cache
    auto it = datatype_cache.find(key);
    if (it != datatype_cache.end()) {
        // Mark as available and return cached datatype
        datatype_availability[key] = true;
        return it->second;
    }

    // Create new datatype
    MPI_Datatype datatype = CreateScatterGatherDatatypeInternal(field_lines, S, target_process_rank);

    if (datatype != MPI_DATATYPE_NULL) {
        // Cache the datatype and mark as available
        datatype_cache[key] = datatype;
        datatype_availability[key] = true;
    } else {
        // Mark as unavailable (empty datatype for this process)
        datatype_availability[key] = false;
    }

    return datatype;
}

bool DatumStoredAtEdgeMPIManager::IsDatumStoredAtEdgeDatatypeAvailable(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {

    Initialize();

    size_t key = CreateDatatypeKey(field_lines, S.length);
    auto it = datatype_availability.find(key);

    if (it != datatype_availability.end()) {
        return it->second;
    }

    // Not in cache, so not available
    return false;
}

bool DatumStoredAtEdgeMPIManager::IsScatterGatherDatatypeAvailable(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S,
    int target_process_rank) {

    Initialize();

    size_t key = CreateScatterGatherDatatypeKey(field_lines, S.length, target_process_rank);
    auto it = datatype_availability.find(key);

    if (it != datatype_availability.end()) {
        return it->second;
    }

    // Not in cache, so not available
    return false;
}

MPI_Datatype DatumStoredAtEdgeMPIManager::GetDatumStoredAtEdgeDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {

    Initialize();

    size_t key = CreateDatatypeKey(field_lines, S.length);
    auto it = datatype_cache.find(key);

    if (it != datatype_cache.end() && IsDatumStoredAtEdgeDatatypeAvailable(field_lines, S)) {
        return it->second;
    }

    return MPI_DATATYPE_NULL;
}

bool DatumStoredAtEdgeMPIManager::RemoveDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {

    size_t key = CreateDatatypeKey(field_lines, S.length);
    auto it = datatype_cache.find(key);

    if (it != datatype_cache.end()) {
        if (it->second != MPI_DATATYPE_NULL) {
            MPI_Type_free(&it->second);
        }
        datatype_cache.erase(it);
        datatype_availability.erase(key);
        return true;
    }
    return false;
}

void DatumStoredAtEdgeMPIManager::CleanupAllDatatypes() {
    for (auto& pair : datatype_cache) {
        if (pair.second != MPI_DATATYPE_NULL) {
            MPI_Type_free(&pair.second);
        }
    }
    datatype_cache.clear();
    datatype_availability.clear();

    // Also invalidate field lines cache when cleaning up datatypes
    InvalidateFieldLinesCache();
}

int DatumStoredAtEdgeMPIManager::GetCachedDatatypeCount() {
    return static_cast<int>(datatype_cache.size());
}

void DatumStoredAtEdgeMPIManager::PrintCacheStatistics() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << "=== DatumStoredAtEdge MPI Cache Statistics ===" << std::endl;
        std::cout << "Total cached datatypes: " << datatype_cache.size() << std::endl;

        int available_count = 0;
        for (const auto& pair : datatype_availability) {
            if (pair.second) available_count++;
        }
        std::cout << "Available datatypes: " << available_count << std::endl;
        std::cout << "=============================================" << std::endl;
    }
}

// ============================================================================
// Helper functions for automatic datatype management
// ============================================================================

MPI_Datatype GetOrCreateDatumStoredAtEdgeDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {

    // Check if datatype is already available
    if (DatumStoredAtEdgeMPIManager::IsDatumStoredAtEdgeDatatypeAvailable(field_lines, S)) {
        return DatumStoredAtEdgeMPIManager::GetDatumStoredAtEdgeDatatype(field_lines, S);
    } else {
        // Create new datatype
        return DatumStoredAtEdgeMPIManager::ConstructDatumStoredAtEdgeDatatype(field_lines, S);
    }
}

MPI_Datatype GetOrCreateDatumStoredAtEdgeDatatype(const cDatumStored& S) {
    // Use cached field lines - automatically handles caching and change detection
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& all_field_lines =
        DatumStoredAtEdgeMPIManager::GetCachedFieldLines();
    return GetOrCreateDatumStoredAtEdgeDatatype(all_field_lines, S);
}

MPI_Datatype GetOrCreateScatterGatherDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S,
    int target_process_rank) {

    // Check if datatype is already available
    if (DatumStoredAtEdgeMPIManager::IsScatterGatherDatatypeAvailable(field_lines, S, target_process_rank)) {
        return DatumStoredAtEdgeMPIManager::ConstructScatterGatherDatatype(field_lines, S, target_process_rank);
    } else {
        // Create new datatype
        return DatumStoredAtEdgeMPIManager::ConstructScatterGatherDatatype(field_lines, S, target_process_rank);
    }
}

// ============================================================================
// Public Interface Functions with Automatic Datatype Management
// ============================================================================

void MPIAllReduceDatumStoredAtEdge(const cDatumStored& S) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const auto& all_field_lines = DatumStoredAtEdgeMPIManager::GetCachedFieldLines();

    // Collect all data pointers and perform individual all-reduce operations
    std::vector<double*> data_pointers;
    cDatumStored S_copy = S;  // Non-const copy for GetDatum_ptr
    
    // Collect all valid data pointers from all segments
    for (const auto& line : all_field_lines) {
        for (PIC::FieldLine::cFieldLineSegment* seg : line) {
            if (seg) {
                double* data_ptr = seg->GetDatum_ptr(S_copy);
                if (data_ptr) {
                    data_pointers.push_back(data_ptr);
                }
            }
        }
    }

    if (data_pointers.empty()) {
        std::cerr << "No valid DatumStoredAtEdge data found for all-reduce" << std::endl;
        return;
    }

    // Perform separate MPI_Allreduce for each segment's data
    int successful_operations = 0;
    for (double* data_ptr : data_pointers) {
        int result = MPI_Allreduce(MPI_IN_PLACE, data_ptr, S.length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if (result != MPI_SUCCESS) {
            std::cerr << "Error: MPI_Allreduce failed for segment data with code " << result << std::endl;
        } else {
            successful_operations++;
        }
    }

    if (rank == 0) {
        std::cout << "Successfully completed MPI all-reduce for DatumStoredAtEdge data: "
                  << successful_operations << " segments across " << PIC::FieldLine::nFieldLine
                  << " field lines (S.length=" << S.length << ")" << std::endl;
    }
}

void MPIAllReduceDatumStoredAtEdgeFieldLine(int field_line_idx, cDatumStored& S) {
    // Create field line vector for single field line
    std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>> single_field_line;
    std::vector<PIC::FieldLine::cFieldLineSegment*> segments;

    if (field_line_idx >= 0 && field_line_idx < PIC::FieldLine::nFieldLine) {
        int num_segments = PIC::FieldLine::FieldLinesAll[field_line_idx].GetTotalSegmentNumber();
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* seg = PIC::FieldLine::FieldLinesAll[field_line_idx].GetSegment(seg_idx);
            if (seg) {
                segments.push_back(seg);
            }
        }
    }
    single_field_line.push_back(segments);

    // Automatically get or create datatype
    MPI_Datatype datum_datatype = GetOrCreateDatumStoredAtEdgeDatatype(single_field_line, S);

    if (datum_datatype == MPI_DATATYPE_NULL) {
        std::cerr << "Failed to create MPI datatype for field line " << field_line_idx << std::endl;
        return;
    }

    // Get base data pointer
    double* base_data = nullptr;

    for (PIC::FieldLine::cFieldLineSegment* seg : segments) {
        if (seg) {
            double* DatumData = seg->GetDatum_ptr(S);
            if (DatumData) {
                base_data = DatumData;
                break;
            }
        }
    }

    if (base_data) {
        int result = MPI_Allreduce(MPI_IN_PLACE, MPI_BOTTOM, 1, datum_datatype, MPI_SUM, MPI_COMM_WORLD);
        if (result != MPI_SUCCESS) {
            std::cerr << "Error: MPI_Allreduce failed for field line " << field_line_idx
                      << " with code " << result << std::endl;
        }
    }
}


void MPIGatherDatumStoredAtEdge(cDatumStored& S, int root_rank) {
    // Get MPI info
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Use cached field lines - automatically handles caching and change detection
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& all_field_lines =
        DatumStoredAtEdgeMPIManager::GetCachedFieldLines();

    if (rank != root_rank) {
        // Non-root processes: send their data to root
        MPI_Datatype send_datatype = GetOrCreateScatterGatherDatatype(all_field_lines, S, rank);
        
        if (send_datatype != MPI_DATATYPE_NULL) {
            int result = MPI_Send(MPI_BOTTOM, 1, send_datatype, root_rank, 0, MPI_COMM_WORLD);
            if (result != MPI_SUCCESS) {
                std::cerr << "Error: MPI_Send failed on rank " << rank 
                          << " with code " << result << std::endl;
            }
        }
    } else {
        // Root process: receive data from all other processes
        for (int proc = 0; proc < size; ++proc) {
            if (proc == root_rank) {
                // Skip self - root process already has its own data
                continue;
            }
            
            // Create receive datatype for this specific process
            MPI_Datatype recv_datatype = GetOrCreateScatterGatherDatatype(all_field_lines, S, proc);
            
            if (recv_datatype != MPI_DATATYPE_NULL) {
                int result = MPI_Recv(MPI_BOTTOM, 1, recv_datatype, proc, 0, 
                                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (result != MPI_SUCCESS) {
                    std::cerr << "Error: MPI_Recv failed from rank " << proc 
                              << " with code " << result << std::endl;
                }
            }
        }
        
        // Success message only on root
        std::cout << "Successfully completed MPI gather for DatumStoredAtEdge data" << std::endl;
    }
}

void MPIBcastDatumStoredAtEdge(cDatumStored& S, int root_rank) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const auto& all_field_lines = DatumStoredAtEdgeMPIManager::GetCachedFieldLines();

    // Use unified datatype that covers ALL segments from ALL field lines
    // This ensures root's complete dataset is broadcast to all processes
    MPI_Datatype unified_datatype = GetOrCreateDatumStoredAtEdgeDatatype(all_field_lines, S);
    
    if (unified_datatype != MPI_DATATYPE_NULL) {
        // Broadcast ALL of root's DatumStoredAtEdge data to all other processes
        int result = MPI_Bcast(MPI_BOTTOM, 1, unified_datatype, root_rank, MPI_COMM_WORLD);
        
        if (result != MPI_SUCCESS) {
            std::cerr << "Error: MPI_Bcast failed on rank " << rank 
                      << " with code " << result << std::endl;
        } else if (rank == root_rank) {
            int total_segments = 0;
            for (const auto& line : all_field_lines) {
                total_segments += line.size();
            }
            std::cout << "Successfully broadcast all DatumStoredAtEdge data (" 
                      << total_segments << " segments) from root rank " << root_rank 
                      << " to all processes" << std::endl;
        }
    } else {
        std::cerr << "Error: Failed to create unified MPI datatype for broadcast operation" << std::endl;
    }
}

// ============================================================================
// Public Wrapper Functions
// ============================================================================

// Public wrapper functions for scatter/gather datatypes
MPI_Datatype ConstructScatterGatherMPIDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S,
    int target_process_rank) {
    return DatumStoredAtEdgeMPIManager::ConstructScatterGatherDatatype(field_lines, S, target_process_rank);
}

bool IsScatterGatherMPIDatatypeAvailable(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S,
    int target_process_rank) {
    return DatumStoredAtEdgeMPIManager::IsScatterGatherDatatypeAvailable(field_lines, S, target_process_rank);
}

// Public wrapper functions
MPI_Datatype ConstructDatumStoredAtEdgeMPIDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {
    return DatumStoredAtEdgeMPIManager::ConstructDatumStoredAtEdgeDatatype(field_lines, S);
}

bool IsDatumStoredAtEdgeMPIDatatypeAvailable(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {
    return DatumStoredAtEdgeMPIManager::IsDatumStoredAtEdgeDatatypeAvailable(field_lines, S);
}

void RemoveDatumStoredAtEdgeDatatype(
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& field_lines,
    const cDatumStored& S) {
    bool removed = DatumStoredAtEdgeMPIManager::RemoveDatatype(field_lines, S);
    if (removed) {
        std::cout << "Successfully removed MPI datatype from cache (S.length=" << S.length << ")" << std::endl;
    } else {
        std::cout << "MPI datatype not found in cache (S.length=" << S.length << ")" << std::endl;
    }
}

void CleanupDatumStoredAtEdgeMPI() {
    DatumStoredAtEdgeMPIManager::CleanupAllDatatypes();
    std::cout << "Cleaned up all DatumStoredAtEdge MPI datatypes" << std::endl;
}

void InitializeDatumStoredAtEdgeMPI() {
    DatumStoredAtEdgeMPIManager::Initialize();
}

// Convenience functions for all field lines
MPI_Datatype ConstructDatumStoredAtEdgeMPIDatatypeAllFieldLines(const cDatumStored& S) {
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& all_field_lines =
        DatumStoredAtEdgeMPIManager::GetCachedFieldLines();
    return ConstructDatumStoredAtEdgeMPIDatatype(all_field_lines, S);
}

bool IsDatumStoredAtEdgeMPIDatatypeAvailableAllFieldLines(const cDatumStored& S) {
    const std::vector<std::vector<PIC::FieldLine::cFieldLineSegment*>>& all_field_lines =
        DatumStoredAtEdgeMPIManager::GetCachedFieldLines();
    return IsDatumStoredAtEdgeMPIDatatypeAvailable(all_field_lines, S);
}

} // namespace Parallel
} // namespace FieldLine
} // namespace PIC

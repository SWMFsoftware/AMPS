/*
PURPOSE:
--------
MPI communication functions for DatumStoredAtVertex data across field line vertices.
Uses direct field line traversal with no state management for maximum simplicity.
Mirrors the functionality of the Edge-based functions but operates on vertex data.

KEY FEATURES:
-------------
1. NO STATE: No initialization, cleanup, or caching - just pure functions
2. DIRECT TRAVERSAL: Uses standard loops over field lines and vertices
3. BUFFER-BASED: Uses std::vector<double> for all MPI operations
4. ZERO SETUP: No setup/teardown required - just call the functions
5. CONSISTENT PATTERN: All functions use same pack/MPI/unpack strategy

CORE FUNCTIONS:
---------------
- MPIAllReduceDatumStoredAtVertex(S)           // Sum across all processes
- MPIReduceDatumStoredAtVertex(S, root_rank)   // Reduce to root process
- MPIBcastDatumStoredAtVertex(S, root_rank)    // Broadcast from root
- MPIAllReduceDatumStoredAtVertexFieldLine(field_line_idx, S)    // All-reduce single field line
- MPIReduceDatumStoredAtVertexFieldLine(field_line_idx, S, root_rank)   // Reduce single field line
- MPIBcastDatumStoredAtVertexFieldLine(field_line_idx, S, root_rank)    // Broadcast single field line

USAGE EXAMPLES:
===============

EXAMPLE 1: Sum magnetic field values across all processes
---------------------------------------------------------
```cpp
// Each process has computed local magnetic field data at vertices
PIC::Datum::cDatum* magneticFieldData;
// ... local computation fills magneticFieldData ...

// Sum all values across all processes (result on all processes)
PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(magneticFieldData);

// Now each process has the global sum at all vertices
std::cout << "Global magnetic field data computed" << std::endl;
```

EXAMPLE 2: Distribute plasma parameters from master process
-----------------------------------------------------------
```cpp
int rank;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

PIC::Datum::cDatum* plasmaParameters;

if (rank == 0) {
    // Master process sets up plasma parameters at vertices
    std::cout << "Setting up plasma parameters..." << std::endl;
    // ... fill plasmaParameters with setup data ...
}

// Broadcast from rank 0 to all other processes
PIC::FieldLine::Parallel::MPIBcastDatumStoredAtVertex(plasmaParameters, 0);

// Now all processes have the plasma parameters
std::cout << "Process " << rank << " received plasma parameters" << std::endl;
```

EXAMPLE 3: Work with specific field line vertex data
----------------------------------------------------
```cpp
// Sum temperature data for field line 2 only across all processes
PIC::Datum::cDatum* temperatureData;
PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertexFieldLine(2, temperatureData);

// Reduce field line 1 density data to rank 0
PIC::Datum::cDatum* densityData;
PIC::FieldLine::Parallel::MPIReduceDatumStoredAtVertexFieldLine(1, densityData, 0);

// Broadcast field line 0 velocity data from rank 0 to all processes
PIC::Datum::cDatum* velocityData;
PIC::FieldLine::Parallel::MPIBcastDatumStoredAtVertexFieldLine(0, velocityData, 0);
```
*/

#include "pic.h"
#include <iostream>
#include <vector>

namespace PIC {
namespace FieldLine {
namespace Parallel {

// ============================================================================
// Helper functions for vertex buffer operations - no caching, direct traversal
// ============================================================================

// Count total vertices across all field lines
int CountTotalVertices() {
    int total_vertices = 0;
    
    for (int iFieldLine = 0; iFieldLine < nFieldLine; iFieldLine++) {
        for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[iFieldLine].GetFirstSegment();
             Segment != NULL; 
             Segment = Segment->GetNext()) {
            
            // Count begin vertex (end vertex will be counted as begin of next segment)
            total_vertices++;
        }
        
        // Add the final vertex of the last segment
        PIC::FieldLine::cFieldLineSegment* LastSegment = FieldLinesAll[iFieldLine].GetLastSegment();
        if (LastSegment != NULL) {
            total_vertices++;
        }
    }
    
    return total_vertices;
}

// Count vertices for a specific field line
int CountFieldLineVertices(int field_line_idx) {
    if (field_line_idx < 0 || field_line_idx >= nFieldLine) {
        return 0;
    }
    
    int vertices = 0;
    for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[field_line_idx].GetFirstSegment();
         Segment != NULL; 
         Segment = Segment->GetNext()) {
        
        // Count begin vertex
        vertices++;
    }
    
    // Add the final vertex of the last segment
    PIC::FieldLine::cFieldLineSegment* LastSegment = FieldLinesAll[field_line_idx].GetLastSegment();
    if (LastSegment != NULL) {
        vertices++;
    }
    
    return vertices;
}

// Count vertices for a specific process (operates on all segments regardless of Thread)
int CountProcessVertices(int target_process) {
    int process_vertices = 0;
    
    for (int iFieldLine = 0; iFieldLine < nFieldLine; iFieldLine++) {
        for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[iFieldLine].GetFirstSegment();
             Segment != NULL; 
             Segment = Segment->GetNext()) {
            
            // Count begin vertex (regardless of Segment->Thread)
            process_vertices++;
            
            // If this is the last segment, also count end vertex
            if (Segment->GetNext() == NULL) {
                process_vertices++;
            }
        }
    }
    
    return process_vertices;
}

// Pack all field line vertex data into a contiguous buffer
int PackAllFieldLinesVertexData(PIC::Datum::cDatum* S, std::vector<double>& buffer) {
    int total_vertices = CountTotalVertices();
    int total_elements = total_vertices * S->length;
    
    buffer.resize(total_elements, 0.0);
    
    if (total_vertices == 0) {
        return 0;
    }
    
    // Pack data using direct traversal
    int element_index = 0;
    
    for (int iFieldLine = 0; iFieldLine < nFieldLine; iFieldLine++) {
        for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[iFieldLine].GetFirstSegment();
             Segment != NULL; 
             Segment = Segment->GetNext()) {
            
            // Pack begin vertex data
            PIC::FieldLine::cFieldLineVertex* BeginVertex = Segment->GetBegin();
            if (BeginVertex) {
                double* vertex_data = BeginVertex->GetDatum_ptr(S,BeginVertex->GetCompletedSamplingOffset());
                if (vertex_data) {
                    for (int i = 0; i < S->length; i++) {
                        buffer[element_index++] = vertex_data[i];
                        
                        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                            validate_numeric(vertex_data[i], __LINE__, __FILE__);
                        }
                    }
                } else {
                    // Leave as zeros if no data
                    element_index += S->length;
                }
            } else {
                element_index += S->length;
            }
            
            // If this is the last segment, also pack end vertex
            if (Segment->GetNext() == NULL) {
                PIC::FieldLine::cFieldLineVertex* EndVertex = Segment->GetEnd();
                if (EndVertex) {
                    double* vertex_data = EndVertex->GetDatum_ptr(S,EndVertex->GetCompletedSamplingOffset());
                    if (vertex_data) {
                        for (int i = 0; i < S->length; i++) {
                            buffer[element_index++] = vertex_data[i];
                            
                            if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                                validate_numeric(vertex_data[i], __LINE__, __FILE__);
                            }
                        }
                    } else {
                        element_index += S->length;
                    }
                } else {
                    element_index += S->length;
                }
            }
        }
    }
    
    return total_vertices;
}

// Unpack buffer data back to field line vertices
void UnpackAllFieldLinesVertexData(PIC::Datum::cDatum* S, const std::vector<double>& buffer) {
    int element_index = 0;
    
    for (int iFieldLine = 0; iFieldLine < nFieldLine; iFieldLine++) {
        for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[iFieldLine].GetFirstSegment();
             Segment != NULL; 
             Segment = Segment->GetNext()) {
            
            // Unpack begin vertex data
            PIC::FieldLine::cFieldLineVertex* BeginVertex = Segment->GetBegin();
            if (BeginVertex) {
                double* vertex_data = BeginVertex->GetDatum_ptr(S,BeginVertex->GetCompletedSamplingOffset());
                if (vertex_data) {
                    for (int i = 0; i < S->length; i++) {
                        vertex_data[i] = buffer[element_index++];
                        
                        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                            validate_numeric(vertex_data[i], __LINE__, __FILE__);
                        }
                    }
                } else {
                    // Skip if no data pointer
                    element_index += S->length;
                }
            } else {
                element_index += S->length;
            }
            
            // If this is the last segment, also unpack end vertex
            if (Segment->GetNext() == NULL) {
                PIC::FieldLine::cFieldLineVertex* EndVertex = Segment->GetEnd();
                if (EndVertex) {
                    double* vertex_data = EndVertex->GetDatum_ptr(S,EndVertex->GetCompletedSamplingOffset());
                    if (vertex_data) {
                        for (int i = 0; i < S->length; i++) {
                            vertex_data[i] = buffer[element_index++];
                            
                            if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                                validate_numeric(vertex_data[i], __LINE__, __FILE__);
                            }
                        }
                    } else {
                        element_index += S->length;
                    }
                } else {
                    element_index += S->length;
                }
            }
        }
    }
}

// Pack data for specific field line vertices
int PackFieldLineVertexData(int field_line_idx, PIC::Datum::cDatum* S, std::vector<double>& buffer) {
    if (field_line_idx < 0 || field_line_idx >= nFieldLine) {
        buffer.clear();
        return 0;
    }
    
    int vertices = CountFieldLineVertices(field_line_idx);
    int total_elements = vertices * S->length;
    
    buffer.resize(total_elements, 0.0);
    
    if (vertices == 0) {
        return 0;
    }
    
    // Pack data for this field line
    int element_index = 0;
    
    for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[field_line_idx].GetFirstSegment();
         Segment != NULL; 
         Segment = Segment->GetNext()) {
        
        // Pack begin vertex data
        PIC::FieldLine::cFieldLineVertex* BeginVertex = Segment->GetBegin();
        if (BeginVertex) {
            double* vertex_data = BeginVertex->GetDatum_ptr(S);
            if (vertex_data) {
                for (int i = 0; i < S->length; i++) {
                    buffer[element_index++] = vertex_data[i];
                }
            } else {
                element_index += S->length;  // Leave as zeros
            }
        } else {
            element_index += S->length;
        }
        
        // If this is the last segment, also pack end vertex
        if (Segment->GetNext() == NULL) {
            PIC::FieldLine::cFieldLineVertex* EndVertex = Segment->GetEnd();
            if (EndVertex) {
                double* vertex_data = EndVertex->GetDatum_ptr(S);
                if (vertex_data) {
                    for (int i = 0; i < S->length; i++) {
                        buffer[element_index++] = vertex_data[i];
                    }
                } else {
                    element_index += S->length;
                }
            } else {
                element_index += S->length;
            }
        }
    }
    
    return vertices;
}

// Unpack data for specific field line vertices
void UnpackFieldLineVertexData(int field_line_idx, PIC::Datum::cDatum* S, const std::vector<double>& buffer) {
    if (field_line_idx < 0 || field_line_idx >= nFieldLine) {
        return;
    }
    
    int element_index = 0;
    
    for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[field_line_idx].GetFirstSegment();
         Segment != NULL; 
         Segment = Segment->GetNext()) {
        
        // Unpack begin vertex data
        PIC::FieldLine::cFieldLineVertex* BeginVertex = Segment->GetBegin();
        if (BeginVertex) {
            double* vertex_data = BeginVertex->GetDatum_ptr(S);
            if (vertex_data) {
                for (int i = 0; i < S->length; i++) {
                    vertex_data[i] = buffer[element_index++];
                }
            } else {
                element_index += S->length;  // Skip
            }
        } else {
            element_index += S->length;
        }
        
        // If this is the last segment, also unpack end vertex
        if (Segment->GetNext() == NULL) {
            PIC::FieldLine::cFieldLineVertex* EndVertex = Segment->GetEnd();
            if (EndVertex) {
                double* vertex_data = EndVertex->GetDatum_ptr(S);
                if (vertex_data) {
                    for (int i = 0; i < S->length; i++) {
                        vertex_data[i] = buffer[element_index++];
                    }
                } else {
                    element_index += S->length;
                }
            } else {
                element_index += S->length;
            }
        }
    }
}

// ============================================================================
// MPI Operations for Vertex Data - Ultra-simplified implementation
// ============================================================================

void MPIAllReduceDatumStoredAtVertex(PIC::Datum::cDatum* S) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Pack all vertex data into buffer
    std::vector<double> buffer;
    int total_vertices = PackAllFieldLinesVertexData(S, buffer);
    
    if (total_vertices == 0) {
        if (rank == 0) {
            std::cerr << "Warning: No vertices found for all-reduce operation" << std::endl;
        }
        return;
    }
    
    // Single MPI_Allreduce call
    int result = MPI_Allreduce(MPI_IN_PLACE, buffer.data(), buffer.size(), 
                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Allreduce failed with code " << result << std::endl;
        return;
    }
    
    // Unpack results back to vertices
    UnpackAllFieldLinesVertexData(S, buffer);
    
    if (rank == 0) {
        std::cout << "Completed MPI all-reduce for " << total_vertices 
                  << " vertices (" << buffer.size() << " elements)" << std::endl;
    }
}

void MPIReduceDatumStoredAtVertex(PIC::Datum::cDatum* S, int root_rank) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (root_rank < 0 || root_rank >= size) {
        std::cerr << "Error: Invalid root_rank " << root_rank << std::endl;
        return;
    }

    // Pack send data
    std::vector<double> send_buffer;
    int total_vertices = PackAllFieldLinesVertexData(S, send_buffer);
    
    if (total_vertices == 0) {
        if (rank == 0) {
            std::cerr << "Warning: No vertices found for reduce operation" << std::endl;
        }
        return;
    }
    
    // Prepare receive buffer (only on root)
    std::vector<double> recv_buffer;
    if (rank == root_rank) {
        recv_buffer.resize(send_buffer.size(), 0.0);
    }
    
    // MPI_Reduce operation
    int result = MPI_Reduce(send_buffer.data(), 
                           (rank == root_rank) ? recv_buffer.data() : nullptr,
                           send_buffer.size(), MPI_DOUBLE, MPI_SUM, 
                           root_rank, MPI_COMM_WORLD);
    
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Reduce failed with code " << result << std::endl;
        return;
    }
    
    // Unpack results (only on root)
    if (rank == root_rank) {
        UnpackAllFieldLinesVertexData(S, recv_buffer);
        std::cout << "Completed MPI reduce for " << total_vertices 
                  << " vertices to root rank " << root_rank << std::endl;
    }
}

void MPIBcastDatumStoredAtVertex(PIC::Datum::cDatum* S, int root_rank) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (root_rank < 0 || root_rank >= size) {
        std::cerr << "Error: Invalid root_rank " << root_rank << std::endl;
        return;
    }

    // Pack data into buffer
    std::vector<double> buffer;
    int total_vertices = PackAllFieldLinesVertexData(S, buffer);
    
    if (total_vertices == 0) {
        if (rank == 0) {
            std::cerr << "Warning: No vertices found for broadcast operation" << std::endl;
        }
        return;
    }
    
    // Broadcast buffer
    int result = MPI_Bcast(buffer.data(), buffer.size(), MPI_DOUBLE, 
                          root_rank, MPI_COMM_WORLD);
    
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Bcast failed with code " << result << std::endl;
        return;
    }
    
    // Unpack data back to vertices
    UnpackAllFieldLinesVertexData(S, buffer);
    
    if (rank == root_rank) {
        std::cout << "Completed MPI broadcast for " << total_vertices 
                  << " vertices from root rank " << root_rank << std::endl;
    }
}

void MPIAllReduceDatumStoredAtVertexFieldLine(int field_line_idx, PIC::Datum::cDatum* S) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (field_line_idx < 0 || field_line_idx >= nFieldLine) {
        std::cerr << "Error: Invalid field_line_idx " << field_line_idx << std::endl;
        return;
    }

    // Pack data for this field line's vertices
    std::vector<double> buffer;
    int vertex_count = PackFieldLineVertexData(field_line_idx, S, buffer);
    
    if (vertex_count == 0) {
        if (rank == 0) {
            std::cerr << "Warning: No vertices found for field line " << field_line_idx << std::endl;
        }
        return;
    }
    
    // All-reduce for this field line's vertices
    int result = MPI_Allreduce(MPI_IN_PLACE, buffer.data(), buffer.size(), 
                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Allreduce failed for field line " << field_line_idx 
                  << " with code " << result << std::endl;
        return;
    }
    
    // Unpack results
    UnpackFieldLineVertexData(field_line_idx, S, buffer);
    
    if (rank == 0) {
        std::cout << "Completed MPI all-reduce for field line " << field_line_idx 
                  << " (" << vertex_count << " vertices)" << std::endl;
    }
}

void MPIReduceDatumStoredAtVertexFieldLine(int field_line_idx, PIC::Datum::cDatum* S, int root_rank) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (field_line_idx < 0 || field_line_idx >= nFieldLine) {
        std::cerr << "Error: Invalid field_line_idx " << field_line_idx << std::endl;
        return;
    }
    
    if (root_rank < 0 || root_rank >= size) {
        std::cerr << "Error: Invalid root_rank " << root_rank << std::endl;
        return;
    }

    // Pack send data for this field line's vertices
    std::vector<double> send_buffer;
    int vertex_count = PackFieldLineVertexData(field_line_idx, S, send_buffer);
    
    if (vertex_count == 0) {
        if (rank == 0) {
            std::cerr << "Warning: No vertices found for field line " << field_line_idx << std::endl;
        }
        return;
    }
    
    // Prepare receive buffer (only on root)
    std::vector<double> recv_buffer;
    if (rank == root_rank) {
        recv_buffer.resize(send_buffer.size(), 0.0);
    }
    
    // MPI_Reduce operation for this field line's vertices
    int result = MPI_Reduce(send_buffer.data(), 
                           (rank == root_rank) ? recv_buffer.data() : nullptr,
                           send_buffer.size(), MPI_DOUBLE, MPI_SUM, 
                           root_rank, MPI_COMM_WORLD);
    
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Reduce failed for field line " << field_line_idx 
                  << " with code " << result << std::endl;
        return;
    }
    
    // Unpack results (only on root)
    if (rank == root_rank) {
        UnpackFieldLineVertexData(field_line_idx, S, recv_buffer);
        std::cout << "Completed MPI reduce for field line " << field_line_idx 
                  << " (" << vertex_count << " vertices) to root rank " << root_rank << std::endl;
    }
}

void MPIBcastDatumStoredAtVertexFieldLine(int field_line_idx, PIC::Datum::cDatum* S, int root_rank) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (field_line_idx < 0 || field_line_idx >= nFieldLine) {
        std::cerr << "Error: Invalid field_line_idx " << field_line_idx << std::endl;
        return;
    }
    
    if (root_rank < 0 || root_rank >= size) {
        std::cerr << "Error: Invalid root_rank " << root_rank << std::endl;
        return;
    }

    // Pack data for this field line's vertices
    std::vector<double> buffer;
    int vertex_count = PackFieldLineVertexData(field_line_idx, S, buffer);
    
    if (vertex_count == 0) {
        if (rank == 0) {
            std::cerr << "Warning: No vertices found for field line " << field_line_idx << std::endl;
        }
        return;
    }
    
    // Broadcast buffer for this field line's vertices
    int result = MPI_Bcast(buffer.data(), buffer.size(), MPI_DOUBLE, 
                          root_rank, MPI_COMM_WORLD);
    
    if (result != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Bcast failed for field line " << field_line_idx 
                  << " with code " << result << std::endl;
        return;
    }
    
    // Unpack data back to this field line's vertices
    UnpackFieldLineVertexData(field_line_idx, S, buffer);
    
    if (rank == root_rank) {
        std::cout << "Completed MPI broadcast for field line " << field_line_idx 
                  << " (" << vertex_count << " vertices) from root rank " << root_rank << std::endl;
    }
}

// ============================================================================
// Utility functions for setting vertex data across all field lines
// ============================================================================

void SetDatumStoredAtVertex(double val, PIC::Datum::cDatum* datum) {
    // Set a single value to all vertices across all field lines
    
    for (int iFieldLine = 0; iFieldLine < nFieldLine; iFieldLine++) {
        for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[iFieldLine].GetFirstSegment();
             Segment != NULL; 
             Segment = Segment->GetNext()) {
            
            // Set begin vertex data
            PIC::FieldLine::cFieldLineVertex* BeginVertex = Segment->GetBegin();
            if (BeginVertex) {
                double* vertex_data = BeginVertex->GetDatum_ptr(datum);
                if (vertex_data) {
                    for (int i = 0; i < datum->length; i++) {
                        vertex_data[i] = val;
                        
                        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                            validate_numeric(vertex_data[i], __LINE__, __FILE__);
                        }
                    }
                }
            }
            
            // If this is the last segment, also set end vertex
            if (Segment->GetNext() == NULL) {
                PIC::FieldLine::cFieldLineVertex* EndVertex = Segment->GetEnd();
                if (EndVertex) {
                    double* vertex_data = EndVertex->GetDatum_ptr(datum);
                    if (vertex_data) {
                        for (int i = 0; i < datum->length; i++) {
                            vertex_data[i] = val;
                            
                            if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                                validate_numeric(vertex_data[i], __LINE__, __FILE__);
                            }
                        }
                    }
                }
            }
        }
    }
}

void SetDatumStoredAtVertex(double* val, PIC::Datum::cDatum* datum) {
    // Set an array of values to all vertices across all field lines
    if (val == nullptr) {
        std::cerr << "Error: SetDatumStoredAtVertex received null pointer for val array" << std::endl;
        return;
    }
    
    for (int iFieldLine = 0; iFieldLine < nFieldLine; iFieldLine++) {
        for (PIC::FieldLine::cFieldLineSegment* Segment = FieldLinesAll[iFieldLine].GetFirstSegment();
             Segment != NULL; 
             Segment = Segment->GetNext()) {
            
            // Set begin vertex data
            PIC::FieldLine::cFieldLineVertex* BeginVertex = Segment->GetBegin();
            if (BeginVertex) {
                double* vertex_data = BeginVertex->GetDatum_ptr(datum);
                if (vertex_data) {
                    for (int i = 0; i < datum->length; i++) {
                        vertex_data[i] = val[i];
                        
                        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                            validate_numeric(vertex_data[i], __LINE__, __FILE__);
                        }
                    }
                }
            }
            
            // If this is the last segment, also set end vertex
            if (Segment->GetNext() == NULL) {
                PIC::FieldLine::cFieldLineVertex* EndVertex = Segment->GetEnd();
                if (EndVertex) {
                    double* vertex_data = EndVertex->GetDatum_ptr(datum);
                    if (vertex_data) {
                        for (int i = 0; i < datum->length; i++) {
                            vertex_data[i] = val[i];
                            
                            if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                                validate_numeric(vertex_data[i], __LINE__, __FILE__);
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace Parallel
} // namespace FieldLine
} // namespace PIC

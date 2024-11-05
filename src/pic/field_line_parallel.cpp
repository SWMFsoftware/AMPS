
#include "pic.h"

#include <unordered_set>
#include <stdexcept>
#include <vector>
#include <map>


PIC::ParallelFieldLines::cThreadSegmentTable* PIC::ParallelFieldLines::ThreadSegmentTable;

void PIC::ParallelFieldLines::StaticDecompositionSegmentNumber() {
namespace FL = PIC::FieldLine;
  int iFieldLine;

  if (FL::FieldLinesAll==NULL) return;

  for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    StaticDecompositionSegmentNumber(FL::FieldLinesAll[iFieldLine].GetFirstSegment());
  }

  GenerateThreadSegmentTable();
} 



void PIC::ParallelFieldLines::StaticDecompositionSegmentNumber(PIC::FieldLine::cFieldLineSegment* SegmentIn) { 
    // Count segments
    int nSegments = 0;
    auto Segment = SegmentIn;
    while (Segment != NULL) {
        nSegments++;
        Segment = Segment->GetNext();
    }
    
    // Handle edge case
    if (nSegments == 0) return;
    
    // Calculate segments per thread
    int nSegmentsPerThread = std::max(1, nSegments / PIC::nTotalThreads);
    
    // Assign threads
    int cnt = 0;
    Segment = SegmentIn;
    while (Segment != NULL) {
        Segment->Thread = std::min(cnt / nSegmentsPerThread, 
                                 PIC::nTotalThreads - 1);
        cnt++;
        Segment = Segment->GetNext();
    }
}

void PIC::ParallelFieldLines::StaticDecompositionFieldLineLength() {
namespace FL = PIC::FieldLine;
  int iFieldLine;

  if (FL::FieldLinesAll==NULL) return;
  
  for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    StaticDecompositionFieldLineLength(FL::FieldLinesAll[iFieldLine].GetFirstSegment());
  }

  GenerateThreadSegmentTable();
} 


void PIC::ParallelFieldLines::StaticDecompositionFieldLineLength(PIC::FieldLine::cFieldLineSegment* SegmentIn) {
    // First pass: calculate total length
    double totalLength = 0.0;
    auto Segment = SegmentIn;
    while (Segment != NULL) {
        totalLength += Segment->GetLength();
        Segment = Segment->GetNext();
    }

    // Handle edge cases
    if (totalLength == 0.0 || SegmentIn == NULL) return;

    // Calculate target length per thread
    double lengthPerThread = totalLength / PIC::nTotalThreads;

    // Second pass: assign threads based on cumulative length
    double currentLength = 0.0;
    Segment = SegmentIn;
    while (Segment != NULL) {
        // Add half of current segment's length to better distribute boundary segments
        double decisionPoint = currentLength + (Segment->GetLength() / 2.0);

        // Calculate thread number based on accumulated length
        int threadNum = static_cast<int>(decisionPoint / lengthPerThread);

        // Ensure thread number is within bounds
        Segment->Thread = std::min(threadNum, PIC::nTotalThreads - 1);

        // Update accumulated length
        currentLength += Segment->GetLength();
        Segment = Segment->GetNext();
    }
}


void PIC::ParallelFieldLines::GenerateThreadSegmentTable() {
namespace FL = PIC::FieldLine;
  int i,iFieldLine,*cnt;

  if (ThreadSegmentTable!=NULL) exit(__LINE__,__FILE__,"Error: re-initialization of ThreadSegmentTable");

  if (FL::FieldLinesAll==NULL) return;
 
  ThreadSegmentTable=new cThreadSegmentTable[FL::nFieldLine];
  cnt=new int [PIC::nTotalThreads];


  //count the number of segments
  for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    for (i=0;i<PIC::nTotalThreads;i++) cnt[i]=0; 

    for (auto Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();Segment!=NULL;Segment=Segment->GetNext()) {
      cnt[Segment->Thread]++;
    }

    ThreadSegmentTable[iFieldLine].TableLength=new int [PIC::nTotalThreads];
    ThreadSegmentTable[iFieldLine].Table=new FL::cFieldLineSegment** [PIC::nTotalThreads];

    for (int thread=0;thread<PIC::nTotalThreads;thread++) {
      ThreadSegmentTable[iFieldLine].TableLength[thread]=cnt[thread];
      ThreadSegmentTable[iFieldLine].Table[thread]=(cnt[thread]!=0) ? new FL::cFieldLineSegment*[cnt[thread]] : NULL; 
    }

    //populate the table 
    for (i=0;i<PIC::nTotalThreads;i++) cnt[i]=0;

    for (auto Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();Segment!=NULL;Segment=Segment->GetNext()) {
      ThreadSegmentTable[iFieldLine].Table[Segment->Thread][cnt[Segment->Thread]]=Segment;
      cnt[Segment->Thread]++;
    }
  }

  delete [] cnt;
}

//===================================   EXCHANGE PARTICLES BETWEEN PROCESSES ============================================
void PIC::ParallelFieldLines::ExchangeFieldLineParticles(cThreadSegmentTable& ThreadSegmentTable) {
    namespace FL = PIC::FieldLine;
    MPI_Request *SendParticleDataRequest = new MPI_Request[PIC::nTotalThreads];
    MPI_Request *RecvParticleDataRequest = new MPI_Request[PIC::nTotalThreads];
    int SendRequestLength = 0, RecvRequestLength = 0;


    PIC::ParallelFieldLines::GetFieldLinePopulationStat();

    // Enhanced descriptor that keeps track of particles to be deleted
    struct cSegmentParticleDescriptor {
        int segmentIndex;
        int nTotalParticles;
        std::vector<long int> particleList;
        FL::cFieldLineSegment* segment;
    };

    std::map<int, std::vector<cSegmentParticleDescriptor>> ProcessParticleMap;

    // First phase: Count particles and prepare send lists
    // This needs to be done by all processes before communication starts
    for (int thread = 0; thread < PIC::nTotalThreads; thread++) {
        if (thread == PIC::ThisThread) continue;
        
        std::vector<cSegmentParticleDescriptor>& segmentList = ProcessParticleMap[thread];
        
        for (int i = 0; i < ThreadSegmentTable.TableLength[thread]; i++) {
            FL::cFieldLineSegment* segment = ThreadSegmentTable.Table[thread][i];
            
            if (segment && segment->FirstParticleIndex != -1) {
                cSegmentParticleDescriptor descriptor;
                descriptor.segmentIndex = i;
                descriptor.nTotalParticles = 0;
                descriptor.segment = segment;
                
                long int particlePtr = segment->FirstParticleIndex;
                while (particlePtr != -1) {
                    descriptor.particleList.push_back(particlePtr);
                    descriptor.nTotalParticles++;
                    particlePtr = PIC::ParticleBuffer::GetNext(particlePtr);
                }
                
                if (descriptor.nTotalParticles > 0) {
                    segmentList.push_back(descriptor);
                }
            }
        }
    }

    // Second phase: Exchange particle counts with all processes
    std::vector<int> sendCounts(PIC::nTotalThreads, 0);
    std::vector<int> recvCounts(PIC::nTotalThreads, 0);
    
    // Prepare send counts
    for (const auto& pair : ProcessParticleMap) {
        int totalParticles = 0;
        for (const auto& desc : pair.second) {
            totalParticles += desc.nTotalParticles;
        }
        sendCounts[pair.first] = totalParticles;
    }

    // Exchange counts using MPI_Alltoall
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT,
                 recvCounts.data(), 1, MPI_INT,
                 MPI_GLOBAL_COMMUNICATOR);

    // Third phase: Exchange particle data
    for (int thread = 0; thread < PIC::nTotalThreads; thread++) {
        if (thread == PIC::ThisThread) continue;

        // Post receives first
        if (recvCounts[thread] > 0) {
            // Allocate space for incoming particles
            std::vector<MPI_Aint> recvOffsets;
            recvOffsets.reserve(recvCounts[thread]);
            
            for (int i = 0; i < recvCounts[thread]; i++) {
                long int newParticle = PIC::ParticleBuffer::GetNewParticle();
                if (newParticle == -1) {
                    throw std::runtime_error("Failed to allocate new particle");
                }
                
                MPI_Aint offset;
                MPI_Get_address(PIC::ParticleBuffer::ParticleDataBuffer + 
                              PIC::ParticleBuffer::GetParticleDataOffset(newParticle),
                              &offset);
                recvOffsets.push_back(offset);
            }

            MPI_Datatype particle_recv_type;
            MPI_Type_create_hindexed_block(recvCounts[thread],
                                         PIC::ParticleBuffer::ParticleDataLength - 2*sizeof(long int),
                                         recvOffsets.data(),
                                         MPI_BYTE,
                                         &particle_recv_type);
            MPI_Type_commit(&particle_recv_type);

            MPI_Irecv(MPI_BOTTOM, 1, particle_recv_type, thread, 0,
                      MPI_GLOBAL_COMMUNICATOR, &RecvParticleDataRequest[RecvRequestLength++]);

            MPI_Type_free(&particle_recv_type);
        }
    }

    // Post sends after all receives are posted
    for (const auto& pair : ProcessParticleMap) {
        int thread = pair.first;
        if (sendCounts[thread] > 0) {
            std::vector<MPI_Aint> sendOffsets;
            sendOffsets.reserve(sendCounts[thread]);
            
            for (const auto& desc : pair.second) {
                for (long int particlePtr : desc.particleList) {
                    MPI_Aint offset;
                    MPI_Get_address(PIC::ParticleBuffer::ParticleDataBuffer + 
                                  PIC::ParticleBuffer::GetParticleDataOffset(particlePtr),
                                  &offset);
                    sendOffsets.push_back(offset);
                }
            }

            MPI_Datatype particle_send_type;
            MPI_Type_create_hindexed_block(sendCounts[thread],
                                         PIC::ParticleBuffer::ParticleDataLength - 2*sizeof(long int),
                                         sendOffsets.data(),
                                         MPI_BYTE,
                                         &particle_send_type);
            MPI_Type_commit(&particle_send_type);

            MPI_Isend(MPI_BOTTOM, 1, particle_send_type, thread, 0,
                      MPI_GLOBAL_COMMUNICATOR, &SendParticleDataRequest[SendRequestLength++]);

            MPI_Type_free(&particle_send_type);
        }
    }

    // Wait for all communication to complete
    if (SendRequestLength > 0) {
        MPI_Waitall(SendRequestLength, SendParticleDataRequest, MPI_STATUSES_IGNORE);
    }
    if (RecvRequestLength > 0) {
        MPI_Waitall(RecvRequestLength, RecvParticleDataRequest, MPI_STATUSES_IGNORE);
    }

    // Final phase: Delete sent particles
    for (const auto& pair : ProcessParticleMap) {
        for (const auto& desc : pair.second) {
            std::unordered_set<long int> deleteSet;
            for (const auto& particleId : desc.particleList) {
                deleteSet.insert(particleId);
            }
            
            long int newFirst = -1;
            long int ptr = desc.segment->FirstParticleIndex;
            
            while (ptr != -1) {
                long int next = PIC::ParticleBuffer::GetNext(ptr);
                
                if (deleteSet.find(ptr) == deleteSet.end()) {
                    PIC::ParticleBuffer::SetNext(ptr, newFirst);
                    newFirst = ptr;
                } else {
                    PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr, true);
                }
                
                ptr = next;
            }
            
            desc.segment->FirstParticleIndex = newFirst;
        }
    }

    delete[] SendParticleDataRequest;
    delete[] RecvParticleDataRequest;
}

void PIC::ParallelFieldLines::ExchangeFieldLineParticles() {
namespace FL = PIC::FieldLine;
  if (FL::FieldLinesAll==NULL) return;

  if (ThreadSegmentTable==NULL) exit(__LINE__,__FILE__,"Error: ThreadSegmentTable is not initialized");

  PIC::ParallelFieldLines::GetFieldLinePopulationStat();

  //loop through field lines  
  for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    ExchangeFieldLineParticles(ThreadSegmentTable[iFieldLine]);
  }
}

//==================================   output the number of particles attached to each of the field lines ======================================
void PIC::ParallelFieldLines::GetFieldLinePopulationStat(PIC::FieldLine::cFieldLineSegment* FirstSegment) {
    // Statistics for this process
    struct ParticleStats {
        long int totalParticles;         // Total particles in field line
        long int localParticles;         // Particles in segments owned by this thread
        long int remoteParticles;        // Particles in segments owned by other threads
        int nLocalSegments;              // Number of segments owned by this thread
        int nRemoteSegments;             // Number of segments owned by other threads
    };
    
    // Collect local statistics
    ParticleStats localStats = {0, 0, 0, 0, 0};
    
    // Loop through all segments in the field line
    auto currentSegment = FirstSegment;
    while (currentSegment != nullptr) {
        int particleCount = 0;
        long int ptr = currentSegment->FirstParticleIndex;
        
        // Count particles in this segment
        while (ptr != -1) {
            particleCount++;
            ptr = PIC::ParticleBuffer::GetNext(ptr);
        }
        
        // Update statistics
        localStats.totalParticles += particleCount;
        
        if (currentSegment->Thread == PIC::ThisThread) {
            localStats.localParticles += particleCount;
            localStats.nLocalSegments++;
        } else {
            localStats.remoteParticles += particleCount;
            localStats.nRemoteSegments++;
        }
        
        currentSegment = currentSegment->GetNext();
    }
    
    // Structure to hold statistics from all processes
    struct GlobalParticleStats {
        long int totalParticles;
        long int localParticles;
        long int remoteParticles;
        int nLocalSegments;
        int nRemoteSegments;
        int processId;
    };
    
    std::vector<GlobalParticleStats> allStats;
    
    if (PIC::ThisThread == 0) {
        allStats.resize(PIC::nTotalThreads);
    }
    
    // Create a custom MPI type for ParticleStats
    MPI_Datatype mpi_stats_type;
    int blocklengths[] = {1, 1, 1, 1, 1, 1};
    MPI_Aint displacements[6];
    MPI_Datatype types[] = {MPI_LONG, MPI_LONG, MPI_LONG, MPI_INT, MPI_INT, MPI_INT};
    
    GlobalParticleStats sample;
    MPI_Aint base_address;
    MPI_Get_address(&sample, &base_address);
    MPI_Get_address(&sample.totalParticles, &displacements[0]);
    MPI_Get_address(&sample.localParticles, &displacements[1]);
    MPI_Get_address(&sample.remoteParticles, &displacements[2]);
    MPI_Get_address(&sample.nLocalSegments, &displacements[3]);
    MPI_Get_address(&sample.nRemoteSegments, &displacements[4]);
    MPI_Get_address(&sample.processId, &displacements[5]);
    
    for (int i = 0; i < 6; i++) {
        displacements[i] = MPI_Aint_diff(displacements[i], base_address);
    }
    
    MPI_Type_create_struct(6, blocklengths, displacements, types, &mpi_stats_type);
    MPI_Type_commit(&mpi_stats_type);
    
    // Prepare local stats for gathering
    GlobalParticleStats myStats = {
        localStats.totalParticles,
        localStats.localParticles,
        localStats.remoteParticles,
        localStats.nLocalSegments,
        localStats.nRemoteSegments,
        PIC::ThisThread
    };
    
    // Gather statistics from all processes
    MPI_Gather(&myStats, 1, mpi_stats_type,
               allStats.data(), 1, mpi_stats_type,
               0, MPI_GLOBAL_COMMUNICATOR);
    
    MPI_Type_free(&mpi_stats_type);
    
    // Process 0 prints the results
    if (PIC::ThisThread == 0) {
        // Calculate totals
        GlobalParticleStats totals = {0, 0, 0, 0, 0, -1};
        for (const auto& stats : allStats) {
            totals.totalParticles += stats.totalParticles;
            totals.localParticles += stats.localParticles;
            totals.remoteParticles += stats.remoteParticles;
            totals.nLocalSegments += stats.nLocalSegments;
            totals.nRemoteSegments += stats.nRemoteSegments;
        }
        
        // Print header
        printf("\n=== Field Line Particle Distribution Analysis ===\n");
        printf("Total Processes: %d\n\n", PIC::nTotalThreads);
        
        // Print per-process statistics
        printf("Per-Process Statistics:\n");
        printf("%-6s | %-12s | %-12s | %-12s | %-12s | %-12s\n",
               "Proc", "Total Parts", "Local Parts", "Remote Parts", "Local Segs", "Remote Segs");
        printf("----------------------------------------------------------------------\n");
        
        for (const auto& stats : allStats) {
            printf("%-6d | %-12ld | %-12ld | %-12ld | %-12d | %-12d\n",
                   stats.processId,
                   stats.totalParticles,
                   stats.localParticles,
                   stats.remoteParticles,
                   stats.nLocalSegments,
                   stats.nRemoteSegments);
        }
        
        // Print summary
        printf("\nGlobal Summary:\n");
        printf("Total Particles: %ld\n", totals.totalParticles);
        printf("Total Segments: %d\n", totals.nLocalSegments + totals.nRemoteSegments);
        printf("Average Particles per Process: %.2f\n", 
               static_cast<double>(totals.totalParticles) / PIC::nTotalThreads);
        printf("Average Segments per Process: %.2f\n", 
               static_cast<double>(totals.nLocalSegments + totals.nRemoteSegments) / PIC::nTotalThreads);
        
        // Print load balance metrics
        printf("\nLoad Balance Metrics:\n");
        long int maxParticles = 0, minParticles = LONG_MAX;
        int maxSegments = 0, minSegments = INT_MAX;
        
        for (const auto& stats : allStats) {
            maxParticles = std::max(maxParticles, stats.localParticles);
            minParticles = std::min(minParticles, stats.localParticles);
            maxSegments = std::max(maxSegments, stats.nLocalSegments);
            minSegments = std::min(minSegments, stats.nLocalSegments);
        }
        
        double particleImbalance = static_cast<double>(maxParticles - minParticles) / ((1>maxParticles) ? 1 : maxParticles)  * 100.0;
        double segmentImbalance = static_cast<double>(maxSegments - minSegments) / maxSegments * 100.0;
        
        printf("Particle Load Imbalance: %.2f%%\n", particleImbalance);
        printf("Segment Load Imbalance: %.2f%%\n", segmentImbalance);
        printf("================================================\n\n");
    }
}

void PIC::ParallelFieldLines::GetFieldLinePopulationStat() {
namespace FL = PIC::FieldLine;
  if (FL::FieldLinesAll==NULL) return;

  if (ThreadSegmentTable==NULL) exit(__LINE__,__FILE__,"Error: ThreadSegmentTable is not initialized");

  //loop through field lines
  for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    if (PIC::ThisThread==0) cout << "Field line particle statistic: field line " << iFieldLine << endl << flush;
    GetFieldLinePopulationStat(FL::FieldLinesAll[iFieldLine].GetFirstSegment());
  }
}

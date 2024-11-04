
#include "pic.h"

PIC::ParallelFieldLines::cThreadSegmentTable* PIC::ParallelFieldLines::ThreadSegmentTable;

void PIC::ParallelFieldLines::StaticDecompositionSegmentNumber() {
namespace FL = PIC::FieldLine;
  int iFieldLine;

  if (FL::FieldLinesAll==NULL) return;

  for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    StaticDecompositionSegmentNumber(FL::FieldLinesAll[iFieldLine].GetFirstSegment());
  }
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

  // Simplified descriptor using segment index instead of ID
  struct cSegmentParticleDescriptor {
    int segmentIndex;  // Index in ThreadSegmentTable.Table[thread]
    int nTotalParticles;
    std::vector<long int> particleList;
  };

  std::map<int, std::vector<cSegmentParticleDescriptor>> ProcessParticleMap;

  // Modified to use array indices
  auto PrepareParticleLists = [&](int targetThread) {
    std::vector<cSegmentParticleDescriptor>& segmentList = ProcessParticleMap[targetThread];
    
    // Loop through segments using array indices
    for (int i = 0; i < ThreadSegmentTable.TableLength[targetThread]; i++) {
      FL::cFieldLineSegment* segment = ThreadSegmentTable.Table[targetThread][i];
      
      if (segment->FirstParticleIndex != -1) {
        cSegmentParticleDescriptor descriptor;
        descriptor.segmentIndex = i;  // Using array index directly
        descriptor.nTotalParticles = 0;
        
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
    
    return segmentList.size();
  };

  auto InitSendParticleData = [&](int To, const std::vector<cSegmentParticleDescriptor>& descriptors) {
    int totalParticles = 0;
    for (const auto& desc : descriptors) {
      totalParticles += desc.nTotalParticles;
    }
    
    if (totalParticles == 0) return;

    MPI_Aint* offsets = new MPI_Aint[totalParticles];
    int offsetIndex = 0;
    
    for (const auto& desc : descriptors) {
      for (long int particlePtr : desc.particleList) {
        MPI_Get_address(PIC::ParticleBuffer::ParticleDataBuffer + 
                       PIC::ParticleBuffer::GetParticleDataOffset(particlePtr),
                       &offsets[offsetIndex++]);
      }
    }

    MPI_Datatype particle_send_type;
    MPI_Type_create_hindexed_block(totalParticles,
                                  PIC::ParticleBuffer::ParticleDataLength - 2*sizeof(long int),
                                  offsets,
                                  MPI_BYTE,
                                  &particle_send_type);
    MPI_Type_commit(&particle_send_type);

    MPI_Isend(MPI_BOTTOM, 1, particle_send_type, To, 0, MPI_GLOBAL_COMMUNICATOR, 
              &SendParticleDataRequest[SendRequestLength++]);

    MPI_Type_free(&particle_send_type);
    delete[] offsets;
  };

  // Modified to use array index
  auto InitNewParticles = [&](FL::cFieldLineSegment* segment, int nParticles) {
    std::vector<MPI_Aint> particleOffsets;
    particleOffsets.reserve(nParticles);
    
    for (int i = 0; i < nParticles; i++) {
      long int newParticle = PIC::ParticleBuffer::GetNewParticle(segment->FirstParticleIndex);
      
      MPI_Aint offset;
      MPI_Get_address(PIC::ParticleBuffer::ParticleDataBuffer + 
                     PIC::ParticleBuffer::GetParticleDataOffset(newParticle),
                     &offset);
      particleOffsets.push_back(offset);
    }
    
    return particleOffsets;
  };

  auto RemoveSentParticles = [&](FL::cFieldLineSegment* segment) {
    long int ptr = segment->FirstParticleIndex;
    
    while (ptr != -1) {
      long int next = PIC::ParticleBuffer::GetNext(ptr);
      PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr, true);
      ptr = next;
    }
    
    segment->FirstParticleIndex = -1;
  };

  // Prepare particle lists for each target process
  for (int thread = 0; thread < PIC::nTotalThreads; thread++) {
    if (thread != PIC::ThisThread) {
      PrepareParticleLists(thread);
    }
  }

  // Simplified exchange info structure using segment indices
  struct ParticleExchangeInfo {
    std::vector<int> segmentIndices;  // Array indices instead of IDs
    std::vector<int> particleCounts;
  };

  std::vector<ParticleExchangeInfo> recvInfo(PIC::nTotalThreads);
  
  // Exchange segment indices and particle counts
  for (const auto& pair : ProcessParticleMap) {
    int targetThread = pair.first;
    ParticleExchangeInfo info;
    
    for (const auto& desc : pair.second) {
      info.segmentIndices.push_back(desc.segmentIndex);
      info.particleCounts.push_back(desc.nTotalParticles);
    }
    
    MPI_Send(info.segmentIndices.data(), info.segmentIndices.size(), MPI_INT,
             targetThread, 1, MPI_GLOBAL_COMMUNICATOR);
    MPI_Send(info.particleCounts.data(), info.particleCounts.size(), MPI_INT,
             targetThread, 2, MPI_GLOBAL_COMMUNICATOR);
  }

  // Receive counts and prepare receiving buffers
  for (int thread = 0; thread < PIC::nTotalThreads; thread++) {
    if (thread == PIC::ThisThread) continue;

    MPI_Status status;
    int count;
    
    MPI_Probe(thread, 1, MPI_GLOBAL_COMMUNICATOR, &status);
    MPI_Get_count(&status, MPI_INT, &count);
    
    recvInfo[thread].segmentIndices.resize(count);
    recvInfo[thread].particleCounts.resize(count);
    
    MPI_Recv(recvInfo[thread].segmentIndices.data(), count, MPI_INT,
             thread, 1, MPI_GLOBAL_COMMUNICATOR, MPI_STATUS_IGNORE);
    MPI_Recv(recvInfo[thread].particleCounts.data(), count, MPI_INT,
             thread, 2, MPI_GLOBAL_COMMUNICATOR, MPI_STATUS_IGNORE);
  }

  // Send particles
  for (const auto& pair : ProcessParticleMap) {
    InitSendParticleData(pair.first, pair.second);
  }

  // Prepare receives using array indices
  for (int thread = 0; thread < PIC::nTotalThreads; thread++) {
    if (thread == PIC::ThisThread) continue;
    
    const auto& info = recvInfo[thread];
    std::vector<MPI_Aint> allOffsets;
    
    // Access segments directly using array indices
    for (size_t i = 0; i < info.segmentIndices.size(); i++) {
      FL::cFieldLineSegment* segment = ThreadSegmentTable.Table[thread][info.segmentIndices[i]];
      auto offsets = InitNewParticles(segment, info.particleCounts[i]);
      allOffsets.insert(allOffsets.end(), offsets.begin(), offsets.end());
    }
    
    if (!allOffsets.empty()) {
      MPI_Datatype particle_recv_type;
      MPI_Type_create_hindexed_block(allOffsets.size(),
                                    PIC::ParticleBuffer::ParticleDataLength - 2*sizeof(long int),
                                    allOffsets.data(),
                                    MPI_BYTE,
                                    &particle_recv_type);
      MPI_Type_commit(&particle_recv_type);

      MPI_Irecv(MPI_BOTTOM, 1, particle_recv_type, thread, 0,
                MPI_GLOBAL_COMMUNICATOR, &RecvParticleDataRequest[RecvRequestLength++]);

      MPI_Type_free(&particle_recv_type);
    }
  }

  // Wait for completion
  MPI_Waitall(SendRequestLength, SendParticleDataRequest, MPI_STATUSES_IGNORE);
  MPI_Waitall(RecvRequestLength, RecvParticleDataRequest, MPI_STATUSES_IGNORE);

  // Clean up sent particles using array indices
  for (const auto& pair : ProcessParticleMap) {
    for (const auto& desc : pair.second) {
      FL::cFieldLineSegment* segment = ThreadSegmentTable.Table[pair.first][desc.segmentIndex];
      RemoveSentParticles(segment);
    }
  }

  delete[] SendParticleDataRequest;
  delete[] RecvParticleDataRequest;
}



void PIC::ParallelFieldLines::ExchangeFieldLineParticles() {
namespace FL = PIC::FieldLine;
  if (FL::FieldLinesAll==NULL) return;

  if (ThreadSegmentTable==NULL) exit(__LINE__,__FILE__,"Error: ThreadSegmentTable is not initialized");

  //loop through field lines  
  for (int iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    ExchangeFieldLineParticles(ThreadSegmentTable[iFieldLine]);
  }
}

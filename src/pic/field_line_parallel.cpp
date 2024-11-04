
#include "pic.h"

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

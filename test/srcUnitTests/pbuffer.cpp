

#include <gtest/gtest.h>
#include "pic.h"

void pbuffer_test_for_linker() {}
 
TEST(InlineTest, IntegerEquality) {
    EXPECT_EQ(5, 5);  // Passes
    EXPECT_NE(5, 3);  // Passes
}


/*
// Test particle allocation and deletion
TEST_F(ParticleBufferTest, ParticleAllocationAndDeletionTest) {
    // Allocate new particle
    long int particlePtr = PIC::ParticleBuffer::GetNewParticle(true);
    EXPECT_GE(particlePtr, 0);
    EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 1);

    // Verify particle is marked as allocated
    EXPECT_TRUE(PIC::ParticleBuffer::IsParticleAllocated(particlePtr));

    // Delete particle
    PIC::ParticleBuffer::DeleteParticle(particlePtr);
    EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 0);
    EXPECT_FALSE(PIC::ParticleBuffer::IsParticleAllocated(particlePtr));
}

void ddd() {
	double d=2;
}
*/ 

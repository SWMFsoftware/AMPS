

#include <gtest/gtest.h>
#include "pic.h"

void idf_test_for_linker() {}
 
class IDFTest : public ::testing::Test {
protected:
   PIC::ParticleBuffer::byte* ParticleDataBuffer; 
   long int MaxNPart,NAllPart,FirstPBufferParticle;
   int *ParticleNumberTable,*ParticleOffsetTable;
   PIC::ParticleBuffer::cParticleTable *ParticlePopulationTable;

    void SetUp() override {
      //set the initial state of the particle buffer
      ParticleDataBuffer=PIC::ParticleBuffer::ParticleDataBuffer;
      MaxNPart=PIC::ParticleBuffer::MaxNPart;
      NAllPart=PIC::ParticleBuffer::NAllPart;
      FirstPBufferParticle=PIC::ParticleBuffer::FirstPBufferParticle;
      ParticleNumberTable=PIC::ParticleBuffer::ParticleNumberTable;
      ParticleOffsetTable=PIC::ParticleBuffer::ParticleOffsetTable;
      ParticlePopulationTable=PIC::ParticleBuffer::ParticlePopulationTable;


      //set the default values for the partice buffer
      PIC::ParticleBuffer::ParticleDataBuffer=NULL;
      PIC::ParticleBuffer::MaxNPart=0;
      PIC::ParticleBuffer::NAllPart=0;
      PIC::ParticleBuffer::FirstPBufferParticle=-1;
      PIC::ParticleBuffer::ParticleNumberTable=NULL;
      PIC::ParticleBuffer::ParticleOffsetTable=NULL;
      PIC::ParticleBuffer::ParticlePopulationTable=NULL; 
        	
      // Initialize buffer with 100 particles
      PIC::ParticleBuffer::Init(100);

      ASSERT_NE(nullptr, PIC::ParticleBuffer::ParticleDataBuffer)
            << "Failed to initialize ParticleDataBuffer";
      ASSERT_EQ(100, PIC::ParticleBuffer::GetMaxNPart())
            << "Failed to set MaxNPart";
    }

    void TearDown() override {
      // Cleanup
      if (PIC::ParticleBuffer::ParticleNumberTable!=NULL) amps_free_managed(PIC::ParticleBuffer::ParticleNumberTable);
      PIC::ParticleBuffer::ParticleNumberTable=ParticleNumberTable;

      if (PIC::ParticleBuffer::ParticlePopulationTable!=NULL) amps_free_managed(PIC::ParticleBuffer::ParticlePopulationTable);
      PIC::ParticleBuffer::ParticlePopulationTable=ParticlePopulationTable;

      if (PIC::ParticleBuffer::ParticleOffsetTable!=NULL) amps_free_managed(PIC::ParticleBuffer::ParticleOffsetTable);
      PIC::ParticleBuffer::ParticleOffsetTable=ParticleOffsetTable;

      if (PIC::ParticleBuffer::ParticleDataBuffer!=NULL) amps_free_managed(PIC::ParticleBuffer::ParticleDataBuffer);
      PIC::ParticleBuffer::ParticleDataBuffer=ParticleDataBuffer;

      PIC::ParticleBuffer::MaxNPart=MaxNPart;
      PIC::ParticleBuffer::NAllPart=NAllPart;
    }

};




// Test particle state vector 
TEST_F(IDFTest, IDFDataAccessTest) {
    long int ptr = PIC::ParticleBuffer::GetNewParticle(true);
    double pos[3] = {1.0, 2.0, 3.0};
    double vel[3] = {4.0, 5.0, 6.0};
    double RotE=2.5;

    auto ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

    // Verify data
    if (_PIC_FIELD_LINE_MODE_!=_PIC_MODE_ON_) {
      double* readPos = PIC::ParticleBuffer::GetX(ptr);
      double* readVel = PIC::ParticleBuffer::GetV(ptr);
      double readRotE,readVibE,VibE=RotE; 

      // Set position and velocity
      PIC::ParticleBuffer::SetX(pos, ptr);
      PIC::ParticleBuffer::SetV(vel, ptr);

      //set rotational and vibrational energy
      if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) { 
         PIC::IDF::LB::SetRotE(RotE,ParticleData); 

	 for (int nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[0];nmode++)  {
           ++VibE;
           PIC::IDF::LB::SetVibE(VibE,nmode,ParticleData);
	 } 

      }

      for(int i = 0; i < 3; i++) {
          EXPECT_DOUBLE_EQ(pos[i], readPos[i]);
          EXPECT_DOUBLE_EQ(vel[i], readVel[i]);
      }

      if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
        readRotE=PIC::IDF::LB::GetRotE(ParticleData);
	EXPECT_DOUBLE_EQ(RotE, readRotE);

        VibE=RotE;

	for (int nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[0];nmode++)  {
           ++VibE;
           readVibE=PIC::IDF::LB::GetVibE(nmode,ParticleData);
	   EXPECT_DOUBLE_EQ(VibE, readVibE);
         }
      }
      
    }
    else {
      double readPos[3],readVel[3];
      double S=4.5,readS,v_parallel,v_normal;
      double readRotE,readVibE,VibE=RotE;

      // Set position and velocity
      PIC::ParticleBuffer::SetFieldLineCoord(S,ptr); 
      PIC::ParticleBuffer::SetV(vel, ptr);

      //set rotational and vibrational energy
      if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
         PIC::IDF::LB::SetRotE(RotE,ParticleData);

	 for (int nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[0];nmode++)  {
           ++VibE;
           PIC::IDF::LB::SetVibE(VibE,nmode,ParticleData);
         }
      }

      readS=PIC::ParticleBuffer::GetFieldLineCoord(ptr);
      PIC::ParticleBuffer::GetV(readVel,ptr);

      for(int i = 0; i < 2; i++) {
        EXPECT_DOUBLE_EQ(vel[i], readVel[i]);
      }

      EXPECT_DOUBLE_EQ(S,readS);

      
      v_parallel=PIC::ParticleBuffer::GetVParallel(ptr);
      EXPECT_DOUBLE_EQ(v_parallel, readVel[0]);

      
      v_normal=PIC::ParticleBuffer::GetVNormal(ptr);
      EXPECT_DOUBLE_EQ(v_normal, readVel[1]); 

      EXPECT_DOUBLE_EQ(readVel[2],0.0);

      if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
        readRotE=PIC::IDF::LB::GetRotE(ParticleData);
        EXPECT_DOUBLE_EQ(RotE, readRotE);

	VibE=RotE;

	for (int nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[0];nmode++)  {
          ++VibE;
          readVibE=PIC::IDF::LB::GetVibE(nmode,ParticleData);
	  EXPECT_DOUBLE_EQ(VibE, readVibE);
        }
      }
    }
    
}



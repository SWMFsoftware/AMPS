

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

    EXPECT_NE(PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_,-1);
    EXPECT_NE(PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_,-1);

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


// Test particle state vector
TEST_F(IDFTest, IDFRotEnregyGenerationTest) {
  int ntest;
  double RotEtheory,SampledRotE=0.0;
  long int ptr = PIC::ParticleBuffer::GetNewParticle(true);
  auto ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

  double Temp=300.0;
  int nTotalTests=1000000;

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    PIC::ParticleBuffer::SetI(s,ptr);
    SampledRotE=0.0;

    for (ntest=0;ntest<nTotalTests;ntest++) {
      PIC::IDF::LB::InitRotTemp(Temp,ParticleData);
      SampledRotE+=PIC::IDF::LB::GetRotE(ParticleData);
    }

    SampledRotE/=nTotalTests;
    RotEtheory=PIC::IDF::nTotalRotationalModes[s]/2.0*Kbol*Temp;

    if (PIC::IDF::nTotalRotationalModes[s]==0) {
      EXPECT_LT(fabs(SampledRotE),0.01) << "RotE distribution: s=" << s << endl;
    }
    else {
      EXPECT_LT(fabs(SampledRotE-RotEtheory)/(SampledRotE+RotEtheory),0.01) << "RotE distribution: s=" << s << ", Rotational Modes=" << PIC::IDF::nTotalRotationalModes[s] << endl;
    }
  }
}

TEST_F(IDFTest, IDFVibEnregyGenerationTest) {
  int nmode,ntest;
  double VibEtheory,SampledVibE[10];
  long int ptr = PIC::ParticleBuffer::GetNewParticle(true);
  auto ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

  double Temp=300.0;
  int nTotalTests=200000;

  for (int iTempTest=0;iTempTest<3;iTempTest++) {
    for (int s=0;s<PIC::nTotalSpecies;s++) {
      Temp=0.5*(1+iTempTest)*PIC::IDF::CharacteristicVibrationalTemperature[s];
      PIC::ParticleBuffer::SetI(s,ptr);
      for (nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[s];nmode++) SampledVibE[nmode]=0.0;

      for (ntest=0;ntest<nTotalTests;ntest++) {
        PIC::IDF::LB::InitVibTemp(Temp,ParticleData);

        for (nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[s];nmode++) {
          SampledVibE[nmode]+=PIC::IDF::LB::GetVibE(nmode,ParticleData);
        }
      }

      for (nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[s];nmode++) {
        double ThetaVib,c,EtaVib; 

	ThetaVib=PIC::IDF::CharacteristicVibrationalTemperature[nmode+s*PIC::IDF::nSpeciesMaxVibrationalModes]; 
	c=ThetaVib/Temp;
        EtaVib=2.0*c/(exp(c)-1.0);

        SampledVibE[nmode]/=nTotalTests;
        VibEtheory=EtaVib/2.0*Kbol*Temp;

        EXPECT_LT(fabs(SampledVibE[nmode]-VibEtheory)/(SampledVibE[nmode]+VibEtheory),0.01) << "VibE distribution: s=" << s << ", nmode=" << nmode << ", Temp=" << Temp <<  endl;
      }
    }
  }
}

TEST_F(IDFTest, IDFEnergyDistributionTest) {
  int nmode,ntest;
  double VibEtheory,SampledVibE[10];
  long int ptr = PIC::ParticleBuffer::GetNewParticle(true);
  auto ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

  double Temp=300.0;
  int nTotalTests=200000;
  double CollisionEnergy_init=6.0*Kbol*Temp;
  double CollisionEnergy,InternalEnergy,InternalEnergy_sum=0.0,TranslationalEnergy,TranslationalEnergy_sum=0.0;

  for (int nInternalDF=1;nInternalDF<=3;nInternalDF++) {
    for (int nTranslationalDF=1;nTranslationalDF<=3;nTranslationalDF++) {
      for (int iTest=0;iTest<nTotalTests;iTest++) {
        CollisionEnergy=CollisionEnergy_init;
        PIC::IDF::LB::DistributeEnergy(InternalEnergy,CollisionEnergy,nInternalDF,nTranslationalDF);

        TranslationalEnergy_sum+=TranslationalEnergy;
        InternalEnergy_sum+=InternalEnergy;
      }

      TranslationalEnergy_sum/=nTotalTests*nTranslationalDF;
      InternalEnergy_sum/=nTotalTests*nInternalDF;

      EXPECT_LT(fabs(TranslationalEnergy_sum-InternalEnergy_sum)/(TranslationalEnergy_sum+InternalEnergy_sum),0.01) << 
      "EnergyDistributionTest: nInternalDF=" << nInternalDF << ", nTranslationalDF=" << nTranslationalDF <<  endl;
    }

  }
}

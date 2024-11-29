#include <gtest/gtest.h>
#include "pic.h"

void collisions_test_for_linker() {}


namespace AMPS_COLLISON_TEST  {
  double mv2sum,mv2sumAfter,mvSumAfter[3];
  double CollsionFrequentcyTheory[2][2];
  double CollsionFrequentcy[2][2];
  double MeanRotEnergy[2],MeanRotEnergyAfter[2],MeanRotEnergyTheory[2];
  double MeanVibEnergy[2],MeanVibEnergyAfter[2],MeanVibEnergyTheory[2];

  void HeatBath(int s0,int s1,void (*fCellCollision)(int, int, int, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)) {
    namespace MD=PIC::MolecularData;
    namespace MC=PIC::MolecularCollisions;
    namespace PB=PIC::ParticleBuffer; 

    int i,j,icell=0,jcell=0,kcell=0;
    long int ptr;
    PIC::Mesh::cDataCenterNode *cell;
    PIC::Mesh::cDataBlockAMR *block;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];    

    int CollFreqOffset=MC::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset+    
      sizeof(double)*MC::ParticleCollisionModel::CollsionFrequentcySampling::Offset(0,0);

    for (i=0;i<2;i++) {
      MeanVibEnergy[i]=0.0,MeanVibEnergyAfter[i]=0.0,MeanVibEnergyTheory[i]=0.0;
      MeanRotEnergy[i]=0.0,MeanRotEnergyAfter[i]=0.0,MeanRotEnergyTheory[i]=0.0;
    }

    //1. Identify the cell
    block=node->block;
    if (block==NULL) exit(__LINE__,__FILE__,"Block is NULL");

    cell=block->GetCenterNode(_getCenterNodeLocalNumber(icell,jcell,kcell));
    if (cell==NULL) exit(__LINE__,__FILE__,"Cell is NULL");

    char *SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;
  
    //2. Populate the particle lists
    double Temp=300.0; //the temperature of the gas
    int nTotalInjectedParticles=50;
    double mmax,v[3],vbulk[3]={0.0,0.0,0.0},mvsum[3]={0.0,0.0,0.0},msum=0.0;

    mv2sum=0.0;

    if (s0==s1) {
      CollFreqOffset=MC::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset+
        sizeof(double)*MC::ParticleCollisionModel::CollsionFrequentcySampling::Offset(s0,s0);

      *((double*)(SamplingData+CollFreqOffset))=0.0;
    }
    else {
      int map[4][2]={{s0,s0},{s0,s1},{s1,s0},{s1,s1}};

      for (i=0;i<4;i++) { 
        CollFreqOffset=MC::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset+
          sizeof(double)*MC::ParticleCollisionModel::CollsionFrequentcySampling::Offset(map[i][0],map[i][1]);

	*((double*)(SamplingData+CollFreqOffset))=0.0;
      }
    }

    mmax=max(MD::GetMass(s0),MD::GetMass(s1));
 
    for (int ispec=0;ispec<2;ispec++) {
      int s=(ispec==0) ? s0 : s1;
      double m=MD::GetMass(s)/mmax; 

      for (int i=0;i<nTotalInjectedParticles;i++) {
        ptr=PB::GetNewParticle(block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)]);
        PIC::Distribution::MaxwellianVelocityDistribution(v,vbulk,Temp,s);

        PB::SetV(v,ptr);
        PB::SetI(s,ptr);

        if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
          PB::SetIndividualStatWeightCorrection(1.0,ptr);
        }

	if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
          PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr); 

          switch (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_) {
          case _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_:
            PIC::IDF::LB::InitRotTemp(Temp,ParticleData); 
            PIC::IDF::LB::InitVibTemp(Temp,ParticleData); 
	    break;

          default:
            exit(__LINE__,__FILE__,"Error: not implemented");
	  }
        }

        for (int j=0;j<3;j++) mvsum[j]+=m*v[j];
        msum+=m;
      }
    }

    //3. Correct builk velocity
    double dv[3];

    for (int j=0;j<3;j++) dv[j]=mvsum[j]/msum;

    for (int i=0;i<2;i++) MeanVibEnergy[i]=0.0,MeanVibEnergyAfter[i]=0.0,MeanRotEnergy[i]=0.0,MeanRotEnergyAfter[i]=0.0;

    ptr=block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)];

    while (ptr!=-1) {
      double *v=PB::GetV(ptr);
      int s=PB::GetI(ptr);
      double m=MD::GetMass(s)/mmax;
      int ispec=(s==s0) ? 0 : 1;


      for (int j=0;j<3;j++) {
        v[j]-=dv[j];
        mv2sum+=m*v[j]*v[j];
      }

      if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
        PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

        switch (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_) {
        case _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_:
          MeanVibEnergy[ispec]+=PIC::IDF::LB::GetVibE(-1,ParticleData);
          MeanRotEnergy[ispec]+=PIC::IDF::LB::GetRotE(ParticleData);
          break;
        default:
          exit(__LINE__,__FILE__,"Error: not implemented");
	}	  
      }

      ptr=PB::GetNext(ptr);
    }

    //Calculate initial rotational and vibrational evergy
    if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
      if (s0==s1) {
        MeanVibEnergy[0]/=2.0*nTotalInjectedParticles; 
        MeanRotEnergy[0]/=2.0*nTotalInjectedParticles; 
      }
      else {
        for (int i=0;i<2;i++) {
          MeanVibEnergy[i]/=nTotalInjectedParticles;
          MeanRotEnergy[i]/=nTotalInjectedParticles;
        }
      }
    }


    //4. Simulate collisions: deermine the particle weight such that for dt=1, there will be 1 collision per particle
    double dt=1.0,Measure=1.0;
    int nIterations=100000;

    //set the weight, time step, and measure
    double InitWeight[2],InitMeasure,InitTimeStep[2];
    double m_s0,m_s1,m,vv[3]={0.0,0.0,0.0};
    double MeanRelativeVelocity,sigma,eanRelativeVelocity,StatWeight[2],n[2];
    int ispec,s;

    m_s0=MD::GetMass(s0);
    m_s1=MD::GetMass(s1); 
    m=m_s0*m_s1/(m_s0+m_s1);
    MeanRelativeVelocity=sqrt(8.0*Kbol*Temp/(Pi*m));

    InitMeasure=cell->Measure;
    cell->Measure=Measure;

    for (ispec=0;ispec<2;ispec++) {
      s=(ispec==0) ? s0 : s1;

      sigma=MD::MolecularModels::GetTotalCrossSection(s,vv,s,vv);
      m=MD::GetMass(s);

      MeanRelativeVelocity=sqrt(16.0*Kbol*Temp/(Pi*m));
      StatWeight[ispec]=Measure/(2.0*nTotalInjectedParticles*sigma*MeanRelativeVelocity*dt);

      if (s0!=s1) {
        InitTimeStep[ispec]=block->GetLocalTimeStep(s);
        InitWeight[ispec]=block->GetLocalParticleWeight(s);
      }
      else {
        InitTimeStep[ispec]=(ispec==0) ? block->GetLocalTimeStep(s) : InitTimeStep[0];
        InitWeight[ispec]=(ispec==0) ? block->GetLocalParticleWeight(s) : InitWeight[0];
      }

      block->SetLocalTimeStep(dt,s);
      block->SetLocalParticleWeight(StatWeight[ispec],s);

      n[ispec]=StatWeight[ispec]*nTotalInjectedParticles/Measure; 
    }

    for (i=0;i<((s0==s1) ? 1 : 2);i++) for (j=0;j<((s0==s1) ? 1 : 2);j++) {
      if (i==j) {
        MeanRelativeVelocity=sqrt(16.0*Kbol*Temp/(Pi*m_s0));
        sigma=MD::MolecularModels::GetTotalCrossSection(s0,vv,s0,vv); 
          
        CollsionFrequentcyTheory[i][j]=pow(n[0]+n[1],2)*sigma*MeanRelativeVelocity;
      }
      else {
        m=m_s0*m_s1/(m_s0+m_s1);
        MeanRelativeVelocity=sqrt(8.0*Kbol*Temp/(Pi*m));
        sigma=MD::MolecularModels::GetTotalCrossSection(s0,vv,s1,vv); 
        CollsionFrequentcyTheory[i][j]=n[0]*n[1]*sigma*MeanRelativeVelocity;
      }   
    }


     if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
        for (ispec=0;ispec<2;ispec++) {
          int n,s=(ispec==0) ? s0 : s1;
          double RotDF,EtaVib,c; 

          switch (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_) {
          case _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_:
            RotDF=PIC::IDF::nTotalRotationalModes[s];
	    MeanRotEnergyTheory[ispec]=0.5*RotDF*Kbol*Temp;

	    MeanVibEnergyTheory[ispec]=0.0;

            for (n=0;n<PIC::IDF::nTotalVibtationalModes[s];n++) {
              double c,ThetaVib,EtaVib;
	    
	      ThetaVib=PIC::IDF::CharacteristicVibrationalTemperature[n+s*PIC::IDF::nSpeciesMaxVibrationalModes];
	      c=ThetaVib/Temp;
	      EtaVib=2.0*c/(exp(c)-1.0);
	      MeanVibEnergyTheory[ispec]+=0.5*EtaVib*Kbol*Temp;
	    }

            break;
          default:
            exit(__LINE__,__FILE__,"Error: not implemented");
          }
        }
      }


    for (int i=0;i<nIterations;i++) {
      fCellCollision(icell,jcell,kcell,node);
    }

    //5. Get the collision frequency
    CollsionFrequentcy[0][0]=(*SamplingData)/nIterations;

    if (s0==s1) {
      CollFreqOffset=MC::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset+
        sizeof(double)*MC::ParticleCollisionModel::CollsionFrequentcySampling::Offset(s0,s0);

      CollsionFrequentcy[0][0]=(*((double*)(SamplingData+CollFreqOffset)))/nIterations;
    }
    else {
      int map[4][2]={{s0,s0},{s0,s1},{s1,s0},{s1,s1}};

      for (i=0;i<4;i++) {
        CollFreqOffset=MC::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset+
          sizeof(double)*MC::ParticleCollisionModel::CollsionFrequentcySampling::Offset(map[i][0],map[i][1]);

        CollsionFrequentcy[map[i][0]][map[i][1]]=(*((double*)(SamplingData+CollFreqOffset)))/nIterations; 
      }
    }



    //6. Verify the results
    mv2sumAfter=0.0;
    for (int idim=0;idim<3;idim++) mvSumAfter[idim]=0.0;

    ptr=block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)];

    while (ptr!=-1) {
      double *v=PB::GetV(ptr);
      int s=PB::GetI(ptr);
      double m=MD::GetMass(s)/mmax;
      int ispec=(s==s0) ? 0 : 1;

      for (int j=0;j<3;j++) {
        mvSumAfter[j]+=m*v[j];
        mv2sumAfter+=m*v[j]*v[j];
      }

      if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
        PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

        switch (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_) {
        case _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_:
          MeanVibEnergyAfter[ispec]+=PIC::IDF::LB::GetVibE(-1,ParticleData);
          MeanRotEnergyAfter[ispec]+=PIC::IDF::LB::GetRotE(ParticleData);
          break;
        default:
          exit(__LINE__,__FILE__,"Error: not implemented");
        }
      }

      long int ptr_next=PB::GetNext(ptr);
      PB::DeleteParticle(ptr);
      ptr=ptr_next;
    }

    //Calculate final rotational and vibrational evergy
    if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
      if (s0==s1) {
        MeanVibEnergyAfter[0]/=2.0*nTotalInjectedParticles;
        MeanRotEnergyAfter[0]/=2.0*nTotalInjectedParticles;
      }
      else {
        for (int i=0;i<2;i++) {
          MeanVibEnergyAfter[i]/=nTotalInjectedParticles;
          MeanRotEnergyAfter[i]/=nTotalInjectedParticles;
        }
      }
    }

    //7. Restore the initial state
    for (ispec=0;ispec<2;ispec++) {
      s=(ispec==0) ? s0 : s1;

      block->SetLocalTimeStep(InitTimeStep[ispec],s); 
      block->SetLocalParticleWeight(InitWeight[ispec],s);
    }


    block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)]=-1;
  }
}

struct ParticleCollisionTestCase {
    void (*fCellCollision)(int, int, int, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*); // Function pointer
    int s0,s1;
    string name;
};


class ParticleCollisionTest : public ::testing::TestWithParam<ParticleCollisionTestCase> {
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
        	

      // Initialize buffer with 200 particles
      PIC::ParticleBuffer::Init(200);

      ASSERT_NE(nullptr, PIC::ParticleBuffer::ParticleDataBuffer)
            << "Failed to initialize ParticleDataBuffer";
      ASSERT_EQ(200, PIC::ParticleBuffer::GetMaxNPart())
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


//Test NTC collision model: tested conservations of momentum and energy, and colliaiio frequency
TEST_P(ParticleCollisionTest, MyHandlesInputs) {
  using namespace AMPS_COLLISON_TEST;
  int s0,s1;

  // Print the function name
  ParticleCollisionTestCase test_case=GetParam();
  s0=test_case.s0;
  s1=test_case.s1;

  std::cout << "\033[1m" << "Testing function: " << test_case.name << " (s0=" << s0 << ", s1=" << s1 <<") \033[0m" << std::endl;

  HeatBath(s0,s1,test_case.fCellCollision);

  if (s0==s1) {
    EXPECT_LT(fabs(mv2sum-mv2sumAfter)/(mv2sum+mv2sumAfter),1.0E-10); 
    EXPECT_LT(fabs(CollsionFrequentcy[0][0]-CollsionFrequentcyTheory[0][0])/(CollsionFrequentcy[0][0]+CollsionFrequentcyTheory[0][0]),1.0E-5) << endl;

    if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
      if (MeanVibEnergyTheory[0]+MeanVibEnergy[0]!=0.0) 
        EXPECT_LT(fabs(MeanVibEnergyTheory[0]-MeanVibEnergy[0])/(MeanVibEnergyTheory[0]+MeanVibEnergy[0]),1.0E-5);

      if (MeanRotEnergyTheory[0]+MeanRotEnergy[0]!=0.0)
        EXPECT_LT(fabs(MeanRotEnergyTheory[0]-MeanRotEnergy[0])/(MeanRotEnergyTheory[0]+MeanRotEnergy[0]),1.0E-5); 

      if (MeanVibEnergyAfter[0]+MeanVibEnergy[0]!=0.0)
        EXPECT_LT(fabs(MeanVibEnergyAfter[0]-MeanVibEnergy[0])/(MeanVibEnergyAfter[0]+MeanVibEnergy[0]),1.0E-5);

      if (MeanRotEnergyAfter[0]+MeanRotEnergy[0]!=0.0)
        EXPECT_LT(fabs(MeanRotEnergyAfter[0]-MeanRotEnergy[0])/(MeanRotEnergyAfter[0]+MeanRotEnergy[0]),1.0E-5);
    }
  }
  else {
    EXPECT_LT(fabs(mv2sum-mv2sumAfter)/(mv2sum+mv2sumAfter),1.0E-10) << test_case.name << ", s0=" << s0 << ", s1=" << s1 << endl;

    for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
      EXPECT_LT(fabs(CollsionFrequentcy[i][j]-CollsionFrequentcyTheory[i][j])/(CollsionFrequentcy[i][j]+CollsionFrequentcyTheory[i][j]),1.0E-5) << " i=" << i << ", j=" << j << endl;
    }

    if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_==_PIC_MODE_ON_) {
      for (int i=0;i<2;i++) {
        if (MeanVibEnergyAfter[i]+MeanVibEnergy[i]!=0.0)
          EXPECT_LT(fabs(MeanVibEnergyAfter[i]-MeanVibEnergy[i])/(MeanVibEnergyAfter[i]+MeanVibEnergy[i]),1.0E-5) << " i=" << i << endl;

        if (MeanRotEnergyAfter[i]+MeanRotEnergy[i]!=0.0)
          EXPECT_LT(fabs(MeanRotEnergyAfter[i]-MeanRotEnergy[i])/(MeanRotEnergyAfter[i]+MeanRotEnergy[i]),1.0E-5) << " i=" << i << endl;

        if (MeanVibEnergyTheory[i]+MeanVibEnergy[i]!=0.0)
          EXPECT_LT(fabs(MeanVibEnergyTheory[i]-MeanVibEnergy[i])/(MeanVibEnergyTheory[i]+MeanVibEnergy[i]),1.0E-5) << " i=" << i << endl;

        if (MeanRotEnergyTheory[i]+MeanRotEnergy[i]!=0.0)
          EXPECT_LT(fabs(MeanRotEnergyTheory[i]-MeanRotEnergy[i])/(MeanRotEnergyTheory[i]+MeanRotEnergy[i]),1.0E-5) << " i=" << i << endl;
      }
    } 
  }

  for (int jj=0;jj<3;jj++) EXPECT_LT(fabs(mvSumAfter[jj])/sqrt(mv2sumAfter),1.0E-10) << ", jj=" << jj <<  endl;
}




INSTANTIATE_TEST_SUITE_P(
    ParticleCollisionTest,             // Test suite name
    ParticleCollisionTest,             // Test fixture name
    ::testing::Values(                 // Test cases
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,0,0,"mf_Yinsi"}, 
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,0,1,"mf_Yinsi"},
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,1,0,"mf_Yinsi"},
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,1,1,"mf_Yinsi"},

        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,0,0,"mf"}, 
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,0,1,"mf"},
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,1,0,"mf"},
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,1,1,"mf"},

        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,0,0,"ntc"}, 
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,0,1,"ntc"},
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,1,0,"ntc"},
        ParticleCollisionTestCase{PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,1,1,"ntc"}
      )
);

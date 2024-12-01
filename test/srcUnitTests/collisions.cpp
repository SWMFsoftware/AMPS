#include <gtest/gtest.h>
#include "pic.h"

#include "ParticleTestBase.h"

void collisions_test_for_linker() {}

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
namespace AMPS_COLLISON_TEST  {
  double mv2sum,mv2sumAfter,mvSumAfter[3];
  double CollsionFrequentcyTheory[2][2];
  double CollsionFrequentcy[2][2];
  double MeanRotEnergy[2],MeanRotEnergyAfter[2],MeanRotEnergyTheory[2];
  double MeanVibEnergy[2],MeanVibEnergyAfter[2],MeanVibEnergyTheory[2];

  void  VelocityAfterCollision(double &RelativeErrorMomentum,double &RelativeErrorEnergy,double &RelativeErrorVel,int s0,int s1,void (*fVelocityAfterCollision)(double*,double,double*,double)) {
    double v_s0[3],v_s1[3],m_s0,m_s1,am,vrel[3],vbulk[3]={0.0,0.0,0.0},mv2sum_init,mv2sum_after,mvsum_init,mvsum_after,vcm[3];
    double vrel_init,vrel_after;
    int nTest,i;

    const int nTotalTests=100000;
    const double Temp=300.0;

    m_s0=PIC::MolecularData::GetMass(s0);
    m_s1=PIC::MolecularData::GetMass(s1);

    RelativeErrorMomentum=0.0;
    RelativeErrorEnergy=0.0;
    RelativeErrorVel=0.0;

    for (nTest=0;nTest<nTotalTests;nTest++) {
      //generate particle velocity and calculate relative speed
      PIC::Distribution::MaxwellianVelocityDistribution(v_s0,vbulk,Temp,s0);
      PIC::Distribution::MaxwellianVelocityDistribution(v_s1,vbulk,Temp,s1);

      //calculating initial parameters
      mv2sum_init=0.0,mvsum_init=0.0;
      for (i=0;i<3;i++) {
        double t;
        
        t=m_s0*v_s0[i];
        mv2sum_init+=t*v_s0[i];
        mvsum_init+=t;

        t=m_s1*v_s1[i];
        mv2sum_init+=t*v_s1[i];
        mvsum_init+=t; 

        vrel[i]=v_s0[i]-v_s1[i];
      }

      vrel_init=Vector3D::Length(vrel);

      //redistribute the particle velocity
      fVelocityAfterCollision(v_s0,m_s0,v_s1,m_s1);

      //Verify the result 
      mv2sum_after=0.0,mvsum_after=0.0;
      for (i=0;i<3;i++) {
        double t;
        
        t=m_s0*v_s0[i];
        mv2sum_after+=t*v_s0[i];
        mvsum_after+=t;

        t=m_s1*v_s1[i];
        mv2sum_after+=t*v_s1[i];
        mvsum_after+=t; 

        vrel[i]=v_s0[i]-v_s1[i];
      }   

      vrel_after=Vector3D::Length(vrel);
      RelativeErrorVel+=fabs(vrel_init-vrel_after)/(vrel_init+vrel_after);

      RelativeErrorEnergy+=fabs(mv2sum_init-mv2sum_after)/(mv2sum_init+mv2sum_after);  
      RelativeErrorMomentum+=fabs(mvsum_init-mvsum_after)/(mvsum_init+mvsum_after);    
    }

  
  }

  void HeatBath(int s0,int s1,bool UseCommonStatWeight,bool UseCommonTimeStep,void (*fCellCollision)(int, int, int, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)) {
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
    int nTotalInjectedParticles=100;
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
    double dt[2]={1.0,1.0},Measure=1.0;
    int nIterations=400000;

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
      StatWeight[ispec]=Measure/(2.0*nTotalInjectedParticles*sigma*MeanRelativeVelocity*dt[ispec]);

      if (s0!=s1) {
        InitTimeStep[ispec]=block->GetLocalTimeStep(s);
        InitWeight[ispec]=block->GetLocalParticleWeight(s);
      }
      else {
        InitTimeStep[ispec]=(ispec==0) ? block->GetLocalTimeStep(s) : InitTimeStep[0];
        InitWeight[ispec]=(ispec==0) ? block->GetLocalParticleWeight(s) : InitWeight[0];
      }

      if (ispec==1) {
        if (UseCommonTimeStep==true) {
          dt[1]=dt[0];
	} else if (dt[0]==dt[1]) {
          dt[1]=0.33*dt[0];
	}

	if (UseCommonStatWeight==true) {
          StatWeight[1]=StatWeight[0];
	}
	else if (StatWeight[0]==StatWeight[1]) {
          StatWeight[1]=0.33*StatWeight[0];
	}
      }

      block->SetLocalTimeStep(dt[ispec],s);
      block->SetLocalParticleWeight(StatWeight[ispec],s);

      n[ispec]=StatWeight[ispec]*nTotalInjectedParticles/Measure; 
    }

    if (s0==s1) {
      MeanRelativeVelocity=sqrt(16.0*Kbol*Temp/(Pi*m_s0));
      sigma=MD::MolecularModels::GetTotalCrossSection(s0,vv,s0,vv);

      CollsionFrequentcyTheory[0][0]=pow(n[0]+n[1],2)*sigma*MeanRelativeVelocity;     
      CollsionFrequentcyTheory[1][1]=CollsionFrequentcyTheory[0][0];
    }
    else {
      for (i=0;i<2;i++) for (j=0;j<2;j++) {
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


// Derived class for particle collision tests
class ParticleCollisionTest :
    public ParticleTestBase<>,
    public ::testing::TestWithParam<ParticleCollisionTestCase> {
protected:
    void SetUp() override {
        ParticleTestBase<>::SetUp();
    }

    void TearDown() override {
        ParticleTestBase<>::TearDown();
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

  HeatBath(s0,s1,true,true,test_case.fCellCollision);

  if (s0==s1) {
    EXPECT_LT(fabs(mv2sum-mv2sumAfter)/(mv2sum+mv2sumAfter),1.0E-10); 
    EXPECT_LT(fabs(CollsionFrequentcy[0][0]-CollsionFrequentcyTheory[0][0])/(CollsionFrequentcy[0][0]+CollsionFrequentcyTheory[0][0]),1.0E-2) << endl;

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
      EXPECT_LT(fabs(CollsionFrequentcy[i][j]-CollsionFrequentcyTheory[i][j])/(CollsionFrequentcy[i][j]+CollsionFrequentcyTheory[i][j]),1.0E-2) << " i=" << i << ", j=" << j << endl;
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

//particle velocity redistribution test 
struct VelocityAfterCollisionTestCase {
  void (*fVelocityAfterCollision)(double*,double,double*,double); // Function pointer
  int s0,s1;
  string name;
};

// Derived class for velocity after collision tests
class VelocityAfterCollisionTest :
    public ParticleTestBase<>,
    public ::testing::TestWithParam<VelocityAfterCollisionTestCase> {
protected:
    void SetUp() override {
        ParticleTestBase<>::SetUp();
    }

    void TearDown() override {
        ParticleTestBase<>::TearDown();
    }
};


TEST_P(VelocityAfterCollisionTest, MyHandlesInputs) {
  using namespace AMPS_COLLISON_TEST;
  double RelativeErrorEnergy,RelativeErrorMomentum,RelativeErrorVel;
  int s0,s1;

  // Print the function name
  VelocityAfterCollisionTestCase test_case=GetParam();
  s0=test_case.s0;
  s1=test_case.s1;

  std::cout << "\033[1m" << "Testing function: " << test_case.name << " (s0=" << s0 << ", s1=" << s1 <<") \033[0m" << std::endl;

  VelocityAfterCollision(RelativeErrorMomentum,RelativeErrorEnergy,RelativeErrorVel,s0,s1,test_case.fVelocityAfterCollision);
  EXPECT_LT(RelativeErrorMomentum,1.0E-10);
  EXPECT_LT(RelativeErrorEnergy,1.0E-10);
  EXPECT_LT(RelativeErrorVel,1.0E-10);
}

INSTANTIATE_TEST_SUITE_P(
    VelocityAfterCollisionTest,        // Test suite name
    VelocityAfterCollisionTest,        // Test fixture name
    ::testing::Values(                 // Test cases
        VelocityAfterCollisionTestCase{PIC::MolecularCollisions::VelocityScattering::HS::VelocityAfterCollision,0,0,"VelocityScattering::HS"}, 
        VelocityAfterCollisionTestCase{PIC::MolecularCollisions::VelocityScattering::HS::VelocityAfterCollision,0,1,"VelocityScattering::HS"},
        VelocityAfterCollisionTestCase{PIC::MolecularCollisions::VelocityScattering::HS::VelocityAfterCollision,1,0,"VelocityScattering::HS"},
        VelocityAfterCollisionTestCase{PIC::MolecularCollisions::VelocityScattering::HS::VelocityAfterCollision,1,1,"VelocityScattering::HS"}
      )
);

#endif
  

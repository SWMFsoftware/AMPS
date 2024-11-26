#include <gtest/gtest.h>
#include "pic.h"

void collisions_test_for_linker() {}

class ParticleCollisionTest : public ::testing::Test {
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
TEST_F(ParticleCollisionTest, CollisionTestNTC) {
  int i,icell=0,jcell=0,kcell=0;
  long int ptr;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];
  
  //1. Identify the cell
  block=node->block;
  if (block==NULL) exit(__LINE__,__FILE__,"Block is NULL");

  cell=block->GetCenterNode(_getCenterNodeLocalNumber(icell,jcell,kcell));
  if (cell==NULL) exit(__LINE__,__FILE__,"Cell is NULL");
  
  //2. Populate the particle lists
  double Temp=300.0; //the temperature of the gas
  int nTotalInjectedParticles=100;
  double v2sum=0.0,v[3],vbulk[3]={0.0,0.0,0.0},vsum[3]={0.0,0.0,0.0};

  for (int i=0;i<nTotalInjectedParticles;i++) {
    ptr=PIC::ParticleBuffer::GetNewParticle(block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)]);
    PIC::Distribution::MaxwellianVelocityDistribution(v,vbulk,Temp,0);

    PIC::ParticleBuffer::SetV(v,ptr);
    PIC::ParticleBuffer::SetI(0,ptr);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,ptr);

    for (int j=0;j<3;j++) vsum[j]+=v[j];
  }

  //3. Correct builk velocity
  double dv[3];

  for (int j=0;j<3;j++) dv[j]=vsum[j]/nTotalInjectedParticles;
 
  ptr=block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)];

  while (ptr!=-1) {
    double *v=PIC::ParticleBuffer::GetV(ptr);

    for (int j=0;j<3;j++) {
      v[j]-=dv[j];
      v2sum+=v[j]*v[j];
    }

    ptr=PIC::ParticleBuffer::GetNext(ptr);
  }

  //4. Simulate collisions: deermine the particle weight such that for dt=1, there will be 1 collision per particle
  double dt=1.0,Measure=1.0;
  int nIterations=1000;

  double m=PIC::MolecularData::GetMass(0),vv[3]={0.0,0.0,0.0};
  double sigma=PIC::MolecularData::MolecularModels::GetTotalCrossSection(0,vv,0,vv);
  double MeanRelativeVelocity=sqrt(16.0*Kbol*Temp/Pi*m);
  double StatWeight=Measure/(nTotalInjectedParticles*sigma*MeanRelativeVelocity*dt); 


  //set the weight, time step, and measure
  double InitWeight,InitMeasure,InitTimeStep;

  InitMeasure=cell->Measure;
  cell->Measure=Measure;

  InitTimeStep=block->GetLocalTimeStep(0);
  block->SetLocalTimeStep(dt,0);

  InitWeight=block->GetLocalParticleWeight(0);
  block->SetLocalParticleWeight(StatWeight,0);

  for (int i=0;i<nIterations;i++) {
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc(icell,jcell,kcell,node);
  }


  //5. Verify the results
  double v2sumAfter=0.0,vSumAfter[3]={0.0,0.0,0.0};

  ptr=block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)];

  while (ptr!=-1) {
    double *v=PIC::ParticleBuffer::GetV(ptr);

    for (int j=0;j<3;j++) {
      vSumAfter[j]+=v[j];
      v2sumAfter+=v[j]*v[j];
    }

    long int ptr_next=PIC::ParticleBuffer::GetNext(ptr);
    PIC::ParticleBuffer::DeleteParticle(ptr);
    ptr=ptr_next;
  }

  EXPECT_DOUBLE_EQ(v2sum,v2sumAfter);
  for (int j=0;j<3;j++) EXPECT_DOUBLE_EQ(0.0,vSumAfter[j]);


  //6. Restore the initial state
  block->SetLocalTimeStep(InitTimeStep,0);
  block->SetLocalParticleWeight(InitWeight,0);
  block->FirstCellParticleTable[icell+_BLOCK_CELLS_X_*(jcell+_BLOCK_CELLS_Y_*kcell)]=-1;
}




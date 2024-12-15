#include <gtest/gtest.h>
#include "pic.h"

#include "ParticleTestBase.h"

void distribution_test_for_linker() {}

// Derived class for particle distribution tests
class ParticleDistributionTest :
    public ParticleTestBase<10000>,
    public ::testing::Test {
protected:
    void SetUp() override {
        ParticleTestBase<10000>::SetUp();
    }

    void TearDown() override {
        ParticleTestBase<10000>::TearDown();
    }
};



TEST_F(ParticleDistributionTest, MexwellianDistributionTest) {
  double v[3],speed_sum=0.0,vel_sum[3]={0.0},vbulk[3]={0.0};
  int spec,iTest;

  const int nTests=100000;
  const double Temp=300.0;

  for (iTest=0;iTest<nTests;iTest++) {
    PIC::Distribution::MaxwellianVelocityDistribution(v,vbulk,Temp,0);

    speed_sum+=Vector3D::Length(v);
    for (int i=0;i<3;i++) vel_sum[i]+=v[i];
  } 

  //verify the results 
  speed_sum/=nTests; 

  for (int i=0;i<3;i++) {
    vel_sum[i]/=nTests;

    EXPECT_LT(fabs(vel_sum[i])/speed_sum,1.0E-2); 
  }

  double SpeedTheory,m;

  m=PIC::MolecularData::GetMass(0);
  SpeedTheory=sqrt(8.0*Kbol*Temp/(Pi*m));

  EXPECT_LT(fabs(speed_sum-SpeedTheory)/(speed_sum+SpeedTheory),1.0E-2);
}



#include "test_harness.h"
#include "pic_deposition_kernels.h"

#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

AMPS_TEST_CASE("PIC-DEPOSITION-ANALYTIC-001",
               "one CIC particle reproduces independent tensor-product corner weights") {
  const double xi[3]={0.25,0.5,0.75};
  double weight[8];
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::TrilinearWeights(xi,weight)==
                   PIC::DepositionKernels::Status::Success);
  // Independent products of {3/4,1/4} x {1/2,1/2} x {1/4,3/4}.
  const double expected[8]={3.0/32.0,1.0/32.0,3.0/32.0,1.0/32.0,
                            9.0/32.0,3.0/32.0,9.0/32.0,3.0/32.0};
  for (int c=0;c<8;++c) AMPS_EXPECT_NEAR(weight[c],expected[c],0.0);
}

AMPS_TEST_CASE("PIC-DEPOSITION-TOTALS-001",
               "CIC deposition preserves total represented charge and q times velocity") {
  const double xi[3]={0.17,0.61,0.43};
  const double q=-7.25;
  const double velocity[3]={120.0,-35.0,0.125};
  double charge[8],qv[8][3];
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::DepositChargeAndVelocityCurrent(
      xi,q,velocity,charge,qv)==PIC::DepositionKernels::Status::Success);
  double charge_sum=0.0,current_sum[3]={0.0,0.0,0.0};
  for (int c=0;c<8;++c) {
    charge_sum+=charge[c];
    for (int d=0;d<3;++d) current_sum[d]+=qv[c][d];
  }
  AMPS_EXPECT_NEAR(charge_sum,q,4.0e-15);
  for (int d=0;d<3;++d)
    AMPS_EXPECT_NEAR(current_sum[d],q*velocity[d],2.0e-12);
}

AMPS_TEST_CASE("PIC-DEPOSITION-CONTINUITY-001",
               "segment current satisfies nodal continuity across multiple cells") {
  const double q=3.75,dt=0.2;
  std::vector<double> old_charge,new_charge,edge_current;
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::DepositChargeConservingCurrent1D(
      0.25,4.6,0.0,1.0,7,q,dt,old_charge,new_charge,edge_current)==
      PIC::DepositionKernels::Status::Success);
  AMPS_EXPECT_EQ(old_charge.size(),std::size_t(7));
  AMPS_EXPECT_EQ(edge_current.size(),std::size_t(6));

  for (std::size_t node=0;node<old_charge.size();++node) {
    const double incoming=node==0 ? 0.0 : edge_current[node-1];
    const double outgoing=node==edge_current.size() ? 0.0 : edge_current[node];
    const double residual=(new_charge[node]-old_charge[node])/dt+
                          outgoing-incoming;
    AMPS_EXPECT_NEAR(residual,0.0,3.0e-14);
  }
  AMPS_EXPECT_NEAR(std::accumulate(old_charge.begin(),old_charge.end(),0.0),
                   q,4.0e-15);
  AMPS_EXPECT_NEAR(std::accumulate(new_charge.begin(),new_charge.end(),0.0),
                   q,4.0e-15);
}

AMPS_TEST_CASE("PIC-DEPOSITION-INVALID-001",
               "deposition rejects nonfinite coordinates, invalid steps, and off-grid segments") {
  double weight[8];
  const double nan_xi[3]={0.5,std::numeric_limits<double>::quiet_NaN(),0.5};
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::TrilinearWeights(nan_xi,weight)==
                   PIC::DepositionKernels::Status::InvalidArgument);
  const double outside[3]={0.5,1.1,0.5};
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::TrilinearWeights(outside,weight)==
                   PIC::DepositionKernels::Status::OutsideDomain);
  std::vector<double> a,b,j;
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::DepositChargeConservingCurrent1D(
      0.5,1.5,0.0,1.0,4,1.0,0.0,a,b,j)==
      PIC::DepositionKernels::Status::InvalidArgument);
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::DepositChargeConservingCurrent1D(
      0.5,5.0,0.0,1.0,4,1.0,1.0,a,b,j)==
      PIC::DepositionKernels::Status::OutsideDomain);
}

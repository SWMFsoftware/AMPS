#include "test_harness.h"
#include "pic_field_line_kernels.h"
#include <cmath>
#include <numeric>
#include <vector>
namespace F=PIC::FieldLineKernels;

AMPS_TEST_CASE("PIC-FIELDLINE-STRAIGHT-INTERP-001","straight segment interpolation reproduces affine vector data") {
  const double a[3]={1.0,-2.0,4.0},b[3]={5.0,6.0,-8.0}; double out[3];
  AMPS_EXPECT_TRUE(F::LinearInterpolate(a,b,3,0.25,out)==F::Status::Success);
  AMPS_EXPECT_NEAR(out[0],2.0,0.0); AMPS_EXPECT_NEAR(out[1],0.0,0.0); AMPS_EXPECT_NEAR(out[2],1.0,0.0);
}
AMPS_TEST_CASE("PIC-FIELDLINE-CURVED-INTERP-001","arc interpolation follows a manufactured piecewise-curved polyline") {
  const std::vector<double> xyz={0,0,0, 3,0,0, 3,4,0, 3,4,12}; double out[3];
  AMPS_EXPECT_TRUE(F::InterpolatePolylineByArc(xyz,5.0,out)==F::Status::Success);
  AMPS_EXPECT_NEAR(out[0],3.0,0.0); AMPS_EXPECT_NEAR(out[1],2.0,0.0); AMPS_EXPECT_NEAR(out[2],0.0,0.0);
  AMPS_EXPECT_TRUE(F::InterpolatePolylineByArc(xyz,13.0,out)==F::Status::Success);
  AMPS_EXPECT_NEAR(out[2],6.0,0.0);
}
AMPS_TEST_CASE("PIC-FIELDLINE-EXCHANGE-001","vertex and edge payload permutations unpack each item exactly once") {
  const std::vector<std::size_t> order={2,0,3,1};
  const std::vector<double> packed={20,21, 0,1, 30,31, 10,11},expected={0,1,10,11,20,21,30,31};
  std::vector<double> unpacked;
  AMPS_EXPECT_TRUE(F::ExchangeIndexedValues(order,packed,4,2,unpacked)==F::Status::Success);
  AMPS_EXPECT_TRUE(unpacked==expected);
}
AMPS_TEST_CASE("PIC-FIELDLINE-PARTICLE-SPLIT-001","field-line splitting preserves represented number, momentum, and energy") {
  const F::WeightedParticle parent{0.35,7.0,-2.0,12.0}; std::vector<F::WeightedParticle> child;
  AMPS_EXPECT_TRUE(F::SplitParticle(parent,5,child)==F::Status::Success);
  double w=0,p=0,e=0,s=0;
  for (const auto& c:child) { w+=c.weight; p+=c.weight*c.v_parallel; e+=0.5*c.weight*(c.v_parallel*c.v_parallel+c.v_normal*c.v_normal); s+=c.weight*c.s; }
  AMPS_EXPECT_NEAR(w,parent.weight,2e-15); AMPS_EXPECT_NEAR(p,parent.weight*parent.v_parallel,2e-14);
  AMPS_EXPECT_NEAR(e,0.5*parent.weight*(53.0),2e-13); AMPS_EXPECT_NEAR(s,parent.weight*parent.s,2e-15);
}
AMPS_TEST_CASE("PIC-FIELDLINE-REPARTITION-001","overlap remap conserves field-line segment totals under repartition") {
  const std::vector<double> old_edges={0,1,3,6},old_total={2,8,18},new_edges={0,0.5,2,4,6};
  std::vector<double> remapped;
  AMPS_EXPECT_TRUE(F::ConservativeRepartition1D(old_edges,old_total,new_edges,remapped)==F::Status::Success);
  AMPS_EXPECT_NEAR(std::accumulate(remapped.begin(),remapped.end(),0.0),28.0,4e-15);
  AMPS_EXPECT_NEAR(remapped[0],1.0,0.0); AMPS_EXPECT_NEAR(remapped[1],5.0,0.0);
  AMPS_EXPECT_NEAR(remapped[2],10.0,0.0); AMPS_EXPECT_NEAR(remapped[3],12.0,0.0);
}
AMPS_TEST_CASE("PIC-FIELDLINE-INVALID-001","field-line kernels reject bad coordinates, topology, and duplicate ownership") {
  const double a[1]={0},b[1]={1}; double out[3];
  AMPS_EXPECT_TRUE(F::LinearInterpolate(a,b,1,1.1,out)==F::Status::OutsideDomain);
  std::vector<double> unpacked;
  AMPS_EXPECT_TRUE(F::ExchangeIndexedValues({0,0},{1,2},2,1,unpacked)==F::Status::DuplicateIndex);
  AMPS_EXPECT_TRUE(F::InterpolatePolylineByArc({0,0,0,0,0,0},0.0,out)==F::Status::InvalidArgument);
}

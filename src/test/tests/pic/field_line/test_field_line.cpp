#include "test_harness.h"
#include "pic_field_line_kernels.h"
#include "pic_field_line_population.h"
#include <cmath>
#include <numeric>
#include <vector>
namespace F=PIC::FieldLineKernels;
namespace S=PIC::ParticleSplitting::FledLine::ExactSplit;

namespace {
// Exercise the ordering/clone/hook kernel used by the native controller. This
// record carries opaque history whose meaning the core must not infer. Actual
// PIC buffer/MPI behavior remains a separate native integration gate.
struct Record {
  double s,parallel,normal,weight;
  std::uint64_t identity,events;
  double cachedNormal,hazard;
  std::vector<unsigned char> opaque;
};
struct Operations {
  std::vector<Record> records;
  std::vector<long int> clonedParents;
  Record copiedBeforeHook;
  void SetWeight(long int particle,double weight) {
    records.at(particle).weight=weight;
  }
  long int Allocate() {
    records.push_back(Record{});
    return static_cast<long int>(records.size()-1);
  }
  void Clone(long int child,long int parent) {
    records.at(child)=records.at(parent);
    copiedBeforeHook=records.at(child);
    clonedParents.push_back(parent);
  }
};
Operations* activeOperations=nullptr;
struct OperationsScope {
  Operations* saved;
  explicit OperationsScope(Operations& operations):saved(activeOperations) {
    activeOperations=&operations;
  }
  ~OperationsScope(){activeOperations=saved;}
};
std::uint64_t Identity(long int particle) {
  return activeOperations->records.at(particle).identity;
}
void ExactSplitHook(long int parent,long int child) {
  // This arbitrary metadata rule demonstrates explicit ownership only; it is
  // not a stochastic/physical prescription for an application.
  AMPS_EXPECT_TRUE(parent!=child);
  AMPS_EXPECT_TRUE(activeOperations->records.at(parent).weight==
                   activeOperations->records.at(child).weight);
  activeOperations->records.at(child).identity=
      activeOperations->records.at(parent).identity+1000;
}
void RejectExactSplit(long int) {
  throw std::runtime_error("controlled application lineage preflight rejection");
}
double Energy(const Record& record) {
  const double c=299792458.0,mass=1.67262192369e-27;
  const double v2=record.parallel*record.parallel+record.normal*record.normal;
  return record.weight*mass*c*c*(1.0/std::sqrt(1.0-v2/(c*c))-1.0);
}
void CheckExactLifecycle() {
  Operations operations;
  const Record original{0.35,1.2e8,0.8e8,12.0,71,83,0.125,0.47,
                        {0,17,255,31,92}};
  operations.records.push_back(original);
  OperationsScope scope(operations);
  std::vector<S::Candidate> candidates{{0,original.weight,0}};
  AMPS_EXPECT_TRUE(S::Prepare(candidates,Identity)==S::Status::Success);
  AMPS_EXPECT_TRUE(S::Apply(candidates.front(),operations,ExactSplitHook));
  const Record& parent=operations.records.at(0);
  const Record& child=operations.records.at(1);
  AMPS_EXPECT_TRUE(operations.clonedParents==std::vector<long int>{0});
  AMPS_EXPECT_TRUE(parent.identity==original.identity);
  AMPS_EXPECT_TRUE(child.identity==original.identity+1000);
  AMPS_EXPECT_TRUE(operations.copiedBeforeHook.identity==original.identity);
  AMPS_EXPECT_TRUE(parent.events==original.events&&child.events==original.events);
  AMPS_EXPECT_TRUE(parent.cachedNormal==original.cachedNormal&&
                   child.cachedNormal==original.cachedNormal);
  AMPS_EXPECT_TRUE(parent.hazard==original.hazard&&child.hazard==original.hazard);
  AMPS_EXPECT_TRUE(parent.opaque==original.opaque&&child.opaque==original.opaque);
  AMPS_EXPECT_TRUE(parent.s==original.s&&child.s==original.s&&
                   parent.parallel==original.parallel&&child.parallel==original.parallel&&
                   parent.normal==original.normal&&child.normal==original.normal);
  AMPS_EXPECT_TRUE(parent.weight+child.weight==original.weight);
  AMPS_EXPECT_NEAR(Energy(parent)+Energy(child),Energy(original),0.0);

  // Equal-weight selection follows physical identity under two different
  // allocation/list permutations instead of following the numeric handles.
  for(const std::vector<std::uint64_t>& order:
      {std::vector<std::uint64_t>{9,2,7,1,8},
       std::vector<std::uint64_t>{8,7,1,9,2}}) {
    operations.records.clear();operations.clonedParents.clear();candidates.clear();
    for(std::uint64_t id:order) {
      Record record=original;record.identity=id;record.weight=2.0;
      const long int pointer=static_cast<long int>(operations.records.size());
      operations.records.push_back(record);
      candidates.push_back({pointer,record.weight,0});
    }
    AMPS_EXPECT_TRUE(S::Prepare(candidates,Identity)==S::Status::Success);
    AMPS_EXPECT_TRUE(candidates[0].identity==1&&candidates[1].identity==2);
    AMPS_EXPECT_TRUE(S::Apply(candidates[0],operations,ExactSplitHook));
    AMPS_EXPECT_TRUE(operations.records.at(operations.clonedParents.front()).identity==1);
  }
}
void CheckInvalidLifecycle() {
  Operations operations;
  operations.records={{0,1,2,3,7,11,0.5,0.3,{91}},
                      {1,3,4,5,7,13,-0.5,0.7,{92}}};
  OperationsScope scope(operations);
  std::vector<S::Candidate> candidates{{0,3,0},{1,5,0}};
  AMPS_EXPECT_TRUE(S::Prepare(candidates,Identity)==S::Status::DuplicateIdentity);
  AMPS_EXPECT_TRUE(operations.clonedParents.empty()&&operations.records.size()==2);
  AMPS_EXPECT_TRUE(operations.records[0].weight==3&&operations.records[1].weight==5);
  operations.records[1].identity=0;
  AMPS_EXPECT_TRUE(S::Prepare(candidates,Identity)==S::Status::InvalidIdentity);
  candidates[1].weight=-1;
  AMPS_EXPECT_TRUE(S::Prepare(candidates,nullptr)==S::Status::InvalidWeight);
  candidates[1].weight=std::numeric_limits<double>::infinity();
  AMPS_EXPECT_TRUE(S::Prepare(candidates,nullptr)==S::Status::InvalidWeight);
  candidates[1].weight=std::numeric_limits<double>::min();
  const std::size_t size=operations.records.size();
  AMPS_EXPECT_TRUE(!S::Apply(candidates[1],operations,ExactSplitHook));
  AMPS_EXPECT_TRUE(operations.records.size()==size&&operations.clonedParents.empty());
  // A lineage/key/counter rejection reaches the actual pre-mutation hook,
  // leaving all physical weights, opaque history and allocation untouched.
  bool rejected=false;
  try {
    S::Apply(candidates[0],operations,ExactSplitHook,RejectExactSplit);
  } catch(const std::runtime_error&) { rejected=true; }
  AMPS_EXPECT_TRUE(rejected&&operations.records.size()==size&&
                   operations.clonedParents.empty());
  AMPS_EXPECT_TRUE(operations.records[0].weight==3&&
                   operations.records[0].events==11&&
                   operations.records[0].hazard==0.3&&
                   operations.records[0].opaque==std::vector<unsigned char>{91});
}
} // namespace

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
  CheckExactLifecycle();
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
  CheckInvalidLifecycle();
}

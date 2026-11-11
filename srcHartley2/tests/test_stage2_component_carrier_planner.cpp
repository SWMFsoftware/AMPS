#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <vector>

#include "srcHartley2/ComponentCarrierPlanner.h"

namespace Planner=Comet::RadialSimulation::ComponentCarrierPlanner;

namespace {

bool Equal(long double left,long double right) {
  return std::fabs(left-right)<=1.0e-13L*
      std::max(1.0L,std::max(std::fabs(left),std::fabs(right)));
}

void Require(bool condition,const char *message) {
  if (!condition) throw std::runtime_error(message);
}

template<class Exception,class Function>
void RequireThrows(Function function,const char *message) {
  bool threw=false;
  try { function(); }
  catch (const Exception&) { threw=true; }
  Require(threw,message);
}

std::vector<Planner::cEdge> Chain(std::size_t size) {
  std::vector<Planner::cEdge> edges;
  for (std::size_t node=1;node<size;node++)
    edges.push_back(Planner::cEdge(node-1,node));
  return edges;
}

void TestOneBlockAndExactBoundaries() {
  const Planner::cPolicy policy(10,40,false);
  Planner::cPlan plan=Planner::BuildPlan(
      {Planner::cNodeState(40,8,1.0)}, {},policy);
  Require(plan.Components.size()==1,"one block must form one component");
  Require(plan.WeightPower[0]==-1,"sparse one-block component must split");
  Require(Equal(plan.Components[0].ExpectedParticlesPerCell,10.0L),
      "one-block split must reach the lower boundary exactly");
  Require(plan.RequiredAdditionalParticles==40,
      "one-block split clone count is not exact");

  plan=Planner::BuildPlan({Planner::cNodeState(80,8,1.0)}, {},policy);
  Require(plan.WeightPower[0]==0,
      "the exact lower boundary must not change weight");
  plan=Planner::BuildPlan({Planner::cNodeState(320,8,1.0)}, {},policy);
  Require(plan.WeightPower[0]==0,
      "the exact upper boundary must not change weight");
  plan=Planner::BuildPlan({Planner::cNodeState(321,8,1.0)}, {},policy);
  Require(plan.WeightPower[0]==1,
      "a population above the upper boundary must roulette once");
  Require(Equal(plan.Components[0].ExpectedParticlesPerCell,20.0625L),
      "roulette prediction must use the exact dyadic factor");
}

void TestLongRadialChainAndEmptyNodes() {
  std::vector<Planner::cNodeState> nodes(17,Planner::cNodeState(0,8,2.0));
  nodes[8].ParticleCount=34;
  Planner::cPlan plan=Planner::BuildPlan(nodes,Chain(nodes.size()),
      Planner::cPolicy(1,4,false));
  Require(plan.Components.size()==1,"radial chain must be connected");
  Require(plan.Components[0].Nodes.size()==17,
      "radial traversal did not reach the full chain");
  for (std::size_t node=0;node<nodes.size();node++)
    Require(plan.WeightPower[node]==-2,
        "one component must receive one common split power");
  Require(Equal(plan.Components[0].TargetCarrier,0.5L),
      "radial component target carrier is wrong");
  Require(plan.RequiredAdditionalParticles==102,
      "only the populated block may contribute split allocations");

  // A wholly empty component cannot be brought into a positive population
  // band without seeding it.  Equalize it to the largest existing dyadic
  // carrier using weight-only changes, which allocate no particles.
  nodes={Planner::cNodeState(0,8,1.0),Planner::cNodeState(0,8,2.0),
      Planner::cNodeState(0,8,4.0)};
  plan=Planner::BuildPlan(nodes,Chain(nodes.size()),
      Planner::cPolicy(1,4,false));
  Require(plan.WeightPower[0]==2 && plan.WeightPower[1]==1 &&
      plan.WeightPower[2]==0,"empty component did not equalize safely");
  Require(plan.RequiredAdditionalParticles==0,
      "empty component must never allocate particles");
  Require(Equal(plan.Components[0].ExpectedParticlesPerCell,0.0L),
      "empty component must remain empty");
}

void TestPopulatedEndsAndDisconnectedComponents() {
  std::vector<Planner::cNodeState> nodes(9,Planner::cNodeState(0,8,1.0));
  nodes.front().ParticleCount=18;
  nodes.back().ParticleCount=18;
  Planner::cPlan plan=Planner::BuildPlan(nodes,Chain(nodes.size()),
      Planner::cPolicy(1,4,false));
  Require(plan.WeightPower.front()==-1 && plan.WeightPower.back()==-1,
      "both populated ends must use the same component decision");

  nodes={Planner::cNodeState(8,8,1.0),Planner::cNodeState(0,8,1.0),
      Planner::cNodeState(320,8,1.0),Planner::cNodeState(320,8,1.0)};
  const std::vector<Planner::cEdge> edges={Planner::cEdge(0,1),
      Planner::cEdge(2,3)};
  plan=Planner::BuildPlan(nodes,edges,Planner::cPolicy(1,4,false));
  Require(plan.Components.size()==2,
      "disconnected graph must retain independent components");
  Require(plan.WeightPower[0]==-1 && plan.WeightPower[1]==-1,
      "sparse component selected the wrong power");
  Require(plan.WeightPower[2]==4 && plan.WeightPower[3]==4,
      "dense component selected the wrong power");
  Require(plan.ComponentIndex[0]!=plan.ComponentIndex[2],
      "disconnected nodes were assigned the same component");
}

std::vector<Planner::cEdge> GridEdges(int nx,int ny,int nz) {
  std::vector<Planner::cEdge> edges;
  const auto index=[=](int i,int j,int k) {
    return static_cast<std::size_t>(i+nx*(j+ny*k));
  };
  for (int k=0;k<nz;k++) for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) {
    if (i+1<nx) edges.push_back(Planner::cEdge(index(i,j,k),index(i+1,j,k)));
    if (j+1<ny) edges.push_back(Planner::cEdge(index(i,j,k),index(i,j+1,k)));
    if (k+1<nz) edges.push_back(Planner::cEdge(index(i,j,k),index(i,j,k+1)));
  }
  return edges;
}

void TestSyntheticHigherDimensionalFaceGraphs() {
  // A 3x3 DIM2 face graph.  Diagonal/corner adjacency is deliberately absent.
  std::vector<Planner::cNodeState> nodes2(9,Planner::cNodeState(0,8,1.0));
  nodes2[4].ParticleCount=72;
  Planner::cPlan plan=Planner::BuildPlan(nodes2,GridEdges(3,3,1),
      Planner::cPolicy(1,4,false));
  Require(plan.Components.size()==1 && plan.Components[0].Nodes.size()==9,
      "DIM2 face graph is not one component");
  for (int power:plan.WeightPower)
    Require(power==0,"DIM2 exact lower boundary must not adapt");

  // A 2x2x2 DIM3 face graph exercises all three face directions.
  std::vector<Planner::cNodeState> nodes3(8,Planner::cNodeState(0,8,1.0));
  nodes3[7].ParticleCount=32;
  plan=Planner::BuildPlan(nodes3,GridEdges(2,2,2),
      Planner::cPolicy(1,4,false));
  Require(plan.Components.size()==1 && plan.Components[0].Nodes.size()==8,
      "DIM3 face graph is not one component");
  for (int power:plan.WeightPower)
    Require(power==-1,"DIM3 component-wide split decision is wrong");
}

void TestDyadicInputsLocksAndCapacityPreflight() {
  std::vector<Planner::cNodeState> nodes={
      Planner::cNodeState(0,8,1.0),Planner::cNodeState(0,8,2.0),
      Planner::cNodeState(8,8,4.0)};
  Planner::cPlan plan=Planner::BuildPlan(nodes,Chain(nodes.size()),
      Planner::cPolicy(1,4,false));
  Require(plan.WeightPower[0]==0 && plan.WeightPower[1]==-1 &&
      plan.WeightPower[2]==-2,
      "unequal dyadic inputs did not converge on the common lower carrier");
  Require(plan.RequiredAdditionalParticles==24,
      "dyadic input split capacity is not exact");

  nodes={Planner::cNodeState(1000,8,1.0),
      Planner::cNodeState(1000,8,1.0)};
  plan=Planner::BuildPlan(nodes,Chain(nodes.size()),
      Planner::cPolicy(1,4,true));
  Require(plan.WeightPower[0]==0 && plan.WeightPower[1]==0,
      "locked direct source must remain bitwise unchanged");

  RequireThrows<std::invalid_argument>([]() {
    Planner::BuildPlan({Planner::cNodeState(1,8,1.0),
        Planner::cNodeState(1,8,2.0)},Chain(2),Planner::cPolicy(1,4,true));
  },"unequal locked carriers must be rejected transactionally");
  RequireThrows<std::invalid_argument>([]() {
    Planner::BuildPlan({Planner::cNodeState(1,8,1.0),
        Planner::cNodeState(1,8,3.0)},Chain(2),Planner::cPolicy(1,4,false));
  },"non-dyadic carriers must be rejected");

  // BuildPlan is pure: a capacity failure is reported before the runtime can
  // mutate any carrier or particle list.
  RequireThrows<std::length_error>([]() {
    Planner::BuildPlan({Planner::cNodeState(40,8,1.0)}, {},
        Planner::cPolicy(10,40,false,39));
  },"insufficient split capacity must reject the complete plan");
}

} // namespace

int main() {
  try {
    TestOneBlockAndExactBoundaries();
    TestLongRadialChainAndEmptyNodes();
    TestPopulatedEndsAndDisconnectedComponents();
    TestSyntheticHigherDimensionalFaceGraphs();
    TestDyadicInputsLocksAndCapacityPreflight();
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,"Stage-2 component carrier planner test failed: %s\n",
        exception.what());
    return 1;
  }
  std::printf("PASS Stage-2 component carrier planner\n");
  return 0;
}

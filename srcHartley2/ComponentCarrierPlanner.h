#ifndef _COMET_COMPONENT_CARRIER_PLANNER_H_
#define _COMET_COMPONENT_CARRIER_PLANNER_H_

// Geometry- and MPI-independent component carrier planning.
//
// AMPS owns particles and block-associated weights on MPI ranks, but the
// numerical decision made at a safe adaptation boundary is a pure graph
// problem.  For one species, every calculation block is a graph node and
// every traversable block face is an undirected edge.  A face transfer
// multiplicity of one requires a single W/dt carrier lattice value over each
// connected component.  Keeping this calculation free of AMPS types makes
// its DIM1, DIM2, and DIM3 behavior testable without manufacturing a mesh or
// invoking MPI.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <climits>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace Comet {
namespace RadialSimulation {
namespace ComponentCarrierPlanner {

struct cNodeState {
  // Actual model-particle population for this species on the authoritative
  // owner.  Zero is meaningful: changing an empty node's carrier must not
  // create a particle.
  long long ParticleCount;
  // Number of physical mesh cells represented by this block.  The requested
  // population policy is stated per cell, so components with unequal block
  // sizes must use the sum of these values rather than a block average.
  long long CellCount;
  // Current W/dt.  All carriers for a species must belong to one dyadic
  // lattice because every supported adaptation changes W by a power of two.
  double Carrier;

  cNodeState(long long particleCount=0,long long cellCount=0,
      double carrier=0.0) :
      ParticleCount(particleCount),CellCount(cellCount),Carrier(carrier) {}
};

struct cEdge {
  std::size_t First,Second;
  cEdge(std::size_t first=0,std::size_t second=0) :
      First(first),Second(second) {}
};

struct cPolicy {
  long long MinimumParticlesPerCell;
  long long MaximumParticlesPerCell;
  bool Locked;
  // The runtime supplies the globally safe number of additional particle
  // slots.  Planning is pure and throws before the caller mutates any block
  // when the exact split requirement exceeds this value.
  long long MaximumAdditionalParticles;

  cPolicy(long long minimum=0,long long maximum=0,bool locked=false,
      long long maximumAdditionalParticles=LLONG_MAX) :
      MinimumParticlesPerCell(minimum),
      MaximumParticlesPerCell(maximum),Locked(locked),
      MaximumAdditionalParticles(maximumAdditionalParticles) {}
};

struct cComponentPlan {
  std::vector<std::size_t> Nodes;
  // Representative target on the first node's carrier lattice.  Each node
  // is changed by WeightPower[node]; small roundoff inherited from W/dt may
  // leave the independently calculated node targets different by a few ulps,
  // which the runtime face check already treats with its explicit tolerance.
  double TargetCarrier;
  long double ExpectedParticlesPerCell;
  long long RequiredAdditionalParticles;

  cComponentPlan() : TargetCarrier(0.0),ExpectedParticlesPerCell(0.0L),
      RequiredAdditionalParticles(0) {}
};

struct cPlan {
  // Target carrier for node i is exactly ldexp(oldCarrier,
  // WeightPower[i]).  A negative value is a split; a positive value is
  // roulette/coarsening; zero leaves the node unchanged.
  std::vector<int> WeightPower;
  std::vector<std::size_t> ComponentIndex;
  std::vector<cComponentPlan> Components;
  long long RequiredAdditionalParticles;

  cPlan() : RequiredAdditionalParticles(0) {}
};

namespace Detail {

inline bool NearlyEqual(double left,double right) {
  const double scale=std::max(std::fabs(left),std::fabs(right));
  return std::fabs(left-right)<=
      512.0*std::numeric_limits<double>::epsilon()*scale;
}

inline int DyadicExponent(double value,double anchor) {
  if (!(value>0.0) || !std::isfinite(value) ||
      !(anchor>0.0) || !std::isfinite(anchor))
    throw std::invalid_argument("component carrier must be finite and positive");

  // Subtract logarithms instead of forming value/anchor directly; the latter
  // can overflow even when the finite dyadic exponent is representable.
  const long double raw=std::log2(static_cast<long double>(value))-
      std::log2(static_cast<long double>(anchor));
  if (!std::isfinite(raw) || raw<-1024.0L || raw>1024.0L)
    throw std::invalid_argument("component carriers exceed the supported dyadic lattice");
  const long long rounded=static_cast<long long>(std::llround(raw));
  if (rounded<-1024 || rounded>1024)
    throw std::invalid_argument("component carrier exponent is out of range");
  const long double expected=std::ldexp(
      static_cast<long double>(anchor),static_cast<int>(rounded));
  const long double scale=std::max(std::fabs(expected),
      std::fabs(static_cast<long double>(value)));
  if (std::fabs(expected-static_cast<long double>(value))>
      512.0L*std::numeric_limits<double>::epsilon()*scale)
    throw std::invalid_argument("component carriers do not share a dyadic lattice");
  return static_cast<int>(rounded);
}

inline long long RequiredClones(long long count,int power) {
  if (count<0) throw std::invalid_argument("negative component particle count");
  if (count==0 || power>=0) return 0;

  long long multiplier=1;
  for (int step=0;step < -power;step++) {
    if (multiplier>LLONG_MAX/2)
      throw std::overflow_error("component split multiplier overflows long long");
    multiplier*=2;
  }
  const long long additionalPerParticle=multiplier-1;
  if (additionalPerParticle>0 && count>LLONG_MAX/additionalPerParticle)
    throw std::overflow_error("component split population overflows long long");
  return count*additionalPerParticle;
}

inline void AddChecked(long long value,long long &total,
    const char *message) {
  if (value<0 || total>LLONG_MAX-value) throw std::overflow_error(message);
  total+=value;
}

} // namespace Detail

inline cPlan BuildPlan(const std::vector<cNodeState> &nodes,
    const std::vector<cEdge> &edges,const cPolicy &policy) {
  if (nodes.empty())
    throw std::invalid_argument("component planner received no nodes");
  if (policy.MinimumParticlesPerCell<=0 ||
      policy.MaximumParticlesPerCell<policy.MinimumParticlesPerCell ||
      policy.MinimumParticlesPerCell>policy.MaximumParticlesPerCell/2)
    throw std::invalid_argument("component population range is invalid");
  if (policy.MaximumAdditionalParticles<0)
    throw std::invalid_argument("component clone capacity is negative");

  std::vector<std::vector<std::size_t> > neighbors(nodes.size());
  for (std::size_t edgeIndex=0;edgeIndex<edges.size();edgeIndex++) {
    const cEdge &edge=edges[edgeIndex];
    if (edge.First>=nodes.size() || edge.Second>=nodes.size())
      throw std::invalid_argument("component edge endpoint is out of range");
    if (edge.First==edge.Second) continue;
    neighbors[edge.First].push_back(edge.Second);
    neighbors[edge.Second].push_back(edge.First);
  }
  for (std::size_t node=0;node<nodes.size();node++) {
    if (nodes[node].ParticleCount<0 || nodes[node].CellCount<=0 ||
        !(nodes[node].Carrier>0.0) || !std::isfinite(nodes[node].Carrier))
      throw std::invalid_argument("component node state is invalid");
    std::sort(neighbors[node].begin(),neighbors[node].end());
    neighbors[node].erase(std::unique(neighbors[node].begin(),
        neighbors[node].end()),neighbors[node].end());
  }

  cPlan result;
  result.WeightPower.assign(nodes.size(),0);
  result.ComponentIndex.assign(nodes.size(),nodes.size());
  std::vector<unsigned char> visited(nodes.size(),0);

  for (std::size_t start=0;start<nodes.size();start++) {
    if (visited[start]) continue;
    cComponentPlan component;
    std::queue<std::size_t> work;
    work.push(start);
    visited[start]=1;
    while (!work.empty()) {
      const std::size_t node=work.front();
      work.pop();
      component.Nodes.push_back(node);
      for (std::size_t neighborIndex=0;
           neighborIndex<neighbors[node].size();neighborIndex++) {
        const std::size_t neighbor=neighbors[node][neighborIndex];
        if (visited[neighbor]) continue;
        visited[neighbor]=1;
        work.push(neighbor);
      }
    }

    const std::size_t componentIndex=result.Components.size();
    const double anchor=nodes[component.Nodes.front()].Carrier;
    std::vector<int> exponent(component.Nodes.size(),0);
    int maximumExponent=-1024;
    long long totalCells=0,totalParticles=0;
    for (std::size_t entry=0;entry<component.Nodes.size();entry++) {
      const std::size_t node=component.Nodes[entry];
      exponent[entry]=Detail::DyadicExponent(nodes[node].Carrier,anchor);
      maximumExponent=std::max(maximumExponent,exponent[entry]);
      Detail::AddChecked(nodes[node].CellCount,totalCells,
          "component cell count overflows long long");
      Detail::AddChecked(nodes[node].ParticleCount,totalParticles,
          "component particle count overflows long long");
      result.ComponentIndex[node]=componentIndex;
    }

    int targetExponent=maximumExponent;
    long double expectedMean=0.0L;
    if (totalParticles>0) {
      // Predict the component population after adopting the highest existing
      // carrier.  A particle in node i represents Carrier[i]/TargetCarrier
      // target particles, so this sum conserves represented physical number
      // even when the input carriers have not yet been equalized.
      long double expectedParticles=0.0L;
      for (std::size_t entry=0;entry<component.Nodes.size();entry++) {
        const std::size_t node=component.Nodes[entry];
        expectedParticles+=static_cast<long double>(nodes[node].ParticleCount)*
            std::ldexp(1.0L,exponent[entry]-targetExponent);
      }
      expectedMean=expectedParticles/static_cast<long double>(totalCells);

      if (!policy.Locked) {
        while (expectedMean<policy.MinimumParticlesPerCell) {
          if (targetExponent<=-1024)
            throw std::overflow_error("component split exponent is out of range");
          targetExponent--;
          expectedMean*=2.0L;
        }
        while (expectedMean>policy.MaximumParticlesPerCell) {
          if (targetExponent>=1024)
            throw std::overflow_error("component roulette exponent is out of range");
          targetExponent++;
          expectedMean*=0.5L;
        }
      }
    }

    if (policy.Locked) {
      // A locked source species must remain bitwise unchanged.  Face
      // multiplicity one is therefore possible only when the component was
      // already equal (within the runtime's W/dt tolerance).  Returning zero
      // powers guarantees that planning never silently edits a source weight.
      for (std::size_t entry=1;entry<component.Nodes.size();entry++)
        if (!Detail::NearlyEqual(nodes[component.Nodes[entry]].Carrier,anchor))
          throw std::invalid_argument(
              "locked component has unequal carriers and cannot be equalized");
      targetExponent=0;
      expectedMean=static_cast<long double>(totalParticles)/totalCells;
    }

    const long double targetLong=std::ldexp(
        static_cast<long double>(anchor),targetExponent);
    if (!(targetLong>0.0L) || !std::isfinite(targetLong) ||
        targetLong>std::numeric_limits<double>::max())
      throw std::overflow_error("component target carrier is not representable");
    component.TargetCarrier=static_cast<double>(targetLong);
    component.ExpectedParticlesPerCell=expectedMean;

    for (std::size_t entry=0;entry<component.Nodes.size();entry++) {
      const std::size_t node=component.Nodes[entry];
      const int power=policy.Locked ? 0 : targetExponent-exponent[entry];
      if (power<-1024 || power>1024)
        throw std::overflow_error("component node weight power is out of range");
      result.WeightPower[node]=power;
      const long long clones=Detail::RequiredClones(
          nodes[node].ParticleCount,power);
      Detail::AddChecked(clones,component.RequiredAdditionalParticles,
          "component clone requirement overflows long long");
      Detail::AddChecked(clones,result.RequiredAdditionalParticles,
          "global clone requirement overflows long long");
    }
    if (result.RequiredAdditionalParticles>
        policy.MaximumAdditionalParticles)
      throw std::length_error(
          "component plan exceeds the available particle capacity");
    result.Components.push_back(component);
  }

  return result;
}

} // namespace ComponentCarrierPlanner
} // namespace RadialSimulation
} // namespace Comet

#endif

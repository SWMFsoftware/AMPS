#include "test_harness.h"
#include "interpolation/complete_stencil.h"
#include <cstring>

namespace {
struct Node {double x[3]; const double* GetX() const {return x;}};
struct Stencil {int Length=0; double Weight[8]{}; Node* cell[8]{}; int LocalCellID[8]{};};
struct Fixture {
  Node nodes[8]; Node* pointers[8]; int ids[8]; double weights[8];
  double lower[3]={-1,-1,-1},upper[3]={1,1,1};
  Fixture() {
    // Deliberately nonmonotone IDs and unequal/zero weights detect reordering;
    // node coordinates are independent of their serialized local index.
    for (int i=0;i<8;i++) {
      for (int d=0;d<3;d++) nodes[i].x[d]=(i&(1<<d)) ? .5:-.5;
      pointers[i]=nodes+i; ids[i]=37-3*i; weights[i]=double(i)/28.0;
    }
  }
  bool Assign(Stencil& s,bool periodic=false) {
    return amps_complete_stencil::Assign(pointers,ids,weights,periodic,lower,upper,3,s);
  }
};
}

AMPS_TEST_CASE("PERF-STENCIL-EXACT-001","complete stencil preserves nonmonotone IDs, zero weights and insertion order") {
  Fixture f; Stencil s;
  AMPS_EXPECT_TRUE(f.Assign(s)); AMPS_EXPECT_EQ(s.Length,8);
  for (int i=0;i<8;i++) {
    AMPS_EXPECT_EQ(s.LocalCellID[i],37-3*i);
    AMPS_EXPECT_TRUE(s.cell[i]==f.nodes+i);
    AMPS_EXPECT_NEAR(s.Weight[i],double(i)/28.0,0.0);
  }
}

AMPS_TEST_CASE("PERF-STENCIL-FALLBACK-001","missing, outside or negative support rejects without partial writes") {
  Fixture f; Stencil s; s.Length=3;
  for (int i=0;i<8;i++) s.Weight[i]=91.0;
  f.pointers[7]=nullptr;
  AMPS_EXPECT_TRUE(!f.Assign(s)); AMPS_EXPECT_EQ(s.Length,3);
  for (int i=0;i<8;i++) AMPS_EXPECT_NEAR(s.Weight[i],91.0,0.0);
  f.pointers[7]=f.nodes+7; f.nodes[7].x[0]=1.01;
  AMPS_EXPECT_TRUE(!f.Assign(s)); AMPS_EXPECT_EQ(s.Length,3);
  AMPS_EXPECT_TRUE(f.Assign(s,true)); // Periodic support intentionally ignores global-box clipping.
  f.weights[2]=-.01; s.Length=3;
  AMPS_EXPECT_TRUE(!f.Assign(s,true)); AMPS_EXPECT_EQ(s.Length,3);
}

AMPS_TEST_CASE("PERF-STENCIL-REFRESH-001","successive queries read new nodes and weights without stale topology") {
  Fixture first,second; Stencil s;
  AMPS_EXPECT_TRUE(first.Assign(s));
  for (int i=0;i<8;i++) {second.ids[i]=i;second.weights[i]=double(7-i)/28.0;}
  AMPS_EXPECT_TRUE(second.Assign(s));
  for (int i=0;i<8;i++) {
    AMPS_EXPECT_TRUE(s.cell[i]==second.nodes+i);
    AMPS_EXPECT_EQ(s.LocalCellID[i],i);
    AMPS_EXPECT_NEAR(s.Weight[i],double(7-i)/28.0,0.0);
  }
}

#ifndef AMPS_CORE_MESH_TEST_FIXTURE_H
#define AMPS_CORE_MESH_TEST_FIXTURE_H

#include "pic.h"
#include "src/test/native/test_harness.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace AMPS {
namespace CoreMeshTests {

using Node = cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;

#if _MESH_DIMENSION_ == 1
#define AMPS_MESH_PROFILE "1D"
#elif _MESH_DIMENSION_ == 2
#define AMPS_MESH_PROFILE "2D"
#elif _MESH_DIMENSION_ == 3
#define AMPS_MESH_PROFILE "3D"
#else
#error "Core mesh tests require DIM=1, DIM=2, or DIM=3"
#endif

#define AMPS_MESH_CASE_ID(CONTRACT) \
  "CORE-MESH-" AMPS_MESH_PROFILE "-" CONTRACT "-001"

// The fixture is intentionally small and dyadic: split the unit root and then
// split its all-lower child.  Every coordinate is exactly representable in
// binary, so topology failures cannot be hidden by a floating-point tolerance.
// The leaves cover unit D-volume and contain both same-level and 2:1 AMR
// interfaces.  Tests derive their expected topology independently from D.
void EnsureTopologyMesh();
void EnsureAllocatedMesh();
void EnsureInterpolationRoutines();
Node* Root();
std::vector<Node*> Leaves();
std::vector<Node*> AllNodes();
double DyadicLeafVolume(const Node* node);
std::uint64_t CanonicalTopologyChecksum();
double UnitLoad(Node* node);
double LevelLoad(Node* node);

// Static-archive anchors.  main.cpp references each symbol so every object
// containing registrations is retained by the application archive linker.
void LinkMeshTopologyTests();
void LinkMeshAdaptationTests();
void LinkMeshCoarsenDataTests();
void LinkPICInterpolationTests();
void LinkMeshCutCellTests();
void RegisterMeshCoarsenData();

}  // namespace CoreMeshTests
}  // namespace AMPS

#endif

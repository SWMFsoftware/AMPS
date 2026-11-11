#ifndef SEP_CORONA_SWCME_CROSS_RAY_MESH_H
#define SEP_CORONA_SWCME_CROSS_RAY_MESH_H

#include "sep_corona_swcme/spherical_sector.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

inline constexpr std::size_t NoMeshIndex=std::numeric_limits<std::size_t>::max();
enum class MeshAxis { Radial,Mu,Phi };

struct StaticMeshInput {
  std::vector<double> radiusEdgesM,muEdges,phiEdgesRad;
  SectorBasis basis;
  bool periodicPhi=false;
  // Explicit caller resource budget, not a physical resolution or mesh floor.
  std::size_t maximumCells=0;
};

struct CurvilinearPoint {
  double radiusM=0,mu=0,phiRad=0;
};

struct OrientedMeshEdge {
  std::size_t index=NoMeshIndex;
  int orientation=0;
};

struct MeshEdge {
  // One canonical direction from the smaller node ID to the larger node ID.
  // A curved mu/phi edge is identified by topology, not a straight chord.
  std::size_t firstNode=NoMeshIndex,lastNode=NoMeshIndex;
  MeshAxis axis=MeshAxis::Radial;
  CurvilinearPoint start,end;
  double arcLengthM=0;
};

struct MeshFace {
  // Normal points in increasing r, mu or phi, from left to right. A missing
  // neighbor is an explicit domain boundary; no ambient/wall law is inferred.
  std::size_t leftCell=NoMeshIndex,rightCell=NoMeshIndex;
  MeshAxis axis=MeshAxis::Radial;
  SectorFaceMeasure measure;
  std::vector<std::size_t> positiveNodeLoop;
  std::vector<OrientedMeshEdge> positiveEdgeLoop;
};

struct MeshFaceUse {
  std::size_t index=NoMeshIndex;
  int outwardOrientation=0;
  bool collapsedGeometricFace=false;
};

struct MeshCell {
  SectorBounds bounds;
  SectorMeasure measure;
  std::array<MeshFaceUse,6> faces;
  CoronalCME::Vec3 sharedAreaClosureResidualM2;
};

struct EdgeSample {
  CoronalCME::Vec3 positionHciM;
  // dX/ds, metres, for dimensionless s in [0,1]. Integrating E.dX/ds ds
  // yields volts, the oriented EMF required by Faraday's law for CT.
  CoronalCME::Vec3 derivativeM;
};

// One immutable STATIC shared-face/edge authority. It verifies exact integer
// boundary-of-boundary cancellation for every cell, pole/seam canonicalization,
// physical measures and all face incidences before publishing a candidate.
// Motion/remapping, legacy-ray conversion, piston cuts and MPI publication
// remain unimplemented and cannot be inferred from a static mesh PASS.
class StaticSphericalMesh final {
 public:
  static Core::Result<std::shared_ptr<const StaticSphericalMesh>> Create(
      StaticMeshInput input,std::string eventFingerprint,
      double epochS,std::uint64_t generation);

  const std::vector<MeshCell>& Cells() const noexcept {return cells_;}
  const std::vector<MeshFace>& Faces() const noexcept {return faces_;}
  const std::vector<MeshEdge>& Edges() const noexcept {return edges_;}
  const std::vector<CoronalCME::Vec3>& NodesHciM() const noexcept {return nodes_;}
  const std::string& Fingerprint() const noexcept {return fingerprint_;}
  const std::string& EventFingerprint() const noexcept {return eventFingerprint_;}
  double EpochS() const noexcept {return epochS_;}
  std::uint64_t Generation() const noexcept {return generation_;}
  Core::Result<EdgeSample> SampleEdge(std::size_t edge,double parameter) const;

 private:
  // Only a validated factory may create an authority; an empty public default
  // object would bypass mesh/epoch/positivity checks before a later query.
  StaticSphericalMesh()=default;
  StaticMeshInput input_;
  std::string eventFingerprint_,fingerprint_;
  double epochS_=0;
  std::uint64_t generation_=0;
  std::vector<MeshCell> cells_;
  std::vector<MeshFace> faces_;
  std::vector<MeshEdge> edges_;
  std::vector<CoronalCME::Vec3> nodes_;
};

} } } // namespace SEP::CoronaSwcme::CrossRay

#endif

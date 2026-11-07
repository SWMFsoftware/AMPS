#include "mesh_test_fixture.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits.h>
#include <string>
#include <vector>

#include <unistd.h>

namespace AMPS {
namespace CoreMeshTests {

#if _MESH_DIMENSION_ == 3

namespace {

const double kPi = 3.141592653589793238462643383279502884;

void SetTriangle(CutCell::cTriangleFace& face,
                 double x0, double y0, double z0,
                 double x1, double y1, double z1,
                 double x2, double y2, double z2) {
  double a[3] = {x0, y0, z0};
  double b[3] = {x1, y1, z1};
  double c[3] = {x2, y2, z2};
  face.SetFaceNodes(a, b, c);
}

double TriangleArea(const double* x0, const double* x1, const double* x2) {
  const double e0[3] = {x1[0]-x0[0], x1[1]-x0[1], x1[2]-x0[2]};
  const double e1[3] = {x2[0]-x0[0], x2[1]-x0[1], x2[2]-x0[2]};
  const double cross[3] = {
      e0[1]*e1[2]-e0[2]*e1[1],
      e0[2]*e1[0]-e0[0]*e1[2],
      e0[0]*e1[1]-e0[1]*e1[0]};
  return 0.5*std::sqrt(cross[0]*cross[0]+cross[1]*cross[1]+
                       cross[2]*cross[2]);
}

// A closed rectangular solid supplies an independent planar-volume oracle.
// The 12 triangles are wound outward.  In the plane-volume case its near face
// is x=0.30 m and the other five faces lie outside the sampled [-1,1]^3 m
// cell, so the exact remaining volume is the elementary prism integral
// integral_{-1}^{0.30} dx integral_{-1}^{1}dy integral_{-1}^{1}dz = 5.2 m^3.
struct OrientedBox {
  CutCell::cNASTRANnode nodes[8];
  CutCell::cTriangleFace faces[12];

  OrientedBox(double xmin, double xmax, double ymin, double ymax,
              double zmin, double zmax) {
    const double coordinates[8][3] = {
        {xmin,ymin,zmin}, {xmax,ymin,zmin},
        {xmax,ymax,zmin}, {xmin,ymax,zmin},
        {xmin,ymin,zmax}, {xmax,ymin,zmax},
        {xmax,ymax,zmax}, {xmin,ymax,zmax}};
    for (int i=0;i<8;i++) {
      std::memcpy(nodes[i].x,coordinates[i],3*sizeof(double));
      nodes[i].id=i;
    }

    // z-, z+, y-, y+, x-, x+; right-hand-rule normals point outward.
    const int connectivity[12][3] = {
        {0,2,1}, {0,3,2}, {4,5,6}, {4,6,7},
        {0,1,5}, {0,5,4}, {3,7,6}, {3,6,2},
        {0,4,7}, {0,7,3}, {1,2,6}, {1,6,5}};
    for (int i=0;i<12;i++) {
      faces[i].SetFaceNodes(nodes[connectivity[i][0]].x,
                            nodes[connectivity[i][1]].x,
                            nodes[connectivity[i][2]].x);
      for (int j=0;j<3;j++) faces[i].node[j]=nodes+connectivity[i][j];
    }
  }
};

// Surface grid resolution is static production state.  Preserve all derived
// angular increments so reverse/shuffled case execution has the same state as
// normal order.
class SphereGridScope {
 public:
  SphereGridScope()
      : nz_(cInternalSphericalData::nZenithSurfaceElements),
        na_(cInternalSphericalData::nAzimuthalSurfaceElements),
        da_(cInternalSphericalData::dAzimuthalAngle),
        dc_(cInternalSphericalData::dCosZenithAngle),
        dz_(cInternalSphericalData::dZenithAngle) {}

  ~SphereGridScope() {
    cInternalSphericalData::nZenithSurfaceElements=nz_;
    cInternalSphericalData::nAzimuthalSurfaceElements=na_;
    cInternalSphericalData::dAzimuthalAngle=da_;
    cInternalSphericalData::dCosZenithAngle=dc_;
    cInternalSphericalData::dZenithAngle=dz_;
  }

 private:
  long int nz_,na_;
  double da_,dc_,dz_;
};

double PlanarSphereArea(cInternalSphericalData& sphere,
                        int nZenith, int nAzimuth) {
  cInternalSphericalData::SetGeneralSurfaceMeshParameters(nZenith,nAzimuth);
  double area=0.0;
  for (int i=0;i<nZenith;i++) for (int j=0;j<nAzimuth;j++) {
    double x00[3],x01[3],x11[3],x10[3];
    sphere.GetSurfaceCoordinate(x00,i,j);
    sphere.GetSurfaceCoordinate(x01,i,j+1);
    sphere.GetSurfaceCoordinate(x11,i+1,j+1);
    sphere.GetSurfaceCoordinate(x10,i+1,j);
    area+=TriangleArea(x00,x01,x11)+TriangleArea(x00,x11,x10);
  }
  return area;
}

class GlobalSurfaceScope {
 public:
  GlobalSurfaceScope()
      : saved_faces_(CutCell::BoundaryTriangleFaces),
        saved_nodes_(CutCell::BoundaryTriangleNodes),
        saved_face_count_(CutCell::nBoundaryTriangleFaces),
        saved_node_count_(CutCell::nBoundaryTriangleNodes) {
    CutCell::BoundaryTriangleFaces=NULL;
    CutCell::BoundaryTriangleNodes=NULL;
    CutCell::nBoundaryTriangleFaces=0;
    CutCell::nBoundaryTriangleNodes=0;
  }

  ~GlobalSurfaceScope() {
    Reset();
    CutCell::BoundaryTriangleFaces=saved_faces_;
    CutCell::BoundaryTriangleNodes=saved_nodes_;
    CutCell::nBoundaryTriangleFaces=saved_face_count_;
    CutCell::nBoundaryTriangleNodes=saved_node_count_;
  }

  void Reset() {
    delete [] CutCell::BoundaryTriangleFaces;
    delete [] CutCell::BoundaryTriangleNodes;
    CutCell::BoundaryTriangleFaces=NULL;
    CutCell::BoundaryTriangleNodes=NULL;
    CutCell::nBoundaryTriangleFaces=0;
    CutCell::nBoundaryTriangleNodes=0;
  }

 private:
  CutCell::cTriangleFace* saved_faces_;
  CutCell::cNASTRANnode* saved_nodes_;
  int saved_face_count_,saved_node_count_;
};

std::string FindFixture(const char* relative) {
  const char* explicit_root=std::getenv("AMPS_TEST_REPO_ROOT");
  if (explicit_root!=NULL) {
    const std::string candidate=std::string(explicit_root)+"/"+relative;
    if (access(candidate.c_str(),R_OK)==0) return candidate;
  }

  char cwd[PATH_MAX];
  if (getcwd(cwd,sizeof(cwd))==NULL) return std::string();
  std::string root(cwd);
  for (int level=0;level<16;level++) {
    const std::string candidate=root+"/"+relative;
    if (access(candidate.c_str(),R_OK)==0) return candidate;
    root+="/..";
  }
  return std::string();
}

class TemporaryPath {
 public:
  TemporaryPath() {
    char pattern[]="/tmp/amps-cutcell-nastran-XXXXXX";
    const int descriptor=mkstemp(pattern);
    AMPS_EXPECT_TRUE(descriptor>=0);
    close(descriptor);
    path_=pattern;
  }
  ~TemporaryPath() { if (!path_.empty()) unlink(path_.c_str()); }
  const char* c_str() const { return path_.c_str(); }
 private:
  std::string path_;
};

struct FaceSnapshot {
  double x[3][3];
  double normal[3];
  int attribute;
};

std::vector<FaceSnapshot> SnapshotSurface() {
  std::vector<FaceSnapshot> result(CutCell::nBoundaryTriangleFaces);
  for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
    const double* coordinates[3] = {
        CutCell::BoundaryTriangleFaces[i].x0Face,
        CutCell::BoundaryTriangleFaces[i].x1Face,
        CutCell::BoundaryTriangleFaces[i].x2Face};
    for (int n=0;n<3;n++)
      std::memcpy(result[i].x[n],coordinates[n],3*sizeof(double));
    std::memcpy(result[i].normal,
                CutCell::BoundaryTriangleFaces[i].ExternalNormal,
                3*sizeof(double));
    result[i].attribute=CutCell::BoundaryTriangleFaces[i].attribute;
  }
  return result;
}

}  // namespace

AMPS_TEST_CASE("MESH-CUTCELL-TRIANGLE-METRICS-001",
               "triangle area, normal, winding, and edge length are analytic") {
  CutCell::cTriangleFace ccw,cw;
  SetTriangle(ccw,0.0,0.0,0.0,3.0,0.0,0.0,0.0,4.0,0.0);
  SetTriangle(cw,0.0,0.0,0.0,0.0,4.0,0.0,3.0,0.0,0.0);
  AMPS_EXPECT_NEAR(ccw.SurfaceArea,6.0,2.0e-14);
  AMPS_EXPECT_NEAR(cw.SurfaceArea,6.0,2.0e-14);
  AMPS_EXPECT_NEAR(ccw.ExternalNormal[0],0.0,2.0e-15);
  AMPS_EXPECT_NEAR(ccw.ExternalNormal[1],0.0,2.0e-15);
  AMPS_EXPECT_NEAR(ccw.ExternalNormal[2],1.0,2.0e-15);
  AMPS_EXPECT_NEAR(cw.ExternalNormal[2],-1.0,2.0e-15);

  // Translation-invariant Euclidean length for endpoints symmetric about the
  // origin is 2 m.  A sum-of-endpoints implementation incorrectly gives 0.
  CutCell::cNASTRANnode endpoint[2];
  endpoint[0].x[0]=-1.0;
  endpoint[1].x[0]=1.0;
  std::list<CutCell::cConnectivityListTriangleNode> nodes;
  nodes.push_back(CutCell::cConnectivityListTriangleNode());
  nodes.back().node=endpoint;
  nodes.push_back(CutCell::cConnectivityListTriangleNode());
  nodes.back().node=endpoint+1;
  CutCell::cConnectivityListTriangleEdge edge;
  std::list<CutCell::cConnectivityListTriangleNode>::iterator node_it=nodes.begin();
  edge.node[0]=node_it;
  ++node_it;
  edge.node[1]=node_it;
  AMPS_EXPECT_NEAR(edge.GetLength(),2.0,2.0e-15);
}

AMPS_TEST_CASE("MESH-CUTCELL-RAY-001",
               "ray hits, misses, origin tolerance, and backward planes are distinguished") {
  CutCell::cTriangleFace face;
  SetTriangle(face,0.0,0.0,0.0,2.0,0.0,0.0,0.0,2.0,0.0);
  double origin[3]={0.5,0.5,1.0},direction[3]={0.0,0.0,-2.0};
  double intersection[3],t;
  AMPS_EXPECT_TRUE(face.RayIntersection(origin,direction,t,intersection,0.0));
  AMPS_EXPECT_NEAR(t,0.5,2.0e-15);
  AMPS_EXPECT_NEAR(intersection[0],0.5,2.0e-15);
  AMPS_EXPECT_NEAR(intersection[1],0.5,2.0e-15);
  AMPS_EXPECT_NEAR(intersection[2],0.0,2.0e-15);

  double outside[3]={1.5,1.5,1.0};
  AMPS_EXPECT_FALSE(face.RayIntersection(outside,direction,t,0.0));
  double backward[3]={0.0,0.0,1.0};
  AMPS_EXPECT_FALSE(face.RayIntersection(origin,backward,t,0.0));

  double near_edge[3]={-5.0e-10,0.5,1.0};
  AMPS_EXPECT_FALSE(face.RayIntersection(near_edge,direction,t,0.0));
  AMPS_EXPECT_TRUE(face.RayIntersection(near_edge,direction,t,1.0e-9));

  // The ray begins 0.5 nm past the plane and travels away from it.  With no
  // tolerance the plane is behind the origin; a 1 nm spatial tolerance maps
  // that near-origin contact to t=0 without accepting more distant planes.
  double just_past[3]={0.5,0.5,-5.0e-10};
  AMPS_EXPECT_FALSE(face.RayIntersection(just_past,direction,t,0.0));
  AMPS_EXPECT_TRUE(face.RayIntersection(just_past,direction,t,1.0e-9));
  AMPS_EXPECT_NEAR(t,0.0,0.0);
}

AMPS_TEST_CASE("MESH-CUTCELL-TANGENT-DEGENERATE-001",
               "parallel, near-tangent, zero ray, and degenerate triangle are typed misses") {
  CutCell::cTriangleFace face,degenerate;
  SetTriangle(face,0.0,0.0,0.0,2.0,0.0,0.0,0.0,2.0,0.0);
  SetTriangle(degenerate,0.0,0.0,0.0,1.0,0.0,0.0,2.0,0.0,0.0);
  double origin[3]={0.25,0.25,1.0};
  double tangent[3]={1.0,0.0,0.0};
  double near_tangent[3]={1.0,0.0,-1.0e-16};
  double zero[3]={0.0,0.0,0.0};
  double down[3]={0.0,0.0,-1.0};
  double t=std::numeric_limits<double>::quiet_NaN();
  AMPS_EXPECT_FALSE(face.RayIntersection(origin,tangent,t,1.0e-12));
  AMPS_EXPECT_NEAR(t,0.0,0.0);
  AMPS_EXPECT_FALSE(face.RayIntersection(origin,near_tangent,t,1.0e-12));
  AMPS_EXPECT_FALSE(face.RayIntersection(origin,zero,t,1.0e-12));
  AMPS_EXPECT_FALSE(degenerate.RayIntersection(origin,down,t,1.0e-12));
  AMPS_EXPECT_NEAR(t,0.0,0.0);
}

AMPS_TEST_CASE("MESH-CUTCELL-CLOSEST-001",
               "closest forward surface hit is independent of face array order") {
  CutCell::cTriangleFace faces[2];
  SetTriangle(faces[0],0.0,0.0,2.0,2.0,0.0,2.0,0.0,2.0,2.0);
  SetTriangle(faces[1],0.0,0.0,0.75,2.0,0.0,0.75,0.0,2.0,0.75);
  double origin[3]={0.25,0.25,0.0},direction[3]={0.0,0.0,1.0};
  double intersection[3],t;
  CutCell::cTriangleFace* selected=NULL;
  AMPS_EXPECT_TRUE(CutCell::GetClosestSurfaceIntersectionPoint(
      origin,direction,intersection,t,selected,faces,2,0.0));
  AMPS_EXPECT_TRUE(selected==faces+1);
  AMPS_EXPECT_NEAR(t,0.75,2.0e-15);
  AMPS_EXPECT_NEAR(intersection[2],0.75,2.0e-15);
}

AMPS_TEST_CASE("MESH-CUTCELL-PLANE-VOLUME-001",
               "closed planar cuts match exact and requested-error prism volumes") {
  OrientedBox box(0.25,2.0,-2.0,2.0,-2.0,2.0);
  CutCell::cTriangleFaceDescriptor descriptors[12];
  for (int i=0;i<12;i++) {
    descriptors[i].TriangleFace=box.faces+i;
    descriptors[i].prev=(i==0) ? NULL : descriptors+i-1;
    descriptors[i].next=(i==11) ? NULL : descriptors+i+1;
  }
  double xmin[3]={-1.0,-1.0,-1.0};
  double xmax[3]={ 1.0, 1.0, 1.0};
  const double dyadic_volume=CutCell::GetRemainedBlockVolume(
      xmin,xmax,1.0e-12,1.0e-3,box.faces,12,descriptors);
  AMPS_EXPECT_NEAR(dyadic_volume,5.0,2.0e-10);

  // A non-dyadic plane cannot align with the bounded recursive grid.  Its
  // independent prism oracle must nevertheless be met within the requested
  // relative volume error, which also tests the level-selection formula.
  OrientedBox non_dyadic_box(0.30,2.0,-2.0,2.0,-2.0,2.0);
  for (int i=0;i<12;i++) descriptors[i].TriangleFace=non_dyadic_box.faces+i;
  const double non_dyadic_volume=CutCell::GetRemainedBlockVolume(
      xmin,xmax,1.0e-12,1.0e-3,non_dyadic_box.faces,12,descriptors);
  AMPS_EXPECT_NEAR(non_dyadic_volume,5.2,5.2e-3);
}

AMPS_TEST_CASE("MESH-CUTCELL-SPHERE-VOLUME-001",
               "translated scaled sphere volumes match exact sphere-box intersections") {
  cInternalSphericalData sphere;
  double center[3]={3.0,-2.0,5.0};
  const double radius=2.0;
  sphere.SetSphereGeometricalParameters(center,radius);
  int status=-1;

  double inside_min[3]={2.75,-2.25,4.75};
  double inside_max[3]={3.25,-1.75,5.25};
  AMPS_EXPECT_NEAR(sphere.GetRemainedBlockVolume(
      inside_min,inside_max,0.0,&status),0.0,2.0e-14);
  AMPS_EXPECT_EQ(status,_AMR_BLOCK_OUTSIDE_DOMAIN_);

  double outside_min[3]={6.0,-2.5,4.5};
  double outside_max[3]={7.0,-1.5,5.5};
  AMPS_EXPECT_NEAR(sphere.GetRemainedBlockVolume(
      outside_min,outside_max,0.0,&status),1.0,2.0e-14);
  AMPS_EXPECT_EQ(status,_AMR_BLOCK_INSIDE_DOMAIN_);

  // One positive octant of a radius-R cube contains pi R^3/6 of sphere.
  double octant_min[3]={3.0,-2.0,5.0};
  double octant_max[3]={5.0,0.0,7.0};
  const double octant_reference=radius*radius*radius*(1.0-kPi/6.0);
  AMPS_EXPECT_NEAR(sphere.GetRemainedBlockVolume(
      octant_min,octant_max,0.0,&status),octant_reference,2.0e-12);
  AMPS_EXPECT_EQ(status,_AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_);

  // A cube of half-width 2R encloses the entire sphere.
  double enclosing_min[3]={center[0]-2.0*radius,
                           center[1]-2.0*radius,
                           center[2]-2.0*radius};
  double enclosing_max[3]={center[0]+2.0*radius,
                           center[1]+2.0*radius,
                           center[2]+2.0*radius};
  const double enclosing_reference=std::pow(4.0*radius,3)-
      4.0*kPi*std::pow(radius,3)/3.0;
  AMPS_EXPECT_NEAR(sphere.GetRemainedBlockVolume(
      enclosing_min,enclosing_max,0.0,&status),enclosing_reference,2.0e-10);
  AMPS_EXPECT_EQ(status,_AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_);
}

AMPS_TEST_CASE("MESH-CUTCELL-SPHERE-AREA-CONVERGENCE-001",
               "planar sphere triangulation area converges at second order") {
  SphereGridScope state;
  cInternalSphericalData sphere;
  double center[3]={4.0,-3.0,2.0};
  const double radius=2.5;
  sphere.SetSphereGeometricalParameters(center,radius);
  const double exact=4.0*kPi*radius*radius;
  const double coarse=PlanarSphereArea(sphere,8,16);
  const double fine=PlanarSphereArea(sphere,16,32);
  const double coarse_error=std::fabs(coarse-exact);
  const double fine_error=std::fabs(fine-exact);
  AMPS_EXPECT_TRUE(fine_error<0.32*coarse_error);
  AMPS_EXPECT_TRUE(fine_error/exact<1.2e-2);
}

AMPS_TEST_CASE("MESH-CUTCELL-NASTRAN-ROUNDTRIP-001",
               "sparse out-of-order NASTRAN IDs and documented write contract round trip") {
  GlobalSurfaceScope global_surface;
  const std::string fixture=FindFixture(
      "test/srcCoreMeshTests/fixtures/sparse_out_of_order_tetrahedron.bdf");
  AMPS_EXPECT_TRUE(!fixture.empty());
  CutCell::ReadNastranSurfaceMeshLongFormat(fixture.c_str(),1.0);
  AMPS_EXPECT_EQ(CutCell::nBoundaryTriangleNodes,4);
  AMPS_EXPECT_EQ(CutCell::nBoundaryTriangleFaces,4);
  const int sparse_ids[4]={42,900,7,305};
  for (int i=0;i<4;i++)
    AMPS_EXPECT_EQ(CutCell::BoundaryTriangleNodes[i].id,sparse_ids[i]);
  AMPS_EXPECT_EQ(CutCell::BoundaryTriangleFaces[0].node[0]-
                 CutCell::BoundaryTriangleNodes,0);
  AMPS_EXPECT_EQ(CutCell::BoundaryTriangleFaces[0].node[1]-
                 CutCell::BoundaryTriangleNodes,3);
  AMPS_EXPECT_EQ(CutCell::BoundaryTriangleFaces[0].node[2]-
                 CutCell::BoundaryTriangleNodes,2);
  AMPS_EXPECT_NEAR(CutCell::BoundaryTriangleFaces[0].ExternalNormal[2],
                   -1.0,2.0e-15);
  const std::vector<FaceSnapshot> reference=SnapshotSurface();

  TemporaryPath serialized;
  CutCell::SaveNastranSurfaceMeshLongFormat(serialized.c_str());
  global_surface.Reset();
  CutCell::ReadNastranSurfaceMeshLongFormat(serialized.c_str(),1.0);
  AMPS_EXPECT_EQ(CutCell::nBoundaryTriangleNodes,4);
  AMPS_EXPECT_EQ(CutCell::nBoundaryTriangleFaces,4);
  const std::vector<FaceSnapshot> roundtrip=SnapshotSurface();
  AMPS_EXPECT_EQ(roundtrip.size(),reference.size());
  for (std::size_t i=0;i<reference.size();i++) {
    AMPS_EXPECT_EQ(roundtrip[i].attribute,reference[i].attribute);
    for (int node=0;node<3;node++) for (int d=0;d<3;d++)
      AMPS_EXPECT_NEAR(roundtrip[i].x[node][d],reference[i].x[node][d],0.0);
    for (int d=0;d<3;d++)
      AMPS_EXPECT_NEAR(roundtrip[i].normal[d],reference[i].normal[d],2.0e-15);
  }
}

#endif  // _MESH_DIMENSION_ == 3

void LinkMeshCutCellTests() {}

}  // namespace CoreMeshTests
}  // namespace AMPS

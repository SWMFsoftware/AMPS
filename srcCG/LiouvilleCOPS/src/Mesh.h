//closed triangulated surface of the nucleus
#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>

#include "Vec3.h"

struct cTriangle {
  int node[3];
  Vec3 Normal;       //unit external normal
  Vec3 Centroid;
  double Area;
  double MaxEdge;    //the longest edge length (used to choose the integration subdivision level)
};

class cMesh {
public:
  std::vector<Vec3> Nodes;
  std::vector<cTriangle> Faces;

  double TotalArea=0.0,Volume=0.0;

  //read the surface triangulation (CTRIA3 elements) of a NASTRAN bulk data file;
  //tetrahedra and all other element cards are ignored
  void ReadNastran(const std::string& fname);

  //read a SPICE DSK type 2 (plate model) file; vertices are converted from km to m; the DSK reference frame
  //must be 'expectedFrame'; the DSK must already be loadable (no kernels are required to read the plates)
  void ReadDSK(const std::string& fname,const std::string& expectedFrame);

  //generate a sphere (subdivided icosahedron) of radius R; used for the analytic verification
  void MakeIcosphere(double R,int nRefinement);

  //verify that the surface is closed and consistently oriented; orient all normals outward
  //requireClosed=false: open/non-manifold edges are reported but accepted (the orientation must still be consistent)
  void FinalizeGeometry(bool verbose,bool requireClosed=true);

  long nOpenEdges=0;   //edges not shared by exactly two triangles (0 for a closed surface)

  Vec3 Vertex(int iface,int i) const {return Nodes[Faces[iface].node[i]];}
};

//closest distance from point p to the triangle (a,b,c)
double PointTriangleDistance(const Vec3& p,const Vec3& a,const Vec3& b,const Vec3& c);

#endif

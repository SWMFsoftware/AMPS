/*
 * Triangulated test bodies whose geometry is known exactly, and the exact convex geometry used to check the
 * decomposition of cut cells into tetrahedra (test/srcCutCellTetraTest/main.cpp).
 *
 * A body is a closed triangulated surface (outward orientation) and, independently, the union of convex polyhedra
 * ('components', each the intersection of half-spaces n.x<=d) that the surface bounds exactly.  The components of
 * a body touch at most along faces, so vol(X ∩ body) = sum over the components of vol(X ∩ component).
 * (This code is local to the test; the shape generator of the mover tests will be merged with it.)
 */
#ifndef _CUTCELL_TETRA_TEST_SHAPES_H_
#define _CUTCELL_TETRA_TEST_SHAPES_H_

#include <string>
#include <vector>

namespace CutCellTetraTest {

struct cPlane {double n[3],d;};                          //the half-space n.x<=d (|n|=1)
struct cComponent {                                      //a convex polyhedron
  std::vector<cPlane> planes;
  bool hasBox;                                           //a bounding box of the component (to skip it quickly)
  double lo[3],hi[3];
  cComponent() : hasBox(false) {}
};

struct cAnalyticSphere {double c[3],r;};                 //an analytic sphere (an internal boundary of the mesh)

struct cBody {
  std::string name;
  std::vector<double> vertex;                            //x,y,z of the vertices
  std::vector<int> face;                                 //three vertices per face, counterclockwise seen from outside
  std::vector<cComponent> components;
  double volume;                                         //the exact volume of the polyhedron

  //analytic spheres: the body is the union of the spheres (disjoint), registered as internal spheres of the mesh; the
  //surface and the components above are then their triangulated approximations of the cut-cell decomposition
  //(CutCell::cSphereApproximation, BuildSphereApproximation())
  std::vector<cAnalyticSphere> spheres;

  double Phi(const double *x) const;                     //>0 outside; the distance to the nearest component face plane (approximately)
};

//"sphere1", "sphere2", "sphere3": icospheres (80/320/1280 faces); "cone", "cylinder": 32 segments; "box": axis-aligned, its faces
//lie on planes of the mesh; "lshape": the union of two boxes (an L-shaped prism), rotated.
//Multiple crossings of cell edges (the mesh: cells of 0.05 at the surface for x>0, their faces on the planes x,y,z=0.05*k):
//"plate": a tilted plate thinner than a cell (an edge enters and leaves the body between its corners);
//"bumps": three small icospheres: one whose cap pokes through a cell edge (the edge is crossed twice, its corners are
//outside), one through a cell face only (no edge is crossed), one inside a cell;
//"lcorner": an L-shaped prism rotated by 45 degrees about z, its concave edge inside a cell column: x-edges near the
//concave edge cross both faces that meet there;
//"wedge": a prism whose sharp edge lies along cell edges, with its end vertices on cell corners and its caps on cell faces.
//Analytic spheres (internal boundaries, no triangulated surface): "asphere": one sphere (radius 0.3, center off the origin);
//"aspheres": two spheres with a gap smaller than a cell (an edge crosses both) and a small sphere through a cell edge.
bool BuildBody(const std::string& name,cBody& body);

//the triangulated approximations of the analytic spheres of the body (N[i]: the resolution of sphere i used by the
//mesh, cMeshAMRgeneric::GetSphereApproximationResolution()): the 'cube sphere' of CutCell::cSphereApproximation,
//computed here independently; the components are the tetrahedra (center, facet)
void BuildSphereApproximation(cBody& body,const std::vector<int>& N);
void WriteNastranFile(const cBody& body,const char *fname);

//---------------------------------- exact convex geometry ----------------------------------
typedef std::vector<std::vector<double> > cPolygon3D;    //the vertices of a planar convex polygon
struct cPolyhedron {std::vector<cPolygon3D> faces;};     //a convex polyhedron

cPolyhedron Box(const double *xmin,const double *xmax);
cPolyhedron Tetrahedron(const double x[4][3]);
cPolyhedron Clip(const cPolyhedron& p,const cPlane& h);  //p ∩ {n.x<=d}
double Volume(const cPolyhedron& p);
cPolygon3D ClipPolygon(const cPolygon3D& p,const cPlane& h);
double Area(const cPolygon3D& p);

//the body restricted to the box [lo,hi] (for the exact geometry inside the box): the components that do not reach the
//box are dropped, and the faces that do not cut the box (the box is inside their half-space) are dropped
cBody RestrictToBox(const cBody& body,const double *lo,const double *hi);

//vol(p ∩ body) and area(polygon ∩ body)
double IntersectionVolume(const cPolyhedron& p,const cBody& body);
double IntersectionArea(const cPolygon3D& p,const cBody& body);

}

#endif

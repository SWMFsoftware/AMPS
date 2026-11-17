/*
 * AnalyticPolyhedra.h: closed triangulated bodies with exact inside tests (header only)
 *
 * A body (cBody) is a closed, consistently oriented triangulation (each face is counter-clockwise seen from
 * outside: (v1-v0)x(v2-v0) points out of the solid) together with the solid it bounds, given as a union of
 * convex pieces.  The inside test is a test of the triangulated polyhedron itself, not of the smooth shape it
 * approximates:
 *   - a convex body has one piece whose planes are the planes of its own triangles (computed from the
 *     final vertex coordinates): a point is inside when it is on the inner side of every face plane;
 *   - a non-convex body is a union of (overlapping) convex boxes whose faces carry its triangles: a point is
 *     inside when it is inside any piece.  The pieces overlap so that every interior point of the union is
 *     strictly inside at least one piece.
 * Depth(x) is a lower bound of the distance of x from the surface when x is inside (>0), and <=0 outside or on
 * the surface; Inside(x,tol) is Depth(x)>tol.
 *
 * Generators (local coordinates, then Transform() with a rotation and a shift): Icosphere (level 0..5), Cone
 * (sharp tip), Cylinder, Box (subdivided faces), voxel unions (an L-shaped body, a hollow cube with a closed
 * cavity).  The voxel unions have their vertices on an integer grid, so faces that share an edge share its
 * vertices exactly and there are no T-junctions.
 *
 * cScene collects several disjoint bodies into one surface and writes it in the NASTRAN long format that
 * CutCell::ReadNastranSurfaceMeshLongFormat() reads (GRID* and CTRIA3 records, the vertices with %.17E so that they are read
 * back bit for bit).
 *
 * Used by the test of the force-free cut-surface mover (test/srcCutSurfaceMoverTest); it has no AMPS dependency.
 */

#ifndef _MESHAMR_ANALYTIC_POLYHEDRA_H_
#define _MESHAMR_ANALYTIC_POLYHEDRA_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace AnalyticPolyhedra {

inline void Cross(const double *a,const double *b,double *c) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

inline double Dot(const double *a,const double *b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}

//a plane n.x=d with the unit normal n pointing out of the solid
struct cPlane {double n[3],d;};

//a convex piece of the solid: the intersection of the inner half-spaces of its planes
struct cConvexPiece {
  std::vector<cPlane> Plane;
  double xmin[3],xmax[3];   //bounding box (of its vertices)

  //min over the planes of the distance to the plane on the inner side (>0: strictly inside)
  double Depth(const double *x) const {
    for (int i=0;i<3;i++) if ((x[i]<=xmin[i])||(x[i]>=xmax[i])) return -1.0;

    double depth=HUGE_VAL;
    for (size_t n=0;n<Plane.size();n++) {
      double s=Plane[n].d-Dot(Plane[n].n,x);
      if (s<depth) {
        depth=s;
        if (depth<=0.0) return depth;
      }
    }
    return depth;
  }
};

class cBody {
public:
  std::string Name;
  std::vector<double> Vertex;        //3 coordinates per vertex
  std::vector<int> Face;             //3 vertex indices per face (counter-clockwise seen from outside)
  std::vector<cConvexPiece> Piece;   //the solid is the union of the pieces
  bool Convex;                       //one piece made of the planes of the faces
  double xmin[3],xmax[3];

  //for the pieces of the non-convex bodies: their corners in local coordinates (box pieces), and the
  //accumulated transformation x -> RTotal x + ShiftTotal from the local coordinates
  std::vector<std::vector<double> > PieceBoxLocal;  //6 numbers per box: xmin[3],xmax[3]
  double RTotal[3][3],ShiftTotal[3];

  cBody() : Convex(true) {
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) RTotal[i][j]=(i==j) ? 1.0 : 0.0;
    for (int i=0;i<3;i++) ShiftTotal[i]=0.0;
  }

  int nVertices() const {return Vertex.size()/3;}
  int nFaces() const {return Face.size()/3;}
  const double* V(int i) const {return &Vertex[3*i];}

  void FaceNormal(int f,double *n) const {
    const double *a=V(Face[3*f]),*b=V(Face[3*f+1]),*c=V(Face[3*f+2]);
    double e0[3],e1[3];
    for (int i=0;i<3;i++) e0[i]=b[i]-a[i],e1[i]=c[i]-a[i];
    Cross(e0,e1,n);
    double l=std::sqrt(Dot(n,n));
    for (int i=0;i<3;i++) n[i]/=l;
  }

  double FaceArea(int f) const {
    const double *a=V(Face[3*f]),*b=V(Face[3*f+1]),*c=V(Face[3*f+2]);
    double e0[3],e1[3],n[3];
    for (int i=0;i<3;i++) e0[i]=b[i]-a[i],e1[i]=c[i]-a[i];
    Cross(e0,e1,n);
    return 0.5*std::sqrt(Dot(n,n));
  }

  double Area() const {
    double s=0.0;
    for (int f=0;f<nFaces();f++) s+=FaceArea(f);
    return s;
  }

  //volume of the polyhedron (divergence theorem, relative to its first vertex)
  double Volume() const {
    double s=0.0;
    const double *o=V(0);
    for (int f=0;f<nFaces();f++) {
      const double *a=V(Face[3*f]),*b=V(Face[3*f+1]),*c=V(Face[3*f+2]);
      double ra[3],rb[3],rc[3],t[3];
      for (int i=0;i<3;i++) ra[i]=a[i]-o[i],rb[i]=b[i]-o[i],rc[i]=c[i]-o[i];
      Cross(rb,rc,t);
      s+=Dot(ra,t)/6.0;
    }
    return s;
  }

  //>0: inside the solid (a lower bound of the distance from the surface); <=0: outside or on the surface
  double Depth(const double *x) const {
    for (int i=0;i<3;i++) if ((x[i]<=xmin[i])||(x[i]>=xmax[i])) return -1.0;

    double depth=-HUGE_VAL;
    for (size_t n=0;n<Piece.size();n++) {
      double d=Piece[n].Depth(x);
      if (d>depth) depth=d;
      if (depth>0.0) break;
    }
    return depth;
  }

  bool Inside(const double *x,double tol=0.0) const {return Depth(x)>tol;}

  //update the bounding boxes and, for a convex body, the planes of the faces from the present vertices
  void Finalize() {
    for (int i=0;i<3;i++) xmin[i]=HUGE_VAL,xmax[i]=-HUGE_VAL;
    for (int v=0;v<nVertices();v++) for (int i=0;i<3;i++) xmin[i]=std::min(xmin[i],V(v)[i]),xmax[i]=std::max(xmax[i],V(v)[i]);

    if (Convex) {
      cConvexPiece p;
      for (int f=0;f<nFaces();f++) {
        cPlane pl;
        FaceNormal(f,pl.n);
        pl.d=Dot(pl.n,V(Face[3*f]));
        p.Plane.push_back(pl);
      }
      for (int i=0;i<3;i++) p.xmin[i]=xmin[i],p.xmax[i]=xmax[i];
      Piece.assign(1,p);
    }
  }

  //x -> R x + shift for the vertices; the box pieces of a non-convex body are built from the local boxes with
  //the accumulated transformation
  void Transform(const double R[3][3],const double *shift) {
    for (int v=0;v<nVertices();v++) {
      double x[3];
      for (int i=0;i<3;i++) x[i]=R[i][0]*Vertex[3*v]+R[i][1]*Vertex[3*v+1]+R[i][2]*Vertex[3*v+2]+shift[i];
      for (int i=0;i<3;i++) Vertex[3*v+i]=x[i];
    }

    {
      double T[3][3],s[3];
      for (int i=0;i<3;i++) {
        s[i]=shift[i];
        for (int j=0;j<3;j++) {
          T[i][j]=0.0;
          for (int k=0;k<3;k++) T[i][j]+=R[i][k]*RTotal[k][j];
          s[i]+=R[i][j]*ShiftTotal[j];
        }
      }
      for (int i=0;i<3;i++) {ShiftTotal[i]=s[i]; for (int j=0;j<3;j++) RTotal[i][j]=T[i][j];}
    }

    if (!Convex) {
      Piece.clear();
      for (size_t b=0;b<PieceBoxLocal.size();b++) {
        const std::vector<double>& box=PieceBoxLocal[b];
        cConvexPiece p;

        for (int i=0;i<3;i++) p.xmin[i]=HUGE_VAL,p.xmax[i]=-HUGE_VAL;
        for (int c=0;c<8;c++) {
          double xl[3]={box[(c&1) ? 3 : 0],box[(c&2) ? 4 : 1],box[(c&4) ? 5 : 2]},x[3];
          for (int i=0;i<3;i++) x[i]=RTotal[i][0]*xl[0]+RTotal[i][1]*xl[1]+RTotal[i][2]*xl[2]+ShiftTotal[i];
          for (int i=0;i<3;i++) p.xmin[i]=std::min(p.xmin[i],x[i]),p.xmax[i]=std::max(p.xmax[i],x[i]);
        }

        for (int axis=0;axis<3;axis++) for (int side=0;side<2;side++) {
          double nl[3]={0.0,0.0,0.0},pl[3]={box[0],box[1],box[2]};
          nl[axis]=(side==0) ? -1.0 : 1.0;
          pl[axis]=box[axis+3*side];
          cPlane plane;
          double xp[3];
          for (int i=0;i<3;i++) {
            plane.n[i]=RTotal[i][0]*nl[0]+RTotal[i][1]*nl[1]+RTotal[i][2]*nl[2];
            xp[i]=RTotal[i][0]*pl[0]+RTotal[i][1]*pl[1]+RTotal[i][2]*pl[2]+ShiftTotal[i];
          }
          double l=std::sqrt(Dot(plane.n,plane.n));
          for (int i=0;i<3;i++) plane.n[i]/=l;
          plane.d=Dot(plane.n,xp);
          p.Plane.push_back(plane);
        }
        Piece.push_back(p);
      }
    }

    Finalize();
  }

  //every edge is shared by exactly two faces that use it in opposite directions; no degenerate faces
  bool CheckClosedManifold(std::string& message) const {
    std::map<std::pair<int,int>,int> edge;
    char buf[400];

    for (int f=0;f<nFaces();f++) {
      for (int e=0;e<3;e++) {
        int a=Face[3*f+e],b=Face[3*f+(e+1)%3];
        if (a==b) {sprintf(buf,"%s: face %i is degenerate",Name.c_str(),f); message=buf; return false;}
        edge[std::make_pair(a,b)]++;
      }
      if (FaceArea(f)<=0.0) {sprintf(buf,"%s: face %i has no area",Name.c_str(),f); message=buf; return false;}
    }

    for (std::map<std::pair<int,int>,int>::const_iterator it=edge.begin();it!=edge.end();it++) {
      std::map<std::pair<int,int>,int>::const_iterator op=edge.find(std::make_pair(it->first.second,it->first.first));
      if ((it->second!=1)||(op==edge.end())||(op->second!=1)) {
        sprintf(buf,"%s: edge %i-%i is used %i times, the opposite edge %i times",Name.c_str(),it->first.first,it->first.second,
            it->second,(op==edge.end()) ? 0 : op->second);
        message=buf;
        return false;
      }
    }

    int chi=nVertices()-(int)edge.size()/2+nFaces();
    sprintf(buf,"%s: %i vertices, %i edges, %i faces, Euler characteristic %i",Name.c_str(),nVertices(),(int)edge.size()/2,nFaces(),chi);
    message=buf;
    return true;
  }
};

//----------------------------------- generators (local coordinates) -----------------------------------

//icosphere: vertices on the sphere of radius r centered at the origin; 20*4^level faces
inline cBody Icosphere(int level,double r) {
  cBody b;
  const double t=(1.0+std::sqrt(5.0))/2.0;
  double v0[12][3]={{-1,t,0},{1,t,0},{-1,-t,0},{1,-t,0},{0,-1,t},{0,1,t},{0,-1,-t},{0,1,-t},{t,0,-1},{t,0,1},{-t,0,-1},{-t,0,1}};
  int f0[20][3]={{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
                 {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};
  std::vector<double> v;
  std::vector<int> f;

  for (int i=0;i<12;i++) {
    double l=std::sqrt(v0[i][0]*v0[i][0]+v0[i][1]*v0[i][1]+v0[i][2]*v0[i][2]);
    for (int d=0;d<3;d++) v.push_back(v0[i][d]/l);
  }
  for (int i=0;i<20;i++) for (int d=0;d<3;d++) f.push_back(f0[i][d]);

  for (int s=0;s<level;s++) {
    std::map<std::pair<int,int>,int> middle;
    std::vector<int> fnew;

    auto Middle=[&](int a,int c) {
      std::pair<int,int> key(std::min(a,c),std::max(a,c));
      std::map<std::pair<int,int>,int>::iterator it=middle.find(key);
      if (it!=middle.end()) return it->second;
      double m[3],l=0.0;
      for (int d=0;d<3;d++) m[d]=0.5*(v[3*a+d]+v[3*c+d]),l+=m[d]*m[d];
      l=std::sqrt(l);
      for (int d=0;d<3;d++) v.push_back(m[d]/l);
      int id=v.size()/3-1;
      middle[key]=id;
      return id;
    };

    for (size_t n=0;n<f.size();n+=3) {
      int a=f[n],bb=f[n+1],c=f[n+2];
      int ab=Middle(a,bb),bc=Middle(bb,c),ca=Middle(c,a);
      int add[12]={a,ab,ca, bb,bc,ab, c,ca,bc, ab,bc,ca};
      fnew.insert(fnew.end(),add,add+12);
    }
    f.swap(fnew);
  }

  b.Vertex.resize(v.size());
  for (size_t n=0;n<v.size();n++) b.Vertex[n]=r*v[n];
  b.Face=f;

  //orientation: outward
  for (int n=0;n<b.nFaces();n++) {
    double nn[3],c[3];
    b.FaceNormal(n,nn);
    for (int d=0;d<3;d++) c[d]=b.V(b.Face[3*n])[d];
    if (Dot(nn,c)<0.0) std::swap(b.Face[3*n+1],b.Face[3*n+2]);
  }

  char name[100];
  sprintf(name,"icosphere%i",level);
  b.Name=name;
  b.Convex=true;
  b.Finalize();
  return b;
}

//cone: base circle (n sides, radius r) in the plane z=0 around the origin, apex (sharp tip) at (0,0,h)
inline cBody Cone(int n,double r,double h) {
  cBody b;
  b.Vertex.insert(b.Vertex.end(),{0.0,0.0,h});      //0: apex
  b.Vertex.insert(b.Vertex.end(),{0.0,0.0,0.0});    //1: base center
  for (int i=0;i<n;i++) {
    double a=2.0*M_PI*i/n;
    b.Vertex.insert(b.Vertex.end(),{r*std::cos(a),r*std::sin(a),0.0});
  }
  for (int i=0;i<n;i++) {
    int p=2+i,q=2+(i+1)%n;
    b.Face.insert(b.Face.end(),{p,q,0});   //side
    b.Face.insert(b.Face.end(),{1,q,p});   //base
  }
  b.Name="cone";
  b.Convex=true;
  b.Finalize();
  return b;
}

//cylinder: axis z, radius r, 0<=z<=h, n sides, nz rings of side faces
inline cBody Cylinder(int n,double r,double h,int nz) {
  cBody b;
  for (int k=0;k<=nz;k++) for (int i=0;i<n;i++) {
    double a=2.0*M_PI*i/n;
    b.Vertex.insert(b.Vertex.end(),{r*std::cos(a),r*std::sin(a),h*k/nz});
  }
  int c0=b.nVertices(); b.Vertex.insert(b.Vertex.end(),{0.0,0.0,0.0});
  int c1=b.nVertices(); b.Vertex.insert(b.Vertex.end(),{0.0,0.0,h});

  auto id=[&](int k,int i) {return k*n+(i%n);};
  for (int k=0;k<nz;k++) for (int i=0;i<n;i++) {
    int a=id(k,i),bb=id(k,i+1),c=id(k+1,i+1),d=id(k+1,i);
    b.Face.insert(b.Face.end(),{a,bb,c});
    b.Face.insert(b.Face.end(),{a,c,d});
  }
  for (int i=0;i<n;i++) {
    b.Face.insert(b.Face.end(),{c0,id(0,i+1),id(0,i)});
    b.Face.insert(b.Face.end(),{c1,id(nz,i),id(nz,i+1)});
  }
  b.Name="cylinder";
  b.Convex=true;
  b.Finalize();
  return b;
}

//a union of unit voxels (cell (i,j,k) is the cube [i,i+1]x[j,j+1]x[k,k+1] times 'size'); every voxel face
//on the boundary is split into m x m squares, two triangles each; 'pieces' are the overlapping boxes (in
//voxel units, 6 numbers each) whose union is the solid
inline cBody VoxelUnion(const std::vector<std::tuple<int,int,int> >& voxels,double size,int m,
    const std::vector<std::vector<double> >& pieces,const char *name) {
  cBody b;
  std::map<std::tuple<int,int,int>,int> filled,vertexId;
  for (size_t n=0;n<voxels.size();n++) filled[voxels[n]]=1;

  auto Vertex=[&](int i,int j,int k) {     //grid point in units of size/m
    std::tuple<int,int,int> key(i,j,k);
    std::map<std::tuple<int,int,int>,int>::iterator it=vertexId.find(key);
    if (it!=vertexId.end()) return it->second;
    int id=b.nVertices();
    b.Vertex.insert(b.Vertex.end(),{size*i/m,size*j/m,size*k/m});
    vertexId[key]=id;
    return id;
  };

  for (size_t n=0;n<voxels.size();n++) {
    int vi=std::get<0>(voxels[n]),vj=std::get<1>(voxels[n]),vk=std::get<2>(voxels[n]);

    for (int axis=0;axis<3;axis++) for (int side=0;side<2;side++) {
      int nb[3]={vi,vj,vk};
      nb[axis]+=(side==0) ? -1 : 1;
      if (filled.count(std::make_tuple(nb[0],nb[1],nb[2]))) continue;

      int a1=(axis+1)%3,a2=(axis+2)%3;
      for (int p=0;p<m;p++) for (int q=0;q<m;q++) {
        int c[4][3];
        int base[3]={vi*m,vj*m,vk*m};
        int du[4]={0,1,1,0},dv[4]={0,0,1,1};
        for (int s=0;s<4;s++) {
          c[s][axis]=base[axis]+side*m;
          c[s][a1]=base[a1]+p+du[s];
          c[s][a2]=base[a2]+q+dv[s];
        }
        int id[4];
        for (int s=0;s<4;s++) id[s]=Vertex(c[s][0],c[s][1],c[s][2]);
        //(a1,a2,axis) is right-handed: the quad 0-1-2-3 is counter-clockwise seen from +axis
        if (side==1) {
          b.Face.insert(b.Face.end(),{id[0],id[1],id[2]});
          b.Face.insert(b.Face.end(),{id[0],id[2],id[3]});
        }
        else {
          b.Face.insert(b.Face.end(),{id[0],id[2],id[1]});
          b.Face.insert(b.Face.end(),{id[0],id[3],id[2]});
        }
      }
    }
  }

  for (size_t n=0;n<pieces.size();n++) {
    std::vector<double> box(6);
    for (int i=0;i<6;i++) box[i]=size*pieces[n][i];
    b.PieceBoxLocal.push_back(box);
  }

  b.Name=name;
  b.Convex=false;
  double R[3][3]={{1,0,0},{0,1,0},{0,0,1}},zero[3]={0,0,0};
  b.Transform(R,zero);   //builds the pieces
  return b;
}

//an L-shaped prism: three voxels (0,0,0),(1,0,0),(0,1,0); the pieces [0,2]x[0,1]x[0,1] and [0,1]x[0,2]x[0,1] overlap
inline cBody LShape(double size,int m) {
  std::vector<std::tuple<int,int,int> > v={std::make_tuple(0,0,0),std::make_tuple(1,0,0),std::make_tuple(0,1,0)};
  std::vector<std::vector<double> > p={{0,0,0,2,1,1},{0,0,0,1,2,1}};
  return VoxelUnion(v,size,m,p,"L-shape");
}

//a hollow cube: 3x3x3 voxels without the central one (a closed cavity); the pieces are the six slabs of
//thickness one voxel at the faces of the cube
inline cBody HollowCube(double size,int m) {
  std::vector<std::tuple<int,int,int> > v;
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) if (!((i==1)&&(j==1)&&(k==1))) v.push_back(std::make_tuple(i,j,k));
  std::vector<std::vector<double> > p={{0,0,0,1,3,3},{2,0,0,3,3,3},{0,0,0,3,1,3},{0,2,0,3,3,3},{0,0,0,3,3,1},{0,0,2,3,3,3}};
  return VoxelUnion(v,size,m,p,"hollow-cube");
}

//box [0,a]x[0,b]x[0,c] with every face split into m x m rectangles
inline cBody Box(double a,double bb,double c,int m) {
  std::vector<std::tuple<int,int,int> > v={std::make_tuple(0,0,0)};
  std::vector<std::vector<double> > p={{0,0,0,1,1,1}};
  cBody b=VoxelUnion(v,1.0,m,p,"box");
  for (int n=0;n<b.nVertices();n++) b.Vertex[3*n]*=a,b.Vertex[3*n+1]*=bb,b.Vertex[3*n+2]*=c;
  b.PieceBoxLocal.clear();
  b.Convex=true;
  b.Finalize();
  return b;
}

//a closed shell between the boxes [-a,a]^3 and [-b,b]^3 (a<b): the outer surface points out, the inner one into
//the box [-a,a]^3 (a container for a gas); every face of both boxes is split into m x m squares; the pieces are the
//six slabs of thickness b-a
inline cBody HollowBox(double a,double b,int m) {
  cBody outer=Box(2.0*b,2.0*b,2.0*b,m),inner=Box(2.0*a,2.0*a,2.0*a,m),s;
  int offset=outer.nVertices();

  for (int v=0;v<outer.nVertices();v++) for (int i=0;i<3;i++) s.Vertex.push_back(outer.Vertex[3*v+i]-b);
  for (int v=0;v<inner.nVertices();v++) for (int i=0;i<3;i++) s.Vertex.push_back(inner.Vertex[3*v+i]-a);
  s.Face=outer.Face;
  for (int f=0;f<inner.nFaces();f++) s.Face.insert(s.Face.end(),{inner.Face[3*f]+offset,inner.Face[3*f+2]+offset,inner.Face[3*f+1]+offset});

  for (int axis=0;axis<3;axis++) for (int side=0;side<2;side++) {
    std::vector<double> box={-b,-b,-b,b,b,b};
    if (side==0) box[axis+3]=-a; else box[axis]=a;
    s.PieceBoxLocal.push_back(box);
  }

  s.Name="shell";
  s.Convex=false;
  double R[3][3]={{1,0,0},{0,1,0},{0,0,1}},zero[3]={0,0,0};
  s.Transform(R,zero);
  return s;
}

//rotation matrix from Euler angles (z-x-z)
inline void RotationMatrix(double alpha,double beta,double gamma,double R[3][3]) {
  double ca=std::cos(alpha),sa=std::sin(alpha),cb=std::cos(beta),sb=std::sin(beta),cg=std::cos(gamma),sg=std::sin(gamma);
  double A[3][3]={{ca,-sa,0},{sa,ca,0},{0,0,1}},B[3][3]={{1,0,0},{0,cb,-sb},{0,sb,cb}},G[3][3]={{cg,-sg,0},{sg,cg,0},{0,0,1}},T[3][3];
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {T[i][j]=0.0; for (int k=0;k<3;k++) T[i][j]+=A[i][k]*B[k][j];}
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {R[i][j]=0.0; for (int k=0;k<3;k++) R[i][j]+=T[i][k]*G[k][j];}
}

//------------------------------------------- scene -------------------------------------------
//several disjoint bodies; the surface is the union of their triangulations
class cScene {
public:
  std::vector<cBody> Body;

  void Add(const cBody& b) {Body.push_back(b);}

  int nFaces() const {int n=0; for (size_t i=0;i<Body.size();i++) n+=Body[i].nFaces(); return n;}

  //the body with the largest depth of x and that depth (>0: x is inside that body)
  double Depth(const double *x,int *body=NULL) const {
    double depth=-HUGE_VAL;
    if (body!=NULL) *body=-1;
    for (size_t i=0;i<Body.size();i++) {
      double d=Body[i].Depth(x);
      if (d>depth) {depth=d; if (body!=NULL) *body=i;}
    }
    return depth;
  }

  bool Inside(const double *x,double tol=0.0) const {return Depth(x)>tol;}

  double Area() const {double s=0.0; for (size_t i=0;i<Body.size();i++) s+=Body[i].Area(); return s;}
  double Volume() const {double s=0.0; for (size_t i=0;i<Body.size();i++) s+=Body[i].Volume(); return s;}

  //NASTRAN long format (GRID*, CTRIA3), the vertices numbered over all bodies from 1
  bool WriteNastran(const char *fname,const char *title) const {
    FILE *f=fopen(fname,"w");
    if (f==NULL) return false;

    //one title line: the reader expects the GRID records from the second line on (no comments between them);
    //the body is the property field of CTRIA3 (the face attribute): 1, 2, ...
    std::string names;
    for (size_t i=0;i<Body.size();i++) names+=std::string((i==0) ? "" : ",")+Body[i].Name;
    fprintf(f,"$ %s: bodies %s\n",title,names.c_str());
    int id=1,offset=0;
    for (size_t i=0;i<Body.size();i++) {
      for (int v=0;v<Body[i].nVertices();v++,id++) {
        const double *x=Body[i].V(v);
        fprintf(f,"GRID* %i 0 %.17E %.17E\n* %.17E\n",id,x[0],x[1],x[2]);
      }
    }
    id=1;
    for (size_t i=0;i<Body.size();i++) {
      for (int n=0;n<Body[i].nFaces();n++,id++) {
        fprintf(f,"CTRIA3 %i %i %i %i %i\n",id,(int)i+1,Body[i].Face[3*n]+1+offset,Body[i].Face[3*n+1]+1+offset,Body[i].Face[3*n+2]+1+offset);
      }
      offset+=Body[i].nVertices();
    }
    fprintf(f,"ENDDATA\n");
    fclose(f);
    return true;
  }
};

}  // namespace AnalyticPolyhedra

#endif

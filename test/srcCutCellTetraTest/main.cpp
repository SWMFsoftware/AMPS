/*
 * Test of the decomposition of cut cells into tetrahedra (CutCell::GetCutCellTetrahedra(),
 * cMeshAMRgeneric::GetCutcellTetrahedronMesh() and the tetrahedral zones of the volume output).
 *
 * Usage: amps sphere1|sphere2|sphere3|cone|cylinder|box|lshape|plate|bumps|lcorner|wedge|asphere|aspheres
 *        amps nastran <file> <scale>   (a diagnostic of any closed surface, see below)
 *   The body (shapes.cpp) is a closed triangulated surface whose geometry is known exactly: icospheres (80, 320 and
 *   1280 faces), a cone and a cylinder (32 segments, caps of coplanar faces), an axis-aligned box whose faces and
 *   edges lie on planes, edges and corners of cells, and the L-shaped union of two boxes (non-convex), rotated.
 *   Bodies built so that cell edges are crossed more than once (the cells at the surface are 0.05 for x>0): a tilted
 *   plate thinner than a cell; small icospheres whose caps poke through one cell edge or one cell face only, and one
 *   inside a cell; an L-shaped prism whose concave edge lies inside a cell, so that edges near it cross both faces of
 *   the notch; a wedge whose sharp edge lies along cell edges with its vertices on cell corners (degenerate cases).
 *   Analytic spheres (internal spheres of the mesh, no triangulated surface): one sphere off the origin; two spheres
 *   with a gap smaller than a cell, so that an edge crosses both, and a small one through a cell edge.  Their exact
 *   geometry in the checks is their triangulated approximation that the decomposition uses (computed here
 *   independently), and the "analytic" check compares the outside volume of every cell with the analytic one
 *   (cInternalSphericalData::GetRemainedBlockVolume()): the remaining approximation error, at most
 *   CutCell::SphereApproximationTolerance (1e-2) V_cell per sphere that crosses the cell.
 *   The domain is [-1,1]^3; the mesh is refined at the surface, twice as fine for x>0 as for x<0, so the surface
 *   crosses jumps of the resolution, and it is distributed over the processes (blocks of other processes).
 *   For every cut cell the crossings of its 12 edges with the surface are counted (INFO line of the report: cells with
 *   an edge crossed 2+ times, with no edge crossed, with a face crossed but none of its edges, with a corner on the
 *   surface, with an edge in the surface, with an edge crossed by two spheres); the cases built for a crossing must
 *   produce it ("crossings" check).
 *
 * For every interior cell of the blocks of this process, with the exact volume Vb=vol(cell ∩ body) (clipping of the
 * cell by the faces of the convex components of the body):
 *   outside    every tetrahedron of the cell in the output (GetCutcellTetrahedronMesh(): the outside tetrahedra of a
 *              cut cell, the full cell when it is not cut, nothing when it is inside) lies outside the body:
 *              vol(T ∩ body) <= 1e-9 V_cell;
 *   inside     every inside tetrahedron of GetCutCellTetrahedra() lies inside: vol(T)-vol(T ∩ body) <= 1e-9 V_cell;
 *   overlap    the tetrahedra of a cell (outside and inside) do not overlap (pairwise vol(T1 ∩ T2) <= 1e-9 V_cell)
 *              and fill the cell (sum = V_cell to 1e-9);
 *   volume     the outside tetrahedra sum to the exact outside volume V_cell-Vb within 1e-9 V_cell;
 *   measure    ... and agree with the cell measure of InitCellMeasure() within its accuracy 1e-2 V_cell (not with
 *              two spheres in a cell: InitCellMeasure() permits one internal boundary per cell);
 *   faces      on each of the six faces of the cell the outside tetrahedra cover exactly the part of the face
 *              outside the body (area to 1e-9 of the face area); as this holds for every cell, the traces seen from
 *              the two cells that share a face (of any size, in any block or process) are the same set;
 *   file       the tetrahedral zones of the volume output file (OutputDistributedDataTECPLOT()) hold exactly the
 *              outside tetrahedra of the checks: the same number, the same vertices (to the printed precision),
 *              valid connectivity, and the same number of values per point as the variable list.
 * One line "PASS|FAIL <check> ..." each is written to cutcell-tetra.<shape>.np<N>.report.
 * cutcell-tetra.<shape>.np<N>.cells lists the cut cells in a canonical order (level, position, number of
 * tetrahedra, volumes): the files of runs with different numbers of processes must be identical.
 * AMPS_CUTCELL_TETRA_LEGACY=1 tests the old procedure of GetCutcellTetrahedronMesh() instead (no inside tetrahedra).
 *
 * "amps nastran <file> <scale>" (not part of the nightly test) reads a closed surface from a NASTRAN file (coordinates
 * multiplied by <scale>), whose geometry is not known in closed form: the tetrahedra of each cut cell must fill it,
 * and the volume of the tetrahedra whose centroid is on the wrong side of the surface (the parity of the crossings of
 * three rays with all faces) must not exceed 1e-6 V_cell.  AMPS_CUTCELL_TETRA_PROBE="xmin ymin zmin xmax ymax zmax px py pz"
 * additionally decomposes that box with all faces of the surface and reports the side of the point p (to examine a cell
 * of another mesh, e.g. of an application run).
 */
#include "pic.h"
#include "shapes.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>
#include <type_traits>
#include <string>
#include <vector>

using namespace CutCellTetraTest;

namespace {

typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cNode;
typedef std::remove_pointer<decltype(PIC::Mesh::mesh)>::type cMesh;

//the return codes of GetCutcellTetrahedronMesh()
const int StatusRemoved=0,StatusCut=1,StatusComplete=2;

cBody Body;
FILE *Report=NULL;

//"nastran <file> <scale>": a surface read from a file (diagnostic mode, not in the nightly test). Its geometry is not
//known in closed form: a tetrahedron is checked by the side of the surface of its centroid (parity of the crossings
//of three rays with all faces, a majority vote)
bool General=false;
std::vector<double> GTri;
double GLo[3],GHi[3],GH=0.0;

bool GeneralInside(const double *x) {
  static const double dir[3][3]={{0.5773502692,0.5773502692,0.5773502692},{0.8017837257,-0.2672612419,0.5345224838},{-0.3015113446,0.9045340337,0.3015113446}};
  int votes=0;
  for (int r=0;r<3;r++) {
    const double *d=dir[r];
    int c=0;
    for (size_t n=0;n<GTri.size();n+=9) {
      const double *a=&GTri[n],*b=&GTri[n+3],*q=&GTri[n+6];
      double e1[3],e2[3],h[3],sv[3],qq[3],det,u,v,t;
      for (int k=0;k<3;k++) e1[k]=b[k]-a[k],e2[k]=q[k]-a[k];
      h[0]=d[1]*e2[2]-d[2]*e2[1],h[1]=d[2]*e2[0]-d[0]*e2[2],h[2]=d[0]*e2[1]-d[1]*e2[0];
      det=e1[0]*h[0]+e1[1]*h[1]+e1[2]*h[2];
      if (fabs(det)<1.0E-300) continue;
      for (int k=0;k<3;k++) sv[k]=x[k]-a[k];
      u=(sv[0]*h[0]+sv[1]*h[1]+sv[2]*h[2])/det;
      if ((u<0.0)||(u>1.0)) continue;
      qq[0]=sv[1]*e1[2]-sv[2]*e1[1],qq[1]=sv[2]*e1[0]-sv[0]*e1[2],qq[2]=sv[0]*e1[1]-sv[1]*e1[0];
      v=(d[0]*qq[0]+d[1]*qq[1]+d[2]*qq[2])/det;
      if ((v<0.0)||(u+v>1.0)) continue;
      t=(e2[0]*qq[0]+e2[1]*qq[1]+e2[2]*qq[2])/det;
      if (t>0.0) c++;
    }
    if (c%2==1) votes++;
  }
  return votes>=2;
}
int nFailed=0;
double SurfaceCellSize=0.12;

//analytic spheres ("asphere", "aspheres"): the internal spheres of the mesh
bool Analytic=false;
std::vector<cInternalSphericalData*> Spheres;

//the crossings of the cell edges with the surface (the triangulated surface of the body; for analytic spheres, their
//triangulated approximation): the triangles with their bounding boxes and the part of the body they belong to
struct cSurfaceTriangle {double v[3][3],lo[3],hi[3]; int part;};
std::vector<cSurfaceTriangle> SurfaceTriangles;

void InitSurfaceTriangles() {
  SurfaceTriangles.clear();

  for (size_t n=0;n<Body.face.size();n+=3) {
    cSurfaceTriangle t;
    double c[3]={0.0,0.0,0.0};

    for (int k=0;k<3;k++) for (int d=0;d<3;d++) t.v[k][d]=Body.vertex[3*Body.face[n+k]+d],c[d]+=t.v[k][d]/3.0;
    for (int d=0;d<3;d++) t.lo[d]=std::min(t.v[0][d],std::min(t.v[1][d],t.v[2][d])),t.hi[d]=std::max(t.v[0][d],std::max(t.v[1][d],t.v[2][d]));

    //the part: the sphere whose surface is the nearest (analytic spheres), otherwise 0
    t.part=0;
    double best=1.0E30;
    for (size_t ns=0;ns<Body.spheres.size();ns++) {
      double r[3]={c[0]-Body.spheres[ns].c[0],c[1]-Body.spheres[ns].c[1],c[2]-Body.spheres[ns].c[2]};
      double dist=fabs(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])-Body.spheres[ns].r);
      if (dist<best) best=dist,t.part=ns;
    }

    SurfaceTriangles.push_back(t);
  }
}

//the crossings of the 12 edges of the cell [lo,hi]: nCross[dir][s1][s2] the number of distinct points where the surface
//crosses the edge along 'dir' at the sides s1, s2 of the axes (dir+1)%3, (dir+2)%3; 'corner': the surface passes
//through a corner of the cell; 'along': an edge of the cell lies in the surface; 'parts': an edge is crossed by two parts
//of the body (two spheres)
void EdgeCrossings(const double *lo,const double *hi,int nCross[3][2][2],bool& corner,bool& along,bool& parts) {
  double h=0.0;
  for (int d=0;d<3;d++) h=std::max(h,hi[d]-lo[d]);
  const double tol=1.0E-9*h;

  corner=false,along=false,parts=false;

  for (int dir=0;dir<3;dir++) for (int s1=0;s1<2;s1++) for (int s2=0;s2<2;s2++) {
    int a1=(dir+1)%3,a2=(dir+2)%3;
    double p0[3],p1[3],elo[3],ehi[3];

    p0[dir]=lo[dir],p1[dir]=hi[dir];
    p0[a1]=p1[a1]=s1 ? hi[a1] : lo[a1];
    p0[a2]=p1[a2]=s2 ? hi[a2] : lo[a2];
    for (int d=0;d<3;d++) elo[d]=std::min(p0[d],p1[d])-tol,ehi[d]=std::max(p0[d],p1[d])+tol;

    std::vector<std::pair<double,int> > t;

    for (size_t n=0;n<SurfaceTriangles.size();n++) {
      const cSurfaceTriangle& f=SurfaceTriangles[n];
      bool overlap=true;
      for (int d=0;d<3;d++) if ((f.lo[d]>ehi[d])||(f.hi[d]<elo[d])) overlap=false;
      if (!overlap) continue;

      double e1[3],e2[3],nn[3],l,d0,d1;
      for (int d=0;d<3;d++) e1[d]=f.v[1][d]-f.v[0][d],e2[d]=f.v[2][d]-f.v[0][d];
      nn[0]=e1[1]*e2[2]-e1[2]*e2[1],nn[1]=e1[2]*e2[0]-e1[0]*e2[2],nn[2]=e1[0]*e2[1]-e1[1]*e2[0];
      l=sqrt(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
      if (!(l>0.0)) continue;
      for (int d=0;d<3;d++) nn[d]/=l;
      d0=nn[0]*(p0[0]-f.v[0][0])+nn[1]*(p0[1]-f.v[0][1])+nn[2]*(p0[2]-f.v[0][2]);
      d1=nn[0]*(p1[0]-f.v[0][0])+nn[1]*(p1[1]-f.v[0][1])+nn[2]*(p1[2]-f.v[0][2]);

      //is the point q of the plane inside the triangle (to the tolerance)
      auto InTriangle=[&](const double *q) {
        for (int k=0;k<3;k++) {
          const double *a=f.v[k],*b=f.v[(k+1)%3];
          double ed[3],r[3],cr[3],len;
          for (int d=0;d<3;d++) ed[d]=b[d]-a[d],r[d]=q[d]-a[d];
          cr[0]=ed[1]*r[2]-ed[2]*r[1],cr[1]=ed[2]*r[0]-ed[0]*r[2],cr[2]=ed[0]*r[1]-ed[1]*r[0];
          len=sqrt(ed[0]*ed[0]+ed[1]*ed[1]+ed[2]*ed[2]);
          if (cr[0]*nn[0]+cr[1]*nn[1]+cr[2]*nn[2]<-tol*len) return false;
        }
        return true;
      };

      if ((fabs(d0)<=tol)&&(fabs(d1)<=tol)) {
        //the edge lies in the plane of the triangle: along the surface when its middle is in the triangle
        double m[3];
        for (int d=0;d<3;d++) m[d]=0.5*(p0[d]+p1[d]);
        if (InTriangle(m)) along=true;
        continue;
      }

      if (((fabs(d0)<=tol)&&InTriangle(p0))||((fabs(d1)<=tol)&&InTriangle(p1))) corner=true;
      if (!(((d0<-tol)&&(d1>tol))||((d0>tol)&&(d1<-tol)))) continue;

      double s=d0/(d0-d1),q[3];
      for (int d=0;d<3;d++) q[d]=p0[d]+s*(p1[d]-p0[d]);
      if (InTriangle(q)) t.push_back(std::make_pair(s,f.part));
    }

    //distinct crossing points (a point on an edge of the surface is found in both triangles that share it)
    std::sort(t.begin(),t.end());
    int cnt=0;
    for (size_t n=0;n<t.size();n++) {
      if ((n>0)&&(t[n].first-t[n-1].first<=1.0E-9)) continue;
      cnt++;
      if ((n>0)&&(t[n].second!=t[0].second)) parts=true;
    }
    nCross[dir][s1][s2]=cnt;
  }
}

//the volume of the cell [lo,hi] inside the analytic spheres (disjoint): cInternalSphericalData::GetRemainedBlockVolume()
double AnalyticInsideVolume(double *lo,double *hi,double eps) {
  double V=(hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2]),v=0.0;

  for (size_t ns=0;ns<Spheres.size();ns++) {
    int status;
    v+=V-Spheres[ns]->GetRemainedBlockVolume(lo,hi,eps,&status);
  }
  return v;
}

void Check(bool pass,const char *name,const char *format,...) {
  char message[2000];
  va_list args;

  va_start(args,format);
  vsnprintf(message,sizeof(message),format,args);
  va_end(args);

  if (PIC::ThisThread==0) {
    fprintf(Report,"%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
    printf("%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
    fflush(Report);
  }

  if (!pass) nFailed++;
}

//the mesh: fine at the surface, twice as fine for x>0; coarse elsewhere
double Resolution(double *x) {
  if (General) {
    for (int d=0;d<3;d++) if ((x[d]<GLo[d]-2.0*GH)||(x[d]>GHi[d]+2.0*GH)) return 8.0*GH;
    return GH;
  }

  double h=(x[0]>0.0) ? SurfaceCellSize : 2.0*SurfaceCellSize;
  if (fabs(Body.Phi(x))<2.0*h) return h;
  return 0.4;
}

double SphereResolution(double *x) {return (x[0]>0.0) ? SurfaceCellSize : 2.0*SurfaceCellSize;}

double SurfaceResolution(CutCell::cTriangleFace *f) {
  if (General) return GH;
  double c[3];
  for (int d=0;d<3;d++) c[d]=(f->x0Face[d]+f->x1Face[d]+f->x2Face[d])/3.0;
  return (c[0]>0.0) ? SurfaceCellSize : 2.0*SurfaceCellSize;
}

double UnitLoad(cNode*) {return 1.0;}

double TetVolume(const double x[4][3]) {
  double e1[3],e2[3],e3[3];
  for (int d=0;d<3;d++) e1[d]=x[1][d]-x[0][d],e2[d]=x[2][d]-x[0][d],e3[d]=x[3][d]-x[0][d];
  return (e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])-e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])+e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]))/6.0;
}

void TetBox(const double x[4][3],double *lo,double *hi) {
  for (int d=0;d<3;d++) {
    lo[d]=hi[d]=x[0][d];
    for (int i=1;i<4;i++) lo[d]=std::min(lo[d],x[i][d]),hi[d]=std::max(hi[d],x[i][d]);
  }
}

//is one of the face planes of a or b a separating plane of the two tetrahedra (to the tolerance)
bool TetSeparated(const double a[4][3],const double b[4][3],double tol) {
  static const int fc[4][4]={{0,1,2,3},{0,1,3,2},{0,2,3,1},{1,2,3,0}};

  for (int t=0;t<2;t++) {
    const double (*p)[3]=(t==0) ? a : b;
    const double (*q)[3]=(t==0) ? b : a;

    for (int f=0;f<4;f++) {
      const double *p0=p[fc[f][0]],*p1=p[fc[f][1]],*p2=p[fc[f][2]],*r=p[fc[f][3]];
      double e1[3],e2[3],n[3],l,sr;
      for (int d=0;d<3;d++) e1[d]=p1[d]-p0[d],e2[d]=p2[d]-p0[d];
      n[0]=e1[1]*e2[2]-e1[2]*e2[1],n[1]=e1[2]*e2[0]-e1[0]*e2[2],n[2]=e1[0]*e2[1]-e1[1]*e2[0];
      l=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      if (!(l>0.0)) continue;
      sr=((r[0]-p0[0])*n[0]+(r[1]-p0[1])*n[1]+(r[2]-p0[2])*n[2])/l;   //the side of the tetrahedron p

      bool sep=true;
      for (int v=0;(v<4)&&sep;v++) {
        double s=((q[v][0]-p0[0])*n[0]+(q[v][1]-p0[1])*n[1]+(q[v][2]-p0[2])*n[2])/l;
        if (s*((sr>0.0) ? 1.0 : -1.0)>tol) sep=false;
      }
      if (sep) return true;
    }
  }
  return false;
}

//the volume of the intersection of two tetrahedra: one clipped by the four face planes of the other
double TetIntersection(const double a[4][3],const double b[4][3]) {
  cPolyhedron p=Tetrahedron(a);
  static const int fc[4][4]={{0,1,2,3},{0,1,3,2},{0,2,3,1},{1,2,3,0}};   //three vertices of a face, the fourth

  for (int f=0;(f<4)&&!p.faces.empty();f++) {
    const double *p0=b[fc[f][0]],*p1=b[fc[f][1]],*p2=b[fc[f][2]],*q=b[fc[f][3]];
    double e1[3],e2[3],n[3],l;
    cPlane h;

    for (int d=0;d<3;d++) e1[d]=p1[d]-p0[d],e2[d]=p2[d]-p0[d];
    n[0]=e1[1]*e2[2]-e1[2]*e2[1],n[1]=e1[2]*e2[0]-e1[0]*e2[2],n[2]=e1[0]*e2[1]-e1[1]*e2[0];
    l=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    if ((q[0]-p0[0])*n[0]+(q[1]-p0[1])*n[1]+(q[2]-p0[2])*n[2]>0.0) l=-l;   //outward: the fourth vertex is inside
    for (int d=0;d<3;d++) h.n[d]=n[d]/l;
    h.d=h.n[0]*p0[0]+h.n[1]*p0[1]+h.n[2]*p0[2];
    p=Clip(p,h);
  }

  return Volume(p);
}

//the cut cells for the canonical list
struct cCellRecord {
  long int key[4];
  int nOut,nIn,status;
  double vOut,vIn,vExact,measure;
  bool operator<(const cCellRecord& c) const {return std::lexicographical_compare(key,key+4,c.key,c.key+4);}
};

//a tetrahedron in the canonical form: its vertices in lexicographic order
struct cTetKey {
  double x[12];
  bool operator<(const cTetKey& t) const {return std::lexicographical_compare(x,x+12,t.x,t.x+12);}
};

cTetKey MakeTetKey(const double x[4][3]) {
  std::vector<std::vector<double> > v(4,std::vector<double>(3));
  cTetKey k;
  for (int i=0;i<4;i++) for (int d=0;d<3;d++) v[i][d]=x[i][d];
  std::sort(v.begin(),v.end());
  for (int i=0;i<4;i++) for (int d=0;d<3;d++) k.x[3*i+d]=v[i][d];
  return k;
}

//read the tetrahedral zones of a TECPLOT file written by OutputDistributedDataTECPLOT()
bool ReadTetraZones(const char *fname,std::vector<cTetKey>& tets,int& nVariables,long& nBadColumns,long& nBadConnectivity) {
  std::ifstream in(fname);
  std::string line;

  tets.clear(),nVariables=0,nBadColumns=0,nBadConnectivity=0;
  if (!in) return false;

  while (std::getline(in,line)) {
    if (line.compare(0,9,"VARIABLES")==0) {
      nVariables=std::count(line.begin(),line.end(),',')+1;
      continue;
    }
    if ((line.compare(0,4,"ZONE")!=0)||(line.find("TETRAHEDRON")==std::string::npos)) continue;

    long N=0,E=0;
    sscanf(strstr(line.c_str(),"N="),"N=%ld",&N);
    sscanf(strstr(line.c_str(),"E="),"E=%ld",&E);

    std::vector<double> x(3*N);
    for (long n=0;n<N;n++) {
      std::getline(in,line);
      std::istringstream s(line);
      double v;
      int cnt=0;
      while (s>>v) {if (cnt<3) x[3*n+cnt]=v; cnt++;}
      if (cnt!=nVariables) nBadColumns++;
    }

    for (long e=0;e<E;e++) {
      long id[4];
      std::getline(in,line);
      if (sscanf(line.c_str(),"%ld %ld %ld %ld",id,id+1,id+2,id+3)!=4) {nBadConnectivity++; continue;}
      bool ok=true;
      for (int i=0;i<4;i++) if ((id[i]<1)||(id[i]>N)) ok=false;
      for (int i=0;i<4;i++) for (int j=i+1;j<4;j++) if (id[i]==id[j]) ok=false;
      if (!ok) {nBadConnectivity++; continue;}
      double t[4][3];
      for (int i=0;i<4;i++) for (int d=0;d<3;d++) t[i][d]=x[3*(id[i]-1)+d];
      tets.push_back(MakeTetKey(t));
    }
  }

  return true;
}

}  // namespace


int main(int argc,char **argv) {
  PIC::Init_BeforeParser();
  PIC::Mesh::initCellSamplingDataBuffer();

  std::string shape=(argc>1) ? argv[1] : "sphere1";
  General=(shape=="nastran");
  if ((General==false)&&(!BuildBody(shape,Body))) {
    if (PIC::ThisThread==0) printf("usage: %s sphere1|sphere2|sphere3|cone|cylinder|box|lshape|plate|bumps|lcorner|wedge|asphere|aspheres | nastran <file> <scale>\n",argv[0]);
    MPI_Finalize();
    return 2;
  }

  const bool legacy=(getenv("AMPS_CUTCELL_TETRA_LEGACY")!=NULL)&&(atoi(getenv("AMPS_CUTCELL_TETRA_LEGACY"))!=0);
  char fname[400];

  //the body
  double xmin[3]={-1.0,-1.0,-1.0},xmax[3]={1.0,1.0,1.0};

  if (General) {
    if (argc<4) {
      if (PIC::ThisThread==0) printf("usage: %s nastran <file> <scale>\n",argv[0]);
      MPI_Finalize();
      return 2;
    }
    CutCell::ReadNastranSurfaceMeshLongFormat(argv[2],atof(argv[3]));
    for (int d=0;d<3;d++) GLo[d]=1.0E300,GHi[d]=-1.0E300;
    for (int n=0;n<CutCell::nBoundaryTriangleFaces;n++) {
      CutCell::cTriangleFace *f=CutCell::BoundaryTriangleFaces+n;
      double *v[3]={f->x0Face,f->x1Face,f->x2Face};
      for (int k=0;k<3;k++) for (int d=0;d<3;d++) {GTri.push_back(v[k][d]); GLo[d]=std::min(GLo[d],v[k][d]),GHi[d]=std::max(GHi[d],v[k][d]);}
    }
    double L=0.0,c[3];
    for (int d=0;d<3;d++) L=std::max(L,GHi[d]-GLo[d]),c[d]=0.5*(GLo[d]+GHi[d]);
    GH=L/24.0;
    for (int d=0;d<3;d++) xmin[d]=c[d]-0.8*L,xmax[d]=c[d]+0.8*L;
    shape=std::string("nastran-")+std::string(strrchr(argv[2],'/') ? strrchr(argv[2],'/')+1 : argv[2]);
  }
  else if (Body.spheres.empty()==false) {
    //analytic spheres: internal spheres of the mesh (no triangulated surface)
    Analytic=true;
    PIC::BC::InternalBoundary::Sphere::Init();

    for (size_t ns=0;ns<Body.spheres.size();ns++) {
      cInternalBoundaryConditionsDescriptor descriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
      cInternalSphericalData *sphere=(cInternalSphericalData*)descriptor.BoundaryElement;
      double x0[3]={Body.spheres[ns].c[0],Body.spheres[ns].c[1],Body.spheres[ns].c[2]};

      sphere->SetSphereGeometricalParameters(x0,Body.spheres[ns].r);
      sphere->localResolution=SphereResolution;
      Spheres.push_back(sphere);
    }
  }
  else {
    sprintf(fname,"cutcell-tetra.%s.bdf",shape.c_str());
    if (PIC::ThisThread==0) WriteNastranFile(Body,fname);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    CutCell::ReadNastranSurfaceMeshLongFormat(fname);
  }

  if (Analytic==false) {
    PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
    PIC::Mesh::mesh->CutCellSurfaceLocalResolution=SurfaceResolution;
  }

  sprintf(fname,"cutcell-tetra.%s%s.np%i.report",shape.c_str(),legacy ? ".legacy" : "",PIC::nTotalThreads);
  if (PIC::ThisThread==0) Report=fopen(fname,"w");

  //the mesh: generated by process 0 and read from the file by the others, as in the applications

  double tStart=MPI_Wtime(),tMesh,tMeasure,tCells,tOutput;
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,Resolution);

  sprintf(fname,"cutcell-tetra.%s.mesh",shape.c_str());
  if (PIC::ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile(fname);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  else {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    PIC::Mesh::mesh->readMeshFile(fname);
  }

  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  tMesh=MPI_Wtime();
  //the orientation of the face normals (ray tracing, as in the applications)
  if (Analytic==false) PIC::Mesh::IrregularSurface::InitExternalNormalVector();

  //the cell measures (a cache of an earlier run is removed). InitCellMeasure() permits one internal boundary per cell: with
  //two spheres in a cell ("aspheres") it stops, so the cell measure is not checked there (the analytic check replaces it)
  const bool MeasureAvailable=(shape!="aspheres");
  if (MeasureAvailable) {
    char cache[600];
    PIC::Mesh::mesh->GetCenterNodeMeasureFileName(cache,NULL,PIC::Mesh::mesh->getMeshSignature());
    if (PIC::ThisThread==0) remove(cache);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    PIC::Mesh::mesh->InitCellMeasure();
  }

  tMeasure=MPI_Wtime();
  //the external normals point out of the body: the bodies of the test are wound counterclockwise seen from outside, so the
  //normal of the vertex order is outward; a surface read from a file must have the same sign everywhere (it is closed and
  //consistently wound) and the outward one (the volume enclosed by the normals is positive)
  if (Analytic) {
    //the triangulated approximations of the spheres that the mesh uses (the resolution from the finest cell size at each
    //sphere): the exact geometry of the checks below, and the crossings of the cell edges
    std::vector<int> N;
    double maxGap=0.0;

    for (size_t ns=0;ns<Spheres.size();ns++) N.push_back(PIC::Mesh::mesh->GetSphereApproximationResolution(Spheres[ns]));
    BuildSphereApproximation(Body,N);

    //the largest distance between the sphere and its approximation (at the centroids of the facets, where it is the
    //largest), relative to the finest cell size of the leaf blocks that the spheres cross
    double hmin=1.0E30;
    std::function<void(cNode*)> scan=[&] (cNode* nd) {
      bool crosses=false;
      for (size_t ns=0;ns<Body.spheres.size();ns++) {
        double dmin2=0.0,dmax2=0.0;
        for (int d=0;d<3;d++) {
          double a=nd->xmin[d]-Body.spheres[ns].c[d],b=nd->xmax[d]-Body.spheres[ns].c[d];
          if (a>0.0) dmin2+=a*a; else if (b<0.0) dmin2+=b*b;
          dmax2+=std::max(a*a,b*b);
        }
        if ((dmin2<=pow(Body.spheres[ns].r,2))&&(dmax2>=pow(Body.spheres[ns].r,2))) crosses=true;
      }
      if (!crosses) return;
      if (nd->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {hmin=std::min(hmin,(nd->xmax[0]-nd->xmin[0])/_BLOCK_CELLS_X_); return;}
      for (int i=0;i<8;i++) if (nd->downNode[i]!=NULL) scan(nd->downNode[i]);
    };
    scan(PIC::Mesh::mesh->rootTree);

    for (size_t n=0;n<Body.face.size();n+=3) {
      double c[3]={0.0,0.0,0.0},best=1.0E30;
      for (int k=0;k<3;k++) for (int d=0;d<3;d++) c[d]+=Body.vertex[3*Body.face[n+k]+d]/3.0;
      for (size_t ns=0;ns<Body.spheres.size();ns++) {
        double r[3]={c[0]-Body.spheres[ns].c[0],c[1]-Body.spheres[ns].c[1],c[2]-Body.spheres[ns].c[2]};
        best=std::min(best,fabs(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])-Body.spheres[ns].r));
      }
      maxGap=std::max(maxGap,best);
    }

    if (PIC::ThisThread==0) {
      double vAnalytic=0.0;
      fprintf(Report,"INFO analytic spheres:");
      for (size_t ns=0;ns<Spheres.size();ns++) {
        fprintf(Report," [center %g %g %g, radius %g, N=%i]",Body.spheres[ns].c[0],Body.spheres[ns].c[1],Body.spheres[ns].c[2],Body.spheres[ns].r,N[ns]);
        vAnalytic+=4.0/3.0*Pi*pow(Body.spheres[ns].r,3);
      }
      fprintf(Report,"; %li facets, the volume of the approximation %.9e (spheres %.9e), the largest distance between a sphere and its facets %.3e (%.3e of the finest cell size %g at the spheres; tolerance %.1e)\n",
          (long)Body.face.size()/3,Body.volume,vAnalytic,maxGap,maxGap/hmin,hmin,CutCell::SphereApproximationTolerance);
    }
  }
  else {
    long nPositive=0,nNegative=0;
    double volume=0.0;

    for (int n=0;n<CutCell::nBoundaryTriangleFaces;n++) {
      CutCell::cTriangleFace *f=CutCell::BoundaryTriangleFaces+n;
      double e1[3],e2[3],w[3];
      for (int d=0;d<3;d++) e1[d]=f->x1Face[d]-f->x0Face[d],e2[d]=f->x2Face[d]-f->x0Face[d];
      w[0]=e1[1]*e2[2]-e1[2]*e2[1],w[1]=e1[2]*e2[0]-e1[0]*e2[2],w[2]=e1[0]*e2[1]-e1[1]*e2[0];
      double dot=w[0]*f->ExternalNormal[0]+w[1]*f->ExternalNormal[1]+w[2]*f->ExternalNormal[2],len=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
      (dot>0.0) ? nPositive++ : nNegative++;
      //the volume by the divergence theorem with the external normals: sum over the faces of area*(x0.n)/3
      volume+=0.5*len*(f->x0Face[0]*f->ExternalNormal[0]+f->x0Face[1]*f->ExternalNormal[1]+f->x0Face[2]*f->ExternalNormal[2])/3.0;
    }

    long bad=General ? std::min(nPositive,nNegative) : nNegative;
    Check((bad==0)&&(volume>0.0),"normals","%li of the %i external normals of the surface are inverted; the volume enclosed by the normals is %.12e%s",
        bad,CutCell::nBoundaryTriangleFaces,volume,General ? "" : (std::string(" (exact ")+std::to_string(Body.volume)+")").c_str());
  }


  //AMPS_CUTCELL_TETRA_PROBE="xmin ymin zmin xmax ymax zmax px py pz" (general mode): the decomposition of that box with all
  //faces of the surface as candidates, and the side of the point p in it (a diagnostic of a cell of another mesh)
  if (General&&(getenv("AMPS_CUTCELL_TETRA_PROBE")!=NULL)&&(PIC::ThisThread==0)) {
    double b[9];
    if (sscanf(getenv("AMPS_CUTCELL_TETRA_PROBE"),"%lf %lf %lf %lf %lf %lf %lf %lf %lf",b,b+1,b+2,b+3,b+4,b+5,b+6,b+7,b+8)==9) {
      std::vector<CutCell::cTriangleFaceDescriptor> d(CutCell::nBoundaryTriangleFaces);
      for (int n=0;n<CutCell::nBoundaryTriangleFaces;n++) d[n].TriangleFace=CutCell::BoundaryTriangleFaces+n,d[n].next=(n+1<CutCell::nBoundaryTriangleFaces) ? &d[n+1] : NULL;
      std::vector<CutCell::cCutCellTetrahedron> out,in;
      double vo,vi;
      int st=CutCell::GetCutCellTetrahedra(b,b+3,&d[0],PIC::Mesh::mesh->EPS,&out,&in,&vo,&vi);
      int inOut=0,inIn=0;
      for (int side=0;side<2;side++) {
        std::vector<CutCell::cCutCellTetrahedron>& t=(side==0) ? out : in;
        for (size_t n=0;n<t.size();n++) {
          double m[3][3],r[3],det,bc[4];
          for (int k=0;k<3;k++) m[k][0]=t[n].x[1][k]-t[n].x[0][k],m[k][1]=t[n].x[2][k]-t[n].x[0][k],m[k][2]=t[n].x[3][k]-t[n].x[0][k],r[k]=b[6+k]-t[n].x[0][k];
          det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
          for (int c=0;c<3;c++) {
            double a[3][3];
            for (int i=0;i<3;i++) for (int j=0;j<3;j++) a[i][j]=(j==c) ? r[i] : m[i][j];
            bc[c+1]=(a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))/det;
          }
          bc[0]=1.0-bc[1]-bc[2]-bc[3];
          if ((bc[0]>=0)&&(bc[1]>=0)&&(bc[2]>=0)&&(bc[3]>=0)) (side==0) ? inOut++ : inIn++;
        }
      }
      printf("cutcell-tetra probe: status %i, outside %.6e inside %.6e of the box volume, the point is in %i outside and %i inside tetrahedra; parity: %s\n",
          st,vo/((b[3]-b[0])*(b[4]-b[1])*(b[5]-b[2])),vi/((b[3]-b[0])*(b[4]-b[1])*(b[5]-b[2])),inOut,inIn,GeneralInside(b+6) ? "inside" : "outside");
    }
  }

  //the checks of the cells of this process
  long nCells=0,nCut=0,nTetOut=0,nTetIn=0;
  long badOutside=0,badInside=0,badOverlap=0,badFill=0,badVolume=0,badMeasure=0,badFaces=0;
  double maxOutside=0.0,maxInside=0.0,maxOverlap=0.0,maxFill=0.0,maxVolume=0.0,maxMeasure=0.0,maxFaces=0.0;
  long badClass=0,badClassMeasure=0,badBrick=0,badRemoved=0,badBox=0,badFaces2=0,nRefCut=0,nOutCut=0;   //(general mode)
  long nEdge2=0,nNoEdge=0,nCorner=0,nAlong=0,nParts=0,nPoke=0,badAnalytic=0;   //the crossings of the cut cells; analytic spheres
  double maxAnalytic=0.0;
  double maxBrick=0.0,maxRemoved=0.0,maxBox=0.0;
  std::vector<cCellRecord> records;
  std::vector<double> outputTets;   //the outside tetrahedra of the cut cells (12 doubles each)
  const double eps=PIC::Mesh::mesh->EPS;

  if (!General) InitSurfaceTriangles();

  for (cNode* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    if ((node->block==NULL)||(node->IsUsedInCalculationFlag==false)) continue;

    double dx[3];
    for (int d=0;d<3;d++) dx[d]=(node->xmax[d]-node->xmin[d])/((d==0) ? _BLOCK_CELLS_X_ : ((d==1) ? _BLOCK_CELLS_Y_ : _BLOCK_CELLS_Z_));

    for (int k=0;k<_BLOCK_CELLS_Z_;k++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
      int idx[3]={i,j,k};
      double lo[3],hi[3],V=1.0;

      if (cell==NULL) continue;
      for (int d=0;d<3;d++) lo[d]=node->xmin[d]+idx[d]*dx[d],hi[d]=node->xmin[d]+(idx[d]+1)*dx[d],V*=dx[d];
      nCells++;

      //the general mode. The reference geometry of the cell is the exact decomposition GetCutCellTetrahedra() (its own
      //checks: the tetrahedra fill the cell and lie on their side of the surface by the parity test); the output of
      //GetCutcellTetrahedronMesh() (the new procedure, or the old one with AMPS_CUTCELL_TETRA_LEGACY=1) is checked against it
      if (General) {
        std::vector<CutCell::cCutCellTetrahedron> out,in;
        double vo=0.0,vi=0.0;
        int st=CutCell::GetCutCellTetrahedra(lo,hi,node->FirstTriangleCutFace,eps,&out,&in,&vo,&vi);
        bool cut=(vo>1.0E-6*V)&&(vi>1.0E-6*V);

        //the output
        list<cMesh::cTetrahedron> output;
        int status=PIC::Mesh::mesh->GetCutcellTetrahedronMesh(output,i,j,k,node);
        std::vector<CutCell::cCutCellTetrahedron> outT;
        for (list<cMesh::cTetrahedron>::iterator it=output.begin();it!=output.end();it++) {
          CutCell::cCutCellTetrahedron t;
          for (int v=0;v<4;v++) for (int d=0;d<3;d++) t.x[v][d]=it->x[v][d];
          outT.push_back(t);
        }
        double vOut=(status==StatusComplete) ? V : 0.0;
        if (status==StatusCut) for (size_t n=0;n<outT.size();n++) vOut+=fabs(TetVolume(outT[n].x));

        if ((st!=CutCell::_cut_cell__cut_)&&(status!=StatusCut)) {
          //a cell that is not cut: written in full when it is outside, not at all when it is inside
          bool outside=(st==CutCell::_cut_cell__outside_);
          if (outside!=(status==StatusComplete)) {
            //(the cell measure of InitCellMeasure() is an independent classification of the cell)
            if (outside==(cell->Measure>0.5*V)) badClassMeasure++;
            if (badClass<5) printf("cutcell-tetra: cell [%e %e %e]-[%e %e %e]: not cut, %s the body, written as %s\n",lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],
                outside ? "outside" : "inside",(status==StatusComplete) ? "a full cell" : "nothing");
            badClass++;
          }
          continue;
        }

        nCut++,nTetOut+=outT.size(),nTetIn+=in.size();
        if (cut) nRefCut++;
        if (status==StatusCut) nOutCut++;

        //the reference decomposition: its inside tetrahedra on the inner side, and it fills the cell
        {
          double wrong=0.0;
          for (size_t n=0;n<in.size();n++) {
            double c[3];
            for (int d=0;d<3;d++) c[d]=0.25*(in[n].x[0][d]+in[n].x[1][d]+in[n].x[2][d]+in[n].x[3][d]);
            if (GeneralInside(c)==false) wrong+=fabs(TetVolume(in[n].x));
          }
          wrong/=V,maxInside=std::max(maxInside,wrong);
          if (wrong>1.0E-6) badInside++;
          double errFill=fabs(vo+vi-V)/V;
          maxFill=std::max(maxFill,errFill);
          if (errFill>1.0E-9) badFill++;
        }

        //a cut cell written in full, or not at all
        if (cut&&(status==StatusComplete)) {
          maxBrick=std::max(maxBrick,vi/V);
          if (badBrick<5) printf("cutcell-tetra: cell [%e %e %e]-[%e %e %e]: cut (%.3e V_cell inside the body) but written as a full cell\n",lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],vi/V);
          badBrick++;
        }
        if (cut&&(status==StatusRemoved)) {
          maxRemoved=std::max(maxRemoved,vo/V);
          badRemoved++;
        }

        //outside: the output tetrahedra whose centroid is inside the body (parity; not reliable for slivers within the
        //round-off of the surface, so a cell fails when their volume exceeds 1e-6 V_cell)
        {
          double wrong=0.0;   //(a cut cell written in full is counted by the 'brick' check)
          if (status==StatusCut) for (size_t n=0;n<outT.size();n++) {
            double c[3];
            for (int d=0;d<3;d++) c[d]=0.25*(outT[n].x[0][d]+outT[n].x[1][d]+outT[n].x[2][d]+outT[n].x[3][d]);
            if (GeneralInside(c)) wrong+=fabs(TetVolume(outT[n].x));
          }
          wrong/=V,maxOutside=std::max(maxOutside,wrong);
          if (wrong>1.0E-6) {
            if (badOutside<5) printf("cutcell-tetra: cell [%e %e %e]-[%e %e %e]: output tetrahedra of %.3e V_cell inside the body\n",lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],wrong);
            badOutside++;
          }
        }

        //volume: against the reference decomposition and against the cell measure
        {
          double errVolume=fabs(vOut-vo)/V,errMeasure=fabs(vOut-std::max(cell->Measure,0.0))/V;
          maxVolume=std::max(maxVolume,errVolume);
          if (errVolume>1.0E-2) badVolume++;
          if (status==StatusCut) maxMeasure=std::max(maxMeasure,errMeasure);
          if ((status==StatusCut)&&(errMeasure>1.0E-2)) {
            if (badMeasure<5) printf("cutcell-tetra: cell [%.12e %.12e %.12e]-[%.12e %.12e %.12e] (level %i, process %i): output outside volume %.6e V_cell, cell measure %.6e V_cell\n",
                lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],node->RefinmentLevel,PIC::ThisThread,vOut/V,cell->Measure/V);
            badMeasure++;
          }
        }

        //the output tetrahedra lie in the cell and do not overlap
        if (status==StatusCut) {
          double errBox=0.0,errOverlap=0.0,h=0.0;
          for (int d=0;d<3;d++) h=std::max(h,hi[d]-lo[d]);
          for (size_t n=0;n<outT.size();n++) for (int v=0;v<4;v++) for (int d=0;d<3;d++) errBox=std::max(errBox,std::max(lo[d]-outT[n].x[v][d],outT[n].x[v][d]-hi[d])/h);
          if (errBox>1.0E-9) badBox++;
          maxBox=std::max(maxBox,errBox);

          for (size_t n=0;n<outT.size();n++) {
            double lo1[3],hi1[3];
            TetBox(outT[n].x,lo1,hi1);
            for (size_t m=n+1;m<outT.size();m++) {
              double lo2[3],hi2[3];
              TetBox(outT[m].x,lo2,hi2);
              bool separated=false;
              for (int d=0;d<3;d++) if ((hi1[d]<=lo2[d]+1.0E-13)||(hi2[d]<=lo1[d]+1.0E-13)) separated=true;
              if (!separated) separated=TetSeparated(outT[n].x,outT[m].x,1.0E-12);
              if (!separated) errOverlap=std::max(errOverlap,TetIntersection(outT[n].x,outT[m].x));
            }
          }
          errOverlap/=V,maxOverlap=std::max(maxOverlap,errOverlap);
          if (errOverlap>1.0E-9) badOverlap++;
        }

        //faces: the part of each cell face covered by the output equals that of the reference decomposition
        {
          double errFaces=0.0;
          for (int axis=0;axis<3;axis++) for (int sd=0;sd<2;sd++) {
            int a1=(axis+1)%3,a2=(axis+2)%3;
            double c=sd ? hi[axis] : lo[axis],faceArea=(hi[a1]-lo[a1])*(hi[a2]-lo[a2]);
            double area[2]={0.0,0.0};

            for (int which=0;which<2;which++) {
              const std::vector<CutCell::cCutCellTetrahedron>& t=(which==0) ? out : outT;
              if ((which==1)&&(status!=StatusCut)) {area[1]=(status==StatusComplete) ? faceArea : 0.0; continue;}
              for (size_t n=0;n<t.size();n++) {
                static const int tf[4][3]={{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
                for (int f=0;f<4;f++) {
                  bool onFace=true;
                  for (int v=0;v<3;v++) if (fabs(t[n].x[tf[f][v]][axis]-c)>1.0E-12) onFace=false;
                  if (!onFace) continue;
                  cPolygon3D tri(3,std::vector<double>(3));
                  for (int v=0;v<3;v++) for (int d=0;d<3;d++) tri[v][d]=t[n].x[tf[f][v]][d];
                  area[which]+=Area(tri);
                }
              }
            }
            errFaces=std::max(errFaces,fabs(area[1]-area[0])/faceArea);
          }
          maxFaces=std::max(maxFaces,errFaces);
          if (errFaces>1.0E-6) badFaces++;
          if (errFaces>1.0E-2) badFaces2++;
        }
        continue;
      }

      //the exact geometry of the cell (the body restricted to the cell: only its faces that cut the cell; a cell outside
      //all analytic spheres meets none of their components)
      bool farOutside=Analytic;
      for (size_t ns=0;(ns<Body.spheres.size())&&farOutside;ns++) {
        double d2=0.0;
        for (int d=0;d<3;d++) {
          double a=lo[d]-Body.spheres[ns].c[d],b=hi[d]-Body.spheres[ns].c[d];
          if (a>0.0) d2+=a*a;
          else if (b<0.0) d2+=b*b;
        }
        if (d2<=pow(Body.spheres[ns].r*(1.0+1.0E-9),2)) farOutside=false;
      }
      const cBody body=farOutside ? cBody() : RestrictToBox(Body,lo,hi);
      double vBody=IntersectionVolume(Box(lo,hi),body),vExact=V-vBody;

      //the output: GetCutcellTetrahedronMesh() (the outside tetrahedra, or the full cell, or nothing)
      list<cMesh::cTetrahedron> output;
      int status=StatusComplete;
      std::vector<std::vector<double> > outT;   //12 doubles per tetrahedron

      status=PIC::Mesh::mesh->GetCutcellTetrahedronMesh(output,i,j,k,node);

      for (list<cMesh::cTetrahedron>::iterator it=output.begin();it!=output.end();it++) {
        std::vector<double> t(12);
        for (int v=0;v<4;v++) for (int d=0;d<3;d++) t[3*v+d]=it->x[v][d];
        outT.push_back(t);
      }

      //the exact decomposition: the outside and inside tetrahedra
      std::vector<CutCell::cCutCellTetrahedron> out,in;
      PIC::Mesh::mesh->GetCutCellDecomposition(i,j,k,node,&out,&in);

      //analytic spheres: the outside volume of the output against the analytic one (the approximation of the spheres:
      //at most SphereApproximationTolerance*V_cell for each sphere that crosses the cell)
      auto CheckAnalytic=[&] (double vOut) {
        if (Analytic==false) return;
        double err=fabs(vOut-(V-AnalyticInsideVolume(lo,hi,eps)))/V;
        int nCrossing=0;
        for (size_t ns=0;ns<Spheres.size();ns++) {
          int st;
          double v=Spheres[ns]->GetRemainedBlockVolume(lo,hi,eps,&st);
          if ((v>0.0)&&(v<V)) nCrossing++;
        }
        err/=std::max(1,nCrossing);
        maxAnalytic=std::max(maxAnalytic,err);
        if (err>CutCell::SphereApproximationTolerance) {
          if (badAnalytic<5) printf("cutcell-tetra: cell [%.6e %.6e %.6e]-[%.6e %.6e %.6e] (level %i): outside volume %.6e V_cell, analytic %.6e V_cell\n",
              lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],node->RefinmentLevel,vOut/V,1.0-AnalyticInsideVolume(lo,hi,eps)/V);
          badAnalytic++;
        }
      };

      //a cell that is not cut is written in full (outside the body) or not at all (inside)
      bool isCut=(status==StatusCut)||((vBody>1.0E-9*V)&&(vBody<(1.0-1.0E-9)*V));
      if (!isCut) {
        double vOut=(status==StatusComplete) ? V : 0.0,errVolume=fabs(vOut-vExact)/V;
        maxVolume=std::max(maxVolume,errVolume);
        CheckAnalytic(vOut);
        if (errVolume>1.0E-9) {
          if (badVolume<5) printf("cutcell-tetra: cell [%e %e %e]-[%e %e %e] (process %i) status %i: vol(cell ∩ body)/V_cell=%e\n",lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],PIC::ThisThread,status,vBody/V);
          badVolume++;
        }
        continue;
      }
      nCut++;

      //the crossings of the edges of the cell with the surface
      int nCross[3][2][2];
      bool corner,along,parts,edge2=false,noEdge=true,poke=false;

      EdgeCrossings(lo,hi,nCross,corner,along,parts);
      for (int dir=0;dir<3;dir++) for (int s1=0;s1<2;s1++) for (int s2=0;s2<2;s2++) {
        if (nCross[dir][s1][s2]>=2) edge2=true;
        if (nCross[dir][s1][s2]>0) noEdge=false;
      }

      //outside: every output tetrahedron lies outside the body; a full cell is outside, a removed one inside
      double vOut=0.0,errOutside=0.0;
      if (status==StatusComplete) errOutside=vBody,vOut=V;
      else if (status==StatusRemoved) vOut=0.0;
      else for (size_t n=0;n<outT.size();n++) {
        double t[4][3];
        for (int v=0;v<4;v++) for (int d=0;d<3;d++) t[v][d]=outT[n][3*v+d];
        vOut+=fabs(TetVolume(t));
        errOutside=std::max(errOutside,IntersectionVolume(Tetrahedron(t),body));
      }
      errOutside/=V;
      maxOutside=std::max(maxOutside,errOutside);
      if (errOutside>1.0E-9) badOutside++;

      //volume: the outside volume of the output, and the cell measure
      double errVolume=fabs(vOut-vExact)/V,errMeasure=MeasureAvailable ? fabs(vOut-std::max(cell->Measure,0.0))/V : 0.0;
      maxVolume=std::max(maxVolume,errVolume),maxMeasure=std::max(maxMeasure,errMeasure);
      CheckAnalytic(vOut);
      if (errVolume>1.0E-9) badVolume++;
      if (errMeasure>1.0E-2) {
        if (badMeasure<5) printf("cutcell-tetra: cell [%.12e %.12e %.12e]-[%.12e %.12e %.12e] (level %i, process %i): outside volume %.6e V_cell (exact %.6e), cell measure %.6e V_cell\n",
            lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],node->RefinmentLevel,PIC::ThisThread,vOut/V,vExact/V,cell->Measure/V);
        badMeasure++;
      }

      //faces: the outside part of each face of the cell is covered by the faces of the outside tetrahedra
      double errFaces=0.0;
      for (int axis=0;axis<3;axis++) for (int s=0;s<2;s++) {
        int a1=(axis+1)%3,a2=(axis+2)%3;
        double c=s ? hi[axis] : lo[axis],area=0.0,exactArea;
        cPolygon3D rect(4,std::vector<double>(3));

        for (int m=0;m<4;m++) {
          rect[m][axis]=c;
          rect[m][a1]=((m==1)||(m==2)) ? hi[a1] : lo[a1];
          rect[m][a2]=(m>=2) ? hi[a2] : lo[a2];
        }
        double faceArea=(hi[a1]-lo[a1])*(hi[a2]-lo[a2]);
        exactArea=faceArea-IntersectionArea(rect,body);

        if (status==StatusComplete) area=faceArea;
        else if (status==StatusCut) for (size_t n=0;n<outT.size();n++) {
          static const int tf[4][3]={{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
          for (int f=0;f<4;f++) {
            bool onFace=true;
            for (int v=0;v<3;v++) if (fabs(outT[n][3*tf[f][v]+axis]-c)>1.0E-12) onFace=false;
            if (!onFace) continue;
            cPolygon3D tri(3,std::vector<double>(3));
            for (int v=0;v<3;v++) for (int d=0;d<3;d++) tri[v][d]=outT[n][3*tf[f][v]+d];
            area+=Area(tri);
          }
        }
        errFaces=std::max(errFaces,fabs(area-exactArea)/faceArea);

        //the surface crosses the face, but none of its edges
        if ((exactArea>1.0E-9*faceArea)&&(exactArea<(1.0-1.0E-9)*faceArea)&&
            (nCross[a1][0][s]+nCross[a1][1][s]+nCross[a2][s][0]+nCross[a2][s][1]==0)) poke=true;
      }

      if (edge2) nEdge2++;
      if (noEdge) nNoEdge++;
      if (corner) nCorner++;
      if (along) nAlong++;
      if (parts) nParts++;
      if (poke) nPoke++;
      maxFaces=std::max(maxFaces,errFaces);
      if (errFaces>1.0E-9) badFaces++;

      //inside, overlap (the exact decomposition only)
      if (!legacy) {
        double errInside=0.0,errOverlap=0.0,vAll=0.0;
        std::vector<CutCell::cCutCellTetrahedron> all(out);
        all.insert(all.end(),in.begin(),in.end());

        for (size_t n=0;n<in.size();n++) errInside=std::max(errInside,fabs(TetVolume(in[n].x))-IntersectionVolume(Tetrahedron(in[n].x),body));
        for (size_t n=0;n<all.size();n++) vAll+=fabs(TetVolume(all[n].x));

        for (size_t n=0;n<all.size();n++) {
          double lo1[3],hi1[3];
          TetBox(all[n].x,lo1,hi1);
          for (size_t m=n+1;m<all.size();m++) {
            double lo2[3],hi2[3];
            TetBox(all[m].x,lo2,hi2);
            bool separated=false;
            for (int d=0;d<3;d++) if ((hi1[d]<=lo2[d]+1.0E-13)||(hi2[d]<=lo1[d]+1.0E-13)) separated=true;
            if (!separated) separated=TetSeparated(all[n].x,all[m].x,1.0E-12);
            if (!separated) errOverlap=std::max(errOverlap,TetIntersection(all[n].x,all[m].x));
          }
        }

        errInside/=V,errOverlap/=V;
        double errFill=fabs(vAll-V)/V;
        maxInside=std::max(maxInside,errInside),maxOverlap=std::max(maxOverlap,errOverlap),maxFill=std::max(maxFill,errFill);
        if (errInside>1.0E-9) badInside++;
        if (errOverlap>1.0E-9) badOverlap++;
        if (errFill>1.0E-9) badFill++;
        nTetIn+=in.size();
      }

      nTetOut+=outT.size();
      for (size_t n=0;n<outT.size();n++) outputTets.insert(outputTets.end(),outT[n].begin(),outT[n].end());

      cCellRecord r;
      r.key[0]=node->RefinmentLevel;
      for (int d=0;d<3;d++) r.key[d+1]=lround((lo[d]+1.0)/dx[d]);
      r.nOut=outT.size(),r.nIn=in.size(),r.status=status,r.vOut=vOut,r.vIn=0.0,r.vExact=vExact,r.measure=cell->Measure;
      for (size_t n=0;n<in.size();n++) r.vIn+=fabs(TetVolume(in[n].x));
      records.push_back(r);
    }
  }

  tCells=MPI_Wtime();
  long lv[12]={nCells,nCut,nTetOut,nTetIn,badOutside,badInside,badOverlap,badFill,badVolume,badMeasure,badFaces,0},glv[12];
  double dv[7]={maxOutside,maxInside,maxOverlap,maxFill,maxVolume,maxMeasure,maxFaces},gdv[7];
  MPI_Allreduce(lv,glv,12,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(dv,gdv,7,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) printf("cutcell-tetra %s%s np=%i: %li cells, %li cut, %li outside and %li inside tetrahedra, surface %i faces, body volume %.12e\n",
      shape.c_str(),legacy ? " (legacy)" : "",PIC::nTotalThreads,glv[0],glv[1],glv[2],glv[3],CutCell::nBoundaryTriangleFaces,Body.volume);

  //the faces of the box lie on cell faces: no cell is cut
  if (!General) Check(((glv[1]>0)||(shape=="box"))&&(glv[4]==0),"outside","%li cut cells, %li with output tetrahedra inside the body (max vol(T ∩ body)/V_cell %.3e)",glv[1],glv[4],gdv[0]);
  if ((!legacy)&&(!General)) {
    Check(glv[5]==0,"inside","%li cut cells with inside tetrahedra not inside the body (max %.3e V_cell)",glv[5],gdv[1]);
    Check((glv[6]==0)&&(glv[7]==0),"overlap","%li cut cells with overlapping tetrahedra (max %.3e V_cell), %li not filled (max |sum-V_cell| %.3e V_cell)",glv[6],gdv[2],glv[7],gdv[3]);
  }
  if (General) {
    long lg[8]={badClass,badBrick,badRemoved,badBox,badClassMeasure,badFaces2,nRefCut,nOutCut},glg[8];
    double dg[3]={maxBrick,maxRemoved,maxBox},gdg[3];
    MPI_Allreduce(lg,glg,8,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(dg,gdg,3,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    const char *procedure=legacy ? "the old procedure (AMPS_CUTCELL_TETRA_LEGACY=1)" : "the new procedure";

    if (PIC::ThisThread==0) fprintf(Report,"INFO output of GetCutcellTetrahedronMesh(): %s; the reference geometry: GetCutCellTetrahedra() (exact decomposition)\n"
        "INFO %li cells, %li cut by the reference decomposition (more than 1e-6 V_cell on each side), %li written as cut cells (outside tetrahedra)\n",procedure,glv[0],glg[6],glg[7]);
    //the reference decomposition
    Check(glv[7]==0,"fill","reference: %li cut cells whose tetrahedra do not fill the cell (max %.3e V_cell)",glv[7],gdv[3]);
    Check(glv[5]==0,"reference-inside","reference: %li cut cells with more than 1e-6 V_cell of inside tetrahedra on the outer side of the surface (max %.3e)",glv[5],gdv[1]);
    //the output
    Check(glv[4]==0,"outside","output: %li of the %li cells written as cut cells with more than 1e-6 V_cell of output (outside) tetrahedra whose centroid is inside the body (max %.3e V_cell)",glv[4],glg[7],gdv[0]);
    Check(glv[9]==0,"measure","output: %li cut cells whose output outside tetrahedra differ from the InitCellMeasure() measure by more than 1e-2 V_cell (max %.3e)",glv[9],gdv[5]);
    Check(glv[8]==0,"volume","output: %li cut cells whose output outside volume (tetrahedra, full cell or nothing) differs from the exact decomposition by more than 1e-2 V_cell (max %.3e)",glv[8],gdv[4]);
    Check(glg[1]==0,"brick","output: %li cut cells written as a full (uncut) cell (max volume inside the body %.3e V_cell)",glg[1],gdg[0]);
    Check(glg[2]==0,"removed","output: %li cut cells not written at all (max volume outside the body %.3e V_cell)",glg[2],gdg[1]);
    Check(glg[0]==0,"uncut","output: %li cells that are not cut written on the wrong side (in full inside the body, or not at all outside); the InitCellMeasure() measure agrees with the reference for %li of them",glg[0],glg[4]);
    Check(glv[6]==0,"overlap","output: %li cut cells with overlapping output tetrahedra (max %.3e V_cell)",glv[6],gdv[2]);
    Check(glg[3]==0,"box","output: %li cut cells with output tetrahedra outside the cell (max %.3e of the cell size)",glg[3],gdg[2]);
    Check(glv[10]==0,"faces","output: %li cut cells whose outside part of a cell face differs from the exact decomposition by more than 1e-6 of the face area, %li by more than 1e-2 (max %.3e)",glv[10],glg[5],gdv[6]);
    if (PIC::ThisThread==0) {fprintf(Report,"%s %s: %i checks failed\n",(nFailed==0) ? "PASS" : "FAIL",shape.c_str(),nFailed); fclose(Report);}
    MPI_Finalize();
    return (nFailed==0) ? 0 : 1;
  }
  {
    long lc[7]={nEdge2,nNoEdge,nCorner,nAlong,nParts,nPoke,badAnalytic},glc[7];
    double gma;
    MPI_Allreduce(lc,glc,7,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&maxAnalytic,&gma,1,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) fprintf(Report,"INFO crossings of the %li cut cells: %li with an edge crossed 2+ times, %li with no edge crossed, %li where the surface crosses a face but none of its edges, "
        "%li with a corner on the surface, %li with an edge in the surface, %li with an edge crossed by two spheres\n",glv[1],glc[0],glc[1],glc[5],glc[2],glc[3],glc[4]);

    //the cases built for these crossings must produce them
    if (shape=="plate") Check(glc[0]>0,"crossings","%li cut cells with an edge that enters and leaves the plate",glc[0]);
    if (shape=="bumps") Check((glc[0]>0)&&(glc[5]>0)&&(glc[1]>0),"crossings","%li cut cells with an edge crossed twice, %li with a face crossed but none of its edges, %li with no edge crossed",glc[0],glc[5],glc[1]);
    if (shape=="lcorner") Check(glc[0]>0,"crossings","%li cut cells with an edge crossed by both faces of the concave corner",glc[0]);
    if (shape=="wedge") Check((glc[2]>0)&&(glc[3]>0),"crossings","%li cut cells with a corner on the surface, %li with an edge in the surface",glc[2],glc[3]);
    if (shape=="aspheres") Check((glc[4]>0)&&(glc[0]>0),"crossings","%li cut cells with an edge crossed by both spheres, %li with an edge crossed 2+ times",glc[4],glc[0]);

    if (Analytic) Check(glc[6]==0,"analytic","%li cells whose outside volume differs from the analytic one (cInternalSphericalData::GetRemainedBlockVolume()) by more than %.1e V_cell per sphere that crosses the cell "
        "(max %.3e): the approximation of the spheres",glc[6],CutCell::SphereApproximationTolerance,gma);
  }

  Check(glv[8]==0,"volume","%li cells whose outside tetrahedra (or full cell) differ from the exact outside volume by more than 1e-9 V_cell (max %.3e)",glv[8],gdv[4]);
  if (MeasureAvailable) Check(glv[9]==0,"measure","%li cut cells whose outside tetrahedra differ from the InitCellMeasure() measure by more than 1e-2 V_cell (max %.3e)",glv[9],gdv[5]);
  else if (PIC::ThisThread==0) fprintf(Report,"INFO measure: not checked (InitCellMeasure() permits one internal boundary per cell; the analytic check replaces it)\n");
  Check(glv[10]==0,"faces","%li cut cells whose outside tetrahedra do not cover the outside part of a cell face exactly (max area error %.3e of the face)",glv[10],gdv[6]);

  //the output file
  {
    sprintf(fname,"cutcell-tetra.%s%s.np%i.dat",shape.c_str(),legacy ? ".legacy" : "",PIC::nTotalThreads);
    PIC::Mesh::mesh->SetAssembleDistributedOutputFileFlag(true);
    PIC::Mesh::mesh->OutputDistributedDataTECPLOT(fname,true,0);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    //all outside tetrahedra of the checks on process 0
    int nLocal=outputTets.size(),*counts=NULL,*displs=NULL,nTotal=0;
    std::vector<double> all;
    if (PIC::ThisThread==0) counts=new int[PIC::nTotalThreads],displs=new int[PIC::nTotalThreads];
    MPI_Gather(&nLocal,1,MPI_INT,counts,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread==0) {
      for (int t=0;t<PIC::nTotalThreads;t++) displs[t]=nTotal,nTotal+=counts[t];
      all.resize(std::max(nTotal,1));
    }
    MPI_Gatherv(outputTets.data(),nLocal,MPI_DOUBLE,all.data(),counts,displs,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      std::vector<cTetKey> checked,file;
      int nVariables;
      long nBadColumns,nBadConnectivity,nMismatch=0;
      double maxDiff=0.0;
      bool found=ReadTetraZones(fname,file,nVariables,nBadColumns,nBadConnectivity);
      if (getenv("AMPS_CUTCELL_TETRA_KEEP_OUTPUT")==NULL) remove(fname);   //(large; AMPS_CUTCELL_TETRA_KEEP_OUTPUT=1 keeps it)

      for (int n=0;n<nTotal/12;n++) {
        double t[4][3];
        for (int v=0;v<4;v++) for (int d=0;d<3;d++) t[v][d]=all[12*n+3*v+d];
        //as printed (%e)
        for (int v=0;v<4;v++) for (int d=0;d<3;d++) {char s[64]; sprintf(s,"%e",t[v][d]); t[v][d]=atof(s);}
        checked.push_back(MakeTetKey(t));
      }

      std::sort(checked.begin(),checked.end()),std::sort(file.begin(),file.end());
      if (checked.size()!=file.size()) nMismatch=std::max(checked.size(),file.size());
      else for (size_t n=0;n<checked.size();n++) for (int d=0;d<12;d++) {
        double diff=fabs(checked[n].x[d]-file[n].x[d]);
        maxDiff=std::max(maxDiff,diff);
        if (diff>1.0E-6) {nMismatch++; break;}
      }

      Check(found&&(nMismatch==0)&&(nBadColumns==0)&&(nBadConnectivity==0)&&((file.size()>0)||(shape=="box")),"file",
          "%s: %li tetrahedra in the file, %li checked; %li differ (max coordinate difference %.3e); %li points without the %i values of the variable list; %li invalid elements",
          fname,(long)file.size(),(long)checked.size(),nMismatch,maxDiff,nBadColumns,nVariables,nBadConnectivity);
      delete [] counts; delete [] displs;
    }
    else Check(true,"file","");
  }

  tOutput=MPI_Wtime();
  if (PIC::ThisThread==0) printf("cutcell-tetra %s: time: mesh %.1f s, normals and cell measures %.1f s, checks of the cells %.1f s, output file %.1f s\n",
      shape.c_str(),tMesh-tStart,tMeasure-tMesh,tCells-tMeasure,tOutput-tCells);

  //the cut cells in the canonical order
  {
    std::vector<long> keys;
    std::vector<double> values;
    std::vector<int> ints;
    for (size_t n=0;n<records.size();n++) {
      keys.insert(keys.end(),records[n].key,records[n].key+4);
      values.push_back(records[n].vOut),values.push_back(records[n].vIn),values.push_back(records[n].measure);
      ints.push_back(records[n].nOut),ints.push_back(records[n].nIn),ints.push_back(records[n].status);
    }

    int nLocal=records.size(),*counts=NULL,*displs=NULL,*c3=NULL,*d3=NULL,*c4=NULL,*d4=NULL,nTotal=0;
    if (PIC::ThisThread==0) {
      counts=new int[PIC::nTotalThreads],displs=new int[PIC::nTotalThreads],c3=new int[PIC::nTotalThreads],d3=new int[PIC::nTotalThreads];
      c4=new int[PIC::nTotalThreads],d4=new int[PIC::nTotalThreads];
    }
    MPI_Gather(&nLocal,1,MPI_INT,counts,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    std::vector<long> allKeys;
    std::vector<double> allValues;
    std::vector<int> allInts;
    if (PIC::ThisThread==0) {
      for (int t=0;t<PIC::nTotalThreads;t++) displs[t]=nTotal,c3[t]=3*counts[t],d3[t]=3*nTotal,c4[t]=4*counts[t],d4[t]=4*nTotal,nTotal+=counts[t];
      allKeys.resize(4*std::max(nTotal,1)),allValues.resize(3*std::max(nTotal,1)),allInts.resize(3*std::max(nTotal,1));
    }
    MPI_Gatherv(keys.data(),4*nLocal,MPI_LONG,allKeys.data(),c4,d4,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gatherv(values.data(),3*nLocal,MPI_DOUBLE,allValues.data(),c3,d3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gatherv(ints.data(),3*nLocal,MPI_INT,allInts.data(),c3,d3,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      std::vector<cCellRecord> all(nTotal);
      for (int n=0;n<nTotal;n++) {
        for (int d=0;d<4;d++) all[n].key[d]=allKeys[4*n+d];
        all[n].vOut=allValues[3*n],all[n].vIn=allValues[3*n+1],all[n].measure=allValues[3*n+2];
        all[n].nOut=allInts[3*n],all[n].nIn=allInts[3*n+1],all[n].status=allInts[3*n+2];
      }
      std::sort(all.begin(),all.end());

      sprintf(fname,"cutcell-tetra.%s%s.np%i.cells",shape.c_str(),legacy ? ".legacy" : "",PIC::nTotalThreads);
      FILE *f=fopen(fname,"w");
      for (int n=0;n<nTotal;n++) fprintf(f,"%li %li %li %li %i %i %i %.17e %.17e %.17e\n",all[n].key[0],all[n].key[1],all[n].key[2],all[n].key[3],
          all[n].status,all[n].nOut,all[n].nIn,all[n].vOut,all[n].vIn,all[n].measure);
      fclose(f);
      delete [] counts; delete [] displs; delete [] c3; delete [] d3; delete [] c4; delete [] d4;
    }
  }

  if (PIC::ThisThread==0) {
    fprintf(Report,"%s %s: %i checks failed\n",(nFailed==0) ? "PASS" : "FAIL",shape.c_str(),nFailed);
    fclose(Report);
  }

  MPI_Finalize();
  return (nFailed==0) ? 0 : 1;
}

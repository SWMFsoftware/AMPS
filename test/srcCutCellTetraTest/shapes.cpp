//Triangulated test bodies with exactly known geometry and the exact convex geometry of the checks (see shapes.h)
#include "shapes.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>

namespace CutCellTetraTest {

namespace {
const double kPi=3.14159265358979323846;
const double kCenter[3]={0.0137,-0.0211,0.0173};

double Dot(const double *a,const double *b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}

//a rotation that aligns no axis of the body with the mesh
void Rotation(double R[3][3]) {
  const double a=0.43,b=0.71,c=0.29;
  double Rx[3][3]={{1,0,0},{0,cos(a),-sin(a)},{0,sin(a),cos(a)}};
  double Ry[3][3]={{cos(b),0,sin(b)},{0,1,0},{-sin(b),0,cos(b)}};
  double Rz[3][3]={{cos(c),-sin(c),0},{sin(c),cos(c),0},{0,0,1}};
  double T[3][3];
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {T[i][j]=0.0; for (int k=0;k<3;k++) T[i][j]+=Ry[i][k]*Rx[k][j];}
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {R[i][j]=0.0; for (int k=0;k<3;k++) R[i][j]+=Rz[i][k]*T[k][j];}
}

//map local coordinates: x = center + R*u (rotate==false: x = center + u)
void Map(const double *u,double *x,bool rotate,const double *center) {
  double R[3][3];
  Rotation(R);
  for (int i=0;i<3;i++) x[i]=center[i]+(rotate ? R[i][0]*u[0]+R[i][1]*u[1]+R[i][2]*u[2] : u[i]);
}

//the half-space n_local.u<=d_local in global coordinates
cPlane MapPlane(const double *nl,double dl,bool rotate,const double *center) {
  double R[3][3],l=sqrt(Dot(nl,nl));
  cPlane p;
  Rotation(R);
  for (int i=0;i<3;i++) p.n[i]=(rotate ? R[i][0]*nl[0]+R[i][1]*nl[1]+R[i][2]*nl[2] : nl[i])/l;
  p.d=dl/l+Dot(p.n,center);
  return p;
}

int AddVertex(cBody& b,const double *x) {
  for (int d=0;d<3;d++) b.vertex.push_back(x[d]);
  return b.vertex.size()/3-1;
}

bool InsideBody(const cBody& b,const double *x) {
  for (size_t c=0;c<b.components.size();c++) {
    bool in=true;
    for (size_t n=0;(n<b.components[c].planes.size())&&in;n++) if (Dot(b.components[c].planes[n].n,x)>b.components[c].planes[n].d) in=false;
    if (in) return true;
  }
  return false;
}

//orient every face outward (the point just in front of the face center is outside) and get the volume
void Finish(cBody& b) {
  b.volume=0.0;

  for (size_t n=0;n<b.face.size();n+=3) {
    double *a=&b.vertex[3*b.face[n]],*p=&b.vertex[3*b.face[n+1]],*q=&b.vertex[3*b.face[n+2]];
    double e1[3],e2[3],nn[3],c[3],x[3],l;

    for (int d=0;d<3;d++) e1[d]=p[d]-a[d],e2[d]=q[d]-a[d],c[d]=(a[d]+p[d]+q[d])/3.0;
    nn[0]=e1[1]*e2[2]-e1[2]*e2[1],nn[1]=e1[2]*e2[0]-e1[0]*e2[2],nn[2]=e1[0]*e2[1]-e1[1]*e2[0];
    l=sqrt(Dot(nn,nn));
    for (int d=0;d<3;d++) x[d]=c[d]+1.0E-7*nn[d]/l;
    if (InsideBody(b,x)) std::swap(b.face[n+1],b.face[n+2]);

    a=&b.vertex[3*b.face[n]],p=&b.vertex[3*b.face[n+1]],q=&b.vertex[3*b.face[n+2]];
    b.volume+=(a[0]*(p[1]*q[2]-p[2]*q[1])-a[1]*(p[0]*q[2]-p[2]*q[0])+a[2]*(p[0]*q[1]-p[1]*q[0]))/6.0;
  }
}

void Icosphere(cBody& b,int nSubdivisions,double radius) {
  const double t=(1.0+sqrt(5.0))/2.0;
  double v0[12][3]={{-1,t,0},{1,t,0},{-1,-t,0},{1,-t,0},{0,-1,t},{0,1,t},{0,-1,-t},{0,1,-t},{t,0,-1},{t,0,1},{-t,0,-1},{-t,0,1}};
  int f0[20][3]={{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
                 {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};
  std::vector<double> v;
  std::vector<int> f;

  for (int i=0;i<12;i++) {double l=sqrt(Dot(v0[i],v0[i])); for (int d=0;d<3;d++) v.push_back(v0[i][d]/l);}
  for (int i=0;i<20;i++) for (int d=0;d<3;d++) f.push_back(f0[i][d]);

  for (int s=0;s<nSubdivisions;s++) {
    std::map<std::pair<int,int>,int> middle;
    std::vector<int> fnew;
    auto Middle=[&](int a,int c) {
      std::pair<int,int> key(std::min(a,c),std::max(a,c));
      std::map<std::pair<int,int>,int>::iterator it=middle.find(key);
      if (it!=middle.end()) return it->second;
      double m[3],l;
      for (int d=0;d<3;d++) m[d]=0.5*(v[3*a+d]+v[3*c+d]);
      l=sqrt(Dot(m,m));
      for (int d=0;d<3;d++) v.push_back(m[d]/l);
      return middle[key]=v.size()/3-1;
    };

    for (size_t n=0;n<f.size();n+=3) {
      int a=f[n],c=f[n+1],e=f[n+2],ac=Middle(a,c),ce=Middle(c,e),ea=Middle(e,a);
      int add[12]={a,ac,ea, c,ce,ac, e,ea,ce, ac,ce,ea};
      fnew.insert(fnew.end(),add,add+12);
    }
    f.swap(fnew);
  }

  for (size_t n=0;n<v.size();n+=3) {double u[3]={radius*v[n],radius*v[n+1],radius*v[n+2]},x[3]; Map(u,x,false,kCenter); AddVertex(b,x);}
  b.face=f;

  //the component: the planes of the faces (outward: the center is inside)
  cComponent c;
  for (size_t n=0;n<f.size();n+=3) {
    const double *a=&b.vertex[3*f[n]],*p=&b.vertex[3*f[n+1]],*q=&b.vertex[3*f[n+2]];
    double e1[3],e2[3],nn[3];
    for (int d=0;d<3;d++) e1[d]=p[d]-a[d],e2[d]=q[d]-a[d];
    nn[0]=e1[1]*e2[2]-e1[2]*e2[1],nn[1]=e1[2]*e2[0]-e1[0]*e2[2],nn[2]=e1[0]*e2[1]-e1[1]*e2[0];
    double l=sqrt(Dot(nn,nn)),r[3]={a[0]-kCenter[0],a[1]-kCenter[1],a[2]-kCenter[2]};
    if (Dot(nn,r)<0.0) l=-l;
    cPlane pl;
    for (int d=0;d<3;d++) pl.n[d]=nn[d]/l;
    pl.d=Dot(pl.n,a);
    c.planes.push_back(pl);
  }
  b.components.push_back(c);
}

//a cone (apex up) or a cylinder with nSeg segments, rotated; the caps are fans around a center vertex (coplanar faces)
void ConeCylinder(cBody& b,bool cone,int nSeg) {
  const double R=cone ? 0.45 : 0.38,zlo=cone ? -0.35 : -0.42,zhi=cone ? 0.45 : 0.42;
  std::vector<int> lo,hi;
  double u[3],x[3];

  for (int i=0;i<nSeg;i++) {
    double phi=2.0*kPi*i/nSeg;
    u[0]=R*cos(phi),u[1]=R*sin(phi),u[2]=zlo; Map(u,x,true,kCenter); lo.push_back(AddVertex(b,x));
    if (!cone) {u[2]=zhi; Map(u,x,true,kCenter); hi.push_back(AddVertex(b,x));}
  }

  u[0]=0.0,u[1]=0.0,u[2]=zlo; Map(u,x,true,kCenter); int clo=AddVertex(b,x);
  u[2]=zhi; Map(u,x,true,kCenter); int chi=AddVertex(b,x);

  for (int i=0;i<nSeg;i++) {
    int j=(i+1)%nSeg;
    int add[3]={clo,lo[j],lo[i]};
    b.face.insert(b.face.end(),add,add+3);

    if (cone) {int s[3]={lo[i],lo[j],chi}; b.face.insert(b.face.end(),s,s+3);}
    else {
      int t[3]={chi,hi[i],hi[j]},s1[3]={lo[i],lo[j],hi[j]},s2[3]={lo[i],hi[j],hi[i]};
      b.face.insert(b.face.end(),t,t+3),b.face.insert(b.face.end(),s1,s1+3),b.face.insert(b.face.end(),s2,s2+3);
    }
  }

  cComponent c;
  double n[3]={0,0,-1};
  c.planes.push_back(MapPlane(n,-zlo,true,kCenter));
  if (!cone) {double m[3]={0,0,1}; c.planes.push_back(MapPlane(m,zhi,true,kCenter));}

  for (int i=0;i<nSeg;i++) {
    double p0=2.0*kPi*i/nSeg,p1=2.0*kPi*(i+1)/nSeg,a[3]={R*cos(p0),R*sin(p0),zlo},q[3]={R*cos(p1),R*sin(p1),zlo};
    double nl[3],dl;

    if (cone) {
      double t[3]={0,0,zhi},e1[3],e2[3];
      for (int d=0;d<3;d++) e1[d]=q[d]-a[d],e2[d]=t[d]-a[d];
      nl[0]=e1[1]*e2[2]-e1[2]*e2[1],nl[1]=e1[2]*e2[0]-e1[0]*e2[2],nl[2]=e1[0]*e2[1]-e1[1]*e2[0];
    }
    else nl[0]=cos(0.5*(p0+p1)),nl[1]=sin(0.5*(p0+p1)),nl[2]=0.0;

    if (Dot(nl,a)<0.0) for (int d=0;d<3;d++) nl[d]=-nl[d];
    dl=Dot(nl,a);
    c.planes.push_back(MapPlane(nl,dl,true,kCenter));
  }

  b.components.push_back(c);
}

//an axis-aligned box whose faces lie on mesh planes; each face is split into 2x2 squares (8 triangles), so faces
//and edges of the surface pass through corners and edges of cells
void AlignedBox(cBody& b) {
  const double lo[3]={-0.4,-0.3,-0.2},hi[3]={0.4,0.3,0.2};
  const double zero[3]={0.0,0.0,0.0};
  std::map<std::vector<long>,int> index;

  auto Vertex=[&](double x,double y,double z) {
    std::vector<long> key(3);
    double p[3]={x,y,z};
    for (int d=0;d<3;d++) key[d]=lround(p[d]*1.0E6);
    std::map<std::vector<long>,int>::iterator it=index.find(key);
    if (it!=index.end()) return it->second;
    return index[key]=AddVertex(b,p);
  };

  for (int axis=0;axis<3;axis++) for (int s=0;s<2;s++) {
    int a1=(axis+1)%3,a2=(axis+2)%3;

    for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
      double c[4][3];
      for (int k=0;k<4;k++) {
        int ii=i+((k==1)||(k==2)),jj=j+(k>=2);
        c[k][axis]=s ? hi[axis] : lo[axis];
        c[k][a1]=lo[a1]+0.5*ii*(hi[a1]-lo[a1]);
        c[k][a2]=lo[a2]+0.5*jj*(hi[a2]-lo[a2]);
      }
      int v[4];
      for (int k=0;k<4;k++) v[k]=Vertex(c[k][0],c[k][1],c[k][2]);
      int t1[3]={v[0],v[1],v[2]},t2[3]={v[0],v[2],v[3]};
      b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);
    }
  }

  cComponent c;
  for (int axis=0;axis<3;axis++) {
    double n[3]={0,0,0};
    n[axis]=1.0; c.planes.push_back(MapPlane(n,hi[axis],false,zero));
    n[axis]=-1.0; c.planes.push_back(MapPlane(n,-lo[axis],false,zero));
  }
  b.components.push_back(c);
}

//the L-shaped prism: the union of the boxes A=[-a,a]x[-h,h]x[-a,0] and B=[-a,0]x[-h,h]x[0,a] (local u,v,w), rotated
void LShape(cBody& b) {
  const double a=0.45,h=0.3;
  //the L polygon in the (u,w) plane, counterclockwise, with the vertex (-a,0) on the left edge
  const double L[7][2]={{-a,-a},{a,-a},{a,0.0},{0.0,0.0},{0.0,a},{-a,a},{-a,0.0}};
  int vlo[7],vhi[7];

  for (int i=0;i<7;i++) {
    double u[3]={L[i][0],-h,L[i][1]},x[3];
    Map(u,x,true,kCenter); vlo[i]=AddVertex(b,x);
    u[1]=h; Map(u,x,true,kCenter); vhi[i]=AddVertex(b,x);
  }

  //the caps: the rectangle A (vertices 0,1,2,3,6) and B (6,3,4,5)
  const int capA[3][3]={{0,1,2},{0,2,3},{0,3,6}},capB[2][3]={{6,3,4},{6,4,5}};
  for (int k=0;k<3;k++) {
    int t1[3]={vlo[capA[k][0]],vlo[capA[k][1]],vlo[capA[k][2]]},t2[3]={vhi[capA[k][0]],vhi[capA[k][1]],vhi[capA[k][2]]};
    b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);
  }
  for (int k=0;k<2;k++) {
    int t1[3]={vlo[capB[k][0]],vlo[capB[k][1]],vlo[capB[k][2]]},t2[3]={vhi[capB[k][0]],vhi[capB[k][1]],vhi[capB[k][2]]};
    b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);
  }

  //the sides
  for (int i=0;i<7;i++) {
    int j=(i+1)%7;
    int t1[3]={vlo[i],vlo[j],vhi[j]},t2[3]={vlo[i],vhi[j],vhi[i]};
    b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);
  }

  double ex[3]={1,0,0},ey[3]={0,1,0},ez[3]={0,0,1},mx[3]={-1,0,0},my[3]={0,-1,0},mz[3]={0,0,-1};
  cComponent A,B;
  A.planes.push_back(MapPlane(ex,a,true,kCenter)),A.planes.push_back(MapPlane(mx,a,true,kCenter));
  A.planes.push_back(MapPlane(ey,h,true,kCenter)),A.planes.push_back(MapPlane(my,h,true,kCenter));
  A.planes.push_back(MapPlane(ez,0.0,true,kCenter)),A.planes.push_back(MapPlane(mz,a,true,kCenter));
  B.planes.push_back(MapPlane(ex,0.0,true,kCenter)),B.planes.push_back(MapPlane(mx,a,true,kCenter));
  B.planes.push_back(MapPlane(ey,h,true,kCenter)),B.planes.push_back(MapPlane(my,h,true,kCenter));
  B.planes.push_back(MapPlane(ez,a,true,kCenter)),B.planes.push_back(MapPlane(mz,0.0,true,kCenter));
  b.components.push_back(A),b.components.push_back(B);
}

//------------- bodies with multiple crossings of cell edges -------------
//a frame: x = c + R u (R = Rz(az) Ry(ay) Rx(ax))
struct cFrame {double R[3][3],c[3];};

cFrame MakeFrame(double ax,double ay,double az,double cx,double cy,double cz) {
  double Rx[3][3]={{1,0,0},{0,cos(ax),-sin(ax)},{0,sin(ax),cos(ax)}};
  double Ry[3][3]={{cos(ay),0,sin(ay)},{0,1,0},{-sin(ay),0,cos(ay)}};
  double Rz[3][3]={{cos(az),-sin(az),0},{sin(az),cos(az),0},{0,0,1}};
  double T[3][3];
  cFrame f;

  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {T[i][j]=0.0; for (int k=0;k<3;k++) T[i][j]+=Ry[i][k]*Rx[k][j];}
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {f.R[i][j]=0.0; for (int k=0;k<3;k++) f.R[i][j]+=Rz[i][k]*T[k][j];}
  if ((ax==0.0)&&(ay==0.0)&&(az==0.0)) for (int i=0;i<3;i++) for (int j=0;j<3;j++) f.R[i][j]=(i==j) ? 1.0 : 0.0;
  f.c[0]=cx,f.c[1]=cy,f.c[2]=cz;
  return f;
}

void FMap(const cFrame& f,const double *u,double *x) {
  for (int i=0;i<3;i++) x[i]=f.c[i]+f.R[i][0]*u[0]+f.R[i][1]*u[1]+f.R[i][2]*u[2];
}

//the half-space nl.u<=dl (local) in global coordinates
cPlane FPlane(const cFrame& f,const double *nl,double dl) {
  double l=sqrt(Dot(nl,nl));
  cPlane p;
  for (int i=0;i<3;i++) p.n[i]=(f.R[i][0]*nl[0]+f.R[i][1]*nl[1]+f.R[i][2]*nl[2])/l;
  p.d=dl/l+Dot(p.n,f.c);
  return p;
}

//a convex polygon (local u,v, counterclockwise) extruded along w over [w0,w1]: the surface and one component
void AddPrism(cBody& b,const cFrame& f,const std::vector<std::vector<double> >& poly,double w0,double w1) {
  int n=poly.size();
  std::vector<int> lo(n),hi(n);
  cComponent c;

  for (int i=0;i<n;i++) {
    double u[3]={poly[i][0],poly[i][1],w0},x[3];
    FMap(f,u,x); lo[i]=AddVertex(b,x);
    u[2]=w1; FMap(f,u,x); hi[i]=AddVertex(b,x);
  }

  for (int i=1;i+1<n;i++) {
    int t1[3]={lo[0],lo[i],lo[i+1]},t2[3]={hi[0],hi[i],hi[i+1]};
    b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);
  }

  for (int i=0;i<n;i++) {
    int j=(i+1)%n;
    int t1[3]={lo[i],lo[j],hi[j]},t2[3]={lo[i],hi[j],hi[i]};
    b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);

    //the side plane: outward normal of the edge (counterclockwise polygon)
    double nl[3]={poly[j][1]-poly[i][1],-(poly[j][0]-poly[i][0]),0.0};
    c.planes.push_back(FPlane(f,nl,nl[0]*poly[i][0]+nl[1]*poly[i][1]));
  }

  double ez[3]={0,0,1},mz[3]={0,0,-1};
  c.planes.push_back(FPlane(f,ez,w1)),c.planes.push_back(FPlane(f,mz,-w0));
  b.components.push_back(c);
}

std::vector<std::vector<double> > Polygon(const double (*p)[2],int n) {
  std::vector<std::vector<double> > r(n,std::vector<double>(2));
  for (int i=0;i<n;i++) r[i][0]=p[i][0],r[i][1]=p[i][1];
  return r;
}

//an icosphere (nSubdivisions, radius, center) added to the body, with its component
void AddIcosphere(cBody& b,int nSubdivisions,double radius,const double *center) {
  cBody s;
  Icosphere(s,nSubdivisions,1.0);
  int offset=b.vertex.size()/3;

  for (size_t n=0;n<s.vertex.size();n+=3) {
    double x[3];
    for (int d=0;d<3;d++) x[d]=center[d]+radius*(s.vertex[n+d]-kCenter[d]);
    AddVertex(b,x);
  }
  for (size_t n=0;n<s.face.size();n++) b.face.push_back(s.face[n]+offset);

  cComponent c;
  for (size_t n=0;n<s.face.size();n+=3) {
    const double *a=&b.vertex[3*(s.face[n]+offset)],*p=&b.vertex[3*(s.face[n+1]+offset)],*q=&b.vertex[3*(s.face[n+2]+offset)];
    double e1[3],e2[3],nn[3];
    for (int d=0;d<3;d++) e1[d]=p[d]-a[d],e2[d]=q[d]-a[d];
    nn[0]=e1[1]*e2[2]-e1[2]*e2[1],nn[1]=e1[2]*e2[0]-e1[0]*e2[2],nn[2]=e1[0]*e2[1]-e1[1]*e2[0];
    double l=sqrt(Dot(nn,nn)),r[3]={a[0]-center[0],a[1]-center[1],a[2]-center[2]};
    if (Dot(nn,r)<0.0) l=-l;
    cPlane pl;
    for (int d=0;d<3;d++) pl.n[d]=nn[d]/l;
    pl.d=Dot(pl.n,a);
    c.planes.push_back(pl);
  }
  b.components.push_back(c);
}

//a plate thinner than a cell (0.04 < 0.1), tilted
void Plate(cBody& b) {
  const double p[4][2]={{-0.45,-0.45},{0.45,-0.45},{0.45,0.45},{-0.45,0.45}};
  AddPrism(b,MakeFrame(0.12,0.19,0.07,kCenter[0],kCenter[1],kCenter[2]),Polygon(p,4),-0.02,0.02);
}

//small spheres (radius 0.015, 0.015, 0.01; the cells are 0.05 there): the first is 0.0075 from the z-edge x=y=0.3,
//0.25<z<0.3 and spans 0.2615<z<0.2915 (its cap pokes through the edge, which is crossed twice; the corners of the edge
//are outside); the second is 0.009 from the face z=0.5 of the cell [0.35,0.4]^2x[0.45,0.5] and at least 0.0215 from
//the edges of that face (no edge is crossed); the third is inside the cell [0.65,0.7]x[-0.35,-0.3]x[0.05,0.1]
void Bumps(cBody& b) {
  const double c1[3]={0.306,0.3045,0.2765},c2[3]={0.3785,0.373,0.509},c3[3]={0.6762,-0.3243,0.0758};
  AddIcosphere(b,2,0.015,c1);
  AddIcosphere(b,2,0.015,c2);
  AddIcosphere(b,2,0.01,c3);
}

//the L-shaped prism A=[-a,a]x[-a,0] u B=[-a,0]x[0,a] (local u,v) extruded along z, rotated by 45 degrees about z: the
//concave edge is the line x=0.376, y=0.381 (inside the cell column [0.35,0.4]^2), the notch opens towards +y, so the
//x-edge y=0.4 of that column crosses both faces of the notch at x=0.376-+0.019 (inside the edge)
void LCorner(cBody& b) {
  const double a=0.4,w0=-0.33,w1=0.31;
  const cFrame f=MakeFrame(0.0,0.0,0.25*kPi,0.376,0.381,0.0);
  const double L[6][2]={{-a,-a},{a,-a},{a,0.0},{0.0,0.0},{0.0,a},{-a,a}};
  int vlo[6],vhi[6];

  for (int i=0;i<6;i++) {
    double u[3]={L[i][0],L[i][1],w0},x[3];
    FMap(f,u,x); vlo[i]=AddVertex(b,x);
    u[2]=w1; FMap(f,u,x); vhi[i]=AddVertex(b,x);
  }

  //the caps: the L polygon 0..5 (counterclockwise) as the triangles (0,1,2),(0,2,3),(0,3,5),(5,3,4)
  const int cap[4][3]={{0,1,2},{0,2,3},{0,3,5},{5,3,4}};
  for (int k=0;k<4;k++) {
    int t1[3]={vlo[cap[k][0]],vlo[cap[k][1]],vlo[cap[k][2]]},t2[3]={vhi[cap[k][0]],vhi[cap[k][1]],vhi[cap[k][2]]};
    b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);
  }

  for (int i=0;i<6;i++) {
    int j=(i+1)%6;
    int t1[3]={vlo[i],vlo[j],vhi[j]},t2[3]={vlo[i],vhi[j],vhi[i]};
    b.face.insert(b.face.end(),t1,t1+3),b.face.insert(b.face.end(),t2,t2+3);
  }

  double ex[3]={1,0,0},ey[3]={0,1,0},ez[3]={0,0,1},mx[3]={-1,0,0},my[3]={0,-1,0},mz[3]={0,0,-1};
  cComponent A,B;
  A.planes.push_back(FPlane(f,ex,a)),A.planes.push_back(FPlane(f,mx,a)),A.planes.push_back(FPlane(f,ey,0.0)),A.planes.push_back(FPlane(f,my,a));
  B.planes.push_back(FPlane(f,ex,0.0)),B.planes.push_back(FPlane(f,mx,a)),B.planes.push_back(FPlane(f,ey,a)),B.planes.push_back(FPlane(f,my,0.0));
  A.planes.push_back(FPlane(f,ez,w1)),A.planes.push_back(FPlane(f,mz,-w0));
  B.planes.push_back(FPlane(f,ez,w1)),B.planes.push_back(FPlane(f,mz,-w0));
  b.components.push_back(A),b.components.push_back(B);
}

//a prism of the triangle (0.3,0.3), (0.3,0.3)+0.33(cos 20,sin 20), (0.3,0.3)+0.33(cos 70,sin 70) over 0.3<z<0.6 (not
//rotated): its sharp edge lies on the cell edges x=y=0.3, its end vertices on cell corners, its caps on cell faces
void Wedge(cBody& b) {
  const double r=0.33,a1=20.0*kPi/180.0,a2=70.0*kPi/180.0;
  const double p[3][2]={{0.3,0.3},{0.3+r*cos(a1),0.3+r*sin(a1)},{0.3+r*cos(a2),0.3+r*sin(a2)}};
  AddPrism(b,MakeFrame(0.0,0.0,0.0,0.0,0.0,0.0),Polygon(p,3),0.3,0.6);
}
}

double cBody::Phi(const double *x) const {
  double res=1.0E30;

  if (!spheres.empty()) {
    for (size_t n=0;n<spheres.size();n++) {
      double r[3]={x[0]-spheres[n].c[0],x[1]-spheres[n].c[1],x[2]-spheres[n].c[2]};
      res=std::min(res,sqrt(Dot(r,r))-spheres[n].r);
    }
    return res;
  }

  for (size_t c=0;c<components.size();c++) {
    double m=-1.0E30;
    for (size_t n=0;n<components[c].planes.size();n++) m=std::max(m,Dot(components[c].planes[n].n,x)-components[c].planes[n].d);
    res=std::min(res,m);
  }
  return res;
}

bool BuildBody(const std::string& name,cBody& b) {
  b=cBody();
  b.name=name;

  if (name=="sphere1") Icosphere(b,1,0.55);
  else if (name=="sphere2") Icosphere(b,2,0.55);
  else if (name=="sphere3") Icosphere(b,3,0.55);
  else if (name=="cone") ConeCylinder(b,true,32);
  else if (name=="cylinder") ConeCylinder(b,false,32);
  else if (name=="box") AlignedBox(b);
  else if (name=="lshape") LShape(b);
  else if (name=="plate") Plate(b);
  else if (name=="bumps") Bumps(b);
  else if (name=="lcorner") LCorner(b);
  else if (name=="wedge") Wedge(b);
  else if ((name=="asphere")||(name=="aspheres")) {
    cAnalyticSphere s;

    if (name=="asphere") {
      for (int d=0;d<3;d++) s.c[d]=kCenter[d];
      s.r=0.3,b.spheres.push_back(s);
    }
    else {
      //the gap between the first two (0.0175<y<0.0325 on the line of the centers) is inside the y-edges 0<y<0.05 near
      //x=0.3, z=0 (the cells are 0.05 there): those edges cross both spheres; the third (radius 0.015) is 0.0075 from the
      //z-edge x=y=0.3, -0.3<z<-0.25 and its cap pokes through that edge
      const double c[3][4]={{0.304,-0.1025,0.012,0.12},{0.304,0.1525,0.012,0.12},{0.306,0.3045,-0.2765,0.015}};
      for (int n=0;n<3;n++) {for (int d=0;d<3;d++) s.c[d]=c[n][d]; s.r=c[n][3],b.spheres.push_back(s);}
    }

    b.volume=0.0;
    for (size_t n=0;n<b.spheres.size();n++) b.volume+=4.0/3.0*kPi*pow(b.spheres[n].r,3);
    return true;   //the surface is built once the mesh is known (BuildSphereApproximation())
  }
  else return false;

  Finish(b);
  return true;
}

void WriteNastranFile(const cBody& b,const char *fname) {
  FILE *f=fopen(fname,"w");

  fprintf(f,"$ AMPS cut-cell tetrahedra test: %s\n",b.name.c_str());
  for (size_t n=0;n<b.vertex.size();n+=3) fprintf(f,"GRID* %i 0 %.17E %.17E\n* %.17E\n",(int)(n/3+1),b.vertex[n],b.vertex[n+1],b.vertex[n+2]);
  for (size_t n=0;n<b.face.size();n+=3) fprintf(f,"CTRIA3 %i 1 %i %i %i\n",(int)(n/3+1),b.face[n]+1,b.face[n+1]+1,b.face[n+2]+1);
  fprintf(f,"ENDDATA\n");
  fclose(f);
}

void BuildSphereApproximation(cBody& b,const std::vector<int>& N) {
  std::map<std::vector<long>,int> index;

  b.vertex.clear(),b.face.clear(),b.components.clear();

  for (size_t ns=0;ns<b.spheres.size();ns++) {
    const cAnalyticSphere& s=b.spheres[ns];
    const int n=N[ns];

    //the vertex of the integer point (ix,iy,iz) of the cube [0,n]^3 (the same expression as CutCell::cSphereApproximation)
    auto Vertex=[&](int ix,int iy,int iz) {
      std::vector<long> key(4);
      key[0]=ns,key[1]=ix,key[2]=iy,key[3]=iz;
      std::map<std::vector<long>,int>::iterator it=index.find(key);
      if (it!=index.end()) return it->second;

      double p[3]={(double)(2*ix-n)/n,(double)(2*iy-n)/n,(double)(2*iz-n)/n},x[3];
      double l=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
      for (int d=0;d<3;d++) x[d]=s.c[d]+s.r*p[d]/l;
      return index[key]=AddVertex(b,x);
    };

    for (int face=0;face<6;face++) for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
      int k=face/2,a=(k+1)%3,c=(k+2)%3,v[4];
      const int di[4]={0,1,1,0},dj[4]={0,0,1,1};

      for (int m=0;m<4;m++) {
        int I[3];
        I[k]=(face%2==1) ? n : 0,I[a]=i+di[m],I[c]=j+dj[m];
        v[m]=Vertex(I[0],I[1],I[2]);
      }

      //the diagonal (i,j)-(i+1,j+1), as in the library
      const int tri[2][3]={{v[0],v[1],v[2]},{v[0],v[2],v[3]}};

      for (int t=0;t<2;t++) {
        b.face.insert(b.face.end(),tri[t],tri[t]+3);

        //the component: the tetrahedron (center, facet)
        const double *q[3]={&b.vertex[3*tri[t][0]],&b.vertex[3*tri[t][1]],&b.vertex[3*tri[t][2]]};
        cComponent comp;

        for (int e=0;e<4;e++) {
          //face e of the tetrahedron: e<3 the side through the center and the edge (q[e],q[e+1]), e=3 the facet
          const double *p0=(e<3) ? s.c : q[0],*p1=(e<3) ? q[e] : q[1],*p2=(e<3) ? q[(e+1)%3] : q[2],*opp=(e<3) ? q[(e+2)%3] : s.c;
          double e1[3],e2[3],nn[3],l;
          cPlane pl;

          for (int d=0;d<3;d++) e1[d]=p1[d]-p0[d],e2[d]=p2[d]-p0[d];
          nn[0]=e1[1]*e2[2]-e1[2]*e2[1],nn[1]=e1[2]*e2[0]-e1[0]*e2[2],nn[2]=e1[0]*e2[1]-e1[1]*e2[0];
          l=sqrt(Dot(nn,nn));
          if (nn[0]*(opp[0]-p0[0])+nn[1]*(opp[1]-p0[1])+nn[2]*(opp[2]-p0[2])>0.0) l=-l;
          for (int d=0;d<3;d++) pl.n[d]=nn[d]/l;
          pl.d=Dot(pl.n,p0);
          comp.planes.push_back(pl);
        }

        comp.hasBox=true;
        for (int d=0;d<3;d++) {
          comp.lo[d]=std::min(s.c[d],std::min(q[0][d],std::min(q[1][d],q[2][d])));
          comp.hi[d]=std::max(s.c[d],std::max(q[0][d],std::max(q[1][d],q[2][d])));
        }
        b.components.push_back(comp);
      }
    }
  }

  //the orientation (outward: the components are the star tetrahedra) and the volume of the polyhedra
  Finish(b);
}

//---------------------------------- exact convex geometry ----------------------------------
namespace {
const double kTol=1.0E-13;   //relative to the size of the coordinates (|x|<=1 in the test)

std::vector<double> P(double x,double y,double z) {std::vector<double> p(3); p[0]=x,p[1]=y,p[2]=z; return p;}

double SDist(const std::vector<double>& x,const cPlane& h) {return h.n[0]*x[0]+h.n[1]*x[1]+h.n[2]*x[2]-h.d;}

//order the points of a convex planar polygon with the normal n
void Order(cPolygon3D& p,const double *n) {
  if (p.size()<3) return;
  double c[3]={0,0,0},e0[3],e1[3],l;
  for (size_t i=0;i<p.size();i++) for (int d=0;d<3;d++) c[d]+=p[i][d]/p.size();
  if (fabs(n[0])<0.9) e0[0]=0.0,e0[1]=n[2],e0[2]=-n[1];
  else e0[0]=-n[2],e0[1]=0.0,e0[2]=n[0];
  l=sqrt(Dot(e0,e0));
  for (int d=0;d<3;d++) e0[d]/=l;
  e1[0]=n[1]*e0[2]-n[2]*e0[1],e1[1]=n[2]*e0[0]-n[0]*e0[2],e1[2]=n[0]*e0[1]-n[1]*e0[0];
  std::vector<std::pair<double,std::vector<double> > > a;
  for (size_t i=0;i<p.size();i++) {
    double r[3]={p[i][0]-c[0],p[i][1]-c[1],p[i][2]-c[2]};
    a.push_back(std::make_pair(atan2(Dot(r,e1),Dot(r,e0)),p[i]));
  }
  std::sort(a.begin(),a.end());
  for (size_t i=0;i<a.size();i++) p[i]=a[i].second;
}
}

cPolyhedron Box(const double *lo,const double *hi) {
  cPolyhedron b;
  static const int fc[6][4]={{0,2,6,4},{1,5,7,3},{0,4,5,1},{2,3,7,6},{0,1,3,2},{4,6,7,5}};
  for (int i=0;i<6;i++) {
    cPolygon3D f;
    for (int j=0;j<4;j++) {int k=fc[i][j]; f.push_back(P((k&1) ? hi[0] : lo[0],(k&2) ? hi[1] : lo[1],(k&4) ? hi[2] : lo[2]));}
    b.faces.push_back(f);
  }
  return b;
}

cPolyhedron Tetrahedron(const double x[4][3]) {
  cPolyhedron t;
  static const int fc[4][3]={{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
  for (int i=0;i<4;i++) {
    cPolygon3D f;
    for (int j=0;j<3;j++) f.push_back(P(x[fc[i][j]][0],x[fc[i][j]][1],x[fc[i][j]][2]));
    t.faces.push_back(f);
  }
  return t;
}

cPolygon3D ClipPolygon(const cPolygon3D& p,const cPlane& h) {
  cPolygon3D r;
  size_t n=p.size();
  for (size_t i=0;i<n;i++) {
    const std::vector<double>& a=p[i];
    const std::vector<double>& b=p[(i+1)%n];
    double sa=SDist(a,h),sb=SDist(b,h);
    if (fabs(sa)<kTol) sa=0.0;
    if (fabs(sb)<kTol) sb=0.0;
    if (sa<=0.0) r.push_back(a);
    if (((sa<0.0)&&(sb>0.0))||((sa>0.0)&&(sb<0.0))) {
      double t=sa/(sa-sb);
      r.push_back(P(a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1]),a[2]+t*(b[2]-a[2])));
    }
  }
  if (r.size()<3) r.clear();
  return r;
}

cPolyhedron Clip(const cPolyhedron& p,const cPlane& h) {
  cPolyhedron r;
  cPolygon3D cap;
  double smin=1.0E300,smax=-1.0E300;

  //entirely on one side (a face may lie in the plane)
  for (size_t i=0;i<p.faces.size();i++) for (size_t j=0;j<p.faces[i].size();j++) {
    double s=SDist(p.faces[i][j],h);
    smin=std::min(smin,s),smax=std::max(smax,s);
  }
  if (smax<=kTol) return p;
  if (smin>=-kTol) return r;

  for (size_t i=0;i<p.faces.size();i++) {
    cPolygon3D f=ClipPolygon(p.faces[i],h);
    if (f.size()>=3) r.faces.push_back(f);
    for (size_t j=0;j<p.faces[i].size();j++) {
      const std::vector<double>& a=p.faces[i][j];
      const std::vector<double>& b=p.faces[i][(j+1)%p.faces[i].size()];
      double sa=SDist(a,h),sb=SDist(b,h);
      if (fabs(sa)<kTol) {cap.push_back(a); continue;}
      if (((sa<-kTol)&&(sb>kTol))||((sa>kTol)&&(sb<-kTol))) {
        double t=sa/(sa-sb);
        cap.push_back(P(a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1]),a[2]+t*(b[2]-a[2])));
      }
    }
  }

  //unique points of the cap
  cPolygon3D u;
  for (size_t i=0;i<cap.size();i++) {
    bool found=false;
    for (size_t j=0;(j<u.size())&&!found;j++) {
      double d2=0.0;
      for (int d=0;d<3;d++) d2+=pow(cap[i][d]-u[j][d],2);
      if (d2<kTol*kTol) found=true;
    }
    if (!found) u.push_back(cap[i]);
  }

  if (u.size()>=3) {Order(u,h.n); r.faces.push_back(u);}
  if (r.faces.size()<4) r.faces.clear();
  return r;
}

double Area(const cPolygon3D& p) {
  double a[3]={0,0,0};
  for (size_t i=1;i+1<p.size();i++) {
    double e1[3],e2[3];
    for (int d=0;d<3;d++) e1[d]=p[i][d]-p[0][d],e2[d]=p[i+1][d]-p[0][d];
    a[0]+=e1[1]*e2[2]-e1[2]*e2[1],a[1]+=e1[2]*e2[0]-e1[0]*e2[2],a[2]+=e1[0]*e2[1]-e1[1]*e2[0];
  }
  return 0.5*sqrt(Dot(a,a));
}

double Volume(const cPolyhedron& p) {
  if (p.faces.empty()) return 0.0;
  double c[3]={0,0,0},v=0.0;
  int n=0;
  for (size_t i=0;i<p.faces.size();i++) for (size_t j=0;j<p.faces[i].size();j++,n++) for (int d=0;d<3;d++) c[d]+=p.faces[i][j][d];
  for (int d=0;d<3;d++) c[d]/=n;
  for (size_t i=0;i<p.faces.size();i++) {
    const cPolygon3D& f=p.faces[i];
    for (size_t j=1;j+1<f.size();j++) {
      double e1[3],e2[3],e3[3];
      for (int d=0;d<3;d++) e1[d]=f[0][d]-c[d],e2[d]=f[j][d]-c[d],e3[d]=f[j+1][d]-c[d];
      v+=fabs(e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])-e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])+e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]))/6.0;
    }
  }
  return v;
}

cBody RestrictToBox(const cBody& body,const double *lo,const double *hi) {
  cBody r;
  r.name=body.name,r.volume=body.volume;

  for (size_t c=0;c<body.components.size();c++) {
    cComponent k;
    bool reaches=true;

    if (body.components[c].hasBox) {
      for (int d=0;d<3;d++) if ((body.components[c].lo[d]>hi[d]+kTol)||(body.components[c].hi[d]<lo[d]-kTol)) reaches=false;
      if (!reaches) continue;
    }

    for (size_t n=0;(n<body.components[c].planes.size())&&reaches;n++) {
      const cPlane& h=body.components[c].planes[n];
      double smin=0.0,smax=0.0;
      for (int d=0;d<3;d++) {
        double a=h.n[d]*lo[d],b=h.n[d]*hi[d];
        smin+=std::min(a,b),smax+=std::max(a,b);
      }
      smin-=h.d,smax-=h.d;
      if (smin>kTol) reaches=false;           //the box is outside this face plane
      else if (smax>-kTol) k.planes.push_back(h);   //the plane cuts (or touches) the box
    }

    if (reaches) {
      k.hasBox=body.components[c].hasBox;
      for (int d=0;d<3;d++) k.lo[d]=body.components[c].lo[d],k.hi[d]=body.components[c].hi[d];
      r.components.push_back(k);
    }
  }

  return r;
}

namespace {
//the component against the points x (the vertices of a convex set): 0 the set is outside it (the result of the clipping
//is empty), 1 inside (clipping changes nothing), 2 it must be clipped by the planes in 'cut'
int Classify(const cComponent& k,const std::vector<const std::vector<double>*>& x,const double *lo,const double *hi,std::vector<int>& cut) {
  cut.clear();
  if (k.hasBox) for (int d=0;d<3;d++) if ((k.lo[d]>hi[d]+kTol)||(k.hi[d]<lo[d]-kTol)) return 0;

  for (size_t n=0;n<k.planes.size();n++) {
    double smin=1.0E300,smax=-1.0E300;
    for (size_t i=0;i<x.size();i++) {double s=SDist(*x[i],k.planes[n]); smin=std::min(smin,s),smax=std::max(smax,s);}
    if (smax<=kTol) continue;       //inside the half-space (a set in the plane is kept, as in Clip() and ClipPolygon())
    if (smin>=-kTol) return 0;
    cut.push_back(n);
  }

  return cut.empty() ? 1 : 2;
}
}

double IntersectionVolume(const cPolyhedron& p,const cBody& body) {
  std::vector<const std::vector<double>*> x;
  std::vector<int> cut;
  double v=0.0,vp=-1.0,lo[3]={1.0E300,1.0E300,1.0E300},hi[3]={-1.0E300,-1.0E300,-1.0E300};

  for (size_t i=0;i<p.faces.size();i++) for (size_t j=0;j<p.faces[i].size();j++) {
    x.push_back(&p.faces[i][j]);
    for (int d=0;d<3;d++) lo[d]=std::min(lo[d],p.faces[i][j][d]),hi[d]=std::max(hi[d],p.faces[i][j][d]);
  }
  if (x.empty()) return 0.0;

  for (size_t c=0;c<body.components.size();c++) {
    switch (Classify(body.components[c],x,lo,hi,cut)) {
    case 0:
      break;
    case 1:
      if (vp<0.0) vp=Volume(p);
      v+=vp;
      break;
    default: {
      cPolyhedron q=p;
      for (size_t n=0;(n<cut.size())&&!q.faces.empty();n++) q=Clip(q,body.components[c].planes[cut[n]]);
      v+=Volume(q);
    }
    }
  }
  return v;
}

double IntersectionArea(const cPolygon3D& p,const cBody& body) {
  std::vector<const std::vector<double>*> x;
  std::vector<int> cut;
  double a=0.0,ap=-1.0,lo[3]={1.0E300,1.0E300,1.0E300},hi[3]={-1.0E300,-1.0E300,-1.0E300};

  if (p.size()<3) return 0.0;
  for (size_t i=0;i<p.size();i++) {
    x.push_back(&p[i]);
    for (int d=0;d<3;d++) lo[d]=std::min(lo[d],p[i][d]),hi[d]=std::max(hi[d],p[i][d]);
  }

  for (size_t c=0;c<body.components.size();c++) {
    switch (Classify(body.components[c],x,lo,hi,cut)) {
    case 0:
      break;
    case 1:
      if (ap<0.0) ap=Area(p);
      a+=ap;
      break;
    default: {
      cPolygon3D q=p;
      for (size_t n=0;(n<cut.size())&&!q.empty();n++) q=ClipPolygon(q,body.components[c].planes[cut[n]]);
      if (q.size()>=3) a+=Area(q);
    }
    }
  }
  return a;
}

}

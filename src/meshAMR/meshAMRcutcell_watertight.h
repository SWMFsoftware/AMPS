/*
 * meshAMRcutcell_watertight.h: watertight intersection of a line with a triangle (header only)
 *
 * CutCell::cTriangleFace::RayIntersection() decides whether a hit point lies in a face from barycentric
 * coordinates that every face computes for itself, with a tolerance (EPS converted to a barycentric margin
 * and xLocalMargin).  Two faces that share an edge therefore round differently: a line through the edge can be
 * rejected by both (a gap through which a particle enters the body) or accepted by both, and the result
 * depends on EPS.  The movers built on it add further tolerances (hits closer than EPS ignored, the last face
 * excluded, particles pushed off the surface).
 *
 * The test here is the one of Woop, Benthin and Wald (2013): the triangle is transformed into a frame in which
 * the line is the z axis (a shear by the line direction), and the line passes through the triangle when the
 * three 2D edge functions E(P,Q)=Px*Qy-Py*Qx of its edges have the same sign.  A vertex has the same transformed
 * coordinates in every face that uses it, and the edge function of a shared edge is evaluated with the same
 * operands in both faces with the order of the endpoints swapped.  The signs are computed exactly
 * (Orient2DSign: a floating-point filter and, when it cannot decide, an exact expansion of the 2x2 determinant
 * of the rounded coordinates), so E(Q,P)=-E(P,Q) holds exactly and the decisions of all faces are those of one
 * consistent planar arrangement: a line that crosses a closed surface passes through at least one of its faces
 * - also through an edge or a vertex - so there are no gaps.  A line through an edge or a vertex (an edge
 * function exactly zero) is accepted by every face that contains the point (closed triangles); callers that need
 * the first crossing take the smallest parameter.
 *
 * No AMPS dependency besides <cmath>.
 */

#ifndef _MESHAMR_CUTCELL_WATERTIGHT_H_
#define _MESHAMR_CUTCELL_WATERTIGHT_H_

#include <cmath>
#include <utility>

namespace CutCell {
namespace Watertight {

//error-free transformations (Shewchuk 1997); they require IEEE double arithmetic without -ffast-math
inline void TwoSum(double a,double b,double& x,double& y) {
  x=a+b;
  double bv=x-a,av=x-bv,br=b-bv,ar=a-av;
  y=ar+br;
}

inline void TwoDiff(double a,double b,double& x,double& y) {
  x=a-b;
  double bv=a-x,av=x+bv,br=bv-b,ar=a-av;
  y=ar+br;
}

inline void TwoProduct(double a,double b,double& x,double& y) {
  x=a*b;
  y=std::fma(a,b,-x);
}

//exact sign of px*qy-py*qx for doubles: +1, -1 or 0
inline int Orient2DSign(double px,double py,double qx,double qy) {
  const double l=px*qy,r=py*qx,det=l-r;
  const double sum=std::fabs(l)+std::fabs(r);
  //(3+16 eps) eps, eps=2^-53: |det-exact| <= errbound (Shewchuk's ccwerrboundA)
  const double errbound=3.3306690738754716e-16*sum;

  if (det>errbound) return 1;
  if (-det>errbound) return -1;
  if (sum==0.0) return 0;

  //exact: (lh+ll)-(rh+rl) as a nonoverlapping expansion x3>x2>x1>x0 (Two_Two_Diff)
  double lh,ll,rh,rl;
  TwoProduct(px,qy,lh,ll);
  TwoProduct(py,qx,rh,rl);

  double i,j,k0,x0,x1,x2,x3;
  TwoDiff(ll,rl,i,x0);
  TwoSum(lh,i,j,k0);
  TwoDiff(k0,rh,i,x1);
  TwoSum(j,i,x3,x2);

  if (x3!=0.0) return (x3>0.0) ? 1 : -1;
  if (x2!=0.0) return (x2>0.0) ? 1 : -1;
  if (x1!=0.0) return (x1>0.0) ? 1 : -1;
  if (x0!=0.0) return (x0>0.0) ? 1 : -1;
  return 0;
}

//The line x(t)=x0+t*l and the closed triangle (a,b,c).  Returns true when the line passes through the
//triangle and is not parallel to its plane; t is then the parameter of the crossing (any real value: the
//caller restricts the range).  The decision is watertight (see above); t is the usual floating-point estimate.
inline bool LineTriangleIntersection(const double *x0,const double *l,const double *a,const double *b,const double *c,double& t) {
  const double al[3]={std::fabs(l[0]),std::fabs(l[1]),std::fabs(l[2])};
  int kz=(al[0]>=al[1]) ? ((al[0]>=al[2]) ? 0 : 2) : ((al[1]>=al[2]) ? 1 : 2);

  if (!(al[kz]>0.0)) return false;

  int kx=(kz+1)%3,ky=(kx+1)%3;
  if (l[kz]<0.0) std::swap(kx,ky);

  const double Sx=l[kx]/l[kz],Sy=l[ky]/l[kz],Sz=1.0/l[kz];
  const double A[3]={a[0]-x0[0],a[1]-x0[1],a[2]-x0[2]};
  const double B[3]={b[0]-x0[0],b[1]-x0[1],b[2]-x0[2]};
  const double C[3]={c[0]-x0[0],c[1]-x0[1],c[2]-x0[2]};

  const double Ax=A[kx]-Sx*A[kz],Ay=A[ky]-Sy*A[kz];
  const double Bx=B[kx]-Sx*B[kz],By=B[ky]-Sy*B[kz];
  const double Cx=C[kx]-Sx*C[kz],Cy=C[ky]-Sy*C[kz];

  //edge functions of the edges b->c, c->a, a->b
  const int sU=Orient2DSign(Cx,Cy,Bx,By);
  const int sV=Orient2DSign(Ax,Ay,Cx,Cy);
  const int sW=Orient2DSign(Bx,By,Ax,Ay);

  if (((sU<0)||(sV<0)||(sW<0))&&((sU>0)||(sV>0)||(sW>0))) return false;
  if ((sU==0)&&(sV==0)&&(sW==0)) return false;  //the line lies in the plane of the triangle

  const double U=Cx*By-Cy*Bx,V=Ax*Cy-Ay*Cx,W=Bx*Ay-By*Ax;
  const double det=U+V+W;
  if (!(det!=0.0)||!std::isfinite(det)) return false;

  const double T=U*(Sz*A[kz])+V*(Sz*B[kz])+W*(Sz*C[kz]);
  t=T/det;
  return std::isfinite(t);
}

//the point q of the closed triangle (a,b,c) nearest to p (Ericson, Real-Time Collision Detection, 5.1.5)
inline void ClosestPointOnTriangle(const double *p,const double *a,const double *b,const double *c,double *q) {
  double ab[3],ac[3],ap[3];
  for (int i=0;i<3;i++) ab[i]=b[i]-a[i],ac[i]=c[i]-a[i],ap[i]=p[i]-a[i];
  auto dot=[](const double *u,const double *v) {return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];};

  const double d1=dot(ab,ap),d2=dot(ac,ap);
  if ((d1<=0.0)&&(d2<=0.0)) {for (int i=0;i<3;i++) q[i]=a[i]; return;}

  double bp[3];
  for (int i=0;i<3;i++) bp[i]=p[i]-b[i];
  const double d3=dot(ab,bp),d4=dot(ac,bp);
  if ((d3>=0.0)&&(d4<=d3)) {for (int i=0;i<3;i++) q[i]=b[i]; return;}

  const double vc=d1*d4-d3*d2;
  if ((vc<=0.0)&&(d1>=0.0)&&(d3<=0.0)) {const double v=d1/(d1-d3); for (int i=0;i<3;i++) q[i]=a[i]+v*ab[i]; return;}

  double cp[3];
  for (int i=0;i<3;i++) cp[i]=p[i]-c[i];
  const double d5=dot(ab,cp),d6=dot(ac,cp);
  if ((d6>=0.0)&&(d5<=d6)) {for (int i=0;i<3;i++) q[i]=c[i]; return;}

  const double vb=d5*d2-d1*d6;
  if ((vb<=0.0)&&(d2>=0.0)&&(d6<=0.0)) {const double w=d2/(d2-d6); for (int i=0;i<3;i++) q[i]=a[i]+w*ac[i]; return;}

  const double va=d3*d6-d5*d4;
  if ((va<=0.0)&&(d4-d3>=0.0)&&(d5-d6>=0.0)) {
    const double w=(d4-d3)/((d4-d3)+(d5-d6));
    for (int i=0;i<3;i++) q[i]=b[i]+w*(c[i]-b[i]);
    return;
  }

  const double denom=1.0/(va+vb+vc),v=vb*denom,w=vc*denom;
  for (int i=0;i<3;i++) q[i]=a[i]+ab[i]*v+ac[i]*w;
}

}  // namespace Watertight
}  // namespace CutCell

#endif

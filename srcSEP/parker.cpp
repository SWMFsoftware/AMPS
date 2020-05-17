/*
 * parker.cpp
 *
 *  Created on: May 16, 2020
 *      Author: vtenishe
 */


#include "sep.h"

void SEP::ParkerSpiral::GetB(double* B,double *x,double u_sw) {
  double B0=1.83E-6,R0=10.0*_RADIUS_(_SUN_);
  double omega=2*Pi/(25.4*24*3600.0);

  double e_r[3],e_phi[3],e_z[3]={0.0,0.0,1.0};
  int idim;

  for (idim=0;idim<3;idim++) e_r[idim]=x[idim];

  Vector3D::CrossProduct(e_phi,e_r,e_z);

  Vector3D::Normalize(e_r);
  Vector3D::Normalize(e_phi);


  double aa,bb,r,theta;


  r=Vector3D::Length(x);
  aa=B0*pow(R0/r,2);

  theta=acos(Vector3D::DotProduct(e_z,e_r));
  bb=aa*(r-R0)*omega*sin(theta)/u_sw;

  for (idim=0;idim<3;idim++) B[idim]=aa*e_r[idim]-bb*e_phi[idim];
}


void SEP::ParkerSpiral::CreateFileLine(list<SEP::cFieldLine> *field_line,double *xstart,double length_rsun) {
  double l[3],dl;
  int idim;
  SEP::cFieldLine p;

  const int npoints=4000;

  double u_sw=400.0E3;

  dl=length_rsun/npoints;

  for (idim=0;idim<3;idim++) p.x[idim]=xstart[3];

  field_line->push_back(p);


  for (int ipoint=0;ipoint<npoints;ipoint++) {
    GetB(l,p.x,u_sw);
    Vector3D::Normalize(l,dl);

    for (idim=0;idim<3;idim++) p.x[idim]+=l[idim];

    field_line->push_back(p);
  }
}

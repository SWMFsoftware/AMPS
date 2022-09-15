#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "constants.h"
#include "sep.h"
#include "array_1d.h"
#include "specfunc.h"

#include "tests.h"

using namespace std;

void TestManager() {
  DxxTest();
  ParkerModelMoverTest_const_plasma_field();
}

void DiffusionCoefficient_const(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment) {
   D=1.0,dD_dmu=0.0;
}


void DxxTest() {
  double c,D,dDxx_dx,FieldLineCoord,v;
  int spec,iFieldLine;
  PIC::FieldLine::cFieldLineSegment *Segment; 

  const bool _pass=true;
  const bool _fail=false;

  bool res=_pass;
  
  //1. constant D_mu_mu 
  v=1.0E6; 
  spec=0;
  FieldLineCoord=1.5;

  Segment=PIC::FieldLine::FieldLinesAll[0].GetFirstSegment();  
  iFieldLine=0; 

  SEP::Diffusion::GetPitchAngleDiffusionCoefficient=DiffusionCoefficient_const;
  SEP::Diffusion::GetDxx(D,dDxx_dx,v,spec,FieldLineCoord,Segment,iFieldLine); 

  if ((fabs(dDxx_dx)>1.0E-5)||(c=fabs(1.0-D/(v*v/8.0*16.0/15.0))>1.0E-5)) {
    res=_fail;

    //debugging:
    c=v*v/8.0*16.0/15.0;
    c-=D; 
  }
}
  
//====================================================================================
void ParkerModelMoverTest_const_plasma_field() {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;
  
  //generate a particle and call the mover; 
  //because of the scattering dX should exibit "normal" distribution
  const int nTotalTests=4000000;
  double v[3]={1.0E6,0.0,0.0},x_new[3],dx;
  double v_comp=sqrt(2.0)*Vector3D::Length(v);

  double dxSampleMin=-1.0E7;
  double dxSampleMax=1.0E7;
  double nSampleIntervals=100;
  double dxSampleStep=(dxSampleMax-dxSampleMin)/nSampleIntervals;
  int iSample,i;

  array_1d<double> SamplingBuffer(nSampleIntervals);
  SamplingBuffer=0.0;

  double dtTotal=3.0;

  //determine the particle location and the starting node 
  double xParticleCoordinate=67.5;
  double x0[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

  long int ptr;
  double xLocal;

  FL::FieldLinesAll[0].GetCartesian(x0,xParticleCoordinate);
  node=PIC::Mesh::mesh->findTreeNode(x0);

  ptr=PB::GetNewParticle();
  PB::SetI(0,ptr);

  for (int ntest=0;ntest<nTotalTests;ntest++) {
    PB::SetVParallel(v_comp,ptr);
    PB::SetVNormal(v_comp,ptr);
    PB::SetFieldLineCoord(xParticleCoordinate,ptr); 

    SEP::ParticleMover_ParkerEquation(ptr,dtTotal,node); 

    xLocal=PB::GetFieldLineCoord(ptr);
    FL::FieldLinesAll[0].GetCartesian(x_new,xLocal); 
    dx=0.0;

    for (int i=0;i<3;i++) {
      double t=x_new[i]-x0[i]; 

      dx+=t*t; 
    }

    dx=sqrt(dx);
    if (xLocal<xParticleCoordinate) dx*=-1.0;

    iSample=(int)((dx-dxSampleMin)/dxSampleStep);

    if ((iSample>=0)&&(iSample<nSampleIntervals-1)) {
      SamplingBuffer(iSample)+=1.0;
    } 
  }

  //output sampled data
  ofstream fout("dxParker.dat");
  
  for (i=0;i<nSampleIntervals;i++) {
    fout << dxSampleMin+(i+0.5)*dxSampleStep << "   " << SamplingBuffer(i) << endl;
  } 

  fout.close();
}  





















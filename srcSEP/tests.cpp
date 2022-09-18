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

  FTE_Convectoin();


  DxxTest();
  ParkerModelMoverTest_convection();
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


//====================================================================================
void ParkerModelMoverTest_convection() {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  struct cVertexData {
    double Vsw,DensityOld,DensityCurrent,v[3];
  };
  
  const bool _pass=true;
  const bool _fail=false;

  bool res=_pass;

  list <cVertexData> VertexData;
  double DensityOld=1.0,DensityCurrent=4.0;
  double SolarWindVelocityOld[3]={1.0E3,0.0,0.0};
  double SolarWindVelocityCurrent[3]={4.0E3,0.0,0.0};
  double dtTotal=1.0;
  
  auto DiffusionCoeffcient=SEP::Diffusion::GetPitchAngleDiffusionCoefficient;
  SEP::Diffusion::GetPitchAngleDiffusionCoefficient=NULL;   
  
  //determine the particle location and the starting node 
  double xTestSegment=67.5;
  int iSegment=(int)xTestSegment;
  double s0,x0[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  auto Segment=FL::FieldLinesAll[0].GetSegment(xTestSegment);; 

  auto Vertex0=Segment->GetBegin(); 
  auto Vertex1=Vertex0->GetNext();

  int nTotalTests=10000;
  double logEmin=log(100.0*KeV2J);
  double logEmax=log(10.0*MeV2J); 

  long int ptr;
  double mu,vNorm,vParallel,e,speed,s1,vNormInit,vParallelInit;
  double mass=PIC::MolecularData::GetMass(0);

  ptr=PB::GetNewParticle();
  PB::SetI(0,ptr);

  bool shock_reached=false;

  for (auto Vertex=FL::FieldLinesAll[0].GetFirstVertex();Vertex!=NULL;Vertex=Vertex->GetNext()) {
    cVertexData t;

    Vertex->GetDatum(FL::DatumAtVertexPlasmaDensity,&t.DensityCurrent);
    Vertex->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&t.DensityOld);
    Vertex->GetPlasmaVelocity(t.v);

    VertexData.push_back(t);

    if (shock_reached==false) {
      Vertex->SetDatum(FL::DatumAtVertexPlasmaDensity,DensityCurrent);
      Vertex->SetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,DensityOld);
      Vertex->SetPlasmaVelocity(SolarWindVelocityCurrent);

      if (Vertex==Vertex0) shock_reached=true;
    }
    else {
      Vertex->SetDatum(FL::DatumAtVertexPlasmaDensity,DensityOld);
      Vertex->SetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,DensityOld);
      Vertex->SetPlasmaVelocity(SolarWindVelocityOld);
    }
  }

  for (int ntest=0;ntest<nTotalTests;ntest++) {
    mu=-1.0+2.0*rnd();
    e=exp(logEmin+rnd()*(logEmax-logEmin));

    speed=Relativistic::E2Speed(e,mass); 
    vParallel=speed*mu;
    vNorm=speed*sqrt(1.0-mu*mu);
    
    vParallelInit=vParallel,vNormInit=vNorm;

    PB::SetVParallel(vParallel,ptr);
    PB::SetVNormal(vNorm,ptr);

    s0=iSegment+rnd();
    PB::SetFieldLineCoord(s0,ptr);
    
    FL::FieldLinesAll[0].GetCartesian(x0,s0);
    node=PIC::Mesh::mesh->findTreeNode(x0);

    SEP::ParticleMover_ParkerEquation(ptr,dtTotal,node);

    //check the new particle location: it sould no change 
    s1=PB::GetFieldLineCoord(ptr);
    
    if (s1!=s0) {
      res=_fail;
    }

    //check the new particle velocity: it sould changes as in Sokolov-2004-AJ
    double p0,p1,p1_theory;

    p0=Relativistic::Speed2Momentum(speed,mass);

    vParallel=PB::GetVParallel(ptr);
    vNorm=PB::GetVNormal(ptr); 
    p1=Relativistic::Speed2Momentum(sqrt(vParallel*vParallel+vNorm*vNorm),mass);

    //calculate div(vSW) : Dln(Rho)=-div(vSW)*dt
    double w0,w1,d_ln_rho_dt;

    w1=s0-((int)s0);
    w0=1.0-w1;

    d_ln_rho_dt=log((w0*DensityCurrent+w1*DensityOld)/DensityOld)/dtTotal;  
    p1_theory=p0*exp(d_ln_rho_dt*dtTotal/3.0);

    if (fabs(1.0-p1_theory/p1)>1.0E-5) {
      res=_fail;
      
      PB::SetVParallel(vParallel,ptr);
      PB::SetVNormal(vNorm,ptr);
      SEP::ParticleMover_ParkerEquation(ptr,dtTotal,node);
    } 
  }

  //return the original parameters of the field line
  auto p=VertexData.begin();
  
  for (auto Vertex=FL::FieldLinesAll[0].GetFirstVertex();Vertex!=NULL;p++,Vertex=Vertex->GetNext()) {
    Vertex->SetDatum(FL::DatumAtVertexPlasmaDensity,p->DensityCurrent);
    Vertex->SetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,p->DensityOld);
    Vertex->SetPlasmaVelocity(p->v);
  }
  
  SEP::Diffusion::GetPitchAngleDiffusionCoefficient=DiffusionCoeffcient;
}

   

void FTE_Convectoin() {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  struct cVertexData {
    double Vsw,DensityOld,DensityCurrent,v[3];
  };

  const bool _pass=true;
  const bool _fail=false;

  bool res=_pass;

  list <cVertexData> VertexData;
  double DensityOld=1.0,DensityCurrent=4.0;
  double SolarWindVelocityOld[3]={1.0E3,0.0,0.0};
  double SolarWindVelocityCurrent[3]={4.0E3,0.0,0.0};
  double dtTotal=1.0;

  auto DiffusionCoeffcient=SEP::Diffusion::GetPitchAngleDiffusionCoefficient;
  SEP::Diffusion::GetPitchAngleDiffusionCoefficient=NULL;

  //determine the particle location and the starting node
  double xTestSegment=67.5;
  int iSegment=(int)xTestSegment;
  double s0,x0[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  auto Segment=FL::FieldLinesAll[0].GetSegment(xTestSegment);;

  auto Vertex0=Segment->GetBegin();
  auto Vertex1=Vertex0->GetNext();

  int nTotalTests=10000;
  double logEmin=log(100.0*KeV2J);
  double logEmax=log(10.0*MeV2J);

  long int ptr;
  double mu,vNorm,vParallel,e,speed,s1,vNormInit,vParallelInit;
  double mass=PIC::MolecularData::GetMass(0);

  ptr=PB::GetNewParticle();
  PB::SetI(0,ptr);
  PB::SetIndividualStatWeightCorrection(1.0,ptr);

  for (auto Vertex=FL::FieldLinesAll[0].GetFirstVertex();Vertex!=NULL;Vertex=Vertex->GetNext()) {
    cVertexData t;

    Vertex->GetDatum(FL::DatumAtVertexPlasmaDensity,&t.DensityCurrent);
    Vertex->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&t.DensityOld);
    Vertex->GetPlasmaVelocity(t.v);

    VertexData.push_back(t);

    Vertex->SetDatum(FL::DatumAtVertexPlasmaDensity,DensityCurrent);
    Vertex->SetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,DensityOld);
    Vertex->SetPlasmaVelocity(SolarWindVelocityCurrent);
  }


  for (int ntest=0;ntest<nTotalTests;ntest++) {
    mu=-1.0+2.0*rnd();
    e=exp(logEmin+rnd()*(logEmax-logEmin));

    speed=Relativistic::E2Speed(e,mass);
    vParallel=speed*mu;
    vNorm=speed*sqrt(1.0-mu*mu);

    vParallelInit=vParallel,vNormInit=vNorm;

    PB::SetVParallel(vParallel,ptr);
    PB::SetVNormal(vNorm,ptr);

    s0=iSegment+rnd();
    PB::SetFieldLineCoord(s0,ptr);

    FL::FieldLinesAll[0].GetCartesian(x0,s0);
    node=PIC::Mesh::mesh->findTreeNode(x0);

    SEP::ParticleMover_He_2011_AJ(ptr,dtTotal,node);

    //check the new particle location: it sould no change
    s1=PB::GetFieldLineCoord(ptr);

    double vParallel_new=PB::GetVParallel(ptr);
    double vNorm_new=PB::GetVNormal(ptr);
    double s,x1[3];
    int idim;

    FL::FieldLinesAll[0].GetCartesian(x1,s1);

    for (s=0.0,idim=0;idim<3;idim++) {
      double t=x0[idim]-x1[idim];
   
      s+=t*t;
    }

    s=sqrt(s); 

    double diff=1.0-s/(dtTotal*fabs(vParallel));

    if (mu!=1.0) if ((fabs(1.0-s/(dtTotal*fabs(vParallel)))>1.0E-2)||(fabs(vParallel_new-vParallel)>1.0E-5)||(fabs(vNorm_new-vNorm)>1.0E-5)) { 
      res=_fail;
    }
  }


  //return the original parameters of the field line
  auto p=VertexData.begin();

  for (auto Vertex=FL::FieldLinesAll[0].GetFirstVertex();Vertex!=NULL;p++,Vertex=Vertex->GetNext()) {
    Vertex->SetDatum(FL::DatumAtVertexPlasmaDensity,p->DensityCurrent);
    Vertex->SetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,p->DensityOld);
    Vertex->SetPlasmaVelocity(p->v);
  }

  SEP::Diffusion::GetPitchAngleDiffusionCoefficient=DiffusionCoeffcient;
}










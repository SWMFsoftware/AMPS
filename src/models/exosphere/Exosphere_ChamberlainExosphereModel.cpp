//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================================================
//$Id$
//===================================================================================
//Chamberlain Exospheric model
//the model is taken from Shen-1963-JAS

#include "pic.h"
#include "Exosphere.h"

double *Exosphere::ChamberlainExosphere::SpecieExobaseEscapeRate=NULL,*Exosphere::ChamberlainExosphere::SpecieExobaseTemperature=NULL;
bool Exosphere::ChamberlainExosphere::ModelInitFlag=false;

//init the Chemberlain Exosphere model
void Exosphere::ChamberlainExosphere::Init(double *ExobaseEscapeRate,double *ExobaseTemperature) {
  SpecieExobaseEscapeRate=new double[PIC::nTotalSpecies];
  SpecieExobaseTemperature=new double[PIC::nTotalSpecies];

  ModelInitFlag=true;
  memcpy(SpecieExobaseEscapeRate,ExobaseEscapeRate,PIC::nTotalSpecies*sizeof(double));
  memcpy(SpecieExobaseTemperature,ExobaseTemperature,PIC::nTotalSpecies*sizeof(double));
}

void Exosphere::ChamberlainExosphere::PrintVariableList(FILE* fout) {
  fprintf(fout,", \"Chamberlain Exosphere Number Density [m^{-3}]\"");
}

void Exosphere::ChamberlainExosphere::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
//Shen-1963-JAS

//  exit(__LINE__,__FILE__,"Error: the coefficients are not consistent with Shen-1963-JAS and Arno's thetis. Chack other implementation of the model");

//  if (PIC::ThisThread==0) {
    double R,H0,g0,E;
    double nT,nS,H,nH,nC,nA;
    double SurfaceNumberDensity;
    double res=0.0;

    if (SpecieExobaseEscapeRate[DataSetNumber]>0.0) {
      //the parameters on the exobase level
      g0=GravityConstant*_MASS_(_TARGET_)/pow(_RADIUS_(_TARGET_),2);
      H0=Kbol*SpecieExobaseTemperature[DataSetNumber]/(PIC::MolecularData::GetMass(DataSetNumber)*g0);
      E=_RADIUS_(_TARGET_)/H0;


      //get the surface density
      SurfaceNumberDensity=SpecieExobaseEscapeRate[DataSetNumber]/Pi/pow(_RADIUS_(_TARGET_),2)/sqrt(8.0*Kbol*SpecieExobaseTemperature[DataSetNumber]/Pi/PIC::MolecularData::GetMass(DataSetNumber));

      //calculate local density
      double x[3];

      CenterNode->GetX(x);
      R=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

      //H=1.0+2.0/E;
      H=3.0/(2.0*E)+1.0/(1.0+sqrt(Pi/E)*exp(E)/_RADIUS_(_TARGET_)*R*(1.0-erf(sqrt(E))));

      if (R>=_RADIUS_(_TARGET_)) {
        nH=exp(-(1.0-_RADIUS_(_TARGET_)/R)*E);
        nC=sqrt(1.0-pow(_RADIUS_(_TARGET_)/R,2))*exp(-E/(1.0+_RADIUS_(_TARGET_)/R));
        nT=nH-nC;

        nS=pow(_RADIUS_(_TARGET_)/R,2)/(2.0*sqrt(Pi*E))*(1.0+E)*exp(-E)/(H-1.0+sqrt(max(1.0-2.0/3.0*H*_RADIUS_(_TARGET_)/R,0.0)));
        nA=nT-nS;
      }
      else nA=0.0;

      res=SurfaceNumberDensity*nA;

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
        if (!isfinite(res)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
  #endif
    }


    fprintf(fout,"  %e  ",res);
//  }

}

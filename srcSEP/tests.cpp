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

#include "tests.h"


void DiffusionCoefficient_const(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment) {
   D=1.0,dD_dmu=0.0;
}


void TestManager() {
  //===============================================
  //Test Dxx 
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
  


























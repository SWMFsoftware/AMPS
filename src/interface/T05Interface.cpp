//$Id$
//Interface to the T05 model

/*
 * T05Interface.cpp
 *
 *  Created: June 2021
 * 	Developed by: Anthony Knighton
 */

#include "pic.h"
#include "constants.h"
#include "constants.PlanetaryData.h"
#include "T05Interface.h"

/* 05::PARMOD:
 * c--------------------------------------------------------------------
C DATA-BASED MODEL CALIBRATED BY (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
C	DISTURBANCE STORM TIME DST (NANOTESLA), (3) Y-COMPONENT OF
C	INTERNAL MAGNETIC FIELD (BYIMF), (4) Z-COMPONENT OF INTERNAL
C       MAGNETIC FIELD (BZIMF), (5)-(10) TIME INTEGRALS FROM
C       BEGINNING OF STORM: W1-W6. THESE ELEMENTS MAKE UP THE 10 ELEMENTS 
C       OF ARRAY PARMOD(10).
 *
 */

double T05::PS=0.170481;
double T05::PARMOD[11]={2.0,-50.0,1.0,-3.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

extern "C"{
  void t05_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

void T05::GetMagneticField(double *B,double *x) {
  int idim,IOPT=0;
  //Calculate global magnetic field in the magnetosphere
  double xLocal[3],bT05[3];

  for (idim=0;idim<3;idim++) xLocal[idim]=x[idim]/_EARTH__RADIUS_;

/* if (pow(xLocal[0],2)+pow(xLocal[1],2)+pow(xLocal[2],2)<1.0) {
   for (idim=0;idim<3;idim++) B[idim]=0.0;
   return;
}*/

  t05_(&IOPT,PARMOD,&PS,xLocal+0,xLocal+1,xLocal+2,bT05+0,bT05+1,bT05+2);
  
  // Calculate Earth's internal magentic field
  IGRF::GetMagneticField(B,x);

  // Sum contributions from Earth's internal magnetic field and global magnetic field
  for (idim=0;idim<3;idim++) B[idim]+=bT05[idim]*_NANO_;

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
  PIC::Debugger::CatchOutLimitValue(xLocal,DIM,__LINE__,__FILE__);
  PIC::Debugger::CatchOutLimitValue(bT05,DIM,__LINE__,__FILE__);
  PIC::Debugger::CatchOutLimitValue(B,DIM,__LINE__,__FILE__);
#endif
#endif
}
  
void T05::SetSolarWindPressure(double PDYN) {PARMOD[0]=PDYN/_NANO_;}
void T05::SetDST(double DST) {PARMOD[1]=DST/_NANO_;}
void T05::SetBYIMF(double BYIMF) {PARMOD[2]=BYIMF/_NANO_;}
void T05::SetBZIMF(double BZIMF) {PARMOD[3]=BZIMF/_NANO_;}

void T05::SetW1(double W1) {PARMOD[4]=W1/_NANO_;}
void T05::SetW2(double W2) {PARMOD[5]=W2/_NANO_;}
void T05::SetW3(double W3) {PARMOD[6]=W3/_NANO_;}
void T05::SetW4(double W4) {PARMOD[7]=W4/_NANO_;}
void T05::SetW5(double W5) {PARMOD[8]=W5/_NANO_;}
void T05::SetW6(double W6) {PARMOD[9]=W6/_NANO_;}

void T05::SetW(double W1,double W2,double W3,double W4,double W5,double W6) {
  SetW1(W1);
  SetW2(W2);
  SetW3(W3);
  SetW4(W4);
  SetW5(W5);
  SetW6(W6);
}


//$Id$
//Interface to Tsyganenko (2016) RBF model (TA16 / RBF_MODEL_2016)

/*
 * TA16Interface.cpp
 *
 * Added to AMPS: 2026-02-26
 *
 * Key practical detail:
 *   TA16 needs the coefficient file TA16_RBF.par (or equivalent).
 *   Without it the model will fail on first call (file open/read).
 *
 * AMPS extension (in TA16.for):
 *   - we introduced TA16_SET_COEFF_FILE(...) which updates the filename in a COMMON block
 *   - default remains TA16_RBF.par for backward compatibility
 */

#include "pic.h"
#include "constants.h"
#include "constants.PlanetaryData.h"
#include "TA16Interface.h"

// Default driver values

double TA16::PS=0.170481;

// PARMOD(10) as documented in TA16_RBF.f
// Only the first 4 are used by the official code; the rest are ignored.
double TA16::PARMOD[10]={
  2.0,  // PDYN [nPa]
  0.0,  // SymHc [nT]
  0.0,  // XIND [dimensionless]
  0.0,  // BYIMF [nT]
  0.0,0.0,0.0,0.0,0.0,0.0
};

std::string TA16::UserFrameName="";

double TA16::UserFrame2GSM[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double TA16::GSM2UserFrame[3][3]={{1,0,0},{0,1,0},{0,0,1}};
bool TA16::Rotate2GSM=false;

double TA16::UserFrame2GSE[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double TA16::GSE2UserFrame[3][3]={{1,0,0},{0,1,0},{0,0,1}};
bool TA16::Rotate2GSE=false;

extern "C"{
  void rbf_model_2016_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
  void ta16_set_coeff_file_(char*,int);
}

void TA16::SetSolarWindPressure(double PDYN) {PARMOD[0]=PDYN/_NANO_;}
void TA16::SetSymHc(double SymHc) {PARMOD[1]=SymHc/_NANO_;}
void TA16::SetXIND(double XIND) {PARMOD[2]=XIND;}
void TA16::SetBYIMF(double BY) {PARMOD[3]=BY/_NANO_;}

void TA16::SetCoeffFileName(const std::string &fname) {
  // Fortran expects a fixed-length CHARACTER argument; the common pattern is
  // to pass (char*, int len) from C/C++.
  //
  // We do *not* pad here; the Fortran side copies the bytes and will include
  // trailing spaces if present. So pass the raw string.
  ta16_set_coeff_file_(const_cast<char*>(fname.c_str()), (int)fname.size());
}

void TA16::GetMagneticField(double *B,double *x) {
  int idim;
  int IOPT=0;

  double bTA16[3];

  double xLocal[3];
  for (idim=0;idim<3;idim++) xLocal[idim]=x[idim]/_EARTH__RADIUS_;

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__TA16_

  if (Rotate2GSM==false) {
    rbf_model_2016_(&IOPT,PARMOD,&PS,xLocal+0,xLocal+1,xLocal+2,bTA16+0,bTA16+1,bTA16+2);
  }
  else {
    double xLocal_GSM[3],bTA16_GSM[3];
    mxv_c(UserFrame2GSM,xLocal,xLocal_GSM);
    rbf_model_2016_(&IOPT,PARMOD,&PS,xLocal_GSM+0,xLocal_GSM+1,xLocal_GSM+2,bTA16_GSM+0,bTA16_GSM+1,bTA16_GSM+2);
    mxv_c(GSM2UserFrame,bTA16_GSM,bTA16);
  }

  #else
  exit(__LINE__,__FILE__,"Error: TA16 is not setup for this run. Use option 'CouplerMode=TA16' in the input file to use TA16");
  #endif

  IGRF::GetMagneticField(B,x);
  for (idim=0;idim<3;idim++) B[idim]+=bTA16[idim]*_NANO_;

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  #if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
  PIC::Debugger::CatchOutLimitValue(xLocal,DIM,__LINE__,__FILE__);
  PIC::Debugger::CatchOutLimitValue(bTA16,DIM,__LINE__,__FILE__);
  PIC::Debugger::CatchOutLimitValue(B,DIM,__LINE__,__FILE__);
  #endif
  #endif
}

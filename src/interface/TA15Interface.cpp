//$Id$
//Interface to the Tsyganenko & Andreeva (2015) forecasting model (TA15)

/*
 * TA15Interface.cpp
 *
 * Added to AMPS: 2026-02-26
 *
 * This file is intentionally verbose (heavy comments) because TA15 is often
 * the first "modern" Tsyganenko-family model users integrate.
 *
 * The interface mirrors T96Interface.cpp:
 *   - inputs are in SI and user frame
 *   - position converted to Re for the Fortran call
 *   - optional frame rotation into GSM
 *   - adds IGRF internal field
 *   - converts external nT -> Tesla
 */

#include <cmath>
#include "pic.h"
#include "constants.h"
#include "constants.PlanetaryData.h"
#include "TA15Interface.h"

// Default driver values (quiet-ish)
TA15::cVersion TA15::Version = TA15::Version_B;

double TA15::PS=0.170481;

// PARMOD(10) as documented in TA_2015_*.f
// Only the first 4 are used by the official code; the rest are ignored.
double TA15::PARMOD[10]={
  2.0,  // PDYN [nPa]
  0.0,  // BYIMF [nT]
  0.0,  // BZIMF [nT]
  0.0,  // XIND  [dimensionless]
  0.0,0.0,0.0,0.0,0.0,0.0
};

std::string TA15::UserFrameName="";

double TA15::UserFrame2GSM[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double TA15::GSM2UserFrame[3][3]={{1,0,0},{0,1,0},{0,0,1}};
bool TA15::Rotate2GSM=false;

double TA15::UserFrame2GSE[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double TA15::GSE2UserFrame[3][3]={{1,0,0},{0,1,0},{0,0,1}};
bool TA15::Rotate2GSE=false;

extern "C"{
  void ta_2015_b_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
  void ta_2015_n_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

void TA15::SetSolarWindPressure(double PDYN) {PARMOD[0]=PDYN/_NANO_;}
void TA15::SetBYIMF(double BY) {PARMOD[1]=BY/_NANO_;}
void TA15::SetBZIMF(double BZ) {PARMOD[2]=BZ/_NANO_;}
void TA15::SetXIND(double XIND) {PARMOD[3]=XIND;}

void TA15::SetVersion(cVersion v) {Version=v;}

static inline void CallTA15_Fortran(int *IOPT,double *PARMOD,double *PS,
                                   double *x,double *y,double *z,
                                   double *bx,double *by,double *bz,
                                   TA15::cVersion Version) {
  // Helper so the call site stays readable.
  if (Version==TA15::Version_B) {
    ta_2015_b_(IOPT,PARMOD,PS,x,y,z,bx,by,bz);
  }
  else {
    ta_2015_n_(IOPT,PARMOD,PS,x,y,z,bx,by,bz);
  }
}

// TA15 describes the field of the magnetospheric current systems INSIDE the
// magnetopause only: its shielding expansions are fitted on the model magnetopause
// and are not bounded outside it.  Upstream of the dayside magnetopause the returned
// "external field" grows without limit (for PDYN=3.7 nPa, BZIMF=-4 nT, XIND=0.18:
// |B| = 134 nT at X=16 Re, 1.5e3 nT at 18 Re, 4.6e9 nT at 20 Re on the Sun-Earth
// line), so backtraced cutoff trajectories that leave sunward through a +X domain
// face got trapped in this unphysical field and were classified forbidden even at
// 20 GV (vertical Rc=-1 at the equator).  T96 and T01 confine their fields to a model
// magnetopause and stay bounded there (3-30 nT at X=8..29 Re for the same event); for
// TA15 the external part is therefore dropped outside the
// Shue et al. (1998) magnetopause for the same PDYN and BZIMF (the internal IGRF field
// is kept everywhere).  The surface depends only on PARMOD, so the test is reentrant.
static inline bool TA15_InsideShue1998Magnetopause(const double xGSM_Re[3],
                                                   const double *PARMOD) {
  const double pdyn=PARMOD[0],bz=PARMOD[2];
  if (!(pdyn>0.0) || !std::isfinite(pdyn) || !std::isfinite(bz)) return true;

  const double r=std::sqrt(xGSM_Re[0]*xGSM_Re[0]+xGSM_Re[1]*xGSM_Re[1]+xGSM_Re[2]*xGSM_Re[2]);
  if (r<=0.0) return true;

  const double r0=(10.22+1.29*std::tanh(0.184*(bz+8.14)))*std::pow(pdyn,-1.0/6.6);
  const double alpha=(0.58-0.007*bz)*(1.0+0.024*std::log(pdyn));
  const double denominator=1.0+xGSM_Re[0]/r;
  if (denominator<=1.0E-12) return true;   // anti-sunward axis: the surface is open there

  return r<r0*std::pow(2.0/denominator,alpha);
}

void TA15::GetMagneticField(double *B,double *x) {
  int idim;
  int IOPT=0;

  double bTA15[3];

  // Convert SI position to Earth radii for Tsyganenko.
  double xLocal[3];
  for (idim=0;idim<3;idim++) xLocal[idim]=x[idim]/_EARTH__RADIUS_;

  if (Rotate2GSM==false) {
    if (TA15_InsideShue1998Magnetopause(xLocal,PARMOD)) {
      CallTA15_Fortran(&IOPT,PARMOD,&PS,xLocal+0,xLocal+1,xLocal+2,bTA15+0,bTA15+1,bTA15+2,Version);
    }
    else {
      for (idim=0;idim<3;idim++) bTA15[idim]=0.0;
    }
  }
  else {
    double xLocal_GSM[3],bTA15_GSM[3];
    Geopack::MatrixVectorMultiply(UserFrame2GSM,xLocal,xLocal_GSM);
    if (TA15_InsideShue1998Magnetopause(xLocal_GSM,PARMOD)) {
      CallTA15_Fortran(&IOPT,PARMOD,&PS,xLocal_GSM+0,xLocal_GSM+1,xLocal_GSM+2,bTA15_GSM+0,bTA15_GSM+1,bTA15_GSM+2,Version);
    }
    else {
      for (idim=0;idim<3;idim++) bTA15_GSM[idim]=0.0;
    }
    Geopack::MatrixVectorMultiply(GSM2UserFrame,bTA15_GSM,bTA15);
  }

  // Internal field
  IGRF::GetMagneticField(B,x);

  // External in nT -> Tesla
  for (idim=0;idim<3;idim++) B[idim]+=bTA15[idim]*_NANO_;

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  #if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
  PIC::Debugger::CatchOutLimitValue(xLocal,DIM,__LINE__,__FILE__);
  PIC::Debugger::CatchOutLimitValue(bTA15,DIM,__LINE__,__FILE__);
  PIC::Debugger::CatchOutLimitValue(B,DIM,__LINE__,__FILE__);
  #endif
  #endif
}

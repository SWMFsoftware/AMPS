//$Id$
//Interface to Tsyganenko (2016) RBF model (TA16 / RBF_MODEL_2016)

/*
 * TA16Interface.h
 *
 * Added to AMPS: 2026-02-26
 *
 * The official TA16 distribution (RBF_MODEL_2016) requires an external coefficient
 * file (TA16_RBF.par) which contains the linear coefficients for the RBF expansions.
 *
 * AMPS modification:
 *   We patched the Fortran source (TA16.for) to allow overriding the coefficient
 *   file location via a helper subroutine TA16_SET_COEFF_FILE(...).
 *
 * This wrapper:
 *   - exposes TA16::SetCoeffFileName(...)
 *   - follows the same frame-rotation and internal-field addition patterns as T96/T01/TA15
 */

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <utility>
#include <map>
#include <cmath>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <string>

#include "SpiceUsr.h"

#ifndef _INTERFACE_TA16INTERFACE_H_
#define _INTERFACE_TA16INTERFACE_H_

#include "GeopackInterface.h"

namespace TA16 {
  using namespace Geopack;

  extern std::string UserFrameName;

  extern double UserFrame2GSM[3][3],GSM2UserFrame[3][3];
  extern bool Rotate2GSM;

  extern double UserFrame2GSE[3][3],GSE2UserFrame[3][3];
  extern bool Rotate2GSE;

  inline void Init(const char* Epoch, std::string FrameNameIn) {
    Geopack::Init(Epoch,FrameNameIn);

    UserFrameName=FrameNameIn;

    if (UserFrameName!="GSM") {
      double et;
      Rotate2GSM=true;
      utc2et_c(Epoch,&et);
      pxform_c(UserFrameName.c_str(),"GSM",et,UserFrame2GSM);
      pxform_c("GSM",UserFrameName.c_str(),et,GSM2UserFrame);
    }

    if (UserFrameName!="GSE") {
      double et;
      Rotate2GSE=true;
      utc2et_c(Epoch,&et);
      pxform_c(UserFrameName.c_str(),"GSE",et,UserFrame2GSE);
      pxform_c("GSE",UserFrameName.c_str(),et,GSE2UserFrame);
    }
  }

  extern double PS;

  // TA16 uses PARMOD(10) (first 4 required):
  //   (1) PDYN   [nPa]
  //   (2) SymHc  [nT]  (corrected SymH index)
  //   (3) XIND   [dimensionless]
  //   (4) BYIMF  [nT]
  extern double PARMOD[10];

  void SetSolarWindPressure(double PDYN);
  void SetSymHc(double SymHc);
  void SetXIND(double XIND);
  void SetBYIMF(double BY);

  // TA16 coefficient file management (AMPS extension).
  void SetCoeffFileName(const std::string &fname);

  void GetMagneticField(double *B,double *x);
}

#endif /* _INTERFACE_TA16INTERFACE_H_ */

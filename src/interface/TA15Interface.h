//$Id$
//Interface to the Tsyganenko & Andreeva (2015) forecasting model (TA15)

/*
 * TA15Interface.h
 *
 * Added to AMPS: 2026-02-26
 *
 * TA15 has two official variants in the distribution:
 *   - TA_2015_B : "B" fit (entire dataset) version
 *   - TA_2015_N : "N" (alternative fit) version
 *
 * Both share the same calling convention:
 *   SUBROUTINE TA_2015_*(IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
 * with PARMOD(10) where the first 4 elements are mandatory:
 *   PARMOD(1) PDYN  [nPa]
 *   PARMOD(2) BYIMF [nT]
 *   PARMOD(3) BZIMF [nT]
 *   PARMOD(4) XIND  [dimensionless coupling index, typical 0..2]
 *
 * This wrapper follows the same pattern as T96Interface/T05Interface:
 *   - optional SPICE rotations for user frames
 *   - adds internal field via IGRF
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

#ifndef _INTERFACE_TA15INTERFACE_H_
#define _INTERFACE_TA15INTERFACE_H_

#include "GeopackInterface.h"

namespace TA15 {
  using namespace Geopack;

  // Which TA15 variant to call.
  enum cVersion {
    Version_B = 0,
    Version_N = 1
  };

  extern cVersion Version;

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

  // dipole tilt
  extern double PS;

  // TA15 uses PARMOD(10)
  extern double PARMOD[10];

  // Common driver setters
  void SetSolarWindPressure(double PDYN);
  void SetBYIMF(double BY);
  void SetBZIMF(double BZ);
  void SetXIND(double XIND);

  // Select variant (B or N)
  void SetVersion(cVersion v);

  // Total field in SI (internal + external)
  void GetMagneticField(double *B,double *x);
}

#endif /* _INTERFACE_TA15INTERFACE_H_ */

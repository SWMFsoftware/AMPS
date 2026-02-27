//$Id$
//Interface to the Tsyganenko T01 model (T01_01)

/*
 * T01Interface.h
 *
 *  Added to AMPS: 2026-02-26
 *
 *  This header mirrors the existing T96Interface.h and T05Interface.h patterns.
 *  The goal is that AMPS can treat T01 just like other analytic magnetospheric
 *  background field models:
 *
 *    - a light C++ wrapper around the official Fortran routine
 *    - optional SPICE frame rotation (User frame -> GSM/GSE)
 *    - (optionally) add Earth's internal field via Geopack/IGRF, consistent with T96
 *
 *  NOTE ON NOMENCLATURE:
 *    "T01" typically refers to the Tsyganenko (2001/2002) near-magnetosphere model.
 *    The official field routine in the distribution is named T01_01.
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

#ifndef _INTERFACE_T01INTERFACE_H_
#define _INTERFACE_T01INTERFACE_H_

#include "GeopackInterface.h"

namespace T01 {
  using namespace Geopack;

  // Frame handling follows the exact pattern used in T96Interface.h:
  //  - User supplies an Epoch string (UTC) and a coordinate frame name
  //  - We build SPICE rotation matrices to/from GSM and GSE
  extern std::string UserFrameName;

  extern double UserFrame2GSM[3][3],GSM2UserFrame[3][3];
  extern bool Rotate2GSM;

  extern double UserFrame2GSE[3][3],GSE2UserFrame[3][3];
  extern bool Rotate2GSE;

  inline void Init(const char* Epoch, std::string FrameNameIn) {
    // Initializes Geopack (dipole tilt, etc.) and caches frame transforms.
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

  // Dipole tilt angle (radians) used by Tsyganenko models.
  extern double PS;

  // T01 uses a PARMOD array. In the official code it is documented as PARMOD(10)
  // (the first entries are PDYN/DST/BY/BZ etc., depending on the exact T01 variant).
  //
  // For consistency with other interfaces, we keep a fixed-size array.
  extern double PARMOD[11];

  // Minimal setters for common solar-wind/geomagnetic drivers.
  // Units follow the conventions in the Fortran implementation (nPa, nT).
  void SetSolarWindPressure(double SolarWindPressure);
  void SetDST(double DST);
  void SetBYIMF(double BYIMF);
  void SetBZIMF(double BZIMF);

  // Main entry: returns TOTAL B (internal IGRF + external T01) in SI units.
  // Input position x is in SI (meters) in the user frame specified in Init().
  void GetMagneticField(double *B,double *x);
}

#endif /* _INTERFACE_T01INTERFACE_H_ */

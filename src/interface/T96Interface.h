//$Id$
//interface to T96 model

/*
 * T96Interface.h
 *
 *  Created on: Jan 5, 2017
 *      Author: vtenishe
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


#ifndef _INTERFACE_T96INTERFACE_H_
#define _INTERFACE_T96INTERFACE_H_

#include "GeopackInterface.h"

namespace T96 {
  using namespace Geopack;

  extern std::string UserFrameName;

  extern double UserFrame2GSM[3][3],GSM2UserFrame[3][3];
  extern bool Rotate2GSM;

  extern double UserFrame2GSE[3][3],GSE2UserFrame[3][3];
  extern bool Rotate2GSE;

  void inline Init(const char* Epoch, std::string FrameNameIn) {
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


  //dipole tilt angle
  extern double PS;

  //the set of themodel paramters
  extern double PARMOD[11];
  void SetSolarWindPressure(double SolarWindPressure);
  void SetDST(double DST);
  void SetBYIMF(double BYIMF);
  void SetBZIMF(double BZIMF);


  void GetMagneticField(double *B,double *x);


}



#endif /* SRC_INTERFACE_T96INTERFACE_H_ */

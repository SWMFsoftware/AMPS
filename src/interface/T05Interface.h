//$Id$
//Interface to the T05 model

/*
 * T05Interface.h
 *
 *  Created: June 2021
 * 	Developed by: Anthony Knighton
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

#ifndef _INTERFACE_T05INTERFACE_H_
#define _INTERFACE_T05INTERFACE_H_

#include "GeopackInterface.h"


namespace T05 {
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

  //set of model parameters
  //PARMOD=array containing all parameters
  extern double PARMOD[11];

  //PDYN=solar wind pressure
  void SetSolarWindPressure(double SolarWindPressure); 
  //DST=disturbance storm time
  void SetDST(double DST);
  //BYIMF=y-component of IMF
  void SetBYIMF(double BYIMF);
  //BZIMF=z-component of IMF
  void SetBZIMF(double BZIMF);
  //W1-W6=time integrals from start of storm. These are defined in detiai in Tsyganenko 2005.
  void SetW1(double W1);
  void SetW2(double W2);
  void SetW3(double W3);
  void SetW4(double W4);
  void SetW5(double W5);
  void SetW6(double W6);
  void SetW(double W1,double W2,double W3,double W4,double W5,double W6); 



  void GetMagneticField(double *B,double *x);


}



#endif /* SRC_INTERFACE_T05INTERFACE_H_ */


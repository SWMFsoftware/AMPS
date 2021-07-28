//Interface to the Geopack package
//$Id$

/*
 * GeopackInterface.h
 *
 *  Created on: Jan 4, 2017
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



#ifndef _SRC_INTERFACE_GEOPACKINTERFACE_H_
#define _SRC_INTERFACE_GEOPACKINTERFACE_H_

namespace Geopack {

  extern std::string UserFrameName;

  extern double UserFrame2GSM[3][3],GSM2UserFrame[3][3];
  extern bool Rotate2GSM;

  extern double UserFrame2GSE[3][3],GSE2UserFrame[3][3];
  extern bool Rotate2GSE;


  //IGRF mgnetoc field model
  namespace IGRF {
    void GetMagneticField(double *B,double *x);
  }

  void Init(const char* Epoch,std::string FrameNameIn);
}



#endif /* SRC_INTERFACE_GEOPACKINTERFACE_H_ */

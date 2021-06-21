//$Id$
//Interface to the T05 model

/*
 * T05Interface.h
 *
 *  Created: June 2021
 * 	Developed by: Anthony Knighton
 */

#ifndef _INTERFACE_T05INTERFACE_H_
#define _INTERFACE_T05INTERFACE_H_

#include "GeopackInterface.h"

namespace T05 {
  using namespace Geopack;

  //dipole tilt angle
  extern double PS;

  //set of model parameters
  //PARMOD=array containing all parameters
  extern double PARMOD[11];
  //PDYN=solar wind pressure
  void SetPDYN(double PDYN);
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


  void GetMagneticField(double *B,double *x);


}



#endif /* SRC_INTERFACE_T05INTERFACE_H_ */


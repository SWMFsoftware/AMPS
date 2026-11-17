//$Id$

/*
 * Surface_temperature.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: vtenishe
 */

#include "Surface.h"


//surface temeprature in the model of the isothremal temeprature; 0 means that the wall temperature is not set: the
//models that use it stop the run (Surface::GetSurfaceTemperature)
double Surface::Temperature::Isothremal::Temp=0.0;




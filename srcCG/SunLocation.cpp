//$Id$
//contains variables and functions used to update the location of the Sun during a model run

/*
 * SunLocation.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "Dust.h"
#include "Comet.h"
#include "RosinaMeasurements.h"

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif


char Comet::SunLocationUpdate::StartTime[_MAX_STRING_LENGTH_PIC_]="";

double Comet::SunLocationUpdate::StartTime_et=NAN;
double Comet::SunLocationUpdate::TimeIncrement=0.0;

int Comet::SunLocationUpdate::OutputCycleCounter=0;
int Comet::SunLocationUpdate::OutputCycleStep=0;
int Comet::SunLocationUpdate::FirstUpdateOutputCycleNumber=0;

extern double HeliocentricDistance,subSolarPointAzimuth,subSolarPointZenith;

void Comet::SunLocationUpdate::Processor() {

  #ifndef _NO_SPICE_CALLS_
  //increment the call counter; proceed with the location only the reqired numer of the prior outputs is reached
  if (FirstUpdateOutputCycleNumber>++OutputCycleCounter) return;

  //if StartTime_et is not initialized yet -> initialize it
  if (std::isnan(StartTime_et)==true) {
    utc2et_c(StartTime,&StartTime_et);
  }

  //update the location of the Sun
  if ((OutputCycleCounter-FirstUpdateOutputCycleNumber)%OutputCycleStep==0) {
    SpiceDouble lt,xSun[3];

    //update the location of the Sun and the nucleus shadowing
    spkpos_c("SUN",StartTime_et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",xSun,&lt);
    reclat_c(xSun,&HeliocentricDistance,&subSolarPointAzimuth,&subSolarPointZenith);
    subSolarPointZenith = Pi/2 - subSolarPointZenith;
    HeliocentricDistance*=1.0E3;

    double xLightSource[3]={HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith),
        HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith),
        HeliocentricDistance*cos(subSolarPointZenith)};

    PIC::RayTracing::SetCutCellShadowAttribute(xLightSource,false);

    //update the output time stamp usen in the file names, and print the triangulation with the shadowing
    char fname[100];

    et2utc_c (StartTime_et,"ISOC",4,100,RosinaSample::SamplingTimeInterval::CurrentTimeWindowStamp);
    sprintf(fname,"SurfaceTriangulation-shadow.time=%s.dat",RosinaSample::SamplingTimeInterval::CurrentTimeWindowStamp);

    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname,PIC::OutputDataFileDirectory);

    //reset the sampling point flags for sampling Rosina measurements
    if (_COMET_SAMPLE_ROSINA_DATA_==_PIC_MODE_ON_) {
      RosinaSample::Init(StartTime_et,StartTime_et+TimeIncrement);
    }

    //increment the time counter
    StartTime_et+=TimeIncrement;
  }
  #endif //_NO_SPICE_CALLS_
}



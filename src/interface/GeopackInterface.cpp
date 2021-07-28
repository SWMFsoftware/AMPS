//interface to Geopack 2008
//$Id$

#include <stdio.h>
#include <time.h>
#include <strings.h>
#include <math.h>



#include "GeopackInterface.h"

#include "pic.h"
#include "constants.h"
#include "constants.PlanetaryData.h"
#include "ifileopr.h"
#include "specfunc.h"


/*
C  IMPORTANT: IF ONLY QUESTIONABLE INFORMATION (OR NO INFORMATION AT ALL) IS AVAILABLE
C             ON THE SOLAR WIND SPEED, OR, IF THE STANDARD GSM AND/OR SM COORDINATES ARE
C             INTENDED TO BE USED, THEN SET VGSEX=-400.0 AND VGSEY=VGSEZ=0. IN THIS CASE,
C             THE GSW COORDINATE SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM.
*/

std::string Geopack::UserFrameName;
double Geopack::UserFrame2GSM[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double Geopack::GSM2UserFrame[3][3]={{1,0,0},{0,1,0},{0,0,1}};
bool Geopack::Rotate2GSM=false;

double Geopack::UserFrame2GSE[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double Geopack::GSE2UserFrame[3][3]={{1,0,0},{0,1,0},{0,0,1}};
bool Geopack::Rotate2GSE=false;

extern "C"{
  void recalc_08_(int*,int*,int*,int*,int*,double*,double*,double*);
  void sphcar_08_(double*,double*,double*,double*,double*,double*,int*);
  void bspcar_08_(double*,double*,double*,double*,double*,double*,double*,double*);
  void igrf_geo_08_(double*,double*,double*,double*,double*,double*);
  void igrf_gsw_08_(double*,double*,double*,double*,double*,double*);

  void gswgse_08_(double*,double*,double*,double*,double*,double*,int*);
}


void Geopack::Init(const char* Epoch,std::string FrameNameIn) {
  CiFileOperations Parser;

  UserFrameName=FrameNameIn;

  if (UserFrameName!="GSE") {
    double et;

    Rotate2GSE=true;

    utc2et_c(Epoch,&et);

    pxform_c(UserFrameName.c_str(),"GSE",et,UserFrame2GSE);
    pxform_c("GSE",UserFrameName.c_str(),et,GSE2UserFrame);
  }

  if (UserFrameName!="GSM") {
    double et;

    Rotate2GSM=true;

    utc2et_c(Epoch,&et);

    pxform_c(UserFrameName.c_str(),"GSM",et,UserFrame2GSM);
    pxform_c("GSM",UserFrameName.c_str(),et,GSM2UserFrame);
  }



  //conver the epoch string to the Geopack format
  int cnt=0;
  char tEpoch[300],str[100],*endptr;
  int Year,Month,Day,Hour,Minute,Second,DayOfYear;

  do {
    tEpoch[cnt]=Epoch[cnt];

    if ((tEpoch[cnt]=='T')||(tEpoch[cnt]=='-')||(tEpoch[cnt]==':')) tEpoch[cnt]=' ';
    cnt++;
  }
  while ((Epoch[cnt]!='\n')&&(Epoch[cnt]!=0));

  Parser.CutInputStr(str,tEpoch);
  Year=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Month=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Day=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Hour=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Minute=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Second=(int)strtol(str,&endptr,10);

  struct tm calendar = {0};
  calendar.tm_year = Year - 1900;
  calendar.tm_mon  = Month - 1;
  calendar.tm_mday = Day;
  time_t date = mktime ( &calendar );
  DayOfYear=calendar.tm_yday + 1;

  //init Geopack
  double VGSE[3]={-400.0,0.0,0.0};

  recalc_08_(&Year,&DayOfYear,&Hour,&Minute,&Second,VGSE+0,VGSE+1,VGSE+2);
}


void Geopack::IGRF::GetMagneticField(double *B,double *x) {
  /*
   *      SUBROUTINE IGRF_GSW_08 (XGSW,YGSW,ZGSW,HXGSW,HYGSW,HZGSW)
   *
   *      SUBROUTINE GSWGSE_08 (XGSW,YGSW,ZGSW,XGSE,YGSE,ZGSE,J)
   *      C                    J>0                       J<0
   *      C-----INPUT:   J,XGSW,YGSW,ZGSW          J,XGSE,YGSE,ZGSE
   *      C-----OUTPUT:    XGSE,YGSE,ZGSE            XGSW,YGSW,ZGSW
   *
   *
   */

  int idim;
  double xLocal[3],xLocalGSM[3],bGSM[3];

  for (idim=0;idim<3;idim++) xLocal[idim]=x[idim]/_EARTH__RADIUS_;

  //convert xLocal in Usr Frame to xLocalGSE
  if (Rotate2GSM==true) {
    mxv_c(UserFrame2GSM,xLocal,xLocalGSM);
  }
  else {
    memcpy(xLocalGSM,xLocal,3*sizeof(double));
  }
    
  //extract the magnetic field vector
  igrf_gsw_08_(xLocalGSM+0,xLocalGSM+1,xLocalGSM+2,bGSM+0,bGSM+1,bGSM+2);
  
  //cover bGSM to b in the User frame
  if (Rotate2GSM==true) {
    mxv_c(GSM2UserFrame,bGSM,B);
    for (idim=0;idim<3;idim++) B[idim]*=_NANO_;
  }
  else {
    for (idim=0;idim<3;idim++) B[idim]=bGSM[idim]*_NANO_;
  }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
  PIC::Debugger::CatchOutLimitValue(B,DIM,__LINE__,__FILE__);
#endif
#endif
}























//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury_OrbitalMotion.cpp
 *
 *  Created on: Mar 24, 2012
 *      Author: vtenishe
 */

//Orbital motion of Mercury
//$Id$


#include "pic.h"
#include "Exosphere.h"

double Exosphere::OrbitalMotion::AccumulatedPlanetRotation=0.0,Exosphere::OrbitalMotion::TotalSimulationTime=0.0,Exosphere::OrbitalMotion::TAA=0.0;
double Exosphere::OrbitalMotion::CoordinateFrameRotationRate=0.0;
double Exosphere::OrbitalMotion::RotationRateVector_SO_J2000[3];

//SPICE ephemeris time
SpiceDouble Exosphere::OrbitalMotion::et=0.0,Exosphere::OrbitalMotion::lt=0.0;

//direction to the Sun and the angle of the rotaio berween planetary axises and the direction to the Sun on the Z-plane
double Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[3]={0.0,0.0,0.0};

//matrixes for tranformation SO->IAU and IAU->SO coordinate frames
SpiceDouble Exosphere::OrbitalMotion::SO_to_IAU_TransformationMartix[6][6]={
    {1.0,0.0,0.0, 0.0,0.0,0.0},{0.0,1.0,0.0, 0.0,0.0,0.0},{0.0,0.0,1.0, 0.0,0.0,0.0},
    {0.0,0.0,0.0, 1.0,0.0,0.0},{0.0,0.0,0.0, 0.0,1.0,0.0},{0.0,0.0,0.0, 0.0,0.0,1.0},
};

SpiceDouble Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[6][6]={
    {1.0,0.0,0.0, 0.0,0.0,0.0},{0.0,1.0,0.0, 0.0,0.0,0.0},{0.0,0.0,1.0, 0.0,0.0,0.0},
    {0.0,0.0,0.0, 1.0,0.0,0.0},{0.0,0.0,0.0, 0.0,1.0,0.0},{0.0,0.0,0.0, 0.0,0.0,1.0},
};

//the number intervals of orbit points printed for the time interval between two outputs of the data file
int Exosphere::OrbitalMotion::nOrbitalPositionOutputMultiplier=1;

//update the transformation matrixes
void Exosphere::OrbitalMotion::UpdateTransformationMartix() {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  sxform_c(Exosphere::SO_FRAME,Exosphere::IAU_FRAME,et,SO_to_IAU_TransformationMartix);
  sxform_c(Exosphere::IAU_FRAME,Exosphere::SO_FRAME,et,IAU_to_SO_TransformationMartix);
#else
  int i,j;

  for (i=0;i<6;i++) for (j=0;j<6;j++) {
    SO_to_IAU_TransformationMartix[i][j]=(i==j) ? 1.0 : 0.0;
    IAU_to_SO_TransformationMartix[i][j]=(i==j) ? 1.0 : 0.0;
  }
#endif
}

void Exosphere::OrbitalMotion::UpdateRotationVector_SO_J2000() {
  int idim;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble xform[6][6],rot[3][3],av[3];

  sxform_c ("J2000",Exosphere::SO_FRAME,OrbitalMotion::et,xform);
  xf2rav_c (xform,rot,av);

  for (CoordinateFrameRotationRate=0.0,idim=0;idim<3;idim++) {
    RotationRateVector_SO_J2000[idim]=av[idim];
    CoordinateFrameRotationRate+=pow(av[idim],2);
  }

  CoordinateFrameRotationRate=sqrt(CoordinateFrameRotationRate);
#else
  for (idim=0;idim<3;idim++) RotationRateVector_SO_J2000[idim]=0.0;
  CoordinateFrameRotationRate=0.0;
#endif
}

double Exosphere::OrbitalMotion::GetTAA(const char* TargetName, const char* CenterBodyName, double CenterBodyMass, SpiceDouble EphemerisTime) {
  double res=0.0;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble ltlocal,TargetState[6];
  double EccentricityVector[3];
  double v2,r,a,c,absEccentricity;
  int idim;

  const double GravitationalParameter=1.0E-9*GravityConstant*CenterBodyMass; //the unit is km^2/s^2

  spkezr_c(TargetName,EphemerisTime,"J2000","none",CenterBodyName,TargetState,&ltlocal);

  r=sqrt(TargetState[0]*TargetState[0]+TargetState[1]*TargetState[1]+TargetState[2]*TargetState[2]);
  v2=TargetState[3+0]*TargetState[3+0]+TargetState[3+1]*TargetState[3+1]+TargetState[3+2]*TargetState[3+2];
  c=TargetState[0]*TargetState[3+0]+TargetState[1]*TargetState[3+1]+TargetState[2]*TargetState[3+2];

  for (idim=0,absEccentricity=0.0,a=0.0;idim<3;idim++) {
    EccentricityVector[idim]=v2*TargetState[idim]/GravitationalParameter - c*TargetState[3+idim]/GravitationalParameter - TargetState[idim]/r;
    absEccentricity+=EccentricityVector[idim]*EccentricityVector[idim];
    a+=EccentricityVector[idim]*TargetState[idim];
  }

  absEccentricity=sqrt(absEccentricity);
  res=acos(a/(absEccentricity*r));

  if (c<0.0) res=2.0*Pi-res;
#endif

  return res;

}


double Exosphere::OrbitalMotion::GetTAA(const char* UTC) {
  double res=0.0;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble EphemerisTime;

  utc2et_c(UTC,&EphemerisTime);
  res=GetTAA(EphemerisTime);
#endif

  return res;
}

double Exosphere::OrbitalMotion::GetPhaseAngle(SpiceDouble EphemerisTime) {
  double res=0.0;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble StateSun[6],StateEarth[6],ltlocal;

  spkezr_c("Sun",EphemerisTime,"J2000","none",ObjectName,StateSun,&ltlocal);
  spkezr_c("Earth",EphemerisTime,"J2000","none",ObjectName,StateEarth,&ltlocal);

  double c=0.0,l0=0.0,l1=0.0;
  int idim;

  for (idim=0;idim<3;idim++) {
    l0+=pow(StateSun[idim],2);
    l1+=pow(StateEarth[idim],2);
    c+=StateSun[idim]*StateEarth[idim];
  }

  res=acos(c/sqrt(l0*l1));
#endif

  return res;
}

double Exosphere::OrbitalMotion::GetPhaseAngle(const char* UTC) {
  double res=0.0;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble EphemerisTime;

  utc2et_c(UTC,&EphemerisTime);
  res=GetPhaseAngle(EphemerisTime);
#endif

  return res;
}

//Get Rotations of the frame
//the calcualtion is described in Zhuravlev-VF-Osnovy-teoreticheskoi-mehaniki-2-e-izdanie, chapter 2, page 30
void Exosphere::OrbitalMotion::FrameRotation::GetRotationAxis(double *RotationAxis,double &RotationAngle,const char *FrameName,double etStartRotation,double etFinishRotation) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble T1[6][6],T2[6][6];
  double t,T[3][3];
  int i,j,k;

  sxform_c(FrameName,"J2000",etFinishRotation,T1);
  sxform_c("J2000",FrameName,etStartRotation,T2);

  for (i=0;i<3;i++) for (j=0;j<3;j++) {
    T[i][j]=0.0;

    for (k=0;k<3;k++) T[i][j]+=T2[i][k]*T1[k][j];
  }

  //determine the rate and the vector of the rotation
  RotationAngle=acos((T[0][0]+T[1][1]+T[2][2]-1.0)/2.0);

  t=2.0*sin(RotationAngle);
  RotationAxis[0]=(T[2][1]-T[1][2])/t;
  RotationAxis[1]=(T[0][2]-T[2][0])/t;
  RotationAxis[2]=(T[1][0]-T[0][1])/t;
#else
  RotationAngle=0.0;
  for (int idim=0;idim<3;idim++) RotationAxis[idim]=0.0;
#endif
}


double Exosphere::OrbitalMotion::FrameRotation::GetRotationVector(double *RotationVector,const char *FrameName,double etStartRotation,double etFinishRotation) {
  double Angle;
  int i;

  GetRotationAxis(RotationVector,Angle,FrameName,etStartRotation,etFinishRotation);
  for (i=0;i<3;i++) RotationVector[i]*=Angle/(etFinishRotation-etStartRotation);

  return Angle/(etFinishRotation-etStartRotation);
}


//the calcualtion is described in Korn-GA-Korn-TM-Spravochnik-po-matematike-dlya-nauchnyh-rabotnikov-i-inzhenerov, chapter 14.10 (Matematicheskoe opisanie vrachenii), page 447
void Exosphere::OrbitalMotion::FrameRotation::GetRotationMatrix(double RotationMatrix[3][3],double *RotationAxis,double RotationAngle) {
  double r,l[3];
  int i,j;

  r=sqrt(RotationAxis[0]*RotationAxis[0]+RotationAxis[1]*RotationAxis[1]+RotationAxis[2]*RotationAxis[2]);

  for (i=0;i<3;i++) {
    l[i]=RotationAxis[i]/r;

    RotationMatrix[0][i]=0.0,RotationMatrix[1][i]=0.0,RotationMatrix[2][i]=0.0;
  }

  //compose the transformation matrix
  double cosRotationAngle,sinRotationAngle;

  cosRotationAngle=cos(RotationAngle);
  sinRotationAngle=sin(RotationAngle);

  for (i=0;i<3;i++) RotationMatrix[i][i]+=cosRotationAngle;
  for (i=0;i<3;i++) for (j=0;j<3;j++) RotationMatrix[i][j]+=(1.0-cosRotationAngle)*l[i]*l[j];

  RotationMatrix[0][1]-=sinRotationAngle*l[2],RotationMatrix[0][2]+=sinRotationAngle*l[1];
  RotationMatrix[1][0]+=sinRotationAngle*l[2],RotationMatrix[1][2]-=sinRotationAngle*l[0];
  RotationMatrix[2][0]-=sinRotationAngle*l[1],RotationMatrix[2][1]+=sinRotationAngle*l[0];
}




















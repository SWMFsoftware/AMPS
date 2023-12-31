/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */


#include "pic.h"

//the object name and the names of the frames
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Mercury";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MERCURY";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="MSGR_MSO";

int Mercury::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=-1;
int Mercury::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=-1;
int Mercury::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=-1;

/*//parameters of the soruce processes
double Exosphere::SourceProcesses::ThermalDesorption::uThermal=1.85*eV2J;
double Exosphere::SourceProcesses::ThermalDesorption::VibrationalFrequency=1.0E13;

double Exosphere::SourceProcesses::PhotonStimulatedDesorption::PhotonFlux_1AU=3.0E15*1.0E+4;
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::CrossSection=2.0E-20*1.0E-4;
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::minInjectionEnergy=pow(10.0,2)*_NA__MASS_/2.0;
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::maxInjectionEnergy=pow(10.0E3,2)*_NA__MASS_/2.0;

double Exosphere::SourceProcesses::ImpactVaporization::SourceRate=2.6E23;
double Exosphere::SourceProcesses::ImpactVaporization::HeliocentricDistance=0.387098*_AU_;
double Exosphere::SourceProcesses::ImpactVaporization::SourceRatePowerIndex=1.9;
double Exosphere::SourceProcesses::ImpactVaporization::SourceTemperature=2500.0;

double Exosphere::SourceProcesses::SolarWindSputtering::Yield=0.1;
double Exosphere::SourceProcesses::SolarWindSputtering::minInjectionEnergy=pow(10.0,2)*_NA__MASS_/2.0;
double Exosphere::SourceProcesses::SolarWindSputtering::maxInjectionEnergy=pow(10.0E3,2)*_NA__MASS_/2.0;


//typical parameters of solar wind
const double Exosphere::swVelocity_Typical[3]={-420.0E3,0.0,0.0};
const double Exosphere::swB_Typical[3]={-12.9E-9,4.71E-9,10.29E-9};
const double Exosphere::swTemperature_Typical=0.174e6,Exosphere::swNumberDensity_Typical=60.0E6;


//interaction with the surface
const double Exosphere::SurfaceInteraction::AccomodationCoefficient=0.1;*/

//energy distribution of particles injected via ion sputtering
double Exosphere::SourceProcesses::SolarWindSputtering::EnergyDistributionFunction(double e,int *spec) {
  static const double Ee=0.015*eV2J;

  return e/pow(Ee+e,2);
}

//sticking probability of sodium atoms
double Exosphere::SurfaceInteraction::StickingProbability(int spec, double& ReemissionParticleFraction,double Temp) {
  double res;

  ReemissionParticleFraction=1.0;

  switch (spec) {
  case _NA_SPEC_: case _NA_PLUS_SPEC_:
    if (Temp<300.0) res=1.0;
    else if (Temp<650.0) res=1.0-(Temp-300.0)/350.0;
    else res=0.0;

    break;
  case _H_PLUS_SPEC_: case _HE_2PLUS_SPEC_:
    res=1.0;
    ReemissionParticleFraction=0.0;

    break;
  default:
    exit(__LINE__,__FILE__,"Error: unknown species");
  }

  return res;
}


//surface temeprature of the planet
double Exosphere::GetSurfaceTemperature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {
  static const double Tn=110.0;
  static const double Td0_Aphelion=590.0,Td0_Perihelion=725.0;
  static const double TAA_Aphelion=Pi,TAA_Perihelion=0.0;
  static const double alpha=(Td0_Aphelion-Td0_Perihelion)/(TAA_Aphelion-TAA_Perihelion);

  double Td,Angle;

  Angle=(OrbitalMotion::TAA<Pi) ? OrbitalMotion::TAA : 2.0*Pi-OrbitalMotion::TAA;
  Td=Td0_Perihelion+alpha*(Angle-TAA_Perihelion);

  return (CosSubSolarAngle>0.0) ? Tn+(Td-Tn)*pow(CosSubSolarAngle,0.25) : Tn;
}

//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\",  \"Mean Speed Along the Line of Sight(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec),PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=2;
  }

  if (_NA_SPEC_>=0) Mercury::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=2*_NA_SPEC_;

  //sodium emission
  if (_NA_SPEC_>=0) {
    Mercury::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=nVariables;
    Mercury::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=nVariables+1;

    if (vlist!=NULL) sprintf(vlist,"%s,  \"Sodium Emission(5891_58A)\",  \"Sodium Emission(5897_56A)\"",vlist);
    nVariables+=2;
  }

  return nVariables;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  int spec,cnt=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (res[cnt]>0.0) res[cnt+1]/=res[cnt];
    cnt+=2;
  }
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd,cnt=0,spec;
  double NumberDensity;

  static unsigned long int nCallCounter=0;

  //DEBUGGING: BEGIN 
  ++nCallCounter;  
  
  if (nCallCounter==18854) {
    double a=0.0;

    a+=2.454;
  }


  //DEBUGGING: END 
   

  nd=PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node);
  for (i=0;i<resLength;i++) res[i]=0.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //get the local density number
    NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(spec);
    res[cnt++]=NumberDensity;
    res[cnt++]=NumberDensity*node->block->GetCenterNode(nd)->GetMeanParticleSpeed(spec);
  }

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    double BulkVelocity_SO[3],v_LOCAL_SO_FROZEN[3],rHeliocentric,vHeliocentric;

    node->block->GetCenterNode(nd)->GetBulkVelocity(BulkVelocity_SO,spec);

    v_LOCAL_SO_FROZEN[0]=Exosphere::vObject_SO_FROZEN[0]+BulkVelocity_SO[0]+
        Exosphere::RotationVector_SO_FROZEN[1]*x[2]-Exosphere::RotationVector_SO_FROZEN[2]*x[1];

    v_LOCAL_SO_FROZEN[1]=Exosphere::vObject_SO_FROZEN[1]+BulkVelocity_SO[1]-
        Exosphere::RotationVector_SO_FROZEN[0]*x[2]+Exosphere::RotationVector_SO_FROZEN[2]*x[0];

    v_LOCAL_SO_FROZEN[2]=Exosphere::vObject_SO_FROZEN[2]+BulkVelocity_SO[2]+
        Exosphere::RotationVector_SO_FROZEN[0]*x[1]-Exosphere::RotationVector_SO_FROZEN[1]*x[0];

    rHeliocentric=sqrt(pow(x[0]-Exosphere::xObjectRadial,2)+(x[1]*x[1])+(x[2]*x[2]));
    vHeliocentric=(
        (v_LOCAL_SO_FROZEN[0]*(x[0]-Exosphere::xObjectRadial))+
        (v_LOCAL_SO_FROZEN[1]*x[1])+(v_LOCAL_SO_FROZEN[2]*x[2]))/rHeliocentric;


    //brightness of the exospheric sodium
    if (spec==_NA_SPEC_) {
      //check if the point is outside of Mercury's shadow
      if ( /*(Moon::EarthShadowCheck(x)==false) &&*/ ((x[0]>0.0)||(x[1]*x[1]+x[2]*x[2]>_RADIUS_(_MERCURY_)*_RADIUS_(_MERCURY_))) ) {
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5891_58A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5897_56A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
      }
      else res[cnt++]=0.0,res[cnt++]=0.0;
    }

  }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
   PIC::Debugger::CatchOutLimitValue(res,cnt,__LINE__,__FILE__);
#endif
#endif 

  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
}


//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble EphemerisTime) {
  return GetTAA("Mercury","Sun",_MASS_(_SUN_),EphemerisTime);
}

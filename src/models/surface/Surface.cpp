//$Id$

/*
 * Surface.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: vtenishe
 */

#include <cmath>
#include <cstdio>

#include "Surface.h"

//the model of each species: initially the compile-time choice of Surface.dfn
namespace {
  struct cInitSelectedModels {
    cInitSelectedModels() {
      for (int s=0;s<PIC::nTotalSpecies;s++) Surface::Model::Selected[s]=_SURFACE__GAS_SURFACE_INTERACTION_MODEL_;
    }
  };
}

int Surface::Model::Selected[PIC::nTotalSpecies];
static cInitSelectedModels InitSelectedModels;

const char* Surface::Model::Name(int model) {
  switch (model) {
  case Diffuse: return "diffuse";
  case Specular: return "specular";
  case CLL: return "cll";
  case Maxwell: return "maxwell";
  case QuasiSpecularCLL: return "quasi_specular_cll";
  case Adsorption: return "adsorption";
  case QuasiSpecularNocilla: return "nocilla";
  case QuasiSpecularSchamberg: return "schamberg";
  }

  return "unknown";
}

int Surface::Model::Code(const std::string& name) {
  for (int model=0;model<nTotalModels;model++) if (name==Name(model)) return model;
  return -1;
}

bool Surface::Model::Set(int spec,int model,std::string *error) {
  if ((model<0)||(model>=nTotalModels)) {
    if (error!=NULL) *error="unknown gas-surface interaction model code "+std::to_string(model);
    return false;
  }

  if (spec>=PIC::nTotalSpecies) {
    if (error!=NULL) *error="species "+std::to_string(spec)+" does not exist";
    return false;
  }

  for (int s=0;s<PIC::nTotalSpecies;s++) if ((spec<0)||(s==spec)) Selected[s]=model;
  return true;
}

bool Surface::Model::UsesWallTemperature(int spec) {
  switch (Selected[spec]) {
  case Diffuse: case CLL: case QuasiSpecularCLL:
    return true;
  case Maxwell:
    return MaxwellReflection::AccommodationCoefficient[spec]>0.0;
  case QuasiSpecularSchamberg:
    return QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[spec]>0.0;
  case QuasiSpecularNocilla:
    return QuasiSpecularReflection::Nocilla::Mode[spec]==QuasiSpecularReflection::Nocilla::AccommodationCoefficients;
  }

  return false;
}

//a value is in [lo,hi] (false for NaN, i.e. for an unset parameter)
static bool InRange(double x,double lo,double hi) {return (x>=lo)&&(x<=hi);}

static bool Fail(std::string *error,const std::string& message) {
  if (error!=NULL) *error=message;
  return false;
}

bool Surface::CheckParameters(int spec,int model,std::string *error) {
  //the message is built only when a value is not valid (this is called at every interaction of some models)
  auto bad=[&](const char *name,double value,const char *range) {
    char buf[400];

    snprintf(buf,sizeof(buf),"model '%s': %s of species %s is %s (allowed: %s)",Model::Name(model),name,
      PIC::MolecularData::GetChemSymbol(spec),std::isnan(value) ? "not set" : std::to_string(value).c_str(),range);
    return Fail(error,buf);
  };

  switch (model) {
  case Model::Diffuse: case Model::Specular: case Model::Adsorption:
    break;

  case Model::Maxwell:
    if (!InRange(MaxwellReflection::AccommodationCoefficient[spec],0.0,1.0)) return bad("the accommodation coefficient",MaxwellReflection::AccommodationCoefficient[spec],"0..1");
    break;

  case Model::CLL:
    if (!InRange(CLL::NormalPartKineticEnergyAccommodationCoefficient[spec],0.0,1.0)) return bad("the normal energy accommodation coefficient alpha_n",CLL::NormalPartKineticEnergyAccommodationCoefficient[spec],"0..1");
    if (!InRange(CLL::TangentialMomentumAccommodationCoefficient[spec],0.0,2.0)) return bad("the tangential momentum accommodation coefficient sigma_t",CLL::TangentialMomentumAccommodationCoefficient[spec],"0..2");
    break;

  case Model::QuasiSpecularCLL:
    if (!InRange(QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient[spec],0.0,1.0)) return bad("the normal energy accommodation coefficient alpha_n",QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient[spec],"0..1");
    if (!InRange(QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient[spec],0.0,1.0)) return bad("the tangential momentum accommodation coefficient sigma_t",QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient[spec],"0..1");
    break;

  case Model::QuasiSpecularNocilla:
    if (QuasiSpecularReflection::Nocilla::Mode[spec]==QuasiSpecularReflection::Nocilla::AccommodationCoefficients) {
      if (!InRange(QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[spec],0.0,1.0)) return bad("the normal momentum accommodation coefficient sigma'",QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[spec],"0..1");
      if (!InRange(QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation[spec],0.0,1.0)) return bad("the tangential momentum accommodation coefficient sigma",QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation[spec],"0..1");
      if (!InRange(QuasiSpecularReflection::Nocilla::EnergyAccommodation[spec],0.0,1.0)) return bad("the energy accommodation coefficient epsilon",QuasiSpecularReflection::Nocilla::EnergyAccommodation[spec],"0..1");
      break;
    }

    if (!InRange(QuasiSpecularReflection::Nocilla::DriftSpeed[spec],0.0,HUGE_VAL)) return bad("the drift speed U_r [m/s]",QuasiSpecularReflection::Nocilla::DriftSpeed[spec],">= 0");
    if (!InRange(QuasiSpecularReflection::Nocilla::DriftAngle[spec],0.0,0.5*Pi)) return bad("the drift angle theta_r [rad]",QuasiSpecularReflection::Nocilla::DriftAngle[spec],"0..pi/2");
    if (!InRange(QuasiSpecularReflection::Nocilla::Temperature[spec],1.0E-300,HUGE_VAL)) return bad("the re-emission temperature T_r [K]",QuasiSpecularReflection::Nocilla::Temperature[spec],"> 0");
    break;

  case Model::QuasiSpecularSchamberg:
    if (!InRange(QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[spec],0.0,1.0)) return bad("the energy accommodation coefficient alpha",QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[spec],"0..1");
    if (!(QuasiSpecularReflection::Schamberg::Nu[spec]>0.0)) return bad("the exponent nu",QuasiSpecularReflection::Schamberg::Nu[spec],"> 0, inf allowed");
    if (!InRange(QuasiSpecularReflection::Schamberg::BeamHalfAngle[spec],0.0,0.5*Pi)) return bad("the beam half-angle phi0 [rad]",QuasiSpecularReflection::Schamberg::BeamHalfAngle[spec],"0..pi/2");
    break;

  default: {
      char buf[200];

      snprintf(buf,sizeof(buf),"unknown gas-surface interaction model code %i of species %s",model,PIC::MolecularData::GetChemSymbol(spec));
      return Fail(error,buf);
    }
  }

  return true;
}

int Surface::ParticleInteractionProcessor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int res;

  //call the particle/surface interation model of the species
  switch (Model::Selected[PIC::ParticleBuffer::GetI(ptr)]) {
  case Model::Specular:
    res=SpecularReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  case Model::Diffuse:
    res=DiffuseReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  case Model::CLL:
    res=CLL::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  case Model::Maxwell:
    res=MaxwellReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  case Model::QuasiSpecularCLL:
    res=QuasiSpecularReflection::CLL::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  case Model::Adsorption:
    res=Adsorption::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  case Model::QuasiSpecularNocilla:
    res=QuasiSpecularReflection::Nocilla::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  case Model::QuasiSpecularSchamberg:
    res=QuasiSpecularReflection::Schamberg::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the gas-surface interaction model of the species is unknown");
  }

  return res;
}

double Surface::GetSurfaceTemperature(CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double res;

#if _SURFACE__TEMPERATURE_MODEL_ == _SURFACE__TEMPERATURE_MODEL__ISOTHERMAL_
  res=Surface::Temperature::Isothremal::Temp;

  //a wall at 0 K (the initial value: not set) would re-emit every particle at rest
  if (!(res>0.0)||std::isinf(res)) {
    char msg[300];

    snprintf(msg,sizeof(msg),"Error: the wall temperature is not set or invalid (Surface::Temperature::Isothremal::Temp=%e K); "
        "set it, e.g. with 'wall_temperature_k' in the runtime block 'physics: gas_surface_interaction'",res);
    exit(__LINE__,__FILE__,msg);
  }
#else
  exit(__LINE__,__FILE__,"Error: the option is unknown");
#endif

  return res;
}



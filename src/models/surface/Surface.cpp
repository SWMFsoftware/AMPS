//$Id$

/*
 * Surface.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: vtenishe
 */

#include "Surface.h"


int Surface::ParticleInteractionProcessor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int res;

  //call approptiate particle/surface interation model
  #if _SURFACE__GAS_SURFACE_INTERACTION_MODEL_ == _SURFACE__GAS_SURFACE_INTERACTION_MODEL__SPECULAR_REFLECTION_
  res=SpecularReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
  #elif _SURFACE__GAS_SURFACE_INTERACTION_MODEL_ == _SURFACE__GAS_SURFACE_INTERACTION_MODEL__DIFFUSE_REFLECTION_
  res=DiffuseReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
  #elif _SURFACE__GAS_SURFACE_INTERACTION_MODEL_ == _SURFACE__GAS_SURFACE_INTERACTION_MODEL__CCL_
  res=CLL::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
  #else
  exit(__LINE__,__FILE__,"Error: the option is unknown");
  #endif

  return res;
}

double Surface::GetSurfaceTemperature(CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double res;

#if _SURFACE__TEMPERATURE_MODEL_ == _SURFACE__TEMPERATURE_MODEL__ISOTHERMAL_
  res=Surface::Temperature::Isothremal::Temp;
#else
  exit(__LINE__,__FILE__,"Error: the option is unknown");
#endif

  return res;
}



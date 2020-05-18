/*
 * mover.cpp
 *
 *  Created on: May 16, 2020
 *      Author: vtenishe
 */

#include "sep.h"

int SEP::ParticleMover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double xInit[3];
  int res;

  PIC::ParticleBuffer::GetX(xInit,ptr);

  res=PIC::Mover::Relativistic::Boris(ptr,dtTotal,startNode);

  if (res==_PARTICLE_MOTION_FINISHED_) {
    //appply particle scatering model if needed 
    double x[3],v[3],speed; 
    double mean_free_path,l2,e,r;
    int idim,spec;

    PIC::ParticleBuffer::GetX(x,ptr);
    PIC::ParticleBuffer::GetV(v,ptr);
    spec=PIC::ParticleBuffer::GetI(ptr);
  
    speed=Vector3D::Length(v);
    e=Relativistic::Speed2E(speed,PIC::MolecularData::GetMass(spec));
    r=Vector3D::Length(x); 
  
    mean_free_path=3.4E9; 

    for (l2=0.0,idim=0;idim<3;idim++) {
      double t;

      t=xInit[idim]-x[idim];
      l2+=t*t;
    }

    if ((1.0-exp(-sqrt(l2)/mean_free_path))>rnd()) {
      //the scattering even occured
      Vector3D::Distribution::Uniform(v,speed);
   
      PIC::ParticleBuffer::SetV(v,ptr);
    }
  }

  return res;
} 



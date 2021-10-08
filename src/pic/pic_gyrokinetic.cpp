//gyrokinetic approaximation 

#include "pic.h"

int PIC::GYROKINETIC::DriftVelocityOffset=-1;

//===================================================================
//Init the model/reserve the space in the particle state vector to stor ethe drift velocity so it is no neer to recalcualte it
void PIC::GYROKINETIC::Init() {
  //request data in the particle state vector
  long int offset;

  PIC::ParticleBuffer::RequestDataStorage(offset,3*sizeof(double));
  DriftVelocityOffset=offset;
}

//===================================================================
//Particle mover 
int PIC::GYROKINETIC::Mover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int res,spec;

  spec=PIC::ParticleBuffer::GetI(ptr);

  switch (spec) {
  case _ELECTRON_SPEC_: 
    //use the gyrokinetic approximation
    res=PIC::Mover::GuidingCenter::Mover_FirstOrder(ptr,dtTotal,startNode); 
    break;

  default: 
    //Lapenta2017 is used for the ions 
    res=PIC::Mover::Lapenta2017(ptr,dtTotal,startNode);
  } 

  return res;    
} 






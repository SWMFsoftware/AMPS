
#include "sep.h"

int SEP::Offset::Momentum=-1;
int SEP::Offset::CosPitchAngle=-1;
int SEP::Offset::p_par=-1;
int SEP::Offset::p_norm=-1;

cInternalSphericalData* SEP::InnerBoundary=NULL;

void SEP::RequestParticleData() {
  int offset;

  switch (_SEP_MOVER_) {
  case _SEP_MOVER_HE_2019_AJL_:
    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::Momentum=offset;

    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::CosPitchAngle=offset;
    break;
  case _SEP_MOVER_BOROVIKOV_2019_ARXIV_:
    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::p_par=offset;

    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::p_norm=offset;
    break;
  }
}

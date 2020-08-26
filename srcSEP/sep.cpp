
#include "sep.h"

int SEP::Offset::Momentum=-1;
int SEP::Offset::CosPitchAngle=-1;

void SEP::RequestParticleData() {
  long int offset;

  PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
  Offset::Momentum=offset;

  PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
  Offset::CosPitchAngle=offset;
}

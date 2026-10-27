#include "pic.h"
#include "Earth.h"
#include "specfunc.h"


double Geospace::OutsideParticleFlux(double ParticleEnergy) {
  return 4.15*pow(ParticleEnergy,1.72); //Kress-2013-ASR
}

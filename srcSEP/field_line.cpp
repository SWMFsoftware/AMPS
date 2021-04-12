
#include "pic.h"
#include "sep.h"
#include "amps2swmf.h"


long int SEP::FieldLine::InjectParticlesSingleFieldLine(int spec,int iFieldLine) {
  int iShockFieldLine,npart;
  double anpart,p[3],ParticleWeightCorrectionFactor;
  int nInjectedParticles=0;


  //determine the filed line to inject particles
  if ((iShockFieldLine=AMPS2SWMF::iShockWaveSegmentTable[iFieldLine])==-1) return 0; 

  //determine the number of particles to inject 
  anpart=23.0;
  npart=(int)anpart;
  if (anpart-npart>rnd()) npart++; 
  
  
  for (int i=0;i<npart;i++) {
    //determine the momentum of the particles 
    Vector3D::Distribution::Uniform(p,1000.0E3*PIC::MolecularData::GetMass(spec));
    ParticleWeightCorrectionFactor=1.0;
  
    if (PIC::FieldLine::InjectParticle_default(spec,p,ParticleWeightCorrectionFactor,iFieldLine,iShockFieldLine)!=-1) nInjectedParticles++;  
  }

  return nInjectedParticles;
} 
   
long int SEP::FieldLine::InjectParticles() {
  long int res=0;

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (int iFieldLine=0;iFieldLine<PIC::FieldLine::nFieldLine;iFieldLine++) {
    res+=InjectParticlesSingleFieldLine(spec,iFieldLine);
  }

  return res;
}
 


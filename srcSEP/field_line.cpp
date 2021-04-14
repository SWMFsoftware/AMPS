
#include "pic.h"
#include "sep.h"
#include "amps2swmf.h"


long int SEP::FieldLine::InjectParticlesSingleFieldLine(int spec,int iFieldLine) {
  int iShockFieldLine,npart;
  double anpart,p[3],ParticleWeightCorrectionFactor;
  int nInjectedParticles=0;


  //determine the filed line to inject particles
  iShockFieldLine=0; 

  #ifdef _SEP_SHOCK_LOCATION_COUPLER_TABLE_
  #if _SEP_SHOCK_LOCATION_COUPLER_TABLE_ == _PIC_MODE_ON_ 
  if ((iShockFieldLine=AMPS2SWMF::iShockWaveSegmentTable[iFieldLine])==-1) return 0; 
  #endif
  #endif

  //determine the number of particles to inject 
  anpart=23.0;
  npart=(int)anpart;
  if (anpart-npart>rnd()) npart++; 
  
//    for (int i=0;i<npart;i++) {
//    //determine the momentum of the particles 
//    Vector3D::Distribution::Uniform(p,1000.0E3*PIC::MolecularData::GetMass(spec));



  double emin=0.1*MeV2J;
  double emax=100.0*MeV2J;

  double s=4.0;
  double q=3.0*s/(3-1.0);

  double pAbs,pmin,pmax,speed,pvect[3];
  double mass=PIC::MolecularData::GetMass(spec);

  pmin=Relativistic::Energy2Momentum(emin,mass);
  pmax=Relativistic::Energy2Momentum(emax,mass);

  double cMin=pow(pmin,-q);

  speed=Relativistic::E2Speed(emin,PIC::MolecularData::GetMass(spec));
  pmin=Relativistic::Speed2Momentum(speed,mass);

  speed=Relativistic::E2Speed(emax,PIC::MolecularData::GetMass(spec));
  pmax=Relativistic::Speed2Momentum(speed,mass);


  double A0=pow(pmin,-q+1.0);
  double A=pow(pmax,-q+1.0)-A0;

  double WeightNorm=pow(pmin,-q);

  for (int i=0;i<npart;i++) {

    pAbs=pmin+rnd()*(pmax-pmin);
    ParticleWeightCorrectionFactor=pow(pAbs,-q)/WeightNorm;

    Vector3D::Distribution::Uniform(p,pAbs);

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
 


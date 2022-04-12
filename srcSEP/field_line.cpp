
#include "pic.h"
#include "sep.h"
#include "amps2swmf.h"


int SEP::FieldLine::InjectionParameters::nParticlesPerIteration=30;
double SEP::FieldLine::InjectionParameters::PowerIndex=4.0;
double SEP::FieldLine::InjectionParameters::emin=0.1,SEP::FieldLine::InjectionParameters::emax=500;

int SEP::FieldLine::InjectionParameters::InjectLocation=SEP::FieldLine::InjectionParameters::_InjectInputFileAMPS;


long int SEP::FieldLine::InjectParticleFieldLineBeginning(int spec,int iFieldLine) {
  namespace FL = PIC::FieldLine;

  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  int nInjectedParticles=0;
  int npart;
  double l[3],pAbs,p[3],ParticleWeightCorrectionFactor=1.0;

  npart=100;
  pAbs=Relativistic::Energy2Momentum(100.0*MeV2J,PIC::MolecularData::GetMass(spec));

  FL::FieldLinesAll[iFieldLine].GetSegment(0)->GetDir(l); 

  for (int i=0;i<npart;i++) {
    //generate a particle
    Vector3D::Distribution::Uniform(p,pAbs);

    if (Vector3D::DotProduct(p,l)<0.0) for (int idim=0;idim<3;idim++) p[idim]=-p[idim];
    
    if (PIC::FieldLine::InjectParticle_default(spec,p,ParticleWeightCorrectionFactor,iFieldLine,0)!=-1) nInjectedParticles++;
  }
   
  return nInjectedParticles;
}

long int SEP::FieldLine::InjectParticlesSingleFieldLine(int spec,int iFieldLine) {
  namespace FL = PIC::FieldLine;

  int iShockFieldLine,npart;
  double anpart,p[3],ParticleWeightCorrectionFactor;
  int nInjectedParticles=0;


  //determine the filed line to inject particles
  iShockFieldLine=0; 

  if (InjectionParameters::InjectLocation==InjectionParameters::_InjectInputFileAMPS) {
    #ifdef _SEP_SHOCK_LOCATION_COUPLER_TABLE_
    #if _SEP_SHOCK_LOCATION_COUPLER_TABLE_ == _PIC_MODE_ON_ 
    if ((iShockFieldLine=AMPS2SWMF::ShockData[iFieldLine].iSegmentShock)==-1) return 0; 
    #endif
    #endif
  }
  else {
    exit(__LINE__,__FILE__);
  }


  //determine the radiaus of the magnetic tube at the middle of the magnetic tube
  FL::cFieldLineSegment* Segment=FL::FieldLinesAll[iFieldLine].GetSegment(iShockFieldLine); 

  if (Segment==NULL) return 0;
 
  //determine the volume swept by the shock wave during the time step 
  double xBegin[3],xEnd[3],xMiddle[3],rMiddle,rMiddleTube,xFirstFieldLine[3];

  Segment->GetBegin()->GetX(xBegin);
  Segment->GetEnd()->GetX(xEnd);

  for (int idim=0;idim<3;idim++) xMiddle[idim]=0.5*(xBegin[idim]+xEnd[idim]);

  rMiddle=Vector3D::Length(xMiddle);

  FL::FieldLinesAll[iFieldLine].GetSegment(0)->GetBegin()->GetX(xFirstFieldLine);  
  rMiddleTube=pow(rMiddle/Vector3D::Length(xFirstFieldLine),2); 

  //velocity of the shock wave
  double *v_sw,vol;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::Search::FindBlock(xMiddle);

  v_sw=Segment->GetBegin()->GetDatum_ptr(FL::DatumAtVertexPlasmaVelocity);
//  vol=node->block->GetLocalTimeStep(spec)*Vector3D::Length(v_sw)*rMiddleTube;

vol=node->block->GetLocalTimeStep(spec)*AMPS2SWMF::ShockData[iFieldLine].ShockSpeed*rMiddleTube;


  //determine the number of particles to inject 
  const double InjectionEfficiency=3.4E-4; //Sokolov-2004-AJ 


  double t_sw_begin,t_sw_end; //=Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaTemperature); 
  double n_sw_begin,n_sw_end; //=Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaDensity); 
  double p_inj=sqrt(2.0*_AMU_*1.0E4*ElectronCharge);  

  Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_begin);
  Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_begin); 

  Segment->GetEnd()->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_end);
  Segment->GetEnd()->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_end);

//  double t=sqrt(2.0*_AMU_*t_sw);
//  anpart=vol*InjectionEfficiency/(8.0*Pi)*n_sw/pow(t,3)*pow(t/p_inj,5);  


  n_sw_end=AMPS2SWMF::ShockData[iFieldLine].DownStreamDensity;

  anpart=vol*InjectionEfficiency*n_sw_end;
  anpart/=node->block->GetLocalParticleWeight(spec);


  double GlobalWeightCorrectionFactor=1.0;

  if (anpart<InjectionParameters::nParticlesPerIteration) {
    GlobalWeightCorrectionFactor=anpart/InjectionParameters::nParticlesPerIteration;
    anpart=InjectionParameters::nParticlesPerIteration;
  }

 
  npart=(int)anpart;
  if (anpart-npart>rnd()) npart++; 
  
  double emin=InjectionParameters::emin*MeV2J;
  double emax=InjectionParameters::emax*MeV2J;

  double s=InjectionParameters::PowerIndex;
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
    ParticleWeightCorrectionFactor=pow(pAbs,-q)/WeightNorm*GlobalWeightCorrectionFactor;

    Vector3D::Distribution::Uniform(p,pAbs);

    if (PIC::FieldLine::InjectParticle_default(spec,p,ParticleWeightCorrectionFactor,iFieldLine,iShockFieldLine)!=-1) nInjectedParticles++;  
  }

  return nInjectedParticles;
} 
   
long int SEP::FieldLine::InjectParticles() {
  long int res=0;

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (int iFieldLine=0;iFieldLine<PIC::FieldLine::nFieldLine;iFieldLine++) {
    switch (_SEP_FIELD_LINE_INJECTION_) {
    case _SEP_FIELD_LINE_INJECTION__SHOCK_:
      res+=InjectParticlesSingleFieldLine(spec,iFieldLine);
      break;
    case _SEP_FIELD_LINE_INJECTION__BEGINNIG_:
      res+=InjectParticleFieldLineBeginning(spec,iFieldLine);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
  }

  return res;
}
 



#include "pic.h"
#include "sep.h"
#include "amps2swmf.h"


int SEP::FieldLine::InjectionParameters::nParticlesPerIteration=30;
double SEP::FieldLine::InjectionParameters::PowerIndex=4.0;
double SEP::FieldLine::InjectionParameters::emin=0.1,SEP::FieldLine::InjectionParameters::emax=500;
double SEP::FieldLine::InjectionParameters::InjectionEfficiency=3.4E-4; //Sokolov-2004-AJ 

double SEP::FieldLine::InjectionParameters::ConstEnergyInjectionValue=0.0;
double SEP::FieldLine::InjectionParameters::ConstSpeedInjectionValue=0.0;
double SEP::FieldLine::InjectionParameters::ConstMuInjectionValue=0.5;

#if _SEP_FIELD_LINE_INJECTION_ == _SEP_FIELD_LINE_INJECTION__SHOCK_
int SEP::FieldLine::InjectionParameters::InjectLocation=SEP::FieldLine::InjectionParameters::_InjectShockLocations;
#else 
int SEP::FieldLine::InjectionParameters::InjectLocation=SEP::FieldLine::InjectionParameters::_InjectBegginingFL;
#endif



int SEP::FieldLine::InjectionParameters::InjectionMomentumModel=SEP::FieldLine::InjectionParameters::_tenishev2005aiaa;


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

long int InjectSolarWindIons(int spec,int iFieldLine) {
  namespace FL = PIC::FieldLine;

  double InjectionArea,n_sw,t_sw,v_sw[3];
  auto Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();
  FL::cFieldLineVertex* FirstVertex=Segment->GetBegin();

  //determine the parameters of of the solar wind at the beginning of the field line  
  FirstVertex->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw);
  FirstVertex->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw);
  FirstVertex->GetPlasmaVelocity(v_sw);
  
  InjectionArea=Pi*pow(SEP::FieldLine::MagneticTubeRadius(FirstVertex->GetX(),iFieldLine),2); 

  //inject model partiles 
  return PIC::FieldLine::InjectMaxwellianLineBeginning(spec,n_sw,t_sw,v_sw,InjectionArea,iFieldLine,200);
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
    switch (InjectionParameters::InjectLocation) {
    case InjectionParameters::_InjectShockLocations:
      #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
      if (AMPS2SWMF::ShockData==NULL) {
        exit(__LINE__,__FILE__,"Error: the shock location table is not allocated");
      }
      else {
        if ((iShockFieldLine=AMPS2SWMF::ShockData[iFieldLine].iSegmentShock)==-1) return 0;
      }
      #else 
      iShockFieldLine=0;
      #endif
  
      break;
    case  InjectionParameters::_InjectBegginingFL:
      iShockFieldLine=0;
      break;
    }
  }

  //determine the radiaus of the magnetic tube at the middle of the magnetic tube
  FL::cFieldLineSegment* Segment=FL::FieldLinesAll[iFieldLine].GetSegment(iShockFieldLine); 

  if (Segment==NULL) return 0;
 
  //determine the volume swept by the shock wave during the time step 
  double xBegin[3],xEnd[3],xMiddle[3],rMiddle,xFirstFieldLine[3];

  Segment->GetBegin()->GetX(xBegin);
  Segment->GetEnd()->GetX(xEnd);

  for (int idim=0;idim<3;idim++) xMiddle[idim]=0.5*(xBegin[idim]+xEnd[idim]);

  //velocity of the shock wave
  double vol;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::Search::FindBlock(xMiddle);

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  if (AMPS2SWMF::ShockData[iFieldLine].ShockSpeed>AMPS2SWMF::MinShockSpeed) {
    vol=node->block->GetLocalTimeStep(spec)*AMPS2SWMF::ShockData[iFieldLine].ShockSpeed*SEP::FieldLine::MagneticTubeRadius(xMiddle,iFieldLine);
  }
  else {
    if (AMPS2SWMF::MinShockSpeed==0.0) exit(__LINE__,__FILE__,"Error: AMPS2SWMF::MinShockSpeed is not set");

    vol=node->block->GetLocalTimeStep(spec)*AMPS2SWMF::MinShockSpeed*SEP::FieldLine::MagneticTubeRadius(xMiddle,iFieldLine);
  }
  #else 
    vol=node->block->GetLocalTimeStep(spec)*SEP::FieldLine::MagneticTubeRadius(xMiddle,iFieldLine);
  #endif


  //determine the number of particles to inject 
  double t_sw_begin,t_sw_end; //=Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaTemperature); 
  double n_sw_begin,n_sw_end; //=Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaDensity); 
  double p_inj=sqrt(2.0*_AMU_*1.0E4*ElectronCharge);  

  Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_begin);
  Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_begin); 

  Segment->GetEnd()->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_end);
  Segment->GetEnd()->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_end);


#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  n_sw_end=AMPS2SWMF::ShockData[iFieldLine].DownStreamDensity;

  anpart=vol*SEP::FieldLine::InjectionParameters::InjectionEfficiency*n_sw_end;
  anpart/=node->block->GetLocalParticleWeight(spec);
#else 
  n_sw_end=1.0;
  anpart=InjectionParameters::nParticlesPerIteration;
#endif

  double GlobalWeightCorrectionFactor=1.0;

  if (anpart<InjectionParameters::nParticlesPerIteration) {
    GlobalWeightCorrectionFactor=anpart/InjectionParameters::nParticlesPerIteration;
    anpart=InjectionParameters::nParticlesPerIteration;
  }

  //in case particle are injected at the beginning of the field line, the actual plasma density is not used -> set the particle weight == 1
  if (InjectionParameters::InjectLocation==InjectionParameters::_InjectBegginingFL) {
    GlobalWeightCorrectionFactor=1.0;
  }
 
  npart=(int)anpart;
  if (anpart-npart>rnd()) npart++; 
  
  auto GetMomentum_Tenishev2005AIAA = [&] (double *pAbsTable,double *WeightCorrectionTable,int nParticles) {
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

    for (int i=0;i<nParticles;i++) {
      pAbsTable[i]=pmin+rnd()*(pmax-pmin);
      WeightCorrectionTable[i]=pow(pAbsTable[i],-q)/WeightNorm*GlobalWeightCorrectionFactor;
    }
  }; 

  auto GetMomentum_Sokolov2004AJ = [&] (double *pAbsTable,double *WeightCorrectionTable,int nParticles) { 
    double pmin=sqrt(10.0*KeV2J*2.0*PIC::MolecularData::GetMass(spec));
    double e,r;

    double p_injection_min=Relativistic::Energy2Momentum(SEP::FieldLine::InjectionParameters::emin,PIC::MolecularData::GetMass(spec));
    double p_injection_max=Relativistic::Energy2Momentum(SEP::FieldLine::InjectionParameters::emax,PIC::MolecularData::GetMass(spec));

    double log_p_injection_min=log(p_injection_min);
    double log_p_injection_max=log(p_injection_max);

    for (int i=0;i<nParticles;i++) {
      pAbsTable[i]=exp(log_p_injection_min+rnd()*(log_p_injection_max-log_p_injection_min));  
      WeightCorrectionTable[i]=pow(p_injection_min/pAbsTable[i],InjectionParameters::PowerIndex);
      WeightCorrectionTable[i]*=1.0/pAbsTable[i];
    } 
  };

  double *pAbsTable=new double [npart];
  double *WeightCorrectionTable=new double [npart];
  double p_const;

  switch (InjectionParameters::InjectionMomentumModel) {
  case InjectionParameters::_tenishev2005aiaa: 
    GetMomentum_Tenishev2005AIAA(pAbsTable,WeightCorrectionTable,npart);
    break;
  case InjectionParameters::_sokolov2004aj:
    GetMomentum_Sokolov2004AJ(pAbsTable,WeightCorrectionTable,npart);
    break;
  case InjectionParameters::_const_speed:
    p_const=Relativistic::Speed2Momentum(SEP::FieldLine::InjectionParameters::ConstSpeedInjectionValue,PIC::MolecularData::GetMass(spec)); 

    for (int i=0;i<npart;i++) pAbsTable[i]=p_const,WeightCorrectionTable[i]=1.0;  
    break;
  case InjectionParameters::_const_energy:
     p_const=Relativistic::Energy2Momentum(SEP::FieldLine::InjectionParameters::ConstEnergyInjectionValue,PIC::MolecularData::GetMass(spec));

    for (int i=0;i<npart;i++) pAbsTable[i]=p_const,WeightCorrectionTable[i]=1.0;
    break;
  defaut:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }  

  for (int i=0;i<npart;i++) {
    if ((InjectionParameters::InjectionMomentumModel==InjectionParameters::_const_speed)||(InjectionParameters::InjectionMomentumModel==InjectionParameters::_const_energy)) {
      double p_parallel[3],p_norm[3];
      double l[3];
      double e0[3],e1[3],c,mu;

      Segment->GetDir(l);
      Vector3D::GetRandomNormFrame(e0,e1,l);

      mu=SEP::FieldLine::InjectionParameters::ConstMuInjectionValue;
      c=sqrt(1.0-mu*mu);

      for (int idim=0;idim<3;idim++) p[idim]=pAbsTable[i]*(mu*l[idim]+c*e0[idim]);
    }
    else {
      Vector3D::Distribution::Uniform(p,pAbsTable[i]);
    } 


    if (PIC::FieldLine::InjectParticle_default(spec,p,GlobalWeightCorrectionFactor*WeightCorrectionTable[i],iFieldLine,iShockFieldLine)!=-1) nInjectedParticles++;  
  }

  delete [] pAbsTable;
  delete [] WeightCorrectionTable;

  return nInjectedParticles;
} 
   
long int SEP::FieldLine::InjectParticles() {
  long int res=0;

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (int iFieldLine=0;iFieldLine<PIC::FieldLine::nFieldLine;iFieldLine++) {
    switch (SEP::FieldLine::InjectionParameters::InjectLocation) {
    case SEP::FieldLine::InjectionParameters::_InjectShockLocations:
      res+=InjectParticlesSingleFieldLine(spec,iFieldLine);
      break;
    
    case SEP::FieldLine::InjectionParameters::_InjectBegginingFL: 
      if (InjectionParameters::InjectionMomentumModel==SEP::FieldLine::InjectionParameters::_background_sw_temperature) {
        res+=InjectSolarWindIons(spec,iFieldLine);
      }
      else {
        res+=InjectParticleFieldLineBeginning(spec,iFieldLine);
      }
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
  }

  return res;
}
 


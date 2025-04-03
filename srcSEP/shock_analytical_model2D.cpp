#include "sep.h"
//analytic shock for 2D case 

//the shoch is represented with a schere 
cInternalSphericalData SEP::ParticleSource::ShockWaveSphere::ShockSurface;
double *SEP::ParticleSource::ShockWaveSphere::CompressionRatioTable=NULL;
double *SEP::ParticleSource::ShockWaveSphere::SourceRateTable=NULL;
int SEP::ParticleSource::ShockWaveSphere::nSurfaceElements=0;
cSingleVariableDiscreteDistribution<int> SEP::ParticleSource::ShockWaveSphere::ShockInjectionDistribution; 
bool SEP::ParticleSource::ShockWaveSphere::InitGenerationSurfaceElement=false; 


double SEP::ParticleSource::ShockWaveSphere::SphericalShockOpeningAngleLimit=Pi/4.0;

void SEP::ParticleSource::ShockWaveSphere::Init() {
  nSurfaceElements=ShockSurface.GetTotalSurfaceElementsNumber(); 

  ShockSurface.Radius=1.0E-5*_SUN__RADIUS_; 

  ShockSurface.OriginPosition[0]=_SUN__RADIUS_;
  ShockSurface.OriginPosition[1]=0.0;
  ShockSurface.OriginPosition[2]=0.0;

  CompressionRatioTable=new double [nSurfaceElements];
  SourceRateTable=new double [nSurfaceElements];
}

//flush elements of the source rate table 
void SEP::ParticleSource::ShockWaveSphere::Flush() {
  for (int i = 0; i < nSurfaceElements; ++i) {
    SourceRateTable[i]= 0.0; // Setting the z element to 0
  }

  InitGenerationSurfaceElement=false;
}

//solar wind density model 
double SEP::ParticleSource::ShockWaveSphere::GetSolarWindDensity(double *x) {
  return 5.0E6*_AU_*_AU_/Vector3D::DotProduct(x,x); 
}

//calcualte the source rate table 
double SEP::ParticleSource::ShockWaveSphere::GetTotalSourceRate() {
  double x[3],res=0.0,n,s,efficientcy;

  //Init the data buggers if needed 
  if (SourceRateTable==NULL) {
    Init();
  }

  Flush();

  SEP::ParticleSource::ShockWave::Tenishev2005::UpdateShockLocation();
  ShockSurface.Radius=SEP::ParticleSource::ShockWave::Tenishev2005::rShock-_SUN__RADIUS_;
  if (ShockSurface.Radius<1.0E-5*_SUN__RADIUS_) ShockSurface.Radius=1.0E-5*_SUN__RADIUS_;

  ShockSurface.Radius/=2.0;
  ShockSurface.OriginPosition[0]=_SUN__RADIUS_+ShockSurface.Radius;
  ShockSurface.OriginPosition[1]=0.0;
  ShockSurface.OriginPosition[2]=0.0;


  const double cosThetaLimit=cos(SphericalShockOpeningAngleLimit);


  if (PIC::ThisThread==0) cout << ShockSurface.Radius/_SUN__RADIUS_ << endl;

  for (int i = 0; i < nSurfaceElements; ++i) {
    ShockSurface.GetSurfaceElementMiddlePoint(x,i); 

    //limit the fraction of the sphere where the particle injection can occur
    if ((x[0]-ShockSurface.OriginPosition[0])/sqrt(pow(x[0]-ShockSurface.OriginPosition[0],2)+x[1]*x[1]+x[2]*x[2])<cosThetaLimit) {
      SourceRateTable[i]=0.0;
      continue;
    }

    if (Vector3D::DotProduct(x,x)<_SUN__RADIUS_*_SUN__RADIUS_) {
      SourceRateTable[i]=0.0;
    }
    else {
      n=GetSolarWindDensity(x);     
      s=SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio();
      efficientcy=(s-1.0)/s; 
     
      SourceRateTable[i]=n*efficientcy*SEP::ParticleSource::ShockWave::Tenishev2005::GetShockSpeed()*ShockSurface.GetSurfaceElementArea(i); 
      res+=SourceRateTable[i];
    }
  } 
  
  return res;
}


int SEP::ParticleSource::ShockWaveSphere::GetInjectionSurfaceElement(double *x) {
  int el;

  //set up the generator of the element index 
  if (InitGenerationSurfaceElement==false) {
    InitGenerationSurfaceElement=true;
     
    ShockInjectionDistribution.InitArray(SourceRateTable,ShockSurface.GetTotalSurfaceElementsNumber(),ShockSurface.GetTotalSurfaceElementsNumber()); 
  } 

  el=ShockInjectionDistribution.DistributeVariable(); 
  ShockSurface.GetSurfaceElementRandomPoint(x,el);

  return el;
}

long int SEP::ParticleSource::ShockWaveSphere::InjectionModel() {
  int res=0;
  double TotalInjectionRate,CompressionRatio;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;

  InitGenerationSurfaceElement=false;

  TotalInjectionRate=GetTotalSourceRate();
  TotalInjectionRate/=PIC::ParticleWeightTimeStep::GlobalParticleWeight[0];

  double emin=SEP::FieldLine::InjectionParameters::emin*MeV2J;
  double emax=SEP::FieldLine::InjectionParameters::emax*MeV2J;
  double mass=PIC::MolecularData::GetMass(0);

  double speed,pmin,pmax;
  speed=Relativistic::E2Speed(emin,PIC::MolecularData::GetMass(0));
  pmin=Relativistic::Speed2Momentum(speed,mass);

  //the loop to generate new particles 
  double s,dtCounter=0.0,dtTotal,x[3],v[3],p,m,q;
  int el;

  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  dtTotal=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
  #elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  dtTotal=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
  #else
  exit(__LINE__,__FILE__,"Error: the global time counter cannot be applied for this case");
  #endif


  //limit the number of injected limit the number of injected particle  
  const int nMaxInjectedParticle=100;
  double WeightCorrectionFactor=1.0;

  if (TotalInjectionRate*dtTotal>nMaxInjectedParticle) {
    double c=nMaxInjectedParticle/(TotalInjectionRate*dtTotal);


    TotalInjectionRate*=c;
    WeightCorrectionFactor=1.0/c;

    if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
      exit(__LINE__,__FILE__,"Error: set _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_");
    }
  }


  while ((dtCounter-=log(rnd())/TotalInjectionRate)<dtTotal) {
    //1. Generate new particle position 
    el=GetInjectionSurfaceElement(x);
    
    if (Vector3D::DotProduct(x,x)<_SUN__RADIUS_*_SUN__RADIUS_) continue; 

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode=PIC::Mesh::mesh->findTreeNode(x);
    if (startNode->block==NULL) continue;
    if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) continue;

    //2. generate particle velocity  
    s=CompressionRatioTable[el]; 
    q=3.0*s/(s-1.0);

    double misc1=pow(pmax,1-q);
    double misc2=pow(pmin,1-q);
    p=pow(rnd()*(misc1-misc2)+misc2,1.0/(1.0-q));  

    //convert to particle speed 
    speed=Relativistic::Momentum2Speed(p,_H__MASS_);
    Vector3D::Distribution::Uniform(v,speed);

    //check the direction of the new particle; it sould move outside of the shock center 
    double d=0.0;

    for (int i=0;i<3;i++) d+=v[i]*(x[i]-ShockSurface.OriginPosition[i]);  
    if (d<0.0) for (int i=0;i<3;i++) v[i]*=-1.0;

    //generate a new particle 
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle); 

    PIC::ParticleBuffer::SetX(x,newParticleData);
    PIC::ParticleBuffer::SetV(v,newParticleData);
    PIC::ParticleBuffer::SetI(0,newParticleData);

    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrectionFactor,newParticleData);
    #endif

    //apply condition of tracking the particle
    #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    PIC::ParticleTracker::InitParticleID(newParticleData);
    PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xInit,xFinal,spec,newParticleData,(void*)node);
    #endif


    SEP::ParticleMoverPtr(newParticle,rnd()*dtTotal,startNode); 
    res++;
  } 

  return res;
}









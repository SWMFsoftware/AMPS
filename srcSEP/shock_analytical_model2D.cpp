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
  double res;

  switch ( SEP::ParticleSource::ShockWaveSphere::SolarWindDensityMode) {
  case SEP::ParticleSource::ShockWaveSphere::SolarWindDensityMode_analytic: 
    res=5.0E6*_AU_*_AU_/Vector3D::DotProduct(x,x); 
    break;
  case SEP::ParticleSource::ShockWaveSphere::SolarWindDensityMode_swmf: 
  default:
    exit(__LINE__,__FILE__,"Error: not implemented");
  }

  return res;
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
      CompressionRatioTable[i]=0.0;
      SourceRateTable[i]=0.0;
      continue;
    }

    if (Vector3D::DotProduct(x,x)<_SUN__RADIUS_*_SUN__RADIUS_) {
      CompressionRatioTable[i]=0.0;
      SourceRateTable[i]=0.0;
    }
    else {
      n=GetSolarWindDensity(x);     
      s=SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio();
      efficientcy=(s-1.0)/s; 
     
      CompressionRatioTable[i]=s;
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

//__attribute__((optimize("O0")))
long int SEP::ParticleSource::ShockWaveSphere::InjectionModel() {
  int res=0;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;

  InitGenerationSurfaceElement=false;
  const double PhysicalInjectionRate=GetTotalSourceRate();

  // GetTotalSourceRate() is a physical rate in particles/s.  Reject corrupt
  // model output before dividing by any species' macroparticle weight.  A
  // zero rate is a valid quiet-shock state and must return without evaluating
  // the exponential waiting-time expression below (which divides by rate).
  if (!isfinite(PhysicalInjectionRate) || PhysicalInjectionRate<0.0)
    exit(__LINE__,__FILE__,"shock source returned an invalid physical injection rate");
  if (PhysicalInjectionRate==0.0) return 0;


  auto GetMomentum_Tenishev2005AIAA = [&] (double& pAbs,double& WeightCorrection,double CompressionRatio,int spec) {
    double emin=SEP::FieldLine::InjectionParameters::emin*MeV2J;
    double emax=SEP::FieldLine::InjectionParameters::emax*MeV2J;

    if (CompressionRatio>SEP::ParticleSource::ShockWave::MaxLimitCompressionRatio) CompressionRatio=SEP::ParticleSource::ShockWave::MaxLimitCompressionRatio;

    double q=3.0*CompressionRatio/(3-1.0);
    double pmin,pmax,speed,pvect[3];
    double mass=PIC::MolecularData::GetMass(spec);

    if (q<1.0) q=1.0;

    pmin=Relativistic::Energy2Momentum(emin,mass);
    pmax=Relativistic::Energy2Momentum(emax,mass);

    double cMin=pow(pmin,-q);

    speed=Relativistic::E2Speed(emin,PIC::MolecularData::GetMass(spec));
    pmin=Relativistic::Speed2Momentum(speed,mass);

    speed=Relativistic::E2Speed(emax,PIC::MolecularData::GetMass(spec));
    pmax=Relativistic::Speed2Momentum(speed,mass);

    double WeightNorm=pow(pmin,1.0-q);

    //to cover the entire range of the particle momentum, the momentum will be generated in the log(p) space 
    //that will be accounted for with a statistical weight correction

    double log_pmin=log(pmin);
    double log_pmax=log(pmax);

    pAbs=pmin*exp(rnd()*(log_pmax-log_pmin));
    WeightCorrection=pow(pAbs,1.0-q)/WeightNorm;

    if (isfinite(WeightCorrection)==false) exit(__LINE__,__FILE__,"Error: nan is found");
  };

  // Apply the physical source independently to every species selected in the
  // build-time AMPS SpeciesList.  Each species has its own mass, base weight,
  // and (in species-global mode) timestep, so its computational Poisson rate
  // and momentum conversion must be calculated inside this loop.
  for (int spec=0;spec<PIC::nTotalSpecies;++spec) {
    const double baseWeight=
        PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
    if (!(baseWeight>0.0) || !isfinite(baseWeight))
      exit(__LINE__,__FILE__,"compiled species has an invalid particle weight");
    double TotalInjectionRate=PhysicalInjectionRate/baseWeight;
    double dtTotal;

    #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    dtTotal=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
    #elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
    dtTotal=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
    #else
    exit(__LINE__,__FILE__,"Error: the global time counter cannot be applied for this case");
    #endif

    if (!(dtTotal>0.0) || !isfinite(dtTotal))
      exit(__LINE__,__FILE__,"compiled species has an invalid particle time step");

    // Limit each compiled species independently.  The individual correction
    // preserves its represented physical rate when the cap is active.
    const int nMaxInjectedParticle=100;
    double WeightCorrectionFactor=1.0;
    if (TotalInjectionRate*dtTotal>nMaxInjectedParticle) {
      const double c=nMaxInjectedParticle/(TotalInjectionRate*dtTotal);
      TotalInjectionRate*=c;
      WeightCorrectionFactor=1.0/c;
      if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
        exit(__LINE__,__FILE__,"Error: set _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_");
      }
    }

    double dtCounter=0.0;
    while ((dtCounter-=log(rnd())/TotalInjectionRate)<dtTotal) {
    double s,x[3],v[3],p;
    int el;
    //1. Generate new particle position 
    el=GetInjectionSurfaceElement(x);
    
    if (Vector3D::DotProduct(x,x)<_SUN__RADIUS_*_SUN__RADIUS_) continue; 

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode=PIC::Mesh::mesh->findTreeNode(x);
    if (startNode->block==NULL) continue;
    if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) continue;

    //2. generate particle velocity  
    double WeightTotalCorrectionFactor;

    s=CompressionRatioTable[el]; 
    GetMomentum_Tenishev2005AIAA(p,WeightTotalCorrectionFactor,s,spec);
    WeightTotalCorrectionFactor*=WeightCorrectionFactor;

    //convert to particle speed 
    const double speed=Relativistic::Momentum2Speed(
        p,PIC::MolecularData::GetMass(spec));
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
    PIC::ParticleBuffer::SetI(spec,newParticleData);

    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightTotalCorrectionFactor,newParticleData);
    #endif

    //apply condition of tracking the particle
    #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    PIC::ParticleTracker::InitParticleID(newParticleData);
    // The tracker contract takes the particle's actual position, velocity,
    // generated AMPS species index, and owning start node.  The former call
    // named nonexistent xInit/xFinal/node variables and failed whenever
    // trajectory tracking was compiled in.
    PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(
        x,v,spec,newParticleData,(void*)startNode);
    #endif


    // Newly injected particles must pass through the same field-line state
    // validation and three-mover dispatch boundary as timestep particles.
    SEP::Mover::DispatchProductionMover(
        newParticle,rnd()*dtTotal,startNode);
    res++;
    }
  }

  return res;
}






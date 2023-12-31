/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"

//the object name and the names of the frames
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Rosetta";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MOON";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="LSO";

int Rosetta::GravityFieldOffset=-1;

/*
int Rosetta::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=-1;
int Rosetta::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=-1;
int Rosetta::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=-1;
*/

static double productionDistributionNASTRAN[30000],cumulativeProductionDistributionNASTRAN[30000];
static double angle;
static double azimuthCenter;
static double zenithCenter;
static cInternalRotationBodyData* Nucleus;

double subSolarPointAzimuth=0.0; //53.0*Pi/180; //0.0;

double DustSizeMin=1.0e-7;
double DustSizeMax=1.0e-2;
double DustTotalMassProductionRate=0.0;
int DustSampleIntervals=10;
double DustSizeDistribution=0.0;

static bool probabilityFunctionDefinedNASTRAN = false;

void Rosetta::Init_BeforeParser() {
#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestDataBuffer);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
#endif

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  //init the dust model                                                                                                                                                                           
  ElectricallyChargedDust::minDustRadius=DustSizeMin; //0.1*_MICROMETER_;
  ElectricallyChargedDust::maxDustRadius=DustSizeMax; //1.0e4*_MICROMETER_;
  ElectricallyChargedDust::Sampling::SetDustSamplingIntervals(DustSampleIntervals);
  ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity=0.01;
  ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity=100.0;
  ElectricallyChargedDust::TotalMassDustProductionRate=DustTotalMassProductionRate;
  ElectricallyChargedDust::SizeDistribution::PowerIndex=DustSizeDistribution;
  ElectricallyChargedDust::Init_BeforeParser();
#endif

}

void Rosetta::Init_AfterParser() {

  //set up the Chamberlen model
  double ExosphereEscapeRate[PIC::nTotalSpecies],ExospehreTemsprature[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) { //ExosphereEscapeRate[spec]=0.0,ExospehreTemsprature[spec]=1000.0;
    ExosphereEscapeRate[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceRate[spec];
    ExospehreTemsprature[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceTemperature[spec];
  }

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  //init the dust model                                                                                                                                                                           
  ElectricallyChargedDust::Init_AfterParser();
#endif
  
  //init Gravity
#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
  InitGravityData();
#endif

  /*
  Exosphere::ChamberlainExosphere::Init(ExosphereEscapeRate,ExospehreTemsprature);



  //set up the model that collected the column integrals at the subsolar point of the limn as a function of the phase angle
  Sampling::SubsolarLimbColumnIntegrals::init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Sampling::SubsolarLimbColumnIntegrals::CollectSample);

  //set up the model that prints the column integrals in the anti-solar direction
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,AntiSolarDirectionColumnMap::Print);

  //set up sampling of velocity distribution functions
  Rosetta::Sampling::VelocityDistribution::Init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Rosetta::Sampling::VelocityDistribution::Sampling,Rosetta::Sampling::VelocityDistribution::OutputSampledData);

  //call init function of the exospheric model
  Exosphere::Init_AfterParser();

  //init the model of calcualting the integrals that correspond to the Kaguya's TVIS observations
  Rosetta::Sampling::Kaguya::Init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Rosetta::Sampling::Kaguya::TVIS::OutputModelData);
  */
}

int Rosetta::RequestDataBuffer(int offset) {
  int TotalDataLength;

  GravityFieldOffset=offset;
  TotalDataLength=3;

  return TotalDataLength*sizeof(double);
}

void Rosetta::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"Gx\",\"Gy\",\"Gz\"");
}

void Rosetta::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;
  double t;

  //Gravity Field
  for (idim=0;idim<3;idim++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+GravityFieldOffset+idim*sizeof(double)));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
}

void Rosetta::InitGravityData(){
  /*  int thread,i,j,k,idim,offset,cnt=0;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  double gravityAccl[3],*position;

  //get coordinated of the center points
  for (node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
	  cell=block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
	  if (cell!=NULL) {  
	    position=cell->GetX();
	    nucleusGravity::gravity(gravityAccl,position);
	    for (idim=0;idim<3;idim++) {
	      *((double*)(cell->GetAssociatedDataBufferPointer()+GravityFieldOffset+idim*sizeof(double)))=gravityAccl[idim];
	      gravityAccl[idim]=0.0;
	    }
	  }
	}
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh->nTotalThreads;thread++) for (node=PIC::Mesh::mesh->DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) {
      block=node->block;

      for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
	for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
	  for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
	    cell=block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
	    if (cell!=NULL) {
	      position=cell->GetX();
	      nucleusGravity::gravity(gravityAccl,position);
	      for (idim=0;idim<3;idim++) {
		*((double*)(cell->GetAssociatedDataBufferPointer()+GravityFieldOffset+idim*sizeof(double)))=gravityAccl[idim];
		gravityAccl[idim]=0.0;
	      }
	    }
	  }
      }
    }

  */
  return ;
}

void Rosetta::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double G[3]={0.0,0.0,0.0};
  int i,idim;
  char *SamplingBuffer;

  for (i=0;i<nInterpolationCoeficients;i++) {

    for (idim=0,SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+GravityFieldOffset;idim<3;idim++) G[idim]+=(*((double*)(SamplingBuffer+idim*sizeof(double))))*InterpolationCoeficients[i];
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+GravityFieldOffset,G,3*sizeof(double));
}

void Rosetta::GetGravityAcceleration(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  register int idim;
  register double *offset=(double*)(GravityFieldOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

  for (idim=0;idim<3;idim++) x[idim]=offset[idim];

}				   

double SodiumStickingProbability(double& ReemissionParticleFraction,double Temp) {
  double res=-1.0;


  //define parametes of the sticking functon
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_           0
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__TEMPERATURE_LIMIT_  1
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__YAKSHINSKY2005SS_   2

#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY__REEMISSION_FRACTION_     1.0
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY__CONSTANT_VALUE_          1.0
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY__TEMPARATURE_LIMIT_       200.0


#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_



  ReemissionParticleFraction=_EXOSPHERE_SODIUM_STICKING_PROBABILITY__REEMISSION_FRACTION_;

#if _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ == _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_
  res= _EXOSPHERE_SODIUM_STICKING_PROBABILITY__CONSTANT_VALUE_;

#elif _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ == _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__TEMPERATURE_LIMIT_
  res= (Temp<_EXOSPHERE_SODIUM_STICKING_PROBABILITY__TEMPARATURE_LIMIT_) ? 1.0 : 0.0;

#elif _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ == _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__YAKSHINSKY2005SS_
  static const double tmin=100.0;
  static const double tmax=500.0;
  static const double dt=5.0;

  static const int nPoints=80;

  struct cDataStruct {
    double t,s;
  };


  static const cDataStruct data[nPoints]={    //digitized from  Yakshinskiy-2005-SS
      {100.00000,0.99830}, {105.00000,0.96768}, {110.00002,0.94512}, {114.99999,0.92417}, {119.99999,0.90322}, {125.00000,0.88065}, {130.00000,0.86293}, {135.00000,0.84681}, {140.00000,0.82586}, {145.00000,0.80813},
      {150.00000,0.78718}, {154.99998,0.76623}, {160.00000,0.75011}, {165.00000,0.73077}, {170.00000,0.71305}, {175.00000,0.69371}, {179.99998,0.67759}, {184.99998,0.65986}, {189.99998,0.64375}, {195.00000,0.62763},
      {200.00000,0.61151}, {205.00000,0.59701}, {210.00002,0.57928}, {214.99998,0.56639}, {220.00000,0.55188}, {225.00000,0.53899}, {229.99998,0.52287}, {235.00000,0.51159}, {240.00000,0.49870}, {245.00000,0.48903},
      {249.99997,0.47936}, {255.00000,0.47130}, {260.00000,0.46325}, {265.00000,0.45035}, {270.00000,0.43907}, {275.00000,0.42618}, {280.00000,0.41651}, {285.00000,0.40523}, {290.00000,0.39717}, {295.00000,0.38911},
      {300.00000,0.38266}, {305.00000,0.36977}, {309.99997,0.35688}, {315.00000,0.34882}, {319.99997,0.33915}, {325.00000,0.33432}, {329.99997,0.32465}, {334.99997,0.31659}, {340.00000,0.31014}, {345.00003,0.30370},
      {350.00000,0.29564}, {354.99997,0.28919}, {360.00000,0.28274}, {365.00000,0.27791}, {370.00003,0.27146}, {375.00000,0.26502}, {379.99997,0.26018}, {384.99997,0.25535}, {390.00000,0.25051}, {395.00003,0.24568},
      {400.00000,0.24084}, {405.00000,0.23762}, {409.99997,0.23278}, {414.99997,0.22956}, {420.00000,0.22473}, {425.00000,0.22150}, {430.00000,0.21828}, {435.00003,0.21506}, {440.00000,0.21022}, {444.99997,0.20700},
      {449.99994,0.20377}, {454.99997,0.20377}, {460.00000,0.20055}, {465.00000,0.19894}, {470.00000,0.19572}, {475.00000,0.19411}, {479.99997,0.19088}, {485.00000,0.19088}, {490.00000,0.18766}, {494.99997,0.18605}
  };

  if (Temp<tmin) res=data[0].s;
  else if (Temp>tmax) res=data[nPoints-1].s;
  else {
    int n;
    double c;

    n=(int)((Temp-tmin)/dt);
    c=(Temp-tmin-dt*n)/dt;

    res=(1.0-c)*data[n].s+c*data[n+1].s;
  }


  return res;
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

  return res;
}


double Exosphere::SurfaceInteraction::StickingProbability(int spec, double& ReemissionParticleFraction,double Temp) {
  double res=0.0;

   switch (spec) {
   case _NA_SPEC_: //case _NAPLUS_SPEC_:
     res=SodiumStickingProbability(ReemissionParticleFraction,Temp);
     break;
   default:
     exit(__LINE__,__FILE__,"the option is not implemented");
   }

   return res;
}


//surface temeprature of the planet
double Exosphere::GetSurfaceTemperature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {
  /*
  //determine if the point on the night side of the Rosetta
  if (CosSubSolarAngle<0.0) return 100.0;

  //determine if the point is within the shadow of the Earth
  if (Rosetta::EarthShadowCheck(x_LOCAL_SO_OBJECT)==true) return 100.0;

  //return the day-side temeprature
  return 280*pow(CosSubSolarAngle,0.25)+100.0;

  const double minTemp[]={172.0,163.0,150.0,145.0,139.0};
  double res,r,zenith,azimuth;
  int angle;

  //determine if the point on the night side of the Rosetta
  if (CosSubSolarAngle<0.0) return minTemp[Rosetta::ndist];
  
  //return the day-side temeprature
  angle=(int) (acos(CosSubSolarAngle)*180.0/Pi);
  if(angle>89) angle=89;
  res=(SurfaceTemp[angle][Rosetta::ndist+1]>minTemp[Rosetta::ndist]) ?  SurfaceTemp[angle][Rosetta::ndist+1] : minTemp[Rosetta::ndist];

  return res;*/
  return 180;
}

//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int nVariables=0;
  /*
  int spec;
  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\",  \"Mean Speed Along the Line of Sight(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec),PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=2;
  }
  
  if (_NA_SPEC_>=0) Rosetta::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=2*_NA_SPEC_;

  //sodium emission
  if (_NA_SPEC_>=0) {
    Rosetta::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=nVariables;
    Rosetta::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=nVariables+1;

    if (vlist!=NULL) sprintf(vlist,"%s,  \"Sodium Emission(5891_58A)\",  \"Sodium Emission(5897_56A)\"",vlist);
    nVariables+=2;
  }
  */
  return nVariables;
  }

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  int spec,cnt=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (res[cnt]>0.0) res[cnt+1]/=res[cnt];
    cnt+=2;
  }
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  /*  int i,j,k,nd,cnt=0,spec;
  double NumberDensity;

  nd=PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node);
  for (i=0;i<resLength;i++) res[i]=0.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //get the local density number
    NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(spec);
    res[cnt++]=NumberDensity;
    res[cnt++]=NumberDensity*node->block->GetCenterNode(nd)->GetMeanParticleSpeed(spec);
  }

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    double BulkVelocity_SO[3],v_LOCAL_SO_FROZEN[3],rHeliocentric,vHeliocentric;

    node->block->GetCenterNode(nd)->GetBulkVelocity(BulkVelocity_SO,spec);

    v_LOCAL_SO_FROZEN[0]=Exosphere::vObject_SO_FROZEN[0]+BulkVelocity_SO[0]+
        Exosphere::RotationVector_SO_FROZEN[1]*x[2]-Exosphere::RotationVector_SO_FROZEN[2]*x[1];

    v_LOCAL_SO_FROZEN[1]=Exosphere::vObject_SO_FROZEN[1]+BulkVelocity_SO[1]-
        Exosphere::RotationVector_SO_FROZEN[0]*x[2]+Exosphere::RotationVector_SO_FROZEN[2]*x[0];

    v_LOCAL_SO_FROZEN[2]=Exosphere::vObject_SO_FROZEN[2]+BulkVelocity_SO[2]+
        Exosphere::RotationVector_SO_FROZEN[0]*x[1]-Exosphere::RotationVector_SO_FROZEN[1]*x[0];

    rHeliocentric=sqrt(pow(x[0]-Exosphere::xObjectRadial,2)+(x[1]*x[1])+(x[2]*x[2]));
    vHeliocentric=(
        (v_LOCAL_SO_FROZEN[0]*(x[0]-Exosphere::xObjectRadial))+
        (v_LOCAL_SO_FROZEN[1]*x[1])+(v_LOCAL_SO_FROZEN[2]*x[2]))/rHeliocentric;


    //brightness of the exospheric sodium
    if (spec==_NA_SPEC_) {
      //check if the point is outside of the Rosetta's and Earth's shadows
      if ( (Rosetta::EarthShadowCheck(x)==false) && ((x[0]>0.0)||(x[1]*x[1]+x[2]*x[2]>_RADIUS_(_MOON_)*_RADIUS_(_MOON_))) ) {
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5891_58A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5897_56A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
      }
      else res[cnt++]=0.0,res[cnt++]=0.0;
    }
  }


  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
*/}


//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble EphemerisTime) {
  double res=0.0;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble State[6],ltlocal;
  double EccentricityVector[3];
  double Speed2,a,c,absEccentricity;
  const double GravitationalParameter=GravityConstant*_MASS_(_EARTH_);
  double vRosetta[3],xRosetta[3],rGeocentric=0.0;
  int idim;


  spkezr_c("Rosetta",EphemerisTime,"MSGR_HCI","none","Earth",State,&ltlocal);

  for (idim=0;idim<3;idim++) {
    xRosetta[idim]=State[idim]*1.0E3;
    vRosetta[idim]=State[idim+3]*1.0E3;

    rGeocentric+=pow(xRosetta[idim],2);
  }

  rGeocentric=sqrt(rGeocentric);
  Speed2=vRosetta[0]*vRosetta[0]+vRosetta[1]*vRosetta[1]+vRosetta[2]*vRosetta[2];
  c=xRosetta[0]*vRosetta[0]+xRosetta[1]*vRosetta[1]+xRosetta[2]*vRosetta[2];

  for (idim=0,absEccentricity=0.0,a=0.0;idim<3;idim++) {
    EccentricityVector[idim]=Speed2/GravitationalParameter*xRosetta[idim] - c/GravitationalParameter*vRosetta[idim] - xRosetta[idim]/rGeocentric;
    absEccentricity+=EccentricityVector[idim]*EccentricityVector[idim];
    a+=EccentricityVector[idim]*xRosetta[idim];
  }

  absEccentricity=sqrt(absEccentricity);
  res=acos(a/(absEccentricity*rGeocentric));

  if (c<0.0) res=2.0*Pi-res;
#endif

  return res;
}


long int Rosetta::InjectionBoundaryModel_Limited() {
  int spec;
  long int res=0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=InjectionBoundaryModel_Limited(spec);

  return res;
}

long int Rosetta::InjectionBoundaryModel_Limited(int spec) {
  cInternalSphericalData *Sphere;
  double ModelParticlesInjectionRate,ParticleWeight,LocalTimeStep,TimeCounter=0.0,x_SO_OBJECT[3],x_IAU_OBJECT[3],v_SO_OBJECT[3],v_IAU_OBJECT[3],*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  double ParticleWeightCorrection=1.0;
  bool flag=false;
  int SourceProcessID;

  double totalProductionRate=Rosetta::GetTotalProductionRateBjornNASTRAN(spec,Sphere);

  const int nMaxInjectedParticles=10*PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber;


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[spec];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  ModelParticlesInjectionRate=totalProductionRate/ParticleWeight;

  //  if (ModelParticlesInjectionRate*ParticleWeight*LocalTimeStep<1.0E-10) return 0;

  if (ModelParticlesInjectionRate*LocalTimeStep>nMaxInjectedParticles) {
    ParticleWeightCorrection=ModelParticlesInjectionRate*LocalTimeStep/nMaxInjectedParticles;
    ModelParticlesInjectionRate/=ParticleWeightCorrection;
  }

  //definition of indexes TEMPORARY!!!!!
  //  int _EXOSPHERE__SOURCE_MAX_ID_VALUE_=0;
  int _EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_=0;

  //calcualte probabilities of each source processes                                                                 
  double TotalFlux,FluxSourceProcess[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_]; //,ProbabilitySourceProcess[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];                                                                                               
int iSource;

for (iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) FluxSourceProcess[iSource]=0.0; //,ProbabilitySourceProcess[iSource]=0.0;                                                                                         

TotalFlux=totalProductionRate;

//only Used defined source here since we only want the Bjorn model so far
//calculate the source rate due to user defined source functions                                                   
FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]=Rosetta::GetTotalProductionRateBjornNASTRAN(spec,Sphere);

 double CalculatedSourceRate[PIC::nTotalSpecies][1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];
 CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]=0.0;


#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
 if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
   static double GrainInjectedMass=0.0;
   PIC::Mesh::cDataBlockAMR *block;
   double GrainRadius,GrainMass,GrainWeightCorrection;
   int GrainVelocityGroup;

   GrainInjectedMass+=ElectricallyChargedDust::TotalMassDustProductionRate*LocalTimeStep;

   while (GrainInjectedMass>0.0) {
     startNode=NULL;

     //generate a particle                                                                                             
     PIC::ParticleBuffer::byte tempParticleData[PIC::ParticleBuffer::ParticleDataLength];

     for (int ii=0;ii<PIC::ParticleBuffer::ParticleDataLength;ii++) tempParticleData[ii]=0;
     
     PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

     flag=Rosetta::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,tempParticleData);
     ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
     GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
     GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
     SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
     if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;

     if (flag==false) continue;
     if ((block=startNode->block)->GetLocalTimeStep(_DUST_SPEC_)/LocalTimeStep<rnd()) continue;

     //determine the velocity group of the injected grain;                                                                                                                                                                                                                     
     //calculate additional particle weight correction because the particle will be placed in a different weight group                                                                                                                                                         
     GrainVelocityGroup=ElectricallyChargedDust::GrainVelocityGroup::GetGroupNumber(v_SO_OBJECT);
     GrainWeightCorrection*=block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup)/block->GetLocalTimeStep(_DUST_SPEC_);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
     GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+GrainVelocityGroup];
#else
     exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

    //determine the surface element of the particle origin                                                            
     PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

     PIC::ParticleBuffer::SetX(x_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
     PIC::ParticleBuffer::SetV(v_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
     PIC::ParticleBuffer::SetI(_DUST_SPEC_+GrainVelocityGroup,(PIC::ParticleBuffer::byte*)tempParticleData);
     
     ElectricallyChargedDust::SetGrainCharge(0.0,(PIC::ParticleBuffer::byte*)tempParticleData);
     ElectricallyChargedDust::SetGrainMass(GrainMass,(PIC::ParticleBuffer::byte*)tempParticleData);
     ElectricallyChargedDust::SetGrainRadius(GrainRadius,(PIC::ParticleBuffer::byte*)tempParticleData);
     
     PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
     
     
     newParticle=PIC::ParticleBuffer::GetNewParticle();
     newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
//     memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);
     
     PIC::ParticleBuffer::CloneParticle(newParticleData,tempParticleData);
     nInjectedParticles++;
     
     //inject the particle into the system                                                                             
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
     if (startNode==NULL) exit(__LINE__,__FILE__,"Error: the node is not defined");
     if ((startNode->Thread!=PIC::ThisThread)||(startNode->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is n\
ot defined");
#endif
     
     _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode,true);
   }
   
 }else{
#endif
   while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
  //determine the source process to generate a particle's properties                                               
  do {
    SourceProcessID=(int)(rnd()*(1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_));
  }
  while (FluxSourceProcess[SourceProcessID]/TotalFlux<rnd());
  
  //generate a particle                                                                                             
  PIC::ParticleBuffer::byte tempParticleData[PIC::ParticleBuffer::ParticleDataLength];

  for (int ii=0;ii<PIC::ParticleBuffer::ParticleDataLength;ii++) tempParticleData[ii]=0;

  PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);
  
  //to satisfy the compiler and fit the while structure                                                             
  if (false) {}
  
  //Add the user defined particle gineration                                                                        
  flag=Rosetta::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,tempParticleData);
  SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
  if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
  
  else {
    continue;
  }
  
  if (flag==false) continue;
  
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  if (startNode->block->GetLocalTimeStep(spec)/LocalTimeStep<rnd()) continue;
#endif
 
  //determine the surface element of the particle origin                                                            
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  PIC::ParticleBuffer::SetX(x_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetV(v_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);


  newParticle=PIC::ParticleBuffer::GetNewParticle();
  newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
//  memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

  PIC::ParticleBuffer::CloneParticle(newParticleData,tempParticleData);
  nInjectedParticles++;

  //inject the particle into the system                                                                             
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (startNode==NULL) exit(__LINE__,__FILE__,"Error: the node is not defined");
  if ((startNode->Thread!=PIC::ThisThread)||(startNode->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is n\
ot defined");
#endif
  
  _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode,true);
   }
   
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
 }
#endif
 
 return nInjectedParticles;
}

double Rosetta::GetTotalProductionRateBjornNASTRAN(int spec, cInternalSphericalData* Sphere){
  return 2.0e15;
}


bool Rosetta::GenerateParticlePropertiesBjornNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere,char* tempParticleData) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  const double NightSideProduction[5]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0};
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=3.3*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i;
  double rSphere=1980.0;
  double area = 0.0;
  double totalArea = 0.0;
  int faceAttribute = 666;
  
  if (probabilityFunctionDefinedNASTRAN==false) {
    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;
    
    for (i=0;i<totalSurfaceElementsNumber;i++)
    {
      faceAttribute = CutCell::BoundaryTriangleFaces[i].attribute;
      if ((faceAttribute == 1) || (faceAttribute == 0)){
        area = CutCell::BoundaryTriangleFaces[i].SurfaceArea;
        totalArea = totalArea + area;
        productionDistributionNASTRAN[i] = area;
      }else{
        productionDistributionNASTRAN[i] = 0.0;
      }
    }
    cout << "total Area: " << totalArea << endl;
    for (i=0; i<totalSurfaceElementsNumber; i++)
    {
      productionDistributionNASTRAN[i] /= totalArea;
    }
    
    cumulativeProductionDistributionNASTRAN[0]=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      if (i==0) {
        cumulativeProductionDistributionNASTRAN[i]+=productionDistributionNASTRAN[i];
      }else{
        cumulativeProductionDistributionNASTRAN[i]=cumulativeProductionDistributionNASTRAN[i-1]+productionDistributionNASTRAN[i];
      }
      
    }
    probabilityFunctionDefinedNASTRAN=true;
  }
  
  //Computation of the segment where the particle will be created
  gamma=rnd();
  i=0;
  while (gamma>cumulativeProductionDistributionNASTRAN[i]){
    i++;
  }
    
  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
  //  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,PIC::Mesh::mesh->EPS);
  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,1e-4);
  //  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,1.0e-1);
  for (idim=0;idim<3;idim++) ExternalNormal[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
  
  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
    c+=ExternalNormal[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
    X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
  }
  cosSubSolarAngle=c/sqrt(X);
  
  for (idim=0;idim<3;idim++) ExternalNormal[idim]*=-1;

  //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
  x_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*x_LOCAL_IAU_OBJECT[2]);
    

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh->findTreeNode(x_LOCAL_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) return false;
  
  //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};

  SurfaceTemperature=GetSurfaceTemperature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-10.0*ExternalNormal[idim];
  else PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#endif
  
  //init the internal degrees of freedom if needed
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  PIC::IDF::InitRotTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
  PIC::IDF::InitVibTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif

  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  v_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_LOCAL_IAU_OBJECT[2]);
  
  memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
  memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));
  
  return true;
}

void Rosetta::AntiSolarDirectionColumnMap::Print(int DataOutputFileNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  FILE *fout=NULL;
  char fname[300];
  SpiceDouble xform[6][6],EarthState[6],lt;
  double xEarthLSO[3];
  SpiceDouble lGSE[6]={0,0,0,0,0,0},lLSO[6]={0,0,0,0,0,0};  //only 3 first components of the vectors are used. all 6 components are needed in the definition in order SPICE routines work correctly
  int idim;

  //calculate the currect position of the Earth
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,"LSO","none","Rosetta",EarthState,&lt);
  for (idim=0;idim<3;idim++) xEarthLSO[idim]=1.0E3*EarthState[idim];

  //calculate the rotation matrix from 'GSE' to 'LSO'
  sxform_c("GSE","LSO",Exosphere::OrbitalMotion::et,xform);

  //determine the number of angular points
  int nZenithPoints;
  double dZ,rr,ZenithAngle,AzimuthAngle,dZenithAngle;

  dZ=dZenithAngleMin;
  rr=(maxZenithAngle+dZenithAngleMax)/(maxZenithAngle+dZenithAngleMin);
  nZenithPoints=(long int)(log(dZenithAngleMax/dZenithAngleMin)/log(rr)-2.0);
  rr=pow(dZenithAngleMax/dZenithAngleMin,1.0/(nZenithPoints+2.0));

  nZenithPoints=0,ZenithAngle=dZenithAngleMin,dZenithAngle=dZenithAngleMin;

  while (ZenithAngle<maxZenithAngle) {
    ZenithAngle+=dZenithAngle;
    dZenithAngle*=rr;
    nZenithPoints++;
  }

  if (PIC::ThisThread==0) {
    const SpiceInt lenout = 35;
    SpiceChar utcstr[lenout+2];
    char vlist[_MAX_STRING_LENGTH_PIC_]="";

    //open data file
    sprintf(fname,"%s/pic.Rosetta.Anti-sunwardColumnIntegrals.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);

    fout=fopen(fname,"w");

    et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
    fprintf(fout,"TITLE=\"UTC=%s\"\n",utcstr);

//    fprintf(fout,"VARIABLES=\"Angle from the anti-solar direction [degree]\", \"Angle Out of Ecpliptic Plane [degree]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");

    ColumnIntegral::GetVariableList(vlist);
    fprintf(fout,"VARIABLES=\"l[0]\", \"l[1]\" %s \n",vlist);


    fprintf(fout,"ZONE T=\"Column Density Map\"\n");
    fprintf(fout,"I=%i, J=%i, K=1, ZONETYPE=Ordered\n",nAzimuthPoints+1,nZenithPoints);
    fprintf(fout,"DATAPACKING=POINT\n");
//    fprintf(fout,"DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n");
  }


  //calculate the integrals
  //calcualte the column integrals
  int StateVectorLength=ColumnIntegral::GetVariableList(NULL);
  double StateVector[StateVectorLength];


  nZenithPoints=0,ZenithAngle=dZenithAngleMin,dZenithAngle=dZenithAngleMin;

  while (ZenithAngle<maxZenithAngle) {

  for (int iAzimuthPoint=0;iAzimuthPoint<nAzimuthPoints+1;iAzimuthPoint++) {
    AzimuthAngle=2.0*Pi*double(iAzimuthPoint)/double(nAzimuthPoints);

      //get the pointing vector of integration in 'GSE' frame
      lGSE[0]=-cos(ZenithAngle);
      lGSE[1]=sin(ZenithAngle)*sin(AzimuthAngle);
      lGSE[2]=sin(ZenithAngle)*cos(AzimuthAngle);

      //convert the pointing vector from 'GSE' to 'LSO'
      mxvg_c(xform,lGSE,6,6,lLSO);

      //get the integrals
      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,xEarthLSO,lLSO,ColumnIntegral::CoulumnDensityIntegrant);
      ColumnIntegral::ProcessColumnIntegrationVector(StateVector,StateVectorLength);

//      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",ZenithAngle/Pi*180.0,AzimuthAngle/Pi*180.0,StateVector[0],StateVector[1],StateVector[2]);
//      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",lGSE[1]/sqrt(lGSE[1]*lGSE[1]+lGSE[2]*lGSE[2]),lGSE[2]/sqrt(lGSE[1]*lGSE[1]+lGSE[2]*lGSE[2]),StateVector[0],StateVector[1],StateVector[2]);


      if (PIC::ThisThread==0) {
        fprintf(fout,"%e   %e",lGSE[1],lGSE[2]);

        for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",StateVector[i]);
        fprintf(fout,"\n");
      }

  }

    ZenithAngle+=dZenithAngle;
    dZenithAngle*=rr;
    nZenithPoints++;
  }


  if (PIC::ThisThread==0) fclose(fout);
#endif
}

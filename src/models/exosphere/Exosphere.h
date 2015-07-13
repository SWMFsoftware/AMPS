//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury.h
 *
 *  Created on: Feb 10, 2012
 *      Author: vtenishe
 */

//$Id$


#ifndef _EXOSPHERE_
#define _EXOSPHERE_



#include "SingleVariableDistribution.h"
#include "SingleVariableDiscreteDistribution.h"
#include "constants.h"

#include "Exosphere.dfn"






//user defined settings of the exospheric model
//#include "UserDefinition.Exosphere.h"
#include "Na.h"

//define the symbolic id of source processes
static const char _EXOSPHERE__SOURCE_SYMBOLIC_ID_[][100]={"ExternalBoundaryInjection","ImpactVaposization","PhotonStimulatedDesorption","ThermalDesorption","SolarWindSputtering"};

//the default value for for the list of the SPICE kernels that will be loaded
static const int nFurnishedSPICEkernels=0;
static const char SPICE_Kernels[][_MAX_STRING_LENGTH_PIC_]={""};

//the default values of the list of the referenced ground based observations
static const int nReferenceGroundBasedObservations=0;
static const char ReferenceGroundBasedObservationTime[][_MAX_STRING_LENGTH_PIC_]={""};

//the default location of the SPICE kernels
const char SPICE_Kernels_PATH[_MAX_STRING_LENGTH_PIC_]="";

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif


namespace Exosphere {

  //the name os the simulated object and related coordinate frames
  extern char ObjectName[_MAX_STRING_LENGTH_PIC_],IAU_FRAME[_MAX_STRING_LENGTH_PIC_],SO_FRAME[_MAX_STRING_LENGTH_PIC_];


  //simulation date and position of Object at the time of simulations
  //const char SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_]="2011-04-13T00:00:00"; //"2001-03-01T00:00:00";  ////"2011-01-01T00:00:00";
  //extern char SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_];

  static const char SimulationStartTimeString[]="2000-01-01T00:00:00";

  extern double xObject_HCI[3],vObject_HCI[3],xEarth_HCI[3],vEarth_HCI[3],xEarth_SO[3],vEarth_SO[3],xSun_SO[3],vSun_SO[3];
  extern double vObjectRadial,xObjectRadial;

  extern double vObject_SO_FROZEN[3]; //Velocity of Object in inertial frame, which axis coinsides with that of SO at a given time
  extern double RotationVector_SO_FROZEN[3],RotationRate_SO_FROZEN; //the vectors of the direction of rotation and the rotation rate of the SO in SO_FROZEN


  //typical solar wind conditions far from the planet
  static const double swVelocity_Typical[]={0.0,0.0,0.0};
  static const double swB_Typical[]={0.0,0.0,0.0};
  static const double swTemperature_Typical=0.0;
  static const double swNumberDensity_Typical=0.0;
  extern double swE_Typical[3];

  //the total number of source processes
  extern int nTotalSourceProcesses;

  //the array of flags that defines wether the source process change the surface aboundance of the volatile
  static const bool Source_DeplitSurfaceSpeciesAbundance_Flag[]={true};


  //the sphere that representd the planet
  extern cInternalSphericalData *Planet;

  //init the model
  void Init_BeforeParser();
  void Init_AfterParser();
  void Init_SPICE();
 // void Init_AfterMesh();

  //ICES data preprocessor -> set up typical values of the solar wind in the regions where the SWMF values have not been found
  void SWMFdataPreProcessor(double *x,PIC::CPLR::DATAFILE::ICES::cDataNodeSWMF& data);


  //make coulumn integration
  namespace ColumnIntegral {
    void CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
    int GetVariableList(char *vlist=NULL);
    void ProcessColumnIntegrationVector(double *res,int resLength);



    void GetSubsolarPointDirection(double *LimbDirection,double *EarthPosition);
    void Tail(char *name);
    void Limb(char *name);
//    void Map(char *fname,double dXmax,double dZmax,int nXpoints);
    void CircularMap(char *fname,double rmax,double dRmin,double dRmax,int nAzimuthPoints,SpiceDouble EphemerisTime);
  }

  //Chamberlain Exosphere Model (Shen-1963-JAS,Valeille-2009)
  namespace ChamberlainExosphere {
    extern double *SpecieExobaseEscapeRate,*SpecieExobaseTemperature;
    extern bool ModelInitFlag;

    void Init(double *ExobaseEscapeRate,double *ExobaseTemperature);

    void PrintVariableList(FILE* fout);
    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  }



  //Sampling
  namespace Sampling {

    //reference ground based observations
    namespace ReferenceGroundBasedObservations {

       struct cObservationTag {
         char TimeStamp[_MAX_STRING_LENGTH_PIC_];
         double TAA,PhaseAngle;
         SpiceDouble et;
         int nOutputFile;
       };

       extern cObservationTag RemoteObservationList[nReferenceGroundBasedObservations];
       void init();
       void OutputSampledData(SpiceDouble etStartInterval,SpiceDouble etFinishInterval,int nObjectOutputFile);
    }

    //sample the source rates
    extern double CalculatedSourceRate[PIC::nTotalSpecies][1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];

    //offsets for sampling densities that are due to different sampling processes
    extern int SamplingDensityOffset[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];
    extern int CellSamplingDataOffset;

    //the field in the particle's data that keeps the id of the source process due to which the particle has beed produced
    extern long int ParticleData_SourceProcessID_Offset;
    extern long int ParticleData_OriginSurfaceElementNumber_Offset;

    //sample the planet's night side return flux
    extern double **PlanetNightSideReturnFlux;

    //total return flux; the rate of sticking to the planet's surface
    extern double *TotalPlanetReturnFlux,*PlanetSurfaceStickingRate;

    namespace OutputDataFile {
      //matrix for transformation SO->HCI coordinate frame (used only when prepare data files)
      extern SpiceDouble SO_to_HCI_TransformationMartix[6][6];

      void PrintVariableList(FILE* fout,int DataSetNumber);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
    }

    namespace OutputSurfaceDataFile {
      void PrintVariableList(FILE* fout);
      void PrintTitle(FILE* fout);
      void PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
      void flushCollectingSamplingBuffer(cInternalSphericalData* Sphere);
    }

    //clean the model sampling buffers after output of the data file
    inline void FlushSamplingDataBuffers() {OutputSurfaceDataFile::flushCollectingSamplingBuffer(Planet);}

    //request the sampling buffer
    int RequestSamplingData(int offset);

    //set and get the source id
    inline int GetParticleSourceID(PIC::ParticleBuffer::byte *pData) {return *(int*)(pData+ParticleData_SourceProcessID_Offset);}
    inline void SetParticleSourceID(int id,PIC::ParticleBuffer::byte *pData) {*(int*)(pData+ParticleData_SourceProcessID_Offset)=id;}

    //set and get the surface elemetn number of the particle's origine
    inline int GetParicleOriginSurfaceElementNumber(PIC::ParticleBuffer::byte *pData) {return *(int*)(pData+ParticleData_OriginSurfaceElementNumber_Offset);}
    inline void SetParicleOriginSurfaceElementNumber(int el,PIC::ParticleBuffer::byte *pData) {*(int*)(pData+ParticleData_OriginSurfaceElementNumber_Offset)=el;}

    //sample particle's data
    void SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec);
    void SampleModelData();
    void OutputSampledModelData(int DataOutputFileNumber);

    typedef void (*fUserDefinedAdditionalData_VariableList_OutputSampledModelData)(FILE*);
    typedef void (*fUserDefinedAdditionalData_OutputSampledModelData)(FILE*,int);
    extern fUserDefinedAdditionalData_VariableList_OutputSampledModelData UserDefinedAdditionalData_VariableList_OutputSampledModelData;
    extern fUserDefinedAdditionalData_OutputSampledModelData UserDefinedAdditionalData_OutputSampledModelData;

    void inline SetUserDefinedAdditionalOutputSampledModelDataFunctions(fUserDefinedAdditionalData_VariableList_OutputSampledModelData t0,fUserDefinedAdditionalData_OutputSampledModelData t1) {
      UserDefinedAdditionalData_VariableList_OutputSampledModelData=t0;
      UserDefinedAdditionalData_OutputSampledModelData=t1;
    }
  }


  //orbital motion of the Planet
  namespace OrbitalMotion {
    extern double AccumulatedPlanetRotation,TotalSimulationTime,TAA;

    //SPICE ephemeris time
    extern SpiceDouble et,lt;

    //direction to the Sun and the angle of the rotaio berween planetary axises and the direction to the Sun on the Z-plane
    extern double SunDirection_IAU_OBJECT[3];

    //matrixes for tranformation SO->IAU and IAU->SO coordinate frames
    extern SpiceDouble SO_to_IAU_TransformationMartix[6][6],IAU_to_SO_TransformationMartix[6][6];


    //parameters of orbital motion of Object
    extern double CoordinateFrameRotationRate;
    extern double RotationRateVector_SO_J2000[3];

    //the number of sub-intervals printed for a single interval between outputs of model data files in the pic.OrbitalData.dat
    extern int nOrbitalPositionOutputMultiplier;

    inline double GetCosineSubsolarAngle(double *x) { return x[0]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);}

    //update the Transformation Matrixes and calculate the rotation vector
    void UpdateTransformationMartix();
    void UpdateRotationVector_SO_J2000();

    //calculate TAA
    double GetTAA(SpiceDouble);
    double GetTAA(const char*);
    double GetTAA(const char* TargetName, const char* CenterBodyName, double CenterBodyMass, SpiceDouble EphemerisTime);

    //calculate the phase angle
    double GetPhaseAngle(SpiceDouble);
    double GetPhaseAngle(const char*);

    namespace FrameRotation {
       void GetRotationAxis(double *RotationAxis,double &RotationAngle,const char *FrameName,double etStartRotation,double etFinishRotation);
       double GetRotationVector(double *RotationVector,const char *FrameName,double etStartRotation,double etFinishRotation);
       void GetRotationMatrix(double RotationMatrix[3][3],double *RotationAxis,double RotationAngle);
    }
  }

  //surface temeprature of the planet
  double GetSurfaceTemeprature(double CosSubSolarAngle,double *x_SO);



  //sources
  namespace SourceProcesses {

    //init the model
    void Init();

    double GetInjectionEnergy(int spec,int SourceProcessID);

    //generate particle position and velocity
  //generate particle properties
  inline bool GenerateParticleProperties(int spec,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere,cSingleVariableDiscreteDistribution<int> *SurfaceInjectionDistribution,cSingleVariableDistribution<int> *EnergyDistribution,int SourceProcessID) {
    double ExternalNormal[3];

    //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
    double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
    int nZenithElement,nAzimuthalElement;
    unsigned int idim,el;

    el=SurfaceInjectionDistribution->DistributeVariable();
    Exosphere::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
    Exosphere::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

    x_LOCAL_IAU_OBJECT[0]=sphereRadius*ExternalNormal[0];
    x_LOCAL_IAU_OBJECT[1]=sphereRadius*ExternalNormal[1];
    x_LOCAL_IAU_OBJECT[2]=sphereRadius*ExternalNormal[2];

    /*
    ExternalNormal[0]*=-1.0;
    ExternalNormal[1]*=-1.0;
    ExternalNormal[2]*=-1.0;
    */

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
    startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
    if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

    //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
    double c=0.0,rVel=0.0,lVel[3];
#if _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION_ == _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION__NUMERIC_
    double Speed=sqrt(EnergyDistribution->DistributeVariable()*2.0/PIC::MolecularData::GetMass(spec));
#elif _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION_ == _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION__USER_DEFINED_
    double Speed=sqrt(GetInjectionEnergy(spec,SourceProcessID)*2.0/PIC::MolecularData::GetMass(spec));
#else
    exit(__LINE__, __FILE__, "ERROR: _ENERGY_DISTRIBUTION_INVERSION_ is not defined ")
#endif

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    if (isfinite(Speed)==false) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

#if _EXOSPHERE__INJECTION_ANGLE_DISTRIBUTION_ == _EXOSPHERE__INJECTION_ANGLE_DISTRIBUTION__UNIFORM_
    for (idim=0;idim<3;idim++) {
      lVel[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
      rVel+=pow(lVel[idim],2);

      c+=ExternalNormal[idim]*lVel[idim];
    }

    rVel=Speed/sqrt(rVel);

    if (c>0.0) {
      //the distributed velocity vector is directed into the domain
      for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=lVel[idim]*rVel;
    }
    else {
      //the distributed velocity vector is directed into the planet -> redirect it
      for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=(lVel[idim]-2.0*c*ExternalNormal[idim])*rVel;
    }
#elif _EXOSPHERE__INJECTION_ANGLE_DISTRIBUTION_ ==  _EXOSPHERE__INJECTION_ANGLE_DISTRIBUTION__KNUDSEN_COSINE_
    // first, distribute velocity assuming normal vector to be {0,0,1}
    double Phi      = PiTimes2*rnd();
    double SinTheta = pow(rnd(),0.5);
    lVel[0] = cos(Phi) * SinTheta;
    lVel[1] = sin(Phi) * SinTheta;
    lVel[2] = pow(1.0 - SinTheta*SinTheta, 0.5);
    // actual normal is obtained from {0,0,1} by 2 rotations:
    //   1: around initial z axis by angle A
    //   2: around new     y axis by angle B
    //
    // rotation matrix:  / CosB*CosA -SinA SinB*CosA \
    //                  |  CosB*SinA  CosA SinB*SinA  |
    //                   \-SinB         0  CosB      /
    //
    //         => normal = {SinB*CosA, SinB*SinA, CosB}
    double CosB = ExternalNormal[2];
    if(CosB < 1.0 && CosB > -1.0) {
      double SinB = pow(1.0 - CosB*CosB, 0.5);
      double CosA, SinA;
      CosA = ExternalNormal[0] / SinB;
      SinA = ExternalNormal[1] / SinB;
      v_LOCAL_IAU_OBJECT[0]=Speed*( lVel[0]*CosB*CosA - lVel[1]*SinA + lVel[2]*SinB*CosA);
      v_LOCAL_IAU_OBJECT[1]=Speed*( lVel[0]*CosB*SinA + lVel[1]*CosA + lVel[2]*SinB*SinA);
      v_LOCAL_IAU_OBJECT[2]=Speed*(-lVel[0]*SinB                     + lVel[2]*CosB     );
    }
    else    // if abs(CosB)==1 no rotation is needed
      for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=lVel[idim]*Speed;
    
#else
    exit(__LINE__, __FILE__, "ERROR: _EXOSPHERE__INJECTION_ANGLE_DISTRIBUTION_ is not defined")
#endif

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



/*    {  //test: set vertical injection velocity
      double l=0.0;
      int i;

      for (i=0;i<3;i++) l+=pow(x_LOCAL_IAU_OBJECT[i],2);

      l=sqrt(l);

      for (i=0;i<3;i++) {
        x_LOCAL_SO_OBJECT[i]=x_LOCAL_IAU_OBJECT[i];

        v_LOCAL_IAU_OBJECT[i]=1500.0*x_LOCAL_IAU_OBJECT[i]/l;
        v_LOCAL_SO_OBJECT[i]=v_LOCAL_IAU_OBJECT[i];
      }

    }*/



    memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
    memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
    memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
    memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      for (int idim=0;idim<3;idim++) {
        if ((isfinite(x_SO_OBJECT[idim])==false) || (isfinite(x_IAU_OBJECT[idim])==false) || (isfinite(v_SO_OBJECT[idim])==false) || (isfinite(v_IAU_OBJECT[idim])==false)) {
          exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
        }
      }
#endif
#endif

    return true;
  }

    namespace ThermalDesorption {
      static const double ThermalDesorption_uThermal[]={0.0};
      static const double ThermalDesorption_VibrationalFrequency[]={0.0};

      extern double SourceRate[PIC::nTotalSpecies],maxLocalSourceRate[PIC::nTotalSpecies];
//      extern double CalculatedTotalSodiumSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution[PIC::nTotalSpecies];
      double GetSurfaceElementProductionRate(int nElement,int *spec);

      inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return SourceRate[spec];}

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],temp,SodiumSurfaceElementPopulation;

        SodiumSurfaceElementPopulation=((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[spec][SurfaceElement];
        if (SodiumSurfaceElementPopulation<0.0) return 0.0;

        memcpy(sun,Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT,3*sizeof(double));
        memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];

        //calculate local temeprature
        //get the local position of the surface element in SO_FRAME
        double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3];

        x_LOCAL_IAU_OBJECT[0]=_RADIUS_(_TARGET_)*norm[0];
        x_LOCAL_IAU_OBJECT[1]=_RADIUS_(_TARGET_)*norm[1];
        x_LOCAL_IAU_OBJECT[2]=_RADIUS_(_TARGET_)*norm[2];

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


        temp=GetSurfaceTemeprature(res,x_LOCAL_SO_OBJECT);
        res=ThermalDesorption_VibrationalFrequency[spec]*SodiumSurfaceElementPopulation*exp(-ThermalDesorption_uThermal[spec]/(Kbol*temp));

        return res;
      }


      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double* v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
        double ExternalNormal[3],CosSubSolarAngle;

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
#endif


        //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
        double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
        int nZenithElement,nAzimuthalElement;
        unsigned int el;

        el=SurfaceInjectionDistribution[spec].DistributeVariable();
        Exosphere::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
        Exosphere::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

        CosSubSolarAngle=(Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[0]*ExternalNormal[0])+
            (Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[1]*ExternalNormal[1])+
            (Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[2]*ExternalNormal[2]);

        x_LOCAL_IAU_OBJECT[0]=sphereRadius*ExternalNormal[0];
        x_LOCAL_IAU_OBJECT[1]=sphereRadius*ExternalNormal[1];
        x_LOCAL_IAU_OBJECT[2]=sphereRadius*ExternalNormal[2];

        ExternalNormal[0]*=-1.0;
        ExternalNormal[1]*=-1.0;
        ExternalNormal[2]*=-1.0;

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
        startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
        if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

        //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
        double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};

        SurfaceTemperature=GetSurfaceTemeprature(CosSubSolarAngle,x_LOCAL_SO_OBJECT);
        PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);


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



    }

    namespace PhotonStimulatedDesorption {
      static const double PhotonStimulatedDesorption_PhotonFlux_1AU=0.0;
      static const double PhotonStimulatedDesorption_CrossSection[]={0.0};

      static const double PhotonStimulatedDesorption_minInjectionEnergy[]={0.0};
      static const double PhotonStimulatedDesorption_maxInjectionEnergy[]={0.0};

      extern double SourceRate[PIC::nTotalSpecies],maxLocalSourceRate[PIC::nTotalSpecies];

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution[PIC::nTotalSpecies];
      double GetSurfaceElementProductionRate(int nElement,int *spec);

      //evaluate nemerically the source rate
//      extern double CalculatedTotalSodiumSourceRate;

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],SurfaceElementPopulation;

        SurfaceElementPopulation=((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[spec][SurfaceElement];
        if (SurfaceElementPopulation<0.0) return 0.0;

        memcpy(sun,Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT,3*sizeof(double));
        memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];


        res*=(res>0.0) ? SurfaceElementPopulation*PhotonStimulatedDesorption_CrossSection[spec]*PhotonStimulatedDesorption_PhotonFlux_1AU*pow(_AU_/Exosphere::xObjectRadial,2.0) : 0.0;
        return res;
      }

      inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return SourceRate[spec];}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution[PIC::nTotalSpecies];
      double EnergyDistributionFunction(double e,int *spec);

      //generate particle properties
      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
#endif

        if (BoundaryElementType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_) exit(__LINE__,__FILE__,"Error: particle ejection from a non-spehtical body is not implemented");

        cInternalSphericalData* Sphere=(cInternalSphericalData*)BoundaryElement;

        return Exosphere::SourceProcesses::GenerateParticleProperties(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,SurfaceInjectionDistribution+spec,EnergyDistribution+spec,_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_);
      }
    }


    namespace VerticalInjection {
      static const double VerticalInjection_SourceRate[]={0.0};
      static const double VerticalInjection_BulkVelocity[]={0.0};

      double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer);

      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
        int idim;
        double ExternalNormal[3],r=0.0;

        //generate the new particle position and velocity
        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
          r+=pow(ExternalNormal[idim],2);
        }

        r=sqrt(r);

        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]/=r;
          x_IAU_OBJECT[idim]=sphereRadius*ExternalNormal[idim]+sphereX0[idim];
          x_SO_OBJECT[idim]=x_IAU_OBJECT[idim];

          v_SO_OBJECT[idim]=VerticalInjection_BulkVelocity[spec]*ExternalNormal[idim];
          v_IAU_OBJECT[idim]=v_SO_OBJECT[idim];
        }

        //determine if the particle belongs to this processor
        startNode=PIC::Mesh::mesh.findTreeNode(x_SO_OBJECT,startNode);
        return (startNode->Thread==PIC::Mesh::mesh.ThisThread) ? true : false;
      }


    }

    namespace BackgroundPlasmaBoundaryIonInjection {

      //the number density fraction of particular ion species in the background plasma flow
      static const double IonNumberDensityFraction[]={1.0}; //the number density fraction of the particular ion species in the total background plasma flow (IonNumberDensityFraction[spec])
      static const double vmax=1000.0E3; //the maximum velocity of the injected ions

      //a cartesian grid is introduces on each face. the number of the cells on that grid corresponds to that is cells in the block to improve modeling of the variatino of the injectino rate
      static const int nFaceInjectionIntervals=std::max(_BLOCK_CELLS_X_,std::max(_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_));
      extern long int nTotalBoundaryInjectionFaces;  //the number of the computationsl mesh faces at the boundary of the domain

      extern double **BoundaryFaceTotalInjectionRate; //the total production rate that is due to a particular block (BoundaryBlockProductionFraction[spec][nInjectionFace]
      extern double *maxBoundaryFaceTotalInjectionRate; //the maximum value of the boundary face injectino rate

      extern double ***BoundaryFaceLocalInjectionRate; //flux at defferent locations across the face BoundaryFaceLocalInjectionRate[spec][nInjectionFace][SurfaceElement]
      extern double **maxBoundaryFaceLocalInjectionRate; //max production rate through a boundary face (a cartisian mesh is introduced on block faces) [spec][nInjectionFace]

      extern double *maxLocalTimeStep; //the maximum value of the time step across the boundary of the computational domain (maxLocalTimeStep[spec])
      extern double *minParticleWeight; //the minimum value of the particle weight across the computational domain (minParticleWeight[spec])

      extern double *TotalInjectionRateTable; //the rate of the ion injection through the boundary of the computational domain; the table is used to save the data calculated by GetTotalProductionRate for a steady state case. In a time-dependent case the value of the source rate is re-calculated each time when that function is acceced

      class cBoundaryFaceDescriptor {
      public:
        int nface;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

        cBoundaryFaceDescriptor() {
          nface=-1,node=NULL;
        }
      };

      extern cBoundaryFaceDescriptor *BoundaryFaceDescriptor;

      void getMinMaxLimits();
      double GetTotalProductionRate(int spec);

      long int ParticleInjection(int spec);
      long int ParticleInjection();

    }

    namespace ImpactVaporization {
      static const double ImpactVaporization_SourceRate[]={0.0};
      static const double ImpactVaporization_HeliocentricDistance=1.0*_AU_;
      static const double ImpactVaporization_SourceRatePowerIndex=0.0;
      static const double ImpactVaporization_SourceTemeprature[]={0.0};

      double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer);

      //evaluate nemerically the source rate
//      extern double CalculatedTotalSodiumSourceRate;

      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
        unsigned int idim;
        double r=0.0,vbulk[3]={0.0,0.0,0.0},ExternalNormal[3];


        //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
        double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
        SpiceDouble xform[6][6];

        memcpy(xform,OrbitalMotion::IAU_to_SO_TransformationMartix,36*sizeof(double));

        //Geenrate new particle position
        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
          r+=pow(ExternalNormal[idim],2);
        }

        r=sqrt(r);

        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]/=r;
          x_LOCAL_IAU_OBJECT[idim]=sphereX0[idim]-sphereRadius*ExternalNormal[idim];
        }

        //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
        x_LOCAL_SO_OBJECT[0]=xform[0][0]*x_LOCAL_IAU_OBJECT[0]+xform[0][1]*x_LOCAL_IAU_OBJECT[1]+xform[0][2]*x_LOCAL_IAU_OBJECT[2];
        x_LOCAL_SO_OBJECT[1]=xform[1][0]*x_LOCAL_IAU_OBJECT[0]+xform[1][1]*x_LOCAL_IAU_OBJECT[1]+xform[1][2]*x_LOCAL_IAU_OBJECT[2];
        x_LOCAL_SO_OBJECT[2]=xform[2][0]*x_LOCAL_IAU_OBJECT[0]+xform[2][1]*x_LOCAL_IAU_OBJECT[1]+xform[2][2]*x_LOCAL_IAU_OBJECT[2];


        //determine if the particle belongs to this processor
        startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
        if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

        //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
        PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,ImpactVaporization_SourceTemeprature[spec],ExternalNormal,spec);


/*
//DEBUG -> injected velocity is normal to the surface

for (int i=0;i<3;i++)  v_LOCAL_IAU_OBJECT[i]=-ExternalNormal[i]*4.0E3;
//END DEBUG
*/


        //transform the velocity vector to the coordinate frame 'MSGR_SO'
        v_LOCAL_SO_OBJECT[0]=xform[3][0]*x_LOCAL_IAU_OBJECT[0]+xform[3][1]*x_LOCAL_IAU_OBJECT[1]+xform[3][2]*x_LOCAL_IAU_OBJECT[2]+
            xform[3][3]*v_LOCAL_IAU_OBJECT[0]+xform[3][4]*v_LOCAL_IAU_OBJECT[1]+xform[3][5]*v_LOCAL_IAU_OBJECT[2];

        v_LOCAL_SO_OBJECT[1]=xform[4][0]*x_LOCAL_IAU_OBJECT[0]+xform[4][1]*x_LOCAL_IAU_OBJECT[1]+xform[4][2]*x_LOCAL_IAU_OBJECT[2]+
            xform[4][3]*v_LOCAL_IAU_OBJECT[0]+xform[4][4]*v_LOCAL_IAU_OBJECT[1]+xform[4][5]*v_LOCAL_IAU_OBJECT[2];

        v_LOCAL_SO_OBJECT[2]=xform[5][0]*x_LOCAL_IAU_OBJECT[0]+xform[5][1]*x_LOCAL_IAU_OBJECT[1]+xform[5][2]*x_LOCAL_IAU_OBJECT[2]+
            xform[5][3]*v_LOCAL_IAU_OBJECT[0]+xform[5][4]*v_LOCAL_IAU_OBJECT[1]+xform[5][5]*v_LOCAL_IAU_OBJECT[2];

        memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
        memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
        memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
        memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));

        //set up the intermal energy if needed
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_

  #if _PIC_INTERNAL_DEGREES_OF_FREEDOM__TR_RELAXATION_MODE_  == _PIC_MODE_ON_
        PIC::IDF::InitRotTemp(ImpactVaporization_SourceTemeprature[spec],tempParticleData);
  #endif

  #if _PIC_INTERNAL_DEGREES_OF_FREEDOM__VT_RELAXATION_MODE_  == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
  #endif

#endif

        return true;
      }

    }

    namespace SolarWindSputtering {
      static const double SolarWindSputtering_Yield[]={0.0};
      static const double SolarWindSputtering_UserRequestedSourceRate[]={0.0};
      static const double SolarWindSputtering_minInjectionEnergy[]={0.0};
      static const double SolarWindSputtering_maxInjectionEnergy[]={0.0};

      extern double SourceRate[PIC::nTotalSpecies],maxLocalSourceRate[PIC::nTotalSpecies];

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution[PIC::nTotalSpecies];
      double GetSurfaceElementProductionRate(int nElement,int *spec);
      void Init();

      //evaluate the typical total source rate based on the typical parameters of the solar wind
      double TypicalIonFluxSputteringRate(int spec);

      //evaluate nemerically the source rate
//      extern double CalculatedTotalSodiumSourceRate;

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double norm_IAU_OBJECT[3],norm_SO_OBJECT[3];

        if (SphereDataPointer==NULL) return 0.0;

        //the sputtering source rate does not depend on the surface population and depends only on the yeild
//        if (((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[spec][SurfaceElement]<=0.0) return 0.0;

        memcpy(norm_IAU_OBJECT,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        //convert the normal vector to the 'SO' frame
        norm_SO_OBJECT[0]=
            (OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*norm_IAU_OBJECT[0])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*norm_IAU_OBJECT[1])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*norm_IAU_OBJECT[2]);

        norm_SO_OBJECT[1]=
            (OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*norm_IAU_OBJECT[0])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*norm_IAU_OBJECT[1])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*norm_IAU_OBJECT[2]);

        norm_SO_OBJECT[2]=
            (OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*norm_IAU_OBJECT[0])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*norm_IAU_OBJECT[1])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*norm_IAU_OBJECT[2]);

        //get the surface element that is pointer by the vectorm norm_SO_OBJECT
        long int nZenithElement,nAzimuthalElement,nd;

        ((cInternalSphericalData*)SphereDataPointer)->GetSurfaceElementProjectionIndex(norm_SO_OBJECT,nZenithElement,nAzimuthalElement);
        nd=((cInternalSphericalData*)SphereDataPointer)->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

        //return the local source rate
        double flux=((cInternalSphericalData*)SphereDataPointer)->SolarWindSurfaceFlux[nd];

        #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE_ == _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE__YIELD_
        return (flux>0.0) ? flux*SolarWindSputtering_Yield[spec]*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementArea[SurfaceElement] : 0.0;
        #elif _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE_ == _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE__USER_SOURCE_RATE_
        return (flux>0.0) ? SolarWindSputtering_UserRequestedSourceRate[spec]*flux*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementArea[SurfaceElement]/((cInternalSphericalData*)SphereDataPointer)->TotalSolarWindSurfaceFlux : 0.0;
        #elif _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE_ == _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE__UNIFORM_USER_SOURCE_RATE_
        return SolarWindSputtering_UserRequestedSourceRate[spec]*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementArea[SurfaceElement]/(4.0*Pi*pow(((cInternalSphericalData*)SphereDataPointer)->Radius,2));
        #else
        exit(__LINE__,__FILE__,"Error: the option is unknown");
        #endif
      }

      inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
        #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE_ == _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE__YIELD_
        return SourceRate[spec];
        #elif _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE_ == _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE__USER_SOURCE_RATE_
        return SolarWindSputtering_UserRequestedSourceRate[spec];
        #elif _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE_ == _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE__UNIFORM_USER_SOURCE_RATE_
        return SolarWindSputtering_UserRequestedSourceRate[spec];
        #else
        exit(__LINE__,__FILE__,"Error: the option is unknown");
        #endif
      }

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution[PIC::nTotalSpecies];
      double DefaultEnergyDistributionFunction(double e,int *spec);
      double EnergyDistributionFunction(double e,int *spec);

      //generate particle properties
      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
#endif

        if (BoundaryElementType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_) exit(__LINE__,__FILE__,"Error: particle ejection from a non-spehtical body is not implemented");

        cInternalSphericalData* Sphere=(cInternalSphericalData*)BoundaryElement;

        return Exosphere::SourceProcesses::GenerateParticleProperties(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,SurfaceInjectionDistribution+spec,EnergyDistribution+spec,_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_);
      }
    }


    double totalProductionRate(int spec,int BoundaryElementType,void *BoundaryElement);
    long int InjectionBoundaryModel(int BoundaryElementType,void *BoundaryElement);
    long int InjectionBoundaryModel(int spec,int BoundaryElementType,void *BoundaryElement);

    long int InjectionBoundaryModelLimited(void *SphereDataPointer);
    long int InjectionBoundaryModelLimited(int spec,void *SphereDataPointer);
  }



  //combine together the surface area densities and recalcualte the
  void ExchangeSurfaceAreaDensity();


  //Interaction of the particles with the surface
  namespace SurfaceInteraction {


    //sampling the surface area density of the sticking species
#define _OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_             3

#define _OBJECT_SURFACE_STICKING_SAMPLING_OFFSET__FLUX_DOWN_           0
#define _OBJECT_SURFACE_STICKING_SAMPLING_OFFSET__FLUX_UP_             1
#define _OBJECT_SURFACE_STICKING_SAMPLING_OFFSET__AREA_NUMBER_DENSITY_ 2


    //sodium/surface interaction model
    static const double AccomodationCoefficient[]={0.0};

    //sticking probability of sodium atoms
    double StickingProbability(int spec,double& ReemissionParticleFraction,double Temp);

    //model of the interaction between particles and the planetary surface
    int ParticleSphereInteraction_SurfaceAccomodation(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);
  }










}



#endif /* OBJECT_H_ */

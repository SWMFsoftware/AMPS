//$Id$

#ifndef _OH_H_
#define _OH_H_

// AMPS core code
#include "pic.h"

// Exosphere model
#include "Exosphere.h"

// self-explanatory
#include "constants.h"

// charge exchange physical model
#include "ChargeExchange.h"

// SPICE is not used for this application
// use empty functions to fill in
#include "SpiceEmptyDefinitions.h"

namespace OH {
  using namespace Exosphere;
  
  // user defined global time step
  extern double UserGlobalTimeStep;

  //  injection boundary condition
  extern double InjectionVelocity[3];
  extern double InjectionNDensity;
  extern double InjectionTemperature;
  
  // computational domain size
  extern double DomainXMin[3];
  extern double DomainXMax[3];
  extern double DomainDXMin;
  extern double DomainDXMax;

  //variables to contain and initialize the physical species number and index within the model
  extern int nPhysSpec;
  extern int PhysSpec[PIC::nTotalSpecies];
  extern bool DoPrintSpec[PIC::nTotalSpecies];
  void InitPhysicalSpecies();
  
  // creating a variable that holds location designator of where it was created
  extern long int OffsetOriginTag;

  // gets origin tag from each particle for output
  inline int GetOriginTag(PIC::ParticleBuffer::byte* ParticleData){
    return *((int*) (ParticleData + OffsetOriginTag));
  }

  // setting the origin tag for each particle
  // iRegion is the variable that holds the integer for each population
  inline void SetOriginTag(int iRegion, PIC::ParticleBuffer::byte* ParticleData){
    *((int*) (ParticleData + OffsetOriginTag)) = iRegion;
  }

  // determines what the origin tag for each particle should be based on loacl plasma parameters where it was created
  int GetEnaOrigin(double PlasmaNumberDensity, double PlasmaPressure, double *PlasmaBulkVelocity);
  void InitializeParticleWithEnaOriginTag(long int ptr, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode, int iInjectionMode);
  

  void Init_BeforeParser();
  
  int user_set_face_boundary(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
  
  //---------------------------------------------------------------------------
  namespace Sampling{
    using namespace Exosphere::Sampling;

    //sample the particle density separately for each source location
    namespace OriginLocation {
      extern int nSampledOriginLocations;

      //offset of the sampled density and velocity data in the AMPS' sampling vector
      extern int OffsetDensitySample;
      extern int OffsetVelocitySample;
      extern int OffsetV2Sample;

      //init the sampling module
      int RequestSamplingData(int offset);
      void SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec);

    }

    namespace DistributionFunctionSample{
      extern const int _LINEAR_SAMPLING_SCALE_,_LOGARITHMIC_SAMPLING_SCALE_;
      extern const bool Use;
      extern int v2SamplingMode,speedSamplingMode;
      extern double vMin;
      extern double vMax;
      extern long int nSampledFunctionPoints;
      extern double** SamplingBuffer;
      extern double SamplingLocations[][3];
      extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
      extern double dV,dV2,dSpeed;
      extern long int *SampleLocalCellNumber;
      extern int nSampleLocations;
      extern bool SamplingInitializedFlag;
      extern int Sample_Velocity_Offset,Sample_Speed_Offset;
      extern int Sample_V2_Offset,SampleDataLength;
      
      void Init();
      void flushSamplingBuffers();
      long int GetSampleDataOffset(int spec,int OriginID,int SampleVariableOffset);    
      void SampleDistributionFnction();
      void printDistributionFunction(int DataOutputFileNumber);
    }  

  }

  //---------------------------------------------------------------------------
  namespace Coupling {
    extern double TimeAfterCoupling[PIC::nTotalSpecies];
    void Send(char *NameVar, int *nVarIn, int *nDimIn, int *nPoint, double *Xyz_DI, double *Data_VI);
  }

  //---------------------------------------------------------------------------
  namespace Output{

    extern int TotalDataLength;
    extern int ohSourceDensityOffset; 
    extern int ohSourceMomentumOffset;
    extern int ohSourceEnergyOffset;
    
    void Init();

    void PrintVariableList(FILE* fout,int DataSetNumber);

    void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);

    int RequestDataBuffer(int offset);
  }

  //---------------------------------------------------------------------------
  namespace Loss {
    double LifeTime(double *x, int spec, long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    void ReactionProcessor(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

  }

  //---------------------------------------------------------------------------
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {

    accl[0]=0.0; accl[1]=0.0;  accl[2]=0.0; 

  }
  

}

#endif

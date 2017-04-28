//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===============================================
//$Id$
//===============================================
//the header for the new particle model

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <utility>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "global.h"

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#include <omp.h>
#endif //_PIC_COMPILATION_MODE_ == _PIC_COMPILATION_MODE__HYBRID_

using namespace std;

#ifndef _PIC_
#define _PIC_

//the global model settings
#include "picGlobal.dfn"

//if the code is compiler for the nightly tests than a) turn the debug model ON, and b) turn on the dynamic domain decomposition based on the particle number
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  #undef _PIC_DEBUGGER_MODE_
  #define _PIC_DEBUGGER_MODE_ _PIC_DEBUGGER_MODE_ON_

  #undef _PIC_DYNAMIC_LOAD_BALANCING_MODE_
  #define _PIC_DYNAMIC_LOAD_BALANCING_MODE_ _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
#endif


//load the macros that defined symbolic references to the species
#include "picSpeciesMacro.dfn"

//load the user defined settings
/*
#if _PIC_USER_DEFINITION_MODE_ == _PIC_USER_DEFINITION_MODE__ENABLED_
#include "UserDefinition.PIC.h"
#endif
*/



#if _PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE_ == _PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE__ON_
$MARKER:SPECIES-MACRO-DEFINIETION-USED-IN-SIMULATION$
#endif

#include "ifileopr.h"
#include "specfunc.h"
#include "constants.h"
#include "rnd.h"
#include "stack.h"

//include the appropriate mesh header
#if DIM == 3
#include "meshAMR3d.h"
#elif DIM == 2
#include "meshAMR2d.h"
#else 
#include "meshAMR1d.h"
#endif

#include "meshAMRinternalSurface.h"


namespace PIC {

  //Global constants of the PIC solver
  //extern int nTotalSpecies;
  static const int nTotalSpecies=1;

  //The currect and total number of processors used in the simulation
  extern int ThisThread,nTotalThreads;

  //the total number of the OpenMP threads (when OpneMP is used in the model run)
  extern int nTotalThreadsOpenMP;

  //The path to the input data of the user-defined physical models
  extern char UserModelInputDataPath[_MAX_STRING_LENGTH_PIC_];

  //the output and prefix for the diagnostic information
  extern FILE* DiagnospticMessageStream;
  extern char DiagnospticMessageStreamName[_MAX_STRING_LENGTH_PIC_];

  //the directory for the input and output files
  extern char OutputDataFileDirectory[_MAX_STRING_LENGTH_PIC_];
  extern char InputDataFileDirectory[_MAX_STRING_LENGTH_PIC_];


  //the variables that controls the sampling procedure
  //LastSampleLength: the length of the last collected sample
  //CollectingSampleCounter: the number of iterations passed in the current sample loop
  //RequireSampleLength: the maximum number of iterations required in the current sample loop
  //DataOutputFileNumber: the number of the NEXT output file
  //SamplingMode: if _RESTART_SAMPLING_MODE_ -> after each output of the data file, the sampling buffer is flushed, _ACCUMULATE_SAMPLING_MODE_ -> the sampling data are saved and used for the next output of the flow file

  #define _RESTART_SAMPLING_MODE_    0
  #define _ACCUMULATE_SAMPLING_MODE_ 1
  extern long int LastSampleLength,CollectingSampleCounter,RequiredSampleLength,DataOutputFileNumber;
  extern int SamplingMode;

  //the tags for the data exchenge between processors
  #define _PIC_SUBDOMAIN_BOUNDARY_LAYER_SAMPLING_DATA_EXCHANGE_TAG_   0
  #define _PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_     1

  //handle run time signals and exeptions
  void SignalHandler(int);

  //perform one time step
  int TimeStep();
//  void Sampling();

  //init the particle solver
  void InitMPI();
  void Init_BeforeParser();
  void Init_AfterParser();

  //the list of functions used to exchenge the execution statiscics
  typedef void (*fExchangeExecutionStatistics) (CMPI_channel*,long int);
  extern vector<fExchangeExecutionStatistics> ExchangeExecutionStatisticsFunctions;

  //structural elements of the solver
  namespace Parser {
    void Run(char*);
    void readMain(CiFileOperations&);
    void readGeneral(CiFileOperations&);

  }

  namespace Datum {
    class cDatum{
      // class with information about data being stored/sampled/printed:
      // - offset in the node's buffer
      // - length of the datum in units sizeof(double)
      // - name of the physical parameter
      // - type of datum (sampled with given averaging or derived)
      // - flag whether to print to the output file
      // also contains generic methods:
      // - activation of datum
      // - printing name to a file
    public:
      static const int Unset_    = 0; // basic datum
      static const int Stored_   = 1; // stored, NOT sampled
      static const int Sampled_  = 2; // sampled, averaging is not yet defined
      static const int Timed_    = 3; // sampled, averaged over time
      static const int Weighted_ = 4; // sampled, averaged over particle weight
      static const int Derived_  = 5; // derived from other types

      long int offset;
      int length;
      char name[_MAX_STRING_LENGTH_PIC_];
      int type;
      bool doPrint;
      
      // print variables' name to file
      //......................................................................
      inline void PrintName(FILE* fout){fprintf(fout, ", %s", name);}

      // constructor
      //......................................................................
      cDatum(int lengthIn, const char* nameIn, bool doPrintIn = true) {
	      // nameIn must be in the format acceptable for printing output:
	      //   "\"var_1\", \"var_2\""
	      length = lengthIn;
	      sprintf(name, "%s", nameIn);

	      // mark as inactive by default
	      offset = -1;
	      type = Unset_;
	      doPrint= doPrintIn;
      }
    };
    // class cDatum -----------------------------------------------------------

    //-------------------------------------------------------------------------
    // the following classes are treated differently at the calls/sampling
    // and activation:
    // - cDatumStored    is stored, but NOT sampled
    // - cDatumSampled   is a parent of the following 2:
    // -- cDatumTimed    is averaged over time
    // -- cDatumWeighted is averaged over particle weight on this node
    //-------------------------------------------------------------------------
    class cDatumStored : public cDatum {
    public:
      // activation procedures:
      // appart from the offset in the data buffer,
      // a storage for data info to keep track of data being sampled
      // at different buffers
      //......................................................................
      inline bool is_active() {return offset >= 0;}

      inline void activate(long int& offsetInOut, vector<cDatumStored*>* DatumVector) {
	      if(is_active()) exit(__LINE__,__FILE__,"ERROR: trying to activate datum a second time");

        // set offset to the variable
	      offset=offsetInOut;

       	// return info about length of the variable
	      offsetInOut+=length*sizeof(double);

	      // add this datum to the provided cDatum vector
	      DatumVector->push_back(this);
      }

      // constructor is inherited
      cDatumStored(int lengthIn, const char* nameIn, bool doPrintIn = true) : cDatum(lengthIn, nameIn, doPrintIn) {
        type = Stored_;
      }

    };

    class cDatumSampled : public cDatum {
    public:
      // activation procedures:
      // appart from the offset in the data buffer,
      // a storage for data info to keep track of data being sampled
      // at different buffers
      //......................................................................
      inline bool is_active() {return offset >= 0;}

      inline void activate(long int& offsetInOut, vector<cDatumSampled*>* DatumVector) {
	      if(is_active()) exit(__LINE__,__FILE__,"ERROR: trying to activate datum a second time");

        // set offset to the variable
	      offset=offsetInOut;

      	// return info about length of the variable
	      // RESERVE SPACE FOR ALL SPECIES
	      offsetInOut += length*sizeof(double)*PIC::nTotalSpecies;

	      // add this datum to the provided cDatum vector
	      DatumVector->push_back(this);
      }

      // constructor is inherited
      cDatumSampled(int lengthIn, const char* nameIn, bool doPrintIn = true) : cDatum(lengthIn, nameIn, doPrintIn) {
        type = Sampled_;
      }
    };

    class cDatumTimed : public cDatumSampled {
    public:
      // constructor is inherited as well
      cDatumTimed(int lengthIn, const char* nameIn, bool doPrintIn = true) : cDatumSampled(lengthIn, nameIn, doPrintIn) {
        type = Timed_;
      }
    };

    class cDatumWeighted : public cDatumSampled {
    public:
      // constructor is inherited as well
      cDatumWeighted(int lengthIn, const char* nameIn, bool doPrintIn = true) : cDatumSampled(lengthIn, nameIn, doPrintIn) {
        type = Weighted_;
      }
    };
    //-------------------------------------------------------------------------

  } 
  // namespace Datum ----------------------------------------------------------

  //field line
  namespace FieldLine{

    //constants for house-keeping--------------------------------------
    //self-explanatory: everything                Vertex  Segment  Line
    const int OK_        = 0;                   // yes     yes      yes
    //isn't physically set, e.g. no coordinates    
    const int Unset_     =-1;                   // yes     yes      yes
    //no neighbors
    const int Hanging_   =-2;                   // yes     yes      no  
    //zero length
    const int Collapsed_ =-3;                   // no      yes      no
    //unpredictable error
    const int Error_     =-4;                   // no      yes      yes
    //error in connectivity
    const int Broken_    =-5;                   // no      no       yes
    //-----------------------------------------------------------------

    // alias
    typedef PIC::Datum::cDatum         cDatum;
    typedef PIC::Datum::cDatumStored   cDatumStored;
    typedef PIC::Datum::cDatumSampled  cDatumSampled;
    typedef PIC::Datum::cDatumTimed    cDatumTimed;
    typedef PIC::Datum::cDatumWeighted cDatumWeighted;

    // standard set of data that is stored/sampled
    extern cDatumStored   DatumAtVertexElectricField;
    extern cDatumStored   DatumAtVertexMagneticField;
    extern cDatumStored   DatumAtVertexPlasmaVelocity;
    extern cDatumStored   DatumAtVertexPlasmaDensity;
    extern cDatumStored   DatumAtVertexPlasmaTemperature;
    extern cDatumStored   DatumAtVertexPlasmaPressure;
    extern cDatumStored   DatumAtVertexMagneticFluxFunction;
    extern cDatumTimed    DatumAtVertexParticleWeight;
    extern cDatumTimed    DatumAtVertexParticleNumber;
    extern cDatumTimed    DatumAtVertexNumberDensity;
    extern cDatumWeighted DatumAtVertexParticleEnergy;
    extern cDatumWeighted DatumAtGridParticleEnergy;



    // vectors with active data
    extern vector<cDatumStored*> DataStoredAtVertex;
    extern vector<cDatumSampled*> DataSampledAtVertex;

    class cFieldLineVertex;
    class cFieldLineSegment;
    class cFieldLine;
  
    extern cFieldLine* FieldLinesAll;
    extern cAssociatedDataAMRstack<cFieldLineVertex> VerticesAll;
    extern cAssociatedDataAMRstack<cFieldLineSegment> SegmentsAll;
    

  class cFieldLineVertex{
    private:
      //flag whether coords of vertex have been set
      char IsSet;
      //coordinates of the vertex
      double x[DIM];
      // data buffer length
      static int totalAssociatedDataLength;
      static int sampleDataLength;
      // sampling offsets
      static int CollectingSamplingOffset;
      static int CompletedSamplingOffset;
      // pinter to the data buffer itself
      char* AssociatedDataPointer;
      //neighboring vertices
      cFieldLineVertex* prev;
      cFieldLineVertex* next;

    public:
      long int Temp_ID;
      bool ActiveFlag; //used to prevent multiple deallocation of the Vertix 

      cFieldLineVertex() {
        IsSet=0;

        for (int idim=0;idim<DIM;idim++) x[idim]=0;

        prev=NULL,next=NULL;
        Temp_ID=0;
        ActiveFlag=false;
      }

      //.......................................................................
      // interface with stack functionality
      inline int AssociatedDataLength() {
	      return totalAssociatedDataLength;
      }

      inline char* GetAssociatedDataBufferPointer() {
	      return AssociatedDataPointer;
      }

      inline void SetAssociatedDataBufferPointer(char* ptr) {
	      AssociatedDataPointer=ptr;
      }

      inline void cleanDataBuffer() {
        int i,length=totalAssociatedDataLength/sizeof(double);
        double *p;
        for (i=0,p=(double*)AssociatedDataPointer;i<length;i++,p++) *p=0.0;
      }

      //.......................................................................
      // offset manipulation
      inline static void SetDataOffsets(int SamplingOffset, int SampleDataLengthIn) {
	      CollectingSamplingOffset = SamplingOffset;
	      CompletedSamplingOffset  = SamplingOffset +   SampleDataLengthIn;
	      totalAssociatedDataLength= SamplingOffset + 2*SampleDataLengthIn;
	      sampleDataLength         = SampleDataLengthIn;
      }

      inline void flushCompletedSamplingBuffer() {
	      int i,length=sampleDataLength/sizeof(double);
	      double *ptr;

	      for(i=0,ptr=(double*)(AssociatedDataPointer+CompletedSamplingOffset);i<length;i++,ptr++) *ptr=0.0;
      }

      inline void flushCollectingSamplingBuffer() {
	      int i,length=sampleDataLength/sizeof(double);
	      double *ptr;

	      for(i=0,ptr=(double*)(AssociatedDataPointer+CollectingSamplingOffset);i<length;i++,ptr++) *ptr=0.0;
      }

      inline static void swapSamplingBuffers() {
	      int tempOffset = CollectingSamplingOffset;
	      CollectingSamplingOffset = CompletedSamplingOffset;
	      CompletedSamplingOffset = tempOffset;
      }

      //.......................................................................
      //check status of vertex
      inline int status() {
	      if (IsSet == 0) return Unset_;
        if (prev == NULL && next == NULL) return Hanging_;
	      return OK_;
      }

      //.......................................................................
      //status of vertex as a string
      inline void status(char* res) {
	      if (IsSet == 0)                  {sprintf(res,"%s","Unset");  return;}
	      if (prev == NULL && next == NULL){sprintf(res,"%s","Hanging");return;}
	      sprintf(res,"%s","OK"); return;
      }

      //.......................................................................
      //access to coordinates
      inline void SetX(double* xIn) {
	      IsSet = 1;
	      for(int idim=0; idim<DIM; idim++) x[idim]=xIn[idim];
      }

      inline void GetX(double* xOut) {
	      for(int idim=0; idim<DIM; idim++) xOut[idim]=x[idim];
      }

      //.......................................................................
      //set individual stored variables
      inline void SetDatum(cDatumStored Datum, double* In) {
	      memcpy(AssociatedDataPointer+Datum.offset, In, Datum.length * sizeof(double));
      }

      inline void SetDatum(cDatumStored Datum, double In) {
	      memcpy(AssociatedDataPointer+Datum.offset, &In, sizeof(double));
      }

      inline void SetElectricField(double* ElectricFieldIn) {
	      SetDatum(DatumAtVertexElectricField, ElectricFieldIn);
      }

      inline void SetMagneticField(double* MagneticFieldIn) {
	      SetDatum(DatumAtVertexMagneticField, MagneticFieldIn);
      }

      inline void SetPlasmaVelocity(double* PlasmaVelocityIn) {
	      SetDatum(DatumAtVertexPlasmaVelocity, PlasmaVelocityIn);
      }

      inline void SetPlasmaDensity(double  PlasmaDensityIn) {
	      SetDatum(DatumAtVertexPlasmaDensity, PlasmaDensityIn);
      }

      inline void SetPlasmaTemperature(double  PlasmaTemperatureIn) {
	      SetDatum(DatumAtVertexPlasmaTemperature, PlasmaTemperatureIn);
      }

      inline void SetPlasmaPressure(double  PlasmaPressureIn) {
	      SetDatum(DatumAtVertexPlasmaPressure, PlasmaPressureIn);
      }

      //.......................................................................
      //set individual sampled variables
      inline void SetDatum(cDatumSampled Datum, double* In, int spec) {
	      memcpy(AssociatedDataPointer + CompletedSamplingOffset + Datum.offset + Datum.length * spec * sizeof(double),In, Datum.length * sizeof(double));
      }

      inline void SetDatum(cDatumSampled Datum, double In, int spec) {
	      memcpy(AssociatedDataPointer + CompletedSamplingOffset + Datum.offset + Datum.length * spec * sizeof(double),&In, sizeof(double));
      }

      //.......................................................................
      // sample data
      inline void SampleDatum(PIC::Datum::cDatumSampled Datum,double* In, int spec, double weight=1.0) {
        for(int i=0; i<Datum.length; i++)  *(i + Datum.length * spec + (double*)(AssociatedDataPointer + CollectingSamplingOffset+Datum.offset))+= In[i] * weight;
      }

      inline void SampleDatum(Datum::cDatumSampled Datum, double In, int spec,double weight=1.0) {
        *(spec + (double*)(AssociatedDataPointer + CollectingSamplingOffset+Datum.offset))+= In * weight;
      }

      //.......................................................................
      //get individual stored variables
      inline void GetDatum(cDatumStored Datum, double* Out) {
	       memcpy(Out, AssociatedDataPointer+Datum.offset, Datum.length * sizeof(double));
      }

      inline void GetDatum(cDatumStored Datum, double& Out) {
        Out = *(double*)(AssociatedDataPointer+Datum.offset);
      }

      inline void GetElectricField(double* ElectricFieldOut) {
	      GetDatum(DatumAtVertexElectricField, ElectricFieldOut);
      }

      inline void GetMagneticField(double* MagneticFieldOut) {
	      GetDatum(DatumAtVertexMagneticField, MagneticFieldOut);
      }

      inline void GetPlasmaVelocity(double* PlasmaVelocityOut) {
	      GetDatum(DatumAtVertexPlasmaVelocity, PlasmaVelocityOut);
      }

      inline void GetPlasmaDensity(double& PlasmaDensityOut) {
	      GetDatum(DatumAtVertexPlasmaDensity, &PlasmaDensityOut);
      }

      inline void GetPlasmaTemperature(double& PlasmaTemperatureOut) {
	      GetDatum(DatumAtVertexPlasmaTemperature, &PlasmaTemperatureOut);
      }

      inline void GetPlasmaPressure(double& PlasmaPressureOut) {
	      GetDatum(DatumAtVertexPlasmaPressure, &PlasmaPressureOut);
      }

      //get accumulated data
      //.......................................................................
      inline void GetDatumCumulative(Datum::cDatumSampled Datum, double* Out, int spec) {
	      for (int i=0; i<Datum.length; i++) Out[i] = *(i + Datum.length * spec + (double*)(AssociatedDataPointer + CompletedSamplingOffset+Datum.offset));
      }

      inline double GetDatumCumulative(Datum::cDatumSampled Datum, int spec) {
	      return *(spec + (double*)(AssociatedDataPointer + CompletedSamplingOffset+Datum.offset));
      }

      //get data averaged over time
      //.......................................................................
      inline void GetDatumAverage(cDatumTimed Datum, double* Out, int spec) {
	      if (PIC::LastSampleLength > 0) {
	        for (int i=0; i<Datum.length; i++) Out[i] = *(i + Datum.length * spec + (double*)(AssociatedDataPointer + CompletedSamplingOffset+Datum.offset)) / PIC::LastSampleLength;
	      }
	      else for(int i=0; i<Datum.length; i++) Out[i] = 0.0;
      }

      inline double GetDatumAverage(cDatumTimed Datum, int spec) {
	      return (PIC::LastSampleLength > 0) ? *(spec + (double*)(AssociatedDataPointer + CompletedSamplingOffset+Datum.offset)) / PIC::LastSampleLength : 0.0;
      }

      //get data averaged over sampled weight
      //.......................................................................
      inline void GetDatumAverage(cDatumWeighted Datum, double* Out, int spec) {
	      double TotalWeight=0.0;

	      GetDatumCumulative(DatumAtVertexParticleWeight, &TotalWeight, spec);

	      if(TotalWeight > 0) for(int i=0; i<Datum.length; i++) Out[i] = *(i + Datum.length * spec + (double*)(AssociatedDataPointer + CompletedSamplingOffset+ Datum.offset)) / TotalWeight;
	      else for(int i=0; i<Datum.length; i++) Out[i] = 0.0;
      }

      inline double GetDatumAverage(cDatumWeighted Datum, int spec) {
	      double TotalWeight=0.0;

	      GetDatumCumulative(DatumAtVertexParticleWeight, &TotalWeight, spec);

	      return (TotalWeight > 0) ? *(spec +(double*)(AssociatedDataPointer + CompletedSamplingOffset+Datum.offset)) / TotalWeight : 0.0;
      }

      //.......................................................................
      //access to neighbors
      inline void SetPrev(cFieldLineVertex* prevIn) {prev=prevIn;}
      inline void SetNext(cFieldLineVertex* nextIn) {next=nextIn;}
      inline void GetPrev(cFieldLineVertex* prevOut) {prevOut=prev;}
      inline void GetNext(cFieldLineVertex* nextOut) {nextOut=next;}
      inline cFieldLineVertex* GetPrev() {return prev;}
      inline cFieldLineVertex* GetNext() {return next;}
    };

    //class cFieldLineVertex --------------------------------------------------
    class cFieldLineSegment{
    private:
      //flag segment has been set (i.e. both vertices are set)
      char IsSet;

      //length of this segment
      double length;

      // weight of the segment for each species,
      // needed for injection of particles
      double weight[PIC::nTotalSpecies];

      //direction of this segment
      double Dir[DIM];

      //neighboring segments
      cFieldLineSegment* prev;
      cFieldLineSegment* next;

      //segment's vertices
      cFieldLineVertex* begin;
      cFieldLineVertex* end;
    public:
      long int Temp_ID;

      cFieldLineSegment(){
	      Temp_ID = 0;
	      IsSet=0, length=0.0;
	      prev  = (next = NULL);
	      begin = (end  = NULL);
      }

      //.......................................................................
      // interface with stack functionality
      inline void cleanDataBuffer(){}
      inline int AssociatedDataLength(){return 0;}
      inline char* GetAssociatedDataBufferPointer(){return NULL;}
      inline void SetAssociatedDataBufferPointer(char* ptr){}

      //.......................................................................
      //check status of segment
      inline int status() {
	      if (IsSet == 0)                return Unset_;
	      if (prev==NULL && next==NULL)  return Hanging_;
	      if (begin==end || length==0.0) return Collapsed_;
	      if (length < 0.0)              return Error_;
	      return OK_;
      }

      //.......................................................................
      //status of segment as string
      inline void status(char* res) {
	      if (IsSet == 0)               {sprintf(res,"%s","Unset");    return;}
	      if (prev==NULL && next==NULL) {sprintf(res,"%s","Hanging");  return;}
	      if (begin==end || length==0.0){sprintf(res,"%s","Collapsed");return;}
	      if (length < 0.0)             {sprintf(res,"%s","Error");    return;}
	      sprintf(res,"%s","OK"); return;
      }

      //.......................................................................
      //set the segment
      inline void SetVertices(cFieldLineVertex* beginIn,cFieldLineVertex* endIn) {

        #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
	      if (beginIn->status() != OK_ || endIn->status() != OK_) {
	        char msg[600];
	        char statusBegin[10],statusEnd[10];

	        beginIn->status(statusBegin), endIn->status(statusEnd);
	        sprintf(msg,"ERROR:: trying to set a segment with invalid vertices, beginIn->status()=%s, endIn->status()=%s",statusBegin, statusEnd);
	        exit(__LINE__,__FILE__, msg);
	      }
        #endif//_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_

	      IsSet=1;
	      if (beginIn!=NULL) begin = beginIn;
	      if (endIn  !=NULL) end   = endIn;

	      double xBegin[DIM], xEnd[DIM];

	      length = 0.0;

	      if (begin!=NULL && end!=NULL) {
	        begin->GetX(xBegin); end->GetX(xEnd);

	        for(int idim=0; idim<DIM; idim++) {
	          Dir[idim] = xEnd[idim]-xBegin[idim];
	          length   += pow(Dir[idim], 2);
	        }

	        length = sqrt(length);
	        Dir[0]/= length; Dir[1]/= length; Dir[2]/= length;
	      }
	      else IsSet = 0;
      }

      //.......................................................................
      //access segment's vertices
      inline void GetBegin(cFieldLineVertex* beginOut) {beginOut = begin;}
      inline void GetEnd(  cFieldLineVertex* endOut  ) {endOut   = end;}
      inline cFieldLineVertex* GetBegin() {return begin;}
      inline cFieldLineVertex* GetEnd() {return end;}

      //.......................................................................
      //access segment's length and and directoon at its beginning
      inline double GetLength() {return length;}
      inline void GetDir(double* DirOut) {memcpy(DirOut,Dir,DIM*sizeof(double));}

      //.......................................................................
      //access segment's statistical weight
      inline double GetWeight(int spec) {return weight[spec];}

      inline void   SetWeight(double* weightIn) {
	      memcpy(weight, weightIn, PIC::nTotalSpecies*sizeof(double));
      }

      inline void   SetWeight(double weightIn, int spec) {
	      weight[spec] = weightIn;
      }

      //.......................................................................
      //access segment's neighbors
      inline void SetPrev(cFieldLineSegment* prevIn){prev = prevIn;}
      inline void SetNext(cFieldLineSegment* nextIn){next = nextIn;}
      inline void GetPrev(cFieldLineSegment* prevOut){prevOut = prev;}
      inline void GetNext(cFieldLineSegment* nextOut){nextOut = next;}
      inline cFieldLineSegment* GetPrev(){return prev;}
      inline cFieldLineSegment* GetNext(){return next;}

      //.......................................................................
      //interpolate individual variables from vertices to point on the segment
      //position 0<=s<=1 on the segment
      inline void GetElectricField(double  S, double* ElectricFieldOut) {
	      double tBegin[DIM], tEnd[DIM];

	      begin->GetElectricField( tBegin), end->GetElectricField(tEnd);

	      for (int idim=0; idim<DIM; idim++) {
	        ElectricFieldOut[idim]  = tBegin[idim] * (1 - S) + tEnd[idim] * S;
	      }
      }

      //position 0<=s<=1 on the segment
      inline void GetMagneticField(double  S, double* MagneticFieldOut) {
	      double tBegin[DIM], tEnd[DIM];

	      begin->GetMagneticField( tBegin), end->GetMagneticField(tEnd);

	      for (int idim=0; idim<DIM; idim++) {
	        MagneticFieldOut[idim]  = tBegin[idim] * (1 - S) + tEnd[idim] * S;
	      }
      }

      //position 0<=s<=1 on segment
      inline void GetPlasmaVelocity(double  S, double* PlasmaVelocityOut) {
	      double tBegin[DIM], tEnd[DIM];

	      begin->GetPlasmaVelocity(tBegin), end->GetPlasmaVelocity(tEnd);

	      for(int idim=0; idim<DIM; idim++) {
	        PlasmaVelocityOut[idim] = tBegin[idim] * (1 - S) + tEnd[idim] * S;
	      }
      }

      //position 0<=s<=1 on the segment
      inline void GetPlasmaDensity(double  S, double& PlasmaDensityOut) {
	      double tBegin, tEnd;

	      //plasma density
	      begin->GetPlasmaDensity(tBegin), end->GetPlasmaDensity(tEnd);
	      PlasmaDensityOut = tBegin * (1 - S) + tEnd * S;
      }

      //position 0<=s<=1 on segment
      inline void GetPlasmaTemperature(double  S, double& PlasmaTemperatureOut) {
	      double tBegin, tEnd;

	      //plasma temperature
	      begin->GetPlasmaTemperature(tBegin), end->GetPlasmaTemperature(tEnd);
	      PlasmaTemperatureOut = tBegin * (1 - S) + tEnd * S;
      }

      //position 0<=s<=1 on segment
      inline void GetPlasmaPressure(double  S, double& PlasmaPressureOut) {
	      double tBegin, tEnd;

	      //plasma pressure
	      begin->GetPlasmaPressure(tBegin), end->GetPlasmaPressure(tEnd);
	      PlasmaPressureOut = tBegin * (1 - S) + tEnd * S;
      }
    };

    //class cFieldLineSegment -------------------------------------------------
    class cFieldLine {
    private:
      //flag whether line has been set (1 or more valid segments)
      char IsSet;
      //total number of segments (for house-keeping)
      int nSegment;
      //total length
      double TotalLength;
      //total statistical weight
      double TotalWeight[PIC::nTotalSpecies];
      //1st and last segments, vertices of the field line
      cFieldLineSegment *FirstSegment, *LastSegment;
      cFieldLineVertex  *FirstVertex,  *LastVertex;
      //check whether the line is broken
      bool is_broken();
    public:

      cFieldLine() {
	      IsSet = 0, nSegment = -1, TotalLength=-1.0;
	      FirstSegment = (LastSegment = NULL);
	      FirstVertex  = (LastVertex  = NULL);
      }

      //check status of line
      inline int status() {
	      if (IsSet == 0)                       return Unset_;
	      if (is_broken())                      return Broken_;
	      if (TotalLength < 0.0 || nSegment < 0)return Error_;
	      return OK_;
      }

      //status of segment as string
      inline void status(char* res) {
	      if(IsSet == 0)                   {sprintf(res,"%s","Unset"); return;}
	      if(is_broken())                  {sprintf(res,"%s","Broken");return;}
	      if(TotalLength<0.0 || nSegment<0){sprintf(res,"%s","Error"); return;}
	      sprintf(res,"%s","OK"); return;
      }

      //check whether the line is a loop
      inline void close_loop() {
	      VerticesAll.deleteElement(LastVertex);

	      LastVertex = FirstVertex;
	      LastVertex->SetPrev(LastSegment->GetBegin());
	      LastSegment->GetBegin()->SetNext(LastVertex);
	      TotalLength -= LastSegment->GetLength();
	      LastSegment->SetVertices(LastSegment->GetBegin(), LastVertex);
	      TotalLength += LastSegment->GetLength();
	      LastSegment->SetNext(FirstSegment);
	      FirstSegment->SetPrev(LastSegment);
      }

      inline bool is_loop() {return LastVertex==FirstVertex;}

      inline void fix_coord(double& Coord) {
	      while(Coord < 0)       Coord += nSegment;
	      while(Coord >=nSegment)Coord -= nSegment;
      }

      inline double move(double SInit, double Increment) {
	      double res = SInit;
	      cFieldLineSegment *Segment = GetSegment(SInit);
	      double Length = Segment->GetLength();
	      double remain;

	      if (Increment>0) {
	        remain = (int)(SInit+1) - SInit;

	        for (Segment = GetSegment(SInit); true; Segment = Segment->GetNext()) {
	          Length = Segment->GetLength();

	          if (Increment < remain*Length) {
	            res += Increment / Length;
	            break;
	          }

	          Increment -= remain*Length;
	          res += remain;
	          remain = 1.0;
	        }
	      }
	      else {
	        remain = SInit - (int)(SInit) ;

	        for (Segment = GetSegment(SInit);true; Segment = Segment->GetPrev()) {
	          Length = Segment->GetLength();

	          if (-Increment < remain*Length) {
	            res += Increment / Length;
	            break;
	          }

	          Increment += remain*Length;
	          res -= remain;
	          remain = 1.0;
	        }
	      }

	      return res;
      }

      // Segment access
      //-----------------------------------------------------------------------
      //access first/last segment
      inline cFieldLineSegment* GetFirstSegment() {return FirstSegment;}
      inline cFieldLineSegment* GetLastSegment() { return LastSegment;}
      inline void GetFirstSegment(cFieldLineSegment* Out) {Out=FirstSegment;}
      inline void GetLastSegment( cFieldLineSegment* Out) {Out=LastSegment;}

      //access an arbitrary segment
      inline cFieldLineSegment* GetSegment(int iSegment) {
	      cFieldLineSegment* Segment=NULL;

	      if (iSegment > 0.5*nSegment && iSegment < nSegment) {
          Segment = LastSegment;

          for(int i=nSegment-1; i > iSegment; i--) Segment = Segment->GetPrev();
        }

        if (iSegment >= 0) {
          Segment = FirstSegment;
          for(int i=0; i < iSegment; i++) Segment = Segment->GetNext();
        }

        if (Segment==NULL) exit(__LINE__,__FILE__, "ERROR: invalid index of a segment");

        return Segment;
      }

      inline cFieldLineSegment* GetSegment(double S) {
        // check correctness
        if(S < 0.0 || S > nSegment) return NULL;

        //  floor(S) is the index of the segment
        int iSegment = (int) S;
        return GetSegment(iSegment);
      };

      inline double GetSegmentLength(double S){
        return GetSegment(S)->GetLength();
      }

      inline void GetSegmentDirection(double* Dir, double S){
        GetSegment(S)->GetDir(Dir);
      }
      //-----------------------------------------------------------------------

      // Vertex access
      //-----------------------------------------------------------------------
      // access first/last vertex
      inline cFieldLineVertex* GetFirstVertex() {return FirstVertex;}
      inline cFieldLineVertex* GetLastVertex() {return LastVertex;}
      inline void GetFirstVertex(cFieldLineVertex* Out) {Out=FirstVertex;}
      inline void GetLastVertex( cFieldLineVertex* Out) {Out=LastVertex;}

      // access an arbitrary vertex
      inline cFieldLineVertex* GetVertex(int iVertex) {
        cFieldLineVertex* Vertex;

        if(iVertex > 0.5*nSegment && iVertex <= nSegment) {
          Vertex = LastVertex;

          for(int i=nSegment; i > iVertex; i--) Vertex = Vertex->GetPrev();
          return Vertex;
        }

        if(iVertex >= 0) {
          Vertex = FirstVertex;

          for(int i=0; i < iVertex; i++) Vertex = Vertex->GetNext();
          return Vertex;
        }

        exit(__LINE__,__FILE__, "ERROR: invalid index of a vertex");
      }

      //-----------------------------------------------------------------------
      //get cartesian coordinats of the location
      inline void GetCartesian(double* xOut, double S) {
        cFieldLineSegment* Segment = GetSegment(S);
        double w = S - (int)S;
        double xBegin[3], xEnd[3];

        Segment->GetBegin()->GetX(xBegin);
        Segment->GetEnd()  ->GetX(xEnd);
        for(int i=0; i<3; i++) xOut[i] = (1-w) * xBegin[i] + w * xEnd[i];
      }

      // add vertex with given coordinates
      cFieldLineVertex* Add(double* xIn);

      // assign statistical weights to segments and normalize them
      void ResetSegmentWeights();

      // get random segment
      //cFieldLineSegment* SegmentOut,
      void GetSegmentRandom(int& iSegment,double& WeightCorrectionFactor, int spec);

      // set magnetic field at a given vertex (last by default)
      void SetMagneticField(double* BIn, int iVertex=-1);

      // get background data at a given 1D coordinate along the field line
      //-----------------------------------------------------------------------
      // get magnetic field
      void   GetMagneticField(double* BOut, double S);
      // get directional derivative of magnetic field
      //-----------------------------------------------------------------------

      // print data stored on the field line
      void Output(FILE* fout, bool GeometryOnly);
    };
    //class cFieldLine --------------------------------------------------------


    // max number of field line in the simulation
    const long int nFieldLineMax=100;

    //current number of field lines
    extern long int nFieldLine;

    //time of last update
    extern double TimeLastUpdate;

    // initialize field line data structure
    void Init();

    //create parker spiral starting at point xStart
    void InitSimpleParkerSpiral(double *xStart);

    //create a 2D field-line loop based on the background field
    void InitLoop2D(double *xStart, double DArc, double DMin, double DMax);
    void Update();
    

    // output data
    void Output(char* fname, bool GeometryOnly);

    //functions for computing field-line segment weight
    void FieldLineWeight_Uniform(double* Weight, cFieldLineSegment* Segment);

    // sample data from a particle
    void Sampling(long int ptr, double Weight, char* SamplingBuffer);

    // inject particle onto the field line:
    // 1st function is a wrapper in the case
    // there is a user-defined procedure
    long int InjectParticle(int spec);
    long int InjectParticle_default(int spec);
  }

  //the first part of the namespace Debugger difinition
  namespace Debugger {
    //save a sequence of the particle data checksums into a file
    void SaveParticleDataIntoDebuggerDataStream(void*,int,int,const char*);
    void SaveParticleDataIntoDebuggerDataStream(void*,int,const char*);
  }

  //the particle trajectory tracing
  namespace ParticleTracker {
    extern long int ParticleDataRecordOffset; //the offset of the trajecotry specific informaion within the particle data
    extern int maxSampledTrajectoryNumber; //tha maximum number of the particle trajectories that will be sampled for each species
    extern int **threadSampledTrajectoryNumber; //the number of trajectories sampled by the current processor for each species. Format: [spec][threadOpenMP]
    extern int *totalSampledTrajectoryNumber; //the number of trajectories sampled by ALL PROCESSORS for each species. Format: [spec]
    extern unsigned long int *SampledTrajectoryCounter; //the total number of traced trajectories originate on the current processor -> used as a part of the trajecotry ID
    extern unsigned long int *SampledTrajectoryPointCounter; //the total number of sampled trajectory points

    extern int nMaxSavedSignleTrajectoryPoints; //the maximum number of the trajectory points saved in the output trajectory file for each particle trajectory

    extern bool AllowRecordingParticleTrajectoryPoints[PIC::nTotalSpecies]; //the array of flags that controls recording of the particle trajectory points

    struct cTrajectoryID {
      unsigned int StartingThread; //the thread where the trajectory has been originated
      unsigned int id; //the counting number of the trajecory on the processor where it was originated
    };


    struct cParticleData {
      bool TrajectoryTrackingFlag; //the trajectory of the praticle is sampled only when TrajectoryTrackingFlag==true; the default value is  TrajectoryTrackingFlag==false

      cTrajectoryID Trajectory;
      unsigned int nSampledTrajectoryPoints; //the number of the points of the particle trajectory
    };

    struct cTrajectoryPhysicalData {
      double x[3],Speed;
      int spec;

      #if _PIC_PARTICLE_TRACKER__TRAJECTORY_TIME_STAMP_MODE_ == _PIC_MODE_ON_
      double TimeStamp;
      #endif

      #if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
      double ElectricCharge,ParticleSize;
      #endif

      #if _PIC_MOVER_INTEGRATOR_MODE_ == _PIC_MOVER_INTEGRATOR_MODE__GUIDING_CENTER_
      double KineticEnergy;
      #endif

      #if _PIC_PARTICLE_TRACKER__INJECTION_FACE_MODE_ ==  _PIC_MODE_ON_
      int InjectionFaceNumber;
      #endif

      #if _PIC_PARTICLE_TRACKER__PARTICLE_WEIGHT_OVER_LOCAL_TIME_STEP_MODE_ == _PIC_MODE_ON_
      double ParticleWeightOverLocalTimeStepRatio;
      #endif
    };

    struct cTrajectoryDataRecord {
      cTrajectoryPhysicalData data;
      cTrajectoryID Trajectory;
      unsigned int offset; //the point number in the trajectory
    };

    struct cTrajectoryListRecord {
      cTrajectoryID Trajectory;
      unsigned int nSampledTrajectoryPoints; //the number of the points of the particle trajectory
    };

    struct cTrajectoryData {
      PIC::ParticleTracker::cTrajectoryDataRecord *buffer;
      unsigned long int CurrentPosition,nfile;
      static unsigned long int Size;

      void flush(); //save in a file the trajecotry information
      inline void clean() {CurrentPosition=0;nfile=0;} //reset to the initial set all counters

      cTrajectoryData() {
        buffer=NULL,CurrentPosition=0,nfile=0;
      }
    };

    extern cTrajectoryData *TrajectoryDataTable;

    struct cTrajectoryList {
      unsigned long int CurrentPosition,nfile;
      static unsigned long int Size;
      PIC::ParticleTracker::cTrajectoryListRecord *buffer;

      void flush(); //save in a file the list of the trajectories
      inline void clean() {CurrentPosition=0;nfile=0;} //reset to the initial set all counters

      cTrajectoryList() {
        buffer=NULL,CurrentPosition=0,nfile=0;
      }
    };

    extern cTrajectoryList *TrajectoryListTable;

    void Init();
    void UpdateTrajectoryCounter(); //calculate the total number of the trajectories sampled on all processors

    void InitParticleID(void *ParticleData);
    void RecordTrajectoryPoint(double *x,double *v,int spec,void *ParticleData,void *node);
    void FinilazeParticleRecord(void *ParticleData);

    //create the output trigectory file
    void OutputTrajectory(const char *fname);
    void CreateTrajectoryOutputFiles(const char *fname,const char *OutputDataFileDirectory,int TrajectoryPointBufferLength);

    void StartParticleTrajectoryTracking(void *ParticleData);
    void StopParticleTrajectoryTracking(void *ParticleData);

    void ApplyTrajectoryTrackingCondition(void *StartNodeVoid=NULL); //search all particles to start tracking those that met the search condition
    void ApplyTrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData,void *nodeIn); //apply the particle tracking conditions to the particle 'ParticleData'
    bool TrajectoryTrackingCondition_default(double *x,double *v,int spec,void *ParticleData); //the default tracking criterion -> return false for all particles
    void SetDefaultParticleTrackingFlag(void *StartNodeVoid=NULL); //set the value of 'TrajectoryTrackingFlag' to 'false' for all particles in a current simulation; used to re-set the trajectory sampling after output of the trajectory data file
  }


  //the interface used to use AMPS as a particle tracker in CCMC
  namespace CCMC {

  extern char BackgroundDataFileName[_MAX_STRING_LENGTH_PIC_];

  //maximum trajectory integraion time
  extern double MaxTrajectoryIntegrationTime;

  //internal boundary surface
  namespace InternalBoundary {
    extern int ParticleProcessingMode;

    //codes of the particle crossing internal boundary processing models
    namespace ParticleBoundaryInteractionCode {
      const int NoBoundary=-1;
      const int Sphere=0;
    }

    namespace Sphere {
      extern double Radius;
      extern double x0[3];
      extern int SamplingMode;
    }
  }


  //size of the computational domain
  namespace Domain {
    extern double xmin[3],xmax[3];

    //Local Cell Resolution
    namespace Resolution {

      //types of the resolution mode
      namespace TYPE {
        const int Constant=0;
        const int BackgroundFieldVariation=1;
        const int Logarithmic=2;
      }

      extern int mode;
      extern double BackgroundParameterVariationLimit;

      //the limits of the computational domain relative to its size
      extern double dxMin,dxMax;

      //location of dxMin,dxMax in case Logarithmic resolution is used
      extern double rXmin,rXmax;

      //get local resolution
      double GetlocalResolution(double *x);
    }
  }

  //characteristic speed of traced particles
  extern double *ParticleCharacteristicSpeedTable;
  double GetParticleCharacteristicSpeed(int spec);

  //constants
  namespace DEF {
    namespace SOURCE {
      namespace TYPE {
        const int Sphere=0;
        const int Table=1;
        const int Quadrilateral=2;
        const int Circle=3;
      }

      namespace SHPERE {
        namespace TYPE {
          const int Uniform=0;
          const int Gaussian=1;
        }
      }
    }

    namespace VELOCITY_DISTRIBUTION {
      namespace TYPE {
        const int Maxwellian=0;
        const int Table=1;
        const int Constant=2;
      }
    }
  }


    namespace ParticleInjection {

      //define the generic parameters that controls injection of the tracked particles
      class cInjectionControl {
      public:
        double StartTime;
        int nTestParticles;

        cInjectionControl() {
          StartTime=0.0;
          nTestParticles=0;
        }
      };

      //define the types of the particle generation for the tracking
      class cVelocityDistributionMaxwellian {
      public:
        double Temeprature,BulkVelocity[3];
      };

      class cVelocityDistributionTable {
      public:
        double v[3];
      };

      class cVelocityDistributionConstant {
      public:
        double v[3];
      };

      class cVelocityDistribution {
      public:
        int Type; //Table, maxwellian

        cVelocityDistributionMaxwellian Maxwellian;
        vector<cVelocityDistributionTable> Table;
        cVelocityDistributionConstant Constant;
      };

      class cInjectionRegionSpherical {
      public:
        double Radius,Origin[3];
        int SpatialDistributionType;  //uniform, gaussian
      };

      class cInjectionRegionCircle : public cInjectionRegionSpherical {
      public:
        double Normal[3],e0[3],e1[3];
      };

      class cInjectionRegionTable {
      public:
        double x[3];
      };

      class cInjectionRegionQuadrilateral {
      public:
        double xCenter[3],dX0[3],dX1[3]; //parameters of the input file: ceneter of the Quadrilateral, displacement of Quadrilateral's point 0 and 1 from the center
      };

      class cInjectionRegion {
      public:
        int Type; //spherical, Table

        cInjectionRegionSpherical Spherical;
        vector<cInjectionRegionTable> Table;
        cInjectionRegionQuadrilateral Quadrilateral;
        cInjectionRegionCircle Circle;
      };

      class cInjectionDescriptor : public cInjectionControl {
      public:
        cInjectionRegion SpatialDistribution;
        cVelocityDistribution VelocityDistribution;
      };

      extern vector<cInjectionDescriptor> InjectionDescriptorList;
    }

    namespace Parser {
      extern char ControlFileName[_MAX_STRING_LENGTH_PIC_];

      //read the file that describes the injection of the particle for the tracking
      void LoadControlFile();

      //read different sectoins of the input file
      namespace Read {
         namespace SourceRegion {
           void Sphere(PIC::CCMC::ParticleInjection::cInjectionDescriptor&,CiFileOperations&);
           void Circle(PIC::CCMC::ParticleInjection::cInjectionDescriptor&,CiFileOperations&);
           void Table(PIC::CCMC::ParticleInjection::cInjectionDescriptor&,CiFileOperations&);
           void Quadrilateral(PIC::CCMC::ParticleInjection::cInjectionDescriptor&,CiFileOperations&);
         }

         namespace VelocityDistribution {
           void Maxwellian(PIC::CCMC::ParticleInjection::cInjectionDescriptor&,CiFileOperations&);
           void Table(PIC::CCMC::ParticleInjection::cInjectionDescriptor&,CiFileOperations&);
           void Constant(PIC::CCMC::ParticleInjection::cInjectionDescriptor&,CiFileOperations&);
         }

         //read limits of the computational domain
         void DomainLimits(CiFileOperations&);

         //read chracteristic speed table
         void CharacteristicSpeedTable(CiFileOperations&);

         //read parameters of the internal boundary sphere
         void InternalBoundarySphere(CiFileOperations&);
      }
    }

    //trace particles
    void LoadParticles();
    int TraceParticles();
  }

  //ray tracing and calculation of the shadow regions on the NASTRAN surfaces
  namespace RayTracing {

    #if _PIC__RAY_TRACING__FACE_ACCESS_COUNTER_BYTE_LENGTH_ == 1
    extern unsigned char *FaceRayTracingOperationIDTable; //the number that is used to avoid checking of the same face twice
    extern unsigned char *FaceAccessCounterTable; //the conter of the Ray Tracking operations

    #elif _PIC__RAY_TRACING__FACE_ACCESS_COUNTER_BYTE_LENGTH_ == 4
    extern unsigned int *FaceRayTracingOperationIDTable; //the number that is used to avoid checking of the same face twice
    extern unsigned int *FaceAccessCounterTable; //the conter of the Ray Tracking operations

    #else
    #error _PIC__RAY_TRACING__FACE_ACCESS_COUNTER_BYTE_LENGTH_ is out of range
    #endif //_PIC__RAY_TRACING__FACE_ACCESS_COUNTER_BYTE_LENGTH_

    void Init();

    bool GetBlockExitPoint(double *xBlockMin,double *xBlockMax,double *x0Ray,double *lRay,double *xBlockExit, double *xFaceExitLocal, int &nExitFace);
    bool TestDirectAccess(double *xStart,double *xTarget);
    int CountFaceIntersectionNumber(double *xStart,double *xTarget,int MeshFileID,void* ExeptionFace=NULL);
    int FindFistIntersectedFace(double *x0Ray,double *lRay,double *xIntersection,void* ExeptionFace=NULL);

    void SetCutCellShadowAttribute(double *xLightSource, bool ParallelExecution=false);
    void FlushFaceRayTracingOperationIDTable(int iThreadOpenMP);
  }

  //define the test-run parameters
  namespace ModelTestRun {
    extern bool mode;
    extern int nTotalIteraction;
  }

  //run tame calculation of the check sums
  namespace RunTimeSystemState {
    void GetParticleFieldCheckSum(const char *msg=NULL);
    void GetParticleFieldCheckSum(long int nline,const char *fname);
    void GetParticleFieldCheckSum_CallCounter(const char *msg=NULL);
    void GetParticleFieldCheckSum_CallCounter(long int nline,const char *fname);

    void GetIndividualParticleFieldCheckSum_CallCounter(long int ptr,const char *msg=NULL);
    void GetIndividualParticleFieldCheckSum_CallCounter(void *ParticleData,const char *msg=NULL);
    void GetIndividualParticleFieldCheckSum_CallCounter(long int ptr,long int nline,const char *fname,const char *msg=NULL);
    void GetIndividualParticleFieldCheckSum_CallCounter(void *ParticleData,long int nline,const char *fname,const char *msg=NULL);

    void GetDomainDecompositionCheckSum(const char *msg=NULL);
    void GetDomainDecompositionCheckSum(long int nline,const char *fname);

    void GetMeanParticleMicroscopicParameters(FILE* fout,const char *msg=NULL);
    void GetMeanParticleMicroscopicParameters(FILE* fout,long int nline,const char *fname);
    void GetMeanParticleMicroscopicParameters(const char *fname);
  }


  namespace MolecularData {

    //Available molecular models
    #define _HS_MOL_MODEL_   0
    #define _VHS_MOL_MODEL_  1
    #define _VSS_MOL_MODEL_  2
    extern int MolModelCode;
    void SetMolType(int);
    int GetMolType();

    //modeling of external species
    #define _EXTERNAL_SPECIES_ON_  true
    #define _EXTERNAL_SPECIES_OFF_ false
    extern bool ExternalSpeciesModelingFlag;

    //modeling of unimolecular reactions
    #define _UNIMOLECULAR_REACTIONS_ON_   true
    #define _UNIMOLECULAR_REACTIONS_OFF_  false
    extern bool UnimolecularReactionFlag;

    //modeling of internal degrees of freedom
    #define _INTERNAL_DEGREES_OF_FREEDOM_ON_  true
    #define _INTERNAL_DEGRESS_OF_FREEDOM_OFF_ false
    extern bool InternalDegreesOfFreedomModelingFlag;

    //molecular models
    namespace MolecularModels {
      namespace HS {
        //the table of the constant collsion cross sections and reference diameter
        static const double ConstantCollisionCrossSectionTable[1][1]={{0.0}};
        static const double ConstantReferenceDiameter[1][1]={{0.0}};


        inline double GetTotalCrossSection(int s0,int s1) {return ConstantCollisionCrossSectionTable[s0][s1];}
        inline double GetDiam(int s0,int s1) {return ConstantReferenceDiameter[s0][s1];}
        inline double GetRefDiam(int s0,int s1) {return ConstantReferenceDiameter[s0][s1];}

      }

      namespace VHS {
      using namespace HS;

      }

      namespace VSS {
      using namespace VHS;

      }
    }

    inline double GetRefDiam(int s0,int s1) {
#if _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION_ == _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION__HS_
      return MolecularModels::HS::GetRefDiam(s0,s1);
#else
      exit(__LINE__,__FILE__,"not implemented");
      return 1.0;
#endif
    }

    inline double GetDiam(int s0,int s1) {
#if _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION_ == _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION__HS_
      return MolecularModels::HS::GetRefDiam(s0,s1);
#else
      exit(__LINE__,__FILE__,"not implemented");
      return 1.0;
#endif
    }

    inline double GetTotalCrossSect(double Vrel,int ptr0,int s0,int ptr1,int s1,long int ncell) {
#if _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION_ == _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION__HS_
      return MolecularModels::HS::GetTotalCrossSection(s0,s1);
#else
      exit(__LINE__,__FILE__,"not implemented");
      return 1.0;
#endif
    }


    //init the molecular data buffers
    void Init();

    //mass of particles
    static const double MolMass[]={0.0};
    static const double ElectricChargeTable[]={0.0};

//    extern double *MolMass;
//    void SetMass(double,int);
    inline double GetMass(int spec) {return MolMass[spec];}
    inline double GetElectricCharge(int spec) {return ElectricChargeTable[spec];}


    //get and set the value of the electric charge for a species
/*    extern double *ElectricCharge;
    int GetElectricCharge(int);
    void SetElectricCharge(double,int);*/

    //get and set the species numbers and chemical symbols
    static const char ChemTable[][_MAX_STRING_LENGTH_PIC_]={"nothin is defined"};


//    extern char **ChemTable; // <- The table of chemical symbols used in the simulation
    extern char **LoadingSpeciesList; // <- the list of species that CAN BE locased in the simualtion


    int GetSpecieNumber(char*);

    //get the chemical name of the species
    void GetChemSymbol(char*,int);
    const char* GetChemSymbol(int);

    //dust could have multiple groups defined as difference species in format DUST:0, DUST:1.....
    //GetChemBaseSymbol returns the species symbol till ':' - DUST
    //GetChemSymbol returns the full species symbol - DUST:0
    void GetChemBaseSymbol(char*,int);
    char* GetChemBaseSymbol(int);

    //set and get the specie type (gas, external, background)
    static const int SpcecieTypeTable[]={-1};

    inline int GetSpecieType(int spec) {
      return SpcecieTypeTable[spec];
    }


    namespace Parser {
      void InitChemTable(CiFileOperations&);
      int GetSpeciesNumber(int&,CiFileOperations&);
      void SpeciesBlock(char*,int,CiFileOperations&);
      void run(CiFileOperations&);
    }

  }

  //===========================================================================
  namespace ParticleBuffer {
    typedef unsigned char byte;
    
    // macro definition for particle data offsets
    #include "picParticleDataMacro.h"

    //the total length of a data allocated for a particle
    extern long int ParticleDataLength;

    //The particle buffer's internal data
    extern byte *ParticleDataBuffer;
    extern long int MaxNPart,NAllPart,FirstPBufferParticle;

    //Request additional data for a particle
    void RequestDataStorage(long int &offset,int TotalDataLength);

    //the basic data access functions for a particle
    byte *GetParticleDataPointer(long int);

    //check the total particles number
    void CheckParticleList();

    //the namespace contains data used in case when OpenMP is used
    namespace Thread {
      extern int NTotalThreads;
      extern long int *AvailableParticleListLength,*FirstPBufferParticle;
      void RebalanceParticleList();
    }

    //add the new particle to the simulation
    //argument 'node' should be type of cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*. but this make the compiler unhappy.
    //so, in the defineition of the function 'node' is void* but in pic_buffer.cpp this argument is transformed to cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*
    #define _PIC_INIT_PARTICLE_MODE__ADD2LIST_  0
    #define _PIC_INIT_PARTICLE_MODE__MOVE_      1
    int InitiateParticle(double *x,double *v,double* WeightCorrectionFactor,int *spec,byte* ParticleData,int InitMode,void *node);


    // Operations related to species ID
    //-------------------------------------------------------------------------
    // the first 7 bits will be used for specie ID, 
    // the last 8th bit will be used to control whether the particle 
    // has been allocated
    inline unsigned int GetI(byte* ParticleDataStart) {
      return ((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x7f);
    }
    //.........................................................................
    inline unsigned int GetI(long int ptr) {
      return ((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x7f);
    }
    //.........................................................................
    inline void SetI(int spec,byte* ParticleDataStart) {
      unsigned char flag,t=spec;

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&spec,sizeof(int),__LINE__,__FILE__);
      #endif

      flag=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x80);
      t|=flag;
      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))=t;
    }
    //.........................................................................
    inline void SetI(int spec,long int ptr) {
      unsigned char flag,t=spec;

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&spec,sizeof(int),__LINE__,__FILE__);
      #endif

      flag=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x80);
      t|=flag;
      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))=t;
    }
    //-------------------------------------------------------------------------

    // Operations related to the next particle in the stack
    //-------------------------------------------------------------------------
    inline long int GetNext(long int ptr) {
      return *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__NEXT_OFFSET_));
    }
    //.........................................................................
    inline long int GetNext(byte* ParticleDataStart) {
      return *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA__NEXT_OFFSET_));
    }
    //.........................................................................
    inline void SetNext(long int next,long int ptr) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&next,sizeof(long int),__LINE__,__FILE__);
      #endif

      *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__NEXT_OFFSET_))=next;
    }
    //.........................................................................
    inline void SetNext(long int next,byte* ParticleDataStart) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&next,sizeof(long int),__LINE__,__FILE__);
      #endif

      *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA__NEXT_OFFSET_))=next;
    }
    //-------------------------------------------------------------------------

    // Operations related to the previous particle in the stack
    //-------------------------------------------------------------------------
    inline long int GetPrev(long int ptr) {
      return *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__PREV_OFFSET_));
    }
    //.........................................................................
    inline long int GetPrev(byte* ParticleDataStart) {
      return *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA__PREV_OFFSET_));
    }
    //.........................................................................
    inline void SetPrev(long int prev,long int ptr) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&prev,sizeof(long int),__LINE__,__FILE__);
      #endif

      *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__PREV_OFFSET_))=prev;
    }
    //.........................................................................
    inline void SetPrev(long int prev,byte* ParticleDataStart) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&prev,sizeof(long int),__LINE__,__FILE__);
      #endif

      *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA__PREV_OFFSET_))=prev;
    }
    //-------------------------------------------------------------------------

    // Operations related to the particle velocity
    //-------------------------------------------------------------------------
    inline double *GetV(long int ptr) {
      return (double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__VELOCITY_OFFSET_);
    }
    //.........................................................................
    inline double *GetV(byte *ParticleDataStart) {
      return (double*) (ParticleDataStart+_PIC_PARTICLE_DATA__VELOCITY_OFFSET_);
    }
    //.........................................................................
    inline void GetV(double* v,long int ptr) {
      memcpy(v,ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__VELOCITY_OFFSET_,3*sizeof(double));
    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      for (int idim=0;idim<3;idim++) if (isfinite(v[idim])==false) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
    #endif
    #endif
    }
    //.........................................................................
    inline void GetV(double* v,byte *ParticleDataStart) {
      memcpy(v,ParticleDataStart+_PIC_PARTICLE_DATA__VELOCITY_OFFSET_,3*sizeof(double));

    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      for (int idim=0;idim<3;idim++) if (isfinite(v[idim])==false) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
    #endif
    #endif
    }
    //.........................................................................
    inline void SetV(double* v,long int ptr) {
/*      if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>1.0e9) {
        exit(__LINE__,__FILE__,"the velocity is too large");
      }*/

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      for (int idim=0;idim<3;idim++) if (isfinite(v[idim])==false) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(v,3*sizeof(double),__LINE__,__FILE__);
      #endif

      memcpy(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__VELOCITY_OFFSET_,v,3*sizeof(double));
    }
    //.........................................................................
    inline void SetV(double* v,byte *ParticleDataStart) {
/*      if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>1.0e9) {
        exit(__LINE__,__FILE__,"the velocity is too large");
      }*/

    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      for (int idim=0;idim<3;idim++) if (isfinite(v[idim])==false) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
    #endif
    #endif

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(v,3*sizeof(double),__LINE__,__FILE__);
      #endif

      memcpy(ParticleDataStart+_PIC_PARTICLE_DATA__VELOCITY_OFFSET_,v,3*sizeof(double));
    }
    //-------------------------------------------------------------------------

    // Operations related to the particle position
    //-------------------------------------------------------------------------
    inline double *GetX(long int ptr) {
      return (double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__POSITION_OFFSET_);
    }
    //.........................................................................
    inline double *GetX(byte *ParticleDataStart) {
      return (double*) (ParticleDataStart+_PIC_PARTICLE_DATA__POSITION_OFFSET_);
    }
    //.........................................................................
    inline void GetX(double* x,long int ptr) {
      memcpy(x,ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__POSITION_OFFSET_,DIM*sizeof(double));
    }
    //.........................................................................
    inline void GetX(double* x,byte *ParticleDataStart) {
      memcpy(x,ParticleDataStart+_PIC_PARTICLE_DATA__POSITION_OFFSET_,DIM*sizeof(double));
    }
    //.........................................................................
    inline void SetX(double* x,long int ptr) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(x,3*sizeof(double),__LINE__,__FILE__);
      #endif

      memcpy(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__POSITION_OFFSET_,x,DIM*sizeof(double));
    }
    //.........................................................................
    inline void SetX(double* x,byte *ParticleDataStart) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(x,3*sizeof(double),__LINE__,__FILE__);
      #endif

      memcpy(ParticleDataStart+_PIC_PARTICLE_DATA__POSITION_OFFSET_,x,DIM*sizeof(double));
    }
    //-------------------------------------------------------------------------

    // Operations related to the individual particle weight correction    
    //-------------------------------------------------------------------------
    inline double GetIndividualStatWeightCorrection(long int ptr) {
    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      return *((double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_));
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      return 1;
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }
    //.........................................................................
    inline double GetIndividualStatWeightCorrection(byte *ParticleDataStart) {
    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      return *((double*) (ParticleDataStart+_PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_));
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      return 1;
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }
    //.........................................................................
    inline void SetIndividualStatWeightCorrection(double WeightCorrectionFactor,long int ptr) {
    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(WeightCorrectionFactor)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
    #endif
    #endif

    #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
    PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&WeightCorrectionFactor,sizeof(double),__LINE__,__FILE__);
    #endif

    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      *((double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_)) =WeightCorrectionFactor;
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      exit(__LINE__,__FILE__,"Error: SetIndividualStatWeightCorrection cannot be used with _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_");
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }
    //.........................................................................
    inline void SetIndividualStatWeightCorrection(double WeightCorrectionFactor,byte *ParticleDataStart) {
    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(WeightCorrectionFactor)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
    #endif
    #endif

    #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
    PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&WeightCorrectionFactor,sizeof(double),__LINE__,__FILE__);
    #endif

    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      *((double*) (ParticleDataStart+_PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_)) =WeightCorrectionFactor;
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      exit(__LINE__,__FILE__,"Error: SetIndividualStatWeightCorrection cannot be used with _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_");
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }
    //-------------------------------------------------------------------------

    // Operations related to the particle magnetic moment
    //-------------------------------------------------------------------------
    inline double GetMagneticMoment(long int ptr) {
      return *(double*)(ParticleDataBuffer + ptr*ParticleDataLength + _PIC_PARTICLE_DATA__MAGNETIC_MOMENT_OFFSET_);
    }

    //.........................................................................
    inline double GetMagneticMoment(byte *ParticleDataStart) {
      return *(double*)(ParticleDataStart + _PIC_PARTICLE_DATA__MAGNETIC_MOMENT_OFFSET_);
    }

    //.........................................................................
    inline void SetMagneticMoment(double MagneticMoment, long int ptr) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&MagneticMoment,sizeof(double),__LINE__,__FILE__);
      #endif

      *(double*)(ParticleDataBuffer +ptr*ParticleDataLength + _PIC_PARTICLE_DATA__MAGNETIC_MOMENT_OFFSET_) = MagneticMoment;
    }

    //.........................................................................
    inline void SetMagneticMoment(double MagneticMoment, byte* ParticleDataStart) {

      #if _PIC_DEBUGGER__SAVE_DATA_STREAM_MODE_ == _PIC_MODE_ON_
      PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(&MagneticMoment,sizeof(double),__LINE__,__FILE__);
      #endif

      *(double*)(ParticleDataStart + _PIC_PARTICLE_DATA__MAGNETIC_MOMENT_OFFSET_) = MagneticMoment;
    }
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    //save/read the face number at which the particle was injected into the simulation
    inline void SetInjectionFaceNumber(int nface,byte *ParticleDataStart) {
      if (_USE_SAVE_INJECTION_FACE_ == _PIC_MODE_ON_) *((int*) (ParticleDataStart+_PIC_PARTICLE_DATA__INJECTION_FACE_OFFSET_))=nface;
    }

    inline int GetInjectionFaceNumber(byte *ParticleDataStart) {
      if (_USE_SAVE_INJECTION_FACE_ == _PIC_MODE_ON_) return *((int*) (ParticleDataStart+_PIC_PARTICLE_DATA__INJECTION_FACE_OFFSET_));

      return -1;
    }
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    //save/read the initial value of the total particle weight over the local time step at the moment of the particle creation
    inline void SetParticleWeightOverTimeStepRatio(double ratio,byte *ParticleDataStart) {
      if (_USE_SAVE_PARTICLE_WEIGHT_OVER_LOCAL_TIME_STEP_==_PIC_MODE_ON_) *((double*) (ParticleDataStart+_PIC_PARTICLE_DATA__WEIGHT_OVER_TIME_STEP_OFFSET_))=ratio;
    }

    inline double GetParticleWeightOverTimeStepRatio(byte *ParticleDataStart) {
      if (_USE_SAVE_PARTICLE_WEIGHT_OVER_LOCAL_TIME_STEP_==_PIC_MODE_ON_) return *((double*) (ParticleDataStart+_PIC_PARTICLE_DATA__WEIGHT_OVER_TIME_STEP_OFFSET_));

      return 0.0;
    }
    //-------------------------------------------------------------------------


    // Operations related to the field line particle sits on                    
    //-------------------------------------------------------------------------
    inline int GetFieldLineId(long int ptr) {
      return *(int*)(ParticleDataBuffer + ptr*ParticleDataLength + _PIC_PARTICLE_DATA__FIELD_LINE_ID_OFFSET_);
    }

    //.........................................................................
    inline int GetFieldLineId(byte *ParticleDataStart) {
      return *(int*)(ParticleDataStart + _PIC_PARTICLE_DATA__FIELD_LINE_ID_OFFSET_);
    }

    //.........................................................................
    inline void SetFieldLineId(int FieldLineId, long int ptr) {
      *(int*)(ParticleDataBuffer +ptr*ParticleDataLength + _PIC_PARTICLE_DATA__FIELD_LINE_ID_OFFSET_) = FieldLineId;
    }

    //.........................................................................
    inline void SetFieldLineId(int FieldLineId, byte* ParticleDataStart) {
      *(int*)(ParticleDataStart + _PIC_PARTICLE_DATA__FIELD_LINE_ID_OFFSET_) = FieldLineId;
    }
    //.........................................................................
    inline double GetFieldLineCoord(long int ptr) {
      return *(double*)(ParticleDataBuffer + ptr*ParticleDataLength +_PIC_PARTICLE_DATA__FIELD_LINE_COORD_OFFSET_);
    }

    //.........................................................................
    inline double GetFieldLineCoord(byte *ParticleDataStart) {
      return *(double*)(ParticleDataStart +_PIC_PARTICLE_DATA__FIELD_LINE_COORD_OFFSET_);
    }

    //.........................................................................
    inline void SetFieldLineCoord(double FieldLineCoord, long int ptr) {
      *(double*)(ParticleDataBuffer +ptr*ParticleDataLength +_PIC_PARTICLE_DATA__FIELD_LINE_COORD_OFFSET_)= FieldLineCoord;
    }

    //.........................................................................
    inline void SetFieldLineCoord(double FieldLineCoord,byte* ParticleDataStart){
      *(double*)(ParticleDataStart +_PIC_PARTICLE_DATA__FIELD_LINE_COORD_OFFSET_)= FieldLineCoord;
    }

    //-------------------------------------------------------------------------

    inline bool IsParticleAllocated(byte* ParticleDataStart) {
      unsigned char flag;

      flag=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x80);
      return (flag==0) ? false : true;
    }

    inline bool IsParticleAllocated(long int ptr) {
      unsigned char flag;

      flag=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x80);
      return (flag==0) ? false : true;
    }

    inline void SetParticleDeleted(byte* ParticleDataStart) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x7f);
      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))=t;
    }

    inline void SetParticleDeleted(long int ptr) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x7f);
      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))=t;
    }

    inline void SetParticleAllocated(byte* ParticleDataStart) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x7f);
      t|=0x80;

      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))=t;
    }

    inline void SetParticleAllocated(long int ptr) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))) & 0x7f);
      t|=0x80;

      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_))=t;
    }


    //========================================================

    //the particle buffer procedure
    void Init(long int);
    long int GetMaxNPart();
    long int GetAllPartNum();
    long int GetParticleDataLength();

    long int GetNewParticle(bool RandomThreadOpenMP=false);
    long int GetNewParticle(long int&,bool RandomThreadOpenMP=false);

    /*DeleteParticle_withoutTrajectoryTermination() acts as  DeleteParticle() when _PIC_PARTICLE_TRACKER_MODE_  == _PIC_MODE_OFF_;
     if _PIC_PARTICLE_TRACKER_MODE_  == _PIC_MODE_ON_ DeleteParticle_withoutTrajectoryTermination() does not terminate sampling of the particle trajectory; the function should be used only
     from PIC::Parallel::ExchangeParticleData() when particles are moved between processors
    */
    void DeleteParticle(long int);
    void DeleteParticle(long int,long int&);
    void DeleteParticle_withoutTrajectoryTermination(long int,bool RandomThreadOpenMP=false);

    void CloneParticle(long int,long int);
    void CloneParticle(byte*,byte*);

    void ExcludeParticleFromList(long int,long int&);

    void SaveImageFile(int);
    void LoadImageFile(int);

    void PackParticleData(char*,long int);
    void UnPackParticleData(char*,long int);

    unsigned long GetChecksum();
    unsigned long GetChecksum(const char *msg);
    unsigned long GetChecksum(int nline,const char *fname);
    unsigned long GetChecksum(int code,int nline,const char *fname);
  }



  namespace Mesh {
    class cDataCenterNode;

    // aliases
    typedef PIC::Datum::cDatumTimed    cDatumTimed;
    typedef PIC::Datum::cDatumWeighted cDatumWeighted;

    // the following class isn't sampled directly
    // values are derived from sampled variables
    //-------------------------------------------------------------------------
    class cDatumDerived : public Datum::cDatumSampled {
    public:
      // function to find an averaged value derived from other variables;
      // CURRENTLY WORKS ONLY FOR CELL-CENTERED DATA
      void (cDataCenterNode::*GetAverage)(double*, int);
      //.......................................................................
      inline bool is_active(){return GetAverage != NULL;}

      inline void activate(void (cDataCenterNode::*GetAverageIn)(double*, int),vector<cDatumDerived*>* DatumVector) {
	      if (is_active()) exit(__LINE__,__FILE__,"ERROR: trying to activate datum a second time");

	      GetAverage = GetAverageIn;
	      // add this datum to the provided cDatumSampled vector
	      DatumVector->push_back(this);
      }

      // constructor is inherited as well
      //.......................................................................
      cDatumDerived(int lengthIn, const char* nameIn, bool doPrintIn = true) : Datum::cDatumSampled(lengthIn, nameIn, doPrintIn) {
	      type = Derived_; GetAverage=NULL;
      }
    };

    // class cDatumDerived ----------------------------------------------------
    //vector of active sampling data
    extern vector<Datum::cDatumSampled*> DataSampledCenterNodeActive;
    //vector of active derived data
    extern vector<cDatumDerived*> DataDerivedCenterNodeActive;
    
    //basic macroscopic parameters sampled in the simulation
    extern cDatumTimed    DatumParticleWeight;
    extern cDatumTimed    DatumParticleNumber;
    extern cDatumTimed    DatumNumberDensity;
    extern cDatumWeighted DatumParticleVelocity;
    extern cDatumWeighted DatumParticleVelocity2;
    extern cDatumWeighted DatumParticleSpeed;
    extern cDatumWeighted DatumParticleParallelVelocity;
    extern cDatumWeighted DatumParticleParallelVelocity2;
    extern cDatumDerived  DatumTranslationalTemperature;
    extern cDatumDerived  DatumParallelTranslationalTemperature;
    extern cDatumDerived  DatumTangentialTranslationalTemperature;

    //the limiting size of the domain 
    extern double xmin[3],xmax[3];
    
    //the function controlling the local mesh resolution
    typedef double (*fLocalMeshResolution) (double*);
    extern fLocalMeshResolution LocalMeshResolution;

    //the offset of the sampled infomation that is stored in 'center nodes'
    extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;
    
    //the data and order that the data are saved in the associated data buffer of 'center nodes'
    //3. The offset of the data buffer for 'completed sample'
    //4. The offset of the data buffer for 'collecting sample'
    
    //sampling data each specie
    //a. for each specie
    //1. particle weight
    //2. particle number
    //3. particle velocity[3]
    //4. particle pow(velocity[3],2)
    //b. sampling data requested for involved physical models and external species
    extern int sampledParticleWeghtRelativeOffset,sampledParticleNumberRelativeOffset,sampledParticleNumberDensityRelativeOffset;
    extern int sampledParticleVelocityRelativeOffset,sampledParticleVelocity2RelativeOffset,sampledParticleSpeedRelativeOffset;
    extern int sampledParticleNormalParallelVelocityRelativeOffset,sampledParticleNormalParallelVelocity2RelativeOffset;
    extern int sampledExternalDataRelativeOffset;
    extern int sampleSetDataLength;
    
    
    
    //user defiend functions for printing the 'center node' data into an output file
    typedef void (*fPrintVariableListCenterNode)(FILE* fout,int DataSetNumber);
    typedef void (*fPrintDataCenterNode)(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,cDataCenterNode *CenterNode);
    typedef void (*fInterpolateCenterNode)(cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,cDataCenterNode *CenterNode);
    
    extern vector<fPrintVariableListCenterNode> PrintVariableListCenterNode;
    extern vector<fPrintDataCenterNode> PrintDataCenterNode;
    extern vector<fInterpolateCenterNode> InterpolateCenterNode;

    void AddVaraibleListFunction(fPrintVariableListCenterNode f);
    
    //the class defining the 'central node' that contains the sampling data
    class cDataCenterNode : public cBasicCenterNode {
    public:
      //parameters that defines the parameters of the associated data used for sampling and code running
      static int totalAssociatedDataLength,LocalParticleVolumeInjectionRateOffset;
      
      //      long int FirstCellParticle,tempParticleMovingList;
      
      char *associatedDataPointer;
      bool ActiveFlag;
      
      inline int AssociatedDataLength() {
        return totalAssociatedDataLength;
      }
      
      void SetAssociatedDataBufferPointer(char* ptr) {
        associatedDataPointer=ptr;
      }
      
      inline char* GetAssociatedDataBufferPointer() {
        return associatedDataPointer;
      }
          
      //clean the sampling buffers
      void cleanDataBuffer() {
        cBasicCenterNode::cleanDataBuffer();
	
        int i,length=totalAssociatedDataLength/sizeof(double);
        double *ptr;
        for (i=0,ptr=(double*)associatedDataPointer;i<length;i++,ptr++) *ptr=0.0;
	
        if (totalAssociatedDataLength%sizeof(double)) exit(__LINE__,__FILE__,"Error: the cell internal buffers contains data different from double");
      }
      
      //init the buffers
      cDataCenterNode() : cBasicCenterNode() {
        associatedDataPointer=NULL;
        ActiveFlag=false;
      }
        
      // access sampled macroscopic parameters;
      // some function have 2 different interfaces:
      // - for array
      // - for scalars
      //-----------------------------------------------------------------------

      void SampleDatum(Datum::cDatumSampled Datum, double* In, int spec, double weight=1.0);
      void SampleDatum(Datum::cDatumSampled Datum, double In, int spec,  double weight=1.0);
      void SetDatum(Datum::cDatumSampled Datum, double* In, int spec);
      void GetDatumCumulative(Datum::cDatumSampled Datum, double* Out, int spec);
      double GetDatumCumulative(Datum::cDatumSampled Datum, int spec);
      void GetDatumAverage(cDatumTimed Datum, double* Out, int spec);
      double GetDatumAverage(cDatumTimed Datum, int spec);
      void GetDatumAverage(cDatumWeighted Datum, double* Out, int spec);
      double GetDatumAverage(cDatumWeighted Datum, int spec);
      void GetDatumAverage(cDatumDerived Datum, double* Out, int spec);

      //backward compatible access
      double GetNumberDensity(int spec);
      void GetBulkVelocity(double* vOut, int spec);
      double GetMeanParticleSpeed(int spec);
      double GetCompleteSampleCellParticleWeight(int spec);

      // data interpolation
      void InterpolateDatum(Datum::cDatumSampled Datum, cDataCenterNode** InterpolationList,double *InterpolationCoefficients,int nInterpolationCoefficients, int spec);


      /*
      inline void SampleDatum(Datum::cDatumSampled Datum, double* In, int spec, double weight=1.0) {
	      for(int i=0; i<Datum.length; i++) {
	        *(i + Datum.length * spec + (double*)(associatedDataPointer + collectingCellSampleDataPointerOffset+Datum.offset))+= In[i] * weight;
	      }
      }

      inline void SampleDatum(Datum::cDatumSampled Datum, double In, int spec,  double weight=1.0) {
	      *(spec + (double*)(associatedDataPointer + collectingCellSampleDataPointerOffset+Datum.offset))+= In * weight;
      }

      //.......................................................................
      inline void SetDatum(Datum::cDatumSampled Datum, double* In, int spec) {
        for(int i=0; i<Datum.length; i++) *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) = In[i];
      }

      //get accumulated data
      //.......................................................................
      inline void GetDatumCumulative(Datum::cDatumSampled Datum, double* Out, int spec) {
	      for(int i=0; i<Datum.length; i++) Out[i] = *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset));
      }

      inline double GetDatumCumulative(Datum::cDatumSampled Datum, int spec) {
	      return *(spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset));
      }

      //get data averaged over time
      //.......................................................................
      inline void GetDatumAverage(cDatumTimed Datum, double* Out, int spec) {
	      if (PIC::LastSampleLength > 0) for (int i=0; i<Datum.length; i++) {
	        Out[i] = *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / PIC::LastSampleLength;
	      }
	      else for (int i=0; i<Datum.length; i++) Out[i] = 0.0;
      }

      inline double GetDatumAverage(cDatumTimed Datum, int spec) {
	      return (PIC::LastSampleLength > 0) ?  *(spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / PIC::LastSampleLength : 0.0;
      }

      //get data averaged over sampled weight
      //.......................................................................
      inline void GetDatumAverage(cDatumWeighted Datum, double* Out, int spec) {
	      double TotalWeight=0.0;

	      GetDatumCumulative(DatumParticleWeight, &TotalWeight, spec);

	      if (TotalWeight > 0) {
	        for(int i=0; i<Datum.length; i++) Out[i] = *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / TotalWeight;
	      }
	      else for(int i=0; i<Datum.length; i++) Out[i] = 0.0;
      }

      inline double GetDatumAverage(cDatumWeighted Datum, int spec) {
	      double TotalWeight=0.0;

	      GetDatumCumulative(DatumParticleWeight, &TotalWeight, spec);
	      return (TotalWeight > 0) ? *(spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / TotalWeight : 0.0;
      }

      //get average for derived data
      //.......................................................................
      inline void GetDatumAverage(cDatumDerived Datum, double* Out, int spec) {
        (this->*Datum.GetAverage)(Out,spec);
      }
      //-----------------------------------------------------------------------

      //backward compatible access
      //-----------------------------------------------------------------------
      inline double GetNumberDensity(int spec) {
	      return GetDatumAverage(DatumNumberDensity, spec);
      }

      inline void GetBulkVelocity(double* vOut, int spec) {
	      GetDatumAverage(DatumParticleVelocity, vOut, spec);
      }

      inline double GetMeanParticleSpeed(int spec) {
	      return GetDatumAverage(DatumParticleSpeed, spec);
      }

      inline double GetCompleteSampleCellParticleWeight(int spec) {
	      return GetDatumAverage(DatumParticleWeight, spec);
      }
      //-----------------------------------------------------------------------

      // data interpolation
      //-----------------------------------------------------------------------
      inline void InterpolateDatum(Datum::cDatumSampled Datum, cDataCenterNode** InterpolationList,double *InterpolationCoefficients,int nInterpolationCoefficients, int spec) {
	      // container for the interpolated value; set it to be zero
	      double value[Datum.length], interpolated[Datum.length];

	      for (int i=0; i<Datum.length; i++) {
	        value[i]=0.0; interpolated[i]=0.0;
	      }

	      // interpolation loop
	      for (int i=0; i<nInterpolationCoefficients; i++) {
	        InterpolationList[i]->GetDatumCumulative(Datum, value, spec);

	        for(int j=0; j<Datum.length; j++) interpolated[j] += InterpolationCoefficients[i] * value[j];
	      }

	      SetDatum(Datum, interpolated, spec);

        #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
        #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
	        for(int i=0; i<Datum.length; i++) if (isfinite(interpolated[i])==false) exit(__LINE__,__FILE__,"Error: Floating Point Exception");
        #endif
        #endif
      }
      //-----------------------------------------------------------------------
*/
    
      //get the sampled macroscopic parameter of the flow (for derived data)
      //-----------------------------------------------------------------------
      void GetTranslationalTemperature(double* T, int s) {
        int idim;
        double res=0.0,v[3]={0.0,0.0,0.0},v2[3]={0.0,0.0,0.0};

        GetDatumAverage(DatumParticleVelocity, v, s);
        GetDatumAverage(DatumParticleVelocity2,v2,s);

        for (idim=0;idim<3;idim++) res+=v2[idim]-v[idim]*v[idim];

        //return 
	      T[0]=PIC::MolecularData::GetMass(s)*res/(3.0*Kbol);
      }

      inline double GetTranslationalTemperature(int s) {
        double res=0.0;

        GetTranslationalTemperature(&res, s);
        return res;
      }

      void GetParallelTranslationalTemperature(double* T, int s) {
        #if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
        //return 
	      GetTranslationalTemperature(T, s);
        #else
        double v=0.0,v2=0.0,w=0.0,res=0.0;

        v = GetDatumAverage(DatumParticleParallelVelocity, s);
	      v2 = GetDatumAverage(DatumParticleParallelVelocity2, s);
	      res=PIC::MolecularData::GetMass(s)*(v2-v*v)/Kbol;
	      T[0]=res;
        #endif
      }

      void GetTangentialTranslationalTemperature(double* T, int s) {
        #if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
	      GetTranslationalTemperature(T,s);
        #else
	      double t,tp;


	      GetTranslationalTemperature(&t,s);
	      GetParallelTranslationalTemperature(&tp,s);
	      T[0]=0.5*(3.0*t-tp);
        #endif
      }
      //-----------------------------------------------------------------------

      //print the sampled data into a file
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread);

      void PrintFileDescriptior(FILE* fout,int DataSetNumber) {
        char sym[_MAX_STRING_LENGTH_PIC_];

        PIC::MolecularData::GetChemSymbol(sym,DataSetNumber);
        fprintf(fout,"TITLE=\"specie=%s\"",sym);
      }


      void PrintVariableList(FILE* fout,int DataSetNumber);
      void Interpolate(cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients);
    };

  
    //the class that contains the run information for the cell's corners
    class cDataCornerNode : public cBasicCornerNode {
    public:
      bool ActiveFlag;

      cDataCornerNode() {
        ActiveFlag=false;
      }
    };
  
    //the data stored in a block
    //1. Local Time Step [NS]: depending on the model mode there will be a 'global' time step for the simulation, 'global' time step for the cell or individual time step for each simulated species
    //2. Local particle weight [NS]: depending on the model mode there will be a 'global' weight for the simulation, 'global' weight for the cell or individual weigh for each simulated species

    class cDataBlockAMR : public cBasicBlockAMR<cDataCornerNode,cDataCenterNode> {
    public:
      static int LocalTimeStepOffset,LocalParticleWeightOffset;
      char *associatedDataPointer;

    private:
      static int tempParticleMovingListTableThreadOffset,tempParticleMovingListTableThreadLength; //the offset and length of the tempParticleMovingListTable for each
      static int totalAssociatedDataLength;

    public:
      static int LoadBalancingMeasureOffset;
      static int UserAssociatedDataOffset;
      long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

      bool ActiveFlag; //used for debugging to prevent repeatable de-allocation of the block

      //the flag defined whether all internal data is allocated
      static bool InternalDataInitFlag;


      //get pointer to element of tempParticleMovingListTableThread when OpenMP is in use
      long int *GetTempParticleMovingListTableThread(int thread,int i,int j,int k) {
        return ((long int*) (associatedDataPointer+tempParticleMovingListTableThreadOffset))+
            thread*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
      }

      //tempParticleMovingListTable is used only when _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
      //when _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_ tempParticleMovingListTableThreadOffset is used instead
      #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
      long int tempParticleMovingListTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
      #endif

      int AssociatedDataLength() {
        return totalAssociatedDataLength;
      }

      void SetAssociatedDataBufferPointer(char* ptr) {
        associatedDataPointer=ptr;
      }


      char* GetAssociatedDataBufferPointer() {
        return associatedDataPointer;
      }


      static void InitInternalData() {
        if (InternalDataInitFlag==false) {
          totalAssociatedDataLength=0;

          //init the tempParticleMovingListTable when OpenMP is in use
          #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
          int nThreadsOpenMP;

          #pragma omp parallel default (none) shared(nThreadsOpenMP)
          {
            #pragma omp single
            {
              nThreadsOpenMP=omp_get_num_threads();
            }
          }

          tempParticleMovingListTableThreadOffset=totalAssociatedDataLength;
          tempParticleMovingListTableThreadLength=nThreadsOpenMP*sizeof(long int)*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;

          totalAssociatedDataLength+=tempParticleMovingListTableThreadLength;
          UserAssociatedDataOffset+=tempParticleMovingListTableThreadLength;

          LoadBalancingMeasureOffset=totalAssociatedDataLength;
          totalAssociatedDataLength+=nThreadsOpenMP*sizeof(double);
          UserAssociatedDataOffset+=nThreadsOpenMP*sizeof(double);
          #endif //_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_

          //set up the init flag
          InternalDataInitFlag=true;
        }
      }

      static int RequestInternalBlockData(int length) {
        if (InternalDataInitFlag==false) InitInternalData();

        int res=totalAssociatedDataLength;
        totalAssociatedDataLength+=length;

        return res;
      }


      cDataBlockAMR () : cBasicBlockAMR<cDataCornerNode,cDataCenterNode> () {
        if (InternalDataInitFlag==false) InitInternalData();
        ActiveFlag=false;
      }

    
      //exchenge of the data between processors
      void sendBoundaryLayerBlockData(CMPI_channel *pipe);
      void recvBoundaryLayerBlockData(CMPI_channel *pipe,int From);

      //send the block to abother processor
      void sendMoveBlockAnotherProcessor(CMPI_channel *pipe);
      void recvMoveBlockAnotherProcessor(CMPI_channel *pipe,int From);

      //clean the sampling buffers
      void cleanDataBuffer() {
        int i,length=(totalAssociatedDataLength-UserAssociatedDataOffset)/sizeof(double);
        double *ptr;

        //clean the associated data buffers
        for (i=0,ptr=(double*)(associatedDataPointer+UserAssociatedDataOffset);i<length;i++,ptr++) *ptr=0.0;

        //clean the base class' data
        cBasicBlockAMR<cDataCornerNode,cDataCenterNode>::cleanDataBuffer();

        //clean the Particle Tables
        length=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
        for (i=0;i<length;i++) FirstCellParticleTable[i]=-1;

 #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
        for (i=0;i<length;i++) tempParticleMovingListTable[i]=-1;
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
        length*=omp_get_max_threads();
        for (i=0;i<length;i++) *(i+(long int *)(associatedDataPointer+tempParticleMovingListTableThreadOffset))=-1;
#else
#error The option is unknown
#endif
      }


      //set and get the local time step
      void SetLocalTimeStep(double dt, int spec);
      double GetLocalTimeStep(int spec);
      void SetLocalParticleWeight(double weight, int spec);
      double GetLocalParticleWeight(int spec);

      //print into a output file the blocks' parameters: the local time step, the local weight
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int BlockThread) {

        struct cOutputData {
          double dtLocal,wLocal;
        } OutputData;



        if (pipe->ThisThread==BlockThread) {
          OutputData.dtLocal=GetLocalTimeStep(DataSetNumber);
          OutputData.wLocal=GetLocalParticleWeight(DataSetNumber);
        }

        if (pipe->ThisThread==0) {
           if (BlockThread!=0) pipe->recv((char*)&OutputData,sizeof(OutputData),BlockThread);

           fprintf(fout,"%e  %e  ",OutputData.dtLocal,OutputData.wLocal);
         }
         else pipe->send((char*)&OutputData,sizeof(OutputData));
      }


      void PrintVariableList(FILE* fout) {
        fprintf(fout,", \"Local Time Step\" ");
        fprintf(fout,", \"log10(Local Particle Weight)\" ");
      }


    };

  
    //init the sampling buffers of the cell's data
    void initCellSamplingDataBuffer();
    void SetCellSamplingDataRequest();

    //return time step and the particle's weights
    inline double GetLocalTimeStep(int spec,cDataBlockAMR* block) {
      return *(spec+(double*)(cDataBlockAMR::LocalTimeStepOffset+block->GetAssociatedDataBufferPointer()));
    }

    inline void SetLocalTimeStep(double dt,int spec,cDataBlockAMR* block) {
      *(spec+(double*)(cDataBlockAMR::LocalTimeStepOffset+block->GetAssociatedDataBufferPointer()))=dt;
    }

    double GetLocalParticleWeight(int,cDataBlockAMR*);
    void SetLocalParticleWeight(double,int,cDataBlockAMR*);

    void flushCompletedSamplingBuffer(cDataCenterNode*);
    void flushCollectingSamplingBuffer(cDataCenterNode*);
    void switchSamplingBuffers();

    //the computational mesh
    #if DIM == 3
    extern cMeshAMR3d<cDataCornerNode,cDataCenterNode,cDataBlockAMR > mesh;
    #elif DIM == 2
    extern cMeshAMR2d<cDataCornerNode,cDataCenterNode,cDataBlockAMR > mesh;
    #else
    extern cMeshAMR1d<cDataCornerNode,cDataCenterNode,cDataBlockAMR > mesh;
    #endif

    //init the computational mesh
    void Init(double*,double*,fLocalMeshResolution);
    void buildMesh();
    void loadMesh(char*);

    //tratment of the cut-cells
    namespace IrregularSurface {
      using namespace CutCell;

      //propagate the information of the cut faces to the neibbouring nodes
      extern int nCutFaceInformationCopyAttempts;
      void CopyCutFaceInformation(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=PIC::Mesh::mesh.rootTree);

      //init the vectors of the external normals at the cut-faces
      void InitExternalNormalVector();
      bool CheckPointInsideDomain_default(double *x,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,bool ParallelCheck,double EPS);

      //determine the closes distance to the triangulated surface
      //if the point is inside the body -> the distance is -1
      double GetClosestDistance(double *x,bool ParallelExecutionModeMPI=false);
      double GetClosestDistance(double *x,double *xClosestPoint,int& iClosestTriangularFace,bool ParallelExecutionModeMPI=false);

    }

    //the namespace for the block/cell search functions
    //-------------------------------------------------------------------------
    namespace Search {
      extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* ***HashTable;
      extern double dx[3],xMinDomain[3]; //left corner of the domain, the interval in each dimension

      //the resolution level of the hash table
      const int HashTableLevel=7;

      //set up the hash table
      void Init();
      void ScanTree(int imin,int jmin,int kmin,int IndexInterval,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

      //find cell nad block
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *FindBlock(double *x);
      PIC::Mesh::cDataCenterNode *FindCell(double *x);

    }
    //namespace Search --------------------------------------------------------

    //get the interpolation stencil (used only when the lineat interpolation is set)
    int GetCenterNodesInterpolationCoefficients(double *x,double *CoefficientsList,PIC::Mesh::cDataCenterNode **InterpolationStencil,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,int nMaxCoefficients);
  } 
  // namespace Mesh ===========================================================

  //volume injection of model particles
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  namespace VolumeParticleInjection {

    //init the model
    void Init();

    //Generte new particle internal properties
    typedef void (*fGenerateInternalParticleProperties)(long int ptr,int spec,int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    extern fGenerateInternalParticleProperties GenerateInternalParticleProperties;

    //generate a random position in a cell
    void GetRandomCellPosition(double *x,int iCell,int jCell,int kCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //The list of volume injection processes
    //the injection rate of particle due to a specific injection process
    typedef void (*fSpeciesInjectionRate)(bool *InjectionFlag,double *Rate, int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //the function that process that injection reaction (generate the particles)
    typedef long int (*fInjectionProcessor)(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //the limit of the local time step due to the injection process
    typedef double (*fLocalTimeStepLimit)(int spec,bool& TimeStepLimitationImposed, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //the array that stores the total injection rate by the volume injection
    extern double *SourceRate;

    struct cVolumeInjectionDescriptor {
      fSpeciesInjectionRate SpeciesInjectionRate;
      fInjectionProcessor InjectionProcessor;
      fLocalTimeStepLimit LocalTimeStepLimit;
    };

    const int nMaxInjectionProcessEntries=128;
    extern int nRegistratedInjectionProcesses;

    extern cVolumeInjectionDescriptor VolumeInjectionDescriptor[nMaxInjectionProcessEntries];

    void inline RegisterVolumeInjectionProcess(fSpeciesInjectionRate f1, fInjectionProcessor f2, fLocalTimeStepLimit f3) {
      if (nRegistratedInjectionProcesses==nMaxInjectionProcessEntries-1) {
        exit(__LINE__,__FILE__,"Error: the volume injection processes buffer is overflow, increase the length of the buffer ('PIC::VolumeParticleInjection::nMaxInjectionProcessEntries'");
      }

      if ((f1==NULL)||(f2==NULL)||(f3==NULL)) exit(__LINE__,__FILE__,"Error: one of the functions is not defined");

      if (SourceRate==NULL) {
        if (PIC::nTotalSpecies==0) exit(__LINE__,__FILE__,"Error: neede to know the total number of species before initializing the sampling buffer");

        SourceRate=new double[PIC::nTotalSpecies];
        for (int s=0;s<PIC::nTotalSpecies;s++) SourceRate[s]=0.0;
      }

      VolumeInjectionDescriptor[nRegistratedInjectionProcesses].SpeciesInjectionRate=f1;
      VolumeInjectionDescriptor[nRegistratedInjectionProcesses].InjectionProcessor=f2;
      VolumeInjectionDescriptor[nRegistratedInjectionProcesses].LocalTimeStepLimit=f3;

      nRegistratedInjectionProcesses++;
    }


    long int  InjectParticle();

    //paericle weight injection rates
    void InitTotalInjectionRate();
    double GetTotalInjectionRate(int);
    double GetTotalTimeStepInjection(int spec);
    double GetBlockInjectionRate(int spec,PIC::Mesh::cDataBlockAMR *block);
    double GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell);
  }
#endif

  //sampling functions
  namespace Sampling {

    //the minimum number of iteration for output of the datafile
    extern int minIterationNumberForDataOutput;

    //'SaveOutputDataFile' determines weather the data file will be created for a particular species (output of data files can be suppress for particular species, such as external, dust....)
    extern bool *SaveOutputDataFile;

    //sample the normal and tangential kinetic temperatures: constant origin of the direction of the normal
    static const double constNormalDirection__SampleParallelTangentialTemperature[3]={0.0,0.0,0.0};

    //the number of the first output file that is printed
    static const int FirstPrintedOutputFile=-1;

    namespace ExternalSamplingLocalVariables {

      //the external procedures for sampling particle data
      typedef void (*fSamplingProcessor) ();

      //the procedure that prints the sampled data into a file
      typedef void (*fPrintOutputFile)(int);

      extern const int nMaxSamplingRoutines;
      extern int SamplingRoutinesRegistrationCounter;

      extern fSamplingProcessor *SamplingProcessor;
      extern fPrintOutputFile *PrintOutputFile;

      inline void RegisterSamplingRoutine(fSamplingProcessor p,fPrintOutputFile f) {
        if (SamplingRoutinesRegistrationCounter==0) {
          SamplingProcessor=new fSamplingProcessor[nMaxSamplingRoutines];
          PrintOutputFile=new fPrintOutputFile[nMaxSamplingRoutines];

          for (int i=0;i<nMaxSamplingRoutines;i++) SamplingProcessor[i]=NULL,PrintOutputFile[i]=NULL;
        }
        else if (SamplingRoutinesRegistrationCounter==nMaxSamplingRoutines-1) {
          exit(__LINE__,__FILE__,"Error: SamplingRoutinesRegistrationCounter exeeds its maximum value: increse the value of PIC::Sampling::ExternalSamplingProcedures::nMaxSamplingRoutines");
        }

        SamplingProcessor[SamplingRoutinesRegistrationCounter]=p;
        PrintOutputFile[SamplingRoutinesRegistrationCounter]=f;

        SamplingRoutinesRegistrationCounter++;
      }


    }



    void Sampling();
    void CatchOutLimitSampledValue();

    //sample the particle data
    inline void SampleParticleData(char* ParticleData,char *SamplingBuffer,PIC::Mesh::cDataBlockAMR *block,PIC::Mesh::cDataCenterNode *cell,double TimeStepFraction) {
      double Speed2,*v,LocalParticleWeight,v2;
      int s,idim;
      double *sampledVelocityOffset,*sampledVelocity2Offset;


      Speed2=0.0;

      s=PIC::ParticleBuffer::GetI((PIC::ParticleBuffer::byte*)ParticleData);
      v=PIC::ParticleBuffer::GetV((PIC::ParticleBuffer::byte*)ParticleData);

      LocalParticleWeight=block->GetLocalParticleWeight(s);
      LocalParticleWeight*=TimeStepFraction*PIC::ParticleBuffer::GetIndividualStatWeightCorrection((PIC::ParticleBuffer::byte*)ParticleData);

      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleWeghtRelativeOffset))+=LocalParticleWeight;
      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberRelativeOffset))+=1;
      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberDensityRelativeOffset))+=LocalParticleWeight/cell->Measure;


      sampledVelocityOffset=3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocityRelativeOffset);
      sampledVelocity2Offset=3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocity2RelativeOffset);

      for (idim=0;idim<3;idim++) {
        v2=v[idim]*v[idim];
        Speed2+=v2;

        *(idim+sampledVelocityOffset)+=v[idim]*LocalParticleWeight;
        *(idim+sampledVelocity2Offset)+=v2*LocalParticleWeight;
      }

      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleSpeedRelativeOffset))+=sqrt(Speed2)*LocalParticleWeight;

    }
  }

  //colecular collisions
  namespace MolecularCollisions {

    //collisions between model particles
    namespace ParticleCollisionModel {

      //the limit of the collision number per a single model particles
      extern double CollisionLimitingThrehold;

      //sample collision statistics
      namespace CollsionFrequentcySampling {
        extern int SamplingBufferOffset;

        //offset to the sampling data for collision frequentcy of specie "CollidingSpecies" when it collides with specie "CollisionPartnerSpecies"
        inline int Offset(int CollidingSpecies, int CollisionPartnerSpecies) {return CollidingSpecies*PIC::nTotalSpecies+CollisionPartnerSpecies;}
      }

      //the namespace containes definitions of the user defined functions
      namespace UserDefined {
        double GetTotalCrossSection(double *v0,double *v1,int s0,int s1,PIC::Mesh::cDataBlockAMR* block,PIC::Mesh::cDataCenterNode *cell);
      }

      int RequestSamplingData(int offset);
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void Init();

      //model of the particle collisions
      void ntc();
      void mf();
    }

    //models for calculation of the relative velocity after a collision
    namespace VelocityScattering {
      namespace HS {
        inline void VelocityAfterCollision(double *vrel,int s0,int s1) {
          double Vrc,V[3];
          double CosKsi,SinKsi,CosEps,SinEps,D,c;

          CosKsi=2.0*rnd()-1.0;
          SinKsi=sqrt(1.0-CosKsi*CosKsi);

          c=2*Pi*rnd();
          SinEps=sin(c);
          CosEps=cos(c);

          D=sqrt(vrel[1]*vrel[1]+vrel[2]*vrel[2]);
          if (D>1.0E-6) {
            Vrc=sqrt(vrel[0]*vrel[0]+vrel[1]*vrel[1]+vrel[2]*vrel[2]);
            V[0]=CosKsi*vrel[0]+SinKsi*SinEps*D;
            V[1]=CosKsi*vrel[1]+SinKsi*(Vrc*vrel[2]*CosEps-vrel[0]*vrel[1]*SinEps)/D;
            V[2]=CosKsi*vrel[2]-SinKsi*(Vrc*vrel[1]*CosEps+vrel[0]*vrel[2]*SinEps)/D;
          }
          else {
            V[0]=CosKsi*vrel[0];
            V[1]=SinKsi*CosEps*vrel[0];
            V[2]=SinKsi*SinEps*vrel[0];
          }

          memcpy(vrel,V,3*sizeof(double));
        }
      }
    }

    //collisions with the background atmosphere
#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_
    namespace BackgroundAtmosphere {

      //the total number of the background species, the mass table, the table of cosntant collision cross sections with the model species
      static const int nTotalBackgroundSpecies=1;
      static const double BackgroundSpeciesMassTable[]={0.0};
//      static const double BackgroundAtmosphereConstantCrossSectionTable[PIC::nTotalSpecies][nTotalBackgroundSpecies];
      static const int Background2ModelSpeciesConversionTable[]={-1};

      inline int GetTotalNumberBackgroundSpecies() {return nTotalBackgroundSpecies;}
      inline double GetBackgroundMolecularMass(int spec) {return BackgroundSpeciesMassTable[spec];}


/*
      //the default value of the user-defined function that calculates the collision cross section
      #define _PIC_BACKGROUND_ATMOSPHERE__COLLISION_CROSS_SECTION_FUNCTION_(spec,BackgroundSpecieNumber,modelParticleData,BackgroundAtmosphereParticleData,TranslationalEnergy,cr2) (0.0)

      //get local, cell mean, and cell maximum density of the background species
      #define _PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_SPECIES_LOCAL_DENSITY_(x,BackgroundSpecieNumber,cell,node) (0.0)
      #define _PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_SPECIES_CELL_MEAN_DENSITY_(BackgroundSpecieNumber,cell,node) (0.0)
      #define _PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_SPECIES_CELL_MAXIMUM_DENSITY_(BackgroundSpecieNumber,cell,node) (0.0)

      //the user-defined function for calcualtion of the scattering angle
      #define _PIC_BACKGROUND_ATMOSPHERE__COLLISION_SCATTERING_ANGLE_(Vrel,TranslationalEnergy,spec,BackgroundSpecieNumber) (0.0)

      //the condition the remove the model particle from the simulation after the collision with the background particle
      #define _PIC_BACKGROUND_ATMOSPHERE__REMOVE_CONDITION_MODEL_PARTICLE_(modelParticleData) (true)

      //the conditions to inject the background atmosphere particle after a collision with the model particle
      #define _PIC_BACKGROUND_ATMOSPHERE__INJECT_CONDITION_BACKGROUND_PARTICLE_(BackgroundAtmosphereParticleData) (false)

      //Evaluate GetSigmaCrMax in a cell
      #define _PIC_BACKGROUND_ATMOSPHERE__GET_SIGMA_CR_MAX(spec,BackgroundSpecieNumber,modelParticleData) (0.0)

      //generate the background atmosphere particle
      #define _PIC_BACKGROUND_ATMOSPHERE__GENERATE_BACKGROUND_PARTICLE_(BackgroundAtmosphereParticleData,BackgroundSpecieNumber,cell,node) (exit(__LINE__,__FILE__,"Error: _PIC_BACKGROUND_ATMOSPHERE__GENERATE_BACKGROUND_PARTICLE_ was not set"))

      //define the mode for loading the definition file
      #define _PIC_BACKGROUND_ATMOSPHERE__LOAD_USER_DEFINITION__MODE_ _PIC_MODE_OFF_
      #define _PIC_BACKGROUND_ATMOSPHERE__UDER_DEFINITION_ "UserDefinition.PIC.BackgroundAtmosphere.h"
*/

      //define functions for calculation of the properties of the background atmosphere species
      double GetCollisionCrossSectionBackgoundAtmosphereParticle(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData,PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,double TranslationalEnergy,double cr2);
      double GetSigmaCrMax(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData);
      double GetCollisionScatteringAngle(double* Vrel,double TranslationalEnergy,int spec,int BackgroundSpecieNumber);

      void GenerateBackgoundAtmosphereParticle(PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

      double GetBackgroundLocalNumberDensity(int BackgroundSpecieNumber,double *x);
      double GetCellMeanBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
      double GetCellMaximumBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
      double GetCellLocalBackgroundNumberDensity(double x[3],int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

      //the conditions to keep the background particle or remove a model particle after a collision
      bool KeepConditionModelParticle(PIC::ParticleBuffer::byte *ModelParticleData);

      //use stopping power for determeinig energy loss of simulated species interactinf with the background atmosphere
      double GetStoppingPower(double *x,double *v,int spec,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
      void StoppingPowerProcessor();



      //include the user defined properties of the background atmosphere
      #if _PIC_BACKGROUND_ATMOSPHERE__LOAD_USER_DEFINITION__MODE_ ==  _PIC_MODE_ON_
      #include _PIC_BACKGROUND_ATMOSPHERE__UDER_DEFINITION_
      #endif


      //Sampling of the model data
      extern int LocalEnergyTransferRateSamplingOffset;
      extern int LocalTotalCollisionFreqSamplingOffset;
      int RequestSamplingData(int offset);
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

      //souce and loss rate of the model particles in interactions with the background atmosphere
      extern double **TotalCollisionModelParticleSourceRate,**TotalCollisionModelParticleLossRate;

      //init the model
      void Init_BeforeParser();
      void Init_AfterParser();

      //output sampled model parameters
      void OutputSampledModelData(int);
      void SampleModelData();

      //processor of the collisions
      void CollisionProcessor();
      void RemoveThermalBackgroundParticles();
    }
#endif

  }

  namespace IndividualModelSampling {

    //reserve memory to store sampling data
    typedef int (*fRequestSamplingData)(int);
    extern vector<fRequestSamplingData> RequestSamplingData;

    //reserve memoty in a cell associated data buffer for non-sampling data
    typedef int (*fRequestStaticCellData)(int);
    extern vector<fRequestStaticCellData> RequestStaticCellData;

    //the list of user defined sampling procedures
    typedef void (*fSamplingProcedure)();
    extern vector<fSamplingProcedure> SamplingProcedure;

    //print the variable list
    typedef void (*fPrintVariableList)(FILE* fout,int nDataSet);
    extern vector<fPrintVariableList> PrintVariableList;

    extern vector<PIC::Datum::cDatumSampled*> DataSampledList;

    //interpolate center node data
    typedef void (*fInterpolateCenterNodeData)(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode* cell);
    extern vector<fInterpolateCenterNodeData> InterpolateCenterNodeData;

    //print the sampled node data
    typedef void (*fPrintSampledData)(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode* cell);
    extern vector<fPrintSampledData> PrintSampledData;
  }

  //sample and output the particle's distribution function
  //the coordinates of the sampled locations are set by ampsConfig.pl (#Sampling -> VelocityDistributionSampling)
  namespace DistributionFunctionSample {
    //the init flag
    extern bool SamplingInitializedFlag;

    //the modes for sampling of the v^2 and the absolute value of velocity
    extern const int _LINEAR_SAMPLING_SCALE_,_LOGARITHMIC_SAMPLING_SCALE_;
    extern int v2SamplingMode,speedSamplingMode;

    //the range of the velocity scale and the number of nodes in the sample
    extern double vMax,vMin;
    extern long int nSampledFunctionPoints;
    extern double dV,dV2,dSpeed;

    //the sampling buffers
    extern double **SamplingBuffer;

    //sampling data offsets
    extern int Sample_Velocity_Offset,Sample_Speed_Offset,Sample_V2_Offset,SampleDataLength;

    //get the offset to the beginig of the sampling data for a particular samplePoint, spec,.....
    long int GetSampleDataOffset(int spec,int nInterval);

    //sampling  locations
    extern double SamplingLocations[][3];
    extern int nSamleLocations;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
    extern long int *SampleLocalCellNumber;

    void Init();

    void SampleDistributionFnction();
    void flushSamplingBuffers();

    void printDistributionFunction(char *fname,int spec);
  }

  //sample energy distribution in the relativistic case
  namespace EnergyDistributionSampleRelativistic {
    //the init flag
    extern bool SamplingInitializedFlag;

    //the modes for sampling of the v^2 and the absolute value of velocity
    extern const int _LINEAR_SAMPLING_SCALE_,_LOGARITHMIC_SAMPLING_SCALE_;
    extern int EnergySamplingMode;

    //combine sampled distribution for several species (in the relativistic mode the same chemical species can have a group of the model particle numbers to fix the problem of
    //substantial velocity different of partilces that belong to the same chemical component
    extern vector<vector<int> > CombinedSpeciesDistributionTable;

    //the range of the velocity scale and the number of nodes in the sample
    extern double eMax,eMin,log10eMax,log10eMin;
    extern long int nSampledFunctionPoints;
    extern double dE,log10dE;

    //the sampling buffers
    extern double ***SamplingBuffer;

    //sampling  locations
    extern double SamplingLocations[][3];
    extern int nSamleLocations;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
    extern long int *SampleLocalCellNumber;

    void AddCombinedCombinedParticleDistributionList(vector<int> CombinedSpecies);

    void Init();

    void SampleDistributionFnction();
    void flushSamplingBuffers();

    void printDistributionFunction(char *fname,int spec);
  }

  //sample and output the particle's pitch angle distribution function
  namespace PitchAngleDistributionSample {
    //the init flag
    extern bool SamplingInitializedFlag;

    //the modes for sampling of the v^2 and the absolute value of velocity
    //extern const int _LINEAR_SAMPLING_SCALE_,_LOGARITHMIC_SAMPLING_SCALE_;
    //extern int v2SamplingMode,speedSamplingMode;

    //the range of the sine of pitch angle scale and the number of nodes in the sample
    extern double CosPAMax, CosPAMin;
    extern long int nSampledFunctionPoints;
    extern double dCosPA;

    //the sampling buffers
    extern double **SamplingBuffer;

    //sampling data offsets
    //extern int Sample_Velocity_Offset,Sample_Speed_Offset,Sample_V2_Offset,SampleDataLength;
    extern int Sample_PitchAngle_Offset,SampleDataLength;

    //get the offset to the beginig of the sampling data for a particular samplePoint, spec,.....
    long int GetSampleDataOffset(int spec,int nInterval);

    //sampling  locations
    extern double SamplingLocations[][3];
    extern int nSampleLocations;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
    extern long int *SampleLocalCellNumber;

    void Init();//double ProbeLocations[][DIM],int nProbeLocations);

    void SampleDistributionFnction();
    void flushSamplingBuffers();

    void printDistributionFunction(char *fname,int spec);
  }


  //sample and output the particle's distribution of particle's flux
  namespace ParticleFluxDistributionSample {
    //the init flag
    extern bool SamplingInitializedFlag;

    //the modes for sampling of the v^2 and the absolute value of velocity
    extern const int _LINEAR_SAMPLING_SCALE_,_LOGARITHMIC_SAMPLING_SCALE_;
    extern int v2SamplingMode,speedSamplingMode;

    //the range of the velocity scale and the number of nodes in the sample
    extern double vMax,vMin;
    extern long int nSampledFunctionPoints;
    extern double dV2,dSpeed;

    //the sampling buffers
    extern double **SamplingBuffer;
    extern double **SamplingFlux;

    //sampling data offsets
    extern int Sample_Speed_Offset,Sample_V2_Offset,SampleDataLength;

    //get the offset to the beginig of the sampling data for a particular samplePoint, spec,.....
    long int GetSampleDataOffset(int spec,int nInterval);

    //sampling  locations
    extern double **SamplingLocations;
    extern double **SamplingPointingDirections;
    extern double maxSamplingConeAngle,cosMaxSamplingConeAngle;

    extern int nSamleLocations;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
    extern long int *SampleLocalCellNumber;

    void Init(double ProbeLocations[][DIM],double ProbeDirections[][DIM],double maxConeAngles,int nProbeLocations);

    void SampleDistributionFnction();
    void flushSamplingBuffers();

    void printDistributionFunction(char *fname,int spec);
    void printMacroscopicParameters(char *fname,int spec);
  }


  //procedures for distribution of particle velocities
  namespace Distribution {

    //the macroses defined the modes for generation of the particle parameters
    //_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_ -> the particle properties are generated with the requested distributuion, no weight correction are needed
    //_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_ -> to recover the requasted distribution, a particle weght correction is needed

    #define _PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_         0
    #define _PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_   1


    void MaxwellianVelocityDistribution(double *v,const double *BulkFlowVelocity,const double Temp,const int spec);
    double InjectMaxwellianDistribution(double *v,const double *BulkFlowVelocity,const double Temp,double *ExternalNormal,int spec,int WeightCorrectionMode=_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_);

    void InjectRingDistribution(double *v,double Energy,double *ExternalNormal,int spec);
  }

  //the namespace contained the domain decomposition information used for loop parallelization in OpenMP
  namespace DomainBlockDecomposition {
    extern unsigned int nLocalBlocks;
    extern int LastMeshModificationID;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> **BlockTable;

    void UpdateBlockTable();
  }

  //the mode of the internal degrees of freedom
  namespace IDF {

    static const int nTotalVibtationalModes[]={-1};
    static const int nTotalRotationalModes[]={-1};
    static const int nSpeciesMaxVibrationalModes=0;
    static const double CharacteristicVibrationalTemperature[]={0.0};    //the rule of access CharacteristicVibrationalTemperature[nmode+s*nSpeciesMaxVibrationalModes]
    static const double RotationZnumber[]={0.0};

    extern int _ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_;
    extern int _VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[PIC::nTotalSpecies];
    extern int _TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_;

    namespace LB {
      extern int _ROTATIONAL_ENERGY_OFFSET_,_VIBRATIONAL_ENERGY_OFFSET_;


      inline double GetRotE(PIC::ParticleBuffer::byte *ParticleDataStart) {
        return *((double*)(ParticleDataStart+_ROTATIONAL_ENERGY_OFFSET_));
      }

      inline void SetRotE(double e,PIC::ParticleBuffer::byte *ParticleDataStart) {
         *((double*)(ParticleDataStart+_ROTATIONAL_ENERGY_OFFSET_))=e;
      }

      inline double GetVibE(int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
        if (nmode>=0) return *(nmode+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_));

        double res=0.0;
        int n,nVibModes,s;

        s=PIC::ParticleBuffer::GetI(ParticleDataStart);
        nVibModes=nTotalVibtationalModes[s];

        for (n=0;n<nVibModes;n++) res+=*(n+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_));
        return res;
      }

      inline void SetVibE(double e,int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
        *(nmode+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_))=e;
      }

      void InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart);
      void InitRotTemp(double RotTemp,PIC::ParticleBuffer::byte *ParticleDataStart);

      double GetCellRotTemp(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellVibTemp(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellVibTemp(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);

      double GetCellMeanRotE(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellMeanVibE(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);

      void RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell);

      //request data for the model
      int RequestSamplingData(int offset);

      //init
      void Init();

      //output the model data
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintVariableList(FILE* fout,int DataSetNumber);

      //calcualte the temperature index
      //get the temperature index
      inline double GetTempIndex(int s0,int s1) {
        static const double TemepratureIndex[1][1]={0.0};

        return TemepratureIndex[s0][s1];
      }
    }

    namespace qLB {
    using namespace LB;

//      extern int NmaxVibModes; //the maximum number of the vibrational modes for the set of used species
      extern int _VIBRATIONAL_GROUND_LEVEL_SAMPLE_DATA_OFFSET_,_VIBRATIONAL_FIRST_EXITED_LEVEL_SAMPLE_DATA_OFFSET_;

      int RequestSamplingData(int offset);
      void Init();

      double GetVibE(int nmode_in,PIC::ParticleBuffer::byte *ParticleDataStart);
      void InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart);
      void SetVibLevel(double VibQuantumNumber,int nmode,PIC::ParticleBuffer::byte *ParticleDataStart);
      int GetVibLevel(int nmode,PIC::ParticleBuffer::byte *ParticleDataStart);

      //the ground and the first exited vibrational states
      double GetGroundVibLevelPopulationFraction(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);
      double GetFirstExitedVibLevelPopulationFraction(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);

      double GetCellVibTemp(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellVibTemp(int nmode_in,int s,PIC::Mesh::cDataCenterNode* cell);

      void RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell);
    }

    inline double GetRotE(PIC::ParticleBuffer::byte *ParticleDataStart) {
      return LB::GetRotE(ParticleDataStart);
    }

    inline double GetVibE(int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
      return LB::GetVibE(nmode,ParticleDataStart);
    }

    inline void RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell) {
      LB::RedistributeEnergy(ptr0,ptr1,vrel,ChangeParticlePropertiesFlag,cell);
    }

    inline void InitRotTemp(double RotTemp,PIC::ParticleBuffer::byte *ParticleDataStart) {
      LB::InitRotTemp(RotTemp,ParticleDataStart);
    }

    inline void InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart) {
      LB::InitVibTemp(VibTemp,ParticleDataStart);
    }

    inline void Init() {
      LB::Init();
    }
  }


  namespace Mover {

//    #include "UserDefinition.PIC.Mover.h"

    //the return codes of the moving procedures
    #define _PARTICLE_REJECTED_ON_THE_FACE_ -1
    #define _PARTICLE_DELETED_ON_THE_FACE_   0
    #define _PARTICLE_CROSS_THE_FACE_        1
    #define _PARTICLE_LEFT_THE_DOMAIN_       2
    #define _PARTICLE_MOTION_FINISHED_       3

    namespace Sampling {
      namespace Errors {
         extern std::vector<int> RemovedModelParticles;
         extern std::vector<double> ModelParticleRemovingRate;
         extern std::vector<std::string> ErrorLineID;

         extern int ErrorDetectionFlag;

         void Init();
         void PrintData();
         void AddRemovedParticleData(double Rate, int spec, int line,const char *fname);
         void AddRemovedParticleData(double Rate, int nRemovedParticles, int spec, std::string &FullLineID);
         void AddRemovedParticleData(double Rate, int spec, std::string &ErrorID);
      }
    }

    namespace FieldLine {

      // procedure that returns parameters of the guiding center motion
      void GuidingCenterMotion(double *Vguide_perp,    double &ForceParal,
                               double &BAbsoluteValue, double *BDirection,
                               double *PParal,
			       int spec,long int ptr,double *x,double *v);
      //mover itself
      int Mover_SecondOrder(long int ptr, double Dt,
                            cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    }

    namespace GuidingCenter{

      namespace Sampling {
	      extern PIC::Mesh::cDatumWeighted DatumTotalKineticEnergy;
	      void SampleParticleData(char* ParticleData, double LocalParticleWeight, char* SamplingBuffer, int spec);
      }

      void Init_BeforeParser();
      void Init();
      void InitiateMagneticMoment(int spec, double *x, double *v, 
				  long int ptr, void *node);

      //mover
      void GuidingCenterMotion_default(double *Vguide, double &ForceParal, 
				       double &BAbsValue, double *BDirection, 
				       double *PParal,
				       int spec,long int ptr,
				       double *x, double *v,
				       cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
      int Mover_SecondOrder(long int ptr,double dtTotal,
			    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    }


    namespace Relativistic {

      int Boris(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    }


//    typedef void (*fTotalParticleAcceleration)(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
    typedef int (*fSpeciesDependentParticleMover) (long int,double,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    typedef int (*fSpeciesDependentParticleMover_BoundaryInjection) (long int,double,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,bool);


    //the vector containing the species specific particle moving procedures
//    extern fSpeciesDependentParticleMover *MoveParticleTimeStep;
//    extern fTotalParticleAcceleration TotalParticleAcceleration;
//    extern fSpeciesDependentParticleMover_BoundaryInjection *MoveParticleBoundaryInjection;

    //process a particle when it leaves the boundary of the computational domain
    typedef int (*fProcessOutsideDomainParticles) (long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
    extern fProcessOutsideDomainParticles ProcessOutsideDomainParticles;

    typedef int (*fProcessTriangleCutFaceIntersection) (long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    extern fProcessTriangleCutFaceIntersection ProcessTriangleCutFaceIntersection;

    void Init_BeforeParser();
    void Init();
    void MoveParticles();

    int UniformWeight_UniformTimeStep_noForce(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,bool FirstBoundaryFlag);

    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_SecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,bool FirstBoundaryFlag,CutCell::cTriangleFace *lastIntersectedTriangleFace=NULL);

    void TotalParticleAcceleration_default(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);

    int Boris(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    void BorisSplitAcceleration_default(double *accl, double *rotation, int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
  }


  namespace ParticleWeightTimeStep {
    typedef double (*fSetFunction) (int,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    extern fSetFunction LocalParticleWeight,LocalTimeStep,LocalBlockInjectionRate;

    //a user-defined function for calculation of an extra source rate
    typedef double (*fUserDefinedExtraSourceRate)(int);
    extern fUserDefinedExtraSourceRate UserDefinedExtraSourceRate;

    //the extra source rate calcualted by the exosphere model (src/models/exosphere)
    typedef double (*fExosphereModelExtraSourceRate)(int);
    extern fExosphereModelExtraSourceRate ExosphereModelExtraSourceRate;

    extern double maxReferenceInjectedParticleNumber;


    //when the global particle weight/time step are used, the following are the buffers where these parameters are stored
    extern double *GlobalParticleWeight,*GlobalTimeStep;

    double GetMaximumBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
    double GetTotalBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

    void initParticleWeight(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
    void SetGlobalParticleWeight(int spec,double GlobalParticleWeight,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

    double GetGlobalTimeStep(int spec);

    void initParticleWeight_ConstantWeight();
    void initParticleWeight_ConstantWeight(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
    void initParticleWeight_ConstantDensity(int spec,double NumberDensity,double TotalModelParticleNumber);


    void initTimeStep(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

    void copyLocalParticleWeightDistribution(int specTarget,int specSource,double ProportionaltyCoefficient=1.0);
    void copyLocalTimeStepDistribution(int specTarger,int specSource,double ProportionaltyCoefficient=1.0);

    //adjust particle weight so Weight/dT=const in all blocks (need to be called after the time step and particle weight are initialized
    void AdjustParticleWeight_ConstantWeightOverTimeStep(int spec,double WeightOverTimeStepRatio,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
    void AdjustParticleWeight_ConstantWeightOverTimeStep_KeepMinParticleWeight(int spec);

    double GetMinLocalParticleWeightValue(int spec,double &WeightOverTimeStepRatio,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
  }

  //in case of time-dependent model runs with glocal time stop for all specis -> count the physical time of the simulation
  namespace SimulationTime {
    extern double TimeCounter;
    void SetInitialValue(double InitalTimeCounterValue);
    double Get();
    void Update();
  }

  namespace Parallel {


     //count the number of particles that were send and recieve by the thread
     extern long int sendParticleCounter,recvParticleCounter,IterationNumberAfterRebalancing;
     extern double RebalancingTime,CumulativeLatency;

     //the factor the trrigeres the emergency load rebalancing. The condition for the rebalancing:
     //(PIC::Parallel::CumulativeLatency>PIC::Parallel::EmergencyLoadRebalancingFactor*PIC::Parallel::RebalancingTime)
     extern double EmergencyLoadRebalancingFactor;

     //exchenge paricles between iterations
     void ExchangeParticleData();

     //Latency of the run
     extern double Latency;
  }

  namespace Debugger {
    //contains functions that are used for debugging the code

    //catch variation of a variable located at a particular address
    //MatchMode == true -> print when the match is found; MatchMode == false -> print when the variables are not match
    template <typename  T>
    inline void Catch(T val,T* ptr,int nline,const char *fname,bool MatchMode) {
      if (ptr!=NULL) {
        if (MatchMode==true) {
          if (*ptr==val) {
            printf("PIC::Debugger::Catch: Match is found (line=%i,file=%s)\n",nline,fname);
          }
        }

        if (MatchMode==false) {
          if (*ptr!=val) {
            cout << "PIC::Debugger::Catch: NO Match is found. [Val=" << *ptr << ", ref=" << val << " (line=" << nline << ", file=" << fname << ")\n";
          }
        }
      }
    }

    //InfiniteLoop==false ==> no problem found; InfiniteLoop==true => the actual number of particles does not consider the that in teh particle buffer
    bool InfiniteLoop(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=NULL);
    void FindDoubleReferencedParticle(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=NULL);

    //check is a variable value is within an allowed range
    const double minAllowedValue=1.0E-80;
    const double maxAllowedValue=1.0E80;

    inline void CatchOutLimitValue(double val,int line,const char *fname) {
      double t=fabs(val);

      if ( (isfinite(val)==false)|| ((t>0.0)&&((t<minAllowedValue)||(t>maxAllowedValue))) ) {
        char msg[600];

        sprintf(msg,"Error: value of limits: val=%e, (line=%i,file=%s)",val,line,fname);
        exit(line,fname,msg);
      }
    }

    inline void CatchOutLimitValue(double* valArray,int lengthArray, int line,const char *fname) {
      int i;

      for (i=0;i<lengthArray;i++) CatchOutLimitValue(valArray[i],line,fname);
    }
  }

  namespace Alarm {
    //finish the execution of the code at a particular value of the walltime
    extern bool AlarmInitialized,WallTimeExeedsLimit;
    extern double StartTime;
    extern double RequestedExecutionWallTime;

    inline void SetAlarm(double requestedWallTime) {
      AlarmInitialized=true;
      StartTime=MPI_Wtime();
      RequestedExecutionWallTime=requestedWallTime;
    }

    inline void FinishExecution() {
      MPI_Finalize();
      printf("$PREFIX:!!!!! Execution is finished by the alarm (PIC::Alarm) !!!!!!\n");
      exit(__LINE__,__FILE__,"!!!!! Execution is finished by the alarm !!!!!!");
    }
  }

  namespace ColumnIntegration {

    //define 3 nodes on the surface of a bounding plane; index value: 0 -> xmin component of the coordinate, 1 -> xmax component of the coordinate
    static const int x0PlaneNodeIndex[6][3]={ {0,0,0},{1,0,0},       {0,0,0},{0,1,0},           {0,0,0},{0,0,1}};
    static const int x1PlaneNodeIndex[6][3]={ {0,1,0},{1,1,0},       {1,0,0},{1,1,0},           {1,0,0},{1,0,1}};
    static const int x2PlaneNodeIndex[6][3]={ {0,0,1},{1,0,1},       {0,0,1},{0,1,1},           {0,1,0},{0,1,1}};
    static const int PlaneNormal[6][3]=     { {1,0,0},{1,0,0},       {0,1,0},{0,1,0},           {0,0,1},{0,0,1}};

    struct cBoundingBoxFace {
      double x[3];
      double e0[3];
      double e1[3];
      double Normal[3];
      double e0Length;
      double e1Length;
    };

    extern cBoundingBoxFace BoundingBoxFace[6];
    extern bool InitializedFlag;

    //control initialization of the model
    extern bool ModelInitFlag;

    void Init();
    bool FindIntegrationLimits(double *x0,double *l,double& IntegrationPathLength,double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &xStartNode,double *xFinish,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &xFinishNode);

//    //get a single integrated value
//    double GetCoulumnIntegral(double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xStartNode,double *l,double IntegrationPathLength,double (*Integrand)(double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));
//    double GetCoulumnIntegral(double *x0,double *l,double (*Integrand)(double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));

    //get values for multiple integrals
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xStartNode,double *l,double IntegrationPathLength,void (*Integrand)(double*,int,double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *x0,double *l,void (*Integrand)(double*,int,double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));
  }


  //interpolation routines
  namespace InterpolationRoutines {

    namespace CellCentered {
      //the maximum number of the elements in the interpolation stencil
      const int nMaxStencilLength=64;

      class cStencil {
      public:
        int Length;
        double Weight[nMaxStencilLength];
        PIC::Mesh::cDataCenterNode* cell[nMaxStencilLength];

        void flush() {
          Length=0;
          for (int i=0;i<nMaxStencilLength;i++) Weight[i]=0.0,cell[i]=NULL;
        }

        void print() {
          printf("Length=%i, weight:",Length);
          for (int i=0;i<Length;i++) printf (" %e[%i]",Weight[i],i);
          printf("\n");  
        }

        void AddCell(double w,PIC::Mesh::cDataCenterNode* c) {
          if (Length==nMaxStencilLength) exit(__LINE__,__FILE__,"Error: the stencil length exeeds 'nMaxStencilLength'. Need to increase 'nMaxStencilLength'");

          Weight[Length]=w;
          cell[Length]=c;
          Length++;
        }

        cStencil() {flush();}
      };

      extern cStencil* StencilTable;

      //types of the cell ceneterd interpolating rourines implemented in AMPS
      namespace Constant {
        PIC::InterpolationRoutines::CellCentered::cStencil *InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL);
      }

      namespace Linear {
        //interface to the fortran written interpolation library
        namespace INTERFACE {
          //list of pointers to nodes to identify them as integers
          //less than 2*2^3 integers might be needed to get an interpolation stencils
          const  int nBlockFoundMax=16;
          extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* BlockFound[nBlockFoundMax];
          //index of the current position in the list
          extern int iBlockFoundCurrent;

          //the last node where the stencil point has been found by the 'find' function
          extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*  last;
        }

        // precision for taking shortcut and moving point to the cell center
        const double PrecisionCellCenter = 1.0e-3;

        //interpolation functions
        PIC::InterpolationRoutines::CellCentered::cStencil *InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL);
        PIC::InterpolationRoutines::CellCentered::cStencil *GetTriliniarInterpolationStencil(double iLoc,double jLoc,double kLoc,double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
        PIC::InterpolationRoutines::CellCentered::cStencil *GetTriliniarInterpolationMutiBlockStencil(double *x,double *xStencilMin,double *xStencilMax,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
      }

    }

    void Init();
  }

  //namespace CPLR contains definitions of all couplers used in AMPS
  namespace CPLR {

    //Set the interpolation stencil that is used for the interpolation in the coupler
   void InitInterpolationStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL);

    //coupling with SWMF
    namespace SWMF {
      extern int MagneticFieldOffset,TotalDataLength,BulkVelocityOffset,PlasmaPressureOffset;
      extern int PlasmaNumberDensityOffset,PlasmaTemperatureOffset;

      //the mean mass of the plasma speies atoms/molecules (needed to conver mass density into number density)
      extern double MeanPlasmaAtomicMass;

      //the flug if 'false; by default and is teruned to 'true' after the first coupling procedure (used to pospond initialization of AMPS till the backround field information is exported to AMPS)
      extern bool FirstCouplingOccured;

      //init the coupler
      void init();
      void ConvertMpiCommunicatorFortran2C(signed int* iComm,signed int* iProc,signed int* nProc);

      //output the interpolated data into a file
      int RequestDataBuffer(int offset);
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);

      //prepare the list of the coordinates for the interpolation
      void ResetCenterPointProcessingFlag();
      void GetCenterPointNumber(int *nCenterPoints);
      void GetCenterPointCoordinates(double *x);
      void RecieveCenterPointData(char* ValiableList, int nVarialbes,double *data,int *index);

      //send AMPS's data to SWMF: multiple user-defined send routines could be implemeted
/*
      Description of the parameters: ARRAYS are fortran style - second index is first!!!!!!!!!
      character(len=*), intent(in):: NameVar ! List of variables
      integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
      integer,          intent(in):: nDimIn  ! Dimensionality of positions
      integer,          intent(in):: nPoint  ! Number of points in Xyz_DI
      real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
      real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array*/

      typedef void (*fSendCenterPointData)(char *NameVar, int *nVarIn, int *nDimIn, int *nPoint, double *Xyz_DI, double *Data_VI);
      extern list <fSendCenterPointData> SendCenterPointData;


      //calcualte physical parameters
      inline void GetBackgroundMagneticField(double *B,PIC::Mesh::cDataCenterNode *cell) {
        int idim;
        double *offset=(double*)(MagneticFieldOffset+cell->GetAssociatedDataBufferPointer());

        for (idim=0;idim<3;idim++) B[idim]=offset[idim];
      }

      inline void GetBackgroundPlasmaVelocity(double *v,PIC::Mesh::cDataCenterNode *cell) {
        int idim;
        double *offset=(double*)(BulkVelocityOffset+cell->GetAssociatedDataBufferPointer());

        for (idim=0;idim<3;idim++) v[idim]=offset[idim];
      }



      inline double GetBackgroundPlasmaPressure(PIC::Mesh::cDataCenterNode *cell) {
        return *((double*)(PlasmaPressureOffset+cell->GetAssociatedDataBufferPointer()));
      }

      inline double GetBackgroundPlasmaNumberDensity(PIC::Mesh::cDataCenterNode *cell) {
        return *((double*)(PlasmaNumberDensityOffset+cell->GetAssociatedDataBufferPointer()));
      }

      inline double GetBackgroundPlasmaTemperature(PIC::Mesh::cDataCenterNode *cell) {
        return *((double*)(PlasmaTemperatureOffset+cell->GetAssociatedDataBufferPointer()));
      }

      inline void GetBackgroundElectricField(double *E,PIC::Mesh::cDataCenterNode *cell) {
        double B[3],v[3];

        GetBackgroundMagneticField(B,cell);
        GetBackgroundPlasmaVelocity(v,cell);

        E[0]=-(v[1]*B[2]-B[1]*v[2]);
        E[1]=-(-v[0]*B[2]+B[0]*v[2]);
        E[2]=-(v[0]*B[1]-B[0]*v[1]);
      }

      inline void GetBackgroundFieldsVector(double *E,double *B,PIC::Mesh::cDataCenterNode *cell) {
        GetBackgroundMagneticField(B,cell);
        GetBackgroundElectricField(E,cell);
      }

    }


    //coupling thrugh a file
    namespace DATAFILE {

      namespace MULTIFILE {
	      //-----------------------------------
	      //number of file to be loaded
	      extern int nFile;

        // number of the next data file in the schedule
        extern int iFileLoadNext;

	      //schedule for loading multiple data files
	      struct cScheduleItem {
	        double Time;
	        char FileName[_MAX_STRING_LENGTH_PIC_];
	      };

	      extern vector<cScheduleItem> Schedule;

	      //comparison function for schedule items: for sorting
	      inline bool _compare(cScheduleItem Item1, cScheduleItem Item2) {
	        return Item1.Time < Item2.Time;
	      }

	      // name of file with table defining schedule
	      extern char FileTable[_MAX_STRING_LENGTH_PIC_];

	      //offsets to the data at next/current datafile
	      //used for time interpolation
	      extern int NextDataFileOffset;
	      extern int CurrDataFileOffset;



	      //variable to track whether to break simulation at the last datafile
	      extern bool BreakAtLastFile;

	      //variable to track whether the last datafile has been reached
	      extern bool ReachedLastFile;

        //check whether it is time to load the next file
        inline bool IsTimeToUpdate() {
          #if     _PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_
          return PIC::SimulationTime::Get() >= Schedule[iFileLoadNext-1].Time;
          #else
          return PIC::SimulationTime::Get() >= Schedule[iFileLoadNext].  Time;
          #endif
        }

        //initialize
        void Init(bool BreakAtLastFileIn = true, int  FileNumber        = 0);

	      //load the schedule
	      void GetSchedule();

	      //extract time from datafile itself
	      double GetFileTime(const char* FileName);

	      //update the datafile
        void UpdateDataFile();
      }

      //path to the location of the datafiles
      extern char path[_MAX_STRING_LENGTH_PIC_];



      //the offset from the cell->AssociatedData()
      extern int CenterNodeAssociatedDataOffsetBegin;
      extern int nTotalBackgroundVariables;

      //Physical quantaties offsets that could be read and stored
      struct cOffsetElement {
        bool active;
        bool allocate;

        int nVars;
        char VarList[_MAX_STRING_LENGTH_PIC_];
        int offset;
      };

      namespace Offset {
        extern bool InitFlag;

        extern cOffsetElement PlasmaNumberDensity;
        extern cOffsetElement PlasmaBulkVelocity;
        extern cOffsetElement PlasmaTemperature;
        extern cOffsetElement PlasmaIonPressure;
        extern cOffsetElement PlasmaElectronPressure;
        extern cOffsetElement MagneticField;
        extern cOffsetElement ElectricField;
        extern cOffsetElement MagneticFieldGradient;
	      extern cOffsetElement MagneticFluxFunction;


        inline void SetAllocate(bool flag,cOffsetElement* offset) {
          if (InitFlag==false) {
            if (flag==false) exit(__LINE__,__FILE__,"Error: the background data fields can be only added not removed");

            offset->allocate=flag;
          }
          else exit(__LINE__,__FILE__,"Error: changing of the allocation flags is not allowed after initialization of the offsets is completed");
        }

        inline void SetActive(bool flag,cOffsetElement* offset) {
          if ((flag==true)&&(offset->allocate==false)) exit(__LINE__,__FILE__,"Error: the offset cannot be set 'active' if it is not allocated before during the initialization");
          offset->active=flag;
        }

        inline void SetAllocateAll(bool flag) {
          SetAllocate(flag,&PlasmaNumberDensity);
          SetAllocate(flag,&PlasmaBulkVelocity);
          SetAllocate(flag,&PlasmaTemperature);
          SetAllocate(flag,&PlasmaIonPressure);
          SetAllocate(flag,&PlasmaElectronPressure);
          SetAllocate(flag,&MagneticField);
          SetAllocate(flag,&ElectricField);
          SetAllocate(flag,&MagneticFieldGradient);
        }

        inline void SetActiveAll(bool flag) {
          SetActive(flag,&PlasmaNumberDensity);
          SetActive(flag,&PlasmaBulkVelocity);
          SetActive(flag,&PlasmaTemperature);
          SetActive(flag,&PlasmaIonPressure);
          SetActive(flag,&PlasmaElectronPressure);
          SetActive(flag,&MagneticField);
          SetActive(flag,&ElectricField);
          SetActive(flag,&MagneticFieldGradient);
        }
      }

      //load the datafile
      void ImportData(const char *fname);

      //routines to generate additional data
      void GenerateMagneticFieldGradient(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

      //return the interpolated value of the background data
      inline void GetBackgroundData(double *DataVector,int DataVectorLength,int BackgroundDataOffset,PIC::Mesh::cDataCenterNode *CenterNode) {
        int i;
        double *offset=(double*)(BackgroundDataOffset+CenterNode->GetAssociatedDataBufferPointer());

        for (i=0;i<DataVectorLength;i++) DataVector[i]=offset[i];
      }

      inline void GetBackgroundData(double *DataVector,int DataVectorLength,int BackgroundDataOffset,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
        int i;
        double *offset=(double*)(BackgroundDataOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

        for (i=0;i<DataVectorLength;i++) DataVector[i]=offset[i];
      }

      //save/read the background data binary file
      bool BinaryFileExists(const char *fNameBase);
      void SaveBinaryFile(const char *fNameBase,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
      void LoadBinaryFile(const char *fNameBase,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);


      //print the background variables into AMPS' output file
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);

      //test the data file reader by comparing the reading results with the reference data
      void SaveTestReferenceData(const char* fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

      //initialize the data file reader
      void Init();

      //print the ion flux at a sphere
      void PrintSphereSurfaceIonFlux(char const* fname,double SphereRadius);
      void EvaluateSurfaceIonFlux(double ShiftFactor);

      //calculate the values of the located parameters
      inline void GetBackgroundValue(double *DataVector,int DataVectorLength,int DataOffsetBegin,PIC::Mesh::cDataCenterNode *cell, double Time) {
        double *offset = (double*)(DataOffsetBegin + MULTIFILE::CurrDataFileOffset + cell->GetAssociatedDataBufferPointer());

        for (int i=0;i<DataVectorLength;i++) DataVector[i]=offset[i];

#if  _PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_
        if (isnan(Time)) Time = PIC::SimulationTime::Get();

        offset = (double*)(DataOffsetBegin+MULTIFILE::NextDataFileOffset+cell->GetAssociatedDataBufferPointer());

         //interpolation weight
        double alpha = (MULTIFILE::Schedule[MULTIFILE::iFileLoadNext-1].Time - Time) /
        (MULTIFILE::Schedule[MULTIFILE::iFileLoadNext-1].Time - MULTIFILE::Schedule[MULTIFILE::iFileLoadNext-2].Time);

        for (int i=0;i<DataVectorLength;i++)  DataVector[i] = DataVector[i] * alpha + offset[i] * (1-alpha);
#endif//_PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_

      #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
        PIC::Debugger::CatchOutLimitValue(DataVector,DataVectorLength,__LINE__,__FILE__);
      #endif
      #endif

      }

      inline void GetBackgroundElectricField(double *E,PIC::Mesh::cDataCenterNode *cell, double Time) {
      	GetBackgroundValue(E, Offset::ElectricField.nVars,Offset::ElectricField.offset, cell, Time);
      }

      inline void GetBackgroundMagneticField(double *B,PIC::Mesh::cDataCenterNode *cell, double Time) {
      	GetBackgroundValue(B,Offset::MagneticField.nVars,Offset::MagneticField.offset, cell, Time);
      }

      inline void GetBackgroundMagneticFieldGradient(double *gradB,PIC::Mesh::cDataCenterNode *cell, double Time) {
	      GetBackgroundValue(gradB,Offset::MagneticFieldGradient.nVars,Offset::MagneticFieldGradient.offset, cell, Time);
      }

      inline void GetBackgroundPlasmaVelocity(double *vel,PIC::Mesh::cDataCenterNode *cell, double Time) {
	      GetBackgroundValue(vel,Offset::PlasmaBulkVelocity.nVars,Offset::PlasmaBulkVelocity.offset, cell, Time);
      }

      inline double GetBackgroundMagneticFluxFunction(PIC::Mesh::cDataCenterNode *cell, double Time) {
	      double ff;

	      GetBackgroundValue(&ff,Offset::MagneticFluxFunction.nVars,Offset::MagneticFluxFunction.offset, cell, Time);
	      return ff;
      }

      inline double GetBackgroundPlasmaPressure(PIC::Mesh::cDataCenterNode *cell, double Time) {
	      double p;

	      GetBackgroundValue(&p,Offset::PlasmaIonPressure.nVars,Offset::PlasmaIonPressure.offset, cell, Time);
	      return p;
      }

      inline double GetBackgroundElectronPlasmaPressure(PIC::Mesh::cDataCenterNode *cell, double Time) {
        double p;
        int nVars, offset;

        if(Offset::PlasmaElectronPressure.offset!=-1) {
          nVars  = Offset::PlasmaElectronPressure.nVars;
          offset = Offset::PlasmaElectronPressure.offset;
        }
        else {
          nVars  = Offset::PlasmaIonPressure.nVars;
          offset = Offset::PlasmaIonPressure.offset;
        }

        GetBackgroundValue(&p,nVars,offset, cell, Time);
        return p;
      }

      inline double GetBackgroundPlasmaNumberDensity(PIC::Mesh::cDataCenterNode *cell, double Time) {
	      double n;

	       GetBackgroundValue(&n,Offset::PlasmaNumberDensity.nVars,Offset::PlasmaNumberDensity.offset, cell, Time);
         return n;
      }	

      inline double GetBackgroundPlasmaTemperature(PIC::Mesh::cDataCenterNode *cell, double Time) {
	       double T;

	       GetBackgroundValue(&T, Offset::PlasmaTemperature.nVars,Offset::PlasmaTemperature.offset, cell, Time);
         return T;
      }

      inline void GetBackgroundFieldsVector(double *E,double *B,PIC::Mesh::cDataCenterNode *cell, double Time) {
         GetBackgroundValue(E, Offset::ElectricField.nVars,Offset::ElectricField.offset, cell, Time);
         GetBackgroundValue(B, Offset::MagneticField.nVars,Offset::MagneticField.offset, cell, Time);
      }

      //individual reading procedures
      namespace ICES {
        extern char locationICES[_MAX_STRING_LENGTH_PIC_]; //location of the data and the dace cases
        extern char ModeCaseSWMF[_MAX_STRING_LENGTH_PIC_]; //the name of the model case that will be used for interpolation with ICES

        void Init();
        void SetLocationICES(const char*);

        //the total number of bytes used to store the ICES data vector; the offset of the data associated with the ICES data vector
        extern int TotalAssociatedDataLength,AssociatedDataOffset;

        //the offsets for the plasma parameters loaded with ICES are stores in PIC::CPLR::DATAFILE::Offset
        //the offsets for parameters loaded from the DSMC model
        extern int NeutralBullVelocityOffset,NeutralNumberDensityOffset,NeutralTemperatureOffset,DataStatusOffsetDSMC;

        //calcualte the total number of cells in the mesh
        long int getTotalCellNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

        //create the trajectory file
        void createCellCenterCoordinateList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

        //evaluate the ion flux at the surface of the spehrical interval body
        void EvaluateSurfaceIonFlux(double ShiftFactor=1.0);


        //retrive the SWMF data file
        #define _PIC_ICES__STATUS_OK_   0

        class cDataNodeSWMF {
        public:
          double swNumberDensity,swTemperature,swPressure,E[3],B[3],swVel[3];
          int status;

          void flush() {
            swNumberDensity=0.0,swTemperature=0.0,swPressure=0.0;
            for (int i=0;i<3;i++) E[i]=0.0,B[i]=0.0,swVel[i]=0.0;
          }
        };

        class cDataNodeDSMC {
        public:
          double neutralNumberDensity,neutralTemperature,neutralVel[3];
          int status;

          void flush() {
            neutralNumberDensity=0.0,neutralTemperature=0.0;

            for (int i=0;i<3;i++) neutralVel[i]=0.0;
          }
        };

        void retriveSWMFdata(const char *DataFile=ModeCaseSWMF);
        void retriveDSMCdata(const char *Case,const char *DataFile,const char *MeshFile);


        void readSWMFdata(const double MeanIonMass,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree); //MeanIonMass -> the mean ion mass of the plasma flow in [amu]
        void readDSMCdata(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

        //user defined pre-processor of the data that is readed by ICES
        typedef void (*fDSMCdataPreProcessor)(double *x,cDataNodeDSMC& data);
        typedef void (*fSWMFdataPreProcessor)(double *x,cDataNodeSWMF& data);

        extern fDSMCdataPreProcessor DSMCdataPreProcessor;
        extern fSWMFdataPreProcessor SWMFdataPreProcessor;

        void PrintVariableList(FILE* fout,int DataSetNumber);
        void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
        void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

        //print the ion flux at a sphere
        void PrintSphereSurfaceIonFlux(char const* fname,double SphereRadius);
      }

      namespace KAMELEON {
        //the mass of the dominant background plasma ion
        extern double PlasmaSpeciesAtomicMass;

        //the conversion factor from units used in the external model (MATSRUS,LFM, etc....) into meters: [m]=cdfDataFile2m_ConversionFactor*[external model units]
        extern double cdfDataFile2m_ConversionFactor;

        //init the reader
        void Init();

	      void LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
	      void GetDomainLimits(double *xmin,double *xmax,const char *fname);

        namespace LFM {
          void LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
          void GetDomainLimits(double *xmin,double *xmax,const char *fname);
        }

      }

      namespace ARMS {
        void Init();
	      double GetFileTime(const char *fname);
        void LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
      }


      namespace BATSRUS {
        extern double PlasmaSpeciesAtomicMass; //the mass of the dominant background plasma ion
        extern double UnitLength; //the spatial units used in the BATSRUS' output file
        extern char filename[_MAX_STRING_LENGTH_PIC_];
        extern bool InitFlag;

        //reserve memory for reading the backdround data
        void Init();

        //open the file and init BATL libraty
        void Init(const char *fname);

        void GetDomainLimits(double *xmin,double *xmax);
        void LoadDataFile(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

	      //the offsets of the physical variables in the .idl file
        extern int rhoBATSRUS2AMPS;
        extern int mxBATSRUS2AMPS,myBATSRUS2AMPS,mzBATSRUS2AMPS;
        extern int uxBATSRUS2AMPS,uyBATSRUS2AMPS,uzBATSRUS2AMPS;
        extern int bxBATSRUS2AMPS,byBATSRUS2AMPS,bzBATSRUS2AMPS;
        extern int pBATSRUS2AMPS;
      }

      namespace TECPLOT {
        extern double xDataMin[3],xDataMax[3]; //the size of the domain. used when DataMode==DataMode_XYZ
        extern double rDataMin,rDataMax; //the radial size of the domain. used when DataMode==DataMode_SPERICAL

        extern double UnitLength; //the spatial units used in the BATSRUS' output file
        extern int maxScriptPointNumber; //the maximum number of point that one single script can have
        extern int nTotalVarlablesTECPLOT; //the total number of the variabled in the TECPLOT output file (including thoses that are not needed)

        const int DataMode_SPHERICAL=0;
        const int DataMode_XYZ=1;
        extern int DataMode;

        //rotation matrix
        extern double RotationMatrix_LocalFrame2DATAFILE[3][3];
        extern double RotationMatrix_DATAFILE2LocalFrame[3][3];
        void SetRotationMatrix_LocalFrame2DATAFILE(const double Matrix[3][3]);
        void SetRotationMatrix_DATAFILE2LocalFrame(const double Matrix[3][3]);
        void CalculateInverseRotationMatrix(double Matrix[3][3],double InverseMatrix[3][3]);

        class cLoadedVariableData {
        public:
          int offset;
          double ScaleFactor;

          cLoadedVariableData() {
            offset=-1,ScaleFactor=0.0;
          }
        };

        extern cLoadedVariableData Velocity,IonPressure,ElectronPressure,MagneticField,Density;
        inline void SetLoadedDensityVariableData(int offset,double ScaleFactor) {Density.offset=offset-1,Density.ScaleFactor=ScaleFactor;}
        inline void SetLoadedVelocityVariableData(int offset,double ScaleFactor) {Velocity.offset=offset-1,Velocity.ScaleFactor=ScaleFactor;}
        inline void SetLoadedIonPressureVariableData(int offset,double ScaleFactor) {IonPressure.offset=offset-1,IonPressure.ScaleFactor=ScaleFactor;}
        inline void SetLoadedElectronPressureVariableData(int offset,double ScaleFactor) {ElectronPressure.offset=offset-1,ElectronPressure.ScaleFactor=ScaleFactor;}
        inline void SetLoadedMagneticFieldVariableData(int offset,double ScaleFactor) {MagneticField.offset=offset-1,MagneticField.ScaleFactor=ScaleFactor;}

        void Init();

        void SetDomainLimitsXYZ(double *xmin,double *xmax);
        void SetDomainLimitsSPHERICAL(double rmin,double rmax);
        void ExtractData(const char *fname);

        void ResetCellProcessingFlag(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

        //function CreatePointList: 1. calculates the number of the points the will be interpolated and 2. is fScript!=NULL save tham into fScript
        int CountInterpolatedPointNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

        int CreateScript(const char *ScriptBaseName,const char* DataFileTECPLOT,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
        void LoadDataFile(const char *fname,int nTotalOutputFiles,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

        //the function call all nessesary methods of the TECPLOT namespace to export the data
        void ImportData(const char* fname);
      }
    }

    //save and load the center node associated data from the AMPS' data buffers
    void SaveCenterNodeAssociatedData(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
    void LoadCenterNodeAssociatedData(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

    inline void GetBackgroundElectricField(double *E, double Time = NAN) {
      double t[3];
      int idim,iStencil;
      PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

      for (idim=0;idim<3;idim++) E[idim]=0.0;

      #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
      #else
      memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
      #endif

      for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
        #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
        SWMF::GetBackgroundElectricField(t,Stencil.cell[iStencil]);
        #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__T96_
        for (int i=0;i<3;i++) t[i]=0.0;
        #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
        DATAFILE::GetBackgroundElectricField(t,Stencil.cell[iStencil], Time);
        #else
        exit(__LINE__,__FILE__,"not implemented");
        #endif

        for (idim=0;idim<3;idim++) E[idim]+=Stencil.Weight[iStencil]*t[idim];
      }
    }

     inline void GetBackgroundMagneticField(double *B, double Time = NAN) {
       double t[3];
       int idim,iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       for (idim=0;idim<3;idim++) B[idim]=0.0;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         SWMF::GetBackgroundMagneticField(t,Stencil.cell[iStencil]);
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__T96_
         DATAFILE::GetBackgroundData(t,3,DATAFILE::Offset::MagneticField.offset,Stencil.cell[iStencil]);
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         DATAFILE::GetBackgroundMagneticField(t,Stencil.cell[iStencil], Time);
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif

         for (idim=0;idim<3;idim++) B[idim]+=Stencil.Weight[iStencil]*t[idim];
       }
     }

     inline double GetAbsBackgroundMagneticField(double Time = NAN) {
       double B[3];

       GetBackgroundMagneticField(B,Time);
       return sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
     }

     inline void GetBackgroundMagneticFieldGradient(double *gradB, double Time = NAN) {
       double t[DATAFILE::Offset::MagneticFieldGradient.nVars];
       int idim,iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       // structure of gradB is the following
       //   gradB[0:2] = {d/dx, d/dy, d/dz} B_x
       //   gradB[3:5] = {d/dx, d/dy, d/dz} B_y
       //   gradB[6:8] = {d/dx, d/dy, d/dz} B_z

       for (idim=0;idim<DATAFILE::Offset::MagneticFieldGradient.nVars;idim++) gradB[idim]=0.0;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         exit(__LINE__,__FILE__,"not implemented");
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         DATAFILE::GetBackgroundMagneticFieldGradient(t,Stencil.cell[iStencil], Time);
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif

         for (idim=0;idim<DATAFILE::Offset::MagneticFieldGradient.nVars;idim++) gradB[idim]+=Stencil.Weight[iStencil]*t[idim];
       }
     }

     inline void GetAbsBackgroundMagneticFieldGradient(double *gradAbsB, double Time = NAN) {
       double gradB[9],B[3],AbsB,b[3];

       // structure of gradB is the following
       //   gradB[0:2] = {d/dx, d/dy, d/dz} B_x
       //   gradB[3:5] = {d/dx, d/dy, d/dz} B_y
       //   gradB[6:8] = {d/dx, d/dy, d/dz} B_z
       GetBackgroundMagneticFieldGradient(gradB,Time);
       PIC::CPLR::GetBackgroundMagneticField(B);
       AbsB=pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2],0.5);

       if (fabs(AbsB)>0.0) {
         b[0] = B[0]/AbsB; b[1] = B[1]/AbsB; b[2] = B[2]/AbsB;

         gradAbsB[0]= b[0] * gradB[0] + b[1] * gradB[3] + b[2] * gradB[6];
         gradAbsB[1]= b[0] * gradB[1] + b[1] * gradB[4] + b[2] * gradB[7];
         gradAbsB[2]= b[0] * gradB[2] + b[1] * gradB[5] + b[2] * gradB[8];
       }
       else for (int i=0;i<3;i++) gradAbsB[i]=0.0;
     }

     inline void GetCurlBackgroundMagneticField(double *CurlB, double Time = NAN) {
       double gradB[9];

       // structure of gradB is the following
       //   gradB[0:2] = {d/dx, d/dy, d/dz} B_x
       //   gradB[3:5] = {d/dx, d/dy, d/dz} B_y
       //   gradB[6:8] = {d/dx, d/dy, d/dz} B_z
       GetBackgroundMagneticFieldGradient(gradB,Time);

       CurlB[0]=+(gradB[7]-gradB[5]);
       CurlB[1]=-(gradB[6]-gradB[2]);
       CurlB[2]=+(gradB[3]-gradB[1]);
     }

     //calculate a particle drift velocity (Elkington-2002-JASTP)
      void GetDriftVelocity(double *vDrift,double *ParticleVelocity,double ParticleMass,double ParticleCharge,double Time = NAN);

      /*{
       double E[3],B[3],absB2,absB,absB4,t[3],c;
       double M,gamma,gradAbsB_perp[3],ParticleMomentum[3],ParticleMomentum_normB[3],pParallel;
       int idim;

       //get the particle momentum
       gamma=1.0/sqrt(1.0-sqrt(ParticleVelocity[0]*ParticleVelocity[0]+ParticleVelocity[1]*ParticleVelocity[1]+ParticleVelocity[2]*ParticleVelocity[2])/SpeedOfLight);
       for (idim=0;idim<3;idim++) ParticleMomentum[idim]=gamma*ParticleMass*ParticleVelocity[idim],vDrift[idim]=0.0;

       //get the background fields
       PIC::CPLR::GetBackgroundMagneticField(B);
       PIC::CPLR::GetBackgroundElectricField(E);
       absB2=B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
       absB=sqrt(absB2);
       absB4=absB2*absB2;

       //E cross B drift (avarage the drift velocities directly)
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;
       int iStencil;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         double cellB[3],cellE[3];

         //loop through all cells of the stencil
         switch (_PIC_COUPLER_MODE_) {
         case _PIC_COUPLER_MODE__SWMF_:
           SWMF::GetBackgroundElectricField(cellE,Stencil.cell[iStencil]);
           SWMF::GetBackgroundMagneticField(cellB,Stencil.cell[iStencil]);

           break;
         case _PIC_COUPLER_MODE__T96_:
           for (int i=0;i<3;i++) cellE[i]=0.0;
           DATAFILE::GetBackgroundData(cellB,3,DATAFILE::Offset::MagneticField.offset,Stencil.cell[iStencil]);

           break;
         case _PIC_COUPLER_MODE__DATAFILE_:
           DATAFILE::GetBackgroundElectricField(cellE,Stencil.cell[iStencil],Time);
           DATAFILE::GetBackgroundMagneticField(cellB,Stencil.cell[iStencil],Time);

           break;
         default:
           exit(__LINE__,__FILE__,"not implemented");
         }

         Vector3D::CrossProduct(t,cellE,cellB);
         c=Stencil.Weight[iStencil]*SpeedOfLight/(cellB[0]*cellB[0]+cellB[1]*cellB[1]+cellB[2]*cellB[2]);

         for (idim=0;idim<3;idim++) vDrift[idim]+=c*t[idim];
       }

       //next
       memcpy(ParticleMomentum_normB,ParticleMomentum,3*sizeof(double));
       Vector3D::Orthogonalize(B,ParticleMomentum_normB);
       M=pow(Vector3D::Length(ParticleMomentum_normB),2)/(2.0*ParticleMass*absB);

       GetAbsBackgroundMagneticFieldGradient(gradAbsB_perp,Time);
       Vector3D::Orthogonalize(B,gradAbsB_perp);
       Vector3D::CrossProduct(t,B,gradAbsB_perp);

       c=M*SpeedOfLight/(ParticleCharge*gamma*absB2);
       for (idim=0;idim<3;idim++) vDrift[idim]+=c*t[idim];

       //next drift coeffecient
       double gradB[9],t1[3];

       // structure of gradB is the following
       //   gradB[0:2] = {d/dx, d/dy, d/dz} B_x
       //   gradB[3:5] = {d/dx, d/dy, d/dz} B_y
       //   gradB[6:8] = {d/dx, d/dy, d/dz} B_z
       GetBackgroundMagneticFieldGradient(gradB,Time);

       //t1=(b\dot)b
       t1[0]=(B[0]*gradB[0]+B[1]*gradB[1]+B[2]*gradB[2])/absB4;
       t1[1]=(B[0]*gradB[3]+B[1]*gradB[4]+B[2]*gradB[5])/absB4;
       t1[2]=(B[0]*gradB[6]+B[1]*gradB[7]+B[2]*gradB[8])/absB4;

       Vector3D::CrossProduct(t,B,t1);
       pParallel=Vector3D::ParallelComponentLength(ParticleMomentum,B);

       c=SpeedOfLight*pow(pParallel,2)/(ParticleCharge*gamma*ParticleMass);
       for (idim=0;idim<3;idim++) vDrift[idim]+=c*t[idim];

     }*/

     inline void GetBackgroundPlasmaVelocity(double *vel, double Time = NAN) {
       double t[3];
       int idim,iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       for (idim=0;idim<3;idim++) vel[idim]=0.0;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         SWMF::GetBackgroundPlasmaVelocity(t,Stencil.cell[iStencil]);
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         DATAFILE::GetBackgroundPlasmaVelocity(t,Stencil.cell[iStencil], Time);
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif

         for (idim=0;idim<3;idim++) vel[idim]+=Stencil.Weight[iStencil]*t[idim];
       }
     }

     inline double GetBackgroundMagneticFluxFunction(double Time = NAN) {
       double res=0.0;
       int iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         res+=DATAFILE::GetBackgroundMagneticFluxFunction(Stencil.cell[iStencil], Time)*Stencil.Weight[iStencil];
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif
       }
       return res;
     }

     inline double GetBackgroundPlasmaPressure(double Time = NAN) {
       double res=0.0;
       int iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         res+=SWMF::GetBackgroundPlasmaPressure(Stencil.cell[iStencil])*Stencil.Weight[iStencil];
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         res+=DATAFILE::GetBackgroundPlasmaPressure(Stencil.cell[iStencil], Time)*Stencil.Weight[iStencil];
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif
       }

       return res;
     }

     inline double GetBackgroundElectronPlasmaPressure(double Time = NAN) {
       double res=0.0;
       int iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         exit(__LINE__,__FILE__,"Error: not implemented");
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         res+=DATAFILE::GetBackgroundElectronPlasmaPressure(Stencil.cell[iStencil], Time)*Stencil.Weight[iStencil];
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif
       }

       return res;
     }

     inline double GetBackgroundPlasmaNumberDensity(double Time = NAN) {
       double res=0.0;
       int iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         res+=SWMF::GetBackgroundPlasmaNumberDensity(Stencil.cell[iStencil])*Stencil.Weight[iStencil];
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         res+=DATAFILE::GetBackgroundPlasmaNumberDensity(Stencil.cell[iStencil], Time)*Stencil.Weight[iStencil];
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif
       }

       return res;
     }

     inline double GetBackgroundPlasmaTemperature(double Time = NAN) {
       double res=0.0;
       int iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         res+=SWMF::GetBackgroundPlasmaTemperature(Stencil.cell[iStencil])*Stencil.Weight[iStencil];
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         res+=DATAFILE::GetBackgroundPlasmaTemperature(Stencil.cell[iStencil], Time)*Stencil.Weight[iStencil];
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif
       }

       return res;
     }

     inline void GetBackgroundFieldsVector(double *E,double *B, double Time = NAN) {
       double tE[3],tB[3];
       int idim,iStencil;
       PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #else
       memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
       #endif

       for (idim=0;idim<3;idim++) E[idim]=0.0,B[idim]=0.0;

       for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
         #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
         SWMF::GetBackgroundFieldsVector(tE,tB,Stencil.cell[iStencil]);
         #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
         DATAFILE::GetBackgroundFieldsVector(tE,tB,Stencil.cell[iStencil], Time);
         #else
         exit(__LINE__,__FILE__,"not implemented");
         #endif

         for (idim=0;idim<3;idim++) {
           E[idim]+=Stencil.Weight[iStencil]*tE[idim];
           B[idim]+=Stencil.Weight[iStencil]*tB[idim];
         }
       }
     }


  }

  //prepopulate the domain
  namespace InitialCondition {
    //constant number density
    long int PrepopulateDomain(int spec,double NumberDensity,double *Velocity,double Temperature);
    // put a single particle (for each thread)
    long int PutParticle(int spec, double *x, double *v);
  }

  //Save and read the restart files
  namespace Restart {
    //sampled data
    namespace SamplingData {
      extern char RestartFileName[_MAX_STRING_LENGTH_PIC_]; //<- SAVE INTO THIS FILE: the name of the sampled restart file that will be used in the 'read' procedures to recoved the sampled data. The name is set in AMPS' input file

      void Save(const char*);
      void SaveBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,CMPI_channel*,FILE*);

      void Read(const char*);
      void ReadBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,FILE*);

      //Read Sampling data range
      //the new function and variables are placed into the namespace. Old variables are
      extern int minReadFileNumber,maxReadFileNumber;

      //preplot the recovered data file
      extern bool PreplotRecoveredData;

      //pointer to a user-defined function that defined the all combinations of the species number/output file number that will be recovered
      typedef void (*fDataRecoveryManager) (list<pair<string,list<int> > >&,int,int);    //(list<pair<string FileName,list<int> Species > >&, min output file number, max output file number);
      extern fDataRecoveryManager DataRecoveryManager;
    }

    //Particle data
    extern int ParticleRestartAutosaveIterationInterval;
    extern char saveParticleDataRestartFileName[_MAX_STRING_LENGTH_PIC_]; //the file name in which the particle data WILL BE WRITTEN
    extern char recoverParticleDataRestartFileName[_MAX_STRING_LENGTH_PIC_]; //the file name FROM WHICH the particle data will be read
    extern bool ParticleDataRestartFileOverwriteMode;

    //save used-defined data in the partcle data restart file
    const int UserAdditionalRestartDataCompletedMarkerLength=43;
    extern char UserAdditionalRestartDataCompletedMarker[UserAdditionalRestartDataCompletedMarkerLength]; //="PARTICLE-RESTART-FILE--END-USER-ADDITIONAL";

    typedef void (*fUserAdditionalRestartData)(FILE*);
    extern fUserAdditionalRestartData UserAdditionalRestartDataSave;
    extern fUserAdditionalRestartData UserAdditionalRestartDataRead;
    void SetUserAdditionalRestartData(fUserAdditionalRestartData fRead,fUserAdditionalRestartData fSave);

    void SaveParticleData(const char*);
    void SaveParticleDataBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,CMPI_channel*,FILE*);

    void ReadParticleData(const char*);
    void ReadParticleDataBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,FILE*);

    //calcualte the check sum of the save/read particle data
    unsigned long GetParticleDataCheckSum();
    void GetParticleDataBlockCheckSum(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,CRC32* CheckSum,int &PrevNodeThread);
  }

  //chemical reactions
  namespace ChemicalReactions {

    namespace PhotolyticReactions {
      //the photolytic reactions: the model is initialized if _PIC_PHOTOLYTIC_REACTIONS_MODE_ = _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
      //the total photolytic lifetime of the species;
      extern double *ConstantTotalLifeTime;

/*
      typedef double (*fTotalLifeTime)(double *x,int spec,long int ptr,bool &ReactionAllowedFlag);
      extern fTotalLifeTime *TotalLifeTime;

      typedef int (*fReactionProcessor)(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData);
      extern fReactionProcessor *ReactionProcessorTable;
*/

      void Init();
//      void SetReactionProcessor(fReactionProcessor f,int spec);
//      void SetSpeciesTotalPhotolyticLifeTime(fTotalLifeTime f,int spec);

      //The return codes of the photolytic model
      #define _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_      0
      #define _PHOTOLYTIC_REACTION_OCCURES_                  1

      #define _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_        2
      #define _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_ 3

      //the default function that returns the constant life time value
      double TotalLifeTime_default(double *x,int spec,long int ptr,bool &ReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);


      //multiply the lifetime by the following constant (use it for example for adjectment to variation of a heliocentric distance)
      static const double ReactionLifetimeMultiplier=1.0;

      //the total number of the reactions that are modeled in the current model run
      static const int nTotalUnimolecularReactions=1;

      //the maximum number of the reaction products in the simulated reactions
      static const int nMaxUnimolecularReactionProducts=1;

      //the descriptor of the Unimolecular reactions
      struct cUnimoleculecularReactionDescriptor {
//      public:
//        cUnimoleculecularReactionDescriptor() {}

        int ReactionId;
        double LifeTime;
        double ReactionRate;
        int SourceSpecie;
        int nProducts;
        int ProductSpecies[nMaxUnimolecularReactionProducts];
      };

      //the array of descriptors of the reactions
      static const cUnimoleculecularReactionDescriptor UnimoleculecularReactionDescriptor[nTotalUnimolecularReactions]={0,0.0,0.0,0,0,{0}};

      //the maximum number of reactions in which a species can particiepate
      static const int nMaxSpeciesUnimolecularReactionNumber=1;

      //the list of the reactions in wich a space can particcipate
      static const int SpeciesUnimolecularReactionList[PIC::nTotalSpecies][nMaxSpeciesUnimolecularReactionNumber]={{0}};

      //the total reaction rate for individual species
      static const double TotalSpecieUnimolecularReactionRate[PIC::nTotalSpecies]={0.0};


      void InitProductStatWeight();

      //the default function for processing the photolytic reactions -> the particle is removed if the reaction occured
      void PhotolyticReactionProcessor_default(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
      /*{
        return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
      }
*/
      //the manager of the photolytic reaction module
      //if the reaction has occured-> spes returns the species number of the transformed particles, TimeInterval return the time when the reaction had occured

       int PhotolyticReaction(double *x,long int ptr,int &spec,double &TimeInterval,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
/*      inline int PhotolyticReaction(double *x,long int ptr,int &spec,double &TimeInterval) {
        int code=_PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
        register double p,lifetime,c;
        bool flag;

        lifetime=_PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(x,spec,ptr,flag);
        if (flag==false) return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;


        c=exp(-TimeInterval/lifetime);
        p=1.0-c; //the probability for reaction to occur

        if (rnd()<p) {
          //determine the time offset when the reaction had occured
          TimeInterval=-lifetime*log(1.0+rnd()*(c-1.0));

          return _PHOTOLYTIC_REACTION_OCCURES_;
        }

        return code;
      }*/

       //execute the photolytic chemistry model
       void ExecutePhotochemicalModel();

    }

    namespace GenericParticleTranformation {
      //contains functions that are used to describe transformations (changing of internal parameters of a particle) that are not due to chemical reactions
      //dt <- can be limited by the function
      typedef int (*fTransformationIndicator)(double *x,double *v,int spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double &dt,bool &TransformationTimeStepLimitFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
//      extern fTransformationIndicator *TransformationIndicator;

      typedef int (*fTransformationProcessor)(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
//      extern fTransformationProcessor *TransformationProcessor;

      inline void Init() {
        /*
        TransformationIndicator=new fTransformationIndicator[PIC::nTotalSpecies];
        TransformationProcessor=new fTransformationProcessor[PIC::nTotalSpecies];

        //only one transformation model can be used!!!
        if (_PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_) exit(__LINE__,__FILE__,"Error: only one particle transformation model can be used");

        for (int s=0;s<PIC::nTotalSpecies;s++) TransformationIndicator[s]=NULL,TransformationProcessor[s]=NULL;
        */
      }

      inline void SetSpeciesModel(fTransformationIndicator Indicator, fTransformationProcessor Processor,int spec) {
        /*
         if (TransformationIndicator==NULL) Init();

         TransformationIndicator[spec]=Indicator;
         TransformationProcessor[spec]=Processor;
         */
      }

/*
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_       0
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_  1
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_        2
*/
    }
  }


  //process model particles with a user-defined functions
  namespace UserParticleProcessing {
    void Processing_default(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

    //call the processing function
    void Processing();
  }


  namespace BC {

    //the list of blocks that are connected to the bounding box, where the injection boundary conditions are applied
    extern list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > boundingBoxInjectionBlocksList;

    typedef bool (*fBlockInjectionIndicator)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    extern fBlockInjectionIndicator BlockInjectionBCindicatior;

    typedef long int (*fBlockInjectionBC)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    extern fBlockInjectionBC userDefinedBoundingBlockInjectionFunction;

    //the number of injected particles and injection rate
    extern long int nTotalInjectedParticles;
    extern long int *nInjectedParticles;
    extern double *ParticleProductionRate,*ParticleMassProductionRate;

    void InitBoundingBoxInjectionBlockList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=PIC::Mesh::mesh.rootTree);

    //model the particle injection for the current time step
    void InjectionBoundaryConditions();

    //calculate of the injection rate of particles distributed with Maxwellian distribution
    double CalculateInjectionRate_MaxwellianDistribution(const double NumberDesnity,const double Temp,const double *BulkVelocity,double *ExternalNormal,const int spec);

    //user-defined generic particle-injection function
    typedef long int (*fUserDefinedParticleInjectionFunction)();
    extern fUserDefinedParticleInjectionFunction UserDefinedParticleInjectionFunction;

    //the extra injection process by the exosphere model (src/models/exosphere)
    typedef long int (*fExosphereModelExtraInjectionFunction)();
    extern fExosphereModelExtraInjectionFunction ExosphereModelExtraInjectionFunction;

    namespace ExternalBoundary {

       //get direction of the gas flow at the boundary of the computational domain
      int ExternalBoundaryFlowDirection(int spec,int nface, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

      //free flow boundary (particle injection from the external boundary with the macroscopic flow conditions at the internal side of the boundary)
      namespace OpenFlow {
        extern bool BoundaryInjectionFlag[PIC::nTotalSpecies];
        extern list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > InjectionBlocksList;
        extern bool InjectionBlocksListInitFlag;

        //init the list of the injection blocks
        void InitBlockList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=PIC::Mesh::mesh.rootTree);

        //particle injection functions
        int InjectBlock(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,int nInjectionFace=-1);
        void Inject();
      }

    }

    namespace InternalBoundary {

      namespace Sphere {
        extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;
        extern int sampledFluxDownRelativeOffset,sampledFluxUpRelativeOffset;
        extern int sampledMeanVelocityDownRelativeOffset,sampledMeanVelocityUpRelativeOffset;
        extern int sampledMeanEnergyDownRelativeOffset,sampledMeanEnergyUpRelativeOffset;
        extern int sampledSurfaceNumberDensityRelativeOffset;

        extern long int TotalSampleSetLength;
        extern long int *SpeciesSampleDataOffset;
        extern long int *SpeciesSampleUserDefinedDataOffset;
        extern long int TotalSurfaceElementNumber;

        extern bool UserDefinedSamplingProcedureFlag;


        extern cAMRheap<cInternalSphericalData> InternalSpheres;

        //default sampling:
        //1. particle flux Down/Up
        //2. mean particles' velocity Down/Up
        //3. mean particles' energy Down/Up
        //4. surface number density

        //the spherical body as a internal body
        /*
        class cSurfaceDataSphere  {
        public:
          unsigned int faceat;
          double *SamplingBuffer;

          cSurfaceDataSphere() {
            faceat=-1,SamplingBuffer=NULL;
          }
        };
        */




        void Init();
        void Init(long int *RequestedSamplingSetDataLength,long int *UserDefinedSampleDataRelativeOffset);

        cInternalBoundaryConditionsDescriptor RegisterInternalSphere();
//        cSurfaceDataSphere* GetSphereSurfaceData(cInternalBoundaryConditionsDescriptor);



        inline double* GetCompletedSamplingBuffer(cInternalSphericalData* Sphere) {return Sphere->SamplingBuffer+completedCellSampleDataPointerOffset;}
        inline double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {return GetCompletedSamplingBuffer((cInternalSphericalData*)Descriptor.BoundaryElement);}

        inline double* GetCollectingSamplingBuffer(cInternalSphericalData* Sphere) {return Sphere->SamplingBuffer+collectingCellSampleDataPointerOffset;}
        inline double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {return GetCollectingSamplingBuffer((cInternalSphericalData*)Descriptor.BoundaryElement);}


        //====================================================
        //the offset of sampling data for a particular specie
        inline int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
          return 2*SurfaceElement*TotalSampleSetLength+completedCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
        }

        inline int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
          return 2*SurfaceElement*TotalSampleSetLength+collectingCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
        }




//        double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor);
//        double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor);

        //the offset of the sampling data for a particular specie
//        int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement);
//        int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int SurfaceElementSamplingSetLength();
        void switchSamplingBuffers();

        //print the 'USER DEFINED' surface data
        typedef void (*fPrintVariableList)(FILE*);
        extern fPrintVariableList PrintUserDefinedVariableList;

        typedef void (*fPrintTitle)(FILE*);
        extern fPrintTitle PrintUserDefinedTitle;

        typedef void (*fPrintDataStateVector)(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
        extern fPrintDataStateVector PrintUserDefinedDataStateVector;

        //print default surface data (3D)
        void PrintDefaultVariableList(FILE*);
        void PrintDefaultTitle(FILE*);
        void PrintDefaultDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);

        //clear the sampling buffers
        void flushCollectingSamplingBuffer(cInternalSphericalData* Sphere);

        //particle-spherical surface interaction
        typedef int (*fParticleSphereInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);
        int ParticleSphereInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);


        //Sampling of the particles data
        typedef void (*fSampleParticleData)(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        extern fSampleParticleData SampleParticleData;
      }

      namespace RotationBody {
      using namespace Sphere;
        extern cAMRheap<cInternalRotationBodyData> InternalRotationBody;

        cInternalBoundaryConditionsDescriptor RegisterInternalRotationBody();
        void PrintDefaultDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalRotationBodyData *RotationBody,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
      }

      namespace NastranSurface {
      using namespace Sphere;
        extern cAMRheap<cInternalNastranSurfaceData> InternalNastranSurface;

        cInternalBoundaryConditionsDescriptor RegisterInternalNastranSurface();
        void PrintDefaultDataStateVector(FILE* fout,long int nElement,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalNastranSurfaceData *NastranSurface,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
      }

      namespace Circle {
        extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;
        extern int sampledFluxDownRelativeOffset,sampledFluxUpRelativeOffset;
        extern int sampledMeanVelocityDownRelativeOffset,sampledMeanVelocityUpRelativeOffset;
        extern int sampledMeanEnergyDownRelativeOffset,sampledMeanEnergyUpRelativeOffset;
        extern int sampledSurfaceNumberDensityRelativeOffset;

        extern long int TotalSampleSetLength;
        extern long int *SpeciesSampleDataOffset;
        extern long int *SpeciesSampleUserDefinedDataOffset;
        extern long int TotalSurfaceElementNumber;

        extern bool UserDefinedSamplingProcedureFlag;


        extern cAMRheap<cInternalCircleData> InternalCircles;

        //default sampling:
        //1. particle flux Down/Up
        //2. mean particles' velocity Down/Up
        //3. mean particles' energy Down/Up
        //4. surface number density

        void Init();
        void Init(long int *RequestedSamplingSetDataLength,long int *UserDefinedSampleDataRelativeOffset);

        cInternalBoundaryConditionsDescriptor RegisterInternalCircle();
        double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor);
        double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor);

        //the offset of the sampling data for a particular specie
        int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int SurfaceElementSamplingSetLength();
        void switchSamplingBuffers();

        //print the 'USER DEFINED' surface data
        typedef void (*fPrintVariableList)(FILE*);
        extern fPrintVariableList PrintUserDefinedVariableList;

        typedef void (*fPrintTitle)(FILE*);
        extern fPrintTitle PrintUserDefinedTitle;

        typedef void (*fPrintDataStateVector)(FILE* fout,long int nPolarPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalCircleData *Circle,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
        extern fPrintDataStateVector PrintUserDefinedDataStateVector;

        //print default surface data (3D)
        void PrintDefaultVariableList(FILE*);
        void PrintDefaultTitle(FILE*);
        void PrintDefaultDataStateVector(FILE* fout,long int nPolarPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalCircleData *Circle,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);

        //clear the sampling buffers
        void flushCollectingSamplingBuffer(cInternalCircleData* Sphere);

        //particle-spherical surface interaction
        typedef int (*fParticleCircleInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);
        int ParticleCircleInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);


        //Sampling of the particles data
        typedef void (*fSampleParticleData)(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        extern fSampleParticleData SampleParticleData;
      }

      namespace Sphere_1D {
        extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;
        extern int sampledFluxDownRelativeOffset,sampledFluxUpRelativeOffset;
        extern int sampledMeanVelocityDownRelativeOffset,sampledMeanVelocityUpRelativeOffset;
        extern int sampledMeanEnergyDownRelativeOffset,sampledMeanEnergyUpRelativeOffset;
        extern int sampledSurfaceNumberDensityRelativeOffset;

        extern long int TotalSampleSetLength;
        extern long int *SpeciesSampleDataOffset;
        extern long int *SpeciesSampleUserDefinedDataOffset;
        extern long int TotalSurfaceElementNumber;

        extern bool UserDefinedSamplingProcedureFlag;


        extern cAMRheap<cInternalSphere1DData> InternalSpheres;

        //default sampling:
        //1. particle flux Down/Up
        //2. mean particles' velocity Down/Up
        //3. mean particles' energy Down/Up
        //4. surface number density

        //the spherical body as a internal body
        /*
        class cSurfaceDataSphere  {
        public:
          unsigned int faceat;
          double *SamplingBuffer;

          cSurfaceDataSphere() {
            faceat=-1,SamplingBuffer=NULL;
          }
        };
        */




        void Init();
        void Init(long int *RequestedSamplingSetDataLength,long int *UserDefinedSampleDataRelativeOffset);

        cInternalBoundaryConditionsDescriptor RegisterInternalSphere();
//        cSurfaceDataSphere* GetSphereSurfaceData(cInternalBoundaryConditionsDescriptor);
        double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor);
        double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor);

        //the offset of the sampling data for a particular specie
        int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int SurfaceElementSamplingSetLength();
        void switchSamplingBuffers();

        //print the 'USER DEFINED' surface data
        typedef void (*fPrintVariableList)(FILE*);
        extern fPrintVariableList PrintUserDefinedVariableList;

        typedef void (*fPrintTitle)(FILE*);
        extern fPrintTitle PrintUserDefinedTitle;

        typedef void (*fPrintDataStateVector)(FILE* fout,cInternalSphere1DData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
        extern fPrintDataStateVector PrintUserDefinedDataStateVector;

        //print default surface data
        void PrintDefaultVariableList(FILE*);
        void PrintDefaultTitle(FILE*);
        void PrintDefaultDataStateVector(FILE* fout,cInternalSphere1DData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);

        //clear the sampling buffers
        void flushCollectingSamplingBuffer(cInternalSphericalData* Sphere);

        //particle-spherical surface interaction
        typedef int (*fParticleSphereInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);
        int ParticleSphereInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);

        //Sampling of the particles data
        typedef void (*fSampleParticleData)(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        extern fSampleParticleData SampleParticleData;
      }
    }
  }

}


#endif

//include headers for individual physical models
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_ 
#include "Dust.h"
#endif

/* the headers are added by the configuration scripts (exosphere.pl) during compiling of the code
//inlude headers for the user defined models
#if _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_ == _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE__ON_
#include "UserDefinition.PIC.PhysicalModelHeaderList.h"
#elif _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_ == _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE__OFF_
//do nothing
#else
!!!! ERROR:  The option is unknown !!!!!!
#endif
*/



/*
 * ProtostellarNebula.h
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */
//$Id$

#ifndef _PROTOSTELLARNEBULA_H_
#define _PROTOSTELLARNEBULA_H_

#include <math.h>

#define _SEP_MOVER_DEFUALT_              0
#define _SEP_MOVER_BOROVIKOV_2019_ARXIV_ 1
#define _SEP_MOVER_HE_2019_AJL_          2 
#define _SEP_MOVER_KARTAVYKH_2016_AJ_    3 

#ifndef _SEP_MOVER_
#define _SEP_MOVER_ _SEP_MOVER_DEFUALT_ 
#endif

#ifndef _SEP_MOVER_DRIFT_ 
#define _SEP_MOVER_DRIFT_ _PIC_MODE_OFF_
#endif

#define _DOMAIN_GEOMETRY_PARKER_SPIRAL_ 0
#define _DOMAIN_GEOMETRY_BOX_    1

#ifndef _DOMAIN_GEOMENTRY_
#define _DOMAIN_GEOMENTRY_ _DOMAIN_GEOMETRY_PARKER_SPIRAL_ 
#endif 

#ifndef _DOMAIN_SIZE_
#define _DOMAIN_SIZE_ 250.0*_RADIUS_(_SUN_)/_AU_; 
#endif


#include "pic.h"
#include "specfunc.h"

#include "Exosphere.h"
#include "constants.h"


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif


//define which diffution model is used in the simulation
#define _DIFFUSION_NONE_                 0
#define _DIFFUSION_ROUX2004AJ_           1 
#define _DIFFUSION_BOROVIKOV_2019_ARXIV_ 2
#define _DIFFUSION_JOKIPII1966AJ_        3


#ifndef _SEP_DIFFUSION_MODEL_
#define _SEP_DIFFUSION_MODEL_  _DIFFUSION_NONE_
#endif 



namespace SEP {
  using namespace Exosphere;

  //functions used to sample and output macroscopic somulated data into the AMPS' output file
  namespace OutputAMPS {
    namespace SamplingParticleData {
      extern int DriftVelocityOffset;
      extern int absDriftVelocityOffset;
      extern int NumberDensity_PlusMu,NumberDensity_MinusMu;

      void PrintVariableList(FILE* fout,int DataSetNumber);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node);
      int RequestSamplingData(int offset);

      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

      void Init();
    }
  }



  //functions related to tracing SEPs along field lines 
  namespace FieldLine {
    long int InjectParticlesSingleFieldLine(int spec,int iFieldLine);
    long int InjectParticles();
  }

  //the namespace contains the diffution models
  namespace Diffusion {
    //LeRoux-2004-AJ
    namespace Roux2004AJ {
      void GetPitchAngleDiffusionCoefficient(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment); 
    }

    //Borovikov-2019-ARXIV (Eq. 6.11)_
    namespace Borovokov_2019_ARXIV {
      void GetPitchAngleDiffusionCoefficient(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment);
    }     

    namespace Jokopii1966AJ {
      void GetPitchAngleDiffusionCoefficient(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment);
    }
 
  }
    
  //functions used for the particle samplein
  namespace Sampling {

    const int SamplingHeliocentricDistanceTableLength=6;
    const double SamplingHeliocentricDistanceTable[]={16.0*_RADIUS_(_SUN_),0.2*_AU_,0.4*_AU_,0.6*_AU_,0.8*_AU_,1.0*_AU_};
    const double MinSampleEnergy=1.0*MeV2J;
    const double MaxSampleEnergy=500.0*MeV2J;
    const int nSampleIntervals=7;


    class cSamplingBuffer {
    public:
      int nEnergyBins;
      double MinEnergy,MaxEnergy,dLogEnergy;
     
      double *DensitySamplingTable;
      double *FluxSamplingTable,*ReturnFluxSamplingTable;
      int SamplingCounter;
      double SamplingTime;

      int iFieldLine;
      double HeliocentricDisctance;
      FILE *foutDensity,*foutFlux,*foutReturnFlux;

      PIC::FieldLine::cFieldLineSegment* GetFieldLineSegment() {
        namespace FL=PIC::FieldLine;

        double xBegin[3],xEnd[3]; 
        double rBegin,rEnd;
        int iSegment=0; //used for debugging only 

        if (iFieldLine>=FL::nFieldLine) exit(__LINE__,__FILE__,"Error: the filed line is out of range");

        FL::cFieldLineSegment* Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment(); 

        while (Segment!=NULL) {
          Segment->GetBegin()->GetX(xBegin);
          Segment->GetEnd()->GetX(xEnd);

          rBegin=Vector3D::Length(xBegin);
          rEnd=Vector3D::Length(xEnd);

          if ( ((rBegin<=HeliocentricDisctance)&&(HeliocentricDisctance<=rEnd)) || ((rBegin>=HeliocentricDisctance)&&(HeliocentricDisctance>=rEnd)) ) {
            break;
          } 
        
          iSegment++;
          Segment=Segment->GetNext();
        }

        return Segment;
      }

      void Sampling () {
        namespace FL=PIC::FieldLine;
        namespace PB=PIC::ParticleBuffer;
    
        FL::cFieldLineSegment* Segment=GetFieldLineSegment();
        if (Segment==NULL) return;

        //find cell;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
        int i,j,k;
        double xBegin[3],xEnd[3],xMiddle[3];
        long int ptr;

        Segment->GetBegin()->GetX(xBegin);
        Segment->GetEnd()->GetX(xEnd);

        for (int idim=0;idim<3;idim++) xMiddle[idim]=0.5*(xBegin[idim]+xEnd[idim]); 
   
        node=PIC::Mesh::Search::FindBlock(xMiddle);
        
        if (node==NULL) return;
        if (node->block==NULL) return;

        double vol=0.0;
        double xFirstFieldLine[3];

        switch (_PIC_PARTICLE_LIST_ATTACHING_) {
        case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
          PIC::Mesh::mesh->fingCellIndex(xMiddle,i,j,k,node);
          ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          vol=(node->xmax[0]-node->xmin[0])*(node->xmax[1]-node->xmin[1])*(node->xmax[2]-node->xmin[2])/(_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);

          break;
        case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
          ptr=Segment->FirstParticleIndex;

          FL::FieldLinesAll[iFieldLine].GetFirstSegment()->GetBegin()->GetX(xFirstFieldLine);
          vol=pow(Vector3D::Length(xMiddle)/Vector3D::Length(xFirstFieldLine),2)*Segment->GetLength(); 

          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is unknown");
        }
 
        SamplingTime+=node->block->GetLocalTimeStep(0);
        SamplingCounter++;

        while (ptr!=-1) {
          double v[3],e,m0;
          int spec,ibin;
          PB::byte* ParticleData;

          ParticleData=PB::GetParticleDataPointer(ptr); 
          spec=PB::GetI(ParticleData);

          //v=PB::GetV(ParticleData);

          v[0]=PB::GetVParallel(ParticleData);
          v[1]=PB::GetVNormal(ParticleData);
          v[2]=0.0;


          m0=PIC::MolecularData::GetMass(spec);

          if (PB::GetFieldLineId(ParticleData)!=iFieldLine) {
            ptr=PB::GetNext(ParticleData);
            continue;
          } 
          
          switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
          case _PIC_MODE_OFF_: 
            e=m0*Vector3D::DotProduct(v,v)/2.0;
            break;

          case _PIC_MODE_ON_: 
            e=Relativistic::Speed2E(Vector3D::Length(v),m0);
            break;
          }

          ibin=(int)(log(e/MinEnergy)/dLogEnergy);  

          if ((ibin>=0)&&(ibin<nEnergyBins)) {  
            double ParticleWeight;

            ParticleWeight=node->block->GetLocalParticleWeight(spec); 
            ParticleWeight*=PB::GetIndividualStatWeightCorrection(ParticleData);

            DensitySamplingTable[ibin]+=ParticleWeight/vol;
            FluxSamplingTable[ibin]+=ParticleWeight/vol*Vector3D::Length(v); 

            if (v[0]<0.0) ReturnFluxSamplingTable[ibin]+=ParticleWeight/vol*Vector3D::Length(v);
          } 

          ptr=PB::GetNext(ParticleData);
        }
      }

      void Clear() {
        for (int i=0;i<nEnergyBins;i++) DensitySamplingTable[i]=0.0,FluxSamplingTable[i]=0.0,ReturnFluxSamplingTable[i]=0.0;
        
        SamplingCounter=0;
      }

      void Output() {
        fprintf(foutDensity,"%e ",SamplingTime);

        for (int i=0;i<nEnergyBins;i++) fprintf(foutDensity,"  %e", DensitySamplingTable[i]/((SamplingCounter!=0) ? SamplingCounter : 1)); 
        
        fprintf(foutDensity,"\n");
        fflush(foutDensity);

        fprintf(foutFlux,"%e ",SamplingTime);

        for (int i=0;i<nEnergyBins;i++) fprintf(foutFlux,"  %e", FluxSamplingTable[i]/((SamplingCounter!=0) ? SamplingCounter : 1));

        fprintf(foutFlux,"\n");
        fflush(foutFlux);

        fprintf(foutReturnFlux,"%e ",SamplingTime);

        for (int i=0;i<nEnergyBins;i++) fprintf(foutReturnFlux,"  %e", ReturnFluxSamplingTable[i]/((SamplingCounter!=0) ? SamplingCounter : 1));

        fprintf(foutReturnFlux,"\n");
        fflush(foutReturnFlux);



        Clear();
      }


      //full name of the output file is saved here for debugging purposes 
      char full_name[200];

      void Init(const char *fname,double e_min,double e_max,int n,double r,int l) {
        nEnergyBins=n;
        MinEnergy=e_min,MaxEnergy=e_max;
        dLogEnergy=log(MaxEnergy/MinEnergy)/nEnergyBins; 
        HeliocentricDisctance=r;
        iFieldLine=l; 

        sprintf(full_name,"%s.density.field-line=%ld.r=%e.dat",fname,l,r/_AU_);
        foutDensity=fopen(full_name,"w"); 
        if (foutDensity==NULL) exit(__LINE__,__FILE__,"Error: cannot open a file for writting"); 

        fprintf(foutDensity,"VARIABLES=\"time\"");
        for (int i=0;i<nEnergyBins;i++) fprintf(foutDensity,", \"E(%e MeV - %e MeV)\"",MinEnergy*exp(i*dLogEnergy)*J2MeV,MinEnergy*exp((i+1)*dLogEnergy)*J2MeV);
        fprintf(foutDensity,"\n");   

        sprintf(full_name,"%s.flux.field-line=%ld.r=%e.dat",fname,l,r/_AU_);
        foutFlux=fopen(full_name,"w");
        if (foutFlux==NULL) exit(__LINE__,__FILE__,"Error: cannot open a file for writting");

        fprintf(foutFlux,"VARIABLES=\"time\"");
        for (int i=0;i<nEnergyBins;i++) fprintf(foutFlux,", \"E(%e MeV - %e MeV)\"",MinEnergy*exp(i*dLogEnergy)*J2MeV,MinEnergy*exp((i+1)*dLogEnergy)*J2MeV);
        fprintf(foutFlux,"\n");

        sprintf(full_name,"%s.return_flux.field-line=%ld.r=%e.dat",fname,l,r/_AU_);
        foutReturnFlux=fopen(full_name,"w");
        if (foutReturnFlux==NULL) exit(__LINE__,__FILE__,"Error: cannot open a file for writting");

        fprintf(foutReturnFlux,"VARIABLES=\"time\"");
        for (int i=0;i<nEnergyBins;i++) fprintf(foutReturnFlux,", \"E(%e MeV - %e MeV)\"",MinEnergy*exp(i*dLogEnergy)*J2MeV,MinEnergy*exp((i+1)*dLogEnergy)*J2MeV);
        fprintf(foutReturnFlux,"\n");

        SamplingTime=0.0;
        DensitySamplingTable=new double[nEnergyBins];
        FluxSamplingTable=new double[nEnergyBins];
        ReturnFluxSamplingTable=new double[nEnergyBins];
        Clear(); 
      }
    };

    extern int SamplingBufferTableLength;
    extern cSamplingBuffer **SamplingBufferTable; 

    //Manager is called by AMPS to perform sampling procedure 
    void Manager();

    //Init the samping module
    void Init();
    void InitSingleFieldLineSampling(int iFieldLine); 
  }

  //sphere describing the inner boundary of the domain 
  extern cInternalSphericalData* InnerBoundary;

  //parameters controlling the model execution
  const int DomainType_ParkerSpiral=0;
  const int DomainType_MultipleParkerSpirals=1; 
  const int DomainType_FLAMPA_FieldLines=2;
  extern int DomainType;
  extern int Domain_nTotalParkerSpirals;

  const int ParticleTrajectoryCalculation_GuidingCenter=0;
  const int ParticleTrajectoryCalculation_RelativisticBoris=1;
  const int ParticleTrajectoryCalculation_IgorFieldLine=2;
  const int ParticleTrajectoryCalculation_FieldLine=3;

  extern int ParticleTrajectoryCalculation;

  //calcualtion of the drift velocity
  extern int b_times_grad_absB_offset;
  extern int CurlB_offset;
  extern int b_b_Curl_B_offset;

  //offsets of the momentum and cos(pitch angle) in particle state vector
  namespace Offset {
    extern int Momentum,CosPitchAngle;
    extern int p_par,p_norm;
  }

  //request data in the particle state vector
  void RequestParticleData();

  int RequestStaticCellData(int);
  void GetDriftVelocity(double *v_drift,double *x,double v_parallel,double v_perp,double ElectricCharge,double mass,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* Node);
  void GetDriftVelocity(double *v_drift,double *x,double v_parallel,double v_perp,double ElectricCharge,double mass,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* Node,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil);
  void InitDriftVelData();


  int ParticleMover_HE_2019_AJL(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
  int ParticleMover_BOROVIKOV_2019_ARXIV(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
  int ParticleMover_default(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
  int ParticleMover__He_2019_AJL(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
  int ParticleMover_Kartavykh_2016_AJ(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node); 
  int ParticleMover_Droge_2009_AJ(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

  //particle mover
  int inline ParticleMover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    int res;


    //init drift velocity data if needed 
    if (_SEP_MOVER_DRIFT_==_PIC_MODE_ON_) {
      static double last_swmf_coupling_time=-10.0;

      if (last_swmf_coupling_time!=PIC::CPLR::SWMF::CouplingTime) {
        last_swmf_coupling_time=PIC::CPLR::SWMF::CouplingTime; 

        SEP::InitDriftVelData();
      }
    } 


    switch(_SEP_MOVER_) {
    case _SEP_MOVER_DEFUALT_:
      res=ParticleMover_default(ptr,dtTotal,startNode);
      break;
    case _SEP_MOVER_BOROVIKOV_2019_ARXIV_:
      res=ParticleMover_BOROVIKOV_2019_ARXIV(ptr,dtTotal,startNode);
      break;
    case _SEP_MOVER_HE_2019_AJL_:
      res=ParticleMover__He_2019_AJL(ptr,dtTotal,startNode);
      break; 
    case _SEP_MOVER_KARTAVYKH_2016_AJ_:
      res=ParticleMover_Kartavykh_2016_AJ(ptr,dtTotal,startNode);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not found");
    }

    return res;
  }

  namespace Sampling {
    using namespace Exosphere::Sampling;
  }

  namespace ParticleSource {
    namespace ShockWave {
      bool IsShock(PIC::Mesh::cDataCenterNode *CenterNode);

      extern int ShockStateFlag_offset;

      namespace Output {
        void PrintVariableList(FILE* fout,int DataSetNumber);
        void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
        void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      }
    }

    namespace OuterBoundary {
      bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
      long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
      long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

      double BoundingBoxInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    }

    namespace InnerBoundary {
      double sphereInjectionRate(int spec,int BoundaryElementType,void *BoundaryElement);
      long int sphereParticleInjection(int spec,int BoundaryElementType,void *SphereDataPointer);
      long int sphereParticleInjection(int BoundaryElementType,void *BoundaryElement);
      long int sphereParticleInjection(void *SphereDataPointer);
    }
  }


  class cFieldLine {
  public:
    double x[3],U[3],B[3],Wave[2];

    cFieldLine() {
      for (int idim=0;idim<3;idim++) x[idim]=0.0,U[idim]=0.0,B[idim]=0.0;
      for (int i=0;i<2;i++) Wave[i]=0.0;
    }
  };

    namespace ParkerSpiral {
      void GetB(double *B,double *x,double u_sw=400.0E3);
      void CreateFileLine(list<SEP::cFieldLine> *field_line,double *xstart,double length_rsun);
    }

    namespace Mesh {
      double localSphericalSurfaceResolution(double *x);
      double localResolution(double *x);
      bool NodeSplitCriterion(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

      void ImportFieldLine(list<SEP::cFieldLine> *field_line);
      void PrintFieldLine(list<SEP::cFieldLine> *field_line,const char *fname);
      void LoadFieldLine_flampa(list<SEP::cFieldLine> *field_line,const char *fname);
      
      double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

      //magnetic field line data
      extern double **FieldLineTable;
      extern int FieldLineTableLength;

    //init field line in AMPS
    void InitFieldLineAMPS(list<SEP::cFieldLine> *field_line);
    }


    namespace Scattering {
      namespace AIAA2005 {
        double MeanFreePath(PIC::ParticleBuffer::byte *ParticleData);
        void Process(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
      }
    }

    namespace Mover {

    }
}


/*
namespace ProtostellarNebula {
  using namespace Exosphere;

  namespace Sampling {
    using namespace Exosphere::Sampling;

  }

  namespace UserdefinedSoruce {
    inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
      return 1.0E20;
    }

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
      startNode=PIC::Mesh::mesh->findTreeNode(x_LOCAL_SO_OBJECT,startNode);
      if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) return false;

      //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
//      PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,ImpactVaporization_SourceTemperature[spec],ExternalNormal,spec);



//DEBUG -> injected velocity is normal to the surface

for (int i=0;i<3;i++)  v_LOCAL_IAU_OBJECT[i]=-ExternalNormal[i]*1.0E3;
//END DEBUG



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
      PIC::IDF::InitRotTemp(ImpactVaporization_SourceTemperature[spec],tempParticleData);
#endif

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM__VT_RELAXATION_MODE_  == _PIC_MODE_ON_
      exit(__LINE__,__FILE__,"Error: not implemented");
#endif

#endif

      return true;
    }

  }



//defined the forces that acts upon a particle on
#define _FORCE_GRAVITY_MODE_ _PIC_MODE_OFF_
#define _FORCE_LORENTZ_MODE_ _PIC_MODE_OFF_
#define _FORCE_FRAMEROTATION_MODE_ _PIC_MODE_OFF_

  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));

    accl[0]=0.0; accl[1]=0.0;  accl[2]=0.0; 

#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_

    exit(__LINE__,__FILE__,"ERROR: Lorentz force not implemented");

    //********************************************************
    // the following code is copied from srcEuropa/Europa.h
    //********************************************************
    long int nd;
    char *offset;
    int i,j,k;
    PIC::Mesh::cDataCenterNode *CenterNode;
    double E[3],B[3];

    if ((nd=PIC::Mesh::mesh->fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
      startNode=PIC::Mesh::mesh->findTreeNode(x_LOCAL,startNode);

      if ((nd=PIC::Mesh::mesh->fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }
    }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
#endif

    CenterNode=startNode->block->GetCenterNode(nd);
    offset=CenterNode->GetAssociatedDataBufferPointer();
    
    PIC::CPLR::GetBackgroundMagneticField(B,x_LOCAL,nd,startNode);
    PIC::CPLR::GetBackgroundElectricField(E,x_LOCAL,nd,startNode);

    double ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
    double mass=PIC::MolecularData::GetMass(spec);

    accl_LOCAL[0]+=ElectricCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
    accl_LOCAL[1]+=ElectricCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
    accl_LOCAL[2]+=ElectricCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;

#endif


#if _FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_
  //the gravity force
  double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
  double r=sqrt(r2);
  int idim;

  for (idim=0;idim<DIM;idim++) {
    accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x_LOCAL[idim]/r;
  }
#endif

#if _FORCE_FRAMEROTATION_MODE_ == _PIC_MODE_ON_
  // by default rotation period is ~25 days (freq = 4.63E-7 sec^-1)
  // frame angular velocity
  static const double Omega        = 2.0*Pi*4.63E-7;//rad/sec 
  // frame angular velocity x 2
  static const double TwoOmega     = 2.0*Omega;
  // frame angular velocity squared
  static const double SquaredOmega = Omega*Omega;
  double aCen[3],aCorr[3];

  aCen[0] = SquaredOmega * x_LOCAL[0];
  aCen[1] = SquaredOmega * x_LOCAL[1];
  aCen[2] = 0.0;


  aCorr[0] =   TwoOmega * v_LOCAL[1];
  aCorr[1] = - TwoOmega * v_LOCAL[0];
  aCorr[2] =   0.0;

  accl_LOCAL[0]+=aCen[0]+aCorr[0];
  accl_LOCAL[1]+=aCen[1]+aCorr[1];
  accl_LOCAL[2]+=aCen[2]+aCorr[2];

#endif

    //copy the local value of the acceleration to the global one
    memcpy(accl,accl_LOCAL,3*sizeof(double));
  }


  inline double ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag) {
    static const double LifeTime=3600.0*5.8/pow(0.4,2);

    // no photoionization for now
    PhotolyticReactionAllowedFlag=false;
    return -1.0;

  }

  inline int ExospherePhotoionizationReactionProcessor(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData) {

    // no photoionization for now
    return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
  }

}
*/

#endif /* PROTOSTELLARNEBULA_H_ */

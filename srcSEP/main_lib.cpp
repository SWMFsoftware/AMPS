

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//$Id$


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "constants.h"
#include "sep.h"

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
#include "amps2swmf.h"
#endif



//the parameters of the domain and the sphere

const double DebugRunMultiplier=4.0;
double rSphere=_RADIUS_(_TARGET_);


const double xMaxDomain=_DOMAIN_SIZE_*_AU_/_RADIUS_(_SUN_);

const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*4.0*1.0/100/2.5,dxMaxSphere=DebugRunMultiplier*2.0/10.0;

const double MarkNotUsedRadiusLimit=100.0;


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

  if (node->IsUsedInCalculationFlag==false) return 0.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
   //delete all particles that was not reflected on the surface
   //PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;
}



void amps_init_mesh() {
  PIC::InitMPI();

  //set the function for the particle sampling 
  SEP::Sampling::Init();
  SEP::Init();

  //set the particle injection function used in case magnetic field lines are used 
  if (_PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_) {
    PIC::BC::UserDefinedParticleInjectionFunction=SEP::FieldLine::InjectParticles;
  }

  //in case the magnetic field line is used -> reserve the nessesaty memory 
  if (_PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_) {
      using namespace PIC::FieldLine;
   
      VertexAllocationManager.MagneticField=true;
      VertexAllocationManager.ElectricField=true;
      VertexAllocationManager.PlasmaVelocity=true;
      VertexAllocationManager.PlasmaDensity=true;
      VertexAllocationManager.PlasmaTemperature=true;
      VertexAllocationManager.PlasmaPressure=true;
      VertexAllocationManager.MagneticFluxFunction=true;
      VertexAllocationManager.PlasmaWaves=true;
      VertexAllocationManager.Fluence=true;


      VertexAllocationManager.PreviousVertexData.MagneticField=true;
      VertexAllocationManager.PreviousVertexData.ElectricField=true;
      VertexAllocationManager.PreviousVertexData.PlasmaVelocity=true;
      VertexAllocationManager.PreviousVertexData.PlasmaDensity=true;
      VertexAllocationManager.PreviousVertexData.PlasmaTemperature=true;
      VertexAllocationManager.PreviousVertexData.PlasmaPressure=true;
      VertexAllocationManager.PreviousVertexData.PlasmaWaves=true;


      PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;
   }

  PIC::Init_BeforeParser();
  SEP::RequestParticleData();

  //request storage for calculating the drift velocity
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(SEP::RequestStaticCellData);

  //SetUp the alarm
//  PIC::Alarm::SetAlarm(2000);


  rnd_seed();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


  if (PIC::CPLR::SWMF::BlCouplingFlag==false) {
    //reserve data for magnetic filed
    PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;
    //PIC::CPLR::DATAFILE::Offset::MagneticField.active=true;
    //PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    //PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=PIC::CPLR::DATAFILE::Offset::MagneticField.nVars*sizeof(double);
  }

  //init the Mercury model
 ////::Init_BeforeParser();
//  PIC::Init_BeforeParser();

//  ProtostellarNebula::OrbitalMotion::nOrbitalPositionOutputMultiplier=10;
///  ProtostellarNebula::Init_AfterParser();



  //register the sphere
  {
    double sx0[3]={0.0,0.0,0.0};
    cInternalBoundaryConditionsDescriptor SphereDescriptor;
    cInternalSphericalData *Sphere;

    //correct radiust of the  sphere to be consistent with the location of the inner boundary of the SWMF/SC
    //use taht in case of coupling to the SWMF
    #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    if (AMPS2SWMF::Heliosphere::rMin>0.0) rSphere=AMPS2SWMF::Heliosphere::rMin;
    #endif

    //reserve memory for sampling of the surface balance of sticking species
    long int ReserveSamplingSpace[PIC::nTotalSpecies];

    for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


    cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);

    PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
    SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
    Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
    Sphere->SetSphereGeometricalParameters(sx0,rSphere);

    //set the innber bounday sphere 
    SEP::InnerBoundary=Sphere;

    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/Sphere.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceMesh(fname);

    sprintf(fname,"%s/SpheraData.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceData(fname,0);


    Sphere->localResolution=SEP::Mesh::localSphericalSurfaceResolution;
    Sphere->InjectionRate=SEP::ParticleSource::InnerBoundary::sphereInjectionRate;
    Sphere->faceat=0;
    Sphere->ParticleSphereInteraction=ParticleSphereInteraction;

    if (_DOMAIN_GEOMETRY_!=_DOMAIN_GEOMETRY_BOX_) {
      Sphere->InjectionBoundaryCondition=SEP::ParticleSource::InnerBoundary::sphereParticleInjection;
    }

    Sphere->PrintTitle=SEP::Sampling::OutputSurfaceDataFile::PrintTitle;
    Sphere->PrintVariableList=SEP::Sampling::OutputSurfaceDataFile::PrintVariableList;
    Sphere->PrintDataStateVector=SEP::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

    //set up the planet pointer in Mercury model
    SEP::Planet=Sphere;
    Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);
  }

  //init the solver
  PIC::Mesh::initCellSamplingDataBuffer();

  //init the mesh
  cout << "Init the mesh" << endl;

  int idim;
  double xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

  for (idim=0;idim<DIM;idim++) {
    xmax[idim]=xMaxDomain*_RADIUS_(_SUN_);
    xmin[idim]=-xMaxDomain*_RADIUS_(_SUN_);
  }

  //generate the magneric field line
  list<SEP::cFieldLine> field_line,field_line_old,field_line_new;
  double xStart[3]={1.1,0.0,0.0};


  if (PIC::CPLR::SWMF::BlCouplingFlag==false) switch (SEP::DomainType) {
  case SEP::DomainType_ParkerSpiral:
    PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;

    PIC::FieldLine::VertexAllocationManager.PlasmaWaves=true;
    PIC::FieldLine::VertexAllocationManager.MagneticField=true;
    PIC::FieldLine::VertexAllocationManager.PlasmaVelocity=true;
  
    if (SEP::Domain_nTotalParkerSpirals==1) { 
      SEP::ParkerSpiral::CreateFileLine(&field_line,xStart,250.0);
      SEP::Mesh::ImportFieldLine(&field_line);

      if (SEP::ParticleTrajectoryCalculation==SEP::ParticleTrajectoryCalculation_FieldLine) {
        PIC::FieldLine::Init();
        SEP::Mesh::InitFieldLineAMPS(&field_line);
      }
    }
    else {
      if (SEP::ParticleTrajectoryCalculation==SEP::ParticleTrajectoryCalculation_FieldLine) PIC::FieldLine::Init();

      rnd_seed(10);

      for (int iline=0;iline<SEP::Domain_nTotalParkerSpirals;iline++) {
        double x0[3],r,phi;  
        double phi_max=15.0*Pi/180.0;

        r=Vector3D::Length(xStart);

        phi=phi_max*rnd();
        if (rnd()<0.5) phi=-phi;

        x0[0]=r*sin(phi);
        x0[1]=r*cos(phi);
        x0[2]=0.0; 

        field_line.clear();

        SEP::ParkerSpiral::CreateFileLine(&field_line,x0,250.0);
        SEP::Mesh::ImportFieldLine(&field_line);

        if (SEP::ParticleTrajectoryCalculation==SEP::ParticleTrajectoryCalculation_FieldLine) SEP::Mesh::InitFieldLineAMPS(&field_line);
      } 

      if (SEP::ParticleTrajectoryCalculation==SEP::ParticleTrajectoryCalculation_FieldLine) PIC::FieldLine::Output("all-field-lines.dat",false);
    }

    break;
  case SEP::DomainType_FLAMPA_FieldLines:
    PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true; 

    SEP::Mesh::LoadFieldLine_flampa(&field_line_old,"FieldLineOld.in");
    SEP::Mesh::ImportFieldLine(&field_line_old);
    SEP::Mesh::PrintFieldLine(&field_line_old,"FieldLineOld.dat");

    SEP::Mesh::LoadFieldLine_flampa(&field_line_new,"FieldLineNew.in");
    SEP::Mesh::ImportFieldLine(&field_line_new);
    SEP::Mesh::PrintFieldLine(&field_line_new,"FieldLineNew.dat");


    //init the field line library in AMPS
    PIC::FieldLine::VertexAllocationManager.PlasmaWaves=true;
    PIC::FieldLine::VertexAllocationManager.MagneticField=true;
    PIC::FieldLine::VertexAllocationManager.PlasmaVelocity=true;

    PIC::FieldLine::DatumAtVertexPlasmaWaves.length=2;

    PIC::FieldLine::Init();

    SEP::Mesh::InitFieldLineAMPS(&field_line_old);
    SEP::Mesh::InitFieldLineAMPS(&field_line_new);
  }


  //refining the mesh along a set of magnetic field lines: use onle wher model SEP
  if (_MODEL_CASE_==_MODEL_CASE_SEP_TRANSPORT_) { 
    PIC::Mesh::mesh->UserNodeSplitCriterion=SEP::Mesh::NodeSplitCriterion;
  }


  //generate only the tree
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,SEP::Mesh::localResolution);
  PIC::Mesh::mesh->memoryAllocationReport();


/*
  if (PIC::Mesh::mesh->ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile("mesh.msh");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    MPI_Barrier(MPI_COMM_WORLD);
    PIC::Mesh::mesh->readMeshFile("mesh.msh");
  }
*/

  PIC::Mesh::mesh->buildMesh();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh->ThisThread << endl;

  //PIC::Mesh::mesh->outputMeshTECPLOT("mesh.dat");

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);
#endif

  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();



  //mark no-used black that are far from the magnetic filed line 
  list <cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> not_used_list; 
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,list <cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*>*)> MarkNotUsed; 
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> GetMaxBlockRefinmentLevel;

  int cnt=0,MaxRefinmentLevel=-1;


  GetMaxBlockRefinmentLevel=[&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) -> void {
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      double x[3];

      for (int idim=0;idim<3;idim++) x[idim]=0.5*(startNode->xmin[idim]+startNode->xmax[idim]);

       if (Vector3D::Length(x)>10.0*_RADIUS_(_SUN_)) {
         if (MaxRefinmentLevel<startNode->RefinmentLevel) MaxRefinmentLevel=startNode->RefinmentLevel;
       }
    }
    else {
      int iDownNode;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

      for (iDownNode=0;iDownNode<(1<<DIM);iDownNode++) if ((downNode=startNode->downNode[iDownNode])!=NULL) {
        GetMaxBlockRefinmentLevel(downNode);
      }
    }
  };  

  
  MarkNotUsed=[&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,list <cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> *not_used_list) -> void {
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      //if (startNode->xmin[0]<0.0) not_used_list->push_back(startNode);

      if (cnt%PIC::nTotalThreads!=PIC::ThisThread) return;

      double r,x[3],d2min=-1.0,dmax;
      int idim;

      for (idim=0;idim<3;idim++) x[idim]=0.5*(startNode->xmin[idim]+startNode->xmax[idim]);

      r=Vector3D::Length(x);

      if (r>MarkNotUsedRadiusLimit*_RADIUS_(_SUN_)) {
        list<SEP::cFieldLine>::iterator it;
        double d2,t;
        
        for (it=field_line.begin();it!=field_line.end();it++) {
          for (idim=0,d2=0.0;idim<3;idim++) {
            t=min(fabs(startNode->xmin[idim]-it->x[idim]),fabs(startNode->xmax[idim]-it->x[idim])); 
            d2+=t*t;
          }

          if ((d2min<0.0)||(d2min>d2)) d2min=d2;
        }

        dmax=MarkNotUsedRadiusLimit*_RADIUS_(_SUN_)+MarkNotUsedRadiusLimit/200.0*(r-50.0*_RADIUS_(_SUN_));

        if (d2min>dmax*dmax) {
          not_used_list->push_back(startNode); 
        }
        else {
          if (startNode->RefinmentLevel<MaxRefinmentLevel-2) not_used_list->push_back(startNode);
        }
      }
    }
    else {
      int iDownNode;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

      for (iDownNode=0;iDownNode<(1<<DIM);iDownNode++) if ((downNode=startNode->downNode[iDownNode])!=NULL) {
        MarkNotUsed(downNode,not_used_list);
      }
    }
  }; 

  GetMaxBlockRefinmentLevel(PIC::Mesh::mesh->rootTree);

  if (_DOMAIN_GEOMETRY_!= _DOMAIN_GEOMETRY_BOX_)  {
    MarkNotUsed(PIC::Mesh::mesh->rootTree,&not_used_list);

    PIC::Mesh::mesh->SetTreeNodeActiveUseFlag(&not_used_list,NULL,false,NULL);
  }

  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();

  //PIC::Mesh::mesh->outputMeshTECPLOT("mesh-reduced.dat");

  //initialize the blocks
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);
#endif

  //init the volume of the cells'
  PIC::Mesh::mesh->InitCellMeasure();


}

void amps_init() {



//init the PIC solver
  PIC::Init_AfterParser ();
	PIC::Mover::Init();


  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=SEP::Mesh::localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //create the list of mesh nodes where the injection boundary conditinos are applied
  if (_DOMAIN_GEOMETRY_==_DOMAIN_GEOMETRY_BOX_) {
    PIC::BC::BlockInjectionBCindicatior=SEP::BoundingBoxInjection::InjectionIndicator;
    PIC::BC::userDefinedBoundingBlockInjectionFunction=SEP::BoundingBoxInjection::InjectionProcessor;
    PIC::BC::InitBoundingBoxInjectionBlockList();
  }

  //set up the particle weight
  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=SEP::ParticleSource::OuterBoundary::BoundingBoxInjectionRate;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H_PLUS_SPEC_);


  //init magnetic filed
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> InitMagneticField; 

  InitMagneticField =[&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) -> void {
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      PIC::Mesh::cDataBlockAMR *block;

      if ((block=startNode->block)!=NULL) {
        double B[3],x[3];
        int idim,i,j,k,LocalCellNumber;
        PIC::Mesh::cDataCenterNode *cell;
        double *data;

        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);

          if ((cell=block->GetCenterNode(LocalCellNumber))!=NULL) {
            cell->GetX(x);
            SEP::ParkerSpiral::GetB(B,x);

            if (cell->Measure==0.0) {
              PIC::Mesh::mesh->InitCellMeasureBlock(startNode); 

              if (cell->Measure==0.0) {
                PIC::Mesh::mesh->CenterNodes.deleteElement(cell);
                startNode->block->SetCenterNode(NULL,LocalCellNumber);
                continue;
              }
            }

            data=(double*)(cell->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset); 

            for (idim=0;idim<3;idim++) data[idim]=B[idim]; 
          } 
        }
      }
    }
    else {
      int iDownNode;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

      for (iDownNode=0;iDownNode<(1<<DIM);iDownNode++) if ((downNode=startNode->downNode[iDownNode])!=NULL) {
        InitMagneticField(downNode);
      }
    }
  };

  if (_PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_) {
    InitMagneticField(PIC::Mesh::mesh->rootTree);
    //PIC::Mesh::mesh->outputMeshDataTECPLOT("magnetic-field.dat",0);
  }


  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::Mesh::mesh->ThisThread==0) cout << "The mesh is generated" << endl;

  //init the particle buffer
  PIC::ParticleBuffer::Init(10000000);
//  double TimeCounter=time(NULL);
//  int LastDataOutputFileNumber=-1;


/*  //init the sampling of the particls' distribution functions: THE DECLARATION IS MOVED INTO THE INPUT FILE
  const int nSamplePoints=3;
  double SampleLocations[nSamplePoints][DIM]={{7.6E5,6.7E5,0.0}, {2.8E5,5.6E5,0.0}, {-2.3E5,3.0E5,0.0}};

  PIC::DistributionFunctionSample::vMin=-40.0E3;
  PIC::DistributionFunctionSample::vMax=40.0E3;
  PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

  PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);*/
}







  //time step
// for (long int niter=0;niter<100000001;niter++) {
void amps_time_step(){

    //make the time advance
    static int LastDataOutputFileNumber=0;

    //make the time advance
     PIC::TimeStep();

     PIC::ParticleSplitting::Split::SplitWithVelocityShift_FL(SEP::MinParticleLimit,SEP::MaxParticleLimit);


     // write output file
     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
//       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;
       
       
       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }
     
}




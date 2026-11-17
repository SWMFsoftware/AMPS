//$Id$


#include "pic.h"
#include "constants.h"
#include "../models/surface/Surface.h"

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


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"
#include "OrbiterRuntimeInput.h"

double BulletLocalResolution(double *x) {
  int idim;
  double res,r;



return  ((fabs(x[0])<100.0)||(x[1]*x[1]+x[2]*x[2]<40.0*40.0)) ? 5.0 : 100.0;

  r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

  if (r<0.5) res=10.0;
  else if (res<5.0) res=0.24;
  else res=10.0;

  return res;
}



namespace ParticleSurfaceInterationModel {
  //diffuse re-emission from the walls at 293 K: the flux-weighted half-range Maxwellian (Surface::DiffuseReflection;
  //until 2026-11 a full Maxwellian was sampled and the particles with v.n<=0 were rejected, which gives the density, not
  //the flux, distribution: <v.n> was 1/sqrt(pi) instead of sqrt(pi)/2 times sqrt(2kT/m))
  int COPS(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    const double Twall=293.0;

    Surface::DiffuseReflection::Reemit(vInit,Twall,PIC::ParticleBuffer::GetI(ptr),TriangleCutFace);
    return _PARTICLE_REJECTED_ON_THE_FACE_; 
  }

  //the same for the CYCNSS spacecraft: the faces with the attribute 8 are at 373 K, the others at 323 K
  int ParticleSurfaceInteraction_CYCNSS(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    double Twall=50.0+273.0;

    if (TriangleCutFace->attribute==8) Twall=100.0+273.0;

    Surface::DiffuseReflection::Reemit(vInit,Twall,PIC::ParticleBuffer::GetI(ptr),TriangleCutFace);
    return _PARTICLE_REJECTED_ON_THE_FACE_;
  }
}


//the drag-coefficient test (MakefileTest/Table, test/Sphere-drag-coefficient): one case of the runtime input
//(input/test/sphere-drag-coefficient.in, "amps -i FILE -case NAME") is simulated: the free stream about the sphere with
//the gas-surface interaction model of the case (Surface::ParticleInteractionProcessor). The run ends at the first output
//after run.iterations iterations; the drag coefficient of that last sample, its statistical standard deviation and the
//input are written to PT/plots/drag-coefficient-test.dat (compared with the independent free-molecular value by
//srcOrbiter/test/sphere_drag_check.py).
namespace DragCoefficientTest {
  int ParticleSurfaceInteractionProcessor_TEST(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    return Surface::ParticleInteractionProcessor(ptr,xInit,vInit,TriangleCutFace,startNode);
  }

  void Execute() {
    int LastDataOutputFileNumber=-1,nTotalIterations=Orbiter::RuntimeInput::nTotalIterations;

    if (Orbiter::RuntimeInput::Loaded==false) exit(__LINE__,__FILE__,"Error: the drag-coefficient test requires the runtime input (amps -i FILE -case NAME)");
    if (PIC::nTotalSpecies!=1) exit(__LINE__,__FILE__,"Error: the test is set to be exdcuted only with a single model speceis");

    PIC::DataOutputFileNumber=0;

    //time step: the run ends at the first output after nTotalIterations iterations
    for (long int niter=0;;niter++) {
      PIC::TimeStep();

      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        if (niter>nTotalIterations) break;

        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>1000) PIC::RequiredSampleLength=1000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }

      if (PIC::Mesh::mesh->ThisThread==0) {
        time_t TimeValue=time(NULL);
        tm *ct=localtime(&TimeValue);

        printf(": (%i/%i %i:%i:%i), Iteration: %ld  (case %s, current sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,Orbiter::RuntimeInput::SelectedCase.c_str(),PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
      }
    }

    //output the particle statistics for the nightly tests
    char fname[400];

    sprintf(fname,"%s/test_Orbiter.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

    //the result of the case: the input and the drag coefficient of the last sample
    if (PIC::ThisThread==0) {
      sprintf(fname,"%s/drag-coefficient-test.dat",PIC::OutputDataFileDirectory);

      FILE *fout=fopen(fname,"w");
      if (fout==NULL) exit(__LINE__,__FILE__,"Error: cannot open the test result file");

      std::string gsi=Surface::RuntimeInput::Describe();

      fprintf(fout,"# The drag-coefficient test (srcOrbiter/main.cpp, DragCoefficientTest): input and result of one case\n");
      fprintf(fout,"%s",Orbiter::RuntimeInput::Describe().c_str());
      for (int s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,"model[%s]=%s\n",PIC::MolecularData::GetChemSymbol(s),Surface::Model::Name(Surface::Model::Selected[s]));
      fprintf(fout,"wall_temperature_k=%.9e\n",Surface::Temperature::Isothremal::Temp);
      fprintf(fout,"%s",Surface::RuntimeInput::ParameterLines().c_str());
      for (size_t i=0,j;i<gsi.size();i=j+1) {
        j=gsi.find('\n',i);
        if (j==std::string::npos) j=gsi.size();
        fprintf(fout,"# %s\n",gsi.substr(i,j-i).c_str());
      }
      fprintf(fout,"mpi_processes=%i\n",PIC::nTotalThreads);
      fprintf(fout,"projection_area_m2=%.9e\n",Orbiter::Sampling::DragCoefficient::LastProjectionArea);
      fprintf(fout,"sample_length=%ld\n",Orbiter::Sampling::DragCoefficient::LastSampleLength);
      fprintf(fout,"hits=%ld\n",Orbiter::Sampling::DragCoefficient::LastnHits);
      fprintf(fout,"Cd=%.9e\n",Orbiter::Sampling::DragCoefficient::LastCd);
      fprintf(fout,"sigma_Cd=%.9e\n",Orbiter::Sampling::DragCoefficient::LastCdSigma);
      fclose(fout);

      printf("$PREFIX: drag-coefficient test, case %s: Cd=%e, statistical standard deviation %e\n",Orbiter::RuntimeInput::SelectedCase.c_str(),
          Orbiter::Sampling::DragCoefficient::LastCd,Orbiter::Sampling::DragCoefficient::LastCdSigma);
    }

    MPI_Finalize();
    cout << "End of the run:" << PIC::nTotalSpecies << endl;

    exit(EXIT_SUCCESS); 
  }
}  


//add to a face data value: atomically when the OpenMP threads of the hybrid mode can update the same face
inline void AddToFaceData(double *p,double x) {
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp atomic
#endif
  *p+=x;
}

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double dt,c,vInitUnchanged[3],ParticleWeight,RateFactor,mass;
  int code,spec;

  spec=PIC::ParticleBuffer::GetI(ptr);
  mass=PIC::MolecularData::GetMass(spec);
  ParticleWeight=startNode->block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr); 

  memcpy(vInitUnchanged,vInit,3*sizeof(double));

  //process particle/surface interaction event
  code=_ORBITER__PARTICLE_SURFACE_INTERACTION_PROCESSOR_(ptr,xInit,vInit,TriangleCutFace,startNode);

  //check for adsorption of the particle on teh face 
  if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
    if (_ORBITER_ADSORPTION_MODEL__MODE_ == _ORBITER_ADSORPTION_MODEL__MODE__OFF_) {
      exit(__LINE__,__FILE__,"Error: with the surrent settings only refrection on the surface is permitted. Change the model settings");
    }

    //the particle has been adsorbed on the surface
    vInit[0]=0.0,vInit[1]=0.0,vInit[2]=0.0;  //this is needed to calculating of the momentum transfer to the surface

    //update the adsorption flux
    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    double *p=TriangleCutFace->UserData.AdsorptionFlux;

    #pragma omp atomic
    p[spec]+=ParticleWeight;

    #else
    TriangleCutFace->UserData.AdsorptionFlux[spec]+=ParticleWeight;
    #endif
  }


  //sample the energy and momentum transfer rates: one code for the MPI and the hybrid (MPI+OpenMP) modes. Each OpenMP
  //thread adds to its own element of the drag sums (they are summed in Orbiter::Sampling::DragCoefficient::
  //PrintOutputFile()), and the face data are updated atomically. Until 2026-11 the hybrid mode had its own copy with
  //the opposite sign of the momentum and twice the energy.
  dt=startNode->block->GetLocalTimeStep(spec);
  RateFactor=mass*ParticleWeight/dt; 

  int thread=0;
  double dp[3],de,dpProjection=0.0,speed;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  thread=omp_get_thread_num();
#elif _COMPILATION_MODE_ != _COMPILATION_MODE__MPI_
  exit(__LINE__,__FILE__,"Error: the option is unknown");
#endif

  Orbiter::Sampling::DragCoefficient::Flux[thread]+=ParticleWeight/dt;
  Orbiter::Sampling::DragCoefficient::ModelParticleCounter[thread]++;

  //the rate of the momentum gain by the orbiter's surface: -(v_out-v_in) m w/dt
  for (int idim=0;idim<3;idim++) dp[idim]=-(vInit[idim]-vInitUnchanged[idim])*RateFactor;

  Orbiter::Sampling::DragCoefficient::dpX[thread]+=dp[0];
  AddToFaceData(TriangleCutFace->UserData.MomentumTransferRateX+spec,dp[0]/TriangleCutFace->SurfaceArea);

  Orbiter::Sampling::DragCoefficient::dpY[thread]+=dp[1];
  AddToFaceData(TriangleCutFace->UserData.MomentumTransferRateY+spec,dp[1]/TriangleCutFace->SurfaceArea);

  Orbiter::Sampling::DragCoefficient::dpZ[thread]+=dp[2];
  AddToFaceData(TriangleCutFace->UserData.MomentumTransferRateZ+spec,dp[2]/TriangleCutFace->SurfaceArea);

  //the squared contribution to the drag (along the upstream flow) for its statistical uncertainty
  speed=sqrt(Orbiter::UpstreamBC::Velocity[0]*Orbiter::UpstreamBC::Velocity[0]+Orbiter::UpstreamBC::Velocity[1]*Orbiter::UpstreamBC::Velocity[1]+
      Orbiter::UpstreamBC::Velocity[2]*Orbiter::UpstreamBC::Velocity[2]);

  if (speed>0.0) {
    for (int idim=0;idim<3;idim++) dpProjection+=dp[idim]*Orbiter::UpstreamBC::Velocity[idim]/speed;
    Orbiter::Sampling::DragCoefficient::dpProjection2[thread]+=dpProjection*dpProjection;
  }

  //the rate of the kinetic energy gain by the gas
  de=(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2] -
      vInitUnchanged[0]*vInitUnchanged[0]-vInitUnchanged[1]*vInitUnchanged[1]-vInitUnchanged[2]*vInitUnchanged[2])*RateFactor/2.0;
  AddToFaceData(TriangleCutFace->UserData.EnergyTransferRate+spec,de/TriangleCutFace->SurfaceArea);

  return code;
}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  double res,size;

  size=t->CharacteristicSize();

  if (size<0.01) res=0.01;
  else if (size<1.0) res=0.01*pow(10.0,size);
  else res=0.25;

  //reduce the mesh resolution when run tests
  if (_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_) { 
    if (_ORBITER__NIGHTLY_TEST_REDUCE_RESOLUTION_MODE_==_PIC_MODE_ON_) {  
      res*=_PIC_NIGHTLY_TEST__GRID_RESOLUTION_MULTIPLIER_;
    }

    //limit the grid resolution with the floor value 
    if (res<_PIC_NIGHTLY_TEST__GRID_RESOLUTION_FLOOR_VALUE_) res=_PIC_NIGHTLY_TEST__GRID_RESOLUTION_FLOOR_VALUE_; 
  }

  return res;
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;

  //Orbiter::CharacteristicSpeed is 50 km/s unless the runtime input sets it
  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/Orbiter::CharacteristicSpeed;
}



double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}




/*  double Orbiter::GetTotalProduction(int spec,void *BoundaryElement) {
    return 1.0E20;
  }

  double Orbiter::GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement) {
    return GetTotalProduction(spec,BoundaryElement);
  }*/



//the command line: "amps" or, for the drag-coefficient test, "amps -i FILE -case NAME" (OrbiterRuntimeInput.h)
static void PrintUsage(const char *program) {
  printf("usage: %s [-i|--input FILE [-case|--case NAME]]\n"
      "  -i FILE    the runtime input (src/INPUT_FILE.md; the drag-coefficient test: input/test/sphere-drag-coefficient.in)\n"
      "  -case NAME the case of FILE; required when FILE has cases\n",program);
}

int main(int argc,char **argv) {
  char fname[_MAX_STRING_LENGTH_PIC_];
  std::string RuntimeInputPath,RuntimeCaseName;
  bool CommandLineError=false,HelpRequested=false;

  //init the particle solver
  PIC::InitMPI();

  //the command line
  for (int i=1;i<argc;i++) {
    std::string option=argv[i];

    if ((option=="-i")||(option=="--input")||(option=="-input")) {
      if ((++i>=argc)||(RuntimeInputPath.empty()==false)) {
        if (PIC::ThisThread==0) fprintf(stderr,"%s requires exactly one file\n",option.c_str());
        CommandLineError=true;
      }
      else RuntimeInputPath=argv[i];
    }
    else if ((option=="-case")||(option=="--case")) {
      if ((++i>=argc)||(RuntimeCaseName.empty()==false)) {
        if (PIC::ThisThread==0) fprintf(stderr,"%s requires exactly one case name\n",option.c_str());
        CommandLineError=true;
      }
      else RuntimeCaseName=argv[i];
    }
    else if ((option=="-h")||(option=="--help")) HelpRequested=true;
    else {
      if (PIC::ThisThread==0) fprintf(stderr,"unknown argument '%s'\n",option.c_str());
      CommandLineError=true;
    }
  }

  if ((CommandLineError==false)&&(RuntimeInputPath.empty())&&(RuntimeCaseName.empty()==false)) {
    if (PIC::ThisThread==0) fprintf(stderr,"-case requires -i FILE\n");
    CommandLineError=true;
  }

  if ((CommandLineError==false)&&(HelpRequested==false)) {
    if ((_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_)&&(_ORBITER_EXECUTE_NEW_TEST_==_PIC_MODE_ON_)) {
      if (RuntimeInputPath.empty()) {
        if (PIC::ThisThread==0) fprintf(stderr,"the drag-coefficient test requires its runtime input: -i FILE -case NAME\n");
        CommandLineError=true;
      }
    }
    else if (RuntimeInputPath.empty()==false) {
      if (PIC::ThisThread==0) fprintf(stderr,"-i: this executable takes no runtime input (only the drag-coefficient test does)\n");
      CommandLineError=true;
    }
  }

  if (CommandLineError||HelpRequested) {
    if (PIC::ThisThread==0) PrintUsage(argv[0]);
    MPI_Finalize();
    return HelpRequested&&(CommandLineError==false) ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  //seed the random number generatore
  if (_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_) {
    rnd_seed(100);
  }


  PIC::Init_BeforeParser();
  PIC::Alarm::SetAlarm(48*3600-30*60);

  Orbiter::Init_BeforeParser();

  //the runtime input: it is read and checked before the surface and the mesh are read or built
  if (RuntimeInputPath.empty()==false) {
    std::string error;

    if (Orbiter::RuntimeInput::LoadAndApply(RuntimeInputPath,RuntimeCaseName,&error)==false) {
      if (PIC::ThisThread==0) fprintf(stderr,"invalid runtime input: %s\n",error.c_str());
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    if (PIC::ThisThread==0) {
      printf("$PREFIX: runtime input=%s case=%s\n",RuntimeInputPath.c_str(),RuntimeCaseName.c_str());
      printf("%s",Orbiter::RuntimeInput::Describe().c_str());
      printf("%s",Surface::RuntimeInput::Describe().c_str());
      fflush(stdout);
    }
  }


  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};



  //load the NASTRAN mesh
  list<CutCell::cSurfaceMeshFile> SurfaceMeshFileList;
  CutCell::cSurfaceMeshFile MeshFile;

  //load the surface model
  for (int iSurfaceModelPart=0;iSurfaceModelPart<Orbiter::SurfaceModel::nTotalSurfaceModelFiles;iSurfaceModelPart++) {
    MeshFile.faceat=Orbiter::SurfaceModel::SurfaceModelSet[iSurfaceModelPart].faceat;
    sprintf(MeshFile.MeshFileName,"%s/%s",PIC::UserModelInputDataPath,Orbiter::SurfaceModel::SurfaceModelSet[iSurfaceModelPart].FileName);

    SurfaceMeshFileList.push_back(MeshFile);
  }

  switch (Orbiter::SurfaceModel::MeshFileFormat) {
  case Orbiter::SurfaceModel::MeshFileFormat_CEA:
    PIC::Mesh::IrregularSurface::ReadCEASurfaceMeshLongFormat(SurfaceMeshFileList,Orbiter::SurfaceModel::ScalingFactor);
    break;

  case Orbiter::SurfaceModel::MeshFileFormat_NASTRAN:
    PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat(SurfaceMeshFileList,Orbiter::SurfaceModel::ScalingFactor);
    break;

  default:
    exit(__LINE__,__FILE__,"Error: the format is unknown");
  }


  PIC::Mesh::IrregularSurface::GetSurfaceSizeLimits(xmin,xmax);
  PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh("SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);

  //calculate the projection are when the drag coefficient it to be calculated 
  if (Orbiter::Sampling::DragCoefficient::SamplingMode==true) {
    Orbiter::CalculateProjectionArea();
  }


//  PIC::Mesh::IrregularSurface::SmoothRefine(0.5);

  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }


  //set up the s/c object
  //init the nucleus
  cInternalBoundaryConditionsDescriptor OrbiterSurfaceDescriptor;
  cInternalNastranSurfaceData *Orbiter;

//  PIC::BC::InternalBoundary::RotationBody::Init(ReserveSamplingSpace,NULL);
  OrbiterSurfaceDescriptor=PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  Orbiter=(cInternalNastranSurfaceData*) OrbiterSurfaceDescriptor.BoundaryElement;

//  Orbiter->PrintTitle=Comet::Sampling::OutputSurfaceDataFile::PrintTitle;
//  Orbiter->PrintVariableList=Comet::Sampling::OutputSurfaceDataFile::PrintVariableList;

//  Nucleus->localResolution=Comet::localSphericalSurfaceResolution;
//  Orbiter->InjectionRate=Orbiter::GetTotalProduction;
  Orbiter->faceat=0;
//  Nucleus->ParticleSphereInteraction=Comet::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
//  Orbiter->InjectionBoundaryCondition=Orbiter::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;
//  Nucleus->InjectionBoundaryCondition=Comet::InjectionBoundaryModel_Limited; ///sphereParticleInjection;
  //PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;


  for (int idim=0;idim<3;idim++) {
    double OrbiterSize=xmax[idim]-xmin[idim];

    xmin[idim]-=OrbiterSize*Orbiter::DomainSize::xMinOffset[idim];
    xmax[idim]+=OrbiterSize*Orbiter::DomainSize::xMaxOffset[idim];
  }

  PIC::Mesh::mesh->CutCellSurfaceLocalResolution=SurfaceResolution;
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,BulletLocalResolution);

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();


  //PIC::Mesh::mesh->buildMesh();

  //generate mesh or read from file
  char mesh[400];
  bool NewMeshGeneratedFlag=false;

  FILE *fmesh=NULL;

  sprintf(mesh,"%s/amr.sig=%s.mesh.bin",PIC::UserModelInputDataPath,Orbiter::Mesh::sign);
  fmesh=fopen(mesh,"r");

  if (fmesh!=NULL) {
    if (PIC::ThisThread==0) printf("$PREFIX: mesh %s is found: loading.... ",mesh);

    fclose(fmesh);
    PIC::Mesh::mesh->readMeshFile(mesh);

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread==0) printf("done. \n");
  }
  else {
    NewMeshGeneratedFlag=true;
    if (PIC::ThisThread==0) {
      printf("$PREFIX: mesh %s is not found: generating.... ",mesh);
      fflush(stdout);
    }

    if (PIC::Mesh::mesh->ThisThread==0) {
      PIC::Mesh::mesh->buildMesh();
      PIC::Mesh::mesh->saveMeshFile("mesh.msh");
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Mesh::mesh->readMeshFile("mesh.msh");
    }

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      printf("  done.\n");
      fflush(stdout);
    }
  }


  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();



  //initialize the blocks
  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();


  //init the external normals of the cut faces
  double xLightSource[3]={6000.0e3,0,0.0};

  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }

  PIC::Mesh::IrregularSurface::InitExternalNormalVector();

  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }

  //test the shadow procedure
  PIC::Mesh::IrregularSurface::CutFaceAccessCounter::Init();
  PIC::RayTracing::SetCutCellShadowAttribute(xLightSource,true);

  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation-shadow.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }

  //output the volume mesh
  sprintf(fname,"%s/VolumeMesh.dat",PIC::OutputDataFileDirectory);
  PIC::Mesh::mesh->outputMeshTECPLOT(fname);

  //init the volume of the cells'
  PIC::Mesh::mesh->InitCellMeasure();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh->getMeshSignature();

    if (PIC::Mesh::mesh->ThisThread==0) {
      char command[300];

      printf("$PREFIX: Renaming the mesh file to amr.sig=0x%lx.mesh.bin.....  ",MeshSignature);
      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
      printf(" done.\n");
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);



//  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=2000; //0; //00; //*10;
//  PIC::RequiredSampleLength=500; //00; //0; //0;


  PIC::Init_AfterParser ();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();



  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=Orbiter::UpstreamBC::BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=Orbiter::UpstreamBC::BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();


  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=Orbiter::UpstreamBC::BoundingBoxInjectionRate;
  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);

  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::AdjustParticleWeight_ConstantWeightOverTimeStep_KeepMinParticleWeight(s);

  //init the particle buffer
//  PIC::ParticleBuffer::Init(1000000);

  //set the model of the boundary conditinos
  PIC::Mover::ProcessTriangleCutFaceIntersection=SurfaceBoundaryCondition;


  //determine the total number of the iterations to perform
  //in the test-mode run 100 iterations and than output the particle data statistics
  int nIterations,nTotalIterations=100000001;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {

    if (_ORBITER_EXECUTE_NEW_TEST_==_PIC_MODE_ON_) {
      DragCoefficientTest::Execute();
    }

    nTotalIterations=_PIC_NIGHTLY_TEST__TOTAL_ITERATION_NUMBER_;   
  } 

  //time step
  for (long int niter=0;niter<nTotalIterations;niter++) {
    static int LastDataOutputFileNumber=-1;

    PIC::TimeStep();

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
      PIC::RequiredSampleLength*=2;
      if (PIC::RequiredSampleLength>1000) PIC::RequiredSampleLength=1000;


      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
      if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
    }

    if (PIC::Mesh::mesh->ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (current sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }
  }

  //output the particle statistics for the nightly tests
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];

    sprintf(fname,"%s/test_Orbiter.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }

  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}

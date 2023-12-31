//$Id$


#include "pic.h"
#include "constants.h"

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

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  double res,size;

  size=t->CharacteristicSize();

  if (size<0.01) res=0.01;
  else if (size<1.0) res=0.01*pow(10.0,size);
  else res=0.25;

  //reduce the mesh resolution when run tests
  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  if (res<0.20) res=0.20;
  #endif

  return res;
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;
  double CharacteristicSpeed=1.0E3;

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}

double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  bool ExternalFaces[6];
  double res=0.0,ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  static double v[3]={2.0e3,000.0,000.0},n=5.0E6,temp=20.0;


  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      if (spec!=_H2O_SPEC_) return 0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n,temp,v,ExternalNormal,_H2O_SPEC_);

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }

  return res;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={2.0e3,0.0,0.0},nNA=5.0E6,tempNA=20.0;

  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_H2O_SPEC_);

      if (ModelParticlesInjectionRate>0.0) return true;
    }
  }

  return false;
}

//injection of model particles through the faces of the bounding box
long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  long int nInjectedParticles=0;

  if (spec!=_H2O_SPEC_) return 0; //inject only spec=0

  static double vNA[3]={2.0e3,000.0,000.0},nNA=5.0E6,tempNA=20.0;
  double v[3];


  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_H2O_SPEC_);


      if (ModelParticlesInjectionRate>0.0) {
        ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;

        PIC::Mesh::mesh->GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

        while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
          //generate the new particle position on the face
          for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

          //generate a particle
          newParticle=PIC::ParticleBuffer::GetNewParticle();
          newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
          nInjectedParticles++;

          //generate particles' velocity
          PIC::Distribution::InjectMaxwellianDistribution(v,vNA,tempNA,ExternalNormal,_H2O_SPEC_,-1);

          PIC::ParticleBuffer::SetX(x,newParticleData);
          PIC::ParticleBuffer::SetV(v,newParticleData);
          PIC::ParticleBuffer::SetI(spec,newParticleData);
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

          //inject the particle into the system
          _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);
        }
      }


    }
  }

  return nInjectedParticles;
}

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  long int nInjectedParticles=0;

  for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

  return nInjectedParticles;
}

double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}




  double Rosetta::GetTotalProduction(int spec,void *BoundaryElement) {
    return 1.0E20;
  }

  double Rosetta::GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement) {
    return GetTotalProduction(spec,BoundaryElement);
  }


  bool Rosetta::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
//  bool Rosetta::GenerateParticleProperties(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, char *tempParticleData,int BoundaryElementType,void *BoundaryElement) {

    static long int ProbabilityTableLengh=0;
    static double ProbabilityTableIncrement=0.0;

    class cFaceDescriptor {
    public:
      double weight;
      int nFace;
      cFaceDescriptor* next;

      cFaceDescriptor() {
        weight=0.0,nFace=-1,next=NULL;
      }
    };

    class cProbabilityTable {
    public:
      int nTotalFaces;
      cFaceDescriptor *firstFaceDescriptor;

      cProbabilityTable() {
        nTotalFaces=0;
        firstFaceDescriptor=NULL;
      }
    };


    static cFaceDescriptor *FaceDescriptorTable=NULL;
    static cProbabilityTable *ProbabilityTable=NULL;

    static bool initflag=false;

    if (initflag==false) {
      initflag=true;

      double t,minSurfaceArea=-1.0,totalSurfaceArea=0.0;
      int nt,cnt,next;

      const int defaultProbabilityTableLengh=10000;

      //calculate the length of the probability table
      for (nt=0;nt<CutCell::nBoundaryTriangleFaces;nt++) {
        if ((minSurfaceArea<0.0)||(minSurfaceArea>CutCell::BoundaryTriangleFaces[nt].SurfaceArea)) minSurfaceArea=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
        totalSurfaceArea+=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
      }

      ProbabilityTableIncrement=totalSurfaceArea/defaultProbabilityTableLengh;
      ProbabilityTableLengh=defaultProbabilityTableLengh;

      ProbabilityTable=new cProbabilityTable [ProbabilityTableLengh];

      //calculate the number of the face descriptors that is needed for the mesh
      int nFaceDescriptor=0,iStart=0,iFinish=0;

      for (nt=0,t=0.0;nt<CutCell::nBoundaryTriangleFaces;nt++) {
        t+=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
        iFinish=(int)(t/ProbabilityTableIncrement);

        nFaceDescriptor+=iFinish-iStart+1;
        iStart=iFinish;
      }

      FaceDescriptorTable=new cFaceDescriptor [nFaceDescriptor];

      //init the face descriptor table
      for (cnt=0,nt=0,t=0.0,iStart=0;nt<CutCell::nBoundaryTriangleFaces;nt++) {
        t+=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
        iFinish=(int)(t/ProbabilityTableIncrement);

        if (iFinish>=ProbabilityTableLengh) iFinish=ProbabilityTableLengh-1;

        for (int ii=iStart;ii<=iFinish;ii++) {
          if (iStart==iFinish) FaceDescriptorTable[cnt].weight=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
          else {
        	  if (ii==iStart) FaceDescriptorTable[cnt].weight=(iStart+1)*ProbabilityTableIncrement-(t-CutCell::BoundaryTriangleFaces[nt].SurfaceArea);
        	  else if (ii==iFinish) FaceDescriptorTable[cnt].weight=t-iFinish*ProbabilityTableIncrement;
        	  else FaceDescriptorTable[cnt].weight=ProbabilityTableIncrement;
          }

          if (cnt<nFaceDescriptor) {
            FaceDescriptorTable[cnt].nFace=nt;
            FaceDescriptorTable[cnt].next=ProbabilityTable[ii].firstFaceDescriptor;

            ProbabilityTable[ii].firstFaceDescriptor=FaceDescriptorTable+cnt;
            ProbabilityTable[ii].nTotalFaces++;
          }

          cnt++;
        }

        iStart=iFinish;
      }

      //normalize weights
      for (int np=0;np<ProbabilityTableLengh;np++) if (ProbabilityTable[np].firstFaceDescriptor!=NULL) {
        double summ=0.0;
        cFaceDescriptor *face,*prev;

        for (face=ProbabilityTable[np].firstFaceDescriptor;face!=NULL;face=face->next) summ+=face->weight;
        for (face=ProbabilityTable[np].firstFaceDescriptor;face!=NULL;face=face->next) face->weight/=summ;

        //Convert the weight into a cumulative distribution
        for (summ=0.0,face=ProbabilityTable[np].firstFaceDescriptor;face!=NULL;face=face->next) {
          summ+=face->weight;
          face->weight=summ;
        }
      }
    }


    //Determine the face number
    int nt,nface=-1;
    double x[3],v[3];
    double weight=rnd();
    cFaceDescriptor *face;

    nt=(int)(rnd()*ProbabilityTableLengh);

    for (face=ProbabilityTable[nt].firstFaceDescriptor;face!=NULL;face=face->next) if (face->weight>=weight){
  	  nface=face->nFace;
	    break;
    }

    bool PositionGenerated;

    do {
    	PositionGenerated=true;

		  CutCell::BoundaryTriangleFaces[nface].GetRandomPosition(x);

	  	//place the point inside the domain
		  for (int idim=0;idim<3;idim++) x[idim]+=0.001*PIC::Mesh::mesh->EPS*CutCell::BoundaryTriangleFaces[nface].ExternalNormal[idim];
    }
    while (PositionGenerated==false);


    //determine if the particle belongs to this processor
    startNode=PIC::Mesh::mesh->findTreeNode(x,startNode);
    if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) return false;

    //get the velocity vector
//    for (int idim=0;idim<3;idim++) v[idim]=500.0*CutCell::BoundaryTriangleFaces[nface].ExternalNormal[idim];

    const static double vbulk[3]={0.0,0.0,0.0};
    const static double SurfaceTemperature=200.0;

    PIC::Distribution::InjectMaxwellianDistribution(v,vbulk,SurfaceTemperature,CutCell::BoundaryTriangleFaces[nface].ExternalNormal,_H2O_SPEC_);



    memcpy(x_SO_OBJECT,x,3*sizeof(double));
    memcpy(x_IAU_OBJECT,x,3*sizeof(double));
    memcpy(v_SO_OBJECT,v,3*sizeof(double));
    memcpy(v_IAU_OBJECT,v,3*sizeof(double));

    return true;

  }



int main(int argc,char **argv) {
  char fname[_MAX_STRING_LENGTH_PIC_];

  //init the particle solver
  PIC::InitMPI();
  PIC::Init_BeforeParser();

  Rosetta::Init_BeforeParser();


  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};



  //load the NASTRAN mesh
  sprintf(fname,"%s/rosetta.surface.reduced.nas",PIC::UserModelInputDataPath);
  PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat(fname); //("C-G_MOC_original.bdf"); //("rosetta.surface.reduced.nas");
  PIC::Mesh::IrregularSurface::GetSurfaceSizeLimits(xmin,xmax);
  PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh("SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);


//  PIC::Mesh::IrregularSurface::SmoothRefine(0.5);

  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }


  //set up the s/c object
  //init the nucleus
  cInternalBoundaryConditionsDescriptor RosettaSurfaceDescriptor;
  cInternalNastranSurfaceData *Rosetta;

//  PIC::BC::InternalBoundary::RotationBody::Init(ReserveSamplingSpace,NULL);
  RosettaSurfaceDescriptor=PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  Rosetta=(cInternalNastranSurfaceData*) RosettaSurfaceDescriptor.BoundaryElement;

//  Rosetta->PrintTitle=Comet::Sampling::OutputSurfaceDataFile::PrintTitle;
//  Rosetta->PrintVariableList=Comet::Sampling::OutputSurfaceDataFile::PrintVariableList;

//  Nucleus->localResolution=Comet::localSphericalSurfaceResolution;
  Rosetta->InjectionRate=Rosetta::GetTotalProduction;
  Rosetta->faceat=0;
//  Nucleus->ParticleSphereInteraction=Comet::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
  Rosetta->InjectionBoundaryCondition=Rosetta::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;
//  Nucleus->InjectionBoundaryCondition=Comet::InjectionBoundaryModel_Limited; ///sphereParticleInjection;
  //PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;




  for (int i=0;i<3;i++) xmin[i]*=2.0,xmax[i]*=2.0;


  PIC::Mesh::mesh->CutCellSurfaceLocalResolution=SurfaceResolution;
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh->memoryAllocationReport();


  //PIC::Mesh::mesh->buildMesh();

  //generate mesh or read from file
  char mesh[200]="amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  bool NewMeshGeneratedFlag=false;

  FILE *fmesh=NULL;

  fmesh=fopen(mesh,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh->readMeshFile(mesh);
  }
  else {
    NewMeshGeneratedFlag=true;

    if (PIC::Mesh::mesh->ThisThread==0) {
       PIC::Mesh::mesh->buildMesh();
       PIC::Mesh::mesh->saveMeshFile("mesh.msh");
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
       PIC::Mesh::mesh->readMeshFile("mesh.msh");
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

  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh->getMeshSignature();

    if (PIC::Mesh::mesh->ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);



  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=2000; //0; //00; //*10;
  PIC::RequiredSampleLength=500; //00; //0; //0;


  PIC::Init_AfterParser ();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=localParticleInjectionRate;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H2O_SPEC_);

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();

  //init the particle buffer
  PIC::ParticleBuffer::Init(1000000);

  //set the model of the boundary conditinos
  PIC::Mover::ProcessTriangleCutFaceIntersection=SurfaceBoundaryCondition;


  //determine the total number of the iterations to perform
  //in the test-mode run 100 iterations and than output the particle data statistics
  int nIterations,nTotalIterations=100000001;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=550;

  //time step
  for (long int niter=0;niter<nTotalIterations;niter++) {
    static int LastDataOutputFileNumber=-1;

    PIC::TimeStep();

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
      PIC::RequiredSampleLength*=2;
      if (PIC::RequiredSampleLength>40000) PIC::RequiredSampleLength=40000;


      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
      if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
    }

  }

  //output the particle statistics for the nightly tests
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];

    sprintf(fname,"%s/test_Rosetta.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }


  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}

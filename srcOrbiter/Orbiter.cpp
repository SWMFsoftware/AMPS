/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"
#include "global.h"

//the object name and the names of the frames
/*char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Orbiter";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MOON";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="LSO";*/

double Orbiter::UpstreamBC::Velocity[3]={-30.0E3,0.0,0.0};
double Orbiter::UpstreamBC::NumberDensity[PIC::nTotalSpecies];
double Orbiter::UpstreamBC::Temperature=293.0;
bool Orbiter::UpstreamBC::UpstreamSourceMode=true;

//scaling factor of the surface model
double Orbiter::SurfaceModel::ScalingFactor=1.0;

bool Orbiter::Sampling::DragCoefficient::SamplingMode=false;
double Orbiter::ProjectionOrbiterSurfaceArea=0.0;

Orbiter::SurfaceModel::cSurfaceModelSet Orbiter::SurfaceModel::SurfaceModelSet[]={{0,"Orbiter"}};
int Orbiter::SurfaceModel::MeshFileFormat=Orbiter::SurfaceModel::MeshFileFormat_NASTRAN;

double Orbiter::DomainSize::xMinOffset[3]={0.1,0.1,0.1};
double Orbiter::DomainSize::xMaxOffset[3]={0.1,0.1,0.1};

char Orbiter::Mesh::sign[_MAX_STRING_LENGTH_PIC_]="";

double *Orbiter::Sampling::DragCoefficient::dpX=NULL;
double *Orbiter::Sampling::DragCoefficient::dpY=NULL;
double *Orbiter::Sampling::DragCoefficient::dpZ=NULL;
double *Orbiter::Sampling::DragCoefficient::Flux=NULL;
int *Orbiter::Sampling::DragCoefficient::ModelParticleCounter=NULL;

Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::cSampledData *Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::SampledData=NULL;
double Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::EndFlowDirection[3]={0.0,0.0,0.0};
double Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::StartFlowDirection[3]={0.0,0.0,0.0};
bool Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::ResetUpstreamFlowMode=false;
int Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::ResetOutputNumberInterval=1;
double Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::BulkFlowSpeed=0.0;
int Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::nTotalDirectionResets=0;


void Orbiter::Sampling::DragCoefficient::Init() {
  dpX=new double [PIC::nTotalThreadsOpenMP];
  dpY=new double [PIC::nTotalThreadsOpenMP];
  dpZ=new double [PIC::nTotalThreadsOpenMP];
  Flux=new double [PIC::nTotalThreadsOpenMP];
  ModelParticleCounter=new int [PIC::nTotalThreadsOpenMP];

  for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
    dpX[thread]=0.0;
    dpY[thread]=0.0;
    dpZ[thread]=0.0;
    Flux[thread]=0.0;
    ModelParticleCounter[thread]=0;
  }

  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(SamplingProcessor,PrintOutputFile);

  //init the module that veried the direction of the upstream flow
  if (ResetUpstreamFlowDirection::ResetUpstreamFlowMode==true) {
    int i;

    memcpy(ResetUpstreamFlowDirection::StartFlowDirection,Orbiter::UpstreamBC::Velocity,3*sizeof(double));

    ResetUpstreamFlowDirection::BulkFlowSpeed=Vector3D::Length(Orbiter::UpstreamBC::Velocity);
    Vector3D::Normalize(ResetUpstreamFlowDirection::StartFlowDirection);
    Vector3D::Normalize(ResetUpstreamFlowDirection::EndFlowDirection);

    ResetUpstreamFlowDirection::SampledData=new ResetUpstreamFlowDirection::cSampledData[ResetUpstreamFlowDirection::nTotalDirectionResets+1];

    for (i=0;i<ResetUpstreamFlowDirection::nTotalDirectionResets+1;i++) {
      ResetUpstreamFlowDirection::SampledData[i].Cd=-1.0;
      ResetUpstreamFlowDirection::SampledData[i].AngleStartingDirection=0.0;
    }
  }
}

//=======================================================================================
//remove all model particles to restart particle sampling for a new boundary condition case
void Orbiter::Sampling::DragCoefficient::ResetUpstreamFlowDirection::RemoveAllParticles() {
  int Ptr,nextPrt;

  //reset particle lists in all cells
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  int i,j,k;

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

    if ((block=node->block)!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) if ((Ptr=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)])!=-1) {
      block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;

      do {
        nextPrt=PIC::ParticleBuffer::GetNext(Ptr);
        PIC::ParticleBuffer::DeleteParticle(Ptr);

        Ptr=nextPrt;
      }
      while (Ptr!=-1);
    }

  }
}

//=======================================================================================
//particle/surface interaction model
int Orbiter::ParticleSurfaceInteractionProcessor_default(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double c;

  c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];
  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//=======================================================================================
//sampling functions that are used in calculating of the drag coefficient 
void Orbiter::Sampling::DragCoefficient::SamplingProcessor() {}

void Orbiter::Sampling::DragCoefficient::PrintOutputFile(int nfile) {
  double Cd,dpTotal[3],FluxTotal;
  int thread,ModelParticleCounterTotal;

  static int iUpstreamFlowDirectionResetCycle=0;

  for (thread=1;thread<PIC::nTotalThreadsOpenMP;thread++) {
    dpX[0]+=dpX[thread];
    dpY[0]+=dpY[thread];
    dpZ[0]+=dpZ[thread];
    Flux[0]+=Flux[thread];
    ModelParticleCounter[0]+=ModelParticleCounter[thread];
  }

  //collect the momentum exchange rate from all processors
  MPI_Reduce(dpX,dpTotal+0,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(dpY,dpTotal+1,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(dpZ,dpTotal+2,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR); 
  MPI_Reduce(Flux,&FluxTotal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(ModelParticleCounter,&ModelParticleCounterTotal,1,MPI_INT,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    //calculate the drag coefficient 
    double speed=0.0,dpProjection=0.0,TotalMassDensity=0.0,TotalNumberDensity=0.0,TotalMassFlux=0.0,TotalNumberDensityFlux=0.0;
    int spec,idim;

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      TotalMassDensity+=Orbiter::UpstreamBC::NumberDensity[spec]*PIC::MolecularData::GetMass(spec);
      TotalMassFlux+=Orbiter::UpstreamBC::NumberDensity[spec]*PIC::MolecularData::GetMass(spec)*Vector3D::Length(Orbiter::UpstreamBC::Velocity);

      TotalNumberDensity+=Orbiter::UpstreamBC::NumberDensity[spec];
      TotalNumberDensityFlux+=Orbiter::UpstreamBC::NumberDensity[spec]*Vector3D::Length(Orbiter::UpstreamBC::Velocity);
    }

    for (idim=0;idim<3;idim++) speed+=pow(Orbiter::UpstreamBC::Velocity[idim],2);
    speed=sqrt(speed);
    for (idim=0;idim<3;idim++) dpProjection+=Orbiter::UpstreamBC::Velocity[idim]*dpTotal[idim]/(speed*PIC::LastSampleLength);

    Cd=2.0*dpProjection/(TotalMassDensity*ProjectionOrbiterSurfaceArea*speed*speed);
    FluxTotal/=PIC::LastSampleLength;

    printf("$PREFIX: Calculated Drag Coefficient:\n");
    printf("$PREFIX: Mass Density=%e; Mass Density Flux=%e\n",TotalMassDensity,TotalMassFlux);
    printf("$PREFIX: Number Density=%e; Number Density Flux=%e\n",TotalNumberDensity,TotalNumberDensityFlux);
    printf("$PREFIX: dp/dt=%e, %e, %e\n",dpTotal[0],dpTotal[1],dpTotal[2]);
    printf("$PREFIX: dpProjection=%e, Speed=%e\n",dpProjection,speed);
    printf("$PREFIX: ProjectionOrbiterSurfaceArea=%e\n",ProjectionOrbiterSurfaceArea);
    printf("$PREFIX: Calculated Drag Coefficient = %e\n",Cd);

    printf("$PREFIX: Total flux of the particle interesected with the surface = %e\n",FluxTotal);
    printf("$PREFIX: Total flux of the particle interesected with the surface / ProjectionOrbiterSurfaceArea = %e\n",FluxTotal/ProjectionOrbiterSurfaceArea);

    printf("$PREFIX: Total number of the model particle that have intersected the surface = %i\n",ModelParticleCounterTotal);
    printf("$PREFIX: Total number of the model particle that have intersected the surface / ProjectionOrbiterSurfaceArea= %e\n",ModelParticleCounterTotal/ProjectionOrbiterSurfaceArea);

    
    //output sampeled drag coefficients calculated for different directions of the upstream flow
    if (ResetUpstreamFlowDirection::ResetUpstreamFlowMode==true) {
      if (iUpstreamFlowDirectionResetCycle<=ResetUpstreamFlowDirection::nTotalDirectionResets) ResetUpstreamFlowDirection::SampledData[iUpstreamFlowDirectionResetCycle].Cd=Cd;
      printf("Drag Coefficient Sampled for Multiple Upstream Directions\ni\tAngle\t\tCd\n");

      for (int i=0;i<ResetUpstreamFlowDirection::nTotalDirectionResets+1;i++) {
        printf("%i\t%e\t%e\n",i,ResetUpstreamFlowDirection::SampledData[i].AngleStartingDirection,ResetUpstreamFlowDirection::SampledData[i].Cd);
      }

      printf("\n");
    }

    //output the particle statistics for the nightly tests
    if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
      char fname[400];
      FILE *fout;

      sprintf(fname,"%s/test_Orbiter_Drag_Coefficient.dat",PIC::OutputDataFileDirectory);
      fout=fopen(fname,"w");

      fprintf(fout,"Mass Density:\n%e\ndp/dt:\n%e %e %e\n",TotalMassDensity,dpTotal[0],dpTotal[1],dpTotal[2]);
      fprintf(fout,"dpProjection\n%e\nSpeed\n%e\n",dpProjection,speed);
      fprintf(fout,"ProjectionOrbiterSurfaceArea\n%e\n",ProjectionOrbiterSurfaceArea);
      fprintf(fout,"Calculated Drag Coefficient\n%e\n",Cd);

      fclose(fout);
    }

  }

  //reset the sampling buffer 
  for (thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
    dpX[thread]=0.0;
    dpY[thread]=0.0;
    dpZ[thread]=0.0;
    Flux[thread]=0.0;
    ModelParticleCounter[thread]=0;
  }  

  //reset the direction of the unstrean flow if needed
  if ((ResetUpstreamFlowDirection::ResetUpstreamFlowMode==true)&&(nfile!=0)&&(nfile%ResetUpstreamFlowDirection::ResetOutputNumberInterval==0)) {
    iUpstreamFlowDirectionResetCycle=nfile/ResetUpstreamFlowDirection::ResetOutputNumberInterval;

     if (iUpstreamFlowDirectionResetCycle<=ResetUpstreamFlowDirection::nTotalDirectionResets) {
       //save the previously lates drag coefficient value
//       ResetUpstreamFlowDirection::SampledData[iResetCycle-1].Cd=Cd;

       //remove all particles and then reset the direcion of the upstream flow
       ResetUpstreamFlowDirection::RemoveAllParticles();

       //the new direction of the upstream flow
       double l[3],c;
       int idim;

       c=double(iUpstreamFlowDirectionResetCycle)/double(ResetUpstreamFlowDirection::nTotalDirectionResets);

       for (idim=0;idim<3;idim++) l[idim]=(1.0-c)*ResetUpstreamFlowDirection::StartFlowDirection[idim]+c*ResetUpstreamFlowDirection::EndFlowDirection[idim];
       Vector3D::Normalize(l);
       for (idim=0;idim<3;idim++) Orbiter::UpstreamBC::Velocity[idim]=ResetUpstreamFlowDirection::BulkFlowSpeed*l[idim];

       ResetUpstreamFlowDirection::SampledData[iUpstreamFlowDirectionResetCycle].AngleStartingDirection=180.0/Pi*acos(Vector3D::DotProduct(l,ResetUpstreamFlowDirection::StartFlowDirection));

       //recalcualte the new projection  area
       Orbiter::CalculateProjectionArea();
     }
  }

} 
  

//=======================================================================================
//calculate the projection area of the spacecraft
double Orbiter::CalculateProjectionArea() {
  unsigned char *map=NULL;
  double *x,de0,de1,c0,c1,e0[3],e1[3],e2[3],e0Min=0.0,e0Max=0.0,e1Min=0.0,e1Max=0.0;
  int idim,iSurfaceElement,iNode;
  int i,j,iOffset,jOffset;

  int nPixels=2000;
  int ByteMapLength=1+(nPixels*nPixels)/8;

  //get the frame of reference with e0 and e1 normal to the ambient gas flow
  for (idim=0;idim<3;idim++) e2[idim]=UpstreamBC::Velocity[idim];

  Vector3D::Normalize(e2);
  Vector3D::GetNormFrame(e0,e1,e2);

  //get the limits of the projected coordinates
  for (iNode=0;iNode<CutCell::nBoundaryTriangleNodes;iNode++) {
    x=CutCell::BoundaryTriangleNodes[iNode].x;

    c0=Vector3D::DotProduct(e0,x);
    c1=Vector3D::DotProduct(e1,x);

    if (iNode==0) e0Min=c0,e0Max=c0;
    else {
      if (c0>e0Max) e0Max=c0;
      if (c0<e0Min) e0Min=c0;
    }

    if (iNode==0) e1Min=c1,e1Max=c1;
    else {
      if (c1>e1Max) e1Max=c1;
      if (c1<e1Min) e1Min=c1;
    }
  }


  de0=(e0Max-e0Min)/nPixels;
  de1=(e1Max-e1Min)/nPixels;

  //allocate the projection map
  map=new unsigned char[ByteMapLength];

  for (i=0;i<ByteMapLength;i++) map[i]=0;

  //loop though all surface elements
  int iSurfaceElementStart,iSurfaceElementStartEnd,nSurfaceElementThread;

  nSurfaceElementThread=CutCell::nBoundaryTriangleFaces/PIC::nTotalThreads;
  iSurfaceElementStart=nSurfaceElementThread*PIC::ThisThread;
  iSurfaceElementStartEnd=iSurfaceElementStart+nSurfaceElementThread;
  if (PIC::ThisThread==PIC::nTotalThreads-1) iSurfaceElementStartEnd=CutCell::nBoundaryTriangleFaces;


  for (iSurfaceElement=iSurfaceElementStart;iSurfaceElement<iSurfaceElementStartEnd;iSurfaceElement++) {
    int nTestPoints;
    int t,iLocal,jLocal;
    double xLocal[2],dxLocal;

    double e0TriangleMax,e0TriangleMin,e1TriangleMax,e1TriangleMin;

    c0=Vector3D::DotProduct(e0,CutCell::BoundaryTriangleFaces[iSurfaceElement].x0Face);
    c1=Vector3D::DotProduct(e1,CutCell::BoundaryTriangleFaces[iSurfaceElement].x0Face);

    e0TriangleMin=c0,e0TriangleMax=c0;
    e1TriangleMin=c1,e1TriangleMax=c1;

    for (i=1;i<3;i++) {
      switch (i) {
      case 1:
        c0=Vector3D::DotProduct(e0,CutCell::BoundaryTriangleFaces[iSurfaceElement].x1Face);
        c1=Vector3D::DotProduct(e1,CutCell::BoundaryTriangleFaces[iSurfaceElement].x1Face);
        break;
      case 2:
        c0=Vector3D::DotProduct(e0,CutCell::BoundaryTriangleFaces[iSurfaceElement].x2Face);
        c1=Vector3D::DotProduct(e1,CutCell::BoundaryTriangleFaces[iSurfaceElement].x2Face);
        break;
      default:
        exit(__LINE__,__FILE__,"Error: something is wrong with counting");
      }

      if (c0<e0TriangleMin) e0TriangleMin=c0;
      if (c0>e0TriangleMax) e0TriangleMax=c0;

      if (c1<e1TriangleMin) e1TriangleMin=c1;
      if (c1>e1TriangleMax) e1TriangleMax=c1;
    }

    nTestPoints=max(e0TriangleMax-e0TriangleMin,e1TriangleMax-e1TriangleMin)/min(de0,de1)*4;
    dxLocal=1.0/nTestPoints;

    for (iLocal=0;iLocal<nTestPoints;iLocal++) {
      xLocal[0]=iLocal*dxLocal;

      if (xLocal[0]<1.0) {
        for (jLocal=0;jLocal<nTestPoints;jLocal++) {
          xLocal[1]=jLocal*dxLocal;

          if ((xLocal[1]<1.0)&&(xLocal[0]+xLocal[1]<1.0)) {
            double xx[3];

            //the point belong to the surface triangle -> register it
            for (idim=0;idim<3;idim++) xx[idim]=CutCell::BoundaryTriangleFaces[iSurfaceElement].x0Face[idim]+
                xLocal[0]*CutCell::BoundaryTriangleFaces[iSurfaceElement].e0[idim]+xLocal[1]*CutCell::BoundaryTriangleFaces[iSurfaceElement].e1[idim];

            c0=Vector3D::DotProduct(e0,xx)-e0Min;
            c1=Vector3D::DotProduct(e1,xx)-e1Min;

            int iPixel,jPixel,iGlobalPixelIndex,iByte,iBit;

            iPixel=c0/de0;
            jPixel=c1/de1;

            iGlobalPixelIndex=iPixel+jPixel*nPixels;

            iByte=iGlobalPixelIndex/8;
            iBit=iGlobalPixelIndex%8;

            //mark the shadow bit
            if ((iByte<ByteMapLength)&&(iBit<8)) {
              map[iByte]|=(1<<iBit);
            }
          }
        }
      }
    }
  }

  //combine maps from all processors, and calculate the total projection area
  int nShadowPixel=0;

  if (PIC::ThisThread==0) {
    MPI_Status status;
    int thread;
    unsigned char *buffer=new unsigned char[ByteMapLength];

    for (thread=1;thread<PIC::nTotalThreads;thread++) {
      MPI_Recv(buffer,ByteMapLength,MPI_UNSIGNED_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

      for (i=0;i<ByteMapLength;i++) map[i]|=buffer[i];
    }

    delete [] buffer;

    for (i=0;i<ByteMapLength;i++) for (j=0;j<8;j++) if ((map[i]&(1<<j)) != 0) nShadowPixel++;
  }
  else {
    MPI_Send(map,ByteMapLength,MPI_UNSIGNED_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);
  }

  MPI_Bcast(&nShadowPixel,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

  ProjectionOrbiterSurfaceArea=(e0Max-e0Min)*(e1Max-e1Min)*double(nShadowPixel)/double(nPixels*nPixels);
  if (PIC::ThisThread==0) printf("$PREFIX: Orbiter Projection Surface Area: %e\n",ProjectionOrbiterSurfaceArea);

  //prepare the projection map of the orbiter
  if (PIC::ThisThread==0) {
    FILE *fout;
    char fname[100];

    sprintf(fname,"%s/OrbiterProjectionArea.dat",PIC::OutputDataFileDirectory);
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLE=\"X\",\"Y\",\"Surface Projection Flag\"\n");
    fprintf(fout,"ZONE I=%i, J=%i, DATAPACKING=POINT\n",nPixels,nPixels);

    for (j=0;j<nPixels;j++) for (i=0;i<nPixels;i++) {
      int flag,iGlobalPixelIndex,iByte,iBit;

      iGlobalPixelIndex=i+j*nPixels;

      iByte=iGlobalPixelIndex/8;
      iBit=iGlobalPixelIndex%8;

      if ((iByte<ByteMapLength)&&(iBit<8)) {
        flag=((map[iByte]&(1<<iBit))!=0) ? 1 : 0;
      }
      else {
        flag=0;
      }

      fprintf(fout,"%e %e %i\n",e0Min+(i+0.5)*de0,e1Min+(j+0.5)*de1,flag);
    }

    fclose(fout);
  }

  delete [] map;
  return ProjectionOrbiterSurfaceArea;
}

//====================================================================
//init the Orbiter model
void Orbiter::Init_BeforeParser() {
  //init the drag coefficient sampling 
  if (Orbiter::Sampling::DragCoefficient::SamplingMode==true) Orbiter::Sampling::DragCoefficient::Init();

  //set the injection boundary processor
  PIC::BC::UserDefinedParticleInjectionFunction=InjectionModel::InjectParticles;
  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate=InjectionModel::GetTotalInjectionRate;

  //init the exosphere model
  Exosphere::Init_BeforeParser();
}  

//====================================================================
//exchange model data
void Orbiter::ExchangeModelData() {
  int iSurfaceElement,spec;


  //prepare and exchange the adsotprion/desorption fluxes, and update the surface aboundance
  double *AdsorptionDesorptionExchangeBuffer=new double [CutCell::nBoundaryTriangleFaces];
  double *sumAdsorptionDesorptionExchangeBuffer=new double [CutCell::nBoundaryTriangleFaces];

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //exchange the adsorption flux
    for (iSurfaceElement=0;iSurfaceElement<CutCell::nBoundaryTriangleFaces;iSurfaceElement++) {
      AdsorptionDesorptionExchangeBuffer[iSurfaceElement]=CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.AdsorptionFlux[spec];
    }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    PIC::Debugger::CatchOutLimitValue(AdsorptionDesorptionExchangeBuffer,CutCell::nBoundaryTriangleFaces,__LINE__,__FILE__);
#endif
#endif

    MPI_Allreduce(AdsorptionDesorptionExchangeBuffer,sumAdsorptionDesorptionExchangeBuffer,CutCell::nBoundaryTriangleFaces,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);


    //update the adsorption flux, exchange the desorption flux
    for (iSurfaceElement=0;iSurfaceElement<CutCell::nBoundaryTriangleFaces;iSurfaceElement++) {
      CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.AdsorptionFlux[spec]=sumAdsorptionDesorptionExchangeBuffer[iSurfaceElement];

      AdsorptionDesorptionExchangeBuffer[iSurfaceElement]=CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.DesorptionFlux[spec];
    }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    PIC::Debugger::CatchOutLimitValue(AdsorptionDesorptionExchangeBuffer,CutCell::nBoundaryTriangleFaces,__LINE__,__FILE__);
#endif
#endif

    MPI_Allreduce(AdsorptionDesorptionExchangeBuffer,sumAdsorptionDesorptionExchangeBuffer,CutCell::nBoundaryTriangleFaces,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    //update the desorption flux, then update the surface aboundance, and then set the default vaelues to the Adsorption/Desorption buffers
    for (iSurfaceElement=0;iSurfaceElement<CutCell::nBoundaryTriangleFaces;iSurfaceElement++) {
      CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.DesorptionFlux[spec]=sumAdsorptionDesorptionExchangeBuffer[iSurfaceElement];

      CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.SpeciesSurfaceAboundance[spec]+=
          CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.AdsorptionFlux[spec]-CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.DesorptionFlux[spec];

      CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.AdsorptionFlux[spec]=0.0;
      CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.DesorptionFlux[spec]=0.0;
    }
  }

  delete [] AdsorptionDesorptionExchangeBuffer;
  delete [] sumAdsorptionDesorptionExchangeBuffer;
}






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
#include "Comet.h"


static double SampleFluxDown[200000];

void FlushElementSampling(double *sample) {
  for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
    sample[i]=0.0;
  }
}

void SampleSurfaceElement(double *x,double *sample) {
  //find cell of interest
  int i,j,k;
  long int LocalCellNumber,nface;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=NULL;  
  startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
  if (startNode->Thread==PIC::Mesh::mesh.ThisThread) {
    PIC::Mesh::cDataBlockAMR *block=startNode->block;
    LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false);
    long int FirstCellParticle=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)],ptr;
    PIC::ParticleBuffer::byte *ParticleData;
    double normalization=0.0;
    
    if (FirstCellParticle!=-1) {
      for (ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {    
	ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
	nface=Comet::GetParticleSurfaceElement(ParticleData);
	sample[nface]+=block->GetLocalParticleWeight(_H2O_SPEC_)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData)/CutCell::BoundaryTriangleFaces[nface].SurfaceArea;
      }
    }
  }
}

void PrintSampledSurfaceElementSurfaceTriangulationMesh(const char *fname,double * x,double * sample) {
#if _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  long int nface,nnode,pnode;
  class cTempNodeData {
  public:
    double Probability;
  };

  cTempNodeData *TempNodeData=new cTempNodeData[CutCell::nBoundaryTriangleNodes];

  //initialization
  for (nnode=0;nnode<CutCell::nBoundaryTriangleNodes;nnode++) TempNodeData[nnode].Probability=0.0;
  
  int i,j,k;
  long int LocalCellNumber;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=NULL;  
  startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);

  if (startNode->Thread==PIC::Mesh::mesh.ThisThread) {
    PIC::Mesh::cDataBlockAMR *block=startNode->block;
    LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false);
    long int FirstCellParticle=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)],ptr;
    PIC::ParticleBuffer::byte *ParticleData;
    double normalization=0.0;

    for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
      for (pnode=0;pnode<3;pnode++) {
	nnode=CutCell::BoundaryTriangleFaces[nface].node[pnode]-CutCell::BoundaryTriangleNodes;
	if ((nnode<0)||(nnode>=CutCell::nBoundaryTriangleNodes)) exit(__LINE__,__FILE__,"Error: out of range");
	
	TempNodeData[nnode].Probability+=sample[nface];
	normalization+=sample[nface];
      }
    }
    
    if (normalization!=0.0) for (nnode=0;nnode<CutCell::nBoundaryTriangleNodes;nnode++) TempNodeData[nnode].Probability=TempNodeData[nnode].Probability*3.0/normalization;
    
    //print the mesh
    FILE *fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\",\"Probability\"");
    fprintf(fout,"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",CutCell::nBoundaryTriangleNodes,CutCell::nBoundaryTriangleFaces);
    
    for (nnode=0;nnode<CutCell::nBoundaryTriangleNodes;nnode++) {
      fprintf(fout,"%e %e %e %e\n",CutCell::BoundaryTriangleNodes[nnode].x[0],CutCell::BoundaryTriangleNodes[nnode].x[1],CutCell::BoundaryTriangleNodes[nnode].x[2],TempNodeData[nnode].Probability);
    }
    
    for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
      fprintf(fout,"%ld %ld %ld\n",1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[0]-CutCell::BoundaryTriangleNodes),1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[1]-CutCell::BoundaryTriangleNodes),1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[2]-CutCell::BoundaryTriangleNodes));
    }
    
    fclose(fout);
    delete [] TempNodeData;
  }
#endif
}


void FlushBackfluxSampling() {
  for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
    SampleFluxDown[i]=0.0;
  }
}


void PrintBackFluxSurfaceTriangulationMesh(const char *fname) {
  long int nface,nnode,pnode;

  int rank;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);
  if (rank!=0) return;
  printf("pass MPI \n");

  class cTempNodeData {
  public:
    double NodeBackflux;
  };

  cTempNodeData *TempNodeData=new cTempNodeData[CutCell::nBoundaryTriangleNodes];

  for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
    for (pnode=0;pnode<3;pnode++) {
      nnode=CutCell::BoundaryTriangleFaces[nface].node[pnode]-CutCell::BoundaryTriangleNodes;
      if ((nnode<0)||(nnode>=CutCell::nBoundaryTriangleNodes)) exit(__LINE__,__FILE__,"Error: out of range");

      TempNodeData[nnode].NodeBackflux=SampleFluxDown[nface]/PIC::LastSampleLength;
    }
  }

   printf("Nodebackflux Done \n");
  //print the mesh
  FILE *fout=fopen(fname,"w");
  fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\",\"BackFlux\"");
  fprintf(fout,"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",CutCell::nBoundaryTriangleNodes,CutCell::nBoundaryTriangleFaces);

  for (nnode=0;nnode<CutCell::nBoundaryTriangleNodes;nnode++) {
    fprintf(fout,"%e %e %e %e\n",CutCell::BoundaryTriangleNodes[nnode].x[0],CutCell::BoundaryTriangleNodes[nnode].x[1],CutCell::BoundaryTriangleNodes[nnode].x[2],TempNodeData[nnode].NodeBackflux);
  }

  for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
    fprintf(fout,"%ld %ld %ld\n",1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[0]-CutCell::BoundaryTriangleNodes),1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[1]-CutCell::BoundaryTriangleNodes),1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[2]-CutCell::BoundaryTriangleNodes));
  }

  fclose(fout);
  delete [] TempNodeData;
}


double BulletLocalResolution(double *x) {
  int idim,i;
  double res,r,l[3];
  double SubsolarAngle;



  for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

  //3.3 AU
  if (r<1600) return 300*12.0; //*4.0 for OSIRIS

  SubsolarAngle=acos(l[0])*180.0/Pi;
  if (SubsolarAngle>=89.5) return 2000.0*pow(r/1600.0,2.0);

 //  return 300*pow(r/1600.0,2.0)*(1+5*SubsolarAngle/180.0);

  return 300*pow(r/1600.0,2.0)*12.0; //only *4.0 for OSIRIS
 
 /*  //2.7 AU
  if (r<1600) return 300;

  SubsolarAngle=acos(l[0])*180.0/Pi;
  if (SubsolarAngle>=89.5) return 700.0*pow(r/1600.0,2.0);

 //  return 300*pow(r/1600.0,2.0)*(1+5*SubsolarAngle/180.0);

 return 70*pow(r/1600.0,2.0);
 */

//  return 1200.0; //Hartley 2
//  return 200.0;
}

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace) {
  int spec=PIC::ParticleBuffer::GetI(ptr);

#if _SAMPLE_BACKFLUX_MODE_ == _SAMPLE_BACKFLUX_MODE__OFF_
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)  return _PARTICLE_DELETED_ON_THE_FACE_;
  else {
  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
  //  return _PARTICLE_DELETED_ON_THE_FACE_;
  }
#else

  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;

#endif
#else

  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  double scalar=0.0,X=0.0;
  double x[3],positionSun[3];
  int idim,i=0,j=0;
  double subSolarPointAzimuth=0.0;
  double subSolarPointZenith=0.0;
  double HeliocentricDistance=3.3*_AU_;

  double ParticleWeight,wc,LocalTimeStep;

  PIC::ParticleBuffer::byte *ParticleData;

  TriangleCutFace->GetCenterPosition(x);

  positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
  positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
  positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

  for (scalar=0.0,X=0.0,idim=0;idim<3;idim++){
    scalar+=TriangleCutFace->ExternalNormal[idim]*(positionSun[idim]-x[idim]);
  }

  if (scalar<0.0 ||  TriangleCutFace->pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) {
    while (CutCell::BoundaryTriangleFaces[j].node[0]!=TriangleCutFace->node[0] || CutCell::BoundaryTriangleFaces[j].node[1]!=TriangleCutFace->node[1] || CutCell::BoundaryTriangleFaces[j].node[2]!=TriangleCutFace->node[2]) j++;

   
    ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
   
    ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
    wc=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
   
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL;
    double x[3];
    TriangleCutFace->GetCenterPosition(x);
    node=PIC::Mesh::mesh.findTreeNode(x,node);
    LocalTimeStep=node->block->GetLocalTimeStep(spec);
 //    LocalTimeStep=Comet::CG->maxIntersectedNodeTimeStep[spec];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif
    SampleFluxDown[j]+=ParticleWeight*wc/TriangleCutFace->SurfaceArea/LocalTimeStep;
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }else{
    vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
    vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
    vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

    return _PARTICLE_REJECTED_ON_THE_FACE_;
  }
#endif

}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  return max(1.0,t->CharacteristicSize()*18.0)/1.5; //4.5
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    double CellSize;
    double CharacteristicSpeed;
    double dt;

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
    if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
      ElectricallyChargedDust::EvaluateLocalTimeStep(spec,dt,startNode); //CharacteristicSpeed=3.0;
      return 5.0* 0.3*startNode->GetCharacteristicCellSize()/50.0;
    } else {
      CharacteristicSpeed=5.0e2*sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/PIC::MolecularData::GetMass(spec));
    }
#else
    CharacteristicSpeed=5.0e2*sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/PIC::MolecularData::GetMass(spec));
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    if (spec==_H_SPEC_) CharacteristicSpeed*=30.0;
    if (spec==_O_SPEC_) CharacteristicSpeed*=10.0;
    if (spec==_H2_SPEC_) CharacteristicSpeed*=10.0;
    if (spec==_OH_SPEC_) CharacteristicSpeed*=5.0;
#endif

#endif

    CellSize=startNode->GetCharacteristicCellSize();
    return 0.3*CellSize/CharacteristicSpeed;
}

double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double res=0.0;
  /*
  bool ExternalFaces[6];
  double ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  static double v[3]={000.0,2.0e3,000.0},n=5.0E6,temp=20.0;


  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      if (spec!=_O2_SPEC_) return 0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n,temp,v,ExternalNormal,_O2_SPEC_);

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }
  */
  return res;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  /*  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={000.0,2.0e3,0.0},nNA=5.0E6,tempNA=20.0;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O2_SPEC_);

      if (ModelParticlesInjectionRate>0.0) return true;
    }
  }
  */
  return false;
}

//injection of model particles through the faces of the bounding box
long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  long int nInjectedParticles=0;
  /*  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;

  if (spec!=_O2_SPEC_) return 0; //inject only spec=0

  static double vNA[3]={000.0,2.0e3,000.0},nNA=5.0E6,tempNA=20.0;
  double v[3];


  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O2_SPEC_);


      if (ModelParticlesInjectionRate>0.0) {
        ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;

        PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

        while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
          //generate the new particle position on the face
          for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

          //generate a particle
          newParticle=PIC::ParticleBuffer::GetNewParticle();
          newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
          nInjectedParticles++;

          //generate particles' velocity
          PIC::Distribution::InjectMaxwellianDistribution(v,vNA,tempNA,ExternalNormal,_O2_SPEC_,-1);

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
  */
  return nInjectedParticles;
}

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    long int nInjectedParticles=0;

    //for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

  return nInjectedParticles;
}


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}


double Comet::GetTotalProduction(int spec,void *BoundaryElement) {
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__OFF_
  return Comet::GetTotalProductionRateBjornNASTRAN(spec)+Comet::GetTotalProductionRateUniformNASTRAN(spec)+Comet::GetTotalProductionRateJetNASTRAN(spec);
#else
  if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
   static bool InitFlag=false;
   static double MeanDustGrainMass=0.0;

   if (InitFlag==false) {
     InitFlag=true;

     //get mean mass of a dust grain
     const int nTotalTest=10000;
     double r,c,cTotal=0.0,m=0.0;


     for (int ntest=0;ntest<nTotalTest;ntest++) {
       ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(r,c);
       m+=pow(r,3.0)*c;
       cTotal+=c;
     }

     MeanDustGrainMass=m*4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity/cTotal;
   }
   return ElectricallyChargedDust::TotalMassDustProductionRate/MeanDustGrainMass;
  }
  else return Comet::GetTotalProductionRateBjornNASTRAN(spec)+Comet::GetTotalProductionRateUniformNASTRAN(spec)+Comet::GetTotalProductionRateJetNASTRAN(spec);   
#endif
}

double Comet::GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement) {
  return GetTotalProduction(spec,BoundaryElement);
}

int main(int argc,char **argv) {

  //init reading of the electron pressure data from the BATSRUS TECPLOT data file
  PIC::CPLR::DATAFILE::Offset::PlasmaElectronPressure.allocate=true;

  //init the particle solver
  PIC::InitMPI();
  PIC::Init_BeforeParser();
  Comet::Init_BeforeParser();

  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};

  //load the NASTRAN mesh
  //  CutCell::ReadNastranSurfaceMeshLongFormat("surface_Thomas_elements.nas",CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,xmin,xmax,1.0E-8);
  //CutCell::ReadNastranSurfaceMeshLongFormat("cg.Lamy-surface.nas",xmin,xmax,1.0E-8);
  // CutCell::ReadNastranSurfaceMeshLongFormat("Sphere_3dCode.nas",xmin,xmax,1.0E-8);
  //  CutCell::ReadNastranSurfaceMeshLongFormat("CG2.bdf",xmin,xmax,1.0E-8);


  //PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat("cg.RMOC.bdf");
  PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat("SHAP5_stefano.bdf",PIC::UserModelInputDataPath);
  //PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat("cg.Sphere.nas");
  PIC::Mesh::IrregularSurface::GetSurfaceSizeLimits(xmin,xmax);
  PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh("SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);

      /*
  //refine the surface mesh
  {
    char fname[_MAX_STRING_LENGTH_PIC_];
    CutCell::SmoothRefine(0.75);
    sprintf(fname,"%s/NucleusSurface-L1.dat",PIC::OutputDataFileDirectory);
    if (PIC::ThisThread==0) CutCell::PrintSurfaceTriangulationMesh(fname,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,1.0E-8);

    CutCell::SmoothRefine(0.75);
    sprintf(fname,"%s/NucleusSurface-L2.dat",PIC::OutputDataFileDirectory);
    if (PIC::ThisThread==0) CutCell::PrintSurfaceTriangulationMesh(fname,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,1.0E-8);
  }
      */

#if _3DGRAVITY__MODE_ == _3DGRAVITY__MODE__ON_
  //Computation of the gravity field for an irregular nucleus shape
  nucleusGravity::readMesh_longformat("SHAP5_Volume.bdf",PIC::UserModelInputDataPath);
  nucleusGravity::setDensity(430);
#endif   

  //set up the nastran object
  //init the nucleus
  cInternalBoundaryConditionsDescriptor CGSurfaceDescriptor;
  cInternalNastranSurfaceData *CG;

  CGSurfaceDescriptor=PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  CG=(cInternalNastranSurfaceData*) CGSurfaceDescriptor.BoundaryElement;

//  CG->PrintTitle=Comet::Sampling::OutputSurfaceDataFile::PrintTitle;
//  CG->PrintVariableList=Comet::Sampling::OutputSurfaceDataFile::PrintVariableList;

  CG->InjectionRate=Comet::GetTotalProduction;
  CG->faceat=0;
//  Nucleus->ParticleSphereInteraction=Comet::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
//  CG->InjectionBoundaryCondition=CG::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;
//  Nucleus->InjectionBoundaryCondition=Comet::InjectionBoundaryModel_Limited; ///sphereParticleInjection;
  //PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;

  Comet::GetNucleusNastranInfo(CG);

  for (int i=0;i<3;i++) {
    if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
      xmin[i]=-100.0e3,xmax[i]=100.0e3;
    }
    else {
      xmin[i]=-450.0e3,xmax[i]=450.0e3;
    }
  }

  PIC::Mesh::mesh.CutCellSurfaceLocalResolution=SurfaceResolution;
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh.memoryAllocationReport();


  //generate mesh or read from file
  char mesh[200]="amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  bool NewMeshGeneratedFlag=false;

  char fullname[STRING_LENGTH];
  sprintf(fullname,"%s/%s",PIC::UserModelInputDataPath,mesh);

  FILE *fmesh=NULL;

  fmesh=fopen(fullname,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh.readMeshFile(fullname);
  }
  else {
    NewMeshGeneratedFlag=true;

    if (PIC::Mesh::mesh.ThisThread==0) {
       PIC::Mesh::mesh.buildMesh();
       PIC::Mesh::mesh.saveMeshFile("mesh.msh");
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
       PIC::Mesh::mesh.readMeshFile("mesh.msh");
    }
  }

  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();

//  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");




  //initialize the blocks
  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();


  //init the volume of the cells'
  PIC::Mesh::IrregularSurface::CheckPointInsideDomain=PIC::Mesh::IrregularSurface::CheckPointInsideDomain_default;
  PIC::Mesh::mesh.InitCellMeasure(PIC::UserModelInputDataPath);

  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh.getMeshSignature();

    if (PIC::Mesh::mesh.ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  /*  //read the background data
    if (PIC::CPLR::DATAFILE::BinaryFileExists("CG-BATSRUS")==true)  {
      PIC::CPLR::DATAFILE::LoadBinaryFile("CG-BATSRUS");
    }
    else {
      double xminTECPLOT[3]={-32,-32,-32},xmaxTECPLOT[3]={16,32,32};

      double RotationMatrix_BATSRUS2AMPS[3][3]={ { 0.725, 0.000, 0.689}, {0.000, 1.000, 0.000}, {-0.689, 0.000, 0.725}};

    //  0.725  0.000  0.689
    //  0.000  1.000  0.000
    // -0.689  0.000  0.725

      PIC::CPLR::DATAFILE::TECPLOT::SetRotationMatrix_DATAFILE2LocalFrame(RotationMatrix_BATSRUS2AMPS);

      PIC::CPLR::DATAFILE::TECPLOT::UnitLength=1000.0;
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsXYZ(xminTECPLOT,xmaxTECPLOT);
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsSPHERICAL(0.0,500.0);

      PIC::CPLR::DATAFILE::TECPLOT::DataMode=PIC::CPLR::DATAFILE::TECPLOT::DataMode_SPHERICAL;
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedVelocityVariableData(5,1.0E3);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedIonPressureVariableData(11,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedElectronPressureVariableData(11,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedMagneticFieldVariableData(8,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedDensityVariableData(4,1.0E6);
      PIC::CPLR::DATAFILE::TECPLOT::nTotalVarlablesTECPLOT=12;
      PIC::CPLR::DATAFILE::TECPLOT::ImportData("/Users/vtenishe/Debugger/eclipse-workspace/MERCURYAMPS/AMPS/data/input/CG/3d__var_4_n00230000-extracted.plt"); //data/input/Mercury/040915-Jia/3d__var_7_t00000200_n0300072.plt

      PIC::CPLR::DATAFILE::SaveBinaryFile("CG-BATSRUS");
    }
  */

  //test the shadow procedure
  double subSolarPointAzimuth=0.0;
  double subSolarPointZenith=0.0;
  double HeliocentricDistance=3.3*_AU_;  

  double xLightSource[3]={HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith),HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith),HeliocentricDistance*cos(subSolarPointZenith)}; //{6000.0e3,1.5e6,0.0};                                                                                                                           
  PIC::Mesh::IrregularSurface::InitExternalNormalVector();
  PIC::RayTracing::SetCutCellShadowAttribute(xLightSource,false);
  PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh("SurfaceTriangulation-shadow.dat",PIC::OutputDataFileDirectory);

  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=60000; //700000;
  PIC::RequiredSampleLength=10;


  PIC::Init_AfterParser();
  PIC::Mover::Init();
  Comet::Init_AfterParser(PIC::UserModelInputDataPath);

  //init the dust particle tracing condition
  Comet::TrajectoryTracking::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=localParticleInjectionRate;

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_
  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);
#else  
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H2O_SPEC_);
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_CO2_SPEC_);
  // PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O2_SPEC_);
  // PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_CO_SPEC_);

  //init weight of the daugter products of the photolytic and electron impact reactions 
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) if (PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]<0.0) {
      double yield=0.0;
      
      /*      yield+=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_]*
	      (PhotolyticReactions::H2O::GetSpeciesReactionYield(spec)+ElectronImpact::H2O::GetSpeciesReactionYield(spec,20.0));*/
      
      yield+=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_]*
  	(PhotolyticReactions::H2O::GetSpeciesReactionYield(spec));
      
      yield/=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_];
      PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(spec,_H2O_SPEC_,yield);
    }
#endif  

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();

  //init the particle buffer
  PIC::ParticleBuffer::Init(20000000);
 

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  const int nSamplingPoints=1;

  double ProbeLocations[nSamplingPoints][DIM]={
    {2.3E3,0.0,0.0},
  };

  ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::Init(ProbeLocations,nSamplingPoints,200);
#endif


  //set the model of the boundary conditinos
  PIC::Mover::ProcessTriangleCutFaceIntersection=SurfaceBoundaryCondition;

  //set the User Definted injection function
  PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;

  //output the volume mesh
  char fname[_MAX_STRING_LENGTH_PIC_];
  sprintf(fname,"%s/VolumeMesh.dat",PIC::OutputDataFileDirectory);
//  PIC::Mesh::mesh.outputMeshTECPLOT(fname);

#if  _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  const int nLocations=2;
  static double SampleElement[nLocations][200000];

  double SampleLocations[nLocations][DIM]={
    {0.0,0.0,10.0e3},
    {10.0e3,10.0e3,10.0e3}};

  //  if (PIC::Mesh::mesh.ThisThread==0) {
    for (int i=0;i<nLocations;i++) FlushElementSampling(SampleElement[i]);
    //}
#endif
  

    PIC::Mesh::mesh.outputMeshDataTECPLOT("loaded.data.dat",0);


  int LastDataOutputFileNumber=-1;

  //the total number of iterations 
  int nTotalIterations=5400;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=50;


  for (long int niter=0;niter<nTotalIterations;niter++) {
    PIC::TimeStep();

    //update the particle tracer counter
    Comet::TrajectoryTracking::UpdateParticleCounter();

    if (PIC::ThisThread==0 && niter==1) {
      //      char fname[_MAX_STRING_LENGTH_PIC_];
      
      //sprintf(fname,"%s/SurfaceTriangulation_Jet.dat",PIC::OutputDataFileDirectory);
      //Comet::PrintSurfaceTriangulationMesh(fname,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,1.0E-8);


#if _COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE_ == _COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE__ON_
      char fname3[_MAX_STRING_LENGTH_PIC_];

      sprintf(fname3,"%s/SurfaceTriangulation_MaxLiftableSize.dat",PIC::OutputDataFileDirectory);
      Comet::PrintMaxLiftableSizeSurfaceTriangulationMesh(fname3);
#endif
    }

#if  _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
    //    if (PIC::Mesh::mesh.ThisThread==0) {
      for (int i=0;i<nLocations;i++) SampleSurfaceElement(SampleLocations[i],SampleElement[i]);
      // }
#endif

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
#if _SAMPLE_BACKFLUX_MODE_ == _SAMPLE_BACKFLUX_MODE__ON_      
      char fname2[_MAX_STRING_LENGTH_PIC_];

      if (PIC::Mesh::mesh.ThisThread==0) cout << "Printing backflux sampling output " << endl;
      
      sprintf(fname2,"%s/SurfaceTriangulation_Backflux.dat",PIC::OutputDataFileDirectory);
      if (PIC::Mesh::mesh.ThisThread==0) PrintBackFluxSurfaceTriangulationMesh(fname2);

      FlushBackfluxSampling();
#endif      

#if  _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
      for (int counter=0;counter<nLocations;counter++) {
	char fnameElement[_MAX_STRING_LENGTH_PIC_];

	if (PIC::Mesh::mesh.ThisThread==0) cout << "Printing Sampled Surface Element output " << endl;
	
	sprintf(fnameElement,"%s/SurfaceTriangulation_SampledElement_%i.dat",PIC::OutputDataFileDirectory,counter);
	PrintSampledSurfaceElementSurfaceTriangulationMesh(fnameElement,SampleLocations[counter],SampleElement[counter]);
	
	FlushElementSampling(SampleElement[counter]);
      }
#endif

      PIC::RequiredSampleLength*=2;
      if (PIC::RequiredSampleLength>700) PIC::RequiredSampleLength=700;

      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
      if (PIC::Mesh::mesh.ThisThread==0) cout << "The new lample length is " << PIC::RequiredSampleLength << endl;
    }
    if (PIC::ThisThread==0)   cout << "niter:" << niter << endl;
#if _PIC_MODEL__RADIATIVECOOLING__MODE_ == _PIC_MODEL__RADIATIVECOOLING__MODE__CROVISIER_
	Comet::StepOverTime();
#endif
  }

  //output the particle statistics for the nightly tests
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];

    sprintf(fname,"%s/test_CG.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }


  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}

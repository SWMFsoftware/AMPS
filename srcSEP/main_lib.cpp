

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
#include "pic.h"
#include "constants.h"
#include "sep.h"




//the parameters of the domain and the sphere

const double DebugRunMultiplier=4.0;
const double rSphere=_RADIUS_(_TARGET_);


const double xMaxDomain=250;

const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*4.0*1.0/100/2.5,dxMaxSphere=DebugRunMultiplier*2.0/10.0;




double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
   //delete all particles that was not reflected on the surface
   //PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;
}



void amps_init() {
  PIC::InitMPI();


  //SetUp the alarm
//  PIC::Alarm::SetAlarm(2000);


  rnd_seed();
  MPI_Barrier(MPI_COMM_WORLD);


  //reserve data for magnetic filed
  PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;
  //PIC::CPLR::DATAFILE::Offset::MagneticField.active=true;
  //PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  //PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=PIC::CPLR::DATAFILE::Offset::MagneticField.nVars*sizeof(double);


  //init the Mercury model
 ////::Init_BeforeParser();
  PIC::Init_BeforeParser();

//  ProtostellarNebula::OrbitalMotion::nOrbitalPositionOutputMultiplier=10;
///  ProtostellarNebula::Init_AfterParser();



  //register the sphere
  {
    double sx0[3]={0.0,0.0,0.0};
    cInternalBoundaryConditionsDescriptor SphereDescriptor;
    cInternalSphericalData *Sphere;


    //reserve memory for sampling of the surface balance of sticking species
    long int ReserveSamplingSpace[PIC::nTotalSpecies];

    for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


    cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);

    PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
    SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
    Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
    Sphere->SetSphereGeometricalParameters(sx0,rSphere);


    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/Sphere.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceMesh(fname);

    sprintf(fname,"%s/SpheraData.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceData(fname,0);


    Sphere->localResolution=SEP::Mesh::localSphericalSurfaceResolution;
    Sphere->InjectionRate=SEP::ParticleSource::InnerBoundary::sphereInjectionRate;
    Sphere->faceat=0;
    Sphere->ParticleSphereInteraction=ParticleSphereInteraction;
    Sphere->InjectionBoundaryCondition=SEP::ParticleSource::InnerBoundary::sphereParticleInjection;

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
  list<SEP::cFieldLine> field_line;
  double xStart[3]={1.1,0.0,0.0};

  SEP::ParkerSpiral::CreateFileLine(&field_line,xStart,250.0);
  SEP::Mesh::LoadFieldLine(&field_line);
  PIC::Mesh::mesh.UserNodeSplitCriterion=SEP::Mesh::NodeSplitCriterion;

  //generate only the tree
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,SEP::Mesh::localResolution);
  PIC::Mesh::mesh.memoryAllocationReport();


  if (PIC::Mesh::mesh.ThisThread==0) {
    PIC::Mesh::mesh.buildMesh();
    PIC::Mesh::mesh.saveMeshFile("mesh.msh");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    MPI_Barrier(MPI_COMM_WORLD);
    PIC::Mesh::mesh.readMeshFile("mesh.msh");
  }

  cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh.ThisThread << endl;

  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif

  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();

  //initialize the blocks
  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif

  //init the volume of the cells'
  PIC::Mesh::mesh.InitCellMeasure();





//init the PIC solver
  PIC::Init_AfterParser ();
	PIC::Mover::Init();


  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=SEP::Mesh::localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=SEP::ParticleSource::OuterBoundary::BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=SEP::ParticleSource::OuterBoundary::BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();

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

  InitMagneticField(PIC::Mesh::mesh.rootTree);
  PIC::Mesh::mesh.outputMeshDataTECPLOT("magnetic-field.dat",0);


  MPI_Barrier(MPI_COMM_WORLD);
  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;

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

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  
  int idim;

    //determine the parameters of the orbital motion of Mercury
    SpiceDouble StateBegin[6],StateEnd[6],lt,StateSun[6],StateMiddle[6];
    double lBegin[3],rBegin,lEnd[3],rEnd,vTangentialBegin=0.0,vTangentialEnd=0.0,c0=0.0,c1=0.0;

    SpiceDouble HCI_to_MSO_TransformationMartix[6][6];

    spkezr_c("Mercury",Mercury::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateBegin,&lt);
    spkezr_c("SUN",Mercury::OrbitalMotion::et,"MSGR_MSO","none","Mercury",StateSun,&lt);

    //calculate Mercury's velocity in an itertial frame, which have dirtectional vectors that coinsides with that of MSO
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et+0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],HCI_to_MSO_TransformationMartix);
    spkezr_c("Mercury",Mercury::OrbitalMotion::et+0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],"MSGR_HCI","none","SUN",StateMiddle,&lt);


    Mercury::OrbitalMotion::et+=PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];
    spkezr_c("Mercury",Mercury::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateEnd,&lt);


    for (rBegin=0.0,rEnd=0.0,idim=0;idim<3;idim++) {
      StateBegin[idim]*=1.0E3,StateBegin[3+idim]*=1.0E3;
      StateEnd[idim]*=1.0E3,StateEnd[3+idim]*=1.0E3;
      StateMiddle[idim]*=1.0E3,StateMiddle[3+idim]*=1.0E3;

      rBegin+=pow(StateBegin[idim],2);
      rEnd+=pow(StateEnd[idim],2);

      Mercury::xObject_HCI[idim]=StateBegin[idim];
      Mercury::vObject_HCI[idim]=StateBegin[3+idim];

      Mercury::xSun_SO[idim]=1.0E3*StateSun[idim];
      Mercury::vSun_SO[idim]=1.0E3*StateSun[3+idim];
    }

    //calculate parameters of SO_FROZEN
    //velocity of the coordinate frame
    for (idim=0;idim<3;idim++) {
      Mercury::vObject_SO_FROZEN[idim]=
          HCI_to_MSO_TransformationMartix[idim][0]*StateMiddle[3+0]+
          HCI_to_MSO_TransformationMartix[idim][1]*StateMiddle[3+1]+
          HCI_to_MSO_TransformationMartix[idim][2]*StateMiddle[3+2];
    }

    //the axis of rotation of the MSO fraim in MSO_FROZEN during the next time step
    //get pointing direction to the Sun at the end of the current iteration in MSO_FROZEN
    SpiceDouble fmatrix[6][6];
    double SunPointingDirectionEnd[3],SunPointingDirectionEnd_MSO_FROZEN[3];

    //calculate Sun pointing at the end of the iteration in HCI frame (et is already incremented!!!!!!)
    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et,fmatrix);

    SunPointingDirectionEnd[0]=fmatrix[0][0];
    SunPointingDirectionEnd[1]=fmatrix[1][0];
    SunPointingDirectionEnd[2]=fmatrix[2][0];

    //convert the pointing direction vector into MSO_FROZEN frame
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],fmatrix);

    for (idim=0;idim<3;idim++) {
      SunPointingDirectionEnd_MSO_FROZEN[idim]=
          fmatrix[idim][0]*SunPointingDirectionEnd[0]+
          fmatrix[idim][1]*SunPointingDirectionEnd[1]+
          fmatrix[idim][2]*SunPointingDirectionEnd[2];
    }

    //calculate the rate of rotation in MSO_FROZEN
    Mercury::RotationRate_SO_FROZEN=acos(SunPointingDirectionEnd_MSO_FROZEN[0])/PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];


    //calculate the direction of rotation
    double c=sqrt(pow(SunPointingDirectionEnd_MSO_FROZEN[1],2)+pow(SunPointingDirectionEnd_MSO_FROZEN[2],2));

    if (c>0.0) {
      Mercury::RotationVector_SO_FROZEN[0]=0.0;
      Mercury::RotationVector_SO_FROZEN[1]=-SunPointingDirectionEnd_MSO_FROZEN[2]/c*Mercury::RotationRate_SO_FROZEN;
      Mercury::RotationVector_SO_FROZEN[2]=SunPointingDirectionEnd_MSO_FROZEN[1]/c*Mercury::RotationRate_SO_FROZEN;
    }
    else {
      Mercury::RotationVector_SO_FROZEN[0]=0.0;
      Mercury::RotationVector_SO_FROZEN[1]=0.0;
      Mercury::RotationVector_SO_FROZEN[2]=0.0;
    }


    //RECALCUALTE THE ROTATION VECTOR USING THE TRANSOFRMATON MARTICX FROM MSO_FROSEN at the time step (n) to the MSO_FROZEN at the time step (n+1)
    //the rotation vector is the eigrnvector of the transformation matrix
    //Zhuravlev, Osnovy teoreticheskoi mehaniki, Chapter 2, paragraph 6.2 (sposoby zadaniya orientacii tverdogo tela)

    //get the transformation matrix T(LSO[n]->LSO[n+1])=T1(LSO[n]->MSGR_HCI)*T2(MSGR_HCI->LSO[n+1])
    SpiceDouble T1[6][6],T2[6][6];
    double T[3][3];
    int i,j,k;



    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et,T1);
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T2);




    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      T[i][j]=0.0;

      for (k=0;k<3;k++) T[i][j]+=T2[i][k]*T1[k][j];
    }

    //determine the rate and the vectrot of the rotation
    double RotationAngle,t,RotationVector[3],RotationRate;

    RotationAngle=acos((T[0][0]+T[1][1]+T[2][2]-1.0)/2.0);

    t=2.0*sin(RotationAngle);
    RotationVector[0]=(T[2][1]-T[1][2])/t;
    RotationVector[1]=(T[0][2]-T[2][0])/t;
    RotationVector[2]=(T[1][0]-T[0][1])/t;

    RotationRate=RotationAngle/PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];

//    t=RotationRate/sqrt(RotationVector[0]*RotationVector[0]+RotationVector[1]*RotationVector[1]+RotationVector[2]*RotationVector[2]);
//    RotationVector[0]*=t,RotationVector[1]*=t,RotationVector[2]*=t;

t=1.0;

    //TEST THE ROTATION RATE AND THE ROTATION VECTOR
    double testRoptationMatrix[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    double cosRotationAngle,sinRotationAngle;

    cosRotationAngle=cos(RotationAngle);
    sinRotationAngle=sin(RotationAngle);

    for (i=0;i<3;i++) testRoptationMatrix[i][i]+=cosRotationAngle;
    for (i=0;i<3;i++) for (j=0;j<3;j++) testRoptationMatrix[i][j]+=(1.0-cosRotationAngle)*RotationVector[i]*RotationVector[j]/pow(t,2);

    testRoptationMatrix[0][1]-=sinRotationAngle*RotationVector[2]/t,testRoptationMatrix[0][2]+=sinRotationAngle*RotationVector[1]/t;
    testRoptationMatrix[1][0]+=sinRotationAngle*RotationVector[2]/t,testRoptationMatrix[1][2]-=sinRotationAngle*RotationVector[0]/t;
    testRoptationMatrix[2][0]-=sinRotationAngle*RotationVector[1]/t,testRoptationMatrix[2][1]+=sinRotationAngle*RotationVector[0]/t;


    //CALCULATE THE EFECT OF THE TRANSFORMATION AND TRANSTER THE TRANSFORMED VEWCTROS TO HCI
    double T3[3][3];

    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T2);

    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      T3[i][j]=0.0;

      for (k=0;k<3;k++) T3[i][j]+=T2[i][k]*testRoptationMatrix[k][j];
    }


    //GET THE ROTATION MATRIX FROM SPICE
//    SpiceDouble rVect1[3],rVect2[3],rVect[3],rot[3][3],xform[6][6];

//    sxform_c (  "MSGR_HCI","MSGR_MSO", et, tsipm ) ;

    //RECALCULATE THE MATRIX AGIN
    double newRotationVector[3],newRate;

    newRate=Exosphere::OrbitalMotion::FrameRotation::GetRotationVector(newRotationVector,"MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],Mercury::OrbitalMotion::et);


    //END OF TRANSFORMATIUON TESTS ------------------

    rBegin=sqrt(rBegin);
    rEnd=sqrt(rEnd);

    for (idim=0;idim<3;idim++) {
      lBegin[idim]=StateBegin[idim]/rBegin;
      lEnd[idim]=StateEnd[idim]/rEnd;

      c0+=StateBegin[3+idim]*lBegin[idim];
      c1+=StateEnd[3+idim]*lEnd[idim];
    }

    Mercury::xObjectRadial=0.5*(rBegin+rEnd);
    Mercury::vObjectRadial=0.5*(c0+c1);

    //calculate TAA
    Mercury::OrbitalMotion::TAA=Mercury::OrbitalMotion::GetTAA(Mercury::OrbitalMotion::et);

    for (idim=0;idim<3;idim++) {
      vTangentialBegin+=pow(StateBegin[3+idim]-c0*lBegin[idim],2);
      vTangentialEnd+=pow(StateEnd[3+idim]-c1*lEnd[idim],2);
    }

    vTangentialBegin=sqrt(vTangentialBegin);
    vTangentialEnd=sqrt(vTangentialEnd);

    Mercury::OrbitalMotion::CoordinateFrameRotationRate=0.5*(vTangentialBegin/rBegin+vTangentialEnd/rEnd);


    //determine direction to the Sun and rotation angle in the coordiname frame related to Mercury
    SpiceDouble state[6],l=0.0;

    spkezr_c("SUN",Mercury::OrbitalMotion::et,"IAU_MERCURY","none","MERCURY",state,&lt);

    for (idim=0;idim<3;idim++) l+=pow(state[idim],2);

    for (l=sqrt(l),idim=0;idim<3;idim++) {
      Mercury::OrbitalMotion::SunDirection_IAU_OBJECT[idim]=state[idim]/l;
    }
    
    //matrixes for tranformation MSO->IAU and IAU->MSO coordinate frames
    sxform_c("MSGR_MSO","IAU_MERCURY",Mercury::OrbitalMotion::et,Mercury::OrbitalMotion::SO_to_IAU_TransformationMartix);
    sxform_c("IAU_MERCURY","MSGR_MSO",Mercury::OrbitalMotion::et,Mercury::OrbitalMotion::IAU_to_SO_TransformationMartix);
 
#endif

    //make the time advance
    static int LastDataOutputFileNumber=0;

    //make the time advance
     PIC::TimeStep();

     // write output file
     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;
       
       
       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }
     
}




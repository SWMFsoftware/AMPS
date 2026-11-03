//$Id: main.cpp,v 1.12 2016/09/11 20:39:17 yunilee Exp $



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

//the species
//unsigned int _C_SPEC_=0;

//forward scattering cross section
#include "ForwardScatteringCrossSection.h"


//the particle class
#include "rnd.h"
#include "pic.h"
#include "Mars.h"
#include "MTGCM.h"



//the parameters of the domain and the sphere

const double DebugRunMultiplier=2.0;


const double rSphere=_RADIUS_(_TARGET_);
const double xMaxDomain=6.0;
const double dxMinGlobal=2.0,dxMaxGlobal=5.0;
const double dxMinSphere=DebugRunMultiplier*2.0/100,dxMaxSphere=DebugRunMultiplier*4.0/100.0;


//the total acceleration of the particles
//#include "species/Na.h"


//double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
/*
void TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {


  accl[0]=0.0,accl[1]=0.0,accl[2]=0.0;

//  return;

  //the gravity force
  double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  double r=sqrt(r2);
  int idim;

  for (idim=0;idim<DIM;idim++) {
    accl[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x[idim]/r;
  }



}
*/



//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
  double res,r,l[3];
  int idim;
  double SubsolarAngle;

  for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

  SubsolarAngle=acos(l[0]);


  SubsolarAngle=0.0;

  res=dxMinSphere+(dxMaxSphere-dxMinSphere)/Pi*SubsolarAngle;



//  return rSphere*res;

  return DebugRunMultiplier*800.0E3;
}

double localResolution(double *x) {
  int idim;
  double lnR,res=0.0,r=0.0;

  for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

  r=sqrt(r);

  if ((3.5E6<r)&&(r<3.8E6)) return DebugRunMultiplier*100.0E3;
  else return DebugRunMultiplier*rSphere*dxMaxGlobal;


  /*
  if (r>dxMinGlobal*rSphere) {
    lnR=log(r);
    res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
  }
  else res=dxMinGlobal;

*/

  return DebugRunMultiplier*rSphere*res;
}

//set up the local time step
//determine the limits on the "geometrical (flight time) times step

void GetTimeStepLimits(double &minTimeStep,double &maxTimeStep,int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double dt;

  const double CharacteristicSpeed_NA=6.0E3;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    dt=0.3*startNode->GetCharacteristicCellSize()/CharacteristicSpeed_NA;

    if (dt>maxTimeStep) maxTimeStep=dt;
    if ((minTimeStep<0.0)||(dt<minTimeStep)) minTimeStep=dt;
  }
  else for (int i=0;i<(1<<_MESH_DIMENSION_);i++) if (startNode->downNode[i]!=NULL) GetTimeStepLimits(minTimeStep,maxTimeStep,spec,startNode->downNode[i]);
}



double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize,dt;

  const double CharacteristicSpeed_NA=6.0E3;
  const double minParticleMovingInterval=0.1E3;

  const double TimeStepMinValue=minParticleMovingInterval/CharacteristicSpeed_NA;



  static double *minTimeStep=NULL,*maxTimeStep=NULL;

  if (minTimeStep==NULL) {
    //init the buffers

    if (PIC::ThisThread==0) cout << "Determine the time step limits:" << endl;

    minTimeStep=new double[PIC::nTotalSpecies];
    maxTimeStep=new double[PIC::nTotalSpecies];

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      minTimeStep[s]=-1.0,maxTimeStep[s]=-1.0;
      GetTimeStepLimits(minTimeStep[s],maxTimeStep[s],s,PIC::Mesh::mesh->rootTree);

      if (PIC::ThisThread==0){
        cout << "spec=" << s << ", minTimeStep=" << minTimeStep[s] << ", maxTimeStep=" << maxTimeStep[s] << ": The time step range that will be actually used is (" << TimeStepMinValue << "," << maxTimeStep[s] << ")" << endl;
      }
    }
  }

  //get the local "geometrical" time step
  CellSize=startNode->GetCharacteristicCellSize();
  dt=0.3*CellSize/CharacteristicSpeed_NA;

  //modify the "geometrical" time step
  dt=TimeStepMinValue+(maxTimeStep[spec]-TimeStepMinValue)/(maxTimeStep[spec]-minTimeStep[spec])*(dt-minTimeStep[spec]);

  if (dt<=0.0) {
    exit(__LINE__,__FILE__,"Error: the time step is out of range");
  }

  return dt;
}


bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {


  return false;
}








double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {

   //delete all particles that was not reflected on the surface
//   PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;
}



void OxigenTGCM() {
/*cDataSetMTGCM O,COp;

 O.PlanetRadius=3388.25E3;
 O.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
 O.ReadDataFile("../data/input/Mars/MTGCM_equinox_SL/O.h");
    COp.PlanetRadius=3388.25E3;
    COp.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    COp.ReadDataFile("../data/input/Mars/MTGCM_equinox_SL/COp.h");

 double t,x[3]={0.0,550.2E3+O.PlanetRadius,0.0};
// t=O.Interpolate(x);


 const int nPoints=300;
 double R=135.0E3+3388.25E3;
 const double dLon=2.0*Pi/(nPoints-1),dLat=Pi/(nPoints-1);

 if (PIC::ThisThread==0) {
 FILE *fout=fopen("interpolation.dat","w");
 fprintf(fout,"VARIABLES = \"Alt\", \"t\"\nZONE I=%i, J=%i, DATAPACKING=POINT\n",nPoints,nPoints);

 int i,j;
 double Lon,Lat;

 //for (i=0;i<nPoints;i++) {
   //Lat=-Pi/2.0+i*dLat;

   for (j=0;j<10;j++) {
     Lat=-Pi+j*dLat;
       Lat=60.0*Pi/180;

     x[0]=-R*cos(Lat);
     x[1]=0;//R*cos(Lat)*sin(Lon);
     x[2]=R*sin(Lat);

       t=COp.Interpolate(x);
     
   //    cout <<(R-3388.25E3)/1000<<"  "<<t<<"\n"<< endl;

     
     fprintf(fout,"%e  %e\n",(R-3388.25E3)/1000,t);
     
    R+=1E3;
   }
 

 fclose(fout);
}
  //  exit(__LINE__,__FILE__,"test");*/
cDataSetMGITM Te,Tn,Ti,O2p,E,O,CO2,CO,N2;
    const std::string fname=newMars::GetMGITMInputFilePath();
    
    Te.PlanetRadius=_RADIUS_(_TARGET_);
    Te.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
    Te.ReadDataFile(_Te_MGITM_,fname.c_str());

    Tn.PlanetRadius=_RADIUS_(_TARGET_);
    Tn.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
    Tn.ReadDataFile(_Tn_MGITM_,fname.c_str());

    Ti.PlanetRadius=_RADIUS_(_TARGET_);
    Ti.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
    Ti.ReadDataFile(_Ti_MGITM_,fname.c_str());

    O2p.PlanetRadius=_RADIUS_(_TARGET_);
    O2p.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    //O2p.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT_;
    //sprintf(fname,"%s/MTGCM_equinox_SL/MGITM_170910_153001MAMPS.dat",PIC::UserModelInputDataPath);
    O2p.ReadDataFile(_nO2P_MGITM_,fname.c_str());
        
    E.PlanetRadius=_RADIUS_(_TARGET_);
    E.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    //E.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT_;
    //sprintf(fname,"%s/MTGCM_equinox_SL/MGITM_170910_153001MAMPS.dat",PIC::UserModelInputDataPath);
    //sprintf(fname,"%s/MTGCM_equinox_SL/MGITM_170910_180001MAMPS.dat",PIC::UserModelInputDataPath);
    E.ReadDataFile(_Ne_MGITM_,fname.c_str());

    CO2.PlanetRadius=_RADIUS_(_TARGET_);
    CO2.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    //sprintf(fname,"%s/MTGCM_equinox_SL/MGITM_170910_153001MAMPS.dat",PIC::UserModelInputDataPath);
    //sprintf(fname,"%s/MTGCM_equinox_SL/MGITM_170910_180001MAMPS.dat",PIC::UserModelInputDataPath);
    CO2.ReadDataFile(_nCO2_MGITM_,fname.c_str());
        
    O.PlanetRadius=_RADIUS_(_TARGET_);
    O.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    //sprintf(fname,"%s/MTGCM_equinox_SL/MGITM_170910_153001MAMPS.dat",PIC::UserModelInputDataPath);
    //sprintf(fname,"%s/MTGCM_equinox_SL/MGITM_170910_180001MAMPS.dat",PIC::UserModelInputDataPath);
    O.ReadDataFile(_nO_MGITM_,fname.c_str());

    N2.PlanetRadius=_RADIUS_(_TARGET_);
    N2.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    N2.ReadDataFile(_nN2_MGITM_,fname.c_str());
        
    CO.PlanetRadius=_RADIUS_(_TARGET_);
    CO.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    CO.ReadDataFile(_nCO_MGITM_,fname.c_str());
        
 double t,t2,x[3]={0.0,0.0,0.0},t3;
// t=O.Interpolate(x);

 const int nPoints1=72, nPoints2=72;
 double R=98.75E3+3388.25E3;
 const double dLon=2.0*Pi/(nPoints1-1),dLat=Pi/(nPoints2-1);

 if (PIC::ThisThread==0) {
 FILE *fout=fopen("interpolation.dat","w");
 fprintf(fout,"VARIABLES = \"Lon\",\"Lat\",\"Alt\",\"Te\",\"Tn\",\"Ti\",\"O2p\",\"E\",\"E-ori\",\"nO\",\"nCO2\",\"nCO\",\"nN2\",\"production rate\"  \nZONE I=%i, J=%i, K=81 DATAPACKING=POINT\n",nPoints1,nPoints2);
 //fprintf(fout,"VARIABLES = \"Alt\", \"O2p\",\"E\",\"E-ori\"\n");

 int i,j,k;
 double Lon,Lat;
     double tE,tO2p,production_rate;
     double r;//=sqrt(r2);
     double Altitude;// = r-_RADIUS_(_TARGET_);
     double cutoffAlt=225.0E3;
     double delZ=30.0E3;

//for (k=0;k<121;k++) {
 //for (i=0;i<nPoints2;i++) {
  // Lat=-Pi/2.0+i*dLat;
  // Lat=(-87.5+i*5.0)*Pi/180;

   for (k=0;k<161;k++) {
    for (i=0;i<nPoints1;i++) {
      for (j=0;j<nPoints2;j++) {
        Lat=-Pi/2+i*dLat;
        Lon=(2.5+j*5.0)*Pi/180;
       // Lon=(-177.5+j*5.0)*Pi/180;
       // Lon=206*Pi/180;
       // Lat=-71*Pi/180;

        x[0]=R*cos(Lat)*cos(Lon);
        x[1]=R*cos(Lat)*sin(Lon);
        x[2]=R*sin(Lat);

        r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        Altitude=r-_RADIUS_(_TARGET_);

        if (Altitude>=cutoffAlt) {
           
           
               double vecMGITMtop[3]={0.0,0.0,0.0},veccutoff[3]={0.0,0.0,0.0},vectop[3]={0.0,0.0,0.0};
               int idim;
               
               for (idim=0;idim<3;idim++) {
                   vectop[idim] = (_RADIUS_(_TARGET_)+(Altitude))*(x[idim]/r);
                   veccutoff[idim]= (_RADIUS_(_TARGET_)+(Altitude+delZ))*(x[idim]/r);
                   vecMGITMtop[idim]= (_RADIUS_(_TARGET_)+(298750))*(x[idim]/r);
               }
               
               double localHO2p=-delZ/(log(O2p.Interpolate(veccutoff)/O2p.Interpolate(x)));
               double localHE=-delZ/(log(E.Interpolate(veccutoff)/E.Interpolate(x)));
               
               //bool isnan(localHE);
               
               if (localHO2p>0.0) {
                 if (Altitude<298750){
                    tE=E.Interpolate(x)*exp(-(Altitude-cutoffAlt)/localHO2p);
                   }
                  else if (Altitude >=298750){
                    tE=E.Interpolate(vecMGITMtop)*exp(-(Altitude-cutoffAlt)/localHO2p);
                   }
    
                }
        
               else if (localHO2p<=0.0) {tE=E.Interpolate(x);}

               
               tO2p=O2p.Interpolate(x);
               
               // if (Altitude>=298750){
                if (tE<tO2p) {tE=tO2p;}
            //       else if ((tE/tO2p)>1.0E1) {tE=tO2p;}
                else if ((tE/tO2p)>2.0E1 && Altitude>240000) {tE=tO2p;}
            //  }
        }

  else{
  tE=E.Interpolate(x);}

     
        if (Te.Interpolate(x)<=1200.0) {production_rate=1.95E-7*pow(300/Te.Interpolate(x),0.7)*(O2p.Interpolate(x)/1.0e6)*(tE/1.0e6);}
        else if (Te.Interpolate(x)>1200.0) {production_rate=7.39E-8*pow(1200/Te.Interpolate(x),0.56)*(O2p.Interpolate(x)/1.0e6)*(tE/1.0e6);}
       
       t=O2p.Interpolate(x);
       t2=tE;//E.Interpolate(x);
       t3=E.Interpolate(x);
     
   //    cout <<(R-3388.25E3)/1000<<"  "<<t<<"\n"<< endl;

     
     // fprintf(fout,"%e  %e  %e  %e\n",(R-3388.25E3)/1000,t,t2,t3);
     fprintf(fout,"%e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",Lon*180/Pi,Lat*180/Pi,Altitude/1000,Te.Interpolate(x),Tn.Interpolate(x),Ti.Interpolate(x),t,t2,t3,O.Interpolate(x),CO2.Interpolate(x),CO.Interpolate(x),N2.Interpolate(x),production_rate);
     }
  }

    R+=2.5E3;
   }
 

 fclose(fout);
}
  //  exit(__LINE__,__FILE__,"test");
}

int main(int argc,char **argv) {
//  MPI_Init(&argc,&argv);
PIC::InitMPI();

  rnd_seed();


 newMars::HotOxygen::init();



  //==================  set up the PIC solver

//  char inputFile[]="mars.input";


 // MPI_Barrier(MPI_COMM_WORLD);
MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //set up the alarm
 // PIC::Alarm::SetAlarm(8.0*3600.0-15*60);

  //init the particle solver
  PIC::Init_BeforeParser();
 // PIC::Parser::Run(inputFile);
		//  PIC::Init_AfterParser();


 // PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=100;
 // PIC::RequiredSampleLength=200; //0;


  //================== print TGCM solution
   OxigenTGCM();
    
    
    
   newMars::Init_AfterParser();




  //register the sphere
  double sx0[3]={0.0,0.0,0.0};
  cInternalBoundaryConditionsDescriptor SphereDescriptor;
  cInternalSphericalData *Sphere;

  PIC::BC::InternalBoundary::Sphere::Init();
  SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
  Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
  Sphere->SetSphereGeometricalParameters(sx0,rSphere);


  Sphere->PrintSurfaceMesh("Sphere.dat");
  Sphere->PrintSurfaceData("SpheraData.dat",0);
  Sphere->localResolution=localSphericalSurfaceResolution;
  Sphere->InjectionRate=NULL;
  Sphere->faceat=0;
  Sphere->ParticleSphereInteraction=ParticleSphereInteraction;
  Sphere->InjectionBoundaryCondition=NULL;




  //init the solver
  PIC::Mesh::initCellSamplingDataBuffer();

  //init the mesh
  cout << "Init the mesh" << endl;

  int maxBlockCellsnumber,minBlockCellsnumber,idim;

  maxBlockCellsnumber=_BLOCK_CELLS_X_;
  if (DIM>1) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Y_);
  if (DIM>2) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Z_);

  minBlockCellsnumber=_BLOCK_CELLS_X_;
  if (DIM>1) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Y_);
  if (DIM>2) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Z_);

  double DomainLength[3],DomainCenterOffset[3],xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

  if (maxBlockCellsnumber==minBlockCellsnumber) {
    for (idim=0;idim<DIM;idim++) {
      DomainLength[idim]=2.0*xMaxDomain*rSphere;
      DomainCenterOffset[idim]=-xMaxDomain*rSphere;
    }
  }
  else {
    if (maxBlockCellsnumber!=_BLOCK_CELLS_X_) exit(__LINE__,__FILE__);
    if (minBlockCellsnumber!=_BLOCK_CELLS_Y_) exit(__LINE__,__FILE__);
    if (minBlockCellsnumber!=_BLOCK_CELLS_Z_) exit(__LINE__,__FILE__);

    DomainLength[0]=xMaxDomain*rSphere*(1.0+double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_);
    DomainLength[1]=DomainLength[0]*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
    DomainLength[2]=DomainLength[0]*double(_BLOCK_CELLS_Z_)/_BLOCK_CELLS_X_;

    DomainCenterOffset[0]=-xMaxDomain*rSphere*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
    DomainCenterOffset[1]=-DomainLength[1]/2.0;
    DomainCenterOffset[2]=-DomainLength[2]/2.0;
  }

  for (idim=0;idim<DIM;idim++) {
    xmin[idim]=DomainCenterOffset[idim];
    xmax[idim]=DomainLength[idim]+DomainCenterOffset[idim];
  }



  //generate only the tree
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,localResolution);
  PIC::Mesh::mesh->memoryAllocationReport();



  if (PIC::Mesh::mesh->ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile("mesh.msh");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    MPI_Barrier(MPI_COMM_WORLD);
    PIC::Mesh::mesh->readMeshFile("mesh.msh");
  }


  cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh->ThisThread << endl;

//  PIC::Mesh::mesh->outputMeshTECPLOT("mesh.dat");

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

//  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);

  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();

  //initialize the blocks
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

//  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);

  //init the volume of the cells'
  PIC::Mesh::mesh->InitCellMeasure();

  //PIC::Sampling::minIterationNumberForDataOutput=15000;


  //============================== DEBUG =========================
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];

    while (node!=NULL) {
      if (node->Temp_ID==2478) {
        cout << __FILE__<< "@" << __LINE__ << endl;

        PIC::Mesh::mesh->InitCellMeasure(node);
      }

      node=node->nextNodeThisThread;
    }

  //============================== END DEBUG ====================================


  //set up the volume particle injection
//  PIC::VolumeParticleInjection::VolumeInjectionRate=newMars::ProductionRateCaluclation;



  //  PIC::VolumeParticleInjection::RegisterVolumeInjectionProcess(newMars::ProductionRateCaluclation,newMars::HotOxygen::HotOProduction,newMars::LocalTimeStep);
//    PIC::VolumeParticleInjection::RegisterVolumeInjectionProcess(newMars::ProductionRateCaluclation,newMars::HotCarbon::HotCProduction,newMars::LocalTimeStep);
    PIC::VolumeParticleInjection::RegisterVolumeInjectionProcess(newMars::ProductionRateCaluclation,newMars::HotAtomProduction_wrapper,newMars::LocalTimeStep);

//  PIC::ChemicalReactions::PhotolyticReactions::Init();
//PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
  
  /*
  //init the interpolation procedure
  newMars::ReadMTGCM();
  */



  //set up the output of the mars model production rate
  PIC::Mesh::PrintVariableListCenterNode.push_back(newMars::PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(newMars::PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(newMars::Interpolate);



//init the PIC solver
  PIC::Init_AfterParser();
  PIC::Mover::Init();
//  PIC::Mover::TotalParticleAcceleration=TotalParticleAcceleration;

  /*
  for (int s=0;s<PIC::nTotalSpecies;s++) {
    PIC::Mover::MoveParticleTimeStep[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder; ///UniformWeight_UniformTimeStep_SecondOrder;
    PIC::Mover::MoveParticleBoundaryInjection[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder;
  }
  */


  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //set up the particle weight
  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=3000;//*PIC::nTotalThreads;
//  if (PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber>5000) PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=500; //50000;

  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=NULL;
    if (_C_SPEC_>=0) {
        PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_C_SPEC_);
    }
   
    if (_O_SPEC_>=0) {
        PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O_SPEC_);
    }



//  PIC::Mesh::mesh->outputMeshTECPLOT("mesh.dat");
//  PIC::Mesh::mesh->outputMeshDataTECPLOT("mesh.data.dat",_C_SPEC_);

  MPI_Barrier(MPI_COMM_WORLD);
  if (PIC::Mesh::mesh->ThisThread==0) cout << "The mesh is generated" << endl;





  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=NULL;
  PIC::BC::InitBoundingBoxInjectionBlockList();


  //init the particle buffer
  PIC::ParticleBuffer::Init(2000000);
//  double TimeCounter=time(NULL);
  int LastDataOutputFileNumber=-1;



  //the total theoretical injection rate of hot oxigen
    if (_C_SPEC_>=0) {
        double rate=PIC::VolumeParticleInjection::GetTotalInjectionRate(_C_SPEC_);
        if (PIC::ThisThread==0) {
            printf("Total theoretical injection rate of hot %s: %e (%s@%i)\n",PIC::MolecularData::GetChemSymbol(_C_SPEC_),rate,__FILE__,__LINE__);
            printf("Integrated DR rate from Fox modes: %e\n",MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalO2PlusDissociativeRecombinationRate());
        }
    }
    if (_O_SPEC_>=0) {
        double rate=PIC::VolumeParticleInjection::GetTotalInjectionRate(_O_SPEC_);
        if (PIC::ThisThread==0) {
            printf("Total theoretical injection rate of hot %s: %e (%s@%i)\n",PIC::MolecularData::GetChemSymbol(_O_SPEC_),rate,__FILE__,__LINE__);
            printf("Integrated DR rate from Fox modes: %e\n",MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalO2PlusDissociativeRecombinationRate());
        }
    }



//============================== DEBUG =========================
  /* cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * */  node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];

  while (node!=NULL) {
    if (node->Temp_ID==2478) {
      cout << __FILE__<< "@" << __LINE__ << endl;


    }

    node=node->nextNodeThisThread;
  }

//============================== END DEBUG ====================================

  //determine the total number of the iterations to perform 
  //in the test-mode run 100 iterations and than output the particle data statistics
  int nIterations,nTotalIterations=100000001;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=100;  

  //time step
  for (long int niter=0;niter<nTotalIterations;niter++) {
     PIC::TimeStep();

//     PIC::MolecularCollisions::BackgroundAtmosphere::CollisionProcessor();

//  PIC::ChemicalReactions::PhotolyticReactions::Init();
//  PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
  
   if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh->ThisThread==0) cout << "The new lample length is " << PIC::RequiredSampleLength << endl;
     }
  }

 
  //output the particle statistics for the nightly tests 
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];
          sprintf(fname,"%s/test_Mars.dat",PIC::OutputDataFileDirectory);
          PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }


  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}






//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the function for particle's data sampling

#include "pic.h"
namespace PIC{
  namespace Mesh{
    //basic macroscopic parameters sampled in the simulation
    cDatumTimed    DatumParticleWeight(1,"\"Particle Weight\"",      false);
    cDatumTimed    DatumParticleNumber(1,"\"Particle Number\"",       true);
    cDatumTimed    DatumNumberDensity( 1,"\"Number Density[1/m^3]\"", true);

    cDatumWeighted DatumParticleVelocity(3, "\"Vx [m/s]\", \"Vy [m/s]\", \"Vz [m/s]\"", true);
    cDatumWeighted DatumParticleVelocity2(3,"\"Vx^2 [(m/s)^2]\", \"Vy^2 [(m/s)^2]\", \"Vz^2 [(m/s)^2]\"", false);

    cDatumWeighted DatumParticleSpeed(1,            "\"|V| [m/s]\"",     true);
    cDatumWeighted DatumParticleParallelVelocity(1, "\"Vpar [m/s]\"",   false);
    cDatumWeighted DatumParticleParallelVelocity2(1,"\"Vpar^2 [(m/s)^2]\"",false);

    //-------------------------------------------------------------------------
    // IMPORTANT: some data may be oncluded only for certain species!!!!
    //if(GetSpecieType(DataSetNumber)==_PIC_SPECIE_TYPE__GAS_)
    //-------------------------------------------------------------------------
    cDatumDerived DatumTranslationalTemperature(1, "\"Translational Temperature [K]\"", true);
    cDatumDerived DatumParallelTranslationalTemperature(1, "\"Parallel Translational Temperature [K]\"", true);
    cDatumDerived DatumTangentialTranslationalTemperature(1, "\"Tangential Translational Temperature [K]\"", true);

    // vector of active sampling data
    vector<cDatumSampled*> DataSampledCenterNodeActive;
    // vector of active derived data
    vector<cDatumDerived*> DataDerivedCenterNodeActive;
  }
}

//init the block's global data
int PIC::Mesh::cDataBlockAMR::LocalTimeStepOffset=0;
int PIC::Mesh::cDataBlockAMR::LocalParticleWeightOffset=0;
int PIC::Mesh::cDataBlockAMR::totalAssociatedDataLength=0;

//init the cells' global data
int PIC::Mesh::cDataCenterNode::totalAssociatedDataLength=0;
int PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset=0;

//the offsets to the sampled data stored in 'center nodes'
int PIC::Mesh::completedCellSampleDataPointerOffset=0,PIC::Mesh::collectingCellSampleDataPointerOffset=0;

int PIC::Mesh::sampledParticleWeghtRelativeOffset=0,PIC::Mesh::sampledParticleNumberRelativeOffset=0,PIC::Mesh::sampledParticleNumberDensityRelativeOffset=0;
int PIC::Mesh::sampledParticleVelocityRelativeOffset=0,PIC::Mesh::sampledParticleVelocity2RelativeOffset=0,PIC::Mesh::sampledParticleSpeedRelativeOffset=0;
int PIC::Mesh::sampledParticleNormalParallelVelocityRelativeOffset=0,PIC::Mesh::sampledParticleNormalParallelVelocity2RelativeOffset=0;
int PIC::Mesh::sampledExternalDataRelativeOffset=0;
int PIC::Mesh::sampleSetDataLength=0;

//the mesh parameters
double PIC::Mesh::xmin[3]={0.0,0.0,0.0},PIC::Mesh::xmax[3]={0.0,0.0,0.0};
PIC::Mesh::fLocalMeshResolution PIC::Mesh::LocalMeshResolution=NULL;
#if DIM == 3
cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR> PIC::Mesh::mesh;
#elif DIM == 2
cMeshAMR2d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR>  PIC::Mesh::mesh;
#else
cMeshAMR1d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR>  PIC::Mesh::mesh;
#endif


//the user defined functions for output of the 'ceneter node' data into a data file
vector<PIC::Mesh::fPrintVariableListCenterNode> PIC::Mesh::PrintVariableListCenterNode;
vector<PIC::Mesh::fPrintDataCenterNode> PIC::Mesh::PrintDataCenterNode;
vector<PIC::Mesh::fInterpolateCenterNode> PIC::Mesh::InterpolateCenterNode;



void PIC::Mesh::SetCellSamplingDataRequest() {
  exit(__LINE__,__FILE__,"not implemented yet");
}

void PIC::Mesh::cDataCenterNode::PrintVariableList(FILE* fout,int DataSetNumber) {
  // printe sampled data names
  vector<cDatumSampled*>::iterator ptrDatumSampled;
  for(ptrDatumSampled = DataSampledCenterNodeActive.begin();
      ptrDatumSampled!= DataSampledCenterNodeActive.end(); ptrDatumSampled++)
    if((*ptrDatumSampled)->doPrint) (*ptrDatumSampled)->PrintName(fout);
  // print derived data names
  vector<cDatumDerived*>::iterator ptrDatumDerived;
  for(ptrDatumDerived = DataDerivedCenterNodeActive.begin();
      ptrDatumDerived!= DataDerivedCenterNodeActive.end(); ptrDatumDerived++)
    if((*ptrDatumDerived)->doPrint) (*ptrDatumDerived)->PrintName(fout);
  
  //print the user defind 'center node' data
  vector<fPrintVariableListCenterNode>::iterator fptr;
  for (fptr=PrintVariableListCenterNode.begin();fptr!=PrintVariableListCenterNode.end();fptr++) (*fptr)(fout,DataSetNumber);
  
  //print varialbes sampled by the user defined sampling procedures
  if (PIC::IndividualModelSampling::PrintVariableList.size()!=0)
    for (unsigned int i=0;
	 i<PIC::IndividualModelSampling::PrintVariableList.size();i++) 
      PIC::IndividualModelSampling::PrintVariableList[i](fout,DataSetNumber);
}

void PIC::Mesh::cDataCenterNode::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread) {
  int idim;

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  static unsigned long int nCallCounter=0;
  ++nCallCounter;
#endif

  static int  nOutput=0;
  static bool IsFirstCall=true;
  static double* OutputData;

  // find size of message at the first call
  if(IsFirstCall){
    // include sampled data
    vector<cDatumSampled*>::iterator ptrDatumSampled;
    for(ptrDatumSampled = DataSampledCenterNodeActive.begin();
	ptrDatumSampled!= DataSampledCenterNodeActive.end(); ptrDatumSampled++)
      if((*ptrDatumSampled)->doPrint) nOutput+=(*ptrDatumSampled)->length;
    // include derived data
    vector<cDatumDerived*>::iterator ptrDatumDerived;
    for(ptrDatumDerived = DataDerivedCenterNodeActive.begin();
	ptrDatumDerived!= DataDerivedCenterNodeActive.end(); ptrDatumDerived++)
      if((*ptrDatumDerived)->doPrint) nOutput+=(*ptrDatumDerived)->length;
    OutputData = new double[nOutput];
    IsFirstCall=false;
  }
  
  if (pipe->ThisThread==CenterNodeThread) {
    // compose a message
    int iOutput=0;
    cDatumTimed* pDatumTimed;
    cDatumWeighted* pDatumWeighted;
    // sampled data values
    vector<cDatumSampled*>::iterator ptrDatumSampled;
    for(ptrDatumSampled = DataSampledCenterNodeActive.begin();
	ptrDatumSampled!= DataSampledCenterNodeActive.end(); ptrDatumSampled++)
      if((*ptrDatumSampled)->doPrint) {
	// cDatumSampled has 2 derived types: cDatumTimed & cDatumWeighted:
	// each is averaged differently => need to downcast them first
	if((*ptrDatumSampled)->type == PIC::Mesh::Timed_){
	  pDatumTimed = static_cast<cDatumTimed*> ((*ptrDatumSampled));
	  GetDatumAverage(*pDatumTimed,&OutputData[iOutput], DataSetNumber);
	}
	else{
	  pDatumWeighted = static_cast<cDatumWeighted*> ((*ptrDatumSampled));
	  GetDatumAverage(*pDatumWeighted,&OutputData[iOutput], DataSetNumber);
	}
	iOutput += (*ptrDatumSampled)->length;
      }
    //derived data values
    vector<cDatumDerived*>::iterator ptrDatumDerived;
    for(ptrDatumDerived = DataDerivedCenterNodeActive.begin();
	ptrDatumDerived!= DataDerivedCenterNodeActive.end(); ptrDatumDerived++)
      if((*ptrDatumDerived)->doPrint) {
	(this->*(*ptrDatumDerived)->GetValue)(&OutputData[iOutput],DataSetNumber);
	iOutput += (*ptrDatumDerived)->length;
      }
  }

  
  if (pipe->ThisThread==0) {
    //print values to the output file
    if (CenterNodeThread!=0)pipe->recv((char*)OutputData,nOutput*sizeof(double),CenterNodeThread);
    for(int iOutput=0; iOutput<nOutput; iOutput++)
      fprintf(fout, "%e ", OutputData[iOutput]);
  }
  else pipe->send((char*)OutputData,nOutput*sizeof(double));

  //print the user defind 'center node' data
  vector<fPrintDataCenterNode>::iterator fptr;

  for (fptr=PrintDataCenterNode.begin();fptr!=PrintDataCenterNode.end();fptr++) (*fptr)(fout,DataSetNumber,pipe,CenterNodeThread,this);

  //print data sampled by the user defined sampling functions
  if (PIC::IndividualModelSampling::PrintSampledData.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::PrintSampledData.size();i++) PIC::IndividualModelSampling::PrintSampledData[i](fout,DataSetNumber,pipe,CenterNodeThread,this);
  }

}

void PIC::Mesh::cDataCenterNode::Interpolate(cDataCenterNode** InterpolationList,double *InterpolationCoefficients,int nInterpolationCoefficients) {
  int i,s,idim;
  double c;



  //==============================  DEBUGGER ===============
           if (nInterpolationCoefficients!=0) Measure=InterpolationList[0]->Measure;

           static long int nCallCounter=0;
           nCallCounter++;

  //============================== END DEBUGGER ============


  #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
  if (associatedDataPointer==NULL) exit(__LINE__,__FILE__,"Error: The associated data buffer is not initialized");
  #endif

  double InterpolatedParticleWeight=0.0,InterpolatedParticleNumber=0.0,InterpolatedParticleNumberDeinsity=0.0,InterpolatedBulkVelocity[3]={0.0,0.0,0.0},InterpolatedBulk2Velocity[3]={0.0,0.0,0.0};
  double InterpolatedParticleSpeed=0.0;
  double pWeight;

#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
#else
  double InterpolatedBulkParallelVelocity,InterpolatedBulk2ParallelVelocity;
#endif

  for (s=0;s<PIC::nTotalSpecies;s++) {
    InterpolatedParticleWeight=0.0,InterpolatedParticleNumber=0.0,InterpolatedParticleNumberDeinsity=0.0,InterpolatedParticleSpeed=0.0;
    for (idim=0;idim<3;idim++) InterpolatedBulkVelocity[idim]=0.0,InterpolatedBulk2Velocity[idim]=0.0;

#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
#else
    InterpolatedBulkParallelVelocity=0.0,InterpolatedBulk2ParallelVelocity=0.0;
#endif

    // interpolate all sampled data
    vector<cDatumSampled*>::iterator ptrDatum;
    for(ptrDatum = DataSampledCenterNodeActive.begin();
    	ptrDatum!= DataSampledCenterNodeActive.end(); ptrDatum++)
      InterpolateDatum(**ptrDatum,
    		       InterpolationList,
    		       InterpolationCoefficients,
    		       nInterpolationCoefficients, s);
  }

  //print the user defind 'center node' data
  vector<fInterpolateCenterNode>::iterator fptr;

  for (fptr=InterpolateCenterNode.begin();fptr!=InterpolateCenterNode.end();fptr++) (*fptr)(InterpolationList,InterpolationCoefficients,nInterpolationCoefficients,this);

  //interpolate data sampled by user defiend sampling procedures
  if (PIC::IndividualModelSampling::InterpolateCenterNodeData.size()!=0) {
    for (unsigned int ifunc=0;ifunc<PIC::IndividualModelSampling::PrintVariableList.size();ifunc++) PIC::IndividualModelSampling::InterpolateCenterNodeData[ifunc](InterpolationList,InterpolationCoefficients,nInterpolationCoefficients,this);
  }
}

void PIC::Mesh::initCellSamplingDataBuffer() {

  if (cDataBlockAMR::totalAssociatedDataLength!=0) exit(__LINE__,__FILE__,"Error: reinitialization of the blocks associated data offsets");

  //local time step
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  if (PIC::ThisThread==0) cout << "$PREFIX:Time step mode: specie dependent local time step" << endl;
  cDataBlockAMR::LocalTimeStepOffset=0;
  cDataBlockAMR::totalAssociatedDataLength+=sizeof(double)*PIC::nTotalSpecies;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  //do nothing for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  #elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  //do nothing for _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif

  //local particle weight
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  if (PIC::ThisThread==0) cout << "$PREFIX:Particle weight mode: specie dependent local weight" << endl;
  cDataBlockAMR::LocalParticleWeightOffset=cDataBlockAMR::totalAssociatedDataLength;
  cDataBlockAMR::totalAssociatedDataLength+=sizeof(double)*PIC::nTotalSpecies;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ ==_SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  //do nothing for _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ ==_SINGLE_GLOBAL_PARTICLE_WEIGHT_
  //do nothing for _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SINGLE_GLOBAL_PARTICLE_WEIGHT_
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif

  //set up the offsets for 'center node' sampled data
  long int offset=0;
  DatumParticleWeight.activate(offset, &DataSampledCenterNodeActive);
  DatumNumberDensity.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleNumber.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleVelocity.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleVelocity2.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleSpeed.activate(offset, &DataSampledCenterNodeActive);
  DatumTranslationalTemperature.activate(&cDataCenterNode::GetTranslationalTemperature, &DataDerivedCenterNodeActive);

  //sampling the 'parallel' and 'tangential' kinetic temperatures
#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
  //do nothing
#else
  DatumParticleParallelVelocity.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleParallelVelocity2.activate(offset,&DataSampledCenterNodeActive);
  DatumParallelTranslationalTemperature.activate(&cDataCenterNode::GetParallelTranslationalTemperature, &DataDerivedCenterNodeActive);
  DatumTangentialTranslationalTemperature.activate(&cDataCenterNode::GetTangentialTranslationalTemperature, &DataDerivedCenterNodeActive);
#endif

  //check if user defined sampling data is requested
  sampledExternalDataRelativeOffset=offset;

  if (PIC::IndividualModelSampling::DataSampledList.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::DataSampledList.size();i++) PIC::IndividualModelSampling::DataSampledList[i]->activate(offset, &DataSampledCenterNodeActive);
  }

  if (PIC::IndividualModelSampling::RequestSamplingData.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::RequestSamplingData.size();i++) offset+=PIC::IndividualModelSampling::RequestSamplingData[i](offset);
  }




  PIC::Mesh::sampleSetDataLength=offset;

  PIC::Mesh::completedCellSampleDataPointerOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::collectingCellSampleDataPointerOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+PIC::Mesh::sampleSetDataLength;

  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=2*PIC::Mesh::sampleSetDataLength;

  //the volume partilce injection: save the volume particle injection rate
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

  //allocate the model requested static (not sampling) cell data
  if (PIC::IndividualModelSampling::RequestStaticCellData.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::RequestStaticCellData.size();i++) {
      PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=PIC::IndividualModelSampling::RequestStaticCellData[i](PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
    }
  }
}


//flush and switch the sampling buffers in 'center' nodes
void PIC::Mesh::flushCompletedSamplingBuffer(cDataCenterNode* node) {
  register int i,length=PIC::Mesh::sampleSetDataLength/sizeof(double);
  register double *ptr;

  for (i=0,ptr=(double*)(node->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset);i<length;i++,ptr++) *ptr=0.0;
}

void PIC::Mesh::flushCollectingSamplingBuffer(cDataCenterNode* node) {
  register int i,length=PIC::Mesh::sampleSetDataLength/sizeof(double);
  register double *ptr;

  for (i=0,ptr=(double*)(node->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset);i<length;i++,ptr++) *ptr=0.0;
}

void PIC::Mesh::switchSamplingBuffers() {
  int tempOffset;

  //exchange the CellSampleData offsets
  tempOffset=PIC::Mesh::completedCellSampleDataPointerOffset;
  PIC::Mesh::completedCellSampleDataPointerOffset=PIC::Mesh::collectingCellSampleDataPointerOffset;
  PIC::Mesh::collectingCellSampleDataPointerOffset=tempOffset;

  //switch the offsets for the internal spherical surfaces installed into the mesh
  PIC::BC::InternalBoundary::Sphere::switchSamplingBuffers();

}

//==============================================================
//init set and set up the computational mesh
void PIC::Mesh::Init(double* xMin,double* xMax,fLocalMeshResolution ResolutionFunction) {

  for (int idim=0;idim<DIM;idim++) xmin[idim]=xMin[idim],xmax[idim]=xMax[idim];

  LocalMeshResolution=ResolutionFunction;
  mesh.init(xMin,xMax,LocalMeshResolution);
}

void PIC::Mesh::buildMesh() {
  mesh.buildMesh();
}



void PIC::Mesh::cDataBlockAMR::sendBoundaryLayerBlockData(CMPI_channel *pipe) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=GetCenterNode(LocalCellNumber);

    pipe->send(cell->associatedDataPointer,cell->totalAssociatedDataLength);
    pipe->send(cell->Measure);
  }

  pipe->send(associatedDataPointer,totalAssociatedDataLength);
}

void PIC::Mesh::cDataBlockAMR::sendMoveBlockAnotherProcessor(CMPI_channel *pipe) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
//  PIC::Mesh::cDataCenterNode *cell=NULL;

  sendBoundaryLayerBlockData(pipe);

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  //send all blocks' data when the blocks is moved to another processor
  long int Particle,NextParticle;
  char *buffer=new char[PIC::ParticleBuffer::ParticleDataLength];

  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
//    LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
//    cell=GetCenterNode(LocalCellNumber);

    //    Particle=cell->FirstCellParticle;
    Particle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

    if  (Particle!=-1) {
      LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
      pipe->send(_CENTRAL_NODE_NUMBER_SIGNAL_);
      pipe->send(LocalCellNumber);

      while (Particle!=-1) {
        PIC::ParticleBuffer::PackParticleData(buffer,Particle);
        pipe->send(_NEW_PARTICLE_SIGNAL_);
        pipe->send(buffer,PIC::ParticleBuffer::ParticleDataLength);

        NextParticle=PIC::ParticleBuffer::GetNext(Particle);
//        PIC::ParticleBuffer::DeleteParticle(Particle);
        PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(Particle);

        Particle=NextParticle;
      }

     // cell->FirstCellParticle=-1;
      FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=-1;
    }
  }

  pipe->send(_END_COMMUNICATION_SIGNAL_);
  delete [] buffer;
}

void PIC::Mesh::cDataBlockAMR::recvBoundaryLayerBlockData(CMPI_channel *pipe,int From) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=GetCenterNode(LocalCellNumber);

    pipe->recv(cell->associatedDataPointer,cell->totalAssociatedDataLength,From);
    pipe->recv(cell->Measure,From);
  }

  pipe->recv(associatedDataPointer,totalAssociatedDataLength,From);
}

//recieve all blocks' data when the blocks is moved to another processo
void PIC::Mesh::cDataBlockAMR::recvMoveBlockAnotherProcessor(CMPI_channel *pipe,int From) {
  long int LocalCellNumber=-1;
//  PIC::Mesh::cDataCenterNode *cell=NULL;
  int i=-10,j=-10,k=-10;

  recvBoundaryLayerBlockData(pipe,From);

  long int Particle;
  char buffer[PIC::ParticleBuffer::ParticleDataLength];

//  char *buffer=new char[PIC::ParticleBuffer::ParticleDataLength];

  int Signal;
  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;


  pipe->recv(Signal,From);
  //LocalCellNumber=-1,cell=NULL;

  while (Signal!=_END_COMMUNICATION_SIGNAL_) {
    switch (Signal) {
    case _CENTRAL_NODE_NUMBER_SIGNAL_ :
      pipe->recv(LocalCellNumber,From);
//      cell=GetCenterNode(LocalCellNumber);

      PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(LocalCellNumber,i,j,k);
      break;
    case _NEW_PARTICLE_SIGNAL_ :
      pipe->recv(buffer,PIC::ParticleBuffer::ParticleDataLength,From);

//      Particle=PIC::ParticleBuffer::GetNewParticle(cell->FirstCellParticle);

      Particle=PIC::ParticleBuffer::GetNewParticle(FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
      PIC::ParticleBuffer::UnPackParticleData(buffer,Particle);
      break;
    default :
      exit(__LINE__,__FILE__,"Error: unknown option");
    }

    pipe->recv(Signal,From);
  }

//  delete [] buffer;
}


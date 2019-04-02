
//$Id$
//the field solver routines

/*
 * pic_field_solver_ecsim.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: vtenishe
 */

#include "pic.h"
//#include "LinearSystemCornerNode.h"
//#include "linear_solver_wrapper_c.h"

//using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

PIC::FieldSolver::Electromagnetic::ECSIM::fSetIC PIC::FieldSolver::Electromagnetic::ECSIM::SetIC=PIC::FieldSolver::Electromagnetic::ECSIM::SetIC_default;

PIC::FieldSolver::Electromagnetic::ECSIM::fUserDefinedParticleBC  PIC::FieldSolver::Electromagnetic::ECSIM::setParticle_BC;
PIC::FieldSolver::Electromagnetic::ECSIM::fUserDefinedFieldBC PIC::FieldSolver::Electromagnetic::ECSIM::setE_half_BC, 
  PIC::FieldSolver::Electromagnetic::ECSIM::setE_curr_BC, PIC::FieldSolver::Electromagnetic::ECSIM::setB_center_BC,
  PIC::FieldSolver::Electromagnetic::ECSIM::setB_corner_BC;

cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1> PIC::FieldSolver::Electromagnetic::ECSIM::Solver;
cLinearSystemCenterNode<PIC::Mesh::cDataCenterNode,1,7,0,1,1,0>  PIC::FieldSolver::Electromagnetic::ECSIM::PoissonSolver;

bool PIC::FieldSolver::Electromagnetic::ECSIM::DoDivECorrection = false;
int PIC::FieldSolver::Electromagnetic::ECSIM::CurrentEOffset=-1;
int PIC::FieldSolver::Electromagnetic::ECSIM::OffsetE_HalfTimeStep=-1;
int PIC::FieldSolver::Electromagnetic::ECSIM::CurrentBOffset=-1;
int PIC::FieldSolver::Electromagnetic::ECSIM::PrevBOffset=-1;
int PIC::FieldSolver::Electromagnetic::ECSIM::OffsetB_corner;

int PIC::FieldSolver::Electromagnetic::ECSIM::ExOffsetIndex=0;
int PIC::FieldSolver::Electromagnetic::ECSIM::EyOffsetIndex=1;
int PIC::FieldSolver::Electromagnetic::ECSIM::EzOffsetIndex=2;
int PIC::FieldSolver::Electromagnetic::ECSIM::JxOffsetIndex;
int PIC::FieldSolver::Electromagnetic::ECSIM::JyOffsetIndex;
int PIC::FieldSolver::Electromagnetic::ECSIM::JzOffsetIndex;
int PIC::FieldSolver::Electromagnetic::ECSIM::BxOffsetIndex=0;
int PIC::FieldSolver::Electromagnetic::ECSIM::ByOffsetIndex=1;
int PIC::FieldSolver::Electromagnetic::ECSIM::BzOffsetIndex=2;
int PIC::FieldSolver::Electromagnetic::ECSIM::MassMatrixOffsetIndex;

int PIC::FieldSolver::Electromagnetic::ECSIM::netChargeOldIndex;
int PIC::FieldSolver::Electromagnetic::ECSIM::netChargeNewIndex;
int PIC::FieldSolver::Electromagnetic::ECSIM::divEIndex;
int PIC::FieldSolver::Electromagnetic::ECSIM::phiIndex;

int PIC::FieldSolver::Electromagnetic::ECSIM::Rho_=0;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUx_=1;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUy_=2;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUz_=3;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUxUx_=4;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUyUy_=5;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUzUz_=6;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUxUy_=7;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUyUz_=8;
int PIC::FieldSolver::Electromagnetic::ECSIM::RhoUxUz_=9;
int * PIC::FieldSolver::Electromagnetic::ECSIM::SpeciesDataIndex=NULL;


double PIC::FieldSolver::Electromagnetic::ECSIM::CumulativeTiming::SolveTime=0.0;
double PIC::FieldSolver::Electromagnetic::ECSIM::CumulativeTiming::UpdateBTime=0.0;
double PIC::FieldSolver::Electromagnetic::ECSIM::CumulativeTiming::UpdateETime=0.0;
double PIC::FieldSolver::Electromagnetic::ECSIM::CumulativeTiming::UpdateJMassMatrixTime=0.0;
double PIC::FieldSolver::Electromagnetic::ECSIM::CumulativeTiming::TotalRunTime=0.0;

void PIC::FieldSolver::Electromagnetic::ECSIM::CumulativeTiming::Print() {
  if (PIC::ThisThread==0) {
    printf("\n$PREFIX: Electromagnetic::ECSIM timing:\n");
    printf("$PREFIX: SolveTime=%e\n",SolveTime);
    printf("$PREFIX: UpdateBTime=%e\n",UpdateBTime);
    printf("$PREFIX: UpdateETime=%e\n",UpdateETime);
    printf("$PREFIX: UpdateJMassMatrixTime=%e\n",UpdateJMassMatrixTime);
    printf("$PREFIX: TotalRunTime=%e\n",TotalRunTime);
  }
}



//location of the solver's data in the corner node associated data vector
int PIC::FieldSolver::Electromagnetic::ECSIM::CornerNodeAssociatedDataOffsetBegin=-1,PIC::FieldSolver::Electromagnetic::ECSIM::CornerNodeAssociatedDataOffsetLast=-1;

double dtTotal = 0.0;
double PIC::FieldSolver::Electromagnetic::ECSIM::cDt=0.0;
double PIC::FieldSolver::Electromagnetic::ECSIM::theta =0.5;
double epsilon0=8.85418782e-12;
double mu0=1.25663706e-6;

#if _PIC_FIELD_SOLVER_INPUT_UNIT_== _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_  
double PIC::FieldSolver::Electromagnetic::ECSIM::LightSpeed =1;
#elif _PIC_FIELD_SOLVER_INPUT_UNIT_== _PIC_FIELD_SOLVER_INPUT_UNIT_SI_
double PIC::FieldSolver::Electromagnetic::ECSIM::LightSpeed =1/sqrt(epsilon0*mu0)*1e2;//in cm/s
#endif

double TotalParticleEnergy=0.0;
double TotalWaveEnergy=0.0;

#if _PIC_FIELD_SOLVER_INPUT_UNIT_== _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_
double E_conv = 1;
double B_conv = 1;
double mass_conv =1.0/_AMU_;
double charge_conv=1.0/ElectronCharge;
double length_conv=1;
#elif _PIC_FIELD_SOLVER_INPUT_UNIT_== _PIC_FIELD_SOLVER_INPUT_UNIT_SI_
double E_conv = 1e6/PIC::FieldSolver::Electromagnetic::ECSIM::LightSpeed; //convert from SI to cgs
double B_conv = 1e4;
double mass_conv = 1e3;
double charge_conv=0.1*PIC::FieldSolver::Electromagnetic::ECSIM::LightSpeed;
double length_conv=1e2; 
#endif


inline double interp2D(double vmm, double vpm, double vpp, double vmp,
                                               double dx, double dy) {

  /*Interpolate the value at p from the surounding 4 points. 'p' means plus, and
    'm' means minus.
        vmp----------------vpp
        |                   |
        |                   |
        |----dx-----p       |
        |           |       |
        |           dy      |
        |           |       |
        vmm----------------vpm
   */

    return vmm * (1 - dx) * (1 - dy) + vpm * dx * (1 - dy) + vpp * dx * dy +
           vmp * (1 - dx) * dy;
}



//magnetic field
int PackBlockData_B(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
  int dataLengthByte=6*sizeof(double);

  return PIC::Mesh::PackBlockData_Internal(NodeTable,NodeTableLength,NodeDataLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      SendDataBuffer,
      NULL,NULL,0,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0);
}



//magnetic fieled
int UnpackBlockData_B(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
  int dataLengthByte=6*sizeof(double);

  return PIC::Mesh::UnpackBlockData_Internal(NodeTable,NodeTableLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      RecvDataBuffer,
      NULL,NULL,0,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0);
}

//net charge
int PackBlockData_netCharge(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer) {
    int ibegin=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset
        +PIC::FieldSolver::Electromagnetic::ECSIM::netChargeNewIndex*sizeof(double);
    int dataLengthByte=1*sizeof(double);

    return PIC::Mesh::PackBlockData_Internal(NodeTable,NodeTableLength,NodeDataLength,
                                                    BlockCenterNodeMask,BlockCornerNodeMask,
                                                    SendDataBuffer,
                                                    NULL,NULL,0,
                                                    &ibegin,&dataLengthByte,1,
                                                    NULL,NULL,0);
}

//net charge
int UnpackBlockData_netCharge(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset
    +PIC::FieldSolver::Electromagnetic::ECSIM::netChargeNewIndex*sizeof(double);
  int dataLengthByte=1*sizeof(double);
  
  return PIC::Mesh::UnpackBlockData_Internal(NodeTable,NodeTableLength,
                                             BlockCenterNodeMask,BlockCornerNodeMask,
                                             RecvDataBuffer,
                                             NULL,NULL,0,
                                             &ibegin,&dataLengthByte,1,
                                             NULL,NULL,0);
}

//phi
int PackBlockData_phi(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset
    +PIC::FieldSolver::Electromagnetic::ECSIM::phiIndex*sizeof(double);
  int dataLengthByte=1*sizeof(double);
  
  return PIC::Mesh::PackBlockData_Internal(NodeTable,NodeTableLength,NodeDataLength,
                                           BlockCenterNodeMask,BlockCornerNodeMask,
                                           SendDataBuffer,
                                           NULL,NULL,0,
                                           &ibegin,&dataLengthByte,1,
                                           NULL,NULL,0);
}


//phi
int UnpackBlockData_phi(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset
    +PIC::FieldSolver::Electromagnetic::ECSIM::phiIndex*sizeof(double);
  int dataLengthByte=1*sizeof(double);
  
  return PIC::Mesh::UnpackBlockData_Internal(NodeTable,NodeTableLength,
                                             BlockCenterNodeMask,BlockCornerNodeMask,
                                             RecvDataBuffer,
                                             NULL,NULL,0,
                                             &ibegin,&dataLengthByte,1,
                                             NULL,NULL,0);
}


//electric field
int PackBlockData_E(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  int dataLengthByte=6*sizeof(double);

  return PIC::Mesh::PackBlockData_Internal(NodeTable,NodeTableLength,NodeDataLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      SendDataBuffer,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0,
      NULL,NULL,0);
}


//electric fieled
int UnpackBlockData_E(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  int dataLengthByte=6*sizeof(double);

  return PIC::Mesh::UnpackBlockData_Internal(NodeTable,NodeTableLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      RecvDataBuffer,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0,
      NULL,NULL,0);
}

//current and massmatrix
int PackBlockData_JMassMatrix(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+6*sizeof(double);
  int dataLengthByte=246*sizeof(double);

  return PIC::Mesh::PackBlockData_Internal(NodeTable,NodeTableLength,NodeDataLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      SendDataBuffer,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0,
      NULL,NULL,0);
}


//current and massmatrix
int UnpackBlockData_JMassMatrix(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+6*sizeof(double);
  int dataLengthByte=246*sizeof(double);

  return PIC::Mesh::UnpackBlockData_Internal(NodeTable,NodeTableLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      RecvDataBuffer,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0,
      NULL,NULL,0);
}

int PackBlockData_JMassMatrixSpeciesData(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+6*sizeof(double);
  int dataLengthByte=(246+10*PIC::nTotalSpecies)*sizeof(double);

  return PIC::Mesh::PackBlockData_Internal(NodeTable,NodeTableLength,NodeDataLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      SendDataBuffer,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0,
      NULL,NULL,0);
}

int UnpackBlockData_JMassMatrixSpeciesData(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer) {
  int ibegin=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+6*sizeof(double);
  int dataLengthByte=(246+10*PIC::nTotalSpecies)*sizeof(double);

  return PIC::Mesh::UnpackBlockData_Internal(NodeTable,NodeTableLength,
      BlockCenterNodeMask,BlockCornerNodeMask,
      RecvDataBuffer,
      &ibegin,&dataLengthByte,1,
      NULL,NULL,0,
      NULL,NULL,0);
}



bool IsInit=false;

//The global initialization procedure
void PIC::FieldSolver::Init() {
  if (IsInit) exit(__LINE__,__FILE__,"Error: The field solver already initialized");

  switch (_PIC_FIELD_SOLVER_MODE_) {
  case _PIC_FIELD_SOLVER_MODE__ELECTROMAGNETIC__ECSIM_:
    PIC::FieldSolver::Electromagnetic::ECSIM::Init();
    break;
  default:
    exit(__LINE__,__FILE__,"Error: The field solver Init() has been called with an unknown _PIC_FIELD_SOLVER_MODE_ value");
  }

  IsInit = true;
}

//Field Solver of IPIC3D
void PIC::FieldSolver::Electromagnetic::ECSIM::Init() {
  //init the electric and magnetic field offsets
  //Magnetic field is in the center nodes
  //Electric field is in the corner nodes
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    

  CornerNodeAssociatedDataOffsetBegin=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

  //register the linear system solver with the core 
  PIC::RegisterLinearSolver(&Solver); 

  //set function for timing of the field solver
  PIC::RunTimeSystemState::CumulativeTiming::PrintTimingFunctionTable.push_back(CumulativeTiming::Print);

  if (PIC::CPLR::DATAFILE::Offset::MagneticField.active==true) {
    exit(__LINE__,__FILE__,"Error: reinitialization of the magnetic field offset");
  }
  else {
    PIC::CPLR::DATAFILE::Offset::MagneticField.active=true;
    PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;   
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=2*PIC::CPLR::DATAFILE::Offset::MagneticField.nVars*sizeof(double);
    CurrentBOffset =0;
    PrevBOffset = 3*sizeof(double);
  }

  netChargeOldIndex = 6;
  netChargeNewIndex = 7;
  divEIndex = 8;
  phiIndex = 9;
  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=4*sizeof(double);

  if (PIC::CPLR::DATAFILE::Offset::ElectricField.active==true) {
    exit(__LINE__,__FILE__,"Error: reinitialization of the electric field offset");
  }
  else {
    PIC::CPLR::DATAFILE::Offset::ElectricField.active=true;
    PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset=PIC::Mesh::cDataCornerNode::totalAssociatedDataLength;
    CurrentEOffset=0;
    PIC::Mesh::cDataCornerNode::totalAssociatedDataLength+=2*PIC::CPLR::DATAFILE::Offset::ElectricField.nVars*sizeof(double);
    OffsetE_HalfTimeStep=3*sizeof(double);
  }

  //allocate memory for Jx,Jy,Jz
  PIC::Mesh::cDataCornerNode::totalAssociatedDataLength+=3*sizeof(double);
  //J starts from PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset
  JxOffsetIndex = 6;
  JyOffsetIndex = 7;
  JzOffsetIndex = 8;
  
  //allocate memory for 81 mass matrix elements
  PIC::Mesh::cDataCornerNode::totalAssociatedDataLength+=243*sizeof(double);
  MassMatrixOffsetIndex = 9;

#if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_  
  PIC::Mesh::cDataCornerNode::totalAssociatedDataLength+=10*sizeof(double)*PIC::nTotalSpecies;
  SpeciesDataIndex = new int[PIC::nTotalSpecies];
  for (int iSp=0; iSp<PIC::nTotalSpecies; iSp++) 
    SpeciesDataIndex[iSp]=9+243+iSp*10; // 10 vars for each species
#endif
  
#if _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CORNER_BASED_
  OffsetB_corner = PIC::Mesh::cDataCornerNode::totalAssociatedDataLength-
    PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  PIC::Mesh::cDataCornerNode::totalAssociatedDataLength += 6*sizeof(double);
#endif

  PIC::Mesh::mesh.GetCenterNodesInterpolationCoefficients=PIC::Mesh::GetCenterNodesInterpolationCoefficients;
     
  PIC::Mesh::AddVaraibleListFunction(PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCenterNodeVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCenterNodeData);
  PIC::Mesh::InterpolateCenterNode.push_back(PIC::FieldSolver::Electromagnetic::ECSIM::output::InterpolateCenterNode);

  PIC::Mesh::PrintVariableListCornerNode.push_back(PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCornerNodeVariableList);
  PIC::Mesh::PrintDataCornerNode.push_back(PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCornerNodeData);
  
  //  PIC::FieldSolver::Electromagnetic::ECSIM::Init_IC(); 

  CornerNodeAssociatedDataOffsetLast=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength-1; //CornerNodeAssociatedDataOffsetLast still belongs to the solver
}

void PIC::FieldSolver::Electromagnetic::ECSIM::Init_IC() {
  //set the initial conditions
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  if (_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_)
    dtTotal=PIC::CPLR::FLUID::dt; 
  //global time step with deallocated blocks in fluid coupler mode will have some issue
  else 
    dtTotal=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
  PIC::FieldSolver::Electromagnetic::ECSIM::cDt=LightSpeed*dtTotal;
  
  PIC::FieldSolver::Electromagnetic::ECSIM::BuildMatrix();
  SetIC();
}

//set default initial conditions
void  PIC::FieldSolver::Electromagnetic::ECSIM::SetIC_default() {
  int i,j,k,iNode,idim;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *CenterNode;
  PIC::Mesh::cDataCornerNode *CornerNode;
  PIC::Mesh::cDataBlockAMR *block;
  double *E,*B;
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    
  if (PIC::CPLR::DATAFILE::Offset::ElectricField.active==false) exit(__LINE__,__FILE__,"Error: the electric field offset is not active");
  if (PIC::CPLR::DATAFILE::Offset::MagneticField.active==false) exit(__LINE__,__FILE__,"Error: the magnetic field offset is not active");

  //loop through all blocks
  for (int iNode=0;iNode<DomainBlockDecomposition::nLocalBlocks;iNode++) {
    node=DomainBlockDecomposition::BlockTable[iNode];
    block=node->block;

    if (block!=NULL) {      
      //set the electric field (corner nodes)
      // the loop index is changed
      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));

        if (CornerNode!=NULL) {
          E=(double*)(CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+CurrentEOffset);
          for (idim=0;idim<3;idim++) E[idim]=0.0;

          E=(double*)(CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+OffsetE_HalfTimeStep);
          for (idim=0;idim<3;idim++) E[idim]=0.0;
        }
      }

      //set the magnetic field (center nodes)
      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        CenterNode=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

        if (CenterNode!=NULL) {
          B=(double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset);
          for (idim=0;idim<3;idim++) B[idim]=0.0;
        }
      }
    }
  }

  //update the 'ghost' cells and 'ghost' blocks
  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData();
    break;
  default:
    PIC::Mesh::mesh.ParallelBlockDataExchange();
  }
}

int MassMatrixOffsetTable[3][81];
bool initMassMatrixOffsetTable=false;

void computeMassMatrixOffsetTable(){

 int indexAddition[3] = {0,-1,1};
  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){
    for (int ii=0;ii<3;ii++){
      for (int jj=0;jj<3;jj++){
        for (int kk=0;kk<3;kk++){
          int iElement = iVarIndex*27+ii+jj*3+kk*9;

          MassMatrixOffsetTable[0][iElement]=(ii+jj*3+kk*9)*9+0*3+iVarIndex;
          MassMatrixOffsetTable[1][iElement]=(ii+jj*3+kk*9)*9+1*3+iVarIndex;
          MassMatrixOffsetTable[2][iElement]=(ii+jj*3+kk*9)*9+2*3+iVarIndex;
        }
      }
    }
  }
  
  initMassMatrixOffsetTable=true;
}

static const double LaplacianCoeff[3][27]={{-0.5,0.25,0.25,-0.25,0.125,0.125,-0.25,0.125,0.125,
			       -0.25,0.125,0.125,-0.125,0.0625,0.0625,-0.125,0.0625,0.0625,
			       -0.25,0.125,0.125,-0.125,0.0625,0.0625,-0.125,0.0625,0.0625},
			      {-0.5,-0.25,-0.25,0.25,0.125,0.125,0.25,0.125,0.125,
			       -0.25,-0.125,-0.125,0.125,0.0625,0.0625,0.125,0.0625,0.0625,
			       -0.25,-0.125,-0.125,0.125,0.0625,0.0625,0.125,0.0625,0.0625},
			      {-0.5,-0.25,-0.25,-0.25,-0.125,-0.125,-0.25,-0.125,-0.125,
			       0.25,0.125,0.125,0.125,0.0625,0.0625,0.125,0.0625,0.0625,
			       0.25,0.125,0.125,0.125,0.0625,0.0625,0.125,0.0625,0.0625}};

static const double graddiv[3][3][27]={{{-0.5,0.25,0.25,-0.25,0.125,0.125,-0.25,0.125,0.125,
				-0.25,0.125,0.125,-0.125,0.0625,0.0625,-0.125,0.0625,0.0625,
				-0.25,0.125,0.125,-0.125,0.0625,0.0625,-0.125,0.0625,0.0625},
			       {0,0,0,0,0.125,-0.125,0,-0.125,0.125,0,0,0,0,0.0625,-0.0625,0,
				-0.0625,0.0625,0,0,0,0,0.0625,-0.0625,0,-0.0625,0.0625},
			       {0,0,0,0,0,0,0,0,0,0,0.125,-0.125,0,0.0625,-0.0625,
				0,0.0625,-0.0625,0,-0.125,0.125,0,-0.0625,0.0625,0,-0.0625,0.0625}},
			      {{0,0,0,0,0.125,-0.125,0,-0.125,0.125,0,0,0,0,0.0625,-0.0625,0,-0.0625,0.0625,
				0,0,0,0,0.0625,-0.0625,0,-0.0625,0.0625},
			       {-0.5,-0.25,-0.25,0.25,0.125,0.125,0.25,0.125,0.125,-0.25,-0.125,-0.125,
				0.125,0.0625,0.0625,0.125,0.0625,0.0625,-0.25,-0.125,-0.125,0.125,0.0625,
				0.0625,0.125,0.0625,0.0625},
			       {0,0,0,0,0,0,0,0,0,0,0,0,0.125,0.0625,0.0625,-0.125,-0.0625,-0.0625,0,0,0,
				-0.125,-0.0625,-0.0625,0.125,0.0625,0.0625}},
			      {{0,0,0,0,0,0,0,0,0,0,0.125,-0.125,0,0.0625,-0.0625,0,0.0625,-0.0625,0,
				-0.125,0.125,0,-0.0625,0.0625,0,-0.0625,0.0625},
			       {0,0,0,0,0,0,0,0,0,0,0,0,0.125,0.0625,0.0625,-0.125,-0.0625,-0.0625,
				0,0,0,-0.125,-0.0625,-0.0625,0.125,0.0625,0.0625},
			       {-0.5,-0.25,-0.25,-0.25,-0.125,-0.125,-0.25,-0.125,-0.125,0.25,0.125,
				0.125,0.125,0.0625,0.0625,0.125,0.0625,0.0625,0.25,0.125,0.125,0.125,
				0.0625,0.0625,0.125,0.0625,0.0625}}};



void PIC::FieldSolver::Electromagnetic::ECSIM::PoissonGetStencil(int i, int j, int k, int iVar,
                       cLinearSystemCenterNode<PIC::Mesh::cDataCenterNode,1,7,0,1,1,0>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,cLinearSystemCenterNode<PIC::Mesh::cDataCenterNode,1,7,0,1,1,0>::cRhsSupportTable* RhsSupportTable_CornerNodes,int& RhsSupportLength_CornerNodes,cLinearSystemCenterNode<PIC::Mesh::cDataCenterNode,1,7,0,1,1,0>::cRhsSupportTable* RhsSupportTable_CenterNodes,int& RhsSupportLength_CenterNodes, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node){
 
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  double x[3];
  //power of 3 array created
  //for some pgi compilers that cannot convert result of pow() from double to int correctly
  int pow3_arr[3]={1,3,9};
  int index[3] = {i,j,k};
  int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
  double dx[3],dx_no[3]; 
  double CoeffSqr[3];

  for (int iDim=0; iDim<3; iDim++){
    dx_no[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim];
    //convert length
    dx[iDim] =dx_no[iDim]* length_conv;
    CoeffSqr[iDim] = 1./(dx[iDim]*dx[iDim]);
    x[iDim]=node->xmin[iDim]+(index[iDim]+0.5)*dx_no[iDim];
  }

  if (isBoundaryCell(x,dx_no,node)) {  
    MatrixRowNonZeroElementTable[0].i=i;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[0].iVar=iVar;
    MatrixRowNonZeroElementTable[0].MatrixElementParameterTable[0]=0.0;
    MatrixRowNonZeroElementTable[0].MatrixElementParameterTableLength=1;
    MatrixRowNonZeroElementTable[0].MatrixElementSupportTableLength = 0;
    MatrixRowNonZeroElementTable[0].Node=node;
    MatrixRowNonZeroElementTable[0].MatrixElementSupportTable[0]=NULL;
    NonZeroElementsFound =1;
    rhs =0.0;
    RhsSupportLength_CenterNodes = 0;
    RhsSupportLength_CornerNodes = 0;
    
    return;
  }

  int iElement =0;
  MatrixRowNonZeroElementTable[iElement].i=i;
  MatrixRowNonZeroElementTable[iElement].j=j;
  MatrixRowNonZeroElementTable[iElement].k=k;
  MatrixRowNonZeroElementTable[iElement].MatrixElementValue=0.0;
  MatrixRowNonZeroElementTable[iElement].iVar=0;
  MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]=0.0;
  MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength=1;
  MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength =0;

  for (int idim=0;idim<3;idim++)
    MatrixRowNonZeroElementTable[iElement].MatrixElementValue -=2*CoeffSqr[idim];
  
  iElement++;
  int indexAddition[2] = {-1,1};
  //  char * CenterNodeDataOffset = node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

  for (int ii=0;ii<3;ii++){
    for (int jj=0;jj<2;jj++){
      int iNode = i+(ii==0?indexAddition[jj]:0);
      int jNode = j+(ii==1?indexAddition[jj]:0);
      int kNode = k+(ii==2?indexAddition[jj]:0);
        
      MatrixRowNonZeroElementTable[iElement].i=iNode;
      MatrixRowNonZeroElementTable[iElement].j=jNode;
      MatrixRowNonZeroElementTable[iElement].k=kNode;
   
      MatrixRowNonZeroElementTable[iElement].MatrixElementValue=CoeffSqr[ii];
      MatrixRowNonZeroElementTable[iElement].iVar=0;
      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]=0.0;
      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength=1;
      MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength =0;

      iElement++;
    }
  }

  NonZeroElementsFound=iElement;

  RhsSupportTable_CenterNodes[0].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();
  RhsSupportTable_CenterNodes[0].Coefficient=1.0;
 
  RhsSupportLength_CenterNodes=1;   
   
  
  rhs=0.0;
  /*
  double * CenterOffset = ((double*)(RhsSupportTable_CornerNodes[0].AssociatedDataPointer+
                                     PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset));

  rhs = -4*Pi*CenterOffset[netChargeIndex]+CenterOffset[divEIndex];
  */
  return;
}


void PIC::FieldSolver::Electromagnetic::ECSIM::GetStencil(int i,int j,int k,int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,
			     cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1>::cRhsSupportTable* RhsSupportTable_CornerNodes,int& RhsSupportLength_CornerNodes,
			     cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1>::cRhsSupportTable* RhsSupportTable_CenterNodes,int& RhsSupportLength_CenterNodes, 
			     cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    
  
  // No.0-No.26  stencil Ex
  // No.27-No.53 stencil Ey
  // No.54-No.80 stencil Ez
  
  // i+indexadd[ii](ii:0,1,2), j+indexAdd[jj](jj:0,1,2), k+indexAdd[kk](kk:0,1,2)
  // index number: ii+3*jj+9*kk 
  // No.0: i,j,k No.1 i-1,j,k No.2 i+1,j,k
  // No.3: i,j-1,k No.4:i-1,j-1,k No.5 i+1,j-1,k 

  
  // double cLighgt, dt;
  double x[3];
  //power of 3 array created
  //for some pgi compilers that cannot convert result of pow() from double to int correctly
  int pow3_arr[3]={1,3,9};
  int index[3] = {i,j,k};
  int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
  double dx[3],coeff[3],coeffSqr[3],coeff4[3]; 

  for (int iDim=0; iDim<3; iDim++){
    dx[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim];
    //convert length
    dx[iDim] *= length_conv;

    x[iDim]=node->xmin[iDim]+index[iDim]*(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim];

    coeff[iDim] = cDt/dx[iDim]*theta; // for  test purpose
    coeffSqr[iDim] = coeff[iDim]*coeff[iDim];

    coeff4[iDim] = coeff[iDim]*0.25; //coefficients for curl calculation

  }

#if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ 
 
  if (isBoundaryCorner(x,node)) {  
    MatrixRowNonZeroElementTable[0].i=i;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=0.0;
    MatrixRowNonZeroElementTable[0].iVar=iVar;
    MatrixRowNonZeroElementTable[0].MatrixElementParameterTable[0]=1.0;
    MatrixRowNonZeroElementTable[0].MatrixElementParameterTableLength=1;
    MatrixRowNonZeroElementTable[0].MatrixElementSupportTableLength = 0;
    MatrixRowNonZeroElementTable[0].Node=node;
    MatrixRowNonZeroElementTable[0].MatrixElementSupportTable[0]=NULL;
    NonZeroElementsFound =1;
    rhs =0.0;
    RhsSupportLength_CenterNodes = 0;
    RhsSupportLength_CornerNodes = 0;
    if  (isRightBoundaryCorner(x,node)) NonZeroElementsFound =0;   
    
    return;
  }

#endif


  if (!initMassMatrixOffsetTable) computeMassMatrixOffsetTable(); 

  int indexAddition[3] = {0,-1,1};
  char * NodeDataOffset = node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){
    for (int ii=0;ii<3;ii++){
      for (int jj=0;jj<3;jj++){
        for (int kk=0;kk<3;kk++){
          int iNode = i+indexAddition[ii];
          int jNode = j+indexAddition[jj];
          int kNode = k+indexAddition[kk];
          int iElement = iVarIndex*27+ii+jj*3+kk*9;

          MatrixRowNonZeroElementTable[iElement].i=iNode;
          MatrixRowNonZeroElementTable[iElement].j=jNode;
          MatrixRowNonZeroElementTable[iElement].k=kNode;
            // already initialized in LinearSystemSolver.h
          MatrixRowNonZeroElementTable[iElement].MatrixElementValue=0.0;
          MatrixRowNonZeroElementTable[iElement].iVar=iVarIndex;
          MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]=0.0;
          MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength=1;
          MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength = 1;
          //question question
          MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTable[0]=(double*)NodeDataOffset+MassMatrixOffsetIndex+MassMatrixOffsetTable[iVar][iElement];
        }
      }
    }
  }
    
  //laplacian
  
  for (int ii=0;ii<3;ii++){
    for (int jj=0;jj<3;jj++){
      for (int kk=0;kk<3;kk++){
        int index=ii+jj*3+kk*9;
        int iElement = index + iVar*27;
        //minus laplacian
        MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]-=
          LaplacianCoeff[0][index]*coeffSqr[0]+LaplacianCoeff[1][index]*coeffSqr[1]+
          LaplacianCoeff[2][index]*coeffSqr[2];
      }
    }
  }
  

  //plus self
  MatrixRowNonZeroElementTable[27*iVar].MatrixElementParameterTable[0]+=1;


  //plus graddiv E
  
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    for (int ii=0;ii<3;ii++){
      for (int jj=0;jj<3;jj++){
        for (int kk=0;kk<3;kk++){
          int nodeIndex=ii+jj*3+kk*9;
          int iElement = nodeIndex + iVarIndex*27;
          MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]+=
            graddiv[iVar][iVarIndex][nodeIndex]*coeff[iVar]*coeff[iVarIndex];
        }
      }
    }
  }
  
    
  //find corners outside the boundary
  vector<int> pointLeft;
  int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;
  for (int ii=0;ii<81;ii++) {
    MatrixRowNonZeroElementTable[ii].Node=node;

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * nodeTemp = node;
    int i0,j0,k0; //for test
    i0 = MatrixRowNonZeroElementTable[ii].i;
    j0 = MatrixRowNonZeroElementTable[ii].j;
    k0 = MatrixRowNonZeroElementTable[ii].k;

    if ((i0>=iMax) && (nodeTemp!=NULL)) {
      i0-=_BLOCK_CELLS_X_;
      nodeTemp=nodeTemp->GetNeibFace(1,0,0);
    }
    else if (i0<0 && nodeTemp!=NULL) {
      i0+=_BLOCK_CELLS_X_;
      nodeTemp=nodeTemp->GetNeibFace(0,0,0);
    }

    if ((j0>=jMax) && (nodeTemp!=NULL)) {
      j0-=_BLOCK_CELLS_Y_;
      nodeTemp=nodeTemp->GetNeibFace(3,0,0);
    }
    else if (j0<0 && nodeTemp!=NULL) {
      j0+=_BLOCK_CELLS_Y_;
      nodeTemp=nodeTemp->GetNeibFace(2,0,0);
    }

    if ((k0>=kMax) && (nodeTemp!=NULL)) {
      k0-=_BLOCK_CELLS_Z_;
      nodeTemp=nodeTemp->GetNeibFace(5,0,0);
    }
    else if (k0<0 && nodeTemp!=NULL) {
      k0+=_BLOCK_CELLS_Z_;
      nodeTemp=nodeTemp->GetNeibFace(4,0,0);
    }

    double xlocal[3];
    int indlocal[3]={MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k};
    int indexG_local[3];
    bool isFixed=false;
    for (int idim=0; idim<3; idim++) {
      xlocal[idim]=MatrixRowNonZeroElementTable[ii].Node->xmin[idim]+indlocal[idim]*dx[idim];
    }
     
    pointLeft.push_back(ii);
    if (MatrixRowNonZeroElementTable[ii].Node==NULL){
      pointLeft.pop_back();
      continue;
    }else if (MatrixRowNonZeroElementTable[ii].Node->Thread==-1){
      pointLeft.pop_back();
      continue;
    }
    
    if (nodeTemp==NULL){
      pointLeft.pop_back();
      continue;
    }else if (nodeTemp->Thread==-1){
      pointLeft.pop_back();
      continue;
    }

    if (isBoundaryCorner(xlocal,node)) pointLeft.pop_back();
  }

  for (int ii=0; ii<pointLeft.size();ii++){

    int copyFrom = pointLeft[ii];

    if (ii!=copyFrom){
      MatrixRowNonZeroElementTable[ii].i=
        MatrixRowNonZeroElementTable[copyFrom].i;
      MatrixRowNonZeroElementTable[ii].j=
        MatrixRowNonZeroElementTable[copyFrom].j;
      MatrixRowNonZeroElementTable[ii].k=
        MatrixRowNonZeroElementTable[copyFrom].k;
      MatrixRowNonZeroElementTable[ii].MatrixElementValue=
        MatrixRowNonZeroElementTable[copyFrom].MatrixElementValue;
      MatrixRowNonZeroElementTable[ii].iVar=
        MatrixRowNonZeroElementTable[copyFrom].iVar;
      MatrixRowNonZeroElementTable[ii].MatrixElementParameterTable[0]=
        MatrixRowNonZeroElementTable[copyFrom].MatrixElementParameterTable[0];
      MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength=
        MatrixRowNonZeroElementTable[copyFrom].MatrixElementParameterTableLength;
      MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength = 
        MatrixRowNonZeroElementTable[copyFrom].MatrixElementSupportTableLength;
         
      MatrixRowNonZeroElementTable[ii].MatrixElementSupportTable[0]=
        MatrixRowNonZeroElementTable[copyFrom].MatrixElementSupportTable[0];
    }

  }
  
  NonZeroElementsFound=pointLeft.size();
  
  //NonZeroElementsFound=81;

  //  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){

    // fill first 27 elements
  for (int ii=0;ii<3;ii++){
    for (int jj=0;jj<3;jj++){
      for (int kk=0;kk<3;kk++){
        int iNode = i+indexAddition[ii];
        int jNode = j+indexAddition[jj];
        int kNode = k+indexAddition[kk];
        int iElement = ii+jj*3+kk*9;

        RhsSupportTable_CornerNodes[iElement].Coefficient= 0.0;
        RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer=node->block->GetCornerNode(_getCornerNodeLocalNumber(iNode,jNode,kNode))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
      }
    }
  }


  for (int iVarIndex=1; iVarIndex<3; iVarIndex++){
    // fill next 54 elements
    for (int ii=0;ii<3;ii++){
      for (int jj=0;jj<3;jj++){
        for (int kk=0;kk<3;kk++){
          int iNode = i+indexAddition[ii];
          int jNode = j+indexAddition[jj];
          int kNode = k+indexAddition[kk];
          int iElement = iVarIndex*27+ii+jj*3+kk*9;
          int jOldElement = ii+jj*3+kk*9;

          RhsSupportTable_CornerNodes[iElement].Coefficient= 0.0;
          RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer=RhsSupportTable_CornerNodes[jOldElement].AssociatedDataPointer;
        }
      }
    }
  }

  //laplacian
  for (int ii=0;ii<3;ii++){
    for (int jj=0;jj<3;jj++){
      for (int kk=0;kk<3;kk++){
        int index=ii+jj*3+kk*9;
        int iElement = index + iVar*27;
        //plus laplacian
        RhsSupportTable_CornerNodes[iElement].Coefficient +=
          LaplacianCoeff[0][index]*coeffSqr[0]+LaplacianCoeff[1][index]*coeffSqr[1]+
          LaplacianCoeff[2][index]*coeffSqr[2];
      }
    }
  }

  
  //minus graddiv E
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    for (int ii=0;ii<3;ii++){
      for (int jj=0;jj<3;jj++){
        for (int kk=0;kk<3;kk++){
          int nodeIndex=ii+jj*3+kk*9;
          int iElement = nodeIndex + iVarIndex*27;
          RhsSupportTable_CornerNodes[iElement].Coefficient -=
            graddiv[iVar][iVarIndex][nodeIndex]*coeff[iVar]*coeff[iVarIndex];
        }
      }
    }
  }



  RhsSupportTable_CornerNodes[81].AssociatedDataPointer=RhsSupportTable_CornerNodes[iVar*27].AssociatedDataPointer;
  RhsSupportTable_CornerNodes[81].Coefficient=-4*Pi*dtTotal*theta;
 
  RhsSupportLength_CornerNodes=82;

  //Ex^n,Ey^n,Ez^n
  rhs=0.0;
  
  
  int indexAdditionB[2] = {-1,0};
  
  int iElement = 0;
  
  double curlB = 0.0;
    //Ex  rhs+= d Bz/dy - d By/dz
  if (iVar==0){
          
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[1]; //c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j,k+indexAdditionB[jj]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //  rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;

        }
      }

      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[1]; //-c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j-1,k+indexAdditionB[jj]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
      
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[2]; //-c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //  rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[2]; //c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k-1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
    }

     //Ey  rhs+= d Bx/dz - d Bz/dx
    if (iVar==1){
     
      
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[2]; //c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          // curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;

        }
      }

      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[2]; //-c(dt)/dz
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[ii],j+indexAdditionB[jj],k-1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
      
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[0]; //-c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[0]; //c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i-1,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          // rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BzOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
      
      // double analytic = -1000*3.14159265/2*cos((x[0]+1)*3.14159265/2)*0.2;
      //printf("Ey,curlB:%f,analytic:%f\n", curlB, analytic);
      //rhs+=curlB;
    }
    
    //Ez  rhs+= d By/dx - d Bx/dy
    if (iVar==2){
     
     
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[0]; //c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[ii].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[ii].Coefficient;
          iElement++;
        }
      }

      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[0]; //-c(dt)/dx
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i-1,j+indexAdditionB[jj],k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[ByOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
      
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=-coeff4[1]; //-c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[jj],j,k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }
    
      for (int ii=0;ii<2;ii++){
        for (int jj=0;jj<2;jj++){
          RhsSupportTable_CenterNodes[iElement].Coefficient=coeff4[1]; //c(dt)/dy
          RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+indexAdditionB[jj],j-1,k+indexAdditionB[ii]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          //rhs+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          //curlB+=((double*)(RhsSupportTable_CenterNodes[iElement].AssociatedDataPointer+CurrentCenterNodeOffset))[BxOffsetIndex]*RhsSupportTable_CenterNodes[iElement].Coefficient;
          iElement++;
        }
      }

      //double analytic = -1000*3.14159265/2*cos((x[0]+1)*3.14159265/2)*0.2;
      //printf("Ez,curlB:%f,analytic:%f\n", curlB, analytic);
      //rhs+=curlB;
    }
   
    RhsSupportLength_CenterNodes = iElement;     
}



int IndexMatrix[8][8]={{0,2,8,6,18,20,26,24},{1,0,6,7,19,18,24,25},{4,3,0,1,22,21,18,19},
			  {3,5,2,0,21,23,20,18},{9,11,17,15,0,2,8,6},{10,9,15,16,1,0,6,7},
			  {13,12,9,10,4,3,0,1},{12,14,11,9,3,5,2,0}};



void PIC::FieldSolver::Electromagnetic::ECSIM::ProcessJMassMatrix(char * realData, char * ghostData){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  
  realData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  ghostData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

  for (int ii=0;ii<3;ii++)   {
    ((double*)realData)[JxOffsetIndex+ii]+=((double*)ghostData)[JxOffsetIndex+ii];
  }

  for (int ii=0;ii<243;ii++) {
    ((double*)realData)[MassMatrixOffsetIndex+ii]+=((double*)ghostData)[MassMatrixOffsetIndex+ii];
  }

}


void PIC::FieldSolver::Electromagnetic::ECSIM::ProcessJMassMatrixSpeciesData(char * realData, char * ghostData){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  
  realData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  ghostData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

  for (int ii=0;ii<3;ii++)   {
    ((double*)realData)[JxOffsetIndex+ii]+=((double*)ghostData)[JxOffsetIndex+ii];
  }

  for (int ii=0;ii<243;ii++) {
    ((double*)realData)[MassMatrixOffsetIndex+ii]+=((double*)ghostData)[MassMatrixOffsetIndex+ii];
  }

  for (int ii=0;ii<10*PIC::nTotalSpecies;ii++) {
    ((double*)realData)[SpeciesDataIndex[0]+ii]+=((double*)ghostData)[SpeciesDataIndex[0]+ii];
  }


}

void PIC::FieldSolver::Electromagnetic::ECSIM::ProcessNetCharge(char * realData, char * ghostData){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  
  realData+=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
  ghostData+=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

  ((double*)realData)[netChargeNewIndex]+=((double*)ghostData)[netChargeNewIndex];


}


void PIC::FieldSolver::Electromagnetic::ECSIM::CopyNetCharge(char * ghostData, char * realData){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  
  realData+=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
  ghostData+=PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

  ((double*)ghostData)[netChargeNewIndex]=((double*)realData)[netChargeNewIndex];
 
}





void PIC::FieldSolver::Electromagnetic::ECSIM::CopyJMassMatrix(char * ghostData, char * realData){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  
  realData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  ghostData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

  for (int ii=0;ii<3;ii++)   {
    ((double*)ghostData)[JxOffsetIndex+ii]=((double*)realData)[JxOffsetIndex+ii];
  }

  for (int ii=0;ii<243;ii++) {
    ((double*)ghostData)[MassMatrixOffsetIndex+ii]=((double*)realData)[MassMatrixOffsetIndex+ii];

  }

}


void PIC::FieldSolver::Electromagnetic::ECSIM::CopyJMassMatrixSpeciesData(char * ghostData, char * realData){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  
  realData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
  ghostData+=PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

  for (int ii=0;ii<3;ii++)   {
    ((double*)ghostData)[JxOffsetIndex+ii]=((double*)realData)[JxOffsetIndex+ii];
  }

  for (int ii=0;ii<243;ii++) {
    ((double*)ghostData)[MassMatrixOffsetIndex+ii]=((double*)realData)[MassMatrixOffsetIndex+ii];
  }

  for (int ii=0;ii<10*PIC::nTotalSpecies;ii++) {
    ((double*)ghostData)[SpeciesDataIndex[0]+ii]=((double*)realData)[SpeciesDataIndex[0]+ii];
  }


}

void PIC::FieldSolver::Electromagnetic::ECSIM::ComputeDivE(){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM; 


  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
    
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;

    
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    PIC::Mesh::cDataBlockAMR * block=node->block;
    
    double CellVolume=1;
    double dx[3];
    for (int iDim=0; iDim<3;iDim++) dx[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim]*length_conv;    
    for (int iDim=0; iDim<3;iDim++) CellVolume*=dx[iDim];
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          
          char * offset[8];
          double *CornerECurr[8];
                               
          offset[0]=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          offset[1]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          offset[2]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,j+1,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          offset[3]=block->GetCornerNode(_getCornerNodeLocalNumber(i,  j+1,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          offset[4]=block->GetCornerNode(_getCornerNodeLocalNumber(i,    j,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          offset[5]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,  j,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          offset[6]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,j+1,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          offset[7]=block->GetCornerNode(_getCornerNodeLocalNumber(i,  j+1,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
          
          for (int ii=0; ii<8; ii++) {
            CornerECurr[ii]=((double*)(offset[ii]+CurrentEOffset));
          }                

              
          double divE = 0.0;
          int positiveCorner[3][4]={{1,2,5,6},{2,3,6,7},{4,5,6,7}};
          int negativeCorner[3][4]={{0,3,4,7},{0,1,4,5},{0,1,2,3}};
          
          for (int iDim=0; iDim<3; iDim++){
            for (int iCorner=0;iCorner<4;iCorner++){
              divE += (CornerECurr[positiveCorner[iDim][iCorner]][iDim]
                       - CornerECurr[negativeCorner[iDim][iCorner]][iDim])/dx[iDim];
              
            }
          }
                
          double * CenterOffset = (double *)(block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset);
                
          CenterOffset[divEIndex] = divE * CellVolume*0.25;

        }
      }
    }

  }

  
}


void PIC::FieldSolver::Electromagnetic::ECSIM::UpdateJMassMatrix(){
  // update J and MassMatrix 
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM; 
  //the table of cells' particles
  //long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  long int *FirstCellParticleTable;
  PIC::ParticleBuffer::byte *ParticleData,*ParticleDataNext;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber,ptr,ptrNext;    

  double ParticleEnergy=0.0;

  PIC::Mesh::SetCornerNodeAssociatedDataValue(0.0,3,JxOffsetIndex*sizeof(double)+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset);
  PIC::Mesh::SetCornerNodeAssociatedDataValue(0.0,243,MassMatrixOffsetIndex*sizeof(double)+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset);

#if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_  
  PIC::Mesh::SetCornerNodeAssociatedDataValue(0.0,10*PIC::nTotalSpecies,SpeciesDataIndex[0]*sizeof(double)+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset);
#endif

  double qom[PIC::nTotalSpecies];
  for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++) 
    qom[iSp] = (PIC::MolecularData::GetElectricCharge(iSp)*charge_conv)
      /(PIC::MolecularData::GetMass(iSp)*mass_conv); 


  int nparticle=0;
  // update J and MassMatrix
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
    double StartTime=MPI_Wtime();

    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;
     
#if  _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CENTER_BASED_  
    double B_Center[_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_][3];
    for (int k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
      for (int j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)  {
        for (int i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
          int LocalCenterId = _getCenterNodeLocalNumber(i,j,k);
          if (!node->block->GetCenterNode(LocalCenterId)) continue; 
          char *offset=node->block->GetCenterNode(LocalCenterId)->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
          double * ptr =  (double*)(offset+CurrentBOffset);
	  //for (int idim=0;idim<3;idim++) B_Center[LocalCenterId][idim]=ptr[idim];
          memcpy(B_Center[LocalCenterId],ptr,3*sizeof(double));
        }
      }
    }
#endif

#if  _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CORNER_BASED_  
    double B_corner[(_TOTAL_BLOCK_CELLS_X_+1)*(_TOTAL_BLOCK_CELLS_Y_+1)*(_TOTAL_BLOCK_CELLS_Z_+1)][3];
    for (int k=0;k<_BLOCK_CELLS_Z_+1;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_+1;j++)  {
        for (int i=0;i<_BLOCK_CELLS_X_+1;i++) {
          int LocalCornerId = _getCornerNodeLocalNumber(i,j,k);
          if (!node->block->GetCornerNode(LocalCornerId)) continue; 
          char *offset=node->block->GetCornerNode(LocalCornerId)->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+OffsetB_corner;
          double * ptr =  (double*)(offset+CurrentBOffset);
	  //for (int idim=0;idim<3;idim++) B_Center[LocalCenterId][idim]=ptr[idim];
          memcpy(B_corner[LocalCornerId],ptr,3*sizeof(double));
        }
      }
    }
#endif

    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    
    block=node->block;
    
    //memcpy(FirstCellParticleTable,block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
    FirstCellParticleTable=block->FirstCellParticleTable;
    double CellVolume=1;
    double dx[3];
    for (int iDim=0; iDim<3;iDim++) dx[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim]*length_conv;      
    for (int iDim=0; iDim<3;iDim++) CellVolume*=dx[iDim];
  
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
	  
          if (ptr!=-1) {

	    // printf("particle, i,j,k,ptr:%d,%d,%d,%ld\n",i,j,k,ptr);	   
	    //iPar=i;jPar=j;kPar=k;
	    //ParticleNode = node;

	    //  LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            //cell=block->GetCenterNode(LocalCellNumber);
	    double vInit[3]={0.0,0.0,0.0},xInit[3]={0.0,0.0,0.0};
	    int spec;
	    double Jg[8][3];

	    for (int ii=0; ii<8; ii++){
	      for (int jj=0; jj<3; jj++){
                Jg[ii][jj]=0.0;
	      }
	    }
            
      double MassMatrix_GGD[8][8][9];
      for (int iCorner=0;iCorner<8;iCorner++){
        for(int jCorner=0;jCorner<8;jCorner++){
          for (int idim=0;idim<9;idim++){
            MassMatrix_GGD[iCorner][jCorner][idim] = 0.0;
          }
        }
      }

#if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_            
      double SpeciesData_GI[8][PIC::nTotalSpecies*10];
      for (int ii=0; ii<8; ii++){
        for (int kk=0; kk<10*PIC::nTotalSpecies; kk++){
          SpeciesData_GI[ii][kk]=0.0;
        }
      }
#endif            
	    ptrNext=ptr;
	    ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
	  
	    while (ptrNext!=-1) {
	      nparticle++;
	      double LocalParticleWeight;
	      ptr=ptrNext;
	      ParticleData=ParticleDataNext;	  	    
	     
              spec=PIC::ParticleBuffer::GetI(ParticleData);
	      PIC::ParticleBuffer::GetV(vInit,ParticleData);
	      PIC::ParticleBuffer::GetX(xInit,ParticleData);
              LocalParticleWeight=block->GetLocalParticleWeight(spec);
              LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
              
              ptrNext=PIC::ParticleBuffer::GetNext(ParticleData);

        if (ptrNext!=-1) {
          ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptrNext);
          PIC::ParticleBuffer::PrefertchParticleData_Basic(ParticleDataNext);
        }


	      double temp[3], B[3]={0.0,0.0,0.0};
#if _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CENTER_BASED_
	      PIC::InterpolationRoutines::CellCentered::cStencil MagneticFieldStencil(false);
	      //interpolate the magnetic field from center nodes to particle location
	      MagneticFieldStencil=*(PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xInit,node));
#endif

#if _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CORNER_BASED_
              PIC::InterpolationRoutines::CornerBased::cStencil MagneticFieldStencil(false);
	      //interpolate the magnetic field from center nodes to particle location
	      MagneticFieldStencil=*(PIC::InterpolationRoutines::CornerBased::InitStencil(xInit,node));             
#endif

	      for (int iStencil=0;iStencil<MagneticFieldStencil.Length;iStencil++) {
		//memcpy(temp,MagneticFieldStencil.cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+CurrentBOffset,3*sizeof(double));
#if _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CENTER_BASED_
		double * B_temp = B_Center[MagneticFieldStencil.LocalCellID[iStencil]];
#endif
#if _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CORNER_BASED_
                double * B_temp = B_corner[MagneticFieldStencil.LocalCellID[iStencil]];
#endif
		for (int idim=0;idim<3;idim++) B[idim]+=MagneticFieldStencil.Weight[iStencil]*B_temp[idim];
	      }

	      //convert from SI to cgs
	      for (int idim=0; idim<3; idim++){
                B[idim] *= B_conv;
                vInit[idim] *= length_conv;
	      }

	      double QdT_over_m,QdT_over_2m,alpha[9],chargeQ;
	      double * WeightPG;
	      double c0,QdT_over_2m_squared;
	      double mass;

	      chargeQ = PIC::MolecularData::GetElectricCharge(spec)*charge_conv; 
	      mass= PIC::MolecularData::GetMass(spec)*mass_conv;
	      //effect of particle weight
	     
	      chargeQ *= LocalParticleWeight;
	      mass *= LocalParticleWeight;

	      QdT_over_m=chargeQ*dtTotal/mass;
	      QdT_over_2m=0.5*QdT_over_m;
	      QdT_over_2m_squared=QdT_over_2m*QdT_over_2m;


	      //to calculate alpha, mdv/dt = q(E+v cross B/c)
	      for (int idim=0; idim<3; idim++){
                B[idim] /= LightSpeed; //divided by the speed of light
              }
	      

	      double BB[3][3],P[3];
	      
	      for (int ii=0;ii<3;ii++) {
                P[ii]=-QdT_over_2m*B[ii];

                for (int jj=0;jj<=ii;jj++) {
                  BB[ii][jj]=QdT_over_2m_squared*B[ii]*B[jj];
                  BB[jj][ii]=BB[ii][jj];
                }
	      }

	      c0=1.0/(1.0+QdT_over_2m_squared*(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));
	      
	      //compute alpha
              /*
	      alpha[0][0]=c0*(1.0+BB[0][0]);
	      alpha[0][1]=c0*(-P[2]+BB[0][1]);
	      alpha[0][2]=c0*(P[1]+BB[0][2]);
	      
	      alpha[1][0]=c0*(P[2]+BB[1][0]);
	      alpha[1][1]=c0*(1.0+BB[1][1]);
	      alpha[1][2]=c0*(-P[0]+BB[1][2]);

	      alpha[2][0]=c0*(-P[1]+BB[2][0]);
	      alpha[2][1]=c0*(P[0]+BB[2][1]);
	      alpha[2][2]=c0*(1.0+BB[2][2]);
              */
              alpha[0]=c0*(1.0+BB[0][0]);
	      alpha[1]=c0*(-P[2]+BB[0][1]);
	      alpha[2]=c0*(P[1]+BB[0][2]);
	      
	      alpha[3]=c0*(P[2]+BB[1][0]);
	      alpha[4]=c0*(1.0+BB[1][1]);
	      alpha[5]=c0*(-P[0]+BB[1][2]);

	      alpha[6]=c0*(-P[1]+BB[2][0]);
	      alpha[7]=c0*(P[0]+BB[2][1]);
	      alpha[8]=c0*(1.0+BB[2][2]);


	      //get weight for each corner
	      // PIC::InterpolationRoutines::CornerBased::cStencil WeightStencil(false);
	      PIC::InterpolationRoutines::CornerBased::InitStencil(xInit,node);
	    
	      /*
	      double xMinCell[3],xMaxCell[3];
	     
	      xMinCell[0]= node->xmin[0]+(node->xmax[0]-node->xmin[0])/nCell[0]*i;
	      xMinCell[1]= node->xmin[1]+(node->xmax[1]-node->xmin[1])/nCell[1]*j;
	      xMinCell[2]= node->xmin[2]+(node->xmax[2]-node->xmin[2])/nCell[2]*k;
	      
	      xMaxCell[0]=node->xmin[0]+(node->xmax[0]-node->xmin[0])/nCell[0]*(i+1);
	      xMaxCell[1]=node->xmin[1]+(node->xmax[1]-node->xmin[1])/nCell[1]*(j+1);
	      xMaxCell[2]=node->xmin[2]+(node->xmax[2]-node->xmin[2])/nCell[2]*(k+1);
	      */

	      WeightPG=PIC::InterpolationRoutines::CornerBased::InterpolationCoefficientTable_LocalNodeOrder;
	      
            /*
        double tempWeight[3][2];
        for (int idim=0;idim<3;idim++){
          printf("xinit:%e,xmin:%e,xmax:%e\n",xInit[idim],xMinCell[idim],xMaxCell[idim]);
          tempWeight[idim][1]= (xInit[idim]-xMinCell[idim])/(xMaxCell[idim]-xMinCell[idim]);
          tempWeight[idim][0]= 1-tempWeight[idim][1];
          printf("weight1:%e,weight2:%e\n", tempWeight[idim][0],tempWeight[idim][1]);
        }

        for (int ii=0;ii<2;ii++){
          for (int jj=0;jj<2;jj++){
            for (int kk=0;kk<2;kk++){
              printf("weight product:%e\n", tempWeight[0][ii]*tempWeight[1][jj]*tempWeight[2][kk]);
            }
          }
        }
            */

	      
	      ParticleEnergy += 0.5*mass*(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2]);  
		//compute alpha*vInit
	      double vRot[3]={0.0,0.0,0.0};

	      for (int iDim =0; iDim<3; iDim++){
                for (int jj=0; jj<3; jj++){
                  vRot[iDim]+=alpha[3*iDim+jj]*vInit[jj];
                }
	      }

	      // printf("vRot[iDim]:%e,%e,%e\n",vRot[0],vRot[1],vRot[2]);
	      
	      for (int iCorner=0; iCorner<8; iCorner++){
                for (int iDim=0; iDim<3; iDim++){
                  Jg[iCorner][iDim]+=chargeQ*vRot[iDim]*WeightPG[iCorner];

            //  printf("Jg[iCorner][iDim]:%e\n",Jg[iCorner][iDim]);

                }
	      }
      
#if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_        
              for (int ii=0; ii<8; ii++){
                int tempOffset = 10*spec;
                SpeciesData_GI[ii][tempOffset+Rho_]+=mass*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUx_]+=mass*vInit[0]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUy_]+=mass*vInit[1]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUz_]+=mass*vInit[2]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUxUx_]+=mass*vInit[0]*vInit[0]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUyUy_]+=mass*vInit[1]*vInit[1]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUzUz_]+=mass*vInit[2]*vInit[2]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUxUy_]+=mass*vInit[0]*vInit[1]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUyUz_]+=mass*vInit[1]*vInit[2]*WeightPG[ii];
                SpeciesData_GI[ii][tempOffset+RhoUxUz_]+=mass*vInit[0]*vInit[2]*WeightPG[ii];
              }
#endif
              
              double matrixConst = chargeQ*QdT_over_2m/CellVolume;
	      for (int iCorner=0; iCorner<8; iCorner++){
                double tempWeightConst = matrixConst*WeightPG[iCorner];
                for (int jCorner=0; jCorner<=iCorner; jCorner++){
                  double tempWeightProduct = WeightPG[jCorner]*tempWeightConst;
                  double * tmpPtr =MassMatrix_GGD[iCorner][jCorner];
                    //for (int ii=0; ii<3; ii++){
                    // for (int jj=0; jj<3; jj++){
                    for (int ii=0; ii<9; ii++){   
                    //double tmp = alpha[ii][jj]*tempWeightProduct;
                    double tmp = alpha[ii]*tempWeightProduct;   
                    //  CornerMassMatrixPtr[iCorner][9*IndexMatrix[iCorner][jCorner]+3*ii+jj] += tmp;
                        //CornerMassMatrixPtr[jCorner][9*IndexMatrix[jCorner][iCorner]+3*ii+jj] += tmp;
                        tmpPtr[ii] +=tmp;
                    }
                }//jCorner
	      }//iCorner
	      
            
	      if (ptrNext!=-1) {
// do nothing;ParticleDataNext is determined earlier in the loop; ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptrNext);
	      } else {

                double * CornerMassMatrixPtr[8];
                double * CornerJPtr[8]; 
                char * offset[8];
                double * specDataPtr[8];
                                             
                offset[0]=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                offset[1]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                offset[2]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,j+1,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                offset[3]=block->GetCornerNode(_getCornerNodeLocalNumber(i,  j+1,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                offset[4]=block->GetCornerNode(_getCornerNodeLocalNumber(i,    j,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                offset[5]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,  j,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                offset[6]=block->GetCornerNode(_getCornerNodeLocalNumber(i+1,j+1,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                offset[7]=block->GetCornerNode(_getCornerNodeLocalNumber(i,  j+1,k+1))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                

                for (int ii=0; ii<8; ii++) {
                  CornerMassMatrixPtr[ii] = ((double*)offset[ii])+MassMatrixOffsetIndex;
                  CornerJPtr[ii]=((double*)offset[ii])+JxOffsetIndex;
                  #if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_
                  specDataPtr[ii]=((double*)offset[ii])+SpeciesDataIndex[0];
                  #endif
                }                
                
                //collect current
                for (int iCorner=0; iCorner<8; iCorner++){
                  for (int ii=0; ii<3; ii++){
                    CornerJPtr[iCorner][ii] += (Jg[iCorner][ii])/CellVolume;
                  }                  
                }

                #if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_
                //collect species data
                for (int iCorner=0; iCorner<8; iCorner++){
                  for (int ii=0; ii<10*PIC::nTotalSpecies; ii++){
                    specDataPtr[iCorner][ii] += SpeciesData_GI[iCorner][ii]/CellVolume;
                  }                  
                }
                #endif
                
                //collect massmatrix
                for (int iCorner=0; iCorner<8; iCorner++){
                  for (int jCorner=0; jCorner<=iCorner; jCorner++){
                 
                    if (iCorner==jCorner){
                      for (int ii=0; ii<3; ii++){
                        for (int jj=0; jj<3; jj++){
                          
                          CornerMassMatrixPtr[iCorner][3*ii+jj] += MassMatrix_GGD[iCorner][iCorner][3*ii+jj];
                          // MassMatrix_GGD[iCorner][iCorner][3*ii+jj] +=tmp;
                          //printf("CornerMassMatrix:%e\n", *(CornerMassMatrixPtr[iCorner]+3*ii+jj));
                        }
                      }
                    } else {
                      for (int ii=0; ii<3; ii++){
                        for (int jj=0; jj<3; jj++){
                       
                          CornerMassMatrixPtr[iCorner][9*IndexMatrix[iCorner][jCorner]+3*ii+jj] += MassMatrix_GGD[iCorner][jCorner][3*ii+jj];
                          CornerMassMatrixPtr[jCorner][9*IndexMatrix[jCorner][iCorner]+3*ii+jj] += MassMatrix_GGD[iCorner][jCorner][3*ii+jj];
                          //MassMatrix_GGD[iCorner][jCorner][3*ii+jj] +=tmp;
                        }
                      }
                    }

                  }//jCorner
                }//iCorner
                
                

	      }

	    }// while (ptrNext!=-1)
	  }//if (ptr!=-1)
	}// for i
      }//for j   
    }//for k

    //increment the time counter
   if (_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_) {
     node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
   }
    
  }
#if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_
  PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData=PIC::FieldSolver::Electromagnetic::ECSIM::ProcessJMassMatrixSpeciesData;
  PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData=PIC::FieldSolver::Electromagnetic::ECSIM::CopyJMassMatrixSpeciesData;
#elif  _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_OFF_
  PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData=PIC::FieldSolver::Electromagnetic::ECSIM::ProcessJMassMatrix;
  PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData=PIC::FieldSolver::Electromagnetic::ECSIM::CopyJMassMatrix;
#endif

  PIC::Parallel::CornerBlockBoundaryNodes::SetActiveFlag(true);
#if _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_ON_
  PIC::BC::ExternalBoundary::UpdateData(PackBlockData_JMassMatrixSpeciesData,UnpackBlockData_JMassMatrixSpeciesData);
#elif  _PIC_FIELD_SOLVER_SAMPLE_SPECIES_ON_CORNER_== _PIC_MODE_OFF_
  PIC::BC::ExternalBoundary::UpdateData(PackBlockData_JMassMatrix,UnpackBlockData_JMassMatrix);
#endif
  PIC::Parallel::CornerBlockBoundaryNodes::SetActiveFlag(false);
  
  MPI_Reduce(&ParticleEnergy, &TotalParticleEnergy, 1, MPI_DOUBLE, MPI_SUM, 0,
	     MPI_GLOBAL_COMMUNICATOR);
  // TotalParticleEnergy *= 1e-7; //in SI
  if (PIC::ThisThread==0) {
    printf("Total Particle Energy:%e\n",TotalParticleEnergy); 
    printf("Total Energy:%.20e,%f\n",TotalParticleEnergy+TotalWaveEnergy,TotalParticleEnergy+TotalWaveEnergy);
    std::cout.precision(20);
    std::cout<<"total energy: "<<TotalParticleEnergy+TotalWaveEnergy<<std::endl;  
  }
  
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ && _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_ )
    PIC::CPLR::FLUID::fix_plasma_node_boundary();
}



void PIC::FieldSolver::Electromagnetic::ECSIM::divECorrection(){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  ComputeNetCharge(true);
  ComputeDivE();
  SetBoundaryChargeDivE();
  PoissonSolver.UpdateRhs(PoissonUpdateRhs);
  linear_solver_matvec_c = PoissonMatvec;
  PoissonSolver.Solve(PoissonSetInitialGuess,PoissonProcessFinalSolution,1e-2,
                      PIC::CPLR::FLUID::EFieldIter,PackBlockData_phi,UnpackBlockData_phi);
  SetBoundaryPHI();
  CorrectParticleLocation();
  PIC::Parallel::ExchangeParticleData();
  ComputeNetCharge(false);
}

void exchangeParticleLocal(){
  
   long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  /*
  for (k=0;k<_BLOCK_CELLS_Z_;k++) {
    for (j=0;j<_BLOCK_CELLS_Y_;j++) {
      for (i=0;i<_BLOCK_CELLS_X_;i++) {
        
        //LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        
          FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;
      }
    }
  }
  */
  for (int ii=0;ii<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;ii++)
    FirstCellParticleTable[ii]=-1;  

  for (int thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *  node=(thread==PIC::Mesh::mesh.ThisThread) ? PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread] : PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];

    if (node==NULL) continue;


    for (;node!=NULL;node=node->nextNodeThisThread) {

      PIC::Mesh::cDataBlockAMR *block=node->block;
      if (!block) continue;
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
      memcpy(block->FirstCellParticleTable,block->tempParticleMovingListTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
      memcpy(block->tempParticleMovingListTable,FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      int thread_OpenMP;

      memcpy(block->FirstCellParticleTable,FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

      //link the lists created by each OpenMP threads
      long int FirstParticle,LastParticle=-1,pNext,*LastParticlePtr;

      for (thread_OpenMP=0;thread_OpenMP<PIC::nTotalThreadsOpenMP;thread_OpenMP++) {
        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          LastParticlePtr=block->GetTempParticleMovingListTableThread(thread_OpenMP,i,j,k);

          LastParticle=(*LastParticlePtr);

          if (LastParticle!=-1) {
            FirstParticle=LastParticle;
            pNext=PIC::ParticleBuffer::GetNext(LastParticle);

            while (pNext!=-1) {
              LastParticle=pNext;
              pNext=PIC::ParticleBuffer::GetNext(LastParticle);
            }

            //link patricle list
            long int *FirstCellParticlePtr=block->FirstCellParticleTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);

            PIC::ParticleBuffer::SetNext(*FirstCellParticlePtr,LastParticle);
            if (*FirstCellParticlePtr!=-1) PIC::ParticleBuffer::SetPrev(LastParticle,*FirstCellParticlePtr);

            *FirstCellParticlePtr=FirstParticle;
          }

          *LastParticlePtr=-1;
        }
      }

#else
#error The option is unknown
#endif


//      node=node->nextNodeThisThread;
    }
  }//for (int thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++)
  

}

void PIC::FieldSolver::Electromagnetic::ECSIM::CorrectParticleLocation(){
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM; 
  //the table of cells' particles
  //long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  long int *FirstCellParticleTable;
  PIC::ParticleBuffer::byte *ParticleData,*ParticleDataNext;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber,ptr,ptrNext;    

  double qom[PIC::nTotalSpecies];
  for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++) 
    qom[iSp] = (PIC::MolecularData::GetElectricCharge(iSp)*charge_conv)
      /(PIC::MolecularData::GetMass(iSp)*mass_conv); 

  
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
  
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;
     
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    
    block=node->block;
    
    //memcpy(FirstCellParticleTable,block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
    FirstCellParticleTable=block->FirstCellParticleTable;
    double CellVolume=1;
    double dx[3],dx_no[3];
    for (int iDim=0; iDim<3;iDim++) dx[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim]*length_conv;  
    for (int iDim=0; iDim<3;iDim++) dx_no[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim];

    for (int iDim=0; iDim<3;iDim++) CellVolume*=dx[iDim];
    
    
    double Phi[_BLOCK_CELLS_X_+2][_BLOCK_CELLS_Y_+2][_BLOCK_CELLS_Z_+2];
    //calculate grad phi for every corner node
    for (int k=-1;k<_BLOCK_CELLS_Z_+1;k++) {
      for (int j=-1;j<_BLOCK_CELLS_Y_+1;j++) {
        for (int i=-1;i<_BLOCK_CELLS_X_+1;i++) {
          int LocalCenterId = _getCenterNodeLocalNumber(i,j,k);
          
          if (!node->block->GetCenterNode(LocalCenterId)) continue;
          Phi[i+1][j+1][k+1] =  ((double *) (block->GetCenterNode(LocalCenterId)->
                   GetAssociatedDataBufferPointer()+
                   PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset))[phiIndex];
                  
        }//i
      }//j
    }//k

    int nparticle=0;
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
	  
          if (ptr!=-1) {
            
	    double xInit[3]={0.0,0.0,0.0};
	    int spec;
            double xNode[3];
            double xCell[3];
            int index[3]={i,j,k};
            for (int iDim=0;iDim<3;iDim++){ 
              xNode[iDim] = node->xmin[iDim]*length_conv+dx[iDim]*index[iDim];
              xCell[iDim] =  node->xmin[iDim]+dx_no[iDim]*(index[iDim]+0.5);
            }
            
            bool atBoundary = isBoundaryCell(xCell,dx_no,node);
            
            
            //if (isBoundaryCell(xCell,dx_no,node)) continue;
            
        
	    ptrNext=ptr;
	    ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
	  
	    while (ptrNext!=-1) {
              nparticle++;
	      double LocalParticleWeight;
	      ptr=ptrNext;
	      ParticleData=ParticleDataNext;	  	    
	     
              spec=PIC::ParticleBuffer::GetI(ParticleData);
              PIC::ParticleBuffer::GetX(xInit,ParticleData);
              LocalParticleWeight=block->GetLocalParticleWeight(spec);
              LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
              ptrNext=PIC::ParticleBuffer::GetNext(ParticleData);
              double xFinal[3];
              
              /*
              bool isTest=false;
              if (fabs(xInit[0]-13)<1.0 && fabs(xInit[1]-1)<1.0 && fabs(xInit[2]-3)<1.0) isTest=true; 
              */
            
              if (spec!=0 || atBoundary){
                for (int idim=0;idim<3;idim++) xFinal[idim]=xInit[idim];
                
              }else{
                double xRel[3];
                for (int iDim=0;iDim<3;iDim++) xRel[iDim] = (xInit[iDim]*length_conv-xNode[iDim])/dx[iDim];
              
                int iClosestNode[3];
                for (int iDim=0;iDim<3;iDim++) iClosestNode[iDim] = int(index[iDim]+round(xRel[iDim]));
              
                int ix=iClosestNode[0]+1,iy=iClosestNode[1]+1,iz=iClosestNode[2]+1;//index of phi array
                
                for (int iDim=0;iDim<3;iDim++) xRel[iDim] = xRel[iDim]>=0.5?xRel[iDim]-0.5:xRel[iDim]+0.5;

              
                double GradPhi[3];
                GradPhi[0] = 
                  interp2D(Phi[ix][iy - 1][iz - 1] - Phi[ix - 1][iy - 1][iz - 1],
                           Phi[ix][iy][iz - 1] - Phi[ix - 1][iy][iz - 1],
                           Phi[ix][iy][iz] - Phi[ix - 1][iy][iz],
                           Phi[ix][iy - 1][iz] - Phi[ix - 1][iy - 1][iz], xRel[1], xRel[2]);
                
                GradPhi[1] =
                  interp2D(Phi[ix - 1][iy][iz - 1] - Phi[ix - 1][iy - 1][iz - 1],
                           Phi[ix][iy][iz - 1] - Phi[ix][iy - 1][iz - 1],
                           Phi[ix][iy][iz] - Phi[ix][iy - 1][iz],
                           Phi[ix - 1][iy][iz] - Phi[ix - 1][iy - 1][iz], xRel[0], xRel[2]);
              
                GradPhi[2] =
                  interp2D(Phi[ix - 1][iy - 1][iz] - Phi[ix - 1][iy - 1][iz - 1],
                           Phi[ix][iy - 1][iz] - Phi[ix][iy - 1][iz - 1],
                           Phi[ix][iy][iz] - Phi[ix][iy][iz - 1],
                           Phi[ix - 1][iy][iz] - Phi[ix - 1][iy][iz - 1], xRel[0], xRel[1]);

                for (int iDim=0;iDim<3;iDim++) GradPhi[iDim] /= dx[iDim];
                
                double eChargeDens,gammaTmp=0.51,eps=0.9;
                
                int localCornerId = _getCornerNodeLocalNumber(iClosestNode[0],iClosestNode[1],iClosestNode[2]);
                eChargeDens = 
                  ((double *) (block->GetCornerNode(localCornerId)->
                               GetAssociatedDataBufferPointer()+
                               PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))[SpeciesDataIndex[0]]*(-25);
                
                
                double displacement[3],temp;
                if (eChargeDens!=0) temp = 1./(4.*Pi*eChargeDens);
                else temp = 0;
                for (int iDim=0; iDim<3; iDim++) {
                  displacement[iDim] = -eps*GradPhi[iDim]*temp;
                  //xFinal[iDim]=xInit[iDim]+displacement[iDim];
                }
                                
                double epsLimit=0.1;
                
                if (fabs( displacement[0] / dx_no[0]) > epsLimit ||
                    fabs( displacement[1] / dx_no[1]) > epsLimit ||
                    fabs( displacement[2] / dx_no[2]) > epsLimit) {
                  double dl =
                    sqrt(pow(displacement[0], 2) + pow(displacement[1], 2) + pow(displacement[2], 2));
                  for (int iDim = 0; iDim < 3; iDim++)
                    displacement[iDim] *= epsLimit * dx_no[0] / dl;
                }
                
                for (int iDim=0; iDim<3; iDim++) xFinal[iDim]=xInit[iDim]+displacement[iDim];  
                //for (int iDim=0; iDim<3; iDim++) xFinal[iDim]=xInit[iDim];  
               
              }
              //bool isTest=false;
              //if (fabs(xInit[0]-15.5)<0.5 && fabs(xInit[1]-7.5)<0.5 && fabs(xInit[2]-3.5)<0.5) isTest=true; 
              long int tempFirstCellParticle,*tempFirstCellParticlePtr;
              int ip, jp, kp;
              
              /*
              if (spec==0)
                printf("particle correction xFinal:    %e  %e  %e\n",xFinal[0],xFinal[1],xFinal[2]);
              */

              cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * newNode=PIC::Mesh::mesh.findTreeNode(xFinal,node);
              if (newNode==NULL || newNode->block==NULL){
                PIC::ParticleBuffer::DeleteParticle(ptr);
              }else{
                  if (PIC::Mesh::mesh.fingCellIndex(xFinal,ip,jp,kp,newNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");
                  
                  PIC::Mesh::cDataBlockAMR * block=newNode->block;
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
                  tempFirstCellParticlePtr=block->tempParticleMovingListTable+ip+_BLOCK_CELLS_X_*(jp+_BLOCK_CELLS_Y_*kp);
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
                  tempFirstCellParticlePtr=block->GetTempParticleMovingListTableThread(omp_get_thread_num(),ip,jp,kp);
#else
#error The option is unknown
#endif
                  tempFirstCellParticle=(*tempFirstCellParticlePtr);
                  
                  PIC::ParticleBuffer::SetX(xFinal,ParticleData);
                  
                  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
                  PIC::ParticleBuffer::SetPrev(-1,ParticleData);
                  
                  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
                  *tempFirstCellParticlePtr=ptr;
              }

              if (ptrNext!=-1) ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptrNext);
              
	    }// while (ptrNext!=-1)
	  }//if (ptr!=-1)
	}// for i
      }//for j   
    }//for k


  }//for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++)

 
  exchangeParticleLocal();

  
}



void PIC::FieldSolver::Electromagnetic::ECSIM::ComputeNetCharge(bool doUpdateOld){
  // update J and MassMatrix 
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM; 
  //the table of cells' particles
  //long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  long int *FirstCellParticleTable;
  PIC::ParticleBuffer::byte *ParticleData,*ParticleDataNext;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber,ptr,ptrNext;    

  double q_I[PIC::nTotalSpecies];
  if (doUpdateOld)
    UpdateOldNetCharge();
  PIC::Mesh::SetCenterNodeAssociatedDataValue(0.0,1,netChargeNewIndex*sizeof(double)+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset);

  for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++) 
    q_I[iSp] = PIC::MolecularData::GetElectricCharge(iSp)*charge_conv;
   
  
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
  
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;
     
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    
    block=node->block;
    
    //memcpy(FirstCellParticleTable,block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
    FirstCellParticleTable=block->FirstCellParticleTable;
    double CellVolume=1;
    double dx[3];
    for (int iDim=0; iDim<3;iDim++) dx[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim]*length_conv;  
    for (int iDim=0; iDim<3;iDim++) CellVolume*=dx[iDim];
    
    
    double q_Center[_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_];
    /*
    for (int k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
      for (int j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++) {
        for (int i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
    */
    for (int k=-1;k<_BLOCK_CELLS_Z_+1;k++) {
      for (int j=-1;j<_BLOCK_CELLS_Y_+1;j++) {
        for (int i=-1;i<_BLOCK_CELLS_X_+1;i++) {
          int LocalCenterId = _getCenterNodeLocalNumber(i,j,k);
          if (!node->block->GetCenterNode(LocalCenterId)) continue;
          q_Center[LocalCenterId]=0.0;
        }
      }
    }

    int nparticle=0;
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
	  
          if (ptr!=-1) {
            
	    double xInit[3]={0.0,0.0,0.0};
	    int spec;

            
	    ptrNext=ptr;
	    ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
	  
	    while (ptrNext!=-1) {
              nparticle++;
	      double LocalParticleWeight;
	      ptr=ptrNext;
	      ParticleData=ParticleDataNext;	  	    
	     
              spec=PIC::ParticleBuffer::GetI(ParticleData);
              PIC::ParticleBuffer::GetX(xInit,ParticleData);
              LocalParticleWeight=block->GetLocalParticleWeight(spec);
              LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
              ptrNext=PIC::ParticleBuffer::GetNext(ParticleData);


              double chargeQ = q_I[spec]*LocalParticleWeight;
	    

              PIC::InterpolationRoutines::CellCentered::cStencil NetChargeStencil(false);
	      //interpolate the magnetic field from center nodes to particle location
	      NetChargeStencil=*(PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xInit,node));

	      for (int iStencil=0;iStencil<NetChargeStencil.Length;iStencil++) {
                q_Center[NetChargeStencil.LocalCellID[iStencil]]+=NetChargeStencil.Weight[iStencil]*chargeQ;
	      }

              if (ptrNext!=-1) ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptrNext);
              
	    }// while (ptrNext!=-1)
	  }//if (ptr!=-1)
	}// for i
      }//for j   
    }//for k

    
 
    for (int k=-1;k<_BLOCK_CELLS_Z_+1;k++) {
      for (int j=-1;j<_BLOCK_CELLS_Y_+1;j++) {
        for (int i=-1;i<_BLOCK_CELLS_X_+1;i++) {
          int LocalCenterId = _getCenterNodeLocalNumber(i,j,k);
          if (!node->block->GetCenterNode(LocalCenterId)) continue;
          ((double *) (node->block->GetCenterNode(LocalCenterId)->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset))[netChargeNewIndex] += q_Center[LocalCenterId];
        
        }
      }
    }
  }//for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++)

  
  PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData=PIC::FieldSolver::Electromagnetic::ECSIM::ProcessNetCharge;
  PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData=PIC::FieldSolver::Electromagnetic::ECSIM::CopyNetCharge;
  
  PIC::Parallel::CenterBlockBoundaryNodes::SetActiveFlag(true);
  PIC::BC::ExternalBoundary::UpdateData(PackBlockData_netCharge,UnpackBlockData_netCharge);
  PIC::Parallel::CenterBlockBoundaryNodes::SetActiveFlag(false);
  
}


void PIC::FieldSolver::Electromagnetic::ECSIM::UpdateOldNetCharge(){

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM; 
  //the table of cells' particles
  //long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  long int *FirstCellParticleTable;
  PIC::ParticleBuffer::byte *ParticleData,*ParticleDataNext;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber,ptr,ptrNext;    

  
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
  
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;
     
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    
    block=node->block;
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          int LocalCenterId = _getCenterNodeLocalNumber(i,j,k);
          if (!node->block->GetCenterNode(LocalCenterId)) continue;
          double * PtrTemp =  ((double *) (block->GetCenterNode(LocalCenterId)->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset));
          PtrTemp[netChargeOldIndex] = PtrTemp[netChargeNewIndex]; 
          
        }// for i
      }//for j   
    }//for k

  }//for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++)


  
}


void PIC::FieldSolver::Electromagnetic::ECSIM::SetBoundaryPHI(){

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM; 


  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber;    

  
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
  
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;
     
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    
    block=node->block;

    double dx[3];
    for (int idim=0;idim<3;idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCell[idim];
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          int index[3]={i,j,k};
          double x[3];
          for (int idim=0;idim<3;idim++) x[idim]=node->xmin[idim]+dx[idim]*(index[idim]+0.5);
          
          int ind = isBoundaryCell(x,dx,node);
          if (!ind) continue;
          int LocalCenterId = _getCenterNodeLocalNumber(i,j,k);
          if (!node->block->GetCenterNode(LocalCenterId)) continue;
          int indNeib[3]={i,j,k};
          
          for (int ii=0; ii<3; ii++){
            int flag = ind%2;
            ind = ind/2;
            if (flag==1) {
              if (indNeib[ii]==0) indNeib[ii]+=1;
              else if (indNeib[ii]==nCell[ii]-1) indNeib[ii]-=1;
              else exit(__LINE__,__FILE__,"Error: Something is wrong!");
            }
          }
          
          
          double * PtrTemp =  ((double *) (block->GetCenterNode(LocalCenterId)->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset));
          double * PtrNeib =  ((double *) (block->GetCenterNode(_getCenterNodeLocalNumber(indNeib[0],indNeib[1],indNeib[2]))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset));
          PtrTemp[phiIndex] = PtrNeib[phiIndex]; 

        }// for i
      }//for j   
    }//for k
    
  }//for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++)

  
  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
    PIC::Mesh::mesh.ParallelBlockDataExchange(PackBlockData_phi,UnpackBlockData_phi);
    break;
    
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData(PackBlockData_phi,UnpackBlockData_phi);
    break;
  }


}
  


void PIC::FieldSolver::Electromagnetic::ECSIM::SetBoundaryChargeDivE(){

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM; 


  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber;    

  
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
  
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;
     
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    
    block=node->block;

    double dx[3];
    for (int idim=0;idim<3;idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCell[idim];
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          int index[3]={i,j,k};
          double x[3];
          for (int idim=0;idim<3;idim++) x[idim]=node->xmin[idim]+dx[idim]*(index[idim]+0.5);
          if (!isBoundaryCell(x,dx,node)) continue;
          
          int LocalCenterId = _getCenterNodeLocalNumber(i,j,k);
          if (!node->block->GetCenterNode(LocalCenterId)) continue;
          double * PtrTemp =  ((double *) (block->GetCenterNode(LocalCenterId)->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset));
          PtrTemp[netChargeNewIndex] = 0.0; 
          PtrTemp[divEIndex]=0.0;// to compare with ipic3d
        }// for i
      }//for j   
    }//for k
    
  }//for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++)

}
  




//void UpdateBWrong(){ 
void PIC::FieldSolver::Electromagnetic::ECSIM::UpdateB(){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;  
  
  //compute B^(n+1) from B^(n) and E^(n+theta)
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }
    
    if (node->Thread!=PIC::ThisThread) continue;
     
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    double dx[3],coeff[3],coeff4[3],x[3]; 
    for (int iDim=0; iDim<3; iDim++){
      dx[iDim]=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim];
      dx[iDim]*=length_conv;
     
      coeff[iDim] = cDt/dx[iDim];
      coeff4[iDim] = coeff[iDim]*0.25; //coefficients for curl calculation
    }
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      char * offset;
      double Ex[2][2][2], Ey[2][2][2], Ez[2][2][2];
      //	  int index[3]={i,j,k};

      for (int kk=0;kk<2;kk++) for (int jj=0;jj<2;jj++) for (int ii=0;ii<2;ii++){
        offset=node->block->GetCornerNode(_getCornerNodeLocalNumber(i+ii,j+jj,k+kk))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
        double * ptr =  (double*)(offset+OffsetE_HalfTimeStep);
        Ex[ii][jj][kk]=ptr[ExOffsetIndex]*E_conv;
        Ey[ii][jj][kk]=ptr[EyOffsetIndex]*E_conv;
        Ez[ii][jj][kk]=ptr[EzOffsetIndex]*E_conv;
      }

      offset=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

      double * CurrentPtr = (double*)(offset+CurrentBOffset);
      double * PrevPtr = (double*)(offset+PrevBOffset);
      //store next B at prevptr
      double tempB[3]={0.0,0.0,0.0};

      for (int ii=0;ii<2;ii++){
        for (int jj=0; jj<2;jj++){
          tempB[BxOffsetIndex] += (-coeff4[1]*(Ez[ii][1][jj]-Ez[ii][0][jj])+coeff4[2]*(Ey[ii][jj][1]-Ey[ii][jj][0]));
          tempB[ByOffsetIndex] += (-coeff4[2]*(Ex[ii][jj][1]-Ex[ii][jj][0])+coeff4[0]*(Ez[1][ii][jj]-Ez[0][ii][jj]));
          tempB[BzOffsetIndex] += (-coeff4[0]*(Ey[1][ii][jj]-Ey[0][ii][jj])+coeff4[1]*(Ex[ii][1][jj]-Ex[ii][0][jj]));
        }
      }

      PrevPtr[BxOffsetIndex] = CurrentPtr[BxOffsetIndex]+tempB[BxOffsetIndex]/B_conv;
      PrevPtr[ByOffsetIndex] = CurrentPtr[ByOffsetIndex]+tempB[ByOffsetIndex]/B_conv;
      PrevPtr[BzOffsetIndex] = CurrentPtr[BzOffsetIndex]+tempB[BzOffsetIndex]/B_conv;
    }
  }


  //swap current and prev pointer
  int tempInt;
  tempInt=PrevBOffset;
  PrevBOffset=CurrentBOffset;
  CurrentBOffset=tempInt;
 

  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
    PIC::Mesh::mesh.ParallelBlockDataExchange(PackBlockData_B,UnpackBlockData_B);
    break;
      
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData(PackBlockData_B,UnpackBlockData_B);
    break;
  }
}

void PIC::FieldSolver::Electromagnetic::ECSIM::InterpolateB_C2N() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (!node->block) continue;
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
          //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
          BoundaryBlock=true;
          break;
        }
      
      if (BoundaryBlock==true) continue;
    }
    
    if (node->Thread!=PIC::ThisThread) continue;
    
    for (int k=0;k<=_BLOCK_CELLS_Z_;k++) for (int j=0;j<=_BLOCK_CELLS_Y_;j++) for (int i=0;i<=_BLOCK_CELLS_X_;i++) {
          char * offset=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+OffsetB_corner;
          
          double tempB_curr[3]={0,0,0}, tempB_prev[3]={0,0,0};
          for (int kk=-1;kk<=0; kk++) for (int jj=-1; jj<=0; jj++) for (int ii=-1; ii<=0; ii++){
            char *  offsetTmp=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+ii,j+jj,k+kk))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
            double * CurrentPtr = (double*)(offsetTmp+CurrentBOffset);
            //   double * PrevPtr = (double*)(offsetTmp+PrevBOffset);
            
            for (int idim=0;idim<3;idim++) {
              tempB_curr[idim] += CurrentPtr[idim];
              // tempB_prev[idim] += PrevPtr[idim];
            }
              }

          double * nodeCurrPtr = (double*)(offset+CurrentBOffset);
          // double * nodePrevPtr = (double*)(offset+PrevBOffset);
          
          for (int idim=0;idim<3;idim++) {
            nodeCurrPtr[idim] = tempB_curr[idim]/8.0; 
            // nodePrevPtr[idim] = tempB_prev[idim]/8.0;
          }
          
        }
  }
  
}


void PIC::FieldSolver::Electromagnetic::ECSIM::InterpolateB_N2C() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (!node->block) continue;
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
          //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
          BoundaryBlock=true;
          break;
        }
      
      if (BoundaryBlock==true) continue;
    }
    
    if (node->Thread!=PIC::ThisThread) continue;
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          //  char * offset=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+OffsetB_corner;
          char * offset=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

          double tempB_curr[3]={0,0,0}, tempB_prev[3]={0,0,0};
          for (int kk=0;kk<=1; kk++) for (int jj=0; jj<=1; jj++) for (int ii=0; ii<=1; ii++){
                //char *  offsetTmp=node->block->GetCenterNode(_getCenterNodeLocalNumber(i+ii,j+jj,k+kk))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
            char *  offsetTmp=node->block->GetCornerNode(_getCornerNodeLocalNumber(i+ii,j+jj,k+kk))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+OffsetB_corner;
            
            double * CurrentPtr = (double*)(offsetTmp+CurrentBOffset);
            double * PrevPtr = (double*)(offsetTmp+PrevBOffset);
            
            for (int idim=0;idim<3;idim++) {
              tempB_curr[idim] += CurrentPtr[idim];
              tempB_prev[idim] += PrevPtr[idim];
            }
              }

          double * nodeCurrPtr = (double*)(offset+CurrentBOffset);
          double * nodePrevPtr = (double*)(offset+PrevBOffset);
          
          for (int idim=0;idim<3;idim++) {
            nodeCurrPtr[idim] = tempB_curr[idim]/8.0; 
            nodePrevPtr[idim] = tempB_prev[idim]/8.0;
          }
          /*
          printf("interpolationB localnumber:%d, ind:%d,%d,%d,prev:%e,%e,%e,curr:%e,%e,%e\n",
                 _getCornerNodeLocalNumber(i,j,k),i,j,k,nodePrevPtr[0],nodePrevPtr[1],nodePrevPtr[2],
                 nodeCurrPtr[0],nodeCurrPtr[1],nodeCurrPtr[2]);
          */
        }
  }


}




void PIC::FieldSolver::Electromagnetic::ECSIM::UpdateE() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    
  //compute E^(n+1)  from E^(n+theta) and E^n
  
  double WaveEnergySum =0.0;

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }
      
      if (BoundaryBlock==true) continue;
    }
      
    if (node->Thread!=PIC::ThisThread) continue;
    
    double CellVolume=1.0;
    int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    for (int iDim=0; iDim<3;iDim++) CellVolume*=(node->xmax[iDim]-node->xmin[iDim])/nCell[iDim]*length_conv; 
    
    for (int k=0;k<_BLOCK_CELLS_Z_;k++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      char * offset;

      offset=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
      char * centerOffset =node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+ PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
      double Bx0,By0,Bz0;
      double * CurrentB_Ptr =  (double*)(centerOffset+CurrentBOffset);
      double * HalfStepPtr = (double*)(offset+OffsetE_HalfTimeStep);
      double * CurrentPtr = (double*)(offset+CurrentEOffset);
      double Ex,Ey,Ez;

      Ex = (HalfStepPtr[ExOffsetIndex]-(1.0-theta)*CurrentPtr[ExOffsetIndex])/theta;
      Ey = (HalfStepPtr[EyOffsetIndex]-(1.0-theta)*CurrentPtr[EyOffsetIndex])/theta;
      Ez = (HalfStepPtr[EzOffsetIndex]-(1.0-theta)*CurrentPtr[EzOffsetIndex])/theta;

      CurrentPtr[ExOffsetIndex] = Ex;
      CurrentPtr[EyOffsetIndex] = Ey;
      CurrentPtr[EzOffsetIndex] = Ez;

      Bx0=CurrentB_Ptr[BxOffsetIndex];
      By0=CurrentB_Ptr[ByOffsetIndex];
      Bz0=CurrentB_Ptr[BzOffsetIndex];
      WaveEnergySum += ((Ex*Ex+Ey*Ey+Ez*Ez)*E_conv*E_conv+(Bx0*Bx0+By0*By0+Bz0*Bz0)*B_conv*B_conv)*0.125/Pi*CellVolume;

    }
  }
    
  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
    PIC::Mesh::mesh.ParallelBlockDataExchange(PackBlockData_E,UnpackBlockData_E);
    break;
    
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData(PackBlockData_E,UnpackBlockData_E);
    break;
  }
 
  MPI_Reduce(&WaveEnergySum, &TotalWaveEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    printf("Total Wave Energy:%f\n",TotalWaveEnergy);
    // printf("Total Energy:%f\n",TotalParticleEnergy+TotalWaveEnergy);
  }
  
}
 

void PIC::FieldSolver::Electromagnetic::ECSIM::UpdateMatrixElement(cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1>::cMatrixRow* row){
  double fourPiDtTheta=4*Pi*dtTotal*theta;

  for (int iElement=0; iElement<row->nNonZeroElements;iElement++){
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1>::cStencilElement* el=row->Elements+iElement;
    el->MatrixElementValue=el->MatrixElementParameterTable[0];
    if (el->MatrixElementSupportTable[0]!=NULL)
      el->MatrixElementValue+=*((double *)el->MatrixElementSupportTable[0])*fourPiDtTheta;
    
    //printf("iElement:%d,const:%f,matrixvalue:%f\n",iElement,el->MatrixElementParameterTable[0],el->MatrixElementValue);
  }  
  
}



double PIC::FieldSolver::Electromagnetic::ECSIM::PoissonUpdateRhs(int iVar,
			      cLinearSystemCenterNode<PIC::Mesh::cDataCenterNode,1,7,0,1,1,0>::cRhsSupportTable* RhsSupportTable_CornerNodes,int RhsSupportLength_CornerNodes,
			      cLinearSystemCenterNode<PIC::Mesh::cDataCenterNode,1,7,0,1,1,0>::cRhsSupportTable* RhsSupportTable_CenterNodes,int RhsSupportLength_CenterNodes) {

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  double res=0.0;

  double * CenterOffset = ((double*)(RhsSupportTable_CenterNodes[0].AssociatedDataPointer+
             PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset));
  double gammaTmp =0.51;
  res = (-4*Pi*(gammaTmp*CenterOffset[netChargeNewIndex]+
               (1-gammaTmp)*CenterOffset[netChargeOldIndex])
         +CenterOffset[divEIndex])/gammaTmp;

  //the equation solves phi/gammaTmp
  return res;
}
 

 //update the RHS vector
double PIC::FieldSolver::Electromagnetic::ECSIM::UpdateRhs(int iVar,
			      cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1>::cRhsSupportTable* RhsSupportTable_CornerNodes,int RhsSupportLength_CornerNodes,
			      cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,81,82,16,1,1>::cRhsSupportTable* RhsSupportTable_CenterNodes,int RhsSupportLength_CenterNodes) {
  int i;
  double res=0.0;

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  double fourPiDtTheta=4*Pi*dtTotal*theta;
   
 
  for (int ii=0;ii<27;ii++) {
    double * tempPtr = (double*)(RhsSupportTable_CornerNodes[ii].AssociatedDataPointer+CurrentEOffset);

    res+=(tempPtr[ExOffsetIndex]*RhsSupportTable_CornerNodes[ii].Coefficient+
      tempPtr[EyOffsetIndex]*RhsSupportTable_CornerNodes[ii+27].Coefficient+
      tempPtr[EzOffsetIndex]*RhsSupportTable_CornerNodes[ii+54].Coefficient)*E_conv;
   }

   double * tempMassMatrixPtr = ((double*)RhsSupportTable_CornerNodes[0].AssociatedDataPointer)+MassMatrixOffsetIndex;
    
  //mass matrix part
  for (int ii=0;ii<27;ii++) {
    double * tempPtr = (double*)(RhsSupportTable_CornerNodes[ii].AssociatedDataPointer+CurrentEOffset);

    res+=(tempPtr[ExOffsetIndex]*tempMassMatrixPtr[MassMatrixOffsetTable[iVar][ii]]+
      tempPtr[EyOffsetIndex]*tempMassMatrixPtr[MassMatrixOffsetTable[iVar][ii+27]]+
      tempPtr[EzOffsetIndex]*tempMassMatrixPtr[MassMatrixOffsetTable[iVar][ii+54]])*(-fourPiDtTheta)*E_conv;
  }
    
    
  // current effect
  res+=((double*)(RhsSupportTable_CornerNodes[81].AssociatedDataPointer))[JxOffsetIndex+iVar]*
    RhsSupportTable_CornerNodes[81].Coefficient;
    
  //contribution from center nodes
  for (i=0; i<8;i++){
    res+=((double*)(RhsSupportTable_CenterNodes[i].AssociatedDataPointer+CurrentBOffset))[(iVar+2)%3]*RhsSupportTable_CenterNodes[i].Coefficient*B_conv;
  }//E=iVar,B=((iVar+2)%3) Ex:Bz, Ey:Bx, Ez:By
    
  for (i=8; i<16;i++){
    res+=((double*)(RhsSupportTable_CenterNodes[i].AssociatedDataPointer+CurrentBOffset))[(iVar+4)%3]*RhsSupportTable_CenterNodes[i].Coefficient*B_conv;
  }//E=iVar,B=((iVar+4)%3)  Ex:By, Ey:Bz, Ez:Bx
    
  //if (fabs(res)>1e-3) printf("rhs:%f\n",res);

  return res;
}





void PIC::FieldSolver::Electromagnetic::ECSIM::BuildMatrix() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    
  Solver.Reset();
  Solver.BuildMatrix(GetStencil);
  if (DoDivECorrection){
    PoissonSolver.Reset();
    PoissonSolver.BuildMatrix(PoissonGetStencil);
  }
}

void PIC::FieldSolver::Electromagnetic::ECSIM::TimeStep() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    
  
  //perform the rest of the field solver calculstions
  double t0,t1,StartTime=MPI_Wtime();
  
  if (PIC::CPLR::FLUID::iCycle==0){  
    UpdateJMassMatrix();    
    ComputeNetCharge(true);
    SetBoundaryChargeDivE();

    {// Output
      double timeNow = 0.0;  
      PIC::CPLR::FLUID::write_output(timeNow);
    }    
  }
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ ){
    setB_center_BC();
    setB_corner_BC();
    setE_curr_BC();
  }
  
  //  PIC::BC::ExternalBoundary::UpdateData();

  Solver.UpdateRhs(UpdateRhs); 
  Solver.UpdateMatrixNonZeroCoefficients(UpdateMatrixElement);

  t0=MPI_Wtime();
  linear_solver_matvec_c = matvec;
  if (PIC::ThisThread==0) printf("---------------Solving E field-----------\n");
  Solver.Solve(SetInitialGuess,ProcessFinalSolution,PIC::CPLR::FLUID::EFieldTol,
	       PIC::CPLR::FLUID::EFieldIter,PackBlockData_E,UnpackBlockData_E);
  t1=MPI_Wtime();
  CumulativeTiming::SolveTime+=t1-t0;

  t0=t1;
  
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ )
    setE_half_BC();
  UpdateB();
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ )
    setB_center_BC();
  PIC::Mesh::mesh.ParallelBlockDataExchange(PackBlockData_B,UnpackBlockData_B);
  InterpolateB_C2N();
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ )
    setB_corner_BC();
  
  t1=MPI_Wtime();
  CumulativeTiming::UpdateBTime+=t1-t0;

  t0=t1;
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ )
    setE_half_BC();
  UpdateE();
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_ )
    setE_curr_BC();
  t1=MPI_Wtime();
  CumulativeTiming::UpdateETime+=t1-t0;
  CumulativeTiming::TotalRunTime+=MPI_Wtime()-StartTime;
 }


//set the initial guess
void PIC::FieldSolver::Electromagnetic::ECSIM::SetInitialGuess(double* x,PIC::Mesh::cDataCornerNode* CornerNode) {
  //  x[0]=*((double*)(CornerNode->GetAssociatedDataBufferPointer()+CurrentCornerNodeOffset));
  x[0]=0.0;
  x[1]=0.0;
  x[2]=0.0;
}


void PIC::FieldSolver::Electromagnetic::ECSIM::PoissonSetInitialGuess(double* x,PIC::Mesh::cDataCenterNode* CenterNode) {
  //  x[0]=*((double*)(CornerNode->GetAssociatedDataBufferPointer()+CurrentCornerNodeOffset));
  x[0]=0.0;

}

//process the solution vector
void PIC::FieldSolver::Electromagnetic::ECSIM::ProcessFinalSolution(double* x,PIC::Mesh::cDataCornerNode* CornerNode) {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    
  char *offset=CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

  ((double*)(offset+OffsetE_HalfTimeStep))[0] =x[0]/E_conv+((double*)(offset+CurrentEOffset))[0];
  ((double*)(offset+OffsetE_HalfTimeStep))[1] =x[1]/E_conv+((double*)(offset+CurrentEOffset))[1];
  ((double*)(offset+OffsetE_HalfTimeStep))[2] =x[2]/E_conv+((double*)(offset+CurrentEOffset))[2];

}


void PIC::FieldSolver::Electromagnetic::ECSIM::PoissonProcessFinalSolution(double* x,PIC::Mesh::cDataCenterNode* CenterNode) {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    
  double *offset=(double *)(CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset);

  offset[phiIndex] = x[0];
 
}


void PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCenterNodeVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"Bx (center node)\",\"By (center node)\",\"Bz (center node)\"");
}

void PIC::FieldSolver::Electromagnetic::ECSIM::output::InterpolateCenterNode(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {

  double Wave[3];
  int i,iDim;
  char *SamplingBuffer;
  
  for (iDim =0;iDim<3; iDim++) Wave[iDim]=0.0;

  for (i=0;i<nInterpolationCoeficients;i++) {
   SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+CurrentBOffset;

   // SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
    Wave[0] += (((double*)SamplingBuffer)[0])*InterpolationCoeficients[i];
    Wave[1] += (((double*)SamplingBuffer)[1])*InterpolationCoeficients[i];
    Wave[2] += (((double*)SamplingBuffer)[2])*InterpolationCoeficients[i];
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+CurrentBOffset,Wave,3*sizeof(double));
}

void PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCenterNodeData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;
  double * t;

  if (pipe->ThisThread==CenterNodeThread) {
    t= (double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+CurrentBOffset);
    // t= (double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset);
  }
  
  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) t=pipe->recvPointer<double>(3,CenterNodeThread);
    fprintf(fout,"%e %e %e ",t[0],t[1],t[2]);
  }
  else pipe->send(t,3);

}

void PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCornerNodeVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"Ex (corner node)\",\"Ey (corner node)\",\"Ez (corner node)\"");
}

void PIC::FieldSolver::Electromagnetic::ECSIM::output::PrintCornerNodeData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CornerNodeThread,PIC::Mesh::cDataCornerNode *CornerNode) {
  int idim;
  double * t;
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;    

  if (pipe->ThisThread==CornerNodeThread) {
    t= ((double*)(CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+CurrentEOffset));
  }

  if (pipe->ThisThread==0) {
    if (CornerNodeThread!=0) t=pipe->recvPointer<double>(3,CornerNodeThread);
    fprintf(fout,"%e %e %e ",t[0],t[1],t[2]);
  }
  else pipe->send(t,3);
}

void PIC::FieldSolver::Electromagnetic::ECSIM::matvec(double* VecIn, double * VecOut, int n){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  Solver.MultiplyVector(VecOut,VecIn,n);
}


void PIC::FieldSolver::Electromagnetic::ECSIM::PoissonMatvec(double* VecIn, double * VecOut, int n){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  PoissonSolver.MultiplyVector(VecOut,VecIn,n);
}

int isFaceBoundary(int sum, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){

  switch (sum) {

  case   1:
    if (!node->GetNeibFace(0,0,0)){
      return 1;
    }else if (node->GetNeibFace(0,0,0)->Thread==-1){
      return 1;
    }
    else return 0;
    
  case   1000:
    if (!node->GetNeibFace(1,0,0)){
      return 1;
    }else if (node->GetNeibFace(1,0,0)->Thread==-1){
      return 1;
    }
    else return 0;
    
  case   10:
    if (!node->GetNeibFace(2,0,0)){
      return 2;
    }else if (node->GetNeibFace(2,0,0)->Thread==-1){
      return 2;
    }
    else return 0;

  case 10000:
    if (!node->GetNeibFace(3,0,0)){
      return 2;
    }else if (node->GetNeibFace(3,0,0)->Thread==-1){
      return 2;
    }
    else return 0;

  case  100:
    if (!node->GetNeibFace(4,0,0)){
      return 4;
    }else if (node->GetNeibFace(4,0,0)->Thread==-1){
      return 4;
    }
    else return 0;
    
  case  100000:
    if (!node->GetNeibFace(5,0,0)){
      return 4;
    }else if (node->GetNeibFace(5,0,0)->Thread==-1){
      return 4;
    }
    else return 0;
    
  default:
    return 0;
  }

}

int isEdgeBoundary(int sum, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){
  
  int temp,a,b;
  switch (sum){
    
  case 110:
    a = isFaceBoundary(100,node);
    b = isFaceBoundary(10,node);
    if (a && b){
      if (a!=b){
        return 6;
      }else{
        return a;
      }
    }
    
    if (a+b) { return a+b;
    }else {
      if (!node->GetNeibEdge(0,0)) return 6;
      else if (node->GetNeibEdge(0,0)->Thread==-1) return 6;
      else return 0;
    }
    
  case 10100:
    a= isFaceBoundary(10000,node);
    b= isFaceBoundary(100,node);
    if (a && b){
      if (a!=b){
        return 6;
      }else{
        return a;
      }
    }
    
    if (a+b) { return a+b;
    }else {
      if (!node->GetNeibEdge(1,0)) return 6;
      else if (node->GetNeibEdge(1,0)->Thread==-1) return 6;
      else return 0;
    }
    
  case 110000:
    a = isFaceBoundary(100000,node);
    b = isFaceBoundary(10000,node);

    if (a && b){
      if (a!=b){
        return 6;
      }else{
        return a;
      }
    }
    
    if (a+b) { return a+b;
    } else {
      if (!node->GetNeibEdge(2,0)) return 6;
      else if (node->GetNeibEdge(2,0)->Thread==-1) return 6;
      else return 0;
    }

  case 100010:
    a = isFaceBoundary(100000,node);
    b = isFaceBoundary(10,node);
    if (a && b){
      if (a!=b) {
        return 6;
      }else{
        return a;
      }
    }
    
    if (a+b) { return a+b;
    }else{
      if (!node->GetNeibEdge(3,0)) return 6;
      else if (node->GetNeibEdge(3,0)->Thread==-1) return 6;
      else return 0;
    }

  case 101:
    a = isFaceBoundary(100,node);
    b = isFaceBoundary(1,node);
    if (a && b){
      if (a!=b){
        return 5;
      }else{
        return a;
      }
    }
      
    if (a+b) { return a+b;
    }else{
      if (!node->GetNeibEdge(4,0)) return 5;
      else if (node->GetNeibEdge(4,0)->Thread==-1) return 5;
      else return 0;
    }
    
  case 1100:
    a = isFaceBoundary(1000,node);
    b = isFaceBoundary(100,node);

    if (a && b){
      if (a!=b){
        return 5;
      }else{
        return a;
      }
    }
   
    if (a+b){ return a+b;
    }else{
      if (!node->GetNeibEdge(5,0)) return 5;
      else if (node->GetNeibEdge(5,0)->Thread==-1) return 5;
      else return 0;
    }
  
  case 101000:
    a = isFaceBoundary(100000,node);
    b = isFaceBoundary(1000,node);

    if (a && b){
      if (a!=b){
        return 5;
      }else{
        return a;
      }
    }
    
    if (a+b){ return a+b;
    }else{
      if (!node->GetNeibEdge(6,0)) return 5;
      else if (node->GetNeibEdge(6,0)->Thread==-1) return 5;
      else return 0;
    }
    
  case 100001:
    a = isFaceBoundary(100000,node);
    b = isFaceBoundary(1,node);
    
    if (a && b){
      if (a!=b){
        return 5;
      }else{
        return a;
      }
    }

    if (a+b){ return a+b;
    }else{
      if (!node->GetNeibEdge(7,0)) return 5;
      else if (node->GetNeibEdge(7,0)->Thread==-1) return 5;
      else return 0;
    }
    
  case 11:
    a = isFaceBoundary(10,node);
    b = isFaceBoundary(1,node);
  
    if (a && b){
      if (a!=b){
        return 3;
      }else{
        return a;
      }
    }
    
    if (a+b){ return a+b;
    }else{
      if (!node->GetNeibEdge(8,0)) return 3;
      else if (node->GetNeibEdge(8,0)->Thread==-1) return 3;
      else return 0;
    }
   
  case 1010:
    a = isFaceBoundary(1000,node);
    b = isFaceBoundary(10,node);
    
    if (a && b){
      if (a!=b){
        return 3;
      }else{
        return a;
      }
    }
    
    if (a+b){ return a+b;
    }else{
      if (!node->GetNeibEdge(9,0)) return 3;
      else if (node->GetNeibEdge(9,0)->Thread==-1) return 3;
      else return 0;
    }
    
  case 11000:
    a = isFaceBoundary(10000,node);
    b = isFaceBoundary(1000,node);
    
    if (a && b){
      if (a!=b){
        return 3;
      }else{
        return a;
      }
    }

    if (a+b){ return a+b;
    }else{
      if (!node->GetNeibEdge(10,0)) return 3;
      else if (node->GetNeibEdge(10,0)->Thread==-1) return 3;
      else return 0;
    }
        
  case 10001:
    a = isFaceBoundary(10000,node);
    b = isFaceBoundary(1,node);
   
    if (a && b){
      if (a!=b){
        return 3;
      }else{
        return a;
      }
    }
    
    if (a+b){ return a+b;
    }else{    
      if (!node->GetNeibEdge(11,0)) return 3;
      else if (node->GetNeibEdge(11,0)->Thread==-1) return 3;
      else return 0;
    }
    
  default:
    return 0;
  }
}

int computeRes(int *a){
  
  int iface[3]={1,2,4};
  int iedge[3]={6,5,3};
 
  for (int ii=0;ii<3; ii++){
    for (int jj=0;jj<3;jj++){
      if (iface[ii]==a[jj]) {
        a[3]+=iface[ii];
        break;
      }
    }
  }
  
  for (int ii=0;ii<3; ii++){
    for (int jj=0;jj<4;jj++){
      if (iedge[ii]==a[jj]) {
        a[4]+=iedge[ii];
        break;
      }
    }
  }

  if (a[4]>=7) {
    return 7;
  }else if (a[4]>0){
    return a[4];
  }else if (a[3]>0){
    return a[3];
  }
  
  return 0;
 
}

int isCornerBoundary(int sum, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){
 
  int temp,a[5],res;
  a[3]=0,a[4]=0;

  switch (sum){
    
  case 111:
    a[0]=isEdgeBoundary(11,node);
    a[1]=isEdgeBoundary(101,node);
    a[2]=isEdgeBoundary(110,node);

    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(0)) {
        return 7;
      }else if (node->GetNeibCorner(0)->Thread==-1){
        return 7;
      }
      else return 0;
    }

  case 1110:
    a[0]=isEdgeBoundary(110,node);
    a[1]=isEdgeBoundary(1010,node);
    a[2]=isEdgeBoundary(1100,node);
    
    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(1)) {
        return 7;
      }else if (node->GetNeibCorner(1)->Thread==-1){
        return 7;
      }
      else return 0;
    }
    
  case 10101:
    a[0]=isEdgeBoundary(101,node);
    a[1]=isEdgeBoundary(10001,node);
    a[2]=isEdgeBoundary(10100,node);
    
    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(2)) {
        return 7;
      }else if (node->GetNeibCorner(2)->Thread==-1){
        return 7;
      }
      else return 0;
    }
       
  case 11100:
    a[0]=isEdgeBoundary(1100,node);
    a[1]=isEdgeBoundary(10100,node);
    a[2]=isEdgeBoundary(11000,node);
    
    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(3)) {
        return 7;
      }else if (node->GetNeibCorner(3)->Thread==-1){
        return 7;
      }
      else return 0;
    }
    
  case 100011:
    a[0]=isEdgeBoundary(11,node);
    a[1]=isEdgeBoundary(100010,node);
    a[2]=isEdgeBoundary(100001,node);

    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(4)) {
        return 7;
      }else if (node->GetNeibCorner(4)->Thread==-1){
        return 7;
      }
      else return 0;
    }
    
  case 101010:
    a[0]=isEdgeBoundary(1010,node);
    a[1]=isEdgeBoundary(101000,node);
    a[2]=isEdgeBoundary(100010,node);

    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(5)) {
        return 7;
      }else if (node->GetNeibCorner(5)->Thread==-1){
        return 7;
      }
      else return 0;
    }

  case 110001:
    a[0]=isEdgeBoundary(10001,node);
    a[1]=isEdgeBoundary(100001,node);
    a[2]=isEdgeBoundary(110000,node);

    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(6)) {
        return 7;
      }else if (node->GetNeibCorner(6)->Thread==-1){
        return 7;
      }
      else return 0;
    }
    
  case 111000:
    a[0]=isEdgeBoundary(11000,node);
    a[1]=isEdgeBoundary(101000,node);
    a[2]=isEdgeBoundary(110000,node);

    res = computeRes(a);
    if (res>0){
      return res;
    }else{
      if (!node->GetNeibCorner(7)) {
        return 7;
      }else if (node->GetNeibCorner(7)->Thread==-1){
        return 7;
      }
      else return 0;
    }

  default:
    return 0;
  }
}


int  PIC::FieldSolver::Electromagnetic::ECSIM::isBoundaryCell(double * x, double *dx, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){
  // 0:    not boundary cell
  // 1:001 x-direction face  
  // 2:010 y-direction face
  // 4:100 z-direction face
  // 6:110 x-direction edge
  // 5:101 y-direction edge
  // 3:011 z-direction edge
  // 7:111 corner
  // 8: outside block and outside domain
  /*
  bool isTest=false;
  if (fabs(x[0]-16.5)<0.1 && fabs(x[1]-8.5)<0.1 && fabs(x[2]-4.5)<0.1) isTest=true;
  */
  for (int idim=0; idim<3 && node; idim++){
    if (x[idim]<(node->xmin[idim]-PIC::Mesh::mesh.EPS)) node=node->GetNeibFace(idim*2,0,0);
  }

  for (int idim=0; idim<3 && node; idim++){
    if (x[idim]>(node->xmax[idim]+PIC::Mesh::mesh.EPS)) node=node->GetNeibFace(idim*2+1,0,0);
  }

  
  if (node==NULL || node->Thread==-1) return 8;


  int addition =1, sum=0;//sum used to indicate the location of the corner
  //0 not at the boundary, 111000 at the right most corner...
  for (int idim=0;idim<3;idim++) {
    if (fabs(x[idim]-0.5*dx[idim]-node->xmin[idim])<PIC::Mesh::mesh.EPS) sum+=addition;
    addition *=10;
  }
  
  for (int idim=0;idim<3;idim++) {
    if (fabs(x[idim]+0.5*dx[idim]-node->xmax[idim])<PIC::Mesh::mesh.EPS) sum+=addition;
    addition *=10;
  }

  if (sum==0) return 0;

  int val =  isCornerBoundary(sum,node)+isEdgeBoundary(sum,node)+isFaceBoundary(sum,node);

  return val;
  
}


bool PIC::FieldSolver::Electromagnetic::ECSIM::isBoundaryCorner(double * x, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){
  
  for (int idim=0; idim<3 && node; idim++){
    if (x[idim]<(node->xmin[idim]-PIC::Mesh::mesh.EPS)) node=node->GetNeibFace(idim*2,0,0);
  }

  for (int idim=0; idim<3 && node; idim++){
    if (x[idim]>(node->xmax[idim]+PIC::Mesh::mesh.EPS)) node=node->GetNeibFace(idim*2+1,0,0);
  }

  
  if (node==NULL || node->Thread==-1) return true;

  for (int idim=0;idim<3;idim++) if (node->GetNeibFace(idim*2,0,0)==NULL || node->GetNeibFace(idim*2,0,0)->Thread==-1) {

      if (fabs(x[idim]-node->xmin[idim])<PIC::Mesh::mesh.EPS){
        return true;
      }
    }
  
  for (int idim=0;idim<3;idim++) if (node->GetNeibFace(idim*2+1,0,0)==NULL || node->GetNeibFace(idim*2+1,0,0)->Thread==-1) {
      if (fabs(x[idim]-node->xmax[idim])<PIC::Mesh::mesh.EPS)
        return true;
    }

  int addition =1, sum=0;//sum used to indicate the location of the corner
  //0 not at the boundary, 111000 at the right most corner...
  for (int idim=0;idim<3;idim++) {
    if (fabs(x[idim]-node->xmin[idim])<PIC::Mesh::mesh.EPS) sum+=addition;
    addition *=10;
  }
  
  for (int idim=0;idim<3;idim++) {
    if (fabs(x[idim]-node->xmax[idim])<PIC::Mesh::mesh.EPS) sum+=addition;
    addition *=10;
  }

  if (sum==0) return false;
  
  switch (sum) {
  case 111:
    if (!node->GetNeibCorner(0)) {
      return true;
    }else if (node->GetNeibCorner(0)->Thread==-1){
      return true;
    }
    else return false;
    
  case 1110:
    if (!node->GetNeibCorner(1)) {
      return true;
    }else if (node->GetNeibCorner(1)->Thread==-1){
      return true;
    }
    else return false;
    
  case 10101:
    if (!node->GetNeibCorner(2)) {
      return true;
    }else if (node->GetNeibCorner(2)->Thread==-1){
      return true;
    }
    else return false;
   
  case 11100:
    if (!node->GetNeibCorner(3)) {
      return true;
    }else if (node->GetNeibCorner(3)->Thread==-1){
      return true;
    }
    else return false;

  case 100011:
    if (!node->GetNeibCorner(4)) {
      return true;
    }else if (node->GetNeibCorner(4)->Thread==-1){
      return true;
    }
    else return false;

  case 101010:
    if (!node->GetNeibCorner(5)) {
      return true;
    }else if (node->GetNeibCorner(5)->Thread==-1){
      return true;
    }
    else return false;


  case 110001:
    if (!node->GetNeibCorner(6)) {
      return true;
    }else if (node->GetNeibCorner(6)->Thread==-1){
      return true;
    }
    else return false;


  case 111000:
    if (!node->GetNeibCorner(7)) {
      return true;
    }else if (node->GetNeibCorner(7)->Thread==-1){
      return true;
    }
    else return false;
   
  case 110:
    if (!node->GetNeibEdge(0,0)) return true;
    else if (node->GetNeibEdge(0,0)->Thread==-1) return true;
    else return false;

  case 10100:
    if (!node->GetNeibEdge(1,0)) return true;
    else if (node->GetNeibEdge(1,0)->Thread==-1) return true;
    else return false;

  case 110000:
    if (!node->GetNeibEdge(2,0)) return true;
    else if (node->GetNeibEdge(2,0)->Thread==-1) return true;
    else return false;

  case 100010:
    if (!node->GetNeibEdge(3,0)) return true;
    else if (node->GetNeibEdge(3,0)->Thread==-1) return true;
    else return false;

  case 101:
    if (!node->GetNeibEdge(4,0)) return true;
    else if (node->GetNeibEdge(4,0)->Thread==-1) return true;
    else return false;

  case 1100:
    if (!node->GetNeibEdge(5,0)) return true;
    else if (node->GetNeibEdge(5,0)->Thread==-1) return true;
    else return false;

  case 101000:
    if (!node->GetNeibEdge(6,0)) return true;
    else if (node->GetNeibEdge(6,0)->Thread==-1) return true;
    else return false;

  case 100001:
    if (!node->GetNeibEdge(7,0)) return true;
    else if (node->GetNeibEdge(7,0)->Thread==-1) return true;
    else return false;

  case 11:
    if (!node->GetNeibEdge(8,0)) return true;
    else if (node->GetNeibEdge(8,0)->Thread==-1) return true;
    else return false;

  case 1010:
    if (!node->GetNeibEdge(9,0)) return true;
    else if (node->GetNeibEdge(9,0)->Thread==-1) return true;
    else return false;

  case 11000:
    if (!node->GetNeibEdge(10,0)) return true;
    else if (node->GetNeibEdge(10,0)->Thread==-1) return true;
    else return false;

  case 10001:
    if (!node->GetNeibEdge(11,0)) return true;
    else if (node->GetNeibEdge(11,0)->Thread==-1) return true;
    else return false;
    
  default:
    return false;
  }
  

}




bool PIC::FieldSolver::Electromagnetic::ECSIM::isRightBoundaryCorner(double * x, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){

  for (int idim=0;idim<3;idim++) if (node->GetNeibFace(idim*2+1,0,0)==NULL || node->GetNeibFace(idim*2+1,0,0)->Thread==-1) {
      if (fabs(x[idim]-node->xmax[idim])<PIC::Mesh::mesh.EPS)
        return true;
    }  
  return false;
  
}

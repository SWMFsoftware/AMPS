//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//the interface between AMPS and SWMF

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "pic.h"

//using namespace std;

int PIC::CPLR::SWMF::MagneticFieldOffset=-1;
int PIC::CPLR::SWMF::PlasmaNumberDensityOffset=-1;
int PIC::CPLR::SWMF::BulkVelocityOffset=-1;
int PIC::CPLR::SWMF::PlasmaPressureOffset=-1;
int PIC::CPLR::SWMF::PlasmaTemperatureOffset=-1;
int PIC::CPLR::SWMF::AlfvenWaveI01Offset=-1;

int PIC::CPLR::SWMF::MagneticFieldOffset_last=-1;
int PIC::CPLR::SWMF::PlasmaNumberDensityOffset_last=-1;
int PIC::CPLR::SWMF::BulkVelocityOffset_last=-1;
int PIC::CPLR::SWMF::PlasmaPressureOffset_last=-1;
int PIC::CPLR::SWMF::PlasmaTemperatureOffset_last=-1;
int PIC::CPLR::SWMF::AlfvenWaveI01Offset_last=-1;


bool PIC::CPLR::SWMF::OhCouplingFlag=false;
bool PIC::CPLR::SWMF::IhCouplingFlag=false;
bool PIC::CPLR::SWMF::BlCouplingFlag=false;

int PIC::CPLR::SWMF::TotalDataLength=0;
double PIC::CPLR::SWMF::MeanPlasmaAtomicMass=1.0*_AMU_;
bool PIC::CPLR::SWMF::FirstCouplingOccured=false;
list<PIC::CPLR::SWMF::fSendCenterPointData> PIC::CPLR::SWMF::SendCenterPointData;
int PIC::CPLR::SWMF::nCommunicatedIonFluids=1;


//the SWMF simulation time when the last two couplings have occured
double PIC::CPLR::SWMF::CouplingTime=-1.0; 
double PIC::CPLR::SWMF::CouplingTime_last=-1.0; 

//set the interpolation stencil that is used for interpolation in the coupler
void PIC::CPLR::InitInterpolationStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  switch( _PIC_COUPLER__INTERPOLATION_MODE_) {
  case _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_:
    PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node);
    break;
  case  _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_:
    PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node);
    break;
  case _PIC_COUPLER__INTERPOLATION_MODE__CORNER_BASED_LINEAR_:
    PIC::InterpolationRoutines::CornerBased::InitStencil(x,node);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }
}


int PIC::CPLR::SWMF::RequestDataBuffer(int offset) {
  MagneticFieldOffset=offset;
  TotalDataLength=3;

  BulkVelocityOffset=offset+TotalDataLength*sizeof(double);
  TotalDataLength+=3*nCommunicatedIonFluids;

  PlasmaPressureOffset=offset+TotalDataLength*sizeof(double);
  TotalDataLength+=nCommunicatedIonFluids;

  PlasmaNumberDensityOffset=offset+TotalDataLength*sizeof(double);
  TotalDataLength+=nCommunicatedIonFluids;

  PlasmaTemperatureOffset=offset+TotalDataLength*sizeof(double);
  TotalDataLength+=nCommunicatedIonFluids;

  if (IhCouplingFlag==true) {
    AlfvenWaveI01Offset=offset+TotalDataLength*sizeof(double);
    TotalDataLength+=2;
  }

  //keep two data sets exported from the SWMF for calculating time-derivatives
  if (_PIC_SWMF_COUPLER__SAVE_TWO_DATA_SETS_== _PIC_MODE_ON_) {
    MagneticFieldOffset_last=offset+TotalDataLength*sizeof(double);
    TotalDataLength+=3;

    BulkVelocityOffset_last=offset+TotalDataLength*sizeof(double);
    TotalDataLength+=3*nCommunicatedIonFluids;

    PlasmaPressureOffset_last=offset+TotalDataLength*sizeof(double);
    TotalDataLength+=nCommunicatedIonFluids;

    PlasmaNumberDensityOffset_last=offset+TotalDataLength*sizeof(double);
    TotalDataLength+=nCommunicatedIonFluids;

    PlasmaTemperatureOffset_last=offset+TotalDataLength*sizeof(double);
    TotalDataLength+=nCommunicatedIonFluids;

    if (IhCouplingFlag==true) {
      AlfvenWaveI01Offset_last=offset+TotalDataLength*sizeof(double);
      TotalDataLength+=2;
    }
  }
  else {
    MagneticFieldOffset_last=MagneticFieldOffset;
    BulkVelocityOffset_last=BulkVelocityOffset;

    PlasmaPressureOffset_last=PlasmaPressureOffset;
    PlasmaNumberDensityOffset_last=PlasmaNumberDensityOffset;

    PlasmaTemperatureOffset_last=PlasmaTemperatureOffset;

    if (IhCouplingFlag==true) {
      AlfvenWaveI01Offset_last=AlfvenWaveI01Offset;
    }
  } 
     


  return TotalDataLength*sizeof(double);
}

void PIC::CPLR::SWMF::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"gmN\",\"gmP\",\"gmVx\",\"gmVy\",\"gmVz\",\"gmBx\",\"gmBy\",\"gmBz\"");

  if (IhCouplingFlag==true) {
     fprintf(fout,",\"AlfvenWaveI01\", \"AlfvenWaveI02\"");
  } 

}

void PIC::CPLR::SWMF::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double B[3]={0.0,0.0,0.0},V[3]={0.0,0.0,0.0},P=0.0,Rho=0.0,i01=0.0,i02=0.0;
  int i,idim;
  char *SamplingBuffer;

  for (i=0;i<nInterpolationCoeficients;i++) {

    for (idim=0,SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+MagneticFieldOffset;idim<3;idim++) B[idim]+=(*((double*)(SamplingBuffer+idim*sizeof(double))))*InterpolationCoeficients[i];
    for (idim=0,SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+BulkVelocityOffset;idim<3;idim++) V[idim]+=(*((double*)(SamplingBuffer+idim*sizeof(double))))*InterpolationCoeficients[i];

    P+=(*((double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+PlasmaPressureOffset)))*InterpolationCoeficients[i];
    Rho+=(*((double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+PlasmaNumberDensityOffset)))*InterpolationCoeficients[i];

    if (IhCouplingFlag==true) {
      i01+=((double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+AlfvenWaveI01Offset))[0]*InterpolationCoeficients[i];
      i02+=((double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+AlfvenWaveI01Offset))[1]*InterpolationCoeficients[i];
    }
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+MagneticFieldOffset,B,3*sizeof(double));
  memcpy(CenterNode->GetAssociatedDataBufferPointer()+BulkVelocityOffset,V,3*sizeof(double));
  memcpy(CenterNode->GetAssociatedDataBufferPointer()+PlasmaPressureOffset,&P,sizeof(double));
  memcpy(CenterNode->GetAssociatedDataBufferPointer()+PlasmaNumberDensityOffset,&Rho,sizeof(double));

  if (IhCouplingFlag==true) {
    ((double*)(CenterNode->GetAssociatedDataBufferPointer()+AlfvenWaveI01Offset))[0]=i01;
    ((double*)(CenterNode->GetAssociatedDataBufferPointer()+AlfvenWaveI01Offset))[1]=i02;
  }
}

void PIC::CPLR::SWMF::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;
  double t;

  //Density
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+PlasmaNumberDensityOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

  //Pressure
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+PlasmaPressureOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

  //Bulk Velocity
  for (idim=0;idim<3;idim++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+BulkVelocityOffset+idim*sizeof(double)));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }

  //Magnetic Field
  for (idim=0;idim<3;idim++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+MagneticFieldOffset+idim*sizeof(double)));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }

  //AlfvenWave
  if (IhCouplingFlag==true) {
    double tt[2];

    if (pipe->ThisThread==CenterNodeThread) {
      for (int i=0;i<2;i++) tt[i]=((double*)(CenterNode->GetAssociatedDataBufferPointer()+AlfvenWaveI01Offset))[i];
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(tt,2,CenterNodeThread);

      fprintf(fout,"%e %e ",tt[0],tt[1]);
    }
    else pipe->send(tt,2);
  }

}


void PIC::CPLR::SWMF::init() {
  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestDataBuffer);

  //print out of the otuput file
  PIC::Mesh::AddVaraibleListFunction(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
}



void PIC::CPLR::SWMF::ConvertMpiCommunicatorFortran2C(signed int* iComm,signed int* iProc,signed int* nProc) {
  MPI_GLOBAL_COMMUNICATOR=MPI_Comm_f2c(*iComm);

 /* In case AMPS is used to trace particles along magnetic field lines, 
 *  it could be beneficial to run multiple copies of AMPS that are not connected via MPI. 
 *  For that, the communicator that AMPS receives from the SWMF is split such that the resulted communicator has only one MPI process.
 */

  if (_PIC_DISCONNECTED_MPI_PROCESSES_==_PIC_MODE_ON_) {
    //set the output directories
    char cmd[300]; 
    int rank;

    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);

    if (strcmp(PIC::OutputDataFileDirectory,".")!=0) {
      sprintf(PIC::OutputDataFileDirectory,"%s.thread=%ld",PIC::OutputDataFileDirectory,rank);
    }
    else {
      sprintf(PIC::OutputDataFileDirectory,"amps-out.thread=%ld",rank);

      sprintf(cmd,"mkdir -p %s",PIC::OutputDataFileDirectory);
      system(cmd);

      sprintf(cmd,"rm -rf %s/*",PIC::OutputDataFileDirectory);
      system(cmd);
    }

    //define the communicator
    MPI_GLOBAL_COMMUNICATOR=MPI_COMM_SELF;
  }


  PIC::InitMPI();

  if (PIC::ThisThread==0) {
    printf("AMPS: MPI Communicatior is imported from SWMF, size=%i\n",PIC::nTotalThreads);
  }
}

void PIC::CPLR::SWMF::ResetCenterPointProcessingFlag() {
  int thread,i,j,k;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  //init the cell processing flags
  for (node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL) {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
          if (cell!=NULL) cell->nodeDescriptor.nodeProcessedFlag=_OFF_AMR_MESH_;
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh->nTotalThreads;thread++) for (node=PIC::Mesh::mesh->DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL)  {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
          if (cell!=NULL) cell->nodeDescriptor.nodeProcessedFlag=_OFF_AMR_MESH_;
        }
    }
  }
}


void PIC::CPLR::SWMF::GetCenterPointNumber(int *nCenterPoints,fTestPointInsideDomain TestPointInsideDomain) {
  int thread,i,j,k;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  //init the cell processing flags
  ResetCenterPointProcessingFlag();
  *nCenterPoints=0;

  //count the number of the center points
  for (node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL) {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            bool inside_domain_flag;

            inside_domain_flag=(TestPointInsideDomain!=NULL) ? TestPointInsideDomain(cell->GetX()) : true;

            if (inside_domain_flag==true) {
              (*nCenterPoints)++;
              cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
            }
          }
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh->nTotalThreads;thread++) for (node=PIC::Mesh::mesh->DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL) {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            bool inside_domain_flag;

            inside_domain_flag=(TestPointInsideDomain!=NULL) ? TestPointInsideDomain(cell->GetX()) : true;

            if (inside_domain_flag==true) {
              (*nCenterPoints)++;
              cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
            }
          }
        }
    }
  }

}

void PIC::CPLR::SWMF::GetCenterPointCoordinates(double *x,fTestPointInsideDomain TestPointInsideDomain) {
  int thread,i,j,k,cnt=0;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  //init the cell processing flags
  ResetCenterPointProcessingFlag();

  //get coordinated of the center points
  for (node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL) {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            bool inside_domain_flag;

            inside_domain_flag=(TestPointInsideDomain!=NULL) ? TestPointInsideDomain(cell->GetX()) : true;

            if (inside_domain_flag==true) {
              cell->GetX(x+3*cnt);
              cnt++;
              cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
            }
          }
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh->nTotalThreads;thread++) for (node=PIC::Mesh::mesh->DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL) {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            bool inside_domain_flag;

            inside_domain_flag=(TestPointInsideDomain!=NULL) ? TestPointInsideDomain(cell->GetX()) : true;

            if (inside_domain_flag==true) {
              cell->GetX(x+3*cnt);
              cnt++;
              cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
            }
          }
        }
    }
  }
}

void PIC::CPLR::SWMF::RecieveCenterPointData(char* ValiableList, int nVarialbes, double *data,int *index,double SimulationTime,fTestPointInsideDomain TestPointInsideDomain) {
  int thread,i,j,k,idim,offset,cnt=0;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;


  //if the 'simulation time' is different from 'CouplingTimeCurrent' then exchange 'current' with 'last' datasets 
  if (CouplingTime!=SimulationTime) {
    //update the time counters 
    CouplingTime_last=CouplingTime;
    CouplingTime=SimulationTime;

    //flip the offsets
    int t;

    t=MagneticFieldOffset_last;
    MagneticFieldOffset_last=MagneticFieldOffset;
    MagneticFieldOffset=t;

    t=BulkVelocityOffset_last;
    BulkVelocityOffset_last=BulkVelocityOffset;
    BulkVelocityOffset=t;

    t=PlasmaPressureOffset_last;
    PlasmaPressureOffset_last=PlasmaPressureOffset;
    PlasmaPressureOffset=t;

    t=PlasmaNumberDensityOffset_last;
    PlasmaNumberDensityOffset_last=PlasmaNumberDensityOffset;
    PlasmaNumberDensityOffset=t;

    t=PlasmaTemperatureOffset_last;
    PlasmaTemperatureOffset_last=PlasmaTemperatureOffset;
    PlasmaTemperatureOffset=t;

    t=AlfvenWaveI01Offset_last;
    AlfvenWaveI01Offset_last=AlfvenWaveI01Offset;
    AlfvenWaveI01Offset=t;
  }

  //set up the 'first coupling occuerd' flag
  FirstCouplingOccured=true;

  //init the cell processing flags
  ResetCenterPointProcessingFlag();

  //determine the relation between SWMF's AMPS' variables
  int Rho_SWMF2AMPS=-1,Vx_SWMF2AMPS=-1,Bx_SWMF2AMPS=-1,P_SWMF2AMPS=-1,I01_SWMF2AMPS=-1;
  int i0=0,i1=0,n=0;
  char vname[200];


  while ((ValiableList[i0]!=0)&&(ValiableList[i0]==' ')) i0++;

  while (n<nVarialbes) {
    i1=i0;
    while (ValiableList[i1]!=' ') {
      vname[i1-i0]=tolower(ValiableList[i1]);
      i1++;
    }

    vname[i1-i0]=0;

    if ((strcmp(vname,"rho")==0)||(strcmp(vname,"swhrho")==0)) Rho_SWMF2AMPS=n;
    if ((strcmp(vname,"mx")==0)||(strcmp(vname,"swhmx")==0))  Vx_SWMF2AMPS=n;
    if ((strcmp(vname,"bx")==0)||(strcmp(vname,"swhbx")==0))  Bx_SWMF2AMPS=n;
    if ((strcmp(vname,"p")==0)||(strcmp(vname,"swhp")==0))   P_SWMF2AMPS=n;
    if (strcmp(vname,"i01")==0)   I01_SWMF2AMPS=n;

    n++;
    i0=i1;

    while ((ValiableList[i0]!=0)&&(ValiableList[i0]==' ')) i0++;
  }

  if ((Rho_SWMF2AMPS==-1)||(Vx_SWMF2AMPS==-1)||(Bx_SWMF2AMPS==-1)||(P_SWMF2AMPS==-1)) exit(__LINE__,__FILE__,"Error: background plasma macroscopic parameter is not found"); 

  //get coordinated of the center points
  for (node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL) {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            bool inside_domain_flag;

            inside_domain_flag=(TestPointInsideDomain!=NULL) ? TestPointInsideDomain(cell->GetX()) : true;

            if (inside_domain_flag==false) continue; 


            offset=nVarialbes*(index[cnt++]-1);
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;

            //convert momentum into velocity
            if ((Vx_SWMF2AMPS!=-1)&&(offset>=0)) {
              if (Rho_SWMF2AMPS!=-1) {
                for (idim=0;idim<3;idim++) data[offset+Vx_SWMF2AMPS+idim]=(data[offset+Rho_SWMF2AMPS]>0.0) ? data[offset+Vx_SWMF2AMPS+idim]/data[offset+Rho_SWMF2AMPS] : 0.0;
              }
              else for (idim=0;idim<3;idim++) data[offset+Vx_SWMF2AMPS+idim]=0.0;
            }

            //the order of the state vector: number density, temperature
            if ((offset>=0)&&(Rho_SWMF2AMPS>=0)) {
              *((double*)(cell->GetAssociatedDataBufferPointer()+PlasmaNumberDensityOffset))=data[offset+Rho_SWMF2AMPS]/MeanPlasmaAtomicMass;
              *((double*)(cell->GetAssociatedDataBufferPointer()+PlasmaTemperatureOffset))=(P_SWMF2AMPS>=0) ? data[offset+P_SWMF2AMPS]/(Kbol*data[offset+Rho_SWMF2AMPS]/MeanPlasmaAtomicMass) : 0.0;
            }
            else {
              *((double*)(cell->GetAssociatedDataBufferPointer()+PlasmaNumberDensityOffset))=0.0;
              *((double*)(cell->GetAssociatedDataBufferPointer()+PlasmaTemperatureOffset))=0.0;
            }

            //get pressure
            *((double*)(cell->GetAssociatedDataBufferPointer()+PlasmaPressureOffset))=((offset>=0)&&(P_SWMF2AMPS>=0)) ? data[offset+P_SWMF2AMPS] : 0.0;

            //AlfvenWave
            if ((offset>=0)&&(I01_SWMF2AMPS>=0)) { 
              ((double*)(cell->GetAssociatedDataBufferPointer()+AlfvenWaveI01Offset))[0]=data[offset+I01_SWMF2AMPS+0]; 
              ((double*)(cell->GetAssociatedDataBufferPointer()+AlfvenWaveI01Offset))[1]=data[offset+I01_SWMF2AMPS+1];
            }


            //bulk velocity and magnetic field
            for (idim=0;idim<3;idim++) {
              *((double*)(cell->GetAssociatedDataBufferPointer()+BulkVelocityOffset+idim*sizeof(double)))=((offset>=0)&&(Vx_SWMF2AMPS>=0)) ? data[offset+Vx_SWMF2AMPS+idim] : 0.0;
              *((double*)(cell->GetAssociatedDataBufferPointer()+MagneticFieldOffset+idim*sizeof(double)))=((offset>=0)&&(Bx_SWMF2AMPS>=0)) ? data[offset+Bx_SWMF2AMPS+idim] : 0.0;
            }

          
            //in case there are mode then just one fluid -> process other fluids
            for (int ifluid=1;ifluid<nCommunicatedIonFluids;ifluid++) {
              *(ifluid+(double*)(cell->GetAssociatedDataBufferPointer()+PlasmaNumberDensityOffset))=data[offset+8+5*(ifluid-1)]/MeanPlasmaAtomicMass;

              for (idim=0;idim<3;idim++) {
                *(idim+3*ifluid+(double*)(cell->GetAssociatedDataBufferPointer()+BulkVelocityOffset))=(data[offset+8+5*(ifluid-1)]>0.0) ? data[offset+9+idim+5*(ifluid-1)]/data[offset+8+5*(ifluid-1)] : 0.0; 
              }

              *(ifluid+(double*)(cell->GetAssociatedDataBufferPointer()+PlasmaPressureOffset))=data[offset+12+5*(ifluid-1)];
            }
          }
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh->nTotalThreads;thread++) for (node=PIC::Mesh::mesh->DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) if ((block=node->block)!=NULL) {
    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            bool inside_domain_flag;

            inside_domain_flag=(TestPointInsideDomain!=NULL) ? TestPointInsideDomain(cell->GetX()) : true;

            if (inside_domain_flag==false) continue;


            offset=nVarialbes*(index[cnt++]-1);
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;

            //convert momentum into velocity
            if ((Vx_SWMF2AMPS!=-1)&&(offset>=0)) {
              if (Rho_SWMF2AMPS!=-1) {
                for (idim=0;idim<3;idim++) (data[offset+Rho_SWMF2AMPS]>0.0) ? data[offset+Vx_SWMF2AMPS+idim]/=data[offset+Rho_SWMF2AMPS] : 0.0;
              }
              else for (idim=0;idim<3;idim++) data[offset+Vx_SWMF2AMPS+idim]=0.0;
            }

            //the order of the state vector: rho, V, B, p
            *((double*)(cell->GetAssociatedDataBufferPointer()+PlasmaNumberDensityOffset))=((offset>=0)&&(Rho_SWMF2AMPS>=0)) ? data[offset+Rho_SWMF2AMPS] : 0.0;
            *((double*)(cell->GetAssociatedDataBufferPointer()+PlasmaPressureOffset))=((offset>=0)&&(P_SWMF2AMPS>=0)) ? data[offset+P_SWMF2AMPS] : 0.0;

            for (idim=0;idim<3;idim++) {
              *((double*)(cell->GetAssociatedDataBufferPointer()+BulkVelocityOffset+idim*sizeof(double)))=((offset>=0)&&(Vx_SWMF2AMPS>=0)) ? data[offset+Vx_SWMF2AMPS+idim] : 0.0;
              *((double*)(cell->GetAssociatedDataBufferPointer()+MagneticFieldOffset+idim*sizeof(double)))=((offset>=0)&&(Bx_SWMF2AMPS>=0)) ? data[offset+Bx_SWMF2AMPS+idim] : 0.0;
            }

          }
        }
    }
  }
}



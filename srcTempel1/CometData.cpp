//$Id$

/* 
 * CometData.cpp
 *          
 *  Created on: Sept 21, 2015 
 *      Author: fougere 
 */

#include <stdio.h>
#include <stdlib.h>

#include "pic.h"

int Comet::CometData::NeutralsFromBinaryOffset=-1;
int Comet::CometData::nNeutrals=0;
int Comet::CometData::nMaxLoadedSpecies=0;

void Comet::CometData::WriteBinaryOutput(const char *fNameBase,int s,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,FILE *fout) {
  static CMPI_channel pipe;

  if (startNode==PIC::Mesh::mesh->rootTree) {
    pipe.init(1000000);

    if (PIC::Mesh::mesh->ThisThread==0) {
      pipe.openRecvAll();
    }
    else pipe.openSend(0); 
  }
  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;
  int idim;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  char SendCellFlag;

  const int nTotalVariables=4;
  
  double data[nTotalVariables];
  double GasBulkVelocity[3];

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((PIC::ThisThread==0)||(startNode->Thread==PIC::ThisThread)) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      SendCellFlag=true;

      //determine whether the cell data neede to be saved
      if (startNode->Thread==PIC::ThisThread) {
        //locate the cell
        if (startNode->block==NULL) {
	  SendCellFlag=false; 
	  for (idim=0;idim<4;idim++) data[idim]=0.0;
	}
        nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k);
        if (SendCellFlag==true) {
	  if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) {
	    data[0]=startNode->block->GetCenterNode(nd)->GetNumberDensity(s)*PIC::MolecularData::GetMass(s);

	    startNode->block->GetCenterNode(nd)->GetBulkVelocity(GasBulkVelocity,s);
	    
	    data[1]=GasBulkVelocity[0];
	    data[2]=GasBulkVelocity[1];
	    data[3]=GasBulkVelocity[2];
	  }
	  else {
	    SendCellFlag=false;
	    for (idim=0;idim<4;idim++) data[idim]=0.0;
	  }
        }

        if (startNode->Thread!=0) pipe.send(SendCellFlag);
      }
      else {
        pipe.recv(SendCellFlag,startNode->Thread);
      }

      //save the cell data saving flag 
      if (PIC::ThisThread==0) fwrite(&SendCellFlag,sizeof(char),1,fout);

      //save the cell data
      if (SendCellFlag==true) {
	if (startNode->Thread==PIC::ThisThread) {
	  if (startNode->Thread==0) {
	    fwrite(data,nTotalVariables*sizeof(double),1,fout);
	  }
	  else {
            pipe.send((double*)(data),nTotalVariables);
          }
	}
	else {
	  pipe.recv(data,nTotalVariables,startNode->Thread);
          fwrite(data,nTotalVariables*sizeof(double),1,fout);
	}
      }
	  }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) Comet::CometData::WriteBinaryOutput(NULL,s,startNode->downNode[nDownNode],fout);
  }
  

  if (startNode==PIC::Mesh::mesh->rootTree) {
    if (PIC::Mesh::mesh->ThisThread==0) {
      pipe.closeRecvAll();
    }
    else pipe.closeSend();

    pipe.remove();
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }

}

void Comet::CometData::LoadBinaryFile(const char *fNameBase) {
  char fname[400];
  FILE *fData;

  sprintf(fname,"%s/amr.sig=0x%lx.f=%s.CenterNodeOutputData.bin",PIC::CPLR::DATAFILE::path,PIC::Mesh::mesh->getMeshSignature(),fNameBase);
  fData=fopen(fname,"r");

  if (fData==NULL) {
    char msg[400];

    sprintf(msg,"Error: The neutral binary file (%s) does not exist.",fname);
    exit(__LINE__,__FILE__,msg);
  }

  //read the number of the saved species
  fread(&nNeutrals,sizeof(int),1,fData);

  if (nNeutrals>nMaxLoadedSpecies) exit(__LINE__,__FILE__,"Error: nNeutrals>nMaxLoadedSpecies. Solution: increase Comet::CometData::nMaxLoadedSpecies to accomodate all speces that are loaded from a binary file");

  //load the background data for each species
  for (int iSpecies=0;iSpecies<nNeutrals;iSpecies++) LoadBinaryFile_Internal(iSpecies,fData,PIC::Mesh::mesh->rootTree);

  //close the data file
  fclose(fData);
}


void Comet::CometData::LoadBinaryFile_Internal(int iSpecies,FILE* fData,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static int CenterNodeAssociatedLength;
  
  //loop through all points                                                                                                                
  //create the list of the points                                                                                                          
  //perform the interpolation loop                                                                                                         
  int i,j,k,nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  
  if (startNode==PIC::Mesh::mesh->rootTree) {
    PIC::Mesh::cDataCenterNode cell;
    CenterNodeAssociatedLength=cell.AssociatedDataLength();
  }

  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  char savedLoadCellFlag;

  const int nTotalVariables=4;

  int LoadVariableOffset[nTotalVariables]={NeutralsFromBinaryOffset,NeutralsFromBinaryOffset+sizeof(double),NeutralsFromBinaryOffset+sizeof(double)*2,NeutralsFromBinaryOffset+sizeof(double)*3};

  double data[nTotalVariables];
  int CenterNodeAssociatedDataOffsetBegin=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->block==NULL) {
      //the block belongs to a other processor -> skip all data                                                                           
      for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)  {
	    fread(&savedLoadCellFlag,sizeof(char),1,fData);

	    if (savedLoadCellFlag==true) {
	      //the cell data is saved -> skip it                                                                                              
	      fseek(fData,nTotalVariables*sizeof(double),SEEK_CUR);
	    }
	  }

    }
    else for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
	    fread(&savedLoadCellFlag,sizeof(char),1,fData);
	    
	    if (savedLoadCellFlag==true) {
	      //determine whether the cell data needed to be read                                                                                
	      //locate the cell                                                                                                                  
	      nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k);
	      
	      if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) {
		//read the data                                                                                                                  
		offset=CenterNode->GetAssociatedDataBufferPointer();

		//read center cells' associated data                                                                                           
		fread(data,nTotalVariables*sizeof(double),1,fData);

		//copy the data                       
		memcpy(4*iSpecies*sizeof(double)+offset+NeutralsFromBinaryOffset,data,4*sizeof(double));

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
		PIC::Debugger::CatchOutLimitValue((double*)(offset+CenterNodeAssociatedDataOffsetBegin),nTotalVariables,__LINE__,__FILE__);
#endif
	      }
	      else {
		//the cell data is saved -> skip it                                                                                              
		fseek(fData,nTotalVariables*sizeof(double),SEEK_CUR);
	      }
	    }

	  }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) Comet::CometData::LoadBinaryFile_Internal(iSpecies,fData,startNode->downNode[nDownNode]);
  }
}

int Comet::CometData::RequestDataBuffer(int offset) {
  int TotalDataLength;

  NeutralsFromBinaryOffset=offset;
  TotalDataLength=4*nMaxLoadedSpecies;

  return TotalDataLength*sizeof(double);
}

void Comet::CometData::PrintVariableList(FILE* fout,int DataSetNumber) {
  for (int iSpecies=0;iSpecies<nNeutrals;iSpecies++) fprintf(fout,",\"MassDensityNeutral_%i\",\"Vneutralx_%i\",\"Vneutraly_%i\",\"Vneutralz_%i\"",iSpecies,iSpecies,iSpecies,iSpecies);
}

void Comet::CometData::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;
  double t;

  for (idim=0;idim<4*nNeutrals;idim++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+ NeutralsFromBinaryOffset+idim*sizeof(double)));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
}

void Comet::CometData::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double Neutrals[4*nNeutrals];
  int i,idim;
  char *SamplingBuffer;

  for (i=0;i<4*nNeutrals;i++) Neutrals[i]=0.0;

  for (i=0;i<nInterpolationCoeficients;i++) {

    for (idim=0,SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+NeutralsFromBinaryOffset;idim<4*nNeutrals;idim++) Neutrals[idim]+=(*((double*)(SamplingBuffer+idim*sizeof(double))))*InterpolationCoeficients[i];
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+NeutralsFromBinaryOffset,Neutrals,4*nNeutrals*sizeof(double));
}

double Comet::CometData::GetNeutralsMassDensity(int s,PIC::Mesh::cDataCenterNode *cell) {
  double res=0.0;

#if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
  if ((s<0)||(s>=nNeutrals)) exit(__LINE__,__FILE__,"Error: 's' is out of the range");
#endif

  res=*((double*)(cell->GetAssociatedDataBufferPointer()+NeutralsFromBinaryOffset+4*s*sizeof(double)));
  return res;
}

double Comet::CometData::GetNeutralsMassDensity(int s,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  return GetNeutralsMassDensity(s,node->block->GetCenterNode(nd));
}


void  Comet::CometData::GetNeutralsVelocity(double *x, int s,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  register int idim;
  register double *offset=(double*)(4*s*sizeof(double)+sizeof(double)+NeutralsFromBinaryOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());
  
  for (idim=0;idim<3;idim++) x[idim]=offset[idim];
}

int Comet::CometData::GetnNeutrals() {
  return nNeutrals;
}

void Comet::CometData::SetnNeutrals(int n) {
  nNeutrals=n;
}

//determine the ckecksum of the background density data
void Comet::CometData::PrintCheckSum(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  static CRC32 checksum;
  static CMPI_channel pipe(1000000);
  static bool FirstCallFlag=true;

  PIC::Mesh::cDataCenterNode *cell;

  if (FirstCallFlag==true) {
    checksum.clear();

    node=PIC::Mesh::mesh->rootTree;
    FirstCallFlag=false;

    if (PIC::Mesh::mesh->ThisThread==0) pipe.openRecvAll();
    else pipe.openSend(0);
  }

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    int s,nd;
    double t;
    PIC::Mesh::cDataCenterNode *cell;

    for (s=0;s<nNeutrals;s++) for (nd=0;nd<_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_;nd++) {
      if (PIC::ThisThread==node->Thread) {
        cell=node->block->GetCenterNode(nd);
        t=(cell!=NULL) ? GetNeutralsMassDensity(s,cell) : -1.0;

        if (PIC::ThisThread!=0) pipe.send(t);
      }
      else {
        if (PIC::ThisThread==0) pipe.recv(t,node->Thread);
      }

      if (PIC::ThisThread==0) checksum.add(t);
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (node->downNode[nDownNode]!=NULL) PrintCheckSum(node->downNode[nDownNode]);
  }

  if (node==PIC::Mesh::mesh->rootTree) {
    if (PIC::ThisThread==0) pipe.closeRecvAll();
    else pipe.closeSend();

    checksum.PrintChecksum(__LINE__,__FILE__);
    FirstCallFlag=true;
  }
}

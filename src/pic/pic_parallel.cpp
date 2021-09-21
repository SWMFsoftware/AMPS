//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the functions that control the interprocessor communication of the code

#include <map>
#include "pic.h"

long int PIC::Parallel::sendParticleCounter=0,PIC::Parallel::recvParticleCounter=0,PIC::Parallel::IterationNumberAfterRebalancing=0;
double PIC::Parallel::RebalancingTime=0.0,PIC::Parallel::CumulativeLatency=0.0;
double PIC::Parallel::EmergencyLoadRebalancingFactor=3.0;
double PIC::Parallel::Latency=0.0;

//processing 'corner' and 'center' node associated data vectors when perform syncronization
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData=NULL;
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData=NULL;
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::CopyCenterNodeAssociatedData=NULL;

PIC::Parallel::CenterBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData=NULL;
PIC::Parallel::CenterBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData=NULL;
 
bool PIC::Parallel::CornerBlockBoundaryNodes::ActiveFlag=false;
void PIC::Parallel::CornerBlockBoundaryNodes::SetActiveFlag(bool flag) {ActiveFlag=flag;}

bool PIC::Parallel::CenterBlockBoundaryNodes::ActiveFlag=false;
void PIC::Parallel::CenterBlockBoundaryNodes::SetActiveFlag(bool flag) {ActiveFlag=flag;}

PIC::Parallel::BoundaryProcessManager PIC::Parallel::BPManager;

//default function forprocessing of the corner node associated data
void PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength);
}

void PIC::Parallel::CornerBlockBoundaryNodes::CopyCenterNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength);
}

//-----------------------------------------------------------------


/*
//processing 'corner' and 'center' node associated data vectors when perform syncronization
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::ProcessCenterNodeAssociatedData=NULL,PIC::Parallel::ProcessCornerNodeAssociatedData=NULL;
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CopyCenterNodeAssociatedData=PIC::Parallel::CopyCenterNodeAssociatedData_default;
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CopyCornerNodeAssociatedData=PIC::Parallel::CopyCornerNodeAssociatedData_default;

//default function forprocessing of the corner node associated data
void PIC::Parallel::CopyCornerNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength);
}

void PIC::Parallel::CopyCenterNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength);
}
*/ 

//====================================================
//Exchange particles between Processors
void PIC::Parallel::ExchangeParticleData() {
  int From,To,i,iFrom,flag;
  long int Particle,NextParticle,newParticle,LocalCellNumber=-1;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *sendNode=NULL,*recvNode=NULL;

  int RecvMessageSizeRequestTableLength=0;
  int SendMessageSizeRequestTableLength=0;

  MPI_Request RecvPaticleDataRequestTable[PIC::nTotalThreads];
  int RecvPaticleDataRequestTableLength=0;
  int RecvPaticleDataProcessTable[PIC::nTotalThreads];

  int SendParticleDataRequestTableLength=0;

  MPI_Request SendParticleDataRequestTable[PIC::nTotalThreads];

  MPI_Request SendMessageSizeRequestTable[PIC::nTotalThreads];
  MPI_Request RecvMessageSizeRequestTable[PIC::nTotalThreads];

  int SendMessageLengthTable[PIC::nTotalThreads];
  int RecvMessageLengthTable[PIC::nTotalThreads];

  int SendMessageLengthProcessTable[PIC::nTotalThreads];
  int RecvMessageLengthProcessTable[PIC::nTotalThreads];

  for (int thread=0;thread<PIC::nTotalThreads;thread++) {
    SendMessageLengthTable[thread]=0,RecvMessageLengthTable[thread]=0;
  }

  //set the default value inthe counters
  for (i=0;i<PIC::nTotalThreads;i++) SendMessageLengthTable[i]=0,RecvMessageLengthTable[i]=0;

  #if DIM == 3
  //  cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
  cAMRnodeID nodeid;
  #elif DIM == 2
  cMeshAMR2d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
  #else
  cMeshAMR1d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
  #endif


  //The signals
  int Signal;
  const int _NEW_BLOCK_ID_SIGNAL_=       0;
  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  struct cMessageDescriptorList {
    int iCell,jCell,kCell;
    int nTotalParticles;
    cAMRnodeID node_id;
  };

  

  auto PrepareMessageDescriptorTable = [&] (list<cMessageDescriptorList> *MessageDescriptorList,list<MPI_Aint> *SendParticleList,int& nTotalSendParticles,int To) {
    int npart,SendSellNumber=0;
    long int *FirstCellParticleTable;
    bool block_header_saved,cell_header_saved; 

    cMessageDescriptorList p;
    MPI_Aint particle_data;

    nTotalSendParticles=0;

    if ((PIC::ThisThread!=To)&&(PIC::Mesh::mesh->ParallelSendRecvMap[PIC::ThisThread][To]==true)) {
      for (sendNode=PIC::Mesh::mesh->DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) if (sendNode->block!=NULL) {
        FirstCellParticleTable=sendNode->block->FirstCellParticleTable;

        for (int kCell=0;kCell<kCellMax;kCell++) for (int jCell=0;jCell<jCellMax;jCell++) for (int iCell=0;iCell<iCellMax;iCell++) {
          Particle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

          if  (Particle!=-1) {
            p.nTotalParticles=0;
            p.iCell=iCell,p.jCell=jCell,p.kCell=kCell;
            p.node_id=sendNode->AMRnodeID;  
          
            SendSellNumber++;

            while (Particle!=-1) {
              MPI_Get_address(PIC::ParticleBuffer::ParticleDataBuffer+PIC::ParticleBuffer::GetParticleDataOffset(Particle),&particle_data); 
              SendParticleList->push_back(particle_data);

              p.nTotalParticles++;
              nTotalSendParticles++;

              Particle=PIC::ParticleBuffer::GetNext(Particle);
            }

            MessageDescriptorList->push_back(p);
          }
        }
      }
    }

    return SendSellNumber;
  };

 

  class cSendDataInfo {
  public:
    int nTotalSendCells,nTotalSendParticles;

    cSendDataInfo() {
      nTotalSendCells=0,nTotalSendParticles=0;
    }
  };


  auto InitSendParticleData = [&] (int To,int nTotalSendParticles,list<MPI_Aint> *SendParticleList,MPI_Request *SendParticleDataRequest) { 
    MPI_Aint *offset_table=new MPI_Aint [nTotalSendParticles];
    int iptr;
    list<MPI_Aint>::iterator it;

    for (iptr=0,it=SendParticleList->begin();it!=SendParticleList->end();it++,iptr++) {
      offset_table[iptr]=*it; 
    }
       
    MPI_Datatype particle_send_type;

    MPI_Type_create_hindexed_block(nTotalSendParticles,PIC::ParticleBuffer::ParticleDataLength-2*sizeof(long int),offset_table,MPI_BYTE,&particle_send_type);
    MPI_Type_commit(&particle_send_type);

    MPI_Isend(MPI_BOTTOM,1,particle_send_type,To,12,MPI_GLOBAL_COMMUNICATOR,SendParticleDataRequest);

    MPI_Type_free(&particle_send_type);
    delete [] offset_table;
  };

  auto RemoveSentParticles = [&] (cMessageDescriptorList *SendMessageDescriptor,int nSendCells) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL;
    int ptr,iCell,jCell,kCell,next,cnt;
    cAMRnodeID node_id;
    cMessageDescriptorList *it;

    for (cnt=0;cnt<nSendCells;cnt++) {
      it=SendMessageDescriptor+cnt; 

      if ((node==NULL)||(node_id!=it->node_id)) {
        node_id=it->node_id;
        node=PIC::Mesh::mesh->findAMRnodeWithID(node_id);
      }

      ptr=node->block->FirstCellParticleTable[it->iCell+_BLOCK_CELLS_X_*(it->jCell+_BLOCK_CELLS_Y_*it->kCell)]; 

      while (ptr!=-1) {
        next=PIC::ParticleBuffer::GetNext(ptr);

        PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr,true);
        ptr=next;
      }

      node->block->FirstCellParticleTable[it->iCell+_BLOCK_CELLS_X_*(it->jCell+_BLOCK_CELLS_Y_*it->kCell)]=-1;
    }
  };

  auto InitNewParticles = [&] (int nTotalSendCells,cMessageDescriptorList *MessageDescriptorTable,list<MPI_Aint> *RecvParticleList) {
    int icell,i,j,k;
    int ptr_cnt,new_particle;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=NULL;
    cAMRnodeID node_id;
    cMessageDescriptorList *p;
    int particle_data_buffer_offset=0;

    MPI_Aint particle_data;

    for (icell=0;icell<nTotalSendCells;icell++) {
      if ((node==NULL)||(node_id!=MessageDescriptorTable[icell].node_id)) { 
        node_id=MessageDescriptorTable[icell].node_id;
        node=PIC::Mesh::mesh->findAMRnodeWithID(node_id);
      }

      p=MessageDescriptorTable+icell;      

      for (ptr_cnt=0;ptr_cnt<p->nTotalParticles;ptr_cnt++) {  
        new_particle=PIC::ParticleBuffer::GetNewParticle(node->block->FirstCellParticleTable[p->iCell+_BLOCK_CELLS_X_*(p->jCell+_BLOCK_CELLS_Y_*p->kCell)],true);  

        MPI_Get_address(PIC::ParticleBuffer::ParticleDataBuffer+PIC::ParticleBuffer::GetParticleDataOffset(new_particle),&particle_data);
        RecvParticleList->push_back(particle_data);
      }
   }
 }; 


  cMessageDescriptorList *SendMessageDescriptorTable[PIC::nTotalThreads];
  cSendDataInfo SendDataInfo[PIC::nTotalThreads];

  MPI_Request SendRequestTable[3*PIC::nTotalThreads];
  int SendRequestTableLength=0;

  for (int thread=0;thread<PIC::Mesh::mesh->nTotalThreads;thread++) SendMessageDescriptorTable[thread]=NULL; 

  //send send_descriptor and particle data
  for (To=0;To<PIC::Mesh::mesh->nTotalThreads;To++) if ((PIC::ThisThread!=To)&&(PIC::Mesh::mesh->ParallelSendRecvMap[PIC::ThisThread][To]==true)) {
    int i,n_sent_cells,nTotalSendParticles;
    list<MPI_Aint> SendParticleList;
    list<cMessageDescriptorList> MessageDescriptorList;

    n_sent_cells=PrepareMessageDescriptorTable(&MessageDescriptorList,&SendParticleList,nTotalSendParticles,To); 

    SendDataInfo[To].nTotalSendCells=n_sent_cells;
    SendDataInfo[To].nTotalSendParticles=nTotalSendParticles;

    MPI_Isend(SendDataInfo+To,sizeof(cSendDataInfo),MPI_BYTE,To,9,MPI_GLOBAL_COMMUNICATOR,SendRequestTable+SendRequestTableLength);
    SendRequestTableLength++;

    //send particle data
    if (n_sent_cells!=0) {
      list<cMessageDescriptorList>::iterator it;

      SendMessageDescriptorTable[To]=new cMessageDescriptorList [n_sent_cells]; 

      for (i=0,it=MessageDescriptorList.begin();it!=MessageDescriptorList.end();it++,i++) SendMessageDescriptorTable[To][i]=*it;
      
      MPI_Isend(SendMessageDescriptorTable[To],n_sent_cells*sizeof(cMessageDescriptorList),MPI_BYTE,To,10,MPI_GLOBAL_COMMUNICATOR,SendRequestTable+SendRequestTableLength);
      SendRequestTableLength++;

      //send the particle data
      InitSendParticleData(To,nTotalSendParticles,&SendParticleList,SendRequestTable+SendRequestTableLength);
      SendRequestTableLength++;
    }
  }


  //recv send_secriptor and particle data
  MPI_Request RecvRequestDescriptorTable[PIC::nTotalThreads];
  int RecvRequestDescriptorTableLength=0;

  cSendDataInfo RecvDataInfo[PIC::nTotalThreads];
  int RecvDataInfoMap[PIC::nTotalThreads];

  MPI_Request RecvRequestDataInfoTable[PIC::nTotalThreads];
  int RecvRequestDataInfoTableLength=0;

  for (int thread=0;thread<PIC::nTotalThreads;thread++) RecvDataInfoMap[thread]=-1; //,RecvCellNumber[thread]=0;

  for (From=0;From<PIC::nTotalThreads;From++) if ((PIC::ThisThread!=From)&&(PIC::Mesh::mesh->ParallelSendRecvMap[From][PIC::ThisThread]==true)) {
    MPI_Irecv(RecvDataInfo+From,sizeof(cSendDataInfo),MPI_BYTE,From,9,MPI_GLOBAL_COMMUNICATOR,RecvRequestDataInfoTable+RecvRequestDataInfoTableLength);
    RecvDataInfoMap[RecvRequestDataInfoTableLength]=From;

    RecvRequestDataInfoTableLength++;
  }

  //recive the particle data
  MPI_Request RecvRequestMessageDescriptorTable[PIC::nTotalThreads];
  int RecvRequestMessageDescriptorTableLength=0;

  MPI_Request RecvRequestParticleDataTable[PIC::nTotalThreads];
  int RecvRequestParticleDataTableLength=0;


  int RecvRequestTableLength=0;
  int nRecvPackages=0;

  bool RecvMessageDescriptorFlagTable[PIC::nTotalThreads];
  bool RecvParticleDataFlagTable[PIC::nTotalThreads]; 

  cMessageDescriptorList *RecvMessageDescriptorTable[PIC::nTotalThreads];
  PIC::ParticleBuffer::byte *RecvParticleDataBuffer[PIC::nTotalThreads]; 
  
  int RecvParticleDataMap[PIC::nTotalThreads];

  for (int thread=0;thread<PIC::nTotalThreads;thread++) {
    RecvParticleDataMap[thread]=-1;
    RecvMessageDescriptorTable[thread]=NULL,RecvParticleDataBuffer[thread]=NULL;
    RecvMessageDescriptorFlagTable[thread]=false,RecvParticleDataFlagTable[thread]=false;
  }

  int MessageDescriptorWaiting=0;

  while ((RecvRequestDataInfoTableLength!=0)||(MessageDescriptorWaiting!=0)) {
    if (RecvRequestDataInfoTableLength!=0) {
      MPI_Testany(RecvRequestDataInfoTableLength,RecvRequestDataInfoTable,&iFrom,&flag,MPI_STATUS_IGNORE);

      if ((flag==true)&&(iFrom!=MPI_UNDEFINED)) {
        int From=RecvDataInfoMap[iFrom];

        MPI_Wait(RecvRequestDataInfoTable+iFrom,MPI_STATUS_IGNORE); 
        RecvRequestDataInfoTable[iFrom]=RecvRequestDataInfoTable[RecvRequestDataInfoTableLength-1];
        RecvDataInfoMap[iFrom]=RecvDataInfoMap[RecvRequestDataInfoTableLength-1];

        RecvRequestDataInfoTableLength--;

        if (RecvDataInfo[From].nTotalSendCells!=0) {
          //init recieving the message descriptor
          RecvMessageDescriptorTable[From]=new cMessageDescriptorList[RecvDataInfo[From].nTotalSendCells];
          RecvParticleDataBuffer[From]=new PIC::ParticleBuffer::byte[PIC::ParticleBuffer::ParticleDataLength*RecvDataInfo[From].nTotalSendParticles];

          RecvParticleDataMap[RecvRequestMessageDescriptorTableLength]=From;

          MPI_Irecv(RecvMessageDescriptorTable[From],RecvDataInfo[From].nTotalSendCells*sizeof(cMessageDescriptorList),MPI_BYTE,From,10,MPI_GLOBAL_COMMUNICATOR,RecvRequestMessageDescriptorTable+RecvRequestMessageDescriptorTableLength);
          RecvRequestMessageDescriptorTableLength++;
          MessageDescriptorWaiting++;
        } 
      }
    }


    if (RecvRequestMessageDescriptorTableLength!=0) {
      MPI_Testany(RecvRequestMessageDescriptorTableLength,RecvRequestMessageDescriptorTable,&iFrom,&flag,MPI_STATUS_IGNORE);

      if ((flag==true)&&(iFrom!=MPI_UNDEFINED)) {
        int From=RecvParticleDataMap[iFrom];
        list<MPI_Aint> RecvParticleList;

        MessageDescriptorWaiting--;

        InitNewParticles(RecvDataInfo[From].nTotalSendCells,RecvMessageDescriptorTable[From],&RecvParticleList);  

        MPI_Aint *offset_table=new MPI_Aint [RecvDataInfo[From].nTotalSendParticles];
        int iptr;
        list<MPI_Aint>::iterator it;

        for (iptr=0,it=RecvParticleList.begin();it!=RecvParticleList.end();it++,iptr++) {
          offset_table[iptr]=*it;
        }

        MPI_Datatype particle_send_type;

        MPI_Type_create_hindexed_block(RecvDataInfo[From].nTotalSendParticles,PIC::ParticleBuffer::ParticleDataLength-2*sizeof(long int),offset_table,MPI_BYTE,&particle_send_type);
        MPI_Type_commit(&particle_send_type);

        MPI_Irecv(MPI_BOTTOM,1,particle_send_type,From,12,MPI_GLOBAL_COMMUNICATOR,RecvRequestParticleDataTable+RecvRequestParticleDataTableLength);
        RecvRequestParticleDataTableLength++;  

        MPI_Type_free(&particle_send_type);
        delete [] offset_table;
      } 
    }
  }
  

  //wait for completing recv operation (clean memory used by MPI) 
  MPI_Waitall(RecvRequestParticleDataTableLength,RecvRequestParticleDataTable,MPI_STATUSES_IGNORE);
  MPI_Waitall(RecvRequestMessageDescriptorTableLength,RecvRequestMessageDescriptorTable,MPI_STATUSES_IGNORE); 

  //wait for completing send operations 
  MPI_Waitall(SendRequestTableLength,SendRequestTable,MPI_STATUSES_IGNORE);

  //delete particle that were send out 
  for (int thread=0;thread<PIC::Mesh::mesh->nTotalThreads;thread++) if (SendMessageDescriptorTable[thread]!=NULL) {
    RemoveSentParticles(SendMessageDescriptorTable[thread],SendDataInfo[thread].nTotalSendCells); 

    delete [] SendMessageDescriptorTable[thread];
  }


  //delete temp recv buffers
  for (int thread=0;thread<PIC::nTotalThreads;thread++) {
    if (RecvMessageDescriptorTable[thread]!=NULL) delete [] RecvMessageDescriptorTable[thread];
    if (RecvParticleDataBuffer[thread]!=NULL) delete [] RecvParticleDataBuffer[thread]; 
  }



  //in case AMPS is compiled in the debugger mode, source particles to eliminate the randomness of the particles in the particle lists
  //caused by the randomness of the particle data exchange time


  auto SortParticleCellList = [&] (int iCell,int jCell,int kCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node) {

    class cParticleDescriptor {
    public:
      long int ptr;
      unsigned long int checksum;

      bool operator < (const cParticleDescriptor& p) const {
        return checksum<p.checksum;
      }
    };

    //put together the particle list
    vector <cParticleDescriptor> ParticleList;
    cParticleDescriptor p;
    long int Particle;
    CRC32 sig;


    Particle=Node->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

    while (Particle!=-1) {
      p.ptr=Particle;
      p.checksum=PIC::ParticleBuffer::GetParticleSignature(Particle,&sig,false);
      sig.clear();

      ParticleList.push_back(p);
      Particle=PIC::ParticleBuffer::GetNext(Particle);
    }


    //sort the list
    std::sort(ParticleList.begin(),ParticleList.end(),
        [](const cParticleDescriptor& a,const cParticleDescriptor& b) {return a.checksum>b.checksum;});

    //create the sorted list
    int ip,ipmax=ParticleList.size();

    Node->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=ParticleList[0].ptr;

    for (ip=0;ip<ipmax;ip++) {
      PIC::ParticleBuffer::SetPrev((ip!=0) ? ParticleList[ip-1].ptr : -1,ParticleList[ip].ptr);
      PIC::ParticleBuffer::SetNext((ip!=ipmax-1) ? ParticleList[ip+1].ptr : -1,ParticleList[ip].ptr);
    }
  };

  auto SortParticleList = [&] () {
    int nLocalNode;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node;
    int i,j,k;


    for (nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      Node=DomainBlockDecomposition::BlockTable[nLocalNode];

      if (Node->block!=NULL) {
        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          if (Node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]!=-1) {
            SortParticleCellList(i,j,k,Node);
          }
        }
      }
    }
  };

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    SortParticleList();
  }

 }






//void PIC::Parallel::ProcessCornerBlockBoundaryNodes() {
//  int iThread,i,j,k,iface;
//  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
//
//  if (CornerBlockBoundaryNodes::ActiveFlag==false) return;
//
//
//  const int iFaceMin[6]={0,_BLOCK_CELLS_X_,0,0,0,0};
//  const int iFaceMax[6]={0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_};
//
//  const int jFaceMin[6]={0,0,0,_BLOCK_CELLS_Y_,0,0};
//  const int jFaceMax[6]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_};
//
//  const int kFaceMin[6]={0,0,0,0,0,_BLOCK_CELLS_Z_};
//  const int kFaceMax[6]={_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,0,_BLOCK_CELLS_Z_};
//
//  struct cStencilElement {
//    int StencilLength; //represent distinct thread number
//    int StencilThreadTable[80];
//    char *AssociatedDataPointer;
//    std::vector<char *> localDataPntArr; //save local data buffer pointers
//  };
//
//  static int StencilTableLength=0;
//  static cStencilElement *St,**StencilTable=NULL;
//
//  static cStack<cStencilElement> StencilElementStack;
//
//  //generate a new stencil table
//  static int nMeshModificationCounter=-1;
//
//  //determine whether the mesh/domain decomposition have been changed
//  int localMeshChangeFlag,globalMeshChangeFlag;
//
//  localMeshChangeFlag=(nMeshModificationCounter==PIC::Mesh::mesh->nMeshModificationCounter) ? 0 : 1;
//  MPI_Allreduce(&localMeshChangeFlag,&globalMeshChangeFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
//
//  static int * nPointsToComm=NULL;
//  static MPI_Request * SendPntNumberList=NULL;
//  static MPI_Request * RecvPntNumberList=NULL;
//  static double * SendCoordBuff = NULL;
//  static MPI_Request * SendCoordList=NULL;
//  static MPI_Request * RecvCoordList=NULL;
//  static MPI_Request * SendDataBufferPtrList=NULL;
//  static MPI_Request * RecvDataBufferPtrList=NULL;
//  static char **  SendDataBuffPointerBuff =NULL;
//  static double ** RecvCoordBuff = NULL;
//  static char *** RecvDataBuffPointerBuff=NULL;
//
//  if (globalMeshChangeFlag!=0) {
//    //the mesh or the domain decomposition has been modified. Need to create a new communucation table
//    bool meshModifiedFlag_CountMeshElements=PIC::Mesh::mesh->meshModifiedFlag_CountMeshElements;
//    double StartTime=MPI_Wtime();
//
//
//    if (nPointsToComm==NULL) nPointsToComm = new int [PIC::nTotalThreads];
//    if (SendPntNumberList==NULL) SendPntNumberList = new MPI_Request[PIC::nTotalThreads-1];
//    if (RecvPntNumberList==NULL) RecvPntNumberList = new MPI_Request[PIC::nTotalThreads-1];
//
//    if (!SendCoordList) SendCoordList = new MPI_Request[PIC::nTotalThreads-1];
//    if (!RecvCoordList) RecvCoordList = new MPI_Request[PIC::nTotalThreads-1];
//    if (!SendDataBufferPtrList) SendDataBufferPtrList = new MPI_Request[PIC::nTotalThreads-1];
//    if (!RecvDataBufferPtrList) RecvDataBufferPtrList = new MPI_Request[PIC::nTotalThreads-1];
//
//    int  nPointsToCommThisThread = 0;
//    //determine the new length of the table
//    //    for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
//    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
//      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
//      double eps = 0.3*PIC::Mesh::mesh->EPS;
//      double dx[3];
//      int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
//      for (int idim=0; idim<3; idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];
//
//      if (node->Thread==PIC::ThisThread && node->block) {
//        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
//
//          double x[3];
//          int ind[3]={i,j,k};
//          bool outBoundary=false;// atBoundary=false;
//          for (int idim=0; idim<3; idim++) {
//            x[idim]=node->xmin[idim]+dx[idim]*ind[idim];
//
//            #if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_
//            if ((x[idim]-PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim])<-eps || (x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim])>eps){
//              outBoundary=true;
//              break;
//            }
//            #endif
//          }
//
//          if (outBoundary==false) nPointsToCommThisThread++;
//
//              }// for (iface=0;iface<6;iface++)  for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)
//      }//if (node->Thread==PIC::ThisThread)
//    }//for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)
//
//    nPointsToComm[PIC::ThisThread] = nPointsToCommThisThread;
//
//    int iProcessCnt=0;
//    //send the number of points on this proc to all other proc
//    for (int iProc=0; iProc<PIC::nTotalThreads; iProc++){
//      if (iProc!=PIC::ThisThread){
//        MPI_Isend(nPointsToComm+PIC::ThisThread,1,MPI_INT,iProc,0,MPI_GLOBAL_COMMUNICATOR,SendPntNumberList+iProcessCnt);
//        iProcessCnt++;
//      }
//    }
//
//    iProcessCnt=0;
//    //recv the number of points on other proc
//    for (int iProc=0; iProc<PIC::nTotalThreads; iProc++){
//      if (iProc!=PIC::ThisThread){
//        MPI_Irecv(nPointsToComm+iProc,1,MPI_INT,iProc,0,MPI_GLOBAL_COMMUNICATOR,RecvPntNumberList+iProcessCnt);
//        iProcessCnt++;
//      }
//    }
//
//    //printf("thisthread:%d, nPointsToComm:%d\n",PIC::ThisThread,nPointsToComm[PIC::ThisThread]);
//
//    if (SendCoordBuff!=NULL) delete [] SendCoordBuff;
//    SendCoordBuff = new double [3*nPointsToComm[PIC::ThisThread]];
//    if (RecvCoordBuff == NULL) {
//      RecvCoordBuff = new double * [PIC::nTotalThreads-1];
//      for (i=0; i<PIC::nTotalThreads-1; i++) RecvCoordBuff[i] = NULL;
//    }
//    if (SendDataBuffPointerBuff!=NULL) delete [] SendDataBuffPointerBuff;
//    SendDataBuffPointerBuff =  new char * [nPointsToComm[PIC::ThisThread]];
//    if  (RecvDataBuffPointerBuff == NULL) {
//      RecvDataBuffPointerBuff= new char ** [PIC::nTotalThreads-1];
//      for (i=0; i<PIC::nTotalThreads-1; i++) RecvDataBuffPointerBuff[i] = NULL;
//    }
//
//    double * tempPnt = SendCoordBuff;
//    char ** tempDataBuffPnt = SendDataBuffPointerBuff;
//    //    for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
//    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
//      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
//
//      if (node->Thread==PIC::ThisThread && node->block) {
//        double dx[3];
//        int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
//        double eps = 0.3*PIC::Mesh::mesh->EPS;
//        for (int idim=0; idim<3; idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];
//
//        for (iface=0;iface<6;iface++)  for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
//
//          double x[3];
//          int ind[3]={i,j,k};
//          bool outBoundary=false;
//
//          // if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_){
//          for (int idim=0; idim<3; idim++) {
//            x[idim]=node->xmin[idim]+dx[idim]*ind[idim];
//
//            #if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_
//            if ((x[idim]-PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim])<-eps || (x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim])>eps) {
//              outBoundary=true;
//              break;
//            }
//            #endif
//          }
//
//
//          if (!outBoundary){
//            for (int idim=0; idim<3; idim++) x[idim]=node->xmin[idim]+dx[idim]*ind[idim];
//
//            for (int idim=0; idim<3; idim++)  {
//              #if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_
//              if (fabs(x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim])<eps)
//                x[idim] = PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim];
//              #endif
//
//              tempPnt[idim]=x[idim];
//            }
//
//            *tempDataBuffPnt = node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();
//            tempDataBuffPnt++;
//            tempPnt += 3;
//          }
//
//
//        }// for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)
//      }//if (node->Thread==PIC::ThisThread)
//    }//for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)
//
//
//    iProcessCnt=0;
//    //send the number of points on this proc to all other proc
//    for (int iProc=0; iProc<PIC::nTotalThreads; iProc++){
//      if (iProc!=PIC::ThisThread){
//        MPI_Isend(SendCoordBuff,nPointsToComm[PIC::ThisThread]*3,MPI_DOUBLE,iProc,1,MPI_GLOBAL_COMMUNICATOR,SendCoordList+iProcessCnt);
//        MPI_Isend(SendDataBuffPointerBuff,nPointsToComm[PIC::ThisThread]*sizeof(char *),MPI_BYTE,iProc,2,MPI_GLOBAL_COMMUNICATOR,SendDataBufferPtrList+iProcessCnt);
//        iProcessCnt++;
//      }
//    }
//
//    int nRecvBuffDone=0;
//
//    while (nRecvBuffDone<PIC::nTotalThreads-1){
//      int q, flag, iProc;
//      MPI_Testany(PIC::nTotalThreads-1, RecvPntNumberList, &q, &flag, MPI_STATUS_IGNORE);
//
//      if (flag!=0){
//        //release memoty used by MPI
//        MPI_Wait(RecvPntNumberList+q,MPI_STATUS_IGNORE);
//
//        iProc = q<PIC::ThisThread?q:q+1;
//        if (RecvCoordBuff[q]!=NULL) delete [] RecvCoordBuff[q];
//        RecvCoordBuff[q] = new double [3*nPointsToComm[iProc]];
//        if (RecvDataBuffPointerBuff[q]!=NULL) delete [] RecvDataBuffPointerBuff[q];
//        RecvDataBuffPointerBuff[q] = new char * [nPointsToComm[iProc]];
//        MPI_Irecv(RecvCoordBuff[q],nPointsToComm[iProc]*3,MPI_DOUBLE,iProc,1,MPI_GLOBAL_COMMUNICATOR,RecvCoordList+q);
//        MPI_Irecv(RecvDataBuffPointerBuff[q],nPointsToComm[iProc]*sizeof(char *),MPI_BYTE,iProc,2,MPI_GLOBAL_COMMUNICATOR,RecvDataBufferPtrList+q);
//
//        nRecvBuffDone++;
//      }
//    }
//
//    PIC::Parallel::XYZTree Tree;
//
//    PIC::Parallel::XYZTree::leafNode ** leafNodeArr = new  PIC::Parallel::XYZTree::leafNode * [PIC::nTotalThreads];
//
//    for (iThread=0; iThread<PIC::nTotalThreads; iThread++) {
//      // printf("nPointsToComm[%d]:%d\n",iThread,nPointsToComm[iThread]);
//      leafNodeArr[iThread] = new PIC::Parallel::XYZTree::leafNode [nPointsToComm[iThread]];
//    }
//    //create tree from local thread
//    for (int iPnt=0; iPnt<nPointsToComm[PIC::ThisThread]; iPnt++){
//      leafNodeArr[PIC::ThisThread][iPnt].iThread = PIC::ThisThread;
//      leafNodeArr[PIC::ThisThread][iPnt].DataBuffer = SendDataBuffPointerBuff[iPnt];
//    }
//
//    //build tree from send buffer
//    for (int iPnt=0; iPnt<nPointsToComm[PIC::ThisThread]; iPnt++){
//      Tree.addXNode(SendCoordBuff[3*iPnt],SendCoordBuff[3*iPnt+1],SendCoordBuff[3*iPnt+2],leafNodeArr[PIC::ThisThread]+iPnt);
//    }
//
//    int counterArr[PIC::nTotalThreads-1];
//    for (i=0; i<PIC::nTotalThreads-1; i++) counterArr[i]=0;
//    int counter[2]={0,0};
//    //process other threads
//
//    int nProcDone=0;
//
//    //build tree from recv buffers
//    while (nProcDone<PIC::nTotalThreads-1){
//      int p[2], flag[2];
//
//      if (counter[0]<PIC::nTotalThreads-1) {
//        MPI_Testany(PIC::nTotalThreads-1, RecvCoordList, p, flag, MPI_STATUS_IGNORE);
//
//        //release memory used by MPI
//        if ((flag[0]==true)&&(p[0]!=MPI_UNDEFINED)) MPI_Wait(RecvCoordList+p[0],MPI_STATUS_IGNORE);
//      }
//
//      if (counter[1]<PIC::nTotalThreads-1) {
//        MPI_Testany(PIC::nTotalThreads-1, RecvDataBufferPtrList, p+1, flag+1, MPI_STATUS_IGNORE);
//
//        //release memory used by MPI
//        if ((flag[1]==true)&&(p[1]!=MPI_UNDEFINED)) MPI_Wait(RecvDataBufferPtrList+p[1],MPI_STATUS_IGNORE);
//      }
//
//      for (i=0;i<2; i++){
//        if (flag[i]!=0 && counter[i]<PIC::nTotalThreads-1){
//          int localInd = p[i];
//          counterArr[localInd]++;
//          counter[i]++;
//          if (counterArr[localInd]==2){
//            int iProc=localInd<PIC::ThisThread?localInd:localInd+1;
//
//            for (int iPnt=0; iPnt<nPointsToComm[iProc]; iPnt++){
//              leafNodeArr[iProc][iPnt].iThread = iProc;
//              leafNodeArr[iProc][iPnt].DataBuffer = RecvDataBuffPointerBuff[localInd][iPnt];
//            }
//
//            for (int iPnt=0; iPnt<nPointsToComm[iProc]; iPnt++){
//              Tree.addXNode(RecvCoordBuff[localInd][3*iPnt],RecvCoordBuff[localInd][3*iPnt+1],RecvCoordBuff[localInd][3*iPnt+2],leafNodeArr[iProc]+iPnt);
//            }
//
//            nProcDone++;
//          }
//        }
//      }
//    }//while (nProcDone<PIC::nTotalThreads-1)
//
//    // printf("nProcDone:%d\n",nProcDone);
//    // printf("\n");
//    // printf("thisthread:%d ,counterarr: ",PIC::ThisThread);
//    // for (i=0;i<PIC::nTotalThreads-1;i++) printf(" %d",counterArr[i]);
//    //printf("\n");
//    //Tree.printTree();
//
//    int totalStencil=0;
//    std::list<PIC::Parallel::XYZTree::xNode*>::iterator itX;
//    for (itX=Tree.xList.begin(); itX!=Tree.xList.end(); itX++){
//      std::list<PIC::Parallel::XYZTree::yNode*>::iterator itY;
//      for (itY=(*itX)->yList.begin(); itY!=(*itX)->yList.end(); itY++){
//        //totalStencil += (*itY)->zList.size();
//        std::list<PIC::Parallel::XYZTree::zNode*>::iterator itZ;
//        for (itZ=(*itY)->zList.begin(); itZ!=(*itY)->zList.end(); itZ++){
//          std::list<PIC::Parallel::XYZTree::leafNode*>::iterator itLeaf;
//          if ((*itZ)->leafList.size()>1) totalStencil +=1;
//        }
//      }
//    }
//
//    //allocate the new Stencile Table
//    if (StencilTableLength!=0) {
//      delete [] StencilTable;
//      StencilElementStack.clear();
//    }
//
//    StencilTableLength=totalStencil;
//    StencilTable=new cStencilElement*[totalStencil];
//    //printf("totalStencil:%d, thisthread:%d\n", totalStencil, PIC::ThisThread);
//
//    for (int i=0;i<totalStencil;i++) StencilTable[i]=NULL;
//
//
//    int totalNodes=0;
//    int iSt =0;
//    for (itX=Tree.xList.begin(); itX!=Tree.xList.end(); itX++){
//      std::list<PIC::Parallel::XYZTree::yNode*>::iterator itY;
//      for (itY=(*itX)->yList.begin(); itY!=(*itX)->yList.end(); itY++){
//        std::list<PIC::Parallel::XYZTree::zNode*>::iterator itZ;
//        for (itZ=(*itY)->zList.begin(); itZ!=(*itY)->zList.end(); itZ++){
//          std::list<PIC::Parallel::XYZTree::leafNode*>::iterator itLeaf;
//          int iTh=0;
//          if ((*itZ)->leafList.size()<2) continue;
//
//          StencilTable[iSt]=StencilElementStack.newElement();
//          St=StencilTable[iSt];
//
//          /*
//          if (fabs((*itX)->x)<1e-3 && fabs((*itY)->y)<1e-3 && fabs((*itZ)->z)<1e-3)
//            printf("test ist:%d,x:%e,y:%e,z:%e,StencilLength:%d\n",iSt,(*itX)->x,(*itY)->y,(*itZ)->z,StencilTable[iSt].StencilLength);
//          */
//          int tempThreadId=-1;
//          St->StencilLength=0;
//          St->AssociatedDataPointer=NULL;
//          for (itLeaf=(*itZ)->leafList.begin(); itLeaf!=(*itZ)->leafList.end(); itLeaf++){
//            /*
//            if (fabs((*itX)->x)<1e-3 && fabs((*itY)->y)<1e-3 && fabs((*itZ)->z)<1e-3){
//              printf("test databuffer thisthread:%d, 0,0,0,iThread:%d, DataBuffer:%p\n", PIC::ThisThread,(*itLeaf)->iThread,(void*)((*itLeaf)->DataBuffer));
//            }
//            */
//            if (tempThreadId!=(*itLeaf)->iThread) {
//              St->StencilLength++;
//              St->StencilThreadTable[iTh]=(*itLeaf)->iThread;
//              tempThreadId = (*itLeaf)->iThread;
//              iTh++;
//            }
//
//            if ((*itLeaf)->iThread==PIC::ThisThread) {
//              if (St->AssociatedDataPointer==NULL) {
//                St->AssociatedDataPointer=(*itLeaf)->DataBuffer;
//                if ((*itLeaf)->DataBuffer==NULL) exit(__LINE__,__FILE__,"Error: NULL pointer");
//              }
//              else{
//                //more than 1 data buffer for one point
//                St->localDataPntArr.push_back((*itLeaf)->DataBuffer);
//                if ((*itLeaf)->DataBuffer==NULL) exit(__LINE__,__FILE__,"Error: NULL pointer");
//              }
//            }
//
//          }
//          //determine the center-processing thread
//          int t = iSt%St->StencilLength;
//          int temp = St->StencilThreadTable[0];
//          St->StencilThreadTable[0]=St->StencilThreadTable[t];
//          St->StencilThreadTable[t]=temp;
//          totalNodes+=St->StencilLength;
//          /*
//          if (PIC::ThisThread==0){
//            printf("iSt:%d, stencilLength:%d, threadtable:",iSt, StencilTable[iSt].StencilLength);
//            for (int tt=0; tt<StencilTable[iSt].StencilLength; tt++) printf("%d ",StencilTable[iSt].StencilThreadTable[tt]);
//            printf("x,y,z:%e,%e,%e\n",(*itX)->x,(*itY)->y,(*itZ)->z);
//            printf("\n");
//          }
//          */
//          iSt++;
//        }
//      }
//    }
//
//    //test
//
//    /*
//    char * testPtr = StencilTable[714].AssociatedDataPointer;
//     for (itX=Tree.xList.begin(); itX!=Tree.xList.end(); itX++){
//      std::list<PIC::Parallel::XYZTree::yNode*>::iterator itY;
//      for (itY=(*itX)->yList.begin(); itY!=(*itX)->yList.end(); itY++){
//        std::list<PIC::Parallel::XYZTree::zNode*>::iterator itZ;
//        for (itZ=(*itY)->zList.begin(); itZ!=(*itY)->zList.end(); itZ++){
//          std::list<PIC::Parallel::XYZTree::leafNode*>::iterator itLeaf;
//
//          if ((*itZ)->leafList.size()<2) continue;
//
//          for (itLeaf=(*itZ)->leafList.begin(); itLeaf!=(*itZ)->leafList.end(); itLeaf++){
//            if (fabs((*itX)->x)>2+1e-3 || fabs((*itY)->y)>2+1e-3 || fabs((*itZ)->z)>2+1e-3 ) printf("points out of range!\n");
//            if (fabs((*itX)->x-2)<1e-3 || fabs((*itY)->y-2)<1e-3 || fabs((*itZ)->z-2)<1e-3 ) printf("wrong points !\n");
//            if (testPtr==(*itLeaf)->DataBuffer){
//              printf("repeated:x:%e,y:%e,z:%e, thread id:%d, Jx:%e\n", (*itX)->x,(*itY)->y, (*itZ)->z, PIC::ThisThread,((double *)(StencilTable[714].AssociatedDataPointer+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))[6] );
//            }
//          }
//
//        }
//
//      }
//
//     }//itX
//     */
//
//    // Tree.printTree();
//    //printf("thisthread:%d,iSt:%d,totalNodes:%d\n",PIC::ThisThread,iSt,totalNodes);
//    //update the coundater
//    nMeshModificationCounter=PIC::Mesh::mesh->nMeshModificationCounter;
//    PIC::Mesh::mesh->meshModifiedFlag_CountMeshElements=meshModifiedFlag_CountMeshElements;
//    PIC::Parallel::RebalancingTime+=MPI_Wtime()-StartTime;
//  }
//
//
//  //performe the data exchange session
//  static int * LoadList=NULL;
//  static std::vector<int> ProcessStencilBufferNumber;
//  static std::vector<char*> CopyDataBufferDest;
//  static std::vector<int> CopyDataBufferStencilIndex;
//  static std::vector<int> WholeIndexOfProcessStencil;
//  static char ** processRecvDataBuffer;
//  static char ** copyRecvDataBuffer;
//  static int * wholeStencilIndexOfBuffer=NULL,* processStencilIndexOfBuffer=NULL;
//  static int * StencilFinishedBufferNumber=NULL;
//  static int sumProcessBufferNumber=0;
//  static int totalProcessStencil;
//  static int totalCopyStencil;
//
//  static MPI_Request * processRecvList=NULL, * copyRecvList=NULL;
//  static MPI_Request * processSendList=NULL, * copySendList=NULL;
//
//
//  //==============================================   BEGINING OF THE DATA EXCHANGE LOOP ========================
//
//  int iStencilStep=1000;
//  int iStencilStart=0,iStencilFinish=iStencilStep;
//
//  for (iStencilStart=0;iStencilStart<StencilTableLength;iStencilStart=((iStencilStart+iStencilStep<StencilTableLength) ? (iStencilStart+iStencilStep) : StencilTableLength)) {
//    iStencilFinish=iStencilStart+iStencilStep;
//    if (iStencilFinish>StencilTableLength) iStencilFinish=StencilTableLength;
//
//  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
//
//    if (LoadList!=NULL) delete [] LoadList;
//    LoadList=new int [PIC::Mesh::mesh->nTotalThreads];
//    for (int i=0;i<PIC::Mesh::mesh->nTotalThreads;i++) LoadList[i]=0;
//
//    sumProcessBufferNumber = 0;
//    //print out the thread distribution of process stencils
//    for (int iStencil=iStencilStart;iStencil<iStencilFinish;iStencil++) if ((St=StencilTable[iStencil])!=NULL) {
//      //if (StencilTable[iStencil].StencilLength>1) {
//      int center=St->StencilThreadTable[0];
//      if (St->StencilLength>1) LoadList[center]++;
//      if (center==PIC::ThisThread) sumProcessBufferNumber += St->StencilLength-1;
//      //}
//    }
//
//    int sumLoad=0;
//    for (int i=0;i<PIC::Mesh::mesh->nTotalThreads;i++) {
//      //  printf("thread:%d, loadlist[%d],%d",PIC::ThisThread,i,LoadList[i]);
//      sumLoad+=LoadList[i];
//    }
//
//    //printf("thread:%d, StencilTableLength:%d\n",PIC::ThisThread,StencilTableLength);
//    //printf("\n");
//    //    if (sumLoad!=StencilTableLength) exit(__LINE__,__FILE__,"Error: Stencils containing only one thread/proc are in the StencilTable");
//
//    if (wholeStencilIndexOfBuffer!=NULL) delete [] wholeStencilIndexOfBuffer;
//    if (processStencilIndexOfBuffer!=NULL) delete [] processStencilIndexOfBuffer;
//    wholeStencilIndexOfBuffer=new int[sumProcessBufferNumber];
//    processStencilIndexOfBuffer =new int[sumProcessBufferNumber];
//
//    int processBufferIndex=0;
//
//    ProcessStencilBufferNumber.clear();
//    CopyDataBufferDest.clear();
//    CopyDataBufferStencilIndex.clear();
//    WholeIndexOfProcessStencil.clear();
//
//    for (int iStencil=iStencilStart;iStencil<iStencilFinish;iStencil++) if ((St=StencilTable[iStencil])!=NULL) {
//      if (St->AssociatedDataPointer){// means this thread is involved
//        if (St->StencilThreadTable[0]==PIC::ThisThread) {
//          if (St->StencilLength>1){
//
//            ProcessStencilBufferNumber.push_back(St->StencilLength-1);
//            WholeIndexOfProcessStencil.push_back(iStencil);
//          }
//          for (int iThread=0; iThread<St->StencilLength-1;iThread++){
//            //process buffer number equals StencilLength-1
//            wholeStencilIndexOfBuffer[processBufferIndex]=iStencil; //index in the whole stensil list
//            processStencilIndexOfBuffer[processBufferIndex]=ProcessStencilBufferNumber.size()-1;
//            processBufferIndex++;
//          }
//        }else{
//          CopyDataBufferDest.push_back(St->AssociatedDataPointer);
//          CopyDataBufferStencilIndex.push_back(iStencil);
//        }//else
//      }
//    }
//    //printf("thread id:%d, processBufferIndex:%d\n", PIC::ThisThread, processBufferIndex);
//
//    if (processRecvDataBuffer!=NULL) {
//      delete [] processRecvDataBuffer[0];
//      delete [] processRecvDataBuffer;
//    }
//
//    processRecvDataBuffer= new char * [sumProcessBufferNumber];
//    processRecvDataBuffer[0] = new char [sumProcessBufferNumber*PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength];
//    for (int i=1;i<sumProcessBufferNumber;i++) {
//      processRecvDataBuffer[i]=i*PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength+processRecvDataBuffer[0];
//    }
//
//    if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
//
//      if (copyRecvDataBuffer!=NULL){
//        delete [] copyRecvDataBuffer[0];
//        delete [] copyRecvDataBuffer;
//      }
//
//      copyRecvDataBuffer=new char * [CopyDataBufferDest.size()];
//      copyRecvDataBuffer[0]=new char [CopyDataBufferDest.size()*PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength];
//      for (int i=1;i<CopyDataBufferDest.size();i++) copyRecvDataBuffer[i]=copyRecvDataBuffer[0]+i*PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength;
//    }
//
//
//    if (processSendList!=NULL) delete[] processSendList;
//    processSendList = new MPI_Request[CopyDataBufferDest.size()];
//    if (processRecvList!=NULL) delete[] processRecvList;
//    processRecvList = new MPI_Request[sumProcessBufferNumber];
//    if (copyRecvList!=NULL) delete [] copyRecvList;
//    copyRecvList = new MPI_Request[CopyDataBufferDest.size()];
//    if (copySendList!=NULL) delete [] copySendList;
//    copySendList = new MPI_Request[sumProcessBufferNumber];
//
//    //printf("sumProcessBufferNumber:%d, copyBufferSize:%d\n",sumProcessBufferNumber,CopyDataBufferDest.size());
//
//    if (StencilFinishedBufferNumber!=NULL) delete [] StencilFinishedBufferNumber;
//    StencilFinishedBufferNumber=new int [LoadList[PIC::ThisThread]];
//    totalProcessStencil= LoadList[PIC::ThisThread];
//    totalCopyStencil = CopyDataBufferDest.size();
//
//  }//if (globalMeshChangeFlag!=0 && PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL)
//
//  //printf("thread id:%d, sumProcessBufferNumber:%d,totalProcessStencil:%d,WholeIndexOfProcessStencil size:%d\n",PIC::ThisThread,sumProcessBufferNumber,totalProcessStencil,WholeIndexOfProcessStencil.size());
//
//  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
//    //clear counter
//    for (int i=0; i<totalProcessStencil;i++) StencilFinishedBufferNumber[i]=0;
//
//    //1. combine 'corner' node data
//    int iProcessBuffer=0, nProcessBufferDone=0, nProcessStencilDone=0;
//    int iCopyBuffer=0, nCopyBufferDone=0;
//    int iProcessSendBuffer=0;
//
//
//
//    for (int iStencil=iStencilStart;iStencil<iStencilFinish;iStencil++) if ((St=StencilTable[iStencil])!=NULL) {//loop through stencil table
//      if (St->AssociatedDataPointer) { // this thread is involved
//        int iThread;
//        //there are more that one MPI processes that contributed to the state vector of the corner node
//
//        if (PIC::ThisThread==St->StencilThreadTable[0]) {
//          //if this thread will do the center processing role and  combine data from all other MPI processes
//          //recv data request start
//          //process local points
//          std::vector<char *>::iterator it;
//          for (it=St->localDataPntArr.begin();it!=St->localDataPntArr.end();it++){
//            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(St->AssociatedDataPointer,*it);
//          }
//
//          for (iThread=1;iThread<St->StencilLength;iThread++) {
//            MPI_Irecv(processRecvDataBuffer[iProcessBuffer],PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,St->StencilThreadTable[iThread],iStencil,MPI_GLOBAL_COMMUNICATOR,processRecvList+iProcessBuffer);
//            iProcessBuffer++;
//          }
//
//        }else{
//          //this thread will be copyThread, send data to the processThread,
//          //and recv the processed data from the processThread
//          std::vector<char *>::iterator it;
//          for (it=St->localDataPntArr.begin();it!=St->localDataPntArr.end();it++){
//            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(St->AssociatedDataPointer,*it);
//          }
//          MPI_Isend(St->AssociatedDataPointer,PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,St->StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,processSendList+iCopyBuffer);
//
//          //start the recv data request
//          if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
//            MPI_Irecv(copyRecvDataBuffer[iCopyBuffer],PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,St->StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,copyRecvList+iCopyBuffer);
//          }else{
//            MPI_Irecv(CopyDataBufferDest[iCopyBuffer],PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,St->StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,copyRecvList+iCopyBuffer);
//          }
//          iCopyBuffer++;
//          // }//for
//        }//else
//
//
//        for (int i=nProcessBufferDone; i<iProcessBuffer;i++) {
//          int q, flag;
//          // printf("test7 pass\n");
//          //if (nProcessBufferDone<sumProcessBufferNumber)
//          MPI_Testany(iProcessBuffer, processRecvList, &q, &flag, MPI_STATUS_IGNORE);
//          // printf("test8 pass\n");
//          //printf("nProcessBufferDone:%d, sumProcessBufferNumber:%d\n",nProcessBufferDone,sumProcessBufferNumber);
//
//          if (flag!=0 && q!=MPI_UNDEFINED){
//            //process buffer of index q is done
//            int processStencilIndex = processStencilIndexOfBuffer[q];
//            int wholeStencilIndex = wholeStencilIndexOfBuffer[q];
//
//            //release memoty used by MPI
//            MPI_Wait(processRecvList+q,MPI_STATUS_IGNORE);
//
//            //printf("thread id:%d,bufferid:%d processed, nProcessBufferDone:%d/total buffer:%d,\n",PIC::ThisThread,q,nProcessBufferDone,sumProcessBufferNumber);
//
//            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[wholeStencilIndex]->AssociatedDataPointer,processRecvDataBuffer[q]);
//
//            nProcessBufferDone++;
//            StencilFinishedBufferNumber[processStencilIndex]++;
//            //if  finished one stencil processing, send it back
//            if (StencilFinishedBufferNumber[processStencilIndex]==ProcessStencilBufferNumber[processStencilIndex]){
//              for (iThread=1;iThread<StencilTable[wholeStencilIndex]->StencilLength;iThread++) {
//                MPI_Isend(StencilTable[wholeStencilIndex]->AssociatedDataPointer,PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[wholeStencilIndex]->StencilThreadTable[iThread],wholeStencilIndex,MPI_GLOBAL_COMMUNICATOR,copySendList+iProcessSendBuffer);
//                iProcessSendBuffer++;
//              }
//              nProcessStencilDone++;
//            }
//
//          }
//        }//for (int i=nProcessBufferDone; i<iProcessBuffer;i++)
//
//        for (int i=nCopyBufferDone; i<iCopyBuffer;i++) {
//          int q, flag;
//
//          MPI_Testany(iCopyBuffer, copyRecvList, &q, &flag, MPI_STATUS_IGNORE);
//
//          if (flag!=0 && q!=MPI_UNDEFINED){
//            //release memoty used by MPI
//            MPI_Wait(copyRecvList+q,MPI_STATUS_IGNORE);
//
//            //int wholeStencilIndex=CopyDataBufferStencilIndex[q];
//            if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
//              PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(CopyDataBufferDest[q],copyRecvDataBuffer[q]);
//            }
//            //if CopyCornerNodeAssociatedData is not defined, the MPI_Irecv already sent the data to the destination state vector
//            nCopyBufferDone++;
//          }
//        }//for (int i=nCopyBufferDone; i<iCopyBuffer;i++)
//      }
//    }//for (iStencil=0;iStencil<StencilTableLength;iStencil++)
//
//
//    //printf("out of for-loop\n");
//
//    //make sure communication is finished
//    while (nProcessStencilDone<totalProcessStencil||nCopyBufferDone<totalCopyStencil){
//
//      for (int i=nProcessBufferDone; i<iProcessBuffer;i++) {
//        int q, flag;
//
//        MPI_Testany(iProcessBuffer, processRecvList, &q, &flag, MPI_STATUS_IGNORE);
//
//        if (flag!=0 &&  q!=MPI_UNDEFINED){
//          //process buffer of index q is done
//          int processStencilIndex = processStencilIndexOfBuffer[q];
//          int wholeStencilIndex = wholeStencilIndexOfBuffer[q];
//
//          //release memoty used by MPI
//          MPI_Wait(processRecvList+q,MPI_STATUS_IGNORE);
//
//          // printf("thread id:%d,bufferid:%d processed, nProcessBufferDone:%d/total buffer:%d,\n",PIC::ThisThread,q,nProcessBufferDone,sumProcessBufferNumber);
//          /*
//            if (wholeStencilIndex==714){
//            printf("ist:%d,in processing thread id: %d, size:%d, (while) before Jx:%e\n",wholeStencilIndex,PIC::ThisThread,StencilTable[wholeStencilIndex].localDataPntArr.size(), ((double *)(StencilTable[wholeStencilIndex].AssociatedDataPointer+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))[6]);
//            }
//          */
//          PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[wholeStencilIndex]->AssociatedDataPointer,processRecvDataBuffer[q]);
//
//          nProcessBufferDone++;
//          StencilFinishedBufferNumber[processStencilIndex]++;
//
//          if (StencilFinishedBufferNumber[processStencilIndex]==ProcessStencilBufferNumber[processStencilIndex]){
//            for (iThread=1;iThread<StencilTable[wholeStencilIndex]->StencilLength;iThread++) {
//              MPI_Isend(StencilTable[wholeStencilIndex]->AssociatedDataPointer,PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[wholeStencilIndex]->StencilThreadTable[iThread],wholeStencilIndex,MPI_GLOBAL_COMMUNICATOR,copySendList+iProcessSendBuffer);
//              iProcessSendBuffer++;
//            }
//            nProcessStencilDone++;
//          }
//        }
//      }//for (int i=nProcessBufferDone; i<iProcessBuffer;i++)
//
//      for (int i=nCopyBufferDone; i<iCopyBuffer;i++) {
//        int q, flag;
//
//        MPI_Testany(iCopyBuffer, copyRecvList, &q, &flag, MPI_STATUS_IGNORE);
//
//        if (flag!=0 && q!=MPI_UNDEFINED){
//          //release memoty used by MPI
//          MPI_Wait(copyRecvList+q,MPI_STATUS_IGNORE);
//
//          if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
//            PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(CopyDataBufferDest[q],copyRecvDataBuffer[q]);
//          }
//          //if CopyCornerNodeAssociatedData is not defined, the MPI_Irecv already sent the data to the destination state vector
//          nCopyBufferDone++;
//        }
//      }//for (int i=nCopyBufferDone; i<iCopyBuffer;i++)
//    }//while (nProcessStencilDone<totalProcessStencil||nCopyBufferDone<totalCopyStencil)
//
//
//
//    //sync the local data buffers
//    for (int iSt=iStencilStart;iSt<iStencilFinish;iSt++) if ((St=StencilTable[iSt])!=NULL) {
//      if (St->AssociatedDataPointer) {
//
//        std::vector<char *>::iterator it;
//        if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
//          for (it=St->localDataPntArr.begin();it!=St->localDataPntArr.end();it++){
//            PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(*it,St->AssociatedDataPointer);
//          }
//        }else{
//          for (it=St->localDataPntArr.begin();it!=St->localDataPntArr.end();it++)
//            memcpy(*it,St->AssociatedDataPointer,PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength);
//        }
//      }
//    }
//
//    MPI_Waitall(iProcessBuffer,processRecvList,MPI_STATUSES_IGNORE);
//    MPI_Waitall(iCopyBuffer,copyRecvList,MPI_STATUSES_IGNORE);
//    MPI_Waitall(iCopyBuffer,processSendList,MPI_STATUSES_IGNORE);
//    MPI_Waitall(iProcessBuffer,copySendList,MPI_STATUSES_IGNORE);
//
//    }//if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL)
//
//  }
//  //===================================  END OF THE DATA EXCHANGE LOOP ========================================================
//
//}
  
void PIC::Parallel::ProcessCornerBlockBoundaryNodes_new() {
  if (CornerBlockBoundaryNodes::ActiveFlag==false) return;
  // Creating a BPManager for each type of communications and passing
  // it as an argument to this function would be a better choice,
  // instead of using BPManager as a global variable. --Yuxi
  ProcessBlockBoundaryNodes(BPManager);
}
//-------------------------------------

void PIC::Parallel::ProcessCenterBlockBoundaryNodes_new() {
   if (CenterBlockBoundaryNodes::ActiveFlag==false) return;
   ProcessBlockBoundaryNodes(BPManager);   
}
//-------------------------------------

void PIC::Parallel::ProcessBlockBoundaryNodes(BoundaryProcessManager &mgr) {

  /*
    1. The algorithm:
    On each MPI, this function 1) loops through all the local blocks, counts
    the corner/center that will be sent to other MPIs, 2) allocates the send
    and receive buffer, 3) packs the data, 4) sends the data, 5) receives the
    data and 6) adds the data to local corners/centers.

    2. If the nearby blocks are on the same MPI, they share the edge/corner 
    memory. This function 'pretend' they do not share the memory. The 
    information of these points may be sent/received several times. Function
    get_n_share() is used to counts how many times the information is sent
    and received, and the data is 'split' before they are added to 
    the local nodes. 

    3. The data in the buffer:
    receive_node_ID, iDir, jDir, kDir, coef, data.. iDir, jDir, kDir, 
    coef, data... receive_node_ID, iDir....

   */
  
  const bool isCorner = mgr.isCorner;
  const int pointBufferSize = mgr.pointBufferSize;
  const int nByteID = sizeof(cAMRnodeID) + 3*sizeof(int);

  // For a face, we divided the face points into pure face points,
  // pure edge points ( x 4), and pure corner points ( x 4); 
  // For a dege: pure dege + pure corners x 2		
  const int nSubDomain[3]={1,3,9};
		
  //----- lambda begin------------------------------
  auto get_neib_node = [](cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node, int i, int j, int k){
    // i/j/k = -1, 0, or 1
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * nodeNeib = NULL;
    int nZeros = 0; 
    if(i==0) nZeros++;
    if(j==0) nZeros++;
    if(k==0) nZeros++; 

    if(nZeros == 0){
      // Corners 
      int idx = -1; 

      if (i==-1) i=0;
      if (j==-1) j=0;
      if (k==-1) k=0;
      idx = i+2*(j+2*k);
      
      nodeNeib = node->GetNeibCorner(idx,PIC::Mesh::mesh);      
    }else if(nZeros == 1){
      // Edges
      int idx = -1; 
      if(i==0){
	if(j == -1 && k == -1) idx = 0; 
	if(j ==  1 && k == -1) idx = 1; 
	if(j ==  1 && k ==  1) idx = 2; 
	if(j == -1 && k ==  1) idx = 3; 
      }else if(j==0){
	if(i == -1 && k == -1) idx = 4; 
	if(i ==  1 && k == -1) idx = 5; 
	if(i ==  1 && k ==  1) idx = 6; 
	if(i == -1 && k ==  1) idx = 7;       
      }else{
	if(i == -1 && j == -1) idx = 8; 
	if(i ==  1 && j == -1) idx = 9; 
	if(i ==  1 && j ==  1) idx = 10; 
	if(i == -1 && j ==  1) idx = 11; 
      }

      nodeNeib = node->GetNeibEdge(idx,0,PIC::Mesh::mesh); 
      
    }else if(nZeros == 2){
      // Faces
      if(i == -1) nodeNeib = node->GetNeibFace(0,0,0,PIC::Mesh::mesh);
      if(i ==  1) nodeNeib = node->GetNeibFace(1,0,0,PIC::Mesh::mesh);

      if(j == -1) nodeNeib = node->GetNeibFace(2,0,0,PIC::Mesh::mesh);
      if(j ==  1) nodeNeib = node->GetNeibFace(3,0,0,PIC::Mesh::mesh);

      if(k == -1) nodeNeib = node->GetNeibFace(4,0,0,PIC::Mesh::mesh);
      if(k ==  1) nodeNeib = node->GetNeibFace(5,0,0,PIC::Mesh::mesh);
    }else{
      nodeNeib = NULL; 
    }

    return nodeNeib; 
  };
  //---------------------------------------------------------------

  auto get_points_number_all = [](int i, int j,int k, bool isCorner=true){    
    // Find the number points at the (i,j,k) direction.
    int nI = _BLOCK_CELLS_X_; 
    int nJ = _BLOCK_CELLS_Y_;
    int nK = _BLOCK_CELLS_Z_;

    if(isCorner){
      nI++; nJ++; nK++; 
    }
    
    if(i != 0) nI=1;
    if(j != 0) nJ=1; 
    if(k != 0) nK=1; 

    return nI*nJ*nK; 
  };
  //---------------------------------------------------------------

  auto get_points_number = [](int i, int j,int k, bool isCorner=true){    
    // Find the number points at the (i,j,k) direction.
    // For a face, this function only counts the 'pure' face points,
    // and excludes the edge and corner points. 
    // For a edge, this function only counts the pure edge points,
    // and excludes the corner points. 
    
    int nI = _BLOCK_CELLS_X_ - 2; 
    int nJ = _BLOCK_CELLS_Y_ - 2;
    int nK = _BLOCK_CELLS_Z_ - 2;

    if(isCorner){
      nI++; nJ++; nK++; 
    }
    
    if(i != 0) nI=1;
    if(j != 0) nJ=1; 
    if(k != 0) nK=1; 

    return nI*nJ*nK; 
  };
  //---------------------------------------------------------------
      
  auto find_loop_range = [] (int const i, int const j, int const k, int& iMin, 
			     int& jMin, int& kMin, int& iMax, int& jMax, int& kMax, 
			     bool isCorner=true, bool isRecv = true, 
			     int const iNeib=-2, int const jNeib=-2, int const kNeib=-2){
   /* Find the loop range for the 'pure' faces, 'pure' edges and corners.
      
      1. For the cell centers, the indices for send and receive are different. 
      For example, for i=1, j=0, and k=0, (a) iMin = iMax = _BLOCK_CELL_X_, which
      represents the ghost cells, for send, (b) but iMin = iMax = _BLOCK_CELL_X_-1
      for receive. 

      2. Case 1: i=1, j=1, k=0; iNeib=1, jNeib=1, kNeib=0
         Case 2: i=1, j=1, k=0; iNeib=1, jNeib=0, kNeib=0
	 The send indices for the two cases above are different. Case 1 sends
	 the 'real edges' of the sending block, while case 2 actually sends the
	 face ghost cells of the sending block. 	 
   */
			   
    bool isCenterSend = (!isCorner) && (!isRecv);

    iMin = 0;
    jMin = 0;
    kMin = 0; 
    iMax = _BLOCK_CELLS_X_ - 1; 
    jMax = _BLOCK_CELLS_Y_ - 1;  
    kMax = _BLOCK_CELLS_Z_ - 1;

    if(isCorner){
      iMax++; jMax++; kMax++;
    }
    
    //------------i--------------------
    if(i ==  1){
      if(isCenterSend && fabs(iNeib)!=0) iMax++; 
	iMin = iMax;
    }

    if(i == -1){
      if(isCenterSend && fabs(iNeib) !=0) iMin--; 
      iMax = iMin; 
    }

    if(i ==  0){
      iMin++;
      iMax--;
    }

    //------------j--------------------
    if(j ==  1){
      if(isCenterSend && fabs(jNeib) !=0) jMax++; 
      jMin = jMax;
    }
    if(j == -1){
      if(isCenterSend && fabs(jNeib) !=0) jMin--; 
      jMax = jMin; 
    }
    if(j ==  0){
      jMin++;
      jMax--;
    }
    
    //------------k--------------------
    if(k ==  1){
      if(isCenterSend && fabs(kNeib) !=0) kMax++; 
      kMin = kMax;
    }
    if(k == -1){
      if(isCenterSend && fabs(kNeib) !=0) kMin--; 
      kMax = kMin; 
    }
    if(k ==  0){
      kMin++; 
      kMax--;
    }
    
  };
  //---------------------------------------------------------------
  
  auto add_map_val = [](std::map<int, int> &mapInt2Int, int key, int val){
    // If the key does not exist: map[key] = val; 
    // If the key already exists: map[key] += val; 
    if( mapInt2Int.find(key) != mapInt2Int.end() ){
      mapInt2Int[key] += val;
    }else{
      mapInt2Int[key] = val; 
    }    
  };
  //---------------------------------------------------------------

  auto is_neib_on_this_thread = [&get_neib_node](cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node, 
				   int i, int j, int k){
    bool isOn = false; 
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode; 
    neibNode = get_neib_node(node, i, j, k);
    if(neibNode != NULL && neibNode->IsUsedInCalculationFlag ){
      int neibThread = neibNode->Thread; 		
      if(neibThread == PIC::ThisThread){
	isOn = true;
      }
    }
    
    return isOn; 
  };
  //---------------------------------------------------------------

  auto get_n_share = [&is_neib_on_this_thread](cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node, int i, int j, int k){
    int nShare=1;
    
    const int nZeros = 3 - fabs(i) - fabs(j) - fabs(k);       

    if(nZeros==2){
      // face
      if(is_neib_on_this_thread(node, i, j, k)) nShare++;      
    }else if(nZeros==1){
      // edge
      // Two face neighbors and one corner neiber may share this edge. 
      if(is_neib_on_this_thread(node, 0, j, k)) nShare++;      
      if(is_neib_on_this_thread(node, i, 0, k)) nShare++;      
      if(is_neib_on_this_thread(node, i, j, 0)) nShare++;      
      
    }else{
      // corner
      // Three face neighbors, three edge neighbors and one corner neighbor
      // may share this corner. 

      // Corner neighbor. 
      if(is_neib_on_this_thread(node, i, j, k)) nShare++;      

      // Edge neighbors. 
      if(is_neib_on_this_thread(node, 0, j, k)) nShare++;      
      if(is_neib_on_this_thread(node, i, 0, k)) nShare++;      
      if(is_neib_on_this_thread(node, i, j, 0)) nShare++;      
      

      // Face neighbors. 
      if(is_neib_on_this_thread(node, i, 0, 0)) nShare++;            
      if(is_neib_on_this_thread(node, 0, j, 0)) nShare++;            
      if(is_neib_on_this_thread(node, 0, 0, k)) nShare++;            
    }

    return nShare; 
  };

  //----------lambda end---------------------------------------------


  std::map<int, int> mapProc2BufferSize; 
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node, *neibNode;
  
  //---- count send/recv buffer size ------------------
  for (int iNodeIdx=0;iNodeIdx<PIC::DomainBlockDecomposition::nLocalBlocks;iNodeIdx++) {
    node=PIC::DomainBlockDecomposition::BlockTable[iNodeIdx];
    if (node->Thread==PIC::ThisThread && node->block) { 
      for(int iNeib = -1; iNeib <= 1; iNeib++)
	for(int jNeib = -1; jNeib <= 1; jNeib++)
	  for(int kNeib = -1; kNeib <= 1; kNeib++){
	    neibNode = get_neib_node(node, iNeib, jNeib, kNeib);
	    if(neibNode != NULL && neibNode->IsUsedInCalculationFlag){
	      int neibThread = neibNode->Thread; 		
	      if(neibThread != PIC::ThisThread){
		// At most 2 zeros. 
		const int nZeros = 3 - fabs(iNeib) - fabs(jNeib) - fabs(kNeib);

		// In the data buffer: recvNodeId, -iNeib, -jNeib, -kNeib, coeff, data.....
		const int nByteData = get_points_number_all(iNeib, jNeib, kNeib, isCorner)*pointBufferSize;
		const int nByteSend = nByteData + nByteID + sizeof(double)*nSubDomain[nZeros];
		add_map_val(mapProc2BufferSize, neibNode->Thread, nByteSend);
	      }
	    }
	  }// for kNeib	      
    }// if node
  }// for iNodeIdx


  // -------------prepare buffer begin-------------------------------
  int nProcSendTo, *procSendToList=nullptr, *sendBufferSizePerProc=nullptr;
  int nProcRecvFrom, *procRecvFromList=nullptr, *recvBufferSizePerProc=nullptr;
  char **sendBuffer = nullptr;  
  char **recvBuffer = nullptr;
  MPI_Request *reqList;
  MPI_Status *statusList;

  nProcSendTo = mapProc2BufferSize.size(); 
  procSendToList = new int[nProcSendTo];
  sendBufferSizePerProc = new int[nProcSendTo];
  sendBuffer = new char*[nProcSendTo];

  // For uniform grid, send and receive is asymmetric. 
  nProcRecvFrom = nProcSendTo; 
  procRecvFromList = new int[nProcRecvFrom];
  recvBufferSizePerProc = new int[nProcRecvFrom];
  recvBuffer = new char*[nProcRecvFrom];

  reqList = new MPI_Request[nProcSendTo+nProcRecvFrom];
  statusList = new MPI_Status[nProcSendTo+nProcRecvFrom];

  int iCount = 0; 
  for(auto const& x : mapProc2BufferSize ){
    procSendToList[iCount] = x.first; 
    sendBufferSizePerProc[iCount] = x.second; 

    procRecvFromList[iCount] = x.first;
    recvBufferSizePerProc[iCount] = x.second; 
    
    iCount++;
  }

  for(int i=0; i<nProcSendTo; i++){
    sendBuffer[i] = new char[sendBufferSizePerProc[i]];
  }
  
  for(int i=0; i<nProcRecvFrom; i++){
    recvBuffer[i] = new char[recvBufferSizePerProc[i]];
  }
  // -------------prepare buffer end-------------------------------
    
  int nRequest=0; 
  for(int iRecv=0; iRecv<nProcRecvFrom; iRecv++){
    MPI_Irecv(recvBuffer[iRecv], recvBufferSizePerProc[iRecv], MPI_CHAR, 
	      procRecvFromList[iRecv],0, MPI_GLOBAL_COMMUNICATOR, &reqList[nRequest++]);
  }
  
  for(int ii = 0; ii < nProcSendTo; ii++ ){
    // Packing data for each target MPI. 
    int iTarget = procSendToList[ii];
    int offset=0; 
    for (int iNodeIdx=0;iNodeIdx<PIC::DomainBlockDecomposition::nLocalBlocks;iNodeIdx++) {
      node=PIC::DomainBlockDecomposition::BlockTable[iNodeIdx];
      if (node->Thread==PIC::ThisThread && node->block) { 

	for(int iNeib = -1; iNeib <= 1; iNeib++) 
	  for(int jNeib = -1; jNeib <= 1; jNeib++)
	    for(int kNeib = -1; kNeib <= 1; kNeib++){
	      neibNode = get_neib_node(node, iNeib, jNeib, kNeib);
	      if(neibNode != NULL && neibNode->IsUsedInCalculationFlag){
		int neibThread = neibNode->Thread; 		
		if(neibThread == iTarget){		  
		  *((cAMRnodeID*)(sendBuffer[ii]+offset)) = neibNode->AMRnodeID;
		  offset += sizeof(cAMRnodeID);		  

		  *((int*)(sendBuffer[ii]+offset)) =  -iNeib; offset +=sizeof(int);
		  *((int*)(sendBuffer[ii]+offset)) =  -jNeib; offset +=sizeof(int);
		  *((int*)(sendBuffer[ii]+offset)) =  -kNeib; offset +=sizeof(int);


		  const int iNeibMin = (iNeib==0)? -1:iNeib;
		  const int jNeibMin = (jNeib==0)? -1:jNeib;
		  const int kNeibMin = (kNeib==0)? -1:kNeib;

		  const int iNeibMax = (iNeib==0)?  1:iNeib;
		  const int jNeibMax = (jNeib==0)?  1:jNeib;
		  const int kNeibMax = (kNeib==0)?  1:kNeib;
		    
		  const bool isRecv = false; 
		  for(int iDir=iNeibMin; iDir<=iNeibMax; iDir++)
		    for(int jDir=jNeibMin; jDir<=jNeibMax; jDir++)
		      for(int kDir=kNeibMin; kDir<=kNeibMax; kDir++){
			int iMin, jMin, kMin, iMax, jMax, kMax;	 
			find_loop_range(iDir, jDir, kDir, iMin, jMin, kMin, iMax, jMax, kMax, isCorner, isRecv, iNeib, jNeib, kNeib);

			double coef = 1./get_n_share(node, iDir, jDir, kDir); 
			
			*((double*)(sendBuffer[ii]+offset)) = coef; offset +=sizeof(double);
		  		 
			for(int i=iMin; i<=iMax; i++)
			  for(int j=jMin; j<=jMax; j++)
			    for(int k=kMin; k<=kMax; k++){
			      mgr.copy_node_to_buffer(node, i, j, k, sendBuffer[ii]+offset);

			      offset += pointBufferSize;
			    }
		  
		      }
		 		 		
		}
	      }
	    }// for kNeib	      
      }// if node
    }// for iNodeIdx

  }// for ii

  


  for(int iSend=0; iSend<nProcSendTo; iSend++){
    MPI_Isend(sendBuffer[iSend], sendBufferSizePerProc[iSend], MPI_CHAR, procSendToList[iSend],0, MPI_GLOBAL_COMMUNICATOR, &reqList[nRequest++]);
  }

  MPI_Waitall(nRequest, &reqList[0], &statusList[0]);

  // Unpacking data
  for(int ii=0; ii<nProcRecvFrom; ii++){
    int offset=0;
    cAMRnodeID nodeId;
    int iNeib, jNeib, kNeib;
    double coef; 
    while(offset < recvBufferSizePerProc[ii]){
      nodeId = *((cAMRnodeID*)(recvBuffer[ii]+offset)); offset += sizeof(cAMRnodeID);
      node = PIC::Mesh::mesh->findAMRnodeWithID(nodeId);

      iNeib = *((int*)(recvBuffer[ii]+offset)); offset += sizeof(int);
      jNeib = *((int*)(recvBuffer[ii]+offset)); offset += sizeof(int);
      kNeib = *((int*)(recvBuffer[ii]+offset)); offset += sizeof(int);      

      const int iNeibMin = (iNeib==0)? -1:iNeib;
      const int jNeibMin = (jNeib==0)? -1:jNeib;
      const int kNeibMin = (kNeib==0)? -1:kNeib;
      
      const int iNeibMax = (iNeib==0)?  1:iNeib;
      const int jNeibMax = (jNeib==0)?  1:jNeib;
      const int kNeibMax = (kNeib==0)?  1:kNeib;
      
      for(int iDir=iNeibMin; iDir<=iNeibMax; iDir++)
	for(int jDir=jNeibMin; jDir<=jNeibMax; jDir++)
	  for(int kDir=kNeibMin; kDir<=kNeibMax; kDir++){
	    coef = *((double*)(recvBuffer[ii]+offset)); offset +=sizeof(double);

	    // For cell centers, the receive cells are all physical cells, and
	    // they can NOT be shared by nearby blocks. 
	    if(isCorner) coef = coef/get_n_share(node, iDir, jDir, kDir);	    
      
	    const bool isRecv=true;
	    int iMin, jMin, kMin, iMax, jMax, kMax;		 
	    find_loop_range(iDir, jDir, kDir, iMin, jMin, kMin, iMax, jMax, kMax, isCorner, isRecv, iNeib, jNeib, kNeib);

	    for(int i=iMin; i<=iMax; i++)
	      for(int j=jMin; j<=jMax; j++)
		for(int k=kMin; k<=kMax; k++){
		  mgr.add_buffer_to_node(node, i, j, k, recvBuffer[ii]+offset, coef);
		  offset += pointBufferSize;
		}
	  }

    }// while 
  } // for ii


  {
    // Delete arrays.
    for(int i=0; i<nProcSendTo; i++){
      delete [] sendBuffer[i];
    }    
    for(int i=0; i<nProcRecvFrom; i++){
      delete [] recvBuffer[i];
    }
    delete [] procSendToList;
    delete [] sendBufferSizePerProc;
    delete [] sendBuffer;

    delete [] procRecvFromList;
    delete [] recvBufferSizePerProc;
    delete [] recvBuffer;
    
    delete [] reqList; 
    delete [] statusList;
  }
  
}


void PIC::Parallel::ProcessCenterBlockBoundaryNodes() {
  int iThread,i,j,k,iface;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

  if (CenterBlockBoundaryNodes::ActiveFlag==false) return;

  
  const int iFaceMin[6]={-1,_BLOCK_CELLS_X_-1,             -1,             -1,             -1,             -1};
  const int iFaceMax[6]={ 0,  _BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_};

  const int jFaceMin[6]={             -1,             -1,-1,_BLOCK_CELLS_Y_-1,             -1,             -1};
  const int jFaceMax[6]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_, 0,  _BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_};

  const int kFaceMin[6]={             -1,             -1,             -1,             -1,-1,_BLOCK_CELLS_Z_-1};
  const int kFaceMax[6]={_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_, 0,  _BLOCK_CELLS_Z_};
  
  struct cStencilElement {
    int StencilLength; //represent distinct thread number
    int StencilThreadTable[80];
    char *AssociatedDataPointer;
    std::vector<char *> localDataPntArr; //save local data buffer pointers
  };
  
  static int StencilTableLength=0;
  static cStencilElement *StencilTable=NULL;

  //generate a new stencil table
  static int nMeshModificationCounter=-1;

  //determine whether the mesh/domain decomposition have been changed
  int localMeshChangeFlag,globalMeshChangeFlag;

  localMeshChangeFlag=(nMeshModificationCounter==PIC::Mesh::mesh->nMeshModificationCounter) ? 0 : 1;
  MPI_Allreduce(&localMeshChangeFlag,&globalMeshChangeFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  static int * nPointsToComm=NULL;
  static MPI_Request * SendPntNumberList=NULL;
  static MPI_Request * RecvPntNumberList=NULL;
  static double * SendCoordBuff = NULL;
  static MPI_Request * SendCoordList=NULL;
  static MPI_Request * RecvCoordList=NULL;
  static MPI_Request * SendDataBufferPtrList=NULL;
  static MPI_Request * RecvDataBufferPtrList=NULL;
  static char **  SendDataBuffPointerBuff =NULL;
  static double ** RecvCoordBuff = NULL;
  static char *** RecvDataBuffPointerBuff=NULL;

  if (globalMeshChangeFlag!=0) {
    //the mesh or the domain decomposition has been modified. Need to create a new communucation table
    bool meshModifiedFlag_CountMeshElements=PIC::Mesh::mesh->meshModifiedFlag_CountMeshElements;
    double StartTime=MPI_Wtime();

   
    if (nPointsToComm==NULL) nPointsToComm = new int [PIC::nTotalThreads];
    if (SendPntNumberList==NULL) SendPntNumberList = new MPI_Request[PIC::nTotalThreads-1]; 
    if (RecvPntNumberList==NULL) RecvPntNumberList = new MPI_Request[PIC::nTotalThreads-1];

    if (!SendCoordList) SendCoordList = new MPI_Request[PIC::nTotalThreads-1];
    if (!RecvCoordList) RecvCoordList = new MPI_Request[PIC::nTotalThreads-1];
    if (!SendDataBufferPtrList) SendDataBufferPtrList = new MPI_Request[PIC::nTotalThreads-1];
    if (!RecvDataBufferPtrList) RecvDataBufferPtrList = new MPI_Request[PIC::nTotalThreads-1]; 
 
    int  nPointsToCommThisThread = 0;
    //determine the new length of the table
    //    for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
      double eps = 0.3*PIC::Mesh::mesh->EPS;
      double dx[3];
      int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
      for (int idim=0; idim<3; idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];
        
      if (node->Thread==PIC::ThisThread && node->block) { 
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {

          double x[3];
          int ind[3]={i,j,k};
          bool outBoundary=false;// atBoundary=false;
          for (int idim=0; idim<3; idim++) {
            x[idim]=node->xmin[idim]+dx[idim]*(ind[idim]+0.5);

            #if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_
            if ((x[idim]-PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim]+0.5*dx[idim])<-eps || (x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim]-0.5*dx[idim])>eps){
              outBoundary=true;
              break;
            }
            #endif
          }

          if (outBoundary==false) nPointsToCommThisThread++;
        
        }// for (iface=0;iface<6;iface++)  for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)
      }//if (node->Thread==PIC::ThisThread)
    }//for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)

    nPointsToComm[PIC::ThisThread] = nPointsToCommThisThread;

    int iProcessCnt=0;
    //send the number of points on this proc to all other proc
    for (int iProc=0; iProc<PIC::nTotalThreads; iProc++){
      if (iProc!=PIC::ThisThread){
        MPI_Isend(nPointsToComm+PIC::ThisThread,1,MPI_INT,iProc,0,MPI_GLOBAL_COMMUNICATOR,SendPntNumberList+iProcessCnt);
        iProcessCnt++;
      }
    }
    
    iProcessCnt=0;
    //recv the number of points on other proc
    for (int iProc=0; iProc<PIC::nTotalThreads; iProc++){
      if (iProc!=PIC::ThisThread){
        MPI_Irecv(nPointsToComm+iProc,1,MPI_INT,iProc,0,MPI_GLOBAL_COMMUNICATOR,RecvPntNumberList+iProcessCnt);
        iProcessCnt++;
      }
    }
   
    //printf("thisthread:%d, nPointsToComm:%d\n",PIC::ThisThread,nPointsToComm[PIC::ThisThread]);

    if (SendCoordBuff!=NULL) delete [] SendCoordBuff;

    SendCoordBuff = new double [3*nPointsToComm[PIC::ThisThread]];

    if (RecvCoordBuff == NULL) {
      RecvCoordBuff = new double * [PIC::nTotalThreads-1];
      for (i=0; i<PIC::nTotalThreads-1; i++) RecvCoordBuff[i] = NULL;
    }

    if (SendDataBuffPointerBuff!=NULL) delete [] SendDataBuffPointerBuff;

    SendDataBuffPointerBuff =  new char * [nPointsToComm[PIC::ThisThread]];

    if  (RecvDataBuffPointerBuff == NULL) {
      RecvDataBuffPointerBuff= new char ** [PIC::nTotalThreads-1];
      for (i=0; i<PIC::nTotalThreads-1; i++) RecvDataBuffPointerBuff[i] = NULL;
    }

    double * tempPnt = SendCoordBuff;
    char ** tempDataBuffPnt = SendDataBuffPointerBuff;

    //    for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

      if (node->Thread==PIC::ThisThread && node->block) { 
        double dx[3];
        int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
        double eps = 0.3*PIC::Mesh::mesh->EPS;
        for (int idim=0; idim<3; idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];

        for (iface=0;iface<6;iface++)  for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {          
          double x[3];
          int ind[3]={i,j,k};
          bool outBoundary=false;
          // if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_){
          for (int idim=0; idim<3; idim++) {
            x[idim]=node->xmin[idim]+dx[idim]*(ind[idim]+0.5);

            #if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_
            if ((x[idim]-PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim]+0.5*dx[idim])<-eps || (x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim]-0.5*dx[idim])>eps) {
              outBoundary=true;
              break;
            }
            #endif
          }

          if (!outBoundary){
            for (int idim=0; idim<3; idim++) x[idim]=node->xmin[idim]+dx[idim]*(ind[idim]+0.5);

            for (int idim=0; idim<3; idim++)  {

              #if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_
              if (fabs(x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim]-0.5*dx[idim])<eps) x[idim] = PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim]-0.5*dx[idim];
              #endif
              tempPnt[idim]=x[idim];
            }

            *tempDataBuffPnt = node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();
            tempDataBuffPnt++;
            tempPnt += 3;
          }


        }// for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)
      }//if (node->Thread==PIC::ThisThread)
    }//for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)
    
    
    iProcessCnt=0;
    //send the number of points on this proc to all other proc
    for (int iProc=0; iProc<PIC::nTotalThreads; iProc++){
      if (iProc!=PIC::ThisThread){
        MPI_Isend(SendCoordBuff,nPointsToComm[PIC::ThisThread]*3,MPI_DOUBLE,iProc,1,MPI_GLOBAL_COMMUNICATOR,SendCoordList+iProcessCnt);
        MPI_Isend(SendDataBuffPointerBuff,nPointsToComm[PIC::ThisThread]*sizeof(char *),MPI_BYTE,iProc,2,MPI_GLOBAL_COMMUNICATOR,SendDataBufferPtrList+iProcessCnt);
        iProcessCnt++;
      }
    }
    
    int nRecvBuffDone=0;      
                      
    while (nRecvBuffDone<PIC::nTotalThreads-1){
      int q, flag, iProc;
      MPI_Testany(PIC::nTotalThreads-1, RecvPntNumberList, &q, &flag, MPI_STATUS_IGNORE);

      if ((flag!=0)&&(q!=MPI_UNDEFINED)) {
        //release memoty used by MPI
        MPI_Wait(RecvPntNumberList+q,MPI_STATUS_IGNORE);

        iProc = q<PIC::ThisThread?q:q+1;
        if (RecvCoordBuff[q]!=NULL) delete [] RecvCoordBuff[q];
        RecvCoordBuff[q] = new double [3*nPointsToComm[iProc]];
        if (RecvDataBuffPointerBuff[q]!=NULL) delete [] RecvDataBuffPointerBuff[q];
        RecvDataBuffPointerBuff[q] = new char * [nPointsToComm[iProc]];
        MPI_Irecv(RecvCoordBuff[q],nPointsToComm[iProc]*3,MPI_DOUBLE,iProc,1,MPI_GLOBAL_COMMUNICATOR,RecvCoordList+q);
        MPI_Irecv(RecvDataBuffPointerBuff[q],nPointsToComm[iProc]*sizeof(char *),MPI_BYTE,iProc,2,MPI_GLOBAL_COMMUNICATOR,RecvDataBufferPtrList+q);
        
        nRecvBuffDone++;
      }
    }     

    PIC::Parallel::XYZTree Tree;

    PIC::Parallel::XYZTree::leafNode ** leafNodeArr = new  PIC::Parallel::XYZTree::leafNode * [PIC::nTotalThreads];
   
    for (iThread=0; iThread<PIC::nTotalThreads; iThread++) {
      // printf("nPointsToComm[%d]:%d\n",iThread,nPointsToComm[iThread]);
      leafNodeArr[iThread] = new PIC::Parallel::XYZTree::leafNode [nPointsToComm[iThread]];
    }
    //create tree from local thread
    for (int iPnt=0; iPnt<nPointsToComm[PIC::ThisThread]; iPnt++){
      leafNodeArr[PIC::ThisThread][iPnt].iThread = PIC::ThisThread;
      leafNodeArr[PIC::ThisThread][iPnt].DataBuffer = SendDataBuffPointerBuff[iPnt];
    }
    
    //build tree from send buffer
    for (int iPnt=0; iPnt<nPointsToComm[PIC::ThisThread]; iPnt++){
      Tree.addXNode(SendCoordBuff[3*iPnt],SendCoordBuff[3*iPnt+1],SendCoordBuff[3*iPnt+2],leafNodeArr[PIC::ThisThread]+iPnt);
    }

    int counterArr[PIC::nTotalThreads-1];
    for (i=0; i<PIC::nTotalThreads-1; i++) counterArr[i]=0;
    int counter[2]={0,0};
    //process other threads
        
    int nProcDone=0;      
                      
    //build tree from recv buffers
    while (nProcDone<PIC::nTotalThreads-1){
      int p[2], flag[2];
      
      if (counter[0]<PIC::nTotalThreads-1) {
        MPI_Testany(PIC::nTotalThreads-1, RecvCoordList, p, flag, MPI_STATUS_IGNORE);

        //release memory used by MPI
        if ((flag[0]==true)&&(p[0]!=MPI_UNDEFINED)) MPI_Wait(RecvCoordList+p[0],MPI_STATUS_IGNORE);
      }


      if (counter[1]<PIC::nTotalThreads-1) {
        MPI_Testany(PIC::nTotalThreads-1, RecvDataBufferPtrList, p+1, flag+1, MPI_STATUS_IGNORE);

        //release memory used by MPI
        if ((flag[1]==true)&&(p[1]!=MPI_UNDEFINED)) MPI_Wait(RecvDataBufferPtrList+p[1],MPI_STATUS_IGNORE);
      }
      
      for (i=0;i<2; i++){
        if (flag[i]!=0 && counter[i]<PIC::nTotalThreads-1){
          int localInd = p[i];
          counterArr[localInd]++;
          counter[i]++;
          if (counterArr[localInd]==2){
            int iProc=localInd<PIC::ThisThread?localInd:localInd+1;
           
            for (int iPnt=0; iPnt<nPointsToComm[iProc]; iPnt++){
              leafNodeArr[iProc][iPnt].iThread = iProc;
              leafNodeArr[iProc][iPnt].DataBuffer = RecvDataBuffPointerBuff[localInd][iPnt];
            }
            
            for (int iPnt=0; iPnt<nPointsToComm[iProc]; iPnt++){
              Tree.addXNode(RecvCoordBuff[localInd][3*iPnt],RecvCoordBuff[localInd][3*iPnt+1],RecvCoordBuff[localInd][3*iPnt+2],leafNodeArr[iProc]+iPnt);
            }
                        
            nProcDone++;
          }          
        }
      }
    }//while (nProcDone<PIC::nTotalThreads-1)
    
    // printf("nProcDone:%d\n",nProcDone);
    // printf("\n");
    // printf("thisthread:%d ,counterarr: ",PIC::ThisThread);
    // for (i=0;i<PIC::nTotalThreads-1;i++) printf(" %d",counterArr[i]);
    //printf("\n");
    //if (PIC::ThisThread==0) Tree.printTree();

    int totalStencil=0;
    std::list<PIC::Parallel::XYZTree::xNode*>::iterator itX;

    for (itX=Tree.xList.begin(); itX!=Tree.xList.end(); itX++){
      std::list<PIC::Parallel::XYZTree::yNode*>::iterator itY;
      for (itY=(*itX)->yList.begin(); itY!=(*itX)->yList.end(); itY++){
        //totalStencil += (*itY)->zList.size();
        std::list<PIC::Parallel::XYZTree::zNode*>::iterator itZ;
        for (itZ=(*itY)->zList.begin(); itZ!=(*itY)->zList.end(); itZ++){
          std::list<PIC::Parallel::XYZTree::leafNode*>::iterator itLeaf;
          if ((*itZ)->leafList.size()>1) totalStencil +=1;
        }
      }
    }

    //allocate the new Stencile Table
    if (StencilTableLength!=0) delete [] StencilTable;

    StencilTableLength=totalStencil;
    StencilTable=new cStencilElement[totalStencil];
    //printf("totalStencil:%d, thisthread:%d\n", totalStencil, PIC::ThisThread);

    int totalNodes=0;
    int iSt =0;

    for (itX=Tree.xList.begin(); itX!=Tree.xList.end(); itX++){
      std::list<PIC::Parallel::XYZTree::yNode*>::iterator itY;

      for (itY=(*itX)->yList.begin(); itY!=(*itX)->yList.end(); itY++){
        std::list<PIC::Parallel::XYZTree::zNode*>::iterator itZ;

        for (itZ=(*itY)->zList.begin(); itZ!=(*itY)->zList.end(); itZ++){
          std::list<PIC::Parallel::XYZTree::leafNode*>::iterator itLeaf;
          int iTh=0;
          if ((*itZ)->leafList.size()<2) continue;
          /*
          if (fabs((*itX)->x)<1e-3 && fabs((*itY)->y)<1e-3 && fabs((*itZ)->z)<1e-3)
            printf("test ist:%d,x:%e,y:%e,z:%e,StencilLength:%d\n",iSt,(*itX)->x,(*itY)->y,(*itZ)->z,StencilTable[iSt].StencilLength);
          */
          int tempThreadId=-1;
          StencilTable[iSt].StencilLength=0;
          StencilTable[iSt].AssociatedDataPointer=NULL;
          for (itLeaf=(*itZ)->leafList.begin(); itLeaf!=(*itZ)->leafList.end(); itLeaf++){            
            /*
            if (fabs((*itX)->x)<1e-3 && fabs((*itY)->y)<1e-3 && fabs((*itZ)->z)<1e-3){
              printf("test databuffer thisthread:%d, 0,0,0,iThread:%d, DataBuffer:%p\n", PIC::ThisThread,(*itLeaf)->iThread,(void*)((*itLeaf)->DataBuffer));
            }
            */            
            if (tempThreadId!=(*itLeaf)->iThread) {
              StencilTable[iSt].StencilLength++;
              StencilTable[iSt].StencilThreadTable[iTh]=(*itLeaf)->iThread;
              tempThreadId = (*itLeaf)->iThread;
              iTh++;
            }

            if ((*itLeaf)->iThread==PIC::ThisThread) {
              if (StencilTable[iSt].AssociatedDataPointer==NULL) {
                StencilTable[iSt].AssociatedDataPointer=(*itLeaf)->DataBuffer;
                if ((*itLeaf)->DataBuffer==NULL) exit(__LINE__,__FILE__,"Error: NULL pointer");
              }
              else{
                //more than 1 data buffer for one point
                StencilTable[iSt].localDataPntArr.push_back((*itLeaf)->DataBuffer);
                if ((*itLeaf)->DataBuffer==NULL) exit(__LINE__,__FILE__,"Error: NULL pointer");
              }
            }
        
          }
          //determine the center-processing thread
          int t = iSt%StencilTable[iSt].StencilLength;
          int temp = StencilTable[iSt].StencilThreadTable[0];
          StencilTable[iSt].StencilThreadTable[0]=StencilTable[iSt].StencilThreadTable[t];
          StencilTable[iSt].StencilThreadTable[t]=temp;         
          totalNodes+=StencilTable[iSt].StencilLength;
          /*
          if (PIC::ThisThread==0){
            printf("iSt:%d, stencilLength:%d, threadtable:",iSt, StencilTable[iSt].StencilLength);
            for (int tt=0; tt<StencilTable[iSt].StencilLength; tt++) printf("%d ",StencilTable[iSt].StencilThreadTable[tt]);
            printf("x,y,z:%e,%e,%e\n",(*itX)->x,(*itY)->y,(*itZ)->z);
            printf("\n");
          }
          */
          iSt++;
        }
      }
    }

    //test

    /*
    char * testPtr = StencilTable[714].AssociatedDataPointer;
     for (itX=Tree.xList.begin(); itX!=Tree.xList.end(); itX++){
      std::list<PIC::Parallel::XYZTree::yNode*>::iterator itY;
      for (itY=(*itX)->yList.begin(); itY!=(*itX)->yList.end(); itY++){
        std::list<PIC::Parallel::XYZTree::zNode*>::iterator itZ;
        for (itZ=(*itY)->zList.begin(); itZ!=(*itY)->zList.end(); itZ++){
          std::list<PIC::Parallel::XYZTree::leafNode*>::iterator itLeaf;
          
          if ((*itZ)->leafList.size()<2) continue;
   
          for (itLeaf=(*itZ)->leafList.begin(); itLeaf!=(*itZ)->leafList.end(); itLeaf++){
            if (fabs((*itX)->x)>2+1e-3 || fabs((*itY)->y)>2+1e-3 || fabs((*itZ)->z)>2+1e-3 ) printf("points out of range!\n");
            if (fabs((*itX)->x-2)<1e-3 || fabs((*itY)->y-2)<1e-3 || fabs((*itZ)->z-2)<1e-3 ) printf("wrong points !\n");
            if (testPtr==(*itLeaf)->DataBuffer){
              printf("repeated:x:%e,y:%e,z:%e, thread id:%d, Jx:%e\n", (*itX)->x,(*itY)->y, (*itZ)->z, PIC::ThisThread,((double *)(StencilTable[714].AssociatedDataPointer+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))[6] );
            }
          }

        }
     
      }
    
     }//itX
     */

    // Tree.printTree();
    //printf("thisthread:%d,iSt:%d,totalNodes:%d\n",PIC::ThisThread,iSt,totalNodes);  
    //update the coundater
    nMeshModificationCounter=PIC::Mesh::mesh->nMeshModificationCounter;
    PIC::Mesh::mesh->meshModifiedFlag_CountMeshElements=meshModifiedFlag_CountMeshElements;
    PIC::Parallel::RebalancingTime+=MPI_Wtime()-StartTime;
  }


  
  static int * LoadList=NULL;
  static std::vector<int> ProcessStencilBufferNumber;
  static std::vector<char*> CopyDataBufferDest;
  static std::vector<int> CopyDataBufferStencilIndex;
  static std::vector<int> WholeIndexOfProcessStencil;
  static char ** processRecvDataBuffer;
  static char ** copyRecvDataBuffer;
  static int * wholeStencilIndexOfBuffer=NULL,* processStencilIndexOfBuffer=NULL;
  static int * StencilFinishedBufferNumber=NULL;
  static int sumProcessBufferNumber=0;
  static int totalProcessStencil;
  static int totalCopyStencil;
 
  static MPI_Request * processRecvList=NULL, * copyRecvList=NULL;
  static MPI_Request * processSendList=NULL, * copySendList=NULL;

  if (globalMeshChangeFlag!=0 && PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData!=NULL){

    if (LoadList!=NULL) delete [] LoadList;
    LoadList=new int [PIC::Mesh::mesh->nTotalThreads];
    for (int i=0;i<PIC::Mesh::mesh->nTotalThreads;i++) LoadList[i]=0;
    
    sumProcessBufferNumber = 0;
    //print out the thread distribution of process stencils 
    for (int iStencil=0;iStencil<StencilTableLength;iStencil++) {
      //if (StencilTable[iStencil].StencilLength>1) {
      int center=StencilTable[iStencil].StencilThreadTable[0];
      if (StencilTable[iStencil].StencilLength>1) LoadList[center]++;
      if (center==PIC::ThisThread) sumProcessBufferNumber += StencilTable[iStencil].StencilLength-1;
      //}
    }
    
    int sumLoad=0;
    for (int i=0;i<PIC::Mesh::mesh->nTotalThreads;i++) {
      //  printf("thread:%d, loadlist[%d],%d",PIC::ThisThread,i,LoadList[i]);
      sumLoad+=LoadList[i];
    }
    
    //printf("thread:%d, StencilTableLength:%d\n",PIC::ThisThread,StencilTableLength);
    //printf("\n");
    //    if (sumLoad!=StencilTableLength) exit(__LINE__,__FILE__,"Error: Stencils containing only one thread/proc are in the StencilTable");
        
    if (wholeStencilIndexOfBuffer!=NULL) delete [] wholeStencilIndexOfBuffer;
    if (processStencilIndexOfBuffer!=NULL) delete [] processStencilIndexOfBuffer;
    wholeStencilIndexOfBuffer=new int[sumProcessBufferNumber];
    processStencilIndexOfBuffer =new int[sumProcessBufferNumber];

    int processBufferIndex=0;

    ProcessStencilBufferNumber.clear();
    CopyDataBufferDest.clear();
    CopyDataBufferStencilIndex.clear();
    WholeIndexOfProcessStencil.clear();
    for (int iStencil=0;iStencil<StencilTableLength;iStencil++) {        
      if (StencilTable[iStencil].AssociatedDataPointer){// means this thread is involved
        if (StencilTable[iStencil].StencilThreadTable[0]==PIC::ThisThread) {
          if (StencilTable[iStencil].StencilLength>1){
            
            ProcessStencilBufferNumber.push_back(StencilTable[iStencil].StencilLength-1);
            WholeIndexOfProcessStencil.push_back(iStencil);
          }
          for (int iThread=0; iThread<StencilTable[iStencil].StencilLength-1;iThread++){
            //process buffer number equals StencilLength-1
            wholeStencilIndexOfBuffer[processBufferIndex]=iStencil; //index in the whole stensil list
            processStencilIndexOfBuffer[processBufferIndex]=ProcessStencilBufferNumber.size()-1; 
            processBufferIndex++;
          }
        }else{ 
          CopyDataBufferDest.push_back(StencilTable[iStencil].AssociatedDataPointer);
          CopyDataBufferStencilIndex.push_back(iStencil);
        }//else    
      }
    }
    //printf("thread id:%d, processBufferIndex:%d\n", PIC::ThisThread, processBufferIndex);

    if (processRecvDataBuffer!=NULL) {
      delete [] processRecvDataBuffer[0]; 
      delete [] processRecvDataBuffer;
    }
    
    processRecvDataBuffer= new char * [sumProcessBufferNumber];
    processRecvDataBuffer[0] = new char [sumProcessBufferNumber*PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength];
    for (int i=1;i<sumProcessBufferNumber;i++) {
      processRecvDataBuffer[i]=i*PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+processRecvDataBuffer[0];
    }

    if (PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData!=NULL) {

      if (copyRecvDataBuffer!=NULL){
        delete [] copyRecvDataBuffer[0];
        delete [] copyRecvDataBuffer;
      }
      
      copyRecvDataBuffer=new char * [CopyDataBufferDest.size()];
      copyRecvDataBuffer[0]=new char [CopyDataBufferDest.size()*PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength];
      for (int i=1;i<CopyDataBufferDest.size();i++) copyRecvDataBuffer[i]=copyRecvDataBuffer[0]+i*PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength;
    }
    
    
    if (processSendList!=NULL) delete[] processSendList;
    processSendList = new MPI_Request[CopyDataBufferDest.size()];
    if (processRecvList!=NULL) delete[] processRecvList;
    processRecvList = new MPI_Request[sumProcessBufferNumber];
    if (copyRecvList!=NULL) delete [] copyRecvList; 
    copyRecvList = new MPI_Request[CopyDataBufferDest.size()];
    if (copySendList!=NULL) delete [] copySendList; 
    copySendList = new MPI_Request[sumProcessBufferNumber];
   
    //printf("sumProcessBufferNumber:%d, copyBufferSize:%d\n",sumProcessBufferNumber,CopyDataBufferDest.size());
    
    if (StencilFinishedBufferNumber!=NULL) delete [] StencilFinishedBufferNumber;
    StencilFinishedBufferNumber=new int [LoadList[PIC::ThisThread]];
    totalProcessStencil= LoadList[PIC::ThisThread];
    totalCopyStencil = CopyDataBufferDest.size();
  
  }//if (globalMeshChangeFlag!=0 && PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL)

  //printf("thread id:%d, sumProcessBufferNumber:%d,totalProcessStencil:%d,WholeIndexOfProcessStencil size:%d\n",PIC::ThisThread,sumProcessBufferNumber,totalProcessStencil,WholeIndexOfProcessStencil.size());
  
  if (PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData!=NULL) {

    //clear counter
    for (int i=0; i<totalProcessStencil;i++) StencilFinishedBufferNumber[i]=0;

    //1. combine 'corner' node data                                                                                                            
    int iProcessBuffer=0, nProcessBufferDone=0, nProcessStencilDone=0;
    int iCopyBuffer=0, nCopyBufferDone=0;
    int iProcessSendBuffer=0;



    for (int iStencil=0;iStencil<StencilTableLength;iStencil++) {//loop through stencil table
      if (StencilTable[iStencil].AssociatedDataPointer) { // this thread is involved
        int iThread;
        //there are more that one MPI processes that contributed to the state vector of the corner node
      
        if (PIC::ThisThread==StencilTable[iStencil].StencilThreadTable[0]) {
          //if this thread will do the center processing role and  combine data from all other MPI processes
          //recv data request start
          //process local points
          std::vector<char *>::iterator it;
          for (it=StencilTable[iStencil].localDataPntArr.begin();it!=StencilTable[iStencil].localDataPntArr.end();it++){
            PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,*it);
          }
        
          for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
            MPI_Irecv(processRecvDataBuffer[iProcessBuffer],PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[iThread],iStencil,MPI_GLOBAL_COMMUNICATOR,processRecvList+iProcessBuffer);
            iProcessBuffer++;
          }
          
        }else{
          //this thread will be copyThread, send data to the processThread, 
          //and recv the processed data from the processThread
          std::vector<char *>::iterator it;
          for (it=StencilTable[iStencil].localDataPntArr.begin();it!=StencilTable[iStencil].localDataPntArr.end();it++){
            PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,*it);                  
          }
          MPI_Isend(StencilTable[iStencil].AssociatedDataPointer,PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,processSendList+iCopyBuffer);
            
          //start the recv data request
          if (PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData!=NULL) {
            MPI_Irecv(copyRecvDataBuffer[iCopyBuffer],PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,copyRecvList+iCopyBuffer);
          }else{
            MPI_Irecv(CopyDataBufferDest[iCopyBuffer],PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,copyRecvList+iCopyBuffer);
          }
          iCopyBuffer++;
          // }//for
        }//else

        
        for (int i=nProcessBufferDone; i<iProcessBuffer;i++) {
          int q, flag;
          // printf("test7 pass\n");
          //if (nProcessBufferDone<sumProcessBufferNumber)
          MPI_Testany(iProcessBuffer, processRecvList, &q, &flag, MPI_STATUS_IGNORE);
          // printf("test8 pass\n");
          //printf("nProcessBufferDone:%d, sumProcessBufferNumber:%d\n",nProcessBufferDone,sumProcessBufferNumber);

          if (flag!=0 && q!=MPI_UNDEFINED){
            //process buffer of index q is done 
            int processStencilIndex = processStencilIndexOfBuffer[q];
            int wholeStencilIndex = wholeStencilIndexOfBuffer[q];
            
            //release memoty used by MPI
            MPI_Wait(processRecvList+q,MPI_STATUS_IGNORE);

            //printf("thread id:%d,bufferid:%d processed, nProcessBufferDone:%d/total buffer:%d,\n",PIC::ThisThread,q,nProcessBufferDone,sumProcessBufferNumber);

            PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData(StencilTable[wholeStencilIndex].AssociatedDataPointer,processRecvDataBuffer[q]);
          
            nProcessBufferDone++;
            StencilFinishedBufferNumber[processStencilIndex]++;
            //if  finished one stencil processing, send it back
            if (StencilFinishedBufferNumber[processStencilIndex]==ProcessStencilBufferNumber[processStencilIndex]){
              for (iThread=1;iThread<StencilTable[wholeStencilIndex].StencilLength;iThread++) {        
                MPI_Isend(StencilTable[wholeStencilIndex].AssociatedDataPointer,PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[wholeStencilIndex].StencilThreadTable[iThread],wholeStencilIndex,MPI_GLOBAL_COMMUNICATOR,copySendList+iProcessSendBuffer);
                iProcessSendBuffer++;
              }
              nProcessStencilDone++;
            }
            
          }
        }//for (int i=nProcessBufferDone; i<iProcessBuffer;i++)

        for (int i=nCopyBufferDone; i<iCopyBuffer;i++) {
          int q, flag;

          MPI_Testany(iCopyBuffer, copyRecvList, &q, &flag, MPI_STATUS_IGNORE);
          
          if (flag!=0 && q!=MPI_UNDEFINED) {
            //release memoty used by MPI
            MPI_Wait(copyRecvList+q,MPI_STATUS_IGNORE);

            //int wholeStencilIndex=CopyDataBufferStencilIndex[q];
            if (PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData!=NULL) {
              PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData(CopyDataBufferDest[q],copyRecvDataBuffer[q]);
            }
            //if CopyCornerNodeAssociatedData is not defined, the MPI_Irecv already sent the data to the destination state vector 
            nCopyBufferDone++;
          }
        }//for (int i=nCopyBufferDone; i<iCopyBuffer;i++)
      }
    }//for (iStencil=0;iStencil<StencilTableLength;iStencil++)
    
    
    //printf("out of for-loop\n");
    
    //make sure communication is finished
    while (nProcessStencilDone<totalProcessStencil||nCopyBufferDone<totalCopyStencil){
      
      for (int i=nProcessBufferDone; i<iProcessBuffer;i++) {
        int q, flag;
            
        MPI_Testany(iProcessBuffer, processRecvList, &q, &flag, MPI_STATUS_IGNORE);
        
        if (flag!=0 &&  q!=MPI_UNDEFINED){
          //process buffer of index q is done 
          int processStencilIndex = processStencilIndexOfBuffer[q];
          int wholeStencilIndex = wholeStencilIndexOfBuffer[q];

          //release memoty used by MPI
          MPI_Wait(processRecvList+q,MPI_STATUS_IGNORE);
            
          // printf("thread id:%d,bufferid:%d processed, nProcessBufferDone:%d/total buffer:%d,\n",PIC::ThisThread,q,nProcessBufferDone,sumProcessBufferNumber);
          /*
            if (wholeStencilIndex==714){
            printf("ist:%d,in processing thread id: %d, size:%d, (while) before Jx:%e\n",wholeStencilIndex,PIC::ThisThread,StencilTable[wholeStencilIndex].localDataPntArr.size(), ((double *)(StencilTable[wholeStencilIndex].AssociatedDataPointer+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))[6]);
            }
          */
          PIC::Parallel::CenterBlockBoundaryNodes::ProcessCenterNodeAssociatedData(StencilTable[wholeStencilIndex].AssociatedDataPointer,processRecvDataBuffer[q]);
          
          nProcessBufferDone++;
          StencilFinishedBufferNumber[processStencilIndex]++;
                
          if (StencilFinishedBufferNumber[processStencilIndex]==ProcessStencilBufferNumber[processStencilIndex]){
            for (iThread=1;iThread<StencilTable[wholeStencilIndex].StencilLength;iThread++) {        
              MPI_Isend(StencilTable[wholeStencilIndex].AssociatedDataPointer,PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength,MPI_BYTE,StencilTable[wholeStencilIndex].StencilThreadTable[iThread],wholeStencilIndex,MPI_GLOBAL_COMMUNICATOR,copySendList+iProcessSendBuffer);
              iProcessSendBuffer++;
            }
            nProcessStencilDone++;
          }
        }
      }//for (int i=nProcessBufferDone; i<iProcessBuffer;i++)
      
      for (int i=nCopyBufferDone; i<iCopyBuffer;i++) {
        int q, flag;
        
        MPI_Testany(iCopyBuffer, copyRecvList, &q, &flag, MPI_STATUS_IGNORE);
        
        if (flag!=0 && q!=MPI_UNDEFINED) {
          //release memoty used by MPI
          MPI_Wait(copyRecvList+q,MPI_STATUS_IGNORE);

          if (PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData!=NULL) {
            PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData(CopyDataBufferDest[q],copyRecvDataBuffer[q]);
          }
          //if CopyCornerNodeAssociatedData is not defined, the MPI_Irecv already sent the data to the destination state vector 
          nCopyBufferDone++;
        }
      }//for (int i=nCopyBufferDone; i<iCopyBuffer;i++)
    }//while (nProcessStencilDone<totalProcessStencil||nCopyBufferDone<totalCopyStencil)
    
    //sync the local data buffers 
    for (int iSt=0;iSt<StencilTableLength;iSt++) {
      if (StencilTable[iSt].AssociatedDataPointer) {
        
        std::vector<char *>::iterator it;
        if (PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData!=NULL) {
          for (it=StencilTable[iSt].localDataPntArr.begin();it!=StencilTable[iSt].localDataPntArr.end();it++){
            PIC::Parallel::CenterBlockBoundaryNodes::CopyCenterNodeAssociatedData(*it,StencilTable[iSt].AssociatedDataPointer);  
          }
        }else{
          for (it=StencilTable[iSt].localDataPntArr.begin();it!=StencilTable[iSt].localDataPntArr.end();it++)
            memcpy(*it,StencilTable[iSt].AssociatedDataPointer,PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength);
        }
      }
    }
   
    MPI_Waitall(iProcessBuffer,processRecvList,MPI_STATUSES_IGNORE);
    MPI_Waitall(iCopyBuffer,copyRecvList,MPI_STATUSES_IGNORE);
    MPI_Waitall(iCopyBuffer,processSendList,MPI_STATUSES_IGNORE);
    MPI_Waitall(iProcessBuffer,copySendList,MPI_STATUSES_IGNORE);
    
  }//if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL)
 
}

//------------------------------------------   ProcessCornerBlockBoundaryNodes: begin  ----------------------------------------
void PIC::Parallel::ProcessCornerBlockBoundaryNodes() {
  const int iFaceMin[6]={0,_BLOCK_CELLS_X_,0,0,0,0};
  const int iFaceMax[6]={0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_};
  
  const int jFaceMin[6]={0,0,0,_BLOCK_CELLS_Y_,0,0};
  const int jFaceMax[6]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_};
  
  const int kFaceMin[6]={0,0,0,0,0,_BLOCK_CELLS_Z_};
  const int kFaceMax[6]={_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,0,_BLOCK_CELLS_Z_};
  
  
  int periodic_bc_pair_real_block=-1;
  int periodic_bc_pair_ghost_block=-1;
  
  auto ReleasePeriodicBCFlags = [&] () {
    PIC::Mesh::mesh->rootTree->ReleaseFlag(periodic_bc_pair_real_block);
    PIC::Mesh::mesh->rootTree->ReleaseFlag(periodic_bc_pair_ghost_block);
  };
  
  
  //-----------------------------------   ACCUMULATION THE REQUESTS: begin  ------------------------------
  class cRequestData {
  public:
    int iSource,jSource,kSource;
    int iTarget,jTarget,kTarget;
    int iSourceThread,iTargetThread;
    cAMRnodeID SourceNodeId,TargetNodeId;
    char *Buffer;
    char *CornerAssociatedData;
    
    cRequestData() {
      iSource=-1,jSource=-1,kSource=-1;
      iTarget=-1,jTarget=-1,kTarget=-1;
      iSourceThread=-1,iTargetThread=-1;
      Buffer=NULL,CornerAssociatedData=NULL;
    }
  };
  
  list <cRequestData> RequestDataList;
  
  //create the request data list 
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,list <cRequestData> *)> BuildRequestDataList;
  
  auto ProcessCornerPoint = [&] (double *xCorner,double *dx,cRequestData* RequestTable,int &RequestTableLength,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* Node, int ExcludeThread) {
    int iSearch,jSearch,kSearch;
    int iCell,jCell,kCell;
    double xSearch[3];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;
    cRequestData *Request;

    RequestTableLength=0;

    for (iSearch=-1;iSearch<=1;iSearch+=2) {
      xSearch[0]=xCorner[0]+iSearch*dx[0];

      for (jSearch=-1;jSearch<=1;jSearch+=2) {
        xSearch[1]=xCorner[1]+jSearch*dx[1];

        for (kSearch=-1;kSearch<=1;kSearch+=2) {
          xSearch[2]=xCorner[2]+kSearch*dx[2];

          neibNode=PIC::Mesh::mesh->findTreeNode(xSearch,Node);

          if (neibNode!=NULL) if (neibNode->TestFlag(periodic_bc_pair_ghost_block)==false) {
            //a MPI process can contribute only ones 
            bool flag=true;

            if (neibNode->Thread==ExcludeThread) {
              flag=false;
            }
            else {
              for (int ii=0;ii<RequestTableLength;ii++) if (neibNode->Thread==RequestTable[ii].iSourceThread) {
                flag=false;
                break;
              }
            }

            if (flag==false) continue;

            PIC::Mesh::mesh->fingCellIndex(xSearch,iCell,jCell,kCell,neibNode,false);

            Request=RequestTable+RequestTableLength;
            RequestTableLength++;

            Request->iSourceThread=neibNode->Thread;     
            Request->SourceNodeId=neibNode->AMRnodeID;

            Request->iSource=(iSearch==-1) ? iCell+1 : iCell;
            Request->jSource=(jSearch==-1) ? jCell+1 : jCell;
            Request->kSource=(kSearch==-1) ? kCell+1 : kCell;
          }
        }
      }
    }
  };
  
  
  //lambda to create the list 
  BuildRequestDataList = [&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* Node,list <cRequestData>* RequestDataList) -> void {
    int i,j,k,iNeib,jNeib,kNeib;
    PIC::Mesh::cDataBlockAMR *block;
    PIC::Mesh::cDataCornerNode *CornerNode;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;
    
    
    if (Node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      //the block at the botton of the graph
      
      //verify that the node is not a ghost
      if (Node->Thread!=PIC::ThisThread) return;
      if (Node->TestFlag(periodic_bc_pair_ghost_block)==true) return;
      if ((block=Node->block)==NULL) return;
      
      double dx[3],xSearch[3],xCorner[3],L[3];
      int iface,ii,iSearch,jSearch,kSearch,RequestTableLength;
      cRequestData LocalRequestTable[8];
      
      dx[0]=0.5*(Node->xmax[0]-Node->xmin[0])/_BLOCK_CELLS_X_;
      dx[1]=0.5*(Node->xmax[1]-Node->xmin[1])/_BLOCK_CELLS_Y_;
      dx[2]=0.5*(Node->xmax[2]-Node->xmin[2])/_BLOCK_CELLS_Z_;
      
      if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_ ) {
        for (int idim=0;idim<3;idim++) L[idim]=PIC::BC::ExternalBoundary::Periodic::L[idim]; //PIC::BC::ExternalBoundary::Periodic::xmaxDomain[idim]-PIC::BC::ExternalBoundary::Periodic::xminDomain[idim];
      }
      
      //loop through all corners
      for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if ((CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) if (CornerNode->TestProcessedFlag()==false) {
          CornerNode->SetProcessedFlag(true);
          CornerNode->GetX(xCorner);
          
          //create the request table for the corner node
          ProcessCornerPoint(xCorner,dx,LocalRequestTable,RequestTableLength,Node,PIC::ThisThread);
          
          for (ii=0;ii<RequestTableLength;ii++) {
            LocalRequestTable[ii].iTarget=i;
            LocalRequestTable[ii].jTarget=j;
            LocalRequestTable[ii].kTarget=k;
            LocalRequestTable[ii].iTargetThread=PIC::ThisThread;
            LocalRequestTable[ii].TargetNodeId=Node->AMRnodeID;
            LocalRequestTable[ii].CornerAssociatedData=CornerNode->GetAssociatedDataBufferPointer();
            
            RequestDataList->push_back(LocalRequestTable[ii]);
          }
          
          //in the case periodic BC are used -> create data requests that include data from 'ghost' blocks
          if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_ ) {
            double xPair[3];
            int iPair,jPair,kPair;
            
            for (iPair=-1;iPair<=1;iPair++) {
              xPair[0]=xCorner[0]+iPair*L[0];
              
              if ((xPair[0]>=PIC::Mesh::mesh->xGlobalMin[0])&&(xPair[0]<=PIC::Mesh::mesh->xGlobalMax[0])) {
                for (jPair=-1;jPair<=1;jPair++) {
                  xPair[1]=xCorner[1]+jPair*L[1];
                   
                  if ((xPair[1]>=PIC::Mesh::mesh->xGlobalMin[1])&&(xPair[1]<=PIC::Mesh::mesh->xGlobalMax[1])) {
                    for (kPair=-1;kPair<=1;kPair++) if ((iPair!=0)||(jPair!=0)||(kPair!=0)) { 
                      xPair[2]=xCorner[2]+kPair*L[2];
                      
                      if ((xPair[2]>=PIC::Mesh::mesh->xGlobalMin[2])&&(xPair[2]<=PIC::Mesh::mesh->xGlobalMax[2])) {
                        //create the request table for the corner node
                        ProcessCornerPoint(xPair,dx,LocalRequestTable,RequestTableLength,Node,-1);
                        
                        for (ii=0;ii<RequestTableLength;ii++) {
                          LocalRequestTable[ii].iTarget=i;
                          LocalRequestTable[ii].jTarget=j;
                          LocalRequestTable[ii].kTarget=k;
                          LocalRequestTable[ii].iTargetThread=PIC::ThisThread;
                          LocalRequestTable[ii].TargetNodeId=Node->AMRnodeID;
                          LocalRequestTable[ii].CornerAssociatedData=CornerNode->GetAssociatedDataBufferPointer();
                          
                          RequestDataList->push_back(LocalRequestTable[ii]);
                        }             
                      }
                    }
                  }
                }
              }
            }          
          }
        }
      }
    }
    else {
      int iDownNode;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;
      
      for (iDownNode=0;iDownNode<(1<<DIM);iDownNode++) if ((downNode=Node->downNode[iDownNode])!=NULL) {
        BuildRequestDataList(downNode,RequestDataList);
      }
    }
  };
  
  //-----------------------------------   ACCUMULATION THE REQUESTS: end  --------------------------------
  
  
  //-----------------------------------   EXCHANGE THE REQUESTS: begin  --------------------------------
  static int *GlobalRequestNumberTable=NULL;
  static cRequestData **SendRequestTable=NULL;
  static cRequestData **RecvRequestTable=NULL;
  
  auto GetRecvRequestNumber = [&] (int iSource) {return GlobalRequestNumberTable[iSource+PIC::ThisThread*PIC::nTotalThreads];};
  auto GetSendRequestNumber = [&] (int iTarget) {return GlobalRequestNumberTable[PIC::ThisThread+iTarget*PIC::nTotalThreads];};
 
  auto CreateGlobalRequestTable = [&] () {
    
    if (GlobalRequestNumberTable==NULL) {
      GlobalRequestNumberTable=new int [PIC::nTotalThreads*PIC::nTotalThreads];
      
      SendRequestTable=new cRequestData* [PIC::nTotalThreads];
      RecvRequestTable=new cRequestData* [PIC::nTotalThreads];
    }
    else {
      for (int thread=0;thread<PIC::nTotalThreads;thread++) {
        if (SendRequestTable[thread]!=NULL) delete [] SendRequestTable[thread];
        if (RecvRequestTable[thread]!=NULL) delete [] RecvRequestTable[thread];
      }
    }
    
    for (int ii=0;ii<PIC::nTotalThreads*PIC::nTotalThreads;ii++) GlobalRequestNumberTable[ii]=0;
    for (int ii=0;ii<PIC::nTotalThreads;ii++) SendRequestTable[ii]=NULL,RecvRequestTable[ii]=NULL;
    
    //create the Recv request table
    int *LocalRecvRequestTableSize=new int [PIC::nTotalThreads];
    int *LocalSendRequestTableSize=new int [PIC::nTotalThreads];
    
    for (int ii=0;ii<PIC::nTotalThreads;ii++) LocalRecvRequestTableSize[ii]=0;
    
    for (auto& it : RequestDataList) LocalRecvRequestTableSize[it.iSourceThread]++;
    
    for (int ii=0;ii<PIC::nTotalThreads;ii++) {
      RecvRequestTable[ii]=new cRequestData [LocalRecvRequestTableSize[ii]];
      LocalRecvRequestTableSize[ii]=0;
    }
                          
    for (auto& it : RequestDataList) {
      RecvRequestTable[it.iSourceThread][LocalRecvRequestTableSize[it.iSourceThread]]=it;
      LocalRecvRequestTableSize[it.iSourceThread]++;
    }
                                                                                                                    
    //gather the global requst number
    MPI_Allgather(LocalRecvRequestTableSize,PIC::nTotalThreads,MPI_INT,GlobalRequestNumberTable,PIC::nTotalThreads,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
    
    MPI_Request *RequestTable=new MPI_Request[2*PIC::nTotalThreads];
    int RequestTableLength=0;
    
    //allocate the request tables and initiate send/recive
    for (int ii=0;ii<PIC::nTotalThreads;ii++) {
      if (GetRecvRequestNumber(ii)!=0) {        
        MPI_Isend(RecvRequestTable[ii],LocalRecvRequestTableSize[ii]*sizeof(cRequestData),MPI_BYTE,ii,9,MPI_GLOBAL_COMMUNICATOR,RequestTable+RequestTableLength);
        RequestTableLength++;
      }
      
      
      if (GetSendRequestNumber(ii)!=0) {
        SendRequestTable[ii]=new cRequestData [GetSendRequestNumber(ii)];

        MPI_Irecv(SendRequestTable[ii],GetSendRequestNumber(ii)*sizeof(cRequestData),MPI_BYTE,ii,9,MPI_GLOBAL_COMMUNICATOR,RequestTable+RequestTableLength);
        RequestTableLength++;
      }
    }
    
    //waite completing the send recive 
    MPI_Waitall(RequestTableLength,RequestTable,MPI_STATUSES_IGNORE);
    
    delete [] RequestTable;
    delete [] LocalSendRequestTableSize;
    delete [] LocalRecvRequestTableSize;
  };
  
  
  //-----------------------------------   EXCHANGE THE REQUESTS: end  ----------------------------------
  
  //-----------------------------------   Allocate Send/Recv Buffers: begin  ----------------------------------
  char *GlobalRecvBuffer=NULL,*GlobalSendBuffer=NULL;
  
  auto DeallocateSendRecvBuffers = [&] () {
    delete [] GlobalRecvBuffer;
    delete [] GlobalSendBuffer;
  };
  
  auto InitSendRecvCornerAssociatedData = [&] {
    int To,ii;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* Node;
    PIC::Mesh::cDataCornerNode *CornerNode;
    
    for (To=0;To<PIC::nTotalThreads;To++) {
      for (ii=0;ii<GetSendRequestNumber(To);ii++) {
        Node=PIC::Mesh::mesh->findAMRnodeWithID(SendRequestTable[To][ii].SourceNodeId);
        CornerNode=Node->block->GetCornerNode(_getCornerNodeLocalNumber(SendRequestTable[To][ii].iSource,SendRequestTable[To][ii].jSource,SendRequestTable[To][ii].kSource));
              
        SendRequestTable[To][ii].CornerAssociatedData=CornerNode->GetAssociatedDataBufferPointer();
      }
    }
  };
  
  auto AllocateSendBuffer = [&] () {
    int offset,ii,To,total_size=0;
  
    for (To=0;To<PIC::nTotalThreads;To++) total_size+=GetSendRequestNumber(To);
    
    total_size*=PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength;
    GlobalSendBuffer=new char [total_size];
    
    for (offset=0,To=0;To<PIC::nTotalThreads;To++) {
      for (ii=0;ii<GetSendRequestNumber(To);ii++) {
        SendRequestTable[To][ii].Buffer=GlobalSendBuffer+offset;
        offset+=PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength;
      }
    }
  };
  
  auto AllocateRecvBuffer = [&] () {
    int offset,ii,From,total_size=0;
  
    for (From=0;From<PIC::nTotalThreads;From++) total_size+=GetRecvRequestNumber(From);
    
    total_size*=PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength;
    GlobalRecvBuffer=new char [total_size];    
    
    for (offset=0,From=0;From<PIC::nTotalThreads;From++) {
      for (ii=0;ii<GetRecvRequestNumber(From);ii++) {
        RecvRequestTable[From][ii].Buffer=GlobalRecvBuffer+offset;
        offset+=PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength;
      }
    }
  };
  
  
  //-----------------------------------   AllocateSend/Recv Buffers: end  ----------------------------------
  
  
  
  
  //-----------------------------------   Pack/Unpack data: begin  -----------------------------------------
  auto PackSendData = [&] () {
    int ii,To,nRequests;
    
    for (To=0;To<PIC::nTotalThreads;To++) {
      nRequests=GetSendRequestNumber(To);

      for (ii=0;ii<nRequests;ii++) {
        memcpy(SendRequestTable[To][ii].Buffer,SendRequestTable[To][ii].CornerAssociatedData,
            PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength);
      }
    }
  };
  
  auto UnpackRecvData = [&] (int From) {
    int ii,nRequests=GetRecvRequestNumber(From);
    
    for (ii=0;ii<nRequests;ii++) {
      PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(
          RecvRequestTable[From][ii].CornerAssociatedData,RecvRequestTable[From][ii].Buffer);
    }
  };
  
  
  //-----------------------------------   Pack/Unpack data: end     -----------------------------------------
  
  
  //-----------------------------------   MPI Send/Recv: begin      -----------------------------------------
  MPI_Request *mpiSendRequestTable=NULL;
  int mpiSendRequestTableLength=0;
  
  
  MPI_Request *mpiRecvRequestTable=NULL;
  int *mpiRecvIndex2ThreadTable=NULL;
  int mpiRecvRequestTableLength=0;
  
  
  auto InitSendData = [&] () {
    int To;
    
    mpiSendRequestTable=new MPI_Request[PIC::nTotalThreads];
    mpiSendRequestTableLength=0;
    
    for (To=0;To<PIC::nTotalThreads;To++) if (GetSendRequestNumber(To)!=0) {    
      MPI_Isend(SendRequestTable[To][0].Buffer,GetSendRequestNumber(To)*PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,To,9,
          MPI_GLOBAL_COMMUNICATOR,mpiSendRequestTable+mpiSendRequestTableLength);
      
      mpiSendRequestTableLength++;
    }
  };
  
  auto WaitSendData = [&] () {
    MPI_Waitall(mpiSendRequestTableLength,mpiSendRequestTable,MPI_STATUS_IGNORE);
    
    delete [] mpiSendRequestTable;
    mpiSendRequestTableLength=0;
  };
  
  auto InitRecvData = [&] () {
    int From;
    
    mpiRecvRequestTable=new MPI_Request[PIC::nTotalThreads];
    mpiRecvIndex2ThreadTable=new int [PIC::nTotalThreads];
    mpiRecvRequestTableLength=0;
    
    for (From=0;From<PIC::nTotalThreads;From++) if (GetRecvRequestNumber(From)!=0) {         
      MPI_Irecv(RecvRequestTable[From][0].Buffer,GetRecvRequestNumber(From)*PIC::Mesh::cDataCornerNode_static_data::totalAssociatedDataLength,MPI_BYTE,From,9,
          MPI_GLOBAL_COMMUNICATOR,mpiRecvRequestTable+mpiRecvRequestTableLength);
      
      mpiRecvIndex2ThreadTable[mpiRecvRequestTableLength]=From;
      mpiRecvRequestTableLength++;
    }
  };
  
  auto WaitRecvData = [&] () {
    int res=-1;
    
    if (mpiRecvRequestTableLength!=0) {
      int index;
      
      MPI_Waitany(mpiRecvRequestTableLength,mpiRecvRequestTable,&index,MPI_STATUS_IGNORE); 
      MPI_Wait(mpiRecvRequestTable+index,MPI_STATUS_IGNORE);

      res=mpiRecvIndex2ThreadTable[index];
      
      mpiRecvRequestTable[index]=mpiRecvRequestTable[mpiRecvRequestTableLength-1];
      mpiRecvIndex2ThreadTable[index]=mpiRecvIndex2ThreadTable[mpiRecvRequestTableLength-1];
      --mpiRecvRequestTableLength; 
      
      if (mpiRecvRequestTableLength==0) {
        delete [] mpiRecvRequestTable;
        delete [] mpiRecvIndex2ThreadTable;
      }
    }
        
    return res;
  };
  
  //-----------------------------------   MPI Send/Recv: end        -----------------------------------------
  
  //function to reset the corner node 'processed' flag
  auto ResetCornerNodeProcessedFlag = [&] () {
    int i,j,k,iface;
    
    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      PIC::Mesh::cDataBlockAMR *block=node->block;
      
      if (block!=NULL) for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        PIC::Mesh::cDataCornerNode *CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));
        
        if (CornerNode!=NULL) CornerNode->SetProcessedFlag(false);
      }
    }
  };
  

  auto ReservePeriodicBCFlags = [&] () {
    periodic_bc_pair_real_block=PIC::Mesh::mesh->rootTree->CheckoutFlag();
    periodic_bc_pair_ghost_block=PIC::Mesh::mesh->rootTree->CheckoutFlag();
    
    if ((periodic_bc_pair_real_block==-1)||(periodic_bc_pair_ghost_block==-1)) exit(__LINE__,__FILE__,"Error: cannot reserve a flag");
  };
  
  
  
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> ResetPeriodicBCFlags;
  
  ResetPeriodicBCFlags = [&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) -> void {
    startNode->SetFlag(false,periodic_bc_pair_real_block);
    startNode->SetFlag(false,periodic_bc_pair_ghost_block);
    
    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;
    
    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      ResetPeriodicBCFlags(downNode);
    }
    
    if (startNode==PIC::Mesh::mesh->rootTree) {
      //loop throught the list of the block pairs used to impose the periodic BC
      int iBlockPair,RealBlockThread,GhostBlockThread;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *RealBlock,*GhostBlock;
      
      //loop through all block pairs
      for (iBlockPair=0;iBlockPair<PIC::BC::ExternalBoundary::Periodic::BlockPairTableLength;iBlockPair++) {
        GhostBlock=PIC::BC::ExternalBoundary::Periodic::BlockPairTable[iBlockPair].GhostBlock;
        RealBlock=PIC::BC::ExternalBoundary::Periodic::BlockPairTable[iBlockPair].RealBlock;
        
        GhostBlock->SetFlag(true,periodic_bc_pair_ghost_block);
        RealBlock->SetFlag(true,periodic_bc_pair_real_block);
      }
    }
  };
  
  //------------------------------------  Init Recv Complete processed CornerNode data:  ----------------
  if (CornerBlockBoundaryNodes::ActiveFlag==false) return; 
  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData==NULL) exit(__LINE__,__FILE__,"Error: the pointer is not allocated"); 
  
  //determine whether the mesh/domain decomposition have been changed
  static int nMeshModificationCounter=-1;
  int localMeshChangeFlag,globalMeshChangeFlag;

  localMeshChangeFlag=(nMeshModificationCounter==PIC::Mesh::mesh->nMeshModificationCounter) ? 0 : 1;
  MPI_Allreduce(&localMeshChangeFlag,&globalMeshChangeFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  nMeshModificationCounter=PIC::Mesh::mesh->nMeshModificationCounter;

  //========================================   main body of the function: communication and processing the data =========================
  //init the flags
  if (globalMeshChangeFlag!=0) {
    ResetCornerNodeProcessedFlag();
  
    ReservePeriodicBCFlags();
    ResetPeriodicBCFlags(PIC::Mesh::mesh->rootTree);
    ReleasePeriodicBCFlags();
    
    ResetCornerNodeProcessedFlag();
    
    
    BuildRequestDataList(PIC::Mesh::mesh->rootTree,&RequestDataList);
    CreateGlobalRequestTable();
  }
  
  AllocateSendBuffer();
  AllocateRecvBuffer();
  
  if (globalMeshChangeFlag!=0) {
    InitSendRecvCornerAssociatedData();
  }
  
  
  PackSendData();
  
  InitSendData();
  InitRecvData();
  
  while (mpiRecvRequestTableLength!=0) {
    int From;
    
    From=WaitRecvData();
    UnpackRecvData(From);
  }
  
  WaitSendData();
  DeallocateSendRecvBuffers();
}  
  
  

//================================================================================
void PIC::Parallel::ProcessBlockBoundaryNodes() {
  #if _PIC_BC__PERIODIC_MODE_== _PIC_BC__PERIODIC_MODE_ON_
  //update associated data accounting for the periodic boundary conditions
  PIC::Parallel::ProcessCornerBlockBoundaryNodes();
  PIC::Parallel::ProcessCenterBlockBoundaryNodes();
  #else

  switch (_PIC_PROCESS_NODE_ASSSOCIATED_DATA_MODE_) {
  case _PIC_PROCESS_NODE_ASSSOCIATED_DATA_MODE__Yuxi_:
    ProcessCornerBlockBoundaryNodes_new();
    ProcessCenterBlockBoundaryNodes_new();
    break;
  case _PIC_PROCESS_NODE_ASSSOCIATED_DATA_MODE__DEFAULT_:
    ProcessCornerBlockBoundaryNodes();
    ProcessCenterBlockBoundaryNodes();
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  #endif
}



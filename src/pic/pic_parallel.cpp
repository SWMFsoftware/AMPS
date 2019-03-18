//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the functions that control the interprocessor communication of the code

#include "pic.h"

long int PIC::Parallel::sendParticleCounter=0,PIC::Parallel::recvParticleCounter=0,PIC::Parallel::IterationNumberAfterRebalancing=0;
double PIC::Parallel::RebalancingTime=0.0,PIC::Parallel::CumulativeLatency=0.0;
double PIC::Parallel::EmergencyLoadRebalancingFactor=3.0;
double PIC::Parallel::Latency=0.0;

//processing 'corner' and 'center' node associated data vectors when perform syncronization
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData=NULL;
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData=NULL;
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::CopyCenterNodeAssociatedData=NULL; 
bool PIC::Parallel::CornerBlockBoundaryNodes::ActiveFlag=false;
void PIC::Parallel::CornerBlockBoundaryNodes::SetActiveFlag(bool flag) {ActiveFlag=flag;}

//default function forprocessing of the corner node associated data
void PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
}

void PIC::Parallel::CornerBlockBoundaryNodes::CopyCenterNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
}

//-----------------------------------------------------------------


/*
//processing 'corner' and 'center' node associated data vectors when perform syncronization
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::ProcessCenterNodeAssociatedData=NULL,PIC::Parallel::ProcessCornerNodeAssociatedData=NULL;
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CopyCenterNodeAssociatedData=PIC::Parallel::CopyCenterNodeAssociatedData_default;
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CopyCornerNodeAssociatedData=PIC::Parallel::CopyCornerNodeAssociatedData_default;

//default function forprocessing of the corner node associated data
void PIC::Parallel::CopyCornerNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
}

void PIC::Parallel::CopyCenterNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
}
*/ 

//====================================================
//Exchange particles between Processors
void PIC::Parallel::ExchangeParticleData() {
  int From,To,i,iFrom,flag;
  long int Particle,NextParticle,newParticle,LocalCellNumber=-1;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *sendNode=NULL,*recvNode=NULL;

  static MPI_Request *SendMessageSizeRequestTable=NULL,*RecvMessageSizeRequestTable=NULL;
  static int *SendMessageLengthTable=NULL,*RecvMessageLengthTable=NULL,*SendMessageLengthProcessTable=NULL,*RecvMessageLengthProcessTable=NULL;
  int SendMessageSizeRequestTableLength=0,RecvMessageSizeRequestTableLength=0;

  static char **RecvParticleDataBuffer=NULL;
  static int *RecvParticleDataBufferLengthTable=NULL;

  MPI_Request RecvPaticleDataRequestTable[PIC::nTotalThreads];
  int RecvPaticleDataRequestTableLength=0;
  int RecvPaticleDataProcessTable[PIC::nTotalThreads];

  static char **SendParticleDataBuffer=NULL;
  static int *SendParticleDataBufferLengthTable=NULL;

  static MPI_Request *SendParticleDataRequestTable=NULL;
  static int SendParticleDataRequestTableLength=0;


  if (SendMessageSizeRequestTable==NULL) {
    SendParticleDataRequestTable=new MPI_Request[PIC::nTotalThreads];

    SendMessageSizeRequestTable=new MPI_Request[PIC::nTotalThreads];
    RecvMessageSizeRequestTable=new MPI_Request[PIC::nTotalThreads];

    SendMessageLengthTable=new int [PIC::nTotalThreads];
    RecvMessageLengthTable=new int [PIC::nTotalThreads];

    SendMessageLengthProcessTable=new int [PIC::nTotalThreads];
    RecvMessageLengthProcessTable=new int [PIC::nTotalThreads];

    RecvParticleDataBuffer=new char* [PIC::nTotalThreads];
    RecvParticleDataBufferLengthTable=new int [PIC::nTotalThreads];

    SendParticleDataBuffer=new char* [PIC::nTotalThreads];
    SendParticleDataBufferLengthTable=new int [PIC::nTotalThreads];

    for (int thread=0;thread<PIC::nTotalThreads;thread++) RecvParticleDataBufferLengthTable[thread]=0,SendParticleDataBufferLengthTable[thread]=0;

  }


  //finish previous send operations in case nor finished yet
  if (SendParticleDataRequestTableLength!=0) {
    MPI_Waitall(SendParticleDataRequestTableLength,SendParticleDataRequestTable,MPI_STATUSES_IGNORE);
    SendParticleDataRequestTableLength=0;
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


  //calculate the number of bytes that will be send
  long int *FirstCellParticleTable;
  int thread;

  for (To=0;To<PIC::Mesh::mesh.nTotalThreads;To++) if ((PIC::ThisThread!=To)&&(PIC::Mesh::mesh.ParallelSendRecvMap[PIC::ThisThread][To]==true)) {
    bool CommunicationInitialed_BLOCK_;
    int iCell,jCell,kCell;

    for (sendNode=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) {
      CommunicationInitialed_BLOCK_=false;
      if (!sendNode->block) continue;
      FirstCellParticleTable=sendNode->block->FirstCellParticleTable;

      for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
        Particle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

        if  (Particle!=-1) {
          if (CommunicationInitialed_BLOCK_==false) {
            SendMessageLengthTable[To]+=sizeof(nodeid)+sizeof(int);

            CommunicationInitialed_BLOCK_=true;
          }

          SendMessageLengthTable[To]+=sizeof(int)+sizeof(LocalCellNumber);

          while (Particle!=-1) {
            SendMessageLengthTable[To]+=sizeof(int)+PIC::ParticleBuffer::ParticleDataLength;

            Particle=PIC::ParticleBuffer::GetNext(Particle);
          }
        }
      }
    }

    //end the part of the sender
    if (SendMessageLengthTable[To]!=0) {
      SendMessageLengthTable[To]+=sizeof(int);
    }

    //the total length of the message to be send
    MPI_Isend(SendMessageLengthTable+To,1,MPI_INT,To,10,MPI_GLOBAL_COMMUNICATOR,SendMessageSizeRequestTable+SendMessageSizeRequestTableLength);
    SendMessageLengthProcessTable[SendMessageSizeRequestTableLength]=To;
    SendMessageSizeRequestTableLength++;
  }

  //recieve the length of the message to be recieved
  for (From=0;From<PIC::Mesh::mesh.nTotalThreads;From++) if ((PIC::ThisThread!=From)&&(PIC::Mesh::mesh.ParallelSendRecvMap[PIC::ThisThread][From]==true)) {
    MPI_Irecv(RecvMessageLengthTable+From,1,MPI_INT,From,10,MPI_GLOBAL_COMMUNICATOR,RecvMessageSizeRequestTable+RecvMessageSizeRequestTableLength);
    RecvMessageLengthProcessTable[RecvMessageSizeRequestTableLength]=From;
    RecvMessageSizeRequestTableLength++;
  }

  //remove unused send buffers
  for (thread=0;thread<PIC::nTotalThreads;thread++) {
    if ((SendMessageLengthTable[thread]==0)&&(SendParticleDataBufferLengthTable[thread]!=0)) {
      //the data buffer is not used
      delete [] SendParticleDataBuffer[thread];
      SendParticleDataBufferLengthTable[thread]=0;
    }
  }

  To=0;

  while ((RecvMessageSizeRequestTableLength!=0)||(RecvPaticleDataRequestTableLength!=0)||(To<PIC::nTotalThreads)) {

    if (RecvMessageSizeRequestTableLength!=0) {
      //check whether sise of the incoming message has been recieved
      MPI_Testany(RecvMessageSizeRequestTableLength,RecvMessageSizeRequestTable,&iFrom,&flag,MPI_STATUS_IGNORE);

      if (flag==true) {
        From=RecvMessageLengthProcessTable[iFrom];
        RecvMessageLengthProcessTable[iFrom]=RecvMessageLengthProcessTable[RecvMessageSizeRequestTableLength-1];
        RecvMessageSizeRequestTable[iFrom]=RecvMessageSizeRequestTable[RecvMessageSizeRequestTableLength-1];
        RecvMessageSizeRequestTableLength--;

        if (RecvMessageLengthTable[From]!=0) {
          //schedule recieving the message
          if (RecvParticleDataBufferLengthTable[From]<RecvMessageLengthTable[From]) {
            if (RecvParticleDataBufferLengthTable[From]!=0) delete [] RecvParticleDataBuffer[From];

            RecvParticleDataBuffer[From]=new char [RecvMessageLengthTable[From]];
            RecvParticleDataBufferLengthTable[From]=RecvMessageLengthTable[From];
          }

          MPI_Irecv(RecvParticleDataBuffer[From],RecvMessageLengthTable[From],MPI_CHAR,From,11,MPI_GLOBAL_COMMUNICATOR,RecvPaticleDataRequestTable+RecvPaticleDataRequestTableLength);
          RecvPaticleDataProcessTable[RecvPaticleDataRequestTableLength]=From;
          RecvPaticleDataRequestTableLength++;
        }
        else {
          //the data buffer is not used. it should be removed in case allocated
          if (RecvParticleDataBufferLengthTable[From]!=0) {
            //the data buffer is not used
            delete [] RecvParticleDataBuffer[From];
            RecvParticleDataBufferLengthTable[From]=0;
          }
        }
      }
    }

    if (RecvPaticleDataRequestTableLength!=0) {
      //check whether new particle data has been recieved
      MPI_Testany(RecvPaticleDataRequestTableLength,RecvPaticleDataRequestTable,&iFrom,&flag,MPI_STATUS_IGNORE);

      if (flag==true) {
        //the new particle data message has been recieve => unpach the particle data
        From=RecvPaticleDataProcessTable[iFrom];

        RecvPaticleDataProcessTable[iFrom]=RecvPaticleDataProcessTable[RecvPaticleDataRequestTableLength-1];
        RecvPaticleDataRequestTable[iFrom]=RecvPaticleDataRequestTable[RecvPaticleDataRequestTableLength-1];
        RecvPaticleDataRequestTableLength--;

        //unpacking the particle data
        if (_PIC_PARTICLE_EXCHANGE_ENFORCE_RECIEVING_ORDER_MODE_==_PIC_MODE_OFF_) {
          int offset=0;
          char *buffer=RecvParticleDataBuffer[From];
          int iCell,jCell,kCell;

          //recieve the data
          Signal=*((int*)(buffer+offset));
          offset+=sizeof(int);

           while (Signal!=_END_COMMUNICATION_SIGNAL_) {
             switch (Signal) {
             case _NEW_BLOCK_ID_SIGNAL_ :
               #if DIM == 3
               nodeid=*((cAMRnodeID*)(buffer+offset));
               #else
               exit(__LINE__,__FILE__,"Error: not implemetned");
               #endif

               offset+=sizeof(nodeid);
               recvNode=PIC::Mesh::mesh.findAMRnodeWithID(nodeid);

               FirstCellParticleTable=recvNode->block->FirstCellParticleTable;

               if (recvNode->block==NULL) exit(__LINE__,__FILE__,"Error: the node is not allocated");
               break;
             case _CENTRAL_NODE_NUMBER_SIGNAL_ :
               //pipe.recv(LocalCellNumber,From);
               LocalCellNumber=*((long int*)(buffer+offset));

               PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(LocalCellNumber,iCell,jCell,kCell);
               offset+=sizeof(long int);

               break;
             case _NEW_PARTICLE_SIGNAL_ :
               newParticle=PIC::ParticleBuffer::GetNewParticle(FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)],true);

               PIC::ParticleBuffer::UnPackParticleData(buffer+offset,newParticle);
               recvParticleCounter++;

               offset+=PIC::ParticleBuffer::ParticleDataLength;
               break;
             default:
               exit(__LINE__,__FILE__,"Error: the option is not recognized");
             }

             Signal=*((int*)(buffer+offset));
             offset+=sizeof(int);
           }

            //end the part of the receiver
           if (offset!=RecvMessageLengthTable[From]) exit(__LINE__,__FILE__,"Error: the amount of recieved data is not consistent");
        }
      }
    }

    //send the particle data message if there is something to send
    if (To<PIC::nTotalThreads) if ((To!=PIC::ThisThread)&&(SendMessageLengthTable[To]!=0)) {
      //there is something to send => pack the data and send it
      bool CommunicationInitialed_BLOCK_;
      int iCell,jCell,kCell;
      bool CellParticleTableModified;

      if (SendParticleDataBufferLengthTable[To]<SendMessageLengthTable[To]) {
        if (SendParticleDataBufferLengthTable[To]!=0) delete [] SendParticleDataBuffer[To];

        SendParticleDataBuffer[To]=new char [SendMessageLengthTable[To]];
        SendParticleDataBufferLengthTable[To]=SendMessageLengthTable[To];
      }

      int offset=0;
      char *buffer=SendParticleDataBuffer[To];

      //reset the proceesed flaf for the blocks to be send
      //send the nodes' data
      for (sendNode=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) {
        CommunicationInitialed_BLOCK_=false;
        if (!sendNode->block) continue;
        FirstCellParticleTable=sendNode->block->FirstCellParticleTable;
        CellParticleTableModified=false;

        for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
          Particle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

          if  (Particle!=-1) {
            LocalCellNumber=_getCenterNodeLocalNumber(iCell,jCell,kCell);

            if (CommunicationInitialed_BLOCK_==false) {
              nodeid=sendNode->AMRnodeID;

              //pipe.send(_NEW_BLOCK_ID_SIGNAL_);
              *((int*)(buffer+offset))=_NEW_BLOCK_ID_SIGNAL_;
              offset+=sizeof(int);

              #if DIM == 3
              *((cAMRnodeID*)(buffer+offset))=nodeid;
              #else
              exit(__LINE__,__FILE__,"Error: not implemetned");
              #endif

              offset+=sizeof(nodeid);

              CommunicationInitialed_BLOCK_=true;
            }

            *((int*)(buffer+offset))=_CENTRAL_NODE_NUMBER_SIGNAL_;
            offset+=sizeof(int);

            *((long int*)(buffer+offset))=LocalCellNumber;
            offset+=sizeof(long int);

            while (Particle!=-1) {
              *((int*)(buffer+offset))=_NEW_PARTICLE_SIGNAL_;
              offset+=sizeof(int);

              PIC::ParticleBuffer::PackParticleData(buffer+offset,Particle);
              offset+=PIC::ParticleBuffer::ParticleDataLength;
              sendParticleCounter++;

              NextParticle=PIC::ParticleBuffer::GetNext(Particle);
              PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(Particle,true);
              Particle=NextParticle;
            }

            FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=-1;
            CellParticleTableModified=true;
          }
        }

      }

      *((int*)(buffer+offset))=_END_COMMUNICATION_SIGNAL_;
      offset+=sizeof(int);

      if (offset!=SendMessageLengthTable[To]) exit(__LINE__,__FILE__,"Error: the data anount to be send is not consistent");

      //end the part of the sender - initiate non-blocks send
      MPI_Isend(buffer,offset,MPI_CHAR,To,11,MPI_GLOBAL_COMMUNICATOR,SendParticleDataRequestTable+SendParticleDataRequestTableLength);
      SendParticleDataRequestTableLength++;
    }

    if (To<PIC::nTotalThreads) To++;
  }

  //in case the order that particles are recieved has to be conserved
  if (_PIC_PARTICLE_EXCHANGE_ENFORCE_RECIEVING_ORDER_MODE_==_PIC_MODE_ON_)  {
    for (From=0;From<PIC::nTotalThreads;From++) if ((From!=PIC::ThisThread)&&(RecvMessageLengthTable[From]!=0)) {
      int offset=0;
      char *buffer=RecvParticleDataBuffer[From];
      int iCell,jCell,kCell;

      //recieve the data
      Signal=*((int*)(buffer+offset));
      offset+=sizeof(int);

       while (Signal!=_END_COMMUNICATION_SIGNAL_) {
         switch (Signal) {
         case _NEW_BLOCK_ID_SIGNAL_ :
           #if DIM == 3
           nodeid=*((cAMRnodeID*)(buffer+offset));
           #else
           exit(__LINE__,__FILE__,"Error: not implemetned");
           #endif

           offset+=sizeof(nodeid);
           recvNode=PIC::Mesh::mesh.findAMRnodeWithID(nodeid);

           FirstCellParticleTable=recvNode->block->FirstCellParticleTable;

           if (recvNode->block==NULL) exit(__LINE__,__FILE__,"Error: the node is not allocated");
           break;
         case _CENTRAL_NODE_NUMBER_SIGNAL_ :
           //pipe.recv(LocalCellNumber,From);
           LocalCellNumber=*((long int*)(buffer+offset));

           PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(LocalCellNumber,iCell,jCell,kCell);
           offset+=sizeof(long int);

           break;
         case _NEW_PARTICLE_SIGNAL_ :
           newParticle=PIC::ParticleBuffer::GetNewParticle(FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)],true);

           PIC::ParticleBuffer::UnPackParticleData(buffer+offset,newParticle);
           recvParticleCounter++;

           offset+=PIC::ParticleBuffer::ParticleDataLength;
           break;
         default:
           exit(__LINE__,__FILE__,"Error: the option is not recognized");
         }

         Signal=*((int*)(buffer+offset));
         offset+=sizeof(int);
       }

        //end the part of the receiver
       if (offset!=RecvMessageLengthTable[From]) exit(__LINE__,__FILE__,"Error: the amount of recieved data is not consistent");
    }
  }
}

void PIC::Parallel::ProcessCornerBlockBoundaryNodes_old() {
  int thread,iThread,i,j,k,iface;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  PIC::Mesh::cDataCornerNode *CornerNode;
  char *CornerNodeAssociatedData;
  PIC::Mesh::cDataBlockAMR *block;
  MPI_Status status;

  if (CornerBlockBoundaryNodes::ActiveFlag==false) return;

  const int iFaceMin[6]={0,_BLOCK_CELLS_X_,0,0,0,0};
  const int iFaceMax[6]={0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_};

  const int jFaceMin[6]={0,0,0,_BLOCK_CELLS_Y_,0,0};
  const int jFaceMax[6]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_};

  const int kFaceMin[6]={0,0,0,0,0,_BLOCK_CELLS_Z_};
  const int kFaceMax[6]={_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,0,_BLOCK_CELLS_Z_};

  static const int BlockCornerTable[8][3]={
      {0,0,0}, {_BLOCK_CELLS_X_,0,0}, {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0}, {0,_BLOCK_CELLS_Y_,0},
      {0,0,_BLOCK_CELLS_Z_}, {_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_}, {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}, {0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}
  };

  static const int BlockCornerOffsetTable[8][3]={
      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},

      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_}
  };

  static const int BlockEdgeOffsetTable[12][3]={
      {0,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {0,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {0,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {0,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},

      {-_BLOCK_CELLS_X_,0,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,0,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,0,+_BLOCK_CELLS_Z_},
      {-_BLOCK_CELLS_X_,0,+_BLOCK_CELLS_Z_},

      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,0},
      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,0},
      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,0},
      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,0}
  };

  static const int BlockFaceOffsetTable[6][3]={
      {-_BLOCK_CELLS_X_,0,0},{+_BLOCK_CELLS_X_,0,0},
      {0,-_BLOCK_CELLS_Y_,0},{0,+_BLOCK_CELLS_Y_,0},
      {0,0,-_BLOCK_CELLS_Z_},{0,0,+_BLOCK_CELLS_Z_}
  };

  static const int iMinEdgeTable[12]={
      0,0,0,0,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0
  };

  static const int iMaxEdgeTable[12]={
      _BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0
  };

  static const int jMinEdgeTable[12]={
      0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,
      0,0,0,0,
      0,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_
  };

  static const int jMaxEdgeTable[12]={
      0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,
      _BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,
      0,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_
  };

  static const int kMinEdgeTable[12]={
    0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
    0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
    0,0,0,0
  };

  static const int kMaxEdgeTable[12]={
      0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
      0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
      _BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_
  };

  struct cStencilElement {
    int StencilLength;
    int StencilThreadTable[80];
    int iCornerNode,jCornerNode,kCornerNode;
    cAMRnodeID nodeid;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
    PIC::Mesh::cDataCornerNode *CornerNode;
    char *AssociatedDataPointer;
  };

  int iStencil;
  static int StencilTableLength=0;
  static cStencilElement *StencilTable=NULL;

  //SendTablePBC,RecvTablePBC,SendTableLengthPBC,RecvTableLengthPBC used for data exchange in case the periodic boundary conditions are enforsed
  struct cSendRecvCornerNodeDataElement {
    int i,j,k;
    cAMRnodeID SendBlockID,RecvBlockID;
    PIC::Mesh::cDataCornerNode *CornerNode;
    char *AssociatedDataPointer;
    int OperationID;
  };

  //pocesseing of the associated data of the nodes located at the boundary of the "real" computational domain in case the periodic boundary conditions are in use
  struct cNodeSetElement {
    cAMRnodeID nodeid;
    int i,j,k;
  };

  struct cNodeSet {
    cNodeSetElement NodeTable[80];
    int NodeTableLength;
  };

  list<cNodeSet> NodeSetList;


  struct cStencilPoint {
    int ThreadTable[80];
    cNodeSetElement NodeSetTable[80];
    int ThreadTableLength;

    char *AssociatedDataPointer;
  };

  struct cStencilPBC {
    cStencilPoint PointTable[80];
    int nStencilPoints;
    int ProcessingThread;
    int SourceThreadTable[80];

    //threads to which the processed associated vector will be send out back
    int InvolvedThreadTable[80],InvolvedThreadTableLength;
    bool InvolvedFlag;
  };

  struct cBlockTable {
    cAMRnodeID GhostBloks[80],RealBlockPair[80];
    int iCorner[80],jCorner[80],kCorner[80];
    int BlockTableLength;
  };

  static int StencilTablePBCLength=0;
  static cStencilPBC *StencilTablePBC=NULL;




  //generate a new stencil table
  static int nMeshModificationCounter=-1;


  //function to reset the corner node 'processed' flag
  auto ResetProcessedFlag = [&] () {
    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      PIC::Mesh::cDataBlockAMR *block=node->block;

      if (block!=NULL) for (i=0;i<_BLOCK_CELLS_X_+1;i++) for (j=0;j<_BLOCK_CELLS_Y_+1;j++) for (k=0;k<_BLOCK_CELLS_Z_+1;k++) {
        PIC::Mesh::cDataCornerNode *CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));

        if (CornerNode!=NULL) CornerNode->SetProcessedFlag(false);
      }
    }
  };



  auto VerifyAllocatedNeighbourBlock = [&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    bool FoundFlag=false;
    int iNeighbour;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;

    //check faces
    for (iNeighbour=0;iNeighbour<6*4;iNeighbour++) if ((neibNode=node->neibNodeFace[iNeighbour])!=NULL) if (neibNode->block!=NULL)  {
      FoundFlag=true;
      break;
    }

    //check corners
    if (FoundFlag==false) for (iNeighbour=0;iNeighbour<8;iNeighbour++) if ((neibNode=node->neibNodeCorner[iNeighbour])!=NULL) if (neibNode->block!=NULL) {
      FoundFlag=true;
      break;
    }

    //check edges
    if (FoundFlag==false) for (iNeighbour=0;iNeighbour<12*2;iNeighbour++) if ((neibNode=node->neibNodeEdge[iNeighbour])!=NULL) if (neibNode->block!=NULL) {
      FoundFlag=true;
      break;
    }

    return FoundFlag;
  };


  //find a corner node in neighbour blocks
  auto FindCornerNbN = [&] (int &iNode,int &jNode,int &kNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &NeibNode,PIC::Mesh::cDataCornerNode *CornerNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    bool res=false;
    int i,j,k,iface,iNeighbour;
    PIC::Mesh::cDataBlockAMR *block;

    //check faces
    for (iNeighbour=0;iNeighbour<6*4;iNeighbour++) if (node->neibNodeFace[iNeighbour]!=NULL) {
      block=node->neibNodeFace[iNeighbour]->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
          res=true;
          iNode=i,jNode=j,kNode=k;
          NeibNode=node->neibNodeFace[iNeighbour];
          return res;
        }
      }
    }

    //check corners
    for (iNeighbour=0;iNeighbour<8;iNeighbour++) if (node->neibNodeCorner[iNeighbour]!=NULL) {
      block=node->neibNodeCorner[iNeighbour]->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
          res=true;
          iNode=i,jNode=j,kNode=k;
          NeibNode=node->neibNodeCorner[iNeighbour];
          return res;
        }
      }
    }

    //check edges
    for (iNeighbour=0;iNeighbour<12*2;iNeighbour++) if (node->neibNodeEdge[iNeighbour]!=NULL) {
      block=node->neibNodeEdge[iNeighbour]->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
          res=true;
          iNode=i,jNode=j,kNode=k;
          NeibNode=node->neibNodeEdge[iNeighbour];
          return res;
        }
      }
    }

    NeibNode=NULL;
    return res;
  };

  auto GetCornerNode = [&] (int &iNode,int &jNode,int &kNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &NodeOut,cAMRnodeID NodeidIn) {
    PIC::Mesh::cDataCornerNode *res=NULL;
    int iface,iedge,icorner,i,j,k,iNeighbour;
    PIC::Mesh::cDataBlockAMR *block;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode,*node=PIC::Mesh::mesh.findAMRnodeWithID(NodeidIn);


    //check self
    if ((block=node->block)!=NULL) {
      res=block->GetCornerNode(_getCornerNodeLocalNumber(iNode,jNode,kNode));
      NodeOut=node;
      return res;
    }
    else {
      //search connections through faces
      for (iface=0;iface<6;iface++) if ((neibNode=node->GetNeibFace(iface,0,0))!=NULL) if (neibNode->block!=NULL) {
        if (node->RefinmentLevel!=neibNode->RefinmentLevel) exit(__LINE__,__FILE__,"Error: not implemented for the case when neibours has different resolution levels");

        i=iNode-BlockFaceOffsetTable[iface][0];

        if ((0<=i)&&(i<=_BLOCK_CELLS_X_)) {
          j=jNode-BlockFaceOffsetTable[iface][1];

          if ((0<=j)&&(j<=_BLOCK_CELLS_Y_)) {
            k=kNode-BlockFaceOffsetTable[iface][2];

            if ((0<=k)&&(k<=_BLOCK_CELLS_Z_)) {
              if ((res=neibNode->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                NodeOut=neibNode;
                iNode=i,jNode=j,kNode=k;
                return res;
              }
            }
          }
        }
      }

      //search connections through edges
      for (iedge=0;iedge<12;iedge++) if ((neibNode=node->GetNeibEdge(iedge,0))!=NULL) if (neibNode->block!=NULL) {
        if (node->RefinmentLevel!=neibNode->RefinmentLevel) exit(__LINE__,__FILE__,"Error: not implemented for the case when neibours has different resolution levels");

        i=iNode-BlockEdgeOffsetTable[iedge][0];

        if ((0<=i)&&(i<=_BLOCK_CELLS_X_)) {
          j=jNode-BlockEdgeOffsetTable[iedge][1];

          if ((0<=j)&&(j<=_BLOCK_CELLS_Y_)) {
            k=kNode-BlockEdgeOffsetTable[iedge][2];

            if ((0<=k)&&(k<=_BLOCK_CELLS_Z_)) {
              if ((res=neibNode->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                NodeOut=neibNode;
                iNode=i,jNode=j,kNode=k;
                return res;
              }
            }
          }
        }
      }

      //search connections through corners
      for (icorner=0;icorner<8;icorner++) if ((neibNode=node->GetNeibCorner(icorner))!=NULL) if (neibNode->block!=NULL) {
        if (node->RefinmentLevel!=neibNode->RefinmentLevel) exit(__LINE__,__FILE__,"Error: not implemented for the case when neibours has different resolution levels");

        i=iNode-BlockCornerOffsetTable[icorner][0];

        if ((0<=i)&&(i<=_BLOCK_CELLS_X_)) {
          j=jNode-BlockCornerOffsetTable[icorner][1];

          if ((0<=j)&&(j<=_BLOCK_CELLS_Y_)) {
            k=kNode-BlockCornerOffsetTable[icorner][2];

            if ((0<=k)&&(k<=_BLOCK_CELLS_Z_)) {
              if ((res=neibNode->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                NodeOut=neibNode;
                iNode=i,jNode=j,kNode=k;
                return res;
              }
            }
          }
        }
      }

    }

    res=NULL;
    return res;
  };


  auto VerifyBoundaryCornerNode = [&] (PIC::Mesh::cDataCornerNode *CornerNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    bool res=false;
    int iNeighbour,i,j,k,iface;
    bool EnternalBlockFlag,StartBlockExternalFlag;
    PIC::Mesh::cDataBlockAMR *block;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode,*neibNode_LastProcessed=NULL;

    StartBlockExternalFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(node)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

    //check faces
    for (iNeighbour=0;iNeighbour<6*4;iNeighbour++) if ((neibNode=node->neibNodeFace[iNeighbour])!=NULL) if (neibNode!=neibNode_LastProcessed) if ((block=neibNode->block)!=NULL) {
      neibNode_LastProcessed=neibNode;
      EnternalBlockFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if ( ((EnternalBlockFlag==true)&&(StartBlockExternalFlag==false)) || ((EnternalBlockFlag==false)&&(StartBlockExternalFlag==true)) ) {
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
            //the 'corner' node is at the boundary of the 'real' domain
            res=true;
            return res;
          }
        }
      }
    }

    //check corners
    for (iNeighbour=0;iNeighbour<8;iNeighbour++) if ((neibNode=node->neibNodeCorner[iNeighbour])!=NULL) if (neibNode!=neibNode_LastProcessed) if ((block=neibNode->block)!=NULL) {
      neibNode_LastProcessed=neibNode;
      EnternalBlockFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if ( ((EnternalBlockFlag==true)&&(StartBlockExternalFlag==false)) || ((EnternalBlockFlag==false)&&(StartBlockExternalFlag==true)) ) {
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
            //the 'corner' node is at the boundary of the 'real' domain
            res=true;
            return res;
          }
        }
      }
    }

    //check edges
    for (iNeighbour=0;iNeighbour<12*2;iNeighbour++) if ((neibNode=node->neibNodeEdge[iNeighbour])!=NULL) if (neibNode!=neibNode_LastProcessed) if ((block=neibNode->block)!=NULL) {
      neibNode_LastProcessed=neibNode;
      EnternalBlockFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if ( ((EnternalBlockFlag==true)&&(StartBlockExternalFlag==false)) || ((EnternalBlockFlag==false)&&(StartBlockExternalFlag==true)) ) {
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
            //the 'corner' node is at the boundary of the 'real' domain
            res=true;
            return res;
          }
        }
      }
    }

    return res;
  };

  auto TestCheckFace = [&] (int iface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode;
    int i,j;

    //check face neighbors
    for (i=0;i<2;i++) for (j=0;j<2;j++) {
      if ((neibNode=node->GetNeibFace(iface,i,j))!=NULL) {
        if (node->Thread!=neibNode->Thread) {
         //the check need to be performed
         return true;
        }

        if (node->RefinmentLevel==neibNode->RefinmentLevel) break;
      }
    }

    //check edge neighbors
    static const int EdgeFaceTable[6][4]={{4,11,7,8},{5,10,6,9}, {0,9,3,8},{1,10,2,11}, {0,5,1,4},{3,6,2,7}};

    for (int iedge=0;iedge<4;iedge++) {
      for (i=0;i<2;i++) {
        if ((neibNode=node->GetNeibEdge(EdgeFaceTable[iface][iedge],i))!=NULL) {
          if (node->Thread!=neibNode->Thread) {
            //the check need to be performed
            return true;
          }

          if (node->RefinmentLevel==neibNode->RefinmentLevel) break;
        }
      }
    }

    //check corner neighbors
    static const int CornerFaceTable[6][4]={{0,3,7,4},{1,2,6,5}, {0,1,5,4},{3,2,6,7}, {0,1,2,3},{4,5,6,7}};

    for (i=0;i<4;i++) if ((neibNode=node->GetNeibCorner(CornerFaceTable[iface][i]))!=NULL) if (node->Thread!=neibNode->Thread) {
      //the check need to be performed
      return true;
    }

    return false;
  };

  auto TestRealBoundaryFace = [&] (int iface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    int i,j;
    PIC::Mesh::cDataCornerNode *CornerNode;
    bool StartBlockExternalFlag,NeibBlockExternalFlag;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode;

    StartBlockExternalFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(node)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

    for (i=0;i<2;i++) for (j=0;j<2;j++) if ((neibNode=node->GetNeibFace(iface,i,j))!=NULL) {
      NeibBlockExternalFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if (((StartBlockExternalFlag==true)&&(NeibBlockExternalFlag==false)) || ((StartBlockExternalFlag==false)&&(NeibBlockExternalFlag==true))) {
        return true;
      }
      else {
        return false;
      }
    }

    return false;
  };


  //determine whether the mesh/domain decomposition have been changed
  int localMeshChangeFlag,globalMeshChangeFlag;

  localMeshChangeFlag=(nMeshModificationCounter==PIC::Mesh::mesh.nMeshModificationCounter) ? 0 : 1;
  MPI_Allreduce(&localMeshChangeFlag,&globalMeshChangeFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);


  if (globalMeshChangeFlag!=0) {
    //the mesh or the domain decomposition has been modified. Need to create a new communucation table
    int NewTableLength=0;
    int iNode,jNode,kNode,FlagSum;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* ActualNode;
    bool meshModifiedFlag_CountMeshElements=PIC::Mesh::mesh.meshModifiedFlag_CountMeshElements;
    double StartTime=MPI_Wtime();

    //reset the flags
    ResetProcessedFlag();

    //determine the new length of the table
    for (iStencil=0,node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      if ((_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_)||(PIC::Mesh::mesh.ExternalBoundaryBlock(node)!=_EXTERNAL_BOUNDARY_BLOCK_)) {
        int flag;
        bool TemporaralyAllocatedBlock=false;
        cAMRnodeID nodeid=node->AMRnodeID;

        //for (int i=0;i<_BLOCK_CELLS_X_+1;i++) for (int j=0;j<_BLOCK_CELLS_Y_+1;j++)  for (int k=0;k<_BLOCK_CELLS_Z_+1;k++) {
        for (iface=0;iface<6;iface++) if (TestCheckFace(iface,node)==true) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          //determine whether the corner node has been already processes
          if (node->Thread==PIC::ThisThread) {
            CornerNode=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));
            flag=(CornerNode->TestProcessedFlag()==true) ? 0 : 1;
          }

          MPI_Bcast(&flag,1,MPI_INT,node->Thread,MPI_GLOBAL_COMMUNICATOR);

          if (flag==0) continue;

          //the corner node has not been processed yet
          iNode=i,jNode=j,kNode=k;

          if ((CornerNode=GetCornerNode(iNode,jNode,kNode,ActualNode,nodeid))!=NULL) {
            //the corner node exists at the current MPI process
            CornerNode->SetProcessedFlag(true);
            flag=1;
          }
          else flag=0;

          MPI_Allreduce(&flag,&FlagSum,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

          if (FlagSum>1) {
            NewTableLength++;
          }
        }
      }
    }

    //allocate the new Stencile Table
    if (StencilTableLength!=0) delete [] StencilTable;

    StencilTableLength=NewTableLength;
    StencilTable=new cStencilElement[NewTableLength];

    //reset the 'processed' flag
    ResetProcessedFlag();

    //populate the Stencil Table
    for (iStencil=0,node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      if ((_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_)||(PIC::Mesh::mesh.ExternalBoundaryBlock(node)!=_EXTERNAL_BOUNDARY_BLOCK_)) {
        int flag;
        bool TemporaralyAllocatedBlock=false;
        cAMRnodeID nodeid=node->AMRnodeID;

        //for (int i=0;i<_BLOCK_CELLS_X_+1;i++) for (int j=0;j<_BLOCK_CELLS_Y_+1;j++)  for (int k=0;k<_BLOCK_CELLS_Z_+1;k++) {
        for (iface=0;iface<6;iface++) if (TestCheckFace(iface,node)==true) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          //determine whether the corner node has been already processes
          if (node->Thread==PIC::ThisThread) {
            CornerNode=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));
            flag=(CornerNode->TestProcessedFlag()==true) ? 0 : 1;
          }

          MPI_Bcast(&flag,1,MPI_INT,node->Thread,MPI_GLOBAL_COMMUNICATOR);

          if (flag==0) continue;

          //the corner node has not been processed yet
          iNode=i,jNode=j,kNode=k;

          if ((CornerNode=GetCornerNode(iNode,jNode,kNode,ActualNode,nodeid))!=NULL) {
            //the corner node exists at the current MPI process
            CornerNode->SetProcessedFlag(true);
            flag=1;
          }
          else flag=0;

          MPI_Allreduce(&flag,&FlagSum,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

          if (FlagSum>1) {
            StencilTable[iStencil].StencilLength=0;

            StencilTable[iStencil].iCornerNode=iNode,StencilTable[iStencil].jCornerNode=jNode,StencilTable[iStencil].kCornerNode=kNode;
            StencilTable[iStencil].nodeid=ActualNode->AMRnodeID;
            StencilTable[iStencil].node=ActualNode;
            StencilTable[iStencil].CornerNode=CornerNode;
            StencilTable[iStencil].AssociatedDataPointer=(CornerNode!=NULL) ? CornerNode->GetAssociatedDataBufferPointer() : NULL;

            //combine information of which MPI processes has a copy of the current corner node
            int FlagTable[PIC::nTotalThreads];

            MPI_Allgather(&flag,1,MPI_INT,FlagTable,1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

            for (thread=0;thread<PIC::nTotalThreads;thread++) if (FlagTable[thread]==1) {
              StencilTable[iStencil].StencilThreadTable[StencilTable[iStencil].StencilLength++]=thread;
            }

            iStencil++;
          }
        }
      }
    }

    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      //In case when the periodic boundary conditions are inforced
      //additional information exchange table need to be generated to link points located at the boundary of the "real" domain
      cBlockTable BlockTable; //contains information of all blocks that containes a given corner node


      //pass twice throug the next loop: during the first pass the value of 'StencilTablePBCLength' will be evaluated, and during the second pass the table will be populated
      for (int iMainPass=0;iMainPass<2;iMainPass++) {

        if (iMainPass==0) {
          if (StencilTablePBCLength!=0) {
            delete [] StencilTablePBC;

            StencilTablePBC=NULL;
            StencilTablePBCLength=0;
          }
        }
        else if (PIC::ThisThread==0) {
          StencilTablePBC=new cStencilPBC[StencilTablePBCLength];
          StencilTablePBCLength=0;
        }


        ResetProcessedFlag();

        for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)  {
          if (PIC::Mesh::mesh.ExternalBoundaryBlock(node)==_EXTERNAL_BOUNDARY_BLOCK_)  {
            for (iface=0;iface<6;iface++) if (TestRealBoundaryFace(iface,node)==true) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
              //the analysis will be performed by the MPI process that the block belongs to
              int ii,jj,kk,ff,iNeighbour,BoundaryNodeFlag=false;
              cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;
              PIC::Mesh::cDataBlockAMR *neibBlock;
              bool found;

              if (node->Thread==PIC::ThisThread) {
                BlockTable.BlockTableLength=0;
                CornerNode=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));

                //determine whether the point is at the boundary
                BoundaryNodeFlag=VerifyBoundaryCornerNode(CornerNode,node);

                if ((BoundaryNodeFlag==true)&&(CornerNode->TestProcessedFlag()==false)) {
                  //the corner node has not been used yet. The set of the corner nodes that are connected due to enforcing the periodic boundary conditions is not determined yet

                  BlockTable.iCorner[0]=i;
                  BlockTable.jCorner[0]=j;
                  BlockTable.kCorner[0]=k;
                  BlockTable.GhostBloks[0]=node->AMRnodeID;
                  BlockTable.BlockTableLength=1;


                  //find all 'ghost' nodes that contains that corner node
                  //check faces
                  static const int FaceNeibBlockOffset[6][3]={
                      {-_BLOCK_CELLS_X_,0,0}, {_BLOCK_CELLS_X_,0,0},
                      {0,-_BLOCK_CELLS_Y_,0}, {0,_BLOCK_CELLS_Y_,0},
                      {0,0,-_BLOCK_CELLS_Z_}, {0,0,_BLOCK_CELLS_Z_}
                  };


                  for (int iiface=0;iiface<6;iiface++) if ((neibNode=node->GetNeibFace(iiface,0,0))!=NULL) if (PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
                    neibBlock=neibNode->block;

                    if (node->RefinmentLevel!=neibNode->RefinmentLevel) {
                      exit(__LINE__,__FILE__,"Error: blocks has different refinement levels");
                    }

                    for (found=false,ii=0;ii<BlockTable.BlockTableLength;ii++) if (BlockTable.GhostBloks[ii]==neibNode->AMRnodeID) {
                      found=true;
                      break;
                    }

                    if (found==true) continue;

                    if (neibBlock!=NULL) {
                      ii=i-FaceNeibBlockOffset[iiface][0];

                      if ((0<=ii)&&(ii<=_BLOCK_CELLS_X_)) {
                        jj=j-FaceNeibBlockOffset[iiface][1];

                        if ((0<=jj)&&(jj<=_BLOCK_CELLS_Y_)) {
                          kk=k-FaceNeibBlockOffset[iiface][2];

                          if ((0<=kk)&&(kk<=_BLOCK_CELLS_Z_)) {
                            if (CornerNode==neibBlock->GetCornerNode(_getCornerNodeLocalNumber(ii,jj,kk))) {
                              //a new block has been found that contains tested 'corner' node
                              BlockTable.iCorner[BlockTable.BlockTableLength]=ii;
                              BlockTable.jCorner[BlockTable.BlockTableLength]=jj;
                              BlockTable.kCorner[BlockTable.BlockTableLength]=kk;
                              BlockTable.GhostBloks[BlockTable.BlockTableLength]=neibNode->AMRnodeID;
                              BlockTable.BlockTableLength++;
                            }
                          }
                        }
                      }

                    }
                  }

                  //check corners
                  static const int CornerNeibBlockOffset[8][3]={
                      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
                      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
                      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
                      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},

                      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
                      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
                      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
                      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_}
                  };


                  for (int icorner=0;icorner<8;icorner++) if ((neibNode=node->GetNeibCorner(icorner))!=NULL) if (PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
                    neibBlock=neibNode->block;

                    if (node->RefinmentLevel!=neibNode->RefinmentLevel) {
                      exit(__LINE__,__FILE__,"Error: blocks has different refinement levels");
                    }
                      
                    for (found=false,ii=0;ii<BlockTable.BlockTableLength;ii++) if (BlockTable.GhostBloks[ii]==neibNode->AMRnodeID) {
                      found=true;
                      break;
                    }

                    if (found==true) continue;

                    if (neibBlock!=NULL) {
                      ii=i-CornerNeibBlockOffset[icorner][0];

                      if ((0<=ii)&&(ii<=_BLOCK_CELLS_X_)) {
                        jj=j-CornerNeibBlockOffset[icorner][1];

                        if ((0<=jj)&&(jj<=_BLOCK_CELLS_Y_)) {
                          kk=k-CornerNeibBlockOffset[icorner][2];

                          if ((0<=kk)&&(kk<=_BLOCK_CELLS_Z_)) {
                            if (CornerNode==neibBlock->GetCornerNode(_getCornerNodeLocalNumber(ii,jj,kk))) {
                              //a new block has been found that contains tested 'corner' node
                              BlockTable.iCorner[BlockTable.BlockTableLength]=ii;
                              BlockTable.jCorner[BlockTable.BlockTableLength]=jj;
                              BlockTable.kCorner[BlockTable.BlockTableLength]=kk;
                              BlockTable.GhostBloks[BlockTable.BlockTableLength]=neibNode->AMRnodeID;
                              BlockTable.BlockTableLength++;
                            }
                          }
                        }
                      }
                    }
                  }


                  //check edges
                  static const int EdgeNeibBlockOffset[12][3]={
                      {0,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},{0,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},{0,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},{0,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
                      {-_BLOCK_CELLS_X_,0,-_BLOCK_CELLS_Z_},{+_BLOCK_CELLS_X_,0,-_BLOCK_CELLS_Z_},{+_BLOCK_CELLS_X_,0,+_BLOCK_CELLS_Z_},{-_BLOCK_CELLS_X_,0,+_BLOCK_CELLS_Z_},
                      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,0},{+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,0},{+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,0},{-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,0}
                  };

                  for (int iedge=0;iedge<12;iedge++) if ((neibNode=node->GetNeibEdge(iedge,0))!=NULL) if (PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
                    neibBlock=neibNode->block;

                    if (node->RefinmentLevel!=neibNode->RefinmentLevel) {
                      exit(__LINE__,__FILE__,"Error: blocks has different refinement levels");
                    }

                    for (found=false,ii=0;ii<BlockTable.BlockTableLength;ii++) if (BlockTable.GhostBloks[ii]==neibNode->AMRnodeID) {
                      found=true;
                      break;
                    }

                    if (found==true) continue;

                    if (neibBlock!=NULL) {
                      ii=i-EdgeNeibBlockOffset[iedge][0];

                      if ((0<=ii)&&(ii<=_BLOCK_CELLS_X_)) {
                        jj=j-EdgeNeibBlockOffset[iedge][1];

                        if ((0<=jj)&&(jj<=_BLOCK_CELLS_Y_)) {
                          kk=k-EdgeNeibBlockOffset[iedge][2];

                          if ((0<=kk)&&(kk<=_BLOCK_CELLS_Z_)) {
                            if (CornerNode==neibBlock->GetCornerNode(_getCornerNodeLocalNumber(ii,jj,kk))) {
                              //a new block has been found that contains tested 'corner' node
                              BlockTable.iCorner[BlockTable.BlockTableLength]=ii;
                              BlockTable.jCorner[BlockTable.BlockTableLength]=jj;
                              BlockTable.kCorner[BlockTable.BlockTableLength]=kk;
                              BlockTable.GhostBloks[BlockTable.BlockTableLength]=neibNode->AMRnodeID;
                              BlockTable.BlockTableLength++;
                            }
                          }
                        }
                      }
                    }
                  }

                  //determine 'real' blocks tha correspond to those in 'BlockTable'
                  for (ii=0;ii<BlockTable.BlockTableLength;ii++) {
                    BlockTable.RealBlockPair[ii]=PIC::BC::ExternalBoundary::Periodic::findCorrespondingRealBlock(PIC::Mesh::mesh.findAMRnodeWithID(BlockTable.GhostBloks[ii]))->AMRnodeID;
                  }
                }
              }

              //Broadcast the 'BlockTable' to other MPI processes
              int iBlock;
              int DataRequestFlagTable[PIC::nTotalThreads];
              cNodeSetElement Set;
              cStencilPBC NewStencilElementPBC;

              MPI_Bcast(&BlockTable,sizeof(cBlockTable),MPI_BYTE,node->Thread,MPI_GLOBAL_COMMUNICATOR);

              NewStencilElementPBC.nStencilPoints=0;

              //loop through the 'BlockTable'
              for (int ipass=0;ipass<2;ipass++) for (iBlock=0;iBlock<((ipass==0) ? 1: BlockTable.BlockTableLength);iBlock++) {
                switch (ipass) {
                case 0:
                  Set.i=BlockTable.iCorner[iBlock];
                  Set.j=BlockTable.jCorner[iBlock];
                  Set.k=BlockTable.kCorner[iBlock];
                  Set.nodeid=BlockTable.GhostBloks[iBlock];
                  break;
                default:
                  Set.i=BlockTable.iCorner[iBlock];
                  Set.j=BlockTable.jCorner[iBlock];
                  Set.k=BlockTable.kCorner[iBlock];
                  Set.nodeid=BlockTable.RealBlockPair[iBlock];
                }

                //determine whether a corner node descrived by 'Set' exist in the current MPI process
                //if the corner node exists -> set the processed flag state "true"
                int CornerNodeExistFlag;
                int CornerNodeRecvFlagTable[PIC::nTotalThreads],CornerNodeSendFlagTable[PIC::nTotalThreads];
                PIC::Mesh::cDataCornerNode *CornerNode;
                cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* NodeOut;

                CornerNode=GetCornerNode(Set.i,Set.j,Set.k,NodeOut,Set.nodeid);

                if (CornerNode!=NULL) {
                  //the corner exists in the currect MPI process
                  if (CornerNode->TestProcessedFlag()==false) {
                    CornerNode->SetProcessedFlag(true);
                    CornerNodeExistFlag=true;

                    Set.nodeid=NodeOut->AMRnodeID;
                  }
                  else CornerNodeExistFlag=false;
                }
                else {
                  CornerNodeExistFlag=false;
                }

                //gather information from all MPI processes at the root MPI process
                MPI_Gather(&CornerNodeExistFlag,1,MPI_INT,CornerNodeRecvFlagTable,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

                if (PIC::ThisThread==0) {
                  cStencilPoint StencilPoint;

                  StencilPoint.ThreadTableLength=0;

                  if (CornerNodeRecvFlagTable[0]==true) {
                    //root thread has the corner node
                    StencilPoint.NodeSetTable[0]=Set;
                    StencilPoint.ThreadTable[0]=PIC::ThisThread;
                    StencilPoint.ThreadTableLength=1;
                  }

                  for (thread=1;thread<PIC::nTotalThreads;thread++) if (CornerNodeRecvFlagTable[thread]==true) {
                    MPI_Status status;

                    MPI_Recv(&Set,sizeof(cNodeSetElement),MPI_BYTE,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
                    StencilPoint.NodeSetTable[StencilPoint.ThreadTableLength]=Set;
                    StencilPoint.ThreadTable[StencilPoint.ThreadTableLength]=thread;
                    StencilPoint.ThreadTableLength++;
                  }

                  if (StencilPoint.ThreadTableLength!=0) NewStencilElementPBC.PointTable[NewStencilElementPBC.nStencilPoints++]=StencilPoint;
                }
                else if (CornerNodeExistFlag==true) {
                  MPI_Send(&Set,sizeof(cNodeSetElement),MPI_BYTE,0,0,MPI_GLOBAL_COMMUNICATOR);
                }
              }

              //create a new entry to the stencil table
              if ((PIC::ThisThread==0)&&(NewStencilElementPBC.nStencilPoints!=0)) {
                //determine thread the will process the data
                NewStencilElementPBC.ProcessingThread=NewStencilElementPBC.PointTable[0].ThreadTable[0];

                //determine thread that is the source of the data for each point of the stencil
                for (ii=0;ii<NewStencilElementPBC.nStencilPoints;ii++) NewStencilElementPBC.SourceThreadTable[ii]=NewStencilElementPBC.PointTable[ii].ThreadTable[0];

                //determine all threads that are involved into the communications
                NewStencilElementPBC.InvolvedThreadTableLength=0;

                for (ii=0;ii<NewStencilElementPBC.nStencilPoints;ii++) {
                  for (jj=0;jj<NewStencilElementPBC.PointTable[ii].ThreadTableLength;jj++) {
                    bool found=false;

                    for (kk=0;kk<NewStencilElementPBC.InvolvedThreadTableLength;kk++) if (NewStencilElementPBC.PointTable[ii].ThreadTable[jj]==NewStencilElementPBC.InvolvedThreadTable[kk]) {
                      found=true;
                      break;
                    }

                    if (found==false) NewStencilElementPBC.InvolvedThreadTable[NewStencilElementPBC.InvolvedThreadTableLength++]=NewStencilElementPBC.PointTable[ii].ThreadTable[jj];
                  }
                }

                //add the new stencil to the stencil list
                if (iMainPass==0) {
                  StencilTablePBCLength++;
                }
                else {
                  StencilTablePBC[StencilTablePBCLength++]=NewStencilElementPBC;
                }
              }
            }
          }
        }
      }

    if (PIC::ThisThread==0) {
      //create the stencil table
      int i;

      MPI_Bcast(&StencilTablePBCLength,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

      CMPI_channel pipe;

      pipe.init(5000000);
      pipe.openBcast(0);

      for (i=0;i<StencilTablePBCLength;i++) pipe.send(StencilTablePBC[i]);

      pipe.closeBcast();
    }
    else {
      MPI_Bcast(&StencilTablePBCLength,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

      StencilTablePBC=new cStencilPBC[StencilTablePBCLength];

      CMPI_channel pipe;

      pipe.init(5000000);
      pipe.openBcast(0);

      for (i=0;i<StencilTablePBCLength;i++) pipe.recv(StencilTablePBC[i],0);

      pipe.closeBcast();
    }

    //prepare data that is needed for a 'fast' performing of the data exchange operation
    for (iStencil=0;iStencil<StencilTablePBCLength;iStencil++) {
      //determine whether the current MPI processes in involved into the communication
      StencilTablePBC[iStencil].InvolvedFlag=false;

      for (i=0;i<StencilTablePBC[iStencil].InvolvedThreadTableLength;i++) if (StencilTablePBC[iStencil].InvolvedThreadTable[i]==PIC::ThisThread) {
        StencilTablePBC[iStencil].InvolvedFlag=true;
        break;
      }

      //if the current thread is involved into the communication -> initialize pointer to the associated data vector
      if (StencilTablePBC[iStencil].InvolvedFlag==true) {
        //loop through all points in the stencil
        for (i=0;i<StencilTablePBC[iStencil].nStencilPoints;i++) {
          //loop through all threads that have this point
          for (j=0;j<StencilTablePBC[iStencil].PointTable[i].ThreadTableLength;j++) {
            if (StencilTablePBC[iStencil].PointTable[i].ThreadTable[j]==PIC::ThisThread) {
              //the current MPI process has point 'j' of the stencil
              cNodeSetElement Set=StencilTablePBC[iStencil].PointTable[i].NodeSetTable[j];
              cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.findAMRnodeWithID(Set.nodeid);


              if (node->block!=NULL) {
                StencilTablePBC[iStencil].PointTable[i].AssociatedDataPointer=node->block->GetCornerNode(_getCornerNodeLocalNumber(Set.i,Set.j,Set.k))->GetAssociatedDataBufferPointer();
              }
              else {
                //the block is not allocated => the corner associated data vectors are not definied
                exit(__LINE__,__FILE__,"Error: the block is not allocated");
              }
            }
          }
        }
      }

    }
  }

    //update the coundater
    nMeshModificationCounter=PIC::Mesh::mesh.nMeshModificationCounter;
    PIC::Mesh::mesh.meshModifiedFlag_CountMeshElements=meshModifiedFlag_CountMeshElements;
    PIC::Parallel::RebalancingTime+=MPI_Wtime()-StartTime;
  }


  //performe the data exchange session
  char recvAssociatedDataBuffer[PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];

  //1. combine 'corner' node data
  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) for (iStencil=0;iStencil<StencilTableLength;iStencil++) {
    int iThread;

    if (StencilTable[iStencil].StencilLength>1) {
      //there are more that one MPI processes that contributed to the state vector of the corner node
      if (PIC::ThisThread==StencilTable[iStencil].StencilThreadTable[0]) {
        //this thread will combine data from all other MPI processes
        //recieve the data
        for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
          MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[iThread],iStencil,MPI_GLOBAL_COMMUNICATOR,&status);

          if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,recvAssociatedDataBuffer);
          }
          else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData is not defined");
        }

        //send out the state vector to MPI processes that have contributed to it
        for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
          MPI_Send(StencilTable[iStencil].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[iThread],iStencil,MPI_GLOBAL_COMMUNICATOR);
        }
      }

      else {
        for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) if (PIC::ThisThread==StencilTable[iStencil].StencilThreadTable[iThread]) {
          //this thread will contribute to the colelcted corner node associated data
          MPI_Send(StencilTable[iStencil].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR);

          //recieve the associated data
          MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,&status);

          if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
            PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,recvAssociatedDataBuffer);
          }
          else{
            memcpy(StencilTable[iStencil].AssociatedDataPointer,recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
          }

          break;
        }
      }
    }
  }

  //2. processes the 'real' domain boundary in case periodic boundary conditions are in use
  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) for (iStencil=0;iStencil<StencilTablePBCLength;iStencil++) {
     int iThread,thread,ipoint;

     if (StencilTablePBC[iStencil].InvolvedFlag==true) {
       char *AssociateDataVector=StencilTablePBC[iStencil].PointTable[0].AssociatedDataPointer;

       for (ipoint=1;ipoint<StencilTablePBC[iStencil].nStencilPoints;ipoint++) { //ipoint=1 is correct because the processing node is selected such that point=0 is already accounted for
         if (StencilTablePBC[iStencil].ProcessingThread==PIC::ThisThread) {
           //the current thread is the processing manager of the stencil
           if (StencilTablePBC[iStencil].SourceThreadTable[ipoint]==PIC::ThisThread) {
             //the data are located with the same MPI process
             if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
               PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(AssociateDataVector,StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer);
             }
             else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData is not defined");
           }
           else {
             //the data need to be recieved before processing
             MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTablePBC[iStencil].SourceThreadTable[ipoint],iStencil,MPI_GLOBAL_COMMUNICATOR,&status);

             if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
               PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(AssociateDataVector,recvAssociatedDataBuffer);
             }
             else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData is not defined");
           }

         }
         else if (StencilTablePBC[iStencil].SourceThreadTable[ipoint]==PIC::ThisThread) {
           //the MPI process is not "Processing" but serves as a source of the associated data vector for the 'ipoint' of the stencil
           MPI_Send(StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTablePBC[iStencil].ProcessingThread,iStencil,MPI_GLOBAL_COMMUNICATOR);
         }
       }

       //processing of the data is complete -> send it out to the involved MPI processes
       if (StencilTablePBC[iStencil].ProcessingThread==PIC::ThisThread) {
         int cnt=0;

         //copy processes associated data vector into a temporary buffer
         if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
           PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(recvAssociatedDataBuffer,AssociateDataVector);
         }
         else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData is not defined");


         //copy/send thecontent of the processed associated data buffer to all MPI processes that are involved in the stencil
         for (ipoint=0;ipoint<StencilTablePBC[iStencil].nStencilPoints;ipoint++) for (iThread=0;iThread<StencilTablePBC[iStencil].PointTable[ipoint].ThreadTableLength;iThread++) {
           thread=StencilTablePBC[iStencil].PointTable[ipoint].ThreadTable[iThread];

           if (thread==PIC::ThisThread) {
             //copy the processes associated data
             PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer,recvAssociatedDataBuffer);
           }
           else {
             //send processes associated data vector
             MPI_Send(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,thread,0,MPI_GLOBAL_COMMUNICATOR);
           }
         }
       }
       else {
         //recieve the processes assoviated vactor
         for (ipoint=0;ipoint<StencilTablePBC[iStencil].nStencilPoints;ipoint++) for (iThread=0;iThread<StencilTablePBC[iStencil].PointTable[ipoint].ThreadTableLength;iThread++) {
           thread=StencilTablePBC[iStencil].PointTable[ipoint].ThreadTable[iThread];

           if (thread==PIC::ThisThread) {
             //recieve the processed associated data
             MPI_Status status;

             MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTablePBC[iStencil].ProcessingThread,0,MPI_GLOBAL_COMMUNICATOR,&status);
             PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer,recvAssociatedDataBuffer);
           }
         }

       }
     }
  }
}


void PIC::Parallel::ProcessCornerBlockBoundaryNodes() {
  int thread,iThread,i,j,k,iface;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  PIC::Mesh::cDataCornerNode *CornerNode;
  char *CornerNodeAssociatedData;
  PIC::Mesh::cDataBlockAMR *block;
  MPI_Status status;

  if (CornerBlockBoundaryNodes::ActiveFlag==false) return;

  
  const int iFaceMin[6]={0,_BLOCK_CELLS_X_,0,0,0,0};
  const int iFaceMax[6]={0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_};

  const int jFaceMin[6]={0,0,0,_BLOCK_CELLS_Y_,0,0};
  const int jFaceMax[6]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_};

  const int kFaceMin[6]={0,0,0,0,0,_BLOCK_CELLS_Z_};
  const int kFaceMax[6]={_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,0,_BLOCK_CELLS_Z_};
  
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

  localMeshChangeFlag=(nMeshModificationCounter==PIC::Mesh::mesh.nMeshModificationCounter) ? 0 : 1;
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
    int NewTableLength=0;
    int iNode,jNode,kNode,FlagSum;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* ActualNode;
    bool meshModifiedFlag_CountMeshElements=PIC::Mesh::mesh.meshModifiedFlag_CountMeshElements;
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
    //    for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
      int flag;
      double eps = 0.3*PIC::Mesh::mesh.EPS;
      double dx[3];
      int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
      for (int idim=0; idim<3; idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];
        
      if (node->Thread==PIC::ThisThread) { 
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {

                double x[3];
                int ind[3]={i,j,k};
                bool outBoundary=false;// atBoundary=false;
                for (int idim=0; idim<3; idim++) {
                  x[idim]=node->xmin[idim]+dx[idim]*ind[idim];
#if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_                   
                  if ((x[idim]-PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim])<-eps || (x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim])>eps){
                    outBoundary=true;
                    break;
                  }
#endif                    
                }
                  
                if (outBoundary==false) nPointsToCommThisThread++;
        
              }// for (iface=0;iface<6;iface++)  for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)
      }//if (node->Thread==PIC::ThisThread)
    }//for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)

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
    //    for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
      int flag;
      
      if (node->Thread==PIC::ThisThread) { 
        double dx[3];
        int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
        double eps = 0.3*PIC::Mesh::mesh.EPS;
        for (int idim=0; idim<3; idim++) dx[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];
        
        for (iface=0;iface<6;iface++)  for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {          
                
                double x[3];
                int ind[3]={i,j,k};
                bool IsBoundaryNode=false, outBoundary=false;
                // if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_){ 
                for (int idim=0; idim<3; idim++) {
                  x[idim]=node->xmin[idim]+dx[idim]*ind[idim];
#if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_   
                  if ((x[idim]-PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim])<-eps || (x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim])>eps) {
                    outBoundary=true;
                    break;
                  }                  
#endif
                }
                    

                if (!outBoundary){
                  
                  for (int idim=0; idim<3; idim++)
                    x[idim]=node->xmin[idim]+dx[idim]*ind[idim];
                                                      
                  for (int idim=0; idim<3; idim++)  {
                                
#if _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_
                    if (fabs(x[idim]-PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[idim])<eps)
                      x[idim] = PIC::BC::ExternalBoundary::Periodic::xminOriginal[idim];
#endif
                    tempPnt[idim]=x[idim];
                  }
                  
                  *tempDataBuffPnt = node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer(); 
                  tempDataBuffPnt++;
                  tempPnt += 3;
                }


              }// for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)
      }//if (node->Thread==PIC::ThisThread)
    }//for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)
    
    
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
      if (flag!=0){
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
      
      if (counter[0]<PIC::nTotalThreads-1)
        MPI_Testany(PIC::nTotalThreads-1, RecvCoordList, p, flag, MPI_STATUS_IGNORE);
      if (counter[1]<PIC::nTotalThreads-1)
        MPI_Testany(PIC::nTotalThreads-1, RecvDataBufferPtrList, p+1, flag+1, MPI_STATUS_IGNORE);
      
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
    //Tree.printTree();

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
    nMeshModificationCounter=PIC::Mesh::mesh.nMeshModificationCounter;
    PIC::Mesh::mesh.meshModifiedFlag_CountMeshElements=meshModifiedFlag_CountMeshElements;
    PIC::Parallel::RebalancingTime+=MPI_Wtime()-StartTime;
  }


  //performe the data exchange session
  char recvAssociatedDataBuffer[PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];

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

  if (globalMeshChangeFlag!=0 && PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL){

    if (LoadList!=NULL) delete [] LoadList;
    LoadList=new int [PIC::Mesh::mesh.nTotalThreads];
    for (int i=0;i<PIC::Mesh::mesh.nTotalThreads;i++) LoadList[i]=0;
    
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
    for (int i=0;i<PIC::Mesh::mesh.nTotalThreads;i++) {
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
    processRecvDataBuffer[0] = new char [sumProcessBufferNumber*PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];
    for (int i=1;i<sumProcessBufferNumber;i++) {
      processRecvDataBuffer[i]=i*PIC::Mesh::cDataCornerNode::totalAssociatedDataLength+processRecvDataBuffer[0];
    }

    if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {

      if (copyRecvDataBuffer!=NULL){
        delete [] copyRecvDataBuffer[0];
        delete [] copyRecvDataBuffer;
      }
      
      copyRecvDataBuffer=new char * [CopyDataBufferDest.size()];
      copyRecvDataBuffer[0]=new char [CopyDataBufferDest.size()*PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];
      for (int i=1;i<CopyDataBufferDest.size();i++) copyRecvDataBuffer[i]=copyRecvDataBuffer[0]+i*PIC::Mesh::cDataCornerNode::totalAssociatedDataLength;
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
  
  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {

    //clear counter
    bool isTestOn=false;
    
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
            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,*it);
          }
        
          for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
            MPI_Irecv(processRecvDataBuffer[iProcessBuffer],PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[iThread],iStencil,MPI_GLOBAL_COMMUNICATOR,processRecvList+iProcessBuffer);
            iProcessBuffer++;
          }
          
        }else{
          //this thread will be copyThread, send data to the processThread, 
          //and recv the processed data from the processThread
          std::vector<char *>::iterator it;
          for (it=StencilTable[iStencil].localDataPntArr.begin();it!=StencilTable[iStencil].localDataPntArr.end();it++){
            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,*it);                  
          }
          MPI_Isend(StencilTable[iStencil].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,processSendList+iCopyBuffer);
            
          //start the recv data request
          if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
            MPI_Irecv(copyRecvDataBuffer[iCopyBuffer],PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,copyRecvList+iCopyBuffer);
          }else{
            MPI_Irecv(CopyDataBufferDest[iCopyBuffer],PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,copyRecvList+iCopyBuffer);
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
            
            //printf("thread id:%d,bufferid:%d processed, nProcessBufferDone:%d/total buffer:%d,\n",PIC::ThisThread,q,nProcessBufferDone,sumProcessBufferNumber);

            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[wholeStencilIndex].AssociatedDataPointer,processRecvDataBuffer[q]);
          
            nProcessBufferDone++;
            StencilFinishedBufferNumber[processStencilIndex]++;
            //if  finished one stencil processing, send it back
            if (StencilFinishedBufferNumber[processStencilIndex]==ProcessStencilBufferNumber[processStencilIndex]){
              for (iThread=1;iThread<StencilTable[wholeStencilIndex].StencilLength;iThread++) {        
                MPI_Isend(StencilTable[wholeStencilIndex].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[wholeStencilIndex].StencilThreadTable[iThread],wholeStencilIndex,MPI_GLOBAL_COMMUNICATOR,copySendList+iProcessSendBuffer);
                iProcessSendBuffer++;
              }
              nProcessStencilDone++;
            }
            
          }
        }//for (int i=nProcessBufferDone; i<iProcessBuffer;i++)

        for (int i=nCopyBufferDone; i<iCopyBuffer;i++) {
          int q, flag;

          MPI_Testany(iCopyBuffer, copyRecvList, &q, &flag, MPI_STATUS_IGNORE);
          
          if (flag!=0 && q!=MPI_UNDEFINED){
            //int wholeStencilIndex=CopyDataBufferStencilIndex[q];
            if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
              PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(CopyDataBufferDest[q],copyRecvDataBuffer[q]);
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
            
          // printf("thread id:%d,bufferid:%d processed, nProcessBufferDone:%d/total buffer:%d,\n",PIC::ThisThread,q,nProcessBufferDone,sumProcessBufferNumber);
          /*
            if (wholeStencilIndex==714){
            printf("ist:%d,in processing thread id: %d, size:%d, (while) before Jx:%e\n",wholeStencilIndex,PIC::ThisThread,StencilTable[wholeStencilIndex].localDataPntArr.size(), ((double *)(StencilTable[wholeStencilIndex].AssociatedDataPointer+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))[6]);
            }
          */
          PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[wholeStencilIndex].AssociatedDataPointer,processRecvDataBuffer[q]);
          
          nProcessBufferDone++;
          StencilFinishedBufferNumber[processStencilIndex]++;
                
          if (StencilFinishedBufferNumber[processStencilIndex]==ProcessStencilBufferNumber[processStencilIndex]){
            for (iThread=1;iThread<StencilTable[wholeStencilIndex].StencilLength;iThread++) {        
              MPI_Isend(StencilTable[wholeStencilIndex].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[wholeStencilIndex].StencilThreadTable[iThread],wholeStencilIndex,MPI_GLOBAL_COMMUNICATOR,copySendList+iProcessSendBuffer);
              iProcessSendBuffer++;
            }
            nProcessStencilDone++;
          }
        }
      }//for (int i=nProcessBufferDone; i<iProcessBuffer;i++)
      
      for (int i=nCopyBufferDone; i<iCopyBuffer;i++) {
        int q, flag;
        
        MPI_Testany(iCopyBuffer, copyRecvList, &q, &flag, MPI_STATUS_IGNORE);
        
        if (flag!=0 && q!=MPI_UNDEFINED){
          int wholeStencilIndex=CopyDataBufferStencilIndex[q];
          if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
            PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(CopyDataBufferDest[q],copyRecvDataBuffer[q]);
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
        if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
          for (it=StencilTable[iSt].localDataPntArr.begin();it!=StencilTable[iSt].localDataPntArr.end();it++){
            PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(*it,StencilTable[iSt].AssociatedDataPointer);  
          }
        }else{
          for (it=StencilTable[iSt].localDataPntArr.begin();it!=StencilTable[iSt].localDataPntArr.end();it++)
            memcpy(*it,StencilTable[iSt].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
        }
      }
    }
    
    }//if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL)
  
}



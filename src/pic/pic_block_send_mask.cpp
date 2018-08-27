
//function for generation Center/Corner node masks used in the data syncronization between MPI processes


/*
 * pic_block_send_mask.cpp
 *
 *  Created on: Aug 24, 2018
 *      Author: vtenishe
 */

#include "pic.h"

void PIC::Mesh::BlockElementSendMask::Set(bool flag,unsigned char* CenterNodeMask,unsigned char* CornerNodeMask) {
  unsigned char m=(flag==true) ? 0xff : 0;
  int i,imax;

  //1. Set the center node mask
  for (i=0,imax=CenterNode::GetSize();i<imax;i++) CenterNodeMask[i]=m;

  //2. St the corner node mask
  for (i=0,imax=CornerNode::GetSize();i<imax;i++) CornerNodeMask[i]=m;
}

int PIC::Mesh::BlockElementSendMask::CornerNode::GetSize() {return 1+((_TOTAL_BLOCK_CELLS_X_+1)*(_TOTAL_BLOCK_CELLS_Y_+1)*(_TOTAL_BLOCK_CELLS_Z_+1))/8;}
int PIC::Mesh::BlockElementSendMask::CenterNode::GetSize() {return 1+(_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_)/8;}

//Init Mask for the Layer Block
void PIC::Mesh::BlockElementSendMask::InitLayerBlockBasic(cTreeNodeAMR<cDataBlockAMR>* Node,int To,unsigned char* CenterNodeMask,unsigned char* CornerNodeMask) {
  int i,j,k;
  cTreeNodeAMR<cDataBlockAMR>* neibNode;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  //set the default value
  Set(false,CenterNodeMask,CornerNodeMask);

  for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
    CenterNode::Set(true,i,j,k,CenterNodeMask);
    CornerNode::Set(true,i,j,k,CornerNodeMask);
  }

  //Loop through the "right" boundary of the block
  int iface,iFaceTable[3]={1,3,5};

   for (int ipass=0;ipass<3;ipass++) {
    iface=iFaceTable[ipass];

    if ((neibNode=Node->GetNeibFace(iface,0,0))!=NULL) if (neibNode->RefinmentLevel<Node->RefinmentLevel) {
      //the current block has more points than the neibour -> need to send the point that exist in the current block but not exist in the neib block

      switch (iface) {
      case 1:
        //the plane normal to 'x' and at the maximum 'x'
        for (i=iCellMax,k=1;k<kCellMax+1;k+=2) for (j=1;j<jCellMax+1;j+=2) {
          CornerNode::Set(true,i,j,k,CornerNodeMask);
        }
        break;

      case 3:
        //the plane is normal to the 'y' direction, and is at maximum 'y'
        for (j=jCellMax,k=1;k<kCellMax+1;k+=2) for (i=1;i<iCellMax+1;i+=2) {
          CornerNode::Set(true,i,j,k,CornerNodeMask);
        }
        break;

      case 5:
        //the plane normal to 'z' and at the maximum 'z'
        for (k=kCellMax,j=0;j<jCellMax;j+=2) for (i=0;i<iCellMax;i+=2) {
          CornerNode::Set(true,i,j,k,CornerNodeMask);
        }
        break;
      }
    }
  }
}

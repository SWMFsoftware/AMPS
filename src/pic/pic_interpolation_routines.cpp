//$Id$
//interpolation routines

/*
 * pic_interpolation_routines.cpp
 *
 *  Created on: Jul 18, 2015
 *      Author: vtenishe
 */

#include "pic.h"

PIC::InterpolationRoutines::CellCentered::cStencil _TARGET_DEVICE_*PIC::InterpolationRoutines::CellCentered::StencilTable=NULL;
PIC::InterpolationRoutines::CornerBased::cStencil _TARGET_DEVICE_ *PIC::InterpolationRoutines::CornerBased::StencilTable=NULL;

int PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent=0;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::nBlockFoundMax];
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::last=NULL;

//the locally ordered interpolation coeffcients for the corner based interpolation procedure
thread_local double PIC::InterpolationRoutines::CornerBased::InterpolationCoefficientTable_LocalNodeOrder[8]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

// macro switch is needed in the case some other interpolation is used
// and interface function is not compiled
#if _PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_
//definition of the interface functions that are used to access the FORTRAN written linear cell cenetered interpolation library
extern "C"{
  void interface__cell_centered_linear_interpolation__init_stencil_(int* nDim, double* XyzIn_D, int* nIndexes, int* nCell_D,  int* nGridOut, double* Weight_I, int* iIndexes_II, int* IsSecondOrder, int* UseGhostCell);
}
#endif//_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_


//initialize the interpolation module
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::Init() {

  //init the stencil table
  CellCentered::StencilTable=new CellCentered::cStencil[PIC::nTotalThreadsOpenMP];
  CornerBased::StencilTable=new CornerBased::cStencil[PIC::nTotalThreadsOpenMP];
}

//determine stencil for the cell centered piecewise constant interpolation
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil) {
  int i,j,k;
  int nd;
  PIC::Mesh::cDataCenterNode *cell;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  //flush the stencil
  Stencil.flush();

  //find the block
  if (node==NULL) {
    node=PIC::Mesh::mesh->findTreeNode(x);

    if (node==NULL) exit(__LINE__,__FILE__,"Error: the location is outside of the computational domain");
    if (node->block==NULL) {
      char msg[200];

      #ifndef __CUDA_ARCH__
      sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
      exit(__LINE__,__FILE__,msg);
      #else
      printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
      exit(__LINE__,__FILE__);
      #endif
    }
  }
  else if (node->block==NULL) {
    char msg[200];

    #ifndef __CUDA_ARCH__
    sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
    exit(__LINE__,__FILE__,msg);
    #else 
    printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
    exit(__LINE__,__FILE__,msg); 
    #endif
  }

  //find cell
  nd = PIC::Mesh::mesh->fingCellIndex(x,i,j,k,node,false);
  cell=node->block->GetCenterNode(nd);//getCenterNodeLocalNumber(i,j,k));

  //add the cell to the stencil
  if (cell!=NULL) {
    Stencil.AddCell(1.0,cell,nd);
  }
  else exit(__LINE__,__FILE__,"Error: cell is not initialized");
}

//determine the stencil for the cell centered linear interpolation using interpolation library from ../share/Library/src/
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(double *XyzIn_D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil) {

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();

  if  (_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__SWMF_) {
    exit(__LINE__,__FILE__,"Error: the procedure is not adapted fro using with OpenMP: INTERFACE::BlockFound... need to be converted into an array as PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP] ");
  }
  #else
  int ThreadOpenMP=0;
  #endif

  Stencil.flush();

  // macro switch is needed in the case some other interpolation is used
  // and interface function is not compiled
  if (_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
    //find the block if needed
    if (node==NULL) {
      node=PIC::Mesh::mesh->findTreeNode(XyzIn_D);

      if (node==NULL) exit(__LINE__,__FILE__,"Error: the location is outside of the computational domain");
      if (node->block==NULL) {
        char msg[200];

        #ifndef __CUDA_ARCH__
        sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
        exit(__LINE__,__FILE__,msg);
        #else 
        printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
        exit(__LINE__,__FILE__);
        #endif
      }
    }
    else if (node->block==NULL) {
      char msg[200];

      #ifndef __CUDA_ARCH__
      sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
      exit(__LINE__,__FILE__,msg);
      #else 
      printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
      exit(__LINE__,__FILE__,msg);
      #endif
    }

    //check whether the point is located deep in the block -> use three linear interpolation
    double iLoc,jLoc,kLoc;
    double *xmin,*xmax;

    xmin=node->xmin;
    xmax=node->xmax;

    iLoc=(XyzIn_D[0]-xmin[0])/(xmax[0]-xmin[0])*_BLOCK_CELLS_X_;
    jLoc=(XyzIn_D[1]-xmin[1])/(xmax[1]-xmin[1])*_BLOCK_CELLS_Y_;
    kLoc=(XyzIn_D[2]-xmin[2])/(xmax[2]-xmin[2])*_BLOCK_CELLS_Z_;

    #if _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
    //in case the mesh is uniform (_AMR_MESH_TYPE_ == _AMR_MESH_TYPE__UNIFORM_) -> use simple trilinear interpolation
    if (_AMR_MESH_TYPE_ == _AMR_MESH_TYPE__UNIFORM_) {
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil);
      return;
    }

    //if the point of interest is deep inside the block or all neighbors of the block has the same resolution level -> use simple trilinear interpolation
    if ((node->RefinmentLevel==node->minNeibRefinmentLevel)&&(node->RefinmentLevel==node->maxNeibRefinmentLevel)) {
      //all blocks around has the same level of the resolution
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil);
      return;
    }
    else if ((1.0<iLoc)&&(iLoc<_BLOCK_CELLS_X_-1) && (1.0<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-1) && (1.0<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-1)) {
      //the point is deep inside the block
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil);
      return;
    }
    else if (node->RefinmentLevel==node->minNeibRefinmentLevel) {
      if  ((0.5<iLoc)&&(iLoc<_BLOCK_CELLS_X_-0.5) && (0.5<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-0.5) && (0.5<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-0.5)) {
        //1. the point is deeper than a half cell size insode the block
        //2. the block is geometrically largest among the neibours
        //3. -> it is only the intermal points of the block that will be used in the stencil
        GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil);
        return;
      }
      else  {
        //the block is largest between the neibours
        //getermine the size limit for the interpolation stencil
        double xStencilMin[3],xStencilMax[3];
        double dxCell[3];

        dxCell[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
        dxCell[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
        dxCell[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

        #pragma ivdep
        for (int idim=0;idim<3;idim++) {
          int iInterval;

          iInterval=(int)((XyzIn_D[idim]-(xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
          xStencilMin[idim]=(xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
          xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
        }

        GetTriliniarInterpolationMutiBlockStencil(XyzIn_D,xStencilMin,xStencilMax,node,Stencil);
        return;
      }
    }
    else {
      //1. find a coarser block that is close to the point, and can be used for interpolation
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *CoarserBlock=NULL;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *NeibNode;
      int idim,iFace;
      double dxCell[3];

      dxCell[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
      dxCell[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
      dxCell[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

      int nBlockCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
      double dmin=10.0*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;

      bool CornerTestFlagTable[8]={false,false,false,false,  false,false,false,false};
      bool EdgeTestFlagTable[12]={false,false,false,false, false,false,false,false, false,false,false,false};

      for (idim=0;idim<3;idim++) {
        switch (idim) {
        case 0:
          if (iLoc<=1.0) iFace=0;
          else if (iLoc>=_BLOCK_CELLS_X_-1.0) iFace=1;
          else continue;

          break;
        case 1:
          if (jLoc<=1.0) iFace=2;
          else if (jLoc>=_BLOCK_CELLS_Y_-1.0) iFace=3;
          else continue;

          break;
        case 2:
          if (kLoc<=1.0) iFace=4;
          else if (kLoc>=_BLOCK_CELLS_Z_-1.0) iFace=5;
          else continue;
        }

        //check blocks connected through the face
        NeibNode=node->GetNeibFace(iFace,0,0,PIC::Mesh::mesh);

        double xLoc[3]={iLoc,jLoc,kLoc};

        if (NeibNode!=NULL) if ((NeibNode->RefinmentLevel<node->RefinmentLevel)&&(NeibNode->IsUsedInCalculationFlag==true)) {
          //found a coarser block
          int cnt=0;

          for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

          if (cnt==3) {
            //the block can be used for interpolation
            switch(iFace) {
            case 0:
              if ((xLoc[0]<1.0)&&(xLoc[0]<dmin)) dmin=xLoc[0],CoarserBlock=NeibNode;
              break;

            case 1:
              if ((xLoc[0]>nBlockCells[0]-1)&&(nBlockCells[0]-xLoc[0]<dmin)) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
              break;

            case 2:
              if ((xLoc[1]<1.0)&&(xLoc[1]<dmin)) dmin=xLoc[1],CoarserBlock=NeibNode;
              break;

            case 3:
              if ((xLoc[1]>nBlockCells[1]-1)&&(nBlockCells[1]-xLoc[1]<dmin)) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
              break;

            case 4:
              if ((xLoc[2]<1.0)&&(xLoc[2]<dmin)) dmin=xLoc[2],CoarserBlock=NeibNode;
              break;

            case 5:
              if ((xLoc[2]>nBlockCells[2]-1)&&(nBlockCells[2]-xLoc[2]<dmin)) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
              break;

            default:
              exit(__LINE__,__FILE__,"error: something is wrong");
            }
          }
        }

        //check blocks connected though the edges
        static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};


        for (int iEdge=0;iEdge<4;iEdge++) if (EdgeTestFlagTable[faceEdges[iFace][iEdge]]==false) {
          EdgeTestFlagTable[faceEdges[iFace][iEdge]]=true;
          NeibNode=node->GetNeibEdge(faceEdges[iFace][iEdge],0,PIC::Mesh::mesh);

          if (NeibNode!=NULL) if ((NeibNode->RefinmentLevel<node->RefinmentLevel)&&(NeibNode->IsUsedInCalculationFlag==true)) {
            //found a coarser block
            int cnt=0;

            for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

            if (cnt==3) {
              //the block can be used for interpolation
              switch (faceEdges[iFace][iEdge]) {
              case 0:
                if ((xLoc[1]<1.0)&&(xLoc[2]<1.0)) {
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 1:
                if ((xLoc[1]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 2:
                if ((xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 3:
                if ((xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;


              case 4:
                if ((xLoc[0]<1.0)&&(xLoc[2]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 5:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 6:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 7:
                if ((xLoc[0]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 8:
                if ((xLoc[0]<1.0)&&(xLoc[1]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                }

                break;

              case 9:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                }

                break;

              case 10:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]>nBlockCells[1]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                }

                break;

              case 11:
                if ((xLoc[0]<1.0)&&(xLoc[1]>nBlockCells[1]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                }


                break;

              default:
                exit(__LINE__,__FILE__,"error: something is wrong");
              }
            }
          }
        }

        //check blocks connected through the corners
        static const int FaceNodeMap[6][4]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};

        for (int iCorner=0;iCorner<4;iCorner++) if (CornerTestFlagTable[FaceNodeMap[iFace][iCorner]]==false) {
          CornerTestFlagTable[FaceNodeMap[iFace][iCorner]]=true;
          NeibNode=node->GetNeibCorner(FaceNodeMap[iFace][iCorner],PIC::Mesh::mesh);

          if (NeibNode!=NULL) if ((NeibNode->RefinmentLevel<node->RefinmentLevel)&&(NeibNode->IsUsedInCalculationFlag==true)) {
            //found a coarser block
            int cnt=0;

            for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

            if (cnt==3) {
              //the block can be used for interpolation
              switch (FaceNodeMap[iFace][iCorner]) {
              case 0:
                if ((xLoc[0]<1.0)&&(xLoc[1]<1.0)&&(xLoc[2]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 1:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]<1.0)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 2:
                if ((xLoc[0]<1.0)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 3:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;


              case 4:
                if ((xLoc[0]<1.0)&&(xLoc[1]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 5:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 6:
                if ((xLoc[0]<1.0)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 7:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              default:
                exit(__LINE__,__FILE__,"error: something is wrong");
              }
            }
          }
        }
      }

      if (CoarserBlock!=NULL) {
        //getermine the size limit for the interpolation stencil
        double xStencilMin[3],xStencilMax[3];
        double dxCell[3];

        dxCell[0]=(CoarserBlock->xmax[0]-CoarserBlock->xmin[0])/_BLOCK_CELLS_X_;
        dxCell[1]=(CoarserBlock->xmax[1]-CoarserBlock->xmin[1])/_BLOCK_CELLS_Y_;
        dxCell[2]=(CoarserBlock->xmax[2]-CoarserBlock->xmin[2])/_BLOCK_CELLS_Z_;

        #pragma ivdep
        for (idim=0;idim<3;idim++) {
          int iInterval;

          iInterval=(int)((XyzIn_D[idim]-(CoarserBlock->xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
          xStencilMin[idim]=(CoarserBlock->xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
          xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
        }

        GetTriliniarInterpolationMutiBlockStencil(XyzIn_D,xStencilMin,xStencilMax,CoarserBlock,Stencil);


        /* the following is used to get continuous interpolsation
        1. if the point of interpolation is in the "coarse" block -> use GetTriliniarInterpolationMutiBlockStencil()
        2. if the point of interpolation is in the "fine" block ->
           a) def "dmin" is the min length to the boundary of the block
           b) if dmin<0.5 -> "coarse" stencil interpolation
           c) if 0.5<dmin<1-> combination of "fine" and "coarse" stencils:
             interpolation=(dmin-0.5)/0.5*fine stencil+[1-(dmin-0.5)/0.5)]*corse stencil
           d) is 1<dmin -> "fine" stencil interpolation
        */

        if ((0.5<dmin)&&(dmin<=1.0)) {
          PIC::InterpolationRoutines::CellCentered::cStencil FineStencil;

          GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,FineStencil);

          Stencil.MultiplyScalar(1.0-(dmin-0.5)/0.5);
          FineStencil.MultiplyScalar((dmin-0.5)/0.5);

          Stencil.Add(&FineStencil);
        }

        return;
      }

      //2.threre is no a coarser block that can be used for interpolation
      //The current block is the coarstest that can be used to construct the interpolation stencil
      double xStencilMin[3],xStencilMax[3];

      dxCell[0]=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
      dxCell[1]=(node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_;
      dxCell[2]=(node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_;

      #pragma ivdep
      for (idim=0;idim<3;idim++) {
        int iInterval;

        iInterval=(int)((XyzIn_D[idim]-(node->xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
        xStencilMin[idim]=(node->xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
        xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
      }

      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil);
      return;
    }

    #elif _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__SWMF_
    if ( ((node->RefinmentLevel==node->minNeibRefinmentLevel)&&(node->RefinmentLevel==node->maxNeibRefinmentLevel)) ||
        (0.5<iLoc)&&(iLoc<_BLOCK_CELLS_X_-0.5) && (0.5<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-0.5) && (0.5<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-0.5) ) {
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil);
      return;
    }

    // if the point of interest is very close to the cell center
    // then truncated stencil to the single point
    if( fabs(iLoc-0.5-(long)iLoc) < PrecisionCellCenter &&
        fabs(jLoc-0.5-(long)jLoc) < PrecisionCellCenter &&
        fabs(kLoc-0.5-(long)kLoc) < PrecisionCellCenter   ){
          PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node,Stencil); 
          return;
    }

    //re-init variables in the INTERFACE and flush the Stencil
    INTERFACE::iBlockFoundCurrent=0;
    for(int iBlock = 0; iBlock < INTERFACE::nBlockFoundMax; iBlock++ ) INTERFACE::BlockFound[iBlock] = NULL;

    Stencil.flush();
    INTERFACE::last=node;

    //prepare the call of FORTRAN interpolate_amr subroutine
    int nDim       = DIM;
    // number of indices to identify cell ON A GIVEN PROCESSOR:
    // block id + DIM indices of cell in a block
    int nIndexes   = DIM + 1;
    int nCell_D[3] = {_BLOCK_CELLS_X_, _BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

    // IMPORTANT NOTE: iIndexes_II has an extra index per cell which is
    // a processor ID => need one more index in addition to nIndexes
    int iIndexes_II[(DIM+1+1)*PIC::InterpolationRoutines::nMaxStencilLength];
    double WeightStencil[PIC::InterpolationRoutines::nMaxStencilLength];
    int IsSecondOrder;
    int UseGhostCell=1;
    int nCellStencil;

    // call the interpolation subroutine
    interface__cell_centered_linear_interpolation__init_stencil_(&nDim,XyzIn_D,&nIndexes, nCell_D, &nCellStencil, WeightStencil, iIndexes_II, &IsSecondOrder, &UseGhostCell);

    // size of the stencil and weights are known
    // need to identify cells in the stencil
    // NOTE: FORTRAN is column-major for 2D arrays like iIndexes_II
    for (int iCellStencil = 0; iCellStencil < nCellStencil; iCellStencil++) {
      int ind[3]={0,0,0};
      PIC::Mesh::cDataBlockAMR  *block;
      PIC::Mesh::cDataCenterNode *cell;

  //    int iThread = iIndexes_II[0    +iCellStencil*(nIndexes+1)];
      for(int i = 0; i < DIM; i++)
        //cell indices are 1-based in FORTRAN
        ind[i]    = iIndexes_II[1+i  +iCellStencil*(nIndexes+1)]-1;
      int iBlock  = iIndexes_II[1+DIM+iCellStencil*(nIndexes+1)];

      //retrieve the pointer to the current cell
      if ((block=INTERFACE::BlockFound[iBlock]->block)==NULL) {
        PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node,Stencil);
        return;
      }

      cell=block->GetCenterNode(_getCenterNodeLocalNumber(ind[0],ind[1],ind[2]));

      if (cell==NULL) {
        //there is no enough information to build up linear interpolation -> use constant interpolation
        PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node,Stencil);
        return;
      }

      int nd=_getCenterNodeLocalNumber(ind[0],ind[1],ind[2]);
      Stencil.AddCell(WeightStencil[iCellStencil],INTERFACE::BlockFound[iBlock]->block->GetCenterNode(nd),nd);
    }
    #else  //_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
    exit(__LINE__,__FILE__,"Error: the option is unknown");
    #endif //_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
  }
  else exit(__LINE__,__FILE__,"ERROR: cell centered linear interpolation is currently available only through interface, add corresponding block to the input file!");
}

//triliniar interpolation used inside blocks
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Linear::GetTriliniarInterpolationStencil(double iLoc,double jLoc,double kLoc,double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil &Stencil) {
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;

  block=node->block;
  if (block==NULL) exit(__LINE__,__FILE__,"Error: the block is node allocated");

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  //flush the stencil
  Stencil.flush();

  //determine the aray of the cell's pointer that will be used in the interpolation stencil
  double w[3],InterpolationWeight,totalInterpolationWeight=0.0;
  int i,j,k,i0,j0,k0,nd;

  i0=(iLoc<0.5) ? -1 : (int)(iLoc-0.50);
  j0=(jLoc<0.5) ? -1 : (int)(jLoc-0.50);
  k0=(kLoc<0.5) ? -1 : (int)(kLoc-0.50);

  //get coefficients used in determening of the interpolation weights
  w[0]=iLoc-(i0+0.5);
  w[1]=jLoc-(j0+0.5);
  w[2]=kLoc-(k0+0.5);

  for (i=0;i<2;i++) for (j=0;j<2;j++) for (k=0;k<2;k++) {
    nd=_getCenterNodeLocalNumber(i0+i,j0+j,k0+k);
    cell=block->GetCenterNode(nd);

    if (cell!=NULL) {
      switch (i+2*j+4*k) {
      case 0+0*2+0*4:
        InterpolationWeight=(1.0-w[0])*(1.0-w[1])*(1.0-w[2]);
        break;
      case 1+0*2+0*4:
        InterpolationWeight=w[0]*(1.0-w[1])*(1.0-w[2]);
        break;
      case 0+1*2+0*4:
        InterpolationWeight=(1.0-w[0])*w[1]*(1.0-w[2]);
        break;
      case 1+1*2+0*4:
        InterpolationWeight=w[0]*w[1]*(1.0-w[2]);
        break;

      case 0+0*2+1*4:
        InterpolationWeight=(1.0-w[0])*(1.0-w[1])*w[2];
        break;
      case 1+0*2+1*4:
        InterpolationWeight=w[0]*(1.0-w[1])*w[2];
        break;
      case 0+1*2+1*4:
        InterpolationWeight=(1.0-w[0])*w[1]*w[2];
        break;
      case 1+1*2+1*4:
        InterpolationWeight=w[0]*w[1]*w[2];
        break;

      default:
        exit(__LINE__,__FILE__,"Error: the option is not defined");
      }

      Stencil.AddCell(InterpolationWeight,cell,nd);
      totalInterpolationWeight+=InterpolationWeight;
    }
  }

  if (Stencil.Length==0) {
    //no cell have been found -> use a canstarnt interpolation stencil
    PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node,Stencil);
    return;
  }
  else if (StencilTable->Length!=8) {
    //the interpolated stencil containes less that 8 elements -> the interpolation weights have to be renormalized
    #pragma ivdep
    for (int i=0;i<Stencil.Length;i++) {
      Stencil.Weight[i]/=totalInterpolationWeight;
    }
  }
}



_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Linear::GetTriliniarInterpolationMutiBlockStencil(double *x,double *xStencilMin,double *xStencilMax,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil &Stencil) {
  int nd;
  PIC::Mesh::cDataCenterNode *cell;

  Stencil.flush();

  class cMultiBlockStencil {
  public:
    int ijk[3];
    double x[3];
    PIC::InterpolationRoutines::CellCentered::cStencil Stencil;
  };

  cMultiBlockStencil StencilTable[2][2][2];

  int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
  double dxCell[3];
  int idim;


  for (idim=0;idim<3;idim++) dxCell[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];

  //determine coordinates of the stencil nodes
  int di,dj,dk,ijkStencilMin[3];

  for (idim=0;idim<3;idim++) {
    ijkStencilMin[idim]=(x[idim]-node->xmin[idim]<0.5*dxCell[idim]) ? -1 : (int)((x[idim]-node->xmin[idim]-0.5*dxCell[idim])/dxCell[idim]);
  }

  for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) {
    StencilTable[di][dj][dk].ijk[0]=ijkStencilMin[0]+di;
    StencilTable[di][dj][dk].ijk[1]=ijkStencilMin[1]+dj;
    StencilTable[di][dj][dk].ijk[2]=ijkStencilMin[2]+dk;

    for (idim=0;idim<3;idim++) {
      StencilTable[di][dj][dk].x[idim]=node->xmin[idim]+(StencilTable[di][dj][dk].ijk[idim]+0.5)*dxCell[idim];
    }

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *StencilNode=PIC::Mesh::mesh->findTreeNode(StencilTable[di][dj][dk].x,node);
    bool return_const_stencil=false;

    if (StencilNode==NULL) return_const_stencil=true;
    else if (StencilNode->IsUsedInCalculationFlag==false) return_const_stencil=true;

    if (return_const_stencil==true) {
      node=PIC::Mesh::mesh->findTreeNode(x,node);
      PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node,Stencil);
      return;
    }

    //init the corner of the stencil
    if (StencilNode->RefinmentLevel==node->RefinmentLevel) {
      //both blocks have the save refinment levels
      nd=_getCenterNodeLocalNumber(StencilTable[di][dj][dk].ijk[0],StencilTable[di][dj][dk].ijk[1],StencilTable[di][dj][dk].ijk[2]);
      cell=node->block->GetCenterNode(nd);//getCenterNodeLocalNumber(i,j,k));

      if (cell!=NULL) {
        StencilTable[di][dj][dk].Stencil.AddCell(1.0,cell,nd);
      }
      else {
        node=PIC::Mesh::mesh->findTreeNode(x,node);
        PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node,Stencil);
        return;
      }
    }
    else {
     //the block has a higher resolution level
     int iNeib[3];

     for (idim=0;idim<3;idim++) iNeib[idim]=2*((int)((StencilTable[di][dj][dk].x[idim]-StencilNode->xmin[idim])/dxCell[idim]));

     for (int ii=0;ii<2;ii++) for (int jj=0;jj<2;jj++) for (int kk=0;kk<2;kk++) {
       nd=_getCenterNodeLocalNumber(iNeib[0]+ii,iNeib[1]+jj,iNeib[2]+kk);
       cell=StencilNode->block->GetCenterNode(nd);

       if (cell!=NULL) {
         StencilTable[di][dj][dk].Stencil.AddCell(1.0/8.0,cell,nd);
       }
       else {
         node=PIC::Mesh::mesh->findTreeNode(x,node);
         PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node,Stencil);
         return;
       }
     }
    }
  }

  //calcualte local coordinates
  double xLoc[3];

  for (idim=0;idim<3;idim++) {
    xLoc[idim]=(x[idim]-StencilTable[0][0][0].x[idim])/dxCell[idim];
    if (xLoc[idim]<0.0) xLoc[idim]=0.0;
    if (xLoc[idim]>1.0) xLoc[idim]=1.0;
  }

  //determine the interpolation coefficients
  for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) {
    double StencilElementWeight;

    switch (di+2*dj+4*dk) {
    case 0+0*2+0*4:
      StencilElementWeight=(1.0-xLoc[0])*(1.0-xLoc[1])*(1.0-xLoc[2]);
      break;
    case 1+0*2+0*4:
      StencilElementWeight=xLoc[0]*(1.0-xLoc[1])*(1.0-xLoc[2]);
      break;
    case 0+1*2+0*4:
      StencilElementWeight=(1.0-xLoc[0])*xLoc[1]*(1.0-xLoc[2]);
      break;
    case 1+1*2+0*4:
      StencilElementWeight=xLoc[0]*xLoc[1]*(1.0-xLoc[2]);
      break;

    case 0+0*2+1*4:
      StencilElementWeight=(1.0-xLoc[0])*(1.0-xLoc[1])*xLoc[2];
      break;
    case 1+0*2+1*4:
      StencilElementWeight=xLoc[0]*(1.0-xLoc[1])*xLoc[2];
      break;
    case 0+1*2+1*4:
      StencilElementWeight=(1.0-xLoc[0])*xLoc[1]*xLoc[2];
      break;
    case 1+1*2+1*4:
      StencilElementWeight=xLoc[0]*xLoc[1]*xLoc[2];
      break;

    default:
      exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    StencilTable[di][dj][dk].Stencil.MultiplyScalar(StencilElementWeight);
    Stencil.Add(&StencilTable[di][dj][dk].Stencil);
  }
}

//init stencil for the corner based interpolation
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CornerBased::InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CornerBased::cStencil &Stencil,double *InterpolationCoefficientTable) {
  int iStencil,jStencil,kStencil,iX[3],nd,idim;
  double w,xLoc[3],dx[3],*xMinNode,*xMaxNode;
  PIC::Mesh::cDataCornerNode* CornerNode;
  PIC::Mesh::cDataBlockAMR *block;

  if (node==NULL) node=PIC::Mesh::mesh->findTreeNode(x);

  xMinNode=node->xmin;
  xMaxNode=node->xmax;

  //the size of the cells
  dx[0]=(xMaxNode[0]-xMinNode[0])/_BLOCK_CELLS_X_;
  dx[1]=(xMaxNode[1]-xMinNode[1])/_BLOCK_CELLS_Y_;
  dx[2]=(xMaxNode[2]-xMinNode[2])/_BLOCK_CELLS_Z_;

  //get the local coordinate for the interpolation point location
  #pragma ivdep
  for (idim=0;idim<3;idim++) {
    if ((x[idim]<xMinNode[idim])||(x[idim]>xMaxNode[idim])) exit(__LINE__,__FILE__,"Error: the point is out of block");
    if (fabs(x[idim]-xMaxNode[idim])<1e-10*dx[idim]) x[idim]= xMaxNode[idim]-1e-10*dx[idim];    
    xLoc[idim]=(x[idim]-xMinNode[idim])/dx[idim];
    iX[idim]=(int)(xLoc[idim]);
    xLoc[idim]-=iX[idim];
  }
  
  //build interpolation stencil
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  Stencil.flush();

  //return an empty stencil if the block is not allocated
  if ((block=node->block)==NULL) return; // PIC::InterpolationRoutines::CornerBased::StencilTable+ThreadOpenMP;

  for (iStencil=0;iStencil<2;iStencil++) for (jStencil=0;jStencil<2;jStencil++) for (kStencil=0;kStencil<2;kStencil++) {
    //get the local ID of the 'corner node'
    nd=_getCornerNodeLocalNumber(iStencil+iX[0],jStencil+iX[1],kStencil+iX[2]);
    CornerNode=block->GetCornerNode(nd);

    if (CornerNode!=NULL) {
      switch (iStencil+2*jStencil+4*kStencil) {
      case 0+0*2+0*4:
        w=(1.0-xLoc[0])*(1.0-xLoc[1])*(1.0-xLoc[2]);
        InterpolationCoefficientTable[0]=w;
        break;
      case 1+0*2+0*4:
        w=xLoc[0]*(1.0-xLoc[1])*(1.0-xLoc[2]);
        InterpolationCoefficientTable[1]=w;
        break;
      case 0+1*2+0*4:
        w=(1.0-xLoc[0])*xLoc[1]*(1.0-xLoc[2]);
        InterpolationCoefficientTable[3]=w;
        break;
      case 1+1*2+0*4:
        w=xLoc[0]*xLoc[1]*(1.0-xLoc[2]);
        InterpolationCoefficientTable[2]=w;
        break;

      case 0+0*2+1*4:
        w=(1.0-xLoc[0])*(1.0-xLoc[1])*xLoc[2];
        InterpolationCoefficientTable[4]=w;
        break;
      case 1+0*2+1*4:
        w=xLoc[0]*(1.0-xLoc[1])*xLoc[2];
        InterpolationCoefficientTable[5]=w;
        break;
      case 0+1*2+1*4:
        w=(1.0-xLoc[0])*xLoc[1]*xLoc[2];
        InterpolationCoefficientTable[7]=w;
        break;
      case 1+1*2+1*4:
        w=xLoc[0]*xLoc[1]*xLoc[2];
        InterpolationCoefficientTable[6]=w;
        break;

      default:
        exit(__LINE__,__FILE__,"Error: the option is not defined");
      }

      Stencil.AddCell(w,CornerNode,nd);
    }
    else {
      switch (iStencil+2*jStencil+4*kStencil) {
      case 0+0*2+0*4:
        InterpolationCoefficientTable[0]=0.0;
        break;
      case 1+0*2+0*4:
        InterpolationCoefficientTable[1]=0.0;
        break;
      case 0+1*2+0*4:
        InterpolationCoefficientTable[3]=0.0;
        break;
      case 1+1*2+0*4:
        InterpolationCoefficientTable[2]=0.0;
        break;

      case 0+0*2+1*4:
        InterpolationCoefficientTable[4]=0.0;
        break;
      case 1+0*2+1*4:
        InterpolationCoefficientTable[5]=0.0;
        break;
      case 0+1*2+1*4:
        InterpolationCoefficientTable[7]=0.0;
        break;
      case 1+1*2+1*4:
        InterpolationCoefficientTable[6]=0.0;
        break;

      default:
        exit(__LINE__,__FILE__,"Error: the option is not defined");
      }
    }
  }

  Stencil.Normalize();
}

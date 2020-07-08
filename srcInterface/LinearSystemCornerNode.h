//$Id$
//interface to SWMF's GMRES solver

/*
 * LinearSystemCornerNode.h
 *
 *  Created on: Nov 20, 2017
 *      Author: vtenishe
 */

#include "linear_solver_wrapper_c.h"
#include <functional>
#include <iostream>

#if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__ON_
#include <immintrin.h>
#endif

#ifndef __PGI
#include <immintrin.h>
#endif

#ifndef _LINEARSYSTEMCORNERNODE_H_
#define _LINEARSYSTEMCORNERNODE_H_

class cLinearSystemCornerNodeDataRequestListElement {
public:
  int CornerNodeID;
  cAMRnodeID NodeID;
};

class cRebuildMatrix {
public:
  virtual void RebuildMatrix() = 0;
};

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
class cLinearSystemCornerNode : public cRebuildMatrix {
public:
  double *SubdomainPartialRHS,*SubdomainPartialUnknownsVector;
  int SubdomainPartialUnknownsVectorLength;

  class cStencilElementData {
  public: 
    int UnknownVectorIndex,iVar,Thread;
    double MatrixElementValue;
  
    cStencilElementData() {
      UnknownVectorIndex=-1,iVar=-1,Thread=-1;
      MatrixElementValue=0.0;
    }
  };

  class cStencilElement {
  public:
    cCornerNode* CornerNode;
    int CornerNodeID;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

    //the user-model parameters neede for recalculating of the matrix coefficients
    double MatrixElementParameterTable[MaxMatrixElementParameterTableLength];
    int MatrixElementParameterTableLength;

    //the user-defined support used for recalculating the matrix coefficients
    void *MatrixElementSupportTable[MaxMatrixElementSupportTableLength];
    int MatrixElementSupportTableLength;

    cStencilElement()  {
      node=NULL,CornerNode=NULL;
      CornerNodeID=-1;

      MatrixElementParameterTableLength=0,MatrixElementSupportTableLength=0;

      for (int i=0;i<MaxMatrixElementParameterTableLength;i++) MatrixElementParameterTable[i]=0.0;
      for (int i=0;i<MaxMatrixElementSupportTableLength;i++) MatrixElementSupportTable[i]=NULL;
    }
  };

  class cRhsSupportTable {
  public:
    double Coefficient;
    char *AssociatedDataPointer;
  };

  class cMatrixRow {
  public:
    cMatrixRow* next;
    double Rhs;

    cStencilElement Elements[MaxStencilLength];
    cStencilElementData ElementDataTable[MaxStencilLength];
    int nNonZeroElements;

    //pointed to the beginitg of the associated data vectros used in calculating of the right-hand side vectors
    cRhsSupportTable RhsSupportTable_CornerNodes[MaxRhsSupportLength_CornerNodes];
    int RhsSupportLength_CornerNodes;

    cRhsSupportTable RhsSupportTable_CenterNodes[MaxRhsSupportLength_CenterNodes];
    int RhsSupportLength_CenterNodes;

    int i,j,k;
    int iVar;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    cCornerNode* CornerNode;

    cMatrixRow() {
      next=NULL,Rhs=0.0,node=NULL,CornerNode=NULL;
      i=0,j=0,k=0,iVar=0,nNonZeroElements=0,RhsSupportLength_CornerNodes=0,RhsSupportLength_CenterNodes=0;
    }
  };

  cStack<cMatrixRow> MatrixRowStack;

  //Table of the rows local to the current MPI process
  cMatrixRow *MatrixRowListFirst,*MatrixRowListLast;
  cMatrixRow **MatrixRowTable;
  int MatrixRowTableLength;

  //exchange buffers
  int* RecvExchangeBufferLength;  //the number of elementf in the recv list
  int* SendExchangeBufferLength;
  int **SendExchangeBufferElementIndex;

  //add new row to the matrix
  //for that a user function is called that returns stenal (elements of the matrix), corresponding rhs, and the number of the elements in the amtrix' line
  class cMatrixRowNonZeroElementTable {
  public:
    int i,j,k; //coordintes of the corner block used in the equation
    int iVar; //the index of the variable used in the stencil
    double MatrixElementValue;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node;

    //the following are needed only in case when the periodic bounday conditinos are imposed
    bool BoundaryNodeFlag;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *OriginalNode;

    //the user-model parameters neede for recalculating of the matrix coefficients
    double MatrixElementParameterTable[MaxMatrixElementParameterTableLength];
    int MatrixElementParameterTableLength;

    //the user-defined support used for recalculating the matrix coefficients
    void *MatrixElementSupportTable[MaxMatrixElementSupportTableLength];
    int MatrixElementSupportTableLength;

    cMatrixRowNonZeroElementTable() {
      i=0,j=0,k=0,iVar=0,MatrixElementValue=0.0,Node=NULL;
      BoundaryNodeFlag=false,OriginalNode=NULL;

      MatrixElementParameterTableLength=0,MatrixElementSupportTableLength=0;

      for (int i=0;i<MaxMatrixElementParameterTableLength;i++) MatrixElementParameterTable[i]=0.0;
      for (int i=0;i<MaxMatrixElementSupportTableLength;i++) MatrixElementSupportTable[i]=NULL;
    }
  };

  cMatrixRowNonZeroElementTable MatrixRowNonZeroElementTable[MaxStencilLength];

  //provcess the right boundary
  int nVirtualBlocks,nRealBlocks;
  void ProcessRightDomainBoundary(int *RecvDataPointCounter,void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node));


  //void reset the data of the obsect to the default state
  void DeleteDataBuffers();
  void Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
  void Reset();

  //reset indexing of the nodes
  void ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

  //build the matrix
  void(*GetStencil)(int,int,int,int,cMatrixRowNonZeroElementTable*,int&,double&,cRhsSupportTable*,int&,cRhsSupportTable*,int&,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*); 
  void BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node));

  //re-build the matrix after the domain re-decomposition operation
  void RebuildMatrix() {
    if (GetStencil!=NULL) BuildMatrix(GetStencil);
  }

  //constructor
  cLinearSystemCornerNode() {
    MatrixRowListFirst=NULL,MatrixRowListLast=NULL;
    MatrixRowTable=NULL,MatrixRowTableLength=0;

    //exchange buffers
    RecvExchangeBufferLength=NULL;
    SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;

    //partial data of the linear system to be solved
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
    SubdomainPartialUnknownsVectorLength=0;

    nVirtualBlocks=0,nRealBlocks=0;
    GetStencil=NULL;
  }

  //exchange the data
  void ExchangeIntermediateUnknownsData(double *x,double** &);

  //matrix/vector multiplication
  void MultiplyVector(double *p,double *x,int length);

  //call the linear system solver, and unpack the solution afterward
  void Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),void (*fUnpackSolution)(double* x,cCornerNode* CornerNode),double Tol,int nMaxIter,
      int (*fPackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer),
      int (*fUnpackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer));

  //update the RHS vector
  void UpdateRhs(double (*fSetRhs)(int,cRhsSupportTable*,int,cRhsSupportTable*,int));

  //update non-zero coefficients of the matrix
  void UpdateMatrixNonZeroCoefficients(void (*UpdateMatrixRow)(cMatrixRow*));

  //destructor
  ~cLinearSystemCornerNode() {
    if (MatrixRowTable!=NULL) delete [] MatrixRowTable;
  }

  //calculate signature of the matrix
  void GetSignature(long int nline,const char* fname) {
    CRC32 Signature,SignatureRhs;
    cMatrixRow* row;
    int cnt,iElementMax,iElement;
    cStencilElementData *data,*ElementDataTable;

    for (row=MatrixRowListFirst,cnt=0;row!=NULL;row=row->next,cnt++) {
      Signature.add(cnt);

      SignatureRhs.add(cnt);
      SignatureRhs.add(row->Rhs);

      iElementMax=row->nNonZeroElements;
      ElementDataTable=row->ElementDataTable;

      for (iElement=0;iElement<iElementMax;iElement++) {
        data=ElementDataTable+iElement;

        Signature.add(data->MatrixElementValue);
        Signature.add(data->Thread);
        Signature.add(data->iVar);
        Signature.add(data->UnknownVectorIndex);
      }
    }

    Signature.PrintChecksum(nline,fname);
    SignatureRhs.PrintChecksum(nline,fname);
  }

};

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    PIC::Mesh::cDataBlockAMR *block=NULL;
    PIC::Mesh::cDataCornerNode *CornerNode=NULL;

    if ((block=node->block)!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
        if ((CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k)))!=NULL) {
          CornerNode->LinearSolverUnknownVectorIndex=-1;
        }
      }
    }

  }
  else for (int i=0;i<(1<<DIM);i++) if (node->downNode[i]!=NULL) ResetUnknownVectorIndex(node->downNode[i]);
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::ProcessRightDomainBoundary(int* RecvDataPointCounter,void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)) {
  int i,j,k,iface,iblock;
  cMatrixRow* NewRow;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  cStencilElement* el;
  cStencilElementData* el_data;
//  int Index=PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;
  int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

//  int Debug=0;

  struct cOffset {
    int di,dj,dk;
  };


  nVirtualBlocks=0;

//  return ;

  cOffset xOffset[1]={{_BLOCK_CELLS_X_,0,0}};
  cOffset yOffset[1]={{0,_BLOCK_CELLS_Y_,0}};
  cOffset zOffset[1]={{0,0,_BLOCK_CELLS_Z_}};

  cOffset xyOffset[3]={{_BLOCK_CELLS_X_,0,0},{0,_BLOCK_CELLS_Y_,0},{_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0}};
  cOffset xzOffset[3]={{_BLOCK_CELLS_X_,0,0},{0,0,_BLOCK_CELLS_Z_},{_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_}};
  cOffset yzOffset[3]={{0,_BLOCK_CELLS_Y_,0},{0,0,_BLOCK_CELLS_Z_},{0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}};

  cOffset xyzOffset[7]={{_BLOCK_CELLS_X_,0,0},{0,_BLOCK_CELLS_Y_,0},{0,0,_BLOCK_CELLS_Z_},
      {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0},{_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_},{0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_},
      {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}};

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    bool flag=false;

    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

    if (node->Thread==PIC::ThisThread) for (iface=1;iface<6;iface+=2) if (node->GetNeibFace(iface,0,0,&PIC::Mesh::mesh)==NULL) {
      flag=true;
      break;
    }

    if (flag==false) continue;

/*Debug++;

//debug -> count the number of rows and blocks
int ntotbl=0;
for ( cMatrixRow* Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) ntotbl++;*/

    //the block has at least one face at the 'right' boundary of the domain
    //determine the combination of the faces that are at the boundary
    int xBoundaryFace=0,yBoundaryFace=0,zBoundaryFace=0;

    if (node->GetNeibFace(1,0,0,&PIC::Mesh::mesh)==NULL) xBoundaryFace=10;
    if (node->GetNeibFace(3,0,0,&PIC::Mesh::mesh)==NULL) yBoundaryFace=100;
    if (node->GetNeibFace(5,0,0,&PIC::Mesh::mesh)==NULL) zBoundaryFace=1000;

    //build "virtual" blocks for the combination of the boundaries
    cOffset *OffsetTable;
    int OffsetTableLength;

    switch (xBoundaryFace+yBoundaryFace+zBoundaryFace) {
     case 10+0+0:
      OffsetTable=xOffset;
      OffsetTableLength=1;
      break;

     case 0+100+0:
      OffsetTable=yOffset;
      OffsetTableLength=1;
      break;

     case 0+0+1000:
      OffsetTable=zOffset;
      OffsetTableLength=1;
      break;

     case 10+100+0:
      OffsetTable=xyOffset;
      OffsetTableLength=3;
      break;

     case 10+0+1000:
      OffsetTable=xzOffset;
      OffsetTableLength=3;
      break;

     case 0+100+1000:
      OffsetTable=yzOffset;
      OffsetTableLength=3;
      break;

     case 10+100+1000:
      OffsetTable=xyzOffset;
      OffsetTableLength=7;
      break;
    }

    nVirtualBlocks+=OffsetTableLength;

    //loop through all 'virtual' blocks
    for (iblock=0;iblock<OffsetTableLength;iblock++) {
      //loop through all 'internal' corner nodes of the block
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        if ( (i+OffsetTable[iblock].di<=_BLOCK_CELLS_X_) && (j+OffsetTable[iblock].dj<=_BLOCK_CELLS_Y_) && (k+OffsetTable[iblock].dk<=_BLOCK_CELLS_Z_) ) {
          //the point is located at the boundary of the domain and the user-defined funtion is needed to be called to get the interpolation stencil
          int NonZeroElementsFound;
          double rhs;

          NewRow=MatrixRowStack.newElement();

          //call the user defined function to determine the non-zero elements of the matrix
          NewRow->RhsSupportLength_CornerNodes=0;
          NewRow->RhsSupportLength_CenterNodes=0;
          rhs=0.0;

          for (int ii=0;ii<MaxStencilLength;ii++) MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength=0,MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength=0;

          f(i+OffsetTable[iblock].di,j+OffsetTable[iblock].dj,k+OffsetTable[iblock].dk,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,NewRow->RhsSupportTable_CornerNodes,NewRow->RhsSupportLength_CornerNodes,NewRow->RhsSupportTable_CenterNodes,NewRow->RhsSupportLength_CenterNodes,node);
          if (NonZeroElementsFound>MaxStencilLength) exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");

          //populate the new row
          NewRow->i=i+OffsetTable[iblock].di,NewRow->j=j+OffsetTable[iblock].dj,NewRow->k=k+OffsetTable[iblock].dk;
          NewRow->iVar=iVar;
          NewRow->node=node;
          NewRow->Rhs=rhs;
          NewRow->CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(NewRow->i,NewRow->j,NewRow->k));
          NewRow->nNonZeroElements=NonZeroElementsFound;

          //add the new row to the matrix
          if (MatrixRowListFirst==NULL) {
            MatrixRowListLast=NewRow;
            MatrixRowListFirst=NewRow;
            NewRow->next=NULL;
          }
          else {
            MatrixRowListLast->next=NewRow;
            MatrixRowListLast=NewRow;
            NewRow->next=NULL;
          }

          //add to the row non-zero elements
          for (int ii=0;ii<NonZeroElementsFound;ii++) {
            MatrixRowNonZeroElementTable[ii].Node=node;

            if ((MatrixRowNonZeroElementTable[ii].i>=iMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,&PIC::Mesh::mesh)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,&PIC::Mesh::mesh);
            }
            else if (MatrixRowNonZeroElementTable[ii].i<0) {
              MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0,&PIC::Mesh::mesh);
            }

            if ((MatrixRowNonZeroElementTable[ii].j>=jMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,&PIC::Mesh::mesh)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,&PIC::Mesh::mesh);
            }
            else if (MatrixRowNonZeroElementTable[ii].j<0) {
              MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0,&PIC::Mesh::mesh);
            }

            if ((MatrixRowNonZeroElementTable[ii].k>=kMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,&PIC::Mesh::mesh)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].k-=_BLOCK_CELLS_Z_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,&PIC::Mesh::mesh);
            }
            else if (MatrixRowNonZeroElementTable[ii].k<0) {
              MatrixRowNonZeroElementTable[ii].k+=_BLOCK_CELLS_Z_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(4,0,0,&PIC::Mesh::mesh);
            }

            if (MatrixRowNonZeroElementTable[ii].Node==NULL) {
              exit(__LINE__,__FILE__,"Error: the block is not found");
            }

            cCornerNode* CornerNode=MatrixRowNonZeroElementTable[ii].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k));

            if (CornerNode==NULL) exit(__LINE__,__FILE__,"Error: something is wrong");

            el=NewRow->Elements+ii;
            el_data=NewRow->ElementDataTable+ii;

            el->CornerNode=CornerNode;
            el->CornerNodeID=PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k);
            el_data->UnknownVectorIndex=-1;
            el_data->MatrixElementValue=MatrixRowNonZeroElementTable[ii].MatrixElementValue;
            el->node=MatrixRowNonZeroElementTable[ii].Node;
            el_data->iVar=MatrixRowNonZeroElementTable[ii].iVar;

            el->MatrixElementParameterTableLength=MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength;
            memcpy(el->MatrixElementParameterTable,MatrixRowNonZeroElementTable[ii].MatrixElementParameterTable,el->MatrixElementParameterTableLength*sizeof(double));

            el->MatrixElementSupportTableLength=MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength;
            memcpy(el->MatrixElementSupportTable,MatrixRowNonZeroElementTable[ii].MatrixElementSupportTable,el->MatrixElementSupportTableLength*sizeof(void*));

            //count the number of the element that are needed
            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[MatrixRowNonZeroElementTable[ii].Node->Thread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
        }
        else {
          //the point is outside of the domain
          NewRow=MatrixRowStack.newElement();

          NewRow->i=-1000000,NewRow->j=-1000000,NewRow->k=-1000000;
          NewRow->iVar=iVar;
          NewRow->node=node;
          NewRow->Rhs=1.0;
          NewRow->CornerNode=NULL;
          NewRow->nNonZeroElements=1;
          NewRow->RhsSupportLength_CornerNodes=0;
          NewRow->RhsSupportLength_CenterNodes=0;

          NewRow->Elements[0].CornerNode=NULL;
          NewRow->Elements[0].CornerNodeID=-1;
          NewRow->ElementDataTable[0].UnknownVectorIndex=-1; //Index;
          NewRow->ElementDataTable[0].MatrixElementValue=1.0;
          NewRow->Elements[0].node=node;
          NewRow->ElementDataTable[0].iVar=iVar;

          //add the new row to the matrix
          if (MatrixRowListFirst==NULL) {
            MatrixRowListLast=NewRow;
            MatrixRowListFirst=NewRow;
            NewRow->next=NULL;
          }
          else {
            MatrixRowListLast->next=NewRow;
            MatrixRowListLast=NewRow;
            NewRow->next=NULL;
          }

        }
      }
    }
  }
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable_CornerNodes,int &RhsSupportLength_CornerNodes,cRhsSupportTable* RhsSupportTable_CenterNodes,int &RhsSupportLength_CenterNodes,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k,thread,nLocalNode;
  int iRow=0;
  cMatrixRow* Row;
  cStencilElement* el;
  cStencilElementData* el_data;
  cCornerNode* CornerNode;

//  int Debug=0;

  nRealBlocks=0;

  //reset indexing of the nodes
  Reset();

  //save the user-defined build stencil function
  GetStencil=f;

  //allocate the counter of the data points to be recieve
  int RecvDataPointCounter[PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) RecvDataPointCounter[thread]=0;

  list<cLinearSystemCornerNodeDataRequestListElement> *DataRequestList=new list<cLinearSystemCornerNodeDataRequestListElement> [PIC::nTotalThreads];

  //build the matrix
  for (nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block==NULL) continue;
    //in case of the periodic boundary condition it is only the points that are inside the "real" computational domain that are considered
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;

      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0,&PIC::Mesh::mesh)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }

      if (BoundaryBlock==true) continue;
    }

    nRealBlocks++;

    //the limits are correct: the point i==_BLOCK_CELLS_X_ belongs to the onother block
    int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

    for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
      for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        int NonZeroElementsFound;
        double rhs;

        //create the new Row entry
        cMatrixRow* NewRow=MatrixRowStack.newElement();

        //call the user defined function to determine the non-zero elements of the matrix
        NewRow->RhsSupportLength_CornerNodes=0;
        NewRow->RhsSupportLength_CenterNodes=0;
        rhs=0.0;

        for (int ii=0;ii<MaxStencilLength;ii++) MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength=0,MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength=0;

        f(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,NewRow->RhsSupportTable_CornerNodes,NewRow->RhsSupportLength_CornerNodes,NewRow->RhsSupportTable_CenterNodes,NewRow->RhsSupportLength_CenterNodes,node);

        if (NonZeroElementsFound==0) {
          //the point is not included into the matrix
          MatrixRowStack.deleteElement(NewRow);
          continue;
        }

        if (NonZeroElementsFound>MaxStencilLength) {
          exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");
        }

        if (NewRow->RhsSupportLength_CenterNodes>MaxRhsSupportLength_CenterNodes) {
          exit(__LINE__,__FILE__,"Error: NewRow->RhsSupportLength_CenterNodes>MaxRhsSupportLength_CenterNodes; Need to increase the value of MaxRhsSupportLength_CenterNodes");
        }

        if (NewRow->RhsSupportLength_CornerNodes>MaxRhsSupportLength_CornerNodes) {
          exit(__LINE__,__FILE__,"Error: NewRow->RhsSupportLength_CornerNodes>MaxRhsSupportLength_CornerNodes; Need to increase the value of MaxRhsSupportLength_CornerNodes");
        }

        //scan through the found stencil and correct blocks and indexing is needed
        for (int ii=0;ii<NonZeroElementsFound;ii++) {
          MatrixRowNonZeroElementTable[ii].Node=node;

          if ((MatrixRowNonZeroElementTable[ii].i>=iMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,&PIC::Mesh::mesh)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0,&PIC::Mesh::mesh);
          }
          else if (MatrixRowNonZeroElementTable[ii].i<0) {
            MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0,&PIC::Mesh::mesh);
          }

          if ((MatrixRowNonZeroElementTable[ii].j>=jMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,&PIC::Mesh::mesh)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0,&PIC::Mesh::mesh);
          }
          else if (MatrixRowNonZeroElementTable[ii].j<0) {
            MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0,&PIC::Mesh::mesh);
          }

          if ((MatrixRowNonZeroElementTable[ii].k>=kMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,&PIC::Mesh::mesh)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].k-=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0,&PIC::Mesh::mesh);
          }
          else if (MatrixRowNonZeroElementTable[ii].k<0) {
            MatrixRowNonZeroElementTable[ii].k+=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(4,0,0,&PIC::Mesh::mesh);
          }

          if (MatrixRowNonZeroElementTable[ii].Node==NULL) {
            exit(__LINE__,__FILE__,"Error: the block is not found");
          }

          //check whether the periodic boundary conditions are in use, and the new block is on the boundary of the domain
          if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
            bool BoundaryBlock=false;

            for (int iface=0;iface<6;iface++) if (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(iface,0,0,&PIC::Mesh::mesh)==NULL) {
              //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
              BoundaryBlock=true;
              break;
            }

            MatrixRowNonZeroElementTable[ii].BoundaryNodeFlag=BoundaryBlock;
            MatrixRowNonZeroElementTable[ii].OriginalNode=MatrixRowNonZeroElementTable[ii].Node;

            if (BoundaryBlock==true) {
              //the block is at the domain boundary -> find the point located in the "real" part of the computational domain
              MatrixRowNonZeroElementTable[ii].Node=PIC::BC::ExternalBoundary::Periodic::findCorrespondingRealBlock(MatrixRowNonZeroElementTable[ii].Node);
            }
          }

        }

        //populate the new row
        NewRow->i=i,NewRow->j=j,NewRow->k=k;
        NewRow->iVar=iVar;
        NewRow->node=node;
        NewRow->Rhs=rhs;
        NewRow->CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
        NewRow->nNonZeroElements=NonZeroElementsFound;


        if (MatrixRowListFirst==NULL) {
          MatrixRowListLast=NewRow;
          MatrixRowListFirst=NewRow;
          NewRow->next=NULL;
        }
        else {
          MatrixRowListLast->next=NewRow;
          MatrixRowListLast=NewRow;
          NewRow->next=NULL;
        }

        //add to the row non-zero elements
        for (int iElement=0;iElement<NonZeroElementsFound;iElement++) {
          el=NewRow->Elements+iElement;
          el_data=NewRow->ElementDataTable+iElement;

          if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_) {
            CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
          }
          else  {
            if (MatrixRowNonZeroElementTable[iElement].BoundaryNodeFlag==false) {
              CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
            }
            else {
              if (MatrixRowNonZeroElementTable[iElement].Node->Thread!=PIC::ThisThread) {
                CornerNode=MatrixRowNonZeroElementTable[iElement].OriginalNode->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
              }
              else {
                CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));
              }
            }
          }

          el->CornerNode=CornerNode;
          el->CornerNodeID=PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k);
          el_data->UnknownVectorIndex=-1;
          el_data->MatrixElementValue=MatrixRowNonZeroElementTable[iElement].MatrixElementValue;
          el->node=MatrixRowNonZeroElementTable[iElement].Node;
          el_data->iVar=MatrixRowNonZeroElementTable[iElement].iVar;

          el->MatrixElementParameterTableLength=MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength;
          memcpy(el->MatrixElementParameterTable,MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable,el->MatrixElementParameterTableLength*sizeof(double));

          el->MatrixElementSupportTableLength=MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength;
          memcpy(el->MatrixElementSupportTable,MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTable,el->MatrixElementSupportTableLength*sizeof(void*));

          //count the number of the element that are needed
          if (MatrixRowNonZeroElementTable[iElement].Node->Thread==PIC::ThisThread) {
            //the point is located at the current MPI process

            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[PIC::ThisThread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
          else {
            //the data point is on another MPI process

            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[MatrixRowNonZeroElementTable[iElement].Node->Thread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
        }


      }
    }
  }


  //add 'virtual' blocks at the right boundary of the domain
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_) {
    //ProcessRightDomainBoundary(RecvDataPointCounter,f);
  }


  //allocate the data buffers for the partial vectors and link the pointers points that are inside the subdomian
  //allocate the exchange buffer for the data the needs to be recieved from other MPI processess
  int iElement;

  int *DataExchangeTableCounter=new int [PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) DataExchangeTableCounter[thread]=0;

  for (iRow=0,Row=MatrixRowListFirst;Row!=NULL;iRow++,Row=Row->next) {

    for (iElement=0;iElement<Row->nNonZeroElements;iElement++) {
      el=Row->Elements+iElement;
      el_data=Row->ElementDataTable+iElement;

      if (el->CornerNode!=NULL) {
        if (el->CornerNode->LinearSolverUnknownVectorIndex==-2) {
          //the node is still not linked to the 'partial data vector'
          el->CornerNode->LinearSolverUnknownVectorIndex=DataExchangeTableCounter[el->node->Thread]++;

          if (el->node->Thread!=PIC::ThisThread) {
            //add the point to the data request list
            cLinearSystemCornerNodeDataRequestListElement DataRequestListElement;
            cAMRnodeID NodeID;

            PIC::Mesh::mesh.GetAMRnodeID(NodeID,el->node);
            DataRequestListElement.CornerNodeID=el->CornerNodeID;
            DataRequestListElement.NodeID=NodeID;

            DataRequestList[el->node->Thread].push_back(DataRequestListElement);
          }
        }

        el_data->UnknownVectorIndex=el->CornerNode->LinearSolverUnknownVectorIndex;
        el_data->Thread=el->node->Thread;
      }
      else {
        //the point is within the 'virtual' block
        el_data->UnknownVectorIndex=DataExchangeTableCounter[el->node->Thread]++;
        el_data->Thread=el->node->Thread;
      }
    }
  }

  //re-count the 'internal' corner nodes so they follow the i,j,k order as well as rows
  for (Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) if (Row->CornerNode!=NULL) Row->CornerNode->LinearSolverUnknownVectorIndex=-1;

  for (iRow=0,Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) {
    if (Row->CornerNode!=NULL) {
      if (Row->iVar==0) {
        if (Row->CornerNode->LinearSolverUnknownVectorIndex>=0) exit(__LINE__,__FILE__,"Error: a courner node is double counted");
        Row->CornerNode->LinearSolverUnknownVectorIndex=iRow;
      }
    }
    else Row->ElementDataTable[0].UnknownVectorIndex=iRow; //rows corresponding to the virtual blocks has only one non-zero element

    if (Row->iVar==NodeUnknownVariableVectorLength-1) iRow++; //all rows pointing to the same 'CornerNode' have the same 'iRow'
  }

  //verify that all all corner nodes have been counted
  for (Row=MatrixRowListFirst;Row!=NULL;iRow++,Row=Row->next) for (iElement=0;iElement<Row->nNonZeroElements;iElement++) {
    int n;

    if (Row->Elements[iElement].CornerNode!=NULL) {
      if ((n=Row->Elements[iElement].CornerNode->LinearSolverUnknownVectorIndex)<0) {
        exit(__LINE__,__FILE__,"Error: an uncounted corner node has been found");
      }

      Row->ElementDataTable[iElement].UnknownVectorIndex=n;
    }
  }


  //exchange the request lists
  int From,To;
  int *nGlobalDataPointTable=new int [PIC::nTotalThreads*PIC::nTotalThreads];
  cLinearSystemCornerNodeDataRequestListElement* ExchangeList;

  MPI_Allgather(DataExchangeTableCounter,PIC::nTotalThreads,MPI_INT,nGlobalDataPointTable,PIC::nTotalThreads,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

  SendExchangeBufferLength=new int [PIC::nTotalThreads];
  SendExchangeBufferElementIndex=new int* [PIC::nTotalThreads];

  RecvExchangeBufferLength=new int [PIC::nTotalThreads];

  for (thread=0;thread<PIC::nTotalThreads;thread++) {
    SendExchangeBufferLength[thread]=0,SendExchangeBufferElementIndex[thread]=0;
    RecvExchangeBufferLength[thread]=0;
  }


  //create the lists
  for (To=0;To<PIC::nTotalThreads;To++) for (From=0;From<PIC::nTotalThreads;From++) if ( (To!=From) && (nGlobalDataPointTable[From+To*PIC::nTotalThreads]!=0) && ((To==PIC::ThisThread)||(From==PIC::ThisThread)) ) {
    ExchangeList=new cLinearSystemCornerNodeDataRequestListElement [nGlobalDataPointTable[From+To*PIC::nTotalThreads]];

    if (PIC::ThisThread==To) {
      list<cLinearSystemCornerNodeDataRequestListElement>::iterator ptr;

      for (i=0,ptr=DataRequestList[From].begin();ptr!=DataRequestList[From].end();ptr++,i++) ExchangeList[i]=*ptr;
      MPI_Send(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cLinearSystemCornerNodeDataRequestListElement),MPI_CHAR,From,0,MPI_GLOBAL_COMMUNICATOR);

      //create the recv list
      RecvExchangeBufferLength[From]=nGlobalDataPointTable[From+To*PIC::nTotalThreads];
    }
    else {
      MPI_Status status;

      MPI_Recv(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cLinearSystemCornerNodeDataRequestListElement),MPI_CHAR,To,0,MPI_GLOBAL_COMMUNICATOR,&status);

      //unpack the SendDataList
      SendExchangeBufferLength[To]=nGlobalDataPointTable[From+To*PIC::nTotalThreads];
      SendExchangeBufferElementIndex[To]=new int[SendExchangeBufferLength[To]];

      for (int ii=0;ii<SendExchangeBufferLength[To];ii++) {
        node=PIC::Mesh::mesh.findAMRnodeWithID(ExchangeList[ii].NodeID);
        CornerNode=node->block->GetCornerNode(ExchangeList[ii].CornerNodeID);

        SendExchangeBufferElementIndex[To][ii]=CornerNode->LinearSolverUnknownVectorIndex;
      }
    }

    delete [] ExchangeList;
  }

  //deallocate 'nGlobalDataPointTable'
  delete [] nGlobalDataPointTable;
  delete [] DataExchangeTableCounter;

  delete [] DataRequestList;

  //cleate a RowTable;
  for (MatrixRowTableLength=0,Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) MatrixRowTableLength++;

  if (MatrixRowTable!=NULL) delete [] MatrixRowTable;
  MatrixRowTable=new cMatrixRow *[MatrixRowTableLength];

  for (MatrixRowTableLength=0,Row=MatrixRowListFirst;Row!=NULL;Row=Row->next) MatrixRowTable[MatrixRowTableLength++]=Row;
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::ExchangeIntermediateUnknownsData(double *x,double** &LocalRecvExchangeBufferTable) {
  int To,From;

  LocalRecvExchangeBufferTable=new double *[PIC::nTotalThreads];
  
  MPI_Datatype *unknown_vector_type_table=new MPI_Datatype[PIC::nTotalThreads];
  
  for (int thread=0;thread<PIC::nTotalThreads;thread++) {
    if (NodeUnknownVariableVectorLength*RecvExchangeBufferLength[thread]>0) {
      LocalRecvExchangeBufferTable[thread]=new double [NodeUnknownVariableVectorLength*RecvExchangeBufferLength[thread]];
    }
    else {
      LocalRecvExchangeBufferTable[thread]=NULL;
    }  
  } 


  MPI_Request RecvRequestTable[PIC::nTotalThreads],SendRequestTable[PIC::nTotalThreads];
  int RecvRequestTableLength=0,SendRequestTableLength=0;

  //Initiate Recv operations
  for (From=0;From<PIC::nTotalThreads;From++) if ((From!=PIC::ThisThread)&&(RecvExchangeBufferLength[From]!=0)) {
    MPI_Irecv(LocalRecvExchangeBufferTable[From],NodeUnknownVariableVectorLength*RecvExchangeBufferLength[From],MPI_DOUBLE,From,0,MPI_GLOBAL_COMMUNICATOR,RecvRequestTable+RecvRequestTableLength);
    RecvRequestTableLength++;
  }

  //Init Send operations
  for (To=0;To<PIC::nTotalThreads;To++) if ((To!=PIC::ThisThread)&&(SendExchangeBufferLength[To]!=0)) {
    //this if the sending thread
    int i,offset=0;
    int Length=SendExchangeBufferLength[To];
    int *IndexTable=SendExchangeBufferElementIndex[To];

    int length_table[Length];
    int offset_table[Length];

    for (i=0;i<Length;i++) {
      length_table[i]=NodeUnknownVariableVectorLength;
      offset_table[i]=NodeUnknownVariableVectorLength*IndexTable[i]; 
    }

    MPI_Type_indexed(Length, length_table, offset_table, MPI_DOUBLE, unknown_vector_type_table+To);
    MPI_Type_commit(unknown_vector_type_table+To);

    for (i=0;i<Length;i++) {
      offset+=NodeUnknownVariableVectorLength;
    }

    if (offset!=NodeUnknownVariableVectorLength*SendExchangeBufferLength[To]) exit(__LINE__,__FILE__,"Error: out of limit");

    MPI_Isend(x,1,unknown_vector_type_table[To],To,0,MPI_GLOBAL_COMMUNICATOR,SendRequestTable+SendRequestTableLength);

    SendRequestTableLength++;
  }

  //wait for completing the send/rect operations
  MPI_Waitall(RecvRequestTableLength,RecvRequestTable,MPI_STATUSES_IGNORE);
  MPI_Waitall(SendRequestTableLength,SendRequestTable,MPI_STATUSES_IGNORE);

  //delete temporary buffers
  delete [] unknown_vector_type_table;
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::DeleteDataBuffers() {
  if (RecvExchangeBufferLength!=NULL) {
    //clear the row stack
    MatrixRowStack.resetStack();

    //deallocate the allocated data buffers
    for (int thread=0;thread<PIC::nTotalThreads;thread++) {

      if (SendExchangeBufferElementIndex!=NULL) if (SendExchangeBufferElementIndex[thread]!=NULL) {
        delete [] SendExchangeBufferElementIndex[thread];
      }
    }

    if (RecvExchangeBufferLength!=NULL) {
      delete [] RecvExchangeBufferLength;
    }

    RecvExchangeBufferLength=NULL;

    if (SendExchangeBufferElementIndex!=NULL) {
      delete [] SendExchangeBufferLength;
      delete [] SendExchangeBufferElementIndex;
    }

    SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;

    if (SubdomainPartialRHS!=NULL) {
      delete [] SubdomainPartialRHS;
      delete [] SubdomainPartialUnknownsVector;
    }

    MatrixRowListFirst=NULL,MatrixRowListLast=NULL;
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
  }
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::Reset() {
  Reset(PIC::Mesh::mesh.rootTree);
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  if ((startNode==PIC::Mesh::mesh.rootTree)&&(RecvExchangeBufferLength!=NULL)) DeleteDataBuffers();

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    PIC::Mesh::cDataBlockAMR *block=NULL;
    PIC::Mesh::cDataCornerNode *CornerNode=NULL;

    if ((block=startNode->block)!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
        if ((CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k)))!=NULL) {
          CornerNode->LinearSolverUnknownVectorIndex=-1;
        }
      }
    }

  }
  else for (int i=0;i<(1<<DIM);i++) if (startNode->downNode[i]!=NULL) Reset(startNode->downNode[i]);
}


//update non-zero coefficients of the matrix
template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::UpdateMatrixNonZeroCoefficients(void (*UpdateMatrixRow)(cMatrixRow*)) {

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int irow=0;irow<MatrixRowTableLength;irow++) {
    UpdateMatrixRow(MatrixRowTable[irow]);
  }
}


template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::UpdateRhs(double (*fSetRhs)(int,cRhsSupportTable*,int,cRhsSupportTable*,int)) {

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp parallel for schedule(dynamic,1) default(none) shared(fSetRhs)
  #endif
  for (int irow=0;irow<MatrixRowTableLength;irow++) {
    cMatrixRow* row=MatrixRowTable[irow];

    if ((row->RhsSupportLength_CornerNodes!=0)||(row->RhsSupportLength_CenterNodes!=0)) {
      row->Rhs=fSetRhs(row->iVar,row->RhsSupportTable_CornerNodes,row->RhsSupportLength_CornerNodes,row->RhsSupportTable_CenterNodes,row->RhsSupportLength_CenterNodes);
    }
  }
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::MultiplyVector(double *p,double *x,int length) {


  double **RecvExchangeBufferTable;

  ExchangeIntermediateUnknownsData(x,RecvExchangeBufferTable);
  RecvExchangeBufferTable[PIC::ThisThread]=x;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel default(none) shared (RecvExchangeBufferTable,PIC::ThisThread,PIC::nTotalThreads,p) firstprivate(length)
#endif
{
  double **LocalRecvExchangeBufferTable;
 
  LocalRecvExchangeBufferTable=new double*[PIC::nTotalThreads];
  for (int thread=0;thread<PIC::nTotalThreads;thread++) LocalRecvExchangeBufferTable[thread]=RecvExchangeBufferTable[thread];


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp for schedule(dynamic,4) 
#endif
  for (int irow=0;irow<MatrixRowTableLength;irow++) {

    int iElementMax=MatrixRowTable[irow]->nNonZeroElements;
    cStencilElementData *ElementDataTable=MatrixRowTable[irow]->ElementDataTable;
    
    double res=0.0;
    int iElement=0;
    cStencilElementData *data,*data_next,*data_next_next=NULL;
    double *u_vect,*u_vect_next;

    #if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__256_
    //add most of the vector
    for (;iElement+3<iElementMax;iElement+=4) {
      alignas(32) double a[4],b[4],*r;
      __m256d av,bv,cv,rv;

      data=ElementDataTable+iElement;
      /*
      a[0]=data->MatrixElementValue,b[0]=LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];

      data=ElementDataTable+iElement+1;
      a[1]=data->MatrixElementValue,b[1]=LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];

      data=ElementDataTable+iElement+2;
      a[2]=data->MatrixElementValue,b[2]=LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];

      data=ElementDataTable+iElement+3;
      a[3]=data->MatrixElementValue,b[3]=LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];
      */


      av=_mm256_set_pd((data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue);

      bv=_mm256_set_pd(
		      LocalRecvExchangeBufferTable[(data+3)->Thread][(data+3)->iVar+NodeUnknownVariableVectorLength*(data+3)->UnknownVectorIndex], 
		      LocalRecvExchangeBufferTable[(data+2)->Thread][(data+2)->iVar+NodeUnknownVariableVectorLength*(data+2)->UnknownVectorIndex], 
		      LocalRecvExchangeBufferTable[(data+1)->Thread][(data+1)->iVar+NodeUnknownVariableVectorLength*(data+1)->UnknownVectorIndex], 
		      LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex]
		      );


   //   av=_mm256_load_pd(a);
   //   bv=_mm256_load_pd(b);
      cv=_mm256_mul_pd(av,bv);
      rv=_mm256_hadd_pd(cv,cv);

      r=(double*)&rv;
      res+=r[1]+r[2];
    }

    #elif _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__512_
    //add most of the vector
    for (;iElement+7<iElementMax;iElement+=8) {
      __m512d av,bv,cv;

      data=ElementDataTable+iElement;

      av=_mm512_set_pd(
		      (data+7)->MatrixElementValue,(data+6)->MatrixElementValue,(data+5)->MatrixElementValue,(data+4)->MatrixElementValue,
		      (data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue);  

      bv=_mm512_set_pd(
		      LocalRecvExchangeBufferTable[(data+7)->Thread][(data+7)->iVar+NodeUnknownVariableVectorLength*(data+7)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[(data+6)->Thread][(data+6)->iVar+NodeUnknownVariableVectorLength*(data+6)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[(data+5)->Thread][(data+5)->iVar+NodeUnknownVariableVectorLength*(data+5)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[(data+4)->Thread][(data+4)->iVar+NodeUnknownVariableVectorLength*(data+4)->UnknownVectorIndex], 

                      LocalRecvExchangeBufferTable[(data+3)->Thread][(data+3)->iVar+NodeUnknownVariableVectorLength*(data+3)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[(data+2)->Thread][(data+2)->iVar+NodeUnknownVariableVectorLength*(data+2)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[(data+1)->Thread][(data+1)->iVar+NodeUnknownVariableVectorLength*(data+1)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex]
                      );

      cv=_mm512_mul_pd(av,bv);
      res+=_mm512_reduce_add_pd(cv);
    }

    //add the rest of the vector
    for (;iElement+3<iElementMax;iElement+=4) {
       alignas(32) double a[4],b[4],*r;
      __m256d av,bv,cv,rv;

      data=ElementDataTable+iElement;

      av=_mm256_set_pd((data+3)->MatrixElementValue,(data+2)->MatrixElementValue,(data+1)->MatrixElementValue,data->MatrixElementValue);

      bv=_mm256_set_pd(
                      LocalRecvExchangeBufferTable[(data+3)->Thread][(data+3)->iVar+NodeUnknownVariableVectorLength*(data+3)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[(data+2)->Thread][(data+2)->iVar+NodeUnknownVariableVectorLength*(data+2)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[(data+1)->Thread][(data+1)->iVar+NodeUnknownVariableVectorLength*(data+1)->UnknownVectorIndex],
                      LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex]
                      );

      cv=_mm256_mul_pd(av,bv);
      rv=_mm256_hadd_pd(cv,cv);

      r=(double*)&rv;
      res+=r[1]+r[2];
    }  
    #endif

    //add the rest of the vector
    for (;iElement<iElementMax;iElement++) {
      data=ElementDataTable+iElement;
      res+=data->MatrixElementValue*LocalRecvExchangeBufferTable[data->Thread][data->iVar+NodeUnknownVariableVectorLength*data->UnknownVectorIndex];
    }

     #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
     if (irow>=length) exit(__LINE__,__FILE__,"Error: out of bound");
     #endif

     p[irow]=res;
  }

  delete [] LocalRecvExchangeBufferTable;

}

  RecvExchangeBufferTable[PIC::ThisThread]=NULL;

  for (int thread=0;thread<PIC::nTotalThreads;thread++) if (RecvExchangeBufferTable[thread]!=NULL) { 
    delete [] RecvExchangeBufferTable[thread]; 
  }

  delete [] RecvExchangeBufferTable;

}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength,
int MaxRhsSupportLength_CornerNodes,int MaxRhsSupportLength_CenterNodes,
int MaxMatrixElementParameterTableLength,int MaxMatrixElementSupportTableLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength,MaxRhsSupportLength_CornerNodes,MaxRhsSupportLength_CenterNodes,MaxMatrixElementParameterTableLength,MaxMatrixElementSupportTableLength>::Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),
      void (*fUnpackSolution)(double* x,cCornerNode* CornerNode),
      double Tol, //the residual tolerance. The recommended value is 1e-5;
      int nMaxIter, //the max iteration error allowed. The recommended value is 100
      int (*fPackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,int* NodeDataLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* SendDataBuffer),
      int (*fUnpackBlockData)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** NodeTable,int NodeTableLength,unsigned char* BlockCenterNodeMask,unsigned char* BlockCornerNodeMask,char* RecvDataBuffer)
      ) {
  cMatrixRow* row;
  int cnt;

  //count the total number of the variables, and create the data buffers
  if (SubdomainPartialRHS==NULL) {
    for (SubdomainPartialUnknownsVectorLength=0,row=MatrixRowListFirst;row!=NULL;row=row->next) SubdomainPartialUnknownsVectorLength+=NodeUnknownVariableVectorLength;

    SubdomainPartialRHS=new double [SubdomainPartialUnknownsVectorLength];
    SubdomainPartialUnknownsVector=new double [SubdomainPartialUnknownsVectorLength];
  }

  //calculate the mean value of the Rhs that could be used later for the point located in the 'virtual' blocks
  int i;
  double MeanRhs[NodeUnknownVariableVectorLength];
  int MeanRhsCounter[NodeUnknownVariableVectorLength];

  for (i=0;i<NodeUnknownVariableVectorLength;i++) MeanRhs[i]=0.0,MeanRhsCounter[i]=0;

  for (row=MatrixRowListFirst;row!=NULL;cnt++,row=row->next) if (row->CornerNode!=NULL) {
    MeanRhs[row->iVar]+=row->Rhs;
    MeanRhsCounter[row->iVar]++;
  }

  for (i=0;i<NodeUnknownVariableVectorLength;i++) {
    if (MeanRhsCounter[i]==0) MeanRhs[i]=0.0;
    else MeanRhs[i]/=MeanRhsCounter[i];
  }


  //populate the buffers
  for (cnt=0,row=MatrixRowListFirst;row!=NULL;cnt++,row=row->next) {
    if (row->CornerNode!=NULL) {
      if (cnt%NodeUnknownVariableVectorLength==0) {
        fInitialUnknownValues(SubdomainPartialUnknownsVector+NodeUnknownVariableVectorLength*row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);
      }

      SubdomainPartialRHS[cnt]=row->Rhs;
    }
    else {
      //the point is located in the 'virtual' block
      SubdomainPartialUnknownsVector[row->iVar+NodeUnknownVariableVectorLength*row->ElementDataTable->UnknownVectorIndex]=MeanRhs[row->iVar];
      SubdomainPartialRHS[cnt]=MeanRhs[row->iVar];
    }
  }

  //call the iterative solver
  int nVar=NodeUnknownVariableVectorLength; //variable number
  int nDim = 3; //dimension
  int nI=_BLOCK_CELLS_X_;
  int nJ=_BLOCK_CELLS_Y_;
  int nK=_BLOCK_CELLS_Z_;
  int nBlock = nRealBlocks+nVirtualBlocks;

  MPI_Fint iComm = MPI_Comm_c2f(MPI_GLOBAL_COMMUNICATOR);

  //double Rhs_I(nVar*nVar*nI*nJ*nK*nBlock)// RHS of the equation
  //double Sol_I(nVar*nVar*nI*nJ*nK*nBlock)// vector for solution
  double * Rhs_I=SubdomainPartialRHS;
  double * Sol_I=SubdomainPartialUnknownsVector;
  double PrecondParam=0; // not use preconditioner
//  double ** precond_matrix_II;// pointer to precondition matrix; us  null if no preconditioner
  int lTest=(PIC::ThisThread==0);//1: need test output; 0: no test statement

  linear_solver_wrapper("GMRES", &Tol,&nMaxIter, &nVar, &nDim,&nI, &nJ, &nK, &nBlock, &iComm, Rhs_I,Sol_I, &PrecondParam, NULL, &lTest);

  //unpack the solution
  for (row=MatrixRowListFirst;row!=NULL;row=row->next) if (row->CornerNode!=NULL)  {
    fUnpackSolution(SubdomainPartialUnknownsVector+NodeUnknownVariableVectorLength*row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);
  }

  //execute the data exchange between 'ghost' blocks
  if (fPackBlockData!=NULL) {
    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh.ParallelBlockDataExchange(fPackBlockData,fUnpackBlockData);
      break;

    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::BC::ExternalBoundary::UpdateData(fPackBlockData,fUnpackBlockData);
      break;
    }
  }
  else {
    exit(__LINE__,__FILE__,"Error: fPackBlockData and fUnpackBlockData have to be defined");
  }

}

#endif /* _LINEARSYSTEMCORNERNODE_H_ */

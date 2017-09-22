//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf

//$Id$
//the AMR mesh in 1,2 and 3D without cut cells

#ifndef _AMR_MESH_GENERIC_ 
#define _AMR_MESH_GENERIC_ 

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>

#include <sys/time.h>
#include <sys/resource.h>


#include <string>
#include <iterator>
#include <utility>
#include <algorithm>
#include <list>
#include <vector>




#include "meshAMRdef.h"


//include the user defined data for the internal boundaries
/*
#if _AMR__LOAD_USER_DEFINITION__MODE_  == _AMR__LOAD_USER_DEFINITION__MODE__ON_
#include "UserDefinition.meshAMR.h"
#endif
*/

#include "meshAMR_UserDefinitions.h"


#include "meshAMRinternalSurface.h"
#include "meshAMRcutcell.h"
#include "meshNASTRAN.h"


#include "specfunc.h"
#include "mpichannel.h"


//define the boolian macro variables used in the mesh 
#define _AMR_FALSE_ 0
#define _AMR_TRUE_ 1




//the limits of the comlutational domain
extern double _MESH_AMR_XMAX_[3],_MESH_AMR_XMIN_[3];


class cBasicNode : public cAMRexit {
public:
  //the place holder for the structure that contained the associated data
  inline int AssociatedDataLength() {return 0;}
  void SetAssociatedDataBufferPointer(char* ptr) {}
  inline char* GetAssociatedDataBufferPointer() {return NULL;}

  //placeholder for the printing procedured
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int ThisThread) {}
  void PrintVariableList(FILE* fout,int DataSetNumber) {}
  void PrintFileDescriptior(FILE* fout,int DataSetNumber) {}

protected:
  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  double x[_MESH_DIMENSION_];
  #endif

public:
  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  struct cNodeDescriptor {
    unsigned long int nodeno : _MESH_ELEMENTS_NUMBERING_BITS_;
    unsigned nodeProcessedFlag : 1;
    unsigned maxRefinmentLevel : _MAX_REFINMENT_LEVEL_BITS_;
    unsigned internalMeshNode : 1;
    unsigned int nNodeConnections : _MAX_CORNER_NODE_CONNECTION_BITS_;
  } nodeDescriptor;


  //get and set the nodes' positions
  inline double *GetX() {
    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_OFF_
    exit(__LINE__,__FILE__,"The operation is allowed only in the debugger mode");
    #endif

    return x;
  }

  inline void GetX(double *l) {
    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_OFF_
    exit(__LINE__,__FILE__,"The operation is allowed only in the debugger mode");
    #endif

    for (int idim=0;idim<_MESH_DIMENSION_;idim++) l[idim]=x[idim];
  }

  inline void SetX(double *l) {


//########  DEBUG ##########
/*	  if (Temp_ID==88861) {
		  *DiagnospticMessageStream << __LINE__ << std::endl;
	  }
*/
//######## END DEBUG #######


    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_OFF_
    exit(__LINE__,__FILE__,"The operation is allowed only in the debugger mode");
    #endif

    for (int idim=0;idim<_MESH_DIMENSION_;idim++) x[idim]=l[idim];
  }


  //clean the data buffers
  void cleanDataBuffer() {
    nodeDescriptor.nodeno=0;
    nodeDescriptor.nodeProcessedFlag=_AMR_FALSE_;
    nodeDescriptor.maxRefinmentLevel=0;
    nodeDescriptor.internalMeshNode=_AMR_TRUE_;
    nodeDescriptor.nNodeConnections=0;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  //increment and decrement the node conenction's counter
  void incrementConnectionCounter() {
    if (nodeDescriptor.nNodeConnections==_MAX_CORNER_NODE_CONNECTION_) exit(__LINE__,__FILE__,"the node's connections exeeds _MAX_CORNER_NODE_CONNECTION_");

//    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_OFF_
    if (nodeDescriptor.nNodeConnections+1==(1<<_MAX_CORNER_NODE_CONNECTION_BITS_)) exit(__LINE__,__FILE__,"Error: the number of connections exeeds the limit. Increase '_MAX_CORNER_NODE_CONNECTION_BITS__MAX_CORNER_NODE_CONNECTION_BITS_'"); 
    nodeDescriptor.nNodeConnections++;
//    #endif
  }

  int getNodeConnectionNumber() {

//    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_OFF_
    return nodeDescriptor.nNodeConnections;
//    #else
//    return 1;
//    #endif
  }

  int decrementConnectionCounter() {

//    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_OFF_
    if (nodeDescriptor.nNodeConnections==0) exit(__LINE__,__FILE__,"the node's connections becomes negative");
    nodeDescriptor.nNodeConnections--;
//    #endif

    return getNodeConnectionNumber();
  }

}; 

class cBasicCornerNode : public cBasicNode {
public: 
  //double x[_MESH_DIMENSION_];
};


class cBasicCenterNode : public cBasicNode {
public: 
  double Measure;

  void cleanDataBuffer() {
    Measure=-1.0;

    cBasicNode::cleanDataBuffer();
  }

  //placeholder for the interpolation procedure
  void Interpolate(cBasicCenterNode **InterpolationStencil,double *InterpolationCoefficients,int StencilLength) {}

};



template <class cBlockAMR> class cTreeNodeAMR;
template <class cCornerNode,class cCenterNode> class cBasicBlockAMR; 

//===================================================================

#define _BOTTOM_BRANCH_TREE_ false

#define _GLOBAL_POSITION_REAL_NODE_UNKNOWN_ 0
#define _GLOBAL_POSITION_REAL_NODE_UP_      1
#define _GLOBAL_POSITION_REAL_NODE_CURRENT_ 2
#define _GLOBAL_POSITION_REAL_NODE_DOWN_    3

class cAMRnodeID {
public:
  int ResolutionLevel;
  unsigned char id[1+3*_MAX_REFINMENT_LEVEL_/8];

  cAMRnodeID() {
    ResolutionLevel=0;
    for (int i=0;i<1+3*_MAX_REFINMENT_LEVEL_/8;i++) id[i]=0;
  }

  bool operator == (cAMRnodeID ID) {
    if (ResolutionLevel!=ID.ResolutionLevel) return false;

    //compare the bit's fields
    int nbytes,i,nbits;

    nbits=3*ResolutionLevel;
    nbytes=nbits/8;
    nbits-=8*ResolutionLevel;

    //compare bytes
    for (i=0;i<nbytes;i++) if (id[i]!=ID.id[i]) return false;

    //compare bits
    unsigned char ComparisonMask=0;

//      for (i=0;i<nbits;i++) ComparisonMask+=(1<<i); substituted by the following
    for (i=0;i<nbits;i++) ComparisonMask=(unsigned char)(ComparisonMask|(1<<i));

    if ((id[nbytes]&ComparisonMask)!=(ID.id[nbytes]&ComparisonMask)) return false;

    return true;
  }

  bool operator != (cAMRnodeID ID) {
    return ((*this)==ID) ? false : true;
  }
};

template <class cBlockAMR>
class cTreeNodeAMR : public cAMRexit {
public:
  cTreeNodeAMR *upNode,*downNode[1<<_MESH_DIMENSION_];
  cBlockAMR *block;
  cAMRnodeID AMRnodeID;

  //Neighbors of the nodes
  #if _MESH_DIMENSION_ == 1
  cTreeNodeAMR *neibNodeFace[2];
  #elif _MESH_DIMENSION_ == 2
  cTreeNodeAMR *neibNodeFace[4*2],*neibNodeCorner[4];
  #elif _MESH_DIMENSION_ == 3
  cTreeNodeAMR *neibNodeFace[6*4],*neibNodeCorner[8],*neibNodeEdge[12*2];
  #endif


  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  bool ActiveFlag; //used to prevent repeatable de-allocation of the block from the stack

  //descriptor of the cut-face list
  struct cCutFaceListDescriptor {
    cCutFaceListDescriptor* next;
    cTreeNodeAMR<cBlockAMR>* node;
    CutCell::cTriangleFaceDescriptor* FirstTriangleCutFace;
  };

  cCutFaceListDescriptor* neibCutFaceListDescriptorList;
  cCutFaceListDescriptor*  neibCutFaceListDescriptorList_temp;

  cTreeNodeAMR() {
    ActiveFlag=false;
    block=NULL,upNode=NULL;
    Temp_ID=-1;
    neibCutFaceListDescriptorList=NULL,neibCutFaceListDescriptorList_temp=NULL;

    for (int i=0;i<1<<_MESH_DIMENSION_;i++) downNode[i]=NULL;

    #if _MESH_DIMENSION_ == 1
    for (int i=0;i<2;i++) neibNodeFace[i]=NULL;

    #elif _MESH_DIMENSION_ == 2
    for (int i=0;i<4*2;i++) neibNodeFace[i]=NULL;
    for (int i=0;i<4;i++) neibNodeCorner[i]=NULL;

    #elif _MESH_DIMENSION_ == 3
    for (int i=0;i<6*4;i++) neibNodeFace[i]=NULL;
    for (int i=0;i<8;i++) neibNodeCorner[i]=NULL;
    for (int i=0;i<12*2;i++) neibNodeEdge[i]=NULL;
    #endif //_MESH_DIMENSION_
  }


  //the pointers to members of the list of the nodes that "belongs" to this processor (or some other lists)
  cTreeNodeAMR *nextNodeThisThread,*prevNodeThisThread;

  //the pointers to the previous and the next nodes in the list of node's that has intersection with the surfaces of the computational domain
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  cTreeNodeAMR *DomainSurfaceBoundaryList_Prev, *DomainSurfaceBoundaryList_Next;
  #endif


  //the list of the cut-face descriptors
#if _AMR__CUT_CELL__MODE_ ==  _AMR__CUT_CELL__MODE__ON_
  CutCell::cTriangleFaceDescriptor *FirstTriangleCutFace,*neibFirstTriangleCutFace,*neibFirstTriangleCutFace_temp;
  double (*CutCellSurfaceLocalResolution)(CutCell::cTriangleFaceDescriptor*);
#endif

  double xmin[_MESH_DIMENSION_],xmax[_MESH_DIMENSION_];
  int RefinmentLevel,minNeibRefinmentLevel,maxNeibRefinmentLevel; //min/max values are used for calculation of the interpolation stencils

  struct cNodeDescriptor {
    unsigned NodeProcessingFlag : 1;
  } nodeDescriptor;

  //the list of the descriptors of the internal boundaries installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptorList;
  #endif

  //the data uswed for the parallel implementation of the mesh; next node in the space filling curve
  #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
  int Thread;
  double ParallelLoadMeasure;
  cTreeNodeAMR *FillingCurveNextNode,*FillingCurvePrevNode;
  #endif


  //calculate the position of the corner nodes
  void GetCornerNodePosition(double *xProbe,int i,int j,int k)  {
    xProbe[0]=xmin[0]+i*((xmax[0]-xmin[0])/_BLOCK_CELLS_X_);
    if (_MESH_DIMENSION_>1) xProbe[1]=xmin[1]+j*((xmax[1]-xmin[1])/_BLOCK_CELLS_Y_);
    if (_MESH_DIMENSION_>2) xProbe[2]=xmin[2]+k*((xmax[2]-xmin[2])/_BLOCK_CELLS_Z_);
  }

  //get the external (directed outward of the computational domain) normal
  void GetExternalNormal(double *norm,int nface) {
    switch(nface) {
    case 0: case 1:
      norm[0]=(nface==0) ? -1.0 : 1.0;
      norm[1]=0.0,norm[2]=0.0;
      break;

    case 2: case 3:
      norm[1]=(nface==2) ? -1.0 : 1.0;
      norm[0]=0.0,norm[2]=0.0;
      break;

    case 4: case 5:
      norm[2]=(nface==4) ? -1.0 : 1.0;
      norm[0]=0.0,norm[1]=0.0;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: 'nface' is out of the range");
    }
  }

  double GetCharacteristicCellSize() {
    double CellSize;

    CellSize=pow((xmax[0]-xmin[0])/_BLOCK_CELLS_X_,2);
    if (_MESH_DIMENSION_>1) CellSize+=pow((xmax[1]-xmin[1])/_BLOCK_CELLS_Y_,2);
    if (_MESH_DIMENSION_>2) CellSize+=pow((xmax[2]-xmin[2])/_BLOCK_CELLS_Z_,2);

    return sqrt(CellSize);
  }

  void ConvertGlobal2LocalCoordinates(double *LocalCoordinates,double *GlobalCoordinates) {
    int idim;

    for (idim=0;idim<_MESH_DIMENSION_;idim++) LocalCoordinates[idim]=(GlobalCoordinates[idim]-xmin[idim])/(xmax[idim]-xmin[idim]);
  }

  double GetBlockFaceSurfaceArea(int nface) {
    double res;

    if (nface>=2*_MESH_DIMENSION_) exit(__LINE__,__FILE__,"Error: 'nface' is out of range");

    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 3

    static const int faceNodeCoordinateFrame[6][3][3]={ {{0,0,0},{0,1,0},{0,0,1}},{{1,0,0},{1,1,0},{1,0,1}},
                                                        {{0,0,0},{1,0,0},{0,0,1}},{{0,1,0},{1,1,0},{0,1,1}},
                                                        {{0,0,0},{1,0,0},{0,1,0}},{{0,0,1},{1,0,1},{0,1,1}} };



    double x0[3],x1[3],x2[3],l0,l1;
    int idim;

    x0[0]=(faceNodeCoordinateFrame[nface][0][0]==0) ? xmin[0] : xmax[0];
    x0[1]=(faceNodeCoordinateFrame[nface][0][1]==0) ? xmin[1] : xmax[1];
    x0[2]=(faceNodeCoordinateFrame[nface][0][2]==0) ? xmin[2] : xmax[2];

    x1[0]=(faceNodeCoordinateFrame[nface][1][0]==0) ? xmin[0] : xmax[0];
    x1[1]=(faceNodeCoordinateFrame[nface][1][1]==0) ? xmin[1] : xmax[1];
    x1[2]=(faceNodeCoordinateFrame[nface][1][2]==0) ? xmin[2] : xmax[2];

    x2[0]=(faceNodeCoordinateFrame[nface][2][0]==0) ? xmin[0] : xmax[0];
    x2[1]=(faceNodeCoordinateFrame[nface][2][1]==0) ? xmin[1] : xmax[1];
    x2[2]=(faceNodeCoordinateFrame[nface][2][2]==0) ? xmin[2] : xmax[2];

    for (idim=0,l0=0.0,l1=0.0;idim<_MESH_DIMENSION_;idim++) l0+=pow(x1[idim]-x0[idim],2),l1+=pow(x2[idim]-x0[idim],2);
    res=sqrt(l0*l1);

    #endif

    return res;
  }

  double GetCellFaceSurfaceArea(int nface) {
    #if _MESH_DIMENSION_ == 1
    const static double CellSurfaceAreaMultiplyer[1]={1.0,1.0};
    #elif _MESH_DIMENSION_ == 2
    const static double CellSurfaceAreaMultiplyer[4]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,
                                                      _BLOCK_CELLS_X_,_BLOCK_CELLS_X_};
    #elif _MESH_DIMENSION_ == 3

    const static double CellSurfaceAreaMultiplyer[6]={_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,
                                                      _BLOCK_CELLS_X_*_BLOCK_CELLS_Z_,_BLOCK_CELLS_X_*_BLOCK_CELLS_Z_,
                                                      _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_};
    #endif

    return GetBlockFaceSurfaceArea(nface)/CellSurfaceAreaMultiplyer[nface];
  }

  struct cTreeNodeDescriptor {
    unsigned int GlobalPositionRealTreeNode : 2;
  } treeNodeDescriptor; 


  void cleanDataBuffer() {
    int i;

    block=NULL,upNode=NULL;

    for (i=0;i<(1<<_MESH_DIMENSION_);i++) downNode[i]=NULL;
  
    treeNodeDescriptor.GlobalPositionRealTreeNode=_GLOBAL_POSITION_REAL_NODE_UNKNOWN_; 

    RefinmentLevel=-1,minNeibRefinmentLevel=-1,maxNeibRefinmentLevel=-1;
    nextNodeThisThread=NULL,prevNodeThisThread=NULL;

    nodeDescriptor.NodeProcessingFlag=_AMR_FALSE_;

#if _AMR__CUT_CELL__MODE_ ==  _AMR__CUT_CELL__MODE__ON_
    FirstTriangleCutFace=NULL,neibFirstTriangleCutFace=NULL,neibFirstTriangleCutFace_temp=NULL;
#endif

    //Neighbors of the nodes
    #if _MESH_DIMENSION_ == 1
    for (i=0;i<2;i++) neibNodeFace[i]=NULL;
    #elif _MESH_DIMENSION_ == 2
    for (i=0;i<4*2;i++) neibNodeFace[i]=NULL;
    for (i=0;i<4;i++) neibNodeCorner[i]=NULL;
    #elif _MESH_DIMENSION_ == 3
    for (i=0;i<6*4;i++) neibNodeFace[i]=NULL;
    for (i=0;i<8;i++) neibNodeCorner[i]=NULL;
    for (i=0;i<12*2;i++) neibNodeEdge[i]=NULL;
    #endif


    #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
    FillingCurveNextNode=NULL;
    #endif

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif

    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    InternalBoundaryDescriptorList=NULL;
    DomainSurfaceBoundaryList_Prev=NULL,DomainSurfaceBoundaryList_Next=NULL;
    #endif

    #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
    Thread=0,ParallelLoadMeasure=0;
    #endif
  }

  //determine the minimum and maximum resolution levels of the neighbor blocks
  void SetNeibRefinmentLevelLimits() {
    cTreeNodeAMR* node;
    int i;

    minNeibRefinmentLevel=-1,maxNeibRefinmentLevel=-1;

    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 3
    //connections through faces
    for (i=0;i<6*4;i++) if ((node=neibNodeFace[i])!=NULL) {
      if ((minNeibRefinmentLevel==-1)||(minNeibRefinmentLevel>node->RefinmentLevel)) minNeibRefinmentLevel=node->RefinmentLevel;
      if (maxNeibRefinmentLevel<node->RefinmentLevel) maxNeibRefinmentLevel=node->RefinmentLevel;
    }

    //connection through corners
    for (i=0;i<8;i++) if ((node=neibNodeCorner[i])!=NULL) {
      if ((minNeibRefinmentLevel==-1)||(minNeibRefinmentLevel>node->RefinmentLevel)) minNeibRefinmentLevel=node->RefinmentLevel;
      if (maxNeibRefinmentLevel<node->RefinmentLevel) maxNeibRefinmentLevel=node->RefinmentLevel;
    }

    //connection through edges
    for (i=0;i<12*2;i++) if ((node=neibNodeEdge[i])!=NULL) {
      if ((minNeibRefinmentLevel==-1)||(minNeibRefinmentLevel>node->RefinmentLevel)) minNeibRefinmentLevel=node->RefinmentLevel;
      if (maxNeibRefinmentLevel<node->RefinmentLevel) maxNeibRefinmentLevel=node->RefinmentLevel;
    }
    #endif
  }


  //find the neighbor of the tree node: this version of the function is used only in 1D case
  /*
  cTreeNodeAMR *GetNeib(int i) {
    cTreeNodeAMR *res=NULL;

    #if _MESH_DIMENSION_ == 1 
    res=neibNode[i];
    #else
    exit(__LINE__,__FILE__,"The wrong dimension");
    #endif
  }

  cTreeNodeAMR *GetNeib(int i,int j,int k) {
    #if _MESH_DIMENSION_ == 1 
    return GetNeib(i);
    #else 
    exit(__LINE__,__FILE__,"not implemented yet");
    #endif
  }
  */


  cTreeNodeAMR *GetNeibCorner(int nCornerNode) {
    cTreeNodeAMR *res;

#if _MESH_DIMENSION_ == 1
//    exit(__LINE__,__FILE__,"Error: out of range");
    res=NULL;
#else
    if ((nCornerNode<0)||(nCornerNode>=1<<_MESH_DIMENSION_)) exit(__LINE__,__FILE__,"Error: out of range");
    res=neibNodeCorner[nCornerNode];
#endif

    return res;
  }

  void SetNeibCorner(cTreeNodeAMR* neibNode,int nCornerNode) {



//################## DEBUG #####################
/*    if ((this->Temp_ID==222)||(this->Temp_ID==107)) {
      *DiagnospticMessageStream << __LINE__ << std::endl;


    }
*/
//################## END DEBUG #####################



#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    //check if the neibNode is an upper branch
    if (neibNode!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) exit(__LINE__,__FILE__,"Error: neibNode is not a bottom tree node");
#endif


#if _MESH_DIMENSION_ == 1
//    exit(__LINE__,__FILE__,"Error: out of range");
// do nothing
#else
    if ((nCornerNode<0)||(nCornerNode>=1<<_MESH_DIMENSION_)) exit(__LINE__,__FILE__,"Error: out of range");
    neibNodeCorner[nCornerNode]=neibNode;
#endif
  }

  inline cTreeNodeAMR *GetNeibFace(int nface,int iFace,int jFace) {
    cTreeNodeAMR* res;

#if _MESH_DIMENSION_ == 1
    res=neibNodeFace[nface];
#elif _MESH_DIMENSION_ == 2
    res=neibNodeFace[iFace+2*nface];
#else
    res=neibNodeFace[iFace+2*(jFace+2*nface)];
#endif

    return res;
  }

  inline void SetNeibFace(cTreeNodeAMR* neibNode,int nface,int iFace,int jFace) {

#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    //check if the neibNode is an upper branch
    if (neibNode!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) exit(__LINE__,__FILE__,"Error: neibNode is not a bottom tree node");
#endif

#if _MESH_DIMENSION_ == 1
    neibNodeFace[nface]=neibNode;
#elif _MESH_DIMENSION_ == 2
    neibNodeFace[iFace+2*nface]=neibNode;
#else
    neibNodeFace[iFace+2*(jFace+2*nface)]=neibNode;
#endif
  }

  inline cTreeNodeAMR *GetNeibEdge(int nedge,int iEdge) {
    cTreeNodeAMR *res;

#if _MESH_DIMENSION_ == 1
    res=NULL;
#elif _MESH_DIMENSION_ == 2
    res=NULL;
#elif _MESH_DIMENSION_ == 3
    res=neibNodeEdge[iEdge+2*nedge];
#endif

    return res;
  }

  inline void SetNeibEdge(cTreeNodeAMR* neibNode,int nedge,int iEdge) {

#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    //check if the neibNode is an upper branch
    if (neibNode!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) exit(__LINE__,__FILE__,"Error: neibNode is not a bottom tree node");
#endif

#if _MESH_DIMENSION_ == 1
    //do nothing
#elif _MESH_DIMENSION_ == 2
    //do nothing
#elif _MESH_DIMENSION_ == 3
    neibNodeEdge[iEdge+2*nedge]=neibNode;
#endif
  }

  inline cTreeNodeAMR *GetNeibNode(int i,int j,int k) {
    cTreeNodeAMR *res=NULL;

#if _MESH_DIMENSION_ == 1
    if (i==-1) res=neibNodeFace[0];
    else if (i==1) res=neibNodeFace[1];
    else if (i!=0) exit(__LINE__,__FILE__,"Error: out of range");

#elif _MESH_DIMENSION_ == 2
     //determine weather the neibNode is a face or a corner node
     register int code=abs(i)+abs(j);

     if (code==2) { //corner node
       if (i==-1) res=(j==-1) ? neibNodeCorner[0] : neibNodeCorner[2];
       else if (i==1) res=(j==-1) ? neibNodeCorner[1] : neibNodeCorner[3];
       else exit(__LINE__,__FILE__,"Error: out of range");
     }
     else if (code==1) { //face node
       if (i==-1) res=neibNodeFace[0*2+0];
       else if (i==1) res=neibNodeFace[1*2+0];
       else if (j==-1) res=neibNodeFace[2*2+0];
       else if (j==1) res=neibNodeFace[3*2+0];
       else exit(__LINE__,__FILE__,"Error: out of range");
     }
     else if (code==0) res=this;
     else exit(__LINE__,__FILE__,"Error: parameters are out or range");

#elif _MESH_DIMENSION_ == 3
     //determine weather the neibNode is a face or a corner node
     register int code=abs(i)+abs(j)+abs(k);

     if (code==3) { //corner node
       if (i==-1) i=0;
       if (j==-1) j=0;
       if (k==-1) k=0;

       res=GetNeibCorner(i+2*(j+2*k));
     } else if (code==2) { //edge node
       int nedge=-1;

       if (i==-1) {
         if (j!=0) nedge=(j==-1) ? 8 : 11;
         else nedge=(k==-1) ? 4 : 7;
       }
       else if (i==1) {
         if (j!=0) nedge=(j==-1) ? 9 : 10;
         else if (k!=0) nedge=(k==-1) ? 5 : 6;
       }
       else if (j==-1) nedge=(k==-1) ? 0 : 3;
       else if (j==1) nedge=(k==-1) ? 1 : 2;
       else exit(__LINE__,__FILE__,"Error: out of range");

       res=GetNeibEdge(nedge,0);
     } else if (code==1) { //face node
       int nface=-1;

       if (i==-1) nface=0;
       else if (i==1) nface=1;
       else if (j==-1) nface=2;
       else if (j==1) nface=3;
       else if (k==-1) nface=4;
       else if (k==1) nface=5;
       else exit(__LINE__,__FILE__,"Error: out of range");

       res=GetNeibFace(nface,0,0);

     } else if (code==0) res=this;
     else exit(__LINE__,__FILE__,"Error: parameters are out or range");
#else
     exit(__LINE__,__FILE__,"Error: wrong option");
#endif

     return res;
  }


  bool allocated() {return (block==NULL) ? false : true;};

  int GetGlobalPositinoRealNode() {return treeNodeDescriptor.GlobalPositionRealTreeNode;}
  void SettGlobalPositinoRealNode(int RealNodePosition) {treeNodeDescriptor.GlobalPositionRealTreeNode=RealNodePosition;} 


  cTreeNodeAMR *GetDownNode(int i,int j,int k) {return downNode[i+2*(j+k*2)];} 
  void SetDownNode(cTreeNodeAMR* node,int i,int j,int k) {downNode[i+2*(j+k*2)]=node;}  

  bool lastBranchFlag() {
    //if _BOTTOM_BRANCH_TREE_ -> the node is on the bottom of the tree
    for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (downNode[nDownNode]!=NULL) return !_BOTTOM_BRANCH_TREE_;

    return _BOTTOM_BRANCH_TREE_;
  }

  //check if a node is a neighbor of the 'this node'
  bool CheckNeibNode(cTreeNodeAMR* neibNode) {
    //check the connections through the corner nodes
    for (register int nd=0;nd<(1<<_MESH_DIMENSION_);nd++) if (GetNeibCorner(nd)==neibNode) {
      return true;
    }

    //check the connection through the faces
#if _MESH_DIMENSION_ == 2
    register int nface,i;

    for (nface=0;nface<4;nface++) for (i=0;i<2;i++) if (GetNeibFace(nface,i,0)==neibNode) {
      return true;
    }
#elif _MESH_DIMENSION_ == 3
    register int nface,i,j;

    for (nface=0;nface<6;nface++) for (i=0;i<2;i++) for (j=0;j<2;j++) if (GetNeibFace(nface,i,j)==neibNode) {
      return true;
    }
#endif


    //check connections through the edges
#if _MESH_DIMENSION_ == 3
    for (register int nedge=0;nedge<12;nedge++) for (i=0;i<2;i++) if (GetNeibEdge(nedge,i)==neibNode) {
      return true;
    }
#endif

    return false;
  }



};





//=======================================================================
template <class cCornerNode,class cCenterNode> 
class cBasicBlockAMR : public cAMRexit {
public:
  //the place holder for the structure that contained the associated data
  int AssociatedDataLength() {return 0;}
  void SetAssociatedDataBufferPointer(char* ptr) {}
  char* GetAssociatedDataBufferPointer() {return NULL;}

  //place holder for print function that outputs the general tree node/blocks values into a file
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int ThisThread) {}
  void PrintVariableList(FILE* fout) {}

protected:
  //nodes of the block

  #if _MESH_DIMENSION_ == 1
  cCornerNode *cornerNodes[1+_TOTAL_BLOCK_CELLS_X_];

  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  cCenterNode *centerNodes[_TOTAL_BLOCK_CELLS_X_];
  #endif

  #elif _MESH_DIMENSION_ == 2
  cCornerNode *cornerNodes[(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)];

  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  cCenterNode *centerNodes[_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_];
  #endif

  #elif _MESH_DIMENSION_ == 3
  cCornerNode *cornerNodes[(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_)]; 

  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  cCenterNode *centerNodes[_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_];
  #endif


  #endif

public:

  //set and get the pointers to the center nodes of the block
  cCenterNode **GetCenterNodeBuffer() {
    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    return centerNodes;
    #else
    return NULL;
    #endif
  }

  inline cCenterNode *GetCenterNode(long int nd) {

    #if _MESH_DIMENSION_ == 1
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_;
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_;
    #endif

    if ((nd<0)||(nd>=nMaxCenterNodes)) {
      printf("$PREFIX:nd=%ld\n",nd);
      exit(__LINE__,__FILE__,"The value is outside of the limit");
    }

    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    return centerNodes[nd];
    #else
    return NULL;
    #endif
  }


  inline void SetCenterNode(cCenterNode* nodeptr,long int nd) {

    #if _MESH_DIMENSION_ == 1
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_;
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_;
    #endif

    if ((nd<0)||(nd>=nMaxCenterNodes)) exit(__LINE__,__FILE__,"The value is outside of the limit");

    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    centerNodes[nd]=nodeptr;
    #endif
  }

  inline cCenterNode **GetCenterNodeList() {return centerNodes;}

  inline int GetCenterNodeListLength() {
    #if _MESH_DIMENSION_ == 1
    return _TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    return _TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_;
    #elif _MESH_DIMENSION_ == 3
    return _TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_;
    #endif
  }

  //set and get the pointers to the corner nodes of the block
  cCornerNode **GetCornerNodeBuffer() {
	 return cornerNodes;
  }

  inline cCornerNode *GetCornerNode(long int nd) {

    #if _MESH_DIMENSION_ == 1
    static const int nMaxCornerNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    if ((nd<0)||(nd>=nMaxCornerNodes)) exit(__LINE__,__FILE__,"The value is outside of the limit");

    return cornerNodes[nd];
  }


  inline void SetCornerNode(cCornerNode* nodeptr,long int nd) {
    #if _MESH_DIMENSION_ == 1
    static const int nMaxCornerNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    if ((nd<0)||(nd>=nMaxCornerNodes)) exit(__LINE__,__FILE__,"The value is outside of the limit");

    cornerNodes[nd]=nodeptr;
  }

  inline cCornerNode **GetCornerNodeList() {return cornerNodes;}

  inline static int GetCornerNodeListLength() {
    #if _MESH_DIMENSION_ == 1
    return 1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    return (1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    return (1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif
  }

  inline static long int getCornerNodeLocalNumber(int i,int j,int k) {
    long int nd;

    #if _MESH_DIMENSION_ == 1
    nd=i+_GHOST_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nd=i+_GHOST_CELLS_X_+(1+_TOTAL_BLOCK_CELLS_X_)*(j+_GHOST_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    nd=i+_GHOST_CELLS_X_+(1+_TOTAL_BLOCK_CELLS_X_)*(j+_GHOST_CELLS_Y_+(k+_GHOST_CELLS_Z_)*(1+_TOTAL_BLOCK_CELLS_Y_));
    #endif

    return nd;
  }

  inline static long int getCenterNodeLocalNumber(int i,int j,int k) {
    long int nd;

    #if _MESH_DIMENSION_ == 1
    nd=i+_GHOST_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nd=i+_GHOST_CELLS_X_+_TOTAL_BLOCK_CELLS_X_*(j+_GHOST_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    nd=i+_GHOST_CELLS_X_+_TOTAL_BLOCK_CELLS_X_*(j+_GHOST_CELLS_Y_+(k+_GHOST_CELLS_Z_)*_TOTAL_BLOCK_CELLS_Y_);
    #endif

    return nd;
  }

  inline void getCenterNodeCoordinate(int& i,int& j,int& k,long int nd) {

    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 3
    k=nd/(_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_);
    j=(nd-k*_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_)/_TOTAL_BLOCK_CELLS_X_;
    i=nd-k*_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_-j*_TOTAL_BLOCK_CELLS_X_;

    i-=_GHOST_CELLS_X_;
    j-=_GHOST_CELLS_Y_;
    k-=_GHOST_CELLS_Z_;
    #endif
  }

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

//  void *nextBlock; //the next block in the list of allocated blocks of the mesh


  struct cBlockDescriptor {
    unsigned RefinmentLevel : _MAX_REFINMENT_LEVEL_BITS_;
    unsigned ghostBlock : 1;
  } blockDescriptor; 

  double CellMeasure(int i,int j, int k) {
    register double vol;

    vol=(_MESH_AMR_XMAX_[0]-_MESH_AMR_XMIN_[0])*(_MESH_AMR_XMAX_[1]-_MESH_AMR_XMIN_[1])*(_MESH_AMR_XMAX_[2]-_MESH_AMR_XMIN_[2]);

    #if _MESH_DIMENSION_ == 1
    vol/=_TOTAL_BLOCK_CELLS_X_*(1<<blockDescriptor.RefinmentLevel);
    #elif _MESH_DIMENSION_ == 2
    vol/=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*(1<<(1+blockDescriptor.RefinmentLevel));
    #elif _MESH_DIMENSION_ == 3
    vol/=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_*(1<<(2+blockDescriptor.RefinmentLevel));
    #endif

    return vol; 
  } 

  double CellCharacteristicSize() {
    double res;

    res=pow((_MESH_AMR_XMAX_[0]-_MESH_AMR_XMIN_[0])/_BLOCK_CELLS_X_,2);
    if (_MESH_DIMENSION_>1) res+=pow((_MESH_AMR_XMAX_[1]-_MESH_AMR_XMIN_[1])/_BLOCK_CELLS_Y_,2);
    if (_MESH_DIMENSION_>2) res+=pow((_MESH_AMR_XMAX_[2]-_MESH_AMR_XMIN_[2])/_BLOCK_CELLS_Z_,2);

    res=sqrt(res)/(1<<blockDescriptor.RefinmentLevel);

    return res;
  }


  void cleanDataBuffer() {
    blockDescriptor.RefinmentLevel=0;
    blockDescriptor.ghostBlock=0;

//    nextBlock=NULL;

#if _MESH_DIMENSION_ == 1
    for (int i=0;i<1+_TOTAL_BLOCK_CELLS_X_;i++) cornerNodes[i]=NULL;

   #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    for (int i=0;i<_TOTAL_BLOCK_CELLS_X_;i++) centerNodes[i]=NULL;
   #endif

#elif _MESH_DIMENSION_ == 2
   for (int i=0;i<(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);i++) cornerNodes[i]=NULL;

   #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
   for (int i=0;i<_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_;i++) centerNodes[i]=NULL;
   #endif

#elif _MESH_DIMENSION_ == 3
   for (int i=0;i<(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);i++) cornerNodes[i]=NULL;

  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
   for (int i=0;i<_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_;i++) centerNodes[i]=NULL;
  #endif
#endif

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  } 

  cBasicBlockAMR() {
    cleanDataBuffer();
  }


  bool GetGhostFlag() {return (blockDescriptor.ghostBlock==0) ? _REAL_BLOCK_ : _GHOST_BLOCK_;}
  void SetGhostFlag(bool flag) {blockDescriptor.ghostBlock=(flag==_REAL_BLOCK_) ? 0 : 1;} //if false -> the block is real, if true -> the block is a ghost block 

  int GetRefinmentLevel() {return blockDescriptor.RefinmentLevel;}
  void SetRefinmentLevel(int level) {blockDescriptor.RefinmentLevel=level;}

};



  

//======================================================================================
template <class cCornerNode,class cCenterNode,class cBlockAMR>
class cMeshAMRgeneric : public cAMRexit {
public:
  cTreeNodeAMR<cBlockAMR>  *rootTree;

  //the stream for output of the mesh diagnisotic information
  FILE *DiagnospticMessageStream;

  //parallel mesh generation flag;
  bool ParallelMeshGenerationFlag;

  //the function that calculates the interpolation coefficients to get an interpolated value for the block's nodes
  //return the number of the interpolation coefficients that was used in the stencil. if the return value <=0 -> the operation is not succesful
  typedef int (*cGetCornerNodesInterpolationCoefficients)(double *x,double *CoefficientsList,cCornerNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients);
  cGetCornerNodesInterpolationCoefficients GetCornerNodesInterpolationCoefficients;

  typedef int (*cGetCenterNodesInterpolationCoefficients)(double *x,double *CoefficientsList,cCenterNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients);
  cGetCenterNodesInterpolationCoefficients GetCenterNodesInterpolationCoefficients;

  //accept tree node function
  typedef bool (*cAcceptBlockFunc)(double*,double*);
  cAcceptBlockFunc accepltTreeNodeFunction;
  void SetAcceptBlockFunction(cAcceptBlockFunc t) {accepltTreeNodeFunction=t;}

  //the counter of any mesh modifications or rebalancing 
  unsigned long int nMeshModificationCounter; 

  //limit of the accuracy of calculation of the mesh parameters
  double EPS;

  //the name and fignature of the mesh: the name contains the data when the mesh has been created; all files associated with the mesh will starts with the mesh name
  char MeshName[STRING_LENGTH];
  unsigned long int MeshSignature;

  //the total number of processores and the current processor
  #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
  int nTotalThreads,ThisThread;
  bool **ParallelSendRecvMap;
  #endif

  //the Local mesh resolution at the surface defined by the cut-cells
  #if _AMR__CUT_CELL__MODE_ ==  _AMR__CUT_CELL__MODE__ON_
  double (*CutCellSurfaceLocalResolution)(CutCell::cTriangleFace*);
  #endif

  //for testing the parallel mesh generation
  int nMpiBarrierCalls;

  void callMpiBarrier() {
    nMpiBarrierCalls++;
//    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }


  //the internal surface triangulation (NASTRAN MESH)
  cNASTRANmesh NsatranSurfaceMesh;

  //mesh modyfied flag -> is set to true each time the mesh is modified, the flag is set to false when the number of elements and the connectivity list are prepared
  bool meshModifiedFlag;

  //when the mesh is modified -> create the new space filling curve
  bool meshModifiedFlag_CreateNewSpaceFillingCurve;
  bool meshModifiedFlag_CountMeshElements;

  //deallocate blocks of the nodes that are not on the  bottom of the tree
  bool DeallocateUnusedBlocks;

  //allow automatic allocation of blocks during the mesh generation (if AllowBlockAllocation==false -> only the tree will be generated);:/
  bool AllowBlockAllocation;

  //the total number of mesh nodes and blocks 
  long int meshNodesNumber,meshBlocksNumber;

  //the upper (maximum) refinement level used on the mesh
  int meshMaximumRefinmentLevel;

  //the limits of the mesh
  double xGlobalMin[_MESH_DIMENSION_],xGlobalMax[_MESH_DIMENSION_];

  //the user defined function for the local mesh resolution 
  double (*localResolution)(double*);

  //the stacks for the tree, blocks and nodes
  cAssociatedDataAMRstack<cCenterNode> CenterNodes;
  cAssociatedDataAMRstack<cCornerNode> CornerNodes;
  cAssociatedDataAMRstack<cBlockAMR> blocks;
  cAMRstack<cTreeNodeAMR<cBlockAMR> > treeNodes;

  //the stack containing the descriptors for all internal surfaces installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  cAMRstack<cInternalBoundaryConditionsDescriptor> InternalBoundaryDescriptors;
  #endif



  //the root block;
  cBlockAMR *rootBlock;
  double dxRootBlock[3];

  //the list of blocks of the mesh
//  void *blockList;

  //the list of the nodes that "belongs" to each processor
  cTreeNodeAMR<cBlockAMR> **ParallelNodesDistributionList;

  #if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
  cTreeNodeAMR<cBlockAMR> **DomainBoundaryLayerNodesList;
  #endif

  //the list of the descriptors of all internal boundaries installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  list<cInternalBoundaryConditionsDescriptor> InternalBoundaryList;
  cTreeNodeAMR<cBlockAMR>* DomainSurfaceBoundaryList;
  #endif

  //calculate the checksumm of the distribution of the tree
  class cTreeCheckSum {
    CRC32 checksum;

    void addblock(cTreeNodeAMR<cBlockAMR>* bl) {
      checksum.add(bl->Temp_ID);

      for (int idim=0;idim<_MESH_DIMENSION_;idim++) {
        checksum.add(bl->xmax[idim]);
        checksum.add(bl->xmin[idim]);
      }

      for (int nd=0;nd<(1<<_MESH_DIMENSION_);nd++) if (bl->downNode[nd]!=NULL) addblock(bl->downNode[nd]);
    }

  public:
    unsigned long Get(cTreeNodeAMR<cBlockAMR>* StartNode,char const *msg,bool printflag) {
      unsigned long res;
      char message[STRING_LENGTH];
      int rank;

      checksum.clear();
      addblock(StartNode);

      MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);

      if (printflag==true) {
        if (msg!=NULL) sprintf(message,", %s[rank=%i,l=%ld,f=%s]",msg,rank,__LINE__,__FILE__);
        else sprintf(message,"rank=%i,l=%ld,f=%s",rank,__LINE__,__FILE__);

        printf("AMR tree checksum:\n%s\nckecksum=0x%lx\n",message,checksum.checksum());
      }

      return checksum.checksum();
    }

  } TreeCheckSum;


  //generate the mesh signeture: the signature contained the time of the mesh creeation, the user name and the computer name where the lesh is created
  void generateMeshSignature(cTreeNodeAMR<cBlockAMR> *node=NULL) {
    static CRC32 Signature;
    static CMPI_channel pipe;

    if (node==NULL) {
      node=rootTree;
      Signature.clear();
      pipe.init(1000000);

      //add the size of the domain to the signature
      Signature.add(rootTree->xmax,3);
      Signature.add(rootTree->xmin,3);

      //add to the signature the number of the cells and ghost cells
      int t[2];

      t[0]=_GHOST_CELLS_X_,t[1]=_BLOCK_CELLS_X_;
      Signature.add(t,2);

      if (_MESH_DIMENSION_>=2) {
        t[0]=_GHOST_CELLS_Y_,t[1]=_BLOCK_CELLS_Y_;
        Signature.add(t,2);
      }

      if (_MESH_DIMENSION_==3) {
        t[0]=_GHOST_CELLS_Z_,t[1]=_BLOCK_CELLS_Z_;
        Signature.add(t,2);
      }

      //open mpi channel
      if (ThisThread==0) pipe.openRecvAll();
      else pipe.openSend(0);
    }

    //scan through the tree and calculate the signature component the is due to the distribution of the cells
    //determine whether the node is at the bottom of the tree
    int i,j,k,nd;

    if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      //add the local number of the allocated center nodes

      #if _MESH_DIMENSION_ == 3
      const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
      const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
      const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
      #else
      exit(__LINE__,__FILE__,"Error: not implemented");
      #endif


      if ((ThisThread==node->Thread)||(ThisThread==0)) {
        int BlockAllocatedFlag;

        if (ThisThread==node->Thread) {
          BlockAllocatedFlag=(node->block==NULL) ? false : true;
          if (ThisThread!=0) pipe.send(BlockAllocatedFlag);
        }
        else pipe.recv(BlockAllocatedFlag,node->Thread);

        if (BlockAllocatedFlag==true) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
          //locate the cell
          if (ThisThread==node->Thread) {
            nd=getCenterNodeLocalNumber(i,j,k);
            if (node->block->GetCenterNode(nd)==NULL) nd=-1;

            if (ThisThread!=0) pipe.send(nd);
          }
          else pipe.recv(nd,node->Thread);

          if (ThisThread==0) Signature.add(nd);
        }
      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
        cTreeNodeAMR<cBlockAMR> *downNode=node->downNode[nDownNode];

        if (downNode!=NULL) {
          if (ThisThread==0) {
            Signature.add(downNode->Temp_ID);
            Signature.add(downNode->xmax,3);
            Signature.add(downNode->xmin,3);
          }

          generateMeshSignature(downNode);
        }
      }
    }


    //finish calculation of the mesh signature and send it to all processors
    if (node==rootTree) {
      //add the information of the internal bodies into the signature
      #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
      if (ThisThread==0) {
        list<cInternalBoundaryConditionsDescriptor>::iterator InternalBoundaryDescriptor;

        cInternalSphericalData* Sphere;
        cInternalCircleData* Circle;
        cInternalSphere1DData* Sphere1D;
        cInternalRotationBodyData* RotationBody;

        for (InternalBoundaryDescriptor=InternalBoundaryList.begin();InternalBoundaryDescriptor!=InternalBoundaryList.end();InternalBoundaryDescriptor++) {
          switch(InternalBoundaryDescriptor->BondaryType) {
          case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
            Sphere=(cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement);

            Signature.add(Sphere->Radius);
            Signature.add(Sphere->OriginPosition,3);
            break;
          case _INTERNAL_BOUNDARY_TYPE_CIRCLE_:
            Circle=(cInternalCircleData*)(InternalBoundaryDescriptor->BoundaryElement);

            Signature.add(Circle->Radius);
            Signature.add(Circle->OriginPosition,3);
            break;
          case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:
            Sphere1D=(cInternalSphere1DData*)(InternalBoundaryDescriptor->BoundaryElement);

            Signature.add(Sphere1D->Radius);
            Signature.add(Sphere1D->OriginPosition,3);
            break;
          case _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_:
            RotationBody=(cInternalRotationBodyData*)(InternalBoundaryDescriptor->BoundaryElement);

            Signature.add(RotationBody->OriginPosition,3);
            Signature.add(RotationBody->AxisOfSymmetry,3);
            Signature.add(RotationBody->xAxisMin);
            Signature.add(RotationBody->xAxisMax);
            Signature.add(RotationBody->rSurfaceMax);

            break;
          case _INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_:
            for (int nt=0;nt<CutCell::nBoundaryTriangleFaces;nt++) for (int np=0;np<3;np++) Signature.add(CutCell::BoundaryTriangleFaces[nt].node[np]->x,3);

            break;
          default:
            exit(__LINE__,__FILE__,"Error: The internal boundary type is not recognized");
          }
        }
      }
      #endif


      //close the pipe and set the signature
      if (ThisThread==0) {
        pipe.closeRecvAll();
        MeshSignature=Signature.checksum();
      }
      else pipe.closeSend();

      //de-allocate all buffers associated with the pipe
      pipe.remove();

      //distribute the signature between all processors
      MPI_Bcast(&MeshSignature,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    }
  }

  unsigned long getMeshSignature() {
    if (meshModifiedFlag==true) generateMeshSignature();

    return MeshSignature;
  }

  //set the mesh name
  void setMeshName(const char *str) {
    nMeshModificationCounter++,meshModifiedFlag=true,meshModifiedFlag_CreateNewSpaceFillingCurve=true,meshModifiedFlag_CountMeshElements=true;

    sprintf(MeshName,"%s",str);
//    generateMeshSignature();
  }

  void getMeshName(char *str) {
    sprintf(str,"%s",MeshName);
  }

  void generateMeshName() {
    char mname[STRING_LENGTH];
    time_t TimeValue=time(0);
    tm *ct=localtime(&TimeValue);

    nMeshModificationCounter++,meshModifiedFlag=true,meshModifiedFlag_CreateNewSpaceFillingCurve=true,meshModifiedFlag_CountMeshElements=true;

    sprintf(mname,"noname.%i-%i-%i.AMR.mesh",ct->tm_mon+1,ct->tm_mday,ct->tm_year+1900);
    setMeshName(mname);
  }

  inline long int getCornerNodeLocalNumber(int i,int j,int k) {
    long int nd;

    #if _MESH_DIMENSION_ == 1
    nd=i+_GHOST_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nd=i+_GHOST_CELLS_X_+(1+_TOTAL_BLOCK_CELLS_X_)*(j+_GHOST_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    nd=i+_GHOST_CELLS_X_+(1+_TOTAL_BLOCK_CELLS_X_)*(j+_GHOST_CELLS_Y_+(k+_GHOST_CELLS_Z_)*(1+_TOTAL_BLOCK_CELLS_Y_));
    #endif

    return nd;
  }

  inline static long int getCenterNodeLocalNumber(int i,int j,int k) {
    long int nd;

    #if _MESH_DIMENSION_ == 1
    nd=i+_GHOST_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nd=i+_GHOST_CELLS_X_+_TOTAL_BLOCK_CELLS_X_*(j+_GHOST_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    nd=i+_GHOST_CELLS_X_+_TOTAL_BLOCK_CELLS_X_*(j+_GHOST_CELLS_Y_+(k+_GHOST_CELLS_Z_)*_TOTAL_BLOCK_CELLS_Y_);
    #endif

    return nd;
  }

  inline static void convertCenterNodeLocalNumber2LocalCoordinates(int LocalNumber,int &i,int &j, int &k) {
    k=LocalNumber/(_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_);
    LocalNumber-=k*_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_;

    j=LocalNumber/_TOTAL_BLOCK_CELLS_X_;
    i=LocalNumber-j*_TOTAL_BLOCK_CELLS_X_;

    i-=_GHOST_CELLS_X_;
    j-=_GHOST_CELLS_Y_;
    k-=_GHOST_CELLS_Z_;
  }

  long int findCornerNodeIndex(double *x,int &i,int &j,int &k,cTreeNodeAMR<cBlockAMR>* startNode) { 
    double dx;
    
    //if ((x[0]<startNode->xmin[0])||(startNode->xmax[0]<x[0])) return -1;
    dx=dxRootBlock[0]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
    i=(x[0]-startNode->xmin[0])/dx;

    if (fabs(startNode->xmin[0]+dx*i-x[0])>EPS) {
      i++; 
      if (fabs(startNode->xmin[0]+dx*i-x[0])>EPS) return -1;
    }

    if (_MESH_DIMENSION_>=2) {
      //if ((x[1]<startNode->xmin[1])||(startNode->xmax[1]<x[1])) return -1;
      dx=dxRootBlock[1]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_);
      j=(x[1]-startNode->xmin[1])/dx;

      if (fabs(startNode->xmin[1]+dx*j-x[1])>EPS) {
        j++;
        if (fabs(startNode->xmin[1]+dx*j-x[1])>EPS) return -1;
      }
    }
    else j=0;

    if (_MESH_DIMENSION_==3) {
      //if ((x[2]<startNode->xmin[2])||(startNode->xmax[2]<x[2])) return -1;
      dx=dxRootBlock[2]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_);
      k=(x[2]-startNode->xmin[2])/dx;

      if (fabs(startNode->xmin[2]+dx*k-x[2])>EPS) {
        k++;
        if (fabs(startNode->xmin[2]+dx*k-x[2])>EPS) return -1;
      }
    }
    else k=0;

    return getCornerNodeLocalNumber(i,j,k);
  }  


  long int findCenterNodeIndex(double *x,int &i,int &j,int &k,cTreeNodeAMR<cBlockAMR>* startNode) {
    double dx,dx2;
  
    //if ((x[0]<startNode->xmin[0])||(startNode->xmax[0]<x[0])) return -1;
    dx=dxRootBlock[0]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
    dx2=dx/2.0;
    i=(int)((x[0]-dx2-startNode->xmin[0])/dx);

    if (fabs(startNode->xmin[0]+dx2+dx*i-x[0])>EPS) {
      i++;
      if (fabs(startNode->xmin[0]+dx2+dx*i-x[0])>EPS) return -1;
    }

    if (_MESH_DIMENSION_>=2) {
      //if ((x[1]<startNode->xmin[1])||(startNode->xmax[1]<x[1])) return -1;
      dx=dxRootBlock[1]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_);
      dx2=dx/2.0;
      j=(int)((x[1]-dx2-startNode->xmin[1])/dx);

      if (fabs(startNode->xmin[1]+dx2+dx*j-x[1])>EPS) {
        j++;
        if (fabs(startNode->xmin[1]+dx2+dx*j-x[1])>EPS) return -1;
      }
    }
    else j=0;

    if (_MESH_DIMENSION_==3) {
      //if ((x[2]<startNode->xmin[2])||(startNode->xmax[2]<x[2])) return -1;
      dx=dxRootBlock[2]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_);
      dx2=dx/2.0;
      k=(int)((x[2]-dx2-startNode->xmin[2])/dx);

      if (fabs(startNode->xmin[2]+dx2+dx*k-x[2])>EPS) {
        k++;
        if (fabs(startNode->xmin[2]+dx2+dx*k-x[2])>EPS) return -1;
      }
    }
    else k=0;

    return getCenterNodeLocalNumber(i,j,k);
  }

  //find the index of the cell where the point 'x' is located
  long int fingCellIndex(double *x,int &i,int &j,int &k,cTreeNodeAMR<cBlockAMR>* startNode,bool ExitFlag=true) {
    double dx;

    if ((x[0]<startNode->xmin[0])||(startNode->xmax[0]<x[0])) {
      int idim=0;
      dx=0.0; //add the place for break point for debbuger purposes

      fprintf(DiagnospticMessageStream,"$PREFIX:Error: a point is out pf the box (file=%s, line=%i)!!!!!!\n",__FILE__,__LINE__);
      for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"$PREFIX:idim=%i, x[idim]=%e, xmin[idim]=%e, xmax[idim]=%e\n",idim,x[idim],startNode->xmin[idim],startNode->xmax[idim]);

      if (ExitFlag==true) exit(__LINE__,__FILE__,"x is outside of the block");
      else return -1;
    }

    dx=dxRootBlock[0]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
    i=(int)((x[0]-startNode->xmin[0])/dx);
    if (i==_BLOCK_CELLS_X_) i=_BLOCK_CELLS_X_-1;

    if (_MESH_DIMENSION_>=2) {
      if ((x[1]<startNode->xmin[1])||(startNode->xmax[1]<x[1])) {

        fprintf(DiagnospticMessageStream,"$PREFIX:Error: a point is out pf the box (file=%s, line=%i)!!!!!!\n",__FILE__,__LINE__);
        for (int idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"$PREFIX:idim=%i, x[idim]=%e, xmin[idim]=%e, xmax[idim]=%e\n",idim,x[idim],startNode->xmin[idim],startNode->xmax[idim]);


        if (ExitFlag==true) exit(__LINE__,__FILE__,"x is outside of the block");
        else return -1;
      }

      dx=dxRootBlock[1]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_);
      j=(int)((x[1]-startNode->xmin[1])/dx);
      if (j==_BLOCK_CELLS_Y_) j=_BLOCK_CELLS_Y_-1;
    }
    else j=0;

    if (_MESH_DIMENSION_==3) {
      if ((x[2]<startNode->xmin[2])||(startNode->xmax[2]<x[2])) {
        fprintf(DiagnospticMessageStream,"$PREFIX:Error: a point is out pf the box (file=%s, line=%i)!!!!!!\n",__FILE__,__LINE__);
        for (int idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"$PREFIX:idim=%i, x[idim]=%e, xmin[idim]=%e, xmax[idim]=%e\n",idim,x[idim],startNode->xmin[idim],startNode->xmax[idim]);

        if (ExitFlag==true) exit(__LINE__,__FILE__,"x is outside of the block");
        else return -1;
      }

      dx=dxRootBlock[2]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_);
      k=(int)((x[2]-startNode->xmin[2])/dx);
      if (k==_BLOCK_CELLS_Z_) k=_BLOCK_CELLS_Z_-1;
    }
    else k=0;

    return getCenterNodeLocalNumber(i,j,k);
  }
 

  //the constructor
  void init(double *xMin,double *xMax,double (*localResolutionFunction)(double*)) {
    int i,idim;

    //check the number of bits reserved to store the number of connection's of a corner node
    if ((1<<_MAX_CORNER_NODE_CONNECTION_BITS_)<_MAX_CORNER_NODE_CONNECTION_+1) exit(__LINE__,__FILE__,"The value of _MAX_CORNER_NODE_CONNECTION_BITS_ in not sufficient to store _MAX_CORNER_NODE_CONNECTION_");

    //check if the number of bits reserved to store the blocks refinment level is sufficient to store all levels' numbers
    if ((1<<_MAX_REFINMENT_LEVEL_BITS_)<=_MAX_REFINMENT_LEVEL_) exit(__LINE__,__FILE__,"The value of _MAX_REFINMENT_LEVEL_BITS_ is not sufficient to store all refinment levels");
    if ((_MESH_DIMENSION_<1)||(_MESH_DIMENSION_>3)) exit(__LINE__,__FILE__,"The mesh dimension is wrong");

    xGlobalMin[0]=xMin[0],xGlobalMax[0]=xMax[0],dxRootBlock[0]=(xMax[0]-xMin[0]);
    if (2*_GHOST_CELLS_X_>=_BLOCK_CELLS_X_) exit(__LINE__,__FILE__,"The mesh dimension is wrong");
    EPS=0.0001*dxRootBlock[0]/double(_BLOCK_CELLS_X_)/(1<<_MAX_REFINMENT_LEVEL_); 

    if (_MESH_DIMENSION_>1) {
      xGlobalMin[1]=xMin[1],xGlobalMax[1]=xMax[1],dxRootBlock[1]=(xMax[1]-xMin[1]); 
      if (2*_GHOST_CELLS_Y_>=_BLOCK_CELLS_Y_) exit(__LINE__,__FILE__,"The mesh dimension is wrong");
      if (EPS>0.0001*dxRootBlock[1]/double(_BLOCK_CELLS_Y_)/(1<<_MAX_REFINMENT_LEVEL_)) EPS=0.0001*dxRootBlock[1]/double(_BLOCK_CELLS_Y_)/(1<<_MAX_REFINMENT_LEVEL_); 
    }

    if (_MESH_DIMENSION_>2) {
      xGlobalMin[2]=xMin[2],xGlobalMax[2]=xMax[2],dxRootBlock[2]=(xMax[2]-xMin[2]);
      if (2*_GHOST_CELLS_Z_>=_BLOCK_CELLS_Z_) exit(__LINE__,__FILE__,"The mesh dimension is wrong");
      if (EPS>0.0001*dxRootBlock[2]/double(_BLOCK_CELLS_Z_)/(1<<_MAX_REFINMENT_LEVEL_)) EPS=0.0001*dxRootBlock[2]/double(_BLOCK_CELLS_Z_)/(1<<_MAX_REFINMENT_LEVEL_);
    }

    for (idim=0;idim<_MESH_DIMENSION_;idim++) _MESH_AMR_XMAX_[idim]=xMax[idim],_MESH_AMR_XMIN_[idim]=xMin[idim]; 

    //check the correctness of the symmetry options
#if _MESH_DIMENSION_ == 3
    if (_AMR_SYMMETRY_MODE_!=_AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_) exit(__LINE__,__FILE__,"The symmetry mode and the mesh dimention are not consistent: only AMR_SYMMETRY_MODE = AMR_SYMMETRY_MODE_PLANAR_SYMMETRY is allowed for MESH_DIMENSION == 3");
#else
    exit(__LINE__,__FILE__,"Error: the option is not found");
#endif


    //set the default value for the 'interpolation functions'
    GetCenterNodesInterpolationCoefficients=NULL;
    GetCornerNodesInterpolationCoefficients=NULL;

    localResolution=localResolutionFunction;

    //set up the tree and the root block
//    rootBlock=blocks.newElement();
    rootTree=treeNodes.newElement();

    //init the ID of the root node
    GetAMRnodeID(rootTree->AMRnodeID,rootTree);


    nMeshModificationCounter++,meshNodesNumber=0,meshBlocksNumber=0,meshModifiedFlag=true,meshModifiedFlag_CreateNewSpaceFillingCurve=true,meshModifiedFlag_CountMeshElements=true,meshMaximumRefinmentLevel=0;

    rootTree->upNode=NULL;
    rootTree->block=rootBlock;
    for (i=0;i<(1<<_MESH_DIMENSION_);i++) rootTree->downNode[i]=NULL; 

    rootTree->RefinmentLevel=0;
    for (i=0;i<_MESH_DIMENSION_;i++) rootTree->xmin[i]=xMin[i],rootTree->xmax[i]=xMax[i];


    //add the boundary faces to the root node
    DistributeBoundaryCutBlocks(rootTree,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces);

    accepltTreeNodeFunction=NULL;
    MeshName[0]='\0',MeshSignature=0;
//    blockList=NULL;

    //set the corner nodes of the block
    /*
    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"Adjust as in 3D");
    #elif _MESH_DIMENSION_ == 2
//    exit(__LINE__,__FILE__,"Adjust as in 3D");
    #elif _MESH_DIMENSION_ == 3
//    int iMax=_BLOCK_CELLS_X_,jMax=_BLOCK_CELLS_Y_,kMax=_BLOCK_CELLS_Z_;
    #else
    exit(__LINE__,__FILE__,"The wrong dimension");
    #endif
    */

    AllocateBlock(rootTree);
    rootBlock=rootTree->block;


    /*
    //set up the corner nodes
    cCornerNode *newCornerNode;
    double xNode[3];

    for (k=0;k<=kMax;k++) for (j=0;j<=jMax;j++) for (i=0;i<=iMax;i++) {
      nd=getCornerNodeLocalNumber(i,j,k);

      newCornerNode=CornerNodes.newElement();
      rootBlock->SetCornerNode(newCornerNode,nd);
      newCornerNode->incrementConnectionCounter();

      xNode[0]=xGlobalMin[0]+i*dxRootBlock[0]/_BLOCK_CELLS_X_;
      if (_MESH_DIMENSION_>=2) xNode[1]=xGlobalMin[1]+j*dxRootBlock[1]/_BLOCK_CELLS_Y_;
      if (_MESH_DIMENSION_==3) xNode[2]=xGlobalMin[2]+k*dxRootBlock[2]/_BLOCK_CELLS_Z_;

      newCornerNode->SetX(xNode);
    }


    //set the center node 
    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    cCenterNode *newCenterNode;

    for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
      nd=getCenterNodeLocalNumber(i,j,k);

      newCenterNode=CenterNodes.newElement();
      rootBlock->SetCenterNode(newCenterNode,nd);
      newCenterNode->incrementConnectionCounter();

      xNode[0]=xGlobalMin[0]+(i+0.5)*dxRootBlock[0]/_BLOCK_CELLS_X_;
      if (_MESH_DIMENSION_>=2) xNode[1]=xGlobalMin[1]+(j+0.5)*dxRootBlock[1]/_BLOCK_CELLS_Y_;
      if (_MESH_DIMENSION_==3) xNode[2]=xGlobalMin[2]+(k+0.5)*dxRootBlock[2]/_BLOCK_CELLS_Z_;

      newCenterNode->SetX(xNode);

    }
    #endif
    */
 
    //install the internal surfaces into the rooBlock
    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    list<cInternalBoundaryConditionsDescriptor>::iterator ptr;
    cInternalBoundaryConditionsDescriptor *newDescriptor;

    for (ptr=InternalBoundaryList.begin();ptr!=InternalBoundaryList.end();ptr++) {
      newDescriptor=InternalBoundaryDescriptors.newElement();

      //copy the content of the descriptor
      newDescriptor->BondaryType=ptr->BondaryType;
      newDescriptor->BoundaryElement=ptr->BoundaryElement;

      //add the descriptor to the rootBlock
      newDescriptor->nextInternalBCelement=rootTree->InternalBoundaryDescriptorList;
      rootTree->InternalBoundaryDescriptorList=newDescriptor;
    }

    //add the 'rootTree' to the list of the boundary nodes
    DomainSurfaceBoundaryList=rootTree;

    #endif

    //init the MPI variables
    int MPIinitFlag;

    switch (_AMR_PARALLEL_MODE_) {
    case _AMR_PARALLEL_MODE_ON_:
      MPI_Initialized(&MPIinitFlag);

      if (MPIinitFlag==true) {
        MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
        MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);
      }
      else exit(__LINE__,__FILE__,"Error: MPI is not initialized");

      ParallelNodesDistributionList=new cTreeNodeAMR<cBlockAMR>*[nTotalThreads];

      #if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
      DomainBoundaryLayerNodesList=new cTreeNodeAMR<cBlockAMR>*[nTotalThreads];
      #endif

      ParallelSendRecvMap=new bool*[nTotalThreads];
      ParallelSendRecvMap[0]=new bool [nTotalThreads*nTotalThreads];


      for (int thread=0;thread<nTotalThreads;thread++) {
        ParallelNodesDistributionList[thread]=NULL;
        DomainBoundaryLayerNodesList[thread]=NULL;

        ParallelSendRecvMap[thread]=ParallelSendRecvMap[0]+thread*nTotalThreads;
        for (i=0;i<nTotalThreads;i++) ParallelSendRecvMap[thread][i]=false;
      }
      break;
    case _AMR_PARALLEL_MODE_OFF_:
      ParallelNodesDistributionList=new cTreeNodeAMR<cBlockAMR>*[1];
      ParallelNodesDistributionList[0]=NULL;

      #if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
      DomainBoundaryLayerNodesList=NULL;
      #endif
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }

  }  


  cMeshAMRgeneric() {
    for (int idim=0;idim<_MESH_DIMENSION_;idim++) _MESH_AMR_XMAX_[idim]=0.0,_MESH_AMR_XMIN_[idim]=0.0;

     //set the default value for the 'interpolation functions'
     GetCenterNodesInterpolationCoefficients=NULL;
     GetCornerNodesInterpolationCoefficients=NULL;
     localResolution=NULL;

     //set the defaul value of the diagnostic stream
     DiagnospticMessageStream=stdout;

     //the counter of the load re-balancing operations
     nParallelListRedistributions=0;

     //the counter of any mesh modifications or rebalancing 
     nMeshModificationCounter=0;

     //default value of the parallel mesh generation flag
     ParallelMeshGenerationFlag=false;

     //set up the tree and the root block
     rootBlock=NULL;
     rootTree=NULL;
     meshNodesNumber=0,meshBlocksNumber=0,meshModifiedFlag=true,meshModifiedFlag_CreateNewSpaceFillingCurve=true,meshModifiedFlag_CountMeshElements=true,meshMaximumRefinmentLevel=0;

     AllowBlockAllocation=true;
     DeallocateUnusedBlocks=true;

     EPS=1.0E-20;

     ThisThread=0,nTotalThreads=1;
     ParallelSendRecvMap=NULL;

     accepltTreeNodeFunction=NULL;
     MeshName[0]='\0',MeshSignature=0;
//     blockList=NULL;
     ParallelNodesDistributionList=NULL;

     #if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
     DomainBoundaryLayerNodesList=NULL;
     #endif

     #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
     DomainSurfaceBoundaryList=NULL;
     #endif

     #if _AMR__CUT_CELL__MODE_ ==  _AMR__CUT_CELL__MODE__ON_
     CutCellSurfaceLocalResolution=NULL;
     #endif

     nMpiBarrierCalls=0;
  }

  //register the 'internal boundary' (the surface determining cut cells)
  void RegisterInternalBoundary(cInternalBoundaryConditionsDescriptor Descriptor) {
    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    if (rootTree!=NULL) exit(__LINE__,__FILE__,"Error: all internal surface must be registered before initialization of the mesh");

    #if _MESH_DIMENSION_ == 1
    if (Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_1D_SPHERE_) exit(__LINE__,__FILE__,"Error: Attempted boundary type is not allowed for _MESH_DIMENSION_ == 1");
    #elif _MESH_DIMENSION_ == 2
    if (Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_CIRCLE_) exit(__LINE__,__FILE__,"Error: Attempted boundary type is not allowed for _MESH_DIMENSION_ == 2");
    #elif _MESH_DIMENSION_ == 3
    if ((Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_)&&(Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_)&&(Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_)) {
      exit(__LINE__,__FILE__,"Error: Attempted boundary type is not allowed for _MESH_DIMENSION_ == 3");
    }
    #else
    exit(__LINE__,__FILE__,"Error: unknown value of _MESH_DIMENSION_");
    #endif


    InternalBoundaryList.push_back(Descriptor);
    #else
    exit(__LINE__,__FILE__,"Error: internal boundary is allowed only when _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_");
    #endif
  }



/*
  inline cTreeNodeAMR<cBlockAMR>* findTreeNode(double *x,cTreeNodeAMR<cBlockAMR>  *startNode=NULL) {
    cTreeNodeAMR<cBlockAMR> *res=NULL,*t=NULL;
    double blockBasedCoordinates[3]={0.0,0.0,0.0};
    int i,j=0,k=0;

    if (startNode==NULL) startNode=rootTree;

Start:
    blockBasedCoordinates[0]=(x[0]-startNode->xmin[0])/dxRootBlock[0]*(1<<startNode->RefinmentLevel);

#if _MESH_DIMENSION_ > 1
    blockBasedCoordinates[1]=(x[1]-startNode->xmin[1])/dxRootBlock[1]*(1<<startNode->RefinmentLevel);
#endif

#if _MESH_DIMENSION_ > 2
    blockBasedCoordinates[2]=(x[2]-startNode->xmin[2])/dxRootBlock[2]*(1<<startNode->RefinmentLevel);
#endif

    if ((0.0<=blockBasedCoordinates[0])&&(blockBasedCoordinates[0]<=1.0)&&(0.0<=blockBasedCoordinates[1])&&(blockBasedCoordinates[1]<=1.0)&&(0.0<=blockBasedCoordinates[2])&&(blockBasedCoordinates[2]<=1.0)) {
      i=(blockBasedCoordinates[0]<0.5) ? 0 : 1;

#if _MESH_DIMENSION_ > 1
      j=(blockBasedCoordinates[1]<0.5) ? 0 : 1;
#endif

#if _MESH_DIMENSION_ > 2
      k=(blockBasedCoordinates[2]<0.5) ? 0 : 1;
#endif

      t=startNode->downNode[i+2*(j+2*k)];

      if (t!=NULL) {
        startNode=t;
        goto Start;
      }
      else res=startNode;
    }
    else {
      if (blockBasedCoordinates[0]<0.0) i=-1;
      else if (blockBasedCoordinates[0]>1.0) i=1;
      else i=0;

      if (blockBasedCoordinates[1]<0.0) j=-1;
      else if (blockBasedCoordinates[1]>1.0) j=1;
      else j=0;

      if (blockBasedCoordinates[2]<0.0) k=-1;
      else if (blockBasedCoordinates[2]>1.0) k=1;
      else k=0;

      t=startNode->GetNeibNode(i,j,k);


      if (t!=NULL) {
        startNode=t;
        goto Start;
      }
      else {
        if ((t=startNode->upNode)!=NULL) {
          startNode=t;
          goto Start;
        }
        else res=NULL;
      }

    }

#if  _AMR_DEBUG_MODE_ == _AMR_DEBUG_MODE_ON_
    if (res!=NULL) for (register int idim=0;idim<_MESH_DIMENSION_;idim++) if ((x[idim]<res->xmin[idim])||(res->xmax[idim]<x[idim])) exit(__LINE__,__FILE__,"Error: did'nt find the tree node");
#endif
    return res;
  }
  */

  inline cTreeNodeAMR<cBlockAMR>* findTreeNode(double *x,cTreeNodeAMR<cBlockAMR>  *startNode=NULL) {
    cTreeNodeAMR<cBlockAMR> *res=NULL,*t=NULL;
    int iState=0,jState=0,kState=0,i,j,k;
    double xmin[3],xmax[3];


    if (startNode==NULL) startNode=rootTree;

Start:
    memcpy(xmin,startNode->xmin,_MESH_DIMENSION_*sizeof(double));
    memcpy(xmax,startNode->xmax,_MESH_DIMENSION_*sizeof(double));

    if (x[0]<xmin[0]) iState=-1;
    else if (x[0]>=xmax[0]) iState=+1;
    else iState=0;

#if _MESH_DIMENSION_ > 1
    if (x[1]<xmin[1]) jState=-1;
    else if (x[1]>=xmax[1]) jState=+1;
    else jState=0;
#endif

#if _MESH_DIMENSION_ > 2
    if (x[2]<xmin[2]) kState=-1;
    else if (x[2]>=xmax[2]) kState=+1;
    else kState=0;
#endif

    if ((iState==0)&&(jState==0)&&(kState==0)) {
      i=((x[0]-xmin[0])/dxRootBlock[0]*(1<<startNode->RefinmentLevel)<0.5) ? 0 : 1;

#if _MESH_DIMENSION_ > 1
      j=((x[1]-xmin[1])/dxRootBlock[1]*(1<<startNode->RefinmentLevel)<0.5) ? 0 : 1;
#endif

#if _MESH_DIMENSION_ > 2
      k=((x[2]-xmin[2])/dxRootBlock[2]*(1<<startNode->RefinmentLevel)<0.5) ? 0 : 1;
#endif

      t=startNode->downNode[i+2*(j+2*k)]; 

      if (t!=NULL) {
        startNode=t;
        goto Start;
      }
      else res=startNode;
    }
    else {
      t=startNode->GetNeibNode(iState,jState,kState);

      if (t!=NULL) {
        startNode=t;
        goto Start;
      }
      else {
        if ((t=startNode->upNode)!=NULL) {
          startNode=t;
          goto Start;
        }
        else res=NULL;
      }

    }

#if  _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if (res!=NULL) for (register int idim=0;idim<_MESH_DIMENSION_;idim++) if ((x[idim]<res->xmin[idim])||(res->xmax[idim]<x[idim])) exit(__LINE__,__FILE__,"Error: did'nt find the tree node");
#endif
    return res;
  }

  //find the node of the tree where the point 'x' is located BUT the nodes' resolution should be less of equal to UpperResolutionLevel (the found free could be not the last in the tree (can be located somewhere in the middle of the tree)
  cTreeNodeAMR<cBlockAMR>*  findTreeNodeLimitedResolutionLevel(double *x,int UpperResolutionLevel,cTreeNodeAMR<cBlockAMR>  *startNode=NULL) {
    cTreeNodeAMR<cBlockAMR> *res=NULL,*t=NULL;
    double blockBasedCoordinates[3]={0.0,0.0,0.0};
    int i,j=0,k=0;

    if (startNode==NULL) startNode=rootTree;

    //go up on the tree
    if (startNode->RefinmentLevel>UpperResolutionLevel) while (startNode->RefinmentLevel==UpperResolutionLevel) startNode=startNode->upNode;


    blockBasedCoordinates[0]=(x[0]-startNode->xmin[0])/dxRootBlock[0]*(1<<startNode->RefinmentLevel);
    if (_MESH_DIMENSION_>1) blockBasedCoordinates[1]=(x[1]-startNode->xmin[1])/dxRootBlock[1]*(1<<startNode->RefinmentLevel);
    if (_MESH_DIMENSION_>2) blockBasedCoordinates[2]=(x[2]-startNode->xmin[2])/dxRootBlock[2]*(1<<startNode->RefinmentLevel);

    if ((0.0<=blockBasedCoordinates[0])&&(blockBasedCoordinates[0]<=1.0)&&(0.0<=blockBasedCoordinates[1])&&(blockBasedCoordinates[1]<=1.0)&&(0.0<=blockBasedCoordinates[2])&&(blockBasedCoordinates[2]<=1.0)) {
      i=(blockBasedCoordinates[0]<0.5) ? 0 : 1;
      if (_MESH_DIMENSION_>1) j=(blockBasedCoordinates[1]<0.5) ? 0 : 1;
      if (_MESH_DIMENSION_>2) k=(blockBasedCoordinates[2]<0.5) ? 0 : 1;

      t=(startNode->RefinmentLevel<UpperResolutionLevel) ? startNode->downNode[i+2*(j+2*k)] : NULL;

      res=(t!=NULL) ? findTreeNodeLimitedResolutionLevel(x,UpperResolutionLevel,t) : startNode;
    }
    else res=(startNode->upNode!=NULL) ? findTreeNodeLimitedResolutionLevel(x,UpperResolutionLevel,startNode->upNode) : NULL;

    return res;
  }


/*
  void reconnectDownTreeNode(cTreeNodeAMR<cBlockAMR> *startNode) {
    cTreeNodeAMR<cBlockAMR> *neibNode,*neibDownNode; 

    #if _MESH_DIMENSION_ == 1
    const static int neibNodeNumber[2][2]={ {0,-1}, {-1,1}};
    const static int neibDownNodeNumber[2][2]={ {1,-1}, {-1,0}};
    const static int reverseFaceMap[2]={1,0};


    #elif _MESH_DIMENSION_ == 2
 
    //the face number of the parent block that the downBlock is connecter through with the surrounding blocks
    const static int neibNodeNumber[4][4]={ {0,-1,2,-1}, {-1,1,2,-1}, {0,-1,-1,3}, {-1,1,-1,3}};

    //the number of the downNodes in the neibNode that the startNode->downNode[?] is connected with
    const static int neibDownNodeNumber[4][4]={ {1,-1,2,-1}, {-1,0,3,-1}, {3,-1,-1,0}, {-1,2,-1,1}};

    //the corresponding face number that neibNode is connected to startNode
    const static int reverseFaceMap[4]={2,0,3,1};

    #else

    //the face number of the parent block that the downBlock is connecter through with the surrounding blocks
    const static int neibNodeNumber[8][6]={ {0,-1,2,-1,4,-1}, {-1,1,2,-1,4,-1}, {0,-1,-1,3,4,-1}, {-1,1,-1,3,4,-1},
                                            {0,-1,2,-1,-1,5}, {-1,1,2,-1,-1,5}, {0,-1,-1,3,-1,5}, {-1,1,-1,3,-1,5},   };

    const static int neibDownNodeNumber[8][6]={ {1,-1,2,-1,4,-1,}, {-1,0,3,-1,5,-1}, {3,-1,-1,0,6,-1}, {-1,2,-1,1,7,-1},
                                                {5,-1,6,-1,-1,0}, {-1,4,7,-1,-1,1,}, {7,-1,-1,4,-1,2}, {-1,6,-1,5,-1,3}  }; 

    const static int reverseFaceMap[6]={1,0,3,2,5,4};



    #endif

    int nDownNode,nface;



    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if (neibNodeNumber[nDownNode][nface]!=-1) {
      neibNode=startNode->neibNode[neibNodeNumber[nDownNode][nface]];  

      if (neibNode!=NULL) {
        neibDownNode=neibNode->downNode[neibDownNodeNumber[nDownNode][nface]];

        startNode->downNode[nDownNode]->neibNode[nface]=neibDownNode;
        if (neibDownNode!=NULL) neibDownNode->neibNode[ reverseFaceMap[nface] ]=startNode->downNode[nDownNode]; 
      }
    }


  }
*/


/*
  cTreeNodeAMR<cBlockAMR> *getNeibNode_2D(int i,int j,cTreeNodeAMR<cBlockAMR>* startNode) {
    double x[3]={0.0,0.0,0.0};

    x[0]=startNode->xmin[0]+0.5*dxRootBlock[0]/(1<<startNode->RefinmentLevel)*(1.0+1.01*i);  
    x[1]=startNode->xmin[1]+0.5*dxRootBlock[1]/(1<<startNode->RefinmentLevel)*(1.0+1.01*j);

    return findTreeNode(x,startNode); 
  }
  */

  cTreeNodeAMR<cBlockAMR> *getNeibNode(int i,int j,int k,cTreeNodeAMR<cBlockAMR>* startNode) {
    return startNode->GetNeibNode(i,j,k);
  }

  cTreeNodeAMR<cBlockAMR> *getNeibNode_DirectTreeSearch(int i,int j,int k,cTreeNodeAMR<cBlockAMR>* startNode) {
    double x[3]={0.0,0.0,0.0};
    cTreeNodeAMR<cBlockAMR> *res;

    x[0]=startNode->xmin[0]+0.5*dxRootBlock[0]/(1<<startNode->RefinmentLevel)*(1.0+1.01*i);
    if (_MESH_DIMENSION_>1) x[1]=startNode->xmin[1]+0.5*dxRootBlock[1]/(1<<startNode->RefinmentLevel)*(1.0+1.01*j);
    if (_MESH_DIMENSION_>2) x[2]=startNode->xmin[2]+0.5*dxRootBlock[2]/(1<<startNode->RefinmentLevel)*(1.0+1.01*k);

    res=findTreeNode(x,startNode);
    
    return res;
  }
 

  /*
 void collectNeibNodeMap_1D(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
   int i,j,nface,iDownBlock,iNeib,iOffset;
   cTreeNodeAMR<cBlockAMR> *neibNode,*downNode;


exit(__LINE__,__FILE__,"update as in 3D collect"); 


   //if the neibNode->RefinmentLevel == startNodes->RefinementLevel -> the 'neibNode' does not contribute to the 'newCornerNodeMap'

   //check if the face neibours can contribute to the lists
   for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=startNode->neibNode[nface])!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
     //transfer the nodes

     static const int _X_FACE_NEIB_OFFSET_[2]={-2*_BLOCK_CELLS_X_,2*_BLOCK_CELLS_X_};

     for (iDownBlock=0;iDownBlock<2;iDownBlock++) { 
       downNode=neibNode->downNode[iDownBlock];
       if (downNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) exit(__LINE__,__FILE__,"the branch is too deap");
       iOffset=_X_FACE_NEIB_OFFSET_[nface]+_BLOCK_CELLS_X_*iDownBlock;

       for (iNeib=-_GHOST_CELLS_X_;iNeib<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;iNeib++) if ((-_GHOST_CELLS_X_<=iNeib+iOffset)&&(iNeib+iOffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
         newCornerNodeMap[iNeib+iOffset+_GHOST_CELLS_X_][0][0]=downNode->block->cornerNodes[getCornerNodeLocalNumber(iNeib,0,0)];
       }

     }
   }
 }


  void collectNeibNodeMap_2D(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
    int i,j,nface,iDownBlock,jDownBlock,iNeib,jNeib,iOffset,jOffset;
    cTreeNodeAMR<cBlockAMR> *neibNode,*downNode;

exit(__LINE__,__FILE__,"update as in 3D collect");

    //if the neibNode->RefinmentLevel == startNodes->RefinementLevel -> the 'neibNode' does not contribute to the 'newCornerNodeMap' 

    //check if the face neibours can contribute to the lists
    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=startNode->neibNode[nface])!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    //transfer the nodes

    static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
    static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;

    
    static const int _X_FACE_NEIB_OFFSET_[4]={-2*_BLOCK_CELLS_X_,0,2*_BLOCK_CELLS_X_,0};
    static const int _Y_FACE_NEIB_OFFSET_[4]={0,-2*_BLOCK_CELLS_Y_,0,2*_BLOCK_CELLS_Y_}; 
 
    for (jDownBlock=0;jDownBlock<2;jDownBlock++) for (iDownBlock=0;iDownBlock<2;iDownBlock++) {
      downNode=neibNode->downNode[iDownBlock+2*jDownBlock]; 

      iOffset=_X_FACE_NEIB_OFFSET_[nface]+_BLOCK_CELLS_X_*iDownBlock;
      jOffset=_Y_FACE_NEIB_OFFSET_[nface]+_BLOCK_CELLS_Y_*jDownBlock; 

      //chech if 'downNode' can overlap 'startNode'
      int nIntersection=0;

      if (iOffset>0) {
        if (iOffset<=2*_BLOCK_CELLS_X_) nIntersection++;
      }
      else if (iOffset<0) {
        if (iOffset==-_BLOCK_CELLS_X_) nIntersection++;
      }
      else nIntersection++;

      if (jOffset>0) {
        if (jOffset<=2*_BLOCK_CELLS_Y_) nIntersection++;
      }
      else if (jOffset<0) {
        if (jOffset==-_BLOCK_CELLS_Y_) nIntersection++;
      }
      else nIntersection++;

      if (nIntersection!=2) continue;

      if (downNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) exit(__LINE__,__FILE__,"the branch is too deap");

      for (iNeib=-_GHOST_CELLS_X_;iNeib<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;iNeib++) if ((-_GHOST_CELLS_X_<=iNeib+iOffset)&&(iNeib+iOffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {   
        for (jNeib=-_GHOST_CELLS_Y_;jNeib<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jNeib++) if ((-_GHOST_CELLS_Y_<=jNeib+jOffset)&&(jNeib+jOffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {  
          newCornerNodeMap[iNeib+iOffset+_GHOST_CELLS_X_][jNeib+jOffset+_GHOST_CELLS_Y_][0]=downNode->block->cornerNodes[getCornerNodeLocalNumber(iNeib,jNeib,0)]; 
        }
      }
    }
  }

  //check the corner blocks
  int iNeibBlock,jNeibBlock;

  for (iNeibBlock=-1;iNeibBlock<2;iNeibBlock+=2) for (jNeibBlock=-1;jNeibBlock<2;jNeibBlock+=2) if ((neibNode=getNeibNode_2D(iNeibBlock,jNeibBlock,startNode))!=NULL) if (neibNode->RefinmentLevel>startNode->RefinmentLevel) { 
    //transfer the nodes
    if ((neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_)||(neibNode->RefinmentLevel!=startNode->RefinmentLevel+1)) exit(__LINE__,__FILE__,"the branch is too deap");
     
    iOffset=_BLOCK_CELLS_X_*((iNeibBlock==-1) ? -1 : 2);
    jOffset=_BLOCK_CELLS_Y_*((jNeibBlock==-1) ? -1 : 2); 

    for (iNeib=-_GHOST_CELLS_X_;iNeib<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;iNeib++) if ((-_GHOST_CELLS_X_<=iNeib+iOffset)&&(iNeib+iOffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
      for (jNeib=-_GHOST_CELLS_Y_;jNeib<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jNeib++) if ((-_GHOST_CELLS_Y_<=jNeib+jOffset)&&(jNeib+jOffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
        newCornerNodeMap[iNeib+iOffset+_GHOST_CELLS_X_][jNeib+jOffset+_GHOST_CELLS_Y_][0]=neibNode->block->cornerNodes[getCornerNodeLocalNumber(iNeib,jNeib,0)];
      }
    }

  }

}
 





  void collectNeibCornerNodeMap_3D(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
    int nIntersection,ii,jj,kk,iNeibNode,jNeibNode,kNeibNode,iDownNode,jDownNode,kDownNode,ioffset,joffset,koffset;
    cTreeNodeAMR<cBlockAMR> *neibNode,*downNode,*upNode;

    //init the node's map
    for (ii=0;ii<1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_);ii++) for (jj=0;jj<1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_);jj++) for (kk=0;kk<1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_);kk++) newCornerNodeMap[ii][jj][kk]=NULL;

    //transfer nodes to the 'newCornerNodeMap' from the 'cCornerNode'
    for (ii=-(_GHOST_CELLS_X_/2);ii<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_/2;ii++) for (jj=-(_GHOST_CELLS_Y_/2);jj<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_/2;jj++) for (kk=-(_GHOST_CELLS_Z_/2);kk<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_/2;kk++) {
      newCornerNodeMap[2*ii+_GHOST_CELLS_X_][2*jj+_GHOST_CELLS_Y_][2*kk+_GHOST_CELLS_Z_]=startNode->block->GetCornerNode(getCornerNodeLocalNumber(ii,jj,kk));
    }

    //transfer nodes from the neibouring blocks
    for (iNeibNode=-1;iNeibNode<=1;iNeibNode++) for (jNeibNode=-1;jNeibNode<=1;jNeibNode++) for (kNeibNode=-1;kNeibNode<=1;kNeibNode++) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
        upNode=neibNode->upNode; 

        for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++)  for (kDownNode=0;kDownNode<2;kDownNode++) {
          downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

          ioffset=_BLOCK_CELLS_X_*(iDownNode+2*iNeibNode);
          joffset=_BLOCK_CELLS_Y_*(jDownNode+2*jNeibNode);
          koffset=_BLOCK_CELLS_Z_*(kDownNode+2*kNeibNode);

          nIntersection=0;

          //chech if 'downNode' can overlap 'startNode'
          if ((ioffset==0)||(ioffset==2*_BLOCK_CELLS_X_)||(ioffset==_BLOCK_CELLS_X_)||(ioffset==-_BLOCK_CELLS_X_)) nIntersection++; 
          if ((joffset==0)||(joffset==2*_BLOCK_CELLS_Y_)||(joffset==_BLOCK_CELLS_Y_)||(joffset==-_BLOCK_CELLS_Y_)) nIntersection++;
          if ((koffset==0)||(koffset==2*_BLOCK_CELLS_Z_)||(koffset==_BLOCK_CELLS_Z_)||(koffset==-_BLOCK_CELLS_Z_)) nIntersection++;

          if (nIntersection!=3) continue;

          for (ii=-_GHOST_CELLS_X_;ii<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
            for (jj=-_GHOST_CELLS_Y_;jj<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
              for (kk=-_GHOST_CELLS_Z_;kk<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=2*_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
                newCornerNodeMap[ii+ioffset+_GHOST_CELLS_X_][jj+joffset+_GHOST_CELLS_Y_][kk+koffset+_GHOST_CELLS_Z_]=downNode->block->GetCornerNode(getCornerNodeLocalNumber(ii,jj,kk));
              }
            }
          }

        }
      }
    }
 
  }




	
  void collectNeibCenterNodeMap_3D(cTreeNodeAMR<cBlockAMR> *startNode,cCenterNode *newCenterNodeMap[2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
    int nIntersection,ii,jj,kk,iNeibNode,jNeibNode,kNeibNode,iDownNode,jDownNode,kDownNode,ioffset,joffset,koffset;
    cTreeNodeAMR<cBlockAMR> *neibNode,*downNode,*upNode;

    //init the node's map
    for (ii=0;ii<2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_);ii++) for (jj=0;jj<2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_);jj++) for (kk=0;kk<2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_);kk++) newCenterNodeMap[ii][jj][kk]=NULL;

    //transfer nodes from the neibouring blocks
    for (iNeibNode=-1;iNeibNode<=1;iNeibNode++) for (jNeibNode=-1;jNeibNode<=1;jNeibNode++) for (kNeibNode=-1;kNeibNode<=1;kNeibNode++) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
        upNode=neibNode->upNode;

       //transfer nodes from the neibouring blocks
       for (iNeibNode=-1;iNeibNode<=1;iNeibNode++) for (jNeibNode=-1;jNeibNode<=1;jNeibNode++) for (kNeibNode=-1;kNeibNode<=1;kNeibNode++) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
        if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
          upNode=neibNode->upNode;

          for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++)  for (kDownNode=0;kDownNode<2;kDownNode++) {
            downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

            ioffset=_BLOCK_CELLS_X_*(iDownNode+2*iNeibNode);
            joffset=_BLOCK_CELLS_Y_*(jDownNode+2*jNeibNode);
            koffset=_BLOCK_CELLS_Z_*(kDownNode+2*kNeibNode);

            nIntersection=0;

            //chech if 'downNode' can overlap 'startNode'
            if ((ioffset==0)||(ioffset==2*_BLOCK_CELLS_X_)||(ioffset==_BLOCK_CELLS_X_)||(ioffset==-_BLOCK_CELLS_X_)) nIntersection++;
            if ((joffset==0)||(joffset==2*_BLOCK_CELLS_Y_)||(joffset==_BLOCK_CELLS_Y_)||(joffset==-_BLOCK_CELLS_Y_)) nIntersection++;
            if ((koffset==0)||(koffset==2*_BLOCK_CELLS_Z_)||(koffset==_BLOCK_CELLS_Z_)||(koffset==-_BLOCK_CELLS_Z_)) nIntersection++;

            if (nIntersection!=3) continue;

            for (ii=-_GHOST_CELLS_X_;ii<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1)) {
              for (jj=-_GHOST_CELLS_Y_;jj<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1)) {
                for (kk=-_GHOST_CELLS_Z_;kk<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_-1;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=2*_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_-1)) {
                  newCenterNodeMap[ii+ioffset+_GHOST_CELLS_X_][jj+joffset+_GHOST_CELLS_Y_][kk+koffset+_GHOST_CELLS_Z_]=downNode->block->GetCenterNode(getCenterNodeLocalNumber(ii,jj,kk));
                }
              }
            }

          }
        }
      }

    }
  }
}

*/


void GetMeshTreeStatistics(cTreeNodeAMR<cBlockAMR> *startNode=NULL) {
  static long int nAllocatedBlocks=0,nAllocatedBoundaryLayerBlocks=0,*nBlocksPerProcessor=NULL,nAllocatedBlocksUpperTreeBranches=0;
  int thread;

  if (startNode==NULL) {
    startNode=rootTree;
    nBlocksPerProcessor=new long int [nTotalThreads];

    for (thread=0;thread<nTotalThreads;thread++) nAllocatedBlocks=0,nBlocksPerProcessor[thread]=0;
    nAllocatedBlocks=0,nAllocatedBlocksUpperTreeBranches=0,nAllocatedBoundaryLayerBlocks=0;
  }

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {

#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
    if ((startNode->Thread<0)||(startNode->Thread>=nTotalThreads)) exit(__LINE__,__FILE__,"Error: the thread number is out of range");

    if (startNode->block!=NULL) {
      nAllocatedBlocks++;
      if (startNode->Thread!=ThisThread) nAllocatedBoundaryLayerBlocks++;
    }


    nBlocksPerProcessor[startNode->Thread]++;
#else
    if (startNode->block!=NULL) nAllocatedBlocks++;
    nBlocksPerProcessor[0]++;
#endif

  }
  else {
    if (startNode->block!=NULL) nAllocatedBlocksUpperTreeBranches++;

    for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) GetMeshTreeStatistics(startNode->downNode[nDownNode]);
  }

  if (startNode==rootTree) {
#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
    long int buffer[nTotalThreads];
    int i;
    MPI_Status status;

    //compare the blocks distribution over the processors
    if (ThisThread==0) {
      for (thread=1;thread<nTotalThreads;thread++) {
        MPI_Recv(buffer,nTotalThreads,MPI_LONG,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        for (i=0;i<nTotalThreads;i++) if (buffer[i]!=nBlocksPerProcessor[i]) exit(__LINE__,__FILE__,"Error: The blocks distribution on different processors is not the same");
       }
    }
    else MPI_Send(nBlocksPerProcessor,nTotalThreads,MPI_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);



    //output the collected information
    //if (ThisThread==0) *DiagnospticMessageStream << "$PREFIX:Tree distribution statistical data:\nThread\t|nTotal Allocated Blocks\t|nAllocated Blocks Upper Tree Branches\t|nAllocated Blocks Domain Boundary Layer\t|nTotal Nodes Per Thread\n";
    if (ThisThread==0) fprintf(DiagnospticMessageStream,"$PREFIX:Tree distribution statistical data:\n$PREFIX:Thread\t|nTotal Allocated Blocks\t|nAllocated Blocks Upper Tree Branches\t|nAllocated Blocks Domain Boundary Layer\t|nTotal Nodes Per Thread\n");


    for (thread=0;thread<nTotalThreads;thread++) {
      if (ThisThread!=0) {
        if (ThisThread==thread) {
          MPI_Send(&nAllocatedBlocks,1,MPI_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
          MPI_Send(&nAllocatedBoundaryLayerBlocks,1,MPI_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
          MPI_Send(&nAllocatedBlocksUpperTreeBranches,1,MPI_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
        }
      }
      else {
        if (thread!=0) {
          MPI_Recv(&nAllocatedBlocks,1,MPI_LONG,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
          MPI_Recv(&nAllocatedBoundaryLayerBlocks,1,MPI_LONG,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
          MPI_Recv(&nAllocatedBlocksUpperTreeBranches,1,MPI_LONG,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
        }

//        *DiagnospticMessageStream << "$PREFIX:" << thread << "\t" << nAllocatedBlocks << "\t\t\t\t\t" << nAllocatedBlocksUpperTreeBranches << "\t\t\t\t\t" << nAllocatedBoundaryLayerBlocks << "\t\t\t\t\t" << nBlocksPerProcessor[thread] <<  std::endl;
        fprintf(DiagnospticMessageStream,"$PREFIX: %i\t%ld\t\t\t\t\t%ld\t\t\t\t\t%ld\t\t\t\t\t%ld\n",thread,nAllocatedBlocks,nAllocatedBlocksUpperTreeBranches,nAllocatedBoundaryLayerBlocks,nBlocksPerProcessor[thread]);


      }
    }

    delete [] nBlocksPerProcessor;

    nBlocksPerProcessor=NULL;
#else
    int NodesBoundaryLayer;

    *DiagnospticMessageStream << "$PREFIX:Tree distribution statistical data:\n$PREFIX:Thread\t\n$PREFIX:Total Allocated Blocks\tnAllocated Blocks Upper Tree Branches\tnTotal Nodes Per Thread\tnNodes in the Boundary Layer\n";
    for (i=0;i<nTotalThreads;i++) if (i!=thread) NodesBoundaryLayer+=nBlocksPerProcessor[i];
    *DiagnospticMessageStream << "$PREFIX:" << thread << "\t" << nAllocatedBlocks << "\t" << nAllocatedBlocksUpperTreeBranches << "\t" << nBlocksPerProcessor[thread] << "\t" << NodesBoundaryLayer << std::endl;
#endif

  }

}
	

void checkMeshConsistency(cTreeNodeAMR<cBlockAMR> *startNode) {
  int i,j,k,iDownBlock,jDownBlock,kDownBlock,iStartNode,iNeibNode;
  cTreeNodeAMR<cBlockAMR> *neibNode,*downNode;
  bool flag;

  #if _MESH_DIMENSION_ == 1
  static const int iDownNodeMax=1,jDownNodeMax=0,kDownNodeMax=0;
  static const int iMin=-1,iMax=2,jMin=0,jMax=1,kMin=0,kMax=1;
  static const int nCornerNodes=1+_TOTAL_BLOCK_CELLS_X_;
  #elif _MESH_DIMENSION_ == 2
  static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=0;
  static const int iMin=-1,iMax=2,jMin=-1,jMax=2,kMin=0,kMax=1;
  static const int nCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);  
  #else 
  static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=1;
  static const int iMin=-1,iMax=2,jMin=-1,jMax=2,kMin=-1,kMax=2;
  static const int nCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
  #endif



  if (startNode==rootTree) {
    //reset the 'processed flag' for the nodes
    resetNodeProcessedFlag();
  }

  //check the 'xmax' and 'xmin' value of the node
  if (startNode!=rootTree) {
    double *xminUpNode,*xmaxUpNode,*xmin,*xmax;
//    int kDownBlock,jDownBlock,iDownBlock;
//    bool flag;

    flag=false;
    xminUpNode=startNode->upNode->xmin,xmaxUpNode=startNode->upNode->xmax;
    xmin=startNode->xmin,xmax=startNode->xmax;

    for (kDownBlock=0;kDownBlock<=kDownNodeMax;kDownBlock++) for (jDownBlock=0;jDownBlock<=jDownNodeMax;jDownBlock++) for (iDownBlock=0;iDownBlock<=iDownNodeMax;iDownBlock++) {
      if (startNode->upNode->downNode[iDownBlock+2*(jDownBlock+2*kDownBlock)]==startNode) {
        flag=true; //there is a connection between the upNode and the startNode;

        if (fabs(xminUpNode[0]+(xmaxUpNode[0]-xminUpNode[0])*iDownBlock/2.0-xmin[0])>EPS) exit(__LINE__,__FILE__,"Block's xmin value is not correct");
        if (_MESH_DIMENSION_>=2) if (fabs(xminUpNode[1]+(xmaxUpNode[1]-xminUpNode[1])*jDownBlock/2.0-xmin[1])>EPS) exit(__LINE__,__FILE__,"Block's xmin value is not correct");
        if (_MESH_DIMENSION_==3) if (fabs(xminUpNode[2]+(xmaxUpNode[2]-xminUpNode[2])*kDownBlock/2.0-xmin[2])>EPS) exit(__LINE__,__FILE__,"Block's xmin value is not correct");

        if (fabs(xminUpNode[0]+(xmaxUpNode[0]-xminUpNode[0])*(iDownBlock+1)/2.0-xmax[0])>EPS) exit(__LINE__,__FILE__,"Block's xmax value is not correct");
        if (_MESH_DIMENSION_>=2) if (fabs(xminUpNode[1]+(xmaxUpNode[1]-xminUpNode[1])*(jDownBlock+1)/2.0-xmax[1])>EPS) exit(__LINE__,__FILE__,"Block's xmax value is not correct");
        if (_MESH_DIMENSION_==3) if (fabs(xminUpNode[2]+(xmaxUpNode[2]-xminUpNode[2])*(kDownBlock+1)/2.0-xmax[2])>EPS) exit(__LINE__,__FILE__,"Block's xmax value is not correct");
      }
    }

    if (flag==false) exit(__LINE__,__FILE__,"There is a disconnection between the upNode and the startNode");
  }


  //check the tree node 
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //&&(AllowBlockAllocation==true)) {

//    if (startNode->block==NULL) return;

    //check consistency of the nodes': connection though the corner nodes'
    for (int nd=0;nd<(1<<_MESH_DIMENSION_);nd++) {
      int ndNeib;

      neibNode=startNode->GetNeibCorner(nd);

      if (neibNode!=NULL) {
        if ((neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_)||(abs(startNode->RefinmentLevel-neibNode->RefinmentLevel)>1)) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");

        for (ndNeib=0;ndNeib<(1<<_MESH_DIMENSION_);ndNeib++) if (neibNode->GetNeibCorner(ndNeib)==startNode) {
          break;
        }

        if (ndNeib==(1<<_MESH_DIMENSION_)) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");
      }
    }

    //check consistency of the nodes': connection though the faces
#if _MESH_DIMENSION_ == 1
    for (int nface=0;nface<2;nface++) {
      int nfaceNeib;
      bool found;

      neibNode=startNode->GetNeibFace(nface,0,0);

      if (neibNode!=NULL) {
        found=false;
        if ((neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_)||(abs(startNode->RefinmentLevel-neibNode->RefinmentLevel)>1)) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");

        for (nfaceNeib=0;nfaceNeib<2;nfaceNeib++) if (found==false) if (neibNode->GetNeibFace(nfaceNeib,0,0)==startNode) {
          found=true;
          break;
        }

        if (found==false) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");
      }
    }
#elif _MESH_DIMENSION_ == 2
    for (int nface=0;nface<4;nface++) for (int iFace=0;iFace<2;iFace++) {
      int nfaceNeib,iFaceNeib;
      bool found;

      neibNode=startNode->GetNeibFace(nface,iFace,0);

      if (neibNode!=NULL) {
        found=false;
        if ((neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_)||(abs(startNode->RefinmentLevel-neibNode->RefinmentLevel)>1)) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");

        for (nfaceNeib=0;nfaceNeib<4;nfaceNeib++) if (found==false) for (iFaceNeib=0;iFaceNeib<2;iFaceNeib++) if (neibNode->GetNeibFace(nfaceNeib,iFaceNeib,0)==startNode) {
          found=true;
          break;
        }

        if (found==false) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");
      }
    }

#elif _MESH_DIMENSION_ == 3
    for (int nface=0;nface<6;nface++) for (int iFace=0;iFace<2;iFace++) for (int jFace=0;jFace<2;jFace++) {
      int nfaceNeib,iFaceNeib,jFaceNeib;
      bool found;

      neibNode=startNode->GetNeibFace(nface,iFace,jFace);

      if (neibNode!=NULL) {
        found=false;
        if ((neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_)||(abs(startNode->RefinmentLevel-neibNode->RefinmentLevel)>1)) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");

        for (nfaceNeib=0;nfaceNeib<6;nfaceNeib++) if (found==false) for (iFaceNeib=0;iFaceNeib<2;iFaceNeib++) for (jFaceNeib=0;jFaceNeib<2;jFaceNeib++) {
          if (neibNode->GetNeibFace(nfaceNeib,iFaceNeib,jFaceNeib)==startNode) {
            found=true;
            break;
          }
        }

        if (found==false) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");
      }
    }
#endif


    ////check consistency of the nodes': connection though the edges
#if _MESH_DIMENSION_ == 3
    for (int nedge=0;nedge<12;nedge++) for (int iEdge=0;iEdge<2;iEdge++) {
      int nedgeNeib,iEdgeNeib;
      bool found;

      neibNode=startNode->GetNeibEdge(nedge,iEdge);

      if (neibNode!=NULL) {
        found=false;

        for (nedgeNeib=0;nedgeNeib<12;nedgeNeib++) if (found==false) for (iEdgeNeib=0;iEdgeNeib<2;iEdgeNeib++) if (neibNode->GetNeibEdge(nedgeNeib,iEdgeNeib)==startNode) {
          found=true;
          break;
        }

        if (found==false) exit(__LINE__,__FILE__,"Error: node's connection is not consistent");
      }
    }

#endif

if ((AllowBlockAllocation==true)&&(startNode->block!=NULL)) {
     //check the neibours 
     for (k=kMin;k<kMax;k++) for (j=jMin;j<jMax;j++) for (i=iMin;i<iMax;i++)  if ((neibNode=getNeibNode(i,j,k,startNode))!=NULL) {
       if (abs(neibNode->RefinmentLevel-startNode->RefinmentLevel)>1) {
         exit(__LINE__,__FILE__,"The resolution change exseeds the limit"); 
       }

       if (neibNode->block==NULL) continue;

       //check if the startNode and neibNode share a cornerNode 
       for (flag=false,iStartNode=0;(iStartNode<nCornerNodes)&&(flag==false);iStartNode++) for (iNeibNode=0;iNeibNode<nCornerNodes;iNeibNode++) {
         if (startNode->block->GetCornerNode(iStartNode)==neibNode->block->GetCornerNode(iNeibNode)) {
            flag=true;
            break;
         }
       }  

       if (flag==false) {
         long int nd,ii,jj,kk,idim;

         fprintf(DiagnospticMessageStream,"$PREFIX:ERROR: Two neiubouring blocks doesn't share nodes\n");
         fprintf(DiagnospticMessageStream,"$PREFIX:Block 1: ID=%ld\nCorner Nodes:\n",startNode->Temp_ID);

         for (ii=0;ii<iMax;ii++) for (jj=0;jj<jMax;jj++) for (kk=0;kk<kMax;kk++) {
           nd=getCornerNodeLocalNumber(ii*_BLOCK_CELLS_X_,jj*_BLOCK_CELLS_Y_,kk*_BLOCK_CELLS_Z_);
  
/*           *DiagnospticMessageStream << "$PREFIX:(i,j,k)=" << ii <<"," << jj << "," << kk << " (nd=" << nd << "), Temp_ID=" << startNode->block->GetCornerNode(nd)->Temp_ID << " x=";
           for (idim=0;idim<_MESH_DIMENSION_;idim++) *DiagnospticMessageStream << startNode->block->GetCornerNode(nd)->GetX()[idim] << " ";
           *DiagnospticMessageStream << std::endl;*/

           fprintf(DiagnospticMessageStream,"$PREFIX:(i,j,k)=%ld,%ld,%ld (nd=%ld), Temp_ID=%ld x=",ii,jj,kk,nd,startNode->block->GetCornerNode(nd)->Temp_ID);
           for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"%e ",startNode->block->GetCornerNode(nd)->GetX()[idim]);
           fprintf(DiagnospticMessageStream,"\n");

         }

         fprintf(DiagnospticMessageStream,"$PREFIX:Block 2: ID=%ld\nCorner Nodes:",neibNode->Temp_ID );
         fprintf(DiagnospticMessageStream,"neib coordinates=%i, %i, %i\n",i,j,k);

         for (ii=0;ii<iMax;ii++) for (jj=0;jj<jMax;jj++) for (kk=0;kk<kMax;kk++) {
           nd=getCornerNodeLocalNumber(ii*_BLOCK_CELLS_X_,jj*_BLOCK_CELLS_Y_,kk*_BLOCK_CELLS_Z_);
  
/*           *DiagnospticMessageStream << "$PREFIX:(i,j,k)=" << ii <<"," << jj << "," << kk << " (nd=" << nd << "), Temp_ID=" << neibNode->block->GetCornerNode(nd)->Temp_ID << " x=";
           for (idim=0;idim<_MESH_DIMENSION_;idim++) *DiagnospticMessageStream << neibNode->block->GetCornerNode(nd)->GetX()[idim] << " ";
           *DiagnospticMessageStream << std::endl;*/


           fprintf(DiagnospticMessageStream,"$PREFIX:(i,j,k)=%ld,%ld,%ld (nd=%ld), Temp_ID=%ld x=",ii,jj,kk,nd,neibNode->block->GetCornerNode(nd)->Temp_ID);
           for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"%e ",neibNode->block->GetCornerNode(nd)->GetX()[idim]);
           fprintf(DiagnospticMessageStream,"\n");
         }

         exit(__LINE__,__FILE__,"The blocks do not share a node"); 
       }

     }

     //check the value of the corner nodes in the block
     int ii,jj,kk;
     long int nd;
     double *xmin,*xmax,dx[3]={0.0,0.0,0.0},*xnode;

     #if _MESH_DIMENSION_ == 1
     static const int iiCornerMin=-_GHOST_CELLS_X_,iiCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjCornerMin=0,jjCornerMax=0,kkCornerMin=0,kkCornerMax=0; 
     #elif _MESH_DIMENSION_ == 2
     static const int iiCornerMin=-_GHOST_CELLS_X_,iiCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjCornerMin=-_GHOST_CELLS_Y_,jjCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkCornerMin=0,kkCornerMax=0;
     #else
     static const int iiCornerMin=-_GHOST_CELLS_X_,iiCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjCornerMin=-_GHOST_CELLS_Y_,jjCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkCornerMin=-_GHOST_CELLS_Z_,kkCornerMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
     #endif

     xmin=startNode->xmin,xmax=startNode->xmax;
     dx[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
     if (_MESH_DIMENSION_>=2) dx[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_; 
     if (_MESH_DIMENSION_==3) dx[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

     for (kk=kkCornerMin;kk<=kkCornerMax;kk++) for (jj=jjCornerMin;jj<=jjCornerMax;jj++)  for (ii=iiCornerMin;ii<=iiCornerMax;ii++) {
       nd=getCornerNodeLocalNumber(ii,jj,kk);

       if (startNode->block->GetCornerNode(nd)==NULL) {
         if ((ii>=0)&&(ii<=_BLOCK_CELLS_X_)&&(jj>=0)&&(jj<=_BLOCK_CELLS_Y_)&&(kk>=0)&&(kk<=_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
         else {
           int iNeib=0,jNeib=0,kNeib=0;
 
           if (ii<0) iNeib=-1;
           if (ii>_BLOCK_CELLS_X_) iNeib=1;
 
           if (jj<0) jNeib=-1;
           if (jj>_BLOCK_CELLS_Y_) jNeib=1;
 
           if (kk<0) kNeib=-1;
           if (kk>_BLOCK_CELLS_Z_) kNeib=1; 
 
           if (getNeibNode(iNeib,jNeib,kNeib,startNode)!=NULL) {
             exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
           }
           else continue; 
         }
       } 

       xnode=startNode->block->GetCornerNode(nd)->GetX();

       if (fabs(xmin[0]+dx[0]*ii-xnode[0])>EPS) exit(__LINE__,__FILE__,"The position of the block's corner nodes is not correct"); 
       if (_MESH_DIMENSION_>=2) if (fabs(xmin[1]+dx[1]*jj-xnode[1])>EPS) exit(__LINE__,__FILE__,"The position of the block's corner nodes is not correct"); 
       if (_MESH_DIMENSION_==3) if (fabs(xmin[2]+dx[2]*kk-xnode[2])>EPS) exit(__LINE__,__FILE__,"The position of the block's corner nodes is not correct");
     }  

     //check if the nodes that the block shares with its neibours the have the same coordinates points to the same location in physical memory
     int ioffset,joffset,koffset/*,iNeibNode*/,jNeibNode,kNeibNode,iDownNode,jDownNode,kDownNode;
     cTreeNodeAMR<cBlockAMR> *upNode; //,*downNode;
     double *xStartNode,*xNeibNode; 
     long int ndStartNode,ndNeibNode;
     int iiDownCornerNode=0,jjDownCornerNode=0,kkDownCornerNode=0;

     #if _MESH_DIMENSION_ == 1
     static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
     static const int jjMin=0,jjMax=0;
     static const int kkMin=0,kkMax=0;
//     static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=0;
     static const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=0,jNeibNodeMax=0,kNeibNodeMin=0,kNeibNodeMax=0;
     #elif _MESH_DIMENSION_ == 2 
     static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
     static const int jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
     static const int kkMin=0,kkMax=0;
//     static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=0;
     static const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=-1,jNeibNodeMax=1,kNeibNodeMin=0,kNeibNodeMax=0;
     #else 
     static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
     static const int jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
     static const int kkMin=-_GHOST_CELLS_Z_,kkMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
//     static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=2;
     static const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=-1,jNeibNodeMax=1,kNeibNodeMin=-1,kNeibNodeMax=1;
     #endif 

     for (iNeibNode=iNeibNodeMin;iNeibNode<=iNeibNodeMax;iNeibNode++) for (jNeibNode=jNeibNodeMin;jNeibNode<=jNeibNodeMax;jNeibNode++) for (kNeibNode=kNeibNodeMin;kNeibNode<=kNeibNodeMax;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel==startNode->RefinmentLevel) {
        ioffset=_BLOCK_CELLS_X_*iNeibNode;
        joffset=_BLOCK_CELLS_Y_*jNeibNode;
        koffset=_BLOCK_CELLS_Z_*kNeibNode;

        if (neibNode->block==NULL) continue;

        for (ii=iiMin;ii<=iiMax;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<=kkMax;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
              ndNeibNode=getCornerNodeLocalNumber(ii,jj,kk);
              ndStartNode=getCornerNodeLocalNumber(ii+ioffset,jj+joffset,kk+koffset);

              xStartNode=startNode->block->GetCornerNode(ndStartNode)->GetX();
              xNeibNode=neibNode->block->GetCornerNode(ndNeibNode)->GetX();

              if ((xNeibNode==NULL)&&(xStartNode==NULL)) continue;
              else if ((xNeibNode==NULL)||(xStartNode==NULL)) exit(__LINE__,__FILE__,"Error: one of the nodes is not defined");
              else {
                if (xStartNode!=xNeibNode) exit(__LINE__,__FILE__,"Error: two different corner nodes have found that has the same coordinate 'x'");
              }
            }
          }
        }
      }
      else if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
        upNode=neibNode->upNode;

        ioffset=_BLOCK_CELLS_X_*iNeibNode;
        joffset=_BLOCK_CELLS_Y_*jNeibNode;
        koffset=_BLOCK_CELLS_Z_*kNeibNode;

        for (ii=iiMin;ii<=iiMax;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<=kkMax;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
              ndStartNode=getCornerNodeLocalNumber(ii+ioffset,jj+joffset,kk+koffset);
              xStartNode=startNode->block->GetCornerNode(ndStartNode)->GetX();

              for (iDownNode=0;iDownNode<=iDownNodeMax;iDownNode++) for (jDownNode=0;jDownNode<=jDownNodeMax;jDownNode++)  for (kDownNode=0;kDownNode<=kDownNodeMax;kDownNode++) {
                downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

                if (downNode->block==NULL) continue;

                //find the index of the node neibNode->block->cornerNodes(ii,jj,kk) in the downNode 
                iiDownCornerNode=2*ii-_BLOCK_CELLS_X_*iDownNode; 
                if (_MESH_DIMENSION_>=2) jjDownCornerNode=2*jj-_BLOCK_CELLS_Y_*jDownNode; 
                if (_MESH_DIMENSION_==3) kkDownCornerNode=2*kk-_BLOCK_CELLS_Z_*kDownNode;

                if ((-_GHOST_CELLS_X_<=iiDownCornerNode)&&(iiDownCornerNode<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
                  if ((-_GHOST_CELLS_Y_<=jjDownCornerNode)&&(jjDownCornerNode<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
                    if ((-_GHOST_CELLS_Z_<=kkDownCornerNode)&&(kkDownCornerNode<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
                      ndNeibNode=getCornerNodeLocalNumber(iiDownCornerNode,jjDownCornerNode,kkDownCornerNode);
                      xNeibNode=downNode->block->GetCornerNode(ndNeibNode)->GetX();

                      if ((xNeibNode==NULL)&&(xStartNode==NULL)) continue;
                      else if ((xNeibNode==NULL)||(xStartNode==NULL)) exit(__LINE__,__FILE__,"Error: one of the nodes is not defined");
                      else {
                        if (xStartNode!=xNeibNode) exit(__LINE__,__FILE__,"Error: two different corner nodes have found that has the same coordinate 'x'");
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

    #if _AMR_CENTER_NODE_  == _ON_AMR_MESH_
    cCenterNode *neibCenterNode,*startCenterNode;


    //check if all centeral nodes (including in 'ghost' cells are defined)
    for (k=((_MESH_DIMENSION_==3) ? -_GHOST_CELLS_Z_ : 0);k<((_MESH_DIMENSION_==3) ? _BLOCK_CELLS_Z_+_GHOST_CELLS_Z_ : 1);k++) for (j=((_MESH_DIMENSION_>=2) ? -_GHOST_CELLS_Y_ : 0);j<((_MESH_DIMENSION_>=2) ? _BLOCK_CELLS_Y_+_GHOST_CELLS_Y_ : 1);j++) for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) { 
      nd=getCenterNodeLocalNumber(i,j,k);

      if (startNode->block->GetCenterNode(nd)==NULL) {
        if ((i>=0)&&(i<_BLOCK_CELLS_X_)&&(j>=0)&&(j<_BLOCK_CELLS_Y_)&&(k>=0)&&(k<_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
        else {
          int iNeib=0,jNeib=0,kNeib=0;

          if (i<0) iNeib=-1;
          if (i>=_BLOCK_CELLS_X_) iNeib=1;

          if (j<0) jNeib=-1;
          if (j>=_BLOCK_CELLS_Y_) jNeib=1;

          if (k<0) kNeib=-1;
          if (k>=_BLOCK_CELLS_Z_) kNeib=1;

          if (getNeibNode(iNeib,jNeib,kNeib,startNode)!=NULL) {
            exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
          }
          else continue;
        }
      }
    }


    //check if the internal nodes of the block are used only ones on the mesh
    for (k=0;k<((_MESH_DIMENSION_==3) ? _BLOCK_CELLS_Z_ : 1);k++) for (j=0;j<((_MESH_DIMENSION_>=2) ? _BLOCK_CELLS_Y_ : 1);j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      nd=getCenterNodeLocalNumber(i,j,k);
      startCenterNode=startNode->block->GetCenterNode(nd);

      if (startCenterNode==NULL) exit(__LINE__,__FILE__,"Error: a center node is not defined"); 

      if (startCenterNode->nodeDescriptor.nodeProcessedFlag==_AMR_TRUE_) exit(__LINE__,__FILE__,"Error: a center node is used twice on the mesh");
      else startCenterNode->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_; 
    }

    //check the location of the center nodes: only blocks that have the same resolution level can share 'center nodes'
     for (iNeibNode=iNeibNodeMin;iNeibNode<=iNeibNodeMax;iNeibNode++) for (jNeibNode=jNeibNodeMin;jNeibNode<=jNeibNodeMax;jNeibNode++) for (kNeibNode=kNeibNodeMin;kNeibNode<=kNeibNodeMax;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel==startNode->RefinmentLevel) {
        ioffset=_BLOCK_CELLS_X_*iNeibNode;
        joffset=_BLOCK_CELLS_Y_*jNeibNode;
        koffset=_BLOCK_CELLS_Z_*kNeibNode;

        if (neibNode->block==NULL) continue;

        for (ii=iiMin;ii<iiMax;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<jjMax;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<kkMax;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
              ndNeibNode=getCenterNodeLocalNumber(ii,jj,kk);
              ndStartNode=getCenterNodeLocalNumber(ii+ioffset,jj+joffset,kk+koffset);

              startCenterNode=startNode->block->GetCenterNode(ndStartNode);
              neibCenterNode=neibNode->block->GetCenterNode(ndNeibNode);


              if ((neibCenterNode==NULL)&&(startCenterNode==NULL)) continue;
              else if ((neibCenterNode==NULL)||(startCenterNode==NULL)) exit(__LINE__,__FILE__,"Error: one of the nodes is not defined");
              else {
                if (startCenterNode!=neibCenterNode) {
//                	long int nd;
                	int iii[3]={0,0,0},idim;

                	//recalculate the position of the neibBlock
                	neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode);

                	//recalculate the indexes of the 'center nodes'
                	nd=findCenterNodeIndex(startCenterNode->GetX(),iii[0],iii[1],iii[2],startNode);
                	fprintf(DiagnospticMessageStream,"$PREFIX:The center node is not found:x=");
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"%e  ",startCenterNode->GetX()[idim]);
                	fprintf(DiagnospticMessageStream,"; (i,j,k,nd)=");
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"%i  ",iii[idim]);
                	fprintf(DiagnospticMessageStream," Local node number=%ld",nd);
                	fprintf(DiagnospticMessageStream," Center node Temp_ID=%ld",startCenterNode->Temp_ID);
                	fprintf(DiagnospticMessageStream," Block Temp_ID=%ld",startNode->Temp_ID);
                	fprintf(DiagnospticMessageStream,"\n");

                	//recalculate the indexes of the 'center nodes'
                	nd=findCenterNodeIndex(neibCenterNode->GetX(),iii[0],iii[1],iii[2],neibNode);
                	fprintf(DiagnospticMessageStream,"$PREFIX:The center node is not found:x=");
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"%e  ",neibCenterNode->GetX()[idim]);
                	fprintf(DiagnospticMessageStream,"; (i,j,k,nd)=");
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"%i  ",iii[idim]);
                	fprintf(DiagnospticMessageStream," Local node number=%ld",nd);
                	fprintf(DiagnospticMessageStream," Center node Temp_ID=%ld",neibCenterNode->Temp_ID);
                	fprintf(DiagnospticMessageStream," Block Temp_ID=%ld",neibNode->Temp_ID);
                	fprintf(DiagnospticMessageStream,"\n");


                	exit(__LINE__,__FILE__,"Error: two different center nodes have found that has the same coordinate 'x'");
                }



              }
            }
          }
        }

      }
    }

    #endif

  }
  }
  else {
    //check the downBlocks
    for (kDownBlock=0;kDownBlock<=kDownNodeMax;kDownBlock++) for (jDownBlock=0;jDownBlock<=jDownNodeMax;jDownBlock++) for (iDownBlock=0;iDownBlock<=iDownNodeMax;iDownBlock++) {
      if ((downNode=startNode->downNode[iDownBlock+2*(jDownBlock+2*kDownBlock)])!=NULL) checkMeshConsistency(downNode);  
    }
  }
}

  

//=================================================================================



/*

void collectneibCornerNodes_1D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)]) {
  int nface,ii,nDownNode;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode; 
  long int nd;

  //collect the nodes from the all upNode->downNodes[:]  
  static const int iiOffset[2]={0,_TOTAL_BLOCK_CELLS_X_}; 
  static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_; 
  
  upNode=startNode->upNode;

  for (nDownNode=0,nDownNode<2;nDownNode++) {
    downNode=upNode->downNode[nDownNode];

    for (ii=iiMin;ii<=iiMax;ii++) {
      nd=getCornerNodeLocalNumber(ii,0,0);
      newCornerNodeMap[ii+iiOffset[nDownNode]-iiMin][0][0]=downNode->cornerNodes[nd];    
    }
  }

  //search the neibours
  static int iiNeibOffset[2]=[-2*_BLOCK_CELLS_X_,2*_BLOCK_CELLS_X_];

  for (nface=0;nface=2*_MESH_DIMENSION_;nface++) if ((neibNode=upNode->neibNode[nface])!=NULL) {
    if (neibNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      int ioffset=iiNeibOffset[nface]; 

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=2*ii+ioffset)&&(2*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) { 
        nd=getCornerNodeLocalNumber(ii,0,0);
        newCornerNodeMap[2*ii+ioffset-iiMin][0][0]=downNode->cornerNodes[nd];
      }
    }
    else {
      int ioffset=iiNeibOffset[nface]+nDownNode*_BLOCK_CELLS_X_;

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {                        
        nd=getCornerNodeLocalNumber(ii,0,0);
        newCornerNodeMap[ii+ioffset-iiMin][0][0]=downNode->cornerNodes[nd];
      }
    }
  }
}

*/

/*
void collectNeibCornerNodes_2D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)]) {
  int nface,ii,jj,iDownNode,jDownNode,ioffset,joffset;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode,*neibNode;
  long int nd;

  //collect the nodes from the all upNode->downNodes[:]
  static const int iiOffset[2]={0,_BLOCK_CELLS_X_};
  static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;

  static const int jjOffset[2]={0,_BLOCK_CELLS_Y_};
  static const int jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
  
  upNode=startNode->upNode;

  for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++) {
    downNode=upNode->downNode[iDownNode+2*jDownNode];

    for (ii=iiMin;ii<=iiMax;ii++) for (jj=jjMin;jj<=jjMax;jj++) {
      nd=getCornerNodeLocalNumber(ii,jj,0);
      newCornerNodeMap[ii+iiOffset[iDownNode]+2*_GHOST_CELLS_X_][jj+jjOffset[jDownNode]+2*_GHOST_CELLS_Y_][0]=downNode->block->cornerNodes[nd];
    }
  }

  //search the neibours
  static const int _X_FACE_NEIB_OFFSET_[4]={-2*_BLOCK_CELLS_X_,0,2*_BLOCK_CELLS_X_,0};
  static const int _Y_FACE_NEIB_OFFSET_[4]={0,-2*_BLOCK_CELLS_Y_,0,2*_BLOCK_CELLS_Y_};

  for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=upNode->neibNode[nface])!=NULL) { 
    if (neibNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      ioffset=_X_FACE_NEIB_OFFSET_[nface],joffset=_Y_FACE_NEIB_OFFSET_[nface]; 

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=2*ii+ioffset)&&(2*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
        for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=2*jj+joffset)&&(2*jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
          nd=getCornerNodeLocalNumber(ii,jj,0);
          newCornerNodeMap[2*ii+ioffset+2*_GHOST_CELLS_X_][2*jj+joffset+2*_GHOST_CELLS_Y_][0]=neibNode->block->cornerNodes[nd];
        }
      }
    }
    else {
      for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++) {
        downNode=neibNode->downNode[iDownNode+2*jDownNode];
        ioffset=_X_FACE_NEIB_OFFSET_[nface]+iDownNode*_BLOCK_CELLS_X_;
        joffset=_Y_FACE_NEIB_OFFSET_[nface]+jDownNode*_BLOCK_CELLS_Y_; 

        for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
            nd=getCornerNodeLocalNumber(ii,jj,0);
            newCornerNodeMap[ii+ioffset+2*_GHOST_CELLS_X_][jj+joffset+2*_GHOST_CELLS_Y_][0]=downNode->block->cornerNodes[nd];
          }
        }
      }
    }
  } 

  
  //search the corner neibours 
  int stepFactor,iNeibNode,jNeibNode; 

  for (iNeibNode=-1;iNeibNode<2;iNeibNode+=2) for (jNeibNode=-1;jNeibNode<2;jNeibNode+=2) if (abs(iNeibNode)+abs(jNeibNode)>1) if ((neibNode=getNeibNode_2D(iNeibNode,jNeibNode,upNode))!=NULL) {
    if (neibNode->RefinmentLevel==upNode->RefinmentLevel) {
      stepFactor=2;
      ioffset=_BLOCK_CELLS_X_*((iNeibNode==-1) ? -2 : 2);
      joffset=_BLOCK_CELLS_Y_*((jNeibNode==-1) ? -2 : 2);
    }
    else {
      stepFactor=1;
      ioffset=_BLOCK_CELLS_X_*((iNeibNode==-1) ? -1 : 2);
      joffset=_BLOCK_CELLS_Y_*((jNeibNode==-1) ? -1 : 2);
    } 

    for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=stepFactor*ii+ioffset)&&(stepFactor*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
      for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=stepFactor*jj+joffset)&&(stepFactor*jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
        nd=getCornerNodeLocalNumber(ii,jj,0);
        newCornerNodeMap[stepFactor*ii+ioffset+2*_GHOST_CELLS_X_][stepFactor*jj+joffset+2*_GHOST_CELLS_Y_][0]=neibNode->block->cornerNodes[nd];
      }
    }
  }

  
}

*/

void collectNeibCornerNodes_3D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)]) {
  int ii,jj,kk,iDownNode,jDownNode,kDownNode,ioffset,joffset,koffset;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode,*neibNode;
  long int nd;

  for (ii=0;ii<1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_);ii++) for (jj=0;jj<1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_);jj++) for (kk=0;kk<1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_);kk++) newCornerNodeMap[ii][jj][kk]=NULL;

  //collect the nodes from the all upNode->downNodes[:]
  static const int iiOffset[2]={0,_BLOCK_CELLS_X_};
  static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;

  static const int jjOffset[2]={0,_BLOCK_CELLS_Y_};
  static const int jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;

  static const int kkOffset[2]={0,_BLOCK_CELLS_Z_};
  static const int kkMin=-_GHOST_CELLS_Z_,kkMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;

  upNode=startNode->upNode;

  for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++) for (kDownNode=0;kDownNode<2;kDownNode++) {
    downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

    for (ii=iiMin;ii<=iiMax;ii++) for (jj=jjMin;jj<=jjMax;jj++) for (kk=kkMin;kk<=kkMax;kk++) {
      nd=getCornerNodeLocalNumber(ii,jj,kk);
      newCornerNodeMap[ii+iiOffset[iDownNode]+2*_GHOST_CELLS_X_][jj+jjOffset[jDownNode]+2*_GHOST_CELLS_Y_][kk+kkOffset[kDownNode]+2*_GHOST_CELLS_Z_]=downNode->block->GetCornerNode(nd);
    }
  }

  //search the neibours
  int iNeibNode,jNeibNode,kNeibNode;

  for (iNeibNode=-1;iNeibNode<2;iNeibNode++) for (jNeibNode=-1;jNeibNode<2;jNeibNode++) for (kNeibNode=-1;kNeibNode<2;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,upNode))!=NULL) {
    if (neibNode->RefinmentLevel==upNode->RefinmentLevel) {
      ioffset=_BLOCK_CELLS_X_*2*iNeibNode;
      joffset=_BLOCK_CELLS_Y_*2*jNeibNode;
      koffset=_BLOCK_CELLS_Z_*2*kNeibNode;

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=2*ii+ioffset)&&(2*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
        for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=2*jj+joffset)&&(2*jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
          for (kk=kkMin;kk<=kkMax;kk++) if ((-2*_GHOST_CELLS_Z_<=2*kk+koffset)&&(2*kk+koffset<=2*_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)) {
            nd=getCornerNodeLocalNumber(ii,jj,kk);
            newCornerNodeMap[2*ii+ioffset+2*_GHOST_CELLS_X_][2*jj+joffset+2*_GHOST_CELLS_Y_][2*kk+koffset+2*_GHOST_CELLS_Z_]=neibNode->block->GetCornerNode(nd);
          }
        }
      }
    }
    else {
      neibNode=neibNode->upNode;

      for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++)  for (kDownNode=0;kDownNode<2;kDownNode++) {
        downNode=neibNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

        ioffset=_BLOCK_CELLS_X_*(iDownNode+2*iNeibNode);
        joffset=_BLOCK_CELLS_Y_*(jDownNode+2*jNeibNode);
        koffset=_BLOCK_CELLS_Z_*(kDownNode+2*kNeibNode);

        for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<=kkMax;kk++) if ((-2*_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=2*_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)) {
              nd=getCornerNodeLocalNumber(ii,jj,kk);
              newCornerNodeMap[ii+ioffset+2*_GHOST_CELLS_X_][jj+joffset+2*_GHOST_CELLS_Y_][kk+koffset+2*_GHOST_CELLS_Z_]=downNode->block->GetCornerNode(nd);
            }
          }
        }
      }
    }
  }


}


void collectNeibCenterNodes_3D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCenterNode *newCenterNodeMap[_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_][_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_][_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_]) {
  int ii,jj,kk,ioffset,joffset,koffset,iNeibNode,jNeibNode,kNeibNode;
  cTreeNodeAMR<cBlockAMR> *upNode,*neibNode;
  long int nd;

  for (ii=0;ii<_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_;ii++) for (jj=0;jj<_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_;jj++) for (kk=0;kk<_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_;kk++) newCenterNodeMap[ii][jj][kk]=NULL; 

  upNode=startNode->upNode;

  for (iNeibNode=-1;iNeibNode<2;iNeibNode++) for (jNeibNode=-1;jNeibNode<2;jNeibNode++) for (kNeibNode=-1;kNeibNode<2;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,upNode))!=NULL) {
    if (neibNode->RefinmentLevel==upNode->RefinmentLevel) {
      ioffset=_BLOCK_CELLS_X_*iNeibNode;
      joffset=_BLOCK_CELLS_Y_*jNeibNode;
      koffset=_BLOCK_CELLS_Z_*kNeibNode;

      for (ii=-_GHOST_CELLS_X_;ii<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
        for (jj=-_GHOST_CELLS_Y_;jj<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
          for (kk=-_GHOST_CELLS_Z_;kk<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
            nd=getCenterNodeLocalNumber(ii,jj,kk);
            newCenterNodeMap[ii+ioffset+_GHOST_CELLS_X_][jj+joffset+_GHOST_CELLS_Y_][kk+koffset+_GHOST_CELLS_Z_]=neibNode->block->GetCenterNode(nd);
          }
        }
      }
    }
  }
}
  
//===================================================================


class cFraction {
public:
  long int Nominator,Denominator;

  cFraction () {
    Nominator=0.0,Denominator=1.0;
  }

  inline void Simplify() {
    int n,t;

    static const int nPrimeNumberList=168;
    static const int PrimeNumberList[nPrimeNumberList]={2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
        149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
        349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569,
        571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
        809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};


    //check if the fraction can be simplified
    for (n=0;n<nPrimeNumberList;n++) {
      t=PrimeNumberList[n];

      while ((Nominator%t==0)&&(Denominator%t==0)) Nominator/=t,Denominator/=t;

      if ((Nominator<t)||(Denominator<t)) break;
    }
  }



  cFraction operator + (cFraction op2) {
    cFraction temp;

    temp.Nominator=Nominator*op2.Denominator+op2.Nominator*Denominator;
    temp.Denominator=Denominator*op2.Denominator;

    temp.Simplify();
    return temp;
  }

  cFraction operator + (int op2) {
    cFraction temp;

    temp.Nominator=Nominator+op2*Denominator;
    temp.Denominator=Denominator;

    temp.Simplify();
    return temp;
  }

  cFraction operator - (cFraction op2) {
    cFraction temp;

    temp.Nominator=Nominator*op2.Denominator-op2.Nominator*Denominator;
    temp.Denominator=Denominator*op2.Denominator;

    temp.Simplify();
    return temp;
  }

  cFraction operator = (cFraction op2) {
    Nominator=op2.Nominator;
    Denominator=op2.Denominator;

    Simplify();
    return *this;
  }

  cFraction operator += (int op2) {
    Nominator+=op2*Denominator;

    Simplify();
    return *this;
  }

  cFraction operator -= (int op2) {
    Nominator-=op2*Denominator;

    Simplify();
    return *this;
  }

  cFraction operator = (int op2) {
    Nominator=op2;
    Denominator=1;

    Simplify();
    return *this;
  }

  cFraction operator * (int op) {
    cFraction temp;

    temp.Nominator=Nominator*op;
    temp.Denominator=Denominator;

    temp.Simplify();
    return temp;
  }

  cFraction operator / (int op) {
    cFraction temp;

    temp.Nominator=Nominator;
    temp.Denominator=Denominator*op;

    temp.Simplify();
    return temp;
  }

  cFraction operator /= (int op) {
    Denominator*=op;

    if (Nominator%Denominator==0) {
      Nominator/=Denominator;
      Denominator=1;
    }

    Simplify();
    return *this;
  }


  friend bool operator >= (const cFraction& op1, const int& op2) {
    return (op1.Nominator>=op2*op1.Denominator) ? true : false;
  }

  friend bool operator <= (const cFraction& op1, const int& op2) {
    return (op1.Nominator<=op2*op1.Denominator) ? true : false;
  }

};


struct cNeibDescriptor {
//  int iOffset,jOffset,kOffset,rLevelStep;

  cFraction iOffset,jOffset,kOffset,BlockSize;
  cTreeNodeAMR<cBlockAMR>* node;
};

void AddNodeNeighborList(cTreeNodeAMR<cBlockAMR>* neibNode,cNeibDescriptor *NeibList,int& nNeibListCounter, int nNeibListMax, cFraction iOffset, cFraction jOffset, cFraction kOffset,cFraction BlockSize) {
  cTreeNodeAMR<cBlockAMR>* downNode;




//==================   DEBUG ===============

/*
  if (neibNode->Temp_ID==7021) {
    *DiagnospticMessageStream << __LINE__ << std::endl;
  }
*/

  //===============  END DEBUG ============






  if (neibNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    int nIntersection;

    //check if the nodes can intersect the 'startNode' and add it to the list if it can

    nIntersection=0;

    if (iOffset.Nominator==0) nIntersection++;
    if ((iOffset.Nominator<0)&&(iOffset+BlockSize*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)>=-_GHOST_CELLS_X_)) nIntersection++;
    if ((iOffset.Nominator>0)&&(iOffset-BlockSize*_GHOST_CELLS_X_<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) nIntersection++;

    if (_MESH_DIMENSION_>1) {
      if (jOffset.Nominator==0) nIntersection++;
      if ((jOffset.Nominator<0)&&(jOffset+BlockSize*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)>=-_GHOST_CELLS_Y_)) nIntersection++;
      if ((jOffset.Nominator>0)&&(jOffset-BlockSize*_GHOST_CELLS_Y_<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) nIntersection++;
    } else nIntersection++;

    if (_MESH_DIMENSION_>2) {
      if (kOffset.Nominator==0) nIntersection++;
      if ((kOffset.Nominator<0)&&(kOffset+BlockSize*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)>=-_GHOST_CELLS_Z_)) nIntersection++;
      if ((kOffset.Nominator>0)&&(kOffset-BlockSize*_GHOST_CELLS_Z_<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) nIntersection++;
    } else nIntersection++;

    if ((nIntersection==3)&&(neibNode->block!=NULL)) {
      //the neibNode can intersect the 'startNode' - check if the node is not in the list
      bool found=false;
      int i;

      for (i=0;i<nNeibListCounter;i++) if (NeibList[i].node==neibNode) {
        found=true;
        break;
      }

      if (found==false) {
        NeibList[nNeibListCounter].iOffset=iOffset;
        NeibList[nNeibListCounter].jOffset=jOffset;
        NeibList[nNeibListCounter].kOffset=kOffset;
        NeibList[nNeibListCounter].BlockSize=BlockSize;
        NeibList[nNeibListCounter].node=neibNode;
        ++nNeibListCounter;
        if (nNeibListCounter>=nNeibListMax) exit(__LINE__,__FILE__,"Error: nNeibListCounter exeeds nNeibListMax -> inxrease the value of nNeibListMax");
      }
    }
  }
  else {
    int nDownNode,iDownNode,jDownNode,kDownNode;
    cFraction iDownOffset,jDownOffset,kDownOffset;

    BlockSize/=2;

#if _MESH_DIMENSION_ == 3
    static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=2;
#elif _MESH_DIMENSION_ == 2
    static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=1;
#elif _MESH_DIMENSION_ == 1
    static const int iDownNodeMax=2,jDownNodeMax=1,kDownNodeMax=1;
#endif

    for (kDownNode=0;kDownNode<kDownNodeMax;kDownNode++) for (jDownNode=0;jDownNode<jDownNodeMax;jDownNode++) for (iDownNode=0;iDownNode<iDownNodeMax;iDownNode++) {
      nDownNode=iDownNode+2*(jDownNode+2*kDownNode);

      if ((downNode=neibNode->downNode[nDownNode])!=NULL) {
//        if (downNode->RefinmentLevel>startNodeRefinmentLevel+1) exit(__LINE__,__FILE__,"Error: some of the nodes can be missed -> the procedure needs to be rewritten!!!!!");

        iDownOffset=iOffset+BlockSize*(iDownNode*_BLOCK_CELLS_X_);
        if (_MESH_DIMENSION_>1) jDownOffset=jOffset+BlockSize*(jDownNode*_BLOCK_CELLS_Y_);
        if (_MESH_DIMENSION_>2) kDownOffset=kOffset+BlockSize*(kDownNode*_BLOCK_CELLS_Z_);

        AddNodeNeighborList(downNode,NeibList,nNeibListCounter,nNeibListMax,iDownOffset,jDownOffset,kDownOffset,BlockSize);
      }
    }
  }


}


void AllocateBlock(cTreeNodeAMR<cBlockAMR> *startNode) {
  int i,j,k,idim; //nDownNode,nDownNodeTemp,idim;
  int ioffset,joffset,koffset;




//#####################  DEBUG #################


static long int nCallCounter=0;

++nCallCounter;

/*
if (startNode->Temp_ID==77) {
// if (nCallCounter==34) {
  *DiagnospticMessageStream << __LINE__ << std::endl;
}
*/

//################### END DEBUG  #############3



  //the block cannot be allocated twice AND the block's allocation must be permitted
  if ((AllowBlockAllocation==false)||(startNode->block!=NULL)) return;
  nMeshModificationCounter++,meshModifiedFlag=true,meshModifiedFlag_CountMeshElements=true;


  //collect all nodes that can intersect 'startNode'
  const int nNeibListMax=1250;
  int nNeibListCounter=0,nlist;



  cNeibDescriptor *neibptr,NeibList[nNeibListMax];


  cTreeNodeAMR<cBlockAMR> *t,*searchNode=startNode;
  int SearchNodeSize=1; //,nConditionsMet; //nConditionsMet,iSearchOffset=0,jSearchOffset=0,kSearchOffset=0;
  cFraction iOffset,jOffset,kOffset,iSearchOffset,jSearchOffset,kSearchOffset;

  iOffset=0;
  jOffset=0;
  kOffset=0;

  //define the search node

  //1. determine ht enumebr of the resolution levels up where to start the search
  static int StartSearchLevelOffset=-1;

  if (StartSearchLevelOffset==-1) {
    int TotalOffset_X_=0,TotalOffset_Y_=0,TotalOffset_Z_=0,nIntersections;

    for (nIntersections=0,StartSearchLevelOffset=1;StartSearchLevelOffset<=_MAX_REFINMENT_LEVEL_;nIntersections=0,StartSearchLevelOffset++) {
      TotalOffset_X_-=_BLOCK_CELLS_X_*(1<<StartSearchLevelOffset);
      if (TotalOffset_X_+(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)*(1<<StartSearchLevelOffset)>=-_GHOST_CELLS_X_) nIntersections++;

      if (_MESH_DIMENSION_>1) {
        TotalOffset_Y_-=_BLOCK_CELLS_Y_*(1<<StartSearchLevelOffset);
        if (TotalOffset_Y_+(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)*(1<<StartSearchLevelOffset)>=-_GHOST_CELLS_Y_) nIntersections++;
      }

      if (_MESH_DIMENSION_>2) {
        TotalOffset_Z_-=_BLOCK_CELLS_Z_*(1<<StartSearchLevelOffset);
        if (TotalOffset_Z_+(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)*(1<<StartSearchLevelOffset)>=-_GHOST_CELLS_Z_) nIntersections++;
      }

      if (nIntersections==0) {
        if (StartSearchLevelOffset>1) StartSearchLevelOffset--;
        break;
      }
    }
  }

  //determine the node to start the search
//  while (searchNode->upNode!=NULL) {

  for (int nSearchLevelUp=0;(nSearchLevelUp<StartSearchLevelOffset)&&(searchNode->upNode!=NULL);nSearchLevelUp++) {
//    nConditionsMet=0;

    /*
    if (_BLOCK_CELLS_X_*(1<<(startNode->RefinmentLevel-searchNode->RefinmentLevel))>_BLOCK_CELLS_X_+6*_GHOST_CELLS_X_) nConditionsMet++;
    if ((_MESH_DIMENSION_<2)||(_BLOCK_CELLS_Y_*(1<<(startNode->RefinmentLevel-searchNode->RefinmentLevel))>_BLOCK_CELLS_Y_+6*_GHOST_CELLS_Y_)) nConditionsMet++;
    if ((_MESH_DIMENSION_<3)||(_BLOCK_CELLS_Z_*(1<<(startNode->RefinmentLevel-searchNode->RefinmentLevel))>_BLOCK_CELLS_Z_+6*_GHOST_CELLS_Z_)) nConditionsMet++;

    if (nConditionsMet==3) break;
*/

    t=searchNode;
    searchNode=searchNode->upNode;

    //determine the downNode coordinates
    for (k=0;k<((_MESH_DIMENSION_==3) ? 2 : 1);k++) for (j=0;j<((_MESH_DIMENSION_>1) ? 2 : 1);j++) for (i=0;i<2;i++) if (searchNode->downNode[i+2*(j+2*k)]==t) {
      iOffset-=i*_BLOCK_CELLS_X_*SearchNodeSize;
      if (_MESH_DIMENSION_>1) jOffset-=j*_BLOCK_CELLS_Y_*SearchNodeSize;
      if (_MESH_DIMENSION_>2) kOffset-=k*_BLOCK_CELLS_Z_*SearchNodeSize;
    }

    SearchNodeSize*=2;
  }


  double xSearchNode[_MESH_DIMENSION_];
  cFraction SearchStartBlockSize;

  SearchStartBlockSize=SearchNodeSize;

#if _MESH_DIMENSION_ == 3
  static const int iSearchMin=-1,iSearchMax=1,jSearchMin=-1,jSearchMax=1,kSearchMin=-1,kSearchMax=1;
#elif _MESH_DIMENSION_ == 2
  static const int iSearchMin=-1,iSearchMax=1,jSearchMin=-1,jSearchMax=1,kSearchMin=0,kSearchMax=0;
#else
  static const int iSearchMin=-1,iSearchMax=1,jSearchMin=0,jSearchMax=0,kSearchMin=0,kSearchMax=0;
#endif

  for (k=kSearchMin;k<=kSearchMax;k++) {
    if (_MESH_DIMENSION_==3) {
      xSearchNode[2]=searchNode->xmin[2]+(searchNode->xmax[2]-searchNode->xmin[2])*(0.5+k);
      kSearchOffset=kOffset+_BLOCK_CELLS_Z_*k*SearchNodeSize;
    }

    for (j=jSearchMin;j<=jSearchMax;j++) {
      if (_MESH_DIMENSION_>=2) {
        xSearchNode[1]=searchNode->xmin[1]+(searchNode->xmax[1]-searchNode->xmin[1])*(0.5+j);
        jSearchOffset=jOffset+_BLOCK_CELLS_Y_*j*SearchNodeSize;
      }

      for (i=iSearchMin;i<=iSearchMax;i++) {
        xSearchNode[0]=searchNode->xmin[0]+(searchNode->xmax[0]-searchNode->xmin[0])*(0.5+i);
        iSearchOffset=iOffset+_BLOCK_CELLS_X_*i*SearchNodeSize;

        t=findTreeNodeLimitedResolutionLevel(xSearchNode,searchNode->RefinmentLevel,searchNode);
        if (t==NULL) continue;
        if (t->RefinmentLevel>searchNode->RefinmentLevel) exit(__LINE__,__FILE__,"Error: something is wrong in the serach procedure");

        AddNodeNeighborList(t,NeibList,nNeibListCounter,nNeibListMax,iSearchOffset,jSearchOffset,kSearchOffset,SearchStartBlockSize);
      }
    }
  }



  /*
  //add the upNode to the list
  cTreeNodeAMR<cBlockAMR> *upNode=startNode->upNode;

  if (upNode!=NULL) if (upNode->block!=NULL) {
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (upNode->downNode[nDownNode]==startNode) break;

    if (nDownNode==1<<_MESH_DIMENSION_) exit(__LINE__,__FILE__,"Error: the mesh is inconsistent");

    k=nDownNode/4;
    nDownNode-=4*k;

    j=nDownNode/2;
    nDownNode-=2*j;

    i=nDownNode;

    NeibList[nNeibListCounter].iOffset=-2*i*_BLOCK_CELLS_X_;
    NeibList[nNeibListCounter].jOffset=-2*j*_BLOCK_CELLS_Y_;
    NeibList[nNeibListCounter].kOffset=-2*k*_BLOCK_CELLS_Z_;
    NeibList[nNeibListCounter].rLevelStep=4;
    NeibList[nNeibListCounter].node=upNode;
    ++nNeibListCounter;
  }

  //add the down nodes to the list
  for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) if (startNode->downNode[nDownNode]->block!=NULL) {
    nDownNodeTemp=nDownNode;
    k=nDownNodeTemp/4;
    nDownNodeTemp-=4*k;

    j=nDownNodeTemp/2;
    nDownNodeTemp-=2*j;

    i=nDownNodeTemp;

    NeibList[nNeibListCounter].iOffset=i*_BLOCK_CELLS_X_;
    NeibList[nNeibListCounter].jOffset=j*_BLOCK_CELLS_Y_;
    NeibList[nNeibListCounter].kOffset=k*_BLOCK_CELLS_Z_;
    NeibList[nNeibListCounter].rLevelStep=1;
    NeibList[nNeibListCounter].node=startNode->downNode[nDownNode];
    ++nNeibListCounter;
  }

  //add neighbors to the list
  int iNeibNode,jNeibNode,kNeibNode;
  cTreeNodeAMR<cBlockAMR> *neibNode;
  bool found;

#if _MESH_DIMENSION_ == 3
  const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=-1,jNeibNodeMax=1,kNeibNodeMin=-1,kNeibNodeMax=1;
#elif _MESH_DIMENSION_ == 2
  const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=-1,jNeibNodeMax=1,kNeibNodeMin=0,kNeibNodeMax=0;
#elif _MESH_DIMENSION_ == 1
  const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=0,jNeibNodeMax=0,kNeibNodeMin=0,kNeibNodeMax=0;
#endif

  for (kNeibNode=kNeibNodeMin;kNeibNode<=kNeibNodeMax;kNeibNode++) {
    for (jNeibNode=jNeibNodeMin;jNeibNode<=jNeibNodeMax;jNeibNode++) {
      for (iNeibNode=iNeibNodeMin;iNeibNode<=iNeibNodeMax;iNeibNode++) if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {

        //check if 'neibNode' is not in the list already
        for (found=false,i=0;i<nNeibListCounter;i++) if (NeibList[i].node==neibNode) {
          found=true;
          break;
        }

        //add the 'neibNode' to the list
        if (found==false) {
          if (neibNode->RefinmentLevel==startNode->RefinmentLevel) {
            if (neibNode->block!=NULL) {
              NeibList[nNeibListCounter].iOffset=2*iNeibNode*_BLOCK_CELLS_X_;
              NeibList[nNeibListCounter].jOffset=2*jNeibNode*_BLOCK_CELLS_Y_;
              NeibList[nNeibListCounter].kOffset=2*kNeibNode*_BLOCK_CELLS_Z_;
              NeibList[nNeibListCounter].rLevelStep=2;
              NeibList[nNeibListCounter].node=neibNode;
              ++nNeibListCounter;
              if (nNeibListCounter>=nNeibListMax) exit(__LINE__,__FILE__,"Error: nNeibListCounter exeeds nNeibListMax -> inxrease the value of nNeibListMax");
            }
          }
          else if (neibNode->RefinmentLevel==startNode->RefinmentLevel+1) {
            upNode=neibNode->upNode;
            if (upNode==NULL) continue;

            for (k=0;k<=kNeibNodeMax;k++) { //if ((kNeibNode==0) ||((kNeibNode==-1)&&(k==1)) || ((kNeibNode==1)&&(k==0)) ) {
              for (j=0;j<=jNeibNodeMax;j++) { //if ((jNeibNode==0) ||((jNeibNode==-1)&&(j==1)) || ((jNeibNode==1)&&(j==0)) ) {
                for (i=0;i<=iNeibNodeMax;i++) { //if ((iNeibNode==0) ||((iNeibNode==-1)&&(i==1)) || ((iNeibNode==1)&&(i==0)) ) {
                  if (upNode->downNode[i+2*(j+2*k)]!=NULL) if (upNode->downNode[i+2*(j+2*k)]->block!=NULL) {
                    NeibList[nNeibListCounter].iOffset=_BLOCK_CELLS_X_*(i+2*iNeibNode);
                    NeibList[nNeibListCounter].jOffset=_BLOCK_CELLS_Y_*(j+2*jNeibNode);
                    NeibList[nNeibListCounter].kOffset=_BLOCK_CELLS_Z_*(k+2*kNeibNode);
                    NeibList[nNeibListCounter].rLevelStep=1;
                    NeibList[nNeibListCounter].node=upNode->downNode[i+2*(j+2*k)];
                    ++nNeibListCounter;
                    if (nNeibListCounter>=nNeibListMax) exit(__LINE__,__FILE__,"Error: nNeibListCounter exeeds nNeibListMax -> inxrease the value of nNeibListMax");
                  }

                }
              }
            }

          }
          else if (neibNode->RefinmentLevel==startNode->RefinmentLevel-1) {
            //find all blocks with the lower resolution level that can intersect 'startNode'
            //1. find the coordinates and offset of 'startNode' in 'startNode->upNode'
            //2. Go thought all neib's of 'startNode->upNode' and recover all its neib's that have the same resolution level
            //3. Add to the list those of the neib's that can intersect 'startNode'

            int nIntersection,iStartNodeOffset=-1,jStartNodeOffset=-1,kStartNodeOffset=-1;

            upNode=startNode->upNode;
            if (upNode==NULL) continue;

            for (k=0;k<=kNeibNodeMax;k++) for (j=0;j<=jNeibNodeMax;j++) for (i=0;i<=iNeibNodeMax;i++) {
              nDownNode=i+2*(j+2*k);

              if (upNode->downNode[nDownNode]==startNode) {
                iStartNodeOffset=-2*i*_BLOCK_CELLS_X_;
                jStartNodeOffset=-2*j*_BLOCK_CELLS_Y_;
                kStartNodeOffset=-2*k*_BLOCK_CELLS_Z_;
                break;
              }
            }

            if (iStartNodeOffset==-1) exit(__LINE__,__FILE__,"Error: the mesh is not consistent");

            for (k=kNeibNodeMin;k<=kNeibNodeMax;k++) {
               for (j=jNeibNodeMin;j<=jNeibNodeMax;j++) {
                 for (i=iNeibNodeMin;i<=iNeibNodeMax;i++) if (abs(i)+abs(j)+abs(k)!=0) if ((neibNode=getNeibNode(i,j,k,upNode))!=NULL) if (neibNode->block!=NULL) {
                   if (neibNode->RefinmentLevel==startNode->RefinmentLevel-1) {
                     nIntersection=0;
                     ioffset=4*i*_BLOCK_CELLS_X_+iStartNodeOffset;
                     joffset=4*j*_BLOCK_CELLS_Y_+jStartNodeOffset;
                     koffset=4*k*_BLOCK_CELLS_Z_+kStartNodeOffset;

 //                    if ((ioffset==0)||(ioffset==-4*_BLOCK_CELLS_X_)||(ioffset==2*_BLOCK_CELLS_X_)||(ioffset==-2*_BLOCK_CELLS_X_)) nIntersection++;
 //                    if ((joffset==0)||(joffset==-4*_BLOCK_CELLS_Y_)||(joffset==2*_BLOCK_CELLS_Y_)||(joffset==-2*_BLOCK_CELLS_Y_)) nIntersection++;
 //                    if ((koffset==0)||(koffset==-4*_BLOCK_CELLS_Z_)||(koffset==2*_BLOCK_CELLS_Z_)||(koffset==-2*_BLOCK_CELLS_Z_)) nIntersection++;

 //                    if (nIntersection==3) {
                       //add the block to the list
                       NeibList[nNeibListCounter].iOffset=ioffset;
                       NeibList[nNeibListCounter].jOffset=joffset;
                       NeibList[nNeibListCounter].kOffset=koffset;
                       NeibList[nNeibListCounter].rLevelStep=4;
                       NeibList[nNeibListCounter].node=neibNode;
                       ++nNeibListCounter;
                       if (nNeibListCounter>=nNeibListMax) exit(__LINE__,__FILE__,"Error: nNeibListCounter exeeds nNeibListMax -> inxrease the value of nNeibListMax");
//                     }
                   }

                 }
               }
            }


          }
          else exit(__LINE__,__FILE__,"Error: the mesh is not consistent");

        }

      }
    }
  }



*/
  if (nNeibListCounter>=nNeibListMax) exit(__LINE__,__FILE__,"Error: nNeibListCounter exeeds nNeibListMax -> inxrease the value of nNeibListMax");


  //allocate the new block
  startNode->block=blocks.newElement();
  startNode->block->SetRefinmentLevel(startNode->RefinmentLevel);

  //determine the distance between nodes in the block
  double dxBlock[3]={0.0,0.0,0.0},*xminBlock=startNode->xmin,*xmaxBlock=startNode->xmax,x[3]={0.0,0.0,0.0};

  //init dxBlock
  dxBlock[0]=(xmaxBlock[0]-xminBlock[0])/_BLOCK_CELLS_X_;
  if (_MESH_DIMENSION_>1) dxBlock[1]=(xmaxBlock[1]-xminBlock[1])/_BLOCK_CELLS_Y_; 
  if (_MESH_DIMENSION_>2) dxBlock[2]=(xmaxBlock[2]-xminBlock[2])/_BLOCK_CELLS_Z_; 

  //create the map of the corner nodes
#if _MESH_DIMENSION_ == 3
  static const int iCornerMin=-_GHOST_CELLS_X_,iCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
  static const int jCornerMin=-_GHOST_CELLS_Y_,jCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
  static const int kCornerMin=-_GHOST_CELLS_Z_,kCornerMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
#elif _MESH_DIMENSION_ == 2
  static const int iCornerMin=-_GHOST_CELLS_X_,iCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
  static const int jCornerMin=-_GHOST_CELLS_Y_,jCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
  static const int kCornerMin=0,kCornerMax=0;
#elif _MESH_DIMENSION_ == 1
  static const int iCornerMin=-_GHOST_CELLS_X_,iCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
  static const int jCornerMin=0,jCornerMax=0;
  static const int kCornerMin=0,kCornerMax=0;
#endif

  long int nd;
  cCornerNode *CornerNodeMap[1+iCornerMax-iCornerMin][1+jCornerMax-jCornerMin][1+kCornerMax-kCornerMin],*ptrCornerNode=NULL;

  for (k=0;k<1+kCornerMax-kCornerMin;k++) for (j=0;j<1+jCornerMax-jCornerMin;j++) for (i=0;i<1+iCornerMax-iCornerMin;i++) CornerNodeMap[i][j][k]=NULL;


  for (nlist=0;nlist<nNeibListCounter;nlist++) {
    neibptr=NeibList+nlist;

    for (k=kCornerMin;k<=kCornerMax;k++) {
      kOffset=neibptr->kOffset+neibptr->BlockSize*k;

      if ((kOffset>=kCornerMin)&&(kOffset<=kCornerMax)&&(kOffset.Nominator%kOffset.Denominator==0)) for (j=jCornerMin;j<=jCornerMax;j++) {
        jOffset=neibptr->jOffset+neibptr->BlockSize*j;

        if ((jOffset>=jCornerMin)&&(jOffset<=jCornerMax)&&(jOffset.Nominator%jOffset.Denominator==0)) for (i=iCornerMin;i<=iCornerMax;i++) {
          iOffset=neibptr->iOffset+neibptr->BlockSize*i;

          if ((iOffset>=iCornerMin)&&(iOffset<=iCornerMax)&&(iOffset.Nominator%iOffset.Denominator==0)) {
            ptrCornerNode=neibptr->node->block->GetCornerNode(getCornerNodeLocalNumber(i,j,k));

            ioffset=iOffset.Nominator/iOffset.Denominator;
            joffset=jOffset.Nominator/jOffset.Denominator;
            koffset=kOffset.Nominator/kOffset.Denominator;

            if ((CornerNodeMap[ioffset-iCornerMin][joffset-jCornerMin][koffset-kCornerMin]!=NULL)&&(CornerNodeMap[ioffset-iCornerMin][joffset-jCornerMin][koffset-kCornerMin]!=ptrCornerNode)) exit(__LINE__,__FILE__,"Error: redifinition of the node");

            CornerNodeMap[ioffset-iCornerMin][joffset-jCornerMin][koffset-kCornerMin]=ptrCornerNode;

#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
#if _CHECK_MESH_CONSISTANCY_ == _ON_AMR_MESH_
            if (ptrCornerNode!=NULL) {
              double *xprobe=ptrCornerNode->GetX();

              if ((xprobe[0]+EPS<startNode->xmin[0]-dxBlock[0]*_GHOST_CELLS_X_)||(xprobe[0]-EPS>startNode->xmax[0]+dxBlock[0]*_GHOST_CELLS_X_)) exit(__LINE__,__FILE__,"Error: new node is outside of the block");
              if (_MESH_DIMENSION_>=2) if ((xprobe[1]+EPS<startNode->xmin[1]-dxBlock[1]*_GHOST_CELLS_Y_)||(xprobe[1]-EPS>startNode->xmax[1]+dxBlock[1]*_GHOST_CELLS_Y_)) exit(__LINE__,__FILE__,"Error: new node is outside of the block");
              if (_MESH_DIMENSION_==3) if ((xprobe[2]+EPS<startNode->xmin[2]-dxBlock[2]*_GHOST_CELLS_Z_)||(xprobe[2]-EPS>startNode->xmax[2]+dxBlock[2]*_GHOST_CELLS_Z_)) exit(__LINE__,__FILE__,"Error: new node is outside of the block");
            }
#endif
#endif


          }
        }
      }
    }
  }

  //populate the new block's corner nodes
  bool OutsideDomainFlag;

  for (k=kCornerMin;k<=kCornerMax;k++) for (j=jCornerMin;j<=jCornerMax;j++) for (i=iCornerMin;i<=iCornerMax;i++) {
    nd=getCornerNodeLocalNumber(i,j,k);
    ptrCornerNode=CornerNodeMap[i-iCornerMin][j-jCornerMin][k-kCornerMin];

    if (ptrCornerNode==NULL) {
      //generate the coordinates ofthe new node
      x[0]=xminBlock[0]+i*dxBlock[0];
      if (_MESH_DIMENSION_>1) x[1]=xminBlock[1]+j*dxBlock[1];
      if (_MESH_DIMENSION_>2) x[2]=xminBlock[2]+k*dxBlock[2];

      //check if the node is within the domain
      for (OutsideDomainFlag=false,idim=0;idim<_MESH_DIMENSION_;idim++) if ((x[idim]<xGlobalMin[idim]-EPS)||(xGlobalMax[idim]+EPS<x[idim])) {
        OutsideDomainFlag=true;
        break;
      }

      if (OutsideDomainFlag==true) {
        startNode->block->SetCornerNode(NULL,nd);
        continue;
      }



//==========================    DEBUG ======================
/*
      if ((fabs(x[0]+225)<EPS)&&(fabs(x[1]+225)<EPS)&&(fabs(x[2]+12.5)<EPS)) {
  *DiagnospticMessageStream << __LINE__ << std::endl;
}
*/
//========================== END DEBUG ===========================



      ptrCornerNode=CornerNodes.newElement();
      ptrCornerNode->SetX(x);
    }

    startNode->block->SetCornerNode(ptrCornerNode,nd);
    ptrCornerNode->incrementConnectionCounter();
  }

  //collect the center nodes
#if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
#if _MESH_DIMENSION_ == 3
  static const int iCenterMin=-_GHOST_CELLS_X_,iCenterMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;
  static const int jCenterMin=-_GHOST_CELLS_Y_,jCenterMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1;
  static const int kCenterMin=-_GHOST_CELLS_Z_,kCenterMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_-1;
#elif _MESH_DIMENSION_ == 2
  static const int iCenterMin=-_GHOST_CELLS_X_,iCenterMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;
  static const int jCenterMin=-_GHOST_CELLS_Y_,jCenterMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1;
  static const int kCenterMin=0,kCenterMax=0;
#elif _MESH_DIMENSION_ == 1
  static const int iCenterMin=-_GHOST_CELLS_X_,iCenterMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;
  static const int jCenterMin=0,jCenterMax=0;
  static const int kCenterMin=0,kCenterMax=0;
#endif

  cCenterNode *CenterNodeMap[1+iCenterMax-iCenterMin][1+jCenterMax-jCenterMin][1+kCenterMax-kCenterMin],*ptrCenterNode;

  for (k=0;k<1+kCenterMax-kCenterMin;k++) for (j=0;j<1+jCenterMax-jCenterMin;j++) for (i=0;i<1+iCenterMax-iCenterMin;i++) CenterNodeMap[i][j][k]=NULL;

  for (nlist=0;nlist<nNeibListCounter;nlist++) {
    neibptr=NeibList+nlist;
    if (neibptr->node->RefinmentLevel!=startNode->RefinmentLevel) continue;

    for (k=kCenterMin;k<=kCenterMax;k++) {
      koffset=neibptr->kOffset.Nominator/neibptr->kOffset.Denominator+k;

      if ((koffset>=kCenterMin)&&(koffset<=kCenterMax)) for (j=jCenterMin;j<=jCenterMax;j++) {
        joffset=neibptr->jOffset.Nominator/neibptr->jOffset.Denominator+j;

        if ((joffset>=jCenterMin)&&(joffset<=jCenterMax)) for (i=iCenterMin;i<=iCenterMax;i++) {
          ioffset=neibptr->iOffset.Nominator/neibptr->iOffset.Denominator+i;

          if ((ioffset>=iCenterMin)&&(ioffset<=iCenterMax)) {
            ptrCenterNode=neibptr->node->block->GetCenterNode(getCenterNodeLocalNumber(i,j,k));

            if ((CenterNodeMap[ioffset-iCenterMin][joffset-jCenterMin][koffset-kCenterMin]!=NULL)&&(CenterNodeMap[ioffset-iCenterMin][joffset-jCenterMin][koffset-kCenterMin]!=ptrCenterNode)) exit(__LINE__,__FILE__,"Error: redifinition of the node");

            CenterNodeMap[ioffset-iCenterMin][joffset-jCenterMin][koffset-kCenterMin]=ptrCenterNode;
          }
        }
      }
    }
  }

  //populate the center nodes of the new block
  for (k=kCenterMin;k<=kCenterMax;k++) for (j=jCenterMin;j<=jCenterMax;j++) for (i=iCenterMin;i<=iCenterMax;i++) {
    nd=getCenterNodeLocalNumber(i,j,k);
    ptrCenterNode=CenterNodeMap[i-iCenterMin][j-jCenterMin][k-kCenterMin];
    if (ptrCenterNode==NULL) {
      //generate the coordinates of the new node
      x[0]=xminBlock[0]+(i+0.5)*dxBlock[0];
      if (_MESH_DIMENSION_>1) x[1]=xminBlock[1]+(j+0.5)*dxBlock[1];
      if (_MESH_DIMENSION_>2) x[2]=xminBlock[2]+(k+0.5)*dxBlock[2];

      //check if the node is within the domain
      for (OutsideDomainFlag=false,idim=0;idim<_MESH_DIMENSION_;idim++) if ((x[idim]<xGlobalMin[idim]-EPS)||(xGlobalMax[idim]+EPS<x[idim])) {
        OutsideDomainFlag=true;
        break;
      }

      if (OutsideDomainFlag==true) {
        startNode->block->SetCenterNode(NULL,nd);
        continue;
      }

      ptrCenterNode=CenterNodes.newElement();
      ptrCenterNode->SetX(x);
      ptrCenterNode->Measure=-1.0;
    }

    startNode->block->SetCenterNode(ptrCenterNode,nd);
    ptrCenterNode->incrementConnectionCounter();
  }
#endif
}


void DeallocateBlock(cTreeNodeAMR<cBlockAMR> *startNode) {
  int i,j,k;
  cCornerNode *ptrCornerNode;

  //the block cannot be deallocated twice
  if (startNode->block==NULL) return;
  nMeshModificationCounter++,meshModifiedFlag=true,meshModifiedFlag_CountMeshElements=true;

#if _MESH_DIMENSION_ == 3
  static const int iCornerMin=-_GHOST_CELLS_X_,iCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
  static const int jCornerMin=-_GHOST_CELLS_Y_,jCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
  static const int kCornerMin=-_GHOST_CELLS_Z_,kCornerMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
#elif _MESH_DIMENSION_ == 2
  static const int iCornerMin=-_GHOST_CELLS_X_,iCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
  static const int jCornerMin=-_GHOST_CELLS_Y_,jCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
  static const int kCornerMin=0,kCornerMax=0;
#elif _MESH_DIMENSION_ == 1
  static const int iCornerMin=-_GHOST_CELLS_X_,iCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
  static const int jCornerMin=0,jCornerMax=0;
  static const int kCornerMin=0,kCornerMax=0;
#endif

  for (k=kCornerMin;k<=kCornerMax;k++) for (j=jCornerMin;j<=jCornerMax;j++) for (i=iCornerMin;i<=iCornerMax;i++) {
    ptrCornerNode=startNode->block->GetCornerNode(getCornerNodeLocalNumber(i,j,k));
    if (ptrCornerNode!=NULL) if (ptrCornerNode->decrementConnectionCounter()==0) CornerNodes.deleteElement(ptrCornerNode);
  }

#if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
#if _MESH_DIMENSION_ == 3
  static const int iCenterMin=-_GHOST_CELLS_X_,iCenterMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;
  static const int jCenterMin=-_GHOST_CELLS_Y_,jCenterMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1;
  static const int kCenterMin=-_GHOST_CELLS_Z_,kCenterMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_-1;
#elif _MESH_DIMENSION_ == 2
  static const int iCenterMin=-_GHOST_CELLS_X_,iCenterMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;
  static const int jCenterMin=-_GHOST_CELLS_Y_,jCenterMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1;
  static const int kCenterMin=0,kCenterMax=0;
#elif _MESH_DIMENSION_ == 1
  static const int iCenterMin=-_GHOST_CELLS_X_,iCenterMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;
  static const int jCenterMin=0,jCenterMax=0;
  static const int kCenterMin=0,kCenterMax=0;
#endif

  cCenterNode *ptrCenterNode;

  for (k=kCenterMin;k<=kCenterMax;k++) for (j=jCenterMin;j<=jCenterMax;j++) for (i=iCenterMin;i<=iCenterMax;i++) {
    ptrCenterNode=startNode->block->GetCenterNode(getCenterNodeLocalNumber(i,j,k));

    if (ptrCenterNode!=NULL) if (ptrCenterNode->decrementConnectionCounter()==0) {
      ptrCenterNode->Measure=-1.0;
      CenterNodes.deleteElement(ptrCenterNode);
    }
  }
#endif

  blocks.deleteElement(startNode->block);
  startNode->block=NULL;

}



//===================================================================
bool deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode) {
  int i,j,k,iMax,jMax,kMax,iMin,jMin,kMin,nDownNode;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode,*neibNode; 


//######## DEBUG #############

  /*
if (startNode->Temp_ID==1169) {
*DiagnospticMessageStream << __LINE__ << __FILE__ << std::endl;
}
*/

//######## END DEBUG #########


  upNode=startNode->upNode;
  if (upNode==NULL) return false;

  //reset the mesh modified flag 
  nMeshModificationCounter++,meshModifiedFlag=true,meshModifiedFlag_CreateNewSpaceFillingCurve=true,meshModifiedFlag_CountMeshElements=true;

  //check if the neib nodes need to be removed first
  //1. check the face neibbours
  //2. check the neibours connected throught edges and corners

  //check the face neibours of all downNodes of hte 'upNode' 
  for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    downNode=upNode->downNode[nDownNode];

    /*
    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=downNode->neibNode[nface])!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      deleteTreeNode(neibNode->downNode[0]); 
    }
    */
    exit(__LINE__,__FILE__,"search is not implemented");

    if (downNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) deleteTreeNode(downNode->downNode[0]);
  

    //check the neibours connected throught edges and corners
    #if _MESH_DIMENSION_ == 1
    iMin=-1,iMax=1,jMin=0,jMax=0,kMin=0,kMax=0; 
    #elif _MESH_DIMENSION_ == 2
    iMin=-1,iMax=1,jMin=-1,jMax=1,kMin=0,kMax=0; 
    #elif _MESH_DIMENSION_ == 3 
    iMin=-1,iMax=1,jMin=-1,jMax=1,kMin=-1,kMax=1;
    #else 
    exit(__LINE__,__FILE__,"The wrong dimension");
    #endif

    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) if (abs(i)+abs(j)+abs(k)>1) {     
      neibNode=getNeibNode(i,j,k,downNode);

      if (neibNode!=NULL) if (neibNode->RefinmentLevel>startNode->RefinmentLevel) deleteTreeNode(neibNode);
    }
  }


  //create a new block for the upNode
  cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)];
  long int ii,iiMin,iiMax,jj,jjMin,jjMax,kk,kkMin,kkMax,nd;

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_OFF_ 
  upNode->block=blocks.newElement();
  #endif


  //transfer the corner nodes of the block to the upNode->block->cornerNodes[]
  #if _MESH_DIMENSION_ == 1 
  exit(__LINE__,__FILE__,"not implemented");
  #elif _MESH_DIMENSION_ == 2 
  iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkMin=0,kkMax=0;
  //collectNeibCornerNodes_2D_deleteTreeNode(startNode,newCornerNodeMap);
  exit(__LINE__,__FILE__,"function comented");
  #elif _MESH_DIMENSION_ == 3 
  iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkMin=-_GHOST_CELLS_Z_,kkMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
  collectNeibCornerNodes_3D_deleteTreeNode(startNode,newCornerNodeMap);
  #else 
  exit(__LINE__,__FILE__,"wrong dimension");
  #endif 

  for (kk=kkMin;kk<=kkMax;kk++) for (jj=jjMin;jj<=jjMax;jj++) for (ii=iiMin;ii<=iiMax;ii++) {  
    nd=getCornerNodeLocalNumber(ii,jj,kk);
    upNode->block->SetCornerNode(newCornerNodeMap[2*(ii-iiMin)][2*(jj-jjMin)][2*(kk-kkMin)],nd);
    if (upNode->block->GetCornerNode(nd)!=NULL) upNode->block->GetCornerNode(nd)->incrementConnectionCounter();

    //check if the node is defined
    if (upNode->block->GetCornerNode(nd)==NULL) {

      if ((ii>=0)&&(ii<=_BLOCK_CELLS_X_)&&(jj>=0)&&(jj<=_BLOCK_CELLS_Y_)&&(kk>=0)&&(kk<=_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
      else {
        int iNeib=0,jNeib=0,kNeib=0;

        if (ii<0) iNeib=-1;
        if (ii>_BLOCK_CELLS_X_) iNeib=1;

        if (jj<0) jNeib=-1;
        if (jj>_BLOCK_CELLS_Y_) jNeib=1;

        if (kk<0) kNeib=-1;
        if (kk>_BLOCK_CELLS_Z_) kNeib=1; 

        if ((neibNode=getNeibNode(iNeib,jNeib,kNeib,upNode))!=NULL) {
          double xMissing[3]={0.0,0.0,0.0};
          cTreeNodeAMR<cBlockAMR>* xMissingNode;

/*          *DiagnospticMessageStream << "$PREFIX:neibNode->Temp_ID=" << neibNode->Temp_ID << ", neibNode->RefinmentLevel=" << neibNode->RefinmentLevel << std::endl;

          //coordinates of the missing node
          *DiagnospticMessageStream << "$PREFIX:Coordinates of the missing node: x=" << (xMissing[0]=upNode->xmin[0]+(upNode->xmax[0]-upNode->xmin[0])/double(_BLOCK_CELLS_X_)*ii);
          if (_MESH_DIMENSION_>=2) *DiagnospticMessageStream << ", " << (xMissing[1]=upNode->xmin[1]+(upNode->xmax[1]-upNode->xmin[1])/double(_BLOCK_CELLS_Y_)*jj);
          if (_MESH_DIMENSION_==3) *DiagnospticMessageStream << ", " << (xMissing[2]=upNode->xmin[2]+(upNode->xmax[2]-upNode->xmin[2])/double(_BLOCK_CELLS_Z_)*kk);
          *DiagnospticMessageStream << std::endl;*/

          fprintf(DiagnospticMessageStream,"$PREFIX:neibNode->Temp_ID=%i, neibNode->RefinmentLevel=%i\n",neibNode->Temp_ID,neibNode->RefinmentLevel);

          //coordinates of the missing node
          fprintf(DiagnospticMessageStream,"$PREFIX:Coordinates of the missing node: x=%e",(xMissing[0]=upNode->xmin[0]+(upNode->xmax[0]-upNode->xmin[0])/double(_BLOCK_CELLS_X_)*ii));
          if (_MESH_DIMENSION_>=2) fprintf(DiagnospticMessageStream,", %e",(xMissing[1]=upNode->xmin[1]+(upNode->xmax[1]-upNode->xmin[1])/double(_BLOCK_CELLS_Y_)*jj));
          if (_MESH_DIMENSION_==3) fprintf(DiagnospticMessageStream,", %e",(xMissing[2]=upNode->xmin[2]+(upNode->xmax[2]-upNode->xmin[2])/double(_BLOCK_CELLS_Z_)*kk));
          fprintf(DiagnospticMessageStream,"\n");

          //find the node that contains xMissing and compare it with neibNode
          xMissingNode=findTreeNode(xMissing);

          if ((xMissing[0]<xMissingNode->xmin[0])||(xMissingNode->xmax[0]<xMissing[0])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");  
          if (_MESH_DIMENSION_>=2) if ((xMissing[1]<xMissingNode->xmin[1])||(xMissingNode->xmax[1]<xMissing[1])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          if (_MESH_DIMENSION_==3) if ((xMissing[2]<xMissingNode->xmin[2])||(xMissingNode->xmax[2]<xMissing[2])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");

          if (neibNode!=xMissingNode) { //recalculate the value of 'neibBlock' for debugging perposes
            cTreeNodeAMR<cBlockAMR>* recalculatedNeibNode;

            recalculatedNeibNode=getNeibNode(iNeib,jNeib,kNeib,upNode); 
            //exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          }
   
          //the position of the missed node in the 'xMissingNode' 
          long int ndMissingNode;
          int iMissingNode=0,jMissingNode=0,kMissingNode=0;

          ndMissingNode=findCornerNodeIndex(xMissing,iMissingNode,jMissingNode,kMissingNode,xMissingNode); 
//          *DiagnospticMessageStream << "$PREFIX:The (i,j,k) index of the missing node in 'xMissingNode': (i,j,k)=" << iMissingNode << ", " << jMissingNode << ", " << kMissingNode << ", ndMissingNode=" << ndMissingNode << std::endl;
          fprintf(DiagnospticMessageStream,"$PREFIX:The (i,j,k) index of the missing node in 'xMissingNode': (i,j,k)=%i, %i, %i, ndMissingNode=%i\n",iMissingNode,jMissingNode,kMissingNode,ndMissingNode);



          exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
        }
      }
    }
  }  

 
  //set up central nodes of the upNode->block
  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  cCenterNode *newCenterNode,*newCenterNodeMap[_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_][_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_][_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_];
  double xNode[3];
  int idim;


  //###########  DEBUG  #######
  /*
  if (startNode->Temp_ID==729) {
  	*DiagnospticMessageStream << __LINE__ << __FILE__ << std::endl;
  }
  */

  //########### END DEBUG ######

 
  #if _MESH_DIMENSION_ == 1
  exit(__LINE__,__FILE__,"not implemented");
  #elif _MESH_DIMENSION_ == 2
  exit(__LINE__,__FILE__,"not implemented");
  #elif _MESH_DIMENSION_ == 3 
  collectNeibCenterNodes_3D_deleteTreeNode(startNode,newCenterNodeMap); 
  #else 
  exit(__LINE__,__FILE__,"The value is not defined");
  #endif


  //create the missing nodes
  for (kk=kkMin;kk<kkMax;kk++) for (jj=jjMin;jj<jjMax;jj++) for (ii=iiMin;ii<iiMax;ii++) {
    nd=getCenterNodeLocalNumber(ii,jj,kk); 

    if (newCenterNodeMap[ii-iiMin][jj-jjMin][kk-kkMin]==NULL) {
      newCenterNode=CenterNodes.newElement();
      newCenterNodeMap[ii-iiMin][jj-jjMin][kk-kkMin]=newCenterNode;

//###########  DEBUG  #######
      /*
if (newCenterNode->Temp_ID==88861) {
	*DiagnospticMessageStream << __LINE__ << __FILE__ << std::endl;
}
*/

//########### END DEBUG ######

      xNode[0]=upNode->xmin[0]+(ii+0.5)*dxRootBlock[0]/(1<<upNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
      if (_MESH_DIMENSION_>=2) xNode[1]=upNode->xmin[1]+(jj+0.5)*dxRootBlock[1]/(1<<upNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_); 
      if (_MESH_DIMENSION_==3) xNode[2]=upNode->xmin[2]+(kk+0.5)*dxRootBlock[2]/(1<<upNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_); 

      #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
      newCenterNode->SetX(xNode);
      #endif

      for (idim=0;idim<_MESH_DIMENSION_;idim++) if ((xNode[idim]<xGlobalMin[idim]-EPS)||(xNode[idim]>xGlobalMax[idim]+EPS)) {
        if (newCenterNode->getNodeConnectionNumber()==0) CenterNodes.deleteElement(newCenterNode);
        newCenterNodeMap[i-iiMin][j-jjMin][k-kkMin]=NULL;
        break;
      }
    }
  }


  //transfer the node table 
  for (kk=kkMin;kk<kkMax;kk++) for (jj=jjMin;jj<jjMax;jj++) for (ii=iiMin;ii<iiMax;ii++) {
    nd=getCenterNodeLocalNumber(ii,jj,kk);
    upNode->block->SetCenterNode(newCenterNodeMap[ii-iiMin][jj-jjMin][kk-kkMin],nd);
    if (upNode->block->GetCenterNode(nd)!=NULL) upNode->block->GetCenterNode(nd)->incrementConnectionCounter();

    //check if the node is defined
    if (upNode->block->GetCenterNode(nd)==NULL) {
      if ((ii>=0)&&(ii<_BLOCK_CELLS_X_)&&(jj>=0)&&(jj<_BLOCK_CELLS_Y_)&&(kk>=0)&&(kk<_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
      else {
        int iNeib=0,jNeib=0,kNeib=0;

        if (ii<0) iNeib=-1;
        if (ii>=_BLOCK_CELLS_X_) iNeib=1;

        if (jj<0) jNeib=-1;
        if (jj>=_BLOCK_CELLS_Y_) jNeib=1;

        if (kk<0) kNeib=-1;
        if (kk>=_BLOCK_CELLS_Z_) kNeib=1;

        if ((neibNode=getNeibNode(iNeib,jNeib,kNeib,upNode))!=NULL) {
          double xMissing[3]={0.0,0.0,0.0};
          cTreeNodeAMR<cBlockAMR>* xMissingNode;

/*          *DiagnospticMessageStream << "$PREFIX:neibNode->Temp_ID=" << neibNode->Temp_ID << ", neibNode->RefinmentLevel=" << neibNode->RefinmentLevel << std::endl;

          //coordinates of the missing node
          *DiagnospticMessageStream << "$PREFIX:Coordinates of the missing node: x=" << (xMissing[0]=upNode->xmin[0]+(upNode->xmax[0]-upNode->xmin[0])/double(_BLOCK_CELLS_X_)*(ii+0.5));
          if (_MESH_DIMENSION_>=2) *DiagnospticMessageStream << ", " << (xMissing[1]=upNode->xmin[1]+(upNode->xmax[1]-upNode->xmin[1])/double(_BLOCK_CELLS_Y_)*(jj+0.5));
          if (_MESH_DIMENSION_==3) *DiagnospticMessageStream << ", " << (xMissing[2]=upNode->xmin[2]+(upNode->xmax[2]-upNode->xmin[2])/double(_BLOCK_CELLS_Z_)*(kk+0.5));
          *DiagnospticMessageStream << std::endl;*/

          fprintf(DiagnospticMessageStream,"$PREFIX:neibNode->Temp_ID=%i, neibNode->RefinmentLevel=%i\n",neibNode->Temp_ID,neibNode->RefinmentLevel);

          //coordinates of the missing node
          fprintf(DiagnospticMessageStream,"$PREFIX:Coordinates of the missing node: x=%e",(xMissing[0]=upNode->xmin[0]+(upNode->xmax[0]-upNode->xmin[0])/double(_BLOCK_CELLS_X_)*(ii+0.5)));
          if (_MESH_DIMENSION_>=2) fprintf(DiagnospticMessageStream,", %e",(xMissing[1]=upNode->xmin[1]+(upNode->xmax[1]-upNode->xmin[1])/double(_BLOCK_CELLS_Y_)*(jj+0.5)));
          if (_MESH_DIMENSION_==3) fprintf(DiagnospticMessageStream,", %e",(xMissing[2]=upNode->xmin[2]+(upNode->xmax[2]-upNode->xmin[2])/double(_BLOCK_CELLS_Z_)*(kk+0.5)));
          fprintf(DiagnospticMessageStream,"\n");

          //find the node that contains xMissing and compare it with neibNode
          xMissingNode=findTreeNode(xMissing);

          if ((xMissing[0]<xMissingNode->xmin[0])||(xMissingNode->xmax[0]<xMissing[0])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          if (_MESH_DIMENSION_>=2) if ((xMissing[1]<xMissingNode->xmin[1])||(xMissingNode->xmax[1]<xMissing[1])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          if (_MESH_DIMENSION_==3) if ((xMissing[2]<xMissingNode->xmin[2])||(xMissingNode->xmax[2]<xMissing[2])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");

          if (neibNode!=xMissingNode) { //recalculate the value of 'neibBlock' for debugging perposes
            cTreeNodeAMR<cBlockAMR>* recalculatedNeibNode;

            recalculatedNeibNode=getNeibNode(iNeib,jNeib,kNeib,upNode);
            //exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          }
  
          //the position of the missed node in the 'xMissingNode'
          long int ndMissingNode;
          int iMissingNode=0,jMissingNode=0,kMissingNode=0;

          ndMissingNode=findCenterNodeIndex(xMissing,iMissingNode,jMissingNode,kMissingNode,xMissingNode);
//          *DiagnospticMessageStream << "$PREFIX:The (i,j,k) index of the missing node in 'xMissingNode': (i,j,k)=" << iMissingNode << ", " << jMissingNode << ", " << kMissingNode << ", ndMissingNode=" << ndMissingNode << std::endl;
          fprintf(DiagnospticMessageStream,"$PREFIX:The (i,j,k) index of the missing node in 'xMissingNode': (i,j,k)=%i, %i, %i, ndMissingNode=%i\n",iMissingNode,jMissingNode,kMissingNode,ndMissingNode);


          exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
        }
      }
    }
  }
  #endif 


  //delete the tree node
  for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    downNode=upNode->downNode[nDownNode];

    //remove the corner nodes
    for (kk=kkMin;kk<=kkMax;kk++) for (jj=jjMin;jj<=jjMax;jj++) for (ii=iiMin;ii<=iiMax;ii++) {
      nd=getCornerNodeLocalNumber(ii,jj,kk);
      if (downNode->block->GetCornerNode(nd)->decrementConnectionCounter()==0) CornerNodes.deleteElement(downNode->block->GetCornerNode(nd));
    }

    //remove the center nodes
    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    for (kk=kkMin;kk<kkMax;kk++) for (jj=jjMin;jj<jjMax;jj++) for (ii=iiMin;ii<iiMax;ii++) {
      nd=getCenterNodeLocalNumber(ii,jj,kk);
      if (downNode->block->GetCenterNode(nd)->decrementConnectionCounter()==0) CenterNodes.deleteElement(downNode->block->GetCenterNode(nd));
    }
    #endif


    //remove the connections with the neighbors
    /*
    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=downNode->neibNode[nface])!=NULL) {
      for (int nNeibFace=0;nNeibFace<2*_MESH_DIMENSION_;nNeibFace++) if (neibNode->neibNode[nNeibFace]==downNode) {
        neibNode->neibNode[nNeibFace]=NULL;
        break;
      }

      downNode->neibNode[nface]=NULL;  
    }
*/
    exit(__LINE__,__FILE__,"removing of the connections is not implemented");

    //remove the connection with the parent tree node
    upNode->downNode[nDownNode]=NULL;

    //empty the memory occupied by the block
    if (downNode->block!=NULL) {
      //remove the nodes of the block
      exit(__LINE__,__FILE__,"not implemnted");

      //remove the block
      blocks.deleteElement(downNode->block);
      downNode->block=NULL;
    }

    //remove the list of the descriptors to the internal surfaces installed into the mesh
    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    cInternalBoundaryConditionsDescriptor *nextDescriptor;

    while (downNode->InternalBoundaryDescriptorList!=NULL) {
      nextDescriptor=downNode->InternalBoundaryDescriptorList->nextInternalBCelement;
      InternalBoundaryDescriptors.deleteElement(downNode->InternalBoundaryDescriptorList);
      downNode->InternalBoundaryDescriptorList=nextDescriptor;
    }

    downNode->InternalBoundaryDescriptorList=NULL;

    cTreeNodeAMR<cBlockAMR> *next,*prev;
    next=downNode->DomainSurfaceBoundaryList_Next,prev=downNode->DomainSurfaceBoundaryList_Prev;
    if (next!=NULL) next->DomainSurfaceBoundaryList_Prev=prev;
    if (prev!=NULL) prev->DomainSurfaceBoundaryList_Next=next;

    if (DomainSurfaceBoundaryList==downNode) DomainSurfaceBoundaryList=next;

    downNode->DomainSurfaceBoundaryList_Next=NULL,downNode->DomainSurfaceBoundaryList_Prev=NULL;
    #endif

    //remove the node
    treeNodes.deleteElement(downNode);

  }

  return true; 
}


    
//=================================================================================================

/*
 bool splitTreeNode(cTreeNodeAMR<cBlockAMR> *startNode) {
   int idim,i,j,k;
   long int nd,rLevel,nDownNode;
   double dxBlock[3];


if (startNode->Temp_ID==44) {
*DiagnospticMessageStream << __LINE__ << std::endl;
}



   if (startNode->RefinmentLevel>=_MAX_REFINMENT_LEVEL_) return false;
   nMeshModificationCounter++,meshModifiedFlag=true,meshModifiedFlag_CreateNewSpaceFillingCurve=true,meshModifiedFlag_CountMeshElements=true;

   //create the daugher blocks
   if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //there is no daugher blocks
     cTreeNodeAMR<cBlockAMR> *newTreeNode,*neibNode;
     cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)];
     cCornerNode *newCornerNode;

     //check if the neibNodes need to be refined  
     for (i=-1;i<2;i++) for (j=-1;j<2;j++) for (k=-1;k<2;k++) if ((neibNode=getNeibNode(i,j,k,startNode))!=NULL) if (neibNode->RefinmentLevel<startNode->RefinmentLevel) { //refine the neibNode   
       splitTreeNode(neibNode);
     }

     #if _MESH_DIMENSION_ == 1
     collectNeibNodeMap_1D(startNode,newCornerNodeMap);
     #elif _MESH_DIMENSION_ == 2
     collectNeibNodeMap_2D(startNode,newCornerNodeMap); 
     #else 
     collectNeibCornerNodeMap_3D(startNode,newCornerNodeMap);
     #endif 

     //transfer the nodes from 'startNode' to 'newCornerNodeMap' and create the mission nodes 
     #if _MESH_DIMENSION_ == 1
     static const int iCellMax=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jCellMax=0,kCellMax=0;
     static const int iCellMin=-_GHOST_CELLS_X_,jCellMin=0,kCellMin=0;
     #elif _MESH_DIMENSION_ == 2
     static const int iCellMax=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jCellMax=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kCellMax=0;
     static const int iCellMin=-_GHOST_CELLS_X_,jCellMin=-_GHOST_CELLS_Y_,kCellMin=0;
     #else
     static const int iCellMax=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jCellMax=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kCellMax=2*_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
     static const int iCellMin=-_GHOST_CELLS_X_,jCellMin=-_GHOST_CELLS_Y_,kCellMin=-_GHOST_CELLS_Z_; 
     #endif

      rLevel=startNode->RefinmentLevel+1;
      dxBlock[0]=(xGlobalMax[0]-xGlobalMin[0])/(1<<rLevel)/_BLOCK_CELLS_X_;
      if (_MESH_DIMENSION_>1) dxBlock[1]=(xGlobalMax[1]-xGlobalMin[1])/(1<<rLevel)/_BLOCK_CELLS_Y_;
      if (_MESH_DIMENSION_>2) dxBlock[2]=(xGlobalMax[2]-xGlobalMin[2])/(1<<rLevel)/_BLOCK_CELLS_Z_;

      for (k=kCellMin;k<=kCellMax;k++) for (j=jCellMin;j<=jCellMax;j++) for (i=iCellMin;i<=iCellMax;i++) if (newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]==NULL) {
        //check if the new node has to be created or transfered from the startBlock
        if ((2*(i/2)!=i)||(2*(j/2)!=j)||(2*(k/2)!=k)) { // create new nodes
          if (newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]==NULL) {
            newCornerNode=CornerNodes.newElement();
            newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=newCornerNode;  

            double xNode[3];

            xNode[0]=startNode->xmin[0]+i*dxBlock[0];
            if (_MESH_DIMENSION_>=2) xNode[1]=startNode->xmin[1]+j*dxBlock[1];
            if (_MESH_DIMENSION_==3) xNode[2]=startNode->xmin[2]+k*dxBlock[2];

            newCornerNode->SetX(xNode);

            for (idim=0;idim<_MESH_DIMENSION_;idim++) if ((xNode[idim]<xGlobalMin[idim]-0.001*dxBlock[idim])||(xNode[idim]>xGlobalMax[idim]+0.001*dxBlock[idim])) {
              if (newCornerNode->getNodeConnectionNumber()==0) CornerNodes.deleteElement(newCornerNode);
              newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=NULL; 
              break;
            }
          }
        }
        else {
          nd=getCornerNodeLocalNumber(i/2,j/2,k/2);
          newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=startNode->block->GetCornerNode(nd);
        }
      }

      //collect and create the block's center nodes
      #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_ 
      cCenterNode *newCenterNodeMap[2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)];
      cCenterNode *newCenterNode;

      #if _MESH_DIMENSION_ == 1
      exit(__LINE__,__FILE__,"not implemented"); 
      #elif _MESH_DIMENSION_ == 2
      exit(__LINE__,__FILE__,"not implemented");
      #else
      collectNeibCenterNodeMap_3D(startNode,newCenterNodeMap);
      #endif 

      //create new 'center nodes' when its needed
      double xNode[3];

      for (k=kCellMin;k<kCellMax;k++) for (j=jCellMin;j<jCellMax;j++) for (i=iCellMin;i<iCellMax;i++) if (newCenterNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]==NULL) { //create bew new node 
        newCenterNode=CenterNodes.newElement();
        newCenterNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=newCenterNode;


//################  DEBUG ##############

if (newCenterNode->Temp_ID==88861) {
	*DiagnospticMessageStream << __LINE__ << std::endl;
}

//################ END DEBUG ###########


        xNode[0]=startNode->xmin[0]+(i+0.5)*dxBlock[0];
        if (_MESH_DIMENSION_>=2) xNode[1]=startNode->xmin[1]+(j+0.5)*dxBlock[1];
        if (_MESH_DIMENSION_==3) xNode[2]=startNode->xmin[2]+(k+0.5)*dxBlock[2];

        #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
        newCenterNode->SetX(xNode);
        #endif

        for (idim=0;idim<_MESH_DIMENSION_;idim++) if ((xNode[idim]<xGlobalMin[idim]-0.001*dxBlock[idim])||(xNode[idim]>xGlobalMax[idim]+0.001*dxBlock[idim])) {
          if (newCenterNode->getNodeConnectionNumber()==0) CenterNodes.deleteElement(newCenterNode);
          newCenterNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=NULL;
          break;
        }
      } 
      #endif


      //create the downNodes
      int iOffset,jOffset,kOffset,kDownNode,jDownNode,iDownNode;

      #if _MESH_DIMENSION_ == 1
      static const int iDownNodeMax=2,jDownNodeMax=1,kDownNodeMax=1;
      #elif _MESH_DIMENSION_ == 2
      static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=1;
      #else
      static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=2;
      static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kMin=-_GHOST_CELLS_Z_,kMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
      #endif


      for (kDownNode=0;kDownNode<kDownNodeMax;kDownNode++) for (jDownNode=0;jDownNode<jDownNodeMax;jDownNode++) for (iDownNode=0;iDownNode<iDownNodeMax;iDownNode++) if (startNode->GetDownNode(iDownNode,jDownNode,kDownNode)==NULL) {
        newTreeNode=treeNodes.newElement();
        startNode->SetDownNode(newTreeNode,iDownNode,jDownNode,kDownNode);

        //set block into the tree
        newTreeNode->block=blocks.newElement();
        newTreeNode->block->SetRefinmentLevel(rLevel);
        newTreeNode->RefinmentLevel=rLevel; 

        //upBlock
        newTreeNode->upNode=startNode;

        //init the corner nodes of the new blocks
        iOffset=iDownNode*_BLOCK_CELLS_X_,jOffset=jDownNode*_BLOCK_CELLS_Y_,kOffset=kDownNode*_BLOCK_CELLS_Z_;

        for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) { 
          nd=getCornerNodeLocalNumber(i,j,k); 

          newTreeNode->block->SetCornerNode(newCornerNodeMap[iOffset+i-iMin][j+jOffset-jMin][k+kOffset-kMin],nd);
          if (newTreeNode->block->GetCornerNode(nd)!=NULL) newTreeNode->block->GetCornerNode(nd)->incrementConnectionCounter();
        }

        //init the center nodes of the block 
        #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
        for (k=kMin;k<kMax;k++) for (j=jMin;j<jMax;j++) for (i=iMin;i<iMax;i++) {
          nd=getCenterNodeLocalNumber(i,j,k);

          newTreeNode->block->SetCenterNode(newCenterNodeMap[iOffset+i-iMin][j+jOffset-jMin][k+kOffset-kMin],nd);
          if (newTreeNode->block->GetCenterNode(nd)!=NULL) newTreeNode->block->GetCenterNode(nd)->incrementConnectionCounter();
        }
        #endif 
  


        //init newTreeNode->xmin
        double *x,*xmin,*xmax;

        nd=getCornerNodeLocalNumber(0,0,0);
        xmin=newTreeNode->xmin; 
        x=newTreeNode->block->GetCornerNode(nd)->GetX();
        for (idim=0;idim<_MESH_DIMENSION_;idim++) xmin[idim]=x[idim]; 

        nd=getCornerNodeLocalNumber(_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_);
        xmax=newTreeNode->xmax;
        x=newTreeNode->block->GetCornerNode(nd)->GetX();
        for (idim=0;idim<_MESH_DIMENSION_;idim++) xmax[idim]=x[idim];

        //determine if the newTreeNode is intersected by any of the internal surface installed into the mesh
        #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
        cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptor,*newDescriptor;
        int IntersectionCode=-1;

        for (InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;InternalBoundaryDescriptor!=NULL;InternalBoundaryDescriptor=InternalBoundaryDescriptor->nextInternalBCelement) {
          switch(InternalBoundaryDescriptor->BondaryType) {
          case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
            IntersectionCode=((cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement))->BlockIntersection(xmin,xmax,EPS);
            break;
          default:
            exit(__LINE__,__FILE__,"Error: The internal boundary type is not recognized");
          }

          switch (IntersectionCode) {
          case _AMR_BLOCK_INSIDE_DOMAIN_:
            //do nothing
            break;
          case _AMR_BLOCK_OUTSIDE_DOMAIN_:
            //the block (node) that is outside of the computational domain can be removed from the tree
            #if _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_ == _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_REMOVE_
            exit(__LINE__,__FILE__,"not implemented");
            #elif _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_ ==  _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_KEEP_
            //do nothing
            #else
            exit(__LINE__,__FILE__,"Error: the code is not recognized");
            #endif

            break;
          case _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_:
            //add the descriptor to the newTreeNode
            newDescriptor=InternalBoundaryDescriptors.newElement();

            //copy the content of the descriptor
            newDescriptor->BondaryType=InternalBoundaryDescriptor->BondaryType;
            newDescriptor->BoundaryElement=InternalBoundaryDescriptor->BoundaryElement;

            //add the descriptor to newTreeNode
            newDescriptor->nextInternalBCelement=newTreeNode->InternalBoundaryDescriptorList;
            newTreeNode->InternalBoundaryDescriptorList=newDescriptor;

            break;
          default:
            exit(__LINE__,__FILE__,"Error: The intersection code is not recognized");
          }



        }

        #endif
      }



      //connect the newly created elements
      #if _MESH_DIMENSION_ == 1
      startNode->downNode[0]->neibNode[1]=startNode->downNode[1];
      startNode->downNode[1]->neibNode[0]=startNode->downNode[0];
      #else

      #if _MESH_DIMENSION_ == 2

      //conenction table
      static const int innerNodeConnectionMap[4][4]={ {-1,-1,1,2}, {0,-1,-1,3}, {-1,0,3,-1}, {2,1,-1,-1}};

      #elif _MESH_DIMENSION_ == 3

      static const int innerNodeConnectionMap[8][6]={ {-1,1,-1,2,-1,4}, {0,-1,-1,3,-1,5}, {-1,3,0,-1,-1,6}, {2,-1,1,-1,-1,7},
                                                      {-1,5,-1,6,0,-1}, {4,-1,-1,7,1,-1}, {-1,7,4,-1,2,-1,}, {6,-1,5,-1,3,-1}   };

      #endif

      for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) for (i=0;i<2*_MESH_DIMENSION_;i++) {
        startNode->downNode[nDownNode]->neibNode[i]=(innerNodeConnectionMap[nDownNode][i]!=-1) ? startNode->downNode[innerNodeConnectionMap[nDownNode][i]] : NULL;
      }
      #endif

      reconnectDownTreeNode(startNode);
   }

   return true;
}


*/

bool splitTreeNode(cTreeNodeAMR<cBlockAMR> *startNode) {
  int i,j,k;


/*
if (startNode->Temp_ID==15) {
*DiagnospticMessageStream << __LINE__ << std::endl;
}
*/



  if ((startNode->RefinmentLevel>=_MAX_REFINMENT_LEVEL_)||(startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_)) return false;
  nMeshModificationCounter++,meshModifiedFlag=true,meshModifiedFlag_CreateNewSpaceFillingCurve=true,meshModifiedFlag_CountMeshElements=true;

  //create the daugher blocks
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //there is no daugher blocks
    cTreeNodeAMR<cBlockAMR> *newTreeNode,*neibNode;

    //check if the neibNodes need to be refined
    for (i=-1;i<2;i++) for (j=-1;j<2;j++) for (k=-1;k<2;k++) if ((neibNode=getNeibNode(i,j,k,startNode))!=NULL) if (neibNode->RefinmentLevel<startNode->RefinmentLevel) { //refine the neibNode
      if (neibNode->RefinmentLevel<startNode->RefinmentLevel-1) exit(__LINE__,__FILE__,"Error: the mesh is not consistent");

      splitTreeNode(neibNode);



//################## DEBUG #####################

//checkMeshConsistency(rootTree);
//#################  END DEBUG #################


    }



     //create the downNodes
     int kDownNode,jDownNode,iDownNode;

     #if _MESH_DIMENSION_ == 1
     static const int iDownNodeMax=2,jDownNodeMax=1,kDownNodeMax=1;
     #elif _MESH_DIMENSION_ == 2
     static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=1;
     #else
     static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=2;
     #endif


     double *xmin,*xmax,dxMddle[3]={0.0,0.0,0.0};

     dxMddle[0]=(startNode->xmax[0]+startNode->xmin[0])/2.0;
     if (_MESH_DIMENSION_>1) dxMddle[1]=(startNode->xmax[1]+startNode->xmin[1])/2.0;
     if (_MESH_DIMENSION_>2) dxMddle[2]=(startNode->xmax[2]+startNode->xmin[2])/2.0;

     for (kDownNode=0;kDownNode<kDownNodeMax;kDownNode++) for (jDownNode=0;jDownNode<jDownNodeMax;jDownNode++) for (iDownNode=0;iDownNode<iDownNodeMax;iDownNode++) if (startNode->GetDownNode(iDownNode,jDownNode,kDownNode)==NULL) {
       newTreeNode=treeNodes.newElement();
       startNode->SetDownNode(newTreeNode,iDownNode,jDownNode,kDownNode);
       newTreeNode->RefinmentLevel=startNode->RefinmentLevel+1;

       //upBlock
       newTreeNode->upNode=startNode;

       //init the node ID
       GetAMRnodeID(newTreeNode->AMRnodeID,newTreeNode);

       //init newTreeNode->xmin,xmax 
       xmin=newTreeNode->xmin;
       xmax=newTreeNode->xmax;

       xmin[0]=(iDownNode==0) ? startNode->xmin[0] : dxMddle[0];
       xmax[0]=(iDownNode==0) ? dxMddle[0] : startNode->xmax[0];

       if (_MESH_DIMENSION_>1) {
         xmin[1]=(jDownNode==0) ? startNode->xmin[1] : dxMddle[1];
         xmax[1]=(jDownNode==0) ? dxMddle[1] : startNode->xmax[1];
       } 

       if (_MESH_DIMENSION_>2) {
         xmin[2]=(kDownNode==0) ? startNode->xmin[2] : dxMddle[2];
         xmax[2]=(kDownNode==0) ? dxMddle[2] : startNode->xmax[2];
       }
      
       //allocate the blocks:THE ORDER IS IMPORTANT
       AllocateBlock(newTreeNode);

       //determine if the newTreeNode is intersected by any of the internal surface installed into the mesh
       #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
       cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptor,*newDescriptor;
       int IntersectionCode=-1;
       bool SurfaceBoundaryListAdded=false;

       for (InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;InternalBoundaryDescriptor!=NULL;InternalBoundaryDescriptor=InternalBoundaryDescriptor->nextInternalBCelement) {
         switch(InternalBoundaryDescriptor->BondaryType) {
         case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
           IntersectionCode=((cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement))->BlockIntersection(xmin,xmax,EPS);
           break;
         case _INTERNAL_BOUNDARY_TYPE_CIRCLE_:
           IntersectionCode=((cInternalCircleData*)(InternalBoundaryDescriptor->BoundaryElement))->BlockIntersection(xmin,xmax,EPS);
           break;
         case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:
           IntersectionCode=((cInternalSphere1DData*)(InternalBoundaryDescriptor->BoundaryElement))->BlockIntersection(xmin,xmax,EPS);
           break;
         case _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_:
           IntersectionCode=((cInternalRotationBodyData*)(InternalBoundaryDescriptor->BoundaryElement))->BlockIntersection(xmin,xmax,EPS);
           break;
         case _INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_:
           //distribute the boundary faces
           DistributeBoundaryCutBlocks(newTreeNode,startNode->FirstTriangleCutFace);

           //!!!! nastran surface faces are installed into the mesh by function DistributeBoundaryCutBlocks(newTreeNode,startNode->FirstTriangleCutFace) above!!!
           IntersectionCode=(newTreeNode->FirstTriangleCutFace==NULL) ? _AMR_BLOCK_INSIDE_DOMAIN_ : _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
           break;
         default:
           exit(__LINE__,__FILE__,"Error: The internal boundary type is not recognized");
         }

         switch (IntersectionCode) {
         case _AMR_BLOCK_INSIDE_DOMAIN_:
           //do nothing
           break;
         case _AMR_BLOCK_OUTSIDE_DOMAIN_:
           //the block (node) that is outside of the computational domain can be removed from the tree
           #if _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_ == _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_REMOVE_
           exit(__LINE__,__FILE__,"not implemented");
           #elif _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_ ==  _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_KEEP_
           //do nothing
           #else
           exit(__LINE__,__FILE__,"Error: the code is not recognized");
           #endif

           break;
         case _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_:
           //add the descriptor to the newTreeNode
           newDescriptor=InternalBoundaryDescriptors.newElement();

           //copy the content of the descriptor
           newDescriptor->BondaryType=InternalBoundaryDescriptor->BondaryType;
           newDescriptor->BoundaryElement=InternalBoundaryDescriptor->BoundaryElement;

           //add the descriptor to newTreeNode
           newDescriptor->nextInternalBCelement=newTreeNode->InternalBoundaryDescriptorList;
           newTreeNode->InternalBoundaryDescriptorList=newDescriptor;

           //add newTreeNode to the lsit of nodes that intersect the surface boundary of the computational domain
           if (SurfaceBoundaryListAdded==false) {
             newTreeNode->DomainSurfaceBoundaryList_Prev=NULL,newTreeNode->DomainSurfaceBoundaryList_Next=DomainSurfaceBoundaryList;
             DomainSurfaceBoundaryList=newTreeNode;
           }


           break;
         default:
           exit(__LINE__,__FILE__,"Error: The intersection code is not recognized");
         }



       }

       #endif
     }



     //connect the newly created elements

     /*
     #if _MESH_DIMENSION_ == 1
     startNode->downNode[0]->neibNode[1]=startNode->downNode[1];
     startNode->downNode[1]->neibNode[0]=startNode->downNode[0];
     #else

     #if _MESH_DIMENSION_ == 2

     //conenction table
     static const int innerNodeConnectionMap[4][4]={ {-1,1,-1,2}, {0,-1,-1,3}, {-1,3,0,-1}, {2,-1,1,-1}};

     #elif _MESH_DIMENSION_ == 3

     static const int innerNodeConnectionMap[8][6]={ {-1,1,-1,2,-1,4}, {0,-1,-1,3,-1,5}, {-1,3,0,-1,-1,6}, {2,-1,1,-1,-1,7},
                                                     {-1,5,-1,6,0,-1}, {4,-1,-1,7,1,-1}, {-1,7,4,-1,2,-1,}, {6,-1,5,-1,3,-1}   };

     #endif

     for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) for (i=0;i<2*_MESH_DIMENSION_;i++) {
       startNode->downNode[nDownNode]->neibNode[i]=(innerNodeConnectionMap[nDownNode][i]!=-1) ? startNode->downNode[innerNodeConnectionMap[nDownNode][i]] : NULL;
     }
     #endif

     reconnectDownTreeNode(startNode);
     */

     //connect
     // 1. connect the newly created nodes (startNode->downNode[:])
     // 2. connect the created nodes with the surrounding blocks (startNode->downNode[:] connect to startNode->neibNode[:])
     // 3. set the connections of the startNode to default value (startNode->neibNode[:]=NULL)



#if _MESH_DIMENSION_ == 1
     int nface,nUpFace,nDownNode;

     //connect internal nodes
     startNode->downNode[0]->SetNeibFace(startNode->downNode[1],1,0,0);
     startNode->downNode[1]->SetNeibFace(startNode->downNode[0],0,0,0);

     //connect external nodes
     struct cExternalNodesFace {
       int nUpFace;
       int nNeibFace;
     };

     static const cExternalNodesFace ExternalNodesFaceConnectionMap_upNode[2][2]={
       {{0,1},{-1,-1}},{{-1,-1},{1,0}} };

     for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
       //make connection through the faces
       for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((nUpFace=ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].nUpFace)!=-1) {
         neibNode=startNode->GetNeibFace(nUpFace,0,0);

         startNode->downNode[nDownNode]->SetNeibFace(neibNode,nface,0,0);

         //connect the neibNode to startNode->downNode
         if (neibNode!=NULL) neibNode->SetNeibFace(startNode->downNode[nDownNode],ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].nNeibFace,0,0);
       }
     }


     //check consistency of the node's reconnection
#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
#if _CHECK_MESH_CONSISTANCY_ == _ON_AMR_MESH_

     for (register int nface=0;nface<2*_MESH_DIMENSION_;nface++) {
       neibNode=startNode->GetNeibFace(nface,0,0);
       if (neibNode!=NULL) if (neibNode->CheckNeibNode(startNode)==true) exit(__LINE__,__FILE__,"Error: the mesh is not consistent");
     }

#endif
#endif

#elif _MESH_DIMENSION_ == 2
     int nface,nNeibNode,nUpCornerNode,nUpFace,iFace,nd,nDownNode;
//     cTreeNodeAMR<cBlockAMR> *neibNode;


     struct cExternalNodesFace {
       int nUpFace;
       int iFace;
       int nNeibFace;
     };

     struct cExternalNodesCorner_upNode {
       int nUpCornerNode;
       int nNeibCornerNode;
     };

     struct cExternalNodesCorner_downNode {
       int nUpFace;
       int iFace;
       int nNeibCornerNode;
     };

     static const int InternalNodesFaceConnectionMap[4][4]={ {-1,1,-1,2},  {0,-1,-1,3},  {-1,3,0,-1},  {2,-1,1,-1}};
     static const int InternalNodesNodeConnectionMap[4][4]={ {-1,-1,-1,3}, {-1,-1,2,-1}, {-1,1,-1,-1}, {0,-1,-1,-1}};


     static const cExternalNodesFace ExternalNodesFaceConnectionMap_upNode[4][4]={
         {{0,0,1},{-1,-1,-1},{2,0,3},{-1,-1,-1}},{{-1,-1,-1},{1,0,0},{2,1,3},{-1,-1,-1}},
         {{0,1,1},{-1,-1,-1},{-1,-1,-1},{3,0,2}},{{-1,-1,-1},{1,1,0},{-1,-1,-1},{3,1,2}}};

     static const cExternalNodesCorner_upNode ExternalNodesCornerNodeConnectionMap_upNode[4][4]={
         {{0,3},{-1,-1},{-1,-1},{-1,-1}},
         {{-1,-1},{1,2},{-1,-1},{-1,-1}},
         {{-1,-1},{-1,-1},{2,1},{-1,-1}},
         {{-1,-1},{-1,-1},{-1,-1},{3,0}} };

     static const cExternalNodesCorner_downNode ExternalNodesCornerNodeConnectionMap_downNode[4][4]={
         {{-1,-1,-1},{2,1,2},{0,1,1},{-1,-1,-1}},
         {{2,0,3},{-1,-1,-1},{-1,-1,-1},{1,1,0}},
         {{0,0,3},{-1,-1,-1},{-1,-1,-1},{3,1,0}},
         {{-1,-1,-1},{1,0,2},{3,0,1},{-1,-1,-1}} };


     //connect internal nodes
     for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
       for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((nNeibNode=InternalNodesFaceConnectionMap[nDownNode][nface])!=-1) {
         for (i=0;i<2;i++) startNode->downNode[nDownNode]->SetNeibFace(startNode->downNode[nNeibNode],nface,i,0);
       }

       for (nd=0;nd<4;nd++) if ((nNeibNode=InternalNodesNodeConnectionMap[nDownNode][nd])!=-1) {
         startNode->downNode[nDownNode]->SetNeibCorner(startNode->downNode[nNeibNode],nd);
       }
     }

     //connect with the "external" neighbors
     for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
       //make connection through the faces
       for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((nUpFace=ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].nUpFace)!=-1) {
         iFace=ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].iFace;
         neibNode=startNode->GetNeibFace(nUpFace,iFace,0);

         for (i=0;i<2;i++) startNode->downNode[nDownNode]->SetNeibFace(neibNode,nface,i,0);

         //connect the neibNode to startNode->downNode
         if (neibNode!=NULL) {
           if (neibNode->RefinmentLevel==startNode->RefinmentLevel) {
             neibNode->SetNeibFace(startNode->downNode[nDownNode],ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].nNeibFace,ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].iFace,0);
           }
           else if (neibNode->RefinmentLevel==startNode->RefinmentLevel+1) {
             for (i=0;i<2;i++) neibNode->SetNeibFace(startNode->downNode[nDownNode],ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].nNeibFace,i,0);
           }
           else exit(__LINE__,__FILE__,"Error: the mesh is not consistent");
         }
       }

       //make connection through the corner nodes of startNode
       for (nd=0;nd<4;nd++) if ((nUpCornerNode=ExternalNodesCornerNodeConnectionMap_upNode[nDownNode][nd].nUpCornerNode)!=-1) {
         neibNode=startNode->GetNeibCorner(nUpCornerNode);

         startNode->downNode[nDownNode]->SetNeibCorner(neibNode,nUpCornerNode);
         if (neibNode!=NULL) neibNode->SetNeibCorner(startNode->downNode[nDownNode],ExternalNodesCornerNodeConnectionMap_upNode[nDownNode][nd].nNeibCornerNode);
       }

       //make connection through the corner nodes of startNode->downNode[:]
       for (nd=0;nd<4;nd++) if ((nUpFace=ExternalNodesCornerNodeConnectionMap_downNode[nDownNode][nd].nUpFace)!=-1) {
         if (startNode->GetNeibFace(nUpFace,0,0)!=startNode->GetNeibFace(nUpFace,1,0)) {
           neibNode=startNode->GetNeibFace(nUpFace,ExternalNodesCornerNodeConnectionMap_downNode[nDownNode][nd].iFace,0);

           startNode->downNode[nDownNode]->SetNeibCorner(neibNode,nd);
           if (neibNode!=NULL) neibNode->SetNeibCorner(startNode->downNode[nDownNode],ExternalNodesCornerNodeConnectionMap_downNode[nDownNode][nd].nNeibCornerNode);
         }
       }
     }


     //check consistency of the node's reconnection
#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
#if _CHECK_MESH_CONSISTANCY_ == _ON_AMR_MESH_

     exit(__LINE__,__FILE__,"not implemented");
#endif
#endif
     //clean connections of startNode
     for (nd=0;nd<(1<<_MESH_DIMENSION_);nd++) startNode->SetNeibCorner(NULL,nd);
     for (nface=0;nface<2*_MESH_DIMENSION_;nface++) for (i=0;i<2;i++) startNode->SetNeibFace(NULL,nface,i,0);

#elif _MESH_DIMENSION_ == 3

     register int nDownNode,nList;
     int nd;
     cTreeNodeAMR<cBlockAMR>* t,*downNode;




//===========================   DEBUG =========================
/*
     if (startNode->Temp_ID==3) {
       *DiagnospticMessageStream << __LINE__ << std::endl;
     }
*/
//=========================== END DEBUG =======================



     //INTERNAL CONENCTIONS: connection of the downNodes



     //connection through the new faces created within the 'startNode'
     static const int InternalNodesFaceConnectionMap[8][6]={
         {-1,1,-1,2,-1,4}, {0,-1,-1,3,-1,5}, {-1,3,0,-1,-1,6}, {2,-1,1,-1,-1,7},
         {-1,5,-1,6,0,-1}, {4,-1,-1,7,1,-1}, {-1,7,4,-1,2,-1}, {6,-1,5,-1,3,-1} };


     for (nDownNode=0;nDownNode<8;nDownNode++) {
       register int nface,iFace,jFace;

       downNode=startNode->downNode[nDownNode];

       for (nface=0;nface<6;nface++) if (InternalNodesFaceConnectionMap[nDownNode][nface]!=-1) {
         t=startNode->downNode[InternalNodesFaceConnectionMap[nDownNode][nface]];

         for (iFace=0;iFace<2;iFace++) for (jFace=0;jFace<2;jFace++) downNode->SetNeibFace(t,nface,iFace,jFace);
       }
     }



     //connection through the central node created in the meddle of 'startNode'
     static const int InternalNodesNodeConnectionMap[8][8]={ {-1,-1,-1,-1,-1,-1,-1,7}, {-1,-1,-1,-1,-1,-1,6,-1}, {-1,-1,-1,-1,-1,5,-1,-1}, {-1,-1,-1,-1,4,-1,-1,-1},
                                                             {-1,-1,-1,3,-1,-1,-1,-1}, {-1,-1,2,-1,-1,-1,-1,-1}, {-1,1,-1,-1,-1,-1,-1,-1}, {0,-1,-1,-1,-1,-1,-1,-1} };

     for (nDownNode=0,nd=7;nDownNode<8;nDownNode++,nd--) {
       startNode->downNode[nDownNode]->SetNeibCorner(startNode->downNode[InternalNodesNodeConnectionMap[nDownNode][nd]],nd);
     }



     //connection through the new edges created within the node 'startNode'
     static const int InternalNodesEdgeConenctionMap[8][12]={ {-1,-1,6,-1,    -1,-1,5,-1,   -1,-1,3,-1},    //nDownNode=0
                                                              {-1,-1,7,-1,    -1,-1,-1,4,   -1,-1,-1,2},    //nDownNode=1
                                                              {-1,-1,-1,4,    -1,-1,7,-1,   -1,1,-1,-1},    //nDownNode=2
                                                              {-1,-1,-1,5,    -1,-1,-1,6,   0,-1,-1,-1},    //nDownNode=3
                                                              {-1,2,-1,-1,    -1,1,-1,-1,   -1,-1,7,-1},    //nDownNode=4
                                                              {-1,3,-1,-1,    0,-1,-1,-1,   -1,-1,-1,6},    //nDownNode=5
                                                              {0,-1,-1,-1,    -1,3,-1,-1,   -1,5,-1,-1},    //nDownNode=6
                                                              {1,-1,-1,-1,    2,-1,-1,-1,   4,-1,-1,-1} };   //nDownNode=7


     static const int InternalNodesEdgeConenctionMap_List[8][3]= { {2,6,10},    //nDownNode=0
                                                                   {2,7,11},    //nDownNode=1
                                                                   {3,6,9},     //nDownNode=2
                                                                   {3,7,8},     //nDownNode=3

                                                                   {1,5,10},    //nDownNode=4
                                                                   {1,4,11},    //nDownNode=5
                                                                   {0,5,9},    //nDownNode=6
                                                                   {0,4,8} };   //nDownNode=7

     static const int  InternalNodesEdgeConenctionMap_List_Length=3;



     for (nDownNode=0;nDownNode<8;nDownNode++) {
       register int iEdge,nedge;

       downNode=startNode->downNode[nDownNode];

       for (nList=0;nList<InternalNodesEdgeConenctionMap_List_Length;nList++) {
         nedge=InternalNodesEdgeConenctionMap_List[nDownNode][nList];
         t=startNode->downNode[InternalNodesEdgeConenctionMap[nDownNode][nedge]];

         for (iEdge=0;iEdge<2;iEdge++) downNode->SetNeibEdge(t,nedge,iEdge);
       }
     }



     //CONNECT THROUGH the FACES OF THE upNODE
     struct cExternalNodesFace {
       int nUpFace;
       int iFace;
       int jFace;
       int nNeibFace;
     };

     static const cExternalNodesFace ExternalNodesFaceConnectionMap_upNode[8][6]={
         {{0,0,0,1},{-1,-1,-1,-1,},{2,0,0,3},{-1,-1,-1,-1},{4,0,0,5},{-1,-1,-1,-1}},    //nDownNode=0
         {{-1,-1,-1,-1}, {1,0,0,0},{2,1,0,3},{-1,-1,-1,-1},{4,1,0,5},{-1,-1,-1,-1}},    //nDownNode=1
         {{0,1,0,1},{-1,-1,-1,-1},{-1,-1,-1,-1},{3,0,0,2},{4,0,1,5},{-1,-1,-1,-1}},     //nDownNode=2
         {{-1,-1,-1,-1},{1,1,0,0},{-1,-1,-1,-1},{3,1,0,2},{4,1,1,5},{-1,-1,-1,-1}},     //nDownNode=3

         {{0,0,1,1},{-1,-1,-1,-1},{2,0,1,3},{-1,-1,-1,-1},{-1,-1,-1,-1},{5,0,0,4}},     //nDownNode=4
         {{-1,-1,-1,-1},{1,0,1,0},{2,1,1,3},{-1,-1,-1,-1},{-1,-1,-1,-1},{5,1,0,4}},     //nDownNode=5
         {{0,1,1,1},{-1,-1,-1,-1},{-1,-1,-1,-1},{3,0,1,2},{-1,-1,-1,-1},{5,0,1,4}},     //nDownNode=6
         {{-1,-1,-1,-1},{1,1,1,0},{-1,-1,-1,-1},{3,1,1,2},{-1,-1,-1,-1},{5,1,1,4}} };   //nDownNode=7

     static const int ExternalNodesFaceConnectionMap_upNode_List[8][3]= { {0,2,4},    //nDownNode=0
                                                                          {1,2,4},    //nDownNode=1
                                                                          {0,3,4},    //nDownNode=2
                                                                          {1,3,4},    //nDownNode=3

                                                                          {0,2,5},    //nDownNode=4
                                                                          {1,2,5},    //nDownNode=5
                                                                          {0,3,5},    //nDownNode=6
                                                                          {1,3,5} };  //nDownNode=7

     static const int ExternalNodesFaceConnectionMap_upNode_List_Length=3;


     for (nDownNode=0;nDownNode<8;nDownNode++) {
       register int nUpFace,iFace,jFace,nNeibFace,nface;

       downNode=startNode->downNode[nDownNode];

       for (nList=0;nList<ExternalNodesFaceConnectionMap_upNode_List_Length;nList++) {
         nface=ExternalNodesFaceConnectionMap_upNode_List[nDownNode][nList];

         nUpFace=ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].nUpFace;
         iFace=ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].iFace;
         jFace=ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].jFace;
         nNeibFace=ExternalNodesFaceConnectionMap_upNode[nDownNode][nface].nNeibFace;

         t=startNode->GetNeibFace(nUpFace,iFace,jFace);

         //connect node 't' to the 'downNode'
         if (t!=NULL) {
           if (t->RefinmentLevel==downNode->RefinmentLevel) {
             for (iFace=0;iFace<2;iFace++) for (jFace=0;jFace<2;jFace++) t->SetNeibFace(downNode,nNeibFace,iFace,jFace);
           }
           else if (t->RefinmentLevel==downNode->RefinmentLevel-1) {
             t->SetNeibFace(downNode,nNeibFace,iFace,jFace);
           }
           else exit(__LINE__,__FILE__,"Error: the mesh is not consistent");
         }

         //connect 'downNode' to the nodes 't'
         for (iFace=0;iFace<2;iFace++) for (jFace=0;jFace<2;jFace++) downNode->SetNeibFace(t,nface,iFace,jFace);
       }
     }



     //CONNECT THROUGH NODES


     //Connect through the corner nodes of the upNode
     struct cExternalNodesCorner_upNode {
       int nUpCornerNode;
       int nNeibCornerNode;
     };


     static const cExternalNodesCorner_upNode ExternalNodesCornerNodeConnectionMap_upNode[8][8]={
         {{0,7},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}},    //nDpwnNode=0
         {{-1,-1},{1,6},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}},    //nDpwnNode=1
         {{-1,-1},{-1,-1},{2,5},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}},    //nDpwnNode=2
         {{-1,-1},{-1,-1},{-1,-1},{3,4},{-1,-1},{-1,-1},{-1,-1},{-1,-1}},    //nDpwnNode=3

         {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{4,3},{-1,-1},{-1,-1},{-1,-1}},    //nDpwnNode=4
         {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{5,2},{-1,-1},{-1,-1}},    //nDpwnNode=5
         {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{6,1},{-1,-1}},    //nDpwnNode=6
         {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{7,0}} };  //nDpwnNode=7


     for (nDownNode=0,nd=0;nDownNode<8;nDownNode++,nd++) {
       downNode=startNode->downNode[nDownNode];
       t=startNode->GetNeibCorner(nd);

       downNode->SetNeibCorner(t,nd);
       if (t!=NULL) t->SetNeibCorner(downNode,ExternalNodesCornerNodeConnectionMap_upNode[nDownNode][nd].nNeibCornerNode);
     }

     //Connect through the nodes created in the middle of the faces of upNode
     struct cExternalNodesCorner_downNode {
       int nUpFace;
       int iFace;
       int jFace;
       int nNeibCornerNode;
     };


     static const cExternalNodesCorner_downNode ExternalNodesCornerNodeConnectionMap_downNode[8][8] = {
         {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{4,1,1,4},  {-1,-1,-1,-1},{2,1,1,2},{0,1,1,1},{-1,-1,-1,-1}},           //nDownNode=0
         {{-1,-1,-1,-1},{-1,-1,-1,-1},{4,0,1,5},{-1,-1,-1,-1},  {2,0,1,3},{-1,-1,-1,-1},{-1,-1,-1,-1},{1,1,1,0}},           //nDownNode=1
         {{-1,-1,-1,-1},{4,1,0,6},{-1,-1,-1,-1},{-1,-1,-1,-1},  {0,0,1,3},{-1,-1,-1,-1},{-1,-1,-1,-1},{3,1,1,0}},           //nDownNode=2
         {{4,0,0,7},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},  {-1,-1,-1,-1},{1,0,1,2},{3,0,1,1},{-1,-1,-1,-1}},           //nDownNode=3

         {{-1,-1,-1,-1},{2,1,0,6},{0,1,0,5},{-1,-1,-1,-1},      {-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{5,1,1,0}},       //nDownNode=4
         {{2,0,0,7},{-1,-1,-1,-1},{-1,-1,-1,-1},{1,1,0,4},      {-1,-1,-1,-1},{-1,-1,-1,-1},{5,0,1,1},{-1,-1,-1,-1}},       //nDownNode=5
         {{0,0,0,7},{-1,-1,-1,-1},{-1,-1,-1,-1},{3,1,0,4},      {-1,-1,-1,-1},{5,1,0,2},{-1,-1,-1,-1},{-1,-1,-1,-1}},       //nDownNode=6
         {{-1,-1,-1,-1},{1,0,0,6},{3,0,0,5},{-1,-1,-1,-1},      {5,0,0,3},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}  };    //nDownNode=7


     static const int ExternalNodesCornerNodeConnectionMap_downNode_List[8][3]= { {3,5,6},    //nDownNode=0
                                                                              {2,4,7},    //nDownNode=1
                                                                              {1,4,7},    //nDownNode=2
                                                                              {0,5,6},    //nDownNode=3

                                                                              {1,2,7},    //nDownNode=4
                                                                              {0,3,6},    //nDownNode=5
                                                                              {0,3,5},    //nDownNode=6
                                                                              {1,2,4}};   //nDownNode=7

     static const int ExternalNodesCornerNodeConnectionMap_downNode_List_Length=3;




     for (nDownNode=0;nDownNode<8;nDownNode++) {
       int iFace,jFace,nUpFace; //,nd;

       downNode=startNode->downNode[nDownNode];

       for (nList=0;nList<ExternalNodesCornerNodeConnectionMap_downNode_List_Length;nList++) {
         nd=ExternalNodesCornerNodeConnectionMap_downNode_List[nDownNode][nList];

         nUpFace=ExternalNodesCornerNodeConnectionMap_downNode[nDownNode][nd].nUpFace;
         iFace=ExternalNodesCornerNodeConnectionMap_downNode[nDownNode][nd].iFace;
         jFace=ExternalNodesCornerNodeConnectionMap_downNode[nDownNode][nd].jFace;

         t=startNode->GetNeibFace(nUpFace,iFace,jFace);

         if (t!=NULL) if (t->RefinmentLevel==downNode->RefinmentLevel) {
           downNode->SetNeibCorner(t,nd);
           t->SetNeibCorner(downNode,ExternalNodesCornerNodeConnectionMap_downNode[nDownNode][nd].nNeibCornerNode);
         }
       }
     }


     //Connect through the nodes created in the middle of the edges of upNode
     struct cExternalNodesEdgeCorner_upNode {
       int nUpEdge;
       int iEdge;
       int nCornerNode;
       int nNeibCornerNode;
     };


     static const cExternalNodesEdgeCorner_upNode ExternalNodesEdgeCornerConnectionMap_upNode[8][12] = {
         {{0,1,1,6},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},    {4,1,2,5},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},   {8,1,4,3},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}},     //nDpwnNode=0
         {{0,0,0,7},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},    {-1,-1,-1,-1},{5,1,3,4},{-1,-1,-1,-1},{-1,-1,-1,-1},   {-1,-1,-1,-1},{9,1,5,2},{-1,-1,-1,-1},{-1,-1,-1,-1}},     //nDpwnNode=1
         {{-1,-1,-1,-1},{1,1,3,4},{-1,-1,-1,-1},{-1,-1,-1,-1},    {4,0,0,7},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},   {-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{11,1,6,1}},    //nDpwnNode=2
         {{-1,-1,-1,-1},{1,0,2,5},{-1,-1,-1,-1},{-1,-1,-1,-1},    {-1,-1,-1,-1},{5,0,1,6},{-1,-1,-1,-1},{-1,-1,-1,-1},   {-1,-1,-1,-1},{-1,-1,-1,-1},{10,1,7,0},{-1,-1,-1,-1}},    //nDpwnNode=3

         {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{3,1,5,2},    {-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{7,1,6,1},   {8,0,0,7},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}},       //nDpwnNode=4
         {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{3,0,4,3},    {-1,-1,-1,-1},{-1,-1,-1,-1},{6,1,7,0},{-1,-1,-1,-1},   {-1,-1,-1,-1},{9,0,1,6},{-1,-1,-1,-1},{-1,-1,-1,-1}},      //nDpwnNode=5
         {{-1,-1,-1,-1},{-1,-1,-1,-1},{2,1,7,0},{-1,-1,-1,-1},    {-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{7,0,4,3},   {-1,-1,-1,-1,},{-1,-1,-1,-1},{-1,-1,-1,-1},{11,0,2,5}},    //nDpwnNode=6
         {{-1,-1,-1,-1},{-1,-1,-1,-1},{2,0,6,1},{-1,-1,-1,-1},    {-1,-1,-1,-1},{-1,-1,-1,-1},{6,0,5,2},{-1,-1,-1,-1},   {-1,-1,-1,-1},{-1,-1,-1,-1},{10,0,3,4},{-1,-1,-1,-1}}  };    //nDpwnNode=7


     static const int ExternalNodesEdgeCornerConnectionMap_upNode_List[8][3]= { {0,4,8},    //nDpwnNode=0
                                                                                {0,5,9},    //nDpwnNode=1
                                                                                {1,4,11},   //nDpwnNode=2
                                                                                {1,5,10},   //nDpwnNode=3

                                                                                {3,7,8},    //nDpwnNode=4
                                                                                {3,6,9},    //nDpwnNode=5
                                                                                {2,7,11},   //nDpwnNode=6
                                                                                {2,6,10}};  //nDpwnNode=7

     static const int ExternalNodesEdgeCornerConnectionMap_upNode_List_Length=3;


     for (nDownNode=0;nDownNode<8;nDownNode++) {
       register int nUpEdge,iEdge,nNeibCornerNode,nedge,nCornerNode;

       downNode=startNode->downNode[nDownNode];

       for (nList=0;nList<ExternalNodesEdgeCornerConnectionMap_upNode_List_Length;nList++) {
         nedge=ExternalNodesEdgeCornerConnectionMap_upNode_List[nDownNode][nList];

         nUpEdge=ExternalNodesEdgeCornerConnectionMap_upNode[nDownNode][nedge].nUpEdge;
         iEdge=ExternalNodesEdgeCornerConnectionMap_upNode[nDownNode][nedge].iEdge;
         nNeibCornerNode=ExternalNodesEdgeCornerConnectionMap_upNode[nDownNode][nedge].nNeibCornerNode;
         nCornerNode=ExternalNodesEdgeCornerConnectionMap_upNode[nDownNode][nedge].nCornerNode;

         t=startNode->GetNeibEdge(nUpEdge,iEdge);

         if (t!=NULL) if (t->RefinmentLevel==downNode->RefinmentLevel) {
           downNode->SetNeibCorner(t,nCornerNode); //connect 't' to 'downNode
           t->SetNeibCorner(downNode,nNeibCornerNode); //connect 'downNode' to 't'
         }
       }
     }



     //CONNECT THROUGH EDGES
      //Connect through the edges formed by the splitting of the edges of upNode

     struct cExternalNodeEdge_upNode {
       int nUpEdge;
       int iEdge;
       int nNeibEdge;
     };


     static const cExternalNodeEdge_upNode ExternalNodesEdgeConnectionMap_upNode[8][12]= {
         {{0,0,2},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},    {4,0,6},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},    {8,0,10},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},     //nDpwnNode=0
         {{0,1,2},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},    {-1,-1,-1},{5,0,7},{-1,-1,-1},{-1,-1,-1},    {-1,-1,-1},{9,0,11},{-1,-1,-1},{-1,-1,-1}},     //nDpwnNode=1
         {{-1,-1,-1},{1,0,3},{-1,-1,-1},{-1,-1,-1},    {4,1,6},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},    {-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{11,0,9}},    //nDpwnNode=2
         {{-1,-1,-1},{1,1,3},{-1,-1,-1},{-1,-1,-1},    {-1,-1,-1},{5,1,7},{-1,-1,-1},{-1,-1,-1},    {-1,-1,-1},{-1,-1,-1},{10,0,8},{-1,-1,-1}},   //nDpwnNode=3

         {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{3,0,1},    {-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{7,0,5},    {8,1,10},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},    //nDpwnNode=4
         {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{3,1,1},    {-1,-1,-1},{-1,-1,-1},{6,0,4},{-1,-1,-1},    {-1,-1,-1},{9,1,11},{-1,-1,-1},{-1,-1,-1}},    //nDpwnNode=5
         {{-1,-1,-1},{-1,-1,-1},{2,0,0},{-1,-1,-1},    {-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{7,1,5},    {-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{11,1,9}},     //nDpwnNode=6
         {{-1,-1,-1},{-1,-1,-1},{2,1,0},{-1,-1,-1},    {-1,-1,-1},{-1,-1,-1},{6,1,4},{-1,-1,-1},    {-1,-1,-1},{-1,-1,-1},{10,1,8},{-1,-1,-1}} };  //nDpwnNode=7


     static const int ExternalNodesEdgeConnectionMap_upNode_List[8][3]={ {0,4,8},   //nDpwnNode=0
                                                                         {0,5,9},   //nDpwnNode=1
                                                                         {1,4,11},  //nDpwnNode=2
                                                                         {1,5,10},  //nDpwnNode=3

                                                                         {3,7,8},   //nDpwnNode=4
                                                                         {3,6,9},   //nDpwnNode=5
                                                                         {2,7,11},  //nDpwnNode=6
                                                                         {2,6,10}}; //nDpwnNode=7

     static const int ExternalNodesEdgeConnectionMap_upNode_List_Length=3;

     for (nDownNode=0;nDownNode<8;nDownNode++) {
       int /*nUpEdge,*/ iEdge,nedge,nNeibEdge;

       downNode=startNode->downNode[nDownNode];

       for (nList=0;nList<ExternalNodesEdgeConnectionMap_upNode_List_Length;nList++) {
         nedge=ExternalNodesEdgeConnectionMap_upNode_List[nDownNode][nList];

//         nUpEdge=ExternalNodesEdgeConnectionMap_upNode[nDownNode][nedge].nUpEdge;
         iEdge=ExternalNodesEdgeConnectionMap_upNode[nDownNode][nedge].iEdge;
         nNeibEdge=ExternalNodesEdgeConnectionMap_upNode[nDownNode][nedge].nNeibEdge;

         t=startNode->GetNeibEdge(nedge,iEdge);

         //connect 'downNode' to 't'
         if (t!=NULL) {
           if (t->RefinmentLevel==startNode->RefinmentLevel) t->SetNeibEdge(downNode,nNeibEdge,iEdge);
           else if (t->RefinmentLevel==downNode->RefinmentLevel) for (iEdge=0;iEdge<2;iEdge++) t->SetNeibEdge(downNode,nNeibEdge,iEdge);
           else exit(__LINE__,__FILE__,"Error: the mesh is not consistent");
         }

         //connect 't' to 'downNode'
         for (iEdge=0;iEdge<2;iEdge++) downNode->SetNeibEdge(t,nedge,iEdge);
       }

     }


     //Connect through the new edges created on the faces of the upNode
     struct cExternalNodesEdgeCorner_downNode {
       int nUpFace;
       int iFace;
       int jFace;
       int nNeibEdge;
     };

     static const cExternalNodesEdgeCorner_downNode ExternalNodesEdgeCornerConnectionMap_downNode[8][12]=
     { {{-1,-1,-1,-1},{4,0,1,3},{-1,-1,-1,-1},{2,0,1,1},    {-1,-1,-1,-1},{4,1,0,7},{-1,-1,-1,-1},{0,0,1,5},    {-1,-1,-1,-1},{2,1,0,11},{-1,-1,-1},{0,1,0,9}},     //nDpwnNode=0
       {{-1,-1,-1,-1},{4,1,1,3},{-1,-1,-1,-1},{2,1,1,1},    {4,0,0,6},{-1,-1,-1,-1},{1,0,1,4},{-1,-1,-1,-1},    {2,0,0,10},{-1,-1,-1,-1},{1,1,0,8},{-1,-1,-1,-1}},  //nDpwnNode=1
       {{4,0,0,2},{-1,-1,-1,-1},{3,0,1,0},{-1,-1,-1,-1},    {-1,-1,-1,-1},{4,1,1,7},{-1,-1,-1,-1},{0,1,1,5},    {0,0,0,10},{-1,-1,-1,-1},{3,1,0,8},{-1,-1,-1,-1}},    //nDpwnNode=2
       {{4,1,0,2},{-1,-1,-1,-1},{3,1,1,0},{-1,-1,-1,-1},    {4,0,1,6},{-1,-1,-1,-1},{1,1,1,4},{-1,-1,-1,-1},    {-1,-1,-1,-1},{1,0,0,11},{-1,-1,-1,-1},{3,0,0,9}},    //nDpwnNode=3

       {{2,0,0,2},{-1,-1,-1,-1},{5,0,1,0},{-1,-1,-1,-1},    {0,0,0,6},{-1,-1,-1,-1},{5,1,0,4},{-1,-1,-1,-1},    {-1,-1,-1,-1},{2,1,1,11},{-1,-1,-1,-1},{0,1,1,9}},    //nDpwnNode=4
       {{2,1,0,2},{-1,-1,-1,-1},{5,1,1,0},{-1,-1,-1,-1},    {-1,-1,-1,-1},{1,0,0,7},{-1,-1,-1,-1},{5,0,0,5},    {2,0,1,10},{-1,-1,-1,-1},{1,1,1,8},{-1,-1,-1,-1}},    //nDpwnNode=5
       {{-1,-1,-1,-1},{3,0,0,3},{-1,-1,-1,-1},{5,0,0,1},    {0,1,0,6},{-1,-1,-1,-1},{5,1,1,4},{-1,-1,-1,-1},    {0,0,1,10,},{-1,-1,-1,-1},{3,1,1,8},{-1,-1,-1,-1}},   //nDpwnNode=6
       {{-1,-1,-1,-1},{3,1,0,3},{-1,-1,-1,-1},{5,1,0,1},    {-1,-1,-1,-1},{1,1,0,7},{-1,-1,-1,-1},{5,0,1,5},    {-1,-1,-1,-1},{1,0,1,11},{-1,-1,-1,-1},{3,0,1,9}} };  //nDpwnNode=7


     static const int ExternalNodesEdgeCornerConnectionMap_downNode_List[8][6]= { {1,3,5,7,9,11},    //nDpwnNode=0
                                                                                  {1,3,4,6,8,10},     //nDpwnNode=1
                                                                                  {0,2,5,7,8,10},    //nDpwnNode=2
                                                                                  {0,2,4,6,9,11},    //nDpwnNode=3

                                                                                  {0,2,4,6,9,11},    //nDpwnNode=4
                                                                                  {0,2,5,7,8,10},    //nDpwnNode=5
                                                                                  {1,3,4,6,8,10},    //nDpwnNode=6
                                                                                  {1,3,5,7,9,11}};   //nDpwnNode=7

     static const int ExternalNodesEdgeCornerConnectionMap_downNode_List_Length=6;

     for (nDownNode=0;nDownNode<8;nDownNode++) {
       register int nUpFace,iFace,jFace,nNeibEdge,nedge,iEdge;

       downNode=startNode->downNode[nDownNode];

       for (nList=0;nList<ExternalNodesEdgeCornerConnectionMap_downNode_List_Length;nList++) {
         nedge=ExternalNodesEdgeCornerConnectionMap_downNode_List[nDownNode][nList];

         nUpFace=ExternalNodesEdgeCornerConnectionMap_downNode[nDownNode][nedge].nUpFace;
         iFace=ExternalNodesEdgeCornerConnectionMap_downNode[nDownNode][nedge].iFace;
         jFace=ExternalNodesEdgeCornerConnectionMap_downNode[nDownNode][nedge].jFace;
         nNeibEdge=ExternalNodesEdgeCornerConnectionMap_downNode[nDownNode][nedge].nNeibEdge;

         t=startNode->GetNeibFace(nUpFace,iFace,jFace);

         //connect 'downNode' to 't' together
         if (t!=NULL) if (t->RefinmentLevel==downNode->RefinmentLevel)  for (iEdge=0;iEdge<2;iEdge++) {
           downNode->SetNeibEdge(t,nedge,iEdge);
           t->SetNeibEdge(downNode,nNeibEdge,iEdge);
         }

       }
     }


     //check consistency of the node's reconnection
#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
#if _CHECK_MESH_CONSISTANCY_ == _ON_AMR_MESH_
     for (nd=0;nd<(1<<_MESH_DIMENSION_);nd++) if ((neibNode=startNode->neibNodeCorner[nd])!=NULL) {
       if (neibNode->CheckNeibNode(startNode)==true) exit(__LINE__,__FILE__,"Error: the mesh is not consistent");
     }

     for (register int nface=0;nface<2*_MESH_DIMENSION_;nface++) for (i=0;i<2;i++) for (j=0;j<2;j++) {
       neibNode=startNode->GetNeibFace(nface,i,j);
       if (neibNode!=NULL) if (neibNode->CheckNeibNode(startNode)==true) exit(__LINE__,__FILE__,"Error: the mesh is not consistent");
     }

     for (register int nedge=0;nedge<12;nedge++) for (register int iEdge=0;iEdge<2;iEdge++) {
       neibNode=startNode->GetNeibEdge(nedge,iEdge);
       if (neibNode!=NULL) if (neibNode->CheckNeibNode(startNode)==true) exit(__LINE__,__FILE__,"Error: the mesh is not consistent");
     }
#endif
#endif

     /*
     //clean connections of startNode
     for (nd=0;nd<(1<<_MESH_DIMENSION_);nd++) startNode->SetNeibCorner(NULL,nd);
     for (register int nface=0;nface<2*_MESH_DIMENSION_;nface++) for (i=0;i<2;i++) for (j=0;j<2;j++) startNode->SetNeibFace(NULL,nface,i,j);
     for (register int nedge=0;nedge<12;nedge++) for (register int iEdge=0;iEdge<2;iEdge++) startNode->SetNeibEdge(NULL,nedge,iEdge);
*/

#endif





     //deallocate the block of 'startNode'
     if (DeallocateUnusedBlocks==true) DeallocateBlock(startNode);

     //remove 'startNode' from the list of nodes that intersects the surface of the computational domain
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
     if (startNode->InternalBoundaryDescriptorList!=NULL) {
       cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptor,*nextDesecriptor;

       for (InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;InternalBoundaryDescriptor!=NULL;InternalBoundaryDescriptor=nextDesecriptor) {
         nextDesecriptor=InternalBoundaryDescriptor->nextInternalBCelement;
         InternalBoundaryDescriptors.deleteElement(InternalBoundaryDescriptor);
       }

       startNode->InternalBoundaryDescriptorList=NULL;

       cTreeNodeAMR<cBlockAMR> *next,*prev;

       next=startNode->DomainSurfaceBoundaryList_Next,prev=startNode->DomainSurfaceBoundaryList_Prev;
       if (next!=NULL) next->DomainSurfaceBoundaryList_Prev=prev;
       if (prev!=NULL) prev->DomainSurfaceBoundaryList_Next=next;

       if (DomainSurfaceBoundaryList==startNode) DomainSurfaceBoundaryList=next;

       startNode->DomainSurfaceBoundaryList_Next=NULL,startNode->DomainSurfaceBoundaryList_Prev=NULL;
     }
#endif



  }

  return true;
}


  //reset the node's 'nodeProcessedFlag' 
  void resetNodeProcessedFlag(bool resetMaxRefinmentLevel=true) {
    long int nMemoryBank,nTotalMemoryBanks,nnode;


    //reset the flag for the 'corner nodes'
    cCornerNode *cornerNodeDataBuffer;
    nTotalMemoryBanks=CornerNodes.dataBufferListPointer;   

    for (nMemoryBank=0;nMemoryBank<nTotalMemoryBanks;nMemoryBank++) {
      cornerNodeDataBuffer=CornerNodes.dataBufferList[nMemoryBank];

      for (nnode=0;nnode<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;nnode++) {
        (cornerNodeDataBuffer+nnode)->nodeDescriptor.nodeProcessedFlag=_AMR_FALSE_;
        if (resetMaxRefinmentLevel==true) (cornerNodeDataBuffer+nnode)->nodeDescriptor.maxRefinmentLevel=0;
      }
    }

    //reset the flag for the 'center nodes'
    #if _AMR_CENTER_NODE_ ==  _ON_AMR_MESH_
    //reset the flag for the 'corner nodes'
    cCenterNode *centerNodeDataBuffer;
    nTotalMemoryBanks=CenterNodes.dataBufferListPointer;

    for (nMemoryBank=0;nMemoryBank<nTotalMemoryBanks;nMemoryBank++) {
      centerNodeDataBuffer=CenterNodes.dataBufferList[nMemoryBank];

      for (nnode=0;nnode<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;nnode++) {
        (centerNodeDataBuffer+nnode)->nodeDescriptor.nodeProcessedFlag=_AMR_FALSE_;
        if (resetMaxRefinmentLevel==true) (centerNodeDataBuffer+nnode)->nodeDescriptor.maxRefinmentLevel=0;
      }
    }
    #endif


  } 

  void CompareTreeStructure(cTreeNodeAMR<cBlockAMR> *startNode) {
    static cAMRnodeID id,idtemp;
    static CMPI_channel pipe(1000000);

    if (startNode==rootTree) {
      if (ThisThread==0)  pipe.openRecvAll();
      else pipe.openSend(0);
    }

    GetAMRnodeID(id,startNode);

    if (ThisThread==0) {
      for (int thread=1;thread<nTotalThreads;thread++) {
        pipe.recv(idtemp,thread);
        if (id!=idtemp) exit(__LINE__,__FILE__,"Error: the tree is not consistent");
      }
    }
    else pipe.send(id);

    for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) CompareTreeStructure(startNode->downNode[nDownNode]);

    if (startNode==rootTree) {
      if (ThisThread==0)  pipe.closeRecvAll();
      else pipe.closeSend();

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
  }

  void CompareGlobalNodeNumbering(cTreeNodeAMR<cBlockAMR> *startNode) {
    static cAMRnodeID id,idtemp;
    static CMPI_channel pipe(1000000);

    if (startNode==rootTree) {
      if (ThisThread==0)  pipe.openRecvAll();
      else pipe.openSend(0);
    }

    GetAMRnodeID(id,startNode);


    if (ThisThread==0) {
      for (int thread=1;thread<nTotalThreads;thread++) {
        pipe.recv(idtemp,thread);
        if (id!=idtemp) exit(__LINE__,__FILE__,"Error: the tree is not consistent");
      }
    }
    else pipe.send(id);


    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      int i,j,k,thread;
      long int nodeno,t,nd;
      cCornerNode *node;


      if (ThisThread==0) {
        for (thread=1;thread<nTotalThreads;thread++) {
          pipe.recv(t,thread);

          if (t!=startNode->Thread) exit(__LINE__,__FILE__,"Error: the parallel threads distribution is not consistent");
        }
      }
      else {
        t=startNode->Thread;
        pipe.send(t);
      }



#if _MESH_DIMENSION_ == 1
      static const int iMax=1+_BLOCK_CELLS_X_,jMax=1,kMax=1;
#elif _MESH_DIMENSION_ == 2
      static const int iMax=1+_BLOCK_CELLS_X_,jMax=1+_BLOCK_CELLS_Y_,kMax=1;
#else
      static const int iMax=1+_BLOCK_CELLS_X_,jMax=1+_BLOCK_CELLS_Y_,kMax=1+_BLOCK_CELLS_Z_;
#endif

      for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
        if (startNode->block==NULL) nodeno=-1;
        else {
          nd=getCornerNodeLocalNumber(i,j,k);
          node=startNode->block->GetCornerNode(nd);
          nodeno=node->nodeDescriptor.nodeno;
        }

        if (ThisThread==0) {
          for (thread=1;thread<nTotalThreads;thread++) {
            pipe.recv(t,thread);

            if (nodeno==-1) nodeno=t;
            else if ((t!=-1)&&(nodeno!=t)) exit(__LINE__,__FILE__,"Error: the global node's numbering is not consistent");
          }
        }
        else pipe.send(nodeno);


//================  DEBUG ==========================
//if (ThisThread!=0) pipe.flush();
//MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


//================= END DEBUG ========================


      }

    }
    else for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) CompareGlobalNodeNumbering(startNode->downNode[nDownNode]);

    if (startNode==rootTree) {
      if (ThisThread==0)  pipe.closeRecvAll();
      else pipe.closeSend();

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
  }



  void countMeshElements(cTreeNodeAMR<cBlockAMR> *startNode,int level) {
    static long int nDownNodes,nd;
    int i,j,k;
    cBlockAMR *block=NULL;
    cCornerNode *ndptr=NULL;
    
    static int SendRequest,*SendRequestVector=NULL;
    static long int *BlockNodeGlabalNumber=NULL;

//======================  DEBUG =======================
//cCornerNode *foundNodePtr;
static long int CallsCounter=0;

CallsCounter++;
/*
if (CallsCounter==83) {
  *DiagnospticMessageStream << __LINE__ << std::endl;
}
*/
//====================== END DEBUG ==============



#if _MESH_DIMENSION_ == 1
    static const int iMax=1+_BLOCK_CELLS_X_,jMax=1,kMax=1;
    static int nBlockCornerNodes=1+_TOTAL_BLOCK_CELLS_X_;
#elif _MESH_DIMENSION_ == 2
    static const int iMax=1+_BLOCK_CELLS_X_,jMax=1+_BLOCK_CELLS_Y_,kMax=1;
    static int nBlockCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
#else
    static const int iMax=1+_BLOCK_CELLS_X_,jMax=1+_BLOCK_CELLS_Y_,kMax=1+_BLOCK_CELLS_Z_;
    static int nBlockCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
#endif


    if (startNode==rootTree) {
      meshBlocksNumber=0,meshNodesNumber=0,meshMaximumRefinmentLevel=0;
      resetNodeProcessedFlag();

      //compare the tree
      CompareTreeStructure(rootTree);

      SendRequestVector=new int [nTotalThreads];
      BlockNodeGlabalNumber=new long int [nBlockCornerNodes];


      //the procedure is developed only for the case where the domain is covered by the layer of boundary blocks
#if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
      //do nothing
#else
      exit(__LINE__,__FILE__,"The procedure is implemented only for the case _AMR_PARALLEL_DATA_EXCHANGE_MODE_ = _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_");
#endif

      nDownNodes=1<<_MESH_DIMENSION_;
    }

    //check the value of the maximum refinment level
    if (level>meshMaximumRefinmentLevel) meshMaximumRefinmentLevel=level;

    //count the number of blocks on the mesh 
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      meshBlocksNumber++;
      SendRequest=false;

#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
      MPI_Bcast(&meshNodesNumber,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
#endif



//======================  DEBUG =======================
/*
      if (startNode->Temp_ID==1684) {
        *DiagnospticMessageStream << __LINE__ << std::endl;
      }

      {
        double xNodeFind[3]={-175.0,-300.0,-187.5};
        double *xNodeTest;
        int thread;


        block=startNode->block;
        SendRequest=false;

        foundNodePtr=NULL;

        if (block==NULL) SendRequest=false;
        else {
          for (i=0;i<6;i++) for (j=0;j<6;j++) for (k=0;k<6;k++) {

            ndptr=block->GetCornerNode(getCornerNodeLocalNumber(i,j,k));
            xNodeTest=ndptr->GetX();

            if ((fabs(xNodeTest[0]-xNodeFind[0])<EPS)&&(fabs(xNodeTest[1]-xNodeFind[1])<EPS)&&(fabs(xNodeTest[2]-xNodeFind[2])<EPS)) SendRequest=true,foundNodePtr=ndptr;
          }
        }

        MPI_Gather(&SendRequest,1,MPI_INT,SendRequestVector,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

        if (ThisThread==0) for (thread=0;thread<nTotalThreads;thread++) if (SendRequestVector[thread]==true) SendRequest=true;

        MPI_Bcast(&SendRequest,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

        if (SendRequest==true) {
          if (ThisThread==0) *DiagnospticMessageStream << __LINE__ << std::endl;
        }

      }

      SendRequest=false;
      */
//====================== END DEBUG ==============

      //count the mesh nodes associated with the block: the nodes must be counted in the order they will be output inth the data file
      if ((startNode->Thread==ThisThread)&&(startNode->block!=NULL)) {
        block=startNode->block;

        #if _MESH_DIMENSION_ == 1
        static const int isubBlockMax=_BLOCK_CELLS_X_;
        static const int jsubBlockMax=1,ksubBlockMax=1;
        static const int nBasicBlockNodes=2;
        static const int nodeOutputOrder_I[2]={0,1};
        static const int nodeOutputOrder_J[2]={0,0};
        static const int nodeOutputOrder_K[2]={0,0};

        #elif _MESH_DIMENSION_ == 2
        static const int isubBlockMax=_BLOCK_CELLS_X_;
        static const int jsubBlockMax=_BLOCK_CELLS_Y_;
        static const int ksubBlockMax=1;
        static const int nBasicBlockNodes=4;
        static const int nodeOutputOrder_I[4]={0,1,1,0};
        static const int nodeOutputOrder_J[4]={0,0,1,1};
        static const int nodeOutputOrder_K[4]={0,0,0,0};
        #else
        static const int isubBlockMax=_BLOCK_CELLS_X_;
        static const int jsubBlockMax=_BLOCK_CELLS_Y_;
        static const int ksubBlockMax=_BLOCK_CELLS_Z_;
        static const int nBasicBlockNodes=8;
        static const int nodeOutputOrder_I[8]={0,1,1,0,0,1,1,0};
        static const int nodeOutputOrder_J[8]={0,0,1,1,0,0,1,1};
        static const int nodeOutputOrder_K[8]={0,0,0,0,1,1,1,1};
        #endif


//        long int nd;
        int ksubBlock,jsubBlock,isubBlock,iNode,jNode,kNode,nnode;

        for (ksubBlock=0;ksubBlock<ksubBlockMax;ksubBlock++) for (jsubBlock=0;jsubBlock<jsubBlockMax;jsubBlock++) for (isubBlock=0;isubBlock<isubBlockMax;isubBlock++) {
          for (nnode=0;nnode<nBasicBlockNodes;nnode++) {
            iNode=nodeOutputOrder_I[nnode]+isubBlock;
            jNode=nodeOutputOrder_J[nnode]+jsubBlock;
            kNode=nodeOutputOrder_K[nnode]+ksubBlock;

            nd=getCornerNodeLocalNumber(iNode,jNode,kNode);
            if ((ndptr=block->GetCornerNode(nd))==NULL) continue;

            if (ndptr->nodeDescriptor.nodeProcessedFlag==_AMR_FALSE_) {
              ++meshNodesNumber;
              if (_MAX_MESH_ELEMENT_NUMBER_==meshNodesNumber) exit(__LINE__,__FILE__,"Error: the maximum values of the mesh element counter is reached. Increase the value of _MESH_ELEMENTS_NUMBERING_BITS_");

              ndptr->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_;
              ndptr->nodeDescriptor.nodeno=meshNodesNumber;
            }

            if (ndptr->nodeDescriptor.maxRefinmentLevel<level) ndptr->nodeDescriptor.maxRefinmentLevel=level;
          }
        }
      }
      else {
        //prepare the send request -> request recievinf of the block if on the currect processor startNode->Neib[:]->block!=NULL
        cTreeNodeAMR<cBlockAMR>* neibNode;

#if _MESH_DIMENSION_ == 1
        if (startNode->block!=NULL) SendRequest=true;

        if (SendRequest==false) {
          for (int nface=0;nface<2;nface++) if (SendRequest==false) {
            if ((neibNode=startNode->GetNeibFace(nface,0,0))!=NULL) if (neibNode->block!=NULL) {
              SendRequest=true;
              break;
            }
          }
        }
#elif _MESH_DIMENSION_ == 2
        if (startNode->block!=NULL) SendRequest=true;

        if (SendRequest==false) {
          for (nd=0;nd<4;nd++) if ((neibNode=startNode->GetNeibCorner(nd))!=NULL) if (neibNode->block!=NULL) {
            SendRequest=true;
            break;
          }
        }

        if (SendRequest==false) {
          for (int nface=0;nface<4;nface++) if (SendRequest==false) {
            for (int iFace=0;iFace<2;iFace++) if ((neibNode=startNode->GetNeibFace(nface,iFace,0))!=NULL) if (neibNode->block!=NULL) {
              SendRequest=true;
              break;
            }
          }
        }
#else

        if (startNode->block!=NULL) SendRequest=true;

        if (SendRequest==false) {
          for (nd=0;nd<8;nd++) if ((neibNode=startNode->GetNeibCorner(nd))!=NULL) if (neibNode->block!=NULL) {
            SendRequest=true;
            break;
          }
        }

        if (SendRequest==false) {
          for (int nedge=0;nedge<12;nedge++) if (SendRequest==false) for (int iEdge=0;iEdge<2;iEdge++) if ((neibNode=startNode->GetNeibEdge(nedge,iEdge))!=NULL) if (neibNode->block!=NULL) {
            SendRequest=true;
            break;
          }
        }

        if (SendRequest==false) {
          for (int nface=0;nface<6;nface++) if (SendRequest==false) {
            for (int iFace=0;iFace<2;iFace++) for (int jFace=0;jFace<2;jFace++) if ((neibNode=startNode->GetNeibFace(nface,iFace,jFace))!=NULL) if (neibNode->block!=NULL) {
              SendRequest=true;
              break;
            }
          }
        }

#endif
      }

#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
      if (startNode->Thread!=0) {
        if (ThisThread==0) {
          MPI_Status status;

          //update the total node's number
          MPI_Recv(&meshNodesNumber,1,MPI_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
        }
        else if (ThisThread==startNode->Thread) {
          MPI_Send(&meshNodesNumber,1,MPI_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
        }
      }

      //collect the SendRequest vector and exchenge the data
      MPI_Gather(&SendRequest,1,MPI_INT,SendRequestVector,1,MPI_INT,startNode->Thread,MPI_GLOBAL_COMMUNICATOR);

      if (ThisThread==startNode->Thread) {
        //compose the vector of global node's number (v)
        int thread;
        bool flag=false;

        for (thread=0;thread<nTotalThreads;thread++) if (SendRequestVector[thread]==true) {
          flag=true;
          break;
        }


        if (flag==true) {
          for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
            nd=getCornerNodeLocalNumber(i,j,k);
            BlockNodeGlabalNumber[nd]=((ndptr=block->GetCornerNode(nd))!=NULL) ? ndptr->nodeDescriptor.nodeno : -1;
          }

          //send the node's global number's vector
          for (thread=0;thread<nTotalThreads;thread++) if (SendRequestVector[thread]==true) {
            MPI_Send(BlockNodeGlabalNumber,nBlockCornerNodes,MPI_LONG,thread,0,MPI_GLOBAL_COMMUNICATOR);
          }
        }
      }
      else if (SendRequest==true) {
        MPI_Status status;
        bool TempBlockAllocationFlag;

        //recieve the vector of node's numbers
        MPI_Recv(BlockNodeGlabalNumber,nBlockCornerNodes,MPI_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        //allocate the block
        TempBlockAllocationFlag=(startNode->block==NULL) ? true : false;

        if (TempBlockAllocationFlag==true) AllocateBlock(startNode);
        block=startNode->block;

        //set up the nodes' numbers
        for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
          nd=getCornerNodeLocalNumber(i,j,k);
          if ((ndptr=block->GetCornerNode(nd))!=NULL) {
            if ((ndptr->nodeDescriptor.nodeProcessedFlag==_AMR_TRUE_)&&(ndptr->nodeDescriptor.nodeno!=BlockNodeGlabalNumber[nd])) exit(__LINE__,__FILE__,"Error: reinitialization of nodes' number");

            ndptr->nodeDescriptor.nodeno=BlockNodeGlabalNumber[nd];
            ndptr->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_;
          }
        }

        //deallocate the block
        if (TempBlockAllocationFlag==true) DeallocateBlock(startNode);
      }


//============================ DEBUG ======================
//      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
//      CompareGlobalNodeNumbering(rootTree);
//      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
//============================ END DEBUG ====================

#endif

    }
    else {
      //*DiagnospticMessageStream the nodes in the downNodes
      for (i=0;i<nDownNodes;i++) if (startNode->downNode[i]!=NULL) countMeshElements(startNode->downNode[i],level+1);
    }

    //reset the 'meshModifiedFlag'
    if (startNode==rootTree) {
      delete [] SendRequestVector;
      delete [] BlockNodeGlabalNumber;

      SendRequestVector=NULL,BlockNodeGlabalNumber=NULL;

      CompareGlobalNodeNumbering(rootTree);

      meshModifiedFlag_CountMeshElements=false;

      MPI_Bcast(&meshNodesNumber,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    }
  } 


  long int getMeshNodesNumber() {
    if (meshModifiedFlag_CountMeshElements==true) countMeshElements(rootTree,0);

    return meshNodesNumber;
  }  

  int getMeshMaxRefinmentLevel() {
    if (meshModifiedFlag_CountMeshElements==true) countMeshElements(rootTree,0);

    return meshMaximumRefinmentLevel;
  } 


  //build the mesh
  bool buildMesh_OneLevelRefinment(cTreeNodeAMR<cBlockAMR> *startNode,int level,int startLevel) {
    bool res=false,flag;
    double c,blockMiddlePoint[3]={0.0,0.0,0.0},requredResolution,characteristicBlockSize;
    int idim,i,j,k;
    double xProbe[3];

    #if _MESH_DIMENSION_ == 1 
    static const double characteristicBlockSize_max=fabs(dxRootBlock[0]/_BLOCK_CELLS_X_);
    static const double characteristicBlockSize_min=characteristicBlockSize_max/(1<<_MAX_REFINMENT_LEVEL_);
    #elif _MESH_DIMENSION_ == 2 
    static const double characteristicBlockSize_max=sqrt(pow(dxRootBlock[0]/_BLOCK_CELLS_X_,2)+pow(dxRootBlock[1]/_BLOCK_CELLS_Y_,2));
    static const double characteristicBlockSize_min=characteristicBlockSize_max/(1<<_MAX_REFINMENT_LEVEL_);
    #else 
    static const double characteristicBlockSize_max=sqrt(pow(dxRootBlock[0]/_BLOCK_CELLS_X_,2)+pow(dxRootBlock[1]/_BLOCK_CELLS_Y_,2)+pow(dxRootBlock[2]/_BLOCK_CELLS_Z_,2));
    static const double characteristicBlockSize_min=characteristicBlockSize_max/(1<<_MAX_REFINMENT_LEVEL_);
    #endif 


    callMpiBarrier();

    if (level==startLevel) { //refine the block if needed
      //get the characterist size of the block and the position of the block's middle point 
      characteristicBlockSize=characteristicBlockSize_max/(1<<level);

      for (idim=0;idim<_MESH_DIMENSION_;idim++) blockMiddlePoint[idim]=startNode->xmin[idim]+0.5*(xGlobalMax[idim]-xGlobalMin[idim])/(1<<level);

      //get the requested mesh resolution 
      requredResolution=localResolution(blockMiddlePoint);

      #if _MESH_DIMENSION_ == 1
      static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
      static const int jMin=0,jMax=0,kMin=0,kMax=0;
      #elif _MESH_DIMENSION_ == 2
      static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
      static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
      static const int kMin=0,kMax=0;
      #elif _MESH_DIMENSION_ == 3
      static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
      static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
      static const int kMin=-_GHOST_CELLS_Z_,kMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
      #endif


      for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)  {
        //evaluate the required resolution within the computational domain
        startNode->GetCornerNodePosition(xProbe,i,j,k);
        c=localResolution(xProbe);
        if (c<requredResolution) requredResolution=c;

        //evaluate the required resolution on the surfaces installed into the mesh
        #if _INTERNAL_BOUNDARY_MODE_  == _INTERNAL_BOUNDARY_MODE_ON_
        if (startNode->InternalBoundaryDescriptorList!=NULL) {
          cInternalBoundaryConditionsDescriptor *SurfaceDescriptor;
          double (*SurfaceLocalResolution)(double*);

          SurfaceLocalResolution=NULL;

          for (SurfaceDescriptor=startNode->InternalBoundaryDescriptorList;SurfaceDescriptor!=NULL;SurfaceDescriptor=SurfaceDescriptor->nextInternalBCelement) {
            if (SurfaceDescriptor->BondaryType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) SurfaceLocalResolution=((cInternalSphericalData*)SurfaceDescriptor->BoundaryElement)->localResolution;
            else if (SurfaceDescriptor->BondaryType==_INTERNAL_BOUNDARY_TYPE_CIRCLE_) SurfaceLocalResolution=((cInternalCircleData*)SurfaceDescriptor->BoundaryElement)->localResolution;
            else if (SurfaceDescriptor->BondaryType==_INTERNAL_BOUNDARY_TYPE_1D_SPHERE_) SurfaceLocalResolution=((cInternalSphere1DData*)SurfaceDescriptor->BoundaryElement)->localResolution;
            else if (SurfaceDescriptor->BondaryType==_INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_) SurfaceLocalResolution=((cInternalRotationBodyData*)SurfaceDescriptor->BoundaryElement)->localResolution;
            else if (SurfaceDescriptor->BondaryType==_INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_) SurfaceLocalResolution=((cInternalNastranSurfaceData*)SurfaceDescriptor->BoundaryElement)->localResolution;


            else exit(__LINE__,__FILE__,"Error: unknown boundary type");

            if (SurfaceLocalResolution!=NULL) {
              c=SurfaceLocalResolution(xProbe);
              if (c<requredResolution) requredResolution=c;
            }

          }
        }
        #endif

        //evaluate the required local resolution are the surface defined by the cut-faces
        #if _AMR__CUT_CELL__MODE_ ==  _AMR__CUT_CELL__MODE__ON_
        if ((CutCellSurfaceLocalResolution!=NULL)&&(startNode->FirstTriangleCutFace!=NULL)) {
          CutCell::cTriangleFaceDescriptor* t;
          int cnt=0;

          #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
          static int nThreadsOpenMP=-1;
          static double *requredResolutionTable=NULL;
          int iThreadOpenMP;

          if (requredResolutionTable==NULL) {
            #pragma omp parallel default (none) shared (nThreadsOpenMP)
            {
              #pragma omp single
              {
              nThreadsOpenMP=omp_get_num_threads();
              }
            }

            requredResolutionTable=new double [nThreadsOpenMP];
          }

          for (iThreadOpenMP=0;iThreadOpenMP<nThreadsOpenMP;iThreadOpenMP++) requredResolutionTable[iThreadOpenMP]=requredResolution;

          #pragma omp parallel default(none) private (t,c,iThreadOpenMP,cnt) firstprivate(cnt,startNode) shared (requredResolutionTable,nThreadsOpenMP)
          {
            #pragma omp single
            {
              for (cnt=0,t=startNode->FirstTriangleCutFace;t!=NULL;t=t->next,cnt++) if ((cnt%nTotalThreads==ThisThread)||(ParallelMeshGenerationFlag==false)) {
                 #pragma omp task default (none) firstprivate (node)
                {
                  iThreadOpenMP=omp_get_thread_num();
                  c=CutCellSurfaceLocalResolution(t->TriangleFace);
                  if (c<requredResolutionTable[iThreadOpenMP]) requredResolutionTable[iThreadOpenMP]=c;
                }
              }
            }
          }

          //get the minimum requested resolution
          for (iThreadOpenMP=1,requredResolution=requredResolutionTable[0];iThreadOpenMP<nThreadsOpenMP;iThreadOpenMP++) {
            if (requredResolution>requredResolutionTable[iThreadOpenMP]) requredResolution=requredResolutionTable[iThreadOpenMP];
          }
          #else //_COMPILATION_MODE__HYBRID_
          for (cnt=0,t=startNode->FirstTriangleCutFace;t!=NULL;t=t->next,cnt++) if ((cnt%nTotalThreads==ThisThread)||(ParallelMeshGenerationFlag==false)) {
            c=CutCellSurfaceLocalResolution(t->TriangleFace);
            if (c<requredResolution) requredResolution=c;
          }
          #endif //_COMPILATION_MODE__HYBRID_

          //combine min requsted resolution determined by all MPI processes
          if (ParallelMeshGenerationFlag==true) {
            double Table[nTotalThreads];

            MPI_Gather(&requredResolution,1,MPI_DOUBLE,Table,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
            if (ThisThread==0) for (int thread=0;thread<nTotalThreads;thread++) if (Table[thread]<requredResolution) requredResolution=Table[thread];
            MPI_Bcast(&requredResolution,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
          }
        }
        #endif //_AMR__CUT_CELL__MODE_


        #if _AMR_ENFORCE_CELL_RESOLUTION_MODE_ == _AMR_ENFORCE_CELL_RESOLUTION_MODE_ON_
        if (requredResolution<characteristicBlockSize_min) exit(__LINE__,__FILE__,"The required resolution is smaller than the minimum resolution allowed for the mesh. Increase the value of _MAX_REFINMENT_LEVEL_");
        #endif
      }  



      //the block is split if 1. the cell size exceeds the required one and 2. the root block has to be split (startNode->upNode==NULL)
      if ((requredResolution<characteristicBlockSize)||(startNode->upNode==NULL)) {
        res=splitTreeNode(startNode);

      } 
    }
    else { // go to the downNodes  
      int nDownNode;

      for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) {
        flag=buildMesh_OneLevelRefinment(startNode->downNode[nDownNode],level+1,startLevel);
        if (flag==true) res=true;
      }
    }

    return res;
  } 

  void buildMesh() {
    int level;
    bool flag;

    //begin function message
    int rank;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);

    if (rank==0) {
      std::cout << "$PREFIX: Building mesh.....  " << std::endl << std::flush;
    }
 
    for (level=0;level<=_MAX_REFINMENT_LEVEL_;level++) {
      flag=buildMesh_OneLevelRefinment(rootTree,0,level);

      #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_ 
      #if _CHECK_MESH_CONSISTANCY_ == _ON_AMR_MESH_
      checkMeshConsistency(rootTree); 
      #endif
      #endif

      if (rank==0) {
        printf("\r$PREFIX: Building mesh [completed level: %i]\n",level);
        fflush(stdout);
      }

      if (flag==false) break;
    }

    //determine the resolution limits for each block
    SetNodeNeibResolutionLevelLimit();

    //end function message
    if (rank==0) {
      printf("$PREFIX: Building mesh.....   done\n");
      fflush(stdout);
    }
  }  
  
//==============================================================
//calcualte the interpolation coefficients
  int CenterNodesInterpolationCoefficients_3D_linear(double *x,double *CoefficientsList,cCenterNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients) {

    //if the length of the coefficient list is not enough -> exist with an error message
    if (nMaxCoefficients<8) {
	  exit(__LINE__,__FILE__,"The length of the interpolation stencil is too short");
  	  return -1;
    }

    //determine the local coordinates of the point
    double xLocal[3],xProbe[3],dx,dy,dz,*xmin=startNode->xmin,*xmax=startNode->xmax;


    if ((x[0]<xmin[0]-EPS)||(x[0]>xmax[0]+EPS) || (x[1]<xmin[1]-EPS)||(x[1]>xmax[1]+EPS) || (x[2]<xmin[2]-EPS)||(x[2]>xmax[2]+EPS)) {
      exit(__LINE__,__FILE__,"The point is outside of the block");
      return -1;
    }

    //determine the offset of the interpolating stencil related to the origin of the block
    int counter,i,j,k,iProbeIndex,jProbeIndex,kProbeIndex;
    long int nd;
    cCenterNode *cell;
    double InterpolationWeight,WeightNorm=0.0;
    cTreeNodeAMR<cBlockAMR>* cellNode;

    dx=0.1*(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
    dy=0.1*(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
    dz=0.1*(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;


    for (i=-1,counter=0;i<2;i+=2) {
      xProbe[0]=x[0]+i*dx;

      for (j=-1;j<2;j+=2) {
        xProbe[1]=x[1]+j*dy;

    	  for (k=-1;k<2;k+=2) {
    	    xProbe[2]=x[2]+k*dz;

    	    if ((cellNode=findTreeNode(xProbe,startNode))==NULL) continue;

    	    //determine the indexes of the center node that corresponds to the point 'xProbe'
    	    cellNode->ConvertGlobal2LocalCoordinates(xLocal,xProbe);

    	    iProbeIndex=(xLocal[0]>=0.0) ? (int)(xLocal[0]*_BLOCK_CELLS_X_) : -1+(int)(xLocal[0]*_BLOCK_CELLS_X_);
    	    jProbeIndex=(xLocal[1]>=0.0) ? (int)(xLocal[1]*_BLOCK_CELLS_Y_) : -1+(int)(xLocal[1]*_BLOCK_CELLS_Y_);
    	    kProbeIndex=(xLocal[2]>=0.0) ? (int)(xLocal[2]*_BLOCK_CELLS_Z_) : -1+(int)(xLocal[2]*_BLOCK_CELLS_Z_);

    	    if ((nd=getCenterNodeLocalNumber(iProbeIndex,jProbeIndex,kProbeIndex))==-1) continue;
    	    cell=cellNode->block->GetCenterNode(nd);

    	    if (cell!=NULL) if (cell->Measure>0.0) {
    	      //calculate the interpolation weight
    	      InterpolationWeight=1.0;


//============  DEBUG ==================
/*
//if cell->r == 0 -> check if the nodes is in the domain boundary list
if (_MESH_DIMENSION_ == 3)  if ((cell->r<0.0001)&&(fabs(cell->GetX()[0])+fabs(cell->GetX()[1])+fabs(cell->GetX()[2])>0.0001)) {
  bool found=false;

  cTreeNodeAMR<cBlockAMR> *ptr;

  for (ptr=DomainBoundaryLayerNodesList[cellNode->Thread];ptr!=NULL;ptr=ptr->nextNodeThisThread) if (ptr==cellNode) found=true;

//  exit(__LINE__,__FILE__,"Error: not initialized");
}
*/
//============  END DEBUG ==============

    	      //update the interpolation stencil
      	    CoefficientsList[counter]=InterpolationWeight;
    	      InterpolationStencil[counter]=cell;
    	      counter+=1;
    	      WeightNorm+=InterpolationWeight;
          }
      	}
      }
    }


    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    //do nothing
    #elif  _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
    if (counter==0) exit(__LINE__,__FILE__,"There is no nodes defined");
    #else 
    exit(__LINE__,__FILE__,"Error: the option is not defined");
    #endif

    if (counter!=0) for (i=0;i<counter;i++) CoefficientsList[i]/=WeightNorm;

    return counter;
  }


  int CenterNodesInterpolationCoefficients_2D_linear(double *x,double *CoefficientsList,cCenterNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients) {

    //if the length of the coefficient list is not enough -> exist with an error message
    if (nMaxCoefficients<4) {
    exit(__LINE__,__FILE__,"The length of the interpolation stencil is too short");
      return -1;
    }

    //determine the local coordinates of the point
    double xLocal[2],xProbe[2],dx,dy,*xmin=startNode->xmin,*xmax=startNode->xmax;


    if ((x[0]<xmin[0]-EPS)||(x[0]>xmax[0]+EPS) || (x[1]<xmin[1]-EPS)||(x[1]>xmax[1]+EPS)) {
      exit(__LINE__,__FILE__,"The point is outside of the block");
      return -1;
    }

    //determine the offset of the interpolating stencil related to the origin of the block
    int counter,i,j,iProbeIndex,jProbeIndex;
    long int nd;
    cCenterNode *cell;
    double InterpolationWeight,WeightNorm=0.0;
    cTreeNodeAMR<cBlockAMR>* cellNode;

    dx=0.1*(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
    dy=0.1*(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;


    for (i=-1,counter=0;i<2;i+=2) {
      xProbe[0]=x[0]+i*dx;

      for (j=-1;j<2;j+=2) {
        xProbe[1]=x[1]+j*dy;

        if ((cellNode=findTreeNode(xProbe,startNode))==NULL) continue;

        //determine the indexes of the center node that corresponds to the point 'xProbe'
        cellNode->ConvertGlobal2LocalCoordinates(xLocal,xProbe);

        iProbeIndex=(xLocal[0]>=0.0) ? (int)(xLocal[0]*_BLOCK_CELLS_X_) : -1+(int)(xLocal[0]*_BLOCK_CELLS_X_);
        jProbeIndex=(xLocal[1]>=0.0) ? (int)(xLocal[1]*_BLOCK_CELLS_Y_) : -1+(int)(xLocal[1]*_BLOCK_CELLS_Y_);

        if ((nd=getCenterNodeLocalNumber(iProbeIndex,jProbeIndex,0))==-1) continue;
        cell=cellNode->block->GetCenterNode(nd);

        if (cell!=NULL) if (cell->Measure>0.0) {
          //calculate the interpolation weight
          InterpolationWeight=1.0;

          //update the interpolation stencil
          CoefficientsList[counter]=InterpolationWeight;
          InterpolationStencil[counter]=cell;
          counter+=1;
          WeightNorm+=InterpolationWeight;
        }
      }
    }


    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    //do nothing
    #elif  _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
    if (counter==0) exit(__LINE__,__FILE__,"There is no nodes defined");
    #else
    exit(__LINE__,__FILE__,"Error: the option is not defined");
    #endif

    if (counter!=0) for (i=0;i<counter;i++) CoefficientsList[i]/=WeightNorm;

    return counter;
  }

  int CenterNodesInterpolationCoefficients_1D_linear(double *x,double *CoefficientsList,cCenterNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients) {

    //if the length of the coefficient list is not enough -> exist with an error message
    if (nMaxCoefficients<2) {
    exit(__LINE__,__FILE__,"The length of the interpolation stencil is too short");
      return -1;
    }

    //determine the local coordinates of the point
    double xLocal[3]={0.0,0.0,0.0},xProbe[3]={0.0,0.0,0.0},dx,*xmin=startNode->xmin,*xmax=startNode->xmax;


    if ((x[0]<xmin[0]-EPS)||(x[0]>xmax[0]+EPS)) {
      exit(__LINE__,__FILE__,"The point is outside of the block");
      return -1;
    }

    //determine the offset of the interpolating stencil related to the origin of the block
    int counter,i,iProbeIndex;
    long int nd;
    cCenterNode *cell;
    double InterpolationWeight,WeightNorm=0.0;
    cTreeNodeAMR<cBlockAMR>* cellNode;

    dx=0.1*(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;

    for (i=-1,counter=0;i<2;i+=2) {
      xProbe[0]=x[0]+i*dx;

      if ((cellNode=findTreeNode(xProbe,startNode))==NULL) continue;

      //determine the indexes of the center node that corresponds to the point 'xProbe'
      cellNode->ConvertGlobal2LocalCoordinates(xLocal,xProbe);

        iProbeIndex=(xLocal[0]>=0.0) ? (int)(xLocal[0]*_BLOCK_CELLS_X_) : -1+(int)(xLocal[0]*_BLOCK_CELLS_X_);

        if ((nd=getCenterNodeLocalNumber(iProbeIndex,0,0))==-1) continue;
        cell=cellNode->block->GetCenterNode(nd);

        if (cell!=NULL) if (cell->Measure>0.0) {
          //calculate the interpolation weight
          InterpolationWeight=1.0;

          //update the interpolation stencil
          CoefficientsList[counter]=InterpolationWeight;
          InterpolationStencil[counter]=cell;
          counter+=1;
          WeightNorm+=InterpolationWeight;
        }

    }


    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    //do nothing
    #elif  _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
    if (counter==0) exit(__LINE__,__FILE__,"There is no nodes defined");
    #else
    exit(__LINE__,__FILE__,"Error: the option is not defined");
    #endif

    if (counter!=0) for (i=0;i<counter;i++) CoefficientsList[i]/=WeightNorm;

    return counter;
  }


  /*
     int CenterNodesInterpolationCoefficients_3D_linear(double *x,double *CoefficientsList,cCenterNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients) {

    //if the length of the coefficient list is not enough -> exist with an error message
    if (nMaxCoefficients<8) {
    exit(__LINE__,__FILE__,"The length of the interpolation stencil is too short");
      return -1;
    }

    //determine the local coordinates of the point
    double xLocal,yLocal,zLocal,*xmin=startNode->xmin,*xmax=startNode->xmax;

    xLocal=_BLOCK_CELLS_X_*(x[0]-xmin[0])/(xmax[0]-xmin[0]);
    yLocal=_BLOCK_CELLS_Y_*(x[1]-xmin[1])/(xmax[1]-xmin[1]);
    zLocal=_BLOCK_CELLS_Z_*(x[2]-xmin[2])/(xmax[2]-xmin[2]);

    if ((xLocal<0.0)||(xLocal>_BLOCK_CELLS_X_)||(yLocal<0.0)||(yLocal>_BLOCK_CELLS_Y_)||(zLocal<0.0)||(zLocal>_BLOCK_CELLS_Z_)) {
      exit(__LINE__,__FILE__,"The point is outside of the block");
      return -1;
    }

    //determine the offset of the interpolating stencil related to the origin of the block
    int i,j,k,ioffset,joffset,koffset,counter;
    long int nd;
    cCenterNode *ptr;

    ioffset=int(xLocal+0.5)-1;
    joffset=int(yLocal+0.5)-1;
    koffset=int(zLocal+0.5)-1;

    //calculate the local coordinated of 'x' within the stencil
    double c,cx,cy,dx,dy,dz;
    double cTotal=0.0;

    dx=xLocal+0.5-(ioffset+1);
    dy=yLocal+0.5-(joffset+1);
    dz=zLocal+0.5-(koffset+1);

    for (i=0,counter=0;i<2;i++) {
      cx=(i==0) ? 1.0-dx : dx;

      for (j=0;j<2;j++) {
      cy=(j==0) ? 1.0-dy : dy;

      for (k=0;k<2;k++) {
        c=((k==0) ? 1.0-dz : dz)*cx*cy;
        nd=getCenterNodeLocalNumber(i+ioffset,j+joffset,k+koffset);
          ptr=startNode->block->GetCenterNode(nd);

        if (ptr!=NULL) if (ptr->Measure>0.0) {
          CoefficientsList[counter]=c;
          InterpolationStencil[counter]=ptr;
          counter+=1;
          cTotal+=c;
          }
      }
      }
    }

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((counter==8)&&(fabs(cTotal-1.0)>1.0E-6)) exit(__LINE__,__FILE__,"The summ of the coefficients is different from one");
    #endif


    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    //do nothing
    #elif  _INTERNAL_BOUNDARY_MODE_ _INTERNAL_BOUNDARY_MODE_OFF_
    if (counter==0) exit(__LINE__,__FILE__,"There is no nodes defined");
    #else
    exit(__LINE__,__FILE__,"Error: the option is not defined");
    #endif

    if (counter!=0) for (int i=0;i<counter;i++) CoefficientsList[i]/=cTotal;

    return counter;
  }
   */
//==============================================================
  //if printCoordinateVector == true  -> print the coordinate vector, else -> printf the connectovity list  
  void outputMeshTECPLOT_BlockCornerNode_BlockConnectivityList(cTreeNodeAMR<cBlockAMR> *startNode,FILE* fout,bool printCoordinateVector,bool PrintMeshData,int DataSetNumber,int ConnectivityListMode) {
    int isubBlock,jsubBlock,ksubBlock; ///,isubBlockMax,jsubBlockMax,ksubBlockMax;
    int iNode,jNode,kNode,idim,nnode; //,nBasicBlockNodes;
    double  xNode[3];
    bool flag=true;

    static CMPI_channel pipe(1000000);
    static long int nGlobalNodeNumber=0;


    if (startNode==rootTree) {
      nGlobalNodeNumber=1;
      ResetAMRnodeProcessingFlag();
      resetNodeProcessedFlag();

      if ((ConnectivityListMode!=_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_)&&(ConnectivityListMode!=_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_)) exit(__LINE__,__FILE__,"Error: unknown option");

      if (ThisThread==0) pipe.openRecvAll();
      else pipe.openSend(0);
    }

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {

      flag=false;
      if (ThisThread==0) flag=true;
      if ((ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_)&&(ThisThread==startNode->Thread)) flag=true;
      if ((ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_)&&(startNode->block!=NULL)) flag=true;

      if (flag==true) {
        //the node is on the bottom of the tree -> printf the node
        //get the limits for the number of the sub-blocks in the block
        #if _MESH_DIMENSION_ == 1
        static const int isubBlockMax=_BLOCK_CELLS_X_;
        static const int jsubBlockMax=1,ksubBlockMax=1;
        static const int nBasicBlockNodes=2;
        static const int nodeOutputOrder_I[2]={0,1};
        static const int nodeOutputOrder_J[2]={0,0};
        static const int nodeOutputOrder_K[2]={0,0};

        #elif _MESH_DIMENSION_ == 2
        static const int isubBlockMax=_BLOCK_CELLS_X_;
        static const int jsubBlockMax=_BLOCK_CELLS_Y_;
        static const int ksubBlockMax=1;
        static const int nBasicBlockNodes=4;
        static const int nodeOutputOrder_I[4]={0,1,1,0};
        static const int nodeOutputOrder_J[4]={0,0,1,1};
        static const int nodeOutputOrder_K[4]={0,0,0,0};
        #else
        static const int isubBlockMax=_BLOCK_CELLS_X_;
        static const int jsubBlockMax=_BLOCK_CELLS_Y_;
        static const int ksubBlockMax=_BLOCK_CELLS_Z_;
        static const int nBasicBlockNodes=8;
        static const int nodeOutputOrder_I[8]={0,1,1,0,0,1,1,0};
        static const int nodeOutputOrder_J[8]={0,0,1,1,0,0,1,1};
        static const int nodeOutputOrder_K[8]={0,0,0,0,1,1,1,1};
        #endif

        cBlockAMR *block=startNode->block;
        cCornerNode *cornerNode;
        long int nd;

        for (ksubBlock=0;ksubBlock<ksubBlockMax;ksubBlock++) for (jsubBlock=0;jsubBlock<jsubBlockMax;jsubBlock++) for (isubBlock=0;isubBlock<isubBlockMax;isubBlock++) {
          for (nnode=0;nnode<nBasicBlockNodes;nnode++) {
            iNode=nodeOutputOrder_I[nnode]+isubBlock;
            jNode=nodeOutputOrder_J[nnode]+jsubBlock;
            kNode=nodeOutputOrder_K[nnode]+ksubBlock;

            startNode->GetCornerNodePosition(xNode,iNode,jNode,kNode);
            nd=getCornerNodeLocalNumber(iNode,jNode,kNode);
            cornerNode=(block!=NULL) ? block->GetCornerNode(nd) : NULL;

            //if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_): Check if the node is already printed, skip printed nodes
            if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_) {
              if ((cornerNode!=NULL)&&(printCoordinateVector==true)) {
                flag=(cornerNode->nodeDescriptor.nodeProcessedFlag==_AMR_TRUE_) ? true : false;
                cornerNode->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_;
              }
              else flag=false;

              if (startNode->Thread==0) {
                if ((ThisThread!=0)||(flag==true)) continue;
              }
              else {
                if ((ThisThread!=0)&&(ThisThread!=startNode->Thread)) continue;

                if (ThisThread!=0) pipe.send(flag);
                else pipe.recv(flag,startNode->Thread);

                if (flag==true) continue;
              }
            }


            if (printCoordinateVector==true) {  //print the nodes' locations /and the data stored on the mesh
              if (ThisThread==0) for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(fout,"%e  ",xNode[idim]);


//====================================== DEBUG ===================

              /*
              if (_MESH_DIMENSION_==3) if ((fabs(xNode[0]+300.0)<EPS)&&(fabs(xNode[1]-1000.0)<EPS)&&(fabs(xNode[2]+00.0)<EPS)) {
  *DiagnospticMessageStream << __FILE__ << __LINE__ << std::endl;
}

              if (startNode->Temp_ID==1116) {
                *DiagnospticMessageStream << __FILE__ << __LINE__ << std::endl;
              }

*/
//=========================== END DEBUG ====================



              if (PrintMeshData==true) { //print the data stored on the mesh
                //print basic parameters
                long int MaxRefinmentLevel,NodeTempID;


//=========================   DEBUG ==============================
/*
                if (cornerNode!=NULL) if (cornerNode->Temp_ID==58786) {
                  *DiagnospticMessageStream << __FILE__ << __LINE__ << std::endl;
                }
*/
//========================= END DEBUG =========================


                if (ThisThread!=0) {
                  MaxRefinmentLevel=cornerNode->nodeDescriptor.maxRefinmentLevel,NodeTempID=cornerNode->Temp_ID;
                  pipe.send(MaxRefinmentLevel);
                  pipe.send(NodeTempID);
                }
                else {
                  if (startNode->Thread==0) MaxRefinmentLevel=cornerNode->nodeDescriptor.maxRefinmentLevel,NodeTempID=cornerNode->Temp_ID;
                  else {
                    pipe.recv(MaxRefinmentLevel,startNode->Thread);
                    pipe.recv(NodeTempID,startNode->Thread);
                  }

                  fprintf(fout,"%ld  %ld %i  ",MaxRefinmentLevel,NodeTempID,startNode->Thread);
                }


                //print the "corner" nodes
#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
                if (cornerNode!=NULL) {
                  cornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                }
                else if (ThisThread==0) { //the root processor
                  cCornerNode *tempCornerNode=CornerNodes.newElement();

                  tempCornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                  CornerNodes.deleteElement(tempCornerNode);
                }
                else exit(__LINE__,__FILE__,"Error: something is wrong");
#else
                cornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
#endif




                //print the data stored in the 'center' nodes
                #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
                const int nMaxCenterInterpolationCoefficients=64;
                cCenterNode *tempCenterNode,*CenterNodeInterpolationStencil[nMaxCenterInterpolationCoefficients];
                double CenterNodeInterpolationCoefficients[nMaxCenterInterpolationCoefficients];
                int centerNodeInterpolationStencilLength;

                if (ThisThread==0) fprintf(fout,"  ");

                if (ThisThread==startNode->Thread) {
                  tempCenterNode=CenterNodes.newElement();
                  tempCenterNode->SetX(xNode);

                  if (GetCenterNodesInterpolationCoefficients==NULL) {
#if _MESH_DIMENSION_ == 1
                    centerNodeInterpolationStencilLength=CenterNodesInterpolationCoefficients_1D_linear(xNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
#elif _MESH_DIMENSION_ == 2
                    centerNodeInterpolationStencilLength=CenterNodesInterpolationCoefficients_2D_linear(xNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
#elif _MESH_DIMENSION_ == 3
                    centerNodeInterpolationStencilLength=CenterNodesInterpolationCoefficients_3D_linear(xNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
#endif
                  }
                  else {
                    centerNodeInterpolationStencilLength=GetCenterNodesInterpolationCoefficients(xNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
                  }

                  if (centerNodeInterpolationStencilLength==-1) exit(__LINE__,__FILE__,"Error in interpolation of the 'center' node's data");

                  tempCenterNode->Interpolate(CenterNodeInterpolationStencil,CenterNodeInterpolationCoefficients,centerNodeInterpolationStencilLength);
                  tempCenterNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                  CenterNodes.deleteElement(tempCenterNode);
                }
                else if (ThisThread==0) {
                  tempCenterNode=CenterNodes.newElement();
                  tempCenterNode->SetX(xNode);

                  tempCenterNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                  CenterNodes.deleteElement(tempCenterNode);
                }
                else exit(__LINE__,__FILE__,"Error: something is wrong");
                #endif



                //print data stored by the tree node (startNode)
                if (cornerNode!=NULL) {
                   startNode->block->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                }
                else if (ThisThread==0) {
                   cBlockAMR *tempBlock=blocks.newElement();

                   startNode->block->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                   blocks.deleteElement(tempBlock);
                }
                else exit(__LINE__,__FILE__,"Error: something is wrong");




                //end of the printing of the data stored on the mesh
              }




              if (ThisThread==0) fprintf(fout,"\n"); //end of the line
            }
            else { // print the connectivity list

              if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_) {
                if (ThisThread==0) {
                  fprintf(fout,"%ld ",nGlobalNodeNumber);
                  ++nGlobalNodeNumber;
                  if (nnode==nBasicBlockNodes-1) fprintf(fout,"\n");
                }
              }
              else if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_) {
                long int nodeno=0;

                if (cornerNode!=NULL) nodeno=cornerNode->nodeDescriptor.nodeno;

                if (startNode->Thread!=0) {
                  if (ThisThread!=0) pipe.send(nodeno);
                  else pipe.recv(nodeno,startNode->Thread);
                }

                if (ThisThread==0) {
                  fprintf(fout,"%ld ",nodeno);



//###################  DEBUG #################
static long int cnttemp=0;

cnttemp++;
//###################  END DEBUG #################

                  if (nnode==nBasicBlockNodes-1) {

//###################  DEBUG #################
cnttemp=0;
//###################  END DEBUG #################

                    fprintf(fout,"\n");
                  }
                }

                }
              else exit(__LINE__,__FILE__,"Error: unknown option");

            }
          }

        }
      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) {
        outputMeshTECPLOT_BlockCornerNode_BlockConnectivityList(startNode->downNode[nDownNode],fout,printCoordinateVector,PrintMeshData,DataSetNumber,ConnectivityListMode);
      }
    }

    if (startNode==rootTree) {
      if (ThisThread==0) {
        pipe.closeRecvAll();
      }
      else {
        pipe.closeSend();
      }

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }

  } 
 
/*
   void outputMeshTECPLOT_BlockCornerNode(cTreeNodeAMR<cBlockAMR> *startNode,FILE* fout,bool printCoordinateVector,bool PrintMeshData,int DataSetNumber) {
    int isubBlock,jsubBlock,ksubBlock,isubBlockMax,jsubBlockMax,ksubBlockMax;
    int iNode=0,jNode=0,kNode=0,nd,idim,nnode,nBasicBlockNodes;

    static CMPI_channel pipe(1000000);
    static long int nGlobalNodeNumber=0;


    if (startNode==rootTree) {
      nGlobalNodeNumber=1;

      if (ThisThread==0) pipe.openRecvAll();
      else pipe.openSend(0);
    }


//##########################   DEBUG #######################
    static long int nMPIops=0;
//##########################  END DEBUG ####################





    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      if ((ThisThread==0)||(startNode->block!=NULL)){ //the node is on the bottom of the tree -> printf the node
        cBlockAMR *block=startNode->block;
        cCornerNode *cornerNode;

        if (startNode->Thread!=0) {
          if (ThisThread==0) MPI_Send(&nGlobalNodeNumber,1,MPI_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR);
          else if (ThisThread==startNode->Thread) {
            MPI_Status status;

            MPI_Recv(&nGlobalNodeNumber,1,MPI_LONG,0,0,MPI_GLOBAL_COMMUNICATOR,&status);
          }
        }

        //get the limits for the number of the sub-blocks in the block
#if _MESH_DIMENSION_ == 1
        isubBlockMax=_BLOCK_CELLS_X_;
        jsubBlockMax=1,ksubBlockMax=1;
        nBasicBlockNodes=2;
        static int nodeOutputOrder_I[2]={0,1};
        static int nodeOutputOrder_J[2]={0,0};
        static int nodeOutputOrder_K[2]={0,0};
 #elif _MESH_DIMENSION_ == 2
        isubBlockMax=_BLOCK_CELLS_X_;
        jsubBlockMax=_BLOCK_CELLS_Y_;
        ksubBlockMax=1;
        nBasicBlockNodes=4;
        static int nodeOutputOrder_I[4]={0,1,1,0};
        static int nodeOutputOrder_J[4]={0,0,1,1};
        static int nodeOutputOrder_K[4]={0,0,0,0};
#else
        isubBlockMax=_BLOCK_CELLS_X_;
        jsubBlockMax=_BLOCK_CELLS_Y_;
        ksubBlockMax=_BLOCK_CELLS_Z_;
        nBasicBlockNodes=8;
        static int nodeOutputOrder_I[8]={0,1,1,0,0,1,1,0};
        static int nodeOutputOrder_J[8]={0,0,1,1,0,0,1,1};
        static int nodeOutputOrder_K[8]={0,0,0,0,1,1,1,1};
#endif


        for (ksubBlock=0;ksubBlock<ksubBlockMax;ksubBlock++) for (jsubBlock=0;jsubBlock<jsubBlockMax;jsubBlock++) for (isubBlock=0;isubBlock<isubBlockMax;isubBlock++) {
          for (nnode=0;nnode<nBasicBlockNodes;nnode++) {
            iNode=nodeOutputOrder_I[nnode]+isubBlock;
            jNode=nodeOutputOrder_J[nnode]+jsubBlock;
            kNode=nodeOutputOrder_K[nnode]+ksubBlock;

            nd=getCornerNodeLocalNumber(iNode,jNode,kNode);
            cornerNode=(block!=NULL) ? block->GetCornerNode(nd) : NULL;


//############  DEBUG ################



//############ ENS DEBUG ##############



            if (printCoordinateVector==true) {
              bool printflag=false;
              double xNode[3];
              long int MaxRefinmentLevel=0,NodeTempID=0;

              if (cornerNode!=NULL) if ((cornerNode->nodeDescriptor.nodeProcessedFlag==_AMR_FALSE_)&&(cornerNode->nodeDescriptor.maxRefinmentLevel==startNode->RefinmentLevel)) printflag=true;

              if (startNode->Thread!=0) {
                if (ThisThread==startNode->Thread) pipe.send(printflag);
                else if (ThisThread==0) pipe.recv(printflag,startNode->Thread);
                else exit(__LINE__,__FILE__,"Error: something is wrong");
              }



              if (printflag==true) {
                if (cornerNode!=NULL) {
                  cornerNode->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_;
                  cornerNode->nodeDescriptor.nodeno=nGlobalNodeNumber++;
                }

                if (ThisThread==0) {
                  startNode->GetCornerNodePosition(xNode,iNode,jNode,kNode);
                  for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(fout,"%e  ",((fabs(xNode[idim])>EPS) ? xNode[idim] : 0.0));
                }

                if (cornerNode!=NULL) MaxRefinmentLevel=cornerNode->nodeDescriptor.maxRefinmentLevel,NodeTempID=cornerNode->Temp_ID;

#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
                if (true) { //(nTotalThreads!=1) {
                  if (ThisThread!=0) {
                    if (startNode->Thread==ThisThread) {
                      pipe.send(MaxRefinmentLevel);
                      pipe.send(NodeTempID);

//##########################   DEBUG #######################
nMPIops++;
//##########################  END DEBUG ####################

                    }
                  }
                  else { //the root processor
                    if (startNode->Thread!=0) {
                      pipe.recv(MaxRefinmentLevel,startNode->Thread);
                      pipe.recv(NodeTempID,startNode->Thread);

//##########################   DEBUG #######################
nMPIops++;
//##########################  END DEBUG ####################


                    }

                    fprintf(fout,"%ld  %ld %i  ",MaxRefinmentLevel,NodeTempID,startNode->Thread);
                  }
                }
#else
                for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(fout,"%e  ",((fabs(cornerNode->GetX()[idim])>EPS) ? cornerNode->GetX()[idim] : 0.0));
                fprintf(fout,"%i  %i  %ld ",cornerNode->nodeDescriptor.maxRefinmentLevel, ((startNode->block->GetGhostFlag()==_GHOST_BLOCK_) ? 1 : 0),cornerNode->Temp_ID);
#endif



                //print the data stared in the block's nodes
                if ((PrintMeshData==true)&&((startNode->Thread==ThisThread)||(ThisThread==0))) {


                  //print the data stared in the corner nodes
#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
                  if (true) { //(nTotalThreads!=1) {
                    if (ThisThread==ThisThread) {
                      if (startNode->Thread!=0) {
                        cornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                      }
                    }
                    else { //the root processor
                      if (startNode->Thread!=0) {
                        cCornerNode *tempCornerNode=CornerNodes.newElement();

                        tempCornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                        CornerNodes.deleteElement(tempCornerNode);
                      }
                      else cornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                    }
                  }


                  if (cornerNode!=NULL) cornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
#else
                  cornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
#endif


                  //print the data stored in the 'center' nodes
                  #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
                  const int nMaxCenterInterpolationCoefficients=64;
                  cCenterNode *tempCenterNode,*CenterNodeInterpolationStencil[nMaxCenterInterpolationCoefficients];
                  double CenterNodeInterpolationCoefficients[nMaxCenterInterpolationCoefficients],xCenterNode[3];
                  int centerNodeInterpolationStencilLength;

                  fprintf(fout,"  ");
                  cornerNode->GetX(xCenterNode);
                  tempCenterNode=CenterNodes.newElement();

                  if (GetCenterNodesInterpolationCoefficients==NULL) {
                    centerNodeInterpolationStencilLength=CenterNodesInterpolationCoefficients_3D_linear(xCenterNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
                  }
                  else {
                    centerNodeInterpolationStencilLength=GetCenterNodesInterpolationCoefficients(xCenterNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
                  }

                  if (centerNodeInterpolationStencilLength==-1) exit(__LINE__,__FILE__,"Errror in interpolation of the 'center' node's data");


                  tempCenterNode->Interpolate(CenterNodeInterpolationStencil,CenterNodeInterpolationCoefficients,centerNodeInterpolationStencilLength);
                  tempCenterNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                  #endif

                  CenterNodes.deleteElement(tempCenterNode);

                  //print data stored by the tree node (startNode)
#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
                  if (true) { //(nTotalThreads!=1) {
                    if (ThisThread!=0) {
                      if (startNode->Thread==ThisThread) {
                        startNode->block->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                      }
                    }
                    else { //the root processor
                      if (startNode->Thread!=0) {
                        cBlockAMR *tempBlock=blocks.newElement();

                        startNode->block->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                        blocks.deleteElement(tempBlock);
                      }
                      else startNode->block->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
                    }
                  }


                  if (cornerNode!=NULL) cornerNode->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
#else
                  startNode->block->PrintData(fout,DataSetNumber,&pipe,startNode->Thread);
#endif
                }

                //finish the line
                if (ThisThread==0) fprintf(fout,"\n");






              }
            }
            else {
              long int nodeno=-1;

              if (startNode->Thread==ThisThread) nodeno=cornerNode->nodeDescriptor.nodeno;


              if (startNode->Thread!=0) {
                if (startNode->Thread==ThisThread) pipe.send(nodeno);
                else if (ThisThread==0) pipe.recv(nodeno,startNode->Thread);
              }


              if (ThisThread==0) fprintf(fout,"%ld ",nodeno);
            }
          }

          if (printCoordinateVector==false) if (ThisThread==0) fprintf(fout,"\n");
        }

        //update the nodes counter
        if (startNode->Thread!=0) {
          if (ThisThread==startNode->Thread) {
            pipe.send(nGlobalNodeNumber);
            pipe.flush();
          }
          else if (ThisThread==0) {
            pipe.recv(nGlobalNodeNumber,startNode->Thread);
          }
        }

      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) outputMeshTECPLOT_BlockCornerNode(startNode->downNode[nDownNode],fout,printCoordinateVector,PrintMeshData,DataSetNumber);
    }


    if (startNode==rootTree) {
      if (ThisThread==0) pipe.closeRecvAll();
      else pipe.closeSend();

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }

  }


 */

  void outputMeshTECPLOTinternal(const char *fname,const bool PrintMeshData, int DataSetNumber,int ConnectivityListMode) {
    FILE *fout=NULL;

    //the procedure is developed only for the case where the domain is covered by the layer of boundary blocks
#if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
    //do nothing
#else
    exit(__LINE__,__FILE__,"The procedure is implemented only for the case _AMR_PARALLEL_DATA_EXCHANGE_MODE_ = _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_");
#endif


    //Count the number of the elements of the mesh
    if (meshModifiedFlag_CountMeshElements==true) {
      countMeshElements(rootTree,0);
      CompareGlobalNodeNumbering(rootTree);
    }

    if ((ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_)&&(meshNodesNumber==0)) {
      PrintErrorLog(__LINE__,__FILE__,"Warning: no mesh blocks allocated. Changed to ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_");
      ConnectivityListMode=_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_;
    }

    if (ThisThread==0) {
      fout=fopen(fname,"w");

#if _MESH_DIMENSION_ == 1
      fprintf(fout,"VARIABLES=\"X\"");

      if (PrintMeshData==true) {
        fprintf(fout,", \"Maximum Refinment Level\", \"Temp_ID\"");

        #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
        fprintf(fout,", \"Thread\"");
        #endif

        if (CornerNodes.usedElements()==0) exit(__LINE__,__FILE__,"Error: CornerNodes are not allocated");
        CornerNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);

        #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
        if (CenterNodes.usedElements()==0) exit(__LINE__,__FILE__,"Error: CenterNodes are not allocated");
        CenterNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);
        #endif

        rootTree->block->PrintVariableList(fout);
      }

      fprintf(fout,"\n");

#elif _MESH_DIMENSION_ == 2
      fprintf(fout,"VARIABLES=\"X\", \"Y\"");

      if (PrintMeshData==true) {
        fprintf(fout,", \"Maximum Refinment Level\", \"Temp_ID\"");

        #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
        fprintf(fout,", \"Thread\"");
        #endif

        if (CornerNodes.usedElements()==0) exit(__LINE__,__FILE__,"Error: CornerNodes are not allocated");
        CornerNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);

        #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
        if (CenterNodes.usedElements()==0) exit(__LINE__,__FILE__,"Error: CenterNodes are not allocated");
        CenterNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);
        #endif

        rootTree->block->PrintVariableList(fout);
      }

      if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_) {
        fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",meshNodesNumber,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_);
      }
      else if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_) {
        fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",4*meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_);
      }
      else exit(__LINE__,__FILE__,"Error: unknown value");

      #else
      fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\"");

      if (PrintMeshData==true) {
        fprintf(fout,", \"Maximum Refinment Level\", \"Temp_ID\"");

        #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
        fprintf(fout,", \"Thread\"");
        #endif

        if (CornerNodes.usedElements()==0) exit(__LINE__,__FILE__,"Error: CornerNodes are not allocated");
        CornerNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);

        #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
        if (CenterNodes.usedElements()==0) exit(__LINE__,__FILE__,"Error: CenterNodes are not allocated");
        CenterNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);
        #endif

        rootTree->block->PrintVariableList(fout);
      }

      if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_) {
        fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",meshNodesNumber,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);
      }
      else if (ConnectivityListMode==_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_) {
        fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",8*meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);
      }
      else exit(__LINE__,__FILE__,"Error: unknown value");
#endif
    }

    //print the node's list
    resetNodeProcessedFlag(false);

    outputMeshTECPLOT_BlockCornerNode_BlockConnectivityList(rootTree,fout,true,PrintMeshData,DataSetNumber,ConnectivityListMode);

    //print the connectivity list 
    #if _MESH_DIMENSION_ != 1
    outputMeshTECPLOT_BlockCornerNode_BlockConnectivityList(rootTree,fout,false,PrintMeshData,DataSetNumber,ConnectivityListMode);
    #endif

    //close the output file 
    if (ThisThread==0) fclose(fout);
  }


/*
  //print the mesh when blocks are not allocated
  void outputMeshTECPLOTinternalNoBlocks_PrintNodes(FILE* fout,cTreeNodeAMR<cBlockAMR> *startNode,bool PrintNodeCoordinates,bool PrintConnectivityList) {
    int isubBlock,jsubBlock,ksubBlock; ///,isubBlockMax,jsubBlockMax,ksubBlockMax;
    int iNode,jNode,kNode,idim,nnode; //,nBasicBlockNodes;
    double  xNode[3];

    static long int nGlobalNodeNumber=-1;

    if (startNode==rootTree) nGlobalNodeNumber=1;

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //the node is on the bottom of the tree -> printf the node
      //get the limits for the number of the sub-blocks in the block
      #if _MESH_DIMENSION_ == 1
      static const int isubBlockMax=_BLOCK_CELLS_X_;
      static const int jsubBlockMax=1,ksubBlockMax=1;
      static const int nBasicBlockNodes=2;
      static const int nodeOutputOrder_I[2]={0,1};
      static const int nodeOutputOrder_J[2]={0,0};
      static const int nodeOutputOrder_K[2]={0,0};

      #elif _MESH_DIMENSION_ == 2
      static const int isubBlockMax=_BLOCK_CELLS_X_;
      static const int jsubBlockMax=_BLOCK_CELLS_Y_;
      static const int ksubBlockMax=1;
      static const int nBasicBlockNodes=4;
      static const int nodeOutputOrder_I[4]={0,1,1,0};
      static const int nodeOutputOrder_J[4]={0,0,1,1};
      static const int nodeOutputOrder_K[4]={0,0,0,0};
      #else
      static const int isubBlockMax=_BLOCK_CELLS_X_;
      static const int jsubBlockMax=_BLOCK_CELLS_Y_;
      static const int ksubBlockMax=_BLOCK_CELLS_Z_;
      static const int nBasicBlockNodes=8;
      static const int nodeOutputOrder_I[8]={0,1,1,0,0,1,1,0};
      static const int nodeOutputOrder_J[8]={0,0,1,1,0,0,1,1};
      static const int nodeOutputOrder_K[8]={0,0,0,0,1,1,1,1};
      #endif


      for (ksubBlock=0;ksubBlock<ksubBlockMax;ksubBlock++) for (jsubBlock=0;jsubBlock<jsubBlockMax;jsubBlock++) for (isubBlock=0;isubBlock<isubBlockMax;isubBlock++) {
        for (nnode=0;nnode<nBasicBlockNodes;nnode++) {
          iNode=nodeOutputOrder_I[nnode]+isubBlock;
          jNode=nodeOutputOrder_J[nnode]+jsubBlock;
          kNode=nodeOutputOrder_K[nnode]+ksubBlock;

          startNode->GetCornerNodePosition(xNode,iNode,jNode,kNode);

          if (PrintNodeCoordinates==true) {
            for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(fout,"%e  ",xNode[idim]);
            fprintf(fout,"\n");
          } else if (PrintConnectivityList) {
            fprintf(fout,"%ld ",nGlobalNodeNumber);
            ++nGlobalNodeNumber;
            if (nnode==nBasicBlockNodes-1) fprintf(fout,"\n");
          }
        }

      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) outputMeshTECPLOTinternalNoBlocks_PrintNodes(fout,startNode->downNode[nDownNode],PrintNodeCoordinates,PrintConnectivityList);
    }
  }


  void outputMeshTECPLOTinternalNoBlocks(const char *fname) {
    FILE *fout;

    //  *DiagnospticMessageStream the number of the elements of the mesh
    if (meshModifiedFlag_CountMeshElements==true) countMeshElements(rootTree,0);

    fout=fopen(fname,"w");

    #if _MESH_DIMENSION_ == 1
    fprintf(fout,"VARIABLES=\"X\"");

    exit(__LINE__,__FILE__,"not implemented");

    if (PrintMeshData==true) {
      CornerNodes.elementStackList[0][0]->PrintVariableList(fout);

      #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
      CenterNodes.elementStackList[0][0]->PrintVariableList(fout);
      #endif

      rootTree->block->PrintVariableList(fout);
    }

    fprintf(fout,"\n");

    #elif _MESH_DIMENSION_ == 2
    fprintf(fout,"VARIABLES=\"X\", \"Y\"");
    fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",4*meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_);
    #else
    fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\"");
    fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",8*meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);
    #endif


    //print nodes coordinates
    outputMeshTECPLOTinternalNoBlocks_PrintNodes(fout,rootTree,true,false);

    //print the connectivity list
    if (_MESH_DIMENSION_>1) outputMeshTECPLOTinternalNoBlocks_PrintNodes(fout,rootTree,false,true);

    fclose(fout);
  }

*/

//#define _AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_    0
//#define _AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_ 1


  void outputMeshTECPLOT(const char *fname) {
    outputMeshTECPLOTinternal(fname,false,-1,_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_);
  }

  void outputMeshDataTECPLOT(const char *fname,int DataSetnumber=-1) {
    outputMeshTECPLOTinternal(fname,true,DataSetnumber,_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_);
  }

  //==============================================================================
  //save and read the mesh from a binary file 
  void saveCornerBlockNodes(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) { 
    static long int countingNumber;
    long nnode,nBlockNodes,nDownNode; 
    cCornerNode **blockNode;

    if (startNode==rootTree) countingNumber=0; 

    #if _MESH_DIMENSION_ == 1
    nBlockNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #else
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    for (nnode=0,blockNode=startNode->block->cornerNodes;nnode<nBlockNodes;nnode++,blockNode++) if ((*blockNode)->nodeDescriptor.nodeProcessedFlag==_AMR_FALSE_) {
      (*blockNode)->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_;
      (*blockNode)->nodeDescriptor.nodeno=countingNumber++;

      fwrite(*blockNode,sizeof(cCornerNode),1,fout);
    }

    //save the nodes in the downNodes 
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) saveCornerBlockNodes(startNode->downNode[nDownNode],fout);
  } 


  void readCornerBlockNodes(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) {
    long nnode;
    cCornerNode *blockNode;

    //clean the 'nodes' stack
    CornerNodes.clear(); 

    //fill up the nodes in a sequential order
    for (nnode=0;nnode<meshNodesNumber;nnode++) {
      blockNode=CornerNodes.newElement();  

      if (fread(blockNode,sizeof(cCornerNode),1,fout)!=1) exit(__LINE__,"Error while reading the mesh nodes list file"); 
    }
  }
    
  void saveTreeStructure(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) {
    //long int nnode,nBlockNodes,nDownNode,countingNumber,i;
    //cCornerNode **blockNode;

    long int nDownNode,countingNumber,i;

    /*
    #if _MESH_DIMENSION_ == 1
    nBlockNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #else
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif
    */

    //save the tree node's pointers 
    countingNumber=treeNodes.GetEntryCountingNumber(startNode->upNode);
    fwrite(&countingNumber,sizeof(long int),1,fout);

    //save downNode[:]
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->downNode[nDownNode]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }

    //save the node ID
    fwrite(&startNode->AMRnodeID,sizeof(cAMRnodeID),1,fout);


    //save the node neighbors
#if _MESH_DIMENSION_ == 1
//    cTreeNodeAMR *neibNodeFace[2];
    for (i=0;i<2;i++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNodeFace[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }
#elif _MESH_DIMENSION_ == 2
    for (i=0;i<4*2;i++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNodeFace[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }

    for (i=0;i<4;i++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNodeCorner[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }
#elif _MESH_DIMENSION_ == 3
//    cTreeNodeAMR *neibNodeFace[6*4],*neibNodeCorner[8],*neibNodeEdge[12*2];
    for (i=0;i<6*4;i++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNodeFace[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }

    for (i=0;i<8;i++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNodeCorner[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }

    for (i=0;i<12*2;i++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNodeEdge[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }
#endif

/*
    for (i=0;i<2*_MESH_DIMENSION_;i++) {
//      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNode[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);

      exit(__LINE__,__FILE__,"writing of neibNodes is node implemeted");
   }
   */

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    fwrite(&startNode->Temp_ID,sizeof(long int),1,fout);
    #endif    

    //save other tree nodes' parameters: xmin,xmax,RefinmentLevel,nodeDescriptor,Thread;
    fwrite(startNode->xmin,_MESH_DIMENSION_*sizeof(double),1,fout);
    fwrite(startNode->xmax,_MESH_DIMENSION_*sizeof(double),1,fout);
    fwrite(&startNode->RefinmentLevel,sizeof(int),1,fout);
    fwrite(&startNode->nodeDescriptor,sizeof(startNode->nodeDescriptor),1,fout);

#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
    fwrite(&startNode->Thread,sizeof(int),1,fout);
#endif


#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    //save the list of nodes that intersect surfaces of the computational domain
    countingNumber=treeNodes.GetEntryCountingNumber(startNode->DomainSurfaceBoundaryList_Next);
    fwrite(&countingNumber,sizeof(long int),1,fout);

    countingNumber=treeNodes.GetEntryCountingNumber(startNode->DomainSurfaceBoundaryList_Prev);
    fwrite(&countingNumber,sizeof(long int),1,fout);

    //save the list of the boundary surface's descriptors
    int InternalSurfaceBoundaryPresentFlag=(startNode->InternalBoundaryDescriptorList==NULL) ? false : true;
    cInternalBoundaryConditionsDescriptor *SurfaceDescriptor;
    list<cInternalBoundaryConditionsDescriptor>::iterator ptrDescriptorList,ptrDescriptorListEnd;

    fwrite(&InternalSurfaceBoundaryPresentFlag,sizeof(int),1,fout);

    if (InternalSurfaceBoundaryPresentFlag==true) {
      bool found;

      for (SurfaceDescriptor=startNode->InternalBoundaryDescriptorList;SurfaceDescriptor!=NULL;SurfaceDescriptor=SurfaceDescriptor->nextInternalBCelement) {
        ptrDescriptorListEnd=InternalBoundaryList.end();
        found=false;

        //save the marker for the beginig of the new record
        fwrite(&InternalSurfaceBoundaryPresentFlag,sizeof(int),1,fout);

        //find the counting number for the intersected surface
        for (countingNumber=0,ptrDescriptorList=InternalBoundaryList.begin();ptrDescriptorList!=ptrDescriptorListEnd;countingNumber++,ptrDescriptorList++) {
          if (ptrDescriptorList->BoundaryElement==SurfaceDescriptor->BoundaryElement) {
            found=true;
            fwrite(SurfaceDescriptor,sizeof(cInternalBoundaryConditionsDescriptor),1,fout);
            fwrite(&countingNumber,sizeof(long int),1,fout);
            break;
          }
        }

        if (found==false) exit(__LINE__,__FILE__,"Error: cannot find the interesected surface data");
      }

      InternalSurfaceBoundaryPresentFlag=false;
      fwrite(&InternalSurfaceBoundaryPresentFlag,sizeof(int),1,fout);
    }



#endif


    /*
    //save the blocks' information
    bool BlockExist=(startNode->block!=NULL) ? true : false;

    fwrite(&BlockExist,sizeof(bool),1,fout);

    if (BlockExist==true) {
      cBlockAMR *block=startNode->block;

      #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
      fwrite(&block->Temp_ID,sizeof(long int),1,fout);
      #endif 

      fwrite(&block->blockDescriptor,sizeof(typename cBlockAMR::cBlockDescriptor),1,fout);

      for (nnode=0,blockNode=startNode->block->cornerNodes;nnode<nBlockNodes;nnode++,blockNode++) {
        countingNumber=(*blockNode)->nodeDescriptor.nodeno;
        fwrite(&countingNumber,sizeof(long int),1,fout);
      }
    }
    */

    //save the downNodes
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) saveTreeStructure(startNode->downNode[nDownNode],fout);
  }    

  void readTreeStructure(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) {
    //long int nnode,nBlockNodes,nDownNode,countingNumber,i;
    long int nDownNode,countingNumber,i;

    /*
    #if _MESH_DIMENSION_ == 1
    nBlockNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #else
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif
    */

    startNode->cleanDataBuffer();

    //read the tree node's pointers
    fread(&countingNumber,sizeof(long int),1,fout);
    startNode->upNode=treeNodes.GetEntryPointer(countingNumber);

    //read downNode[:]
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->downNode[nDownNode]=treeNodes.GetEntryPointer(countingNumber);      

      if (startNode->downNode[nDownNode]!=NULL) startNode->downNode[nDownNode]->upNode=startNode;
    }

    //read the node ID
    fread(&startNode->AMRnodeID,sizeof(cAMRnodeID),1,fout);

    //read the node neighbors
#if _MESH_DIMENSION_ == 1
//    cTreeNodeAMR *neibNodeFace[2];
    for (i=0;i<2;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->neibNodeFace[i]=treeNodes.GetEntryPointer(countingNumber);
    }
#elif _MESH_DIMENSION_ == 2
    for (i=0;i<4*2;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->neibNodeFace[i]=treeNodes.GetEntryPointer(countingNumber);
    }

    for (i=0;i<4;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->neibNodeCorner[i]=treeNodes.GetEntryPointer(countingNumber);
    }
#elif _MESH_DIMENSION_ == 3
//    cTreeNodeAMR *neibNodeFace[6*4],*neibNodeCorner[8],*neibNodeEdge[12*2];
    for (i=0;i<6*4;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->neibNodeFace[i]=treeNodes.GetEntryPointer(countingNumber);
    }

    for (i=0;i<8;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->neibNodeCorner[i]=treeNodes.GetEntryPointer(countingNumber);
    }

    for (i=0;i<12*2;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->neibNodeEdge[i]=treeNodes.GetEntryPointer(countingNumber);
    }

#endif

    /*
    for (i=0;i<2*_MESH_DIMENSION_;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
//      startNode->neibNode[i]=treeNodes.GetEntryPointer(countingNumber);
      exit(__LINE__,__FILE__,"reading of neibNodes is not implemented");
    }
    */

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    fread(&startNode->Temp_ID,sizeof(long int),1,fout);
    #endif

    //read other tree nodes' parameters: xmin,xmax,RefinmentLevel,nodeDescriptor,Thread;
    fread(startNode->xmin,_MESH_DIMENSION_*sizeof(double),1,fout);
    fread(startNode->xmax,_MESH_DIMENSION_*sizeof(double),1,fout);
    fread(&startNode->RefinmentLevel,sizeof(int),1,fout);
    fread(&startNode->nodeDescriptor,sizeof(startNode->nodeDescriptor),1,fout);

#if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
    fread(&startNode->Thread,sizeof(int),1,fout);
#endif

    //save the list of nodes that intersect surfaces of the computational domain
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    fread(&countingNumber,sizeof(long int),1,fout);
    startNode->DomainSurfaceBoundaryList_Next=treeNodes.GetEntryPointer(countingNumber);

    fread(&countingNumber,sizeof(long int),1,fout);
    startNode->DomainSurfaceBoundaryList_Prev=treeNodes.GetEntryPointer(countingNumber);

    //read the list of the boundary surface's descriptors
    int InternalSurfaceBoundaryPresentFlag;
    cInternalBoundaryConditionsDescriptor *SurfaceDescriptor;
    list<cInternalBoundaryConditionsDescriptor>::iterator ptrDescriptorList,ptrDescriptorListEnd;

    fread(&InternalSurfaceBoundaryPresentFlag,sizeof(int),1,fout);

    if (InternalSurfaceBoundaryPresentFlag==true) {
//      int i;

      ptrDescriptorListEnd=InternalBoundaryList.end();
      fread(&InternalSurfaceBoundaryPresentFlag,sizeof(int),1,fout);

      while (InternalSurfaceBoundaryPresentFlag==true) {
        SurfaceDescriptor=InternalBoundaryDescriptors.newElement();

        fread(SurfaceDescriptor,sizeof(cInternalBoundaryConditionsDescriptor),1,fout);
        fread(&countingNumber,sizeof(long int),1,fout);

        for (i=0,SurfaceDescriptor->BoundaryElement=NULL,ptrDescriptorList=InternalBoundaryList.begin();ptrDescriptorList!=ptrDescriptorListEnd;i++,ptrDescriptorList++) {
          if (i==countingNumber) {
            SurfaceDescriptor->BoundaryElement=ptrDescriptorList->BoundaryElement;
            break;
          }
        }

        if (SurfaceDescriptor->BoundaryElement==NULL) exit(__LINE__,__FILE__,"Error: cannot find approproate surface element");

        SurfaceDescriptor->nextInternalBCelement=startNode->InternalBoundaryDescriptorList;
        startNode->InternalBoundaryDescriptorList=SurfaceDescriptor;

        fread(&InternalSurfaceBoundaryPresentFlag,sizeof(int),1,fout);
      }
    }
#endif

    /*
    //save the blocks' information
    bool BlockExist;

    fread(&BlockExist,sizeof(bool),1,fout);

    if (BlockExist==true) {
      cBlockAMR *block;

      startNode->block=blocks.newElement(); 
      block=startNode->block; 

      #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
      fread(&block->Temp_ID,sizeof(long int),1,fout);
      #endif

      fread(&block->blockDescriptor,sizeof(typename cBlockAMR::cBlockDescriptor),1,fout);

      for (nnode=0;nnode<nBlockNodes;nnode++) {
        fread(&countingNumber,sizeof(long int),1,fout);
        startNode->block->cornerNodes[nnode]=CornerNodes.GetEntryPointer(countingNumber);  
      }
    }
    else startNode->block=NULL;
    */

    //read the downNodes
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) readTreeStructure(startNode->downNode[nDownNode],fout); 
  }



  void saveMeshFile(const char *MeshFileName=NULL) {
    FILE *fout=NULL;
    long int countingNumber,ThisThread=0;
    char fname[STRING_LENGTH];
 
    //get the name of the mesh file 
    if (MeshName[0]=='\0') generateMeshName();
//    if (meshModifiedFlag==true) generateMeshSignature();   

    if (MeshFileName!=NULL) sprintf(fname,"%s",MeshFileName);
    else getMeshName(fname);
 


    //the 'global' parameters of the mesh is saved in the file by the root processor only 
    if (ThisThread==0) {
      fout=fopen(fname,"w");

      //save the mesh parameters: in the beginig of the mesh file keep the name of the mesh!
      fwrite(MeshName,sizeof(char),STRING_LENGTH,fout) ;
//      fwrite(&MeshSignature,sizeof(unsigned long),1,fout);


      //the global definition of the domain limits 
      fwrite(_MESH_AMR_XMAX_,sizeof(double),3,fout);
      fwrite(_MESH_AMR_XMIN_,sizeof(double),3,fout); 

      //the global parameters of the mesh  
      fwrite(&EPS,sizeof(double),1,fout);
//      fwrite(&meshModifiedFlag,sizeof(bool),1,fout);

      fwrite(&meshNodesNumber,sizeof(long int),1,fout);
      fwrite(&meshBlocksNumber,sizeof(long int),1,fout);
      fwrite(&meshMaximumRefinmentLevel,sizeof(int),1,fout);

      fwrite(&xGlobalMin,sizeof(double),_MESH_DIMENSION_,fout);
      fwrite(&xGlobalMax,sizeof(double),_MESH_DIMENSION_,fout);

      fwrite(&AllowBlockAllocation,sizeof(bool),1,fout);

      //save allocation for the stacks
      fwrite("AMR-MESH-FILE-MARKER:MEMORY_ALLOCATION",sizeof(char),STRING_LENGTH,fout);

//      CornerNodes.saveAllocationParameters(fout);
//      blocks.saveAllocationParameters(fout);
      treeNodes.saveAllocationParameters(fout);

      //save the pointer to the rootTree
      countingNumber=treeNodes.GetEntryCountingNumber(rootTree);
      fwrite(&countingNumber,sizeof(long int),1,fout); 


      //save the tree structure
      fwrite("AMR-MESH-FILE-MARKER:MESH-TREE",sizeof(char),STRING_LENGTH,fout);
      saveTreeStructure(rootTree,fout);


      //save the list of nodes that intersect surfaces of the computational domain
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
      fwrite("AMR-MESH-FILE-MARKER:BOUNDARY_NODES",sizeof(char),STRING_LENGTH,fout);
      countingNumber=treeNodes.GetEntryCountingNumber(DomainSurfaceBoundaryList);
      fwrite(&countingNumber,sizeof(long int),1,fout);

      //save the distribution of the cut faces
      if (CutCell::nBoundaryTriangleFaces!=0) {

        struct {
          void SaveBlock(cTreeNodeAMR<cBlockAMR> *node,FILE *fout) {
            if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
              //the node is at the bottom of the tree -> save the list of the NASTRAN cut cells;
              CutCell::cTriangleFaceDescriptor *t;
              bool NextTriangleFaceExist=false;
              std::ptrdiff_t TriangleFaceID;


              for (t=node->FirstTriangleCutFace;t!=NULL;t=t->next) {
                NextTriangleFaceExist=true;
                fwrite(&NextTriangleFaceExist,sizeof(bool),1,fout);

                TriangleFaceID=t->TriangleFace-CutCell::BoundaryTriangleFaces;
                fwrite(&TriangleFaceID,sizeof(std::ptrdiff_t),1,fout);
              }

              NextTriangleFaceExist=false;
              fwrite(&NextTriangleFaceExist,sizeof(bool),1,fout);


            }
            else for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) SaveBlock(node->downNode[nDownNode],fout);
          }
        } SaveBlockCutCellDistribution;

        SaveBlockCutCellDistribution.SaveBlock(rootTree,fout);
      }
#endif




    /*
    //save corner block nodes
    fwrite("AMR-MESH-FILE-MARKER:CORNER-BLOCK-NODES",sizeof(char),STRING_LENGTH,fout);
    resetNodeProcessedFlag();
    saveCornerBlockNodes(rootTree,fout); 
    */

      //save the 'end' marker
      fwrite("AMR-MESH-FILE-MARKER:END",sizeof(char),STRING_LENGTH,fout);
      fclose(fout);
    }
  }
   
  void readMeshFile(const char *MeshFileName) {
    FILE *fout=NULL;
    long int countingNumber;
    char marker[STRING_LENGTH];

    //the 'global' parameters of the mesh is saved in the file by the root processor only
//    if (ThisThread==0) {
      fout=fopen(MeshFileName,"r");

      //save the mesh parameters: in the beginig of the mesh file keep the name of the mesh!
      fread(MeshName,sizeof(char),STRING_LENGTH,fout) ;
//      fread(&MeshSignature,sizeof(unsigned long),1,fout);


      //the global definition of the domain limits
      fread(_MESH_AMR_XMAX_,sizeof(double),3,fout);
      fread(_MESH_AMR_XMIN_,sizeof(double),3,fout);

      //the global parameters of the mesh
      fread(&EPS,sizeof(double),1,fout);
//      fread(&meshModifiedFlag,sizeof(bool),1,fout);

      fread(&meshNodesNumber,sizeof(long int),1,fout);
      fread(&meshBlocksNumber,sizeof(long int),1,fout);
      fread(&meshMaximumRefinmentLevel,sizeof(int),1,fout);

      fread(&xGlobalMin,sizeof(double),_MESH_DIMENSION_,fout);
      fread(&xGlobalMax,sizeof(double),_MESH_DIMENSION_,fout);

      fread(&AllowBlockAllocation,sizeof(bool),1,fout);

      //read allocation for the stacks
      fread(marker,sizeof(char),STRING_LENGTH,fout);
      if (strcmp("AMR-MESH-FILE-MARKER:MEMORY_ALLOCATION",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong");

//      CornerNodes.readAllocationParameters(fout);
//      blocks.readAllocationParameters(fout);
      treeNodes.readAllocationParameters(fout);

      //read the pointer to the rootTree
      fread(&countingNumber,sizeof(long int),1,fout);
      rootTree=treeNodes.GetEntryPointer(countingNumber);
//    }

    //read the tree structure
    fread(marker,sizeof(char),STRING_LENGTH,fout);
    if (strcmp("AMR-MESH-FILE-MARKER:MESH-TREE",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong"); 
    readTreeStructure(rootTree,fout);


    //read the list of nodes that intersect surfaces of the computational domain
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    fread(marker,sizeof(char),STRING_LENGTH,fout);
    if (strcmp("AMR-MESH-FILE-MARKER:BOUNDARY_NODES",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong");
    fread(&countingNumber,sizeof(long int),1,fout);
    DomainSurfaceBoundaryList=treeNodes.GetEntryPointer(countingNumber);


    //read the distribution of the cut faces
    if (CutCell::nBoundaryTriangleFaces!=0) {

      struct {
        void LoadBlock(cTreeNodeAMR<cBlockAMR> *node,FILE *fout) {
          if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
            //the node is at the bottom of the tree -> save the list of the NASTRAN cut cells;
            CutCell::cTriangleFaceDescriptor *t;
            bool NextTriangleFaceExist=false;
            std::ptrdiff_t TriangleFaceID;

            fread(&NextTriangleFaceExist,sizeof(bool),1,fout);

            while (NextTriangleFaceExist==true) {
              //set up the description of the cut face
              CutCell::cTriangleFaceDescriptor *t=CutCell::BoundaryTriangleFaceDescriptor.newElement();

              fread(&TriangleFaceID,sizeof(std::ptrdiff_t),1,fout);
              t->TriangleFace=CutCell::BoundaryTriangleFaces+TriangleFaceID;

              t->prev=NULL,t->next=node->FirstTriangleCutFace;
              if (node->FirstTriangleCutFace!=NULL) node->FirstTriangleCutFace->prev=t;
              node->FirstTriangleCutFace=t;

              fread(&NextTriangleFaceExist,sizeof(bool),1,fout);
            }
          }
          else for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) LoadBlock(node->downNode[nDownNode],fout);
        }
      } LoadBlockCutCellDistribution;

      LoadBlockCutCellDistribution.LoadBlock(rootTree,fout);
    }
#endif

    /*
    //read corner block nodes
    fread(marker,sizeof(char),STRING_LENGTH,fout);
    if (strcmp("AMR-MESH-FILE-MARKER:CORNER-BLOCK-NODES",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong");
    readCornerBlockNodes(rootTree,fout);
    */

    //read the 'end' marker 
    fread(marker,sizeof(char),STRING_LENGTH,fout);
    if (strcmp("AMR-MESH-FILE-MARKER:END",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong");
    fclose(fout);

    //determine the range of a block neighbor's relution level
    SetNodeNeibResolutionLevelLimit();
  }

  //create the memory allocation report
  void countTreeNodes(cTreeNodeAMR<cBlockAMR> *startNode,long int *Counter,int level) {
    int nDownNode;

    if (level>_MAX_REFINMENT_LEVEL_) exit(__LINE__,__FILE__,"The level of refinment exeeds the value of _MAX_REFINMENT_LEVEL_");

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      Counter[level]++;
    }
    else for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) countTreeNodes(startNode->downNode[nDownNode],Counter,level+1);
  }
 
  //Determine resolution limits of the neib nodes
  void SetNodeNeibResolutionLevelLimit(cTreeNodeAMR<cBlockAMR> *startNode=NULL) {
    int nDownNode;

    if (startNode==NULL) startNode=rootTree;

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) startNode->SetNeibRefinmentLevelLimits();
    else for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) SetNodeNeibResolutionLevelLimit(startNode->downNode[nDownNode]);
  }

  //distribute the cut-faces among the blocks
  void DistributeBoundaryCutBlocks(cTreeNodeAMR<cBlockAMR> *startNode=NULL) {
    int nDownNode;
    static bool initflag=false;


    if ((startNode==NULL)||(startNode==rootBlock)) {
      if (initflag!=false) exit(__LINE__,__FILE__,"Error: repeated call of 'DistributeBoundaryCutBlocks'");
      initflag=true;
    }

    if (startNode==NULL) {
      startNode=rootTree;
    }

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      //check intersection of the cut-faces with the block
      //if the list is already allocated -> remove it

      if (startNode->FirstTriangleCutFace!=NULL) {
        CutCell::cTriangleFaceDescriptor *tnext,*t=startNode->FirstTriangleCutFace;

        while (t!=NULL) {
          tnext=t->next;
          CutCell::BoundaryTriangleFaceDescriptor.deleteElement(t);
          t=tnext;
        }
      }

      //check all boundary faces for intersection with the block
      int nface;

      double xmin[3],xmax[3];

      xmin[0]=startNode->xmin[0]-(startNode->xmax[0]-startNode->xmin[0])/double(_BLOCK_CELLS_X_)*double(_GHOST_CELLS_X_);
      xmin[1]=startNode->xmin[1]-(startNode->xmax[1]-startNode->xmin[1])/double(_BLOCK_CELLS_Y_)*double(_GHOST_CELLS_Y_);
      xmin[2]=startNode->xmin[2]-(startNode->xmax[2]-startNode->xmin[2])/double(_BLOCK_CELLS_Z_)*double(_GHOST_CELLS_Z_);

      xmax[0]=startNode->xmax[0]+(startNode->xmax[0]-startNode->xmin[0])/double(_BLOCK_CELLS_X_)*double(_GHOST_CELLS_X_);
      xmax[1]=startNode->xmax[1]+(startNode->xmax[1]-startNode->xmin[1])/double(_BLOCK_CELLS_Y_)*double(_GHOST_CELLS_Y_);
      xmax[2]=startNode->xmax[2]+(startNode->xmax[2]-startNode->xmin[2])/double(_BLOCK_CELLS_Z_)*double(_GHOST_CELLS_Z_);

      for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
        if (CutCell::BoundaryTriangleFaces[nface].BlockIntersection(xmin,xmax,EPS)==true) {
          //the block is intersected by the face
          CutCell::cTriangleFaceDescriptor *t=CutCell::BoundaryTriangleFaceDescriptor.newElement();

          t->TriangleFace=CutCell::BoundaryTriangleFaces+nface;

          t->prev=NULL,t->next=startNode->FirstTriangleCutFace;
          if (startNode->FirstTriangleCutFace!=NULL) startNode->FirstTriangleCutFace->prev=t;
          startNode->FirstTriangleCutFace=t;
        }
      }
    }
    else for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) DistributeBoundaryCutBlocks(startNode->downNode[nDownNode]);

  }


  void DistributeBoundaryCutBlocks(cTreeNodeAMR<cBlockAMR> *startNode,CutCell::cTriangleFaceDescriptor *BoundaryFaces) {
    //check intersection of the cut-faces with the block
    //if the list is already allocated -> remove it
    CutCell::cTriangleFaceDescriptor *tnext,*t=startNode->FirstTriangleCutFace;

    if (startNode->FirstTriangleCutFace!=NULL) {
      while (t!=NULL) {
        tnext=t->next;
        CutCell::BoundaryTriangleFaceDescriptor.deleteElement(t);
        t=tnext;
      }
    }

    //check all boundary faces for intersection with the block
    //check all boundary faces for intersection with the block
        double xmin[3],xmax[3];

        xmin[0]=startNode->xmin[0]-(startNode->xmax[0]-startNode->xmin[0])/double(_BLOCK_CELLS_X_)*double(_GHOST_CELLS_X_);
        xmin[1]=startNode->xmin[1]-(startNode->xmax[1]-startNode->xmin[1])/double(_BLOCK_CELLS_Y_)*double(_GHOST_CELLS_Y_);
        xmin[2]=startNode->xmin[2]-(startNode->xmax[2]-startNode->xmin[2])/double(_BLOCK_CELLS_Z_)*double(_GHOST_CELLS_Z_);

        xmax[0]=startNode->xmax[0]+(startNode->xmax[0]-startNode->xmin[0])/double(_BLOCK_CELLS_X_)*double(_GHOST_CELLS_X_);
        xmax[1]=startNode->xmax[1]+(startNode->xmax[1]-startNode->xmin[1])/double(_BLOCK_CELLS_Y_)*double(_GHOST_CELLS_Y_);
        xmax[2]=startNode->xmax[2]+(startNode->xmax[2]-startNode->xmin[2])/double(_BLOCK_CELLS_Z_)*double(_GHOST_CELLS_Z_);

    for (t=BoundaryFaces;t!=NULL;t=t->next) {
      if (t->TriangleFace->BlockIntersection(xmin,xmax,EPS)==true) {
        //the block is intersected by the face
        CutCell::cTriangleFaceDescriptor *startNodeBoundaryDescriptor=CutCell::BoundaryTriangleFaceDescriptor.newElement();

        startNodeBoundaryDescriptor->TriangleFace=t->TriangleFace;

        startNodeBoundaryDescriptor->prev=NULL,startNodeBoundaryDescriptor->next=startNode->FirstTriangleCutFace;
        if (startNode->FirstTriangleCutFace!=NULL) startNode->FirstTriangleCutFace->prev=startNodeBoundaryDescriptor;
        startNode->FirstTriangleCutFace=startNodeBoundaryDescriptor;



/*
        //DEBUG BEGIN
        t->TriangleFace->BlockIntersection(startNode->xmin,startNode->xmax,EPS);
        t->TriangleFace->BlockIntersection(startNode->xmin,startNode->xmax,EPS);
        t->TriangleFace->BlockIntersection(startNode->xmin,startNode->xmax,EPS);
        t->TriangleFace->BlockIntersection(startNode->xmin,startNode->xmax,EPS);
        t->TriangleFace->BlockIntersection(startNode->xmin,startNode->xmax,EPS);
        //DEBUG END
*/





      }
    }
  }

  void DistributeBoundaryCutBlocks(cTreeNodeAMR<cBlockAMR> *startNode,CutCell::cTriangleFace *BoundaryFaces, int nBoundaryFaces) {
    //check intersection of the cut-faces with the block
    //if the list is already allocated -> remove it
    CutCell::cTriangleFaceDescriptor *tnext,*t=startNode->FirstTriangleCutFace;

    if (startNode->FirstTriangleCutFace!=NULL) {
      while (t!=NULL) {
        tnext=t->next;
        CutCell::BoundaryTriangleFaceDescriptor.deleteElement(t);
        t=tnext;
      }
    }

    //check all boundary faces for intersection with the block
    double xmin[3],xmax[3];

    xmin[0]=startNode->xmin[0]-(startNode->xmax[0]-startNode->xmin[0])/double(_BLOCK_CELLS_X_)*double(_GHOST_CELLS_X_);
    xmin[1]=startNode->xmin[1]-(startNode->xmax[1]-startNode->xmin[1])/double(_BLOCK_CELLS_Y_)*double(_GHOST_CELLS_Y_);
    xmin[2]=startNode->xmin[2]-(startNode->xmax[2]-startNode->xmin[2])/double(_BLOCK_CELLS_Z_)*double(_GHOST_CELLS_Z_);

    xmax[0]=startNode->xmax[0]+(startNode->xmax[0]-startNode->xmin[0])/double(_BLOCK_CELLS_X_)*double(_GHOST_CELLS_X_);
    xmax[1]=startNode->xmax[1]+(startNode->xmax[1]-startNode->xmin[1])/double(_BLOCK_CELLS_Y_)*double(_GHOST_CELLS_Y_);
    xmax[2]=startNode->xmax[2]+(startNode->xmax[2]-startNode->xmin[2])/double(_BLOCK_CELLS_Z_)*double(_GHOST_CELLS_Z_);

    for (int nface=0;nface<nBoundaryFaces;nface++)  {
      if (BoundaryFaces[nface].BlockIntersection(xmin,xmax,EPS)==true) {
        //the block is intersected by the face
        CutCell::cTriangleFaceDescriptor *startNodeBoundaryDescriptor=CutCell::BoundaryTriangleFaceDescriptor.newElement();

        startNodeBoundaryDescriptor->TriangleFace=BoundaryFaces+nface;

        startNodeBoundaryDescriptor->prev=NULL,startNodeBoundaryDescriptor->next=startNode->FirstTriangleCutFace;
        if (startNode->FirstTriangleCutFace!=NULL) startNode->FirstTriangleCutFace->prev=startNodeBoundaryDescriptor;
        startNode->FirstTriangleCutFace=startNodeBoundaryDescriptor;

      }
    }
  }

  /*
  //generate the list of the blocks of the mesh
  void createMeshBlockList(cTreeNodeAMR<cBlockAMR>  *startNode=NULL) {
    int iDownNode,jDownNode,kDownNode;  

    #if _MESH_DIMENSION_ == 1
    static const int iDownNodeMax=1,jDownNodeMax=0,kDownNodeMax=0; 
    #elif _MESH_DIMENSION_ == 2
    static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=0;
    #elif _MESH_DIMENSION_ == 3
    static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=1;
    #endif

    if (startNode==NULL) blockList=NULL,startNode=rootTree; 

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      startNode->block->nextBlock=blockList;
      blockList=(void*)startNode->block;
    }   
    else {
      for (kDownNode=0;kDownNode<=kDownNodeMax;kDownNode++) for (jDownNode=0;jDownNode<=jDownNodeMax;jDownNode++) for (iDownNode=0;iDownNode<=iDownNodeMax;iDownNode++) {
        createMeshBlockList(startNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)]);
      } 
    }
  }
  */


  //determine the memory allocated by the mesh
  void memoryAllocationReport() {
    long int Counter[1+_MAX_REFINMENT_LEVEL_];
    int level,thread;

    CMPI_channel pipe(100000);

    for (level=0;level<=_MAX_REFINMENT_LEVEL_;level++) Counter[level]=0;

    //count the number of tree nodes for each refinment level 
    countTreeNodes(rootTree,Counter,0);

    if (ThisThread==0) {
/*      *DiagnospticMessageStream << "$PREFIX:Mesh blocks report:" << std::endl;
      for (level=0;level<=_MAX_REFINMENT_LEVEL_;level++) *DiagnospticMessageStream << "$PREFIX:refinment level=" << level << ", blocks=" << Counter[level] << std::endl;

      *DiagnospticMessageStream << "$PREFIX:Memory usage:" << std::endl;
      *DiagnospticMessageStream << "$PREFIX:Thread \t Tree \t Blocks \t Nodes" << std::endl;*/


      fprintf(DiagnospticMessageStream,"$PREFIX:Mesh blocks report:\n");
      for (level=0;level<=_MAX_REFINMENT_LEVEL_;level++) fprintf(DiagnospticMessageStream,"$PREFIX:refinment level=%i, blocks=%ld\n",level,Counter[level]);

      fprintf(DiagnospticMessageStream,"$PREFIX:Memory usage:\n");
      fprintf(DiagnospticMessageStream,"$PREFIX:Thread \t Tree \t Blocks \t Nodes\n");

      pipe.openRecvAll();
    }
    else pipe.openSend(0);

    long int TreeNodesAllocation,BlocksAllocation,CornerNodesAllocation;

    TreeNodesAllocation=treeNodes.getAllocatedMemory();
    BlocksAllocation=blocks.getAllocatedMemory();
    CornerNodesAllocation=CornerNodes.getAllocatedMemory();

    for (thread=0;thread<nTotalThreads;thread++) {
      if (thread!=0) {
        if (ThisThread==0) {
          pipe.recv(TreeNodesAllocation,thread);
          pipe.recv(BlocksAllocation,thread);
          pipe.recv(CornerNodesAllocation,thread);
        }
        else if (ThisThread==thread) {
          pipe.send(TreeNodesAllocation);
          pipe.send(BlocksAllocation);
          pipe.send(CornerNodesAllocation);
        }
      }

//      if (ThisThread==0) *DiagnospticMessageStream << "$PREFIX:" << thread << "\t" << TreeNodesAllocation/1.0E6 << "MB\t" << BlocksAllocation/1.0E6 << "MB\t" << CornerNodesAllocation/1.0E6 << "MB" << std::endl;
      if (ThisThread==0) fprintf(DiagnospticMessageStream,"$PREFIX:%i\t%eMB\t%eMB\t%eMB\n",thread,TreeNodesAllocation/1.0E6,BlocksAllocation/1.0E6,CornerNodesAllocation/1.0E6);
    }



    rusage ResourceUsage;
    if (getrusage(RUSAGE_SELF,&ResourceUsage)!=-1)  {
      long int ru_maxrss,ru_ixrss,ru_idrss;

      ru_maxrss=ResourceUsage.ru_maxrss;
      ru_ixrss=ResourceUsage.ru_ixrss;
      ru_idrss=ResourceUsage.ru_idrss;

//      if (ThisThread==0) *DiagnospticMessageStream << "$PREFIX:Thread \t \"maximum resident set size\"\t \"integral shared memory size\"\t  \"integral unshared data size\"" <<std::endl;
      if (ThisThread==0) fprintf(DiagnospticMessageStream,"$PREFIX:Thread \t \"maximum resident set size\"\t \"integral shared memory size\"\t  \"integral unshared data size\"\n");

      for (thread=0;thread<nTotalThreads;thread++) {
        if (thread!=0) {
          if (ThisThread==0) {
            pipe.recv(ru_maxrss,thread);
            pipe.recv(ru_ixrss,thread);
            pipe.recv(ru_idrss,thread);
          }
          else if (ThisThread==thread) {
            pipe.send(ru_maxrss);
            pipe.send(ru_ixrss);
            pipe.send(ru_idrss);
          }
        }

//        if (ThisThread==0) *DiagnospticMessageStream << "$PREFIX:" << thread << "\t" << ru_maxrss/1.0E6 << "MB\t" << ru_ixrss/1.0E6 << "MB\t" << ru_idrss/1.0E6 << "MB" << std::endl;
        if (ThisThread==0) fprintf(DiagnospticMessageStream,"$PREFIX:%i\t%eMB\t%eMB\t%eMB\n",thread,ru_maxrss/1.0E6,ru_ixrss/1.0E6,ru_idrss/1.0E6);
      }
    } 

    if (ThisThread==0) pipe.closeRecvAll();
    else pipe.closeSend();
  }
    
  //determine weather 'startNode' is attached to the boundary of the bounding box (external boundary)
  #define _EXTERNAL_BOUNDARY_BLOCK_  true
  #define _INTERNAL_BLOCK_           false

  bool ExternalBoundaryBlock(cTreeNodeAMR<cBlockAMR>* startNode,bool *BoundaryFaceFlag=NULL) {
    int i,j,k,counter,nface;
    double x,*xmin,*xmax;
    bool found=false,ExternalFaceNodes[8]={false,false,false,false,false,false,false,false};

    #if _MESH_DIMENSION_ == 1
    static const int nFaceNodes=1;
    static const int FaceNodeMap[2*_MESH_DIMENSION_][nFaceNodes]={ {0}, {1}};
    static const int iMax=2,jMax=1,kMax=1;
    #elif _MESH_DIMENSION_ == 2
    static const int nFaceNodes=2;
    static const int FaceNodeMap[2*_MESH_DIMENSION_][nFaceNodes]={ {0,2}, {1,3}, {0,1}, {2,3}};
    static const int iMax=2,jMax=2,kMax=1;
    #elif _MESH_DIMENSION_ == 3
    static const int nFaceNodes=4;
    static const int FaceNodeMap[2*_MESH_DIMENSION_][nFaceNodes]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};
    static const int iMax=2,jMax=2,kMax=2;
    #endif

    xmin=startNode->xmin;
    xmax=startNode->xmax;

    for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
      x=(i==0) ? xmin[0] : xmax[0];
      if ((fabs(x-xGlobalMin[0])<EPS)||(fabs(x-xGlobalMax[0])<EPS)) {
        ExternalFaceNodes[i+2*(j+2*k)]=true;
        found=true;
        continue;
      }

      if (_MESH_DIMENSION_>1) {
        x=(j==0) ? xmin[1] : xmax[1];
        if ((fabs(x-xGlobalMin[1])<EPS)||(fabs(x-xGlobalMax[1])<EPS)) {
          ExternalFaceNodes[i+2*(j+2*k)]=true;
          found=true;
          continue;
        }
      }

      if (_MESH_DIMENSION_>2) {
        x=(k==0) ? xmin[2] : xmax[2];
        if ((fabs(x-xGlobalMin[2])<EPS)||(fabs(x-xGlobalMax[2])<EPS)) {
          ExternalFaceNodes[i+2*(j+2*k)]=true;
          found=true;
          continue;
        }
      }
    }

    if (found==false) return _INTERNAL_BLOCK_;

    //determine which face is connected to the external boundary
    if (BoundaryFaceFlag==NULL) return _EXTERNAL_BOUNDARY_BLOCK_;

    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) {
      for (i=0,counter=0;i<nFaceNodes;i++) if (ExternalFaceNodes[ FaceNodeMap[nface][i] ]==true) counter++;

      BoundaryFaceFlag[nface]=(counter==nFaceNodes) ? true : false;
    }

    return _EXTERNAL_BOUNDARY_BLOCK_;
  }

  //calcualte the coorfinate frame on a face of a block
  void GetBlockFaceCoordinateFrame_3D(double *x0,double *e0,double *e1,int nface,cTreeNodeAMR<cBlockAMR>* startNode) {
    int nd0,nd1,nd2,idim;
    double *xmin,*xmax;

    if (nface>=6) exit(__LINE__,__FILE__"Error: nface is out of the range");

    static const int nFaceNodes=4;
    static const int FaceNodeMap[2*3][nFaceNodes]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};  //the array dimension: FaceNodeMap[2*_MESH_DIMENSION_][nFaceNodes]
    static const int CornerNodesMap[8][3]={ {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},   {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};

    nd0=FaceNodeMap[nface][0];
    nd1=FaceNodeMap[nface][1];
    nd2=FaceNodeMap[nface][2];

    xmin=startNode->xmin;
    xmax=startNode->xmax;

    for (idim=0;idim<_MESH_DIMENSION_;idim++) {
      x0[idim]=(CornerNodesMap[nd0][idim]==1) ? xmax[idim] : xmin[idim];

      e0[idim]=((CornerNodesMap[nd1][idim]==1) ? xmax[idim] : xmin[idim]) - x0[idim];
      e1[idim]=((CornerNodesMap[nd2][idim]==1) ? xmax[idim] : xmin[idim]) - x0[idim];
    }
  }

  //collect into a list all nodes that "belongs" to this processor
  void InitDomainBoundaryLayer(cTreeNodeAMR<cBlockAMR>* startNode) {
    long int i,neibThread,nCounter=0;
    cTreeNodeAMR<cBlockAMR> *node,*prev,*next,*neibNode,*tNode;
    bool found;

    for (i=0;i<nTotalThreads;i++) DomainBoundaryLayerNodesList[i]=NULL;
    node=ParallelNodesDistributionList[ThisThread];

    while (node!=NULL) {


      //Add the corner nodes to the list
#if _MESH_DIMENSION_ == 2
#define _AMR_InitDomainBoundaryLayer_ADD_CORNER_NODES_
#elif _MESH_DIMENSION_ == 3
#define _AMR_InitDomainBoundaryLayer_ADD_CORNER_NODES_
#endif

#ifdef _AMR_InitDomainBoundaryLayer_ADD_CORNER_NODES_
      for (i=0;i<(1<<_MESH_DIMENSION_);i++) if ((neibNode=node->neibNodeCorner[i])!=NULL) if ((neibThread=neibNode->Thread)!=ThisThread) {
        //check if the node is not in the list already
        found=false,tNode=DomainBoundaryLayerNodesList[neibThread];

        while (tNode!=NULL) {
          if (tNode==neibNode) {
            found=true;
            break;
          }

          tNode=tNode->nextNodeThisThread;
        }

        if (found==true) continue;

        nCounter++;
        prev=neibNode->prevNodeThisThread;
        next=neibNode->nextNodeThisThread;

        if (prev!=NULL) prev->nextNodeThisThread=next;
        if (next!=NULL) next->prevNodeThisThread=prev;

        neibNode->prevNodeThisThread=NULL;
        neibNode->nextNodeThisThread=DomainBoundaryLayerNodesList[neibThread];

        if (DomainBoundaryLayerNodesList[neibThread]!=NULL) DomainBoundaryLayerNodesList[neibThread]->prevNodeThisThread=neibNode;
        DomainBoundaryLayerNodesList[neibThread]=neibNode;
      }
#endif


      //Add faces to the list
//#if _MESH_DIMENSION_ == 2
//#define _AMR_InitDomainBoundaryLayer_ADD_FACES_
//#elif _MESH_DIMENSION_ == 3
//#define _AMR_InitDomainBoundaryLayer_ADD_FACES_
//#endif

//#ifdef _AMR_InitDomainBoundaryLayer_ADD_FACES_
      for (i=0;i<_MESH_DIMENSION_*(1<<_MESH_DIMENSION_);i++) if ((neibNode=node->neibNodeFace[i])!=NULL) if ((neibThread=neibNode->Thread)!=ThisThread) {
        //check if the node is not in the list already
        found=false,tNode=DomainBoundaryLayerNodesList[neibThread];

        while (tNode!=NULL) {
          if (tNode==neibNode) {
            found=true;
            break;
          }

          tNode=tNode->nextNodeThisThread;
        }

        if (found==true) continue;

        nCounter++;
        prev=neibNode->prevNodeThisThread;
        next=neibNode->nextNodeThisThread;

        if (prev!=NULL) prev->nextNodeThisThread=next;
        if (next!=NULL) next->prevNodeThisThread=prev;

        neibNode->prevNodeThisThread=NULL;
        neibNode->nextNodeThisThread=DomainBoundaryLayerNodesList[neibThread];

        if (DomainBoundaryLayerNodesList[neibThread]!=NULL) DomainBoundaryLayerNodesList[neibThread]->prevNodeThisThread=neibNode;
        DomainBoundaryLayerNodesList[neibThread]=neibNode;
      }

//#endif


      //Add edges to the list
#if _MESH_DIMENSION_ == 3
//cTreeNodeAMR *neibNodeFace[6*4],*neibNodeCorner[8],*neibNodeEdge[12*2];
      for (i=0;i<12*2;i++) if ((neibNode=node->neibNodeEdge[i])!=NULL) if ((neibThread=neibNode->Thread)!=ThisThread) {
        //check if the node is not in the list already
        found=false,tNode=DomainBoundaryLayerNodesList[neibThread];

        while (tNode!=NULL) {
          if (tNode==neibNode) {
            found=true;
            break;
          }

          tNode=tNode->nextNodeThisThread;
        }

        if (found==true) continue;

        nCounter++;
        prev=neibNode->prevNodeThisThread;
        next=neibNode->nextNodeThisThread;

        if (prev!=NULL) prev->nextNodeThisThread=next;
        if (next!=NULL) next->prevNodeThisThread=prev;

        neibNode->prevNodeThisThread=NULL;
        neibNode->nextNodeThisThread=DomainBoundaryLayerNodesList[neibThread];

        if (DomainBoundaryLayerNodesList[neibThread]!=NULL) DomainBoundaryLayerNodesList[neibThread]->prevNodeThisThread=neibNode;
        DomainBoundaryLayerNodesList[neibThread]=neibNode;
      }
#endif



      node=node->nextNodeThisThread;
    }


    long int *buffer=new long int [nTotalThreads];
    int thread;

    buffer[0]=nCounter;

    long int bufferRecv[nTotalThreads];
    MPI_Gather(buffer,1,MPI_LONG,bufferRecv,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    memcpy(buffer,bufferRecv,nTotalThreads*sizeof(long int));

    if (ThisThread==0) {
      fprintf(DiagnospticMessageStream,"$PREFIX:The length of the domain's boundary nodes list:\n$PREFIX:Thread\tThe number of the domain's boundary nodes\n");

      for (thread=0;thread<nTotalThreads;thread++) fprintf(DiagnospticMessageStream,"$PREFIX:%i\t%ld\n",thread,buffer[thread]);
    }

    delete [] buffer;


    //test the real length of the list (check the consisntency of the lists)
    long int nTest;

    for (thread=0;thread<nTotalThreads;thread++) {
      nTest=0;
      node=DomainBoundaryLayerNodesList[thread];

      if (node!=NULL) if (node->prevNodeThisThread!=NULL) exit(__LINE__,__FILE__,"Error: the list is not consistent");

      while (node!=NULL) {
        nTest++;
        if (nTest>nCounter) { //the list is locked somewhere. find the lock up of the list
          cTreeNodeAMR<cBlockAMR> *startSearch=DomainBoundaryLayerNodesList[thread];

          for (i=0;i<nCounter;i++) {
            node=startSearch->nextNodeThisThread;

            while (node!=NULL) {
              if (node==startSearch) {
                fprintf(DiagnospticMessageStream,"$PREFIX:Error: have found a repearting node in the list. nodeid=%ld\n",node->Temp_ID);
              }

              node=node->nextNodeThisThread;
            }

            startSearch=startSearch->nextNodeThisThread;
          }

          exit(__LINE__,__FILE__,"Error: the list is locked-up somewhere");
        }

        node=node->nextNodeThisThread;
      }
    }

  }


  /*
  void InitLocalNodesList(cTreeNodeAMR<cBlockAMR>* startNode=NULL) {
    int thread;

    if (startNode==NULL) {
      startNode=rootTree;

      for (thread=0;thread<nTotalThreads;thread++) ParallelNodesDistributionList[thread]=NULL;
    }

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //add the node to the list
      thread=startNode->Thread;

      startNode->prevNodeThisThread=NULL,startNode->nextNodeThisThread=ParallelNodesDistributionList[thread];
      if (ParallelNodesDistributionList[thread]!=NULL) ParallelNodesDistributionList[thread]->prevNodeThisThread=startNode;
      ParallelNodesDistributionList[thread]=startNode;
    }
    else for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InitLocalNodesList(startNode->downNode[nDownNode]);


#if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
    if (startNode==rootTree) InitDomainBoundaryLayer(rootTree);
#elif _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__GHOST_CELLS_
    //do nothing
#else
    exit(__LINE__,__FILE__,"Error: wrong option");
#endif


  }

*/

  void InitCellMeasure_ResetToZero(cTreeNodeAMR<cBlockAMR>* startNode) {
#if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    int i,j,k;
    long int nd;
    cCenterNode *centerNode;

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InitCellMeasure_ResetToZero(startNode->downNode[nDownNode]);
    }
    else if (startNode->block!=NULL) {
#if _MESH_DIMENSION_ == 1
      static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
      static const int jMin=0,jMax=1;
      static const int kMin=0,kMax=1;
#elif _MESH_DIMENSION_ == 2
       static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
       static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
       static const int kMin=0,kMax=1;
#elif _MESH_DIMENSION_ == 3
       static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
       static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
       static const int kMin=-_GHOST_CELLS_Z_,kMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
#else
       exit(__LINE__,__FILE__,"Error: unknown option");
#endif

       for (k=kMin;k<kMax;k++) for (j=jMin;j<jMax;j++) for (i=iMin;i<iMax;i++) {
         nd=getCenterNodeLocalNumber(i,j,k);

         if ((centerNode=startNode->block->GetCenterNode(nd))!=NULL) centerNode->Measure=-1.0;
       }
    }
#endif
  }


  //save/load the binary file containing the center node measures
  void SaveCenterNodeMeasure(cTreeNodeAMR<cBlockAMR> *node=NULL) {
    static CMPI_channel pipe;
    static FILE *fout=NULL;

    if (node==NULL) {
      unsigned long Signature=getMeshSignature();

      node=rootTree;
      pipe.init(1000000);

      if (ThisThread==0) {
        pipe.openRecvAll();

        //open the file
        char fname[400];
        sprintf(fname,"amr.sig=0x%lx.CenterNodeMeasure.bin",Signature);

        fout=fopen(fname,"w");
      }
      else pipe.openSend(0);
    }

    //loop through all points
    //create the list of the points
    //perform the interpolation loop
    int i,j,k,nd;
    cCenterNode *CenterNode;
    double Measure;


    #if _MESH_DIMENSION_ == 3
    const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
    const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
    const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
    #else
    exit(__LINE__,__FILE__,"Error: not implemented");
    #endif



    if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      if ((ThisThread==0)||(node->Thread==ThisThread)) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
        //get cell's measure
        if (node->Thread==ThisThread) {
          Measure=-1.0;

          //locate the cell
          if (node->block!=NULL) {
            nd=getCenterNodeLocalNumber(i,j,k);
            if ((CenterNode=node->block->GetCenterNode(nd))!=NULL) Measure=CenterNode->Measure;
          }

          if (node->Thread!=0) pipe.send(Measure);
        }
        else {
          pipe.recv(Measure,node->Thread);
        }

        //save cells' measure
        if (ThisThread==0) fwrite(&Measure,sizeof(double),1,fout);

      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (node->downNode[nDownNode]!=NULL) SaveCenterNodeMeasure(node->downNode[nDownNode]);
    }

    if (node==rootTree) {
      if (ThisThread==0) {
        pipe.closeRecvAll();
        fclose(fout);
      }
      else pipe.closeSend();

      pipe.remove();
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
  }

  bool LoadCenterNodeMeasure(const char *CellMeasureFilePath, cTreeNodeAMR<cBlockAMR> *node=NULL) {
    static FILE *fData=NULL;

    if (node==NULL) {
      node=rootTree;

      //open the file
      char fname[600];

      if (CellMeasureFilePath==NULL) {
        sprintf(fname,"amr.sig=0x%lx.CenterNodeMeasure.bin",getMeshSignature());
      }
      else {
        sprintf(fname,"%s/amr.sig=0x%lx.CenterNodeMeasure.bin",CellMeasureFilePath,getMeshSignature());
      }

      fData=fopen(fname,"r");

      if (fData==NULL) {
        //the file does not exists -> the measure has to be calculated directly
        if (ThisThread==0) printf("$PREFIX: the calculated volume data file %s is not found\n",fname);
        return false;
      }
      else {
        if (ThisThread==0) printf("$PREFIX: the calculated volume data file %s is found: loading....  ",fname);
      }
    }

    //loop through all points
    //create the list of the points
    //perform the interpolation loop
    int i,j,k,nd;
    cCenterNode *CenterNode;
    double Measure;


    #if _MESH_DIMENSION_ == 3
    const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
    const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
    const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
    #else
    exit(__LINE__,__FILE__,"Error: not implemented");
    #endif


    if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      if (node->block!=NULL) {
        for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
          fread(&Measure,sizeof(double),1,fData);

          nd=getCenterNodeLocalNumber(i,j,k);
          if ((CenterNode=node->block->GetCenterNode(nd))!=NULL) CenterNode->Measure=Measure;
        }
      }
      else {
        //skip the measure data for the block
        fseek(fData,(iMax-iMin+1)*(jMax-jMin+1)*(kMax-kMin+1)*sizeof(double),SEEK_CUR);
      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (node->downNode[nDownNode]!=NULL) LoadCenterNodeMeasure(CellMeasureFilePath,node->downNode[nDownNode]);
    }

    if (node==rootTree) {
      fclose(fData);

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      if (ThisThread==0) printf("done. \n");
    }

    return true;
  }


  //save/load the center node associated data. nSaveVariables==-1 -> save the whole associated data vector
  void SaveCenterNodeAssociatedData(const char *fNameBase, int *SaveVariableOffset, int nSaveVariables,cTreeNodeAMR<cBlockAMR>  *node=NULL) {
    static CMPI_channel pipe;
    static FILE *fout=NULL;

    if (node==NULL) {
      unsigned long Signature;

      Signature=getMeshSignature();
      node=rootTree;
      pipe.init(1000000);

      if (ThisThread==0) {
        pipe.openRecvAll();

        //open the file
        char fname[400];
        sprintf(fname,"amr.sig=0x%lx.f=%s.CenterNodeAssociatedData.bin",Signature,fNameBase);

        fout=fopen(fname,"w");

        //save the number of the variables
        fwrite(&nSaveVariables,sizeof(int),1,fout);
      }
      else pipe.openSend(0);
    }

    //loop through all points
    //create the list of the points
    //perform the interpolation loop
    int i,j,k,nd;
    cCenterNode *CenterNode;
    char *offset;

    #if _MESH_DIMENSION_ == 3
    const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
    const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
    const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
    #else
    exit(__LINE__,__FILE__,"Error: not implemented");
    #endif

    char SendCellFlag;

    if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      if ((ThisThread==0)||(node->Thread==ThisThread)) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
        SendCellFlag=true;

        //determine whether the cell data neede to be saved
        if (node->Thread==ThisThread) {
          //locate the cell
          if (node->block==NULL) SendCellFlag=false,offset=NULL;

          nd=getCenterNodeLocalNumber(i,j,k);

          if (SendCellFlag==true) {
            if ((CenterNode=node->block->GetCenterNode(nd))!=NULL) offset=CenterNode->GetAssociatedDataBufferPointer();
            else SendCellFlag=false,offset=NULL;
          }

          if (node->Thread!=0) pipe.send(SendCellFlag);
        }
        else {
          pipe.recv(SendCellFlag,node->Thread);
        }

        //save the cell data saving flag
        if (ThisThread==0) fwrite(&SendCellFlag,sizeof(char),1,fout);

        //save the cell data
        if (SendCellFlag==true) {
          if (node->Thread==ThisThread) {
            if (node->Thread==0) {
              if (nSaveVariables==-1) fwrite(offset,sizeof(char),CenterNode->AssociatedDataLength(),fout);
              else for (int nvar=0;nvar<nSaveVariables;nvar++) {
                fwrite(offset+SaveVariableOffset[nvar],sizeof(double),1,fout);
              }
            }
            else {
              pipe.send(offset,CenterNode->AssociatedDataLength());
            }
          }
          else {
            char data[CenterNode->AssociatedDataLength()];

            //recieve the data
            pipe.recv(data,CenterNode->AssociatedDataLength(),node->Thread);

            //save the data
            if (nSaveVariables==-1) fwrite(data,sizeof(char),CenterNode->AssociatedDataLength(),fout);
            else for (int nvar=0;nvar<nSaveVariables;nvar++) {
              fwrite(data+SaveVariableOffset[nvar],sizeof(double),1,fout);
            }

          }
        }
      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (node->downNode[nDownNode]!=NULL) SaveCenterNodeAssociatedData(fNameBase,SaveVariableOffset,nSaveVariables,node->downNode[nDownNode]);
    }

    if (node==rootTree) {
      if (ThisThread==0) {
        pipe.closeRecvAll();
        fclose(fout);
      }
      else pipe.closeSend();

      pipe.remove();
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
  }

  void SaveCenterNodeAssociatedData(const char *fNameBase) {
    SaveCenterNodeAssociatedData(fNameBase,NULL,-1,NULL);
  }

  void LoadCenterNodeAssociatedData(const char *fNameBase, const char* DataFilePath, int *LoadVariableOffset, int nLoadVariables,cTreeNodeAMR<cBlockAMR>  *node=NULL) {
    static FILE *fData=NULL;
    static int SavedDataLength;

    if (node==NULL) {
      node=rootTree;

      //open the file
      char fname[600];

      if (DataFilePath==NULL) {
        sprintf(fname,"amr.sig=0x%lx.f=%s.CenterNodeAssociatedData.bin",getMeshSignature(),fNameBase);
      }
      else {
        sprintf(fname,"%s/amr.sig=0x%lx.f=%s.CenterNodeAssociatedData.bin",DataFilePath,getMeshSignature(),fNameBase);
      }

      fData=fopen(fname,"r");

      if (fData==NULL) {
        //the file does not exists -> the measure has to be calculated directly
        exit(__LINE__,__FILE__,"Error: cannot open file");
      }

      //compare the number of the saved variables with 'nLoadVariables'
      int nSaveVariables;

      fread(&nSaveVariables,sizeof(int),1,fData);
      if (nSaveVariables!=nLoadVariables) exit(__LINE__,__FILE__,"Error: the number of load varibles is different from the number of saved variables");

      //determine the size of the saved data
      cCenterNode cell;
      SavedDataLength=(nSaveVariables!=-1) ? nSaveVariables*sizeof(double) : cell.AssociatedDataLength();
    }


    //loop through all points
    //create the list of the points
    //perform the interpolation loop
    int i,j,k,nd;
    cCenterNode *CenterNode;
    char *offset;

    #if _MESH_DIMENSION_ == 3
    const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
    const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
    const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
    #else
    exit(__LINE__,__FILE__,"Error: not implemented");
    #endif

    char savedLoadCellFlag,data[SavedDataLength];

    if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      if (node->block==NULL) {
         //the block belongs to a other processor -> skip all data
        for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)  {
          fread(&savedLoadCellFlag,sizeof(char),1,fData);

          if (savedLoadCellFlag==true) {
            //the cell data is saved -> skip it
            fseek(fData,SavedDataLength,SEEK_CUR);
          }
        }

      }
      else for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
        fread(&savedLoadCellFlag,sizeof(char),1,fData);

        if (savedLoadCellFlag==true) {
          //determine whether the cell data needed to be read
          //locate the cell
          nd=getCenterNodeLocalNumber(i,j,k);

          if ((CenterNode=node->block->GetCenterNode(nd))!=NULL) {
            offset=CenterNode->GetAssociatedDataBufferPointer();

            //read center cells' associated data
            fread(data,sizeof(char),SavedDataLength,fData);

            //copy the data
            if (nLoadVariables==-1) memcpy(offset,data,CenterNode->AssociatedDataLength());
            else for (int nvar=0;nvar<nLoadVariables;nvar++) {
              memcpy(offset+LoadVariableOffset[nvar],data+nvar*sizeof(double),sizeof(double));
            }

          }
          else fseek(fData,SavedDataLength,SEEK_CUR);
        }

      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (node->downNode[nDownNode]!=NULL) LoadCenterNodeAssociatedData(fNameBase,DataFilePath,LoadVariableOffset,nLoadVariables,node->downNode[nDownNode]);
    }

    if (node==rootTree) {
      fclose(fData);
    }
  }

  void LoadCenterNodeAssociatedData(const char *fNameBase,const char* DataFilePath) {
    LoadCenterNodeAssociatedData(fNameBase,DataFilePath,NULL,-1,NULL);
  }

  void LoadCenterNodeAssociatedData(const char *fNameBase) {
    LoadCenterNodeAssociatedData(fNameBase,NULL);
  }

  //check whether the associated data file exists
  bool AssociatedDataFileExists(const char *fNameBase,const char* DataFilePath) {
    FILE *fData=NULL;
    char fname[400];

    if (DataFilePath==NULL) {
      sprintf(fname,"amr.sig=0x%lx.f=%s.CenterNodeAssociatedData.bin",getMeshSignature(),fNameBase);
    }
    else {
      sprintf(fname,"%s/amr.sig=0x%lx.f=%s.CenterNodeAssociatedData.bin",DataFilePath,getMeshSignature(),fNameBase);
    }


    fData=fopen(fname,"r");

    if (fData!=NULL) {
      fclose(fData);
      return true;
    }

    return false;
  }

  bool AssociatedDataFileExists(const char *fNameBase) {
    return AssociatedDataFileExists(fNameBase,NULL);
  }

  //evaluate the total volume of the computational domain
  double GetTotalDomainMeasure(cTreeNodeAMR<cBlockAMR>* startNode=NULL) {
#if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    int i,j,k;
    long int nd;
    cCenterNode *centerNode;
    double res=0.0;

    if (startNode==NULL) startNode=rootTree;

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) res+=GetTotalDomainMeasure(startNode->downNode[nDownNode]);
    }
    else if ((startNode->block!=NULL)&&(startNode->Thread==ThisThread)) {
#if _MESH_DIMENSION_ == 1
      static const int iMin=0,iMax=_BLOCK_CELLS_X_;
      static const int jMin=0,jMax=1;
      static const int kMin=0,kMax=1;
#elif _MESH_DIMENSION_ == 2
       static const int iMin=0,iMax=_BLOCK_CELLS_X_;
       static const int jMin=0,jMax=_BLOCK_CELLS_Y_;
       static const int kMin=0,kMax=1;
#elif _MESH_DIMENSION_ == 3
       static const int iMin=0,iMax=_BLOCK_CELLS_X_;
       static const int jMin=0,jMax=_BLOCK_CELLS_Y_;
       static const int kMin=0,kMax=_BLOCK_CELLS_Z_;
#else
       exit(__LINE__,__FILE__,"Error: unknown option");
#endif

       for (k=kMin;k<kMax;k++) for (j=jMin;j<jMax;j++) for (i=iMin;i<iMax;i++) {
         nd=getCenterNodeLocalNumber(i,j,k);

         if ((centerNode=startNode->block->GetCenterNode(nd))!=NULL) res+=centerNode->Measure;
       }
    }
#endif

    //collect measure from all processors
    if (startNode==rootTree) {
      double TotalMeasure;

      MPI_Allreduce(&res,&TotalMeasure,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      res=TotalMeasure;
    }

    return res;
  }


  void InitCellMeasureBlock(cTreeNodeAMR<cBlockAMR>* startNode) {
       double *xNodeMin,*xNodeMax,Measure,xCellMin[3],xCellMax[3],dx=0.0,dy=0.0,dz=0.0,xTotalMin[3]={0.0,0.0,0.0},xTotalMax[3]={0.0,0.0,0.0};
       long int i,j,k,nd;
       cCenterNode *centerNode;
       int IntersectionCode=_AMR_BLOCK_INSIDE_DOMAIN_;


#if _MESH_DIMENSION_ == 1
       static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
       static const int jMin=0,jMax=1;
       static const int kMin=0,kMax=1;
#elif _MESH_DIMENSION_ == 2
       static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
       static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
       static const int kMin=0,kMax=1;
#elif _MESH_DIMENSION_ == 3
       static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
       static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
       static const int kMin=-_GHOST_CELLS_Z_,kMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
#else
       exit(__LINE__,__FILE__,"Error: unknown option");
#endif


       if (startNode->block==NULL)  exit(__LINE__,__FILE__,"Error: the block is not allocated");

       xNodeMin=startNode->xmin,xNodeMax=startNode->xmax;

       dx=(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_;
       Measure=dx;
       xTotalMin[0]=xNodeMin[0]-dx*_GHOST_CELLS_X_,xTotalMax[0]=xNodeMax[0]+dx*_GHOST_CELLS_X_;

       if (_MESH_DIMENSION_>1) {
         dy=(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_;
         Measure*=dy;
         xTotalMin[1]=xNodeMin[1]-dy*_GHOST_CELLS_Y_,xTotalMax[1]=xNodeMax[1]+dy*_GHOST_CELLS_Y_;
       }

       if (_MESH_DIMENSION_>2) {
         dz=(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_;
         Measure*=dz;
         xTotalMin[2]=xNodeMin[2]-dz*_GHOST_CELLS_Z_,xTotalMax[2]=xNodeMax[2]+dz*_GHOST_CELLS_Z_;
       }

#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
       if (startNode->InternalBoundaryDescriptorList!=NULL) {
         int IntersectionCodeTemp=-1;

         for (cInternalBoundaryConditionsDescriptor *bc=startNode->InternalBoundaryDescriptorList;bc!=NULL;bc=bc->nextInternalBCelement) {
            if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) {
              IntersectionCodeTemp=((cInternalSphericalData*)(bc->BoundaryElement))->BlockIntersection(xTotalMin,xTotalMax,EPS);
            }
            else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_CIRCLE_) {
              IntersectionCodeTemp=((cInternalCircleData*)(bc->BoundaryElement))->BlockIntersection(xTotalMin,xTotalMax,EPS);
            }
            else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_1D_SPHERE_) {
              IntersectionCodeTemp=((cInternalSphere1DData*)(bc->BoundaryElement))->BlockIntersection(xTotalMin,xTotalMax,EPS);
            }
            else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_) {
              IntersectionCodeTemp=((cInternalRotationBodyData*)(bc->BoundaryElement))->BlockIntersection(xTotalMin,xTotalMax,EPS);
            }
            else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_) {
              IntersectionCodeTemp=((cInternalNastranSurfaceData*)(bc->BoundaryElement))->BlockIntersection(xTotalMin,xTotalMax,EPS);
            }
            else exit(__LINE__,__FILE__,"Error: unknown boundary type");

            if (IntersectionCodeTemp==_AMR_BLOCK_OUTSIDE_DOMAIN_) IntersectionCode=_AMR_BLOCK_OUTSIDE_DOMAIN_;
            else if ((IntersectionCodeTemp==_AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_)&&(IntersectionCode!=_AMR_BLOCK_OUTSIDE_DOMAIN_)) IntersectionCode=_AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
         }


       }
#endif

       if (IntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) Measure=-1.0;

       for (k=kMin;k<kMax;k++) {
         if (_MESH_DIMENSION_>2) xCellMin[2]=xNodeMin[2]+k*dz,xCellMax[2]=xNodeMin[2]+(k+1)*dz;

         for (j=jMin;j<jMax;j++) {
           if (_MESH_DIMENSION_>1) xCellMin[1]=xNodeMin[1]+j*dy,xCellMax[1]=xNodeMin[1]+(j+1)*dy;

           for (i=iMin;i<iMax;i++) {
             xCellMin[0]=xNodeMin[0]+i*dx,xCellMax[0]=xNodeMin[0]+(i+1)*dx;
             nd=getCenterNodeLocalNumber(i,j,k);

             if ((centerNode=startNode->block->GetCenterNode(nd))!=NULL) if (centerNode->Measure<=0.0) {
               if ((IntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_)||(IntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_)) {
                 //set the cell measure with the accounting for possible symmetry cases
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
                 centerNode->Measure=Measure;
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
                 centerNode->Measure=4.0*Pi/3.0*fabs(pow(fabs(xCellMax[0]),3)-pow(fabs(xCellMin[0]),3));
#else
                 exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif
               }

               else if (startNode->FirstTriangleCutFace!=NULL) {


                 //TEST
                 double *a=&centerNode->Measure;

                 if ((a==(double*)0x14f374938) || (a==(double*)0x14ab5f938) || (centerNode->Temp_ID==9444) )  {
                   double test=0;

                   test+=1;
                 }

                 //END TEST


                 //the block is cutted by a list of triangulat faces
                 centerNode->Measure=GetRemainedBlockVolume(xCellMin,xCellMax,EPS,1.0E-2,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,startNode->FirstTriangleCutFace);
               }



               else {
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
                 int BoundarySurfaceCounter;
                 cInternalBoundaryConditionsDescriptor *bc;

                 if (centerNode->Measure<=0.0) {
                   for (BoundarySurfaceCounter=0,bc=startNode->InternalBoundaryDescriptorList;bc!=NULL;bc=bc->nextInternalBCelement,BoundarySurfaceCounter++) {
                     int cellIntersectionStatus;

                     if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) {
                       centerNode->Measure=((cInternalSphericalData*)(bc->BoundaryElement))->GetRemainedBlockVolume(xCellMin,xCellMax,EPS,&cellIntersectionStatus);
                     }
                     else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_CIRCLE_) {
                       centerNode->Measure=((cInternalCircleData*)(bc->BoundaryElement))->GetRemainedBlockVolume(xCellMin,xCellMax,EPS,&cellIntersectionStatus);
                     }
                     else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_1D_SPHERE_) {
                       centerNode->Measure=((cInternalSphere1DData*)(bc->BoundaryElement))->GetRemainedBlockVolume(xCellMin,xCellMax,EPS,&cellIntersectionStatus);
                     }
                     else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_) {
                       centerNode->Measure=((cInternalRotationBodyData*)(bc->BoundaryElement))->GetRemainedBlockVolume(xCellMin,xCellMax,EPS,&cellIntersectionStatus);
                     }
                     else if (bc->BondaryType==_INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_) {
                       centerNode->Measure=((cInternalNastranSurfaceData*)(bc->BoundaryElement))->GetRemainedBlockVolume(xCellMin,xCellMax,EPS,&cellIntersectionStatus);
                     }
                     else exit(__LINE__,__FILE__,"Error: unknown boundary type");


                   }

                   if (BoundarySurfaceCounter!=1) exit(__LINE__,__FILE__,"Error: only one internal surface is permitted for a cell");
                 }

#else
                 exit(__LINE__,__FILE__,"Error: unknown option");
#endif
               }

             }


           }
         }
       }


  }



  void InitCellMeasure(const char* CellMeasureFilePath,cTreeNodeAMR<cBlockAMR>* startNode=NULL) {
   #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    if (startNode==NULL) {
      //try to load the saved file with the cell measures
      if (_AMR_READ_SAVE_CENTER_NODE_MEASURE__MODE_==_ON_AMR_MESH_) {
        if (ThisThread==0) printf("$PREFIX: check existance of the calculated volume data file.... \n");

        if (LoadCenterNodeMeasure(CellMeasureFilePath,NULL)==true) {
          return;
        }
        else if (ThisThread==0) printf("$PREFIX: the calculated volume data file is not found: generating....\n");
      }

      startNode=rootTree;
    }


    if (startNode==rootTree) InitCellMeasure_ResetToZero(rootTree);




//=============   DEBUG ==================
/*
    if (startNode!=NULL) if (startNode->Temp_ID==5096) {
  *DiagnospticMessageStream << __FILE__ << "@" << __LINE__ << std::endl;
}
*/
//=============  END DEBUG ===============


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  //allocate the Thread::FirstPBufferParticle array and distribute particles between processes
  #pragma omp parallel
  {
    #pragma omp single
    {
#endif

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InitCellMeasure(startNode->downNode[nDownNode]);
    }
    else if (startNode->block!=NULL) { //set up the cell's volume
      static long int nProcessedBlocks=0;
      static  time_t LastTimeValue=time(NULL);
      time_t TimeValue=time(NULL);

      #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      #pragma omp task default (none) firstprivate (startNode)
      #endif
      InitCellMeasureBlock(startNode);

      nProcessedBlocks++;

      if (TimeValue-LastTimeValue>120.0) {
        tm *ct=localtime(&TimeValue);

        LastTimeValue=TimeValue;
        printf("MESH: InitMeasure: Processed %li blocks (thread=%i; %i/%i %i:%i:%i)\n",nProcessedBlocks,ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
      }
    }
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    }}
#endif


    if (startNode==rootTree) {
      //try to save the file with the cell measures
      if (_AMR_READ_SAVE_CENTER_NODE_MEASURE__MODE_==_ON_AMR_MESH_) {
        if (ThisThread==0) printf("$PREFIX: saving the calculated volume data file.... ");
        SaveCenterNodeMeasure(NULL);
        if (ThisThread==0) printf("   done\n");
      }
    }
    #endif
  }


  void InitCellMeasure(cTreeNodeAMR<cBlockAMR>* startNode) {
    InitCellMeasure(NULL,startNode);
  }

  void InitCellMeasure() {
    InitCellMeasure(NULL,NULL);
  }

  //Calculate ID of an AMR node,  and find an AMR node by an ID
  /*
  class cAMRnodeID {
  public:
    int ResolutionLevel;
    unsigned char id[1+3*_MAX_REFINMENT_LEVEL_/8];

    cAMRnodeID() {
      ResolutionLevel=0;
      for (int i=0;i<1+3*_MAX_REFINMENT_LEVEL_/8;i++) id[i]=0;
    }

    bool operator == (cAMRnodeID ID) {
      if (ResolutionLevel!=ID.ResolutionLevel) return false;

      //compare the bit's fields
      int nbytes,i,nbits;

      nbits=3*ResolutionLevel;
      nbytes=nbits/8;
      nbits-=8*ResolutionLevel;

      //compare bytes
      for (i=0;i<nbytes;i++) if (id[i]!=ID.id[i]) return false;

      //compare bits
      unsigned char ComparisonMask=0;

//      for (i=0;i<nbits;i++) ComparisonMask+=(1<<i); substituted by the following
      for (i=0;i<nbits;i++) ComparisonMask=(unsigned char)(ComparisonMask|(1<<i));

      if ((id[nbytes]&ComparisonMask)!=(ID.id[nbytes]&ComparisonMask)) return false;

      return true;
    }

    bool operator != (cAMRnodeID ID) {
      return ((*this)==ID) ? false : true;
    }
  };
  */

  cTreeNodeAMR<cBlockAMR>* findAMRnodeWithID(cAMRnodeID node) {
    int Level,i,j,k,nDownNode;
    cTreeNodeAMR<cBlockAMR>* res=rootTree;
    unsigned char mask=node.id[0];
    int byteOffset=0,bitOffset=0;


    for (Level=0;Level<node.ResolutionLevel;Level++) {
      i=((mask&(1<<bitOffset))==0) ? 0 : 1;

      if (++bitOffset==8) {
        bitOffset=0,byteOffset++;
        mask=node.id[byteOffset];
      }

      j=((mask&(1<<bitOffset))==0) ? 0 : 1;

      if (++bitOffset==8) {
        bitOffset=0,byteOffset++;
        mask=node.id[byteOffset];
      }

      k=((mask&(1<<bitOffset))==0) ? 0 : 1;

      if (++bitOffset==8) {
        bitOffset=0,byteOffset++;
        mask=node.id[byteOffset];
      }

      nDownNode=i+2*(j+2*k);
      res=res->downNode[nDownNode];

      if (res==NULL) exit(__LINE__,__FILE__,"Error: the mesh tree is not consistent");
    }

    return res;
  }


  void GetAMRnodeID(cAMRnodeID& node,cTreeNodeAMR<cBlockAMR>* startNode) {
    int Level,i,j,k,nDownNode,nDownNodes;
    cTreeNodeAMR<cBlockAMR>* upNode;
    unsigned char mask=0;
    int byteOffset=0,bitOffset=0;

    node.ResolutionLevel=startNode->RefinmentLevel;

    bitOffset=3*node.ResolutionLevel-1;
    byteOffset=bitOffset/8;
    bitOffset-=8*byteOffset;

    nDownNodes=1<<_MESH_DIMENSION_;
    upNode=startNode->upNode;

    for (Level=node.ResolutionLevel;Level>0;Level--) {
      for (nDownNode=0;nDownNode<nDownNodes;nDownNode++) if (upNode->downNode[nDownNode]==startNode) break;
      if (nDownNode==nDownNodes) exit(__LINE__,__FILE__,"Error: cannot find the down node");

      //Determine the coordinates of the down node
      k=nDownNode/4;
      nDownNode-=4*k;

      j=nDownNode/2;
      nDownNode-=2*j;

      i=nDownNode;

      //save the coordinates
      if (k==1) mask=(unsigned char)(mask|(1<<bitOffset));
      if (--bitOffset==-1) {
        node.id[byteOffset]=mask;
        mask=0,bitOffset=7,byteOffset--;
      }

      if (j==1) mask=(unsigned char)(mask|(1<<bitOffset));
      if (--bitOffset==-1) {
        node.id[byteOffset]=mask;
        mask=0,bitOffset=7,byteOffset--;
      }

      if (i==1) mask=(unsigned char)(mask|(1<<bitOffset));
      if (--bitOffset==-1) {
        node.id[byteOffset]=mask;
        mask=0,bitOffset=7,byteOffset--;
      }

      startNode=upNode;
      upNode=startNode->upNode;
    }
  }





  //The MPI routines used in the mesh
  double CalculateTotalParallelLoadMeasure(cTreeNodeAMR<cBlockAMR>* startNode,double *ThreadLoad=NULL) {
    double res=0.0;
    static CMPI_channel pipe(100000);

    if (startNode==rootTree) {
      if (ThisThread==0) {
        pipe.openRecvAll();
        if (ThreadLoad!=NULL) for (int thread=0;thread<nTotalThreads;thread++) ThreadLoad[thread]=0.0;
      }
      else pipe.openSend(0);
    }

    startNode->nodeDescriptor.NodeProcessingFlag=_AMR_FALSE_;

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
      startNode->ParallelLoadMeasure=0.0;
      #endif

      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) res+=CalculateTotalParallelLoadMeasure(startNode->downNode[nDownNode],ThreadLoad);
    }
    else {
      #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
      res=startNode->ParallelLoadMeasure;

      //check whether the node load measure is normalized
      if (isfinite(res)==false) {
        char msg[200];
   
        sprintf(msg,"Error: the parallel load measure is not normalized (res=%e)",res);
        exit(__LINE__,__FILE__,msg);
      }


      if (ThisThread==0) {
        double t;

        if (ThreadLoad!=NULL) ThreadLoad[0]+=res;

        for (int thread=1;thread<nTotalThreads;thread++) {
          pipe.recv(t,thread);
          res+=t;
          if (ThreadLoad!=NULL) ThreadLoad[thread]+=t;
        }

        startNode->ParallelLoadMeasure=res;
      }
      else pipe.send(res);
      #elif _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_OFF_
      res=0.0;
      #else
      exit(__LINE__,__FILE__,"Error: unknown option");
      #endif
    }

    if (startNode==rootTree) {
      if (ThisThread==0) pipe.closeRecvAll();
      else pipe.closeSend();
    }

    return res;
  }

  //normalize the load, calculate the total number of blocks, and the number of blocks on each resolution level
  long int  NormalizeParallelLoadMeasure(double Norm,long int *nResolutionLevelBlocks,cTreeNodeAMR<cBlockAMR>* startNode) {
    long int res=0;

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) res+=NormalizeParallelLoadMeasure(Norm,nResolutionLevelBlocks,startNode->downNode[nDownNode]);
    }
    else {
      nResolutionLevelBlocks[startNode->RefinmentLevel]++;
      res=1;

      #if _AMR_PARALLEL_MODE_ == _AMR_PARALLEL_MODE_ON_
      startNode->ParallelLoadMeasure/=Norm;

      //check whether the node load measure is normalized
      if (isfinite(startNode->ParallelLoadMeasure)==false) {
        char msg[200];

        sprintf(msg,"Error: the parallel load measure is not normalized (startNode->ParallelLoadMeasure=%e)\n",startNode->ParallelLoadMeasure);
        exit(__LINE__,__FILE__,msg);
      }

      //add the load to all parent nodes
      cTreeNodeAMR<cBlockAMR>* upNode=startNode->upNode;

      while (upNode!=NULL) {
        upNode->ParallelLoadMeasure+=startNode->ParallelLoadMeasure;
        upNode=upNode->upNode;
      }

      #endif
    }

    return res;
  }

  //attach all down nodes to the new distribution list
  void AttachDownNodesToDistributionList(cTreeNodeAMR<cBlockAMR> **DistributionList,cTreeNodeAMR<cBlockAMR> *startNode) {
    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) AttachDownNodesToDistributionList(DistributionList,startNode->downNode[nDownNode]);
    }
    else {
      startNode->nextNodeThisThread=(*DistributionList);
      startNode->prevNodeThisThread=NULL;

      if (*DistributionList!=NULL) (*DistributionList)->prevNodeThisThread=startNode;
      *DistributionList=startNode;
    }
  }

  /*
  void ResetAMRnodeProcessingFlag(cTreeNodeAMR<cBlockAMR>* startNode) {
    startNode->nodeDescriptor.NodeProcessingFlag=_AMR_FALSE_;

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) ResetAMRnodeProcessingFlag(startNode->downNode[nDownNode]);
    }
  }
  */
  void ResetAMRnodeProcessingFlag() {
    long int nMemoryBank,nTotalMemoryBanks,nnode;
    cTreeNodeAMR<cBlockAMR>* TreeNode;

    nTotalMemoryBanks=treeNodes.dataBufferListPointer;

    for (nMemoryBank=0;nMemoryBank<nTotalMemoryBanks;nMemoryBank++) {
      TreeNode=treeNodes.dataBufferList[nMemoryBank];

      for (nnode=0;nnode<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;nnode++) {
        (TreeNode+nnode)->nodeDescriptor.NodeProcessingFlag=_AMR_FALSE_;
      }
    }
  }

  void ResetParallelLoadMeasure() {
    long int nMemoryBank,nTotalMemoryBanks,nnode;
    cTreeNodeAMR<cBlockAMR>* TreeNode;

    nTotalMemoryBanks=treeNodes.dataBufferListPointer;

    for (nMemoryBank=0;nMemoryBank<nTotalMemoryBanks;nMemoryBank++) {
      TreeNode=treeNodes.dataBufferList[nMemoryBank];

      for (nnode=0;nnode<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;nnode++) {
        (TreeNode+nnode)->ParallelLoadMeasure=0.0;
      }
    }
  }

  void SetConstantParallelLoadMeasure(double ParallelLoad,cTreeNodeAMR<cBlockAMR>* startNode=NULL) {
    if (startNode==NULL) startNode=rootTree;
    startNode->ParallelLoadMeasure=0.0;

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) SetConstantParallelLoadMeasure(ParallelLoad,startNode->downNode[nDownNode]);
    }
    else startNode->ParallelLoadMeasure=ParallelLoad;
  }

  void SetParallelLoadMeasure(double(*ParallelLoad)(cTreeNodeAMR<cBlockAMR>*),cTreeNodeAMR<cBlockAMR>* startNode=NULL) {
    if (startNode==NULL) startNode=rootTree;
    startNode->ParallelLoadMeasure=0.0;

    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) SetParallelLoadMeasure(ParallelLoad,startNode->downNode[nDownNode]);
    }
    else startNode->ParallelLoadMeasure=ParallelLoad(startNode);
  }


#define _BLOCK_INSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_         0
#define _BLOCK_INTERSECTS_DOMAIN__ALLOCATE_TREE_BLOCKS_     1
#define _BLOCK_OUTINSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_     -1

  void AllocateTreeBlocks() {
    int flag=_BLOCK_INSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_;

    int rank;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);

    if (rank==0) {
      printf("$PREFIX: Allocating tree blocks.....  ");
      fflush(stdout);
    }

#if _AMR__CUT_CELL__MODE_ == _AMR__CUT_CELL__MODE__ON_
    if (rootTree->FirstTriangleCutFace!=NULL) flag=_BLOCK_INTERSECTS_DOMAIN__ALLOCATE_TREE_BLOCKS_;
#endif

    AllocateTreeBlocks(rootTree,flag);

    if (rank==0) {
      printf("done.\n");
      fflush(stdout);
    }
  }


  void AllocateTreeBlocks(cTreeNodeAMR<cBlockAMR>* startNode,int DomainIntersectionFlag) {
    static bool ThisThreadBlockFound;
    static long int nAllocatedBlocks;

    if (startNode==rootTree) {
      ThisThreadBlockFound=false,nAllocatedBlocks=0;

      //init the ray tracing module if needed
      CutCell::InitRayTracingModule();
    }


    if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      startNode->Thread=-1;
      if ((startNode->block!=NULL)&&(DeallocateUnusedBlocks==true)) DeallocateBlock(startNode);

      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) {
        bool downNodeDomainIntersectionFlag=DomainIntersectionFlag;

#if _AMR__CUT_CELL__MODE_ == _AMR__CUT_CELL__MODE__ON_
        switch (DomainIntersectionFlag) {
        case _BLOCK_OUTINSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_ :
        case _BLOCK_INSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_ :
          downNodeDomainIntersectionFlag=DomainIntersectionFlag;

          break;
        case _BLOCK_INTERSECTS_DOMAIN__ALLOCATE_TREE_BLOCKS_:
          if (startNode->FirstTriangleCutFace!=NULL) {
            //getermine if the block inside, cross, or outside of the domain
            //determine if the block is intersected by a cut-face
            bool IntersectionFound=false;



/*            for (int nCutFace=0;nCutFace<CutCell::nBoundaryTriangleFaces;nCutFace++) {
              if (CutCell::BoundaryTriangleFaces[nCutFace].BlockIntersection(startNode->downNode[nDownNode]->xmin,startNode->downNode[nDownNode]->xmax,EPS)==true) {
                IntersectionFound=true;
                downNodeDomainIntersectionFlag=_BLOCK_INTERSECTS_DOMAIN__ALLOCATE_TREE_BLOCKS_;
                break;
              }
            }*/

            for (CutCell::cTriangleFaceDescriptor* Descriptor=startNode->FirstTriangleCutFace;Descriptor!=NULL;Descriptor=Descriptor->next) {
              if (Descriptor->TriangleFace->BlockIntersection(startNode->downNode[nDownNode]->xmin,startNode->downNode[nDownNode]->xmax,EPS)==true) {
                IntersectionFound=true;
                downNodeDomainIntersectionFlag=_BLOCK_INTERSECTS_DOMAIN__ALLOCATE_TREE_BLOCKS_;
                break;
              }
            }


            //if intersection with the cut-face is not found, check the middle point of the block is within the domain
            if (IntersectionFound==false) {
              int idim;
              double xMiddle[3];

              for (idim=0;idim<DIM;idim++) xMiddle[idim]=0.5*(startNode->downNode[nDownNode]->xmin[idim]+startNode->downNode[nDownNode]->xmax[idim]);

              if (CutCell::CheckPointInsideDomain(xMiddle,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,EPS)==true) {
                downNodeDomainIntersectionFlag=_BLOCK_INSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_;
              }
              else downNodeDomainIntersectionFlag=_BLOCK_OUTINSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_;
            }
          }

          break;
        default:
          exit(__LINE__,__FILE__,"error: unknown option");
        }
#endif

        if ((downNodeDomainIntersectionFlag==_BLOCK_INSIDE_DOMAIN__ALLOCATE_TREE_BLOCKS_)||(downNodeDomainIntersectionFlag==_BLOCK_INTERSECTS_DOMAIN__ALLOCATE_TREE_BLOCKS_)) {
          AllocateTreeBlocks(startNode->downNode[nDownNode],downNodeDomainIntersectionFlag);
        }
      }
    }
    else {
      if (startNode->Thread==ThisThread) {
        ThisThreadBlockFound=true,nAllocatedBlocks++;


//====================== DEBUG =======================

static long int TmpAllocationCounter=0;

TmpAllocationCounter++;

/*
if (TmpAllocationCounter==2437) {
  *DiagnospticMessageStream << __LINE__ << std::endl;
}
*/

//checkMeshConsistency(rootTree);
//====================== END DEBUG =============



        if (startNode->block==NULL) AllocateBlock(startNode);


//====================== DEBUG =======================

//checkMeshConsistency(rootTree);
//====================== END DEBUG =============



      }
      else if ((startNode->Thread<0)||(startNode->Thread>=nTotalThreads)) exit(__LINE__,__FILE__,"Error: the thread number is out of range");
      else if (startNode->block!=NULL) exit(__LINE__,__FILE__,"Error: a block is allocated for another thread");
    }

    if (startNode==rootTree) {
      if (_AMR__NO_BLOCKS_FOUND__EXIT_MODE_ == _ON_AMR_MESH_) {
        if (ThisThreadBlockFound==false) exit(__LINE__,__FILE__,"Error: no blocks found that belong to thins thread");
      }

      int mpiInitFlag;
      MPI_Initialized(&mpiInitFlag);

      if (mpiInitFlag==true) {
        long int thread=0,*buffer=NULL;

        if (ThisThread==0) {
          buffer=new long int [nTotalThreads];
          buffer[0]=nAllocatedBlocks;

          long int bufferRecv[nTotalThreads];
          MPI_Gather(buffer,1,MPI_LONG,bufferRecv,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
          memcpy(buffer,bufferRecv,nTotalThreads*sizeof(long int));

          fprintf(DiagnospticMessageStream,"$PREFIX:Blocks Allocation Report:\n$PREFIX: Thread\tAllocatedBlocks\n");
          for (thread=0;thread<nTotalThreads;thread++) fprintf(DiagnospticMessageStream,"$PREFIX:%ld\t%ld\n",thread,buffer[thread]);

          delete [] buffer;
        }
        else MPI_Gather(&nAllocatedBlocks,1,MPI_LONG,buffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);

        MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

        //allocate the domain's boundary blocks
#if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
        cTreeNodeAMR<cBlockAMR> *node;

        nAllocatedBlocks=0;

        for (thread=0;thread<nTotalThreads;thread++) {
          node=DomainBoundaryLayerNodesList[thread];

          while (node!=NULL) {
            nAllocatedBlocks++;
            if (node->block==NULL) AllocateBlock(node);

            node=node->nextNodeThisThread;
          }
        }

        if (ThisThread==0) {
          buffer=new long int [nTotalThreads];
          buffer[0]=nAllocatedBlocks;

          long int bufferRecv[nTotalThreads];
          MPI_Gather(buffer,1,MPI_LONG,bufferRecv,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
          memcpy(buffer,bufferRecv,nTotalThreads*sizeof(long int));

          fprintf(DiagnospticMessageStream,"$PREFIX:Blocks Allocation Report:\n$PREFIX: Thread\tAllocated Domain's Boundary Blocks\n");
          for (thread=0;thread<nTotalThreads;thread++) fprintf(DiagnospticMessageStream,"$PREFIX:%ld\t%ld\n",thread,buffer[thread]);

          delete [] buffer;
        }
        else MPI_Gather(&nAllocatedBlocks,1,MPI_LONG,buffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);

        MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
#endif

      }
    }

  }


  void CreateMortonSpaceFillingCurve(cTreeNodeAMR<cBlockAMR> **startNodeFillingCurve,cTreeNodeAMR<cBlockAMR>* upLevelNode) {
    int nd;

    if (upLevelNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      //there are levels down the tree - they should be inserted into the curve

      //connect downNode's
      for (nd=0;nd<(1<<_MESH_DIMENSION_)-1;nd++) upLevelNode->downNode[nd]->FillingCurveNextNode=upLevelNode->downNode[nd+1];
      for (nd=1;nd<(1<<_MESH_DIMENSION_);nd++) upLevelNode->downNode[nd]->FillingCurvePrevNode=upLevelNode->downNode[nd-1];

      //reconnect the nodes with the elements of the curve
      cTreeNodeAMR<cBlockAMR>* t;

      upLevelNode->downNode[0]->FillingCurvePrevNode=(t=upLevelNode->FillingCurvePrevNode);
      if (t!=NULL) t->FillingCurveNextNode=upLevelNode->downNode[0];

      upLevelNode->downNode[(1<<_MESH_DIMENSION_)-1]->FillingCurveNextNode=(t=upLevelNode->FillingCurveNextNode);
      if (t!=NULL) t->FillingCurvePrevNode=upLevelNode->downNode[(1<<_MESH_DIMENSION_)-1];

      upLevelNode->FillingCurveNextNode=NULL;
      upLevelNode->FillingCurvePrevNode=NULL;

      //update the start point of the curve
      if (startNodeFillingCurve!=NULL) *startNodeFillingCurve=upLevelNode->downNode[0];

      //insert into the curve nodes of higher levels
      for (nd=0;nd<(1<<_MESH_DIMENSION_);nd++) {
        CreateMortonSpaceFillingCurve(startNodeFillingCurve,upLevelNode->downNode[nd]);

        //startNodeFillingCurve=NULL for other 'downNode's
        startNodeFillingCurve=NULL;
      }
    }
  }


/*
  void CreateSpaceFillingCurve(cTreeNodeAMR<cBlockAMR> **startNodeFillingCurve,int UpperResolutionLevel) {
    int idim,i,j,k,iMax,jMax,kMax;
    double xProbe[3]={0.0,0.0,0.0},dxCurve[3]={0.0,0.0,0.0};
    cTreeNodeAMR<cBlockAMR> *TreeNode=NULL;
    int iStart=0,iIncrement=-1,jStart=0,jIncrement=-1;

    ResetAMRnodeProcessingFlag();
    *startNodeFillingCurve=NULL;

#if _MESH_DIMENSION_ == 1
    iMax=1<<UpperResolutionLevel,jMax=1,kMax=1;
#elif _MESH_DIMENSION_ == 2
    iMax=1<<UpperResolutionLevel,jMax=1<<UpperResolutionLevel,kMax=1;
#elif _MESH_DIMENSION_ == 3
    iMax=1<<UpperResolutionLevel,jMax=1<<UpperResolutionLevel,kMax=1<<UpperResolutionLevel;
#else
    exit(__LINE__,__FILE__,"Error: unknown option");
#endif


    for (idim=0;idim<_MESH_DIMENSION_;idim++) dxCurve[idim]=(xGlobalMax[idim]-xGlobalMin[idim])/(1<<UpperResolutionLevel);


    for (k=0;k<kMax;k++) {
      if (_MESH_DIMENSION_==3) xProbe[2]=xGlobalMin[2]+(k+0.5)*dxCurve[2];

      if (jIncrement==1) jIncrement=-1,jStart=jMax-1;
      else jIncrement=1,jStart=0;

      for (j=jStart;(j<jMax)&&(j>=0);j+=jIncrement) {
        if (_MESH_DIMENSION_>=2) xProbe[1]=xGlobalMin[1]+(j+0.5)*dxCurve[1];

        if (iIncrement==1) iIncrement=-1,iStart=iMax-1;
        else iIncrement=1,iStart=0;

        for (i=iStart;(i<iMax)&&(i>=0);i+=iIncrement) {
          xProbe[0]=xGlobalMin[0]+(i+0.5)*dxCurve[0];
          TreeNode=findTreeNodeLimitedResolutionLevel(xProbe,UpperResolutionLevel,TreeNode);
          if (TreeNode==NULL) exit(__LINE__,__FILE__,"Error: cannot find the location");

          if (TreeNode->nodeDescriptor.NodeProcessingFlag==_AMR_FALSE_) {
            TreeNode->nodeDescriptor.NodeProcessingFlag=_AMR_TRUE_;
            TreeNode->FillingCurveNextNode=(*startNodeFillingCurve);
            *startNodeFillingCurve=TreeNode;
          }
        }
      }
    }
  }
*/


  long int nParallelListRedistributions;

  void CreateNewParallelDistributionLists(int userDefinedCodeForSendingBlockData=-1) {
    double LoadMeasureNormal;
    long int nTotalBlocks,nResolutionLevelBlocks[_MAX_REFINMENT_LEVEL_+1];
    int i,nLevel;

    static cTreeNodeAMR<cBlockAMR> *startNodeFillingCurve=NULL;

    //increment the counter of the load redistributions
    nParallelListRedistributions++;
    nMeshModificationCounter++;

    //calculate and normalized the load measure
    for (nLevel=0;nLevel<=_MAX_REFINMENT_LEVEL_;nLevel++) nResolutionLevelBlocks[nLevel]=0;

    double InitialProcessorLoad[nTotalThreads];
    LoadMeasureNormal=CalculateTotalParallelLoadMeasure(rootTree,InitialProcessorLoad)/nTotalThreads;

    MPI_Bcast(&LoadMeasureNormal,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

    if (LoadMeasureNormal<=0.0) {
      SetConstantParallelLoadMeasure(1.0,rootTree);
      LoadMeasureNormal=CalculateTotalParallelLoadMeasure(rootTree,InitialProcessorLoad)/nTotalThreads;
    }

    if (ThisThread==0) {
      cTreeNodeAMR<cBlockAMR> *ptr;
      long int nblocks;

      fprintf(DiagnospticMessageStream,"$PREFIX:Initial Cumulative Parallel Load Distribution\n$PREFIX:Thread\tLoad\tNormalized Load\tBlock's Number\n");

      for (int t=0;t<nTotalThreads;t++) {
        for (nblocks=0,ptr=ParallelNodesDistributionList[t];ptr!=NULL;ptr=ptr->nextNodeThisThread) nblocks++;

        fprintf(DiagnospticMessageStream,"$PREFIX:%i\t%8.2e\t%8.2e\t%ld\n",t,InitialProcessorLoad[t],InitialProcessorLoad[t]/LoadMeasureNormal,nblocks);
      }
    }




    //normalize the load
    nTotalBlocks=NormalizeParallelLoadMeasure(LoadMeasureNormal,nResolutionLevelBlocks,rootTree);

    //Recalculate the space filling curve if its needed
    if (meshModifiedFlag_CreateNewSpaceFillingCurve==true) {
//      long int nBalancingLevel,nBalancigElements,nBalancigElementsUpLevels;

      meshModifiedFlag_CreateNewSpaceFillingCurve=false;

      rootTree->FillingCurveNextNode=NULL;
      rootTree->FillingCurvePrevNode=NULL;

      startNodeFillingCurve=rootTree;

      CreateMortonSpaceFillingCurve(&startNodeFillingCurve,rootTree);

      /*
      for (nBalancingLevel=0;nBalancingLevel<=_MAX_REFINMENT_LEVEL_;nBalancingLevel++) {
        for (nBalancigElements=0,i=0;i<=nBalancingLevel;i++) nBalancigElements+=nResolutionLevelBlocks[i];

        for (nBalancigElementsUpLevels=nBalancingLevel+1;i<=_MAX_REFINMENT_LEVEL_;i++) {
          nBalancigElementsUpLevels+=nResolutionLevelBlocks[i];
 //         nBalancigElements+=(long int)(double(nResolutionLevelBlocks[i])/(1<<(_MESH_DIMENSION_*(i-nBalancingLevel))));
        }

        if ((nBalancigElementsUpLevels==0)||((nBalancigElementsUpLevels+nBalancigElements)*0.1>nBalancigElementsUpLevels)) break;
      }

      if (nBalancingLevel>_MAX_REFINMENT_LEVEL_) nBalancingLevel=_MAX_REFINMENT_LEVEL_;

      CreateSpaceFillingCurve(&startNodeFillingCurve,nBalancingLevel);
      */
    }


    //redistribute the parallel load
    double CumulativeProcessorLoad=0.0;
    cTreeNodeAMR<cBlockAMR>* CurveNode=startNodeFillingCurve;
    int nCurrentProcessorBalancing=0;
    CMPI_channel pipe(100000);
    cTreeNodeAMR<cBlockAMR> *node;

    for (i=0;i<nTotalThreads;i++) {
      ParallelNodesDistributionList[i]=NULL;

      //deallocate the boundary layer
#if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
      if (DomainBoundaryLayerNodesList[i]!=NULL) { //exit(__LINE__,__FILE__,"not implemented");
        node=DomainBoundaryLayerNodesList[i];

        while (node!=NULL) {
          DeallocateBlock(node);
          node=node->nextNodeThisThread;
        }

        DomainBoundaryLayerNodesList[i]=NULL;
      }
#endif
    }

    //redistribute the load on the root processor
    double *CumulativeThreadLoad=new double [nTotalThreads];
    cTreeNodeAMR<cBlockAMR>** ThreadStartNode=new cTreeNodeAMR<cBlockAMR>* [nTotalThreads];
    int thread;

    for (thread=0;thread<nTotalThreads;thread++) CumulativeThreadLoad[thread]=0.0,ThreadStartNode[thread]=NULL;

    if (ThisThread==0) {
      ThreadStartNode[0]=CurveNode;
      nCurrentProcessorBalancing=0;

      while (CurveNode!=NULL) {
        CumulativeProcessorLoad+=CurveNode->ParallelLoadMeasure;
        CumulativeThreadLoad[nCurrentProcessorBalancing]+=CurveNode->ParallelLoadMeasure;

        if ((isfinite(CumulativeProcessorLoad)==false)||(isfinite(CumulativeThreadLoad[nCurrentProcessorBalancing])==false)) {
          char msg[500];

          sprintf(msg,"Error: the parallel load measure is not finite (CumulativeProcessorLoad=%e, CumulativeThreadLoad[nCurrentProcessorBalancing]=%e, nCurrentProcessorBalancing=%i)",CumulativeProcessorLoad,CumulativeThreadLoad[nCurrentProcessorBalancing],nCurrentProcessorBalancing); 
          exit(__LINE__,__FILE__,msg);
        }  

        //increment the processor number if needed
        if ((CumulativeProcessorLoad>1.0+nCurrentProcessorBalancing)&&(nCurrentProcessorBalancing!=nTotalThreads-1)) {
          nCurrentProcessorBalancing++;
          ThreadStartNode[nCurrentProcessorBalancing]=CurveNode->FillingCurveNextNode;
        }

        CurveNode=CurveNode->FillingCurveNextNode;
      }

      //fine tuning of the processor's load
      for (thread=0;thread<nTotalThreads;thread++) {
         const double UpperLoadLimit=1.05;

         if (CumulativeThreadLoad[thread]>UpperLoadLimit) for (int nLoadOptimizationPass=0;nLoadOptimizationPass<2;nLoadOptimizationPass++) {
           //check if the load can be exchanged with one of the neighboring
           int smallLoadNeighbour=-2;
           double neibLoad=CumulativeThreadLoad[thread];

           if (thread!=0) if (CumulativeThreadLoad[thread-1]<neibLoad) neibLoad=CumulativeThreadLoad[thread-1],smallLoadNeighbour=thread-1;
           if (thread!=nTotalThreads-1) if (CumulativeThreadLoad[thread+1]<neibLoad) neibLoad=CumulativeThreadLoad[thread+1],smallLoadNeighbour=thread+1;

           if (smallLoadNeighbour!=-1) {
             //try to transfer some load to the neighboring processor
             if ((smallLoadNeighbour==thread-1)&&(ThreadStartNode[thread-1]!=NULL)) {
               while (CumulativeThreadLoad[thread-1]+ThreadStartNode[thread]->ParallelLoadMeasure<max(1.0,CumulativeThreadLoad[thread]-ThreadStartNode[thread]->ParallelLoadMeasure)) {

                 //ckeck whether there are any blocks to exchange
                 if (thread<nTotalThreads-1) {
                   if (ThreadStartNode[thread]->FillingCurveNextNode==ThreadStartNode[thread+1]) {
                     //'thread' has only one block
                     break;
                   }
                 } 
                 else {
                   if (ThreadStartNode[thread]->FillingCurveNextNode==NULL) {
                     break;
                   }
                 }

                 CumulativeThreadLoad[thread-1]+=ThreadStartNode[thread]->ParallelLoadMeasure;
                 CumulativeThreadLoad[thread]-=ThreadStartNode[thread]->ParallelLoadMeasure;

                 ThreadStartNode[thread]=ThreadStartNode[thread]->FillingCurveNextNode;
               }
             }
             else if ((smallLoadNeighbour==thread+1)&&(ThreadStartNode[thread+1]!=NULL)) {
               while (CumulativeThreadLoad[thread+1]+ThreadStartNode[thread+1]->FillingCurvePrevNode->ParallelLoadMeasure<max(1.0,CumulativeThreadLoad[thread]-ThreadStartNode[thread+1]->FillingCurvePrevNode->ParallelLoadMeasure)) {

                 //ckeck whether there are any blocks to exchange
                 if (thread<nTotalThreads-1) {
                   if (ThreadStartNode[thread]->FillingCurveNextNode==ThreadStartNode[thread+1]) {
                     //'thread' has only one block
                     break;
                   }
                 }
                 else {
                   if (ThreadStartNode[thread]->FillingCurveNextNode==NULL) {
                     break;
                   }
                 }

                 CumulativeThreadLoad[thread+1]+=ThreadStartNode[thread+1]->FillingCurvePrevNode->ParallelLoadMeasure;
                 CumulativeThreadLoad[thread]-=ThreadStartNode[thread+1]->FillingCurvePrevNode->ParallelLoadMeasure;

                 ThreadStartNode[thread+1]=ThreadStartNode[thread+1]->FillingCurvePrevNode;
               }
             }
           }

         }
      }
    }

    //Exchange the load distribution
    cTreeNodeAMR<cBlockAMR>* nextThreadStartNode=startNodeFillingCurve;

    CurveNode=startNodeFillingCurve;
    thread=-1;

    pipe.openBcast(0);

    while (CurveNode!=NULL) {
      if (CurveNode==nextThreadStartNode) {
        thread+=1;
        nextThreadStartNode=(thread!=nTotalThreads-1) ? ThreadStartNode[thread+1] : NULL;
      }

      if (ThisThread==0) pipe.send(thread);
      else pipe.recv(thread,0);

      //add the node to the processor list
      CurveNode->nextNodeThisThread=ParallelNodesDistributionList[thread];
      CurveNode->prevNodeThisThread=NULL;

      if (ParallelNodesDistributionList[thread]!=NULL) ParallelNodesDistributionList[thread]->prevNodeThisThread=CurveNode;
      ParallelNodesDistributionList[thread]=CurveNode;

      CurveNode=CurveNode->FillingCurveNextNode;
    }

    pipe.closeBcast();

    //check if the nodes' distribution is the same on all processors


    if (ThisThread==0) {
      int t;
      cTreeNodeAMR<cBlockAMR> *ptr;
      cAMRnodeID id,idtemp;

      //cumulutive load distribution
      double newCumulativeParallelLoadMeasure[nTotalThreads];
      double TotalParallelLoadMeasure=0.0;

      for (t=0;t<nTotalThreads;t++) newCumulativeParallelLoadMeasure[t]=0.0;

      pipe.openRecvAll();

      for (t=0;t<nTotalThreads;t++) for (ptr=ParallelNodesDistributionList[t];ptr!=NULL;ptr=ptr->nextNodeThisThread) {
        GetAMRnodeID(id,ptr);

        TotalParallelLoadMeasure+=ptr->ParallelLoadMeasure;
        newCumulativeParallelLoadMeasure[t]+=ptr->ParallelLoadMeasure;

        for (thread=1;thread<nTotalThreads;thread++) {
          pipe.recv(idtemp,thread);

          if (id!=idtemp) exit(__LINE__,__FILE__,"Error: the processor's distribution of the tree is not consistent");
        }
      }

      pipe.closeRecvAll();

      fprintf(DiagnospticMessageStream,"$PREFIX:Cumulative Parallel Load Distribution\n$PREFIX:Thread\tLoad\tNormalized Load\tBlock's Number\n");

      int minThreadBlockNumber=-1,maxThreadBlockNumber=-1;
      double minLoadMeasure=-1.0,maxLoadMeasure=-1.0;
      int minLoadThread=-1,maxLoadThread=-1;
      double minBlockLoad=-1.0,maxBlockLoad=-1.0;
      cTreeNodeAMR<cBlockAMR> *maxLoadBlock=NULL;

      for (t=0;t<nTotalThreads;t++) {
        long int nblocks;

        for (nblocks=0,ptr=ParallelNodesDistributionList[t];ptr!=NULL;ptr=ptr->nextNodeThisThread) {
          nblocks++;

          if ((minBlockLoad<0.0)||(minBlockLoad>ptr->ParallelLoadMeasure)) minBlockLoad=ptr->ParallelLoadMeasure;
          if ((maxBlockLoad<0.0)||(maxBlockLoad<ptr->ParallelLoadMeasure)) maxBlockLoad=ptr->ParallelLoadMeasure,maxLoadBlock=ptr;
        }

        fprintf(DiagnospticMessageStream,"$PREFIX:%i\t%8.2e\t%8.2e\t%ld\n",t,newCumulativeParallelLoadMeasure[t],nTotalThreads*newCumulativeParallelLoadMeasure[t]/TotalParallelLoadMeasure,nblocks);

        if ((minThreadBlockNumber==-1)||(minThreadBlockNumber>nblocks)) minThreadBlockNumber=nblocks;
        if ((maxThreadBlockNumber==-1)||(maxThreadBlockNumber<nblocks)) maxThreadBlockNumber=nblocks;

        if ((minLoadMeasure<0.0)||(minLoadMeasure>nTotalThreads*newCumulativeParallelLoadMeasure[t]/TotalParallelLoadMeasure)) {
          minLoadMeasure=nTotalThreads*newCumulativeParallelLoadMeasure[t]/TotalParallelLoadMeasure;
          minLoadThread=t;
        }

        if ((maxLoadMeasure<0.0)||(maxLoadMeasure<nTotalThreads*newCumulativeParallelLoadMeasure[t]/TotalParallelLoadMeasure)) {
          maxLoadMeasure=nTotalThreads*newCumulativeParallelLoadMeasure[t]/TotalParallelLoadMeasure;
          maxLoadThread=t;
        }
      }

      fprintf(DiagnospticMessageStream,"$PREFIX:Min number of blocks per processor: %i\n$PREFIX:Max number of blocks per processor: %i\n",minThreadBlockNumber,maxThreadBlockNumber);
      fprintf(DiagnospticMessageStream,"$PREFIX:Min processor load: %e, thread=%i\n",minLoadMeasure,minLoadThread);
      fprintf(DiagnospticMessageStream,"$PREFIX:Max processor load: %e, thread=%i\n",maxLoadMeasure,maxLoadThread);
      fprintf(DiagnospticMessageStream,"$PREFIX:Individual block load: min=%e, max=%e\n",minBlockLoad,maxBlockLoad);

      fprintf(DiagnospticMessageStream,"$PREFIX:Parameters of the block with the maximum load\n");
      fprintf(DiagnospticMessageStream,"$PREFIX:Block->Temp_ID=%ld\n",maxLoadBlock->Temp_ID);

      //the position of the block's nodes
      double middleX[3]={0.0,0.0,0.0},xnode[3];
      int j,k,idim;

      for (i=0;i<2;i++) for (j=0;j<((_MESH_DIMENSION_>1) ? 2 : 1);j++) for (k=0;k<((_MESH_DIMENSION_>2) ? 2 : 1);k++) {
        maxLoadBlock->GetCornerNodePosition(xnode,i*_BLOCK_CELLS_X_,j*_BLOCK_CELLS_Y_,k*_BLOCK_CELLS_Z_);

        fprintf(DiagnospticMessageStream,"$PREFIX:(i,j,k)=(%i,%i,%i): x=",i,j,k);
        for (idim=0;idim<_MESH_DIMENSION_;idim++) {
          fprintf(DiagnospticMessageStream,"%e ",xnode[idim]);
          middleX[idim]+=xnode[idim];
        }

        fprintf(DiagnospticMessageStream,"\n");
      }


      fprintf(DiagnospticMessageStream,"$PREFIX:Middle block's coordinates=");
      for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(DiagnospticMessageStream,"%e ",middleX[idim]/(1<<_MESH_DIMENSION_));
      fprintf(DiagnospticMessageStream,"\n\n");

    }
    else {
      int t;
      cTreeNodeAMR<cBlockAMR> *ptr;
      cAMRnodeID id;

      pipe.openSend(0);

      for (t=0;t<nTotalThreads;t++) for (ptr=ParallelNodesDistributionList[t];ptr!=NULL;ptr=ptr->nextNodeThisThread) {
        GetAMRnodeID(id,ptr);
        pipe.send(id);
      }

      pipe.closeSend();
    }


    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    //de-allocate the temporary data buffers
    delete [] CumulativeThreadLoad;
    delete [] ThreadStartNode;

    CumulativeThreadLoad=NULL,ThreadStartNode=NULL;

    //check that all blocks are presented in the new blocks' distribution lists
    long int nDistributedNodes=0;
    //cTreeNodeAMR<cBlockAMR> *node;

    ResetAMRnodeProcessingFlag();

    for (/*int*/ thread=0;thread<nTotalThreads;thread++) {
      node=ParallelNodesDistributionList[thread];

      while (node!=NULL) {
        nDistributedNodes++;

        if (node->nodeDescriptor.NodeProcessingFlag==_AMR_TRUE_) exit(__LINE__,__FILE__,"Error: the same node in present in the nodes' distribution list twice");
        node->nodeDescriptor.NodeProcessingFlag=_AMR_TRUE_;

        if ((node->nextNodeThisThread!=NULL)&&(node->nextNodeThisThread==node->prevNodeThisThread)) exit(__LINE__,__FILE__,"Error: prev==next");

        node=node->nextNodeThisThread;
      }
    }

    if (nTotalBlocks!=nDistributedNodes) exit(__LINE__,__FILE__,"Error: some nodes were lost during the load rebalancing");

//    delete [] nResolutionLevelBlocks;

    //clear the node's load sample
    SetConstantParallelLoadMeasure(0.0,rootTree);

    /*
    //deallocate blocks from the subdomain boundary layer
    node=DomainBoundaryLayerNodesList[ThisThread];

    while (node!=NULL) {
      DeallocateBlock(node);
      node=node->nextNodeThisThread;
    }
    */

    //move the cells between processes according the new distribution
    pipe.openSend(0);
    pipe.openRecv(0);
    int pipeLastRecvThread=0;

    for (thread=0;thread<nTotalThreads;thread++) {
      node=ParallelNodesDistributionList[thread];

      while (node!=NULL) {
        if ((userDefinedCodeForSendingBlockData!=-1)&&(node->Thread!=thread)&&((ThisThread==node->Thread)||(ThisThread==thread))) { //the block is moved from processor=node->Thread to processor=thread
//          exit(__LINE__,__FILE__,"blocks' exchange is not implemented"); //userDefinedCodeForSendingBlockData

          if (ThisThread==thread) { //the block is moved to processor 'thread'
            if (pipeLastRecvThread!=node->Thread) {
              pipe.RedirectRecvBuffer(node->Thread);
              pipeLastRecvThread=node->Thread;
            }

            AllocateBlock(node);
//            InitCellMeasure(node);
            node->block->recvMoveBlockAnotherProcessor(&pipe,node->Thread);
          }
          else { //the block is moved out from the processor
            if (pipe.sendThread!=thread) pipe.RedirectSendBuffer(thread);

            node->block->sendMoveBlockAnotherProcessor(&pipe);
            DeallocateBlock(node);
          }

        }

        //update the node's thread
        node->Thread=thread;
        node=node->nextNodeThisThread;
      }

      pipe.flush();

#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
//      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
#endif

    }

    pipe.closeSend();
    pipe.closeRecv(pipeLastRecvThread);

    //create the Send/Recv flags vectors
    for (thread=0;thread<nTotalThreads;thread++) for (i=0;i<nTotalThreads;i++) ParallelSendRecvMap[thread][i]=false;




    for (thread=0;thread<nTotalThreads;thread++) {
      node=ParallelNodesDistributionList[thread];

      while (node!=NULL) {
        //search through the neighbors of the blocks
#if _MESH_DIMENSION_ == 1
//for (i=0;i<2;i++) neibNodeFace[i]=NULL;
        for (i=0;i<2;i++) if (node->neibNodeFace[i]!=NULL) if (node->neibNodeFace[i]->Thread!=node->Thread) {
          ParallelSendRecvMap[node->neibNodeFace[i]->Thread][node->Thread]=true;
          ParallelSendRecvMap[node->Thread][node->neibNodeFace[i]->Thread]=true;
        }

#elif _MESH_DIMENSION_ == 2
       for (i=0;i<4*2;i++) if (node->neibNodeFace[i]!=NULL) if (node->neibNodeFace[i]->Thread!=node->Thread) {
         ParallelSendRecvMap[node->neibNodeFace[i]->Thread][node->Thread]=true;
         ParallelSendRecvMap[node->Thread][node->neibNodeFace[i]->Thread]=true;
       }

       for (i=0;i<4;i++) if (node->neibNodeCorner[i]!=NULL) if (node->neibNodeCorner[i]->Thread!=node->Thread) {
         ParallelSendRecvMap[node->neibNodeCorner[i]->Thread][node->Thread]=true;
         ParallelSendRecvMap[node->Thread][node->neibNodeCorner[i]->Thread]=true;
       }
#elif _MESH_DIMENSION_ == 3

//for (i=0;i<6*4;i++) neibNodeFace[i]=NULL;
//for (i=0;i<8;i++) neibNodeCorner[i]=NULL;
//for (i=0;i<12*2;i++) neibNodeEdge[i]=NULL;

       for (i=0;i<6*4;i++) if (node->neibNodeFace[i]!=NULL) if (node->neibNodeFace[i]->Thread!=node->Thread) {
         ParallelSendRecvMap[node->neibNodeFace[i]->Thread][node->Thread]=true;
         ParallelSendRecvMap[node->Thread][node->neibNodeFace[i]->Thread]=true;
       }

       for (i=0;i<8;i++) if (node->neibNodeCorner[i]!=NULL) if (node->neibNodeCorner[i]->Thread!=node->Thread) {
         ParallelSendRecvMap[node->neibNodeCorner[i]->Thread][node->Thread]=true;
         ParallelSendRecvMap[node->Thread][node->neibNodeCorner[i]->Thread]=true;
       }

       for (i=0;i<12*2;i++) if (node->neibNodeEdge[i]!=NULL) if (node->neibNodeEdge[i]->Thread!=node->Thread) {
         ParallelSendRecvMap[node->neibNodeEdge[i]->Thread][node->Thread]=true;
         ParallelSendRecvMap[node->Thread][node->neibNodeEdge[i]->Thread]=true;
       }

#else
        exit(__LINE__,__FILE__,"wrong option");
#endif

        node=node->nextNodeThisThread;
      }
    }

    //create the list of boundary layer, allocated and init the corresponding blocks
#if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_
    for (i=0;i<nTotalThreads;i++) DomainBoundaryLayerNodesList[i]=NULL;

    InitDomainBoundaryLayer(rootTree);

    //allocate blocks from the subdomain boundary layer
    for (i=0;i<nTotalThreads;i++) if (i!=ThisThread) {
      node=DomainBoundaryLayerNodesList[i];

      while (node!=NULL) {
        AllocateBlock(node);
        node=node->nextNodeThisThread;
      }
    }

    //exchenge data from the boundary layer blocks
    if (blocks.usedElements()!=0) ParallelBlockDataExchange();

 //   if (blocks.usedElements()!=0)  exit(__LINE__,__FILE__,"not implemented"); //the blocks are not exists -> simply create the boundary layer list

#endif

#if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
#endif

    {//broadcast to everyone if changes have been made
      int send_flag, recv_flag;
      send_flag = meshModifiedFlag_CountMeshElements;
      MPI_Allreduce(&send_flag,&recv_flag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      if(recv_flag > 0) meshModifiedFlag_CountMeshElements=true;
    }
  }

  void ParallelBlockDataExchange() {
    int From,To,i,pipeLastRecvThread;
    CMPI_channel pipe(100000);

    cTreeNodeAMR<cBlockAMR> *sendNode,*recvNode;
    cAMRnodeID nodeid;

    pipe.openSend(0);
    pipe.openRecv(0);
    pipeLastRecvThread=0;

    //communication signals
    const int _Next_Node_SIGNAL_=0;
    const int _End_Communication_SIGNAL_=1;
    int Signal;

    /*
    //the list of nodes to be send
    const int SendNodelListCounterMax=60;
    int SendNodeListCounter;
    cTreeNodeAMR<cBlockAMR>* SendNodeList[SendNodelListCounterMax];
    bool found;
*/

    //the data exchange loop
    for (From=0;From<nTotalThreads;From++) for (To=0;To<nTotalThreads;To++) if ((From!=To)&&(ParallelSendRecvMap[From][To]==true)) {
//      SendNodeListCounter=0;


      //the part of the sender
      if (ThisThread==From) {
        //redirect the send pipe buffers
        pipe.RedirectSendBuffer(To);

        //reset the proceesed flaf for the blocks to be send
        //send the nodes' data
        for (recvNode=DomainBoundaryLayerNodesList[To];recvNode!=NULL;recvNode=recvNode->nextNodeThisThread) {



//================  DEBUG ========================

          /*
if (ThisThread==1) if ((pow(recvNode->xmin[0]+500.0,2)+pow(recvNode->xmin[1]+1000.0,2)+pow(recvNode->xmin[2]+500.0,2)<0.00001)||(recvNode->Temp_ID==14)) {
   *DiagnospticMessageStream << __LINE__ << std::endl;
}
*/
//================ END DEBUG =====================


#if _MESH_DIMENSION_ == 2
#define _AMR_ParallelBlockDataExchange_SEND_CORNER_NODES_
#elif _MESH_DIMENSION_ == 3
#define _AMR_ParallelBlockDataExchange_SEND_CORNER_NODES_
#endif

#ifdef _AMR_ParallelBlockDataExchange_SEND_CORNER_NODES_
          for (i=0;i<(1<<_MESH_DIMENSION_);i++) if ((sendNode=recvNode->neibNodeCorner[i])!=NULL) if (sendNode->Thread==From) sendNode->nodeDescriptor.NodeProcessingFlag=_AMR_FALSE_;
#endif



          //add face-neighbors into the send list
//#if _MESH_DIMENSION_ == 2
//#define _AMR_ParallelBlockDataExchange_SEND_FACES_
//#elif _MESH_DIMENSION_ == 3
//#define _AMR_ParallelBlockDataExchange_SEND_FACES_
//#endif

//#ifdef _AMR_ParallelBlockDataExchange_SEND_FACES_
         for (i=0;i<_MESH_DIMENSION_*(1<<_MESH_DIMENSION_);i++) if ((sendNode=recvNode->neibNodeFace[i])!=NULL) if (sendNode->Thread==From) sendNode->nodeDescriptor.NodeProcessingFlag=_AMR_FALSE_;
         /*
           //Search the send-list if the node is already there
           for (found=false,k=0;k<SendNodeListCounter;k++) {
             if (sendNode==SendNodeList[k]) {
               found=true;
               break;
             }
           }

           //add the sendBlock to the send-list
           if (found==true) {
             if (SendNodeListCounter+1==SendNodelListCounterMax) exit(__LINE__,__FILE__,"Error: increase the value of SendNodelListCounterMax");
             SendNodeList[SendNodeListCounter++]=sendNode;
           }
         }
         */
//#endif


         //add the edge-neighbors to the send list
#if _MESH_DIMENSION_ == 3
         for (i=0;i<12*2;i++) if ((sendNode=recvNode->neibNodeEdge[i])!=NULL) if (sendNode->Thread==From) sendNode->nodeDescriptor.NodeProcessingFlag=_AMR_FALSE_;
         /*
           //Search the send-list if the node is already there
           for (found=false,k=0;k<SendNodeListCounter;k++) {
             if (sendNode==SendNodeList[k]) {
               found=true;
               break;
             }
           }

           //add the sendBlock to the send-list
           if (found==true) {
             if (SendNodeListCounter+1==SendNodelListCounterMax) exit(__LINE__,__FILE__,"Error: increase the value of SendNodelListCounterMax");
             SendNodeList[SendNodeListCounter++]=sendNode;
           }
         }
         */
#endif

//         recvNode=recvNode->nextNodeThisThread;
        }

/*
        for (k=0;k<SendNodeListCounter;k++) {
          sendNode=SendNodeList[k];
          GetAMRnodeID(nodeid,sendNode);

          pipe.send(_Next_Node_SIGNAL_);
          pipe.send((char*)(&nodeid),sizeof(nodeid));

           //send the data
           sendNode->block->SendNodeData(&pipe,DataSetTag);
         }*/


        //send the data
        for (recvNode=DomainBoundaryLayerNodesList[To];recvNode!=NULL;recvNode=recvNode->nextNodeThisThread) {


#if _MESH_DIMENSION_ == 2
#define _AMR_ParallelBlockDataExchange_SEND_CORNER_NODES_
#elif _MESH_DIMENSION_ == 3
#define _AMR_ParallelBlockDataExchange_SEND_CORNER_NODES_
#endif

#ifdef _AMR_ParallelBlockDataExchange_SEND_CORNER_NODES_

          for (i=0;i<(1<<_MESH_DIMENSION_);i++) if ((sendNode=recvNode->neibNodeCorner[i])!=NULL) if ((sendNode->Thread==From)&&(sendNode->nodeDescriptor.NodeProcessingFlag==_AMR_FALSE_)) {
            sendNode->nodeDescriptor.NodeProcessingFlag=_AMR_TRUE_;
            GetAMRnodeID(nodeid,sendNode);

            pipe.send(_Next_Node_SIGNAL_);
            pipe.send((char*)(&nodeid),sizeof(nodeid));

            //send the data
            sendNode->block->sendBoundaryLayerBlockData(&pipe);
          }
#endif

//#if _MESH_DIMENSION_ == 2
//#define _AMR_ParallelBlockDataExchange_SEND_FACES_
//#elif _MESH_DIMENSION_ == 3
//#define _AMR_ParallelBlockDataExchange_SEND_FACES_
//#endif

//#ifdef _AMR_ParallelBlockDataExchange_SEND_FACES_
         for (i=0;i<_MESH_DIMENSION_*(1<<_MESH_DIMENSION_);i++) if ((sendNode=recvNode->neibNodeFace[i])!=NULL) if ((sendNode->Thread==From)&&(sendNode->nodeDescriptor.NodeProcessingFlag==_AMR_FALSE_)) {
           sendNode->nodeDescriptor.NodeProcessingFlag=_AMR_TRUE_;
           GetAMRnodeID(nodeid,sendNode);

           pipe.send(_Next_Node_SIGNAL_);
           pipe.send((char*)(&nodeid),sizeof(nodeid));

           //send the data
           sendNode->block->sendBoundaryLayerBlockData(&pipe);
         }


//#endif

#if _MESH_DIMENSION_ == 3
         for (i=0;i<12*2;i++) if ((sendNode=recvNode->neibNodeEdge[i])!=NULL) if ((sendNode->Thread==From)&&(sendNode->nodeDescriptor.NodeProcessingFlag==_AMR_FALSE_)) {
           sendNode->nodeDescriptor.NodeProcessingFlag=_AMR_TRUE_;
           GetAMRnodeID(nodeid,sendNode);

           pipe.send(_Next_Node_SIGNAL_);
           pipe.send((char*)(&nodeid),sizeof(nodeid));

           //send the data
           sendNode->block->sendBoundaryLayerBlockData(&pipe);
         }
#endif
//         recvNode=recvNode->nextNodeThisThread;
        }


         pipe.send(_End_Communication_SIGNAL_);
         pipe.flush();
        //end the part of the sender
      }
      else if (ThisThread==To) {
        //the part of the receiver

        //redirect the recv's pipe buffers
        if (pipeLastRecvThread!=From) pipe.RedirectRecvBuffer(From);
        pipeLastRecvThread=From;

        pipe.recv(Signal,From);

        while (Signal==_Next_Node_SIGNAL_) {
          pipe.recv((char*)(&nodeid),sizeof(nodeid),From);
          recvNode=findAMRnodeWithID(nodeid);

          if (recvNode->block==NULL) {
            //chack the consistence ofthe mesh: 1. weather the 'recvNode' have neiborous block from 'ThisThread' and 2. weather 'recvNode' is in the DomainBoundaryLayerNodesList[From]
            bool found=false;
            cTreeNodeAMR<cBlockAMR> *searchNode;

            fprintf(DiagnospticMessageStream,"$PREFIX:Error: the node is not allocated:\n");

#ifdef _AMR_ParallelBlockDataExchange_SEND_CORNER_NODES_
            for (i=0;i<(1<<_MESH_DIMENSION_);i++) if (recvNode->neibNodeCorner[i]!=NULL) if (recvNode->neibNodeCorner[i]->Thread==ThisThread) {
              found=true;
              fprintf(DiagnospticMessageStream,"$PREFIX:'recvNode' has neibours on 'ThisThread': (file=%s, line=%i)\n",__FILE__,__LINE__);
            }
#endif

//#ifdef _AMR_ParallelBlockDataExchange_SEND_FACES_
           for (i=0;i<_MESH_DIMENSION_*(1<<_MESH_DIMENSION_);i++) if (recvNode->neibNodeFace[i]!=NULL) if (recvNode->neibNodeFace[i]->Thread==ThisThread) {
             found=true;
             fprintf(DiagnospticMessageStream,"$PREFIX:'recvNode' has neibours on 'ThisThread': (file=%s, line=%i)\n",__FILE__,__LINE__);
           }
//#endif

#if _MESH_DIMENSION_ == 3
           for (i=0;i<12*2;i++) if (recvNode->neibNodeEdge[i]!=NULL) if (recvNode->neibNodeEdge[i]->Thread==ThisThread) {
             found=true;
             fprintf(DiagnospticMessageStream,"$PREFIX:'recvNode' has neibours on 'ThisThread': (file=%s, line=%i)\n",__FILE__,__LINE__);
           }
#endif

           if (found==false) fprintf(DiagnospticMessageStream,"$PREFIX:'recvNode' doesn't have neibours at ThisThread (file=%s, line=%i)\n",__FILE__,__LINE__);

           //2. check if the 'recvNode' in the list 'DomainBoundaryLayerNodesList[From]'
           found=false;
           searchNode=DomainBoundaryLayerNodesList[From];

           while (searchNode!=NULL) {
             if (searchNode==recvNode) {
               found=true;
               fprintf(DiagnospticMessageStream,"$PREFIX:recvNode in the 'DomainBoundaryLayerNodesList[From]' (file=%s, line=%i)\n",__FILE__,__LINE__);
             }

             searchNode=searchNode->nextNodeThisThread;
           }

           if (found==false) fprintf(DiagnospticMessageStream,"$PREFIX:'recvNode' is not in 'DomainBoundaryLayerNodesList[From]' (file=%s, line=%i)\n",__FILE__,__LINE__);


            exit(__LINE__,__FILE__,"Error: the node is not allocated");
          }


//================  DEBUG ========================

          /*
if (ThisThread==2) if (pow(recvNode->xmin[0]+250.0,2)+pow(recvNode->xmin[1]+500.0,2)+pow(recvNode->xmin[0]+0.0,2)<0.00001) {
  *DiagnospticMessageStream << __LINE__ << std::endl;
}
*/

//================ END DEBUG =====================




          recvNode->block->recvBoundaryLayerBlockData(&pipe,From);
          pipe.recv(Signal,From);
        }

        //end the part of the receiver
      }

    }


    pipe.closeSend();
    pipe.closeRecv(pipeLastRecvThread);
  }

};










#endif 


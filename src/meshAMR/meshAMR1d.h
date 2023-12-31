//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf

//$Id$ 
//1D version of the AMR mesh

#ifndef _MESH_AMR_
#define _MESH_AMR_

//define the dimension of the mesh 
#define _MESH_DIMENSION_ 1

#include "meshAMRdef.h"




#include "meshAMRgeneric.h"

template <class cCornerNode,class cCenterNode,class cBlockAMR>
class cMeshAMR1d : public cMeshAMRgeneric<cCornerNode,cCenterNode,cBlockAMR> {
public:

  cMeshAMR1d() : cMeshAMRgeneric<cCornerNode,cCenterNode,cBlockAMR> () {
  } 

  

};


#endif

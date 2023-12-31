//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf

//$Id$ 
//2D version of the AMR mesh

#ifndef _MESH_AMR_
#define _MESH_AMR_

//define the dimension of the mesh 
#define _MESH_DIMENSION_ 3

#include "meshAMRdef.h"




#include "meshAMRgeneric.h"

template <class cCornerNode,class cCenterNode,class cBlockAMR> 
class cMeshAMR3d : public cMeshAMRgeneric<cCornerNode,cCenterNode,cBlockAMR> {
public:

_TARGET_HOST_ _TARGET_DEVICE_
  cMeshAMR3d() : cMeshAMRgeneric<cCornerNode,cCenterNode,cBlockAMR> () {
  }

};


#endif

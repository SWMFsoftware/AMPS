//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//define the global variables for the whole execution code

#include "mpi.h"

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "global.dfn"

#ifndef _GLOBAL_VARIABLES_
#define _GLOBAL_VARIABLES_

extern MPI_Comm MPI_GLOBAL_COMMUNICATOR;

#endif

//compilation mode
#define _COMPILATION_MODE_ _COMPILATION_MODE__MPI_

//using AVX instructions in calculations
#define _AVX_INSTRUCTIONS_USAGE_MODE_  _AVX_INSTRUCTIONS_USAGE_MODE__OFF_

//macros used for CUDA
#define _TARGET_GLOBAL_
#define _TARGET_HOST_ 
#define _TARGET_DEVICE_
#define _CUDA_MODE_ _OFF_

//inlcude settings of the general block
#include "../../.general.conf"


#endif


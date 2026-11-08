#ifndef AMPS_T7_SURFACE_EXOSPHERE_STUB_SPECFUNC_H
#define AMPS_T7_SURFACE_EXOSPHERE_STUB_SPECFUNC_H

// Exosphere_DefaultSurfaceDataStructure.h only needs the AMPS three-argument
// fatal hook from specfunc.h.  The test translation unit defines it to throw a
// harness failure, keeping this bounded buffer test independent of MPI/PIC.
#include <cstdio>

void exit(long int line, const char* file, const char* message);

#endif

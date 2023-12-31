//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//use the header when SPICE libraty is not available

#ifndef _SPICE_EMPTY_DEFINITIONS_
#define _SPICE_EMPTY_DEFINITIONS_

typedef double SpiceDouble;
typedef int SpiceInt;
typedef char SpiceChar;
typedef const char ConstSpiceChar;
typedef int SpiceBoolean;

inline void furnsh_c(ConstSpiceChar *file) {}
inline void utc2et_c(ConstSpiceChar *utcstr,SpiceDouble *et) {*et=0.0;}
inline void et2utc_c(SpiceDouble et, ConstSpiceChar *format, SpiceInt prec, SpiceInt lenout, SpiceChar *utcstr) {}

inline void spkezr_c(ConstSpiceChar *targ,SpiceDouble et,ConstSpiceChar *ref,ConstSpiceChar *abcorr,ConstSpiceChar *obs,SpiceDouble starg[6],SpiceDouble *lt) {
  *lt=0.0;
  for (int i=0;i<6;i++) starg[i]=1.0;
}

#endif

#ifndef _CONSTFIELDBC_MAIN_LIB_H_
#define _CONSTFIELDBC_MAIN_LIB_H_

#include <string>

#define _UNIFORM_MESH_ 1
#define _NONUNIFORM_MESH_ 2
#define _TEST_MESH_MODE_ _UNIFORM_MESH_

#include "cli.h"
#include "pic.h"
//#include "pic_units_normalization.h" 

// defined in main.cpp
extern TestConfig cfg;

// Convert any physical-unit inputs in cfg (SI) into the unit system expected by the
// field solver / particle routines. This must be called once after
// ConfigureTestFromArgsWithInput() and before any IC/setup routines.
void FinalizeConfigUnits(TestConfig& cfg);
extern double xmin[3];
extern double xmax[3];
extern int g_TestStencilOrder;


double BulletLocalResolution(double *x);
void InitGlobalParticleWeight_TargetPPC(const TestConfig& cfg);
void CleanParticles();
long int PrepopulateDomain();
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
void SetIC();
#endif

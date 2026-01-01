#ifndef _CONSTFIELDBC_MAIN_LIB_H_
#define _CONSTFIELDBC_MAIN_LIB_H_

#include <string>

#define _UNIFORM_MESH_ 1
#define _NONUNIFORM_MESH_ 2
#define _TEST_MESH_MODE_ _UNIFORM_MESH_

#include "cli.h"
#include "pic.h"

// defined in main.cpp
extern TestConfig cfg;
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




#ifndef _CONSTFIELDBC_MAIN_LIB_H_
#define _CONSTFIELDBC_MAIN_LIB_H_

#include <string>

#define _UNIFORM_MESH_ 1
#define _NONUNIFORM_MESH_ 2
#define _TEST_MESH_MODE_ _UNIFORM_MESH_

#include "cli.h"
#include "pic.h"
#include "../pic/units/pic_units_normalization.h"

// defined in main.cpp
extern TestConfig cfg;

// Convert any physical-unit inputs in cfg (SI) into the unit system expected by the
// field solver / particle routines. This must be called once after
// ConfigureTestFromArgsWithInput() and before any IC/setup routines.
picunits::Factors FinalizeConfigUnits(TestConfig& cfg);
// Global domain extents used by the test driver. These are filled in main.cpp
// either from legacy xmin/xmax defaults or from -L (domain size) option.
// Units: the mesh coordinate units for this test (same as used everywhere else).
extern double xmin[3];
extern double xmax[3];
extern int g_TestStencilOrder;


// Local surface resolution callback used by AMPS internal boundaries.
// The framework calls this to decide the surface-element target size near x.
// In this test we keep a uniform surface resolution, so BulletLocalResolution()
// returns a constant (see main_lib.cpp).
double BulletLocalResolution(double *x);
void InitGlobalParticleWeight_TargetPPC(const picunits::Factors& F,const TestConfig& cfg);
void CleanParticles();
long int PrepopulateDomain(int spec,picunits::Factors,const TestConfig& cfg);
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
void SetIC();

double EvaluateGlobalParticleWeightForTargetPPC(const picunits::Factors& F,const TestConfig& cfg); 
double EvaluateCFLTimeStepForSpecies(int spec, double CFL);

int InjectBoundaryParticles(const picunits::Factors& F,
                                 const TestConfig& cfg,
                                 int spec,
                                 double dt_no); 
//particle mover 
int MoverTestConstBC(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

// Species-dependent mover dispatch table. Populated from cfg (CLI / input file)
// on first use (lazy-init) and can also be explicitly re-initialized.
extern PIC::Mover::fSpeciesDependentParticleMover g_SpeciesParticleMoverTable[_TOTAL_SPECIES_NUMBER_];
void InitSpeciesParticleMoverTable(const TestConfig& cfg);

// ---------------------------------------------------------------------------
// Optional internal spherical boundary ("Enceladus" placeholder)
// ---------------------------------------------------------------------------
// When enabled, we register an internal sphere through the AMPS internal-boundary
// module (PIC::BC::InternalBoundary::Sphere). The intent is to treat the sphere
// as a solid obstacle inside the domain (a simple Enceladus stand-in).
//
// IMPORTANT ordering:
//   * Mesh must already exist: PIC::Mesh::mesh is created and xmin/xmax are final.
//   * Sphere must be registered before particles are created, so:
//       - PrepopulateDomain() can avoid putting particles inside the body.
//       - The mover sees the sphere as an internal boundary from the first step.
void InitInternalSphericalBoundary(const TestConfig& cfg);

// Utility used by particle initialization (PrepopulateDomain) to reject samples
// inside the internal sphere. If the sphere is disabled, this always returns false.
bool IsPointInsideInternalSphere(const double x[3]);

//==============================================================================
// Optional initialization of reduced (aligned) velocity state at particle birth
//------------------------------------------------------------------------------
// See the detailed discussion in bc.cpp. The key points are:
//   - When _USE_PARTICLE_V_PARALLEL_NORM_ is ON, particles have dedicated storage
//     for V_parallel and V_normal (in addition to the full 3D velocity vector).
//   - For guiding-center / gyrokinetic / magnetic-moment workflows it is useful
//     to also compute the magnetic moment at birth:
//         mu = (gamma^2 m v_perp^2) / (2|B|)
//   - This test driver uses a spatially uniform background magnetic field B0
//     (cfg.B0 in solver units), so we can compute these values in the generic
//     PIC::ParticleBuffer::InitiateParticle() callback without x/node.
//==============================================================================
namespace VparVnormMu {
  extern double gUniformB0_no[3];
  void InitParticle(PIC::ParticleBuffer::byte* ParticleDataStart); 
}
#endif

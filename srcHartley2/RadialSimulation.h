#ifndef _GENERIC_COMET_RADIAL_SIMULATION_H_
#define _GENERIC_COMET_RADIAL_SIMULATION_H_

#include <string>

template <class T> class cTreeNodeAMR;
namespace PIC { namespace Mesh { class cDataBlockAMR; } }

namespace Comet {
namespace RadialSimulation {

// These are numerical controls, not comet properties. CharacteristicSpeed is
// the maximum resolved transport speed used by the cell-crossing timestep
// criterion; it is indexed by the compiled AMPS species number.
extern bool Enabled;
extern double CourantNumber;
extern double CharacteristicSpeed[];
extern double InnerBoundaryTemperature;
// MaxIterations is a hard safety ceiling. SamplingLength defines the number
// of particle steps in each statistically independent convergence window.
extern long int MaxIterations;
extern long int SamplingLength;
// ConvergenceTolerance applies to the volume-weighted relative L1 change of
// each species' complete radial number-density profile. The run stops only
// after every species passes for RequiredConvergedWindows successive windows.
extern double ConvergenceTolerance;
extern int RequiredConvergedWindows;

bool Validate(std::string *error=NULL);

// Complete the AMPS initialization that follows radial mesh construction,
// install the spherical source and boundary callbacks, and advance the run.
void Initialize();
void Run();

// The following functions match AMPS callback signatures.  They are public so
// the core can store function pointers to them; case code should normally call
// only Validate, Initialize, and Run.
double LocalTimeStep(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
double LocalBlockInjectionRate(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
bool InjectionBlockIndicator(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
long int InjectBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
int ProcessDomainBoundary(long int ptr,double *position,double *velocity,
    int intersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

}
}

#endif

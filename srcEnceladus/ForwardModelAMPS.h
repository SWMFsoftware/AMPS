#ifndef SRCENCELADUS_FORWARD_MODEL_AMPS_H
#define SRCENCELADUS_FORWARD_MODEL_AMPS_H

#include "ForwardModel.h"
#include "pic.h"

#include <string>

namespace Enceladus {
namespace ForwardModel {
namespace AMPS {

// Parse and transactionally install the immutable runtime/source receipts.
// This must be called before PIC::InitMPI().
bool InstallFromInput(const std::string& inputFile, std::string* error);
bool IsInstalled();
const Configuration& InstalledConfiguration();
const SourceDistribution& InstalledDistribution();
const NumericalPlan& InstalledNumericalPlan();

// Called after the final distributed AMR blocks exist. It scans the actual
// allocated minimum cell width, computes the single global Courant step and
// species weights, validates the generated AMPS species table, and installs
// the generic split/merge selection.
void FinalizeNumericsAfterMeshAllocation();

// Initialize an explicit zero background for plasma_background_mode=none,
// load the configured cache, or invoke ICES only when explicitly selected.
void LoadCoupledPlasmaBackground();

// AMPS callback installed through PIC::BC::UserDefinedParticleInjectionFunction.
long int InjectSurfaceParticles();

// AMPS user-particle-processing callback. It is selected in enceladus.input
// so reaction choice remains runtime data while generic PIC owns traversal.
void ProcessParticle(long int particle, long int& firstParticleCell,
                     cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
int MoveParticle(long int particle, double timeStep,
                 cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

int AbsorbParticleAtSurface(int species, long int particle, double* position,
    double* velocity, double& remainingTime, void* node, void* sphere);
void TotalParticleAcceleration(double* acceleration, int species,
    long int particle, double* position, double* velocity,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
int TransformationIndicator(double* position, double* velocity, int species,
    long int particle, PIC::ParticleBuffer::byte* particleData, double& dt,
    bool& limitTimeStep, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
int TransformationProcessor(double* initialPosition, double* finalPosition,
    double* velocity, int& species, long int particle,
    PIC::ParticleBuffer::byte* particleData, double dt,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

double GlobalTimeStep(int species,
                      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
double PlannedSourceRate(int species);
void PrintReceipt();
void PrintFinalReceipt(long completedIterations);

}  // namespace AMPS
}  // namespace ForwardModel
}  // namespace Enceladus

#endif

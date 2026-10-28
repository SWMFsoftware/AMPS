#include "pic.h"
#include "CometChemistry.h"
#include "CometModel.h"
#include "CometPhysics.h"
#include "RadialMesh.h"
#include "RadialSimulation.h"

#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc,char **argv) {
  (void)argc;
  (void)argv;

#if DIM != 1
  std::fprintf(stderr,"Comet RadialMain requires DIM=1\n");
  return EXIT_FAILURE;
#else
  // AMPS input preprocessing has already compiled case values into the source
  // tree.  This call initializes MPI, molecular data, and core services that
  // the validators use; no runtime text input is parsed here.
  PIC::Init_BeforeParser();

  std::string validationError;
  if (!Comet::Model::Validate(&validationError,
      Comet::RadialSimulation::Enabled)) {
    std::fprintf(stderr,"Invalid comet model: %s\n",validationError.c_str());
    return EXIT_FAILURE;
  }
  if (!Comet::Chemistry::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid comet chemistry: %s\n",validationError.c_str());
    return EXIT_FAILURE;
  }
  if (!Comet::Physics::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid comet physics database: %s\n",
        validationError.c_str());
    return EXIT_FAILURE;
  }
  if (!Comet::RadialSimulation::Validate(&validationError)) {
    std::fprintf(stderr,"Invalid radial simulation: %s\n",
        validationError.c_str());
    return EXIT_FAILURE;
  }

  // The computational coordinate is cometocentric radius.  DIM1 spherical
  // symmetry in the AMR core supplies shell volumes and radial face areas.
  Comet::RadialModel::InitMeshProfile();

  double xmin[1]={Comet::RadialModel::NucleusRadius};
  double xmax[1]={Comet::RadialModel::GetDomainOuterRadius()};

  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,Comet::RadialModel::LocalResolution);

  // One rank constructs and serializes the AMR tree.  Other ranks read the
  // same tree before AMPS creates its load-distribution lists.
  const char *meshFile="comet.radial.mesh.msh";
  if (PIC::ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile(meshFile);
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread!=0) PIC::Mesh::mesh->readMeshFile(meshFile);

  PIC::Mesh::mesh->SetParallelLoadMeasure([](
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) -> double {return 1.0;});
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  Comet::RadialModel::PrintMeshDiagnostics();

  if (Comet::RadialSimulation::Enabled) {
    Comet::RadialSimulation::Initialize();
    Comet::RadialSimulation::Run();
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
#endif
}

#include "pic.h"
#include "RadialMesh.h"

#include <cstdio>
#include <cstdlib>

int main(int argc,char **argv) {
  (void)argc;
  (void)argv;

#if DIM != 1
  std::fprintf(stderr,"Hartley2 RadialMain requires DIM=1\n");
  return EXIT_FAILURE;
#else
  PIC::Init_BeforeParser();
  Comet::RadialModel::InitMeshProfile();

  double xmin[1]={Comet::RadialModel::NucleusRadius};
  double xmax[1]={Comet::RadialModel::GetDomainOuterRadius()};

  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,Comet::RadialModel::LocalResolution);

  const char *meshFile="hartley2.radial.mesh.msh";
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

  MPI_Finalize();
  return EXIT_SUCCESS;
#endif
}

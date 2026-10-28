#include "RadialMesh.h"

#include <cstdio>

#include "pic.h"

double Comet::RadialModel::NucleusRadius=-1.0;
double Comet::RadialModel::DomainRadialExtent=-1.0;
double Comet::RadialModel::DomainOuterRadius=-1.0;
double Comet::RadialModel::SurfaceResolution=-1.0;
int Comet::RadialModel::RequestedCellNumber=-1;

Comet::RadialModel::cLogarithmicResolutionProfile
    Comet::RadialModel::ResolutionProfile;

namespace {
long int CountCells(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) return 0;
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) return _BLOCK_CELLS_X_;

  long int result=0;
  for (int child=0;child<(1<<DIM);child++) result+=CountCells(node->downNode[child]);
  return result;
}

}

void Comet::RadialModel::InitMeshProfile() {
  const double outerRadius=GetDomainOuterRadius();
  const char *errorMessage=NULL;
  if (!ResolutionProfile.Init(NucleusRadius,outerRadius,SurfaceResolution,
      RequestedCellNumber,&errorMessage)) {
    exit(__LINE__,__FILE__,errorMessage);
  }
}

double Comet::RadialModel::GetDomainOuterRadius() {
  const bool extentIsSet=DomainRadialExtent>0.0;
  const bool outerRadiusIsSet=DomainOuterRadius>0.0;

  if (extentIsSet==outerRadiusIsSet) {
    exit(__LINE__,__FILE__,
        "set exactly one of DomainRadialExtent and DomainOuterRadius");
  }

  return extentIsSet ? NucleusRadius+DomainRadialExtent : DomainOuterRadius;
}

double Comet::RadialModel::LocalResolution(double *x) {
  return ResolutionProfile.GetResolution(x[0]);
}

long int Comet::RadialModel::CountAllocatedRadialCells() {
  return CountCells(PIC::Mesh::mesh->rootTree);
}

void Comet::RadialModel::PrintMeshDiagnostics() {
  if (PIC::ThisThread!=0) return;

  const double outerRadius=ResolutionProfile.GetOuterRadius();
  std::printf("$PREFIX: comet radial mesh: nucleus radius=%e m, outer radius=%e m\n",
      NucleusRadius,outerRadius);
  std::printf("$PREFIX: comet radial mesh: radial extent=%e m\n",
      outerRadius-NucleusRadius);
  std::printf("$PREFIX: comet radial mesh: requested surface resolution=%e m, ideal outer-cell resolution=%e m\n",
      SurfaceResolution,ResolutionProfile.GetOuterResolution());
  std::printf("$PREFIX: comet radial mesh: geometric cell-width ratio=%.12e\n",
      ResolutionProfile.GetExpansionRatio());
  std::printf("$PREFIX: comet radial mesh: requested ideal cells=%d, allocated dyadic AMR cells=%ld\n",
      RequestedCellNumber,
      CountAllocatedRadialCells());
}

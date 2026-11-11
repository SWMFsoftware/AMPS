#include "CometRadiativeCooling.h"

#include <cmath>
#include <cstring>
#include <sstream>

#include "pic.h"

int Comet::RadiativeCooling::SpeciesModel[PIC::nTotalSpecies]={0};

namespace {
bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
}

double CometocentricDistance(PIC::Mesh::cDataCenterNode *cell) {
  const double *position=cell->GetX();
  double radiusSquared=0.0;
  for (int idim=0;idim<DIM;idim++)
    radiusSquared+=position[idim]*position[idim];
  return std::sqrt(radiusSquared);
}
}

bool Comet::RadiativeCooling::Validate(std::string *error) {
  bool enabled=false;
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    const int model=SpeciesModel[species];
    if (model==static_cast<int>(CometaryPhysics::RadiativeCooling::Model::Disabled))
      continue;
    if (model!=static_cast<int>(
        CometaryPhysics::RadiativeCooling::Model::Gombosi1986Water)) {
      std::ostringstream out;
      out << "species " << species << " has an unknown radiative-cooling model "
          << model;
      return Fail(out.str(),error);
    }
    enabled=true;
  }

  if (!enabled) return true;

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ != _PIC_MODE_ON_
  return Fail("radiative cooling requires the AMPS internal-degrees-of-freedom model",error);
#elif _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_ != _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_
  return Fail("radiative cooling currently requires the AMPS Larsen-Borgnakke internal-energy model",error);
#else
  return true;
#endif
}

void Comet::RadiativeCooling::Apply() {
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_ && \
    _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_ == _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__LB_
  // Each rank owns the blocks in its parallel distribution list.  The cell
  // particle list is updated in place; no MPI exchange is required because
  // cooling does not move or create particles.
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=
      PIC::Mesh::mesh->ParallelNodesDistributionList[
          PIC::Mesh::mesh->ThisThread];
      node!=NULL;node=node->nextNodeThisThread) {
    PIC::Mesh::cDataBlockAMR *block=node->block;
    if (block==NULL) continue;

    for (int i=0;i<_BLOCK_CELLS_X_;i++)
      for (int j=0;j<_BLOCK_CELLS_Y_;j++)
        for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
          const int localCell=_getCenterNodeLocalNumber(i,j,k);
          PIC::Mesh::cDataCenterNode *cell=block->GetCenterNode(localCell);
          if (cell==NULL) continue;

          long int ptr=block->FirstCellParticleTable[
              i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
          if (ptr==-1) continue;

          // Evaluate a cell rate once per configured species.  Both accessors
          // read the completed sample.  Before the first complete window they
          // return zero, so cooling begins only when the requested prior-cycle
          // state exists.
          double lossRate[PIC::nTotalSpecies];
          for (int species=0;species<PIC::nTotalSpecies;species++) {
            lossRate[species]=0.0;
            const int modelCode=SpeciesModel[species];
            if (modelCode==static_cast<int>(
                CometaryPhysics::RadiativeCooling::Model::Disabled)) continue;

            const double rotationalTemperature=
                PIC::IDF::LB::GetCellRotTemp(species,cell);
            const double numberDensity=cell->GetNumberDensity(species);
            lossRate[species]=CometaryPhysics::RadiativeCooling::EnergyLossRate(
                static_cast<CometaryPhysics::RadiativeCooling::Model>(modelCode),
                rotationalTemperature,numberDensity,
                CometocentricDistance(cell));
          }

          for (;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
            PIC::ParticleBuffer::byte *data=
                PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            const int species=PIC::ParticleBuffer::GetI(data);
            if (!(lossRate[species]>0.0)) continue;

            const double timeStep=block->GetLocalTimeStep(species);
            double rotationalEnergy=PIC::IDF::LB::GetRotE(data)-
                lossRate[species]*timeStep;
            if (rotationalEnergy<0.0) rotationalEnergy=0.0;
            PIC::IDF::LB::SetRotE(rotationalEnergy,data);
          }
        }
  }
#endif
}

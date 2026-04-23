#include "Mode3D.h"
#include "ElectricField.h"

#include <cstdio>
#include <cmath>
#include <string>

#include "pic.h"
#include "../../interface/T96Interface.h"
#include "../../interface/T05Interface.h"
#include "../../interface/TA16Interface.h"

void amps_init_mesh();
void amps_init();

namespace Earth {
namespace Mode3D {

bool ParsedDomainActive=false;
double ParsedDomainMin[3]={0.0,0.0,0.0};
double ParsedDomainMax[3]={0.0,0.0,0.0};

namespace {

void ApplyParsedDomain(const EarthUtil::AmpsParam& prm) {
  ParsedDomainActive=true;
  ParsedDomainMin[0]=prm.domain.xMin*1000.0;
  ParsedDomainMin[1]=prm.domain.yMin*1000.0;
  ParsedDomainMin[2]=prm.domain.zMin*1000.0;
  ParsedDomainMax[0]=prm.domain.xMax*1000.0;
  ParsedDomainMax[1]=prm.domain.yMax*1000.0;
  ParsedDomainMax[2]=prm.domain.zMax*1000.0;
}

void ConfigureBackgroundFieldModel(const EarthUtil::AmpsParam& prm) {
  Earth::T96::active_flag=false;
  Earth::T05::active_flag=false;
  Earth::BackgroundMagneticFieldModelType=Earth::_undef;

  const std::string model=EarthUtil::ToUpper(prm.field.model);
  if (model=="T96") {
    Earth::BackgroundMagneticFieldModelType=Earth::_t96;
    Earth::T96::active_flag=true;
    Earth::T96::solar_wind_pressure=prm.field.pdyn_nPa*_NANO_;
    Earth::T96::dst=prm.field.dst_nT*_NANO_;
    Earth::T96::by=prm.field.imfBy_nT*_NANO_;
    Earth::T96::bz=prm.field.imfBz_nT*_NANO_;
    ::T96::SetSolarWindPressure(Earth::T96::solar_wind_pressure);
    ::T96::SetDST(Earth::T96::dst);
    ::T96::SetBYIMF(Earth::T96::by);
    ::T96::SetBZIMF(Earth::T96::bz);
    ::T96::Init(Exosphere::SimulationStartTimeString,Exosphere::SO_FRAME);
  }
  else if (model=="T05") {
    Earth::BackgroundMagneticFieldModelType=Earth::_t05;
    Earth::T05::active_flag=true;
    Earth::T05::solar_wind_pressure=prm.field.pdyn_nPa*_NANO_;
    Earth::T05::dst=prm.field.dst_nT*_NANO_;
    Earth::T05::by=prm.field.imfBy_nT*_NANO_;
    Earth::T05::bz=prm.field.imfBz_nT*_NANO_;
    for (int i=0;i<6;i++) Earth::T05::W[i]=prm.field.w[i];
    ::T05::SetSolarWindPressure(Earth::T05::solar_wind_pressure);
    ::T05::SetDST(Earth::T05::dst);
    ::T05::SetBXIMF(prm.field.imfBx_nT*_NANO_);
    ::T05::SetBYIMF(Earth::T05::by);
    ::T05::SetBZIMF(Earth::T05::bz);
    ::T05::SetW(Earth::T05::W[0],Earth::T05::W[1],Earth::T05::W[2],Earth::T05::W[3],Earth::T05::W[4],Earth::T05::W[5]);
    ::T05::Init(Exosphere::SimulationStartTimeString,Exosphere::SO_FRAME);
  }
  else if (model=="TA16") {
    // TA16 does not use BackgroundMagneticFieldModelType — it is driven
    // entirely through _PIC_COUPLER_MODE__TA16_ compile-time guards,
    // consistent with TA15 and T01.
    if (!prm.field.ta16CoeffFile.empty())
      ::TA16::SetCoeffFileName(prm.field.ta16CoeffFile);
    // TA16 PARMOD: [PDYN, SymHc, XIND, BYIMF, W1..W6]
    // SetSolarWindPressure/SetSymHc accept SI values (Pa / T); the _NANO_
    // factor converts from nPa / nT to SI, matching the T05 convention.
    ::TA16::SetSolarWindPressure(prm.field.pdyn_nPa*_NANO_);
    ::TA16::SetSymHc(prm.field.dst_nT*_NANO_);
    ::TA16::SetXIND(prm.field.xind);
    ::TA16::SetBYIMF(prm.field.imfBy_nT*_NANO_);
    ::TA16::Init(Exosphere::SimulationStartTimeString,Exosphere::SO_FRAME);
  }
}

void WriteTecplotMesh(const EarthUtil::AmpsParam& prm,const char* fnameBase) {
  // In MPI runs each rank writes its own Tecplot file. This avoids the need for
  // manual gather logic in this patch while still preserving all initialized
  // cell-center data from the distributed AMPS mesh.
  char fname[256];
  if (PIC::nTotalThreads>1) std::sprintf(fname,"%s.thread=%04d.dat",fnameBase,PIC::ThisThread);
  else std::sprintf(fname,"%s",fnameBase);

  FILE* fout=std::fopen(fname,"w");
  if (fout==nullptr) return;

  std::fprintf(fout,
    "TITLE=\"AMPS 3D Mesh Field Initialization\"\n"
    "VARIABLES=\"x_Re\",\"y_Re\",\"z_Re\",\"CellSize_Re\",\"Bx_nT\",\"By_nT\",\"Bz_nT\",\"Bmag_nT\",\"Ex_mVm\",\"Ey_mVm\",\"Ez_mVm\",\"Emag_mVm\",\"PhiE_kV\"\n"
    "ZONE T=\"AMPS_CELL_CENTERS\", F=POINT\n");

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread]; node!=NULL; node=node->nextNodeThisThread) {
    if (node->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_ || node->block==NULL) continue;

    for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
          const int nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k);
          PIC::Mesh::cDataCenterNode* center=node->block->GetCenterNode(nd);
          if (center==NULL) continue;

          double xCell[3];
          xCell[0]=node->xmin[0]+(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_*(0.5+i);
          xCell[1]=node->xmin[1]+(node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_*(0.5+j);
          xCell[2]=node->xmin[2]+(node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_*(0.5+k);

          double B[3],E[3];
          EvaluateBackgroundMagneticFieldSI(B,xCell,prm);
          EvaluateElectricFieldSI(E,xCell,prm);

          const double Bmag=std::sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
          const double Emag=std::sqrt(E[0]*E[0]+E[1]*E[1]+E[2]*E[2]);
          const double phiV=EvaluateElectricPotential_V(xCell,prm);

          std::fprintf(fout,
            "%e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            xCell[0]/_EARTH__RADIUS_,xCell[1]/_EARTH__RADIUS_,xCell[2]/_EARTH__RADIUS_,
            node->GetCharacteristicCellSize()/_EARTH__RADIUS_,
            B[0]/_NANO_,B[1]/_NANO_,B[2]/_NANO_,Bmag/_NANO_,
            E[0]/1.0e-3,E[1]/1.0e-3,E[2]/1.0e-3,Emag/1.0e-3,
            phiV/1000.0);
        }
      }
    }
  }

  std::fclose(fout);
}

} // namespace

int Run(const EarthUtil::AmpsParam& prm) {
  // Reuse the standard AMPS initialization path so that the 3D mode generates
  // the same tree, internal sphere, block allocation, and cell metrics as the
  // rest of the code rather than creating a separate standalone mesh.
  Earth::ModelMode=Earth::CutoffRigidityMode;
  ApplyParsedDomain(prm);

  PIC::InitMPI();
  Exosphere::Init_SPICE();
  amps_init_mesh();
  amps_init();

  ConfigureBackgroundFieldModel(prm);

  // The output filename is passed as a base name; in parallel runs each rank
  // appends its own suffix inside WriteTecplotMesh().
  WriteTecplotMesh(prm,"amps_3d_initialized_mesh.dat");
  PIC::Mesh::mesh->outputMeshDataTECPLOT("amps_3d_initialized.data.dat",0);
  return 0;
}

} // namespace Mode3D
} // namespace Earth

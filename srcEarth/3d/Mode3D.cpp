#include "Mode3D.h"
#include "ElectricField.h"
#include "CutoffRigidityMode3D.h"
#include "GlobalMagneticField.h"

#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>

#include "pic.h"

#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
#include "../../interface/T96Interface.h"
#include "../../interface/T05Interface.h"
#include "../../interface/TA16Interface.h"
#endif

void amps_init_mesh();
void amps_init();

// Sphere surface-mesh resolution parameters and the per-surface-element
// resolution function are defined in main_lib.cpp and shared with Mode3D.
extern int nZenithElements;
extern int nAzimuthalElements;
double localSphericalSurfaceResolution(double *x);

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
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  // In the live AMPS--SWMF coupling mode the background magnetic field is not
  // selected from the standalone Tsyganenko/TA16 wrappers.  SWMF supplies the
  // MHD state to AMPS through PIC::CPLR, and the field that should be used by
  // particle movers / diagnostics is exposed by the standard AMPS coupler
  // accessors:
  //
  //   PIC::CPLR::InitInterpolationStencil(...)
  //   PIC::CPLR::GetBackgroundMagneticField(...)
  //
  // Therefore this setup routine must be a no-op in SWMF builds: do not set
  // Earth::T96/Earth::T05 active flags and do not call ::T96::Init(),
  // ::T05::Init(), or ::TA16::Init().  Those model interfaces may not even be
  // linked in the coupled executable.
  (void)prm;
#else
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
#endif
}

// Traverse the full AMR tree and write B and E into every cell's data buffer,
// following the same ghost-cell-inclusive iteration pattern used by
// Earth::InitMagneticField in Earth.cpp.  Unlike that function, this version:
//   - uses EvaluateBackgroundMagneticFieldSI / EvaluateElectricFieldSI so all
//     models (T96, T05, TA16, DIPOLE) are handled uniformly, and
//   - initialises E from the configured electric-field model instead of
//     unconditionally writing zero.
void InitMeshFields(const EarthUtil::AmpsParam& prm,
                    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  // In SWMF-coupled builds the cell-centered B and E values are owned by the
  // live coupler data structures, not by the DATAFILE/MULTIFILE buffers that
  // this standalone initializer fills.  Overwriting DATAFILE fields here would
  // create a second, stale field source and would bypass the AMPS/SWMF access
  // pattern used elsewhere:
  //
  //   PIC::CPLR::GetBackgroundMagneticField(...)
  //   PIC::CPLR::GetBackgroundElectricField(...)
  //
  // Leave the mesh field buffers untouched; downstream field evaluation must
  // obtain the coupled fields through PIC::CPLR.
  (void)prm;
  (void)startNode;
  return;
#endif

  const int iMin=-_GHOST_CELLS_X_, iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_, jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_, kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->block!=NULL) {
      const int S=(kMax-kMin+1)*(jMax-jMin+1)*(iMax-iMin+1);

      for (int ii=0;ii<S;ii++) {
        int S1=ii;
        const int i=iMin+S1/((kMax-kMin+1)*(jMax-jMin+1));
        S1=S1%((kMax-kMin+1)*(jMax-jMin+1));
        const int j=jMin+S1/(kMax-kMin+1);
        const int k=kMin+S1%(kMax-kMin+1);

        const int nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k);
        PIC::Mesh::cDataCenterNode* CenterNode=startNode->block->GetCenterNode(nd);
        if (CenterNode==NULL) continue;

        char* offset=CenterNode->GetAssociatedDataBufferPointer()
                    +PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin
                    +PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset;

        double xCell[3];
        xCell[0]=startNode->xmin[0]+(startNode->xmax[0]-startNode->xmin[0])/_BLOCK_CELLS_X_*(0.5+i);
        xCell[1]=startNode->xmin[1]+(startNode->xmax[1]-startNode->xmin[1])/_BLOCK_CELLS_Y_*(0.5+j);
        xCell[2]=startNode->xmin[2]+(startNode->xmax[2]-startNode->xmin[2])/_BLOCK_CELLS_Z_*(0.5+k);

        double B[3],E[3];
        EvaluateBackgroundMagneticFieldSI(B,xCell,prm);
        EvaluateElectricFieldSI(E,xCell,prm);

//for (int i=0;i<3;i++) B[i]=xCell[i];

        for (int idim=0;idim<3;idim++) {
          if (PIC::CPLR::DATAFILE::Offset::MagneticField.active==true) {
            *((double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+idim*sizeof(double)))=B[idim];
          }
          if (PIC::CPLR::DATAFILE::Offset::ElectricField.active==true) {
            *((double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+idim*sizeof(double)))=E[idim];
          }
        }
      }
    }
  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* downNode;
    for (int i=0;i<(1<<DIM);i++) {
      if ((downNode=startNode->downNode[i])!=NULL) InitMeshFields(prm,downNode);
    }
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

// ---------------------------------------------------------------------------
// InitSphere — explicitly initialises and configures the inner Earth boundary
// sphere following the same pattern as the SphereInsideDomain block in
// main_lib.cpp::amps_init_mesh().
//
// amps_init_mesh() registers the sphere and sets Earth::Planet; we retrieve
// that handle here and re-apply every property so that the sphere setup is
// self-contained and auditable in the Mode3D flow.
//
// One important difference from main_lib.cpp: the cutoff-rigidity output
// callbacks are wired unconditionally.  In main_lib.cpp they are only set
// when (RigidityCalculationMode == Earth::_sphere &&
//       CutoffRigidity::SampleRigidityMode == true).
// Mode3D::Run() always operates in CutoffRigidityMode, so both conditions
// are implicitly satisfied and the guard can be dropped.
// ---------------------------------------------------------------------------
void InitSphere() {
  // amps_init_mesh() has already registered the sphere and assigned
  // Earth::Planet.  Retrieve the pointer and bail if it is somehow null.
  cInternalSphericalData* Sphere =
      static_cast<cInternalSphericalData*>(Earth::Planet);
  if (Sphere == nullptr) return;

  // ---- Geometry -----------------------------------------------------------
  // Sphere centred at the origin, radius = Earth radius (matches main_lib.cpp
  // where sx0={0,0,0} and rSphere=_EARTH__RADIUS_).
  double sx0[3] = {0.0, 0.0, 0.0};
  Sphere->SetSphereGeometricalParameters(sx0, _EARTH__RADIUS_);
  Sphere->Radius = _RADIUS_(_EARTH_);

  // ---- Surface mesh -------------------------------------------------------
  // Surface-mesh discretisation is shared with main_lib.cpp via the externs
  // declared at the top of this file.
  cInternalSphericalData::SetGeneralSurfaceMeshParameters(
      nZenithElements, nAzimuthalElements);

  // ---- Callbacks ----------------------------------------------------------
  // Particle–sphere interaction and injection (same values as main_lib.cpp).
  Sphere->ParticleSphereInteraction  = Earth::BC::ParticleSphereInteraction;
  Sphere->InjectionRate              = Exosphere::SourceProcesses::totalProductionRate;
  Sphere->InjectionBoundaryCondition = Exosphere::SourceProcesses::InjectionBoundaryModel;

  // Per-surface-element resolution function and face offset.
  Sphere->localResolution = localSphericalSurfaceResolution;
  Sphere->faceat          = 0;

  // ---- Diagnostic surface files -------------------------------------------
  // Re-emit the surface-mesh and initial surface-data files so that the output
  // on disk reflects the final Mode3D sphere configuration (geometry and
  // callbacks set above), not just the partial state left by amps_init_mesh().
  // The filenames match those written by main_lib.cpp::amps_init_mesh() so
  // that downstream post-processing scripts find the expected files.
  Sphere->PrintSurfaceMesh("Sphere.dat");
  Sphere->PrintSurfaceData("SpheraData.dat", 0);

  // ---- Cutoff-rigidity output callbacks -----------------------------------
  // In Mode3D, Run() always sets Earth::ModelMode = CutoffRigidityMode, so
  // the cutoff-rigidity output callbacks must always be connected.  Wire them
  // unconditionally here instead of repeating the conditional from
  // main_lib.cpp (which only fires when RigidityCalculationMode == _sphere &&
  // SampleRigidityMode == true).
  Earth::CutoffRigidity::AllocateCutoffRigidityTable();
  Sphere->PrintDataStateVector =
      Earth::CutoffRigidity::OutputDataFile::PrintDataStateVector;
  Sphere->PrintVariableList =
      Earth::CutoffRigidity::OutputDataFile::PrintVariableList;
}

} // namespace

int Run(const EarthUtil::AmpsParam& prm) {
  //============================================================================
  // Mode3D::Run — entry point for the standalone 3-D cutoff-rigidity workflow.
  //
  // Initialization sequence
  // -----------------------
  // 1. Mark execution mode and configure the physical domain bounds from prm.
  //    main_lib.cpp::amps_init_mesh() reads ParsedDomainMin/Max while building the
  //    AMR tree, and RunCutoffRigidity() later uses the same values for the outer
  //    escape box.  This keeps mesh geometry and trajectory classification aligned.
  //
  // 2. Use the normal AMPS MPI initialization, not independentDomainMode.
  //    Earlier standalone cutoff code called
  //
  //       PIC::InitMPI(/*independentDomainMode=*/true)
  //
  //    so that every MPI process built a private complete copy of the AMR mesh.
  //    That is no longer needed.  The current path initializes AMPS in the standard
  //    distributed-domain mode used elsewhere in the model.
  //
  // 3. Fill the selected standalone B/E field into the owner-rank DATAFILE cell
  //    buffers.  At this point only the normal distributed set of blocks is
  //    allocated on each rank.
  //
  // 4. Materialize a global read-only magnetic-field snapshot for cutoff tracing.
  //    This reuses the same data-management strategy developed for the SWMF-coupled
  //    cutoff path:
  //
  //       assign deterministic node->Temp_ID values,
  //       allocate missing leaf blocks on every rank,
  //       pack only owner-rank interior-cell B values,
  //       MPI_Allreduce the dense B cache and presence mask,
  //       populate all allocated interior and ghost cells from that cache.
  //
  //    After this step each MPI rank can trace its assigned cutoff trajectories
  //    through the whole AMR domain using the existing cell-centered interpolation
  //    code.  The ranks are no longer independent private AMPS universes; they are
  //    normal MPI ranks sharing one MPI_GLOBAL_COMMUNICATOR.
  //
  // 5. RunCutoffRigidity performs MPI × OpenMP work distribution:
  //    - static point distribution across MPI ranks;
  //    - OpenMP parallelism over points/directions inside each rank;
  //    - rank 0 gathers scalar Rc/Emin products and writes Tecplot output.
  //============================================================================

  Earth::ModelMode = Earth::CutoffRigidityMode;
  ApplyParsedDomain(prm);

  // ---- Standard distributed-domain initialization --------------------------
  //
  // Do not request independentDomainMode here.  We want the usual AMPS MPI domain
  // decomposition so owner blocks are distributed across ranks.  The global B-field
  // materialization below will later allocate nonlocal blocks and gather the owner
  // cell-centered magnetic field into them specifically for cutoff tracing.
  // -------------------------------------------------------------------------
  PIC::InitMPI();
  Exosphere::Init_SPICE();

  amps_init_mesh();   // build the replicated AMR tree topology and distributed blocks
  amps_init();        // data-offset registration and cutoff-mode particle settings

  //----------------------------------------------------------------------------
  // Sphere and field-model configuration
  //----------------------------------------------------------------------------
  InitSphere();
  ConfigureBackgroundFieldModel(prm);

  //----------------------------------------------------------------------------
  // Mesh field initialisation (diagnostic Tecplot output)
  //
  // Standalone/non-SWMF builds:
  //   Fill DATAFILE-style B/E buffers with the selected analytic/Tsyganenko field.
  //
  // SWMF-coupled builds:
  //   InitMeshFields() intentionally returns without writing DATAFILE buffers.
  //   The magnetic/electric fields are the live SWMF-coupled fields exposed by
  //   PIC::CPLR::GetBackgroundMagneticField/ElectricField.
  //----------------------------------------------------------------------------
  InitMeshFields(prm, PIC::Mesh::mesh->rootTree);

  // Diagnostic mesh output of the distributed owner/local field immediately after
  // initialization.  This is intentionally written BEFORE the global cutoff snapshot
  // allocates nonlocal blocks on every rank; at this point outputMeshDataTECPLOT sees
  // the normal AMPS parallel ownership layout and can assemble a clean file.
  if (prm.mode3d.outputInitializedFile) {
    PIC::Mesh::mesh->outputMeshDataTECPLOT("amps_3d_initialized.data.dat", 0);
  }

#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  // ---- Replicate the standalone magnetic field for cutoff tracing -----------
  //
  // InitMeshFields() just filled B on the owner/local DATAFILE cells that exist in
  // the normal distributed AMPS mesh.  The cutoff tracer on each rank can move a
  // particle through any AMR block, so before RunCutoffRigidity() starts we create
  // a read-only global B snapshot on every rank.  The helper below is the same
  // generic data-management operation used by the SWMF-coupled cutoff path; only
  // the cell-data offset differs.
  //
  // The electric field is not gathered here because the 3-D cutoff mover uses only
  // B for rigidity classification.  If a future cutoff formulation needs E, the
  // same helper pattern can be extended to another vector slot.
  if (!PIC::CPLR::DATAFILE::Offset::MagneticField.active) {
    exit(__LINE__,__FILE__,
         "[Mode3D] DATAFILE magnetic-field offset is inactive before standalone 3-D cutoff materialization.");
  }

  Earth::Mode3D::GlobalMagneticField::MaterializeCellCenteredMagneticFieldForCutoff(
      "Mode3D",
      Earth::Mode3D::GlobalMagneticField::DataFileMagneticFieldDataOffset(),
      true);
#endif

  //----------------------------------------------------------------------------
  // Cutoff rigidity computation (MPI_COMM_WORLD × OpenMP)
  // Each normal MPI rank processes its static point slice; OpenMP threads
  // within each rank parallelise over individual trajectory points.
  //----------------------------------------------------------------------------
  std::cout << "[Mode3D] Starting cutoff rigidity calculation...\n";
  std::cout.flush();

  const int status = RunCutoffRigidity(prm);

  std::cout << "[Mode3D] Cutoff rigidity calculation complete (status="
            << status << ").\n";
  std::cout.flush();

  return status;
}

} // namespace Mode3D
} // namespace Earth

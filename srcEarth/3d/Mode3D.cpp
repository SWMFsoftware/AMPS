#include "Mode3D.h"
#include "ElectricField.h"
#include "CutoffRigidityMode3D.h"
#include "DensityMode3D.h"
#include "GlobalMagneticField.h"

#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cctype>
#include <stdexcept>
#include <algorithm>

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

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

// User-defined Mode3D AMR mesh-resolution profile.  These globals are read by
// main_lib.cpp::localResolution() while AMPS builds the tree.  The inactive default
// preserves the historical hard-coded resolution function exactly.
bool MeshResolutionProfileActive=false;
double MeshResolutionEarth_m=0.0;
double MeshResolutionBoundary_m=0.0;
double MeshResolutionOuterRadius_Re=0.0;
double MeshResolutionExponent=1.0;
int MeshResolutionCoarseningCode=0; // 0=LINEAR, 1=LOG/EXPONENTIAL, 2=POWER, 3=CONSTANT

namespace {

//--------------------------------------------------------------------------------------
// Mode3D calculation-target helpers
//--------------------------------------------------------------------------------------
//
// The parser stores CALC_TARGET as one string for backward compatibility.  For Mode3D we
// intentionally accept both historical single-product tokens and compact combined tokens
// so one input file can request:
//   CALC_TARGET  CUTOFF_RIGIDITY                 -> cutoff only
//   CALC_TARGET  DENSITY_SPECTRUM                -> density + flux only
//   CALC_TARGET  CUTOFF_RIGIDITY+DENSITY_SPECTRUM -> both products from one field snapshot
//
// We use substring tests rather than a rigid enum so separators such as '+', ',', or
// whitespace (if preserved by an input reader) all work.  The solver still fails below if
// neither recognized product is requested.
static bool Mode3DTargetRequestsCutoff(const EarthUtil::AmpsParam& prm) {
  const std::string t = EarthUtil::ToUpper(prm.calc.target);
  return t.find("CUTOFF") != std::string::npos || t=="ALL" || t=="BOTH";
}

static bool Mode3DTargetRequestsDensityFlux(const EarthUtil::AmpsParam& prm) {
  const std::string t = EarthUtil::ToUpper(prm.calc.target);
  return t.find("DENSITY") != std::string::npos || t.find("FLUX") != std::string::npos ||
         t=="ALL" || t=="BOTH";
}


void ApplyParsedDomain(const EarthUtil::AmpsParam& prm) {
  ParsedDomainActive=true;
  ParsedDomainMin[0]=prm.domain.xMin*1000.0;
  ParsedDomainMin[1]=prm.domain.yMin*1000.0;
  ParsedDomainMin[2]=prm.domain.zMin*1000.0;
  ParsedDomainMax[0]=prm.domain.xMax*1000.0;
  ParsedDomainMax[1]=prm.domain.yMax*1000.0;
  ParsedDomainMax[2]=prm.domain.zMax*1000.0;
}

static double MaxAbsDomainFaceRadiusRe(const EarthUtil::AmpsParam& prm) {
  const double rEarth_km = _EARTH__RADIUS_ / 1000.0;
  double rmax_km = 0.0;
  rmax_km = std::max(rmax_km, std::fabs(prm.domain.xMin));
  rmax_km = std::max(rmax_km, std::fabs(prm.domain.xMax));
  rmax_km = std::max(rmax_km, std::fabs(prm.domain.yMin));
  rmax_km = std::max(rmax_km, std::fabs(prm.domain.yMax));
  rmax_km = std::max(rmax_km, std::fabs(prm.domain.zMin));
  rmax_km = std::max(rmax_km, std::fabs(prm.domain.zMax));
  return (rEarth_km > 0.0) ? (rmax_km / rEarth_km) : 0.0;
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
    ::T96::Init(prm.field.epoch.c_str(),Exosphere::SO_FRAME);
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
    ::T05::Init(prm.field.epoch.c_str(),Exosphere::SO_FRAME);
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
    ::TA16::Init(prm.field.epoch.c_str(),Exosphere::SO_FRAME);
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
      // In the original standalone cutoff implementation every MPI rank used
      // independentDomainMode=true, so the fact that a leaf had an allocated
      // block meant that the current rank was responsible for initializing it.
      //
      // The new standalone cutoff path deliberately uses the normal distributed
      // AMPS domain decomposition during field initialization and then builds a
      // replicated read-only B-field cache immediately before tracing.  After
      // that cache/materialization step, all used leaf blocks may remain
      // allocated on every MPI rank.  This becomes important for time-dependent
      // Tsyganenko runs: on the second and later snapshots InitMeshFields() is
      // called again while those replicated non-owner blocks already exist.
      //
      // Only the AMPS owner rank (startNode->Thread) is allowed to recompute the
      // cell-centered field here.  The global-field materialization helper packs
      // owner interior cells, MPI-reduces them, and then overwrites every
      // replicated block/ghost cell from that authoritative cache.  If this
      // routine also initialized non-owner replicated blocks, it would waste a
      // large amount of CPU time and could leave confusing stale values in
      // buffers before the materialization pass.
      if (startNode->Thread!=PIC::ThisThread) return;

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


// ---------------------------------------------------------------------------
// Mode3D time-snapshot helpers
// ---------------------------------------------------------------------------
// The AMPS_PARAM parser already knows how to load a Tsyganenko driver table from
// either
//
//   #TEMPORAL
//   TS_INPUT_MODE  FILE
//   TS_INPUT_FILE  ...
//
// or the more compact
//
//   #BACKGROUND_FIELD
//   DRIVER_FILE    ...
//
// and stores that table in prm.temporal.driverTable.  The helpers below are the
// runtime counterpart for the standalone mesh-backed 3-D cutoff path: they build
// the list of magnetic-field snapshots, interpolate the driver table at each
// snapshot epoch, update prm.field, and create a unique output suffix so that a
// multi-snapshot run does not overwrite its own Tecplot files.
// ---------------------------------------------------------------------------

bool Mode3DTimeSeriesRequested(const EarthUtil::AmpsParam& prm) {
  return EarthUtil::ToUpper(prm.temporal.mode) == "TIME_SERIES";
}

std::string Mode3DSanitizeSuffixToken(const std::string& in) {
  // Output suffixes become part of filenames such as
  //   cutoff_3d_shells_snapshot_000003_2024_05_10T12_15_00.dat
  // Keep only portable filename characters.  This intentionally converts ':',
  // '-', '.', whitespace, and any SPICE-specific decorations to underscores.
  std::string out;
  out.reserve(in.size());
  for (std::size_t i=0;i<in.size();++i) {
    const unsigned char c = static_cast<unsigned char>(in[i]);
    if (std::isalnum(c)) out.push_back(static_cast<char>(c));
    else if (in[i]=='T' || in[i]=='Z') out.push_back(in[i]);
    else out.push_back('_');
  }

  // Avoid extremely long filenames if a user supplies a verbose time string.
  if (out.size()>48) out.resize(48);
  return out;
}

std::string Mode3DSnapshotSuffix(int iSnapshot,const std::string& epochUTC) {
  std::ostringstream os;
  os << "_snapshot_" << std::setw(6) << std::setfill('0') << iSnapshot
     << "_" << Mode3DSanitizeSuffixToken(epochUTC);
  return os.str();
}

#ifndef _NO_SPICE_CALLS_
double Mode3DEpochToEtOrExit(const std::string& epochUTC,const char* context) {
  SpiceDouble et = 0.0;
  try {
    str2et_c(epochUTC.c_str(),&et);
  }
  catch (...) {
    std::ostringstream msg;
    msg << "[Mode3D] Cannot convert UTC epoch '" << epochUTC
        << "' to SPICE ET while " << context << ".";
    exit(__LINE__,__FILE__,msg.str().c_str());
  }
  return static_cast<double>(et);
}

std::string Mode3DEtToUtc(double et) {
  // ISOC gives a compact sortable UTC label: YYYY-MM-DDTHH:MM:SS.sss.
  // Three decimals are enough for field-update cadences specified in minutes.
  char utc[96];
  et2utc_c(static_cast<SpiceDouble>(et),"ISOC",3,static_cast<SpiceInt>(sizeof(utc)),utc);
  return std::string(utc);
}
#endif

std::vector<std::string> Mode3DBuildSnapshotEpochs(const EarthUtil::AmpsParam& prm) {
  std::vector<std::string> epochs;

  // Legacy/snapshot behavior: one field realization at #BACKGROUND_FIELD/EPOCH.
  // A driver table may still be present in this mode; it is sampled once at this
  // epoch by Mode3DBuildSnapshotParam() below.
  if (!Mode3DTimeSeriesRequested(prm)) {
    epochs.push_back(prm.field.epoch);
    return epochs;
  }

  // TIME_SERIES is meaningful only when an external driver table is loaded.  The
  // parser loads that table from #TEMPORAL/TS_INPUT_FILE or #BACKGROUND_FIELD/DRIVER_FILE.
  if (prm.temporal.driverTable.empty()) {
    exit(__LINE__,__FILE__,
         "[Mode3D] TEMPORAL_MODE=TIME_SERIES requires a Tsyganenko driver table: "
         "set #TEMPORAL TS_INPUT_MODE=FILE and TS_INPUT_FILE, or set "
         "#BACKGROUND_FIELD DRIVER_FILE.");
  }

  if (prm.temporal.eventStart.empty() || prm.temporal.eventEnd.empty()) {
    exit(__LINE__,__FILE__,
         "[Mode3D] TEMPORAL_MODE=TIME_SERIES requires EVENT_START and EVENT_END "
         "in the #TEMPORAL section.");
  }

  if (!(prm.temporal.fieldUpdateDt_min>0.0)) {
    exit(__LINE__,__FILE__,
         "[Mode3D] FIELD_UPDATE_DT must be positive for TEMPORAL_MODE=TIME_SERIES.");
  }

#ifdef _NO_SPICE_CALLS_
  exit(__LINE__,__FILE__,
       "[Mode3D] TEMPORAL_MODE=TIME_SERIES requires SPICE to convert EVENT_START/END "
       "and driver-table timestamps. Rebuild without _NO_SPICE_CALLS_.");
#else
  const double etStart = Mode3DEpochToEtOrExit(prm.temporal.eventStart,"building the Mode3D snapshot list");
  const double etEnd   = Mode3DEpochToEtOrExit(prm.temporal.eventEnd,  "building the Mode3D snapshot list");
  const double dt      = prm.temporal.fieldUpdateDt_min * 60.0;

  if (etEnd < etStart) {
    exit(__LINE__,__FILE__,
         "[Mode3D] EVENT_END must not be earlier than EVENT_START for TIME_SERIES.");
  }

  // Include both endpoints when they lie on the requested cadence.  The small
  // tolerance protects against roundoff in floating-point ET arithmetic.
  const double tol = std::max(1.0e-6,1.0e-9*dt);
  for (double et=etStart; et<=etEnd+tol; et+=dt) {
    epochs.push_back(Mode3DEtToUtc(std::min(et,etEnd)));
  }

  if (epochs.empty()) epochs.push_back(Mode3DEtToUtc(etStart));
  return epochs;
#endif
}

EarthUtil::AmpsParam Mode3DBuildSnapshotParam(const EarthUtil::AmpsParam& base,
                                              const std::string& epochUTC) {
  EarthUtil::AmpsParam snap = base;
  snap.field.epoch = epochUTC;

  // If a driver table is loaded, interpolate it at this snapshot epoch and copy
  // the resulting Pdyn/Dst/IMF/W/G/BZ/XIND values into snap.field.  Downstream
  // code then treats this snapshot exactly like a normal one-time input-file
  // configuration; no field evaluator needs to know whether the parameters came
  // from a static block or from a time table.
  if (!snap.temporal.driverTable.empty()) {
#ifdef _NO_SPICE_CALLS_
    exit(__LINE__,__FILE__,
         "[Mode3D] A Tsyganenko driver table was loaded, but this build has "
         "_NO_SPICE_CALLS_ defined. Cannot convert snapshot UTC to ET.");
#else
    const double et = Mode3DEpochToEtOrExit(epochUTC,"sampling the Mode3D Tsyganenko driver table");
    const EarthUtil::TsDriverRecord rec = snap.temporal.driverTable.Lookup(et);
    EarthUtil::TsDriverTable::ApplyToField(rec,snap.field);
#endif
  }

  return snap;
}

void Mode3DPrepareMagneticFieldSnapshot(const EarthUtil::AmpsParam& snap,
                                         const std::string& suffix,
                                         bool verbose) {
  // Reinitialize the selected standalone field backend for this snapshot.  For
  // Tsyganenko models this updates the model parameters and Geopack epoch before
  // any cell-centered field values are evaluated.  For DIPOLE this is cheap and
  // simply resets the internal Mode3D model flags.
  ConfigureBackgroundFieldModel(snap);

  // Fill the DATAFILE B/E buffers on all currently allocated blocks.  During the
  // first snapshot this is the normal distributed AMPS mesh; after the global-B
  // materialization step below it also includes replicated nonlocal blocks.  That
  // is safe: the subsequent gather packs only node->Thread==PIC::ThisThread owner
  // cells as authoritative values, so nonowner recalculations are ignored.
  InitMeshFields(snap,PIC::Mesh::mesh->rootTree);

  // Optional diagnostic output of the initialized cell-centered field.  Append the
  // same snapshot suffix used by cutoff outputs so multiple snapshots do not
  // overwrite each other.
  if (snap.mode3d.outputInitializedFile) {
    const std::string fname = "amps_3d_initialized" + suffix + ".data.dat";
    PIC::Mesh::mesh->outputMeshDataTECPLOT(fname.c_str(),0);
  }

#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  // Replicate the owner-rank cell-centered B field to every MPI process for this
  // snapshot.  The same reusable helper is called for every snapshot, so the code
  // path is identical to the single-snapshot standalone cutoff and to the SWMF
  // coupled data-management strategy.
  if (!PIC::CPLR::DATAFILE::Offset::MagneticField.active) {
    exit(__LINE__,__FILE__,
         "[Mode3D] DATAFILE magnetic-field offset is inactive before standalone "
         "3-D cutoff materialization.");
  }

  Earth::Mode3D::GlobalMagneticField::MaterializeCellCenteredMagneticFieldForCutoff(
      "Mode3D",
      Earth::Mode3D::GlobalMagneticField::DataFileMagneticFieldDataOffset(),
      verbose);
#else
  (void)verbose;
#endif
}

} // namespace

void ConfigureMeshResolutionProfile(const EarthUtil::AmpsParam& prm) {
  MeshResolutionProfileActive = prm.mode3d.meshResolutionProfileActive;

  if (!MeshResolutionProfileActive) {
    MeshResolutionEarth_m=0.0;
    MeshResolutionBoundary_m=0.0;
    MeshResolutionOuterRadius_Re=0.0;
    MeshResolutionExponent=1.0;
    MeshResolutionCoarseningCode=0;
    return;
  }

  MeshResolutionEarth_m = prm.mode3d.meshResolutionEarth_km * 1000.0;
  MeshResolutionBoundary_m = prm.mode3d.meshResolutionBoundary_km * 1000.0;
  MeshResolutionExponent = prm.mode3d.meshResolutionExponent;

  if (prm.mode3d.meshResolutionOuterRadius_km > 0.0) {
    MeshResolutionOuterRadius_Re = (prm.mode3d.meshResolutionOuterRadius_km * 1000.0) / _EARTH__RADIUS_;
  }
  else {
    MeshResolutionOuterRadius_Re = MaxAbsDomainFaceRadiusRe(prm);
  }
  if (!(MeshResolutionOuterRadius_Re > 1.0)) MeshResolutionOuterRadius_Re = 1.01;

  const std::string c = EarthUtil::ToUpper(prm.mode3d.meshResolutionCoarsening);
  if (c=="LOG" || c=="EXP" || c=="EXPONENTIAL" || c=="GEOMETRIC" || c=="LOGARITHMIC")
    MeshResolutionCoarseningCode=1;
  else if (c=="POWER" || c=="POW" || c=="EXPONENT" || c=="POLYNOMIAL")
    MeshResolutionCoarseningCode=2;
  else if (c=="CONSTANT" || c=="CONST" || c=="UNIFORM")
    MeshResolutionCoarseningCode=3;
  else
    MeshResolutionCoarseningCode=0;

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3D mesh] User-defined AMR resolution profile: "
              << "res_earth=" << MeshResolutionEarth_m/_EARTH__RADIUS_ << " Re, "
              << "res_boundary=" << MeshResolutionBoundary_m/_EARTH__RADIUS_ << " Re, "
              << "r_boundary=" << MeshResolutionOuterRadius_Re << " Re, "
              << "coarsening=" << c << ", exponent=" << MeshResolutionExponent
              << std::endl;
  }
}

double ConfiguredMeshResolutionSI(const double *x) {
  if (!MeshResolutionProfileActive) return -1.0;

  double r2=0.0;
  for (int idim=0; idim<DIM; ++idim) r2 += x[idim]*x[idim];
  const double r = std::sqrt(r2);
  const double r_Re = r / _EARTH__RADIUS_;

  // Preserve the legacy AMR behavior inside the loss sphere.  The Mode3D mesh
  // exists in a Cartesian box that includes the inactive Earth interior, but
  // particle tracing stops at the inner boundary and does not need a finely
  // resolved interior volume.  Applying the user requested near-Earth surface
  // resolution throughout r < Re can create a very large number of unnecessary
  // AMR blocks and may make validation cases such as C5 fail before they ever
  // test the mesh-interpolated field.  The old localResolution() path used a
  // coarse cell size for r < 0.98 Re; keep that protection when the optional
  // user-defined radial profile is active.
  if (r < 0.98*_EARTH__RADIUS_) return _EARTH__RADIUS_;

  double t = 0.0;
  if (MeshResolutionOuterRadius_Re > 1.0) {
    t = (r_Re - 1.0) / (MeshResolutionOuterRadius_Re - 1.0);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
  }

  const double a = MeshResolutionEarth_m;
  const double b = MeshResolutionBoundary_m;

  if (!(a > 0.0) || !(b > 0.0)) return -1.0;

  switch (MeshResolutionCoarseningCode) {
    case 1: { // LOG/EXPONENTIAL/GEOMETRIC interpolation in resolution.
      return std::exp(std::log(a) + t*(std::log(b)-std::log(a)));
    }
    case 2: { // POWER profile in normalized altitude.
      const double tt = std::pow(t,MeshResolutionExponent);
      return a + (b-a)*tt;
    }
    case 3: { // CONSTANT/UNIFORM resolution.
      return a;
    }
    case 0:
    default: { // LINEAR profile.
      return a + (b-a)*t;
    }
  }
}

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
  ConfigureMeshResolutionProfile(prm);
  Exosphere::Init_SPICE();

  amps_init_mesh();   // build the replicated AMR tree topology and distributed blocks
  amps_init();        // data-offset registration and cutoff-mode particle settings

  //----------------------------------------------------------------------------
  // Sphere and static mesh geometry
  //----------------------------------------------------------------------------
  InitSphere();

  // Build the list of magnetic-field snapshots for this run.
  //
  //   SNAPSHOT mode (default): one entry, prm.field.epoch.
  //   TIME_SERIES mode       : EVENT_START..EVENT_END in steps of FIELD_UPDATE_DT.
  //
  // The parser has already loaded the Tsyganenko driver file, if requested, into
  // prm.temporal.driverTable.  For every epoch below we interpolate that table into
  // a temporary AmpsParam copy and then reuse the normal single-snapshot field
  // initialization + global-B materialization + cutoff solver.
  const std::vector<std::string> snapshotEpochs = Mode3DBuildSnapshotEpochs(prm);

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3D] Magnetic-field snapshot count: "
              << snapshotEpochs.size();
    if (Mode3DTimeSeriesRequested(prm)) {
      std::cout << " (TIME_SERIES, FIELD_UPDATE_DT="
                << prm.temporal.fieldUpdateDt_min << " min)";
    }
    std::cout << "\n";
    std::cout.flush();
  }

  int finalStatus = 0;

  for (std::size_t iSnapshot=0; iSnapshot<snapshotEpochs.size(); ++iSnapshot) {
    // Convert the base run parameters into this snapshot's effective magnetic
    // configuration.  If a driver table is present, this is where Pdyn/Dst/IMF/W
    // are interpolated at snapshotEpochs[iSnapshot].
    EarthUtil::AmpsParam snap = Mode3DBuildSnapshotParam(prm,snapshotEpochs[iSnapshot]);

    const std::string suffix = (snapshotEpochs.size()>1 || Mode3DTimeSeriesRequested(prm))
        ? Mode3DSnapshotSuffix(static_cast<int>(iSnapshot),snap.field.epoch)
        : std::string("");

    // Tell the cutoff writer to append the snapshot suffix before ".dat".  This
    // prevents outputs from later snapshots from overwriting earlier snapshots.
    SetCutoffOutputFileSuffix(suffix);
    SetDensityOutputFileSuffix(suffix);

    if (PIC::ThisThread == 0) {
      std::cout << "[Mode3D] Preparing magnetic-field snapshot "
                << (iSnapshot+1) << "/" << snapshotEpochs.size()
                << ": epoch=" << snap.field.epoch;
      if (!suffix.empty()) std::cout << ", suffix=" << suffix;
      std::cout << "\n";
      std::cout.flush();
    }

    // Reuse the existing single-snapshot data-management path for every time.
    // This call reinitializes the selected Tsyganenko/dipole field, writes it into
    // the mesh DATAFILE buffers, and materializes a global read-only B snapshot on
    // every MPI process for the cutoff tracer.
    Mode3DPrepareMagneticFieldSnapshot(snap,suffix,/*verbose=*/true);

    //------------------------------------------------------------------------
    // Requested physics products for this snapshot
    //------------------------------------------------------------------------
    // Both products below use the same prepared mesh-field snapshot.  Running them
    // consecutively here guarantees that, when the input requests
    // CUTOFF_RIGIDITY+DENSITY_SPECTRUM, cutoff, directional maps, density, spectra,
    // and flux are all derived from the identical magnetic-field state and carry the
    // same snapshot suffix in their file names.
    //------------------------------------------------------------------------
    const bool doCutoff      = Mode3DTargetRequestsCutoff(snap);
    const bool doDensityFlux = Mode3DTargetRequestsDensityFlux(snap);

    if (!doCutoff && !doDensityFlux) {
      std::ostringstream msg;
      msg << "Unsupported CALC_TARGET for -mode 3d: '" << snap.calc.target
          << "'. Supported targets include CUTOFF_RIGIDITY, DENSITY_SPECTRUM, "
          << "and CUTOFF_RIGIDITY+DENSITY_SPECTRUM.";
      throw std::runtime_error(msg.str());
    }

    if (doCutoff) {
      if (PIC::ThisThread == 0) {
        std::cout << "[Mode3D] Starting cutoff rigidity calculation for snapshot "
                  << (iSnapshot+1) << "/" << snapshotEpochs.size() << "...\n";
        std::cout.flush();
      }

      const int status = RunCutoffRigidity(snap,/*showProgressBar=*/true);
      if (status != 0) finalStatus = status;

      if (PIC::ThisThread == 0) {
        std::cout << "[Mode3D] Cutoff rigidity calculation complete for snapshot "
                  << (iSnapshot+1) << "/" << snapshotEpochs.size()
                  << " (status=" << status << ").\n";
        std::cout.flush();
      }
    }

    if (doDensityFlux) {
      if (PIC::ThisThread == 0) {
        std::cout << "[Mode3D] Starting density/flux calculation for snapshot "
                  << (iSnapshot+1) << "/" << snapshotEpochs.size() << "...\n";
        std::cout.flush();
      }

      const int status = RunDensityAndFlux(snap);
      if (status != 0) finalStatus = status;

      if (PIC::ThisThread == 0) {
        std::cout << "[Mode3D] Density/flux calculation complete for snapshot "
                  << (iSnapshot+1) << "/" << snapshotEpochs.size()
                  << " (status=" << status << ").\n";
        std::cout.flush();
      }
    }
  }

  // Reset the global suffix so a future in-process caller starts from the legacy
  // filename convention unless it explicitly installs another suffix.
  SetCutoffOutputFileSuffix("");
  SetDensityOutputFileSuffix("");

  return finalStatus;
}

} // namespace Mode3D
} // namespace Earth

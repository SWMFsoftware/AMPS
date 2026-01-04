/*
 * main.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: fougere and vtenishe
 */

/***************************************************************************************
 * Test driver for PIC with/without particles and constant field boundary conditions
 *
 * OVERVIEW
 * --------
 * This test can run in four modes controlled via the command line:
 *
 *   1) -particles
 *        • PIC with particles enabled.
 *        • Field boundary conditions are held constant.
 *        • Particle boundary condition is *specular reflection* at domain faces.
 *        • Initial fields: by default E = 0, B = 0 unless you also set a background.
 *
 *   2) -no-particles
 *        • Field-only (no particles).
 *        • Field boundary conditions are held constant.
 *        • Initial fields default to E = 0, B = 0 (unless -B/-E are provided).
 *
 *   3) -B [Bx By Bz]
 *        • Field-only (implies -no-particles).
 *        • Initialize B to a uniform vector; E = 0.
 *        • If the three numbers are omitted, defaults to B = (0, 1, 0).
 *
 *   4) -E [Ex Ey Ez]
 *        • Field-only (implies -no-particles).
 *        • Initialize E to a uniform vector; B = 0.
 *        • If the three numbers are omitted, defaults to E = (1, 0, 0).
 *
 * Additionally:
 *   -stencil-order=N
 *        • Sets the finite-difference stencil order used by divergence/curl/Poisson
 *          operators in this test (typical choices: 2,4,6,8). Default is 2.
 *
 * BOUNDARY CONDITIONS
 * -------------------
 * • Fields: “constant BC” — the field values at the domain boundary are held fixed to
 *   the values chosen by the mode (e.g., the uniform E or B you set). In practice this
 *   can be enforced via your project’s boundary callback or by reapplying the uniform
 *   values to halo layers after each step in field-only runs.
 * • Particles: “specular reflection” — when a particle hits a domain face, its velocity
 *   component normal to that face is flipped (tangential components are preserved).
 *
 * HOW SetIC() WORKS NOW
 * ---------------------
 * SetIC(cfg) zeros all fields, then applies the requested initial condition:
 *   • -particles: leaves E and B at your default (often zeros) unless your test adds a
 *     background; particles are (pre)populated later as usual.
 *   • -no-particles: leaves both fields zero.
 *   • -B: sets cell-centered B uniformly to (Bx,By,Bz), keeps corner E = 0.
 *   • -E: sets corner E uniformly to (Ex,Ey,Ez), keeps cell-centered B = 0.
 *
 * EXAMPLES
 * --------
 *   # PIC with particles, specular reflection, 4th-order stencil
 *   ./amps_test -particles -stencil-order=4
 *
 *   # Field-only, uniform B=(0,1,0), default 2nd-order
 *   ./amps_test -B
 *
 *   # Field-only, E=(0.2,0,0), 8th-order stencil
 *   ./amps_test -E 0.2 0.0 0.0 -stencil-order 8
 *
 *   # Field-only, all-zero fields, 6th-order stencil
 *   ./amps_test -no-particles -stencil-order=6
 *
 ***************************************************************************************/


#include "pic.h"
#include "constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"
#include "meshAMRgeneric.h"

#include "../../srcInterface/LinearSystemCornerNode.h"
#include "linear_solver_wrapper_c.h"

//#include "PeriodicBCTest.dfn"

#if _CUDA_MODE_ == _ON_
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

//for lapenta mover

#include "pic.h"
#include "Exosphere.dfn"
#include "Exosphere.h"



#include "main_lib.h"
#include "pic_units_normalization.h"

// -----------------------------------------------------------------------------
// InitGlobalParticleWeight_TargetPPC()
//   Choose a global particle weight so that PrepopulateDomain() injects roughly
//   cfg.target_ppc macro-particles per cell *per species* for the uniform
//   solar-wind-like IC.
//
//   PrepopulateDomain() injects particles with:
//       npart ≈ NumberDensity * CellVolume / ParticleWeight
//   where npart is per species (ions and electrons are injected separately).
//
//   We estimate a representative cell volume using the MIN cell volume across
//   all MPI ranks (and skip periodic halo blocks) to avoid under-resolving the
//   smallest cells when AMR/nonuniform meshes are used.
//
//   NOTE: This routine uses the same normalization logic as PrepopulateDomain()
//         (rho_conv and the definition of NumberDensity_ref) to ensure the ppc
//         target matches the injection formula used by the test.
// -----------------------------------------------------------------------------
void InitGlobalParticleWeight_TargetPPC(const TestConfig& cfg) {
  if (cfg.mode != TestConfig::Mode::WithParticles) return;

  const double Nppc_target = cfg.target_ppc;
  if (Nppc_target <= 0.0) return;

  // Species indices used by this test: ion=0, electron=1
  const int ionSpec=0, electronSpec=1;

  // Mass normalization consistent with PrepopulateDomain() in this driver.
  const double ionMass      = PIC::MolecularData::GetMass(ionSpec)/_AMU_;
  const double electronMass = PIC::MolecularData::GetMass(electronSpec)/_AMU_;

  // Ensure BlockTable is up to date.
  PIC::DomainBlockDecomposition::UpdateBlockTable();

  // Conservative representative cell volume: min cell volume over all ranks.
  double localMinCellVolume = 1.0e100;
  const int nBlock[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  for (int nLocalNode=0; nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks; nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->Thread!=PIC::ThisThread) continue;

    // In periodic mode, skip the periodic halo blocks (blocks with missing neighbors).
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      for (int iface=0; iface<6; iface++) {
        if (node->GetNeibFace(iface,0,0,PIC::Mesh::mesh)==NULL) { BoundaryBlock=true; break; }
      }
      if (BoundaryBlock==true) continue;
    }

    double CellVolume=1.0;
    for (int idim=0; idim<3; idim++) {
      const double dx=(node->xmax[idim]-node->xmin[idim])/(double)nBlock[idim];
      CellVolume*=dx;
    }
    if (CellVolume < localMinCellVolume) localMinCellVolume = CellVolume;
  }

  double globalMinCellVolume=0.0;
  MPI_Allreduce(&localMinCellVolume,&globalMinCellVolume,1,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);

  // If for some reason all local blocks were skipped (should not happen in normal runs),
  // fall back to a safe positive volume to avoid division by zero.
  if (!(globalMinCellVolume>0.0) || globalMinCellVolume>1.0e90) {
    globalMinCellVolume = 1.0;
  }

  // Reference number density derived from the configured mass density.
  // IMPORTANT: In this test driver, cfg.sw_* fields are finalized (including
  // any SI->normalized conversions) by FinalizeConfigUnits(cfg) in main.cpp.
  // Therefore, the injection formula uses cfg.sw_rho0/cfg.sw_p0 directly, with
  // no additional ad-hoc scaling.
  const double rho = cfg.sw_rho0;
  const double NumberDensity_ref = rho / (ionMass + electronMass);

  const double ParticleWeight_ref = NumberDensity_ref * globalMinCellVolume / Nppc_target;

  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(ionSpec,ParticleWeight_ref);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(electronSpec,ParticleWeight_ref);

  if (PIC::ThisThread==0) {
    std::printf("[ConstFieldBC] ParticleWeight set for ~%.1f ppc/spec: weight=%e, n_ref=%e, Vcell(min)=%e\n",
                Nppc_target, ParticleWeight_ref, NumberDensity_ref, globalMinCellVolume);
  }
}




//------------------------------------------------------------------------
void FinalizeConfigUnits(TestConfig& cfg) {
  // This routine converts *physical* (SI) user inputs (if provided) into the
  // unit system expected by the compiled ECSIM field solver.
  //
  // The AMPS/ECSIM build can be configured either to store fields in SI units
  // (_PIC_FIELD_SOLVER_INPUT_UNIT_SI_) or in a normalized "code-unit" system
  // (_PIC_FIELD_SOLVER_INPUT_UNIT_NORM_). For the latter, we use the
  // normalization in pic_units_normalization.h:
  //   - choose normalization scales (lSI, uSI, mSI)
  //   - compute factors
  //   - map SI -> normalized for {rho, p, v, B, E}
  //
  // IMPORTANT:
  //   * This is a *driver-level* convenience. The test logic (SetIC,
  //     PrepopulateDomain) expects cfg.B0/E0 and cfg.sw_* already expressed in
  //     the solver's input units after this routine runs.
  //   * For the solar-wind option, the convective field is E = -u x B.

  // Build normalization factors once (used only if the solver expects NORM units).
  // NOTE: picunits::build() takes a NormScalesSI struct (see pic_units_normalization.h).
  picunits::NormScalesSI norm_scales{cfg.units_lSI_m, cfg.units_uSI_mps, cfg.units_mSI_kg};
  picunits::Factors F = picunits::build(norm_scales);

  auto si_to_solver_B = [&](const double B_T[3], double B_out[3]) {
#if _PIC_FIELD_SOLVER_INPUT_UNIT_ == _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_
    picunits::si2no_B3(B_T, B_out, F);
#else
    // SI build: store Tesla directly
    for (int d=0; d<3; d++) B_out[d] = B_T[d];
#endif
  };

  auto si_to_solver_E = [&](const double E_Vm[3], double E_out[3]) {
#if _PIC_FIELD_SOLVER_INPUT_UNIT_ == _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_
    picunits::si2no_E3(E_Vm, E_out, F);
#else
    // SI build: store V/m directly
    for (int d=0; d<3; d++) E_out[d] = E_Vm[d];
#endif
  };

  auto si_to_solver_v = [&](const double u_mps[3], double u_out[3]) {
#if _PIC_FIELD_SOLVER_INPUT_UNIT_ == _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_
    picunits::si2no_v3(u_mps, u_out, F);
#else
    for (int d=0; d<3; d++) u_out[d] = u_mps[d];
#endif
  };

  auto si_to_solver_rho = [&](double rho_kgm3) -> double {
#if _PIC_FIELD_SOLVER_INPUT_UNIT_ == _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_
    return picunits::si2no_rho(rho_kgm3, F);
#else
    return rho_kgm3;
#endif
  };

  auto si_to_solver_p = [&](double p_Pa) -> double {
#if _PIC_FIELD_SOLVER_INPUT_UNIT_ == _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_
    return picunits::si2no_P(p_Pa, F);
#else
    return p_Pa;
#endif
  };

  // ---------------- Background fields (non-solar-wind mode too) ----------------
  if (cfg.userB_SI) {
    si_to_solver_B(cfg.B0_SI_T, cfg.B0);
    cfg.userB = true;
  }

  if (cfg.userE_SI) {
    si_to_solver_E(cfg.E0_SI_Vm, cfg.E0);
    cfg.userE = true;
  }

  // ---------------- Solar wind: convert physical IC to solver units ----------------
  if (cfg.mode == TestConfig::Mode::WithParticles) {
    // --- Velocity ---
    bool have_u_SI = false;
    double u_SI[3] = {0.0,0.0,0.0};
    if (cfg.sw_has_u_kms) {
      for (int d=0; d<3; d++) u_SI[d] = cfg.sw_u_kms[d] * 1.0e3; // km/s -> m/s
      have_u_SI = true;
    }
    else if (cfg.sw_has_u_mps) {
      for (int d=0; d<3; d++) u_SI[d] = cfg.sw_u_mps[d];
      have_u_SI = true;
    }

    if (have_u_SI) {
      si_to_solver_v(u_SI, cfg.sw_u0);
    }

    // --- Magnetic field ---
    // If the user provided SW B in nT, store it also as background B in SI so it
    // gets converted above (cfg.userB_SI).
    if (cfg.sw_has_BnT) {
      for (int d=0; d<3; d++) cfg.B0_SI_T[d] = cfg.sw_BnT[d] * 1.0e-9; // nT -> Tesla
      cfg.userB_SI = true;
      si_to_solver_B(cfg.B0_SI_T, cfg.B0);
      cfg.userB = true;
    }

    // --- Density & pressure ---
    // If n and T are provided in physical units, we build SI rho,p and convert.
    if (cfg.sw_has_ncm3) {
      const double kB_SI = 1.380649e-23; // J/K
      const double n_m3 = cfg.sw_n_cm3 * 1.0e6; // cm^-3 -> m^-3

      // NOTE: we keep the legacy assumption of one ion species (species=0) and
      //       quasi-neutral plasma with Ti=Te=T.
      const double mi_kg = PIC::MolecularData::GetMass(0);
      const double me_kg = PIC::MolecularData::GetMass(PIC::nTotalSpecies-1);

      const double rho_SI = n_m3 * (mi_kg + me_kg); // kg/m^3
      cfg.sw_rho0 = si_to_solver_rho(rho_SI);

      if (cfg.sw_has_TK) {
        const double T = cfg.sw_TK;
        const double p_SI = n_m3 * kB_SI * (T + T); // Pi+Pe with Ti=Te=T
        cfg.sw_p0 = si_to_solver_p(p_SI);
      }
    }

    // --- Electric field ---
    // Explicit physical E overrides E=-u×B.
    if (cfg.sw_has_EmVm || cfg.sw_has_EVm) {
      double E_SI[3] = {0.0,0.0,0.0};
      if (cfg.sw_has_EmVm) {
        for (int d=0; d<3; d++) E_SI[d] = cfg.sw_E_mVm[d] * 1.0e-3; // mV/m -> V/m
      }
      else {
        for (int d=0; d<3; d++) E_SI[d] = cfg.sw_E_Vm[d];
      }

      for (int d=0; d<3; d++) cfg.E0_SI_Vm[d] = E_SI[d];
      cfg.userE_SI = true;
      si_to_solver_E(cfg.E0_SI_Vm, cfg.E0);
      cfg.userE = true;
    }
    else {
      // If sw-evxb is unset, default to computing E=-u×B when particles are enabled.
      const bool do_evxb = (cfg.sw_evxb == 1) || (cfg.sw_evxb < 0);

      if (do_evxb && cfg.userB && (cfg.sw_has_u_kms || cfg.sw_has_u_mps || (cfg.sw_u0[0]!=0.0 || cfg.sw_u0[1]!=0.0 || cfg.sw_u0[2]!=0.0))) {
        // Compute in solver units (works for both SI and normalized builds as long as
        // u and B are in consistent units):  E = -u × B
        const double ux = cfg.sw_u0[0], uy = cfg.sw_u0[1], uz = cfg.sw_u0[2];
        const double Bx = cfg.B0[0],    By = cfg.B0[1],    Bz = cfg.B0[2];

        cfg.E0[0] = -(uy*Bz - uz*By);
        cfg.E0[1] = -(uz*Bx - ux*Bz);
        cfg.E0[2] = -(ux*By - uy*Bx);
        cfg.userE = true;
      }
    }
  }

  // Optional sanity print (rank 0 only) to help validate unit conversions.
  {
    int do_print = 1;
#ifdef MPI_ON
    int mpi_init = 0;
    MPI_Initialized(&mpi_init);
    if (mpi_init) {
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      do_print = (rank == 0);
    }
#endif
    if (do_print) {
#if _PIC_FIELD_SOLVER_INPUT_UNIT_ == _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_
      printf("[ConstFieldBC] ECSIM input units: NORM (code units). Using normalization scales: lSI=%g m, uSI=%g m/s, mSI=%g kg\n",
             cfg.units_lSI_m, cfg.units_uSI_mps, cfg.units_mSI_kg);
#else
      printf("[ConstFieldBC] ECSIM input units: SI. Background fields interpreted as SI (Tesla, V/m).\n");
#endif
    }
  }
}

// Uniform setters used by SetIC()
// NOTE: These use the same buffer/offset accessors that your legacy SetIC() uses.
void SetUniformCornerE(const double E0[3]) {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node; node = node->nextNodeThisThread) {
    if (!node->block) continue;

    char *offset;
    for (int k=0; k<=_BLOCK_CELLS_Z_; ++k)
    for (int j=0; j<=_BLOCK_CELLS_Y_; ++j)
    for (int i=0; i<=_BLOCK_CELLS_X_; ++i) {
      auto* cn = node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(i,j,k));
      if (!cn) continue;
      offset = cn->GetAssociatedDataBufferPointer() + PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

      ((double*)(offset + CurrentEOffset))[ExOffsetIndex] = E0[0];
      ((double*)(offset + CurrentEOffset))[EyOffsetIndex] = E0[1];
      ((double*)(offset + CurrentEOffset))[EzOffsetIndex] = E0[2];

      ((double*)(offset + OffsetE_HalfTimeStep))[ExOffsetIndex] = E0[0];
      ((double*)(offset + OffsetE_HalfTimeStep))[EyOffsetIndex] = E0[1];
      ((double*)(offset + OffsetE_HalfTimeStep))[EzOffsetIndex] = E0[2];
    }
  }
}

void SetUniformCenterB(const double B0[3]) {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node; node = node->nextNodeThisThread) {
    if (!node->block) continue;

    char *offset;
    for (int k=0; k<_BLOCK_CELLS_Z_; ++k)
    for (int j=0; j<_BLOCK_CELLS_Y_; ++j)
    for (int i=0; i<_BLOCK_CELLS_X_; ++i) {
      auto* cc = node->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
      if (!cc) continue;
      offset = cc->GetAssociatedDataBufferPointer() + PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

      ((double*)(offset + CurrentBOffset))[BxOffsetIndex] = B0[0];
      ((double*)(offset + CurrentBOffset))[ByOffsetIndex] = B0[1];
      ((double*)(offset + CurrentBOffset))[BzOffsetIndex] = B0[2];

      ((double*)(offset + PrevBOffset))[BxOffsetIndex] = B0[0];
      ((double*)(offset + PrevBOffset))[ByOffsetIndex] = B0[1];
      ((double*)(offset + PrevBOffset))[BzOffsetIndex] = B0[2];
    }
  }
}

// New SetIC based on CLI
void SetIC() {
  const double Z[3] = {0.0,0.0,0.0};
  // Zero fields
  SetUniformCornerE(Z);
  SetUniformCenterB(Z);

  // Apply chosen mode
  switch (cfg.mode) {
    case TestConfig::Mode::FieldOnlyB:
      SetUniformCenterB(cfg.B0);
      break;
    case TestConfig::Mode::FieldOnlyE:
      SetUniformCornerE(cfg.E0);
      break;
    case TestConfig::Mode::WithParticles:
      // Particle mode: keep fields uniform/constant if provided via -B/-E, otherwise leave as zero.
      if (cfg.userB) SetUniformCenterB(cfg.B0);
      if (cfg.userE) SetUniformCornerE(cfg.E0);
      break;
    case TestConfig::Mode::NoParticles:
      // Field-only run with no particles: honor user-specified backgrounds if provided.
      if (cfg.userB) SetUniformCenterB(cfg.B0);
      if (cfg.userE) SetUniformCornerE(cfg.E0);
      break;
    default:
      break;
  }

  PIC::Mesh::mesh->ParallelBlockDataExchange();
}

void CleanParticles(){
  
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

  for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) if (node->block!=NULL) {
   
     long int *  FirstCellParticleTable=node->block->FirstCellParticleTable;
     if (FirstCellParticleTable==NULL) continue;
     for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
	 for (int i=0;i<_BLOCK_CELLS_X_;i++) {
	   long int * ptr=FirstCellParticleTable+(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k));
	   while ((*ptr)!=-1) PIC::ParticleBuffer::DeleteParticle(*ptr,*ptr);
	 }   
       }
     }
     
  }

}


long int PrepopulateDomain() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  int iCell,jCell,kCell;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  long int nd,nGlobalInjectedParticles,nLocalInjectedParticles=0;
  double Velocity[3];
  /*
  //local copy of the block's cells
  int cellListLength=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread]->block->GetCenterNodeListLength();
  PIC::Mesh::cDataCenterNode *cellList[cellListLength];
  */
  //particle ejection parameters
  double ParticleWeight;//beta=PIC::MolecularData::GetMass(spec)/(2*Kbol*Temperature);
  double waveNumber[3]={0.0,0.0,0.0};
  double lambda=32.0;
 
  waveNumber[0]=2*Pi/lambda;

  double *ParticleDataTable=NULL,*ParticleDataTable_dev=NULL;
  int ParticleDataTableIndex=0,ParticleDataTableLength=0;

  
  int nBlock[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  //the boundaries of the block and middle point of the cell
  double *xminBlock,*xmaxBlock;
  double v[3],anpart;
  int npart;
  char * offset=NULL;
  int ionSpec=0, electronSpec=1;
  double ionMass = PIC::MolecularData::GetMass(ionSpec)/_AMU_;
  double electronMass = PIC::MolecularData::GetMass(electronSpec)/_AMU_;

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0,PIC::Mesh::mesh)==NULL) {
	  //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
	  BoundaryBlock=true;
	  break;
	}
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;


    // }

    // PIC::Mesh::cDataCenterNode *cellList[cellListLength];
  
    //memcpy(cellList,node->block->GetCenterNodeList(),cellListLength*sizeof(PIC::Mesh::cDataCenterNode*));

    xminBlock=node->xmin,xmaxBlock=node->xmax;
    double dx[3];
    double CellVolume=1;
    for (int idim=0;idim<3;idim++) {
      dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nBlock[idim];
      CellVolume *= dx[idim];
    }
    //particle stat weight
#ifndef _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
#error ERROR: _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_ is used but not defined
#endif
#ifndef _SIMULATION_PARTICLE_WEIGHT_MODE_
#error ERROR: _SIMULATION_PARTICLE_WEIGHT_MODE_ is used but not defined
#endif

    //assume ion and electron have the same particle weight
    #if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[ionSpec];
    #else
    ParticleWeight=node->block->GetLocalParticleWeight(ionSpec);
    #endif

//    double *ParticleDataTable=NULL,*ParticleDataTable_dev=NULL;
//    int ParticleDataTableIndex=0,ParticleDataTableLength=0;

    for (kCell=0;kCell<nBlock[2];kCell++) for (jCell=0;jCell<nBlock[1];jCell++) for (iCell=0;iCell<nBlock[0];iCell++) {
	  //      nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(iCell,jCell,kCell);

      // cell=cellList[nd];
      //  xMiddle=cell->GetX();
      //offset = cell->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
	  int ind[3]={iCell,jCell,kCell};
	  double x[3];
	  for (int idim=0; idim<3; idim++) x[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];          // --- Uniform solar-wind-like plasma (no waves) ---
          // The solar-wind IC is specified through cfg.sw_* and has already been
          // converted into the ECSIM unit system by FinalizeConfigUnits(cfg).
          // Here we use those values directly to build Maxwellian particle samples.
          const double rho = cfg.sw_rho0;
          const double p   = cfg.sw_p0;

          const double NumberDensity = rho/(ionMass+electronMass);
          const double rho_i = NumberDensity*ionMass;
          const double rho_e = NumberDensity*electronMass;

          // Split total scalar pressure between ions and electrons (isotropic).
          const double pi = 0.5*p;
          const double pe = 0.5*p;

          double ionBulkVelocity[3]      = {cfg.sw_u0[0], cfg.sw_u0[1], cfg.sw_u0[2]};
          double electronBulkVelocity[3] = {cfg.sw_u0[0], cfg.sw_u0[1], cfg.sw_u0[2]};

          // Component thermal speeds for isotropic Maxwellians: <v_x'^2> = p_species / rho_species
          const double uth_i = sqrt(pi / rho_i);
          const double uth_e = sqrt(pe / rho_e);


          //inject particles into the cell
          anpart=NumberDensity*CellVolume/ParticleWeight;
          //std::cout<<"CellLoc:"<<x[0]<<" "<<x[1]<<" "<<x[2]<<" NumberDensity: "<<NumberDensity<<"cell volume: "<<CellVolume<<"anpart: "<<anpart<<std::endl;
          npart=(int)(anpart);
          if (cfg.sw_use_rounding && (rnd() < anpart - npart)) npart++;
          nLocalInjectedParticles+=npart*2;
          //std::cout<<"need to inject npart: "<<npart<<std::endl;
          
          #if _CUDA_MODE_ == _ON_
          if (ParticleDataTableLength<npart) {
            if (ParticleDataTable!=NULL) {
              delete [] ParticleDataTable;
              cudaFree(ParticleDataTable_dev);
            }

            ParticleDataTable=new double [9*npart];
            cudaMalloc(&ParticleDataTable_dev,9*npart*sizeof(double));
          }
          
          ParticleDataTableIndex=0;
          #endif 

          while (npart-->0) {
            double xPar[3];
            xPar[0]=x[0]+dx[0]*(rnd()-0.5);
            xPar[1]=x[1]+dx[1]*(rnd()-0.5);
            
            // xPar[0]=x[0];
            // xPar[1]=x[1];
            xPar[2]=x[2];

            
            double electronVelocity[3],ionVelocity[3];
                        for (int idim=0; idim<3; idim++) {
              // Box-Muller: Gaussian(0,1) * uth + bulk
              const double g1 = sqrt(-2.0 * log(1.0 - 0.999999999 * rnd())) * cos(2.0*Pi*rnd());
              const double g2 = sqrt(-2.0 * log(1.0 - 0.999999999 * rnd())) * cos(2.0*Pi*rnd());
              electronVelocity[idim] = uth_e * g1 + electronBulkVelocity[idim];
              ionVelocity[idim]      = uth_i * g2 + ionBulkVelocity[idim];
            }
            
            /*  
            for (int idim=0;idim<3;idim++) {
              //in this test case B field is in y-direction
              double ElectronTemp= idim!=1?kTemp_perp/electronMass:kTemp_par/electronMass; 
              double IonTemp= idim!=1?kTemp_perp/ionMass:kTemp_par/ionMass; 
              
              electronVelocity[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())*(2*ElectronTemp))+electronBulkVelocity[idim];
              ionVelocity[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())*(2*IonTemp))+ionBulkVelocity[idim];       
              }
            */      
            //initiate the new particle
            
            #if _CUDA_MODE_ == _OFF_ 
            PIC::ParticleBuffer::InitiateParticle(xPar, electronVelocity,NULL,&electronSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            PIC::ParticleBuffer::InitiateParticle(xPar, ionVelocity,NULL,&ionSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            #else 
           
            memcpy(ParticleDataTable+0+9*ParticleDataTableIndex,xPar,3*sizeof(double));
            memcpy(ParticleDataTable+3+9*ParticleDataTableIndex,electronVelocity,3*sizeof(double));
            memcpy(ParticleDataTable+6+9*ParticleDataTableIndex,ionVelocity,3*sizeof(double));

            ParticleDataTableIndex++;
            #endif
            
          }
      //end of the particle injection block
      //std::cout<<"finished injecting npart: "<<npart<<std::endl;
      
     #if _CUDA_MODE_ == _ON_ 
          auto InitParticle = [=] _TARGET_DEVICE_ (double *ParticleDataTable, int ParticleDataTableIndex,int electronSpec, int ionSpec, void *node) {
            int id=blockIdx.x*blockDim.x+threadIdx.x;
            int increment=gridDim.x*blockDim.x;

            for (int i=id;i<ParticleDataTableIndex;i+=increment) {
              PIC::ParticleBuffer::InitiateParticle(ParticleDataTable+0+9*i,ParticleDataTable+3+9*i,NULL,&electronSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
              PIC::ParticleBuffer::InitiateParticle(ParticleDataTable+0+9*i,ParticleDataTable+6+9*i,NULL,&ionSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            }
          };


          cudaMemcpy(ParticleDataTable_dev,ParticleDataTable,9*ParticleDataTableIndex*sizeof(double),cudaMemcpyHostToDevice);

          kernel_5<<<1,1>>>(InitParticle,ParticleDataTable_dev,ParticleDataTableIndex,electronSpec,ionSpec,node);
          cudaDeviceSynchronize();
      #endif
      
      
        }
        }

    #if _CUDA_MODE_ == _ON_
    delete [] ParticleDataTable;
    cudaFree(ParticleDataTable_dev);
    #endif

  MPI_Allreduce(&nLocalInjectedParticles,&nGlobalInjectedParticles,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  printf("particles prepopulated!\n");
  return nGlobalInjectedParticles;
}



double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    double CellSize;
    double CharacteristicSpeed;
    double dt;


    CellSize=startNode->GetCharacteristicCellSize();
    //return 0.3*CellSize/CharacteristicSpeed;

    //return 0.05;
    return 1;
}


double BulletLocalResolution(double *x) {                                                                                           
  double dist = xmax[0]-xmin[0];

#ifndef _UNIFORM_MESH_
#error ERROR: _UNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_UNIFORM_MESH_  
  double res = 3;
#endif

#ifndef _NONUNIFORM_MESH_
#error ERROR: _NONUNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  double highRes = dist/32.0, lowRes= dist/2.0;     
  double res =(5-1)/dist*(x[0]-xmin[0])+1;  
#endif

  res=sqrt(3)+0.1;
  return res;
}
                       



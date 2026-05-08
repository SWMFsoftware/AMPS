//======================================================================================
// Mode3DForward.cpp
//======================================================================================
//
// Full PIC-backed 3D forward particle transport solver.
// See Mode3DForward.h for the design overview and initialisation sequence.
//
//======================================================================================

#include "Mode3DForward.h"
#include "Density3D.h"
#include "BoundaryDistribution.h"

#include "../boundary/spectrum.h"
#include "../3d/ElectricField.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <memory>
#include <algorithm>
#include <functional>

#include "pic.h"
#include "../../interface/T96Interface.h"
#include "../../interface/T05Interface.h"
#include "../../interface/TA16Interface.h"

void amps_init_mesh();
void amps_init();
void amps_time_step();

// Sphere surface-mesh resolution (defined in main_lib.cpp; shared with Mode3D)
extern int nZenithElements;
extern int nAzimuthalElements;
double localSphericalSurfaceResolution(double* x);

namespace Earth {
namespace Mode3DForward {

// ============================================================================
//  Constants / physical constants
// ============================================================================
// NOTE: Pi is already defined as a #define macro by the AMPS framework
// (build/general/constants.h). Do NOT redefine it here.
static constexpr double MeV_in_J     = 1.602176634e-13;  // 1 MeV in Joules
// Fraction of cell-size used for time-step calculation
static constexpr double DtCellFrac   = 0.25;

// ============================================================================
//  Module-level state (set in Run, read by injection/sphere callbacks)
// ============================================================================
static double sParticleWeight      = 1.0;   // physical particles per sim particle
static double sDt                  = 1.0;   // time step [s]
static int    sSpecies             = 0;     // AMPS species index for injection

// Pointer to the inner absorbing sphere (set by InitAbsorptionSphere)
static cInternalSphericalData* sAbsorptionSphere = nullptr;

// ============================================================================
//  Helper: configure background field (identical to Mode3D path)
// ============================================================================
static void ConfigureBackgroundFieldModel(const EarthUtil::AmpsParam& prm) {
  Earth::T96::active_flag = false;
  Earth::T05::active_flag = false;
  Earth::BackgroundMagneticFieldModelType = Earth::_undef;

  const std::string model = EarthUtil::ToUpper(prm.field.model);
  if (model == "T96") {
    Earth::BackgroundMagneticFieldModelType = Earth::_t96;
    Earth::T96::active_flag        = true;
    Earth::T96::solar_wind_pressure = prm.field.pdyn_nPa  * _NANO_;
    Earth::T96::dst                 = prm.field.dst_nT    * _NANO_;
    Earth::T96::by                  = prm.field.imfBy_nT  * _NANO_;
    Earth::T96::bz                  = prm.field.imfBz_nT  * _NANO_;
    ::T96::SetSolarWindPressure(Earth::T96::solar_wind_pressure);
    ::T96::SetDST(Earth::T96::dst);
    ::T96::SetBYIMF(Earth::T96::by);
    ::T96::SetBZIMF(Earth::T96::bz);
    ::T96::Init(Exosphere::SimulationStartTimeString, Exosphere::SO_FRAME);
  }
  else if (model == "T05") {
    Earth::BackgroundMagneticFieldModelType = Earth::_t05;
    Earth::T05::active_flag        = true;
    Earth::T05::solar_wind_pressure = prm.field.pdyn_nPa  * _NANO_;
    Earth::T05::dst                 = prm.field.dst_nT    * _NANO_;
    Earth::T05::by                  = prm.field.imfBy_nT  * _NANO_;
    Earth::T05::bz                  = prm.field.imfBz_nT  * _NANO_;
    for (int i = 0; i < 6; i++) Earth::T05::W[i] = prm.field.w[i];
    ::T05::SetSolarWindPressure(Earth::T05::solar_wind_pressure);
    ::T05::SetDST(Earth::T05::dst);
    ::T05::SetBXIMF(prm.field.imfBx_nT * _NANO_);
    ::T05::SetBYIMF(Earth::T05::by);
    ::T05::SetBZIMF(Earth::T05::bz);
    ::T05::SetW(Earth::T05::W[0], Earth::T05::W[1], Earth::T05::W[2],
                Earth::T05::W[3], Earth::T05::W[4], Earth::T05::W[5]);
    ::T05::Init(Exosphere::SimulationStartTimeString, Exosphere::SO_FRAME);
  }
  else if (model == "TA16") {
    if (!prm.field.ta16CoeffFile.empty())
      ::TA16::SetCoeffFileName(prm.field.ta16CoeffFile);
    ::TA16::SetSolarWindPressure(prm.field.pdyn_nPa * _NANO_);
    ::TA16::SetSymHc(prm.field.dst_nT * _NANO_);
    ::TA16::SetXIND(prm.field.xind);
    ::TA16::SetBYIMF(prm.field.imfBy_nT * _NANO_);
    ::TA16::Init(Exosphere::SimulationStartTimeString, Exosphere::SO_FRAME);
  }
  // DIPOLE is handled implicitly (no active_flag needed; the field evaluator
  // falls through to the IGRF/analytic dipole path when all model flags are false).
}

// ============================================================================
//  Helper: initialise B/E field values in all AMR cells
// ============================================================================
static void InitMeshFields(const EarthUtil::AmpsParam& prm,
                           cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  const int iMin = -_GHOST_CELLS_X_, iMax = _GHOST_CELLS_X_ + _BLOCK_CELLS_X_ - 1;
  const int jMin = -_GHOST_CELLS_Y_, jMax = _GHOST_CELLS_Y_ + _BLOCK_CELLS_Y_ - 1;
  const int kMin = -_GHOST_CELLS_Z_, kMax = _GHOST_CELLS_Z_ + _BLOCK_CELLS_Z_ - 1;

  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    if (node->block == nullptr) return;

    const int S = (kMax-kMin+1)*(jMax-jMin+1)*(iMax-iMin+1);
    for (int ii = 0; ii < S; ii++) {
      int S1 = ii;
      const int i = iMin + S1 / ((kMax-kMin+1)*(jMax-jMin+1));
      S1 %= (kMax-kMin+1)*(jMax-jMin+1);
      const int j = jMin + S1 / (kMax-kMin+1);
      const int k = kMin + S1 % (kMax-kMin+1);

      const int nd = PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k);
      PIC::Mesh::cDataCenterNode* cn = node->block->GetCenterNode(nd);
      if (cn == nullptr) continue;

      char* offset = cn->GetAssociatedDataBufferPointer()
                   + PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin
                   + PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset;

      double xCell[3];
      xCell[0] = node->xmin[0] + (node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_*(0.5+i);
      xCell[1] = node->xmin[1] + (node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_*(0.5+j);
      xCell[2] = node->xmin[2] + (node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_*(0.5+k);

      double B[3], E[3];
      Earth::Mode3D::EvaluateBackgroundMagneticFieldSI(B, xCell, prm);
      Earth::Mode3D::EvaluateElectricFieldSI(E, xCell, prm);

      for (int idim = 0; idim < 3; idim++) {
        if (PIC::CPLR::DATAFILE::Offset::MagneticField.active)
          *((double*)(offset + PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset + idim*sizeof(double))) = B[idim];
        if (PIC::CPLR::DATAFILE::Offset::ElectricField.active)
          *((double*)(offset + PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset + idim*sizeof(double))) = E[idim];
      }
    }
  }
  else {
    for (int i = 0; i < (1<<DIM); i++)
      if (node->downNode[i] != nullptr)
        InitMeshFields(prm, node->downNode[i]);
  }
}

// ============================================================================
//  Helper: evaluate time step
// ============================================================================
//  dt = DtCellFrac * min_cell_size / v_max(DENS_EMAX)
//
static double EvaluateTimeStep(const EarthUtil::AmpsParam& prm) {
  // Maximum particle speed: relativistic speed at DENS_EMAX
  const double mass   = PIC::MolecularData::GetMass(sSpecies);
  const double E_max  = prm.density3d.Emax_MeV * MeV_in_J;
  const double v_max  = Relativistic::E2Speed(E_max, mass);
  if (!(v_max > 0.0)) return prm.numerics.dtTrace_s;

  // Minimum cell size across the local mesh
  double min_cell = 1.0e30;
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> TraverseMin;
  TraverseMin = [&](cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
      const double cs = node->GetCharacteristicCellSize();
      if (cs < min_cell) min_cell = cs;
    }
    else {
      for (int i = 0; i < (1<<DIM); i++)
        if (node->downNode[i]) TraverseMin(node->downNode[i]);
    }
  };
  TraverseMin(PIC::Mesh::mesh->rootTree);

  // MPI_Allreduce to get global minimum
  double global_min = min_cell;
  MPI_Allreduce(&min_cell, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_GLOBAL_COMMUNICATOR);

  const double dt = DtCellFrac * global_min / v_max;

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForward] Time step: dt=" << dt << " s"
              << "  (v_max=" << v_max/3e8 << " c, min_cell=" << global_min/_EARTH__RADIUS_ << " Re)\n";
  return dt;
}

// ============================================================================
//  BoundaryInjectionSourceRate
// ============================================================================
//
// PURPOSE
// -------
// Wrapper function whose address is assigned to
//   PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate
// before calling
//   PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(spec, startNode)
//
// This follows the same pattern used in srcMOP/main.cpp:
//   PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate = MOP::SourceRate;
//   PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber = N;
//   PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);
//
// The framework then derives the global particle weight internally as:
//   W = BoundaryInjectionSourceRate(spec)
//         x GlobalTimeStep[spec]
//         / maxReferenceInjectedParticleNumber
//     = (pi x integral_{Emin}^{Emax} J(E) dE x A_boundary)
//         x dt / nParticlesPerIter
//
// which is identical to the previous hand-computed formula in EvaluateParticleWeight,
// but now goes through the standard AMPS weight-initialisation machinery so that
// all internal accounting, diagnostic counters, and normalisation paths see a
// consistent value from the very start of the run.
//
// PRECONDITIONS
// -------------
//   - gSpectrum must already be initialised (InitGlobalSpectrumFromKeyValueMap called).
//   - PIC::Mesh::mesh->xGlobalMin/xGlobalMax must be set (amps_init_mesh called).
//   - sSpecies must be set to the forward-mode particle species index.
//
// PHYSICS
// -------
// For an isotropic external differential intensity J(E) [m^-2 s^-1 sr^-1 J^-1],
// the one-way particle flux incident on a surface element is:
//   dN/dA/dt = pi x integral_{Emin}^{Emax} J(E) dE     [particles m^-2 s^-1]
//
// Integrated over the total outer rectangular boundary area A_boundary:
//   R_total = pi x integral J(E) dE x A_boundary        [particles s^-1]
//
static double BoundaryInjectionSourceRate(int spec) {
  // Only the species used by this forward run contributes a non-zero rate.
  if (spec != sSpecies) return 0.0;

  // ---- Log-spaced trapezoidal integration of J(E) over the spectrum range ----
  const int    nIntegPts = 500;
  const double logEmin   = std::log(gSpectrum.Emin_MeV() * MeV_in_J);
  const double logEmax   = std::log(gSpectrum.Emax_MeV() * MeV_in_J);
  const double dlogE     = (logEmax - logEmin) / nIntegPts;

  double integralJ = 0.0;
  for (int i = 0; i < nIntegPts; i++) {
    const double E0 = std::exp(logEmin + i       * dlogE);
    const double E1 = std::exp(logEmin + (i+1.0) * dlogE);
    integralJ += 0.5 * (gSpectrum.GetSpectrum(E0) + gSpectrum.GetSpectrum(E1)) * (E1 - E0);
  }
  // pi sr * integral J dE = one-way isotropic flux  [particles m^-2 s^-1]
  const double onewayFlux = Pi * integralJ;

  // ---- Total outer rectangular boundary area [m^2] ----
  const double Lx = PIC::Mesh::mesh->xGlobalMax[0] - PIC::Mesh::mesh->xGlobalMin[0];
  const double Ly = PIC::Mesh::mesh->xGlobalMax[1] - PIC::Mesh::mesh->xGlobalMin[1];
  const double Lz = PIC::Mesh::mesh->xGlobalMax[2] - PIC::Mesh::mesh->xGlobalMin[2];
  const double totalArea = 2.0 * (Lx*Ly + Ly*Lz + Lz*Lx);

  // Total physical injection rate [particles / s]
  return onewayFlux * totalArea;
}
// ============================================================================
//  Inner absorption sphere: particle-sphere interaction callback
// ============================================================================
//  Called by AMPS when a particle reaches the inner sphere surface.
//  Action: delete the particle (absorbed by Earth / inner boundary).
//
static int ForwardModeParticleSphereInteraction(
    int              /*spec*/,
    long int         ptr,
    double*          /*x*/,
    double*          /*v*/,
    double&          /*dtReturned*/,
    void*            /*sphereData*/,
    void*            /*nodeData*/) {

  // Particles reaching the inner sphere are absorbed by Earth.
  PIC::ParticleBuffer::DeleteParticle(ptr);
  return _PARTICLE_LEFT_THE_DOMAIN_;
}

// ============================================================================
//  Helper: initialise the inner absorbing sphere
// ============================================================================
static void InitAbsorptionSphere(const EarthUtil::AmpsParam& prm) {
  // Earth::Planet is set by amps_init_mesh() when RigidityCalculationMode == _sphere.
  // In Mode3DForward we always use a sphere, so re-configure it here unconditionally.
  sAbsorptionSphere = static_cast<cInternalSphericalData*>(Earth::Planet);
  if (sAbsorptionSphere == nullptr) return;

  // Geometry: sphere at origin, radius from R_INNER (default _EARTH__RADIUS_)
  double sx0[3] = {0.0, 0.0, 0.0};
  const double rInner = prm.domain.rInner * 1000.0;  // km → m
  const double rSphere = (rInner > 0.0) ? rInner : _EARTH__RADIUS_;

  sAbsorptionSphere->SetSphereGeometricalParameters(sx0, rSphere);
  sAbsorptionSphere->Radius = rSphere;

  // Surface mesh
  cInternalSphericalData::SetGeneralSurfaceMeshParameters(nZenithElements, nAzimuthalElements);
  sAbsorptionSphere->localResolution = localSphericalSurfaceResolution;
  sAbsorptionSphere->faceat           = 0;

  // Absorption callback — simply deletes arriving particles
  sAbsorptionSphere->ParticleSphereInteraction = ForwardModeParticleSphereInteraction;

  // No injection from the sphere (forward mode: injection is from outer boundary)
  sAbsorptionSphere->InjectionRate              = nullptr;
  sAbsorptionSphere->InjectionBoundaryCondition = nullptr;

  // Diagnostic surface mesh files
  sAbsorptionSphere->PrintSurfaceMesh("ForwardModeSphere.dat");
  sAbsorptionSphere->PrintSurfaceData("ForwardModeSphereData.dat", 0);

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForward] Absorption sphere: r=" << rSphere/_EARTH__RADIUS_
              << " Re\n";
}

// ============================================================================
//  Helper: inject particles from boundary faces
// ============================================================================
//  For each of the 6 domain-boundary faces, randomly place particles on the face
//  and give them an inward direction sampled from the cosine-weighted hemisphere.
//  Energy is sampled from the boundary spectrum using rejection sampling.
//
//  The distribution object allows future non-isotropic modes (see BoundaryDistribution.h).
//
static void InjectBoundaryParticles(int nParticles,
                                    double weight,
                                    const BoundaryDistributionBase* dist) {
  const double xmin[3] = {PIC::Mesh::mesh->xGlobalMin[0],
                          PIC::Mesh::mesh->xGlobalMin[1],
                          PIC::Mesh::mesh->xGlobalMin[2]};
  const double xmax[3] = {PIC::Mesh::mesh->xGlobalMax[0],
                          PIC::Mesh::mesh->xGlobalMax[1],
                          PIC::Mesh::mesh->xGlobalMax[2]};

  const double Lx = xmax[0] - xmin[0];
  const double Ly = xmax[1] - xmin[1];
  const double Lz = xmax[2] - xmin[2];

  // Face areas and their cumulative sum for random face selection
  const double areas[6] = { Ly*Lz, Ly*Lz,   // ±X faces
                             Lz*Lx, Lz*Lx,   // ±Y faces
                             Lx*Ly, Lx*Ly }; // ±Z faces
  double cumArea[7];
  cumArea[0] = 0.0;
  for (int f = 0; f < 6; f++) cumArea[f+1] = cumArea[f] + areas[f];
  const double totalArea = cumArea[6];

  // Pre-compute species mass and spectrum bounds for energy sampling
  const double mass     = PIC::MolecularData::GetMass(sSpecies);
  const double logEmin  = std::log(gSpectrum.Emin_MeV() * MeV_in_J);
  const double logEmax  = std::log(gSpectrum.Emax_MeV() * MeV_in_J);
  // Upper envelope for rejection sampling: J at Emin (max for typical falling spectra)
  const double Jmax     = gSpectrum.GetSpectrum(gSpectrum.Emin_MeV() * MeV_in_J);

  for (int ip = 0; ip < nParticles; ip++) {
    // ---- Select face proportional to area ----
    const double aRand = rnd() * totalArea;
    int face = 5;
    for (int f = 0; f < 6; f++)
      if (aRand < cumArea[f+1]) { face = f; break; }

    // face 0 = -X, 1 = +X, 2 = -Y, 3 = +Y, 4 = -Z, 5 = +Z
    // Inward normals: (-X face) → (+1,0,0), etc.
    //   face 0: n=(+1,0,0); face 1: n=(-1,0,0)
    //   face 2: n=(0,+1,0); face 3: n=(0,-1,0)
    //   face 4: n=(0,0,+1); face 5: n=(0,0,-1)

    // ---- Random position on chosen face ----
    double xInj[3];
    switch (face) {
      case 0: // -X face
        xInj[0] = xmin[0];
        xInj[1] = xmin[1] + rnd() * Ly;
        xInj[2] = xmin[2] + rnd() * Lz;
        break;
      case 1: // +X face
        xInj[0] = xmax[0];
        xInj[1] = xmin[1] + rnd() * Ly;
        xInj[2] = xmin[2] + rnd() * Lz;
        break;
      case 2: // -Y face
        xInj[0] = xmin[0] + rnd() * Lx;
        xInj[1] = xmin[1];
        xInj[2] = xmin[2] + rnd() * Lz;
        break;
      case 3: // +Y face
        xInj[0] = xmin[0] + rnd() * Lx;
        xInj[1] = xmax[1];
        xInj[2] = xmin[2] + rnd() * Lz;
        break;
      case 4: // -Z face
        xInj[0] = xmin[0] + rnd() * Lx;
        xInj[1] = xmin[1] + rnd() * Ly;
        xInj[2] = xmin[2];
        break;
      default: // +Z face
        xInj[0] = xmin[0] + rnd() * Lx;
        xInj[1] = xmin[1] + rnd() * Ly;
        xInj[2] = xmax[2];
        break;
    }

    // ---- Sample kinetic energy by rejection sampling ----
    // Reject if position is inside the inner sphere (pre-check to save work).
    const double r2 = xInj[0]*xInj[0] + xInj[1]*xInj[1] + xInj[2]*xInj[2];
    const double rInner = sAbsorptionSphere ? sAbsorptionSphere->Radius : _EARTH__RADIUS_;
    if (r2 < rInner * rInner * 0.99 * 0.99) continue; // extremely unlikely but safe

    double E_J = 0.0;
    {
      bool accepted = false;
      for (int attempt = 0; attempt < 200; attempt++) {
        const double Etry = std::exp(logEmin + rnd() * (logEmax - logEmin));
        const double Jtry = gSpectrum.GetSpectrum(Etry);
        if (Jmax > 0.0 && rnd() * Jmax <= Jtry) {
          E_J = Etry;
          accepted = true;
          break;
        }
      }
      if (!accepted) E_J = std::exp(logEmin); // fallback: inject at Emin
    }

    // ---- Sample inward direction (cosine-weighted hemisphere) ----
    // Local frame: e_n is the inward normal, e1 and e2 are tangents.
    double e_n[3] = {0,0,0}, e1[3] = {0,0,0}, e2[3] = {0,0,0};
    switch (face) {
      case 0: e_n[0]= 1; e1[1]= 1; e2[2]= 1; break;
      case 1: e_n[0]=-1; e1[1]= 1; e2[2]= 1; break;
      case 2: e_n[1]= 1; e1[0]= 1; e2[2]= 1; break;
      case 3: e_n[1]=-1; e1[0]= 1; e2[2]= 1; break;
      case 4: e_n[2]= 1; e1[0]= 1; e2[1]= 1; break;
      default: e_n[2]=-1; e1[0]= 1; e2[1]= 1; break;
    }

    // Cosine-weighted hemisphere: cos(θ) = sqrt(ξ₁), φ = 2π ξ₂
    const double xi1 = rnd();
    const double xi2 = rnd();
    const double cosTheta = std::sqrt(xi1);
    const double sinTheta = std::sqrt(1.0 - xi1);
    const double phi = 2.0 * Pi * xi2;
    const double cosPhi = std::cos(phi);
    const double sinPhi = std::sin(phi);

    double v_hat[3];
    for (int idim = 0; idim < 3; idim++)
      v_hat[idim] = cosTheta * e_n[idim]
                  + sinTheta * cosPhi * e1[idim]
                  + sinTheta * sinPhi * e2[idim];

    // ---- Distribution weighting factor ----
    // Retrieve B at injection point for anisotropic modes.
    // For ISOTROPIC this call always returns 1.0.
    double B_inj[3] = {0.0, 0.0, 0.0};
    // (B is not needed for isotropic; skip evaluation for performance)

    const double angularWeight = dist->GetAngularWeight(xInj, v_hat, B_inj, sSpecies);
    if (!(angularWeight > 0.0)) continue;

    // ---- Particle speed from kinetic energy ----
    const double speed  = Relativistic::E2Speed(E_J, mass);
    double v[3];
    for (int idim = 0; idim < 3; idim++) v[idim] = v_hat[idim] * speed;

    // ---- Find cell and inject ----
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode =
        PIC::Mesh::mesh->findTreeNode(xInj);
    if (startNode == nullptr) continue;
    if (startNode->Thread != PIC::ThisThread) continue;

    int iCell, jCell, kCell;
    if (PIC::Mesh::mesh->FindCellIndex(xInj, iCell, jCell, kCell, startNode, false) == -1)
      continue;

    long int newP = PIC::ParticleBuffer::GetNewParticle(
        startNode->block->FirstCellParticleTable[
            iCell + _BLOCK_CELLS_X_ * (jCell + _BLOCK_CELLS_Y_ * kCell)]);

    PIC::ParticleBuffer::byte* pData = PIC::ParticleBuffer::GetParticleDataPointer(newP);

    PIC::ParticleBuffer::SetV(v, pData);
    PIC::ParticleBuffer::SetX(xInj, pData);
    PIC::ParticleBuffer::SetI(sSpecies, pData);

    // Apply angular weighting to individual particle weight correction
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(angularWeight, pData);

    // Per-particle physical weight stored in the block
    startNode->block->SetLocalParticleWeight(weight, sSpecies);

    // Optional: trajectory tracking
    if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
      PIC::ParticleTracker::InitParticleID(pData);
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xInj, v, sSpecies, pData,
                                                              static_cast<void*>(startNode));
    }
  }
}

// ============================================================================
//  Run — main entry point
// ============================================================================
int Run(const EarthUtil::AmpsParam& prm) {
  //--------------------------------------------------------------------------
  // 1. Set model mode
  //--------------------------------------------------------------------------
  Earth::ModelMode = Earth::BoundaryInjectionMode;

  // Domain bounds (convert km → m; same convention as Mode3D)
  // Domain bounds are encoded directly in the PIC AMR mesh
  // (xGlobalMin/xGlobalMax), set by amps_init_mesh() from the .dfn file.
  // The Mode3D::ParsedDomain* variables belong to the cutoff-rigidity path
  // and are not used here.

  //--------------------------------------------------------------------------
  // 2. Register callbacks and configure the field model, THEN init mesh + PIC.
  //--------------------------------------------------------------------------
  // All three assignments below must happen BEFORE amps_init_mesh() / amps_init()
  // because both of those calls reach into main_lib.cpp::amps_init() which, for
  // BoundaryInjectionMode, calls Earth::BoundingBoxInjection::InitDirectionIMF().
  // InitDirectionIMF() needs:
  //   (a) UserDefinedExtraSourceRate -- for the particle weight path (MOP pattern)
  //   (b) s_prm (via SetPrm)         -- so InitDirectionIMF()'s default case can
  //                                    call EvaluateBackgroundMagneticFieldSI
  //   (c) the field model state      -- ConfigureBackgroundFieldModel sets the
  //                                    T96/T05/TA16 active flags and calls their
  //                                    Init(); EvaluateBackgroundMagneticFieldSI
  //                                    dispatches on those flags at runtime.
  //
  // (a) cf. srcMOP/main.cpp:238 -- UserDefinedExtraSourceRate set before amps_init.
  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate = BoundaryInjectionSourceRate;

  // (b) register prm so InitDirectionIMF()'s default branch can call
  //     EvaluateBackgroundMagneticFieldSI(b, x, prm).
  Earth::BoundingBoxInjection::SetPrm(prm);

  // (c) set active flags and initialise T96/T05/TA16 before amps_init() fires
  //     InitDirectionIMF(), which evaluates the background field.
  //     ConfigureBackgroundFieldModel does not touch the mesh, so it is safe
  //     to call here before amps_init_mesh().
  ConfigureBackgroundFieldModel(prm);

  PIC::InitMPI();
  Exosphere::Init_SPICE();
  amps_init_mesh();
  amps_init();

  //--------------------------------------------------------------------------
  // 3. Background field model (mesh-dependent setup, e.g. InitMeshFields)
  //--------------------------------------------------------------------------
  // ConfigureBackgroundFieldModel was already called above (step 2c). The
  // mesh cells are populated here now that amps_init() has allocated blocks.

  //--------------------------------------------------------------------------
  // 4. Initialise B/E fields in all mesh cells
  //--------------------------------------------------------------------------
  InitMeshFields(prm, PIC::Mesh::mesh->rootTree);

  if (prm.mode3dForward.outputInitializedFile) {
    if (PIC::nTotalThreads != 1) {
      // Multi-rank PIC run: outputMeshDataTECPLOT is a collective operation --
      // all PIC ranks must call it together (it gathers cell data via the PIC
      // communicator before rank 0 writes the file).
      PIC::Mesh::mesh->outputMeshDataTECPLOT("amps_3dforward_initialized.data.dat", 0);
    }
    else {
      // Single-rank PIC run: MPI may be configured so each process runs
      // independently (e.g. embarrassingly parallel over parameter space).
      // Only the process that owns MPI_COMM_WORLD rank 0 should write the file.
      int worldRank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
      if (worldRank == 0)
        PIC::Mesh::mesh->outputMeshDataTECPLOT("amps_3dforward_initialized.data.dat", 0);
    }
  }

  //--------------------------------------------------------------------------
  // 5. Initialise global spectrum from #SPECTRUM section
  //    Must happen before step 7 (particle weight) because
  //    BoundaryInjectionSourceRate() integrates gSpectrum over [Emin, Emax].
  //--------------------------------------------------------------------------
  InitGlobalSpectrumFromKeyValueMap(prm.spectrum);

  //--------------------------------------------------------------------------
  // 6. Time step: dt = DtCellFrac * min_cell / v_max(DENS_EMAX)
  //--------------------------------------------------------------------------
  // Identify species index (use species 0 as default; all share the same mass
  // in single-species runs).
  sSpecies = 0;
  sDt = EvaluateTimeStep(prm);
  PIC::ParticleWeightTimeStep::GlobalTimeStep[sSpecies] = sDt;
  cDensity3D::dt_s = sDt;

  //--------------------------------------------------------------------------
  // 7. Particle weight — via initParticleWeight_ConstantWeight (MOP pattern)
  //--------------------------------------------------------------------------
  // BoundaryInjectionSourceRate was registered with UserDefinedExtraSourceRate
  // before amps_init_mesh() (step 2), mirroring the MOP pattern:
  //   srcMOP/main.cpp:238  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate=MOP::SourceRate;
  //   (mesh + amps_init)
  //   srcMOP/main.cpp:400  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=...;
  //   srcMOP/main.cpp:405  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);
  //
  // initParticleWeight_ConstantWeight now derives GlobalParticleWeight as:
  //   W = BoundaryInjectionSourceRate(spec) x GlobalTimeStep[spec]
  //         / maxReferenceInjectedParticleNumber
  //     = (pi x integral J(E) dE x A_boundary) x dt / nParticlesPerIter
  const int nSimPerIter = prm.mode3dForward.nParticlesPerIter;
  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber = nSimPerIter;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(sSpecies);

  // Cache the weight for use in InjectBoundaryParticles and progress logging.
  sParticleWeight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[sSpecies];

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForward] Particle weight: W=" << sParticleWeight
              << " phys/sim  (nParticlesPerIter=" << nSimPerIter
              << ", injectionRate=" << BoundaryInjectionSourceRate(sSpecies)
              << " particles/s)\n";

  //--------------------------------------------------------------------------
  // 8. Inner absorption sphere
  //--------------------------------------------------------------------------
  InitAbsorptionSphere(prm);

  //--------------------------------------------------------------------------
  // 9. 3D density sampling (CG/Moon ExternalSamplingLocalVariables pattern)
  //--------------------------------------------------------------------------
  // Init() registers SampleParticleData and OutputSampledData callbacks.
  cDensity3D::Init(prm);

  //--------------------------------------------------------------------------
  // 10. Select boundary distribution (ISOTROPIC now; extensible for future modes)
  //--------------------------------------------------------------------------
  const BoundaryDistributionType distType =
      ParseBoundaryDistributionType(prm.mode3dForward.boundaryDistType);
  std::unique_ptr<BoundaryDistributionBase> boundaryDist(MakeBoundaryDistribution(distType));

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForward] Boundary distribution: " << boundaryDist->TypeName() << "\n"
              << "[Mode3DForward] Starting " << prm.mode3dForward.nIterations
              << " forward iterations...\n";
  std::cout.flush();

  //--------------------------------------------------------------------------
  // 11. Main loop
  //--------------------------------------------------------------------------
  PIC::Mover::BackwardTimeIntegrationMode = _PIC_MODE_OFF_; // forward in time
  PIC::SamplingMode = _SINGLE_OUTPUT_FILE_SAMPING_MODE_;

  const int nIter = prm.mode3dForward.nIterations;

  for (int iter = 0; iter < nIter; iter++) {
    // ---- Inject particles from boundary ----
    InjectBoundaryParticles(nSimPerIter, sParticleWeight, boundaryDist.get());

    // ---- Advance one time step (moves particles, applies sphere BC,
    //      triggers sampling callbacks registered via ExternalSamplingLocalVariables) ----
    amps_time_step();

    // ---- Progress report every 100 iterations (rank 0 only) ----
    if (PIC::ThisThread == 0 && iter % 100 == 0) {
      const long int nParticles = PIC::ParticleBuffer::GetAllPartNum();
      std::printf("[Mode3DForward] iter=%d/%d  nParticles=%ld\n",
                  iter, nIter, nParticles);
      std::fflush(stdout);
    }
  }

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForward] Forward integration complete.\n";

  //--------------------------------------------------------------------------
  // 12. Force final output (trigger output callbacks at end of run)
  //--------------------------------------------------------------------------
  cDensity3D::OutputSampledData(PIC::DataOutputFileNumber);

  return EXIT_SUCCESS;
}

} // namespace Mode3DForward
} // namespace Earth

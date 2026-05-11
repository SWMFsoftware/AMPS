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
static int    sNParticlesPerIter   = 1000;  // sim particles injected per iteration
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
// ============================================================================
// Boundary injection table (follows BoundaryInjection_SEP.cpp pattern)
// ============================================================================
//
// PURPOSE
// -------
// Pre-computes per-face injection data once from the loaded spectrum gSpectrum
// and the IMF direction b (Earth::BoundingBoxInjection::b set by InitDirectionIMF).
// The design mirrors the SEP model:
//
//   SEP::InitEnergySpectrum()   -> cBoundaryInjectionTable::Init()
//   SEP::EnergyDistributor[]    -> cBoundaryInjectionTable::energyCDF[face][bin]
//   SEP::IntegratedSpectrum[]   -> cBoundaryInjectionTable::fluxPerArea[face]
//   SEP::GetNewParticle()       -> the (mu, phi, b) direction-sampling loop in
//                                  InjectBoundaryParticles()
//
// PHYSICS
// -------
// For an isotropic external differential intensity J(E) [m^-2 s^-1 sr^-1 J^-1]
// the one-way particle flux incident on a surface element is the same for every
// face orientation:
//
//   Phi = pi * integral_{Emin}^{Emax} J(E) dE   [particles m^-2 s^-1]
//
// The per-face flux (and therefore the per-face energy CDF) is identical for
// all six domain faces. Face selection is thus proportional to face area, and
// the energy CDF is built once and shared.
//
// The velocity direction is sampled in the inward hemisphere relative to the IMF
// vector b exactly as in SEP::GetNewParticle():
//
//   mu  = rnd()                       pitch-angle cosine (uniform in [0,1])
//   phi = 2*pi * rnd()               azimuth (uniform)
//   v_hat = mu*b + sqrt(1-mu^2)*(sin(phi)*ee0 + cos(phi)*ee1)
//   accept if v_hat . inwardNormal < 0   (particle enters the domain)
//
// This sampling produces the correct cos(theta) marginal distribution for the
// one-way flux at a planar surface for an isotropic external field.

static constexpr int kNEnergyBins = 500;

struct cBoundaryInjectionTable {
  bool   initialized{false};

  // Log-uniform energy bin edges [J], size kNEnergyBins+1
  double eBinEdge[kNEnergyBins + 1]{};

  // Energy CDF shared across all faces (isotropic: same spectrum for each face).
  // energyCDF[k] = P(E < eBinEdge[k+1]).
  // Populated from the cumulative trapezoid of J(E_mid) * dE.
  // Mirrors SEP::EnergyDistributor (cSingleVariableDiscreteDistribution).
  double energyCDF[kNEnergyBins]{};

  // Per-face one-way flux per unit area [particles m^-2 s^-1].
  // For isotropic injection: fluxPerArea[f] = pi * integralJ for all f.
  // Mirrors SEP::IntegratedSpectrum[6].
  double fluxPerArea[6]{};

  // Face areas [m^2]: face 0=-X,1=+X -> Ly*Lz; 2=-Y,3=+Y -> Lz*Lx; 4=-Z,5=+Z -> Lx*Ly
  double faceArea[6]{};

  // Cumulative face selection weights (fluxPerArea[f]*faceArea[f], normalised).
  // cumFaceWeight[0]=0, cumFaceWeight[6]=1.
  double cumFaceWeight[7]{};
};

static cBoundaryInjectionTable sBndTable;

// ---------------------------------------------------------------------------
// InitBoundaryInjectionTable — analogous to SEP::InitEnergySpectrum()
// ---------------------------------------------------------------------------
static void InitBoundaryInjectionTable() {
  if (sBndTable.initialized) return;

  // ---- 1. Log-uniform energy bins over the spectrum range ----
  const double Emin_J  = gSpectrum.Emin_MeV() * MeV_in_J;
  const double Emax_J  = gSpectrum.Emax_MeV() * MeV_in_J;
  const double logEmin = std::log(Emin_J);
  const double logEmax = std::log(Emax_J);
  const double dlogE   = (logEmax - logEmin) / kNEnergyBins;

  for (int k = 0; k <= kNEnergyBins; k++)
    sBndTable.eBinEdge[k] = std::exp(logEmin + k * dlogE);

  // ---- 2. Raw bin flux values (trapezoid rule) and energy CDF ----
  // rawBin[k] = 0.5*(J(E0)+J(E1))*(E1-E0) — the flux contribution of bin k.
  // Mirrors SEP::ProbabilityTable[] before normalisation.
  double rawBin[kNEnergyBins]{};
  double integralJ = 0.0;

  for (int k = 0; k < kNEnergyBins; k++) {
    const double E0 = sBndTable.eBinEdge[k];
    const double E1 = sBndTable.eBinEdge[k + 1];
    rawBin[k]  = 0.5 * (gSpectrum.GetSpectrum(E0) + gSpectrum.GetSpectrum(E1)) * (E1 - E0);
    integralJ += rawBin[k];
  }

  // Normalised CDF (same as SEP's ProbabilityTable[] after dividing by IntegratedSpectrum)
  double running = 0.0;
  for (int k = 0; k < kNEnergyBins; k++) {
    running += (integralJ > 0.0) ? rawBin[k] / integralJ : 1.0 / kNEnergyBins;
    sBndTable.energyCDF[k] = running;
  }

  // ---- 3. Per-face one-way flux per unit area ----
  // For isotropic: Phi = pi * integralJ  [particles m^-2 s^-1] for every face.
  // This is the 3d_forward analogue of SEP::IntegratedSpectrum[iface].
  const double isoFluxPerArea = Pi * integralJ;
  for (int f = 0; f < 6; f++) sBndTable.fluxPerArea[f] = isoFluxPerArea;

  // ---- 4. Face areas and cumulative selection CDF ----
  const double Lx = PIC::Mesh::mesh->xGlobalMax[0] - PIC::Mesh::mesh->xGlobalMin[0];
  const double Ly = PIC::Mesh::mesh->xGlobalMax[1] - PIC::Mesh::mesh->xGlobalMin[1];
  const double Lz = PIC::Mesh::mesh->xGlobalMax[2] - PIC::Mesh::mesh->xGlobalMin[2];

  sBndTable.faceArea[0] = sBndTable.faceArea[1] = Ly * Lz;
  sBndTable.faceArea[2] = sBndTable.faceArea[3] = Lz * Lx;
  sBndTable.faceArea[4] = sBndTable.faceArea[5] = Lx * Ly;

  double totalWeight = 0.0;
  for (int f = 0; f < 6; f++)
    totalWeight += sBndTable.fluxPerArea[f] * sBndTable.faceArea[f];

  sBndTable.cumFaceWeight[0] = 0.0;
  for (int f = 0; f < 6; f++)
    sBndTable.cumFaceWeight[f + 1] = sBndTable.cumFaceWeight[f]
        + sBndTable.fluxPerArea[f] * sBndTable.faceArea[f] / totalWeight;

  sBndTable.initialized = true;

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForward] BoundaryInjectionTable: integralJ=" << integralJ
              << " m^-2 s^-1 sr^-1, isoFluxPerArea=" << isoFluxPerArea
              << " m^-2 s^-1\n";
}

// ============================================================================
// InjectBoundaryParticles
// ============================================================================
//
// Inject exactly nParticles simulation particles at the domain boundary for
// this iteration, drawing energy and direction from the input-file spectrum.
//
// Implementation follows the SEP pattern from BoundaryInjection_SEP.cpp:
//
//  1. Face selection:   proportional to fluxPerArea[f] * faceArea[f]
//                       (= area for isotropic; analogous to SEP::InjectionRate)
//
//  2. Position:         uniform random on the chosen face
//                       (cf. BoundaryInjection.cpp::InjectionProcessor,
//                            GetBlockFaceCoordinateFrame_3D + random (c0,c1))
//
//  3. Energy:           inverse-CDF sampling from per-face energyCDF[]
//                       (cf. SEP::GetNewParticle -> EnergyDistributor[nface])
//
//  4. Direction:        inward half-hemisphere relative to IMF vector b
//                       (cf. SEP::GetNewParticle, default InjectionMode)
//                         mu  = rnd()          pitch-angle cosine (uniform)
//                         phi = 2*pi*rnd()     azimuth (uniform)
//                         v = mu*b + sqrt(1-mu^2)*(sin(phi)*ee0 + cos(phi)*ee1)
//                         accept if v . inwardNormal < 0
//
//  5. Inject:           GetNewParticle / SetV / SetX / SetI
//                       (same as BoundaryInjection.cpp::InjectionProcessor)
//
// ============================================================================
// InjectBoundaryParticles
// ============================================================================
//
// Inject exactly nParticles simulation particles at the domain boundary for
// this iteration, drawing energy and direction from the input-file spectrum.
//
// Follows the SEP injection pattern (BoundaryInjection_SEP.cpp):
//
//  1. Face selection:  proportional to fluxPerArea[f] * faceArea[f]
//                      (mirrors SEP::InjectionRate returning IntegratedSpectrum[nface])
//
//  2. Position:        uniform random on the chosen face
//                      (cf. BoundaryInjection.cpp::InjectionProcessor:
//                           x = x0 + c0*e0 + c1*e1 from GetBlockFaceCoordinateFrame_3D)
//
//  3. Energy:          inverse-CDF sampling from energyCDF[]
//                      (cf. SEP::GetNewParticle -> EnergyDistributor[nface].DistributeVariable()
//                           then e = e0 + rnd()*(e1-e0) within the selected bin)
//
//  4. Direction:       inward half-hemisphere relative to IMF vector b
//                      (cf. SEP::GetNewParticle, default InjectionMode:
//                           mu=rnd(), phi=2*pi*rnd(),
//                           v = mu*b + sqrt(1-mu^2)*(sin(phi)*ee0 + cos(phi)*ee1),
//                           while (v.ExternalNormal >= 0) repeat)
//
//  5. Inject:          SetV / SetX / SetI + particle tracker hook
//                      (same as BoundaryInjection.cpp::InjectionProcessor)
//
// ---------------------------------------------------------------------------
// InjectParticles — UserDefinedParticleInjectionFunction callback
// ---------------------------------------------------------------------------
// Signature matches PIC::BC::fUserDefinedParticleInjectionFunction:
//   long int (*)()
// Assigned to PIC::BC::UserDefinedParticleInjectionFunction in Run() before
// amps_init_mesh(), following the pattern in srcMOP/main.cpp line 239:
//   PIC::BC::UserDefinedParticleInjectionFunction = MOP::InjectParticles;
// PIC::TimeStep() (via amps_time_step()) then calls it automatically each
// iteration — no explicit call is needed in the main loop.
static long int InjectParticles() {
  const int    nParticles = sNParticlesPerIter;
  const double weight     = sParticleWeight;

  // Pre-compute injection tables from gSpectrum on first call (lazy, idempotent).
  // Analogous to SEP::InjectionRate() calling InitEnergySpectrum() on first use.
  InitBoundaryInjectionTable();

  const double* xmin  = PIC::Mesh::mesh->xGlobalMin;
  const double* xmax  = PIC::Mesh::mesh->xGlobalMax;
  const double  Lx    = xmax[0] - xmin[0];
  const double  Ly    = xmax[1] - xmin[1];
  const double  Lz    = xmax[2] - xmin[2];
  const double  mass  = PIC::MolecularData::GetMass(sSpecies);

  // Per-face inward normals and fixed tangent-vector pairs.
  // Each face has an orthonormal frame (n_in, t1, t2) used for direction sampling.
  // n_in points into the domain; t1 and t2 span the face plane.
  static const double inwardNormal[6][3] = {
    { 1, 0, 0}, {-1, 0, 0},   // face 0=-X, face 1=+X
    { 0, 1, 0}, { 0,-1, 0},   // face 2=-Y, face 3=+Y
    { 0, 0, 1}, { 0, 0,-1}    // face 4=-Z, face 5=+Z
  };
  // Tangent vectors completing the right-handed frame for each face.
  static const double tangent1[6][3] = {
    {0, 1, 0}, {0, 1, 0},     // face 0,1: t1 = +Y
    {0, 0, 1}, {0, 0, 1},     // face 2,3: t1 = +Z
    {1, 0, 0}, {1, 0, 0}      // face 4,5: t1 = +X
  };
  static const double tangent2[6][3] = {
    {0, 0, 1}, {0, 0, 1},     // face 0,1: t2 = +Z
    {1, 0, 0}, {1, 0, 0},     // face 2,3: t2 = +X
    {0, 1, 0}, {0, 1, 0}      // face 4,5: t2 = +Y
  };

  for (int ip = 0; ip < nParticles; ip++) {

    // ------------------------------------------------------------------
    // 1. Select face from pre-computed CDF (proportional to area for isotropic)
    //    SEP analogue: InjectionProcessor loops over nface with ExternalFaces[nface];
    //    here the face is drawn probabilistically rather than processed per-block.
    // ------------------------------------------------------------------
    const double fRand = rnd();
    int face = 5;
    for (int f = 0; f < 6; f++)
      if (fRand < sBndTable.cumFaceWeight[f + 1]) { face = f; break; }

    // ------------------------------------------------------------------
    // 2. Uniform random position on the chosen face, offset slightly inward.
    //
    //    The injection point is placed ONE CELL WIDTH inward from the outer
    //    boundary face along the inward normal.  This is necessary because
    //    PIC::Mesh::mesh->findTreeNode() treats the domain as a half-open
    //    interval: a coordinate exactly equal to xmin[d] or xmax[d] lies
    //    outside the valid range and findTreeNode returns NULL, which causes
    //    the particle to be silently dropped.
    //
    //    The reference BoundaryInjection.cpp avoids this entirely by iterating
    //    over pre-built boundary-block lists (InitBoundingBoxInjectionBlockList)
    //    and never calling findTreeNode on an outer-boundary coordinate.  Here
    //    we replicate the same guarantee with a small offset equal to 1e-6 of
    //    the relevant domain dimension -- small enough that the injection point
    //    is still effectively on the boundary face for sampling purposes.
    // ------------------------------------------------------------------
    double xInj[3];
    const double offX = (xmax[0]-xmin[0]) * 1.0e-6;
    const double offY = (xmax[1]-xmin[1]) * 1.0e-6;
    const double offZ = (xmax[2]-xmin[2]) * 1.0e-6;
    switch (face) {
      case 0: xInj[0]=xmin[0]+offX; xInj[1]=xmin[1]+rnd()*Ly; xInj[2]=xmin[2]+rnd()*Lz; break;
      case 1: xInj[0]=xmax[0]-offX; xInj[1]=xmin[1]+rnd()*Ly; xInj[2]=xmin[2]+rnd()*Lz; break;
      case 2: xInj[0]=xmin[0]+rnd()*Lx; xInj[1]=xmin[1]+offY; xInj[2]=xmin[2]+rnd()*Lz; break;
      case 3: xInj[0]=xmin[0]+rnd()*Lx; xInj[1]=xmax[1]-offY; xInj[2]=xmin[2]+rnd()*Lz; break;
      case 4: xInj[0]=xmin[0]+rnd()*Lx; xInj[1]=xmin[1]+rnd()*Ly; xInj[2]=xmin[2]+offZ; break;
      default: xInj[0]=xmin[0]+rnd()*Lx; xInj[1]=xmin[1]+rnd()*Ly; xInj[2]=xmax[2]-offZ; break;
    }

    // Skip positions inside the inner absorption sphere (Earth's surface)
    {
      const double r2     = xInj[0]*xInj[0] + xInj[1]*xInj[1] + xInj[2]*xInj[2];
      const double rInner = sAbsorptionSphere ? sAbsorptionSphere->Radius : _EARTH__RADIUS_;
      if (r2 < rInner * rInner) continue;
    }

    // ------------------------------------------------------------------
    // 3. Sample kinetic energy using inverse-CDF on sBndTable.energyCDF[]
    //    SEP: EnergyDistributor[nface].DistributeVariable() -> iInterval
    //         e = e0 + rnd()*(e1-e0) within the selected bin (eV, then *eV2J)
    //    Here bins are already in Joules.
    // ------------------------------------------------------------------
    double E_J;
    {
      const double u = rnd();
      int k = kNEnergyBins - 1;
      for (int kk = 0; kk < kNEnergyBins; kk++)
        if (u <= sBndTable.energyCDF[kk]) { k = kk; break; }
      E_J = sBndTable.eBinEdge[k] + rnd() * (sBndTable.eBinEdge[k+1] - sBndTable.eBinEdge[k]);
    }

    // ------------------------------------------------------------------
    // 4. Sample inward velocity direction: isotropic 4pi, cosine-weighted.
    //
    // For an isotropic external intensity f(Omega) = 1/(4*pi) the one-way
    // flux through a surface element is weighted by v_norm = v * cos(theta).
    // The conditional PDF of an inward direction, given that the particle
    // crosses the surface, is:
    //
    //   p(theta, phi) = (1/pi) * cos(theta)  for theta in [0, pi/2], phi in [0, 2*pi)
    //
    // This is sampled exactly (no rejection loop) by cosine-weighted hemisphere
    // sampling:
    //
    //   cos(theta) = sqrt(xi_1)      <- accounts for v_norm * f(v) factor
    //   phi        = 2*pi * xi_2
    //   v_hat = cos(theta)*n_in + sin(theta)*(cos(phi)*t1 + sin(phi)*t2)
    //
    // where (n_in, t1, t2) is the face-local orthonormal frame defined above.
    // ------------------------------------------------------------------
    double v_hat[3];
    {
      const double cosTheta = std::sqrt(rnd());   // samples v_norm*f(v) exactly
      const double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
      const double phi      = 2.0 * Pi * rnd();
      const double cosPhi   = std::cos(phi);
      const double sinPhi   = std::sin(phi);
      for (int d = 0; d < 3; d++)
        v_hat[d] = cosTheta * inwardNormal[face][d]
                 + sinTheta * (cosPhi * tangent1[face][d] + sinPhi * tangent2[face][d]);
    }

    // ------------------------------------------------------------------
    // 5. Physical speed from relativistic kinetic energy
    //    (same as SEP::GetNewParticle: Speed = Relativistic::E2Speed(e, mass))
    // ------------------------------------------------------------------
    const double speed = Relativistic::E2Speed(E_J, mass);
    double v[3];
    for (int d = 0; d < 3; d++) v[d] = v_hat[d] * speed;

    // ------------------------------------------------------------------
    // 6. Locate AMR cell and inject particle
    //    (mirrors BoundaryInjection.cpp::InjectionProcessor:
    //     GetNewParticle + CloneParticle + SetX/SetI + MoveParticle)
    // ------------------------------------------------------------------
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode =
        PIC::Mesh::mesh->findTreeNode(xInj);
    if (!startNode || startNode->Thread != PIC::ThisThread) continue;

    int iCell, jCell, kCell;
    if (PIC::Mesh::mesh->FindCellIndex(xInj, iCell, jCell, kCell, startNode, false) == -1)
      continue;

    long int newP = PIC::ParticleBuffer::GetNewParticle(); 
    PIC::ParticleBuffer::byte* pData = PIC::ParticleBuffer::GetParticleDataPointer(newP);

    PIC::ParticleBuffer::SetV(v,      pData);
    PIC::ParticleBuffer::SetX(xInj,   pData);
    PIC::ParticleBuffer::SetI(sSpecies, pData);

    if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_)
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0, pData);

    startNode->block->SetLocalParticleWeight(weight, sSpecies);

#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    PIC::ParticleTracker::InitParticleID(pData);
    PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xInj, v, sSpecies, pData,
                                                            static_cast<void*>(startNode));
#endif

    // Move the freshly injected particle for a random fraction of the local
    // time step, placing it at a statistically uniform position along its
    // trajectory within the first iteration interval.
    // This is identical to BoundaryInjection.cpp line 219:
    //   _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle, rnd()*LocalTimeStep, startNode)
    // Without this step all injected particles would start on the domain face,
    // producing a spurious density spike at the boundary in the first sample.
    const double localDt = startNode->block->GetLocalTimeStep(sSpecies);
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newP, rnd() * localDt, startNode);
  }  // end for (int ip = 0; ip < nParticles; ip++)
  return nParticles;
}  // end InjectParticles

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

  // Register injection callback — mirrors srcMOP/main.cpp:239:
  //   PIC::BC::UserDefinedParticleInjectionFunction = MOP::InjectParticles;
  // PIC::TimeStep() calls this automatically each iteration; the direct
  // InjectBoundaryParticles() call in the main loop is therefore removed.
  PIC::BC::UserDefinedParticleInjectionFunction = InjectParticles;

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
  // Store in module-level static so InjectParticles() (called by
  // PIC::TimeStep via UserDefinedParticleInjectionFunction) can read it.
  sNParticlesPerIter = prm.mode3dForward.nParticlesPerIter;
  const int nSimPerIter = sNParticlesPerIter;
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
  // 10. Pre-compute boundary injection table from spectrum + IMF direction.
  //--------------------------------------------------------------------------
  // InitBoundaryInjectionTable() is also called lazily on the first injection,
  // but calling it here ensures any startup diagnostics are printed before the
  // main loop begins. Analogous to SEP::Init() -> InitEnergySpectrum().
  InitBoundaryInjectionTable();

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForward] Boundary injection: SEP-style spectrum CDF sampling\n"
              << "[Mode3DForward] IMF direction: b=(" << Earth::BoundingBoxInjection::b[0]
              << ", " << Earth::BoundingBoxInjection::b[1]
              << ", " << Earth::BoundingBoxInjection::b[2] << ")\n"
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
    // ---- Advance one time step ----
    // InjectParticles() is called automatically by PIC::TimeStep() via
    // PIC::BC::UserDefinedParticleInjectionFunction (set in step 2, above).
    // No explicit injection call is needed here.
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

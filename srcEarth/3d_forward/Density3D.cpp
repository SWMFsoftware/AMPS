//======================================================================================
// Density3D.cpp
//======================================================================================
//
// Implementation of per-cell 3D volumetric density sampling for Mode3DForward.
// See Density3D.h for the full design description and indexing conventions.
//
// PARALLEL STRATEGY
// -----------------
// Each MPI rank accumulates densityBuffer_ independently for its local blocks.
// At output time OutputSampledData() performs an MPI_Reduce(SUM) across all ranks
// so that rank 0 receives the complete domain-wide sum. The normalised result is
// then broadcast and rank 0 writes the Tecplot file.
// In replicated-domain mode (PIC::nTotalThreads == 1 per rank, Mode3D style) every
// rank owns every block, so MPI_Reduce still produces the correct aggregate.
//
//======================================================================================

#include "Density3D.h"
#include "../util/amps_param_parser.h"


#include "pic.h"
//#include "../constants.h"
#include "../Earth.h"

#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>

namespace Earth {
namespace Mode3DForward {

//======================================================================================
// Static member definitions
//======================================================================================
int    cDensity3D::nEnergyBins      = 30;
double cDensity3D::Emin_J           = 1.0 * 1.602176634e-13; // 1 MeV in J
double cDensity3D::Emax_J           = 2.0e7 * 1.602176634e-13; // 20000 MeV in J
bool   cDensity3D::logSpacing        = true;

int    cDensity3D::nLocalBlocks_    = 0;
int    cDensity3D::nCellsPerBlock_  = 0;
int    cDensity3D::nTotalCells_     = 0;
long int cDensity3D::nSampleSteps_  = 0;

double cDensity3D::dt_s             = 1.0;

std::vector<double> cDensity3D::densityBuffer_;
std::vector<double> cDensity3D::densitySampled_;

//======================================================================================
// Energy-bin helpers
//======================================================================================

double cDensity3D::GetBinLowJ(int iE) {
  if (logSpacing) {
    const double logMin = std::log(Emin_J);
    const double logMax = std::log(Emax_J);
    const double dlog   = (logMax - logMin) / nEnergyBins;
    return std::exp(logMin + iE * dlog);
  }
  else {
    const double dE = (Emax_J - Emin_J) / nEnergyBins;
    return Emin_J + iE * dE;
  }
}

double cDensity3D::GetBinHighJ(int iE) {
  return GetBinLowJ(iE + 1);
}

double cDensity3D::GetBinCentreJ(int iE) {
  return 0.5 * (GetBinLowJ(iE) + GetBinHighJ(iE));
}

double cDensity3D::GetBinWidthJ(int iE) {
  return GetBinHighJ(iE) - GetBinLowJ(iE);
}

//======================================================================================
// Init
//======================================================================================

// Internal helper: map kinetic energy to energy bin index.
// Returns -1 if outside [Emin, Emax].
static int EnergyToBin(double E_J,
                       int nBins, double Emin, double Emax, bool logSp) {
  if (E_J < Emin || E_J > Emax) return -1;
  int iE;
  if (logSp) {
    const double logMin = std::log(Emin);
    const double logMax = std::log(Emax);
    const double dlog   = (logMax - logMin) / nBins;
    iE = static_cast<int>((std::log(E_J) - logMin) / dlog);
  }
  else {
    const double dE = (Emax - Emin) / nBins;
    iE = static_cast<int>((E_J - Emin) / dE);
  }
  if (iE < 0)      iE = 0;
  if (iE >= nBins) iE = nBins - 1;
  return iE;
}

void cDensity3D::Init(const EarthUtil::AmpsParam& prm) {
  // ---- Energy grid from #DENSITY_3D section ----
  constexpr double MeV_in_J = 1.602176634e-13;
  nEnergyBins = prm.density3d.nEnergyBins;
  Emin_J      = prm.density3d.Emin_MeV * MeV_in_J;
  Emax_J      = prm.density3d.Emax_MeV * MeV_in_J;
  logSpacing  = (prm.density3d.spacing == EarthUtil::Density3DParam::Spacing::LOG);

  // ---- Cell buffer dimensions ----
  nCellsPerBlock_ = _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_;
  nLocalBlocks_   = PIC::DomainBlockDecomposition::nLocalBlocks;
  nTotalCells_    = nLocalBlocks_ * nCellsPerBlock_;

  // Allocate and zero-initialise
  densityBuffer_.assign(static_cast<std::size_t>(nTotalCells_) * nEnergyBins, 0.0);
  densitySampled_.assign(static_cast<std::size_t>(nTotalCells_) * nEnergyBins, 0.0);
  nSampleSteps_ = 0;

  // ---- Register sampling callbacks (Moon/CG pattern) ----
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(
      SampleParticleData,
      OutputSampledData);

  if (PIC::ThisThread == 0) {
    std::cout << "[Density3D] Initialized:"
              << " nBlocks=" << nLocalBlocks_
              << " nCellsPerBlock=" << nCellsPerBlock_
              << " nEnergyBins=" << nEnergyBins
              << " Emin=" << prm.density3d.Emin_MeV << " MeV"
              << " Emax=" << prm.density3d.Emax_MeV << " MeV"
              << " spacing=" << (logSpacing ? "LOG" : "LINEAR")
              << "\n";
  }
}

//======================================================================================
// SampleParticleData — registered sampling function (CG/Moon pattern)
//======================================================================================
// Called once per PIC iteration. Traverses all local blocks and accumulates
// the differential number-density contribution from every particle.
//
// Contribution of one particle with weight W in a cell of volume V:
//   Δn(iE) += W / (V * ΔE_iE)    [m^-3 J^-1]
//
// After N_sample iterations the time-averaged density is:
//   n_avg(iE) = sum(Δn) / N_sample
//
void cDensity3D::SampleParticleData() {
  // Guard: if the block table has grown (mesh refinement), the buffer is stale.
  // In practice Mode3DForward runs on a frozen mesh, so this should never trigger.
  if (PIC::DomainBlockDecomposition::nLocalBlocks != nLocalBlocks_) {
    if (PIC::ThisThread == 0)
      std::cerr << "[Density3D] WARNING: block count changed after Init; "
                   "density buffer may be inconsistent.\n";
    return;
  }

  const int nBins  = nEnergyBins;
  const double Emin = Emin_J;
  const double Emax = Emax_J;
  const bool   logSp = logSpacing;

  for (int iBlock = 0; iBlock < nLocalBlocks_; iBlock++) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[iBlock];
    if (node == nullptr || node->block == nullptr) continue;

    // Cell volume: use cell-size³ (valid for cubic cells in the AMR tree).
    const double cellSize = node->GetCharacteristicCellSize();
    const double cellVol  = cellSize * cellSize * cellSize;
    if (!(cellVol > 0.0)) continue;

    for (int i = 0; i < _BLOCK_CELLS_X_; i++)
    for (int j = 0; j < _BLOCK_CELLS_Y_; j++)
    for (int k = 0; k < _BLOCK_CELLS_Z_; k++) {
      const int localIdx = i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k);
      const int base     = (iBlock * nCellsPerBlock_ + localIdx) * nBins;

      long int ptr = node->block->FirstCellParticleTable[localIdx];
      while (ptr != -1) {
        PIC::ParticleBuffer::byte* pData =
            PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        const int spec = PIC::ParticleBuffer::GetI(pData);

        // Particle weight (physical particles per simulation particle)
        double weight = node->block->GetLocalParticleWeight(spec);
        weight *= PIC::ParticleBuffer::GetIndividualStatWeightCorrection(pData);

        // Kinetic energy from velocity
        double v[3];
        PIC::ParticleBuffer::GetV(v, pData);
        const double vSq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        const double mass = PIC::MolecularData::GetMass(spec);
        const double E_J  = Relativistic::Speed2E(std::sqrt(vSq), mass);

        // Find energy bin and accumulate
        const int iE = EnergyToBin(E_J, nBins, Emin, Emax, logSp);
        if (iE >= 0) {
          // n contribution [m^-3 J^-1]
          const double dE = cDensity3D::GetBinWidthJ(iE);
          if (dE > 0.0)
            densityBuffer_[base + iE] += weight / (cellVol * dE);
        }

        ptr = PIC::ParticleBuffer::GetNext(pData);
      }
    }
  }

  nSampleSteps_++;
}

//======================================================================================
// OutputSampledData — registered output function (CG/Moon pattern)
//======================================================================================
// 1. MPI_Reduce sum across all ranks → rank 0 has global totals.
// 2. Normalise by nSampleSteps_ to get time-averaged density.
// 3. Rank 0 writes Tecplot ASCII file.
// 4. Reset buffer for the next output window.
//
void cDensity3D::OutputSampledData(int dataOutputFileNumber) {
  const int bufSize = nTotalCells_ * nEnergyBins;
  if (bufSize <= 0) return;

  // ---- MPI gather ----
  std::vector<double> globalBuf(bufSize, 0.0);

  MPI_Reduce(densityBuffer_.data(), globalBuf.data(),
             bufSize, MPI_DOUBLE, MPI_SUM, 0,
             MPI_GLOBAL_COMMUNICATOR);

  // ---- Normalise and cache ----
  if (PIC::ThisThread == 0 && nSampleSteps_ > 0) {
    const double norm = 1.0 / static_cast<double>(nSampleSteps_);
    for (int i = 0; i < bufSize; i++)
      densitySampled_[i] = globalBuf[i] * norm;
  }

  // ---- Tecplot output (rank 0 only) ----
  if (PIC::ThisThread == 0) {
    char fname[512];
    std::snprintf(fname, sizeof(fname),
                  "%s/density3d.out=%04d.dat",
                  PIC::OutputDataFileDirectory, dataOutputFileNumber);

    FILE* fout = std::fopen(fname, "w");
    if (fout == nullptr) {
      std::cerr << "[Density3D] Cannot open output file: " << fname << "\n";
      return;
    }

    // ---- Header ----
    std::fprintf(fout,
        "TITLE=\"AMPS 3D Forward: Volumetric Particle Density\"\n"
        "VARIABLES=\"x_Re\" \"y_Re\" \"z_Re\" \"CellSize_Re\"");

    constexpr double MeV_in_J = 1.602176634e-13;
    const double Re = _EARTH__RADIUS_;
    for (int iE = 0; iE < nEnergyBins; iE++) {
      const double E_lo_MeV = GetBinLowJ(iE)    / MeV_in_J;
      const double E_hi_MeV = GetBinHighJ(iE)   / MeV_in_J;
      std::fprintf(fout, " \"n[%g-%g_MeV]_m3J\"", E_lo_MeV, E_hi_MeV);
    }
    std::fprintf(fout, " \"n_total_m3\"\n");
    std::fprintf(fout,
        "ZONE T=\"AMPS_3DFORWARD_DENSITY\", F=POINT\n");

    // ---- Data rows ----
    for (int iBlock = 0; iBlock < nLocalBlocks_; iBlock++) {
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
          PIC::DomainBlockDecomposition::BlockTable[iBlock];
      if (node == nullptr || node->block == nullptr) continue;

      const double sz = node->GetCharacteristicCellSize();

      for (int i = 0; i < _BLOCK_CELLS_X_; i++)
      for (int j = 0; j < _BLOCK_CELLS_Y_; j++)
      for (int k = 0; k < _BLOCK_CELLS_Z_; k++) {
        // Cell-centre position [m]
        const double x = node->xmin[0] + (node->xmax[0]-node->xmin[0])
                         / _BLOCK_CELLS_X_ * (0.5 + i);
        const double y = node->xmin[1] + (node->xmax[1]-node->xmin[1])
                         / _BLOCK_CELLS_Y_ * (0.5 + j);
        const double z = node->xmin[2] + (node->xmax[2]-node->xmin[2])
                         / _BLOCK_CELLS_Z_ * (0.5 + k);

        const int localIdx = i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k);
        const int base = (iBlock * nCellsPerBlock_ + localIdx) * nEnergyBins;

        std::fprintf(fout, "%.6e %.6e %.6e %.6e",
                     x/Re, y/Re, z/Re, sz/Re);

        double n_total = 0.0;
        for (int iE = 0; iE < nEnergyBins; iE++) {
          const double n_perJ = densitySampled_[base + iE];
          std::fprintf(fout, " %.6e", n_perJ);
          n_total += n_perJ * GetBinWidthJ(iE);
        }
        std::fprintf(fout, " %.6e\n", n_total);
      }
    }

    std::fclose(fout);
    std::cout << "[Density3D] Output written: " << fname << "\n";
  }

  // ---- Reset accumulator for next window ----
  std::fill(densityBuffer_.begin(), densityBuffer_.end(), 0.0);
  nSampleSteps_ = 0;
}

} // namespace Mode3DForward
} // namespace Earth

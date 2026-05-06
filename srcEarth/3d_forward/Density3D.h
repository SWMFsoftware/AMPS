#ifndef _SRC_EARTH_3D_FORWARD_DENSITY3D_H_
#define _SRC_EARTH_3D_FORWARD_DENSITY3D_H_

//======================================================================================
// Density3D.h
//======================================================================================
//
// PURPOSE
// -------
// Per-cell volumetric particle number-density sampler for Mode3DForward.
// Accumulates differential number density n(E) [m^-3 J^-1] in each AMR cell
// across energy bins defined by the #DENSITY_3D input section.
//
// DESIGN FOLLOWS CG/MOON PATTERN
// --------------------------------
// Sampling is registered with PIC::Sampling::ExternalSamplingLocalVariables
// (same pattern as Comet::Sampling::InjectedDustSizeDistribution and
//  Moon::Sampling::VelocityDistribution).
//
// Two callbacks are registered:
//   SampleParticleData()    — called every PIC iteration to accumulate density
//   OutputSampledData(int)  — called each output cycle to gather MPI contributions
//                             and write Tecplot
//
// INDEXING
// --------
// After amps_init(), the mesh block table is stable. A sequential block index
// iBlock ∈ [0, nLocalBlocks) maps to PIC::DomainBlockDecomposition::BlockTable[iBlock].
// Each block contains nCellsPerBlock = _BLOCK_CELLS_X_ × _BLOCK_CELLS_Y_ × _BLOCK_CELLS_Z_ cells.
//
// The flat density buffer is indexed as:
//
//   densityBuffer[(iBlock * nCellsPerBlock + localCellIdx) * nEnergyBins + iE]
//
// where localCellIdx = i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k).
//
// UNITS
// -----
// densitySampled[cell][iE] — time-averaged differential number density
//   [particles m^-3 J^-1]
//
// To convert to [particles m^-3 MeV^-1]:
//   n_perMeV = n_perJ * MeV_in_J   (where MeV_in_J = 1.602176634e-13)
//
// To convert to total number density in bin iE:
//   N_iE [m^-3] = n_perJ * GetBinWidthJ(iE)
//
//======================================================================================

#include <vector>
#include <string>

// Forward-declare PIC types so this header remains independent of pic.h
template <class T> class cTreeNodeAMR;
class cDataBlockAMR;

namespace EarthUtil { struct AmpsParam; }

namespace Earth {
namespace Mode3DForward {

class cDensity3D {
public:
  // -------------------------------------------------------------------
  // Initialisation
  // -------------------------------------------------------------------
  // Must be called once after amps_init_mesh() and amps_init() have
  // completed so that the block table is stable and cell measures are set.
  //
  // Allocates the density buffer, sets up energy bins, and registers
  // SampleParticleData / OutputSampledData with
  //   PIC::Sampling::ExternalSamplingLocalVariables.
  //
  static void Init(const EarthUtil::AmpsParam& prm);

  // -------------------------------------------------------------------
  // Registered callbacks
  // -------------------------------------------------------------------
  // SampleParticleData — traverse every local block×cell, sum weighted
  //   particle contributions into densityBuffer_.
  static void SampleParticleData();

  // OutputSampledData — MPI_Reduce across ranks, normalise by sample
  //   count × dt, write Tecplot file "density3d_<N>.dat".
  static void OutputSampledData(int dataOutputFileNumber);

  // -------------------------------------------------------------------
  // Energy-bin accessors
  // -------------------------------------------------------------------
  static int    nEnergyBins;          ///< number of energy bins (from DENS_NENERGY)
  static double Emin_J;               ///< lower bound [J]  (= DENS_EMIN * MeV_in_J)
  static double Emax_J;               ///< upper bound [J]  (= DENS_EMAX * MeV_in_J)
  static bool   logSpacing;           ///< true for LOG, false for LINEAR

  // Return energy bounds and width of bin iE (0-based) in Joules.
  static double GetBinLowJ(int iE);
  static double GetBinHighJ(int iE);
  static double GetBinCentreJ(int iE);
  static double GetBinWidthJ(int iE);

  // -------------------------------------------------------------------
  // Buffer sizes (set by Init)
  // -------------------------------------------------------------------
  static int nLocalBlocks_;     ///< PIC::DomainBlockDecomposition::nLocalBlocks snapshot
  static int nCellsPerBlock_;   ///< _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_
  static int nTotalCells_;      ///< nLocalBlocks_ * nCellsPerBlock_
  static long int nSampleSteps_;///< iterations accumulated since last reset

  // Flat sampling accumulator; reset to zero at start of each output window.
  // Size: nTotalCells_ * nEnergyBins
  static std::vector<double> densityBuffer_;

  // Time-averaged result written to Tecplot; populated by OutputSampledData.
  static std::vector<double> densitySampled_;

  // Physical time step (set externally by Mode3DForward::Run before Init).
  static double dt_s;

private:
  // Not instantiable; all members are static.
  cDensity3D() = delete;
};

} // namespace Mode3DForward
} // namespace Earth

#endif // _SRC_EARTH_3D_FORWARD_DENSITY3D_H_

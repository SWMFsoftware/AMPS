/*
================================================================================
        WAVE-NUMBER-RESOLVED ALFVÉN TURBULENCE MODEL FOR SEP TRANSPORT
================================================================================

This file implements a spectral extension of the existing Kolmogorov turbulence
model.  The legacy model carries only two branch-integrated segment energies,
E+ and E-.  That is efficient and remains the default.  The model implemented
here carries a log-uniform wave-number distribution for each branch:

    E_+(k_j,s,t),  E_-(k_j,s,t),  j=0,...,NK-1 .

The spectral model is deliberately integrated with the existing AMPS/SEP code in
a minimally invasive way:

  * The existing particle movers continue to accumulate the k-binned streaming
    arrays G_plus_streaming[j] and G_minus_streaming[j].

  * The new coupling manager applies the j-th growth/damping rate only to the
    j-th wave-energy bin, so particles affect the wave-number bin selected by
    their resonant k.

  * Existing output and scattering diagnostics continue to use the compact
    branch-integrated datum CellIntegratedWaveEnergy and the derived datum
    WaveEnergyDensity.  Therefore, after every spectral operation, this file
    refreshes CellIntegratedWaveEnergy by summing the spectral bins.

  * Advection, particle coupling, reflection, and nonlinear cascade are all
    implemented directly on E_±(k_j).  The compact integrated E_± datum is kept
    only as a diagnostic/compatibility view obtained by summing the spectral
    bins after each spectral operator.

  * The shock-turbulence source is still formulated as an integrated energy
    increment.  ProjectIntegratedEnergyToSpectrum() is therefore retained for
    that one operator: it maps the shock-updated integrated energy back to the
    spectrum while preserving the previous spectral shape whenever possible.

The code below is intentionally verbose.  The detailed comments explain the
units and data layout at every modified location because these arrays are easy
to misinterpret: the stored values are cell-integrated energies per logarithmic
wave-number bin, not pointwise spectral densities.
================================================================================
*/

#include "turbulence_wave_number_resolved.h"
#include "../sep.h"
#include "../wave_particle_coupling_kolmogorov.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <vector>

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {
namespace WaveNumberResolved {

using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::NK;
using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::K_MIN;
using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::K_MAX;
using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::DLNK;
using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::Q;
using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::M;
using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::GetKBinIndex;
using SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::RedistributeWaveEnergyToParticles;

// Default remains the legacy branch-integrated model.  The CLI switches this to
// WaveNumberResolved only when explicitly requested by the user.
ModelMode TurbulenceModelMode = ModelMode::Integrated;

// The datum is hidden from the generic Tecplot field-line output.  It stores
// 2*NK values per field-line segment: first all outward-wave bins, then all
// inward-wave bins.  The string name is empty because doPrint=false.
PIC::Datum::cDatumStored SpectralWaveEnergy(2*NK,"",false);

// Hidden diagnostic datum for wave-energy exchange rates.  It is intentionally
// not printed by the generic AMPS field-line writer because the number of
// columns would be large: two branches times several source terms times NK
// wave-number bins.  OutputSpectrumTecplot2D() writes these rates in a dedicated
// 2-D Tecplot file together with W_+(s,k), W_-(s,k), and sigma_c(k).
PIC::Datum::cDatumStored SpectralWaveEnergyExchangeRate(
    2*nSpectralEnergyExchangeRateComponents*NK,"",false);

namespace {

const double kTinyEnergy = 1.0e-300;
const double kMu0 = 4.0e-7 * M_PI;
const double kProtonMass = 1.67262192e-27;

// -----------------------------------------------------------------------------
// Cached right-boundary W-(k) density.
// -----------------------------------------------------------------------------
// For the integrated model the right boundary stores one pre-existing W- value.
// For the spectral model the boundary must preserve a complete W-(k_j) density.
// We store density rather than energy because the plotted and physically useful
// quantity is W=E/V.  If the boundary-cell volume is recomputed, the energy that
// corresponds to the fixed density must be E_j = W_j * V_current.
// -----------------------------------------------------------------------------
std::vector<std::vector<double> > RightBoundaryInitialWminusSpectrum;
std::vector<char> RightBoundaryInitialSpectrumIsSet;
int RightBoundaryInitialSpectrumNFieldLine = -1;

void EnsureRightBoundarySpectrumStorage() {
  if (RightBoundaryInitialSpectrumNFieldLine != PIC::FieldLine::nFieldLine) {
    RightBoundaryInitialWminusSpectrum.assign(
        PIC::FieldLine::nFieldLine, std::vector<double>(NK,0.0));
    RightBoundaryInitialSpectrumIsSet.assign(PIC::FieldLine::nFieldLine,0);
    RightBoundaryInitialSpectrumNFieldLine = PIC::FieldLine::nFieldLine;
  }
}

// Return the center wave number for a logarithmic bin.
inline double KCenter(int j) {
  return K_MIN * std::exp(j * DLNK);
}

// Kolmogorov initialization weights for a log-uniform grid.
//
// If P(k) = C k^{-5/3} is the one-dimensional wave-power spectrum and the grid
// is uniform in ln k, then the energy in one log bin is approximately
//
//   E_j ≈ P(k_j) Δk_j = C k_j^{-5/3} (k_j Δln k)
//       ∝ k_j^{-2/3}.
//
// The common Δln k cancels when the weights are normalized, so we use
// k_j^{-2/3}.  This allows the sum over bins to exactly match the existing
// branch-integrated energy stored in CellIntegratedWaveEnergy.
void BuildKolmogorovLogBinWeights(std::vector<double>& weights) {
  weights.assign(NK,0.0);

  double sum = 0.0;
  for (int j=0; j<NK; ++j) {
    const double k = KCenter(j);
    weights[j] = std::pow(k,-2.0/3.0);
    sum += weights[j];
  }

  if (sum <= 0.0 || !std::isfinite(sum)) {
    const double uniform = 1.0/static_cast<double>(NK);
    std::fill(weights.begin(),weights.end(),uniform);
    return;
  }

  for (int j=0; j<NK; ++j) weights[j] /= sum;
}

// Geometry helper copied from the integrated advection model.  It returns the
// Alfvén speed and flux-tube cross-sectional area at a field-line vertex.  Every
// spectral bin uses the same wave propagation speed because the current model is
// nondispersive Alfvénic transport along the mean field.
bool GetFaceAlfvenSpeedAndArea(
    PIC::FieldLine::cFieldLineVertex* vertex,
    int field_line_idx,
    double& V_A,
    double& area) {
  namespace FL = PIC::FieldLine;

  if (!vertex) return false;

  double* B0 = vertex->GetDatum_ptr(FL::DatumAtVertexMagneticField);
  double n_sw = 0.0;
  vertex->GetDatum(FL::DatumAtVertexPlasmaDensity, &n_sw);

  if (!B0 || n_sw <= 0.0) return false;

  const double B_mag = std::sqrt(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2]);
  if (B_mag <= 0.0) return false;

  const double rho = n_sw * kProtonMass;
  V_A = B_mag / std::sqrt(kMu0 * rho);

  double* x = vertex->GetX();
  const double radius = SEP::FieldLine::MagneticTubeRadius(x, field_line_idx);
  area = M_PI * radius * radius;

  return std::isfinite(V_A) && std::isfinite(area) && V_A > 0.0 && area > 0.0;
}

inline double ComputeAdvectiveFlux(
    double source_energy,
    double source_volume,
    PIC::FieldLine::cFieldLineVertex* face_vertex,
    int field_line_idx,
    double dt,
    bool limit_to_source_energy) {
  if (source_energy <= 0.0 || source_volume <= 0.0 || dt <= 0.0) return 0.0;

  double V_A = 0.0, area = 0.0;
  if (!GetFaceAlfvenSpeedAndArea(face_vertex, field_line_idx, V_A, area)) return 0.0;

  const double source_density = source_energy / source_volume;
  double flux = V_A * source_density * area * dt;
  if (!std::isfinite(flux) || flux <= 0.0) return 0.0;

  if (limit_to_source_energy) flux = std::min(flux,source_energy);
  return flux;
}

inline double ComputeAdvectiveFluxFromDensity(
    double source_density,
    PIC::FieldLine::cFieldLineVertex* face_vertex,
    int field_line_idx,
    double dt) {
  if (source_density <= 0.0 || dt <= 0.0) return 0.0;

  double V_A = 0.0, area = 0.0;
  if (!GetFaceAlfvenSpeedAndArea(face_vertex, field_line_idx, V_A, area)) return 0.0;

  const double flux = V_A * source_density * area * dt;
  return (std::isfinite(flux) && flux > 0.0) ? flux : 0.0;
}

// Return a small, always-positive denominator for divisions that are allowed to
// approach zero only because a physical quantity is absent.  The value is small
// enough to be irrelevant for normal heliospheric parameters, but it prevents a
// zero density, zero length, or zero wave energy from producing NaN values that
// can contaminate all spectral bins after MPI synchronization.
inline double PositiveFloor(double x, double floor_value) {
  return (std::isfinite(x) && x > floor_value) ? x : floor_value;
}

// Extract the proton mass density at the center of a field-line segment.  The
// spectral cascade operator uses the same Elsasser-energy normalization as the
// integrated cascade operator, W^±=(rho/4) Z_±^2.  Therefore a local mass density
// is needed to convert wave-energy density into the counter-propagating Elsasser
// amplitude that controls nonlinear decorrelation.
bool GetSegmentMassDensity(
    PIC::FieldLine::cFieldLineSegment* segment,
    double& rho) {
  if (!segment) return false;

  double n_sw = 0.0;
  segment->GetPlasmaDensity(0.5,n_sw);

  // In the rest of the SEP turbulence code GetPlasmaDensity() is interpreted as
  // a proton number density.  Keep the same convention here so the spectral and
  // integrated cascade operators use the same physical normalization.
  rho = n_sw*_H__MASS_;

  return std::isfinite(rho) && rho > 0.0;
}

// Branch helper: the first NK entries are E+(k), the second NK entries are E-(k).
// The helpers are defined before all spectral operators because advection,
// reflection, cascade, diagnostics, and particle coupling all use the same
// memory layout.  Keeping the indexing in one place reduces the chance of
// accidentally mixing E+ and E- bins when new spectral operators are added.
inline int OffsetPlus(int j) { return j; }
inline int OffsetMinus(int j) { return NK+j; }

// Index helper for SpectralWaveEnergyExchangeRate.  All exchange-rate terms use
// the same branch/k layout as SpectralWaveEnergy, but include an additional
// component index that identifies the physical process.  The helper is kept in
// one place so future diagnostics do not accidentally mix branch, component,
// and wave-number offsets.
inline int RateOffset(int branch, int component, int j) {
  return ((branch*nSpectralEnergyExchangeRateComponents + component)*NK + j);
}

inline int BranchPlus() { return 0; }
inline int BranchMinus() { return 1; }

// Add a signed integrated exchange rate to the hidden diagnostic datum.  The
// operators work with cell-integrated bin energies E_\pm(k_j) [J per log-k bin],
// so the natural rate accumulated here is dE/dt [J/s per log-k bin].  The
// Tecplot writer later divides by the segment volume to obtain the plotted
// density rate dW/dt [J m^{-3} s^{-1} per log-k bin].
inline void AddIntegratedEnergyExchangeRate(
    double* rates,
    int branch,
    SpectralEnergyExchangeRateComponent component,
    int j,
    double dE,
    double dt_normalization) {

  if (!rates) return;
  if (branch != BranchPlus() && branch != BranchMinus()) return;
  if (component < 0 || component >= nSpectralEnergyExchangeRateComponents) return;
  if (j < 0 || j >= NK) return;
  if (!(dt_normalization > 0.0) || !std::isfinite(dt_normalization)) return;
  if (!std::isfinite(dE)) return;

  rates[RateOffset(branch,component,j)] += dE/dt_normalization;
}

// Compute a safe normalized rate (1/W) dW/dt.  This quantity is useful because
// it is the local e-folding rate of the wave-energy density.  The floor avoids
// extremely large diagnostic values in bins where the wave energy is exactly
// zero; those bins carry no physical turbulence and a normalized rate would be
// ill-defined.
inline double NormalizedRate(double density_rate, double wave_density) {
  const double floor = 1.0e-300;
  if (!(wave_density > floor) || !std::isfinite(wave_density)) return 0.0;
  if (!std::isfinite(density_rate)) return 0.0;
  return density_rate/wave_density;
}

// Compute the large-scale, non-WKB reflection rate for a single segment.  This
// is the same physical rate used by the legacy integrated reflection operator,
// but this helper returns it to the wave-number-resolved code so the rate can be
// applied to every E_+(k_j),E_-(k_j) bin independently.
//
// The estimate is based on the gradient of the Alfvén speed along the magnetic
// field line:
//
//   G_R = C_R/2 * | V_A * d ln(V_A) / ds | .
//
// Reflection does not move energy in k-space; it changes the propagation branch
// at the same wave number.  Therefore the caller uses this scalar G_R to mix
// E_+(k_j) and E_-(k_j) inside each bin j.
bool GetSegmentReflectionRate(
    PIC::FieldLine::cFieldLineSegment* segment,
    double C_reflection,
    double grad_floor,
    double& G_reflection) {
  namespace FL = PIC::FieldLine;

  G_reflection = 0.0;
  if (!segment || C_reflection <= 0.0) return false;

  FL::cFieldLineVertex* v0 = segment->GetBegin();
  FL::cFieldLineVertex* v1 = segment->GetEnd();
  if (!v0 || !v1) return false;

  const double ds = segment->GetLength();
  if (!(ds > 0.0) || !std::isfinite(ds)) return false;

  double* B0 = v0->GetDatum_ptr(FL::DatumAtVertexMagneticField);
  double* B1 = v1->GetDatum_ptr(FL::DatumAtVertexMagneticField);
  if (!B0 || !B1) return false;

  double n0 = 0.0, n1 = 0.0;
  v0->GetDatum(FL::DatumAtVertexPlasmaDensity,&n0);
  v1->GetDatum(FL::DatumAtVertexPlasmaDensity,&n1);
  if (!(n0 > 0.0) || !(n1 > 0.0)) return false;

  const double B2_0 = B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2];
  const double B2_1 = B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2];
  if (!(B2_0 > 0.0) || !(B2_1 > 0.0)) return false;

  const double rho0 = n0*kProtonMass;
  const double rho1 = n1*kProtonMass;
  if (!(rho0 > 0.0) || !(rho1 > 0.0)) return false;

  const double ratio = (B2_1*rho0)/(B2_0*rho1);
  if (!(ratio > 0.0) || !std::isfinite(ratio)) return false;

  const double dlnVA_ds = 0.5*std::log(ratio)/ds;
  if (grad_floor > 0.0 && std::fabs(dlnVA_ds) < grad_floor) return false;

  const double VA2_0 = B2_0/(kMu0*rho0);
  const double VA2_1 = B2_1/(kMu0*rho1);
  if (!(VA2_0 > 0.0) || !(VA2_1 > 0.0)) return false;

  const double VA_c = std::sqrt(0.5*(VA2_0+VA2_1));
  if (!(VA_c > 0.0) || !std::isfinite(VA_c)) return false;

  G_reflection = 0.5*C_reflection*std::fabs(VA_c*dlnVA_ds);
  return std::isfinite(G_reflection) && G_reflection > 0.0;
}

// Apply one explicit, positivity-preserving cascade sweep to one local spectrum.
//
// The spectral cascade is represented as an energy flux in logarithmic wave
// number.  For each branch and each bin j, a fraction of E_±(k_j) is transferred
// to the next higher-k bin.  The highest bin has no resolved neighbor, so the
// outgoing flux from that bin is interpreted as unresolved dissipation/heating.
// This is a simple finite-volume model in ln(k), not a diagnostic rescaling of
// the integrated spectrum.
//
// The nonlinear rate uses the same counter-propagating-wave idea as the legacy
// integrated operator,
//
//   a_±(k_j) = C_nl f_sigma sqrt(W_∓(k_j)) / [sqrt(rho) lambda_perp]
//              * (k_j/k_min)^(2/3) .
//
// The factor (k/k_min)^(2/3) shortens the nonlinear time toward smaller scales,
// as expected for a Kolmogorov-like cascade.  The exponential fraction
// 1-exp(-a dt) keeps the update positive even when dt is larger than the local
// cascade time.
void CascadeOneSpectralSweep(
    double* spectrum,
    double* exchange_rates,
    double volume,
    double rho,
    double dt,
    double rate_normalization_dt,
    double C_nl,
    double lambda_perp_m,
    bool enable_cross_helicity_modulation,
    double& transferred_energy,
    double& dissipated_energy) {

  if (!spectrum) return;
  if (!(volume > 0.0) || !(rho > 0.0) || !(dt > 0.0)) return;
  if (!(C_nl > 0.0) || !(lambda_perp_m > 0.0)) return;

  // For diagnostics, an operator called with multiple sub-sweeps should report
  // the average exchange rate over the full physical operator time step, not
  // the instantaneous rate within one sub-sweep.  rate_normalization_dt is
  // therefore the total CascadeSpectrumAllFieldLines() dt, while dt is the
  // sub-sweep time used by this explicit update.
  const double rate_dt = (rate_normalization_dt > 0.0 && std::isfinite(rate_normalization_dt))
                         ? rate_normalization_dt : dt;

  const double sqrt_rho = std::sqrt(PositiveFloor(rho,1.0e-60));
  const double lambda = PositiveFloor(lambda_perp_m,1.0);

  std::vector<double> transfer_plus(NK,0.0);
  std::vector<double> transfer_minus(NK,0.0);

  for (int j=0; j<NK; ++j) {
    const double Eplus = std::max(0.0,spectrum[OffsetPlus(j)]);
    const double Eminus = std::max(0.0,spectrum[OffsetMinus(j)]);

    const double Wplus = Eplus/volume;
    const double Wminus = Eminus/volume;

    // Per-bin normalized cross helicity.  When one branch is absent, nonlinear
    // interaction should be strongly suppressed because Alfvénic cascade needs
    // counter-propagating wave packets.  This mirrors the integrated cascade
    // option, but applies the suppression separately at each k_j.
    double f_sigma = 1.0;
    if (enable_cross_helicity_modulation) {
      const double Wsum = Wplus + Wminus;
      if (Wsum > 0.0) {
        const double sigma_c = (Wplus-Wminus)/Wsum;
        f_sigma = std::sqrt(std::max(0.0,1.0-sigma_c*sigma_c));
      }
      else {
        f_sigma = 0.0;
      }
    }

    const double k_ratio = PositiveFloor(KCenter(j)/K_MIN,1.0);
    const double spectral_rate_factor = std::pow(k_ratio,2.0/3.0);

    const double a_plus = C_nl*f_sigma*std::sqrt(std::max(0.0,Wminus))/(sqrt_rho*lambda)*spectral_rate_factor;
    const double a_minus = C_nl*f_sigma*std::sqrt(std::max(0.0,Wplus))/(sqrt_rho*lambda)*spectral_rate_factor;

    const double frac_plus = (a_plus > 0.0 && std::isfinite(a_plus)) ? (1.0-std::exp(-a_plus*dt)) : 0.0;
    const double frac_minus = (a_minus > 0.0 && std::isfinite(a_minus)) ? (1.0-std::exp(-a_minus*dt)) : 0.0;

    transfer_plus[j] = Eplus*std::min(1.0,std::max(0.0,frac_plus));
    transfer_minus[j] = Eminus*std::min(1.0,std::max(0.0,frac_minus));
  }

  for (int j=0; j<NK; ++j) {
    spectrum[OffsetPlus(j)] = std::max(0.0,spectrum[OffsetPlus(j)]-transfer_plus[j]);
    spectrum[OffsetMinus(j)] = std::max(0.0,spectrum[OffsetMinus(j)]-transfer_minus[j]);

    if (j+1 < NK) {
      spectrum[OffsetPlus(j+1)] += transfer_plus[j];
      spectrum[OffsetMinus(j+1)] += transfer_minus[j];

      // Cascade is an exchange in wave-number space.  At a given k-bin it is a
      // signed source/sink: the donor bin loses energy and the receiving bin
      // gains the same amount.  Recording both signs lets the Tecplot output show
      // where the cascade is depleting power and where it is depositing power in
      // the resolved spectrum.
      AddIntegratedEnergyExchangeRate(exchange_rates,BranchPlus(),RateCascade,j,
          -transfer_plus[j],rate_dt);
      AddIntegratedEnergyExchangeRate(exchange_rates,BranchMinus(),RateCascade,j,
          -transfer_minus[j],rate_dt);
      AddIntegratedEnergyExchangeRate(exchange_rates,BranchPlus(),RateCascade,j+1,
          +transfer_plus[j],rate_dt);
      AddIntegratedEnergyExchangeRate(exchange_rates,BranchMinus(),RateCascade,j+1,
          +transfer_minus[j],rate_dt);

      transferred_energy += transfer_plus[j] + transfer_minus[j];
    }
    else {
      // Energy that leaves the resolved k-range at k_max is not placed into an
      // unresolved spectral bin.  It is counted as turbulent dissipation/heating.
      // In the diagnostic rate, this is a true sink from the highest resolved bin.
      AddIntegratedEnergyExchangeRate(exchange_rates,BranchPlus(),RateCascade,j,
          -transfer_plus[j],rate_dt);
      AddIntegratedEnergyExchangeRate(exchange_rates,BranchMinus(),RateCascade,j,
          -transfer_minus[j],rate_dt);

      dissipated_energy += transfer_plus[j] + transfer_minus[j];
    }
  }
}

// Sum one branch of the spectral datum.
inline double SumBranch(const double* spectrum, bool plus_branch) {
  double sum = 0.0;
  const int offset = plus_branch ? 0 : NK;
  for (int j=0; j<NK; ++j) sum += std::max(0.0,spectrum[offset+j]);
  return sum;
}

// Redistribute a wave-energy change for a single k-bin to particles that are
// resonant with that same k-bin and wave branch.
//
// This is the part that makes the new model genuinely wave-number resolved for
// particle coupling.  The legacy code sums the wave-energy change over all k and
// redistributes it to all particles resonant with the branch.  Here the wave
// energy change of bin j is redistributed only to particles whose local
// resonant wave number falls into that same bin.
void RedistributeWaveNumberBinEnergyToParticles(
    PIC::FieldLine::cFieldLineSegment* segment,
    double particle_energy_change,
    int sigma,
    int kbin) {

  if (!segment || segment->Thread != PIC::ThisThread) return;
  if (sigma != +1 && sigma != -1) return;
  if (kbin < 0 || kbin >= NK) return;
  if (std::fabs(particle_energy_change) < 1.0e-50) return;

  double rho_tmp = 0.0;
  segment->GetPlasmaDensity(0.5,rho_tmp);
  const double rho = rho_tmp * _H__MASS_;
  if (rho <= 0.0) return;

  double B[3];
  segment->GetMagneticField(0.5,B);
  const double B0 = Vector3D::Length(B);
  if (B0 <= 0.0) return;

  const double vA = B0/std::sqrt(VacuumPermeability*rho);
  const double sigma_vA = sigma*vA;
  const double c2 = SpeedOfLight*SpeedOfLight;

  std::vector<long int> particle_idx;
  std::vector<double> resonance_weight;
  std::vector<double> current_speed;
  std::vector<double> current_vParallel;
  std::vector<double> current_vNormal;
  std::vector<double> stat_weight;
  std::vector<double> particle_mass;

  double Wsum = 0.0;

  long int p = segment->FirstParticleIndex;
  while (p != -1) {
    const double vParallel = PIC::ParticleBuffer::GetVParallel(p);
    const double vNormal = PIC::ParticleBuffer::GetVNormal(p);
    const double v2 = vParallel*vParallel + vNormal*vNormal;

    if (v2 <= 0.0 || v2 >= c2) {
      p = PIC::ParticleBuffer::GetNext(p);
      continue;
    }

    const double v_mag = std::sqrt(v2);
    const double mu = vParallel/v_mag;

    // Keep the same branch-selection convention as the legacy redistribution:
    // sigma=+1 (outward waves) exchanges energy with inward-moving particles,
    // sigma=-1 (inward waves) exchanges energy with outward-moving particles.
    if ((sigma > 0 && mu > 0.0) || (sigma < 0 && mu < 0.0) || std::fabs(mu) < 1.0e-12) {
      p = PIC::ParticleBuffer::GetNext(p);
      continue;
    }

    const double Omega = std::fabs(Q)*B0/M;
    const double kRes = Omega/(std::fabs(mu)*v_mag);
    const int particle_kbin = GetKBinIndex(kRes);

    if (particle_kbin != kbin) {
      p = PIC::ParticleBuffer::GetNext(p);
      continue;
    }

    const int species = PIC::ParticleBuffer::GetI(p);
    if (species < 0 || species >= PIC::nTotalSpecies) {
      p = PIC::ParticleBuffer::GetNext(p);
      continue;
    }

    const double mass = PIC::MolecularData::GetMass(species);
    const double gamma = 1.0/std::sqrt(1.0-v2/c2);
    const double p_momentum = gamma*mass*v_mag;
    const double w_cnt = PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] *
                         PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);

    if (w_cnt <= 0.0) {
      p = PIC::ParticleBuffer::GetNext(p);
      continue;
    }

    const double K_sigma = v_mag*mu - sigma_vA;
    const double wi_s = w_cnt*p_momentum*p_momentum*mu*K_sigma;

    if (std::fabs(wi_s) > 0.0 && std::isfinite(wi_s)) {
      particle_idx.push_back(p);
      resonance_weight.push_back(wi_s);
      current_speed.push_back(v_mag);
      current_vParallel.push_back(vParallel);
      current_vNormal.push_back(vNormal);
      stat_weight.push_back(w_cnt);
      particle_mass.push_back(mass);
      Wsum += wi_s;
    }

    p = PIC::ParticleBuffer::GetNext(p);
  }

  // If no particles are resonant with this exact bin, fall back to the legacy
  // branch redistribution.  This prevents wave-energy changes from being lost
  // in very sparse Monte-Carlo cells while still using the bin-resolved path
  // whenever resonant particles are present.
  if (particle_idx.empty() || std::fabs(Wsum) < 1.0e-60) {
    RedistributeWaveEnergyToParticles(segment,particle_energy_change,sigma);
    return;
  }

  const double scale = particle_energy_change/Wsum;
  for (std::size_t i=0; i<particle_idx.size(); ++i) {
    const double dE_i = scale*resonance_weight[i];
    const double dE_physical = dE_i/stat_weight[i];

    const double v_current = current_speed[i];
    const double gamma_current = 1.0/std::sqrt(1.0-v_current*v_current/c2);
    const double Ek_current = (gamma_current-1.0)*particle_mass[i]*c2;
    double Ek_new = Ek_current + dE_physical;
    if (Ek_new < 0.0) Ek_new = 0.0;

    if (Ek_new > 0.0) {
      const double gamma_new = Ek_new/(particle_mass[i]*c2) + 1.0;
      const double v_new = SpeedOfLight*std::sqrt(1.0 - 1.0/(gamma_new*gamma_new));
      const double scale_v = v_new/v_current;
      PIC::ParticleBuffer::SetVParallel(current_vParallel[i]*scale_v,particle_idx[i]);
      PIC::ParticleBuffer::SetVNormal(current_vNormal[i]*scale_v,particle_idx[i]);
    }
    else {
      PIC::ParticleBuffer::SetVParallel(0.0,particle_idx[i]);
      PIC::ParticleBuffer::SetVNormal(0.0,particle_idx[i]);
    }
  }
}

} // anonymous namespace

void ResetRightBoundarySpectrumInitialCondition() {
  RightBoundaryInitialWminusSpectrum.clear();
  RightBoundaryInitialSpectrumIsSet.clear();
  RightBoundaryInitialSpectrumNFieldLine = -1;
}

void ResetSpectralEnergyExchangeRates() {
  // The exchange-rate datum is a per-iteration diagnostic.  It must be cleared
  // before the wave-particle coupling, reflection, and cascade operators are
  // applied; otherwise the Tecplot file would show an accumulated history rather
  // than the rate associated with the current main-loop update.
  //
  // Only the owner rank writes locally-owned segments.  The MPI gather at the
  // end makes the zeroed values visible in the replicated FieldLinesAll arrays
  // used by rank 0 for output.
  if (!IsActive()) return;

  const int n = SpectralWaveEnergyExchangeRate.length;

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment || segment->Thread != PIC::ThisThread) continue;

      double* rates = segment->GetDatum_ptr(SpectralWaveEnergyExchangeRate);
      if (!rates) continue;

      for (int q=0; q<n; ++q) rates[q] = 0.0;
    }
  }

  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergyExchangeRate);
}

void InitializeSpectrumFromIntegratedEnergy(PIC::Datum::cDatumStored& IntegratedWaveEnergy) {
  std::vector<double> weights;
  BuildKolmogorovLogBinWeights(weights);

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment || segment->Thread != PIC::ThisThread) continue;

      double* integrated = segment->GetDatum_ptr(IntegratedWaveEnergy);
      double* spectrum = segment->GetDatum_ptr(SpectralWaveEnergy);
      if (!integrated || !spectrum) continue;

      const double Eplus = std::max(0.0,integrated[0]);
      const double Eminus = std::max(0.0,integrated[1]);

      for (int j=0; j<NK; ++j) {
        spectrum[OffsetPlus(j)] = Eplus*weights[j];
        spectrum[OffsetMinus(j)] = Eminus*weights[j];
      }
    }
  }

  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergy);
  UpdateIntegratedEnergyFromSpectrum(IntegratedWaveEnergy);

  if (PIC::ThisThread == 0) {
    std::cout << "Initialized wave-number-resolved turbulence spectrum from integrated E+/E- "
              << "using NK=" << NK << ", k_min=" << K_MIN << " m^-1, k_max=" << K_MAX
              << " m^-1 and Kolmogorov log-bin weights.\n";
  }
}

void UpdateIntegratedEnergyFromSpectrum(PIC::Datum::cDatumStored& IntegratedWaveEnergy) {
  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment || segment->Thread != PIC::ThisThread) continue;

      double* integrated = segment->GetDatum_ptr(IntegratedWaveEnergy);
      double* spectrum = segment->GetDatum_ptr(SpectralWaveEnergy);
      if (!integrated || !spectrum) continue;

      integrated[0] = SumBranch(spectrum,true);
      integrated[1] = SumBranch(spectrum,false);

      if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
        validate_numeric(integrated[0],__LINE__,__FILE__);
        validate_numeric(integrated[1],__LINE__,__FILE__);
      }
    }
  }
}

void ProjectIntegratedEnergyToSpectrum(PIC::Datum::cDatumStored& IntegratedWaveEnergy) {
  std::vector<double> weights;
  BuildKolmogorovLogBinWeights(weights);

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment || segment->Thread != PIC::ThisThread) continue;

      double* integrated = segment->GetDatum_ptr(IntegratedWaveEnergy);
      double* spectrum = segment->GetDatum_ptr(SpectralWaveEnergy);
      if (!integrated || !spectrum) continue;

      const double target_plus = std::max(0.0,integrated[0]);
      const double target_minus = std::max(0.0,integrated[1]);
      const double current_plus = SumBranch(spectrum,true);
      const double current_minus = SumBranch(spectrum,false);

      if (current_plus > kTinyEnergy) {
        const double scale = target_plus/current_plus;
        for (int j=0; j<NK; ++j) spectrum[OffsetPlus(j)] = std::max(0.0,spectrum[OffsetPlus(j)]*scale);
      }
      else {
        for (int j=0; j<NK; ++j) spectrum[OffsetPlus(j)] = target_plus*weights[j];
      }

      if (current_minus > kTinyEnergy) {
        const double scale = target_minus/current_minus;
        for (int j=0; j<NK; ++j) spectrum[OffsetMinus(j)] = std::max(0.0,spectrum[OffsetMinus(j)]*scale);
      }
      else {
        for (int j=0; j<NK; ++j) spectrum[OffsetMinus(j)] = target_minus*weights[j];
      }
    }
  }

  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergy);
}

void CaptureRightBoundarySpectrumInitialCondition() {
  EnsureRightBoundarySpectrumStorage();

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    if (nseg <= 0) continue;

    PIC::FieldLine::cFieldLineSegment* segment_last = field_line->GetSegment(nseg-1);
    if (!segment_last || segment_last->Thread != PIC::ThisThread) continue;

    double* spectrum = segment_last->GetDatum_ptr(SpectralWaveEnergy);
    if (!spectrum) continue;

    const double volume_last = SEP::FieldLine::GetSegmentVolume(segment_last,field_line_idx);
    if (volume_last <= 0.0) continue;

    for (int j=0; j<NK; ++j) {
      RightBoundaryInitialWminusSpectrum[field_line_idx][j] =
          std::max(0.0,spectrum[OffsetMinus(j)])/volume_last;
    }
    RightBoundaryInitialSpectrumIsSet[field_line_idx] = 1;
  }
}

void EnforceRightBoundarySpectrumInitialCondition() {
  EnsureRightBoundarySpectrumStorage();

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    if (!RightBoundaryInitialSpectrumIsSet[field_line_idx]) continue;

    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    if (nseg <= 0) continue;

    PIC::FieldLine::cFieldLineSegment* segment_last = field_line->GetSegment(nseg-1);
    if (!segment_last || segment_last->Thread != PIC::ThisThread) continue;

    double* spectrum = segment_last->GetDatum_ptr(SpectralWaveEnergy);
    if (!spectrum) continue;

    const double volume_last = SEP::FieldLine::GetSegmentVolume(segment_last,field_line_idx);
    if (volume_last <= 0.0) continue;

    for (int j=0; j<NK; ++j) {
      spectrum[OffsetMinus(j)] = RightBoundaryInitialWminusSpectrum[field_line_idx][j]*volume_last;
    }
  }
}

void AdvectSpectrumAllFieldLines(double dt, double TurbulenceLevelBeginning, double TurbulenceLevelEnd) {
  (void)TurbulenceLevelEnd;

  if (dt <= 0.0) return;

  std::vector<double> kolmogorov_weights;
  BuildKolmogorovLogBinWeights(kolmogorov_weights);

  int processed_field_lines = 0;

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    if (nseg < 1) continue;

    bool has_local_segments = false;
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (segment && segment->Thread == PIC::ThisThread) {
        has_local_segments = true;
        break;
      }
    }
    if (!has_local_segments) continue;

    std::vector<double> delta_plus(static_cast<std::size_t>(nseg)*NK,0.0);
    std::vector<double> delta_minus(static_cast<std::size_t>(nseg)*NK,0.0);

    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment_i = field_line->GetSegment(i);
      if (!segment_i || segment_i->Thread != PIC::ThisThread) continue;

      double* spectrum_i = segment_i->GetDatum_ptr(SpectralWaveEnergy);
      if (!spectrum_i) continue;

      const double volume_i = SEP::FieldLine::GetSegmentVolume(segment_i,field_line_idx);
      if (volume_i <= 0.0) continue;

      for (int j=0; j<NK; ++j) {
        const std::size_t idx = static_cast<std::size_t>(i)*NK + j;

        const double Eplus_i = std::max(0.0,spectrum_i[OffsetPlus(j)]);
        const double Eminus_i = std::max(0.0,spectrum_i[OffsetMinus(j)]);

        // E+(k_j): outward propagation, increasing segment index.
        delta_plus[idx] -= ComputeAdvectiveFlux(
            Eplus_i, volume_i, segment_i->GetEnd(), field_line_idx, dt,
            /*limit_to_source_energy=*/true);

        if (i > 0) {
          PIC::FieldLine::cFieldLineSegment* segment_left = field_line->GetSegment(i-1);
          if (segment_left) {
            double* spectrum_left = segment_left->GetDatum_ptr(SpectralWaveEnergy);
            const double volume_left = SEP::FieldLine::GetSegmentVolume(segment_left,field_line_idx);
            if (spectrum_left && volume_left > 0.0) {
              delta_plus[idx] += ComputeAdvectiveFlux(
                  std::max(0.0,spectrum_left[OffsetPlus(j)]), volume_left,
                  segment_i->GetBegin(), field_line_idx, dt,
                  /*limit_to_source_energy=*/true);
            }
          }
        }
        else {
          // Inner-boundary source for E+(k).  The total injected W+ follows the
          // same turbulence-level prescription as the integrated model, then is
          // distributed over k with the Kolmogorov log-bin weights.  This keeps
          // the boundary condition compatible with the legacy scalar parameter.
          double V_A = 0.0, area = 0.0;
          if (GetFaceAlfvenSpeedAndArea(segment_i->GetBegin(),field_line_idx,V_A,area)) {
            double* B_inner = segment_i->GetBegin()->GetDatum_ptr(PIC::FieldLine::DatumAtVertexMagneticField);
            if (B_inner) {
              const double Bmag = std::sqrt(B_inner[0]*B_inner[0]+B_inner[1]*B_inner[1]+B_inner[2]*B_inner[2]);
              const double W_total_bc = TurbulenceLevelBeginning*TurbulenceLevelBeginning*Bmag*Bmag/(2.0*kMu0);
              const double W_plus_bc_bin = 0.5*W_total_bc*kolmogorov_weights[j];
              delta_plus[idx] += V_A*W_plus_bc_bin*area*dt;
            }
          }
        }

        // E-(k_j): inward propagation, decreasing segment index.  The last cell
        // is a fixed W-(k) reservoir and is not depleted by its own inward flux.
        if (i != nseg-1) {
          delta_minus[idx] -= ComputeAdvectiveFlux(
              Eminus_i, volume_i, segment_i->GetBegin(), field_line_idx, dt,
              /*limit_to_source_energy=*/true);
        }

        if (i < nseg-1) {
          PIC::FieldLine::cFieldLineSegment* segment_right = field_line->GetSegment(i+1);
          if (segment_right) {
            if (i+1 == nseg-1 &&
                field_line_idx < static_cast<int>(RightBoundaryInitialSpectrumIsSet.size()) &&
                RightBoundaryInitialSpectrumIsSet[field_line_idx]) {
              delta_minus[idx] += ComputeAdvectiveFluxFromDensity(
                  RightBoundaryInitialWminusSpectrum[field_line_idx][j],
                  segment_i->GetEnd(), field_line_idx, dt);
            }
            else {
              double* spectrum_right = segment_right->GetDatum_ptr(SpectralWaveEnergy);
              const double volume_right = SEP::FieldLine::GetSegmentVolume(segment_right,field_line_idx);
              if (spectrum_right && volume_right > 0.0) {
                delta_minus[idx] += ComputeAdvectiveFlux(
                    std::max(0.0,spectrum_right[OffsetMinus(j)]), volume_right,
                    segment_i->GetEnd(), field_line_idx, dt,
                    /*limit_to_source_energy=*/true);
              }
            }
          }
        }
      }
    }

    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment || segment->Thread != PIC::ThisThread) continue;

      double* spectrum = segment->GetDatum_ptr(SpectralWaveEnergy);
      if (!spectrum) continue;

      for (int j=0; j<NK; ++j) {
        const std::size_t idx = static_cast<std::size_t>(i)*NK + j;
        spectrum[OffsetPlus(j)] = std::max(0.0,spectrum[OffsetPlus(j)] + delta_plus[idx]);
        spectrum[OffsetMinus(j)] = std::max(0.0,spectrum[OffsetMinus(j)] + delta_minus[idx]);
      }
    }

    processed_field_lines++;
  }

  EnforceRightBoundarySpectrumInitialCondition();
  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergy);

  if (PIC::ThisThread == 0) {
    std::cout << "Wave-number-resolved turbulence advection completed for "
              << processed_field_lines << " field lines; each of " << NK
              << " k-bins was advected independently.\n";
  }
}

double GetGlobalMaxStableTimeStep() {
  // The CFL condition does not depend on k in the present nondispersive Alfvénic
  // model, so the spectral solver can reuse the integrated solver's global CFL
  // estimate exactly.
  return SEP::AlfvenTurbulence_Kolmogorov::GetGlobalMaxStableTimeStep();
}

void ReflectSpectrumAllFieldLines(
    double dt,
    double C_reflection,
    double grad_floor,
    bool enable_logging) {

  if (dt <= 0.0 || C_reflection <= 0.0) return;

  int processed_segments = 0;
  int updated_bins = 0;
  double total_abs_exchanged = 0.0;

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment || segment->Thread != PIC::ThisThread) continue;

      double* spectrum = segment->GetDatum_ptr(SpectralWaveEnergy);
      double* exchange_rates = segment->GetDatum_ptr(SpectralWaveEnergyExchangeRate);
      if (!spectrum) continue;

      double G_reflection = 0.0;
      if (!GetSegmentReflectionRate(segment,C_reflection,grad_floor,G_reflection)) continue;

      // Reflection is a branch-conversion process, not a k-space cascade.  For
      // each wave number, keep the bin energy E_tot(k_j)=E_+(k_j)+E_-(k_j)
      // fixed and relax the branch imbalance toward zero:
      //
      //   d/dt [E_+ - E_-] = -2 G_R [E_+ - E_-] .
      //
      // The exact exponential update below is positive and conservative for
      // every bin separately.  It avoids the previous approximation used in the
      // wave-number-resolved mode, where the integrated E_± were reflected first
      // and the result was projected back onto the spectrum.
      const double factor = std::exp(-2.0*G_reflection*dt);

      for (int j=0; j<NK; ++j) {
        const double Eplus_old = std::max(0.0,spectrum[OffsetPlus(j)]);
        const double Eminus_old = std::max(0.0,spectrum[OffsetMinus(j)]);
        const double Etot = Eplus_old + Eminus_old;
        if (Etot <= 0.0) continue;

        const double imbalance_old = Eplus_old - Eminus_old;
        const double imbalance_new = imbalance_old*factor;

        const double Eplus_new = std::max(0.0,0.5*(Etot+imbalance_new));
        const double Eminus_new = std::max(0.0,Etot-Eplus_new);

        spectrum[OffsetPlus(j)] = Eplus_new;
        spectrum[OffsetMinus(j)] = Eminus_new;

        // Reflection changes the propagation branch at fixed wave number.  The
        // diagnostic records the signed change of each branch in the same bin.
        // For a perfectly conservative reflection update, the plus and minus
        // rates at a given k-bin should sum to zero, apart from roundoff and any
        // boundary enforcement applied after the operator.
        AddIntegratedEnergyExchangeRate(exchange_rates,BranchPlus(),RateReflection,j,
            Eplus_new-Eplus_old,dt);
        AddIntegratedEnergyExchangeRate(exchange_rates,BranchMinus(),RateReflection,j,
            Eminus_new-Eminus_old,dt);

        total_abs_exchanged += std::fabs(Eplus_new-Eplus_old) + std::fabs(Eminus_new-Eminus_old);
        updated_bins++;
      }

      processed_segments++;
    }
  }

  // The right boundary is prescribed as a pre-existing W-(k) reservoir.  Keep it
  // fixed even after reflection, because otherwise reflection inside the last
  // segment could slowly change the imposed boundary condition itself.
  EnforceRightBoundarySpectrumInitialCondition();
  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergy);

  if (enable_logging && PIC::ThisThread == 0) {
    std::cout << "Wave-number-resolved reflection processed " << processed_segments
              << " local segments and " << updated_bins
              << " spectral bins; total absolute branch exchange = "
              << total_abs_exchanged << " J.\n";
  }
}

void CascadeSpectrumAllFieldLines(
    double dt,
    double C_nl,
    double lambda_perp_m,
    bool enable_cross_helicity_modulation,
    bool two_sweep_imex,
    bool enable_logging) {

  if (dt <= 0.0 || C_nl <= 0.0 || lambda_perp_m <= 0.0) return;

  const int nsweep = two_sweep_imex ? 2 : 1;
  const double dt_sweep = dt/static_cast<double>(nsweep);

  int processed_segments = 0;
  double transferred_energy = 0.0;
  double dissipated_energy = 0.0;

  for (int sweep=0; sweep<nsweep; ++sweep) {
    for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
      PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
      if (!field_line) continue;

      const int nseg = field_line->GetTotalSegmentNumber();
      for (int i=0; i<nseg; ++i) {
        PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
        if (!segment || segment->Thread != PIC::ThisThread) continue;

        double* spectrum = segment->GetDatum_ptr(SpectralWaveEnergy);
        double* exchange_rates = segment->GetDatum_ptr(SpectralWaveEnergyExchangeRate);
        if (!spectrum) continue;

        const double volume = SEP::FieldLine::GetSegmentVolume(segment,field_line_idx);
        if (!(volume > 0.0)) continue;

        double rho = 0.0;
        if (!GetSegmentMassDensity(segment,rho)) continue;

        // The sweep is local to one segment.  It moves energy from low k to high
        // k within the same segment and branch.  It does not advect energy along
        // the field line; spatial transport has already been handled by
        // AdvectSpectrumAllFieldLines().  This operator therefore represents the
        // nonlinear spectral transfer/dissipation part of the turbulence model.
        CascadeOneSpectralSweep(
            spectrum,exchange_rates,volume,rho,dt_sweep,dt,C_nl,lambda_perp_m,
            enable_cross_helicity_modulation,transferred_energy,dissipated_energy);

        if (sweep == 0) processed_segments++;
      }
    }
  }

  // The last-segment W-(k) boundary condition represents imposed pre-existing
  // inward turbulence.  Do not let the local cascade operator deplete or reshape
  // that boundary spectrum.
  EnforceRightBoundarySpectrumInitialCondition();
  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergy);

  if (enable_logging && PIC::ThisThread == 0) {
    std::cout << "Wave-number-resolved spectral cascade processed "
              << processed_segments << " local segments with " << nsweep
              << " sweep(s); energy moved to higher k = " << transferred_energy
              << " J, unresolved dissipation at k_max = " << dissipated_energy
              << " J.\n";
  }
}


void OutputSpectrumTecplot2D(long int iteration, double simulation_time) {
  // The writer is meaningful only for the spectral model.  Keeping the guard in
  // the routine itself makes accidental calls harmless when the legacy
  // branch-integrated turbulence model is selected from the CLI.
  if (!IsActive()) return;

  // Synchronize the per-operator exchange-rate diagnostic before rank 0 writes
  // the file.  Unlike the compact spectral energy datum, this rate datum is not
  // used by the physics operators on other ranks; it exists only for output.
  // Doing the gather here keeps the output routine self-contained and avoids
  // requiring every call site to remember to gather both the spectrum and the
  // diagnostic-rate arrays.
  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergyExchangeRate);

  // Only one MPI rank should create the diagnostic file.  The spectral datum and
  // the exchange-rate datum have now been synchronized, so rank 0 can traverse
  // FieldLinesAll and print a complete Tecplot view without every rank opening
  // the same file.
  if (PIC::ThisThread != 0) return;

  // Store spectral diagnostics in a dedicated subdirectory of the usual AMPS
  // output directory.  That keeps the potentially large 2-D spectrum files
  // separate from the compact field-line output amps.FieldLines.out=*.dat.
  const char* base_dir = (PIC::OutputDataFileDirectory && PIC::OutputDataFileDirectory[0] != '\0')
                         ? PIC::OutputDataFileDirectory : ".";

  char spectrum_dir[1024];
  std::snprintf(spectrum_dir,sizeof(spectrum_dir),"%s/turbulence_wave_number_resolved",base_dir);

  // POSIX mkdir is used for consistency with the rest of the AMPS/SEP source
  // tree.  Ignore EEXIST-like failures here: if the directory already exists,
  // the following file open will succeed; if creation truly failed, the open
  // below will report the problem.
  mkdir(spectrum_dir,0777);

  char fname[1200];
  std::snprintf(fname,sizeof(fname),"%s/spectrum.iter=%08ld.dat",spectrum_dir,iteration);

  std::ofstream fout(fname);
  if (!fout.is_open()) {
    std::cerr << "WARNING: could not open wave-number-resolved turbulence spectrum output file '"
              << fname << "' for writing.\n";
    return;
  }

  // Put the shock location directly into the Tecplot title.  This mirrors the
  // field-line title hook used for amps.FieldLines.out=*.dat and lets a Tecplot
  // animation identify where the CME/shock is relative to the spectral
  // turbulence features at the same simulation time.
  double r_shock = 0.0;
  const char* shock_model_name = "unknown";

  switch (SEP::ShockModelType) {
  case SEP::cShockModelType::Analytic1D:
    r_shock = SEP::ParticleSource::ShockWave::Tenishev2005::rShock;
    shock_model_name = "analytic-1D";
    break;
  case SEP::cShockModelType::SwCme1d:
    r_shock = SEP::SW1DAdapter::gState.r_sh_m;
    shock_model_name = "SW-CME-1D";
    break;
  }

  fout.setf(std::ios::scientific);
  fout << std::setprecision(16);

  fout << "TITLE=\"Wave-number-resolved Alfven turbulence spectrum; "
       << "time=" << simulation_time << " s; "
       << "shock model=" << shock_model_name << "; "
       << "R_sh=" << r_shock << " m = " << (r_shock/_AU_)
       << " AU = " << (r_shock/_SUN__RADIUS_) << " R_s\"\n";

  fout << "VARIABLES=\"s[m]\",\"s[AU]\",\"k[m^-1]\",\"log10(k)\","
       << "\"W+[J/m^3 per log-k bin]\",\"W-[J/m^3 per log-k bin]\",\"sigma_c(k)\","
       << "\"dW+dt_cascade[J/m^3/s per log-k bin]\",\"dW-dt_cascade[J/m^3/s per log-k bin]\","
       << "\"dW+dt_particle_excitation[J/m^3/s per log-k bin]\",\"dW-dt_particle_excitation[J/m^3/s per log-k bin]\","
       << "\"dW+dt_particle_damping[J/m^3/s per log-k bin]\",\"dW-dt_particle_damping[J/m^3/s per log-k bin]\","
       << "\"dW+dt_particle_net[J/m^3/s per log-k bin]\",\"dW-dt_particle_net[J/m^3/s per log-k bin]\","
       << "\"dW+dt_reflection[J/m^3/s per log-k bin]\",\"dW-dt_reflection[J/m^3/s per log-k bin]\","
       << "\"dW+dt_total[J/m^3/s per log-k bin]\",\"dW-dt_total[J/m^3/s per log-k bin]\","
       << "\"dWdt_total_both[J/m^3/s per log-k bin]\","
       << "\"nu_total_plus[s^-1]\",\"nu_total_minus[s^-1]\",\"nu_total_both[s^-1]\"\n";

  // Each field line is written as a separate ordered 2-D zone.  The horizontal
  // Tecplot index I corresponds to segment-center distance from the beginning of
  // the field line.  The vertical index J corresponds to the logarithmic k-grid.
  // Data are written with J as the outer loop and I as the inner loop, so the
  // first row of a zone is W(s,k_min) along the field line, followed by the next
  // wave-number row.
  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    if (nseg <= 0) continue;

    std::vector<double> s_center(nseg,0.0);
    std::vector<double> volume(nseg,0.0);
    std::vector<double*> spectrum_ptr(nseg,nullptr);
    std::vector<double*> rate_ptr(nseg,nullptr);

    double s_begin = 0.0;
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment) {
        s_center[i] = s_begin;
        volume[i] = 0.0;
        spectrum_ptr[i] = nullptr;
        rate_ptr[i] = nullptr;
        continue;
      }

      double segment_length = 0.0;
      PIC::FieldLine::cFieldLineVertex* begin_vertex = segment->GetBegin();
      PIC::FieldLine::cFieldLineVertex* end_vertex = segment->GetEnd();
      if (begin_vertex && end_vertex) {
        double* xb = begin_vertex->GetX();
        double* xe = end_vertex->GetX();
        if (xb && xe) {
          const double dx = xe[0]-xb[0];
          const double dy = xe[1]-xb[1];
          const double dz = xe[2]-xb[2];
          segment_length = std::sqrt(dx*dx + dy*dy + dz*dz);
        }
      }

      // The spectrum is stored as segment-integrated energy per logarithmic
      // wave-number bin.  To output a density, divide each bin by exactly the
      // same magnetic-tube segment volume used for the compact W+,W- diagnostics.
      s_center[i] = s_begin + 0.5*segment_length;
      volume[i] = SEP::FieldLine::GetSegmentVolume(segment,field_line_idx);
      spectrum_ptr[i] = segment->GetDatum_ptr(SpectralWaveEnergy);
      rate_ptr[i] = segment->GetDatum_ptr(SpectralWaveEnergyExchangeRate);

      s_begin += segment_length;
    }

    fout << "ZONE T=\"FieldLine " << field_line_idx << "\", I=" << nseg
         << ", J=" << NK << ", DATAPACKING=POINT\n";

    for (int j=0; j<NK; ++j) {
      const double k = KCenter(j);
      const double log10k = (k > 0.0) ? std::log10(k) : -300.0;

      for (int i=0; i<nseg; ++i) {
        double Wplus = 0.0;
        double Wminus = 0.0;

        if (spectrum_ptr[i] && volume[i] > 0.0) {
          Wplus  = std::max(0.0,spectrum_ptr[i][OffsetPlus(j)]) / volume[i];
          Wminus = std::max(0.0,spectrum_ptr[i][OffsetMinus(j)]) / volume[i];
        }

        const double Wsum = Wplus + Wminus;
        const double sigma_c = (Wsum > 0.0) ? (Wplus-Wminus)/Wsum : 0.0;

        // Exchange-rate diagnostics are stored as integrated dE/dt values.
        // Convert them to density rates dW/dt by dividing by the same segment
        // volume used for W_+ and W_-.  If a segment is absent or has no rate
        // datum, all diagnostic rates are printed as zero so the Tecplot grid
        // remains rectangular.
        auto rate_density = [&](int branch, SpectralEnergyExchangeRateComponent component) -> double {
          if (!rate_ptr[i] || !(volume[i] > 0.0)) return 0.0;
          const double value = rate_ptr[i][RateOffset(branch,component,j)]/volume[i];
          return std::isfinite(value) ? value : 0.0;
        };

        const double cascade_plus = rate_density(BranchPlus(),RateCascade);
        const double cascade_minus = rate_density(BranchMinus(),RateCascade);

        const double particle_excitation_plus = rate_density(BranchPlus(),RateParticleExcitation);
        const double particle_excitation_minus = rate_density(BranchMinus(),RateParticleExcitation);

        const double particle_damping_plus = rate_density(BranchPlus(),RateParticleDamping);
        const double particle_damping_minus = rate_density(BranchMinus(),RateParticleDamping);

        const double particle_net_plus = rate_density(BranchPlus(),RateParticleNet);
        const double particle_net_minus = rate_density(BranchMinus(),RateParticleNet);

        const double reflection_plus = rate_density(BranchPlus(),RateReflection);
        const double reflection_minus = rate_density(BranchMinus(),RateReflection);

        // The total source/sink rate reported here intentionally includes only
        // local exchange processes requested for this diagnostic: spectral
        // cascade/dissipation, particle excitation/damping, and wave reflection.
        // Spatial advection and boundary fluxes are not included because they
        // redistribute turbulence along the field line rather than represent a
        // local physical exchange between waves, particles, and unresolved scales.
        const double total_plus = cascade_plus + particle_net_plus + reflection_plus;
        const double total_minus = cascade_minus + particle_net_minus + reflection_minus;
        const double total_both = total_plus + total_minus;

        const double nu_plus = NormalizedRate(total_plus,Wplus);
        const double nu_minus = NormalizedRate(total_minus,Wminus);
        const double nu_both = NormalizedRate(total_both,Wsum);

        fout << s_center[i] << ' ' << (s_center[i]/_AU_) << ' '
             << k << ' ' << log10k << ' '
             << Wplus << ' ' << Wminus << ' ' << sigma_c << ' '
             << cascade_plus << ' ' << cascade_minus << ' '
             << particle_excitation_plus << ' ' << particle_excitation_minus << ' '
             << particle_damping_plus << ' ' << particle_damping_minus << ' '
             << particle_net_plus << ' ' << particle_net_minus << ' '
             << reflection_plus << ' ' << reflection_minus << ' '
             << total_plus << ' ' << total_minus << ' ' << total_both << ' '
             << nu_plus << ' ' << nu_minus << ' ' << nu_both << "\n";
      }
    }
  }

  fout.close();

  std::cout << "Wrote wave-number-resolved turbulence spectrum Tecplot file: "
            << fname << "\n";
}

void WaveParticleCouplingManager(PIC::Datum::cDatumStored& IntegratedWaveEnergy, double dt) {
  if (dt <= 0.0) return;

  int processed_segments = 0;
  double total_wave_energy_change = 0.0;

  for (int field_line_idx=0; field_line_idx<PIC::FieldLine::nFieldLine; ++field_line_idx) {
    PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
    if (!field_line) continue;

    const int nseg = field_line->GetTotalSegmentNumber();
    for (int i=0; i<nseg; ++i) {
      PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(i);
      if (!segment || segment->Thread != PIC::ThisThread) continue;

      double* spectrum = segment->GetDatum_ptr(SpectralWaveEnergy);
      double* exchange_rates = segment->GetDatum_ptr(SpectralWaveEnergyExchangeRate);
      double* G_plus = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
      double* G_minus = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
      double* gamma_plus = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::gamma_plus_array);
      double* gamma_minus = segment->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::gamma_minus_array);

      if (!spectrum || !G_plus || !G_minus || !gamma_plus || !gamma_minus) continue;

      for (int j=0; j<NK; ++j) {
        gamma_plus[j] = G_plus[j];
        gamma_minus[j] = G_minus[j];

        const double Eplus_old = std::max(0.0,spectrum[OffsetPlus(j)]);
        const double Eminus_old = std::max(0.0,spectrum[OffsetMinus(j)]);

        const double Eplus_new = Eplus_old * std::exp(2.0*gamma_plus[j]*dt);
        const double Eminus_new = Eminus_old * std::exp(2.0*gamma_minus[j]*dt);

        const double dEplus_wave = Eplus_new - Eplus_old;
        const double dEminus_wave = Eminus_new - Eminus_old;

        spectrum[OffsetPlus(j)] = std::max(0.0,Eplus_new);
        spectrum[OffsetMinus(j)] = std::max(0.0,Eminus_new);

        // Record wave-particle exchange rates before the streaming arrays are
        // cleared.  Positive dE_wave means particle streaming has excited/grown
        // the wave bin.  Negative dE_wave means the wave bin has been damped and
        // the corresponding energy has been transferred to particles.  The
        // Tecplot output contains excitation and damping separately, plus their
        // signed net sum, so users can distinguish where particles amplify waves
        // from where waves energize/damp particles.
        if (dEplus_wave >= 0.0) {
          AddIntegratedEnergyExchangeRate(exchange_rates,BranchPlus(),RateParticleExcitation,j,
              dEplus_wave,dt);
        }
        else {
          AddIntegratedEnergyExchangeRate(exchange_rates,BranchPlus(),RateParticleDamping,j,
              dEplus_wave,dt);
        }
        AddIntegratedEnergyExchangeRate(exchange_rates,BranchPlus(),RateParticleNet,j,
            dEplus_wave,dt);

        if (dEminus_wave >= 0.0) {
          AddIntegratedEnergyExchangeRate(exchange_rates,BranchMinus(),RateParticleExcitation,j,
              dEminus_wave,dt);
        }
        else {
          AddIntegratedEnergyExchangeRate(exchange_rates,BranchMinus(),RateParticleDamping,j,
              dEminus_wave,dt);
        }
        AddIntegratedEnergyExchangeRate(exchange_rates,BranchMinus(),RateParticleNet,j,
            dEminus_wave,dt);

        total_wave_energy_change += dEplus_wave + dEminus_wave;

        // Transfer the equal-and-opposite energy change to particles resonant
        // with the same wave-number bin.  This is the key difference from the
        // integrated model, where all k-bin changes are summed before particle
        // redistribution.
        RedistributeWaveNumberBinEnergyToParticles(segment,-dEplus_wave,+1,j);
        RedistributeWaveNumberBinEnergyToParticles(segment,-dEminus_wave,-1,j);

        // Streaming arrays are one-time accumulators for the current time step.
        G_plus[j] = 0.0;
        G_minus[j] = 0.0;
      }

      processed_segments++;
    }
  }

  EnforceRightBoundarySpectrumInitialCondition();
  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SpectralWaveEnergy);
  UpdateIntegratedEnergyFromSpectrum(IntegratedWaveEnergy);
  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(IntegratedWaveEnergy);

  if (PIC::ThisThread == 0) {
    std::cout << "Wave-number-resolved wave-particle coupling processed "
              << processed_segments << " local segments; total spectral wave-energy change = "
              << total_wave_energy_change << " J.\n";
  }
}

} // namespace WaveNumberResolved
} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP

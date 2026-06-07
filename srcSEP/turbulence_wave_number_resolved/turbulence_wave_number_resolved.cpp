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

  * Operators that have not yet been written in spectral form in the legacy code
    (shock injection, reflection, nonlinear cascade) may still be called on the
    integrated energy.  ProjectIntegratedEnergyToSpectrum() then maps the new
    integrated value back to the spectrum while preserving the previous spectral
    shape if one exists.  This preserves compatibility without pretending that
    those operators are already full spectral cascade/reflection solvers.

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

// Branch helper: the first NK entries are E+(k), the second NK entries are E-(k).
inline int OffsetPlus(int j) { return j; }
inline int OffsetMinus(int j) { return NK+j; }

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

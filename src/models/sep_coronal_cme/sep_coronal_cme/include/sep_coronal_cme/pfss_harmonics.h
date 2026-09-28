#ifndef SEP_CORONAL_CME_PFSS_HARMONICS_H
#define SEP_CORONAL_CME_PFSS_HARMONICS_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <cstdint>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

// Real, orthonormal, Condon--Shortley spherical-harmonic coefficient.  cosine
// multiplies N_lm P_l^m(cos theta) cos(m phi), and sine multiplies the sine
// counterpart. For m=0, sine must be zero. Both values are photospheric radial
// magnetic-field amplitudes in tesla, not potential coefficients.
struct HarmonicCoefficient {
  int degree = 0;
  int order = 0;
  double cosineT = 0.0;
  double sineT = 0.0;
};

enum class MonopolePolicy { Reject, Remove };
enum class FieldLineTopology { OpenToOuterBoundary, ClosedBelowOuterBoundary,
                               Separatrix, Invalid };

struct SphericalField {
  double brT = 0.0, bThetaT = 0.0, bPhiT = 0.0;
  // Analytic partial derivatives. First index is field component (r,theta,phi)
  // and second is coordinate (r,theta,phi); angular derivatives are per radian.
  double derivative[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
};

struct MapSample {
  double thetaRad = 0.0, phiRad = 0.0, radialFieldT = 0.0;
  double solidAngleWeightSr = 0.0;
};

struct ReconstructionMetrics {
  double weightedRmsT = 0.0;
  double unsignedFluxChangeFraction = 0.0;
  double removedMonopoleT = 0.0;
};

struct TopologyPair {
  FieldLineTopology pfssOpenToRb = FieldLineTopology::Invalid;
  FieldLineTopology compositeOpenToRi = FieldLineTopology::Invalid;
  std::uint64_t generation = 0;
  std::string plasmaAuthority;
  std::string targetSpeedTopologyAuthority;
};

class PfssHarmonics {
 public:
  static Core::Result<PfssHarmonics> Create(double solarRadiusM,
      double sourceSurfaceRadiusM, std::vector<HarmonicCoefficient> coefficients,
      MonopolePolicy monopolePolicy = MonopolePolicy::Reject,
      double monopoleToleranceT = 1.0e-15);

  Core::Result<SphericalField> Evaluate(double radiusM, double thetaRad,
                                        double phiRad) const;
  Core::Result<Vec3> EvaluateCartesian(Vec3 positionM) const;
  Core::Result<FieldLineTopology> Classify(Vec3 positionM, double stepM,
                                           int maximumSteps = 200000) const;

  // Pure transforms always return a new set. The retained raw coefficients
  // are immutable, preventing a Stage-3 calibration pass from scaling twice.
  PfssHarmonics Scaled(double factor) const;
  PfssHarmonics HeatKernelFiltered(int apodizationDegree) const;
  const std::vector<HarmonicCoefficient>& Coefficients() const noexcept {
    return coefficients_;
  }
  double SolarRadiusM() const noexcept { return solarRadiusM_; }
  double SourceSurfaceRadiusM() const noexcept { return sourceSurfaceRadiusM_; }

  static Core::Result<std::vector<HarmonicCoefficient>> ProjectMap(
      const std::vector<MapSample>& samples, int maximumDegree,
      MonopolePolicy policy, double monopoleToleranceT,
      ReconstructionMetrics* metrics = nullptr);
  static Core::Result<ReconstructionMetrics> CompareMap(
      const std::vector<MapSample>& samples,
      const std::vector<HarmonicCoefficient>& coefficients);

 private:
  double solarRadiusM_ = 0.0, sourceSurfaceRadiusM_ = 0.0;
  std::vector<HarmonicCoefficient> coefficients_;
};

// Keeps both topology authorities explicit when R_i<R_b. This helper does not
// infer composite topology from the unused PFSS continuation.
Core::Result<TopologyPair> RouteTopology(FieldLineTopology pfss,
    FieldLineTopology composite, std::uint64_t generation,
    const std::string& targetSpeedTopologyAuthority);

} }  // namespace SEP::CoronalCME
#endif

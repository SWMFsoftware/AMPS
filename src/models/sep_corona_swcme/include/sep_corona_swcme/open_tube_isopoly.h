#ifndef SEP_CORONA_SWCME_OPEN_TUBE_ISOPOLY_H
#define SEP_CORONA_SWCME_OPEN_TUBE_ISOPOLY_H

#include "sep_corona_swcme/ambient_state.h"
#include "sep_coronal_cme/tube_isopoly_wind.h"

namespace SEP::CoronaSwcme {

struct OpenTubeThermalBand {
  double polytropicIndex=1, blendBeginRadiusM=0, blendEndRadiusM=0;
};

struct OpenTubeIsoPolyInputs {
  // Mass loading is rho*u/|B|, in kg/(s Wb), from one declared anchor.
  // The reference density anchors the outer species temperature laws. The
  // temperatures at that density are read from the event's composition.
  double massLoadingKgPerSWb=0, referenceDensityKgM3=0;
  OpenTubeThermalBand electron, proton, alpha;
  double angularTraceToleranceRad=0;
  unsigned criticalScanIntervals=256;
  double relativeWindTolerance=1e-9;
};

struct OpenTubeIsoPolyRequest : OpenTubeIsoPolyInputs {
  // A tube is labelled by its source-surface colatitude in radians.
  double sourceColatitudeRad=0;
};

struct OpenTubeIsoPolySample {
  AmbientPrimitive primitive;
  // Scalar spatial derivatives are with respect to heliocentric radius.
  // Net heat uses the physical radial advection rate for those derivatives.
  CoronalCME::IsoPolyTubeSample balance;
  CoronalCME::Vec3 positionM;
};

// Immutable resolved tube for an axisymmetric PFSS sector and its rotating
// Parker continuation. It can be evaluated at any local inertial azimuth at
// the solved radius/colatitude. A single resolved tube does not supply a
// spatial atlas, transverse derivatives, or a native material lifecycle.
class OpenTubeIsoPolyWind {
 public:
  static Core::Result<OpenTubeIsoPolyWind> Prepare(std::shared_ptr<const AmbientModel> ambient,
      const OpenTubeIsoPolyRequest& request);
  Core::Result<OpenTubeIsoPolySample> Evaluate(double radiusM,
      double azimuthRad,double epochS) const;
  double CriticalRadiusM() const { return wind_.CriticalCoordinateM(); }
  bool UsesAmbient(const std::shared_ptr<const AmbientModel>& ambient) const {
    return ambient&&ambient==ambient_;
  }

 private:
  friend class OpenTubeIsoPolyFamily;
  static Core::Result<OpenTubeIsoPolyWind> PrepareThrough(
      std::shared_ptr<const AmbientModel> ambient,const OpenTubeIsoPolyInputs&,
      double radiusM,double colatitudeRad);
  struct Geometry;
  std::shared_ptr<const AmbientModel> ambient_;
  std::shared_ptr<const Geometry> geometry_;
  CoronalCME::IsoPolyTubeWind wind_;
};

struct OpenTubeIsoPolyFamilyInput {
  // A declared single open population shares these thermal and mass-loading
  // inputs. The band is measured at the source surface, not at the query
  // radius, so its boundaries remain actual magnetic/material surfaces.
  OpenTubeIsoPolyInputs population;
  double lowerSourceColatitudeRad=0, upperSourceColatitudeRad=0;
  std::string inputIdentity;
};

struct OpenTubeIsoPolySpatialSample {
  AmbientState ambient;
  // Species order is electron, proton and (when present) alpha, W/m^3.
  std::vector<double> speciesNetHeatWPerM3;
};

// Controlled continuous field on a monotone, nonpolar, axisymmetric open
// band. A query traces its own tube and solves the regular wind branch.
// There is no interpolation across a neutral line and no fallback provider.
class OpenTubeIsoPolyFamily final {
 public:
  static Core::Result<std::shared_ptr<const OpenTubeIsoPolyFamily>> Prepare(
      std::shared_ptr<const AmbientModel>,const OpenTubeIsoPolyFamilyInput&);
  const std::string& Identity() const { return identity_; }
  Core::Result<OpenTubeIsoPolySample> Evaluate(CoronalCME::Vec3,double epochS) const;
  Core::Result<OpenTubeIsoPolySpatialSample> EvaluateWithDerivatives(
      CoronalCME::Vec3,double epochS,std::uint64_t generation) const;

 private:
  OpenTubeIsoPolyFamily() = default;
  std::shared_ptr<const AmbientModel> ambient_;
  OpenTubeIsoPolyFamilyInput input_;
  std::string identity_;
  int sector_=0;
};

}  // namespace SEP::CoronaSwcme

#endif

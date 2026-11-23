#ifndef SEP_CORONA_SWCME_OPEN_TUBE_ISOPOLY_H
#define SEP_CORONA_SWCME_OPEN_TUBE_ISOPOLY_H

#include "sep_corona_swcme/ambient_state.h"
#include "sep_coronal_cme/tube_isopoly_wind.h"

namespace SEP::CoronaSwcme {

struct OpenTubeThermalBand {
  double polytropicIndex=1, blendBeginRadiusM=0, blendEndRadiusM=0;
};

struct OpenTubeIsoPolyRequest {
  // A tube is labelled by its source-surface colatitude. Angles are radians.
  // Mass loading is rho*u/|B|, in kg/(s Wb), from one declared anchor.
  // The reference density anchors the outer species temperature laws. The
  // temperatures at that density are read from the event's composition.
  double sourceColatitudeRad=0, massLoadingKgPerSWb=0, referenceDensityKgM3=0;
  OpenTubeThermalBand electron, proton, alpha;
  double angularTraceToleranceRad=0;
  unsigned criticalScanIntervals=256;
  double relativeWindTolerance=1e-9;
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
  struct Geometry;
  std::shared_ptr<const AmbientModel> ambient_;
  std::shared_ptr<const Geometry> geometry_;
  CoronalCME::IsoPolyTubeWind wind_;
};

}  // namespace SEP::CoronaSwcme

#endif

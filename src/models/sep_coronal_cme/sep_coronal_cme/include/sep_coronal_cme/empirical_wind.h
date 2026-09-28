#ifndef SEP_CORONAL_CME_EMPIRICAL_WIND_H
#define SEP_CORONAL_CME_EMPIRICAL_WIND_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <map>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

struct ProfileNode { double x=0.0,value=0.0,first=0.0,second=0.0; };

// Positive physical profiles are represented by an exact quintic Hermite
// polynomial in z=ln(y/y_ref), including supplied first/second derivatives.
class CertifiedPositiveProfile {
 public:
  static Core::Result<CertifiedPositiveProfile> Create(
      std::vector<ProfileNode> nodes,double referenceValue,
      bool requireNoOvershoot=true);
  Core::Result<ProfileNode> Evaluate(double x) const;
 private:
  double reference_=1.0;
  std::vector<ProfileNode> nodes_;
  std::vector<std::vector<double>> coefficients_;
};

struct TwoZoneState { double densityKgM3=0.0,speedMPerS=0.0,logMismatch=0.0; };
Core::Result<TwoZoneState> BlendTwoZoneWind(double radiusM,double joinInnerM,
    double joinOuterM,double innerDensityKgM3,double outerSpeedMPerS,
    double magneticFieldT,double massLoadingKgPerSWb,double referenceDensityKgM3);

enum class VelocityComponent { Radial, FieldAligned };
enum class VelocityFrame { Inertial, Corotating };
enum class ProfileAbscissa { HeliocentricRadius, OrientedArcLength };
Core::Result<double> ConvertToCorotatingFieldAlignedSpeed(double storedSpeedMPerS,
    VelocityComponent component,VelocityFrame frame,double radialProjection,
    double minimumRadialProjection,Vec3 omegaRadPerS,Vec3 positionM,Vec3 tangent);

struct ConsumerMeasure {
  std::string stableId;
  double openFluxWb=0.0,openAreaM2=0.0,sourceNumberRate=0.0;
  double observerExposureM2S=0.0,exportLengthM=0.0;
  bool covered=false;
  std::string rejectionReason;
};
struct CoverageCensus {
  double coveredFluxFraction=0.0,coveredAreaFraction=0.0,
      coveredSourceFraction=0.0,coveredObserverFraction=0.0,
      coveredExportFraction=0.0;
  std::vector<std::string> rejectedIds;
};
Core::Result<CoverageCensus> BuildCoverageCensus(
    const std::vector<ConsumerMeasure>& consumers,bool eventNominal);

} }  // namespace SEP::CoronalCME
#endif

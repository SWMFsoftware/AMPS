#ifndef _GENERIC_COMET_BRIGHTNESS_MATH_H_
#define _GENERIC_COMET_BRIGHTNESS_MATH_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Comet {
namespace Brightness {
namespace Math {

// Rotate a velocity stored by the DIM1 spherical mover back into a Cartesian
// frame at the supplied unit radial direction.  This is the exact transpose
// (inverse) of NormalizeSphericalParticleState() in src/pic/pic_mover.cpp.
// Keeping the operation here makes the brightness estimator testable without
// a configured PIC build and prevents a subtly different LOS basis from being
// invented in the application sampler.
inline bool RotateLocalVelocityToCartesian(const double localVelocity[3],
    const double radialDirection[3],double cartesianVelocity[3]) {
  if (localVelocity==NULL || radialDirection==NULL ||
      cartesianVelocity==NULL) return false;
  for (int i=0;i<3;i++)
    if (!std::isfinite(localVelocity[i]) ||
        !std::isfinite(radialDirection[i])) return false;

  const double norm=std::sqrt(radialDirection[0]*radialDirection[0]+
      radialDirection[1]*radialDirection[1]+
      radialDirection[2]*radialDirection[2]);
  if (!(norm>0.0) || !std::isfinite(norm)) return false;
  const double nx=radialDirection[0]/norm;
  const double ny=radialDirection[1]/norm;
  const double nz=radialDirection[2]/norm;
  const double transverse=std::sqrt(nx*nx+ny*ny);
  double cosZ=1.0,sinZ=0.0,cosY=0.0,sinY=nz;
  if (transverse>1.0e-20) {
    cosZ=nx/transverse;
    sinZ=ny/transverse;
    cosY=transverse;
  }

  // Columns of this matrix are the local radial and two tangential basis
  // vectors.  Multiplication by those columns is the inverse of the mover's
  // Cartesian-to-local dot products.
  cartesianVelocity[0]=cosY*cosZ*localVelocity[0]-
      sinZ*localVelocity[1]-sinY*cosZ*localVelocity[2];
  cartesianVelocity[1]=cosY*sinZ*localVelocity[0]+
      cosZ*localVelocity[1]-sinY*sinZ*localVelocity[2];
  cartesianVelocity[2]=sinY*localVelocity[0]+
      cosY*localVelocity[2];
  return std::isfinite(cartesianVelocity[0]) &&
      std::isfinite(cartesianVelocity[1]) &&
      std::isfinite(cartesianVelocity[2]);
}

inline double Dot3(const double a[3],const double b[3]) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

// One DIM1 AMPS cell represents a spherical shell with a completed-window
// number density.  Keeping this small value type independent of PIC and MPI
// makes the geometry and unit conversions testable without a configured AMPS
// build.
struct cShell {
  double InnerRadius;
  double OuterRadius;
  double NumberDensity;
};

// Convert an integrated radiated power per molecule to photon events per
// molecule per second at an application-selected effective wavelength.  This
// operation is separate from every core fluorescence model because a band
// energy g factor does not uniquely define a photon count.  Invalid or
// non-positive inputs fail visibly with NaN.
inline double EnergyRateToPhotonRate(double emittedPowerW,
    double effectiveWavelengthM) {
  constexpr double PlanckJS=6.62607015e-34;
  constexpr double LightSpeedMPS=299792458.0;
  if (!std::isfinite(emittedPowerW) || !(emittedPowerW>0.0) ||
      !std::isfinite(effectiveWavelengthM) ||
      !(effectiveWavelengthM>0.0))
    return std::numeric_limits<double>::quiet_NaN();
  return emittedPowerW/(PlanckJS*LightSpeedMPS/effectiveWavelengthM);
}

// Return the visible path length through one shell for a ray whose projected
// cometocentric distance is impactParameter.  An opaque spherical nucleus
// blocks the rear half of the ray only when b<Rn.  At b=Rn the ray is tangent
// to the nucleus, so the measure-zero contact does not hide the rear coma.
// Shells are clipped to the physical nucleus surface as a defensive measure;
// normal radial meshes already start exactly at Rn.
inline double VisibleShellChordLength(double innerRadius,double outerRadius,
    double impactParameter,double nucleusRadius) {
  if (!std::isfinite(innerRadius) || !std::isfinite(outerRadius) ||
      !std::isfinite(impactParameter) || !std::isfinite(nucleusRadius) ||
      innerRadius<0.0 || !(outerRadius>innerRadius) ||
      impactParameter<0.0 || nucleusRadius<0.0 ||
      impactParameter>=outerRadius) {
    return 0.0;
  }

  const double clippedInner=std::max(innerRadius,nucleusRadius);
  if (!(outerRadius>clippedInner)) return 0.0;

  const double innerAlongRay=std::max(clippedInner,impactParameter);
  const double outerSquared=outerRadius*outerRadius-
      impactParameter*impactParameter;
  const double innerSquared=innerAlongRay*innerAlongRay-
      impactParameter*impactParameter;
  // max(0,...) protects a valid tangent from a small negative roundoff.
  const double halfChord=std::sqrt(std::max(0.0,outerSquared))-
      std::sqrt(std::max(0.0,innerSquared));
  if (!(halfChord>0.0)) return 0.0;

  return (impactParameter<nucleusRadius ? 1.0 : 2.0)*halfChord;
}

// Piecewise-constant Abel sum used by the production writer.  The summation
// order is the deterministic radial mesh order.  A negative/nonfinite sample
// is rejected with NaN rather than being hidden in a brightness product.
inline double ColumnDensity(const std::vector<cShell> &shells,
    double impactParameter,double nucleusRadius) {
  if (!std::isfinite(impactParameter) || impactParameter<0.0 ||
      !std::isfinite(nucleusRadius) || nucleusRadius<0.0)
    return std::numeric_limits<double>::quiet_NaN();

  double column=0.0;
  for (std::size_t shell=0;shell<shells.size();shell++) {
    if (!std::isfinite(shells[shell].InnerRadius) ||
        !std::isfinite(shells[shell].OuterRadius) ||
        shells[shell].InnerRadius<0.0 ||
        !(shells[shell].OuterRadius>shells[shell].InnerRadius) ||
        !std::isfinite(shells[shell].NumberDensity) ||
        shells[shell].NumberDensity<0.0)
      return std::numeric_limits<double>::quiet_NaN();
    column+=shells[shell].NumberDensity*VisibleShellChordLength(
        shells[shell].InnerRadius,shells[shell].OuterRadius,
        impactParameter,nucleusRadius);
  }
  return column;
}

// A Rayleigh is 10^6 photons cm^-2 s^-1 integrated over 4*pi sr.  For an
// optically thin isotropic fluorescence rate g [photon molecule^-1 s^-1],
// B[R]=g*N[cm^-2]/10^6=g*N[m^-2]/10^10.
inline double OpticallyThinRayleigh(double columnDensityM2,
    double photonGFactorPerSecond) {
  if (!std::isfinite(columnDensityM2) || columnDensityM2<0.0 ||
      !std::isfinite(photonGFactorPerSecond) ||
      photonGFactorPerSecond<0.0)
    return std::numeric_limits<double>::quiet_NaN();
  return photonGFactorPerSecond*columnDensityM2/1.0e10;
}

// This is a diagnostic only.  The caller supplies a conservative line-center
// cross section and the code reports tau_0=sigma_0*N and a threshold flag;
// it never invents an escape probability or attenuates the thin brightness.
inline double LineCenterOpticalDepth(double columnDensityM2,
    double lineCenterCrossSectionM2) {
  if (!std::isfinite(columnDensityM2) || columnDensityM2<0.0 ||
      !std::isfinite(lineCenterCrossSectionM2) ||
      lineCenterCrossSectionM2<0.0)
    return std::numeric_limits<double>::quiet_NaN();
  return lineCenterCrossSectionM2*columnDensityM2;
}

inline bool OpticalDepthThresholdExceeded(double opticalDepth,
    double threshold) {
  return std::isfinite(opticalDepth) && std::isfinite(threshold) &&
      threshold>0.0 && opticalDepth>=threshold;
}

} // namespace Math
} // namespace Brightness
} // namespace Comet

#endif

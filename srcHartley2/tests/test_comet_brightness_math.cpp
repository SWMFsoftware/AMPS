#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "srcHartley2/CometBrightnessMath.h"

namespace BrightnessMath=Comet::Brightness::Math;

namespace {
bool Close(double actual,double expected,double relative=1.0e-12) {
  return std::fabs(actual-expected)<=relative*std::max(1.0,
      std::max(std::fabs(actual),std::fabs(expected)));
}
}

int main() {
  bool pass=true;
  const double nucleus=2.0;
  const double outer=10.0;

  // The LOS sampler must invert the mover's Cartesian-to-local spherical
  // rotation exactly. Exercise a general direction and both tangential
  // components, then reconstruct the mover dot products independently.
  const double localVelocity[3]={430.0,-71.0,29.0};
  const double radialDirection[3]={2.0,-3.0,6.0};
  double cartesianVelocity[3];
  pass=BrightnessMath::RotateLocalVelocityToCartesian(localVelocity,
      radialDirection,cartesianVelocity) && pass;
  const double norm=7.0;
  const double nx=radialDirection[0]/norm;
  const double ny=radialDirection[1]/norm;
  const double nz=radialDirection[2]/norm;
  const double transverse=std::sqrt(nx*nx+ny*ny);
  const double cosZ=nx/transverse,sinZ=ny/transverse;
  const double reconstructed[3]={
      nx*cartesianVelocity[0]+ny*cartesianVelocity[1]+
          nz*cartesianVelocity[2],
      -sinZ*cartesianVelocity[0]+cosZ*cartesianVelocity[1],
      -nz*cosZ*cartesianVelocity[0]-nz*sinZ*cartesianVelocity[1]+
          transverse*cartesianVelocity[2]};
  for (int i=0;i<3;i++) pass=Close(reconstructed[i],localVelocity[i]) && pass;

  // A central ray sees only the front gas because the opaque nucleus blocks
  // the complete rear segment.
  pass=Close(BrightnessMath::VisibleShellChordLength(nucleus,outer,0.0,
      nucleus),outer-nucleus) && pass;

  // Outside the projected nucleus both halves are visible.  At exact
  // tangency the point contact has zero measure and must not occult the rear.
  const double b=5.0;
  pass=Close(BrightnessMath::VisibleShellChordLength(nucleus,outer,b,nucleus),
      2.0*std::sqrt(outer*outer-b*b)) && pass;
  pass=Close(BrightnessMath::VisibleShellChordLength(nucleus,outer,nucleus,
      nucleus),2.0*std::sqrt(outer*outer-nucleus*nucleus)) && pass;
  pass=BrightnessMath::VisibleShellChordLength(nucleus,outer,outer,nucleus)==
      0.0 && pass;

  // Adjacent piecewise-constant shells telescope exactly for equal density.
  std::vector<BrightnessMath::cShell> shells;
  shells.push_back(BrightnessMath::cShell{2.0,4.0,3.0});
  shells.push_back(BrightnessMath::cShell{4.0,10.0,3.0});
  pass=Close(BrightnessMath::ColumnDensity(shells,0.0,nucleus),
      3.0*(outer-nucleus)) && pass;
  pass=Close(BrightnessMath::ColumnDensity(shells,b,nucleus),
      3.0*2.0*std::sqrt(outer*outer-b*b)) && pass;
  // A requested observational coordinate may lie outside the finite model
  // domain.  It has no represented shell contribution and must be zero, not
  // an implicit analytic continuation.
  pass=BrightnessMath::ColumnDensity(shells,2.0*outer,nucleus)==0.0 && pass;

  // Exercise SI-to-Rayleigh conversion and the deliberately separate tau
  // flag.  No attenuation factor is part of this API.
  const double column=2.5e16;
  const double gFactor=4.0e-5;
  // Case-1 Lowell power at 1.75 AU and the explicitly selected 308.5-nm
  // convention reproduce the independently archived photon rate.
  pass=Close(BrightnessMath::EnergyRateToPhotonRate(
      1.1591836734693877e-22,308.5e-9),1.8002411793823956e-4) && pass;
  pass=Close(BrightnessMath::OpticallyThinRayleigh(column,gFactor),100.0) &&
      pass;
  const double tau=BrightnessMath::LineCenterOpticalDepth(column,2.0e-17);
  pass=Close(tau,0.5) &&
      BrightnessMath::OpticalDepthThresholdExceeded(tau,0.5) &&
      !BrightnessMath::OpticalDepthThresholdExceeded(tau,0.6) && pass;

  shells[0].NumberDensity=-1.0;
  pass=std::isnan(BrightnessMath::ColumnDensity(shells,0.0,nucleus)) && pass;
  pass=std::isnan(BrightnessMath::OpticallyThinRayleigh(-1.0,gFactor)) && pass;
  pass=std::isnan(BrightnessMath::EnergyRateToPhotonRate(1.0,0.0)) && pass;

  std::cout << (pass ? "PASS" : "FAIL")
      << " comet-brightness rotation=1 exact-chord=1 occultation=1 units=1 tau-flag=1 invalid=1\n";
  return pass ? 0 : 1;
}

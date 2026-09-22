#include "../../util/FluxNumerics.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace FN=Earth::FluxNumerics;
using Earth::GridlessMode::TrajectoryTermination;

static int failures=0;

static void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

static bool Near(double actual,double expected,double relative=1.0e-12,double absolute=1.0e-12) {
  return std::fabs(actual-expected)<=absolute+relative*std::max(std::fabs(actual),std::fabs(expected));
}

int main() {
  const double protonMass=1.007276466621*FN::kAtomicMassUnit_kg;
  const double alphaMass=4.001506179127*FN::kAtomicMassUnit_kg;

  // U-F03: conversions must be mutual inverses over nonrelativistic and relativistic
  // energies and for more than one charge-to-mass ratio.
  const double energies[]={0.01,1.0,100.0,1000.0,1.0e5};
  for (double E:energies) {
    const double Rp=FN::RigidityFromEnergyGV(E*FN::kMeVToJ,FN::kElementaryCharge_C,protonMass);
    Check(Near(FN::EnergyFromRigidityMeV(Rp,FN::kElementaryCharge_C,protonMass),E,2.0e-12,1.0e-11),
          "proton energy-rigidity round trip");
    const double Ra=FN::RigidityFromEnergyGV(E*FN::kMeVToJ,2.0*FN::kElementaryCharge_C,alphaMass);
    Check(Near(FN::EnergyFromRigidityMeV(Ra,2.0*FN::kElementaryCharge_C,alphaMass),E,2.0e-12,1.0e-11),
          "alpha energy-rigidity round trip");
  }
  Check(FN::RelativisticSpeed(1.0e9*FN::kMeVToJ,protonMass)<FN::kSpeedOfLight_m_s,
        "relativistic speed remains below c");

  // U-F04: exact endpoints, monotonicity, and the LINEAR midpoint regression.  The
  // midpoint assertion specifically fails for the removed Mode3D a*a expression.
  const std::vector<double> linear=FN::BuildEnergyGridMeV(
      10.0,50.0,5,FN::EnergySpacing::Linear,false,0,0,FN::kElementaryCharge_C,protonMass);
  Check(linear.size()==5 && Near(linear[1],20.0) && Near(linear[2],30.0),
        "linear grid uses one interpolation factor");
  const std::vector<double> logarithmic=FN::BuildEnergyGridMeV(
      1.0,10000.0,5,FN::EnergySpacing::Log,false,0,0,FN::kElementaryCharge_C,protonMass);
  Check(Near(logarithmic[1],10.0) && Near(logarithmic[3],1000.0),"log grid nodes");
  const std::vector<double> rigidity=FN::BuildEnergyGridMeV(
      1.0,10000.0,5,FN::EnergySpacing::Log,true,9,7,FN::kElementaryCharge_C,protonMass);
  Check(rigidity.size()==7 && rigidity.front()==1.0 && rigidity.back()==10000.0,
        "rigidity scan cap and exact endpoints");
  for (std::size_t i=1;i<rigidity.size();++i)
    Check(rigidity[i]>rigidity[i-1],"rigidity scan is strictly monotone");

  // U-F05: quadrature, channel clipping, and density units.  A constant integrand is
  // exact under trapezoidal quadrature, so these checks isolate unit factors and limits.
  const std::vector<double> E{1.0,2.0,3.0};
  const std::vector<double> T{0.5,0.5,0.5};
  const auto unitSpectrum=[](double) { return 1.0; };
  const double expectedFlux=4.0*FN::kPi*0.5*2.0*FN::kMeVToJ;
  Check(Near(FN::IntegrateFlux(E,T,1.0,3.0,unitSpectrum),expectedFlux),"total flux quadrature");
  const double expectedChannel=4.0*FN::kPi*0.5*1.5*FN::kMeVToJ;
  Check(Near(FN::IntegrateFlux(E,T,0.5,2.5,unitSpectrum),expectedChannel),"channel clipping");
  const double density=FN::IntegrateDensity(E,T,protonMass,unitSpectrum);
  Check(std::isfinite(density) && density>0.0,"density quadrature is finite and positive");

  // U-F06: equal-solid-angle weights and deterministic subsampling.
  const std::vector<FN::DirectionSample> directions=FN::BuildEqualSolidAngleDirections(24,48);
  double solidAngle=0.0;
  for (const FN::DirectionSample& d:directions) {
    solidAngle+=d.solidAngleWeight_sr;
    Check(Near(d.x*d.x+d.y*d.y+d.z*d.z,1.0,1.0e-13),"direction is unit length");
  }
  Check(Near(solidAngle,4.0*FN::kPi,2.0e-14),"angular weights sum to 4*pi");
  const std::vector<int> ids{0,1,2,3,4,5,6,7,8,9};
  const std::vector<int> selected=FN::SelectDeterministic(ids,4);
  Check(selected==std::vector<int>({0,3,6,9}),"deterministic sky subsampling");

  // U-F12: unresolved outcomes produce bounds, not false shielding.  Two allowed, one
  // forbidden, and one unresolved direction give nominal 2/3 and bounds [1/2,3/4].
  FN::AccessAccumulator access;
  access.Record(TrajectoryTermination::OuterBoundaryAllowed,1.0);
  access.Record(TrajectoryTermination::OuterBoundaryAllowed,1.0,true);
  access.Record(TrajectoryTermination::InnerBoundaryForbidden,0.0);
  access.Record(TrajectoryTermination::TimeLimit,0.0);
  const FN::AccessEstimate estimate=FN::ResolveAccess(access);
  Check(Near(estimate.nominal,2.0/3.0) && Near(estimate.lower,0.5) && Near(estimate.upper,0.75),
        "resolved estimate and conservative unresolved bounds");
  Check(Near(estimate.unresolvedFraction,0.25) && access.retried==1,
        "unresolved and retry accounting");
  Check(FN::ExceedsUnresolvedTolerance(access,0.20) && !FN::ExceedsUnresolvedTolerance(access,0.25),
        "validation tolerance boundary");
  FN::AccessAccumulator mergedLeft,mergedRight;
  mergedLeft.Record(TrajectoryTermination::OuterBoundaryAllowed,1.0);
  mergedRight.Record(TrajectoryTermination::NumericalFailure,0.0,true);
  mergedLeft.Merge(mergedRight);
  Check(mergedLeft.sampled==2 && mergedLeft.resolved==1 && mergedLeft.retried==1 &&
        mergedLeft.terminationCounts[static_cast<std::size_t>(
            static_cast<int>(TrajectoryTermination::NumericalFailure))]==1,
        "block merge preserves structured accounting");

  FN::AccessAccumulator noneResolved;
  noneResolved.Record(TrajectoryTermination::StepLimit,0.0);
  const FN::AccessEstimate none=FN::ResolveAccess(noneResolved);
  Check(std::isnan(none.nominal) && Near(none.lower,0.0) && Near(none.upper,1.0),
        "no-resolved case is NaN with [0,1] bounds");

  const FN::TransmissionDiagnostics diagnostic=FN::ComputeTransmissionDiagnostics(
      std::vector<double>({10.0,100.0,1000.0}),std::vector<double>({0.0,0.5,1.0}),
      FN::kElementaryCharge_C,protonMass);
  Check(diagnostic.RcUpper_GV>=diagnostic.RcLower_GV && diagnostic.PenumbraWidth_GV>=0.0,
        "transmission diagnostics are ordered");

  if (failures!=0) {
    std::cerr << failures << " FluxNumerics test(s) failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: U-F03 U-F04 U-F05 U-F06 U-F12\n";
  return EXIT_SUCCESS;
}
